#!/usr/bin/env python
from __future__ import print_function
import math

from simtk import unit as u

from simtk import openmm as mm
from simtk.openmm import app

from argparse import ArgumentParser
import sys

parser = ArgumentParser()
group = parser.add_argument_group('Input File Options')
group.add_argument('--xml', dest='xml', default='system.xml', metavar='FILE',
                   help='''OpenMM System XML file. Default is %(default)s''')
group.add_argument('-p', '--pdb', dest='pdb', metavar='<PDB_FILE>', required=True,
                   help='''PDB file with coordinates for all atoms. Is also the
                   reference coordinates''')
group.add_argument('-s', '--state', dest='state', metavar='FILE', default=None,
                   help='''Restart file (any format)''')
group = parser.add_argument_group('Positional Restraint Options')
group.add_argument('-k', '--force-constant', dest='force_constant', type=float,
                   metavar='FLOAT', help='''Force constant for cartesian
                   constraints. Default 10 kcal/mol/A^2''', default=10)
group = parser.add_argument_group('Output File Options')
group.add_argument('-l' , '--log', dest='log', default='run.out',
                   metavar='FILE', help='''Log file''')
#group.add_argument('-x', '--trajectory', dest='trajectory', default='md.nc',
#                   metavar='FILE', help='''NetCDF trajectory to generate.
#                   Snapshots written every --interval steps.''')
group.add_argument('--checkpoint', dest='checkpoint', metavar='FILE',
                   default=None, help='''Name of a checkpoint file to write
                   periodically throughout the simulation. Primarily useful for
                   debugging intermittent and rare errors.''')
group.add_argument('--interval', dest='interval', default=500, metavar='INT',
                   help='Interval between printing state data. Default 500',
                   type=int)
group = parser.add_argument_group('Simulation Options')
group.add_argument('-n', '--num-steps', dest='num_steps', required=True, type=int,
                   help='Number of MD steps to run. Required', metavar='INT')
group.add_argument('--ntp', dest='ntp', default=False, action='store_true',
                   help='Do constant pressure simulation')
group.add_argument('--aniso', dest='aniso', default=False, action='store_true',
                   help='Do anisotropic pressure scaling')
group.add_argument('--dt', dest='timestep', type=float,
                   metavar='FLOAT', help='''time step for integrator (outer
                   time-step for RESPA integrator) Default 1 fs''', default=1.0)
group.add_argument('--nrespa', dest='nrespa', type=int, metavar='INT',
                   default=1, help='''Number of inner time steps to run
                   (evaluating fast forces) for every outer timestep. Default is
                   1 (no RESPA). Best value to use for AMOEBA is at most 2. Only
                   AMOEBA is supported.''')
group.add_argument('--gamma_ln', dest='gamma_ln', type=float,
                   metavar='FLOAT', help='''collision frequency for Langevin
                   integrator. Default %(default)s ps-1''', default=0.0)
group.add_argument('--temp', dest='temp', type=float,
                   metavar='FLOAT', help='''target temperature for NVT
                   simulation. Default %(default)s K''', default=300.0)
group.add_argument('--nve', dest='nve', default=False, action='store_true',
                    help='Do constant energy simulation')
group.add_argument('--watch-for-errors', dest='hawkeye', default=False,
                   action='store_true',
                   help='''Watch energy every step and if energy becomes large
                   or NaN, print out each component to find where the energy is
                   blowing up. This may slow the simulation down considerably,
                   so it is primarily of use when debugging.''')
group.add_argument('--nframes', dest='frame_num', type=int,
                   metavar='INT', help='''Number of frames to be generated 
                   during the md run''', default=100)

opt = parser.parse_args()

print('Command line:\n\t%s' % ' '.join(sys.argv)); sys.stdout.flush()

print('Parsing XML file %s' % opt.xml); sys.stdout.flush()
with open(opt.xml, 'r') as f:
    system = mm.XmlSerializer.deserialize(f.read())

if opt.hawkeye:
    if opt.nrespa > 1:
        raise ValueError('Cannot use MTS integrator and watch for errors')
    groups_and_names = []
    for i, force in enumerate(system.getForces()):
#       if isinstance(force, mm.AmoebaMultipoleForce):
#           # Skip the multipole force, since it's so expensive
#           force.setForceGroup(20)
#           continue
        groups_and_names.append((type(force).__name__, i))
        force.setForceGroup(i)
    
    class ErrorDetectionReporter(app.StateDataReporter):
        def __init__(self, groups_and_names):
            self._reportInterval = 1
            self.groups_and_names = groups_and_names

        def describeNextReport(self, simulation):
            """Get information about the next report this object will generate.
    
            Parameters
            ----------
            simulation : Simulation
                The Simulation to generate a report for
    
            Returns
            -------
            tuple
                A five element tuple. The first element is the number of steps
                until the next report. The remaining elements specify whether
                that report will require positions, velocities, forces, and
                energies respectively.
            """
            return 1, False, False, False, True # only need energies

        def report(self, simulation, state):
            ene = state.getPotentialEnergy().value_in_unit(u.kilojoules_per_mole)
            if not math.isnan(ene) and ene < 1e5:
                return
            print('%30s %.6f kcal/mol' % ('Total Energy', ene))
            for name, i in self.groups_and_names:
                ene = simulation.context.getState(getEnergy=True, groups=1<<i).getPotentialEnergy().value_in_unit(u.kilocalories_per_mole)
                print('%30s %.6f kcal/mol' % (name, ene), file=sys.stderr)
            sys.exit('Bad energy! Look at components')
    
# Remove the andersen thermostat that might exist (e.g. from Tinker)
for i in range(system.getNumForces()):
    if isinstance(system.getForce(i), mm.AndersenThermostat):
        print('Deleting the Andersen thermostat'); sys.stdout.flush()
        system.removeForce(i)
        break
else:
    print('No Andersen thermostat found in the system.xml file'); sys.stdout.flush()

print('Parsing PDB file %s' % opt.pdb); sys.stdout.flush()
pdb =  app.PDBFile(opt.pdb)

if opt.ntp:
    if opt.aniso:
        print('Using anisotropic barostat'); sys.stdout.flush()
        baro = mm.MonteCarloAnisotropicBarostat(1*u.bar, opt.temp*u.kelvin)
    else:
        print('Using isotropic barostat'); sys.stdout.flush()
        baro = mm.MonteCarloBarostat(1*u.bar, opt.temp*u.kelvin)
    system.addForce(baro)

if opt.gamma_ln == 0.0 and not opt.nve:
    print('Adding Anderson thermostat at %sK, 0.1 psec^-1' % opt.temp); sys.stdout.flush()
    thermo = mm.AndersenThermostat(opt.temp*u.kelvin, 0.1/u.picosecond)
    system.addForce(thermo)

# Create the simulation
if opt.gamma_ln > 0.0:
    if opt.nrespa > 1:
        raise ValueError('nrespa and Langevin dynamics are incompatible')
    integrator = mm.LangevinIntegrator(opt.temp*u.kelvin,
                 opt.gamma_ln/u.picosecond, opt.timestep*u.femtoseconds)
    print('Langevin: %8.2fK, %8.2f ps-1, %8.2f fs' %
       (opt.temp, opt.gamma_ln, opt.timestep) ); sys.stdout.flush()
elif opt.nrespa > 1:
    slow = (mm.AmoebaMultipoleForce, mm.AmoebaVdwForce,
            mm.AmoebaGeneralizedKirkwoodForce, mm.AmoebaWcaDispersionForce)
    found_slow = False
    for force in system.getForces():
        if isinstance(force, slow):
            found_slow = True
            force.setForceGroup(1)
        else:
            force.setForceGroup(0)
    if not found_slow:
        raise ValueError('No slow AMOEBA forces found for MTS integrator!')
    # The list given to MTSIntegrator defining the time steps and force
    # decompositions is a list of 2-element tuples where the first element is
    # a force group and the second element is how many times to evaluate that
    # force group each "outer" time-step. So [(0, opt.nrespa), (1, 1)] means
    # force group 0 is executed nrespa times each time step and force group 1 is
    # executed only once. The slow forces are defined above as force group 1 and
    # all others as force group 0.
    integrator = mm.MTSIntegrator(opt.timestep*u.femtoseconds,
                                  [(0, opt.nrespa), (1, 1)])
    print('RESPA MTS Integrator: %8.2f fs outer time-step with %d inner steps' %
          (opt.timestep, opt.nrespa))
else:
    integrator = mm.VerletIntegrator(opt.timestep*u.femtoseconds)
    print('Verlet: %8.2f fs' % opt.timestep )

sim = app.Simulation(pdb.topology, system, integrator,
                     mm.Platform.getPlatformByName('CUDA'),
                     dict(CudaPrecision='mixed') )

if opt.hawkeye:
    # Watch every step... slow!
    sim.reporters.append(ErrorDetectionReporter(groups_and_names))

sim.reporters.append(
        app.StateDataReporter(opt.log, reportInterval=opt.interval,
                         step=True, volume=True, density=True, temperature=True, speed=True, separator='\t')
)

if opt.checkpoint is not None:
    sim.reporters.append(
            app.CheckpointReporter(opt.checkpoint, opt.interval*10)
    )

if opt.state is not None:
    print('Setting coordinates and velocities from restart file %s' %
        opt.state); sys.stdout.flush()

    sim.loadState(opt.state)
    previous_state = sim.context.getState(getVelocities=True, getPositions=True, getForces = True, getEnergy = True, getParameters = True, groups = -1)
    positions = previous_state.getPositions()
    velocities = previous_state.getVelocities()
    sim.context.setPositions(positions)
    sim.context.setVelocities(velocities)

else:
    print('Setting coordinates from PDB file %s' % opt.pdb); sys.stdout.flush()
    sim.context.setPositions(pdb.positions)
    sim.context.setVelocitiesToTemperature(opt.temp)

frameInterval = int(opt.num_steps / opt.frame_num)
print("Running the simulation for %i steps and producing %i frames %i steps separated." % (opt.num_steps, opt.frame_num, frameInterval)); sys.stdout.flush()

for i in range(0, opt.frame_num):
     sim.step(frameInterval)
     fileName = "state_number_"+str(i)+".xml"    
     print(fileName); sys.stdout.flush()
     with open(fileName, 'w') as f:
        f.write(
            mm.XmlSerializer.serialize(
                sim.context.getState(
                    getPositions=True, getVelocities=True, getForces=True,
                    enforcePeriodicBox=system.usesPeriodicBoundaryConditions(),
                    getEnergy=True
                )
            )
    )

# The last step may not have resulted in a restart file being written. Force it
# here
state = sim.context.getState(getPositions=True, getVelocities=True,
        getEnergy=True, getForces=True,
        enforcePeriodicBox=system.usesPeriodicBoundaryConditions())

print('Finished')


