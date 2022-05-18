#!/usr/bin/env python
from __future__ import division, print_function

from argparse import ArgumentParser

from simtk.openmm.app import *
from simtk.openmm import *
import simtk.openmm as mm
import simtk.unit as u
from sys import stdout, argv
from mdtraj.reporters import DCDReporter

parser = ArgumentParser()
group = parser.add_argument_group('Input/Output file options')
group.add_argument('-p', '--pdb', metavar='<PDB FILE>', required=True,
                   dest='pdb', help='PDB file with target system')
parser.add_argument('-x', '--xml', dest='output', metavar='FILE',
                    default='minimized.xml', help='''Output XML with minimized
                    coordinates. Default is %(default)s''')
group.add_argument('-v', '--vdw-cutoff', metavar='FLOAT', default=9,
                   type=float, help='''Cutoff in Angstroms to use for van der
                   Waals interactions. This is only applicable to AMOEBA force
                   fields and will be ignored otherwise. Default is
                   %(default)s''', dest='vdwcut')
group.add_argument('-e', '--epsilon', metavar='FLOAT', default=1e-5, dest='eps',
                   help='''Convergence criteria for mutually induced polarizable
                   dipoles in the AMOEBA force field. Default is %(default)g.
                   This option is ignored for fixed-charge FFs''', type=float)

opt = parser.parse_args()

# Parse the PDB file
print('Parsing the PDB file [%s]...' % opt.pdb)
pdb = app.PDBFile(opt.pdb)

forcefield = app.ForceField('amoeba2013.xml')
print('Creating system...')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=0.9*u.nanometer, constraints=None, rigidWater=False, ewaldErrorTolerance=0.0005)
for force in system.getForces():
    if isinstance(force, mm.AmoebaVdwForce):
        print('Adjusting the vdW cutoff to %g Angstroms...' % opt.vdwcut)
        force.setCutoff(opt.vdwcut*u.angstroms)
    elif isinstance(force, mm.AmoebaMultipoleForce):
        print('Setting the induced dipole convergence criteria to %g' % opt.eps)
        force.setMutualInducedTargetEpsilon(opt.eps)
        
print('Creating the Simulation...')
simulation = app.Simulation(pdb.topology, system, mm.VerletIntegrator(0.001),
                     platform=mm.Platform.getPlatformByName('CUDA'),
                     platformProperties=dict(CudaPrecision='mixed'))

simulation.context.setPositions(pdb.positions)

e=simulation.context.getState(getEnergy=True).getPotentialEnergy()
print(' Initial energy = %10.4f kcal/mol' % e.value_in_unit(u.kilocalories_per_mole))
simulation.minimizeEnergy()
e=simulation.context.getState(getEnergy=True).getPotentialEnergy()
print('   Final energy = %10.4f kcal/mol' % e.value_in_unit(u.kilocalories_per_mole))
# Now write a serialized state that has coordinates
print('Finished. Writing serialized XML restart file...')
with open(opt.output, 'w') as f:
    f.write(
            mm.XmlSerializer.serialize(
                simulation.context.getState(
                    getPositions=True, getVelocities=True, getForces=True,
                    enforcePeriodicBox=system.usesPeriodicBoundaryConditions(),
                    getEnergy=True
                )
            )
    )

print('Done!')
