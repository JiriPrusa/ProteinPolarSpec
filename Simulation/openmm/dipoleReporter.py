"""
statedatareporter.py: Outputs data about a simulation

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012-2013 Stanford University and the Authors.
Authors: Peter Eastman
Contributors: Robert McGibbon

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from __future__ import absolute_import
from __future__ import print_function
__author__ = "Peter Eastman"
__version__ = "1.0"

try:
    import bz2
    have_bz2 = True
except: have_bz2 = False

try:
    import gzip
    have_gzip = True
except: have_gzip = False

import simtk.openmm as mm
from simtk.openmm.vec3 import Vec3
import simtk.unit as unit
import math
import time

class DipoleReporter(object):
    """DipoleReporter outputs information about a dipole moments to a file.

    To use it, create a DipoleReporter, then add it to the Simulation's list of reporters.  The set of
    data to write is configurable using boolean flags passed to the constructor.  By default the data is
    written in comma-separated-value (CSV) format, but you can specify a different separator to use.
    """

    def __init__(self, file, reportInterval, step=False, time=False, totalDipole = False, inducedDipoles = False, makeSum = False,  separator=',', indexList=None):
        """Create a StateDataReporter.

        Parameters
        ----------
        file : string or file
            The file to write to, specified as a file name or file object
        reportInterval : int
            The interval (in time steps) at which to write frames
        step : bool=False
            Whether to write the current step index to the file
        time : bool=False
            Whether to write the current time to the file
        permanentDipoles : bool=True
        	Whether to write permanent dipoles
        inducedDipoles : bool=True
        	Whether to write induced dipoles
        makeSum : bool=False
		Whether to output the sum of induced dipoles per frame. If false induced dipoles are written for each atom.
        volume : bool=False
            Whether to write the periodic box volume to the file
        density : bool=False
            Whether to write the system density to the file
        separator : string=','
            The separator to use between columns in the file
        """
        self._reportInterval = reportInterval
        self._openedFile = isinstance(file, str)
        
        if self._openedFile:
            # Detect the desired compression scheme from the filename extension
            # and open all files unbuffered
            if file.endswith('.gz'):
                if not have_gzip:
                    raise RuntimeError("Cannot write .gz file because Python could not import gzip library")
                self._out = gzip.GzipFile(fileobj=open(file, 'wb', 0))
            elif file.endswith('.bz2'):
                if not have_bz2:
                    raise RuntimeError("Cannot write .bz2 file because Python could not import bz2 library")
                self._out = bz2.BZ2File(file, 'w', 0)
            else:
                self._out = open(file, 'w')
        else:
            self._out = file
        self._step = step
        self._time = time
        self._reportIndexes = indexList
        self._totalDipole = totalDipole
        self._inducedDipoles = inducedDipoles
        self._makeSum = makeSum
        self._separator = separator
        self._hasInitialized = False
        self._needsPositions = False
        self._needsVelocities = False
        self._needsForces = False
        self._needEnergy = False

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
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, self._needsPositions, self._needsVelocities, self._needsForces, self._needEnergy)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        sep=self._separator
        if not self._hasInitialized:
            self._initializeConstants(simulation)
            headers = self._constructHeaders()
            print('#"%s"' % ('\t\t\t').join(headers), file=self._out)
            try:
                self._out.flush()
            except AttributeError:
                pass
            self._initialClockTime = time.time()
            self._initialSimulationTime = state.getTime()
            self._initialSteps = simulation.currentStep
            self._hasInitialized = True

        # Check for errors.
        self._checkForErrors(simulation, state)

        # Query for the values
        values = self._constructReportValues(simulation, state, simulation.context.getSystem(), simulation.context)

        # Write the values.
        
        #print("frame: " + str(self._initialSteps), file =self._out)
        
        #print(self._separator.join(', '.join(map(str,v)) for v in values), file=self._out)
        if self._step or self._time:
            print(values.pop(0), end=sep, file=self._out)    

        print(sep.join(sep.join(map(str,v)) for v in values), file=self._out)

        try:
            self._out.flush()
        except AttributeError:
            pass

    def _constructReportValues(self, simulation, state, system, context):
        """Query the simulation for the current state of our observables of interest.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation

        Returns
        -------
        A list of values summarizing the current state of
        the simulation, to be printed or saved. Each element in the list
        corresponds to one of the columns in the resulting CSV file.
        """
        debye = 0.02081943
        values = []
        permSum = Vec3(0,0,0)
        indSum = Vec3(0,0,0)
        box = state.getPeriodicBoxVectors()
        volume = box[0][0]*box[1][1]*box[2][2]
        clockTime = time.time()
        separator = self._separator
        makeSum = self._makeSum
        totalDipole = self._totalDipole
        inducedDipoles = self._inducedDipoles 
        if self._step:
            values.append(simulation.currentStep)
        if self._time:
            values.append(state.getTime().value_in_unit(unit.picosecond))
        #if self._volume:
        #    values.append(volume.value_in_unit(unit.nanometer**3))
        #if self._density:
        #    values.append((self._totalMass/volume).value_in_unit(unit.gram/unit.item/unit.milliliter))
        if inducedDipoles or totalDipole:
            # check if AmoebaMultipoleForce exists since charges needed
            # if it has not been created, raise an error
            existing = [system.getForce(i) for i in range(system.getNumForces())]
            amoebaMultipoleForceList = [f for f in existing if type(f) == mm.AmoebaMultipoleForce]
            if (len(amoebaMultipoleForceList) > 0):
                amoebaMultipoleForce = amoebaMultipoleForceList[0]
            else:
                raise RuntimeError("Cannot get induced dipoles! Non-polarizable FF used.")
        
            if totalDipole:
                permanent_all = amoebaMultipoleForce.getTotalDipoles(context)
                if self._reportIndexes is not None:
                    # make a subset based on selected chain
                    permanent = [permanent_all[i] for i in self._reportIndexes]
                else:
                    permanent = permanent_all
                # convert to debye
                per_asDebye = [x.__div__(debye) for x in permanent]
        
            if inducedDipoles:
                induced_all = amoebaMultipoleForce.getInducedDipoles(context)
                if self._reportIndexes is not None:
                    # make a subset based on selected chain
                    induced = [induced_all[i] for i in self._reportIndexes]
                else:
                    induced = induced_all
                # convert to debye
                ind_asDebye = [x.__div__(debye) for x in induced]   
 

            if inducedDipoles and totalDipole:
                for perdip, indip in zip(per_asDebye, ind_asDebye):
                    if makeSum:
                        indSum=indSum.__add__(indip)
                        permSum=permSum.__add__(perdip) 
                    else:
                        values.append(perdip[0:3] + indip[0:3])
                if makeSum:
                    values.append(permSum[0:3] + indSum[0:3])

            elif totalDipole:
                for dip in per_asDebye:
                    if makeSum:
                        permSum=permSum.__add__(dip)
                    else:
                        values.append(dip[0:3])
                if makeSum:
                    values.append(permSum)
            else:        
                for dip in ind_asDebye:
                    if makeSum:
                        indSum=indSum.__add__(dip)
                    else:
                        values.append(dip[0:3])
                if makeSum:
                    values.append(indSum)
        return values

    def _initializeConstants(self, simulation):
        """Initialize a set of constants required for the reports

        Parameters
        - simulation (Simulation) The simulation to generate a report for
        """
        system = simulation.system

        #if self._density:
        #    if self._totalMass is None:
        #        # Compute the total system mass.
        #        self._totalMass = 0*unit.dalton
        #        for i in range(system.getNumParticles()):
        #            self._totalMass += system.getParticleMass(i)
        #    elif not unit.is_quantity(self._totalMass):
        #        self._totalMass = self._totalMass*unit.dalton

    def _constructHeaders(self):
        """Construct the headers for the CSV output

        Returns: a list of strings giving the title of each observable being reported on.
        """
        sep = self._separator
        headers = []
        if self._totalDipole:
            headers.append('x'+sep+'y'+sep+'z')
        if self._inducedDipoles:
            headers.append('x'+sep+'y'+sep+'z')
        #if self._volume:
        #    headers.append('Box Volume (nm^3)')
        #if self._density:
        #    headers.append('Density (g/mL)')
        return headers
    

    def _checkForErrors(self, simulation, state):
        """Check for errors in the current state of the simulation

        Parameters
         - simulation (Simulation) The Simulation to generate a report for
         - state (State) The current state of the simulation
        """
        pass

    def __del__(self):
        if self._openedFile:
            self._out.close()


