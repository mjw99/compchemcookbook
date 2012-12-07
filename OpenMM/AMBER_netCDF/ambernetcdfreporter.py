from ambernetcdffile import AmberNetCDFFile

class AmberNetCDFReporter(object):
    def __init__(self, filename, reportInterval):
        """Create a AmberNetCDFReporter.

        Parameters:
         - file (string) The file to write to
         - reportInterval (int) The interval (in time steps) at which to write frames
        """
        self._reportInterval = reportInterval
        self._netcdf = None
        self._filename = filename

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters:
         - simulation (Simulation) The Simulation to generate a report for
        Returns: A five element tuple.  The first element is the number of steps until the
        next report.  The remaining elements specify whether that report will require
        positions, velocities, forces, and energies respectively.
        """
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, True, False, False, False)

    def report(self, simulation, state):
        """Generate a report.

        Parameters:
         - simulation (Simulation) The Simulation to generate a report for
         - state (State) The current state of the simulation
        """
        if self._netcdf is None:
            self._netcdf = AmberNetCDFFile(self._filename, simulation.topology)

        self._netcdf.writeModel(state.getPositions(), state.getTime() )

    def __del__(self):
        self._netcdf.close()
