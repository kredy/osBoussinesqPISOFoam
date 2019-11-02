/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         |
   \\    /   O peration     | Copyright (C) 2019 Krishna Reddy Maryada
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this software.  If not, see <http://www.gnu.org/licenses/>.

Application
    osBoussinesqPISOFoam

Description
    Transient solver for incompressible, buoyant flow, using the PISO 
    algorithm. Boussinesq approximation in its original formulation
    is used. The source term is added as a flux at the faces instead of the 
    cell-centres.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"


int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);
  
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // Start of time loop
    Info<< "Starting Time Loop\n" << endl;
    
/*---------------------------------------------------------------------------*/
   
    while (runTime.loop())
    {
        Info<< "Time = "<< runTime.timeName() << nl << endl;

        #include "CourantNo.H"
	
        #include "UEqn.H"
        
        while (piso.correct())
        {
            #include "pEqn.H"
        }

        #include "TEqn.H"
	
	    #include "continuityErrs.H"
	
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

/*---------------------------------------------------------------------------*/
