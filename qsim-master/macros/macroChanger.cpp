#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>


using namespace std;

int main() {
float E=2;

  for(int a=12507; a<512507; a=a+10000){
	
		stringstream stream; // Declares a stream to be saved as string
		stream << "runexample_"<< a << ".mac"; // Saves text into a string
	
		string filename = stream.str(); // declares previous string to be called filename
		
		
		ofstream myfile(filename.c_str()); // opens file which has the name of the string saved in filename
		
		
		if (myfile.is_open()){
			myfile << "/qsim/filename /home/sluddani/showerData/qsim_seed_" << a << ".root" <<  endl << endl;

			myfile << "/qsim/seed " << a << " #change number to change random number seed" << endl << endl << endl;


			myfile << "# Set run mode" <<  endl;
			myfile << "/qsim/fDetMode 3 " << endl;
			myfile << "/qsim/fStandMode 0" << endl;
			myfile << "/qsim/fSourceMode 1" << endl;
			myfile << "/qsim/fQMode 2" << endl << endl;

			myfile << "# Set detector properties" <<  endl; 
			myfile << "/qsim/fQuartzPolish 0.981" <<  endl; 
			myfile << "/qsim/fDetAngle 0 deg " <<  endl; 
			myfile << "/qsim/fDetPosX 0 cm" << endl;
			myfile << "/qsim/fDetPosY 0 cm" << endl << endl;

			myfile << "# set beam to \"raster\" " <<  endl; 
			myfile << " /qsim/xmin -5.25 cm" << endl;
			myfile << " /qsim/xmax 4.65 cm" << endl;
			myfile << " /qsim/ymin -12.3 cm" << endl;
			myfile << " /qsim/ymax 12.3 cm" << endl << endl;

			myfile << "/run/initialize" << endl << endl;
			
			myfile << "/qsim/emin   " << E << " GeV #beam energies" << endl;
			myfile << "/qsim/emax   " << E << " GeV" << endl << endl;
			E = E + 0.14;
			myfile << "/gun/particle e- #type of particle is an electron" << endl << endl;

			myfile << "/run/beamOn 1000" << endl << endl;

			myfile.close(); // closes the file
			
		
		}
		
		else cout << "unable to open file";
	
		
	
	}
	
	return 0;
	
}


