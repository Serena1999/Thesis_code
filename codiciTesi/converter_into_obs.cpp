/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                  converter_into_obs.cpp:                 ****
****  given a column of the file, I compute an observable and ****
****                put it into an output file.               ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

#include "../library.h"

//-----------------------------------------------------------------
//OBSERVABLE TO COMPUTE:

double obs(double x) {
	return x * x;
}

//-----------------------------------------------------------------
//MAIN:

int main() {

	int skip_in_lines = 1;

	int col_index = 5;//CAMBIA SECONDO NECESSITà: 5 per re, 6 oer im
	string directory_in_file = "19_05_2025/build_good_mpi1500_16_out/";//CAMBIA SECONDO NECESSITà
	string name_in_file = "ferm_obs2239384639.txt";//CAMBIA SECONDO NECESSITà
	string obs_name = "reff2";//CAMBIA SECONDO NECESSITà: usa (re/im)P2 per ((re/im)loop di Pol)^2 e (re/im)ff2 per ((re/im)cond

	string path_in_file = directory_in_file + name_in_file;
	string directory_out_file = "results/";
	string path_out_file = directory_out_file + obs_name + "_" + name_in_file;
	
	string line;
	
	ifstream input_file;
	input_file.open(path_in_file);
	if (!input_file) {
		cout << "Error opening file input file" << endl;
		return 1;
	}

	ofstream output_file;
	output_file.open(path_out_file);
	if (!output_file) {
		cout << "Error opening file output file" << endl;
		return 1;
	}

	output_file << setprecision(numeric_limits<double>::max_digits10);

	for (int ii = 0; ii < skip_in_lines; ii++) {
		if (!getline(input_file, line)) {
			cerr << "Error: there are less than " << skip_in_lines << " lines in the file: " << path_in_file << endl;
			return 1;
		}
	}

	while (getline(input_file, line)){
		double conf_id, value;
		istringstream iss(line);
		iss >> conf_id;
		for (int ii = 2; ii < col_index; ii++) {
			double discard;
			iss >> discard;
		}
		if (iss >> value){
			output_file << conf_id << "\t" << obs(value) << endl;
		}
		else {
			cerr << "Skipped line (badly formatted) from" << path_in_file << ": " << line << endl;
			return 1;
		}
	}

	output_file.close();
	input_file.close();

	return 0;
}