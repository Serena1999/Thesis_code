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
	return x;
}

//-----------------------------------------------------------------
//MAIN:

int main() {

	int skip_in_lines = 1;

	//PARAMETRI MODIFICABILI:-------------------------------------------------------------------------------------------------------------------
	
	//CARTELLE PER 800 MEV;
	vector <string> name_directories = {//CAMBIA SECONDO NECESSITà
		"11_05_2025/build_good_017_out/",
		"11_05_2025/build_good_018_out/",
		"11_05_2025/build_good_019_out/",
		"11_05_2025/build_good_020_out/"
	};
	
	//OSSERVABILI DI GAUGE A 800 MEV:
	vector <string> name_files = {//CAMBIA SECONDO NECESSITà
		"gauge_obs3314513751.txt",
		"gauge_obs1615200272.txt",
		"gauge_obs2602939787.txt", 
		"gauge_obs3624614669.txt"
	};

	int col_index = 6;//CAMBIA SECONDO NECESSITà: 5 per re, 6 oer im (numerazione a partire da 1)
	string obs_name = "imP";//CAMBIA SECONDO NECESSITà: usa (re/im)P2 per ((re/im)loop di Pol)^2 e (re/im)ff2 per ((re/im)cond

	//reff->imff->reP->imP


	//PARAMETRI NON MODIFICABILI:-------------------------------------------------------------------------------------------------------------------
	string path_in_file;
	string directory_out_file;
	string path_out_file;
	string line;


	for (int ii = 0; ii < name_files.size(); ii++) {
		path_in_file = name_directories[ii] + name_files[ii];
		directory_out_file = "results/";
		path_out_file = directory_out_file + obs_name + "_" + name_files[ii];

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

		output_file << fixed << setprecision(numeric_limits<double>::max_digits10);

		for (int ii = 0; ii < skip_in_lines; ii++) {
			if (!getline(input_file, line)) {
				cerr << "Error: there are less than " << skip_in_lines << " lines in the file: " << path_in_file << endl;
				return 1;
			}
		}

		while (getline(input_file, line)) {
			double conf_id, value;
			istringstream iss(line);
			iss >> conf_id;
			for (int ii = 2; ii < col_index; ii++) {
				double discard;
				iss >> discard;
			}
			if (iss >> value) {
				output_file << conf_id << "\t" << obs(value) << endl;
			}
			else {
				cerr << "Skipped line (badly formatted) from" << path_in_file << ": " << line << endl;
				return 1;
			}
		}

		output_file.close();
		input_file.close();
	}

	return 0;
}


//READ TO USE:

/*
//CARTELLE PER 1500 MEV;
	vector <string> name_directories = {//CAMBIA SECONDO NECESSITà
		"19_05_2025/build_good_mpi1500_1_out/",
		"19_05_2025/build_good_mpi1500_2_out/",
		"19_05_2025/build_good_mpi1500_3_out/",
		"19_05_2025/build_good_mpi1500_4_out/",
		"19_05_2025/build_good_mpi1500_5_out/",
		"19_05_2025/build_good_mpi1500_6_out/",
		"19_05_2025/build_good_mpi1500_7_out/",
		"19_05_2025/build_good_mpi1500_8_out/",
		"19_05_2025/build_good_mpi1500_9_out/",
		"19_05_2025/build_good_mpi1500_10_out/",
		"19_05_2025/build_good_mpi1500_11_out/",
		"19_05_2025/build_good_mpi1500_12_out/",
		"19_05_2025/build_good_mpi1500_13_out/",
		"19_05_2025/build_good_mpi1500_14_out/",
		"19_05_2025/build_good_mpi1500_15_out/",
		"19_05_2025/build_good_mpi1500_16_out/"
	};

	//OSSERVABILI DI FERM A 1500 MEV:
	vector <string> name_files = {//CAMBIA SECONDO NECESSITà
		"ferm_obs1194794680.txt",
		"ferm_obs958284955.txt",
		"ferm_obs2174723772.txt",
		"ferm_obs2393834332.txt",
		"ferm_obs2698338921.txt",
		"ferm_obs2087301928.txt",
		"ferm_obs1439591720.txt",
		"ferm_obs3844636190.txt",
		"ferm_obs4003409327.txt",
		"ferm_obs1698962147.txt",
		"ferm_obs1287052570.txt",
		"ferm_obs230875197.txt",
		"ferm_obs802462068.txt",
		"ferm_obs14065219.txt",
		"ferm_obs4048334874.txt",
		"ferm_obs2239384639.txt"
	};
 
//OSSERVABILI DI GAUGE A 1500 MEV:
vector <string> name_files = {//CAMBIA SECONDO NECESSITà
	"gauge_obs1194794680.txt",
	"gauge_obs958284955.txt",
	"gauge_obs2174723772.txt",
	"gauge_obs2393834332.txt",
	"gauge_obs2698338921.txt",
	"gauge_obs2087301928.txt",
	"gauge_obs1439591720.txt",
	"gauge_obs3844636190.txt",
	"gauge_obs4003409327.txt",
	"gauge_obs1698962147.txt",
	"gauge_obs1287052570.txt",
	"gauge_obs230875197.txt",
	"gauge_obs802462068.txt",
	"gauge_obs14065219.txt",
	"gauge_obs4048334874.txt",
	"gauge_obs2239384639.txt",
};

	//CARTELLE PER 800 MEV;
	vector <string> name_directories = {//CAMBIA SECONDO NECESSITà
		"11_05_2025/build_good_004_out/",
		"11_05_2025/build_good_006_out/",
		"11_05_2025/build_good_001_out/",
		"11_05_2025/build_good_007_out/",
		"11_05_2025/build_good_008_out/",
		"11_05_2025/build_good_009_out/",
		"11_05_2025/build_good_002_out/",
		"11_05_2025/build_good_010_out/",
		"11_05_2025/build_good_011_out/",
		"11_05_2025/build_good_012_out/",
		"11_05_2025/build_good_013_out/",
		"11_05_2025/build_good_003_out/",
		"11_05_2025/build_good_014_out/",
		"11_05_2025/build_good_015_out/",
		"11_05_2025/build_good_016_out/",
		"11_05_2025/build_good_005_out/",
		"11_05_2025/build_good_017_out/",
		"11_05_2025/build_good_018_out/",
		"11_05_2025/build_good_019_out/",
		"11_05_2025/build_good_020_out/"
	};

	//OSSERVABILI DI GAUGE A 800 MEV:
	vector <string> name_files = {//CAMBIA SECONDO NECESSITà
		"gauge_obs4125601828.txt",
		"gauge_obs3379062322.txt",
		"gauge_obs3725763671.txt",
		"gauge_obs1709664641.txt",
		"gauge_obs3812052727.txt",
		"gauge_obs920250475.txt",
		"gauge_obs1375234971.txt",
		"gauge_obs637002155.txt",
		"gauge_obs1229215857.txt",
		"gauge_obs3483143181.txt",
		"gauge_obs3946243793.txt",
		"gauge_obs443817478.txt",
		"gauge_obs1917658109.txt",
		"gauge_obs4237028932.txt",
		"gauge_obs2847033391.txt",
		"gauge_obs133777317.txt",
		"gauge_obs3314513751.txt",
		"gauge_obs1615200272.txt",
		"gauge_obs2602939787.txt",
		"gauge_obs3624614669.txt"

	};


	//OSSERVABILI DI FERM A 800 MEV:
	vector <string> name_files = {//CAMBIA SECONDO NECESSITà
		"ferm_obs4125601828.txt",
		"ferm_obs3379062322.txt",
		"ferm_obs3725763671.txt",
		"ferm_obs1709664641.txt",
		"ferm_obs3812052727.txt",
		"ferm_obs920250475.txt",
		"ferm_obs1375234971.txt",
		"ferm_obs637002155.txt",
		"ferm_obs1229215857.txt",
		"ferm_obs3483143181.txt",
		"ferm_obs3946243793.txt",
		"ferm_obs443817478.txt",
		"ferm_obs1917658109.txt",
		"ferm_obs4237028932.txt",
		"ferm_obs2847033391.txt",
		"ferm_obs133777317.txt",
		"ferm_obs3314513751.txt",
		"ferm_obs1615200272.txt",
		"ferm_obs2602939787.txt",
		"ferm_obs3624614669.txt"
	};



*/