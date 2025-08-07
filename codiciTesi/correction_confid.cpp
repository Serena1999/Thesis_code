/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                   correction_confid.cpp:                 ****
****   to correct files from confid erroneously duplicated.   ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

//-----------------------------------------------------------------
//HEADERS && LIBRARY:

#include "../library.h"

//-----------------------------------------------------------------
//FUNCTION DECLARATIONS:

void read_file_list_DIRECTORIES_THERM(//SCRIVI SOTTO COME IMPLEMENTARLA, PRENDI COME MINIMO DI TERMALIZZAZIONE IL MASSIMO FRA GLI N_SKIP, COSì SEI SICURA CHE VA BENE;
	const string& name_file_list,
	const int skipLines_file_list,
	vector<string>& directories,
	vector<int>& n_skip,
	int step_sample_gauge,
	int step_sample_fermion
);

//-----------------------------------------------------------------
//MAIN:

int main() {

	string mpi = "800";
	string line;

	string list_file = "11_05_2025/file_list_therm_extended.txt";
	int skipLines_list_file = 1;
	int step_sample_gauge = 1;
	int step_sample_fermion = 10;

	vector <string> directories;//LEGGILI DA FILE
	vector <int> n_skip; //LEGGILO DA FILE;
	string name_input_file = "mon.dat";
	string name_output_file = "corrected_" + name_input_file;
	int skip_lines_input_file = 0;

	read_file_list_DIRECTORIES_THERM(
		list_file,
		skipLines_list_file,
		directories,
		n_skip,
		step_sample_gauge,
		step_sample_fermion
	);

	cout << "READ: GOOD" << endl;

	//MOMENTANEO:
	//directories.clear();
	//directories.push_back("19_05_2025/build_good_mpi1500_16_out/");
	//name_input_file = "ferm_obs2239384639.txt";
	//name_output_file = "corrected_" + name_input_file;

	for (int ii = 0; ii < directories.size(); ++ii) {

		cout << ii << "STARTED" << endl;

		ifstream input_file;
		input_file.open(directories[ii] + name_input_file);
		if (!input_file) {
			cerr << "Error opening input " << ii << "-th file." << endl;
			return 1;
		}

		ofstream output_file;
		output_file.open(directories[ii] + name_output_file);
		if (!output_file) {
			cerr << "Error opening first output file" << endl;
			return 1;
		}

		//MOMENTANEO:
		//output_file << "#conf_id    acc     plq     rect    ReP     ImP" << endl;
		//output_file << "#0.conf     1.icopy 02.Plaquette    03.Rectangle    04.Reff_light_s           05.Imff_light_s           06.ReN_light_s            07.ImN_light_s            08.ReMag_light_s           09.ImMag_light_s           10.ReChSuscConn_light_s  11.ImChSuscConn_light_s  12.ReQNSuscConn1_light_s 13.ImQNSuscConn1_light_s 14.ReQNSuscConn2_light_s  15.ImQNSuscConn2_light_s" << endl;

		//we don't consider the first skip_lines_input_file:
		for (int jj = 0; jj < skip_lines_input_file; ++jj) {
			if (!getline(input_file, line)) {
				cerr << "Error: there are less than " << skip_lines_input_file << " lines in the input " << ii << "-th file." << endl;
				return 1;
			}
		}

		vector <double> conf_id_stored;
		double conf_id = -999, conf_tmp = 999;
		bool flag_init = 0;

		//for each configurazion:
		double discard1, discard2, discard3, discard4, discard5, discard6, discard7, wrap_value;
		double n_accum = 0;
		int n_conf = 0;
		vector <double> value_rhok, value_rhok_rho1, value_rhok_norm, k_array;
		bool flag = 0, flag_output = 0;

		while (getline(input_file, line)) {
			istringstream iss(line);
			if (iss >> conf_id) {
				if (!flag_init) {
					conf_tmp = conf_id;
					flag_init = 1;
				}
				if (conf_id != conf_tmp) {
					flag_output = 0;
					for (int kk = 0; kk < conf_id_stored.size(); ++kk) {
						if(conf_id == conf_id_stored[kk]) {
							flag_output = 1;
							break;
						}
					}
					if (!flag_output) {
						conf_id_stored.push_back(conf_id);
					}
					conf_tmp = conf_id;
				}
				if (!flag_output) {
					output_file << line << endl;
				}
			}
			else {
				cerr << "Skipped line (badly formatted) from" << (ii + 1) << "-th input file: " << line << endl;
			}
		}

		output_file.close();
		input_file.close();

	}

	return 0;
}

//-----------------------------------------------------------------
//FUNCTION DEFINITIONS:

void read_file_list_DIRECTORIES_THERM(
	const string& name_file_list,
	const int skipLines_file_list,
	vector<string>& directories,
	vector<int>& n_skip,
	int step_sample_gauge,
	int step_sample_fermion
) {
	string line;
	ifstream file_list;
	file_list.open(name_file_list);
	if (!file_list) {
		cout << "Error opening file list" << endl;
		return;
	}

	for (int i = 0; i < skipLines_file_list; i++) {
		if (!getline(file_list, line)) {
			cerr << "Error: there are less than " << skipLines_file_list << " lines in the file list." << endl;
			return;
		}
	}

	while (getline(file_list, line)) {
		if (line.empty()) {
			cerr << "Skipped blank/whitespace-only line in file: " << name_file_list << endl;
			continue;
		}
		istringstream iss(line);
		string dir, gauge, ferm;
		int max_n_therm;
		vector <int> n_therm(4, 0);
		if (iss >> dir >> gauge >> n_therm[0] >> n_therm[1] >> ferm >> n_therm[2] >> n_therm[3]) {
			directories.push_back(dir);
			n_therm[0] /= step_sample_gauge;
			n_therm[1] /= step_sample_gauge;
			n_therm[2] /= step_sample_fermion;
			n_therm[3] /= step_sample_fermion;
			max_n_therm = n_therm[0];
			for (int ii = 1; ii < n_therm.size(); ++ii) {
				if (max_n_therm < n_therm[ii]) {
					max_n_therm = n_therm[ii];
				}
			}
			n_skip.push_back(max_n_therm);
		}
		else {
			cerr << "Poorly formatted line: " << line << endl;
		}
	}

	file_list.close();
}