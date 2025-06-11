/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****              susceptibility_with_errors.cpp:              ****
****    this script takes the "data_vale"-style measurement   ****
****  files for each temperature and returns a file with the  ****
****  temperature along with the mean and standard deviation  ****
****        of the suscebility of the wanted variable.        ****
****    When selected, blocking is applied with block sizes   ****
****             selectable from an external file.            ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

#include "../library.h"
#include "../root_include.h"

const double hbar_c = 197.3269804; //MeV * fm

const bool debug_mode = 0;
const string tipology = "reff";//YOU CAN CHOOSE BETWEEN reP, imP, reff, imff;

//-----------------------------------------------------------------
//DECLARATIONS:

double obs_function(double prefactor, double x_mean, double x2_mean);

void read_file_LPC(
	const string& name_file_lpc,
	const int Nt,
	const int skipLines_file_lpc,
	const double mpi,
	vector<double>& aml,
	vector<double>& beta,
	vector<double>& afm,
	vector<double>& temp
);

void read_file_list(
	const string& name_file_list,
	const int skipLines_file_list,
	const int skipLines,
	vector<string>& directories,
	vector<string>& files,
	vector<int>& n_skip,
	vector<int>& dim_block,
	const string tipology
);

void susceptibility_with_errors(
	const int n_steps,
	const string& input_path,
	const string& output_path,
	const bool append_mode,
	const double temp,
	const double prefactor,
	int n_skip,
	int dim_block,
	int n_skip_file
);

template <class T> int blocking_sample(
	vector <T>& original_draws,
	vector <T>& blocked_draws,
	const int dim_block
);

//-----------------------------------------------------------------
//MAIN:

int main() {
	int n_steps = 1e4;//TO CHOOSE: NUMBER OF BOOSTRAP STEPS 
	const int Nt = 8; //BE CAREFUL TO CHOOSE IT WELL;
	const int Ns = 32; //BE CAREFUL TO CHOOSE IT WELL;
	const int Vs = Ns * Ns * Ns;
	int skipLines_file_lpc = 2, skipLines_file_list = 1, skipLines = 0;
	double mpi = 1500; //MeV //BE CAREFUL TO CHOOSE IT WELL;
	bool bool_startFile = 1;//BE CAREFUL TO CHOOSE IT WELL;
	vector<int> append_mode(20, 1); //80 entries with value = 1 (same size of beta); 
	ostringstream mpi_stream;//TO INTRODUCE ALSO IN NUMERICAL METHODS CODE: IT IS USEFUL;
	mpi_stream << std::fixed << std::setprecision(1) << mpi; //set to 1 decimal place
	string mpi_string = mpi_stream.str(); // conversion into string
	string name_output_file = "results/" + mpi_string + "_" + tipology + "_results.txt";
	string name_file_lpc = "19_05_2025/data_value/lcp_data_value.txt";
	string name_file_list = "19_05_2025/data_value/file_list_therm.txt";

	double temp_value;
	vector<int> n_skip, n_skip2, dim_block, dim_block2;
	vector<double> aml, beta, afm, temp;//T = \hbar * c /(Nt * a[fm]) (1.60), Nt = 8; 
	vector<string> directories, directories2, files, files2;

	read_file_LPC(
		name_file_lpc,
		Nt,
		skipLines_file_lpc,
		mpi,
		aml,
		beta,
		afm,
		temp
	);

	//per i file di value;
	read_file_list(
		name_file_list,
		skipLines_file_list,
		skipLines,
		directories,
		files,
		n_skip,
		dim_block,
		tipology
	);

	if (bool_startFile) {

		ofstream output_file; //declaration of output file
		output_file.open(name_output_file);
		if (!output_file) {
			cerr << "Error opening output file" << endl;
			return 1;
		}

		if (tipology == "reP") {
			output_file << "# T \t Chi_reP \t err(Chi_reP)" << endl;
		}
		else if (tipology == "imP") {
			output_file << "# T \t Chi_imP \t err(Chi_imP)" << endl;
		}
		else if (tipology == "reff") {
			output_file << "# T \t Chi_reff \t err(Chi_reff)" << endl;
		}
		else if (tipology == "imff") {
			output_file << "# T \t Chi_imff \t err(Chi_imff)" << endl;
		}
		else {
			cerr << "PROBLEM OR IN FILE LIST OR IN TIPOLOGY";
			return 1;
		}
		output_file.close();
	}

	for (int ii = 0; ii < directories.size(); ii++) {
		susceptibility_with_errors(
			n_steps,
			directories[ii] + files[ii],
			name_output_file,
			append_mode[ii],
			temp[ii],
			Vs,
			n_skip[ii],
			dim_block[ii],
			skipLines
		);
		cout << "file n°" << ii << " DONE! T = " << temp[ii] << endl;
		cout << endl;
	}

	return 0;
}

//-----------------------------------------------------------------
//FUNCTION DEFINITION:

double obs_function(double prefactor, double x_mean, double x2_mean) {
	return prefactor * (x2_mean - x_mean * x_mean); // Vs*<x^2> - <x>^2
}

void read_file_LPC(
	const string& name_file_lpc,
	const int Nt,
	const int skipLines_file_lpc,
	const double mpi,
	vector<double>& aml,
	vector<double>& beta,
	vector<double>& afm,
	vector<double>& temp
) {
	string line;
	ifstream file_lpc;

	file_lpc.open(name_file_lpc);
	if (!file_lpc) {
		cerr << "Error opening file list" << endl;
		return;
	}

	for (int i = 0; i < skipLines_file_lpc; i++) {
		if (!getline(file_lpc, line)) {
			cerr << "Error: there are less than " << skipLines_file_lpc << " lines in the lpc file." << endl;
			return;
		}
	}

	while (getline(file_lpc, line)) {
		if (line.empty()) {
			cerr << "Skipped blank/whitespace-only line in file: " << name_file_lpc << endl;
			continue;
		}
		istringstream iss(line);
		double aml_value, beta_value, afm_value;
		string dir, gauge, ferm;
		if (iss >> aml_value >> beta_value >> afm_value) {
			aml.push_back(aml_value);
			beta.push_back(beta_value);
			afm.push_back(afm_value);
		}
		else {
			cerr << "Poorly formatted line: " << line << endl;
		}
	}

	for (int ii = 0; ii < (beta.size()); ii++) {
		temp.push_back(hbar_c / (Nt * afm[ii]));
	}

	file_lpc.close();
}

void read_file_list(
	const string& name_file_list,
	const int skipLines_file_list,
	const int skipLines,
	vector<string>& directories,
	vector<string>& files,
	vector<int>& n_skip,
	vector<int>& dim_block,
	const string tipology
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
		string dir, group_name, file_name;
		int n_therm_value, step_sample, dim_block_value;
		if (iss >> dir >> file_name >> group_name >> n_therm_value >> step_sample >> dim_block_value) {
			if (group_name == tipology) {
				directories.push_back(dir);
				files.push_back(file_name);
				n_skip.push_back(n_therm_value / step_sample);
				dim_block.push_back(dim_block_value);
			}
		}
		else {
			cerr << "Poorly formatted line: " << line << endl;
		}
	}

	file_list.close();
}

void susceptibility_with_errors(
	const int n_steps,
	const string& input_path,
	const string& output_path,
	const bool append_mode,
	const double temp,
	const double prefactor,
	int n_skip,
	int dim_block,
	int n_skip_file
) {

	int index;
	double value;
	double mean_chi = 0, var_chi = 0, delta_chi;
	double mean = 0, mean_tmp = 0, delta;
	double mean2 = 0, mean_tmp2 = 0, delta2;

	double mean_chi2 = 0;

	sample_gen sampler;
	vector <double> original_draws, blocked_draws;
	vector <double> original_draws2, blocked_draws2;

	ofstream output_file; //declaration of output file
	if (append_mode) {
		output_file.open(output_path, ios::app);
	}
	else
	{
		output_file.open(output_path);
	}
	if (!output_file) {
		cout << "Error opening output file" << endl;
		return;
	}

	read_1from2columns(1, n_skip + n_skip_file, original_draws, input_path);
	
	if (blocking_sample(original_draws, blocked_draws, dim_block)) {
		cout << "Problem in blocking of original_draws" << endl;
		output_file.close();
		return;
	}

	/*
		for (int ii = 0; ii < original_draws.size(); ii++) {
			value = original_draws[ii];
			original_draws2.push_back(value * value);
		}

		if(blocking_sample(original_draws2, blocked_draws2, dim_block)) {
			cout << "Problem in blocking of original_draws2" << endl;
			output_file.close();
			return;
		}
	*/

	uniform_int_distribution<> dist_int(0, blocked_draws.size() - 1);

	for (int ii = 0; ii < n_steps; ii++) {

		mean_tmp = 0;
		mean_tmp2 = 0;
		
		for (int jj = 0; jj < blocked_draws.size(); jj++) {
			
			index = dist_int(sampler.rng);
			//cout << index << "\t";
			value = blocked_draws[index];
			mean_tmp += value;//
			//delta = value - mean_tmp;
			//mean_tmp = mean_tmp + delta / (double) (jj + 1);


			/*
			
			value = blocked_draws2[index];
			mean_tmp2 += value;//

			*/
			
			mean_tmp2 += (value * value);//

			//delta2 = value - mean_tmp2;
			//mean_tmp2 = mean_tmp2 + delta2 / (double)(jj + 1);

		}

		mean_tmp /= (double)blocked_draws.size();
		mean_tmp2 /= (double)blocked_draws.size();

		//cout << endl;

		//delta = mean_tmp - mean;
		//mean = mean + delta / (double) (ii + 1);
		//delta2 = mean_tmp2 - mean2;
		//mean2 = mean2 + delta2 / (double) (ii + 1);

		mean = mean_tmp;//
		mean2 = mean_tmp2;

		value = obs_function(prefactor, mean, mean2);
		//delta_chi = value - mean_chi;
		//mean_chi = mean_chi + delta_chi / (double) (ii + 1);
		mean_chi += value;//
		mean_chi2 += value * value;
		//assert(var_chi >= 0);
		//cout << var_chi << endl;
		//var_chi = var_chi + delta * (value - mean_chi);
		//assert(isnan(var_chi) == false);
		//assert(var_chi < 0);
		//assert(var_chi >= 0);

	}

	mean_chi /= (double)n_steps;
	mean_chi2 /= (double)n_steps;
	var_chi = mean_chi2 - mean_chi * mean_chi;
	assert(var_chi > 0);
	var_chi *= (double)n_steps;
	var_chi /= (double) (n_steps - 1);
	assert(var_chi > 0);

	if (var_chi < 0) {
		cout << "[WARNING] var_chi problem at T = " << temp << ", var_chi = " << var_chi << endl;
		cout << "dim_block = " << dim_block << endl;
		cout << "len original = " << original_draws.size() << endl;
		cout << "len blocked = " << blocked_draws.size() << endl;
	}


	cout << tipology << ": <Chi> = " << mean_chi << " +- " << sqrt(var_chi) << endl;
	
	output_file << setprecision(numeric_limits<double>::max_digits10);
	output_file << temp << "\t" << mean_chi << "\t" << sqrt(var_chi) << endl;

	original_draws.clear();
	original_draws2.clear();
	blocked_draws.clear();
	blocked_draws2.clear();

	output_file.close();
}

template <class T> int blocking_sample(
	vector <T>& original_draws,
	vector <T>& blocked_draws,
	const int dim_block
) {

	int n_blocks = original_draws.size() / dim_block;//=number of blocks
	if (n_blocks < 2) {
		std::cerr << "Warning: not enough blocks (" << n_blocks << ") to estimate variance reliably. Returning NaN." << std::endl;
		return 1;
	}

	int n_max = n_blocks * dim_block;//N_max to consider to compute variance
	int index = 0;
	T mean_tmp, value;
	for (int jj = 0; jj < n_max; jj += dim_block) {
		index++;
		mean_tmp = 0;
		for (int kk = 1; kk <= dim_block; kk++) {
			//delta = original_draws[jj + kk - 1] - mean_tmp;
			//mean_tmp = mean_tmp + delta / (double) kk;
			value = original_draws[jj + kk - 1];
			mean_tmp += value;
		}
		mean_tmp /= (double)dim_block;
		blocked_draws.push_back(mean_tmp);
	}

	return 0;
}
