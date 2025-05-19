/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                        polyff.cpp:                       ****
****  this script takes the fermionic and bosonic measurement ****
****  files for each temperature and returns a file with the  ****
****  temperature along with the mean and standard deviation  ****
****                  for these observables.                  ****
****    When selected, blocking is applied with block sizes   ****
****             selectable from an external file.            ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

#include "../library.h"
#include "../root_include.h"

const double hbar_c = 197.3269804; //MeV * fm

//-----------------------------------------------------------------
//DECLARATIONS:

void read_file_LPC(
	const string& name_file_lpc,
	const int Nt,
	const int skipLines_file_lpc,
	const double mpi,
	vector<double>& aml,
	vector<double>& beta,
	vector<double>& afm,
	vector<double>& temp,
	vector<int>& dim_block_modP,
	vector<int>& dim_block_reP,
	vector<int>& dim_block_imP,
	vector<int>& dim_block_modff,
	vector<int>& dim_block_reff,
	vector<int>& dim_block_imff
);

void read_file_list(
	const string& name_file_list,
	const int skipLines_file_list,
	const int skipLines,
	vector<string>& directories,
	vector<string>& gauge_files,
	vector<string>& fermion_files,
	vector<int>& n_copy,
	vector<int>& n_skip_rep,
	vector<int>& n_skip_imp,
	vector<int>& n_skip_reff,
	vector<int>& n_skip_imff
);

void stats_thesis(
	const string& input_path,
	const string& output_path,
	const string& tipology, //fermion/gauge
	const bool append_mode,
	const double temp,
	int n_skip_re,
	int n_skip_im,
	int dim_block,
	int dim_block_re,
	int dim_block_im,
	int n_copy
);

//-----------------------------------------------------------------
//MAIN:

int main() {
	int Nt = 8; //BE CAREFUL TO CHOOSE IT WELL;
	int skipLines_file_lpc = 2, skipLines_file_list = 1, skipLines = 1;
	double mpi = 800; //MeV //BE CAREFUL TO CHOOSE IT WELL;
	bool bool_startFile_poly = 1, bool_startFile_ff = 1;//BE CAREFUL TO CHOOSE IT WELL;
	double temp_value;
	vector<int> append_mode_poly = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };//20 entries (same size of beta);
	vector<int> append_mode_ff = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 , 1, 1, 1, 1};//20 entries (same size of beta);
	vector<int> n_skip_rep, n_skip_imp, n_skip_reff, n_skip_imff, n_copy;
	vector<int> dim_block_modP, dim_block_reP, dim_block_imP, dim_block_modff, dim_block_reff, dim_block_imff;
	vector<double> aml, beta, afm, temp;//T = \hbar * c /(Nt * a[fm]) (1.60), Nt = 8; 
	vector<string> directories, gauge_files, fermion_files;
	ostringstream mpi_stream;//TO INTRODUCE ALSO IN NUMERICAL METHODS CODE: IT IS USEFUL;
	mpi_stream << std::fixed << std::setprecision(1) << mpi; //set to 1 decimal place
	string mpi_string = mpi_stream.str(); // conversion into string
	string name_output_file_poly = "results/" + mpi_string + "_poly_results.txt";
	string name_output_file_ff = "results/" +  mpi_string + "_ff_results.txt";
	string name_file_lpc = "11_05_2025/LCP_800MeV_dimblock_extended.txt";
	string name_file_list = "11_05_2025/file_list_therm_extended.txt";

	read_file_LPC(
		name_file_lpc,
		Nt,
		skipLines_file_lpc,
		mpi,
		aml,
		beta,
		afm,
		temp,
		dim_block_modP,
		dim_block_reP,
		dim_block_imP,
		dim_block_modff,
		dim_block_reff,
		dim_block_imff
	);

	read_file_list(
		name_file_list,
		skipLines_file_list,
		skipLines,
		directories,
		gauge_files,
		fermion_files,
		n_copy,
		n_skip_rep, 
		n_skip_imp,
		n_skip_reff,
		n_skip_imff
	);

	if (bool_startFile_poly) {
		ofstream output_file; //declaration of output file
		output_file.open(name_output_file_poly);
		if (!output_file) {
			cerr << "Error opening output file" << endl;
			return 1;
		}
		output_file << "# T \t |<P * P^dag>| \t err(|<P* P^dag>|) \t Re{P} \t err(Re{P}) \t Im{P} \t err(Im{P})" << endl;
		output_file.close();
	}
	
	if (bool_startFile_ff) {
		ofstream output_file; //declaration of output file
		output_file.open(name_output_file_ff);
		if (!output_file) {
			cerr << "Error opening output file" << endl;
			return 1;
		}
		output_file << "# T \t |<ff * ff^dag>| \t err(|<ff * ff^dag>|) \t Re{ff} \t err(Re{ff}) \t Im{ff} \t err(Im{ff})" << endl;
		output_file.close();
	}


	for (int ii = 0; ii < temp.size(); ii++) {
		
		stats_thesis(
			directories[ii] + gauge_files[ii],
			name_output_file_poly,
			"gauge",
			append_mode_poly[ii],
			temp[ii],
			n_skip_rep[ii],
			n_skip_imp[ii],
			dim_block_modP[ii],
			dim_block_reP[ii],
			dim_block_imP[ii],
			1
		);
		cout << "poly n°" << ii << " DONE! T = " << temp[ii] << endl;
		cout << endl;

		if (fermion_files[ii] == "NONE") continue;

		stats_thesis(
			directories[ii] + fermion_files[ii],
			name_output_file_ff,
			"fermion",
			append_mode_ff[ii],
			temp[ii],
			n_skip_reff[ii],
			n_skip_imff[ii],
			dim_block_modff[ii],
			dim_block_reff[ii],
			dim_block_imff[ii],
			n_copy[ii]
		);
		cout << "ff n°" << ii << " DONE! T = " << temp[ii] << endl;
		cout << endl;
	}
	
	return 0;
}

//-----------------------------------------------------------------
//FUNCTION DEFINITION:

void read_file_LPC(
	const string& name_file_lpc,
	const int Nt,
	const int skipLines_file_lpc,
	const double mpi,
	vector<double>& aml,
	vector<double>& beta,
	vector<double>& afm,
	vector<double>& temp,
	vector<int>& dim_block_modP,
	vector<int>& dim_block_reP,
	vector<int>& dim_block_imP,
	vector<int>& dim_block_modff,
	vector<int>& dim_block_reff,
	vector<int>& dim_block_imff
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
		int dim_block_modP_value, dim_block_reP_value, dim_block_imP_value;
		int dim_block_modff_value, dim_block_reff_value, dim_block_imff_value;
		string dir, gauge, ferm;
		if (iss >> aml_value >> beta_value >> afm_value >> dim_block_modP_value >> dim_block_reP_value >> dim_block_imP_value >> dim_block_modff_value >> dim_block_reff_value >> dim_block_imff_value) {
			aml.push_back(aml_value);
			beta.push_back(beta_value);
			afm.push_back(afm_value);
			dim_block_modP.push_back(dim_block_modP_value);
			dim_block_reP.push_back(dim_block_reP_value);
			dim_block_imP.push_back(dim_block_imP_value);
			dim_block_modff.push_back(dim_block_modff_value);
			dim_block_reff.push_back(dim_block_reff_value);
			dim_block_imff.push_back(dim_block_imff_value);
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
	vector<string>& gauge_files,
	vector<string>& fermion_files,
	vector<int>& n_copy,
	vector<int>& n_skip_rep,
	vector<int>& n_skip_imp,
	vector<int>& n_skip_reff,
	vector<int>& n_skip_imff
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
		int n_therm_rep, n_therm_imp, n_therm_reff, n_therm_imff, n_copy_value;
		if (iss >> dir >> gauge >> n_therm_rep >> n_therm_imp >> ferm >> n_therm_reff >> n_therm_imff >> n_copy_value) {
			directories.push_back(dir);
			gauge_files.push_back(gauge);
			fermion_files.push_back(ferm);
			n_skip_rep.push_back(skipLines + n_therm_rep);
			n_skip_imp.push_back(skipLines + n_therm_imp);
			n_skip_reff.push_back(skipLines + n_therm_reff * n_copy_value);
			n_skip_imff.push_back(skipLines + n_therm_imff * n_copy_value);
			n_copy.push_back(n_copy_value);
		}
		else {
			cerr << "Poorly formatted line: " << line << endl;
		}
	}

	file_list.close();
}

void stats_thesis(
	const string& input_path,
	const string& output_path,
	const string& tipology, //fermion/gauge
	const bool append_mode,
	const double temp,
	int n_skip_re,
	int n_skip_im,
	int dim_block,
	int dim_block_re,
	int dim_block_im,
	int n_copy
){
	vector <double> y, yr, yi;
	double discard1, discard2, discard3, discard4;
	double obs, obs_re, obs_im;
	double mean, mean_re, mean_im, var_m, var_re, var_im;
	string line;

	ifstream input_file; //declaration of input file
	input_file.open(input_path);
	if (!input_file) {
		cout << "Error opening input file: " << input_path << endl;
		return;
	}

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

	for (int jj = 0; jj < min(n_skip_re, n_skip_im); jj++) {
		if (!getline(input_file, line)) {
			cerr << "Error: there are less than " << min(n_skip_re, n_skip_im) << " lines in the file: " << input_path << endl;
			return;
		}
	}

	if (n_skip_re < n_skip_im) {
		for (int jj = 0; jj < (n_skip_im - n_skip_re); jj++) {
			getline(input_file, line);
			if (line.empty()) {
				cerr << "Skipped blank/whitespace-only line in file: " << input_path << endl;
				continue;
			}
			istringstream iss(line);
			if (iss >> discard1 >> discard2 >> discard3 >> discard4 >> obs_re >> obs_im) {
				yr.push_back(obs_re);
			}
			else {
				cerr << "Skipped line (badly formatted) from" << input_path << ": " << line << endl;
			}
		}
	}
	else if (n_skip_im < n_skip_re) {
		for (int jj = 0; jj < (n_skip_re - n_skip_im); jj++) {
			getline(input_file, line);
			if (line.empty()) {
				cerr << "Skipped blank/whitespace-only line in file: " << input_path << endl;
				continue;
			}
			istringstream iss(line);
			if (iss >> discard1 >> discard2 >> discard3 >> discard4 >> obs_re >> obs_im) {
				yi.push_back(obs_im);
			}
			else {
				cerr << "Skipped line (badly formatted) from" << input_path << ": " << line << endl;
			}
		}
	}

	while (getline(input_file, line)) {
		if (line.empty()) {
			cerr << "Skipped blank/whitespace-only line in file: " << input_path << endl;
			continue;
		}
		istringstream iss(line);
		if (iss >> discard1 >> discard2 >> discard3 >> discard4 >> obs_re >> obs_im) {
			obs = obs_re * obs_re + obs_im * obs_im;//obs = |obs|^2 = obs_re^2 + obs_im^2;
			y.push_back(obs);
			yr.push_back(obs_re);
			yi.push_back(obs_im);
		}
		else {
			cerr << "Skipped line (badly formatted) from" << input_path << ": " << line << endl;
		}
	}

	dim_block *= n_copy;
	dim_block_re *= n_copy;
	dim_block_im *= n_copy;

	if (dim_block == 1) {
		stats_indipendent_unbiased(&mean, &var_m, y);
	}
	else {
		blocking(&mean, &var_m, y, dim_block);
	}
	if (dim_block_re == 1) {
		stats_indipendent_unbiased(&mean_re, &var_re, yr);
	}
	else {
		blocking(&mean_re, &var_re, yr, dim_block_re);
	}
	if (dim_block_im == 1) {
		stats_indipendent_unbiased(&mean_im, &var_im, yi);
	}
	else
	{
		blocking(&mean_im, &var_im, yi, dim_block_im);
	}

	cout << tipology << ": |<Obs * Obs^dag>| = " << mean << " +- " << sqrt(var_m) << endl;
	cout << tipology << ": <Re{Obs}> = " << mean_re << " +- " << sqrt(var_re) << endl;
	cout << tipology << ": <Im{Obs}> = " << mean_im << " +- " << sqrt(var_im) << endl;

	output_file << setprecision(numeric_limits<double>::max_digits10);
	output_file << temp << "\t" << mean << "\t" << sqrt(var_m) << "\t" << mean_re << "\t" << sqrt(var_re) << "\t" << mean_im << "\t" << sqrt(var_im) << endl;

	input_file.close();
	output_file.close();
}
