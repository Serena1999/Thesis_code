/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****            2nd_susceptibility_with_errors.cpp:           ****
****        this script takes the original measurement        ****
****  files for each temperature and returns a file with the  ****
****  temperature along with the mean and standard deviation  ****
****        of the suscebility of the wanted variable.        ****
****    When selected, blocking is applied with block sizes   ****
****             selectable from an external file.            ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

//here we have mod(obs), re(obs), im(obs)

//-----------------------------------------------------------------
//HEADERS && LIBRARY:

#include "../library.h"
#include "../root_include.h"

//-----------------------------------------------------------------
//GLOBAL CONSTANTS:

const double hbar_c = 197.3269804; //MeV * fm

//-----------------------------------------------------------------
//VARIABLES TO SET:

const bool debug_mode = 0;

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
	vector<int>& n_skip_rep,
	vector<int>& n_skip_imp,
	vector<int>& n_skip_reff,
	vector<int>& n_skip_imff,
	int step_sample_gauge,
	int step_sample_fermion
);

void susceptibility_with_errors(
	const string tipology,
	const int n_steps,
	const string& input_path,
	const string& output_path_mod,
	const string& output_path_re,
	const string& output_path_im,
	const double temp,
	const double prefactor,
	const int n_skip_re,
	const int n_skip_im,
	int dim_block_mod,
	int dim_block_re,
	int dim_block_im,
	int n_skip_file
);

template <class T> int blocking_sample(
	vector <T>& original_draws,
	vector <T>& blocked_draws,
	const int dim_block
);

void mean_chi_f(
	sample_gen& sampler,
	const int n_steps,
	vector <double>& original_draws,
	vector <double>& blocked_draws,
	double prefactor,
	double* chi_estimate,
	double* var_chi
);

bool read_measurements(
	const string input_path,
	vector <double>& original_draws_re,
	vector <double>& original_draws_im,
	vector <double>& original_draws_mod,
	const int n_skip_file,
	const int n_skip_re,
	const int n_skip_im
);

//-----------------------------------------------------------------
//MAIN:

int main() {
	//-----------------------------------------------------------------
	// VARIABILES TO SET:
	int n_steps = 1100;//TO CHOOSE: NUMBER OF BOOSTRAP STEPS 
	int Nt = 8; //BE CAREFUL TO CHOOSE IT WELL;
	const int Ns = 32; //BE CAREFUL TO CHOOSE IT WELL;
	const int Vs = Ns * Ns * Ns;
	int skipLines_file_lpc = 2, skipLines_file_list = 1, skipLines = 1;
	int step_sample_fermion = 10;
	int step_sample_gauge = 1;
	double mpi = 1500; //MeV //BE CAREFUL TO CHOOSE IT WELL;
	bool bool_startFile_poly = 1, bool_startFile_ff = 1;//BE CAREFUL TO CHOOSE IT WELL;
	
	string name_file_lpc = "19_05_2025/LCP_1500MeV_dimblock.txt";
	string name_file_list = "19_05_2025/file_list_therm.txt";

	//-----------------------------------------------------------------
	// DO NOT MODIFY THE FOLLOWING:
	double temp_value;
	vector<int> n_skip_rep, n_skip_imp, n_skip_reff, n_skip_imff;
	vector<int> dim_block_modP, dim_block_reP, dim_block_imP, dim_block_modff, dim_block_reff, dim_block_imff;
	vector<double> aml, beta, afm, temp;//T = \hbar * c /(Nt * a[fm]) (1.60), Nt = 8; 
	vector<string> directories, gauge_files, fermion_files;
	ostringstream mpi_stream;//TO INTRODUCE ALSO IN NUMERICAL METHODS CODE: IT IS USEFUL;
	mpi_stream << fixed << setprecision(1) << mpi; //set to 1 decimal place
	string mpi_string = mpi_stream.str(); // conversion into string
	string name_output_file_poly_re = "results/" + mpi_string + "_reP_results.txt";
	string name_output_file_poly_im = "results/" + mpi_string + "_imP_results.txt";
	string name_output_file_poly_mod = "results/" + mpi_string + "_modP_results.txt";
	string name_output_file_ff_re = "results/" + mpi_string + "_reff_results.txt";
	string name_output_file_ff_im = "results/" + mpi_string + "_imff_results.txt";
	string name_output_file_ff_mod = "results/" + mpi_string + "_modff_results.txt";

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
		n_skip_rep,
		n_skip_imp,
		n_skip_reff,
		n_skip_imff,
		step_sample_gauge,
		step_sample_fermion
	);

	if (bool_startFile_poly) {
		ofstream output_file_reP; //declaration of output file
		output_file_reP.open(name_output_file_poly_re);
		if (!output_file_reP) {
			cerr << "Error opening reP output file" << endl;
			return 1;
		}
		output_file_reP << "# T \t <Chi_Re{P}> \t err(Chi_Re{P})" << endl;
		output_file_reP.close();

		ofstream output_file_imP; //declaration of output file
		output_file_imP.open(name_output_file_poly_im);
		if (!output_file_imP) {
			cerr << "Error opening imP output file" << endl;
			return 1;
		}
		output_file_imP << "# T \t <Chi_Im{P}> \t err(Chi_Im{P})" << endl;
		output_file_imP.close();

		ofstream output_file_modP; //declaration of output file
		output_file_modP.open(name_output_file_poly_mod);
		if (!output_file_modP) {
			cerr << "Error opening mod{P} output file" << endl;
			return 1;
		}
		output_file_modP << "# T \t <Chi_|P|> \t err(Chi_|P|)" << endl;
		output_file_modP.close();
	}

	if (bool_startFile_ff) {
		ofstream output_file_reff; //declaration of output file
		output_file_reff.open(name_output_file_ff_re);
		if (!output_file_reff) {
			cerr << "Error opening Re{ff} output file" << endl;
			return 1;
		}
		output_file_reff << "# T \t <Chi_Re{ff}> \t err(Chi_Re{ff})" << endl;
		output_file_reff.close();

		ofstream output_file_imff; //declaration of output file
		output_file_imff.open(name_output_file_ff_im);
		if (!output_file_imff) {
			cerr << "Error opening Im{ff} output file" << endl;
			return 1;
		}
		output_file_imff << "# T \t <Chi_Im{ff}> \t err(Chi_Im{ff})" << endl;
		output_file_imff.close();

		ofstream output_file_modff; //declaration of output file
		output_file_modff.open(name_output_file_ff_mod);
		if (!output_file_modff) {
			cerr << "Error opening mod{ff} output file" << endl;
			return 1;
		}
		output_file_modff << "# T \t <Chi_|ff|> \t err(Chi_|ff|)" << endl;
		output_file_modff.close();
	}

	//cout << "temp size: " << temp.size() << endl;
	//cout << "directories size: " << directories.size() << endl;


	for (int ii = 0; ii < temp.size(); ii++) {

		susceptibility_with_errors(
			"gauge",
			n_steps,
			directories[ii] + gauge_files[ii],
			name_output_file_poly_mod,
			name_output_file_poly_re,
			name_output_file_poly_im,
			temp[ii],
			Vs,
			n_skip_rep[ii],
			n_skip_imp[ii],
			dim_block_modP[ii],
			dim_block_reP[ii],
			dim_block_imP[ii],
			skipLines
		);

		cout << "poly n°" << ii << " DONE! T = " << temp[ii] << endl;
		cout << endl;

		if (fermion_files[ii] == "NONE") continue;

		susceptibility_with_errors(
			"fermion",
			n_steps,
			directories[ii] + fermion_files[ii],
			name_output_file_ff_mod,
			name_output_file_ff_re,
			name_output_file_ff_im,
			temp[ii],
			Vs,
			n_skip_rep[ii],
			n_skip_imp[ii],
			dim_block_modff[ii],
			dim_block_reff[ii],
			dim_block_imff[ii],
			skipLines
		);

		cout << "ff n°" << ii << " DONE! T = " << temp[ii] << endl;
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

	if (debug_mode) {
		cout << "LCP FILE:" << endl;
		cout << "T \t aml \t beta \t a[fm] \t dim_block_modP \t dim_block_reP \t dim_block_imP \t dim_block_modff \t dim_block_reff \t dim_block_imff:" << endl;

		for (int ii = 0; ii < beta.size(); ++ii) {
			cout << temp[ii] << "\t" << aml[ii] << "\t" << beta[ii] << "\t" << afm[ii] << "\t" <<
				dim_block_modP[ii] << "\t" << dim_block_reP[ii] << "\t" << dim_block_imP[ii] << "\t" <<
				dim_block_modff[ii] << "\t" << dim_block_reff[ii] << "\t" << dim_block_imff[ii] << endl;
		}
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
	vector<int>& n_skip_rep,
	vector<int>& n_skip_imp,
	vector<int>& n_skip_reff,
	vector<int>& n_skip_imff,
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
		int n_therm_rep, n_therm_imp, n_therm_reff, n_therm_imff;
		if (iss >> dir >> gauge >> n_therm_rep >> n_therm_imp >> ferm >> n_therm_reff >> n_therm_imff) {
			directories.push_back(dir);
			gauge_files.push_back(gauge);
			fermion_files.push_back(ferm);
			n_skip_rep.push_back(n_therm_rep / step_sample_gauge);
			n_skip_imp.push_back(n_therm_imp / step_sample_gauge);
			n_skip_reff.push_back(n_therm_reff / step_sample_fermion);
			n_skip_imff.push_back(n_therm_imff / step_sample_fermion);
		}
		else {
			cerr << "Poorly formatted line: " << line << endl;
		}
	}

	if (debug_mode) {
		cout << "LIST FILE:" << endl;
		cout << "(step_sample_gauge, step_sample_fermion) = (" << step_sample_gauge << ", " << step_sample_fermion << ")" << endl;
		cout << "directory \t gauge_file \t n_skip_reP \t n_skip_imP \t ferm_file \t n_skip_reff \t n_skip_imff:" << endl;

		for (int ii = 0; ii < directories.size(); ++ii) {
			cout << directories[ii] << "\t" << gauge_files[ii] << "\t" << n_skip_rep[ii] << "\t" << n_skip_imp[ii] << "\t" <<
				fermion_files[ii] << "\t" << n_skip_reff[ii] << "\t" << n_skip_imff[ii] << endl;
		}
	}

	file_list.close();
}

void susceptibility_with_errors(
	const string tipology,
	const int n_steps,
	const string& input_path,
	const string& output_path_mod,
	const string& output_path_re,
	const string& output_path_im,
	const double temp,
	const double prefactor,
	const int n_skip_re,
	const int n_skip_im,
	int dim_block_mod,
	int dim_block_re,
	int dim_block_im,
	int n_skip_file
) {
	sample_gen sampler;
	vector <double> original_draws_re, original_draws_im, original_draws_mod;
	vector <double> blocked_draws_re, blocked_draws_im, blocked_draws_mod;

	read_measurements(
		input_path,
		original_draws_re,
		original_draws_im,
		original_draws_mod,
		n_skip_file,
		n_skip_re,
		n_skip_im
		);

	if (blocking_sample(original_draws_re, blocked_draws_re, dim_block_re)) {
		cout << "Problem in blocking of re_original_draws" << endl;
		return;
	}

	if (blocking_sample(original_draws_im, blocked_draws_im, dim_block_im)) {
		cout << "Problem in blocking of im_original_draws" << endl;
		return;
	}

	if (blocking_sample(original_draws_mod, blocked_draws_mod, dim_block_mod)) {
		cout << "Problem in blocking of mod_original_draws" << endl;
		return;
	}


	//mean computation:

	double estimator_chi_re, estimator_chi_im, estimator_chi_mod;
	double var_chi_re, var_chi_im, var_chi_mod;

	mean_chi_f(sampler, n_steps, original_draws_re, blocked_draws_re, prefactor, &estimator_chi_re, &var_chi_re);


	if (var_chi_re < 0) {
		cout << "[WARNING] var_chi_re problem at T = " << temp << ", var_chi_re = " << var_chi_re << endl;
		cout << "dim_block = " << dim_block_re << endl;
		cout << "len original = " << original_draws_re.size() << endl;
		cout << "len blocked = " << blocked_draws_re.size() << endl;
	}

	ofstream output_file_re; //declaration of output file
	output_file_re.open(output_path_re, ios::app);
	if (!output_file_re) {
		cout << tipology << ": Error opening re output file" << endl;
		return;
	}

	cout << tipology << ": <Chi_re> = " << estimator_chi_re << " +- " << sqrt(var_chi_re) << endl;

	output_file_re << fixed << setprecision(numeric_limits<double>::max_digits10);
	output_file_re << temp << "\t" << estimator_chi_re << "\t" << sqrt(var_chi_re) << endl;

	original_draws_re.clear();
	blocked_draws_re.clear();

	output_file_re.close();


	mean_chi_f(sampler, n_steps, original_draws_im, blocked_draws_im, prefactor, &estimator_chi_im, &var_chi_im);
	
	if (var_chi_im < 0) {
		cout << "[WARNING] var_chi_im problem at T = " << temp << ", var_chi_im = " << var_chi_im << endl;
		cout << "dim_block = " << dim_block_im << endl;
		cout << "len original = " << original_draws_im.size() << endl;
		cout << "len blocked = " << blocked_draws_im.size() << endl;
	}

	ofstream output_file_im; //declaration of output file
	output_file_im.open(output_path_im, ios::app);
	if (!output_file_im) {
		cout << tipology << ": Error opening im output file" << endl;
		return;
	}

	cout << tipology << ": <Chi_im> = " << estimator_chi_im << " +- " << sqrt(var_chi_im) << endl;

	output_file_im << fixed << setprecision(numeric_limits<double>::max_digits10);
	output_file_im << temp << "\t" << estimator_chi_im << "\t" << sqrt(var_chi_im) << endl;

	original_draws_im.clear();
	blocked_draws_im.clear();

	output_file_im.close();


	mean_chi_f(sampler, n_steps, original_draws_mod, blocked_draws_mod, prefactor, &estimator_chi_mod, &var_chi_mod);
	if (var_chi_mod < 0) {
		cout << "[WARNING] var_chi_mod problem at T = " << temp << ", var_chi_mod = " << var_chi_mod << endl;
		cout << "dim_block = " << dim_block_mod << endl;
		cout << "len original = " << original_draws_mod.size() << endl;
		cout << "len blocked = " << blocked_draws_mod.size() << endl;
	}
	
	ofstream output_file_mod; //declaration of output file
	output_file_mod.open(output_path_mod, ios::app);
	if (!output_file_mod) {
		cout << tipology << ": Error opening mod output file" << endl;
		return;
	}

	cout << tipology << ": <Chi_mod> = " << estimator_chi_mod << " +- " << sqrt(var_chi_mod) << endl;

	output_file_mod << fixed << setprecision(numeric_limits<double>::max_digits10);
	output_file_mod << temp << "\t" << estimator_chi_mod << "\t" << sqrt(var_chi_mod) << endl;

	original_draws_mod.clear();
	blocked_draws_mod.clear();

	output_file_mod.close();

}


void mean_chi_f(
	sample_gen& sampler,
	const int n_steps,
	vector <double>& original_draws,
	vector <double>& blocked_draws,
	double prefactor,
	double* chi_estimate,
	double* var_chi
) {
	double mean = 0, mean2 = 0, mean_chi2 = 0, mean_chi = 0, value;
	double mean_tmp, mean_tmp2;
	int N = original_draws.size(), index;

	for (int ii = 0; ii < N; ii++) {
		value = original_draws[ii];
		mean += value;
		mean2 += (value * value);
	}
	mean /= (double)original_draws.size();
	mean2 /= (double)original_draws.size();
	(*chi_estimate) = obs_function(prefactor, mean, mean2);

	//variance computation:

	N = blocked_draws.size();
	uniform_int_distribution<> dist_int(0, N - 1);
	//sampler.init(10);//SEED

	for (int ii = 0; ii < n_steps; ii++) {

		mean_tmp = 0;
		mean_tmp2 = 0;

		for (int jj = 0; jj < N; jj++) {
			index = dist_int(sampler.rng);
			value = blocked_draws[index];
			mean_tmp += value;
			mean_tmp2 += (value * value);
			//cout << value << endl;
		}

		mean_tmp /= (double)blocked_draws.size();
		mean_tmp2 /= (double)blocked_draws.size();

		value = obs_function(prefactor, mean_tmp, mean_tmp2);
		mean_chi += value;
		mean_chi2 += value * value;
	}

	mean_chi /= (double)n_steps;
	mean_chi2 /= (double)n_steps;
	(*var_chi) = mean_chi2 - mean_chi * mean_chi;
	//cout << endl;
	//cout << prefactor << endl;
	//cout << mean_tmp << endl;
	//cout << mean_tmp2 << endl;
	//cout << *chi_estimate << endl;
	//cout << n_steps << endl;
	//cout << N << endl;
	//cout << (*var_chi) << endl;
	assert((*var_chi) > 0);
	(*var_chi) *= (double)n_steps;
	(*var_chi) /= (double)(n_steps - 1);
	assert((*var_chi) > 0);

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

	//cout << endl;
	//cout << n_blocks << endl;
	//cout << endl;

	int n_max = n_blocks * dim_block;//N_max to consider to compute variance
	int index = 0;
	T mean_tmp, value;

	for (int jj = 0; jj < n_max; jj += dim_block) {
		index++;
		mean_tmp = 0;
		for (int kk = 1; kk <= dim_block; kk++) {
			value = original_draws[jj + kk - 1];
			mean_tmp += value;
		}
		mean_tmp /= (double)dim_block;
		blocked_draws.push_back(mean_tmp);

		//cout << mean_tmp << endl;

	}

	//cout << endl;

	return 0;
}

bool read_measurements(
	const string input_path,
	vector <double>& original_draws_re,
	vector <double>& original_draws_im,
	vector <double>& original_draws_mod,
	const int n_skip_file,
	const int n_skip_re,
	const int n_skip_im
) {

	original_draws_re.clear();
	original_draws_im.clear();
	original_draws_mod.clear();

	double obs, obs_re, obs_im;
	string discard1, discard2, discard3, discard4;
	string line;

	ifstream input_file; //declaration of input file
	input_file.open(input_path);
	if (!input_file) {
		cout << "Error opening input file: " << input_path << endl;
		return 1;
	}

	double  value_re_tmp = 0, value_im_tmp = 0, value_mod_tmp = 0;
	int conf_id = -999, conf_tmp = 999, n_copy_accum_r = 0, n_copy_accum_i = 0, n_copy_accum_m = 0, flag = 0, counter = 0;
	int cnt = 0, flag_init = 0;

	for (int jj = 0; jj < n_skip_file; jj++) {
		if (!getline(input_file, line)) {
			cerr << "Error: there are less than " << n_skip_file << " lines in the file: " << input_path << endl;
			return 1;
		}
	}

	for (int jj = 0; jj < min(n_skip_re, n_skip_im); jj++) {
		while (conf_id != conf_tmp) {
			if (!getline(input_file, line)) {
				cerr << "Error: there are less than " << min(n_skip_re, n_skip_im) + n_skip_file << " lines in the file: " << input_path << endl;
				return 1;
			}
			istringstream iss(line);
			if (iss >> conf_id) {
				if (jj == 0) {
					conf_tmp = conf_id;
				}
			}
		}
	}

	if (n_skip_re < n_skip_im) {
		while (counter < (n_skip_im - n_skip_re)) {
			if (!getline(input_file, line)) break;
			istringstream iss(line);
			if (iss >> conf_id >> discard1 >> discard2 >> discard3 >> obs_re >> obs_im) {
				if (!flag_init) {
					conf_tmp = conf_id;
					flag_init = 1;
				}
				if (conf_id != conf_tmp) {
					value_re_tmp /= (double)n_copy_accum_r;
					original_draws_re.push_back(value_re_tmp);

					counter++;
					conf_tmp = conf_id;
					value_re_tmp = 0;
					n_copy_accum_r = 0;
				}
				n_copy_accum_r++;
				value_re_tmp += obs_re;
			}
			else {
				cerr << "Skipped line (badly formatted) from" << input_path << ": " << line << endl;
			}
		}
	}
	else if (n_skip_im < n_skip_re) {
		while (counter < (n_skip_re - n_skip_im)) {
			if (!getline(input_file, line)) break;
			istringstream iss(line);
			if (iss >> conf_id >> discard1 >> discard2 >> discard3 >> obs_re >> obs_im) {
				if (!flag_init) {
					conf_tmp = conf_id;
					flag_init = 1;
				}
				if (conf_id != conf_tmp) {
					value_im_tmp /= (double)n_copy_accum_i;
					original_draws_im.push_back(value_im_tmp);
					
					counter++;
					conf_tmp = conf_id;
					value_im_tmp = 0;
					n_copy_accum_i = 0;
				}
				n_copy_accum_i++;
				value_im_tmp += obs_im;
			}
			else {
				cerr << "Skipped line (badly formatted) from" << input_path << ": " << line << endl;
			}
		}
	}

	while (getline(input_file, line)) {
		istringstream iss(line);
		if (iss >> conf_id >> discard1 >> discard2 >> discard3 >> obs_re >> obs_im) {
			if (!flag_init) {
				conf_tmp = conf_id;
				flag_init = 1;
			}
			flag = 0;
			if ((debug_mode) && (input_path == "19_05_2025/build_good_mpi1500_12_out/ferm_obs230875197.txt")) {
				cout << "(conf_id, conf_tmp) = (" << conf_id << ", " << conf_tmp << ")" << endl;
			}

			if (conf_id != conf_tmp) {
				if (n_copy_accum_i != 0) {

					value_im_tmp /= (double)n_copy_accum_i;
					original_draws_im.push_back(value_im_tmp);
				}
				if (n_copy_accum_r != 0) {

					value_re_tmp /= (double)n_copy_accum_r;					
					original_draws_re.push_back(value_re_tmp);
				}
				if (n_copy_accum_m != 0) {

					value_mod_tmp /= (double)n_copy_accum_m;
					original_draws_mod.push_back(value_mod_tmp);
				}

				cnt++;
				conf_tmp = conf_id;

				if ((debug_mode) && (input_path == "19_05_2025/build_good_mpi1500_12_out/ferm_obs230875197.txt")) {
					cout << "Some acculumulated data in fermion file: " << cnt << "\t" << value_mod_tmp << "\t" << value_re_tmp << "\t" << value_im_tmp << endl;
				}
				value_mod_tmp = 0;
				value_re_tmp = 0;
				value_im_tmp = 0;
				n_copy_accum_m = 0;
				n_copy_accum_r = 0;
				n_copy_accum_i = 0;
			}
			n_copy_accum_m++;
			n_copy_accum_r++;
			n_copy_accum_i++;
			value_mod_tmp += sqrt(obs_re * obs_re + obs_im * obs_im);//obs = |obs| = sqrt(obs_re^2 + obs_im^2);
			value_re_tmp += obs_re;
			value_im_tmp += obs_im;
		}
		else {
			flag = 1;
			cerr << "Skipped line (badly formatted) from" << input_path << ": " << line << endl;
		}
	}
	if ((conf_id == conf_tmp) && (!flag)) {
		value_re_tmp /= n_copy_accum_r;
		value_im_tmp /= n_copy_accum_i;
		value_mod_tmp /= n_copy_accum_m;

		original_draws_mod.push_back(value_mod_tmp);
		original_draws_re.push_back(value_re_tmp);
		original_draws_im.push_back(value_im_tmp);
	}

	return 0;
}
