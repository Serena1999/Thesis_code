/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
**** dependent_data_analysis_autocorrelation_calculated.cpp:  ****
****   implementation of the calculation of autocorrelation   ****
****                                                          ****
****  Autocorrelation is considered by explicit calculation.  ****
****    Autocorrelation vs (difference between indices) are   ****
****                  written on a text file.                 ****
****                                                          ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

#include "../library.h"
#include <cassert>

template <class T> bool renormalize_conf_draw(
	int n_first_conf,
	vector <int>& conf_draws,
	vector <T>& original_draws,
	vector <T>& renormalized_conf_draws
);

//-----------------------------------------------------------------
//MAIN:

int main() {

	int skipLines = 0; //= number of lines to skip while reading input file;
	int n_term = 2000;// 50000; // 200;//=number of data to discard because of thermalization;
	int step_sample = 10; //10 for fermion obs, 1 for gauge ones
	n_term = floor(n_term/step_sample);
	//cout << n_term << endl;
	int N; //number of draws on which we do the statistics;
	double mean = 0, mean2 = 0, var_m = 0, dev_m;
	double tau_int = 0;
	int kmax;//k up to which the summation is done (selectable). 
					//The maximum possible is kmax = N, where N = (number of considered draws). 
	double value;
	vector<double> original_draws, renormalized_conf_draws, conf_draws_d;
	vector <int> conf_draws;

	string name_input_directory = "19_05_2025/data_value/";
	string name_input_file = "reff_ferm_obs2239384639.txt";
	string input_path = name_input_directory + name_input_file;
	string name_output_fileCk = "results/dependent_data_analysis_Ck_" + name_input_file;
	
	ifstream input_file; //declaration of input file
	int index = 0;
	string line;

	read_2columns(skipLines, conf_draws_d, original_draws, input_path);
	
	for (int ii = 0; ii < conf_draws_d.size(); ii++) {
		conf_draws.push_back((int) conf_draws_d[ii]);
	}


	if (renormalize_conf_draw(n_term, conf_draws, original_draws, renormalized_conf_draws)) {
		cerr << "Error in renormalize_conf_draw" << endl;
		conf_draws.clear();
		original_draws.clear();
		renormalized_conf_draws.clear();
		return 1;
	}

	N = renormalized_conf_draws.size();
	//cout << N << endl;

	for (int ii = 0; ii < N; ii++) {
		value = renormalized_conf_draws[ii];
		mean += value;
		mean2 += (value * value);
	}

	mean /= (double) N;
	mean2 /= (double)N;
	
	var_m = mean2 - mean * mean;
	var_m /= (double)(N - 1);

	kmax = N;
	
	vector<double> corr(kmax, 0); // = <dFi * dFj> //create a vector of type double with kmax elements, all initialized to the value 0.
	vector<double> ck(kmax, 0); //array for C(k)
	vector<int> ik(kmax, 0); //numbers of data on which the sample mean of C(k) is taken for given k

	ofstream output_fileCk;

	output_fileCk.open(name_output_fileCk);
	if (!output_fileCk) {
		cout << "Error opening Ck vs k file" << endl;
		return 1;
	}

	for (int i = 0; i < N; i++) {
		corr[0] += (renormalized_conf_draws[i] - mean) * (renormalized_conf_draws[i] - mean);
	}
	//cout << N << endl;//for debug;
	ik[0] = N;
	corr[0] = corr[0] / ((double)ik[0]);
	ck[0] = 1.0;

	output_fileCk << "#k\tC[k]" << endl;
	output_fileCk << fixed << setprecision(numeric_limits<double>::max_digits10);//This way, the output to the file will show the double values ​​with the greatest possible precision.
	output_fileCk << 0 << "\t" << ck[0] << endl;

	for (int k = 1; k < kmax; k++) {
		for (int i = 0; (i + k) < N; i++) {
			corr[k] += (renormalized_conf_draws[i] - mean) * (renormalized_conf_draws[i + k] - mean);
			//cout << renormalized_conf_draws[i] << endl;
			//cout << renormalized_conf_draws[i + k] << endl;
			//cout << mean << endl;
			//cout << endl;
		}
		//cout << k << endl;//for debug;
		ik[k] = N - k;
		corr[k] = corr[k] / ((double)ik[k]);
		ck[k] = corr[k] / corr[0];//normalized in a way in which ck[0] = 1;
		if (ck[k] < 0) break;
		output_fileCk << k << "\t" << ck[k] << endl;
		//cout << k << endl;//for debug;
	};//breaking out of the loop: 
	// -> corr[k] = <d(draws_i) d(draws_{i+k})>, with 0<=k<=kmax;
	// -> ck[k] = <d(draws_i) d(draws_{i+k})>/sigma_F ^2, with 0<=k<=kmax;
	// -> tau_int = \sum_{k=1} ^kmax C(k)

	output_fileCk.close();

	return 0;
}


template <class T> bool renormalize_conf_draw(
	int n_first_conf,
	vector <int>& conf_draws,
	vector <T>& original_draws,
	vector <T>& renormalized_conf_draws
) {
	int index = 0;
	bool flag = 1;
	int cnt = 0, n_copy_accum = 0, conf_id = conf_draws[0], conf_tmp;
	T value_tmp = 0;

	while (flag) {
		if ((conf_draws[index]) > n_first_conf) {
			flag = 0;
		}
		index++;
	}

	conf_tmp = conf_draws[index];

	for (int ii = index; ii < original_draws.size(); ii++) {
		conf_id = conf_draws[ii];
		if (conf_id != conf_tmp) {
			cnt++;
			conf_tmp = conf_id;
			value_tmp /= (double)n_copy_accum;
			renormalized_conf_draws.push_back(value_tmp);
			cout << "conf_id: " << conf_id << ", valore medio: " << value_tmp << ", n_accumulated = " << n_copy_accum << endl;
			value_tmp = 0;
			n_copy_accum = 0;
		}
		n_copy_accum++;
		value_tmp += original_draws[ii];
	}
	if (n_copy_accum != 0) {
		value_tmp /= (double)n_copy_accum;
		renormalized_conf_draws.push_back(value_tmp);
		cout << "conf_id: " << conf_id << ", valore medio: " << value_tmp << endl;
	}

	return flag;
}