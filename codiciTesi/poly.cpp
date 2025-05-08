/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                         poly.cpp:                        ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

#include "../library.h"
#include "../root_include.h"

const double hbar_c = 197.3269804; //MeV * fm

void poly(double temp, string name_output_file, bool append_mode);

//0h/3h done: 
//usa k = 50 per fare il binning: DONE
// fai analisi anche del condensato chirale

//-----------------------------------------------------------------
//MAIN:

int main() {
	int Nt = 8; //BE CAREFUL TO CHOOSE IT WELL;
	double mpi = 800; //MeV //BE CAREFUL TO CHOOSE IT WELL;
	double temp_value;
	bool bool_startFile = 1;//BE CAREFUL TO CHOOSE IT WELL;
	//vector<double> aml = { 0.081869, 0.075858, 0.073104, 0.065886, 0.061722, 0.057988, 0.054603, 0.051498, 0.048612, 0.045891, 0.043293, 0.040783, 0.038333, 0.035928, 0.033558, 0.031653 }; //a*m_l //BE CAREFUL TO CHOOSE IT WELL;

	//vector<double> beta = { 3.80000, 3.82700, 3.84063, 3.88100, 3.90800, 3.93500, 3.96200, 3.98900, 4.01600, 4.04300, 4.07000, 4.09700, 4.12400, 4.15100, 4.17800, 4.20000 }; //BE CAREFUL TO CHOOSE IT WELL;

	//vector<double> afm = { 0.1377037468, 0.1290036669, 0.1249189781, 0.1139049357, 0.1073456174, 0.1013484202, 0.0958468116, 0.0907797794, 0.0860918322, 0.0817329985, 0.0776588275, 0.0738303888, 0.0702142723, 0.0667825883, 0.0635129674, 0.0609567905 }; //a[fm] //BE CAREFUL TO CHOOSE IT WELL;

	
	vector<double> aml = { 0.073104 };//a*m_l //BE CAREFUL TO CHOOSE IT WELL;
	vector<double> beta = { 3.84063 }; //BE CAREFUL TO CHOOSE IT WELL;
	vector<double> afm = { 0.1249189781 }; //a[fm] //BE CAREFUL TO CHOOSE IT WELL;

	vector<double> temp;//T = \hbar * c /(Nt * a[fm]) (1.60), Nt = 8; 
	vector<int> append_mode = { 1 };//BE CAREFUL TO CHOOSE IT WELL;

	ostringstream mpi_stream;//TO INTRODUCE ALSO IN NUMERICAL METHODS CODE: IT IS USEFUL;
	mpi_stream << std::fixed << std::setprecision(1) << mpi; //set to 1 decimal place
	string mpi_string = mpi_stream.str(); // conversion into string

	//string name_output_file = mpi_string + "_poly_results.txt";
	string name_output_file = mpi_string + "_conversionToSeeTemp.txt";

	for (int ii = 0; ii < (beta.size()); ii++) {
		temp_value = hbar_c / (Nt * afm[ii]);
		temp.push_back(temp_value);
		cout << temp_value << endl;
	}

	if (bool_startFile) {
		ofstream output_file; //declaration of output file
		output_file.open("results/" + name_output_file);
		if (!output_file) {
			cout << "Error opening output file" << endl;
		}

		output_file << "# T \t |<L* L^dag>| \t err(|<L* L^dag>|) \t Re{L} \t err(Re{L}) \t Im{L} \t err(Im{L})" << endl;

		output_file.close();
	}
	


	//poly(temp[0], name_output_file, append_mode[0]);
	return 0;
}

void poly(double temp, string name_output_file, bool append_mode) {

	int skipLines = 1; //= number of lines to skip while reading input file;
	int dim_block = 50; //TO CHOOSE BEFORE COMPILATION;
	vector<double> y; //to contain data distributed as gaussians
	vector <double> poly_vec, polyr_vec, polyi_vec;
	double poly, poly_re, poly_im, delta, delta_re, delta_im, value_tmp;
	double mean = 0, var_m = 0, mean_re = 0, var_re = 0, mean_im = 0, var_im = 0;
	int index = 0;
	string line;
	ifstream input_file; //declaration of input file
	ofstream output_file; //declaration of output file
	string name_input_file = "gauge_obs2277865449.txt";

	input_file.open("02_05_2025/" + name_input_file);
	if (!input_file) {
		cout << "Error opening input file" << endl;
	}

	if (append_mode) {
		output_file.open("results/" + name_output_file, ios::app);
	}
	else
	{
		output_file.open("results/" + name_output_file);
	}

	if (!output_file) {
		cout << "Error opening output file" << endl;
	}

	for (int i = 0; i<skipLines; i++) {
		if (!getline(input_file, line)) {
			cerr << "Error: there are less than " << skipLines << " lines in the file." << endl;
		}
	} 

	while(input_file >> value_tmp >> value_tmp >> value_tmp >> value_tmp >> poly_re >> poly_im) {//DEVI CORREGGEERE ERRORI
		index++;
		poly = poly_re * poly_re + poly_im * poly_im;//poly = |poly|^2 = poly_re^2 + poly_im^2;
		poly_vec.push_back(poly);
		polyr_vec.push_back(poly_re);
		polyi_vec.push_back(poly_im);
	}

	blocking_faster(&mean, &var_m, poly_vec, dim_block);
	blocking_faster(&mean_re, &var_re, polyr_vec, dim_block);
	blocking_faster(&mean_im, &var_im, polyi_vec, dim_block);

	cout << "|<L * L^dag>| = " << mean << " +- " << sqrt(var_m) << endl;
	cout << "<Re{L}> = " << mean_re << " +- " << sqrt(var_re) << endl;
	cout << "<Im{L}> = " << mean_im << " +- " << sqrt(var_im) << endl;

	output_file << temp << "\t" << mean << "\t" << sqrt(var_m) << "\t" << mean_re << "\t" << sqrt(var_re) << "\t" << mean_im << "\t" << sqrt(var_im) << endl;

	input_file.close();
	output_file.close();
}
