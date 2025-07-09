/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                  analysis_rhok_vs_T.cpp:                 ****
****      Starting from dedicated files ("mon.dat"), the      ****
****    computes total monopole density per configuration.    ****
****    Averages over configurations at each temperature T.   ****
****    Plots <rho_tot/T^3> with its error vs temperature.    ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

//-----------------------------------------------------------------
//HEADERS && LIBRARY:

#include "../library.h"//DA STUDIARE BENE!!
#include "../root_include.h"

//-----------------------------------------------------------------
//GLOBAL CONSTANTS:

const double hbar_c = 197.3269804; //MeV * fm

//-----------------------------------------------------------------
//ROOT MACRO TO DO FIT AND GRAPH:

void plot_points_errors(
	vector<double>& x,
	vector<double>& y,
	vector<double>& y_err,
	string name_image,
	string title,
	string y_name,
	double pos_title,
	double pos_y,
	double heigh_y,
	double width_canvas
);

//-----------------------------------------------------------------
//FUNCTION DECLARATIONS AND DEFINITIONS:

void minimal_read_file_LPC(
	const string& name_file_lpc,
	const int Nt,
	const int skipLines_file_lpc,
	const double mpi,
	vector<double>& aml,
	vector<double>& beta,
	vector<double>& afm,
	vector<double>& temp
);


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

	string mpi = "1500";
	int mpi_int = 1500;
	string line;

	int Nt = 8;
	int Ns = 32;
	int ad_Vs = Ns * Ns * Ns; //adimensional spatial volume

	string list_file = "19_05_2025/file_list_therm.txt";
	int skipLines_list_file = 1;
	int step_sample_gauge = 1;
	int step_sample_fermion = 10;

	vector <string> directories;//LEGGILI DA FILE
	vector <int> n_skip; //LEGGILO DA FILE;
	string name_input_file = "mon.dat";
	int skip_lines_input_file = 0;

	string name_file_lpc = "19_05_2025/LCP_1500MeV_dimblock.txt";
	int skipLines_file_lpc = 2;
	vector <double> aml, beta, afm, temp;


	minimal_read_file_LPC(
		name_file_lpc,
		Nt,
		skipLines_file_lpc,
		mpi_int,
		aml,
		beta,
		afm,
		temp
	);

	vector <double> mean_rho_T3;
	/*
	-> mean_rho_T3[ii] = mean of rho_tot/T^3 over configurations = <rho_tot>/T^3
	*/
	vector <double> mean_2; //mean of squares
	vector <double> err_rho_T3; //standard deviation from the mean

	string name_image = "results/rho_tot_mpi" + mpi + "_attempt.png";
	string name_output_file = "results/rho_tot_mpi" + mpi + "_attempt.txt";

	string title = "#LT #rho / T^{3} #GT vs temperature:";

	string y_name = "#LT #rho / T^{3} #GT";

	double pos_title = 0.3;

	double pos_y = 0.04;

	double heigh_y = 0.4;

	double width_canvas = 800;

	read_file_list_DIRECTORIES_THERM(
		list_file,
		skipLines_list_file,
		directories,
		n_skip,
		step_sample_gauge,
		step_sample_fermion
	);

	for (int ii = 0; ii < directories.size(); ++ii) {

		mean_rho_T3.push_back(0);
		mean_2.push_back(0);

		ifstream input_file;
		input_file.open(directories[ii] + name_input_file);
		if (!input_file) {
			cerr << "Error opening input " << ii << "-th file." << endl;
			return 1;
		}

		//we don't consider the first skip_lines_input_file:
		for (int jj = 0; jj < skip_lines_input_file; ++jj) {
			if (!getline(input_file, line)) {
				cerr << "Error: there are less than " << skip_lines_input_file << " lines in the input " << ii << "-th file." << endl;
				return 1;
			}
		}

		//we don't condider the first n_skip[ii] configurations:

		double conf_id = -999, conf_tmp = 999;
		bool flag_init = 0;

		n_skip[ii] = 0;//MOMENTANEO

		for (int jj = 0; jj < n_skip[ii]; jj++) {
			while (conf_id != conf_tmp) {
				if (!getline(input_file, line)) {
					cerr << "Error: there are less than " << n_skip[ii] << " configurations in the " << ii << "-th input file. " << endl;
					return 1;
				}
				istringstream iss(line);
				if (iss >> conf_id) {
					if (jj == 0) {
						conf_tmp = conf_id;
						flag_init = 1;
					}
				}
			}
		}

		//for each configurazion, we compute the densities and then we do the mean:
		double discard1, discard2, discard3, discard4, discard5, discard6, discard7, wrap_value;
		int n_conf = 0;
		vector <double> nk;
		double rho_T3_value;

		while (getline(input_file, line)) {
			istringstream iss(line);
			if (iss >> conf_id >> discard1 >> discard2 >> discard3 >> discard4 >> discard5 >> discard6 >> discard7 >> wrap_value) {
				if (!flag_init) {
					conf_tmp = conf_id;
					flag_init = 1;
				}
				if (conf_id != conf_tmp) {
					//cout << conf_id << endl;
					//cout << conf_id << endl;
					//cout << discard1 << endl;
					//cout << discard2 << endl;
					//cout << discard3 << endl;
					//cout << discard4 << endl;
					//cout << discard5 << endl;
					//cout << discard6 << endl;
					//cout << discard7 << endl;
					//cout << wrap_value << endl;
					//cout << endl;
					if (nk.empty()) {
						cout << "No monopoles for conf_id = " << to_string(conf_tmp) << endl;
						cout << "--> so we omit this configuration." << endl;
						conf_tmp = conf_id;
					}
					else {
						++n_conf;
						rho_T3_value = 0;
						for (int kk = 0; kk < nk.size(); ++kk) {
							//rho_T3_value += kk * nk[kk];
							rho_T3_value += nk[kk];
						}
						mean_rho_T3[ii] += rho_T3_value;
						mean_2[ii] += (rho_T3_value * rho_T3_value);
					}
					nk.clear();
					conf_tmp = conf_id;
				}
				if (wrap_value == 0) {
					continue;
				}
				if (abs(wrap_value) > nk.size()) {
					nk.resize(abs(wrap_value), 0.0);
				}
				nk[abs(wrap_value) - 1] += 1;
			}
			else {
				cerr << "Skipped line (badly formatted) from" << (ii + 1) << "-th input file: " << line << endl;
			}
		}
		if (conf_id == conf_tmp) {
			if (nk.empty()) {
				cout << "No monopoles for conf_id = " << to_string(conf_tmp) << endl;
				cout << "--> so we omit this configuration." << endl;
			}
			else {
				++n_conf;
				rho_T3_value = 0;
				for (int kk = 0; kk < nk.size(); ++kk) {
					//rho_T3_value += kk * nk[kk];
					rho_T3_value += nk[kk];
				}
				mean_rho_T3[ii] += rho_T3_value;
				mean_2[ii] += (rho_T3_value * rho_T3_value);
			}
			nk.clear();
			//if (ii == 12) {
			//	cout << endl;
			//	//cout << kk << endl;
			//	cout << "mean_rhok.size() = " << mean_rhok.size() << endl;
			//	cout << endl;
			//}
		}

		if (n_conf == 0) {
			cout << "No configuration considerd for " << to_string(ii + 1) << "-th temperature: T = " << temp[ii] << endl;
			continue;
		}

		double factor = pow(Nt, 3) / ad_Vs;
		/*
		look at formula 32 of page 10 of second article given by professor.
		*/

		mean_rho_T3[ii] /= (double)n_conf;
		mean_2[ii] /= (double)n_conf;

		mean_rho_T3[ii] *= factor;
		mean_2[ii] *= (factor * factor);

		err_rho_T3.push_back(mean_2[ii] - mean_rho_T3[ii] * mean_rho_T3[ii]);

		if (n_conf == 1) {
			cout << "No sufficient configurations to evaluete error for configuration of " << to_string(ii + 1) << "-th temperature: T = " << temp[ii] << endl;
			continue;
		}

		err_rho_T3[ii] /= (double)(n_conf - 1);

		if (err_rho_T3[ii] < 0) {
			cout << "Negative variance for configuration of " << to_string(ii + 1) << "-th temperature: T = " << temp[ii] << endl;
			continue;
		}

		if (isnan(err_rho_T3[ii])) {
			cout << "Nan variance for configuration of " << to_string(ii + 1) << "-th temperature: T = " << temp[ii] << endl;
			continue;
		}

		err_rho_T3[ii] = sqrt(err_rho_T3[ii]);

		input_file.close();

		cout << mean_rho_T3[ii] << endl;
		cout << mean_2[ii] << endl;
		cout << n_conf << endl;
	}

	if (err_rho_T3.size() > 1) {

		plot_points_errors(
			temp,
			mean_rho_T3,
			err_rho_T3,
			name_image,
			title,
			y_name,
			pos_title,
			pos_y,
			heigh_y,
			width_canvas
		);
	}
	else {
		cout << "No image generated since:" << endl;
		cout << "Not enough data points to produce a plot (need at least 2)." << endl;
	}

	//output_file:

	ofstream output_file;
	output_file.open(name_output_file);
	if (!output_file) {
		cerr << "Error opening output file" << endl;
		return 1;
	}

	output_file << "#T[MeV] \t <rho_tot/T^3>[adimensional] \t  \t error(rho_tot/T^3):" << endl;

	output_file << setprecision(numeric_limits<double>::max_digits10);

	for (int ii = 0; ii < mean_rho_T3.size(); ++ii) {
		output_file << temp[ii] << "\t" << mean_rho_T3[ii] << "\t" << err_rho_T3[ii] << endl;
	}

	output_file.close();

	return 0;
}

//-----------------------------------------------------------------
//FUNCTION DEFINITIONS:

void minimal_read_file_LPC(
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
		cerr << "Error opening file lcp list" << endl;
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


//-----------------------------------------------------------------
//ROOT MACRO TO GRAPH:

void plot_points_errors(
	vector<double>& x,
	vector<double>& y,
	vector<double>& y_err,
	string name_image,
	string title,
	string y_name,
	double pos_title,
	double pos_y,
	double heigh_y,
	double width_canvas
) {

	if (x.size() != y.size() || y.size() != y_err.size()) {
		cerr << "Error: mismatched vector sizes in plot_points_errors()." << endl;
		return;
	}

	//CANVAS creation: (to draw the graph)
	TCanvas* canvas = new TCanvas("canvas", "Canvas for Drawing Points", width_canvas, 600);
	canvas->SetGrid();//to set grid

	// 1. Graph with only error bars:
	TGraphErrors* g_errors = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, y_err.data());

	//Using std::max_element: (useful to create my histogram object)
	//	-> std::max_element takes two iterators that define the range of the array to operate on. 
	//	-> returns an iterator that points to the maximum element found.
	//The usage of std::min_element is analogous but with for minimum element.
	auto min_x = *min_element(x.begin(), x.end());
	auto max_x = *max_element(x.begin(), x.end());
	auto min_y = *min_element(y.begin(), y.end()) - *max_element(y_err.begin(), y_err.end());
	auto max_y = *max_element(y.begin(), y.end()) + *max_element(y_err.begin(), y_err.end());

	g_errors->SetLineColor(kBlack);//color of error bars
	g_errors->SetMarkerColor(kBlack);
	g_errors->SetMarkerStyle(20);
	g_errors->SetTitle("");
	g_errors->GetXaxis()->SetLimits(min_x - 10, max_x + 10);
	//g_errors->GetYaxis()->SetRangeUser(min_y - 0.01 * fabs(min_y), max_y + 0.01 * fabs(max_y));
	g_errors->GetYaxis()->SetRangeUser(min_y - 0.01, max_y + 0.01);
	g_errors->Draw("AP");

	// 2. Graph with only the line joining the points:
	TGraph* g_line = new TGraph(x.size(), x.data(), y.data());
	g_line->SetLineColor(kBlack);
	g_line->SetLineStyle(2);
	g_line->Draw("LP SAME");

	gPad->Update();

	//LATEX: I add LaTeX titles and axis labels:
	TLatex latex;
	latex.SetNDC(); //sets the use of Normalized Device Coordinates (NDC).
	latex.SetTextSize(0.05); //changes text size for title
	latex.DrawLatex(pos_title, 0.94, title.c_str());
	latex.SetTextSize(0.04); //changes text size for axis labels
	latex.DrawLatex(0.45, 0.03, "T[MeV]");
	latex.SetTextAngle(90);
	latex.DrawLatex(pos_y, heigh_y, y_name.c_str());


	//SAVE: I save the canvas as an image
	canvas->SaveAs(name_image.c_str());

	//to save also in vectorial pdf form:
	canvas->SaveAs((name_image.substr(0, name_image.find_last_of(".")) + ".pdf").c_str());

	//DELETE:
	delete g_line;
	delete g_errors;
	delete canvas;
}