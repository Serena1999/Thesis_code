/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                    cluster_lenght.cpp:                   ****
****    Starting from dedicated files ("mon_clusters.dat"),   ****
****  for each configuration, it calculates the length of the ****
**** largest cluster (max_l) and the sum of the lengths of all****
****    clusters present (sum_l). With this information, it   ****
**** calculates the parameter r_c = max_l/sum_l for each gauge****
****  configuration. It then estimates the mean <r_c> and the ****
**** standard deviation from the mean over all configurations.****
****  It returns an output .txt file with temperature, <r_c>  ****
****          and standard deviation and an image for         ****
****                   <r_c> vs temperature.                  ****
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


void read_file_list_DIRECTORIES_THERM(//SCRIVI SOTTO COME IMPLEMENTARLA, PRENDI COME MINIMO DI TERMALIZZAZIONE IL MASSIMO FRA GLI N_SKIP, COS� SEI SICURA CHE VA BENE;
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
	int mpi_int = 800;
	string line;

	int Nt = 8;

	string list_file = "11_05_2025/file_list_therm_extended.txt";
	int skipLines_list_file = 1;
	int step_sample_gauge = 1;
	int step_sample_fermion = 10;

	vector <string> directories;
	vector <int> n_skip;
	string name_input_file = "mon_clusters.dat";
	int skip_lines_input_file = 0;

	string name_file_lpc = "11_05_2025/LCP_800MeV_dimblock_extended.txt";
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

	vector <double> mean_rc, err_rc;
	/*
	-> mean_rc[ii] = (max(L)/sum_L mediated over configurations) = <max(L)/sum_L> (eq. (2.49), page 51, Leone thesis)
	-> err_rc[ii] = standard deviation of the mean
	*/
	
	string name_image = "results/rc_mpi" + mpi + ".png";
	string name_output_file = "results/rc_mpi" + mpi + ".txt";

	string title = "#LT r_{c} #GT vs temperature:";

	string y_name = "#LT r_{c} #GT";

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

		mean_rc.push_back(0);
		err_rc.push_back(0);

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

		//n_skip[ii] = 0;//MOMENTANEO

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

		//for each configurazion, we compute rc and then we do the mean:

		//QUIII
		double mon_tipology, l_value;
		int discard1, discard2, discard3, discard4;
		int n_conf = 0;
		bool flag_new_conf = 1;
		double max_l = 0, sum_l = 0, rc_value;

		while (getline(input_file, line)) {
			istringstream iss(line);
			if (iss >> discard1 >> discard2 >> discard3 >> discard4 >> mon_tipology >> l_value) {

				if (mon_tipology == 0) {
					if (flag_new_conf) {
						++n_conf;
						max_l = 0;
						sum_l = 0;
						flag_new_conf = 0;
					}
					if (l_value > max_l) {
						max_l = l_value;
					}
					sum_l += l_value;
				}
				else {
					if (!flag_new_conf) {
						if (sum_l == 0) {
							cout << "sum_l = 0: so conf " << conf_id << " for" << (ii + 1) << "-th temperature T = " << temp[ii] << endl;
							continue;
						}
						rc_value = max_l / sum_l;
						mean_rc[ii] += rc_value;
						err_rc[ii] += (rc_value * rc_value);
						flag_new_conf = 1;
					}
				}
			}
			else {
				cerr << "Skipped line (badly formatted) from" << (ii + 1) << "-th input file: " << line << endl;
			}
		}
		if ((!flag_new_conf) && (sum_l > 0)) {
			if (sum_l == 0) {
				cout << "sum_l = 0: so conf " << conf_id << " for" << (ii + 1) << "-th temperature T = " << temp[ii] << endl;
			}
			else {
				rc_value = max_l / sum_l;
				mean_rc[ii] += rc_value;
				err_rc[ii] += (rc_value * rc_value);
			}
		}

		if (n_conf == 0) {
			cout << "No configuration considerd for " << to_string(ii + 1) << "-th temperature: T = " << temp[ii] << endl;
			continue;
		}

		mean_rc[ii] /= (double)n_conf;
		err_rc[ii] /= (double)n_conf;

		err_rc[ii] -= mean_rc[ii] * mean_rc[ii];

		if (n_conf == 1) {
			cout << "No sufficient configurations to evaluete error for configuration of " << to_string(ii + 1) << "-th temperature: T = " << temp[ii] << endl;
			continue;
		}

		err_rc[ii] /= (double)(n_conf - 1);

		if (err_rc[ii] < 0) {
			cout << "Negative variance for configuration of " << to_string(ii + 1) << "-th temperature: T = " << temp[ii] << endl;
			continue;
		}

		if (isnan(err_rc[ii])) {
			cout << "Nan variance for configuration of " << to_string(ii + 1) << "-th temperature: T = " << temp[ii] << endl;
			continue;
		}

		err_rc[ii] = sqrt(err_rc[ii]);

		input_file.close();

		cout << mean_rc[ii] << endl;
		cout << err_rc[ii] << endl;
		cout << n_conf << endl;
	}

	//images:
	
	if (err_rc.size() > 1) {

		plot_points_errors(
			temp,
			mean_rc,
			err_rc,
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

	output_file << "#T[MeV] \t <r_c>[adimensional] \t  \t error(r_c):" << endl;

	output_file << setprecision(numeric_limits<double>::max_digits10);

	for (int ii = 0; ii < mean_rc.size(); ++ii) {
		output_file << temp[ii] << "\t" << mean_rc[ii] << "\t" << err_rc[ii] << endl;
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