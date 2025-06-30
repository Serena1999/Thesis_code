/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                  analysis_rhok_vs_T.cpp:                 ****
**** Starting from dedicated files ("corrected_mon.dat"), the ****
**** number nk of monopoles is calculated for each number k of****
****  windings, for each gauge configuration. Then, averaging ****
****  over all configurations, the mean of nk (at a given k), ****
****    the mean of the ratio between nk and n1, the ratio    ****
****  between the means of nk and n1 are estimated. Two plots ****
****   are returned for the last 2 quantities mentioned as a  ****
****                      function of k.                      ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

//-----------------------------------------------------------------
//HEADERS && LIBRARY:

#include "../library.h"
#include "../root_include.h"

//IMPORTANT CONSIDERATIONS:
// -> BLOCKING SEEMS UNECESSARY FROM AUTOCORR_CHECK_RHOK.CPP (TO VERIFY BETTER)... SO FOR NOW I IMPLEMENT THE INDEPENDENT ESTIMATE OF VARIANCE FOR <nk/n1> 


//TODO:
// MANDA DI NUOVO UNO DEI FILE USCITI STRANI DI MON, SALVANDO IL MON VECCHIO IN QUALCHE MODO, PER VEDERE SE CAMBIA QUALCOSA: 
// ---> FORSE CE Dà DIVIDERE QUESTA FASE IN PIù PARTI, VERIFICA CIò.
//DEVI VEDERE SE IL BLOCKING AFFLIGGE LE MISURE SUI MONOPOLI... BASTA FARE AUTOCORR_CHECK APPOSITO PER IL SOLO NK: 
// --> fai un file che ti calcolo nk vs k + uno successivo che fa come autocorr_check; 
//PER L'ULTIMO ARRAY, IL CALCOLO DELL'ERRORE NON è BANALE... ESSENDO <>/<>, CREDO DI DOVER RICORRERE AL BOOTSTRAP. TODO
//POI FATTI RESTITUIRE ANCHE GLI ERRORI SUI FILE DI OUTPUT;
//POI FAI ANCHE LINEE A PIù T SOVRAPPOSTE, CON FIT -> 4*

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

void read_file_list_DIRECTORIES_THERM(//SCRIVI SOTTO COME IMPLEMENTARLA, PRENDI COME MINIMO DI TERMALIZZAZIONE IL MASSIMO FRA GLI N_SKIP, COSì SEI SICURA CHE VA BENE;
	const string& name_file_list,
	const int skipLines_file_list,
	vector<string>& directories,
	vector<int>& n_skip,
	int step_sample_gauge,
	int step_sample_fermion
);

double f(double rho_k, double rho_1) {
	return rho_k / rho_1;
}

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
	string name_input_file = "corrected_mon.dat";
	int skip_lines_input_file = 0;

	int n_boot_steps = 100; //for estimate errors of secondary observables with Bootstrap.
	int seed = 1; //for estimate errors of secondary observables with Bootstrap.

	vector <double> mean_rhok, mean_rhok_rho1, mean_rhok_norm;
	/*
	-> mean_rhok[ii] = mean of nk[ii] over configurations = <nk[ii]>
	-> mean_rhok_rho1[ii] = <nk[ii]/nk[0]> (0 corrispondente ad 1 wrapping)
	-> mean_rhok_norm[ii] = <nk[ii]>/<nk[0]>
	*/
	vector <double> mean0_2, mean1_2, mean2_2; //mean of squares
	vector <double> err_rhok, err_rhok_rho1, err_rhok_norm; //standard deviation from the mean

	string name_image1, name_image2, name_output_file1, name_output_file2;

	string title1 = "#LT #rho_{k} / #rho_{1} #GT vs Number of wrapping:";
	string title2 = "#LT #rho_{k} #GT / #LT #rho_{1} #GT vs Number of wrapping:";

	string y_name1 = "#LT #rho_{k} / #rho_{1} #GT";
	string y_name2 = "#LT #rho_{k} #GT / #LT #rho_{1} #GT";

	double pos_title1 = 0.5;
	double pos_title2 = 0.5;

	double pos_y1 = 0.020;
	double pos_y2 = 0.020;

	double heigh_y1 = 0.45;
	double heigh_y2 = 0.45;

	double width_canvas1 = 900;
	double width_canvas2 = 900;

	read_file_list_DIRECTORIES_THERM(
		list_file,
		skipLines_list_file,
		directories,
		n_skip,
		step_sample_gauge,
		step_sample_fermion
	);

	sample_gen sampler;
	sampler.init(seed);

	for (int ii = 0; ii < directories.size(); ++ii) {
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
		double n_accum = 0;
		int n_conf = 0;
		vector <double> value_rhok, value_rhok_rho1, value_rhok_norm, k_array;
		vector <vector<double>> matrix_value_rhok;

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
					if ((value_rhok.empty()) || (value_rhok[0] == 0)) {
						cout << "No monopoles with only one wrapping for conf_id = " << to_string(conf_tmp) << endl;
						cout << "--> so we omit this configuration." << endl;
						conf_tmp = conf_id;
					}
					else {
						++n_conf;
						for (int kk = 0; kk < value_rhok.size(); ++kk) {
							if (value_rhok_rho1.size() <= kk) {
								value_rhok_rho1.resize(kk + 1, 0.0);
							}
							value_rhok_rho1[kk] = value_rhok[kk] / (double)value_rhok[0];
							if (kk >= mean_rhok.size()) {
								mean_rhok.resize(kk + 1, 0.0);
								mean0_2.resize(kk + 1, 0.0);
								mean_rhok_rho1.resize(kk + 1, 0.0);
								mean1_2.resize(kk + 1, 0.0);
							}
							if (kk >= matrix_value_rhok.size()) {
								matrix_value_rhok.resize(kk + 1, vector<double>());
							}
							matrix_value_rhok[kk].push_back(value_rhok[kk]);
							value_rhok[kk] /= n_accum;
							mean_rhok[kk] += value_rhok[kk];
							mean0_2[kk] += (value_rhok[kk] * value_rhok[kk]);
							value_rhok_rho1[kk] /= n_accum;
							mean_rhok_rho1[kk] += value_rhok_rho1[kk];
							mean1_2[kk] += (value_rhok_rho1[kk] * value_rhok_rho1[kk]);
						}
					}
					value_rhok.clear();
					value_rhok_rho1.clear();
					n_accum = 0;
					conf_tmp = conf_id;
				}
				if (wrap_value == 0) {
					continue;
				}
				if (abs(wrap_value) > value_rhok.size()) {
					value_rhok.resize(abs(wrap_value), 0.0);
					value_rhok_rho1.resize(abs(wrap_value), 0.0);
				}
				++n_accum;
				value_rhok[abs(wrap_value) - 1] += 1;
			}
			else {
				cerr << "Skipped line (badly formatted) from" << (ii + 1) << "-th input file: " << line << endl;
			}
		}
		if (conf_id == conf_tmp) {
			if ((value_rhok.empty()) || (value_rhok[0] == 0)) {
				cout << "No monopoles with only one wrapping for conf_id = " << to_string(conf_tmp) << endl;
				cout << "--> so we omit this configuration." << endl;
			}
			else {
				++n_conf;
				for (int kk = 0; kk < value_rhok.size(); ++kk) {
					value_rhok_rho1[kk] = value_rhok[kk] / (double)value_rhok[0];
					//cout << "value_rhok_rho1.size() = " << value_rhok_rho1.size() << endl;
					if (kk >= mean_rhok.size()) {
						mean_rhok.resize(kk + 1, 0.0);
						mean0_2.resize(kk + 1, 0.0);
						mean_rhok_rho1.resize(kk + 1, 0.0);
						mean1_2.resize(kk + 1, 0.0);
					}
					if (kk >= matrix_value_rhok.size()) {
						matrix_value_rhok.resize(kk + 1, vector<double>());
					}
					matrix_value_rhok[kk].push_back(value_rhok[kk]);
					value_rhok[kk] /= n_accum;
					mean_rhok[kk] += value_rhok[kk];
					mean0_2[kk] += (value_rhok[kk] * value_rhok[kk]);
					value_rhok_rho1[kk] /= n_accum;
					mean_rhok_rho1[kk] += value_rhok_rho1[kk];
					mean1_2[kk] += (value_rhok_rho1[kk] * value_rhok_rho1[kk]);
				}
			}
			value_rhok.clear();
			value_rhok_rho1.clear();
			//if (ii == 12) {
			//	cout << endl;
			//	//cout << kk << endl;
			//	cout << "mean_rhok.size() = " << mean_rhok.size() << endl;
			//	cout << endl;
			//}
		}

		input_file.close();

		cout << mean_rhok.size() << endl;
		cout << mean_rhok_rho1.size() << endl;
		cout << n_conf << endl;

		if(mean_rhok.size() <= 1){
			cout << "No monopoles with more than 1 wrapping in " << ii << "-th iteration" << endl;
			cout << "Not considered " << ii << "-th iteration" << endl;
			continue;
		}

		if (mean_rhok[0] == 0) {
			cout << "No monopoles with 1 wrapping in " << ii << "-th iteration" << endl;
			cout << "Not considered " << ii << "-th iteration" << endl;
			continue;
		}

		for (int kk = 0; kk < mean_rhok.size(); ++kk) {
			mean_rhok[kk] /= n_conf;
			mean_rhok_rho1[kk] /= n_conf;
			mean0_2[kk] /= n_conf;
			mean1_2[kk] /= n_conf;

			err_rhok.push_back(mean0_2[kk] - mean_rhok[kk] * mean_rhok[kk]);
			err_rhok_rho1.push_back(mean1_2[kk] - mean_rhok_rho1[kk] * mean_rhok_rho1[kk]);

			err_rhok[kk] /= (n_conf - 1);
			err_rhok_rho1[kk] /= (n_conf - 1);

			if (err_rhok_rho1[kk] < 0) {
				cout << "err_rhok_rho1[" << to_string(kk) << "] < 0" << endl;
				err_rhok_rho1[kk] = 0;
			}
			else {
				err_rhok_rho1[kk] = sqrt(err_rhok_rho1[kk]);
			}
		}

		if (mean_rhok_rho1.size() > 1) {
			mean_rhok_rho1.erase(mean_rhok_rho1.begin()); //non ci interessa rho_1/rho_1 = 1
			err_rhok_rho1.erase(err_rhok_rho1.begin());
		}
		else {
			cout << "No monopoles with more than 1 wrapping in " << ii << "-th iteration" << endl;
			cout << "Not considered " << ii << "-th iteration" << endl;
			continue;
		}

		int n_not_considered = 0;

		for (int kk = 1; kk < mean_rhok.size(); ++kk) {

			if (matrix_value_rhok[kk].empty() || matrix_value_rhok[0].empty()) {
				cerr << "Warning: empty data for Bootstrap at k=" << kk << endl;
				++n_not_considered;
				continue;
			}

			vector <double> x1_vec, x2_vec;

			for (int jj = 0; jj < min(matrix_value_rhok[kk].size(), matrix_value_rhok[0].size()); ++jj) {
				if (matrix_value_rhok[0][jj] != 0) {
					x1_vec.push_back(matrix_value_rhok[kk][jj]);
					x2_vec.push_back(matrix_value_rhok[0][jj]);
				}
				else {
					cout << "Omitted a zero rho_1 due to division issue" << endl;
				}
			}

			mean_rhok_norm.push_back(0);
			err_rhok_norm.push_back(0);
			uniform_int_distribution<> dist_int(0, x1_vec.size() - 1);

			if (bootstrap_2in(n_boot_steps, sampler, dist_int, f, x1_vec, x2_vec, &mean_rhok_norm[kk - 1 - n_not_considered], &err_rhok_norm[kk - 1 - n_not_considered]))
			{
				cout << "Error in Bootstrap for (ii, kk) = (" << to_string(ii) << ", " << to_string(kk) << endl;
				return 1;
			}

			if (isnan(err_rhok_norm[kk - 1 - n_not_considered])) {
				cout << "Encontered NaN in Bootstrap" << endl;
				err_rhok_norm.resize(kk - 1 - n_not_considered, 0.0);;
				break;
			}
			else if (err_rhok_norm[kk - 1 - n_not_considered] >= 0) {
				err_rhok_norm[kk - 1 - n_not_considered] = sqrt(err_rhok_norm[kk - 1 - n_not_considered]);
			}
			else {
				cout << "Encontered negative value in Bootstrap" << endl;
				err_rhok_norm.resize(kk - 1 - n_not_considered, 0.0);;
				break;
			}
			k_array.push_back(kk + 1);
			x1_vec.clear();
			x2_vec.clear();
		}

		//images:
		name_image1 = "results/rhok_rho1_mpi" + mpi + "_" + to_string(ii + 1) + "thTemperature.png"; //<rho_k/rho_1>
		name_image2 = "results/rhok_norm_mpi" + mpi + "_" + to_string(ii + 1) + "thTemperature.png"; //<rho_k>/<rho_1>

		if (mean_rhok_norm.size() > 1) {
			
			plot_points_errors(
				k_array,
				mean_rhok_rho1,
				err_rhok_rho1,
				name_image1,
				title1,
				y_name1,
				pos_title1,
				pos_y1,
				heigh_y1,
				width_canvas1
			);
			
			plot_points_errors(
				k_array,
				mean_rhok_norm,
				err_rhok_norm,
				name_image2,
				title2,
				y_name2,
				pos_title2,
				pos_y2,
				heigh_y2,
				width_canvas2
			);
		}
		else {
			cout << "No image generated for iteration n°" << ii << endl;
			cout << "Not enough data points to produce a plot (need at least 2)." << endl;
		}

		//output_file:
		name_output_file1 = "results/" + mpi + "_TempN" + to_string(ii + 1) + "_mean_rhok_rho1VSk.txt";
		name_output_file2 = "results/" + mpi + "_TempN" + to_string(ii + 1) + "_mean_rhok_mean_rho1VSk.txt";

		ofstream output_file1;
		output_file1.open(name_output_file1);
		if (!output_file1) {
			cerr << "Error opening first output file" << endl;
			return 1;
		}

		ofstream output_file2;
		output_file2.open(name_output_file2);
		if (!output_file2) {
			cerr << "Error opening second output file" << endl;
			return 1;
		}

		output_file1 << "#k \t <rhok/rho1> \t err(rhok/rho1):" << endl;
		output_file2 << "#k \t <rhok>/<rho1> \t err(<rhok>/<rho1>):" << endl;

		output_file1 << setprecision(numeric_limits<double>::max_digits10);
		output_file2 << setprecision(numeric_limits<double>::max_digits10);

		for (int kk = 0; kk < k_array.size(); ++kk) {
			output_file1 << k_array[kk] << "\t" << mean_rhok_rho1[kk] << "\t" << err_rhok_rho1[kk] << endl;
			output_file2 << k_array[kk] << "\t" << mean_rhok_norm[kk] << "\t" << err_rhok_norm[kk] << endl;
		}

		output_file1.close();
		output_file2.close();

		mean_rhok.clear();
		mean_rhok_norm.clear();
		mean_rhok_rho1.clear();
		k_array.clear();
		err_rhok.clear();
		err_rhok_rho1.clear();
		err_rhok_norm.clear();
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
	g_errors->GetXaxis()->SetLimits(0, max_x + 1);
	//g_errors->GetYaxis()->SetRangeUser(min_y - 0.01 * fabs(min_y), max_y + 0.01 * fabs(max_y));
	g_errors->GetYaxis()->SetRangeUser(0, max_y + 0.01 * fabs(max_y));
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
