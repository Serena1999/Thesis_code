/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                 autocorr_check_rhok.cpp:                 ****
****      Starting from dedicated files ("mon.dat"), the      ****
**** number nk of monopoles is calculated for each number k of****
**** windings, for each gauge configuration. Then, we compute ****
****  nk/n1 for each k, estimate its variance by varying the  ****
****        dimension of blocks in blocking procedure.        ****
****   In output: images and file with <nk/n1> vs dim_block.  ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

//-----------------------------------------------------------------
//HEADERS && LIBRARY:

#include "../library.h"
#include "../root_include.h"


//RIUSALO QUANDO AVRAI TOLTO CORRETTAMENTE N_THERM, CIOè QUANDO HAI TUTTE LE CONFIGURAZIONI CORRETTE.
//IMPORTANT! DA QUESTO FILE, SEMBRA CHE IL BLOCKING NON SERVA PER I VALORI RHO_K/RHO_1.

//-----------------------------------------------------------------
//ROOT MACRO TO DO FIT AND GRAPH:

void plot_points(
	vector<double>& x,
	vector<double>& y,
	string name_image,
	string title,
	string y_name,
	double pos_title,
	double pos_y,
	double heigh_y,
	double width_canvas
);

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

	double n_sub_ratio = 0.5;//=0.5*len(data) -> you can modify this number from 0 to 1;

	string mpi = "800";
	string line;

	string list_file = "11_05_2025/file_list_therm_extended.txt";
	int skipLines_list_file = 1;
	int step_sample_gauge = 1;
	int step_sample_fermion = 10;

	vector <string> directories;//LEGGILI DA FILE
	vector <int> n_skip; //LEGGILO DA FILE;
	string name_input_file = "mon.dat";
	int skip_lines_input_file = 0;

	/*
	-> mean_rhok[ii] = mean of nk[ii] over configurations = <nk[ii]>
	-> mean_rhok_rho1[ii] = <nk[ii]/nk[0]> (0 corrispondente ad 1 wrapping)
	-> mean_rhok_norm[ii] = <nk[ii]>/<nk[0]>
	*/
	vector <double> var_rhok_rho1; //standard deviation from the mean
	double mean;

	string name_image, name_output_file;

	string title = "var(#rho_{k} / #rho_{1}) vs dimension of blocks:";

	string y_name = "var(#rho_{k} / #rho_{1})";
	
	double pos_y = 0.020;
	
	double heigh_y = 0.45;
	
	double width_canvas = 900;
	
	read_file_list_DIRECTORIES_THERM(
		list_file,
		skipLines_list_file,
		directories,
		n_skip,
		step_sample_gauge,
		step_sample_fermion
	);

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
		vector <double> value_rhok, value_rhok_rho1;
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
							if (kk >= matrix_value_rhok.size()) {
								matrix_value_rhok.resize(kk + 1, vector<double>());
							}
							matrix_value_rhok[kk].push_back(value_rhok_rho1[kk]);
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
					if (value_rhok_rho1.size() <= kk) {
						value_rhok_rho1.resize(kk + 1, 0.0);
					}
					value_rhok_rho1[kk] = value_rhok[kk] / (double)value_rhok[0];
					if (kk >= matrix_value_rhok.size()) {
						matrix_value_rhok.resize(kk + 1, vector<double>());
					}
					matrix_value_rhok[kk].push_back(value_rhok_rho1[kk]);
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

		for (int kk = 0; kk < matrix_value_rhok.size(); ++kk) {

			name_output_file = "results/" + mpi + "_TempN" + to_string(ii + 1) + "_mean_rho" + to_string(kk + 1) + "_rho1VSk.txt";;

			ofstream output_file; //declaration of output file
			output_file.open(name_output_file);
			if (!output_file) {
				cout << "Error opening output file " << name_output_file << endl;
				return 1;
			}

			output_file << setprecision(numeric_limits<double>::max_digits10);
			output_file << "#dim_block \t var(rho" << to_string(kk + 1) << "/rho1):" << endl;

			int n_sub_max = floor(n_sub_ratio * matrix_value_rhok[kk].size());
			vector <double> dim_block_vec;

			if (dim_block_vec.size() != var_rhok_rho1.size()) {
				cout << "There is a conceptual errors in vector dimensiong." << endl;
				return 1;
			}


			for (int dim_block = 1; dim_block <= n_sub_max; dim_block++) {
				var_rhok_rho1.push_back(0);
				if (!blocking_faster(&mean, &var_rhok_rho1[dim_block - 1], matrix_value_rhok[kk], dim_block)) {
					var_rhok_rho1.resize(dim_block - 1, 0.0);
					break;
				}
				if (!isnan(var_rhok_rho1[dim_block - 1])) {
					output_file << dim_block << "\t" << var_rhok_rho1[dim_block - 1] << endl;
					dim_block_vec.push_back(dim_block);
				}
				else {
					var_rhok_rho1.resize(dim_block - 1, 0.0);
					break;
				}
			}

			output_file.close();

			name_image = "results/rho" + to_string(kk + 1) + "_rho1_mpi" + mpi + "_" + to_string(ii + 1) + "thTemperature.png";

			if (var_rhok_rho1.size() > 1) {
				vector <double> n_sub_d;
				copy(dim_block_vec.begin(), dim_block_vec.end(), std::back_inserter(n_sub_d));

				plot_points(
					n_sub_d,
					var_rhok_rho1,
					name_image,
					title,
					y_name,
					5.0,
					pos_y,
					heigh_y,
					width_canvas
				);
				n_sub_d.clear();
			}
			else {
				cout << "No sufficient point to graph for " << to_string(ii + 1) << " T and " << to_string(kk + 1) << " number of windings" << endl;
			}

			var_rhok_rho1.clear();
			dim_block_vec.clear();
		}

		value_rhok.clear();
		value_rhok_rho1.clear();
		matrix_value_rhok.clear();
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

void plot_points(
	vector<double>& x,
	vector<double>& y,
	string name_image,
	string title,
	string y_name,
	double pos_title,
	double pos_y,
	double heigh_y,
	double width_canvas
) {

	if (x.size() != y.size()) {
		cerr << "Error: mismatched vector sizes in plot_points_errors()." << endl;
		return;
	}

	//CANVAS creation: (to draw the graph)
	TCanvas* canvas = new TCanvas("canvas", "Canvas for Drawing Points", width_canvas, 600);
	canvas->SetGrid();//to set grid

	// 1. Graph with only error bars:
	TGraph* g = new TGraph(x.size(), x.data(), y.data());

	//Using std::max_element: (useful to create my histogram object)
	//	-> std::max_element takes two iterators that define the range of the array to operate on. 
	//	-> returns an iterator that points to the maximum element found.
	//The usage of std::min_element is analogous but with for minimum element.
	auto min_x = *min_element(x.begin(), x.end());
	auto max_x = *max_element(x.begin(), x.end());
	auto min_y = *min_element(y.begin(), y.end());
	auto max_y = *max_element(y.begin(), y.end());

	g->SetLineColor(kBlack);//color of error bars
	g->SetMarkerColor(kBlack);
	g->SetMarkerStyle(20);
	g->SetTitle("");
	g->GetXaxis()->SetLimits(min_x, max_x);
	g->GetYaxis()->SetRangeUser(min_y, max_y);
	g->Draw("AP");

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
	latex.DrawLatex(0.45, 0.03, "dim_block");
	latex.SetTextAngle(90);
	latex.DrawLatex(pos_y, heigh_y, y_name.c_str());


	//SAVE: I save the canvas as an image
	canvas->SaveAs(name_image.c_str());

	//to save also in vectorial pdf form:
	canvas->SaveAs((name_image.substr(0, name_image.find_last_of(".")) + ".pdf").c_str());

	//DELETE:
	delete g_line;
	delete g;
	delete canvas;
}
