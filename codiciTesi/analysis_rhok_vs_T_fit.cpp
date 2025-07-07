/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                analysis_rhok_vs_T_fit.cpp:               ****
****Best fit of output .txt results of analysis_rhok_vs_T.cpp.****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

//GLI ERRORI NULLI MI DANNO PROBLEMI DI CHI2... 
// NON SONO NORMALI, PENSO CHE RISOLVANO QUANDO USERò I DATI GIUSTI.
// ALTRIMENTI VA COMUNQUE MESSO A MANO UN VALORE MINIMO PER L'ERRORE.
// PENSAVO UNA COSA DEL TIPO: min(err != 0)/100.0;

//-----------------------------------------------------------------
//HEADERS && LIBRARY:

#include "../library.h"
#include "../root_include.h"

#define CHOOSE_FIT_FUNCTION 2
/*
 -> 0 for rho_k/rho_1 = exp(-par[0]*(x-1))/pow(x,par[1])
 -> 1 for rho_k/rho_1 = exp(-par[0]*(x-1))/pow(x, 2.5)
 -> 2 for rho_k/rho_1 = x^(-par[0])
*/
const bool bool_choose_at_eye = 0; //0 if you want an automatic set of parameters, 1 if you want to impose them by hand;
// -> if 1, modify the corrisponding if condition in par_estimate function to choose parameters;

const bool only_one_graph = 0;//to choose to focus only on a single graph
const int index_graph = 1; //index of the graph to focus on if only_one_graph = 1.

//-----------------------------------------------------------------
//ROOT MACRO TO DO FIT AND GRAPH:

void fit_plot_points_errors(
	const vector<double>& x,
	const vector<double>& y,
	const vector<double>& y_err,
	const string name_image,
	const string title,
	const string y_name,
	const double pos_title,
	const double pos_y,
	const double heigh_y,
	const int n_par_fit,
	string name_out_file
);

void silly_plot(
	const vector<double>& x,
	const vector<double>& y,
	const vector<double>& y_err,
	const string name_image,
	const string title,
	const string y_name,
	const double pos_title,
	const double pos_y,
	const double heigh_y,
	const int n_par_fit
);


//-----------------------------------------------------------------
//FIT && ESTIMATE OF PARAMETERS FUNCTIONS: 

#if CHOOSE_FIT_FUNCTION == 0

	double fit_function(
		double* x, //number of windings
		double* p //parameters
	) {
		//if ((p[0] < 0) || (p[1] < 0)) return 1e100;
		return exp(-p[0] * (x[0] - 1)) / pow(x[0], p[1]);
	};
	
	const int n_par_fit = 2;
	
	bool par_estimate(
		//return 0 if success, 1 if not;
		const vector<double>& x, //temperatures
		const vector<double>& y, //observables
		vector<double>& p //parameters
	) {
	
		p[0] = 1;
		p[1] = 1;

		if (bool_choose_at_eye || (x.size() < n_par_fit) || (x.size() <= 3)) {//TO CHANGE THE FOLLOWING FOR "CHOOSE BY EYE" SETTING
			if (x.size() < n_par_fit) {
				cerr << "Not enough points for estimate." << endl;
				return 1;
			}
			cout << "Parameters used for fit: p0 = " << p[0]
				<< ", p1 = " << p[1] << endl;
			return 0;
		}

		//By a linear fit: (y = exp(-p[0] * (x - 1)) / pow(x, p[1]) <-> y_log = log(y) =  - p[0] * (x-1) - p[1] log(x) 
		// -> p[0] *a0 + p[1] * a1 = c
		// with a0 = 1 - x, a1 = - log(x), c = log(y) -> we can make linear fit with only p[0] and p[1] unknown)

		vector <double> a0, a1, c;

		TLinearFitter fitter(2, "x0 ++ x1");

		for (int ii = 0; ii < x.size(); ++ii) {
			if ((y[ii] > 0) && (x[ii] > 0)) {
				a0.push_back(1 - x[ii]);
				a1.push_back(-log(x[ii]));
				c.push_back(log(y[ii]));

				//cout << "(a0, a1, c) = (" << a0.back() << ", " << a1.back() << ", " << c.back() << ")" << endl;

				double vars[] = {a0.back(), a1.back() }; // being c[ii] = x0*a0[ii]+x1*a1[ii];
				fitter.AddPoint(vars, c.back()); // being c[ii] = x0*a0[ii]+x1*a1[ii];
			}
			else {
				cout << "zero y in par_estimate: not used in par estimate." << endl;
			}
			//cout << ii << ":\t" << x[ii] << "\t" << y[ii] << endl;
		}

		//cout << "N punti inseriti nel fit: " << a0.size() << endl;

		if (a0.size() < n_par_fit) {//TO CHANGE THE FOLLOWING FOR "CHOOSE BY EYE" SETTING
			if (y[2] != 0) {
				p[0] = log(y[1] / y[2]);
			}

			if (!bool_choose_at_eye) {
				p[1] = p[0] * log(2.0 / 3.0);
			}
		}
		else {
			fitter.Eval();
			
			p[0] = fitter.GetParameter(0);
			p[1] = fitter.GetParameter(1);

		}

		cout << "Parameters used for fit: p0 = " << p[0]
			<< ", p1 = " << p[1] << endl;

		return 0;
	}

#elif CHOOSE_FIT_FUNCTION == 1 

	double fit_function(
		double* x, //number of windings
		double* p //parameters
	) {
		return exp(-p[0] * (x[0] - 1)) / pow(x[0], 2.5);
	};
	
	const int n_par_fit = 1;
	
	bool par_estimate(
		//return 0 if success, 1 if not;
		const vector<double>& x, //temperatures
		const vector<double>& y, //observables
		vector<double>& p //parameters
	) {
	
		p[0] = 1;
	
		if (bool_choose_at_eye || (x.size() < n_par_fit)) {//TO CHANGE THE FOLLOWING FOR "CHOOSE BY EYE" SETTING
			if (x.size() < n_par_fit) {
				cerr << "Not enough points for estimate." << endl;
				return 1;
			}
			cout << "Parameters used for fit: p0 = " << p[0] << endl;
			return 0;
		}
	
		//By a linear fit: (y = exp(-p[0] * (x - 1)) / pow(x, 2.5) <-> z_log = log(y) + 2.5 log(x) = - p[0] * (x-1)  
		// -> p[0] *a0 = c
		// with a0 = 1 - x, c = log(y) + 2.5 log(x) -> we can make linear fit with only p[0] and p[1] unknown)

		vector <double> a0, c;

		TLinearFitter fitter(1, "x0");

		for (int ii = 0; ii < x.size(); ++ii) {
			if ((y[ii] > 0) && (x[ii] > 0)) {
				a0.push_back(1 - x[ii]);
				c.push_back(log(y[ii]) + 2.5 * log(x[ii]));

				//cout << "(a0, c) = (" << a0.back() << ", " << c.back() << ")" << endl;

				double vars[] = { a0.back()}; // being c[ii] = x0*a0[ii];
				fitter.AddPoint(vars, c.back()); // being c[ii] = x0*a0[ii];
			}
			else {
				cout << "zero y in par_estimate: not used in par estimate." << endl;
			}
			//cout << ii << ":\t" << x[ii] << "\t" << y[ii] << endl;
		}

		//cout << "N punti inseriti nel fit: " << a0.size() << endl;

		if (a0.size() < n_par_fit) {//TO CHANGE THE FOLLOWING FOR "CHOOSE BY EYE" SETTING
			p[0] = 0;

			for (int ii = 0; ii < (x.size() - 1); ++ii) {
				p[0] += log(y[ii] * pow(x[ii], 2.5)) / (1 - x[ii]);
			}

			p[0] /= (x.size() - 1);
		}
		else {
			fitter.Eval();

			p[0] = fitter.GetParameter(0);

		}

		cout << "Parameters used for fit: p0 = " << p[0] << endl;
	
		/*
		My considerations:
			y[ii] = exp(-p[0] * (x[ii] - 1)) / pow(x[ii], 2.5)
			-> p[0] = log( y[ii]*pow(x[ii], 2.5) ) /(1 - x[ii])
		*/
	
		return 0;
	}

#elif CHOOSE_FIT_FUNCTION == 2

	double fit_function(
		double* x, //number of windings
		double* p //parameters
	) {
		return pow(x[0], -p[0]);
	};
	
	const int n_par_fit = 1;
	
	bool par_estimate(
		//return 0 if success, 1 if not;
		const vector<double>& x, //temperatures
		const vector<double>& y, //observables
		vector<double>& p //parameters
	) {
	
		p[0] = 1;
	
		if (bool_choose_at_eye || (x.size() < n_par_fit)) {//TO CHANGE THE FOLLOWING FOR "CHOOSE BY EYE" SETTING
			if (x.size() < n_par_fit) {
				cerr << "Not enough points for estimate." << endl;
				return 1;
			}
			cout << "Parameters used for fit: p0 = " << p[0] << endl;
			return 0;
		}
	
		//By a linear fit: (y = pow(x, -p[0]) <-> y_log = log(y) = -p[0] * log(x)
		// -> p[0] *a0 = c
		// with a0 = -log(x), c = log(y) -> we can make linear fit with only p[0] and p[1] unknown)

		vector <double> a0, c;

		TLinearFitter fitter(1, "x0");

		for (int ii = 0; ii < x.size(); ++ii) {
			if ((y[ii] > 0) && (x[ii] > 0)) {
				a0.push_back(-log(x[ii]));
				c.push_back(log(y[ii]));

				//cout << "(a0, c) = (" << a0.back() << ", " << c.back() << ")" << endl;

				double vars[] = { a0.back() }; // being c[ii] = x0*a0[ii];
				fitter.AddPoint(vars, c.back()); // being c[ii] = x0*a0[ii];
			}
			else {
				cout << "zero y in par_estimate: not used in par estimate." << endl;
			}
			//cout << ii << ":\t" << x[ii] << "\t" << y[ii] << endl;
		}

		//cout << "N punti inseriti nel fit: " << a0.size() << endl;

		if (a0.size() < n_par_fit) {//TO CHANGE THE FOLLOWING FOR "CHOOSE BY EYE" SETTING
			p[0] = 0;

			for (int ii = 0; ii < (x.size() - 1); ++ii) {
				p[0] -= log(y[ii] / x[ii]);
			}

			p[0] /= (x.size() - 1);

			cout << "Parameters used for fit: p0 = " << p[0] << endl;

			/*
			My considerations:
				y[ii] = pow(x[ii], 1-par[0])
				-> log(y[ii]) = -par[0] * log(x[ii])
				--> par[0] = -log(y[ii]/x[ii])
			*/
		}
		else {
			fitter.Eval();

			p[0] = fitter.GetParameter(0);

		}

		return 0;
	}

#endif

//-----------------------------------------------------------------
//FUNCTION DECLARATIONS:

double chi2_reduced_estimate(
	const vector<double>& x,
	const vector<double>& y,
	const vector<double>& y_err,
	vector<double>& par,
	ofstream& output_file
);

void read_name_files(
	vector <string>& files,
	const string name_file_list
);

//-----------------------------------------------------------------
//MAIN:
int main() {

	string directory = "19_05_2025/rhok_vs_k/";
	string name_file_list = "file_list.txt";
	int skiplines_input_file = 1;

	vector <string> files;
	read_name_files(files, directory + name_file_list);

	vector <string> y_name = {//BE CAREFUL TO CHOOSE WELL;
		"#LT #rho_{k} #GT / #LT #rho_{1} #GT",
		"#LT #rho_{k} #GT / #LT #rho_{1} #GT",
		"#LT #rho_{k} #GT / #LT #rho_{1} #GT",
		"#LT #rho_{k} #GT / #LT #rho_{1} #GT",
		"#LT #rho_{k} #GT / #LT #rho_{1} #GT",
		"#LT #rho_{k} #GT / #LT #rho_{1} #GT",
		"#LT #rho_{k} #GT / #LT #rho_{1} #GT",
		"#LT #rho_{k} #GT / #LT #rho_{1} #GT",
		"#LT #rho_{k} #GT / #LT #rho_{1} #GT",
		"#LT #rho_{k} #GT / #LT #rho_{1} #GT",
		"#LT #rho_{k} #GT / #LT #rho_{1} #GT",
		"#LT #rho_{k} #GT / #LT #rho_{1} #GT",
		"#LT #rho_{k} #GT / #LT #rho_{1} #GT",
		"#LT #rho_{k} #GT / #LT #rho_{1} #GT",
		"#LT #rho_{k} #GT / #LT #rho_{1} #GT",
		"#LT #rho_{k} #GT / #LT #rho_{1} #GT",
		"#LT #rho_{k} / #rho_{1} #GT",
		"#LT #rho_{k} / #rho_{1} #GT",
		"#LT #rho_{k} / #rho_{1} #GT",
		"#LT #rho_{k} / #rho_{1} #GT",
		"#LT #rho_{k} / #rho_{1} #GT",
		"#LT #rho_{k} / #rho_{1} #GT",
		"#LT #rho_{k} / #rho_{1} #GT",
		"#LT #rho_{k} / #rho_{1} #GT",
		"#LT #rho_{k} / #rho_{1} #GT",
		"#LT #rho_{k} / #rho_{1} #GT",
		"#LT #rho_{k} / #rho_{1} #GT",
		"#LT #rho_{k} / #rho_{1} #GT",
		"#LT #rho_{k} / #rho_{1} #GT",
		"#LT #rho_{k} / #rho_{1} #GT",
		"#LT #rho_{k} / #rho_{1} #GT",
		"#LT #rho_{k} / #rho_{1} #GT"
	};

	vector <double> pos_title = {//BE CAREFUL TO CHOOSE WELL;
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5,
		0.5
	};

	vector <double> pos_y = {//BE CAREFUL TO CHOOSE WELL;
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020,
		0.020
	};

	vector <double> heigh_y = {//BE CAREFUL TO CHOOSE WELL;
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45,
		0.45
	};

	string line;
	string name_image, base_name, title;
	size_t pos;

	for (int ii = 0; ii < files.size(); ++ii) {
		
		if (only_one_graph) {
			if (ii != index_graph) {
				continue;
			}
		}

		vector <double> k, nk, err_nk;

		ifstream input_file;
		input_file.open(directory + files[ii]);
		if (!input_file) {
			cerr << directory + files[ii] << endl;
			cerr << "Error opening " << ii <<"-th input file." << endl;
			return 1;
		}

		for (int jj = 0; jj < skiplines_input_file; ++jj) {
			if (!getline(input_file, line)) {
				cerr << "Error: there are less than " << skiplines_input_file << " lines in the input " << ii << "-th file." << endl;
				return 1;
			}
		}

		while (getline(input_file, line)) {
			if (line.empty()) {
				cerr << "Skipped blank/whitespace-only line in file: " << name_file_list << endl;
				continue;
			}
			istringstream iss(line);
			double k_value, nk_value, err_nk_value;
			if (iss >> k_value >> nk_value >> err_nk_value) {
				k.push_back(k_value);
				nk.push_back(nk_value);
				err_nk.push_back(err_nk_value);
			}
			else {
				cerr << "Poorly formatted line: " << line << endl;
			}
		}

		input_file.close();

		if (nk.empty() || nk.size() <= n_par_fit) {
			cout << "No sufficient points to fit for " << files[ii] << endl;
		}
		else {
			pos = files[ii].find_last_of(".");
			if (pos != std::string::npos) {
				base_name = files[ii].substr(0, pos); //I remove extension using substr
			}

			name_image = "results/FIT_function_" + to_string(CHOOSE_FIT_FUNCTION) + "_" + base_name + ".png"; //<rho_k/rho_1>
			string name_out_file = "results/FIT_function_" + to_string(CHOOSE_FIT_FUNCTION) + "_" + base_name + ".txt";
			
			title = "Fit result:";

			for (int kk = 0; kk < nk.size(); ++kk) {
				if (nk[kk] == 0) {
					cout << "Found 0-error: I remove it: doesn't make sense -> it is associated to the of very little statistics" << endl;
					nk.erase(nk.begin() + kk);
					err_nk.erase(err_nk.begin() + kk);
					k.erase(k.begin() + kk);
					kk--;
				}
			}

			silly_plot(k,
				nk,
				err_nk, 
				"results/silly_plot_"+ to_string(ii+1) + ".png",
				title,
				y_name[ii],
				pos_title[ii],
				pos_y[ii],
				heigh_y[ii],
				n_par_fit
			);


			fit_plot_points_errors(
				k,
				nk,
				err_nk,
				name_image,
				title,
				y_name[ii],
				pos_title[ii],
				pos_y[ii],
				heigh_y[ii],
				n_par_fit,
				name_out_file
			);
		}

		k.clear();
		nk.clear();
		err_nk.clear();
	}
	return 0;
}

//-----------------------------------------------------------------
//SOME FUNCTION DEFINITIONS:

void read_name_files(
	vector <string> & files,
	const string name_file_list
) {
	ifstream file_list;
	file_list.open(name_file_list);
	if (!file_list) {
		cerr << "Error opening file-list file" << endl;
		return;
	}

	string line;
	while (getline(file_list, line)) {
		if (line.empty()) {
			cerr << "Skipped blank/whitespace-only line in file: " << name_file_list << endl;
			continue;
		}
		istringstream iss(line);
		string name_file;
		if (iss >> name_file) {
			files.push_back(name_file);
		}
		else {
			cerr << "Poorly formatted line: " << line << endl;
		}
	}
	file_list.close();
}

double chi2_reduced_estimate(
	const vector<double>& x,
	const vector<double>& y,
	const vector<double>& y_err,
	vector<double>& par,
	ofstream& output_file
) {
	double chi2r = 0, res;
	for (int ii = 0; ii < x.size(); ++ii) {
		double xx[1] = { x[ii] };
		res = (y[ii] - fit_function(xx, par.data())) / y_err[ii];
		res *= res;
		chi2r += res;
		//cout << "x["<< ii<< "]:" << x[ii] << endl;
		//cout << "y[" << ii << "]:" << y[ii] << endl;
		//cout << "y_err[" << ii << "]:" << y_err[ii] << endl;
	}
	cout << "CHI2 = " << chi2r << endl;
	output_file << "CHI2 = " << chi2r << endl;
	chi2r /= (x.size() - par.size());
	cout << "CHI2 reduced = " << chi2r << endl;
	output_file << "CHI2 reduced = " << chi2r << endl;
	return chi2r;
}

//-----------------------------------------------------------------
//ROOT MACRO TO GRAPH:

void fit_plot_points_errors(
	const vector<double>& x,
	const vector<double>& y,
	const vector<double>& y_err,
	const string name_image,
	const string title,
	const string y_name,
	const double pos_title,
	const double pos_y,
	const double heigh_y,
	const int n_par_fit,
	string name_out_file
) {

	if (x.size() != y.size() || y.size() != y_err.size()) {
		cerr << "Error: mismatched vector sizes in plot_points_errors()." << endl;
		return;
	}

	vector <double> par(n_par_fit, 0);

	ofstream output_file;
	output_file.open(name_out_file);
	if (!output_file) {
		cerr << "Error opening output file: " << name_out_file  << endl;
		return;
	}

	output_file << setprecision(numeric_limits<double>::max_digits10);
	output_file << "FIT RESULTS:" << endl;
	output_file << "CHOOSE_FIT_FUNCTION \t" << CHOOSE_FIT_FUNCTION << endl;
	output_file << "array lenght \t" << x.size() << endl;
	output_file << "INITIAL PARAMETERS:" << endl;
	output_file << "\t -> number of initial parameters \t" << n_par_fit << endl;
	output_file << "\t -> number of degrees of freedom \t" << (x.size() - n_par_fit) << endl;

	cout << "array lenght \t" << x.size() << endl;
	cout << "INITIAL PARAMETERS:" << endl;
	cout << "\t -> number of initial parameters \t" << n_par_fit << endl;
	cout << "\t -> number of degrees of freedom \t" << (x.size() - n_par_fit) << endl;

	//CANVAS creation: (to draw the graph)
	TCanvas* canvas = new TCanvas("canvas", "Canvas for Drawing Points", 900, 600);
	canvas->SetGrid();//to set grid

	// 1. Graph with only error bars:
	TGraphErrors* g_errors = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, y_err.data());

	//Using std::max_element: (useful to create my histogram object)
	//	-> std::max_element takes two iterators that define the range of the array to operate on. 
	//	-> returns an iterator that points to the maximum element found.
	//The usage of std::min_element is analogous but with for minimum element.
	auto min_x = *min_element(x.begin(), x.end());
	auto max_x = *max_element(x.begin(), x.end());
	auto min_y = *min_element(y.begin(), y.end());
	auto max_y = *max_element(y.begin(), y.end());

	g_errors->SetLineColor(kBlack);//color of error bars
	g_errors->SetMarkerColor(kBlack);
	g_errors->SetMarkerStyle(20);
	g_errors->SetTitle("");
	g_errors->GetXaxis()->SetLimits(min_x, max_x);
	g_errors->GetYaxis()->SetRangeUser(min_y - 0.01 * fabs(min_y), max_y + 0.01 * fabs(max_y));
	g_errors->Draw("AP");

	// 2. Fit of points:

	gStyle->SetOptFit(1111);

	double fit_min = 0;
	double fit_max = max_x;

	TF1* p_plot = new TF1("fit_function", fit_function, fit_min, fit_max, n_par_fit);

	par_estimate(x, y, par);

	for (int ii = 0; ii < par.size(); ++ii) {
		p_plot->SetParameter(ii, par[ii]); //setting the ii-th parameter of the function
		output_file << "\t -> inital estimate of par[" << ii << "] \t" << par[ii] << endl;
		p_plot->SetParLimits(ii, 0, 1e3);    // par[ii]: solo positivi
	}

	p_plot->GetXaxis()->SetRangeUser(min_x, max_x);
	p_plot->GetYaxis()->SetRangeUser(0, g_errors->GetMaximum() * 1.01);


	//FIT:
	gStyle->SetOptFit(0);// (1111);
	g_errors->Fit(p_plot, "SW"); // S = return result, W = use weights (i.e., y errors)

	TVirtualFitter* fit = TVirtualFitter::GetFitter(); //to get fit results;
	//ATTENTA!! Spesso, ROOT stampa "Chi2" come la media dei residui quadratici(RMS²), NON come la vera somma pesata!
	//Questo è uno di questi casi, dei valori ci si può fidare se il chi2 stimato con l'apposita funzione nel seguito è ragionevole.

	for (int ii = 0; ii < par.size(); ++ii) {
		par[ii] = p_plot->GetParameter(ii);
		output_file << "\t -> FINAL estimate of par[" << ii << "] \t" << par[ii] << "+-" << sqrt(fit->GetCovarianceMatrixElement(ii, ii)) << endl;
	}

	chi2_reduced_estimate(
		x,
		y,
		y_err,
		par,
		output_file
	);

	double cov[10][10];

	output_file << endl;
	cout << endl;
	output_file << "COV_MATRIX" << endl;
	cout << "COV_MATRIX:" << endl;
	for (int ii = 0; ii < n_par_fit; ++ii)
	{
		for (int jj = 0; jj < n_par_fit; ++jj)
		{
			cov[ii][jj] = fit->GetCovarianceMatrixElement(ii, jj);
			cout << cov[ii][jj] << "\t";
			output_file << cov[ii][jj] << "\t";
		}

		cout << endl;
		output_file << endl;

	}

	cout << endl;
	output_file << endl;

	cout << endl;
	output_file << endl;
	cout << "CORR_MATRIX:" << endl;
	output_file << "CORR_MATRIX:" << endl;
	for (int ii = 0; ii < n_par_fit; ++ii)
	{
		for (int jj = 0; jj < n_par_fit; ++jj)
		{
			cout << cov[ii][jj] / sqrt(cov[ii][ii] * cov[jj][jj]) << "\t";
			output_file << cov[ii][jj] / sqrt(cov[ii][ii] * cov[jj][jj]) << "\t";
		}
		output_file << endl;
		cout << endl;

	}

	cout << endl;
	output_file << endl;


	p_plot->SetLineColor(kRed);
	p_plot->Draw("SAME");

	gPad->Update();

	//LATEX: I add LaTeX titles and axis labels:
	TLatex latex;
	latex.SetNDC(); //sets the use of Normalized Device Coordinates (NDC).
	latex.SetTextSize(0.05); //changes text size for title
	latex.DrawLatex(pos_title, 0.92, title.c_str());
	latex.SetTextSize(0.04); //changes text size for axis labels
	latex.DrawLatex(0.45, 0.03, "T[MeV]");
	latex.SetTextAngle(90);
	latex.DrawLatex(pos_y, heigh_y, (y_name).c_str());


	//SAVE: I save the canvas as an image
	canvas->SaveAs(name_image.c_str());

	//to save also in vectorial pdf form:
	canvas->SaveAs((name_image.substr(0, name_image.find_last_of(".")) + ".pdf").c_str());

	//DELETE:
	delete p_plot;
	delete g_errors;
	delete canvas;

	output_file.close();
}

void silly_plot(
	const vector<double>& x,
	const vector<double>& y,
	const vector<double>& y_err,
	const string name_image,
	const string title,
	const string y_name,
	const double pos_title,
	const double pos_y,
	const double heigh_y,
	const int n_par_fit
) {
	if (x.size() != y.size() || y.size() != y_err.size()) {
		cerr << "Error: mismatched vector sizes in plot_points_errors()." << endl;
		return;
	}

	vector <double> par(n_par_fit, 0);
	vector <double> y_attempt;

	//CANVAS creation: (to draw the graph)
	TCanvas* canvas = new TCanvas("canvas", "Canvas for Drawing Points", 900, 600);
	canvas->SetGrid();//to set grid

	// 1. Graph with only error bars:
	TGraphErrors* g_errors = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, y_err.data());

	//Using std::max_element: (useful to create my histogram object)
	//	-> std::max_element takes two iterators that define the range of the array to operate on. 
	//	-> returns an iterator that points to the maximum element found.
	//The usage of std::min_element is analogous but with for minimum element.
	auto min_x = *min_element(x.begin(), x.end());
	auto max_x = *max_element(x.begin(), x.end());
	auto min_y = *min_element(y.begin(), y.end());
	auto max_y = *max_element(y.begin(), y.end());

	g_errors->SetLineColor(kBlack);//color of error bars
	g_errors->SetMarkerColor(kBlack);
	g_errors->SetMarkerStyle(20);
	g_errors->SetTitle("");
	g_errors->GetXaxis()->SetLimits(min_x, max_x);
	g_errors->GetYaxis()->SetRangeUser(min_y - 0.01 * fabs(min_y), max_y + 0.01 * fabs(max_y));
	g_errors->Draw("AP");

	// 2. function plot:

	TF1* f1 = new TF1("f", fit_function, min_x, max_x, n_par_fit);

	par_estimate(x, y, par);

	for (int ii = 0; ii < par.size(); ii++) {
		f1->SetParameter(ii, par[ii]); //setting the ii-th parameter of the function
	}

	f1->SetLineColor(kRed);
	f1->Draw("SAME");

	gPad->Update();

	//LATEX: I add LaTeX titles and axis labels:
	TLatex latex;
	latex.SetNDC(); //sets the use of Normalized Device Coordinates (NDC).
	latex.SetTextSize(0.05); //changes text size for title
	latex.DrawLatex(pos_title, 0.92, title.c_str());
	latex.SetTextSize(0.04); //changes text size for axis labels
	latex.DrawLatex(0.45, 0.03, "T[MeV]");
	latex.SetTextAngle(90);
	latex.DrawLatex(pos_y, heigh_y, y_name.c_str());


	//SAVE: I save the canvas as an image
	canvas->SaveAs(name_image.c_str());

	//to save also in vectorial pdf form:
	canvas->SaveAs((name_image.substr(0, name_image.find_last_of(".")) + ".pdf").c_str());

	//DELETE:
	delete f1;
	delete g_errors;
	delete canvas;
}