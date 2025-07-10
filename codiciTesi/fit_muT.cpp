//-----------------------------------------------------------------
//HEADERS && LIBRARY:

#include "../library.h"
#include "../root_include.h"

const bool bool_choose_at_eye = 0; //0 if you want an automatic set of parameters, 1 if you want to impose them by hand;
// -> if 1, modify the corrisponding if condition in par_estimate function to choose parameters;

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
//FUNCTION DECLARATIONS:

double chi2_reduced_estimate(
	const vector<double>& x,
	const vector<double>& y,
	const vector<double>& y_err,
	vector<double>& par,
	ofstream& output_file
);

//-----------------------------------------------------------------
//FIT && ESTIMATE OF PARAMETERS FUNCTIONS: 

double fit_function(
	double* x, //number of windings
	double* p //parameters
) {
	return p[0] * pow(x[0] - p[1], p[2]);
};

const int n_par_fit = 3;

bool par_estimate(
	//return 0 if success, 1 if not;
	const vector<double>& x, //temperatures
	const vector<double>& y, //observables
	vector<double>& p //parameters
) {

	p[0] = 1;
	p[1] = 250;
	p[2] = 1;

	if (bool_choose_at_eye || (x.size() < n_par_fit) || (x.size() <= 3)) {//TO CHANGE THE FOLLOWING FOR "CHOOSE BY EYE" SETTING
		if (x.size() < n_par_fit) {
			cerr << "Not enough points for estimate." << endl;
			return 1;
		}
		cout << "Parameters used for fit: p0 = " << p[0]
			<< ", p1 = " << p[1]
			<< ", p2 = " << p[2] << endl;
		return 0;
	}

	//SE VEDI CHE SERVE, IMPLEMENTA ROBA PIù ARTICOLATA

	cout << "Parameters used for fit: p0 = " << p[0]
		<< ", p1 = " << p[1]
		<< ", p2 = " << p[2] << endl;

	return 0;
}

//-----------------------------------------------------------------
//MAIN:

int main(int argc, char** argv) {

	//TO CHOOSE:
	string mpi = "1500";
	string quantity = "1"; //1 for <rhok>/<rho1>, 2 for <rhok/rho1> 

	TApplication app("App", &argc, argv);
	
	string name_out_file = "results/FIT_MUT_" + quantity + "_mpi" + mpi + ".txt";
	string name_image = "results/FIT_MUT_" + quantity + "_mpi" + mpi + ".png";

	vector <double> x = {
	244.89 ,
	261.844,
	279.271,
	297.374,
	316.388,
	336.545,
	358.027,
	380.897,
	404.993,
	433.553
	};

	vector <double> y = {
		0.528143,
		1.09037 ,
		1.22604 ,
		1.69116 ,
		1.96273 ,
		2.20709 ,
		2.55058 ,
		2.77996 ,
		3.07537 ,
		3.40656
	};

	vector <double> dy = {
		0.200577 ,
		0.0871921,
		0.0542246,
		0.0411082,
		0.0579388,
		0.0670202,
		0.11059	 ,
		0.0626644,
		0.0697073,
		0.0258094
	};

	fit_plot_points_errors(
		x,
		y,
		dy,
		name_image,
		"Effective chemical potential as a function of T:",
		"#hat{#mu}",
		0.2,
		0.02,
		0.5,
		n_par_fit,
		name_out_file
	);

	return 0;
}

//-----------------------------------------------------------------
//SOME FUNCTION DEFINITIONS:

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
		cerr << "Error opening output file: " << name_out_file << endl;
		return;
	}

	output_file << setprecision(numeric_limits<double>::max_digits10);
	output_file << "FIT RESULTS:" << endl;
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
	auto min_y = *min_element(y.begin(), y.end()) - *max_element(y_err.begin(), y_err.end());
	auto max_y = *max_element(y.begin(), y.end()) + *max_element(y_err.begin(), y_err.end());

	g_errors->SetLineColor(kBlack);//color of error bars
	g_errors->SetMarkerColor(kBlack);
	g_errors->SetMarkerStyle(20);
	g_errors->SetTitle("");

	double fit_min = min_x;
	double fit_max = max_x;

	g_errors->GetXaxis()->SetLimits(min_x, max_x);
	g_errors->GetYaxis()->SetRangeUser(0, max_y + 0.01 * fabs(max_y));
	g_errors->Draw("AP");

	// 2. Fit of points:

	gStyle->SetOptFit(1111);
	TF1* p_plot = new TF1("fit_function", fit_function, fit_min, fit_max, n_par_fit);

	par_estimate(x, y, par);

	for (int ii = 0; ii < par.size(); ++ii) {
		p_plot->SetParameter(ii, par[ii]); //setting the ii-th parameter of the function
		output_file << "\t -> inital estimate of par[" << ii << "] \t" << par[ii] << endl;
		//p_plot->SetParLimits(ii, 0, 1e6);    // par[ii]: solo positivi
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

	gApplication->Run(true); // <--- TRUE = non bloccare il terminale

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


/*
VALUE USED FOR MPI = 800 MEV and <rhok>/<rho1>:
	vector <double> x = {
		243.377,
		257.347,
		271.711,
		286.507,
		301.786,
		317.618,
		334.088,
		351.294,
		369.346//,
//		388.360//,
//		404.645
	};

	vector <double> y = {
		0.284661,//INDECISA SE METTERLO
		0.398901,
		0.545283,
		0.77939,
		0.96586,
		1.09371,
		1.36387,//ho messo 1, se non va bene metti 2
		1.56162,
		1.75081//,
//		1.87469//,
//		2.13651
	};

	vector <double> dy = {
		0.108201,//INDECISA SE METTERLO
		0.193213,
		0.357746,
		0.0816677,
		0.16642,
		0.0754943,
		0.0374763, //ho messo 1, se non va bene metti 2
		0.0200966,
		0.0147706//,
//		0.0160792//,
//		0.0775104//INDECISA SE METTERLO
	};

	USED FOR MPI = 800 MEV AND <RHOK/RHO1>:

		vector <double> x = {
		243.377,
		257.347,
		271.711,
		286.507,
		301.786,
		317.618,
		334.088,
		351.294,
		369.346//,
//		388.360//,
//		404.645
	};

	vector <double> y = {
			0.277787,//INDECISA SE METTERLO
			0.392144,
			0.538699,
			0.768353,
			0.952251,
			1.08448,
			1.34192, //ho messo 1, se non va bene metti 2
			1.54398,
			1.72737//,
//			1.85577//,
//			2.1217
	};

	vector <double> dy = {
		0.117498,//INDECISA SE METTERLO
		0.171818,
		0.471515,
		0.0832165,
		0.180344,
		0.0902923,
		0.0359896, //ho messo 1, se non va bene metti 2
		0.0201032,
		0.0159643//,
//		0.0192587//,
//		0.0775332//INDECISA SE METTERLO
	};

	VALUE USED FOR MPI = 1500 MEV and <rhok>/<rho1>:
		vector <double> x = {
	244.89 ,
	261.844,
	279.271,
	297.374,
	316.388,
	336.545,
	358.027,
	380.897,
	404.993,
	433.553
	};

	vector <double> y = {
		0.543064,
		0.9274  ,
		1.24625 ,
		1.70588 ,
		1.80587 ,
		2.23249 ,
		2.57356,
		2.79353,
		3.09302,
		3.44086
	};

	vector <double> dy = {
		0.22407  ,
		0.114908 ,
		0.0321993,
		0.0421979,
		0.345458 ,
		0.0771824,
		0.111547 ,
		0.0753449,
		0.0787561,
		0.0623854
	};

	VALUE USED FOR MPI = 1500 MEV and <rhok/rho1>:


*/