/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                visual_autocorrelation.cpp:               ****
****    visualization of the correlation function C(k) vs k   ****
****        and of (the sum of C(k) up to kmax) vs kmax       **** 
****        (kmax = the maximum value of k considered).       ****
****                                                          ****
****  The estimate of the integrated time tau_int is done by  ****
****      fitting C(k) vs k with an exponential function      ****
****                  C(k) = exp(-k/tau_int).                 ****
****  Finally, I draw the horizontal line with height=tau_int ****
****         on the same graph of sum_k C(k) vs k_max,        ****
****   so it can be seen that our curve tends asymptotically  ****
****         to tau_int for a suitable sample of data.        ****
****  This is followed by the estimate of the sample variance ****
****     of the sample mean which also takes into account     ****
****             the autocorrelation of the data.             ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

//-----------------------------------------------------------------
//HEADERS && LIBRARY:

#include "../library.h"
#include "../root_include.h"

/*On the ROOT shell, we will have the output:
(for "":, therm = 2000)
****************************************
****************************************
Minimizer is Minuit2 / Migrad
Chi2                      =   0.00329755
NDf                       =            7
Edm                       =  6.04787e-08
NCalls                    =           22
p0                        =      1.40161   +/-   0.0805495
Info in <TCanvas::Print>: png file results/image_Ck_vs_k_reff_ferm_obs2239384639.png has been created
integrated autoccorrelation time = 0.713465
integrated autoccorrelation time standard deviation = 0.0410023
Info in <TCanvas::Print>: png file results/image_sumCk_vs_kmax_reff_ferm_obs2239384639.png has been created

C:\Users\User\Desktop\UNI\TESI_MAGISTRALE\dati_da_remoto\codiciTesi\out\build\x64-debug\codiciTesi\visual_autocorrelation.exe (processo 2576) terminato. Codice restituito: 0 (0x0).
where:
	-> Chi2 = chi square of the fit procedure;
	-> NDf = number of degrees of freedom;
	-> Edm = "estimated distance from minimum", which is very small here, indicating good convergence;
	-> NCalls = number of function calls to find minimum;
	-> p0 = estimated value of the parameter +/- associated standard deviation.
https://root.cern.ch/doc/master/classROOT_1_1Fit_1_1FitResult.html
https://root.cern/manual/fitting/
*/

//-----------------------------------------------------------------
//FUNCTIONS DEFINITIONS:

double f_fit(double* y, double* par) {//par[1] = 1/tau_int
	return exp(-y[0]*par[0]);
}

//-----------------------------------------------------------------
//ROOT MACRO TO DO FIT AND GRAPH:

void graph_points_fit(vector<double>& x, vector<double>& y, int len_y, string name_image, double * result, double * dresult);
void plot_points_orizline(vector<double>& x, vector<double>& y, double y0, string name_image);

//-----------------------------------------------------------------
//MAIN:

int main() {

	int skipLines = 1; //= number of lines to skip while reading input file;
	vector<double> Ck, sumCk, k;
	double k_value, Ck_value;
	double sumCk_value = 0;
	double in_tau_int, tau_int, err_in_tau_int, err_tau_int;
	double mean = 0.0, var_m = 0.0, dev_m = 0.0;
	int n_therm;
	string line;

	string name_input_file = "results/dependent_data_analysis_Ck_reff_ferm_obs2239384639.txt"; //TO CHOOSE BEFORE COMPILATION
	string name_output_file = "results/dependent_data_analysis_output_reff_ferm_obs2239384639.txt"; //TO CHOOSE BEFORE COMPILATION
	string name_image1 = "results/image_Ck_vs_k_reff_ferm_obs2239384639.png"; //TO CHOOSE BEFORE COMPILATION
	string name_image2 = "results/image_sumCk_vs_kmax_reff_ferm_obs2239384639.png"; //TO CHOOSE BEFORE COMPILATION

	ifstream input_file; //declaration of input file

	input_file.open(name_input_file);
	if (!input_file) {
		cout << "Error opening input file" << endl;
		return 1;
	}

	for (int i = 0; i < skipLines; i++) {
		if (!getline(input_file, line)) {
			cerr << "Error: there are less than " << skipLines << " lines in the file." << endl;
			return 1;
		}
	}

	while (input_file >> k_value >> Ck_value) {
		k.push_back(k_value);
		Ck.push_back(Ck_value);
		sumCk_value += Ck_value;
		sumCk.push_back(sumCk_value);
	}

	input_file.close();

	//GRAPHIC REPRESENTATION:
	graph_points_fit(k, Ck, k.size(), name_image1, &in_tau_int, &err_in_tau_int);
	
	tau_int = 1 / in_tau_int;
	err_tau_int = err_in_tau_int / (in_tau_int * in_tau_int);
	
	ofstream output_file; //declaration of input file

	output_file.open(name_output_file);
	if (!output_file) {
		cout << "Error opening output file" << endl;
		return 1;
	}
	output_file << setprecision(numeric_limits<double>::max_digits10);//This way, the output to the file will show the double values ​​with the greatest possible precision.
	output_file << "integrated autoccorrelation time \t" << tau_int << endl;
	output_file << "integrated autoccorrelation time standard deviation \t" << err_tau_int << endl;
	cout << "integrated autoccorrelation time = "<< tau_int << endl;
	cout << "integrated autoccorrelation time standard deviation = " << err_tau_int << endl;
	
	output_file.close();

	plot_points_orizline(k, sumCk, tau_int, name_image2);

	return 0;
}

//-----------------------------------------------------------------
//ROOT MACRO TO DO FIT AND GRAPH:

void graph_points_fit(vector<double>& x, vector<double>& y, int len_y, string name_image, double* result, double* dresult) {

	string func_name = "p_plot_" + name_image;

	//CANVAS creation: (to draw the graph)
	TCanvas* canvas = new TCanvas("canvas", "Canvas for Drawing Points", 800, 600);
	canvas->SetGrid();//to set grid

	TGraph* g = new TGraph(y.size(), x.data(), y.data());

	//Using std::max_element: (useful to create my histogram object)
	//	-> std::max_element takes two iterators that define the range of the array to operate on. 
	//	-> returns an iterator that points to the maximum element found.
	//The usage of std::min_element is analogous but with for minimum element.
	auto min_x = *min_element(x.begin(), x.end());
	auto max_x = *max_element(x.begin(), x.end());
	auto min_y = *min_element(y.begin(), y.end());
	auto max_y = *max_element(y.begin(), y.end());

	g->SetTitle(""); //to clear default title
	g->GetXaxis()->SetTitle(""); //to clear default X axis title
	g->GetYaxis()->SetTitle(""); //to clear default Y axis title

	g->SetLineColor(1);

	g->GetYaxis()->SetRangeUser(min_y - 0.01 * abs(min_y), max_y + 0.01 * abs(max_y));//to set y axis range
	g->GetXaxis()->SetRangeUser(min_x - 0.01 * abs(min_x), max_x + 0.01 * abs(max_x));//to set x axis range

	gStyle->SetOptFit(1111);

	g->Draw("ALP");

	//1D FUNCTION object creation:
	TF1* p_plot = new TF1(func_name.c_str(), f_fit, min_x, max_x, 1); //1 is the number of parameters

	p_plot->SetParameter(0, 1); //setting the first parameter of the function
	//p_plot->SetParameter(1, 100); //setting the second parameter of the function
	p_plot->GetXaxis()->SetRangeUser(min_x, max_x);
	p_plot->GetYaxis()->SetRangeUser(0, g->GetMaximum() * 1.01);


	//FIT:
	gStyle->SetOptFit(0);// (1111);
	g->Fit(p_plot, "S");//With the fit option S, you can access the full result of the fit including the covariance and correlation matrix.
	TVirtualFitter* fit = TVirtualFitter::GetFitter(); //to get fit results;
	// Get covariance matrix and print it
	//int npar = p_plot->GetNpar();
	//TMatrixD covMatrix(npar, npar, fit->GetCovarianceMatrix()); //commented as not useful for now, but may be later
	//covMatrix.Print();
	p_plot->SetParameter(0, p_plot->GetParameter(0));//setting the first parameter of the function, using the value given by the fit procedure
	//p_plot->SetParameter(1, p_plot->GetParameter(1));//setting the second parameter of the function, using the value given by the fit procedure


	//DRAW: I draw hist and p_plot in the same graph:
	g->Draw("");
	p_plot->Draw("SAME");


	//LATEX: I add LaTeX titles and axis labels:
	TLatex latex;
	latex.SetNDC(); //sets the use of Normalized Device Coordinates (NDC).
	latex.SetTextSize(0.05); //changes text size for title
	latex.DrawLatex(0.15, 0.92, "Autocorrelation function vs indexes distance:");
	latex.SetTextSize(0.04); //changes text size for axis labels
	latex.DrawLatex(0.5, 0.03, "k");
	latex.SetTextAngle(90);
	latex.DrawLatex(0.05, 0.47, "C_{x} (k)");

	//SAVE: I save the canvas as an image
	canvas->SaveAs((name_image).c_str());

	(*result) = p_plot->GetParameter(0);
	(*dresult) = p_plot->GetParError(0);
	
	//DELETE:
	delete p_plot;
	delete g;
	delete canvas;
}

//-----------------------------------------------------------------
//ROOT MACRO TO GRAPH:

void plot_points_orizline(vector<double>& x, vector<double>& y, double y0, string name_image) {

	//CANVAS creation: (to draw the graph)
	TCanvas* canvas = new TCanvas("canvas", "Canvas for Drawing Points", 1000, 600);
	canvas->SetGrid();//to set grid

	TGraph* g = new TGraph(y.size(), x.data(), y.data());
	//Using std::max_element: (useful to create my histogram object)
	//	-> std::max_element takes two iterators that define the range of the array to operate on. 
	//	-> returns an iterator that points to the maximum element found.
	//The usage of std::min_element is analogous but with for minimum element.
	auto min_x = *min_element(x.begin(), x.end());
	auto max_x = *max_element(x.begin(), x.end());
	auto min_y = *min_element(y.begin(), y.end());
	auto max_y = *max_element(y.begin(), y.end());

	max_y = max(max_y, y0);

	g->SetTitle(""); //to clear default title
	g->GetXaxis()->SetTitle(""); //to clear default X axis title
	g->GetYaxis()->SetTitle(""); //to clear default Y axis title

	g->SetLineColor(1);

	g->GetYaxis()->SetRangeUser(min_y - 0.01 * abs(min_y), max_y + 0.01 * abs(max_y));//to set y axis range
	g->GetXaxis()->SetRangeUser(min_x - 0.01 * abs(min_x), max_x);//to set x axis range

	g->Draw("ALP");

	TLine* line = new TLine(min_x - 0.01 * abs(min_x), y0, max_x + 0.01 * abs(min_x), y0);
	line->SetLineColor(kRed);
	line->SetLineStyle(2);    //dotted line
	line->SetLineWidth(2);    //line width
	line->Draw("same");

	//LATEX: I add LaTeX titles and axis labels:
	TLatex latex;
	latex.SetNDC(); //sets the use of Normalized Device Coordinates (NDC).
	latex.SetTextSize(0.05); //changes text size for title
	latex.DrawLatex(0.075, 0.92, "(the sum of C(k) up to kmax) vs (the maximum value of k considered):");
	latex.SetTextSize(0.04); //changes text size for axis labels
	latex.DrawLatex(0.5, 0.03, "k_{max}");
	latex.SetTextAngle(90);
	latex.DrawLatex(0.04, 0.47, "#sum_{k=1}^{k_{max}} C_{x}(k)");


	//SAVE: I save the canvas as an image
	canvas->SaveAs((name_image).c_str());

	//DELETE:
	delete g;
	delete canvas;
}
