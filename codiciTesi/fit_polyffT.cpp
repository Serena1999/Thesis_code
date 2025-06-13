/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                     fit_polyffT.cpp:                     ****
****      visualization and fit of the mean values ​​of the     **** 
****    observables (with corresponding standard deviation)   ****
****             as a function of the temperature             ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

#include "../library.h"
#include "../root_include.h"

const string tipology = "gauge"; //gauge/fermion, CHOOSABLE --> to do the gauge/fermion observables graph
#define CHOOSE_FIT_FUNCTION 0 //0 for polynomial, 1 for arctg, 2 for logistic function;

//FIT FUNCTION: 

#if CHOOSE_FIT_FUNCTION == 0
	double fit_function(double* y, double* par) {//
		return par[0] + par[1] * y[0] + par[2] * y[0] * y[0] + par[3] * y[0] * y[0] * y[0];
	};

	void par_estimate(//considering polynomial as y(x) = p[0] + p[1]*x + p[2] * x^2 + p[3] * x^3;
		vector<double>& x,
		vector<double>& y,
		vector<double>& p
	) {

		int step = x.size() / 4;

		//using Cramer method to solve system of equations:
		double detM = (x[0] - x[step]) * (x[0] - x[2 * step]) * (x[step] - x[2 * step]) * (x[0] - x[3 * step]) * (x[step] - x[3 * step]) * (x[2 * step] - x[3 * step]);
		
		if (detM == 0) {
			cout << "detM = 0, so I use {1,1,1,1} as init" << endl;
			p[0] = 1;
			p[1] = 1;
			p[2] = 1;
			p[3] = 1;
			return;
		}
		
		p[0] += -x[0] * x[step] * (x[0] - x[step]) * (x[2 * step] * y[step * 3] * (x[0] - x[2 * step]) * (x[2 * step] - x[step]) + x[3 * step] * y[step * 2] * (x[0] - x[3 * step]) * (x[step] - x[3 * step])) + x[0] * x[2 * step] * x[3 * step] * y[step] * (x[0] - x[2 * step]) * (x[0] - x[3 * step]) * (x[2 * step] - x[3 * step]) - x[step] * x[2 * step] * x[3 * step] * y[0] * (x[step] - x[2 * step]) * (x[step] - x[3 * step]) * (x[2 * step] - x[3 * step]);
		p[1] += y[step] * (pow(x[0], 3) * (pow(x[step * 3], 2) - pow(x[step * 2], 2)) + pow(x[0], 2) * (pow(x[step * 2], 3) - pow(x[step * 3], 3)) + pow(x[step * 2], 2) * pow(x[step * 3], 2) * (x[step * 3] - x[step * 2])) + (x[0] - x[step]) * (y[step * 3] * (x[0] - x[step * 2]) * (x[step * 2] - x[step]) * (x[0] * (x[step] + x[step * 2]) + x[step] * x[step * 2]) + y[step * 2] * (x[0] - x[step * 3]) * (x[step] - x[step * 3]) * (x[0] * (x[step] + x[step * 3]) + x[step] * x[step * 3])) + y[0] * (pow(x[step], 3) * (pow(x[step * 2], 2) - pow(x[step * 3], 2)) + pow(x[step], 2) * (pow(x[step * 3], 3) - pow(x[step * 2], 3)) + pow(x[step * 2], 2) * pow(x[step * 3], 2) * (x[step * 2] - x[step * 3]));
		p[2] += pow(x[0], 3) * (-x[step]) * y[step * 2] + pow(x[0], 3) * x[step] * y[step * 3] + y[step] * (pow(x[0], 3) * (x[step * 2] - x[step * 3]) + x[0] * (pow(x[step * 3], 3) - pow(x[step * 2], 3)) + x[step * 2] * x[step * 3] * (pow(x[step * 2], 2) - pow(x[step * 3], 2))) - pow(x[0], 3) * x[step * 2] * y[step * 3] + pow(x[0], 3) * x[step * 3] * y[step * 2] + x[0] * pow(x[step], 3) * y[step * 2] - x[0] * pow(x[step], 3) * y[step * 3] + x[0] * pow(x[step * 2], 3) * y[step * 3] - x[0] * pow(x[step * 3], 3) * y[step * 2] + y[0] * (pow(x[step], 3) * (x[step * 3] - x[step * 2]) + x[step] * (pow(x[step * 2], 3) - pow(x[step * 3], 3)) - pow(x[step * 2], 3) * x[step * 3] + x[step * 2] * pow(x[step * 3], 3)) + pow(x[step], 3) * x[step * 2] * y[step * 3] - pow(x[step], 3) * x[step * 3] * y[step * 2] - x[step] * pow(x[step * 2], 3) * y[step * 3] + x[step] * pow(x[step * 3], 3) * y[step * 2];
		p[3] += (x[0] - x[step]) * (y[step * 3] * (x[0] - x[step * 2]) * (x[step * 2] - x[step]) + y[step * 2] * (x[0] - x[step * 3]) * (x[step] - x[step * 3])) - y[step] * (x[0] - x[step * 2]) * (x[0] - x[step * 3]) * (x[step * 2] - x[step * 3]) + y[0] * (x[step] - x[step * 2]) * (x[step] - x[step * 3]) * (x[step * 2] - x[step * 3]);
		
		

		p[0] /= detM;
		p[1] /= detM;
		p[2] /= detM;
		p[3] /= detM;

		cout << "Parameters used for fit: [" << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << "]" << endl;
	}


#elif CHOOSE_FIT_FUNCTION == 1
	double fit_function(double* y, double* par) {//
		return par[0] + par[1] * atan(par[2] * (y[0] - par[3]));
	};


	/*
	
	void par_estimate(//considering polynomial as par[0] + par[1] * atan(par[2] * (x[0] - par[3]));
		vector<double>& x,
		vector<double>& y,
		vector<double>& p
	) {
		/
		p[0] = 0.004;
		p[1] = 0.06;
		p[2] = 0.003;
		p[3] = 160;
		/
		p[0] = 0.02;
		p[1] = 0.05;
		p[2] = 0.015;
		p[3] = 200;

		cout << "Parameters used for fit: [" << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << "]" << endl;
	}
	*/
	
	void par_estimate(
		vector<double>& x, // temperature
		vector<double>& y, // observable
		vector<double>& p  // parameters to fill
	) {

		if (0) {
			p[0] = 0.024;// 0.005;
			p[1] = 0.1 /PI;
			p[2] = 0.006;//0.002;//0.005
			p[3] = 250; // 180;
			return;
		}

		size_t n = x.size();
		if (n < 5) {
			cerr << "Not enough points for estimate." << endl;
			p = { 0, 1, 0.01, 250 }; // fallback
			return;
		}

		/*
		// ---- 1. p3: punto della transizione (lo scelgo ad occhio)
		int idx = 6;
		p[3] = x[idx]; // T al massimo gradiente

		// ---- 2. p0: punto del flesso:
		p[0] = y[idx];

		// ---- 3. p1: altezza del gradino/PI
		p[1] = 2.0 * p[0] / PI;

		// ---- 4. p2: uso le stime precedenti:
		p[2] = tan((y[idx+1] - p[0]) / p[1]) / (x[idx+1] - p[3]);

		cout << "Initial guess: p0 = " << p[0]
			<< ", p1 = " << p[1]
			<< ", p2 = " << p[2]
			<< ", p3 = " << p[3] << endl;

		*/

		/*
			// 1. Trova il punto di massimo gradiente (flesso)
		int idx = 1;
		double max_der = -1e10;
		for (int i = 1; i < n - 1; i++) {
			double dy = y[i + 1] - y[i - 1];
			double dx = x[i + 1] - x[i - 1];
			double der = dy / dx;
			if (abs(der) > max_der) {
				max_der = abs(der);
				idx = i;
			}
		}
		p[3] = x[idx]; // centro del gradino

		// 2. p0: valore nel flesso (altezza al centro)
		p[0] = y[idx];

		// 3. p1: ampiezza del gradino / PI
		double y_max = *max_element(y.begin(), y.end());
		double y_min = *min_element(y.begin(), y.end());
		p[1] = (y_max - y_min) / PI;

		// 4. p2: ripidità (1/larghezza tra 20% e 80%)
		double y20 = y_min + 0.2 * (y_max - y_min);
		double y80 = y_min + 0.8 * (y_max - y_min);
		double T_low = x.front(), T_high = x.back();
		for (int i = 1; i < n; i++) {
			if ((y[i - 1] < y20 && y[i] > y20)) T_low = x[i];
			if ((y[i - 1] < y80 && y[i] > y80)) {
				T_high = x[i];
				break;
			}
		}
		double width = T_high - T_low;
		p[2] = (width > 0.0) ? 1.0 / width : 1.0;
		*/

		// ---- 1. p3: punto x della transizione (lo scelgo ad occhio)
		p[3] = 240; // 210; // T al massimo gradiente

		// ---- 2. p0: punto y della transizione (lo scelgo ad occhio):
		p[0] = 0.02; //0.015;

		// ---- 3. p1: altezza del gradino/PI
		p[1] = 2.0 / PI * (y[x.size() - 1] * 1.1 - p[0]);

		// ---- 4. p2: uso le stime precedenti:
		p[2] = tan((y[2] - p[0]) / p[1]) / (x[2] - p[3]);

		cout << "Initial guess: p0 = " << p[0]
			<< ", p1 = " << p[1]
			<< ", p2 = " << p[2]
			<< ", p3 = " << p[3] << endl;

	}
	
#elif CHOOSE_FIT_FUNCTION == 2

	double fit_function(double* y, double* par) {
		return par[0] + par[1] / (1 + exp(-par[2] * (y[0] - par[3])));
	}

	void par_estimate(
		vector<double>& x, // temperature
		vector<double>& y, // observable
		vector<double>& p  // parameters to fill
	) {

		if (1) {
			p[0] = 0;
			p[1] = 1.0;
			p[2] = 1.0;
			p[3] = 200;
			return;
		}

		int n = y.size() - 1;
		p[0] = 0.0; // baseline
		p[1] = y[n - 1] - p[0]; // ampiezza
		// Trova l'indice del massimo gradiente
		int idx = 1;
		double max_der = -1e10;
		for (int i = 1; i < n - 1; i++) {
			double dy = y[i + 1] - y[i - 1];
			double dx = x[i + 1] - x[i - 1];
			double der = dy / dx;
			if (der > max_der) {
				max_der = der;
				idx = i;
			}
		}
		p[3] = x[idx]; // posizione centrale
		// Calcola larghezza tra 20% e 80%
		double y20 = p[0] + 0.2 * p[1];
		double y80 = p[0] + 0.8 * p[1];
		double x20 = x[0], x80 = x[n - 1];
		for (int i = 1; i < n; i++) {
			if ((y[i - 1] < y20 && y[i] > y20)) x20 = x[i];
			if ((y[i - 1] < y80 && y[i] > y80)) { x80 = x[i]; break; }
		}
		double width = x80 - x20;
		p[2] = width > 0 ? 4.0 / width : 1.0;

		cout << "Initial guess: p0 = " << p[0]
			<< ", p1 = " << p[1]
			<< ", p2 = " << p[2]
			<< ", p3 = " << p[3] << endl;

	}

#endif


//-----------------------------------------------------------------
//ROOT MACRO TO DO FIT AND GRAPH:

void fit_plot_points_errors(
	vector<double>& x,
	vector<double>& y,
	vector<double>& y_err,
	string name_image,
	string title,
	string y_name,
	double pos_title,
	double pos_y,
	double heigh_y
);

void silly_plot(
	vector<double>& x,
	vector<double>& y,
	vector<double>& y_err,
	string name_image,
	string title,
	string y_name,
	double pos_title,
	double pos_y,
	double heigh_y
);


//-----------------------------------------------------------------
//MAIN:

int main() {

	int skipLines = 1; //= number of lines to skip while reading input file;
	double temp_value, mod_value, mod_err_value, re_value, re_err_value, im_value, im_err_value;
	vector <double> temp, mod, mod_err, re, re_err, im, im_err;
	size_t pos;
	string line, name_tmp;
	string input_directory = "11_05_2025/polyff_results/";
	string output_directory = "results/";
	string name_input_file;
	string name_image_mod;
	string name_image_re;
	string name_image_im;
	string title_mod;
	string title_re;
	string title_im;
	string y_name_mod;
	string y_name_re;
	string y_name_im;
	double pos_ymod;
	double pos_yre;
	double pos_yim;
	double height_mod;
	double height_re;
	double height_im;
	double pos_title_mod;
	double pos_title_re;
	double pos_title_im;

	if (tipology == "gauge") {
		name_input_file = "800.0_poly_results.txt";
		pos = name_input_file.find_last_of(".");
		if (pos != string::npos) {
			name_tmp = name_input_file.substr(0, pos); //I remove extension using substr
		}
		name_input_file = input_directory + name_input_file;
		name_image_mod = output_directory + "FIT_modPvsT_" + name_tmp + ".png";
		name_image_re = output_directory + "FIT_rePvsT_" + name_tmp + ".png";
		name_image_im = output_directory + "FIT_imPvsT_" + name_tmp + ".png";
		title_mod = "#LT|PP^{+}|#GT vs temperature:";
		title_re = "#LTRe{P}#GT vs temperature :";
		title_im = "#LTIm{P}#GT vs temperature:";
		y_name_mod = "#LT|PP^{+}|#GT";
		y_name_re = "#LTRe{P}#GT";
		y_name_im = "#LTIm{P}#GT";
		pos_ymod = 0.02;
		pos_yre = 0.03;
		pos_yim = 0.03;
		height_mod = 0.4;
		height_re = 0.45;
		height_im = 0.45;
		pos_title_mod = 0.3;
		pos_title_re = 0.3;
		pos_title_im = 0.3;
	}
	else if (tipology == "fermion") {
		name_input_file = "800.0_ff_results.txt";
		pos = name_input_file.find_last_of(".");
		if (pos != string::npos) {
			name_tmp = name_input_file.substr(0, pos); //I remove extension using substr
		}
		name_input_file = input_directory + name_input_file;
		name_image_mod = output_directory + "FIT_modffvsT_" + name_tmp + ".png";
		name_image_re = output_directory + "FIT_reffvsT_" + name_tmp + ".png";
		name_image_im = output_directory + "FIT_imffvsT_" + name_tmp + ".png";
		title_mod = "#LT|(#bar{#psi}#psi)(#bar{#psi}#psi)^{+}|#GT vs temperature:";
		title_re = "#LTRe{#bar{#psi}#psi}#GT vs temperature :";
		title_im = "#LTIm{#bar{#psi}#psi}#GT vs temperature:";
		y_name_mod = "#LT|(#bar{#psi}#psi)(#bar{#psi}#psi)^{+}|#GT";
		y_name_re = "#LTRe{#bar{#psi}#psi}#GT";
		y_name_im = "#LTIm{#bar{#psi}#psi}#GT";
		pos_ymod = 0.03;
		pos_yre = 0.03;
		pos_yim = 0.04;
		height_mod = 0.4;
		height_re = 0.45;
		height_im = 0.45;
		pos_title_mod = 0.3;
		pos_title_re = 0.3;
		pos_title_im = 0.3;
	}
	else {
		cerr << "Invalid tipology: you must write gauge or fermion";
		return 1;
	}

	ifstream input_file;
	input_file.open(name_input_file);
	if (!input_file) {
		cerr << "Error opening input file." << endl;
		return 1;
	}

	for (int i = 0; i < skipLines; i++) {
		if (!getline(input_file, line)) {
			cerr << "Error: there are less than " << skipLines << " lines in the input file." << endl;
			return 1;
		}
	}

	while (getline(input_file, line)) {
		istringstream iss(line);
		if (iss >> temp_value >> mod_value >> mod_err_value >> re_value >> re_err_value >> im_value >> im_err_value) {
			temp.push_back(temp_value);
			mod.push_back(mod_value);
			mod_err.push_back(mod_err_value);
			re.push_back(re_value);
			re_err.push_back(re_err_value);
			im.push_back(im_value);
			im_err.push_back(im_err_value);
		}
		else {
			cerr << "Poorly formatted line: " << line << endl;
		}
	}

	input_file.close();

	//GRAPHIC REPRESENTATION:
	silly_plot(temp, mod, mod_err, "results/silly_plot_1.png", title_mod, y_name_mod, pos_title_mod, pos_ymod, height_mod);
	silly_plot(temp, re, re_err, "results/silly_plot_2.png", title_re, y_name_re, pos_title_re, pos_yre, height_re);
	silly_plot(temp, im, im_err, "results/silly_plot_3.png", title_im, y_name_im, pos_title_im, pos_yim, height_im);


	//GRAPHIC REPRESENTATION:
	fit_plot_points_errors(temp, mod, mod_err, name_image_mod, title_mod, y_name_mod, pos_title_mod, pos_ymod, height_mod);
	fit_plot_points_errors(temp, re, re_err, name_image_re, title_re, y_name_re, pos_title_re, pos_yre, height_re);
	fit_plot_points_errors(temp, im, im_err, name_image_im, title_im, y_name_im, pos_title_im, pos_yim, height_im);

	return 0;
}

//-----------------------------------------------------------------
//ROOT MACRO TO GRAPH:

void fit_plot_points_errors(
	vector<double>& x,
	vector<double>& y,
	vector<double>& y_err,
	string name_image,
	string title,
	string y_name,
	double pos_title,
	double pos_y,
	double heigh_y
) {

	if (x.size() != y.size() || y.size() != y_err.size()) {
		cerr << "Error: mismatched vector sizes in plot_points_errors()." << endl;
		return;
	}

	vector <double> par(4,0);

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

	double fit_min = min_x - 20;
	double fit_max = max_x + 20;

	TF1* p_plot = new TF1("fit_function", fit_function, fit_min, fit_max, 4); //4 is the number of parameters


	par_estimate(x, y, par);

	for (int ii = 0; ii < par.size(); ii++) {
		p_plot->SetParameter(ii, par[ii]); //setting the ii-th parameter of the function
	}

	p_plot->GetXaxis()->SetRangeUser(min_x, max_x);
	p_plot->GetYaxis()->SetRangeUser(0, g_errors->GetMaximum() * 1.01);
	
	
	//FIT:
	gStyle->SetOptFit(0);// (1111);
	g_errors->Fit(p_plot, "SW"); // S = return result, W = use weights (i.e., y errors)
	TVirtualFitter* fit = TVirtualFitter::GetFitter(); //to get fit results;

	p_plot ->SetLineColor(kRed);
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
	latex.DrawLatex(pos_y, heigh_y, y_name.c_str());


	//SAVE: I save the canvas as an image
	canvas->SaveAs(name_image.c_str());

	//to save also in vectorial pdf form:
	canvas->SaveAs((name_image.substr(0, name_image.find_last_of(".")) + ".pdf").c_str());

	//DELETE:
	delete p_plot;
	delete g_errors;
	delete canvas;
}

void silly_plot(
	vector<double>& x,
	vector<double>& y,
	vector<double>& y_err,
	string name_image,
	string title,
	string y_name,
	double pos_title,
	double pos_y,
	double heigh_y
) {
	if (x.size() != y.size() || y.size() != y_err.size()) {
		cerr << "Error: mismatched vector sizes in plot_points_errors()." << endl;
		return;
	}

	vector <double> par(4, 0);
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

	TF1* f1 = new TF1("f", fit_function, min_x, max_x, 4);
	
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