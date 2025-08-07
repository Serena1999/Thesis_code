/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                    visual_polyffT.cpp:                   ****
****    visualization of the mean values ​​of the observables   **** 
****          (with corresponding standard deviation)         ****
****             as a function of the temperature             ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

//-----------------------------------------------------------------
//HEADERS && LIBRARY:

#include "../library.h"
#include "../root_include.h"

//-----------------------------------------------------------------
//VARIABLES TO SET:

const string tipology = "fermionT0"; //gauge/fermion/fermionT0, CHOOSABLE --> to do the gauge/fermion observables graph

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
	double heigh_y
);

//-----------------------------------------------------------------
//MAIN:

int main() {

	int skipLines = 1; //= number of lines to skip while reading input file;
	//vector<double> time1, time2, polyr_vec, polyi_vec, reff_vec, imff_vec; //to contain data
	//int index = 0, conf_id;
	//double value_tmp, poly_re, poly_im, ff_re, ff_im;
	double temp_value, mod_value, mod_err_value, re_value, re_err_value, im_value, im_err_value;
	vector <double> temp, mod, mod_err, re, re_err, im, im_err;
	size_t pos;
	string line, name_tmp;
	string input_directory = "19_05_2025/polyff_results/";
	string output_directory = "results/";
	string name_input_file;
	pos = name_input_file.find_last_of(".");
	if (pos != std::string::npos) {
		name_tmp = name_input_file.substr(0, pos); //I remove extension using substr
	}
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
		name_image_mod = output_directory + "modPvsT_" + name_tmp + ".png";
		name_image_re = output_directory + "rePvsT_" + name_tmp + ".png";
		name_image_im = output_directory + "imPvsT_" + name_tmp + ".png";
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
		name_image_mod = output_directory + "modffvsT_" + name_tmp + ".png";
		name_image_re = output_directory + "reffvsT_" + name_tmp + ".png";
		name_image_im = output_directory + "imffvsT_" + name_tmp + ".png";
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
	else if (tipology == "fermionT0") {
		name_input_file = "0Tsubtracted_1500.0_ff_results.txt";
		pos = name_input_file.find_last_of(".");
		if (pos != string::npos) {
			name_tmp = name_input_file.substr(0, pos); //I remove extension using substr
		}
		name_input_file = input_directory + name_input_file;
		name_image_mod = output_directory + "modffvsT_" + name_tmp + ".png";
		name_image_re = output_directory + "reffvsT_" + name_tmp + ".png";
		name_image_im = output_directory + "imffvsT_" + name_tmp + ".png";
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
	plot_points_errors(temp, mod, mod_err, name_image_mod, title_mod, y_name_mod, pos_title_mod, pos_ymod, height_mod);
	plot_points_errors(temp, re, re_err, name_image_re, title_re, y_name_re, pos_title_re, pos_yre, height_re);
	plot_points_errors(temp, im, im_err, name_image_im, title_im, y_name_im, pos_title_im, pos_yim, height_im);

	return 0;
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
	double heigh_y
) {

	if (x.size() != y.size() || y.size() != y_err.size()) {
		cerr << "Error: mismatched vector sizes in plot_points_errors()." << endl;
		return;
	}
	
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
	delete g_line;
	delete g_errors;
	delete canvas;
}
