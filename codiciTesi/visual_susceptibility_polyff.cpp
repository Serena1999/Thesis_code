/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****             visual_susceptibility_polyff.cpp:            ****
****     visualization of the susceptibilities of the the     ****
****           Polyakov loops and chiral condensates          ****
****             as a function of the temperature             ****
****    (we define a susceptibility Chi_p for an intensive    ****
****                      observable p as                     ****
****        Chi_p = N_s^3 (<p^2> - <p>^2) = V_s*Var(p))       ****
****         (ref.: https://arxiv.org/pdf/1911.07668)         ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

#include "../library.h"
#include "../root_include.h"

const string tipology = "gauge"; //gauge/fermion, CHOOSABLE --> to do the gauge/fermion observables graph
const int N_s = 32; //number of sites along each spatial direction
const int V_s = N_s * N_s * N_s;

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
//MAIN:

int main() {

	int skipLines = 1; //= number of lines to skip while reading input file;
	double temp_value, mod_value, chi_mod_value, re_value, chi_re_value, im_value, chi_im_value;
	vector <double> temp, chi_mod, chi_re, chi_im;
	size_t pos;
	string line, name_tmp;
	//string input_directory = "19_05_2025/polyff_results/";
	string input_directory = "results/";
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
	double width_canvas_mod;
	double width_canvas_re;
	double width_canvas_im;

	if (tipology == "gauge") {
		name_input_file = "1500.0_poly_results.txt";
		pos = name_input_file.find_last_of(".");
		if (pos != string::npos) {
			name_tmp = name_input_file.substr(0, pos); //I remove extension using substr
		}
		name_input_file = input_directory + name_input_file;
		name_image_mod = output_directory + "CHImodPvsT_" + name_tmp + ".png";
		name_image_re = output_directory + "CHIrePvsT_" + name_tmp + ".png";
		name_image_im = output_directory + "CHIimPvsT_" + name_tmp + ".png";
		title_mod = "#chi_{|PP^{+}|} vs temperature:";
		title_re = "#chi_{|Re(P)|} vs temperature:";
		title_im = "#chi_{|Im(P)|} vs temperature:";
		y_name_mod = "#chi_{|PP^{+}|}";
		y_name_re = "#chi_{|Re(P)|}";
		y_name_im = "#chi_{|Im(P)|}";
		pos_ymod = 0.03;
		pos_yre = 0.02;
		pos_yim = 0.04;
		height_mod = 0.45;
		height_re = 0.45;
		height_im = 0.45;
		pos_title_mod = 0.35;
		pos_title_re = 0.35;
		pos_title_im = 0.35;
		width_canvas_mod = 900;
		width_canvas_re = 1100;
		width_canvas_im = 900;

	}
	else if (tipology == "fermion") {
		name_input_file = "1500.0_ff_results.txt";
		pos = name_input_file.find_last_of(".");
		if (pos != string::npos) {
			name_tmp = name_input_file.substr(0, pos); //I remove extension using substr
		}
		name_input_file = input_directory + name_input_file;
		name_image_mod = output_directory + "CHImodffvsT_" + name_tmp + ".png";
		name_image_re = output_directory + "CHIreffvsT_" + name_tmp + ".png";
		name_image_im = output_directory + "CHIimffvsT_" + name_tmp + ".png";
		title_mod = "#chi_{|(#bar{#psi}#psi)(#bar{#psi}#psi)^{+}|} vs temperature:";
		title_re = "#chi_{Re(#bar{#psi}#psi)} vs temperature:";
		title_im = "#chi_{Im(#bar{#psi}#psi)} vs temperature:";
		y_name_mod = "#chi_{|(#bar{#psi}#psi)(#bar{#psi}#psi)^{+}|}";
		y_name_re = "#chi_{Re(#bar{#psi}#psi)}";
		y_name_im = "#chi_{Im(#bar{#psi}#psi)}";
		pos_ymod = 0.04;
		pos_yre = 0.03;
		pos_yim = 0.04;
		height_mod = 0.45;
		height_re = 0.45;
		height_im = 0.45;
		pos_title_mod = 0.35;
		pos_title_re = 0.35;
		pos_title_im = 0.35;
		width_canvas_mod = 900;
		width_canvas_re = 900;
		width_canvas_im = 900;
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
		if (iss >> temp_value >> mod_value >> chi_mod_value >> re_value >> chi_re_value >> im_value >> chi_im_value) {
			temp.push_back(temp_value);
			chi_mod_value = chi_mod_value * chi_mod_value * V_s;
			chi_re_value = chi_re_value * chi_re_value * V_s;
			chi_im_value = chi_im_value * chi_im_value * V_s;
			chi_mod.push_back(chi_mod_value);
			chi_re.push_back(chi_re_value);
			chi_im.push_back(chi_im_value);
		}
		else {
			cerr << "Poorly formatted line: " << line << endl;
		}
	}

	input_file.close();

	//GRAPHIC REPRESENTATION:
	plot_points(temp, chi_mod, name_image_mod, title_mod, y_name_mod, pos_title_mod, pos_ymod, height_mod, width_canvas_mod);
	plot_points(temp, chi_re, name_image_re, title_re, y_name_re, pos_title_re, pos_yre, height_re, width_canvas_re);
	plot_points(temp, chi_im, name_image_im, title_im, y_name_im, pos_title_im, pos_yim, height_im, width_canvas_im);

	return 0;
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
	g->GetYaxis()->SetRangeUser(min_y - 0.01 * fabs(min_y), max_y + 0.01 * fabs(max_y));
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
	latex.DrawLatex(0.45, 0.03, "T[MeV]");
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
