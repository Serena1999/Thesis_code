/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                    visual_polyffT.cpp:                   ****
****    visualization of the mean values ​​of the observables   ****
****          (with corresponding standard deviation)         ****
****             as a function of the temperature             ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

#include "../library.h"
#include "../root_include.h"

const string tipology = "reff"; //YOU CAN CHOOSE BETWEEN reP, imP, reff, imff;
const bool bool_log_scale = 0;

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
//MAIN:

int main() {

	int skipLines = 1; //= number of lines to skip while reading input file;
	double temp_value, value, err_value;
	vector <double> temp, chi, chi_err;
	size_t pos;
	string line, name_tmp;
	string input_directory = "19_05_2025/data_susceptiblities_with_errors/";
	string output_directory = "results/";
	string name_input_file;
	pos = name_input_file.find_last_of(".");
	if (pos != std::string::npos) {
		name_tmp = name_input_file.substr(0, pos); //I remove extension using substr
	}
	string name_image;
	string title;
	string y_name;
	double pos_y;
	double height;
	double pos_title;
	double width_canvas;

	if (tipology == "reP") {
		name_input_file = "1500.0_reP_results.txt";
		pos = name_input_file.find_last_of(".");
		if (pos != string::npos) {
			name_tmp = name_input_file.substr(0, pos); //I remove extension using substr
		}
		name_input_file = input_directory + name_input_file;
		name_image = output_directory + "chi(reP)vsT_" + name_tmp + ".png";
		title = "#chi_{|Re(P)|} vs temperature:";
		y_name = "#chi_{|Re(P)|}";
		pos_y = 0.035;
		height = 0.45;
		pos_title = 0.35;
		width_canvas = 900;
	}
	else if (tipology == "imP") {
		name_input_file = "1500.0_imP_results.txt";
		pos = name_input_file.find_last_of(".");
		if (pos != string::npos) {
			name_tmp = name_input_file.substr(0, pos); //I remove extension using substr
		}
		name_input_file = input_directory + name_input_file;
		name_image = output_directory + "chi(imP)vsT_" + name_tmp + ".png";
		title = "#chi_{|Im(P)|} vs temperature:";
		y_name = "#chi_{|Im(P)|}";
		pos_y = 0.03;
		height = 0.45;
		pos_title = 0.35;
		width_canvas = 900;
	}
	else if (tipology == "reff") {
		name_input_file = "1500.0_reff_results.txt";
		pos = name_input_file.find_last_of(".");
		if (pos != string::npos) {
			name_tmp = name_input_file.substr(0, pos); //I remove extension using substr
		}
		name_input_file = input_directory + name_input_file;
		name_image = output_directory + "chi(reff)vsT_" + name_tmp + ".png";
		title = "#chi_{Re(#bar{#psi}#psi)} vs temperature:";
		y_name = "#chi_{Re(#bar{#psi}#psi)}";
		pos_y = 0.03;
		height = 0.45;
		pos_title = 0.35;
		width_canvas = 900;
	}
	else if (tipology == "imff") {
		name_input_file = "1500.0_imff_results.txt";
		pos = name_input_file.find_last_of(".");
		if (pos != string::npos) {
			name_tmp = name_input_file.substr(0, pos); //I remove extension using substr
		}
		name_input_file = input_directory + name_input_file;
		name_image = output_directory + "chi(imff)vsT_" + name_tmp + ".png";
		title = "#chi_{Im(#bar{#psi}#psi)} vs temperature:";
		y_name = "#chi_{Im(#bar{#psi}#psi)}";
		pos_y = 0.03;
		height = 0.45;
		pos_title = 0.35;
		width_canvas = 900;
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
		if (iss >> temp_value >> value >> err_value) {
			temp.push_back(temp_value);
			if (bool_log_scale) { 
				value = log(value);
				//err_value = log(err_value);
			}
			chi.push_back(value);
			chi_err.push_back(err_value);
		}
		else {
			cerr << "Poorly formatted line: " << line << endl;
		}
	}

	input_file.close();

	//GRAPHIC REPRESENTATION:
	plot_points_errors(temp, chi, chi_err, name_image, title, y_name, pos_title, pos_y, height, width_canvas);

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
	auto min_y = *min_element(y.begin(), y.end()) - *min_element(y_err.begin(), y_err.end());
	auto max_y = *max_element(y.begin(), y.end()) + *max_element(y_err.begin(), y_err.end());

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
