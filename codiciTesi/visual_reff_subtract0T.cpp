/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

//-----------------------------------------------------------------
//HEADERS && LIBRARY:

#include "../library.h"
#include "../root_include.h"

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

	string input_directory = "19_05_2025/polyff_results/";
	string input_directory_0T = "19_05_2025/0T_ff_results/";
	string output_directory = "results/";
	string name_input_file = "1500.0_ff_results.txt";
	string name_input_file_0T = "0T_1500.0_ff_results.txt";
	string name_output_file = "0Tsubtracted_1500.0_ff_results.txt";
	int skipLines = 1; //= number of lines to skip while reading input file;
	int skipLines_0T = 1;
	
	double temp_value, mod_value, mod_err_value, re_value, re_err_value, im_value, im_err_value;
	vector <double> temp, mod, mod_err, re, re_err, im, im_err;
	vector <double> mod_0T, mod_err_0T, re_0T, re_err_0T, im_0T, im_err_0T;
	size_t pos;
	string line, name_tmp;
	pos = name_input_file.find_last_of(".");
	if (pos != string::npos) {
		name_tmp = name_input_file.substr(0, pos); //I remove extension using substr
	}
	name_input_file = input_directory + name_input_file;
	name_input_file_0T = input_directory_0T + name_input_file_0T;

	string name_image_mod = output_directory + "modffvsT_" + name_tmp + "_subracted0T.png";
	string name_image_re = output_directory + "reffvsT_" + name_tmp + "_subracted0T.png";
	string name_image_im = output_directory + "imffvsT_" + name_tmp + "_subracted0T.png";
	string title_mod = "#LT|(#bar{#psi}#psi)(#bar{#psi}#psi)^{+}|#GT vs temperature:";
	string title_re = "#LTRe{#bar{#psi}#psi}#GT vs temperature:";
	string title_im = "#LTIm{#bar{#psi}#psi}#GT vs temperature:";
	string y_name_mod = "#LT|(#bar{#psi}#psi)(#bar{#psi}#psi)^{+}|#GT";
	string y_name_re = "#LTRe{#bar{#psi}#psi}#GT-#LTRe{#bar{#psi}#psi}#GT_{T=0}";
	string y_name_im = "#LTIm{#bar{#psi}#psi}#GT";
	double pos_ymod = 0.03;
	double pos_yre = 0.03;
	double pos_yim = 0.04;
	double height_mod = 0.4;
	double height_re = 0.35;
	double height_im = 0.45;
	double pos_title_mod = 0.3;
	double pos_title_re = 0.3;
	double pos_title_im = 0.3;
	

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


	ifstream input_file_0T;
	input_file_0T.open(name_input_file_0T);
	if (!input_file_0T) {
		cerr << "Error opening 0T-input file." << endl;
		return 1;
	}

	for (int i = 0; i < skipLines_0T; i++) {
		if (!getline(input_file_0T, line)) {
			cerr << "Error: there are less than " << skipLines_0T << " lines in the 0T-input file." << endl;
			return 1;
		}
	}

	while (getline(input_file_0T, line)) {
		istringstream iss(line);
		if (iss >> temp_value >> mod_value >> mod_err_value >> re_value >> re_err_value >> im_value >> im_err_value) {
			mod_0T.push_back(mod_value);
			mod_err_0T.push_back(mod_err_value);
			re_0T.push_back(re_value);
			re_err_0T.push_back(re_err_value);
			im_0T.push_back(im_value);
			im_err_0T.push_back(im_err_value);
		}
		else {
			cerr << "Poorly formatted line: " << line << endl;
		}
	}

	input_file_0T.close();

	if (temp.size() != mod_0T.size()) {
		cerr << "Problem in size" << endl;
		return 1;
	}


	ofstream output_file;
	output_file.open(output_directory + name_output_file);
	if (!output_file) {
		cerr << "Error opening output file." << endl;
		return 1;
	}

	output_file << "# T 	 |<ff * ff^dag>| 	 err(|<ff * ff^dag>|) 	 Re{ff} 	 err(Re{ff}) 	 Im{ff} 	 err(Im{ff})" << endl;

	for (int ii = 0; ii < temp.size(); ++ii) {
		mod[ii] = mod[ii] - mod_0T[ii];
		mod_err[ii] = sqrt(mod_err[ii] * mod_err[ii] + mod_err_0T[ii] * mod_err_0T[ii]);
		re[ii] = re[ii] - re_0T[ii];
		re_err[ii] = sqrt(re_err[ii] * re_err[ii] + re_err_0T[ii] * re_err_0T[ii]);
		im[ii] = im[ii] - im_0T[ii];
		im_err[ii] = sqrt(im_err[ii] * im_err[ii] + im_err_0T[ii] * im_err_0T[ii]);;


		//mod_err[ii] = sqrt(pow(mod_err[ii] / mod_0T[ii], 2) + pow(mod[ii] * mod_err_0T[ii] / pow(mod_0T[ii],2), 2));
		//mod[ii] = mod[ii] / mod_0T[ii];
		//re_err[ii] = sqrt(pow(re_err[ii] / re_0T[ii], 2) + pow(re[ii] * re_err_0T[ii] / pow(re_0T[ii], 2), 2));
		//re[ii] = re[ii] / re_0T[ii];
		//im_err[ii] = sqrt(pow(im_err[ii] / im_0T[ii], 2) + pow(im[ii] * im_err_0T[ii] / pow(im_0T[ii], 2), 2));
		//im[ii] = im[ii] / im_0T[ii];

		//mod[ii] -= 1;
		//re[ii] -= 1;
		//im[ii] -= 1;

		output_file << temp[ii] << "\t" << mod[ii] << "\t" << mod_err[ii] << "\t" << re[ii] << "\t" << re_err[ii] << "\t" << im[ii] << "\t" << im_err[ii] << endl;

		cout << re[ii] << "\t" << re_err[ii] << endl;
	}

	output_file.close();

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
