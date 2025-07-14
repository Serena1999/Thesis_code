/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****          0Tfile_visual_monte_carlo_history.cpp:          **** 
****  visualization of Monte Carlo history of fermionic files **** 
***                    at zero temperature.                   ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

//-----------------------------------------------------------------
//HEADERS && LIBRARY:

#include "../library.h"
#include "../root_include.h"

//-----------------------------------------------------------------
//ROOT MACRO TO DO FIT AND GRAPH:

void plot_points(vector<double>& x, vector<double>& y, string name_image, string title, string y_name, double heigh_y);

//-----------------------------------------------------------------
//MAIN:

int main() {

	int skipLines_file_list = 1, skipLines = 1; //= number of lines to skip while reading input file;
	
	vector<double> time, reff_vec, imff_vec; //to contain data
	int conf_id;
	double value_tmp, ff_re, ff_im;
	string line, y_name1, y_name2, reff_image, imff_image, title;
	vector<string> paths;
	string name_file_list = "19_05_2025/file_list_0T.txt";

	ifstream file_list;
	file_list.open(name_file_list);
	if (!file_list) {
		cout << "Error opening file list" << endl;
		return 1;
	}

	for (int i = 0; i < skipLines; i++) {
		if (!getline(file_list, line)) {
			cerr << "Error: there are less than " << skipLines_file_list << " lines in the file list." << endl;
			return 1;
		}
	}

	while (getline(file_list, line)) {
		istringstream iss(line);
		string single_path;
		if (iss >> single_path) {
			paths.push_back(single_path);
		}
		else {
			cerr << "Poorly formatted line: " << line << endl;
		}
	}

	file_list.close();

	size_t pos, pos_slash;
	string name_tmp;
	
	y_name1 = "Re{ff}";
	y_name2 = "Im{ff}";

	for (int ii = 0; ii < paths.size(); ii++) {
		ifstream input_file; //declaration of input file
		input_file.open(paths[ii]);
		if (!input_file) {
			cerr << "Error opening " << paths[ii] << endl;
			return 1;
		}

		pos = paths[ii].find_last_of(".");
		pos_slash = paths[ii].find_last_of("/");
		if (pos != string::npos) {
			paths[ii] = paths[ii].substr(0, pos); //I remove extension using substr
		}
		if (pos_slash != string::npos) {
			paths[ii] = paths[ii].substr(pos_slash + 1); //I remove extension using substr
		}

		cout << paths[ii] << endl;

		title = "Monte Carlo History " + paths[ii] + ":";
		reff_image = y_name1 + "_" + paths[ii] + ".png";
		imff_image = y_name2 + "_" + paths[ii] + ".png";

		for (int i = 0; i < skipLines; i++) {
			if (!getline(input_file, line)) {
				cerr << "Error: there are less than " << skipLines << " lines in the file" << paths[ii] << endl;
				return 1;
			}
		}

		while (getline(input_file, line)) {
			istringstream iss(line);
			if (iss >> conf_id >> value_tmp >> value_tmp >> value_tmp >> ff_re >> ff_im) { //read the line only until ff_im: if next elements are omitted isn't interesting.
				if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
				time.push_back(static_cast<double>(conf_id));
				reff_vec.push_back(ff_re);
				imff_vec.push_back(ff_im);
			}
			else {
				cerr << "Skipped line (badly formatted): " << line << endl;
			}
		}

		input_file.close();

		//GRAPHIC REPRESENTATION:
		plot_points(time, reff_vec, reff_image, title, y_name1, 0.47);
		cout << ii << ": Re{ff} image DONE!" << endl;
		plot_points(time, imff_vec, imff_image, title, y_name2, 0.47);
		cout << ii << ": Im{ff} image DONE!" << endl;
		time.clear();
		reff_vec.clear();
		imff_vec.clear();
	}

	return 0;
}

//-----------------------------------------------------------------
//ROOT MACRO TO GRAPH:

void plot_points(vector<double>& x, vector<double>& y, string name_image, string title, string y_name, double heigh_y) {

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
	g->GetXaxis()->SetLimits(min_x - 100, max_x + 100);//to set x axis range //SetLimits cannot be ignored by ROOT, while SetRangeUser yes.

	gStyle->SetOptFit(1111);

	g->Draw("ALP");
	gPad->Update();

	//LATEX: I add LaTeX titles and axis labels:
	TLatex latex;
	latex.SetNDC(); //sets the use of Normalized Device Coordinates (NDC).
	latex.SetTextSize(0.05); //changes text size for title
	latex.DrawLatex(0.15, 0.92, title.c_str());
	latex.SetTextSize(0.04); //changes text size for axis labels
	latex.DrawLatex(0.4, 0.03, "Monte Carlo time");
	latex.SetTextAngle(90);
	latex.DrawLatex(0.03, heigh_y, y_name.c_str());


	//SAVE: I save the canvas as an image
	canvas->SaveAs(("results/" + name_image).c_str());

	//DELETE:
	delete g;
	delete canvas;
}
