/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****              visual_monte_carlo_history.cpp:             ****
****           visualization of Monte Carlo history           ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

#include "../library.h"
#include "../root_include.h"


//-----------------------------------------------------------------
//ROOT MACRO TO DO FIT AND GRAPH:

void plot_points(vector<double>& x, vector<double>& y, string name_image, string title, string y_name, double heigh_y);

//-----------------------------------------------------------------
//MAIN:

int main() {

	int skipLines = 1; //= number of lines to skip while reading input file;
	vector<double> time, polyr_vec, polyi_vec; //to contain data
	int index = 0, conf_id;
	double value_tmp, poly_re, poly_im;
	string line;
	
	string directory = "11_05_2025/build_good_001_out/";
	string name_input_file = "gauge_obs3725763671.txt";
	name_input_file = directory + name_input_file;
	string name_image1 = "gauge_obs3725763671.png";
	string title = "Monte Carlo History";
	title = title + " gauge_obs3725763671:";
	string y_name1 = "Re{P}";
	string y_name2 = "Im{P}";
	string name_image2 = y_name2 + "_" + name_image1;
	name_image1 = y_name1 + "_" + name_image1;

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

	while (input_file >> conf_id >> value_tmp >> value_tmp >> value_tmp >> poly_re >> poly_im) {//DEVI CORREGGEERE ERRORI
		time.push_back((double)conf_id);
		polyr_vec.push_back(poly_re);
		polyi_vec.push_back(poly_im);
	}

	input_file.close();

	//GRAPHIC REPRESENTATION:
	plot_points(time, polyr_vec, name_image1, title, y_name1, 0.47);
	plot_points(time, polyi_vec, name_image2, title, y_name1, 0.47);

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
	g->GetXaxis()->SetRangeUser(min_x - 0.01 * abs(min_x), max_x);//to set x axis range

	gStyle->SetOptFit(1111);

	g->Draw("ALP");

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
