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
	double heigh_y,
	double width_canvas
);


int main(int argc, char** argv) {

	TApplication app("App", &argc, argv);

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

	plot_points_errors(
		x,
		y,
		dy,
		"rendering_mpi1500_muT_1.png",
		"rendering",
		"#mu",
		0.4,
		0.02,
		0.5,
		900
	);

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
	g_errors->GetXaxis()->SetLimits(min_x-10, max_x+10);
	g_errors->GetYaxis()->SetRangeUser(0, max_y + 0.01 * fabs(max_y));
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

	gApplication->Run(true); // <--- TRUE = non bloccare il terminale

	//DELETE:
	delete g_line;
	delete g_errors;
	delete canvas;
}
