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
		243.377,
		257.347,
		271.711,
		286.507,
		301.786,
		317.618,
		334.088,
		351.294,
		369.346,
		388.360,
		404.645
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
		1.75081,
		1.87469,
		2.13651
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
		0.0147706,
		0.0160792,
		0.0775104//INDECISA SE METTERLO
	};

	plot_points_errors(
		x,
		y,
		dy,
		"rendering_mpi800_muT.png",
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
