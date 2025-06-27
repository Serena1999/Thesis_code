//-----------------------------------------------------------------
//HEADERS && LIBRARY:

#include "../library.h"
#include "../root_include.h"

//TODO:
//INSERISCI L'OPZIONE DI SALTARE LE CONFIGURAZIONI DI TERMALIZZAZIONE; -> 2*
//INSERISCI MEDIA SU PIù CONFIGURAZIONI -> 1*
//MODIFICA NOMI IN LATEX -> 3*
//SALVA DATI SU FILE, COSì CHE POI FAI ANCHE LINEE A PIù T SOVRAPPOSTE -> 4*

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

int main() {

	string input_dir = "11_05_2025/build_good_015_out/";
	string name_input_file = "mon.dat";
	string name_image = "results/rhok_vs_rho1_mpi800_n15.png";

	string title = "title";
	string y_name = "<rho_k/rho_1>"; //MODIFICA POI 
	double pos_y = 0.020;
	double heigh_y = 0.45;
	double width_canvas = 900;

	string line;
	int Ns = 32;
	int Nt = 8;
	int Vs = Ns * Ns * Ns; //spatial volume (adimensional)
	double discard1, discard2, discard3, discard4, discard5, discard6, discard7, discard8, wrap_value;

	vector <double> n_k; //number of monopoles at given |k| wrapping
	vector <double>  k_array;
	vector <double> rho_k_norm; //= <rho_k/rho_1>; rho_k = n_k/(Vs*a^3)

	ifstream input_file;
	input_file.open(input_dir + name_input_file);
	if (!input_file) {
		cerr << "Error opening input file." << endl;
		return 1;
	}
	
	//reading monopoles vs |k|

	while (getline(input_file, line)) {
		istringstream iss(line);
		if (iss >> discard1 >> discard2 >> discard3 >> discard4 >> discard5 >> discard6 >> discard7 >> discard8 >> wrap_value) {
			if (wrap_value == 0) {
				continue;
			}
			while(abs(wrap_value) > n_k.size()) {
				n_k.push_back(0);
			}
			n_k[abs(wrap_value) - 1] += 1;
		}
		else {
			cerr << "Poorly formatted line: " << line << endl;
		}
	}

	input_file.close();

	if (n_k.size() > 1) {
		for (int ii = 1; ii < n_k.size(); ii++) {//vogliamo da rho_2/rho_1 in su con numeri di avvolgimenti
			if (n_k[0] == 0) {
				cout << "No monopoles with 1 wrapping: you need more statistics";
				return 1;
			}
			rho_k_norm.push_back((double)n_k[ii] / (double)n_k[0]);
			k_array.push_back(ii + 1);
		}
	}
	else {
		cout << "No monopoles exept with n_wrap = 1" << endl;
		return 1;
	}
	ofstream minimal_output;
	minimal_output.open("results/minimal_output.txt");
	if (!minimal_output) {
		cerr << "Error opening minimal output file." << endl;
		return 1;
	}
	
	if (rho_k_norm.size() == 1) {
		cout << endl;
		cout << "File: " << (input_dir + name_input_file) << endl;
		cout << "We have only n_wrap in 1 and 2" << endl;
		cout << "rho_2/rho_1 = " << rho_k_norm[0];
		minimal_output << "File: " << (input_dir + name_input_file) << endl;
		minimal_output << "We have only n_wrap in 1 and 2" << endl;
		minimal_output << "rho_2/rho_1 = " << rho_k_norm[0];
		return 1;
	}

	minimal_output.close();


	plot_points(
		k_array,
		rho_k_norm,
		name_image,
		title,
		y_name,
		5.0,
		pos_y,
		heigh_y,
		width_canvas
	);

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
	latex.DrawLatex(0.45, 0.03, "k");
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
