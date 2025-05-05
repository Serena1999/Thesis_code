/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                    autocorr_check.cpp:                   ****
****               using the blocking technique               ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

#include "../library.h"
#include "../root_include.h"

//1h/3h done: then errors must be implemented more carefully 
//sistema nomi sugli assi, passali come parametri;
//poi controlla se abbastanza inversioni per condensato chirale e fai analisi anche di quello.
//poi scegli bene intervallo di 15? valori e + then other 4 of them in CINECA_SCRATCH (CHECK);

void plot_points(vector<double>& x, vector<double>& y, string name_image);

//-----------------------------------------------------------------
//MAIN:

int main() {
	int Nt = 8; //BE CAREFUL TO CHOOSE IT WELL;
	int skipLines = 1; //= number of lines to skip while reading input file;
	int index = 0, index_new = 0;
	int n_sub_max;
	int N_tmp;
	double n_sub_ratio = 0.03;//0.03;//=0.5*len(data) -> you can modify this number from 0 to 1;
	double mpi = 800; //MeV //BE CAREFUL TO CHOOSE IT WELL;
	double poly, poly_re, poly_im, delta, delta_re, delta_im, value_tmp;
	double mean = 0, var_m = 0, mean_re = 0, var_re = 0, mean_im = 0, var_im = 0;
	double mean_new = 0, meanr_new = 0, meani_new = 0;
	string name_input_directory = "02_05_2025/";
	string name_output_directory = "results/";
	string name_input_file = "gauge_obs2277865449";
	string name_output_file = "autocorr_check_" + name_input_file;
	string line;
	string name_image_poly = "var_poly_" + name_output_file;
	string name_image_polyre = "var_polyre_" + name_output_file;
	string name_image_polyim = "var_polyim_" + name_output_file;
	name_input_file = name_input_directory + name_input_file + ".txt";
	name_output_file = name_output_directory + name_output_file + ".txt";
	name_image_poly = name_image_poly + ".png";
	name_image_polyre = name_image_polyre + ".png";
	name_image_polyim = name_image_polyim + ".png";
	vector <int> n_sub;//will contain N°elements in each subset that we consider
	vector<double> y, yr, yi;
	vector<double> var_new, varr_new, vari_new;
	
	ofstream output_file; //declaration of output file
	output_file.open(name_output_file);
	if (!output_file) {
		cout << "Error opening output file" << endl;
	}
	ifstream input_file; //declaration of output file
	input_file.open(name_input_file);
	if (!input_file) {
		cout << "Error opening input file" << endl;
	}

	for (int i = 0; i < skipLines; i++) {
		if (!getline(input_file, line)) {
			cerr << "Error: there are less than " << skipLines << " lines in the file." << endl;
		}
	}

	output_file << "# N°elements in each subset \t var_m" << endl;
	cout << "# N°elements in each subset \t var_m" << endl;

	while (input_file >> value_tmp >> value_tmp >> value_tmp >> value_tmp >> poly_re >> poly_im) {//DEVI CORREGGEERE ERRORI
		index++;
		poly = poly_re * poly_re + poly_im * poly_im;//poly = |poly|^2 = poly_re^2 + poly_im^2;
		y.push_back(poly);
		yr.push_back(poly_re);
		yi.push_back(poly_im);
	}

	n_sub_max = floor(n_sub_ratio * y.size());

	index = 0;
	for (int dim_block = 1; dim_block <= n_sub_max; dim_block++) {
		blocking_faster(&mean, &var_m, y, dim_block);
		blocking_faster(&mean_re, &var_im, yi, dim_block);
		blocking_faster(&mean_re, &var_re, yr, dim_block);
		n_sub.push_back(dim_block);
		var_new.push_back(var_m);
		varr_new.push_back(var_re);
		vari_new.push_back(var_im);
		output_file << dim_block << "\t" << var_m << "\t" << var_re << "\t" << var_im << endl;
		cout << dim_block << "\t" << var_m << "\t" << var_re << "\t" << var_im << endl;
	}

	vector <double> n_sub_d;
	copy(n_sub.begin(), n_sub.end(), std::back_inserter(n_sub_d));

	plot_points(n_sub_d, var_new, name_image_poly);
	plot_points(n_sub_d, varr_new, name_image_polyre);
	plot_points(n_sub_d, vari_new, name_image_polyim);

	output_file.close();
	input_file.close();

	return 0;
}


//-----------------------------------------------------------------
//ROOT MACRO TO GRAPH:

void plot_points(vector<double>& x, vector<double>& y, string name_image) {
	
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

	g->Draw("ALP");

	//LATEX: I add LaTeX titles and axis labels:
	TLatex latex;
	latex.SetNDC(); //sets the use of Normalized Device Coordinates (NDC).
	latex.SetTextSize(0.05); //changes text size for title
	latex.DrawLatex(0.065, 0.94, "Sample variance from the mean from blocking techinque:");
	latex.SetTextSize(0.04); //changes text size for axis labels
	latex.DrawLatex(0.4, 0.03, "block dimension");
	latex.SetTextAngle(90);
	latex.DrawLatex(0.04, 0.20, "Sample variance from the mean");


	//SAVE: I save the canvas as an image
	canvas->SaveAs(("results/" + name_image).c_str());

	//DELETE:
	delete g;
	delete canvas;
}
