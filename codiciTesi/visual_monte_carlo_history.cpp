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

	int skipLines_file_list = 1, skipLines = 1; //= number of lines to skip while reading input file;
	vector<double> time1, time2, polyr_vec, polyi_vec, reff_vec, imff_vec; //to contain data
	int index = 0, conf_id;
	double value_tmp, poly_re, poly_im, ff_re, ff_im;
	string line;
	vector<string> directories;
	vector<string> gauge_files;
	vector<string> fermion_files;
	string polyre_image;
	string polyim_image;
	string y_name1;
	string y_name2;
	string reff_image;
	string imff_image;
	string title;
	string name_file_list = "11_05_2025/file_list_pt2.txt";
	
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
		string dir, gauge, ferm;
		if (iss >> dir >> gauge >> ferm) {
			directories.push_back(dir);
			gauge_files.push_back(gauge);
			fermion_files.push_back(ferm);
		}
		else {
			cerr << "Poorly formatted line: " << line << endl;
		}
	}

	file_list.close();

	size_t pos;
	string name_tmp;
	for (int ii = 0; ii < directories.size(); ii++) {
		
		y_name1 = "Re{P}";
		y_name2 = "Im{P}";

		ifstream input_file; //declaration of input file
		input_file.open(directories[ii] + gauge_files[ii]);
		if (!input_file) {
			cout << "Error opening " << gauge_files[ii] << endl;
			return 1;
		}
		
		pos = gauge_files[ii].find_last_of(".");
		if (pos != std::string::npos) {
			gauge_files[ii] = gauge_files[ii].substr(0, pos); //I remove extension using substr
		}
		title = "Monte Carlo History " + gauge_files[ii] + ":";
		polyre_image = y_name1 + "_" + gauge_files[ii] + ".png";
		polyim_image = y_name2 + "_" + gauge_files[ii] + ".png";

		for (int i = 0; i < skipLines; i++) {
			if (!getline(input_file, line)) {
				cerr << "Error: there are less than " << skipLines << " lines in the file: " << gauge_files[ii] << endl;
				return 1;
			}
		}

		while (getline(input_file, line)) {
			if (line.find_first_not_of(" \t\r\n") == string::npos) continue;//skips empty lines or lines composed only of spaces/tabs

			istringstream iss(line);
			if (iss >> conf_id >> value_tmp >> value_tmp >> value_tmp >> poly_re >> poly_im) {
				time1.push_back(static_cast<double>(conf_id));
				polyr_vec.push_back(poly_re);
				polyi_vec.push_back(poly_im);
			}
			else {
				cerr << "Skipped gauge line (badly formatted): " << line << endl;
			}
		}

		input_file.close();

		//GRAPHIC REPRESENTATION:
		plot_points(time1, polyr_vec, polyre_image, title, y_name1, 0.47);
		cout << ii << ": Re{P} image DONE!" << endl;
		plot_points(time1, polyi_vec, polyim_image, title, y_name2, 0.47);
		cout << ii << ": Im{P} image DONE!" << endl;
		time1.clear();
		polyr_vec.clear();
		polyi_vec.clear();
	}

	for (int ii = directories.size(); ii < 2*directories.size(); ii++) {

		index = ii - directories.size();
		y_name1 = "Re{ff}";
		y_name2 = "Im{ff}";

		ifstream input_file; //declaration of input file
		input_file.open(directories[index] + fermion_files[index]);
		if (!input_file) {
			cerr << "Error opening " << fermion_files[index] << endl;
			return 1;
		}

		pos = fermion_files[index].find_last_of(".");
		if (pos != std::string::npos) {
			fermion_files[index] = fermion_files[index].substr(0, pos); //I remove extension using substr
		}
		title = "Monte Carlo History " + fermion_files[index] + ":";
		reff_image = y_name1 + "_" + fermion_files[index] + ".png";
		imff_image = y_name2 + "_" + fermion_files[index] + ".png";

		for (int i = 0; i < skipLines; i++) {
			if (!getline(input_file, line)) {
				cerr << "Error: there are less than " << skipLines << " lines in the file" << fermion_files[index] << endl;
				return 1;
			}
		}

		while (getline(input_file, line)) {
			istringstream iss(line);
			if (iss >> conf_id >> value_tmp >> value_tmp >> value_tmp >> ff_re >> ff_im) { //read the line only until ff_im: if next elements are omitted isn't interesting.
				if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
				time2.push_back(static_cast<double>(conf_id));
				reff_vec.push_back(ff_re);
				imff_vec.push_back(ff_im);
			}
			else {
				cerr << "Skipped line (badly formatted): " << line << endl;
			}
		}

		input_file.close();

		//GRAPHIC REPRESENTATION:
		plot_points(time2, reff_vec, reff_image, title, y_name1, 0.47);
		cout << index << ": Re{ff} image DONE!" << endl;
		plot_points(time2, imff_vec, imff_image, title, y_name2, 0.47);
		cout << index << ": Im{ff} image DONE!" << endl;
		time2.clear();
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
	g->GetXaxis()->SetLimits(min_x- 100, max_x + 100);//to set x axis range //SetLimits cannot be ignored by ROOT, while SetRangeUser yes.

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
