/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                    autocorr_check.cpp:                   ****
****               using the blocking technique               ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

#include "../library.h"
#include "../root_include.h"

void plot_points(vector<double>& x, vector<double>& y, string name_image, string title);

//-----------------------------------------------------------------
//MAIN:

int main() {
	bool bool_long = 1; //1 if you want "_long" in the end of images names;
	int Nt = 8; //BE CAREFUL TO CHOOSE IT WELL;
	int skipLines = 1, skipLines_file_list_therm = 1; //= number of lines to skip while reading input file;
	vector <int> n_skip_rep, n_skip_imp, n_skip_reff, n_skip_imff;
	vector <int> n_sub;//will contain N°elements in each subset that we consider
	vector<double> y, yr, yi;
	vector<double> var_new, varr_new, vari_new;
	vector<string> directories;
	vector<string> gauge_files;
	vector<string> fermion_files;
	int n_sub_max;
	int N_tmp;
	double n_sub_ratio = 0.5;//0.03;//=0.5*len(data) -> you can modify this number from 0 to 1;
	double mpi = 800; //MeV //BE CAREFUL TO CHOOSE IT WELL;
	double poly, poly_re, poly_im, ff, ff_re, ff_im, delta, delta_re, delta_im, value_tmp;
	double mean = 0, var_m = 0, mean_re = 0, var_re = 0, mean_im = 0, var_im = 0;
	double mean_new = 0, meanr_new = 0, meani_new = 0;
	string line, word, title1, title2, title3, var_dimblock_poly_image, var_dimblock_polyre_image, var_dimblock_polyim_image, name_output_file;
	string var_dimblock_ff_image, var_dimblock_ffre_image, var_dimblock_ffim_image;
	string name_file_list_therm = "11_05_2025/file_list_therm.txt";
	
	ifstream file_list;
	file_list.open(name_file_list_therm);
	if (!file_list) {
		cout << "Error opening file list" << endl;
		return 1;
	}

	for (int i = 0; i < skipLines; i++) {
		if (!getline(file_list, line)) {
			cerr << "Error: there are less than " << skipLines_file_list_therm << " lines in the file list." << endl;
			return 1;
		}
	}

	while (getline(file_list, line)) {
		istringstream iss(line);
		string dir, gauge, ferm;
		int n_therm_rep, n_therm_imp, n_therm_reff, n_therm_imff, n_copy;
		if (iss >> dir >> gauge >> n_therm_rep >> n_therm_imp  >> ferm >> n_therm_reff >> n_therm_imff >> n_copy) {
			directories.push_back(dir);
			gauge_files.push_back(gauge);
			fermion_files.push_back(ferm);
			n_skip_rep.push_back(skipLines+n_therm_rep);
			n_skip_imp.push_back(skipLines + n_therm_imp);
			n_skip_reff.push_back(skipLines + n_therm_reff * n_copy);
			n_skip_imff.push_back(skipLines + n_therm_imff * n_copy);
		}
		else {
			cerr << "Poorly formatted line: " << line << endl;
		}
	}

	file_list.close();

	size_t pos;
	string name_tmp;
	for (int ii = 0; ii < directories.size(); ii++) {
		title1 = "|<P * P^dag>|";
		title2 = "Var(Re{P})";
		title3 = "Var(Im{P})";

		ifstream input_file; //declaration of input file
		input_file.open(directories[ii] + gauge_files[ii]);
		if (!input_file) {
			cout << "Error opening " << gauge_files[ii] << endl;
			return 1;
		}

		name_output_file = "results/autocorr_check_" + gauge_files[ii];
		pos = gauge_files[ii].find_last_of(".");
		if (pos != std::string::npos) {
			gauge_files[ii] = gauge_files[ii].substr(0, pos); //I remove extension using substr
		}
		
		if (bool_long) {
			var_dimblock_poly_image = "modP_" + gauge_files[ii] + "_long" + ".png";
			var_dimblock_polyre_image = title2 + "_" + gauge_files[ii] + "_long" + ".png";
			var_dimblock_polyim_image = title3 + "_" + gauge_files[ii] + "_long" + ".png";
		}
		else {
			var_dimblock_poly_image = "modP_" + gauge_files[ii] + ".png";
			var_dimblock_polyre_image = title2 + "_" + gauge_files[ii] + ".png";
			var_dimblock_polyim_image = title3 + "_" + gauge_files[ii] + ".png";
		}
		title1 = title1 + " " + gauge_files[ii] + ":";
		title2 = title2 + " " + gauge_files[ii] + ":";
		title3 = title3 + " " + gauge_files[ii] + ":";
		
		for (int jj = 0; jj < min(n_skip_rep[ii], n_skip_imp[ii]); jj++) {
			if (!getline(input_file, line)) {
				cerr << "Error: there are less than " << min(n_skip_rep[ii], n_skip_imp[ii]) << " lines in the file: " << gauge_files[ii] << endl;
				return 1;
			}
		}

		if (n_skip_rep[ii] < n_skip_imp[ii]) {

			for(int jj=0; jj < (n_skip_imp[ii] - n_skip_rep[ii]); jj++){
				getline(input_file, line);
				istringstream iss(line);
				if (iss >> value_tmp >> value_tmp >> value_tmp >> value_tmp >> poly_re >> poly_im) {
					yr.push_back(poly_re);
				}
				else {
					cerr << "Skipped gauge line (badly formatted): " << line << endl;
				}
			}
		}
		else if (n_skip_imp[ii] < n_skip_rep[ii]) {
			
			for (int jj = 0; jj < (n_skip_rep[ii] - n_skip_imp[ii]); jj++) {
				getline(input_file, line);
				istringstream iss(line);
				if (iss >> value_tmp >> value_tmp >> value_tmp >> value_tmp >> poly_re >> poly_im) {
					yi.push_back(poly_im);
				}
				else {
					cerr << "Skipped gauge line (badly formatted): " << line << endl;
				}
			}
		}

		while (getline(input_file, line)) {
			istringstream iss(line);
			if (iss >> value_tmp >> value_tmp >> value_tmp >> value_tmp >> poly_re >> poly_im) {
				poly = poly_re * poly_re + poly_im * poly_im;//poly = |poly|^2 = poly_re^2 + poly_im^2;
				y.push_back(poly);
				yr.push_back(poly_re);
				yi.push_back(poly_im);
			}
			else {
				cerr << "Skipped gauge line (badly formatted): " << line << endl;
			}
		}

		input_file.close();

		ofstream output_file; //declaration of output file
		output_file.open(name_output_file);
		if (!output_file) {
			cout << "Error opening output file " << name_output_file << endl;
		}

		n_sub_max = floor(n_sub_ratio * y.size());

		output_file << "# N°elements in each subset \t var(|<P * P^dag>|) \t var(Re{P}) \t var(Im{P}):" << endl;
		for (int dim_block = 1; dim_block <= n_sub_max; dim_block++) {
			blocking_faster(&mean, &var_m, y, dim_block);
			blocking_faster(&mean_re, &var_re, yr, dim_block);
			blocking_faster(&mean_im, &var_im, yi, dim_block);
			n_sub.push_back(dim_block);
			var_new.push_back(var_m);
			varr_new.push_back(var_re);
			vari_new.push_back(var_im);
			output_file << dim_block << "\t" << var_m << "\t" << var_re << "\t" << var_im << endl;
		}

		vector <double> n_sub_d;
		copy(n_sub.begin(), n_sub.end(), std::back_inserter(n_sub_d));

		//GRAPHIC REPRESENTATION:
		plot_points(n_sub_d, var_new, var_dimblock_poly_image, title1);
		cout << ii << ": P image DONE!" << endl;
		plot_points(n_sub_d, varr_new, var_dimblock_polyre_image, title2);
		cout << ii << ": Re{P} image DONE!" << endl;
		plot_points(n_sub_d, vari_new, var_dimblock_polyim_image, title3);
		cout << ii << ": Re{P} image DONE!" << endl;

		output_file.close();
		

		y.clear();
		yi.clear();
		yr.clear();
		n_sub.clear();
		n_sub_d.clear();
		var_new.clear();
		varr_new.clear();
		vari_new.clear();
	}

	int index = 0;
	for (int ii = directories.size(); ii < 2 * directories.size(); ii++) {
		index = ii - directories.size();
		title1 = "|<ff * ff^dag>|";
		title2 = "Var(Re{ff})";
		title3 = "Var(Im{ff})";

		ifstream input_file; //declaration of input file
		input_file.open(directories[index] + fermion_files[index]);
		if (!input_file) {
			cout << "Error opening " << fermion_files[index] << endl;
			return 1;
		}

		name_output_file = "results/autocorr_check_" + fermion_files[index];
		pos = fermion_files[index].find_last_of(".");
		if (pos != std::string::npos) {
			fermion_files[index] = fermion_files[index].substr(0, pos); //I remove extension using substr
		}
		if (bool_long) {
			var_dimblock_ff_image = "modff_" + fermion_files[index] + "_long" + ".png";
			var_dimblock_ffre_image = title2 + "_" + fermion_files[index] + "_long" + ".png";
			var_dimblock_ffim_image = title3 + "_" + fermion_files[index] + "_long" + ".png";
		}
		else {
			var_dimblock_ff_image = "modff_" + fermion_files[index] + ".png";
			var_dimblock_ffre_image = title2 + "_" + fermion_files[index] + ".png";
			var_dimblock_ffim_image = title3 + "_" + fermion_files[index] + ".png";
		}
		title1 = title1 + " " + fermion_files[index] + ":";
		title2 = title2 + " " + fermion_files[index] + ":";
		title3 = title3 + " " + fermion_files[index] + ":";

		for (int jj = 0; jj < min(n_skip_reff[index], n_skip_imff[index]); jj++) {
			if (!getline(input_file, line)) {
				cerr << "Error: there are less than " << min(n_skip_reff[index], n_skip_imff[index]) << " lines in the file: " << fermion_files[index] << endl;
				return 1;
			}
		}

		if (n_skip_rep[index] < n_skip_imp[index]) {

			for (int jj = 0; jj < (n_skip_imp[index] - n_skip_rep[index]); jj++) {
				getline(input_file, line);
				istringstream iss(line);
				if (iss >> value_tmp >> value_tmp >> value_tmp >> value_tmp >> ff_re >> ff_im) {
					yr.push_back(ff_re);
				}
				else {
					cerr << "Skipped fermion line (badly formatted): " << line << endl;
				}
			}
		}
		else if (n_skip_imff[index] < n_skip_reff[index]) {
			for (int jj = 0; jj < (n_skip_rep[index] - n_skip_imp[index]); jj++) {
				getline(input_file, line);
				istringstream iss(line);
				if (iss >> value_tmp >> value_tmp >> value_tmp >> value_tmp >> ff_re >> ff_im) {
					yi.push_back(ff_im);
				}
				else {
					cerr << "Skipped fermion line (badly formatted): " << line << endl;
				}
			}
		}

		while (getline(input_file, line)) {
			istringstream iss(line);
			if (iss >> value_tmp >> value_tmp >> value_tmp >> value_tmp >> ff_re >> ff_im) {
				ff = ff_re * ff_re + ff_im * ff_im; //ff = |ff|^2 = ff_re^2 + ff_im^2;
				y.push_back(ff);
				yr.push_back(ff_re);
				yi.push_back(ff_im);
			}
			else {
				cerr << "Skipped fermion line (badly formatted): " << line << endl;
			}
		}

		input_file.close();

		ofstream output_file; //declaration of output file
		output_file.open(name_output_file);
		if (!output_file) {
			cout << "Error opening output file " << name_output_file << endl;
		}

		n_sub_max = floor(n_sub_ratio * y.size());

		output_file << "# N°elements in each subset \t var(|<ff * ff^dag>|) \t var(Re{ff}) \t var(Im{ff}):" << endl;
		for (int dim_block = 1; dim_block <= n_sub_max; dim_block++) {
			blocking_faster(&mean, &var_m, y, dim_block);
			blocking_faster(&mean_re, &var_re, yr, dim_block);
			blocking_faster(&mean_im, &var_im, yi, dim_block);
			n_sub.push_back(dim_block);
			var_new.push_back(var_m);
			varr_new.push_back(var_re);
			vari_new.push_back(var_im);
			output_file << dim_block << "\t" << var_m << "\t" << var_re << "\t" << var_im << endl;
		}

		vector <double> n_sub_d;
		copy(n_sub.begin(), n_sub.end(), std::back_inserter(n_sub_d));

		//GRAPHIC REPRESENTATION:
		plot_points(n_sub_d, var_new, var_dimblock_ff_image, title1);
		cout << index << ": ff image DONE!" << endl;
		plot_points(n_sub_d, varr_new, var_dimblock_ffre_image, title2);
		cout << index << ": Re{ff} image DONE!" << endl;
		plot_points(n_sub_d, vari_new, var_dimblock_ffim_image, title3);
		cout << index << ": Im{ff} image DONE!" << endl;

		output_file.close();


		y.clear();
		yi.clear();
		yr.clear();
		n_sub.clear();
		n_sub_d.clear();
		var_new.clear();
		varr_new.clear();
		vari_new.clear();
	}


	return 0;
}


//-----------------------------------------------------------------
//ROOT MACRO TO GRAPH:

void plot_points(vector<double>& x, vector<double>& y, string name_image, string title) {
	
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
	latex.DrawLatex(0.04, 0.20, title.c_str());


	//SAVE: I save the canvas as an image
	canvas->SaveAs(("results/" + name_image).c_str());

	//DELETE:
	delete g;
	delete canvas;
}
