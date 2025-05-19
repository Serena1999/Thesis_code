/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                    autocorr_check.cpp:                   ****
****     implements the blocking technique and ranges over    ****
****             different block sizes dim_block,             **** 
****    returning the plot of variance vs dim_block for all   ****
****            files specified in file_list_therm.           ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

#include "../library.h"
#include "../root_include.h"

const bool debug_mode = 0;
const bool fermion_mode = 1;
const bool gauge_mode = 1;

//-----------------------------------------------------------------
//DECLARATIONS:

void plot_points(vector<double>& x, vector<double>& y, string name_image, string title);

void process_autocorr_block(
	const string& input_path,
	const string& output_file,
	const string& output_image_mod,
	const string& output_image_re,
	const string& output_image_im,
	const string& title_mod,
	const string& title_re,
	const string& title_im,
	const string& tipology, //fermion/gauge
	const string& first_out_line,
	int index_loop,
	int n_skip_file,
	int n_skip_re,
	int n_skip_im,
	double n_sub_ratio,
	int step_sample
);

//-----------------------------------------------------------------
//MAIN:

int main() {
	bool bool_long = 1; //1 if you want "_long" in the end of images names;
	int skipLines_g = 1, skipLines_f = 1, skipLines_file_list_therm = 1; //= number of lines to skip while reading input file;
	int step_sample_fermion = 10;
	int step_sample_gauge = 1;
	vector <int> n_skip_rep, n_skip_imp, n_skip_reff, n_skip_imff;
	vector <int> n_sub;//will contain N°elements in each subset that we consider
	vector<string> directories;
	vector<string> gauge_files;
	vector<string> fermion_files;
	double n_sub_ratio = 0.5;//=0.5*len(data) -> you can modify this number from 0 to 1;
	string line, word, title1, title2, title3, var_dimblock_poly_image, var_dimblock_polyre_image, var_dimblock_polyim_image, name_output_file;
	string var_dimblock_ff_image, var_dimblock_ffre_image, var_dimblock_ffim_image;
	const string name_file_list_therm = "19_05_2025/file_list_therm.txt";
	const string first_out_line_gauge = "# N°elements in each subset \t var(|<P * P^dag>| ) \t var(Re{ P }) \t var(Im{ P }):";
	const string first_out_line_ferm = "# N°elements in each subset \t var(|<ff * ff^dag>|) \t var(Re{ff}) \t var(Im{ff}):";
	
	ifstream file_list;
	file_list.open(name_file_list_therm);
	if (!file_list) {
		cout << "Error opening file list" << endl;
		return 1;
	}

	for (int i = 0; i < skipLines_file_list_therm; i++) {
		if (!getline(file_list, line)) {
			cerr << "Error: there are less than " << skipLines_file_list_therm << " lines in the file list." << endl;
			return 1;
		}
	}

	while (getline(file_list, line)) {
		istringstream iss(line);
		string dir, gauge, ferm;
		int n_therm_rep, n_therm_imp, n_therm_reff, n_therm_imff;
		if (iss >> dir >> gauge >> n_therm_rep >> n_therm_imp  >> ferm >> n_therm_reff >> n_therm_imff) {
			directories.push_back(dir);
			gauge_files.push_back(gauge);
			fermion_files.push_back(ferm);
			n_skip_rep.push_back(n_therm_rep / step_sample_gauge);
			n_skip_imp.push_back(n_therm_imp / step_sample_gauge);
			n_skip_reff.push_back(n_therm_reff / step_sample_fermion);
			n_skip_imff.push_back(n_therm_imff / step_sample_fermion);
		}
		else {
			cerr << "Poorly formatted line: " << line << endl;
		}
	}

	file_list.close();

	size_t pos;
	string name_tmp;
	for (int ii = 0; ii < directories.size(); ii++) {

		if (gauge_mode) {
			title1 = "|<P * P^dag>|";
			title2 = "Var(Re{P})";
			title3 = "Var(Im{P})";

			pos = gauge_files[ii].find_last_of(".");
			if (pos != std::string::npos) {
				name_tmp = gauge_files[ii].substr(0, pos); //I remove extension using substr
			}

			if (bool_long) {
				var_dimblock_poly_image = "modP_" + name_tmp + "_long" + ".png";
				var_dimblock_polyre_image = title2 + "_" + name_tmp + "_long" + ".png";
				var_dimblock_polyim_image = title3 + "_" + name_tmp + "_long" + ".png";
			}
			else {
				var_dimblock_poly_image = "modP_" + name_tmp + ".png";
				var_dimblock_polyre_image = title2 + "_" + name_tmp + ".png";
				var_dimblock_polyim_image = title3 + "_" + name_tmp + ".png";
			}
			title1 = title1 + " " + name_tmp + ":";
			title2 = title2 + " " + name_tmp + ":";
			title3 = title3 + " " + name_tmp + ":";

			process_autocorr_block(
				directories[ii] + gauge_files[ii],
				"results/autocorr_check_" + gauge_files[ii],
				var_dimblock_poly_image,
				var_dimblock_polyre_image,
				var_dimblock_polyim_image,
				title1,
				title2,
				title3,
				"gauge", //fermion/gauge
				first_out_line_gauge,
				ii,
				skipLines_g,
				n_skip_rep[ii],
				n_skip_imp[ii],
				n_sub_ratio,
				step_sample_gauge
			);
		}
		if (fermion_mode) {
			title1 = "|<ff * ff^dag>|";
			title2 = "Var(Re{ff})";
			title3 = "Var(Im{ff})";

			pos = fermion_files[ii].find_last_of(".");
			if (pos != std::string::npos) {
				name_tmp = fermion_files[ii].substr(0, pos); //I remove extension using substr
			}

			if (bool_long) {
				var_dimblock_ff_image = "modff_" + name_tmp + "_long" + ".png";
				var_dimblock_ffre_image = title2 + "_" + name_tmp + "_long" + ".png";
				var_dimblock_ffim_image = title3 + "_" + name_tmp + "_long" + ".png";
			}
			else {
				var_dimblock_ff_image = "modff_" + name_tmp + ".png";
				var_dimblock_ffre_image = title2 + "_" + name_tmp + ".png";
				var_dimblock_ffim_image = title3 + "_" + name_tmp + ".png";
			}
			title1 = title1 + " " + name_tmp + ":";
			title2 = title2 + " " + name_tmp + ":";
			title3 = title3 + " " + name_tmp + ":";

			process_autocorr_block(
				directories[ii] + fermion_files[ii],
				"results/autocorr_check_" + fermion_files[ii],
				var_dimblock_ff_image,
				var_dimblock_ffre_image,
				var_dimblock_ffim_image,
				title1,
				title2,
				title3,
				"fermion", //fermion/gauge
				first_out_line_ferm,
				ii,
				skipLines_f,
				n_skip_reff[ii],
				n_skip_imff[ii],
				n_sub_ratio,
				step_sample_fermion
			);
		}
		
	}

	return 0;
}

//-----------------------------------------------------------------
//FUNCTION DEFINITION:

void process_autocorr_block(
	const string& input_path,
	const string& output_path,
	const string& output_image_mod,
	const string& output_image_re,
	const string& output_image_im,
	const string& title_mod,
	const string& title_re,
	const string& title_im,
	const string& tipology, //fermion/gauge
	const string& first_out_line,
	int index_loop,
	int n_skip_file,
	int n_skip_re,
	int n_skip_im,
	double n_sub_ratio,
	int step_sample
) {
	if (tipology == "fermion") {
		size_t pos;
		string line;
		double obs, obs_re, obs_im;
		vector <double> y, yi, yr, var_new, varr_new, vari_new;
		vector <int> n_sub;
		double mean, mean_re, mean_im, var_m, var_re, var_im, value_re_tmp = 0, value_im_tmp = 0, value_mod_tmp = 0;
		string value_tmp;
		int conf_id = -999, conf_tmp = 999, n_copy_accum_r = 0, n_copy_accum_i = 0, n_copy_accum_m = 0, flag = 0;
		double counter = 0;

		ifstream input_file; //declaration of input file
		input_file.open(input_path);
		if (!input_file) {
			cout << "Error opening " << input_path << endl;
			return;
		}

		for (int jj = 0; jj < n_skip_file; jj++) {
			if (!getline(input_file, line)) {
				cerr << "Error: there are less than " << n_skip_file << " lines in the file: " << input_path << endl;
				return;
			}
		}

		for (int jj = 0; jj < min(n_skip_re, n_skip_im); jj++) {
			while (conf_id != conf_tmp) {
				if (!getline(input_file, line)) {
					cerr << "Error: there are less than " << min(n_skip_re, n_skip_im) + n_skip_file << " lines in the file: " << input_path << endl;
					return;
				}
				istringstream iss(line);
				if (iss >> conf_id) {
					if (jj == 0) {
						conf_tmp = conf_id;
					}
				}
			}
		}

		if (n_skip_re < n_skip_im) {
			while(counter < (n_skip_im - n_skip_re)){
				if(!getline(input_file, line)) break;
				istringstream iss(line);
				if (iss >> conf_id >> value_tmp >> value_tmp >> value_tmp >> obs_re >> obs_im) {
					if (conf_id != conf_tmp) {
						counter++;
						conf_tmp = conf_id;
						value_re_tmp /= (double) n_copy_accum_r;
						yr.push_back(value_re_tmp * step_sample);
						value_re_tmp = 0;
						n_copy_accum_r = 0;
					}
					n_copy_accum_r++;
					value_re_tmp += obs_re;
				}
				else {
					cerr << "Skipped line (badly formatted) from" << input_path << ": " << line << endl;
				}
			}
		}
		else if (n_skip_im < n_skip_re) {
			while(counter < (n_skip_re - n_skip_im)){
				if (!getline(input_file, line)) break;
				istringstream iss(line);
				if (iss >> conf_id >> value_tmp >> value_tmp >> value_tmp >> obs_re >> obs_im) {
					if (conf_id != conf_tmp) {
						counter++;
						conf_tmp = conf_id;
						value_im_tmp /= (double) n_copy_accum_i;
						yi.push_back(value_im_tmp * step_sample);
						value_im_tmp = 0;
						n_copy_accum_i = 0;
					}
					n_copy_accum_i++;
					value_im_tmp += obs_im;
				}
				else {
					cerr << "Skipped line (badly formatted) from" << input_path << ": " << line << endl;
				}
			}
		}

		int cnt = 0, flag_start = 0;
		if (n_skip_im == n_skip_re) flag_start = 1;
		while (getline(input_file, line)) {
			istringstream iss(line);
			if (iss >> conf_id >> value_tmp >> value_tmp >> value_tmp >> obs_re >> obs_im) {
				if (flag_start) {
					conf_tmp = conf_id;
				}
				flag = 0;
				if (debug_mode) {
					cout << "(conf_id, conf_tmp) = (" << conf_id << ", " << conf_tmp << ")" << endl;
				}
				if (conf_id != conf_tmp) {
					cnt++;
					conf_tmp = conf_id;
					value_re_tmp /= (double) n_copy_accum_r;
					value_im_tmp /= (double) n_copy_accum_i;
					value_mod_tmp /= (double) n_copy_accum_m;
					y.push_back(value_mod_tmp * step_sample * step_sample);
					yr.push_back(value_re_tmp * step_sample);
					yi.push_back(value_im_tmp * step_sample);
					if (debug_mode) {
						cout << "Some acculumulated data in fermion file: " << cnt << "\t" << value_mod_tmp << "\t" << value_re_tmp << "\t" << value_im_tmp << endl;
					}
					value_mod_tmp = 0;
					value_re_tmp = 0;
					value_im_tmp = 0;
					n_copy_accum_m = 0;
					n_copy_accum_r = 0;
					n_copy_accum_i = 0;
				}
				n_copy_accum_m++;
				n_copy_accum_r++;
				n_copy_accum_i++;
				value_mod_tmp += obs_re * obs_re + obs_im * obs_im;//obs = |obs|^2 = obs_re^2 + obs_im^2;
				value_re_tmp += obs_re;
				value_im_tmp += obs_im;
			}
			else {
				flag = 1;
				cerr << "Skipped line (badly formatted) from" << input_path << ": " << line << endl;
			}
		}
		if ((conf_id == conf_tmp) && (!flag)) {
			value_re_tmp /= n_copy_accum_r;
			value_im_tmp /= n_copy_accum_i;
			value_mod_tmp /= n_copy_accum_m;
			y.push_back(value_mod_tmp * step_sample * step_sample);
			yr.push_back(value_re_tmp * step_sample);
			yi.push_back(value_im_tmp * step_sample);
		}

		input_file.close();

		ofstream output_file; //declaration of output file
		output_file.open(output_path);
		if (!output_file) {
			cout << "Error opening output file " << output_path << endl;
		}

		int n_sub_max = floor(n_sub_ratio * y.size());
		if (debug_mode) {
			cout << "n_sub_max = " << n_sub_max << endl;
			cout << "size of y array = " << y.size() << endl;
		}

		output_file << first_out_line << endl;
		for (int dim_block = 1; dim_block <= n_sub_max; dim_block++) {
			if(!blocking_faster(&mean, &var_m, y, dim_block)) break;
			if(!blocking_faster(&mean_re, &var_re, yr, dim_block)) break;
			if(!blocking_faster(&mean_im, &var_im, yi, dim_block)) break;
			n_sub.push_back(dim_block);
			var_new.push_back(var_m);
			varr_new.push_back(var_re);
			vari_new.push_back(var_im);
			output_file << dim_block << "\t" << var_m << "\t" << var_re << "\t" << var_im << endl;
		}
		
		vector <double> n_sub_d;
		copy(n_sub.begin(), n_sub.end(), std::back_inserter(n_sub_d));

		//GRAPHIC REPRESENTATION:
		plot_points(n_sub_d, var_new, output_image_mod, title_mod);
		cout << index_loop << ": " << tipology << " Obs image DONE!" << endl;
		plot_points(n_sub_d, varr_new, output_image_re, title_re);
		cout << index_loop << ": " << tipology << " Re{obs} image DONE!" << endl;
		plot_points(n_sub_d, vari_new, output_image_im, title_im);
		cout << index_loop << ": " << tipology << " Im{obs} image DONE!" << endl;

		output_file.close();

		y.clear();
		yi.clear();
		yr.clear();
		n_sub.clear();
		n_sub_d.clear();
		var_new.clear();
		varr_new.clear();
		vari_new.clear();

		return;
	}
	size_t pos;
	string line;
	double value_tmp, obs, obs_re, obs_im;
	vector <double> y, yi, yr, var_new, varr_new, vari_new;
	vector <int> n_sub;
	double mean, mean_re, mean_im, var_m, var_re, var_im;

	ifstream input_file; //declaration of input file
	input_file.open(input_path);
	if (!input_file) {
		cout << "Error opening " << input_path << endl;
		return;
	}

	if (!getline(input_file, line)) {
		cerr << "Error: there are less than " << min(n_skip_re, n_skip_im) << " lines in the file: " << input_path << endl;
		return;
	}

	for (int jj = 0; jj < min(n_skip_re, n_skip_im); jj++) {
		if (!getline(input_file, line)) {
			cerr << "Error: there are less than " << min(n_skip_re, n_skip_im) << " lines in the file: " << input_path << endl;
			return;
		}
	}

	if (n_skip_re < n_skip_im) {
		for (int jj = 0; jj < (n_skip_im - n_skip_re); jj++) {
			getline(input_file, line);
			istringstream iss(line);
			if (iss >> value_tmp >> value_tmp >> value_tmp >> value_tmp >> obs_re >> obs_im) {
				yr.push_back(obs_re * step_sample);
			}
			else {
				cerr << "Skipped line (badly formatted) from" << input_path << ": " << line << endl;
			}
		}
	}
	else if (n_skip_im < n_skip_re) {

		for (int jj = 0; jj < (n_skip_re - n_skip_im); jj++) {
			getline(input_file, line);
			istringstream iss(line);
			if (iss >> value_tmp >> value_tmp >> value_tmp >> value_tmp >> obs_re >> obs_im) {
				yi.push_back(obs_im * step_sample);
			}
			else {
				cerr << "Skipped line (badly formatted) from" << input_path << ": " << line << endl;
			}
		}
	}

	while (getline(input_file, line)) {
		istringstream iss(line);
		if (iss >> value_tmp >> value_tmp >> value_tmp >> value_tmp >> obs_re >> obs_im) {
			obs = obs_re * obs_re + obs_im * obs_im;//obs = |obs|^2 = obs_re^2 + obs_im^2;
			y.push_back(obs * step_sample * step_sample);
			yr.push_back(obs_re * step_sample);
			yi.push_back(obs_im * step_sample);
		}
		else {
			cerr << "Skipped line (badly formatted) from" << input_path << ": " << line << endl;
		}
	}

	input_file.close();

	ofstream output_file; //declaration of output file
	output_file.open(output_path);
	if (!output_file) {
		cout << "Error opening output file " << output_path << endl;
	}

	int n_sub_max = floor(n_sub_ratio * y.size());

	output_file << first_out_line << endl;
	for (int dim_block = 1; dim_block <= n_sub_max; dim_block++) {
		if(!blocking_faster(&mean, &var_m, y, dim_block)) break;
		if(!blocking_faster(&mean_re, &var_re, yr, dim_block)) break;
		if(!blocking_faster(&mean_im, &var_im, yi, dim_block)) break;
		n_sub.push_back(dim_block);
		var_new.push_back(var_m);
		varr_new.push_back(var_re);
		vari_new.push_back(var_im);
		output_file << dim_block << "\t" << var_m << "\t" << var_re << "\t" << var_im << endl;
	}

	vector <double> n_sub_d;
	copy(n_sub.begin(), n_sub.end(), std::back_inserter(n_sub_d));

	//GRAPHIC REPRESENTATION:
	plot_points(n_sub_d, var_new, output_image_mod, title_mod);
	cout << index_loop << ": " << tipology << " Obs image DONE!" << endl;
	plot_points(n_sub_d, varr_new, output_image_re, title_re);
	cout << index_loop << ": " << tipology << " Re{obs} image DONE!" << endl;
	plot_points(n_sub_d, vari_new, output_image_im, title_im);
	cout << index_loop << ": " << tipology << " Im{obs} image DONE!" << endl;

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
	latex.DrawLatex(0.065, 0.94, "Sample variance from the mean from blocking techique:");
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
