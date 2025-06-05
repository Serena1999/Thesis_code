/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                    autocorr_check.cpp:                   ****
****     implements the blocking technique and ranges over    ****
****             different block sizes dim_block,             **** 
****    returning the plot of variance vs dim_block for all   ****
****            files specified in file_list_therm            ****
****          FORMAT FILE: 2 COLUMNS (CONF_ID, VALUE)         ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

#include "../library.h"
#include "../root_include.h"

const bool debug_mode = 0;

//-----------------------------------------------------------------
//DECLARATIONS:

void plot_points(
	vector<double>& x,
	vector<double>& y, 
	string name_image, 
	string title);

void process_autocorr_block(
	const string& input_path,
	const string& output_path,
	const string& output_image,
	const string& title,
	const string& first_out_line,
	int index_loop,
	int n_skip_file,
	int n_skip,
	double n_sub_ratio
);

//-----------------------------------------------------------------
//MAIN:

int main() {
	bool bool_long = 0; //1 if you want "_long" in the end of images names;
	int skipLines_file = 0, skipLines_file_list_therm = 1; //= number of lines to skip while reading input file;
	//int step_sample_fermion = 10;
	//int step_sample_gauge = 1;
	vector <int> n_skip;
	vector <int> n_sub;//will contain N°elements in each subset that we consider
	vector<string> directories, files, titles;
	double n_sub_ratio = 0.05;//=0.5*len(data) -> you can modify this number from 0 to 1;
	string line, word, title, var_dimblock_poly_image, var_dimblock_polyre_image, var_dimblock_polyim_image, name_output_file;
	string image;
	const string name_file_list_therm = "11_05_2025/data_square/file_list_therm.txt";
	const string first_out_line = "# N°elements in each subset \t	mean(value^2) \t var(value^2)";
	
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
//ATTENTA!
	while (getline(file_list, line)) {
		istringstream iss(line);
		string dir, file, title_tmp;//title_tmp = title of the graph;
		int n_therm, step_sample;//n_therm: da scartare per la termalizzazione, step_sample: le configurazioni sono numerate a passi di step_sample
		if (iss >> dir >> file >> title_tmp >> n_therm >> step_sample) {
			directories.push_back(dir);
			files.push_back(file);
			titles.push_back(title_tmp);
			n_skip.push_back(n_therm / step_sample);
		}
		else {
			cerr << "Poorly formatted line: " << line << endl;
		}
	}

	file_list.close();

	size_t pos;
	string name_tmp;
	for (int ii = 0; ii < directories.size(); ii++) {
		title = titles[ii];
		pos = files[ii].find_last_of(".");
		if (pos != std::string::npos) {
			name_tmp = files[ii].substr(0, pos); //I remove extension using substr
		}
		if (bool_long) {
			image = "autocorr_" + name_tmp + "_long" + ".png";
		}
		else {
			image = "autocorr_" + name_tmp + ".png";
		}
		
		title = title + " " + name_tmp + ":";

		process_autocorr_block(
			directories[ii] + files[ii],
			"results/autocorr_check_" + files[ii],
			image,
			title,
			first_out_line,
			ii,
			skipLines_file,
			n_skip[ii],
			n_sub_ratio
			);

	}

	return 0;
}

//-----------------------------------------------------------------
//FUNCTION DEFINITION:

void process_autocorr_block(
	const string& input_path,
	const string& output_path,
	const string& output_image,
	const string& title,
	const string& first_out_line,
	int index_loop,
	int n_skip_file,
	int n_skip,
	double n_sub_ratio
) {
	size_t pos;
	string line;
	double obs;
	vector <double> y, var_new;
	vector <int> n_sub;
	double mean, var_m, value_tmp = 0;
	int conf_id = -999, conf_tmp = 999, n_copy_accum = 0, flag = 0, flag_init = 0;
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

	for (int jj = 0; jj < n_skip; jj++) {
		while (conf_id != conf_tmp) {
			if (!getline(input_file, line)) {
				cerr << "Error: there are less than " << n_skip + n_skip_file << " lines in the file: " << input_path << endl;
				return;
			}
			istringstream iss(line);
			if (iss >> conf_id) {
				if (jj == 0) {
					conf_tmp = conf_id;
					flag_init = 1;
				}
			}
		}
	}

	int cnt = 0;
	while (getline(input_file, line)) {
		istringstream iss(line);
		if (iss >> conf_id >> obs) {
			if (!flag_init) {
				conf_tmp = conf_id;
				flag_init = 1;
			}
			flag = 0;
			if ((debug_mode) && (0)) {
				cout << "(conf_id, conf_tmp) = (" << conf_id << ", " << conf_tmp << ")" << endl;
			}
			if (conf_id != conf_tmp) {
				cnt++;
				conf_tmp = conf_id;
				if (debug_mode) {
					if (n_copy_accum == 0) cout << "NAN HERE 3!" << endl;
				}
				value_tmp /= (double) n_copy_accum;
				y.push_back(value_tmp);
				value_tmp = 0;
				n_copy_accum = 0;
			}
			n_copy_accum++;
			value_tmp += obs;
		}
		else {
			flag = 1;
			cerr << "Skipped line (badly formatted) from" << input_path << ": " << line << endl;
		}
	}
	if ((conf_id == conf_tmp) && (!flag)) {
		value_tmp /= (double) n_copy_accum;
		y.push_back(value_tmp);
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
	output_file << setprecision(numeric_limits<double>::max_digits10);
	for (int dim_block = 1; dim_block <= n_sub_max; dim_block++) {
		if(!blocking_faster(&mean, &var_m, y, dim_block)) break;
		n_sub.push_back(dim_block);
		var_new.push_back(var_m);
		output_file << dim_block << "\t" << mean << "\t" << var_m << endl;
	}
		
	vector <double> n_sub_d;
	copy(n_sub.begin(), n_sub.end(), std::back_inserter(n_sub_d));

	//GRAPHIC REPRESENTATION:
	plot_points(n_sub_d, var_new, output_image, title);
	cout << index_loop << ": " << " Obs image DONE!" << endl;

	output_file.close();

	y.clear();
	n_sub.clear();
	n_sub_d.clear();
	var_new.clear();
}


//-----------------------------------------------------------------
//ROOT MACRO TO GRAPH:

void plot_points(
	vector<double>& x, 
	vector<double>& y, 
	string name_image, 
	string title
) {
	
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
