/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                     fit_polyffT.cpp:                     ****
****      visualization and fit of the mean values ​​of the     **** 
****    observables (with corresponding standard deviation)   ****
****             as a function of the temperature             ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

//-----------------------------------------------------------------
//HEADERS && LIBRARY:

#include "../library.h"
#include "../root_include.h"

//-----------------------------------------------------------------
//VARIABLES TO SET:

const string tipology = "gauge"; //gauge/fermion, CHOOSABLE --> to do the gauge/fermion observables graph
#define CHOOSE_FIT_FUNCTION 6 //0 for polynomial, 1 for arctg, 2 for logistic function, 3 for Hill function, 4 for arctg+linear, 5 for sigmoid, 6 for personalized sigmoid;
const bool bool_choose_at_eye = 0; //0 if you want an automatic set of parameters, 1 if you want to impose them by hand;
// -> if 1, modify the corrisponding if condition in par_estimate function to choose parameters;
bool bool_enlarge = 1;//0 if you want to use the original draws, 1 if you want that y-values are multiplied by a the following factor:
double enlarge_factor = 100;//factor to multiply the draws if bool_enlarge = 0;

//-----------------------------------------------------------------
//FIT && ESTIMATE OF PARAMETERS FUNCTIONS: 

#if CHOOSE_FIT_FUNCTION == 0
	
	double fit_function(//T_c = flex point = - p[2]/(3*p[3])
		double* y, //temperatures
		double* p //parameters
	) {
		return p[0] + p[1] * y[0] + p[2] * y[0] * y[0] + p[3] * y[0] * y[0] * y[0];
	};

	const int n_par_fit = 4;

	bool par_estimate(//considering polynomial as y(x) = p[0] + p[1]*x + p[2] * x^2 + p[3] * x^3;
		//return 0 if success, 1 if not;
		const vector<double>& x, //temperatures
		const vector<double>& y, //observables
		vector<double>& p //parameters
	) {

		if (bool_choose_at_eye || (x.size() < n_par_fit)) {//TO CHANGE THE FOLLOWING FOR "CHOOSE BY EYE" SETTING
			p[0] = 1;
			p[1] = 1;
			p[2] = 1;
			p[3] = 1;

			if (x.size() < n_par_fit) {
				cerr << "Not enough points for estimate." << endl;
				return 1;
			}

			return 0;
		}

		int step = x.size() / 4;

		//using Cramer method to solve system of equations:
		double detM = (x[0] - x[step]) * (x[0] - x[2 * step]) * (x[step] - x[2 * step]) * (x[0] - x[3 * step]) * (x[step] - x[3 * step]) * (x[2 * step] - x[3 * step]);
		
		if (detM == 0) {
			cout << "detM = 0, so I use {1,1,1,1} as init" << endl;
			p[0] = 1;
			p[1] = 1;
			p[2] = 1;
			p[3] = 1;
			return;
		}
		
		for (int ii = 0; ii < p.size(); ++ii) {
			p[ii] = 0;
		}
		
		for (int ii = 0; ii < (x.size() - 4); ++ii) {
			p[0] += -x[ii] * x[ii + 1] * (x[ii] - x[ii + 1]) * (x[ii + 2] * y[ii + 3] * (x[ii] - x[ii + 2]) * (x[ii + 2] - x[ii + 1]) + x[ii + 3] * y[ii + 2] * (x[ii] - x[ii + 3]) * (x[ii + 1] - x[ii + 3])) + x[ii] * x[ii + 2] * x[ii + 3] * y[ii] * (x[ii] - x[ii + 2]) * (x[ii] - x[ii + 3]) * (x[ii + 2] - x[ii + 3]) - x[ii + 1] * x[ii + 2] * x[ii + 3] * y[0] * (x[ii + 1] - x[ii + 2]) * (x[ii + 1] - x[ii + 3]) * (x[ii + 2] - x[ii + 3]);
			p[1] += y[ii] * (pow(x[ii], 3) * (pow(x[ii + 3], 2) - pow(x[ii + 2], 2)) + pow(x[ii], 2) * (pow(x[ii + 2], 3) - pow(x[ii + 3], 3)) + pow(x[ii + 2], 2) * pow(x[ii + 3], 2) * (x[ii + 3] - x[ii + 2])) + (x[ii] - x[ii + 1]) * (y[ii + 3] * (x[ii] - x[ii + 2]) * (x[ii + 2] - x[ii + 1]) * (x[ii] * (x[ii + 1] + x[ii + 2]) + x[ii + 1] * x[ii + 2]) + y[ii + 2] * (x[ii] - x[ii + 3]) * (x[ii + 1] - x[ii + 3]) * (x[ii] * (x[ii + 1] + x[ii + 3]) + x[ii + 1] * x[ii + 3])) + y[0] * (pow(x[ii + 1], 3) * (pow(x[ii + 2], 2) - pow(x[ii + 3], 2)) + pow(x[ii + 1], 2) * (pow(x[ii + 3], 3) - pow(x[ii + 2], 3)) + pow(x[ii + 2], 2) * pow(x[ii + 3], 2) * (x[ii + 2] - x[ii + 3]));
			p[2] += pow(x[ii], 3) * (-x[ii + 1]) * y[ii + 2] + pow(x[ii], 3) * x[ii + 1] * y[ii + 3] + y[ii] * (pow(x[ii], 3) * (x[ii + 2] - x[ii + 3]) + x[ii] * (pow(x[ii + 3], 3) - pow(x[ii + 2], 3)) + x[ii + 2] * x[ii + 3] * (pow(x[ii + 2], 2) - pow(x[ii + 3], 2))) - pow(x[ii], 3) * x[ii + 2] * y[ii + 3] + pow(x[ii], 3) * x[ii + 3] * y[ii + 2] + x[ii] * pow(x[ii + 1], 3) * y[ii + 2] - x[ii] * pow(x[ii + 1], 3) * y[ii + 3] + x[ii] * pow(x[ii + 2], 3) * y[ii + 3] - x[ii] * pow(x[ii + 3], 3) * y[ii + 2] + y[0] * (pow(x[ii + 1], 3) * (x[ii + 3] - x[ii + 2]) + x[ii + 1] * (pow(x[ii + 2], 3) - pow(x[ii + 3], 3)) - pow(x[ii + 2], 3) * x[ii + 3] + x[ii + 2] * pow(x[ii + 3], 3)) + pow(x[ii + 1], 3) * x[ii + 2] * y[ii + 3] - pow(x[ii + 1], 3) * x[ii + 3] * y[ii + 2] - x[ii + 1] * pow(x[ii + 2], 3) * y[ii + 3] + x[ii + 1] * pow(x[ii + 3], 3) * y[ii + 2];
			p[3] += (x[ii] - x[ii + 1]) * (y[ii + 3] * (x[ii] - x[ii + 2]) * (x[ii + 2] - x[ii + 1]) + y[ii + 2] * (x[ii] - x[ii + 3]) * (x[ii + 1] - x[ii + 3])) - y[ii] * (x[ii] - x[ii + 2]) * (x[ii] - x[ii + 3]) * (x[ii + 2] - x[ii + 3]) + y[0] * (x[ii + 1] - x[ii + 2]) * (x[ii + 1] - x[ii + 3]) * (x[ii + 2] - x[ii + 3]);
		}

		for (int ii = 0; ii < p.size(); ++ii) {
			p[ii] /= (x.size() - 4);
		}

		p[0] /= detM;
		p[1] /= detM;
		p[2] /= detM;
		p[3] /= detM;

		cout << "Parameters used for fit: p0 = " << p[0]
			<< ", p1 = " << p[1]
			<< ", p2 = " << p[2]
			<< ", p3 = " << p[3] << endl;

		return 0;
	}


#elif CHOOSE_FIT_FUNCTION == 1 

	double fit_function(//T_c = flex point = p[3]
		double* y, //temperatures
		double* p //parameters
	) {
		if (p[3] < 0 ) return 1e30;
		return p[0] + p[1] * atan(p[2] * (y[0] - p[3]));
	};

	const int n_par_fit = 4;
	
	bool par_estimate(//considering fit function as y(x) = p[0] + p[1] * atan{p[2] * (x - p[3])};
		//return 0 if success, 1 if not;
		const vector<double>& x, //temperatures
		const vector<double>& y, //observables
		vector<double>& p //parameters
	) {

		if (bool_chose_at_eye || (x.size() < n_par_fit)) {//TO CHANGE THE FOLLOWING FOR "CHOOSE BY EYE" SETTING
			p[0] = 0.024;
			p[1] = 0.1 /PI;
			p[2] = 0.006;
			p[3] = 250;

			if (x.size() < n_par_fit) {
				cerr << "Not enough points for estimate." << endl;
				return 1;
			}

			return 0;
		}

		// ---- 1. p3: point x of the transition (I choose it by eye)
		p[3] = 240;//must be about the inflection point

		// ---- 2. p0: y point of the transition (I choose it by eye):
		p[0] = 0.02;//must be about the inflection point

		// ---- 3. p1: step height/PI
		p[1] = 2.0 / PI * (y[x.size() - 1] * 1.1 - p[0]);

		// ---- 4. p2: I use previous estimates and the function form:
		p[2] = tan((y[2] - p[0]) / p[1]) / (x[2] - p[3]);

		cout << "Parameters used for fit: p0 = " << p[0]
			<< ", p1 = " << p[1]
			<< ", p2 = " << p[2]
			<< ", p3 = " << p[3] << endl;
		return 0;
	}
	
#elif CHOOSE_FIT_FUNCTION == 2

	const int n_par_fit = 4;

	double fit_function(//T_c = flex point = p[3]
		double* y, //temperatures
		double* p //parameters
	) {
		return p[0] + p[1] / (1 + exp(-p[2] * (y[0] - p[3])));
	}

	bool par_estimate(//considering fit function as y(x) = p[0] + p[1]/{1 + exp[-p[2]*(x - p[3])]};
		//return 0 if success, 1 if not;
		const vector<double>& x, //temperatures
		const vector<double>& y, //observables
		vector<double>& p //parameters
	) {
		
		if (bool_chose_at_eye || (x.size() < n_par_fit)) {//TO CHANGE THE FOLLOWING FOR "CHOOSE BY EYE" SETTING
			p[0] = 1;
			p[1] = 1;
			p[2] = 1;
			p[3] = 1;

			if (x.size() < n_par_fit) {
				cerr << "Not enough points for estimate." << endl;
				return 1;
			}

			return 0;
		}
		
		int i_infl = 1;
		double max_der = -1e10, dy, dx, der, y20, y80, x20, width;

		size_t n = x.size();

		//1. p[0] = baseline: average of the first 3 values:
		p[0] = (y[0] + y[1] + y[2]) / 3.0;
		
		//2. p[1] = amplitude = (average of last 3) - baseline:
		p[1] = (y[n - 1] + y[n - 2] + y[n - 3]) / 3.0 - p[0];
		
		//3. p[3] = inflection point: maximum numerical gradient:
		for (int ii = 1; ii < n - 1; ++ii) {
			dy = y[ii + 1] - y[ii - 1];
			dx = x[ii + 1] - x[ii - 1];
			der = dy / dx;
			if (der > max_der) {
				max_der = der;
				i_infl = ii;
			}
		}
		p[3] = x[i_infl];

		//4. p[2] = steepness: I assume it to be: 4/(width between 20% and 80%)
		y20 = p[0] + 0.2 * p[1];
		y80 = p[0] + 0.8 * p[1];
		x20 = x[0];
		x80 = x[n - 1];
		for (int ii = 1; ii < n; ++ii) {
			if (((y[ii - 1] < y20) && (y[ii] > y20))) {
				x20 = x[ii];
			}
			if (((y[ii - 1] < y80) && (y[ii] > y80))) {
				x80 = x[ii]; 
				break;
			}
		}
		width = x80 - x20;
		p[2] = width > 0 ? (4.0 / width) : 1.0;

		cout << "Parameters used for fit: p0 = " << p[0]
			<< ", p1 = " << p[1]
			<< ", p2 = " << p[2]
			<< ", p3 = " << p[3] << endl;
		return 0;
	}
	
#elif CHOOSE_FIT_FUNCTION == 3

	double fit_function(//Hill function: //T_c = flex point = 1/p[2]
		double* y, //temperatures
		double* p //parameters
	) {
		return p[0] + p[1] / (1 + pow(p[2] * y[0],p[3]));
	}

	const int n_par_fit = 4;

	bool par_estimate(//considering fit function as y(x) = p[0] + p[1]/{1 + (p[2] * x)^p[3]};
		//return 0 if success, 1 if not;
		const vector<double>& x, //temperatures
		const vector<double>& y, //observables
		vector<double>& p //parameters
	) {
		
		if (bool_chose_at_eye || (x.size() < n_par_fit)) {//TO CHANGE THE FOLLOWING FOR "CHOOSE BY EYE" SETTING
			p[0] = 1;
			p[1] = 1;
			p[2] = 1;
			p[3] = 1;
			
			if (x.size() < n_par_fit) {
				cerr << "Not enough points for estimate." << endl;
				return 1;
			}

			return 0;
		}
		
		double y_half, x_half, min_diff;

		int n = x.size();

		//1. p[0] = baseline = average of the first 3 values:
		p[0] = (y[0] + y[1] + y[2]) / 3.0;

		//2. p[1] = amplitude = (average of last 3) - baseline:
		p[1] = (y[n - 1] + y[n - 2] + y[n - 3]) / 3.0 - p[0];

		//3. p[2] = central_height_point (where y = p0 + 0.5*p1)
		y_half = p[0] + 0.5 * p[1];
		x_half = x[0];
		min_diff = fabs(y[0] - y_half);
		for (int ii = 1; ii < n; ++ii) {
			if (fabs(y[ii] - y_half) < min_diff) {
				min_diff = fabs(y[ii] - y_half);
				x_half = x[ii];
			}
		}
		p[2] = 1.0 / x_half;

		//4. p[3] = Hill coefficient (stepness): (I give a default value = 2)
		p[3] = 2.0;

		cout << "Parameters used for fit: p0 = " << p[0]
			<< ", p1 = " << p[1]
			<< ", p2 = " << p[2]
			<< ", p3 = " << p[3] << endl;
		
		return 0;
		}

#elif CHOOSE_FIT_FUNCTION == 4

	const int n_par_fit = 5;

	double fit_function(//T_c = flex point = p[4]
		double* y, //temperatures
		double* p //parameters
	) {
		// p[0] = (intercept bottom)
		// p[1] = (bottom slope)
		// p[2] = (step width)/PI
		// p[3] = stepness
		// p[4] = transition temperature
		return p[0] + p[1] * y[0] + p[2] * atan(p[3] * (y[0] - p[4]));
	}

	bool par_estimate(//considering fit function as y(x) = p[0] + p[1]*x + p[2] * atan{p[3] * (x - p[4])};
		//return 0 if success, 1 if not;
		const vector<double>& x, //temperatures
		const vector<double>& y, //observables
		vector<double>& p //parameters
	) {

		if (bool_choose_at_eye || (x.size() < n_par_fit)) {//TO CHANGE THE FOLLOWING FOR "CHOOSE BY EYE" SETTING
			p[0] = 1;
			p[1] = 1;
			p[2] = 1;
			p[3] = 1;

			if (x.size() < n_par_fit) {
				cerr << "Not enough points for estimate." << endl;
				return 1;
			}

			return 0;
		}

		int n = x.size(), n_fit = 3, i_infl = 1;
		double sum_x = 0, sum_y = 0, sum_xx = 0, sum_xy = 0, max_der = -1e20;
		double denom, y_top, x_top, fondo_top, y20, x20, y80, x80, width;

		// 1. p[0] and p[1], by a linear fit on the first 3 points: (y =ca p[0] + p[1]*x)
		for (int ii = 0; ii < n_fit; ++ii) {
			sum_x += x[ii];
			sum_y += y[ii];
			sum_xx += x[ii] * x[ii];
			sum_xy += x[ii] * y[ii];
		}
		denom = n_fit * sum_xx - sum_x * sum_x;
		p[1] = (denom != 0) ? (n_fit * sum_xy - sum_x * sum_y) / denom : 0;
		p[0] = (sum_y - p[1] * sum_x) / n_fit;

		// 2. p[2] = amplitude/PI, where: amplitude = (difference between average of last points and bottom line (at maximum x))
		y_top = (y[n - 1] + y[n - 2] + y[n - 3]) / 3.0;
		x_top = (x[n - 1] + x[n - 2] + x[n - 3]) / 3.0;
		fondo_top = p[0] + p[1] * x_top;
		p[2] = (y_top - fondo_top) / PI; //amplitude/PI

		//3. p[3] = inflection point: maximum numerical gradient:
		for (int ii = 1; ii < (n - 1); ++ii) {
			double dy = y[ii + 1] - y[ii - 1];
			double dx = x[ii + 1] - x[ii - 1];
			double der = dy / dx;
			if (abs(der) > max_der) {
				max_der = abs(der);
				i_infl = ii;
			}
		}
		p[4] = x[i_infl];

		// 4. p[3] = steepness: I assume it to be: 4/(width between 20% and 80%)
		y20 = fondo_top + 0.2 * (y_top - fondo_top);
		y80 = fondo_top + 0.8 * (y_top - fondo_top);
		x20 = x[0];
		x80 = x[n - 1];
		for (int ii = 1; ii < n; ++ii) {
			if ((y[ii - 1] < y20 && y[ii] > y20)) {
				x20 = x[ii];
			}
			if ((y[ii - 1] < y80 && y[ii] > y80)) {
				x80 = x[ii]; break; 
			}
		}
		width = x80 - x20;
		p[3] = (width > 0) ? 1.0 / width : 0.01;

		cout << "Parameters used for fit: p0 = " << p[0]
			<< ", p1 = " << p[1]
			<< ", p2 = " << p[2]
			<< ", p3 = " << p[3]
			<< ", p4 = " << p43] << endl;

		return 0;
}

#elif CHOOSE_FIT_FUNCTION == 5//QUIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

double fit_function(//T_c = flex point = p[3]
	double* y, //temperatures
	double* p //parameters
) {
	return p[0] + p[1] * y[0] / (1 + exp(-p[2] * (y[0] - p[3])));
}

const int n_par_fit = 4;

bool par_estimate(//considering fit function as y(x) = p[0] + p[1]*x/{ 1 - exp[ -p[2] * (x - p[3]) ] };
	//return 0 if success, 1 if not;
	const vector<double>& x, //temperatures
	const vector<double>& y, //observables
	vector<double>& p //parameters
) {

	if (bool_choose_at_eye || (x.size() < n_par_fit)) {//TO CHANGE THE FOLLOWING FOR "CHOOSE BY EYE" SETTING
		p[0] = 1;
		p[1] = 1;
		p[2] = 1;
		p[3] = 1;

		if (x.size() < n_par_fit) {
			cerr << "Not enough points for estimate." << endl;
			return 1;
		}

		return 0;
	}

	size_t n = x.size();
	// 1. baseline: media dei primi 3 y
	p[0] = (y[0] + y[1] + y[2]) / 3.0;

	// 2. pendenza (ampiezza “scalante”) fit lineare ultimi 3
	double sum_x = 0, sum_y = 0, sum_xx = 0, sum_xy = 0;
	int n_fit = 3;
	for (int ii = n - n_fit; ii < n; ++ii) {
		sum_x += x[ii];
		sum_y += y[ii];
		sum_xx += x[ii] * x[ii];
		sum_xy += x[ii] * y[ii];
	}
	double denom = n_fit * sum_xx - sum_x * sum_x;
	p[1] = (denom != 0) ? (n_fit * sum_xy - sum_x * sum_y) / denom : 0.0;

	// 3. centro (p3): massimo gradiente numerico
	int idx = 1;
	double max_der = -1e20;
	for (int ii = 1; ii < n - 1; ++ii) {
		double dy = y[ii + 1] - y[ii - 1];
		double dx = x[ii + 1] - x[ii - 1];
		double der = dy / dx;
		if (der > max_der) {
			max_der = der;
			idx = i;
		}
	}
	p[3] = x[idx];

	// 4. ripidità (p2): 4/(larghezza tra 20% e 80%)
	double y_min = *std::min_element(y.begin(), y.end());
	double y_max = *std::max_element(y.begin(), y.end());
	double y20 = y_min + 0.2 * (y_max - y_min);
	double y80 = y_min + 0.8 * (y_max - y_min);
	double x20 = x[0], x80 = x[n - 1];
	for (int ii = 1; ii < n; ++ii) {
		if ((y[ii - 1] < y20 && y[ii] > y20)) x20 = x[ii];
		if ((y[ii - 1] < y80 && y[ii] > y80)) { x80 = x[ii]; break; }
	}
	double width = x80 - x20;
	p[2] = (width > 0) ? 4.0 / width : 1.0;

	cout << "Parameters used for fit: p0 = " << p[0]
		<< ", p1 = " << p[1]
		<< ", p2 = " << p[2]
		<< ", p3 = " << p[3] << endl;
	return 0;
}

#elif CHOOSE_FIT_FUNCTION == 6

double fit_function(//T_c = flex point = p[3]
	double* y, //temperatures
	double* p //parameters
) {
	//return pow(p[0] * y[0] / (1 + exp(-p[2] * (y[0] - p[3]))), p[1]);
	// Proteggi exp
	double expo = -p[2] * (y[0] - p[3]);
	if (expo > 700) expo = 700;
	if (expo < -700) expo = -700;
	double denom = 1.0 + exp(expo);

	// Calcola argomento della pow
	double arg = p[0] * y[0] / denom;
	if (arg <= 1e-12) arg = 1e-12; // evita base negativa/zero

	// Penalità per parametri non fisici
	if (p[0] <= 0 || p[1] <= 0 || p[2] <= 0)
		return 1e30;

	return pow(arg, p[1]);
}


const int n_par_fit = 4;

bool par_estimate(//considering fit function as y(x) = { p[0]*x / [ 1 + exp( -p[2] * (x - p[3]) ) ] }^p[1];
	//return 0 if success, 1 if not;
	const vector<double>& x, //temperatures
	const vector<double>& y, //observables
	vector<double>& p //parameters
) {

	if (bool_choose_at_eye || (x.size() < n_par_fit)) {//TO CHANGE THE FOLLOWING FOR "CHOOSE BY EYE" SETTING
		
		p[0] = 1.0 / 20000;
		p[1] = 0.79;
		p[2] = 1.0 / 10;
		p[3] = 240.0;
		
		if (x.size() < n_par_fit) {
			cerr << "Not enough points for estimate." << endl;
			return 1;
		}

		return 0;
	}

	size_t n = x.size();
	// Fit log-log agli ultimi 3 punti per scala e potenza
	double logx1 = std::log(x[n - 3]), logx2 = std::log(x[n - 2]), logx3 = std::log(x[n - 1]);
	double logy1 = std::log(y[n - 3]), logy2 = std::log(y[n - 2]), logy3 = std::log(y[n - 1]);
	double sx = logx1 + logx2 + logx3;
	double sy = logy1 + logy2 + logy3;
	double sxx = logx1 * logx1 + logx2 * logx2 + logx3 * logx3;
	double sxy = logx1 * logy1 + logx2 * logy2 + logx3 * logy3;
	double denom = 3 * sxx - sx * sx;
	p[1] = (denom != 0) ? (3 * sxy - sx * sy) / denom : 1.0; // esponente potenza
	double loga = (3 * sy - p[1] * sx) / 3.0;
	p[0] = std::exp(loga); // scala

	// Centro della transizione (p3): massimo gradiente numerico
	int idx = 1;
	double max_der = -1e20;
	for (int ii = 1; ii < n - 1; ++ii) {
		double dy = y[ii + 1] - y[ii - 1];
		double dx = x[ii + 1] - x[ii - 1];
		double der = dy / dx;
		if (der > max_der) {
			max_der = der;
			idx = ii;
		}
	}
	p[3] = x[idx];

	// Ripidità (p2): 4/(larghezza tra 20% e 80%)
	double y_min = *std::min_element(y.begin(), y.end());
	double y_max = *std::max_element(y.begin(), y.end());
	double y20 = y_min + 0.2 * (y_max - y_min);
	double y80 = y_min + 0.8 * (y_max - y_min);
	double x20 = x[0], x80 = x[n - 1];
	for (int ii = 1; ii < n; ++ii) {
		if ((y[ii - 1] < y20 && y[ii] > y20)) x20 = x[ii];
		if ((y[ii - 1] < y80 && y[ii] > y80)) { x80 = x[ii]; break; }
	}
	double width = x80 - x20;
	p[2] = (width > 0) ? 4.0 / width : 1.0;

	cout << "Parameters used for fit: p0 = " << p[0]
		<< ", p1 = " << p[1]
		<< ", p2 = " << p[2]
		<< ", p3 = " << p[3] << endl;

	return 0;
}

#endif


//-----------------------------------------------------------------
//ROOT MACRO TO DO FIT AND GRAPH:

void fit_plot_points_errors(
	const vector<double>& x,
	const vector<double>& y,
	const vector<double>& y_err,
	const string name_image,
	const string title,
	const string y_name,
	const double pos_title,
	const double pos_y,
	const double heigh_y,
	const int n_par_fit
);

void fit_plot_inverted_points_errors(
	const vector<double>& x,
	const vector<double>& y,
	const vector<double>& y_err,
	const string name_image,
	const string title,
	const string y_name,
	const double pos_title,
	const double pos_y,
	const double heigh_y
);

void silly_plot(
	const vector<double>& x,
	const vector<double>& y,
	const vector<double>& y_err,
	const string name_image,
	const string title,
	const string y_name,
	const double pos_title,
	const double pos_y,
	const double heigh_y,
	const int n_par_fit
);

double chi2_reduced_estimate(
	const vector<double>& x,
	const vector<double>& y,
	const vector<double>& y_err,
	vector<double>& par
);


//-----------------------------------------------------------------
//MAIN:

int main() {

	int skipLines = 1; //= number of lines to skip while reading input file;
	double temp_value, mod_value, mod_err_value, re_value, re_err_value, im_value, im_err_value;
	vector <double> temp, mod, mod_err, re, re_err, im, im_err;
	size_t pos;
	string line, name_tmp;
	string input_directory = "19_05_2025/polyff_results/";
	string output_directory = "results/";
	string name_input_file;
	string name_image_mod;
	string name_image_re;
	string name_image_im;
	string title_mod;
	string title_re;
	string title_im;
	string y_name_mod;
	string y_name_re;
	string y_name_im;
	double pos_ymod;
	double pos_yre;
	double pos_yim;
	double height_mod;
	double height_re;
	double height_im;
	double pos_title_mod;
	double pos_title_re;
	double pos_title_im;

	if (tipology == "gauge") {
		name_input_file = "1500.0_poly_results.txt";
		pos = name_input_file.find_last_of(".");
		if (pos != string::npos) {
			name_tmp = name_input_file.substr(0, pos); //I remove extension using substr
		}
		name_input_file = input_directory + name_input_file;
		name_image_mod = output_directory + "FIT_modPvsT_" + name_tmp + ".png";
		name_image_re = output_directory + "FIT_rePvsT_" + name_tmp + ".png";
		name_image_im = output_directory + "FIT_imPvsT_" + name_tmp + ".png";
		title_mod = "#LT|PP^{+}|#GT vs temperature:";
		title_re = "#LTRe{P}#GT vs temperature :";
		title_im = "#LTIm{P}#GT vs temperature:";
		y_name_mod = "#LT|PP^{+}|#GT";
		y_name_re = "#LTRe{P}#GT";
		y_name_im = "#LTIm{P}#GT";
		pos_ymod = 0.02;
		pos_yre = 0.03;
		pos_yim = 0.03;
		height_mod = 0.4;
		height_re = 0.45;
		height_im = 0.45;
		pos_title_mod = 0.3;
		pos_title_re = 0.3;
		pos_title_im = 0.3;
	}
	else if (tipology == "fermion") {
		name_input_file = "1500.0_ff_results.txt";
		pos = name_input_file.find_last_of(".");
		if (pos != string::npos) {
			name_tmp = name_input_file.substr(0, pos); //I remove extension using substr
		}
		name_input_file = input_directory + name_input_file;
		name_image_mod = output_directory + "FIT_modffvsT_" + name_tmp + ".png";
		name_image_re = output_directory + "FIT_reffvsT_" + name_tmp + ".png";
		name_image_im = output_directory + "FIT_imffvsT_" + name_tmp + ".png";
		title_mod = "#LT|(#bar{#psi}#psi)(#bar{#psi}#psi)^{+}|#GT vs temperature:";
		title_re = "#LTRe{#bar{#psi}#psi}#GT vs temperature :";
		title_im = "#LTIm{#bar{#psi}#psi}#GT vs temperature:";
		y_name_mod = "#LT|(#bar{#psi}#psi)(#bar{#psi}#psi)^{+}|#GT";
		y_name_re = "#LTRe{#bar{#psi}#psi}#GT";
		y_name_im = "#LTIm{#bar{#psi}#psi}#GT";
		pos_ymod = 0.03;
		pos_yre = 0.03;
		pos_yim = 0.04;
		height_mod = 0.4;
		height_re = 0.45;
		height_im = 0.45;
		pos_title_mod = 0.3;
		pos_title_re = 0.3;
		pos_title_im = 0.3;
	}
	else {
		cerr << "Invalid tipology: you must write gauge or fermion";
		return 1;
	}

	ifstream input_file;
	input_file.open(name_input_file);
	if (!input_file) {
		cerr << "Error opening input file." << endl;
		return 1;
	}

	for (int ii = 0; ii < skipLines; ++ii) {
		if (!getline(input_file, line)) {
			cerr << "Error: there are less than " << skipLines << " lines in the input file." << endl;
			return 1;
		}
	}

	while (getline(input_file, line)) {
		istringstream iss(line);
		if (iss >> temp_value >> mod_value >> mod_err_value >> re_value >> re_err_value >> im_value >> im_err_value) {
			temp.push_back(temp_value);
			mod.push_back(mod_value);
			mod_err.push_back(mod_err_value);
			re.push_back(re_value);
			re_err.push_back(re_err_value);
			im.push_back(im_value);
			im_err.push_back(im_err_value);
		}
		else {
			cerr << "Poorly formatted line: " << line << endl;
		}
	}

	input_file.close();

	if (bool_enlarge) {
		for (int ii = 0; ii < temp.size(); ++ii) {
			mod[ii] *= 100;
			mod_err[ii] *= 100;
			re[ii] *= 100;
			re_err[ii] *= 100;
			im[ii] *= 100;
			im_err[ii] *= 100;
		}
		y_name_mod = y_name_mod + "*100";
		y_name_re = y_name_re + "*100";
		y_name_im = y_name_im + "*100";
	}


	//GRAPHIC REPRESENTATION:
	silly_plot(temp, mod, mod_err, "results/silly_plot_1.png", title_mod, y_name_mod, pos_title_mod, pos_ymod, height_mod, n_par_fit);
	silly_plot(temp, re, re_err, "results/silly_plot_2.png", title_re, y_name_re, pos_title_re, pos_yre, height_re, n_par_fit);
	silly_plot(temp, im, im_err, "results/silly_plot_3.png", title_im, y_name_im, pos_title_im, pos_yim, height_im, n_par_fit);


	//GRAPHIC REPRESENTATION:
	fit_plot_points_errors(temp, mod, mod_err, name_image_mod, title_mod, y_name_mod, pos_title_mod, pos_ymod, height_mod, n_par_fit);
	fit_plot_points_errors(temp, re, re_err, name_image_re, title_re, y_name_re, pos_title_re, pos_yre, height_re, n_par_fit);
	fit_plot_points_errors(temp, im, im_err, name_image_im, title_im, y_name_im, pos_title_im, pos_yim, height_im, n_par_fit);

	//GRAPHIC REPRESENTATION:
	/*
	fit_plot_inverted_points_errors(temp, mod, mod_err, "results/inverted_MOD.png", title_mod, y_name_mod, pos_title_mod, pos_ymod, height_mod);
	fit_plot_inverted_points_errors(temp, re, re_err, "results/inverted_RE.png", title_re, y_name_re, pos_title_re, pos_yre, height_re);
	fit_plot_inverted_points_errors(temp, im, im_err, "results/inverted_IM.png", title_im, y_name_im, pos_title_im, pos_yim, height_im);
	*/

	return 0;
}

//-----------------------------------------------------------------
//SOME FUNCTION DEFINITIONS:


double chi2_reduced_estimate(
	const vector<double>& x,
	const vector<double>& y,
	const vector<double>& y_err,
	vector<double>& par
) {
	double chi2r = 0, res;
	for (int ii = 0; ii < x.size(); ++ii) {
		double xx[1] = { x[ii] };
		res = (y[ii] - fit_function(xx, par.data()))/y_err[ii];
		res *= res;
		chi2r += res;
		//cout << "x["<< ii<< "]:" << x[ii] << endl;
		//cout << "y[" << ii << "]:" << y[ii] << endl;
		//cout << "y_err[" << ii << "]:" << y_err[ii] << endl;
	}
	cout << "CHI2 = " << chi2r << endl;
	chi2r /= (x.size() - par.size());
	cout << "CHI2 reduced = " << chi2r << endl;
	return chi2r;
}

//-----------------------------------------------------------------
//ROOT MACRO TO GRAPH:

void fit_plot_points_errors(
	const vector<double>& x,
	const vector<double>& y,
	const vector<double>& y_err,
	const string name_image,
	const string title,
	const string y_name,
	const double pos_title,
	const double pos_y,
	const double heigh_y,
	const int n_par_fit
) {

	if (x.size() != y.size() || y.size() != y_err.size()) {
		cerr << "Error: mismatched vector sizes in plot_points_errors()." << endl;
		return;
	}

	cout << y_err[4] << endl;

	vector <double> par(n_par_fit,0);

	//CANVAS creation: (to draw the graph)
	TCanvas* canvas = new TCanvas("canvas", "Canvas for Drawing Points", 900, 600);
	canvas->SetGrid();//to set grid

	// 1. Graph with only error bars:
	TGraphErrors* g_errors = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, y_err.data());

	//Using std::max_element: (useful to create my histogram object)
	//	-> std::max_element takes two iterators that define the range of the array to operate on. 
	//	-> returns an iterator that points to the maximum element found.
	//The usage of std::min_element is analogous but with for minimum element.
	auto min_x = *min_element(x.begin(), x.end());
	auto max_x = *max_element(x.begin(), x.end());
	auto min_y = *min_element(y.begin(), y.end());
	auto max_y = *max_element(y.begin(), y.end());

	g_errors->SetLineColor(kBlack);//color of error bars
	g_errors->SetMarkerColor(kBlack);
	g_errors->SetMarkerStyle(20);
	g_errors->SetTitle("");
	g_errors->GetXaxis()->SetLimits(min_x, max_x);
	g_errors->GetYaxis()->SetRangeUser(min_y - 0.01 * fabs(min_y), max_y + 0.01 * fabs(max_y));
	g_errors->Draw("AP");

	// 2. Fit of points:

	gStyle->SetOptFit(1111);

	double fit_min = 0;
	double fit_max = max_x + 200;

	TF1* p_plot = new TF1("fit_function", fit_function, fit_min, fit_max, n_par_fit);

	//p_plot->SetParLimits(0, 1e-10, 10);    // par[0]: solo positivi
	//p_plot->SetParLimits(1, 0.1, 5);       // par[1]: esponente ragionevole
	//p_plot->SetParLimits(2, 1e-10, 10);    // par[2]: ripidità positiva
	//p_plot->SetParLimits(3, min_x, max_x);   // par[3]: centro nella regione dei dati

	par_estimate(x, y, par);

	for (int ii = 0; ii < par.size(); ++ii) {
		p_plot->SetParameter(ii, par[ii]); //setting the ii-th parameter of the function
	}

	p_plot->GetXaxis()->SetRangeUser(min_x, max_x);
	p_plot->GetYaxis()->SetRangeUser(0, g_errors->GetMaximum() * 1.01);
	
	
	//FIT:
	gStyle->SetOptFit(0);// (1111);
	g_errors->Fit(p_plot, "SW"); // S = return result, W = use weights (i.e., y errors)
	
	TVirtualFitter* fit = TVirtualFitter::GetFitter(); //to get fit results;
	//ATTENTA!! Spesso, ROOT stampa "Chi2" come la media dei residui quadratici(RMS²), NON come la vera somma pesata!
	//Questo è uno di questi casi, dei valori ci si può fidare se il chi2 stimato con l'apposita funzione nel seguito è ragionevole.
	
	for (int ii = 0; ii < par.size(); ++ii) {
		par[ii] = p_plot->GetParameter(ii);
	}

	chi2_reduced_estimate(
		x,
		y,
		y_err,
		par
	);
	
	double cov[10][10];

	cout << endl;
	cout << "COV_MATRIX:" << endl;
	for (int ii = 0; ii < n_par_fit; ++ii)
	{
		for (int jj = 0; jj < n_par_fit; ++jj)
		{
			cov[ii][jj] = fit->GetCovarianceMatrixElement(ii, jj);
			cout << cov[ii][jj] << "\t";
		}
	
		cout << endl;

	}

	cout << endl;

	cout << endl;
	cout << "CORR_MATRIX:" << endl;
	for (int ii = 0; ii < n_par_fit; ++ii)
	{
		for (int jj = 0; jj < n_par_fit; ++jj)
		{
			cout << cov[ii][jj] /sqrt(cov[ii][ii] * cov[jj][jj]) << "\t";
		}

		cout << endl;

	}

	cout << endl;



	p_plot ->SetLineColor(kRed);
	p_plot->Draw("SAME");

	gPad->Update();

	//LATEX: I add LaTeX titles and axis labels:
	TLatex latex;
	latex.SetNDC(); //sets the use of Normalized Device Coordinates (NDC).
	latex.SetTextSize(0.05); //changes text size for title
	latex.DrawLatex(pos_title, 0.92, title.c_str());
	latex.SetTextSize(0.04); //changes text size for axis labels
	latex.DrawLatex(0.45, 0.03, "T[MeV]");
	latex.SetTextAngle(90);
	latex.DrawLatex(pos_y, heigh_y, (y_name).c_str());


	//SAVE: I save the canvas as an image
	canvas->SaveAs(name_image.c_str());

	//to save also in vectorial pdf form:
	canvas->SaveAs((name_image.substr(0, name_image.find_last_of(".")) + ".pdf").c_str());

	//DELETE:
	delete p_plot;
	delete g_errors;
	delete canvas;
}


void fit_plot_inverted_points_errors(//CAMBIA NOME AGLI ASSI
	const vector<double>& x,
	const vector<double>& y,
	const vector<double>& y_err,
	const string name_image,
	const string title,
	const string y_name,
	const double pos_title,
	const double pos_y,
	const double heigh_y
) {

	if (x.size() != y.size() || y.size() != y_err.size()) {
		cerr << "Error: mismatched vector sizes in plot_points_errors()." << endl;
		return;
	}

	vector <double> par(4, 0);

	//CANVAS creation: (to draw the graph)
	TCanvas* canvas = new TCanvas("canvas", "Canvas for Drawing Points", 900, 600);
	canvas->SetGrid();//to set grid

	// 1. Graph with only error bars:
	TGraphErrors* g_errors = new TGraphErrors(x.size(), y.data(), x.data(), y_err.data(), nullptr);

	//Using std::max_element: (useful to create my histogram object)
	//	-> std::max_element takes two iterators that define the range of the array to operate on. 
	//	-> returns an iterator that points to the maximum element found.
	//The usage of std::min_element is analogous but with for minimum element.
	auto min_x = *min_element(x.begin(), x.end());
	auto max_x = *max_element(x.begin(), x.end());
	auto min_y = *min_element(y.begin(), y.end());
	auto max_y = *max_element(y.begin(), y.end());

	g_errors->SetLineColor(kBlack);//color of error bars
	g_errors->SetMarkerColor(kBlack);
	g_errors->SetMarkerStyle(20);
	g_errors->SetTitle("");
	g_errors->GetXaxis()->SetRangeUser(min_y - 0.01 * fabs(min_y), max_y + 0.01 * fabs(max_y));
	g_errors->GetYaxis()->SetLimits(min_x, max_x); 
	g_errors->Draw("AP");

	// 2. Fit of points:

	gStyle->SetOptFit(1111);

	double fit_min = min_x - 20;
	double fit_max = max_x + 20;

	TF1* p_plot = new TF1("fit_function", fit_function, fit_min, fit_max, n_par_fit); //4 is the number of parameters

	par_estimate(y, x, par);

	for (int ii = 0; ii < par.size(); ++ii) {
		p_plot->SetParameter(ii, par[ii]); //setting the ii-th parameter of the function
	}

	p_plot->GetXaxis()->SetRangeUser(0, g_errors->GetMaximum() * 1.01);
	p_plot->GetYaxis()->SetRangeUser(min_x, max_x);

	//FIT:
	gStyle->SetOptFit(0);// (1111);
	g_errors->Fit(p_plot, "SW"); // S = return result, W = use weights (i.e., y errors)
	
	TVirtualFitter* fit = TVirtualFitter::GetFitter(); //to get fit results;

	TVirtualFitter* fitter = TVirtualFitter::GetFitter();
	double cov[10][10]; // scegli una dimensione sufficientemente grande

	for (int ii = 0; ii < n_par_fit; ++ii)
	{
		for (int jj = 0; jj < n_par_fit; ++jj)
		{
			cov[ii][jj] = fitter->GetCovarianceMatrixElement(ii, jj);
		}

	}
	p_plot->SetLineColor(kRed);
	p_plot->Draw("SAME");

	gPad->Update();

	//LATEX: I add LaTeX titles and axis labels:
	TLatex latex;
	latex.SetNDC(); //sets the use of Normalized Device Coordinates (NDC).
	latex.SetTextSize(0.05); //changes text size for title
	latex.DrawLatex(pos_title, 0.92, title.c_str());
	latex.SetTextSize(0.04); //changes text size for axis labels
	latex.DrawLatex(0.45, 0.03, "T[MeV]");
	latex.SetTextAngle(90);
	latex.DrawLatex(pos_y, heigh_y, y_name.c_str());


	//SAVE: I save the canvas as an image
	canvas->SaveAs(name_image.c_str());

	//to save also in vectorial pdf form:
	canvas->SaveAs((name_image.substr(0, name_image.find_last_of(".")) + ".pdf").c_str());

	//DELETE:
	delete p_plot;
	delete g_errors;
	delete canvas;
}

void silly_plot(
	const vector<double>& x,
	const vector<double>& y,
	const vector<double>& y_err,
	const string name_image,
	const string title,
	const string y_name,
	const double pos_title,
	const double pos_y,
	const double heigh_y,
	const int n_par_fit
) {
	if (x.size() != y.size() || y.size() != y_err.size()) {
		cerr << "Error: mismatched vector sizes in plot_points_errors()." << endl;
		return;
	}

	vector <double> par(n_par_fit, 0);
	vector <double> y_attempt;

	//CANVAS creation: (to draw the graph)
	TCanvas* canvas = new TCanvas("canvas", "Canvas for Drawing Points", 900, 600);
	canvas->SetGrid();//to set grid

	// 1. Graph with only error bars:
	TGraphErrors* g_errors = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, y_err.data());

	//Using std::max_element: (useful to create my histogram object)
	//	-> std::max_element takes two iterators that define the range of the array to operate on. 
	//	-> returns an iterator that points to the maximum element found.
	//The usage of std::min_element is analogous but with for minimum element.
	auto min_x = *min_element(x.begin(), x.end());
	auto max_x = *max_element(x.begin(), x.end());
	auto min_y = *min_element(y.begin(), y.end());
	auto max_y = *max_element(y.begin(), y.end());

	g_errors->SetLineColor(kBlack);//color of error bars
	g_errors->SetMarkerColor(kBlack);
	g_errors->SetMarkerStyle(20);
	g_errors->SetTitle("");
	g_errors->GetXaxis()->SetLimits(min_x, max_x);
	g_errors->GetYaxis()->SetRangeUser(min_y - 0.01 * fabs(min_y), max_y + 0.01 * fabs(max_y));
	g_errors->Draw("AP");

	// 2. function plot:

	TF1* f1 = new TF1("f", fit_function, min_x, max_x, n_par_fit);
	
	par_estimate(x, y, par);
	
	for (int ii = 0; ii < par.size(); ii++) {
		f1->SetParameter(ii, par[ii]); //setting the ii-th parameter of the function
	}

	f1->SetLineColor(kRed);
	f1->Draw("SAME");

	gPad->Update();

	//LATEX: I add LaTeX titles and axis labels:
	TLatex latex;
	latex.SetNDC(); //sets the use of Normalized Device Coordinates (NDC).
	latex.SetTextSize(0.05); //changes text size for title
	latex.DrawLatex(pos_title, 0.92, title.c_str());
	latex.SetTextSize(0.04); //changes text size for axis labels
	latex.DrawLatex(0.45, 0.03, "T[MeV]");
	latex.SetTextAngle(90);
	latex.DrawLatex(pos_y, heigh_y, y_name.c_str());


	//SAVE: I save the canvas as an image
	canvas->SaveAs(name_image.c_str());

	//to save also in vectorial pdf form:
	canvas->SaveAs((name_image.substr(0, name_image.find_last_of(".")) + ".pdf").c_str());

	//DELETE:
	delete f1;
	delete g_errors;
	delete canvas;
}