/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****                  simple_computation.cpp:                 ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

/*
OUTPUT:
#beta    ml[MeV]         ml_adimensional
3.8     117.317 0.0388878
3.827   116.034 0.0362886
3.84063 115.478 0.0350957
3.881   114.14  0.0319629
3.908   113.46  0.0301512
3.935   112.904 0.0285228
3.962   112.415 0.0270421
3.989   111.941 0.0256782
4.016   111.421 0.0244032
4.043   110.794 0.0231922
4.07    110.005 0.0220253
4.097   109.001 0.020886
4.124   107.729 0.0197607
4.151   106.159 0.0186421
4.178   104.261 0.0175257
4.2     102.466 0.0166178
*/


#include "../library.h"
#include "../root_include.h"

const double hbar_c = 197.3269804; //MeV * fm

//-----------------------------------------------------------------
//MAIN:

int main() {

	int Nt = 8;
	
	vector<double> aml = { 0.081869, 0.075858, 0.073104, 0.065886, 0.061722, 0.057988, 0.054603, 0.051498, 0.048612, 0.045891, 0.043293, 0.040783, 0.038333, 0.035928, 0.033558, 0.031653 }; //a*m_l //BE CAREFUL TO CHOOSE IT WELL;
		
	vector<double> beta = { 3.80000, 3.82700, 3.84063, 3.88100, 3.90800, 3.93500, 3.96200, 3.98900, 4.01600, 4.04300, 4.07000, 4.09700, 4.12400, 4.15100, 4.17800, 4.20000 }; //BE CAREFUL TO CHOOSE IT WELL;

	vector<double> afm = { 0.1377037468, 0.1290036669, 0.1249189781, 0.1139049357, 0.1073456174, 0.1013484202, 0.0958468116, 0.0907797794, 0.0860918322, 0.0817329985, 0.0776588275, 0.0738303888, 0.0702142723, 0.0667825883, 0.0635129674, 0.0609567905 }; //a[fm] //BE CAREFUL TO CHOOSE IT WELL;

	double ml_value;//[in MeV] ml[fm^{-1}] = aml/afm -> ml[fm^{-1}] = hbar_c * aml/afm; 
	double ml_adimensional; //beta = Nt *a -> a = Nt / beta ->ml_adimensional * a = aml->ml_adimensional = aml / a = aml * beta / ((double)Nt)

	cout << "#beta \t ml[MeV] \t ml_adimensional" << endl;

	for (int ii = 0; ii < (beta.size()); ii++) {
		ml_value = hbar_c * aml[ii] / afm[ii];
		ml_adimensional = aml[ii] * beta[ii] / ((double)Nt);
		cout << beta[ii] << "\t" << ml_value << "\t" << ml_adimensional << endl;
	}

	return 0;
}
