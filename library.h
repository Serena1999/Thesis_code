﻿/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****   library.h: LIBRARY THAT CONTAINS HEADERS, DEFINITIONS  ****
****     AND FUNCTIONS FOR CODES FOR THE COURSE "NUMERICAL    ****
****    METHODS FOR PHYSICISTS" (author = Serena Bruzzesi)    ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/


using namespace std;

#include "CONFIG.h"

//-----------------------------------------------------------------
//DEFINITIONS:

# define PI acos(-1) 

struct sample_gen {

//-----------------------------------------------------------------
//RANDOM NUMBER GENERATOR:

#if random_choice_par == 1
	mt19937_64 rng;

	/*mt19937_64 is a Mersenne Twister 19937 generator (64 bit).
	A pseudo-random number generator engine that produces
	unsigned integer numbers in the closed interval [0,2^w - 1],
	where w = Word size: Number of bits of each word in the
	state sequence, that for mt19937_64 is w = 64.

	References used for the choice of pseudo-random generator:
		- https://cplusplus.com/reference/random/
		- https://cplusplus.com/reference/random/mersenne_twister_engine/
		- https://cplusplus.com/reference/random/mt19937_64/
		- https://www.delftstack.com/it/howto/cpp/cpp-timing/
		- https://it.wikipedia.org/wiki/Mersenne_Twister
		- https://en.wikipedia.org/wiki/Mersenne_Twister
	*/

#elif random_choice_par == 2
	pcg64 rng;

	/*
	Permuted congruential generator (pcg) (64 bit)
	("codesFromExternalSources/pcg-cpp-0.98")
	[2014 M.E. O'Neill / pcg-random.org]
	(https://github.com/imneme/pcg-cpp/tree/master)
	*/

#else
	ranlux48 rng;

	/*
		48-bit RANLUX generator by Martin Lüscher and Fred James (1994)
		(a subtract-with-carry pseudo-random generator)
		(it is a default generator in C++)
		(https://cplusplus.com/reference/random/ranlux48/)
	*/
#endif
	//int current_seed = 5489;//to remember the default seed of Mersenne;

	//Function to initialize the seed:
	void init(int seed) {
		if (seed != -1) {//if seed = -1, it means that we DON'T want to modify the seed
			rng.seed((unsigned)seed);//we assign the selected seed in file4Settings
		}
	}

	//Function to extract a random value in the default interval: rng()

	//Function to extract a value according to a given distribution ("distribution"): 
	template <class D> typename D::result_type extract(D& distribution) {
		return distribution(rng);
	}
	//https://en.cppreference.com/w/cpp/named_req/RandomNumberDistribution
	//Distribution_type = D (which is a generic type, assigned whitin main functions)

	//Function to extract N random values in the default interval:
	void extract(int N, uint64_t* extractions) {
		for (int i = 0; i < N; i++) {
			extractions[i] = rng();
		}
	}

	//Function to extract N values according to a given distribution ("distribution"):
	template <class D> void extract(D& distribution, int N, typename D::result_type* extractions) {
		for (int i = 0; i < N; i++) {
			extractions[i] = distribution(rng);
		}
	}

	//Functions to extract N random values in the default interval and put them in files:
	void extract_to_file(int N, int seed) {
		/*
			seed = random number generator seed. Particular cases:
				- if seed = -1: we use the default seed of the generator
			I write on these files:
				- "dataRNG.txt" contains the generated random numbers
				- "dataRNGprop.txt" contains the properties used in random number generator
		*/
		string name_file4data = "dataRNG.txt";
		string name_file4properties = "dataRNGprop.txt";
		////DATA WRITE AND READ FILES:
		ofstream file4data; //declaration of output file
		ofstream file4properties; //declaration of output file properties

		file4data.open(name_file4data);
		if (!file4data) {
			cout << "Error opening data file" << endl;
			return;
		}

		file4properties.open(name_file4properties);
		if (!file4properties) {
			cout << "Error opening properties file" << endl;
			file4data.close();
			return;
		}

		if (seed != -1) {//if seed = -1, it means that we want to use the default seed
			rng.seed((unsigned)seed);//we assign the selected seed in file4Settings
			file4properties << "#random_choice_par \t seed \t number of generated values:" << endl;
			file4properties << random_choice_par << "\t" << seed << "\t" << N;
			//!!
			file4data << "#Generator n." << random_choice_par << " - " << N << " extractions, seed = " << seed << endl;
		}
		else {
			file4properties << "#random_choice_par \t seed \t number of generated values:" << endl;
			file4properties << random_choice_par << "\t" << "default\t" << N;
			//!!
			file4data << "#Generator n." << random_choice_par << " - " << N << " extractions, seed = default" << endl;
		}

		//PSEUDO-RANDOM NUMBERS GENERATION:
		for (int i = 0; i < N; i++) {
			file4data << rng() << endl;
		}

		file4data.close();
		file4properties.close();
	}

	void extract_to_file(int N, int seed, char* name_file4data) {
		/*
			name_file4data = name of the file where I write the data
		*/
		string name_file4properties = "dataRNGprop.txt";
		////DATA WRITE AND READ FILES:
		ofstream file4data; //declaration of output file
		ofstream file4properties; //declaration of output file properties

		file4data.open(name_file4data);
		if (!file4data) {
			cout << "Error opening data file" << endl;
			return;
		}

		file4properties.open(name_file4properties);
		if (!file4properties) {
			cout << "Error opening properties file" << endl;
			file4data.close();
			return;
		}

		if (seed != -1) {//if seed = -1, it means that we want to use the default seed
			rng.seed((unsigned)seed);//we assign the selected seed in file4Settings
			file4properties << "#random_choice_par \t seed \t number of generated values:" << endl;
			file4properties << random_choice_par << "\t" << seed << "\t" << N;
			//!!
			file4data << "#Generator n." << random_choice_par << " - " << N << " extractions, seed = " << seed << endl;
		}
		else {
			file4properties << "#random_choice_par \t seed \t number of generated values:" << endl;
			file4properties << random_choice_par << "\t" << "default\t" << N;
			//!!
			file4data << "#Generator n." << random_choice_par << " - " << N << " extractions, seed = default" << endl;
		}

		for (int i = 0; i < N; i++) {
			file4data << rng() << endl;
		}

		file4data.close();
		file4properties.close();
	}

	void extract_to_file(int N, int seed, char* name_file4data, char* name_file4properties) {
		/*
			name_file4properties = name of the file where I write the properties of our random number extraction
		*/

		////DATA WRITE AND READ FILES:
		ofstream file4data; //declaration of output file
		ofstream file4properties; //declaration of output file properties

		file4data.open(name_file4data);
		if (!file4data) {
			cout << "Error opening data file" << endl;
			return;
		}

		file4properties.open(name_file4properties);
		if (!file4properties) {
			cout << "Error opening properties file" << endl;
			file4data.close();
			return;
		}

		if (seed != -1) {//if seed = -1, it means that we want to use the default seed
			rng.seed((unsigned)seed);//we assign the selected seed in file4Settings
			file4properties << "#random_choice_par \t seed \t number of generated values:" << endl;
			file4properties << random_choice_par << "\t" << seed << "\t" << N;
			//!!
			file4data << "#Generator n." << random_choice_par << " - " << N << " extractions, seed = " << seed << endl;
		}
		else {
			file4properties << "#random_choice_par \t seed \t number of generated values:" << endl;
			file4properties << random_choice_par << "\t" << "default\t" << N;
			//!!
			file4data << "#Generator n." << random_choice_par << " - " << N << " extractions, seed = default" << endl;
		}

		for (int i = 0; i < N; i++) {
			file4data << rng() << endl;
		}

		file4data.close();
		file4properties.close();
	}

	//Functions to extract N values according to a given distribution ("distribution"):
	template <class D> void extract_to_file(D& distribution, int N, int seed) {
		string name_file4data = "dataRNG_differentDistribution.txt";
		string name_file4properties = "dataRNGprop_differentDistribution.txt";
		////DATA WRITE AND READ FILES:
		ofstream file4data; //declaration of output file
		ofstream file4properties; //declaration of output file properties

		file4data.open(name_file4data);
		if (!file4data) {
			cout << "Error opening file" << endl;
			return;
		}

		file4properties.open(name_file4properties);
		if (!file4properties) {
			cout << "Error opening file" << endl;
			file4data.close();
			return;
		}

		if (seed != -1) {//if seed = -1, it means that we want to use the default seed
			rng.seed((unsigned)seed);//we assign the selected seed in file4Settings
			file4properties << "#random_choice_par \t seed \t number of generated values:" << endl;
			file4properties << random_choice_par << "\t" << seed << "\t" << N;
			//!!
			file4data << "#Generator n." << random_choice_par << " - " << N << " extractions, seed = " << seed << endl;
		}
		else {
			file4properties << "#random_choice_par \t seed \t number of generated values:" << endl;
			file4properties << random_choice_par << "\t" << "default\t" << N;
			//!!
			file4data << "#Generator n." << random_choice_par << " - " << N << " extractions, seed = default" << endl;
		}

		if (random_choice_par == 3) {
			file4data << scientific << setprecision(14);
			/*
			From https://cplusplus.com/reference/random/ranlux48_base/,
			we know that wordsize = 48. So the generated random number
			contains information in the form of only the entire sequence
			of digits.
			So, I put it in scientific format and with 15 significant digits
			to keep all the information contained in the numbers
			generated by ranlux48 (which have 48 bits). In fact, the double is
			made up of 64 bits and the information contained
			is not all the digits in full, but:
				-> the sign is contained in one bit;
				-> the exponent is contained in 11 bits;
				-> the mantissa is contained in 52 bits.
			Therefore, the mantissa reaches a maximum of 2^52 = (approximately) 4.5e15.
			So, being 48 digits less than 52, we can represent all of them in the mantissa.
			So, the maximum number generated by our procedure is converted into a double
			with 48 occupied digits in mantissa and so: the maximum mantissa that can
			be generated by this generator is 2^48 = (approximately) 3e14. Consequently,
			the maximum number of significant digits generated by ranlux48 in the mantissa
			is 15: in exponential form, the first of this digit is before the decimal
			point and therefore setprecision(14) must be placed above it.
		*/
		}
		else {
			file4data << scientific << setprecision(15);

			/*
				I put it in scientific format and with 15 significant digits
				to keep all the information contained in the numbers
				generated by pcg64/mersenne (which have 64 bits). In fact, the double is
				also made up of 64 bits, but the information contained
				is not all the digits in full, but:
					-> the sign is contained in a bit;
					-> the exponent is contained in 11 bits;
					-> the mantissa is contained in 52 bits.
				So, the mantissa reaches a maximum of 2^52 = (approximately) 4.5e15.
				Consequently, the maximum number of significant digits in the
				mantissa is 16. However, if we exceed the number of digits, we
				would not make conceptual errors because, when we recall this
				file, the data will be interpreted in double and therefore the
				less significant digits are not considered if the more
				significant ones have already filled all the space available for the
				double.
			*/
		}

		for (int i = 0; i < N; i++) {
			file4data << distribution(rng) << endl;
		}

		file4data.close();
		file4properties.close();
	}

	template <class D> void extract_to_file(D& distribution, int N, int seed, char* name_file4data) {

		string name_file4properties = "dataRNGprop_differentDistribution.txt";
		////DATA WRITE AND READ FILES:
		ofstream file4data; //declaration of output file
		ofstream file4properties; //declaration of output file properties

		file4data.open(name_file4data);
		if (!file4data) {
			cout << "Error opening data file" << endl;
			return;
		}

		file4properties.open(name_file4properties);
		if (!file4properties) {
			cout << "Error opening properties file" << endl;
			file4data.close();
			return;
		}
		if (seed != -1) {//if seed = -1, it means that we want to use the default seed
			rng.seed((unsigned)seed);//we assign the selected seed in file4Settings
			file4properties << "#random_choice_par \t seed \t number of generated values:" << endl;
			file4properties << random_choice_par << "\t" << seed << "\t" << N;
			//!!
			file4data << "#Generator n." << random_choice_par << " - " << N << " extractions, seed = " << seed << endl;
		}
		else {
			file4properties << "#random_choice_par \t seed \t number of generated values:" << endl;
			file4properties << random_choice_par << "\t" << "default\t" << N;
			//!!
			file4data << "#Generator n." << random_choice_par << " - " << N << " extractions, seed = default" << endl;
		}

		if (random_choice_par == 3) {
			file4data << scientific << setprecision(14);
		}
		else {
			file4data << scientific << setprecision(15);
		}

		for (int i = 0; i < N; i++) {
			file4data << distribution(rng) << endl;
		}

		file4data.close();
		file4properties.close();
	}

	template <class D> void extract_to_file(D& distribution, int N, int seed, char* name_file4data, char* name_file4properties) {
		////DATA WRITE AND READ FILES:
		ofstream file4data; //declaration of output file
		ofstream file4properties; //declaration of output file properties

		file4data.open(name_file4data);
		if (!file4data) {
			cout << "Error opening file" << endl;
			return;
		}

		file4properties.open(name_file4properties);
		if (!file4properties) {
			cout << "Error opening file" << endl;
			file4data.close();
			return;
		}
		if (seed != -1) {//if seed = -1, it means that we want to use the default seed
			rng.seed((unsigned)seed);//we assign the selected seed in file4Settings
			file4properties << "#random_choice_par \t seed \t number of generated values:" << endl;
			file4properties << random_choice_par << "\t" << seed << "\t" << N;
			//!!
			file4data << "#Generator n." << random_choice_par << " - " << N << " extractions, seed = " << seed << endl;
		}
		else {
			file4properties << "#random_choice_par \t seed \t number of generated values:" << endl;
			file4properties << random_choice_par << "\t" << "default\t" << N;
			//!!
			file4data << "#Generator n." << random_choice_par << " - " << N << " extractions, seed = default" << endl;
		}

		if (random_choice_par == 3) {
			file4data << scientific << setprecision(14);
		}
		else {
			file4data << scientific << setprecision(15);
		}

		for (int i = 0; i < N; i++) {
			file4data << distribution(rng) << endl;
		}

		file4data.close();
		file4properties.close();
	}

	//-----------------------------------------------------------------
	//SAMPLING:

	//Function to implement Box-Muller algorithm with trigonometric functions:
	template <class D> void box_muller(int seed, int n_steps, D mu, D var, char* name_file4data) {
		/*
		 -> seed = seed del generatore di numeri casuali;
		 -> n_steps = number of extractions;
		 -> mu = average of the gaussian distribution;
		 -> var = variance of the gaussian distribution;
		 -> name_file4data = name for output file;
		 */

		D x, y, g1, g2; //to construct gaussian distributions;

		//uniform distribution in [0, 1]
		uniform_real_distribution <D> dist(0, 1);

		ofstream file4data; //declaration of output file

		file4data.open(name_file4data);
		if (!file4data) {
			cout << "Error opening data file" << endl;
			return;
		}
		file4data << setprecision(numeric_limits<D>::max_digits10);//This way, the output to the file will show the double values ​​with the greatest possible precision.
		file4data << "#BOX MULLER DATA: (average = " << mu << ", variance = " << var << ", number of extractions = " << n_steps << ")" << endl;

		if (seed != -1) {//if seed = -1, it means that we want to use the default seed
			rng.seed((unsigned)seed);//we assign the selected seed in file4Settings
			file4data << "#Generator n." << random_choice_par << ", seed = " << seed << endl;
		}
		else {
			file4data << "#Generator n." << random_choice_par << ", seed = default" << endl;
		}

		file4data << "#N°point \t Value_ditribution_1  \t Value_ditribution_2" << endl;

		for (int ii = 1; ii <= n_steps; ii++) {

			//Extraction of x and y in [0,1]
			x = dist(rng); //implicit conversion int to double
			y = dist(rng);

			//x with a uniform  distribution in [0,2pi] and y with a distribution exp(-y) in [0, +inf], writing:
			x = 2 * acos(-1) * x;
			y = sqrt(-2 * log(1 - y));

			//Gaussian Distributed Points:
			g1 = y * cos(x) * sqrt(var) + mu;
			g2 = y * sin(x) * sqrt(var) + mu;

			//	cout << ii << ") value = " << y << endl;
			file4data << ii << "\t" << g1 << "\t" << g2 << endl;

		}

		file4data.close();

	}

	//Function to implement Box-Muller algorithm without trigonometric functions:
	template <class D> void box_muller_no_trig(int seed, int n_steps, D mu, D var, char* name_file4data) {
		/*
		 -> seed = seed del generatore di numeri casuali;
		 -> n_steps = number of extractions;
		 -> mu = average of the gaussian distribution;
		 -> var = variance of the gaussian distribution;
		 -> name_file4data = name for output file;
		 */

		D x, y, g1, g2; //to construct gaussian distributions;

		//uniform distribution in [0, 1]
		uniform_real_distribution <D> dist(-1, 1);

		ofstream file4data; //declaration of output file

		file4data.open(name_file4data);
		if (!file4data) {
			cout << "Error opening data file" << endl;
			return;
		}

		file4data << fixed << setprecision(numeric_limits<D>::max_digits10);//This way, the output to the file will show the double values ​​with the greatest possible precision.
		file4data << "#BOX MULLER DATA: (average = " << mu << ", variance = " << var << ", number of extractions = " << n_steps << ")" << endl;

		if (seed != -1) {//if seed = -1, it means that we want to use the default seed
			rng.seed((unsigned)seed);//we assign the selected seed in file4Settings
			file4data << "#Generator n." << random_choice_par << ", seed = " << seed << endl;
		}
		else {
			file4data << "#Generator n." << random_choice_par << ", seed = default" << endl;
		}

		file4data << "#N°point \t Value_ditribution_1  \t Value_ditribution_2" << endl;

		int step = 0;
		double r2;
		
		while (step < n_steps) {

			//Extraction of x and y in [0,1]
			x = dist(rng); //implicit conversion int to double
			y = dist(rng);

			r2 = (x * x + y * y);

			if ((r2 < 1) && (r2 > 0)) {
				step += 1;
				g1 = sqrt(-2 * log(1 - r2) / r2) * x * sqrt(var) + mu;
				g2 = sqrt(-2 * log(1 - r2) / r2) * y * sqrt(var) + mu;
				//	cout << ii << ") value = " << y << endl;
				file4data << step << "\t" << g1 << "\t" << g2 << endl;
			}
		}

		file4data.close();

	}
	
	//Function to implement a single step of Box-Muller algorithm with trigonometric functions:
	template <class T, class D> void box_muller_single_step(D& distribution, T* mu, T* var, T* result1, T* result2) {
		/*
		 -> distribution deve essere dist così definita: std::uniform_real_distribution <T> dist(0, 1);
		 -> mu = average of the gaussian distribution;
		 -> var = variance of the gaussian distribution;
		 -> result1 and result2 are the 2 numbers which are distributed with two indipendent gaussian distributions.
		 */

		T x, y; //to construct gaussian distributions;

		//Extraction of x and y in [0,1]
		x = distribution(rng); //implicit conversion int to double
		y = distribution(rng);

		//x with a uniform  distribution in [0,2pi] and y with a distribution exp(-y) in [0, +inf], writing:
		x = 2 * acos(-1) * x;
		y = sqrt(-2 * log(1 - y));

		//Gaussian Distributed Points:
		(*result1) = y * cos(x) * sqrt(*var) + (*mu);
		(*result2) = y * sin(x) * sqrt(*var) + (*mu);

	}

	//Function to implement a single step of Box-Muller algorithm without trigonometric functions:
	template <class T, class D> void box_muller_no_trig_single_step(D& distribution, T* mu, T* var, T* result1, T* result2) {
		/*
		 -> distribution deve essere dist così definita: std::uniform_real_distribution <T> dist(-1, 1);
		 -> mu = average of the gaussian distribution;
		 -> var = variance of the gaussian distribution;
		 -> result1 and result2 are the 2 numbers which are distributed with two indipendent gaussian distributions.
		 */

		T x, y, r2; //to construct gaussian distributions;

		int flag = 1;
		
		while (flag) {
			//Extraction of x and y in [0,1]
			x = distribution(rng); //implicit conversion int to double
			y = distribution(rng);

			r2 = (x * x + y * y);

			if ((r2 < 1) && (r2 > 0)) {
				flag = 0;
				(*result1) = sqrt(-2 * log(1 - r2) / r2) * x * sqrt(*var) + *mu;
				(*result2) = sqrt(-2 * log(1 - r2) / r2) * y * sqrt(*var) + *mu;
			}
		}
	}

	//Function to implement Metropolis algorithm (no Hastings):
	template <class F, class L, class T, class D> int metropolis(D dist_unif, L ptest, F prob, int seed, int n_steps, T value, string name_file4data, bool append_mode) {
		/*
			-> dist_unif = dist defined as std::uniform_real_distribution <double> dist(0, 1);
			-> ptest = probability density for value->(test value) step;
			-> prob = probability density that we want to sample;
			-> seed = seed of the random number generator. If seed = -1, it means that we DON'T want to modify the seed from the previosuly chosen one (or default if no choice was made);
			-> n_steps = number of Metropolis steps;
			-> value = starting value of Metropolis algorithm;
			-> name_file4data = name of output file: it will contain the Markov chain;
			-> append_mode = boolean variable which indicates: 1 if you want to open the output file in append mode, 0 otherwise.
		*/

		int n_acc = 0;
		T x_test, dp;

		ofstream file4data; //declaration of output file

		if (append_mode) {
			file4data.open(name_file4data, ios::app);
		}
		else
		{
			file4data.open(name_file4data);
		}

		if (!file4data) {
			cout << "Error opening file" << endl;
			return 1;
		}

		file4data << setprecision(numeric_limits<T>::max_digits10);//This way, the output to the file will show the double values ​​with the greatest possible precision.

		if (seed != -1) {//if seed = -1, it means that we DON'T want to modify the seed
			rng.seed((unsigned)seed);//we assign the selected seed in file4Settings
		}

		for (int ii = 0; ii < n_steps; ii++) {
			x_test = ptest(value);
			if ((dist_unif(rng)) < (pacc(prob(x_test)/prob(value)))) {
				n_acc++;
				value = x_test;
				file4data << value << endl;
			}
			else {
				file4data << value << endl;
			}
		}

		file4data.close();

		return n_acc;
	}

	//Function to implement a single step of Metropolis algorithm (no Hastings), in append mode:
	template <class F, class L, class T, class D> int metropolis_1s(D dist_unif, L ptest, F prob, T value, string name_file4data) {
	/*
		-> dist_unif = dist defined as std::uniform_real_distribution <double> dist(0, 1);
		-> ptest = probability density for value->(test value) step;
		-> prob = probability density that we want to sample;
		-> value = starting value of Metropolis algorithm;
		-> name_file4data = name of output file: it will contain the Markov chain;
	*/

		int n_acc = 0;
		T x_test, dp;

		ofstream file4data; //declaration of output file

		file4data.open(name_file4data, ios::app);

		if (!file4data) {
			cout << "Error opening file" << endl;
			return 1;
		}

		file4data << setprecision(numeric_limits<T>::max_digits10);//This way, the output to the file will show the double values ​​with the greatest possible precision.

		x_test = ptest(value);
		if ((dist_unif(rng)) < (pacc(prob(x_test) / prob(value)))) {
			n_acc++;
			value = x_test;
			file4data << value << endl;
		}
		else
		{
			file4data << value << endl;
		}

		file4data.close();

		return n_acc;
	}


	//Function to implement Metropolis algorithm (no Hastings), for distributions writable in an exponential form:
	template <class F, class L, class T, class D> int metroexp(D dist_unif, L ptest, F pr, int seed, int n_steps, T value, string name_file4data, bool append_mode) {
		/*
			-> dist_unif = dist defined as std::uniform_real_distribution <double> dist(0, 1);
			-> ptest = probability density for value->(test value) step;
			-> pr = - (exponent of probability ratio);
			-> seed = seed of the random number generator. If seed = -1, it means that we DON'T want to modify the seed from the previosuly chosen one (or default if no choice was made);
			-> n_steps = number of Metropolis steps;
			-> value = starting value of Metropolis algorithm;
			-> name_file4data = name of output file: it will contain the Markov chain;
			-> append_mode = boolean variable which indicates: 1 if you want to open the output file in append mode, 0 otherwise.
		*/
		int n_acc = 0;
		T x_test, dp;

		ofstream file4data; //declaration of output file

		if (append_mode){
			file4data.open(name_file4data, ios::app);
		}
		else
		{
			file4data.open(name_file4data);
		}
		if (!file4data) {
			cout << "Error opening file" << endl;
			return 1;
		}

		file4data << setprecision(numeric_limits<T>::max_digits10);//This way, the output to the file will show the double values ​​with the greatest possible precision.

		if (seed != -1) {//if seed = -1, it means that we DON'T want to modify the seed
			rng.seed((unsigned)seed);//we assign the selected seed in file4Settings
		}

		for (int ii = 0; ii < n_steps; ii++) {
			x_test = ptest(value);
			dp = pr(x_test, value);
			if (dp < 0) // the trial is accepted
			{
				n_acc++;
				value = x_test;
			}
			else if (dist_unif(rng) < exp(-dp)) // the trial is accepted
			{
				n_acc++;
				value = x_test;
			}

			file4data << value << endl;
		}

		file4data.close();

		return n_acc;
	}

	//Function to implement Metropolis algorithm (no Hastings), for distributions writable in an exponential form, WITHOUT SAVING (useful for "visual_acceptance.cpp"):
	template <class F, class L, class T, class D> int metroexp(D dist_unif, L ptest, F pr, int seed, int n_steps, T value) {
		/*
			-> dist_unif = dist defined as std::uniform_real_distribution <double> dist(0, 1);
			-> ptest = probability density for value->(test value) step;
			-> pr = - (exponent of probability ratio);
			-> seed = seed of the random number generator. If seed = -1, it means that we DON'T want to modify the seed from the previosuly chosen one (or default if no choice was made);
			-> n_steps = number of Metropolis steps;
			-> value = starting value of Metropolis algorithm;
			-> name_file4data = name of output file: it will contain the Markov chain;
			-> append_mode = boolean variable which indicates: 1 if you want to open the output file in append mode, 0 otherwise.
		*/
		int n_acc = 0;
		T x_test, dp;

		if (seed != -1) {//if seed = -1, it means that we DON'T want to modify the seed
			rng.seed((unsigned)seed);//we assign the selected seed in file4Settings
		}

		for (int ii = 0; ii < n_steps; ii++) {
			x_test = ptest(value);
			dp = pr(x_test, value);
			if (dp < 0) // the trial is accepted
			{
				n_acc++;
				value = x_test;
			}
			else if (dist_unif(rng) < exp(-dp)) // the trial is accepted
			{
				n_acc++;
				value = x_test;
			}
		}

		return n_acc;
	}

	//Function to implement a single step of Metropolis algorithm (no Hastings), for distributions writable in an exponential form, in append mode:
	template <class F, class L, class T, class D> int metroexp_1s(D dist_unif, L ptest, F pr,  T value, string name_file4data) {
		/*
			-> dist_unif = dist defined as std::uniform_real_distribution <double> dist(0, 1);
			-> ptest = probability density for value->(test value) step;
			-> pr = - (exponent of probability ratio);
			-> value = starting value of Metropolis algorithm;
			-> name_file4data = name of output file: it will contain the Markov chain;
		*/

		int n_acc = 0;
		T x_test, dp;

		ofstream file4data; //declaration of output file

		file4data.open(name_file4data, ios::app);
		
		if (!file4data) {
			cout << "Error opening file" << endl;
			return 1;
		}
		
		file4data << setprecision(numeric_limits<T>::max_digits10);//This way, the output to the file will show the double values ​​with the greatest possible precision.

		x_test = ptest(value);
		dp = pr(x_test, value);
		if (dp < 0) // the trial is accepted
		{
			n_acc++;
			value = x_test;
		}
		else if (dist_unif(rng) < exp(-dp)) // the trial is accepted
		{
			n_acc++;
			value = x_test;
		}

		file4data << value << endl;

		file4data.close();

		return n_acc;
	}

	//Function to implement Metropolis-Hastings algorithm:
	template <class F, class F1, class L, class T, class D> int hastings(D dist_unif, F prob, L ptest, T(*pacc)(T), F1 trasfer_value_to_test, F1 trasfer_test_to_value, int seed, int n_steps, T value, string name_file4data, bool append_mode) {
		/*
			-> dist_unif = dist defined as std::uniform_real_distribution <double> dist(0, 1);
			-> prob = probability density that we want to sample;
			-> ptest = probability density for value->(test value) step;
			-> pacc = probability density for acceptance step;
			-> trasfer_value_to_test = probability density for value->(test_value) step, where value and test_value are respectively our "value" and ptest(value) computed before;
			-> trasfer_test_to_value = probability density for test_value->(value) step, where value and test_value are respectively our "value" and ptest(value) computed before;
			-> seed = seed of the random number generator. If seed = -1, it means that we DON'T want to modify the seed from the previosuly chosen one (or default if no choice was made);
			-> n_steps = number of Metropolis-Hastings steps;
			-> value = starting value of Metropolis-Hastings algorithm;
			-> name_file4data = name of output file: it will contain the Markov chain;
			-> append_mode = boolean variable which indicates: 1 if you want to open the output file in append mode, 0 otherwise.
		*/

		int n_acc = 0;
		T x_test;

		ofstream file4data; //declaration of output file

		if (append_mode) {
			file4data.open(name_file4data, ios::app);
		}
		else
		{
			file4data.open(name_file4data);
		}
		
		if (!file4data) {
			cout << "Error opening file" << endl;
			return 1;
		}

		file4data << setprecision(numeric_limits<T>::max_digits10);//This way, the output to the file will show the double values ​​with the greatest possible precision.

		if (seed != -1) {//if seed = -1, it means that we DON'T want to modify the seed
			rng.seed((unsigned)seed);//we assign the selected seed in file4Settings
		}

		for (int ii = 0; ii < n_steps; ii++) {
			x_test = ptest(value);
			if ((dist_unif(rng)) < (pacc(prob(x_test) * trasfer_test_to_value(x_test, value) / (prob(value) * trasfer_value_to_test(x_test, value))))){
				n_acc++;
				value = x_test;
			}
			file4data << value << endl;
		}

		file4data.close();

		return n_acc;
	}

	//Function to implement a single step of Metropolis-Hastings algorithm (no Hastings) in append mode:
	template <class F, class F1, class L, class T, class D> int hastings_1s(D dist_unif, F prob, L ptest, T(*pacc)(T), F1 trasfer_value_to_test, F1 trasfer_test_to_value, T value, string name_file4data) {
		/*
			-> dist_unif = dist defined as std::uniform_real_distribution <double> dist(0, 1);
			-> prob = probability density that we want to sample;
			-> ptest = probability density for value->(test value) step;
			-> pacc = probability density for acceptance step;
			-> trasfer_value_to_test = probability density for value->(test_value) step, where value and test_value are respectively our "value" and ptest(value) computed before;
			-> trasfer_test_to_value = probability density for test_value->(value) step, where value and test_value are respectively our "value" and ptest(value) computed before;
			-> value = starting value of Metropolis-Hastings algorithm;
			-> name_file4data = name of output file: it will contain the Markov chain;
		*/

		int n_acc = 0;
		T x_test;

		ofstream file4data; //declaration of output file

		if (append_mode) {
			file4data.open(name_file4data, ios::app);
		}
		else
		{
			file4data.open(name_file4data);
		}

		if (!file4data) {
			cout << "Error opening file" << endl;
			return 1;
		}

		file4data << setprecision(numeric_limits<T>::max_digits10);//This way, the output to the file will show the double values ​​with the greatest possible precision.

		x_test = ptest(value);
		if ((dist_unif(rng)) < (pacc(prob(x_test) * trasfer_test_to_value(x_test, value) / (prob(value) * trasfer_value_to_test(x_test, value))))) {
			n_acc++;
			value = x_test;
		}
		
		file4data << value << endl;

		file4data.close();

		return n_acc;
	}
};


//-----------------------------------------------------------------
//ACCEPTANCE PROBABILITY:

//Function to do the accept-reject step in Metropolis(-Hastings): it isn't the only possible choiche but is the most common.
template <class D> D pacc(D x) {
	return min<D>(1, x);
}

//Function to do the accept-reject step in Metropolis(-Hastings): (to provide an alternative to the previous one)
template <class D> D pacc2(D x) {
	return x/(1+x);
}

//-----------------------------------------------------------------
//STATISTICAL ANALYSIS FUNCTIONS:

//Function to compute sample mean and sample variance of the sample mean from a given file of data, with n_skip discarded lines.
template <class T> int stats_indipendent_unbiased(T * mean, T * var_m, int n_skip, string name_input_file) {
	/*
		-> *mean will contain sample mean value;
		-> *var_m will contain sample variance value of the sample mean;
		-> n_skip = number of discarded lines from the document start;
		-> name_input_file = file from which we take data to analyze
		-> return the number of draws on which we have done the statistics 
		   if n_skip < number of lines in the document, 1 otherwise.
		
		NOTE FOR ME: page 115 of statnotes_3.3.0.pdf on Desktop.

		Welford's algorithm was used: B. P. Welford. Note on a method for 
		calculating corrected sums of squares and products. Technometrics, 
		4(3):419–420, 1962.
	*/

	ifstream input_file; //declaration of input file
	int index = 0;
	T value, delta;
	string line;

	input_file.open(name_input_file);
	if (!input_file) {
		cout << "Error opening file" << endl;
		return 1;
	}

	for (int i = 0; i < n_skip; i++) {
		if (!getline(input_file, line)) {
			cerr << "Error: there are less than " << n_skip << " lines in the file." << endl;
			return 1;
		}
	}

	(*mean) = 0;
	(*var_m) = 0;

	while (input_file >> value) {
		index++;
		delta = value - (*mean);
		(*mean) = (*mean) + delta / index;
		(*var_m) = (*var_m) + delta * (value - (*mean));
	}

	(*var_m) = (*var_m) / (index - 1);

	(*var_m) = (*var_m) / index;
	
	input_file.close();

	return index;
}

//Function to compute sample mean and sample variance of the sample mean from a given file of data, with n_skip discarded lines.
template <class T> int stats_indipendent_unbiased(T* mean, T* var_m, int n_skip, vector<T>& draws, string name_input_file) {
	/*
		-> *mean will contain sample mean value;
		-> *var_m will contain sample variance value of the sample mean;
		-> n_skip = number of discarded lines from the document start;
		-> draws will contain all the draws which we are using in the mean and variance computations;
		-> name_input_file = file from which we take data to analyze;
		-> return the number of draws on which we have done the statistics
		   if n_skip < number of lines in the document, 1 otherwise.

		NOTE FOR ME: page 115 of statnotes_3.3.0.pdf on Desktop.

		Welford's algorithm was used: B. P. Welford. Note on a method for
		calculating corrected sums of squares and products. Technometrics,
		4(3):419–420, 1962.
	*/

	ifstream input_file; //declaration of input file
	int index = 0;
	T value, delta;
	string line;

	input_file.open(name_input_file);
	if (!input_file) {
		cout << "Error opening file" << endl;
		return 1;
	}

	for (int i = 0; i < n_skip; i++) {
		if (!getline(input_file, line)) {
			cerr << "Error: there are less than " << n_skip << " lines in the file." << endl;
			return 1;
		}
	}

	(*mean) = 0;
	(*var_m) = 0;

	while (input_file >> value) {
		index++;
		delta = value - (*mean);
		(*mean) = (*mean) + delta / index;
		(*var_m) = (*var_m) + delta * (value - (*mean));
		draws.push_back(value);
	}

	(*var_m) = (*var_m) / (index - 1);

	(*var_m) = (*var_m) / index;

	input_file.close();

	return index;
}


//Function to compute sample mean and sample variance of the sample mean from a given array draws of data.
template <class T> void stats_indipendent_unbiased(T* mean, T* var_m, vector<T>& draws) {
	/*
		-> *mean will contain sample mean value;
		-> *var_m will contain sample variance value of the sample mean;
		-> draws contain all the draws which we are using in the mean and variance computations;
		
		NOTE FOR ME: page 115 of statnotes_3.3.0.pdf on Desktop.

		Welford's algorithm was used: B. P. Welford. Note on a method for
		calculating corrected sums of squares and products. Technometrics,
		4(3):419–420, 1962.
	*/

	int index = draws.size();
	T delta;

	(*mean) = 0;
	(*var_m) = 0;
	for (int ii = 0; ii < index; ii++) {
		delta = draws[ii] - (*mean);
		(*mean) = (*mean) + delta / index;
		(*var_m) = (*var_m) + delta * (draws[ii] - (*mean));
	}
		
	(*var_m) = (*var_m) / (index - 1);
	(*var_m) = (*var_m) / index;
}

template <class T> void stats_indipendent_unbiased_more_faster(T* mean, T* var_m, vector<T>& draws) {
	/*
		-> *mean will contain sample mean value;
		-> *var_m will contain sample variance value of the sample mean;
		-> draws contain all the draws which we are using in the mean and variance computations;

		NOTE FOR ME: page 115 of statnotes_3.3.0.pdf on Desktop.

		Welford's algorithm was used: B. P. Welford. Note on a method for
		calculating corrected sums of squares and products. Technometrics,
		4(3):419–420, 1962.
	*/

	int index = draws.size();
	T delta, mean2 = 0;

	(*mean) = 0;
	(*var_m) = 0;
	for (int ii = 0; ii < index; ii++) {
		
		delta = draws[ii];
		(*mean) += delta;
		mean2 += delta * delta;
	}

	(*mean) /= draws.size();
	mean2 /= draws.size();

	(*var_m) = mean2 - (*mean) * (*mean);

	(*var_m) = (*var_m) / (index - 1);
}



//Function to compute sample mean and sample variance of the sample mean from a given vector of data, with autocorrelation estimated by the blocking techinque.
template <class T> int blocking(T* mean, T* var_m, vector<T>& draws, int dim_block) {
	/*	
		-> *mean will contain sample mean value;
		-> *var_m will contain sample variance value of the sample mean;
		-> draws contain all the draws which we are using in the mean and variance computations;
		-> dim_block = block dimension used in blocking procedure.

	Uses the Welford algorithm to avoid numerical errors due to divisions between large numbers.
	This implies a slowdown. Therefore, a faster version of this function has also been also 
	implemented in "blocking_faster", although it is subject to problems when large divisions
	are present (for blocks of very large dimension).
	If you need precision, use the "blocking" version, especially for numerically sensitive datasets.
	If performance is more important and the data is not extremely sensitive to numerical errors, 
	the "blocking_faster" version is a good choice.
	*/

	(*mean) = 0;
	(*var_m) = 0;

	int n_blocks = draws.size() / dim_block;//=number of blocks
	if (n_blocks < 2) {
		std::cerr << "Warning: not enough blocks (" << n_blocks << ") to estimate variance reliably. Returning NaN." << std::endl;
		*mean = std::numeric_limits<T>::quiet_NaN();
		*var_m = std::numeric_limits<T>::quiet_NaN();
		return 0;
	}
	
	int n_max = n_blocks * dim_block;//N_max to consider to compute variance
	int index = 0;
	T mean_tmp, delta;
	for (int jj = 0; jj < n_max; jj += dim_block) {
		index++;
		mean_tmp = 0;
		for (int kk = 1; kk <= dim_block; kk++) {
			delta = draws[jj + kk - 1] - mean_tmp;
			mean_tmp = mean_tmp + delta / kk;
		}
		delta = mean_tmp - (*mean);
		//(*mean) = (*mean) + delta / (index+1);sbagliato->cancellalo poi
		(*mean) = (*mean) + delta / index;
		(*var_m) = (*var_m) + delta * (mean_tmp - (*mean));
	}
	(*var_m) = ((*var_m) / (n_blocks - 1)) / n_blocks;
	for (int ii = n_max + 1; ii < draws.size(); ii++) {
		delta = draws[ii - 1] - (*mean);
		(*mean) = (*mean) + delta / ii;
	}
	return 1;
}

//Function to compute sample mean and sample variance of the sample mean from a given vector of data, with autocorrelation estimated by the blocking techinque.
//More fast than blocking, but less precise.
template <class T> int blocking_faster(T* mean, T* var_m, vector<T>& draws, int dim_block) {
	/*
		-> *mean will contain sample mean value;
		-> *var_m will contain sample variance value of the sample mean;
		-> draws contain all the draws which we are using in the mean and variance computations;
		-> dim_block = block dimension used in blocking procedure.

	If you need precision, use the "blocking" version, especially for numerically sensitive datasets.
	If performance is more important and the data is not extremely sensitive to numerical errors, 
	the "blocking_faster" version is a good choice.
	*/
	(*mean) = 0;
	(*var_m) = 0;
	int n_blocks = draws.size() / dim_block;//=number of blocks //be careful: it is a division between int and so decimals are discarded
	if (n_blocks < 2) {
		std::cerr << "Warning: not enough blocks (" << n_blocks << ") to estimate variance reliably. Returning NaN." << std::endl;
		*mean = std::numeric_limits<T>::quiet_NaN();
		*var_m = std::numeric_limits<T>::quiet_NaN();
		return 0;
	}
	int n_max = n_blocks * dim_block;//N_max to consider to compute variance
	int index = 0;
	T mean_tmp, delta;
	for (int jj = 0; jj < n_max; jj += dim_block) {
		index++;
		mean_tmp = 0;
		for (int kk = 0; kk < dim_block; kk++) {
			mean_tmp += draws[jj + kk];
		}
		mean_tmp /= dim_block;
		delta = mean_tmp - (*mean);
		(*mean) = (*mean) + delta / index;
		(*var_m) = (*var_m) + delta * (mean_tmp - (*mean));
	}
	(*var_m) = ((*var_m) / (n_blocks - 1)) / n_blocks;
	for (int ii = n_max + 1; ii < draws.size(); ii++) {
		delta = draws[ii - 1] - (*mean);
		(*mean) = (*mean) + delta / ii;
	}
	return 1;
}

template <class T> int blocking_more_faster(T* mean, T* var_m, vector<T>& draws, int dim_block) {
	/*
		-> *mean will contain sample mean value;
		-> *var_m will contain sample variance value of the sample mean;
		-> draws contain all the draws which we are using in the mean and variance computations;
		-> dim_block = block dimension used in blocking procedure.

	If you need precision, use the "blocking" version, especially for numerically sensitive datasets.
	If performance is more important and the data is not extremely sensitive to numerical errors,
	the "blocking_faster" version is a good choice.
	*/
	(*mean) = 0;
	(*var_m) = 0;
	int n_blocks = draws.size() / dim_block;//=number of blocks //be careful: it is a division between int and so decimals are discarded
	if (n_blocks < 2) {
		std::cerr << "Warning: not enough blocks (" << n_blocks << ") to estimate variance reliably. Returning NaN." << std::endl;
		*mean = std::numeric_limits<T>::quiet_NaN();
		*var_m = std::numeric_limits<T>::quiet_NaN();
		return 0;
	}
	int n_max = n_blocks * dim_block;//N_max to consider to compute variance
	int index = 0;
	T mean_tmp, delta, mean2_tmp, mean2 = 0;
	for (int jj = 0; jj < n_max; jj += dim_block) {
		index++;
		mean_tmp = 0;
		mean2_tmp = 0;
		for (int kk = 0; kk < dim_block; kk++) {
			delta = draws[jj + kk];
			mean_tmp += delta;
		}
		mean_tmp /= (double) dim_block;
		(*mean) += mean_tmp;
		mean2 += (mean_tmp * mean_tmp);
	}
	(*mean) /= index;
	mean2 /= index;
	(*var_m) = mean2 - (*mean) * (*mean);
	(*var_m) = (*var_m) / (n_blocks - 1);
	for (int ii = n_max + 1; ii < draws.size(); ii++) {
		delta = draws[ii - 1] - (*mean);
		(*mean) = (*mean) + delta / (double) ii;
	}
	return 1;
}

//Function to compute sample mean and sample variance of an observable, described by a function of two real parameters, 
// implementing Bootstrap technique, without blocking (we are assuming to have non-autocorrelated draws or to have done
// the blocking resample in a previous step;
// Using the same implementation of blocking_more_faster.
template <class S, class D, class F, class T> bool bootstrap_2in(int n_steps, S sampler, D dist_index, F func, vector <T>& x1, vector <T>& x2, T* mean_f, T* var_f) {
	//0 if success, 1 if failure
	/*
		-> n_steps = number of iteration of bootstrap procedure;
		-> dist_index = dist defined as std::uniform_int_distribution<> dist_int(0, x1.size() - 1), useful to choose randomly the draws in our sample;
		-> func = function that defines the observable we want to study. It accepts 2 real parameters: they are sampled from x1 and x2 passed into the
				  function and ordered in f exactly as they are ordered in this call. BE CAREFUL NOT TO SWAP THEIR ORDER.
		-> x1 = vector containing the samples of the first parameter of f;
		-> x2 = vector containing the samples of the first parameter of f;
		-> mean_f = it will contain the sample mean of f;
		-> var_f = it will contain the sample variance of f (BE CAREFUL: IT IS NOT THE VARIANCE OF THE MEAN: IN THIS CONTEXT, THE SAMPLE ONE SHOULD BE USED:
				   in fact, the standard deviation of the mean decreases as the number of repetitions increases, but here the only thing I'm interested in
				   seeing is how much f fluctuates, I don't have to decrease the error)
	*/

	if (x1.size() != x2.size() || x1.empty()) {
		cerr << "Error: vectors x1 and x2 must have the same non-zero size." << endl;
		return 1;
	}

	int index, n_sample = x1.size();
	T mean_x1, mean_x2, mean_f2 = 0, value;

	(*mean_f) = 0;
	(*var_f) = 0;

	for (int ii = 0; ii < n_steps; ++ii) {

		mean_x1 = 0;
		mean_x2 = 0;

		for (int jj = 0; jj < n_sample; ++jj) {
			index = dist_index(sampler.rng);
			mean_x1 += x1[index];
			mean_x2 += x2[index];
		}
		mean_x1 /= (double)n_sample;
		mean_x2 /= (double)n_sample;

		value = func(mean_x1, mean_x2);
		(*mean_f) += value;
		mean_f2 += (value * value);
	}

	(*mean_f) /= (double)n_steps;
	mean_f2 /= (double)n_steps;
	(*var_f) = mean_f2 - (*mean_f) * (*mean_f);
	(*var_f) *= (double)n_steps;
	(*var_f) /= (double)(n_steps - 1);

	return 0;
}

//Function to compute sample means, sample variances of the sample mean and covariance between N samples of different 
// data, given in N input files which names are passed throught a vector,
// with autocorrelation estimated by the blocking techinque.
// Using the same implementation of blocking_more_faster.
template <class T> bool blocking_double_corr_more_faster(
	vector<T>& means, // -> it will refer to sample means of name_input_files data;
	vector<vector<T>>& cov_matrix, // -> it will refer to covariance matrix between name_input_files data; (FOR ME: defined as (148) of page 116 statnotes)
	vector<vector<T>>& corr_matrix, // -> it will refer to correlation matrix between name_input_files data; (FOR ME: defined as (148) of page 116 statnotes)
	vector<string>& name_input_files, // -> names of input files of data
	int dim_block, // -> dim_block = block dimension used in blocking procedure.
	int n_skip, // -> n_skip = number of data to discard from input_file;
	int n_data // -> n_data = number of data in the file; IT MUST BE THE SAME FOR each input_file
) {

	T value;
	vector<T> means_tmp, means2;
	string line;
	vector <ifstream> input_files;
	int n_blocks = n_data / dim_block;//=number of blocks
	int N = name_input_files.size();

	if (n_blocks < 2) {
		cerr << "Warning: not enough blocks (" << n_blocks << ") to estimate variance reliably. Returning NaN." << endl;
		means.assign(N, numeric_limits<T>::quiet_NaN());
		cov_matrix.assign(N, vector<T>(N, numeric_limits<T>::quiet_NaN()));
		corr_matrix.assign(N, vector<T>(N, numeric_limits<T>::quiet_NaN()));
		return 1;
	}

	int n_max = n_blocks * dim_block;//N_max to consider to compute variance

	means.resize(N, 0);
	means_tmp.resize(N, 0);
	means2.resize(N, 0);
	cov_matrix.resize(N, vector<T>(N, 0));
	corr_matrix.resize(N, vector<T>(N, 0));

	for (int ii = 0; ii < N; ++ii) {
		input_files.emplace_back(name_input_files[ii]);
		//emplace_back() creates and adds an element directly to the end of the vector,
		// building it in-place. Here it directly builds an ifstream by opening "file.txt"

		if (!(input_files.back())) {//back() returns a reference to the last element of the vector
			cerr << "Error opening " << (ii + 1) << "-th input file" << endl;
			return 1;
		}

		for (int jj = 0; jj < n_skip; ++jj) {
			if (!getline(input_files[ii], line)) {
				cerr << "Error: there are less than " << n_skip << " lines in the " << (ii + 1) << "-th file." << endl;
				return 1;
			}
		}
	}

	for (int ii = 0; ii < n_max; ii += dim_block) {
		for (int jj = 0; jj < N; ++jj) {//The order of these two for (in ii and jj) is important because we have to refer at the same blocks while sampling the covariances.
			means_tmp[jj] = 0;
			for (int kk = 0; kk < dim_block; ++kk) {
				input_files[jj] >> value;
				means_tmp[jj] += value;
			}
			means_tmp[jj] /= (double)dim_block;
			means[jj] += means_tmp[jj];
			means2[jj] += ((means_tmp[jj]) * (means_tmp[jj]));
			for (int ll = 0; ll < jj; ++ll) {
				cov_matrix[jj][ll] += means_tmp[jj] * means_tmp[ll];
			}
		}
	}

	for (int ii = 0; ii < N; ++ii) {
		means[ii] /= n_blocks;
		means2[ii] /= n_blocks;
		cov_matrix[ii][ii] = means2[ii] - (means[ii] * means[ii]);
		corr_matrix[ii][ii] = 1;
		for (int jj = 0; jj < ii; ++jj) {
			cov_matrix[ii][jj] /= n_blocks;
			cov_matrix[ii][jj] -= (means[ii] * means[jj]);
			if ((cov_matrix[ii][ii] > 0) && (cov_matrix[jj][jj] > 0)) {
				corr_matrix[ii][jj] = cov_matrix[ii][jj] / (sqrt(cov_matrix[ii][ii] * cov_matrix[jj][jj])); // corr(sampl_mean_a, sampl_mean_b) = cov(sampl_mean_a, sampl_mean_b)/sqrt(var_a * var_b)
				corr_matrix[jj][ii] = corr_matrix[ii][jj];
			}
			else {
				corr_matrix[ii][jj] = numeric_limits<T>::quiet_NaN();
				corr_matrix[jj][ii] = numeric_limits<T>::quiet_NaN();
			}
			cov_matrix[ii][jj] /= (n_blocks - 1); //cov(sampl_mean_a, sampl_mean_b) = (<sampl_mean_a*sampl_mean_b> - <sampl_mean_a> * <sampl_mean_b>)/( n - 1 )
			cov_matrix[jj][ii] = cov_matrix[ii][jj];
		}
		cov_matrix[ii][ii] /= (n_blocks - 1); //= variance[ii]
	}

	for (int ii = 0; ii < N; ++ii) {
		input_files[ii].close();
	}

	int residual = n_data % dim_block;
	if (residual != 0) {
		cout << "Note: " << residual << " data points were discarded to fit full blocks." << endl;
	}

	return 0;
}

//-----------------------------------------------------------------
//I/O FUNCTIONS:

//Function to read data from a file, in which datas are alligned as a single column.
template <class T> void read_column(int n_skip, vector<T>& draws, string name_input_file) {
	/*
		-> name_input_file = name of the file from where we read the values that we want to put into draws;
		-> n_skip = number of lines to discard from the start of name_input_file;
		-> draws will contain all the values read from name_input_file.
	*/
	ifstream input_file;
	string line;
	T value;

	input_file.open(name_input_file);
	if (!input_file) {
		cout << "Error opening input file" << endl;
		return;
	}

	for (int i = 0; i < n_skip; i++) {
		if (!getline(input_file, line)) {
			cerr << "Error: there are less than " << n_skip << " lines in the file." << endl;
			return;
		}
	}

	while (input_file >> value) {
		draws.push_back(value);
	}

	input_file.close();
}

//Function to read data one columns (which_column=0,1 in order) from a file, in which datas are alligned as 2 columns.
template <class T> void read_1from2columns(bool which_column, int n_skip, vector<T>& draws, string name_input_file) {
	/*
		-> name_input_file = name of the file from where we read the values that we want to put into draws;
		-> n_skip = number of lines to discard from the start of name_input_file;
		-> draws will contain all the values read from name_input_file.
	*/
	ifstream input_file;
	string line;
	T value, discard;

	input_file.open(name_input_file);
	if (!input_file) {
		cout << "Error opening input file" << endl;
		return;
	}

	for (int i = 0; i < n_skip; i++) {
		if (!getline(input_file, line)) {
			cerr << "Error: there are less than " << n_skip << " lines in the file." << endl;
			return;
		}
	}

	if (which_column) {
		while (input_file >> discard >> value) {
			draws.push_back(value);
		}
	}
	else {
		while (input_file >> value >> discard) {
			draws.push_back(value);
		}
	}

	input_file.close();
}

//Function to read data 2 columns from a file, in which datas are alligned as 2 columns.
template <class T, class T1> void read_2columns(int n_skip, vector<T>& draws1, vector<T1>& draws2, string name_input_file) {
	/*
		-> name_input_file = name of the file from where we read the values that we want to put into draws;
		-> n_skip = number of lines to discard from the start of name_input_file;
		-> draws will contain all the values read from name_input_file.
	*/
	ifstream input_file;
	string line;
	T value1;
	T1 value2;

	input_file.open(name_input_file);
	if (!input_file) {
		cout << "Error opening input file" << endl;
		return;
	}

	for (int i = 0; i < n_skip; i++) {
		if (!getline(input_file, line)) {
			cerr << "Error: there are less than " << n_skip << " lines in the file." << endl;
			return;
		}
	}

	while (input_file >> value1 >> value2) {
		draws1.push_back(value1);
		draws2.push_back(value2);
	}

	input_file.close();
}