/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
****      CONFIG.h: this file contains the configuration      ****
****        variables, which can be modified as needed.       ****
****                (author = Serena Bruzzesi)                ****
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

#include "headers.h"

//random_choice_par = parameter to choose the random number generator:
//	-> 1 = Mersenne Twister 19937 generator (64 bit)
//	-> 2 = Permuted congruential generator (pcg) (64 bit) [2014 M.E. O'Neill / pcg-random.org]
//	-> 3 = 48-bit RANLUX generator by Martin Lüscher and Fred James(1994)
#define random_choice_par 1
