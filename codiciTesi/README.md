## CONTENT DESCRIPTION:

- "autocorr_check.cpp": implements the blocking technique and ranges over different block sizes dim_block, returning the plot of variance vs dim_block for all files specified in file_list.txt.

- "polyff.cpp": this script takes the fermionic and bosonic measurement files for each temperature and returns a file with the temperature along with the mean and standard deviation for these observables. When selected, blocking is applied with block sizes selectable from an external file.

- "simple_computation".cpp: simple script to compute mass of quark in MeV.

- "visual_monte_carlo_history.cpp": visualization of Monte Carlo history.

- "visual_polyffT.cpp": visualization of the mean values ​​of the observables (with corresponding standard deviation) as a function of the temperature.

- "visual_susceptibility_polyff.cpp": this script provides a visualization of the susceptibility of Polyakov loops and chiral condensates, estimated directly from the variances computed in polyff.cpp

- "converter_into_obs.cpp": script to convert value data in observable files to obs(value), arranged in two-column (one for conf_id, one for obs(value)).

- "autocorr_check_simpler_files.cpp":  implements the blocking technique and ranges over different block sizes dim_block, returning the plot of variance vs dim_block for all files specified in file_list_therm.txt, in the 2-columns format.

- "susceptibility_with_errors.cpp": script dedicated to the estimation of the susceptibility of Polyakov loops and chiral condensates, with related error (use of blocking and bootstrap).

- "visual_susceptibility_with_errors.cpp": visualization of the chiral condensate, with errors, as a function of temperature (input file from susceptibility_with_errors.cpp).

- "dependent_data_analysis_autocorrelation_calculated.cpp": direct calculation of the autocorrelation function.

- "visual_autocorrelation.cpp": visualization of autocorrelation as a function of distance between indices, estimation of autocorrelation time;

- "fit_polyffT.cpp": script dedicated to the fit of Polyakov loops and chiral condensates as a function of temperature.

- "temperature_compute.cpp": simple script to compute temperature in MeV, given Nt e a[fm] values.

- "autocorr_check_rhok.cpp": Starting from dedicated files ("corrected_mon.dat"), the number nk of monopoles is calculated for each number k of windings, for each gauge configuration. Then, we compute nk/n1 for each k, estimate its variance by varying the dimension of blocks in blocking procedure. In output: images and file with <nk/n1> vs dim_block.

- "analysis_rhok_vs_T.cpp": Starting from dedicated files ("corrected_mon.dat"), the number nk of monopoles is calculated for each number k of windings, for each gauge configuration. Then, averaging over all configurations, the mean of nk (at a given k), the mean of the ratio between nk and n1, the ratio between the means of nk and n1 are estimated. Two plots are returned for the last 2 quantities mentioned as a function of k.

- "correction_confid.cpp": to correct files from confid erroneously duplicated.

- "analysis_rhok_vs_T_fit.cpp": Best fit of output .txt results of analysis_rhok_vs_T.cpp.

- "analysis_rhok_vs_T.cpp": Starting from dedicated files ("corrected_mon.dat"), computes total monopole density per configuration. Averages over configurations at each temperature T. Plots <rho_tot/T^3> with its error vs temperature.

- "cluster_lenght.cpp": ...