# CMakeList.txt: progetto CMake per codiciTesi. Includere l'origine e definire
# qui la logica specifica del progetto.
#

find_package(ROOT REQUIRED COMPONENTS RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript MathCore Physics Matrix)

# Include ROOT directories
include_directories(${ROOT_INCLUDE_DIRS})

# Aggiungere l'origine all'eseguibile di questo progetto.
add_executable(polyff "polyff.cpp")
add_executable(autocorr_check "autocorr_check.cpp")
add_executable(simple_computation "simple_computation.cpp")
add_executable(visual_monte_carlo_history "visual_monte_carlo_history.cpp")
add_executable(visual_polyffT "visual_polyffT.cpp")
add_executable(visual_susceptibility_polyff "visual_susceptibility_polyff.cpp")
add_executable(converter_into_obs "converter_into_obs.cpp")
add_executable(autocorr_check_simpler_files "autocorr_check_simpler_files.cpp")
add_executable(susceptibility_with_errors "susceptibility_with_errors.cpp")
add_executable(visual_susceptibility_with_errors "visual_susceptibility_with_errors.cpp")
add_executable(dependent_data_analysis_autocorrelation_calculated "dependent_data_analysis_autocorrelation_calculated.cpp")
add_executable(visual_autocorrelation "visual_autocorrelation.cpp")
add_executable(fit_polyffT "fit_polyffT.cpp")
add_executable(temperature_compute "temperature_compute.cpp")
add_executable(analysis_rhok_vs_T "analysis_rhok_vs_T.cpp")
add_executable(correction_confid "correction_confid.cpp")
add_executable(analysis_rhok_vs_T_fit "analysis_rhok_vs_T_fit.cpp")
add_executable(autocorr_check_rhok "autocorr_check_rhok.cpp")
add_executable(rho_vs_T "rho_vs_T.cpp")

target_link_libraries(polyff ${ROOT_LIBRARIES})
target_link_libraries(autocorr_check ${ROOT_LIBRARIES})
target_link_libraries(simple_computation ${ROOT_LIBRARIES})
target_link_libraries(visual_monte_carlo_history ${ROOT_LIBRARIES})
target_link_libraries(visual_polyffT ${ROOT_LIBRARIES})
target_link_libraries(visual_susceptibility_polyff ${ROOT_LIBRARIES})
target_link_libraries(autocorr_check_simpler_files ${ROOT_LIBRARIES})
target_link_libraries(susceptibility_with_errors ${ROOT_LIBRARIES})
target_link_libraries(visual_susceptibility_with_errors ${ROOT_LIBRARIES})
target_link_libraries(visual_autocorrelation ${ROOT_LIBRARIES})
target_link_libraries(fit_polyffT ${ROOT_LIBRARIES})
target_link_libraries(analysis_rhok_vs_T ${ROOT_LIBRARIES})
target_link_libraries(analysis_rhok_vs_T_fit ${ROOT_LIBRARIES})
target_link_libraries(autocorr_check_rhok ${ROOT_LIBRARIES})
target_link_libraries(rho_vs_T ${ROOT_LIBRARIES})