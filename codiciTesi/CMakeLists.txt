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

target_link_libraries(polyff ${ROOT_LIBRARIES})
target_link_libraries(autocorr_check ${ROOT_LIBRARIES})
target_link_libraries(simple_computation ${ROOT_LIBRARIES})
target_link_libraries(visual_monte_carlo_history ${ROOT_LIBRARIES})
target_link_libraries(visual_polyffT ${ROOT_LIBRARIES})