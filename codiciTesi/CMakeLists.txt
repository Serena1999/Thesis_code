# CMakeList.txt: progetto CMake per codiciTesi. Includere l'origine e definire
# qui la logica specifica del progetto.
#

find_package(ROOT REQUIRED COMPONENTS RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript MathCore Physics Matrix)

# Include ROOT directories
include_directories(${ROOT_INCLUDE_DIRS})

# Aggiungere l'origine all'eseguibile di questo progetto.
add_executable(poly "poly.cpp")
add_executable(autocorr_check "autocorr_check.cpp")
add_executable(simple_computation "simple_computation.cpp")

target_link_libraries(poly ${ROOT_LIBRARIES})
target_link_libraries(autocorr_check ${ROOT_LIBRARIES})
target_link_libraries(simple_computation ${ROOT_LIBRARIES})