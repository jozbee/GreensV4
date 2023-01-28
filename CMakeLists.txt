cmake_minimum_required(VERSION 3.24)
project(secomb_green_cpu)

set(CMAKE_CXX_STANDARD 17)

add_executable(greens analyzenet.cpp bicgstab.cpp blood.cpp cmgui.cpp contour.cpp contr_lines.cpp
        contr_shade.cpp convect.cpp eval.cpp gaussj.cpp greens.cpp histogram.cpp initgreens.cpp
        input.cpp ludcmp.cpp nrutil.cpp outboun.cpp picturenetwork.cpp postgreens.cpp
        postgreens.cpp.dat putrank.cpp setuparrays0.cpp setuparrays1.cpp setuparrays2.cpp
        testconvect.cpp tissrate.cpp main.cpp)

# move .dat files into ${PROJECT_BINARY_DIR}
configure_file(ContourParams.dat ContourParams.dat)
configure_file(IntravascRes.dat IntravascRes.dat)
configure_file(Network.dat Network.dat)
configure_file(PostGreensParams.dat PostGreensParams.dat)
configure_file(SoluteParams.dat SoluteParams.dat)
configure_file(VaryParams.dat VaryParams.dat)
configure_file(postgreens.cpp.dat postgreens.cpp.dat)
configure_file(tissrate.cpp.dat tissrate.cpp.dat)
