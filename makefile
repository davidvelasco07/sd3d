FLAGS= -std=c++11 -O2 -lineinfo -I../include 

LINK=-lcublas -lcuda

OBJS= bin/main.o bin/global.o bin/initial_conditions.o bin/variables.o bin/outputs.o bin/polynomials.o bin/time_step.o bin/hydro.o bin/riemann.o bin/transforms.o bin/mesh.o bin/ader.o bin/boundary.o bin/build_comms.o bin/comms.o bin/trouble_detection.o bin/godunov_2O.o bin/correct_fluxes.o bin/split.o

exe: $(OBJS) 
	CC $(OBJS) -o sd3d -I$(CRAY_MPICH2_DIR)/include 

bin/%.o: src/%.cpp
	CC $(FLAGS) -c $< -o $@

clean :
	rm -f bin/*
	rm sd3d