# Usage:
# make        # compile all binary
# make clean  # remove ALL binaries and objects
# $@ --> target , $< first prerequisite , $^ --> prerequisites

.PHONY = all clean

CC = g++ -std=c++23		# compiler to use
LIBS = -lm $(LAPACK_BLAS)	# external libraries
LAPACK_BLAS = -L/home/kingofbroccoli/Programmazione/My_C++/Libraries/ -llapacke -llapack -lrefblas -lgfortran
OPT = -O3 -march=native -mtune=native -fopenmp
CFLAGS = -Wall # -g	# Use -g just for debug			# compiler flags: all warnings + debugger meta-data
EIGEN = /home/kingofbroccoli/Programmazione/My_C++/Libraries/eigen-3.4.0/		# Path to Eigen # I could use -I$(EIGEN) instead of having the path in #include

SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:%.cpp=%.o)
BINS := Ecosystem_Extract_and_Evolve
#BINS := $(SRCS:%.c=%)

all: $(BINS)

# This default rule compiles the executable program
$(BINS): $(OBJS)
	@echo "Checking.."
#	@echo $(BINS)
#	@echo $(OBJS)
#	@echo $(SRCS)
	$(CC) $^ $(OPT) $(CFLAGS) $(LIBS) -o $@
	chmod a+x $@

# This rule compiles each module into its object file
%.o: %.cpp
	@echo "Creating object.."
	$(CC) $(OPT) $(CFLAGS) -c $< -o $@ $(LIBS)

clean:
	@echo "Cleaning up..."
	rm -rvf *.o $(BINS)
