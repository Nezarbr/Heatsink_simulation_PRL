For the Makefile:
make all: Compiles all .c files with the compilation command.
make clean: Deletes txt files and executables (txt files are deleted to avoid cluttering folders after executing codes with checkpoint).

IF THE MAKEFILE DOES NOT WORK ON GRID5000: Open the Makefile on Grid5000 and indent with tabs where needed. (a rather odd issue)

Reservation command: oarsub -l host=4/core=16 -I
Compilation command: mpicc -Wall -O3 -o file_name file_name.c -lm
Execution command: mpiexec -n [number_of_cores] --mca pml ^ucx --hostfile $OAR_NODEFILE ./file_name