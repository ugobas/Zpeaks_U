PRG=ZPeaks_U

# 64 bits
CFLAGS=-g -O2 -march=nocona # for 64 bit compilation 
#CFLAGS=-g 

CC=gcc # c compiler

# Sources and objects

# HEAD= .h    	# header files
CSRC = ZPeaks_U.c optimization.c nrutil.c random3.c allocate.c Genome_profile_aux.c Peaks_aux.c cluster_peaks.c

# Likelihood_score.c optimization.c Get_peaks.c HMM_aux.c  cluster_score.c
COBJ = ZPeaks_U.o optimization.o nrutil.o random3.o allocate.o  Genome_profile_aux.o Peaks_aux.o cluster_peaks.o
# Likelihood_score.o optimization.o Get_peaks.o HMM_aux.o cluster_score.o

$(PRG): $(COBJ)
	$(CC) $(CFLAGS) $(OFLAGS) -o $(PRG) $(COBJ) -lm
#	rm -f *.o
	
	echo "Executable file generated:" $(PRG)
	

$(COBJ) : $(HEAD) $(CSRC)
	$(CC) $(CFLAGS) $(OFLAGS) -c  $(CSRC)

