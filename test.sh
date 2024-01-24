# Pull the docker image
docker pull nostrad1/utk-eos:v2
# Run a test script
docker run nostrad1/utk-eos:v2 sh -c \
	 "mpirun --allow-run-as-root -np 1 ./eos_nuclei -load data/fid_3_5_22.o2 -point-nuclei 0.16 0.4 1.0" \
# Run create MUSES table function
#docker run nostrad1/utk-eos:v2 sh -c \
#     "make mbnuc"
