docker_clean:
	-sudo docker rm \
	`sudo docker ps --all | grep -i -v container | awk '{print $$1}'`
	-sudo docker rmi \
	`sudo docker images --all | grep -i -v container | awk '{print $$3}'`

docker_show:
	- sudo docker ps --all
	- sudo docker images --all

build:
	sudo docker build . --rm -f dockerfile \
		-t awsteiner/eos:alpha2 > build.out 2>&1
	sudo docker run awsteiner/eos:alpha2 sh -c \
		"curl -o fid_3_5_22.o2 https://isospin.roam.utk.edu/public_data/eos_tables/du21/fid_3_5_22.o2 && mpirun --allow-run-as-root -np 1 ./eos_nuclei -load fid_3_5_22.o2 -point-nuclei 0.16 0.4 1.0" \
		> test_alpha2.out 2>&1

push:
	sudo docker push awsteiner/eos:alpha2
