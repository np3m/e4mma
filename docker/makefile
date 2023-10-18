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
		-t awsteiner/eos:alpha > build.out 2>&1 &

push:
	sudo docker push awsteiner/eos:alpha
