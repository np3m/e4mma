Using Docker
============

Docker images are periodically updated which contain a working version
of the code and examples. These images are large, but they contain a
full version of E4MMA together with O2scl, O2sclpy, Scipy, Matplotlib,
and Jupyter (Scikit-learn, Tensorflow, and PyTorch are also included).
To use the images (available at the `E4MMA docker repository
<https://hub.docker.com/repository/docker/awsteiner/e4mma/general>`_),
you can try, for example::

  sudo docker pull awsteiner/e4mma:alpha5_ju_o930a4_u24.04
  sudo docker run --rm -p 8888:8888 -it awsteiner/e4mma:alpha5_ju_o930a4_u24.04

Then, inside the docker image, you can either run the code directly::

  /home/np3m/e4mma/eos_nuclei -help

run the examples::

  cd /home/np3m/e4mma/examples
  ./A_point_nuc.scr

or start a Jupyter server::
  
  source /home/np3m/e4mma/docker/start_server

and then you can start a Jupyter notebook using the link specified by
the server log. If you want to transfer files between the docker
container and your directory on the host machine, you can try a
"bind mount", e.g.::

  sudo docker run -rm -p 8888:8888 --mount \
  type=bind,source=/path/on/host,target=/home/np3m/e4mma/mounted_dir \
  -it awsteiner/e4mma:alpha5_ju_o930a4_u24.04

To get more information about how bind mounts work in docker, see
`this page <https://docs.docker.com/engine/storage/bind-mounts/>`_.
Keep in mind, inside the container, you will need to use ``sudo`` to
create or modify files in the bound directory (no superuser password
is needed).
