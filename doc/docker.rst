Using Docker
============

Docker images are periodically updated which contain a working version
of the code and examples. These images are large, but they contain a
full version of E4MMA together with O2scl, O2sclpy, Scipy, Matplotlib,
and Jupyter (Scikit-learn, Tensorflow, and PyTorch are also included).
To use the images (available at
https://hub.docker.com/repository/docker/awsteiner/e4mma/general), you
can try, for example::

  sudo docker pull awsteiner/e4mma:alpha5_ju_o930a4_u24.04
  sudo docker run -p 8888:8888 -it \
    `sudo docker images | grep "alpha5_ju_o930a4_u24.04" | awk '{print $3}'`

Then, inside the docker image, you can either run the code directly::

  ./eos_nuclei -help

run the examples::

  cd examples
  ./A_point_nuc.scr

or start a Jupyter server::
  
  source docker/start_server

and then you can start a Jupyter notebook using the link specified by
the server log. 


