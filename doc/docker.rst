EOS with Docker
===============

Docker images are periodically updated which contain a working version
of the code and examples. To use the images (available at
https://hub.docker.com/repository/docker/awsteiner/e4mma/general),
just use, for example::

  sudo docker pull awsteiner/e4mma:alpha4_ju_o930a2_u24.04
  sudo docker run -it `sudo docker images | grep "alpha4_ju_o930a2_u24.04" | awk '{print $3}'`

Currently the code reads the EOS table and depending on the baryon
density :math:`n_B`, electron fraction :math:`Y_e` and temperature
:math:`T` computes various nuclear properties.


