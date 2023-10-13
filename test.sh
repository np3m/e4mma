docker run -it --rm --name=eosv2 utk_eos2:latest \
    curl https://isospin.roam.utk.edu/public_data/eos_tables/du21/fid_3_5_22.o2 --output data/fid_3_5_22.o2 && make mbmuses\
    make mbmuses