#!/bin/bash
# Execute this file to recompile locally
c++ -Wall -shared -fPIC -std=c++11 -O2 -I/usr/local/lib/python3.6/dist-packages/ffc/backends/ufc -I/usr/local/include -I/usr/local/slepc-32/include -I/usr/local/petsc-32/include -I/usr/include/mpich -I/usr/include/hdf5/mpich -I/usr/include/eigen3 -I/home/fenics/.cache/dijitso/include ffc_form_74d937094ade8f29f005e97b083181d387b7a9f6.cpp -L/home/fenics/.cache/dijitso/lib -Wl,-rpath,/home/fenics/.cache/dijitso/lib -ldijitso-ffc_element_e8f528befdb49b0b6ffad88d40def4bf76e92d88 -ldijitso-ffc_element_a0486fd8c8124acd09b2ea68e480deaef4ac7177 -ldijitso-ffc_element_e6169480e17fdc10385e914338ba04e90d06d050 -ldijitso-ffc_element_16af6b649fadf61cfbbc905ea52262bc71249c2b -ldijitso-ffc_element_c6d1f8363bada3fbb8b2fdbd5d434c66ca66441f -ldijitso-ffc_coordinate_mapping_bbe5001d60345544f4153023a4a8d8259bf5e4fa -olibdijitso-ffc_form_74d937094ade8f29f005e97b083181d387b7a9f6.so