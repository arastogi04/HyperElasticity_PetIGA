#!/bin/bash

make HyperElasticity

./HyperElasticity -neohook -iga_load thin_film_3d.dat -snes_max_it 50 -snes_atol 1e-10 -snes_type newtonls -snes_monitor -snes_converged_reason -st_type sinvert -eps_target 0.0 -eps_target_magnitude -iga_view -malloc_debug

python create_vtk.py
