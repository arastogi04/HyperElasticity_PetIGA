#!/bin/bash

make HyperElasticity

./HyperElasticity -iga_elements 5 iga_degree 2 -snes_max_it 80 -snes_atol 1e-5 -snes_type newtontr -snes_monitor -snes_converged_reason -iga_view

python create_vtk.py
