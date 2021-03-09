#!/bin/bash

make HyperElasticity

./HyperElasticity -neohook -snes_max_it 50 -snes_atol 1e-10 -snes_type newtonls -snes_monitor -snes_converged_reason -iga_view

python create_vtk.py
