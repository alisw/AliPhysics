#!/bin/bash
#--------------------------------------------------
# CONFIG_BEAMS:  pp or PbPb
# CONFIG_MESON_PDG: 111, 221, 223, 331 
#--------------------------------------------------

export CONFIG_BEAMS=PbPb
export CONFIG_MESON_PDG=221
export CONFIG_SEED=20120609

aliroot -b -q sim.C
