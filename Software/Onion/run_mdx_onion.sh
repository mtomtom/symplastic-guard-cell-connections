#!/bin/bash

# Define the model path and library
MODEL="Model/CCF/01 FEM Membranes"
LIBRARY="usrLibFemMembranes.so"
MDXV_FILE="FemMembranes.mdxv"

# Run the MorphoDynamX command
MorphoDynamX --model "$MODEL" --addlibrary "$LIBRARY" "$MDXV_FILE"
