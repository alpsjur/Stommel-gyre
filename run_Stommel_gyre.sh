#!/bin/bash

# Define the filename variable
FILENAME="Stommel_gyre.jl"

nice julia --project=. $FILENAME 
