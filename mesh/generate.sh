#!/bin/bash

for file in geo/*.geo; do
    filename=$(basename "$file" .geo)
    gmsh "geo/${filename}.geo" -3 -o "${filename}.msh"
done