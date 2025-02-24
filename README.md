# ProgrammingProject_Homework_A_Serrano
**Author:** Álvaro Serrano López de la Vieja

This github is used by Álvaro Serrano to present the resolution of the ProgrammingProject Homework.

The Python script presented aims to satisfy the requeriments of a program to perform structural optimization (geometry optimization) for molecular systems, using both Cartesian and internal coordinates. Here it is important to mention that this script is based on solving saturated hydrocarbon systems and will not work with other systems. The algorithm used to do the geometry optimization is the BFGS Algorithm. The final remarks are that all the steps proposed at the homework layout has been followed and achieved.

**A breakdown of each file contained in the repository can be found below.**

- **geom_opt_notes_2024.pdf:** This is the pdf file followed in order to fulfill the homework. Important ideas inside the pdf file are highlighted.

- **inputs directory:** This directory stores all the inputs used to test the code.

- **outputs directory:** This directory stores all the outputs provided by the professor and used to compare with the results of the code.

- **geometry_optimization_A_Serrano.py:** This is the Python code presented as the resolution of the homework.




# Code Explanation

This section is used to explain the key features of the code.


**Below a breakdown of the code can be found.**

## -Reading in the input-
The first thing to do for our program to work is to read in the number of atoms, coordinates and mass for the specific molecule (input). The corresponding lines that perform these actions can be found next.
