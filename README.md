# ProgrammingProject_Homework_A_Serrano
**Author:** Álvaro Serrano López de la Vieja

This github is used by Álvaro Serrano to present the resolution of the ProgrammingProject Homework.

The Python script presented aims to satisfy the requeriments of a program to perform structural optimization (geometry optimization) for molecular systems, using both Cartesian and internal coordinates. Here it is important to mention that this script is based on solving saturated hydrocarbon systems and will not work with other systems. The algorithm used to do the geometry optimization is the BFGS Algorithm. The final remarks are that all the steps proposed at the homework layout has been followed and achieved.

**A breakdown of each file contained in the repository can be found below.**

- **geom_opt_notes_2024.pdf:** This is the pdf file followed in order to fulfill the homework. Important ideas inside the pdf file are highlighted.

- **inputs directory:** This directory stores all the inputs used to test the code.

- **outputs directory:** This directory stores all the outputs provided by the professor and used to compare with the results of the code.

- **geometry_optimization_A_Serrano.py:** This is the Python code presented as the resolution of the homework.

# How do I use the code?
This section is used to explain how to run the code.

First of all it would be necessary for you to download the Python script and the inputs directory. Once downloaded both files (script and directory) you will need to get them together in the same directory (1 directory containing the inputs directory and also containing the Python script). Then you will run the code by doing:

```console
python3 geometry_optimization_A_Serrano.py
```

Now the code will ask for an input inside the inputs directory, in case you want to study the ethane structure you will need to write ethane.mol2.

The results are printed in the terminal.
