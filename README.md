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

**This printes the results using Internal Coordinates, in case you want to get the results using Cartesian Coordinates it would be necessary to uncomment the lines for Cartesian Coordiantes and now comment all the lines for the Internal Coordinates.**

# Explanation of the code
This section is used to shortly and quickly explain the code.

- **parse_molecule_file funtion:** This function is used to get the atoms, coordinates and bonds arrays from an input. It devides all the lines of the file excluding the first one in parts (columns) and depending the length of columns it saves the data in one of the arrays.

- **calculate_bond_distances funtion:** In the bonds array the bonds are already specified. This function calculates the distance of each bond by going through all the bonds and using the equation:

     ![eq](https://github.com/user-attachments/assets/fe834418-c5a6-4e95-bfdd-4342541c76be)

- **calculate_bond_energies funtion:** This calculates the bond energies by making us of the bonds, atoms and distances arrays. Something to pay attention to is that in order to get the actual atom we need to substract 1 to the bond index. Python starts at 0. The equation used is:

     ![eq2](https://github.com/user-attachments/assets/57ee26da-5f3b-48d0-9b14-72ffab654946)

- **calculate_angles funtion:** This function calculates the present angles making sure they are not repeated and checking for a common atom. The equation used can be found on the pdf.

- **calculate_angle_energies funtion:** This function calculates the energy of the angles using some parameters and the following equation:

     ![eq3](https://github.com/user-attachments/assets/d84cbe24-dd69-44d0-b235-d229fbf170a7)

- **calculate_dihedral_angles funtion:** This function calculates the dihedral angles using a similar logic than before for the angles.

- **calculate_dihedral_energies funtion:** This function calculates the dihedral energy terms as the following equation:




