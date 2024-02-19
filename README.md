# DCOPF

Follow these steps to run the Optimal Power Flow successfully:

1) Download the folders: Input data, Source Code, Results. Save all the folders inside another one, for example "SernaTorre_project".

2) Go to the folder Source code. Open the file directory paths.csv. In the column 3, modify the first row Main folder name according to the path where the folder Source Code, Input data, and Results are located. For example: /home/paul/Documents/SernaTorre_project

3) Go to the folder Input data. Modify all the information according to the characteristics and parameters of the power system you want to run. All the files have the format ERCOT 120-bus 500kV system.

4) Open a Julia terminal from the folder "Source code". Type in the terminal: include("DCOPF.jl"). You can also run the file in VSCode. In VSCode, you may need to click on File, then Open Folder, open the folder Source Code, and finally run DCOPF.jl

5) After the execution, files (.csv) and plots (.pdf) are created automatically and saved inside the folder ``Results". These files contain the optimal solution.
