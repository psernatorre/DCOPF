# DCOPF

Follow these steps to run the Optimal Power Flow successfully:

1) Download the folders: Input data, Source Code, Results. Save all the folders inside another one, for example "SernaTorre_project".

2) Go to the folder ``Source code''. Open the file ``directory paths.csv". In the column 3, modify the first row "Main folder name" according to the path where the folder Source Code, Input data, and Results are located. For example: /home/paul/Documents/SernaTorre_project

3) Go to the folder ``Input data". Modify all the information according to the characteristics and parameters of the power system you want to run. All the files have the format ERCOT 120-bus 500kV system.

4) Open a Julia terminal from the folder "Source code". Execute the command include >include("Lossy DC PowerFlow with BESS.jl"). See the example:

	paul@paul-Inspiron-15-5510:~/Documents/SernaTorre_project/Source code$ julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.6.3 (2021-09-23)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

	julia> include("Lossy DC PowerFlow with BESS.jl")


5) After the execution, files (.csv) and plots (.pdf) are created automatically and saved inside the folder ``Results". These files contain the optimal solution.
