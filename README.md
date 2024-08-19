This branch holds the preliminary MOSES C++ code written in 2023. Very small modifications and a smaller refactor has been made since then. The code is based on pseudocode from 1976.
The model is very suceptible to changes in input data, which can cause extreme outcomes. The model can create exponential increases of values that is unrealistic if inconsistent input data is utilised.
No more time will be allocated to fine-tuning this outdated and preliminary C++ model based on the pseudocode from 1976. Instead, time will be used to create a more stable and realistic MOSES model based on the 1989 pseudocode.
The time spent coding this model is not wasted as it proves the function of MOSES in object code and creates a foundation for the next version. This is since the pseudocode from 1989 is in many ways similar to the pseudocode from 1976.

The "bin" folder holds binaries. The model is executed by running the "moses.sh" script in the "scripts" folder. Input data is in the folder "data". The source code is in the "src" folder. The output is written to the file "out" in the "data" folder. 
The g++ compiler and bash is required to run the model. 
A model run output is examplified in the data/out file. The data used for the run has been hidden. 

The model has only been tested on Linux and may have file-related bugs on Windows. The LF and CRLF line ending difference may cause issues on Windows.
