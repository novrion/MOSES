# MOSES-1976
MOSES-1976 is the preliminary MOSES C++ code written in 2023. Very small modifications and a smaller refactor has been made since then. The code is based on the pseudocode in **_A MICRO-MACRO INTERACTIVE SIMULATION MODEL OF THE SWEDISH ECONOMY_ (Gunnar Eliasson, 1976)**.

The model is susceptible to changes in input data, which can cause extreme outcomes. The model can create exponential increases of values that is unrealistic if inconsistent input data is utilised.

# Running The Model
MOSES-1976 is executed with **scripts/moses.sh**.

The model requires the **g++ compiler** and **bash** to execute.

Results are written to **data/out**. A model run is examplified in **data/out**. Its input data has been hidden.

**moses.sh** compiles input data from **data/** to **data/in**, which the program reads from line by line for each variable. The current model does not have a robust system for handling data. Use the function **INPUT** in **src/main.cpp** and **moses.sh** as a reference for how input data should be formatted. Note the formatting of data at the end of **INPUT**, which is used to facilitate the use of data from the original MOSES APL program.

The model has only been tested on Linux and may have file-related bugs on Windows. The **LF** and **CRLF** line ending difference may cause issues on Windows.
