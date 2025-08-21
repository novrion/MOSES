# MOSES - Model Of The Swedish Economic System

This repository holds C++ models of the **Model Of The Swedish Economic System**.

The aspiration of this project is to replicate the most recent APL MOSES model in object code, and to convert the MOSES household sector to micro. The nature of object-oriented programming should facilitate the micro conversion of the household sector significantly.

### How to run the model
1. **Compilation**: run **_make_** in the main directory of the project in a bash terminal.\
   (Project compilation requires **_g++ 17_** or later and **_GNU Make_**.

2. Run the compiled binary **_.test_run_** using the command: **_./.test_run_**.\
   (Run **_make clean_** before recompilation to ensure changes are represented properly)

### Branches

**MOSES-1989**  
Based on pseudocode in **_MOSES Code_ (James W. Albrecht et al., 1989)** and the current MOSES APL project.  

**MOSES-1976** (old)  
Based on pseudocode in **_A MICRO-MACRO INTERACTIVE SIMULATION MODEL OF THE SWEDISH ECONOMY_ (Gunnar Eliasson, 1976)**.
