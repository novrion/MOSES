#!/bin/bash

# change to correct dir
cd ..

# define files
orderfile="data/in/order"
infile="data/in/in"
outfile="data/out"

# clear files
> $infile
> $outfile

# assemble input files
for line in $(cat $orderfile | tr -d " \r")
do

  #tmp=${line//\r\n}
  #if cmp -s "$line" "$tmp";
  #then
  #  echo $line
  #fi

  cat $line >> $infile

  #printf "$line\n" >> $infile
done

# compile
g++ -o bin/run.exe src/main.cpp src/MOSES.cpp src/DOM.cpp src/LABOR.cpp src/QUARTER.cpp src/YEAR.cpp

# run
./bin/run.exe < $infile > $outfile
