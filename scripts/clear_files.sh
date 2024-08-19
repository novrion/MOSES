#!/bin/bash

cd ..

# define files
orderfile="data/in/order"
infile="data/in/in"

# clear input files
for line in $(cat $orderfile | tr -d " \r")
do
  > $line
done

# clear in file
> $infile
