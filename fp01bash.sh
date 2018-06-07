#!/bin/bash
CPLEX_BIN_DIR=/opt/ibm/ILOG/CPLEX_Studio1271/cplex/bin/x86-64_linux
INPUT_DIR=/home/giorgiograni/Downloads/miplib2010-benchmark
JAR=/home/giorgiograni/IdeaProjects/FP01/out/artifacts/FP01_jar/FP01.jar
OUTPUT=/home/giorgiograni/IdeaProjects/FP01/results/statsTrial.csv
for d in $INPUT_DIR/*; do
		  fromname=${d##*/}
          echo  "Processing  $fromname"
          java   -Djava.library.path=$CPLEX_BIN_DIR -jar $JAR $d $OUTPUT $fromname  
done
