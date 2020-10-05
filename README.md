# Agent based models for simulating infectious diseases

Here are agent based models for simulating infectious diseases. The two main
programs are:

- abm_spaces.cc - a structured agent based model of a Covid-19 epidemic
- abm_hetgen.cc - an agent based model with highly heterogeneous agents to simulate
  short duration infection outbreaks

Details for compiling and running the programs are in the comments section at
the top of each programme. We develop and test under Linux.

Earlier prototypes of both programs are c19p5.R and hetgen.py. Please note
though that the C++ code has changed substantially since the prototypes which
are buggy and seldom if ever maintained. The speed difference between the C++
programs and the prototypes is so massive that we are now committed solely to
maintaining the C++ code.

Other useful files:

- Makefile - contains various useful targets

- analyze5.R, analyzeTaT2.R, analyze10000five.R, k_neighbours_average.R - for
  analysing the output of abm_spaces

- analyze_hetgen.py - for running abm_hetgen or analysing its output or both

- out_hetgen.csv.gz - output of abm_hetgen used in the paper doc_hetgen.tex

- spaces_data.tar.gz - output used in the paper doc_spaces.tex 

The analysis in the ABM-Spaces paper (doc_spaces.tex) can be reproduced by running:

- make analyze10000five (This is the primary analysis which will output five.htm
  and combine all the raw data in abm_spaces_output.csv)

- make analyzeTaT2 (This is the secondary analysis which will output TaT2.htm)
