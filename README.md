SWIFTsimIO Examples
===================

This repository contains a number of example repositories for the
SWIFTsimIO code. They use SWIFT snapshots and where possible note
example simulations that they may be used with. They are catalogued
below:

+ `rho_T_movie`: Creates a rho-T movie (using FuncAnimation; this
  requires ffmpeg to be installed) or a simple plot. Comes with
  a version that reads GADGET/EAGLE data to produce the same
  plot but for older codes. There is a parallel version also available
  if you have `p_tqdm` installed.

+ `temperature_map`: Creates a temperature map of a given snapshot and
  saves it to disk.
