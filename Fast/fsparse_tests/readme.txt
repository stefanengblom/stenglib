FSPARSE
S. Engblom, D. Lukarski 2014-06-02

For a license statement and additional code, see:
http://user.it.uu.se/~stefane/freeware

-Run make.m to compile. You may have to edit this file depending on
your platform.

-Type 'help fsparse' to see how to use FSPARSE.

- To reproduce the results from the paper you need to run:
1) serialmatlab.m compare our serial fsparse and MATLABs sparse
2) parallelmatlab.m compare our parallel fsparse and MATLABs sparse
3) percserial.m shows the percentage of each section in the serial
fsparse code
4) percparallel.m shows the percentage of each section in the parallel
fsparse code
5) perc.m shows the speed-up of each section in fsparse
(serial/parallel)

Note: All corresponding plots are created in .eps format.
