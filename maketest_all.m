%MAKETEST_ALL Make and runs tests for all of stenglib/.

% S. Engblom 2010-09-23

cd Tensor
startup
make
spin

cd ../Fast
startup
make
spin
scrub

cd ../Scicomp
startup
spin

cd ../Utils
startup
spin
scrub

cd ..
