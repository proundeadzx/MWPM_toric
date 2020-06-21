# MWPM_toric
This is a little program that solves defects on the toric code, using the blossom algorithm, the program will solve random defects for a number of generations for a number of error probabilities. The succes rate, p_s, for all error probabilities will be outputed to a txt file in a results directory.

##Prerequisites
After cloning the repo Vladimir Kolmogorov's Blossom V implementation needs to be downloaded from: http://pub.ist.ac.at/~vnk/software/blossom5-v2.05.src.tar.gz.

Then untar the downloaded tar-archive into the repo folder, then the program should be good to go!

## How to use
At the present the program is a bit barebones, the only way to change arguments for the program is at the present to change the values of variables in the program.

### p_errors
p_errors is the probability of an error appearing on a qubit, this variable takes a list of values that the program will use in order.

### system_size
system_size is the width and heigth of the toric code, takes an integer. The integer MUST be odd

### nbr_of_iterations
nbr_of_iterations is the number of syndromes the program will try to solve, takes an integer.

Below is an example of a configuration for error probability 0.05, 0.06 and 0.07, system size 5 and solving 1 000 000 syndromes

p_errors = [0.05, 0.06, 0.07]
system_size = 5
nbr_of_iterations = int(1e6)
