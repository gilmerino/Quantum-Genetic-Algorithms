# Quantum Genetic Algorithms in R
This repository includes the companion software to the Master Thesis "Building Quantum Genetic Algorithms: Basis and a First Implementation" for the Master in Research in Artificial Intelligence at the UIMP-AEPIA.

There are three pieces of code written in the R language. Two of them are quantum genetic algorithms (QGAs): a basic one (BQGA) that includes quantum mutation but not crossover, with a similar structure to the classical versions of GAs; a "reduced" one (RQGA) that it is based on the Grover's algorithm, including an oracle and the diffusion operator. The third piece of code is a simple classical genetic algorithm (SGA) to compare, side-by-side, to the BQGA. Bibliographic references on which the algorithms are based, and other details, are included in the Master Thesis that will be soon available through University repositories.

The problem to solve in all of them is to maximize the 2D function **f(x,y)=x²+y²** as a benchmark function. The fitness function can be modified and adapted to other 2D or 1D problems. Multi-dimensional functions could also be implemented with a small modification of the chromosome encoding.
