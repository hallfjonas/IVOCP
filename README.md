# IVOCP
An initial value optimal control problem with complementarity constraints.
This benchmark was used for demonstration [in this paper](https://arxiv.org/abs/2103.05965).

## Prerequisits
This is a simple MATLAB benchmark created under Ubuntu with the following prerequisits:

- Install CasADi for MATLAB (https://web.casadi.org/).
- Install the qpOASES 3.2.0 MATLAB interface (https://github.com/coin-or/qpOASES).
- Adjust the casadi path in the file `RunBenchmark.m`.
- Adjust the qpOASES MATLAB interface path in the file `RunBenchmark.m`.

## CODE Outline
- A MATLAB implementation of the LCQP algorithm is given in the directory `LCQP`.
- The script `RunBenchmark.m` can be used to run the benchmark.
- Afterwards, the function `CreatePerformancePlots.m` can be used to create performance plots regarding time and solution quality.
