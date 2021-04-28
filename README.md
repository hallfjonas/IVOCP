# IVOCP
An initial value optimal control problem with complementarity constraints.
This benchmark was used for demonstration [in this paper](https://arxiv.org/abs/2103.05965).

## Prerequisits
This is a simple MATLAB benchmark created under Ubuntu with the following prerequisits:


### Obligatory
- Install CasADi for MATLAB (https://web.casadi.org/).
- Install qpOASES (https://github.com/coin-or/qpOASES).
- Adjust the casadi path in the file `RunBenchmark.m`.
- Adjust the qpOASES MATLAB interface path in the file `RunBenchmark.m`.

### Optional
If you want to make use of sparse solution variant within qpOASES you must compile qpOASES with the `SPARSE_SOLVER=MA57` flag in `make_linux.mk`.

## CODE Outline
- A MATLAB implementation of the LCQP algorithm is given in the directory `LCQP`.
- The script `RunBenchmark.m` can be used to run the benchmark.
- Afterwards, the script `CreatePerformancePlots.m` can be used to create performance plots regarding time and solution quality.
