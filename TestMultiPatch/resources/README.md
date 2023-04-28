# Pluto.jl notebook

The `Pluto.jl` notebook [`coordinate_tests.jl`](./resources/coordinate_tests.jl) compares data output by `CarpetX` via the [`run_cake_tests.par`](../../par/run_cake_tests.par) parameter file with the expected output of `global2local` and `local2global`. If the tests pass, no output is produced, otherwise an error message will be printed. Note that in order to reproduce the tests, one must change the directory information near the top of the notebook 
