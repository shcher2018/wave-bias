# wave-bias
This is a collection of MATLAB scripts for simulation of wave-induced biases in ADCP measurements from quasi-Lagrangian platforms

## Quick start
Running 
`
  run_waves_bias_simulation
`
should produce a simulation of wave-induced bias for a test case and plot it alongside the analytical expression. 

## Setting up your problem
Configuration of the semi-analytical model is achieved by populating a few structures before calling the main model routine (`simulate_wave_bias`):
*   WS - wave field parameters
*   PL - platform response parameters
*   ADCP - ADCP configuration

See `run_waves_bias_simulation.m` for an example.
