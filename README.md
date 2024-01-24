# Kolmogorov-Zurbenko Adaptive Filter 

## Build 

Makefile has multiple targets:
- prefix_sum - prefix sum optimization;
- threads_client_server - client-server multithreading implementation;
- threads_loop - loop multithreading implementation;
- prefix_sum-threads_client_server - client-server multithreading implementation with prefix sum optimization;
- prefix_sum-threads_loop - loop multithreading implementation with prefix sum
optimization;
- kza_timer - default target without any optimizations with time
measurement.

Also each target has a prefix -timer, which enables time measurement of
algorithms. 

```console
$ make prefix_sum-threads_client_server-debug
```

## Usage 
Makefile builds shared library. For example, we can use it in python.

```console
$ python kza_example.py
$ python kz_example.py 
```

## Testing 

```console
$ python run_tests.py
```
