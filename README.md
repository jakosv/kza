# Kolmogorov-Zurbenko Adaptive Filter 

## Build 

Makefile has multiple targets:
- prefix_sum - prefix sum optimization;
- kza_timer - time measurement;
- prefix_sum-timer - prefix sum optimization with time measurement.

```console
$ make prefix_sum-timer
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
