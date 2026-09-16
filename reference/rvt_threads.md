# Number of threads used by the C++ kernels

Defaults to all detected cores. Set with `rvt_threads(n)` or by setting
`options(rvtr.threads = n)`.

## Usage

``` r
rvt_threads(n = NULL)
```

## Arguments

- n:

  number of threads, or NULL to just query the current setting

## Value

the number of threads in effect

## Examples

``` r
rvt_threads()
#> [1] 4
```
