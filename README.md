# Functions

![Version](https://img.shields.io/badge/version-0.0.1-blue.svg)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/gabri-aero/fft/blob/main/LICENSE)
[![docs](https://img.shields.io/badge/Doxygen-Documentation-5A7BA7?logo=doxygen&logoColor=white&style=flat)](https://gabri-aero.github.io/functions/)

> A header-only library of functions used in astrodynamics.

## Features
- Associated Legendre functions through standard forward column recursive approach (Holmes & Featherstone, 2002). First and second order derivatives are also supported.
- Inclination function computation through FFT (Wagner, 1983). First derivatives can also be computed similarly. Besides, cross-track inclination functions computation is provided (Balmino et al., 1996).
- Full-normalization constants (Heiskanen & Moritz, 1967). Fully-normalized versions of the functions are also available.

## TO DO
- Eccentricity functions

## Setup

### Makefile

1. Clone [`functions`](https://github.com/gabri-aero/functions) and [`fft`](https://github.com/gabri-aero/fft) repository into your `external` folder
```sh
mkdir -p external
cd external
git clone https://github.com/gabri-aero/functions
git clone https://github.com/gabri-aero/fft
```
2. Link to your project accordingly
```sh
# Define path to external repositories
FUNCTIONS_DIR = $(EXTERNAL_DIR)/functions
FFT_DIR = $(EXTERNAL_DIR)/fft

# Define flags for compiler
INCLUDE_FLAGS = -I$(FUNCTIONS_DIR) -I$(FFT_DIR) # include flags
```

3. Include headers in your code
```cpp
// Include headers
#include <functions>
```

## Tests
The unit tests can be built and run as follows:
```sh
make test
```

## References

Balmino, G., Schrama, E., & Sneeuw, N. (1996). Compatibility of first-order circular orbit perturbations theories; consequences for cross-track inclination functions. _Journal of Geodesy, 70_(9), 554–561. https://doi.org/10.1007/bf00867863

Heiskanen, W., & Moritz, H. (1967). _Physical Geodesy_. W. H. Freeman.  

Holmes, S. A., & Featherstone, W. E. (2002). A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions. _Journal of Geodesy, 76_(5), 279–299. https://doi.org/10.1007/s00190002-0216-2

Kaula, W. M. (1966). _Theory of Satellite Geodesy: Applications of Satellites to Geodesy._ Blaisdell Publishing Company.

Wagner, C. A. (1983). Direct determination of gravitational harmonics from low-low GRAVSAT data. _Journal of Geophysical Research: Solid Earth, 88_(B12), 10309–10321. https://doi.org/10.1029/jb088ib12p10309
