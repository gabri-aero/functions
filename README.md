# Functions

![Version](https://img.shields.io/badge/version-0.0.1-blue.svg)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/gabri-aero/fft/blob/main/LICENSE)
[![docs](https://img.shields.io/badge/Doxygen-Documentation-5A7BA7?logo=doxygen&logoColor=white&style=flat)](https://gabri-aero.github.io/functions/)

> A static library of functions used in astrodynamics.

## Features
- Fully-normalized Associated Legendre functions through standard forward column recursive approach (Holmes & Featherstone, 2002). First and second order derivatives are also supported.
- Fully-normalized inclination function computation through FFT (Wagner, 1983). First derivatives can also be computed similarly. Besides, cross-track inclination functions computation is provided (Balmino et al., 1996).
- Full-normalization constants (Heiskanen & Moritz, 1967).

## TO DO
- Eccentricity functions

## Setup

### Makefile

1. Clone this repository into your `external` folder
```sh
git clone https://github.com/gabri-aero/functions
```
2. Build the library
```sh
cd functions
make 
```
3. Link to your project accordingly
```sh
# Define path to functions repository
FUNCTIONS_DIR = $(EXTERNAL_DIR)/functions

# Define flags for compiler
FUNCTIONS_INCLUDE = -I$(FUNCTIONS_DIR)/include/ -I$(FUNCTIONS_DIR)  # include flags
FUNCTIONS_LIB_FLAGS = -L$(FUNCTIONS_DIR)/lib -lfunctions            # link flags
```

4. Include headers in your code
```cpp
// Include single headers
#include <functions/Flmp.hpp>
#include <functions/Plm.hpp>
#include <functions/Nlm.hpp>
// Alternatively all at once
#include <functions>
```

## References

Balmino, G., Schrama, E., & Sneeuw, N. (1996). Compatibility of first-order circular orbit perturbations theories; consequences for cross-track inclination functions. _Journal of Geodesy, 70_(9), 554–561. https://doi.org/10.1007/bf00867863

Heiskanen, W., & Moritz, H. (1967). _Physical Geodesy_. W. H. Freeman.  

Holmes, S. A., & Featherstone, W. E. (2002). A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions. _Journal of Geodesy, 76_(5), 279–299. https://doi.org/10.1007/s00190002-0216-2

Kaula, W. M. (1966). _Theory of Satellite Geodesy: Applications of Satellites to Geodesy._ Blaisdell Publishing Company.

Wagner, C. A. (1983). Direct determination of gravitational harmonics from low-low GRAVSAT data. _Journal of Geophysical Research: Solid Earth, 88_(B12), 10309–10321. https://doi.org/10.1029/jb088ib12p10309
