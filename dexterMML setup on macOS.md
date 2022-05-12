# dexterMML setup on macOS

## Install *gcc* using [Homebrew](https://brew.sh/)

```sh
brew install gcc
```
## Specify compiler options for R

Create the file `~/.R/Makevars` and copy in the contents depending on your Mac's processor platform:

### Apple Silicon based Mac

```sh
SHLIB_OPENMP_CFLAGS = -fopenmp
SHLIB_OPENMP_CXXFLAGS = -fopenmp
GCC_LOC = /opt/homebrew
#CC = $(GCC_LOC)/bin/gcc-11
#CXX = $(GCC_LOC)/bin/g++-11
CXX11 = $(GCC_LOC)/bin/g++-11
FLIBS = -L$(GCC_LOC)/lib/gcc/11 -lgfortran -lm
LDFLAGS = -L$(GCC_LOC)/lib -Wl,-rpath,$(GCC_LOC)/lib
```

### Intel based Mac

```sh
SHLIB_OPENMP_CFLAGS = -fopenmp
SHLIB_OPENMP_CXXFLAGS = -fopenmp
GCC_LOC = /usr/local/opt/gcc
#CC = $(GCC_LOC)/bin/gcc-11
#CXX = $(GCC_LOC)/bin/g++-11
CXX11 = $(GCC_LOC)/bin/g++-11
FLIBS = -L$(GCC_LOC)/lib/gcc/11 -lgfortran -lm
LDFLAGS = -L$(GCC_LOC)/lib -Wl,-rpath,$(GCC_LOC)/lib
```

## Install R and RStudio if you haven't already

```sh
brew install --cask r
brew install --cask rstudio
```

## Install dexterMML in R/RStudio

```R
install.packages("devtools")
devtools::install_github("dexter-psychometrics/dexterMML")
```

## Notes

The default compiler options for R can be found in `/Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/etc/Makeconf`.
