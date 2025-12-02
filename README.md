#080808
Inherited (and adapted) from github.com/ignabelitzky/gcas
#000000

# Runs
1-dimensional gas with discretized values for linear momenta $p$, with slight uncertainties both in the recipient length $L$ (1m ± 0.1mm) and in $p$,added when bumping on the walls: each $E_j$ is changed to $E_j'$ chosen from a uniform distribution within $[E_j-\Delta E/2, E_j+\Delta E/2]$, where $\Delta E = \alpha\left[(p - p_{\rm min}) \cdot (p_{\rm max}-p)\right]$.
Since the spatial distributions allow smooth ends, 2 $x$-channels are added beyond each wall.

### Instructions for program execution
This program is designed for compilation through the command "make".

### Program compilation
1. Open a terminal
2. Clone this repository
```Bash
git clone git@github.com:gcasgit/g1D-1000
```
3. Move inside the cloned directory
```Bash
cd g1D-1000
```
4. Run the following command to compile the program
```Bash
make
```
An executable (`main`) program will be created.

### Program execution
Once the program has been successfully compiled, it can be executed following these steps:
1. Make sure the input file `datos.in` is contained in the folder containing the executable program.
2. Run the `main` program:
```Bash
./main
```

**Note:** The files .dmp save dynamic info for $2^{21}$  particles.
