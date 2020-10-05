# README

Summary: Import Comsol as FinEtools model, solve the free vibration problem, export  the solution. 

## How to use

1. Download the zip file of the repository.
Expand it, and change the working folder to `FinEtoolsMultimaterialVibEP`.

2. Run Julia.

3. Activate and instantiate the package:
```
using Pkg; Pkg.activate("."); Pkg.instantiate()
```

4. Test that the repository works:
```
using Pkg; Pkg.test()
```

## Creating a new simulation

Copy `"twoblocks.json"` to `"mynewexample.json"`. Change the settings in the file to suit your particular situation. Note that the geometry is expected in a Nastran file (free format, only linear tetrahedra).

Now run
```
using FinEtoolsMultimaterialVibEP: solve_ep
solve_ep("mynewexample.json")
```

If everything went well, there will be a Matlab file `"mynewexample.mat"` with the data in the current directory.


