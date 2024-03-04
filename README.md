# perigrain-v2

Granular media with peridynamics

# Installation along with dependencies

- For linux (tested on Ubuntu 20.04)
```bash
./install_linux.sh
```

- For- Mac
```bash
./install_mac.sh
```

# Dependencies

* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) (Specify the location of the `Eigen` library path in `makefile` if using other versions)
* `gmsh` to generate mesh
* To read `hdf5` files in linux using `h5dump`
* Python packages
    * `matplotlib`, `numpy` as usual
    * mesh related tools `pygmsh`, `gmsh`, `meshio`
    * optimization using `gekko`
    * parallel processing using `multiprocessing`, `pathos`

# Examples

Example simulations are located in `examples/` directory. Call the corresponding `run.sh` files from the root like this:

```
examples/3d-coll-sph/run.sh
```

By default, the outputs are stored in `examples_output/` directory.
