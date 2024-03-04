for shape in circ n4 pertdisk plus
#for shape in n4
do
    ## small
    #python3 examples/wheel/setup.py  --L 500e-3 --hh 300e-3 --shape $shape --max_rad 20e-3 --min_rad 10e-3 --arr_grid_order rand_unif --dotsize 1 --pngfile setup_$shape.png

    #python3 examples/wheel/setup.py  --L 500e-3 --hh 300e-3 --shape $shape --max_rad 20e-3 --min_rad 10e-3 --arr_grid_order rand_unif --dotsize 0.1 --pngfile setup_interior_$shape.png --wheel_meshsize 10e-3 --grain_meshsize 3e-3 --arrangement interior

    ## big
    #python3 examples/wheel/setup.py  --shape $shape --dotsize 0.1 --pngfile setup_big_$shape.png --wheel_meshsize 8e-3 --grain_meshsize 10e-3 --arrangement interior --P_meshsize 170e-3

    # grid
    python3 examples/wheel/setup.py  --shape $shape --dotsize 0.1 --pngfile setup_big_grid_$shape.png --wheel_meshsize 20e-3 --grain_meshsize 10e-3 --arrangement grid --P_meshsize 170e-3 --min_rad 40e-3 --max_rad 45e-3
done
