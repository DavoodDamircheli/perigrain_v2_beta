import argparse

def readargs():
    # Instantiate the parser
    parser = argparse.ArgumentParser(description='Optional app description')
    # Optional argument
    parser.add_argument('--shape', type=str, help='shape of grain', default='pertdisk')
    parser.add_argument('--seed', type=int, help='random see', default=1)


    ## wall
    parser.add_argument('--L', type=float, help='half-width of wall', default=3000e-3)
    parser.add_argument('--hh', type=float, help='half-height of wall', default=1200e-3)

    ## grain distribution
    parser.add_argument('--arrangement', type=str, help='type of particle arrangement', default='interior')
    # parser.add_argument('--arrangement', type=str, help='type of particle arrangement', default='grid')
    parser.add_argument('--P_meshsize', type=float, help='avg grain size when using incircles', default=160e-3)
    parser.add_argument('--manual_P_trim', action='store_true', help='specify lower limit of particle size manually')
    parser.add_argument('--P_trim_val', type=float, help='drop particle size below this', default=25e-3)

    parser.add_argument('--arr_grid_order', type=str, help='order of particles in grid arrangement', default='rand_unif')
    parser.add_argument('--min_rad', type=float, help='min radius of particle in arrangement', default=4e-3)
    parser.add_argument('--max_rad', type=float, help='min radius of particle in arrangement', default=20e-3)

    parser.add_argument('--delta', type=float, help='peridynamic horizon', default=45e-3)
    parser.add_argument('--contact_radius', type=float, help='contact radius', default=12e-3)

    ## wheel
    parser.add_argument('--wheel_rad', type=float, help='radius of wheel', default=0.355)
    parser.add_argument('--wheel_meshsize', type=float, help='wheel meshsize', default=12e-3)
    parser.add_argument('--wheel_rho_scale', type=float, help='wheel density multiplier', default=3)
    parser.add_argument('--wheel_rad_file', type=str, help='output filename containing wheel radius', default='output/wheel_rad')

    parser.add_argument('--fit_wheel', action='store_true', help='scale wheel to fit within the wall, for testing smaller setup')
    parser.add_argument('--fit_wheel_scale', type=float, help='fraction of wall half-height to set as wheel radius', default=3)

    ## grain properties
    parser.add_argument('--grain_meshsize', type=float, help='grain meshsize', default=12e-3)
    parser.add_argument('--grain_rho_scale', type=float, help='grain density multiplier', default=0.8)
    parser.add_argument('--grain_Gnot_scale', type=float, help='grain fracture toughness multiplier', default=1)
    parser.add_argument('--grain_G_scale', type=float, help='grain modulus multiplier', default=1)

    ## shape param
    parser.add_argument('--plus_notch_dist', type=float, help='shape param', default= 0.25)
    parser.add_argument('--ring_inner_ratio', type=float, help='shape param', default= 0.8)
    parser.add_argument('--ring_meshsize_factor', type=float, help='shape param', default= 2)
    parser.add_argument('--circ_steps', type=float, help='shape param', default= 16)
    parser.add_argument('--pertdisk_steps', type=float, help='shape param', default= 16)
    parser.add_argument('--pertdisk_incirc_ratio', type=float, help='shape param', default= 0.7)
    parser.add_argument('--n4_ratio_vol', type=float, help='shape param', default= 1)

    parser.add_argument('--setup_file', type=str, help='output setup directory', default='data/hdf5/all.h5')
    parser.add_argument('--pngfile', type=str, help='png file name with path for exp setup', default='setup.png')
    parser.add_argument('--size_distribution_pngfile', type=str, help='particle size distribution plot', default='size_distribution.png')
    parser.add_argument('--other_output_dir', type=str, help='output directory', default='data/hdf5')
    parser.add_argument('--dotsize', type=float, help='dotsize in plot', default= 0.1)
    parser.add_argument('--plot', action='store_true', help='whether plot or not')
    # finish parsing
    args = parser.parse_args()

    return args
