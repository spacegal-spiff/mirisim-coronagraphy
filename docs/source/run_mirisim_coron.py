# import the configuration file parsers so they can be written to file
from mirisim.config_parser import SimConfig, SimulatorConfig, SceneConfig

# yaml parsing to get the overall parameters for the observation and target
import sys
import yaml

# import scene component generators and simulation tools
from mirisim.skysim import Background, sed, Point
from mirisim.skysim import wrap_pysynphot as wS
from mirisim.skysim import Skycube
from mirisim import MiriSimulation

# more general imports
import numpy as np
import glob                 
import os                   
from astropy.io import fits 
import time
import matplotlib.pyplot as plt    
from matplotlib import colors,cm
from astropy.visualization import ZScaleInterval

# custom tools + webbpsf_ext
from mirisim.imsim.webbPsf_imgtools import rotate
from mirisim.imsim.webbPsf_imgtools import convolve_extended_model
from mirisim.imsim.webbPsf_imgtools import compute_roll_ref_from_file

"""
TODO Need to add extended docstring

command line arguments:
argument 1: yaml file 
argument 2: whether to simulate science 'sci', else None
argument 3: whether to simulate PSF reference 'psf' else None

usage:
# generate both science images and PSF reference images
run_mirisim.py config.yaml sci psf 

# generate only science images 
run_mirisim.py config.yaml sci None 

# generate only PSF reference images
run_mirisim.py config.yaml None psf 


"""


## ---- observation setup ---- ##

yaml_file = sys.argv[1]

# import yaml simulation setup information
config = yaml.load(open(yaml_file), Loader=yaml.Loader)

# set up prefix for simulation to write out filenames
prefix = config['systemparams']['prefix']

# required mirisim coronagraphy dictionaries
cfigpath = {'F1065C':'IMA_CORO1065', 'F1140C':'IMA_CORO1140', 'F1550C':'IMA_CORO1550', 'F2300C':'IMA_COROLYOT'}
detector = {'F1065C': 'MASK1065', 'F1140C': 'MASK1140', 'F1550C': 'MASK1550', 'F2300C':'MASKLYOT'}
coronmsk = {'F1065C': '4QPM_1065', 'F1140C': '4QPM_1140', 'F1550C': '4QPM_1550', 'F2300C':'LYOT_2300'}

# set filtername for all observations (sci + PSF ref)
filtername = config['sci_setup']['filtername']

# set up for science observations
n_exp = config['sci_setup']['n_exp']
n_int = config['sci_setup']['n_int']
n_groups = config['sci_setup']['n_groups']

# set up simulator configuration
simulator_config = SimulatorConfig(config['configfilepaths']['simulator_config'])

 


## -- scene generation -- # 
 
if sys.argv[2] == 'sci':
    
    t_init = time.time()
    
    # set up the background for the science observations
    bg = Background(level = config['sci_setup']['bg_level'], gradient = config['sci_setup']['bg_gradient'], pa = config['sci_setup']['bg_pa'])


            
    # make stars and planets as blackbodies - TODO update later
    def create_star(xpos,ypos):
        star = Point(Cen = (xpos,ypos))
        BBparams = {'Temp': config['starparams']['teff'], # 1e4 for Vega, 4170 for pi Her
                    'wref': config['starparams']['flux_ujy_refwv'], # was 10 um
                    'flux': float(config['starparams']['flux_ujy'])} # units are microjansky, 16.09e6 for Vega, 18.99e6 for pi Her
        Blackbody = sed.BBSed(**BBparams)
        star.set_SED(Blackbody)
        
        return star

    def create_planet(xpos,ypos):
        star = Point(Cen = (xpos,ypos))
        BBparams = {'Temp': config['companion1']['teff'], # 1e4 for Vega, 4170 for pi Her
                    'wref': config['companion1']['flux_ujy_refwv'], # was 10 um
                    'flux': float(config['companion1']['flux_ujy'])} # units are microjansky, 16.09e6 for Vega, 18.99e6 for pi Her
        Blackbody = sed.BBSed(**BBparams)
        star.set_SED(Blackbody)
        
        return star

    # do you want to convolve a disk?
    if config['extendedscene']['disk'].upper() == 'Y':
        
        # check to see if existing convolution can be used; 
        # otherwise re-do convolution (takes more time)
        if config['extendedscene']['use_existing_convolution'].upper() == 'N':
            convolve_extended_model(config['extendedscene']['modelpath'], \
                filtername, dist_pc = config['systemparams']['dist'], pa = config['systemparams']['pa_deg'], pixscale = config['extendedscene']['model_px_scale'], model_wv = config['extendedscene']['model_wv'], model_units = config['extendedscene']['model_units'], central_star = config['extendedscene']['central_star'], model_path_output = config['extendedscene']['model_path_output'], simulator_config = simulator_config)

        if config['extendedscene']['central_star'] == False:   
            # assume star perfectly centered - TODO update for SGD and TA errors
            star = create_star(config['systemparams']['starposx'], config['systemparams']['starposx']) 

        # read in convolved disk cube:
        cube = Skycube(config['extendedscene']['model_path_output'], simulator_config=simulator_config)        
        
        targetlist = []
        targetlist.append(cube)
        if config['extendedscene']['central_star'] == False:
            targetlist.append(star)
        
        # generate scene with convolved model cube
        scene_config = SceneConfig.makeScene(loglevel=0,
                                        background=bg, targets=targetlist)

        os.system(f'rm {prefix}_scene.ini')
        scene_config.write(f'{prefix}_scene.ini')
        print(f"Saved scene config to {prefix}_scene.ini")        


    # no disk, only companion(s) 
    else:
        # assume star perfectly centered - TODO update for SGD and TA errors
        star = create_star(config['systemparams']['starposx'], config['systemparams']['starposx'])

        # read in planet position from yaml config
        planet_ra_offset = config['companion1']['ra_off']
        planet_dec_offset = config['companion1']['dec_off']

        # assume some sort of position angle offset for roll angle of observations... TODO CHECK ROTATION ANGLE, currently PA_V3
        planet_rotate = rotate([(planet_ra_offset, planet_dec_offset)], origin=(0,0), degrees=config['systemparams']['pa_deg'])

        # make that planet
        planet = create_planet(planet_rotate[0], planet_rotate[1])

        # TODO make more planets!

        # make the targetlist for the scene:
        targetlist = [star, planet]

        # make scene itself
        scene_config = SceneConfig.makeScene(loglevel=0, background=bg, targets=targetlist)

        # remove existing scene
        os.system(f'rm {prefix}_scene.ini')
        scene_config.write(f'{prefix}_scene.ini')
        print(f"saved scene config to {prefix}_scene.ini")


    ## ---- simulation setup ---- ##

    # set up simulaTION config:
    sim_config = SimConfig.makeSim(
        name = 'ima_simulation',    # name given to simulation
        scene = f'{prefix}_scene.ini', # name of scene file to input - IMPORTANT - this needs to match name of scene defined earlier
        rel_obsdate = 0.0,          # relative observation date (0 = launch, 1 = end of 5 yrs)
        POP = 'IMA',                # Component on which to center (Imager or MRS)
        ConfigPath = cfigpath[filtername],       # Configure the Optical path (MRS sub-band) - e.g., IMA_CORO1140
        Dither = False,             # Don't Dither
        StartInd = 1,               # start index for dither pattern [NOT USED HERE]
        NDither = 2,                # number of dither positions [NOT USED HERE]
        DitherPat = 'ima_recommended_dither.dat', # dither pattern to use [NOT USED HERE]
        disperser = 'SHORT',        # [NOT USED HERE]
        detector = 'SW',            # [NOT USED HERE]
        mrs_mode = 'SLOW',          # [NOT USED HERE]
        mrs_exposures = 1,          # [NOT USED HERE]
        mrs_integrations = 1,       # [NOT USED HERE]
        mrs_frames = 1,             # [NOT USED HERE]
        ima_exposures = n_exp,          # number of exposures
        ima_integrations = n_int,       # number of integrations
        ima_frames = n_groups,             # number of groups (for MIRI, # Groups = # Frames)
        ima_mode = 'FAST',          # Imager read mode (default is FAST ~ 2.3 s)
        filter = filtername,          # Imager Filter to use - F1140C
        readDetect = detector[filtername]         # Portion of detector to read out - e.g., MASK1550
    )

    # Save these to disk, overwriting any previous version
    os.system(f'rm {prefix}_simulation.ini')
    sim_config.write(f'{prefix}_simulation.ini')

    mysim = MiriSimulation(sim_config, scene_config, simulator_config)


    ## ---- run all the science stuff ---- ##

    t0 = time.time()
    mysim.run()
    t1 = time.time()
    print(f"Total time to run science simulation with PA_V3 roll angle {config['systemparams']['pa_deg']}: {np.round((t1-t0)/60.,3)} minutes")



    ## ---- housekeeping for the pipeline ---- ##

    outputdir = sorted(glob.glob('*_*_mirisim'),key=os.path.getmtime)[-1]    #[-1] takes the last entry found

    print("Output directory is...", outputdir)

    outputDirContents = os.listdir(outputdir)

    directories = [name for name in outputDirContents if os.path.isdir(os.path.join(outputdir,name))]
    files = [name for name in outputDirContents if not os.path.isdir(os.path.join(outputdir,name))]


    infits = glob.glob('{}/det_images/*.fits'.format(outputdir))

    # in case there are multiple exposures
    for idx, exposure in enumerate(infits):

        hdul = fits.open(exposure)

        # add in the fits header keywords that we need to run the pipeline...

        roll_ref = compute_roll_ref_from_file(exposure, config['systemparams']['pa_deg']) # assuming the latter is pa_v3

        hdul[0].header['CORONMSK'] = coronmsk[filtername]
        hdul[0].header['XOFFSET'] = config['systemparams']['starposx']
        hdul[0].header['YOFFSET'] = config['systemparams']['starposy']
        hdul[0].header['TARGPROP'] = config['starparams']['name']

        hdul[0].header['HISTORY'] = 'MIRISim Coron: Manually added keywords: CORONMSK, XOFFSET, YOFFSET, TARGPROP.'         

        hdul[1].header['ROLL_REF'] = roll_ref
        hdul[1].header['PA_V3'] = config['systemparams']['pa_deg']

        hdul[1].header['HISTORY'] = 'MIRISim Coron: Overwrote default keywords PA_V3, ROLL_REF.'    

        hdul.writeto(exposure, overwrite=True)
        print("Overwrote", exposure)

    # rename directory to be something a little more informative
    os.rename(outputdir, outputdir + f'_{prefix}')
    print("Renamed output directory to" + outputdir + f'_{prefix}')

    t_fin = time.time()

    print(f"Total running time including convolution: {np.round((t_fin-t_init)/60.,3)} minutes")

"""
PSF Reference Observations!
"""



## ---- OPTIONAL: run simulations for PSF reference ---- ##
if sys.argv[3] == 'psf':
    prefix = prefix + '_psfref'

    # set up for PSF observations
    n_exp = config['psf_setup']['n_exp']
    n_int = config['psf_setup']['n_int']
    n_groups = config['psf_setup']['n_groups']

    # generate set of SGD values for 5 point dither positions
    from webbpsf_ext.coords import get_idl_offset

    base_offset = get_idl_offset(base_std=None) # where base_std is the TA uncertainty, assumed to be 5.0 mas
    # dith_std is the one-sigma pointing uncertainty for dithers. It's 2.5 mas if SGD and 5.0 mas otherwise     
    # TODO add in functionality for 9 pt SGD pattern

    if config['psf_setup']['sgd_pattern'] == 5:

        dith0 = get_idl_offset(base_offset, dith_offset=(0,0), dith_std=None)
        dith1 = get_idl_offset(base_offset, dith_offset=(-0.01,+0.01), dith_std=None)
        dith2 = get_idl_offset(base_offset, dith_offset=(+0.01,+0.01), dith_std=None)
        dith3 = get_idl_offset(base_offset, dith_offset=(+0.01,-0.01), dith_std=None)
        dith4 = get_idl_offset(base_offset, dith_offset=(-0.01,-0.01), dith_std=None)

        sgd_list = [dith0, dith1, dith2, dith3, dith4]


    # set up the background for the observations
    bg = Background(level = config['sci_setup']['bg_level'], gradient = config['sci_setup']['bg_gradient'], pa = config['sci_setup']['bg_pa'])

    # Make a single star offset from the center by some amount - in this case, centered behind 4QPM
    # Create point source object within the Imager Field of View, using MIRISim built-in specifications
    # TODO: Use a more sophisticated spectrum using webbpsf_ext tools
    def create_star(xpos,ypos):
        star = Point(Cen = (xpos,ypos))
        BBparams = {'Temp': config['psfrefparams']['teff'], # 1e4 for Vega, 4170 for pi Her
                    'wref': config['psfrefparams']['flux_ujy_refwv'], # was 10 um
                    'flux': float(config['psfrefparams']['flux_ujy'])} # units are microjansky, 16.09e6 for Vega, 18.99e6 for pi Her
        Blackbody = sed.BBSed(**BBparams)
        star.set_SED(Blackbody)
        
        return star   

    # first, save the values so we know which were randomly generated for this run
    np.savetxt(f'{prefix}_SGD_values.txt', np.array(sgd_list))

    # time total run
    t_all_start = time.time()

    for idx, dith in enumerate(sgd_list):
        print(f"Simulating small grid dither offset at {dith}")
        
        star = create_star(dith[0], dith[1]) # x, y in arcsec, following RA convention (negative to the west, detector right)

        targetlist = [star] # only central star

        # Now we can turn this all into a configuration file (scene.ini)
        scene_config = SceneConfig.makeScene(loglevel=0,
                                            background=bg, targets=targetlist)

        os.system(f'rm {prefix}_scene.ini')
        scene_config.write(f'{prefix}_scene.ini')
        print(f"saved scene config to {prefix}_scene.ini")
        
        sim_config = SimConfig.makeSim(
            name = 'ima_simulation',    # name given to simulation
            scene = f'{prefix}_scene.ini', # name of scene file to input - IMPORTANT - this needs to match name of scene defined earlier
            rel_obsdate = 0.0,          # relative observation date (0 = launch, 1 = end of 5 yrs)
            POP = 'IMA',                # Component on which to center (Imager or MRS)
            ConfigPath = cfigpath[filtername],       # Configure the Optical path (MRS sub-band) - e.g., IMA_CORO1140
            Dither = False,             # Don't Dither
            StartInd = 1,               # start index for dither pattern [NOT USED HERE]
            NDither = 2,                # number of dither positions [NOT USED HERE]
            DitherPat = 'ima_recommended_dither.dat', # dither pattern to use [NOT USED HERE]
            disperser = 'SHORT',        # [NOT USED HERE]
            detector = 'SW',            # [NOT USED HERE]
            mrs_mode = 'SLOW',          # [NOT USED HERE]
            mrs_exposures = 1,          # [NOT USED HERE]
            mrs_integrations = 1,       # [NOT USED HERE]
            mrs_frames = 1,             # [NOT USED HERE]
            ima_exposures = n_exp,          # number of exposures
            ima_integrations = n_int,       # number of integrations
            ima_frames = n_groups,             # number of groups (for MIRI, # Groups = # Frames)
            ima_mode = 'FAST',          # Imager read mode (default is FAST ~ 2.3 s)
            filter = filtername,          # Imager Filter to use - F1140C
            readDetect = detector[filtername]         # Portion of detector to read out - e.g., MASK1550
        )
        
        # Save these to disk, overwriting any previous version
        os.system(f'rm {prefix}_simulation.ini')
        sim_config.write(f'{prefix}_simulation.ini')
        
        mysim = MiriSimulation(sim_config, scene_config, simulator_config)
        
        # time each of the runs
        t0 = time.time()
        mysim.run()
        t1 = time.time()
        
        print(f"Time taken to run sgd: {np.round((t1-t0)/60.,2)} minutes")
        
        # add in the fits header keywords that we need to run the pipeline... 
        outputdir = sorted(glob.glob('*_*_mirisim'),key=os.path.getmtime)[-1]    #[-1] takes the last entry found

        print("Output directory is...", outputdir)

        outputDirContents = os.listdir(outputdir)

        directories = [name for name in outputDirContents if os.path.isdir(os.path.join(outputdir,name))]
        files = [name for name in outputDirContents if not os.path.isdir(os.path.join(outputdir,name))]

        #print('The subdirectories in the outputdirectory are:\n{}'.format(directories))
        #print('The files in the outputdirectory are:\n{}'.format(files))
        
        infits = glob.glob('{}/det_images/*.fits'.format(outputdir))

        # assuming that the SGD index is the same as the order of the simulated files... I think (hope)
        dith = sgd_list[idx]        
        
        # in case there are multiple exposures
        for exposure in infits:

            hdul = fits.open(exposure)

            # add in the fits header keywords that we need to run the pipeline...

            roll_ref = compute_roll_ref_from_file(exposure, config['systemparams']['pa_deg']) # assuming the latter is pa_v3

            hdul[0].header['CORONMSK'] = coronmsk[filtername]
            hdul[0].header['XOFFSET'] = config['systemparams']['starposx']
            hdul[0].header['YOFFSET'] = config['systemparams']['starposy']
            hdul[0].header['TARGPROP'] = config['starparams']['name']

            hdul[0].header['HISTORY'] = 'MIRISim Coron: Manually added keywords: CORONMSK, XOFFSET, YOFFSET, TARGPROP.'         

            hdul[1].header['ROLL_REF'] = roll_ref
            hdul[1].header['PA_V3'] = config['systemparams']['pa_deg']

            hdul[1].header['HISTORY'] = 'MIRISim Coron: Overwrote default keywords PA_V3, ROLL_REF.'    

            hdul.writeto(exposure, overwrite=True)
            print("Overwrote", exposure)
        
        #rename directory to be something a little more informative
        os.rename(outputdir, outputdir + f'_{prefix}_SGD{idx}')
        print("Renamed output directory to" + outputdir + f'_{prefix}_SGD{idx}')

        
    t_all_end = time.time()

    print(f"Time taken to run full set of sgd: {np.round((t_all_end-t_all_start)/60.,2)} minutes")