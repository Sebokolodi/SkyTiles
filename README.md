# SkyTiles

#### Purpose : Extract tiles IDs for specific SBs and generate tile images.
    NB: The scripts are currently suitable for ASKAP fits data. 

#### Definitions of Scripts: 

	PyMapSkyTiles.py : Used for extracting tile locations (coordinates in degrees and pixels) 
                       and tile HPX IDs for a specific SB, 
                       and outputs csv files needed for creating tiles.  
			   It outputs 3 csv files. One csv contains the coordinates  
                       CRVAL1/2 and CRPIX1/2 values of the tiles (this actually the file
		           used for tilling, others are used for other purposes), 
                       the second file contains HPX IDs of tiles that are 
                       completely filled by one SB (we refer to these as 'SINGLE'),
			   and the last file contains tiles IDs requiring multiple 
                       SBs to be complete. 

	PyTiling.py      : The actual script for generating tile images. This script creates 
                        tiles using Montage. After several test, we wouldn't recommend using 
			    Montage for radio images, instead use CASA tiling (casa_tiling.py).      

    casa_tiling.py   : The actual script for generating tile images. This script uses CASA
	                  (previously known as CASAPY) to generate the tiles. CASA function 
		          called imregrid is used for tiling (interpolation='cubic'). CASA
	                  works well with our radio data.

    renaming_tiles.py: Script used to rename tiles fits data to a format agreed upon by the team:
	                   prefix_cenfreq_resol_RADEC_tileID_Stokes_version.fits
	
	Config_SkyMapTiles.json: Input to PyMapSkyTiles.py. This file contains the configurations 
                             needed for running PyMapSkyTiles.py, such as directory,
				 path to footprints, primary beam size, etc. See below for details. 
	
	Config_Tiling.json : Input to PyTiling.py. Works the same as the above json file but this is
                         is specifically for configuring PyTilling.py.   

#### The scripts are executed as follows:

	#./PyMapSkyTiles.py -j Config_SkyMapTiles.json
							
	#./PyTiling.py -j Config_Tiling.json	

    # casa -c ./casa_tiling.py  -h  (option 1: if you installed CASA all inclusive version) 
    
    #./casa_tiling.py -h (option 2: if you installed CASA modular version) 
  
    #./rename_tiles.py -i image.fits -pf PSM_pilot1 -v v2 -mfs
							

#### Config_SkyMapTiles.json Input Definitions and PyMapSkyTiles.py outputs:

    path_footprints   : Path to footprints.
    run_all_footprints: true/false. If true, look inside path_footprints and extract tile 
                        information for all footprints inside this path/directory. 
    footprint_files   : If run_all_footprints is false, then specify a specific footprint file(s). 
                        If run_all_footprints is false and no footprint_files is given, the code will
			    break.
    HPX_nside         : Healpix NSIDE.  
    tile_naxis        : Tile size (NAXIS of a tile. A value provided will be used for both naxis1 
                        and 2).
    tile_cdelt        : Tile pixel size (in degrees).
    beam_radius       : Size of the observation's primary beam (the radius, in degrees).
    beam_sample_points: Sampling the circumference of the beam. Used for determining pixel 
                        locations (ensures that no tile ID is missed). 
    number_of_beams   : The number of beams per SB. Most relevant for phased arrays. 
                        For POSSUM, we have 36 beams.
    outfile_prefix    : The prefix name to use for the output csvs files (include a path also, 
                        if not the files will be stored in your current directory).
    generate_ds9regions: if true, generate a DS9-suitable region file for the SB. 
     
    
    Output of executing PyMapSkyTiles.py: 
    i) csv file for each SB containing tile IDs, CRPIX1/2 and CRVAL1/2.The output name takes 
		the form: 'output_prefix_TileConfig_SBID.csv'
    ii) csv file containing tile IDs of the tiles completed by a single SB.  The output 
		name takes the form: 'output_prefix_SINGLE.csv'. Note that this file is only derived
            when multiple SBs (footprints) are provided via run_all_footprints or footprints_files. 
	        The file will be empty if only one SB is evaluated or if the SBs do not overlap. 
    iii) csv file containing tile IDs of the tiles needing more than a single SB to be complete. 
            The output name takes the form: 'output_prefix_REPEAT.csv'. This file will 
	        only be generated when multiple SBs are provided. See (ii).

     NB: Note that a file containing footprints (beam centres) is generated separately using ASKAP 
     tools called aces-aps. You may need an account to install this package. The tool repository can
     be found at  ssh://git@bitbucket.csiro.au:7999/aces/aces-apps.git and the documentation here: 
     https://bitbucket.csiro.au/projects/ACES/repos/aces-apps/browse. For assitance in installing this, 
     contact ASKAP directory (or ask POSSUM members involved with ASKAP software).

     Once you have it installed, you can simply run it using for example:
     # ./footprint-plan.py -n ak:square_6x6 -p 0.90 -a 0 -r 12.5,-52.5 -o possum_full_survery -f ds9 -w 0.90

     where -n is Name of footprint, -p is pitch angle, -a is Position angle of footprint (degrees, you need to 
     add 45 degrees to the rotation value), -r is observation RA and Dec, -o output filename, 
     -f format (e.g. ds9), -w Beam FWHM used in overlay files (degrees).  The information about observation
     RA, Dec, PA and name, can we found in a csv file generated for observation scheduling (see TG1 leads).
  
									  

#### Config_Tiling.json Input Definitions:

     input_image: Information about the input image. Insert path and name (below).
         path: 
         name: 
		 
    input_tile_config: CSV file containing tile ID, CRPIX1/2 information. 
                       Specify the path and name to these files (below). Note: 
		       The input csv file must take the format: 'filename_SBID.csv'. 
		       We use this file to also extract SBID. 
        path:
        SB_tile:
   
    output_tile: Output tile images and name. 
        path : Path to save the tiles and headers.
        naxis: Tile size should be the same as tile_naxis in Config_SkyMapTiles.json.
        cdelt: Tile pixel sample. Should be the same as tile_cdelt in Config_SkyMapTiles.json.
        fits_prefix: The prefix to give the resulting tile images.         
			The output tile names take the form: 
		       'fits_prefix-SBID-TileID.fits.     
			The stokes parameter is obtained from the input fits header
                    with crval3/4: 1 is I, 2 is Q, 3 is U and 4 is V.
					 

#### CASA Tiling Input Definitions

    -h, --help   show this help message and exit
    -i OBS_ID    Observation ID (e.g. 2156-54. It will be used on the filename. 
                  You can use anything e.g. SB ID such as 10040).
    -c CUBE      Image cube (or an MFS image) you want to tile. 
    -m MAP       Tiling map for the image cube [csv]. This is an csv contains the coordinates  
                 CRVAL1/2 and CRPIX1/2 values of the tiles.Example is provided, see PoSSUM-TileConfig-2156-54.csv.
    -o OUTPUT    Output write directory for tiles cubes/images.
    -t TEMPLATE  The template fits file. This is required (and is provided see tile-template.fits). This contains 
                 the header for a specific tile setting. In this case, it is for naxis=2048, NSIDE of 32, pixel
		     size of 3.2''. If a different tile settings are required, a user need to generate these themselves,
                 either using PyMapSkyTile and PyTiling.
    -n NAXIS     tile naxis. In this case 2048. If you change this, you may need to change csv file coordinates using
                 PyMapSkyTiles, and template fits using PyTiling.
    -p PREFIX    Prefix for output tile filenames

#### Rename_tiles Input Definitions:

    -i:   user must provide path and name of input image to rename.
    -pf:  user must provide path and prefix name to use for the output image. 
          If no path is provided, the fits file will be saved in the current directory.
    -v:   You can specify the version number of the output e.g v1 for version 1.
    -mfs: If file is MFS, you need to use this paramter indicate this. This ensures that
           the extension 't0' is added to the file name, and that the central frequency is 
           calculated accurately.

