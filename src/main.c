/*****************************************************************************
* 
* The following source code is in the public domain.
* Specifically, we give to the public domain all rights for future licensing
* of the source code, all resale rights, and all publishing rights.
* 
* We ask, but do not require, that the following message be included in all
* derived works:
* 
* Portions developed at the National Center for Supercomputing Applications at
* the University of Illinois at Urbana-Champaign.
* 
* THE UNIVERSITY OF ILLINOIS GIVES NO WARRANTY, EXPRESSED OR IMPLIED, FOR THE
* SOFTWARE AND/OR DOCUMENTATION PROVIDED, INCLUDING, WITHOUT LIMITATION,
* WARRANTY OF MERCHANTABILITY AND WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE
* 
****************************************************************************/

/*
 * Updated - option to now either display the polygons on the SGI and/or
 * write them out in either HDF VSet format, or Wavefront format. Modified by
 * Mike Krogh, NCSA, July 11, 1990 
 *
 */

/*
 * This program implements the marching cubes surface tiler described by
 * Lorensen & Cline in the Siggraph 87 Conference Proceedings. 
 *
 * This program gets its data from a 3D Scientific Data Set in an NCSA HDF file
 * and creates two files, one containing vertices and one containing
 * connectivity.  The user must specify the three file names and the
 * threshold value on the command line. 
 *
 * Written by Mike Krogh, NCSA, Feb.  2, 1990 
 *
 */

#include <stdio.h>

#define MAIN
#include "isovis.h"

main(argc, argv)
int argc;
char *argv[];
{

    float *data;		/* a pointer to the 3d volume of data */
    float threshold;		/* the threshold to look for */
    float max, min;		/* the maximum & minimum data values */
    int xdim, ydim, zdim;	/* the dimensions of the data */
    char inp_file[80];		/* name of input file */
#ifdef SGI
    char *image;		/* the raster image */
    int img_xdim, img_ydim;	/* the dimensions of the raster image */
#endif
    int parse_options();	/* parse the command line options */
    void usage();		/* print out the usage for this program */
    int get_hdf_data();		/* reads in hdf data */
    int get_raw_data();		/* reads in raw data */
    int iso_surface();		/* generates the iso-surface */
    int dump_wft();		/* writes out a wft file from stored polygons */
    int dump_vset();		/* writes out a vset from stored polygons */
    int dump_byu();		/* writes out a movie.byu file from polygons */
#ifdef SGI
    int init_graphics();	/* open the graphics window */
    int get_image();		/* captures the raster image on the screen */
    int write_image();		/* stores the raster image to an HDF file */
    void window_event_loop();	/* refreshes the screen, etc. */
#endif
    double atof();


    strcpy(MY_NAME, argv[0]);

#ifdef NCSA_VERBIAGE
    printf("\n");
    printf("    +-------------------------------------------------------------------+\n");
    printf("    |   isovis, version: 1.00, date: Aug. 14, 1990                      |\n");
    printf("    |   public domain software, NCSA, Urbana-Champaign, IL              |\n");
    printf("    +-------------------------------------------------------------------+\n");
    printf("\n");
#endif

    if (argc == 1) {
	usage();
	exit(1);
    }
    strcpy(inp_file, argv[(argc - 2)]);

    threshold = (float) atof(argv[(argc - 1)]);

    if (parse_options(argc, argv) == -1) {
	fprintf(stderr, "%s: error from parse_options\n", MY_NAME);
	exit(1);
    }
    if (VERBOSE)
	printf("%s: looking for threshold %f in %s\n", MY_NAME, threshold, inp_file);


    if (RAW_INPUT) {
	if (XDIM <= 0 || YDIM <= 0 || ZDIM <= 0) {
	    fprintf(stderr, "%s: bad dimensions for raw input file\n", MY_NAME);
	    exit(1);
	}
	xdim = XDIM; ydim = YDIM; zdim = ZDIM;
	if (get_raw_data(inp_file, &data, xdim, ydim, zdim, &max, &min) == -1) {
	    fprintf(stderr, "%s: error from get_raw_data\n", MY_NAME);
	    exit(1);
	}
    } else {
	if (get_hdf_data(inp_file, &data, &xdim, &ydim, &zdim, &max, &min) == -1) {
	    fprintf(stderr, "%s: error from get_hdf_data\n", MY_NAME);
	    exit(1);
	}
    }

    if (VERBOSE)
	printf("%s: minimum value %f maximum value %f\n", MY_NAME, min, max);

    if ((threshold < min) || (max < threshold)) {
	fprintf(stderr, "%s: minimum value %f maximum value %f\n", MY_NAME, min, max);
	fprintf(stderr, "%s: error, threshold must be between the maximum\n",
		MY_NAME);
	fprintf(stderr, "and minimum data values.\n");
	exit(1);
    }
#ifdef SGI
    if (DISPLAY)
	if (init_graphics(xdim, ydim, zdim) == -1) {
	    fprintf(stderr, "%s: error from init_graphics\n", MY_NAME);
	    exit(1);
	}
#endif

    if (strcmp(WFT_NAME, ""))
	/* If a wft name is specified, set the polygon store flag */
	STORE_POLYGONS = 1;

    if (strcmp(VSET_NAME, ""))
	/* If a VSet name is specified, set the polygon store flag */
	STORE_POLYGONS = 1;

    if (strcmp(BYU_NAME, ""))
	/* If a BYU name is specified, set the polygon store flag */
	STORE_POLYGONS = 1;

    if (DTM_OUTPUT)
	/* If a VSet name is specified, set the polygon store flag */
	STORE_POLYGONS = 1;


    if (iso_surface(data, xdim, ydim, zdim, threshold) == -1) {
	printf("%s: error from iso_surface\n", MY_NAME);
	exit(1);
    }
    if (SMOOTH) {
	if (NORMAL_TYPE) {
	    fprintf(stderr, "%s: warning smoothing gradient normals\n", MY_NAME);
	}
	smooth_norms();
	NORMAL_TYPE = 1;
    }
    if (strcmp(WFT_NAME, "")) {
	/* If a wft name is specified, write out a wft object file */
	if (dump_wft() == -1) {
	    fprintf(stderr, "%s: error from dump_wft\n", MY_NAME);
	    exit(1);
	}
    }
    if (strcmp(VSET_NAME, "")) {
	/* If a VSet name is specified, write out a HDF VSet */
	if (dump_vset() == -1) {
	    fprintf(stderr, "%s: error from dump_vset\n", MY_NAME);
	    exit(1);
	}
    }
    if (strcmp(BYU_NAME, "")) {
	/* If a BYU name is specified, write out a Movie.BYU file */
	if (dump_byu() == -1) {
	    fprintf(stderr, "%s: error from dump_byu\n", MY_NAME);
	    exit(1);
	}
    }
    if (DTM_OUTPUT) {
	/* output to the NCSA DTM protocol is requested */
	if (dump_dtm() == -1) {
	    fprintf(stderr, "%s: error from dump_dtm\n", MY_NAME);
	    exit(1);
	}
    }
#ifdef SGI
    if (DISPLAY)
	dump_sgi();		/* draw the picture */
    if (OUTNAME[0] != '\0') {
	/* save the raster image */
	if (get_image(&image, &img_xdim, &img_ydim) == -1) {
	    fprintf(stderr, "%s: error from get_image\n", MY_NAME);
	    exit(1);
	}
	if (write_image(image, img_xdim, img_ydim) == -1) {
	    fprintf(stderr, "%s: error from write_image\n", MY_NAME);
	    exit(1);
	}
	exit(0);
    }
#ifdef SGI_IMAGE
    if (SOUTNAME[0] != '\0') {
	/* save the raster image */
	if (write_sgi_image(SOUTNAME) == -1) {
	    fprintf(stderr, "%s: error from write_sgi_image\n", MY_NAME);
	    exit(1);
	}
	exit(0);
    }
#endif
    if (VERBOSE)
	printf("Done.\n");
    if (DISPLAY)
	window_event_loop(data, xdim, ydim, zdim, threshold);

#endif

    exit(0);

}
