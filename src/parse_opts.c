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

#include <stdio.h>
#include "isovis.h"

#ifdef DTM
/* This is for the NCSA DTM communications protocol */
#include <dtm.h>
#endif


/**************************** parse_options ****************************/
/**************************** parse_options ****************************/
/**************************** parse_options ****************************/
/**************************** parse_options ****************************/

int parse_options(argc,argv)
int argc;
char *argv[];
{

  int i;
  int found;
  int atoi();
  double atof();
  void defaults();


  if (!strcmp(argv[1],"-def")) {
     defaults();
     exit(1);
  }

  i=1;
  while (i<(argc-2)) {
    found=0;

    if (!strcmp(argv[i],"-r")) {
       if (i >= (argc-5)) {
          fprintf(stderr,"%s: not enough arguments for -r\n",MY_NAME);
          return -1;
       }
       i++;
       found=1;
#ifdef SGI
       ROTATION[0] = (float)atof(argv[i++]);
       ROTATION[1] = (float)atof(argv[i++]);
       ROTATION[2] = (float)atof(argv[i++]);
       if ((ROTATION[0] < -360.0) || (360.0 < ROTATION[0])) {
          fprintf(stderr,"%s: error, x rotation out of bounds\n",
                  MY_NAME);
          return -1;
       }
       if ((ROTATION[1] < -360.0) || (360.0 < ROTATION[1])) {
          fprintf(stderr,"%s: error, y rotation out of bounds\n",
                  MY_NAME);
          return -1;
       }
       if ((ROTATION[2] < -360.0) || (360.0 < ROTATION[2])) {
          fprintf(stderr,"%s: error, z rotation out of bounds\n",
                  MY_NAME);
          return -1;
       }
#else
       fprintf(stderr,"%s: warning, '-r' not a valid option for non-SGI\n",
               MY_NAME);
       i+=3;  /* skip over rotational values */
#endif
    }


    if (!strcmp(argv[i],"-t")) {
       if (i >= (argc-5)) {
          fprintf(stderr,"%s: not enough arguments for -t\n",MY_NAME);
          return -1;
       }
       i++;
       found=1;
#ifdef SGI
       TRANSLATION[0] = (float)atof(argv[i++]);
       TRANSLATION[1] = (float)atof(argv[i++]);
       TRANSLATION[2] = (float)atof(argv[i++]);
#else
       fprintf(stderr,"%s: warning, '-t' not a valid option for non-SGI\n",
               MY_NAME);
       i+=3;  /* skip over translational values */
#endif
    }

    if (!strcmp(argv[i],"-s")) {
       if (i >= (argc-5)) {
          fprintf(stderr,"%s: not enough arguments for -s\n",MY_NAME);
          return -1;
       }
       i++;
       found=1;
#ifdef SGI
       SCALE[0] = (float)atof(argv[i++]);
       SCALE[1] = (float)atof(argv[i++]);
       SCALE[2] = (float)atof(argv[i++]);
#else
       fprintf(stderr,"%s: warning, '-s' not a valid option for non-SGI\n",
               MY_NAME);
       i+=3;  /* skip over scaling values */
#endif
    }

    if (!strcmp(argv[i],"-ll")) {
       if (i >= (argc-5)) {
          fprintf(stderr,"%s: not enough arguments for -ll\n",MY_NAME);
          return -1;
       }
       i++;
       found=1;
#ifdef SGI
       LIGHT_LOC[0] = (float)atof(argv[i++]);
       LIGHT_LOC[1] = (float)atof(argv[i++]);
       LIGHT_LOC[2] = (float)atof(argv[i++]);
#else
       fprintf(stderr,"%s: warning, '-ll' not a valid option for non-SGI\n",
               MY_NAME);
       i+=3;  /* skip over light location values */
#endif
    }

    if (!strcmp(argv[i],"-lc")) {
       if (i >= (argc-5)) {
          fprintf(stderr,"%s: not enough arguments for -lc\n",MY_NAME);
          return -1;
       }
       i++;
       found=1;
#ifdef SGI
       LIGHT_COLOR[0] = atof(argv[i++]);
       LIGHT_COLOR[1] = atof(argv[i++]);
       LIGHT_COLOR[2] = atof(argv[i++]);
       if ((LIGHT_COLOR[0] < 0.0) || (1.0 < LIGHT_COLOR[0])) {
          fprintf(stderr,"%s: error, illegal red light value\n",
                  MY_NAME);
          return -1;
       }
       if ((LIGHT_COLOR[1] < 0.0) || (1.0 < LIGHT_COLOR[1])) {
          fprintf(stderr,"%s: error, illegal green light value\n",
                  MY_NAME);
          return -1;
       }
       if ((LIGHT_COLOR[2] < 0.0) || (1.0 < LIGHT_COLOR[2])) {
          fprintf(stderr,"%s: error, illegal blue light value\n",
                  MY_NAME);
          return -1;
       }
#else
       fprintf(stderr,"%s: warning, '-lc' not a valid option for non-SGI\n",
               MY_NAME);
       i+=3;  /* skip over light color values */
#endif
    }

    if (!strcmp(argv[i],"-k")) {
       if (i >= (argc-5)) {
          fprintf(stderr,"%s: not enough arguments for -k\n",MY_NAME);
          return -1;
       }
       i++;
       found=1;
#ifdef SGI
       K_VALUES[0] = (float)atof(argv[i++]);
       K_VALUES[1] = (float)atof(argv[i++]);
       K_VALUES[2] = (float)atof(argv[i++]);
       if ((K_VALUES[0] < 0.0) || (1.0 < K_VALUES[0])) {
          fprintf(stderr,"%s: error, illegal Ka value\n", MY_NAME);
          return -1;
       }
       if ((K_VALUES[1] < 0.0) || (1.0 < K_VALUES[1])) {
          fprintf(stderr,"%s: error, illegal Kd value\n", MY_NAME);
          return -1;
       }
       if ((K_VALUES[2] < 0.0) || (1.0 < K_VALUES[2])) {
          fprintf(stderr,"%s: error, illegal Ks value\n", MY_NAME);
          return -1;
       }
#else
       fprintf(stderr,"%s: warning, '-k' not a valid option for non-SGI\n",
               MY_NAME);
       i+=3;  /* skip over material values */
#endif
    }

    if (!strcmp(argv[i],"-ns")) {
       if (i >= (argc-3)) {
          fprintf(stderr,"%s: not enough arguments for -ns\n",MY_NAME);
          return -1;
       }
       i++;
       found=1;
#ifdef SGI
       K_EXP = (float)atof(argv[i++]);
       if ((K_EXP < 0.0) || (1000.0 < K_EXP)) {
          fprintf(stderr,"%s: error, illegal Ns value\n",
                  MY_NAME);
          return -1;
       }
#else
       fprintf(stderr,"%s: warning, '-ns' not a valid option for non-SGI\n",
               MY_NAME);
       i++;  /* skip over specular value */
#endif
    }

    if (!strcmp(argv[i],"-ntsc")) {
       i++;
       found=1;
#ifdef SGI
       WIN_POS[0] = 0;
       WIN_POS[1] = 0;
       WIN_SIZE[0] = 640;
       WIN_SIZE[1] = 483;
#else
       fprintf(stderr,"%s: warning, '-ntsc' not a valid option for non-SGI\n",
               MY_NAME);
#endif
    }

    if (!strcmp(argv[i],"-v")) {
       if (i >= (argc-6)) {
          fprintf(stderr,"%s: not enough arguments for -v\n",MY_NAME);
          return -1;
       }
       i++;
       found=1;
#ifdef SGI
       WIN_POS[0] = atoi(argv[i++]);
       WIN_POS[1] = atoi(argv[i++]);
       WIN_SIZE[0] = atoi(argv[i++]);
       WIN_SIZE[1] = atoi(argv[i++]);
       if ((WIN_POS[0] < 0) || (1279 < WIN_POS[0])) {
          fprintf(stderr,"%s: error, illegal window position\n", MY_NAME);
          return -1;
       }
       if ((WIN_POS[1] < 0) || (1023 < WIN_POS[1])) {
          fprintf(stderr,"%s: error, illegal window position\n", MY_NAME);
          return -1;
       }
       if ((WIN_SIZE[0] < 0) || ((1279-WIN_POS[0]) < WIN_SIZE[0])) {
          fprintf(stderr,"%s: error, illegal window size\n", MY_NAME);
          return -1;
       }
       if ((WIN_SIZE[1] < 0) || ((1023-WIN_POS[1]) < WIN_SIZE[1])) {
          fprintf(stderr,"%s: error, illegal window size\n", MY_NAME);
          return -1;
       }
#else
       fprintf(stderr,"%s: warning, '-v' not a valid option for non-SGI\n",
               MY_NAME);
       i+=4;  /* skip over window values */
#endif
    }

    if (!strcmp(argv[i],"-bc")) {
       if (i >= (argc-5)) {
          fprintf(stderr,"%s: not enough arguments for -bc\n",MY_NAME);
          return -1;
       }
       i++;
       found=1;
#ifdef SGI
       BACKG_COLOR[0] = atof(argv[i++]);
       BACKG_COLOR[1] = atof(argv[i++]);
       BACKG_COLOR[2] = atof(argv[i++]);
       if ((BACKG_COLOR[0] < 0.0) || (1.0 < BACKG_COLOR[0])) {
          fprintf(stderr,"%s: error, illegal red background value\n", MY_NAME);
          return -1;
       }
       if ((BACKG_COLOR[1] < 0.0) || (1.0 < BACKG_COLOR[1])) {
          fprintf(stderr,"%s: error, illegal green background value\n", MY_NAME);
          return -1;
       }
       if ((BACKG_COLOR[2] < 0.0) || (1.0 < BACKG_COLOR[2])) {
          fprintf(stderr,"%s: error, illegal blue background value\n", MY_NAME);
          return -1;
       }
#else
       fprintf(stderr,"%s: warning, '-bc' not a valid option for non-SGI\n",
               MY_NAME);
       i+=3;  /* skip over color values */
#endif
    }

    if (!strcmp(argv[i],"-norm")) {
       if (i >= (argc-3)) {
          fprintf(stderr,"%s: not enough arguments for -norm\n",MY_NAME);
          return -1;
       }
       i++;
       found=1;
       NORMAL_TYPE = atoi(argv[i++]);
       if ((NORMAL_TYPE < 0) || (1 < NORMAL_TYPE)) {
          fprintf(stderr,"%s: error, illegal normal type\n", MY_NAME);
          return -1;
       }
    }

    if (!strcmp(argv[i],"-o")) {
       if (i >= (argc-3)) {
          fprintf(stderr,"%s: not enough arguments for -o\n",MY_NAME);
          return -1;
       }
       i++;
       found=1;
#ifdef SGI
       strcpy(OUTNAME,argv[i++]);
#else
       fprintf(stderr,"%s: warning, '-o' not a valid option for non-SGI\n",
               MY_NAME);
       i++;  /* skip over raster name */
#endif
    }

    if (!strcmp(argv[i],"-obj")) {
       if (i >= (argc-3)) {
          fprintf(stderr,"%s: not enough arguments for -obj\n",MY_NAME);
          return -1;
       }
       i++;
       found=1;
       strcpy(WFT_NAME,argv[i++]);
    }

    if (!strcmp(argv[i],"-vset")) {
       if (i >= (argc-3)) {
          fprintf(stderr,"%s: not enough arguments for -vset\n",MY_NAME);
          return -1;
       }
       i++;
       found=1;
       strcpy(VSET_NAME,argv[i++]);
    }

#ifdef DTM 
    if (!strcmp(argv[i],"-DTM")) {
       if (i >= (argc-3)) {
          fprintf(stderr,"%s: not enough arguments for -DTM\n",MY_NAME);
          return -1;
       }
       i++;
       found=1;
       if (DTMmakePort(argv[i++]) == DTMERROR) {
          fprintf(stderr,"%s: error with DTM port specification\n",
                  MY_NAME);
          return -1;
       }
       DTM_OUTPUT=1;
    }
#endif

    if (!strcmp(argv[i],"-p")) {
       i++;
       found=1;
       VERBOSE=1;
    }       

    if (!strcmp(argv[i],"-d")) {
       if (i >= (argc-3)) {
          fprintf(stderr,"%s: not enough arguments for -d\n",MY_NAME);
          return -1;
       }
       i++;
       found=1;
#ifdef SGI
       DISPLAY = atoi(argv[i++]);
       if ((DISPLAY<0) || (DISPLAY>1)) {
          fprintf(stderr,"%s: error, illegal display value\n",MY_NAME);
          return -1;
       }
#else
       fprintf(stderr,"%s: warning, '-d' not a valid option for non-SGI\n",
               MY_NAME);
       i++;  /* skip over display value */
#endif
    }       

    if (!found) {
      printf("%s: error, illegal option %s\n",MY_NAME,argv[i]);
      return -1;
    }

  }

  /* can only select '-o' if image is displayed; check */
  if ((strcmp(OUTNAME,"")) && (!DISPLAY)) {
     fprintf(stderr,"%s: error, display must be on to save as raster image\n",
             MY_NAME);
     return -1;
  }

  return 0;

}




/**************************** usage ****************************/
/**************************** usage ****************************/
/**************************** usage ****************************/
/**************************** usage ****************************/

void usage()
{

  printf("usage: %s [-options] <3d_sds.hdf> <threshold>\n\n",
          MY_NAME);

#ifdef SGI
  printf("           -bc <br bg bb>               Background color\n");
  printf("           -d <0|1>                     Display off or on(default)\n");
  printf("           -k <Ka Kd Ks>                Ambient,Diffuse,Specular\n");
  printf("           -lc <lr lg lb>               Light color\n");
  printf("           -ll <lx ly lz>               Light location\n");
  printf("           -ns <Ns>                     Specular exponent\n");
  printf("           -ntsc                        NTSC window\n");
  printf("           -o <hdf.r24>                 Save raster image\n");
  printf("           -r <rx ry rz>                Rotation\n");
  printf("           -s <sx sy sz>                Scaling\n");
  printf("           -t <tx ty tz>                Translation\n");
  printf("           -v <xorg yorg xsize ysize>   Window location\n");
#endif
  printf("           -def                         Display default values\n");
  printf("           -p                           Print information\n");
  printf("           -vset <hdf.vset>             Save VSet polygons\n");
  printf("           -obj <obj.file>              Save obj polygons\n");
  printf("\n");

  printf("This program generates polygonal isosurfaces, for a given\n");
  printf("threshold, contained in the input file.  The isosurface\n");
#ifdef SGI
  printf("can be displayed on the SGI workstation and optionally.\n");
  printf("saved to a raster image file or to polygon files.\n\n");
#else
  printf("can be saved to polygon files.\n\n");
#endif

  printf("The input file must be an NCSA HDF file containing a\n");
  printf("three dimensional Scientific Data Set (SDS).\n");

}




/**************************** defaults ****************************/
/**************************** defaults ****************************/
/**************************** defaults ****************************/
/**************************** defaults ****************************/

void defaults()
{

  printf("\n\n");
  printf("%s: Default values\n",MY_NAME);
#ifdef SGI
  printf("           -r %f %f %f\n",ROTATION[0],ROTATION[1],ROTATION[2]);
  printf("           -t %f %f %f\n",TRANSLATION[0],
                                    TRANSLATION[1],TRANSLATION[2]);
  printf("           -s %f %f %f\n",SCALE[0],SCALE[1],SCALE[2]);
  printf("           -ll %f %f %f\n",LIGHT_LOC[0],LIGHT_LOC[1],
                                     LIGHT_LOC[2]);
  printf("           -lc %f %f %f\n",LIGHT_COLOR[0],LIGHT_COLOR[1],
                                     LIGHT_COLOR[2]);
  printf("           -k %f %f %f\n",K_VALUES[0],K_VALUES[1],K_VALUES[2]);
  printf("           -ns %f\n",K_EXP);
  printf("           -bc %f %f %f\n",BACKG_COLOR[0],BACKG_COLOR[1],
                                     BACKG_COLOR[2]);
  printf("           -d %d\n",DISPLAY);
#endif
  printf("\n");

}
