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

#ifdef MAIN

float ROTATION[3] = {0.0, 0.0, 0.0};
float TRANSLATION[3] = {0.0, 0.0, 0.0};
float SCALE[3] = {1.0, 1.0, 1.0};
float LIGHT_LOC[3] = {50.0, 50.0, 50.0};
float LIGHT_COLOR[3] = {1.0, 1.0, 1.0};
float MAT_COLOR[3] = {0.9, 0.9, 0.9};
float K_VALUES[3] = {0.3, 0.7, 0.0};
float K_EXP = 3.0;
int WIN_POS[2] = {-1, -1};
int WIN_SIZE[2] = {-1, -1};
int WIN_FOREGROUND = 1;
float BACKG_COLOR[3] = {0.0, 0.0, 0.0};
float FACE_COLOR[3] = {0.065,0.170,0.779};
int NORMAL_TYPE = 1;
char OUTNAME[80] = "";
char MY_NAME[80];
int VERBOSE = 0;
char VSET_NAME[80] = "";
char WFT_NAME[80] = "";
char R3D_NAME[80] = "";
int STORE_POLYGONS = 0;
int DTM_OUTPUT = 0;
char BYU_NAME[80] = "";
int XDIM, YDIM, ZDIM, RAW_INPUT;
float XINC=0.0,YINC=0.0,ZINC=0.0;
float XMIN=0.0,YMIN=0.0,ZMIN=0.0;
int BBOX = 0, AUXLIGHT = 0;
int SMOOTH = 0, NTSC = 0;
#if defined(SGI) && defined(SGI_IMAGE)
char SOUTNAME[80] = "";
#endif

int ORBITAL = 0;

#ifdef SGI
int DISPLAY = 1;
#else
int DISPLAY = 0;
#endif

#else

extern float ROTATION[3];
extern float TRANSLATION[3];
extern float SCALE[3];
extern float LIGHT_LOC[3];
extern float LIGHT_COLOR[3];
extern float MAT_COLOR[3];
extern float K_VALUES[3];
extern float K_EXP;
extern int WIN_POS[2];
extern int WIN_SIZE[2];
extern int WIN_FOREGROUND;
extern float BACKG_COLOR[3];
extern float FACE_COLOR[3];
extern int NORMAL_TYPE;
extern char OUTNAME[80];
extern char MY_NAME[80];
extern int VERBOSE;
extern char VSET_NAME[80];
extern char WFT_NAME[80];
extern char R3D_NAME[80];
extern int STORE_POLYGONS;
extern int ORBITAL;
extern int DISPLAY;
extern int DTM_OUTPUT;
extern char BYU_NAME[80];
extern int XDIM, YDIM, ZDIM, RAW_INPUT;
extern float XINC,YINC,ZINC;
extern float XMIN,YMIN,ZMIN;
extern int BBOX, AUXLIGHT;
extern int SMOOTH, NTSC;

#if defined(SGI) && defined(SGI_IMAGE)
extern char SOUTNAME[80];
#endif

#endif
