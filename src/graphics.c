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
#ifdef SGI
#include <gl.h>
#include <device.h>
#endif
#include <math.h>
#include "isovis.h"



/**************************** init_graphics ****************************/
/**************************** init_graphics ****************************/
/**************************** init_graphics ****************************/
/**************************** init_graphics ****************************/

#ifdef SGI

int init_graphics(xdim,ydim,zdim)
int xdim,ydim,zdim;
{
    int gid;

    void define_lights();
    void use_lights();

    static Matrix idmat = { 1.0, 0.0, 0.0, 0.0,
                            0.0, 1.0, 0.0, 0.0,
                            0.0, 0.0, 1.0, 0.0,
                            0.0, 0.0, 0.0, 1.0};  /* identity matrix */

    if (WIN_POS[0] != -1)
       prefposition (WIN_POS[0], WIN_SIZE[0]+WIN_POS[0],
                     WIN_POS[1], WIN_SIZE[1]+WIN_POS[1]);

    foreground();
    gid = winopen (MY_NAME);
    if (gid == -1) {
       fprintf(stderr,"%s: error, can't open a graphics window\n",
               MY_NAME);
       return -1;
    }

    RGBmode();
    gconfig();
    zbuffer(1);
    lsetdepth(0x0000,0xfffff);

    mmode(MVIEWING);
    perspective(450,1.0,0.1,500);
    loadmatrix(idmat);
    lookat((float)(xdim/2),(float)(ydim/2),(float)(zdim+100),
           0.0,0.0,0.0,0.0);

    define_lights();

    if (NORMAL_TYPE == 0)
       shademodel(FLAT);

    zclear();
    RGBcolor( (int)(BACKG_COLOR[0] * 255.0),
              (int)(BACKG_COLOR[1] * 255.0),
              (int)(BACKG_COLOR[2] * 255.0) );
    clear();


    scale(SCALE[0],SCALE[1],SCALE[2]);

    /* can't use 'rot' due to bug in DGL for the Cray */
    rotate( ((int)ROTATION[0]*10),'x');
    rotate( ((int)ROTATION[1]*10),'y');
    rotate( ((int)ROTATION[2]*10),'z');

    translate(TRANSLATION[0],TRANSLATION[1],TRANSLATION[2]);

    if (VERBOSE) {
       printf("%s: scaling %f %f %f\n",MY_NAME,SCALE[0],
              SCALE[1],SCALE[2]);
       printf("%s: rotation %f %f %f\n",MY_NAME,ROTATION[0],
              ROTATION[1],ROTATION[2]);
       printf("%s: translation %f %f %f\n",MY_NAME,TRANSLATION[0],
              TRANSLATION[1],TRANSLATION[2]);
    }
    /* Queue devices for <esc> key as an exit, and redraw events */
    qdevice(ESCKEY);
    qdevice(REDRAW);

    return 0;

}




/**************************** define_lights ****************************/
/**************************** define_lights ****************************/
/**************************** define_lights ****************************/
/**************************** define_lights ****************************/

void define_lights()
{

  float light[10];
  float material[15];

  light[0] = POSITION;
  light[1] = LIGHT_LOC[0];
  light[2] = LIGHT_LOC[1];
  light[3] = LIGHT_LOC[2];
  light[4] = 0.0;
  light[5] = LCOLOR;
  light[6] = LIGHT_COLOR[0];
  light[7] = LIGHT_COLOR[1];
  light[8] = LIGHT_COLOR[2];
  light[9] = LMNULL;

  material[0] = SPECULAR;
  material[1] = K_VALUES[2];
  material[2] = K_VALUES[2];
  material[3] = K_VALUES[2];
  material[4] = AMBIENT;
  material[5] = K_VALUES[0];
  material[6] = K_VALUES[0];
  material[7] = K_VALUES[0];
  material[8] = DIFFUSE;
  material[9] = K_VALUES[1];
  material[10] = K_VALUES[1];
  material[11] = K_VALUES[1];
  material[12] = SHININESS;
  material[13] = K_EXP;
  material[14] = LMNULL;

  lmdef(DEFMATERIAL,1,15,material);
  lmdef(DEFLIGHT,1,10,light);
  lmdef(DEFLMODEL,1,0,NULL);

  lmbind(MATERIAL,1);
  lmbind(LIGHT0,1);
  lmbind(LMODEL,1);

}




/**************************** draw_triangle ****************************/
/**************************** draw_triangle ****************************/
/**************************** draw_triangle ****************************/
/**************************** draw_triangle ****************************/

void draw_triangle(p1,p2,p3,n1,n2,n3)
float *p1,*p2,*p3;
float *n1,*n2,*n3;
{

  bgnpolygon();
    n3f(n1);
    v3f(p1);
    n3f(n2);
    v3f(p2);
    n3f(n3);
    v3f(p3);
  endpolygon();

}




/**************************** get_image ****************************/
/**************************** get_image ****************************/
/**************************** get_image ****************************/
/**************************** get_image ****************************/

int get_image(image,xdim,ydim)
char **image;
int *xdim, *ydim;
{

  long xsize,ysize;
  register int x,y,t,i,j,y2;
  register long *sgi_image;
  long ret;

  getsize(&xsize,&ysize);
  *xdim = (int)xsize;     *ydim = (int)ysize;

  if ((sgi_image = (long *)malloc((*xdim)*(*ydim)*sizeof(long))) == NULL) {
     fprintf(stderr,"%s: error, not enough memory for sgi_image\n",MY_NAME);
     return -1;
  }

  ret = lrectread(0,0,(xsize-1),(ysize-1),sgi_image);
  if (ret != (xsize*ysize)) {
     fprintf(stderr,"%s: error from lrectread %ld\n",MY_NAME,ret);
     return -1;
  }

  if ((*image = (char *)malloc((*xdim)*(*ydim)*3)) == NULL) {
     fprintf(stderr,"%s: error, not enough memory for image\n",MY_NAME);
     return -1;
  }

  y2 = (*ydim)-1;
  i=0;
  for (y=y2;y>=0;y--) {
    t=y*(*xdim);
    for (x=0;x<(*xdim);x++) {
      j=t+x;
      *(*image+i) = (char)( *(sgi_image+j) & (long)0xff);
      i++;
      *(*image+i) = (char)(( *(sgi_image+j) & (long)0xff00)>>8);
      i++;
      *(*image+i) = (char)(( *(sgi_image+j) & (long)0xff0000)>>16);
      i++;
    }
  }

  return 0;

}


void window_event_loop(data,xdim,ydim,zdim,threshold)
float *data;         /* the volume of data */
int xdim,ydim,zdim;  /* the dimensions of the data */
float threshold;     /* the current threshold of interest */
/* This subroutine is called by main and implements a simple event loop
 * that handles screen refreshes and <esc> key for exiting.
 */
{
  short value;
  int dev;
  int iso_surface();

  /* do a simple event loop */
  while (TRUE) {
    while (qtest()) {        /* look for an event */
      dev=qread(&value);     /* get the event */
      if (dev == ESCKEY)     /* <esc> key was pressed */
         exit(0);
      if (dev == REDRAW) {   /* window needs redrawing */
         reshapeviewport();
         STORE_POLYGONS=0;   /* do not resave any polygons. */
         zclear();           /* clear the z-buffer */
         RGBcolor( (int)(BACKG_COLOR[0] * 255.0),
                   (int)(BACKG_COLOR[1] * 255.0),
                   (int)(BACKG_COLOR[2] * 255.0) );
         clear();            /* reset the background color */
         if (iso_surface(data,xdim,ydim,zdim,threshold) == -1) {
            printf("%s: error from iso_surface\n",MY_NAME);
            exit(1);
         }
      } /* redraw */
    }  /* qtest */
  } /* while true */

}

#endif /* SGI */
