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
 * This program implements the marching cubes surface tiler described by
 * Lorensen & Cline in the Siggraph 87 Conference Proceedings.
 *
 * This program gets its data from a 3D Scientific Data Set in an
 * NCSA HDF file and creates two files, one containing vertices and
 * one containing connectivity.  The user must specify the three
 * file names and the threshold value on the command line.
 *
 * Written by Mike Krogh, NCSA, Feb.  2, 1990
 *
*/

#include <stdio.h>

typedef struct {
    int    nverts;
    int    verts[8];
    int    nedges;
    int    edges[12];
    int    npolys;
    int    polys[30];
    } CELL_ENTRY;
#include "cell_table.h"


#include "isovis.h"




/**************************** Temporary Globals ****************************/
/**************************** Temporary Globals ****************************/
/**************************** Temporary Globals ****************************/
/**************************** Temporary Globals ****************************/

float DATA1,DATA2,DATA3,DATA4,DATA5,DATA6,DATA7,DATA8;
int XDIMYDIM;





/**************************** iso_surface ****************************/
/**************************** iso_surface ****************************/
/**************************** iso_surface ****************************/
/**************************** iso_surface ****************************/

int iso_surface(data,xdim,ydim,zdim,threshold)
float *data;
int xdim,ydim,zdim;
float threshold;
/* This subroutine will generate the polygonal isosurface. */
{

  register int x,y,z;
  register int xdim1,ydim1,zdim1;
  float xtrans,ytrans,ztrans;        /* used for centering the volume */
  int index;
  int npolys;
  int old_verts;
  int edge[13];
  float crossings[13][3];
  void get_cell_verts();
  void get_cell_polys();
  void calc_index_and_temps();

  zdim1 = zdim-1;
  ydim1 = ydim-1;
  xdim1 = xdim-1;

  xtrans = -((float)xdim/2.0);
  ytrans = -((float)ydim/2.0);
  ztrans = -((float)zdim/2.0);

  XDIMYDIM = xdim*ydim;

  npolys=0;     /* keep count of total polygons */

  for (z=0;z<zdim1;z++)      /* process each cell in the volume */
    for (y=0;y<ydim1;y++)
      for (x=0;x<xdim1;x++) {
        calc_index_and_temps(data,x,y,z,xdim,ydim,zdim,threshold,&index);
        if (index) {
          get_cell_verts(index,x,y,z,xtrans,ytrans,ztrans,
                         threshold,crossings);
          get_cell_polys(index,&npolys,crossings);
        }
      }

  /* record some statistics */
  if (VERBOSE)
     printf("%s: %d polygons generated\n",MY_NAME,npolys);

  return 0;

}


/**************************** calc_index_and_temps ****************************/
/**************************** calc_index_and_temps ****************************/
/**************************** calc_index_and_temps ****************************/
/**************************** calc_index_and_temps ****************************/

void calc_index_and_temps(data,x1,y1,z1,xdim,ydim,zdim,threshold,index)
register float *data;
int x1,y1,z1;
int xdim,ydim,zdim;
register float threshold;
int *index;
/* This subroutine calculates the index and creates some global */
/* temporary variables (for speed). */
{

  register float *tmp;

  *index = 0;

  tmp = data + (z1*XDIMYDIM) + (y1*xdim) + x1;

  *index += (threshold <= (DATA1 = *(tmp)));
  *index += (threshold <= (DATA2 = *(tmp + 1))) * 2;

  tmp += xdim;
  *index += (threshold <= (DATA3 = *(tmp + 1))) * 4;
  *index += (threshold <= (DATA4 = *(tmp))) * 8;

  tmp = tmp - xdim + XDIMYDIM;
  *index += (threshold <= (DATA5 = *(tmp))) * 16;
  *index += (threshold <= (DATA6 = *(tmp + 1))) * 32;

  tmp += xdim;
  *index += (threshold <= (DATA7 = *(tmp + 1))) * 64;
  *index += (threshold <= (DATA8 = *(tmp))) * 128;
 
}




/**************************** get_cell_verts ****************************/
/**************************** get_cell_verts ****************************/
/**************************** get_cell_verts ****************************/
/**************************** get_cell_verts ****************************/

void get_cell_verts(index,x1,y1,z1,xtrans,ytrans,ztrans,threshold,crossings)
int index;
int x1,y1,z1;
float xtrans,ytrans,ztrans;
float threshold;
float crossings[13][3];
{

  register int i;
  register int x2,y2,z2;
  int nedges;
  int crnt_edge;
 
#define linterp(a1,a2,a,b1,b2) ((float)(((a-a1) * (float)(b2-b1) / (a2-a1)) + (float)b1))

  x2 = x1+1;
  y2 = y1+1;
  z2 = z1+1;

  nedges = cell_table[index].nedges;
  for (i=0;i<nedges;i++) {
     crnt_edge = cell_table[index].edges[i];
     switch (crnt_edge) {
	case 1:
		crossings[1][0] = linterp(DATA1,DATA2,threshold,x1,x2)+xtrans;
		crossings[1][1] = (float)y1+ytrans;
		crossings[1][2] = (float)z1+ztrans;
		break;

	case 2:
		crossings[2][1] = linterp(DATA2,DATA3,threshold,y1,y2)+ytrans;
		crossings[2][0] = (float)x2+xtrans;
		crossings[2][2] = (float)z1+ztrans;
		break;

	case 3:
		crossings[3][0] = linterp(DATA4,DATA3,threshold,x1,x2)+xtrans;
		crossings[3][1] = (float)y2+ytrans;
		crossings[3][2] = (float)z1+ztrans;
		break;

	case 4:
		crossings[4][1] = linterp(DATA1,DATA4,threshold,y1,y2)+ytrans;
		crossings[4][0] = (float)x1+xtrans;
		crossings[4][2] = (float)z1+ztrans;
		break;

	case 5:
		crossings[5][0] = linterp(DATA5,DATA6,threshold,x1,x2)+xtrans;
		crossings[5][1] = (float)y1+ytrans;
		crossings[5][2] = (float)z2+ztrans;
		break;

	case 6:
		crossings[6][1] = linterp(DATA6,DATA7,threshold,y1,y2)+ytrans;
		crossings[6][0] = (float)x2+xtrans;
		crossings[6][2] = (float)z2+ztrans;
		break;

	case 7:
		crossings[7][0] = linterp(DATA8,DATA7,threshold,x1,x2)+xtrans;
		crossings[7][1] = (float)y2+ytrans;
		crossings[7][2] = (float)z2+ztrans;
		break;

	case 8:
		crossings[8][1] = linterp(DATA5,DATA8,threshold,y1,y2)+ytrans;
		crossings[8][0] = (float)x1+xtrans;
		crossings[8][2] = (float)z2+ztrans;
		break;

	case 9:
		crossings[9][2] = linterp(DATA1,DATA5,threshold,z1,z2)+ztrans;
		crossings[9][1] = (float)y1+ytrans;
		crossings[9][0] = (float)x1+xtrans;
		break;

	case 10:
		crossings[10][2] = linterp(DATA2,DATA6,threshold,z1,z2)+ztrans;
		crossings[10][1] = (float)y1+ytrans;
		crossings[10][0] = (float)x2+xtrans;
		break;

	case 11:
		crossings[11][2] = linterp(DATA4,DATA8,threshold,z1,z2)+ztrans;
		crossings[11][1] = (float)y2+ytrans;
		crossings[11][0] = (float)x1+xtrans;
		break;

	case 12:
		crossings[12][2] = linterp(DATA3,DATA7,threshold,z1,z2)+ztrans;
		crossings[12][1] = (float)y2+ytrans;
		crossings[12][0] = (float)x2+xtrans;
		break;

    } /* end switch */
  } /* end for */
}




/**************************** get_cell_polys ****************************/
/**************************** get_cell_polys ****************************/
/**************************** get_cell_polys ****************************/
/**************************** get_cell_polys ****************************/

void get_cell_polys(index,npolys,crossings)
int index;
int *npolys;
float crossings[13][3];
/* This subroutine will calculate the polygons */
{

  register int num_o_polys;
  register int poly;
  float *p1,*p2,*p3;
  float n1[3],n2[3],n3[3];

  int add_polygon();  /* stores polygons for file output */

#ifdef SGI
  void draw_triangle();
  void calc_normal();
#endif

  num_o_polys = cell_table[index].npolys;
  for (poly=0;poly<num_o_polys;poly++) {

    p1 = &crossings[cell_table[index].polys[(poly*3)]][0];
    p2 = &crossings[cell_table[index].polys[(poly*3)+1]][0];
    p3 = &crossings[cell_table[index].polys[(poly*3)+2]][0];

    if (!NORMAL_TYPE) {
       /* store the polygons if required */
       if (STORE_POLYGONS)
          if (add_polygon(p1,p2,p3) == -1) {
             fprintf(stderr,"%s: error from add_polygon\n",MY_NAME);
             exit(1);
          }

#ifdef SGI
       if (DISPLAY) {
          calc_normal(p1,p2,p3,n1);
          draw_triangle(p1,p2,p3,n1,n1,n1);
       }
#endif

    } else {
       fprintf(stderr,"%s: error, unsupported - gradient normals\n",
               MY_NAME);
       exit(1);
    }

  }

  (*npolys) += num_o_polys;

}

