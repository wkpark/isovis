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
 * This program gets its data from a 3D Scientific Data Set in an NCSA HDF file
 * and creates two files, one containing vertices and one containing
 * connectivity.  The user must specify the three file names and the
 * threshold value on the command line. 
 *
 * Written by Mike Krogh, NCSA, Feb.  2, 1990 
 *
 */

#include <stdio.h>
#include <math.h>

typedef struct {
    int nverts;
    int verts[8];
    int nedges;
    int edges[12];
    int npolys;
    int polys[30];
} CELL_ENTRY;
#include "cell_table.h"

#include "isovis.h"

float DATA1, DATA2, DATA3, DATA4, DATA5, DATA6, DATA7, DATA8;
int XDIMYDIM;

int 
iso_surface(data, xdim, ydim, zdim, threshold)
float *data;
int xdim, ydim, zdim;
float threshold;
/* This subroutine will generate the polygonal isosurface. */
{

    register int x, y, z;
    register int xdim1, ydim1, zdim1;
    float xtrans, ytrans, ztrans;	/* used for centering the volume */
    int index;
    int npolys;
    float crossings[13][3];
    float cell_norms[8][3];
    float normals[13][3];
    void get_cell_verts();
    void get_cell_polys();
    void get_cell_polys_and_norms();
    int calc_index_and_temps();
    void calc_cell_normals();
    void get_cell_norms();

    zdim1 = zdim - 1;
    ydim1 = ydim - 1;
    xdim1 = xdim - 1;

    xtrans = -((float) xdim / 2.0);
    ytrans = -((float) ydim / 2.0);
    ztrans = -((float) zdim / 2.0);

    XDIMYDIM = xdim * ydim;

    npolys = 0;			/* keep count of total polygons */

    if (NORMAL_TYPE) {
	/* gradient normals - scary stuff kids */
	for (z = 0; z < zdim1; z++)	/* process each cell in the volume */
	    for (y = 0; y < ydim1; y++)
		for (x = 0; x < xdim1; x++) {
		    index = calc_index_and_temps(data, x, y, z, xdim, ydim, zdim, threshold);
		    if (index) {
			calc_cell_normals(data, x, y, z, xdim, ydim, zdim, cell_norms);
			get_cell_verts(index, x, y, z, xtrans, ytrans, ztrans,
				       threshold, crossings);
			get_cell_norms(index, x, y, z, threshold, cell_norms, normals);
			get_cell_polys_and_norms(index, &npolys, crossings, normals);
		    }
		}
    } else {
	for (z = 0; z < zdim1; z++)	/* process each cell in the volume */
	    for (y = 0; y < ydim1; y++)
		for (x = 0; x < xdim1; x++) {
		    index = calc_index_and_temps(data, x, y, z, xdim, ydim, zdim, threshold);
		    if (index) {
			get_cell_verts(index, x, y, z, xtrans, ytrans, ztrans,
				       threshold, crossings);
			get_cell_polys(index, &npolys, crossings);
		    }
		}
    }

    /* record some statistics */
    if (VERBOSE)
	printf("%s: %d polygons generated\n", MY_NAME, npolys);

    /* don't do this when its time for multiple iso-surfaces */
//    free((char *)data);
    return 0;

}


calc_index_and_temps(data, x1, y1, z1, xdim, ydim, zdim, thresh)
register float *data;
int x1, y1, z1;
int xdim, ydim, zdim;
float thresh;
/* This subroutine calculates the index and creates some global */
/* temporary variables (for speed). */
{

    register float *tmp;
    int index = 0;
    register float threshold = thresh;	/* thresh is really a double :-( */

    tmp = data + (z1 * XDIMYDIM) + (y1 * xdim) + x1;

    index += (threshold <= (DATA1 = *(tmp)));
    index += (threshold <= (DATA2 = *(tmp + 1))) * 2;

    tmp += xdim;
    index += (threshold <= (DATA3 = *(tmp + 1))) * 4;
    index += (threshold <= (DATA4 = *(tmp))) * 8;

    tmp = tmp - xdim + XDIMYDIM;
    index += (threshold <= (DATA5 = *(tmp))) * 16;
    index += (threshold <= (DATA6 = *(tmp + 1))) * 32;

    tmp += xdim;
    index += (threshold <= (DATA7 = *(tmp + 1))) * 64;
    index += (threshold <= (DATA8 = *(tmp))) * 128;

    return index;

}

void 
get_cell_verts(index, x1, y1, z1, xtrans, ytrans, ztrans, threshold, crossings)
int index;
int x1, y1, z1;
float xtrans, ytrans, ztrans;
float threshold;
float crossings[13][3];
{

    register int i;
    register int x2, y2, z2;
    int nedges;
    int crnt_edge;

#define linterp(a1,a2,a,b1,b2) ((float)(((a-a1) * (float)(b2-b1) / (a2-a1)) + (float)b1))

    x2 = x1 + 1;
    y2 = y1 + 1;
    z2 = z1 + 1;

    nedges = cell_table[index].nedges;
    for (i = 0; i < nedges; i++) {
	crnt_edge = cell_table[index].edges[i];
	switch (crnt_edge) {
	case 1:
	    crossings[1][0] = linterp(DATA1, DATA2, threshold, x1, x2) + xtrans;
	    crossings[1][1] = (float) y1 + ytrans;
	    crossings[1][2] = (float) z1 + ztrans;
	    break;

	case 2:
	    crossings[2][1] = linterp(DATA2, DATA3, threshold, y1, y2) + ytrans;
	    crossings[2][0] = (float) x2 + xtrans;
	    crossings[2][2] = (float) z1 + ztrans;
	    break;

	case 3:
	    crossings[3][0] = linterp(DATA4, DATA3, threshold, x1, x2) + xtrans;
	    crossings[3][1] = (float) y2 + ytrans;
	    crossings[3][2] = (float) z1 + ztrans;
	    break;

	case 4:
	    crossings[4][1] = linterp(DATA1, DATA4, threshold, y1, y2) + ytrans;
	    crossings[4][0] = (float) x1 + xtrans;
	    crossings[4][2] = (float) z1 + ztrans;
	    break;

	case 5:
	    crossings[5][0] = linterp(DATA5, DATA6, threshold, x1, x2) + xtrans;
	    crossings[5][1] = (float) y1 + ytrans;
	    crossings[5][2] = (float) z2 + ztrans;
	    break;

	case 6:
	    crossings[6][1] = linterp(DATA6, DATA7, threshold, y1, y2) + ytrans;
	    crossings[6][0] = (float) x2 + xtrans;
	    crossings[6][2] = (float) z2 + ztrans;
	    break;

	case 7:
	    crossings[7][0] = linterp(DATA8, DATA7, threshold, x1, x2) + xtrans;
	    crossings[7][1] = (float) y2 + ytrans;
	    crossings[7][2] = (float) z2 + ztrans;
	    break;

	case 8:
	    crossings[8][1] = linterp(DATA5, DATA8, threshold, y1, y2) + ytrans;
	    crossings[8][0] = (float) x1 + xtrans;
	    crossings[8][2] = (float) z2 + ztrans;
	    break;

	case 9:
	    crossings[9][2] = linterp(DATA1, DATA5, threshold, z1, z2) + ztrans;
	    crossings[9][1] = (float) y1 + ytrans;
	    crossings[9][0] = (float) x1 + xtrans;
	    break;

	case 10:
	    crossings[10][2] = linterp(DATA2, DATA6, threshold, z1, z2) + ztrans;
	    crossings[10][1] = (float) y1 + ytrans;
	    crossings[10][0] = (float) x2 + xtrans;
	    break;

	case 11:
	    crossings[11][2] = linterp(DATA4, DATA8, threshold, z1, z2) + ztrans;
	    crossings[11][1] = (float) y2 + ytrans;
	    crossings[11][0] = (float) x1 + xtrans;
	    break;

	case 12:
	    crossings[12][2] = linterp(DATA3, DATA7, threshold, z1, z2) + ztrans;
	    crossings[12][1] = (float) y2 + ytrans;
	    crossings[12][0] = (float) x2 + xtrans;
	    break;

	} /* end switch */
    } /* end for */
}

void 
get_cell_polys(index, npolys, crossings)
int index;
int *npolys;
float crossings[13][3];
/* This subroutine will calculate the polygons */
{

    register int num_o_polys;
    register int poly;
    float *p1, *p2, *p3;
    float n1[3], n2[3], n3[3];

    int add_polygon();		/* stores polygons for file output */
    void calc_normal();		/* calculate polygon normal */

    num_o_polys = cell_table[index].npolys;
    for (poly = 0; poly < num_o_polys; poly++) {

	p1 = &crossings[cell_table[index].polys[(poly * 3)]][0];
	p2 = &crossings[cell_table[index].polys[(poly * 3) + 1]][0];
	p3 = &crossings[cell_table[index].polys[(poly * 3) + 2]][0];
#define CULL
#ifdef CULL
	if ((p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2]) ||
	    (p1[0] == p3[0] && p1[1] == p3[1] && p1[2] == p3[2]) ||
	    (p2[0] == p3[0] && p2[1] == p3[1] && p2[2] == p3[2]))  {
	    (*npolys)--;
	    continue;
	}
#endif

	calc_normal(p1, p2, p3, n1);
	if (add_polygon(p1, p2, p3, n1, n1, n1) == -1) {
	    fprintf(stderr, "%s: error from add_polygon\n", MY_NAME);
	    exit(1);
	}
    }

    (*npolys) += num_o_polys;
}

#ifdef ndef
void 
calc_cell_normals(data, x1, y1, z1, xdim, ydim, zdim, norms)
register float *data;
int x1, y1, z1;
int xdim, ydim, zdim;
float norms[8][3];
/* This subroutine calculates the vertex normals using a central difference */
{

    int i;
#define D(x,y,z) *(data+(z)*XDIMYDIM + (y)*xdim + (x))

    /* clamp */
    if (x1 == 0) x1++;
    else if (x1 == xdim - 2) x1 = xdim - 3;
    if (y1 == 0) y1++;
    else if (y1 == ydim - 2) y1 = ydim - 3;
    if (z1 == 0) z1++;
    else if (z1 == zdim - 2) z1 = zdim - 3;

    norms[0][0] = D(x1 + 1, y1, z1) - D(x1 - 1, y1, z1);
    norms[0][1] = D(x1, y1 + 1, z1) - D(x1, y1 - 1, z1);
    norms[0][2] = D(x1, y1, z1 + 1) - D(x1, y1, z1 - 1);

    norms[1][0] = D(x1 + 1 + 1, y1, z1) - D(x1 + 1 - 1, y1, z1);
    norms[1][1] = D(x1 + 1, y1 + 1, z1) - D(x1 + 1, y1 - 1, z1);
    norms[1][2] = D(x1 + 1, y1, z1 + 1) - D(x1 + 1, y1, z1 - 1);

    norms[2][0] = D(x1 + 1 + 1, y1 + 1, z1) - D(x1 + 1 - 1, y1 + 1, z1);
    norms[2][1] = D(x1 + 1, y1 + 1 + 1, z1) - D(x1 + 1, y1 + 1 - 1, z1);
    norms[2][2] = D(x1 + 1, y1 + 1, z1 + 1) - D(x1 + 1, y1 + 1, z1 - 1);

    norms[3][0] = D(x1 + 1, y1 + 1, z1) - D(x1 - 1, y1 + 1, z1);
    norms[3][1] = D(x1, y1 + 1 + 1, z1) - D(x1, y1 + 1 - 1, z1);
    norms[3][2] = D(x1, y1 + 1, z1 + 1) - D(x1, y1 + 1, z1 - 1);

    norms[4][0] = D(x1 + 1, y1, z1 + 1) - D(x1 - 1, y1, z1 + 1);
    norms[4][1] = D(x1, y1 + 1, z1 + 1) - D(x1, y1 - 1, z1 + 1);
    norms[4][2] = D(x1, y1, z1 + 1 + 1) - D(x1, y1, z1 + 1 - 1);

    norms[5][0] = D(x1 + 1 + 1, y1, z1 + 1) - D(x1 + 1 - 1, y1, z1 + 1);
    norms[5][1] = D(x1 + 1, y1 + 1, z1 + 1) - D(x1 + 1, y1 - 1, z1 + 1);
    norms[5][2] = D(x1 + 1, y1, z1 + 1 + 1) - D(x1 + 1, y1, z1 + 1 - 1);

    norms[6][0] = D(x1 + 1 + 1, y1 + 1, z1 + 1) - D(x1 + 1 - 1, y1 + 1, z1 + 1);
    norms[6][1] = D(x1 + 1, y1 + 1 + 1, z1 + 1) - D(x1 + 1, y1 + 1 - 1, z1 + 1);
    norms[6][2] = D(x1 + 1, y1 + 1, z1 + 1 + 1) - D(x1 + 1, y1 + 1, z1 + 1 - 1);

    norms[7][0] = D(x1 + 1, y1 + 1, z1 + 1) - D(x1 - 1, y1 + 1, z1 + 1);
    norms[7][1] = D(x1, y1 + 1 + 1, z1 + 1) - D(x1, y1 + 1 - 1, z1 + 1);
    norms[7][2] = D(x1, y1 + 1, z1 + 1 + 1) - D(x1, y1 + 1, z1 + 1 - 1);
    if (NORMAL_TYPE > 1)
	for (i = 0; i < 8; i++) {
	    norms[i][0] = -norms[i][0];
	    norms[i][1] = -norms[i][1];
	    norms[i][2] = -norms[i][2];
	}
    for (i = 0; i < 8; i++) {
#ifdef mips
	float d = fsqrt(norms[i][0] * norms[i][0] +
		        norms[i][1] * norms[i][1] +
		        norms[i][2] * norms[i][2]);
#else
	float d = sqrt(norms[i][0] * norms[i][0] +
		       norms[i][1] * norms[i][1] +
		       norms[i][2] * norms[i][2]);
#endif
	if (d < 1.e-6) {
	    printf("small normal alert\n");
	} else {
	    norms[i][0] /= d;
	    norms[i][1] /= d;
	    norms[i][2] /= d;
	}
    }
}

#else
void 
calc_cell_normals(data, x1, y1, z1, xdim, ydim, zdim, norms)
register float *data;
int x1, y1, z1;
int xdim, ydim, zdim;
float norms[8][3];
/* This subroutine calculates the vertex normals using a central difference */
{

    int i;
    int x11, y11, z11;
#define D(x,y,z) *(data+(z)*XDIMYDIM + (y)*xdim + (x))

    /* clamp */
    x11 = x1+1; y11 = y1+1; z11 = z1+1;
    if (x1 == 0) x1++;
    else if (x1 == xdim - 2) x11--;
    if (y1 == 0) y1++;
    else if (y1 == ydim - 2) y11--;
    if (z1 == 0) z1++;
    else if (z1 == zdim - 2) z11--;

    norms[0][0] = D(x1 + 1, y1, z1) - D(x1 - 1, y1, z1);
    norms[3][0] = D(x1 + 1, y11, z1) - D(x1 - 1, y11, z1);
    norms[4][0] = D(x1 + 1, y1, z11) - D(x1 - 1, y1, z11);
    norms[7][0] = D(x1 + 1, y11, z11) - D(x1 - 1, y11, z11);

    norms[1][0] = D(x11 + 1, y1, z1) - D(x11 - 1, y1, z1);
    norms[2][0] = D(x11 + 1, y11, z1) - D(x11 - 1, y11, z1);
    norms[5][0] = D(x11 + 1, y1, z11) - D(x11 - 1, y1, z11);
    norms[6][0] = D(x11 + 1, y11, z11) - D(x11 - 1, y11, z11);

    norms[0][1] = D(x1, y1 + 1, z1) - D(x1, y1 - 1, z1);
    norms[1][1] = D(x11, y1 + 1, z1) - D(x11, y1 - 1, z1);
    norms[4][1] = D(x1, y1 + 1, z11) - D(x1, y1 - 1, z11);
    norms[5][1] = D(x11, y1 + 1, z11) - D(x11, y1 - 1, z11);

    norms[2][1] = D(x11, y11 + 1, z1) - D(x11, y11 - 1, z1);
    norms[3][1] = D(x1, y11 + 1, z1) - D(x1, y11 - 1, z1);
    norms[7][1] = D(x1, y11 + 1, z11) - D(x1, y11 - 1, z11);
    norms[6][1] = D(x11, y11 + 1, z11) - D(x11, y11 - 1, z11);

    norms[0][2] = D(x1, y1, z1 + 1) - D(x1, y1, z1 - 1);
    norms[1][2] = D(x11, y1, z1 + 1) - D(x11, y1, z1 - 1);
    norms[2][2] = D(x11, y11, z1 + 1) - D(x11, y11, z1 - 1);
    norms[3][2] = D(x1, y11, z1 + 1) - D(x1, y11, z1 - 1);

    norms[4][2] = D(x1, y1, z11 + 1) - D(x1, y1, z11 - 1);
    norms[5][2] = D(x11, y1, z11 + 1) - D(x11, y1, z11 - 1);
    norms[6][2] = D(x11, y11, z11 + 1) - D(x11, y11, z11 - 1);
    norms[7][2] = D(x1, y11, z11 + 1) - D(x1, y11, z11 - 1);

    if (NORMAL_TYPE > 1)
	for (i = 0; i < 8; i++) {
	    norms[i][0] = -norms[i][0];
	    norms[i][1] = -norms[i][1];
	    norms[i][2] = -norms[i][2];
	}
    for (i = 0; i < 8; i++) {
#ifdef mips
	float d = fsqrt(norms[i][0] * norms[i][0] +
		        norms[i][1] * norms[i][1] +
		        norms[i][2] * norms[i][2]);
#else
	float d = sqrt(norms[i][0] * norms[i][0] +
		       norms[i][1] * norms[i][1] +
		       norms[i][2] * norms[i][2]);
#endif
	if (d != 0.) {
	    norms[i][0] /= d;
	    norms[i][1] /= d;
	    norms[i][2] /= d;
	}
    }
}
#endif

void 
get_cell_norms(index, x1, y1, z1, thresh, cnorm, normals)
int index;
int x1, y1, z1;
float thresh;
float cnorm[8][3];
float normals[13][3];
{

    register int i;
    int nedges;
    int crnt_edge;
    register float threshold = thresh;	/* thresh is double */

#define lerp(a1,a2,a,b1,b2) ((a-a1)/(a2-a1)*(b2-b1) + b1)

    nedges = cell_table[index].nedges;
    for (i = 0; i < nedges; i++) {
	crnt_edge = cell_table[index].edges[i];
	switch (crnt_edge) {
	case 1:
	normals[1][0] = lerp(DATA1, DATA2, threshold, cnorm[0][0], cnorm[1][0]);
	normals[1][1] = lerp(DATA1, DATA2, threshold, cnorm[0][1], cnorm[1][1]);
	normals[1][2] = lerp(DATA1, DATA2, threshold, cnorm[0][2], cnorm[1][2]);
	break;

	case 2:
	normals[2][0] = lerp(DATA2, DATA3, threshold, cnorm[1][0], cnorm[2][0]);
	normals[2][1] = lerp(DATA2, DATA3, threshold, cnorm[1][1], cnorm[2][1]);
	normals[2][2] = lerp(DATA2, DATA3, threshold, cnorm[1][2], cnorm[2][2]);
	break;

	case 3:
	normals[3][0] = lerp(DATA3, DATA4, threshold, cnorm[2][0], cnorm[3][0]);
	normals[3][1] = lerp(DATA3, DATA4, threshold, cnorm[2][1], cnorm[3][1]);
	normals[3][2] = lerp(DATA3, DATA4, threshold, cnorm[2][2], cnorm[3][2]);
	break;

	case 4:
	normals[4][0] = lerp(DATA4, DATA1, threshold, cnorm[3][0], cnorm[0][0]);
	normals[4][1] = lerp(DATA4, DATA1, threshold, cnorm[3][1], cnorm[0][1]);
	normals[4][2] = lerp(DATA4, DATA1, threshold, cnorm[3][2], cnorm[0][2]);
	break;

	case 5:
	normals[5][0] = lerp(DATA5, DATA6, threshold, cnorm[4][0], cnorm[5][0]);
	normals[5][1] = lerp(DATA5, DATA6, threshold, cnorm[4][1], cnorm[5][1]);
	normals[5][2] = lerp(DATA5, DATA6, threshold, cnorm[4][2], cnorm[5][2]);
	break;

	case 6:
	normals[6][0] = lerp(DATA6, DATA7, threshold, cnorm[5][0], cnorm[6][0]);
	normals[6][1] = lerp(DATA6, DATA7, threshold, cnorm[5][1], cnorm[6][1]);
	normals[6][2] = lerp(DATA6, DATA7, threshold, cnorm[5][2], cnorm[6][2]);
	break;

	case 7:
	normals[7][0] = lerp(DATA7, DATA8, threshold, cnorm[6][0], cnorm[7][0]);
	normals[7][1] = lerp(DATA7, DATA8, threshold, cnorm[6][1], cnorm[7][1]);
	normals[7][2] = lerp(DATA7, DATA8, threshold, cnorm[6][2], cnorm[7][2]);
	break;

	case 8:
	normals[8][0] = lerp(DATA8, DATA5, threshold, cnorm[7][0], cnorm[4][0]);
	normals[8][1] = lerp(DATA8, DATA5, threshold, cnorm[7][1], cnorm[4][1]);
	normals[8][2] = lerp(DATA8, DATA5, threshold, cnorm[7][2], cnorm[4][2]);
	break;

	case 9:
	normals[9][0] = lerp(DATA1, DATA5, threshold, cnorm[0][0], cnorm[4][0]);
	normals[9][1] = lerp(DATA1, DATA5, threshold, cnorm[0][1], cnorm[4][1]);
	normals[9][2] = lerp(DATA1, DATA5, threshold, cnorm[0][2], cnorm[4][2]);
	break;

	case 10:
	normals[10][0] = lerp(DATA2, DATA6, threshold, cnorm[1][0], cnorm[5][0]);
	normals[10][1] = lerp(DATA2, DATA6, threshold, cnorm[1][1], cnorm[5][1]);
	normals[10][2] = lerp(DATA2, DATA6, threshold, cnorm[1][2], cnorm[5][2]);
	break;

	case 11:
	normals[11][0] = lerp(DATA4, DATA8, threshold, cnorm[3][0], cnorm[7][0]);
	normals[11][1] = lerp(DATA4, DATA8, threshold, cnorm[3][1], cnorm[7][1]);
	normals[11][2] = lerp(DATA4, DATA8, threshold, cnorm[3][2], cnorm[7][2]);
	break;

	case 12:
	normals[12][0] = lerp(DATA3, DATA7, threshold, cnorm[2][0], cnorm[6][0]);
	normals[12][1] = lerp(DATA3, DATA7, threshold, cnorm[2][1], cnorm[6][1]);
	normals[12][2] = lerp(DATA3, DATA7, threshold, cnorm[2][2], cnorm[6][2]);
	break;

	} /* end switch */
    } /* end for */
}

void 
get_cell_polys_and_norms(index, npolys, crossings, normals)
int index;
int *npolys;
float crossings[13][3];
float normals[13][3];
/* This subroutine will calculate the polygons */
{

    register int num_o_polys;
    register int poly;
    float *p1, *p2, *p3;
    float *n1, *n2, *n3;

    int add_polygon();		/* stores polygons for file output */

    num_o_polys = cell_table[index].npolys;
    for (poly = 0; poly < num_o_polys; poly++) {

	p1 = &crossings[cell_table[index].polys[(poly * 3)]][0];
	p2 = &crossings[cell_table[index].polys[(poly * 3) + 1]][0];
	p3 = &crossings[cell_table[index].polys[(poly * 3) + 2]][0];
	n1 = &normals[cell_table[index].polys[(poly * 3)]][0];
	n2 = &normals[cell_table[index].polys[(poly * 3) + 1]][0];
	n3 = &normals[cell_table[index].polys[(poly * 3) + 2]][0];
#define CULL
#ifdef CULL
	if ((p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2]) ||
	    (p1[0] == p3[0] && p1[1] == p3[1] && p1[2] == p3[2]) ||
	    (p2[0] == p3[0] && p2[1] == p3[1] && p2[2] == p3[2]))  {
	    (*npolys)--;
	    continue;
	}
    /*
    if (poly + *npolys == 4399 || poly + *npolys == 4400 ||
	poly + *npolys == 3801 || poly + *npolys == 3802) {
    */
    if ((*npolys <= 4399 && 4399 < num_o_polys + *npolys) ||
        (*npolys <= 4400 && 4400 < num_o_polys + *npolys) ||
        (*npolys <= 3801 && 3801 < num_o_polys + *npolys) ||
        (*npolys <= 3802 && 3802 < num_o_polys + *npolys) ) {
	printf("n %6d, poly = %3d index = %x\n", *npolys+poly, poly, index);
	printf("p1\t%g %g %g\n", p1[0], p1[1], p1[2]);
	printf("p2\t%g %g %g\n", p2[0], p2[1], p2[2]);
	printf("p3\t%g %g %g\n", p3[0], p3[1], p3[2]);
    }
#endif

	if (add_polygon(p1, p2, p3, n1, n2, n3) == -1) {
	    fprintf(stderr, "%s: error from add_polygon\n", MY_NAME);
	    exit(1);
	}
    }

    (*npolys) += num_o_polys;

}
