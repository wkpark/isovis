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
#include <math.h>
#include <malloc.h>
#include "isovis.h"


void
get_max_min(data, xdim, ydim, zdim, maxp, minp)
register float *data;
int xdim, ydim, zdim;
float *maxp, *minp;
/* This subroutine finds the maximum & minimum data values */
{
    float *enddata;
    float max, min;

    enddata = data + (xdim * ydim * zdim);
    max = min = *(data++);
#pragma ivdep
    for (; data < enddata; data++) {
	if (*data > max) max = *data;
	if (*data < min) min = *data;
    }
    *maxp = max; *minp = min;
}

void
calc_normal(p1, p2, p3, n)
float p1[3], p2[3], p3[3];
float n[3];
/* This subroutine calculates a normal from three vertices */
{
    float u[3], v[3];
    float sum, mag;

    u[0] = p3[0] - p2[0];
    u[1] = p3[1] - p2[1];
    u[2] = p3[2] - p2[2];

    v[0] = p1[0] - p2[0];
    v[1] = p1[1] - p2[1];
    v[2] = p1[2] - p2[2];

    n[0] = u[1] * v[2] - u[2] * v[1];
    n[1] = u[2] * v[0] - u[0] * v[2];
    n[2] = u[0] * v[1] - u[1] * v[0];

    sum = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
#ifdef mips
    mag = fsqrt(sum);
#else
    mag = (float) sqrt((double) sum);
#endif

    if (mag == 0.0)
	mag = 1.0;

    n[0] = n[0] / mag;
    n[1] = n[1] / mag;
    n[2] = n[2] / mag;

}

float *VERTICES = NULL;		/* a pointer to the x coordinates */
float *NORMALS = NULL;		/* a pointer to the x normals */
int NUM_VERTICES = 0;		/* number of vertices currently stored */
int VERT_LIMIT = 0;		/* currently allocated space for vertices */
int VERT_INCR = 20000;		/* allocate space for this many vertices at a
				 * time */

add_polygon(p1, p2, p3, n1, n2, n3)
float *p1, *p2, *p3, *n1, *n2, *n3;
/* This subroutine stores a polygon (triangle) in a list of vertices */
/* and connectivity.  This list can then be written out in different */
/* file formats. */
{
    unsigned size, offset;
    float *ptr;
    int vert_alloc();		/* allocates space for the vertices */

    /* see if we have enough space to store the vertices */
    if (NUM_VERTICES >= (VERT_LIMIT - 3)) {
	/* get more space */
	VERT_LIMIT += VERT_INCR;/* This is the space we need */
	size = VERT_LIMIT * sizeof(float);	/* size for malloc/realloc */
	if (vert_alloc(size) == -1) {
	    fprintf(stderr, "%s: error from vert_alloc\n", MY_NAME);
	    return -1;
	}
    }
    /* store the vertices */
    ptr = VERTICES + NUM_VERTICES * 3;
    *ptr++ = *p1++;		/* x of first vertex */
    *ptr++ = *p1++;		/* y of first vertex */
    *ptr++ = *p1++;		/* z of first vertex */
    *ptr++ = *p2++;		/* x of second vertex */
    *ptr++ = *p2++;		/* y of second vertex */
    *ptr++ = *p2++;		/* z of second vertex */
    *ptr++ = *p3++;		/* x of third vertex */
    *ptr++ = *p3++;		/* y of third vertex */
    *ptr++ = *p3++;		/* z of third vertex */

    ptr = NORMALS + NUM_VERTICES * 3;
    *ptr++ = *n1++; *ptr++ = *n1++; *ptr++ = *n1++;
    *ptr++ = *n2++; *ptr++ = *n2++; *ptr++ = *n2++;
    *ptr++ = *n3++; *ptr++ = *n3++; *ptr++ = *n3++;

    NUM_VERTICES += 3;

    return 0;
}

int 
vert_alloc(size)
int size;
/* This subroutine is for allocating memory for the vertex lists */
{

    /* flag to allocate memory from 'malloc' first time through */
    static int first_alloc = 1;

    if (first_alloc) {
	/* use 'malloc' for the first time */
	if ((VERTICES = (float *) malloc(size * 3)) == NULL) {
	    fprintf(stderr, "%s: error, not enough memory to store vertices\n",
		    MY_NAME);
	    return -1;
	}
	if ((NORMALS = (float *) malloc(size * 3)) == NULL) {
	    fprintf(stderr, "%s: error, not enough memory to store normals\n",
		    MY_NAME);
	    return -1;
	}
	first_alloc = 0;	/* use 'realloc' from now on */
    } else {
	/* use 'realloc' from now on */
	if ((VERTICES = (float *) realloc((char *) VERTICES, size * 3)) == NULL) {
	    fprintf(stderr, "%s: error, not enough memory to store vertices\n",
		    MY_NAME);
	    return -1;
	}
	if ((NORMALS = (float *) realloc((char *) NORMALS, size * 3)) == NULL) {
	    fprintf(stderr, "%s: error, not enough memory to store normals\n",
		    MY_NAME);
	    return -1;
	}
    }

    return 0;
}

void
release_memory()
{

    if (VERTICES) {
	free(VERTICES);
	VERTICES = NULL;
    }
    if (NORMALS) {
	free(NORMALS);
	NORMALS = NULL;
    }
}

dump_vset()
/* This subroutine calls a subroutine to write out a HDF VSet */
{
    int write_vset();

    if (VERBOSE)
	printf("%s: writing VSet output\n", MY_NAME);

    return (write_vset(VERTICES, NUM_VERTICES));
}

dump_wft()
/* This subroutine calls a subroutine to write out a wavefront .obj file */
{
    int write_wft();

    if (VERBOSE)
	printf("%s: writing WFT output\n", MY_NAME);

    return (write_wft(VERTICES, NUM_VERTICES, NORMALS));
}

dump_dtm()
/* This subroutine calls a subroutine to write out a wft file */
{
    int write_dtm();

    if (VERBOSE)
	printf("%s: writing DTM output\n", MY_NAME);

    return (write_dtm(VERTICES, NUM_VERTICES));
}

dump_sgi()
{
#ifdef SGI
    int write_sgi();

    return (write_sgi(VERTICES, NORMALS, NUM_VERTICES));
#endif
}

dump_byu()
/* This subroutine calls a subroutine to write out a MOVIE.BYU file */
{
    int write_byu();

    if (VERBOSE)
	printf("%s: writing BYU output\n", MY_NAME);

    return (write_byu(VERTICES, NUM_VERTICES));
}

smooth_norms() {
    printf("smoothing ...\n");
    if (VERBOSE)
	printf("%s: smoothing ...", MY_NAME);
    smooth(VERTICES, NORMALS, NUM_VERTICES);
    if (VERBOSE)
	printf("%s: done\n", MY_NAME);
}
