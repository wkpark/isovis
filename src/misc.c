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


/**************************** get_max_min ****************************/
/**************************** get_max_min ****************************/
/**************************** get_max_min ****************************/
/**************************** get_max_min ****************************/

void get_max_min(data,xdim,ydim,zdim,max,min)
register float *data;
int xdim,ydim,zdim;
float *max,*min;
/* This subroutine finds the maximum & minimum data values */
{

  float *enddata;

  enddata = data + (xdim * ydim * zdim);
  *max = *min = *(data++);
  for ( ;data<enddata;data++) {
    if (*data > *max)
       *max = *data;
    else
       if (*data < *min)
          *min = *data;
  }

}




/**************************** calc_normal ****************************/
/**************************** calc_normal ****************************/
/**************************** calc_normal ****************************/
/**************************** calc_normal ****************************/

void calc_normal(p1,p2,p3,n)
float p1[3],p2[3],p3[3];
float n[3];
/* This subroutine calculates a normal from three vertices */
{

  float u[3],v[3];
  float sum, mag;

  u[0]=p3[0] - p2[0];
  u[1]=p3[1] - p2[1];
  u[2]=p3[2] - p2[2];

  v[0]=p1[0] - p2[0];
  v[1]=p1[1] - p2[1];
  v[2]=p1[2] - p2[2];

  n[0] = u[1]*v[2] - u[2]*v[1];
  n[1] = u[2]*v[0] - u[0]*v[2];
  n[2] = u[0]*v[1] - u[1]*v[0];

  sum = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  mag = (float)sqrt((double)sum);

  if (mag == 0.0)
     mag = 1.0;

  n[0] = n[0] / mag;
  n[1] = n[1] / mag;
  n[2] = n[2] / mag;

}




/**************************** Polygon Globals ****************************/
/**************************** Polygon Globals ****************************/
/**************************** Polygon Globals ****************************/
/**************************** Polygon Globals ****************************/

float *XVERTICES=NULL;     /* a pointer to the x coordinates */
float *YVERTICES=NULL;     /* a pointer to the y coordinates */
float *ZVERTICES=NULL;     /* a pointer to the z coordinates */
int   NUM_VERTICES=0;      /* number of vertices currently stored */
int   VERT_LIMIT=0;        /* currently allocated space for vertices */
int   VERT_INCR=20000;     /* allocate space for this many vertices at a time */
int   *CONNECTIVITY=NULL;  /* a pointer to the connectivity */
int   NUM_CONN=0;          /* number of connections currently stored */
int   CONN_LIMIT=0;        /* currently allocated space for connections */
int   CONN_INCR=20000;     /* allocate space for this many connections/time */
int   CONN_SIZE=3;         /* number of elements per connection (triangle) */



/**************************** add_polygon ****************************/
/**************************** add_polygon ****************************/
/**************************** add_polygon ****************************/
/**************************** add_polygon ****************************/

int add_polygon(p1,p2,p3)
float *p1,*p2,*p3;
/* This subroutine stores a polygon (triangle) in a list of vertices */
/* and connectivity.  This list can then be written out in different */
/* file formats. */
{

  unsigned size,offset;
  float *ptr;
  int vert_alloc();  /* allocates space for the vertices */
  int conn_alloc();  /* allocates space for the connectivity */

  /* see if we have enough space to store the vertices */
  if (NUM_VERTICES >= (VERT_LIMIT-3)) {
     /* get more space */
     VERT_LIMIT += VERT_INCR; /* This is the space we need */
     size = VERT_LIMIT * sizeof(float);  /* size for malloc/realloc */
     if (vert_alloc(size) == -1) {
        fprintf(stderr,"%s: error from vert_alloc\n",MY_NAME);
        return -1;
     }
  }

  /* store the vertices */
  ptr = XVERTICES + NUM_VERTICES;
  *(ptr++) = *p1;      /* x of first vertex */
  *(ptr++) = *p2;      /* x of second vertex */
  *(ptr) = *p3;        /* x of third vertex */

  ptr = YVERTICES + NUM_VERTICES;
  *(ptr++) = *(p1+1);  /* y of first vertex */
  *(ptr++) = *(p2+1);  /* y of second vertex */
  *(ptr)   = *(p3+1);  /* y of third vertex */

  ptr = ZVERTICES + NUM_VERTICES;
  *(ptr++) = *(p1+2);  /* z of first vertex */
  *(ptr++) = *(p2+2);  /* z of second vertex */
  *(ptr)   = *(p3+2);  /* z of third vertex */


  /* see if we have enough space to store the connectivity */
  if (NUM_CONN >= (CONN_LIMIT-1)) {
     /* get more space */
     CONN_LIMIT += CONN_INCR;  /* this is the space we need */
     size = CONN_LIMIT * CONN_SIZE * sizeof(int);  /* size for malloc/realloc */
     if (conn_alloc(size) == -1) {
        fprintf(stderr,"%s: error from vert_alloc\n",MY_NAME);
        return -1;
     }
  }

  offset = NUM_CONN*CONN_SIZE;
  /* store the connectivity info */
  *(CONNECTIVITY+offset)   = NUM_VERTICES+1;
  *(CONNECTIVITY+offset+1) = NUM_VERTICES+2;
  *(CONNECTIVITY+offset+2) = NUM_VERTICES+3;

  NUM_VERTICES += 3;  /* keep track of how many vertices we have */
  NUM_CONN++;         /* increment the connectivity */

  return 0;

}




/**************************** vert_alloc ****************************/
/**************************** vert_alloc ****************************/
/**************************** vert_alloc ****************************/
/**************************** vert_alloc ****************************/

int vert_alloc(size)
int size;
/* This subroutine is for allocating memory for the vertex lists */
{

  /* flag to allocate memory from 'malloc' first time through */
  static int first_alloc=1;

  if (first_alloc) {
     /* use 'malloc' for the first time */
     if ((XVERTICES=(float *)malloc(size)) == NULL) {
        fprintf(stderr,"%s: error, not enough memory to store x vertices\n",
                MY_NAME);
        return -1;
     }
     if ((YVERTICES=(float *)malloc(size)) == NULL) {
        fprintf(stderr,"%s: error, not enough memory to store y vertices\n",
                MY_NAME);
        return -1;
     }
     if ((ZVERTICES=(float *)malloc(size)) == NULL) {
        fprintf(stderr,"%s: error, not enough memory to store z vertices\n",
                MY_NAME);
        return -1;
     }
     first_alloc=0;  /* use 'realloc' from now on */
  } else {
     /* use 'realloc' from now on */
     if ((XVERTICES=(float *)realloc((char *)XVERTICES,size)) == NULL) {
        fprintf(stderr,"%s: error, not enough memory to store x vertices\n",
                MY_NAME);
        return -1;
     }
     if ((YVERTICES=(float *)realloc((char *)YVERTICES,size)) == NULL) {
        fprintf(stderr,"%s: error, not enough memory to store y vertices\n",
                MY_NAME);
        return -1;
     }
     if ((ZVERTICES=(float *)realloc((char *)ZVERTICES,size)) == NULL) {
        fprintf(stderr,"%s: error, not enough memory to store z vertices\n",
                MY_NAME);
        return -1;
     }
  }

  return 0;

}




/**************************** conn_alloc ****************************/
/**************************** conn_alloc ****************************/
/**************************** conn_alloc ****************************/
/**************************** conn_alloc ****************************/

int conn_alloc(size)
int size;
/* This subroutine is for allocating memory for the connectivity list */
{

  /* flag to allocate memory from 'malloc' first time through */
  static int first_alloc=1;

  if (first_alloc) {
     /* use 'malloc' for the first time thru */
     if ((CONNECTIVITY=(int *)malloc(size)) == NULL) {
        fprintf(stderr,"%s: error, not enough memory to store connectivity\n",
                MY_NAME);
        return -1;
     }
     first_alloc=0;  /* turn off this flag */
  } else {
     /* use 'realloc' from now on */
     if ((CONNECTIVITY=(int *)realloc((char *)CONNECTIVITY,size)) == NULL) {
        fprintf(stderr,"%s: error, not enough memory to store connectivity\n",
                MY_NAME);
        return -1;
     }

  }

  return 0;

}




/**************************** dump_vset ****************************/
/**************************** dump_vset ****************************/
/**************************** dump_vset ****************************/
/**************************** dump_vset ****************************/

int dump_vset()
/* This subroutine calls a subroutine to write out a HDF VSet */
{

  int write_vset();

  if (VERBOSE)
     printf("%s: writing VSet output\n",MY_NAME);

  return (write_vset(XVERTICES,YVERTICES,ZVERTICES,NUM_VERTICES,
             CONNECTIVITY,NUM_CONN));


}




/**************************** dump_wft ****************************/
/**************************** dump_wft ****************************/
/**************************** dump_wft ****************************/
/**************************** dump_wft ****************************/
int dump_wft()
/* This subroutine calls a subroutine to write out a wft file */
{

  int write_wft();

  if (VERBOSE)
     printf("%s: writing WFT output\n",MY_NAME);

  return (write_wft(XVERTICES,YVERTICES,ZVERTICES,NUM_VERTICES,
             CONNECTIVITY,NUM_CONN));


}




/**************************** dump_dtm ****************************/
/**************************** dump_dtm ****************************/
/**************************** dump_dtm ****************************/
/**************************** dump_dtm ****************************/
/**************************** dump_dtm ****************************/
int dump_dtm()
/* This subroutine calls a subroutine to write out a wft file */
{

  int write_dtm();

  if (VERBOSE)
     printf("%s: writing DTM output\n",MY_NAME);

  return (write_dtm(XVERTICES,YVERTICES,ZVERTICES,NUM_VERTICES,
             CONNECTIVITY,NUM_CONN));


}

