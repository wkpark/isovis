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
#include <df.h>
#include <vg.h>
#include "isovis.h"

#ifdef DTM
  /* This is for the NCSA Data Transport Mechanism protocol */
#include <dtm.h>
#include <sdl.h>
#endif

/**************************** get_data ****************************/
/**************************** get_data ****************************/
/**************************** get_data ****************************/
/**************************** get_data ****************************/

int get_data(filename,data,xdim,ydim,zdim)
char *filename;
float **data;
int *xdim,*ydim,*zdim;
/* This subroutine will read in the data file. */
{

  int rank;
  int shape[3];
  int size;
  int ret;

  ret=DFSDgetdims(filename,&rank,shape,3);
  if (ret != 0) {
     fprintf(stderr,"%s: error from DFSDgetdims %d for %s\n",
             MY_NAME,ret,filename);
     return -1;
  }

  if (rank != 3) {
     fprintf(stderr,"%s: error, %s is rank %d; must be 3\n",
             MY_NAME,filename,rank);
     return -1;
  }

  *xdim = shape[2];  *ydim = shape[1];  *zdim = shape[0];
  if (VERBOSE)
     printf("%s: data set size xdim=%d ydim=%d zdim=%d\n",
            MY_NAME,*xdim,*ydim,*zdim);

  size = (*xdim) * (*ydim) * (*zdim) * sizeof(float);
  if ((*data = (float *)malloc(size)) == NULL) {
     fprintf(stderr,"%s: error, not enough memory for the data set\n",MY_NAME);
     return -1;
  }


  ret=DFSDgetdata(filename,rank,shape,*data);
  if (ret != 0) {
     fprintf(stderr,"%s: error from DFSDgetdata %d file %s\n",
             MY_NAME,ret,filename);
     return -1;
  }

  return 0;

}


#ifdef SGI
int write_image(image,xdim,ydim)
char *image;
int xdim,ydim;
{

  if (DF24addimage(OUTNAME,image,xdim,ydim) != 0) {
     fprintf(stderr,"%s: error, DF24addimage returned an error\n",MY_NAME);
     return -1;
  }

  if (VERBOSE)
     printf("%s: saved R24 image to %s\n",MY_NAME,OUTNAME);

  return 0;

}
#endif





/**************************** write_wft ****************************/
/**************************** write_wft ****************************/
/**************************** write_wft ****************************/
/**************************** write_wft ****************************/

int write_wft(px,py,pz,nverts,conn,nconn)
float *px,*py,*pz;
int *conn;
int nverts,nconn;
/* This subroutine writes out a wft object file from the polygons */
{

  FILE *wft_fp;
  register int i,j;

  /* open the ascii file for writing */
  if ((wft_fp=fopen(WFT_NAME,"w")) == NULL) {
     fprintf(stderr,"%s: error, can't open wft file %s\n",MY_NAME,WFT_NAME);
     return -1;
  }

  /* insert the default group 'd' */
  if (fprintf(wft_fp,"g d\n") == -1) {
     fprintf(stderr,"%s: error, can't write to wft file %s\n",MY_NAME,WFT_NAME);
     return -1;
  }

  /* write out all of the vertices */
  for (i=0;i<nverts;i++)
    if (fprintf(wft_fp,"v %f %f %f\n",*(px+i),*(py+i),*(pz+i)) == -1) {
       fprintf(stderr,"%s: error, can't write to wft file %s\n",
               MY_NAME,WFT_NAME);
       return -1;
    }

  /* write out the connectivity */
  j=0;
  for (i=0;i<nconn;i++) {
    if (fprintf(wft_fp,"fo %d %d %d\n",*(conn+j),
                *(conn+j+1),*(conn+j+2)) == -1) {
       fprintf(stderr,"%s: error, can't write to wft file %s\n",
               MY_NAME,WFT_NAME);
       return -1;
    }
    j+=3;
  }

  /* close the file */
  if (fclose(wft_fp) == -1) {
     fprintf(stderr,"%s: error, can't close wft file %s\n",MY_NAME,WFT_NAME);
     return -1;
  }

  return 0;

}




/**************************** write_vset ****************************/
/**************************** write_vset ****************************/
/**************************** write_vset ****************************/
/**************************** write_vset ****************************/

int write_vset(px,py,pz,nverts,conn,nconn)
float *px,*py,*pz;
int *conn;
int nverts,nconn;
/* This subroutine will write out a set of vertices and connectivity */
/* to an HDF file containing a single VSet for programs such as      */
/* NCSA PolyView. */
{

  DF *df_ptr;
  VGROUP *root;
  VDATA *vpx, *vpy, *vpz;
  VDATA *vplist3;

  /* open an HDF file for the VSet polygons */
  if ((df_ptr = DFopen(VSET_NAME, DFACC_ALL, 0)) == NULL) {
     fprintf(stderr,"%s: error, couldn't open %s for vset\n",
             MY_NAME,VSET_NAME);
     return -1;
  }

  /* Create a root node for the Vgroup */
  if ((root = (VGROUP *) Vattach(df_ptr, -1, "w")) == NULL) {
     fprintf(stderr,"%s: error from Vattach in write_vset\n",
             MY_NAME);
     return -1;
  }

  Vsetname(root, "/");

  /* Create the vdatas - actual vertices and connectivity */
  if ((vpx = (VDATA *) VSattach(df_ptr, -1, "w")) == NULL) {
     fprintf(stderr,"%s: error from Vattach in write_vset\n",
             MY_NAME);
     return -1;
  }
  VSsetname(vpx, "px");
  VSsetfields(vpx, "px");
  if (VSwrite(vpx, px, nverts, NO_INTERLACE) == -1) {
     fprintf(stderr,"%s: error from VSwrite in write_vset\n",
             MY_NAME);
     return -1;
  }

  if ((vpy = (VDATA *) VSattach(df_ptr, -1, "w")) == NULL) {
     fprintf(stderr,"%s: error from VSattach in write_vset\n",
             MY_NAME);
     return -1;
  }
  VSsetname(vpy, "py");
  VSsetfields(vpy, "py");
  if (VSwrite(vpy, py, nverts, NO_INTERLACE) == -1) {
     fprintf(stderr,"%s: error from VSwrite in write_vset\n",
             MY_NAME);
     return -1;
  }

  if ((vpz = (VDATA *) VSattach(df_ptr, -1, "w")) == NULL) {
     fprintf(stderr,"%s: error from VSattach in write_vset\n",
             MY_NAME);
     return -1;
  }
  VSsetname(vpz, "pz");
  VSsetfields(vpz, "pz");
  if (VSwrite(vpz, pz, nverts, NO_INTERLACE) == -1) {
     fprintf(stderr,"%s: error from VSwrite in write_vset\n",
             MY_NAME);
     return -1;
  }

  /* Create the connectivity */
  if ((vplist3 = (VDATA *) VSattach(df_ptr, -1, "w")) == NULL) {
     fprintf(stderr,"%s: error from VSattach in write_vset\n",
             MY_NAME);
     return -1;
  }
  VSfdefine(vplist3, "plist3", LOCAL_INTTYPE, 3);
  VSsetname(vplist3, "plist3");
  VSsetfields(vplist3, "plist3");
  if (VSwrite(vplist3, conn, nconn, NO_INTERLACE) == -1) {
     fprintf(stderr,"%s: error from VSwrite in write_vset\n",
             MY_NAME);
     return -1;
  }

  /* insert data in the root */
  if (Vinsert(root,vpx) == -1) {
     fprintf(stderr,"%s: error from Vinsert in write_vset\n",
             MY_NAME);
     return -1;
  }
  if (Vinsert(root,vpy) == -1) {
     fprintf(stderr,"%s: error from Vinsert in write_vset\n",
             MY_NAME);
     return -1;
  }
  if (Vinsert(root,vpz) == -1) {
     fprintf(stderr,"%s: error from Vinsert in write_vset\n",
             MY_NAME);
     return -1;
  }
  if (Vinsert(root,vplist3) == -1) {
     fprintf(stderr,"%s: error from Vinsert in write_vset\n",
             MY_NAME);
     return -1;
  }

  /* Detach Vsets */
  VSdetach(vpx);
  VSdetach(vpy);
  VSdetach(vpz);
  VSdetach(vplist3);
  Vdetach(root);

  /* Close the HDF file */
  if (DFclose(df_ptr) == -1) {
     fprintf(stderr,"%s: error from DFclose in write_vset\n",
             MY_NAME);
     return -1;
  }

  return 0;

}




/**************************** write_dtm ****************************/
/**************************** write_dtm ****************************/
/**************************** write_dtm ****************************/
/**************************** write_dtm ****************************/

int write_dtm(px,py,pz,nverts,conn,nconn)
float *px,*py,*pz;
int *conn;
int nverts,nconn;
/* This subroutine writes out a dtm object file from the polygons */
/* This subroutine will use the DTM communications protocol */
{

#ifdef DTM

  char header[SDLsize];
  float triangle[27];
  float *ptr;
  int i;
  float norm[3];

  void calc_normal();

  SDLsetClass(header);
  SDLsetPrimative(header,SDLtriangle,9);
  DTMbeginWrite(1,header,SDLsize);

  i=0;
  /* set the color at each vertex - white for now */
  triangle[21] = triangle[12] = triangle[3] = 1.0;
  triangle[22] = triangle[13] = triangle[4] = 1.0;
  triangle[23] = triangle[14] = triangle[5] = 1.0;

  while (i<nverts) {
    triangle[0] = *(px+i);
    triangle[1] = *(py+i);
    triangle[2] = *(pz+i);
    i++;

    triangle[9] = *(px+i);
    triangle[10]= *(py+i);
    triangle[11]= *(pz+i);
    i++;

    triangle[18] = *(px+i);
    triangle[19] = *(py+i);
    triangle[20] = *(pz+i);
    i++;

    /* handle the normals */
    calc_normal(&triangle[0],&triangle[9],&triangle[18],&triangle[6]);
    triangle[24] = triangle[15] = triangle[6];
    triangle[25] = triangle[16] = triangle[7];
    triangle[26] = triangle[17] = triangle[8];

    DTMsendDataset(1,triangle,27,DTM_FLOAT);
  }

  DTMend(1);

  return 0;

#endif  /* DTM */

}
