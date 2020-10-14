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
#include <stdlib.h>
#include "isovis.h"

#ifdef HDF
/* This is for the NCSA HDF format */
#include <df.h>
#endif

#ifdef VSET
/* This is for the NCSA Vset format */
#include <vg.h>
#endif

#ifdef DTM
/* This is for the NCSA Data Transport Mechanism protocol */
#include <dtm.h>
#include <sdl.h>
#endif

get_hdf_data(filename, data, xdim, ydim, zdim, max, min)
char *filename;
float **data;
int *xdim, *ydim, *zdim;
float *max, *min;
/* This subroutine will read in the HDF data file. */
{

#ifdef HDF
    int rank;
    int shape[3];
    int size;
    int ret;
    int maxmin;
    void get_max_min();		/* searches for max and min */

    ret = DFSDgetdims(filename, &rank, shape, 3);
    if (ret != 0) {
	fprintf(stderr, "%s: error from DFSDgetdims %d for %s\n",
		MY_NAME, ret, filename);
	return -1;
    }
    if (rank != 3) {
	fprintf(stderr, "%s: error, %s is rank %d; must be 3\n",
		MY_NAME, filename, rank);
	return -1;
    }
    *xdim = shape[2];
    *ydim = shape[1];
    *zdim = shape[0];
    if (VERBOSE)
	printf("%s: data set size xdim=%d ydim=%d zdim=%d\n",
	       MY_NAME, *xdim, *ydim, *zdim);

    size = (*xdim) * (*ydim) * (*zdim) * sizeof(float);
    if ((*data = (float *) malloc(size)) == NULL) {
	fprintf(stderr, "%s: error, not enough memory for the data set\n",
		MY_NAME);
	return -1;
    }
//    maxmin = DFSDgetmaxmin(max, min);
    maxmin = DFSDgetrange(max, min);
    ret = DFSDgetdata(filename, rank, shape, *data);
    if (ret != 0) {
	fprintf(stderr, "%s: error from DFSDgetdata %d file %s\n",
		MY_NAME, ret, filename);
	return -1;
    }
    if (maxmin == -1)
	get_max_min(*data, *xdim, *ydim, *zdim, max, min);
//////////
    {
// read grid dim/scale/inc
	float t[100];
    maxmin = DFSDgetdimscale(1, *xdim,&t);
	if (VERBOSE)
		fprintf(stderr, "%f,%f: dimscale\n",t[0],t[1]);
	XMIN=t[0];
	XINC=t[1]-t[0];
	if (VERBOSE)
		fprintf(stderr, "xmin=%f, xinc=%f\n",XMIN,XINC);
    maxmin = DFSDgetdimscale(2, *ydim,&t);
	if (VERBOSE)
		fprintf(stderr, "%f,%f: dimscale\n",t[0],t[1]);
	YMIN=t[0];
	YINC=t[1]-t[0];
	if (VERBOSE)
		fprintf(stderr, "ymin=%f, yinc=%f\n",YMIN,YINC);
    maxmin = DFSDgetdimscale(3, *xdim,&t);
	if (VERBOSE)
		fprintf(stderr, "%f,%f: dimscale\n",t[0],t[1]);
	ZMIN=t[0];
	ZINC=t[1]-t[0];
	if (VERBOSE)
		fprintf(stderr, "zmin=%f, zinc=%f\n",ZMIN,ZINC);
    }
/////////

    return 0;
#else
    fprintf(stderr, "%s: hdf support not installed\n");
    exit(1);
#endif

}


get_raw_data(filename, data, xdim, ydim, zdim, max, min)
char *filename;
float **data;
int xdim, ydim, zdim;
float *max, *min;
/* This subroutine will read in the raw data file. */
{

    int size;
    int ret;
    FILE *fd;
    void get_max_min();		/* searches for max and min */

    size = xdim*ydim*zdim*sizeof(float);
    if ((*data = (float *) malloc(size)) == NULL) {
	fprintf(stderr, "%s: error, not enough memory for the data set\n",
		MY_NAME);
	return -1;
    }
    if ((fd = fopen(filename, "r")) == NULL) {
	fprintf(stderr, "%s: can't open file %s\n", MY_NAME, filename);
	return -1;
    }
    if (fread(*data, size, 1, fd) != 1) {
	fprintf(stderr, "%s: short read from file %s\n", MY_NAME, filename);
	return -1;
    }
    fclose(fd);
    get_max_min(*data, xdim, ydim, zdim, max, min);
    return 0;
}


#ifdef SGI
write_image(image, xdim, ydim)
char *image;
int xdim, ydim;
{

#ifdef HDF
    if (DF24addimage(OUTNAME, image, xdim, ydim) != 0) {
	fprintf(stderr, "%s: error, DF24addimage returned an error\n", MY_NAME);
	return -1;
    }
    if (VERBOSE)
	printf("%s: saved R24 image to %s\n", MY_NAME, OUTNAME);

    return 0;
#else
    fprintf(stderr, "%s: HDF output support not installed\n", MY_NAME);
    exit(1);
#endif

}
#endif

/* Raster 3D */
write_r3d(p, nverts, n)
float *p, *n;
int nverts;
/* This subroutine writes out a raster3d file from the polygons */
{
    FILE *r3d_fp;
    register int i, j;
    float ***conn, ***tricompact(), ***tricompactn();
    int uniq;
    float **vlist;

    float scalex,scaley,scalez;

    if (XINC==0.0 && YINC==0.0 && ZINC==0.0) {
	   scalex=1.0;
	   scaley=1.0;
	   scalez=1.0;
    } else {
	   scalex=ZINC;
	   scaley=YINC;
	   scalez=XINC;
    }

    /* write out minimal number number of vertices */
    conn = NORMAL_TYPE ? tricompactn(p, n, nverts) : tricompact(p, nverts);
    uniq = conn[nverts] - conn[0]; vlist = conn[0];

    /* open the ascii file for writing */
    if ((r3d_fp = fopen(R3D_NAME, "w")) == NULL) {
	fprintf(stderr, "%s: error, can't open r3d file %s\n",
		MY_NAME, R3D_NAME);
	return -1;
    }
    /* write out all of the vertices */
    for (i = 0; i < uniq; i++) {
	float *v = vlist[i];
	v[0]*=scalex;
	v[1]*=scaley;
	v[2]*=scalez;
    }
    if (1) {
	/* write out all of the normals */
	for (i = 0; i < uniq; i++) {
	    float *v = n+(vlist[i]-p);
	    v[0]*=scaley*scalez;
	    v[1]*=scalex*scalez;
	    v[2]*=scalex*scaley;
	}
	/* write out the connectivity */
	j = 0;
	for (i = 0; i < nverts / 3; i++) {
	    int i1=conn[j]-conn[0],
		i2=conn[j+1]-conn[0],
		i3=conn[j+2]-conn[0];
	    float *v1=vlist[i1],
	          *v2=vlist[i2],
		  *v3=vlist[i3];
	    float *n1 = n+(vlist[i1]-p),
	          *n2 = n+(vlist[i2]-p),
	          *n3 = n+(vlist[i3]-p);
	    if (fprintf(r3d_fp, "# f %d %d %d\n",i1+1,i2+2,i3+3) == -1) {
	        fprintf(stderr, "%s: error, can't write to r3d file %s\n",
		    MY_NAME, R3D_NAME);
		return -1;
	    }
	    if (fprintf(r3d_fp, "1\n%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f  %4.3f %4.3f %4.3f\n",
			v1[0], v1[1], v1[2],
			v2[0], v2[1], v2[2],
			v3[0], v3[1], v3[2],
			FACE_COLOR[0],FACE_COLOR[1],FACE_COLOR[2]) == -1) {
		fprintf(stderr, "%s: error, can't write to r3d file %s\n",
			MY_NAME, R3D_NAME);
		return -1;
	    }
	    if (NORMAL_TYPE) {
	    if (fprintf(r3d_fp, "7\n%13.10f %13.10f %13.10f %13.10f %13.10f %13.10f %13.10f %13.10f %13.10f\n",
			n1[0], n1[1], n1[2],
			n2[0], n2[1], n2[2],
			n3[0], n3[1], n3[2]) == -1) {
		fprintf(stderr, "%s: error, can't write to r3d file %s\n",
			MY_NAME, R3D_NAME);
		return -1;
	    }
	    }
	    j += 3;
	}
    }

    /* close the file */
    if (fclose(r3d_fp) == -1) {
	fprintf(stderr, "%s: error, can't close r3d file %s\n",
		MY_NAME, R3D_NAME);
	return -1;
    }
    free((char *)conn[0]); free((char *)conn);
    return 0;

}

write_wft(p, nverts, n)
float *p, *n;
int nverts;
/* This subroutine writes out a wavefront object file from the polygons */
{

    FILE *wft_fp;
    register int i, j;
    float ***conn, ***tricompact(), ***tricompactn();
    int uniq;
    float **vlist;

    /* write out minimal number number of vertices */
    conn = NORMAL_TYPE ? tricompactn(p, n, nverts) : tricompact(p, nverts);
    uniq = conn[nverts] - conn[0]; vlist = conn[0];

    /* open the ascii file for writing */
    if ((wft_fp = fopen(WFT_NAME, "w")) == NULL) {
	fprintf(stderr, "%s: error, can't open wft file %s\n",
		MY_NAME, WFT_NAME);
	return -1;
    }
    /* insert the default group 'd' */
    if (fprintf(wft_fp, "g d\n") == -1) {
	fprintf(stderr, "%s: error, can't write to wft file %s\n",
		MY_NAME, WFT_NAME);
	return -1;
    }
    /* write out all of the vertices */
    for (i = 0; i < uniq; i++) {
	float *v = vlist[i];
	if (fprintf(wft_fp, "v %f %f %f\n", v[0], v[1], v[2]) == -1) {
	    fprintf(stderr, "%s: error, can't write to wft file %s\n",
		    MY_NAME, WFT_NAME);
	    return -1;
	}
    }
    if (NORMAL_TYPE) {
	/* write out all of the normals */
	for (i = 0; i < uniq; i++) {
	    float *v = n+(vlist[i]-p);
	    if (fprintf(wft_fp, "vn %f %f %f\n", v[0], v[1], v[2]) == -1) {
		fprintf(stderr, "%s: error, can't write to wft file %s\n",
			MY_NAME, WFT_NAME);
		return -1;
	    }
	}
	/* smoothing group */
	if (fprintf(wft_fp, "s 1\n") == -1) {
	    fprintf(stderr, "%s: error, can't write to wft file %s\n",
		    MY_NAME, WFT_NAME);
	}
	/* write out the connectivity */
	j = 0;
	for (i = 0; i < nverts / 3; i++) {
	    int c1=conn[j]-conn[0]+1,
		c2=conn[j+1]-conn[0]+1,
		c3=conn[j+2]-conn[0]+1;
//	    if (fprintf(wft_fp, "fo %d//%d %d//%d %d//%d\n", c1, c1, c2, c2,
	    if (fprintf(wft_fp, "f %d//%d %d//%d %d//%d\n", c1, c1, c2, c2,
			c3, c3) == -1) {
		fprintf(stderr, "%s: error, can't write to wft file %s\n",
			MY_NAME, WFT_NAME);
		return -1;
	    }
	    j += 3;
	}
    } else {

	/* write out the connectivity */
	j = 0;
	for (i = 0; i < nverts / 3; i++) {
	    int c1=conn[j]-conn[0]+1,
		c2=conn[j+1]-conn[0]+1,
		c3=conn[j+2]-conn[0]+1;
//	    if (fprintf(wft_fp, "fo %d %d %d\n", c1, c2, c3) == -1) {
	    if (fprintf(wft_fp, "f %d %d %d\n", c1, c2, c3) == -1) {
		fprintf(stderr, "%s: error, can't write to wft file %s\n",
			MY_NAME, WFT_NAME);
		return -1;
	    }
	    j += 3;
	}
    }

    /* close the file */
    if (fclose(wft_fp) == -1) {
	fprintf(stderr, "%s: error, can't close wft file %s\n",
		MY_NAME, WFT_NAME);
	return -1;
    }
    free((char *)conn[0]); free((char *)conn);
    return 0;

}

write_vset(p, nverts)
float *p;
int nverts;
/* This subroutine will write out a set of vertices and connectivity */
/* to an HDF file containing a single VSet for programs such as      */
/* NCSA PolyView. */
{

#ifdef VSET
    DF *df_ptr;
    VGROUP *root;
    VDATA *vpx, *vpy, *vpz;
    VDATA *vplist3;
    int i, *icon;
    float *pp;
    char *malloc();
    float ***conn, ***tricompact();
    int uniq;
    float **vlist;

    /* write out minimal number number of vertices */
    conn = tricompact(p, nverts);
    uniq = conn[nverts] - conn[0]; vlist = conn[0];


    if ((pp = (float *) malloc(nverts * sizeof(float))) == 0) {
	fprintf(stderr, "%s: error, could allocate buffer for vset\n", MY_NAME);
	return -1;
    }
    /* open an HDF file for the VSet polygons */
    if ((df_ptr = DFopen(VSET_NAME, DFACC_ALL, 0)) == NULL) {
	fprintf(stderr, "%s: error, couldn't open %s for vset\n",
		MY_NAME, VSET_NAME);
	return -1;
    }
    /* Create a root node for the Vgroup */
    if ((root = (VGROUP *) Vattach(df_ptr, -1, "w")) == NULL) {
	fprintf(stderr, "%s: error from Vattach in write_vset\n",
		MY_NAME);
	return -1;
    }
    Vsetname(root, "/");

    /* Create the vdatas - actual vertices and connectivity */
    if ((vpx = (VDATA *) VSattach(df_ptr, -1, "w")) == NULL) {
	fprintf(stderr, "%s: error from Vattach in write_vset\n",
		MY_NAME);
	return -1;
    }
    VSsetname(vpx, "px");
    VSsetfields(vpx, "px");
#pragma ivdep
    for (i = 0; i < uniq; i++)
	pp[i] = vlist[i][0];
    if (VSwrite(vpx, pp, uniq, NO_INTERLACE) == -1) {
	fprintf(stderr, "%s: error from VSwrite in write_vset\n",
		MY_NAME);
	return -1;
    }
    if ((vpy = (VDATA *) VSattach(df_ptr, -1, "w")) == NULL) {
	fprintf(stderr, "%s: error from VSattach in write_vset\n",
		MY_NAME);
	return -1;
    }
    VSsetname(vpy, "py");
    VSsetfields(vpy, "py");
#pragma ivdep
    for (i = 0; i < uniq; i++)
	pp[i] = vlist[i][1];
    if (VSwrite(vpy, pp, uniq, NO_INTERLACE) == -1) {
	fprintf(stderr, "%s: error from VSwrite in write_vset\n",
		MY_NAME);
	return -1;
    }
    if ((vpz = (VDATA *) VSattach(df_ptr, -1, "w")) == NULL) {
	fprintf(stderr, "%s: error from VSattach in write_vset\n",
		MY_NAME);
	return -1;
    }
    VSsetname(vpz, "pz");
    VSsetfields(vpz, "pz");
#pragma ivdep
    for (i = 0; i < uniq; i++)
	pp[i] = vlist[i][2];
    if (VSwrite(vpz, pp, uniq, NO_INTERLACE) == -1) {
	fprintf(stderr, "%s: error from VSwrite in write_vset\n",
		MY_NAME);
	return -1;
    }
    icon = (int *) pp;
    for (i = 0; i < nverts; i++)
	icon[i] = conn[i]-conn[0] + 1;

    /* Create the connectivity */
    if ((vplist3 = (VDATA *) VSattach(df_ptr, -1, "w")) == NULL) {
	fprintf(stderr, "%s: error from VSattach in write_vset\n",
		MY_NAME);
	return -1;
    }

#define LOCAL_INTTYPE 4
    VSfdefine(vplist3, "plist3", LOCAL_INTTYPE, 3);
    VSsetname(vplist3, "plist3");
    VSsetfields(vplist3, "plist3");
    if (VSwrite(vplist3, icon, nverts / 3, NO_INTERLACE) == -1) {
	fprintf(stderr, "%s: error from VSwrite in write_vset\n",
		MY_NAME);
	return -1;
    }
    free((char *) pp); free((char *)conn[0]); free((char *)conn);

    /* insert data in the root */
    if (Vinsert(root, vpx) == -1) {
	fprintf(stderr, "%s: error from Vinsert in write_vset\n",
		MY_NAME);
	return -1;
    }
    if (Vinsert(root, vpy) == -1) {
	fprintf(stderr, "%s: error from Vinsert in write_vset\n",
		MY_NAME);
	return -1;
    }
    if (Vinsert(root, vpz) == -1) {
	fprintf(stderr, "%s: error from Vinsert in write_vset\n",
		MY_NAME);
	return -1;
    }
    if (Vinsert(root, vplist3) == -1) {
	fprintf(stderr, "%s: error from Vinsert in write_vset\n",
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
	fprintf(stderr, "%s: error from DFclose in write_vset\n",
		MY_NAME);
	return -1;
    }
    return 0;
#endif

}

write_dtm(p, nverts)
float *p;
int nverts;
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
    SDLsetPrimative(header, SDLtriangle, 9);
    DTMbeginWrite(1, header, SDLsize);

    i = 0;
    /* set the color at each vertex - white for now */
    triangle[21] = triangle[12] = triangle[3] = 1.0;
    triangle[22] = triangle[13] = triangle[4] = 1.0;
    triangle[23] = triangle[14] = triangle[5] = 1.0;

    while (i < nverts) {
	triangle[0] = *p++;
	triangle[1] = *p++;
	triangle[2] = *p++;
	i++;

	triangle[9] = *p++;
	triangle[10] = *p++;
	triangle[11] = *p++;
	i++;

	triangle[18] = *p++;
	triangle[19] = *p++;
	triangle[20] = *p++;
	i++;

	/* handle the normals */
	calc_normal(&triangle[0], &triangle[9], &triangle[18], &triangle[6]);
	triangle[24] = triangle[15] = triangle[6];
	triangle[25] = triangle[16] = triangle[7];
	triangle[26] = triangle[17] = triangle[8];

	DTMsendDataset(1, triangle, 27, DTM_FLOAT);
    }

    DTMend(1);

    return 0;

#endif				/* DTM */

}

write_byu(p, nverts)
float *p;
int nverts;
/* This subroutine writes out a MOVIE.BYU object file from the polygons */
{

    FILE *byu_fp;
    register int i;
    float ***conn, ***tricompact();
    int uniq;
    float **vlist;

    /* write out minimal number number of vertices */
    conn = tricompact(p, nverts);
    uniq = conn[nverts] - conn[0]; vlist = conn[0];

    /* open the ascii file for writing */
    if ((byu_fp = fopen(BYU_NAME, "w")) == NULL) {
	fprintf(stderr, "%s: error, can't open byu file %s\n",
		MY_NAME, BYU_NAME);
	return -1;
    }
    /* write out nparts nvertices nelements nconnections */
    if (fprintf(byu_fp, "%8d%8d%8d%8d\n", 1, uniq, nverts / 3, nverts) == -1)
	goto write_error;
    /* write out beginning and ending elements of part */
    if (fprintf(byu_fp, "%8d%8d\n", 1, nverts / 3) == -1)
	goto write_error;

    /* write out all of the vertices */
    for (i = 0; i < uniq; i++) {
	float *v = vlist[i];
	if (fprintf(byu_fp, "%12.5e%12.5e%12.5e", v[0], v[1], v[2]) == -1)
	    goto write_error;
	if ((i & 1) && fputc('\n', byu_fp) == -1)
	    goto write_error;
    }
    if ((i & 1) && fputc('\n', byu_fp) == -1)
	goto write_error;

    /* write out the connectivity */
    for (i = 0; i < nverts; i++) {
	int c = conn[i]-conn[0]+1;
	if (i % 3 == 2)
	    c = -c;
	if (fprintf(byu_fp, "%8d", c) == -1)
	    goto write_error;
	if (i % 10 == 9 && fputc('\n', byu_fp) == 1)
	    goto write_error;
    }
    if (i % 10 && fputc('\n', byu_fp) == 1)
	goto write_error;

    /* close the file */
    if (fclose(byu_fp) == -1) {
	fprintf(stderr, "%s: error, can't close byu file %s\n",
		MY_NAME, BYU_NAME);
	return -1;
    }
    free((char *)conn[0]); free((char *)conn);
    return 0;

write_error:
    fprintf(stderr, "%s: error, can't write to byu file %s\n",
	    MY_NAME, BYU_NAME);
    return -1;

}
