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
#if SGI_IMAGE
#include <gl/image.h>
#endif
#endif
#include <math.h>
#include "isovis.h"


#ifdef SGI

#define TRANSPARENCY
#ifdef TRANSPARENCY
static transp;
#endif
static long backc;
static float zdist;
static int have_stencil;
#define TORADIANS(x)	((x)*M_PI/180.)

int 
init_graphics(xdim, ydim, zdim)
int xdim, ydim, zdim;
{
    int gid;

    void define_lights();
    void use_lights();

    static Matrix idmat = {1.0, 0.0, 0.0, 0.0,
			   0.0, 1.0, 0.0, 0.0,
			   0.0, 0.0, 1.0, 0.0,
			   0.0, 0.0, 0.0, 1.0};	/* identity matrix */

    if (getgdesc(GD_BITS_NORM_ZBUFFER) == 0) {
	fprintf(stderr, "%s: error, no zbuffer\n", MY_NAME);
	return -1;
    }
    if (WIN_POS[0] != -1)
	prefposition(WIN_POS[0], WIN_SIZE[0]-1 + WIN_POS[0],
		     WIN_POS[1], WIN_SIZE[1]-1 + WIN_POS[1]);
    else
	keepaspect(1, 1);

    if (WIN_FOREGROUND)
	foreground();
    gid = winopen(MY_NAME);
    subpixel(TRUE);
    if (gid == -1) {
	fprintf(stderr, "%s: error, can't open a graphics window\n",
		MY_NAME);
	return -1;
    }
    RGBmode();

    if (have_stencil = getgdesc(GD_BITS_STENCIL) >= 1)
	/* for hollow fill */
	stensize(1);
    doublebuffer();
    gconfig();
    frontbuffer(1);

    zbuffer(1);
    /*
     * lsetdepth(0x0000,0xfffff); 
     */
    /* use czclear() */
    lsetdepth(0x7fffff, 0);
    zfunction(ZF_GEQUAL);

    mmode(MVIEWING);
    if (WIN_POS[0] != -1)
	perspective(450, ((float)WIN_SIZE[0])/WIN_SIZE[1], 0.1, 1000.);
    else
	perspective(450, 1.0, 0.1, 1000.);
    loadmatrix(idmat);

#define max(a,b) a > b ? a : b;
    zdist = max(xdim,ydim);
    zdist = max(zdist,zdim);
#undef max
    zdist = zdist*.5 + M_SQRT2*zdist*cos(TORADIANS(45.))/sin(TORADIANS(45.));
    /*
    lookat((float) (xdim / 2), (float) (ydim / 2), zdist, 0.0, 0.0, 0.0, 0);
    */
    lookat(0., 0., zdist, 0.0, 0.0, 0.0, 0);
    define_lights();
    pushmatrix();

    if (NORMAL_TYPE == 0 && !SMOOTH)
	shademodel(FLAT);

    backc = ((int) (BACKG_COLOR[2] * 255.0) << 16) |
	    ((int) (BACKG_COLOR[1] * 255.0) << 8) |
	     (int) (BACKG_COLOR[0] * 255.0);
    czclear(backc, 0);


    translate(TRANSLATION[0], TRANSLATION[1], TRANSLATION[2]);

    /* can't use 'rot' due to bug in DGL for the Cray */
    rotate(((int) ROTATION[0] * 10), 'x');
    rotate(((int) ROTATION[1] * 10), 'y');
    rotate(((int) ROTATION[2] * 10), 'z');

    scale(SCALE[0], SCALE[1], SCALE[2]);

    if (VERBOSE) {
	printf("%s: scaling %f %f %f\n", MY_NAME, SCALE[0],
	       SCALE[1], SCALE[2]);
	printf("%s: rotation %f %f %f\n", MY_NAME, ROTATION[0],
	       ROTATION[1], ROTATION[2]);
	printf("%s: translation %f %f %f\n", MY_NAME, TRANSLATION[0],
	       TRANSLATION[1], TRANSLATION[2]);
    }
    return 0;

}

void 
define_lights()
{

    float light[10];
    float material[17];

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
    material[1] = K_VALUES[2] * MAT_COLOR[0];
    material[2] = K_VALUES[2] * MAT_COLOR[1];
    material[3] = K_VALUES[2] * MAT_COLOR[2];
    material[4] = AMBIENT;
    material[5] = K_VALUES[0] * MAT_COLOR[0];
    material[6] = K_VALUES[0] * MAT_COLOR[1];
    material[7] = K_VALUES[0] * MAT_COLOR[2];
    material[8] = DIFFUSE;
    material[9] = K_VALUES[1] * MAT_COLOR[0];
    material[10] = K_VALUES[1] * MAT_COLOR[1];
    material[11] = K_VALUES[1] * MAT_COLOR[2];
    material[12] = SHININESS;
    material[13] = K_EXP;
    material[14] = ALPHA;
#ifdef TRANSPARENCY
    if (transp)
	material[15] = .5;
    else
	material[15] = 1.;
#else
    material[15] = 1.;
#endif
    material[16] = LMNULL;

    lmdef(DEFMATERIAL, 1, 15, material);
    lmdef(DEFLIGHT, 1, 10, light);
    lmdef(DEFLMODEL, 1, 0, NULL);

    lmbind(MATERIAL, 1);
    lmbind(LIGHT0, 1);
    lmbind(LMODEL, 1);
    light[1] = -LIGHT_LOC[0];
    light[2] = -LIGHT_LOC[1];
    light[3] = -LIGHT_LOC[2];
    lmdef(DEFLIGHT, 2, 10, light);
    if (AUXLIGHT)
	lmbind(LIGHT1, 2);

}

int 
get_image(image, xdim, ydim)
char **image;
int *xdim, *ydim;
{
    long xsize, ysize;
    register int x, y, i;
    register long *sgi_image;
    long ret;
    char *malloc();

    getsize(&xsize, &ysize);
    *xdim = (int) xsize; *ydim = (int) ysize;

    if ((sgi_image = (long *) malloc(xsize*ysize*sizeof(long))) == NULL) {
	fprintf(stderr, "%s: error, not enough memory for sgi_image\n", MY_NAME);
	return -1;
    }
    ret = lrectread(0, 0, (xsize - 1), (ysize - 1), sgi_image);
    if (ret != (xsize * ysize)) {
	fprintf(stderr, "%s: error from lrectread %ld\n", MY_NAME, ret);
	return -1;
    }
    if ((*image = (char *) malloc(xsize*ysize*3)) == NULL) {
	fprintf(stderr, "%s: error, not enough memory for image\n", MY_NAME);
	return -1;
    }
    i = 0;
    for (y = ysize-1; y >= 0; y--) {
	for (x = 0; x < (*xdim); x++) {
	    long s = sgi_image[y*xsize+x];
	    *(*image + i) =  s & 0xff;
	    i++;
	    *(*image + i) =  (s >> 8) & 0xff;
	    i++;
	    *(*image + i) =  (s >> 16) & 0xff;
	    i++;
	}
    }
    free((char *)sgi_image);

    return 0;

}

#ifdef SGI_IMAGE
write_sgi_image(name)
char *name;
{
    long xsize, ysize;
    int x, y;
    long *sgi_image;
    short *pbuf;
    IMAGE *image;
    char *malloc();

    getsize(&xsize, &ysize);

    image = iopen(name,"w",RLE(1),3,xsize,ysize,3);

    pbuf = (short *)malloc(xsize*sizeof(short)*3);
    sgi_image = (long *) malloc(xsize*sizeof(long));
    if (!pbuf || !sgi_image) {
	fprintf(stderr, "%s: error, not enough memory for write_sgi_image\n", MY_NAME);
	return -1;
    }
    for (y = 0; y < ysize; y++) {
	lrectread(0, y, xsize-1, y, sgi_image);
	for (x = 0; x < xsize; x++) {
	    long s = sgi_image[x];
	    pbuf[x] =  s & 0xff;
	    pbuf[x+xsize] =  (s >> 8) & 0xff;
	    pbuf[x+2*xsize] =  (s >> 16) & 0xff;
	}
	putrow(image,pbuf,y,0);
	putrow(image,pbuf+xsize,y,1);
	putrow(image,pbuf+2*xsize,y,2);
    }
    free((char *)sgi_image); free((char *)pbuf);
    iclose(image);

    return 0;
}
#endif


int rot_which;
int cur_rot[3];
int trans_which;
int rottrans;
float cur_trans[3];
int shade;
int show_normals;
long menu, rmenu, smenu;
int light = 1;
#ifdef ANTIALIAS
static antialias;
#endif

void 
window_event_loop(data, xdim, ydim, zdim, threshold)
float *data;			/* the volume of data */
int xdim, ydim, zdim;		/* the dimensions of the data */
float threshold;		/* the current threshold of interest */
/*
 * This subroutine is called by main and implements a simple event loop that
 * handles screen refreshes and <esc> key for exiting. 
 */
{
    short value;
    int dev;

    if (AUXLIGHT) light += 2;

    /* Queue devices for <esc> key as an exit, and redraw events */
    qdevice(ESCKEY);
    qdevice(REDRAW);
    qdevice(RIGHTMOUSE);
    qdevice(LEFTMOUSE);
    qdevice(MIDDLEMOUSE);
#ifdef MYPUP
    rmenu = defpup("position %t%s|rotate x %I%r0%x1|rotate y %r0%x2|rotate z %r0%x3|translate x %r0%x4|translate y %r0%x5|translate z %r0%x6");
#if defined(ANTIALIAS) || defined(TRANSPARENCY)
    if (getgdesc(GD_BLEND))
	smenu = defpup("shade %t%s|flat %r0%x7|smooth %I%r0%x8|wire %r0%x9|hidden wire %r0%x10|off %r0%l%x11|normals %i%l%x12|antialias %i%x13|transparent %i%x14");
    else
	smenu = defpup("shade %t%s|flat %r0%x7|smooth %I%r0%x8|wire %r0%x9|hidden wire %r0%x10|off %r0%l%x11|normals %i%x12");
#else
    smenu = defpup("shade %t%s|flat %r0%x7|smooth %I%r0%x8|wire %r0%x9|hidden wire %r0%x10|off %r0%l%x11|normals %i%x12");
#endif
    menu = defpup("isovis %t%s|shading %m|position %m|bbox %i%x15|reset %x16|+ light %I%x17|- light %i%x18|exit %x100", smenu, rmenu);
    if (light&2) changepup(menu,6,"- light %I%x18");
    if (BBOX) changepup(menu,3,"bbox %I%x15");
#else /* !MYPUP */
    rmenu = defpup("position %t|rotate x %x1|rotate y %x2|rotate z %x3|translate x %x4|translate y %x5|translate z %x6");
#if defined(ANTIALIAS) || defined(TRANSPARENCY)
    if (getgdesc(GD_BLEND))
	smenu = defpup("shade %t|flat %x7|smooth %x8|wire %x9|hidden wire %x10|off %l%x11|normals %l%x12|antialias %x13|transparent %x14");
    else
	smenu = defpup("shade %t|flat %x7|smooth %x8|wire %x9|hidden wire %x10|off %l%x11|normals %x12");
#else
    smenu = defpup("shade %t|flat %x7|smooth %x8|wire %x9|hidden wire %x10|off %l%x11|normals %x12");
#endif
    menu = defpup("isovis %t|shading %m|position %m|bbox %x15|reset %x16|+ light %x17|- light %x18|exit %x100", smenu, rmenu);
#endif
    if (NORMAL_TYPE) shade = 1;  /* shade also may have been set in write_sgi */

#ifdef MYPUP
    if (!shade) {
	changepup(smenu, 1, "flat %r0%x7%I");
	changepup(smenu, 2, "smooth %r0%x8");
    }
#endif
    if (NTSC) linewidth(2);
    /* do a simple event loop */
#define DISPLAY_LIST
#ifdef DISPLAY_LIST
    /* try to build the triangle mesh before the user starts to
     * interact with the program so he won't notice a big delay
     */
    if (NORMAL_TYPE)
	make_tmesh();
#endif
    cur_rot[0] = ROTATION[0]*10;
    cur_rot[1] = ROTATION[1]*10;
    cur_rot[2] = ROTATION[2]*10;
    cur_trans[0] = TRANSLATION[0];
    cur_trans[1] = TRANSLATION[1];
    cur_trans[2] = TRANSLATION[2];
    while (1) {
	dev = qread(&value);	/* get the event */
	if (dev == ESCKEY)	/* <esc> key was pressed */
	    exit(0);
	else if (dev == REDRAW) {	/* window needs redrawing */
	    reshapeviewport();
	    redraw(BBOX, xdim, ydim, zdim);
	} else if (dev == RIGHTMOUSE) {
	    int m;
	    switch (m = dopup(menu)) {
	    case 1: case 2: case 3:
		rot_which = m - 1; rottrans = 0; break;
	    case 4: case 5: case 6:
		trans_which = m - 4; rottrans = 1; break;
	    case 7: case 8: case 9: case 10: case 11:
		shade = m - 7;
		if (shade < 2) shademodel(shade ? GOURAUD : FLAT);
		redraw(BBOX, xdim, ydim, zdim);
		break;
	    case 12: show_normals ^= 1; redraw(BBOX, xdim, ydim, zdim); break;
	    case 13: antialias ^= 1; redraw(BBOX, xdim, ydim, zdim); break;
	    case 14: transp ^= 1; popmatrix(); pushmatrix();
		     define_lights(); redraw(BBOX, xdim, ydim, zdim); break;
	    case 15: BBOX ^= 1; redraw(BBOX, xdim, ydim, zdim); break;
	    case 16:
		cur_rot[0] = ROTATION[0]*10;
		cur_rot[1] = ROTATION[1]*10;
		cur_rot[2] = ROTATION[2]*10;
		cur_trans[0] = TRANSLATION[0];
		cur_trans[1] = TRANSLATION[1];
		cur_trans[2] = TRANSLATION[2];
		redraw(BBOX, xdim, ydim, zdim);
		break;
	    case 17: case 18:
		light ^= (m-16); if (!light) light = 3 - (m-16);
#ifdef MYPUP
		changepup(menu, 5, ((light&1) ? "+ light %I%x17" : "+ light %i%x15"), 0);
		changepup(menu, 6, ((light&2) ? "- light %I%x18" : "- light %i%x16"), 0);
#endif
		/* set tranformations for light definition */
		popmatrix(); pushmatrix();
		lmbind(LIGHT0, ((light&1) ? 1 : 0));
		lmbind(LIGHT1, ((light&2) ? 2 : 0));
		redraw(BBOX, xdim, ydim, zdim); break;
	    case 100: exit(0);
	    default:
		break;
	    }
	} else if (dev == LEFTMOUSE) {
	    frontbuffer(0);
	    while (getbutton(LEFTMOUSE)) {
		if (!rottrans) {
		    cur_rot[rot_which] += 50;
		    if (cur_rot[rot_which] > 3600) cur_rot[rot_which] = 0;
		} else {
		    cur_trans[trans_which] +=.05 * xdim;
		}
		redraw(BBOX, xdim, ydim, zdim);
		swapbuffers();
	    }
	    frontbuffer(1);
	} else if (dev == MIDDLEMOUSE) {
	    frontbuffer(0);
	    while (getbutton(MIDDLEMOUSE)) {
		if (!rottrans) {
		    cur_rot[rot_which] -= 50;
		    if (cur_rot[rot_which] < 0) cur_rot[rot_which] = 3600;
		} else {
		    cur_trans[trans_which] -=.05 * xdim;
		}
		redraw(BBOX, xdim, ydim, zdim);
		swapbuffers();
	    }
	    frontbuffer(1);
	}
    }	/* while true */
}

#ifdef DISPLAY_LIST
static float *vlist, *nlist;
static int vsize;
#endif

redraw(bbox, xdim, ydim, zdim)
{
    struct { float x, y, z; } v;
    czclear(backc, 0);
    popmatrix(); pushmatrix();
    translate(cur_trans[0], cur_trans[1], cur_trans[2]);
    rotate(cur_rot[0], 'x');
    rotate(cur_rot[1], 'y');
    rotate(cur_rot[2], 'z');
    scale(SCALE[0], SCALE[1], SCALE[2]);
    if (bbox) {
#ifdef ANTIALIAS
	if (antialias) {
	    linesmooth(SML_ON|SML_SMOOTHER|SML_END_CORRECT);
	    blendfunction(BF_SA, BF_MSA);
	}
#endif
	cpack(0xFFFFFF00);
	bgnclosedline();
	v.x = -xdim/2; v.y = -ydim/2; v.z = -zdim/2;
	v3f(&v);
	v.x = xdim/2; v3f(&v);
	v.y = ydim/2; v3f(&v);
	v.x = -xdim/2; v3f(&v);
	v.y = -ydim/2; /*v3f(&v);*/
	endclosedline();
	bgnclosedline();
	v.z = zdim/2; v3f(&v);
	v.x = xdim/2; v3f(&v);
	v.y = ydim/2; v3f(&v);
	v.x = -xdim/2; v3f(&v);
	v.y = -ydim/2; /*v3f(&v);*/
	endclosedline();
	bgnline(); v3f(&v); v.z = -zdim/2; v3f(&v); endline();
	v.x = xdim/2;
	bgnline(); v3f(&v); v.z = zdim/2; v3f(&v); endline();
	v.y = ydim/2;
	bgnline(); v3f(&v); v.z = -zdim/2; v3f(&v); endline();
	v.x = -xdim/2;
	bgnline(); v3f(&v); v.z = zdim/2; v3f(&v); endline();
#ifdef ANTIALIAS
	if (antialias) {
	    blendfunction(BF_ONE, BF_ZERO);
	    linesmooth(SML_OFF);
	}
#endif
    }
#ifdef TRANSPARENCY
	if (transp) {
	    blendfunction(BF_SA, BF_MSA);
	    zwritemask(0);
	}
#endif
#ifdef DISPLAY_LIST
    if (NORMAL_TYPE)
	display(vlist,nlist,vsize);
    else
	dump_sgi();
#else
    dump_sgi();
#endif
#ifdef TRANSPARENCY
	if (transp) {
	    blendfunction(BF_ONE, BF_ZERO);
	    zwritemask(0xffffffff);
	}
#endif
}

draw_normals(p, n, nv)
float *p, *n;
{
    struct { float x, y, z; } nn;
    int i;

    cpack(0xFF00FFFF);
    for (i = 0; i < nv; i++) {
	nn.x = *n++ + *p++;
	nn.y = *n++ + *p++;
	nn.z = *n++ + *p++;
	bgnline(); v3f(p - 3); v3f(&nn); endline();
    }
}

#ifdef DISPLAY_LIST
static tmesh_obj;

make_tmesh() {
    if (VERBOSE) printf("tmeshing ... ");
    tmesh_obj = tmesh(vlist,nlist,vsize, (NORMAL_TYPE ? 1 : 0));
    if (VERBOSE) printf("done\n");
}

display(p, n, nvert)
float *p;
float *n;
{
    static wire_obj;
    static norm_obj;

    if (!tmesh_obj && shade != 2) {
/*
#define SMOOTH
#ifdef SMOOTH
	if (!NORMAL_TYPE) {
	    printf("smoothing ... ");
	    smooth(p, n, nvert);
	    printf("done\n");
	    shade = 1;
	    shademodel(GOURAUD);
	}
#endif
*/
	if (VERBOSE) printf("tmeshing ... ");
	tmesh_obj = tmesh(p,n,nvert, (NORMAL_TYPE ? 1 : 0));
	if (VERBOSE) printf("done\n");
    }
    if (shade == 2 || shade == 3) {
#ifdef POLYMODE
	if (have_stencil) {
	    if (shade == 2) {
		lmbind(MATERIAL, 0);
		cpack(0xFF000000 |
		      ((int) (255 * MAT_COLOR[2]) << 16) |
		      ((int) (255 * MAT_COLOR[1]) << 8) |
		       (int) (255 * MAT_COLOR[0]));
		polymode(PYM_LINE);
		callobj(tmesh_obj);
		polymode(PYM_FILL);
		lmbind(MATERIAL, 1);
	    } else {
		lmbind(MATERIAL, 0);
		cpack(backc);
		callobj(tmesh_obj);
		/* draw the lines as hollow polygons */
		polymode(PYM_HOLLOW);
		stencil(TRUE,1,SF_EQUAL,1,ST_KEEP,ST_KEEP,ST_KEEP);
		cpack(0xFF000000 |
		      ((int) (255 * MAT_COLOR[2]) << 16) |
		      ((int) (255 * MAT_COLOR[1]) << 8) |
		       (int) (255 * MAT_COLOR[0]));
		callobj(tmesh_obj);
		polymode(PYM_FILL);
		stencil(FALSE,0,0,0,0,0,0);
		lmbind(MATERIAL, 1);
	    }
	    goto out;
	}
#endif
	if (shade == 3) {
	    lmbind(MATERIAL, 0);
	    cpack(backc);
	    callobj(tmesh_obj);
	    lmbind(MATERIAL, 1);
	}
	if (!wire_obj) {
	    int i;
	    if (VERBOSE) printf("creating wire obj\n");
#ifdef FOO
	    makeobj(wire_obj = genobj());
	    for (i = 0; i < nvert / 3; i++, p += 9) {
		bgnclosedline();
		v3f(p); v3f(p + 3); v3f(p + 6);
		endclosedline();
	    }
	    closeobj();
	    p -= nvert*3;
#else
	    wire_obj = trilines(p, nvert);
#endif
	}
	cpack(0xFF000000 |
	      ((int) (255 * MAT_COLOR[2]) << 16) |
	      ((int) (255 * MAT_COLOR[1]) << 8) |
	       (int) (255 * MAT_COLOR[0]));
#ifdef ANTIALIAS
	if (antialias) {
	    linesmooth(SML_ON|SML_SMOOTHER|SML_END_CORRECT);
	    blendfunction(BF_SA, BF_MSA);
	    callobj(wire_obj);
	    blendfunction(BF_ONE, BF_ZERO);
	    linesmooth(SML_OFF);
	} else
	    callobj(wire_obj);
#else
	callobj(wire_obj);
#endif
#ifdef POLYMODE
out:;
#endif
    } else if (shade < 2) {
	callobj(tmesh_obj);
    }
    if (show_normals) {
	if (!norm_obj) {
	    if(VERBOSE) printf("creating normal obj\n");
#ifdef FOO
	    makeobj(norm_obj = genobj());
	    draw_normals(p, n, nvert);
	    closeobj();
#else
	    norm_obj = trinorms(p, n, nvert, 1.0);
#endif
	}
	cpack(0xFF00FFFF);
#ifdef ANTIALIAS
	if (antialias) {
	    linesmooth(SML_ON|SML_SMOOTHER|SML_END_CORRECT);
	    blendfunction(BF_SA, BF_MSA);
	    callobj(norm_obj);
	    blendfunction(BF_ONE, BF_ZERO);
	    linesmooth(SML_OFF);
	} else
	    callobj(norm_obj);
#else
	callobj(norm_obj);
#endif
    }
}
#endif

write_sgi(p, n, nvert)
float *p;
float *n;
{
    int i;
/*
#define SMOOTH
#ifdef SMOOTH
    static smoothed;
    if (!NORMAL_TYPE && !smoothed) {
	printf("smoothing ... ");
	smooth(p, n, nvert);
	printf("done\n");
	shade = 1;
	shademodel(GOURAUD);
	smoothed++;
    }
#endif
*/
    if (BBOX) {
	struct { float x, y, z; } v;
	cpack(0xFFFFFF00);
#ifdef ANTIALIAS
	if (antialias) {
	    linesmooth(SML_ON|SML_SMOOTHER|SML_END_CORRECT);
	    blendfunction(BF_SA, BF_MSA);
	}
#endif
	bgnclosedline();
	v.x = -XDIM/2; v.y = -YDIM/2; v.z = -ZDIM/2;
	v3f(&v);
	v.x = XDIM/2; v3f(&v);
	v.y = YDIM/2; v3f(&v);
	v.x = -XDIM/2; v3f(&v);
	v.y = -YDIM/2; /*v3f(&v);*/
	endclosedline();
	bgnclosedline();
	v.z = ZDIM/2; v3f(&v);
	v.x = XDIM/2; v3f(&v);
	v.y = YDIM/2; v3f(&v);
	v.x = -XDIM/2; v3f(&v);
	v.y = -YDIM/2; /*v3f(&v);*/
	endclosedline();
	bgnline(); v3f(&v); v.z = -ZDIM/2; v3f(&v); endline();
	v.x = XDIM/2;
	bgnline(); v3f(&v); v.z = ZDIM/2; v3f(&v); endline();
	v.y = YDIM/2;
	bgnline(); v3f(&v); v.z = -ZDIM/2; v3f(&v); endline();
	v.x = -XDIM/2;
	bgnline(); v3f(&v); v.z = ZDIM/2; v3f(&v); endline();
#ifdef ANTIALIAS
	if (antialias) {
	    blendfunction(BF_ONE, BF_ZERO);
	    linesmooth(SML_OFF);
	}
#endif
    }
    vlist = p; nlist = n; vsize = nvert;
    if (shade == 2 || shade == 3) {
	if (shade == 3) {
	    cpack(backc);
	    for (i = 0; i < nvert / 3; i++, p += 9) {
#define TMESH
#ifdef TMESH
		bgntmesh();
		v3f(p);
		v3f(p + 3);
		v3f(p + 6);
		endtmesh();
#else
		bgnpolygon();
		v3f(p);
		v3f(p + 3);
		v3f(p + 6);
		endpolygon();
#endif
	    }
	    p -= nvert*3;
	}
	cpack(0xFF000000 |
	      ((int) (255 * MAT_COLOR[2]) << 16) |
	      ((int) (255 * MAT_COLOR[1]) << 8) |
	       (int) (255 * MAT_COLOR[0]));
	for (i = 0; i < nvert / 3; i++, p += 9) {
	    bgnclosedline();
	    v3f(p);
	    v3f(p + 3);
	    v3f(p + 6);
	    v3f(p);
	    endclosedline();
	}
	p -= nvert*3;
    } else if (shade < 2) {
#ifdef SMOOTH
	if (NORMAL_TYPE || smoothed) {
#else
	if (NORMAL_TYPE) {
#endif
	    for (i = 0; i < nvert / 3; i++, p += 9, n += 9) {
#define TMESH
#ifdef TMESH
		bgntmesh();
		n3f(n);
		v3f(p);
		n3f(n + 3);
		v3f(p + 3);
		n3f(n + 6);
		v3f(p + 6);
		endtmesh();
#else
		bgnpolygon();
		n3f(n);
		v3f(p);
		n3f(n + 3);
		v3f(p + 3);
		n3f(n + 6);
		v3f(p + 6);
		endpolygon();
#endif
	    }
	} else {
	    for (i = 0; i < nvert / 3; i++, p += 9, n += 9) {
		bgnpolygon();
		n3f(n);
		v3f(p);
		v3f(p + 3);
		v3f(p + 6);
		endpolygon();
	    }
	}
	p -= nvert*3; n -= nvert*3;
    }
    if (show_normals)
	draw_normals(p, n, nvert);
}

#endif	/* SGI */
