/*
 *  lines.c
 *
 *  Given a sequence of lines, generate a display list which
 *  uses long sequences of bgnline()/endline() and which eliminates
 *  redundant lines (so the display list is fast).  Useful for displaying
 *  wireframe meshes from the polygonal description.
 */

#include <stdio.h>
#ifdef SGI
#include <gl/gl.h>

#define VERBOSE(x)	x

struct hash {
    float *v1, *v2;
    struct hash *next;
    struct hash *vnext;
};

#define TABLE_SIZE 9967

static struct hash *htable[TABLE_SIZE];

char *malloc();

static
hashit(v)
float *v;
{
   unsigned long *p = (unsigned long *)v;
   return (((p[0]*283+p[1])*283)+p[2]) % TABLE_SIZE;
}

static
#ifdef ndef
hinsert(v, l)
float *v;
#else
hinsert(v1, v2)
float *v1, *v2;
#endif
{
#ifdef ndef
    int h = hashit(&v[6*l]);
#else
    int h = hashit(v1);
#endif
    struct hash *hh, *hhh, *newhash();

    /* if line already exists - ignore */
    if (htable[h]) {
	for(hh = htable[h]; hh; hh=hh->next) {
	    if (hh->v1[0] == v1[0] &&
		hh->v1[1] == v1[1] &&
		hh->v1[2] == v1[2]) {
		struct hash *vv;
		for(vv = hh; vv; vv = vv->vnext)
		    if (vv->v2[0] == v2[0] &&
			vv->v2[1] == v2[1] &&
			vv->v2[2] == v2[2]) return 1;
	       goto insert_1;
	    }
	}
    }
#ifdef MALLOC
    hhh = (struct hash *)malloc(2*sizeof *hh);
#else
    hhh = newhash();
#endif
    hhh->next = htable[h]; htable[h] = hhh;
    hhh->vnext = 0;
    hhh->v1 = v1;
    goto other;

insert_1:
#ifdef MALLOC
    hhh = (struct hash *)malloc(2*sizeof *hh);
#else
    hhh = newhash();
#endif
    hhh->next = hh->next;
    hhh->vnext = hh->vnext; hh->vnext = hhh; 
    hhh->v1 = hh->v1;

other:
    hhh++; h = hashit(v2);
    for(hh = htable[h]; hh; hh = hh->next) {
	if (hh->v1[0] == v2[0] &&
	    hh->v1[1] == v2[1] &&
	    hh->v1[2] == v2[2]) goto insert_2;
    }
    hhh->next = htable[h]; htable[h] = hhh;
    hhh->vnext = 0;
    hhh->v1 = v2; hhh[-1].v2 = hhh->v1; hhh->v2 = hhh[-1].v1;
    return 0;

insert_2:
    hhh->next = hh->next;
    hhh->vnext = hh->vnext; hh->vnext = hhh; 
    hhh->v1 = hh->v1;  hhh[-1].v2 = hhh->v1; hhh->v2 = hhh[-1].v1;
    return 0;
}

line_mesh(v, nv)
float **v;
{
    int i;
    int np = 0, vp = 0, dup = 0;
    int obj;
    int nlines;
    int current = 0;
    int nbg = 0, broken = 0, singlets = 0;

    /* build hash table of unique lines */
    for(i = 0; i < nv/2; i++)
	dup += hinsert(v[2*i], v[2*i+1]);

    VERBOSE(printf ("dups = %d\n", dup));

    nlines = nv/2 - dup;

    for(i = 0; i < TABLE_SIZE; i++) {
	struct hash *p;
	if (!htable[i]) continue;
	for(p = htable[i]; p; p = p->next) {
	    vp++;
	}
	np++;
    }
    printf("vp %d np %d  %f\n", vp, np, (float)vp/np);

#ifdef DEBUG
for(i = 0; i < TABLE_SIZE; i++) {
    struct hash *p;
    if (!htable[i]) continue;
    printf("%d\n", i);
    for(p = htable[i]; p; p = p->next) {
	struct hash *pp;
	for(pp = p; pp; pp = pp->vnext)
	    printf("\t%x %x  %x\n", pp->v1, pp->v2, pp->next);
	printf("\n");
    }
}
#endif

    makeobj(obj = genobj());

    /* thread out a display list */
    while(nlines) {
	struct hash *p, *p1, *p2;
	int h;

	/* find a seed vertex consisting of an end vertex */
	for(i = current; i < TABLE_SIZE; i++) {
	    if (!htable[i]) continue;
	    for(p1 = 0, p = htable[i]; p; p1 = p, p = p->next)
		if (!p->vnext) goto found_one;
	}
	for(i = 0; i < current; i++) {
	    if (!htable[i]) continue;
	    for(p1 = 0, p = htable[i]; p; p1 = p, p = p->next)
		if (!p->vnext) goto found_one;
	}
	/* take any vertex */
	for(p1 = 0, i = 0; i < TABLE_SIZE; i++)
	    if (htable[i]) break;
	p = htable[i];
	broken++;

found_one:
	current = i /*== TABLE_SIZE-1 ? 0 : i+1*/;

	/* remove from hash */
	if (!p1) {
	    /* beginning of list */
	    htable[i] = p->vnext ? p->vnext : p->next;
	} else {
	    struct hash *pp, *next;
	    next = p->vnext ? p->vnext : p->next;
	    for(pp = p1; pp; pp = pp->vnext) pp->next = next;
	}

	bgnline();
	v3f(p->v1);

	for(i = 1;; ) {
	    struct hash *p3, *p4;

	    /* find his partner */
	    h = hashit(p->v2);
	    for(p2 = 0, p1 = htable[h]; p1; p2 = p1, p1 = p1->next)
		if (p1->v1 == p->v2) break;
	    for(p4 = 0, p3 = p1; p3; p4 = p3, p3 = p3->vnext)
		if (p3->v2 == p->v1) break;

	    /* remove from sublist */
	    if (!p4) {
		p3 = p1->vnext;
		p4 = p1->vnext ? p1->vnext : p1->next;
	    } else {
		p4->vnext = p3->vnext;
		p3 = p4 = p1;
	    }
	    /* adjust hash table */
	    if (!p2) { /* beginning of list */
		htable[h] = p4;
	    } else {
		struct hash *pp, *next;
		for(pp = p2; pp; pp = pp->vnext) pp->next = p4;
	    }
	    v3f(p->v2);

	    i++; nlines--;

	    /* new candidate ? */

	    if (!(p = p3)) {
		/*
		if (p1->off) p1--; free((char*)p1);
		*/
		break; /* i lose */
	    }
	    /*
	    if (p1->off) p1--; free((char*)p1);
	    */

	    /* remove him */
	    if (!p2) { /* beginning of list */
		htable[h] = p->vnext ? p->vnext : p->next;
	    } else {
		struct hash *pp, *next;
		next = p->vnext ? p->vnext : p->next;
		for(pp = p2; pp; pp = pp->vnext) pp->next = next;
	    }
	}
	endline(); nbg++;
#ifdef DEBUG
	VERBOSE(printf("i = %d\n", i));
#endif
	if (i == 2) singlets++;
    }
    freehash();
    printf("nbg = %d %% %f  broken %d singlets %d\n", nbg, (nv/2-dup)/(float)nbg, broken, singlets);

    for(i = 0; i < TABLE_SIZE; i++) {
	if (htable[i]) printf("hash botch - entry %d not zero\n", i);
    }

    closeobj();

    return obj;
}

/*
 *  take a list of triangle vertices, convert them to line segments
 *  (3 per triangle) and convert them to a mesh of lines, with redundant
 *  lines eliminated.
 */
trilines(p, nv)
float *p;
{
    float **lines = (float **)malloc(nv*2*sizeof(float *));
    int i, j;

    if (!lines) {
	fprintf(stderr, "trilines: out of memory\n");
	exit(1);
    }
    /* 3 lines per triangle */
    for(j = i = 0; i < nv/3; i++, p+=9) {
	lines[2*j+0] = &p[0];
	lines[2*j+1] = &p[3];
	j++;
	lines[2*j+0] = &p[3];
	lines[2*j+1] = &p[6];
	j++;
	lines[2*j+1] = &p[6];
	lines[2*j+0] = &p[0];
	j++;
    }
    printf("lines = %d\n", j);
    j = line_mesh(lines, j*2);
    free((char *)lines);
    return j;
}

/*
 *  take a list of triangle vertices and normal vectors, convert these
 *  to short line segments and generate a display list with redundant
 *  vectors eliminated.
 */
trinorms(v, n, nv, sc)
float *v, *n, sc;
{
    int i, obj;
    float *nn, *norms = (float *)malloc(6*nv*sizeof(float));
    int dup = 0;

    if (!norms) {
	fprintf(stderr, "trinorms: out of memory\n");
	exit(1);
    }
    if (sc <= 0.) sc = 1.0;
    for (nn = norms, i = 0; i < nv; i++) {
	*nn++ = v[0]; *nn++ = v[1]; *nn++ = v[2];
	*nn++ = v[0] + n[0]*sc;
	*nn++ = v[1] + n[1]*sc;
	*nn++ = v[2] + n[2]*sc;
	v += 3;
	n += 3;
    }
    v -= 3*nv;

    /* build hash table of unique lines */
    for(i = 0; i < nv; i++) {
	int h = hashit(&norms[6*i]);
	struct hash *hh;

	/* if line already exists - ignore */
	if (htable[h]) {
	    for(hh = htable[h]; hh; hh=hh->next) {
		if (hh->v1[0] == norms[6*i+0] &&
		    hh->v1[1] == norms[6*i+1] &&
		    hh->v1[2] == norms[6*i+2] &&
		    hh->v2[0] == norms[6*i+3+0] &&
		    hh->v2[1] == norms[6*i+3+1] &&
		    hh->v2[2] == norms[6*i+3+2]) goto next;
	    }
	}
	hh = (struct hash *)malloc(sizeof *hh);
	hh->next = htable[h]; htable[h] = hh;
	hh->v1 = norms+6*i; hh->v2 = norms+6*i+3;
	continue;
next:;
	dup++;
    }
    VERBOSE(printf ("lines = %d dups = %d\n", nv, dup));

    { int vp = 0, np = 0;
    for(i = 0; i < TABLE_SIZE; i++) {
	struct hash *p;
	if (!htable[i]) continue;
	for(p = htable[i]; p; p = p->next) {
	    vp++;
	}
	np++;
    }
    printf("vp %d np %d  %f\n", vp, np, (float)vp/np);
    }

    /* now make display list */
    makeobj(obj = genobj());
    for(i = 0; i < TABLE_SIZE; i++) {
	struct hash *p, *p1;
	for(p = htable[i]; p;) {
	    bgnline();  v3f(p->v1); v3f(p->v2); endline();
	    p1 = p->next;
	    free((char *)p);
	    p = p1;
	}
	htable[i] = 0;
    }
    closeobj();
    free((char *)norms);
    return obj;
}

static struct hash *manage;

struct hash *
newhash() {
#define N 1001
    static avail;
    static struct hash *cur;
    struct hash *h;

    if (!avail) {
	if ((cur = (struct hash *) malloc(N*sizeof h[0])) == NULL) {
	    fprintf(stderr, "line_mesh: out of memory\n");
	    exit(1);
	}
	cur->next = manage;
	manage = cur;
	avail = N-1;
	cur++;
    }
    avail -= 2; h = cur; cur += 2;
    return h;
}

freehash() {
    struct hash *h, *h1;
    for(h = manage; h;) {
	h1 = h->next; free((char *)h); h = h1;
    }
}
#endif
