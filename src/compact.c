#include <stdio.h>

#define VERBOSE(x)	x

struct hash {
    float *v1, *n1;
    float **v;
    struct hash *next;
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

/*
 *  take a list of triangle vertices and normal vectors, produce
 *  a connectivity array with redundant vertices eliminated.
 */
float ***
tricompactn(v, n, nv)
float *v, *n;
{
    int i;
    float ***conn = (float ***)malloc((nv+1)*sizeof(float **));
    float **vert = (float **)malloc(nv*sizeof(float *));
    int nvert = 0;
    int dup = 0;

    if (!conn || !vert) {
	fprintf(stderr, "tricompact: out of memory\n");
	exit(1);
    }
    for(i = 0; i < nv; i++) {
	int h = hashit(v+3*i);
	struct hash *hh, *newhash();

	/* if vertex already exists - ignore */
	for(hh = htable[h]; hh; hh=hh->next) {
	    if (hh->v1[0] == v[3*i+0] &&
		hh->v1[1] == v[3*i+1] &&
		hh->v1[2] == v[3*i+2] &&
		hh->n1[0] == n[3*i+0] &&
		hh->n1[1] == n[3*i+1] &&
		hh->n1[2] == n[3*i+2]) {
		    goto vdup;
	    }
	}
	vert[nvert] = v+3*i;
	hh = newhash();
	hh->next = htable[h]; htable[h] = hh;
	hh->v1 = v+3*i; hh->n1 = n+3*i; hh->v = vert+nvert;
	conn[i] = vert+nvert;
	nvert++;
	continue;
vdup:
	conn[i] = hh->v;
	dup++;
    }
    conn[i] = vert+nvert;
    VERBOSE(printf ("verts = %d dups = %d\n", nv, dup));

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
    freehash();
    return conn;
}

float ***
tricompact(v, nv)
float *v;
{
    int i;
    float ***conn = (float ***)malloc((nv+1)*sizeof(float **));
    float **vert = (float **)malloc(nv*sizeof(float *));
    int nvert = 0;
    int dup = 0;

    if (!conn || !vert) {
	fprintf(stderr, "tricompact: out of memory\n");
	exit(1);
    }
    for(i = 0; i < nv; i++) {
	int h = hashit(v+3*i);
	struct hash *hh, *newhash();

	/* if vertex already exists - ignore */
	for(hh = htable[h]; hh; hh=hh->next) {
	    if (hh->v1[0] == v[3*i+0] &&
		hh->v1[1] == v[3*i+1] &&
		hh->v1[2] == v[3*i+2]) {
		    goto vdup;
	    }
	}
	vert[nvert] = v+3*i;
	hh = newhash();
	hh->next = htable[h]; htable[h] = hh;
	hh->v1 = v+3*i; hh->v = vert+nvert;
	conn[i] = vert+nvert;
	nvert++;
	continue;
vdup:
	conn[i] = hh->v;
	dup++;
    }
    conn[i] = vert+nvert;
    VERBOSE(printf ("verts = %d dups = %d\n", nv, dup));

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
    freehash();
    return conn;
}

static struct hash *manage;

static
struct hash *
newhash() {
#define N 1001
    static avail;
    static struct hash *cur;
    struct hash *h;

    if (!avail) {
	cur = (struct hash *) malloc(N*sizeof h[0]);
	cur->next = manage;
	manage = cur;
	avail = N-1;
	cur++;
    }
    avail--; h = cur++;
    return h;
}

static
freehash() {
    struct hash *h, *h1;
    for(h = manage; h;) {
	h1 = h->next; free((char *)h); h = h1;
    }
}
