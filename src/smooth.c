#define WEIGHTED
#ifdef WEIGHTED
#include <math.h>
#endif

struct hash {
    int tri;	/* triangle # */
    int off;	/* vertex in triangle */
#ifdef WEIGHTED
    float w;
#endif
    struct hash *next;
};

#define TABLE_SIZE 9967

static struct hash *htable[TABLE_SIZE];

static
hashit(v)
float *v;
{
   unsigned long *p = (unsigned long *)v;
   return (((p[0]*283+p[1])*283)+p[2]) % TABLE_SIZE;
}

static
hinsert(v, t, o)
float *v;
{
    int h = hashit(v);
    struct hash *hh = (struct hash *)malloc(sizeof *hh);
    hh->next = htable[h]; htable[h] = hh;
    hh->tri = t; hh->off = o;
}

smooth(v, n, nv)
float *v, *n;
{
    int i;
    int vp = 0;
    int np = 0;
    int csrej = 0;

    for(i = 0; i < nv/3; i++) {
	hinsert(v+i*9, i, 0);
	hinsert(v+i*9+3, i, 1);
	hinsert(v+i*9+6, i, 2);
    }

    for(i = 0; i < TABLE_SIZE; i++) {
	struct hash *p;
	if (!htable[i]) continue;
	for(p = htable[i]; p; p = p->next) {
	    vp++;
	}
	/*
	printf("%d  %d\n", i, vp);
	*/
	np++;
    }
    printf("vp %d np %d  %f\n", vp, np, (float)vp/np);
    vp = 0;

    for(i = 0; i < TABLE_SIZE; i++) {
	while(htable[i]) {
	    struct hash *p, *nequal = 0, *equal = htable[i];
	    int j = 1;
	    float nn[3], vv[3], cs;
	    vv[0] = v[equal->tri*9+equal->off*3];
	    vv[1] = v[equal->tri*9+equal->off*3+1];
	    vv[2] = v[equal->tri*9+equal->off*3+2];

	    nn[0] = n[equal->tri*9+equal->off*3];
	    nn[1] = n[equal->tri*9+equal->off*3+1];
	    nn[2] = n[equal->tri*9+equal->off*3+2];

	    /* make a list of equal vertices */
	    for(p = equal->next, equal->next = 0; p;) {
		struct hash *p1 = p->next;
		if (v[p->tri*9+p->off*3] == vv[0] &&
		    v[p->tri*9+p->off*3+1] == vv[1] &&
		    v[p->tri*9+p->off*3+2] == vv[2]) {
		    cs = n[p->tri*9+p->off*3]*nn[0] +
			 n[p->tri*9+p->off*3+1]*nn[1] +
			 n[p->tri*9+p->off*3+2]*nn[2];
		    if (cs > .707) {
			p->next = equal; equal = p;
			j++;
			p = p1;
			continue;
		    } else csrej++;
		}
		p->next = nequal; nequal = p;
		p = p1;
	    }
	    htable[i] = nequal;
	    vp += j;
	    if (j == 1) {
		free((char *)equal);
		continue;
	    }
#ifdef DEBUG
	    printf("adjacencies = %d\n", j);
	    for(p = equal; p; p = p->next) {
		printf("%d %d %f %f %f   %f %f %f\n", p->tri, p->off,
	v[p->tri*9+p->off*3], v[p->tri*9+p->off*3+1], v[p->tri*9+p->off*3+2],
	n[p->tri*9+p->off*3], n[p->tri*9+p->off*3+1], n[p->tri*9+p->off*3+2]);
	    }
#endif
#ifdef WEIGHTED
	    /* calculate weights == distances from triangle centroids */
	    for(cs = 0., p = equal; p; p = p->next) {
		nn[0] = (v[p->tri*9]+v[p->tri*9+3]+v[p->tri*9+6])/3.
			- v[p->tri*9+p->off*3];
		nn[1] = (v[p->tri*9+1]+v[p->tri*9+3+1]+v[p->tri*9+6+1])/3.
			- v[p->tri*9+p->off*3+1];
		nn[2] = (v[p->tri*9+2]+v[p->tri*9+3+2]+v[p->tri*9+6+2])/3.
			- v[p->tri*9+p->off*3+2];
#ifdef mips
		p->w = fsqrt(nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2]);
#else
		p->w = sqrt(nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2]);
#endif
		cs += p->w;
	    }
#endif
	    /* average the normal */
	    nn[0] = nn[1] = nn[2] = 0.;
#ifdef WEIGHTED
	    for (p = equal; p; p = p->next) {
		nn[0] += n[p->tri*9+p->off*3]*p->w/cs;
		nn[1] += n[p->tri*9+p->off*3+1]*p->w/cs;
		nn[2] += n[p->tri*9+p->off*3+2]*p->w/cs;
	    }
#else
	    for (p = equal; p; p = p->next) {
		nn[0] += n[p->tri*9+p->off*3];
		nn[1] += n[p->tri*9+p->off*3+1];
		nn[2] += n[p->tri*9+p->off*3+2];
	    }
	    nn[0] /= j; nn[1] /= j; nn[2] /= j;
#endif
	    /* put it back */
	    for (p = equal; p;) {
		struct hash *p1 = p->next;
		n[p->tri*9+p->off*3] = nn[0];
		n[p->tri*9+p->off*3+1] = nn[1];
		n[p->tri*9+p->off*3+2] = nn[2];
#ifdef DEBUG
		printf("%d %d %f %f %f   %f %f %f\n", p->tri, p->off,
	v[p->tri*9+p->off*3], v[p->tri*9+p->off*3+1], v[p->tri*9+p->off*3+2],
	n[p->tri*9+p->off*3], n[p->tri*9+p->off*3+1], n[p->tri*9+p->off*3+2]);
#endif
		free((char *)p);
		p = p1;
	    }
	}
    }
    printf("vertices processed = %d, cosine rejects = %d\n", vp, csrej);
}
