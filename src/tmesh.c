/* 
 *	tmesh -
 *		build a display list of triangles using tmesh primitives
 *		from an array of triangles
 *
 *			Kurt Akeley and Paul Haeberli - 1990
 *
 *	Algorithm Description
 *
 *	    Triangles are sorted into 4 groups, based on whether they are
 *	    adjacent to 0, 1, 2, or 3 other triangles.  Triangles in group
 *	    0 are output immediately as independent triangles.  Then a
 *	    seed triangle with minimum adjacencies is chosen (i.e. a triangle
 *	    is selected from group 1 if it is not empty, otherwise from
 *	    group 2 if it is not empty, otherwise from group 3).  A mesh
 *	    is extended from this triangle, with triangles again selected
 *	    based on minimum adjacencies.  Adjacency counts for all triangles
 *	    are updated as each triangle is used.
 *
 *	    Because triangles with minimum adjacencies are selected, the
 *	    mesh consumes triangle patches along their edges, rather than
 *	    through their interiors.  This results in long meshes.
 */

#include <stdio.h>
#ifdef SGI
#include <gl/gl.h>
/*
*/

/* eliminate degenerate and redundant triangles */
/*#define ELIMINATE_DEGENERATE*/
#define ELIMINATE_REDUNDANT

/* test adjacency pointers for consistency */
/*#define CHECK_CONSISTENCY*/
/* more self checking */
/*#define CAUTIOUS*/

/*#define VERBOSE(x)	x*/
#define VERBOSE(x)

tmesh(p, n, nv, nflag)
float *p, *n;
{
    long o = genobj();
    makeobj(o);
    do_tmesh(p, n, nv, nflag);
    closeobj();
    return o;
}

#ifndef FALSE
#define FALSE 0
#define TRUE 1
#endif

#define VERTHASHSIZE 9991
#define TRIHASHSIZE 9991
#define EDGEHASHSIZE 9991

typedef struct Vert {
    struct Vert *next,*prev;
    struct Vert *hnext;
    float x,y,z;
    float nx,ny,nz;
#ifdef TEXTURE
    float tx,ty,tz;
#endif
    int index;
} Vert;

typedef struct Vertlist {
    Vert *first,*last;
} Vertlist;

typedef struct Edge {
    struct Edge *hnext;
    struct Tri *tri;
    int index0, index1;
} Edge;

typedef struct Tri {
    struct Tri *next,*prev;
    struct Tri *hnext;
    struct Tri *(adj[3]);
    Vert *vert[3];
#ifdef ndef
    int drawvert[3];
#endif
    int adjcount;
    int used;
} Tri;

typedef struct Trilist {
    Tri *first,*last;
} Trilist;

Tri *minadj();

static int npolys;
static int connectcount = 0;
static int independentcount = 0;
static Trilist *adjtrilist[4];
#ifdef CAUTIOUS
static Trilist *donetrilist;
#endif
static Vert **verthashlist;
static Tri **trihashlist;
static Edge **edgehashlist;
static Edge *edgearray;
static Edge *freeedges;

/* 
 *	Vertex hashing
 */

hashvertbegin()
{
    int i;

    verthashlist = (Vert**)malloc(VERTHASHSIZE*sizeof(Vert*));
    for(i=0; i<VERTHASHSIZE; i++) 
	verthashlist[i] = 0;
}

hashvert(vert)
Vert *vert;
{
    long *buf = (long *)&vert->x;
    long hash = 0;
#ifdef TEXTURE
    int n = 9;
#else
    int n = 6;
#endif

#ifdef ndef
    while(n--)
	hash = (273*hash)+(*buf++);
#else
    hash = 273*(273*(273*(273*((273*buf[0]) + buf[1]) + buf[2]) + buf[3]) + buf[4]) + buf[5];
#endif
    hash &= 0x7fffffff;

    return hash%VERTHASHSIZE;
}

#ifdef ndef
vequal(p,q)
Vert *p,*q;
{
#define close(a,b) (a == b)

    return (close(p->x,q->x) && close(p->y,q->y) && close(p->z,q->z) &&
#ifdef TEXTURE
	close(p->nx,q->nx) && close(p->ny,q->ny) && close(p->nz,q->nz) &&
	close(p->tx,q->tx) && close(p->ty,q->ty) && close(p->tz,q->tz));
#else
	close(p->nx,q->nx) && close(p->ny,q->ny) && close(p->nz,q->nz));
#endif
}
#else
#define vequal(p,q) (p->x==q->x && p->y==q->y && p->z==q->z && \
	             p->nx==q->nx && p->ny==q->ny && p->nz==q->nz)
#endif

Vert *
hashvertadd(vert)
Vert *vert;
{
    int pos;
    Vert *vptr;

    pos = hashvert(vert);

    /* check if already present */
    for(vptr = verthashlist[pos]; vptr; vptr = vptr->hnext)
	if(vequal(vert,vptr)) return vptr;
    vert->hnext = verthashlist[pos];
    verthashlist[pos] = vert;
    return 0;
}

hashvertend()
{
    free(verthashlist);
}

/* 
 *	Triangle hashing
 */

hashtribegin()
{
    int i;

    trihashlist = (Tri**)malloc(TRIHASHSIZE*sizeof(Tri*));
    for(i=0; i<TRIHASHSIZE; i++) 
	trihashlist[i] = 0;
}

hashtri(tri)
Tri *tri;
{
    long val;

    val = (long)tri->vert[0];
    val ^= (long)tri->vert[1];
    val ^= (long)tri->vert[2];
    return val%VERTHASHSIZE;
}

hashtriadd(tri)
Tri *tri;
{
    int pos;
    Tri *tptr;

    for(tptr = trihashlist[pos=hashtri(tri)]; tptr; tptr = tptr->hnext)
	if(triequal(tri,tptr)) return 1;
    tri->hnext = trihashlist[pos];
    trihashlist[pos] = tri;
    return 0;
}

hashtriend()
{
    free(trihashlist);
}

/* 
 *	Edge hashing
 */

hashedgebegin()
{
    int i;

    edgehashlist = (Edge**)malloc(EDGEHASHSIZE*sizeof(Edge*));
    edgearray = freeedges = (Edge*)malloc(3*npolys*sizeof(Edge));
    for(i=0; i<EDGEHASHSIZE; i++) 
	edgehashlist[i] = 0;
}

hashedge(index0,index1)
int index0, index1;
{
    long val;

    val = index0*index1;
    val = val&0x7fffffff;
    return val%EDGEHASHSIZE;
}

Edge *
hashedgefind(v0,v1)
Vert *v0, *v1;
{
#ifdef ndef
    int pos;
    int index0, index1;
    Edge *tptr;

    index0 = v0->index;
    index1 = v1->index;
    pos = hashedge(index0,index1);
    tptr = edgehashlist[pos];
    return tptr;
#else
    return edgehashlist[hashedge(v0->index,v1->index)];
#endif
}

hashedgeadd(tri,v0,v1)
Tri *tri;
Vert *v0, *v1;
{
    int pos;
    int index0, index1;
    Edge *edge;

    index0 = v0->index;
    index1 = v1->index;
    pos = hashedge(index0,index1);
    edge = freeedges++;
    edge->index0 = index0;
    edge->index1 = index1;
    edge->tri = tri;
    edge->hnext = edgehashlist[pos];
    edgehashlist[pos] = edge;
}

hashedgeend()
{
    free(edgehashlist);
    free(edgearray);
}

Vertlist *
makevertlist() 
{
    /* allocate space for and initialize a vert list */
    Vertlist *tmp;

    if ((tmp=(Vertlist*)(malloc(sizeof(Vertlist))))==0) {
	fprintf(stderr,"makevertlist: out of memory.  abort.\n");
	exit(1);
    }
    tmp->first = tmp->last = 0;
    return tmp;
}

Trilist *
maketrilist() 
{
    /* allocate space for and initialize a tri list */
    Trilist *tmp;

    if ((tmp=(Trilist*)(malloc(sizeof(Trilist))))==0) {
	fprintf(stderr,"maketrilist: out of memory.  abort.\n");
	exit(1);
    }
    tmp->first = tmp->last = 0;
    return tmp;
}

Vert *
makevert() 
{
    /* allocate space for and initialize a vert */
    Vert *tmp;

#ifdef MALLOC
    if ((tmp=(Vert*)(malloc(sizeof(Vert))))==0) {
	fprintf(stderr,"makevert: out of memory.  abort.\n");
	exit(1);
    }
#else
    Vert *newvert();
    tmp = newvert();
#endif
    tmp->prev = tmp->next = 0;
    tmp->index = 0;
    return tmp;
}

Tri *
maketri() 
{
    /* allocate space for and initialize a tri */
    Tri *tmp;

    register i;
#ifdef MALLOC
    if ((tmp=(Tri*)(malloc(sizeof(Tri))))==0) {
	fprintf(stderr,"maketri: out of memory.  abort.\n");
	exit(1);
    }
#else
    Tri *newtri();
    tmp = newtri();
#endif
    tmp->prev = tmp->next = 0;
    for (i=0; i<3; i++) {
	tmp->adj[i] = 0;
	tmp->vert[i] = 0;
    }
#ifdef ndef
    tmp->drawvert[0] = 0;
    tmp->drawvert[1] = 1;
    tmp->drawvert[2] = 2;
#endif
    tmp->adjcount = 0;
    tmp->used = FALSE;
    return tmp;
}

void
insertvert(list,item)
Vertlist *list;
Vert *item;
{
    /* insert the new item at the top of the list */
    Vert *tmp;

    if (list->first) {
	tmp = list->first;
	list->first = item;
	item->next = tmp;
	item->prev = 0;
	tmp->prev = item;
    } else {
	list->first = list->last = item;
	item->prev = item->next = 0;
    }
}

void
inserttri(list,item)
Trilist *list;
Tri *item;
{
    /* insert the new item at the top of the list */
    Tri *tmp;

    if (list->first) {
	tmp = list->first;
	list->first = item;
	item->next = tmp;
	item->prev = 0;
	tmp->prev = item;
    } else {
	list->first = list->last = item;
	item->prev = item->next = 0;
    }
}

void
deletevert(list,item)
Vertlist *list;
Vert *item;
{
    /* delete the item from the list */
    if (list->first == item) {
	if (list->last == item) {
	    /* this is the only item in the list */
	    list->first = list->last = 0;
	} else {
	    /* this is the first item in the list */
	    list->first = item->next;
	    list->first->prev = 0;
	}
    } else if (list->last == item) {
	/* this is the last item in the list */
	list->last = item->prev;
	list->last->next = 0;
    } else {
	item->prev->next = item->next;
	item->next->prev = item->prev;
    }
#ifdef CAUTIOUS
    item->prev = item->next = 0;
#endif
}

void deletetri(list,item)
Trilist *list;
Tri *item;
{
    /* delete the item from the list */
    if (list->first == item) {
	if (list->last == item) {
	    /* this is the only item in the list */
	    list->first = list->last = 0;
	} else {
	    /* this is the first item in the list */
	    list->first = item->next;
	    list->first->prev = 0;
	}
    } else if (list->last == item) {
	/* this is the last item in the list */
	list->last = item->prev;
	list->last->next = 0;
    } else {
	item->prev->next = item->next;
	item->next->prev = item->prev;
    }
#ifdef CAUTIOUS
    item->prev = item->next = 0;
#endif
}

triequal(tri1,tri2) 
Tri *tri1, *tri2;
{
    int i, j, k;

    k = 0;
#ifdef ndef
    for (i=0; i<3; i++) {
	for (j=0; j<3; j++) {
	    k += (tri1->vert[i] == tri2->vert[j]);
	}
    }
#else
    k += (tri1->vert[0] == tri2->vert[0]);
    k += (tri1->vert[0] == tri2->vert[1]);
    k += (tri1->vert[0] == tri2->vert[2]);
    k += (tri1->vert[1] == tri2->vert[0]);
    k += (tri1->vert[1] == tri2->vert[1]);
    k += (tri1->vert[1] == tri2->vert[2]);
    k += (tri1->vert[2] == tri2->vert[0]);
    k += (tri1->vert[2] == tri2->vert[1]);
    k += (tri1->vert[2] == tri2->vert[2]);
#endif
    return k == 3;
}

do_tmesh(p, n, nv, nf)
float *p, *n;
{
    register Vert *vert;
    register Vert *tmpvert;
    Vert *vert0,*vert1;
    register Tri *tri;
    register Tri *tmptri;
    register Tri *nexttri;
    register Edge *edge;
    int i, j, k;
    int vertcount;
    int degeneratecount;
    int equivalentcount;
    Vertlist *vertlist;
    Trilist *trilist;
    Trilist *newtrilist;
    int count;
    int adjcount[4];
    int adjnumber;
    int poly;

    /*** initialize lists ***/
    vertlist = makevertlist();
    trilist = maketrilist();
    newtrilist = maketrilist();
#ifdef CAUTIOUS
    donetrilist = maketrilist();
#endif
    for (i=0; i<4; i++)
	adjtrilist[i] = maketrilist();

    tmpvert = makevert();
    vertcount = 0;

    npolys = nv/3;
    hashvertbegin();
    VERBOSE(tpercentdone(0.0));
    for(poly=0; poly<npolys; poly++) {
	VERBOSE(tpercentdone(100.0*poly/(npolys-1)));
	tri = maketri();
	inserttri(trilist,tri);
	for (i=0; i<3; i++) {
	    tri->vert[i] = 0;
	    tmpvert->x = *p++;
	    tmpvert->y = *p++;
	    tmpvert->z = *p++;
	    tmpvert->nx = *n++;
	    tmpvert->ny = *n++;
	    tmpvert->nz = *n++;
#ifdef TEXTURE
	    tmpvert->tx = 0;
	    tmpvert->ty = 0;
	    tmpvert->tz = 0;
#endif
	    vert = hashvertadd(tmpvert);
	    if(vert)
		tri->vert[i] = vert;
	    else  {
		/* add a new vertex to the list */
		tmpvert->index = vertcount;
		vertcount++;
		tri->vert[i] = tmpvert;
		insertvert(vertlist,tmpvert);
		tmpvert = makevert();
	    }
	}
    }
#ifdef MALLOC
    free((char *)tmpvert);
#endif
    VERBOSE(tpercentdone(100.0));
    hashvertend();
    VERBOSE(printf("%d triangles read\n",npolys));
    VERBOSE(printf("%d vertices created\n",vertcount));

#ifdef ELIMINATE_DEGENERATE
    /*** eliminate degenerate triangles ***/
    VERBOSE(printf("eliminating degenerate triangles\n"));
    degeneratecount = 0;
    for (tri=trilist->first; tri;) {
	if ((tri->vert[0] == tri->vert[1]) ||
	    (tri->vert[0] == tri->vert[2]) ||
	    (tri->vert[1] == tri->vert[2])) {
	    degeneratecount += 1;
	    tmptri = tri->next;
	    deletetri(trilist,tri);
	    tri = tmptri;
	} else 
	    tri = tri->next;
    }
    VERBOSE(printf("%d degenerate triangles eliminated\n",degeneratecount));
#endif

#ifdef ELIMINATE_REDUNDANT
    /*** eliminate equivalent triangles ***/
    VERBOSE(printf("eliminating equivalent triangles\n"));
    count = 0;
    equivalentcount = 0;

    hashtribegin();
    VERBOSE(tpercentdone(0.0));
    for (tri=trilist->first; tri;) {
	VERBOSE(tpercentdone(100.0*count/(npolys-1)));
	count += 1;
	if(hashtriadd(tri)) {
	    equivalentcount += 1;
	    nexttri = tri->next;
	    deletetri(trilist,tri);
	    tri = nexttri;
	} else {
	    tri = tri->next;
	}
    }
    VERBOSE(tpercentdone(100.0));
    hashtriend();
    VERBOSE(printf("%d equivalent triangles eliminated\n",equivalentcount));
#endif

    /*** compute triangle adjacencies ***/
    VERBOSE(printf("computing adjacent triangles\n"));
    hashedgebegin();
    VERBOSE(printf("adding to hash table . . "));
    for (tri=trilist->first; tri; tri=tri->next) {
#ifdef ndef
	for(i=0; i<3; i++) {
	    vert0 = tri->vert[(i+0)%3];
	    vert1 = tri->vert[(i+1)%3];
	    hashedgeadd(tri,vert0,vert1);
	}
#else
	hashedgeadd(tri,tri->vert[0],tri->vert[1]);
	hashedgeadd(tri,tri->vert[1],tri->vert[2]);
	hashedgeadd(tri,tri->vert[2],tri->vert[0]);
#endif
    }
    VERBOSE(printf("done\n"));
    count = 0;
    VERBOSE(tpercentdone(0.0));
    for (tri=trilist->first; tri; tri=tri->next) {
	VERBOSE(tpercentdone(100.0*count/(npolys-1)));
	count += 1;
	for (i=0; i<3; i++) {
	    static suc[3] = {1, 2, 0};
	    if (tri->adj[i]) continue;
	    vert0 = tri->vert[i];
#ifdef ndef
	    vert1 = tri->vert[(i+1)%3];
#else
	    vert1 = tri->vert[suc[i]];
#endif
	    for (edge=hashedgefind(vert0,vert1); edge; edge=edge->hnext) {
		nexttri = edge->tri;
		if(nexttri == tri) continue;
		for (j=0; j<3; j++) {
		    if (vert0 == nexttri->vert[j]) {
			for (k=0; k<3; k++) {
			    if (k==j) continue;
			    if (vert1 == nexttri->vert[k]) {
#ifdef ndef
				switch (j+k) {
				    case 1:
					adjnumber = 0;
					break;
				    case 2:
					adjnumber = 2;
					break;
				    case 3:
					adjnumber = 1;
					break;
				    default:
					fprintf(stderr,
					    "ERROR: bad adjnumber\n");
					break;
				}
#else
				{
				static t[4] = {0, 0, 2, 1};
				adjnumber = t[j+k];
				}
#endif
				if (tri->adj[i]||nexttri->adj[adjnumber]) {
				} else {	
				    tri->adj[i] = nexttri;
				    nexttri->adj[adjnumber] = tri;
				}
			    }
			}
		    }
		}
	    }
	}
    }
    VERBOSE(tpercentdone(100.0));
    hashedgeend();
    VERBOSE(printf(" done\n"));

#ifdef CHECK_CONSISTENCY
    /*** test adjacency pointers for consistency ***/
    for (tri=trilist->first; tri; tri=tri->next) {
	for (i=0; i<3; i++) {
	    if (nexttri = tri->adj[i]) {
		for (j=0,k=0; j<3; j++) {
		    if (tri == nexttri->adj[j])
			k += 1;
		}
		if (k != 1) {
		    fprintf(stderr," ERROR: %x to %x k = %d\n",tri,nexttri,k);
		}
	    }
	}
    }
#endif

    /*** compute adjacency statistics ***/
    for (i=0; i<4; i++)
	adjcount[i] = 0;
    for (tri=trilist->first; tri;) {
#ifdef ndef
	for (i=0,count=0; i<3; i++) {
	    if (tri->adj[i])
		count += 1;
	}
#else
	count = (tri->adj[0] != 0) + (tri->adj[1] != 0) + (tri->adj[2] != 0);
#endif
	tri->adjcount = count;
	adjcount[count] += 1;
	nexttri = tri->next;
	deletetri(trilist,tri);
	inserttri(adjtrilist[count],tri);
	tri = nexttri;
    }
    printf("adjacencies: 0:%d, 1:%d, 2:%d, 3:%d\n",
	adjcount[0],adjcount[1],adjcount[2],adjcount[3]);

    /*** search for connected triangles and output ***/
    while (1) {
	/*** output singular triangles, if any ***/
	while (tri = adjtrilist[0]->first) {
	    deletetri(adjtrilist[0],tri);
	    inserttri(newtrilist,tri);
	    tri->used = TRUE;
	    outmesh(newtrilist);
	}
	/*** choose a seed triangle with the minimum number of adjacencies ***/
	if (tri = adjtrilist[1]->first)
	    deletetri(adjtrilist[1],tri);
	else if (tri = adjtrilist[2]->first)
	    deletetri(adjtrilist[2],tri);
	else if (tri = adjtrilist[3]->first)
	    deletetri(adjtrilist[3],tri);
	else
	    break;
	inserttri(newtrilist,tri);
	removeadjacencies(tri);

	/*** extend in one direction using triangles with min adjacencies ***/
	while (tri = minadj(tri)) {
	    deletetri(adjtrilist[tri->adjcount],tri);
	    inserttri(newtrilist,tri);
	    removeadjacencies(tri);
	}

	/*** if seed has two or more adjacencies, extend in other direction **/
	tri = newtrilist->first;
	nexttri = 0;
	for (i=0; i<3; i++) {
	    if (tri->adj[i] &&
		(tri->adj[i] != tri->next) &&
		(!(tri->adj[i]->used))) {
		nexttri = tri->adj[i];
		break;
	    }
	}
	tri = nexttri;
	while (tri) {
	    deletetri(adjtrilist[tri->adjcount],tri);
	    inserttri(newtrilist,tri);
	    removeadjacencies(tri);
	    tri = minadj(tri);
	}

	/*** output the resulting mesh ***/
	outmesh(newtrilist);
    }
    if (connectcount)
	printf("%d triangle mesh%s output\n",
	    connectcount,connectcount==1?"":"es");
    if (independentcount)
	printf("%d independent triangle%s output\n",
	    independentcount,independentcount==1?"":"s");
}

int ismember(vert,tri)
Vert *vert;
Tri *tri;
{
    /*** return TRUE if vert is one of the vertexes in tri, otherwise FALSE ***/
    register int i;

#ifdef ndef
    for (i=0; i<3; i++)
	if (vert == tri->vert[i])
	    return TRUE;
    return FALSE;
#else
    return (vert==tri->vert[0]) | (vert==tri->vert[1]) | (vert==tri->vert[2]);
#endif
}

int notcommon(tri,tri2)
Tri *tri,*tri2;
{
    /*** returns the index of the vertex in tri that is not in tri2 ***/
    int i;

    for (i=0; i<3; i++)
	if (!ismember(tri->vert[i],tri2))
	    return i;
    return -1;
}

#define INITMESH	replaceB = FALSE
#define SWAPMESH	replaceB = !replaceB
#define REPLACEVERT(v)	{vert[replaceB] = (v); SWAPMESH;}

outmesh(trilist)
Trilist *trilist;
{
    Tri *tri;
    int i;
    Vert *vert[2];
    Vert *nextvert;
    int replaceB;

    /*** output trilist - transfer to donelist ***/
    tri = trilist->first;
    if (tri == trilist->last) {
	/*** there is only one triangle in the mesh - use polygon command ***/
	independentcount += 1;
	bgnpolygon();
	n3f(&tri->vert[0]->nx); v3f(&tri->vert[0]->x);
	n3f(&tri->vert[1]->nx); v3f(&tri->vert[1]->x);
	n3f(&tri->vert[2]->nx); v3f(&tri->vert[2]->x);
	endpolygon();
	deletetri(trilist,tri);
#ifdef CAUTIOUS
	inserttri(donetrilist,tri);
#endif
    } else {
	/*** a real mesh output is required ***/
	connectcount += 1;
	/*** start output with vertex that is not in the second triangle ***/
	i = notcommon(tri,tri->next);
	INITMESH;
	bgntmesh();
	n3f(&tri->vert[i]->nx); v3f(&tri->vert[i]->x);
	REPLACEVERT(tri->vert[i]);
	n3f(&tri->vert[(i+1)%3]->nx); v3f(&tri->vert[(i+1)%3]->x);
	REPLACEVERT(tri->vert[(i+1)%3]);
	n3f(&tri->vert[(i+2)%3]->nx); v3f(&tri->vert[(i+2)%3]->x);
	REPLACEVERT(tri->vert[(i+2)%3]);
	/*** compute vertex of second triangle that is not in the first ***/
	i = notcommon(tri->next,tri);
	nextvert = (tri->next)->vert[i];
	/*** transfer triangle to done list ***/
	deletetri(trilist,tri);
#ifdef CAUTIOUS
	inserttri(donetrilist,tri);
#endif
	tri = trilist->first;
	while (tri->next) {
	    /*** check for errors ***/
#ifdef CAUTIOUS
	    if ((!ismember(vert[0],tri)) || (!ismember(vert[1],tri)) ||
		(!ismember(nextvert,tri))) {
		fprintf(stderr,"ERROR in mesh generation\n");
	    }
	    if ((vert[0] == vert[1]) || (vert[0] == nextvert)) {
		fprintf(stderr,"ERROR in mesh generation\n");
	    }
#endif
	    /*** decide whether to swap or not ***/
	    if (ismember(vert[replaceB],tri->next)) {
		swaptmesh();
		SWAPMESH;
	    }
	    /*** output the next vertex ***/
	    n3f(&nextvert->nx); v3f(&nextvert->x);
	    REPLACEVERT(nextvert);
	    /*** determine the next output vertex ***/
	    i = notcommon(tri->next,tri);
	    nextvert = (tri->next)->vert[i];   
	    /*** transfer tri to the done list ***/
	    deletetri(trilist,tri);
#ifdef CAUTIOUS
	    inserttri(donetrilist,tri);
#endif
	    tri = trilist->first;
	}
	/*** output the last vertex ***/
	n3f(&nextvert->nx); v3f(&nextvert->x);
	REPLACEVERT(nextvert);
	deletetri(trilist,tri);
#ifdef CAUTIOUS
	inserttri(donetrilist,tri);
#endif
	endtmesh();
    }
}

removeadjacencies(tri)
Tri *tri;
{
    register i,j;
    Tri *adjtri;

    tri->used = TRUE;
    for (i=0; i<3; i++) {
	adjtri = tri->adj[i];
	if (adjtri) {
	    deletetri(adjtrilist[adjtri->adjcount],adjtri);
	    adjtri->adjcount -= 1;
	    for (j=0; j<3; j++) {
		if (tri == adjtri->adj[j]) {
		    adjtri->adj[j] = 0;
		    break;
		}
	    }
	    inserttri(adjtrilist[adjtri->adjcount],adjtri);
	}
    }
}

Tri *minadj(tri)
Tri *tri;
{
    int min0,min1;

    switch (tri->adjcount) {
	case 0:
	    return 0;
	case 1:
	    if (tri->adj[0])
		return tri->adj[0];
	    else if (tri->adj[1])
		return tri->adj[1];
	    else
		return tri->adj[2];
	case 2:
	    if (tri->adj[0]) {
		min0 = 0;
		if (tri->adj[1])
		    min1 = 1;
		else
		    min1 = 2;
	    } else {
		min0 = 1;
		min1 = 2;
	    }
	    if ((tri->adj[min0])->adjcount <= (tri->adj[min1])->adjcount)
		return tri->adj[min0];
	    else
		return tri->adj[min1];
	case 3:
	    if ((tri->adj[0])->adjcount <= (tri->adj[1])->adjcount)
		min0 = 0;
	    else
		min0 = 1;
	    min1 = 2;
	    if ((tri->adj[min0])->adjcount <= (tri->adj[min1])->adjcount)
		return tri->adj[min0];
	    else
		return tri->adj[min1];
    }
    /*NOTREACHED*/
}


/*
 *	tpercent - 
 *		Make a row of dots that show percent done . . .
 *
 *				Paul Haeberli - 1989
 *
 */
#include "stdio.h"

static int started = 0;
static int pos;
static FILE *outf;

#define NDOTS	66

tpercentfile(f)
FILE *f;
{
    outf = f;
}

tpercentdone(p)
float p;
{
    int newpos;
    FILE *outfile;

    if(!outf)
	outfile = stderr;
    else
	outfile = outf;
    p = p/100.0;
    if(!started && p <= 0.01) {
	fprintf(outfile,"working: [");
	fflush(outfile);
	started = 1;
	pos = 0;
	return;
    }
    if(started) {
	if(p<0.999) 
	    newpos = NDOTS*p;
	else
	    newpos = NDOTS;
        if(newpos>pos) {
	    while(pos<newpos) {
	        fprintf(outfile,".");
		pos++;
	    }
	    fflush(outfile);
        }
	if(p>0.999) {
	    fprintf(outfile,"]\n");
	    fflush(outfile);
	    started = 0;
	}
    }
}

#ifdef MAIN
#include <sys/param.h>
#include <sys/types.h>
#include <sys/times.h>

main() {
   int verts, norms;
   float *v, *n;
   struct tms t1, t2;

   scanf("%*s %d %d\n", &verts, &norms);
   printf("verts = %d norms = %d\n", verts, norms);
   v = (float *)malloc(sizeof(v[0])*verts*3);
   fread(v, sizeof(v[0]), verts*3, stdin);
   n = (float *)malloc(sizeof(n[0])*verts*3);
   fread(n, sizeof(n[0]), verts*3, stdin);
   times(&t1);
   tmesh(v, n, verts, 0);
   times(&t2);
   printf("time: %f\n", (float)(t2.tms_utime-t1.tms_utime)/HZ);
   free((char *)n);
   free((char *)v);
}

bgnpolygon() { }
endpolygon() { }
bgntmesh() { }
endtmesh() { }
swaptmesh() { }
makeobj() { }
closeobj() { }
genobj() { }
v3f() { }
n3f() { }
#endif

static Vert *vmanage;

Vert *
newvert() {
#define N 1001
    static avail;
    static Vert *cur;
    Vert *v;

    if (!avail) {
	cur = (Vert *) malloc(N*sizeof v[0]);
	cur->next = vmanage;
	vmanage = cur;
	avail = N-1;
	cur++;
    }
    avail--; v = cur++;
    return v;
}

freeverts() {
    Vert *v, *v1;
    for(v = vmanage; v;) {
	v1 = v->next; free((char *)v); v = v1;
    }
}

static Tri *tmanage;

Tri *
newtri() {
#define N 1001
    static avail;
    static Tri *cur;
    Tri *t;

    if (!avail) {
	cur = (Tri *) malloc(N*sizeof t[0]);
	cur->next = tmanage;
	tmanage = cur;
	avail = N-1;
	cur++;
    }
    avail--; t = cur++;
    return t;
}

freetris() {
    Tri *t, *t1;
    for(t = tmanage; t;) {
	t1 = t->next; free((char *)t); t = t1;
    }
}
#endif
