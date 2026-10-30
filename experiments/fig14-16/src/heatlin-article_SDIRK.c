#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "braid.h"
static const double alpha = 0.1;
typedef struct _braid_App_struct {
    MPI_Comm comm;
    int rank, nprocs, nspace, spatial, mgrit_iterations;
    double tstop, eps, xLeft, xRight;
    double *one, *half, *two, *work, *w1, *w2, *w3, *w4, *w5, *predictor;
}
my_App;
typedef struct _braid_Vector_struct {
    int size;
    double *values;
}
my_Vector;
static void die(const char *m) {
    fprintf(stderr, "%s\n", m);
    MPI_Abort(MPI_COMM_WORLD, 1);
}
static void *xm(size_t n) {
    void*p=malloc(n);
    if(!p)die("Allocation failure.");
    return p;
}
static void alloc_work(my_App*a) {
    size_t n=(size_t)a->nspace*sizeof(double);
    a->one=xm(n);
    a->half=xm(n);
    a->two=xm(n);
    a->work=xm(n);
    a->w1=xm(n);
    a->w2=xm(n);
    a->w3=xm(n);
    a->w4=xm(n);
    a->w5=xm(n);
    a->predictor=xm(n);
}
static void free_work(my_App*a) {
    free(a->one);
    free(a->half);
    free(a->two);
    free(a->work);
    free(a->w1);
    free(a->w2);
    free(a->w3);
    free(a->w4);
    free(a->w5);
    free(a->predictor);
}
static void heat_rhs(my_App*a, const double*y, double*f) {
    int n=a->nspace;
    double dx=1.0/(n-1), scale=alpha/(dx*dx);
    f[0]=0.0;
    for(int i=1; i<n-1; i++)f[i]=scale*(y[i-1]-2.0*y[i]+y[i+1]);
    f[n-1]=0.0;
}
static void sdirk_solve(my_App*a, const double*rhs, double h, double gamma, double*z) {
    int n=a->nspace;
    double dx=1.0/(n-1), r=alpha*gamma*h/(dx*dx), lo=-r, diag=1.0+2.0*r, *cp=a->w3;
    z[0]=a->xLeft;
    cp[0]=0.0;
    for(int i=1; i<n-1; i++) {
        double d=diag-lo*cp[i-1];
        cp[i]=lo/d;
        z[i]=(rhs[i]-lo*z[i-1])/d;
    }
    z[n-1]=a->xRight;
    for(int i=n-2; i>=1; i--)z[i]-=cp[i]*z[i+1];
}
static double sdirk_step(my_App*a, const double*y, double h, double*z) {
    const double gamma=(2.0-sqrt(2.0))/2.0;
    double e=0.0;
    sdirk_solve(a,y,h,gamma,a->one);
    heat_rhs(a,a->one,a->w1);
    for(int i=0; i<a->nspace; i++)a->predictor[i]=y[i]+(1.0-gamma)*h*a->w1[i];
    a->predictor[0]=a->xLeft;
    a->predictor[a->nspace-1]=a->xRight;
    sdirk_solve(a,a->predictor,h,gamma,a->half);
    heat_rhs(a,a->half,a->w2);
    for(int i=0; i<a->nspace; i++) {
        double d=h*gamma*(a->w1[i]-a->w2[i]);
        z[i]=y[i]+(1.0-gamma)*h*a->w1[i]+gamma*h*a->w2[i];
        e+=d*d;
    }
    z[0]=a->xLeft;
    z[a->nspace-1]=a->xRight;
    return sqrt(e/(double)a->nspace);
}
static int my_Step(braid_App a, braid_Vector us, braid_Vector fs, braid_Vector u, braid_StepStatus s) {
    double ta, tb;
    int it, l, nr;
    braid_StepStatusGetTstartTstop(s, &ta, &tb);
    braid_StepStatusGetIter(s, &it);
    braid_StepStatusGetLevel(s, &l);
    braid_StepStatusGetNRefine(s, &nr);
    double e=sdirk_step(a, u->values, tb-ta, a->work);
    memcpy(u->values, a->work, (size_t)u->size*sizeof(double));
    if(l==0&&it>3&&(it+nr)%3==1) {
        int rf=(int)ceil(sqrt(e/a->eps)-.4);
        if(rf<1)rf=1;
        braid_StepStatusSetRFactor(s, rf);
    }
    (void)us;
    (void)fs;
    return 0;
}
static int my_Init(braid_App a, double t, braid_Vector*up) {
    my_Vector*u=xm(sizeof(*u));
    u->size=a->nspace;
    u->values=xm((size_t)u->size*sizeof(double));
    u->values[0]=a->xLeft;
    u->values[u->size-1]=a->xRight;
    for(int i=1; i<u->size-1; i++) {
        double x=(double)i/(u->size-1);
        u->values[i]=(1-t/a->tstop)*(10*exp(-(x-.435576)*(x-.435576))-7.27185)+t/a->tstop*(1-x);
    }
    *up=u;
    return 0;
}
static int my_Clone(braid_App a, braid_Vector u, braid_Vector*v) {
    my_Vector*x=xm(sizeof(*x));
    x->size=u->size;
    x->values=xm((size_t)x->size*sizeof(double));
    memcpy(x->values, u->values, (size_t)x->size*sizeof(double));
    *v=x;
    (void)a;
    return 0;
}
static int my_Free(braid_App a, braid_Vector u) {
    free(u->values);
    free(u);
    (void)a;
    return 0;
}
static int my_Sum(braid_App a, double x, braid_Vector u, double y, braid_Vector v) {
    for(int i=0; i<v->size; i++)v->values[i]=x*u->values[i]+y*v->values[i];
    (void)a;
    return 0;
}
static int my_SpatialNorm(braid_App a, braid_Vector u, double*n) {
    double q=0;
    for(int i=0; i<u->size; i++)q+=u->values[i]*u->values[i];
    *n=sqrt(q);
    (void)a;
    return 0;
}
static int my_Access(braid_App a, braid_Vector u, braid_AccessStatus s) {
    braid_Int it;
    braid_AccessStatusGetIter(s, &it);
    if((int)it>a->mgrit_iterations)a->mgrit_iterations=(int)it;
    (void)u;
    return 0;
}
static int my_BufSize(braid_App a, int*n, braid_BufferStatus s) {
    *n=(a->nspace+1)*(int)sizeof(double);
    (void)s;
    return 0;
}
static int my_BufPack(braid_App a, braid_Vector u, void*b, braid_BufferStatus s) {
    double*d=b;
    d[0]=u->size;
    memcpy(d+1, u->values, (size_t)u->size*sizeof(double));
    braid_BufferStatusSetSize(s, (u->size+1)*(int)sizeof(double));
    (void)a;
    return 0;
}
static int my_BufUnpack(braid_App a, void*b, braid_Vector*up, braid_BufferStatus s) {
    double*d=b;
    my_Vector*u=xm(sizeof(*u));
    u->size=(int)d[0];
    u->values=xm((size_t)u->size*sizeof(double));
    memcpy(u->values, d+1, (size_t)u->size*sizeof(double));
    *up=u;
    (void)a;
    (void)s;
    return 0;
}
static int my_Coarsen(braid_App a, braid_Vector f, braid_Vector*cp, braid_CoarsenRefStatus s) {
    int l;
    braid_CoarsenRefStatusGetLevel(s, &l);
    if(l>a->spatial&&f->size>3) {
        my_Vector*c=xm(sizeof(*c));
        c->size=(f->size-1)/2+1;
        c->values=xm((size_t)c->size*sizeof(double));
        c->values[0]=f->values[0];
        c->values[c->size-1]=f->values[f->size-1];
        for(int j=1; j<c->size-1; j++) {
            int i=2*j;
            c->values[j]=.25*f->values[i-1]+.5*f->values[i]+.25*f->values[i+1];
        }
        *cp=c;
        return 0;
    }
    return my_Clone(a, f, cp);
}
static int my_Refine(braid_App a, braid_Vector c, braid_Vector*fp, braid_CoarsenRefStatus s) {
    int l;
    braid_CoarsenRefStatusGetLevel(s, &l);
    if(l>a->spatial&&c->size<a->nspace) {
        my_Vector*f=xm(sizeof(*f));
        f->size=2*(c->size-1)+1;
        f->values=xm((size_t)f->size*sizeof(double));
        for(int j=0; j<c->size; j++)f->values[2*j]=c->values[j];
        for(int j=0; j<c->size-1; j++)f->values[2*j+1]=.5*(c->values[j]+c->values[j+1]);
        *fp=f;
        return 0;
    }
    return my_Clone(a, c, fp);
}
static double max_elapsed(MPI_Comm c, double t) {
    double x=MPI_Wtime()-t, m=0;
    MPI_Reduce(&x, &m, 1, MPI_DOUBLE, MPI_MAX, 0, c);
    return m;
}
int main(int argc, char**argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm c=MPI_COMM_WORLD;
    int rank, np;
    MPI_Comm_rank(c, &rank);
    MPI_Comm_size(c, &np);
    int nt=2049, nx=17, ml=15, mc=64, nu=1, nu0=-1, nuc=1, cf=2, mi=100, skip=0, storage=-1, sc=0, sp=-1;
    double eps=1e-7, tol=1e-9;
    for(int i=1; i<argc; i++) {
        if((!strcmp(argv[i], "-nt")||!strcmp(argv[i], "-ntime"))&&i+1<argc)nt=atoi(argv[++i]);
        else if(!strcmp(argv[i], "-nx")&&i+1<argc)nx=atoi(argv[++i]);
        else if(!strcmp(argv[i], "-eps")&&i+1<argc)eps=atof(argv[++i]);
        else if(!strcmp(argv[i], "-tol")&&i+1<argc)tol=atof(argv[++i]);
        else if(!strcmp(argv[i], "-ml")&&i+1<argc)ml=atoi(argv[++i]);
        else if(!strcmp(argv[i], "-mc")&&i+1<argc)mc=atoi(argv[++i]);
        else if(!strcmp(argv[i], "-nu")&&i+1<argc)nu=atoi(argv[++i]);
        else if(!strcmp(argv[i], "-nu0")&&i+1<argc)nu0=atoi(argv[++i]);
        else if(!strcmp(argv[i], "-nuc")&&i+1<argc)nuc=atoi(argv[++i]);
        else if(!strcmp(argv[i], "-cf")&&i+1<argc)cf=atoi(argv[++i]);
        else if(!strcmp(argv[i], "-mi")&&i+1<argc)mi=atoi(argv[++i]);
        else if(!strcmp(argv[i], "-skip")&&i+1<argc)skip=atoi(argv[++i]);
        else if(!strcmp(argv[i], "-storage")&&i+1<argc)storage=atoi(argv[++i]);
        else if(!strcmp(argv[i], "-sc")&&i+1<argc)sc=atoi(argv[++i]);
        else if(!strcmp(argv[i], "-spatial")&&i+1<argc)sp=atoi(argv[++i]);
        else die("Unknown or incomplete option.");
    }
    my_App*a=calloc(1, sizeof(*a));
    if(!a)die("Allocation failure.");
    a->comm=c;
    a->rank=rank;
    a->nprocs=np;
    a->nspace=nx;
    a->spatial=sp;
    a->tstop=5;
    a->eps=eps;
    a->xLeft=1;
    a->xRight=0;
    alloc_work(a);
    MPI_Barrier(c);
    double total_start=MPI_Wtime();
    braid_Core core;
    braid_Init(c, c, 0, a->tstop, nt, a, my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);
    braid_SetPrintLevel(core, 0);
    braid_SetMaxLevels(core, ml);
    braid_SetNRelax(core, -1, nu);
    if(nu0>=0)braid_SetNRelax(core, 0, nu0);
    braid_SetNRelax(core, ml-1, nuc);
    braid_SetMinCoarse(core, mc);
    braid_SetAbsTol(core, tol);
    braid_SetCFactor(core, -1, cf);
    braid_SetMaxIter(core, mi);
    braid_SetSkip(core, skip);
    braid_SetRefine(core, 1);
    if(storage>=-2)braid_SetStorage(core, storage);
    if(sc) {
        braid_SetSpatialCoarsen(core, my_Coarsen);
        braid_SetSpatialRefine(core, my_Refine);
    }
    MPI_Barrier(c);
    double drive_start=MPI_Wtime();
    braid_Drive(core);
    double drive=max_elapsed(c, drive_start);
    int iterations=0;
    MPI_Reduce(&a->mgrit_iterations, &iterations, 1, MPI_INT, MPI_MAX, 0, c);
    braid_Destroy(core);
    MPI_Barrier(c);
    double total=max_elapsed(c, total_start);
    if(rank==0) {
        printf("%d %d %.14e %.14e %d %d\n", nt, np, drive, total, iterations, sc);
        fflush(stdout);
    }
    free_work(a);
    free(a);
    MPI_Finalize();
    return 0;
}
