#define FILTERED 1
//#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "contact.h
"#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
//#include "log-conform.h"
#include "tension.h"
//#include "view.h"

#include "vtknew_cell.h"

#include "adapt_wavelet_leave_interface.h"

#define Ca 0.6     // Capillary number
#define Re 0.3     // Reynold number
#define We (Ca*Re) // Weber number
#define MUr 1.     // ratio of outer(matrix) to inner(drop) viscosity
#define M 1.       // ratio of outer to inner density
#define Deb 0.4    // Deborah number
#define Beta 0.5   // ratio of the solvent visc. to the total viscoelastic visc.

#define My_PI 3.141592654 // pi
#define xrim 0.0005 // orifice radius
#define xrimsq xrim*xrim // orifice radius squared
#define influx 5.0e-6 // flow rate 

#define avg_vel influx/(2.0*xrim)

int maxlevel = 8;
int    level = 7;
int minlevel = 6;

double tEnd = 50.0;

//scalar mupd[], lam[];

vector h[];
double theta_min = 20.0*pi/180.0;
double theta_max = 160.0*pi/180.0;

h.t[bottom] = contact_angle( (x>-xrim && x<xrim) ? theta_max : theta_min);

p[top] = dirichlet (0.);

u.n[bottom] = ((x>-xrim && x<xrim) ? dirichlet(1.5*avg_vel*(1.0-(x*x)/xrimsq)) : dirichlet(0.0) );
//u.n[bottom] = dirichlet (0.);
u.t[bottom] = dirichlet (0.);
//uf.n[bottom] = 0.;

u.n[left] = dirichlet (0.);
u.n[right] = dirichlet (0.);

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);

f[bottom] = ((x>-xrim && x<xrim) ? 0. : 1.0 );

int main() {
  L0 = 0.01 ;
  size (L0);
  init_grid(1 << maxlevel);
  origin (-L0/2, 0.);

  mu1 = 0.001;
  mu2 = 1.8e-5;
  rho1 = 1000.0;
  rho2 = 1.2;
  f.sigma = 0.072;

  f.height = h;

  TOLERANCE = 1e-4;
  DT = 1.0e-6;
  run();
}

event init (i = 0) {
//  mask (y > 4 ? top : none);

//  refine ( sq(x - 0.01) + sq(y) - sq(0.0022) && level < MAXLEVEL );

//  refine ( (x>0.006 && x<0.014 && y<0.005) && level < MAXLEVEL );
  fraction (f, sq(x) + sq(y) - xrimsq);

//  fraction (f, sq(x) + sq(y-0.0015) - xrimsq);

//  fraction (f, x<0.01 ? sq(x - 0.008) + sq(y-0.0018) - sq(0.002) : sq(x - 0.012) + sq(y-0.0018) - sq(0.002));

//  f.refine = f.prolongation = fraction_refine;
//  restriction ({f}); // for boundary conditions on levels

//  output_ppm (f, linear = true, n=1024, file = "initial_shape.png");
}
/*
event properties (i++) {
  foreach () {
    lam[] = Deb*f[];
    mupd[] = MUr*(1. - Beta)*f[]/Re;
  }
  boundary ({lam, mupd});
}
*/

event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= 9.81;
}

event adapt (i++) {
//  double femax = 5.0e-3;
  double uemax = 5.0e-3;
  adapt_wavelet_leave_interface ({u.x, u.y}, {f},
      (double[]){uemax, uemax}, maxlevel, minlevel, 1);
}


event display_time (i+=100) {
  fprintf(stderr,"%d\t%g\t%g\n", i,dt,t);
}

//event display_image (i+=100) {
//  output_ppm (f, linear = true, n=1024, file = "current_shape.png");
//}

/*
event movie (i+=100)
{
view (quat = {0.000, 0.000, 0.000, 1.000},
      fov = 30, near = 0.01, far = 1000,
      tx = 0.0, ty = -0.454, tz = -2.432,
      width = 888, height = 888);
  clear();
  box ();
  draw_vof ("f");
  save ("movie.mp4");
}

event pictures (i+=1000)
{
  char name[80];
view (quat = {0.000, 0.000, 0.000, 1.000},
      fov = 30, near = 0.01, far = 1000,
      tx = 0.0, ty = -0.454, tz = -2.432,
      width = 888, height = 888);
//  squares ("f", spread = -1, linear = true, map = cool_warm);
  squares ("f", linear = true);
  draw_vof ("f", lc = {1,0,0}, lw = 4);
  box ();
  cells ();
  sprintf (name, "shape-%d.png", i);
  save (name);
}
*/

event write_vtk_new (i+=2000)
{     
  char fname[80];
  sprintf(fname,"fields_%i.vtk",i);
  FILE * fp = fopen(fname,"w");
  output_vtk((scalar*){p, f, u.x, u.y}, fp);  
  fclose(fp);
}

event end (t = tEnd) {}
