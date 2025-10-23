/*****************************************************************
   Soil water dynamics
   Q. Zhuang, February 12, 2003
/***************************************************************** */
using namespace std;
#if !defined(HYDROLOFY_H)
  #include "soilhydrology.hpp"
#endif

// to calculate surface runoff and infiltration occured in the first soil layer
// The implementation is accounting for sub-grid cell resolution
// See Entekhabi and Eagleson, 1989, Bonan, 1996
double HYDM::surfacerunoff(const double& snowmelt, const double& rainthrough, const double& theta1,
       const double& b_layer1, const double& ksat, const double& thetasat, const double& prsat, const double& sthick1)
 {
   // sthick1,  first soil layer thickness (mm)
   // theta1, soil volume moisture at layer 1 (mm3/mm3)
   // ksat, thetasat, presat, are saturation conductivity, voumetric moiture, and wate potential
   // b_layer1, clapp and hornberger constant b

   double dunnerf; // dunne runoff (mm)
   double hortonrf1; // Horton runoff (mm) from the area (sq <s < 1)
   double hortonrf2; // horton runoff (mm) from the area (0<= s < s1)
   double sq; // the value of s at which waters >= fstar
   double v; // d soil potential (pr) / d soil moiture (s) |moisture (s) =1 / d depth (z)
   double s; // soil moisture at first layer
   double a1, a2, a3, a4, a5, a6, a7;
   double srunoff; // surface runoff from the first soil layer (mm day-1 ?)

   // assume the grid cell is homogeneous for the precipitation, snowmelt, and moisture
   s = theta1 / thetasat;
   v = b_layer1 *  prsat / (0.5 * sthick1);

   a1 = snowmelt - ksat * (1.0 - v);
   a2 = ksat * v;

   sq = max( min(a1/a2, 1.0), 0.0);
   a3 = max(rainthrough, powl(1.0,-20));
   a4 = exp(- 1.0 / s);
   a5 = exp(- sq / s);
   a6 = exp(min(80.0, a1/a3));
   a7 = exp(min(80.0, - sq/s - sq*a2/a3));

   if (snowmelt <= ksat) dunnerf = rainthrough * a4 * exp( (snowmelt - ksat)/ a3);
   else dunnerf = (rainthrough + snowmelt - ksat) * a4;

   hortonrf1 = (rainthrough + a1) * (a5 - a4) + a2 * (a4 * (1.0 + s) - a5 * (sq + s) );

  // hortonrf1 = 0.7 * hortonrf1;  // a parameter to control too much runoff

   hortonrf2 = powl(rainthrough,2.0) / (rainthrough + a2*s) * a6 * (1.0 - a7);

   srunoff = dunnerf + hortonrf1 + hortonrf2;

   if (srunoff < 0.0) srunoff = 0.0;
  // adjusted for limiting runoff amount
   srunoff = srunoff * 0.1;

   return srunoff;
 }

// demtermine infiltration from the surface to the soil
 double HYDM::infiltration ( const double& srunoff, const double& snowmelt, const double& rainthrough)
  {
    double infil; // infiltration from surface to the soil (mm day-1)

    infil = (rainthrough + snowmelt) - srunoff;

    return infil;
  }

// Set up soil depth, hydaulic properties
void HYDM::set_soil_pro(const int& cmnt)
{

int n;

/*
 See Clapp and Hornberger, 1978 for the first 11 type of soil
 type 12 and 13 are estimated
 b_const: Clapp and Hornberg function 'b'
 SatK:    Saturated hydraulic conductivity of soil (mm h2o/s)
 SatP:    Saturated matric potential of soil (mm)
 SatV:   Saturated volumetric water content of soil( porosity mm3mm-3)

                                 b_const    SatK          SatP           SatV
 1 Sand                          4.05       0.0176        -121.0         0.395
 2 Loamy sand                    4.38       0.0156        -90.0          0.410
 3 Sandy Loam                    4.90       0.0035        -218.0         0.435
 4 Silt loam                     5.30       0.0007        -786.0         0.485
 5 Loam                          5.39       0.00069       -478.0         0.451
 6 Sandy clay loam               7.12       0.00063       -299.0         0.420
 7 Silty Clay Loam               7.75       0.00017       -356.0         0.477
 8 clay loam                     8.52       0.00025       -630.0         0.476
 9 Sandy clay                    10.4       0.00022       -153.0         0.426
 10 Silty Clay                   10.4       0.00010       -490.0         0.492
 11 Clay                         11.4       0.00013       -405.0         0.482
 12 Peat(organic)                4.00       0.02          -120.0         0.650
 13 Moss,Lichen                  1.00       0.06          -120.0         0.800
 14 Water                        0.0        0.0           0.0            0.00
*/

   nsl = 6;
// define the thickness of each layer (mm)
   dzsoi[0] = thick1[cmnt];
   dzsoi[1] = thick2[cmnt];
   dzsoi[2] = thick3[cmnt];
   dzsoi[3] = thick4[cmnt];
   dzsoi[4] = thick5[cmnt];
   dzsoi[5] = thick6[cmnt];

// initialize soil water content for each layer
   h2osoi[0]= 0.60;
   h2osoi[1]= 0.40;
   h2osoi[2]= 0.40;
   h2osoi[3]= 0.40;
   h2osoi[4]= 0.25;
   h2osoi[5]= 0.25;

 for (n =0; n < 6; n++) {
    // the depth of each layer is defined at the middle of each layer
    if (n == 0) zsoi[n] = 0.5 * dzsoi[n];  // for fist layer
    // the depth for the rest of layers
    else zsoi[n] = zsoi[n-1] + 0.5 * (dzsoi[n-1]+ dzsoi[n]);

    // sign the root bandance for each soil layer
    if (n == 0) root[n] = 0.1; // first layer
    if (n == 1) root[n] = 0.1;
    if (n == 2) root[n] = 0.4;
    if (n == 3) root[n] = 0.3;
    if (n == 4) root[n] = 0.05;
    if (n == 5) root[n] = 0.05; // sixth layer
   }

 // Set up hydraulic properties for each layer according to the constant data array
 // which is from emperical estimates or measuremnt
 for (n =0; n < 6; n++) {
   if (n ==0 ) {
   // first layer will be moss, identification number will be 13, in the array will be [12]
    watsat[n] = thetasat1[cmnt];
    hksat[n]  = ksat1[cmnt];
    smpsat[n] = presat1[cmnt];
    bch[n]    = b1[cmnt];
               }
   if (n ==1 ) {
   // secnd layer we assume as peat identification number will be 12, in the array will be [11];
    watsat[n] = thetasat2[cmnt];
    hksat[n]  = ksat2[cmnt];
    smpsat[n] = presat2[cmnt];
    bch[n]    = b2[cmnt];
        }
   if (n ==2 ) {
   // third layer will be identification number 1, in the array will be [0];
    watsat[n] = thetasat3[cmnt];
    hksat[n]  = ksat3[cmnt];
    smpsat[n] = presat3[cmnt];
    bch[n]    = b3[cmnt];
    }
   if (n ==3 ) {
   // forth layer will be identification number will be 2, in the array will be [1];
    watsat[n] = thetasat4[cmnt];
    hksat[n]  = ksat4[cmnt];
    smpsat[n] = presat4[cmnt];
    bch[n]    = b4[cmnt];
       }
   if (n ==4 ) {
   // fifth layer will be identification number will be 3, in the array will be [2];
    watsat[n] = thetasat5[cmnt];
    hksat[n]  = ksat5[cmnt];
    smpsat[n] = presat5[cmnt];
    bch[n]    = b5[cmnt];
       }
   if (n ==5 ) {
   // sixth layer will be identification number 4, in the array will be [3];
    watsat[n] = thetasat6[cmnt];
    hksat[n]  = ksat6[cmnt];
    smpsat[n] = presat6[cmnt];
    bch[n]    = b6[cmnt];
      }
  }
 // setup time steps
  dthyd = 24 * 3600; // daily time step (sec)
  dtsoi = 60 * 60;    // 60 min. for model runnig step (in seconds)

}

// See Press et al., 1990 for solving tridiagonal system of equations
void HYDM::tridiag(int& nsl, double a[20], double b[20], double c[20], double r[20], double u[20])
 {
  double gam[20];
  double bet;
  int j;

  for (j =0; j < nsl-1; j++)
    {
    if (b[j] == 0.0) b[j] = 0.01;  // for debugging
       bet = b[j];
       u[j] = r[j] / bet;  //printf("hydm 198: j = %i, u = %.3f, r = %.3f, bet = %.3f\n", j, u[j], r[j], bet);
    }
  for ( j = 0; j < nsl; j++) {
      gam[j] = c[j-1] / bet;
      bet = b[j] - a[j] * gam[j];
      u[j] = (r[j] - a[j]*u[j-1]) / bet;
	//printf("hydm 204: j = %i, u = %.3f, r = %.3f, a = %.3f, u-1 = %.3f, bet = %.3f\n", j, u[j], r[j], a[j], u[j-1], bet);
    }
  for (j = nsl-1; j < 0; j--) {
      u[j] = u[j] - gam[j+1] * u[j+1];
	//printf("hydm 208: j = %i, u = %.3f, gam+1 = %.3f, u+1 = %.3f\n", j, u[j], gam[j+1], u[j+1]);
     }

 }

//calculate soil moisture, drainage (sub-surface runoff)
void HYDM::main_soil(int& nsl, double qtran, double qinfl,  double qseva,  double qdrai,
           double dtsoi)
{

int     j,k,n; //  do loop/array indices
int     jimp; //   first impervious soil layer (1 to nsl)
int     iter; //   number of iterations (dthyd/dtsoi)

double  r[20]; // solution matrix for tridiagonal equations
double  a[20]; // "a" vector for tridiagonal matrix
double  b[20]; // "b" vector for tridiagonal matrix
double  c[20]; // "c" vector for tridiagonal matrix
double  dwat[20]; // change in soil water
double  smp[20];  //  soil matrix potential (mm)
double  hk[20];   //   hydraulic conductivity (mm h2o/s)
double  hk2[20];  //  hk**2
double  dsmpdw[20]; // d(smp)/d(h2osoi)
double  dhkdw[20]; //  d(hk)/d(h2osoi)

double  s;      //  h2osoi/watsat
double  hydcon; //  hydraulic conductivity (mm h2o/s)
double  qin;    //  flux of water into soil layer (mm h2o/s)
double  qout;   //  flux of water out of soil layer (mm h2o/s)
double  num;    //  used in calculating qi, qout
double  den;    //  used in calculating qi, qout
double  den2;   //  den**2 used in calculating qi, qout
double  dqidw1; //  d(qin)/d(h2osoi(i-1))
double  dqidw2; //  d(qin)/d(h2osoi(i))
double  dqodw1; //  d(qout)/d(h2osoi(i))
double  dqodw2; //  d(qout)/d(h2osoi(i+1))
double  xs;     //  soil water > sat (mm h2o) or < some minimum (m h2o)
double  axs;    //  amount of xs added to a soil layer
double  x;      //  temporary value of axs (needed for vectorization)
double  newwat; //  temporary value of updated h2osoi

double ftran; // temprary varibles
double finfl;
double fseva;
double fdrai;

// the following three inputs come from canopy and ground surface level processes
  ftran =   qtran / (24.0 * 3600 ) ; // transpiration water flux (mm h2o/s)
  finfl =   qinfl / (24.0 * 3600) ; // infiltration rate (mm h2o/s)
  fseva =   qseva / (24.0 * 3600) ; // ground surface evaporation rate (mm h2o/s)
  fdrai = 4.00 / (24.0 * 3600);  // assume the drainage is 0.0 at the beginning (mm h2o/s)

//printf("hydm 260: ftran = %.3f, finfl = %.3f, qtran = %.3f, qinfl = %.3f\n", ftran, finfl, qtran, qinfl);
// Running for a daily and time step 10 min. total time iteration is 144
for ( iter = 0; iter < int(dthyd/ dtsoi); iter++) {

// initialize jimp and xs
     jimp = nsl-1; // last layer which is the sixth layer
     xs = 0.0;    // water beyond saturation

// make sure the soil water with 1% for debugging purpose
    for ( j = 0; j < nsl; j++) {
      if (h2osoi[j] <= 0.0) h2osoi[j] = 0.01;
      }

// evaluate hydraulic conductivity, soil matrix potential for each layer
// d(smp)/d(h2osoi), and d(hk)/d(h2osoi)
     for (j = 0; j < nsl; j++) {
     	   s = max(h2osoi[j]/watsat[j],0.05);
         smp[j] = smpsat[j] * pow(s,-bch[j]);  // smpsat[j] could be estimated using soil texture
         dsmpdw[j] = -bch[j]*smp[j]/h2osoi[j];

// here should be coupled with soil thermal model
   //      if (tsoil[j] > 0.0)   // no freezing
         hydcon = hksat[j] * powl(s,2.*bch[j]+3.); // hksat[j] shoud be estimated wih soil texture
   //      else hydcon = 0.0;  // freezing prevent the soil recharge
         hk[j] = max(hydcon,1.e-10);
         dhkdw[j] = hydcon*(2.*bch[j]+3.)/h2osoi[j];
         hk2[j] = hk[j]*hk[j];
       }

// find first impervious layer basing on soil thermal model
       for (j = nsl-1; j > -1; j--)
        {
         if ( hk[j] <= 1.e-10) jimp = j;
         else jimp = nsl-1;  // the last soil layer
        }

// set up r, a, b, and c vectors for tridiagonal solution for the first soil layer
       j = 0;
       num = -2.*(smp[j]-smp[j+1]) - (dzsoi[j]+dzsoi[j+1]);
       den = dzsoi[j]/hk[j] + dzsoi[j+1]/hk[j+1];
       den2 = den*den;
       qout = num/den;
       dqodw1 = (-2.*den*dsmpdw[j] + num*dzsoi[j]/hk2[j]*dhkdw[j]) / den2;
       dqodw2 = ( 2.*den*dsmpdw[j+1] + num*dzsoi[j+1]/hk2[j+1]*dhkdw[j+1]) / den2;
       r[j] = (fseva + ftran * root[j]) - finfl - qout;
	//printf("hydm 305: j = %i, r = %.3f, fseva = %.3f, ftran = %.3f, root = %.3f, finfl = %.3f, qout = %.3f\n", j, r[j], fseva, ftran, root[j], finfl, qout);
       a[j] = 0.0;
       b[j] = dqodw1 - dzsoi[j]/dtsoi;
       c[j] = dqodw2;

// for the last soil layer
       j = nsl-1;
       num = -2.*(smp[j-1]-smp[j]) - (dzsoi[j-1]+dzsoi[j]);
       den = dzsoi[j-1]/hk[j-1] + dzsoi[j]/hk[j];
       den2 = den*den;
       qin = num/den;
       dqidw1 = (-2.* den * dsmpdw[j-1] + num * dzsoi[j-1]/hk2[j-1]*dhkdw[j-1]) / den2;
       dqidw2 = ( 2.* den * dsmpdw[j] + num * dzsoi[j]/hk2[j]*dhkdw[j]) / den2;
       qout = -hk[j];
       dqodw1 = -dhkdw[j];
       r[j] = ftran * root[j] + qin - qout;
	//printf("hydm 320: j = %i, r = %.3f, ftran = %.3f, root = %.3f, qin = %.3f, qout = %.3f\n", j, r[j], ftran, root[j], qin, qout);
       a[j] = -dqidw1;
       b[j] = dqodw1 - dqidw2 - dzsoi[j]/dtsoi;
       c[j] = 0.0;

// for layer from second (j=1) to 4th layer
	    for (j = 1; j < nsl-1; j++) {
	      num = -2.*(smp[j-1]-smp[j]) - (dzsoi[j-1]+dzsoi[j]);
         den = dzsoi[j-1]/hk[j-1] + dzsoi[j]/hk[j];
         den2 = den*den;
         qin = num/den;
         dqidw1 = (-2.*den*dsmpdw[j-1] + num*dzsoi[j-1]/hk2[j-1]*dhkdw[j-1]) / den2;
         dqidw2 = ( 2.*den*dsmpdw[j] + num*dzsoi[j]/hk2[j]*dhkdw[j]) / den2;
         num = -2.*(smp[j]-smp[j+1]) - (dzsoi[j]+dzsoi[j+1]);
         den = dzsoi[j]/hk[j] + dzsoi[j+1]/hk[j+1];
         den2 = den*den;
         qout = num/den;
         dqodw1 = (-2.*den*dsmpdw[j] + num*dzsoi[j]/hk2[j]*dhkdw[j]) / den2;
         dqodw2 = ( 2.*den*dsmpdw[j+1] + num*dzsoi[j+1]/hk2[j+1]*dhkdw[j+1])/ pow(den,2.0);
         r[j] = ftran*root[j] + qin - qout;
	//printf("hydm 341: j = %i, r = %.3f, ftran = %.3f, root = %.3f, qin = %.3f, qout = %.3f\n", j, r[j], ftran, root[j], qin, qout);
         a[j] = -dqidw1;
         b[j] = dqodw1 - dqidw2 - dzsoi[j]/dtsoi;
         c[j] = dqodw2;
	    }

// solve for dwat with a, b, c,and r
      tridiag(nsl ,a ,b ,c ,r ,dwat);

// update water constraining h2osoi <= watsat. accumulate excess water
      for ( j = 0; j < nsl; j++) {
         newwat = h2osoi[j] + dwat[j];
	//printf("hydm 348: newwat = %.3f, j = %i, h2osoi = %.3f, dwat = %.3f\n", newwat, j, h2osoi[j], dwat[j]);
         xs = xs + max(newwat-watsat[j],0.0) * dzsoi[j];  // estimate extra water
         h2osoi[j] = min(watsat[j],newwat); // water content was constraining less than saturation
	//printf("hydm 350: j = %i, h2osoi = %.3f, watsat = %.3f, newwat = %.3f\n", j, h2osoi[j], watsat[j], newwat);
       }

// add excess water back to soil and bring soil layers up to saturation.
// only do if soil layer is less than or equal to the impervious layer
      for (j = 0; j< nsl; j++) {
        if (j  <= jimp)
           x = min((watsat[j]-h2osoi[j])*dzsoi[j],xs);
        else  x = 0.0;
        axs = x;
        xs = xs - axs;
        h2osoi[j] = h2osoi[j] + axs/dzsoi[j];
	//printf("hydm 362: j = %i, h2osoi = %.3f, axs = %.3f, dasoi = %.3f\n", j, h2osoi[j], axs, dzsoi[j]);
    	}

// sub-surface drainage (accumulate over dthyd/dtsoi iterations)
      fdrai = fdrai + xs + (hk[nsl-1] + dhkdw[nsl-1] * dwat[nsl-1]) * dtsoi;
} // the end of iterations

// sub-surface drainage over time step = dthyd
       fdrai = fdrai/dthyd;

// limit h2osoi >= 0.01. get water needed to bring h2osoi = 0.01 from lower layer.
      for ( j = 0; j < nsl ; j++) {
         if (h2osoi[j] < 0.01) xs = (0.01-h2osoi[j])*dzsoi[j];
         else  xs = 0.0;
         h2osoi[j] = h2osoi[j] + xs / dzsoi[j];
	//printf("hydm 377: j = %i, h2osoi = %.3f, xs = %.3f, dzsoi = %.3f\n", j, h2osoi[j], xs, dzsoi[j]);
        }

  //printf(" %3.2f %3.2f  %3.2f %3.2f %3.2f %3.2f \n", h2osoi[0], h2osoi[1],h2osoi[2],h2osoi[3],h2osoi[4],h2osoi[5]);

  }


//Interpolating moisture to obtain every 1 cm depth soil moistures with Lagrange method
// Numerical method of computation and FORTRAN lanuage, Xianyi Zheng, 1986
double HYDM::InterpolatSM(const double& vsm1, const double& vsm2,const double& vsm3,
                const double& vsm4, const double& vsm5, const double& vsm6,  double x, const int& cmnt)
{
  double f, p;
  int i, j;
  double x0[6], y0[6];

    x0[0] = thick1[cmnt] /10.0;
    x0[1] = thick2[cmnt] /10.0 + x0[0];
    x0[2] = thick3[cmnt] /10.0 + x0[1];
    x0[3] = thick4[cmnt] /10.0 + x0[2];
    x0[4] = thick5[cmnt] /10.0 + x0[3];
    x0[5] = thick6[cmnt] /10.0 + x0[4];

    y0[0] = vsm1;
    y0[1] = vsm2;
    y0[2] = vsm3;
    y0[3] = vsm4;
    y0[4] = vsm5;
    y0[5] = vsm6;

  f = 0.0;
  for (i = 0; i < 6; i++)
  {
   p = 1.0;
   for (j =0 ; j < 6; j++)
   {
    if (i!=j) {
               p = p * (x - x0[j]) / (x0[i] - x0[j]);
              }

   }
   f = f + p * y0[i];
  }
 return f;

}


// Strictly use Granberg et al., 1999
void HYDM::wetlandTAWA(const double& rainthrough, const double& snowmelt,const double& qtrans,
                      const double& qseva, const int& dm, const int& m, const double& soilt,
                      const double& sand, const double& silt, const double& clay, const double& qdraimax)
{

 int i;
 double A, B, C;
 double wt, wt1, wt2; // water table depth
 double fcoarse;
 double q, q1,z;
 double qdrai;
 double thetas; // water content in the soil profile;
 double az; // the gradient in the linearly decreasing interval

 double phi = 0.9; // porosity of soil; = 0.9 cm-3cm-3 (Frolking 1996)
 double zm = 300; // maximum water table depth, assumed as 30 cm as Frolking, 1996, Granberg et al. 1999
 double thetamin = 0.42; // minimum water content at the unstaturated zone (= 0.25)
 double zmin = 100; // minimum depth with the minimum water content ( = 0.1 m = 100mm)
 const double pvsand = 0.45; // relative volume of coarse pores in sandy soils
 const double pvsilt = 0.20; // silt soils
 const double pvclay = 0.14; // clay soils

 // qdrai = qdraimax * (fcoarse/maximum) See Walter et al., 2001;
 //here modeled as the product of maximum water and water holding capacity as measure by porce space of the soil
 fcoarse = (sand * pvsand + silt * pvsilt + clay * pvclay) * 0.01;
 qdrai = qdraimax * fcoarse;

 q = rainthrough + snowmelt - qtrans - qseva - qdrai;
 az = (phi - thetamin) / zmin;

 wt1 = 1.5 * (phi * zm - q) / (phi - thetamin);
 wt2 =  sqrt( 1.5 * (phi * zm - q) / az);

// if (wt1 <= zmin) wt = wt2;
// else wt = wt1;

 wt  = min (wt1, wt2);

 if (wt > 300.0) wt = 300.0;

 if (soilt < 0.0) wt = 100.;

 watertable[dm][m] = wt;

// deterine the water content above water table, similar to Granberg et al., 1999
// assume water content linearly decrease with depth
// the soil moiture is evaluated for every 10 mm, i.e. every 1 cm
 thetas = max (thetamin, (phi - az * wt));
 int d = int (ceil((300 - wt) /10.0));
 for (i = 0; i < d; i++)
  {
   z = 10 * i;
   q1 = pow( (z / wt), 2.0);
   unsatthetaWL[i] = min(phi,(thetas + (phi - thetas) * q1) );
  }

}

// get soil ecd
void HYDM::getecdsoil(char ecd[80])
{
  getsoilhydpro(ecd);

}

void HYDM::getecdsoil(ofstream& rflog1)
{

  char ecd[80];

  cout << "Enter name of the file with soil hydraulic parameter values (.ECD):" << endl;
  //cin >> ecd;
  fpara >> ecd;
  
  rflog1 << "Enter name of the file with soil hydaulic parameter values (.ECD):";
  rflog1 << ecd << endl << endl;

  getsoilhydpro(ecd);

};

void HYDM::getsoilhydpro(char ecd[80])
{

  const int NUMVAR = 34;
  char dummy[35][40];
  ifstream infile;
  int i;
  int dcmnt;

  long update[35];

  infile.open(ecd, ios::in);

  if (!infile)
  {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (dcmnt = 1; dcmnt < 35; dcmnt++)
  {

  infile >> dummy[dcmnt] >> dummy[dcmnt]>> thick1[dcmnt] >> b1[dcmnt]>> ksat1[dcmnt] >> thetasat1[dcmnt] >> presat1[dcmnt]
         >> thick2[dcmnt] >> b2[dcmnt]>> ksat2[dcmnt] >> thetasat2[dcmnt] >> presat2[dcmnt]
         >> thick3[dcmnt] >> b3[dcmnt]>> ksat3[dcmnt] >> thetasat3[dcmnt] >> presat3[dcmnt]
         >> thick4[dcmnt] >> b4[dcmnt]>> ksat4[dcmnt] >> thetasat4[dcmnt] >> presat4[dcmnt]
         >> thick5[dcmnt] >> b5[dcmnt]>> ksat5[dcmnt] >> thetasat5[dcmnt] >> presat5[dcmnt]
         >> thick6[dcmnt] >> b6[dcmnt]>> ksat6[dcmnt] >> thetasat6[dcmnt] >> presat6[dcmnt] >> qdraimax[dcmnt] >> update[dcmnt];
//  printf("hydm dcmnt = %i\n", dcmnt);
   }

  infile.close();

};


void HYDM::getecdhydpro(ofstream& rflog1) {

  const int NUMVAR = 34;
  char dummy[NUMVAR][8];
  ifstream infile;
  int i;
  int dcmnt;
  char ecd[20];

  int update[12];

  infile.open(ecd, ios::in);

  cout << "Enter name of soil hydraulic (.ECD) data file with parameter values:" << endl;
  //cin >> ecd;

  rflog1 << "Enter name of soil hydraulic(.ECD) data file with parameter values:" << ecd << endl << endl;

  infile.open(ecd, ios::in);

  if (!infile) {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (dcmnt = 1; dcmnt < 35+1; dcmnt++)
  {
  infile >> dummy[0] >> dummy[0]>> thick1[dcmnt] >> b1[dcmnt]>> ksat1[dcmnt] >> thetasat1[dcmnt] >> presat1[dcmnt]
         >> thick2[dcmnt] >> b2[dcmnt]>> ksat2[dcmnt] >> thetasat2[dcmnt] >> presat2[dcmnt]
         >> thick3[dcmnt] >> b3[dcmnt]>> ksat3[dcmnt] >> thetasat3[dcmnt] >> presat3[dcmnt]
         >> thick4[dcmnt] >> b4[dcmnt]>> ksat4[dcmnt] >> thetasat4[dcmnt] >> presat4[dcmnt]
         >> thick5[dcmnt] >> b5[dcmnt]>> ksat5[dcmnt] >> thetasat5[dcmnt] >> presat5[dcmnt]
         >> thick6[dcmnt] >> b6[dcmnt]>> ksat6[dcmnt] >> thetasat6[dcmnt] >> presat6[dcmnt] >> qdraimax[dcmnt] >> update[dcmnt];


  }


  infile.close();

};

void HYDM::wetland_soil(int& nsl, double qtran, double qinfl,  double qseva,  double qdrai,
           double dtsoi)
{

int     j,k,n; //  do loop/array indices
int     jimp; //   first impervious soil layer (1 to nsl)
int     iter; //   number of iterations (dthyd/dtsoi)

double  r[20]; // solution matrix for tridiagonal equations
double  a[20]; // "a" vector for tridiagonal matrix
double  b[20]; // "b" vector for tridiagonal matrix
double  c[20]; // "c" vector for tridiagonal matrix
double  dwat[20]; // change in soil water
double  smp[20];  //  soil matrix potential (mm)
double  hk[20];   //   hydraulic conductivity (mm h2o/s)
double  hk2[20];  //  hk**2
double  dsmpdw[20]; // d(smp)/d(h2osoi)
double  dhkdw[20]; //  d(hk)/d(h2osoi)

double  s;      //  h2osoi/watsat
double  hydcon; //  hydraulic conductivity (mm h2o/s)
double  qin;    //  flux of water into soil layer (mm h2o/s)
double  qout;   //  flux of water out of soil layer (mm h2o/s)
double  num;    //  used in calculating qi, qout
double  den;    //  used in calculating qi, qout
double  den2;   //  den**2 used in calculating qi, qout
double  dqidw1; //  d(qin)/d(h2osoi(i-1))
double  dqidw2; //  d(qin)/d(h2osoi(i))
double  dqodw1; //  d(qout)/d(h2osoi(i))
double  dqodw2; //  d(qout)/d(h2osoi(i+1))
double  xs;     //  soil water > sat (mm h2o) or < some minimum (m h2o)
double  axs;    //  amount of xs added to a soil layer
double  x;      //  temporary value of axs (needed for vectorization)
double  newwat; //  temporary value of updated h2osoi

double ftran; // temprary varibles
double finfl;
double fseva;
double fdrai;

// the following three inputs come from canopy and ground surface level processes
  ftran =   qtran / (24.0 * 3600 ) ; // transpiration water flux (mm h2o/s)
  finfl =   qinfl / (24.0 * 3600) ; // infiltration rate (mm h2o/s)
  fseva =   qseva / (24.0 * 3600) ; // ground surface evaporation rate (mm h2o/s)
  fdrai = 4.00 / (24.0 * 3600);  // assume the drainage is 0.0 at the beginning (mm h2o/s)


// Running for a daily and time step 10 min. total time iteration is 144
for ( iter = 0; iter < int(dthyd/ dtsoi); iter++) {

// initialize jimp and xs
     jimp = nsl-1; // last layer which is the sixth layer
     xs = 0.0;    // water beyond saturation

// make sure the soil water with 1% for debugging purpose
    for ( j = 0; j < nsl; j++) {
      if (h2osoi[j] <= 0.0) h2osoi[j] = 0.01;
      }

// evaluate hydraulic conductivity, soil matrix potential for each layer
// d(smp)/d(h2osoi), and d(hk)/d(h2osoi)
     for (j = 0; j < nsl; j++) {
     	   s = max(h2osoi[j]/watsat[j],0.05);
         smp[j] = smpsat[j] * pow(s,-bch[j]);  // smpsat[j] could be estimated using soil texture
         dsmpdw[j] = -bch[j]*smp[j]/h2osoi[j];

// here should be coupled with soil thermal model
   //      if (tsoil[j] > 0.0)   // no freezing
         hydcon = hksat[j] * powl(s,2.*bch[j]+3.); // hksat[j] shoud be estimated wih soil texture
   //      else hydcon = 0.0;  // freezing prevent the soil recharge
         hk[j] = max(hydcon,1.e-10);
         dhkdw[j] = hydcon*(2.*bch[j]+3.)/h2osoi[j];
         hk2[j] = hk[j]*hk[j];
       }

// find first impervious layer basing on soil thermal model
       for (j = nsl-1; j > -1; j--)
        {
         if ( hk[j] <= 1.e-10) jimp = j;
         else jimp = nsl-1;  // the last soil layer
        }

// set up r, a, b, and c vectors for tridiagonal solution for the first soil layer
       j = 0;
       num = -2.*(smp[j]-smp[j+1]) - (dzsoi[j]+dzsoi[j+1]);
       den = dzsoi[j]/hk[j] + dzsoi[j+1]/hk[j+1];
       den2 = den*den;
       qout = num/den;
       dqodw1 = (-2.*den*dsmpdw[j] + num*dzsoi[j]/hk2[j]*dhkdw[j]) / den2;
       dqodw2 = ( 2.*den*dsmpdw[j+1] + num*dzsoi[j+1]/hk2[j+1]*dhkdw[j+1]) / den2;
       r[j] = (fseva + ftran * root[j]) - finfl - qout;
       a[j] = 0.0;
       b[j] = dqodw1 - dzsoi[j]/dtsoi;
       c[j] = dqodw2;

// for the last soil layer
       j = nsl-1;
       num = -2.*(smp[j-1]-smp[j]) - (dzsoi[j-1]+dzsoi[j]);
       den = dzsoi[j-1]/hk[j-1] + dzsoi[j]/hk[j];
       den2 = den*den;
       qin = num/den;
       dqidw1 = (-2.* den * dsmpdw[j-1] + num * dzsoi[j-1]/hk2[j-1]*dhkdw[j-1]) / den2;
       dqidw2 = ( 2.* den * dsmpdw[j] + num * dzsoi[j]/hk2[j]*dhkdw[j]) / den2;
       qout = -hk[j];
       dqodw1 = -dhkdw[j];
       r[j] = ftran * root[j] + qin - qout;
       a[j] = -dqidw1;
       b[j] = dqodw1 - dqidw2 - dzsoi[j]/dtsoi;
       c[j] = 0.0;

// for layer from second (j=1) to 4th layer
	    for (j = 1; j < nsl-1; j++) {
	      num = -2.*(smp[j-1]-smp[j]) - (dzsoi[j-1]+dzsoi[j]);
         den = dzsoi[j-1]/hk[j-1] + dzsoi[j]/hk[j];
         den2 = den*den;
         qin = num/den;
         dqidw1 = (-2.*den*dsmpdw[j-1] + num*dzsoi[j-1]/hk2[j-1]*dhkdw[j-1]) / den2;
         dqidw2 = ( 2.*den*dsmpdw[j] + num*dzsoi[j]/hk2[j]*dhkdw[j]) / den2;
         num = -2.*(smp[j]-smp[j+1]) - (dzsoi[j]+dzsoi[j+1]);
         den = dzsoi[j]/hk[j] + dzsoi[j+1]/hk[j+1];
         den2 = den*den;
         qout = num/den;
         dqodw1 = (-2.*den*dsmpdw[j] + num*dzsoi[j]/hk2[j]*dhkdw[j]) / den2;
         dqodw2 = ( 2.*den*dsmpdw[j+1] + num*dzsoi[j+1]/hk2[j+1]*dhkdw[j+1])/ pow(den,2.0);
         r[j] = ftran*root[j] + qin - qout;
         a[j] = -dqidw1;
         b[j] = dqodw1 - dqidw2 - dzsoi[j]/dtsoi;
         c[j] = dqodw2;
	    }

// solve for dwat with a, b, c,and r
      tridiag(nsl ,a ,b ,c ,r ,dwat);

// update water constraining h2osoi <= watsat. accumulate excess water
      for ( j = 0; j < nsl; j++) {
         newwat = h2osoi[j] + dwat[j];
         xs = xs + max(newwat-watsat[j],0.0) * dzsoi[j];  // estimate extra water
         h2osoi[j] = min(watsat[j],newwat); // water content was constraining less than saturation
       }

// add excess water back to soil and bring soil layers up to saturation.
// only do if soil layer is less than or equal to the impervious layer
      for (j = 0; j< nsl; j++) {
        if (j  <= jimp)
           x = min((watsat[j]-h2osoi[j])*dzsoi[j],xs);
        else  x = 0.0;
        axs = x;
        xs = xs - axs;
        h2osoi[j] = h2osoi[j] + axs/dzsoi[j];
    	}

// sub-surface drainage (accumulate over dthyd/dtsoi iterations)
      fdrai = fdrai + xs + (hk[nsl-1] + dhkdw[nsl-1] * dwat[nsl-1]) * dtsoi;
      fdrai = 0.0;  // set the drainage to zero
} // the end of iterations

// sub-surface drainage over time step = dthyd
       fdrai = fdrai/dthyd;

// limit h2osoi >= 0.01. get water needed to bring h2osoi = 0.01 from lower layer.
      for ( j = 0; j < nsl ; j++) {
         if (h2osoi[j] < 0.01) xs = (0.01-h2osoi[j])*dzsoi[j];
         else  xs = 0.0;
         h2osoi[j] = h2osoi[j] + xs / dzsoi[j];
        }

   //  printf(" %3.2f %3.2f  %3.2f %3.2f %3.2f %3.2f \n", h2osoi[0], h2osoi[1],h2osoi[2],h2osoi[3],h2osoi[4],h2osoi[5]);

  };
  
  
