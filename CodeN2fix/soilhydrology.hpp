
/*********************************************************
**********************************************************
hydrology.HPP -object describing Hydrology models
changed to 1-functional type Q. Zhuang 01/11/2002

Q. Zhuang, Change the function of transpiration spliting, 05/Jan/2003
**********************************************************
*********************************************************/
using namespace std;
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

class HYDM {

	public:

/******************************************************
			public Functions
******************************************************/
    void    getecdsoil(ofstream& rflog1);
    void    getecdsoil(char ecd[80]);
    void    getecdhydpro(ofstream& rflog1);
    void    getsoilhydpro(char ecd[80]);


//    void    getecdsoiltexture(ofstream& rflog1);

//    void    getecdsoiltexture(char ecd[80]);
//    void    getsoiltexture(ofstream& rflog1);
//    void    getsoiltexture(char ecd[80]);


    double surfacerunoff(const double& snowmelt, const double& rainthrough, const double& theta1,
          const double& b_layer1, const double& ksat, const double& thetasat, const double& prsat, const double& sthick1);
    double infiltration ( const double& srunoff, const double& snowmelt, const double& rainthrough);
    void set_soil_pro(const int& cmnt);
    void tridiag(int& nsl, double a[20], double b[20], double c[20], double r[20], double u[20]);
//    void main_soil(int& nsl, const double& qtran, const double& qinfl, const double& qeva, const double& qdrai, double& dtsoi);
    void main_soil(int& nsl, double qtran,  double qinfl,  double qeva,  double qdrai, double dtsoi);
    void wetland_soil(int& nsl, double qtran,  double qinfl,  double qeva,  double qdrai, double dtsoi);
    double InterpolatSM(const double& vsm1, const double& vsm2,const double& vsm3,
                const double& vsm4, const double& vsm5, const double& vsm6,  double x, const int& cmnt);

    void wetlandTAWA(const double& rainthrough, const double& snowmelt,const double& qtrans,
                      const double& qseva, const int& dm, const int& m, const double& soilt,
                      const double& sand, const double& silt, const double& clay, const double& qdraimax);
/**************************************************************
			Public Variables
**************************************************************/

   double thick1[35], thick2[35],thick3[35],thick4[35],thick5[35],thick6[35]; //layers thickness mm
   double b1[35],b2[35],b3[35],b4[35],b5[35],b6[35];  // b constant
   double ksat1[35],ksat2[35],ksat3[35],ksat4[35],ksat5[35],ksat6[35];  // saturation conductivity
   double thetasat1[35],thetasat2[35],thetasat3[35],thetasat4[35],thetasat5[35],thetasat6[35]; // saturation volumetric moisture
   double presat1[35],presat2[35],presat3[35],presat4[35],presat5[35],presat6[35]; // saturation potential

   double qdraimax[35];

   double  watsat[20]; // saturated volumetric soil water content (porosity)
   double  smpsat[20]; // soil matrix potential at saturation (mm)
   double  bch[20];     // clapp and hornberger "b"
   double  hksat[20];   // hydraulic conductivity at saturation (mm h2o/s)

   double  tsoi[20];    // soil temperature (degree k)
   double  root[20];    // relative root abundance (0 to 1)
   double  dzsoi[20];   // soil layer thickness (m)
   double  zsoi[20];
   double  dtsoi;  // time step (e.g. 10 min)
   double  dthyd; //daily
   double  h2osoi[20];   // volumetric soil water content (mm-3/mm-3)
   int nsl;  // soil layers

   double surrunoff[CYCLE][31]; // run off from surface (mm)
   double qtran; // transpiration water flux (mm h2o/s)
   double surinfl[CYCLE][31]; // infiltration rate (mm h2o/s)
   double qseva; // ground surface evaporation rate (mm h2o/s)
   double bdrai[CYCLE][31]; // assume the drainage is 0.0 at the beginning (mm h2o/s)

   double watertable[CYCLE][31]; // water depth for each day of a month
   double unsatthetaWL[10]; // to store the water content for layers of wetland (every 10 mm)

   double intersoilm[6500]; // to store the soil moisture for each cm
};
