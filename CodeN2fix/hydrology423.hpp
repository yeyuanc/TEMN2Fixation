/* **************************************************************
*****************************************************************
HYDROLOGY423.HPP - object describes physical characteristics of the
	      atmosphere

           - created by QZ 06162000
           - change to 1-plant functional type, 01/11/2002 Q. Zhuang
*****************************************************************
************************************************************** */
using namespace std;
// Class uses global constants CYCLE and MAXRTIME

class Hydrology
{

  public:

/* **************************************************************
		 Public Functions
************************************************************** */
    Hydrology();  // initialization

    void    showpenmon(int& dcmnt);
    void    getecd(char ecd[80]);
    void    getecd(ofstream& rflog1);
    void    getpenmonecd(ofstream& rflog1);
    void    getpenmonecd(char ecd[80]);

    void   drad(const double& LAI, const double& EXT, const double& nirr, const int& dm, const int& m);
	 double LWP(const double& SOILH2O, const double& SOILCAP, const double sand,
              const double clay, const double silt);
    double CanopyConductance(const double& LAI, const double& gmax,
           const double& airt, const double& vpd, const double& leafwaterp,
           const double& vpd_close,const double& vpd_open, const double& psi_close,
           const double& psi_open);
    double densityA(const double& airt);
    double latentV(const double& airt);
    double penmanmonteith_new(const double& airtemp, double drad, const double& rule,
           const double& vpd, const double& lamuda, const double& Gs, const double& LAI);
    double intercept(const double& rain, const double& inter_coef, const int& dm, const int& m);
    void   vpd(double& vaporP, const double& airt);
    double evaporation_layer1(const double& airtemp, const double& pa,  double& lanmuda, double& vapord, const double& soilrad);
    double snowsublimation_ground( double snowwatereq, const double& airtemp, double surfacerad);
    double snowmelt_srad(double snow_rad, const double airt, const double snowthrough);
    double sublimation_canopy(double lai, double canopy_rad, const double airt, const double snowfall,
           const int& dm, const int& m);

/* **************************************************************
		 Public Variables
************************************************************** */

   double yrhevap;
   double yrhtrans;
   double yrsevap;
   double yrsnowsub;
   double yrsubcan;
 //  pressure Variables

   double vpdeficit; // kpa
   double absoluteHD; // g/m3
   double vaporpressure[CYCLE][31];
   double snowpack_canopy[CYCLE][31]; // mm/month
// for soil surface evaporation
   double sub_from_canopy[CYCLE][31]; // sublimation from canopy mm/day
   double snowsub[CYCLE][31]; // snow sublimation mm/day (water equivalent) from ground
   double soil_evap[CYCLE][31]; // actual evaporation of soil surface mm/day
   double ratio;
   double potevap_surface[CYCLE][31]; // potential soil surface evaporation mm/day
   int dayssincerain; // days since rain
   int sublimationdays[CYCLE][31]; // canopy sublimation days
// end of soil surface

// Evaporation Variables // added by QZ 14 June 2000
   double hdrad[CYCLE][31]; // canopy average daily radiation  (MJ m-2 day-1)
   double hLWP[CYCLE][31];  // daily maximum leaf water potential (MPa)
   double hcch[CYCLE][31]; // canopy conductance with humidity reduction (ms-1)
   double hslope[CYCLE][31]; // slope of the saturation vapor pressure curve
   double hpa[CYCLE][31]; // air density
   double hlv[CYCLE][31]; // latent heat of vapor
   double htrans[CYCLE][31]; // transpiration (mmday-1)
   double hevap[CYCLE][31];  // Estimated Actual Evaporation by canopy and bare soil (?)(mm)
   double surfacerad[CYCLE][31]; // soil surface radiation MJ m-2 day-1
   double Gs[CYCLE][31]; // canopy conductance
//   Climatic Variables
   double  rainthrough[CYCLE][31];       // Rainfall to soil ground (mm)
   double  canopy_water[CYCLE][31]; // rain on canopy (mm) // added by Qianlai zhuang 14 June 2000
   double  snowthroughfall[CYCLE][31]; // store snow throughfall mm
   double  prec_to_canopy[CYCLE][31];
   double snowmelt[CYCLE][31];
   double inter_coef[35], hvpd[35], EXT[35],dewpoint[CYCLE][31];
   double CCmax[35],LWPc[35],SLOPEcc[35];
   double hle,LWPmin[35];
   double gmax[35]; // maximum canopy conductance
   double vpd_close[35];
   double vpd_open[35];
   double psi_close[35];
   double psi_open[35];


  // added for hydrology model
   int vapflag;

   double  vap[CYCLE][31];       // Total vap
   double  tvap[MAXRTIME][CYCLE][31];
   long    vapyear[MAXRTIME];       // year of vap  data

   double yrvap;             // Annual sum of total vap
   double yrtvap[MAXRTIME];

   double mxvap;             // Maximum monthly vap
   double mxtvap[MAXRTIME];

   double daze[CYCLE];          // number of days per month

  private:
/* **************************************************************
		      Private Variables
************************************************************** */

};

