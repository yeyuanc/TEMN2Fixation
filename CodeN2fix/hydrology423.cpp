
/* ****************************************************************
HYDROLOGY.CPP - object describes physical and biological processes of water
    dyanmics
 - Created by QZ for hydrology model 06162000
 - changed to 1-plant functional type 01/11/2002
 - Changed to daily version of the model, 05/Jan/2003, Q. Zhuang

************************************************************** */
using namespace std;
#if !defined(HYDROLOFY423_H)
  #include "hydrology423.hpp"
#endif


Hydrology::Hydrology() {

// Initialize number of days per month for each month of year

  daze[0] = daze[2] = daze[4] = daze[6] = daze[7] = 31.0;
  daze[9] = daze[11] = 31.0;
  daze[3] = daze[5] = daze[8] = daze[10] = 30.0;
  daze[1] = 28.0;

};


// Interception at canopy for evaporation
// interception rate is obtained from Helvey (1971) and helvey and Patric (1965)
double Hydrology::intercept(const double& rain, const double& inter_coef,
       const int& dm, const int& m)
{

   double max_int; // maximum daily canopy interception mm rain in the canopy */
   double tempuse;

   tempuse = rain /10.0 ; //convert to  rainfall from mm to cm
 	max_int = tempuse - ((tempuse*0.77-0.05)+ 0.02 * tempuse);
   max_int = max_int * 10; // covert from cm back to mm

   if (rain  <= max_int)          /* all intercepted */
		{
			prec_to_canopy[dm][m] = rain;
			rainthrough[dm][m] = 0.0;
		}
   else                          /* canopy limits interception */
		{
			prec_to_canopy[dm][m] = max_int;
			rainthrough[dm][m] = rain - max_int;
		}

      return rainthrough[dm][m];
 };


// Split the surafce radiation to canopy and ground surface
void Hydrology::drad(const double& LAI, const double& EXT, const double& nirr,
      const int& dm, const int& m)

 {
  // nirr[CYCLE];        // Net Irradiance   (cal/(sq. cm * day))
  // hdrad (MJ m-2 day-1) canopy daily average radiation
  // LAI (m2 m-2)  leaf area index
  // EXT dimentionless  extinction coefficient of radiation through the canopy
  // 2.2 is the coefficiet changig total LAI to projected LAI (dimentionless)
  double drad;
  double nirrt;
  double lai;
  // solar radiation to the canopy hdrad[dm][m]
  // following Running and Coughlan 1988
  // assume 50% is absorbed by canopy See P.31 Table 2.4., Jones,1994
     lai = LAI;
     if (lai < 0.1) {lai = 0.1;}
     drad = nirr * 0.5;
     drad = (drad * (1-exp(-(lai/2.2) * EXT))) / (EXT * lai / 2.2);
     //printf("hyd 75: drad = %.3f, lai = %.3f, ext = %.3f\n", drad, lai, EXT);
     nirrt = drad / 0.2388; // J /cm2 /d
     nirrt = nirrt / pow(10,6); //MJ cm-2 d-1
     nirrt = nirrt * pow(10,4); // MJ M-2 D-1
     drad = nirrt;
     hdrad[dm][m] = drad;

  // solar radiation to the ground (or soil surface layer)
  // page 124, Coughlan and Running 1997, landscape ecology, beer's law
     surfacerad[dm][m] = drad * exp(- EXT * lai /2.2);
     surfacerad[dm][m] = surfacerad[dm][m] /0.2388 /pow(10,2); // converse
                                                               // to MJ M-2 D-1

  }


// Snow interception and sublimation, Coughlan and Running 1997
double Hydrology::sublimation_canopy(double lai, double canopy_rad,
        const double airt, const double snowfall, const int& dm, const int& m)
{
  double sub_canopy; // mm/day
  double latentv; // latent heat of vaporization (MJmm-1M-2)
  static double ISMAX = 0.5; // mmLAI-1, daily snow interception of snow per
                             // unit leaf area
  static double k = 3.5e5; // MJmm-1M-2 latent heat fusion
  double PSI; // potential snow inteception (mm)
  double SI; // snow inteception (mm)
  double psub; // potential sublimation of snow from canopy (mm)

  latentv = 2.5023e6 - 2430.54 * airt; // J/kg
  latentv = (latentv / 10e6) * 10e9 / 10e6; // MJ mm-1M-2;
  PSI = lai * ISMAX / 2.0;

  if (snowfall < PSI) SI = snowfall;
  else SI = PSI;

  psub = canopy_rad * 1000.0 / (latentv + k);
  if (SI <= psub) sub_canopy = snowfall;
  else sub_canopy = psub;

  snowthroughfall[dm][m] = snowfall - SI;

  return sub_canopy;
}


//  Evaporation from the first layer of soil during wet (Maybe moss or
//  bair ground, depeding on the site
//  Follow Monteith and Unsworth, 1990
double Hydrology::evaporation_layer1(const double& airtemp, const double& pa,
        double& lanmuda, double& vapord, const double& soilrad)
{

// airtemp - air temperature (oc)
// pa - density of air  (kg/m3)
//lamuda - latent heat of vaporization (J/kg)
// vapord - vapor pressure deficit  (Kpa)
// soilr - radiation projected to the soil surface (MJ m-2 day-1)

double evap_soil; // soil evaporation rate mm/day
double eeq; // equlibrium evaporation mm/s
double yipucilong; // the change of latent heat
double gb;  //the boundary layer conductance (which is soil surface) mm/s
double zeta; // as function of serveral parameters
double Tk; // air temperature in degrees Kelvin
double secondpart; // evaporation rate resulted from vapord mm/s
double ttt;
double Gv = 0.462; //  m3 kpa kg-1 K-1, as gas constant for water vapor

  // Jm-2s-1 = W m-2
  ttt = soilrad * pow(10,6) / (24 * 3600); // transfer MJ m-2 day-1 to J m-2 s-1
  gb = 0.01; // here assume 0.01 ms-1; typical value = 10mm/s for soil surface
            //about aerodynamics, see waring et al, 1997, page. 27 	Niger:
            // aero_resis_layer1 = 107 s m-1 (Wallace and Holwill, 1997).
  Tk = airtemp + 273.15 ;  // change airtemp (oc) to kelvin, later on need to
                           //replace Tk with soil surface temperature from STM
  yipucilong = 0.7185 * exp(0.0544 * airtemp);
  zeta = pa * (yipucilong + 1.0) * Gv * Tk;
  eeq = ttt * yipucilong / (pa * lanmuda * (yipucilong + 1.0)) ; // transfer
                                                          //rate m/s to mm/s
  vapord = vapord * pow(10,2)/pow(10,3);  // bar = 10^5pa; mbar = 10^2pa,
                                          //trasfer mbar to kpa

  secondpart = vapord * gb /zeta; //mm/s

  evap_soil = eeq + secondpart; // mm/s total evaporation rate,
  evap_soil = evap_soil * 24 * 3600; // daytime evaporation, mm/day

  return evap_soil;
}


//Sublimation from ground snowpack based on Coughlan PhD thesis, 1991
double Hydrology::snowsublimation_ground( double snowwatereq,
        const double& airtemp, double surfacerad)
{
 static double latent_s = 2845.0; // (KJ/kg) latent heat of sublimation
 static double snowabs = 0.6; // radiation absorptivity of snow;
 double incident_rad; //incident radiation to the snowpack KJ/m2/day
 double sublimation; // actual sublimation, mm/day
 double potentialsub; // potential sublimation of snow

  incident_rad = surfacerad * snowabs * pow(10,3); // convert MJ/m2/day
                                                  // to KJ/m2/day
  if ((airtemp < -1.0) && (snowwatereq > 0.0))
  {
    potentialsub = incident_rad / latent_s;

    if (potentialsub > snowwatereq)
    {
     sublimation = snowwatereq;
    }
    else
     {
     sublimation = potentialsub;
     }
  }
  else sublimation = 0.0;

  return sublimation;
}


//Snow melt from ground driven by incident radiation, Coughlan and Running 1997.
double Hydrology::snowmelt_srad(double snow_rad,
           const double airt, const double snowthrough)

{
 double albedo = 0.80; // dimensionless Running and Coughlan 1988
// static double k = 3.5e5; // MJmm-1M-2 latent heat fusion
 double k = 335.0;  /* (kJ/kg) latent heat of fusion */
 double meltmaxday; // melt rate day mm/day
 double smelt;  // snow melt with mm/day
 double tt;

/*
 if ((airt > -1.0))
  {
   tt = ((snow_rad /0.2388) / 10e6) * 10e4; // covert Cal cm-1 day to MJm-2day-1
   meltmaxday = albedo * (tt / k);

   if (meltmaxday < snowthrough) smelt = meltmaxday;
   else smelt = snowthrough;
  }
 else smelt = 0.0;
*/


 //Following approach is based on Rastetter et al., and Brubaker et al., 1996
 double mq = 2.99; // kg MJ-1 = mm m2 MJ-1
 double ar = 2.0; // mm oC-1 day-1 Brubaker et al. 1996,
  tt = ((snow_rad /0.2388) / 10e6) * 10e4; // covert Cal cm-1 day to MJm-2day-1
  smelt = mq * tt + ar * airt;
  if (smelt < 0.0) smelt =0.0;

 return smelt;
}


// Determine the leaf water potential
double Hydrology::LWP(const double& SOILH2O, const double& SOILCAP,
      const double sand, const double clay, const double silt)
   {
    // LWP daily maximum leaf water potential (Mpa)
    // SOILH2O soil water content at time (t) (mm)
    // SOILCAP soil water capacity, here using the field capacity
    // using SOILCAP = 235 mm before
    double lwp;
    double temp_water;

    temp_water = SOILH2O;
    if (temp_water ==0) temp_water = 0.01; // protected
    lwp = 0.2 /(temp_water / SOILCAP); // Running and Coughlan 1988

/* New solution, see:
	Cosby, B.J., G.M. Hornberger, R.B. Clapp, and T.R. Ginn, 1984.  A
	   statistical exploration of the relationships of soil moisture
	   characteristics to the physical properties of soils.  Water Res.
	   Res. 20:682-690.

	Saxton, K.E., W.J. Rawls, J.S. Romberger, and R.I. Papendick, 1986.
		Estimating generalized soil-water characteristics from texture.
		Soil Sci. Soc. Am. J. 50:1031-1036.

   double soil_b = -(3.10 + 0.157*clay - 0.003*sand);
	double vwc_sat = (50.5 - 0.142*sand - 0.037*clay)/100.0;
  double psi_sat = -(exp((1.54 - 0.0095*sand /100 + 0.0063*silt /100)
                    *log(10.0))*9.8e-5);
	LWP = psi_sat * powl((SOILH2O / vwc_sat), soil_b);
	LWP = psi_sat * powl((SOILH2O /vwc_sat) , soil_b);
*/
// We may also look at Bachelet et al., 1998, sensitivity of a biogeography
//model to soil properties
//about the soil water potential

   return lwp;
  }



// Detemine the vapor pressure deficit and humidity deficit

void Hydrology::vpd(double& vaporP, const double& airt)
{

    double Mw = 18; // g/mol
    double R = 8.3143; // J/mol/k

    double eadt;
    double ea;
    double dewpointT;  //oC

  // vaporP ambient vapor pressure Kpa
  // airt oC
  // absoluteHD g/m3
  // detemine depoint temperature first. see Campbell and Norman,
  //Td = c * ln(ea/ a) /(b - ln(ea/a))
  // where c= 240.97 OC, a = 0.611Kpa, b =17.502

    if (vaporP <= 0.0) vaporP = 0.001;   // for debug purpose, Apr/27/2003

    dewpointT = (240.97 * logl(vaporP/0.611)) / (17.502 - logl(vaporP/0.611));

    eadt = 0.61078 * exp((17.269 * airt) / (airt + 237.30)); // Rosenberg et al.
                                                             // 1983, P. 170
//    vpdeficit =  fabs(eadt - vaporP);
    vpdeficit = eadt - vaporP;
    vpdeficit = (vpdeficit * 10); // trasfer to mbar

    absoluteHD = fabs((Mw / R) * ( (eadt * 1000 / (273.2 +dewpointT))
                 - (vaporP * 1000 /  (273.2 + airt))));

  }


//calculate density of air (rho) as a function of air temperature
double Hydrology::densityA(const double& airt)

{
// pa   (kg/m3)   density of air
// airt air temperature oC
  double pa;
  pa = 1.292 - (0.00428 * airt);
  return pa;
}

// calculate latent heat of vaporization as a function of ta
double Hydrology::latentV(const double& airt)
{
  // le  (J/kg)        latent heat of vaporization of water
  //airt air temperature oC
   double le;
   le = 2.5023e6 - 2430.54 * airt;
   return le;
}

//Simplified equation if Waring and Running, 1998
double Hydrology::CanopyConductance(const double& LAI, const double& gmax,
        const double& airt, const double& vpd, const double& leafwaterp,
        const double& vpd_close,const double& vpd_open, const double& psi_close,
        const double& psi_open)
 {

 // -0.6   open       (MPa) leaf water potential: start of conductance reduction
 //-2.3    close      (MPa) leaf water potential: complete conductance reduction
 //930.0   open       (Pa)  = 0.93Kpa = 9.3 Mbar   vapor pressure deficit:
                                               //start of conductance reduction
 //4100.0  close      (Pa)  = 4.1Kpa  = 41.0 Mbar  vapor pressure deficit:
                      // complete conductance reduction
 // maximum canopy conductance gmax mms-1, Waring ans Running, 1998, P313

 // assume CO2 has no effects on canopy conductance, See Pataki et al., 1998,
                                                     // Oecologia,
  double Gs;
  double ft, fvpd, fsl;

 // air temperature effects
  if (airt > 0.0) ft = 1.0;  // no effect
  else
   {
    if (airt < -8.0) ft = 0.0; // ful effect
    else ft = 1 + 0.125 * airt; // partial effect
   }

 // vapor pressure deficit effects Mbar
   if (vpd < vpd_open) fvpd = 1.0;
   else
   {
    if (vpd > vpd_close) fvpd = 0.0;
    else fvpd = (vpd_close - (-vpd)) / (vpd_close - vpd_open); //partial effects
   }

 // Leaf water potential effect Mpa
  if (leafwaterp > psi_open) fsl = 1.0;
  else
    {
     if (leafwaterp <= psi_close) fsl = 0.0;
     else fsl = (psi_close - leafwaterp) / (psi_close - psi_open);
    }

  Gs = gmax * ft * fvpd * fsl;

  if (Gs ==0.0 ) Gs = 0.001;

  return Gs;
 }

// Implementing Penman-Monteith equation as a whole for the canopy
double Hydrology::penmanmonteith_new(const double& airtemp, double drad,
        const double& rule, const double& vpd, const double& lamuda,
        const double& Gs, const double& LAI)

{
//lamuda = latent heat of vaporization of water (J kg-1)from function of LatentV
 // rule = density of air (kgm-3)  from function of DensityA
 // vpd vapor pressure deficit from canopy to air (mbar)
 double Gv = 0.462;  // gas constant for water vapor (0.462 m3 kpa kg-1 k-1)
 long dayl = 86400; // 24*3600  (second day-1)
 double trans; //  canopy transpiration daily or monthly (m3 day-1)
 double epuci;
 double Omiga;
//canopy aerodynamic resistance(fixed at 5.0 sm-1)typical value 5-10 for forests
 double gb;  // the boundary-layer conductance (determine by wind speed) mm s-1
 double Et;
 double Eeq;
 double Eimp;
 double ttt;
 double gs; // the mean of stomatal conductance to vapor water
 double Tk;
 double lai;
  lai = LAI;
  if (lai < 0.1) {lai = 0.1;}
 // air temperture degree kelvin
  Tk =  airtemp + 273.15; // converse oC to kelvin, Oke, 1995

 // radiation to canopy, LAI to indicating the whole canopy (J m-2 s-1 = W m-2)
  ttt = drad * 10e6 / (24 * 3600); // convert MJ/m2/day to J/m2/s

 // calculate epuci, the change of latent heat relative to the change of sensible
 // heat of saturated air which is 1.26 at 10oC and increase s
 //exponetially with temperature
  epuci = 0.7185 * exp(0.0544 * airtemp);

 // calcualte coefficient of Omiga
  gb = 0.02; // ms-1
  Omiga = (1 + epuci) / (1 + epuci + gb / Gs);

 // calculate evaporation rate (ms-1) in equalibrium (McNaughton, 1976)
  Eeq = ttt * epuci / (rule * lamuda * (1 + epuci));

 // caculate evapotranspiration from dry canopy, vpd with kPa
  gs = Gs / lai;
  Eimp = (vpd /10) * gs / (rule * Gv * Tk);

 // calculate evaporation and transpiration (McNaugton and Jarvis, 1983)
  Et = Omiga * Eeq + (1 - Omiga) * Eimp;
//  trans =  LAI * Et * dayl /2.0; // mm day-1
  trans =  LAI * Et * dayl * 0.4; // mm day-1; light time for transpiration
  //nanprintf("hyd 431: ttt = %.3f, epuci = %.3f, Gs = %.3f, omiga= %.3f, rule = %.3f, lamuda = %.3f, eeq = %.3f, gs = %.3f, vpd = %.3f, eimp = %.3f, et = %.3f, lai = %.3f, dayl = %i, trans = %.3f\n", ttt, epuci, Gs, Omiga, rule, lamuda, Eeq, gs, vpd, Eimp, Et, LAI, dayl, trans);
  return trans;

 }


/* **************************************************************
************************************************************** */
// the following is to read ecd files
void Hydrology::getecd(char ecd[80])
{
  getpenmonecd(ecd);

}

void Hydrology::getecd(ofstream& rflog1)
{

  char ecd[80];

  cout << "Enter name of the file with Penman-Monteith parameter values (.ECD):"
      << endl;
 //cin >> ecd;
  fpara >>ecd;

  rflog1 << "Enter name of the file with Penman-Monteith parameter values (.ECD):";
  rflog1 << ecd << endl << endl;

  getpenmonecd(ecd);

};

void Hydrology::getpenmonecd(char ecd[80])
{

  const int NUMVAR = 10;
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

  infile >> dummy[dcmnt] >> dummy[dcmnt]>> inter_coef[dcmnt] >> EXT[dcmnt]>>
  gmax[dcmnt] >> vpd_close[dcmnt] >> vpd_open[dcmnt] >> psi_close[dcmnt] >>
  psi_open[dcmnt] >> update[dcmnt];
//  printf("hyd dcmnt = %i\n", dcmnt);
   }

  infile.close();

};

/* *************************************************************
************************************************************* */


void Hydrology::getpenmonecd(ofstream& rflog1) {

  const int NUMVAR = 10;
  char dummy[NUMVAR][8];
  ifstream infile;
  int i;
  int dcmnt;
  char ecd[20];

  int update[12];

  infile.open(ecd, ios::in);

  cout << "Enter name of Penman Monteith (.ECD) data file with parametervalues:"
       << endl;
  //cin >> ecd;

rflog1 << "Enter name of Penman Monteith(.ECD) data file with parameter values:"
       << ecd << endl << endl;

  infile.open(ecd, ios::in);

  if (!infile) {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (dcmnt = 1; dcmnt < 35+1; dcmnt++)
  {
    infile >> dummy[0] >> dummy[0]>> inter_coef[dcmnt] >> EXT[dcmnt]
           >> gmax[dcmnt] >> vpd_close[dcmnt] >> vpd_open[dcmnt]
           >> psi_close[dcmnt] >> psi_open[dcmnt] >> update[dcmnt];

  }


  infile.close();

};

/* *************************************************************
************************************************************* */

void Hydrology::showpenmon(int& dcmnt)
{

  cout << endl << "         PARAMETERS FOR THE PENMAN-MONTEITH EQUATION";
  cout << endl << endl;
  printf("     INTER_COEF = %7.5lf     EXT = %7.5lf         gmax = %8.5lf\n",
         inter_coef[dcmnt], EXT[dcmnt], gmax[dcmnt]);
  printf("   vpd_close = %4.2lf       vpd_open = %7.4lf   \n",
         vpd_close[dcmnt], vpd_open[dcmnt]);

  printf("   psi_close = %4.2lf  psi_open = %4.2lf \n",
         psi_close[dcmnt], psi_open[dcmnt] );
};


