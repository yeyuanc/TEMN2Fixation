/* **************************************************************
*****************************************************************
TSOIL423.HPP - object describing characteristics of soil used by
	       the Terrestrial Ecosystem Model (TEM)
            - modified by DWK on 20000102
*****************************************************************
************************************************************** */
using namespace std;
#if !defined(BIOMS423_H)
  #include "bioms423.hpp"  // Tsoil4 uses Biomass class
#endif

// Tsoil4 also uses the global constants CYCLE, MAXRTIME, MISSING
//   and NUMVEG

#if !defined(TSOIL423_H)
#define TSOIL423_H

class Tsoil4 {

  public:

     Tsoil4(void);

/* **************************************************************
		 Public Functions
************************************************************** */

     void getecd (ofstream& rflog1);
     void getecd (char ecd[80]);
     void getrootz(ofstream& rflog1);
     void getrootz(char ecd[80]);
     void lake(double& tair, double& prec, double& rain, double& snowfall,
	       double& pet, double& eet, int& dm);
     void percol(double& rain, double& snowinf, double& eet, double& avlh2o, const int& dm);
     double rrunoff(const double& rgrndh2o, const double& rperc);
     void showecd(void);
     double snowmelt(const double& elev, const double& tair, const double& prevtair, const double& snowpack);
     double srunoff(const double& elev, const double& tair, const double& prevtair, const double& prev2tair,
		    const double& sgrndh2o, const double& sperc);
     void xtext(int& ez, double& pctsilt, double& pctclay);
		 void hydm_xtext(int& ez, double& pctsilt, double& pctclay);
		 void docmb(const double& soilorgc, const double& density);
/* **************************************************************
		 Public Variables
************************************************************** */
		 double docm;					// DOC parameter, sorption coefficient
		 double docb;					// DOC parameter, desorption coefficient
			
     int text;              // soil texture (categorical data)
     int wsoil;             // wetland soil type designation
			    // (categorical data)

     double pctsand;        // %sand
     double pctsilt;        // %silt
     double pctclay;        // %clay
     double psiplusc;       // proportion silt and clay
    // 3-layer models

     double dpwbox1 ;        // depth (m) of first-layer soil water box
     double dpwbox2 ;        // depth (m) of first-layer soil water box
     double dpwbox3 ;        // depth (m) of first-layer soil water box

     double awcapmm1 ;        // available water capacity (mm)
	  double awcapmm2 ;        // available water capacity (mm)
	  double awcapmm3 ;        // available water capacity (mm)

     double rtostorage, srtostorage, rrtostorage; // for depositing runoff storage, 07/Jan/2003

     double fldcap1 ;         // volume of water at field capacity (mm)
	  double wiltpt1 ;         // volume of water at wilting point  (mm)
	  double totpor1 ;         // volume of total pore space        (mm)

	  double fldcap2 ;         // volume of water at field capacity (mm)
	  double wiltpt2 ;         // volume of water at wilting point  (mm)
	  double totpor2 ;

	  double fldcap3 ;         // volume of water at field capacity (mm)
	  double wiltpt3 ;         // volume of water at wilting point  (mm)
	  double totpor3 ;

     double awcapmm;        // available water capacity (mm)
     double fldcap;         // volume of water at field capacity (mm)
     double wiltpt;         // volume of water at wilting point  (mm)
     double totpor;         // volume of total pore space        (mm)
		 
		 double doc[CYCLE];
     double avlh2o[CYCLE];  // available water (mm)
     double moist[CYCLE];   // soil moisture (mm)
     double pcfc[CYCLE];    // soil h2o as %field capacity
     double pctp[CYCLE];    // soil h2o as %total pore space
     double vsm[CYCLE];     // soil h2o as %rooting depth

     double rperc[CYCLE];      // rain excess (mm)
     double rgrndh2o[CYCLE];   // rain runoff storage (mm)
     double rrun[CYCLE];       // rain runoff (mm)
     double snowpack[CYCLE];   // snowpack (mm)
     double snowinf[CYCLE];    // snow melt infiltration (mm)
     double sperc[CYCLE];      // snow melt excess (mm)
     double sgrndh2o[CYCLE];   // snowmelt runoff storage (mm)
     double srun[CYCLE];       // snow runoff (mm)
     double h2oyld[CYCLE];     // water yield (mm)
			
		// for hydrology model
		double avlh2o1[CYCLE][31];
		double snowthrough[CYCLE][31];

     double yravlh2o;
     double yrrgrndh2o;
     double yrrperc;
     double yrrrun;
     double yrsnowpack;
     double yrsnowinf;
     double yrsperc;
     double yrsgrndh2o;
     double yrsrun;
     double yrh2oyld;

     
     double yrvsm;
     double yrsmoist;
     double yrpctp;

     double meanvsm;             // mean annual volumetric soil moisture (%)


     Biomass org[CYCLE];         // soil organic matter
     double availn[CYCLE];       // available nitrogen

     // total nitrogen input to soils (g N / (sq. meter * month))
     double ninput[CYCLE];

     // total nitrogen input to soils (g N / (sq. meter * month))
     double tninput[MAXRTIME][CYCLE];
     long ninyear[MAXRTIME];

     // total nitrogen lost from soils (g N / (sq. meter * month))
     double nlost[CYCLE];

     // total nitrogen lost as NH4 from soils (g N / (sq. meter * month))
     double NH4lost[CYCLE];      

     // total nitrogen lost as NO3 from soils (g N / (sq. meter * month))
     double NO3lost[CYCLE];      
            
      // total n gas emission from soils (g N / (sq. meter * month)) N2O N2 NOx
     double ngas[CYCLE];

     int tndepflag;             // flag for transient N deposition data

     double yrorgc;
     double yrorgn;
     double yrc2n;
     double yravln;
     double yrnin;
     double yrnlost;
     double yrngas;


/* **************************************************************
		 Public Parameters
************************************************************** */


     double pctpor;           //porosity of soil (%soil volume)
     double pctpora;
     double pctporb;
     double pcfldcap;         //field capacity (%soil volume)
     double fldcapa;
     double fldcapb;
     double pcwiltpt;         //wilting point (%soil volume)
     double wiltpta;
     double wiltptb;
     double rootz;            //effective rooting depth (m)
     double rootza[MAXCMNT];
     double rootzb[MAXCMNT];
     double rootzc[MAXCMNT];
     double minrootz[MAXCMNT];


// Calibration Parameters

//total nitrogen input into soil (g N / (square meter * month))
     double nintot[MAXCMNT];

//proportion of available nitrogen lost (g N / (square meter * month))
     double nloss[MAXCMNT];

};

#endif
