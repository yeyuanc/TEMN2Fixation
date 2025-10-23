/* **************************************************************
*****************************************************************
TMCRB423.HPP - object describing characteristics of soil microbes
	       used in the Terrestrial Ecosystem Model (TEM)
*****************************************************************
************************************************************** */

// Tmicrobe4 uses the global constants CYCLE, NUMMSAC and NUMVEG

#if !defined(TMCRB423_H)
#define TMCRB423_H

class Tmicrobe4 {

   public:

/* **************************************************************
		 Public Functions
************************************************************** */

      void   getvegecd(std::ofstream& rflog1);
      void   getvegecd(char ecd[80]);
      double nminxclm(const int& dcmnt, const int& dm, double& soilh2o,
                      double& soilorgc, double& soilorgn, double& availn,
                      double& decay, double& rh, double& ksoil);
      double rhxclm(double& soilorgc, double& dq10, double& moist);
      void   showecd(const int& dcmnt);
      // ez changed to dcmnt in yrkd() by DWK on 20000130
      double yrkd(int nfeed, double& yrltrc, double& yrltrn, int dcmnt);

// for n cycle added by YYL, YY
			void nflux(const int& dcmnt,const double& fn, const double& fdn, const double& rn2, const double& rnox, const double& noxp, const int& dm, const int& m);
															//double& dayn2, double& dayn2o, double& daynox);
			double Nitrif(const int& dcmnt,const double& nmin, const double& t, double pctp, const double& psiplusc, const double& ph, const double& NH4);
			double Denitrif(const int& dcmnt,const double& no3, const double& co2, const double& pctp, const double& ph, const double& t);
			double Rn2(double diff,const double& no3,const double& co2,const double& ph, const double& pctp);
			double N2OUptake(const double& dayn2o, const double& t);
      double NH3Vola(const double& t, const double& NH4);
      double N2Fixation(const int& dcmnt, const double& t, const double& pctp, const double& availn, double& vegcb, const double& density);
      double N2FixationFL (const int& dcmnt,const double& t, const double& pctp, const double& availn, double soilc);
			double Rnox(double diff);
			double noxp(const double& p);
			double Diffusivity(const double& fc, const double& pc, const double& vsm);	
			
/* **************************************************************
		 Public Variables
************************************************************** */

      // monthly heterotrophic respiration (g C / (sq. meter * month))

      double rh[CYCLE];         

      // monthly nitrogen uptake by microbes (g N / (sq. meter * month))

      double nuptake[CYCLE];
      
      // monthly net nitrogen mineralization (g N / (sq. meter * month))

      double netnmin[CYCLE];    

      // annual heterotrophic respiration (g C / (sq. meter * year))

      double yrrh;

      // annual net nitrogen mineralization (g N / (sq. meter * year))

      double yrnmin;            

			// for n cycle YYL
			double N2O[CYCLE];
			double N2[CYCLE];
			double NOx[CYCLE];
//			double n2ot,n2t,noxt;
      
      // for N2O uptake by YYuan
      double N2OF[CYCLE];
			double N2F[CYCLE];
      double N2OAir[CYCLE];
			double N2Air[CYCLE];
      double N2OUpt[CYCLE];
      double N2ON[CYCLE];
      double N2ODN[CYCLE];
      double dayn2air[CYCLE][31];  // YY
      double dayn2oair[CYCLE][31];  // YY
      double dayn2of[CYCLE][31];  // YY
      double dayn2f[CYCLE][31];   // YY 
      double dayn2oupt[CYCLE][31];  // YY
      double NH3[CYCLE];
      double daynh3[CYCLE][31];

			double dayn2odn[CYCLE][31];
			double dayn2on[CYCLE][31];      
			double dayn2o[CYCLE][31];
			double dayn2[CYCLE][31];
			double daynox[CYCLE][31];
			
			double NH4[CYCLE];
			double NO3[CYCLE];
			double Fnit[CYCLE];
			double Fdenit[CYCLE];
			double yrnh4,yrno3;
//			double nh4,no3;
			
      double N2Fix[CYCLE];  //YY N2 fix
      double dayn2fix[CYCLE][31];   //YY N2 fix
      double N2Fixs[CYCLE];  //YY N2 fix
      double dayn2fixs[CYCLE][31];   //YY N2 fix
      double N2Fixfl[CYCLE];  //YY N2 fix
      double dayn2fixfl[CYCLE][31];   //YY N2 fix   
            
			double rn2;
			double rnox;
			double diff;
			
			double noxpi;
			double tp[15], pretp, pulsp;
			int pday, qdry;
			
			
			
/* **************************************************************
		 Public Parameters
************************************************************** */

      // Parameter representing the quality of soil organic matter

      double decay;

      // Biome-specific decomposition parameters for function rhxclm

      double kd;
      double kda[MAXCMNT];
      double kdb[MAXCMNT];
      double kdc;

      double kdin[NUMMSAC];     // kd values read in from file
      double kdsave[NUMMSAC];   // kd values saved to a file

      double lcclnc[MAXCMNT];
      double propftos[MAXCMNT];

      // Biome-specific microbial nitrogen uptake parameters for function
      //   nminxclm

      double nup;
      double nupa[MAXCMNT];
      double nupb[MAXCMNT];

      // Biome-specific N fixation parameter

      double nfixpar[MAXCMNT];

      double cnsoil[MAXCMNT];


      double rhq10[MAXCMNT];
      double k1[MAXCMNT];
      double kmax[MAXCMNT];
      double fno3k[MAXCMNT];
      double fco2k[MAXCMNT];
      double kde[MAXCMNT];
      double kn[MAXCMNT];      
      //Biome-specific parameters describing the influence of soil moisture
      // on decomposition (i.e., moist)
      double npot[MAXCMNT];
      double toptl[MAXCMNT];
      double topth[MAXCMNT];
      double phi2[MAXCMNT];      
      double wuph[MAXCMNT];
      double fnup[MAXCMNT];
      double kc[MAXCMNT];
      double fplant[MAXCMNT];
      double rfl[MAXCMNT];

      double moistmin[MAXCMNT];
      double moistopt[MAXCMNT];
      double moistmax[MAXCMNT];

      // Biome-specific half saturation parameter for function nminxclm
      //   describing the effect of available nitrogen on microbial
      //   nitrogen uptake

      double kn2[MAXCMNT];

};

#endif

