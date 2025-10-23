/* **************************************************************
*****************************************************************
TSOIL423.CPP - object describing general characteristics of soil
            - modified by DWK on 20000102
*****************************************************************
************************************************************** */

#if !defined(TSOIL423_H)
  #include "tsoil423.hpp"
#endif

/* **************************************************************
************************************************************** */

Tsoil4::Tsoil4(void)
{

  text  = -99;
  wsoil = -99;

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tsoil4::getecd(std::ofstream& rflog1)
{

  char ecd[80];

  std::cout << "Enter name of the soil (.ECD) data file with parameter values: ";
  std::cout << std::endl;
  //std::cin >> ecd;
  fpara >> ecd;
  
  rflog1 << "Enter name of the soil (.ECD) data file with parameter values: ";
  rflog1 << ecd << std::endl;

  getecd(ecd);

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void Tsoil4::getecd (char ecd[80])
{

  char dummy[12];
  ifstream infile;

  long update;

  infile.open(ecd, ios::in);

  if (!infile)
  {
    std::cerr << "\nCannot open " << ecd << " for data input" << std::endl;
    exit(-1);
  }

  infile >> dummy >> dummy >> dummy;
  infile >> dummy >> pctpora >> update;
  infile >> dummy >> pctporb >> update;
  infile >> dummy >> fldcapa >> update;
  infile >> dummy >> fldcapb >> update;
  infile >> dummy >> wiltpta >> update;
  infile >> dummy >> wiltptb >> update;

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil4::getrootz(std::ofstream& rflog1)
{

  char ecd[80];

  std::cout << "Enter name of the data file containing the rooting depths:";
  std::cout << std::endl;
  std::cout << "               (e.g., ROOTZVEG.ECD)" << std::endl;
  //std::cin >> ecd;
  fpara >> ecd;

  rflog1 << "Enter name of the data file containing the rooting depths:";
  rflog1 << std::endl;
  rflog1 << "               (e.g., ROOTZVEG.ECD)" << std::endl;
  rflog1 << ecd << std::endl;

  getrootz(ecd);

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil4::getrootz(char ecd[80])
{

  const int NUMVAR = 7;
  char dummy[NUMVAR][10];
  ifstream infile;

  int i;
  int dcmnt;
  int  rootveg[MAXCMNT];
  long update[MAXCMNT];
  char vegname[MAXCMNT][31];

  infile.open(ecd, ios::in);

  if (!infile)
  {
    std::cerr << "\nCannot open " << ecd << " for data input" << std::endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (dcmnt = 1; dcmnt < MAXCMNT; dcmnt++)
  {
    infile >> rootveg[dcmnt] >> vegname[dcmnt];
    infile >> rootza[dcmnt] >> rootzb[dcmnt] >> rootzc[dcmnt];
    infile >> minrootz[dcmnt] >> update[dcmnt];
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil4::lake(double& tair,double& prec,double& rain,double& snowfall,
       		  double& pet, double& eet, int& dm)
{

  rgrndh2o[dm] = 0.0;
  sperc[dm] = 0.0;
  snowpack[dm] = 0.0;
  sgrndh2o[dm] = 0.0;
  moist[dm] = 0.0;

  if (tair >= -1.0)
  {
   rain = prec;
    snowfall = 0.0;
  }
  else
  {
    rain = 0.0;
    snowfall = prec;
  }

  eet = pet;
  h2oyld[dm] = prec - pet;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil4::percol(double& rain, double& snowinf, double& eet,
                    double& avlh2o, const int& dm)
{

  double extra;
  double recharge;
  sperc[dm] = 0.0;
  rperc[dm] = 0.0;

  recharge = rain + snowinf;
  if (recharge <= 0.0) { recharge = 0.001; }
  if ((avlh2o + rain + snowinf - eet) > awcapmm)
  {
    extra = rain + snowinf + avlh2o - awcapmm - eet;
    sperc[dm] = snowinf * extra / recharge;
    rperc [dm] = rain * extra / recharge;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

double Tsoil4::rrunoff(const double& rgrndh2o, const double& rperc)
{

  double rrunof;

  rrunof = 0.5 * (rgrndh2o + rperc);

  return rrunof;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void Tsoil4::showecd(void)
{

  std::cout << std::endl << "                   SOIL CHARACTERISTICS OF SITE";
  std::cout << std::endl << std::endl;
  printf("PSAND    = %5.2lf      PSILT = %5.2lf      PCLAY = %5.2lf\n",
         pctsand, pctsilt, pctclay);

  printf("POROSITY = %5.2lf   PCFLDCAP = %5.2lf   PCWILTPT = %5.2lf\n",
         pctpor, pcfldcap, pcwiltpt);

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

double Tsoil4::snowmelt(const double& elev, const double& tair,
                        const double& prevtair, const double& snowpack)
{

  double snowflux = 0.0;

  if (tair >= -1.0)
  {
    if (elev <= 500.0) { snowflux = snowpack;}
    else
    {
      if (prevtair < -1) { snowflux = 0.5 * snowpack; }
      else { snowflux = snowpack; }
    }
  }

  return snowflux;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

double Tsoil4::srunoff(const double& elev, const double& tair,
                       const double& prevtair, const double& prev2tair,
                       const double& sgrndh2o, const double& sperc)
{

  double srunof = 0.0;

  if (tair >= -1.0)
  {
    if (prevtair < -1.0) { srunof = 0.1 * (sgrndh2o + sperc); }
    else
    {
      if (prev2tair < -1)
      {
	if (elev <= 500.0) { srunof = 0.5 * (sgrndh2o + sperc); }
	else { srunof = 0.25 * (sgrndh2o + sperc); }
      }
      else { srunof = 0.5 * (sgrndh2o + sperc); }
    }
  }

  return srunof;

};

/* *************************************************************
************************************************************* */
// To determine the parameters in DOC modelling, refer to (Neff & Asner 2001) in Ecosystems 2001, 29-48.
void Tsoil4::docmb(const double& soilorgc, const double& density)
{
	double soilc; //soil c content in terms of %
	double soilorgc2;
  soilorgc2=soilorgc;
  //printf("tsoil 314:soilorgc=%.3f\n",soilorgc);
  if (soilorgc2 < 0.001) {soilorgc2=0.001;}
  
  soilc = soilorgc2 / (density * 10000)*10; //transfer soilorgc(g/m2) to soil c content(%), with the assumption that soil c pools are based on 1m depth  added *10 by YYuan 011122 (100 in paper, *10 to match the value in paper)
	//printf("tsoil 311: soilc=%.3f soilorgc=%.3f, density = %.3f\n", soilc, soilorgc, density);
  
  docm = 0.15 * log(soilc) + 0.51;
  if (docm < 0.0) {docm = 0.0;}
	if (docm > 1.0) {docm = 1.0;}
	docb = 0.05 * soilc + 0.09;
	if (docb > 1.0) {docb = 1.0;}	
 	if (docb < 0.0) {docb = 0.0;}	
	docb *= soilc / (soilc + 1.0);
};

/* *************************************************************
************************************************************** */

void Tsoil4::xtext(int& cmnt, double& pctsilt, double& pctclay)
{

  totpor = fldcap = wiltpt = MISSING;
  awcapmm =  MISSING;

  psiplusc = (pctsilt + pctclay) * 0.01;  //printf("tsoil338: psiplusc=%.3f, pctsilt =%.3f, pctclay = %.3f\n", psiplusc ,pctsilt, pctclay);
  if (psiplusc < 0.01) { psiplusc = 0.01; }

  rootz = (rootza[cmnt] * pow(psiplusc, 2.0)) + (rootzb[cmnt] * psiplusc)
          + rootzc[cmnt];
  if (rootz < minrootz[cmnt]) { rootz = minrootz[cmnt]; }

  pctpor = (pctpora * psiplusc) + pctporb;
  pcfldcap = (fldcapa * psiplusc) + fldcapb;
  pcwiltpt = (wiltpta * psiplusc) + wiltptb;

  totpor  = rootz * pctpor * 10.0;
  fldcap  = rootz * pcfldcap * 10.0;
  wiltpt  = rootz * pcwiltpt * 10.0;

  awcapmm = fldcap - wiltpt;

};

/* *************************************************************
************************************************************* */
// added for 3-layer hydrological model, Q. Zhuang, 05/Jan/2003

void Tsoil4::hydm_xtext(int& ez, double& pctsilt, double& pctclay)
 {
 double rootmin;
 double rootmx;

  totpor = fldcap = wiltpt = MISSING;
  awcapmm =  MISSING;

  rootmin = 99.9;
  rootmx = -999.9;

  totpor1 = fldcap1 = wiltpt1 = -999.9;
  totpor2 = fldcap2 = wiltpt2 = -999.9;
  totpor3 = fldcap3 = wiltpt3 = -999.9;

  awcapmm1 =  -999.9;
  awcapmm2 =  -999.9;
  awcapmm3 =  -999.9;

  psiplusc = (pctsilt + pctclay) * 0.01;

  if (psiplusc < 0.01) { psiplusc = 0.01; }

// unit is meter

  rootz = (rootza[ez] * pow(psiplusc, 2.0)) + (rootzb[ez] * psiplusc) + rootzc[ez];

  if (rootz < minrootz[ez]) { rootz = minrootz[ez]; }

  if (rootz< rootmin) {rootmin= rootz;}
  if (rootz> rootmx) { rootmx=rootz;}

  dpwbox1 = rootz * 0.1;  // a thin layer (surface laler 0.1m for the first root zone, maybe the moss, Desborough, 1997
  dpwbox2 = rootz * 0.7; // needed to reconsider it
  dpwbox3 = rootz * 0.2; // mineral account for 20%

  pctpor = (pctpora * psiplusc) + pctporb;
  pcfldcap = (fldcapa * psiplusc) + fldcapb;
  pcwiltpt = (wiltpta * psiplusc) + wiltptb;

  totpor1  = rootz * 0.1 * 0.7 * pctpor * 10.0 ; //dpwbox1  //multiple 0.7 as soil thinkness 1+2+3 = 0.7
  fldcap1  = rootz * 0.1 * 0.7 * pcfldcap * 10.0;
  wiltpt1  = rootz * 0.1 * 0.7 * pcwiltpt * 10.0;

  totpor2  = rootz * 0.7 * 0.7 * pctpor * 10.0 ; //dpwbox1 change 0.8 to 0.7
  fldcap2  = rootz * 0.7 * 0.7 * pcfldcap * 10.0;
  wiltpt2  = rootz * 0.7 * 0.7 * pcwiltpt * 10.0;

  totpor3  = rootz * 0.2 * 0.7 * pctpor * 10.0; //dpwbox1 remove *5
  fldcap3  = rootz * 0.2 * 0.7 * pcfldcap * 10.0;
  wiltpt3  = rootz * 0.2 * 0.7 * pcwiltpt * 10.0;

  totpor  = rootz * pctpor * 10.0;
  fldcap  = rootz * pcfldcap * 10.0;
  wiltpt  = rootz * pcwiltpt * 10.0;

  awcapmm1 = fldcap1 - wiltpt1;
  awcapmm2 = fldcap2 - wiltpt2;
  awcapmm3 = fldcap3 - wiltpt3;

  awcapmm = fldcap - wiltpt;  // estimate capacity for daily hydrological model

};

/* *************************************************************
************************************************************** */

