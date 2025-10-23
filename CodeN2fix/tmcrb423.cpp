/* **************************************************************
*****************************************************************
TMCRB423.CPP - object describing characteristics of soil microbes
               used in the Terrestrial Ecosystem Model (TEM)
*****************************************************************
************************************************************** */

#if !defined(TMCRB423_H)
  #include "tmcrb423.hpp"
#endif

/* **************************************************************
************************* Public Functions **********************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tmicrobe4::getvegecd(std::ofstream& rflog1)
{

  char ecd[80];

  std::cout << "Enter name of the file with microbe parameter values (.ECD)";
  std::cout << std::endl;
  std::cout << "dependent upon vegetation : " << std::endl;
  //std::cin >> ecd;
  fpara >> ecd;
  
  rflog1 << "Enter name of the file with microbe parameter values (.ECD)";
  rflog1 << std::endl;
  rflog1 << "dependent upon vegetation: " << ecd << std::endl << std::endl;

  getvegecd(ecd);

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tmicrobe4::getvegecd(char ecd[80])
{

  const int NUMVAR = 23;
  char dummy[MAXCMNT][31];
  ifstream infile;
  int i;
  int dcmnt;

  int vegid[MAXCMNT];
  char vegname[MAXCMNT][31];
  long update[MAXCMNT];

  infile.open(ecd, ios::in);

  if (!infile)
  {
    std::cerr << "\nCannot open " << ecd << " for data input" << std::endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (dcmnt = 1; dcmnt < MAXCMNT; dcmnt++)
  {
     infile >> vegid[dcmnt] >> vegname[dcmnt];
     infile >> rhq10[dcmnt];
     infile >> kn2[dcmnt];
     infile >> moistmin[dcmnt] >> moistopt[dcmnt] >> moistmax[dcmnt];
     infile >> k1[dcmnt];
     infile >> kmax[dcmnt];
     infile >> fno3k[dcmnt];
     infile >> fco2k[dcmnt];
     infile >> kde[dcmnt];     
     infile >> kn[dcmnt];
     infile >> npot[dcmnt];
     infile >> toptl[dcmnt];
     infile >> topth[dcmnt];
     infile >> phi2[dcmnt];
     infile >> wuph[dcmnt];
     infile >> fnup[dcmnt];
     infile >> kc[dcmnt];
     infile >> fplant[dcmnt];
     infile >> rfl[dcmnt];
     infile >> update[dcmnt];
  }

  infile.close();

};

/* **************************************************************
************************************************************** */



/* **************************************************************
************************************************************** */

double Tmicrobe4::nminxclm(const int& dcmnt, const int& dm, double& soilh2o,
                           double& soilorgc, double& soilorgn, double& availn,
                           double& decay, double& rh, double& ksoil)
{

  double nmin;
  double tcnsoil;
  double immo;

  tcnsoil = cnsoil[dcmnt]; 

  nuptake[dm] = 0.0;
  nmin = 0.0; //printf("tmcrb116:soiln=%.3f\n",soilorgn);
  if (soilorgc > 0.0 && soilorgn > 0.0)
  {
    nuptake[dm]  = (availn * ksoil) / soilh2o;  //printf ("tmcrb 119:nuptake=%.5f, availn=%.3f, ksoil=%.3f，soilh2o=%.3f\n", nuptake[dm], availn, ksoil,soilh2o); 
//   if (nuptake[dm] > 0.00002) nuptake[dm] = 0.000018;  // debug by Q. zhuang, 2/Nov/2001
//   if (nuptake[dm] < 0.00001) nuptake[dm] = 0.000018;  // debug by Q. zhuang, 2/Nov/2001
    nuptake[dm] /= (kn2[dcmnt] + nuptake[dm]); //printf("tmcrb125:nuptake=%.5f, kn2=%.3f\n", nuptake[dm],kn2[dcmnt]);
    
    immo = (nup * nuptake[dm] * decay) * rh;  //printf("tmcrb 124: immo=%.3f, nup=%.3f, nuptake=%.3f, decay=%.3f, rh=%.3f\n", immo, nup,nuptake[dm],decay,rh); 
    nmin = (soilorgn / soilorgc) * rh - immo; 
    //nmin = ((soilorgn / soilorgc) - (nup * nuptake[dm] * decay)) * rh; 
    //printf("tmcrb127: nmin=%.3f,soilorgn=%.3f, soilorgc=%.3f,rh=%.3f\n",nmin,soilorgn,soilorgc,rh);
    if (nmin >= 0.0) { nmin *= (soilorgn/soilorgc) * tcnsoil;}
    else { nmin *= (soilorgc/soilorgn) / tcnsoil; }
    //printf("tmcrb 130: nmin=%.3f,soilorgn=%.3f,soilorgc=%.3f,tcnsoil=%.3f\n",nmin,soilorgn,soilorgc, tcnsoil);
  }
    //if (nmin > 39) {nmin=39;} 
    //printf("tmcrb132: nmin=%.3f\n",nmin); 
  return nmin;

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

double Tmicrobe4::rhxclm(double& soilorgc, double& dq10, double& moist)
{
  //printf("tmcrb147: kd = %.3f, soilc = %.3f, moist = %.3f, dq10 = %.3f\n", kd, soilorgc, moist, dq10);
  return kd * soilorgc * moist * dq10;

//  if (soilorgc > 0.001) { return kd * soilorgc * moist * dq10; }
//  else { return 0.0; }
  
};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tmicrobe4::showecd(const int& dcmnt)
{

  std::cout << std::endl << "   MICROBIAL PARAMETERS INFLUENCED BY CLIMATE";
  std::cout << std::endl << std::endl;
  printf("          KN2 = %6.4lf\n", kn2[dcmnt]);
  printf("        RHQ10 = %6.2lf\n", rhq10[dcmnt]);
  printf("     MOISTMIN = %8.6lf\n", moistmin[dcmnt]);
  printf("     MOISTOPT = %8.6lf\n", moistopt[dcmnt]);
  printf("     MOISTMAX = %8.6lf\n", moistmax[dcmnt]);
  printf("        K1 = %6.4lf\n", k1[dcmnt]);
  printf("        KMAX = %6.4lf\n", kmax[dcmnt]);
  printf("        FNO3K = %6.4lf\n", fno3k[dcmnt]);
  printf("        FCO2K = %6.4lf\n", fco2k[dcmnt]);
  printf("        KDE = %6.4lf\n", kde[dcmnt]);  
  printf("        KN = %6.4lf\n", kn[dcmnt]);
  printf("        NPOT = %6.4lf\n", npot[dcmnt]);
  printf("        TOPTL = %6.4lf\n", toptl[dcmnt]);
  printf("        TOPTH = %6.4lf\n", topth[dcmnt]);
  printf("        PHI2 = %6.4lf\n", phi2[dcmnt]);  
  printf("        WUPH = %6.4lf\n", wuph[dcmnt]);
  printf("        FNUP = %6.4lf\n", fnup[dcmnt]);
  printf("        KC = %6.4lf\n", kc[dcmnt]);
  printf("        FPLANT = %6.4lf\n", fplant[dcmnt]);     
  printf("        RFL = %6.4lf\n", rfl[dcmnt]);          
  printf("       CNSOIL = %5.2lf\n", cnsoil[dcmnt]);
};

/* **************************************************************
************************************************************** */

/* **************************************************************
************************************************************** */

double Tmicrobe4::yrkd(int nfeed, double& yrltrc, double& yrltrn, int dcmnt)
{

  double yrkd;

  if (yrltrn < 0.001 || yrltrc < 0.001) { yrkd = kdc;}
  else {
  /*{
    std::cout << "YRLTRC is < 0.0 in microbe.yrkd()" << std::endl;
    exit(-1);
  } */
  if (nfeed == 0) { yrkd = kdc; }
  else
  {
    yrkd = kdc * pow((yrltrc/yrltrn),-0.784) / pow(lcclnc[dcmnt],-0.784);
  } }
  //printf("tmcb 270: dcmnt=%i, yrkd = %.3f, kdc = %.3f, yrltrc = %.3f, yrltrn = %.3f, lcclnc = %.3f\n", dcmnt, yrkd, kdc, yrltrc, yrltrn, lcclnc[dcmnt]);
  return yrkd;
};

/* **************************************************************
************************************************************** */

/* **************************************************************
************************************************************** */
//added by Yanyu Lu 06/03/2009, for N cycle including Nitrification and Denitrification
/* **************************************************************
************************************************************** */
void Tmicrobe4::nflux(const int& dcmnt,const double& fn, const double& fdn, const double& rn2, const double& rnox, const double& noxp, const int& dm, const int& m)
															///double& dayn2, double& dayn2o, double& daynox)
{
	//double n2o,n2,nox; //N gas flux
	double n2o_n, n2o_dn;
	
	//const double kn=0.02; // 0.02 the fraction of nitrificated N lost as N2O flux, number from Parton et al (2001) 0.0006 DNDC Li et al 2000
	n2o_n=kn[dcmnt]*fn;
  //printf("tmcrb 230: kn=%.5f\n", kn[dcmnt]);
	const double kdn=1/(1+rn2); // the fraction of denitrificated N lost as N2O flux modified by YYuan (110121), constant 0.2 in previous version
	n2o_dn=kdn*fdn;
	
  dayn2on[dm][m] = n2o_n; 
  dayn2odn[dm][m] = n2o_dn;
  dayn2o[dm][m] = n2o_n+n2o_dn;
	dayn2[dm][m]  = n2o_dn*rn2;
	daynox[dm][m] = rnox*n2o_dn + rnox*n2o_n*noxp ; // the precipitation events? Parton et al (2001)
	//printf("tmcrb 242: dm = %i, m = %i, dayn2o = %.3f, dayn2 = %.3f, daynox = %.3f\n", dm, m, dayn2o[dm][m], dayn2[dm][m], daynox[dm][m]);
  //printf("tmcrb 243: fn = %.3f, n2o_n = %.3f, rn2 = %.3f, kdn = %.3f, fdn = %.3f, n2o_dn = %.3f\n", fn, n2o_n, rn2, kdn, fdn, n2o_dn);

	//return 1;

}



/* **************************************************************
************************************************************** */
// Nitrification
double Tmicrobe4::Nitrif(const int& dcmnt, const double& nmin, const double& t, double pctp, const double& psiplusc, const double& ph, const double& NH4)
{
	double fn;
	//const double K1 = 0.05; //Different in different sites (0.2 original)
	//const double Kmax = 0.01; //Different in different sites (0.1 original)
	double ft, fm, fph; // the effect of soilt, soil moisture, pH on nitrification
	double wfps; 
  
 	if (t < 0){ft = 0.185;}
  //else if (t < 10){ft = pow((60-t))/25.78,3.503) * exp(3.502*(t-34.22)/25.78) + 0.2;} 
	else{ft = pow((60-t)/25.78,3.503) * exp(3.502*(t-34.22)/25.78);}	 		// the effect of soil t, Li et al 2000 and Huang & Gerber 2015 optimum at 35

	double a,b,c,d;
	wfps = pctp / 100; 

	if (psiplusc < 50){a=0.55; b=1.70; c=-0.007; d=3.22;} // psiplusc is the ratio of silt and clay to determine the type of soil texture
	else{a=0.60; b=1.27; c=0.0012; d=2.84;}
	fm = pow((wfps-b)/(a-b), d*(b-a)/(a-c))	* pow((wfps-c)/(a-c),d); // the effect of soil moisture, Parton et al (1996)
	//printf("tmcrb 233: psiplusc = %.3f, wfps = %.3f, pow11 = %.3f, pow12 = %.3f, pow21 = %.3f, pow22 = %.3f\n", psiplusc, wfps, (wfps-b)/(a-b), d*(b-a)/(a-c), (wfps-c)/(a-c), d);
	fph = 0.56 + atan(3.1415926*0.45*(-5+ph))/3.1415926; // the effect of soil pH, Parton et al (1996)
	
	if (wfps > 0.80) {fn = 0;} 		// when soil wfps > 80, the nitrification will be inhibited. Mosier et al (2002), Henault et al (2005)
	else{ 
		 fn = k1[dcmnt] * nmin + kmax[dcmnt] * NH4 * ft * fm * fph;
	}
	//printf("trmcb 279: fn = %.3f, nmin = %.3f, nh4 = %.3f, ft = %.3f, fm = %.3f, fph = %.3f\n",fn, nmin, NH4, ft, fm, fph);
  //printf("tmcrb 280:k1=%.3f,kmax=%.3f\n",k1[dcmnt],kmax[dcmnt]);
	return fn;
}	

//double Tmicrobe4::N_EffectT(const double& t)
//{
//	double n_ft;
//	
//	if (t < 11){n_ft = exp(((t-11)*log(89)-9*log(2.1))/10);}
//	else{n_ft = exp((t-20)*log(2.1)/10);}	
//	
//	return n_ft;	
//}
//
//double Tmicrobe4::N_EffectM(const double& wfps, const double& psiplusc)
//{
//	double n_fm;
//	double a,b,c,d;
//	
//	if (psiplusc < 50){a=0.55; b=1.70; c=-0.007; d=3.22;}
//	else{a=0.60; b=1.27; c=0.0012; d=2.84;}
//		
//	n_fm = pow((wfps-b)/(a-b), d*(b-a)/(a-c))	* pow((wfps-c)/(a-c),d);
//	return n_fm;
//}	
//
//double Tmicrobe4::N_EffectpH(const double& ph)
//{
//	double n_ph;
//	
//	n_ph = 0.56 + atan(3.1415926*0.45*(-5+ph))/3.1415926;
//	
//	return n_ph;
//
//}
/* **************************************************************
************************************************************** */

/* **************************************************************
************************************************************** */
// Denitrification
//double Tmicrobe4::Denitrif(double dn_fno3, double dn_fco2, double dn_fm, double dn_ft)
double Tmicrobe4::Denitrif(const int& dcmnt,const double& no3, const double& co2, const double& pctp, const double& ph, const double& t)	
{
	double wfps;
  double fdn;   
	double fno3, fco2, fm, ft; // the effect of NO3, C, moisture and temperature on denitrification
	double fphdn;  //the effect of pH on denitrification, added by YYuan 110121
	double no3_2, co2_2;
 
	no3_2 = no3;
  co2_2 =co2;
	if (no3_2 < 0.001) {no3_2 = 0.001;}
  if (co2_2 < 0.001) {co2_2 = 0.001;}
	
  fno3 = fno3k[dcmnt] * pow(no3_2, 0.57); // the unit of no3 is ug N g-1 d-1, Del Grosso(2000) (1.15)
	fco2 = fco2k[dcmnt] * pow(co2_2, 1.3); // the unit of co2 is ug C g-1 d-1, Del Grosso(2000)	(0.1)
	
	wfps = pctp / 100;
  if (dcmnt==9){wfps=wfps*5;}
  //printf("tmcrb 334: wfps=%.3f, pctp=%.3f\n", wfps,pctp);
	if (wfps < 0.4){fm = 0;} //previous 0.62
	else {fm = pow((wfps-0.4)/0.38, 1.74);} //NOE algorithm, Henault(2005)
 
 	if (t < 1){ft = 0.1;}
	else if (t < 11){ft = exp(((t-11)*log(89)-9*log(2.1))/10);}          //original t-11 
	else{ft = exp((t-20)*log(2.1)/10);}	 									//NOE, Henault(2005) 	
	
	if (ph<=3.5){fphdn=0.3;}   //add fphdn by YYuan(110121), SWAT(Moges et al., 2017) 
  else if (ph<=6){fphdn=(ph-3)/3;}
	else {fphdn=1;}

	if (fno3 > fco2){ fdn = kde[dcmnt] * fco2 * fm * ft * fphdn ;}		// the Liebig law of the minimum with either available C or NO3 limiting the denitrification
	else { fdn =  kde[dcmnt] * fno3 * fm * ft *fphdn;} //*fphdn
	//printf("tmcrb 353: fdn = %.3f, kde=%.3f, fco2 = %.3f, fno3 = %.3f, fm = %.3f, ft = %.3f, fphdn = %.3f\n", fdn, kde[dcmnt], fco2, fno3, fm, ft, fphdn);
	//printf("tmcrb312: fdn=%.6f\n",fdn);	
	return fdn;	
}

//double Tmicrobe4::DN_EffectNO3(const double& no3)
//{
//	double dn_fno3;
//	//dn_fno3 = 11000 + 40000 * atan(3.1415926*0.002*(no3-180))/3.1415926; //the unit of no3 is ug N g-1
//	dn_fno3 = 1.15 * pow(no3, 0.57); // the unit of no3 is ug N g-1 d-1, Del Grosso(2000)
//	return dn_fno3;
//}
//
//double Tmicrobe4::DN_EffectC(const double& co2)
//{
//	double dn_co2;
//	//dn_co2 = 24000/(1+200/exp(0.35*co2))-100; //co2 is the soil heterotrophic respiration rate  kg C ha-1 d-1
//	dn_co2 = 0.1 * pow(co2, 1.3); // the unit of co2 is ug C g-1 d-1, Del Grosso(2000)
//	return dn_co2;
//}
//
//double Tmicrobe4::DN_EffectM(const double& wfps)
//{
//	double dn_fm;
//	
//	if (wfps < 0.62){dn_fm = 0;}
//	else {dn_fm = pow((wfps-0.62)/0.38, 1.74);} //NOE algorithm, Henault(2005)
//		
//	return dn_fm;
//}
//
//double Tmicrobe4::DN_EffectT(const double& t)
//{
//	double dn_ft;
//	
//	if (t < 11){dn_ft = exp(((t-11)*log(89)-9*log(2.1))/10);}
//	else{dn_ft = exp((t-20)*log(2.1)/10);}	 									//NOE, Henault(2005) 
//	
//	return dn_ft;	
//}	

/* **************************************************************
************************************************************** */

/* **************************************************************
************************************************************** */

double Tmicrobe4::Rn2(double diff,const double& no3,const double& co2,const double& ph, const double& pctp)
{	
	double rn2; // the ratio of N2/N2O
	double k1;
	double fnc; // the effect of nitrate level and soil respiration
	double fm; // the effect of WFPS
	double fphrn2; //the effect of pH on N2/N2O ratio 
	double wfps;
  double co2_rn2, no3_rn2;
	wfps = pctp; //(%)
  co2_rn2 = co2;
  no3_rn2 = no3;
  
  if (ph<=4.5){fphrn2=0.1;}
  else {fphrn2 =1/(1470*exp(-1.1*ph));}  //added by YYuan (110121), SWAT (Moges et al.,2017)

	if ((38.4-350*diff)<1.7){k1=1.7;}
	else {k1=38.4-350*diff;}
	
  if (co2_rn2 < 0.001) {co2_rn2 = 0.001;}
  //if (no3_rn2 < 0.001) {no3_rn2 = 0.001;}
  
	if (0.16*k1 > k1*exp(-0.8*(no3_rn2/co2_rn2))){fnc=0.16*k1;}
	else{fnc = k1*exp(-0.8*(no3_rn2/co2_rn2));}
	//printf("tmcrb 415: fnc=%.3f,no3=%.3f,co2=%.3f\n",fnc,no3_rn2,co2_rn2);	
	if ((0.015*wfps-0.32) < 0.1){fm=0.1;}
	else{fm=0.015*wfps-0.32;}

	rn2=fnc*fm*fphrn2;   
	//printf("tmcrb 426: rn2 = %.3f, fnc = %.3f, fm = %.3f, fphrn2 = %.3f\n", rn2, fnc, fm, fphrn2);
  if (rn2 < 0.3) {rn2 = 0.3;}
  if (rn2 > 99 ) {rn2 = 99;}  // n2:n2o ranges from 0.3 to 261 from Butterbach-Bahl et al., 2013 and 0-33 from Schlesinger 2009 PNAS
  //printf("tmcrb 429: rn2=%.3f",rn2);
	return rn2;
	
}
/* **************************************************************
************************************************************** */

///// Atmospheric N2O uptake

double Tmicrobe4::N2OUptake(const double& dayn2o,const double& daysnowpack)
{
  double n2oupt;
	double dayn2ocon;
  const double dayn2oconair = 0.596 /1000; //atmospheric n2o concentration 331.1 ppb (0.3311 ppm *44/22.4 * 274/(273+25)mg , convert to g N2O m-3

	dayn2ocon = dayn2o*5.24; //soil n2o concentration in 30 cm soil (g N2O m-3 day-1) dayn2o/0.3/28*44
  //printf("dayn2ocon=%.5f\n",dayn2ocon);
  
  if (dayn2oconair <= dayn2ocon){n2oupt = 0;} 
  else if (daysnowpack>5){n2oupt = 0.1*(dayn2oconair-dayn2ocon)*0.1089*28/44;} 
  else {n2oupt=(dayn2oconair-dayn2ocon)*0.1089*28/44;} //n2o uptake in (mg N m-3 day-1); 30cm or 1m? (van Bochove, E., Bertrand, N., & Caron, J. (1998) 1.26*10-6 m-2 s-1)
  //printf("tmcrb399: daysnowpack=%.2f\n",daysnowpack);
  //printf("tmcrb 410: n2oupt=%.6f\n",n2oupt);
  return n2oupt;
}

/* **************************************************************
************************************************************** */
double Tmicrobe4::NH3Vola(const double& t,const double& NH4)
{
  double NH3V;
  double daynh3con;
  double hm;
  double tk;
  const double daynh3conair = 0.0021/1000;  // atmospheric nh3 concentration 1-5 ppb (0.003 ppm *17/22.4*274/(273+25)) mg/m-3, convert to g N2O m-3
  tk=t+273.15;
  
  if (tk>=273.15){
  daynh3con= NH4*(0.2138/tk)*pow(10,6.123-1826/tk);
  hm=0.000612*pow(t,0.382);}
  
  if (daynh3con <=daynh3conair) NH3V =0;
  else NH3V = hm * 3600 * 12 * (daynh3con-daynh3conair);

  //printf("tmcrb424:t=%.3f, nh4=%.3f, daynh3con=%.9f,hm=%.9f,NH3=%.9f\n",t, NH4, daynh3con,hm,NH3V);
  return NH3V;
}
  
     
/* **************************************************************
************************************************************** */

///// Biological nitrogen fixation

double Tmicrobe4::N2Fixation (const int& dcmnt, const double& t, const double& pctp, const double& availn, double& vegcb, const double& density)
{
  double n2fixs;
	double ft, fm;
  double fncon;
  double fcavail;

  
  if (t >= 40 || t <= 0) {ft = 0;} 
  else if (t <= toptl[dcmnt]) {ft = t/toptl[dcmnt];}
  else if (t <= topth[dcmnt]) {ft = 1;}
  else {ft = (40-t)/(40-topth[dcmnt]);}
  //printf("tmcrb484:t=%.3f, toptl=%.3f, topth=%.3f, ft=%.3f\n",t,toptl[dcmnt],topth[dcmnt],ft);
  
  if (pctp <= 0) {fm = 0;}
  else if (pctp <= wuph[dcmnt]) {fm=phi2[dcmnt]*pctp;}
  else {fm=1;}
  if (fm>1) {fm=1;}
  //printf("tmcrb490: pctp=%.3f, wuph=%.3f, ohi2=%.3f, fm=%.3f\n",pctp,wuph[dcmnt],phi2[dcmnt],fm);
  
  if (availn >= 0.001 ) {fncon = 1 - fnup[dcmnt] * log (1000*availn);}
  else {fncon =1;}
  //printf("tmcrb494: availn=%.3f, fnup=%.3f, fncon=%.3f\n", availn,fnup[dcmnt],fncon);
  
  fcavail = 1/(1+kc[dcmnt]*density*1000000/vegcb);
  //printf ("tmcrb497: kc=%.3f, soilorgc=%.3f,fcavail=%.3f,density=%.3f\n", kc[dcmnt],vegcb,fcavail,density);
  
  n2fixs = npot[dcmnt] * ft * fm * fncon * fcavail * fplant[dcmnt];
  
  //printf ("tmcrb500: n2fixs=%.3f, npot=%.3f,ft=%.3f,fm=%.3f, fncon=%.3f,fcavail=%.3f,fplant=%.3f\n", n2fixs,npot[dcmnt], ft,fm,fncon, fcavail,fplant[dcmnt]);
  if (n2fixs<0) {n2fixs=0;}
  return n2fixs;
  
}


double Tmicrobe4::N2FixationFL (const int& dcmnt,const double& t, const double& pctp, const double& availn, double soilc)
{
  double n2fixfl;
  double ft,fm, fn;
  double wfps;
  wfps = pctp/100;
  
  ft = exp(-2.6+0.21*t*(1-0.5*t/24.4)); //printf("tmcrb514: t=%.4f, ft=%.4f\n",t,ft); 
  
  if (wfps >= 0.67) {fm=1;}
  else {fm=5.0682 * pow(wfps, 2) - 2.4045 * wfps + 0.3495;} 
  //printf("tmcrb529: wfps=%.3f,pctp=%.3f,fm=%.3f\n",wfps,pctp,fm);
  
  if (availn >= 0.055 ) {fn = 0.413 * exp (0.0487/availn);}
  else {fn =1;} //printf("tmcrb517: availn=%.4f, fn=%.6f\n",availn,fn);
  
  n2fixfl= rfl[dcmnt]*ft*fn*soilc*1.5/(1500+soilc);  //*1500soilc/(1500+soilc)
  //printf("tmcrb519: rfl=%.4f, soilc=%.3f, n2fixfl=%.6f\n",rfl[dcmnt],soilc,n2fixfl);
  
  if (n2fixfl<0) {n2fixfl=0;}
  return n2fixfl; 

}

/* **************************************************************
************************************************************** */
double Tmicrobe4::Rnox(double diff)
{	
	double rnox; // the ratio of NOx/N2O
	
	rnox = 15.2 + (35.5* atan(0.68*3.1415926*(10*diff-1.86)))/3.1415926;
	
	return rnox;
	
	
}	

double Tmicrobe4::noxp(const double& p)
{

	int i;
	for (i=1;i<15;i++){tp[i]=tp[i-1];pretp+=tp[i];}
	tp[0]= p;
			
	if (p>0)pday++;
	//	else pday = 0;
			
	if (pday==1){
		if (pretp <=10)qdry=1;
			else qdry=0;
	}
	
	if (qdry==1){
		for(i=0;i<pday;i++) pulsp+=tp[i];
		if (pday <= 3){ 
			if(pulsp >=1 && pulsp < 5) noxpi=11.10*exp(-0.805*pday);
		}
		if (pday <= 7){
			if(pulsp >=5 && pulsp < 15) noxpi=14.68*exp(-0.384*pday); 
		}
		if (pday <= 14){
			if(pulsp >=15) noxpi=18.46*exp(-0.208*pday); 
		}	
	}
	if (p == 0) {qdry = 0; pday = 0;}
	if (qdry != 1) {noxpi = 1.0;}	
	
	return noxpi;
}		
/* **************************************************************
************************************************************** */

/* **************************************************************
************************************************************** */
double Tmicrobe4::Diffusivity(const double& fc, const double& pc, const double& vsm)		
{
//	double diff; // Diffusivity index refer to Potter et al(1996)
	double Wa,Wp,A,P,S,Swa,Swp;
	double m,x,y,z;
	
	A=fc;
	P=pc-A;
	S=1-pc;
	
	m=vsm;
	if(m<A){Wa=m;Wp=0;}
	else{Wa=0.999*A;Wp=m-A;} // for month step, the Wp should be divided by the number of days per month
	if(Wp > P) Wp = P ;
	
	Swa=Wa/A;
	Swp=Wp/P;
	
	x=0.477*pow(P,3)-0.596*pow(P,2)+0.437*P+0.564;
	y=0.477*pow((P-Wp),3)-0.596*pow((P-Wp),2)+0.437*(P-Wp)+0.564;
	z=0.477*pow(((A-Wa)/(A+S)),3)-0.596*pow(((A-Wa)/(A+S)),2)+0.437*((A-Wa)/(A+S))+0.564;
	
	diff=pow((1-Swa),2)*pow(((A-Wa)/(A+S)),2*z)*(1-pow(P,2*x))*(P-Wp-pow(P-Wp,2*y))/
			(pow((1-Swa),2)*pow(((A-Wa)/(A+S)),2)*(1-pow(P,2*x))+(P-Wp)-pow(P-Wp,2*y)) 
			+pow((1-Swp),2)*pow((P-Wp),2*y);
  if (diff < 0.01){diff=0.01;}
	return diff;	
	//return 0.5;
}

