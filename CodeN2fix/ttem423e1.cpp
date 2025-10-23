 /* *************************************************************
****************************************************************
TTEM423.CPP - Terrestrial Ecosystem Model Version 4.2
****************************************************************

Modifications:

19991028 - DWK add bug fixes
20000107 - DWK adds compiler directive
20000107 - DWK renames functions to better describe purpose
20000130 - DWK changes prevleafmx to initleafmx in fourth getsiteecd()
20000130 - DWK changes y[2] to y[I_SOLC] and y[3] to y[I_SOLN]
           at bottom of stepyr(
20000201 - DWK changes pstate[1] to pstate[I_STRN] in delta()
20000207 - DWK added a minimum VEGC and SOLC (0.00001 gC m-2) condition
           to massbal(); Commented out by DWK on 20000210
20000614 - DWK veg.ingpp changed to veg.ingpp[] in delta(), setELMNTflux(), and
           setmonth()
20000614 - DWK veg.inpp changed to veg.innpp[] in delta(), setELMNTflux(), and
           setmonth()
20000614 - DWK resorts function names alphabetically to allow function
           descriptions to be found easier
20000616 - CFD removed Rh and Nmin from boundcon(), added RvMaint to boundcon()
20000616 - CFD restructured massbal()
20000630 - DWK changes missing values from -99.99 to -999.99 in ecdqc()
20010418 - Q. Z. Soil thermal model
20020202 - Kick and Qianlai, add leaf to calculation of LAI in delta()
20020202 - Q. Z change the include from ttem423e.cpp to ttem423e1.cpp

****************************************************************
************************************************************** */

#ifndef TTEM423E1_H
  #include "ttem423e1.hpp"
#endif

/* *********************************************************** */

TTEM::TTEM() : Odeint4()
{

  nfert = 1.00000;
  tol = inittol;
  totyr = -99;

// Identify potential output variables from TEM

// Ecosystem carbon pools ***************************************

  strcpy(predstr[I_VEGC],"VEGC");       // vegetation carbon
  strcpy(predstr[I_SOLC],"SOILORGC");   // soil organic carbon
  strcpy(predstr[I_TOTC],"TOTALC");     // total carbon
  strcpy(predstr[I_VEGCAG],"VEGCAG");     // aboveground vegetation carbon
  strcpy(predstr[I_VEGCBG],"VEGCBG");     // belowground vegetation carbon

// Ecosystem nitrogen pools *************************************

  // total nitrogen stored in vegetation
  strcpy(predstr[I_VEGN],"VEGN");

  // vegetation structural nitrogen
  strcpy(predstr[I_STRN],"VSTRUCTN");

  // vegetation labile nitrogen
  strcpy(predstr[I_STON],"VSTOREN");

  strcpy(predstr[I_SOLN],"SOILORGN");  // soil organic nitrogen
  strcpy(predstr[I_AVLN],"AVAILN");    // soil available nitrogen


// Carbon and nitrogen pools associated with human products *****

  // carbon in agricultural products
  strcpy(predstr[I_AGPRDC],"AGPRODC");

  // nitrogen in agricultural products
  strcpy(predstr[I_AGPRDN],"AGPRODN");

  // carbon pool of products that decompose in 10 years
  strcpy(predstr[I_PROD10C],"PROD10C");

  // nitrogen pool of products that decompose in 10 years
  strcpy(predstr[I_PROD10N],"PROD10N");

  // carbon pool of products that decompose in 100 years
  strcpy(predstr[I_PROD100C],"PROD100C");

  // nitrogen pool of products that decompose in 100 years
  strcpy(predstr[I_PROD100N],"PROD100N");

  // carbon in all product pools
  strcpy(predstr[I_TOTPRDC],"TOTPRODC");

  // nitrogen in all product pools
  strcpy(predstr[I_TOTPRDN],"TOTPRODN");

  // total carbon pool found in ecosystem excluding products
  strcpy(predstr[I_TOTEC],"TOTEC");

  // total carbon pool found in ecosystem including products
  strcpy(predstr[I_TOTGC],"TOTGC");


// Ecosystem water pools ****************************************

  // available soil moisture
  strcpy(predstr[I_AVLW],"AVAILH2O");

  // groundwater pool resulting from rainfall
  strcpy(predstr[I_RGRW],"RGRNDH2O");

  strcpy(predstr[I_SNWPCK],"SNOWPACK");  // snowpack

  // groundwater pool resulting from snow melt
  strcpy(predstr[I_SGRW],"SGRNDH2O");

  strcpy(predstr[I_SM],"SOILH2O");       // soil moisture

   // added on 15/NOV/98 for soil temperature
   strcpy(predstr[I_TSOIL],"TSOIL");
   strcpy(predstr[I_DST5],"DST5");
   strcpy(predstr[I_DST10],"DST10");
   strcpy(predstr[I_DST20],"DST20");
   strcpy(predstr[I_DST50],"DST50");
   strcpy(predstr[I_DST100],"DST100");
   strcpy(predstr[I_DST200],"DST200");
   strcpy(predstr[I_FRONTD],"FRONTD");
   strcpy(predstr[I_THAWBE],"THAWBE");
   strcpy(predstr[I_THAWEND],"THAWEND");
// end of ...



  // soil moisture expressed as percent total porosity
  strcpy(predstr[I_PCTP],"PCTP");

  strcpy(predstr[I_VSM],"VSM");      // volumetric soil moisture


// Carbon fluxes for natural ecosystems *************************

  // GPP not limited by nutrient availability
  strcpy(predstr[I_INGPP],"VEGINGPP");

  strcpy(predstr[I_GPP],"GPP");      // gross primary production

  // NPP not limited by nutrient availability
  strcpy(predstr[I_INNPP],"VEGINNPP");

  strcpy(predstr[I_NPP],"NPP");      // net primary production
  strcpy(predstr[I_GPR],"GPR");      // gross plant respiration

  // vegetation maintenance respiration
  strcpy(predstr[I_RVMNT],"RVMAINT");

  // vegetation growth respiration
  strcpy(predstr[I_RVGRW],"RVGRWTH");

  strcpy(predstr[I_LTRC],"LTRC");   // litterfall carbon
  strcpy(predstr[I_RH],"RH");       // heterotrophic respiration
  strcpy(predstr[I_NEP],"NEP");     // net ecosystem production


// Nitrogen fluxes for natural ecosystems ***********************

  // total nitrogen inputs into ecosystem
  strcpy(predstr[I_NINP],"NINPUT");

  // VEGNUP not limited by carbon availability
  strcpy(predstr[I_INNUP],"VEGINNUP");

  // nitrogen uptake by vegetation
  strcpy(predstr[I_VNUP],"VEGNUP");

  // vegetation nitrogen uptake for structural components
  strcpy(predstr[I_VSUP],"VEGSUP");

  // vegetation nitrogen uptake for labile components
  strcpy(predstr[I_VLUP],"VEGLUP");

  // nitrogen mobilization by vegetation
  strcpy(predstr[I_VNMBL],"VNMOBIL");

  strcpy(predstr[I_VNRSRB],"VNRESORB"); // nitrogen resorption by vegetation
  strcpy(predstr[I_LTRN],"LTRN");       // litterfall nitrogen
  strcpy(predstr[I_MNUP],"MICRONUP");   // nitrogen uptake by microbes
  strcpy(predstr[I_NMIN],"NETNMIN");    // net nitrogen mineralization
  strcpy(predstr[I_NLST],"NLOST");      // total nitrogen losses from ecosystem


// Water fluxes *************************************************

  strcpy(predstr[I_RAIN],"RAIN");        // rainfall

  // percolation of rainwater through soil profile
  strcpy(predstr[I_RPERC],"RPERC");

  strcpy(predstr[I_RRUN],"RRUN");        // runoff of rainwater
  strcpy(predstr[I_SNWFAL],"SNOWFALL");  // snowfall

  // infiltration into the soil of water from snowmelt
  strcpy(predstr[I_SNWINF],"SNOWINF");

  // percolation of snowmelt through soil profile
  strcpy(predstr[I_SPERC],"SPERC");

  strcpy(predstr[I_SRUN],"SRUN");        // runoff of snowmelt
  strcpy(predstr[I_PET],"PET");          // potential evapotranspiration
  strcpy(predstr[I_EET],"EET");          // estimated evapotranspiration
  strcpy(predstr[I_WYLD],"H2OYIELD");    // water yield


// Phenology variables for natural ecosystems *******************

  // un-normalized relative phenology variable
  strcpy(predstr[I_UNRMLF],"UNRMLEAF");

  // normalized relative phenology variable (0 - 1.0)
  strcpy(predstr[I_LEAF],"LEAF");

  strcpy(predstr[I_FPC],"FPC");          // foliar projected cover
  strcpy(predstr[I_LAI],"LAI");          // leaf area index


// Carbon and nitrogen fluxes associated with agricultural
// conversion ***************************************************

  strcpy(predstr[I_CNVRTC],"CONVERTC"); // carbon loss during conversion
  strcpy(predstr[I_CNVRTN],"CONVERTN"); // nitrogen loss during conversion
  strcpy(predstr[I_SCNVRTC],"SCONVRTC");
  strcpy(predstr[I_SCNVRTN],"SCONVRTN");
  strcpy(predstr[I_NVRTNT],"NVRETENT");
  strcpy(predstr[I_NSRTNT],"NSRETENT");
  strcpy(predstr[I_NRETNT],"NRETENT");

  // carbon associated with slash left after conversion
  strcpy(predstr[I_SLASHC],"SLASHC");

  // nitrogen associated with slash left after conversion
  strcpy(predstr[I_SLASHN],"SLASHN");

  // carbon loss to formation of products that decompose in 10 years
  strcpy(predstr[I_PRDF10C],"PRDF10C");

  // nitrogen loss to formation of products that decompose in 10 years
  strcpy(predstr[I_PRDF10N],"PRDF10N");

  // carbon loss to formation of products that decompose in 100 years
  strcpy(predstr[I_PRDF100C],"PRDF100C");

  // nitrogen loss to formation of products that decompose in 100 years
  strcpy(predstr[I_PRDF100N],"PRDF100N");


// Carbon and nitrogen fluxes associated with agriculture *******

  // agricultral net primary production
  strcpy(predstr[I_AGNPPC],"AGNPPC");

  // nitrogen uptake by crops associated with AGNPPC
  strcpy(predstr[I_AGNPPN],"AGNPPN");

  // carbon loss to formation of agricultural products
  strcpy(predstr[I_AGFPRDC],"AGFPRODC");

  // nitrogen loss to formation of agricultural products
  strcpy(predstr[I_AGPRDN],"AGFPRODN");

  strcpy(predstr[I_AGFRTN],"AGFERTN"); // nitrogen fertilization
  strcpy(predstr[I_AGLTRC],"AGLTRFC"); // litterfall carbon from crops
  strcpy(predstr[I_AGLTRN],"AGLTRFN"); // litterfall nitrogen from crops


// Carbon and nitrogen fluxes associated with human products ****

  // carbon loss to the formation of all products

  strcpy(predstr[I_TOTFPRDC],"TOTFPRDC");

  // nitrogen loss to the formation of all products

  strcpy(predstr[I_TOTFPRDN],"TOTFPRDN");

  // carbon loss to resulting from decomposition of agricultural
  //   products

  strcpy(predstr[I_AGPRDFC],"AGPRODFC");

  // nitrogen loss resulting from decomposition of agricultural
  //   products

  strcpy(predstr[I_AGPRDFN],"AGPRODFN");

  // carbon loss resulting from decomposition of PROD10C

  strcpy(predstr[I_PRD10FC],"PRD10FC");

  // nitrogen loss resulting from decomposition of PROD10N

  strcpy(predstr[I_PRD10FN],"PRD10FN");

  // carbon loss resulting from decomposition of PROD100C

  strcpy(predstr[I_PRD100FC],"PRD100FC");

  // nitrogen loss resulting from decomposition of PROD100N

  strcpy(predstr[I_PRD100FN],"PRD100FN");

  // carbon loss resulting from decomposition of all products

  strcpy(predstr[I_TOTPRDFC],"TOTPRDFC");

  // nitrogen loss resulting from decomposition of all products

  strcpy(predstr[I_TOTPRDFN],"TOTPRDFN");


// Integrated carbon fluxes *************************************

  // total net primary production (NPP+AGNPPC)

  strcpy(predstr[I_TOTNPP],"TOTNPP");

  // carbon flux from ecosystem (NEP+CONVERTC)

  strcpy(predstr[I_CFLX],"CFLUX");
  
  // soil NH4 content

  strcpy(predstr[I_NH4],"NH4");

  // soil NO3 content

  strcpy(predstr[I_NO3],"NO3");

  // N2O flux

  strcpy(predstr[I_N2O],"N2O");

  // N2 flux

  strcpy(predstr[I_N2],"N2");

  // NOx flux

  strcpy(predstr[I_NOX],"NOX");
  
    // Total N gas flux

  strcpy(predstr[I_NGAS],"NGAS");
  
    // Nitrification rate

  strcpy(predstr[I_FNIT],"FNIT");  
    
    // Denitrification rate

  strcpy(predstr[I_FDENIT],"FDENIT");
    
    // Final N2O flux
    
  strcpy(predstr[I_N2OF],"N2OF");
  
    // Final N2 production 
    
  strcpy(predstr[I_N2F],"N2F");
  
    // Atmospheric source N2O 
    
  strcpy(predstr[I_N2OAIR],"N2OAIR");
   
    // Atmospheric source N2 
    
  strcpy(predstr[I_N2AIR],"N2AIR");

    // Atmospheric source N2O uptake 
    
  strcpy(predstr[I_N2OUPT],"N2OUPT");
    
    // N2O from nitrification 
    
  strcpy(predstr[I_N2ON],"N2ON");
  
    // N2O from denitrification
    
  strcpy(predstr[I_N2ODN],"N2ODN");
  
    // NH3 volatilization 
    
  strcpy(predstr[I_NH3],"NH3");
    // N2 fixzation 
    
  strcpy(predstr[I_N2FIX],"N2FIX");  
    // N2 fixzation 
    
  strcpy(predstr[I_N2FIXS],"N2FIXS");  
    // N2 fixzation 
    
  strcpy(predstr[I_N2FIXFL],"N2FIXFL");    
    // DOC
  strcpy(predstr[I_DOC],"DOC");
};

/* **************************************************************
************************* Public Functions **********************
************************************************************** */


/* **************************************************************
************************************************************** */
//Rh and Nmin removed by CFD 20000616

int TTEM::boundcon(double ptstate[], double err[], double& ptol)
{

  int test = ACCEPT;

// Check carbon and nitrogen state variables
  if (err[I_VEGC] > fabs(ptol * ptstate[I_VEGC]))
  {
    return test = temkey(I_VEGC)+1;
  }
  if (nfeed == 1 && err[I_STRN] > fabs(ptol * ptstate[I_STRN]))
  {
    return test = temkey(I_STRN)+1;
  }
  if (err[I_SOLC] > fabs(ptol * ptstate[I_SOLC]))
  {
    return test = temkey(I_SOLC)+1;
  }
  if (nfeed == 1 && err[I_SOLN] > fabs(ptol * ptstate[I_SOLN]))  
  {
    return test = temkey(I_SOLN)+1; //printf("ttem434: ptsoln=%.3f\n", ptstate[I_SOLN]);
  }
  if (nfeed == 1 && err[I_AVLN] > fabs(ptol * ptstate[I_AVLN]))
  {
    return test = temkey(I_AVLN)+1;
  }
  if (err[I_GPP] > fabs(ptol * ptstate[I_GPP])) 
  {
    return test = temkey(I_GPP)+1;   printf("ttem442: ptgpp=%.3f\n",ptstate[I_GPP]);
  }
  if (err[I_NPP] > fabs(ptol * ptstate[I_NPP]))
  {
    return test = temkey(I_NPP)+1;
  }
  if (nfeed == 1 && err[I_VNUP] > fabs(ptol * ptstate[I_VNUP]))
  {
    return test = temkey(I_VNUP)+1;
  }
  if (nfeed == 1 && err[I_VSUP] > fabs(ptol * ptstate[I_VSUP]))
  {
    return test = temkey(I_VSUP)+1;
  }
  if (nfeed == 1 && err[I_STON] > fabs(ptol * ptstate[I_STON]))
  {
    return test = temkey(I_STON)+1;
  }
  if (nfeed == 1 && err[I_VNMBL] > fabs(ptol * ptstate[I_VNMBL]))
  {
    return test = temkey(I_VNMBL)+1;
  }
  //20000616 - CFD added RvMaint
  if (nfeed == 1 && err[I_RVMNT] > fabs(ptol * ptstate[I_RVMNT]))
  {
    return test = temkey(I_RVMNT)+1;
  }
  //20000616 - CFD removed Rh and Nmin

  // Check water state variables

  if (err[I_AVLW]  > fabs(ptol * ptstate[I_AVLW]))
  {
    return test = temkey(I_AVLW)+1;
  }
  if (err[I_RGRW]  > fabs(ptol * ptstate[I_RGRW]))
  {
    return test = temkey(I_RGRW)+1;
  }
  if (err[I_SNWPCK]  > fabs(ptol * ptstate[I_SNWPCK]))
  {
    return test = temkey(I_SNWPCK)+1;
  }
  if (err[I_SGRW]  > fabs(ptol * ptstate[I_SGRW]))
  {
    return test = temkey(I_SGRW)+1;
  }
  if (err[I_RPERC]  > fabs((ptol) * ptstate[I_RPERC]))
  {
    return test = temkey(I_RPERC)+1;
  }
  if (err[I_EET]  > fabs((ptol) * ptstate[I_EET]))
  {
    return test = temkey(I_EET)+1;
  }


  return test;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TTEM::delta(const int& dm, double pstate[], double pdstate[])
{

  nfert = 1.0000000;

  atms.eet[dm] = atms.xeet(atms.rain[dm], soil.snowinf[dm], atms.pet[dm],
                           pstate[I_AVLW], soil.awcapmm, dm);

  soil.percol(atms.rain[dm], soil.snowinf[dm], atms.eet[dm],
              pstate[I_AVLW],dm);

  if ((pstate[I_AVLW] + soil.snowinf[dm] + atms.rain[dm] - atms.eet[dm]
     - soil.rperc[dm] - soil.sperc[dm]) < 0.0)
  {
    atms.eet[dm] = pstate[I_AVLW] + soil.snowinf[dm] + atms.rain[dm]
                   - soil.rperc[dm] - soil.sperc[dm];
  }


  soil.rrun[dm] = soil.rrunoff(pstate[I_RGRW], soil.rperc[dm]);
  soil.srun[dm] = soil.srunoff(elev, atms.tair[dm], atms.prevtair,
                               atms.prev2tair, pstate[I_SGRW], soil.sperc[dm]);

  soil.h2oyld[dm] = soil.rrun[dm] + soil.srun[dm];


  if (moistlim == 0)
  {
    veg.unnormleaf[dm] = veg.deltaleaf(veg.cmnt, atms.pet[dm],
                                       atms.prvpetmx, prevy[I_UNRMLF]);
    //printf("ttem540:unnormleaf=%.3f,atms.pet=%.3f, prvpetmx=%.3f,prevy_UNRMLF=%.3f\n",veg.unnormleaf[dm],atms.pet[dm],atms.prvpetmx,prevy[I_UNRMLF]);                           
  }
  else
  {
    veg.unnormleaf[dm] = veg.deltaleaf(veg.cmnt, atms.eet[dm],
                                       atms.prveetmx, prevy[I_UNRMLF]);
    //printf("ttem546:unnormleaf=%.3f,atms.eet=%.3f, prvpetmx=%.3f,prevy_UNRMLF=%.3f\n",veg.unnormleaf[dm],atms.eet[dm],atms.prvpetmx,prevy[I_UNRMLF]);    
  }
  if (veg.prvleafmx[veg.cmnt] <= 0.0) { veg.leaf[dm] = 0.0; }
  else { veg.leaf[dm] = veg.unnormleaf[dm]/veg.prvleafmx[veg.cmnt]; }   
  //printf("ttem547: leaf=%.3f,unnormleaf=%.3f,prvleafmx=%.3f\n", veg.leaf[dm],veg.unnormleaf[dm],veg.prvleafmx[veg.cmnt]);
  if (veg.leaf[dm] < veg.minleaf[veg.cmnt])
  {
    veg.leaf[dm] = veg.minleaf[veg.cmnt];   
  }
  if (veg.leaf[dm] > 1.0) { veg.leaf[dm] = 1.0; }

  veg.alleaf = veg.leafmxc[veg.cmnt]/(1.0 + veg.kleafc[veg.cmnt]
               * exp(veg.cov[veg.cmnt]*pstate[I_VEGC]));
  veg.lai[dm] = veg.sla[veg.cmnt] * veg.alleaf * veg.leaf[dm];  // modified 02/02/2002 Kick and Qianlai
  veg.foliage = veg.alleaf / veg.leafmxc[veg.cmnt];
  veg.fpc[dm] = 1.0 - exp(-0.5*veg.lai[dm]);   // modified by Kick and Qialai 02/02/2002
   
  deltaxclm(veg.cmnt, soil.pcfldcap, dm);
  
  if (ag.state == 0)
  {

      // added veg.thawpercent for thaw-frozen mechanism by qianlai 28/08/2000
    veg.ingpp[dm] = veg.gppxclm(veg.cmnt, atms.co2[dm], atms.par[dm],
                            temp, gv, veg.leaf[dm], veg.foliage,veg.thawpercent[dm]); //,veg.ac 
//    veg.ingpp[dm] = veg.gppxclm(veg.cmnt, atms.co2[dm], atms.par[dm],
//                            temp, gv, veg.leaf[dm], veg.foliage);

    //printf("ttem 571: dm=%i, ingpp = %.3f\n", dm, veg.ingpp[dm]);
    
    if (veg.ingpp[dm] < 0.0) { veg.ingpp[dm] = 0.0; }

    veg.inuptake = veg.nupxclm(veg.cmnt, pstate[I_SM], pstate[I_AVLN],
                               respq10, ksoil, veg.foliage);
    //printf("ttem 577: dm=%i,inuptake = %.3f, Psm = %.3f, Pavln = %.3f, respq10 = %.3f, ksoil = %.3f, foliage = %.3f\n", dm, veg.inuptake, pstate[I_SM], pstate[I_AVLN], respq10, ksoil, veg.foliage);
  }
  else
  {
    veg.ingpp[dm] = 0.0;
    veg.inuptake = 0.0;
  }
  microbe.rh[dm] = microbe.rhxclm(pstate[I_SOLC], dq10, rhmoist);
  //printf("ttem 585: dm=%i, rh = %.3f,soc=%.3f, dq10=%.3f,rhmoist=%.3f\n", dm, microbe.rh[dm],pstate[I_SOLC],dq10,rhmoist);
  if (microbe.rh[dm] < 0.0) { microbe.rh[dm] = 0.0; }

  soil.ninput[dm] = atms.nhxdepo[dm] + atms.noydepo[dm]; //printf("ttem588: dm=%i, ninput=%.3f,nhx=%.3f,noy=%.3f\n",dm,soil.ninput[dm],atms.nhxdepo[dm],atms.noydepo[dm]);

  microbe.netnmin[dm] = microbe.nminxclm(veg.cmnt, dm, pstate[I_SM],
                                         pstate[I_SOLC], pstate[I_SOLN],
                                         pstate[I_AVLN], microbe.decay,
                                         microbe.rh[dm], ksoil);
  //printf("ttem 594: dm=%i, netmin=%.3f, ksoil=%.3f, PSOLN=%.3f, aviln=%.3f\n",dm, microbe.netnmin[dm], ksoil, pstate[I_SOLN],pstate[I_AVLN]);
  //printf("ttem 595: vegcmnt=%i, pSM=%.3f, PSOLC=%.3f, decay=%.3f, rh=%.3f,ksoil=%.3f\n", veg.cmnt, pstate[I_SM], pstate[I_SOLC], microbe.decay,microbe.rh[dm],ksoil);
  if (ag.state == 0)
  {
    veg.ltrfal[dm].carbon = veg.cfall[veg.cmnt] * pstate[I_VEGC];
    if (veg.ltrfal[dm].carbon < 0.0) { veg.ltrfal[dm].carbon = 0.0; }

    veg.ltrfal[dm].nitrogen = veg.nfall[veg.cmnt] * pstate[I_STRN];  //printf("ttem601:dm=%i, ltrfaln=%.3f,nfall=%.3f,pdstrn=%.3f\n",dm,veg.ltrfal[dm].nitrogen,veg.nfall[veg.cmnt],pstate[I_STRN]);
    if (veg.ltrfal[dm].nitrogen < 0.0) { veg.ltrfal[dm].nitrogen = 0.0; }

    veg.rm[dm] = veg.rmxclm(veg.cmnt,pstate[I_VEGC],respq10); 
    if (veg.rm[dm] < 0.0) { veg.rm[dm] = 0.0; }

    veg.innpp[dm] = veg.ingpp[dm] - veg.rm[dm];  //printf("ttem611: innpp=%.3f,ingpp=%.3f,rm=%.3f\n", veg.innpp[dm],veg.ingpp[dm],veg.rm[dm]);

    veg.rg[dm] = 0;
    if (veg.innpp[dm] > 0.0)
    {
      veg.rg[dm]  = 0.2 * veg.innpp[dm];
      veg.innpp[dm] *= 0.8;  //printf("ttem617: innpp=%.3f\n", veg.innpp[dm]);
    }

    if (veg.inuptake > pstate[I_AVLN] + (nfert*soil.ninput[dm]) + microbe.netnmin[dm] + microbe.N2Fix[dm])
    {
      veg.inuptake = pstate[I_AVLN] + (nfert*soil.ninput[dm]) + microbe.netnmin[dm] + microbe.N2Fix[dm];
	    //printf("ttem 623: inuptake = %.3f, pAVLN = %.3f, nfert = %.3f, ninput = %.3f, netnmin = %.3f,nfix=%.3f\n", veg.inuptake, pstate[I_AVLN], nfert, soil.ninput[dm], microbe.netnmin[dm],microbe.N2Fix[dm]);
    }

    if (veg.inuptake < 0.0) { veg.inuptake = 0.0; } 
    //if(fixflag == 0) {veg.nuptake[dm] = veg.inuptake;}     //101922 stable vegetation N uptake
    veg.nuptake[dm] = veg.inuptake;  //printf("ttem 624: dm=%i,vegnuptake = %.3f,inuptake=%.3f\n", dm,veg.nuptake[dm],veg.inuptake);
    veg.suptake[dm] = veg.nuptake[dm];  //printf("ttem 625: suptake = %.3f\n", veg.suptake[dm]);
    veg.luptake[dm] = 0.0;
    veg.gpp[dm] = veg.ingpp[dm]; //printf("ttem 631: dm=%i,veggpp = %.3f\n", dm,veg.gpp[dm]);
    veg.npp[dm] = veg.innpp[dm]; //printf("ttem 632: dm=%i,vegnpp = %.3f\n", dm,veg.npp[dm]);
    veg.nmobil[dm] = 0.0;
    veg.nresorb[dm] = 0.0;

// Nitrogen feedback of GPP (nfeed == 1)

    if (nfeed == 1)
    {
      if (veg.inuptake == 0.0) { veg.inuptake = 0.000001; }
      veg.inprodcn = veg.innpp[dm] / (veg.inuptake + pstate[I_STON]); //
      if (veg.inprodcn <= 0) { veg.inprodcn = veg.cnmin[veg.cmnt]; }
      //printf("ttem638: inprodcn=%.3f,innpp=%.3f,inuptake=%.3f,psston=%.3f\n", veg.inprodcn,veg.innpp[dm],veg.inuptake,pstate[I_STON]);
      if (veg.ltrfal[dm].nitrogen <= veg.ltrfal[dm].carbon/veg.cneven[veg.cmnt]) 
      {
      veg.nresorb[dm] = veg.ltrfal[dm].carbon/veg.cneven[veg.cmnt]- veg.ltrfal[dm].nitrogen;
      }
      else
      {
      veg.ltrfal[dm].nitrogen = veg.ltrfal[dm].carbon/veg.cneven[veg.cmnt];
      veg.nresorb[dm] = 0.0;
      }
      
      if (pstate[I_VEGC] > 0.0) {veg.nresorb[dm] *= (pstate[I_STRN]/pstate[I_VEGC]) * veg.c2n;}

      if (veg.inprodcn > veg.cneven[veg.cmnt])
      { 
        veg.npp[dm] = veg.cneven[veg.cmnt] * (veg.nuptake[dm] + pstate[I_STON]);
        if (veg.npp[dm] < 0.0) { veg.npp[dm] = 0.0; }
 	      veg.rg[dm] = 0.25 * veg.npp[dm];
  	    veg.gpp[dm] = veg.npp[dm] + veg.rg[dm] + veg.rm[dm];  
        //printf("ttem 661:gpp=%.3f, npp=%.3f,veg.rg=%.3f,veg.rm=%.3f\n",veg.gpp[dm],veg.npp[dm],veg.rg[dm],veg.rm[dm]);
   
	      if (veg.gpp[dm] < 0.0) { veg.gpp[dm] = 0.0; } 
	      veg.nmobil[dm] = pstate[I_STON]; //printf("ttem659: nmobile=%.3f,ston=%.3f\n",veg.nmobil[dm],pstate[I_STON]);
      }
      if (veg.inprodcn <= veg.cneven[veg.cmnt]) 
      {
       //if(fixflag == 0)   {                                                         //101922 stable vegetation N uptake
        veg.nuptake[dm] = veg.inuptake * (veg.inprodcn - veg.cnmin[veg.cmnt])
                          * (veg.inprodcn - 2*veg.cneven[veg.cmnt] + veg.cnmin[veg.cmnt]);
	      veg.nuptake[dm] /= ((veg.inprodcn - veg.cnmin[veg.cmnt]) * (veg.inprodcn - 2*veg.cneven[veg.cmnt] + veg.cnmin[veg.cmnt])) 
                           - pow(veg.inprodcn - veg.cneven[veg.cmnt],2.0); 
        //printf("ttem 672: vegnuptake = %.3f,veginuptake=%.3f, inprodcn=%.3f,cneven=%.3f\n", veg.nuptake[dm],veg.inuptake,veg.inprodcn,veg.cneven[veg.cmnt]);
        // }
	     if (veg.nuptake[dm] < 0.0 ) { veg.nuptake[dm] = 0; }//printf("ttem675: vegnuptake = %.3f\n", veg.nuptake[dm]);    //&& fixflag == 0
       if (pstate[I_STON] >= veg.npp[dm]/veg.cneven[veg.cmnt])
       {
	       veg.nmobil[dm] = veg.npp[dm]/veg.cneven[veg.cmnt];
	       if (veg.nmobil[dm] < 0.0 && pstate[I_VEGC] > 0.0){ veg.nmobil[dm] *= (pstate[I_STRN]/pstate[I_VEGC]) * veg.c2n;}
	       veg.suptake[dm] = 0.0;
	     }
     	 else
       {
       veg.nmobil[dm] = pstate[I_STON];
	     veg.suptake[dm] = (veg.npp[dm]/veg.cneven[veg.cmnt]) - veg.nmobil[dm];
		   //printf("ttem 693: suptake = %.3f, npp = %.3f, cneven = %.3f, nmobil = %.3f\n", veg.suptake[dm], veg.npp[dm], veg.cneven[veg.cmnt], veg.nmobil[dm]);
	     if (veg.suptake[dm] < 0.0) { veg.suptake[dm] = 0.0; }
	     if (veg.suptake[dm] > veg.nuptake[dm]){veg.suptake[dm] = veg.nuptake[dm];} //printf("ttem687: vegnuptake = %.3f\n", veg.nuptake[dm]);
       }
		   //printf("ttem 689: suptake = %.3f, nuptake = %.3f\n", veg.suptake[dm], veg.nuptake[dm]);
       // ptate[1] changed to pstate[I_STRN] by DWK on 20000201
	     if ((pstate[I_STON] + veg.nuptake[dm] - veg.suptake[dm]
           + veg.nresorb[dm] - veg.nmobil[dm]) < (veg.labncon[veg.cmnt]
           * (pstate[I_STRN] + veg.suptake[dm] - veg.ltrfal[dm].nitrogen
           - veg.nresorb[dm] + veg.nmobil[dm])))
       {
       veg.luptake[dm] = veg.nuptake[dm] - veg.suptake[dm];
	     }
	     else
       {
	     veg.luptake[dm] = (veg.labncon[veg.cmnt] * (pstate[I_STRN]
                             + veg.suptake[dm] - veg.ltrfal[dm].nitrogen
                             - veg.nresorb[dm] + veg.nmobil[dm]))
                             - (pstate[I_STON] + veg.nresorb[dm]
                             - veg.nmobil[dm]);
	     if (veg.luptake[dm] < 0.0) { veg.luptake[dm] = 0.0; }
	   //if(fixflag == 0) {
       veg.nuptake[dm] = veg.suptake[dm] + veg.luptake[dm]; 
 //}     //101922 stable vegetation N uptake 
	     }
      //printf("ttem 706: vegnupta=%.3f,suptake=%.3f,luptake=%.3f\n", veg.nuptake[dm],veg.suptake[dm], veg.luptake[dm]);  
      if (veg.nuptake[dm] <= 0.0 ) { 
         veg.nuptake[dm] = veg.innpp[dm]/veg.cneven[veg.cmnt]/1.3;
         if (veg.inuptake > pstate[I_AVLN] + (nfert*soil.ninput[dm]) + microbe.netnmin[dm] + microbe.N2Fix[dm])
            {veg.inuptake = pstate[I_AVLN] + (nfert*soil.ninput[dm]) + microbe.netnmin[dm] + microbe.N2Fix[dm];}}
      // printf("ttem 716: vegnupta=%.3f,suptake=%.3f,luptake=%.3f\n", veg.nuptake[dm],veg.suptake[dm], veg.luptake[dm]);     
      }
    }
    //printf("ttem711: vegnuptake = %.3f\n", veg.nuptake[dm]);
    veg.gpr[dm] = veg.rm[dm] + veg.rg[dm];
    nep[dm] = veg.npp[dm] - microbe.rh[dm];  //printf("ttem721:nep=%.3f,npp=%.3f,rh=%.3f\n",nep[dm],veg.npp[dm],microbe.rh[dm]);
    cflux[dm] = nep[dm];

    ag.npp[dm].carbon = 0.0;
    ag.npp[dm].nitrogen = 0.0;
    ag.fertn[dm] = 0.0;
    ag.ltrfal[dm].carbon = 0.0;
    ag.ltrfal[dm].nitrogen = 0.0;
    ag.slash[dm].carbon = 0.0;
    ag.slash[dm].nitrogen = 0.0;
    ag.sconvrtflx[dm].carbon = 0.0;
    ag.sconvrtflx[dm].nitrogen = 0.0;
    ag.nsretent[dm] = 0.0;
  }
  else
  {
    veg.ltrfal[dm].carbon = 0.0;
    veg.ltrfal[dm].nitrogen = 0.0;
    veg.gpr[dm] = 0.0;
    veg.rm[dm] = 0.0;
    veg.rg[dm] = 0.0;
    veg.innpp[dm] = 0.0;
    veg.nuptake[dm] = 0.0;
    veg.suptake[dm] = 0.0;
    veg.luptake[dm] = 0.0;
    veg.gpp[dm] = 0.0;  //printf("ttem 746: gpp=%.3f\n",veg.gpp[dm]);
    veg.npp[dm] = 0.0;
    veg.nmobil[dm] = 0.0;
    veg.nresorb[dm] = 0.0;

    ag.npp[dm].carbon = ag.RAP * ag.potnpp[dm];
    if (ag.npp[dm].carbon < 0.0) { ag.npp[dm].carbon = 0.0; }
    ag.npp[dm].nitrogen = ag.npp[dm].carbon / ag.c2n;
    ag.fertn[dm] = 0.0;
    ag.ltrfal[dm].carbon = ag.npp[dm].carbon * ag.cfall[veg.cmnt];
    ag.ltrfal[dm].nitrogen = ag.npp[dm].nitrogen * ag.nfall[veg.cmnt];
    nep[dm] = ag.npp[dm].carbon - microbe.rh[dm];
    cflux[dm] = nep[dm] - ag.convrtflx[dm].carbon;
  }

// Changes in available nitrogen in soil
double c_no3,c_nh4;
c_no3=pstate[I_NO3];
c_nh4=pstate[I_NH4];

	
  if (avlnflag == 1)
  {
  	if (days[dm]==0) days[dm]=30;
  		 // for N cycle
 //****************************************************************************/
    	int m;
//     double CDM;
//     CDM = 0.0;
//     for (m = 0; m < CYCLE; m++)
//     CDM = CDM + (10.0 - atms.tair[m]);

//   	int pday;
//   	switch (dm){
//   		case 0: pday = 1;  break;
//   		case 1: pday = 32; break;
//   		case 2: pday = 60; break;
//   		case 3: pday = 91; break;
//   		case 4: pday = 121; break;
//   		case 5: pday = 152; break;
//   		case 6: pday = 182; break;
//   		case 7: pday = 213; break;
//   		case 8: pday = 244; break;	
//   		case 9: pday = 274; break;
//   		case 10: pday = 305; break;
//   		case 11: pday = 335; break;
//   	}
   	// daily version for HM, STM and N cycle; the daily environmental condition from getenviron() function;
	microbe.N2O[dm] = 0.0;
	microbe.NOx[dm] = 0.0;
	microbe.N2[dm] = 0.0;
	microbe.Fnit[dm] = 0.0;
	microbe.Fdenit[dm] = 0.0;
  microbe.N2OF[dm] = 0.0;  // YYuan 010122 for N2O uptake
  microbe.N2F[dm] = 0.0;  // YYuan 010122 for N2O uptake
  microbe.N2OAir[dm] = 0.0;   // YYuan 010122 for N2O uptake
  microbe.N2Air[dm] = 0.0;  // YYuan 010122 for N2O uptake
  microbe.N2OUpt[dm] = 0.0;  // YYuan 010122 for N2O uptake
  microbe.N2ON[dm] = 0.0;  // 
  microbe.N2ODN[dm] = 0.0;  // 
  microbe.NH3[dm] = 0.0;  // YYuan 033122 for NH3 volatilization
  microbe.N2Fix[dm]= 0.0;
  microbe.N2Fixs[dm]= 0.0;
  microbe.N2Fixfl[dm]= 0.0;
	daynit = 0.0;
	daydenit = 0.0;
   for (m = 0; m < days[dm]; m++){
   
//     		//use HM daily step added by YYL
//   	   // convert the monthly data to daily
//
//   	    daylai[dm][m]=veg.lai[dm];
//   	    daynirr[dm][m]=atms.nirr[dm];
//   	   
//		    daysnowpack[dm][m]=soil.snowpack[dm];
//		
//   	  // Determine if daily precipitation occurs as rain or snow   
// //  	 atms.precsplt(atms.prec[dm][m],atms.tair[dm][m],atms.rain[dm][m],atms.snowfall[dm][m]);
//			atms.precsplt(dayprec[dm][m],daytair[dm][m],dayrain[dm][m],daysnowfall[dm][m]);
// //  if (atms.GISlai[dm][m] == 0.0) atms.GISlai[dm][m] = 0.01;  // for drbugging Q. Z.
// //  if (atms.GISlai[dm][m] < 0.0) atms.GISlai[dm][m] = 0.5;  // to deal with non-data grid cells
// //  veg.lai[dm][m] = atms.GISlai[dm][m]; // for calibration version
//
//    hyd.vaporpressure[dm][m] = dayvap[dm][m]/10.0;  //use CRU data sets, converse to kpa
// //  hyd.vpd(hyd.vaporpressure[dm][m], atms.tair[dm][m]);
//		hyd.vpd(hyd.vaporpressure[dm][m], daytair[dm][m]);
//
// //  hyd.rainthrough[dm][m] = hyd.intercept(atms.rain[dm][m], hyd.inter_coef[veg.cmnt],dm, m);
//   hyd.rainthrough[dm][m] = hyd.intercept(dayrain[dm][m], hyd.inter_coef[veg.cmnt],dm, m);
//
//   soil.avlh2o1[dm][m] = 2000.0 * (hydm.h2osoi[0] + hydm.h2osoi[1] + hydm.h2osoi[2] + hydm.h2osoi[3] + hydm.h2osoi[4]
//                          + hydm.h2osoi[5])/6.0;
//
////   soil.avlh2o1[dm][m] =  soil.awcapmm + hyd.rainthrough[dm][m];
////   soil.avlh2o2[dm][m] =  soil.awcapmm + hyd.rainthrough[dm][m];
//
//   hyd.canopy_water[dm][m] = hyd.prec_to_canopy[dm][m];
//
//   hyd.drad(daylai[dm][m], hyd.EXT[veg.cmnt], daynirr[dm][m], dm, m);
//
//   hyd.hLWP[dm][m] =  hyd.LWP(soil.avlh2o1[dm][m], soil.fldcap, soil.pctsand, soil.pctclay, soil.pctsilt);
//
//   hyd.hpa[dm][m] = hyd.densityA(daytair[dm][m]);
//   hyd.hlv[dm][m] = hyd.latentV(daytair[dm][m]);
//
//  //rainfall is changed to rainthrough fall due to interception
//  //  atms.rain[dm][m] = hyd.rainthrough[dm][m];
//
//   hyd.Gs[dm][m] = hyd.CanopyConductance(daylai[dm][m], hyd.gmax[veg.cmnt],
//           daytair[dm][m], hyd.vpdeficit, hyd.hLWP[dm][m],
//           hyd.vpd_close[veg.cmnt],hyd.vpd_open[veg.cmnt], hyd.psi_close[veg.cmnt], hyd.psi_open[veg.cmnt]);
//
//   hyd.htrans[dm][m] = hyd.penmanmonteith_new(daytair[dm][m], hyd.hdrad[dm][m], hyd.hpa[dm][m],
//             hyd.vpdeficit, hyd.hlv[dm][m], hyd.Gs[dm][m], daylai[dm][m]);
//
//   hyd.potevap_surface[dm][m] = hyd.evaporation_layer1(daytair[dm][m],
//                                 hyd.hpa[dm][m],hyd.hlv[dm][m],hyd.vpdeficit, hyd.surfacerad[dm][m]);
//
// // consider only the precipitation flux reaching the soil
// // check for precipitation >= potential evaporation
//
////   	if (hyd.rainthrough[dm][m] >= hyd.potevap_surface[dm][m])
////         	hyd.soil_evap[dm][m] = hyd.potevap_surface[dm][m];
////   	else	hyd.soil_evap[dm][m] = 0.0;
//
//  	hyd.soil_evap[dm][m] = hyd.potevap_surface[dm][m];
//
// //snow intercept
//    if (daysnowpack[dm][m] > 0.0) {
//          hyd.sub_from_canopy[dm][m] = hyd.sublimation_canopy(daylai[dm][m],
//              hyd.hdrad[dm][m], daytair[dm][m], daysnowpack[dm][m], dm, m);
//          if (hyd.sub_from_canopy[dm][m] == 0.0)  hyd.sub_from_canopy[dm][m] = 0.001; // to protect
//     }
//
//    soil.snowthrough[dm][m] = hyd.snowthroughfall[dm][m];
//    hyd.snowmelt[dm][m] = hyd.snowmelt_srad(hyd.hdrad[dm][m], daytair[dm][m], soil.snowthrough[dm][m]);
//
// // snowpack is updated to snowthroug[dm][m]
//    daysnowpack[dm][m] = soil.snowthrough[dm][m] - hyd.snowmelt[dm][m];
//
// // sublimation from ground snow
//    hyd.snowsub[dm][m] = hyd.snowsublimation_ground(soil.snowthrough[dm][m], daytair[dm][m], hyd.surfacerad[dm][m]);
//
// // the estimate of actual evapotranspiration
//    dayeet3[dm][m] = hyd.htrans[dm][m] + hyd.sub_from_canopy[dm][m] + hyd.soil_evap[dm][m] + hyd.snowsub[dm][m];
//
// // the following is the original one to estimate the potential evapotranspiration
//    daypet3[dm][m] = atms.new_petjh(daynirr[dm][m], daytair[dm][m], dm);
//    if (daypet3[dm][m] <= 0.0) { daypet3[dm][m] = 0.001; }
//
// // deal with soil hydrological dynamics
//    hydm.surrunoff[dm][m] = hydm.surfacerunoff(hyd.snowmelt[dm][m],hyd.rainthrough[dm][m],hydm.h2osoi[0],
//            hydm.b1[veg.cmnt], hydm.ksat1[veg.cmnt], hydm.thetasat1[veg.cmnt], hydm.presat1[veg.cmnt],
//            hydm.thick1[veg.cmnt]);
//   hydm.surinfl[dm][m] = hydm.infiltration(hydm.surrunoff[dm][m], hyd.snowmelt[dm][m],hyd.rainthrough[dm][m]);
//
// // calculate soil water content and sub-surafce runoff (drainage)
//    hydm.qtran = hyd.htrans[dm][m];
//    hydm.qseva = hyd.soil_evap[dm][m];
// //  hydm.bdrai[dm][m] = hydm.qdrai[dm][m];
//
// // for general soil
//    hydm.main_soil(hydm.nsl, hyd.htrans[dm][m], hydm.surinfl[dm][m], hyd.soil_evap[dm][m], hydm.bdrai[dm][m], hydm.dtsoi);
//
////	for wetland
////			hydm.wetland_soil(hydm.nsl, hyd.htrans[dm][m], hydm.surinfl[dm][m], hyd.soil_evap[dm][m], hydm.bdrai[dm][m], hydm.dtsoi);
//		
//		if (hydm.h2osoi[0] > soil.totpor1) hydm.h2osoi[0] = soil.totpor1;
//		if (hydm.h2osoi[1] > soil.totpor2) hydm.h2osoi[1] = soil.totpor2;					
//		if (hydm.h2osoi[2] > soil.totpor3) hydm.h2osoi[2] = soil.totpor3;		
//			
//		daypctp[dm][m] = 100 * (hydm.h2osoi[0] /(soil.totpor1 / 100.0) +  hydm.h2osoi[1] /(soil.totpor2 / 100.0) 
//											+  hydm.h2osoi[2] /(soil.totpor3 / 100.0))/3;
// //end of daily HM
////***********************************************************************************************************/
// // call soil thermal subroutine
//     switch (pday) {
//     case 1: sthermal.airt19= tempdaily[pday]; // satisfy soil thermal model
//             sthermal.airt29= tempdaily[pday];
//             sthermal.airt39= (tempdaily[pday]+tempdaily[pday+1])/2.0;
//             break;
//     case 365:sthermal.airt19= (tempdaily[pday-1]+tempdaily[pday])/2.0; // satisfy soil thermal model
//             sthermal.airt29= tempdaily[pday];
//             sthermal.airt39= tempdaily[pday];
//            break;
//     default:sthermal.airt19= (tempdaily[pday-1]+tempdaily[pday])/2.0; // satisfy soil thermal model
//             sthermal.airt29= tempdaily[pday];
//             sthermal.airt39= (tempdaily[pday]+tempdaily[pday+1])/2.0;
//            break;
//         }
//
//    switch (dm) {
//     case 0: sthermal.hsnow19=(soil.snowpack[CYCLE-1]+soil.snowpack[0])/2.0/1000.0; // satisfy soil thermal model
//             sthermal.hsnow29= soil.snowpack[0]/1000.0;
//             sthermal.hsnow39= (soil.snowpack[1]+soil.snowpack[0])/2.0/1000.0;
//             break;
//     case 11: sthermal.hsnow19=(soil.snowpack[11]+soil.snowpack[10])/2.0/1000.0; // satisfy soil thermal model
//              sthermal.hsnow29= soil.snowpack[11]/1000.0;
//              sthermal.hsnow39= (soil.snowpack[11]+soil.snowpack[0])/2.0/1000.0;
//            break;
//     default: sthermal.hsnow19=(soil.snowpack[dm-1]+soil.snowpack[dm])/2.0/1000.0; // satisfy soil thermal model
//              sthermal.hsnow29= soil.snowpack[dm]/1000.0;
//              sthermal.hsnow39= (soil.snowpack[dm]+soil.snowpack[dm+1])/2.0/1000.0;
//            break;
//         }
////   cout << airt19 << " " << airt29 << " " << airt39 << " " << hsnow19 << " " << hsnow29 << " " << hsnow39 <<endl;
//
//	sthermal.is9 = 10;
//    sthermal.ism19 = 9;
//    sthermal.smass9 = 0.f;
//
//   int i;
//
//	for (i=1;i <=210; i++) {
//    sthermal.t9[i] = 0;
//  	sthermal.xfa9[i] = -1e10f;
//	sthermal.xfb9[i] = 0.f;
//	sthermal.water9[i] = 1.f;
//	sthermal.dx9[i] = 0.f;
//	sthermal.weight9[i] = 0.f; }
//
//
//   sthermal.calcu_cdsnow = 0.2; // constant for equilibrium run --- an assumpition
//   sthermal.water2 = daypctp[dm][m]/100.0; // assume moss and organic have the same water content, 29/march/2001
//
//  if (sthermal.soiltemp_(&sthermal.water2, &sthermal.calcu_cdsnow, &sthermal.airt19, &sthermal.airt29, &sthermal.airt39, &sthermal.hsnow19, &sthermal.hsnow29,
//    &sthermal.hsnow39,sthermal.weight9, &sthermal.smass9, &sthermal.is9, &sthermal.ism19, sthermal.x9, sthermal.dx9,
//     sthermal.xfb9, sthermal.xfa9,	sthermal.water9, sthermal.t9, &sthermal.tsoil,&sthermal.frontd,&sthermal.thawbegin,
//    &sthermal.thawend,sthermal.diffsoilt, veg.cmnt) !=0)
//    { printf("bad tem"); 
////        getch();
//		};
//
//   //      ++kswitch;
//				daytsoil[dm][m]=sthermal.tsoil;
////      atms.frontd[dm]=sthermal.frontd;
////      atms.thawbe[dm]=sthermal.thawbegin;
////      atms.thawend[dm]=sthermal.thawend;
////
////      atms.tsoil[dm]=sthermal.tsoil;
////      atms.dst5[dm]=sthermal.diffsoilt[0];
////      atms.dst10[dm]=sthermal.diffsoilt[1];
////      atms.dst20[dm]=sthermal.diffsoilt[2];
////      atms.dst50[dm]=sthermal.diffsoilt[3];
////      atms.dst100[dm]=sthermal.diffsoilt[4];
////      atms.dst200[dm]=sthermal.diffsoilt[5];
//
// /*** end of calling soil thermal model ***/
//*********************************************************************************************/	
// for N cycle to determine the N gas flux	

	daynit = microbe.Nitrif(veg.cmnt,microbe.netnmin[dm] /days[dm], daytsoil[dm][m], daypctp[dm][m], soil.psiplusc, ph, c_nh4 );
	daydenit = microbe.Denitrif(veg.cmnt,c_no3, pstate[I_DOC], daypctp[dm][m], ph, daytsoil[dm][m]);
	//printf("ttem 990: daynit = %.3f, netnmin=%.3f, pctp = %.3f, tsoil = %.3f, c_nh4=%.3f\n", daynit, microbe.netnmin[dm],daypctp[dm][m],daytsoil[dm][m], c_nh4 );
  //printf("ttem 991: dm =%i, daydenit =%.6f, c_no3 = %.3f, pDOC = %.3f\n", dm, daydenit,c_no3,pstate[I_DOC]);
	if (daynit <= 0.0) {daynit = 0.0;}
	if (daydenit <= 0.0){daydenit = 0.0;}

	microbe.diff = microbe.Diffusivity(soil.pcfldcap/100, soil.pctpor/100, daypctp[dm][m]*soil.pctpor/10000);
	microbe.rn2 = microbe.Rn2(microbe.diff, c_no3 , pstate[I_DOC], ph, daypctp[dm][m]);
  //printf("ttem974: dm =%i, rn2 =%.3f\n", dm, microbe.rn2);
	microbe.rnox = microbe.Rnox(microbe.diff);
	
  microbe.dayn2fixs[dm][m] = microbe.N2Fixation(veg.cmnt, daytsoil[dm][m], daypctp[dm][m], pstate[I_AVLN], pstate[I_VEGCBG],density);
  microbe.dayn2fixfl[dm][m] = microbe.N2FixationFL(veg.cmnt,daytsoil[dm][m], daypctp[dm][m],pstate[I_AVLN],pstate[I_SOLC]);
  microbe.dayn2fix[dm][m] = microbe.dayn2fixs[dm][m] + microbe.dayn2fixfl[dm][m]; //printf("ttem1002: fixs=%.3f, fixfl=%.3f\n",microbe.dayn2fixs[dm][m],microbe.dayn2fixfl[dm][m]);
  //printf("ttem1003: dm=%i,m=%i,fix=%.3f,veg=%i,tsoil= %.3f,pctp=%.3f\n",dm,m,microbe.dayn2fix[dm][m], veg.cmnt, daytsoil[dm][m],daypctp[dm][m]);
  //printf("ttem1004: availn=%.3f,vegcb=%.3f,,soilc=%.3f, density=%.3f\n",pstate[I_AVLN],pstate[I_VEGCBG],pstate[I_SOLC],density);
  
	microbe.noxpi = 1.0; // when using the monthly data input
	microbe.nflux(veg.cmnt,daynit, daydenit, microbe.rn2, microbe.rnox, microbe.noxpi, dm, m);
							//	microbe.dayn2[dm][m], microbe.dayn2o[dm][m], microbe.daynox[dm][m]);
  
  microbe.dayn2oupt[dm][m] = microbe.N2OUptake(microbe.dayn2o[dm][m],daysnowpack[dm][m]);    // YYuan 010122 for N2O uptake
  microbe.dayn2air[dm][m] = microbe.dayn2oupt[dm][m]*microbe.rn2;    //Atmosphere N2O to N2
  if (microbe.dayn2air[dm][m] > microbe.dayn2oupt[dm][m]) {microbe.dayn2air[dm][m]=microbe.dayn2oupt[dm][m];}
  microbe.dayn2oair[dm][m] = microbe.dayn2oupt[dm][m]-(microbe.dayn2oupt[dm][m]*microbe.rn2);  //Atmospheric N2O exist in soil after consumpion
  //printf("ttem 1014: dayn2oair = %.6f, dayn2oupt = %.6f, rn2 = %.6f\n", microbe.dayn2oair[dm][m], microbe.dayn2oupt[dm][m], microbe.rn2);
  microbe.dayn2f[dm][m] = microbe.dayn2[dm][m] + microbe.dayn2air[dm][m]; //total N2 production (soil bron + atmosphere born) 
  //printf("ttem 1016: dayn2f=%.3f, dayn2=%.3f,dayn2air=%.3f\n", microbe.dayn2f[dm][m], microbe.dayn2[dm][m],microbe.dayn2air[dm][m]);
  microbe.dayn2of[dm][m] = microbe.dayn2o[dm][m] - microbe.dayn2oupt[dm][m]; //N2O flux 
  //printf("ttem 1018: dayn2of=%.3f, dayn2o=%.3f,dayn2oair=%.3f\n", microbe.dayn2of[dm][m],microbe.dayn2o[dm][m],microbe.dayn2oupt[dm][m]);                                                                  
  microbe.daynh3[dm][m] = microbe.NH3Vola(daytsoil[dm][m],c_nh4);
  
	microbe.N2O[dm] += microbe.dayn2o[dm][m]; //printf("ttem1021: n2o=%.3f,dayn2o=%.3f\n",microbe.N2O[dm],microbe.dayn2o[dm][m]);
	microbe.NOx[dm] += microbe.daynox[dm][m];
	microbe.N2[dm] += microbe.dayn2[dm][m]; // printf("ttem1023: n2=%.3f,dayn2=%.3f\n",microbe.N2[dm],microbe.dayn2[dm][m]);
	microbe.Fnit[dm] = microbe.Fnit[dm] + daynit;
	microbe.Fdenit[dm] = microbe.Fdenit[dm] + daydenit;
 
 	microbe.N2OF[dm] += microbe.dayn2of[dm][m];
	microbe.N2F[dm] += microbe.dayn2f[dm][m];
	microbe.N2OAir[dm] += microbe.dayn2oair[dm][m];
	microbe.N2Air[dm] += microbe.dayn2air[dm][m];
  microbe.N2OUpt[dm] = microbe.N2OUpt[dm] + microbe.dayn2oupt[dm][m];
 	microbe.N2ON[dm] += microbe.dayn2on[dm][m];
 	microbe.N2ODN[dm] += microbe.dayn2odn[dm][m];
	microbe.N2Fix[dm] += microbe.dayn2fix[dm][m]; 
	microbe.N2Fixs[dm] += microbe.dayn2fixs[dm][m]; 
	microbe.N2Fixfl[dm] += microbe.dayn2fixfl[dm][m]; 
  microbe.NH3[dm] += microbe.daynh3[dm][m];
//	if (daypctp[dm][m]>microbe.N2[dm]) microbe.N2[dm] = daypctp[dm][m];

      if (c_nh4 + c_no3 == 0)	{
      	rnitrate = 0.5;
      	rammonium = 0.5;
      }	
			else{
				rammonium = c_nh4 / (c_nh4 + c_no3); //printf ("ttem 1044: rammonium=%.3f,c_nh4=%.2f, c_no3=%.3f\n", rammonium, c_nh4, c_no3);
				rnitrate = c_no3 / (c_nh4 + c_no3); //printf ("ttem 1045: rnitrate=%.3f\n", rnitrate);
			}
	c_nh4 += microbe.netnmin[dm]/days[dm] + atms.nhxdepo[dm]/days[dm]- daynit - veg.nuptake[dm]/days[dm] * rammonium + microbe.dayn2fix[dm][m] ; //- microbe.daynh3[dm][m]
	c_no3 += daynit + atms.noydepo[dm]/days[dm] - microbe.dayn2o[dm][m] - microbe.dayn2[dm][m] - veg.nuptake[dm]/days[dm] * rnitrate; // - microbe.daynox[dm][m]
	if (c_nh4 < 0)  c_nh4=0;
	if (c_no3 < 0)  c_no3=0;
  //printf ("ttem1051: c_nh4=%.3f,netmin=%.3f, daynit=%.3f, veg.nuptake=%.3f,rammonium=%.3f,dayfix=%.3f\n", c_nh4, microbe.netnmin[dm], daynit, veg.nuptake[dm], rammonium,microbe.dayn2fix[dm][m]);
  //printf ("ttem1052: c_no3=%.3f, daynit = %.3f, dayn2o=%.3f, daynox=%.3f, dayn2=%.3f, veguptake=%.3f, rnitrate=%.3f\n", c_no3, daynit,microbe.dayn2o[dm][m], microbe.daynox[dm][m], microbe.dayn2[dm][m], veg.nuptake[dm], rnitrate);
//	pday++;
  }//end of daily version
//********************************************************************************/
//  microbe.Fnit[dm] = microbe.Nitrif(microbe.netnmin[dm]*1000/days[dm], atms.tsoil[dm], soil.pctp[dm], soil.psiplusc, ph, pstate[I_NH4]*1000);
//	microbe.Fdenit[dm] = microbe.Denitrif(pstate[I_NO3]*1000, microbe.rh[dm]*1000/days[dm], soil.pctp[dm], atms.tsoil[dm]);
////  microbe.Fnit[dm] = microbe.Nitrif(microbe.netnmin[dm]*1000/30, atms.tsoil[dm], soil.pctp[dm], soil.psiplusc, ph, pstate[I_AVLN]*1000/2);
////	microbe.Fdenit[dm] = microbe.Denitrif(pstate[I_AVLN]*1000/2, microbe.rh[dm]*1000/30, soil.pctp[dm], atms.tsoil[dm]);
//
//		
//	
//	microbe.diff = microbe.Diffusivity(soil.pcfldcap/100, soil.pctpor/100, soil.vsm[dm]);
//	microbe.rn2 = microbe.Rn2(microbe.diff, pstate[I_NO3]*1000, microbe.rh[dm]*1000/days[dm], soil.pctp[dm]);
//	microbe.rnox = microbe.Rnox(microbe.diff);
//	
//	microbe.noxpi = 1.0; // when using the monthly data input
//	microbe.nflux(microbe.Fnit[dm], microbe.Fdenit[dm], microbe.rn2, microbe.rnox, microbe.noxpi);
	
//	microbe.N2[dm] *=  0.001;
//	microbe.N2O[dm] *=  0.001;
//	microbe.NOx[dm] *=  0.001;
//	 
//	microbe.Fnit[dm] *=  0.001;
//	microbe.Fdenit[dm] *=  0.001;
	  soil.ngas[dm] = microbe.N2[dm] + microbe.N2O[dm]; // + microbe.NOx[dm] + microbe.NH3[dm]
	  //printf("ttem 1077: ngas = %.3f, n2 = %.3f, n2o = %.3f\n", soil.ngas[dm], microbe.N2[dm], microbe.N2O[dm]);
//	soil.ngas[dm] = 10.0;
  	//soil.ninput[dm] = 0.01; //0.075 + 0.0028*atms.prec[dm];

    if (ag.state == 0)  
    {  //printf("ttem1083: ag.state=%i", ag.state);
      soil.nlost[dm] = pstate[I_NO3] / pstate[I_SM]; //printf("ttem1084: dm=%i,avln=%.3f,sm=%.3f,nlost=%.3f\n",dm, pstate[I_AVLN],pstate[I_SM],soil.nlost[dm]);
      //soil.nlost[dm] *= ((atms.rain[dm] + soil.snowinf[dm] - atms.eet[dm])
      //                  + (soil.rootz * 1000.0)) / (soil.rootz * 1000.0);
       soil.nlost[dm] *= soil.h2oyld[dm];
      //printf("ttem1087:sm = %.3f, runoff = %.3f, nlost = %.3f\n", pstate[I_SM],soil.h2oyld[dm],soil.nlost[dm]);
//     if  (soil.nlost[dm] > 1000.0 ) soil.nlost[dm] = 0.002771;  // denugged by Q. Z. 04/Nov/2001
//     if  (soil.nlost[dm] < 0.000001 ) soil.nlost[dm] = 0.002771;  // denugged by Q. Z. 04/Nov/2001

      soil.nlost[dm] *= soil.nloss[veg.cmnt];
   	  //printf("ttem 1092: dm=%i, nlost=%.3f.nloss=%.3f\n", dm, soil.nlost[dm],soil.nloss[veg.cmnt]);		 		
      if (soil.nlost[dm] > pstate[I_AVLN] - veg.nuptake[dm] - soil.ngas[dm]+ microbe.netnmin[dm] + microbe.N2Fix[dm] + soil.ninput[dm])
      {
        soil.nlost[dm] = pstate[I_AVLN] - veg.nuptake[dm] + microbe.netnmin[dm] + microbe.N2Fix[dm] - soil.ngas[dm] + soil.ninput[dm];
      //printf("ttem1096: nlost=%.3f, availn=%.3f,nuptake=%.3f, netnmin=%.3f,nfix=%.3f, ngas=%.3f, ninput=%.3f\n",soil.nlost[dm], pstate[I_AVLN], veg.nuptake[dm], microbe.netnmin[dm],  microbe.N2Fix[dm], soil.ngas[dm],  soil.ninput[dm]);         
      } 
      if (soil.nlost[dm] < 0.0)
      {
        soil.nlost[dm]  = 0.0;
        microbe.netnmin[dm] = soil.nlost[dm] + veg.nuptake[dm] - soil.ninput[dm] + soil.ngas[dm] - pstate[I_AVLN] - microbe.N2Fix[dm];
        //printf("ttem1101: dm=%i, netnmin=%.3f, nlost=%.3f, nuptake=%.3f, ninput=%.3f, ngas=%.3f,avail=%.3f,nfix=%.3f\n",dm,microbe.netnmin[dm],soil.nlost[dm], veg.nuptake[dm],soil.ninput[dm],soil.ngas[dm],pstate[I_AVLN],microbe.N2Fix[dm]);
//				 soil.ninput[dm] = soil.nlost[dm] + veg.nuptake[dm] - microbe.netnmin[dm] + soil.ngas[dm]
//                              - pstate[I_AVLN];
      }
      
      // to determint the n lost as NH4 and NO3	
      if (pstate[I_NH4] + pstate[I_NO3] == 0)	{
      	rnitrate = 0.5;
      	rammonium = 0.5;
      }	
			else{
				rammonium = pstate[I_NH4] / (pstate[I_NH4] + pstate[I_NO3]);
				rnitrate = pstate[I_NO3] / (pstate[I_NH4] + pstate[I_NO3]);
			}
				soil.NH4lost[dm] = soil.nlost[dm] * rammonium;  //printf("ttem 1115: soil.NH4lost=%.3f,nlost=%.3f,ramonium=%.3f\n", soil.NH4lost[dm],soil.nlost[dm], rammonium);
				soil.NO3lost[dm] = soil.nlost[dm] * rnitrate; //printf("ttem 1116: soil.no3lost=%.3f,nlost=%.3f,rnitrate=%.3f\n", soil.NO3lost[dm],soil.nlost[dm],rnitrate);
							
      if (soil.NH4lost[dm] > pstate[I_NH4] - veg.nuptake[dm] * rammonium + microbe.netnmin[dm] + microbe.N2Fix[dm] + atms.nhxdepo[dm] - microbe.Fnit[dm])
      {
        soil.NH4lost[dm] = pstate[I_NH4] - veg.nuptake[dm] * rammonium + microbe.netnmin[dm] + microbe.N2Fix[dm] + atms.nhxdepo[dm] - microbe.Fnit[dm];
      }
      if (soil.NH4lost[dm] < 0.0)
      {
      	soil.NH4lost[dm] = 0.0;
//      	microbe.Fnit[dm] = pstate[I_NH4] - veg.nuptake[dm] * rammonium + microbe.netnmin[dm]
//        									 + 0.5 * soil.ninput[dm] - soil.NH4lost[dm];
        microbe.netnmin[dm] = soil.nlost[dm] + veg.nuptake[dm] * rammonium - atms.nhxdepo[dm]	+ microbe.Fnit[dm] - pstate[I_NH4] - microbe.N2Fix[dm]; 
      //printf("ttem1129: dm=%i, netnmin=%.3f, nlost=%.3f, nuptake=%.3f, ninput=%.3f, fnit=%.3f,nh4=%.3f,nfix=%.3f\n",dm,microbe.netnmin[dm],soil.nlost[dm], veg.nuptake[dm],soil.ninput[dm],microbe.Fnit[dm],pstate[I_NH4],microbe.N2Fix[dm]);                               
      }
      if (soil.NO3lost[dm] > pstate[I_NO3] - veg.nuptake[dm] * rnitrate + microbe.Fnit[dm] + atms.noydepo[dm] - soil.ngas[dm])
      {
        soil.NO3lost[dm] = pstate[I_NO3] - veg.nuptake[dm] * rnitrate	+ microbe.Fnit[dm] + atms.noydepo[dm] - soil.ngas[dm];
      } 
      if (soil.NO3lost[dm] < 0.0)
      {
      	soil.NO3lost[dm] = 0.0;
      	microbe.Fnit[dm] = veg.nuptake[dm] * rnitrate + soil.ngas[dm] + soil.NO3lost[dm] - atms.noydepo[dm] - pstate[I_NO3];
      }
     }
    else
    {
      soil.ninput[dm] = ag.nretent[dm];
      soil.nlost[dm] = 0.0;

      // Condition added to limit addition of agfertn to only
      //   periods of crop growth
      if (ag.npp[dm].nitrogen > 0.0)
      {
        if (pstate[I_AVLN] + soil.ninput[dm] + microbe.netnmin[dm] + microbe.N2Fix[dm] - soil.nlost[dm] < ag.npp[dm].nitrogen)
        {
   	    ag.fertn[dm] = ag.npp[dm].nitrogen + soil.nlost[dm] - pstate[I_AVLN] - soil.ninput[dm] - microbe.netnmin[dm] - microbe.N2Fix[dm];
	      if (ag.fertn[dm] < 0.0)
        {
	      ag.fertn[dm] = 0.0;
	      microbe.netnmin[dm] = soil.nlost[dm] + ag.npp[dm].nitrogen - soil.ninput[dm] - pstate[I_AVLN] - microbe.N2Fix[dm];
	      }
	      soil.ninput[dm] += ag.fertn[dm];
        }
      }
    }
  }
  else
  {
    soil.nlost[dm] = soil.ninput[dm] - veg.nuptake[dm] - ag.npp[dm].nitrogen + microbe.netnmin[dm] + microbe.N2Fix[dm] ;
  }

// Describe monthly changes to carbon pools and fluxes for ODE state variables
// (i.e., pdstate)

  // Carbon pools in natural ecosystems
 //determine the DOC
	soil.docmb(pstate[I_SOLC], density);
	pdstate[I_DOC] = soil.docb * pstate[I_SOLC] - soil.docm * pstate[I_DOC];
	//printf("ttem 1176: pdDOC =%.3f soil.docb=%.3f soil.docm=%.3f psSOLC=%.3f psDOC=%.3f\n",pdstate[I_DOC],soil.docb,soil.docm,pstate[I_SOLC],pstate[I_DOC]);
  pdstate[I_VEGC] = veg.gpp[dm] - veg.gpr[dm] - veg.ltrfal[dm].carbon;  
  //printf("ttem1186: pdvegc=%.3f,gpp=%.3f,gpr=%.3f,ltrfalc=%.3f\n", pdstate[I_VEGC],veg.gpp[dm],veg.gpr[dm],veg.ltrfal[dm].carbon);  
  pdstate[I_SOLC] = veg.ltrfal[dm].carbon + ag.ltrfal[dm].carbon
                    + ag.slash[dm].carbon - ag.sconvrtflx[dm].carbon
                    - microbe.rh[dm] ; //  YYuan 122922 add n from permafrost thaw  + atms.frontddiff*14500; 
  //printf("ttem 1180: PDsolc = %.3f, ltrc = %.3f, agc = %.3f\n", pdstate[I_SOLC], veg.ltrfal[dm].carbon, ag.ltrfal[dm].carbon); 
  //printf("ttem 1181: slashc = %.3f, convc = %.3f, rh = %.3f\n", ag.slash[dm].carbon, ag.sconvrtflx[dm].carbon, microbe.rh[dm]);
  pstate[I_VEGCAG] = veg.rsratio(veg.cmnt,atms.tair[dm], pstate[I_VEGC]); 
  pstate[I_VEGCBG] = pstate[I_VEGC]- pstate[I_VEGCAG] ; 
  //printf ("ttem 1168: veg=%.i, tair=%.3f, pvegc=%.3f, pvegca=%.3f\n", veg.cmnt,atms.tair[dm], pstate[I_VEGC],pstate[I_VEGCA]);
  // Nitrogen pools in natural ecosystems
  pdstate[I_STRN] = veg.suptake[dm] - veg.ltrfal[dm].nitrogen - veg.nresorb[dm]
                    + veg.nmobil[dm];
  //printf("ttem 1187 : pdSTRN = %.3f, suptake = %.3f, ltrn = %.3f, nresorb = %.3f, nmobil = %.3f\n", pdstate[I_STRN], veg.suptake[dm], veg.ltrfal[dm].nitrogen, veg.nresorb[dm], veg.nmobil[dm]);
  pdstate[I_SOLN] = veg.ltrfal[dm].nitrogen + ag.ltrfal[dm].nitrogen
                    + ag.slash[dm].nitrogen - ag.sconvrtflx[dm].nitrogen
                    - ag.nsretent[dm] - microbe.netnmin[dm];  //YYuan 122922 add n from permafrost thaw //  + (atms.frontddiff)*1000
  //printf("ttem 1192: dm=%i,netnmin=%.3f,lttern=%.3f,pdsoilorgn=%.3f\n", dm, microbe.netnmin[dm], veg.ltrfal[dm].nitrogen, pdstate[I_SOLN]);
          
 // for NH4 and NO3 pools, soil N input is divided equally into NH4 and NO3. added by YYL             
  pdstate[I_NH4] =  atms.nhxdepo[dm] + microbe.netnmin[dm] + microbe.N2Fix[dm] - microbe.Fnit[dm] 
  							- soil.NH4lost[dm] - veg.nuptake[dm] * rammonium;
  									 
  pdstate[I_NO3] = atms.noydepo[dm] + microbe.Fnit[dm] - soil.ngas[dm]
   						  - soil.NO3lost[dm] - veg.nuptake[dm] * rnitrate;
   									 
  pdstate[I_AVLN] = soil.ninput[dm] - soil.nlost[dm] + microbe.netnmin[dm] + microbe.N2Fix[dm] - soil.ngas[dm]
                    - veg.nuptake[dm] - ag.npp[dm].nitrogen;
 	//if (pdstate[I_NH4]<0)  pdstate[I_NH4]=0;
	//if (pdstate[I_NO3]<0)  pdstate[I_NO3]=0;                  
  //printf("ttem 1204: pdAVLN = %.5f, ninput = %.3f, nlost = %.3f, netnmin = %.3f, nfix=%.3f, ngas = %.3f, nuptake = %.3f, nppn = %.3f\n", pdstate[I_AVLN], soil.ninput[dm], soil.nlost[dm], microbe.netnmin[dm],  microbe.N2Fix[dm], soil.ngas[dm], veg.nuptake[dm], ag.npp[dm].nitrogen);
  //printf("ttem 1170:I_NH4=%.3f, I_NO3=%.3f\n ",pdstate[I_NH4],pdstate[I_NO3]);
  pdstate[I_STON] = veg.luptake[dm] + veg.nresorb[dm] - veg.nmobil[dm];
  if (pdstate[I_STON] <0) {pdstate[I_STON]=0;}

  // Human product pools

  pdstate[I_AGPRDC] = 0.0;
  pdstate[I_AGPRDN] = 0.0;
  pdstate[I_PROD10C] = 0.0;
  pdstate[I_PROD10N] = 0.0;
  pdstate[I_PROD100C] = 0.0;
  pdstate[I_PROD100N] = 0.0;
  pdstate[I_TOTPRDC] = 0.0;
  pdstate[I_TOTPRDN] = 0.0;
  pdstate[I_TOTEC] = 0.0;
  pdstate[I_TOTGC] = 0.0;

  // Water pools

  pdstate[I_AVLW] = soil.snowinf[dm] + atms.rain[dm] - atms.eet[dm]
                    - soil.rperc[dm] - soil.sperc[dm];
  pdstate[I_RGRW] = soil.rperc[dm] - soil.rrun[dm];
  pdstate[I_SNWPCK] = atms.snowfall[dm] - soil.snowinf[dm];
  pdstate[I_SGRW] = soil.sperc[dm] - soil.srun[dm];
  pdstate[I_SM] = pdstate[I_AVLW];
  pdstate[I_PCTP] = 100.0 * pdstate[I_SM]/soil.totpor;
  pdstate[I_VSM] = pdstate[I_SM]/(soil.rootz*1000.0);
  if (pstate[I_VSM]+pdstate[I_VSM] <= 0.0) {
    pdstate[I_VSM] = 0.001 - pstate[I_VSM];
  }

  // Carbon fluxes in natural ecosystems

  pdstate[I_INGPP] = veg.ingpp[dm];  //printf("ttem 1248: pdINGPP = %.3f,psINGPP = %.3f,ingpp=%.3f\n", pdstate[I_INGPP], pstate[I_INGPP],veg.ingpp[dm]);
  pdstate[I_GPP] = veg.gpp[dm];      //printf("ttem 1249: pdGPP = %.3f,psegpp=%.3f,gpp=%.3f\n", pdstate[I_GPP],pstate[I_GPP],veg.gpp[dm]);
  pdstate[I_INNPP] = veg.innpp[dm];
  pdstate[I_NPP] = veg.npp[dm];      //printf("ttem 1252: pdNPP = %.3f,psnpp=%.3f,npp=%.3f\n", pdstate[I_NPP],pstate[I_NPP],veg.npp[dm]);
  pdstate[I_GPR] = veg.gpr[dm];      //printf("ttem 1253: pdGPR = %.3f,psgpr=%.3f,gpr=%.3f\n", pdstate[I_GPR],pstate[I_GPR],veg.gpr[dm]);
  pdstate[I_RVMNT] = veg.rm[dm];
  pdstate[I_RVGRW] = veg.rg[dm];
  pdstate[I_LTRC] = veg.ltrfal[dm].carbon; //printf("ttem 1256: ltrc = %.3f\n", veg.ltrfal[dm].carbon);    
  pdstate[I_RH] = microbe.rh[dm];
  pdstate[I_NEP] = nep[dm];

  // Nitrogen fluxes in natural ecosystems

  pdstate[I_NINP] = soil.ninput[dm];  //printf("ttem1255:pdNINP=%.3f,ninput=%.3f\n",pdstate[I_NINP],soil.ninput[dm]);
  pdstate[I_INNUP] = veg.inuptake;
  pdstate[I_VNUP] = veg.nuptake[dm];   //printf("ttem1253:pdIVNUP=%.3f,vegnupt=%.3f\n",pdstate[I_VNUP] * 1000,veg.nuptake[dm]*1000);
  pdstate[I_VSUP] = veg.suptake[dm];
  pdstate[I_VLUP] = veg.luptake[dm];
  pdstate[I_VNMBL] = veg.nmobil[dm];
  pdstate[I_VNRSRB] = veg.nresorb[dm];
  pdstate[I_LTRN] = veg.ltrfal[dm].nitrogen; //printf("ttem 1257: ltrn = %.3f\n", veg.ltrfal[dm].nitrogen);
  pdstate[I_MNUP] = microbe.nuptake[dm];
  pdstate[I_NMIN] = microbe.netnmin[dm]; //printf("ttem 1254: nmin = %.3f\n", microbe.netnmin[dm]);
  pdstate[I_FNIT] = microbe.Fnit[dm];
  pdstate[I_FDENIT] = microbe.Fdenit[dm]; //printf("ttem 1267: denit = %.3f\n", microbe.Fdenit[dm]);
  pdstate[I_NLST] = soil.nlost[dm];
  pdstate[I_N2] = microbe.N2[dm];
  pdstate[I_N2O] = microbe.N2O[dm];
  pdstate[I_NOX] = microbe.NOx[dm];
  pdstate[I_NGAS] = soil.ngas[dm];
  pdstate[I_N2F] = microbe.N2F[dm]; //YY
  pdstate[I_N2OF] = microbe.N2OF[dm];// YY
  pdstate[I_N2OAIR] = microbe.N2OAir[dm]; // YY
  pdstate[I_N2AIR] = microbe.N2Air[dm]; // YY 
  pdstate[I_N2OUPT] = microbe.N2OUpt[dm];// YY 
  pdstate[I_N2ON] = microbe.N2ON[dm];// YY 
  pdstate[I_N2ODN] = microbe.N2ODN[dm];// YY 
  pdstate[I_N2FIX] = microbe.N2Fix[dm];// YY 
  pdstate[I_N2FIXS] = microbe.N2Fixs[dm];// YY 
  pdstate[I_N2FIXFL] = microbe.N2Fixfl[dm]; //printf("ttem 1277: n2fixflmon = %.3f\n", microbe.N2Fixfl[dm]);
  pdstate[I_NH3] = microbe.NH3[dm]; // YY 
  // Water fluxes

  pdstate[I_RAIN] = atms.rain[dm];
  pdstate[I_RPERC] = soil.rperc[dm];
  pdstate[I_RRUN] = soil.rrun[dm];
  pdstate[I_SNWFAL] = atms.snowfall[dm];
  pdstate[I_SNWINF] = soil.snowinf[dm];
  pdstate[I_SPERC] = soil.sperc[dm];
  pdstate[I_SRUN] = soil.srun[dm];
  pdstate[I_PET] = atms.pet[dm];
  pdstate[I_EET] = atms.eet[dm];
  pdstate[I_WYLD] = soil.rrun[dm] + soil.srun[dm];

       // for soil temperature
  pdstate[I_TSOIL] = atms.tsoil[dm];
  pdstate[I_DST5] = atms.dst5[dm];
  pdstate[I_DST10] = atms.dst10[dm];
  pdstate[I_DST20] = atms.dst20[dm];
  pdstate[I_DST50] = atms.dst50[dm];
  pdstate[I_DST100] = atms.dst100[dm];
  pdstate[I_DST200] = atms.dst200[dm];
  pdstate[I_FRONTD] = atms.frontd[dm];
  pdstate[I_THAWBE] = atms.thawbe[dm];
  pdstate[I_THAWEND] = atms.thawend[dm];
// end of ...


  // Phenology

  pdstate[I_UNRMLF] = veg.unnormleaf[dm]; //printf("ttem1313: dm=%i,pdUNRMLF=%.3f, unnormleaf=%.3f,psUNRMLF=%.3f\n",dm, pdstate[I_UNRMLF],veg.unnormleaf[dm],pstate[I_UNRMLF]);
  pdstate[I_LEAF] = veg.leaf[dm];

  // Carbon and nitrogen fluxes from agricultural conversion

  pdstate[I_CNVRTC] = ag.convrtflx[dm].carbon;
  pdstate[I_CNVRTN] = ag.convrtflx[dm].nitrogen;
  pdstate[I_SCNVRTC] = ag.sconvrtflx[dm].carbon;
  pdstate[I_SCNVRTN] = ag.sconvrtflx[dm].nitrogen;
  pdstate[I_NVRTNT] = ag.nvretent[dm];
  pdstate[I_NSRTNT] = ag.nsretent[dm];
  pdstate[I_NRETNT] = ag.nretent[dm];
  pdstate[I_SLASHC] = ag.slash[dm].carbon;
  pdstate[I_SLASHN] = ag.slash[dm].nitrogen;
  pdstate[I_PRDF10C] = ag.formPROD10.carbon / (double) CYCLE;
  pdstate[I_PRDF10N] = ag.formPROD10.nitrogen / (double) CYCLE;
  pdstate[I_PRDF100C] = ag.formPROD100.carbon / (double) CYCLE;
  pdstate[I_PRDF100N] = ag.formPROD100.nitrogen / (double) CYCLE;

  // Carbon and nitrogen fluxes in agricultural ecosystems

  pdstate[I_AGNPPC] = ag.npp[dm].carbon;
  pdstate[I_AGNPPN] = ag.npp[dm].nitrogen;
  pdstate[I_AGFPRDC] = ag.npp[dm].carbon - ag.ltrfal[dm].carbon;
  pdstate[I_AGFPRDN] = ag.npp[dm].nitrogen - ag.ltrfal[dm].nitrogen;
  pdstate[I_AGLTRC] = ag.ltrfal[dm].carbon;
  pdstate[I_AGLTRN] = ag.ltrfal[dm].nitrogen;
  pdstate[I_AGFRTN] = ag.fertn[dm];

  // Carbon and nitrogen fluxes from human product pools

  pdstate[I_TOTFPRDC] = ag.formTOTPROD.carbon / (double) CYCLE;
  pdstate[I_TOTFPRDN] = ag.formTOTPROD.nitrogen / (double) CYCLE;
  pdstate[I_AGPRDFC] = ag.PROD1decay.carbon / (double) CYCLE;
  pdstate[I_AGPRDFN] = ag.PROD1decay.nitrogen / (double) CYCLE;
  pdstate[I_PRD10FC] = ag.PROD10decay.carbon / (double) CYCLE;
  pdstate[I_PRD10FN] = ag.PROD10decay.nitrogen / (double) CYCLE;
  pdstate[I_PRD100FC] = ag.PROD100decay.carbon / (double) CYCLE;
  pdstate[I_PRD100FN] = ag.PROD100decay.nitrogen / (double) CYCLE;
  pdstate[I_TOTPRDFC] = ag.TOTPRODdecay.carbon / (double) CYCLE;
  pdstate[I_TOTPRDFN] = ag.TOTPRODdecay.nitrogen / (double) CYCLE;

  // Integrated carbon fluxes

  pdstate[I_TOTNPP] = veg.npp[dm] + ag.npp[dm].carbon;
  pdstate[I_CFLX] = cflux[dm];

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TTEM::deltaxclm(const int& dcmnt, const double& pcfldcap, const int& dm)
{

  double vfc;


  vfc = pcfldcap * 0.01;


/* gv: effect of moisture on primary productivity */

  if (moistlim == 0)
  {
    gv = 1.0;
  }
  else
  {
    if (atms.eet[dm]/atms.pet[dm] <= 0.1)
    {
      gv = (-10.0 * pow((atms.eet[dm]/atms.pet[dm]),2.0))
           + (2.9 * (atms.eet[dm]/atms.pet[dm]));
      if (gv < 0.0) { gv = 0.0; }
    }
    else
    {
      gv = 0.1 + (0.9 * atms.eet[dm] / atms.pet[dm]);
    }
  }


/* ksoil: effect of soil moisture on nitrogen uptake by plants
			 and microbes */

  if (moistlim == 0) { ksoil = pow(vfc,3.0); }
  else { ksoil = pow(y[I_VSM],3.0); } //3.0 orignal
  //printf("ttem1363: ksoil=%.3f, moistlim=%.3f,vfc=%3f,y[I_VSM]=%.3f\n", ksoil, moistlim, vfc, y[I_VSM]);

/* rhmoist: effect of moisture on decomposition */

  if (moistlim == 0)
  {
    rhmoist = (vfc - microbe.moistmin[dcmnt])
              * (vfc - microbe.moistmax[dcmnt]);  //printf("ttem1406: rhmoist=%.3f, vfc=%.3f,moistmin=%.3f,moistman=%.3f\n",rhmoist, vfc, microbe.moistmin[dcmnt], microbe.moistmax[dcmnt]);
    rhmoist /= rhmoist - pow((vfc - microbe.moistopt[dcmnt]),2.0);  //printf("ttem1407: rhmoist=%.3f, vfc=%.3f,moistopt=%.3f\n",rhmoist, vfc, microbe.moistopt[dcmnt]);
  }
  else
  {
    rhmoist = (y[I_VSM] - microbe.moistmin[dcmnt])
              * (y[I_VSM] - microbe.moistmax[dcmnt]);  //printf("ttem1412: rhmoist=%.3f,vfc=%.3f,moistmin=%.3f,moistman=%.3f\n",rhmoist,y[I_VSM], microbe.moistmin[dcmnt], microbe.moistmax[dcmnt]);
    rhmoist /= rhmoist - pow((y[I_VSM] - microbe.moistopt[dcmnt]),2.0);//printf("ttem1413: rhmoist=%.3f, yvsm=%.3f,moistopt=%.3f\n",rhmoist,y[I_VSM], microbe.moistopt[dcmnt]);
  }
  if (rhmoist < 0.0) { rhmoist = 0.0; }

};

/* *************************************************************
************************************************************* */



/* *************************************************************
************************************************************* */

int TTEM::ecdqc(const int& dcmnt)
{

  int qc = ACCEPT;

  if (vegca[dcmnt] <= -999.99) { return qc = REJECT; }
  if (vegcb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (strna[dcmnt] <= -999.99) { return qc = REJECT; }
  if (strnb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (solca[dcmnt] <= -999.99) { return qc = REJECT; }
  if (solcb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (solna[dcmnt] <= -999.99) { return qc = REJECT; }
  if (solnb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (avlna[dcmnt] <= -999.99) { return qc = REJECT; }
  if (avlnb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (stona[dcmnt] <= -999.99) { return qc = REJECT; }
  if (stonb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.unleaf12[dcmnt] <= -9.99) { return qc = REJECT; }
  if (veg.prvleafmx[dcmnt] <= -9.99) { return qc = REJECT; }
  if (veg.cmaxcut[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.cmax1a[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.cmax1b[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.cmax2a[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.cmax2b[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.cfall[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.kra[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.krb[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.kda[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.kdb[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.lcclnc[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.propftos[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.nmaxcut[dcmnt] <= -99.99) { return qc = REJECT; }
  // DWK changed missing values from -99.99 to -999.99 for "nmax"s
  if (veg.nmax1a[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.nmax1b[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.nmax2a[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.nmax2b[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.nfall[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.nupa[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.nupb[dcmnt] <= -99.99) { return qc = REJECT; }
  if (soil.nloss[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.nfixpar[dcmnt] <= -99.9) { return qc = REJECT; }
  if (veg.cneven[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.cnmin[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.c2na[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.c2nb[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.c2nmin[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.cnsoil[dcmnt] <= -999.99) { return qc = REJECT; }

  return qc;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void TTEM::ECDsetELMNTstate(const int& dcmnt, const double& psiplusc)
{

  int dyr;
  int dm;

  for (dm = 0; dm < CYCLE; dm++)
  {
    y[I_AVLW] = soil.avlh2o[dm] = soil.awcapmm;
    y[I_RGRW] = soil.rgrndh2o[dm] = 0.0;
    y[I_SNWPCK] = soil.snowpack[dm] = 0.0;
    y[I_SGRW] = soil.sgrndh2o[dm] = 0.0;
    y[I_SM] = soil.moist[dm] = soil.awcapmm + soil.wiltpt;
    y[I_PCTP] = soil.pctp[dm] = 100.0 * soil.moist[dm] / soil.totpor;
    soil.vsm[dm] = soil.moist[dm] / (soil.rootz * 1000.0);
    if (soil.vsm[dm] <= 0.0)
    {
      soil.vsm[dm] = 0.001;
    }
    y[I_VSM] = soil.vsm[dm];

    veg.plant[dm].carbon = vegca[dcmnt] * psiplusc + vegcb[dcmnt];
    y[I_VEGC] = veg.plant[dm].carbon; //printf("ttem1520: dm=%i, yvegc=%.3f, vegplantcarbon=%.3f\n",dm, y[I_VEGC], veg.plant[dm].carbon);
    
    veg.planta[dm].carbon = 0.74 * veg.plant[dm].carbon;  //Mokany et al.,2006; Reich et al 2014
    y[I_VEGCAG] = veg.planta[dm].carbon; //printf("ttem1494: dm=%i, yvegcb=%.3f, vegplantbcarbon=%.3f\n",dm, y[I_VEGCBG], veg.plantb[dm].carbon);

    veg.plantb[dm].carbon = 0.26 * veg.plant[dm].carbon;  //Mokany et al.,2006; Reich et al 2014
    y[I_VEGCBG] = veg.plantb[dm].carbon; //printf("ttem1494: dm=%i, yvegcb=%.3f, vegplantbcarbon=%.3f\n",dm, y[I_VEGCBG], veg.plantb[dm].carbon);
    
    veg.strctrl[dm].nitrogen = strna[dcmnt] * psiplusc + strnb[dcmnt];// printf("ttem1516: dm=%i,dcmnt=%i,vegstrctrln=%.3f,strna=%.3f,psiplusc=%.3f,strnb=%.3f\n",dm,dcmnt,veg.strctrl[dm].nitrogen,strna[dcmnt],psiplusc,strnb[dcmnt]);
    y[I_STRN] = veg.strctrl[dm].nitrogen; //printf("ttem1517: dm=%i, ystern=%.3f,vegstrctrln=%.3f\n", dm,y[I_STRN], veg.strctrl[dm].nitrogen);

    soil.org[dm].carbon = solca[dcmnt] * psiplusc + solcb[dcmnt];
    y[I_SOLC] = soil.org[dm].carbon;

		soil.doc[dm] = 0.1 * soil.org[dm].carbon;
		y[I_DOC] = soil.doc[dm];
		
    soil.org[dm].nitrogen = solna[dcmnt] * psiplusc + solnb[dcmnt]; //printf("ttem1525: dm=i.soln=%.3f, solna=%.3f, psiplusc=%.3f, solnb=%.3f\n",dm,soil.org[dm].nitrogen,solna[dcmnt],psiplusc,solnb[dcmnt]);
    y[I_SOLN] = soil.org[dm].nitrogen; 

    soil.availn[dm] = avlna[dcmnt] * psiplusc + avlnb[dcmnt]; //printf("ttem1529: soil.availn=%.3f,psiplusc=%.3f\n",soil.availn[dm],psiplusc);
    y[I_AVLN] = soil.availn[dm];  //printf("ttem1530:  yravaln=%.3f\n",y[I_AVLN]); 
    
    microbe.NH4[dm] = 0.5 * soil.availn[dm];
    y[I_NH4] = microbe.NH4[dm];
    
    microbe.NO3[dm] = 0.5 * soil.availn[dm];
    y[I_NO3] = microbe.NO3[dm];

    veg.labile[dm].nitrogen = stona[dcmnt] * psiplusc + stonb[dcmnt];
    y[I_STON] = veg.labile[dm].nitrogen;

    y[I_UNRMLF] = veg.unnormleaf[dm] = veg.unleaf12[dcmnt]; //printf("ttem 1546: yUNRMLF=%.3f,unnormleaf=%.3f,unleaf12=%.3f\n",y[I_UNRMLF],veg.unnormleaf[dm],veg.unleaf12[dcmnt]);

    veg.plant[dm].nitrogen = 0.0;

    veg.alleaf = veg.leafmxc[dcmnt]/(1.0 + veg.kleafc[dcmnt]
                 * exp(veg.cov[dcmnt]*y[I_VEGC]));
    veg.lai[dm] = veg.sla[dcmnt] * veg.alleaf;
    veg.foliage = veg.alleaf / veg.leafmxc[dcmnt];
    veg.fpc[dm] = 1.0;

    y[I_PROD10C] = ag.PROD10.carbon = 0.0;
    y[I_PROD10N] = ag.PROD10.nitrogen = 0.0;
    y[I_PROD100C] = ag.PROD100.carbon = 0.0;
    y[I_PROD100N] = ag.PROD100.nitrogen = 0.0;
    y[I_AGPRDC] = ag.PROD1.carbon = 0.0;
    y[I_AGPRDN] = ag.PROD1.nitrogen = 0.0;
    y[I_TOTPRDC] = ag.TOTPROD.carbon = 0.0;
    y[I_TOTPRDN] = ag.TOTPROD.nitrogen = 0.0;
    y[I_TOTEC] = totalc[dm] = veg.plant[dm].carbon + soil.org[dm].carbon;
    y[I_TOTGC] = totalc[dm];
  }

  for (dyr = 0; dyr < 10; dyr++)
  {
    ag.initPROD10[dyr].carbon = 0.0;
    ag.initPROD10[dyr].nitrogen = 0.0;
  }

  for (dyr = 0; dyr < 100; dyr++)
  {
    ag.initPROD100[dyr].carbon = 0.0;
    ag.initPROD100[dyr].nitrogen = 0.0;
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

int TTEM::equilibrium(const int& itype, double& tol)
{

  int dyr = 0;

  setPrevState(prevy,y);

  totyr = 0;
  endeq = 0;
  intflag = 0;
  initFlag = 0;
  while ((dyr < runsize) && (endeq < 2))
  {
    endeq = stepyr(dyr,itype,intflag,tol);
    ++dyr;
    ++totyr;

 // Reset product fluxes and pools to zero

    ag.resetPROD();

// Check to see if steady state conditions have been reached.

    if (dyr >= strteq && endeq == 0)
    {
      if (nfeed == 0 && rheqflag == 0
         && (wtol >= fabs(atms.yrrain + atms.yrsnowfall - atms.yreet
         - soil.yrh2oyld))
         && ((ctol >= fabs(veg.yrnpp - veg.yrltrc))
         || ( 0.001 >= fabs(veg.yrnpp - veg.yrltrc))))
      {
	     endeq = 1;
      }
      if (nfeed == 0 && rheqflag == 1 && (wtol >= fabs(atms.yrrain
         + atms.yrsnowfall - atms.yreet - soil.yrh2oyld))
         && (ctol >= fabs(yrnep))
         && (ctol >= fabs(veg.yrnpp - veg.yrltrc))
         && (ctol >= fabs(veg.yrltrc - microbe.yrrh)))
      {
	     endeq = 1;
      }
      if (nfeed == 1 && rheqflag == 1 && (wtol >= fabs(atms.yrrain
         + atms.yrsnowfall - atms.yreet - soil.yrh2oyld))
         && (ntol >= fabs(soil.yrnin - soil.yrnlost - soil.yrngas))
         && (ntol >= fabs(veg.yrnup - veg.yrltrn))
         && (ntol >= fabs(veg.yrnup - microbe.yrnmin))
         && (ntol >= fabs(veg.yrltrn - microbe.yrnmin))
         && (ctol >= fabs(yrnep)) && (ctol >= fabs(veg.yrnpp - veg.yrltrc))
         && (ctol >= fabs(veg.yrltrc - microbe.yrrh)))
      {
        endeq = 1;
      }
    }
  }

  if (endeq == 2)
  {
    nattempt = maxnrun;
    initFlag = 1;
  }

  if (dyr >= runsize && endeq < 2) { ++nattempt;
  initFlag = 1;} 

  return nattempt;

};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TTEM::getenviron(const int& dm)
{

  // Determine monthly potential evapotranspiration

  atms.pet[dm] = atms.petjh(atms.nirr[dm], atms.tair[dm], dm);
  if (atms.pet[dm] <= 0.0) { atms.pet[dm] = 0.001; }

  // Determine if monthly precipitation occurs as rain or snow

  atms.precsplt(atms.prec[dm],atms.tair[dm],atms.rain[dm],atms.snowfall[dm]);

  // Determine previous two month's air temperatures for following
  //   snowmelt calculations

  switch (dm)
  {
    case 0:  atms.prevtair = atms.tair[CYCLE-1];
             atms.prev2tair = atms.tair[CYCLE-2];
             break;
    case 1:  atms.prevtair = atms.tair[0];
	     atms.prev2tair = atms.tair[CYCLE-1];
	     break;
    default: atms.prevtair = atms.tair[dm-1];
	     atms.prev2tair = atms.tair[dm-2];
	     break;
  }

  // Determine contribution of snowmelt to soil moisture

  soil.snowinf[dm] = soil.snowmelt(elev, atms.tair[dm],
                                   atms.prevtair, y[I_SNWPCK]);

  monthxclm(veg.cmnt, veg.topt, dm);
  
  //get daily soil temperature and moisture; 
  	if (days[dm]==0) days[dm]=30;
  		 // for N cycle
 //****************************************************************************/
    	int m;
//     double CDM;
//     CDM = 0.0;
//     for (m = 0; m < CYCLE; m++)
//     CDM = CDM + (10.0 - atms.tair[m]);

   	int pday;
   	switch (dm){
   		case 0: pday = 1;  break;
   		case 1: pday = 32; break;
   		case 2: pday = 60; break;
   		case 3: pday = 91; break;
   		case 4: pday = 121; break;
   		case 5: pday = 152; break;
   		case 6: pday = 182; break;
   		case 7: pday = 213; break;
   		case 8: pday = 244; break;	
   		case 9: pday = 274; break;
   		case 10: pday = 305; break;
   		case 11: pday = 335; break;
   	}

   for (m = 0; m < days[dm]; m++){
     		//use HM daily step added by YYL
   	   // convert the monthly data to daily

   	    daylai[dm][m]=veg.lai[dm];
   	    daynirr[dm][m]=atms.nirr[dm];
   	   
		    daysnowpack[dm][m]=soil.snowpack[dm];
		
   	  // Determine if daily precipitation occurs as rain or snow   
 //  	 atms.precsplt(atms.prec[dm][m],atms.tair[dm][m],atms.rain[dm][m],atms.snowfall[dm][m]);
			atms.precsplt(dayprec[dm][m],daytair[dm][m],dayrain[dm][m],daysnowfall[dm][m]);
 //  if (atms.GISlai[dm][m] == 0.0) atms.GISlai[dm][m] = 0.01;  // for drbugging Q. Z.
 //  if (atms.GISlai[dm][m] < 0.0) atms.GISlai[dm][m] = 0.5;  // to deal with non-data grid cells
 //  veg.lai[dm][m] = atms.GISlai[dm][m]; // for calibration version

    hyd.vaporpressure[dm][m] = dayvap[dm][m]/10.0;  //use CRU data sets, converse to kpa
 //  hyd.vpd(hyd.vaporpressure[dm][m], atms.tair[dm][m]);
		hyd.vpd(hyd.vaporpressure[dm][m], daytair[dm][m]);

 //  hyd.rainthrough[dm][m] = hyd.intercept(atms.rain[dm][m], hyd.inter_coef[veg.cmnt],dm, m);
   hyd.rainthrough[dm][m] = hyd.intercept(dayrain[dm][m], hyd.inter_coef[veg.cmnt],dm, m);

   soil.avlh2o1[dm][m] = 2000.0 * (hydm.h2osoi[0] + hydm.h2osoi[1] + hydm.h2osoi[2] + hydm.h2osoi[3] + hydm.h2osoi[4]
                          + hydm.h2osoi[5])/6.0;

//   soil.avlh2o1[dm][m] =  soil.awcapmm + hyd.rainthrough[dm][m];
//   soil.avlh2o2[dm][m] =  soil.awcapmm + hyd.rainthrough[dm][m];

   hyd.canopy_water[dm][m] = hyd.prec_to_canopy[dm][m];

   hyd.drad(daylai[dm][m], hyd.EXT[veg.cmnt], daynirr[dm][m], dm, m);
   //printf("ttem 1699: dm = %i, m = %i, dard = %.3f, lai = %.3f, ext = %.3f, nirr = %.3f\n", dm, m, hyd.hdrad[dm][m], daylai[dm][m], hyd.EXT[veg.cmnt], daynirr[dm][m]);
   hyd.hLWP[dm][m] =  hyd.LWP(soil.avlh2o1[dm][m], soil.fldcap, soil.pctsand, soil.pctclay, soil.pctsilt);

   hyd.hpa[dm][m] = hyd.densityA(daytair[dm][m]);
   hyd.hlv[dm][m] = hyd.latentV(daytair[dm][m]);

  //rainfall is changed to rainthrough fall due to interception
  //  atms.rain[dm][m] = hyd.rainthrough[dm][m];

   hyd.Gs[dm][m] = hyd.CanopyConductance(daylai[dm][m], hyd.gmax[veg.cmnt],
           daytair[dm][m], hyd.vpdeficit, hyd.hLWP[dm][m],
           hyd.vpd_close[veg.cmnt],hyd.vpd_open[veg.cmnt], hyd.psi_close[veg.cmnt], hyd.psi_open[veg.cmnt]);

   hyd.htrans[dm][m] = hyd.penmanmonteith_new(daytair[dm][m], hyd.hdrad[dm][m], hyd.hpa[dm][m],
             hyd.vpdeficit, hyd.hlv[dm][m], hyd.Gs[dm][m], daylai[dm][m]);
   //printf("ttem 1714: dm = %i, m = %i, htrans = %.3f, tair = %.3f, hdard = %.3f, hpa = %.3f, vpd = %.3f, hlv = %.3f, gs = %.3f, lai = %.3f\n", dm, m, hyd.htrans[dm][m], daytair[dm][m], hyd.hdrad[dm][m], hyd.hpa[dm][m], hyd.vpdeficit, hyd.hlv[dm][m], hyd.Gs[dm][m], daylai[dm][m]);
   hyd.potevap_surface[dm][m] = hyd.evaporation_layer1(daytair[dm][m],
                                 hyd.hpa[dm][m],hyd.hlv[dm][m],hyd.vpdeficit, hyd.surfacerad[dm][m]);

 // consider only the precipitation flux reaching the soil
 // check for precipitation >= potential evaporation

//   	if (hyd.rainthrough[dm][m] >= hyd.potevap_surface[dm][m])
//         	hyd.soil_evap[dm][m] = hyd.potevap_surface[dm][m];
//   	else	hyd.soil_evap[dm][m] = 0.0;

  	hyd.soil_evap[dm][m] = hyd.potevap_surface[dm][m];

 //snow intercept
    if (daysnowpack[dm][m] > 0.0) {
          hyd.sub_from_canopy[dm][m] = hyd.sublimation_canopy(daylai[dm][m],
              hyd.hdrad[dm][m], daytair[dm][m], daysnowpack[dm][m], dm, m);
          if (hyd.sub_from_canopy[dm][m] == 0.0)  hyd.sub_from_canopy[dm][m] = 0.001; // to protect
     }

    soil.snowthrough[dm][m] = hyd.snowthroughfall[dm][m];
    hyd.snowmelt[dm][m] = hyd.snowmelt_srad(hyd.hdrad[dm][m], daytair[dm][m], soil.snowthrough[dm][m]);

 // snowpack is updated to snowthroug[dm][m]
    daysnowpack[dm][m] = soil.snowthrough[dm][m] - hyd.snowmelt[dm][m];
    //printf("ttem 1744:daysnowpack=%.3f\n ", daysnowpack[dm][m]);
 // sublimation from ground snow
    hyd.snowsub[dm][m] = hyd.snowsublimation_ground(soil.snowthrough[dm][m], daytair[dm][m], hyd.surfacerad[dm][m]);

 // the estimate of actual evapotranspiration
    dayeet3[dm][m] = hyd.htrans[dm][m] + hyd.sub_from_canopy[dm][m] + hyd.soil_evap[dm][m] + hyd.snowsub[dm][m];

 // the following is the original one to estimate the potential evapotranspiration
    daypet3[dm][m] = atms.new_petjh(daynirr[dm][m], daytair[dm][m], dm);
    if (daypet3[dm][m] <= 0.0) { daypet3[dm][m] = 0.001; }

 // deal with soil hydrological dynamics
    hydm.surrunoff[dm][m] = hydm.surfacerunoff(hyd.snowmelt[dm][m],hyd.rainthrough[dm][m],hydm.h2osoi[0],
            hydm.b1[veg.cmnt], hydm.ksat1[veg.cmnt], hydm.thetasat1[veg.cmnt], hydm.presat1[veg.cmnt],
            hydm.thick1[veg.cmnt]);
   hydm.surinfl[dm][m] = hydm.infiltration(hydm.surrunoff[dm][m], hyd.snowmelt[dm][m],hyd.rainthrough[dm][m]);
  // printf("ttem 1755: dm = %i, m = %i, surinfl = %.3f, surroff = %.3f, snowmelt = %.3f, raint = %.3f\n", dm, m,  hydm.surinfl[dm][m], hydm.surrunoff[dm][m], hyd.snowmelt[dm][m],hyd.rainthrough[dm][m]);
 // calculate soil water content and sub-surafce runoff (drainage)
    hydm.qtran = hyd.htrans[dm][m];
    hydm.qseva = hyd.soil_evap[dm][m];
 //  hydm.bdrai[dm][m] = hydm.qdrai[dm][m];

 // for general soil
    hydm.main_soil(hydm.nsl, hyd.htrans[dm][m], hydm.surinfl[dm][m], hyd.soil_evap[dm][m], hydm.bdrai[dm][m], hydm.dtsoi);

//	for wetland
//			hydm.wetland_soil(hydm.nsl, hyd.htrans[dm][m], hydm.surinfl[dm][m], hyd.soil_evap[dm][m], hydm.bdrai[dm][m], hydm.dtsoi);
		
		if (hydm.h2osoi[0] > soil.totpor1) hydm.h2osoi[0] = soil.totpor1;
		if (hydm.h2osoi[1] > soil.totpor2) hydm.h2osoi[1] = soil.totpor2;					
		if (hydm.h2osoi[2] > soil.totpor3) hydm.h2osoi[2] = soil.totpor3;		
    //printf("ttem 1777: totpor1 = %.3f, totpor2 = %.3f, totpor3 = %.3f\n", soil.totpor1, soil.totpor2, soil.totpor3);
		//daypctp[dm][m] = 100 * (hydm.h2osoi[0] /(soil.totpor1 / 100.0) +  hydm.h2osoi[1] /(soil.totpor2 / 100.0) +  hydm.h2osoi[2] /(soil.totpor3 / 100.0))/3;
    daypctp[dm][m] = 100 * (hydm.h2osoi[0] /(soil.totpor1 / 100.0) +  hydm.h2osoi[1] /(soil.totpor2 / 100.0) + hydm.h2osoi[2] /(soil.totpor3 / 100.0))/3;
		//daypctp[dm][m] = 100 * ((hydm.h2osoi[0] + hydm.h2osoi[1] + hydm.h2osoi[2])/(soil.totpor1 / 100.0))/3;                                                                  
    //printf("ttem 1810: daypctp=%.3f, dm =%i, m = %i, h2o0 = %.3f, h2o1 = %.3f, h2o2 = %.3f,h2o3=%.3f\n", daypctp[dm][m], dm, m, hydm.h2osoi[0], hydm.h2osoi[1], hydm.h2osoi[2],hydm.h2osoi[3]);
    //printf("ttem 1811: rootz = %.3f, totpor1 =%.3f, totpor2=%.3f, totpor3=%.3f\n",soil.rootz, soil.totpor1, soil.totpor2, soil.totpor3);
 //end of daily HM
//***********************************************************************************************************/
 // call soil thermal subroutine
     switch (pday) {
     case 1: sthermal.airt19= tempdaily[pday]; // satisfy soil thermal model
             sthermal.airt29= tempdaily[pday];
             sthermal.airt39= (tempdaily[pday]+tempdaily[pday+1])/2.0;
             break;
     case 365:sthermal.airt19= (tempdaily[pday-1]+tempdaily[pday])/2.0; // satisfy soil thermal model
             sthermal.airt29= tempdaily[pday];
             sthermal.airt39= tempdaily[pday];
            break;
     default:sthermal.airt19= (tempdaily[pday-1]+tempdaily[pday])/2.0; // satisfy soil thermal model
             sthermal.airt29= tempdaily[pday];
             sthermal.airt39= (tempdaily[pday]+tempdaily[pday+1])/2.0;
            break;
         }

    switch (dm) {
     case 0: sthermal.hsnow19=(soil.snowpack[CYCLE-1]+soil.snowpack[0])/2.0/1000.0; // satisfy soil thermal model
             sthermal.hsnow29= soil.snowpack[0]/1000.0;
             sthermal.hsnow39= (soil.snowpack[1]+soil.snowpack[0])/2.0/1000.0;
             break;
     case 11: sthermal.hsnow19=(soil.snowpack[11]+soil.snowpack[10])/2.0/1000.0; // satisfy soil thermal model
              sthermal.hsnow29= soil.snowpack[11]/1000.0;
              sthermal.hsnow39= (soil.snowpack[11]+soil.snowpack[0])/2.0/1000.0;
            break;
     default: sthermal.hsnow19=(soil.snowpack[dm-1]+soil.snowpack[dm])/2.0/1000.0; // satisfy soil thermal model
              sthermal.hsnow29= soil.snowpack[dm]/1000.0;
              sthermal.hsnow39= (soil.snowpack[dm]+soil.snowpack[dm+1])/2.0/1000.0;
            break;
         }
//   cout << airt19 << " " << airt29 << " " << airt39 << " " << hsnow19 << " " << hsnow29 << " " << hsnow39 <<endl;

	sthermal.is9 = 10;
    sthermal.ism19 = 9;
    sthermal.smass9 = 0.f;

   int i;

	for (i=1;i <=210; i++) {
    sthermal.t9[i] = 0;
  	sthermal.xfa9[i] = -1e10f;
	sthermal.xfb9[i] = 0.f;
	sthermal.water9[i] = 1.f;
	sthermal.dx9[i] = 0.f;
	sthermal.weight9[i] = 0.f; }


   sthermal.calcu_cdsnow = 0.2; // constant for equilibrium run --- an assumpition
   sthermal.water2 = daypctp[dm][m]/100.0; // assume moss and organic have the same water content, 29/march/2001
    //printf("ttem 1864: dm = %i, m = %i, tmax = %.3f, tmin = %.3f, tave = %.3f\n", dm, m, tmax, tmin, tave);
  if (sthermal.soiltemp_(&sthermal.water2, &sthermal.calcu_cdsnow, &sthermal.airt19, &sthermal.airt29, &sthermal.airt39, &sthermal.hsnow19, &sthermal.hsnow29,
    &sthermal.hsnow39,sthermal.weight9, &sthermal.smass9, &sthermal.is9, &sthermal.ism19, sthermal.x9, sthermal.dx9,
     sthermal.xfb9, sthermal.xfa9,	sthermal.water9, sthermal.t9, &sthermal.tsoil,&sthermal.frontd,&sthermal.thawbegin,
    &sthermal.thawend,sthermal.diffsoilt, veg.cmnt, tmax, tmin, tave, atms.longmean) !=0) //,veg.tmax, veg.tmin, veg.tave, atms.longmean, soil.vsm1[dm],soil.vsm2[dm], soil.vsm3[dm], pthick_m
    { printf("bad tem"); 
//        getch();
		};

   //      ++kswitch;
				daytsoil[dm][m]=sthermal.tsoil;
		pday++;

 /*** end of calling soil thermal model ***/
//*********************************************************************************************/
}
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM::getsitecd(std::ofstream& rflog1)
{

  char ecd[80];

  std::cout << "Enter name of the site (.ECD) data file with the parameter values:";
  std::cout << std::endl;
  //std::cin >> ecd;
  fpara >> ecd;

  rflog1 << "Enter name of the site (.ECD) data file with the parameter values:";
  rflog1 << ecd << std::endl << std::endl;

  getsitecd(ecd);

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM::getsitecd(const int& numcmnt, std::ofstream& rflog1)
{

  int dv;
  char ecd[80];

  std::cout << "Enter name of the site (.ECD) data file with the parameter values  for cmnt" << std::endl;
  rflog1 << "Enter name of the site (.ECD) data file with the parameter values cmnt" << std::endl;
  for (dv = 0; dv < numcmnt; dv++)
  {
    std::cout << (dv+1) << ": ";
    //std::cin >> ecd;
    fpara >> ecd;
    
    rflog1 << (dv+1) << ": " << ecd << std::endl;

    getsitecd(dv, ecd);
  }
  rflog1 << std::endl;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM::getsitecd(char ecd[80])
{
  const int NUMVAR = 52;
  char dummy[NUMVAR][10];
  ifstream infile;
  int i;
  int dcmnt;

  int sitetveg[MAXCMNT];
  int sitewsoil[MAXCMNT];
  long update[MAXCMNT];
  char sitevegtype[MAXCMNT][31];
  char sitename[MAXCMNT][17];
  char sitetext[MAXCMNT][14];
  float sitecol[MAXCMNT];
  float siterow[MAXCMNT];

  infile.open(ecd, ios::in);

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (dcmnt = 0; dcmnt < MAXCMNT; dcmnt++)
  {
    infile >> sitetveg[dcmnt] >> sitevegtype[dcmnt];
    infile >> sitecol[dcmnt] >> siterow[dcmnt];
    infile >> sitename[dcmnt] >> sitetext[dcmnt] >> sitewsoil[dcmnt];
    infile >> vegca[dcmnt] >> vegcb[dcmnt];
    infile >> strna[dcmnt] >> strnb[dcmnt];
    infile >> solca[dcmnt] >> solcb[dcmnt];
    infile >> solna[dcmnt] >> solnb[dcmnt];
    infile >> avlna[dcmnt] >> avlnb[dcmnt];
    infile >> stona[dcmnt] >> stonb[dcmnt];
    infile >> veg.unleaf12[dcmnt];
    infile >> veg.initleafmx[dcmnt];      // Changed by DWK on 19991028
    infile >> veg.cmaxcut[dcmnt];
    infile >> veg.cmax1a[dcmnt] >> veg.cmax1b[dcmnt];
    infile >> veg.cmax2a[dcmnt] >> veg.cmax2b[dcmnt];
    infile >> veg.cfall[dcmnt];
    infile >> veg.kra[dcmnt] >> veg.krb[dcmnt];
    infile >> microbe.kda[dcmnt] >> microbe.kdb[dcmnt];
    infile >> microbe.lcclnc[dcmnt] >> microbe.propftos[dcmnt];
    infile >> veg.nmaxcut[dcmnt];
    infile >> veg.nmax1a[dcmnt] >> veg.nmax1b[dcmnt];
    infile >> veg.nmax2a[dcmnt] >> veg.nmax2b[dcmnt];
    infile >> veg.nfall[dcmnt];
    infile >> microbe.nupa[dcmnt] >> microbe.nupb[dcmnt];
    infile >> soil.nloss[dcmnt];
    infile >> microbe.nfixpar[dcmnt];
    infile >> veg.initcneven[dcmnt] >> veg.cnmin[dcmnt];
    infile >> veg.c2na[dcmnt] >> veg.c2nb[dcmnt] >> veg.c2nmin[dcmnt];
    infile >> microbe.cnsoil[dcmnt];
    infile >> update[dcmnt];

//    veg.initcneven[i] = veg.cneven[i];
//    veg.adjc2n = 1.0 + (veg.dc2n * (atms.co2[11] - atms.initco2));
//    veg.cneven[i] = veg.initcneven[i] * veg.adjc2n;
  }

  infile.close();

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM::getsitecd(const int& dv, char ecd[80])
{
  // Function added by DWK on 20000102

  char dummy[12];
  // string changed to char[80] by DWK on 20000210
  char sitename[80];
  float sitecol;
  float siterow;
  long updated;

  fecd[dv].open(ecd,ios::in);

  if (!fecd[dv])
  {
    std::cerr <<  std::endl << "Cannot open " << ecd << " for data input" <<  std::endl;
    exit(-1);
  }

  fecd[dv] >> dummy >> veg.cmnt;
  fecd[dv] >> dummy >> veg.cmnt_name;
  fecd[dv] >> dummy >> sitename;
  fecd[dv] >> dummy >> sitecol;
  fecd[dv] >> dummy >> siterow;
  fecd[dv] >> dummy >> updated;

  fecd[dv] >> dummy >> vegca[veg.cmnt];
  fecd[dv] >> dummy >> vegcb[veg.cmnt];
  fecd[dv] >> dummy >> strna[veg.cmnt];
  fecd[dv] >> dummy >> strnb[veg.cmnt];
  fecd[dv] >> dummy >> solca[veg.cmnt];
  fecd[dv] >> dummy >> solcb[veg.cmnt];
  fecd[dv] >> dummy >> solna[veg.cmnt];
  fecd[dv] >> dummy >> solnb[veg.cmnt];
  fecd[dv] >> dummy >> avlna[veg.cmnt];
  fecd[dv] >> dummy >> avlnb[veg.cmnt];
  fecd[dv] >> dummy >> stona[veg.cmnt];
  fecd[dv] >> dummy >> stonb[veg.cmnt];

  fecd[dv] >> dummy >> veg.unleaf12[veg.cmnt];
  // veg.prevleafmx changed to veg.initleafmx by DWK on 20000130
  fecd[dv] >> dummy >> veg.initleafmx[veg.cmnt];
  fecd[dv] >> dummy >> veg.cmaxcut[veg.cmnt];
  fecd[dv] >> dummy >> veg.cmax1a[veg.cmnt];
  fecd[dv] >> dummy >> veg.cmax1b[veg.cmnt];
  fecd[dv] >> dummy >> veg.cmax2a[veg.cmnt];
  fecd[dv] >> dummy >> veg.cmax2b[veg.cmnt];
  fecd[dv] >> dummy >> veg.cfall[veg.cmnt];
  fecd[dv] >> dummy >> veg.kra[veg.cmnt];
  fecd[dv] >> dummy >> veg.krb[veg.cmnt];
  fecd[dv] >> dummy >> microbe.kda[veg.cmnt];
  fecd[dv] >> dummy >> microbe.kdb[veg.cmnt];
  fecd[dv] >> dummy >> microbe.lcclnc[veg.cmnt];
  fecd[dv] >> dummy >> microbe.propftos[veg.cmnt];
  fecd[dv] >> dummy >> veg.nmaxcut[veg.cmnt];
  fecd[dv] >> dummy >> veg.nmax1a[veg.cmnt];
  fecd[dv] >> dummy >> veg.nmax1b[veg.cmnt];  
  fecd[dv] >> dummy >> veg.nmax2a[veg.cmnt];
  fecd[dv] >> dummy >> veg.nmax2b[veg.cmnt];
  fecd[dv] >> dummy >> veg.nfall[veg.cmnt];
  fecd[dv] >> dummy >> microbe.nupa[veg.cmnt];
  fecd[dv] >> dummy >> microbe.nupb[veg.cmnt];
  fecd[dv] >> dummy >> soil.nloss[veg.cmnt];
  fecd[dv] >> dummy >> microbe.nfixpar[veg.cmnt];
  fecd[dv] >> dummy >> veg.initcneven[veg.cmnt];
  fecd[dv] >> dummy >> veg.cnmin[veg.cmnt];
  fecd[dv] >> dummy >> veg.c2na[veg.cmnt];
  fecd[dv] >> dummy >> veg.c2nb[veg.cmnt];
  fecd[dv] >> dummy >> veg.c2nmin[veg.cmnt];
  fecd[dv] >> dummy >> microbe.cnsoil[veg.cmnt];

  fecd[dv].close();

};

/* *************************************************************
************************************************************** */
// YY added for MPI
void TTEM::initrun(Temstac& temstac)
{
	avlnflag=temstac.avlnflag;
	nfeed=temstac.nfeed;
	initbase=temstac.initbase;
	baseline=temstac.baseline;
	
	moistlim=temstac.moistlim;
	strteq=temstac.strteq;
	maxyears=temstac.maxyears;
	runsize=temstac.runsize;
	maxnrun=temstac.maxnrun;
	rheqflag=temstac.rheqflag;
		
	
	wtol=temstac.wtol;
	ctol=temstac.ctol;
 	ntol=temstac.ntol;
 	startyr=temstac.startyr;
 	endyr=temstac.endyr;
 	diffyr=temstac.diffyr;
    	
};


/* **************************************************************
************************************************************** */

void TTEM::initrun( std::ofstream& rflog1, const int& equil, Temstac& temstac)
{

  avlnflag = nfeed = rheqflag = 0;

/* **************************************************************
		  Run Model with Nitrogen Limitation?
************************************************************** */

  std::cout <<  std::endl << "Do you want to allow available N to fluctuate?" << std::endl;
  std::cout << "  Enter 0 for No" << std::endl;
  std::cout << "  Enter 1 for Yes: ";
  //std::cin >> avlnflag;
  fpara >> avlnflag;
  temstac.avlnflag=avlnflag; //YY added for MPI

  rflog1 << std::endl << "Do you want to allow available N to fluctuate?" << std::endl;
  rflog1 << "  Enter 0 for No" << std::endl;
  rflog1 << "  Enter 1 for Yes: " << std::endl;
  rflog1 << "avlnflag = " << avlnflag << std::endl << std::endl;

  std::cout << std::endl << "Do you want nitrogen feedback on GPP?" << std::endl;
  std::cout << "  Enter 0 for No" << std::endl;
  std::cout << "  Enter 1 for Yes: ";
  //std::cin >> nfeed;
  fpara >> nfeed;
  temstac.nfeed=nfeed; //YY added for MPI

  rflog1 << std::endl << "Do you want nitrogen feedback on GPP?" << std::endl;
  rflog1 << "  Enter 0 for No" << std::endl;
  rflog1 << "  Enter 1 for Yes: " << std::endl;
  rflog1 << "nfeed = " << nfeed << std::endl << std::endl;

  baseline = initbase = 0;
  if (nfeed == 1)
  {
    std::cout << std::endl << "Do you want to solve for baseline soil nitrogen?" << std::endl;
    std::cout << "  Enter 0 for No" << std::endl;
    std::cout << "  Enter 1 for Yes: ";
    fpara >> initbase; //std::cin >> initbase;
    baseline = initbase;

    rflog1 << std::endl << "Do you want to solve for baseline soil nitrogen?" << std::endl;
    rflog1 << "  Enter 0 for No" << std::endl;
    rflog1 << "  Enter 1 for Yes: " << std::endl;
    rflog1 << "baseline = " << baseline << std::endl << std::endl;
  }
// YY added for mpi
  temstac.initbase=initbase;
  temstac.baseline=baseline;
/* **************************************************************
			 Run Model with Moisture Limitation?
************************************************************** */

  moistlim = 0;
  std::cout << std::endl << "Do you want to run the model with moisture limitation?" << std::endl;
  std::cout << "  Enter 0 for No" << std::endl;
  std::cout << "  Enter 1 for Yes: ";
  //std::cin >> moistlim;
  fpara >> moistlim;
  temstac.moistlim=moistlim;//YY added for MPI

  rflog1 << std::endl << "Do you want to run the model with moisture limitation?" << std::endl;
  rflog1 << "  Enter 0 for No" << std::endl;
  rflog1 << "  Enter 1 for Yes: " << std::endl;
  rflog1 << "moistlim = " << moistlim << std::endl << std::endl;


/* ***************************************************************
	       Details for Steady State Conditions
************************************************************** */


  maxyears = 0;
  maxnrun = 0;

  std::cout << std::endl << "How many years do you want to wait before checking equilibrium conditions? ";
  //std::cin >> strteq;
  fpara >> strteq;
  temstac.strteq=strteq;

  rflog1 << std::endl;
  rflog1 << "How many years do you want to wait before checking equilibrium conditions? ";
  rflog1 << std::endl;
  rflog1 << "strteq = " << strteq << std::endl << std::endl;

  std::cout << std::endl << "Enter the maximum number of years for the model to run: ";
  //std::cin >> maxyears;
  fpara >> maxyears;
  temstac.maxyears=maxyears;

  rflog1 << std::endl << "Enter the maximum number of years for the model to run: ";
  rflog1 << std::endl;
  rflog1 << "maxyears = " << maxyears << std::endl << std::endl;

  runsize = maxyears;
  temstac.runsize=runsize;
  
  std::cout << std::endl << "Enter the maximum number of attempts to reach a solution: ";
  //std::cin >> maxnrun;
  fpara >> maxnrun;
  temstac.maxnrun=maxnrun;

  rflog1 << std::endl;
  rflog1 << "Enter the maximum number of attempts to reach a solution: ";
  rflog1 << std::endl;
  rflog1 << "maxnrun = " << maxnrun << std::endl << std::endl;

  if (nfeed == 0)
  {
    std::cout << std::endl << "Do you want decomposition to come into equilibrium? ";
    std::cout << "  Enter 0 for No" << std::endl;
    std::cout << "  Enter 1 for Yes: ";
    //std::cin >> rheqflag;
    fpara >> rheqflag;
    
    rflog1 << std::endl;
    rflog1 << "Do you want decomposition to come into equilibrium? " << std::endl;
    rflog1 << "  Enter 0 for No" << std::endl;
    rflog1 << "  Enter 1 for Yes: " << std::endl;
    rflog1 << "rheqflag = " << rheqflag << std::endl << std::endl;
  }

  wtol = 1000.0;
  std::cout << std::endl;
  std::cout << "What absolute tolerance do you want to use for checking equilibrium";
  std::cout << std::endl;
  std::cout << "of the water cycle? ";
  //std::cin >> wtol;
  fpara >> wtol;

  rflog1 << std::endl;
  rflog1 << "What absolute tolerance do you want to use for checking equilibrium";
  rflog1 << std::endl;
  rflog1 << "of the water cycle? wtol = " << wtol << std::endl << std::endl;

  ctol = 1000.0;
  std::cout << std::endl;
  std::cout << "What absolute tolerance do you want to use for checking equilibrium";
  std::cout << std::endl;
  std::cout << "of the carbon cycle? ";
  //std::cin >> ctol;
  fpara >> ctol;

  rflog1 << std::endl;
  rflog1 << "What absolute tolerance do you want to use for checking equilibrium";
  rflog1 << std::endl;
  rflog1 << "of the carbon cycle?" << std::endl;
  rflog1 << "ctol = " << ctol << std::endl << std::endl;

  ntol = 1000.0;
  if (nfeed == 1)
  {
    rheqflag = 1;
    std::cout << std::endl;
    std::cout << "What absolute tolerance do you want to use for checking equilibrium";
    std::cout << std::endl;
    std::cout << "of the nitrogen cycle? ";
    //std::cin >> ntol;
    fpara >> ntol;
    
    rflog1 << std::endl;
    rflog1 << "What absolute tolerance do you want to use for checking equilibrium";
    rflog1 << std::endl;
    rflog1 << "of the nitrogen cycle?" << std::endl;
    rflog1 << "ntol = " << ntol << std::endl << std::endl;
  }
    temstac.ntol=ntol;
    
  if (equil == 0)
  {

    std::cout << std::endl << std::endl;
    std::cout << "What year do you want to start collecting output data? ";
    //std::cin >> startyr;
    fpara >> startyr;

    rflog1 << std::endl << std::endl;
    rflog1 << "What year do you want to start collecting output data? ";
    rflog1 << "startyr = " << startyr << std::endl;

    std::cout << std::endl << std::endl;
    std::cout << "What year do you want to stop collecting output data? ";
    //std::cin >> endyr;
    fpara >> endyr;

    rflog1 << std::endl << std::endl;
    rflog1 << "What year do you want to stop collecting output data? ";
    rflog1 << "endyr = " << endyr << std::endl;

    std::cout << "How often (x years) should data be collected after the initial year? ";
    //std::cin >> diffyr;
    fpara >> diffyr;

    rflog1 << "How often (x years) should data be collected after the initial year? ";
    rflog1 << "diffyr = " << diffyr << std::endl;

  }
//YY added for mpi
	temstac.startyr=startyr;
	temstac.endyr=endyr;
	temstac.diffyr=diffyr;	
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TTEM::massbal(double y[NUMEQ], double prevy[NUMEQ])
{


  if (y[I_SM] != y[I_AVLW] + soil.wiltpt)
  {
    y[I_SM] = y[I_AVLW] + soil.wiltpt;
  }

  if (y[I_PCTP] != 100.0 * y[I_SM] / soil.totpor)
  {
    y[I_PCTP] = 100.0 * y[I_SM] / soil.totpor;
  }

  if (y[I_VSM] != y[I_SM] / (soil.rootz * 1000.0))
  {
    y[I_VSM] = y[I_SM] / (soil.rootz * 1000.0);
    if (y[I_VSM] <= 0.0) { y[I_VSM] = 0.001; }
  }

  if ((y[I_SNWPCK] - prevy[I_SNWPCK]) != (y[I_SNWFAL] - y[I_SNWINF]))
  {
    y[I_SNWINF] = y[I_SNWFAL] - y[I_SNWPCK] + prevy[I_SNWPCK];
  }

  if ((y[I_AVLW] - prevy[I_AVLW]) != (y[I_SNWINF] + y[I_RAIN]
       - y[I_RPERC] - y[I_EET] - y[I_SPERC]))
  {
    y[I_SPERC] = y[I_SNWINF] + y[I_RAIN] - y[I_RPERC] - y[I_EET]
                 - y[I_AVLW] + prevy[I_AVLW];
  }

  if ((y[I_RGRW] - prevy[I_RGRW]) != (y[I_RPERC] - y[I_RRUN]))
  {
    y[I_RRUN] = y[I_RPERC] - y[I_RGRW] + prevy[I_RGRW];
  }

  if ((y[I_SGRW] - prevy[I_SGRW]) != (y[I_SPERC] - y[I_SRUN]))
  {
    y[I_SRUN] = y[I_SPERC] - y[I_SGRW] + prevy[I_SGRW];
  }

  if (y[I_WYLD] != y[I_RRUN] + y[I_SRUN])
  {
    y[I_WYLD] = y[I_RRUN] + y[I_SRUN];
  }
/************************* Carbon Cycle Balances **************************/
  //printf("ttem2399: yinnpp=%.3f,ynpp=%.3f\n",y[I_INNPP],y[I_NPP]);
  if (y[I_INNPP] < y[I_NPP]) { y[I_INNPP] = y[I_NPP]; }
  //printf("ttem2401: yingpp=%.3f,ygpp=%.3f\n",y[I_INGPP],y[I_GPP]);
  if (y[I_INGPP] < y[I_GPP]) { y[I_INGPP] = y[I_GPP]; }
  //printf("ttem2403: ygrp=%.3f,ygpp=%.3f,ynpp=%.3f\n",y[I_GPR],y[I_GPP], y[I_NPP]);
  if (y[I_GPR] != y[I_GPP] - y[I_NPP])
  {
    y[I_GPR] = y[I_GPP] - y[I_NPP];
  }

  if (y[I_GPR] != y[I_RVMNT] + y[I_RVGRW])
  {
    y[I_RVGRW] = y[I_GPR] - y[I_RVMNT];
  }

  if (ag.state == 0)
  {
    // Minimum VEGC added by DWK on 20000207
    //  if (y[I_VEGC] < 0.00001) { y[I_VEGC] = 0.00001; }
    if (y[I_VEGC] - prevy[I_VEGC] != y[I_NPP] - y[I_LTRC]) 
    {
      //printf("ttem 2420: vegc=%.3f,prevvegc=%.3f,npp=%.3f,ltrc=%.3f\n",y[I_VEGC],prevy[I_VEGC],y[I_NPP],y[I_LTRC]);
      y[I_LTRC] = y[I_NPP] - y[I_VEGC] + prevy[I_VEGC];  //printf("ttem 2421: vegc=%.3f,prevvegc=%.3f,npp=%.3f,ltrc=%.3f\n",y[I_VEGC],prevy[I_VEGC],y[I_NPP],y[I_LTRC]);
    }
  }
  else
  {
    if (y[I_AGNPPC] != y[I_AGLTRC] + y[I_AGFPRDC])
    {
      y[I_AGLTRC] = y[I_AGNPPC] - y[I_AGFPRDC];
    }
  }

  // Minimum SOLC added by DWK on 20000207
  //  if (y[I_SOLC] < 0.00001) { y[I_SOLC] = 0.00001; }
  //printf("ttem 2434:rh=%.3f,ltrc=%.3f,agltrc=%.3f,slashc=%.3f,scnvrtc=%.3f,ysolc=%.3f,pysolc=%.3f\n", y[I_RH],y[I_LTRC],y[I_AGLTRC],y[I_SLASHC],y[I_SCNVRTC],y[I_SOLC], prevy[I_SOLC]);
  if (y[I_SOLC] - prevy[I_SOLC] != y[I_LTRC] + y[I_AGLTRC]
      + y[I_SLASHC] - y[I_SCNVRTC] - y[I_RH])
  {
    y[I_RH] = y[I_LTRC] + y[I_AGLTRC] + y[I_SLASHC] - y[I_SCNVRTC]
    - y[I_SOLC] + prevy[I_SOLC]  ; // + atms.frontddiff*14500   
  }
  //printf("ttem 2441:nep=%.3f,npp=%.3f,agnppc=%.3f,rh=%.3f\n",y[I_NEP],y[I_NPP],y[I_AGNPPC],y[I_RH]);
  if (y[I_NEP] != y[I_NPP] + y[I_AGNPPC] - y[I_RH])
  {
    y[I_NEP] = y[I_NPP] + y[I_AGNPPC] - y[I_RH];
  }
  //printf("ttem2446: cflx=%.3f,npp=%.3f,cnvertc=%.3f\n",y[I_CFLX],y[I_NEP],y[I_CNVRTC]);
  if (y[I_CFLX] != y[I_NEP] - y[I_CNVRTC])
  {
    y[I_CFLX] = y[I_NEP] - y[I_CNVRTC];
  }

  /*********************Nitrogen Cycle Balances**********************/

  if (y[I_VNUP] < 0.0 ) { y[I_VNUP] = 0.0; }

  //if (y[I_VNUP] > y[I_INNUP]) { y[I_VNUP] = y[I_INNUP]; }
  if (y[I_INNUP] < y[I_VNUP]) { y[I_INNUP] = y[I_VNUP]; }

  if (y[I_VSUP] < 0.0) { y[I_VSUP] = 0.0; }

  if (y[I_VSUP] > y[I_VNUP]) { y[I_VSUP] = y[I_VNUP]; }

  if (y[I_VLUP] != y[I_VNUP] - y[I_VSUP])
  {
    y[I_VLUP] = y[I_VNUP] - y[I_VSUP];
  }

  if (ag.state == 0)
  {
  if ((y[I_STON] + y[I_STRN]) * veg.cneven[veg.cmnt] < y[I_VEGC]) 
    {
    double scaler = y[I_VEGC] / (y[I_STON] + y[I_STRN]) / veg.cneven[veg.cmnt]; 
    y[I_STON] = y[I_STON] * scaler;
    y[I_STRN] = y[I_STRN] * scaler;   
    }
   if ((y[I_STON] + y[I_STRN]) * veg.cnmin[veg.cmnt] > y[I_VEGC]) 
    {
    double scaler = y[I_VEGC] / (y[I_STON] + y[I_STRN]) / veg.cnmin[veg.cmnt];  
    y[I_STON] = y[I_STON] * scaler;
    y[I_STRN] = y[I_STRN] * scaler;   
    //printf("ttem2476: scaler=%.3f,ston=%.3f,strn=%.3f\n", scaler, y[I_STON],y[I_STRN]);
    }   
    if (y[I_STON] - prevy[I_STON] != y[I_VLUP] + y[I_VNRSRB] - y[I_VNMBL])
    {
      y[I_VNRSRB] = y[I_STON] - prevy[I_STON] + y[I_VNMBL] - y[I_VLUP];
    }

    if (y[I_STRN] - prevy[I_STRN] != y[I_VSUP] - y[I_LTRN]
        - y[I_VNRSRB] + y[I_VNMBL])
    {
      y[I_LTRN] = y[I_VSUP] - y[I_STRN] + prevy[I_STRN]
                  - y[I_VNRSRB] + y[I_VNMBL];
	  //printf("ttem 2465: yLTRN = %.3f, yVSUP = %.3f, ySTRN = %.3f, pySTRN = %.3f, yVNRSRB = %.3f, yVNMBL = %.3f\n", y[I_LTRN], y[I_VSUP], y[I_STRN], prevy[I_STRN], y[I_VNRSRB], y[I_VNMBL]);
    }
  }
  else
  {
    if (y[I_AGNPPN] != y[I_AGLTRN] + y[I_AGFPRDN])
    {
      y[I_AGLTRN] = y[I_AGNPPN] - y[I_AGFPRDN];
    }
  }

  if (y[I_SOLN] - prevy[I_SOLN] != y[I_LTRN] + y[I_AGLTRN] + y[I_SLASHN] 
      - y[I_NMIN] - y[I_SCNVRTN] - y[I_NSRTNT])
  {
    y[I_NMIN] = y[I_LTRN] + y[I_AGLTRN] + y[I_SLASHN] - y[I_SCNVRTN] 
                - y[I_NSRTNT] - y[I_SOLN] + prevy[I_SOLN] ;  //  + (atms.frontddiff)*1000
  }
  //printf("ttem 2482: ySOLN = %.3f, prevySOLN = %.3f, yLTRN = %.3f, yAGLTRN = %.3f, ySLASHN = %.3f, yNMIN = %.3f, ySCNVRTN = %.3f, yNSRTNT = %.3f,\n", y[I_SOLN], prevy[I_SOLN], y[I_LTRN], y[I_AGLTRN], y[I_SLASHN], y[I_NMIN], y[I_SCNVRTN],y[I_NSRTNT]);

  if (y[I_NLST] < 0.0) { y[I_NLST] = 0.0; }
  	
 
  if (ag.state == 0)
  {
    if (y[I_NINP] < 0.0) { y[I_NINP] = 0.0; }
    //printf("ttem2499: yavln=%.3f,preavln=%.3f,ninp=%.3f,nlost=%.3f,nmin=%.3f,n2fix=%.3f,ngas=%.3f,vnup=%.3f\n", y[I_AVLN],prevy[I_AVLN],y[I_NINP],y[I_NLST],y[I_NMIN],y[I_N2FIX],y[I_NGAS],y[I_VNUP]);
    if (y[I_AVLN] - prevy[I_AVLN] != y[I_NINP] - y[I_NLST] + y[I_NMIN] + y[I_N2FIX] - y[I_NGAS] - y[I_VNUP])
    {
      if (y[I_NINP] + y[I_NMIN] + y[I_N2FIX] - y[I_VNUP] - y[I_AVLN] + prevy[I_AVLN] - y[I_NGAS] > 0.0)
      {
        y[I_NLST] =  y[I_NINP] + y[I_NMIN] + y[I_N2FIX] - y[I_VNUP] - y[I_AVLN] - y[I_NGAS] + prevy[I_AVLN];
      }
      else
      {
        y[I_NINP] = y[I_NLST] - y[I_NMIN] - y[I_N2FIX]+ y[I_VNUP] + y[I_AVLN] + y[I_NGAS]- prevy[I_AVLN];
        //printf("ttem2509: yavln=%.3f,preavln=%.3f,ninp=%.3f,nlost=%.3f,nmin=%.3f,n2fix=%.3f,ngas=%.3f,vnup=%.3f\n", y[I_AVLN],prevy[I_AVLN],y[I_NINP],y[I_NLST],y[I_NMIN],y[I_N2FIX],y[I_NGAS],y[I_VNUP]);
      }
    }
  }

  if(ag.state == 1)
  {

// Following condition added to correct for negative AGFERTN
// D. Kicklighter 19990721

    if (y[I_NRETNT] != y[I_NVRTNT] + y[I_NSRTNT])
    {
      y[I_NRETNT] = y[I_NVRTNT] + y[I_NSRTNT];
    }

    if (y[I_AGFRTN] < 0.0) { y[I_AGFRTN] = 0.0; }

    if (y[I_NINP] != y[I_AGFRTN] + y[I_NRETNT])
    {
      y[I_NINP] = y[I_AGFRTN] + y[I_NRETNT];
    }

    if (y[I_NINP] < 0.0) { y[I_NINP] = 0.0; }

    if (y[I_AVLN] - prevy[I_AVLN] != y[I_NINP] - y[I_NLST] + y[I_NMIN] + y[I_N2FIX] - y[I_NGAS]- y[I_AGNPPN])
    {
      y[I_NLST] =  prevy[I_AVLN] - y[I_AVLN] + y[I_NINP] + y[I_NMIN] + y[I_N2FIX] - y[I_NGAS]- y[I_AGNPPN];
    }
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM::monthxclm(const int& dcmnt, const double& tgppopt, const int& dm)
{

double raq10;


/* temp: effect of temperature on primary productivity */

  if (atms.tair[dm] <= veg.tmin[dcmnt] || atms.tair[dm] >= veg.tmax[dcmnt])
  {
    temp = 0.0;
  }
  else
  {
    if (atms.tair[dm] >= tgppopt && atms.tair[dm] <= veg.toptmax[dcmnt])
    {
      temp = 1.0;
    }
    else
    {
      if (atms.tair[dm] > veg.tmin[dcmnt] && atms.tair[dm] < tgppopt)
      {
	temp = (atms.tair[dm]-veg.tmin[dcmnt])*(atms.tair[dm]-veg.tmax[dcmnt])
               /((atms.tair[dm]-veg.tmin[dcmnt])*(atms.tair[dm]-veg.tmax[dcmnt])
               - pow((atms.tair[dm]-tgppopt),2.0));
      }
      else
      {
	temp = (atms.tair[dm]-veg.tmin[dcmnt])*(atms.tair[dm]-veg.tmax[dcmnt])
               /((atms.tair[dm]-veg.tmin[dcmnt])*(atms.tair[dm]-veg.tmax[dcmnt])
               - pow((atms.tair[dm]-veg.toptmax[dcmnt]),2.0));
      }
    }
  }


/* respq10: effect of temperature on plant respiration  */

  raq10 = veg.raq10a0[dcmnt] + (veg.raq10a1[dcmnt]*atms.tair[dm])
          + (veg.raq10a2[dcmnt]*pow(atms.tair[dm],2.0))
          + (veg.raq10a3[dcmnt]*pow(atms.tair[dm],3.0));
  respq10 = pow(raq10,atms.tair[dm]/10.0);
  //printf("ttem2581: respq10=%.3f, tair=%.3f, raq10=%.3f,raq10a0=%.3f,raq10a1=%.3f,raq10a2=%.3f,raq10a3=%.3f\n",respq10,atms.tair[dm],raq10,veg.raq10a0[dcmnt],veg.raq10a1[dcmnt],veg.raq10a2[dcmnt],veg.raq10a3[dcmnt]);


/* dq10: effect of temperature on decomposition */

// 19990821-previous year snowpack (soil.snowpack[dm]) changed
//          to previous month snowpack (y[NUMEEQ+2]) by Jim Long



 if (tmflgs.stmflg==1) {
    dq10 = pow(microbe.rhq10[dcmnt],atms.tsoil[dm]/10.0);
    //if (0 < atms.tsoil[dm] < 10) dq10 = pow(microbe.rhq10[dcmnt],(atms.tsoil[dm]/10.0)+1);
    //else dq10 = pow(microbe.rhq10[dcmnt],atms.tsoil[dm]/10.0);  // Replace air Temperature with Soil Temperature
    }
  else
   {
//  if (y[I_SNWPCK] > 0.0) { dq10 = 1.0; }
//  else
 // {
     dq10 = pow(microbe.rhq10[dcmnt],atms.tsoil[dm]/10.0);
     //if (0 < atms.tsoil[dm] < 10) dq10 = pow(microbe.rhq10[dcmnt],(atms.tsoil[dm]/10.0)+1);
     //else dq10 = pow(microbe.rhq10[dcmnt],atms.tair[dm]/10.0);
   }
   //printf("ttem2515: dm=%i, dq10=%.3f, rhq10=%.3f, tsoil=%.3f\n", dm, dq10, microbe.rhq10[dcmnt],atms.tsoil[dm]);
//  if (y[I_SNWPCK] > 0.0) { dq10 = 1.0; }
//  else
//  {
//    dq10 = pow(microbe.rhq10[dcmnt],atms.tair[dm]/10.0);
//  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void TTEM::resetODEflux(double y[])
{
  int i;

  for (i = MAXSTATE; i < NUMEQ; i++) { y[i] = 0.0; }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TTEM::resetYrFlux(void)
{

  int dm;

  // Annual carbon storage
  veg.yrcarbon = 0.0;
  soil.yrorgc = 0.0;
  veg.yrcarbonb = 0.0;
  
  // Annual nitrogen storage
  veg.yrnitrogen = 0.0;
  veg.yrstructn = 0.0;
  veg.yrc2n = 0.0;
  veg.yrstoren = 0.0;
  soil.yrorgn = 0.0;
  soil.yrc2n = 0.0;
  soil.yravln = 0.0;
  microbe.yrnh4 = 0.0;
  microbe.yrno3 = 0.0;

  // Annual carbon & nitrogen storage in agricultural ecosystems
  ag.PROD1.carbon = 0.0;
  ag.PROD1.nitrogen = 0.0;

  // Annual water storage
  soil.yravlh2o = 0.0;
  soil.yrrgrndh2o = 0.0;
  soil.yrsnowpack = 0.0;
  soil.yrsgrndh2o = 0.0;
  soil.yrsmoist = 0.0;
  soil.yrpctp = 0.0;
  soil.yrvsm = 0.0;

  // Annual carbon fluxes
  veg.yringpp = 0.0;
  veg.yrgpp = 0.0;
  veg.yrinnpp = 0.0;
  veg.yrnpp = 0.0;
  veg.yrltrc = 0.0;
  microbe.yrrh = 0.0;
  yrnep = 0.0;

  // Annual nitrogen fluxes
  soil.yrnin = 0.0;
  veg.yrinnup = 0.0;
  veg.yrnup = 0.0;
  veg.yrsup = 0.0;
  veg.yrlup = 0.0;
  veg.yrnmobil = 0.0;
  veg.yrnrsorb = 0.0;
  veg.yrltrn = 0.0;
  microbe.yrnmin = 0.0;
  soil.yrnlost = 0.0;
  soil.yrngas = 0.0;

  // Annual water fluxes
  atms.yrrain = 0.0;
  soil.yrrperc = 0.0;
  soil.yrrrun = 0.0;
  atms.yrsnowfall = 0.0;
  soil.yrsnowinf = 0.0;
  soil.yrsperc = 0.0;
  soil.yrsrun = 0.0;
  atms.yrpet = 0.0;
  atms.yreet = 0.0;
  soil.yrh2oyld = 0.0;

     // for soil temperature
  atms.yrfrontd =0.0;
  atms.yrthawbegin =0.0;
  atms.yrthawend =0.0;
  atms.yrtsoil = 0.0;
  atms.yrdst5 = 0.0;
  atms.yrdst10 = 0.0;
  atms.yrdst20 = 0.0;
  atms.yrdst50 = 0.0;
  atms.yrdst100 = 0.0;
  atms.yrdst200 = 0.0;


  // Phenology
  veg.yrunleaf = 0.0;
  veg.yrleaf = 0.0;
  veg.yrfpc = 0.0;

  // Annual carbon and nitrogen fluxes from agricultural
  // conversion
  ag.yrconvrtC = 0.0;
  ag.yrvconvrtC = 0.0;
  ag.yrsconvrtC = 0.0;
  ag.yrconvrtN = 0.0;
  ag.yrvconvrtN = 0.0;
  ag.yrsconvrtN = 0.0;
  ag.yrnrent = 0.0;
  ag.yrnvrent = 0.0;
  ag.yrnsrent = 0.0;
  ag.yrslashC = 0.0;
  ag.yrslashN = 0.0;
  ag.formPROD10.carbon  = 0.0;
  ag.formPROD10.nitrogen  = 0.0;
  ag.formPROD100.carbon = 0.0;
  ag.formPROD100.nitrogen = 0.0;

  // Annual carbon and nitrogen fluxes from agriculture
  ag.yrnppC = 0.0;
  ag.yrnppN = 0.0;
  ag.formPROD1.carbon = 0.0;
  ag.formPROD1.nitrogen = 0.0;
  ag.yrltrc = 0.0;
  ag.yrltrn = 0.0;
  ag.yrfertn = 0.0;

  // Annual carbon and nitrogen fluxes from human products
  ag.formTOTPROD.carbon = 0.0;
  ag.formTOTPROD.nitrogen = 0.0;

  // Monthly carbon and nitrogen fluxes
  for (dm = 0; dm < CYCLE; dm++)
  {
    soil.ninput[dm] = 0.0;
    soil.nlost[dm] = 0.0;
    soil.ngas[dm] = 0.0;
    

    ag.slash[dm].carbon = 0.0;
    ag.slash[dm].nitrogen = 0.0;
    ag.vconvrtflx[dm].carbon = 0.0;
    ag.sconvrtflx[dm].carbon = 0.0;
    ag.convrtflx[dm].carbon = 0.0;
    ag.vconvrtflx[dm].nitrogen = 0.0;
    ag.sconvrtflx[dm].nitrogen = 0.0;
    ag.convrtflx[dm].nitrogen = 0.0;
    ag.nvretent[dm] = 0.0;
    ag.nsretent[dm] = 0.0;
    ag.nretent[dm] = 0.0;
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TTEM::setELMNTecd(const int& kdinflg, const int& dcmnt,
                       const double& psiplusc)
{


// Initialize TEM parameters dependent upon a grid cell's soil texture

  if (psiplusc <= veg.cmaxcut[dcmnt])
  {
    veg.cmax = (veg.cmax1a[dcmnt] * psiplusc) + veg.cmax1b[dcmnt];
  }
  else
  {
    veg.cmax = (veg.cmax2a[dcmnt] * psiplusc) + veg.cmax2b[dcmnt];
  }

  if (kdinflg == 0)
  {
    microbe.kdc = (microbe.kda[dcmnt] / psiplusc) + microbe.kdb[dcmnt];
  }

  if (psiplusc <= veg.nmaxcut[dcmnt])
  {
    veg.nmax = (veg.nmax1a[dcmnt] * psiplusc) + veg.nmax1b[dcmnt];
  }
  else
  {
    veg.nmax = (veg.nmax2a[dcmnt] * psiplusc) + veg.nmax2b[dcmnt];
  }

  microbe.nup = (microbe.nupa[dcmnt] / psiplusc) + microbe.nupb[dcmnt];  
  //printf("ttem2814:nup=%.3f,nupa=%.3f,psiplusc=%.3f,nupb=%.3f\n", microbe.nup,microbe.nupa[dcmnt],psiplusc, microbe.nupb[dcmnt]);
// Set initial maximum relative leaf area

  veg.prvleafmx[dcmnt] = veg.initleafmx[dcmnt];  // Changed by DWK on 19991028


// Determine the "decay" parameter

  microbe.decay = 0.26299 + (1.14757*microbe.propftos[dcmnt])
                  - (0.42956*pow(microbe.propftos[dcmnt],2.0));
  //printf("ttem2732:decay=.3%f, microbe.propftos.3%f\n", microbe.decay, microbe.propftos[dcmnt]);
  
  // set inital maximum active layer thickness (m), Lei, 7/8/20
	atms.frontdmx = 0;
 	atms.frontdymx = 0;
	atms.frontddiff = 0;

};
 
/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TTEM::setELMNTevap(const int& stateflg, const int& dcmnt,
                        double pet[CYCLE], double tair[CYCLE])
{

  int i;

// Determine initial values for atms.prvpetmx, atms.prveetmx,
//   and veg.topt

  if (stateflg == 1)
  {
    for (i = 0; i < CYCLE; i++)
    {
      if (pet[i] >= atms.prvpetmx)
      {
	     veg.topt = tair[i];
      }
      if (veg.aleaf[dcmnt] == 0.0 && veg.bleaf[dcmnt] == 0.0
         && veg.cleaf[dcmnt] == 1.0)
      {
	     if (tair[i] > veg.topt) { veg.topt = tair[i]; }
      }
      else
      {
	     if (veg.unnormleaf[i] >= veg.prvleafmx[dcmnt])
        {
	       veg.topt = tair[i];
	     }
      }
    }

    if (veg.topt > veg.toptmax[dcmnt]) { veg.topt = veg.toptmax[dcmnt]; }
    if (veg.topt < veg.toptmin[dcmnt]) { veg.topt = veg.toptmin[dcmnt]; }
  }

  else
  {
    for (i = 0; i < CYCLE; i++)
    {
      atms.pet[i] = atms.petjh(atms.nirr[i], tair[i], i);
      if (i == 0)
      {
	atms.prvpetmx = atms.prveetmx = atms.pet[0];
	veg.topt  = tair[0];
      }
      else
     {
	if (atms.pet[i] > atms.prvpetmx)
        {
	  atms.prvpetmx = atms.prveetmx = atms.pet[i];
	  veg.topt = tair[i];
	}
      }
    }

    atms.yrpet = 1.0;
    atms.yreet = 1.0;
  }


// Determine vegetation C/N parameter as a function of vegetation type,
// annual PET, and annual EET (annual EET initially equal to yrpet)

  veg.updateC2N(dcmnt,atms.yreet,atms.yrpet,atms.co2[11],atms.initco2);

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TTEM::setELMNTflux(void)
{

  // Initialize carbon, nitrogen and water fluxes including
  // ODE state variables (i.e., y[])
  int dm;

  for (dm = 0; dm < CYCLE; dm++)
  {
    // Initialize carbon fluxes in natural ecosystems to zero

    y[I_INGPP] = veg.ingpp[dm] = 0.0;
    y[I_GPP] = veg.gpp[dm] = 0.0;
    y[I_INNPP] = veg.innpp[dm] = 0.0;
    y[I_NPP] = veg.npp[dm] = 0.0;
    y[I_GPR] = veg.gpr[dm] = 0.0;
    y[I_RVMNT] = veg.rm[dm] = 0.0;
    y[I_RVGRW] = veg.rg[dm] = 0.0;
    y[I_LTRC] = veg.ltrfal[dm].carbon = 0.0;
    y[I_RH] = microbe.rh[dm] = 0.0;
    y[I_NEP] = nep[dm] = 0.0;

    // Initialize nitrogen fluxes in natural ecosystems to zero

    y[I_NINP] = soil.ninput[dm] = 0.0;
    y[I_INNUP] = veg.inuptake = 0.0;
    y[I_VNUP] = veg.nuptake[dm] = 0.0;
    y[I_VSUP] = veg.suptake[dm] = 0.0;
    y[I_VLUP] = veg.luptake[dm] = 0.0;
    y[I_VNMBL] = veg.nmobil[dm] = 0.0;
    y[I_VNRSRB] = veg.nresorb[dm] = 0.0;
    y[I_LTRN] = veg.ltrfal[dm].nitrogen = 0.0;
    y[I_MNUP] = microbe.nuptake[dm] = 0.0;
    y[I_NMIN] = microbe.netnmin[dm] = 0.0;
    y[I_FNIT] = microbe.Fnit[dm] = 0.0;
    y[I_FDENIT] = microbe.Fdenit[dm] = 0.0;
    y[I_NLST] = soil.nlost[dm] = 0.0;
    y[I_N2] = microbe.N2[dm] = 0.0;
    y[I_N2O] = microbe.N2O[dm] = 0.0;
    y[I_NOX] = microbe.NOx[dm] = 0.0;
    y[I_NGAS] = soil.ngas[dm] = 0.0;
    y[I_N2F] = microbe.N2F[dm] = 0.0;
    y[I_N2OF] = microbe.N2OF[dm] = 0.0;
    y[I_N2AIR] = microbe.N2Air[dm] = 0.0;
    y[I_N2OAIR] = microbe.N2OAir[dm] = 0.0;
    y[I_N2OUPT] = microbe.N2OUpt[dm] = 0.0;
    y[I_N2ON] = microbe.N2ON[dm] = 0.0;
    y[I_N2ODN] = microbe.N2ODN[dm] = 0.0;
    y[I_N2FIX] = microbe.N2Fix[dm] = 0.0;
    y[I_N2FIXS] = microbe.N2Fixs[dm] = 0.0;  
    y[I_N2FIXFL] = microbe.N2Fixfl[dm] = 0.0;
    y[I_NH3] = microbe.NH3[dm] = 0.0;
    // Initialize water fluxes to zero

    y[I_RAIN] = 0.0;
    y[I_RPERC] = soil.rperc[dm] = 0.0;
    y[I_RRUN] = soil.rrun[dm] = 0.0;
    y[I_SNWFAL] = 0.0;
    y[I_SNWINF] = soil.snowinf[dm] = 0.0;
    y[I_SPERC] = soil.sperc[dm] = 0.0;
    y[I_SRUN] = soil.srun[dm] = 0.0;
    y[I_PET] = 0.0;
    y[I_EET] = atms.eet[dm] = 0.0;
    y[I_WYLD] = soil.h2oyld[dm] = 0.0;

// for soil temperature
  y[I_TSOIL] = atms.tsoil[dm]=0.0;
  y[I_DST5] = atms.dst5[dm]=0.0;
  y[I_DST10] = atms.dst10[dm]=0.0;
  y[I_DST20] = atms.dst20[dm]=0.0;
  y[I_DST50] = atms.dst50[dm]=0.0;
  y[I_DST100] = atms.dst100[dm]=0.0;
  y[I_DST200] = atms.dst200[dm]=0.0;
  y[I_FRONTD] = atms.frontd[dm]=0.0;
  y[I_THAWBE] = atms.thawbe[dm]=0.0;
  y[I_THAWEND] = atms.thawend[dm]=0.0;
// end of ...


    // Initialize carbon and nitrogen fluxes during conversion
    //  to zero

    y[I_CNVRTC] = ag.convrtflx[dm].carbon = 0.0;
    y[I_CNVRTN] = ag.convrtflx[dm].nitrogen = 0.0;
    y[I_SCNVRTC] = ag.sconvrtflx[dm].carbon = 0.0;
    y[I_SCNVRTN] = ag.sconvrtflx[dm].nitrogen = 0.0;
    y[I_NVRTNT] = ag.nvretent[dm] = 0.0;
    y[I_NSRTNT] = ag.nsretent[dm] = 0.0;
    y[I_NRETNT] = ag.nretent[dm] = 0.0;
    y[I_SLASHC] = ag.slash[dm].carbon = 0.0;
    y[I_SLASHN] = ag.slash[dm].nitrogen = 0.0;
    y[I_PRDF10C] = ag.formPROD10.carbon = 0.0;
    y[I_PRDF10N] = ag.formPROD10.nitrogen = 0.0;
    y[I_PRDF100C] = ag.formPROD100.carbon = 0.0;
    y[I_PRDF100N] = ag.formPROD100.nitrogen = 0.0;

    // Initialize carbon and nitrogen in agricultural ecosystems
    //   to zero

    y[I_AGNPPC] = ag.npp[dm].carbon = 0.0;
    y[I_AGNPPN] = ag.npp[dm].nitrogen = 0.0;
    y[I_AGFPRDC] = ag.formPROD1.carbon = 0.0;
    y[I_AGFPRDN] = ag.formPROD1.nitrogen = 0.0;
    y[I_AGLTRC] = ag.ltrfal[dm].carbon = 0.0;
    y[I_AGLTRN] = ag.ltrfal[dm].nitrogen = 0.0;
    y[I_AGFRTN] = ag.fertn[dm] = 0.0;

    // Initialize carbon and nitrogen from human products to zero

    y[I_TOTFPRDC] = ag.formTOTPROD.carbon = 0.0;
    y[I_TOTFPRDN] = ag.formTOTPROD.nitrogen = 0.0;
    y[I_AGPRDFC] = ag.PROD1decay.carbon = 0.0;
    y[I_AGPRDFN] = ag.PROD1decay.nitrogen = 0.0;
    y[I_PRD10FC] = ag.PROD10decay.carbon = 0.0;
    y[I_PRD10FN] = ag.PROD10decay.nitrogen = 0.0;
    y[I_PRD100FC] = ag.PROD100decay.carbon = 0.0;
    y[I_PRD100FN] = ag.PROD100decay.nitrogen = 0.0;
    y[I_TOTPRDFC] = ag.TOTPRODdecay.carbon = 0.0;
    y[I_TOTPRDFN] = ag.TOTPRODdecay.nitrogen = 0.0;

    // Initialize integrated carbon fluxes to zero

    y[I_TOTNPP] = ag.totnpp[dm] = 0.0;
    y[I_CFLX] = cflux[dm] = 0.0;

  }
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM::setMonth(int& dm, double y[])
{

  // Carbon pools
  soil.doc[dm] = y[I_DOC];
  veg.plant[dm].carbon = y[I_VEGC];  //printf("ttem 3069: dm=%i, vegc=%.3f, vegplantcarbon=%.3f\n", dm, y[I_VEGC],veg.plant[dm].carbon);
  soil.org[dm].carbon = y[I_SOLC];
  totalc[dm] = veg.plant[dm].carbon + soil.org[dm].carbon;
  veg.planta[dm].carbon = y[I_VEGCAG];  //printf("ttem 3033: dm=%i, vegca=%.3f, vegplantbcarbon=%.3f\n", dm, y[I_VEGCAG],veg.planta[dm].carbon);
  veg.plantb[dm].carbon = y[I_VEGCBG];
  
  // Nitrogen pools
  veg.strctrl[dm].nitrogen = y[I_STRN]; //printf("ttem3061: dm=%i, vegstrctrln=%.3f,ystern=%.3f\n", dm, veg.strctrl[dm].nitrogen, y[I_STRN]);
  veg.labile[dm].nitrogen = y[I_STON];
  veg.plant[dm].nitrogen = veg.strctrl[dm].nitrogen + veg.labile[dm].nitrogen;
  soil.org[dm].nitrogen = y[I_SOLN];
  soil.availn[dm] = y[I_AVLN];   //printf("ttem3066:availn=%.3f,yravln=%.3f\n",soil.availn[dm],y[I_AVLN]);
  microbe.NH4[dm] = y[I_NH4];
  microbe.NO3[dm] = y[I_NO3];

  // Water pools
  soil.avlh2o[dm] = y[I_AVLW];
  soil.rgrndh2o[dm] = y[I_RGRW];
  soil.snowpack[dm] = y[I_SNWPCK];
  soil.sgrndh2o[dm] = y[I_SGRW];
  soil.moist[dm] = y[I_SM];
  soil.pctp[dm] = y[I_PCTP];
  soil.vsm[dm] = y[I_VSM];

  // Monthly carbon fluxes in natural ecosystems
  veg.ingpp[dm] = y[I_INGPP]; //printf("ttem3090: ingpp=%.3f,INGPP=%.3f\n",veg.ingpp[dm],y[I_INGPP]);
  veg.gpp[dm] = y[I_GPP];   //printf("ttem3091: gpp=%.3f,IGPP=%.3f\n",veg.gpp[dm],y[I_GPP]);
  veg.innpp[dm] = y[I_INNPP];
  veg.npp[dm] = y[I_NPP];
  veg.gpr[dm] = y[I_GPR];
  veg.rm[dm] = y[I_RVMNT];
  veg.rg[dm] = y[I_RVGRW];
  veg.ltrfal[dm].carbon = y[I_LTRC];
  microbe.rh[dm] = y[I_RH];
  nep[dm] = y[I_NEP];

  // Monthly nitrogen fluxes in natural ecosystems
  soil.ninput[dm] = y[I_NINP];  //printf("ttem3100: ninput=%.3f,ININP=%.3f\n",soil.ninput[dm],y[I_NINP]);
  veg.inuptake = y[I_INNUP];
  veg.nuptake[dm] = y[I_VNUP];  //printf("ttem3094 vegnupt=%.3f,IVUNP=%.3f\n",veg.nuptake[dm],y[I_VNUP]);
  veg.suptake[dm] = y[I_VSUP];
  veg.luptake[dm] = y[I_VLUP];
  veg.nmobil[dm] = y[I_VNMBL];
  veg.nresorb[dm] = y[I_VNRSRB];
  veg.ltrfal[dm].nitrogen = y[I_LTRN];
  microbe.nuptake[dm] = y[I_MNUP];
  microbe.netnmin[dm] = y[I_NMIN];
  microbe.Fnit[dm] = y[I_FNIT];
  microbe.Fdenit[dm] = y[I_FDENIT];
  soil.nlost[dm] = y[I_NLST];
  microbe.N2[dm] = y[I_N2];
  microbe.N2O[dm] = y[I_N2O];
  microbe.NOx[dm] = y[I_NOX];
  soil.ngas[dm] = y[I_NGAS];
  microbe.N2OF[dm] = y[I_N2OF];  //YY
  microbe.N2F[dm] = y[I_N2F];  //YY
  microbe.N2OAir[dm] = y[I_N2OAIR];  //printf("ttem 2984: microbe.N2OAir = %.3f\n", microbe.N2OAir[dm]);//YY
  microbe.N2Air[dm] = y[I_N2AIR];  //YY
  microbe.N2OUpt[dm] = y[I_N2OUPT];  //YY
  microbe.N2ON[dm] = y[I_N2ON];  //YY
  microbe.N2ODN[dm] = y[I_N2ODN];  //YY
  microbe.N2Fix[dm] = y[I_N2FIX];  //YY
  microbe.N2Fixs[dm] = y[I_N2FIXS];  
  microbe.N2Fixfl[dm] = y[I_N2FIXFL];  //printf("ttem 3118: microbe.N2Fixfl = %.3f\n", microbe.N2Fixfl[dm]);
  microbe.NH3[dm] = y[I_NH3];  //YY
    
  // Monthly water fluxes
  atms.rain[dm] = y[I_RAIN];
  soil.rperc[dm] = y[I_RPERC];
  soil.rrun[dm] = y[I_RRUN];
  atms.snowfall[dm] = y[I_SNWFAL];
  soil.snowinf[dm] = y[I_SNWINF];
  soil.sperc[dm] = y[I_SPERC];
  soil.srun[dm] = y[I_SRUN];
  atms.pet[dm] = y[I_PET];
  atms.eet[dm] = y[I_EET];
  soil.h2oyld[dm] = y[I_WYLD];

  // added for soil thermal model

  atms.tsoil[dm]= y[I_TSOIL];
  atms.dst5[dm]= y[I_DST5];
  atms.dst10[dm]= y[I_DST10];
  atms.dst20[dm]= y[I_DST20];
  atms.dst50[dm]= y[I_DST50];
  atms.dst100[dm]= y[I_DST100];
  atms.dst200[dm]= y[I_DST200];
  atms.frontd[dm]= y[I_FRONTD];
  atms.thawbe[dm]= y[I_THAWBE];
  atms.thawend[dm]= y[I_THAWEND];

// end of adding


  // Monthly phenology
  veg.unnormleaf[dm] = y[I_UNRMLF];  //printf("ttem3155: dm=%i,unnormleaf=%.3f,yUNRMLF=%.3f\n",dm,veg.unnormleaf[dm],y[I_UNRMLF]);
  veg.leaf[dm] = y[I_LEAF];

  // Monthly carbon and nitrogen fluxes associated with
  //  agricultural conversion
  ag.convrtflx[dm].carbon = y[I_CNVRTC];
  ag.convrtflx[dm].nitrogen = y[I_CNVRTN];
  ag.sconvrtflx[dm].carbon = y[I_SCNVRTC];
  ag.sconvrtflx[dm].nitrogen = y[I_SCNVRTN];
  ag.nvretent[dm] = y[I_NVRTNT];
  ag.nsretent[dm] = y[I_NSRTNT];
  ag.nretent[dm] = y[I_NRETNT];

  // Monthly carbon and nitrogen fluxes from agricultural
  //   ecosystems
  ag.npp[dm].carbon = y[I_AGNPPC];
  ag.npp[dm].nitrogen = y[I_AGNPPN];
  ag.ltrfal[dm].carbon = y[I_AGLTRC];
  ag.ltrfal[dm].nitrogen = y[I_AGLTRN];
  ag.fertn[dm] = y[I_AGFRTN];

  // Monthly integrated carbon fluxes
  ag.totnpp[dm] = y[I_TOTNPP];
  cflux[dm] = y[I_CFLX];

  // Update sum of annual carbon storage
  veg.yrcarbon  += y[I_VEGC];  //printf("ttem 3148: vegc=%.3f\n", y[I_VEGC]);
  soil.yrorgc += y[I_SOLC];
  veg.yrcarbona += y[I_VEGCAG];
  veg.yrcarbonb += y[I_VEGCBG]; //printf("ttem 3150: vegyrcb=%.3f, vegcb=%.3f\n", veg.yrcarbonb , y[I_VEGCBG]);
  
  // Update sum of annual nitrogen storage
  veg.yrnitrogen  += y[I_STRN] + y[I_STON];
  veg.yrstructn += y[I_STRN];
  soil.yrorgn += y[I_SOLN]; //printf("ttem 3212: soln=%.3f\n",y[I_SOLN]);
  soil.yravln  += y[I_AVLN];  //printf("ttem 3185:soilyravln=%.3f,Iavln=%.3f\n",soil.yravln,y[I_AVLN]); //printf("ttem 3148: ston=%.3f, strn=%.3f, totn = %.3f\n", y[I_STON], y[I_STRN], y[I_STRN] + y[I_STON]);
  veg.yrstoren += y[I_STON];
  microbe.yrnh4 += y[I_NH4];
  microbe.yrno3 += y[I_NO3];

  // Update sum of annual water storage
  soil.yravlh2o += y[I_AVLW];
  soil.yrrgrndh2o += y[I_RGRW];
  soil.yrsnowpack += y[I_SNWPCK];
  soil.yrsgrndh2o += y[I_SGRW];
  soil.yrsmoist += y[I_SM];
  soil.yrpctp += y[I_PCTP];
  soil.yrvsm += y[I_VSM];

  // Update sum of annual carbon fluxes in natural ecosystems
  veg.yringpp += y[I_INGPP];   //printf("ttem3214: yringpp=%.3f,INGPP=%.3f\n",veg.yringpp,y[I_INGPP]);
  veg.yrgpp   += y[I_GPP];    //printf("ttem3215: yrgpp=%.3f,GPP=%.3f\n",veg.yrgpp,y[I_GPP]);
  veg.yrinnpp += y[I_INNPP];
  veg.yrnpp   += y[I_NPP];   //printf("ttem3217: yrnpp=%.3f,NPP=%.3f\n",veg.yrnpp,y[I_NPP]);
  veg.yrltrc  += y[I_LTRC];  //printf("ttem3218: yrltrc=%.3f,ltrc=%.3f\n",veg.yrltrc ,y[I_LTRC]);
  microbe.yrrh    += y[I_RH];
  yrnep   += y[I_NEP];

 // Update sum of annual nitrogen fluxes in natural ecosystems
  soil.yrnin   += y[I_NINP];
  veg.yrinnup += y[I_INNUP];
  veg.yrnup   += y[I_VNUP];
  veg.yrsup    += y[I_VSUP];
  veg.yrlup    += y[I_VLUP];
  veg.yrnmobil += y[I_VNMBL];
  veg.yrnrsorb += y[I_VNRSRB];
  veg.yrltrn  += y[I_LTRN]; //printf("ttem3215: yrltrn=%.3f,ltrn=%.3f\n",veg.yrltrn ,y[I_LTRN]);
  microbe.yrnmin  += y[I_NMIN];
  soil.yrnlost += y[I_NLST];
  soil.yrngas += y[I_NGAS];

   // Update sum of annual water fluxes
  atms.yrrain += y[I_RAIN];
  soil.yrrperc += y[I_RPERC];
  soil.yrrrun += y[I_RRUN];
  atms.yrsnowfall += y[I_SNWFAL];
  soil.yrsnowinf += y[I_SNWINF];
  soil.yrsperc += y[I_SPERC];
  soil.yrsrun += y[I_SRUN];
  atms.yrpet += y[I_PET];
  atms.yreet += y[I_EET];
  soil.yrh2oyld += y[I_WYLD];

  // Update sum of annual phenology in natural ecosystems
  veg.yrunleaf += y[I_UNRMLF];  //printf("ttem3239: yrunleaf=%.3f,yUNRMLF=%.3f\n",veg.yrunleaf,y[I_UNRMLF]);
  veg.yrleaf += y[I_LEAF];
  veg.yrfpc += veg.fpc[dm];

 // Update sum of annual carbon and nitrogen fluxes from
 //   agricultural conversion
  ag.yrconvrtC += y[I_CNVRTC];
  ag.yrconvrtN += y[I_CNVRTN];
  ag.yrsconvrtC += y[I_SCNVRTC];
  ag.yrsconvrtN += y[I_SCNVRTN];
  ag.yrslashC += y[I_SLASHC];
  ag.yrslashN += y[I_SLASHN];

 // Update sum of annual carbon and nitrogen fluxes from
 //   agricultural ecosystems
  ag.yrnppC += y[I_AGNPPC];
  ag.yrnppN += y[I_AGNPPN];
  ag.formPROD1.carbon += y[I_AGFPRDC];
  ag.formPROD1.nitrogen += y[I_AGFPRDN];
  ag.yrltrc += y[I_AGLTRC];
  ag.yrltrn += y[I_AGLTRN];
  ag.yrfertn += y[I_AGFRTN];

   // Update sum of annual integrated carbon fluxes
  yrcflux += y[I_CFLX];

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TTEM::setPrevState(double prevState[],double currentState[])
{
  for (int i = 0; i < NUMEQ; i++) { prevState[i] = currentState[i]; }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

int TTEM::stepyr(const int& dyr, const int& itype, int& intflag, double& tol)
{

  int dm;
  int mintflag;

  if (dyr == 0) { microbe.kd = microbe.kdc; }
  else
  {
    if (ag.state == 0 && ag.prvstate == 0)
    {
      microbe.kd = microbe.yrkd(nfeed, veg.yrltrc, veg.yrltrn, veg.cmnt); 
      ag.kd = microbe.kd;
      ag.natsoil = soil.org[CYCLE-1].carbon;
    }
    else
    {
      if (soil.org[CYCLE-1].carbon < ag.natsoil)
      {
        microbe.kd = ag.kd * soil.org[CYCLE-1].carbon / ag.natsoil;
      }
      else { microbe.kd = ag.kd; }
    }
  }

  // Reset annual fluxes to zero

  resetYrFlux();

  // Convert natural vegetation to agriculture

  if (ag.state == 1 && ag.prvstate == 0)
  {

    y[I_VEGC] = ag.vrespar[veg.cmnt] * veg.plant[CYCLE-1].carbon;
    y[I_STRN] = ag.vrespar[veg.cmnt] * veg.strctrl[CYCLE-1].nitrogen;
    y[I_STON] = ag.vrespar[veg.cmnt] * veg.labile[CYCLE-1].nitrogen;

    // Create PROD10 and PROD100 from conversion to agriculture

    ag.conversion(veg.cmnt, veg, soil);

  }

  // Revert to natural vegetation after cropland abandonment

  if (ag.state == 0 && ag.prvstate == 1)
  {
    veg.plant[CYCLE-1].carbon = y[I_VEGC];
    veg.strctrl[CYCLE-1].nitrogen = y[I_STRN];
    veg.labile[CYCLE-1].nitrogen = y[I_STON];
  }

  for (dm = 0; dm < CYCLE; dm++)
  {  
    if (initFlag == 1)
    {
      if (atms.tco2flag != 0) { atms.co2[dm] = atms.tco2[dyr][dm]; }
      if (atms.ttairflag != 0) { atms.tair[dm] = atms.ttair[dyr][dm]; } 
      if (atms.tprecflag != 0) { atms.prec[dm] = atms.tprec[dyr][dm]; } 
      if (ag.tlulcflag != 0)
      {
        if (ag.RAP0flag == 0)
        {
	  ag.potnpp[dm] = ag.tpotnpp[itype][dyr][dm];
        }
        else { ag.potnpp[dm] = 0.0; }
      }
    }

    // Get environmental conditions for month "dm"

     // calculate snowpack for the specific year

//    atms.pet[dm] = atms.petjh(atms.nirr[dm], atms.tair[dm], dm);
//    if (atms.pet[dm] <= 0.0) { atms.pet[dm] = 0.001; }

	 atms.precsplt(atms.prec[dm], atms.tair[dm], atms.rain[dm], atms.snowfall[dm]);
    switch (dm) {
      case 0:  atms.prevtair = atms.tair[CYCLE-1];
	       atms.prev2tair = atms.tair[CYCLE-2];
          break;
      case 1:  atms.prevtair = atms.tair[0];
	       atms.prev2tair = atms.tair[CYCLE-1];
	       break;
      default: atms.prevtair = atms.tair[dm-1];
	       atms.prev2tair = atms.tair[dm-2];
	       break;
      }
    soil.snowinf[dm] = soil.snowmelt(elev, atms.tair[dm], atms.prevtair, y[I_SNWPCK]);
    soil.snowpack[dm]= atms.snowfall[dm] - soil.snowinf[dm];
    if (soil.snowpack[dm] < 0.0) { soil.snowpack[dm] = 0.0; }

   if (tmflgs.stmflg == 1) {
     switch (dm) {
     case 0: sthermal.airt19= (atms.tair[CYCLE-1]+atms.tair[0])/2.0; // satisfy soil thermal model
             sthermal.airt29= atms.tair[0];
             sthermal.airt39= (atms.tair[0]+atms.tair[1])/2.0;
             break;
     case 11:sthermal.airt19= (atms.tair[10]+atms.tair[11])/2.0; // satisfy soil thermal model
             sthermal.airt29= atms.tair[11];
             sthermal.airt39= (atms.tair[11]+atms.tair[0])/2.0;
            break;
     default:sthermal.airt19= (atms.tair[dm-1]+atms.tair[dm])/2.0; // satisfy soil thermal model
             sthermal.airt29= atms.tair[dm];
             sthermal.airt39= (atms.tair[dm]+atms.tair[dm+1])/2.0;
            break;
         }

    switch (dm) {
     case 0: sthermal.hsnow19=(soil.snowpack[CYCLE-1]+soil.snowpack[0])/2.0/1000.0; // satisfy soil thermal model
             sthermal.hsnow29= soil.snowpack[0]/1000.0;
             sthermal.hsnow39= (soil.snowpack[1]+soil.snowpack[0])/2.0/1000.0;
             break;
     case 11: sthermal.hsnow19=(soil.snowpack[11]+soil.snowpack[10])/2.0/1000.0; // satisfy soil thermal model
              sthermal.hsnow29= soil.snowpack[11]/1000.0;
              sthermal.hsnow39= (soil.snowpack[11]+soil.snowpack[0])/2.0/1000.0;
            break;
     default: sthermal.hsnow19=(soil.snowpack[dm-1]+soil.snowpack[dm])/2.0/1000.0; // satisfy soil thermal model
              sthermal.hsnow29= soil.snowpack[dm]/1000.0;
              sthermal.hsnow39= (soil.snowpack[dm]+soil.snowpack[dm+1])/2.0/1000.0;
            break;
         }
//   cout << airt19 << " " << airt29 << " " << airt39 << " " << hsnow19 << " " << hsnow29 << " " << hsnow39 <<endl;

	sthermal.is9 = 10;
    sthermal.ism19 = 9;
    sthermal.smass9 = 0.f;

   int i;

	for (i=1;i <=210; i++) {
    sthermal.t9[i] = 0;
  	sthermal.xfa9[i] = -1e10f;
	sthermal.xfb9[i] = 0.f;
	sthermal.water9[i] = 1.f;
	sthermal.dx9[i] = 0.f;
	sthermal.weight9[i] = 0.f; }


   sthermal.calcu_cdsnow = 0.2; // constant for equilibrium run --- an assumpition
   sthermal.water2 = soil.pctp[dm]/100.0; // assume moss and organic have the same water content, 29/march/2001
   tmax = atms.mxttair[0];
   tmin = atms.mittair[0];
   tave = atms.avttair[0];
        //printf("ttem 3406: dyr = %i, tmax = %.3f, tmin = %.3f, tave = %.3f\n",dyr, tmax, tmin, tave);
  if (sthermal.soiltemp_(&sthermal.water2, &sthermal.calcu_cdsnow, &sthermal.airt19, &sthermal.airt29, &sthermal.airt39, &sthermal.hsnow19, &sthermal.hsnow29,
    &sthermal.hsnow39,sthermal.weight9, &sthermal.smass9, &sthermal.is9, &sthermal.ism19, sthermal.x9, sthermal.dx9,
     sthermal.xfb9, sthermal.xfa9,	sthermal.water9, sthermal.t9, &sthermal.tsoil,&sthermal.frontd,&sthermal.thawbegin,
    &sthermal.thawend,sthermal.diffsoilt, veg.cmnt, tmax, tmin, tave, atms.longmean) !=0) //tmax, tmin, tave, atms.longmean, soil.vsm1[dm],soil.vsm2[dm], soil.vsm3[dm], pthick_m
    { printf("bad tem"); 
//        getch();
};

   //      ++kswitch;

      atms.frontd[dm]=sthermal.frontd;  //printf("ttem 3414: frontd = %.3f\n", atms.frontd[dm]);
      atms.thawbe[dm]=sthermal.thawbegin;
      atms.thawend[dm]=sthermal.thawend;

      atms.tsoil[dm]=sthermal.tsoil;
      atms.dst5[dm]=sthermal.diffsoilt[0];
      atms.dst10[dm]=sthermal.diffsoilt[1];
      atms.dst20[dm]=sthermal.diffsoilt[2];
      atms.dst50[dm]=sthermal.diffsoilt[3];
      atms.dst100[dm]=sthermal.diffsoilt[4];
      atms.dst200[dm]=sthermal.diffsoilt[5];

 /*** end of calling soil thermal model ***/
  } // stmflg ===1

 // * calculate the thawing-proportion in early month and later fall

  if (atms.dst5[dm] ==0.0) atms.dst5[dm] =0.001;
 if (atms.dst5[dm-1] ==0.0) atms.dst5[dm-1] =0.002;
 if (atms.dst5[dm+1] ==0.0) atms.dst5[dm+1] =0.003;

 if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] > 0.0) && (atms.dst5[dm+1] > 0.0))
     veg.thawpercent[dm] = 1- (atms.dst5[dm] / (atms.dst5[dm] - atms.dst5[dm-1]));

 if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] < 0.0) && (atms.dst5[dm+1] < 0.0))
     veg.thawpercent[dm] = 0.0;

 if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] < 0.0) && (atms.dst5[dm+1] > 0.0))
     veg.thawpercent[dm] = 1- (atms.dst5[dm+1] / (atms.dst5[dm+1] - atms.dst5[dm]));

 if ((atms.dst5[dm-1] > 0.0) && (atms.dst5[dm] > 0.0) && (atms.dst5[dm+1] < 0.0))
     veg.thawpercent[dm] = atms.dst5[dm] / (atms.dst5[dm] - atms.dst5[dm+1]);

 if ((atms.dst5[dm-1] > 0.0) && (atms.dst5[dm] > 0.0) && (atms.dst5[dm+1] > 0.0))
     veg.thawpercent[dm] = 1.0;

 if ((atms.dst5[dm-1] > 0.0) && (atms.dst5[dm] < 0.0) && (atms.dst5[dm+1] < 0.0))
     veg.thawpercent[dm] = atms.dst5[dm-1] / (atms.dst5[dm-1] - atms.dst5[dm]);

  if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] > 0.0) && (atms.dst5[dm+1] < 0.0))
    veg.thawpercent[dm] = ((1- (atms.dst5[dm] / (atms.dst5[dm] - atms.dst5[dm-1])))
        + (atms.dst5[dm] / (atms.dst5[dm] - atms.dst5[dm+1])))/2.0;

 if ((atms.dst5[dm-1] > 0.0) && (atms.dst5[dm] < 0.0) && (atms.dst5[dm+1] > 0.0))
   veg.thawpercent[dm] = ((atms.dst5[dm-1] / (atms.dst5[dm-1] - atms.dst5[dm]))
       + ( 1- (atms.dst5[dm+1] / (atms.dst5[dm+1] - atms.dst5[dm]))))/2.0;

 if (veg.thawpercent[dm] > 1.0) veg.thawpercent[dm] = 1.0;
 if (veg.thawpercent[dm] < 0.0) veg.thawpercent[dm] = 0.0;


/*
 if (atms.dst5[dm] ==0.0) atms.dst5[dm] =0.001;
 if (atms.dst5[dm-1] ==0.0) atms.dst5[dm-1] =0.002;
 if (atms.dst5[dm+1] ==0.0) atms.dst5[dm+1] =0.003;

 if (atms.dst5[dm] > 0.0)
     {
      if (atms.dst5[dm-1] > 0.0)
         {
          if (atms.dst5[dm+1] > 0.0)
            {
             veg.thawpercent[dm] = 1.0; // summer months
             // To overcome wired JUNE phenomena, Q. Z. Sept. 03 /2001, delete the following dm ==5 sentences
     //        if (dm == 5) {
     //          if (atms.frontd[5] < 1.5) veg.thawpercent[5] = abs(atms.frontd[5] / 1.5);

//              veg.thawpercent[4] = abs(atms.dst5[3] / (atms.dst5[4] - atms.dst5[3])); // may
//              veg.thawpercent[dm] = 1.1 * veg.thawpercent[4]; // Assume Jun thawing has not reached the rooting depth yet
      //        }
            }
          else //atms.dst5[dm+1] < 0.0)
            {
         //   veg.thawpercent[dm] = atms.dst5[dm]/(atms.dst5[dm] - atms.dst5[dm+1]); // sept.
           veg.thawpercent[dm] = 1.0;
            }
         }
         else // atms.dst5[dm-1]<0.0
         {
           if (atms.dst5[dm+1] > 0.0)
           {
//            veg.thawpercent[dm] = abs(atms.dst5[dm-1] / (atms.dst5[dm] - atms.dst5[dm-1])); // may
            veg.thawpercent[dm] =1.0;

           }
           else // //atms.dst5[dm+1] < 0.0)
            {
             veg.thawpercent[dm] =0.0;  // rare to happen
            }
         }
     }
  else // (atms.dst5[dm] < 0.0)
   {
     if (atms.dst5[dm-1] > 0.0)
         {
          if (atms.dst5[dm+1] > 0.0)
            {
             veg.thawpercent[dm] = 1.0; // rare happen
            }
          else
            {
                if (atms.dst5[dm-1] - atms.dst5[dm] == 0)
                {
                   veg.thawpercent[dm] = 1.0;
                }
                else
                {
                   veg.thawpercent[dm] = atms.dst5[dm-1]/(atms.dst5[dm-1] - atms.dst5[dm]); // october
                }
//           veg.thawpercent[dm] = atms.dst5[dm-1]/(atms.dst5[dm-1] - atms.dst5[dm]); // October
//              veg.thawpercent[dm] = 1.0;
             }
         }
         else // atms.dst5[dm-1]<0.0
         {
           if (atms.dst5[dm+1] > 0.0)
           {
                if (atms.dst5[dm] - atms.dst5[dm-1] == 0)
                {
                   veg.thawpercent[dm] = 0.01;
                }
                else
                {
                   veg.thawpercent[dm] = abs(atms.dst5[dm] / (atms.dst5[dm] - atms.dst5[dm-1])); // april
                }
           }
           else // //atms.dst5[dm+1] < 0.0)
            {
             veg.thawpercent[dm] =0.0;  // rare to happen
            }
         }
      }
  // end of adding for thawing-frozen
  */

    getenviron(dm);
    //printf("ttem3574: stepyr: year = %i, month = %i\n", dyr, dm);
    mintflag = adapt(NUMEQ,y,tol,dm);
    
    //if (dyr == 200) {nuptake[dm] = veg.nuptake[dm];}     //YY for stable vegetation uptake 101722
    //if (dyr > 200) {veg.nuptake[dm]=nuptake[dm];}    
    
    if (mintflag == 1) { intflag = 1; }

    if (blackhol != 0) { qualcon[dyr][itype] = 10; }

    massbal(y,prevy);
    setPrevState(prevy,y);

    // Update carbon, nitrogen and water pools and fluxes from
    //  integrator results
    setMonth(dm, y);

    resetODEflux(y);
  }


  //  Update maximum EET, maximum PET, GPP optimum temperature (veg.topt),
  //  and maximum leaf cover (veg.prvleafmx) for the current year
  for (dm = 0; dm < CYCLE; dm++)
  {
    if (dm == 0)
    {
      atms.prveetmx = atms.eet[0];
      atms.prvpetmx = atms.pet[0];
      veg.topt   = atms.tair[0];
      veg.prvleafmx[veg.cmnt] = veg.unnormleaf[0];  
      //printf("ttem 3607: prvleafmx=%.3f,unnormleaf=%.3f\n",veg.prvleafmx[veg.cmnt],veg.unnormleaf[0]);
    }
    else
    {
      if (atms.eet[dm] > atms.prveetmx) { atms.prveetmx = atms.eet[dm]; }
      if (atms.pet[dm] > atms.prvpetmx) { atms.prvpetmx = atms.pet[dm]; }
      if (veg.aleaf[veg.cmnt] == 0.0 && veg.bleaf[veg.cmnt] == 0.0
         && veg.cleaf[veg.cmnt] == 1.0)
      {
	if (atms.tair[dm] > veg.topt) { veg.topt = atms.tair[dm]; }
      }
      else
      {
	if (veg.unnormleaf[dm] > veg.prvleafmx[veg.cmnt])
        {
	  veg.prvleafmx[veg.cmnt] = veg.unnormleaf[dm];   
    //printf("ttem 3622: prvleafmx=%.3f,unnormleaf=%.3f\n",veg.prvleafmx[veg.cmnt],veg.unnormleaf[dm]);
	  veg.topt = atms.tair[dm];
	}
      }
    }
  }

  // Update optimum temperature parameters for GPP

  if (veg.topt > veg.toptmax[veg.cmnt]) { veg.topt = veg.toptmax[veg.cmnt]; }
  if (veg.topt < veg.toptmin[veg.cmnt]) { veg.topt = veg.toptmin[veg.cmnt]; }


// Determine vegetation C/N parameter as a function of vegetation type,
// annual PET, annual EET, CO2 concentration

  veg.updateC2N(veg.cmnt, atms.yreet, atms.yrpet, atms.co2[11], atms.initco2);

  soil.yravlh2o /= 12.0;
  soil.yrrgrndh2o /= 12.0;
  soil.yrsnowpack /= 12.0;
  soil.yrsgrndh2o /= 12.0;
  soil.yrsmoist /= 12.0;
  soil.yrpctp /= 12.0;
  soil.yrvsm /= 12.0;

  atms.yrtsoil /= 12.0; // for soil temperature
  atms.yrfrontd /= 12.0; // for soil temperature
  atms.yrthawbegin /= 12.0; // for soil temperature
  atms.yrthawend /= 12.0; // for soil temperature


  veg.yrcarbon  /= 12.0 ;
  veg.yrnitrogen  /= 12.0;
  veg.yrstructn /= 12.0;

  if (veg.yrstructn != 0.0)
  {
    veg.yrc2n  = veg.yrcarbon / veg.yrstructn;
  }
//  else { veg.yrc2n = MISSING; }

  soil.yrorgc /= 12.0;
  soil.yrorgn /= 12.0;

  if (soil.yrorgn != 0.0)
  {
    soil.yrc2n = soil.yrorgc / soil.yrorgn;
  }
//  else { soil.yrc2n = MISSING; }

  soil.yravln  /= 12.0;
  veg.yrstoren /= 12.0;
  veg.yrunleaf /= 12.0;  //printf("ttem3679:yrunleaf=%.3f\n",veg.yrunleaf);
  veg.yrleaf /= 12.0;
  veg.yrfpc /= 12.0;
	microbe.yrnh4 /= 12.0;
	microbe.yrno3 /= 12.0;

  // y[2] changed to y[I_SOLC] by DWK on 20000130 and
  // y[3] changed to y[I_SOLN] by DWK on 20000130
  if (baseline == 1)
  {
    soil.yrnin = 0.0;
    soil.yrnlost = 0.0;
    if (y[I_SOLC]/microbe.cnsoil[veg.cmnt] >= y[I_SOLN])  
    {
      soil.yrnin = (y[I_SOLC]/microbe.cnsoil[veg.cmnt]) - y[I_SOLN];// + soil.yrngas ; 
    }
    else
    {
      soil.yrnlost = y[I_SOLN] - (y[I_SOLC]/microbe.cnsoil[veg.cmnt]);// - soil.yrngas;
    }
    y[I_SOLN] = y[I_SOLC]/microbe.cnsoil[veg.cmnt];  
  }
  //printf("ttem 3692: soln=%.3f,solc=%.3f,cnsoil=%.3f\n", y[I_SOLN], y[I_SOLC], microbe.cnsoil[veg.cmnt]);
  if (endeq > 0) { ++endeq; }

  return endeq;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */
// modified for Soil Thermal model

//int TTEM::transient(const int& dyr, const int& itype, double& tol)
int TTEM::transient(const int& dyr, const int& itype, double& tol, const int& RTIME)
{

  endeq = 0;

  if (atms.tco2flag == 1) { totyr = atms.co2year[dyr]; }
  else if (atms.ttairflag == 1) { totyr = atms.tairyear[dyr]; }
  else if (atms.tprecflag == 1) { totyr = atms.precyear[dyr]; }

  if (atms.ttairflag != 0) { atms.mxtair = atms.mxttair[dyr]; }
  if (atms.tprecflag != 0) { atms.yrprec = atms.yrtprec[dyr]; }

  if (ag.tlulcflag == 1)
  {
    ag.state = ag.tstate[dyr];
    ag.RAP = ag.tRAP[dyr];
  }

 // stepyr(dyr,itype, intflag, tol);
   transtepyr(dyr,itype, intflag, tol, RTIME);
   //printf("ttem3669");
  // Update annual agricultural product pools and fluxes

  ag.updateyr(dyr);

  if (totyr == startyr) { wrtyr = 0;}
  if (totyr > startyr) {++wrtyr; }

  return wrtyr;

};

// addition in order to
  int dm;
  int mintflag; //use RTIME for soil thermal model
int TTEM::transtepyr(const int& dyr, const int& itype, int& intflag, double& tol, const int& RTIME)
{
  if (dyr == 0) { microbe.kd = microbe.kdc; }
  else
  {
    if (ag.state == 0 && ag.prvstate == 0)
    {
      microbe.kd = microbe.yrkd(nfeed, veg.yrltrc, veg.yrltrn, veg.cmnt);
      ag.kd = microbe.kd;
      ag.natsoil = soil.org[CYCLE-1].carbon;
    }
    else
    {
      if (soil.org[CYCLE-1].carbon < ag.natsoil)
      {
        microbe.kd = ag.kd * soil.org[CYCLE-1].carbon / ag.natsoil;
      }
      else { microbe.kd = ag.kd; }
    }
  }

  // Reset annual fluxes to zero

  resetYrFlux();

  // Convert natural vegetation to agriculture

  if (ag.state == 1 && ag.prvstate == 0)
  {

    y[I_VEGC] = ag.vrespar[veg.cmnt] * veg.plant[CYCLE-1].carbon;
    y[I_STRN] = ag.vrespar[veg.cmnt] * veg.strctrl[CYCLE-1].nitrogen;
    y[I_STON] = ag.vrespar[veg.cmnt] * veg.labile[CYCLE-1].nitrogen;
    // Create PROD10 and PROD100 from conversion to agriculture

    ag.conversion(veg.cmnt, veg, soil);

  }
   
  // Revert to natural vegetation after cropland abandonment

  if (ag.state == 0 && ag.prvstate == 1)
  {
    veg.plant[CYCLE-1].carbon = y[I_VEGC];
    veg.strctrl[CYCLE-1].nitrogen = y[I_STRN];
    veg.labile[CYCLE-1].nitrogen = y[I_STON];
  }

       //calculate CDM
     double CDM;
     CDM = 0.0;
     for (dm = 0; dm < CYCLE; dm++)
     CDM = CDM + (10.0 - atms.ttair[dyr][dm]);
  for (dm = 0; dm < CYCLE; dm++)
  {
      //printf("ttem 3776: transient: year = %i, month = %i\n", dyr, dm);
      //printf("ttem 3777: initFlag = %i, tairflag = %.i, tprecflag = %i\n", initFlag, atms.ttairflag, atms.tprecflag);
    if (initFlag == 1)
    {
      if (atms.tco2flag != 0) { atms.co2[dm] = atms.tco2[dyr][dm]; }
      if (atms.ttairflag != 0) { atms.tair[dm] = atms.ttair[dyr][dm]; } 
      if (atms.tprecflag != 0) { atms.prec[dm] = atms.tprec[dyr][dm]; } 

      if (ag.tlulcflag != 0)
      {
        if (ag.RAP0flag == 0)
        {
         ag.potnpp[dm] = ag.tpotnpp[itype][dyr][dm];
        }
        else { ag.potnpp[dm] = 0.0; }
      }
    }

    // Get environmental conditions for month "dm"
  // calculate snowpack for the specific year

//  for (dm = 0; dm < CYCLE; dm++) {
 //    atms.pet[dm] = atms.petjh(atms.nirr[dm], atms.tair[dm], dm);
 //   if (atms.pet[dm] <= 0.0) { atms.pet[dm] = 0.001; }
	 atms.precsplt(atms.prec[dm], atms.tair[dm], atms.rain[dm], atms.snowfall[dm]);
     switch (dm) {
      case 0:  atms.prevtair = atms.tair[CYCLE-1];
	       atms.prev2tair = atms.tair[CYCLE-2];
          break;
      case 1:  atms.prevtair = atms.tair[0];
	       atms.prev2tair = atms.tair[CYCLE-1];
	       break;
      default: atms.prevtair = atms.tair[dm-1];
	       atms.prev2tair = atms.tair[dm-2];
	       break;
      }
    soil.snowinf[dm] = soil.snowmelt(elev, atms.tair[dm], atms.prevtair, y[I_SNWPCK]);
    soil.snowpack[dm]= atms.snowfall[dm] - soil.snowinf[dm];
     if (soil.snowpack[dm] < 0.0) { soil.snowpack[dm] = 0.0; }


  if (tmflgs.stmflg==1) {
// call soil thermal subroutine
    switch (dm) {
     case 0: if (dyr==0) {
             sthermal.airt19= atms.ttair[dyr][0];; // satisfy soil thermal model
             sthermal.airt29= (atms.ttair[dyr][1]+atms.ttair[dyr][0])/2.0;
             sthermal.airt39= (atms.ttair[dyr][1]+atms.ttair[dyr][2])/2.0; }
             else {
             sthermal.airt19= (atms.ttair[dyr-1][11]+atms.ttair[dyr][0])/2.0; // satisfy soil thermal model
             sthermal.airt29= atms.ttair[dyr][0];
             sthermal.airt39= (atms.ttair[dyr][0]+atms.ttair[dyr][1])/2.0;
              }
             break;
     case 11:if ( dyr<RTIME-3) {
             sthermal.airt19= (atms.ttair[dyr][CYCLE-3]+atms.ttair[dyr][10])/2.0; // satisfy soil thermal model
             sthermal.airt29= atms.ttair[dyr][11];
             sthermal.airt39= (atms.ttair[dyr][11]+atms.ttair[dyr+1][0])/2.0; }
             else {
             sthermal.airt19= (atms.ttair[dyr][CYCLE-2]+atms.ttair[dyr][11])/2.0; // satisfy soil thermal model
             sthermal.airt29= atms.ttair[dyr][11];
             sthermal.airt39= (atms.ttair[dyr][11]+atms.ttair[dyr][0])/2.0; }
             break;
     default: sthermal.airt19= (atms.ttair[dyr][dm-1]+atms.ttair[dyr][dm])/2.0; // satisfy soil thermal model
             sthermal.airt29= atms.ttair[dyr][dm];
             sthermal.airt39= (atms.ttair[dyr][dm]+atms.ttair[dyr][dm+1])/2.0;
             break;
         }

     // using soil.snowpack, which is calculated from the WBM, to drive soil thermal model
     // convert mm to m for snow depth by dividing 1000, also assume the snow density is changed according to the snow classification, March,28, 2001, by Q. Z.
     // According to the vegetation type determine the wind speed, calculate cooling degree month (CDM) and daily averaged precipitation to determine the snow classication
     // snow classification determine the snow density, and therefore the snow depth, and snow thermal conductivity

     // determine the wind speed level
     int wind_speed; // high = 1; low =0.0

     switch (veg.cmnt) {
         case 1: wind_speed =1; break;
         case 2: wind_speed =1; break;
         case 3: wind_speed =0; break;
         case 4: wind_speed =0; break;
         case 5: wind_speed =0; break;
         case 6: wind_speed =1; break;
         case 7: wind_speed =1; break;
         case 8: wind_speed =1; break;
         case 9: wind_speed =0; break;
         case 10: wind_speed =1; break;
         case 11: wind_speed =0; break;
         case 12: wind_speed =1; break;
        }

     // determine the snow classes
     int snow_class;

     if (CDM < 50.0) snow_class = 6; // class Ephemeral
     else
       {
        if (CDM > 125.0) {
                         if (atms.prec[dm] / atms.daze[dm] >2.0)
                           {
                           if (wind_speed ==1) snow_class = 1;
                           else  snow_class = 2;
                           }
                         else
                            {
                           if (wind_speed ==0) snow_class = 2;
                           else snow_class = 1;
                            }
                         }

        else {
              if (atms.prec[dm] / atms.daze[dm] <2.0)
              {
               if (wind_speed ==1) snow_class = 4;
               else snow_class = 3;
              }
              else
               {
                if (wind_speed ==0) snow_class = 5;
                else snow_class = 5;
               }

            }
        }

// determine the snow denisty and snow thermal conductivity
  double snow_dens; // g mm-3

   switch (snow_class) {
   case 1: snow_dens = 0.28;  break;
   case 2: snow_dens = 0.225;  break;
   case 3: snow_dens = 0.25;  break;
   case 4: snow_dens = 0.25;  break;
   case 5: snow_dens = 0.30;  break;
   case 6: snow_dens = 0.35;  break;
   }

  sthermal.calcu_cdsnow = powl(10, (2.650 * snow_dens - 1.652));

  sthermal.water2 = soil.pctp[dm]/100.0; // assume moss and organic have the same water content, 29/march/2001

  // end of the program of snow classification
    //printf("nodepo 3918: snowpack = %.3f, snowdens = %.3f, tair = %.3f, prec = %.3f, rain = %.3f, snow = %.3f\n", soil.snowpack[dm], snow_dens, atms.ttair[dyr][dm], atms.prec[dm], atms.rain[dm], atms.snowfall[dm]);
    switch (dm) {
     case 0: sthermal.hsnow19=(soil.snowpack[CYCLE-1]+soil.snowpack[0])/2.0/1000.0/snow_dens; // satisfy soil thermal model
             sthermal.hsnow29= soil.snowpack[0]/1000.0/snow_dens;
             sthermal.hsnow39= (soil.snowpack[1]+soil.snowpack[0])/2.0/1000.0/snow_dens;
             break;
     case 11: sthermal.hsnow19=(soil.snowpack[11]+soil.snowpack[10])/2.0/1000.0/snow_dens; // satisfy soil thermal model
              sthermal.hsnow29= soil.snowpack[11]/1000.0/snow_dens;
              sthermal.hsnow39= (soil.snowpack[11]+soil.snowpack[0])/2.0/1000.0/snow_dens;
            break;
     default: sthermal.hsnow19=(soil.snowpack[dm-1]+soil.snowpack[dm])/2.0/1000.0/snow_dens; // satisfy soil thermal model
              sthermal.hsnow29= soil.snowpack[dm]/1000.0/snow_dens;
              sthermal.hsnow39= (soil.snowpack[dm]+soil.snowpack[dm+1])/2.0/1000.0/snow_dens;
            break;
               }

    int i;

    sthermal.is9 = 10;
    sthermal.ism19 = 9;
    sthermal.smass9 = 0.f;

	for (i=1;i <=210; i++) {
    sthermal.t9[i] = 0;
  	sthermal.xfa9[i] = -1e10f;
	sthermal.xfb9[i] = 0.f;
	sthermal.water9[i] = 1.f;
	sthermal.dx9[i] = 0.f;
	sthermal.weight9[i] = 0.f; }

// Need to change the soil moisture (water content for organic layer), snow thermal conductivity based on the snow classification ??
    tmax = atms.mxttair[dyr];
    tmin = atms.mittair[dyr];
    tave = atms.avttair[dyr];
        //printf("ttem 3954: dyr = %i, tmax = %.3f, tmin = %.3f, tave = %.3f\n",dyr, tmax, tmin, tave);
        //printf("nodepo ttem 3953: tmax = %.3f, tmin = %.3f, tave = %.3f, longmean = %.3f\n", tmax, tmin, tave, atms.longmean);
         // printf("nodepo ttem 3954: water2 = %.3f, cdsnow = %.3f, airt19 = %.3f, arit29 = %.3f, airt39 = %.3f\n", sthermal.water2, sthermal.calcu_cdsnow, sthermal.airt19, sthermal.airt29, sthermal.airt39);
      //printf("nodepo ttem 3955: hsnow19 = %.3f, hsnow29 = %.3f, hsnow39 = %.3f\n", sthermal.hsnow19, sthermal.hsnow29, sthermal.hsnow39);         
     if (sthermal.soiltemp_(&sthermal.water2, &sthermal.calcu_cdsnow, &sthermal.airt19, &sthermal.airt29, &sthermal.airt39, &sthermal.hsnow19, &sthermal.hsnow29, &sthermal.hsnow39,sthermal.weight9, &sthermal.smass9, &sthermal.is9, &sthermal.ism19, sthermal.x9, sthermal.dx9, sthermal.xfb9, sthermal.xfa9,
		sthermal.water9, sthermal.t9, &sthermal.tsoil,&sthermal.frontd,&sthermal.thawbegin,&sthermal.thawend,sthermal.diffsoilt,veg.cmnt, tmax, tmin, tave, atms.longmean) !=0) //,tmax, tmin, tave, atms.longmean, soil.vsm1[dm],soil.vsm2[dm], soil.vsm3[dm], pthick_m
      { printf("bad tem");
  //      getch();
        };
    //  ++ integer (kswitch);
      atms.frontd[dm]=sthermal.frontd;
      //printf("ttem 3963: dm=%i, atms.frontd = %.3f, sthermal.frontd=%.3f\n", dm, atms.frontd[dm], sthermal.frontd);
      atms.thawbe[dm]=sthermal.thawbegin;
      atms.thawend[dm]=sthermal.thawend;

      atms.tsoil[dm]=sthermal.tsoil;
      atms.dst5[dm]=sthermal.diffsoilt[0];
      atms.dst10[dm]=sthermal.diffsoilt[1];
      atms.dst20[dm]=sthermal.diffsoilt[2];
      atms.dst50[dm]=sthermal.diffsoilt[3];
      atms.dst100[dm]=sthermal.diffsoilt[4];
      atms.dst200[dm]=sthermal.diffsoilt[5];
     } //stmflg ==1
   /*** end of calling soil thermal model ***/
    //printf("ttem 3933: dm =%i, atms.tsoil=%.3f,atms.dst5=%.3f, atms.dst10=%.3f, atms.dst20=%.3f, atms.dst50%.3f\n", dm, atms.tsoil[dm],atms.dst5[dm], atms.dst10[dm],atms.dst20[dm],atms.dst50[dm]);
   // * calculate the thawing-proportion in early month and later fall, 28/Nov/2001

 if (atms.dst5[dm] ==0.0) atms.dst5[dm] =0.001;
 if (atms.dst5[dm-1] ==0.0) atms.dst5[dm-1] =0.002;
 if (atms.dst5[dm+1] ==0.0) atms.dst5[dm+1] =0.003;

 if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] > 0.0) && (atms.dst5[dm+1] > 0.0))
     veg.thawpercent[dm] = 1- (atms.dst5[dm] / (atms.dst5[dm] - atms.dst5[dm-1]));

 if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] < 0.0) && (atms.dst5[dm+1] < 0.0))
     veg.thawpercent[dm] = 0.0;

 if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] < 0.0) && (atms.dst5[dm+1] > 0.0))
     veg.thawpercent[dm] = 1- (atms.dst5[dm+1] / (atms.dst5[dm+1] - atms.dst5[dm]));

 if ((atms.dst5[dm-1] > 0.0) && (atms.dst5[dm] > 0.0) && (atms.dst5[dm+1] < 0.0))
     veg.thawpercent[dm] = atms.dst5[dm] / (atms.dst5[dm] - atms.dst5[dm+1]);

 if ((atms.dst5[dm-1] > 0.0) && (atms.dst5[dm] > 0.0) && (atms.dst5[dm+1] > 0.0))
     veg.thawpercent[dm] = 1.0;

 if ((atms.dst5[dm-1] > 0.0) && (atms.dst5[dm] < 0.0) && (atms.dst5[dm+1] < 0.0))
     veg.thawpercent[dm] = atms.dst5[dm-1] / (atms.dst5[dm-1] - atms.dst5[dm]);

 if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] > 0.0) && (atms.dst5[dm+1] < 0.0))
    veg.thawpercent[dm] = ((1- (atms.dst5[dm] / (atms.dst5[dm] - atms.dst5[dm-1])))
        + (atms.dst5[dm] / (atms.dst5[dm] - atms.dst5[dm+1])))/2.0;

 if ((atms.dst5[dm-1] > 0.0) && (atms.dst5[dm] < 0.0) && (atms.dst5[dm+1] > 0.0))
   veg.thawpercent[dm] = ((atms.dst5[dm-1] / (atms.dst5[dm-1] - atms.dst5[dm]))
       + ( 1- (atms.dst5[dm+1] / (atms.dst5[dm+1] - atms.dst5[dm]))))/2.0;

 if (veg.thawpercent[dm] > 1.0) veg.thawpercent[dm] = 1.0;
 if (veg.thawpercent[dm] < 0.0) veg.thawpercent[dm] = 0.0;

// mark following; added above new algorithm

	// calculate maximum front depth of current year - YYuan, 12/28/22
	//for (dm = 0; dm < CYCLE; dm++)
	//{
		if (dm == 0){atms.frontdmx = atms.frontd[0];} //c frontdmx of last year
		else
		{
      if (atms.frontd[dm] < atms.frontdmx) // atms.frontd is negitave. 
			{
        atms.frontdmx = atms.frontd[dm];				
			}
		}
    //printf("ttem4015: dm=%i, frontd=%.3f, frontd0=%.3f, frontdmx= %.3f\n",dm, atms.frontd[dm], atms.frontd[0], atms.frontdmx);     

    if (atms.frontdmx < -13) {atms.frontdmx = -13;}
    if (dyr <= 215 ) {
    if (atms.frontdmx < atms.frontdymx ) {
        atms.frontdymx = atms.frontdmx;
        atms.frontddiff = 0;}
    }
    //printf("ttem4026: frontdmx= %.3f,frontdymx= %.3f\n",atms.frontdmx,atms.frontdymx);
    if (dyr > 216 ) {
      if (atms.frontdmx < atms.frontdymx ) {//printf("ttem4028: frontdmx= %.3f,frontdymx= %.3f\n",atms.frontdmx,atms.frontdymx);
        atms.frontddiff = atms.frontdymx - atms.frontdmx; 
        //printf("ttem4031: frontdmx= %.5f,frontdymx= %.5f,frontddiff= %.5f\n",atms.frontdmx,atms.frontdymx,atms.frontddiff);
        atms.frontdymx = atms.frontdmx;
        //printf("ttem4033: frontdymx= %.3f,frontmx= %.3f\n",atms.frontdymx,atms.frontdmx);
        }
      else {atms.frontddiff = 0;} //printf("ttem4035: frontddiff= %.3f\n",atms.frontddiff);
    } 
    
    /*if (atms.frontddiff > 0){
       if (atms.frontdmx < 0 && atms.frontdmx > -0.3 ) {cperma = atms.frontddiff * 33898 };
       if (atms.frontdmx <= -0.3 && atms.frontdmx > -1.0 ) 
          if (atms.frontdymx > -0.3) {cperma = (atms.frontddiff-0.3) * 23199 + (0.3-atms.frontdymx)*33898};
          else {cperma = (atms.frontddiff-0.3) * 23199 };
       if (atms.frontdmx <= -1.0 && atms.frontdmx > -3.0 ) 
          if (atms.frontdymx > -0.3) {cperma =  (atms.frontddiff-1) * 14056 + (1-0.3) * 23199 + (0.3-atms.frontdymx)*33898};
          else if (atms.frontdymx > -1.0 && atms.frontdymx )
    }*/
    //if (atms.frontddiff > 0 ) {atms.frontddiff2 = atms.frontddiff ;}
    //else {atms.frontddiff2=0;}
    //printf("ttem4027: frontddiff= %.3f,frontddiff2= %.3f\n",atms.frontddiff,atms.frontddiff2);
	//} 
 
/*
 if (atms.dst5[dm] ==0.0) atms.dst5[dm] =0.001;
 if (atms.dst5[dm-1] ==0.0) atms.dst5[dm-1] =0.002;
 if (atms.dst5[dm+1] ==0.0) atms.dst5[dm+1] =0.003;

 if (atms.dst5[dm] > 0.0)
     {
      if (atms.dst5[dm-1] > 0.0)
         {
          if (atms.dst5[dm+1] > 0.0)
            {
             veg.thawpercent[dm] = 1.0; // summer months
             // To overcome wired JUNE phenomena, Q. Z. Sept. 03 /2001, delete the following dm ==5 sentences
     //        if (dm == 5) {
     //          if (atms.frontd[5] < 1.5) veg.thawpercent[5] = abs(atms.frontd[5] / 1.5);

//              veg.thawpercent[4] = abs(atms.dst5[3] / (atms.dst5[4] - atms.dst5[3])); // may
//              veg.thawpercent[dm] = 1.1 * veg.thawpercent[4]; // Assume Jun thawing has not reached the rooting depth yet
      //        }
            }
          else //atms.dst5[dm+1] < 0.0)
            {
         //   veg.thawpercent[dm] = atms.dst5[dm]/(atms.dst5[dm] - atms.dst5[dm+1]); // sept.
           veg.thawpercent[dm] = 1.0;
            }
         }
         else // atms.dst5[dm-1]<0.0
         {
           if (atms.dst5[dm+1] > 0.0)
           {
//            veg.thawpercent[dm] = abs(atms.dst5[dm-1] / (atms.dst5[dm] - atms.dst5[dm-1])); // may
            veg.thawpercent[dm] =1.0;

           }
           else // //atms.dst5[dm+1] < 0.0)
            {
             veg.thawpercent[dm] =0.0;  // rare to happen
            }
         }
     }
  else // (atms.dst5[dm] < 0.0)
   {
     if (atms.dst5[dm-1] > 0.0)
         {
          if (atms.dst5[dm+1] > 0.0)
            {
             veg.thawpercent[dm] = 1.0; // rare happen
            }
          else
            {
                if (atms.dst5[dm-1] - atms.dst5[dm] == 0)
                {
                   veg.thawpercent[dm] = 1.0;
                }
                else
                {
                   veg.thawpercent[dm] = atms.dst5[dm-1]/(atms.dst5[dm-1] - atms.dst5[dm]); // october
                }
//           veg.thawpercent[dm] = atms.dst5[dm-1]/(atms.dst5[dm-1] - atms.dst5[dm]); // October
//              veg.thawpercent[dm] = 1.0;
             }
         }
         else // atms.dst5[dm-1]<0.0
         {
           if (atms.dst5[dm+1] > 0.0)
           {
                if (atms.dst5[dm] - atms.dst5[dm-1] == 0)
                {
                   veg.thawpercent[dm] = 0.01;
                }
                else
                {
                   veg.thawpercent[dm] = abs(atms.dst5[dm] / (atms.dst5[dm] - atms.dst5[dm-1])); // april
                }
           }
           else // //atms.dst5[dm+1] < 0.0)
            {
             veg.thawpercent[dm] =0.0;  // rare to happen
            }
         }
      }
  */
    getenviron(dm);
    //printf("ttem4106 transient: year=%i, month = %i\n",dyr,dm);
    
    //if (dyr == 262 && dm == 11 ) {printf("n2o = %.2f\n", microbe.N2OF[dm]*1000);}   
    //if (dyr == 263 && dm == 0 ) {printf("n2o = %.2f\n", microbe.N2OF[dm]*1000);}    
    //if (dyr == 263 && dm == 1 ) {printf("n2o = %.2f\n", microbe.N2OF[dm]*1000);}
    //if (dyr == 263 && dm == 2 ) {printf("n2o = %.2f\n", microbe.N2OF[dm]*1000);}   
    //if (dyr == 263 && dm == 3 ) {printf("n2o = %.2f\n", microbe.N2OF[dm]*1000);}                  
    //if (dyr == 263 && dm == 4 ) {printf("n2o = %.2f\n", microbe.N2OF[dm]*1000);}
    //if (dyr == 263 && dm == 5 ) {printf("n2o = %.2f\n", microbe.N2OF[dm]*1000);}    
    //if (dyr == 263 && dm == 6 ) {printf("n2o = %.2f\n", microbe.N2OF[dm]*1000);}    
   // if (dyr == 263 && dm == 7 ) {printf("n2o = %.2f\n", microbe.N2OF[dm]*1000);}    
    //if (dyr == 263 && dm == 8 ) {printf("n2o = %.2f\n", microbe.N2OF[dm]*1000);}    
    //if (dyr == 263 && dm == 9 ) {printf("n2o = %.2f\n", microbe.N2OF[dm]*1000);}    
    //if (dyr == 263 && dm == 10 ) {printf("n2o = %.2f\n", microbe.N2OF[dm]*1000);}    
    //if (dyr == 263 && dm == 11 ) {printf("n2o = %.2f\n", microbe.N2OF[dm]*1000);}     
                    


        
    mintflag = adapt(NUMEQ,y,tol,dm);


    if (mintflag == 1) { intflag = 1; }

    if (blackhol != 0) { qualcon[dyr][itype] = 10; }

    massbal(y,prevy);

    setPrevState(prevy,y);

    // Update carbon, nitrogen and water pools and fluxes from
    //  integrator results

    setMonth(dm, y);
    resetODEflux(y);

  }

  //  Update maximum EET, maximum PET, GPP optimum temperature (veg.topt),
  //  and maximum leaf cover (veg.prvleafmx) for the current year
  for (dm = 0; dm < CYCLE; dm++)
  {
    if (dm == 0)
    {
      atms.prveetmx = atms.eet[0];
      atms.prvpetmx = atms.pet[0];
      veg.topt   = atms.tair[0];
      veg.prvleafmx[veg.cmnt] = veg.unnormleaf[0];
    }
    else
    {
      if (atms.eet[dm] > atms.prveetmx) { atms.prveetmx = atms.eet[dm]; }
      if (atms.pet[dm] > atms.prvpetmx) { atms.prvpetmx = atms.pet[dm]; }
      if (veg.aleaf[veg.cmnt] == 0.0 && veg.bleaf[veg.cmnt] == 0.0
         && veg.cleaf[veg.cmnt] == 1.0)
      {
	if (atms.tair[dm] > veg.topt) { veg.topt = atms.tair[dm]; }
      }
      else
      {
	if (veg.unnormleaf[dm] > veg.prvleafmx[veg.cmnt])
        {
	  veg.prvleafmx[veg.cmnt] = veg.unnormleaf[dm];
	  veg.topt = atms.tair[dm];
	}
      }
    }
  }

  // Update optimum temperature parameters for GPP

  if (veg.topt > veg.toptmax[veg.cmnt]) { veg.topt = veg.toptmax[veg.cmnt]; }
  if (veg.topt < veg.toptmin[veg.cmnt]) { veg.topt = veg.toptmin[veg.cmnt]; }


// Determine vegetation C/N parameter as a function of vegetation type,
// annual PET, annual EET, CO2 concentration

  veg.updateC2N(veg.cmnt, atms.yreet, atms.yrpet, atms.co2[11], atms.initco2);

  soil.yravlh2o /= 12.0;
  soil.yrrgrndh2o /= 12.0;
  soil.yrsnowpack /= 12.0;
  soil.yrsgrndh2o /= 12.0;
  soil.yrsmoist /= 12.0;
  soil.yrpctp /= 12.0;
  soil.yrvsm /= 12.0;

  atms.yrtsoil /= 12.0; // for soil temperature
  atms.yrfrontd /= 12.0; // for soil temperature
  atms.yrthawbegin /= 12.0; // for soil temperature
  atms.yrthawend /= 12.0; // for soil temperature

  veg.yrcarbon  /= 12.0;
  veg.yrnitrogen  /= 12.0;
  veg.yrstructn /= 12.0;

  if (veg.yrstructn != 0.0)
  {
    veg.yrc2n  = veg.yrcarbon / veg.yrstructn;
  }
//  else { veg.yrc2n = MISSING; }

  soil.yrorgc /= 12.0;
  soil.yrorgn /= 12.0;

  if (soil.yrorgn != 0.0)
  {
    soil.yrc2n = soil.yrorgc / soil.yrorgn;
  }
//  else { soil.yrc2n = MISSING; }

  soil.yravln  /= 12.0; 
  veg.yrstoren /= 12.0;
  veg.yrunleaf /= 12.0;
  veg.yrleaf /= 12.0;
  veg.yrfpc /= 12.0;
	microbe.yrnh4 /= 12.0;
	microbe.yrno3 /= 12.0;
	

  // y[2] changed to y[I_SOLC] by DWK on 20000130 and
  // y[3] changed to y[I_SOLN] by DWK on 20000130
  if (baseline == 1)
  {
    soil.yrnin = 0.0;
    soil.yrnlost = 0.0;
    if (y[I_SOLC]/microbe.cnsoil[veg.cmnt] >= y[I_SOLN])
    {
      soil.yrnin = (y[I_SOLC]/microbe.cnsoil[veg.cmnt]) - y[I_SOLN];// + soil.yrngas;
    }
    else
    {
      soil.yrnlost = y[I_SOLN] - (y[I_SOLC]/microbe.cnsoil[veg.cmnt]);// - soil.yrngas;
    }
    y[I_SOLN] = y[I_SOLC]/microbe.cnsoil[veg.cmnt];
  }

  if (endeq > 0) { ++endeq; }

  return endeq;

};


//added by YYL to determine the number of days per month
int TTEM::daysnumber(const int& whichmonth, const int& whichyear)
{
  int days;

    switch (whichmonth)
       {
        case 0: days =31; break;
        case 2: days =31; break;
        case 4: days =31; break;
        case 6: days =31; break;
        case 7: days =31; break;
        case 9: days =31; break;
        case 11: days =31; break;
        case 3: days =30; break;
        case 5: days =30; break;
        case 8: days =30; break;
        case 10: days =30; break;
        case 1:
               if ( fmod((whichyear - 1992), 4) == 0 ) days = 29;
               else days =28; break;
        default: return (-1);
       }

  return days;
  
};

// month to daily -- climate
 int TTEM::monthTOdaily ( double airtemp[], double preci[], double vap[] )
 {
 const int NY = 1;
 const int NM = 12*NY;
 const int NJD=NY*365;
 const int NJ=(NY+1)*365;

 int MMONTH[NM+1];
// double TEMP[NM+1] = {1000,-18.75,-15.50,-3.68,3.91,11.06,12.83,15.54,15.20,8.72,2.33,-8.71,-12.51};
// double PREC[NM+1] = {1000,14.20,1.40,5.40,21.80,10.80,146.60,161.60,37.00,46.20,11.00,23.00,6.80};
 double TEMP[13];
 double PREC[13];
 double VAPO[13];

 double TMP1[13+1],TMP15[NM+1];
 double VAP1[13+1],VAP15[NM+1];

 double B[NM+1], T[NM+1], S[NM+1], RB[NM+1];
 double DURB[NM+1], DURT[NM+1], DURS[NM+1];
 double TEM[NJD+1], RAININTE[NJD+1], RAINDUR[NJD+1];

 double VAPOR[NJD+1];

 int RAINN;
 int JDAY[NJD+1], JJD[NJ+1], DRYD[NJD+1], RAINNUMB[NJD+1];
 int DX[12+1] = {100,0,1,-1,0,0,1,1,2,3,3,4,4};  // the first value 100 is just the space holder
 int DD[12+1] = {100,1,-2,1,0,1,0,1,1,0,1,0,1};  // the first value 100 is just the space holder
 int MONTH, DAY, DT;

 int I, JJ, K, L;

 //  double TMO, PMO;
 double RT, RS, R, TPREC, RTN, ATEMP;
 int NI, NS, NSY, NN;

 double BB, TT;
 int KTT, KTD, KDD, KKTD;

// ofstream fout;

  for (K = 1; K <= NM; K++) {
        TEMP[K] = airtemp[K];
        PREC[K] = preci[K];
        VAPO[K] = vap[K];
  }

//   The integration date begins at MONTH/DAY. i.e. Jan. 1.
      MONTH=1;
      DAY=1;

      RT=1.778;
      RS=0.635;
      R=0.5;
      TPREC=0.0;
      RTN=0.0;
      ATEMP=0.0;

  for (I = 1; I <= NY; I++) {
      for (K=(I-1)*12+1; K <= I*12; K++) {
          MMONTH[K]=MONTH+K-1-(I-1)*12;
          if (MMONTH[K] > 12) MMONTH[K]=MMONTH[K]-12;
      }
   }

  for (K=1; K <= NM; K++) {

      PREC[K]=PREC[K]/10.0/2.54; //comverse mm to cm, then to in.
      DURT[K]=RT/R;
      DURS[K]=RS/R;

//   Case 1, TEMP<0.
      if (TEMP[K] <= 0.0) {
        if (PREC[K] <= 1.0) {
          B[K]=1.0;
          T[K]=0.0;
          S[K]=1.0;
          }
        else {
          B[K]=1.0;
          T[K]=1.0;
          S[K]=1.0;
        }
        RB[K]=(PREC[K]*2.54-RS*S[K]-RT*T[K])/B[K];
        DURB[K]=RB[K]/R;
      }

//   Case 2, PREC<1.0 inch.
     else if (PREC[K] <= 1.0) {
          B[K]=1.0;
          T[K]=0.0;
          S[K]=1.0;
          RB[K]=(PREC[K]*2.54-RS*S[K]-RT*T[K])/B[K];
          DURB[K]=RB[K]/R;
      }

//   Case 3, 1.0<PREC<2.5 inches.
      else if ((2.5 >= PREC[K]) && (PREC[K] > 1.0)) {
           B[K]=1.0;
           T[K]=1.0;
           S[K]=1.0;
           RB[K]=(PREC[K]*2.54-RS*S[K]-RT*T[K])/B[K];
           DURB[K]=RB[K]/R;
      }

//   Case 4, 2.5<PREC<4.0 inches.
      else if ((4.0 >= PREC[K]) && (PREC[K] > 2.5)) {
           B[K]=1.0;
           S[K]=4.0;
           if (PREC[K] < 3.7)
             T[K]=1.0;
           else
             T[K]=2.0;
           RB[K] = (PREC[K]*2.54-RS*S[K]-RT*T[K])/B[K];
           DURB[K] = RB[K]/R;
      }

//   Case 5, 4.0<PREC<5.0 inches.
      else if ((5.0 >= PREC[K]) && (PREC[K] > 4.0)) {
           B[K]=1.0;
           S[K]=4.0;
           if (PREC[K] < 4.43)
             T[K]=1.0;
           else
             T[K]=2.0;
           RB[K]=(PREC[K]*2.54-RS*S[K]-RT*T[K])/B[K];
           DURB[K]=RB[K]/R;
      }

//   Case 6, 5.0<PREC<7.0 inches.
      else if ((7.0 >= PREC[K]) && (PREC[K] > 5.0)) {
           B[K]=2.0;
           S[K]=4.0;
           if (PREC[K] < 5.65)
             T[K]=1.0;
           else
             T[K]=2.0;
           RB[K]= (PREC[K]*2.54-RS*S[K]-RT*T[K])/B[K];
           DURB[K]=RB[K]/R;
      }

//   Case 7, 7.0<PREC<9.0 inches.
      else if ((9.0 >= PREC[K]) && (PREC[K] > 7.0)) {
           B[K]=2.0;
           S[K]=6.0;
           if (PREC[K] < 8.21)
             T[K]=3.0;
           else
             T[K]=4.0;
           RB[K]=(PREC[K]*2.54-RS*S[K]-RT*T[K])/B[K];
           DURB[K]=RB[K]/R;
      }

//   Case 8, 9.0<PREC<11.0 inches.
      else if ((11.0 >= PREC[K]) && (PREC[K] > 9.0)) {
           B[K]=3.0;
           S[K]=6.0;
           if (PREC[K] < 10.0)
              T[K]=4.0;
           else
              T[K]=5.0;
           RB[K]=(PREC[K]*2.54-RS*S[K]-RT*T[K])/B[K];
           DURB[K]=RB[K]/R;
      }

//   Case 9, PREC>11.0 inches.
      else if (PREC[K] > 11.0) {
           B[K]=4.0;
           S[K]=7.0;
           if (PREC[K] < 13.0)
             T[K]=4.0;
            else
             T[K]=5.0;
           RB[K]=(PREC[K]*2.54-RS*S[K]-RT*T[K])/B[K];
           DURB[K]=RB[K]/R;
       }

      if (DURB[K] < 0.0) DURB[K]=0.0;    // !added

      PREC[K]=PREC[K]*2.54 * 10.0;  // converse back to cm, and then to mm

      TPREC=TPREC+PREC[K];
      RTN=RTN+B[K]+T[K];
      ATEMP=ATEMP+TEMP[K];
 }
      ATEMP=ATEMP/NM;

//   Following is to calculate the Julian days for NY years given a
//   specific date [MONTH and DAY]. Be aware that NJD=NY*365 while
//   NJ=[NY+1]*365. The output is in an arrary JDAY in dimension of NJD.
//   NI is the Julian day for the date DAY/MONTH.

//      data DX/0,1,-1,0,0,1,1,2,3,3,4,4/

      NI=(MONTH-1)*30+DAY+DX[MONTH];

//      data DD/1,-2,1,0,1,0,1,1,0,1,0,1/
      NS=0;
      NSY=0;

      for (I=1; I<=(NY+1); I++) {
          NSY=NSY+1;
          for (K=1; K<=12; K++) {
              DT=30+DD[K];
              for (JJ=1; JJ<=DT; JJ++) {
                  NS=NS+1;
//        The calculation of Julian day.
                  JJD[NS]=NS-(NSY-1)*365;
              }
          }
      }

      for (I = 1; I <= (NY*365); I++) JDAY[I]=JJD[I+NI-1];

//   The end of calculating Julian day.

      NN=0;

      for (I = 1; I <= NY; I++) {
          RAINN=0;

          for (K=((I-1)*12+1); K <= I*12; K++) {
              DT=30+DD[K-(I-1)*12];
              KTT=B[K]+T[K];
              KTD=DT/KTT;
              KDD=DT-KTT*KTD;
              BB=B[K];
              TT=T[K];

              for (JJ=1; JJ<=KTT; JJ++) {
                 if (BB > 0.0) {
                    BB=BB-1.0;
                    RAINN=RAINN+1;

                    for (L=1; L<=KTD; L++) {
                       NN=NN+1;
                       DRYD[NN]=KTD;
                       RAINNUMB[NN]=RAINN;
                       TEM[NN]=TEMP[K];

                       VAPOR[NN] = VAPO[K];

                       RAININTE[NN]=0.0;
                       RAINDUR[NN]=0.0;
                       if (L == KTD) {
                           RAININTE[NN]=5.0; // unit with mm /hr
                           RAINDUR[NN]=DURB[K];
                       }
                    }
                  }

                 if (TT > 0.0) {
                    TT=TT-1.0;
                    RAINN=RAINN+1;
                    if (JJ == 1)
                       KKTD=KTD+KDD;
                    else
                       KKTD=KTD;
                    for (L=1; L <= KKTD; L++) {
                        NN=NN+1;
                        DRYD[NN]=KKTD;
                        RAINNUMB[NN]=RAINN;
                        TEM[NN]=TEMP[K];

                        VAPOR[NN] = VAPO[K];

                        RAININTE[NN]=0.0;
                        RAINDUR[NN]=0.0;
                        if (L == KKTD) {
                           RAININTE[NN]=5.0; //unit mm/hr
                           RAINDUR[NN]=DURT[K];
                         }
                    }
                 }

               }  // end of for J

          }   // end of for K
      }   // end of for I

//  The following block is to interpolate the monthly mean tempearature into
//  daily temperature assuming on Day 15 the daily temperature is equal to
//  the monthly mean and Day 1 and Day 30 or 31 the temperatures are equal
//  to the average of the two consective months' means.

      for (K=1; K<=NM; K++) {
         if (K < 12)
           {
           TMP1[K+1]=(TEMP[K]+TEMP[K+1])*0.5;
           VAP1[K+1]=(VAPO[K]+VAPO[K+1])*0.5;
           }
         else {
           TMP1[1]=(TEMP[1]+TEMP[12])*0.5;
           TMP1[13]=(TEMP[1]+TEMP[12])*0.5;
           VAP1[1]=(VAPO[1]+VAPO[12])*0.5;
           VAP1[13]=(VAPO[1]+VAPO[12])*0.5;
         }
      }

      for (K=1; K<=NM; K++)
           {
            TMP15[K]=2.0*TEMP[K]-((DD[K]+15.0)/(DD[K]+30.0))*TMP1[K+1]-15.0/(DD[K]+30.0)*TMP1[K];
            VAP15[K]=2.0*VAPO[K]-((DD[K]+15.0)/(DD[K]+30.0))*VAP1[K+1]-15.0/(DD[K]+30.0)*VAP1[K];
           }

      NS=0;
      for (K=1; K<=NM; K++) {
         DT=DD[K]+30;
         for (JJ=1; JJ<=DT; JJ++) {
            NS=NS+1;
            if (JJ <= 15)
               {
                TEM[NS]=TMP1[K]+(TMP15[K]-TMP1[K])*JJ/15.0;
                VAPOR[NS]=VAP1[K]+(VAP15[K]-VAP1[K])*JJ/15.0;
               }
            else
               {
                TEM[NS]=TMP15[K]+(TMP1[K+1]-TMP15[K])*(JJ-15)/(DT-15.0);
                VAPOR[NS]=VAP15[K]+(VAP1[K+1]-VAP15[K])*(JJ-15)/(DT-15.0);
               }
         }
      }

/*
    fout.open("test.out", ios::out);
    fout << "K     TEMP[C]     PREC[cm]    RB    B    DIRB[hr]    RT    T    DURT[hr]    RS    S    DURS[hr]" << endl;
    for (K=1; K<=NM; K++) {
        fout << MMONTH[K] << "\t" << TEMP[K] <<"\t" << PREC[K] <<"\t" << RB[K] <<"\t";
        fout << "\t" << B[K] << "\t" << DURB[K] <<"\t" << RT <<"\t" << T[K] <<"\t";
        fout << DURT[K] << "\t" <<  RS << "\t" << S[K] << "\t" << DURS[K] << endl;
    }
    fout << endl << endl;
    fout << "Total Rainfall[cm]      Total Rainfall events     Average Temperature[C]" << endl;
    fout << TPREC << "\t" << "\t" << "\t" << "\t" << RTN << "\t" << "\t" << "\t" << "\t" << ATEMP << endl;

    fout << endl << endl;
    fout << "JDAY        DRYD        RAINNUMB       TEM        RAININTE        RAINDUR" <<endl;
 */
    for (K=1; K<=NJD; K++) {
//        fout << JDAY[K] << "\t" << DRYD[K] << "\t" << RAINNUMB[K] << "\t";
//        fout << TEM[K] << "    \t" << RAININTE[K] << "\t" << RAINDUR[K] << endl;
        precdaily[K] =  RAININTE[K] *  RAINDUR[K];
        tempdaily[K] = TEM[K];
        vapdaily[K] = VAPOR[K];
    }
//    fout.close();


 };
 
 
