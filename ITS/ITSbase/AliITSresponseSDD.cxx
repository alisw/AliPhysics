/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//////////////////////////////////////////////////////
//  Base response class forITS                      //
//  It is used to set static data members           //
//  connected to parameters equal for all           //
//  the modules                                     //
//                                                  //
//                                                  //
//////////////////////////////////////////////////////

#include <TMath.h>

#include "AliITSresponseSDD.h"
#include <AliITSgeomTGeo.h>

const Float_t AliITSresponseSDD::fgkTimeOffsetDefault = 54.30;
const Float_t AliITSresponseSDD::fgkADC2keVDefault = 3.34;
const Float_t AliITSresponseSDD::fgkChargevsTimeDefault = 0.00355;
const Float_t AliITSresponseSDD::fgkADCvsDrTimeDefault = 0.0101;
const Float_t AliITSresponseSDD::fgkCarlosRXClockPeriod = 25.;
ClassImp(AliITSresponseSDD)

//_________________________________________________________________________
AliITSresponseSDD::AliITSresponseSDD():
TObject(),
  fTimeOffset(fgkTimeOffsetDefault),
  fADC2keV(fgkADC2keVDefault),
  fChargevsTime(fgkChargevsTimeDefault)
{
  // default constructor
  for(Int_t i=0; i<kNSDDmods;i++){
    fTimeZero[i]=fgkTimeOffsetDefault;
    fDeltaVDrift[i] = fDeltaVDrift[i+kNSDDmods] = 0.;
    fADCtokeV[i]=fgkADC2keVDefault;
    fADCvsDriftTime[i]=fgkADCvsDrTimeDefault;
    fADCvsDriftTimeMC[i]=fgkADCvsDrTimeDefault;
  }  
  SetMCDefaults();
  SetVDCorr2Side(kTRUE); // default for new objects will be separate corrections for 2 sides (bwd compatible)
  //  SetVDCorrMult(kTRUE); // default for new objects will have multiplicative correction v'=(1+corr)*v (bwd compatible)
}
//_________________________________________________________________________
void AliITSresponseSDD::SetMCDefaults(){
  // SetMC defaults -> calibrated on LHC14jd4 (R. Russo, 16-Jan-15)
  fADCvsDriftTimeMC[1]=0.008384;
  fADCvsDriftTimeMC[2]=0.008747;
  fADCvsDriftTimeMC[3]=0.007127;
  fADCvsDriftTimeMC[4]=0.008040;
  fADCvsDriftTimeMC[5]=0.007344;
  fADCvsDriftTimeMC[6]=0.010241;
  fADCvsDriftTimeMC[7]=0.007235;
  fADCvsDriftTimeMC[8]=0.006573;
  fADCvsDriftTimeMC[9]=0.005941;
  fADCvsDriftTimeMC[10]=0.007528;
  fADCvsDriftTimeMC[11]=0.008039;
  fADCvsDriftTimeMC[12]=0.007193;
  fADCvsDriftTimeMC[13]=0.008126;
  fADCvsDriftTimeMC[14]=0.010100;
  fADCvsDriftTimeMC[15]=0.008038;
  fADCvsDriftTimeMC[16]=0.008456;
  fADCvsDriftTimeMC[17]=0.007424;
  fADCvsDriftTimeMC[18]=0.010100;
  fADCvsDriftTimeMC[19]=0.010100;
  fADCvsDriftTimeMC[20]=0.010100;
  fADCvsDriftTimeMC[21]=0.007915;
  fADCvsDriftTimeMC[22]=0.007949;
  fADCvsDriftTimeMC[23]=0.006127;
  fADCvsDriftTimeMC[24]=0.007887;
  fADCvsDriftTimeMC[25]=0.008140;
  fADCvsDriftTimeMC[26]=0.006520;
  fADCvsDriftTimeMC[27]=0.007024;
  fADCvsDriftTimeMC[28]=0.007271;
  fADCvsDriftTimeMC[29]=0.007348;
  fADCvsDriftTimeMC[30]=0.007820;
  fADCvsDriftTimeMC[31]=0.007544;
  fADCvsDriftTimeMC[32]=0.004143;
  fADCvsDriftTimeMC[33]=0.007644;
  fADCvsDriftTimeMC[34]=0.007170;
  fADCvsDriftTimeMC[35]=0.006653;
  fADCvsDriftTimeMC[36]=0.008632;
  fADCvsDriftTimeMC[37]=0.007728;
  fADCvsDriftTimeMC[38]=0.008096;
  fADCvsDriftTimeMC[39]=0.010100;
  fADCvsDriftTimeMC[40]=0.010100;
  fADCvsDriftTimeMC[41]=0.010100;
  fADCvsDriftTimeMC[42]=0.007480;
  fADCvsDriftTimeMC[43]=0.007949;
  fADCvsDriftTimeMC[44]=0.007662;
  fADCvsDriftTimeMC[45]=0.007848;
  fADCvsDriftTimeMC[46]=0.006875;
  fADCvsDriftTimeMC[47]=0.007592;
  fADCvsDriftTimeMC[48]=0.007716;
  fADCvsDriftTimeMC[49]=0.007461;
  fADCvsDriftTimeMC[50]=0.007555;
  fADCvsDriftTimeMC[51]=0.007551;
  fADCvsDriftTimeMC[52]=0.008068;
  fADCvsDriftTimeMC[53]=0.008595;
  fADCvsDriftTimeMC[54]=0.006893;
  fADCvsDriftTimeMC[55]=0.007608;
  fADCvsDriftTimeMC[56]=0.007715;
  fADCvsDriftTimeMC[57]=0.007554;
  fADCvsDriftTimeMC[58]=0.007716;
  fADCvsDriftTimeMC[59]=0.006613;
  fADCvsDriftTimeMC[60]=0.008064;
  fADCvsDriftTimeMC[61]=0.006857;
  fADCvsDriftTimeMC[62]=0.007170;
  fADCvsDriftTimeMC[63]=0.006693;
  fADCvsDriftTimeMC[64]=0.006828;
  fADCvsDriftTimeMC[65]=0.006384;
  fADCvsDriftTimeMC[66]=0.006835;
  fADCvsDriftTimeMC[67]=0.007131;
  fADCvsDriftTimeMC[68]=0.008520;
  fADCvsDriftTimeMC[69]=0.007338;
  fADCvsDriftTimeMC[70]=0.008039;
  fADCvsDriftTimeMC[71]=0.005399;
  fADCvsDriftTimeMC[72]=0.007444;
  fADCvsDriftTimeMC[73]=0.007333;
  fADCvsDriftTimeMC[74]=0.007522;
  fADCvsDriftTimeMC[75]=0.007645;
  fADCvsDriftTimeMC[76]=0.008382;
  fADCvsDriftTimeMC[77]=0.010100;
  fADCvsDriftTimeMC[78]=0.005804;
  fADCvsDriftTimeMC[79]=0.007902;
  fADCvsDriftTimeMC[80]=0.007034;
  fADCvsDriftTimeMC[81]=0.008094;
  fADCvsDriftTimeMC[82]=0.007014;
  fADCvsDriftTimeMC[83]=0.007301;
  fADCvsDriftTimeMC[84]=0.010100;
  fADCvsDriftTimeMC[85]=0.007371;
  fADCvsDriftTimeMC[86]=0.005768;
  fADCvsDriftTimeMC[87]=0.007079;
  fADCvsDriftTimeMC[88]=0.009297;
  fADCvsDriftTimeMC[89]=0.008228;
  fADCvsDriftTimeMC[90]=0.007923;
  fADCvsDriftTimeMC[91]=0.005789;
  fADCvsDriftTimeMC[92]=0.005612;
  fADCvsDriftTimeMC[93]=0.007189;
  fADCvsDriftTimeMC[94]=0.007237;
  fADCvsDriftTimeMC[95]=0.007478;
  fADCvsDriftTimeMC[96]=0.006879;
  fADCvsDriftTimeMC[97]=0.007035;
  fADCvsDriftTimeMC[98]=0.006412;
  fADCvsDriftTimeMC[99]=0.006723;
  fADCvsDriftTimeMC[100]=0.006292;
  fADCvsDriftTimeMC[101]=0.007809;
  fADCvsDriftTimeMC[102]=0.008281;
  fADCvsDriftTimeMC[103]=0.010100;
  fADCvsDriftTimeMC[104]=0.007527;
  fADCvsDriftTimeMC[105]=0.010100;
  fADCvsDriftTimeMC[106]=0.007810;
  fADCvsDriftTimeMC[107]=0.007578;
  fADCvsDriftTimeMC[108]=0.010100;
  fADCvsDriftTimeMC[109]=0.009124;
  fADCvsDriftTimeMC[110]=0.006268;
  fADCvsDriftTimeMC[111]=0.005359;
  fADCvsDriftTimeMC[112]=0.007384;
  fADCvsDriftTimeMC[113]=0.007310;
  fADCvsDriftTimeMC[114]=0.008525;
  fADCvsDriftTimeMC[115]=0.010100;
  fADCvsDriftTimeMC[116]=0.006898;
  fADCvsDriftTimeMC[117]=0.007718;
  fADCvsDriftTimeMC[118]=0.007631;
  fADCvsDriftTimeMC[119]=0.006599;
  fADCvsDriftTimeMC[120]=0.006054;
  fADCvsDriftTimeMC[121]=0.007661;
  fADCvsDriftTimeMC[122]=0.008095;
  fADCvsDriftTimeMC[123]=0.007686;
  fADCvsDriftTimeMC[124]=0.007319;
  fADCvsDriftTimeMC[125]=0.008086;
  fADCvsDriftTimeMC[126]=0.004947;
  fADCvsDriftTimeMC[127]=0.008500;
  fADCvsDriftTimeMC[128]=0.007279;
  fADCvsDriftTimeMC[129]=0.007272;
  fADCvsDriftTimeMC[130]=0.007216;
  fADCvsDriftTimeMC[131]=0.005968;
  fADCvsDriftTimeMC[132]=0.005603;
  fADCvsDriftTimeMC[133]=0.007193;
  fADCvsDriftTimeMC[134]=0.007272;
  fADCvsDriftTimeMC[135]=0.007144;
  fADCvsDriftTimeMC[136]=0.005582;
  fADCvsDriftTimeMC[137]=0.007927;
  fADCvsDriftTimeMC[138]=0.007598;
  fADCvsDriftTimeMC[139]=0.006109;
  fADCvsDriftTimeMC[140]=0.005826;
  fADCvsDriftTimeMC[141]=0.007517;
  fADCvsDriftTimeMC[142]=0.007181;
  fADCvsDriftTimeMC[143]=0.006917;
  fADCvsDriftTimeMC[144]=0.006550;
  fADCvsDriftTimeMC[145]=0.008342;
  fADCvsDriftTimeMC[146]=0.007485;
  fADCvsDriftTimeMC[147]=0.008732;
  fADCvsDriftTimeMC[148]=0.007370;
  fADCvsDriftTimeMC[149]=0.006544;
  fADCvsDriftTimeMC[150]=0.007054;
  fADCvsDriftTimeMC[151]=0.010100;
  fADCvsDriftTimeMC[152]=0.007192;
  fADCvsDriftTimeMC[153]=0.007222;
  fADCvsDriftTimeMC[154]=0.006940;
  fADCvsDriftTimeMC[155]=0.007468;
  fADCvsDriftTimeMC[156]=0.007776;
  fADCvsDriftTimeMC[157]=0.007824;
  fADCvsDriftTimeMC[158]=0.008047;
  fADCvsDriftTimeMC[159]=0.006877;
  fADCvsDriftTimeMC[160]=0.006306;
  fADCvsDriftTimeMC[161]=0.006674;
  fADCvsDriftTimeMC[162]=0.007904;
  fADCvsDriftTimeMC[163]=0.005996;
  fADCvsDriftTimeMC[164]=0.009781;
  fADCvsDriftTimeMC[165]=0.006822;
  fADCvsDriftTimeMC[166]=0.010100;
  fADCvsDriftTimeMC[167]=0.006567;
  fADCvsDriftTimeMC[168]=0.007594;
  fADCvsDriftTimeMC[169]=0.006436;
  fADCvsDriftTimeMC[170]=0.010100;
  fADCvsDriftTimeMC[171]=0.006541;
  fADCvsDriftTimeMC[172]=0.010100;
  fADCvsDriftTimeMC[173]=0.005869;
  fADCvsDriftTimeMC[174]=0.006169;
  fADCvsDriftTimeMC[175]=0.007122;
  fADCvsDriftTimeMC[176]=0.005458;
  fADCvsDriftTimeMC[177]=0.006642;
  fADCvsDriftTimeMC[178]=0.008467;
  fADCvsDriftTimeMC[179]=0.007741;
  fADCvsDriftTimeMC[180]=0.008877;
  fADCvsDriftTimeMC[181]=0.007502;
  fADCvsDriftTimeMC[182]=0.007330;
  fADCvsDriftTimeMC[183]=0.006217;
  fADCvsDriftTimeMC[184]=0.007978;
  fADCvsDriftTimeMC[185]=0.007441;
  fADCvsDriftTimeMC[186]=0.007820;
  fADCvsDriftTimeMC[187]=0.007762;
  fADCvsDriftTimeMC[188]=0.007411;
  fADCvsDriftTimeMC[189]=0.007028;
  fADCvsDriftTimeMC[190]=0.008228;
  fADCvsDriftTimeMC[191]=0.006653;
  fADCvsDriftTimeMC[192]=0.007009;
  fADCvsDriftTimeMC[193]=0.005980;
  fADCvsDriftTimeMC[194]=0.006573;
  fADCvsDriftTimeMC[195]=0.010100;
  fADCvsDriftTimeMC[196]=0.007062;
  fADCvsDriftTimeMC[197]=0.007691;
  fADCvsDriftTimeMC[198]=0.008381;
  fADCvsDriftTimeMC[199]=0.007031;
  fADCvsDriftTimeMC[200]=0.007874;
  fADCvsDriftTimeMC[201]=0.011088;
  fADCvsDriftTimeMC[202]=0.007996;
  fADCvsDriftTimeMC[203]=0.007077;
  fADCvsDriftTimeMC[204]=0.008236;
  fADCvsDriftTimeMC[205]=0.008347;
  fADCvsDriftTimeMC[206]=0.008797;
  fADCvsDriftTimeMC[207]=0.005973;
  fADCvsDriftTimeMC[208]=0.007560;
  fADCvsDriftTimeMC[209]=0.007505;
  fADCvsDriftTimeMC[210]=0.007965;
  fADCvsDriftTimeMC[211]=0.007379;
  fADCvsDriftTimeMC[212]=0.007974;
  fADCvsDriftTimeMC[213]=0.006762;
  fADCvsDriftTimeMC[214]=0.010100;
  fADCvsDriftTimeMC[215]=0.007151;
  fADCvsDriftTimeMC[216]=0.005802;
  fADCvsDriftTimeMC[217]=0.006873;
  fADCvsDriftTimeMC[218]=0.007286;
  fADCvsDriftTimeMC[219]=0.006774;
  fADCvsDriftTimeMC[220]=0.007708;
  fADCvsDriftTimeMC[221]=0.007048;
  fADCvsDriftTimeMC[222]=0.004330;
  fADCvsDriftTimeMC[223]=0.010100;
  fADCvsDriftTimeMC[224]=0.007970;
  fADCvsDriftTimeMC[225]=0.010100;
  fADCvsDriftTimeMC[226]=0.008034;
  fADCvsDriftTimeMC[227]=0.006442;
  fADCvsDriftTimeMC[228]=0.007122;
  fADCvsDriftTimeMC[229]=0.007737;
  fADCvsDriftTimeMC[230]=0.007755;
  fADCvsDriftTimeMC[231]=0.007995;
  fADCvsDriftTimeMC[232]=0.006997;
  fADCvsDriftTimeMC[233]=0.008777;
  fADCvsDriftTimeMC[234]=0.010994;
  fADCvsDriftTimeMC[235]=0.005986;
  fADCvsDriftTimeMC[236]=0.007349;
  fADCvsDriftTimeMC[237]=0.008129;
  fADCvsDriftTimeMC[238]=0.005306;
  fADCvsDriftTimeMC[239]=0.008648;
  fADCvsDriftTimeMC[240]=0.007879;
  fADCvsDriftTimeMC[241]=0.007558;
  fADCvsDriftTimeMC[242]=0.007368;
  fADCvsDriftTimeMC[243]=0.010837;
  fADCvsDriftTimeMC[244]=0.006780;
  fADCvsDriftTimeMC[245]=0.005178;
  fADCvsDriftTimeMC[246]=0.010100;
  fADCvsDriftTimeMC[247]=0.007705;
  fADCvsDriftTimeMC[248]=0.008350;
  fADCvsDriftTimeMC[249]=0.008654;
  fADCvsDriftTimeMC[250]=0.007147;
  fADCvsDriftTimeMC[251]=0.008399;
  fADCvsDriftTimeMC[252]=0.010100;
  fADCvsDriftTimeMC[253]=0.007036;
  fADCvsDriftTimeMC[254]=0.007749;
  fADCvsDriftTimeMC[255]=0.009103;
  fADCvsDriftTimeMC[256]=0.007138;
  fADCvsDriftTimeMC[257]=0.006691;
  fADCvsDriftTimeMC[258]=0.006366;
  fADCvsDriftTimeMC[259]=0.006701;
}
//_________________________________________________________________________
void AliITSresponseSDD::SetHalfLadderATimeZero(Int_t lay, Int_t lad, Float_t tzero){
  // Sets time Zero for all modules of a ladder on side A (Z>0)
  Int_t minMod,maxMod;
  if(lay==3){
    minMod=1; 
    maxMod=3;
    if(lad>kNLaddersLay3){
      AliError(Form("Ladder number %d out of range",lad));
      return;
    }
  }else if(lay==4){
    minMod=1; 
    maxMod=4;
    if(lad>kNLaddersLay4){
      AliError(Form("Ladder number %d out of range",lad));
      return;
    }
  }else{
    AliError(Form("Layer number %d out of range",lay));
    return;
  }
  for(Int_t iMod=minMod; iMod<=maxMod; iMod++){
    Int_t modIndex=AliITSgeomTGeo::GetModuleIndex(lay,lad,iMod);
    SetModuleTimeZero(modIndex,tzero);
  }
}
//_________________________________________________________________________
void AliITSresponseSDD::SetHalfLadderCTimeZero(Int_t lay, Int_t lad, Float_t tzero){
  // Sets time Zero for all modules of a ladder on side C (Z<0)
  Int_t minMod,maxMod;
  if(lay==3){
    minMod=4; 
    maxMod=6;
    if(lad>kNLaddersLay3){
      AliError(Form("Ladder number %d out of range",lad));
      return;
    }
  }else if(lay==4){
    minMod=5; 
    maxMod=8;
    if(lad>kNLaddersLay4){
      AliError(Form("Ladder number %d out of range",lad));
      return;
    }
  }else{
    AliError(Form("Layer number %d out of range",lay));
    return;
  }
  for(Int_t iMod=minMod; iMod<=maxMod; iMod++){
    Int_t modIndex=AliITSgeomTGeo::GetModuleIndex(lay,lad,iMod);
    SetModuleTimeZero(modIndex,tzero);
  }
}
//_________________________________________________________________________
void AliITSresponseSDD::PrintChargeCalibrationParams() const{
  // Dump charge calibration parameters

  printf("ADC vs. drift time corr=%f\n",GetChargevsTime());
  printf("-------------------------------------\n");
  printf("Layer 3\n");
  for(Int_t ilad=1; ilad<=14; ilad++){
    for(Int_t idet=1; idet<=6;idet++){
      Int_t modIndex=AliITSgeomTGeo::GetModuleIndex(3,ilad,idet);
      Float_t tz=GetADCtokeV(modIndex);
      printf("%7.2f   ",tz);
    }
    printf("\n");
  }
  printf("\n");
  printf("Layer 4\n");
  for(Int_t ilad=1; ilad<=22; ilad++){
    for(Int_t idet=1; idet<=8;idet++){
      Int_t modIndex=AliITSgeomTGeo::GetModuleIndex(4,ilad,idet);
      Float_t tz=GetADCtokeV(modIndex);
      printf("%7.2f   ",tz);
    }
    printf("\n");
  }  
}
//_________________________________________________________________________
void AliITSresponseSDD::PrintTimeZeroes() const{
  // Dump time zero values

  printf("Layer 3\n");
  for(Int_t ilad=1; ilad<=14; ilad++){
    for(Int_t idet=1; idet<=6;idet++){
      Int_t modIndex=AliITSgeomTGeo::GetModuleIndex(3,ilad,idet);
      Float_t tz=GetTimeZero(modIndex);
      printf("%7.2f   ",tz);
    }
    printf("\n");
  }
  printf("\n");
  printf("Layer 4\n");
  for(Int_t ilad=1; ilad<=22; ilad++){
    for(Int_t idet=1; idet<=8;idet++){
      Int_t modIndex=AliITSgeomTGeo::GetModuleIndex(4,ilad,idet);
      Float_t tz=GetTimeZero(modIndex);
      printf("%7.2f   ",tz);
    }
    printf("\n");
  }
  
}
//_________________________________________________________________________
void AliITSresponseSDD::PrintVdriftCorerctions() const{
  // Dump corrections to vdrift

  for(Int_t iMod=240; iMod<500; iMod++){
    printf("Module %d   dVleft=%f   dVright=%f\n",iMod,GetDeltaVDrift(iMod,0),GetDeltaVDrift(iMod,1));
  }
}
