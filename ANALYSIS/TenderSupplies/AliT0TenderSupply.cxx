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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  T0 Tender supply    //
//  //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <AliESDEvent.h>
#include <AliESDtrack.h>

#include <AliTender.h>
#include <AliT0TenderSupply.h>


ClassImp(AliT0TenderSupply)

//________________________________________________________________________
AliT0TenderSupply::AliT0TenderSupply():
  AliTenderSupply(),
  fCorrectMeanTime(kFALSE),
  fCorrectStartTimeOnAmplSatur(kFALSE),
  fAmplitudeThreshold(100) 
{
  //
  // default constructor
  //
  for(int i=0; i<3; i++) fTimeOffset[i]=0;
}

//________________________________________________________________________
AliT0TenderSupply::AliT0TenderSupply(const char *name, const AliTender *tender):
  AliTenderSupply(name,tender),
  fCorrectMeanTime(kFALSE),
  fCorrectStartTimeOnAmplSatur(kFALSE),
  fAmplitudeThreshold(100) 
{
  //
  // constructor
  //
  for(int i=0; i<3; i++) fTimeOffset[i]=0;
}

//________________________________________________________________________
AliT0TenderSupply::~AliT0TenderSupply(){
  //
  // destructor
  //
  
}

//________________________________________________________________________
void AliT0TenderSupply::Init(){
  //
  // Init
  //
  Int_t run = fTender->GetRun();
  if (run == 0) return;                // to skip first init, when we don't have yet a run number

  fCorrectMeanTime = kFALSE; //reset
  for(int i=0; i<3; i++) fTimeOffset[i]=0;

  if(run==167706){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=37.636; fTimeOffset[1]=38.2358;  fTimeOffset[2] =37.2348;}
  if(run==167711){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=7.71361; fTimeOffset[1]=15.0231;  fTimeOffset[2] =18.0363;}
  if(run==167713){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=19.4678; fTimeOffset[1]=8.26494;  fTimeOffset[2] =39.7651;}
  if(run==167806){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=32.8404; fTimeOffset[1]=25.6742;  fTimeOffset[2] =43.0574;}
  if(run==167807){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=39.2029; fTimeOffset[1]=15.9077;  fTimeOffset[2] =46.8384;}
  if(run==167808){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=27.2107; fTimeOffset[1]=20.1576;  fTimeOffset[2] =27.9818;}
  if(run==167813){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=22.5833; fTimeOffset[1]=16.1395;  fTimeOffset[2] =31.3737;}
  if(run==167814){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=37.756; fTimeOffset[1]=20.6428;  fTimeOffset[2] =56.2365;}
  if(run==167818){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=38.7765; fTimeOffset[1]=37.3998;  fTimeOffset[2] =35.8552;}
  if(run==167902){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=41.0647; fTimeOffset[1]=31.3968;  fTimeOffset[2] =33.8023;}
  if(run==167903){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=35.2311; fTimeOffset[1]=49.3209;  fTimeOffset[2] =30.5427;}
  if(run==167915){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=33.2731; fTimeOffset[1]=27.4272;  fTimeOffset[2] =35.8073;}
  if(run==167920){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=31.8057; fTimeOffset[1]=21.8993;  fTimeOffset[2] =36.615;}
  if(run==167921){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=30.6407; fTimeOffset[1]=8.60546;  fTimeOffset[2] =53.1829;}
  if(run==167985){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=31.7908; fTimeOffset[1]=25.753;  fTimeOffset[2] =37.3178;}
  if(run==167986){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=22.331; fTimeOffset[1]=4.7895;  fTimeOffset[2] =40.9276;}
  if(run==167987){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=30.4746; fTimeOffset[1]=15.5886;  fTimeOffset[2] =50.3252;}
  if(run==167988){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=45.1788; fTimeOffset[1]=29.2541;  fTimeOffset[2] =67.222;}
  if(run==168066){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=23.5272; fTimeOffset[1]=26.8631;  fTimeOffset[2] =31.2057;}
  if(run==168068){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=32.4972; fTimeOffset[1]=9.76024;  fTimeOffset[2] =48.3329;}
  if(run==168069){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=25.6839; fTimeOffset[1]=11.4422;  fTimeOffset[2] =34.8375;}
  if(run==168076){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=39.2055; fTimeOffset[1]=25.607;  fTimeOffset[2] =52.0409;}
  if(run==168103){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=67.7271; fTimeOffset[1]=50.5392;  fTimeOffset[2] =48.0866;}
  if(run==168104){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=42.6881; fTimeOffset[1]=43.5591;  fTimeOffset[2] =48.7308;}
  if(run==168105){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=29.958; fTimeOffset[1]=8.65483;  fTimeOffset[2] =42.7843;}
  if(run==168107){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=39.0088; fTimeOffset[1]=26.1971;  fTimeOffset[2] =53.5561;}
  if(run==168108){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=47.8417; fTimeOffset[1]=33.5841;  fTimeOffset[2] =60.2913;}
  if(run==168115){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=26.9248; fTimeOffset[1]=16.2735;  fTimeOffset[2] =37.1457;}
  if(run==168171){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=841.144; fTimeOffset[1]=927.182;  fTimeOffset[2] =751.284;}
  if(run==168172){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=868.888; fTimeOffset[1]=972.012;  fTimeOffset[2] =766.496;}
  if(run==168173){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=833.223; fTimeOffset[1]=925.62;  fTimeOffset[2] =732.57;}
  if(run==168175){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=807.481; fTimeOffset[1]=900.059;  fTimeOffset[2] =718.232;}
  if(run==168181){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=78.7363; fTimeOffset[1]=115.034;  fTimeOffset[2] =39.0348;}
  if(run==168203){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=34.9107; fTimeOffset[1]=18.5032;  fTimeOffset[2] =50.9434;}
  if(run==168204){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=800.427; fTimeOffset[1]=875.242;  fTimeOffset[2] =744.428;}
  if(run==168205){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=70.2692; fTimeOffset[1]=57.6408;  fTimeOffset[2] =75.372;}
  if(run==168206){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=754.293; fTimeOffset[1]=846.635;  fTimeOffset[2] =670.418;}
  if(run==168207){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=26.3374; fTimeOffset[1]=-54.9438;  fTimeOffset[2] =76.9511;}
  if(run==168208){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=73.6404; fTimeOffset[1]=20.3787;  fTimeOffset[2] =75.0944;}
  if(run==168212){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=57.4094; fTimeOffset[1]=38.2747;  fTimeOffset[2] =41.9405;}
  if(run==168213){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=35.4503; fTimeOffset[1]=27.5827;  fTimeOffset[2] =39.1856;}
  if(run==168310){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=28.6722; fTimeOffset[1]=33.2269;  fTimeOffset[2] =27.9654;}
  if(run==168311){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=40.203; fTimeOffset[1]=19.1132;  fTimeOffset[2] =56.441;}
  if(run==168318){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=20.5591; fTimeOffset[1]=26.1756;  fTimeOffset[2] =16.0513;}
  if(run==168322){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=56.2993; fTimeOffset[1]=44.9021;  fTimeOffset[2] =61.1342;}
  if(run==168325){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=44.1044; fTimeOffset[1]=24.2727;  fTimeOffset[2] =61.4175;}
  if(run==168341){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=43.254; fTimeOffset[1]=40.1345;  fTimeOffset[2] =44.7359;}
  if(run==168342){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=55.7112; fTimeOffset[1]=48.5852;  fTimeOffset[2] =56.7952;}
  if(run==168356){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=-34.5093; fTimeOffset[1]=-5.55039;  fTimeOffset[2] =-65.6368;}
  if(run==168361){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=-248.811; fTimeOffset[1]=-250.006;  fTimeOffset[2] =-153.319;}
  if(run==168362){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=53.8183; fTimeOffset[1]=48.2561;  fTimeOffset[2] =54.1145;}
  if(run==168458){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=33.4396; fTimeOffset[1]=32.8394;  fTimeOffset[2] =39.4757;}
  if(run==168460){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=43.74; fTimeOffset[1]=49.8685;  fTimeOffset[2] =41.4665;}
  if(run==168461){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=29.7889; fTimeOffset[1]=27.8432;  fTimeOffset[2] =38.4801;}
  if(run==168464){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=33.1291; fTimeOffset[1]=23.5964;  fTimeOffset[2] =38.9655;}
  if(run==168467){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=37.1316; fTimeOffset[1]=48.4763;  fTimeOffset[2] =19.1521;}
  if(run==168511){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=45.4487; fTimeOffset[1]=31.5609;  fTimeOffset[2] =44.3649;}
  if(run==168512){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=34.615; fTimeOffset[1]=24.6349;  fTimeOffset[2] =38.7753;}
  if(run==168514){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=28.8178; fTimeOffset[1]=22.1466;  fTimeOffset[2] =40.6307;}
  if(run==168777){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=32.2621; fTimeOffset[1]=17.7215;  fTimeOffset[2] =44.2911;}
  if(run==168826){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=27.7899; fTimeOffset[1]=19.4315;  fTimeOffset[2] =34.24;}
  if(run==168984){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=21.9251; fTimeOffset[1]=4.93046;  fTimeOffset[2] =27.7368;}
  if(run==168988){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=32.3264; fTimeOffset[1]=19.2958;  fTimeOffset[2] =42.3047;}
  if(run==168992){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=40.2228; fTimeOffset[1]=27.3191;  fTimeOffset[2] =51.0479;}
  if(run==169035){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=16.2383; fTimeOffset[1]=12.7409;  fTimeOffset[2] =21.7943;}
  if(run==169044){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=18.3657; fTimeOffset[1]=2.523;  fTimeOffset[2] =39.3647;}
  if(run==169040){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=30     ; fTimeOffset[1]=-42;     fTimeOffset[2] = 0;}
  if(run==169045){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=29.9529; fTimeOffset[1]=21.3877;  fTimeOffset[2] =35.4287;}
  if(run==169094){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=16.3535; fTimeOffset[1]=12.5777;  fTimeOffset[2] =15.5841;}
  if(run==169099){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=24.34; fTimeOffset[1]=17.5334;  fTimeOffset[2] =28.812;}
  if(run==169143){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=34.3485; fTimeOffset[1]=-2.89818;  fTimeOffset[2] =67.149;}
  if(run==169145){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=30.1335; fTimeOffset[1]=13.647;  fTimeOffset[2] =44.4143;}
  if(run==169148){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=28.0664; fTimeOffset[1]=15.8024;  fTimeOffset[2] =32.7707;}
  if(run==169156){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=20.3655; fTimeOffset[1]=15.1023;  fTimeOffset[2] =28.4647;}
  if(run==169160){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=9.7256; fTimeOffset[1]=3.94123;  fTimeOffset[2] =16.3137;}
  if(run==169167){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=30.4764; fTimeOffset[1]=11.632;  fTimeOffset[2] =42.8243;}
  if(run==169238){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=15.5803; fTimeOffset[1]=17.6367;  fTimeOffset[2] =12.3469;}
  if(run==169411){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=40.604; fTimeOffset[1]=24.3761;  fTimeOffset[2] =56.2715;}
  if(run==169415){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=45.9346; fTimeOffset[1]=27.2209;  fTimeOffset[2] =54.7067;}
  if(run==169417){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=33.8899; fTimeOffset[1]=25.979;  fTimeOffset[2] =45.0143;}
  if(run==169418){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=34.6481; fTimeOffset[1]=31.0951;  fTimeOffset[2] =39.9321;}
  if(run==169419){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=42.138; fTimeOffset[1]=32.447;  fTimeOffset[2] =47.7478;}
  if(run==169420){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=51.7414; fTimeOffset[1]=42.8009;  fTimeOffset[2] =56.2767;}
  if(run==169498){ fCorrectMeanTime = kTRUE; fTimeOffset[0]=47.3734; fTimeOffset[1]=48.0008;  fTimeOffset[2] =45.101;}

  //-----
  /*

     fCorrectStartTimeOnAmplSatur = kFALSE;
     fAmplitudeThreshold = 100; //in mips

     if(167693<= run && run<=170593){  // LHC11h
     fCorrectStartTimeOnAmplSatur = kTRUE;
     fAmplitudeThreshold = 40; //in mips
     }
     */
}

//________________________________________________________________________
void AliT0TenderSupply::ProcessEvent(){
    //
    // loop over all online T0 candidates and flag
    // selected daughter tracks using the status bis of the TObject
    //

    AliESDEvent *event=fTender->GetEvent();
    if (!event) return;

    //Do something when the run number changed, like loading OCDB entries etc.
    if(fTender->RunChanged()){
        Init();
    }


    if(fCorrectStartTimeOnAmplSatur){
        //correct A side ORA on amplitude saturation
        const Double32_t* time = event->GetT0time();
        const Double32_t* amplitude = event->GetT0amplitude();

        Int_t idxOfFirstPmtA = -1;
        Double32_t timeOrA   = 99999;
        for(int ipmt=12; ipmt<24; ipmt++){ //loop over A side
            if( amplitude[ipmt] < fAmplitudeThreshold){
                if( time[ipmt] > -200 && time[ipmt]!=0 && time[ipmt] < timeOrA ){ 
                    timeOrA        = time[ipmt];
                    idxOfFirstPmtA = ipmt;
                }
            }
        }

        if(idxOfFirstPmtA>-1){ //a hit in aside with less than 40 mips
            const Double32_t* mean = event->GetT0TOF();
            Double32_t timeOrC = mean[2];
            Double32_t timeOrAplusOrC = (timeOrA+timeOrC)/2;

            event->SetT0TOF(0, timeOrAplusOrC);
            event->SetT0TOF(1, timeOrA);
        }
    }

    //...........................................
    if(fCorrectMeanTime){
        // correct mean time offsets  
        const Double32_t* mean = event->GetT0TOF();
        for(int it0=0; it0<3; it0++){
            if(-200 < mean[it0]){
                event->SetT0TOF(it0, mean[it0] - fTimeOffset[it0]); 
            }
        }
    }

}
