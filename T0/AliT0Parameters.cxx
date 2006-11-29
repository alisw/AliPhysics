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

/* $Id:  */

//____________________________________________________________________
//                                                                          
// T0 - T0. 
//
// This class is a singleton that handles various parameters of
// the T0 detectors.  
// Eventually, this class will use the Conditions DB to get the
// various parameters, which code can then request from here.
//                                                       
#include "AliLog.h"		  
#include "AliT0Parameters.h"	  
#include "AliT0CalibData.h"   
#include <AliCDBManager.h>        
#include <AliCDBEntry.h>          
#include <AliCDBStorage.h>         
#include <Riostream.h>

AliT0CalibData* AliT0Parameters::fgCalibData = 0;
//====================================================================
ClassImp(AliT0Parameters)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliT0Parameters* AliT0Parameters::fgInstance = 0;
//____________________________________________________________________
AliT0Parameters* 
AliT0Parameters::Instance() 
{
  // Get static instance 
  if (!fgInstance) fgInstance = new AliT0Parameters;
  return fgInstance;
}

//____________________________________________________________________
AliT0Parameters::AliT0Parameters()
  :fIsInit(kFALSE),fPh2Mip(0),fmV2Mip(0),fChannelWidth(0),fmV2Channel(0),fQTmin(0),fQTmax(0),fFixedGain(0),fSlewingLED(),fSlewingRec(),fPMTeff(),fTimeDelayLED(0),fTimeDelayCFD(0),fTimeDelayTVD(0),fCalibentry() 
{
  // Default constructor 

  for (Int_t ipmt=0; ipmt<24; ipmt++)
    {
      SetTimeDelayCablesCFD(ipmt);
      SetTimeDelayCablesLED(ipmt);
      SetTimeDelayElectronicCFD(ipmt);
      SetTimeDelayElectronicLED(ipmt);
      SetTimeDelayPMT(ipmt);
       SetVariableDelayLine(ipmt);
      SetSlewingLED(ipmt);
      SetSlewingRec(ipmt);
      SetPh2Mip();      
      SetmV2Mip();      
      SetChannelWidth();
      SetmV2channel();
      SetGain();
      SetQTmin();
      SetQTmax();
      SetPMTeff(ipmt);
   }
  SetTimeDelayTVD();
  SetZposition();
  
}

//__________________________________________________________________
void
AliT0Parameters::Init()
{
  // Initialize the parameters manager.  We need to get stuff from the
  // CDB here. 
  //   if (fIsInit) return;
  
  AliCDBManager* cdb      = AliCDBManager::Instance();
  //  AliCDBStorage *stor = cdb->GetStorage("local://$ALICE_ROOT");
  fCalibentry  = cdb->Get("T0/Calib/Gain_TimeDelay_Slewing_Walk");
  if (fCalibentry){
   fgCalibData  = (AliT0CalibData*)fCalibentry->GetObject();
  }

  fIsInit = kTRUE;
}


//__________________________________________________________________
Float_t
AliT0Parameters::GetGain(Int_t ipmt) const
{
  // Returns the calibrated gain for each PMT 
  // 

  if (!fCalibentry) 
    return fFixedGain;
   
  return fgCalibData->GetGain(ipmt);
}

//__________________________________________________________________
Float_t
AliT0Parameters::GetTimeDelayLED(Int_t ipmt) 
{
  // return time delay for LED channel
  // 
  if (!fCalibentry) {
    fTimeDelayLED = fTimeDelayCablesLED[ipmt] + fTimeDelayElectronicLED[ipmt] + fTimeDelayPMT[ipmt];
    return  fTimeDelayLED;
  } 
  return fgCalibData ->GetTimeDelayLED(ipmt);
}
//__________________________________________________________________
Float_t
AliT0Parameters::GetTimeDelayCFD(Int_t ipmt) 
{
  // return time delay for CFD channel
   // 
  if (!fCalibentry) 
    {
      fTimeDelayCFD = fTimeDelayCablesCFD[ipmt] + fTimeDelayElectronicCFD[ipmt] + fTimeDelayPMT[ipmt] + fVariableDelayLine[ipmt];
      return fTimeDelayCFD+37;
    }
   
  return fgCalibData->GetTimeDelayCFD(ipmt);
}

//__________________________________________________________________

void 
AliT0Parameters::SetSlewingLED(Int_t ipmt)
{
  //  Set Slweing Correction for LED channel 
     Float_t mv[23] = {25, 30,40,60, 80,100,150,200,250,300,
		       400,500,600,800,1000,1500, 2000, 3000, 4000, 5500,
		       6000, 7000,8000};
      Float_t y[23] = {5044, 4719, 3835, 3224, 2847, 2691,2327, 2067, 1937, 1781,
		       1560, 1456 ,1339, 1163.5, 1027, 819, 650, 520, 370.5, 234,
		       156, 78, 0};
      
      TGraph* gr = new TGraph(23,mv,y);
      fSlewingLED.AddAtAndExpand(gr,ipmt);
  }
//__________________________________________________________________

Float_t AliT0Parameters::GetSlewingLED(Int_t ipmt, Float_t mv) const
{
  if (!fCalibentry) {
    return ((TGraph*)fSlewingLED.At(ipmt))->Eval(mv); 
  } 
  return fgCalibData->GetSlewingLED(ipmt, mv) ;
}


//__________________________________________________________________

TGraph *AliT0Parameters::GetSlew(Int_t ipmt) const
{
  if (!fCalibentry) {
    return  (TGraph*)fSlewingLED.At(ipmt); 
  } 
  return fgCalibData -> GetSlew(ipmt) ;
}

//__________________________________________________________________


void 
AliT0Parameters::SetSlewingRec(Int_t ipmt)
{
  //  Set Slweing Correction for LED channel 
      Float_t mv[23] = {25, 30, 40,60, 80,100,150,200,250,300,
			400,500,600,800,1000,1500, 2000, 3000, 4000, 5500,
			6000, 7000,8000};
      Float_t y[23] = {5044, 4719, 3835, 3224, 2847, 2691,2327, 2067, 1937, 1781, 
		       1560, 1456 ,1339, 1163.5, 1027, 819, 650, 520, 370.5, 234, 
		       156, 78, 0};
      Float_t y1[23], mv1[23];
      for (Int_t i=0; i<23; i++){
	y1[i] = y[22-i]; mv1[i] = mv[22-i];}

      TGraph* gr = new TGraph(23,y1,mv1);
      fSlewingRec.AddAtAndExpand(gr,ipmt);

}
//__________________________________________________________________

Float_t AliT0Parameters::GetSlewingRec(Int_t ipmt, Float_t mv) const
{
  if (!fCalibentry) {
    return ((TGraph*)fSlewingRec.At(ipmt))->Eval(mv); 
  } 
  return fgCalibData -> GetSlewingRec(ipmt, mv) ;
}

//__________________________________________________________________

TGraph *AliT0Parameters::GetSlewRec(Int_t ipmt) const
{
  if (!fCalibentry) {
    return  (TGraph*)fSlewingRec.At(ipmt); 
  } 
  return fgCalibData -> GetSlewRec(ipmt) ;
}

//__________________________________________________________________
void 
AliT0Parameters::SetPMTeff(Int_t ipmt)
{
  Float_t lambda[50];
  Float_t eff[50 ] = {0,        0,       0.23619,  0.202909, 0.177913, 
		    0.175667, 0.17856, 0.190769, 0.206667, 0.230286,
		    0.252276, 0.256267,0.26,     0.27125,  0.281818,
		    0.288118, 0.294057,0.296222, 0.301622, 0.290421, 
		    0.276615, 0.2666,  0.248,    0.23619,  0.227814, 
		    0.219818, 0.206667,0.194087, 0.184681, 0.167917, 
		    0.154367, 0.1364,  0.109412, 0.0834615,0.0725283, 
		    0.0642963,0.05861, 0.0465,   0.0413333,0.032069, 
		    0.0252203,0.02066, 0.016262, 0.012,    0.00590476,
		    0.003875, 0.00190, 0,        0,        0          } ;
  for (Int_t i=0; i<50; i++) lambda[i]=200+10*i; 

  TGraph* gr = new TGraph(50,lambda,eff);
  fPMTeff.AddAtAndExpand(gr,ipmt);
}
