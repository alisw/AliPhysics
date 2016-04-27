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
#include "AliT0.h"
#include "AliLog.h"		  
#include "AliT0Parameters.h"	  
#include "AliT0CalibData.h"   
#include "AliT0CalibWalk.h"   
#include "AliT0CalibTimeEq.h"   
#include "AliT0CalibLatency.h"   
#include "AliT0LookUpKey.h"
#include "AliT0LookUpValue.h"
#include <AliCDBManager.h>        
#include <AliCDBEntry.h>          
#include <AliCDBStorage.h>  
#include <TMath.h>
#include <TSystem.h>
//#include <Riostream.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <TGeoMatrix.h>
#include <AliGeomManager.h>

AliT0CalibTimeEq* AliT0Parameters::fgCalibData = 0;
AliT0CalibData* AliT0Parameters::fgLookUp = 0;
AliT0CalibWalk* AliT0Parameters::fgSlewCorr =0;
AliT0CalibLatency *AliT0Parameters::fgLatency=0;
//====================================================================
using std::cout;
ClassImp(AliT0Parameters)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliT0Parameters* AliT0Parameters::fgInstance = 0;
//____________________________________________________________________
AliT0Parameters* AliT0Parameters::Instance() 
{
  // Get static instance 
  if (!fgInstance) {
    fgInstance = new AliT0Parameters;
  }
  return fgInstance;
}

//____________________________________________________________________
AliT0Parameters::AliT0Parameters()
  :fIsInit(kFALSE),
   fPh2Mip(0),fmV2Mip(0),
   fChannelWidth(0),fmV2Channel(0),
   fQTmin(0),fQTmax(0),
   fAmpLEDRec(0), 
   fPMTeff(),
   fWalk(0),
   fQTC(0),
   fAmpLED(0),
   fTimeDelayCFD(0), 
 //  fTimeV0(0), 
   fTimeDelayTVD(0),
   fMeanT0(512),
   fMeanVertex(0),
   fLatencyHPTDC(0),
   fLatencyL1(0),
   fLatencyL1A(0),
   fLatencyL1C(0),
   fLookUp(0),
   fNumberOfTRMs(2),
   fCalibentry(), 
   fLookUpentry(),
   fSlewCorr(),
   fLatency()
  
{
  // Default constructor 
  for (Int_t ipmt=0; ipmt<24; ipmt++)
    {
      SetPh2Mip();      
      SetmV2Mip();      
      SetChannelWidth();
      SetmV2channel();
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
   if (fIsInit) return;

   AliCDBManager *stor =AliCDBManager::Instance();
   //time equalizing
   fCalibentry  = stor->Get("T0/Calib/TimeDelay");
   if (fCalibentry)
     fgCalibData  = (AliT0CalibTimeEq*)fCalibentry->GetObject();
   else {
         AliFatal(" ALARM !!!! No time delays in CDB "); 
     fIsInit = kFALSE;
     return;
   }
 //slewing correction
  fSlewCorr  = stor->Get("T0/Calib/Slewing_Walk");
  if (fSlewCorr){
    fgSlewCorr  = (AliT0CalibWalk*)fSlewCorr->GetObject();
  }
  else {
      AliFatal(" ALARM !!!! No slewing correction in CDB "); 
    fIsInit = kFALSE;
    return;
  }
  //lookup table
  fLookUpentry  = stor->Get("T0/Calib/LookUp_Table");
  if (fLookUpentry){
    fgLookUp  = (AliT0CalibData*)fLookUpentry->GetObject();
  }
  else {
     AliFatal(" ALARM !!!! No Lookup table  in CDB "); 
    fIsInit = kFALSE;
    return;
  }
  //latency
  
 fLatency  = stor->Get("T0/Calib/Latency");
  if (fLatency){
    fgLatency  = (AliT0CalibLatency*)fLatency->GetObject();
  }
  else {
     AliWarning(" !!! no latency  in CDB "); 
    return;
  }
  
fIsInit = kTRUE;
}


//__________________________________________________________________

void AliT0Parameters::InitIfOnline()
{
// should be used in online
// for switching to this one should write
  // AliT0RawReader myrawreader(rawReader);
//	myrawreader.SetOnlineMode(kTRUE);
  
  if (fIsInit) return;
  //standart configuration (used for simulation)
  //Int_t trm=0; Int_t tdc=0; Int_t chain=0; Int_t channel=0;
  // configuration for test Jun07.
  fgLookUp = new AliT0CalibData("T0");
  
  fNumberOfTRMs = 2;
  fgLookUp-> SetNumberOfTRMs(fNumberOfTRMs);
  Int_t trm=7;
  Int_t tdc=0;  Int_t channel=0;
  Int_t ikey=0; Int_t chain=0;
  for (Int_t ik=0; ik<226; ik++)
    {
     if (ik==56) {trm=7; chain=1; tdc=0; channel=0;}
      if (ik==106) { trm=9; chain=0; tdc=0; channel=0;}
      if (ik==162) { trm=9; chain=1; tdc=0; channel=0;}
      if (ik==211) { trm=7; chain=1; tdc=14; channel=0;}
      if (ik==215) { trm=9; chain=1; tdc=12; channel=0;}
      AliT0LookUpKey * lookkey= new AliT0LookUpKey();
      AliT0LookUpValue * lookvalue= new AliT0LookUpValue();
      lookvalue->SetTRM(trm);
      lookvalue->SetTDC(tdc);
      lookvalue->SetChain(chain);
      lookvalue->SetChannel(channel);
      lookkey->SetKey(ik);
      fgLookUp->GetMapLookup()->Add((TObject*)lookvalue,(TObject*)lookkey);
      //       printf(" LookUp ik %i trm %i chain %i tdc %i  channel %i\n",ik, trm, chain, tdc, channel);	
      if (channel<6) channel +=2;
      else {channel = 0; tdc++;}
       ikey++;     
    }
  
  fIsInit=kTRUE;
}
//__________________________________________________________________
Float_t
AliT0Parameters::GetTimeDelayCFD(Int_t ipmt) 
  {
  // return time delay for CFD channel
   // 
  if (!fCalibentry) 
    {
      fTimeDelayCFD = 1000+ipmt*100;
      return fTimeDelayCFD;
    }
   
  return fgCalibData->GetTimeEq(ipmt);
}
//__________________________________________________________________
Float_t
AliT0Parameters::GetCFD(Int_t ipmt) 
  {
  // return  CFD channel
   
    return fgCalibData->GetCFDvalue(ipmt,0);
}
//__________________________________________________________________
Float_t
AliT0Parameters::GetQT1(Int_t ipmt) 
  {
  // return  CFD channel
   
    return fgCalibData->GetCFDvalue(ipmt,1);
}
//__________________________________________________________________
Float_t
AliT0Parameters::GetPedestalOld(Int_t ipmt) 
  {
  // return  CFD channel
   
    return fgCalibData->GetCFDvalue(ipmt,3);
}
//__________________________________________________________________
Float_t
AliT0Parameters::GetMeanOrA() 
  {
  // return  CFD channel
   
    return fgCalibData->GetCFDvalue(0,2);
}
Float_t
AliT0Parameters::GetMeanOrC() 
  {
  // return  CFD channel
   
    return fgCalibData->GetCFDvalue(1,2);
}
Float_t
AliT0Parameters::GetMeanTVDC() 
  {
  // return  CFD channel
   
    return fgCalibData->GetMeanVertex();
}

//__________________________________________________________________
Float_t
AliT0Parameters::GetLatencyHPTDC() 
  {
  // return LatencyHPTDC for CFD channel
  if (!fLatency) 
    {
      fLatencyHPTDC=9000.;
      return fLatencyHPTDC;
    }
   
  return fgLatency->GetLatencyHPTDC();
}
//__________________________________________________________________
Float_t
AliT0Parameters::GetLatencyL1() 
  {
  // return time delay for CFD channel
   
  return fgLatency->GetLatencyL1();
}

//__________________________________________________________________
Float_t
AliT0Parameters::GetLatencyL1A() 
  {
  // return time delay for CFD channel
   
  return fgLatency->GetLatencyL1A();
}

//__________________________________________________________________
Float_t
AliT0Parameters::GetLatencyL1C() 
  {
  // return time delay for CFD channel
   
  return fgLatency->GetLatencyL1C();
}
//__________________________________________________________________

Float_t
AliT0Parameters:: GetMeanVertex()
{ 
  if (!fCalibentry) 
    {
      fMeanVertex=0;
      return fMeanVertex;
    }
   
  return fgCalibData->GetMeanVertex();
}
//__________________________________________________________________

TGraph *AliT0Parameters::GetAmpLEDRec(Int_t ipmt) const
{
   if (!fSlewCorr) {
     AliError("No slewing correction is available!");
     return  (TGraph*)fAmpLEDRec.At(ipmt); 
  } 
  return fgSlewCorr -> GetAmpLEDRec(ipmt) ;
}

//__________________________________________________________________

TGraph *AliT0Parameters::GetWalk(Int_t ipmt) const
{
  if (!fSlewCorr) {
    AliError("No walk correction is available!");
    return  (TGraph*)fWalk.At(ipmt); 
  } 
  return fgSlewCorr -> GetWalk(ipmt) ;
}

//__________________________________________________________________

TGraph *AliT0Parameters::GetQTC(Int_t ipmt) const
{
  if (!fSlewCorr) {
    AliError("No walk correction is available!");
    //    return  (TGraph*)fQTC.At(ipmt); 
   return  0; 
  } 
  return fgSlewCorr -> GetQTC(ipmt) ;
}

//__________________________________________________________________
TGraph *AliT0Parameters::GetAmpLED(Int_t ipmt) const
{
  if (!fSlewCorr) {
    AliError("No walk correction is available!");
    //    return  (TGraph*)fQTC.At(ipmt); 
   return  0; 
  } 
  return fgSlewCorr -> GetAmpLED(ipmt) ;
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
//________________________________________________________________

Int_t 
AliT0Parameters::GetChannel(Int_t trm,  Int_t tdc, Int_t chain, Int_t channel)
{

  if (fgLookUp) {
    AliT0LookUpValue key(trm,tdc,chain,channel);
      AliT0LookUpKey *val = (AliT0LookUpKey*) fgLookUp->GetMapLookup()->GetValue((TObject*)&key);
      // AliT0LookUpKey *val = (AliT0LookUpKey*) fLookUp.GetValue((TObject*)&key);
    if (val )
      return val->GetKey();
    else {
      AliWarning(Form("No such address (%d %d %d %d)!",trm,tdc,chain,channel));
      return -1;
    }
  }
  else {
    AliError("No look up table has been loader!");
    return -1;
  }

}
//__________________________________________________________________
TMap *AliT0Parameters::GetMapLookup()
{
  if (!fgLookUp){
    cout<<" No look up table in OCDB";
    return 0;
  }
  return   fgLookUp->GetMapLookup();
}
//__________________________________________________________________

Int_t
AliT0Parameters::GetNumberOfTRMs() 
{
  // return number of trms
  // 
  if (!fgLookUp) {
    //  fNumberOfTRMs = 2;
    return  fNumberOfTRMs;
  } 
  return  fgLookUp ->GetNumberOfTRMs();
}
/*
//________________________________________________________________________________
Double_t AliT0Parameters::GetZPosition(const char* symname){
// Get the global z coordinate of the given T0 alignable volume
//
  Double_t *tr = AliGeomManager::GetMatrix(symname)->GetTranslation();

  return tr[2];
}
*/
//________________________________________________________________________________
Double_t AliT0Parameters::GetZPosition(const char* symname){
// Get the global z coordinate of the given T0 alignable volume
//
  Double_t *tr;
  TGeoPNEntry *pne = gGeoManager->GetAlignableEntry(symname);
  if (!pne) return 0;
  

  TGeoPhysicalNode *pnode = pne->GetPhysicalNode();
  if(pnode){
          TGeoHMatrix* hm = pnode->GetMatrix();
           tr = hm->GetTranslation();
  }else{
          const char* path = pne->GetTitle();
          if(!gGeoManager->cd(path)){
                  AliErrorClass(Form("Volume path %s not valid!",path));
                  return 0;
          }
         tr = gGeoManager->GetCurrentMatrix()->GetTranslation();
  }
  return tr[2];

}
//________________________________________________________________________________

Double_t AliT0Parameters::GetZPositionShift(const char* symname)
{
// Get the global z coordinate of the given T0 alignable volume
//
  Double_t *tr = AliGeomManager::GetMatrix(symname)->GetTranslation();

  TGeoHMatrix origmat;
  AliGeomManager::GetOrigGlobalMatrix(symname,origmat);
  Double_t *otr = origmat.GetTranslation();

  return (tr[2]-otr[2]);
}

