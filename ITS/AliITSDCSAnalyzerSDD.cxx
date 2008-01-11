/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the class for SDD DCS data analysis         //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////


#include "AliITSDCSAnalyzerSDD.h"
#include "AliDCSValue.h"



const Int_t AliITSDCSAnalyzerSDD::fgkNcathodes = 291;
const Float_t AliITSDCSAnalyzerSDD::fgkCathodePitch = 0.0120;
const Int_t AliITSDCSAnalyzerSDD::fgkTemperatureStatusOK = 4;

ClassImp(AliITSDCSAnalyzerSDD)

//---------------------------------------------------------------
AliITSDCSAnalyzerSDD::AliITSDCSAnalyzerSDD():
TObject(),
fDCSData(0)
{
  // Default constructor
  Init();
  fDCSData=new AliITSDCSDataSDD*[kNmodules];
}
//---------------------------------------------------------------
AliITSDCSAnalyzerSDD::AliITSDCSAnalyzerSDD(const AliITSDCSAnalyzerSDD& /* dcsa */):
TObject(),
fDCSData(0)
{
  // copy constructor
  // Copies are not allowed. The method is protected to avoid misuse.
  AliError("Copy constructor not allowed");
}
//---------------------------------------------------------------
AliITSDCSAnalyzerSDD& AliITSDCSAnalyzerSDD::operator=(const AliITSDCSAnalyzerSDD& /* dcsa */){
  // assigment operator
  // Assignment is not allowed. The method is protected to avoid misuse.
  AliError("Assignment operator not allowed");
  return *this;
}
//---------------------------------------------------------------
AliITSDCSAnalyzerSDD::~AliITSDCSAnalyzerSDD(){
  // destructor
  for(int j=0; j<kNmodules; j++){
    if(fDCSData[j]) delete fDCSData[j];
  }
  delete [] fDCSData;
}
//---------------------------------------------------------------
void AliITSDCSAnalyzerSDD::AnalyzeData(TMap* dcsMap){
  // Data processing
   
  for(int j=0; j<kNmodules; j++){
    TObjArray* arrHV = (TObjArray*) dcsMap->GetValue(fHVDPNames[j].Data());   
    if(!arrHV){
      AliError(Form("DCS HV alias %s not found!", fHVDPNames[j].Data()));
      continue;
    }
    TObjArray* arrMV = (TObjArray*) dcsMap->GetValue(fMVDPNames[j].Data());   
    if(!arrMV){
      AliError(Form("DCS MV alias %s not found!", fMVDPNames[j].Data()));
      continue;
    }
    TObjArray* arrTL = (TObjArray*) dcsMap->GetValue(fTLDPNames[j].Data());   
    if(!arrTL){
      AliError(Form("DCS TEMP_L alias %s not found!", fTLDPNames[j].Data()));
      continue;
    }
    TObjArray* arrTR = (TObjArray*) dcsMap->GetValue(fTRDPNames[j].Data());   
    if(!arrTR){
      AliError(Form("DCS TEMP_R alias %s not found!", fTRDPNames[j].Data()));
      continue;
    }
    TObjArray* arrStTL = (TObjArray*) dcsMap->GetValue(fTLStDPNames[j].Data());   
    if(!arrStTL){
      AliError(Form("DCS TEMP_L_STATE alias %s not found!", fTLStDPNames[j].Data()));
      continue;
    }
    TObjArray* arrStTR = (TObjArray*) dcsMap->GetValue(fTRStDPNames[j].Data());   
    if(!arrStTR){
      AliError(Form("DCS TEMP_R_STATE alias %s not found!", fTRStDPNames[j].Data()));
      continue;
    }

    Int_t nData=arrHV->GetEntries();
    fDCSData[j]=new AliITSDCSDataSDD(nData);
    for(int i =0; i<arrHV->GetEntries(); i++) {
      AliDCSValue* valHV = (AliDCSValue*) arrHV->At(i);
      AliDCSValue* valMV = (AliDCSValue*) arrMV->At(i);
      AliDCSValue* valTL = (AliDCSValue*) arrTL->At(i);
      AliDCSValue* valTR = (AliDCSValue*) arrTR->At(i);
      AliDCSValue* valStTL = (AliDCSValue*) arrStTL->At(i);
      AliDCSValue* valStTR = (AliDCSValue*) arrStTR->At(i);
      Int_t timeStamp = valHV->GetTimeStamp();
      Float_t hv = valHV->GetFloat();
      Float_t mv = valMV->GetFloat();
      Float_t edrift = (hv-mv)/(Float_t)fgkNcathodes/fgkCathodePitch;
      Int_t statusTL=valStTL->GetInt();
      Int_t statusTR=valStTR->GetInt();
      Float_t tleft = -9999.;
      if(statusTL==fgkTemperatureStatusOK) tleft=valTL->GetFloat();
      Float_t tright = -9999.;
      if(statusTR==fgkTemperatureStatusOK) tright=valTR->GetFloat();
      fDCSData[j]->SetValues(timeStamp,edrift,tleft,tright);
    }
  }
}

//---------------------------------------------------------------
void AliITSDCSAnalyzerSDD::Init(){
  // Initialization of DCS DP names
  Char_t dpName[50];
  Char_t modName[50];
  for(Int_t iLad=0; iLad<kNladders3; iLad++){
    for(Int_t iMod=0; iMod<kNmodLad3;iMod++){
      sprintf(modName,"SDD_LAYER3_LADDER%02d_MODULE%d",iLad,iMod);
      Int_t id=iMod+iLad*kNmodLad3;
      sprintf(dpName,"%s_HV",modName);
      fHVDPNames[id]=dpName;
      sprintf(dpName,"%s_MV",modName);
      fMVDPNames[id]=dpName;
      sprintf(dpName,"%s_TEMP_L",modName);
      fTLDPNames[id]=dpName;
      sprintf(dpName,"%s_TEMP_R",modName);
      fTRDPNames[id]=dpName;
      sprintf(dpName,"%s_TEMP_L_STATE",modName);
      fTLStDPNames[id]=dpName;
      sprintf(dpName,"%s_TEMP_R_STATE",modName);
      fTRStDPNames[id]=dpName;
    }
  }
  for(Int_t iLad=0; iLad<kNladders4; iLad++){
    for(Int_t iMod=0; iMod<kNmodLad4;iMod++){
      sprintf(modName,"SDD_LAYER4_LADDER%02d_MODULE%d",iLad,iMod);
      Int_t id=kNladders3*kNmodLad3+iMod+iLad*kNmodLad4;
      sprintf(dpName,"%s_HV",modName);
      fHVDPNames[id]=dpName;
      sprintf(dpName,"%s_MV",modName);
      fMVDPNames[id]=dpName;
      sprintf(dpName,"%s_TEMP_L",modName);
      fTLDPNames[id]=dpName;
      sprintf(dpName,"%s_TEMP_R",modName);
      fTRDPNames[id]=dpName;
      sprintf(dpName,"%s_TEMP_L_STATE",modName);
      fTLStDPNames[id]=dpName;
      sprintf(dpName,"%s_TEMP_R_STATE",modName);
      fTRStDPNames[id]=dpName;
    }
  } 
}

//---------------------------------------------------------------
void AliITSDCSAnalyzerSDD::PrintDCSDPNames(){
  // Data processing
  for(int j=0; j<kNmodules; j++){
    printf("Module %d      %s   %s   %s   %s\n",j,fHVDPNames[j].Data(),fMVDPNames[j].Data(),fTLDPNames[j].Data(),fTRDPNames[j].Data());
  }
}
