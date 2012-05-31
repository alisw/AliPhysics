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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class to set HMPID reconstruction parameters (normal, HTA, UserCut ...    //
//                                                                           //  
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
//
//Email: domenico.dibari@ba.infn.it
//
#include "AliLog.h"
#include "AliHMPIDRecoParam.h"
#include "AliHMPIDParam.h"

ClassImp(AliHMPIDRecoParam)

//_____________________________________________________________________________
AliHMPIDRecoParam::AliHMPIDRecoParam():AliDetectorRecoParam(),   
  fHmpRecoMode(kTRUE),fHmpFixedDistCut(kTRUE),
  fHmpTrackMatchingDist(1.0)
{
  //
  // ctor
  //
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++) fHmpUserCut[iCh]=3;
  for(Int_t iPol=0;iPol<5;iPol++) fHmpTrackMatchingDistParas[iPol]=0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDRecoParam::AliHMPIDRecoParam(const AliHMPIDRecoParam &p):AliDetectorRecoParam(p),   
    fHmpRecoMode(kTRUE),fHmpFixedDistCut(kTRUE),
    fHmpTrackMatchingDist(1.0)
{ 
   //copy Ctor

   fHmpRecoMode= p.fHmpRecoMode;
   fHmpFixedDistCut=p.fHmpFixedDistCut;
   for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++) fHmpUserCut[iCh]=p.fHmpUserCut[iCh];
   fHmpTrackMatchingDist=p.fHmpTrackMatchingDist;
   for(Int_t iPol=0;iPol<5;iPol++) fHmpTrackMatchingDistParas[iPol]=p.fHmpTrackMatchingDistParas[iPol];
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDRecoParam& AliHMPIDRecoParam::operator=(const AliHMPIDRecoParam &p)
{
//
// assign. operator
//
  if(this!=&p){
    AliDetectorRecoParam::operator=(p);
    this->fHmpRecoMode= p.fHmpRecoMode;
    this->fHmpFixedDistCut=p.fHmpFixedDistCut;
    for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++) this->fHmpUserCut[iCh] = p.fHmpUserCut[iCh];   
    this->fHmpTrackMatchingDist=p.fHmpTrackMatchingDist;
    for(Int_t iPol=0;iPol<5;iPol++) this->fHmpTrackMatchingDistParas[iPol]=p.fHmpTrackMatchingDistParas[iPol];
  }
  return *this;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDRecoParam::~AliHMPIDRecoParam() 
{
  //
  // dtor
  //  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDRecoParam *AliHMPIDRecoParam::GetLowFluxParam(){
  //
  // Set HMPID Reco Params for Low Flux environment
  //
  AliHMPIDRecoParam *hmpParam = new AliHMPIDRecoParam;
  hmpParam->fHmpRecoMode = kTRUE;                                                     //kTRUE = normal reco. kFLASE = HTA
  hmpParam->fHmpUserCut[0] = 3;                                                       //HMPID Module 0 User DAQ Sigma cut
  hmpParam->fHmpUserCut[1] = 3;                                                       //HMPID Module 1 User DAQ Sigma cut
  hmpParam->fHmpUserCut[2] = 3;                                                       //HMPID Module 2 User DAQ Sigma cut
  hmpParam->fHmpUserCut[3] = 3;                                                       //HMPID Module 3 User DAQ Sigma cut
  hmpParam->fHmpUserCut[4] = 3;                                                       //HMPID Module 4 User DAQ Sigma cut
  hmpParam->fHmpUserCut[5] = 3;                                                       //HMPID Module 5 User DAQ Sigma cut
  hmpParam->fHmpUserCut[6] = 3;                                                       //HMPID Module 6 User DAQ Sigma cut
  hmpParam->fHmpFixedDistCut=kTRUE;                                                   //HMPID fixed (kTRUE) or parameterized distance cut (kFALSE)
  hmpParam->fHmpTrackMatchingDist = 3.0;                                              //HMPID Track Matching distance cut
  for(Int_t iPol=0;iPol<5;iPol++) hmpParam->fHmpTrackMatchingDistParas[iPol]=1;       //Prevision for momentum dependen track matching
  hmpParam->SetName("HMP Low Flux");
  hmpParam->SetTitle("HMP Low Flux");
  return hmpParam;
    
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDRecoParam *AliHMPIDRecoParam::GetHighFluxParam(){
  //
  // Set HMPID Reco Params for Low Flux environment
  //
  AliHMPIDRecoParam *hmpParam = new AliHMPIDRecoParam;
  hmpParam->fHmpRecoMode = kTRUE;                                                       //kTRUE = normal reco. kFLASE = HTA
  hmpParam->fHmpUserCut[0] = 3;                                                         //HMPID Module 0 User DAQ Sigma cut
  hmpParam->fHmpUserCut[1] = 3;                                                         //HMPID Module 1 User DAQ Sigma cut
  hmpParam->fHmpUserCut[2] = 3;                                                         //HMPID Module 2 User DAQ Sigma cut
  hmpParam->fHmpUserCut[3] = 3;                                                         //HMPID Module 3 User DAQ Sigma cut
  hmpParam->fHmpUserCut[4] = 3;                                                         //HMPID Module 4 User DAQ Sigma cut
  hmpParam->fHmpUserCut[5] = 3;                                                         //HMPID Module 5 User DAQ Sigma cut
  hmpParam->fHmpUserCut[6] = 3;                                                         //HMPID Module 6 User DAQ Sigma cut
  hmpParam->fHmpFixedDistCut=kTRUE;                                                     //HMPID fixed (kTRUE) or parameterized distance cut (kFALSE)
  hmpParam->fHmpTrackMatchingDist = 1.0;                                                //HMPID Track Matching distance cut
  for(Int_t iPol=0;iPol<5;iPol++) hmpParam->fHmpTrackMatchingDistParas[iPol]=1;         //Prevision for momentum dependen track matching
  hmpParam->SetName("HMP High Flux");
  hmpParam->SetTitle("HMP High Flux");
  return hmpParam;
    
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDRecoParam *AliHMPIDRecoParam::GetCosmicParam(){
  //
  // Set HMPID Reco Params for Low Flux environment
  //
  AliHMPIDRecoParam *hmpParam = new AliHMPIDRecoParam;
  hmpParam->fHmpRecoMode = kTRUE;                                                                   //kTRUE = normal reco. kFLASE = HTA
  hmpParam->fHmpUserCut[0] = 3;                                                                     //HMPID Module 0 User DAQ Sigma cut
  hmpParam->fHmpUserCut[1] = 3;                                                                     //HMPID Module 1 User DAQ Sigma cut
  hmpParam->fHmpUserCut[2] = 3;                                                                     //HMPID Module 2 User DAQ Sigma cut
  hmpParam->fHmpUserCut[3] = 3;                                                                     //HMPID Module 3 User DAQ Sigma cut
  hmpParam->fHmpUserCut[4] = 3;                                                                     //HMPID Module 4 User DAQ Sigma cut
  hmpParam->fHmpUserCut[5] = 3;                                                                     //HMPID Module 5 User DAQ Sigma cut
  hmpParam->fHmpUserCut[6] = 3;                                                                     //HMPID Module 6 User DAQ Sigma cut
  hmpParam->fHmpFixedDistCut=kTRUE;                                                                 //HMPID fixed (kTRUE) or parameterized distance cut (kFALSE)
  hmpParam->fHmpTrackMatchingDist = 4.0;                                                            //HMPID Track Matching distance cut
  for(Int_t iPol=0;iPol<5;iPol++) hmpParam->fHmpTrackMatchingDistParas[iPol]=1;                     //Prevision for momentum dependen track matching
  hmpParam->SetName("HMP Cosmic ");
  hmpParam->SetTitle("HMP Cosmic");
  return hmpParam;
    
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRecoParam::PrintParameters() const
{
  //
  // Printing of the used HMPID reconstruction parameters
  //
   AliInfo(Form("%s",AliHMPIDRecoParam::GetName()));
   AliInfo(Form("IsDefault: %d",AliHMPIDRecoParam::IsDefault()));
   AliInfo(Form(" Running HMPID Reco: %d (1=Standard, 0=HTA)",fHmpRecoMode));
   AliInfo(Form(" Is track matching distance fixed (1) or momentum dependent (0): %d",fHmpFixedDistCut));
   AliInfo(Form(" HMPID track matching distance cut: %.3lf",fHmpTrackMatchingDist));
   for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++)
    AliInfo(Form(" HMPID Chamber: %d User DAQ Sigma cut: %d",iCh,fHmpUserCut[iCh]));
   for(Int_t iPol=0;iPol<5;iPol++) 
     AliInfo(Form(" HMPID momentum dependent distnce parameters: param[%d]=%lf",iPol,fHmpTrackMatchingDistParas[iPol]));
  
}
