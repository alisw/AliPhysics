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

/* $Id: AliPIDResponse.cxx 46193 2010-12-21 09:00:14Z wiechula $ */

//-----------------------------------------------------------------
//        Base class for handling the pid response               //
//        functions of all detectors                             //
//        and give access to the nsigmas                         //
//                                                               //
//   Origin: Jens Wiechula, Uni Tuebingen, jens.wiechula@cern.ch //
//-----------------------------------------------------------------

#include <TList.h>
#include <TObjArray.h>
#include <TPRegexp.h>
#include <TF1.h>
#include <TSpline.h>
#include <TFile.h>

#include <AliVEvent.h>
#include <AliLog.h>
#include <AliPID.h>

#include "AliPIDResponse.h"

ClassImp(AliPIDResponse);

AliPIDResponse::AliPIDResponse(Bool_t isMC/*=kFALSE*/) :
TNamed("PIDResponse","PIDResponse"),
fITSResponse(isMC),
fTPCResponse(),
fTRDResponse(),
fTOFResponse(),
fIsMC(isMC),
fOADBPath(),
fBeamType("PP"),
fLHCperiod(),
fMCperiodTPC(),
fRecoPass(0),
fRun(0),
fOldRun(0),
fArrPidResponseMaster(0x0),
fResolutionCorrection(0x0),
fTOFTimeZeroType(kBest_T0),
fTOFres(100.)
{
  //
  // default ctor
  //
  AliLog::SetClassDebugLevel("AliPIDResponse",10);
  AliLog::SetClassDebugLevel("AliESDpid",10);
  AliLog::SetClassDebugLevel("AliAODpidUtil",10);
}

//______________________________________________________________________________
AliPIDResponse::~AliPIDResponse()
{
  //
  // dtor
  //
  delete fArrPidResponseMaster;
}

//______________________________________________________________________________
AliPIDResponse::AliPIDResponse(const AliPIDResponse &other) :
TNamed(other),
fITSResponse(other.fITSResponse),
fTPCResponse(other.fTPCResponse),
fTRDResponse(other.fTRDResponse),
fTOFResponse(other.fTOFResponse),
fIsMC(other.fIsMC),
fOADBPath(other.fOADBPath),
fBeamType("PP"),
fLHCperiod(),
fMCperiodTPC(),
fRecoPass(0),
fRun(0),
fOldRun(0),
fArrPidResponseMaster(0x0),
fResolutionCorrection(0x0),
fTOFTimeZeroType(AliPIDResponse::kBest_T0),
fTOFres(100.)
{
  //
  // copy ctor
  //
}

//______________________________________________________________________________
AliPIDResponse& AliPIDResponse::operator=(const AliPIDResponse &other)
{
  //
  // copy ctor
  //
  if(this!=&other) {
    delete fArrPidResponseMaster;
    TNamed::operator=(other);
    fITSResponse=other.fITSResponse;
    fTPCResponse=other.fTPCResponse;
    fTRDResponse=other.fTRDResponse;
    fTOFResponse=other.fTOFResponse;
    fOADBPath=other.fOADBPath;
    fIsMC=other.fIsMC;
    fBeamType="PP";
    fLHCperiod="";
    fMCperiodTPC="";
    fRecoPass=0;
    fRun=0;
    fOldRun=0;
    fArrPidResponseMaster=0x0;
    fResolutionCorrection=0x0;
    fTOFTimeZeroType=AliPIDResponse::kBest_T0;
    fTOFres=100.;
  }
  return *this;
}

//______________________________________________________________________________
void AliPIDResponse::InitialiseEvent(AliVEvent *event, Int_t pass)
{
  //
  // Apply settings for the current event
  //
  fRecoPass=pass;

  if (!event) return;
  fRun=event->GetRunNumber();
  
  if (fRun!=fOldRun){
    ExecNewRun();
    fOldRun=fRun;
  }
  
  //TPC resolution parametrisation PbPb
  if ( fResolutionCorrection ){
    Double_t corrSigma=fResolutionCorrection->Eval(GetTPCMultiplicityBin(event));
    fTPCResponse.SetSigma(3.79301e-03*corrSigma, 2.21280e+04);
  }
  
  //TOF resolution
  SetTOFResponse(event, (AliPIDResponse::EStartTimeType_t)fTOFTimeZeroType);
  
}

//______________________________________________________________________________
void AliPIDResponse::ExecNewRun()
{
  //
  // Things to Execute upon a new run
  //
  SetRecoInfo();
  
  SetITSParametrisation();
  
  SetTPCPidResponseMaster();
  SetTPCParametrisation();
  
  fTOFResponse.SetTimeResolution(fTOFres);
}

//_____________________________________________________
Double_t AliPIDResponse::GetTPCMultiplicityBin(const AliVEvent * const event)
{
  //
  // Get TPC multiplicity in bins of 150
  //
  
  const AliVVertex* vertexTPC = event->GetPrimaryVertex();
  Double_t tpcMulti=0.;
  if(vertexTPC){
    Double_t vertexContribTPC=vertexTPC->GetNContributors();
    tpcMulti=vertexContribTPC/150.;
    if (tpcMulti>20.) tpcMulti=20.;
  }
  
  return tpcMulti;
}

//______________________________________________________________________________
void AliPIDResponse::SetRecoInfo()
{
  //
  // Set reconstruction information
  //
  
  //reset information
  fLHCperiod="";
  fMCperiodTPC="";
  
  fBeamType="";
    
  fBeamType="PP";
  
  //find the period by run number (UGLY, but not stored in ESD and AOD... )
  if (fRun>=114737&&fRun<=117223)      { fLHCperiod="LHC10B"; fMCperiodTPC="LHC10D1";  }
  else if (fRun>=118503&&fRun<=121040) { fLHCperiod="LHC10C"; fMCperiodTPC="LHC10D1";  }
  else if (fRun>=122195&&fRun<=126437) { fLHCperiod="LHC10D"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=127719&&fRun<=130850) { fLHCperiod="LHC10E"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=133004&&fRun<=135029) { fLHCperiod="LHC10F"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=135654&&fRun<=136377) { fLHCperiod="LHC10G"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=136851&&fRun<=139517) { fLHCperiod="LHC10H"; fMCperiodTPC="LHC10H8"; fBeamType="PBPB"; }
  else if (fRun>=139699) { fLHCperiod="LHC11A"; fMCperiodTPC="LHC10F6A"; }
  
}

//______________________________________________________________________________
void AliPIDResponse::SetITSParametrisation()
{
  //
  // Set the ITS parametrisation
  //
}

//______________________________________________________________________________
void AliPIDResponse::SetTPCPidResponseMaster()
{
  //
  // Load the TPC pid response functions from the OADB
  //
  
  //reset the PID response functions
  delete fArrPidResponseMaster;
  fArrPidResponseMaster=0x0;
  
  TString fileName(Form("%s/COMMON/PID/data/TPCPIDResponse.root", fOADBPath.Data()));
  
  TFile f(fileName.Data());
  if (f.IsOpen() && !f.IsZombie()){
    fArrPidResponseMaster=dynamic_cast<TObjArray*>(f.Get("TPCPIDResponse"));
    f.Close();
  }
  
  if (!fArrPidResponseMaster){
    AliFatal(Form("Could not retrieve the TPC pid response from: %s",fileName.Data()));
    return;
  }
  fArrPidResponseMaster->SetOwner();
}

//______________________________________________________________________________
void AliPIDResponse::SetTPCParametrisation()
{
  //
  // Change BB parametrisation for current run
  //
  
  if (fLHCperiod.IsNull()) {
    AliFatal("No period set, not changing parametrisation");
    return;
  }
  
  //
  // Set default parametrisations for data and MC
  //
  
  //data type
  TString datatype="DATA";
  //in case of mc fRecoPass is per default 1
  if (fIsMC) {
    datatype="MC";
    fRecoPass=1;
  }
  
  //
  //reset old splines
  //
  for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec){
    fTPCResponse.SetResponseFunction((AliPID::EParticleType)ispec,0x0);
  }
  
  //
  //set the new PID splines
  //
  TString period=fLHCperiod;
  if (fArrPidResponseMaster){
    TObject *grAll=0x0;
    //for MC don't use period information
//     if (fIsMC) period="[A-Z0-9]*";
    //for MC use MC period information
    if (fIsMC) period=fMCperiodTPC;
//pattern for the default entry (valid for all particles)
    TPRegexp reg(Form("TSPLINE3_%s_([A-Z]*)_%s_PASS%d_%s_MEAN",datatype.Data(),period.Data(),fRecoPass,fBeamType.Data()));
    
    //loop over entries and filter them
    for (Int_t iresp=0; iresp<fArrPidResponseMaster->GetEntriesFast();++iresp){
      TObject *responseFunction=fArrPidResponseMaster->At(iresp);
      if (responseFunction==0x0) continue;
      TString responseName=responseFunction->GetName();
      
      if (!reg.MatchB(responseName)) continue;
      
      TObjArray *arr=reg.MatchS(responseName);
      TString particleName=arr->At(1)->GetName();
      delete arr;
      if (particleName.IsNull()) continue;
      if (particleName=="ALL") grAll=responseFunction;
      else {
        //find particle id
        for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec){
          TString particle=AliPID::ParticleName(ispec);
          particle.ToUpper();
          if ( particle == particleName ){
            fTPCResponse.SetResponseFunction((AliPID::EParticleType)ispec,responseFunction);
            fTPCResponse.SetUseDatabase(kTRUE);
            AliInfo(Form("Adding graph: %d - %s",ispec,responseFunction->GetName()));
            break;
          }
        }
      }
    }
    
    //set default response function to all particles which don't have a specific one
    if (grAll){
      for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec){
        if (!fTPCResponse.GetResponseFunction((AliPID::EParticleType)ispec)){
          fTPCResponse.SetResponseFunction((AliPID::EParticleType)ispec,grAll);
          AliInfo(Form("Adding graph: %d - %s",ispec,grAll->GetName()));
        }
      }
    }
  }
  
  //
  // Setup resolution parametrisation
  //
  
  //default
  fTPCResponse.SetSigma(3.79301e-03, 2.21280e+04);
  
  if (fRun>=122195){
    fTPCResponse.SetSigma(2.30176e-02, 5.60422e+02);
  }
  
  fResolutionCorrection=(TF1*)fArrPidResponseMaster->FindObject(Form("TF1_%s_ALL_%s_PASS%d_%s_SIGMA",datatype.Data(),period.Data(),fRecoPass,fBeamType.Data()));
  
  if (fResolutionCorrection) AliInfo(Form("Setting multiplicity correction function: %s",fResolutionCorrection->GetName()));
}

