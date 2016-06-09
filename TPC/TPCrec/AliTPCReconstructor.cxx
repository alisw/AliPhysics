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

//--------------------------------------------------------------------
//          Options for the TPC Reconstruction in rec.C
//
//  4 options can be set to change the input for TPC reconstruction
//  which overwrites the usage of fUseHLTClusters of the AliTPCRecoParam
//
//  1) useRAW        - use RAW, if not present -> do nothing
//  2) useRAWorHLT   - use RAW, if not present -> use HLT clusters
//  3) useHLT        - use HLT clusters, if not present -> do nothing
//  4) useHLTorRAW   - use HLT clusters, if not present -> use RAW
//
//  -> The current default is useHLTorRAW
//--------------------------------------------------------------------

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for TPC reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TFile.h>

#include <AliLog.h>
#include <AliPID.h>
#include <AliESDpid.h>
#include <AliTPCPIDResponse.h>
#include "AliTPCReconstructor.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliRawReader.h"
#include "AliTPCclusterer.h"
#include "AliTPCtracker.h"
#include "AliTPCParam.h"
#include "AliTPCParamSR.h"
#include "AliTPCcalibDB.h"
#include "AliTracker.h"
#include "AliMagF.h"
#include "TTreeStream.h"

ClassImp(AliTPCReconstructor)


Int_t    AliTPCReconstructor::fgStreamLevel     = 0;             // stream (debug) level
AliTPCAltroEmulator *  AliTPCReconstructor::fAltroEmulator=0;    // ALTRO emulator
TTreeSRedirector    *  AliTPCReconstructor::fgDebugStreamer=0;                          // NOTE -  AliTPCReconstructor is not an owner of the streamer
TString                AliTPCReconstructor::fgPIDRespnonsePath="$ALICE_PHYSICS/OADB/COMMON/PID/data/TPCPIDResponse.root";
 
TVectorD *  AliTPCReconstructor::fSystematicErrors=0;
TVectorD *  AliTPCReconstructor::fSystematicErrorClusters=0;
TVectorD *  AliTPCReconstructor::fgExtendedRoads=0;
TVectorD *  AliTPCReconstructor::fgPrimaryDCACut=0;
Double_t    AliTPCReconstructor::fgPrimaryZ2XCut = 0;
Double_t    AliTPCReconstructor::fgZOutSectorCut = 0;
Bool_t      AliTPCReconstructor::fgCompactClusters = kFALSE;

AliTPCReconstructor::AliTPCReconstructor():
AliReconstructor(),
fClusterer(NULL),
fArrSplines(NULL)
{
  //
  // default constructor
  //
  //
  //
  AliTPCcalibDB * calib = AliTPCcalibDB::Instance();
  const AliMagF * field = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  calib->SetExBField(field);
  AliTPCParam* param = GetTPCParam();
  if (!param) {
    AliWarning("Loading default TPC parameters !");
    param = new AliTPCParamSR;
  }
  fClusterer = new AliTPCclusterer(param);
}

AliTPCReconstructor::AliTPCReconstructor(const AliTPCReconstructor& /*rec*/):
AliReconstructor(),
fClusterer(NULL),
fArrSplines(NULL)
{
  //
  // Dummy copu constructor
  //
}

AliTPCReconstructor& AliTPCReconstructor::operator=(const AliTPCReconstructor&){
  //
  // dummy operator
  //
  return *this;
}

//_____________________________________________________________________________
AliTPCReconstructor::~AliTPCReconstructor()
{
  if (fClusterer)   delete fClusterer;
  delete fArrSplines;
  delete AliTPCcalibDB::Instance();
}

//_____________________________________________________________________________
void AliTPCReconstructor::Reconstruct(TTree* digitsTree, TTree* clustersTree) const {
  // single event local reconstruction
  // of TPC data
  fClusterer->SetInput(digitsTree);
  fClusterer->SetOutput(clustersTree);
  fClusterer->Digits2Clusters();
}

//_____________________________________________________________________________
void AliTPCReconstructor::Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const {
  // single event local reconstruction
  // of TPC data starting from raw data

  fClusterer->SetOutput(clustersTree);
  fClusterer->Digits2Clusters(rawReader);
}

//_____________________________________________________________________________
AliTracker* AliTPCReconstructor::CreateTracker() const
{
// create a TPC tracker

  AliTPCParam* param = GetTPCParam();
  if (!param) {
    AliWarning("Loading default TPC parameters !");
    param = new AliTPCParamSR;
  }
  param->ReadGeoMatrices();
  
  AliTPCtracker* tracker = new AliTPCtracker(param);

  ParseOptions(tracker);

  return tracker;
}

//_____________________________________________________________________________
void AliTPCReconstructor::FillESD(TTree */*digitsTree*/, TTree */*clustersTree*/,
				  AliESDEvent* /*esd*/) const
{
// make PID
/*  Now done in AliESDpid
  Double_t parTPC[] = {50., 0.07, 5.};  // MIP nnormalized to channel 50 -MI
  AliTPCpidESD tpcPID(parTPC);
  tpcPID.MakePID(esd);
*/
}


//_____________________________________________________________________________
AliTPCParam* AliTPCReconstructor::GetTPCParam() const
{
// get the TPC parameters

  AliTPCParam* param = AliTPCcalibDB::Instance()->GetParameters();

  return param;
}

//_____________________________________________________________________________
void AliTPCReconstructor::SetSplinesFromOADB(const char* tmplt, AliESDpid *esdPID)
{
  //
  //  load splines from the OADB using 'template'
  //

  // only load splines if not already set
  if (!fArrSplines) {
    fArrSplines=new TObjArray(Int_t(AliPID::kSPECIES));
    fArrSplines->SetOwner();
    TString stemplate(tmplt);

    TFile f(GetPIDRespnonsePath());

    TObjArray *arrPidResponseMaster=0x0;

    if (f.IsOpen() && !f.IsZombie()){
      arrPidResponseMaster=dynamic_cast<TObjArray*>(f.Get("TPCPIDResponse"));
    }
    f.Close();

    if (!arrPidResponseMaster){
      AliError("PID response array not found, cannot assign proper splines");
      return;
    }

    for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec)
    {
      Int_t ispec2=ispec;
      if (ispec==Int_t(AliPID::kMuon)) ispec2=Int_t(AliPID::kPion);

      TString particle=AliPID::ParticleName(ispec2);
      particle.ToUpper();

      TString splineName;
      splineName.Form(stemplate.Data(),particle.Data());
      TObject *spline=arrPidResponseMaster->FindObject(splineName.Data());
      if (!spline) {
        AliError(Form("No spline found for '%s'", splineName.Data()));
        continue;
      };
      AliInfo(Form("Adding Response function %d:%s",ispec,splineName.Data()));
      fArrSplines->AddAt(spline->Clone(), ispec);
    }    
    arrPidResponseMaster->Delete();
    delete arrPidResponseMaster;
    if (fArrSplines->GetEntries()!=Int_t(AliPID::kSPECIES)) {
      AliError("Splines not found for all species, cannot use proper PID");
      delete fArrSplines;
      fArrSplines=NULL;
      return;
    }
  }

  for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec)
  {
    esdPID->GetTPCResponse().SetResponseFunction( (AliPID::EParticleType)ispec, fArrSplines->UncheckedAt(ispec) );
  }

  esdPID->GetTPCResponse().SetUseDatabase(kTRUE);
}

//_____________________________________________________________________________
void AliTPCReconstructor::GetPidSettings(AliESDpid *esdPID)
{
  //
  // Get TPC pid splines. They should be written to the OCDB during the CPass
  // the splines themselves are owned by the OCDB object
  //

  // parse options
  TString allopt(GetOption());
  TObjArray *optArray=allopt.Tokenize(";");

  // defines whether the pid was set via a specific option in the rec.C
  Bool_t pidSetInOptions = kFALSE;
  
  for (Int_t iopt=0; iopt<optArray->GetEntriesFast(); ++iopt){
    if (!optArray->At(iopt)) continue;
    TString option(static_cast<TObjString*>(optArray->At(iopt))->GetString().Strip(TString::kBoth,' '));

    if (!option.BeginsWith("PID.")) continue;

    // remove 'PID.' identifyer
    option.Remove(0,4);

    // parse PID type
    if (option.BeginsWith("Static=")){
      option.Remove(0,option.First('=')+1);
      if (option.Contains("LHC13b2_fix_PID")) {
        esdPID->GetTPCResponse().SetBetheBlochParameters(0.0320981, 19.9768, 2.52666e-16, 2.72123, 6.08092);
        esdPID->GetTPCResponse().SetMip(53.4968);
        pidSetInOptions=kTRUE;
      }
      
    } else if (option.BeginsWith("OADB=")) {
      option.Remove(0,option.First('=')+1);
      AliInfo(Form("Setting splines From OADB using template: '%s'",option.Data()));
      SetSplinesFromOADB(option, esdPID);
      pidSetInOptions=kTRUE;
    } else if (option.BeginsWith("OCDB=")){
      option.Remove(0,option.First('=')+1);
      // not yet implemented
    }
    
  }

  delete optArray;

  //
  // Initialisation of BB parameters from the OCDB.
  // They are stored in the AliTPCParam
  //
  if (!pidSetInOptions) {
    AliTPCParam* param = AliTPCcalibDB::Instance()->GetParameters();
    if (param) {
      TVectorD *paramBB=param->GetBetheBlochParameters();
      if (paramBB){
        esdPID->GetTPCResponse().SetBetheBlochParameters((*paramBB)(0),(*paramBB)(1),(*paramBB)(2),(*paramBB)(3),(*paramBB)(4));
        AliInfo(Form("Setting BB parameters from OCDB (AliTPCParam): %.2g, %.2g, %.2g, %.2g, %.2g",
                     (*paramBB)(0),(*paramBB)(1),(*paramBB)(2),(*paramBB)(3),(*paramBB)(4)));
      } else {
        AliError("Couldn't get BB parameters from OCDB, the old default values will be used instead");
      }
    } else {
      AliError("Couldn't get TPC parameters");
    }
  }
  
/*
  AliTPCcalibDB * calib = AliTPCcalibDB::Instance();
  
  //Get pid splines array
  TObjArray *arrSplines=calib->GetPidResponse();
  if (!arrSplines) return;
  AliTPCPIDResponse &tpcPID=esdPID->GetTPCResponse();
  tpcPID.SetUseDatabase(kTRUE);

  // check if parametrisations are already set.
  // since this is uniq for one run, we don't have to reload them
  if (tpcPID.GetResponseFunction(AliPID::kPion)) return;

  // get the default object
  TObject *defaultPID=arrSplines->At(AliPID::kUnknown);
  
  // loop over all particle species and set the response functions
  for (Int_t ispec=0; ispec<AliPID::kUnknown; ++ispec){
    TObject *pidSpline=arrSplines->At(ispec);
    if (!pidSpline) pidSpline=defaultPID;
    tpcPID.SetResponseFunction((AliPID::EParticleType)ispec,pidSpline);
  }
 */ 
}

//_____________________________________________________________________________
void AliTPCReconstructor::ParseOptions( AliTPCtracker* tracker ) const
{
// parse options from rec.C and set in clusterer and tracker
  
  TString option = GetOption();
  
  Int_t useHLTClusters = 3;

  if (option.Contains("use")) {
    
    AliInfo(Form("Overide TPC RecoParam with option %s",option.Data()));
    
    if (option.Contains("useRAW")) {
      useHLTClusters = 1;
      if (option.Contains("useRAWorHLT"))
	useHLTClusters = 2;
    }
    else if (option.Contains("useHLT")) {
      useHLTClusters = 3;
      if (option.Contains("useHLTorRAW"))
	useHLTClusters = 4;
    }
  }
  else {
    const AliTPCRecoParam* param = GetRecoParam();
    useHLTClusters = param->GetUseHLTClusters();
  }

  AliInfo(Form("Usage of HLT clusters in TPC reconstruction : %d", useHLTClusters));

  fClusterer->SetUseHLTClusters(useHLTClusters);
  tracker->SetUseHLTClusters(useHLTClusters);

  return;
}

//_____________________________________________________________________________
void AliTPCReconstructor::SetSystematicErrorCluster( TVectorD *vec ) 
{ 
  // set clusters syst.errors which will override persistent data member from AliTPCRecoParam::fSystematicErrors
  fSystematicErrorClusters=vec;
  AliTPCRecoParam::SetSystematicErrorClusterCustom(vec);
}

//_____________________________________________________________________________
void AliTPCReconstructor::SetPrimaryDCACut( TVectorD *dcacut ) 
{ 
  // set clusters syst.errors which will override persistent data member from AliTPCRecoParam::fSystematicErrors
  fgPrimaryDCACut=dcacut;
  AliTPCRecoParam::SetPrimaryDCACut(dcacut);
}
