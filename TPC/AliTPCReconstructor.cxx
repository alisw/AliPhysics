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


#include "AliTPCReconstructor.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliRawReader.h"
#include "AliTPCclustererMI.h"
#include "AliTPCtrackerMI.h"
//#include "AliTPCpidESD.h"
#include "AliTPCParam.h"
#include "AliTPCParamSR.h"
#include "AliTPCcalibDB.h"
#include "AliTracker.h"
#include "AliMagF.h"

ClassImp(AliTPCReconstructor)


Int_t    AliTPCReconstructor::fgStreamLevel     = 1;        // stream (debug) level
AliTPCAltroEmulator *  AliTPCReconstructor::fAltroEmulator=0;    // ALTRO emulator

AliTPCReconstructor::AliTPCReconstructor():
AliReconstructor(),
fClusterer(NULL)
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
  fClusterer = new AliTPCclustererMI(param);
}

AliTPCReconstructor::AliTPCReconstructor(const AliTPCReconstructor& /*rec*/):
AliReconstructor(),
fClusterer(NULL)
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
  
  AliTPCtrackerMI* tracker = new AliTPCtrackerMI(param);

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
void AliTPCReconstructor::ParseOptions( AliTPCtrackerMI* tracker ) const
{
// parse options from rec.C and set in clusterer and tracker
  
  TString option = GetOption();
  
  Int_t useHLTClusters = 3;

  if (option.Contains("use")) {
    
    AliInfo(Form("Overide TPC RecoParam with option %s",option.Data()));
    
    if (!option.CompareTo("useRAW"))
      useHLTClusters = 1;
    if (!option.CompareTo("useRAWorHLT"))
      useHLTClusters = 2;
    if (!option.CompareTo("useHLT"))
      useHLTClusters = 3;
    if (!option.CompareTo("useHLTorRAW"))
      useHLTClusters = 4;
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
