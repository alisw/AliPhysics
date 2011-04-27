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
  return new AliTPCtrackerMI(param);
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
