/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerInclusiveJet
// \brief helper class to handle a tree for D0 cut optimisation and MVA analyses
// \authors:
// N. Zardoshti, nima.zardoshti@cern.ch
/////////////////////////////////////////////////////////////

#include <TString.h>
#include "AliHFTreeHandlerInclusiveJet.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandlerInclusiveJet);
/// \endcond

//________________________________________________________________
AliHFTreeHandlerInclusiveJet::AliHFTreeHandlerInclusiveJet():
  AliHFTreeHandler()
{
  //
  // Default constructor
  //

}

//________________________________________________________________
AliHFTreeHandlerInclusiveJet::~AliHFTreeHandlerInclusiveJet()
{
  //
  // Default Destructor
  //
}

//________________________________________________________________
TTree* AliHFTreeHandlerInclusiveJet::BuildTree(TString name, TString title) 
{
  fIsMCGenTree=false;
  
  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=nullptr;
  }
  fTreeVar = new TTree(name.Data(),title.Data());

  //set common variables
  AddCommonDmesonVarBranches();

  //set single-track variables
  AddSingleTrackBranches();
  if (fFillJets) AddJetBranches();
  
  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerInclusiveJet::SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, AliAODRecoDecayHF* cand, float bfield, int masshypo, AliPIDResponse *pidrespo) 
{
  fRunNumber=runnumber;
  fEvID=eventID;
  fEvIDExt=eventID_Ext;
  fEvIDLong=eventID_Long;
  fIsMCGenTree=false;

  return true;
}

//________________________________________________________________
bool AliHFTreeHandlerInclusiveJet::SetMCGenVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long) {

  fRunNumber = runnumber;
  fEvID = eventID;
  fEvIDExt = eventID_Ext;
  fEvIDLong = eventID_Long;
  
  return true;
}
