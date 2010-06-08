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
/* $Id: $ */

//__________________________________________
// Class that uses info from the AliJetCorrelSelector object to set up the
// two-particle correlations to be run in one instance of the analysis module.
//-- Author: Paul Constantin

#include "AliJetCorrelMaker.h"

using namespace std;

ClassImp(AliJetCorrelMaker)

AliJetCorrelMaker::AliJetCorrelMaker() : 
  fNumCorrel(0), fNumTrigg(0), fNumAssoc(0), 
  fCorrelType(NULL), fCorrelStr(NULL),
  fTriggType(NULL), fAssocType(NULL),
  fIdxTrigg(NULL), fIdxAssoc(NULL){
  // (default) constructor
}

AliJetCorrelMaker::~AliJetCorrelMaker(){
  // destructor
  fNumCorrel = 0;
  fNumTrigg  = 0;
  fNumAssoc  = 0;
  if(fCorrelType)  delete [] fCorrelType;
  if(fCorrelStr)   delete [] fCorrelStr;
  if(fTriggType)   delete [] fTriggType;
  if(fAssocType)   delete [] fAssocType;
  if(fIdxTrigg)    delete [] fIdxTrigg;
  if(fIdxAssoc)    delete [] fIdxAssoc;
}

Bool_t AliJetCorrelMaker::Init(UInt_t s, UInt_t * const v){
  // Main method. Returns false on initialisation error
  fNumCorrel = s;
  fCorrelType = new UInt_t[fNumCorrel];
  fCorrelStr  = new TString[fNumCorrel];
  fTriggType  = new cPartType_t[fNumCorrel];
  fAssocType  = new cPartType_t[fNumCorrel];
  for(UInt_t k=0; k<fNumCorrel; k++){
    fCorrelType[k] = v[k];
    switch(fCorrelType[k]){
    case 0:
      fTriggType[k] = t_hadron;     fAssocType[k] = t_hadron; fCorrelStr[k] = "DiHadron";
      break;
    case 1:
      fTriggType[k] = t_diphoton;   fAssocType[k] = t_hadron; fCorrelStr[k] = "Pi0Hadron";
      break;
    case 2:
      fTriggType[k] = t_photon;     fAssocType[k] = t_hadron; fCorrelStr[k] = "PhotHadron";
      break;
    case 3:
      fTriggType[k] = t_dielectron; fAssocType[k] = t_hadron; fCorrelStr[k] = "Z0Hadron";
      break;
    case 4:
      fTriggType[k] = t_jet;        fAssocType[k] = t_jet;    fCorrelStr[k] = "JetJet";
      break;
    case 5:
      fTriggType[k] = t_photon;     fAssocType[k] = t_jet;    fCorrelStr[k] = "PhotJet";
      break;
    case 6:
      fTriggType[k] = t_dielectron; fAssocType[k] = t_jet;    fCorrelStr[k] = "Z0Jet";
      break;
    default:
      std::cerr<<"AliJetCorrelMaker::Initialize - ERROR: unknown correlation type!"<<std::endl;
      return kFALSE;
      break;
    }
  } // loop over correlations

  Bool_t notStored;
  fIdxTrigg  = new UInt_t[fNumCorrel];
  cPartType_t *fTriggUniq = new cPartType_t[fNumCorrel];
  for(UInt_t k=0; k<fNumCorrel; k++){
    notStored = kTRUE;
    for(UInt_t i=0; i<fNumTrigg; i++) 
      if(fTriggUniq[i]==fTriggType[k]) notStored = kFALSE;
    if(notStored){fTriggUniq[fNumTrigg]=fTriggType[k]; fNumTrigg++;}
  }
  for(UInt_t k=0; k<fNumCorrel; k++)
    for(UInt_t i=0; i<fNumTrigg; i++)
      if(fTriggType[k]==fTriggUniq[i]) fIdxTrigg[k] = i;
  delete [] fTriggUniq;

  fIdxAssoc  = new UInt_t[fNumCorrel];
  cPartType_t *fAssocUniq = new cPartType_t[fNumCorrel];
  for(UInt_t k=0; k<fNumCorrel; k++){
    notStored = kTRUE;
    for(UInt_t i=0; i<fNumAssoc; i++) 
      if(fAssocUniq[i]==fAssocType[k]) notStored = kFALSE;
    if(notStored){fAssocUniq[fNumAssoc]=fAssocType[k]; fNumAssoc++;}
  }
  for(UInt_t k=0; k<fNumCorrel; k++)
    for(UInt_t i=0; i<fNumAssoc; i++)
      if(fAssocType[k]==fAssocUniq[i]) fIdxAssoc[k] = i;
  delete [] fAssocUniq;

  if(!Check()){
    std::cerr<<"AliJetCorrelMaker::Initialize - array sanity check failed!"<<std::endl;
    return kFALSE;
  }

  return kTRUE;
}

Bool_t AliJetCorrelMaker::Check() const {
  // performs initialization sanity checks
  if(fNumTrigg<1 || fNumAssoc<1) return kFALSE;
  if(fNumTrigg>fNumCorrel || fNumAssoc>fNumCorrel) return kFALSE;
  for(UInt_t k=0; k<fNumCorrel; k++){
    if(fIdxTrigg[k]>=fNumTrigg) return kFALSE;
    if(fIdxAssoc[k]>=fNumAssoc) return kFALSE;
  }
  return kTRUE;
}

TString AliJetCorrelMaker::Descriptor(UInt_t k) const {
  if(k>=fNumCorrel)
    {std::cerr<<"AliJetCorrelMaker::Descriptor overflow!"<<std::endl; return "?";}
  return fCorrelStr[k];
}

UInt_t AliJetCorrelMaker::IdxTrigg(UInt_t k) const {
  if(k>=fNumCorrel)
    {std::cerr<<"AliJetCorrelMaker::IdxTrigg overflow!"<<std::endl; exit(-1);}
  return fIdxTrigg[k];
}

UInt_t AliJetCorrelMaker::IdxAssoc(UInt_t k) const {
  if(k>=fNumCorrel)
    {std::cerr<<"AliJetCorrelMaker::IdxAssoc overflow!"<<std::endl; exit(-1);}
  return fIdxAssoc[k];
}

cPartType_t AliJetCorrelMaker::TriggType(UInt_t k) const {
  if(k>=fNumCorrel)
    {std::cerr<<"AliJetCorrelMaker::TriggType overflow!"<<std::endl; return t_unknown;}
  return fTriggType[k];
}

cPartType_t AliJetCorrelMaker::AssocType(UInt_t k) const {
  if(k>=fNumCorrel)
    {std::cerr<<"AliJetCorrelMaker::AssocType overflow!"<<std::endl; return t_unknown;}
  return fAssocType[k];
}

Bool_t AliJetCorrelMaker::RecoTrigger(UInt_t k) const {
  if(fTriggType[k]==t_diphoton || fTriggType[k]==t_dielectron)
    return kTRUE;
  return kFALSE;
}

Bool_t AliJetCorrelMaker::RecoTrigger() const {
  for(UInt_t k=0; k<fNumCorrel; k++)
    if(RecoTrigger(k)) return kTRUE;
  return kFALSE;
}

void AliJetCorrelMaker::Show() const {
  // print out whole correlation setup
  std::cout<<"Number of Correlations:"<<fNumCorrel
	   <<" Triggers:"<<fNumTrigg<<" Associated:"<<fNumAssoc<<std::endl;
  for(UInt_t k=0; k<fNumCorrel; k++)
    std::cout<<"Correlation("<<k<<"):"<<fCorrelStr[k]
	     <<" TriggType="<<fTriggType[k]<<" AssocType="<<fAssocType[k]
	     <<" IdxTrigg="<<fIdxTrigg[k]<<" IdxAssoc="<<fIdxAssoc[k]<<std::endl;
}
