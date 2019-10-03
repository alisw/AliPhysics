
/* $Id$ */

#include "AliEventPoolMuon.h"
#include "AliRunTagCuts.h"
#include "AliLHCTagCuts.h"
#include "AliDetectorTagCuts.h"
#include "AliEventTagCuts.h"
#include "AliTagAnalysis.h"

// Realisation of an AliVEventPool via
// on the flight generation of the bin using AliTagAnalysis.
// Created expanding AliEventPoolOTF functionalities
//
// Authors Alessandro De Falco and Antonio Uras, INFN Cagliari
// alessandro.de.falco@ca.infn.it  antonio.uras@ca.infn.it

#define AliEventPoolMuon_CXX

ClassImp(AliEventPoolMuon)

//=====================================================================================================

AliEventPoolMuon::AliEventPoolMuon():
  AliVEventPool(),
  fTagAnalysis(0),
  fRunCuts(0),
  fLHCCuts(0),
  fDetectorCuts(0),
  fEventCuts(0),
  fTagDirectory(0),
  fMultiplicityMin(0),
  fMultiplicityMax(0),
  fMultiplicityStep(0),
  fMultiplicity(0),
  fNFWMuonMin(0),
  fNFWMuonMax(0),
  fNFWMuonStep(0),
  fNFWMuon(0),
  fPrimaryVertexZMin(0),
  fPrimaryVertexZMax(0),
  fPrimaryVertexZStep(0),
  fPrimaryVertexZ(0),
  fBinNumber(0) {
  
  // Default constructor

}

//=====================================================================================================

AliEventPoolMuon::AliEventPoolMuon(const Char_t *name, const Char_t *title):
  AliVEventPool(name, title),
  fTagAnalysis(new AliTagAnalysis(title)),
  fRunCuts(new AliRunTagCuts()),
  fLHCCuts(new AliLHCTagCuts()),
  fDetectorCuts(new AliDetectorTagCuts()),
  fEventCuts(new AliEventTagCuts()),
  fTagDirectory("."),
  fMultiplicityMin(0),
  fMultiplicityMax(0),
  fMultiplicityStep(0),
  fMultiplicity(0),
  fNFWMuonMin(0),
  fNFWMuonMax(0),
  fNFWMuonStep(0),
  fNFWMuon(0),
  fPrimaryVertexZMin(0),
  fPrimaryVertexZMax(0),
  fPrimaryVertexZStep(0),
  fPrimaryVertexZ(0),
  fBinNumber(0) {

  // Constructor

}

//=====================================================================================================

AliEventPoolMuon::AliEventPoolMuon(const AliEventPoolMuon& obj):
  AliVEventPool(obj),
  fTagAnalysis(0),
  fRunCuts(0),
  fLHCCuts(0),
  fDetectorCuts(0),
  fEventCuts(0),
  fTagDirectory(0),
  fMultiplicityMin(0),
  fMultiplicityMax(0),
  fMultiplicityStep(0),
  fMultiplicity(0),
  fNFWMuonMin(0),
  fNFWMuonMax(0),
  fNFWMuonStep(0),
  fNFWMuon(0),
  fPrimaryVertexZMin(0),
  fPrimaryVertexZMax(0),
  fPrimaryVertexZStep(0),
  fPrimaryVertexZ(0),
  fBinNumber(0) {

  // Copy constructor

}

//=====================================================================================================

AliEventPoolMuon& AliEventPoolMuon::operator=(const AliEventPoolMuon& other) {

  // Assignment operator
  AliVEventPool::operator=(other);
  return *this;

}

//=====================================================================================================

void AliEventPoolMuon::Init() {
  
  fTagAnalysis -> ChainLocalTags(fTagDirectory);

  fMultiplicity   = fMultiplicityMin;
  fNFWMuon        = fNFWMuonMin;
  fPrimaryVertexZ = fPrimaryVertexZMin;

}

//=====================================================================================================

TChain* AliEventPoolMuon::GetNextChain() {
  
  TChain *chain = 0;
  fBinNumber++;

  // hierarchic order of variables: multiplicity -> nFWMuons -> primaryVertexZ
  
  Double_t primaryVertexZMax_TMP = fPrimaryVertexZ + fPrimaryVertexZStep;
  if (primaryVertexZMax_TMP > fPrimaryVertexZMax) {
    fPrimaryVertexZ = fPrimaryVertexZMin;
    primaryVertexZMax_TMP = fPrimaryVertexZ + fPrimaryVertexZStep;
    fNFWMuon += fNFWMuonStep;
  }

  Int_t nFWMuonMax_TMP = fNFWMuon + fNFWMuonStep - 1;
  if (nFWMuonMax_TMP > fNFWMuonMax) {
    fNFWMuon = fNFWMuonMin;
    nFWMuonMax_TMP = fNFWMuon + fNFWMuonStep - 1;
    fMultiplicity += fMultiplicityStep;
  }

  Int_t multiplicityMax_TMP = fMultiplicity + fMultiplicityStep - 1;
  if (multiplicityMax_TMP > fMultiplicityMax) return 0;
  else {
    printf("\n");
    printf("mixing events in pool #%02d:  multiplicity    %d -> %d\n",fBinNumber,fMultiplicity,multiplicityMax_TMP);
    printf("                            nFWMuons        %d -> %d\n",fNFWMuon,nFWMuonMax_TMP);
    printf("                            vertexZ         %f -> %f\n\n",fPrimaryVertexZ, primaryVertexZMax_TMP);
    fEventCuts->SetPrimaryVertexZRange(fPrimaryVertexZ, primaryVertexZMax_TMP);
    fEventCuts->SetNFWMuonRange(fNFWMuon, nFWMuonMax_TMP);
    fEventCuts->SetMultiplicityRange(fMultiplicity, multiplicityMax_TMP);
    chain = fTagAnalysis->QueryTags(fRunCuts, fLHCCuts, fDetectorCuts, fEventCuts);
    fPrimaryVertexZ += fPrimaryVertexZStep;     // here the innermost bin-variable has to be increased
    return chain;
  }
  
}

//=====================================================================================================

void  AliEventPoolMuon::GetCurrentBin(Float_t* /*bin*/) {

  //

}

//=====================================================================================================

Int_t AliEventPoolMuon::GetDimension() {

  return (1);

}

//=====================================================================================================
