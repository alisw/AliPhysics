#ifndef ALIANALYSISUTILS_H
#define ALIANALYSISUTILS_H

////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with functions useful for different analyses                      //
// - vertex selection                                                       //
//    * 2013 pA default cuts                                             //
// - identification of the fist event of the chunk                         //
// - identification pileup events                                           //
//                                                                          //
///////////////////////////////////////////////////////////////////

#include "TObject.h"

class AliVEvent;
class AliVVertex;

class AliAnalysisUtils : public TObject {

 public:

  AliAnalysisUtils();
  virtual ~AliAnalysisUtils(){};
  
  Bool_t IsVertexSelected2013pA(AliVEvent *event);
  Bool_t IsFirstEventInChunk(AliVEvent *event);
  
  Bool_t IsPileUpEvent(AliVEvent *event); //to be used in the analysis
  Bool_t IsPileUpMV(AliVEvent *event); //MV pileup selection implemented here
  Bool_t IsPileUpSPD(AliVEvent *event); //this calls IsPileUpFromSPD
  Bool_t IsOutOfBunchPileUp(AliVEvent *event); //out-of-bunch pileup rejection using trigger information
  
  Double_t GetWDist(const AliVVertex* v0, const AliVVertex* v1);
  
  void SetMinVtxContr(Int_t contr=1) {fMinVtxContr=contr;}
  void SetMaxVtxZ(Float_t z=1e6) {fMaxVtxZ=z;}
  void SetCutOnZVertexSPD(Bool_t iscut=true) { fCutOnZVertexSPD = iscut; }
  
  //general pileup selection settings
  void SetUseMVPlpSelection(Bool_t useMVPlpSelection) { fUseMVPlpSelection = useMVPlpSelection;}
  void SetUseOutOfBunchPileUp(Bool_t useOutOfBunchPileUp) { fUseOutOfBunchPileUp = useOutOfBunchPileUp;}
  
  //Multi Vertex pileup selection
  void SetMinPlpContribMV(Int_t minPlpContribMV) { fMinPlpContribMV = minPlpContribMV;}
  void SetMaxPlpChi2MV(Float_t maxPlpChi2MV) { fMaxPlpChi2MV = maxPlpChi2MV;}
  void SetMinWDistMV(Float_t minWDistMV) { fMinWDistMV = minWDistMV;}
  void SetCheckPlpFromDifferentBCMV(Bool_t checkPlpFromDifferentBCMV) { fCheckPlpFromDifferentBCMV = checkPlpFromDifferentBCMV;}
  //SPD Pileup slection
  void SetMinPlpContribSPD(Int_t minPlpContribSPD) { fMinPlpContribSPD = minPlpContribSPD;}
  void SetMinPlpZdistSPD(Float_t minPlpZdistSPD) { fMinPlpZdistSPD = minPlpZdistSPD;}
  void SetnSigmaPlpZdistSPD(Float_t nSigmaPlpZdistSPD) { fnSigmaPlpZdistSPD = nSigmaPlpZdistSPD;}
  void SetnSigmaPlpDiamXYSPD(Float_t nSigmaPlpDiamXYSPD) { fnSigmaPlpDiamXYSPD = nSigmaPlpDiamXYSPD;}
  void SetnSigmaPlpDiamZSPD(Float_t nSigmaPlpDiamZSPD) { fnSigmaPlpDiamZSPD = nSigmaPlpDiamZSPD;}
  void SetUseSPDCutInMultBins(Bool_t useSPDCutInMultBins) { fUseSPDCutInMultBins = useSPDCutInMultBins;}
  
  
 private:
  
  Bool_t fisAOD; // flag for AOD:1 or ESD:0
  
  Int_t fMinVtxContr; // minimum vertex contributors
  Float_t fMaxVtxZ;   // maximum |z| of primary vertex
  
  Bool_t fCutOnZVertexSPD; // 0: no cut, 1: |zvtx-SPD - zvtx-TPC|<0.5cm
  
  Bool_t  fUseMVPlpSelection; //switch between SPD and MV pileup
  Bool_t  fUseOutOfBunchPileUp; //BC information from the trigger is used to tag out-of-bunch pileup
  
  Int_t    fMinPlpContribMV; //minimum contributors to the pilup vertices, multi-vertex
  Float_t  fMaxPlpChi2MV; //minimum value of Chi2perNDF of the pileup vertex, multi-vertex
  Float_t  fMinWDistMV; //minimum of the sqrt of weighted distance between the primary and the pilup vertex, multi-vertex
  Bool_t  fCheckPlpFromDifferentBCMV; //pileup from different BC (Bunch Crossings)
  
  Int_t    fMinPlpContribSPD; //minimum contributors to the pilup vertices, SPD
  Float_t  fMinPlpZdistSPD;  //minimum distance for the SPD pileup vertex
  Float_t  fnSigmaPlpZdistSPD;  //cut on Z for SPD pileup
  Float_t  fnSigmaPlpDiamXYSPD;  //cut on nsigma diamond XY for SPD pileup
  Float_t  fnSigmaPlpDiamZSPD;  //cut on nsigma diamond Z for SPD pileup
  Bool_t  fUseSPDCutInMultBins;  //use IsPileupFromSPDInMultBins instead of IsPileupFromSPD
  
  ClassDef(AliAnalysisUtils,1) // base helper class
};
#endif
 
