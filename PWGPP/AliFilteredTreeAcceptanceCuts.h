#ifndef ALIFILTEREDTREEACCEPTANCECUTS_H
#define ALIFILTEREDTREEACCEPTANCECUTS_H

//------------------------------------------------------------------------------
// Class to keep selection cuts for MC tracks. 
// 
// Author: J.Otwinowski 03/11/2008 
// last change: 2011-04-04 by M.Knichel
// Modifications - usage of DCA cuts  marian.ivanov@cern.ch
//------------------------------------------------------------------------------

class TParticle;
class AliESDtrack;
class AliExternalTrackParam;

#include "AliAnalysisCuts.h"

class AliFilteredTreeAcceptanceCuts : public AliAnalysisCuts
{
public:
  AliFilteredTreeAcceptanceCuts(const Char_t* name ="AliFilteredTreeAcceptanceCuts", const Char_t *title ="");
  virtual ~AliFilteredTreeAcceptanceCuts(); 
 
  // setters 
  void SetEtaRange(Float_t min=-1e99, Float_t max=1e99)  { fMinEta=min; fMaxEta=max; }
  void SetPhiRange(Float_t min=-1e99, Float_t max=1e99)  { fMinPhi=min; fMaxPhi=max;}
  void SetPtRange(Float_t min=-1e99, Float_t max=1e99)   { fMinPt=min;  fMaxPt=max;}
  void SetExcludeEtaPhiRange(Float_t etaMin, Float_t etaMax, Float_t phiMin, Float_t phiMax)
  	{ fExcludeMinEta = etaMin; fExcludeMaxEta = etaMax; fExcludeMinPhi = phiMin; fExcludeMaxPhi = phiMax; fCheckRange=kTRUE; }
  void SetExcludeEtaPhiRange2(Float_t etaMin, Float_t etaMax, Float_t phiMin, Float_t phiMax)
  	{ fExcludeMinEta2 = etaMin; fExcludeMaxEta2 = etaMax; fExcludeMinPhi2 = phiMin; fExcludeMaxPhi2 = phiMax; fCheckRange=kTRUE; }  	

  void SetMaxDCAr(Float_t max=100) { fMaxDCAr=max;}
  void SetMaxDCAz(Float_t max=100) { fMaxDCAz=max;}

  // getters 
  Float_t GetMinEta() const {return fMinEta;}
  Float_t GetMaxEta() const {return fMaxEta;}
  Float_t GetMinPhi() const {return fMinPhi;}
  Float_t GetMaxPhi() const {return fMaxPhi;}
  Float_t GetMinPt() const {return fMinPt;}
  Float_t GetMaxPt() const {return fMaxPt;}
  
  Bool_t  GetCheckRange() const { return fCheckRange; }
  Float_t GetExcludeMinEta() const { return fExcludeMinEta; }
  Float_t GetExcludeMaxEta() const { return fExcludeMaxEta; }
  Float_t GetExcludeMinPhi() const { return fExcludeMinPhi; }
  Float_t GetExcludeMaxPhi() const { return fExcludeMaxPhi; }  

  Float_t GetMaxDCAr() const {return fMaxDCAr;}
  Float_t GetMaxDCAz() const {return fMaxDCAz;}

  // cuts init function
  void Init();

  // check MC tracks
  virtual Bool_t IsSelected(TObject *) {return kTRUE;}
  virtual Bool_t IsSelected(TList *) {return kTRUE;}

  //
  Bool_t AcceptTrack(AliESDtrack *track);
  Bool_t AcceptTrackLocalTPC(AliESDtrack *track);
  Bool_t AcceptTrack(AliExternalTrackParam *track);
  Bool_t AcceptTrack(TParticle *particle);
  
  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* list);

private:
  Float_t fMinEta; // min pseudorapidity 
  Float_t fMaxEta; // max pseudorapidity
  Float_t fMinPhi; // min azimuthal angle (rad)
  Float_t fMaxPhi; // max azimuthal angle (rad)
  Float_t fMinPt;  // min pt
  Float_t fMaxPt;  // max pt
  
  Float_t fExcludeMinEta;
  Float_t fExcludeMaxEta;
  Float_t fExcludeMinPhi;
  Float_t fExcludeMaxPhi;
  Float_t fExcludeMinEta2;
  Float_t fExcludeMaxEta2;
  Float_t fExcludeMinPhi2;
  Float_t fExcludeMaxPhi2;  
  Bool_t  fCheckRange;

  // max DCAr and DCAz with respect
  // to nominal vertex position
  Float_t fMaxDCAr; // min DCAr
  Float_t fMaxDCAz; // max DCAz
 
  AliFilteredTreeAcceptanceCuts(const AliFilteredTreeAcceptanceCuts&); // not implemented
  AliFilteredTreeAcceptanceCuts& operator=(const AliFilteredTreeAcceptanceCuts&); // not implemented

  ClassDef(AliFilteredTreeAcceptanceCuts, 1)
};

#endif // 
