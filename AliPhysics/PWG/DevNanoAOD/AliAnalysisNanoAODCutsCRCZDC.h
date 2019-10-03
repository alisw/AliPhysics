#ifndef _ALIANALYSISNANOAODCUTSANDSETTERSCRCZDC_H_
#define _ALIANALYSISNANOAODCUTSANDSETTERSCRCZDC_H_

#include "AliAnalysisCuts.h"
#include "AliNanoAODCustomSetter.h"
//#include "AliNanoAODCustomSetter.h"

class AliVVertex;

class AliAnalysisNanoAODTrackCutsCRCZDC : public AliAnalysisCuts
{
public:
  AliAnalysisNanoAODTrackCutsCRCZDC();
  virtual ~AliAnalysisNanoAODTrackCutsCRCZDC()  {}
  virtual Bool_t IsSelected(TObject* obj); // TObject should be an AliAODTrack
  virtual Bool_t IsSelected(TList*   /* list */ ) { return kTRUE; }
  UInt_t GetBitMask() { return fBitMask; }
  void  SetBitMask (UInt_t var) { fBitMask = var;}
  Float_t GetMinPt() { return fMinPt; }
  void  SetMinPt (Float_t var) { fMinPt = var;}
  Float_t GetMaxEta() { return fMaxEta; }
  void  SetMaxEta (Float_t var) { fMaxEta = var;}

private:
  UInt_t fBitMask; // Only AOD tracks matching this bit mask are accepted
  Float_t fMinPt; // miminum pt of the tracks
  Float_t fMaxEta; // MaxEta

  ClassDef(AliAnalysisNanoAODTrackCutsCRCZDC,1); // track cut object for nano AOD filtering
};

class AliAnalysisNanoAODEventCutsCRCZDC  : public AliAnalysisCuts
{
public:
  AliAnalysisNanoAODEventCutsCRCZDC();
  virtual ~AliAnalysisNanoAODEventCutsCRCZDC() {}
  virtual Bool_t IsSelected(TObject* obj); // TObject should be an AliAODEvent
  virtual Bool_t IsSelected(TList*   /* list */ ) { return kTRUE; }
  Float_t GetVertexRange() { return fVertexRange; }
  void  SetVertexRange (Float_t var) { fVertexRange = var;}
  void  SetMultiplicityRange(AliAnalysisCuts* cutObject, Int_t minMultiplicity, Int_t maxMultiplicity) { fTrackCut = cutObject; fMinMultiplicity = minMultiplicity; fMaxMultiplicity = maxMultiplicity; }
  
public:
  Float_t fVertexRange; // Only events with primary vertex within this range are accepted (whatever the vertex)
  
  AliAnalysisCuts* fTrackCut; // track cut object for multiplicity cut
  Int_t fMinMultiplicity;   // minimum number of tracks to accept this event
  Int_t fMaxMultiplicity;   // maximal number of tracks to accept this event
  
  ClassDef(AliAnalysisNanoAODEventCutsCRCZDC,1); // event cut object for nano AOD filtering
};

class AliNanoAODSimpleSetterCRCZDC  : public AliNanoAODCustomSetter
{
public:
  AliNanoAODSimpleSetterCRCZDC(){;}
  virtual ~AliNanoAODSimpleSetterCRCZDC(){;}

  virtual void SetNanoAODHeader(const AliAODEvent * event   , AliNanoAODHeader * head, TString varListHeader   );
  virtual void SetNanoAODTrack (const AliAODTrack * /*aodTrack*/, AliNanoAODTrack * /*spTrack*/){;}
  Bool_t       SelectPileup    (AliAODEvent* aod);
  Bool_t       plpMV           (const AliAODEvent* aod);
  Double_t     GetWDist        (const AliVVertex* v0, const AliVVertex* v1);
  virtual void SetRejectPileUpTight(Bool_t pileuptight) { fRejectPileUpTight = pileuptight; }

private:
  Bool_t fRejectPileUpTight;  

  ClassDef(AliNanoAODSimpleSetterCRCZDC, 1)

};




#endif /* _ALIANALYSISNANOAODCUTSANDSETTERSCRCZDC_H_ */
