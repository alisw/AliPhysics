#ifndef _ALIANALYSISNANOAODCUTSANDSETTERS_H_
#define _ALIANALYSISNANOAODCUTSANDSETTERS_H_

#include "AliAnalysisCuts.h"
#include "AliNanoAODCustomSetter.h"
#include "AliNanoAODCustomSetter.h"

class AliAnalysisNanoAODTrackCuts : public AliAnalysisCuts
{
public:
  AliAnalysisNanoAODTrackCuts();
  virtual ~AliAnalysisNanoAODTrackCuts()  {}
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



  ClassDef(AliAnalysisNanoAODTrackCuts,1); // track cut object for nano AOD filtering
};

class AliAnalysisNanoAODEventCuts : public AliAnalysisCuts
{
public:
  AliAnalysisNanoAODEventCuts();
  virtual ~AliAnalysisNanoAODEventCuts() {}
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
  
  ClassDef(AliAnalysisNanoAODEventCuts,2); // event cut object for nano AOD filtering
};

class AliNanoAODSimpleSetter : public AliNanoAODCustomSetter
{
public:
  AliNanoAODSimpleSetter(){;}
  virtual ~AliNanoAODSimpleSetter(){;}

  virtual void SetNanoAODHeader(const AliAODEvent * event   , AliNanoAODHeader * head  );
  virtual void SetNanoAODTrack (const AliAODTrack * /*aodTrack*/, AliNanoAODTrack * /*spTrack*/){;}

  ClassDef(AliNanoAODSimpleSetter, 1)

};




#endif /* _ALIANALYSISNANOAODCUTSANDSETTERS_H_ */
