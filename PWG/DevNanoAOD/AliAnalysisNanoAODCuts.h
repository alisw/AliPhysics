#ifndef _ALIANALYSISNANOAODCUTSANDSETTERS_H_
#define _ALIANALYSISNANOAODCUTSANDSETTERS_H_

#include "AliAnalysisCuts.h"
#include "AliNanoAODCustomSetter.h"
#include "AliAnalysisUtils.h"
#include <map>

class AliEventCuts;

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

class AliAnalysisNanoAODV0Cuts : public AliAnalysisCuts
{
public:
  AliAnalysisNanoAODV0Cuts() {}
  virtual ~AliAnalysisNanoAODV0Cuts() {}
  virtual Bool_t IsSelected(TObject* obj); // TObject should be an AliAODv0
  virtual Bool_t IsSelected(TList*   /* list */ ) { return kTRUE; }

private:

  ClassDef(AliAnalysisNanoAODV0Cuts, 1); // track cut object for nano AOD filtering
};

class AliAnalysisNanoAODEventCuts : public AliAnalysisCuts
{
public:
  AliAnalysisNanoAODEventCuts();
  virtual ~AliAnalysisNanoAODEventCuts() {}
  virtual Bool_t IsSelected(TObject* obj); // TObject should be an AliAODEvent
  virtual Bool_t IsSelected(TList*   /* list */ ) { return kTRUE; }
  
  Float_t GetVertexRange() { return fVertexRange; }
  void SetVertexRange (Float_t var) { fVertexRange = var;}
  
  void SetMultiplicityRange(AliAnalysisCuts* cutObject, Int_t minMultiplicity, Int_t maxMultiplicity) { fTrackCut = cutObject; fMinMultiplicity = minMultiplicity; fMaxMultiplicity = maxMultiplicity; }
  
  AliAnalysisUtils* GetAnalysisUtils() { return fAnalysisUtils; }
  void SetCutPileUpMV(Bool_t flag) { fCutPileUpMV = flag; }
  
  void SetAliEventCuts(AliEventCuts* cuts) { fEventCuts = cuts; }
  
private:
  Float_t fVertexRange; // Only events with primary vertex within this range are accepted (whatever the vertex)
  
  AliAnalysisCuts* fTrackCut; // track cut object for multiplicity cut
  Int_t fMinMultiplicity;   // minimum number of tracks to accept this event
  Int_t fMaxMultiplicity;   // maximal number of tracks to accept this event
  
  AliAnalysisUtils* fAnalysisUtils; // AnalysisUtils object
  Bool_t fCutPileUpMV;      // Use fAnalysisUtils->IsPileUpMV to remove pile up. Customize by using GetAnalysisUtils()
  AliEventCuts* fEventCuts; // AliEventCut object for Run 2
  
  ClassDef(AliAnalysisNanoAODEventCuts, 3); // event cut object for nano AOD filtering
};

class AliNanoAODSimpleSetter : public AliNanoAODCustomSetter
{
public:
  AliNanoAODSimpleSetter() : fInitialized(kFALSE) {;}
  virtual ~AliNanoAODSimpleSetter(){;}

  virtual void SetNanoAODHeader(const AliAODEvent * event   , AliNanoAODHeader * head ,TString varListHeader  );
  virtual void SetNanoAODTrack (const AliAODTrack * /*aodTrack*/, AliNanoAODTrack * /*spTrack*/);
  
protected:
  void Init(AliNanoAODHeader* head, TString varListHeader);
  
  Bool_t fInitialized;
  std::map<TString,int> fMultMap;

  ClassDef(AliNanoAODSimpleSetter, 2)

};

#endif /* _ALIANALYSISNANOAODCUTSANDSETTERS_H_ */
