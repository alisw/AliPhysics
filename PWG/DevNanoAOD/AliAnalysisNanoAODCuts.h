#ifndef _ALIANALYSISNANOAODCUTSANDSETTERS_H_
#define _ALIANALYSISNANOAODCUTSANDSETTERS_H_

#include "AliAnalysisCuts.h"
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

private:
  UInt_t fBitMask; // Only AOD tracks matching this bit mask are accepted

  ClassDef(AliAnalysisNanoAODTrackCuts,1); // Select muon spectrometer tracks
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
public:
  Float_t fVertexRange; // Only events with primary vertex within this range are accepted (whathever the vertex)

  
  ClassDef(AliAnalysisNanoAODEventCuts,1); // Select primary vertices
};





#endif /* _ALIANALYSISNANOAODCUTSANDSETTERS_H_ */
