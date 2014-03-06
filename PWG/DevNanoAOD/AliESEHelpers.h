#ifndef _ALIESEHELPERS_H_
#define _ALIESEHELPERS_H_

#include "AliAnalysisCuts.h"
#include "AliSpectraAODTrackCuts.h"
#include "AliSpectraAODEventCuts.h"
#include "AliNanoAODCustomSetter.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"


// Wrappers for the spectra classes which inherit from TNamed and not from AliAnalysisCuts


class AliESETrkCut : public AliAnalysisCuts
{
public: 
  AliESETrkCut () : AliAnalysisCuts("trkcuts", "Track cuts"), fTrackCuts(0) {;}
  virtual ~AliESETrkCut(){;}
  AliSpectraAODTrackCuts * GetTrackCuts() { return fTrackCuts; }
  void  SetTrackCuts (AliSpectraAODTrackCuts * var) { fTrackCuts = var; }
  virtual Bool_t IsSelected(TObject* obj) {return fTrackCuts->IsSelected(dynamic_cast<AliAODTrack*>(obj), 1);};
  virtual Bool_t IsSelected(TList*   /* list */ ) { return kTRUE; }

private:
  AliSpectraAODTrackCuts * fTrackCuts; // TrackCuts

  AliESETrkCut(const AliESETrkCut&); // not implemented
  AliESETrkCut& operator=(const AliESETrkCut&); // not implemented
  
  ClassDef(AliESETrkCut, 1)
};

class AliESEEvtCut  : public AliAnalysisCuts
{
public:
  AliESEEvtCut() : AliAnalysisCuts("evtcuts", "Event cuts"), fEventCuts(0), fTrackCuts(0) {;}
  virtual ~AliESEEvtCut(){;}
  AliSpectraAODEventCuts * GetEventCuts() { return fEventCuts; }
  void  SetEventCuts (AliSpectraAODEventCuts * var) { fEventCuts = var;}
  virtual Bool_t IsSelected(TObject* obj) { return fEventCuts->IsSelected(dynamic_cast<AliAODEvent*>(obj), fTrackCuts);};
  virtual Bool_t IsSelected(TList*   /* list */ ) { return kTRUE; }
  virtual void Print(Option_t *) const { fEventCuts->PrintCuts(); }
  AliSpectraAODTrackCuts * GetTrackCuts() { return fTrackCuts; }
  void  SetTrackCuts (AliSpectraAODTrackCuts * var) { fTrackCuts = var;}

private:
  AliSpectraAODEventCuts * fEventCuts; // EventCuts, static memeber for copying
  AliSpectraAODTrackCuts * fTrackCuts;  // pointer to track cuts, needed 

  AliESEEvtCut(const AliESEEvtCut&); // not implemented
  AliESEEvtCut& operator=(const AliESEEvtCut&); // not implemented


  ClassDef(AliESEEvtCut,1)
};


class AliAnalysisESESetter  : public AliNanoAODCustomSetter
{
public:
  AliAnalysisESESetter(const char * name = "AliAnalysisNanoAODFilters") : AliNanoAODCustomSetter(name), fEventCuts(0) {;}
  virtual ~AliAnalysisESESetter() {;}

  virtual void SetNanoAODHeader(const AliAODEvent * event   , AliNanoAODHeader * head  );
  virtual void SetNanoAODTrack (const AliAODTrack * aodTrack, AliNanoAODTrack * spTrack);
  AliSpectraAODEventCuts* GetEventCuts() { return fEventCuts; }
  void  SetEventCuts (AliSpectraAODEventCuts* var) { fEventCuts = var;}

private:  
  AliSpectraAODEventCuts* fEventCuts; // EventCuts

  AliAnalysisESESetter(const AliAnalysisESESetter&); // not implemented
  AliAnalysisESESetter& operator=(const AliAnalysisESESetter&); // not implemented

  ClassDef(AliAnalysisESESetter, 1)
};
  

#endif /* _ALIESEHELPERS_H_ */
