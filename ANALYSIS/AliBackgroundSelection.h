#ifndef ALIBACKGROUNDSELECTION_H
#define ALIBACKGROUNDSELECTION_H
 

// Select events which are not flagged as backgroud
// Author Michele Floris
// michele.floris@cern.ch

#include "AliAnalysisCuts.h"
#include "AliLog.h"

class TList;
class TH2F;
class TH1F;
class TCollection;


class AliBackgroundSelection : public AliAnalysisCuts
{
public:
  // Inherited methods
  AliBackgroundSelection();
  AliBackgroundSelection(const char* name, const char* title);
  AliBackgroundSelection(const AliBackgroundSelection& obj);  
  virtual ~AliBackgroundSelection();
  virtual Bool_t IsSelected(TObject* const obj);
  virtual Bool_t IsSelected(TList*  ) {AliFatal("Not implemented");return 0;}
  virtual void   Init();
  virtual TList * GetOutput() {return fOutputHist;}

  // Helper Methods
  void BookClusterVsTrackletsHisto(const char * trigger_name);
  void BookDeltaPhiHisto(const char * trigger_name);

  TH2F * GetClusterVsTrackletsHisto(const char * trigger_name);
  TH2F * GetClusterVsTrackletsHistoAccepted(const char * trigger_name);
  const char *  GetClusterVsTrackletsHistoName(const char * trigger_name);
  const char *  GetClusterVsTrackletsHistoNameAccepted(const char * trigger_name);

  TH1F * GetDeltaPhiHisto(const char * trigger_name);
  TH1F * GetDeltaPhiHistoAccepted(const char * trigger_name);
  const char *  GetDeltaPhiHistoName(const char * trigger_name);
  const char *  GetDeltaPhiHistoNameAccepted(const char * trigger_name);

  Long64_t Merge(TCollection* const list);
  // Cuts Setters & Getters
  void SetCutParameters(Float_t a,Float_t b) {fACut = a; fBCut =b;}
  void SetDeltaPhiCut(Float_t cut) { fDeltaPhiCut = cut;}
  const Float_t GetCutParameterA() const {return fACut;}
  const Float_t GetCutParameterB() const {return fBCut;}
  const Float_t GetDeltaPhiCut()   const {return fDeltaPhiCut;}

  // TODO: implement cut on global vertex DCA?

private:
  TList * fOutputHist; // contains 2 histo Cluster vs Tracklets and delta phiper trigger type (all and accepted)
  Float_t fACut; // Cut on y = ax + b in the Cluster Vs Tracklets correlation. This is the "a" parameter of the cut
  Float_t fBCut; // Cut on y = ax + b in the Cluster Vs Tracklets correlation. This is the "b" parameter of the cut
  Float_t fDeltaPhiCut; // events with vertex from vertexer Z and DeltaPhi>fDeltaPhiCut are rejected

  AliBackgroundSelection& operator=(const AliBackgroundSelection&);

  ClassDef(AliBackgroundSelection, 1); 
};
 
#endif
