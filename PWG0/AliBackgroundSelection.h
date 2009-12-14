#ifndef ALIBACKGROUNDSELECTION_H
#define ALIBACKGROUNDSELECTION_H
 

// Select events which are not flagged as backgroud
// Author Michele Floris
// michele.floris@cern.ch

#include <TNamed.h>
#include "AliAnalysisCuts.h"
#include "AliLog.h"

class TList;
class TH2F;
class TCollection;

class AliBackgroundSelection : public AliAnalysisCuts
{
public:
  AliBackgroundSelection();
  AliBackgroundSelection(const char* name, const char* title);
  AliBackgroundSelection(const AliBackgroundSelection& obj);  
  virtual ~AliBackgroundSelection();
  virtual Bool_t IsSelected(TObject* obj);
  virtual Bool_t IsSelected(TList*  ) {AliFatal("Not implemented");return 0;}
  virtual void   Init();
  virtual TList * GetOutput() {return fOutputHist;}

  void BookClusterVsTrackletsHisto(const char * trigger_name);
  TH2F * GetClusterVsTrackletsHisto(const char * trigger_name);
  TH2F * GetClusterVsTrackletsHistoAccepted(const char * trigger_name);
  const char *  GetClusterVsTrackletsHistoName(const char * trigger_name);
  const char *  GetClusterVsTrackletsHistoNameAccepted(const char * trigger_name);
  Long64_t Merge(TCollection* list);
  void SetCutParameters(Float_t a,Float_t b) {fACut = a; fBCut =b;}
  Float_t GetCutParameterA(){return fACut;}
  Float_t GetCutParameterB(){return fBCut;}

  // TODO: implement cut on global vertex DCA?

private:
  TList * fOutputHist; // contains 2 histo Cluster vs Tracklets per trigger type (all and accepted)
  Float_t fACut; // Cut on y = ax + b in the Cluster Vs Tracklets correlation. This is the "a" parameter of the cut
  Float_t fBCut; // Cut on y = ax + b in the Cluster Vs Tracklets correlation. This is the "b" parameter of the cut
    
  AliBackgroundSelection& operator=(const AliBackgroundSelection&);

  ClassDef(AliBackgroundSelection, 1); 
};
 
#endif
