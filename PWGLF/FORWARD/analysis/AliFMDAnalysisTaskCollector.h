#ifndef ALIFMDANALYSISTASKCOLLECTOR_H
#define ALIFMDANALYSISTASKCOLLECTOR_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTaskSE.h"
#include "TH1F.h"
#include "TObjArray.h"

class AliESDEvent;
class TChain;
class AliAODEvent;
class AliFMDAnaParameters;
//
//Class to fit energy distributions in the FMD
//
class AliFMDAnalysisTaskCollector : public AliAnalysisTaskSE
{
 public:
  AliFMDAnalysisTaskCollector();
  AliFMDAnalysisTaskCollector(const char* name);
  AliFMDAnalysisTaskCollector(const AliFMDAnalysisTaskCollector& o) : 
    AliAnalysisTaskSE(),
    //    fDebug(o.fDebug),
    fOutputList(o.fOutputList),
    fArray(o.fArray),
    fZvtxDist(o.fZvtxDist),
    fEvents(0),
    fEmptyEvents(0),
    fClusters(0),
    fClustersEmpty(0),
    fFirstEvent(kTRUE), 
    fParam(0)
  {}
  
  AliFMDAnalysisTaskCollector& operator=(const AliFMDAnalysisTaskCollector&) { return *this; }
  virtual ~AliFMDAnalysisTaskCollector() {;}
  // Implementation of interface methods
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void SetDebugLevel(Int_t level) {fDebug = level;}
  virtual void Terminate(Option_t */*option*/);
  void ReadFromFile(const Char_t* filename, Bool_t store=kFALSE, Int_t speciesOption = 0);
  static Double_t  TripleLandau(const Double_t *x, Double_t *par);
  TF1* FitEnergyDistribution(TH1F* hEnergy, Int_t speciesOption);

private:
  void          GetVertex(Double_t* vertexXYZ); 
  //Int_t         fDebug;        //  Debug flag
  TList*        fOutputList;     //Output list
  TObjArray*    fArray;          //Array for storage
  TH1F*         fZvtxDist;       //Dist of z vertex
  Int_t         fEvents;         //Number of events
  Int_t         fEmptyEvents;    //Number of events with empty trigger
  Float_t       fClusters;       //Number of clusters
  Float_t       fClustersEmpty;  //Number of clusters in empty events
  Bool_t        fFirstEvent;     //Have we had events yet ?
  AliFMDAnaParameters* fParam;   //The parameters class for IO

  ClassDef(AliFMDAnalysisTaskCollector, 1); // Analysis task for FMD analysis
};
 
#endif
// Local Variables:
//   mode: C++ 
// End:
// EOF
