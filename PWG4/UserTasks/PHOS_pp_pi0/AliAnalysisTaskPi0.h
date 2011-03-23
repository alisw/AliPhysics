#ifndef AliAnalysisTaskPi0_cxx
#define AliAnalysisTaskPi0_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TObjArray;
class TH1F;
class TH2I;
class TH2F;
class TH3F;
class AliESDtrackCuts;
class AliPHOSGeometry;
class AliTriggerAnalysis;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPi0 : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskPi0(const char *name = "AliAnalysisTaskPi0");
  virtual ~AliAnalysisTaskPi0() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  void SetPHOSBadMap(Int_t mod,TH2I * h)
  {
    if(fPHOSBadMap[mod]) delete fPHOSBadMap[mod] ;
    fPHOSBadMap[mod]=new TH2I(*h) ;
    printf("Set %s \n",fPHOSBadMap[mod]->GetName());
  }
  
private:
  AliAnalysisTaskPi0(const AliAnalysisTaskPi0&); // not implemented
  AliAnalysisTaskPi0& operator=(const AliAnalysisTaskPi0&); // not implemented
  Bool_t IsGoodChannel(const char * det, Int_t mod,Int_t ix, Int_t iz);
  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  Bool_t TestLambda(Double_t l1,Double_t l2) ;

 
private:
  AliESDtrackCuts *fESDtrackCuts; // Track cut
  TList * fOutputContainer;       //final histogram container
  TList * fPHOSEvents[10][2] ;    //Container for PHOS photons
  TClonesArray * fPHOSEvent ;     //PHOS photons in current event
 
  Int_t fnCINT1B;           // Number of CINT1B triggers
  Int_t fnCINT1A;           // Number of CINT1A triggers
  Int_t fnCINT1C;           // Number of CINT1C triggers
  Int_t fnCINT1E;           // Number of CINT1E triggers

  TH2I *fPHOSBadMap[6] ;    //Container for PHOS bad channels map

  AliPHOSGeometry  *fPHOSGeo;  // PHOS geometry
  Int_t fEventCounter;         // number of analyzed events
  AliTriggerAnalysis *fTriggerAnalysis; //! Trigger Analysis for Normalisation

  ClassDef(AliAnalysisTaskPi0, 2); // PHOS analysis task
};

#endif
