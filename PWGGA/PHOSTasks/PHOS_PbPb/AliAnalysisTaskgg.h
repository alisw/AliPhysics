#ifndef AliAnalysisTaskgg_cxx
#define AliAnalysisTaskgg_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class THashList ;
class AliPHOSGeometry;
class AliCaloPhoton ;
class AliAODTrack ;
class AliEPFlattener ;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskgg : public AliAnalysisTaskSE {
public:
  
  enum CutList {kDefault, kDisp, kCPV, kBoth, kDistance1, kDistance2, kDistance3} ;  
  
  
  AliAnalysisTaskgg(const char *name = "AliAnalysisTaskgg");
  virtual ~AliAnalysisTaskgg() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
 
  void SelectPbPb(Bool_t isPbPb=kFALSE){fIsPbPb=isPbPb;} 
  
private:
  AliAnalysisTaskgg(const AliAnalysisTaskgg&); // not implemented
  AliAnalysisTaskgg& operator=(const AliAnalysisTaskgg&); // not implemented

  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z, Double_t w) const ; //Fill 3D histogram witn name key

  
  Int_t ConvertRunNumber(Int_t run) ; 
  Bool_t PairCut(const AliCaloPhoton * ph1, const AliCaloPhoton * ph2, Int_t cut) const ; 
  Bool_t SecondaryPi0Cut(const AliCaloPhoton * ph1, const AliCaloPhoton * ph2) const ;

private:
//  AliStack * fStack ;
  THashList *   fOutputContainer;        //final histogram container
  AliAODEvent * fEvent ;        //!
//  TClonesArray * fStack ;  
  TList *       fPHOSEvents[10][10][11] ; //Containers for events with PHOS photons
  TClonesArray* fPHOSEvent ;      //PHOS photons in current event
 
  //Reaction plain for v2
  AliEPFlattener * fV0AFlat ; //!
  AliEPFlattener * fV0CFlat ; //!
  

  Int_t fRunNumber ;    //Current run number
  Float_t fCentrality ; //!Centrality of the currecnt event
  Int_t fCenBin ;       //! Current centrality bin

  AliPHOSGeometry  *fPHOSGeo;  //! PHOS geometry
  Int_t fEventCounter;         // number of analyzed events
  Bool_t fIsPbPb ;             // Switch to PbPb or pp/pPb mode

  ClassDef(AliAnalysisTaskgg, 1); // PHOS analysis task
};

#endif
