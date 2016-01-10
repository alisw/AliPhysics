#ifndef AliAnalysisTaskgg_cxx
#define AliAnalysisTaskgg_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class THashList ;
class AliPHOSGeometry;
class AliCaloPhoton ;
class AliAODTrack ;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskgg : public AliAnalysisTaskSE {
public:
  
  enum CutList {kDefault, kDisp, kCPV, kBoth, kDistance1, kDistance2, kDistance3} ;  
  
  
  AliAnalysisTaskgg(const char *name = "AliAnalysisTaskgg");
  virtual ~AliAnalysisTaskgg() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
 
private:
  AliAnalysisTaskgg(const AliAnalysisTaskgg&); // not implemented
  AliAnalysisTaskgg& operator=(const AliAnalysisTaskgg&); // not implemented

  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z, Double_t w) const ; //Fill 3D histogram witn name key

  Bool_t GetTPCEventPlane(Double_t &epAngle, Double_t &qsubRes) ;
  Bool_t SelectTrack(AliAODTrack * t) ;
  
  Int_t ConvertRunNumber(Int_t run) ; 
  Bool_t PairCut(const AliCaloPhoton * ph1, const AliCaloPhoton * ph2, Int_t cut) const ; 
  Bool_t SecondaryPi0Cut(const AliCaloPhoton * ph1, const AliCaloPhoton * ph2) const ;

private:
//  AliStack * fStack ;
  THashList *   fOutputContainer;        //final histogram container
  AliAODEvent * fEvent ;
  TClonesArray * fStack ;  
  TList *       fPHOSEvents[1][10][11] ; //Containers for events with PHOS photons
  TClonesArray* fPHOSEvent ;      //PHOS photons in current event
 
  //Reaction plain for v2
  Float_t fRP ; //!Reaction plane calculated with full TPC 
  Float_t fRPV0A ;    //!Reaction plain calculated with A-side TPC: eta>0.15 
  Float_t fRPV0C ;    //!Reaction plain calculated with C-side TPC: eta<-0.15
  Bool_t fHaveTPCRP ; //! Is TPC RP defined?
  
  //V0 calibration

  static const Int_t nCentrBinV0 = 9; // # cenrality bins

  TProfile *fMultV0;                  // object containing VZERO calibration information

  Float_t fV0Cpol,fV0Apol;            // loaded by OADB

  Float_t fMeanQ[nCentrBinV0][2][2];    // and recentering

  Float_t fWidthQ[nCentrBinV0][2][2];   // ...


  Int_t fRunNumber ;    //Current run number
  Float_t fCentrality ; //!Centrality of the currecnt event
  Int_t fCenBin ;       //! Current centrality bin

  AliPHOSGeometry  *fPHOSGeo;  //! PHOS geometry
  Int_t fEventCounter;         // number of analyzed events
  TH1D*  fTPCflatC2 ;   //HIstos with flattening parameters
  TH1D*  fTPCflatS2 ;
  TH1D*  fTPCflatC4 ;
  TH1D*  fTPCflatS4 ;

  TH1D*  fV0AflatC2 ;
  TH1D*  fV0AflatS2 ;
  TH1D*  fV0AflatC4 ;
  TH1D*  fV0AflatS4 ;
    
  TH1D*  fV0CflatC2 ;
  TH1D*  fV0CflatS2 ;
  TH1D*  fV0CflatC4 ;
  TH1D*  fV0CflatS4 ;

  ClassDef(AliAnalysisTaskgg, 1); // PHOS analysis task
};

#endif
