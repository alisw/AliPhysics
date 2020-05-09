#ifndef AliAnalysisPhotonDDA_h
#define AliAnalysisPhotonDDA_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

// Analysis task for identified tracks matched to a PHOS cluster
// Authors: Yuri Kharlov
// 14-Oct-2012

class TObjArray;
class TH1F;
class TH2I;
class TH2F;
class TH3F;
class AliPHOSGeometry;
class AliPIDResponse;
class AliAODCaloCluster ;
class AliAODCaloCells ;

#include "TH2I.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisPhotonDDA : public AliAnalysisTaskSE {
public:
  AliAnalysisPhotonDDA(const char *name = "AliAnalysisPhotonDDA");
  virtual ~AliAnalysisPhotonDDA() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetCentralityEstimator(Int_t i){fCentEstimator=i;}
  void SetMultiplicityBins(TArrayI *ar){fNCenBin=ar->GetSize() ; fCenBinEdges.Set(ar->GetSize(),ar->GetArray());} 

  void SetDistanceToBad(Double_t d=2.5){fMinBCDistance=d;} 
  void SetMC(Bool_t a){fIsMC=a;}
  
private:
  AliAnalysisPhotonDDA(const AliAnalysisPhotonDDA&); // not implemented
  AliAnalysisPhotonDDA& operator=(const AliAnalysisPhotonDDA&); // not implemented
  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  Int_t FindTrackMatching(Int_t mod,TVector3 *locpos,
					    Double_t &dx, Double_t &dz,
					    Double_t &pt,Int_t &charge) ;
  Double_t TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge) ;

 
private:
  THashList * fOutputContainer ;       // final histogram container  
  Int_t   fCenBin ;                    // centrality bin
  Int_t   fCentEstimator;              // Centrality estimator: 1: V0A/C, 2: V0M, 3: ZNA/C,  4: CL1
  Int_t   fNCenBin ;                   // Number of centrality bins
  TArrayI fCenBinEdges;                //Centrality binning
  Double_t     fMinBCDistance;         // Cut on distance to bad
  Double_t     fCentrality ;           // centrality
  Int_t fRunNumber ;
  Bool_t fIsMC ;
  AliPHOSGeometry  *fPHOSGeo ;         //! PHOS geometry
  Int_t fEventCounter;                 // number of analyzed events
  AliPIDResponse     *fPIDResponse ;   //! PID response 
  TClonesArray * fPHOSEvent  ;
  TList        * fCurrentMixedList ;
  TList        * fPHOSEvents[5] ;
  TH2F         * fhCont2D[15][4] ;
  ClassDef(AliAnalysisPhotonDDA, 1);   // PHOS DDA analysis task
};

#endif
