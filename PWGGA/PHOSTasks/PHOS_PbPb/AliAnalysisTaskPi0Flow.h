#ifndef AliAnalysisTaskPi0Flow_cxx
#define AliAnalysisTaskPi0Flow_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TObjArray;
class TH1F;
class TH2I;
class TH2F;
class TH3F;
class TF1 ;
class AliStack ;
class AliESDtrackCuts;
class AliPHOSGeometry;
class AliESDEvent ;
class AliPHOSCalibData;
class AliESDtrack ;
class AliESDCaloCluster ;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPi0Flow : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskPi0Flow(const char *name = "AliAnalysisTaskPi0Flow");
  virtual ~AliAnalysisTaskPi0Flow() {}
  
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
  AliAnalysisTaskPi0Flow(const AliAnalysisTaskPi0Flow&); // not implemented
  AliAnalysisTaskPi0Flow& operator=(const AliAnalysisTaskPi0Flow&); // not implemented
  Bool_t IsGoodChannel(const char * det, Int_t mod,Int_t ix, Int_t iz); //Use addisional bad map for PHOS

  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key

  Bool_t TestLambda(Double_t pt,Double_t l1,Double_t l2) ;  //Evaluate Dispersion cuts for photons
  Bool_t TestLambda2(Double_t pt,Double_t l1,Double_t l2) ;  //Evaluate Dispersion cuts for photons
  Double_t TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge);
  Int_t ConvertRunNumber(Int_t run) ; 

  void OpenInfoCalbration(Int_t run); //V0 calibration
  void EvalV0ReactionPlane(AliESDEvent * event) ;
  Double_t ApplyFlattening(Double_t phi, Double_t c) ; //Apply centrality-dependent flattening
  Double_t ApplyFlatteningV0A(Double_t phi, Double_t c) ; //Apply centrality-dependent flattening
  Double_t ApplyFlatteningV0C(Double_t phi, Double_t c) ; //Apply centrality-dependent flattening

  Double_t CoreEnergy(AliESDCaloCluster * clu); 
  void    Reclusterize(AliESDCaloCluster * clu) ;
  Bool_t  AreNeibors(Int_t id1,Int_t id2) ;
private:
  AliESDtrackCuts *fESDtrackCuts; // Track cut
  AliStack * fStack ;
  TList * fOutputContainer;        //final histogram container
  TList * fPHOSEvents[1][10][11] ; //Containers for events with PHOS photons
  TClonesArray * fPHOSEvent ;      //PHOS photons in current event
  AliPHOSCalibData *fPHOSCalibData; // PHOS calibration object
  TF1 *fNonLinCorr;          // Non-linearity correction
 
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


//  TF1 * fRecent[5][12] ;//Recentering corrections

  Int_t fRunNumber ;    //Current run number
  Float_t fCentrality ; //!Centrality of the currecnt event

  Int_t fCenBin ;       //! Current centrality bin

  TH2I *fPHOSBadMap[6] ;    //Container for PHOS bad channels map

  AliPHOSGeometry  *fPHOSGeo;  //! PHOS geometry
  Int_t fEventCounter;         // number of analyzed events

  ClassDef(AliAnalysisTaskPi0Flow, 1); // PHOS analysis task
};

#endif
