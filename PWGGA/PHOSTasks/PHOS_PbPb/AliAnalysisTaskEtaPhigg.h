#ifndef AliAnalysisTaskEtaPhigg_cxx
#define AliAnalysisTaskEtaPhigg_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class THashList ;
class AliPHOSGeometry;
class AliCaloPhoton ;
class AliAODTrack ;
class AliEPFlattener ;
class AliV0ReaderV1 ;
class AliConvEventCuts ;
class AliConversionPhotonCuts ;
class AliAODConversionPhoton ;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEtaPhigg : public AliAnalysisTaskSE {
public:
  
  enum CutList {kDefault, kDisp, kCPV, kBoth, kDistance1, kDistance2, kDistance3} ;  
  
  
  AliAnalysisTaskEtaPhigg(const char *name = "AliAnalysisTaskEtaPhigg");
  virtual ~AliAnalysisTaskEtaPhigg() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
 
  void SetMassWindow(Float_t MinMass, Float_t MaxMass) { fMinMass = MinMass; fMaxMass = MaxMass; }
  void SetKappaWindow(Float_t MinKappa, Float_t MaxKappa) { fMinKappa = MinKappa; fMaxKappa = MaxKappa; }
  void SetEventCutList(Int_t nCuts, TList *CutArray){
        fnCuts = nCuts;
        fEventCutArray = CutArray;
    }
  void SetConversionCutList(Int_t nCuts, TList *CutArray){
        fnCuts = nCuts;
        fCutArray = CutArray;
    }
  void SetDoMesonAnalysis(Bool_t flag){fDoMesonAnalysis = flag;}
  void SetDoMesonQA(Int_t flag){fDoMesonQA = flag;}
  void SetDoPhotonQA(Int_t flag){fDoPhotonQA = flag;}
  void SetIsHeavyIon(Int_t flag){fIsHeavyIon = flag; }
  
  
  
private:
  AliAnalysisTaskEtaPhigg(const AliAnalysisTaskEtaPhigg&); // not implemented
  AliAnalysisTaskEtaPhigg& operator=(const AliAnalysisTaskEtaPhigg&); // not implemented

  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z, Double_t w) const ; //Fill 3D histogram witn name key

  
  Int_t ConvertRunNumber(Int_t run) ; 
  Bool_t PairCut(const AliCaloPhoton * ph1, const AliCaloPhoton * ph2, Int_t cut) const ; 
  Bool_t PHOSCut(const AliCaloPhoton * ph1, Int_t cut) const ; 
  Bool_t SecondaryPi0Cut(const AliCaloPhoton * ph1, const AliCaloPhoton * ph2) const ;
  void ProcessPCMPhotonCandidates() ;

private:
//  AliStack * fStack ;
  THashList *   fOutputContainer;        //final histogram container
  AliAODEvent * fEvent ;        //!
//  TClonesArray * fStack ;  
  TList *       fPHOSEvents[10][10][1] ; //Containers for events with PHOS photons
  TList *       fPCMEvents[10][10][1] ; //Containers for events with PHOS photons
  TClonesArray* fHadrEvents[10][10][1] ; //Containers for events with PHOS photons
  TClonesArray* fPHOSEvent ;      //PHOS photons in current event
  TClonesArray* fPCMEvent ;       //PCM photons in current event
  TClonesArray* fHadrEvent;       //hadrons
  
  //Reaction plain for v2
  AliEPFlattener * fV0AFlat ; //!
  AliEPFlattener * fV0CFlat ; //!
  Double_t         fRP ;      //!

  Int_t fRunNumber ;    //Current run number
  Float_t fCentrality ; //!Centrality of the currecnt event
  Int_t fCenBin ;       //! Current centrality bin

  AliPHOSGeometry  *fPHOSGeo;  //! PHOS geometry
  Int_t fEventCounter;         // number of analyzed events

  AliV0ReaderV1	*fV0Reader;											//
  Bool_t         fIsFromMBHeader;									//
  Bool_t         fDoMesonAnalysis;									//
  Int_t          fDoMesonQA;											//
  Int_t          fDoPhotonQA;										//
  Int_t          fnCuts;												//
  Int_t          fiCut;												//
  Int_t          fIsHeavyIon;	
  Float_t        fPtGamma;											//
  Float_t        fDCAzPhoton;										//
  Float_t        fRConvPhoton;										//
  Float_t        fEtaPhoton;											//
  UChar_t        fiCatPhoton;											//
  UChar_t        fiPhotonMCInfo; 										//
  Float_t               fMinMass;                        //
  Float_t               fMaxMass;                        //
  Float_t               fMinKappa;                       //
  Float_t               fMaxKappa;                       //
  TList                   *fEventCutArray;									//
  AliConvEventCuts	  *fEventCuts;										//
  TList                   *fCutArray;											//
  AliConversionPhotonCuts *fConversionCuts;									//
  TClonesArray            *fGammaCandidates;									//
  TList                  **fCutFolder;										//
  
  
  ClassDef(AliAnalysisTaskEtaPhigg, 1); // PHOS analysis task
};

#endif
