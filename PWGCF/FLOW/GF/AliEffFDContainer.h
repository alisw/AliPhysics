#ifndef ALIEFFFDCONTAINER__H
#define ALIEFFFDCONTAINER__H
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TList.h"
#include "TNamed.h"
#include "AliESDtrackCuts.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliMCSpectraWeights.h"
#include "TString.h"
#include "TF1.h"
#include "AliESDVertex.h"
class AliEffFDContainer: public TNamed {
  public:
    AliEffFDContainer();
    AliEffFDContainer(TString lName, TString lTitle);
    ~AliEffFDContainer() {};
    TList *fOutList;
    TList *fCutList; //! might be interesting to store for reference, but irrelevant otherwise
    Double_t fChi2Cut;
    TF1 *fDCAXYPtCut;
    Bool_t fIsMC;
    //Initialization settings
    void SetCentralityEstimator(TString newEst) { fCentEst=newEst; };
    void SetCentralityBins(Int_t nBins, Double_t *lbins) { StoreBins(nBins,lbins,fNCentBins, fCentBins); };
    void SetPtBins(Int_t nBins, Double_t *lbins) { StoreBins(nBins,lbins,fNPtBins, fPtBins); fPtMin=fPtBins[0]; fPtMax=fPtBins[fNPtBins];};
    void SetEta(Double_t newval) {fEta = newval; };
    void AddCut(AliESDtrackCuts *inCuts);
    void NewEvent(AliVEvent &InputEventAddress, AliMCEvent &MCEventAddress);
    void Fill(AliVEvent &InputEventAddress, AliMCEvent &MCEventAddress);
    void SetTest(TH1D &hadd) { htmp=&hadd; };
    void PrintHName() {printf("%s\n",htmp->GetName()); };
  // private:
    AliMCEvent *flMCEvent; //! fetched runtime, no need to store
    AliESDEvent *flESDEvent; //! cast at NewEvent(), so that we don't need to do that multiple times
    AliMultSelection *flMultSel;//! do not store
    AliMCSpectraWeightsHandler* flmcWeightsHandler;//! do not store this
    AliMCSpectraWeights *flMCSpectraWeights;//! do not store
    //Pointers to histograms -- so that we don't need to look them up in the list all the time:
    TH2D **fEff; //! Stored by TList
    TH3D **fDCA;//! Stored by TList
    TH2D **fWithinDCA;//! Stored by TList
    //Pointers to axes
    Double_t *fPtBins; //!
    Int_t fNPtBins; //!
    Double_t *fCentBins; //!
    Int_t fNCentBins; //!
    //Centrality
    TString fCentEst;
    Double_t fCent;//! irrelevant
    //Kinematic cuts
    Double_t fPtMin;
    Double_t fPtMax;
    Double_t fEta;
    //Helper functions
    void StoreBins(Int_t nBins, Double_t *lBins, Int_t &tNBins, Double_t *&tBins);
    void CreateHistograms();
    Double_t GetChi2TPCConstrained(const AliESDtrack *l_Tr);
    //Testing
    TH1D *htmp; //! written by the list
    ClassDef(AliEffFDContainer,1);
};
#endif
