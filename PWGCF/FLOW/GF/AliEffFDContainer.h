/*
  Container to calculate charged-particle efficiencies with PCC and also record DCAxy distributions for feed-down estimation.
  The PCC framework (AliMCSpectraWeights) used here is written by Patrick Huhn.
  For postprocessing of output to get efficiencies and feed-down corrections, use the
  postprocessing macro at: https://github.com/vvislavi/Feeddown
  If used, please acknowledge the authors of AliMCSpectraWeights as well as myself
  Author: Vytautas Vislavicius
*/
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
#include "TMath.h"
class AliEffFDContainer: public TNamed {
  public:
    AliEffFDContainer();
    AliEffFDContainer(TString lName, TString lTitle, Bool_t lIsMC=kFALSE);
    ~AliEffFDContainer();
    Long64_t Merge(TCollection *collist);
    //Initialization settings
    void SetCentralityEstimator(TString newEst) { fCentEst=newEst; };
    void SetCentralityBins(Int_t nBins, Double_t *lbins) { StoreBins(nBins,lbins,fNCentBins, fCentBins); };
    void SetPtBins(Int_t nBins, Double_t *lbins) { StoreBins(nBins,lbins,fNPtBins, fPtBins); fPtMin=fPtBins[0]; fPtMax=fPtBins[fNPtBins];};
    void SetEta(Double_t newval) {fEta = newval; };
    void SetEta(Double_t etaLow, Double_t etaHigh) {fEtaLow = etaLow, fEta = etaHigh; };
    void AddCut(AliESDtrackCuts *inCuts);
    void AddCut(Int_t lFilterBit);
    void Fill(AliESDEvent &inputESD, AliMCEvent &inputMC);
    void Fill(AliESDEvent &inputESD);
    TList *GetOutList() { return fOutList; };
    Bool_t AddContainer(AliEffFDContainer *target);
    TH1* fetchObj(TString inname) { return (TH1*)fOutList->FindObject(makeName(inname).Data()); };
  // private:
    //Helper functions
    void NewEvent(AliESDEvent &inputESD);
    void NewEvent(AliMCEvent &inputMC);
    void StoreBins(Int_t nBins, Double_t *lBins, Int_t &tNBins, Double_t *&tBins);
    void CreateHistograms(Bool_t forceRecreate=kFALSE);
    Double_t GetChi2TPCConstrained(const AliESDtrack *l_Tr);
    Bool_t CheckEta(Double_t &lEta) { if(fEtaLow>-999) return ((lEta>fEtaLow) && (lEta<fEta)); else return (TMath::Abs(lEta)<fEta);  };
    void SetIdentifier(TString newname) { fIdentifier->SetTitle(newname.Data()); };
    TString makeName(TString pf) { return pf+fIdentifier->GetTitle(); };
    //Members
    TList *fOutList;
    TList *fCutList; //! might be interesting to store for reference, but irrelevant otherwise
    Double_t fChi2Cut;
    TF1 *fDCAXYPtCut;
    Bool_t fIsMC;
    Bool_t fInitialized;
    AliMCEvent *flMCEvent; //! fetched runtime, no need to store
    AliESDEvent *flESDEvent; //! cast at NewEvent(), so that we don't need to do that multiple times
    AliMultSelection *flMultSel;//! do not store
    AliMCSpectraWeightsHandler* flmcWeightsHandler;//! do not store this
    AliMCSpectraWeights *flMCSpectraWeights;//! do not store
    //Pointers to histograms -- so that we don't need to look them up in the list all the time:
    TH2D **fEff; //! Stored by TList
    TH3D **fDCA;//! Stored by TList
    TH2D **fWithinDCA;//! Stored by TList
    TNamed *fIdentifier; //
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
    Double_t fEtaLow;
    ClassDef(AliEffFDContainer,2);
};
#endif
