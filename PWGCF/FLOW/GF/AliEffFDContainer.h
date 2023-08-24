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
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliGFWFilter.h"

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
    void EnablePID(Bool_t newval) { fAddPID = newval; };
    void SetPIDObjects(AliPIDResponse *lPIDResponse, AliPIDCombined *lPIDCombined) {fPIDResponse = lPIDResponse; fBayesPID = lPIDCombined; if(lPIDCombined&&lPIDResponse) fAddPID=kTRUE; };
    void Fill(AliESDEvent &inputESD, AliMCEvent &inputMC);
    void Fill(AliESDEvent &inputESD);
    void Fill(AliAODEvent &inputAOD, AliMCEvent &inputMC);
    TList *GetOutList() { return fOutList; };
    Bool_t AddContainer(AliEffFDContainer *target);
    TH1* fetchObj(TString inname, Int_t iSpecie=0) { return (TH1*)fOutList->FindObject(makeName(inname,iSpecie).Data()); };
    void SetUseGenPt(Bool_t newval) { fUseGenPt = newval; };
    void SetBayesianProbs(std::vector<Double_t> probs) {fMinBayesProb.clear(); for(auto i: probs) fMinBayesProb.push_back(i); };
    void SetAODSelectionFlags(UInt_t evFlag, UInt_t trFlags) { fEvNomFlag=evFlag; fTrNomFlag=trFlags; };
    TH2 *getPureEfficiency2D(Int_t iSpecie, Bool_t weighted=kTRUE) {if(weighted) return                   get2DRatio("nChRec_Weighted","nChGen_Weighted",iSpecie); return get2DRatio("nChRec_Uneighted","nChGen_Uneighted",iSpecie); };
    TH1 *getPureEfficiency1D(Int_t iSpecie, Bool_t weighted=kTRUE, Int_t yb1=-1, Int_t yb2=-1) {if(weighted) return get1DRatio("nChRec_Weighted","nChGen_Weighted",iSpecie,yb1,yb2); return get1DRatio("nChRec_Uneighted","nChGen_Uneighted",iSpecie,yb1,yb2); };
    TH2 *getPurity2D(Int_t iSpecie, Bool_t weighted=kTRUE) {if(weighted) return get2DRatio("nChRec_Weighted","PurityPrimary_Weighted",iSpecie); return get2DRatio("nChRec_Uneighted","PurityPrimary_Uneighted",iSpecie); };
    TH1 *getPurity1D(Int_t iSpecie, Bool_t weighted=kTRUE, Int_t yb1=-1, Int_t yb2=-1) {if(weighted) return get1DRatio("nChRec_Weighted","PurityPrimary_Weighted",iSpecie,yb1,yb2); return get1DRatio("nChRec_Uneighted","PurityPrimary_Uneighted",iSpecie,yb1,yb2); };
    TH2 *getFDvsPhi(Int_t iSpecie, Bool_t ratioToIntegrated=kTRUE, Int_t cBin1=-1, Int_t cBin2=-1);
    TH1 *getPureFeeddown(Int_t iSpecie, Int_t centBin1=-1, Int_t centBin2=-1);
    Int_t getNCentBins() { TH2 *hTemp = (TH2*)fetchObj("nChRec_Weighted"); return hTemp->GetNbinsY(); };
      // private:
    //Helper functions
    void NewEvent(AliESDEvent &inputESD);
    void NewEvent(AliMCEvent &inputMC);
    void NewAODEvent(AliAODEvent &inputAOD, AliMCEvent &inputMC);
    void StoreBins(Int_t nBins, Double_t *lBins, Int_t &tNBins, Double_t *&tBins);
    void CreateHistograms(Bool_t forceRecreate=kFALSE);
    Double_t GetChi2TPCConstrained(const AliESDtrack *l_Tr);
    Bool_t CheckEta(Double_t &lEta) { if(fEtaLow>-999) return ((lEta>fEtaLow) && (lEta<fEta)); else return (TMath::Abs(lEta)<fEta);  };
    void SetIdentifier(TString newname) { fIdentifier->SetTitle(newname.Data()); };
    Int_t GetBayesPIDIndex(AliVTrack *l_track);
    Int_t GetTruePIDIndex(const Int_t &pdgcode);
    Int_t GetMCSWPrimIndex(AliMCParticle *lPart);
    Int_t GetMCSWMotherIndex(AliMCParticle *lPart, Double_t &lMotherPt);
    void SetMCSWeights(TH3D *inh) {fMCSWeights=inh;};
    Int_t CalculateMult(); //Required to pick up the correct weights
    TString getSpecieName(Int_t ind) {if(ind>=(Int_t)fSpNames.size() || ind<0) return "Undefined_"; return fSpNames[ind];};
    TString makeName(TString pf, Int_t spInd=0) { return getSpecieName(spInd)+pf+fIdentifier->GetTitle(); };
    //Getters for "pure" observables:
    TH2 *get2DRatio(TString numID, TString denID, Int_t iSpecie);
    TH1 *get1DRatio(TString numID, TString denID, Int_t iSpecie, Int_t yb1=-1, Int_t yb2=-1);
    TH1 *getPureFeeddown(Int_t iSpecie, TH2 *inh);
    //Members
    TList *fOutList;
    TList *fCutList; //! might be interesting to store for reference, but irrelevant otherwise
    Double_t fChi2Cut;
    TF1 *fDCAXYPtCut;
    UInt_t fEvNomFlag; //Flag for event selection -- relevant for AODs. The whole flag, not bit index
    UInt_t fTrNomFlag; //Flag for track selection -- relevant for AODs. The whole flag, not bit index
    Bool_t fAddPID;
    Int_t fNSpecies;
    AliPIDResponse *fPIDResponse; //! for PID
    AliPIDCombined *fBayesPID; //! for PID
    std::vector<Double_t> fMinBayesProb; //Minimum Bayesian probabilities
    Bool_t fIsMC;
    Bool_t fUseGenPt;
    Bool_t fInitialized;
    AliMCEvent *flMCEvent; //! fetched runtime, no need to store
    AliESDEvent *flESDEvent; //! cast at NewEvent(), so that we don't need to do that multiple times
    AliAODEvent *flAODEvent; //! cast at NewEvent(), so that we don't need to do that multiple times
    AliMultSelection *flMultSel;//! do not store
    AliMCSpectraWeightsHandler* flmcWeightsHandler;//! do not store this
    AliMCSpectraWeights *flMCSpectraWeights;//! do not store
    TH3D *fMCSWeights; //! do not store... I think?
    //Pointers to histograms -- so that we don't need to look them up in the list all the time:
    TH2D ***fEff; //! Stored by TList
    TH3D ***fDCA;//! Stored by TList
    TH3D ***fWithinDCA;//! Stored by TList
    TH2D ***fPurity;//! Stored by TList
    TH3D ***fFDvsPhi;//! Stored by TList
    TNamed *fIdentifier; //
    std::vector<TString> fSpNames;//! No need to store, defined at initialization
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
    ClassDef(AliEffFDContainer,6);
};
#endif
