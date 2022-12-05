/**************************************************************************
 *    Author:       Zuzana Moravcova                                      *
 *    Framework for calculating di-hadron correlation                     *
 *    for extraction of v_n{2} and v_n[2] coefficients.                   *
 *                                                                        *
 *    If used, modified, or distributed,                                  *
 *    please aknowledge the author of this code.                          *
 **************************************************************************/

#ifndef ALIANALYSISTASKCORRFORFLOW_H
#define ALIANALYSISTASKCORRFORFLOW_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliLog.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TRandom.h"
#include "TList.h"
#include "TMath.h"
#include "THnSparse.h"
#include "TString.h"
#include <vector>
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliMultSelection.h"
#include "AliAODInputHandler.h"
#include "AliEventPoolManager.h"
#include "AliVEvent.h"
#include "AliTHn.h"

class AliAnalysisTaskCorrForFlow : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskCorrForFlow();
                                AliAnalysisTaskCorrForFlow(const char *name, Bool_t bUseEff);
        virtual                 ~AliAnalysisTaskCorrForFlow();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        //event and track selection
        void                    SetTrigger(AliVEvent::EOfflineTriggerTypes trigger) { fTrigger = trigger; }
        void                    SetFilterBit(UInt_t filter) { fFilterBit = filter; }
        void                    SetPtRangeTrig(Double_t min, Double_t max) {fPtMinTrig = min; fPtMaxTrig = max; }
        void                    SetPtRangeAss(Double_t min, Double_t max) {fPtMinAss = min; fPtMaxAss = max; }
        void                    SetAbsEta(Double_t etaAbs) {fAbsEtaMax = etaAbs; }
        void                    SetPhiStarCur(Double_t phiStar) {fMergingCut = phiStar; }
        void                    SetCentrality(TString cent, Double_t min = 0.0, Double_t max = 20.0) { fCentEstimator = cent; fCentMin = min; fCentMax = max; }
        void                    SetUseNchRange(Bool_t range, Int_t min, Int_t max) { fUseNch = range; fNchMin = min; fNchMax = max; }
        void                    SetPtBins(std::vector<Double_t> bins) { fPtBinsTrigCharged = bins; }
        void                    SetPtBinsAss(std::vector<Double_t> bins) { fPtBinsAss = bins; }
        void                    SetCentBinsForMixing(Int_t nofBins, std::vector<Double_t> bins) { fNCentBins = nofBins; fCentBins = bins; }
        void                    SetSystematicsFlag(TString flag) { fSystematicsFlag = flag; }
        void                    SetNofSamples(Int_t n) { fNOfSamples = n; } //sampling setter
        void                    SetEtaPolarity(Int_t polarity) {fEtaPolarity = polarity; }
             /* code */

        void                    SetPVZcut(Double_t cut) { fPVzCut = cut; } //systematics setters
        void                    SetutDCAz(Double_t cut) { fCutDCAz = cut; }
        void                    SetutDCAxy(Double_t cut) { fCutDCAxySigma = cut; }
        void                    SetCutTPCchi2pCl(Double_t cut) { fCutTPCchi2pCl = cut; }
        void                    SetTPCclMincut(Double_t cut) { fTPCclMin = cut; }

        void                    SetIsHMpp(Bool_t hm = kTRUE) { fIsHMpp = hm; }
        void                    SetUseEtaDependentEfficiencies(Bool_t ef = kTRUE) { fEfficiencyEtaDependent = ef; }

    private:

        void                    PrintSetup();
        Bool_t                  IsEventSelected();
        Bool_t                  IsEventRejectedAddPileUp() const;
        Bool_t                  IsTrackSelected(const AliAODTrack* track) const;
        Bool_t                 IsTrackSelectedPolarity(const AliAODTrack* track) const;

        Double_t                RangePhi(Double_t DPhi);
        Double_t                GetDPhiStar(Double_t phi1, Double_t pt1, Double_t charge1, Double_t phi2, Double_t pt2, Double_t charge2, Double_t radius);
        void                    FillCorrelations();
        void                    FillCorrelationsMixed();
        Bool_t                  AreEfficienciesLoaded();
        Double_t                GetEff(const Double_t dPt, const Double_t dEta = 0.0);
        Int_t                   GetEtaRegion(const Double_t dEta);


        AliAnalysisTaskCorrForFlow(const AliAnalysisTaskCorrForFlow&); // not implemented
        AliAnalysisTaskCorrForFlow& operator=(const AliAnalysisTaskCorrForFlow&); // not implemented

        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputListCharged;    //! output list
        TList*                  fInputListEfficiency;    //! input list
        TObjArray*              fTracksTrigCharged; //!
        TObjArray*              fTracksAss; //!
        AliEventPoolManager*    fPoolMgr;  //!  event pool manager for Event Mixing
        //output histograms
        TH1D*                   fhEventCounter; //!
        TH1D*                   fhEventMultiplicity; //!
        TH2D*                   fHistPhiEta; //!
        TH3D*                   fhTrigTracks; //!
        AliTHn*                 fhChargedSE; //!
        AliTHn*                 fhChargedME; //!
        TH2D*                   fhEfficiencyEta[8]; //! eta dependent (4 sectors)

        //event and track selection
        AliVEvent::EOfflineTriggerTypes    fTrigger;
        Bool_t                  fIsHMpp; // [kFALSE]
        Bool_t                  fUseNch; // [kFALSE]
        Bool_t                  fUseEfficiency; // [kFALSE]
        Bool_t                  fEfficiencyEtaDependent; // [kFALSE]
        UInt_t                  fFilterBit;
        Int_t                   fbSign;
        Int_t                   fNofTracks;
        Int_t                   fNchMin;
        Int_t                   fNchMax;
        Double_t                fPtMinTrig;
        Double_t                fPtMaxTrig;
        Double_t                fPtMinAss;
        Double_t                fPtMaxAss;
        std::vector<Double_t>   fPtBinsTrigCharged;
        std::vector<Double_t>   fPtBinsAss;
        Double_t                fCentMin;
        Double_t                fCentMax;
        Double_t                fCentrality;
        Double_t                fAbsEtaMax;
        Double_t                fPVz;
        TString                 fCentEstimator;
        AliEventCuts            fEventCuts;

        // mixing
        Int_t                   fPoolMaxNEvents;   // maximum number of events in the pool
        Int_t                   fPoolMinNTracks;   // minimum number of tracks in the pool
        Int_t                   fMinEventsToMix;   // minimum number of events for mixing
        Int_t                   fNzVtxBins; // number of PV z bins
        Int_t                   fNCentBins; // number of centrality bins
        std::vector<Double_t>   fzVtxBins;
        std::vector<Double_t>   fCentBins;
        Double_t                fMergingCut; // [0.02] cut for track spliting/merging
        Int_t                   fNOfSamples;
        Int_t                   fSampleIndex;
        std::vector<Double_t>   fsampleBins; //sampling
        TString                 fSystematicsFlag; // ""
        Int_t                   fEtaPolarity; // trig particles in one half of detector, -1,0,1

        Double_t                fPVzCut;
        Double_t                fCutDCAz;
        Double_t                fCutDCAxySigma;
        Double_t                fCutTPCchi2pCl;
        Double_t                fTPCclMin;


        ClassDef(AliAnalysisTaskCorrForFlow, 8);
};

#endif
