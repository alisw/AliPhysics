
#ifndef AliAnalysisTaskSEDplusCorrelations_H
#define AliAnalysisTaskSEDplusCorrelations_H


#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3D.h>
#include <THnSparse.h>
#include <TArrayD.h>

#include "AliRDHFCutsDplustoKpipi.h"
#include "AliHFAssociatedTrackCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"
#include "AliEventPoolManager.h"
#include "AliNormalizationCounter.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliHFCorrelator.h"
#include "AliCentrality.h"
#include "AliHFOfflineCorrelator.h"

class TParticle ;
class TClonesArray ;
class AliAODMCParticle;
class AliAODEvent;
class AliVParticle;
class TObjArray;
class AliEventPoolManager;
class AliESDEvent;


class AliAnalysisTaskSEDplusCorrelations : public AliAnalysisTaskSE
{
    public :
    AliAnalysisTaskSEDplusCorrelations();
    AliAnalysisTaskSEDplusCorrelations(const Char_t* name, AliRDHFCutsDplustoKpipi* DplusCuts, AliHFAssociatedTrackCuts *AsscCuts);
    virtual ~AliAnalysisTaskSEDplusCorrelations();
    
    // Class Interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    
    // Setters.
    void SetCorrFormPart(Bool_t genMC){fMCParticle=genMC;}
    void SetCorrFormTrack(Bool_t reco){fRecoTrk=reco;}
    void SetDataOrMC(Bool_t readMC){fReadMC=readMC;}
    void SetEventMixing(Bool_t mixing){fMixing=mixing;}
    void SetCorrelator(Int_t number) {fAssoParType = number;} // select 1 for hadrons, 2 for Kaons, 3 for Kzeros
    void SetSystem(Bool_t system){fSystem=system;} // select between pp (kFALSE) or PbPb (kTRUE)
    //void SetEtaRange(Double_t etacorr) {fEtaRange=etacorr;}
    void SetUseBit(Bool_t bits=kTRUE){fUseBit=bits;}
    void SetTCConfig(Bool_t TCcong=kFALSE){fTCconfig=TCcong;}
    void SetTrackEffActive(Bool_t effTrack=kFALSE){fEffTrack=effTrack;}
    void SetDplusEffActive(Bool_t effDplus=kFALSE){fEffDplus=effDplus;}
    void SetUseCentrality(Bool_t flag, Int_t estimator, Bool_t flag2){fEvalCentrality=flag; fCentralityEstimator=estimator; fPoolbyCent=flag2;}
    void SetBinWidth(Float_t BinW){fBinWidth=BinW;}
    void SetMCGevEventType(Bool_t sel1=kFALSE){fMCGenEvType=sel1;}
    void SetPoolByPoolCorr(Bool_t sel2=kFALSE){fPoolByPool=sel2;}
    void SetCheckCutDistandChoice(Bool_t sel3=kFALSE, Bool_t sel4=kFALSE){fCheckCutDist=sel3;fRawCutQA=sel4;}
    void SetAODMismatchProtection(Int_t sel5=1) {fAODProtection=sel5;}
    void SetLeadPartCorrelation(Bool_t Sel){fLeadPartCorr = Sel;}
    void SetAutoSignalSBRange(Bool_t autosignalSBrange){fAutoSignalSBRange = autosignalSBrange;}
    
    //void SetUseDisplacement(Int_t m) {fDisplacement=m;} // select 0 for no displ, 1 for abs displ, 2 for d0/sigma_d0
    
    void SetMinDPt(Double_t minDPt){fMinDPt=minDPt;}
    Bool_t GetFillTrees() const {return fFillTrees;}
    void SetFillTrees(Bool_t fillTrees, Double_t fractAccME) {fFillTrees=fillTrees; fFractAccME=fractAccME;}
    
    //Setting left sideband bin ranges
    void SetLSBLowerUpperLim(Double_t* LSBLowLim, Double_t* LSBUppLim){
        for(int i=0;i<fNPtBins;i++){
            fLSBLowLim.push_back(LSBLowLim[i]);
            fLSBUppLim.push_back(LSBUppLim[i]);
        }
    }
    
    //Setting S+B bin ranges
    void SetSandBLowerUpperLim(Double_t* SandBLowLim, Double_t* SandBUppLim){
        for(int i=0;i<fNPtBins;i++){
            fSandBLowLim.push_back(SandBLowLim[i]);
            fSandBUppLim.push_back(SandBUppLim[i]);
        }
    }
    
    //Setting right sideband bin ranges
    void SetRSBLowerUpperLim(Double_t* RSBLowLim, Double_t* RSBUppLim){
        for(int i=0;i<fNPtBins;i++){
            fRSBLowLim.push_back(RSBLowLim[i]);
            fRSBUppLim.push_back(RSBUppLim[i]);
        }
    }
    
    
    private :
    
    AliAnalysisTaskSEDplusCorrelations(const AliAnalysisTaskSEDplusCorrelations &source);
    AliAnalysisTaskSEDplusCorrelations& operator=(const AliAnalysisTaskSEDplusCorrelations& source);
    
    //correlation methods
    void DoDplusCutDistFill(AliAODRecoDecayHF3Prong *d);
    void HistoNomenclature();
    void HadronCorrelations(AliAODRecoDecayHF3Prong* d, Int_t isDplus);
    void CorrelationNSparsePlots(AliAODRecoDecayHF3Prong *d, AliReducedParticle* track, Int_t iPtBin, Int_t origDplus, Double_t weightEff);
    Int_t CheckOriginPartOfDPlus(TClonesArray* arrayMC, AliAODMCParticle *mcDplus) const;
    
    //Offline
    void OfflineDPlusTree(AliAODRecoDecayHF3Prong* d, AliAODEvent* aod);
    void OfflineAssoTrackTree(AliAODEvent* aod);
    Bool_t AcceptTrackForMEOffline(Double_t TrackPt);
    
    
    //variables..
    Bool_t fSystem; // pp or PbPb
    Bool_t fReadMC; //  MC Switch
    Bool_t fRecoTrk; // Switch to reco track
    Bool_t fMCParticle; // Switch to Montecarlo Gen level
    Bool_t fMCGenEvType; //Gen MC event type
    TClonesArray* farrayMC; //! mcarray
    Bool_t fMixing;// switch for event mixing
    Int_t fAssoParType; // Correlation Option between D+ and (1-chargedtracks,2-chargedkaons,3-k0s )
    AliRDHFCutsDplustoKpipi *fDplusCuts;  // Cuts D+
    AliHFAssociatedTrackCuts *fAssoCuts; // cuts for associated track
    Bool_t fEffTrack; //Track eff ON/OFF
    Bool_t fEffDplus; //Dplus eff ON/OFF
    Int_t  fCentralityEstimator;   // enum from AliRDHFCuts..
    Bool_t    fEvalCentrality; // Switch to ON/OFF the centrality interface
    Double_t  fMinCentrality; // Minimun Centrality Value
    Double_t  fMaxCentrality; // Maximum Centrality Value
    Double_t fCentrOrMult; // Multiplicity of Event for D eff
    Bool_t fTCconfig; //  TC Cuts option
    Bool_t fUseBit; //  filterbit option
    Bool_t fAutoSignalSBRange; // mainly for offline correlation
    AliHFCorrelator  *fCorrelator; //object for correlations
    Int_t fNPtBins; // number of event at different Stages
    TH1F *fHistNEvents; //!hist. for No. of events
    TH1F *fHistNDplus; //!hist. for No. of Dplus
    TH1F *fHistNDTrkOff; //!hist. for No. of Dplus
    TH1F *fHistNDDauOnRemoved; //!hist. for No. of Dplus
    TH1F *fHistNDDauTrigID; //!hist. for No. of Dplus
    AliNormalizationCounter *fCounter; // counter
    Float_t fBinWidth;//width of one bin in output histos
    Bool_t  fPoolByPool;
    Int_t  fWhichPool;
    Bool_t fPoolbyCent;
    Int_t fEvtMult;
    Bool_t fCheckCutDist;
    Int_t fAODProtection; //New by Fabio
    TString fCutSuffix; //suffix for cut
    Bool_t fRawCutQA; //if D cut before sel
    TList *fOutput;     //! user output data
    TList *fOutputCorr; //! user output data
    Bool_t fLeadPartCorr; // Added by shyam Flag for leading particle correlation THnsparse
    //Offline
    AliHFCorrelationBranchD   *fBranchD; //!
    AliHFCorrelationBranchTr  *fBranchTr; //!
    TTree *fTreeD;    //tree for D+ mesons
    TTree *fTreeTr;   //tree for Assoc tracks
    Bool_t    fFillTrees;  //Flag to fill ME offline trees
    Double_t  fFractAccME; //Fraction of tracks to be accepted in the ME offline
    
    Int_t fNtrigDplusInR; // # of D+ filled (for association with decay tracks in TTrees)
    Int_t fNtrigDplusOutR; // # of D+ filled (for association with decay tracks in TTrees)
    Bool_t    fAlreadyFilled;            // D+ in an event already analyzed (for track distribution plots)
    Double_t  fMinDPt;                   // Minimum pT of the D+ to allow selection
    TObjArray *fTrackArray;  // Array with selected tracks for association
    Bool_t fTrackArrayFilled; // Flag to fill fTrackArray or not (if already filled)
    Double_t  fzVtx; // event zVtx position (for track eff)
    
    std::vector<Int_t> fDaughTrackID;       // ID of tagged daughters
    std::vector<Int_t> fDaughTrigNum;       // ID of D-trigger for daughters
    std::vector<Double_t>  fLSBLowLim;      // Left SB lower lim
    std::vector<Double_t>  fLSBUppLim;      // Left SB lower lim
    std::vector<Double_t>  fSandBLowLim;    // Signal +Bkg lower lim
    std::vector<Double_t>  fSandBUppLim;    // Signal _Bkg upper lim
    std::vector<Double_t>  fRSBLowLim;      // Right SB upper lim
    std::vector<Double_t>  fRSBUppLim;      // Right SB upper lim
    
    ClassDef(AliAnalysisTaskSEDplusCorrelations,9); // class for D+ meson correlations
    
};

#endif
