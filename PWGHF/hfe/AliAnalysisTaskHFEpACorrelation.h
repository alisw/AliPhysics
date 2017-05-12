#ifndef AliAnalysisTaskHFEpACorrelation_cxx
#define AliAnalysisTaskHFEpACorrelation_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//      Task for Heavy-flavour electron analysis in pPb collisions    //
//      (+ Electron-Hadron Jetlike Azimuthal Correlation)             //
//																	  //
//		version: Nov 18, 2016.							              //
//                                                                    //
//	    Authors 							                          //
//	    Authors 							                          //
//		Elienos Pereira de Oliveira Filho (epereira@cern.ch)	      //
//		Cristiane Jahnke		(cristiane.jahnke@cern.ch)            //
//      Henrique Zanoli (h.zanoli@cern.ch)                            //
//      Alexis Mas (aleximas@if.usp.br)                               //
////////////////////////////////////////////////////////////////////////

class TH1F;
class TH2F;
class AliESDEvent;
class AliESDtrackCuts;
class AliESDtrack;
class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliCFManager;
class AliPIDResponse;
class AliCentrality;
class AliAODEvent;
class AliVEvent;
class AliAODMCHeader;
class AliSelectNonHFE;
class AliEventPoolManager;
class AliEventPool;
class TObjArray;
class TF1;
//Lucile
class AliCaloTrackAODReader;
class AliCaloTrackReader;
//exotic
class AliAODReader;
class AliCalorimeterUtils;
class AliAnalysisUtils;

// --- ROOT system ---
#include <TObject.h>
#include <TString.h>
#include <TH3F.h>
#include <TObjArray.h>
#include <TArrayF.h>
#include <TGeoMatrix.h>

//--- ANALYSIS system ---
class AliVEvent;
class AliVTrack;


//______________________________________________________________________
//Library
#include "AliAnalysisTaskSE.h"
#include "AliHFEpid.h"
#include "AliLog.h"
//______________________________________________________________________

//______________________________________________________________________
class AliAnalysisTaskHFEpACorrelation : public AliAnalysisTaskSE
{
    //______________________________________________________________________
public:
    
    AliAnalysisTaskHFEpACorrelation();
    AliAnalysisTaskHFEpACorrelation(const char *name);
    virtual ~AliAnalysisTaskHFEpACorrelation();
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    
    //Setters
    void SetAssHadronPtRange(Double_t AssHadronPtMin, Double_t AssHadronPtMax) {fAssHadronPtMin = AssHadronPtMin; fAssHadronPtMax = AssHadronPtMax; };
    void SetHFECuts(AliHFEcuts * const cuts) {fCuts = cuts;};
    void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) {fRejectKinkMother = rejectKinkMother;};
    void SetCorrelationAnalysis(Bool_t CorrelationFlag=kTRUE) {fCorrelationFlag = CorrelationFlag;};
    void SetMCanalysis() {fIsMC = kTRUE;};
    void SetUseKF(Bool_t use) {fUseKF = use;};
    void SetPPanalysis(Bool_t Ispp = kFALSE) {fIspp = Ispp;};
    void SetCentrality(Double_t CentralityMin, Double_t CentralityMax) { fCentralityMin = CentralityMin; fCentralityMax = CentralityMax; fHasCentralitySelection = kTRUE; };
    void SetAODanalysis(Bool_t IsAOD) {fIsAOD = IsAOD;};
    void SetEventMixing(Bool_t EventMixingFlag) { fEventMixingFlag = EventMixingFlag;};
    void SetNonHFEmassCut(Double_t MassCut) { fMassCut = MassCut; fMassCutFlag = kTRUE;};
    void SetEtaCut(Double_t EtaCutMin,Double_t EtaCutMax ) { fEtaCutMin = EtaCutMin; fEtaCutMax = EtaCutMax; };
    void SetpTBins(Int_t n, Float_t* array) { fpTBins.Set(n,array); };
    
    void SetNonHFEangleCut(Double_t AngleCut) { fAngleCut = AngleCut; fAngleCutFlag = kTRUE;};
    void SetNonHFEchi2Cut(Double_t Chi2Cut) { fChi2Cut = Chi2Cut; fChi2CutFlag = kTRUE;};
    void SetNonHFEdcaCut(Double_t DCAcut) { fDCAcut = DCAcut; fDCAcutFlag = kTRUE;};
    
    void SetEfficiencyHadron(TH3F *hMap){if(fEffHadron) delete fEffHadron;fEffHadron = (TH3F*)hMap->Clone();}
    
    void SetBackgroundPi0Weight(TH1F *hBkgPi0W) {if(fBkgPi0Weight) delete fBkgPi0Weight; fBkgPi0Weight = (TH1F*) hBkgPi0W->Clone("fBkgPi0Weight");}
    void SetBackgroundEtaWeight(TH1F *hBkgEtaW) {if(fBkgEtaWeight) delete fBkgEtaWeight; fBkgEtaWeight = (TH1F*) hBkgEtaW->Clone("fBkgEtaWeight");}
    
    void SetBackgroundPi0WeightToData(TH1F *Wpion) {fBkgPi0WeightToData = (TH1F*) Wpion->Clone("PionWToData");}
    void SetBackgroundEtaWeightToData(TH1F *WEta) {fBkgEtaWeightToData = (TH1F*) WEta->Clone("EtaWToData");}

    
    //DCA cut main particle
    void SetdcaCut(Double_t DCAcutr, Double_t DCAcutz) { fDCAcutr = DCAcutr; fDCAcutz = DCAcutz;};
    
    //DCA cut for hadrons
    void SetDCACutHadron(Double_t DCAcutr, Double_t DCAcutz) { fDCAcutrHadron = DCAcutr; fDCAcutzHadron = DCAcutz;};
    void SetUseDCACutHadron() {fUseDCACutforHadrons = kTRUE;};
    
    void UseGlobalTracksHadron() { fUseGlobalTracksHadron = kTRUE;};
    
    void SetUseAlternativeBinning() {fUseAlternativeBinnig = kTRUE;};
    
    void SetCentralityEstimator(Int_t Estimator) { fEstimator=Estimator; }; //0 = ZNA, 1 = V0A
    void SetAdditionalCuts(Double_t PtMinAsso, Int_t TpcNclsAsso) {fPtMinAsso = PtMinAsso; fTpcNclsAsso = TpcNclsAsso;};
    void SetSPDCutForHadrons() {fAssocWithSPD = kTRUE;};
    
    //Getters
    AliHFEpid *GetPID() const {return fPID;};
    //bad channel
    //AliEMCALGeometry * GetEMCALGeometry()              const { return fEMCALGeo; }
    //AliCalorimeterUtils * GetCaloUtils()               const { return fCaloUtils; }
    /*AliCalorimeterUtils * GetCaloUtils()                                { if(!fCaloUtils) fCaloUtils = new AliCalorimeterUtils();
     return fCaloUtils      ; }*/
    
    //______________________________________________________________________
    
    //______________________________________________________________________
private:
    
    enum CocktailType_t{
        kNoCoktail = 0,
        kHijing = 1,
        kHFEnhanced = 2,
        kBackgroundEnhanced = 3,
        kUndefined = 4
    } ;
    
    //Function to process track cuts
    Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);
    //Function to process eh analysis
    void ElectronHadronCorrelation(AliVTrack *track, Int_t trackIndex, AliVParticle *vtrack);
    //Function to find non-HFE and fill histos
    //Selected Hadrons, for mixed event analysis
    TObjArray* SelectedHadrons();
    //DiHadron Correlation Background
    void DiHadronCorrelation(AliVTrack *track, Int_t trackIndex);
    Double_t GetHadronEfficiency(Double_t pT, Double_t eta, Double_t zvtx);
    
    Double_t CalculateWeight(Int_t pdg_particle, Double_t x);
    Double_t CalculateWeightRun2(Int_t pdg_particle, Double_t pT);
    Double_t CalculateWeightRun2ToData(Int_t pdg_particle, Double_t pT);
    void FillHistBkgWtoData(AliAODMCParticle* MCMotheWtoData, AliVTrack* track);
    
    void ComputeWeightInEnhancedSample();
    CocktailType_t FindTrackGenerator(Int_t label, AliAODMCHeader *header,TClonesArray *arrayMC);
    void GetTrackPrimaryGenerator(Int_t lab,AliAODMCHeader *header,TClonesArray *arrayMC,TString &nameGen);
    TString GetGenerator(Int_t label, AliAODMCHeader* header);

    void TaggingEfficiencyCalculation(AliVTrack *track,Bool_t *lIsNHFe,Bool_t *lIsHFe,Bool_t *lIsOther, Bool_t *lHasMother);
    void TaggingEfficiencyCalculationRun1(AliVTrack *track, Bool_t *lIsNHFe,Bool_t *lIsHFe,Bool_t *lIsOther, Bool_t *lHasMother);
    void TaggingEfficiencyCalculationRun2(AliVTrack *track, Bool_t *lIsNHFe,Bool_t *lIsHFe,Bool_t *lIsOther, Bool_t *lHasMother);
    
    //Flags for specifics analysis
    Bool_t 				fCorrelationFlag; //
    Bool_t              fUseKF;

    Bool_t				fIspp; //pp analysis
    Bool_t				fIsMC; //MC flag
    //Check different binning
    Bool_t              fUseAlternativeBinnig; //Smaller binning for trigger electron
    Bool_t				fAssocWithSPD; // hadron with SPD cut
    
    //Used in the function FindMother
    Bool_t				fIsHFE1; //legacy
    Bool_t				fIsHFE2; //legacy
    Bool_t				fIsNonHFE; //legacy
    Bool_t				fIsFromD; //legacy
    Bool_t				fIsFromB; //legacy
    Bool_t				fIsFromPi0; //legacy
    Bool_t				fIsFromEta; //legacy
    Bool_t				fIsFromGamma; //legacy
    
    //General variables
    TArrayF                 fpTBins;
    AliESDEvent 			*fESD; //!
    AliAODEvent 		   	*fAOD;//!
    AliVEvent 		      	*fVevent; //!
    AliESDtrackCuts         *fPartnerCuts;//!
    TList       			*fOutputList;//!
    AliPIDResponse 			*fPidResponse;//!
    AliSelectNonHFE 		*fNonHFE;//!
    
    //For the case of AOD analysis
    Bool_t					fIsAOD;	//flag for AOD analysis
    
    //For Centrality Selection
    AliCentrality			*fCentrality; //!
    Double_t				fCentralityMin;
    Double_t				fCentralityMax;
    Double_t				fCentralityValue;
    Bool_t					fHasCentralitySelection;
    TH1F					*fCentralityHist; //!
    TH1F					*fCentralityHistPass; //!
    Float_t					fZvtx;
    Int_t					fEstimator;
    
    //New Hadron DCA cut and Efficiency dependence
    Bool_t                  fUseDCACutforHadrons;
    
    //Efficiency Maps
    TH3F                    *fEffHadron;
    
    //Histograms
    TH1F				*fNevent; //!
    TH1F				*fNevent2; //!
    TH1F				*fPtElec_Inc;//!
    TH1F				*fPtElec_ULS; //!
    TH1F				*fPtElec_LS; //!
    
    TH1F				*fPtElec_ULS_NoPid; //!
    TH1F				*fPtElec_LS_NoPid; //!
    
    TH1F				*fPtElec_ULS_weight; //!
    TH1F				*fPtElec_LS_weight; //!
    TH1F				*fPtElec_ULS2_weight; //!
    TH1F				*fPtElec_LS2_weight; //!
    
    //PID Histograms
    
    TH2F				*fTOF01; //!
    TH2F				*fTOF02; //!
    TH2F				*fTOF03; //!
    TH1F				*fpid; //!
    
    TH1F				**fTPC_pt; //!
    TH2F				**fTPC_p; //!
    
    TH2F				*fTPC_momentum; //!
    TH2F				*fTPC_eta; //!
    TH2F				*fTPC_momentum1; //!
    TH2F				*fTPC_eta1; //!
    
    TH1F				**fTPCnsigma_pt; //!
    TH2F				   **fTOFTPCnsigma_pt; //!
    TH2F				**fTPCnsigma_p; //!
    
    TH2F				*fTPCnsigma_eta; //!
    TH2F				*fTPCnsigma_phi; //!

    TH1F				**fVtxZ; //!
    TH1F				*fVtxZ_new1; //!
    TH1F				*fVtxZ_new2; //!
    TH1F				*fVtxZ_new3; //!
    TH1F				*fVtxZ_new4;  //!
    
    TH1F		        *fzRes1; //!
    TH1F		    	*fzRes2; //!
    TH1F		    	*fSPD_track_vtx1; //!
    TH1F			    *fSPD_track_vtx2; //!
    
    TH1F				**fEtad; //!
    TH1F				**fNTracks; //!
    TH1F				*fTrack_Multi; //!
    TH2F				**fNTracks_pt; //!
    TH2F				**fNTracks_eta; //!
    TH2F				**fNTracks_phi; //!
    TH2F				**fTPCNcls_pid; //!

    //Electron-Hadron Correlation Histograms
    TH2F				**fCEtaPhi_Inc; //!
    TH2F				**fCEtaPhi_ULS_Weight; //!
    TH2F				**fCEtaPhi_LS_Weight; //!
    TH2F				**fCEtaPhi_ULS_NoP_Weight; //!
    TH2F				**fCEtaPhi_LS_NoP_Weight; //!
    
    TH1F				**fInvMassULS; //!
    TH1F				**fInvMassLS; //!
    TH1F				*fDCA; //!
    TH1F				*fDCABack; //!
    TH1F				*fOpAngle; //!
    TH1F				*fOpAngleBack; //!
    TH1F				*fInvMass2; //!
    TH1F				*fInvMass2_weight; //!
    TH1F				*fInvMassBack2; //!
    TH1F				*fInvMassBack2_weight; //!
    
    TH1F				**fInvMass_pT; //!
    TH1F				**fInvMassBack_pT; //!
    
    
    TH1F				*fDCA2; //!
    TH1F				*fDCABack2; //!
    TH1F				*fOpAngle2; //!
    TH1F				*fOpAngleBack2; //!
    
    Double_t			fMassCut;
    Double_t			fEtaCutMin;
    Double_t			fEtaCutMax;
    Double_t            fMinpTElec;
    Double_t            fMaxpTElec;
    Double_t			fAngleCut;
    Double_t			fChi2Cut;
    Double_t			fDCAcut; //Background DCA Cut
    Double_t			fDCAcutr; //
    Double_t			fDCAcutz;
    Double_t            fDCAcutrHadron;
    Double_t            fDCAcutzHadron;
    
    Bool_t				fMassCutFlag;
    Bool_t				fAngleCutFlag;
    Bool_t				fChi2CutFlag;
    Bool_t				fDCAcutFlag;
    //Correlation Function
    Double_t			fAssHadronPtMin;
    Double_t			fAssHadronPtMax;
    //Non-HFE reconstruction efficiency
    TH1F				*fPtBackgroundBeforeReco;  //!

    Double_t			fPtMinAsso; //
    Int_t			    fTpcNclsAsso; //
    
    //For the HFE package
    AliHFEcuts 			*fCuts; // Cut Collection for HFE
    //Lucile
    AliCFManager 		*fCFM; //! Correction Framework Manager
    AliHFEpid 			*fPID; // PID
    AliHFEpidQAmanager 	*fPIDqa; //! PID QA manager
    
    //Others
    AliStack 			*fMCstack;	//!
    Bool_t              fRejectKinkMother;	//!
    TParticle 			*fMCtrack; //!
    TParticle 			*fMCtrackMother; //!
    TParticle 			*fMCtrackGMother; //!
    TParticle 			*fMCtrackGGMother; //!
    TParticle 			*fMCtrackGGGMother; //!
    TClonesArray 		*fMCarray; //!
    AliAODMCHeader 		*fMCheader; //!
    AliAODMCParticle 	*fMCparticle; //!
    AliAODMCParticle 	*fMCparticle2; //!
    AliAODMCParticle 	*fMCparticleMother; //!
    AliAODMCParticle 	*fMCparticleGMother; //!
    AliAODMCParticle 	*fMCparticleGGMother; //!
    AliAODMCParticle 	*fMCparticleGGGMother; //!
    AliMCEventHandler	*fEventHandler; //!
    AliMCEvent			*fMCevent; //!
    
    //______________________________________________________________________
    //Mixed event analysis
    AliEventPoolManager *fPoolMgr; //!
    AliEventPool		*fPool; //!
    TObjArray			*fTracksClone; //!
    TObjArray			*fTracks; //!
    
    TH2F				**fCEtaPhi_Inc_EM; //!
    
    TH2F				**fCEtaPhi_ULS_Weight_EM; //!
    TH2F				**fCEtaPhi_LS_Weight_EM; //!
    
    TH2F				**fCEtaPhi_Inc_NoULSP; //!
    TH2F				**fCEtaPhi_Back_ULS_NoULSP; //!
    TH2F				**fCEtaPhi_Back_LS_NoULSP; //!
    
    TH1F				*fPoolNevents; //!
    
    Bool_t				fEventMixingFlag; //
    //______________________________________________________________________
    
    //______________________________________________________________________
    //Di-hadron correlation
    TH2F				**fCEtaPhi_Inc_DiHadron;  //!
    TH1F				*fPtTrigger_Inc;  //!
    AliAnalysisUtils *fAnalysisUtils;     //! Analysis Utils for pA pileup cut
    
    
    //MC CT Analysis
    TH1F				**fDCAElectronXY; //!
    TH1F				**fDCAElectronZ; //!
    
    //Generated pT,eta,zvtx distributions
    TH1F                *fNoEtaCutElectronGeneratedSignalPtEtaZvtx; //!
    TH1F                *fEtaCutElectronGeneratedSignalPtEtaZvtx;//!
    TH1F                *fEtaCutElectronInclusiveRecoPtEtaZvtx; //!
    TH1F                *fEtaCutElectronBKNoTag; //!
    TH1F                *fEtaCutElectronBKWithLabelULS; //!
    TH1F                *fEtaCutElectronBKWithLabelLS; //!
    TH1F                *fEtaCutElectronRecoHFEMC; //!
    TH1F                *fEtaCutElectronRecoOtherMC; //!
    TH1F                *fMissIDElectronsReco; //!
    TH3F                *fHadronsReco; //!
    TH3F                *fHadronsRecoPP; //!
    TH3F                *fHadronsGenerated; //!
    TH1F                *fElectronNoLabel; //!
    TH1F                *fElectronNoLabelULS; //!
    TH1F                *fElectronNoLabelLS; //!
    TH1F                *fEtaCutElectronBKULSMainSources; //!
    TH1F                *fEtaCutElectronBKLSMainSources; //!
    TH1F                *fEtaCutElectronBKULSOtherSources; //!
    TH1F                *fEtaCutElectronBKLSOtherSources; //!
    TH1F                *fEtaCutElectronHFEULS; //!
    TH1F                *fEtaCutElectronHFELS; //!
    TH1F                *fEtaCutElectronMissIDULS; //!
    TH1F                *fEtaCutElectronMissIDLS; //!
    
    //Test NHFE weight
    TH1F                *fEtaCutElectronBKULSMainSources_NW; //!
    TH1F                *fEtaCutElectronBKLSMainSources_NW; //!
    
    //Test Background weight (as Cris)
    TH1F                *fEtaCutElectronBKNoTag_WithMotherW; //!
    TH1F                *fEtaCutElectronBKULSMainSources_WithMotherW; //!
    TH1F                *fEtaCutElectronBKULSMainSources_WithMotherW_NW; //!
    TH1F                *fEtaCutElectronBKLSMainSources_WithMotherW; //!
    TH1F                *fEtaCutElectronBKLSMainSources_WithMotherW_NW; //!
    
    TH1F                *fElectronBKGNoEnhULS; //!
    TH1F                *fElectronBKGNoEnhLS; //!
    TH1F                *fElectronBKGNoEnhTotalNumber; //!
    TH1F                *fElectronBKGWToDataTotal; //!
    TH1F                *fElectronBKGWToDataULS; //!
    TH1F                *fElectronBKGWToDataLS; //!
    
    TH1F                *fElectronBKGNoEnhULS_WithW; //!
    TH1F                *fElectronBKGNoEnhLS_WithW; //!
    TH1F                *fElectronBKGNoEnhTotalNumber_WithW; //!
    
    //Background weight calculation
    
    TH1F                *fPtMCpi0_NoMother; //!
    TH1F                *fPtMCpi0_PureHijing; //!
    TH1F                *fPtMCEta_NoMother; //!
    TH1F                *fPtMCEta_PureHijing; //!
    
    TH1F                *fBkgPi0Weight; //
    TH1F                *fBkgEtaWeight; //
    TH1F                *fBkgPi0WeightToData; //
    TH1F                *fBkgEtaWeightToData; //
    
   
    
    //DPhi MC
    TH2F                **fCEtaPhiNoEtaCutInclusive;  //!
    TH2F                **fCEtaPhiNoEtaCutBKG; //!
    TH2F                **fCEtaPhiNoEtaCutHFe; //!
    TH2F                **fCEtaPhiNoEtaCutHFeNoDCA; //!
    TH2F                **fCEtaPhiNoEtaCutOther; //!
    TH2F                **fCEtaPhiNoEtaCutNoMother; //!
    
    TH2F                **fCEtaPhiCutInclusive; //!
    TH2F                **fCEtaPhiCutBKG; //!
    TH2F                **fCEtaPhiCutHFe; //!
    TH2F                **fCEtaPhiCutOther; //!
    TH2F                **fCEtaPhiCutNoMother; //!
    
    //Data With MC information
    TH2F                **fCEtaPhi_Back_MC_Tag; //!
    TH2F                **fCEtaPhi_Other_MC_Tag; //!
    TH2F                **fCEtaPhi_HFe_MC_Tag; //!
    TH2F                **fCEtaPhi_MC_NoMother_Tag; //!
    
    
    Bool_t fUseGlobalTracksHadron;
    
    
    //______________________________________________________________________
    
    AliAnalysisTaskHFEpACorrelation(const AliAnalysisTaskHFEpACorrelation&); 			// not implemented
    AliAnalysisTaskHFEpACorrelation& operator=(const AliAnalysisTaskHFEpACorrelation&); 		// not implemented
    
    ClassDef(AliAnalysisTaskHFEpACorrelation, 7); 								// example of analysis
    //______________________________________________________________________
};

///_________________________________________________________________________________________________
///Class copied from : $ALICE_ROOT/PWGCF/Correlations/DPhi/AliAnalysisTaskLongRangeCorrelations.h
///Author: Christoph Mayer
class AliHFEHCParticle : public TObject {
public:
    AliHFEHCParticle(Double_t eta=0, Double_t phi=0, Double_t pt=0, Double_t weight = 1)
    : fEta(eta), fPhi(phi), fPt(pt) , fWeight(weight){}
    virtual ~AliHFEHCParticle() {}
    
    Double_t Eta() const { return fEta; }
    Double_t Phi() const { return fPhi; }
    Double_t Pt() const { return fPt; }
    Double_t GetWeight() const { return fWeight; }
    
protected:
private:
    AliHFEHCParticle(const AliHFEHCParticle&);
    AliHFEHCParticle& operator=(const AliHFEHCParticle&);
    
    
    
    Double_t fEta;
    Double_t fPhi;
    Double_t fPt;
    Double_t fWeight;
    
    ClassDef(AliHFEHCParticle, 1);
} ;
///_________________________________________________________________________________________________

#endif
