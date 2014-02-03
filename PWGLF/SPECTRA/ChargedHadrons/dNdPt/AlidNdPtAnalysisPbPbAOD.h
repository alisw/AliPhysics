#ifndef AlidNdPtAnalysisPbPbAOD_H
#define AlidNdPtAnalysisPbPbAOD_H


//------------------------------------------------------------------------------
// AlidNdPtAnalysisPbPbAOD class used for dNdPt analysis in PbPb collision
// via AODs 
// 
// Author: P. Luettig, 15.05.2013
// last modified: 08.10.2013
//------------------------------------------------------------------------------



class iostream;

#include "AliAnalysisTaskSE.h"
#include "TObject.h"
#include "TList.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "THn.h"
#include "TClonesArray.h"

#include "TParticlePDG.h"
#include "TDatabasePDG.h"

#include "AliLog.h"
#include "AliCentrality.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"

#include "AliInputEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliExternalTrackParam.h"
#include "AliESDtrack.h"

#include "TSystem.h"
#include "TROOT.h"

class AlidNdPtAnalysisPbPbAOD : public AliAnalysisTaskSE {
  public :
    enum CheckQuantity { cqCrossedRows = 0, cqNcluster = 1, cqChi = 2, cqLength = 3 };
    enum KinematicQuantity { kqPt = 0, kqEta = 1, kqPhi = 2 };
    enum MaxCheckQuantity { cqMax = 4 };
    enum MaxKinematicQuantity { kqMax = 3 };
    
    AlidNdPtAnalysisPbPbAOD(const char *name = "dNdPtPbPbAOD");
    ~AlidNdPtAnalysisPbPbAOD();
    
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);
    
    // Set binning for Histograms (if not set default binning is used)
    void SetBinsMult(Int_t nbins, Double_t* edges) 			{ Printf("[I] Setting Mult Bins"); fMultNbins = nbins; fBinsMult = GetArrayClone(nbins,edges); }
    void SetBinsPt(Int_t nbins, Double_t* edges) 			{ Printf("[I] Setting pT Bins"); fPtNbins = nbins; fBinsPt = GetArrayClone(nbins,edges); }
    void SetBinsPtCorr(Int_t nbins, Double_t* edges) 			{ Printf("[I] Setting pTcorr Bins"); fPtCorrNbins = nbins; fBinsPtCorr = GetArrayClone(nbins,edges); }
    void SetBinsPtCheck(Int_t nbins, Double_t* edges) 			{ Printf("[I] Setting pTcheck Bins"); fPtCheckNbins = nbins; fBinsPtCheck = GetArrayClone(nbins,edges); }
    void SetBinsEta(Int_t nbins, Double_t* edges) 			{ Printf("[I] Setting Eta Bins"); fEtaNbins = nbins; fBinsEta = GetArrayClone(nbins,edges); }
    void SetBinsEtaCheck(Int_t nbins, Double_t* edges) 			{ Printf("[I] Setting EtaCheck Bins"); fEtaCheckNbins = nbins; fBinsEtaCheck = GetArrayClone(nbins,edges); }
    void SetBinsZv(Int_t nbins, Double_t* edges) 			{ Printf("[I] Setting Zv Bins"); fZvNbins = nbins; fBinsZv= GetArrayClone(nbins,edges); }
    void SetBinsCentrality(Int_t nbins, Double_t* edges) 		{ Printf("[I] Setting Cent Bins"); fCentralityNbins = nbins; fBinsCentrality = GetArrayClone(nbins,edges); }
    void SetBinsPhi(Int_t nbins, Double_t* edges) 			{ Printf("[I] Setting Phi Bins"); fPhiNbins = nbins; fBinsPhi = GetArrayClone(nbins,edges); }
    
    // set event cut variables
    void SetCutMaxZVertex( Double_t d)					{ fCutMaxZVertex = d; }
    Double_t GetCutMaxZVertex()						{ return fCutMaxZVertex; }
    
    // set track kinematic cut parameters
    void SetCutPtRange(Double_t ptmin, Double_t ptmax)			{ fCutPtMin = ptmin; fCutPtMax = ptmax; }
    Double_t GetCutPtMin()						{ return fCutPtMin; }
    Double_t GetCutPtMax()						{ return fCutPtMax; }
    
    void SetCutEtaRange(Double_t etamin, Double_t etamax)		{ fCutEtaMin = etamin; fCutEtaMax = etamax; }
    Double_t GetCutEtaMin()						{ return fCutEtaMin; }
    Double_t GetCutEtaMax()						{ return fCutEtaMax; }
    
    void EnableRelativeCuts()						{ Printf("[I] Relative Cuts enabled"); fUseRelativeCuts = kTRUE; }
    Bool_t AreRelativeCutsEnabled()					{ return fUseRelativeCuts; }
    
    // setter and getter track quality cut parameters
    void SetFilterBit(Int_t b)						{ fFilterBit = b; };
    Int_t GetFilterBit()						{ return fFilterBit; }
    
    void SetCutRequireTPCRefit(Bool_t *b) 				{ fCutRequireTPCRefit = b; } 
    Bool_t IsTPCRefitRequired() 					{ return fCutRequireTPCRefit; } 
    
    void SetCutRequireITSRefit(Bool_t *b) 				{ fCutRequireITSRefit = b; } 
    Bool_t IsITSRefitRequired() 					{ return fCutRequireITSRefit; } 
    
    void SetCutMinNClustersTPC(Double_t d)				{ fCutMinNumberOfClusters = d; }
    Double_t GetCutMinNClustersTPC()					{ return fCutMinNumberOfClusters; }
    
    void SetCutPercMinNClustersTPC(Double_t d)				{ Printf("[I] Take only %.2f%% tracks with most clusters", d*100.); fCutPercMinNumberOfClusters = d; }
    Double_t GetCutPercMinNClustersTPC()				{ return fCutPercMinNumberOfClusters; }
    
    void SetCutMinNCrossedRowsTPC(Double_t d) 				{ fCutMinNumberOfCrossedRows = d; }    
    Double_t GetCutMinNCrossedRowsTPC()					{ return fCutMinNumberOfCrossedRows; }
    
    void SetCutPercMinNCrossedRowsTPC(Double_t d) 			{ Printf("[I] Take only %.2f%% tracks with most crossedRows", d*100.); fCutPercMinNumberOfCrossedRows = d; }    
    Double_t GetCutPercMinNCrossedRowsTPC()				{ return fCutPercMinNumberOfCrossedRows; }
    
    void SetCutMinRatioCrossedRowsOverFindableClustersTPC(Double_t d) 	{ fCutMinRatioCrossedRowsOverFindableClustersTPC = d; }
    Double_t GetCutMinRatioCrossedRowsOverFindableClustersTPC()		{ return fCutMinRatioCrossedRowsOverFindableClustersTPC; }
    
    void SetCutLengthInTPCPtDependent()					{ fCutLengthInTPCPtDependent = kTRUE; }
    Bool_t DoCutLengthInTPCPtDependent()				{ return fCutLengthInTPCPtDependent; }
    
    void SetPrefactorLengthInTPCPtDependent(Double_t d)			{ fPrefactorLengthInTPCPtDependent = d; }
    Double_t GetPrefactorLengthInTPCPtDependent()			{ return fPrefactorLengthInTPCPtDependent; }
     
    void SetCutMaxChi2PerClusterTPC(Double_t d) 			{ fCutMaxChi2PerClusterTPC = d; }
    void SetCutMaxFractionSharedTPCClusters(Double_t d) 		{ fCutMaxFractionSharedTPCClusters = d; }
    void SetCutMaxDCAToVertexZ(Double_t d) 				{ fCutMaxDCAToVertexZ = d; }
    void SetCutMaxDCAToVertexXY(Double_t d) 				{ fCutMaxDCAToVertexXY = d; }
    void SetCutMaxChi2PerClusterITS(Double_t d) 			{ fCutMaxChi2PerClusterITS = d; }
    void SetCutDCAToVertex2D(Bool_t *b) 				{ fCutDCAToVertex2D = b; } 
    void SetCutRequireSigmaToVertex(Bool_t *b) 				{ fCutRequireSigmaToVertex = b; } 
    void SetCutMaxDCAToVertexXYPtDep(Double_t d0, Double_t d1, Double_t d2)
    {
      fCutMaxDCAToVertexXYPtDepPar0 = d0;
      fCutMaxDCAToVertexXYPtDepPar1 = d1;
      fCutMaxDCAToVertexXYPtDepPar2 = d2;
    }
    void SetCutAcceptKinkDaughters(Bool_t *b) 				{ fCutAcceptKinkDaughters = b; } 
    void SetCutMaxChi2TPCConstrainedGlobal(Double_t d) 			{ fCutMaxChi2TPCConstrainedGlobal = d; }
       
    // fill function for cross check histos
    Bool_t FillDebugHisto(Double_t *dCrossCheckVar, Double_t *dKineVar, Double_t dCentrality, Bool_t bIsAccepted);
    
    // fill function for cut settings
    void StoreCutSettingsToHistogram();
    
    // getter for DCA
    Bool_t GetDCA(const AliAODTrack *track, AliAODEvent *evt, Double_t d0z0[2]);
    
    THnSparseF * GetHistZvPtEtaCent() const { return fZvPtEtaCent; }
    TH1F * GetHistEventStatistics() const { return fEventStatistics; }
    
    const char * GetParticleName(Int_t pdg);
    
    AliGenHijingEventHeader* GetHijingEventHeader(AliAODMCHeader *header);
    AliGenPythiaEventHeader* GetPythiaEventHeader(AliAODMCHeader *header);
    
    
    Bool_t SetRelativeCuts(AliAODEvent *event);
    
    Bool_t IsTrackAccepted(AliAODTrack *tr, Double_t dCentrality, Double_t bMagZ);
    Bool_t IsMCTrackAccepted(AliAODMCParticle *part);
    
    Bool_t IsHijingParticle(const AliAODMCParticle *part, AliGenHijingEventHeader* hijingGenHeader);
    Bool_t IsPythiaParticle(const AliAODMCParticle *part, AliGenPythiaEventHeader* pythiaGenHeader);
    
    static Double_t* GetArrayClone(Int_t n, Double_t* source);
    
  private :
    
    // Output List
    TList	*fOutputList;
    
    // Histograms
    TH1F	*fPt; // simple pT histogramm
    TH1F	*fMCPt; // simple pT truth histogramm
    THnSparseF 	*fZvPtEtaCent; //-> Zv:Pt:Eta:Cent
    THnSparseF 	*fPhiPtEtaCent; //-> Phi:Pt:Eta:Cent
    THnSparseF 	*fPtResptCent; //-> 1/pt:ResolutionPt:Cent
    THnSparseF 	*fMCRecPrimZvPtEtaCent; //-> MC Zv:Pt:Eta:Cent
    THnSparseF 	*fMCGenZvPtEtaCent; //-> MC Zv:Pt:Eta:Cent
    THnSparseF 	*fMCRecSecZvPtEtaCent; //-> MC Zv:Pt:Eta:Cent, only secondaries
    THnSparseF 	*fMCRecPrimPhiPtEtaCent; //-> MC Phi:Pt:Eta:Cent
    THnSparseF 	*fMCGenPhiPtEtaCent; //-> MC Phi:Pt:Eta:Cent
    THnSparseF 	*fMCRecSecPhiPtEtaCent; //-> MC Phi:Pt:Eta:Cent, only secondaries
    TH1F	*fEventStatistics; // contains statistics of number of events after each cut
    TH1F        *fEventStatisticsCentrality; // contains number of events vs centrality, events need to have a track in kinematic range
    TH1F	*fMCEventStatisticsCentrality; // contains MC number of events vs centrality, events need to have a track in kinematic range
    TH1F	*fAllEventStatisticsCentrality; // contains number of events vs centrality, events need to be triggered
    TH2F	*fEventStatisticsCentralityTrigger; // contains number of events vs centrality in 1% bins vs trigger
    THnSparseF	*fZvMultCent; // Zv:Mult:Cent
    TH1F	*fTriggerStatistics; // contains number of events per trigger
    TH1F	*fMCTrackPdgCode; // contains statistics of pdg codes of tracks
    TH1F	*fMCTrackStatusCode; // contains statistics of status codes of tracks
    TH1F	*fCharge; // charge distribution in data
    TH1F	*fMCCharge; // charge distribution in MC
    TH2F	*fMCPdgPt; // PDGvs PT for MC Particles
    TH1F	*fMCHijingPrim; // number of particles, which are Hijing particles and primaries
    THnSparseF	*fDCAPtAll; //control histo: DCAz:DCAxy:pT:eta:phi for all reconstructed tracks
    THnSparseF	*fDCAPtAccepted; //control histo: DCAz:DCAxy:pT:eta:phi for all accepted reco tracks
    THnSparseF	*fMCDCAPtSecondary; //control histo: DCAz:DCAxy:pT:eta:phi for all accepted reco track, which are secondaries (using MC info)
    THnSparseF	*fMCDCAPtPrimary; //control histo: DCAz:DCAxy:pT:eta:phi for all accepted reco track, which are primaries (using MC info)
    THnF	*fCrossCheckAll[4]; //control histo: {CrossedRows,Ncluster,Chi} vs pT,eta,phi,Centrality for all tracks
    THnF	*fCrossCheckAcc[4]; //control histo: {CrossedRows,Ncluster,Chi} vs pT,eta,phi,Centrality after cuts
    TH1F	*fCutPercClusters; // control histo: number of clusters, where the relative cut has been set e-by-e
    TH1F	*fCutPercCrossed; // control histo: number of crossed rows, where the relative cut has been set e-by-e
    TH2F	*fCrossCheckRowsLength; // control histo: number of crossed rows vs length in TPC
    TH2F	*fCrossCheckClusterLength; // control histo: number of clusters vs length in TPC
    TH2F	*fCrossCheckRowsLengthAcc; // control histo: number of crossed rows vs length in TPC for all accepted tracks
    TH2F	*fCrossCheckClusterLengthAcc; // control histo: number of clusters vs length in TPC for all accepted tracks
    TH1F        *fCutSettings; // control histo: cut settings
    
    
    // global variables
    Bool_t fIsMonteCarlo;
    
    // event cut variables
    Double_t fCutMaxZVertex;
    
    // track kinematic cut variables
    Double_t fCutPtMin;
    Double_t fCutPtMax;
    Double_t fCutEtaMin;
    Double_t fCutEtaMax;
    
    // track quality cut variables
    Int_t	fFilterBit;
    Bool_t 	fUseRelativeCuts;
    Bool_t 	fCutRequireTPCRefit;
    Bool_t 	fCutRequireITSRefit;
    Double_t	fCutMinNumberOfClusters;
    Double_t	fCutPercMinNumberOfClusters;
    Double_t 	fCutMinNumberOfCrossedRows;
    Double_t 	fCutPercMinNumberOfCrossedRows;
    Double_t 	fCutMinRatioCrossedRowsOverFindableClustersTPC;
    Double_t 	fCutMaxChi2PerClusterTPC;
    Double_t 	fCutMaxFractionSharedTPCClusters;
    Double_t 	fCutMaxDCAToVertexZ;
    Double_t 	fCutMaxDCAToVertexXY;
    Double_t 	fCutMaxChi2PerClusterITS;
    Bool_t  	fCutDCAToVertex2D;
    Bool_t 	fCutRequireSigmaToVertex;
    Double_t 	fCutMaxDCAToVertexXYPtDepPar0;
    Double_t 	fCutMaxDCAToVertexXYPtDepPar1;
    Double_t 	fCutMaxDCAToVertexXYPtDepPar2;
    Bool_t 	fCutAcceptKinkDaughters;
    Double_t 	fCutMaxChi2TPCConstrainedGlobal;
    Bool_t	fCutLengthInTPCPtDependent;
    Double_t	fPrefactorLengthInTPCPtDependent;
    
    //binning for THNsparse
    Int_t   fMultNbins;
    Int_t   fPtNbins;
    Int_t   fPtCorrNbins;
    Int_t   fPtCheckNbins;
    Int_t   fEtaNbins;
    Int_t   fEtaCheckNbins;
    Int_t   fZvNbins;
    Int_t   fCentralityNbins;
    Int_t   fPhiNbins;
    Double_t* fBinsMult; //[fMultNbins]
    Double_t* fBinsPt; //[fPtNbins]
    Double_t* fBinsPtCorr; //[fPtCorrNbins]
    Double_t* fBinsPtCheck; //[fPtCheckNbins]
    Double_t* fBinsEta; //[fEtaNbins]
    Double_t* fBinsEtaCheck; //[fEtaCheckNbins]
    Double_t* fBinsZv; //[fZvNbins]
    Double_t* fBinsCentrality; //[fCentralityNbins]
    Double_t* fBinsPhi; //[fPhiNbins]
    
    AlidNdPtAnalysisPbPbAOD(const AlidNdPtAnalysisPbPbAOD&); // not implemented
    AlidNdPtAnalysisPbPbAOD& operator=(const AlidNdPtAnalysisPbPbAOD&); // not implemented  
    
    ClassDef(AlidNdPtAnalysisPbPbAOD,5); // has to be at least 1, otherwise not streamable...
};

#endif
