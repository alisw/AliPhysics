#ifndef AlidNdPtAnalysisPbPbAOD_H
#define AlidNdPtAnalysisPbPbAOD_H


//------------------------------------------------------------------------------
// AlidNdPtAnalysisPbPbAOD class used for dNdPt analysis in PbPb collision
// via AODs 
// 
// Author: P. Luettig, 15.05.2013
//------------------------------------------------------------------------------

class iostream;

class TObject;
class TFile;
class TCint;
class THnSparse;

#include "AliAnalysisTaskSE.h"


#include "TList.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
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

#include "TSystem.h"
#include "TROOT.h"

class AlidNdPtAnalysisPbPbAOD : public AliAnalysisTaskSE {
  public :
    AlidNdPtAnalysisPbPbAOD(); 
    AlidNdPtAnalysisPbPbAOD(const char *name);
    ~AlidNdPtAnalysisPbPbAOD();
    
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);
    
    // Set binning for Histograms (if not set default binning is used)
    void SetBinsMult(Int_t nbins, Double_t* edges) 		{ Printf("[I] Setting Mult Bins"); fMultNbins = nbins; fBinsMult = GetArrayClone(nbins,edges); }
    void SetBinsPt(Int_t nbins, Double_t* edges) 		{ Printf("[I] Setting pT Bins"); fPtNbins = nbins; fBinsPt = GetArrayClone(nbins,edges); }
    void SetBinsPtCorr(Int_t nbins, Double_t* edges) 		{ Printf("[I] Setting pTcorr Bins"); fPtCorrNbins = nbins; fBinsPtCorr = GetArrayClone(nbins,edges); }
    void SetBinsEta(Int_t nbins, Double_t* edges) 		{ Printf("[I] Setting Eta Bins"); fEtaNbins = nbins; fBinsEta = GetArrayClone(nbins,edges); }
    void SetBinsZv(Int_t nbins, Double_t* edges) 		{ Printf("[I] Setting Zv Bins"); fZvNbins = nbins; fBinsZv= GetArrayClone(nbins,edges); }
    void SetBinsCentrality(Int_t nbins, Double_t* edges) 	{ Printf("[I] Setting Cent Bins"); fCentralityNbins = nbins; fBinsCentrality = GetArrayClone(nbins,edges); }
    
    
    // set event cut variables
    void SetCutMaxZVertex( Double_t d)					{ dCutMaxZVertex = d; }
    
    // set track kinematic cut parameters
    void SetCutPtRange(Double_t ptmin, Double_t ptmax)			{ dCutPtMin = ptmin; dCutPtMax = ptmax; }
    void SetCutEtaRange(Double_t etamin, Double_t etamax)		{ dCutEtaMin = etamin; dCutEtaMax = etamax; }
    
    // set track quality cut parameters
    void SetCutRequireTPCRefit(Bool_t *b) 				{ bCutRequireTPCRefit = b; } 
    void SetCutMinNCrossedRowsTPC(Double_t d) 				{ dCutMinNumberOfCrossedRows = d; }    
    void SetCutMinRatioCrossedRowsOverFindableClustersTPC(Double_t d) 	{ dCutMinRatioCrossedRowsOverFindableClustersTPC = d; }
    void SetCutMaxChi2PerClusterTPC(Double_t d) 			{ dCutMaxChi2PerClusterTPC = d; }
    void SetCutMaxFractionSharedTPCClusters(Double_t d) 		{ dCutMaxFractionSharedTPCClusters = d; }
    void SetCutMaxDCAToVertexZ(Double_t d) 				{ dCutMaxDCAToVertexZ = d; }
    void SetCutMaxDCAToVertexXY(Double_t d) 				{ dCutMaxDCAToVertexXY = d; }
    void SetCutRequireITSRefit(Bool_t *b) 				{ bCutRequireITSRefit = b; } 
    void SetCutMaxChi2PerClusterITS(Double_t d) 			{ dCutMaxChi2PerClusterITS = d; }
    void SetCutDCAToVertex2D(Bool_t *b) 				{ dCutDCAToVertex2D = b; } 
    void SetCutRequireSigmaToVertex(Bool_t *b) 				{ dCutRequireSigmaToVertex = b; } 
    void SetCutMaxDCAToVertexXYPtDep(Double_t d0, Double_t d1, Double_t d2)
    {
      dCutMaxDCAToVertexXYPtDepPar0 = d0;
      dCutMaxDCAToVertexXYPtDepPar1 = d1;
      dCutMaxDCAToVertexXYPtDepPar2 = d2;
    }
    void SetCutAcceptKinkDaughters(Bool_t *b) 				{ bCutAcceptKinkDaughters = b; } 
    void SetCutMaxChi2TPCConstrainedGlobal(Double_t d) 			{ dCutMaxChi2TPCConstrainedGlobal = d; }
    
    // getter for qualtiy track cuts
    Double_t GetCutMinNCrossedRowsTPC()					{ return dCutMinNumberOfCrossedRows; }
    
    // getter for DCA
    Bool_t GetDCA(const AliAODTrack *track, AliAODEvent *evt, Double_t d0z0[2]);
    
    THnSparseF *GetHistZvPtEtaCent() const { return hnZvPtEtaCent; }
    TH1F *GetHistEventStatistics() const { return hEventStatistics; }
    
    const char * GetParticleName(Int_t pdg);
    
    AliGenHijingEventHeader* GetHijingEventHeader(AliAODMCHeader *header);
    AliGenPythiaEventHeader* GetPythiaEventHeader(AliAODMCHeader *header);
    
    Bool_t IsTrackAccepted(AliAODTrack *tr);
    Bool_t IsMCTrackAccepted(AliAODMCParticle *part);
    
    Bool_t IsHijingParticle(const AliAODMCParticle *part, AliGenHijingEventHeader* hijingGenHeader);
    Bool_t IsPythiaParticle(const AliAODMCParticle *part, AliGenPythiaEventHeader* pythiaGenHeader);
    
    static Double_t* GetArrayClone(Int_t n, Double_t* source);
    
  private :
    
    // Output List
    TList	*fOutputList;
    
    // Histograms
    TH1F	*hPt; // simple pT histogramm
    TH1F	*hMCPt; // simple pT truth histogramm
    THnSparseF 	*hnZvPtEtaCent; //-> Zv:Pt:Eta:Cent
    THnSparseF 	*hnMCRecPrimZvPtEtaCent; //-> MC Zv:Pt:Eta:Cent
    THnSparseF 	*hnMCGenZvPtEtaCent; //-> MC Zv:Pt:Eta:Cent
    THnSparseF 	*hnMCRecSecZvPtEtaCent; //-> MC Zv:Pt:Eta:Cent, only secondaries
    TH1F	*hEventStatistics; // contains statistics of number of events after each cut
    TH1F	*hEventStatisticsCentrality; // contains number of events vs centrality, events need to have a track in kinematic range
    TH1F	*hMCEventStatisticsCentrality; // contains MC number of events vs centrality, events need to have a track in kinematic range
    TH1F	*hAllEventStatisticsCentrality; // contains number of events vs centrality, events need to be triggered
    TH2F	*hEventStatisticsCentralityTrigger; // contains number of events vs centrality in 1% bins vs trigger
    THnSparseF	*hnZvMultCent; // Zv:Mult:Cent
    TH1F	*hTriggerStatistics; // contains number of events per trigger
    TH1F	*hMCTrackPdgCode; // contains statistics of pdg codes of tracks
    TH1F	*hMCTrackStatusCode; // contains statistics of status codes of tracks
    TH1F	*hCharge; // charge distribution in data
    TH1F	*hMCCharge; // charge distribution in MC
    TH2F	*hMCPdgPt; // PDGvs PT for MC Particles
    TH1F	*hMCHijingPrim; // number of particles, which are Hijing particles and primaries
    THnSparseF	*hDCAPtAll; //control histo: DCAz:DCAxy:pT:eta:phi for all reconstructed tracks
    THnSparseF	*hDCAPtAccepted; //control histo: DCAz:DCAxy:pT:eta:phi for all accepted reco tracks
    THnSparseF	*hMCDCAPtSecondary; //control histo: DCAz:DCAxy:pT:eta:phi for all accepted reco track, which are secondaries (using MC info)
    THnSparseF	*hMCDCAPtPrimary; //control histo: DCAz:DCAxy:pT:eta:phi for all accepted reco track, which are primaries (using MC info)
    THnSparseF	*hnCrossedRowsClustersChiPtEtaPhiAll; // CrossedRows:Cluster:ChiperCluster:pT:Eta:Phi before track cuts
    THnSparseF	*hnCrossedRowsClustersChiPtEtaPhiAcc; // CrossedRows:Cluster:ChiperCluster:pT:Eta:Phi after track cuts
   
    
    // global variables
    Bool_t bIsMonteCarlo;
     
    // event cut variables
    Double_t dCutMaxZVertex;
    
    // track kinematic cut variables
    Double_t dCutPtMin;
    Double_t dCutPtMax;
    Double_t dCutEtaMin;
    Double_t dCutEtaMax;
    
    // track quality cut variables
    Bool_t 	bCutRequireTPCRefit;
    Double_t 	dCutMinNumberOfCrossedRows;
    Double_t 	dCutMinRatioCrossedRowsOverFindableClustersTPC;
    Double_t 	dCutMaxChi2PerClusterTPC;
    Double_t 	dCutMaxFractionSharedTPCClusters;
    Double_t 	dCutMaxDCAToVertexZ;
    Double_t 	dCutMaxDCAToVertexXY;
    Bool_t 	bCutRequireITSRefit;
    Double_t 	dCutMaxChi2PerClusterITS;
    Bool_t 	dCutDCAToVertex2D;
    Bool_t 	dCutRequireSigmaToVertex;
    Double_t 	dCutMaxDCAToVertexXYPtDepPar0;
    Double_t 	dCutMaxDCAToVertexXYPtDepPar1;
    Double_t 	dCutMaxDCAToVertexXYPtDepPar2;
    Bool_t 	bCutAcceptKinkDaughters;
    Double_t 	dCutMaxChi2TPCConstrainedGlobal;
    
    //binning for THNsparse
    Int_t fMultNbins;
    Int_t fPtNbins;
    Int_t fPtCorrNbins;
    Int_t fEtaNbins;
    Int_t fZvNbins;
    Int_t fCentralityNbins;
    Double_t* fBinsMult; //[fMultNbins]
    Double_t* fBinsPt; //[fPtNbins]
    Double_t* fBinsPtCorr; //[fPtCorrNbins]
    Double_t* fBinsEta; //[fEtaNbins]
    Double_t* fBinsZv; //[fZvNbins]
    Double_t* fBinsCentrality; //[fCentralityNbins]
    
    AlidNdPtAnalysisPbPbAOD(const AlidNdPtAnalysisPbPbAOD&); // not implemented
    AlidNdPtAnalysisPbPbAOD& operator=(const AlidNdPtAnalysisPbPbAOD&); // not implemented  
    
    ClassDef(AlidNdPtAnalysisPbPbAOD,3); // has to be at least 1, otherwise not streamable...
};

#endif
