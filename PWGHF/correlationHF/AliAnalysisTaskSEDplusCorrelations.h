#ifndef ALIANALYSISTASKSEDPLUSCORRELATIONS_H
#define ALIANALYSISTASKSEDPLUSCORRELATIONS_H


/* Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisTaskSEDplusCorrelations.h 58883 2012-10-02 09:41:01Z prino $ */

//*************************************************************************
// Class AliAnalysisTaskSEDplusCorrelations
// AliAnalysisTaskSE for Dplus candidates (3Prongs) and hadrons correlations
// Authors: Jitendra


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
    void SetCorrFormPart(Bool_t genMC){fMontecarlo=genMC;}
    void SetCorrFormTrack(Bool_t reco){fReco=reco;}
    void SetDataOrMC(Bool_t readMC){fReadMC=readMC;}
    void SetEventMixing(Bool_t mixing){fMixing=mixing;}
    void SetCorrelator(Int_t number) {fSelect = number;} // select 1 for hadrons, 2 for Kaons, 3 for Kzeros
    void SetSystem(Bool_t system){fSystem=system;} // select between pp (kFALSE) or PbPb (kTRUE)
    void SetEtaRagne(Double_t etacorr) {fEtaRange=etacorr;}
    void SetUseBit(Bool_t bits=kTRUE){fUseBit=bits;}
    void SetTCConfig(Bool_t TCcong=kFALSE){fTCconfig=TCcong;}
    void SetTrackEffActive(Bool_t effTrack=kFALSE){fEffTrack=effTrack;}
    void SetDplusEffActive(Bool_t effDplus=kFALSE){fEffDplus=effDplus;}
    void SetUseCentrality(Bool_t flag, Int_t estimator){fEvalCentrality=flag; fCentralityEstimator=estimator;}
    void SetBinWidth(Float_t BinW){fBinWidth=BinW;}
    void SetMCGevEventType(Bool_t sel1=kFALSE){fMCGenEvType=sel1;}
    void SetPoolByPoolCorr(Bool_t sel2=kFALSE){fPoolByPool=sel2;}
    void SetCheckCutDist(Bool_t sel3=kFALSE){fCheckCutDist=sel3;}
    // void SetUseDisplacement(Int_t m) {fDisplacement=m;} // select 0 for no displ, 1 for abs displ, 2 for d0/sigma_d0
    
    
    private :
    
    AliAnalysisTaskSEDplusCorrelations(const AliAnalysisTaskSEDplusCorrelations &source);
    AliAnalysisTaskSEDplusCorrelations& operator=(const AliAnalysisTaskSEDplusCorrelations& source);
    
    //correlation methods
    void HistoNomenclature();
    void HadronCorrelations(AliAODRecoDecayHF3Prong* d,TClonesArray *arrayMC, Bool_t isDplus);
    void CorrelationNSparsePlots(AliAODRecoDecayHF3Prong *d, AliReducedParticle* track, Int_t iPtBin, Bool_t *origDplus, Double_t weightEff);
    
    
    Int_t fSelect; // Correlation Option between D+ and (1-chargedtracks,2-chargedkaons,3-k0s )
    TList *fOutput;                  //! user output data
    TList *fOutputCorr;                  //! user output data
    Bool_t fReadMC; //  MC Switch
    Bool_t fReco; // Switch to reco track
    Bool_t fMontecarlo; // Switch to Montecarlo Gen level
    Bool_t fMCGenEvType; //Gen MC event type
    Bool_t fMixing;// switch for event mixing
    TClonesArray* farrayMC; //! mcarray
    Bool_t fSystem; // pp or PbPb
    Bool_t fUseBit; //  filterbit option
    Bool_t fTCconfig; //  TC Cuts option
    TH1F *fHistNEvents; //!hist. for No. of events
    TH1F *fHistNDplus; //!hist. for No. of Dplus
    AliNormalizationCounter *fCounter; // counter
    AliRDHFCutsDplustoKpipi *fDplusCuts;  // Cuts D+
    AliHFAssociatedTrackCuts *fAssoCuts; // cuts for associated track
    AliHFCorrelator  *fCorrelator; //object for correlations
    Double_t  fEtaRange;		// cut for Dplus eta to
    Int_t fNPtBins; // number of event at different Stages
    Float_t fBinWidth;//width of one bin in output histos
    Double_t fCentrOrMult; // Multiplicity of Event for D eff
    Double_t fMultiplicity; //Multiplicity for maps
    Bool_t fEffTrack; //Track eff ON/OFF
    Bool_t fEffDplus; //Dplus eff ON/OFF
    
    Int_t  fCentralityEstimator;   // enum from AliRDHFCuts..
    Bool_t    fEvalCentrality; // Switch to ON/OFF the centrality interface
    Double_t  fMinCentrality; // Minimun Centrality Value
    Double_t  fMaxCentrality; // Maximum Centrality Value
    Bool_t  fPoolByPool;
    Int_t  fWhichPool;
    Bool_t fCheckCutDist; //flag to check topological cuts distribuition
    ClassDef(AliAnalysisTaskSEDplusCorrelations,5); // class for D+ meson correlations
    
};

#endif
