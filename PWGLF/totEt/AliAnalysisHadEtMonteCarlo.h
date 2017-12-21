//_________________________________________________________________________
//  Utility Class for transverse energy studies, charged hadrons
//  Base class for MC analysis
//  - MC output
// 
//Created by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville 
//_________________________________________________________________________
#ifndef ALIANALYSISHADETMONTECARLO_H
#define ALIANALYSISHADETMONTECARLO_H

#include "AliAnalysisHadEt.h"
class AliVEvent;
class TRandom;
class AliAnalysisHadEtReconstructed;
class AliMCEvent;
class AliGenEventHeader;
class AliStack;

class AliAnalysisHadEtMonteCarlo : public AliAnalysisHadEt
{

public:
   
  AliAnalysisHadEtMonteCarlo();
  virtual ~AliAnalysisHadEtMonteCarlo();
   
    virtual Int_t AnalyseEvent(AliVEvent* event);
    virtual Int_t AnalyseEvent(AliVEvent* event,AliVEvent* event2);

    //void FillHistograms();
    void CreateHistograms();
    virtual void Init();

    Float_t GetSimulatedHadronicEt() const {return fSimHadEt;}
    Float_t GetSimulatedTotalEt() const {return fSimTotEt;}

    void FillSimTotEtVsRecoTotEtFullAcceptanceTPC(Float_t et){if(fGoodEvent)FillHisto2D("SimTotEtVsRecoTotEtFullAcceptanceTPC",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtFullAcceptanceITS(Float_t et){if(fGoodEvent)FillHisto2D("SimTotEtVsRecoTotEtFullAcceptanceITS",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtEMCALAcceptanceTPC(Float_t et){if(fGoodEvent)FillHisto2D("SimTotEtVsRecoTotEtEMCALAcceptanceTPC",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtEMCALAcceptanceITS(Float_t et){if(fGoodEvent)FillHisto2D("SimTotEtVsRecoTotEtEMCALAcceptanceITS",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtPHOSAcceptanceTPC(Float_t et){if(fGoodEvent)FillHisto2D("SimTotEtVsRecoTotEtPHOSAcceptanceTPC",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtPHOSAcceptanceITS(Float_t et){if(fGoodEvent)FillHisto2D("SimTotEtVsRecoTotEtPHOSAcceptanceITS",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtFullAcceptanceTPCNoPID(Float_t et){if(fGoodEvent)FillHisto2D("SimTotEtVsRecoTotEtFullAcceptanceTPCNoPID",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtFullAcceptanceITSNoPID(Float_t et){if(fGoodEvent)FillHisto2D("SimTotEtVsRecoTotEtFullAcceptanceITSNoPID",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtEMCALAcceptanceTPCNoPID(Float_t et){if(fGoodEvent)FillHisto2D("SimTotEtVsRecoTotEtEMCALAcceptanceTPCNoPID",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtEMCALAcceptanceITSNoPID(Float_t et){if(fGoodEvent)FillHisto2D("SimTotEtVsRecoTotEtEMCALAcceptanceITSNoPID",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtPHOSAcceptanceTPCNoPID(Float_t et){if(fGoodEvent)FillHisto2D("SimTotEtVsRecoTotEtPHOSAcceptanceTPCNoPID",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtPHOSAcceptanceITSNoPID(Float_t et){if(fGoodEvent)FillHisto2D("SimTotEtVsRecoTotEtPHOSAcceptanceITSNoPID",fSimTotEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtFullAcceptanceTPC(Float_t et){if(fGoodEvent)FillHisto2D("SimHadEtVsRecoHadEtFullAcceptanceTPC",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtFullAcceptanceITS(Float_t et){if(fGoodEvent)FillHisto2D("SimHadEtVsRecoHadEtFullAcceptanceITS",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtEMCALAcceptanceTPC(Float_t et){if(fGoodEvent)FillHisto2D("SimHadEtVsRecoHadEtEMCALAcceptanceTPC",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtEMCALAcceptanceITS(Float_t et){if(fGoodEvent)FillHisto2D("SimHadEtVsRecoHadEtEMCALAcceptanceITS",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtPHOSAcceptanceTPC(Float_t et){if(fGoodEvent)FillHisto2D("SimHadEtVsRecoHadEtPHOSAcceptanceTPC",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtPHOSAcceptanceITS(Float_t et){if(fGoodEvent)FillHisto2D("SimHadEtVsRecoHadEtPHOSAcceptanceITS",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtFullAcceptanceTPCNoPID(Float_t et){if(fGoodEvent)FillHisto2D("SimHadEtVsRecoHadEtFullAcceptanceTPCNoPID",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtFullAcceptanceITSNoPID(Float_t et){if(fGoodEvent)FillHisto2D("SimHadEtVsRecoHadEtFullAcceptanceITSNoPID",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtEMCALAcceptanceTPCNoPID(Float_t et){if(fGoodEvent)FillHisto2D("SimHadEtVsRecoHadEtEMCALAcceptanceTPCNoPID",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtEMCALAcceptanceITSNoPID(Float_t et){if(fGoodEvent)FillHisto2D("SimHadEtVsRecoHadEtEMCALAcceptanceITSNoPID",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtPHOSAcceptanceTPCNoPID(Float_t et){if(fGoodEvent)FillHisto2D("SimHadEtVsRecoHadEtPHOSAcceptanceTPCNoPID",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtPHOSAcceptanceITSNoPID(Float_t et){if(fGoodEvent)FillHisto2D("SimHadEtVsRecoHadEtPHOSAcceptanceITSNoPID",fSimHadEt,et,1.0);}
    //PiKP histograms
    void FillSimPiKPEtVsRecoPiKPEtFullAcceptanceTPC(Float_t et){if(fGoodEvent)FillHisto2D("SimPiKPEtVsRecoPiKPEtFullAcceptanceTPC",fSimPiKPEt,et,1.0);}
    void FillSimPiKPEtVsRecoPiKPEtFullAcceptanceITS(Float_t et){if(fGoodEvent)FillHisto2D("SimPiKPEtVsRecoPiKPEtFullAcceptanceITS",fSimPiKPEt,et,1.0);}
    void FillSimPiKPEtVsRecoPiKPEtFullAcceptanceTPCNoPID(Float_t et){if(fGoodEvent)FillHisto2D("SimPiKPEtVsRecoPiKPEtFullAcceptanceTPCNoPID",fSimPiKPEt,et,1.0);}
    void FillSimPiKPEtVsRecoPiKPEtFullAcceptanceITSNoPID(Float_t et){if(fGoodEvent)FillHisto2D("SimPiKPEtVsRecoPiKPEtFullAcceptanceITSNoPID",fSimPiKPEt,et,1.0);}
    //Raw histograms
    void FillSimRawEtVsRecoRawEtFullAcceptanceTPC(Float_t et){if(fGoodEvent)FillHisto2D("SimRawEtVsRecoRawEtFullAcceptanceTPC",fSimRawEtTPC,et,1.0);}
    void FillSimRawEtVsRecoRawEtFullAcceptanceITS(Float_t et){if(fGoodEvent)FillHisto2D("SimRawEtVsRecoRawEtFullAcceptanceITS",fSimRawEtITS,et,1.0);}
    void FillSimRawEtVsRecoRawEtFullAcceptanceTPCNoPID(Float_t et){if(fGoodEvent)FillHisto2D("SimRawEtVsRecoRawEtFullAcceptanceTPCNoPID",fSimRawEtTPC,et,1.0);}
    void FillSimRawEtVsRecoRawEtFullAcceptanceITSNoPID(Float_t et){if(fGoodEvent)FillHisto2D("SimRawEtVsRecoRawEtFullAcceptanceITSNoPID",fSimRawEtITS,et,1.0);}

    void FillSimTotEtMinusRecoTotEtFullAcceptanceTPC(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRecoTotEtFullAcceptanceTPC",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimTotEtMinusRecoTotEtFullAcceptanceITS(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRecoTotEtFullAcceptanceITS",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimTotEtMinusRecoTotEtEMCALAcceptanceTPC(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRecoTotEtEMCALAcceptanceTPC",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimTotEtMinusRecoTotEtEMCALAcceptanceITS(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRecoTotEtEMCALAcceptanceITS",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimTotEtMinusRecoTotEtPHOSAcceptanceTPC(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRecoTotEtPHOSAcceptanceTPC",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimTotEtMinusRecoTotEtPHOSAcceptanceITS(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRecoTotEtPHOSAcceptanceITS",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimTotEtMinusRecoTotEtFullAcceptanceTPCNoPID(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRecoTotEtFullAcceptanceTPCNoPID",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimTotEtMinusRecoTotEtFullAcceptanceITSNoPID(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRecoTotEtFullAcceptanceITSNoPID",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimTotEtMinusRecoTotEtEMCALAcceptanceTPCNoPID(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRecoTotEtEMCALAcceptanceTPCNoPID",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimTotEtMinusRecoTotEtEMCALAcceptanceITSNoPID(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRecoTotEtEMCALAcceptanceITSNoPID",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimTotEtMinusRecoTotEtPHOSAcceptanceTPCNoPID(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRecoTotEtPHOSAcceptanceTPCNoPID",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimTotEtMinusRecoTotEtPHOSAcceptanceITSNoPID(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRecoTotEtPHOSAcceptanceITSNoPID",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimHadEtMinusRecoHadEtFullAcceptanceTPC(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRecoHadEtFullAcceptanceTPC",et,(fSimHadEt-et)/fSimHadEt,1.0);}
    void FillSimHadEtMinusRecoHadEtFullAcceptanceITS(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRecoHadEtFullAcceptanceITS",et,(fSimHadEt-et)/fSimHadEt,1.0);}
    void FillSimHadEtMinusRecoHadEtEMCALAcceptanceTPC(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRecoHadEtEMCALAcceptanceTPC",et,(fSimHadEt-et)/fSimHadEt,1.0);}
    void FillSimHadEtMinusRecoHadEtEMCALAcceptanceITS(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRecoHadEtEMCALAcceptanceITS",et,(fSimHadEt-et)/fSimHadEt,1.0);}
    void FillSimHadEtMinusRecoHadEtPHOSAcceptanceTPC(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRecoHadEtPHOSAcceptanceTPC",et,(fSimHadEt-et)/fSimHadEt,1.0);}
    void FillSimHadEtMinusRecoHadEtPHOSAcceptanceITS(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRecoHadEtPHOSAcceptanceITS",et,(fSimHadEt-et)/fSimHadEt,1.0);}
    void FillSimHadEtMinusRecoHadEtFullAcceptanceTPCNoPID(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRecoHadEtFullAcceptanceTPCNoPID",et,(fSimHadEt-et)/fSimHadEt,1.0);}
    void FillSimHadEtMinusRecoHadEtFullAcceptanceITSNoPID(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRecoHadEtFullAcceptanceITSNoPID",et,(fSimHadEt-et)/fSimHadEt,1.0);}
    void FillSimHadEtMinusRecoHadEtEMCALAcceptanceTPCNoPID(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRecoHadEtEMCALAcceptanceTPCNoPID",et,(fSimHadEt-et)/fSimHadEt,1.0);}
    void FillSimHadEtMinusRecoHadEtEMCALAcceptanceITSNoPID(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRecoHadEtEMCALAcceptanceITSNoPID",et,(fSimHadEt-et)/fSimHadEt,1.0);}
    void FillSimHadEtMinusRecoHadEtPHOSAcceptanceTPCNoPID(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRecoHadEtPHOSAcceptanceTPCNoPID",et,(fSimHadEt-et)/fSimHadEt,1.0);}
    void FillSimHadEtMinusRecoHadEtPHOSAcceptanceITSNoPID(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRecoHadEtPHOSAcceptanceITSNoPID",et,(fSimHadEt-et)/fSimHadEt,1.0);}


    void FillSimTotEtMinusRawEtFullAcceptanceTPC(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRawEtFullAcceptanceTPC",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimTotEtMinusRawEtFullAcceptanceITS(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRawEtFullAcceptanceITS",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimTotEtMinusRawEtEMCALAcceptanceTPC(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRawEtEMCALAcceptanceTPC",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimTotEtMinusRawEtEMCALAcceptanceITS(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRawEtEMCALAcceptanceITS",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimTotEtMinusRawEtPHOSAcceptanceTPC(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRawEtPHOSAcceptanceTPC",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimTotEtMinusRawEtPHOSAcceptanceITS(Float_t et){if(fGoodEvent && fSimTotEt>0.0&&et>0.0)FillHisto2D("SimTotEtMinusRawEtPHOSAcceptanceITS",et,(fSimTotEt-et)/fSimTotEt,1.0);}
    void FillSimHadEtMinusRawEtFullAcceptanceTPC(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRawEtFullAcceptanceTPC",et,(fSimHadEt-et)/fSimHadEt,1.0);}
    void FillSimHadEtMinusRawEtFullAcceptanceITS(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRawEtFullAcceptanceITS",et,(fSimHadEt-et)/fSimHadEt,1.0);}
    void FillSimHadEtMinusRawEtEMCALAcceptanceTPC(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRawEtEMCALAcceptanceTPC",et,(fSimHadEt-et)/fSimHadEt,1.0);}
    void FillSimHadEtMinusRawEtEMCALAcceptanceITS(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRawEtEMCALAcceptanceITS",et,(fSimHadEt-et)/fSimHadEt,1.0);}
    void FillSimHadEtMinusRawEtPHOSAcceptanceTPC(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRawEtPHOSAcceptanceTPC",et,(fSimHadEt-et)/fSimHadEt,1.0);}
    void FillSimHadEtMinusRawEtPHOSAcceptanceITS(Float_t et){if(fGoodEvent && fSimHadEt>0.0&&et>0.0)FillHisto2D("SimHadEtMinusRawEtPHOSAcceptanceITS",et,(fSimHadEt-et)/fSimHadEt,1.0);}

    void FillSimPiKPMinusRecoPiKPFullAcceptanceTPC(Float_t et){if(fGoodEvent && fSimPiKPEt>0.0)FillHisto2D("SimPiKPMinusRecoPiKPFullAcceptanceTPC",et,(fSimPiKPEt-et)/fSimPiKPEt,1.0);}
    void FillSimPiKPMinusRecoPiKPFullAcceptanceITS(Float_t et){if(fGoodEvent && fSimPiKPEt>0.0)FillHisto2D("SimPiKPMinusRecoPiKPFullAcceptanceITS",et,(fSimPiKPEt-et)/fSimPiKPEt,1.0);}
    void FillSimPiKPMinusRecoPiKPFullAcceptanceTPCNoPID(Float_t et){if(fGoodEvent && fSimPiKPEt>0.0)FillHisto2D("SimPiKPMinusRecoPiKPFullAcceptanceTPCNoPID",et,(fSimPiKPEt-et)/fSimPiKPEt,1.0);}
    void FillSimPiKPMinusRecoPiKPFullAcceptanceITSNoPID(Float_t et){if(fGoodEvent && fSimPiKPEt>0.0)FillHisto2D("SimPiKPMinusRecoPiKPFullAcceptanceITSNoPID",et,(fSimPiKPEt-et)/fSimPiKPEt,1.0);}

    void InvestigateSmearing(Bool_t val){fInvestigateSmearing=val;}
    void RunLightweight(Bool_t val){fRunLightweight=val;}//Switches off O(100) histograms that aren't regularly used.
    void InvestigateFull(Bool_t val){fInvestigateFull=val;}
    void InvestigateEMCAL(Bool_t val){fInvestigateEMCal=val;}
    void InvestigatePHOS(Bool_t val){fInvestigatePHOS=val;}
    void InvestigatePiKP(Bool_t val){fInvestigatePiKP=val;}
    void RequireITSHits(Bool_t val){fRequireITSHits=val;}
    void EnhanceBaryons(Bool_t val){fBaryonEnhancement=val;}
    void UseReconstructedPt(Bool_t val){fUseRecoPt=val;}
    Bool_t Full() const {return fInvestigateFull;}
    Bool_t EMCAL() const {return fInvestigateEMCal;}
    Bool_t PHOS()const {return fInvestigatePHOS;}
    Bool_t PiKP() const {return fInvestigatePiKP;}
    Bool_t BaryonEnhancement() const {return fBaryonEnhancement;}
    Bool_t ReconstructedPt() const {return fUseRecoPt;}
    void SetHadEtReconstructed(AliAnalysisHadEtReconstructed *reco){fHadEtReco = reco;}
    void DoTriggerChecks(){kDoTriggerChecks=kTRUE;}
    void DoTriggerChecksOnly(){kDoTriggerChecksOnly=kTRUE;}
    void SetIsOfflineV0AND(Bool_t val){kIsOfflineV0AND = val;}
    void SetIsOfflineMB(Bool_t val){kIsOfflineMB = val;}
    void SetUseOldCentrality(Bool_t val){useOldCentrality = val;}


    void SetGeneratorMinMaxParticles(AliMCEvent *eventMC);
    AliGenEventHeader* GetGenEventHeader(AliMCEvent *eventMC) const;
    Bool_t IsHIJINGLabel(Int_t label,AliMCEvent *eventMC,AliStack *stack);
    void CheckForHIJINGLabel(){checkLabelForHIJING = kTRUE;}

 protected:
    Bool_t checkLabelForHIJING;//Boolean - do we need to check this data set for HIJING labels?  We don't need to do this unless there were signals embedded
    Int_t            fNMCProducedMin;             // In case of cocktail, select particles in the list with label from this value
    Int_t            fNMCProducedMax;             // In case of cocktail, select particles in the list with label up to this value
    Bool_t useOldCentrality;//Use AliCentrality or AliMultSelection

 private:
    //Declare it private to avoid compilation warning
    AliAnalysisHadEtMonteCarlo & operator = (const AliAnalysisHadEtMonteCarlo & g) ;//cpy assignment
    AliAnalysisHadEtMonteCarlo(const AliAnalysisHadEtMonteCarlo & g) ; // cpy ctor

    Float_t fSimPiKPEt;//simulated Et for pi,k,p event by event
    Float_t fSimRawEtTPC;//simulated Et for pi,k,p event by event
    Float_t fSimRawEtITS;//simulated Et for pi,k,p event by event
    Float_t fSimHadEt;//simulated Et event by event
    Float_t fSimTotEt;//total et event by event
    Float_t fSimPiKPEtShouldBeReco;//simulated Et for pi,k,p event by event
    Float_t fSimPiKPEtShouldBeRecoPi;//simulated Et for pi event by event
    Float_t fSimPiKPEtShouldBeRecoK;//simulated Et for k event by event
    Float_t fSimPiKPEtShouldBeRecoP;//simulated Et for p event by event

    Bool_t fRunLightweight;//Turns off O(100) histograms that are not regularly used to reduce the size & complexity of the output
    Bool_t fInvestigateSmearing;//Turns on and off functions and histos for investigating momentum, efficiency, pid smearing
    Bool_t fInvestigateFull;//Turns on and off functions and histos for investigating event-by-event et for the full acceptance
    Bool_t fInvestigateEMCal;//Turns on and off functions and histos for investigating event-by-event et for the full acceptance
    Bool_t fInvestigatePHOS;//Turns on and off functions and histos for investigating event-by-event et for the full acceptance
    Bool_t fInvestigatePiKP;//Turns on and off functions and histos for looking pi/k/p Et event-by-event
    Bool_t fRequireITSHits;//Also investigates Et for track cuts with ITS+TPC hits
    Bool_t fBaryonEnhancement;//Turns on and off baryon enhancement
    Bool_t fUseRecoPt;//Switch to decide whether to use reconstructed or simulated pT for reconstructed tracks
    Bool_t kIsOfflineV0AND;//Boolean to keep track of whether or not this matches the offline trigger
    Bool_t kIsOfflineMB;//Boolean to keep track of whether or not this matches the offline trigger
    Bool_t kDoTriggerChecks;//Boolean to keep track of whether or not I want to let the physics selection do the work
    Bool_t kDoTriggerChecksOnly;//Boolean to keep track of whether or not I want to let the physics selection do the work

    void ResetEventValues();

    //Float_t fSimPiKPEtSmeared[4];//simulated Et for pi,k,p smeared for each event by different momentum resolutions
    static Float_t fgSmearWidths[4];//array with widths for smearing with different momentum resultions
    static Int_t fgNumSmearWidths;//number of entries in the array above
    TRandom *fPtSmearer;//a TRandom used for investigating momentum smearing
    AliAnalysisHadEtReconstructed *fHadEtReco;//corresponding class for data reconstruction.  Used to access corrections directly used here.

    ClassDef(AliAnalysisHadEtMonteCarlo, 1);
};

#endif // ALIANALYSISHADETMONTECARLO_H
