#ifndef ALIANALYSISTASKELECHADRONCORREL_H
#define ALIANALYSISTASKELECHADRONCORREL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  Task for Heavy Flavour Electron-Hadron DeltaPhi Correlation       //
//                                                                    //
//  Author: Deepa Thomas (Utrecht University)                         //
//                                                                    //
////////////////////////////////////////////////////////////////////////

class THnSparse;
class TH2F;
class TLorentzVector;

class AliEMCALTrack;
class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliAnalysisFilter;
class AliESDtrackCuts;
class AliESDtrack;
class AliAODTrack;
class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliCFManager;
//class AliEventPoolManager;

#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"

class AliAnalysisTaskElecHadronCorrel : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskElecHadronCorrel();
    AliAnalysisTaskElecHadronCorrel(const char *name);
    virtual ~AliAnalysisTaskElecHadronCorrel();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };
    Bool_t IsESDanalysis() const { return !TestBit(kAODanalysis); };

    void SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };
    void SetInvariantMassCut (Double_t invmass) {fInvmassCut = invmass;};
    AliHFEpid *GetPID() const { return fPID; }
    void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) { fRejectKinkMother = rejectKinkMother; };
    void SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec);
    void ElectronHadCorrel(Int_t itrack, AliVTrack *track, TH2F *DphiPt, TH2F *DphiPt1,TH2F *DphiPt2,TH2F *DphiPt3,TH2F *DphiPt4);	
    void ElectronHadCorrelEtaFarSide(Int_t itrack, AliVTrack *track, TH2F *DphiPt, TH2F *DphiPt1,TH2F *DphiPt2,TH2F *DphiPt3,TH2F *DphiPt4);	
    void ElectronHadCorrelEtaBins(Int_t itrack, AliVTrack *track, TH2F *DphiPtEta1, TH2F *DphiPtEta11,TH2F *DphiPtEta12,TH2F *DphiPtEta13,TH2F *DphiPtEta14,TH2F *DphiPtEta2, TH2F *DphiPtEta21,TH2F *DphiPtEta22,TH2F *DphiPtEta23,TH2F *DphiPtEta24);	
   // void ElectronHadCorrelEtaBins(Int_t itrack, AliVTrack *track, TH3F *DphiPtEta1, TH3F *DphiPtEta11,TH3F *DphiPtEta12,TH3F *DphiPtEta13,TH3F *DphiPtEta14);	
    void ElectronHadCorrelNoPartner(Int_t itrack,Int_t jtrack, AliVTrack *track, TH2F *DphiPtNew,TH2F *DphiPtNew1,TH2F *DphiPtNew2,TH2F *DphiPtNew3,TH2F *DphiPtNew4);	
    void ElectronHadCorrelEtaBinsNoPartner(Int_t itrack,Int_t jtrack, AliVTrack *track, TH2F *DphiPtEta1, TH2F *DphiPtEta11,TH2F *DphiPtEta12,TH2F *DphiPtEta13,TH2F *DphiPtEta14,TH2F *DphiPtEta2, TH2F *DphiPtEta21,TH2F *DphiPtEta22,TH2F *DphiPtEta23,TH2F *DphiPtEta24);	
   // void ElectronHadCorrelEtaBinsNoPartner(Int_t itrack,Int_t jtrack, AliVTrack *track, TH3F *DphiPtEta1, TH3F *DphiPtEta11,TH3F *DphiPtEta12,TH3F *DphiPtEta13,TH3F *DphiPtEta14);	
    void HadronInfo(Int_t itrack);
    void    SetCentralityParameters(Double_t CentralityMin, Double_t CentralityMax, const char* CentralityMethod); //select centrality
    void    CheckCentrality(AliVEvent *event,Bool_t &centralitypass); //to use only events with the correct centrality....

    void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
    void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };
    /*    void MixedEvent(AliAODTrack *track, TH2F *DphiPt, TH2F *DphiPt1, TH2F *DphiPt2);
    TObjArray* CloneAndReduceTrackList();
  */
    private:

    enum{
      kAODanalysis = BIT(20),
    };

    Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);

    AliVEvent 		      *fVevent;		//V event object
    AliESDEvent 		   *fESD;			//ESD object
    AliAODEvent 		   *fAOD;			//AOD object
    AliEMCALGeometry  	*fGeom; 		// emcal geometry 

    TList       	   	*fOutputList;		//output list

    AliESDtrackCuts		*fTrackCuts1;		//ESD track cuts
    AliESDtrackCuts		*fTrackCuts2;		//ESD track cuts
    AliHFEcuts     		*fCuts;                 //Cut Collection
    Bool_t 			      fIdentifiedAsOutInz;    //Out Of Range in z
    Bool_t 	   	   	fPassTheEventCut;       //Pass The Event Cut
    Bool_t 		   	   fRejectKinkMother;      //Reject Kink Mother
    Double_t 		      fVz;                    //z position of the primary vertex
    AliCFManager 	   	*fCFM;                  //!Correction Framework Manager
    AliHFEpid 		      *fPID;                  //PID
    AliHFEpidQAmanager 	*fPIDqa;		//! PID QA manager
    Double_t	         	fInvmassCut;		//invariant mass cut value
    Double_t             fCentrality; // event centrality for QA
    Double_t             fCentralityMin; // lower bound of cenrality bin
    Double_t             fCentralityMax; // upper bound of centrality bin
    const char           *fkCentralityMethod; // method used to determine centrality (V0 by default)   
//    AliEventPoolManager*     fPoolMgr;         //! event pool manager

    TH1F			*fNoEvents;		//no of events
//    TH1F			*fTrkpt;		//track pt
    TH2F			*fTrkEovPAft;		//track E/p after HFE pid
 //   TH2F			*fTrkEovPBefHad;		//track E/p before HFE pid
//    TH2F			*fdEdxBef;		//track dEdx vs p before HFE pid
    TH2F			*fSemiIncElecDphi;  	//Semi Inclusive elec - had DPhi
    TH2F			*fSemiIncElecDphi1;  	//Semi Inclusive elec - had DPhi
    TH2F			*fSemiIncElecDphi2;  	//Semi Inclusive elec - had DPhi
    TH2F			*fSemiIncElecDphi3;  	//Semi Inclusive elec - had DPhi
    TH2F			*fSemiIncElecDphi4;  	//Semi Inclusive elec - had DPhi
    TH2F			*fPhotElecDphi;  	//Photon elec - had DPhi
    TH2F			*fPhotElecDphi1;  	//Photon elec - had DPhi
    TH2F			*fPhotElecDphi2;  	//Photon elec - had DPhi
    TH2F			*fPhotElecDphi3;  	//Photon elec - had DPhi
    TH2F			*fPhotElecDphi4;  	//Photon elec - had DPhi
    TH2F			*fInclusiveElecDphi;  	//Inclusive elec - had DPhi
    TH2F			*fInclusiveElecDphi1;  	//Inclusive elec - had DPhi
    TH2F			*fInclusiveElecDphi2;  	//Inclusive elec - had DPhi
    TH2F			*fInclusiveElecDphi3;  	//Inclusive elec - had DPhi
    TH2F			*fInclusiveElecDphi4;  	//Inclusive elec - had DPhi
    TH2F       *fInclusiveElecDphiEtaFS;    //Inclusive elec EtaFS- had DPhi
    TH2F       *fInclusiveElecDphiEtaFS1;   //Inclusive elec EtaFS- had DPhi
    TH2F       *fInclusiveElecDphiEtaFS2;   //Inclusive elec EtaFS- had DPhi
    TH2F       *fInclusiveElecDphiEtaFS3;   //Inclusive elec EtaFS- had DPhi
    TH2F       *fInclusiveElecDphiEtaFS4;   //Inclusive elec EtaFS- had DPhi
    TH2F			*fDphiULSMassLow;	//Dphi - ULS, mass< mass cut
    TH2F			*fDphiULSMassLow1;	//Dphi - ULS, mass< mass cut
    TH2F			*fDphiULSMassLow2;	//Dphi - ULS, mass< mass cut
    TH2F			*fDphiULSMassLow3;	//Dphi - ULS, mass< mass cut
    TH2F			*fDphiULSMassLow4;	//Dphi - ULS, mass< mass cut
    TH2F        *fDphiLSMassLow;  //Dphi - LS, mass< mass cut
    TH2F        *fDphiLSMassLow1;  //Dphi - LS, mass< mass cut
    TH2F        *fDphiLSMassLow2;  //Dphi - LS, mass< mass cut
    TH2F        *fDphiLSMassLow3;  //Dphi - LS, mass< mass cut
    TH2F        *fDphiLSMassLow4;  //Dphi - LS, mass< mass cut
    TH2F        *fDphiULSMassLowNoPartner; //Dphi - ULS, mass< mass cut no partner
    TH2F        *fDphiULSMassLowNoPartner1; //Dphi - ULS, mass< mass cut no partner
    TH2F        *fDphiULSMassLowNoPartner2; //Dphi - ULS, mass< mass cut no partner
    TH2F        *fDphiULSMassLowNoPartner3; //Dphi - ULS, mass< mass cut no partner
    TH2F        *fDphiULSMassLowNoPartner4; //Dphi - ULS, mass< mass cut no partner
    TH2F			*fDphiLSMassLowNoPartner;	//Dphi - LS, mass< mass cut
    TH2F			*fDphiLSMassLowNoPartner1;	//Dphi - LS, mass< mass cut
    TH2F			*fDphiLSMassLowNoPartner2;	//Dphi - LS, mass< mass cut
    TH2F			*fDphiLSMassLowNoPartner3;	//Dphi - LS, mass< mass cut
    TH2F			*fDphiLSMassLowNoPartner4;	//Dphi - LS, mass< mass cut
    TH1F			*fPhotoElecPt;		//photonic elec pt 
    TH1F			*fSemiInclElecPt;	//Semi inclusive ele pt
    TH1F        *fInclusiveElecPt; // Inclusive elec pt
    TH1F        *fULSElecPt; //ULS elec Pt
    TH1F        *fLSElecPt;// LS elec pt 
    //Eta bins (Deta < 1)
    TH2F       *fSemiIncElecDphiEta1;   //Semi Inclusive elec - had DPhi
    TH2F       *fSemiIncElecDphiEta11;     //Semi Inclusive elec - had DPhi
    TH2F       *fSemiIncElecDphiEta12;     //Semi Inclusive elec - had DPhi
    TH2F       *fSemiIncElecDphiEta13;     //Semi Inclusive elec - had DPhi
    TH2F       *fSemiIncElecDphiEta14;     //Semi Inclusive elec - had DPhi
    TH2F       *fPhotElecDphiEta1;   //Photon elec - had DPhi
    TH2F       *fPhotElecDphiEta11;     //Photon elec - had DPhi
    TH2F       *fPhotElecDphiEta12;     //Photon elec - had DPhi
    TH2F       *fPhotElecDphiEta13;     //Photon elec - had DPhi
    TH2F       *fPhotElecDphiEta14;     //Photon elec - had DPhi
    TH2F       *fInclusiveElecDphiEta1;    //Inclusive elec - had DPhi
    TH2F       *fInclusiveElecDphiEta11;   //Inclusive elec - had DPhi
    TH2F       *fInclusiveElecDphiEta12;   //Inclusive elec - had DPhi
    TH2F       *fInclusiveElecDphiEta13;   //Inclusive elec - had DPhi
    TH2F       *fInclusiveElecDphiEta14;   //Inclusive elec - had DPhi
    TH2F       *fDphiULSMassLowEta1; //Dphi - ULS, mass< mass cut
    TH2F       *fDphiULSMassLowEta11;   //Dphi - ULS, mass< mass cut
    TH2F       *fDphiULSMassLowEta12;   //Dphi - ULS, mass< mass cut
    TH2F       *fDphiULSMassLowEta13;   //Dphi - ULS, mass< mass cut
    TH2F       *fDphiULSMassLowEta14;   //Dphi - ULS, mass< mass cut
    TH2F        *fDphiLSMassLowEta1;  //Dphi - LS, mass< mass cut
    TH2F        *fDphiLSMassLowEta11;  //Dphi - LS, mass< mass cut
    TH2F        *fDphiLSMassLowEta12;  //Dphi - LS, mass< mass cut
    TH2F        *fDphiLSMassLowEta13;  //Dphi - LS, mass< mass cut
    TH2F        *fDphiLSMassLowEta14;  //Dphi - LS, mass< mass cut
    TH2F        *fDphiULSMassLowNoPartnerEta1; //Dphi - ULS, mass< mass cut no partner
    TH2F        *fDphiULSMassLowNoPartnerEta11; //Dphi - ULS, mass< mass cut no partner
    TH2F        *fDphiULSMassLowNoPartnerEta12; //Dphi - ULS, mass< mass cut no partner
    TH2F        *fDphiULSMassLowNoPartnerEta13; //Dphi - ULS, mass< mass cut no partner
    TH2F        *fDphiULSMassLowNoPartnerEta14; //Dphi - ULS, mass< mass cut no partner
    TH2F       *fDphiLSMassLowNoPartnerEta1;  //Dphi - LS, mass< mass cut
    TH2F       *fDphiLSMassLowNoPartnerEta11; //Dphi - LS, mass< mass cut
    TH2F       *fDphiLSMassLowNoPartnerEta12; //Dphi - LS, mass< mass cut
    TH2F       *fDphiLSMassLowNoPartnerEta13; //Dphi - LS, mass< mass cut
    TH2F       *fDphiLSMassLowNoPartnerEta14; //Dphi - LS, mass< mass cut

    //Eta bins (Deta > 1)
    TH2F       *fSemiIncElecDphiEta2;   //Semi Inclusive elec - had DPhi
    TH2F       *fSemiIncElecDphiEta21;     //Semi Inclusive elec - had DPhi
    TH2F       *fSemiIncElecDphiEta22;     //Semi Inclusive elec - had DPhi
    TH2F       *fSemiIncElecDphiEta23;     //Semi Inclusive elec - had DPhi
    TH2F       *fSemiIncElecDphiEta24;     //Semi Inclusive elec - had DPhi
    TH2F       *fPhotElecDphiEta2;   //Photon elec - had DPhi
    TH2F       *fPhotElecDphiEta21;     //Photon elec - had DPhi
    TH2F       *fPhotElecDphiEta22;     //Photon elec - had DPhi
    TH2F       *fPhotElecDphiEta23;     //Photon elec - had DPhi
    TH2F       *fPhotElecDphiEta24;     //Photon elec - had DPhi
    TH2F       *fInclusiveElecDphiEta2;    //Inclusive elec - had DPhi
    TH2F       *fInclusiveElecDphiEta21;   //Inclusive elec - had DPhi
    TH2F       *fInclusiveElecDphiEta22;   //Inclusive elec - had DPhi
    TH2F       *fInclusiveElecDphiEta23;   //Inclusive elec - had DPhi
    TH2F       *fInclusiveElecDphiEta24;   //Inclusive elec - had DPhi
    TH2F       *fDphiULSMassLowEta2; //Dphi - ULS, mass< mass cut                       
    TH2F       *fDphiULSMassLowEta21;   //Dphi - ULS, mass< mass cut                        
    TH2F       *fDphiULSMassLowEta22;   //Dphi - ULS, mass< mass cut                            
    TH2F       *fDphiULSMassLowEta23;   //Dphi - ULS, mass< mass cut
    TH2F       *fDphiULSMassLowEta24;   //Dphi - ULS, mass< mass cut
    TH2F        *fDphiLSMassLowEta2;  //Dphi - LS, mass< mass cut
    TH2F        *fDphiLSMassLowEta21;  //Dphi - LS, mass< mass cut
    TH2F        *fDphiLSMassLowEta22;  //Dphi - LS, mass< mass cut
    TH2F        *fDphiLSMassLowEta23;  //Dphi - LS, mass< mass cut
    TH2F        *fDphiLSMassLowEta24;  //Dphi - LS, mass< mass cut
    TH2F        *fDphiULSMassLowNoPartnerEta2; //Dphi - ULS, mass< mass cut no partner
    TH2F        *fDphiULSMassLowNoPartnerEta21; //Dphi - ULS, mass< mass cut no partner
    TH2F        *fDphiULSMassLowNoPartnerEta22; //Dphi - ULS, mass< mass cut no partner
    TH2F        *fDphiULSMassLowNoPartnerEta23; //Dphi - ULS, mass< mass cut no partner
    TH2F        *fDphiULSMassLowNoPartnerEta24; //Dphi - ULS, mass< mass cut no partner
    TH2F       *fDphiLSMassLowNoPartnerEta2;  //Dphi - LS, mass< mass cut
    TH2F       *fDphiLSMassLowNoPartnerEta21; //Dphi - LS, mass< mass cut
    TH2F       *fDphiLSMassLowNoPartnerEta22; //Dphi - LS, mass< mass cut
    TH2F       *fDphiLSMassLowNoPartnerEta23; //Dphi - LS, mass< mass cut
    TH2F       *fDphiLSMassLowNoPartnerEta24; //Dphi - LS, mass< mass cut

    //   TH1F			*fTrackPtBefTrkCuts;	//Track pt before track cuts	
    //   TH1F			*fTrackPtAftTrkCuts;	//Track pt after track cuts
    //   TH2F			*fTPCnsigma;		//TPC n sigma vs p	
    //   TH1F			*fNCellv1;		//No of cells in cluster, all EMCAL cluster
    //   TH1F			*fClsEv1;		//Cluster energy, all EMCAL cluster
    //   TH1F			*fNClusv1;		//No of clusters in event, all EMCAL cluster
    TH1F        *fInvmassLS1; //LS Invmass for all rec par
    //   TH1F        *fInvmassLS2; //LS Invmass for all rec par
    //   TH1F        *fInvmassLS3; //LS Invmass for all rec par
    //   TH1F        *fInvmassLS4; //LS Invmass for all rec par
    //   TH1F        *fInvmassLS5; //LS Invmass for all rec par
    TH1F        *fInvmassULS1;//ULS Invmass for all rec par
    //   TH1F        *fInvmassULS2;//ULS Invmass for all rec par
    //   TH1F        *fInvmassULS3;//ULS Invmass for all rec par
    //   TH1F        *fInvmassULS4;//ULS Invmass for all rec par
    //   TH1F        *fInvmassULS5;//ULS Invmass for all rec par
    TH1F        *fcentrality;//
    TH1F        *fElecPhi;//
    TH1F        *fElecPhiTPChalf;//
    TH2F        *fElecPhiPt;//
    //    TH1F        *fElecPhiTPC;//
    //    TH1F        *fElecPhiTPCEovP;//
    TH1F        *fHadronPhi;//
    TH1F        *fHadronPhiTPChalf;//
    TH2F        *fHadronPhiPt;//
    /*    TH1F        *fTrackHFEcuts;//
          TH1F        *fTrakPhiSPD1;//
          TH1F        *fTrakPhiSPD2;//
          TH1F        *fTrakPhiSPDOr;//
          TH1F        *fTrakPhiSPDAnd;//
          TH1F        *fTrackHFEcutsITS;//
     */
    /*    TH1F        *fNoMixedEvents;//
          TH2F        *fMixStat; //no of events in pool vs multplicity
          TH2F        *fMixStat1; //no of events in pool vs zvtx 
          TH2F        *fMixedIncElecDphi; //Mixed event - inclusive elec DPhi
          TH2F        *fMixedIncElecDphi1; //Mixed event - inclusive elec DPhi
          TH2F        *fMixedIncElecDphi2; //Mixed event - inclusive elec DPhi
          TH2F        *fMixedPhotElecDphi; //
          TH2F        *fMixedPhotElecDphi1; //
          TH2F        *fMixedPhotElecDphi2; //
          TH2F        *fMixedSemiIncElecDphi; //
          TH2F        *fMixedSemiIncElecDphi1; //
          TH2F        *fMixedSemiIncElecDphi2; //
          TH2F        *fMixedDphiULSMassLow;//
          TH2F        *fMixedDphiULSMassLow1;//
          TH2F        *fMixedDphiULSMassLow2;//
          TH2F        *fMixedDphiLSMassLow;//
          TH2F        *fMixedDphiLSMassLow1;//
          TH2F        *fMixedDphiLSMassLow2;//
     */
    TH1F        *fHadronPt;//
    TH1F       *fCentralityPass; // ! QA histogram of events that pass centrality cut
    TH1F       *fCentralityNoPass; //! QA histogram of events that do not pass centrality cut

    TH2F       *fHadronDphi;    //Hadron - had DPhi
    TH2F       *fHadronDphi1;   //Hadron - had DPhi
    TH2F       *fHadronDphi2;   //Hadron - had DPhi
    TH2F       *fHadronDphi3;   //Hadron - had DPhi
    TH2F       *fHadronDphi4;   //Hadron - had DPhi
    TH1F       *fPiPt; //TPC nsig < 3.5 pt

    TH2F       *fHadronDphiNoSS;    //Hadron - had DPhi
    TH2F       *fHadronDphiNoSS1;   //Hadron - had DPhi
    TH2F       *fHadronDphiNoSS2;   //Hadron - had DPhi
    TH2F       *fHadronDphiNoSS3;   //Hadron - had DPhi
    TH2F       *fHadronDphiNoSS4;   //Hadron - had DPhi
    TH1F       *fPiPtNoSS; //TPC nsig < 3.5 pt
    TH2F       *fEovPWoSS;//
    TH2F       *fEovPWSS;//
    TH2F       *fEovPHadWoSS;//
    TH2F       *fEovPHadWSS;//

    //Deta < 0.8
    TH2F       *fHadronDphiEta1;   //Hadron - had DPhi
    TH2F       *fHadronDphiEta11;     //Hadron - had DPhi
    TH2F       *fHadronDphiEta12;     //Hadron - had DPhi
    TH2F       *fHadronDphiEta13;     //Hadron - had DPhi
    TH2F       *fHadronDphiEta14;     //Hadron - had DPhi
    TH2F       *fHadronDphiNoSSEta1;   //Hadron - had DPhi NoSS
    TH2F       *fHadronDphiNoSSEta11;     //Hadron - had DPhi NoSS
    TH2F       *fHadronDphiNoSSEta12;     //Hadron - had DPhi NoSS
    TH2F       *fHadronDphiNoSSEta13;     //Hadron - had DPhi NoSS
    TH2F       *fHadronDphiNoSSEta14;     //Hadron - had DPhi NoSS

    //Deta > 0.8
    TH2F       *fHadronDphiEta2;   //Hadron - had DPhi
    TH2F       *fHadronDphiEta21;     //Hadron - had DPhi
    TH2F       *fHadronDphiEta22;     //Hadron - had DPhi
    TH2F       *fHadronDphiEta23;     //Hadron - had DPhi
    TH2F       *fHadronDphiEta24;     //Hadron - had DPhi
    TH2F       *fHadronDphiNoSSEta2;   //Hadron - had DPhi NoSS
    TH2F       *fHadronDphiNoSSEta21;     //Hadron - had DPhi NoSS
    TH2F       *fHadronDphiNoSSEta22;     //Hadron - had DPhi NoSS
    TH2F       *fHadronDphiNoSSEta23;     //Hadron - had DPhi NoSS
    TH2F       *fHadronDphiNoSSEta24;     //Hadron - had DPhi NoSS


    //THnSparse  *fSparseElectron;//!Electron info 
    //Double_t *fvalueElectron;//!Electron info 

    AliAnalysisTaskElecHadronCorrel(const AliAnalysisTaskElecHadronCorrel&); // not implemented
    AliAnalysisTaskElecHadronCorrel& operator=(const AliAnalysisTaskElecHadronCorrel&); // not implemented

    ClassDef(AliAnalysisTaskElecHadronCorrel, 2); //!example of analysis
};
/*
   class AliehDPhiBasicParticle : public AliVParticle
   {
   public:
   AliehDPhiBasicParticle(Float_t eta, Float_t phi, Float_t pt, Short_t charge)
   : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge)
   {
   }
   ~AliehDPhiBasicParticle() {}

// kinematics
virtual Double_t Px() const { AliFatal("Not implemented"); return 0; }
virtual Double_t Py() const { AliFatal("Not implemented"); return 0; }
virtual Double_t Pz() const { AliFatal("Not implemented"); return 0; }
virtual Double_t Pt() const { return fpT; }
virtual Double_t P() const { AliFatal("Not implemented"); return 0; }
virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

virtual Double_t Xv() const { AliFatal("Not implemented"); return 0; }
virtual Double_t Yv() const { AliFatal("Not implemented"); return 0; }
virtual Double_t Zv() const { AliFatal("Not implemented"); return 0; }
virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

virtual Double_t OneOverPt()  const { AliFatal("Not implemented"); return 0; }
virtual Double_t Phi()        const { return fPhi; }
virtual Double_t Theta()      const { AliFatal("Not implemented"); return 0; }


virtual Double_t E()          const { AliFatal("Not implemented"); return 0; }
virtual Double_t M()          const { AliFatal("Not implemented"); return 0; }

virtual Double_t Eta()        const { return fEta; }
virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }

virtual Short_t Charge()      const { return fCharge; }
virtual Int_t   GetLabel()    const { AliFatal("Not implemented"); return 0; }
// PID
virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }
virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }

private:
Float_t fEta;      // eta
Float_t fPhi;      // phi
Float_t fpT;       // pT
Short_t fCharge;   // charge

ClassDef( AliehDPhiBasicParticle, 1); // class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};
 */
#endif


