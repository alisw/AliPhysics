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

    void SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };
    void SetOpeningAngleCut (Double_t openingAngle) {fOpeningAngleCut = openingAngle;};
    void SetInvariantMassCut (Double_t invmass) {fInvmassCut = invmass;};
    AliHFEpid *GetPID() const { return fPID; }
    void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) { fRejectKinkMother = rejectKinkMother; };
    void SelectPhotonicElectron(Int_t itrack, AliESDtrack *track, Bool_t &fFlagPhotonicElec);
    void ElectronHadCorrel(Int_t itrack, AliESDtrack *track, TH2F *DphiPt, TH2F *DphiPt1,TH2F *DphiPt2,TH2F *DphiPt3,TH2F *DphiPt4);	
    void ElectronHadCorrelNoPartner(Int_t itrack,Int_t jtrack, AliESDtrack *track, TH2F *DphiPtNew,TH2F *DphiPtNew1,TH2F *DphiPtNew2,TH2F *DphiPtNew3,TH2F *DphiPtNew4);	
//    void MixedEvent(AliESDtrack *track, TH2F *DphiPt, TH2F *DphiPt1, TH2F *DphiPt2);
    void HadronInfo(Int_t itrack);
    void    SetCentralityParameters(Double_t CentralityMin, Double_t CentralityMax, const char* CentralityMethod); //select centrality
    void    CheckCentrality(AliESDEvent *event,Bool_t &centralitypass); //to use only events with the correct centrality....
//    TObjArray* CloneAndReduceTrackList();
  private:

    Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);

    AliESDEvent 		   *fESD;			//ESD object
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
    Double_t 		      fOpeningAngleCut;	//openingAngle cut value
    Double_t	         	fInvmassCut;		//invariant mass cut value
    Double_t             fCentrality; // event centrality for QA
    Double_t             fCentralityMin; // lower bound of cenrality bin
    Double_t             fCentralityMax; // upper bound of centrality bin
    const char           *fkCentralityMethod; // method used to determine centrality (V0 by default)   
//    AliEventPoolManager*     fPoolMgr;         //! event pool manager

    TH1F			*fNoEvents;		//no of events
//    TH1F			*fTrkpt;		//track pt
    TH2F			*fTrkEovPBef;		//track E/p before HFE pid
    TH2F			*fTrkEovPBefHad;		//track E/p before HFE pid
/*    TH2F			*fTrkEovPAft;		//track E/p after HFE pid
    TH2F			*fTrkEovPAftOwn;		//track E/p after HFE pid
    TH2F			*fdEdxBef;		//track dEdx vs p before HFE pid
    TH2F			*fdEdxAft;		//track dEdx vs p before HFE pid
    TH2F			*fdEdxAftOwn;		//track dEdx vs p before HFE pid
    TH1F			*fOpeningAngleLS;	//opening angle for LS pairs
    TH1F			*fOpeningAngleULS;	//opening angle for ULS pairs
*/    TH2F			*fSemiIncElecDphi;  	//Semi Inclusive elec - had DPhi
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

//    TH1F			*fTrackPtBefTrkCuts;	//Track pt before track cuts	
//    TH1F			*fTrackPtAftTrkCuts;	//Track pt after track cuts
    TH2F			*fTPCnsigma;		//TPC n sigma vs p	
//    TH2F			*fTPCnsigmaAft;		//TPC n sigma vs p	
//    TH2F			*fTPCnsigmaAftOwn;		//TPC n sigma vs p	
//    TH1F			*fNCellv1;		//No of cells in cluster, all EMCAL cluster
//    TH1F			*fClsEv1;		//Cluster energy, all EMCAL cluster
//    TH1F			*fNClusv1;		//No of clusters in event, all EMCAL cluster

//    TH1F        *fKFParticleP; //KFparticle rec P distr
//    TH1F        *fKFParticleE; //KFparticle rec E distr
    TH1F        *fInvmassLS1; //LS Invmass for all rec par
    TH1F        *fInvmassLS2; //LS Invmass for all rec par
    TH1F        *fInvmassLS3; //LS Invmass for all rec par
    TH1F        *fInvmassLS4; //LS Invmass for all rec par
    TH1F        *fInvmassLS5; //LS Invmass for all rec par
    TH1F        *fInvmassULS1;//ULS Invmass for all rec par
    TH1F        *fInvmassULS2;//ULS Invmass for all rec par
    TH1F        *fInvmassULS3;//ULS Invmass for all rec par
    TH1F        *fInvmassULS4;//ULS Invmass for all rec par
    TH1F        *fInvmassULS5;//ULS Invmass for all rec par
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
/*    TH1F        *fNLSminus;//
    TH1F        *fNLSplus;//
    TH1F        *fNULS;//  
*/    TH1F        *fHadronIPxy;//  
    TH1F        *fHadronIPz;//  
    TH1F        *fHadronPt;//
    TH1F       *fCentralityPass; // ! QA histogram of events that pass centrality cut
    TH1F       *fCentralityNoPass; //! QA histogram of events that do not pass centrality cut
    //  THnSparse  *fSparseElectron;//!Electron info 
    //  Double_t *fvalueElectron;//!Electron info 

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


