/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIANALYSISTASKHAHFECORREL_H
#define ALIANALYSISTASKHAHFECORREL_H 

class THnSparse;
class TH2F;
class TLorentzVector;

class AliEMCALTrack;
class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliAODMCParticle;
class AliAODMCHeader;
class AliGenEventHeader;
class AliVEvent;
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
class AliPIDResponse;
class AliMultSelection;
class AliEventPoolManager;


#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskHaHFECorrel : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskHaHFECorrel();
    AliAnalysisTaskHaHFECorrel(const char *name);
    ~AliAnalysisTaskHaHFECorrel();
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    
    //******************** ANALYSIS
    AliAODTrack* FindLPAndHFE(TObjArray* RedTracksHFE, const AliVVertex *pVtx, Int_t nMother, Double_t listMother[]);
    void FindPhotonicPartner(Int_t iTracks, AliAODTrack* track,  const AliVVertex *pVtx, Int_t nMother, Double_t listMother[], Int_t &LSPartner, Int_t&ULSPartner);
    void CorrelateElectron(TObjArray* RedTracksHFE);

    void CorrelateLP(AliAODTrack* LPtrack, TObjArray* RedTracksHFE);
    void CorrelateLPMixedEvent(AliAODTrack* LPtrack, Float_t mult, Float_t zVtx);
    
    void CorrelateHadron(TObjArray* RedTracksHFE,  const AliVVertex* pVtx, Int_t nMother, Double_t listMother[], Float_t mult);
    void CorrelateHadronMixedEvent(AliAODTrack* Htrack, Float_t mult, Float_t zVtx);

    //*********************ANALYSIS Helper
    Bool_t ChargedHadronTrackCuts(const AliVVertex *pVtx,AliAODTrack *Htrack, Int_t nMother, Double_t listMother[]);
    Bool_t ChargedHadronPIDCuts(AliAODTrack *Htrack);

    Bool_t InclElecTrackCuts(const AliVVertex *pVtx,AliAODTrack *ietrack, Int_t nMother, Double_t listMother[]);
    Bool_t InclElecPIDCuts(AliAODTrack *track, Bool_t IsPrimary);

    Bool_t PhotElecPIDCuts(AliAODTrack *track);
    Bool_t PhotElecTrackCuts(const AliVVertex *pVtx,AliAODTrack *aetrack, Int_t nMother, Double_t listMother[]);
    
    Bool_t CloneAndReduceTrackList(TObjArray* RedTracks, AliAODTrack* track, Int_t LSPartner, Int_t ULSPartner);

    void BinLogX(TAxis *axis);



    //**************  SETTINGS
    void SetTender (Bool_t UseTender) {fUseTender = UseTender;};
    void SetPeriod (Double_t period) {fWhichPeriod = period;};

    void SetAssPtCut (Double_t AssPtCut) {fAssPtCut = AssPtCut;};
    void SetITSnCut (Int_t ITSnCut) {fITSnCut = ITSnCut;};
    void SetAssTPCnCut (Int_t AssTPCnCut) {fAssTPCnCut = AssTPCnCut;};
    void SetTPCnCut(Int_t TPCnCut) {fTPCnCut = TPCnCut;};
    void SetHTPCnCut(Int_t HTPCnCut) {fHTPCnCut = HTPCnCut;}
    void SetAssITSrefitCut(Bool_t AssITSrefitCut) {fAssITSrefitCut = AssITSrefitCut;};
    void SetHITSrefitCut(Bool_t HITSrefitCut) {fHITSrefitCut = HITSrefitCut;};
    
    void SetHTPCrefitCut(Bool_t HTPCrefitCut) {fHTPCrefitCut = HTPCrefitCut;};

    void SetSigmaITScut(Double_t SigmaITScut) {fSigmaITScut = SigmaITScut;};
    void SetSigmaTOFcut(Double_t SigmaTOFcut) {fSigmaTOFcut = SigmaTOFcut;};
    void SetSigmaTPCcut(Double_t SigmaTPCcut) {fSigmaTPCcut = SigmaTPCcut;};

    // void SetWeightSyst(Bool_t WeightSyst) {fWeightSyst = WeightSyst;};
    // void SetSystTOFcut(Bool_t SystTOFcut) {fSystTOFcut = SystTOFcut;};
  
    //void SetPileUpCut(Bool_t EnablePileupRejVZEROTPCout){fEnablePileupRejVZEROTPCout = EnablePileupRejVZEROTPCout;};
    void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) { fRejectKinkMother = rejectKinkMother; };

    void SetHadronCorrelation(Bool_t CorrHadron) {fCorrHadron = CorrHadron;};
    void SetLPCorrelation(Bool_t CorrLP) {fCorrLParticle = CorrLP;};
    
    void SetOpeningAngleCut(Bool_t OpeningAngleCut) {fOpeningAngleCut=OpeningAngleCut;};
    void SetInvmassCut(Double_t InvmassCut) {fInvmassCut=InvmassCut;};
 

 private:
    
    Bool_t                ProcessCutStep(Int_t cutStep, AliVParticle *track);
    Double_t              GetDeltaPhi(Double_t phiA,Double_t phiB) const;
    Double_t              GetDeltaEta(Double_t etaA,Double_t etaB) const;
    
    Bool_t                fUseTender;               // Use tender
    Int_t                 fWhichPeriod;             // period
    Double_t              fAssPtCut;                // pt cut for associated electron
    Int_t                 fITSnCut;                 // ITs number of clusters for tagged electrons
    Int_t                 fAssTPCnCut;              // TPC number of clusters for associated electron
    Int_t                 fTPCnCut;                 // TPC number of clusters for tagged electron
    Int_t                 fHTPCnCut;                // TPC number of clusters for trigger hadron
    Bool_t                fAssITSrefitCut;          // ITS refit for associated electron
    Bool_t                fHITSrefitCut;            // ITS refit for trigger hadron
    Bool_t                fHTPCrefitCut;            // TPC refit for trigger hadron
    Double_t              fSigmaITScut;             // ITS nSigma cut
    Double_t              fSigmaTOFcut;             // TOF nSigma cut
    Double_t              fSigmaTPCcut;             // TPC nSigma cut 

    Bool_t                fWeightSyst;              // used?
    Bool_t                fSystTOFcut;              // used?

    Bool_t                fEnablePileupRejVZEROTPCout;   // used?
    Bool_t                fRejectKinkMother;        // Reject Kink Mother
    
    Double_t              fOpeningAngleCut;         // openingAngle cut for non-HFE selection
    Double_t              fInvmassCut;              // invariant mass cut  for non-HFE selection
    Double_t              fChi2Cut;                 //! used?? Chi2 cut  for non-HFE selection
    Double_t              fDCAcut;                  //! used?? DCA cut  for non-HFE selection
  
    // ******* Switch for analysis modes
    Bool_t                fCorrHadron;              // Choste Hadron-HFE Correl
    Bool_t                fCorrLParticle;           // Chose LP-HFE Correl
    Bool_t                fMixedEvent;              // Fill Mixed Event for the cases chosen above
    Bool_t                fLParticle;               //! Is LP found?
    

    AliAODEvent           *fAOD;                    //! AOD object
    AliVEvent             *fVevent;                 //! VEvent
    AliPIDResponse        *fpidResponse;            //! PID response
    AliMultSelection      *fMultSelection;          //! MulSelection
    AliCentrality         *fCentrality;             //! Centrality

    AliEventPoolManager   *fPoolMgr;                //! event pool manager
    
    AliMCEvent            *fMC;                     //! MC object
    AliStack              *fStack;                  //! stack
    AliAODMCParticle      *fMCparticle;             //! MC particle
    TClonesArray          *fMCarray;                //! MC array
    AliAODMCHeader        *fMCheader;               //! MC header

    
    TClonesArray          *fTracks_tender;          //Tender tracks
    TClonesArray          *fCaloClusters_tender;    //Tender clusters
    
    
  
    
    AliHFEcuts            *fCuts;                   //! Cut Collection
    
    Bool_t                fIdentifiedAsOutInz;      //! used?? Out Of Range in z
    Bool_t                fPassTheEventCut;         //! used?? Pass The Event Cut
   
    Bool_t                fIsMC;                    //! flag for MC analysis
    Bool_t                fIsAOD;                   //! flag for AOD analysis
    Bool_t                fSetMassConstraint;       //! used?? set mass constraint
    Double_t              fVz;                      //! used?? z position of the primary vertex
    AliCFManager          *fCFM;                    //! Correction Framework Manager
    AliHFEpid             *fPID;                    //! PID
    AliHFEpidQAmanager    *fPIDqa;                  //! PID QA manager
    
    Int_t                 fWhichDecay;              //! used?? which decay
    Double_t              fPi0EtaWeight;            //! used?? weight for pi0/eta in MC with enhanced signal
    
    TList                 *fOutputList;             //! output list
    TH1F                  *fNoEvents;               //! no of events for different cuts
    TH2F                  *fTrkpt;                  //! track pt for different cuts
    TH1F                  *fMultiplicity;	    //! multiplicity distribution
    
    TH2F                  *fHistITSnSig;            //! ITS sigma vs p
    TH2F                  *fHistTOFnSig;            //! TOF sigma vs p
    TH2F                  *fHistTPCnSig;            //! TPC sigma vs p
    TH2F                  *fHistTPCnSigITScut;      //! TPC sigma vs p (ITS cut)
    TH2F                  *fHistTPCnSigTOFcut;      //! TPC sigma vs p (TOF cut)
    TH2F                  *fHistTPCnSigITSTOFcut;   //! TPC sigma vs p (ITS+TOF cuts)
    

    TH1F                  *fInclElecPt;              //! inclusive electron pt
    TH1F                  *fULSElecPt;              //! ULS electron pt (after IM cut)
    TH1F                  *fLSElecPt;               //! LS electron pt (after IM cut)
    TH1F                  *fElecTrigger;            //! trigger electron vs pt
    TH2F                  *fInclElecPhi;            //! electron (trigger): phi vs pt
    TH2F                  *fULSElecPhi;             //! phi vs pt for electrons from ULS pairs
    TH2F                  *fLSElecPhi;              //! phi vs pt for electrons from LS pairs
    TH2F                  *fElecDphi;               //! inlcusive electron: dPhi vs pt of triggered electron
    TH2F                  *fULSElecDphi;            //! electron from ULS pairs: dPhi vs pt of triggered electron
    TH2F                  *fLSElecDphi;             //! electron from LS pairs: dPhi vs pt of triggered electron
    TH2F                  *fULSElecDphiDiffMethod;   //! electron from ULS pairs: dPhi vs pt of triggered electron
    TH2F                  *fLSElecDphiDiffMethod;   //! electron from LS pairs: dPhi vs pt of triggered electron
    
    TH1F                  *fElecHaTrigger;          //!
    THnSparse             *fElecHa;                 //!
    THnSparse             *fLSElecHa;               //!
    THnSparse             *fULSElecHa;              //!

    THnSparse             *fElecHaLSNoPartner;      //!
    THnSparse             *fElecHaULSNoPartner;     //!

    THnSparse             *fElecHaMixedEvent;       //!
    THnSparse             *fLSElecHaMixedEvent;     //!
    THnSparse             *fULSElecHaMixedEvent;    //!

    TH1F                  *fElecLPTrigger;          //!
   
    THnSparse             *fElecLP;                 //!
    THnSparse             *fLSElecLP;               //!
    THnSparse             *fULSElecLP;              //!

    THnSparse             *fElecLPMixedEvent;       //!
    THnSparse             *fLSElecLPMixedEvent;     //!
    THnSparse             *fULSElecLPMixedEvent;    //!


    TH2F                  *fInvmassLS;              //! Inv mass of LS (e,e)
    TH2F                  *fInvmassULS;             //! Inv mass of ULS (e,e)
    TH2F                  *fOpeningAngleLS;         //! opening angle for LS pairs
    TH2F                  *fOpeningAngleULS;        //! opening angle for ULS pairs
    TH2F                  *fCheckLSULS;              //!
    
    TH2F                  *fPi0Pt;                  //! primary pi0 pt to compute the weight
    TH2F                  *fEtaPt;                  //! primary eta pt to compute the weight
    
   
    
    AliAnalysisTaskHaHFECorrel(const AliAnalysisTaskHaHFECorrel&);
    AliAnalysisTaskHaHFECorrel& operator=(const AliAnalysisTaskHaHFECorrel&);
    
    ClassDef(AliAnalysisTaskHaHFECorrel, 1);
};


// class storing reduced track information for mixed event pool

class AliBasicParticleHaHFE : public AliVParticle
{
 public:
 AliBasicParticleHaHFE() 
   : fID(0), fEta(0), fPhi(0), fpT(0), fCharge(0), fLSpartner(0), fULSpartner(0) 
    {}
 AliBasicParticleHaHFE(Int_t id, Float_t eta, Float_t phi, Float_t pt, Short_t charge, Short_t LS, Short_t ULS)
   : fID(id), fEta(eta), fPhi(phi), fpT(pt), fCharge(charge), fLSpartner(LS), fULSpartner(ULS)
  {
  }
  virtual ~AliBasicParticleHaHFE() {}

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
  
  virtual Short_t LS()          const {return fLSpartner; }
  virtual Short_t ULS()         const {return fULSpartner; }
  virtual Short_t ID()          const {return fID;}
  
 private:
  Float_t fEta;        // eta
  Float_t fPhi;        // phi
  Float_t fpT;         // pT
  Short_t fCharge;     // charge
  Short_t fULSpartner; // no of ULS partner
  Short_t fLSpartner;  // no of LS partner
  Int_t   fID;         // particle id

  ClassDef(AliBasicParticleHaHFE, 1); // class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};



#endif
