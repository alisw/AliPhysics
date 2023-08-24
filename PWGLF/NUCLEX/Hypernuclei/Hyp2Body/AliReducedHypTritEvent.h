/// \class AliReducedHypTritEvent
/// \brief Reduced Event class
///
/// Reduced Event class for saving hypertriton information in a TTree for later analysis.
/// Each Event contains an array of v0, which all have a positive and a negative track.
///
/// \author Lukas Kreis <l.kreis@gsi.de>, GSI
/// \date Feb 17, 2016


#ifndef AliReducedHypTritEVENT_H
#define AliReducedHypTritEVENT_H

#ifndef ROOT_Object
#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <TClonesArray.h>
#endif

class AliReducedHypTritTrack : public TObject {

  friend class AliAnalysisTaskHypTritEventTree;

 public:
  AliReducedHypTritTrack();
  ~AliReducedHypTritTrack();

  // Getters
  TLorentzVector P() const {return fP;}
  Double_t Ptrack() const {return fPtrack;}
  Double_t Dca() const {return fDca;}
  Double_t SignedDca() const {return fDcaSigned;}
  Double_t Phi() const {return fPhi;}
  Double_t Eta() const {return fEta;}
  Double_t Dedx() const {return fDedx;}
  Double_t DedxSigma() const {return fDedxSigma;}
  Double_t DedxSigmaTriton() const {return fDedxSigmaTriton;}
  
  Double_t TpcNcls() const {return fTpcNClusters;}
  Double_t ITSNcls() const {return fITSNClusters;}
  Double_t TpcChi2() const {return fTpcChi2;}
  Int_t    Kink() const {return fKink;}
  Int_t    TPCrefit() const {return fTPCrefit;}
  Int_t    ITSrefit() const {return fITSrefit;}
  Double_t GeoLength() const {return fGeoLength;}
  Int_t    TrkCutsPassed() const {return fTrkCutsPassed;}
  Int_t    TRDvalid() const {return fTRDvalid;}
  Int_t    TRDtrigHNU() const {return fTRDtrigHNU;}
  Int_t    TRDtrigHQU() const {return fTRDtrigHQU;}
  Int_t    TRDPid() const {return fTRDPid;}
  Int_t    TRDnTracklets() const {return fTRDnTracklets;}
  Int_t    TRDPt() const {return fTRDPt;}
  Int_t    TRDLayerMask() const {return fTRDLayerMask;}
  Double_t TRDSagitta() const {return fTRDSagitta;}
  Int_t    TRDStack() const {return fTRDStack;}
  Int_t    TRDSector() const {return fTRDSector;}
  UInt_t   TRDPID0() const {return fTRDPID0;} 
  UInt_t   TRDPID1() const {return fTRDPID1;} 
  UInt_t   TRDPID2() const {return fTRDPID2;} 
  UInt_t   TRDPID3() const {return fTRDPID3;} 
  UInt_t   TRDPID4() const {return fTRDPID4;} 
  UInt_t   TRDPID5() const {return fTRDPID5;}    
private:
  TLorentzVector  fP;               //< 4 momentum of track
  Double_t        fPtrack;          //< Total momentum of Track
  Double_t        fDca;             //< DCA to prim vertex
  Double_t        fDcaSigned;       //< DCA to prim vertex
  Double_t        fDedx;            //< specific energyloss in TPC of track
  Double_t        fDedxSigma;       //< dEdx sigma
  Double_t        fDedxSigmaTriton; //< dEdx sigma of triton hypothesis
  Double_t        fEta;             //< eta of track
  Double_t        fPhi;             //< phi of track
 
  Double_t        fTpcNClusters;    //< number of clusters
  Double_t        fITSNClusters;    //< number of clusters
  Double_t        fTpcChi2;         //< chi2 of TPC fit
  Int_t           fKink;            //< kink doughters
  Int_t           fTPCrefit;        //< TPC refit
  Int_t           fITSrefit;        //< ITS refit
  Double_t        fGeoLength;       //< geometric length cut
  Int_t           fTrkCutsPassed;   //< track passed track cuts
  Int_t           fTRDvalid;        //< has valid TRD track
  Int_t           fTRDtrigHNU;      //< HNU fired by track
  Int_t           fTRDtrigHQU;      //< HQU fired by track
  Int_t           fTRDPid;          //< PID value of TRD track
  Int_t           fTRDnTracklets;   //< number of TRD tracklets
  Int_t           fTRDPt;           //< Pt of TRD track
  Int_t           fTRDLayerMask;    //< TRD track layer mask
  Double_t        fTRDSagitta;      //< sagitta value of TRD track
  Int_t           fTRDStack;
  Int_t           fTRDSector;
  UInt_t          fTRDPID0;
  UInt_t          fTRDPID1;
  UInt_t          fTRDPID2;
  UInt_t          fTRDPID3;
  UInt_t          fTRDPID4;
  UInt_t          fTRDPID5;

AliReducedHypTritTrack(const AliReducedHypTritTrack&);
AliReducedHypTritTrack &operator = (const AliReducedHypTritTrack&);
ClassDef(AliReducedHypTritTrack, 8)
};

class AliReducedHypTritV0 : public TObject {

  friend class AliAnalysisTaskHypTritEventTree;

 public:
  AliReducedHypTritV0();
  ~AliReducedHypTritV0();

  // Getters
  TVector3 Position() const {return fPosition;}
  TVector3 Pvect() const {return fPvect;}
  AliReducedHypTritTrack* Pi() const {return fPiTrack;}
  AliReducedHypTritTrack* He() const {return fHeTrack;}
  Double_t M() const {return fM;}
  Double_t Pt() const {return fPt;}
  Double_t P() const {return fP;}
  Double_t Dca() const {return fDcaV0;}
  Double_t Ct() const {return fDecayLength;}
  Double_t Cos() const {return fCosPointingAngle;}
  Bool_t  Mc() const {return fMcTruth;}
  Double_t Y() const {return fRapidity;}
  Short_t Charge() const {return fCharge;}
  Short_t Species() const {return fParticleSpecies;}
  Bool_t  OnFlyStatus() const {return fOnFlyStatus;}


private:
  TVector3                fPosition;         //< Lorentzvector of v0
  TVector3                fPvect;            //< Momentum vector of mother
  AliReducedHypTritTrack* fPiTrack;          //< positive daughter of v0
  AliReducedHypTritTrack* fHeTrack;          //< negative daughter of v0
  Double_t                fP;                //< momentum of mother
  Double_t                fPt;               //<  transverse momentum of mother
  Double_t                fM;                //< reconstructed invariant mass of mother
  Double_t                fDcaV0;            //< DCA of daughters
  Double_t                fCosPointingAngle; //< cosine of pointing angle of vertex
  Double_t                fDecayLength;      //< decay radius of mother particle
  Bool_t                  fMcTruth;          //< Monte Carlo truth of mother type
  Double_t                fRapidity;         //< Rapidity of V0
  Short_t                 fCharge;           //< anti or particle
  Short_t                 fParticleSpecies;  //< particle species
  Bool_t                  fOnFlyStatus;      //< ontheflyStatus


  AliReducedHypTritV0(const AliReducedHypTritV0&);
  AliReducedHypTritV0 &operator = (const AliReducedHypTritV0&);
  ClassDef(AliReducedHypTritV0, 5);
};

class AliReducedHypTritEvent : public TObject {

  friend class AliAnalysisTaskHypTritEventTree;

 public:
  AliReducedHypTritEvent();
  ~AliReducedHypTritEvent();

  AliReducedHypTritV0* V0(Int_t i) const
      {return (i < fNumberV0s ? (AliReducedHypTritV0*) fV0s->At(i) : 0x0);}
  TVector3 VertexPosition() const {return fVertexPosition;}
  Double_t  Centrality() const {return fCentrality;}
  UShort_t NumberV0s() const {return fNumberV0s;}
 
  UShort_t Trigger() const {return fTrigger;}
  UShort_t IsMBtriggered() const {return fTrigMB;}
  UShort_t IsHNUtriggered() const {return fTrigHNU;}
  UShort_t IsHQUtriggered() const {return fTrigHQU;}
  UShort_t IsHJTtriggered() const {return fTrigHJT;}
  UShort_t IsHSEtriggered() const {return fTrigHSE;}
  UShort_t IsV0triggered() const {return fTrigV0;}
  UShort_t IsSPDtriggered() const {return fTrigSPD;}
  TString  TriggerClasses() const {return fTriggerClasses;}
 
  Int_t    RunNumber() const {return fRunNumber;}
  Double_t MagField()const {return fMagField;}

  Double_t SPDFiredChips0() const {return fSPDFiredChips0;}
  Double_t SPDFiredChips1() const {return fSPDFiredChips1;}
  Double_t SPDTracklets() const {return fSPDTracklets;}
  Double_t SPDCluster() const {return fSPDCluster;}
  Double_t V0Multiplicity() const {return fV0Multiplicity;}

  Double_t MultV0M()const {return fMultV0M;}
  Double_t MultOfV0M()const {return fMultOfV0M;}
  Double_t MultSPDTracklet() const {return fMultSPDTracklet;}
  Double_t MultSPDCluster() const {return fMultSPDCluster;}
  Double_t MultRef05()const {return fMultRef05;}
  Double_t MultRef08()const {return fMultRef08;}
  ULong64_t EventID() const {return fEventID;}
  Int_t EvCutsPassed() const {return fEvCutsPassed;} 
    
  void ClearEvent();

private:
  TVector3      fVertexPosition; //< position of primary vertex
  TClonesArray* fV0s;            //< array of v0s in event
  UShort_t      fNumberV0s;      //< number of v0s in event
  Double_t      fCentrality;     //< centrality of event

  UShort_t      fTrigger;        //< array of Triggers
  UShort_t      fTrigMB;         //< Flag for MB trigger
  UShort_t      fTrigHNU;        //< Flag for TRD-HNU trigger
  UShort_t      fTrigHQU;        //< Flag for TRD-HQU trigger
  UShort_t      fTrigHJT;        //< Flag for TRD-HJT trigger
  UShort_t      fTrigHSE;        //< Flag for TRD-HSE trigger
  UShort_t      fTrigV0;         //< Flag for HM-V0 trigger
  UShort_t      fTrigSPD;        //< Flag for HM-SPD trigger
  TString       fTriggerClasses; //< fired trigger classes

  Int_t         fEvCutsPassed;   //< event passed event cuts
  ULong64_t     fEventID;        //< global event ID
  Int_t         fRunNumber;      //< number of run
  Double_t      fMagField;

  Double_t      fSPDFiredChips0; // multiplicity triggers
  Double_t      fSPDFiredChips1;
  Double_t      fSPDTracklets;
  Double_t      fSPDCluster;
  Double_t      fV0Multiplicity;

  Double_t      fMultV0M;// multiplicity estimators
  Double_t      fMultOfV0M;
  Double_t      fMultSPDTracklet;
  Double_t      fMultSPDCluster;
  Double_t      fMultRef05;
  Double_t      fMultRef08;

  AliReducedHypTritEvent(const AliReducedHypTritEvent&);
  AliReducedHypTritEvent &operator = (const AliReducedHypTritEvent&);
  ClassDef(AliReducedHypTritEvent, 7);
};

#endif
