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
  Float_t Ptrack() const {return fPtrack;}
  Float_t Dca() const {return fDca;}
  Float_t Phi() const {return fPhi;}
  Float_t Eta() const {return fEta;}
  Float_t Dedx() const {return fDedx;}
  Float_t DedxSigma() const {return fDedxSigma;}
  Float_t DedxSigmaTriton() const {return fDedxSigmaTriton;}
  Float_t TpcNcls() const {return fTpcNClusters;}
  Float_t TpcChi2() const {return fTpcChi2;}
  Int_t 	Kink() const {return fKink;}
  Int_t 	TPCrefit() const {return fTPCrefit;}
private:
  TLorentzVector fP;               //< 4 momentum of track
  Float_t        fPtrack;          //< Total momentum of Track
  Float_t        fDca;             //< DCA to prim vertex
  Float_t        fDedx;            //< specific energyloss in TPC of track
  Float_t        fDedxSigma;       //< dEdx sigma
  Float_t        fDedxSigmaTriton; //< dEdx sigma of triton hypothesis
  Float_t        fEta;             //< eta of track
  Float_t        fPhi;             //< phi of track
  Float_t        fTpcNClusters;    //< number of clusters
  Float_t	 			 fTpcChi2;	   		 //< chi2 of TPC fit
  Int_t					 fKink;						 //< kink doughters
  Int_t					 fTPCrefit;				 //< TPC refit
  Float_t        fGeoLength;       //< geometric length cut
	Int_t					 fTRDvalid;	       //< has valid TRD track
	Int_t					 fTRDtrigHNU;	     //< HNU fired by track
	Int_t					 fTRDtrigHQU;			 //< HQU fired by track
	Int_t					 fTRDPid;					 //< PID value of TRD track
	Int_t	 				 fTRDnTracklets;	 //< number of TRD tracklets
	Int_t	         fTRDPt;           //< Pt of TRD track
	Int_t	         fTRDLayerMask;		 //< TRD track layer mask
	Float_t	       fTRDSagitta;			 //< sagitta value of TRD track
	
AliReducedHypTritTrack(const AliReducedHypTritTrack&);
AliReducedHypTritTrack &operator = (const AliReducedHypTritTrack&);
ClassDef(AliReducedHypTritTrack, 6)
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
  Float_t M() const {return fM;}
  Float_t Pt() const {return fPt;}
  Float_t P() const {return fP;}
  Float_t Dca() const {return fDcaV0;}
  Float_t Ct() const {return fDecayLength;}
  Float_t Cos() const {return fCosPointingAngle;}
  Bool_t  Mc() const {return fMcTruth;}
  Float_t Y() const {return fRapidity;}
  Short_t Charge() const {return fCharge;}
  Short_t Species() const {return fParticleSpecies;}
  Bool_t  OnFlyStatus() const {return fOnFlyStatus;}


private:
  TVector3                fPosition;         //< Lorentzvector of v0
  TVector3                fPvect;            //< Momentum vector of mother
  AliReducedHypTritTrack* fPiTrack;          //< positive daughter of v0
  AliReducedHypTritTrack* fHeTrack;          //< negative daughter of v0
  Float_t                 fP;                //< momentum of mother
  Float_t                 fPt;               //<  transverse momentum of mother
  Float_t                 fM;                //< reconstructed invariant mass of mother
  Float_t                 fDcaV0;            //< DCA of daughters
  Float_t                 fCosPointingAngle; //< cosine of pointing angle of vertex
  Float_t                 fDecayLength;      //< decay radius of mother particle
  Bool_t                  fMcTruth;          //< Monte Carlo truth of mother type
  Float_t                 fRapidity;         //< Rapidity of V0
  Short_t                 fCharge;           //< anti or particle
  Short_t                 fParticleSpecies;  //< particle species
  Bool_t                  fOnFlyStatus;       //< ontheflyStatus


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
  Float_t  Centrality() const {return fCentrality;}
  UShort_t NumberV0s() const {return fNumberV0s;}
  UShort_t Trigger() const {return fTrigger;}
  UShort_t IsMBtriggered() const {return fTrigMB;}
  UShort_t IsHNUtriggered() const {return fTrigHNU;}
  UShort_t IsHQUtriggered() const {return fTrigHQU;}
  UShort_t IsHJTtriggered() const {return fTrigHJT;}
  UShort_t IsHSEtriggered() const {return fTrigHSE;}
  UShort_t IsV0triggered() const {return fTrigV0;}
  UShort_t IsSPDtriggered() const {return fTrigSPD;}
  TString TriggerClasses() const {return fTriggerClasses;}
  Int_t    RunNumber() const {return fRunNumber;}
	Float_t	SPDFiredChips0() const {return fSPDFiredChips0;}
	Float_t SPDFiredChips1() const {return fSPDFiredChips1;}
	Float_t SPDTracklets() const {return fSPDTracklets;}
	Float_t SPDCluster() const {return fSPDCluster;}
	Float_t	V0Multiplicity() const {return fV0Multiplicity;}
	Float_t MultV0M()	const {return fMultV0M;}	
	Float_t MultOfV0M()	const {return fMultOfV0M;}		
	Float_t MultSPDTracklet() const {return fMultSPDTracklet;}	
	Float_t MultSPDCluster() const {return fMultSPDCluster;}
	Float_t MultRef05()	const {return fMultRef05;}		
	Float_t MultRef08()	const {return fMultRef08;}
  void ClearEvent();

private:
  TVector3      fVertexPosition; //< position of primary vertex
  TClonesArray* fV0s;            //< array of v0s in event
  UShort_t      fNumberV0s;      //< number of v0s in event
  Float_t       fCentrality;     //< centrality of event
  Int_t         fRunNumber;      //< number of run
  UShort_t      fTrigger;        //< array of Triggers
  UShort_t      fTrigMB;         //< Flag for MB trigger
  UShort_t      fTrigHNU;        //< Flag for TRD-HNU trigger
  UShort_t      fTrigHQU;        //< Flag for TRD-HQU trigger
  UShort_t      fTrigHJT;        //< Flag for TRD-HJT trigger
  UShort_t      fTrigHSE;        //< Flag for TRD-HSE trigger
  UShort_t      fTrigV0;         //< Flag for HM-V0 trigger
  UShort_t      fTrigSPD;        //< Flag for HM-SPD trigger
  
  TString 	    fTriggerClasses; //< fired trigger classes

	Float_t					fSPDFiredChips0;	// multiplicity triggers
	Float_t					fSPDFiredChips1;
	Float_t 				fSPDTracklets;
	Float_t 				fSPDCluster;
	Float_t					fV0Multiplicity;

	Float_t 				fMultV0M;			// multiplicity estimators
	Float_t 				fMultOfV0M;			
	Float_t 				fMultSPDTracklet;	
	Float_t 				fMultSPDCluster;	
	Float_t 				fMultRef05;			
	Float_t 				fMultRef08;		
		
  AliReducedHypTritEvent(const AliReducedHypTritEvent&);
  AliReducedHypTritEvent &operator = (const AliReducedHypTritEvent&);
  ClassDef(AliReducedHypTritEvent, 6);
};

#endif
