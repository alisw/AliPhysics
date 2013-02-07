#ifndef ALIHFENONPHOTONICELECTRON_H
#define ALIHFENONPHOTONICELECTRON_H

 /************************************************************************************
  *											*
  *	Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.	*
  *	See cxx source for full Copyright notice					*
  *											*
  *											*
  *											*
  *	Task for the Selection of Non-photonic Electron study				*
  *											*
  *			Author: R.Bailhache, C.A.Schmidt				*
  *											*
  ************************************************************************************/

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

class AliESDtrackCuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliMCEvent;
class AliVEvent;
class AliVParticle;
class THnSparse;
class TClonesArray;
class TList;

class AliHFENonPhotonicElectron : public TNamed {
 public:

  typedef enum{
    kElectronfromconversion = 0,
    kElectronfromconversionboth = 1,
    kElectronfrompi0 = 2,
    kElectronfrompi0both = 3,
    kElectronfrometa = 4,
    kElectronfrometaboth = 5,
    kElectronfromC = 6,
    kElectronfromB = 7,
    kElectronfromother = 8,
    kNoElectron = 9
  } Source_t;

  typedef enum{
    kS = 0,
    kOp = 1
  } Sign_t;

  AliHFENonPhotonicElectron();
  AliHFENonPhotonicElectron(const char *name, const Char_t *title);
  AliHFENonPhotonicElectron &operator=(const AliHFENonPhotonicElectron &ref);
  virtual ~AliHFENonPhotonicElectron();

  void SetAOD    		(Bool_t isAOD)     		{ fIsAOD = isAOD; };
  void SetMCEvent		(AliMCEvent *mcEvent)		{ fMCEvent = mcEvent; };
  void SetAODArrayMCInfo	(TClonesArray *aodArrayMCInfo) { fAODArrayMCInfo = aodArrayMCInfo; };
  void SetUseFilterAOD		(Bool_t useFilterAOD)		{ fUseFilterAOD = useFilterAOD; };
  void SetFilter		(UInt_t filter)			{ fFilter = filter; };
  void SetHFEBackgroundCuts	(AliHFEcuts * const cuts)	{ fHFEBackgroundCuts = cuts; };

  AliHFEpid		*GetPIDBackground()		const	{ return fPIDBackground; };
  AliHFEpidQAmanager	*GetPIDBackgroundQAManager()	const	{ return fPIDBackgroundQA; };

  void  SetMaxInvMass		(Double_t MaxInvMass)		{ fMaxInvMass		= MaxInvMass; };
  void  SetMaxOpening3D		(Double_t MaxOpening3D)		{ fMaxOpening3D		= MaxOpening3D; };
//  void  SetMaxOpeningTheta	(Double_t MaxOpeningTheta)	{ fMaxOpeningTheta	= MaxOpeningTheta; };
//  void  SetMaxOpeningPhi	(Double_t MaxOpeningPhi)	{ fMaxOpeningPhi	= MaxOpeningPhi; };
  void  SetAlgorithmMA		(Bool_t algorithmMA)	 	{ fAlgorithmMA		= algorithmMA; };
  void  SetMassConstraint	(Bool_t MassConstraint)		{ fSetMassConstraint	= MassConstraint; };

  TList      *GetListOutput()		const	{ return fListOutput; };
  THnSparseF *GetAssElectronHisto()	const	{ return fAssElectron; };
  THnSparseF *GetIncElectronHisto()	const	{ return fIncElectron; };
  THnSparseF *GetUSignHisto()		const	{ return fUSign; };
  THnSparseF *GetLSignHisto()		const	{ return fLSign; };
//  THnSparseF *GetUSignAngleHisto() const { return fUSignAngle; };
//  THnSparseF *GetLSignAngleHisto() const { return fLSignAngle; };

  void     Init				();
  void     InitRun			(const AliVEvent *inputEvent, const AliPIDResponse *pidResponse);
  Int_t    FillPoolAssociatedTracks	(AliVEvent *inputEvent, Int_t binct=-1);
  Int_t    CountPoolAssociated		(AliVEvent *inputEvent, Int_t binct=-1);
  Int_t    LookAtNonHFE			(Int_t iTrack1, AliVTrack *track1, AliVEvent *vEvent, Double_t weight=1., Int_t binct=-1, Double_t deltaphi=-1, Int_t source=-1, Int_t indexmother=-1);

  Int_t    FindMother		(Int_t tr, Int_t &indexmother);
  Int_t    CheckPdg		(Int_t tr);
  Int_t    IsMotherGamma	(Int_t tr);
  Int_t    IsMotherPi0		(Int_t tr);
  Int_t    IsMotherC		(Int_t tr);
  Int_t    IsMotherB		(Int_t tr);
  Int_t    IsMotherEta		(Int_t tr);


 private:
  Bool_t		 fIsAOD;			// Is AOD
  AliMCEvent		*fMCEvent;			//! MC event ESD
  TClonesArray		*fAODArrayMCInfo;		//! MC info particle AOD
  AliHFEcuts		*fHFEBackgroundCuts;		// HFE background cuts
  AliHFEpid		*fPIDBackground;		// PID background cuts
  AliHFEpidQAmanager	*fPIDBackgroundQA;		// QA Manager Background
  const AliPIDResponse	*fkPIDRespons;			// PID response
  Bool_t		 fAlgorithmMA;			// algorithm MA
  Bool_t		 fUseFilterAOD;			// Use the preselected AOD track
  UInt_t		 fFilter;			// filter AOD status
  Double_t		 fChi2OverNDFCut;		// Limit chi2
  Double_t		 fMaxDCA;			// Limit dca
//  Double_t		 fMaxOpeningTheta;		// Limit opening angle in theta
//  Double_t		 fMaxOpeningPhi;		// Limit opening angle in phi
  Double_t		 fMaxOpening3D;			// Limit opening 3D
  Double_t		 fMaxInvMass;			// Limit invariant mass
  Bool_t		 fSetMassConstraint;		// Set mass constraint
  TArrayI		*fArraytrack;			//! list of associated tracks
  Int_t			 fCounterPoolBackground;	// number of associated electrons
  Int_t			 fnumberfound;			// number of inclusive  electrons
  TList			*fListOutput;			// List of histos
  THnSparseF		*fAssElectron;			//! centrality, pt, Source MC, P, TPCsignal
  THnSparseF		*fIncElectron;			//! centrality, pt, Source MC, P, TPCsignal
  THnSparseF		*fUSign;			//! delta phi, c, pt, inv, source
  THnSparseF		*fLSign;			//! delta phi, c, pt, inv, source
//  THnSparseF		*fUSignAngle;			//! angle, c, source
//  THnSparseF		*fLSignAngle;			//! angle, c, source


  AliHFENonPhotonicElectron(const AliHFENonPhotonicElectron &ref); 

  ClassDef(AliHFENonPhotonicElectron, 1); //!example of analysis
};

#endif
