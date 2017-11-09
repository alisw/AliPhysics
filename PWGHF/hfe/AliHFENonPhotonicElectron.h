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

#ifndef ROOT_TArrayD
#include <TArrayD.h>
#endif

class AliESDtrackCuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliMCEvent;
class AliKFVertex;
class AliVEvent;
class AliVParticle;
class AliVTrack;
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
    kElectronfromomega = 6,
    kElectronfromomegaboth = 7,
    kElectronfromC = 8,
    kElectronfromB = 9,
    kElectronfromother = 10,
    kNoElectron = 11
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
  void SetMCEvent		(AliMCEvent *mcEvent);
  void SetAODArrayMCInfo	(TClonesArray *aodArrayMCInfo);
  void SetHFEBackgroundCuts	(AliHFEcuts * const cuts)	{ fHFEBackgroundCuts = cuts; };
  void SetWithWeights(Int_t levelBack);

  AliHFEpid		*GetPIDBackground()		const	{ return fPIDBackground; };
  AliHFEpidQAmanager	*GetPIDBackgroundQAManager()	const	{ return fPIDBackgroundQA; };
  Int_t GetLevelBack()                                  const   { return fLevelBack; };

  void  SetMaxInvMass		(Double_t MaxInvMass)		{ fMaxInvMass		= MaxInvMass; };
  void  SetMaxOpening3D		(Double_t MaxOpening3D)		{ fMaxOpening3D		= MaxOpening3D; };
//  void  SetMaxOpeningTheta	(Double_t MaxOpeningTheta)	{ fMaxOpeningTheta	= MaxOpeningTheta; };
//  void  SetMaxOpeningPhi	(Double_t MaxOpeningPhi)	{ fMaxOpeningPhi	= MaxOpeningPhi; };
  void  SetStudyRadius		(Bool_t studyRadius)	 	{ fStudyRadius		= studyRadius; };
  void  SetAlgorithmMA		(Bool_t algorithmMA)	 	{ fAlgorithmMA		= algorithmMA; };
  void  SetMassConstraint	(Bool_t MassConstraint)		{ fSetMassConstraint	= MassConstraint; };
  void  SetITSMeanShift         (Double_t meanshift)            { fITSmeanShift = meanshift; }
  void  SetITSnSigmaHigh        (Double_t nSigmaHigh)           { fITSnSigmaHigh = nSigmaHigh; }
  void  SetITSnSigmaLow         (Double_t nSigmaLow)            { fITSnSigmaLow = nSigmaLow; }
  void  SetminPt             (Double_t minpt)                   { fminPt = minpt; }
  void  SetEtaDalitzWeightFactor(Double_t etaDalitzWeightFactor){ fEtaDalitzWeightFactor = etaDalitzWeightFactor;}

  void SelectCategory1Tracks(Bool_t doSelect = kTRUE)           { fSelectCategory1tracks = doSelect; }
  void SelectCategory2Tracks(Bool_t doSelect = kTRUE)           { fSelectCategory2tracks = doSelect; }

  void SetAnaPairGen(Bool_t setAna = kTRUE, Int_t nGen = 2)     { fAnaPairGen = setAna; fNumberofGenerations = nGen;};
  void SetNPairGenerations(Int_t nGen)                          { fNumberofGenerations = nGen;};
  void SetDisplayMCStack(Bool_t setDisplay = kTRUE)             { fDisplayMCStack = setDisplay;};

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
  Int_t    LookAtNonHFE			(Int_t iTrack1, AliVTrack *track1, AliVEvent *vEvent, Double_t weight=1., Int_t binct=-1, Double_t deltaphi=-1, Int_t source=-1, Int_t indexmother=-1,Int_t mcQAsource=-1);

  Int_t    FindMother		(Int_t tr, Int_t &indexmother) const;

  void SetPtBinning(const TArrayD &binning) { fPtBinning = binning; }
  void SetPtBinning(Int_t nbins, const Double_t *const binning) { fPtBinning.Set(nbins+1, binning); }
  void SetEtaBinning(const TArrayD &binning) { fEtaBinning = binning; }
  void SetEtaBinning(Int_t nbins, const Double_t *const binning) { fEtaBinning.Set(nbins+1, binning); }
    void SetPhiBinning(const TArrayD &binning) { fPhiBinning = binning; }
    void SetPhiBinning(Int_t nbins, const Double_t *const binning) { fPhiBinning.Set(nbins+1, binning); }
  void SetInvMassBinning(const TArrayD &binning) { fInvMassBinning = binning; }
  void SetInvMassBinning(Int_t nbins, const Double_t *const binning) { fInvMassBinning.Set(nbins+1, binning); }


 private:
  void	   FillMotherArray(Int_t tr, int index, Int_t a[], int NumberofGenerations);
  Int_t    FindGeneration(Int_t a[], Int_t b[], int NumberofGenerations); 
  Int_t    GetMotherPDG(Int_t tr, Int_t &motherIndex) const;
  Int_t    CheckPdg		(Int_t tr) const;
  Double_t Radius               (Int_t tr) const;
  Int_t    IsMotherGamma	(Int_t tr) const;
  Int_t    IsMotherPi0		(Int_t tr) const;
  Int_t    IsMotherC		(Int_t tr) const;
  Int_t    IsMotherB		(Int_t tr) const;
  Int_t    IsMotherEta		(Int_t tr) const;
  Int_t    IsMotherOmega	(Int_t tr) const;
  Bool_t MakePairDCA(const AliVTrack *inclusive, const AliVTrack *associated, AliVEvent *vEvent, Bool_t isAOD, Double_t &invMass, Double_t &angle) const;
  Bool_t MakePairKF(const AliVTrack *inclusive, const AliVTrack *associated, AliKFVertex &primV, Double_t &invMass, Double_t &angle) const;
  Bool_t FilterCategory1Track(const AliVTrack * const track, Bool_t isAOD, Int_t binct);
  Bool_t FilterCategory2Track(const AliVTrack * const track, Bool_t isAOD);

  Bool_t                    fIsAOD;                         // Is AOD
  AliMCEvent                *fMCEvent;                      //! MC event ESD
  TClonesArray              *fAODArrayMCInfo;               //! MC info particle AOD
  Int_t                     fLevelBack;                     // Level Background
  AliHFEcuts                *fHFEBackgroundCuts;            // HFE background cuts
  AliHFEpid                 *fPIDBackground;                // PID background cuts
  AliHFEpidQAmanager        *fPIDBackgroundQA;              // QA Manager Background
  const AliPIDResponse      *fkPIDRespons;                  // PID response
  TArrayD                   fPtBinning;                     // pt binning
  TArrayD                   fEtaBinning;                    // eta binning
    TArrayD                   fPhiBinning;                  // phi binning
  TArrayD                   fInvMassBinning;                // Inv mass binning
  Bool_t                    fStudyRadius;                   // Study radius
  Bool_t                    fAlgorithmMA;                   // algorithm MA
  Double_t                  fChi2OverNDFCut;                // Limit chi2
  Double_t                  fMaxDCA;                        // Limit dca
//  Double_t                fMaxOpeningTheta;               // Limit opening angle in theta
//  Double_t                fMaxOpeningPhi;                 // Limit opening angle in phi
  Double_t                  fMaxOpening3D;                  // Limit opening 3D
  Double_t                  fMaxInvMass;                    // Limit invariant mass
  Bool_t                    fSetMassConstraint;             // Set mass constraint
  Bool_t                    fSelectCategory1tracks;         // Category 1 tracks: Standard track cuts
  Bool_t                    fSelectCategory2tracks;         // Category 2 tracks: tracks below 300 MeV/c
  Double_t                  fITSmeanShift;                  // Shift of the mean in the ITS
  Double_t                  fITSnSigmaHigh;                 // ITS n Sigma electron cut high (>0)
  Double_t                  fITSnSigmaLow;                  // ITS n Sigma electron cut low (<0)
  Double_t                  fminPt;                         // min pT cut for the associated leg
  Double_t                  fEtaDalitzWeightFactor;         // Relative modification for the weighting factor for electrons from Eta Dalitz decays (default = 1);
  TArrayI                   *fArraytrack;                   //! list of associated tracks
  Int_t                     fCounterPoolBackground;         // number of associated electrons
  Int_t                     fnumberfound;                   // number of inclusive  electrons
  TList                     *fListOutput;                   // List of histos
  THnSparseF                *fAssElectron;                  //! centrality, pt, Source MC, P, TPCsignal
  THnSparseF                *fIncElectron;                  //! centrality, pt, Source MC, eta, phi, charge
  THnSparseF                *fUSign;                        //! delta phi, c, pt, inv. mass, source, opening angle, pt_assoc, eta_inc, eta_assoc, phi_inc, charge_inc
  THnSparseF                *fLSign;                        //! delta phi, c, pt, inv. mass, source, opening angle, pt_assoc, eta_inc, eta_assoc, phi_inc, charge_inc
  THnSparseF                *fUSmatches;                    //! number of matched tracks with oposite sign per inclusive track after inv mass cut
  THnSparseF                *fLSmatches;                    //! number of matched tracks with same sign per inclusive track after inv mass cut
  TH2F                      *fHnsigmaITS;                    //! Control histogram for ITS pid of category 2 tracks
  TH2F                      *fWeightsSource;                 //! Control histo for sources for weights  

  THnSparseF                *fIncElectronRadius;            //! For fakes
  THnSparseF                *fRecElectronRadius;            //! For fakes                    
//  THnSparseF              *fUSignAngle;                   //! angle, c, source
//  THnSparseF              *fLSignAngle;                   //! angle, c, source

  Bool_t                    fAnaPairGen;                     // switch on the analysis of the pair generation (switch for performance)
  Int_t                     fNumberofGenerations;            // number of generations stored in pair container variable nGen
  Bool_t                    fDisplayMCStack;                 // display MC stack for true likesign pairs (usually misidentification), for debugging

  AliHFENonPhotonicElectron(const AliHFENonPhotonicElectron &ref); 

  ClassDef(AliHFENonPhotonicElectron, 5); //!example of analysis
};

#endif
