#ifndef ALIANALYSISTPCTOFPA_H
#define ALIANALYSISTPCTOFPA_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// This analysis extracts pT-spectra of charged kaons, protons, and pions.  //
// It is based on particles identifation via the dE/dx signal of the TPC.   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

class TH1;
class TH1F;
class TH2D;
class TH3D;
class TList;
class TObjArray;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliHeader;
class AliESDpid;


#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"

class AliAnalysisTPCTOFpA : public AliAnalysisTaskSE {
 public:
  AliAnalysisTPCTOFpA(const char *name);
  AliAnalysisTPCTOFpA();
  virtual ~AliAnalysisTPCTOFpA() {}
  //
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  //
  Bool_t         SelectOnImpPar(AliESDtrack* t);
  //
  void           SetESDtrackCuts(AliESDtrackCuts * trackCuts){fESDtrackCuts = trackCuts;};
  void           SetAlephParameters(const Double_t * parameters){for(Int_t j=0;j<5;j++) fAlephParameters[j] = parameters[j]; Initialize();};
  void           SetIsMCtrue(Bool_t isMCdata = kTRUE){fMCtrue = isMCdata;};
  void           SetUseHBTmultiplicity(Bool_t useHBTmultiplicity = kTRUE){fUseHBTmultiplicity = useHBTmultiplicity;};
  void		 SetUseTPConlyTracks(Bool_t useTPConlyTracks = kTRUE){fUseTPConlyTracks = useTPConlyTracks;};
  void           SetSaveMotherPDG(Bool_t saveMotherPDG =kTRUE){fSaveMotherPDG = saveMotherPDG;};
  void           SetSmallTHnSparse(Bool_t smallTHnSparse = kTRUE) {fSmallTHnSparse = smallTHnSparse;};
  void           SetTPCnSigmaCuts(Double_t nSigmaTPCLow = -3., Double_t nSigmaTPCHigh = 3.){fTPCnSigmaCutLow = nSigmaTPCLow; fTPCnSigmaCutHigh = nSigmaTPCHigh;};
  void           SetRapidityCuts(Double_t rapidityLow = 0., Double_t rapidityHigh = 0.5){fRapidityCutLow = rapidityLow; fRapidityCutHigh = rapidityHigh;};
  void           SetEvenDCAbinning(Bool_t EvenDCAbinning = kTRUE) {fEvenDCAbinning = EvenDCAbinning;};
  void           SetIspA(Bool_t ispA = kTRUE) {fIspA = ispA;};
  void           SetRapCMS(Bool_t rapCMS = kTRUE) {fRapCMS = rapCMS;};
  void           SetCentEst(TString centEst = "V0M") {fCentEst = centEst.Data();};
  void           SetTOFmisMatch(Int_t TOFmisMatch = 2) {fTOFmisMatch = TOFmisMatch;};
  void           SetTOFwindow(Double_t TOFwindow = 10.) {fTOFwindow = TOFwindow;};
  void           SetCrossedRows(Double_t crossedRows = 70.) {fCrossedRows = crossedRows;};
  void           SetRatioRowsClusters(Double_t ratioRowsClusters = 0.8) {fRatioRowsClusters = ratioRowsClusters;};
  void           SetTRDinReject(Bool_t TRDinReject = kFALSE) {fTRDinReject = TRDinReject;};
  void           SetDCAzCut(Double_t dcaZcut = 2.){fDCAzCut = dcaZcut;};
  void           Initialize();
  //
  
 private:
  //
  void  BinLogAxis(const TH1 *h);
  Int_t GetPythiaEventProcessType(const AliHeader* aHeader, const Bool_t adebug = kFALSE) const;
  //
  AliESDEvent *fESD;                   //! ESD object
  TList       *fListHist;              //! list for histograms
  //
  AliESDtrackCuts * fESDtrackCuts;     // basic cut variables
  AliESDtrackCuts * fESDTrackCutsMult; // cuts for the MULTIPLICITY DETERMINATION
  AliESDpid       * fESDpid;           // basic TPC object for n-sigma cuts
  Bool_t        fMCtrue;               // flag if real data or MC is processed
  Bool_t        fOnlyQA;               // flag if only QA histograms should be filled
  Bool_t        fUseHBTmultiplicity;   // flag if multiplicity determination should be done as in the HBT paper
  Bool_t	fUseTPConlyTracks;     // flag if TPConly-track should be used
  Bool_t        fSaveMotherPDG;        // flag if PDG of mother should be saved (weak decays)
  Bool_t        fSmallTHnSparse;       // flag if to do cuts on TPC n-sigma and rapidity in task or not
  Bool_t        fIspA;                 // flag for pA analysis                                                               
  Bool_t        fRapCMS;               // flag if rapitidy should be shifted by 0.465 do have rap in CMS of pPb
  TString       fCentEst;              // string which contains the string for the centrality estimator
  Int_t         fTOFmisMatch;          // switch for how tof mismatch should be handled. possible options 0,1,2
  Bool_t        fTRDinReject;          // flag to reject all tracks with TRDin flag set
  Double_t      fTOFwindow;            // set cut on dx and dz TOF window
  Double_t      fDCAzCut;              // set cut on DCA z -standard is 2cm
  Double_t      fCrossedRows;          // min. number of crossed rows for track cuts
  Double_t      fRatioRowsClusters;    // ratio of findable clusters over crossed rows
  Double_t      fTPCnSigmaCutLow;      // low border for TPC n-sigma cut
  Double_t      fTPCnSigmaCutHigh;     // high border for TPC n-sigma cut
  Double_t      fRapidityCutLow;       // low border for rapidity cut
  Double_t      fRapidityCutHigh;      // high border for rapidity cut
  Double_t      fEvenDCAbinning;       // same bin width for all bins in DCA xy
  Double_t      fAlephParameters[5];   // Aleph Parameters for Bethe-Bloch

  //
  //
  //
  THnSparseF * fHistRealTracks;        //! histogram with all necessary information for real tracks
  THnSparseF * fHistMCparticles;       //! histogram with all necessary information for MC particles
  //
  TH3D       * fHistPidQA;             //! histogram for the QA of the PID
  TH2D       * fHistMult;              //! control histogram for multiplicity
  TH1D       * fHistCentrality;        //! control histogram for centrality
  TH2D       * fHistTOFwindow;         //! control histogram for TOF window
  //
  AliAnalysisTPCTOFpA(const AliAnalysisTPCTOFpA&); 
  AliAnalysisTPCTOFpA& operator=(const AliAnalysisTPCTOFpA&); 

  ClassDef(AliAnalysisTPCTOFpA, 1); 
};

#endif
