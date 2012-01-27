#ifndef ALIANALYSISTASKCHARGEDHADRONSPECTRA_H
#define ALIANALYSISTASKCHARGEDHADRONSPECTRA_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// This analysis extracts pT-spectra of charged kaons, protons, and pions.  //
// It is based on particles identifation via the dE/dx signal of the TPC.   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

class TH1;
class TH1F;
class TH2F;
class TH3F;
class TList;
class TObjArray;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliHeader;
class AliESDpid;


#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"

class AliAnalysisTaskChargedHadronSpectra : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskChargedHadronSpectra(const char *name);
  AliAnalysisTaskChargedHadronSpectra();
  virtual ~AliAnalysisTaskChargedHadronSpectra() {}
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
  void           SetTrackingMode(Int_t trackingMode = 0){fTrackingMode = trackingMode; Initialize();};
  void           Initialize();
  //
  static TH1D *  AnalyseClassicProton(const TH3F * input, Int_t EtaBinLow,Int_t EtaBinHigh, const Double_t * AlephParams);
  static TH1D *  AnalyseClassicPion(const TH3F * input, Int_t EtaBinLow, Int_t EtaBinHigh, const Double_t * AlephParams);
  static TH1D *  AnalyseClassicKaon(const TH3F * input, Int_t EtaBinLow, Int_t EtaBinHigh, const Double_t * AlephParams);
  //
  static void    Postprocess(const TList * ListOfHistogramsMC,const  TList * ListOfHistogramsData, const Char_t *filename);

 private:
  //
  void  BinLogX(const TH1 *h);
  Int_t GetPythiaEventProcessType(const AliHeader* aHeader, const Bool_t adebug = kFALSE) const;
  //
  AliESDEvent *fESD;                  //! ESD object
  TList       *fListHist;             //! list for histograms
  //
  AliESDtrackCuts * fESDtrackCuts;    // basic cut variables
  AliESDpid       * fESDpid;          // basic TPC object for n-sigma cuts
  Bool_t        fMCtrue;              // flag if real data or MC is processed
  Int_t         fTrackingMode;        // flag which traking mode should be used: 0: TPC only, 1: global tracking
  Double_t      fAlephParameters[5];  // Aleph Parameters for Bethe-Bloch
  //
  // MC histogram
  //
  TH3F        *fHistPtMCKaon;        //! (mult,eta,pT) for Kaons MC truth; neg. x-axis for neg. particles, pos. x-axis for pos. particles
  TH3F        *fHistPtMCProton;      //! (mult,eta,pT) for Protons MC truth; neg. x-axis for neg. particles, pos. x-axis for pos. particles
  TH3F        *fHistPtMCPion;        //! (mult,eta,pT) for Pions MC truth; neg. x-axis for neg. particles, pos. x-axis for pos. particles

  // reconstructed particle histograms
  TH3F        *fHistPtEtaKaon;       //!  (mult,eta,pT) for Kaons; neg. x-axis for neg. particles, pos. x-axis for pos. particles
  TH3F        *fHistPtEtaKaonNoKink; //!  (mult,eta,pT) for Kaons withou accepting the Kink mother; neg. x-axis for neg. particles, pos. x-axis for pos. particles
  TH3F        *fHistPtEtaProton;     //!  (mult,eta,pT) for Protons; neg. x-axis for neg. particles, pos. x-axis for pos. particles
  TH3F        *fHistPtEtaProtonDCA;  //!  (DCA,eta,pT) for Protons; neg. x-axis for neg. particles, pos. x-axis for pos. particles; special histogram for protons
  TH3F        *fHistPtEtaPion;       //!  (mult,eta,pT) for Pions; neg. x-axis for neg. particles, pos. x-axis for pos. particles
  //
  TH3F        *fHistClassicKaon;     //! (Pt,eta,delta dEdx) for Kaons for different eta; neg. x-axis for neg. particles, pos. x-axis for pos. particles
  TH3F        *fHistClassicProton;   //! (Pt,eta,delta dEdx) for Protons for different eta; neg. x-axis for neg. particles, pos. x-axis for pos. particles
  TH3F        *fHistClassicPion;     //! (Pt,eta,delta dEdx) for Pions for different eta; neg. x-axis for neg. particles, pos. x-axis for pos. particles
   
  // histograms of general interest
  TH3F        *fDeDx;                 //! dEdx spectrum
  TH2F        *fHistTrackPerEvent;    //! tracks per event for multiplicity studies; code: (0) all calls; (1) all selected; (2) all selected with vtx. < 10cm; 
  TH3F        *fHistTrackPerEventMC;  //! (code, TrackPerEvent, isSelected) tracks per event, codes represent different event types (non-diffractive,..), is selected according to event selection procedure
  TH2F        *fSecProtons;           //! control histogram for secondary interactions
  TH3F        *fVertexZ;              //! control histogram for the z-position of the vertex
  //
  TH2F        *fHistEtaNcls;          //! 2d histogram (eta, nTPCclusters) which will define our acceptance
  TH2F        *fHistEtaPhi;           //! 2d histogram (eta, phi) which will show dead regions

  // histograms for efficiency studies
  TH3F        *fHistEffProton;       //! (code,eta,pT) special hist. for eff. studies; code 0: true primary p, 1: true sec. p, 2: misidentified, 3: weak decay sec.
  TH3F        *fHistEffProtonDCA;    //! (code,dca,pT) special hist. for eff. studies; code 0: true primary p, 1: true sec. p, 2: misidentified, 3: weak decay sec.
  TH3F        *fHistEffPion;         //! (code,eta,pT) special hist. for eff. studies; code 0: true primary pi, 1: true sec. pi, 2: misidentified, 3: weak decay sec., 4: muons
  TH3F        *fHistEffKaon;         //! (code,eta,pT) special hist. for eff. studies; code 0: true primary K, 1: true sec. K, 2: misidentified, 3: weak decay sec.
  //
  //
  //
  THnSparseS * fHistRealTracks;      //! histogram with all necessary information for real tracks
  THnSparseS * fHistMCparticles;     //! histogram with all necessary information for MC particles


  AliAnalysisTaskChargedHadronSpectra(const AliAnalysisTaskChargedHadronSpectra&); 
  AliAnalysisTaskChargedHadronSpectra& operator=(const AliAnalysisTaskChargedHadronSpectra&); 

  ClassDef(AliAnalysisTaskChargedHadronSpectra, 1); 
};

#endif
