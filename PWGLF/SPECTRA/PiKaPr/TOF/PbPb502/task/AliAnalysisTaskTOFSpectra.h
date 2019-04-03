#ifndef AliAnalysisTaskTOFSpectra_H
#define AliAnalysisTaskTOFSpectra_H

///////////////////////////////////////////////////////////////////////////////
///                                                                          //
/// This analysis extracts pT-spectra of charged kaons, protons, and pions.  //
/// It is based on particles identifation via the TOF signal.                //
///                                                                          //
///                                                                          //
/// Authors:                                                                 //
/// N. Jacazio,  nicolo.jacazio[AROBASe]bo.infn.it                           //
///////////////////////////////////////////////////////////////////////////////

class TH1;
class TH1F;
class TH2F;
class TH3F;
class TH3I;
class TH3S;
class TList;
class TObjArray;
class AliESDEvent;
class AliESDtrack;
class AliMCEvent;
class AliHeader;
class AliESDpid;
class AliAnalysisFilter;
class AliCFContainer;
class TDatabasePDG;
class AliESDVertex;
class AliAnTOFevent;
class TProfile;
class TTree;
class AliTOFT0maker;
class AliTOFcalib;
class AliPIDResponse;
class AliESDTOFCluster;

//Includes
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "AliMultSelection.h"
#include "AliTOFPIDResponse.h"
#include "AliUtilTOFParams.h"
#include "TBenchmark.h"

using namespace AliUtilTOFParams;

// #define TESTFLAG//Flag intendended for testing
#define USETREECLASS //Flag temporary to use a class to be stored in the tree
//#define CHECKPERFORMANCE//Flag to check the time usage

// #define CHECKTRACKCUTS//Flag to build advanced histograms checking track cuts

// #define CHECKCOMPUTEDVALUES//Flag to compute the TOF values and compare to the ones expected

// #define BUILDTOFDISTRIBUTIONS// Flag to prepare distributions of T-Texp-T0 for analysis with and without mismatch with TPC information

// #define BUILDT0PLOTS// Flag to prepare distributions of T0 for analysis of the correlation with the vertex position

class AliAnalysisTaskTOFSpectra : public AliAnalysisTaskSE {
  public:
  //Constructors and destructor
  AliAnalysisTaskTOFSpectra(const TString taskname = "TaskTOFChargedHadron", const Int_t collsys = AliUtilTOFParams::kPbPb, Bool_t mc = kFALSE, Bool_t tree = kTRUE, Bool_t chan = kFALSE, Bool_t cuts = kFALSE, Int_t simplecuts = -1);
  virtual ~AliAnalysisTaskTOFSpectra();

  //////////////////////////////
  // Soft configuration flags //
  //////////////////////////////
  Bool_t fBuilTPCTOF;     ///<  Flag to build the TPC TOF separation
  Bool_t fBuilDCAchi2;    ///<  Flag to build the DCAxy distributions with the cut on the Golden Chi2
  Bool_t fUseTPCShift;    ///<  Flag to use the Shift of the TPC nsigma
  Bool_t fPerformance;    ///<  Flag to fill the performance plots
  Bool_t fMCPerformance;  ///<  Flag to fill the MC performance plots
  Bool_t fRecalibrateTOF; ///<  Flag to require to recalibrate the TOF signal
  Bool_t fCutOnMCImpact;  ///<  Flag to cut the MC on the impact parameter instead of the centrality. This is intendend for the generation of MC predictions
  Bool_t fFineTOFReso;    ///<  Flag to compute a finer TOF resolution as a function of the number of tracks with TOF signal
  Bool_t fFineEfficiency; ///<  Flag to use 3D  histograms with the MC information as a function of pT, eta, phi

  //Standard AnalysisTask functions
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*);

  //Output functions
  ///
  /// Defines all the output objects
  void DefineAllTheOutput()
  {
    DefineOutput(1, TList::Class());
    if (fTreemode) {
      DefineOutput(2, TTree::Class());
    }
  };

  ///
  /// Posts all the data to the output objects
  void PostAllTheData()
  {
    PostData(1, fListHist);
    if (fTreemode) {
      PostData(2, fTreeTrack);
    }
  };

  //Initialization methods

  ///
  /// Method to initalize all the variables to zero
  void Init();

  ///
  /// Method to initalize all the variables that need a reset for each track
  void InitializeTrackVar();

  ///
  /// Method to initalize all the variables that need a reset for each MC track
  void InitializeMCTrackVar();

  ///
  /// Method to initalize all the variables that need a reset for each event
  void InitializeEventVar();

  ///
  /// Method to print out the configuration flags of the task
  void PrintStatus();

  ///
  /// Method to obtain the output list for histograms, this is intended so as to check the data size of output
  TList* GetOutput() { return fListHist; }

  //Utility methods
  ///
  ///  Method to get the event Multiplicity
  void ComputeEvtMultiplicity();

  ///
  ///  Method to get the event Multiplicity (binned)
  void ComputeEvtMultiplicityBin();

  ///
  /// Method to get the information of the MC truth of the generated event
  void GatherEventMCInfo();

  ///
  /// Method to get the information of the MC truth before the event is selected
  void AnalyseMCParticles();

  ///
  /// Method to get the information of the MC truth of the reconstructed tracks and fill the histograms
  Bool_t AnalyseMCTracks();

  ///
  /// Method to get the information of the MC truth of the reconstructed tracks
  Bool_t GatherTrackMCInfo(const AliESDtrack* trk);

  ///
  /// Method to decide if the event selection is selected or not
  Bool_t SelectEvents(Int_t& binstart);

  ///
  /// Method to obtain in a common way the vertex information, this also checks the vertex quality
  const AliESDVertex* ObtainVertex();

  // TOF calibration methods

  ///
  /// Method to prepare to run with a new TOF calibratin, this initialize from the run, therefore should only once for each run
  Bool_t TOFCalibInitRun();

  ///
  /// Method to prepare to run with a new TOF calibratin, this initialize from the event, therefore should once for each event
  Bool_t TOFCalibInitEvent();

  //
  //Mask
  //

  //Indexes to store event masks
  enum fEvtMaskIndex { kIsNewEvent,
    kWithVtx,
    kWithVtxCont,
    kVtxInRange,
    kVtxRes,
    kGoodZVtx,
    kPassPileUp,
    kMCVtxInRange,
    kLimitfEvtMask }; //Event selection bitmask fEvtMask
  enum fMCTrkMaskIndex { kIsNegative,
    kIsPion,
    kIsKaon,
    kIsProton,
    kIsPhysicalPrimary,
    kIsFromStrangeness,
    kIsFromMaterial,
    kLimitfMCTrkMask }; //MC information fMCTrkMask

  //Implementation of mask read/write
  void SetEvtMaskBit(fEvtMaskIndex bit, Bool_t value);
  void ResetEvtMaskBit() { ResetMask(fEvtMask); };

  void SetTrkMaskBit(fTrkMaskIndex bit, Bool_t value);
  void ResetTrkMaskBit() { ResetMask(fTrkMask); };

  void SetTPCPIDMaskBit(fTPCPIDMaskIndex bit, Bool_t value);
  void ResetTPCPIDMaskBit() { ResetMask(fTPCPIDMask); };

  void SetTrkCutMaskBit(fTrkCutMaskIndex bit, Bool_t value);
  void ResetTrkCutMaskBit() { ResetMask(fTrkCutMask); };

  void SetMCTrkMaskBit(fMCTrkMaskIndex bit, Bool_t value);
  void ResetMCTrkMaskBit() { ResetMask(fMCTrkMask); };

  //Particle kinematics
  ///
  ///Converts Transverse momentum into linear momentum
  Double_t Pt2P(Double_t pt, Double_t eta)
  {
    if (pt == 0 && eta == 0)
      return -999;
    Double_t tnh = TMath::TanH(eta);
    return TMath::Sqrt(pt * pt / (1 - tnh * tnh));
  }

  ///
  /// Find the index of the pt bin of the track
  void FindPtBin();

  ///
  /// Method to compute the expected time of a particle from the integrated track length
  Double_t ComputeExpectedTime(Double_t mass, Double_t momentum) { return fLength / CSPEED * TMath::Sqrt(momentum * momentum + mass * mass) / momentum; }

  ///
  /// Method to compute the expected momentum of a particle from the expected time
  Double_t ComputeExpectedMomentum(Double_t mass, Double_t time) { return mass * fLength / TMath::Sqrt(CSPEED * CSPEED * time * time - fLength * fLength); }

  ///
  /// Computes the Y for the given particle mass!
  Double_t ComputeY(Double_t mass) { return TMath::ASinH(fPt / TMath::Sqrt(mass * mass + fPt * fPt) * TMath::SinH(fEta)); }

  ///
  /// Computes the Rapidity of the track in the 3 mass hypothesis: pions kaons and protons
  void ComputeRapidity()
  {
    for (Int_t species = 0; species < 3; species++)
      fRapidity[species] = ComputeY(AliPID::ParticleMass(species + AliPID::kPion));
  }

  ///
  /// Method to fill the TOF per channel histogram
  void RunTOFChannel();

  ///
  /// Method to check which are the track cuts of the cut variation parameters which are passed
  Bool_t AnalyseCutVariation(const AliESDtrack* track);

  ///
  /// Method to fill the histograms containing the values of the variables under cut, the flag is to fill the ones before and after the flag
  void FillCutVariable(const Bool_t pass = kFALSE);

  ///
  /// Method to start the time usage
  void StartTimePerformance(const UInt_t test)
  {
#ifdef CHECKPERFORMANCE
    tr[test] += tb.GetRealTime(tests[test]);
    tc[test] += tb.GetCpuTime(tests[test]);
    tb.Start(tests[test]);
#endif
  }

  ///
  /// Method to start the time usage
  void StopTimePerformance(const UInt_t test)
  {
#ifdef CHECKPERFORMANCE
    tb.Stop(tests[test]);
    tr[test] = tb.GetRealTime(tests[test]) - tr[test];
    tc[test] = tb.GetCpuTime(tests[test]) - tc[test];
#endif
  }

  ///
  /// Method to check the time usage
  void FillTimePerformance()
  {
#ifdef CHECKPERFORMANCE
    for (Int_t i = 0; i < ntests; i++) {
      // hPerformanceTime->Fill(i, tb.GetRealTime(tests[i]));
      // hPerformanceCPUTime->Fill(i, tb.GetCpuTime(tests[i]));
      hPerformanceTime->Fill(i, tr[i]);
      hPerformanceCPUTime->Fill(i, tc[i]);
    }
#endif
  }

  ///
  /// Method to define the histograms the time usage
  void DefineTimePerformance()
  {
#ifdef CHECKPERFORMANCE
    hPerformanceTime = new TH1F("hPerformanceTime", "Real Time", ntests, -.5 + 0, -.5 + ntests);
    fListHist->AddLast(hPerformanceTime);

    hPerformanceCPUTime = new TH1F("hPerformanceCPUTime", "CPU Time", ntests, -.5 + 0, -.5 + ntests);
    fListHist->AddLast(hPerformanceCPUTime);

    for (Int_t i = 1; i <= ntests; i++) {
      hPerformanceTime->GetXaxis()->SetBinLabel(i, tests[i - 1]);
      hPerformanceCPUTime->GetXaxis()->SetBinLabel(i, tests[i - 1]);
    }
#endif
  }

  ///
  /// Method to define the histograms used for the Performance
  void DefinePerformanceHistograms();

  ///
  /// Method to fill the histograms used for the Performance
  void FillPerformanceHistograms(const AliVTrack* track);

  ///
  /// Method to define the histograms used for the MC Performance
  void DefineMCPerformanceHistograms();

  ///
  /// Method to fill the histograms used for the MC Performance
  void FillMCPerformanceHistograms(const AliVTrack* track);

  ///
  /// Method to define the histograms used for the T0 information
  void DefineT0Histograms();

  ///
  /// Method to fill the histograms used for the T0 information
  Bool_t FillT0Histograms();

  //Setter methods
  void SetHeavyIonFlag(CollSys mode) { fCollSysMode = mode; };
  void SetMCFlag(Bool_t mode = kTRUE) { fMCmode = mode; };
  void SetTreeFlag(Bool_t mode = kTRUE) { fTreemode = mode; };
  void SetChannelFlag(Bool_t mode = kTRUE) { fChannelmode = mode; };
  void SetCutFlag(Bool_t mode = kTRUE) { fCutmode = mode; };
  void SetTrigger(UInt_t trg = AliVEvent::kMB) { fSelectBit = trg; };
  void SetMultiplicityBinning(TArrayD bin);

  //
  //Single cuts
  void SetTPCRowsCut(Double_t cut);
  void SetTrkChi2Cut(Double_t cut);
  void SetTrkChi2CutITS(Double_t cut);
  void SetDCAxyCut(Double_t cut);
  void SetDCAzCut(Double_t cut);
  void SetGeoCut(Double_t fDeadZoneWidth, Double_t fCutGeoNcrNclLength, Double_t fCutGeoNcrNclGeom1Pt, Double_t fCutGeoNcrNclFractionNcr, Double_t fCutGeoNcrNclFractionNcl);
  void SetRatioCrossedRowsFindableCls(Double_t cut);

  ///
  /// Method to initialize the standard cuts
  void SetTrackCuts();

  ///
  /// Method to set the simple cut variation via the index
  void SetSimpleCutVar();

  ///
  /// Method to set the tree cut variation
  void SetCutVar();

  ///
  /// Method to print the cut variables for the whole list
  void PrintCutVariables();

  ///
  /// Method to print the cut variables for the tree only
  void PrintCutVariablesForTree();

  ///
  ///Bins numbers
  enum { kEtaBins = 36,
    kPhiBins = 18,
    kEvtMultBins = 12,
    kChannelBins = 3276,
    kTimeBins = 10000 };

  //Standard track and event cuts
  AliEventCuts fEventCut;            //!<! basic cut variables for events
  AliESDtrackCuts* fESDtrackCuts;    /// basic cut variables
  AliESDtrackCuts* fESDtrackCutsPrm; /// basic cut variables for primaries

  private:
  /////////////////////
  // Utility objects //
  /////////////////////
  AliESDEvent* fESD;                  //!<! ESD object
  AliMCEvent* fMCEvt;                 //!<! MC event
  AliESDtrackCuts* fCutVar[nCutVars]; //!<! basic cut variables cut variations
  AliMultSelection* fMultSel;         //!<! Multiplicity selection

  //TOF specific objects
  AliESDTOFCluster* fTOFcls;  //!<! TOF cluster object
  AliTOFcalib* fTOFcalib;     //!<! TOF calibration object
  AliTOFT0maker* fTOFT0maker; //!<! TOF T0 maker object

  //PID utilities
  AliPIDResponse* fPIDResponse;      //!<! PID response object
  AliTOFPIDResponse fTOFPIDResponse; //!<! TOF PID response object
  AliESDpid* fESDpid;                //!<! ESD PID

  //TOF specific parameters
  const Double_t fTimeResolution; ///  TOF time resolutions expected

  ///////////////////////
  // Output containers //
  ///////////////////////
  TList* fListHist;           //!<! TList for histograms
  TTree* fTreeTrack;          //!<! TTree to store information of the single track
  AliAnTOFevent* fAnTOFevent; //!<! Event container containing all tracks from one event in the format of AliAnTOFtrack
  TTree* fTreeTrackMC;        //!<! TTree to store MC information of the single track

  /////////////////////////
  // Configuration Flags //
  /////////////////////////
  CollSys fCollSysMode;       ///<  Flag for the Heavy Ion Mode
  Bool_t fMCmode;             ///<  Flag for the Monte Carlo Mode
  Bool_t fTreemode;           ///<  Flag for the Tree analysis Mode
  Bool_t fChannelmode;        ///<  Flag to set the analysis only on channel TOF
  Bool_t fCutmode;            ///<  Flag to set the cut variation mode, cuts are not the standard cuts but are modified accordingly to the requirements
  const Int_t fSimpleCutmode; ///<  Index to set simple configuration of the track cuts
  const Bool_t fUseAliEveCut; ///<  Index to set the usage of the AliEventCuts class from OADB to select events
  const Float_t fDCAXYshift;  ///<  Shift to the DCAxy of the track to compensate for the bias, to be used only in MC, to be defined in the constructor
  UInt_t fSelectBit;          ///<  Mask for Trigger selection

  /////////////////////
  // Event Variables //
  /////////////////////

  //Event Cut values
  const Double_t fVtxZCut; ///<  Max Z displacement of the vertex position

  //Event variables
  Int_t fRunNumber; ///<  Run number under analysis

  //Event flags
  Bool_t fEvtPhysSelected;   ///<  Event is selected by the Physics Selection
  Bool_t fEvtSelected;       ///<  Event is selected by the Event Selection
  Bool_t fEvtMCSampSelected; ///<  Event is selected by the MC Sample Selection

  //Vertex Info
  Double_t fPrimVertex[3]; ///<  X,Y,Z position of the primary vertex
  Int_t fNContrPrimVertex; ///<  Contributors to the primary vertex
  Short_t fVertStatus;     ///<  vertex status (0) vertex does not exist, (1) not enough Contributors, (2) vertex exists but outside the fiducial volume, (3) Good vertex

  //Multiplicity
  Float_t fEvtMult;         ///<  Event Multiplicity
  Short_t fEvtMultBin;      ///<  Event Multiplicity bin to compact the information in the Tree
  TArrayD fMultiplicityBin; ///<  Array of the Event Multiplicity bins

  /////////////////////
  // Track Variables //
  /////////////////////

  //Track Cut values
  const Double_t fTOFmax;      ///<  Max TOF time mesured for tracks
  const Double_t fTOFmin;      ///<  Min TOF time mesured for tracks
  const Double_t fLengthmin;   ///<  Min length for tracks
  const Double_t fRapidityCut; ///<  Max rapidity for tracks
  const Double_t fEtaRange;    ///<  Max pseudorapidity for tracks coincides with the range in eta of the histograms

  //Track flags
  ///
  ///Method to set the flags of the track
  void SetTrackFlags(const AliESDtrack* track);
  Bool_t fTOFout;         ///<  Track has the TOFout
  Bool_t fTRDout;         ///<  Track has the TRDout
  Bool_t fTime;           ///<  Track has the Time
  Bool_t fSign;           ///<  Track charge sign kTRUE (1) for negative and kFALSE (0) for positive
  Bool_t fMismatch;       ///<  Track is flagged as Mismatch (i.e. track in the TPC wrongly matched to TOF time)
  Bool_t fPassGoldenChi2; ///<  Track pass the Golden Chi2 cut
  Bool_t fPassDCAxy;      ///<  Track pass the DCAxy cut
  Bool_t fPassDCAz;       ///<  Track pass the DCAz cut
  Bool_t fITSTPCMatch;    ///<  Track is matched to the TPC

  //Track values
  Int_t fBinPtIndex; ///<  Index of the Pt bin for the track

  ///
  ///Method to set the values of the track
  void SetTrackValues(const AliESDtrack* track, const AliESDVertex* vertex);
  UShort_t fTPCClusters;      ///<  Track clusters in the TPC
  UShort_t fTPCFindClusters;  ///<  Track findable clusters in the TPC
  UChar_t fITSClusters;       ///<  Track clusters in the ITS
  Float_t fTPCCrossedRows;    ///<  Track crossed rows in the TPC
  Double_t fTPCChi2;          ///<  Track chi2 in TPC
  Double_t fTPCChi2PerNDF;    ///<  Track chi2 per NDF in TPC
  Double_t fITSChi2;          ///<  Track chi2 in ITS
  Double_t fITSChi2PerNDF;    ///<  Track chi2 per NDF in ITS
  Double_t fGoldenChi2;       ///<  Track chi2 constrained to the vertex
  Double_t fLengthActiveZone; ///<  Track length in active zone
  Double_t fLength;           ///<  Track length
  Float_t fLengthRatio;       ///<  Track length ratio between the matched track length and the one coming from another cluster
  Float_t fEta;               ///<  Eta of the track
  Double_t fRapidity[3];      ///<  Rapidity of the track in the three mass hypothesis
  Float_t fPhi;               ///<  Phi of the track
  Double_t fP;                ///<  Momentum of the track
  Double_t fPTPC;             ///<  Momentum of the track in the TPC
  Double_t fPt;               ///<  Transverse momentum of the track
  Float_t fDCAXY;             ///<  DCAxy of the track
  Float_t fDCAZ;              ///<  DCAz of the track

  //TOF values
  Double_t fTOFTime;                        ///<  TOF signal
  Float_t fTOFMismatchTime;                 ///<  TOF signal from unmatched tracks
  Float_t fTOFSignalTot;                    ///<  TOF total signal
  Double_t fTOFImpactDZ;                    ///<  Local difference along z of track's impact on the TOF pad and the extrapolated track from the TPC
  Double_t fTOFImpactDX;                    ///<  Local difference along z of track's impact on the TOF pad and the extrapolated track from the TPC
  Int_t fTOFchan;                           ///<  Channel Index of the TOF Signal
  Float_t fT0TrkTime;                       ///<  Best start time of the track
  Int_t fT0UsedMask;                        ///<  Mask with the T0 used (0x1=T0-TOF,0x2=T0A,0x3=TOC) for p bins
  Float_t fT0TrkSigma;                      ///<  Measured resolution on the T0
  Float_t fTOFExpSigma[kExpSpecies];        ///<  TOF expected Sigma of the track in the hypothesis (0) Electron, (1) Muon, (2) Pion, (3) Kaon, (4) Proton
  Float_t fTOFExpTime[kExpSpecies];         ///<  TOF expected time of the track in the hypothesis (0) Electron, (1) Muon, (2) Pion, (3) Kaon, (4) Proton
  Float_t fTOFSigma[kExpSpecies];           ///<  Measured Number of Sigmas of the TOF signal in the hypothesis (0) Electron, (1) Muon, (2) Pion, (3) Kaon, (4) Proton
  Double_t fTOFPIDProbability[kExpSpecies]; ///<  TOF PID Probability of the track
  Int_t fNTOFClusters;                      ///<  Number of TOF clusters matched to track
  Int_t* fTOFClusters;                      //!<! TOF clusters matched to the track

  //TPC values
  Float_t fTPCSignal;             ///<  TPC signal
  Float_t fTPCSigma[kExpSpecies]; ///<  Measured Sigma of the TPC signal in the three hypothesis
  Double_t fPhiout;               ///<  Phi position of the external Track
  Double_t fXout;                 ///<  X position of the external Track
  Double_t fYout;                 ///<  Y position of the external Track
  Double_t fZout;                 ///<  Z position of the external Track
  Double_t fTPCShift[3][kPtBins]; ///<  Shift of the TPC sigmas in order to be centered to zero

  //Combined PID values                       ///<
  Float_t fCombinedSigma[3]; ///<  Measured Sigma of the combined TOF and TPC signal in the three hypothesis

  //MC info
  Double_t fRapidityMC;    ///<  Rapidity of the track with MC truth
  Double_t fPtMC;          ///<  Transverse momentum of the track with MC truth
  Double_t fPMC;           ///<  Momentum of the track with MC truth
  Bool_t fSignMC;          ///<  MC sign of the track kFAlSE for Positive, kTRUE for Negative
  Double_t fPhiMC;         ///<  Phi information of the track with MC truth
  Double_t fEtaMC;         ///<  eta information of the track with MC truth
  Short_t fProdInfo;       ///<  Information on the origin of the particle, if it is Physical Primary (0), Decay from strangeness (1), or Material (2)
  Int_t fPdgcode;          ///<  PID identity of the track with MC truth -> Real PDG code
  Short_t fPdgIndex;       ///<  PID identity of the track with MC truth -> Just the index (0) Pion, (1) Kaon, (2) Proton
  Int_t fPdgcodeMother;    ///<  PID identity of the track Mother with MC truth -> Real PDG code
  Int_t fFirstMotherIndex; ///<  Index of the mother in the stack the track with MC truth
  Int_t fNMCTracks;        ///<  Number of MC tracks in the event
  Int_t fMCPrimaries;      ///<  Number of MC primary tracks in the event
  Int_t fMCTOFMatch;       ///<  Flag to see if the track has a real match in the TOF detector -> (0) Match, (1) Mismatch, (2) Mother label is compatible with match, (-1) Not matched

  //Masks to store the event and track information in the Tree
  UChar_t fEvtMask;     ///< Event information mask
  UShort_t fTrkMask;    ///< Track information mask
  UChar_t fTPCPIDMask;  ///< TPC PID information mask
  UShort_t fTrkCutMask; ///< Mask with the information for the cut variation
  UChar_t fMCTrkMask;   ///< MC information mask
#ifdef TESTFLAG
  static_assert(kLimitfEvtMask <= 8 * sizeof(fEvtMask), "Wrong dimension of fEvtMask");
  static_assert(kLimitfTrkMask <= 8 * sizeof(fTrkMask), "Wrong dimension of fTrkMask");
  static_assert(kLimitfTPCPIDMask <= 8 * sizeof(fTPCPIDMask), "Wrong dimension of fTPCPIDMask");
  static_assert(kLimitfTrkCutMask <= 8 * sizeof(fTrkCutMask), "Wrong dimension of fTrkCutMask");
  static_assert(kLimitfMCTrkMask <= 8 * sizeof(fMCTrkMask), "Wrong dimension of fMCTrkMask");
#endif

  ////////////////
  // Histograms //
  ////////////////

  //Bin and range definition
  //157248 number of Channels
  //48 * 2 number of Channels per strip
  //1638 number of full strips
  //3276 number of half strips
  //Binning for TOF Channel
  const Double_t fChannelFirst; ///<  First TOF Channel
  const Double_t fChannelLast;  ///<  Last TOF Channel

  //Performance
  enum { ntests = 10 };
  const TString tests[ntests] = {
    "InitializeEventVar",
    "Multiplicity",
    "SelectEvents",
    "ComputeEvtMultiplicityBin",
    "TrackLoop",
    "FillTree",
    "ClearTree",
    "EvtFlags",
    "EvtVariables",
    "EvtMasks"

  };
  Double_t tc[ntests];
  Double_t tr[ntests];
  TBenchmark tb;             ///<  TBenchmark
  TH1F* hPerformanceTime;    ///<  Histogram with the Time used
  TH1F* hPerformanceCPUTime; ///<  Histogram with the CPU Time used

  //Event Info
  TH1F* hNEvt;             ///<  Histogram with the number of events and all the events that passed each cut
  TH1F* hEvtMult;          ///<  Histogram with the event Multiplicity Before any physics selection
  TH1F* hEvtMultAftEvSel;  ///<  Histogram with the event Multiplicity After the physics selection
  TH1F* hEvtVtxXYBefSel;   ///<  Histogram with the linear distance in the XY plane of the primary vertex before any selection on vertex position
  TH1F* hEvtVtxZBefSel;    ///<  Histogram with the linear distance along the Z axis of the primary vertex before any selection on vertex position
  TH1F* hEvtVtxXY;         ///<  Histogram with the linear distance in the XY plane of the primary vertex
  TH1F* hEvtVtxZ;          ///<  Histogram with the linear distance along the Z axis of the primary vertex
  TH1F* hEvtVtxZMCGen;     ///<  Histogram with the linear distance along the Z axis of the MC vertex of generated MC events
  TH1F* hEvtVtxZMCPhysSel; ///<  Histogram with the linear distance along the Z axis of the MC vertex after physics selection
  TH1F* hEvtVtxZMCReco;    ///<  Histogram with the linear distance along the Z axis of the MC vertex of the events with reconstructed vertex

  //Track Info
  TH1F* hNTrk; ///<  Histogram with the number of tracks and all the tracks that passed each cut
  //->Track cuts information divided into before [0] and after [1] the cut
  TH1F* hTrkTPCCls[2];                     ///<  Histogram with the number of TPC crossed Rows of all not accepted and accepted tracks
  TH1F* hTrkTPCRows[2];                    ///<  Histogram with the number of TPC crossed Rows of all not accepted and accepted tracks
  TH1F* hTrkTPCRatioRowsFindCls[2];        ///<  Histogram with the ratio between the number of TPC crossed Rows and the number of findable clusters in TPC of all not accepted and accepted tracks
  TH1F* hTrkTPCChi2NDF[2];                 ///<  Histogram with the reduced chi2 in TPC of all not accepted and accepted tracks
  TH1F* hTrkITSChi2NDF[2];                 ///<  Histogram with the reduced chi2 in ITS of all not accepted and accepted tracks
  TH1F* hTrkActiveLength[2];               ///<  Histogram with the length of the track in the active region
  TH1F* hTrkITSTPCMatch[2];                ///<  Histogram with the match between ITS and TPC
  TH1F* hTrkDCAxy[2];                      ///<  Histogram with the DCA XY of all not accepted and accepted tracks
  TH1F* hTrkDCAz[2];                       ///<  Histogram with the DCA Z of all not accepted and accepted tracks
#ifdef CHECKTRACKCUTS                      //Only if required
  TH2I* hTrkTPCClsCorr[2][3];              ///<  Correlation between pt [0] - eta [1] - phi [2] and the number of TPC clusters before and after the cut
  TH2I* hTrkTPCRowsCorr[2][3];             ///<  Correlation between pt [0] - eta [1] - phi [2] and the number of TPC crossed rows before and after the cut
  TH2I* hTrkTPCRatioRowsFindClsCorr[2][3]; ///<  Correlation between pt [0] - eta [1] - phi [2] and the ratio between the number of TPC crossed rows and findable clusters in TPC before and after the cut
  TH2I* hTrkTPCChi2NDFCorr[2][3];          ///<  Correlation between pt [0] - eta [1] - phi [2] and the TPC reduced chi2 before and after the cut
  TH2I* hTrkITSChi2NDFCorr[2][3];          ///<  Correlation between pt [0] - eta [1] - phi [2] and the ITS reduced chi2 before and after the cut
  TH2I* hTrkActiveLengthCorr[2][3];        ///<  Correlation between pt [0] - eta [1] - phi [2] and the track length in TPC active region before and after the cut
  TH2I* hTrkITSTPCMatchCorr[2][3];         ///<  Correlation between pt [0] - eta [1] - phi [2] and the track matched between ITS and TPC active region before and after the cut
  TH2I* hTrkDCAxyCorr[2][3];               ///<  Correlation between pt [0] - eta [1] - phi [2] and the DCAxy before and after the cut
  TH2I* hTrkDCAzCorr[2][3];                ///<  Correlation between pt [0] - eta [1] - phi [2] and the DCAz before and after the cut
#endif
  TH1F* hCutVariation; ///<  Histogram with the number of tracks which pass each cut
  //->TOF information
  TH1F* hTOFResidualX;                        ///<  Histogram with the Impact Residual X
  TH1F* hTOFResidualZ;                        ///<  Histogram with the Impact Residual Z
  TH2F* hPadDist;                             ///<  Histogram with the Impact Residual X and Residual Z values
  TH2F* hTOFDist;                             ///<  Histogram with the distributions of the TOF strips and sectors matched to tracks
  TH1F* hTOFChannel;                          ///<  Histogram with the Channel in the TOF
  TH1F* hT0;                                  ///<  Histogram with the T0 used for each track
  TH1F* hT0Resolution;                        ///<  Histogram with the resolution on the T0
  TH1F* hTimeOfFlightRes;                     ///<  Histogram to compute the Time Of Flight resolution
  TH1F* hTimeOfFlightTOFRes;                  ///<  Histogram to compute the Time Of Flight resolution for events without the T0 Fill
  TH1F* hTimeOfFlightGoodRes;                 ///<  Histogram to compute the Time Of Flight resolution for tracks with good matching
  TH1F* hTimeOfFlightResNoMismatch;           ///<  Histogram to compute the Time Of Flight resolution for PID consistent with TPC for Pi K P
  TH2F* hTimeOfFlightResFine;                 ///<  Histogram to compute the Time Of Flight resolution as a function of the matched tracks to TOF
  TH1F* hTimeOfFlightResFinePerEvent;         ///<  Histogram to compute the Time Of Flight resolution per event, this particular one should not be added to the output list as it yields no information but it is rather auxiliary to the computation of the TOF resolution as a function of the TOF tracks
  TH2I* hBeta;                                ///<  Histogram with the track beta vs the track momentum
  TProfile* hBetaExpected[kExpSpecies];       ///<  TProfile with the track beta vs the track momentum obtained with the expected time
  TProfile* hBetaExpectedTOFPID[kExpSpecies]; ///<  TProfile with the track beta vs the track momentum obtained with the expected time but with the 3sigma TOF PID on the particle hypothesis
  TProfile* hTOFSepVsP[kExpSpecies - 1];      ///<  TProfile with the difference of the track nsigma vs the track momentum
  TProfile* hTOFSepVsPt[kExpSpecies - 1];     ///<  TProfile with the difference of the track nsigma vs the track transverse momentum
  TH2I* hBetaNoMismatch;                      ///<  Histogram with the track beta vs the track momentum with a cut on the maximum number of clusters to reduce the mismatch
  TH2I* hBetaNoMismatchEtaCut;                ///<  Histogram with the track beta vs the track momentum with a cut on the maximum number of clusters to reduce the mismatch and a cut on the eta range
  TH2I* hBetaNoMismatchEtaCutOut;             ///<  Histogram with the track beta vs the track momentum with a cut on the maximum number of clusters to reduce the mismatch and a lower cut on the eta range
  TH2I* hBetaCentral;                         ///<  Histogram with the track beta vs the track momentum for central events
  TH2I* hBetaNoMismatchCentral;               ///<  Histogram with the track beta vs the track momentum for central events with a cut on the maximum number of clusters to reduce the mismatch
  TH2I* hBetaNoMismatchCentralEtaCut;         ///<  Histogram with the track beta vs the track momentum for central events with a cut on the maximum number of clusters to reduce the mismatch and a cut on the eta range
  TH2I* hBetaNoMismatchCentralEtaCutOut;      ///<  Histogram with the track beta vs the track momentum for central events with a cut on the maximum number of clusters to reduce the mismatch and a cut on the eta range
  TH3F* hBetaMC[3];                           ///<  Histogram with the track beta vs the track momentum with MC information on the PID and on the production mode
  TH3F* hBetaMCMother[kExpSpecies + 1];       ///<  Histogram with the track beta vs the track momentum with MC information on the PID and on the decay mother
  TH2F* hBetaMCMotherMode[kExpSpecies + 1];   ///<  Histogram with the track MC information on the PID and on the decay mother/production mode
  TH1F* hChannelEvents;                       ///<  Histogram with the numbe of events used in the TOF channel-time mode
  TH2I* hChannelTime;                         ///<  Histogram with the measured time at TOF divided into each channel (or strip) -> Used to get the mismatch
  TH1F* hChannelHits;                         ///<  Histogram with the numer of TOF Hits (matched and not matched) per channel
  TH1F* hChannelClusters;                     ///<  Histogram with the numer of TOF Clusters (matched and not matched) per channel
  TH1F* hTOFClusters;                         ///<  Histogram with the number of TOF clusters per track
  TH1F* hTOFClustersDCApass;                  ///<  Histogram with the number of TOF clusters per track, for tracks which passed the DCA cut for primaries
  //->TPC information
  TH2I* hTPCdEdx;                               ///<  Histogram with the track energy loss in the TOC vs the track momentum
  TH2I* hTPCdEdxTPCp;                           ///<  Histogram with the track energy loss in the TOC vs the track momentum for TPC
  TProfile* hdEdxExpected[kExpSpecies + 2];     ///<  TProfile with the track dEdx vs the track momentum
  TProfile* hdEdxExpectedTPCp[kExpSpecies + 2]; ///<  TProfile with the track dEdx vs the track momentum for TPC
  //->TPC + TOF information
  TH2F* hTPCTOFSeparation[kExpSpecies][kPtBins]; ///<  Histogram with the PID separation of the TPC and TOF signal
#ifdef CHECKCOMPUTEDVALUES                       // Only if checks on computed values are required
  TH1F* hTOFExpectedComputed[kExpSpecies];       ///<  Histogram with the difference of the expected TOF calculated with the formula and extracated from the track
  TH1F* hTOFExpectedComputedTPC[kExpSpecies];    ///<  Histogram with the difference of the expected TOF calculated with the formula with the TPC momentum and extracated from the track
  TH1F* hTOFMomComputed[kExpSpecies];            ///<  Histogram with the difference of the track momentum and the one which is calculated from the expected time
  TH1F* hTOFMomComputedTPC[kExpSpecies];         ///<  Histogram with the difference of the track momentum of the TPC and the one which is calculated from the expected time
#endif
#ifdef BUILDTOFDISTRIBUTIONS                              // Build TOF distributions only if requested
  TH1F* hTOF[kPtBins][kCharges][kSpecies];                ///<  Distribution for T-Texp-T0
  TH1F* hTOFNoYCut[kPtBins][kCharges][kSpecies];          ///<  Distribution for T-Texp-T0 without the cut on rapidity
  TH1F* hTOFSigma[kPtBins][kCharges][kSpecies];           ///<  Sigma distributions for T-Texp-T0
  TH1F* hTOFNoMismatch[kPtBins][kCharges][kSpecies];      ///<  Distribution for T-Texp-T0 without Mismatch, removed with the information on the TPC
  TH1F* hTOFSigmaNoMismatch[kPtBins][kCharges][kSpecies]; ///<  Distribution for T-Texp-T0 without Mismatch, removed with the information on the TPC
#endif
#ifdef BUILDT0PLOTS
  TH2F* hT0VsVtxZ[4];     ///<  Distribution for T0 vs Vtx Z position for all the ways of computing T0 (i.e. A&C, A, C and TOF)
  TH2F* hT0VsVtxZbest[4]; ///<  Distribution for T0 vs Vtx Z position for all the ways of computing T0 (i.e. A&C, A, C best)
#endif

  //
  //MC Info
  //
  TH1F* hNumMatchMC[2][3];                   ///<  Matching efficiency numerator with MC information on PID
  TH1F* hDenMatchMC[2][3];                   ///<  Matching efficiency denominator with MC information on PID
  TH1F* hNumMatchPrimMC[2][3];               ///<  Matching efficiency numerator with MC information on PID and on Primary production
  TH1F* hDenMatchPrimMC[2][3];               ///<  Matching efficiency denominator with MC information on PID and on Primary production
  TH2F* hNumMatchPrimMCYCut[2][3];           ///<  Matching efficiency numerator with MC information on PID and on Primary production with a cut on |y| < 0.5
  TH2F* hDenMatchPrimMCYCut[2][3];           ///<  Matching efficiency denominator with MC information on PID and on Primary production with a cut on |y| < 0.5
  TH2F* hNumMatchMultTrkTRDOut[2][3];        ///<  Matching efficiency numerator with kTRDOut flag and MC information on PID
  TH2F* hDenMatchMultTrkTRDOut[2][3];        ///<  Matching efficiency denominator with kTRDOut flag and MC information on PID
  TH2F* hNumMatchMultTrkNoTRDOut[2][3];      ///<  Matching efficiency numerator without kTRDOut flag and with MC information on PID
  TH2F* hDenMatchMultTrkNoTRDOut[2][3];      ///<  Matching efficiency denominator without kTRDOut flag and with MC information on PID
  TH2F* hNumMatchMultTrk[2][3];              ///<  Matching efficiency numerator with kTIME, kTRDOut flags and MC information on PID
  TH2F* hDenMatchMultTrk[2][3];              ///<  Matching efficiency denominator with kTIME, kTRDOut flags and MC information on PID
  TH2F* hNumMatchMultTrkInc[2];              ///<  Matching efficiency numerator with kTIME, kTRDOut flags
  TH2F* hDenMatchMultTrkInc[2];              ///<  Matching efficiency denominator with kTIME, kTRDOut flags
  TH2F* hDenTrkTrigger[2][3];                ///<  Generated particles with MC truth on PID for all events: that passed and didn't pass Physics Selection
  TH1F* hDenTrkMCVertexZ[2][3];              ///<  Generated particles with MC truth on PID for all events: that passed and didn't pass Physics Selection. But that are inside the sample selection in the MC
  TH1F* hDenTrkVertex[2][3];                 ///<  Generated particles with MC truth on PID for events that passed Vertex Cuts
  TH1F* hDenTrkVertexMCVertexZ[2][3];        ///<  Generated particles with MC truth on PID for events that passed Vertex Cuts
  TH2F* hDenPrimMCYCut[2][3];                ///<  Pt Distribution of Primary Particles with MC Truth on PID, that passed Physics Selection and Event Selection with a cut on the max y
  TH2F* hDenPrimMCEtaCut[2][3];              ///<  Pt Distribution of Primary Particles with MC Truth on PID, that passed Physics Selection and Event Selection with a cut on the max eta
  TH2F* hDenPrimMCEtaYCut[2][3];             ///<  Pt Distribution of Primary Particles with MC Truth on PID, that passed Physics Selection and Event Selection with a cut on the max eta and max y
  TH2F* hNumPrimMCTrueMatch[2][3];           ///<  Pt Distribution of Tracks from primary particles with MC Truth on PID, with true match in the TOF detector
  TH2F* hNumPrimMCTrueMatchYCut[2][3];       ///<  Pt Distribution of Tracks from primary particles with MC Truth on PID, with true match in the TOF detector with a Y cut
  TH2F* hNumPrimMCTrueMatchYCutTPC[2][3];    ///<  Pt Distribution of Tracks from primary particles with MC Truth on PID, with true match in the TOF detector with a Y cut and a TPC 5sigma cut on the signal for pi/k/p
  TH2F* hNumPrimMCConsistentMatchYCut[2][3]; ///<  Pt Distribution of Tracks from primary particles with MC Truth on PID, without true match in the TOF detector but compatible PID with a Y cut
  TH2F* hNumPrimMCMotherMatchYCut[2][3];     ///<  Pt Distribution of Tracks from primary particles with MC Truth on PID, with mother with true match in the TOF detector with a Y cut
  TH3S* hDenMatchPrimNoCut[2][3];            ///<  Matching efficiency denominator with MC information on PID and on Primary production with pT eta and phi
  TH3S* hDenPrimMCNoCut[2][3];               ///<  Distribution of Primary Particles with MC Truth on PID, that passed Physics Selection and Event Selection with pT eta and phi

  //histograms for matching efficiency calculation (Also in Real Data)
  //Positive / Negative
  TH1F* hNumMatch[2];            ///<  Matching efficiency numerator Pt Distribution
  TH1F* hDenMatch[2];            ///<  Matching efficiency denominator Pt Distribution
  TH1F* hNumMatchEta[2];         ///<  Matching efficiency numerator Eta Distribution
  TH1F* hDenMatchEta[2];         ///<  Matching efficiency denominator Eta Distribution
  TH1F* hNumMatchphiOut[2];      ///<  Matching efficiency numerator Phi Distribution
  TH1F* hDenMatchphiOut[2];      ///<  Matching efficiency denominator Phi Distribution
  TH1F* hNumMatchEtaPtMa[2];     ///<  Matching efficiency numerator Eta Distribution with a Cut on low pt (< 0.5 GeV)
  TH1F* hDenMatchEtaPtMa[2];     ///<  Matching efficiency denominator Eta Distribution with a Cut on low pt (< 0.5 GeV)
  TH1F* hNumMatchphiOutPtMa[2];  ///<  Matching efficiency numerator Phi Distribution with a Cut on low pt (< 0.5 GeV)
  TH1F* hDenMatchphiOutPtMa[2];  ///<  Matching efficiency denominator Phi Distribution with a Cut on low pt (< 0.5 GeV)
  TH1F* hNumMatchTRDOut[2];      ///<  Matching efficiency numerator Pt Distribution with TRD out
  TH1F* hDenMatchTRDOut[2];      ///<  Matching efficiency denominator Pt Distribution with TRD out
  TH1F* hNumMatchNoTRDOut[2];    ///<  Matching efficiency numerator Pt Distribution without TRD out
  TH1F* hDenMatchNoTRDOut[2];    ///<  Matching efficiency denominator Pt Distribution without TRD out
  TH3I* hNumMatchPtEtaPhiout[2]; ///<  Matching efficiency numerator Pt, Eta and Phi out Distribution
  TH3I* hDenMatchPtEtaPhiout[2]; ///<  Matching efficiency denominator Pt, Eta and Phi out Distribution

  //Positive / Negative
  //Pion / Kaon / Proton
  TH1F* hNumMatchTPC[2][3]; ///<  Matching efficiency numerator Pt Distribution with information of the TPC PID
  TH1F* hDenMatchTPC[2][3]; ///<  Matching efficiency denominator Pt Distribution with information of the TPC PID

  //DCA Histograms
  //-> Data and MC
  TH2F* hDCAxy[2][3][kPtBins];           ///< DCAxy Distribution in Pt and Multiplicity/Centrality bins for all the reconstructed tracks identified via a 2sigma TOF/TPC cut
  TH2F* hDCAxyGoldenChi2[2][3][kPtBins]; ///< DCAxy Distribution in Pt and Multiplicity/Centrality bins for all the reconstructed tracks but with the cut on the golden Chi2
  //-> MC only
  TH1F* hDCAxyPrimMC[2][3][kPtBins];   ///< DCAxy Distribution in Pt bins for Primary reconstructed tracks identified with MC Truth
  TH1F* hDCAxySecStMC[2][3][kPtBins];  ///< DCAxy Distribution in Pt bins for Secondary from Strangeness reconstructed tracks identified with MC Truth
  TH1F* hDCAxySecMatMC[2][3][kPtBins]; ///< DCAxy Distribution in Pt bins for Secondary from Material reconstructed tracks identified with MC Truth

  //
  AliAnalysisTaskTOFSpectra(const AliAnalysisTaskTOFSpectra&);            //! Copy constructor
  AliAnalysisTaskTOFSpectra& operator=(const AliAnalysisTaskTOFSpectra&); //! Not implemented

  ClassDef(AliAnalysisTaskTOFSpectra, 15); //AliAnalysisTaskTOFSpectra used for the Pi/K/p analysis with TOF
};

#endif
