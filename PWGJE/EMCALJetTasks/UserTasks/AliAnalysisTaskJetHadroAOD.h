#ifndef ALIANALYSISTASKJETHADROAOD_H
#define ALIANALYSISTASKJETHADROAOD_H

////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                        //
//                        Analysis for Jet Hadrochemistry                                 //
//                                                                                        //
//    This analysis extracts pT-spectra of charged kaons, protons, and pions              //
//                      for the inclusive event and in jets.                              //
//   It is based on particles identification via the dE/dx signal of the TPC              //
//                   and the time of flight nsigma from the TOF.                          //
//                              This is the AOD version.                                  //
//                                                                                        //
// Author: Sierra Cantway (Weyhmiller) <sierra.lisa.weyhmiller@cern.ch>, Yale University  //
//      Author: Mesut Arslandok <mesut.arslandok@cern.ch>, Yale University                //
//                                                                                        //
////////////////////////////////////////////////////////////////////////////////////////////

class TH1;
class THn;
class TH1F;
class TH2F;
class TH3F;
class TH3F;
class TList;
class TTree;
class TObjArray;
class AliAODEvent;
class AliAODTrack;
class AliPIDResponse;
class AliPIDCombined;
class AliJetContainer;
class AliAnalysisTaskRho;


#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliPIDCombined.h"
#include "AliTPCdEdxInfo.h"
#include "THnSparse.h"
#include "TClonesArray.h"
#include "THn.h"
#include "TVectorF.h"
#include "TTreeStream.h"
#include "AliEventCuts.h"
#include "TRandom3.h"
#include "AliFJWrapper.h"

class AliAnalysisTaskJetHadroAOD : public AliAnalysisTaskEmcalJet {
public:

  AliEventCuts fEventCuts;     /// Event cuts

  // ---------------------------------------------------------------------------------
  //                           Constructor and Destructor
  // ---------------------------------------------------------------------------------

  AliAnalysisTaskJetHadroAOD(const char *name);
  AliAnalysisTaskJetHadroAOD();
  virtual ~AliAnalysisTaskJetHadroAOD();

  enum kPDGpart{
    kPDGel=11,
    kPDGpi=211,
    kPDGka=321,
    kPDGpr=2212,
    kPDGde=1000010020,
    kPDGmu=13,
    kPDGla=3122,
  };
  //
  // ---------------------------------------------------------------------------------
  //                                    Methods
  // ---------------------------------------------------------------------------------
  //
  virtual void   UserCreateOutputObjects();            // create output objects
  virtual Bool_t Run();                                // Equivalent to UserExec --> tuned in
  virtual void   Terminate(Option_t *);                // run only once and terminate

  // ---------------------------------------------------------------------------------
  //                                    Settings
  // ---------------------------------------------------------------------------------

  void   Initialize();
  void   SetFilterBits(UInt_t bits)                                   {fAOD_FilterBits      = bits;}
  void   SetIsMCtrue(Bool_t isMCdata = kTRUE)                         {fMCtrue              = isMCdata;};
  void   SetYear(const Int_t ifYear = 0)                              {fYear                = ifYear;}
  void   SetPeriodName(const TString ifPeriodName = "")               {fPeriodName          = ifPeriodName;}
  void   SetPassIndex(const Int_t ifPassIndex = 0)                    {fPassIndex           = ifPassIndex;}
  // Some boolian settings
  void   SetSmallOut(const Bool_t ifSmallOut = kTRUE)               {fSmallOut           = ifSmallOut;}
  void   SetIncludeITScuts(const Bool_t ifITSCuts = kTRUE)            {fIncludeITS          = ifITSCuts;}
  void   SetFilljetsFJBGTree(const Bool_t ifjetsFJBGTree = kTRUE)     {fFilljetsFJBGTree    = ifjetsFJBGTree;}
  void   SetFillJetsFJBGConst(const Bool_t ifJetsFJBGConst = kTRUE)     {fFillJetsFJBGConst     = ifJetsFJBGConst;}
  void   SetDoIncTracks(const Bool_t ifdoIncTracks = kTRUE)               {fDoIncTracks        = ifdoIncTracks;}
  void   SetDoPerpCone(const Bool_t ifdoPerpCone = kTRUE)               {fDoPerpCone        = ifdoPerpCone;}
  void   SetDoMatchPerpCone(const Bool_t ifdoMatchPerpCone = kTRUE)               {fDoMatchPerpCone        = ifdoMatchPerpCone;}
  void   SetDoRandCone(const Bool_t ifdoRandCone = kTRUE)               {fDoRandCone        = ifdoRandCone;}
  void   SetDoFastJet(const Bool_t ifFastJet = kTRUE)               {fDoFastJet         = ifFastJet;}
  void   SetDoEMCJet(const Bool_t ifEMCJet = kTRUE)               {fDoEMCJet         = ifEMCJet;}
  void   SetDoEventTree(const Bool_t ifEventTree = kTRUE)               {fDoEventTree         = ifEventTree;}
  void   SetDoEventHistos(const Bool_t ifEventHisto = kTRUE)               {fDoEventHistos         = ifEventHisto;}
  void   SetFillFastJet(const Bool_t ifFastJet = kTRUE)               {fFillFastJet         = ifFastJet;}
  void   SetFillEMCJet(const Bool_t ifEMCJet = kTRUE)               {fFillEMCJet         = ifEMCJet;}
  void   SetFillJetsEMCConst(const Bool_t ifJetsEMCConst = kTRUE)     {fFillJetsEMCConst     = ifJetsEMCConst;}
  void   SetFillJetsFJConst(const Bool_t ifJetsFJConst = kTRUE)     {fFillJetsFJConst     = ifJetsFJConst;}
  void   SetFillJetsEMCBG(const Bool_t ifJetsEMCBG = kTRUE)     {fFillJetsEMCBG     = ifJetsEMCBG;}
  void   SetFillJetsEMCBGConst(const Bool_t ifJetsEMCBGConst = kTRUE)     {fFillJetsEMCBGConst     = ifJetsEMCBGConst;}
  void   SetJetMinPtSub(const Double_t jetminptsub = -1000.0)         {fjetMinPtSub         = jetminptsub;}
  void   SetJetMaxPtSub(const Double_t jetmaxptsub = -1000.0)         {fjetMaxPtSub         = jetmaxptsub;}
  void   SetJetMinArea(const Double_t jetminarea = -1000.0)         {fjetMinArea         = jetminarea;}
  void   SetMinCent(const Double_t mincent = 0.0)         {fcent_min         = mincent;}
  void   SetMaxCent(const Double_t maxcent = 100.0)         {fcent_max         = maxcent;}
  void   SetDoRapCut(const Bool_t ifRapCut = kFALSE)         {fDoRapCut         = ifRapCut;}
  void   SetEtaCut(const Double_t etacut = 100.0)         {fEtaCut         = etacut;}
  void   SetYCut(const Double_t ycut = 100.0)         {fYCut         = ycut;}

  void   SetFillTreeMC(const Bool_t ifTreeMC = kFALSE)                {fFillTreeMC= ifTreeMC;}
  void   SetFillIncTracks(const Bool_t ifIncTracks = kTRUE)           {fFillIncTracks       = ifIncTracks;}
  void   SetFill_TPC(const Bool_t if_TPC = kTRUE)           {fFill_TPC       = if_TPC;}
  void   SetFillpTPC_pT(const Bool_t ifpTPC_pT = kTRUE)           {fFillpTPC_pT       = ifpTPC_pT;}
  void   SetFillp_pT(const Bool_t ifp_pT = kTRUE)           {fFillp_pT       = ifp_pT;}
  void   SetFillpTPC_p(const Bool_t ifpTPC_p = kTRUE)           {fFillpTPC_p       = ifpTPC_p;}
  void   SetFill_TOF(const Bool_t if_TOF = kTRUE)           {fFill_TOF       = if_TOF;}
  void   SetFill_TOF_expecs(const Bool_t if_TOF_expecs = kTRUE)           {fFill_TOF_expecs       = if_TOF_expecs;}
  void   SetFill_TPC_expecs(const Bool_t if_TPC_expecs = kTRUE)           {fFill_TPC_expecs       = if_TPC_expecs;}
  void   SetUseCouts(const Bool_t ifUseCouts = kFALSE)                {fUseCouts            = ifUseCouts;}
  void   SetPercentageOfEvents(const Int_t nPercentageOfEvents = 0)   {fPercentageOfEvents = nPercentageOfEvents;}
  void   SetDoEmbedding(const Bool_t ifDoEmbedding = kFALSE)   {fDoEmbedding   =  ifDoEmbedding;}
  void   SetMultTrueFlag(const Bool_t ifMultTrueFlag = kFALSE)   {fMultTrueFlag   =  ifMultTrueFlag;}

  // Setters for the eta, momentum, dEdx, etc
  void   SetDeDxBins(const Int_t ndEdxBins, Float_t dEdxBins[]) {
      fNdEdxBins = ndEdxBins;
      fdEdxBins.resize(fNdEdxBins+1);
      for (Int_t i=0; i<(fNdEdxBins+1); i++) fdEdxBins[i] = dEdxBins[i];
  }

  void   SetBetaBins(const Int_t nbetaBins, Float_t betaBins[]) {
      fNBetaBins = nbetaBins;
      fBetaBins.resize(fNBetaBins+1);
      for (Int_t i=0; i<(fNBetaBins+1); i++) fBetaBins[i] = betaBins[i];
  }

  void   SetTOFNSigmaBins(const Int_t nTOFNSigmaBins, Float_t TOFNSigmaBins[]) {
      fNTOFNSigmaBins = nTOFNSigmaBins;
      fTOFNSigmaBins.resize(fNTOFNSigmaBins+1);
      for (Int_t i=0; i<(fNTOFNSigmaBins+1); i++) fTOFNSigmaBins[i] = TOFNSigmaBins[i];
  }

  void   SetTPCMmom_Bins(const Int_t nTPCMombins, Float_t TPCMombins[]) {
      fNTPCMom_Bins = nTPCMombins;
      fTPCMom_Bins.resize(fNTPCMom_Bins+1);
      for (Int_t i=0; i<(fNTPCMom_Bins+1); i++) fTPCMom_Bins[i] = TPCMombins[i];
  }

  void   SetTOFMom_Bins(const Int_t nTOFMombins, Float_t TOFMombins[]) {
      fTOFMom_NBins = nTOFMombins;
      fTOFMom_Bins.resize(fTOFMom_NBins+1);
      for (Int_t i=0; i<(fTOFMom_NBins+1); i++) fTOFMom_Bins[i] = TOFMombins[i];
  }

  void   SetEta_Bins(const Int_t nEtabins, Float_t Etabins[]) {
      fNEta_Bins = nEtabins;
      fEta_Bins.resize(fNEta_Bins+1);
      for (Int_t i=0; i<(fNEta_Bins+1); i++) fEta_Bins[i] = Etabins[i];
  }

  void   SetTPCmom_choice(const Float_t TPCmom = 0.)         {fSetTPCmom            = TPCmom;}
  void   SetTOFmom_choice(const Float_t TOFmom = 0.)         {fSetTOFmom            = TOFmom;}
  void   SetBetamom_choice(const Float_t Betamom = 0.)         {fSetBetamom            = Betamom;}
  void   SetEta_choice(const Float_t eta = 0.)         {fSetEta            = eta;}

//From jet extractor task
  void   SetJetMatchingRadius(Double_t val)                  { fJetMatchingRadius = val; }
  void   SetTruthLabelRange(Int_t min, Int_t max)            { fTruthMinLabel = min; fTruthMaxLabel = max; }
  void   SetSaveMCInformation(Bool_t val) {fSaveMCInformation = val; fInitialized = kFALSE;}
  void   SetJetMatchingSharedPtFraction(Double_t val)        { fJetMatchingSharedPtFraction = val; }

  //function to add jet task
  void   AddJet(AliJetContainer* jet = 0)                         {fJetContainer        = jet; }

  Double_t  GetDistance(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2)
  {
    Double_t deltaPhi = TMath::Min(TMath::Abs(phi1-phi2),TMath::TwoPi() - TMath::Abs(phi1-phi2));
    return TMath::Sqrt((eta1-eta2)*(eta1-eta2) + deltaPhi*deltaPhi);
  }

private:

  AliAnalysisTaskJetHadroAOD(const AliAnalysisTaskJetHadroAOD&);
  AliAnalysisTaskJetHadroAOD& operator=(const AliAnalysisTaskJetHadroAOD&);

  // ---------------------------------------------------------------------------------
  //                                   Functions
  // ---------------------------------------------------------------------------------

  void FindJetsEMC();                          // Find and Fill Jets with EMCAL framework
  void FindJetsFJ();                          // Find and Fill Jets with FJ framework
  void FillIncTracksReal();                   // Fill all inclusive track information
  void FillTreeMC();
  void FillMCJets();
  void FillEmbJets();
  void GetExpecteds(AliAODTrack* track);
  void FillEventTree();
  void FillEventHistos();
  void GetTrueJetPtFraction(AliEmcalJet* jet, Double_t& truePtFraction, Double_t& truePtFraction_mcparticles);
  bool PerformGeometricalJetMatching(AliJetContainer& contBase, AliJetContainer& contTag, double maxDist);
  void DoJetMatching();
  void GetMatchedJetObservables(AliEmcalJet* jet,  Double_t& detJetPt, Double_t& partJetPt, Double_t& detJetPhi, Double_t& detJetEta, Double_t& partJetPhi, Double_t& partJetEta, Double_t& detJetDistance, Double_t& partJetDistance);

  //
  Bool_t CountEmptyEvents();  // Just count if there is empty events
  void SetCutBitsAndSomeTrackVariables(AliAODTrack* track);
  // ---------------------------------------------------------------------------------
  //                                   Members
  // ---------------------------------------------------------------------------------

  AliPIDResponse   * fPIDResponse;            //!<! PID response object
  AliAODEvent      * fAOD;                    //!<! AOD object
  AliMCEvent       * fMCEvent;                //!<! corresponding MC event
  TList            * fListHist;               //!<! list for histograms
  UInt_t            fAOD_FilterBits;
  AliPIDCombined   * fPIDCombined;            //!<! combined PID object
  AliStack         * fMCStack;                //!<! stack object to get Mc info
  const AliAODVertex * fVertex_AOD;               //!<! primary vertex

  TTreeSRedirector * fTreeSRedirector;        //!<! temp tree to dump output
  TTree            * fTreeMC;                 //!<! tree for mc samples
  TTree            * fTreeCuts;               //!<! tree to save all variables for control plots
  TTree            * fTreejetsEMCconst;       //!<! tree for EMCal signal jet constituents
  TTree            * fTreejetsEMCBGconst;       //!<! tree for EMCal background jet constituents
  TTree            * fTreejetsFJ;             //!<! tree for fastjet signal jets
  TTree            * fTreejetsFJBG;           //!<! tree for fastjet background jets
  TTree            * fTreejetsFJconst;          //!<! tree for fastjet signal jet constituents
  TTree            * fTreejetsFJBGconst;          //!<! tree for fastjet signal jet constituents
  TTree            * fTreejetEvents;          //!<! tree for event level observables
  TTree            * fTreejetsEMC;            //!<! tree for EMCal signal jets
  TTree            * fTreejetsEMCBG;            //!<! tree for EMCal background jets
  TRandom3         * fRandom;                 //!<! random number generator


  TString            fPeriodName;
  Int_t              fYear;
  Int_t              fPassIndex;

  Int_t             fPercentageOfEvents;     // when only a fPercentageOfEvents is enough

  Bool_t            fSmallOut;              // flag for small output
  Bool_t            fMCtrue;                 // flag if real data or MC is processed
  Bool_t            fIncludeITS;             // decide whether to use ITS or not
  Bool_t            fFillIncTracks;          // switch whether to fill tracks tree
  Bool_t            fFill_TPC;               // switch whether to fill TPC histos
  Bool_t            fFillpTPC_pT;           // switch whether to fill pTPC v pT conversion histo
  Bool_t            fFillp_pT;           // switch whether to fill p v pT conversion histo
  Bool_t            fFillpTPC_p;           // switch whether to fill pTPC v p conversion histo
  Bool_t            fFill_TOF;           // switch whether to fill TOF histos
  Bool_t            fFill_TOF_expecs;           // switch whether to fill expected TOF histos
  Bool_t            fFill_TPC_expecs;           // switch whether to fill expected TPC histos
  Bool_t            fFilljetsFJBGTree;         // switch whether to fill BG Jets FJ tree
  Bool_t            fFillJetsFJBGConst;        // switch whether to fill jetsEMCBG constituent tree
  Bool_t            fFillFastJet;         // switch whether to fill FJ tree
  Bool_t            fFillEMCJet;         // switch whether to fill EMC tree
  Bool_t            fDoIncTracks;         // switch whether to use inclusive tracks
  Bool_t            fDoPerpCone;         // switch whether to use perpendicular cone tracks
  Bool_t            fDoMatchPerpCone;         // switch whether to use perpendicular cone tracks for matched jets only
  Bool_t            fDoRandCone;         // switch whether to use random cone tracks
  Bool_t            fDoFastJet;         // switch whether to use FJ jets
  Bool_t            fDoEMCJet;         // switch whether to use EMC jets
  Bool_t            fDoEventTree;         // switch whether to use tree for event variables
  Bool_t            fDoEventHistos;         // switch whether to use histos for event variables
  Bool_t            fFillJetsEMCConst;        // switch whether to fill jetsEMC constituent tree
  Bool_t            fFillJetsFJConst;        // switch whether to fill jetsFJ constituent tree
  Bool_t            fFillJetsEMCBG;        // switch whether to fill jetsEMCBG tree
  Bool_t            fFillJetsEMCBGConst;        // switch whether to fill jetsEMCBG constituent tree
  Double_t          fjetMinPtSub;            // minimium jet pt after subtraction to keep jet
  Double_t          fjetMaxPtSub;            // maximum jet pt after subtraction to keep jet
  Double_t          fjetMinArea;            // minimium jet area to keep jet
  Float_t           fcent_min;            // minimium centrality cut
  Float_t           fcent_max;            // maximium centrality cut
  Bool_t            fDoRapCut;            //employ rap cut or not
  Float_t           fEtaCut;              //cut on absolute eta
  Float_t           fYCut;              //cut on absolute y

  Bool_t            fFillTreeMC;

  Bool_t            fUseCouts;               // for debugging

  Int_t             fSetTPCmom;             // set which momentum to use for TPC plots
  Int_t             fSetTOFmom;             // set which momentum to use for TOF plots
  Int_t             fSetBetamom;            // set which momentum to use for beta plots
  Int_t             fSetEta;                // set which eta or y to use for histos

  Int_t             fNdEdxBins;           //!<! number of bins for dEdx histograms
  Int_t             fNBetaBins;           //!<! number of bins for the beta histogram
  Int_t             fNTOFNSigmaBins;      //!<! number of bins for the TOF nsigma histograms
  Int_t             fNTPCMom_Bins;        //!<! number of TPC momentum bins
  Int_t             fTOFMom_NBins;        //!<! number of TOF momentum bins
  Int_t             fNEta_Bins;           //!<! number of absolute eta bins

  std::vector<float>  fdEdxBins;          //!<! variable bins for dEdx histograms
  std::vector<float>  fBetaBins;          //!<! variable bins for beta histograms
  std::vector<float>  fTOFNSigmaBins;     //!<! variable bins for TOFnsigma histograms
  std::vector<float>  fTPCMom_Bins;       //!<! variable bins for TPC momentum histograms
  std::vector<float>  fTOFMom_Bins;       //!<! variable bins for TOF momentum histograms
  std::vector<float>  fEta_Bins;          //!<! variable bins for eta histograms

  Float_t           fNSigmasElTOF;           // TOF N sigma for Electron
  Float_t           fNSigmasMuTOF;           // TOF N sigma for Muon
  Float_t           fNSigmasPiTOF;           // TOF N sigma for Pion
  Float_t           fNSigmasKaTOF;           // TOF N sigma for Kaon
  Float_t           fNSigmasPrTOF;           // TOF N sigma for Proton
  Float_t           fNSigmasDeTOF;           // TOF N sigma for Deuteron

  Float_t           fDEdxEl;                 // Expected Electron dEdx
  Float_t           fDEdxKa;                 // Expected Kaon dEdx
  Float_t           fDEdxPi;                 // Expected Pion dEdx
  Float_t           fDEdxPr;                 // Expected Proton dEdx
  Float_t           fDEdxDe;                 // Expected Deuteron dEdx

  Float_t           fSigmaEl;                // Expected Electron sigma
  Float_t           fSigmaKa;                // Expected Kaon sigma
  Float_t           fSigmaPi;                // Expected Pion sigma
  Float_t           fSigmaPr;                // Expected Proton sigma
  Float_t           fSigmaDe;                // Expected Deuteron sigma

  Float_t           fNSigmasElTPC;           // TOF N sigma for Electron
  Float_t           fNSigmasPiTPC;           // TOF N sigma for Pion
  Float_t           fNSigmasKaTPC;           // TOF N sigma for Kaon
  Float_t           fNSigmasPrTPC;           // TOF N sigma for Proton
  Float_t           fNSigmasDeTPC;           // TOF N sigma for Deuteron

  Float_t           fTPCSignalMC;
  Float_t           fPtotMC;
  Float_t           fPtMC;
  Float_t           fEtaMC;
  Int_t             fSignMC;

  Double_t          fMCImpactParameter;
  Int_t             fNHardScatters;            // Number of hard scatterings
  Int_t             fNProjectileParticipants;  // Number of projectiles participants
  Int_t             fNTargetParticipants;      // Number of target participants
  Int_t             fNNColl;                   // Number of N-N collisions
  Int_t             fNNwColl;                  // Number of N-Nwounded collisions
  Int_t             fNwNColl;                  // Number of Nwounded-N collisons
  Int_t             fNwNwColl;                 // Number of Nwounded-Nwounded collisions

  Float_t           fPtot;                   // TPC momentum
  Float_t           fPVertex;                // TPC momentum
  Float_t           fPt;                     // Transverse momentum
  Float_t           fY;                      // rapidity

  Float_t            fCentrality;             // centrality information
  Float_t            fCentImpBin;
  Double_t           fVz;                     // Vertex position
  Int_t              fEventCountInFile;       // event count per job

  Float_t            fTPCSignal;              // Measured dE/dx
  Float_t            fEta;                    // pseudo rapidity
  Float_t            fNContributors;          // Ntracks
  Float_t            fPhi;                    // azimuthal angle
  Int_t              fSign;                   // sign of the particle

  AliJetContainer*   fJetContainer;   //!<! signal jet container
  AliJetContainer*   fbgJetContainer;   //!<! background jet container
  AliJetContainer*   fEmbJetContainer;   //!<! signal embedded jet container
  AliJetContainer*   fRecoJetContainer;   //!<! signal reco jet container
  AliJetContainer*   fGenJetContainer;   //!<! signal gen jet container
  Double_t           fJetPt;
  Double_t           fJetEta;
  Double_t           fJetPhi;
  Float_t            fjetRhoVal;
  Float_t            fjetRecoRhoVal;
  Float_t            fjetEmbRhoVal;
  Float_t            fjetGenRhoVal;
  Float_t            frhoFJ;
  Bool_t             fisGoodIncEvent;
  Bool_t             fFilledUECone_Emb;
  Bool_t             fFilledUECone_Gen;
  Bool_t             fFilledUECone_Rec;
  Bool_t             fhasAcceptedFJjet;
  Bool_t             fhasAcceptedEMCjet;
  Bool_t             fhasRealFJjet;
  Bool_t             fhasRealEMCjet;
  Int_t              fNumRealFJJets;
  Int_t              fNumRealEMCJets;
  Bool_t             fhasAcceptedFJjet_Truth;
  Bool_t             fhasRealFJjet_Truth;
  Int_t              fNumRealFJJets_Truth;

  Bool_t              fMultTrueFlag;

//From jet extractor task
  Bool_t              fDoPartLevelMatching;
  Bool_t              fDoDetLevelMatching;
  Bool_t              fIsEmbeddedEvent;
  TString             fMCParticleArrayName;     ///< Array name of MC particles in event (mcparticles)
  TClonesArray*       fMCParticleArray;         //!<! Array of MC particles in event (usually mcparticles)
  Double_t            fJetMatchingRadius;       ///< Matching radius for geometrically matching jets
  Int_t               fTruthMinLabel;           ///< min track label to consider it as true particle
  Int_t               fTruthMaxLabel;           ///< max track label to consider it as true particle
  Bool_t              fSaveMCInformation;       ///< save MC information
  Double_t            fJetMatchingSharedPtFraction; ///< Shared pT fraction required in matching
  Bool_t              fDoEmbedding;             //Flag to enable all embedding code

  // Cut variables
  Double_t fTrackProbElTPC;
  Double_t fTrackProbPiTPC;
  Double_t fTrackProbKaTPC;
  Double_t fTrackProbPrTPC;
  Bool_t   fTrackProbDeTPC;
  Double_t fTrackProbElTOF;
  Double_t fTrackProbPiTOF;
  Double_t fTrackProbKaTOF;
  Double_t fTrackProbPrTOF;
  Bool_t   fTrackProbDeTOF;

  int fall_reco_jets_w_multiple_matches;
  int fall_reco_jets_w_matches;

  //histograms
  TH1F             * fHistCentrality;            //!<! control histogram for centrality
  TH1F             * fHistImpParam;              //!<! control histogram for impact parameter
  TH1F             * fHistVertex;                //!<! control histogram for vertexZ
  TH1F             * fHist_Is_Good_Inc_Event;    //!<! histogram for good inclusive events
  TH1F             * fHist_Has_Acc_Real_FJ_Jet;    //!<! histogram for good FJ jet events
  TH1F             * fHist_Has_Acc_Real_FJ_Jet_Truth;    //!<! histogram for good FJ truth jet events
  TH1F             * fHist_Has_Acc_Real_EMC_Jet;    //!<! histogram for good jet events
  TH1F             * fHist_Num_Real_FJ_Jets;    //!<! histogram for number of real FJ jets
  TH1F             * fHist_Num_Real_FJ_Jets_Truth;    //!<! histogram for number of real FJ jets
  TH1F             * fHist_Num_Real_EMC_Jets;    //!<! histogram for number of real EMC jets
  TH1F             * fHist_RhoFJ;    //!<! histogram for FJ rho values
  TH1F             * fHist_RhoEMC;    //!<! histogram for EMC rho values
  TH1F             * fHist_RhoEMCReco;    //!<! histogram for EMC reco rho values
  TH1F             * fHist_RhoEMCGen;    //!<! histogram for EMC gen rho values
  TH1F             * fHist_RhoEMCEmb;    //!<! histogram for EMC Emb rho values
  TH1F             * fHist_Filled_UEC_Emb;    //!<! histogram for filled Emb UE cones or not
  TH1F             * fHist_Filled_UEC_Gen;    //!<! histogram for filled Gen UE cones or not
  TH1F             * fHist_Filled_UEC_Rec;    //!<! histogram for filled Rec UE cones or not

  TH3F             * fHistIncTracks_dEdx;        //!<! histogram for inclusive tracks dEdx all eta v some momentum form
  TH2F             * fHistIncTracks_moms;        //!<! histogram for inclusive tracks ptpc to pT
  TH2F             * fHistIncTracks_moms_p;        //!<! histogram for inclusive tracks p to pT
  TH2F             * fHistIncTracks_moms_pTPC_p;        //!<! histogram for inclusive tracks pTPC to p
  TH3F             * fHistIncTracks_kin;        //!<! histogram for inclusive tracks dEdx all eta
  TH2F             * fHistIncTracks_beta;        //!<! histogram for inclusive tracks beta all eta v p
  TH2F             * fHistIncTracks_t0;        //!<! histogram for inclusive tracks t0 all eta v p
  TH3F             * fHistIncTracks_TOFpi_nsigma;        //!<! histogram for inclusive tracks TOF nsigma under pion hypothesis vs pT
  TH3F             * fHistIncTracks_TOFka_nsigma;        //!<! histogram for inclusive tracks TOF nsigma under kaon hypothesis vs pT
  TH3F             * fHistIncTracks_TOFpr_nsigma;        //!<! histogram for inclusive tracks TOF nsigma under proton hypothesis vs pT

  TH3F             * fHistJetTracks_dEdx;        //!<! histogram for jet tracks dEdx all eta v some momentum form
  TH2F             * fHistJetTracks_moms;        //!<! histogram for jet tracks ptpc to pT
  TH2F             * fHistJetTracks_moms_p;        //!<! histogram for jet tracks p to Pt
  TH2F             * fHistJetTracks_moms_pTPC_p;        //!<! histogram for jet tracks pTPC to p
  TH3F             * fHistJetTracks_kin;         //!<! histogram for jet tracks dEdx all eta
  TH2F             * fHistJetTracks_beta;        //!<! histogram for jet tracks beta all eta v p
  TH3F             * fHistJetTracks_TOFpi_nsigma;        //!<! histogram for jet tracks TOF nsigma under pion hypothesis vs pT
  TH3F             * fHistJetTracks_TOFka_nsigma;        //!<! histogram for jet tracks TOF nsigma under kaon hypothesis vs pT
  TH3F             * fHistJetTracks_TOFpr_nsigma;        //!<! histogram for jet tracks TOF nsigma under proton hypothesis vs pT

  TH2F             * fHistBetaExpec_pi; //!<! histogram for inc expected pion beta v pT
  TH2F             * fHistBetaExpec_ka; //!<! histogram for inc expected kaon beta v pT
  TH2F             * fHistBetaExpec_pr; //!<! histogram for inc expected proton beta v pT

  TH2F             * fHistjet_BetaExpec_pi; //!<! histogram for jet expected pion beta v pT
  TH2F             * fHistjet_BetaExpec_ka; //!<! histogram for jet expected kaon beta v pT
  TH2F             * fHistjet_BetaExpec_pr; //!<! histogram for jet expected proton beta v pT

  TH3F             * fHist_pi_mismatch; //!<! histogram for inc pion mismatch v pT
  TH3F             * fHist_ka_mismatch; //!<! histogram for inc kaon mismatch v pT
  TH3F             * fHist_pr_mismatch; //!<! histogram for inc proton mismatch v pT

  TH3F             * fHist_jet_pi_mismatch; //!<! histogram for jet pion mismatch v pT
  TH3F             * fHist_jet_ka_mismatch; //!<! histogram for jet kaon mismatch v pT
  TH3F             * fHist_jet_pr_mismatch; //!<! histogram for jet proton mismatch v pT

  TH3F             * fHist_elExpec_pihyp; //!<! histogram for expected inc electron nsigma under the pion hypothesis v pT
  TH3F             * fHist_muExpec_pihyp; //!<! histogram for expected inc muon nsigma under the pion hypothesis v pT
  TH3F             * fHist_kaExpec_pihyp; //!<! histogram for expected inc kaon nsigma under the pion hypothesis v pT
  TH3F             * fHist_prExpec_pihyp; //!<! histogram for expected inc proton nsigma under the pion hypothesis v pT
  TH3F             * fHist_piExpec_kahyp; //!<! histogram for expected inc pion nsigma under the kaon hypothesis v pT
  TH3F             * fHist_prExpec_kahyp; //!<! histogram for expected inc proton nsigma under the kaon hypothesis v pT
  TH3F             * fHist_piExpec_prhyp; //!<! histogram for expected inc pion nsigma under the proton hypothesis v pT
  TH3F             * fHist_kaExpec_prhyp; //!<! histogram for expected inc kaon nsigma under the proton hypothesis v pT
  TH3F             * fHist_deExpec_prhyp; //!<! histogram for expected inc deuteron nsigma under the proton hypothesis v pT

  TH3F             * fHist_jet_elExpec_pihyp; //!<! histogram for expected jet electron nsigma under the pion hypothesis v pT
  TH3F             * fHist_jet_muExpec_pihyp; //!<! histogram for expected jet muon nsigma under the pion hypothesis v pT
  TH3F             * fHist_jet_kaExpec_pihyp; //!<! histogram for expected jet kaon nsigma under the pion hypothesis v pT
  TH3F             * fHist_jet_prExpec_pihyp; //!<! histogram for expected jet proton nsigma under the pion hypothesis v pT
  TH3F             * fHist_jet_piExpec_kahyp; //!<! histogram for expected jet pion nsigma under the kaon hypothesis v pT
  TH3F             * fHist_jet_prExpec_kahyp; //!<! histogram for expected jet proton nsigma under the kaon hypothesis v pT
  TH3F             * fHist_jet_piExpec_prhyp; //!<! histogram for expected jet pion nsigma under the proton hypothesis v pT
  TH3F             * fHist_jet_kaExpec_prhyp; //!<! histogram for expected jet kaon nsigma under the proton hypothesis v pT
  TH3F             * fHist_jet_deExpec_prhyp; //!<! histogram for expected jet deuteron nsigma under the proton hypothesis v pT

  TH3F             * fHistTOFSigmaExpec_pi; //!<! histogram for expected inc pion TOF Sigma v pT
  TH3F             * fHistTOFSigmaExpec_ka; //!<! histogram for expected inc kaon TOF Sigma v pT
  TH3F             * fHistTOFSigmaExpec_pr; //!<! histogram for expected inc proton TOF Sigma v pT

  TH3F             * fHistjet_TOFSigmaExpec_pi; //!<! histogram for expected jet pion TOF Sigma v pT
  TH3F             * fHistjet_TOFSigmaExpec_ka; //!<! histogram for expected jet kaon TOF Sigma v pT
  TH3F             * fHistjet_TOFSigmaExpec_pr; //!<! histogram for expected jet proton TOF Sigma v pT

  TH3F             * fHistIncTracks_mpi;        //!<! intermediate histogram for inclusive tracks dEdx expected pion mean
  TH3F             * fHistIncTracks_spi;        //!<! intermediate histogram for inclusive tracks dEdx expected pion sigma
  TH3F             * fHistIncTracks_mel;        //!<! intermediate histogram for inclusive tracks dEdx expected electron mean
  TH3F             * fHistIncTracks_sel;        //!<! intermediate histogram for inclusive tracks dEdx expected electron sigma
  TH3F             * fHistIncTracks_mka;        //!<! intermediate histogram for inclusive tracks dEdx expected kaon mean
  TH3F             * fHistIncTracks_ska;        //!<! intermediate histogram for inclusive tracks dEdx expected kaon sigma
  TH3F             * fHistIncTracks_mpr;        //!<! intermediate histogram for inclusive tracks dEdx expected proton mean
  TH3F             * fHistIncTracks_spr;        //!<! intermediate histogram for inclusive tracks dEdx expected proton sigma

  TH2F             * fHistJet_ptsub_v_area;     //!<! histogram for before any cuts, jet pt after bg subtraction vs jet area
  TH3F             * fHistJet_kin;     //!<! histogram for jet ptsub, eta, phi after area cut
  TH2F             * fHistJet_moms;     //!<! histogram for jet pt v jet ptsub after area cut*/

  TH2F             * fHist_pi_DCAxy;       //!<! histogram for Data DCAxy distributions pion v pT
  TH2F             * fHist_pr_DCAxy;       //!<! histogram for Data DCAxy distributions pion v pT

  TH2F             * fHistJet_pi_DCAxy;       //!<! histogram for Data DCAxy distributions pion v pT
  TH2F             * fHistJet_pr_DCAxy;       //!<! histogram for jet Data DCAxy distributions pion v pT

  TH2F             * fHistMCTruth_TrackEff_Den_pi;     //!<! histogram for MC Truth Level Tracking Efficiency Denominator pion v pT, eta
  TH2F             * fHistMCTruth_TrackEff_Den_ka;     //!<! histogram for MC Truth Level Tracking Efficiency Denominator kaon v pT, eta
  TH2F             * fHistMCTruth_TrackEff_Den_pr;     //!<! histogram for MC Truth Level Tracking Efficiency Denominator proton v pT, eta
  TH2F             * fHistMCReco_TrackEff_Num_pi;     //!<! histogram for MC Reconstructed Level Tracking Efficiency Numerator pion v pT, eta
  TH2F             * fHistMCReco_TrackEff_Num_ka;     //!<! histogram for MC Reconstructed Level Tracking Efficiency Numerator kaon v pT, eta
  TH2F             * fHistMCReco_TrackEff_Num_pr;     //!<! histogram for MC Reconstructed Level Tracking Efficiency Numerator proton v pT, eta
  TH2F             * fHistMCReco_MatchEff_Num_pi;     //!<! histogram for MC Reconstructed Level TOF Matching Efficiency Numerator pion v pT, eta
  TH2F             * fHistMCReco_MatchEff_Num_ka;     //!<! histogram for MC Reconstructed Level TOF Matching Efficiency Numerator kaon v pT, eta
  TH2F             * fHistMCReco_MatchEff_Num_pr;     //!<! histogram for MC Reconstructed Level TOF Matching Efficiency Numerator proton v pT, eta

  TH2F             * fHistMCReco_pi_prim_DCAxy;       //!<! histogram for MC Reconstructed Level Primary DCAxy distributions pion v pT
  TH2F             * fHistMCReco_pi_scdweak_DCAxy;    //!<! histogram for MC Reconstructed Level Secondary from Weak Decays DCAxy distributions pion v pT
  TH2F             * fHistMCReco_pr_prim_DCAxy;       //!<! histogram for MC Reconstructed Level Primary DCAxy distributions proton v pT
  TH2F             * fHistMCReco_pr_scdweak_DCAxy;    //!<! histogram for MC Reconstructed Level Secondary from Weak Decays DCAxy distributions proton v pT
  TH2F             * fHistMCReco_pr_scdmat_DCAxy;    //!<! histogram for MC Reconstructed Level Secondary from Material DCAxy distributions proton v pT


  TH2F             * fHist_Jet_MCTruth_TrackEff_Den_pi;     //!<! jet histogram for MC Truth Level Tracking Efficiency Denominator pion v pT, eta
  TH2F             * fHist_Jet_MCTruth_TrackEff_Den_ka;     //!<! jet histogram for MC Truth Level Tracking Efficiency Denominator kaon v pT, eta
  TH2F             * fHist_Jet_MCTruth_TrackEff_Den_pr;     //!<! jet histogram for MC Truth Level Tracking Efficiency Denominator proton v pT, eta
  TH2F             * fHist_Jet_MCReco_TrackEff_Num_pi;     //!<! jet histogram for MC Reconstructed Level Tracking Efficiency Numerator pion v pT, eta
  TH2F             * fHist_Jet_MCReco_TrackEff_Num_ka;     //!<! jet histogram for MC Reconstructed Level Tracking Efficiency Numerator kaon v pT, eta
  TH2F             * fHist_Jet_MCReco_TrackEff_Num_pr;     //!<! jet histogram for MC Reconstructed Level Tracking Efficiency Numerator proton v pT, eta
  TH2F             * fHist_Jet_MCReco_MatchEff_Num_pi;     //!<! jet histogram for MC Reconstructed Level TOF Matching Efficiency Numerator pion v pT, eta
  TH2F             * fHist_Jet_MCReco_MatchEff_Num_ka;     //!<! jet histogram for MC Reconstructed Level TOF Matching Efficiency Numerator kaon v pT, eta
  TH2F             * fHist_Jet_MCReco_MatchEff_Num_pr;     //!<! jet histogram for MC Reconstructed Level TOF Matching Efficiency Numerator proton v pT, eta

  TH3F             * fHist_JetMatch_MCTruth_TrackEff_Den_pi;     //!<! jet match histogram for MC Truth Level Tracking Efficiency Denominator pion v pT
  TH3F             * fHist_JetMatch_MCTruth_TrackEff_Den_ka;     //!<! jet match histogram for MC Truth Level Tracking Efficiency Denominator kaon v pT
  TH3F             * fHist_JetMatch_MCTruth_TrackEff_Den_pr;     //!<! jet match histogram for MC Truth Level Tracking Efficiency Denominator proton v pT
  TH3F             * fHist_JetMatch_MCReco_TrackEff_Num_pi;     //!<! jet match histogram for MC Reconstructed Level Tracking Efficiency Numerator pion v pT
  TH3F             * fHist_JetMatch_MCReco_TrackEff_Num_ka;     //!<! jet match histogram for MC Reconstructed Level Tracking Efficiency Numerator kaon v pT
  TH3F             * fHist_JetMatch_MCReco_TrackEff_Num_pr;     //!<! jet match histogram for MC Reconstructed Level Tracking Efficiency Numerator proton v pT
  TH3F             * fHist_JetMatch_MCReco_MatchEff_Num_pi;     //!<! jet match histogram for MC Reconstructed Level TOF Matching Efficiency Numerator pion v pT
  TH3F             * fHist_JetMatch_MCReco_MatchEff_Num_ka;     //!<! jet match histogram for MC Reconstructed Level TOF Matching Efficiency Numerator kaon v pT
  TH3F             * fHist_JetMatch_MCReco_MatchEff_Num_pr;     //!<! jet match histogram for MC Reconstructed Level TOF Matching Efficiency Numerator proton v pT

  TH2F             * fHistMCT_Jet_ptsub_v_area;     //!<! MC Truth histogram for before any cuts, jet pt after bg subtraction vs jet area
  TH2F             * fHistMCR_Jet_ptsub_v_area;     //!<! MC Reco histogram for before any cuts, jet pt after bg subtraction vs jet area

  TH2F             * fHistMC_Jet_shared_pt_frac;     //!<! fHistMC_Jet_shared_pt_frac
  TH2F             * fHistMC_Jet_deltaR;     //!<! fHistMC_Jet_deltaR

  TH3F             * fHist_TrueJetPtFraction; //!<! fHist_TrueJetPtFraction
  TH3F             * fHist_MatchJetPts; //!<! fHist_MatchJetPts
  TH3F             * fHist_MatchJetEtas; //!<! fHist_MatchJetEtas
  TH3F             * fHist_MatchJetDeltaRs; //!<! fHist_MatchJetDeltaRs

  TH3F             * fHist_JERS_PbPbunid_tru; //!<! fHist_JERS_PbPbunid_tru
  TH3F             * fHist_JERS_PbPbunid_det; //!<! fHist_JERS_PbPbunid_det
  TH3F             * fHist_JERS_PbPbunid_truhyb; //!<! fHist_JERS_PbPbunid_truhyb
  TH3F             * fHist_JERS_PbPbunid_dethyb; //!<! fHist_JERS_PbPbunid_dethyb

  TH3F             * fHist_JERS_truPythia; //!<! fHist_JERS_truPythia
  TH3F             * fHist_JERS_truPythia_pi; //!<! fHist_JERS_truPythia_pi
  TH3F             * fHist_JERS_truPythia_ka; //!<! fHist_JERS_truPythia_ka
  TH3F             * fHist_JERS_truPythia_pr; //!<! fHist_JERS_truPythia_pr

  TH3F             * fHist_JERS_detPythia; //!<! fHist_JERS_detPythia
  TH3F             * fHist_JERS_detPythia_pi; //!<! fHist_JERS_detPythia_pi
  TH3F             * fHist_JERS_detPythia_ka; //!<! fHist_JERS_detPythia_ka
  TH3F             * fHist_JERS_detPythia_pr; //!<! fHist_JERS_detPythia_pr

  TH3F             * fHist_JERS_truhybPythia; //!<! fHist_JERS_truhybPythia
  TH3F             * fHist_JERS_truhybPythia_pi; //!<! fHist_JERS_truhybPythia_pi
  TH3F             * fHist_JERS_truhybPythia_ka; //!<! fHist_JERS_truhybPythia_ka
  TH3F             * fHist_JERS_truhybPythia_pr; //!<! fHist_JERS_truhybPythia_pr

  TH3F             * fHist_JERS_dethybPythia; //!<! fHist_JERS_dethybPythia
  TH3F             * fHist_JERS_dethybPythia_pi; //!<! fHist_JERS_dethybPythia_pi
  TH3F             * fHist_JERS_dethybPythia_ka; //!<! fHist_JERS_dethybPythia_ka
  TH3F             * fHist_JERS_dethybPythia_pr; //!<! fHist_JERS_dethybPythia_pr

  TH1F             * fHist_PC_spectra_Pythia; //!<! fHist_PC_spectra_Pythia
  TH1F             * fHist_PC_spectra_Pythia_pi; //!<! fHist_PC_spectra_Pythia_pi
  TH1F             * fHist_PC_spectra_Pythia_ka; //!<! fHist_PC_spectra_Pythia_ka
  TH1F             * fHist_PC_spectra_Pythia_pr; //!<! fHist_PC_spectra_Pythia_pr

  TH2F             * fHist_det_spectra_Pythia; //!<! fHist_det_spectra_Pythia
  TH2F             * fHist_det_spectra_Pythia_pi; //!<! fHist_det_spectra_Pythia_pi
  TH2F             * fHist_det_spectra_Pythia_ka; //!<! fHist_det_spectra_Pythia_ka
  TH2F             * fHist_det_spectra_Pythia_pr; //!<! fHist_det_spectra_Pythia_pr

  TH2F             * fHist_part_spectra_Pythia; //!<! fHist_part_spectra_Pythia
  TH2F             * fHist_part_spectra_Pythia_pi; //!<! fHist_part_spectra_Pythia_pi
  TH2F             * fHist_part_spectra_Pythia_ka; //!<! fHist_part_spectra_Pythia_ka
  TH2F             * fHist_part_spectra_Pythia_pr; //!<! fHist_part_spectra_Pythia_pr

  AliFJWrapper* fFastJetWrapper;            //!<! signal data FJwrapper
  AliFJWrapper* fFastJetWrapperBG;            //!<! background data FJwrapper
  AliFJWrapper* fFastJetWrapper_Rec;            //!<! signal MC reco FJwrapper
  AliFJWrapper* fFastJetWrapperBG_Rec;            //!<! background MC reco FJwrapper
  AliFJWrapper* fFastJetWrapper_Gen;            //!<! signal MC gen FJwrapper
  AliFJWrapper* fFastJetWrapperBG_Gen;            //!<! background MC gen FJwrapper

  ClassDef(AliAnalysisTaskJetHadroAOD, 7);

};

#endif
