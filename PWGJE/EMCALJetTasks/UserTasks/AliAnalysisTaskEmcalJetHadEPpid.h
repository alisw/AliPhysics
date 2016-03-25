#ifndef AliAnalysisTaskEmcalJetHadEPpid_h
#define AliAnalysisTaskEmcalJetHadEPpid_h

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

//-------------------------------------------------------------------------
// 1) Analysis Task to perform Jet-Hadron Correlations
// 2) Event plane dependence task.
// 3) performs event plane resolution calculation
// 4) does PID of the associated pi/k/p hadrons
//
// Author: Joel Mazer (joel.mazer@cern.ch)
//-------------------------------------------------------------------------

// root classes
class TClonesArray;
class TF1;
class TH1;
class TH2;
class TH3;
class THnSparse;
class TProfile;
class TList;
class TLorentzVector;
class TGraph;

// AliROOT classes
class AliEventPoolManager;
class AliLocalRhoParameter;
class AliEMCALTrack;
class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliESDtrack;
class AliESDtrackCuts;

// container classes
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;

// includes
#include <AliAnalysisTaskEmcalJet.h>
#include <AliEmcalJet.h>
#include <AliVEvent.h>
#include <AliVTrack.h>
#include <AliVCluster.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TRandom3.h>
#include <AliLog.h>
#include <TArrayD.h>
#include "AliESDtrackCuts.h"

// Local Rho includes
#include "AliAnalysisTaskLocalRho.h"
#include "AliLocalRhoParameter.h"

// PID includes
#include "AliPIDResponse.h"

#include "AliAnalysisFilter.h"

class AliAnalysisTaskEmcalJetHadEPpid : public AliAnalysisTaskEmcalJet {
 public:
  enum detectorType       { kTPC, kVZEROA, kVZEROC, kVZEROComb, kFixedEP};  // detector that was used for event plane
  AliAnalysisTaskEmcalJetHadEPpid();
  AliAnalysisTaskEmcalJetHadEPpid(const char *name);
  //virtual ~AliAnalysisTaskEmcalJetHadEPpid() {}
  virtual ~AliAnalysisTaskEmcalJetHadEPpid();

  virtual void            UserCreateOutputObjects();
  // THnSparse Setup
  virtual THnSparse*      NewTHnSparseF(const char* name, UInt_t entries);
  virtual void            GetDimParams(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
  virtual THnSparse*      NewTHnSparseFPID(const char* name, UInt_t entries);
  virtual void            GetDimParamsPID(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
  virtual THnSparse*      NewTHnSparseFCorr(const char* name, UInt_t entries);
  virtual void            GetDimParamsCorr(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);

        // note that the cdf of the chisquare distribution is the normalized lower incomplete gamma function
//        /* inline */    static Double_t ChiSquareCDF(Int_t ndf, Double_t x) { return TMath::Gamma(ndf/2., x/2.); }
//        /* inline */    static Double_t ChiSquare(TH1& histo, TF1* func) {
            // evaluate the chi2 using a poissonian error estimate on bins
//            Double_t chi2(0.);
/*
            for(Int_t i(0); i < histo.GetXaxis()->GetNbins(); i++) {
                if(histo.GetBinContent(i+1) <= 0.) continue;
                chi2 += TMath::Power((histo.GetBinContent(i+1)-func->Eval(histo.GetXaxis()->GetBinCenter(1+i))), 2)/histo.GetBinContent(i+1);
            }
           return chi2;
        }
*/

  // set a bunch of histogram switches up
  void                    SetPlotGlobalRho(Bool_t g)            { doPlotGlobalRho = g; } // plot global rho switch
  void                    SetVariableBinning(Bool_t v)          { doVariableBinning = v; } // do variable binning switch
  void	                  SetvarbinTHnSparse(Bool_t vb)         { dovarbinTHnSparse = vb; } // variable THnSparse bin switch
  void	                  SetallpidAXIS(Bool_t allAXIS)		{ allpidAXIS = allAXIS; } // fill all PID sparse axis's
  void	                  SetmakeQAhistos(Bool_t QAhist)        { makeQAhistos = QAhist; } // make QA histos  
  void	                  SetmakeBIAShistos(Bool_t BIAShist)    { makeBIAShistos = BIAShist; } // make bias histos
  void	                  SetmakeextraCORRhistos(Bool_t Xhist)  { makeextraCORRhistos = Xhist; } // make extra correlations histos
  void	                  SetoldJEThadhistos(Bool_t oldJH)      { makeoldJEThadhistos = oldJH; } // make older JH histos for comparison

  // set data, detectors type, and PID and PID w bias switches
  void	                  SetcutType(TString cut)               { fcutType = cut; }    // EMCAL / TPC acceptance cut
  void                    SetdoPID(Bool_t p)                    { doPID = p; }   // do PID switch
  void 	                  SetdoPIDtrackBIAS(Bool_t PIDbias)     { doPIDtrackBIAS = PIDbias; } // do PID track bias switch
  void                    SetdoaltPIDbinning(Bool_t altPIDbin)  { doaltPIDbinning = altPIDbin; } // alternate PID binning (fewer bins - TOF focus)

  // esd track cuts setters
  void                    SetTrackCuts(AliESDtrackCuts *cuts)                      { fesdTrackCuts = cuts; }

  // reference of detector for event plane resolution
  void                    SetReferenceDetector(detectorType type)         {fDetectorType = type; }

  // set soft track min/max
  void                    SetSoftTrackMinMaxPt_ep(Float_t min, Float_t max)          {fSoftTrackMinPt_ep = min; fSoftTrackMaxPt_ep = max;}
  void                    SetExcludeLeadingJetsFromFit(Float_t n)         {fExcludeLeadingJetsFromFit = n; }

  // set centrality classes up
  void                    SetCentralityClasses(TArrayD* c)      {fCentralityClasses = c;}

  // set Chi2 for VZERO A and C
  void                    SetChi2VZEROA(TArrayD* a)                       { fChi2A = a;}
  void                    SetChi2VZEROC(TArrayD* a)                       { fChi2C = a;}
  void                    SetChi3VZEROA(TArrayD* a)                       { fChi3A = a;}
  void                    SetChi3VZEROC(TArrayD* a)                       { fChi3C = a;}
  void                    SetUseChiWeightForVZERO(Bool_t w)               { fUseChiWeightForVZERO = w; }
  void                    ReadVZEROCalibration2011h();
  Int_t                   GetVZEROCentralityBin() const;
  AliEmcalJet*            GetLeadingJet(AliLocalRhoParameter* localRho = 0x0);

  // switch for Event Plane Resolution analysis
  void                    SetdoEventPlaneRes(Bool_t depr)       { doEventPlaneRes = depr; } // do EP res switch

  // give comments setter
  void	                  SetdoComments(Bool_t comm)            { doComments = comm; } // give comment switch

  // setter switch for flavour jet analysis
  void 	                  SetFlavourJetAnalysis(Bool_t flj)     { doFlavourJetAnalysis = flj; } // set on flavour jet analysis
  virtual void	          SetJETFlavourTag(Int_t fltag)         { fJetFlavTag = fltag; } // set manual tag #

  // setter for beamtype (needed for UserCreateObjects section)
  virtual void	          SetCollType(BeamType bm) { fBeam = bm; } // set beamtype 

  // getters
  TString                 GetLocalRhoName() const		{return fLocalRhoName; }

  // set names of some objects
  virtual void            SetLocalRhoName(const char *ln)       { fLocalRhoName = ln; }
  virtual void            SetTracksName(const char *tn)         { fTracksName = tn; }
  virtual void	          SetTracksNameME(const char *MEtn)     { fTracksNameME = MEtn; }
  virtual void            SetJetsName(const char *jn)           { fJetsName = jn; }
  virtual void            SetCaloClustersName(const char *cn)   { fCaloClustersName=cn; }

  // bias and cuts - setters
  virtual void            SetAreaCut(Double_t a)                { fAreacut    = a; }
  virtual void            SetTrkBias(Double_t b)                { fTrkBias    = b; }  //require a track with pt > b in jet
  virtual void            SetClusBias(Double_t b)               { fClusBias   = b; }  //require a cluster with pt > b in jet
  virtual void            SetTrkEta(Double_t e)                 { fTrkEta   = e; }  //eta range of the associated tracks
  virtual void            SetJetPtcut(Double_t jpt)             { fJetPtcut = jpt; } // jet pt cut
  virtual void	          SetJetRad(Double_t jrad)             	{ fJetRad = jrad; } // jet radius 
  virtual void 	          SetConstituentCut(Double_t constCut)  { fConstituentCut = constCut; } // constituent Cut

  // eta and phi limits of jets - setters
  virtual void            SetJetEta(Double_t emin, Double_t emax)  { fEtamin = emin; fEtamax = emax; }
  virtual void            SetJetPhi(Double_t pmin, Double_t pmax)  { fPhimin = pmin; fPhimax = pmax; }

  // event mixing - setters
  virtual void            SetEventMixing(Int_t yesno)	       { fDoEventMixing=yesno; }
  virtual void	          SetMixingTracks(Int_t tracks)	       { fMixingTracks = tracks; }
  virtual void            SetNMixedTr(Int_t nmt)               { fNMIXtracks = nmt; }
  virtual void            SetNMixedEvt(Int_t nme)              { fNMIXevents = nme; }

  // event trigger/mixed selection - setters
  virtual void            SetTriggerEventType(UInt_t te)       { fTriggerEventType = te; }
  virtual void            SetMixedEventType(UInt_t me)         { fMixingEventType = me; }
  virtual void            SetCentBinSize(Int_t centbins)       { fCentBinSize = centbins; }
  virtual void            SetReduceStatsCent(Int_t red)        { fReduceStatsCent = red; }

  // efficiency correction setter
  void                    SetDoEffCorr(Int_t effcorr)          { fDoEffCorr = effcorr; }

  // use local rho to correct jet pt in correlation sparses
  void                    SetCorrectJetPt(Bool_t cpt)           { fcorrJetPt = cpt; }

  // jet container - setters
  void SetContainerAllJets(Int_t c)         { fContainerAllJets      = c;}
  void SetContainerPIDJets(Int_t c)         { fContainerPIDJets      = c;}

protected:
  // functions 
  void	                 ExecOnce();
  virtual Bool_t         Notify();
  Bool_t                 Run();
  virtual void           Terminate(Option_t *); 
  virtual Int_t          AcceptMyJet(AliEmcalJet *jet);   // applies basic jet tests/cuts before accepting
  virtual Int_t          GetCentBin(Double_t cent) const; // centrality bin of event
  Double_t               RelativePhi(Double_t mphi,Double_t vphi) const; // relative jet track angle
  Double_t               RelativeEPJET(Double_t jetAng, Double_t EPAng) const;  // relative jet event plane angle
  virtual Int_t          GetEtaBin(Double_t eta) const;      // eta bins
  virtual Int_t          GetpTjetBin(Double_t pt) const;     // jet pt bins
  virtual Int_t          GetpTtrackBin(Double_t pt) const;   // track pt bins
  virtual Int_t          GetzVertexBin(Double_t zVtx) const; // zVertex bin
  void                   SetfHistPIDcounterLabels(TH1* fHistPID) const;  // PID counter
  void	                 SetfHistQAcounterLabels(TH1* h) const; // QA counter
  void                   SetfHistEvtSelQALabels(TH1* h) const; // Event Selection Counter
  //virtual Int_t                AcceptFlavourJet(AliEmcalJet *jet, Int_t NUM, Int_t NUM2, Int_t NUM3); // flavour jet acceptor
  virtual Int_t	         AcceptFlavourJet(AliEmcalJet *jet, Int_t NUM); // flavour jet acceptor
  Double_t               EffCorrection(Double_t trkETA, Double_t trkPT, Int_t effswitch) const; // efficiency correction function
  static Double_t        CalculateEventPlaneChi(Double_t res);
  void                   CalculateEventPlaneVZERO(Double_t vzero[2][2]) const;
  void                   CalculateEventPlaneCombinedVZERO(Double_t* comb) const;
  void                   CalculateEventPlaneTPC(Double_t* tpc);
  void                   CalculateEventPlaneResolution(Double_t vzero[2][2], Double_t* vzeroComb, Double_t* tpc);
  TH1*                   FillEventTriggerQA(TH1* h, UInt_t t); // filled event trigger QA plots

  // parameters of detector to cut on for event
  Double_t               fPhimin;                  // phi min
  Double_t               fPhimax;                  // phi max
  Double_t               fEtamin;                  // eta min
  Double_t               fEtamax;                  // eta max
  Double_t               fAreacut;                 // area cut
  Double_t               fTrkBias;                 // track bias
  Double_t               fClusBias;                // cluster bias
  Double_t               fTrkEta;                  // eta min/max of tracks
  Double_t               fJetPtcut;                // jet pt to cut on for correlations
  Double_t               fJetRad;                  // jet radius
  Double_t               fConstituentCut;          // jet constituent cut

  // esd track cuts
  AliESDtrackCuts       *fesdTrackCuts;            // esdTrackCuts

  // detector type
  detectorType           fDetectorType;            // type of detector

  // used for event plane resolution calculation
  Float_t                fSoftTrackMinPt_ep;        // min pt for soft tracks
  Float_t                fSoftTrackMaxPt_ep;        // max pt for soft tracks
  Int_t                  fNAcceptedTracks;       //! number of accepted tracks
  AliEmcalJet*           fLeadingJet;            //! leading jet
  Float_t                fExcludeLeadingJetsFromFit;    // exclude n leading jets from fit
  TArrayD*               fCentralityClasses;     //-> centrality classes (maximum 10)
  Int_t                  fInCentralitySelection; //! centrality bin

  // event mixing
  Int_t	         fDoEventMixing;
  Int_t	         fMixingTracks;
  Int_t          fNMIXtracks;
  Int_t          fNMIXevents;
  Int_t          fCentBinSize; // centrality bin size of mixed event pools
  Int_t          fReduceStatsCent; // bins to use for reduced statistics of sparse

  // event selection types
  UInt_t         fTriggerEventType;
  UInt_t         fMixingEventType;

  // efficiency correction
  Int_t          fDoEffCorr;

  // correct jet pt
  Bool_t         fcorrJetPt;

  // switches for plots
  Bool_t         doPlotGlobalRho;
  Bool_t         doVariableBinning;
  Bool_t         dovarbinTHnSparse;
  Bool_t         makeQAhistos;
  Bool_t         makeBIAShistos;
  Bool_t         makeextraCORRhistos; 
  Bool_t         makeoldJEThadhistos;
  Bool_t         allpidAXIS;

  // Cut type (EMCAL/TPC acceptance)
  TString        fcutType;

  // switches for PID
  Bool_t         doPID;
  Bool_t         doPIDtrackBIAS;
  Bool_t         doaltPIDbinning;

  // do EP resolution switch
  Bool_t         doEventPlaneRes;

  // do comment switch
  Bool_t         doComments;

  // do flavour jet analysis switch, and set flavour jet tag
  Bool_t         doFlavourJetAnalysis;
  Int_t	         fJetFlavTag;

  // beam type
  BeamType       fBeam;

  // local rho value
  Double_t       fLocalRhoVal;

  // object names
  TString        fTracksName; // name of track collection (for signal events)
  TString        fTracksNameME; // name of mixed event track collection
  TString        fJetsName; // name of jet collection
  TString        fCaloClustersName; // name of Calo Cluster collection

  // event counter
  Int_t	         event;

  // boolean functions for PID
  Bool_t         isPItpc, isKtpc, isPtpc;
  Bool_t         isPIits, isKits, isPits;
  Bool_t         isPItof, isKtof, isPtof;

  // event pool
  TObjArray      *CloneAndReduceTrackList(TObjArray* tracks);
  AliEventPoolManager   *fPoolMgr;//!  // event pool Manager object

  // PID
  AliPIDResponse	*fPIDResponse;   // PID response object
  AliTPCPIDResponse	*fTPCResponse;   // TPC pid response object

 private:
  // needed for PID, track objects
  AliESDEvent       *fESD;//!         // ESD object
  AliAODEvent	    *fAOD;//!         // AOD object
  AliVEvent         *fVevent;//!      // Vevent object

  TH1F              *fHistEventQA;//!
  TH1F              *fHistEventSelectionQA;//!
  TH1F              *fHistEventSelectionQAafterCuts;//!

  TH2F              *fHistCentZvertGA;//!
  TH2F              *fHistCentZvertJE;//!
  TH2F              *fHistCentZvertMB;//!
  TH2F              *fHistCentZvertAny;//!

  TH2F              *fHistTPCdEdX;//!
  TH2F	            *fHistITSsignal;//!
//  TH2F              *fHistTOFsignal;//!

  TH2F              *fHistRhovsCent;//!
  TH2F              *fHistNjetvsCent;//! number of jets versus Centrality
  TH2F              *fHistJetPtvsTrackPt[6];//!
  TH1F              *fHistTrackPt[6];//!
  TH1F              *fHistEP0[6];//!
  TH1F              *fHistEP0A[6];//!
  TH1F              *fHistEP0C[6];//!
  TH2F              *fHistEPAvsC[6];//!
  TH1F              *fHistJetPtcorrGlRho[6];//!
  TH2F              *fHistJetPtvsdEP[6];//!
  TH2F              *fHistJetPtvsdEPBias[6];//!
  TH2F              *fHistRhovsdEP[6];//!
  TH3               *fHistJetEtaPhiPt[6];//!
  TH3               *fHistJetEtaPhiPtBias[6];//!
  TH2F              *fHistJetPtArea[6];//!
  TH2F              *fHistJetPtAreaBias[6];//!
  TH2F              *fHistJetPtNcon[6];//!
  TH2F              *fHistJetPtNconBias[6];//!
  TH2F              *fHistJetPtNconCh[6];//!
  TH2F              *fHistJetPtNconBiasCh[6];//!
  TH2F              *fHistJetPtNconEm[6];//!
  TH2F              *fHistJetPtNconBiasEm[6];//!
  TH1F              *fHistJetHaddPhiINcent[6];//!
  TH1F              *fHistJetHaddPhiOUTcent[6];//!
  TH1F              *fHistJetHaddPhiMIDcent[6];//!

  TH1               *fHistMult;//!
  TH1               *fHistJetPhi;//!
  TH1               *fHistTrackPhi;//!
  TH1               *fHistLocalRhoJetpt;//!
  TH1               *fHistJetHaddPhiIN;//!
  TH1               *fHistJetHaddPhiOUT;//!
  TH1               *fHistJetHaddPhiMID;//!
  TH1               *fHistJetHaddPhiBias;//!
  TH1               *fHistJetHaddPhiINBias;//!
  TH1               *fHistJetHaddPhiOUTBias;//!
  TH1               *fHistJetHaddPhiMIDBias;//!

  TH1               *fHistMEdPHI;//! // phi distrubtion of mixed events
  TH1               *fHistTrackPtallcent;//!

  TH2               *fHistJetEtaPhi;//!
  TH2               *fHistTrackEtaPhi[4][7];//!

  TH1               *fHistJetHadbindPhi[9];//! 
  TH1               *fHistJetHadbindPhiIN[9];//! 
  TH1               *fHistJetHadbindPhiMID[9];//! 
  TH1               *fHistJetHadbindPhiOUT[9];//! 
  TH2               *fHistJetHEtaPhi;//!

  TH1               *fHistJetPt[6];//!
  TH1               *fHistJetPtBias[6];//!
  TH1               *fHistJetPtTT[6];//!
  TH2               *fHistAreavsRawPt[6];//!
  TH2               *fHistJetH[6][5][3];//!
  TH2               *fHistJetHBias[6][5][3];//!
  TH2               *fHistJetHTT[6][5][3];//!
  TH2F              *fHistSEphieta;//! // single events phi-eta distributions
  TH2F              *fHistMEphieta;//! // mixed events phi-eta distributions
  TH1F              *fHistJetHaddPHI;//!

  // more QA histos
  TH3               *fHistClusEtaPhiEnergy;//!
  TH1               *fHistClusEnergy;//!
  TH1               *fHistClusofJetEnergy;//!
  TH1               *fHistJetNClusterConstit;//!
  TH1               *fHistJetNTrackConstit;//!
  TH1               *fHistJetNConstit;//!

  // PID status histo's
  TH1               *fHistPID;//!

  // THn Sparse's
  THnSparse             *fhnPID;//!          // PID sparse
  THnSparse             *fhnMixedEvents;//!  // mixed events matrix
  THnSparse             *fhnJH;//!           // jet hadron events matrix
  THnSparse             *fhnCorr;//!         // sparse to get # jet triggers

  // EP resoltuion profiles and chi2 array
  TProfile              *fProfV2Resolution[10];//! resolution parameters for v2
  TProfile              *fProfV3Resolution[10];//! resolution parameters for v3
  TProfile              *fProfV4Resolution[10];//! resolution parameters for v4
  TProfile              *fProfV5Resolution[10];//! resolution parameters for v5
  TArrayD*              fChi2A;                         // chi vs cent for vzero A ep_2
  TArrayD*              fChi2C;                         // chi vs cent for vzero C ep_2
  TArrayD*              fChi3A;                         // chi vs cent for vzero A ep_3
  TArrayD*              fChi3C;                         // chi vs cent for vzero C ep_3
  Bool_t                fUseChiWeightForVZERO;          // use chi weight for vzero

  // save containers in clones array (object)
  TClonesArray          *fTracksFromContainer;       //!jets

  // container objects
  AliJetContainer       *fJetsCont;                //!Jets
  AliParticleContainer  *fTracksCont;              //!Tracks
  AliClusterContainer   *fCaloClustersCont;        //!Clusters

  // container specifier
  Int_t	                fContainerAllJets;  // number of container with all full jets
  Int_t	                fContainerPIDJets;  // number of container with full jets meeting Pt cut (for PID)

  //TObjArray                  *fTrgJet;                       //!jets
// ***********************************************************
   
  //Declare it private to avoid compilation warning
  AliAnalysisTaskEmcalJetHadEPpid(const AliAnalysisTaskEmcalJetHadEPpid & g) ; // cpy ctor

  AliAnalysisTaskEmcalJetHadEPpid& operator=(const AliAnalysisTaskEmcalJetHadEPpid&); // not implemented
  ClassDef(AliAnalysisTaskEmcalJetHadEPpid, 4); // Emcal jet hadron PID - Event plane dependence
};
#endif
