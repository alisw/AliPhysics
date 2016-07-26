#ifndef AliAnalysisTaskEmcalTriggerPatchJetMatch_h
#define AliAnalysisTaskEmcalTriggerPatchJetMatch_h
//#ifndef ALIANALYSISTASKEMCALTRIGGERPATCHJETMATCH_H
//#define ALIANALYSISTASKEMCALTRIGGERPATCHJETMATCH_H

//-------------------------------------------------------------------------
// 1) Analysis task to identify the jet whose cluster fired the trigger patch
// 2) perform some QA on the patch / cluster / jet
// 3) and pass the saved out collection of jet(s) to other tasks
//
// currently set up for GA trigger
//
// Author: Joel Mazer (joel.mazer@cern.ch)
//-------------------------------------------------------------------------
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TH1;
class TH2;
class TH3;
class TH3F;
class TProfile;
class TClonesArray;
class TArrayI;
class AliEMCALTriggerPatchInfo;

#include <TRef.h>
#include <TBits.h>
#include <TMath.h>
#include <AliVEvent.h>
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerSetupInfo.h"
#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalTriggerPatchJetMatch : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskEmcalTriggerPatchJetMatch();
  AliAnalysisTaskEmcalTriggerPatchJetMatch(const char *name);
  virtual ~AliAnalysisTaskEmcalTriggerPatchJetMatch();

  enum MainPatchType {
    kManual = 0,    //just select highest energy patch in array
    kEmcalJet = 1   //use functionality of AliAnalysisTaskEmcal
  };

//  void                        ExecOnce();
  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);
  Bool_t                      SelectEvent();              //decides if event is used for analysis
  void                        FillTriggerPatchHistos();
  void                        SetAttachToEvent(Bool_t a)                      {fAttachToEvent = a;}

  //Setters
  void SetDebug(Int_t d)                    { fDebug = d;}
  void SetTriggerClass(const char *n)       { fTriggerClass = n; }
  void SetNFastorPatch(Int_t i)             { fNFastOR = i;}
 
  // containers
  void SetJetContainer(Int_t c)            { fJetContainer      = c;}

  void SetMainPatchType(MainPatchType t)    { fMainPatchType      = t;}//!
  void SetMainTriggerTypeCat(TriggerCategory cat, Bool_t b) {fMainTrigCat = cat; fMainTrigSimple = b;}

  // cuts and biases
  virtual void            SetJetPtcut(Double_t jpt)             { fJetPtCut = jpt; } // jet pt cut
  virtual void            SetPatchEcut(Double_t pE)             { fPatchECut = pE; } // jet pt cut
  virtual void            SetTrkBias(Double_t b)                { fTrkBias    = b; }  //require a track with pt > b in jet
  virtual void            SetClusBias(Double_t b)               { fClusBias   = b; }  //require a cluster with pt > b in jet

  // give comments setter and various switches
  void                    SetdoComments(Bool_t comm)               { doComments = comm; } // give comment switch
  void                    SetUseALLrecalcPatches(Bool_t useall)    { fUseALLrecalcPatches = useall; } // use all firing (offline) patches

  virtual void            SetJetTriggeredEventname(const char *jcol)           { fJetTriggeredEventname = jcol; }
  TString                 GetJetTriggeredEvent() const          { return fJetTriggeredEventname; }
  virtual void            SetCaloClustersName(const char *cn)   { fCaloClustersName=cn; }

  Int_t    GetLeadingCellId(const AliVCluster *clus) const;
  Double_t GetEnergyLeadingCell(const AliVCluster *clus) const;
  Double_t GetECross(Int_t absID) const;

  Double_t GetZ(const AliVParticle *trk, const AliEmcalJet *jet)       const;
  Double_t GetZ(Double_t trkPx, Double_t trkPy, Double_t trkPz, Double_t jetPx, Double_t jetPy, Double_t jetPz) const;

 protected:
  void                        ExecOnce();
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;
  Float_t                     RelativeEP(Double_t objAng, Double_t EPAng) const;
  Bool_t                      TestFilterBit(Int_t trigBit, UInt_t bitJetTrig) const {return (Bool_t) ((trigBit & bitJetTrig) != 0);}

  // Trigger bit - do i need?????
  void                        ExtractMainPatch();//!
  TH1*                        FillTriggerPatchQA(TH1* h, UInt_t t, AliEMCALTriggerPatchInfo* fPatch); // filled trigger patch QA
  TH1*                        FillEventTriggerQA(TH1* h, UInt_t t); // fill event trigger QA
  Bool_t                      CorrelateToTrigger(Double_t etaclust, Double_t phiclust, TList *triggerpatches) const;

 private:
  Bool_t             fDebug;                 // debug level
  Bool_t             fAttachToEvent;         // attach local rho to the event
  TString            fTriggerClass;          // trigger class to analyze EJ or GA patches
  Int_t              fJetContainer;          // number of container of jets
  Double_t           fMaxPatchEnergy;        // energy of patch with largest energy (offline)
  Double_t           fMaxPatchADCEnergy;     // energy of patch with largest energy from online ADC
  Int_t              fTriggerType;           // trigger type
  Int_t              fNFastOR;               // size of trigger patch fNFastORxfNFastOR
  TriggerCategory    fMainTrigCat;           // trigger category for main trigger
  TriggerCategory    fTriggerCategory;       // trigger category
  Bool_t             fMainTrigSimple;        // use offline trigger instead of online from AliAnalysisTaskEmcal::GetMainTriggerPatch
  Double_t	     fJetPtCut;	             // jet pt to cut on for correlations
  Double_t           fPatchECut;             // patch energy cut
  Double_t           fTrkBias;               // track bias
  Double_t           fClusBias;              // cluster bias
  Bool_t             doComments;             // summary (debugging comments)
  Bool_t             fUseALLrecalcPatches;   // use all/ just max recalulated (offline) patches
  TString	     fJetTriggeredEventname; // name of jet that triggered event collection
  TString            fCaloClustersName;      // name of Calo Cluster collection

  AliEMCALTriggerPatchInfo      *fMaxPatch;//!                           main patch
  MainPatchType                 fMainPatchType;//!                       method to select main patch

  TH1F     *fhNEvents;                         //! Histo number of events
  TProfile *fhTriggerbit;                      //! histogram containing the triggerbit (fOfflineTriggerMask)
  TH2F     *fHistRhovsCent;                    //! rho vs. centrality
  TH3F     *fh3PtEtaPhiTracks;                 //! pt,eta,phi of tracks at vertex
  TH3F     *fh3PtEtaPhiTracksOnEmcal;          //! pt,eta,phi of tracks at Emcal surface
  TH3F     *fh3PtEtaPhiTracksToProp;           //! pt,eta,phi of tracks at vertex
  TH3F     *fh3PtEtaPhiTracksProp;             //! pt,eta,phi of tracks at vertex
  TH3F     *fh3PtEtaPhiTracksNoProp;           //! pt,eta,phi of tracks at vertex
  TH2F     *fh2CentPtJet;                      //! cent, pt of jets
  TH3F     *fh3PtEtaPhiJet;                    //! pt,eta,phi of jets
  TH2F     *fh2NJetsPt;                        //! NJets per event vs pT,jet
  TH3F     *fh3PtEtaAreaJet;                   //! pt,eta,area of jet
  TH2F     *fh2PtNConstituentsCharged;         //! pt, # charged jet constituents
  TH2F     *fh2PtNConstituents;                //! pt, # jet constituents
  TH2F     *fh2PtMeanPtConstituentsCharged;    //! pt, <pt> charged constituents
  TH2F     *fh2PtMeanPtConstituentsNeutral;    //! pt, <pt> neutral constituents
  TH2F     *fh2PtNEF;                          //! pt, NEF (neutral energy fraction)
  TH3F     *fh3NEFEtaPhi;                      //! NEF, eta, phi
  TH2F     *fh2NEFNConstituentsCharged;        //! NEF, # charged jet constituents
  TH2F     *fh2NEFNConstituentsNeutral;        //! NEF, # neutral jet constituents
  TH2F     *fh2Ptz;                            //! pt, z=pT,h,proj/p,jet jet
  TH2F     *fh2PtLeadJet1VsLeadJet2;           //! correlation between leading jet of the two branches
  TH3F     *fh3EEtaPhiCluster;                 //! cluster E, eta, phi
  TH3F     *fh3PtLeadJet1VsPatchEnergy;        //! leading jet energy vs leading patch energy vs jet trigger (J1/J2)
  TH3F     *fh3PtLeadJet2VsPatchEnergy;        //! leading jet energy vs leading patch energy vs jet trigger (J1/J2)
  TH3F     *fh3PtLeadJet1PatchEnergyVZEROAmp;  //! leading jet energy vs leading patch energy vs VZERO amplitude
  TH3F     *fh3PtLeadJet1RawPatchEnergyVZEROAmp;  //! leading jet energy vs online leading patch energy vs VZERO amplitude
  TH3F     *fh3PatchEnergyEtaPhiCenterJ1;      //! patch energy vs eta, phi at center of patch, high threshold
  TH3F     *fh3PatchEnergyEtaPhiCenterJ2;      //! patch energy vs eta, phi at center of patch, low threshold
  TH3F     *fh3PatchEnergyEtaPhiCenterJ1J2;    //! patch energy vs eta, phi at center of patch, low + high threshold
  TH3F     *fh3PatchADCEnergyEtaPhiCenterJ1;   //! patch ADC energy vs eta, phi at center of patch, high threshold
  TH3F     *fh3PatchADCEnergyEtaPhiCenterJ2;   //! patch ADC energy vs eta, phi at center of patch, low threshold
  TH3F     *fh3PatchADCEnergyEtaPhiCenterJ1J2; //! patch ADC energy vs eta, phi at center of patch, low + high threshold
  TH3F     *fh3PatchEnergyEtaPhiCenterG1;      //! patch energy vs eta, phi at center of patch, high threshold
  TH3F     *fh3PatchEnergyEtaPhiCenterG2;      //! patch energy vs eta, phi at center of patch, low threshold
  TH3F     *fh3PatchEnergyEtaPhiCenterG1G2;    //! patch energy vs eta, phi at center of patch, low + high threshold
  TH3F     *fh3PatchADCEnergyEtaPhiCenterG1;   //! patch ADC energy vs eta, phi at center of patch, high threshold
  TH3F     *fh3PatchADCEnergyEtaPhiCenterG2;   //! patch ADC energy vs eta, phi at center of patch, low threshold
  TH3F     *fh3PatchADCEnergyEtaPhiCenterG1G2; //! patch ADC energy vs eta, phi at center of patch, low + high threshold
  TH3F     *fh3PatchADCEnergyEtaPhiCenterAll;  //! patch ADC energy vs eta, phi at center of patch, all trigger patches
  TH3F     *fh3EEtaPhiCell;                    //! cell E, eta, phi
  TH2F     *fh2ECellVsCent;                    //! cell E vs centrality
  TH2F     *fh2CellEnergyVsTime;               //! emcal cell energy vs time
  TH3F     *fh3EClusELeadingCellVsTime;        //! cluster energy vs energy of leading cell in cluster vs time of the leading cell
  TH3F     *fh3JetReacCent;                    //! jet energy vs cent vs dphi(jet,event plane)

  TH1F     *fHistClusEnergy;                   //!
  TH1F     *fHistClusofJetEnergy;              //!
  TH1F     *fHistEventSelectionQA;             //! trigger event class QA
  TH1F     *fhQAinfoAllPatchesCounter;         //!
  TH1F     *fhQAinfoCounter;                   //!
  TH1F     *fhQAmaxinfoCounter;                //!
  TH1F     *fhRecalcGammaPatchEnergy;          //! Recalculated Gamma Patch Energy distribution
  TH1F     *fhGammaLowPatchEnergy;             //! Gamma Low Patch Energy distribution
  TH1F     *fhGammaLowSimplePatchEnergy;       //! Gamma Low Simple Patch Energy distribution
  TH1F     *fhRecalcJetPatchEnergy;            //! Recalculated Jet Patch Energy distribution
  TH1F     *fhJetLowPatchEnergy;               //! Jet Low patch energy distribution
  TH1F     *fhJetLowSimplePatchEnergy;         //! Jet Low Simple patch energy distribution
  TH1F     *fhMainTriggerPatchEnergy;          //! Main Trigger patch energy distribution

  TH2F     *fHistdPhidEtaPatchJetCluster[16];  //!

  TClonesArray               *fJetTriggeredEvent;                    //!jets
  TClonesArray               *fRecalcTriggerPatches;                 //!recalculated patches

  THnSparse                  *fhnPatchMaxClus;//!                    // patch-maxclus distributions sparse matrix
  THnSparse                  *fhnPatchMatch;//!                      // QA before matching patch sparse matrix
  THnSparse                  *fhnPatchMatch2;//!                     // QA after matching patch sparse matrix
  THnSparse                  *fhnPatchMatchJetLeadClus;//!           // patch matching sparse matrix

  AliAnalysisTaskEmcalTriggerPatchJetMatch(const AliAnalysisTaskEmcalTriggerPatchJetMatch&);            // not implemented
  AliAnalysisTaskEmcalTriggerPatchJetMatch &operator=(const AliAnalysisTaskEmcalTriggerPatchJetMatch&); // not implemented

  ClassDef(AliAnalysisTaskEmcalTriggerPatchJetMatch, 4)
};
#endif

