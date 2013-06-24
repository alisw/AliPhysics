#ifndef ALIANALYSISTASKEMCALJETTRIGGERQA_H
#define ALIANALYSISTASKEMCALJETTRIGGERQA_H

class TH1;
class TH2;
class TH3;
class TH3F;
class TClonesArray;
class TArrayI;
class AliAnalysisUtils;

#include <TRef.h>
#include <TBits.h>

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalJetTriggerQA : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskEmcalJetTriggerQA();
  AliAnalysisTaskEmcalJetTriggerQA(const char *name);
  virtual ~AliAnalysisTaskEmcalJetTriggerQA();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  void                        InitOnce();
  void                        LoadExtraBranches();

  Bool_t                      SelectEvent();              //decides if event is used for analysis
  void                        FindTriggerPatch();

  Bool_t                      AcceptJet2(const AliEmcalJet *jet) const;

  //Setters
  void SetDebug(Int_t d)                    { fDebug = d;}
  void SetTriggerClass(const char *n)       { fTriggerClass = n; }
  void SetNFastorPatch(Int_t i)             { fNFastOR = i;}
 
  void SetJetsName2(const char *n)          { fJetsName2 = n; }
  void SetRhoChName(const char *n)          { fRhoChName = n; }
  void SetMinEtaJets2(Double_t p)           { fEtaMinJet2 = p;}
  void SetMaxEtaJets2(Double_t p)           { fEtaMaxJet2 = p;}
  void SetMinPhiJets2(Double_t p)           { fPhiMinJet2 = p;}
  void SetMaxPhiJets2(Double_t p)           { fPhiMaxJet2 = p;}

  Double_t GetZ(const AliVParticle *trk, const AliEmcalJet *jet)       const;
  Double_t GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz, const Double_t jetPx, const Double_t jetPy, const Double_t jetPz) const;

 protected:
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;

  Bool_t                      TestFilterBit(Int_t trigBit, UInt_t bitJetTrig) const {return (Bool_t) ((trigBit & bitJetTrig) != 0);}



 private:
  Bool_t            fDebug;                 //  debug level
  Bool_t            fUseAnaUtils;           //  used for LHC13* data
  AliAnalysisUtils *fAnalysisUtils;         //! vertex selection
  TString           fJetsName2;             //  name of charged jet collection
  TClonesArray     *fJets2;                 //! list with charged jets
  Double_t          fEtaMinJet2;            //  min eta of jet axis, jets in 2nd branch
  Double_t          fEtaMaxJet2;            //  max eta of jet axis, jets in 2nd branch
  Double_t          fPhiMinJet2;            //  min phi of jet axis, jets in 2nd branch
  Double_t          fPhiMaxJet2;            //  max phi of jet axis, jets in 2nd branch
  Double_t          fMaxTrackPtJet2;        //  maximum track pT in jet

  TString           fRhoChName;             //  name of charged rho branch
  AliRhoParameter  *fRhoCh;                 //! event rho charged
  Double_t          fRhoChVal;              //  charged rho value

  TString           fTriggerClass;          // trigger class to analyze EJ1 or EJ2    
  UInt_t            fBitJ1;                 // trigger bit of EJE1
  UInt_t            fBitJ2;                 // trigger bit of EJE2
  Double_t          fMaxPatchEnergy;        // energy of patch with largest energy
  Int_t             fTriggerType;           // trigger type
  Int_t             fNFastOR;               // size of trigger patch fNFastORxfNFastOR

  TH1F  *fhNEvents;                         //! Histo number of events
  TH3F  *fh3PtEtaPhiJetFull;                //! pt,eta,phi of full jets
  TH3F  *fh3PtEtaPhiJetCharged;             //! pt,eta,phi of charged jets
  TH2F  *fh2NJetsPtFull;                    //! NJets per event vs pT,jet
  TH2F  *fh2NJetsPtCharged;                 //! NJets per event vs pT,jet
  TH3F  *fh3PtEtaAreaJetFull;               //! pt,eta,area of full jet
  TH3F  *fh3PtEtaAreaJetCharged;            //! pt,eta,area of charged jets
  TH2F  *fh2PtNConstituentsCharged;         //! pt, # charged jet constituents
  TH2F  *fh2PtNConstituents;                //! pt, # jet constituents
  TH2F  *fh2PtMeanPtConstituentsCharged;    //! pt, <pt> charged constituents
  TH2F  *fh2PtMeanPtConstituentsNeutral;    //! pt, <pt> neutral constituents
  TH2F  *fh2PtNEF;                          //! pt, NEF (neutral energy fraction)
  TH2F  *fh2Ptz;                            //! pt, z=pT,h,proj/p,jet
  TH2F  *fh2PtLeadJet1VsLeadJet2;           //! correlation between leading jet of the two branches
  TH3F  *fh3EEtaPhiCluster;                 //! cluster E, eta, phi
  TH3F  *fh3PtLeadJet1VsPatchEnergy;        //! leading jet energy vs leading patch energy vs jet trigger (J1/J2)
  TH3F  *fh3PtLeadJet2VsPatchEnergy;        //! leading jet energy vs leading patch energy vs jet trigger (J1/J2)
  TH3F  *fh3PatchEnergyEtaPhiCenter;        //! patch energy vs eta, phi at center of patch
  TH2F  *fh2CellEnergyVsTime;               //! emcal cell energy vs time


  AliAnalysisTaskEmcalJetTriggerQA(const AliAnalysisTaskEmcalJetTriggerQA&);            // not implemented
  AliAnalysisTaskEmcalJetTriggerQA &operator=(const AliAnalysisTaskEmcalJetTriggerQA&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetTriggerQA, 1) // jet sample analysis task
};
#endif
