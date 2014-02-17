#ifndef ALIANALYSISTASKEMCALJETTRIGGERQA_H
#define ALIANALYSISTASKEMCALJETTRIGGERQA_H

class TH1;
class TH2;
class TH3;
class TH3F;
class TClonesArray;
class TArrayI;

#include <TRef.h>
#include <TBits.h>
#include <TMath.h>

#include <AliVEvent.h>

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalJetTriggerQA : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskEmcalJetTriggerQA();
  AliAnalysisTaskEmcalJetTriggerQA(const char *name);
  virtual ~AliAnalysisTaskEmcalJetTriggerQA();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  Bool_t                      SelectEvent();              //decides if event is used for analysis
  void                        FindTriggerPatch();

  //Setters
  void SetDebug(Int_t d)                    { fDebug = d;}
  void SetTriggerClass(const char *n)       { fTriggerClass = n; }
  void SetNFastorPatch(Int_t i)             { fNFastOR = i;}
 
  void SetContainerFull(Int_t c)            { fContainerFull      = c;}
  void SetContainerCharged(Int_t c)         { fContainerCharged   = c;}

  Int_t    GetLeadingCellId(const AliVCluster *clus) const;
  Double_t GetEnergyLeadingCell(const AliVCluster *clus) const;
  Double_t GetECross(Int_t absID) const;

  Double_t GetZ(const AliVParticle *trk, const AliEmcalJet *jet)       const;
  Double_t GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz, const Double_t jetPx, const Double_t jetPy, const Double_t jetPz) const;

 protected:
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;
  Float_t                     RelativeEP(Double_t objAng, Double_t EPAng) const;
    
  Bool_t                      TestFilterBit(Int_t trigBit, UInt_t bitJetTrig) const {return (Bool_t) ((trigBit & bitJetTrig) != 0);}


 private:
  Bool_t             fDebug;                 //  debug level
  TString            fTriggerClass;          // trigger class to analyze EJ1 or EJ2    
  UInt_t             fBitJ1;                 // trigger bit of EJE1
  UInt_t             fBitJ2;                 // trigger bit of EJE2
  Int_t              fContainerFull;         //  number of container with full jets DET
  Int_t              fContainerCharged;      //  number of container with charged jets DET
  Double_t           fMaxPatchEnergy;        // energy of patch with largest energy
  Int_t              fTriggerType;           // trigger type
  Int_t              fNFastOR;               // size of trigger patch fNFastORxfNFastOR

  TH1F  *fhNEvents;                         //! Histo number of events
  TH2F  *fHistRhovsCentFull;                //!rho vs. centrality
  TH2F  *fHistRhovsCentCharged;             //!rho vs. centrality
  TH3F  *fh3PtEtaPhiTracks;                 //! pt,eta,phi of tracks
  TH3F  *fh3PtEtaPhiTracksOnEmcal;          //! pt,eta,phi of tracks
  TH3F  *fh3PtEtaPhiTracksProp;             //! pt,eta,phi of tracks
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
  TH3F  *fh3NEFEtaPhi;                      //! NEF, eta, phi
  TH2F  *fh2NEFNConstituentsCharged;        //! NEF, # charged jet constituents
  TH2F  *fh2NEFNConstituentsNeutral;        //! NEF, # neutral jet constituents
  TH2F  *fh2Ptz;                            //! pt, z=pT,h,proj/p,jet full jet
  TH2F  *fh2PtzCharged;                     //! pt, z=pT,h,proj/p,jet charged jet
  TH2F  *fh2PtLeadJet1VsLeadJet2;           //! correlation between leading jet of the two branches
  TH3F  *fh3EEtaPhiCluster;                 //! cluster E, eta, phi
  TH3F  *fh3PtLeadJet1VsPatchEnergy;        //! leading jet energy vs leading patch energy vs jet trigger (J1/J2)
  TH3F  *fh3PtLeadJet2VsPatchEnergy;        //! leading jet energy vs leading patch energy vs jet trigger (J1/J2)
  TH3F  *fh3PatchEnergyEtaPhiCenterJ1;      //! patch energy vs eta, phi at center of patch, high threshold
  TH3F  *fh3PatchEnergyEtaPhiCenterJ2;      //! patch energy vs eta, phi at center of patch, low threshold
  TH3F  *fh3PatchEnergyEtaPhiCenterJ1J2;    //! patch energy vs eta, phi at center of patch, low + high threshold
  TH3F  *fh3PatchADCEnergyEtaPhiCenterJ1;   //! patch ADC energy vs eta, phi at center of patch, high threshold
  TH3F  *fh3PatchADCEnergyEtaPhiCenterJ2;   //! patch ADC energy vs eta, phi at center of patch, low threshold
  TH3F  *fh3PatchADCEnergyEtaPhiCenterJ1J2; //! patch ADC energy vs eta, phi at center of patch, low + high threshold
  TH2F  *fh2CellEnergyVsTime;               //! emcal cell energy vs time
  TH3F  *fh3EClusELeadingCellVsTime;        //! cluster energy vs energy of leading cell in cluster vs time of the leading cell

  TH3F  *fh3JetReacCent;
  TH2F  *fh2FullJetCent;
    
  AliAnalysisTaskEmcalJetTriggerQA(const AliAnalysisTaskEmcalJetTriggerQA&);            // not implemented
  AliAnalysisTaskEmcalJetTriggerQA &operator=(const AliAnalysisTaskEmcalJetTriggerQA&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetTriggerQA, 9)
};
#endif
