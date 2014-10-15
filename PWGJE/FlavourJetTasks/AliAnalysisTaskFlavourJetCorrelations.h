#ifndef ALIANALYSISTASKFLAVOURJETCORRELATIONS_H
#define ALIANALYSISTASKFLAVOURJETCORRELATIONS_H
/**************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------
// Author : A. Grelli,  Utrecht University
//          C. Bianchin, Utrecht University
//          X. Zhang, LBNL
//-----------------------------------------------------------------------


#include <TH2F.h>
#include "AliAODEvent.h"
#include "AliPicoTrack.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliRhoParameter.h"

class TH3F;
class TParticle ;
class TClonesArray ;
class THnSparse;
class AliMCParticle;
class AliAODMCParticle;
class AliRDHFCuts;
class AliEmcalJet;
class AliAODRecoDecayHF;
class AliAODRecoCascadeHF;
class AliAODEvent;
class AliParticleContainer;
class AliClusterContainer;
class AliJetContainer;

class AliAnalysisTaskFlavourJetCorrelations : public AliAnalysisTaskEmcalJet 
{
   
public:
   
   enum ECandidateType{ kD0toKpi, kDstartoKpipi };
   
   AliAnalysisTaskFlavourJetCorrelations();
   AliAnalysisTaskFlavourJetCorrelations(const Char_t* name,AliRDHFCuts* cuts, ECandidateType candtype, Bool_t jetOnly=kFALSE);
   virtual ~AliAnalysisTaskFlavourJetCorrelations();
   
   virtual void     UserCreateOutputObjects();
   virtual Bool_t   Run();
   virtual void     Terminate(Option_t *);
   virtual void     Init();
   virtual void     LocalInit() {Init();}
   
   // inizializations
   Bool_t DefineHistoForAnalysis();  
   
   // set MC usage
   void   SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
   Bool_t GetMC() const {return fUseMCInfo;}
   // set usage of reconstructed tracks
   void   SetUseReco(Bool_t reco) {fUseReco=reco;}
   Bool_t GetUseReco() const {return fUseReco;}
   //no setter because needed in the constructor
   Bool_t GetJetOnlyMode() const {return fJetOnlyMode;}
   void SetTypeDJetSelection(Int_t type) {fTypeDInJet=type;} //see IsDInJet for possibilities
   Int_t GetTypeDJetSelection() const {return fTypeDInJet;} 
   
   void SetMassLimits(Double_t range, Int_t pdg);
   void SetMassLimits(Double_t lowlimit, Double_t uplimit);
   
   //jet and track arrays
   void SetJetArrayName(TString jetArrName) {fJetArrName=jetArrName;}
   TString GetJetArrayName() const {return fJetArrName;}
   void SetTrackArrayName(TString trackArrName) {fTrackArrName=trackArrName;}
   TString GetTrackArrayName() const {return fTrackArrName;}
   
   // trigger on jet events
   void SetTriggerOnLeadingJet(Bool_t triggerOnLeadingJet) {fLeadingJetOnly=triggerOnLeadingJet;};
   Bool_t GetTriggerOnLeadingJet() const {return fLeadingJetOnly;}
   
   // resolution parameter used in this task
   void SetRadius(Double_t r){fJetRadius=r;}
   Double_t GetRadius() const {return fJetRadius;}
   
   // Array of D0 width for the Dstar
   Bool_t SetD0WidthForDStar(Int_t nptbins,Float_t* width);
   
   //Bool_t   FillMCDJetInfo(AliPicoTrack *jetTrk,AliEmcalJet* jet, TClonesArray *mcArray,Double_t ptjet);
   void FillHistogramsRecoJetCorr(AliVParticle* candidate, AliEmcalJet *jet, AliAODEvent* aodEvent);
   void FillHistogramsD0JetCorr(AliAODRecoDecayHF* candidate, Double_t dPhi, Double_t z, Double_t ptD, Double_t ptj, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc, AliAODEvent* aodEvent);
   
   void FillHistogramsDstarJetCorr(AliAODRecoCascadeHF* dstar, Double_t dPhi,  Double_t z, Double_t ptD, Double_t ptj, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc);
   void FillHistogramsMCGenDJetCorr(Double_t dPhi, Double_t z,Double_t ptD,Double_t ptjet, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc);
   void SideBandBackground(AliAODRecoCascadeHF *candDstar, AliEmcalJet *jet);
   void MCBackground(AliAODRecoDecayHF *candbg, AliEmcalJet *jet);
   void FillMassHistograms(Double_t mass,Double_t ptD);
   void FlagFlavour(AliEmcalJet* jet);
   Int_t IsDzeroSideBand(AliAODRecoCascadeHF *candDstar);
   Bool_t InEMCalAcceptance(AliVParticle *vpart);
   
   void LightTHnSparse(){
      fSwitchOnSB=0;
      fSwitchOnPhiAxis=0;
      fSwitchOnOutOfConeAxis=0;
   }
   void HeavyTHnSparse(){
      fSwitchOnSB=1;
      fSwitchOnPhiAxis=1;
      fSwitchOnOutOfConeAxis=1;
   }
   
   void TurnOffTHnSparse() {fSwitchOnSparses=kFALSE;}
   
   /*
   void SwitchOnSB(){fSwitchOnSB=1;};
   void SwitchOffSB(){fSwitchOnSB=0;};
   void SwitchOnPhiAxis(){fSwitchOnPhiAxis=1;}
   void SwitchOffPhiAxis(){fSwitchOffPhiAxis=0;}
   void SwitchOnOutOfConeAxis(){fSwitchOnOutOfConeAxis=1;}
   void SwitchOffOutOfConeAxis(){fSwitchOnOutOfConeAxis=0;}
   */
   
private:
   
   AliAnalysisTaskFlavourJetCorrelations(const AliAnalysisTaskFlavourJetCorrelations &source);
   AliAnalysisTaskFlavourJetCorrelations& operator=(const AliAnalysisTaskFlavourJetCorrelations& source); 
   
   Double_t Z(AliVParticle* part,AliEmcalJet* jet, Bool_t transverse=kFALSE) const;
   Double_t Z(Double_t* p, Double_t *pj) const;
   Double_t ZT(Double_t* p, Double_t *pj) const;
   Float_t DeltaR(AliEmcalJet *p1, AliVParticle *p2) const;
   Bool_t AreDaughtersInJet(AliEmcalJet *thejet, AliAODRecoDecayHF* charm, Int_t*& daughOutOfJet, AliAODTrack**& dtrks, Bool_t fillH);
   Bool_t IsDInJet(AliEmcalJet *thejet, AliAODRecoDecayHF* charm, Bool_t fillH=kFALSE);
   void RecalculateMomentum(Double_t* pj, const Double_t* pmissing) const;
   
   
   Bool_t fUseMCInfo;             //  Use MC info
   Bool_t fUseReco;               // use reconstructed tracks when running on MC
   Int_t  fCandidateType;         // Dstar or D0
   Int_t  fPDGmother;             // PDG code of D meson
   Int_t  fNProngs;               // number of prong of the decay channel  
   Int_t  fPDGdaughters[4];       // PDG codes of daughters
   Float_t fSigmaD0[30];          //
   TString fBranchName;           // AOD branch name
   AliRDHFCuts *fCuts;            // Cuts 
   
   Double_t fMinMass;             // mass lower limit histogram
   Double_t fMaxMass;             // mass upper limit histogram
   
   TString  fJetArrName;          // name of the jet array, taken from the task running the jet finder
   TString fTrackArrName;         // name of the array of tracks, default "PicoTracks"
   TString fCandArrName;          // string which correspond to the candidate type
   Bool_t fLeadingJetOnly;        // use only the leading jet in the event to make the correlations
   Double_t fJetRadius;           // jet radius (filled from the JetContainer)
   TClonesArray *fCandidateArray;   //! contains candidates selected by AliRDHFCuts
   TClonesArray *fSideBandArray;    //! contains candidates selected by AliRDHFCuts::IsSelected(kTracks), to be used for side bands (DStar case only!!)
   Bool_t fJetOnlyMode;           //switch to simple version which analyses jets only
   Double_t fPmissing[3];             // jet missing momentum due to D mesons out of cone
   Double_t fPtJet;                   // pt jet, if requeted corrected with the missing pt
   Bool_t   fIsDInJet;                // flag D meson in jet
   Int_t    fTypeDInJet;              // way of selecting D in jets (see IsDInJet in cxx)
   TClonesArray *fTrackArr;      //array of the tracks used in the jet finder 
   Bool_t fSwitchOnSB;          // flag to switch off/on the SB analysis (default is off)
   Bool_t fSwitchOnPhiAxis;     // flag to switch off/on the DeltaPhi axis in THnSparse (to be used in combination with switch OutOfCone)
   Bool_t fSwitchOnOutOfConeAxis; //flag to switch off/on the out of cone axis in THnSparse (to be switch on if DeltaPhi is on)
   Bool_t fSwitchOnSparses;     // turn on/off all THnSparse
    
   Int_t fNAxesBigSparse;      // number of axis
   AliJetContainer      *fJetCont;      //! jets attachedto the task
   AliParticleContainer *fTrackCont;    //! tracks attached to the jet container
   AliClusterContainer  *fClusterCont;  //! Clusters  attached to the jet container 
  
   
   // Histograms
   TH1I* fhstat;                    //!
   //generic jet and jet track distributions
   TH1F* fhPtJetTrks;		    //!
   TH1F* fhPhiJetTrks;              //!
   TH1F* fhEtaJetTrks;              //!
   TH1F* fhEjetTrks;                //!
   TH1F* fhPtJet;                   //!
   TH1F* fhPhiJet;                  //!
   TH1F* fhEtaJet;                  //!
   TH1F* fhEjet;                    //!
   TH1F* fhNtrArr;                  //!
   TH1F* fhNJetPerEv;               //!
   TH1F* fhdeltaRJetTracks;         //!
   THnSparse* fhsJet; //available in jet only mode
   // event characteristics;        
   TH1F* fhNDPerEvNoJet;            //!
   TH1F* fhptDPerEvNoJet;           //!
   TH1F* fhNJetPerEvNoD;            //!
   TH1F* fhPtJetPerEvNoD;           //!
   //D mesons                       
   THnSparse* fhsDstandalone;       //!
   TH2F* fhInvMassptD;              //!
   TH2F* fhDiffSideBand;            //!
   TH2F* fhInvMassptDbg;            //!
   TH1F* fhPtPion;                  //!
                                    //!
   //histograms for checks          
   TH1F* fhztracksinjet;            //!
   TH1F* fhDiffPtTrPtJzNok;         //!
   TH1F* fhDiffPtTrPtJzok;          //!
   TH1F* fhDiffPzTrPtJzok;          //!
   TH1F* fhDiffPzTrPtJzNok;         //!
   TH1I* fhNtrkjzNok;               //!
   TH1F* fhztracksinjetT;           //!
   TH1I* fhControlDInJ;             //!
   TH1I* fhIDddaugh   ;             //!
   TH1I* fhIDddaughOut;             //!
   TH1I* fhIDjetTracks;             //!
   TH1F* fhDRdaughOut ;             //!
   TH1F* fhzDinjet;                 //!
   TH1F* fhmissingp;                //!
   TH1F**fhMissPi;                  //!
   TH1F* fhDeltaPtJet;              //!
   TH1F* fhRelDeltaPtJet;           //!
   // D-jet correlation histograms  
   TH1F* fhzDT;                     //!
   TH1F* fhDeltaRD;                 //!
   TH3F* fhDeltaRptDptj;            //!
   TH3F* fhDeltaRptDptjB;           //!
   //main histograms                
   THnSparse* fhsDphiz;             //!
   

   ClassDef(AliAnalysisTaskFlavourJetCorrelations,7); // class for charm-jet CorrelationsExch
};

#endif
