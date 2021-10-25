/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
//______________________________________________________________________________
// Analysis task for providing various dijet informations
// author: O. Saarimaki, D.J. Kim(dong.jo.kim@cern.ch)
// ALICE Group University of Jyvaskyla
// Finland
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////

// used in local and grid execution

#ifndef ALIJEMCALDIJETANA_H
#define ALIJEMCALDIJETANA_H

#include <TH1D.h>
#include <TClonesArray.h>
#include "AliJBaseTrack.h"
#include "TRandom3.h"
#include "THistManager.h"
#include "AliAnalysisTaskEmcalJet.h"

// Fastjet includes
#include "FJ_includes.h"

class AliJEmcalDijetAna : public AliAnalysisTaskSE
{
    public:
        AliJEmcalDijetAna(); // Default contructor
        virtual ~AliJEmcalDijetAna(); // Destructor
        AliJEmcalDijetAna(const AliJEmcalDijetAna& obj); // Copy constructor
        AliJEmcalDijetAna& operator=(const AliJEmcalDijetAna& obj); // Equal sign operator

#if !defined(__CINT__) && !defined(__MAKECINT__)
        vector<vector<fastjet::PseudoJet>> GetJets() { return jets; }
        vector<vector<vector<fastjet::PseudoJet>>> GetDijets() { return dijets; }
        bool HasDijet() { return bHasDijet; }
        bool HasDeltaPhiDijet() { return bHasDeltaPhiDijet; }
        void InitHistos(bool bIsMC);

        void SetSettings(int    lDebug,
                         double lParticleEtaCut,
                         double lParticlePtCut,
                         double lJetCone,
                         double lktJetCone,
                         int    lktScheme,
                         int    lantiktScheme,
                         bool   lusePionMass,
                         bool   luseDeltaPhiBGSubtr,
                         double lConstituentCut,
                         double lLeadingJetCut,
                         double lSubleadingJetCut,
                         double lMinJetPt,
                         double lDeltaPhiCut,
                         double lmatchingR,
                         double ltrackingIneff);

        int CalculateJets(TClonesArray *inList, int lCBin);
        void SetJets(vector<fastjet::PseudoJet> jetsOutside);
        void SetPtHardBin(double flptHardBin) {fptHardBin = flptHardBin; }
        void SetParticleCollArray(TObjArray arr) {fParticleCollArray = arr; }
        void SetCentBins(int ncents) {fNcentBins = ncents; }
        void FillJetsDijets(int lCBin);
        void CalculateResponse(AliJEmcalDijetAna *anaDetMC, int iJetSetPart, int iJetSetDet);
        void ResetObjects();
        double DeltaR(fastjet::PseudoJet jet1, fastjet::PseudoJet jet2);
        double DeltaR(double eta1, double eta2, double phi1, double phi2);
        bool CheckDeltaPhi(fastjet::PseudoJet leadingJet, fastjet::PseudoJet subleadingJet, double deltaPhiCut);
        double GetDeltaPhi(fastjet::PseudoJet leadingJet, fastjet::PseudoJet subleadingJet);
        //void SetBinning(TH1F *h) {hBinning=h;}
        //void SetBinningALICE(TH1F *h) {hBinningALICE=h;}
        //void SetHistGroupName(TString groupname) { fGroupname = groupname; }

        enum jetClasses {iAcc, iBGSubtr, iBGSubtrConstCut, iConstCut, iktJets, jetClassesSize};
#endif

    protected:
        void AllocateTrackHistograms();
        THistManager fhistos ;///< Histogram manager

    private:
        int fDebug;
        int fNcentBins;
        double fParticleEtaCut;
        double fParticlePtCut;
        bool fusePionMass;
        bool fUseDeltaPhiBGSubtr;
        double fConstituentCut;
        double fLeadingJetCut;
        double fSubleadingJetCut;
        double fDeltaPhiCut;
        double etaMaxCutForJet;
        double etaMaxCutForKtJet;
        double MinJetPt;
        double fJetCone;
        double fktJetCone;
        double pionmass;
        double matchingR;
        double ftrackingIneff;
        bool bEvtHasAreaInfo;

        double phi, eta, pt, pt2, rho, rhom, area, mjj, ptpair, dPhi, deltaRMin, deltaR;
        bool leadingTrackOverThreshold;
        unsigned noTracks;
        bool removed;
        //For loops:
        unsigned utrack, uktjet, ujet, ujet2, uconst, udijet, ujetDetMC;
        bool bHasDijet;
        bool bHasDeltaPhiDijet;
        bool bHasDeltaPhiSubLeadJet;
        TRandom3 *randomGenerator;
        double fDeltaPt;
        double fptHardBin;
        double randConePhi;
        double randConeEta;
        double randConePt;
        double fDeltaM;
        //TH1F *hBinning;
        //TH1F *hBinningALICE;
        TString fGroupname;
        TString fHistname;
        TObjArray fParticleCollArray;
        AliEmcalList *fOutput;                     //!<!output list

#if !defined(__CINT__) && !defined(__MAKECINT__)
        vector<fastjet::PseudoJet> chparticles;
        vector<fastjet::PseudoJet> ktchparticles;
        // This list contains all accepted jets in different categories:
        vector<vector<fastjet::PseudoJet>> jets;
        // These 'raw' lists contain all jets by fastjet:
        vector<fastjet::PseudoJet> rawJets;
        vector<fastjet::PseudoJet> rawKtJets;
        vector<fastjet::PseudoJet> rhoEstJets;
        vector<vector<vector<fastjet::PseudoJet>>> dijets;

        fastjet::RecombinationScheme ktScheme;
        fastjet::RecombinationScheme antiktScheme;
        fastjet::PseudoJet jetAreaVector;
        fastjet::PseudoJet jet_bgSubtracted;
        fastjet::PseudoJet dijet;

        fastjet::JetDefinition jet_def;
        fastjet::JetDefinition jet_def_bge;

        fastjet::GhostedAreaSpec area_spec;
        fastjet::AreaDefinition area_def;
        fastjet::AreaDefinition area_def_bge;

        fastjet::Selector selectorAllButTwo;
        fastjet::Selector selectorEta;
        fastjet::Selector selectorBoth;
        fastjet::JetMedianBackgroundEstimator bge;

        unique_ptr<fastjet::ClusterSequenceArea> cs;
        unique_ptr<fastjet::ClusterSequenceArea> cs_bge;
#endif

        ClassDef(AliJEmcalDijetAna, 1); // ClassDef needed if inheriting from TObject

};

#endif
