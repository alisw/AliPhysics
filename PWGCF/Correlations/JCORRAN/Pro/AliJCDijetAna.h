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

#ifndef ALIJCDIJETANA_H
#define ALIJCDIJETANA_H

#include <TH1D.h>
#include <TClonesArray.h>
#include "AliJBaseTrack.h"
#include "AliJCDijetHistos.h"
#include "TRandom3.h"

// Fastjet includes
#include "FJ_includes.h"

class AliJCDijetHistos;

class AliJCDijetAna : public TObject
{
    public:
        AliJCDijetAna(); // Default contructor
        virtual ~AliJCDijetAna(); // Destructor
        AliJCDijetAna(const AliJCDijetAna& obj); // Copy constructor
        AliJCDijetAna& operator=(const AliJCDijetAna& obj); // Equal sign operator

#if !defined(__CINT__) && !defined(__MAKECINT__)
        vector<vector<fastjet::PseudoJet>> GetJets() { return jets; }
        vector<vector<vector<fastjet::PseudoJet>>> GetDijets() { return dijets; }
        bool HasDijet() { return bHasDijet; }
        bool HasDeltaPhiDijet() { return bHasDeltaPhiDijet; }
        void InitHistos(AliJCDijetHistos *histos, bool bIsMC, int nCentBins);

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

        int CalculateJets(TClonesArray *inList, AliJCDijetHistos *fhistos, int lCBin);
        void SetJets(vector<fastjet::PseudoJet> jetsOutside);
        void SetPtHardBin(double flptHardBin) {fptHardBin = flptHardBin; }
        void FillJetsDijets(AliJCDijetHistos *fhistos, int lCBin);
        void CalculateResponse(AliJCDijetAna *anaDetMC, AliJCDijetHistos *fhistos);
        void ResetObjects();
        double DeltaR(fastjet::PseudoJet jet1, fastjet::PseudoJet jet2);
        double DeltaR(double eta1, double eta2, double phi1, double phi2);
        bool CheckDeltaPhi(fastjet::PseudoJet leadingJet, fastjet::PseudoJet subleadingJet, double deltaPhiCut);
        double GetDeltaPhi(fastjet::PseudoJet leadingJet, fastjet::PseudoJet subleadingJet);
#endif

    private:
        int fDebug;
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

        enum jetClasses {iAcc, iBGSubtr, iBGSubtrConstCut, iConstCut, iktJets, jetClassesSize};
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

        ClassDef(AliJCDijetAna, 1); // ClassDef needed if inheriting from TObject

};

#endif
