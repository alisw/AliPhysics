/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKRHOVNMODULATION_H
#define ALIANALYSISTASKRHOVNMODULATION_H

#include <AliAnalysisTaskEmcalJet.h>
#include <AliEmcalJet.h>
#include <AliVEvent.h>
#include <AliVTrack.h>
#include <AliVCluster.h>
#include <TClonesArray.h>
#include <TMath.h>

class TF1;
class THF1;
class THF2;
class TProfile;
class TRandom3;

class AliAnalysisTaskRhoVnModulation : public AliAnalysisTaskEmcalJet
{
    public:
         // enumerators
        enum fitModulationType  { kNoFit, kV2, kV3, kCombined, kUser, kFourierSeries }; // fit type
        enum runModeType        { kLocal, kGrid };                      // run mode type
        enum dataType           { kESD, kAOD, kESDMC, kAODMC };         // data type
        enum detectorType       { kTPC, kVZEROA, kVZEROC };             // detector that was used
        // constructors, destructor
                                AliAnalysisTaskRhoVnModulation();
                                AliAnalysisTaskRhoVnModulation(const char *name, runModeType type);
        virtual                 ~AliAnalysisTaskRhoVnModulation();
       
        // setting up the task and technical aspects
        Bool_t                  InitializeAnalysis();
        virtual void            UserCreateOutputObjects();
        virtual Bool_t          Run();
        TH1F*                   BookTH1F(const char* name, const char* x, Int_t bins, Double_t min, Double_t max, Int_t c = -1);
        TH2F*                   BookTH2F(const char* name, const char* x, const char* y, Int_t binsx, Double_t minx, Double_t maxx, Int_t binsy, Double_t miny, Double_t maxy, Int_t c = -1);
        /* inline */    Double_t PhaseShift(Double_t x) const {  
            while (x>=TMath::TwoPi())x-=TMath::TwoPi();
            while (x<0.)x+=TMath::TwoPi();
            return x; }
        /* inline */    Double_t ChiSquare(Int_t ndf, Double_t x) const {
            Double_t n(ndf/2.), denom(TMath::Power(2, n)*TMath::Gamma(n));
            if (denom!=0)  return ((1./denom)*TMath::Power(x, n-1)*TMath::Exp(-x/2.)); 
            return -999; }
        /* inline */    Double_t RhoVal() const { return (fRho) ? fRho->GetVal(): -999.;}                 
        /* inline */    Double_t RhoVal(Double_t phi, Double_t r, Double_t n) const {
            if(!fFitModulation) return -999; // coverity
            switch (fFitModulationType) {
                case kNoFit : return RhoVal();
                default : {
                    Double_t denom(2*r*fFitModulation->GetParameter(0));
                    return  (denom <= 0.) ? RhoVal() : n*(fFitModulation->Integral(phi-r, phi+r)/denom); 
                }
            }
        }
        // setters - analysis setup
        void                    SetDebugMode(Int_t d)                           {fDebug = d;}
        void                    SetFillQAHistograms(Bool_t qa)                  {fFillQAHistograms = qa;}
        void                    SetCentralityClasses(TArrayI* c)                {fCentralityClasses = c;}
        void                    SetNameJetClones(const char* name)              {fNameJetClones = name; }
        void                    SetNamePicoTrackClones(const char* name)        {fNamePicoTrackClones = name; }
        void                    SetNameRho(const char* name)                    {fNameRho = name; }
        void                    SetRandomSeed(TRandom3* r)                      {if (fRandom) delete fRandom; fRandom = r; }
        void                    SetModulationFit(TF1* fit)                      {if (fFitModulation) delete fFitModulation;
                                                                                 fFitModulation = fit; }
        void                    SetModulationFitMinimumPvalue(Float_t p)        {fMinPvalue = p;}
        void                    SetModulationFitType(fitModulationType type)    {fFitModulationType = type; }
        void                    SetModulationFitOptions(TString opt)            {fFitModulationOptions = opt; }
        void                    SetModulationFitDetector(detectorType type)     {fDetectorType = type; }
        void                    SetRunModeType(runModeType type)                {fRunModeType = type; }
        void                    SetAbsVertexZ(Float_t v)                        {fAbsVertexZ = v; }
        void                    SetMinDistanceRctoLJ(Float_t m)                 {fMinDisanceRCtoLJ = m; }
        void                    SetRandomConeRadius(Float_t r)                  {fRandomConeRadius = r; }
        // getters
        /* FIXME implement getters */
        // 'trivial' helper calculations
        void                    CalculateEventPlaneVZERO(Double_t vzero[2][2]) const;
        void                    CalculateEventPlaneTPC(Double_t* tpc) const;
        void                    CalculateRandomCone(Float_t &pt, Float_t &eta, Float_t &phi, AliEmcalJet* jet = 0x0, Bool_t randomize = 0) const;
        // analysis details
        void                    CorrectRho(Double_t* params, Double_t psi2, Double_t psi3) const;
        // event and track selection
        /* inline */    Bool_t PassesCuts(const AliVTrack* track) const {
            if(!track) return kFALSE;
            return kTRUE; }
        Bool_t                  PassesCuts(AliVEvent* event);
        Bool_t                  PassesCuts(const AliVCluster* track) const;
        Bool_t                  PassesCuts(const AliEmcalJet* jet) const;
        // filling histograms
        void                    FillHistogramsAfterSubtraction(Double_t vzero[2][2], Double_t* tpc) const;
        void                    FillTrackHistograms() const;
        void                    FillClusterHistograms() const;
        void                    FillCorrectedClusterHistograms() const;
        void                    FillEventPlaneHistograms(Double_t vzero[2][2], Double_t* tpc) const;
        void                    FillRhoHistograms() const;
        void                    FillDeltaPtHistograms(Double_t* tpc) const; 
        void                    FillJetHistograms(Double_t vzero[2][2], Double_t* tpc) const;
        void                    FillDeltaPhiHistograms(Double_t vzero[2][2], Double_t* tpc) const;
        void                    FillQAHistograms(AliVTrack* vtrack) const;
        void                    FillQAHistograms(AliVEvent* vevent);
        virtual void            Terminate(Option_t* option);
    private:
        // analysis flags and settings
        Int_t                   fDebug;                 // debug level (0 none, 1 fcn calls, 2 verbose)
        Bool_t                  fInitialized;           //! is the analysis initialized?
        Bool_t                  fFillQAHistograms;      // fill qa histograms
        TArrayI*                fCentralityClasses;     //-> centrality classes (maximum 10)
        // members
        fitModulationType       fFitModulationType;     // fit modulation type
        detectorType            fDetectorType;          // type of detector used for modulation fit
        TString                 fFitModulationOptions;  // fit options for modulation fit
        runModeType             fRunModeType;           // run mode type 
        dataType                fDataType;              // datatype 
        TRandom3*               fRandom;                //-> dont use gRandom to not interfere with other tasks
        Int_t                   fMappedRunNumber;       //! mapped runnumer (for QA)
        Int_t                   fInCentralitySelection; //! centrality bin
        TF1*                    fFitModulation;         //-> modulation fit for rho
        Float_t                 fMinPvalue;             // minimum value of p
        const char*             fNameJetClones;         //! collection of tclones array with jets
        const char*             fNamePicoTrackClones;   //! collection of tclones with pico tracks
        const char*             fNameRho;               //! name of rho
        // event cuts
        Float_t                 fAbsVertexZ;            // cut on zvertex
        // general qa histograms
        TH1F*                   fHistCentrality;        //! accepted centrality
        TH1F*                   fHistVertexz;           //! accepted verte
        TH2F*                   fHistRunnumbersPhi;     //! run numbers averaged phi
        TH2F*                   fHistRunnumbersEta;     //! run numbers averaged eta
        // general settings
        Float_t                 fMinDisanceRCtoLJ;      //! min distance between rc and leading jet
        Float_t                 fRandomConeRadius;      //! radius of random cone
        // transient object pointers
        TList*                  fOutputList;            //! output list
        TList*                  fOutputListGood;        //! output list for local analysis
        TList*                  fOutputListBad;         //! output list for local analysis
        TH1F*                   fHistAnalysisSummary;   //! analysis summary
        TH1F*                   fHistSwap;              //! swap histogram
        TProfile*               fProfVn;                //! extracted vn
        // qa histograms for accepted pico tracks
        TH1F*                   fHistPicoTrackPt[10];    //! pt of all charged tracks
        TH2F*                   fHistPicoCat1[10];       //! pico tracks spd hit and refit
        TH2F*                   fHistPicoCat2[10];       //! pico tracks wo spd hit w refit, constrained
        TH2F*                   fHistPicoCat3[10];       //! pico tracks wo spd hit wo refit, constrained
        // qa histograms for accepted emcal clusters
        /* TH1F*                   fHistClusterPt[10];      //! pt uncorrected emcal clusters */
        /* TH1F*                   fHistClusterPhi[10];     //! phi uncorrected emcal clusters */
        /* TH1F*                   fHistClusterEta[10];     //! eta uncorrected emcal clusters */
        //// qa histograms for accepted emcal clusters aftehadronic correction
        /* TH1F*                   fHistClusterCorrPt[10];  //! pt corrected emcal clusters */
        /* TH1F*                   fHistClusterCorrPhi[10]; //! phi corrected emcal clusters */
        /* TH1F*                   fHistClusterCorrEta[10]; //! eta corrected emcal clusters */
        // qa event planes
        TProfile*               fHistPsi2;               //! Psi_2 estimates
        TProfile*               fHistPsi2Spread;         //! spread between 
        TH1F*                   fHistPsiVZEROA;          //! psi 2 from vzero a
        TH1F*                   fHistPsiVZEROC;          //! psi 2 from vzero c
        TH1F*                   fHistPsiTPC;             //! psi 2 from tpc
        // background
        TH1F*                   fHistRhoPackage[10];     //! rho as estimated by emcal jet package
        TH1F*                   fHistRho[10];            //! background
        TH2F*                   fHistRhoVsMult;          //! rho versus multiplicity
        TH2F*                   fHistRhoVsCent;          //! rho veruss centrality
        TH2F*                   fHistRhoAVsMult;         //! rho * A vs multiplicity for all jets
        TH2F*                   fHistRhoAVsCent;         //! rho * A vs centrality for all jets
        // delta pt distributions
        TH2F*                   fHistRCPhiEta[10];              //! random cone eta and phi
        TH2F*                   fHistRhoVsRCPt[10];             //! rho * A vs rcpt
        TH1F*                   fHistRCPt[10];                  //! rcpt
        TH2F*                   fHistDeltaPtDeltaPhi[10];       //! dpt vs dphi
        TH2F*                   fHistRCPhiEtaExLJ[10];          //! random cone eta and phi, excl leading jet
        TH2F*                   fHistRhoVsRCPtExLJ[10];         //! rho * A vs rcpt, excl leading jet
        TH1F*                   fHistRCPtExLJ[10];              //! rcpt, excl leading jet
        TH2F*                   fHistDeltaPtDeltaPhiExLJ[10];   //! dpt vs dphi, excl leading jet
        TH2F*                   fHistRCPhiEtaRand[10];          //! random cone eta and phi, randomized
        TH2F*                   fHistRhoVsRCPtRand[10];         //! rho * A vs rcpt, randomized
        TH1F*                   fHistRCPtRand[10];              //! rcpt, randomized
        TH2F*                   fHistDeltaPtDeltaPhiRand[10];   //! dpt vs dphi, randomized
        // jet histograms (after kinematic cuts)
        TH1F*                   fHistJetPtRaw[10];              //! jet pt - no background subtraction
        TH1F*                   fHistJetPt[10];                 //! pt of found jets (background subtracted)
        TH2F*                   fHistJetEtaPhi[10];             //! eta and phi correlation
        TH2F*                   fHistJetPtArea[10];             //! jet pt versus area
        TH2F*                   fHistJetPtConstituents[10];     //! jet pt versus number of constituents
        // in plane, out of plane jet spectra
        TH2F*                   fHistJetPsiTPCPt[10];            //! psi tpc versus pt
        TH2F*                   fHistJetPsiVZEROAPt[10];         //! psi vzeroa versus pt
        TH2F*                   fHistJetPsiVZEROCPt[10];         //! psi vzeroc versus pt
        // phi minus psi 
        TH1F*                   fHistDeltaPhiVZEROA[10]; //! phi minus psi_A
        TH1F*                   fHistDeltaPhiVZEROC[10]; //! phi minus psi_C
        TH1F*                   fHistDeltaPhiTPC[10];    //! phi minus psi_TPC

        AliAnalysisTaskRhoVnModulation(const AliAnalysisTaskRhoVnModulation&);                  // not implemented
        AliAnalysisTaskRhoVnModulation& operator=(const AliAnalysisTaskRhoVnModulation&);       // not implemented

        ClassDef(AliAnalysisTaskRhoVnModulation, 2);
};

#endif
