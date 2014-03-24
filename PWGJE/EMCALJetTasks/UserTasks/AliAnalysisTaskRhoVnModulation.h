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
#include <TArrayD.h>
#include <TRandom3.h>
#include <AliJetContainer.h>
#include <AliParticleContainer.h>

class TF1;
class THF1;
class THF2;
class TProfile;
class AliLocalRhoParameter;

class AliAnalysisTaskRhoVnModulation : public AliAnalysisTaskEmcalJet {
    public:
         // enumerators
        enum fitModulationType  { kNoFit, kV2, kV3, kCombined, kFourierSeries, kIntegratedFlow, kQC2, kQC4 }; // fit type
        enum fitGoodnessTest    { kChi2ROOT, kChi2Poisson, kKolmogorov, kKolmogorovTOY, kLinearFit };
        enum collisionType      { kPbPb, kPythia };                     // collision type
        enum qcRecovery         { kFixedRho, kNegativeVn, kTryFit };    // how to deal with negative cn value for qcn value
        enum runModeType        { kLocal, kGrid };                      // run mode type
        enum dataType           { kESD, kAOD, kESDMC, kAODMC };         // data type
        enum detectorType       { kTPC, kVZEROA, kVZEROC, kVZEROComb};  // detector that was used
        // constructors, destructor
                                AliAnalysisTaskRhoVnModulation();
                                AliAnalysisTaskRhoVnModulation(const char *name, runModeType type);
        virtual                 ~AliAnalysisTaskRhoVnModulation();
        // setting up the task and technical aspects
        void                    ExecOnce();
        Bool_t                  InitializeAnalysis();
        virtual void            UserCreateOutputObjects();
        virtual Bool_t          Run();
        TH1F*                   BookTH1F(const char* name, const char* x, Int_t bins, Double_t min, Double_t max, Int_t c = -1, Bool_t append = kTRUE);
        TH2F*                   BookTH2F(const char* name, const char* x, const char* y, Int_t binsx, Double_t minx, Double_t maxx, Int_t binsy, Double_t miny, Double_t maxy, Int_t c = -1, Bool_t append = kTRUE);
        /* inline */    Double_t PhaseShift(Double_t x) const {  
            while (x>=TMath::TwoPi())x-=TMath::TwoPi();
            while (x<0.)x+=TMath::TwoPi();
            return x; }
        /* inline */    Double_t PhaseShift(Double_t x, Double_t n) const {
            x = PhaseShift(x);
            if(TMath::Nint(n)==2) while (x>TMath::Pi()) x-=TMath::Pi();
            if(TMath::Nint(n)==3) {
                if(x>2.*TMath::TwoPi()/n) x = TMath::TwoPi() - x;
                if(x>TMath::TwoPi()/n) x = TMath::TwoPi()-(x+TMath::TwoPi()/n);
            }
            return x; }
        /* inline */    Double_t ChiSquarePDF(Int_t ndf, Double_t x) const {
            Double_t n(ndf/2.), denom(TMath::Power(2, n)*TMath::Gamma(n));
            if (denom!=0)  return ((1./denom)*TMath::Power(x, n-1)*TMath::Exp(-x/2.)); 
            return -999; }
        // note that the cdf of the chisquare distribution is the normalized lower incomplete gamma function
        /* inline */    Double_t ChiSquareCDF(Int_t ndf, Double_t x) const { return TMath::Gamma(ndf/2., x/2.); }
        /* inline */    Double_t ChiSquare(TH1& histo, TF1* func) const {
            // evaluate the chi2 using a poissonian error estimate on bins
            Double_t chi2(0.);
            for(Int_t i(0); i < histo.GetXaxis()->GetNbins(); i++) {
                if(histo.GetBinContent(i+1) <= 0.) continue;
                chi2 += TMath::Power((histo.GetBinContent(i+1)-func->Eval(histo.GetXaxis()->GetBinCenter(1+i))), 2)/histo.GetBinContent(i+1);
            }
           return chi2;
        }
        /* inline*/ Double_t KolmogorovTest(TH1F& histo, TF1* func) const {
            // return the probability from a Kolmogorov test
            TH1F test(histo);       // stack copy of test statistic
            for(Int_t i(0); i < test.GetXaxis()->GetNbins(); i++) test.SetBinContent(i+1, func->Eval(test.GetXaxis()->GetBinCenter(1+i)));
            if(fFitGoodnessTest == kKolmogorovTOY) return histo.TH1::KolmogorovTest((&test), "X");
            return histo.TH1::KolmogorovTest((&test));
        }
 
        // setters - analysis setup
        void                    SetDebugMode(Int_t d)                           {fDebug = d;}
        void                    SetRunToyMC(Bool_t t)                           {fRunToyMC = t; }
        void                    SetAttachToEvent(Bool_t b)                      {fAttachToEvent = b;}
        void                    SetSemiCentralInclusive(Bool_t b)               {fSemiCentralInclusive = b;}
        void                    SetFillHistograms(Bool_t b)                     {fFillHistograms = b;}
        void                    SetFillQAHistograms(Bool_t qa)                  {fFillQAHistograms = qa;}
        void                    SetReduceBinsXYByFactor(Float_t x, Float_t y)   {fReduceBinsXByFactor = x;
                                                                                 fReduceBinsYByFactor = y;}
        void                    SetNoEventWeightsForQC(Bool_t e)                {fNoEventWeightsForQC = e;}
        void                    SetCentralityClasses(TArrayD* c)                {fCentralityClasses = c;}
        void                    SetPtBinsHybrids(TArrayD* p)                    {fPtBinsHybrids = p;}
        void                    SetPtBinsJets(TArrayD* p)                       {fPtBinsJets = p;}
        void                    SetExpectedRuns(TArrayI* r)                     {fExpectedRuns = r;}
        void                    SetExpectedSemiGoodRuns(TArrayI* r)             {fExpectedSemiGoodRuns = r;}
        void                    SetIntegratedFlow(TH1F* i, TH1F* j)             {fUserSuppliedV2 = i;
                                                                                 fUserSuppliedV3 = j; }
        void                    SetOnTheFlyResCorrection(TH1F* r2, TH1F* r3)    {fUserSuppliedR2 = r2;
                                                                                 fUserSuppliedR3 = r3; }
        void                    SetNameJetClones(const char* name)              {fNameJetClones = name; }
        void                    SetNamePicoTrackClones(const char* name)        {fNamePicoTrackClones = name; }
        void                    SetNameRho(const char* name)                    {fNameRho = name; }
        void                    SetNameRhoSmall(TString name)                   {fNameSmallRho = name; }
        void                    SetUseScaledRho(Bool_t s)                       {fUseScaledRho = s; }
        void                    SetRandomSeed(TRandom3* r)                      {if (fRandom) delete fRandom; fRandom = r; }
        void                    SetModulationFit(TF1* fit);
        void                    SetUseControlFit(Bool_t c);
        void                    SetModulationFitMinMaxP(Float_t m, Float_t n)   {fMinPvalue = m; fMaxPvalue = n; }
        void                    SetModulationFitType(fitModulationType type)    {fFitModulationType = type; }
        void                    SetGoodnessTest(fitGoodnessTest test)           {fFitGoodnessTest = test; }
        void                    SetQCnRecoveryType(qcRecovery type)             {fQCRecovery = type; }
        void                    SetModulationFitOptions(TString opt)            {fFitModulationOptions = opt; }
        void                    SetReferenceDetector(detectorType type)         {fDetectorType = type; }
        void                    SetCollisionType(collisionType type)            {fCollisionType = type; }
        void                    SetUsePtWeight(Bool_t w)                        {
            fUsePtWeight = w; 
            if(!fUsePtWeight) fUsePtWeightErrorPropagation = kFALSE; }
        void                    SetUsePtWeightErrorPropagation(Bool_t w)        {fUsePtWeightErrorPropagation = w; }
        void                    SetRunModeType(runModeType type)                {fRunModeType = type; }
        void                    SetAbsVertexZ(Float_t v)                        {fAbsVertexZ = v; }
        void                    SetMinDistanceRctoLJ(Float_t m)                 {fMinDisanceRCtoLJ = m; }
        void                    SetRandomConeRadius(Float_t r)                  {fRandomConeRadius = r; }
        void                    SetMaxNoRandomCones(Int_t m)                    {fMaxCones = m; }
        void                    SetMinLeadingHadronPt(Double_t m)               {fMinLeadingHadronPt = m; }
        void                    SetSetPtSub(Bool_t s)                           {fSubtractJetPt = s;}
        void                    SetForceAbsVnHarmonics(Bool_t f)                {fAbsVnHarmonics = f; }
        void                    SetExcludeLeadingJetsFromFit(Float_t n)         {fExcludeLeadingJetsFromFit = n; }
        void                    SetRebinSwapHistoOnTheFly(Bool_t r)             {fRebinSwapHistoOnTheFly = r; }
        void                    SetSaveThisPercentageOfFits(Float_t p)          {fPercentageOfFits = p; }
        void                    SetUseV0EventPlaneFromHeader(Bool_t h)          {fUseV0EventPlaneFromHeader = h;}
        void                    SetExplicitOutlierCutForYear(Int_t y)           {fExplicitOutlierCut = y;}
        // getters - these are used as well by AliAnalyisTaskJetFlow, so be careful when changing them
        TString                 GetJetsName() const                             {return GetJetContainer()->GetArrayName(); }
        TString                 GetTracksName() const                           {return GetParticleContainer()->GetArrayName(); }
        TString                 GetLocalRhoName() const                         {return fLocalRhoName; }
        TArrayD*                GetCentralityClasses() const                    {return fCentralityClasses;}
        TArrayD*                GetPtBinsHybrids() const                        {return fPtBinsHybrids; }
        TArrayD*                GetPtBinsJets() const                           {return fPtBinsJets; }
        TProfile*               GetResolutionParameters(Int_t h, Int_t c) const {return (h==2) ? fProfV2Resolution[c] : fProfV3Resolution[c];}
        TList*                  GetOutputList() const                           {return fOutputList;}
        AliLocalRhoParameter*   GetLocalRhoParameter() const                    {return fLocalRho;}
        Double_t                GetJetRadius() const                            {return GetJetContainer()->GetJetRadius();}
        /* inline */    AliEmcalJet* GetLeadingJet() {
            // return pointer to the highest pt jet (before background subtraction) within acceptance
            // only rudimentary cuts are applied on this level, hence the implementation outside of
            // the framework
            Int_t iJets(fJets->GetEntriesFast());
            Double_t pt(0);
            AliEmcalJet* leadingJet(0x0);
            for(Int_t i(0); i < iJets; i++) {
                AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(i));
                if(!PassesSimpleCuts(jet)) continue;
                if(jet->Pt() > pt) {
                   leadingJet = jet;
                   pt = leadingJet->Pt();
                }
            }
            return leadingJet;
        }
        void                    ExecMe()                                        {ExecOnce();}
        AliAnalysisTaskRhoVnModulation* ReturnMe()                              {return this;}
        // local cuts
        void                    SetLocalJetMinMaxEta(Float_t min, Float_t max)  {fLocalJetMinEta = min; fLocalJetMaxEta = max;}
        void                    SetLocalJetMinMaxEta(Float_t R)                 {fLocalJetMinEta = - 0.9 + R; fLocalJetMaxEta = 0.9 - R; }
        void                    SetLocalJetMinMaxPhi(Float_t min, Float_t max)  {fLocalJetMinPhi = min; fLocalJetMaxEta = max;}
        void                    SetSoftTrackMinMaxPt(Float_t min, Float_t max)  {fSoftTrackMinPt = min; fSoftTrackMaxPt = max;}
        void                    SetSemiGoodJetMinMaxPhi(Double_t a, Double_t b) {fSemiGoodJetMinPhi = a; fSemiGoodJetMaxPhi = b;}
        void                    SetSemiGoodTrackMinMaxPhi(Double_t a, Double_t b)       {fSemiGoodTrackMinPhi = a; fSemiGoodTrackMaxPhi = b;}
        // numerical evaluations
        void                    CalculateEventPlaneVZERO(Double_t vzero[2][2]) const;
        void                    CalculateEventPlaneTPC(Double_t* tpc);
        void                    CalculateEventPlaneCombinedVZERO(Double_t* comb) const;
        void                    CalculateEventPlaneResolution(Double_t vzero[2][2], Double_t* vzeroComb, Double_t* tpc);
        Double_t                CalculateEventPlaneChi(Double_t resEP) const;
        void                    CalculateRandomCone(Float_t &pt, Float_t &eta, Float_t &phi, AliEmcalJet* jet = 0x0) const;
        Double_t                CalculateQC2(Int_t harm);
        Double_t                CalculateQC4(Int_t harm);
        // helper calculations for the q-cumulant analysis, also used by AliAnalyisTaskJetFlow
        void                    QCnQnk(Int_t n, Int_t k, Double_t &reQ, Double_t &imQ);
        void                    QCnDiffentialFlowVectors(
            TClonesArray* pois, TArrayD* ptBins, Bool_t vpart, Double_t* repn, Double_t* impn, 
            Double_t *mp, Double_t *reqn, Double_t *imqn, Double_t* mq, Int_t n);
        Double_t                QCnS(Int_t i, Int_t j);
        Double_t                QCnM();
        Double_t                QCnM11();
        Double_t                QCnM1111();
        Bool_t                  QCnRecovery(Double_t psi2, Double_t psi3);
        // analysis details
        Bool_t                  CorrectRho(Double_t psi2, Double_t psi3);
        // event and track selection, also used by AliAnalyisTaskJetFlow
        /* inline */    Bool_t PassesCuts(AliVTrack* track) const { return AcceptTrack(track, 0); }
        /* inline */    Bool_t PassesCuts(AliEmcalJet* jet) { return AcceptJet(jet, 0); }
        /* inline */    Bool_t PassesSimpleCuts(AliEmcalJet* jet) {
            Float_t minPhi(GetJetContainer()->GetJetPhiMin()), maxPhi(GetJetContainer()->GetJetPhiMax());
            return (jet && jet->Pt() > 1 && jet->Eta() < .9-GetJetRadius() && jet->Eta() > -.9+GetJetRadius() && jet->Phi() > minPhi && jet->Phi() < maxPhi && jet->Area() > .557*GetJetRadius()*GetJetRadius()*TMath::Pi());
        }
        Bool_t                  PassesCuts(AliVEvent* event);
        Bool_t                  PassesCuts(Int_t year);
        Bool_t                  PassesCuts(const AliVCluster* track) const;
        // filling histograms
        void                    FillHistogramsAfterSubtraction(Double_t psi2, Double_t psi3, Double_t vzero[2][2], Double_t* vzeroComb, Double_t* tpc);
        void                    FillTrackHistograms() const;
        void                    FillClusterHistograms() const;
        void                    FillCorrectedClusterHistograms() const;
        void                    FillEventPlaneHistograms(Double_t vzero[2][2], Double_t* vzeroComb, Double_t* tpc) const;
        void                    FillRhoHistograms();
        void                    FillDeltaPtHistograms(Double_t psi2, Double_t psi3) const; 
        void                    FillJetHistograms(Double_t psi2, Double_t psi3);
        void                    FillQAHistograms(AliVTrack* vtrack) const;
        void                    FillQAHistograms(AliVEvent* vevent);
        void                    FillAnalysisSummaryHistogram() const;
        virtual void            Terminate(Option_t* option);
        // interface methods for the output file
        void                    SetOutputList(TList* l) {fOutputList = l;}
        TH1F*                   GetResolutionFromOuptutFile(detectorType detector, Int_t h = 2, TArrayD* c = 0x0);
        TH1F*                   CorrectForResolutionDiff(TH1F* v, detectorType detector, TArrayD* cen, Int_t c, Int_t h = 2);
        TH1F*                   CorrectForResolutionInt(TH1F* v, detectorType detector, TArrayD* cen, Int_t h = 2);
        TH1F*                   GetDifferentialQC(TProfile* refCumulants, TProfile* diffCumlants, TArrayD* ptBins, Int_t h);
    private:
        // analysis flags and settings
        Int_t                   fDebug;                 // debug level (0 none, 1 fcn calls, 2 verbose)
        Bool_t                  fRunToyMC;              // run toy mc for fit routine
        Bool_t                  fLocalInit;             //! is the analysis initialized?
        Bool_t                  fAttachToEvent;         // attach local rho to the event
        Bool_t                  fSemiCentralInclusive;  // semi central inclusive event selection
        Bool_t                  fFillHistograms;        // fill histograms
        Bool_t                  fFillQAHistograms;      // fill qa histograms
        Float_t                 fReduceBinsXByFactor;   // reduce the bins on x-axis of histo's by this much
        Float_t                 fReduceBinsYByFactor;   // reduce the bins on y-axis of histo's by this much
        Bool_t                  fNoEventWeightsForQC;   // don't store event weights for qc analysis
        TArrayD*                fCentralityClasses;     //-> centrality classes (maximum 10)
        TArrayD*                fPtBinsHybrids;         //-> pt bins for hybrid track vn anaysis
        TArrayD*                fPtBinsJets;            //-> pt bins for jet vn analysis
        TArrayI*                fExpectedRuns;          //-> array of expected run numbers, used for QA
        TArrayI*                fExpectedSemiGoodRuns;  //-> array of expected semi-good runs, used for cuts and QA
        TH1F*                   fUserSuppliedV2;        // histo with integrated v2
        TH1F*                   fUserSuppliedV3;        // histo with integrated v3
        TH1F*                   fUserSuppliedR2;        // correct the extracted v2 with this r
        TH1F*                   fUserSuppliedR3;        // correct the extracted v3 with this r
        AliParticleContainer*   fTracksCont;            //!tracks
        AliJetContainer*        fJetsCont;              //!jets
        AliEmcalJet*            fLeadingJet;            //! leading jet
        // members
        Bool_t                  fUseScaledRho;          // use scaled rho
        Int_t                   fNAcceptedTracks;       //! number of accepted tracks
        Int_t                   fNAcceptedTracksQCn;    //! accepted tracks for QCn
        fitModulationType       fFitModulationType;     // fit modulation type
        fitGoodnessTest         fFitGoodnessTest;       // fit goodness test type
        qcRecovery              fQCRecovery;            // recovery type for e-by-e qc method
        Bool_t                  fUsePtWeight;           // use dptdphi instead of dndphi
        Bool_t                  fUsePtWeightErrorPropagation;   // recalculate the bin errors in case of pt weighting 
        detectorType            fDetectorType;          // type of detector used for modulation fit
        TString                 fFitModulationOptions;  // fit options for modulation fit
        runModeType             fRunModeType;           // run mode type 
        dataType                fDataType;              // datatype 
        collisionType           fCollisionType;         // collision type
        TRandom3*               fRandom;                //-> dont use gRandom to not interfere with other tasks
        Int_t                   fRunNumber;             //! current runnumber (for QA and jet, track selection)
        Int_t                   fMappedRunNumber;       //! mapped runnumer (for QA)
        Int_t                   fInCentralitySelection; //! centrality bin
        TF1*                    fFitModulation;         //-> modulation fit for rho
        TF1*                    fFitControl;            //-> control fit
        Float_t                 fMinPvalue;             // minimum value of p
        Float_t                 fMaxPvalue;             // maximum value of p
        const char*             fNameJetClones;         //! collection of tclones array with jets
        const char*             fNamePicoTrackClones;   //! collection of tclones with pico tracks
        const char*             fNameRho;               //! name of rho
        TString                 fNameSmallRho;          // name of small rho
        AliRhoParameter*        fCachedRho;             //! temp cache for rho pointer
        // additional jet cuts (most are inherited)
        Float_t                 fLocalJetMinEta;        // local eta cut for jets
        Float_t                 fLocalJetMaxEta;        // local eta cut for jets
        Float_t                 fLocalJetMinPhi;        // local phi cut for jets
        Float_t                 fLocalJetMaxPhi;        // local phi cut for jets
        Float_t                 fSoftTrackMinPt;        // min pt for soft tracks
        Float_t                 fSoftTrackMaxPt;        // max pt for soft tracks
        Double_t                fSemiGoodJetMinPhi;     // min phi for semi good tpc runs
        Double_t                fSemiGoodJetMaxPhi;     // max phi for semi good tpc runs
        Double_t                fSemiGoodTrackMinPhi;   // min phi for semi good tpc runs
        Double_t                fSemiGoodTrackMaxPhi;   // max phi for semi good tpc runs
        // event cuts
        Float_t                 fAbsVertexZ;            // cut on zvertex
        // general qa histograms
        TH1F*                   fHistCentrality;        //! accepted centrality
        TH1F*                   fHistVertexz;           //! accepted verte
        TH2F*                   fHistRunnumbersPhi;     //! run numbers averaged phi
        TH2F*                   fHistRunnumbersEta;     //! run numbers averaged eta
        TH1F*                   fHistPvalueCDFROOT;     //! pdf value of chisquare p
        TH2F*                   fHistPvalueCDFROOTCent; //! p value versus centrlaity from root
        TH2F*                   fHistChi2ROOTCent;      //! reduced chi2 from ROOT, centrality correlation
        TH2F*                   fHistPChi2Root;         //! correlation p value and reduced chi2
        TH1F*                   fHistPvalueCDF;         //! cdf value of chisquare p
        TH2F*                   fHistPvalueCDFCent;     //! p value vs centrality
        TH2F*                   fHistChi2Cent;          //! reduced chi2, centrlaity correlation
        TH2F*                   fHistPChi2;             //! correlation p value and reduced chi2
        TH1F*                   fHistKolmogorovTest;    //! KolmogorovTest value
        TH2F*                   fHistKolmogorovTestCent;//! KolmogorovTest value, centrality correlation
        TH2F*                   fHistPKolmogorov;       //! p value vs kolmogorov value
        TH2F*                   fHistRhoStatusCent;     //! status of rho as function of centrality
        TH1F*                   fHistUndeterminedRunQA; //! undetermined run QA
        // general settings
        Float_t                 fMinDisanceRCtoLJ;      // min distance between rc and leading jet
        Float_t                 fRandomConeRadius;      // radius of random cone
        Int_t                   fMaxCones;              // max number of random cones
        Bool_t                  fAbsVnHarmonics;        // force postive local rho
        Float_t                 fExcludeLeadingJetsFromFit;    // exclude n leading jets from fit
        Bool_t                  fRebinSwapHistoOnTheFly;       // rebin swap histo on the fly
        Float_t                 fPercentageOfFits;      // save this percentage of fits
        Bool_t                  fUseV0EventPlaneFromHeader;             // use the vzero event plane from the header
        Int_t                   fExplicitOutlierCut;    // cut on correlation of tpc and global multiplicity
        Double_t                fMinLeadingHadronPt;    // minimum pt for leading hadron
        Bool_t                  fSubtractJetPt;         // save subtracted jet pt by calling SetPtSub
        // transient object pointers
        TList*                  fOutputList;            //! output list
        TList*                  fOutputListGood;        //! output list for local analysis
        TList*                  fOutputListBad;         //! output list for local analysis
        TH1F*                   fHistAnalysisSummary;   //! analysis summary
        TH1F*                   fHistSwap;              //! swap histogram
        TProfile*               fProfV2;                //! extracted v2
        TProfile*               fProfV2Cumulant;        //! v2 cumulant
        TProfile*               fProfV2Resolution[10];  //! resolution parameters for v2
        TProfile*               fProfV3;                //! extracted v3
        TProfile*               fProfV3Cumulant;        //! v3 cumulant
        TProfile*               fProfV3Resolution[10];  //! resolution parameters for v3
        // qa histograms for accepted pico tracks
        TH1F*                   fHistPicoTrackPt[10];    //! pt of all charged tracks
        TH1F*                   fHistPicoTrackMult[10];  //! multiplicity of accepted pico tracks
        TH2F*                   fHistPicoCat1[10];       //! pico tracks spd hit and refit
        TH2F*                   fHistPicoCat2[10];       //! pico tracks wo spd hit w refit, constrained
        TH2F*                   fHistPicoCat3[10];       //! pico tracks wo spd hit wo refit, constrained
        // qa histograms for accepted emcal clusters
        /* TH1F*                   fHistClusterPt[10];      //! pt uncorrected emcal clusters */
        /* TH1F*                   fHistClusterPhi[10];     //! phi uncorrected emcal clusters */
        /* TH1F*                   fHistClusterEta[10];     //! eta uncorrected emcal clusters */
        // qa histograms for accepted emcal clusters aftehadronic correction
        /* TH1F*                   fHistClusterCorrPt[10];  //! pt corrected emcal clusters */
        /* TH1F*                   fHistClusterCorrPhi[10]; //! phi corrected emcal clusters */
        /* TH1F*                   fHistClusterCorrEta[10]; //! eta corrected emcal clusters */
        // qa event planes
        TProfile*               fHistPsiControl;         //! event plane control histogram
        TProfile*               fHistPsiSpread;          //! event plane spread histogram
        TH1F*                   fHistPsiVZEROA;          //! psi 2 from vzero a
        TH1F*                   fHistPsiVZEROC;          //! psi 2 from vzero c
        TH1F*                   fHistPsiVZERO;           //! psi 2 from combined vzero
        TH1F*                   fHistPsiTPC;             //! psi 2 from tpc
        TH2F*                   fHistPsiVZEROAV0M;      //! psi 2 from vzero a
        TH2F*                   fHistPsiVZEROCV0M;      //! psi 2 from vzero c
        TH2F*                   fHistPsiVZEROVV0M;      //! psi 2 from combined vzero
        TH2F*                   fHistPsiTPCiV0M;        //! psi 2 from tpc
        TH2F*                   fHistPsiVZEROATRK;      //! psi 2 from vzero a
        TH2F*                   fHistPsiVZEROCTRK;      //! psi 2 from vzero c
        TH2F*                   fHistPsiVZEROTRK;       //! psi 2 from combined vzero
        TH2F*                   fHistPsiTPCTRK;         //! psi 2 from tpc
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
        TH2F*                   fHistDeltaPtDeltaPhi2[10];      //! dpt vs dphi (psi2 - phi)
        TH2F*                   fHistDeltaPtDeltaPhi3[10];      //! dpt vs dphi (psi3 - phi)
        TH2F*                   fHistRCPhiEtaExLJ[10];          //! random cone eta and phi, excl leading jet
        TH2F*                   fHistRhoVsRCPtExLJ[10];         //! rho * A vs rcpt, excl leading jet
        TH1F*                   fHistRCPtExLJ[10];              //! rcpt, excl leading jet
        TH2F*                   fHistDeltaPtDeltaPhi2ExLJ[10];  //! dpt vs dphi, excl leading jet
        TH2F*                   fHistDeltaPtDeltaPhi3ExLJ[10];  //! dpt vs dphi, excl leading jet
        /* TH2F*                   fHistRCPhiEtaRand[10];          //! random cone eta and phi, randomized */
        /* TH2F*                   fHistRhoVsRCPtRand[10];         //! rho * A vs rcpt, randomized */
        /* TH1F*                   fHistRCPtRand[10];              //! rcpt, randomized */ 
        /* TH2F*                   fHistDeltaPtDeltaPhi2Rand[10];  //! dpt vs dphi, randomized */
        /* TH2F*                   fHistDeltaPtDeltaPhi3Rand[10];  //! dpt vs dphi, randomized */
        // jet histograms (after kinematic cuts)
        TH1F*                   fHistJetPtRaw[10];              //! jet pt - no background subtraction
        TH1F*                   fHistJetPt[10];                 //! pt of found jets (background subtracted)
        TH2F*                   fHistJetEtaPhi[10];             //! eta and phi correlation
        TH2F*                   fHistJetPtArea[10];             //! jet pt versus area
        TH2F*                   fHistJetPtConstituents[10];     //! jet pt versus number of constituents
        TH2F*                   fHistJetEtaRho[10];             //! jet eta versus jet rho
        // in plane, out of plane jet spectra
        TH2F*                   fHistJetPsi2Pt[10];             //! psi tpc versus pt
        TH2F*                   fHistJetPsi3Pt[10];             //! psi vzeroc versus pt

        AliAnalysisTaskRhoVnModulation(const AliAnalysisTaskRhoVnModulation&);                  // not implemented
        AliAnalysisTaskRhoVnModulation& operator=(const AliAnalysisTaskRhoVnModulation&);       // not implemented

        ClassDef(AliAnalysisTaskRhoVnModulation, 25);
};

#endif
