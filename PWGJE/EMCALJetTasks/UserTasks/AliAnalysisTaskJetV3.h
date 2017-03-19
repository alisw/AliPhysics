/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskJetV3_H
#define AliAnalysisTaskJetV3_H

// uncomment or define externally to enable debug information
// 
// flag for global and e-by-e debug info
// #define ALIANALYSISTASKJETV3_DEBUG_FLAG_1
// flag for debug statements that may be repeated multiple times per event
// #define ALIANALYSISTASKJETV3_DEBUG_FLAG_2

#include <AliAnalysisTaskEmcalJet.h>
#include <AliEmcalJet.h>
#include <AliVEvent.h>
#include <AliVParticle.h>
#include <AliVCluster.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TArrayD.h>
#include <TRandom3.h>
#include <AliJetContainer.h>
#include <AliTrackContainer.h>

class TFile;
class TF1;
class TH1F;
class TH1I;
class TH2F;
class TH3F;
class TProfile;
class AliLocalRhoParameter;
class AliClusterContainer;
class AliVTrack;

class AliAnalysisTaskJetV3 : public AliAnalysisTaskEmcalJet {
    public:
         // enumerators
        enum fitModulationType  { kNoFit, kV2, kV3, kCombined, kFourierSeries, kIntegratedFlow, kQC2, kQC4 }; // fit type
        enum fitGoodnessTest    { kChi2ROOT, kChi2Poisson, kKolmogorov, kKolmogorovTOY, kLinearFit };
        enum collisionType      { kPbPb, kPythia, kPbPb10h, kPbPb11h, kJetFlowMC, kPbPb15o }; // collision type, kPbPb = 11h, kept for backward compatibilitiy
        enum qcRecovery         { kFixedRho, kNegativeVn, kTryFit };    // how to deal with negative cn value for qcn value
        enum runModeType        { kLocal, kGrid };                      // run mode type
        enum dataType           { kESD, kAOD, kESDMC, kAODMC};          // data type
        enum detectorType       { kTPC, kVZEROA, kVZEROC, kVZEROComb, kFixedEP};  // detector that was used for event plane
        enum EPweightType       { kNone, kChi, kSigmaSquared};          // event plane weight type
        enum analysisType       { kCharged, kFull };                    // analysis type
        // constructors, destructor
                                AliAnalysisTaskJetV3();
                                AliAnalysisTaskJetV3(
                                        const char *name,               // task name
                                        runModeType type,               // grid or local mode
                                        Bool_t baseClassHistos = kFALSE);       // book framework histos
        virtual                 ~AliAnalysisTaskJetV3();
        // setting up the task and technical aspects
        void                    ExecOnce();
        virtual Bool_t          Notify();
        Bool_t                  InitializeAnalysis();
        virtual void            UserCreateOutputObjects();
        virtual void            Exec(Option_t *);
        virtual Bool_t          Run();
        TH1F*                   BookTH1F(const char* name, const char* x, Int_t bins, Double_t min, Double_t max, Int_t c = -1, Bool_t append = kTRUE);
        TH2F*                   BookTH2F(const char* name, const char* x, const char* y, Int_t binsx, Double_t minx, Double_t maxx, Int_t binsy, Double_t miny, Double_t maxy, Int_t c = -1, Bool_t append = kTRUE);
        TH3F*                   BookTH3F(const char* name, const char* x, const char* y, const char* z, Int_t binsx, Double_t minx, Double_t maxx, Int_t binsy, Double_t miny, Double_t maxy, Int_t binsz, Double_t minz, Double_t maxz, Int_t c = -1, Bool_t append = kTRUE);
        /* inline */    static Double_t PhaseShift(Double_t x) {  
            while (x>=TMath::TwoPi())x-=TMath::TwoPi();
            while (x<0.)x+=TMath::TwoPi();
            return x; }
        /* inline */    static Double_t PhaseShift(Double_t x, Double_t n) {
            x = PhaseShift(x);
            if(TMath::Nint(n)==2) while (x>TMath::Pi()) x-=TMath::Pi();
            if(TMath::Nint(n)==3) {
                if(x>2.*TMath::TwoPi()/n) x = TMath::TwoPi() - x;
                if(x>TMath::TwoPi()/n) x = TMath::TwoPi()-(x+TMath::TwoPi()/n);
            }
            return x; }
        /* inline */    static Bool_t   IsInPlane(Double_t dPhi) {
            return (dPhi < -1.*TMath::Pi()/4. || dPhi > TMath::Pi()/4.); }
        /* inline */    static Double_t ChiSquarePDF(Int_t ndf, Double_t x) {
            Double_t n(ndf/2.), denom(TMath::Power(2, n)*TMath::Gamma(n));
            if (denom!=0)  return ((1./denom)*TMath::Power(x, n-1)*TMath::Exp(-x/2.)); 
            return -999; }
        // note that the cdf of the chisquare distribution is the normalized lower incomplete gamma function
        /* inline */    static Double_t ChiSquareCDF(Int_t ndf, Double_t x) { return TMath::Gamma(ndf/2., x/2.); }
        /* inline */    static Double_t ChiSquare(TH1& histo, TF1* func) {
            // evaluate the chi2 using a poissonian error estimate on bins
            Double_t chi2(0.);
            for(Int_t i(0); i < histo.GetXaxis()->GetNbins(); i++) {
                if(histo.GetBinContent(i+1) <= 0.) continue;
                chi2 += TMath::Power((histo.GetBinContent(i+1)-func->Eval(histo.GetXaxis()->GetBinCenter(1+i))), 2)/histo.GetBinContent(i+1);
            }
           return chi2;
        }
        /* inline */ Double_t KolmogorovTest(/*TH1F& histo, TF1* func*/) const {
            // return the probability from a Kolmogorov test
            return .5;
            /* this test is disabeled as it eats a lot of resources but kept as a dummty to 
             * ensure compatibility of the output with offline macros
            TH1F test(histo);       // stack copy of test statistic
            for(Int_t i(0); i < test.GetXaxis()->GetNbins(); i++) test.SetBinContent(i+1, func->Eval(test.GetXaxis()->GetBinCenter(1+i)));
            if(fFitGoodnessTest == kKolmogorovTOY) return histo.TH1::KolmogorovTest((&test), "X");
            return histo.TH1::KolmogorovTest((&test));
            */
        }

        // setters - analysis setup
        void                    SetRunToyMC(Bool_t t)                           {fRunToyMC = t; }
        void                    SetAttachToEvent(Bool_t b)                      {fAttachToEvent = b;}
        void                    SetFillHistograms(Bool_t b)                     {fFillHistograms = b;}
        void                    SetFillQAHistograms(Bool_t qa)                  {fFillQAHistograms = qa;}
        void                    SetReduceBinsXYByFactor(Float_t x, Float_t y)   {fReduceBinsXByFactor = x;
                                                                                 fReduceBinsYByFactor = y;}
        void                    SetNoEventWeightsForQC(Bool_t e)                {fNoEventWeightsForQC = e;}
        void                    SetCentralityClasses(TArrayD* c)                {fCentralityClasses = c;}
        void                    SetExpectedRuns(TArrayI* r)                     {fExpectedRuns = r;}
        void                    SetExpectedSemiGoodRuns(TArrayI* r)             {fExpectedSemiGoodRuns = r;}
        void                    SetIntegratedFlow(TH1F* i, TH1F* j)             {fUserSuppliedV2 = i;
                                                                                 fUserSuppliedV3 = j; }
        void                    SetOnTheFlyResCorrection(TH1F* r2, TH1F* r3)    {fUserSuppliedR2 = r2;
                                                                                 fUserSuppliedR3 = r3; }
        void                    SetEventPlaneWeights(TH1F* ep, Int_t c)         {fEventPlaneWeights[c] = ep; }
        void                    SetAcceptanceWeights(Bool_t w)                  {fAcceptanceWeights = w; }
        void                    SetNameRhoSmall(TString name)                   {fNameSmallRho = name; }
        void                    SetRandomSeed(TRandom3* r)                      {if (fRandom) delete fRandom; fRandom = r; }
        void                    SetModulationFit(TF1* fit);
        void                    SetUseControlFit(Bool_t c);
        void                    SetModulationFitMinMaxP(Float_t m, Float_t n)   {fMinPvalue = m; fMaxPvalue = n; }
        void                    SetModulationFitType(fitModulationType type)    {fFitModulationType = type; }
        void                    SetGoodnessTest(fitGoodnessTest test)           {fFitGoodnessTest = test; }
        void                    SetQCnRecoveryType(qcRecovery type)             {fQCRecovery = type; }
        void                    SetModulationFitOptions(TString opt)            {fFitModulationOptions = opt; }
        void                    SetReferenceDetector(detectorType type)         {fDetectorType = type; }
        void                    SetAnalysisType(analysisType type)              {fAnalysisType = type; }
        void                    SetCollisionType(collisionType type)            {fCollisionType = type; }
        void                    SetUsePtWeight(Bool_t w)                        {
            fUsePtWeight = w; 
            if(!fUsePtWeight) fUsePtWeightErrorPropagation = kFALSE; }
        void                    SetUsePtWeightErrorPropagation(Bool_t w)        {fUsePtWeightErrorPropagation = w; }
        void                    SetUse2DIntegration(Bool_t k)                   {fUse2DIntegration= k;}
        void                    SetRunModeType(runModeType type)                {fRunModeType = type; }
        void                    SetMinDistanceRctoLJ(Float_t m)                 {fMinDisanceRCtoLJ = m; }
        void                    SetMaxNoRandomCones(Int_t m)                    {fMaxCones = m; }
        void                    SetExcludeLeadingJetsFromFit(Float_t n)         {fExcludeLeadingJetsFromFit = n; }
        void                    SetExcludeJetsWithTrackPt(Float_t n)            {fExcludeJetsWithTrackPt = n; }
        void                    SetRebinSwapHistoOnTheFly(Bool_t r)             {fRebinSwapHistoOnTheFly = r; }
        void                    SetSaveThisPercentageOfFits(Float_t p)          {fPercentageOfFits = p; }
        void                    SetChi2VZEROA(TArrayD* a)                       { fChi2A = a;}
        void                    SetChi2VZEROC(TArrayD* a)                       { fChi2C = a;}
        void                    SetChi3VZEROA(TArrayD* a)                       { fChi3A = a;}
        void                    SetChi3VZEROC(TArrayD* a)                       { fChi3C = a;}
        void                    SetSigma2VZEROA(TArrayD* a)                     { fSigma2A = a;}
        void                    SetSigma2VZEROC(TArrayD* a)                     { fSigma2C = a;}
        void                    SetSigma3VZEROA(TArrayD* a)                     { fSigma3A = a;}
        void                    SetSigma3VZEROC(TArrayD* a)                     { fSigma3C = a;}
        void                    SetEPWeightForVZERO(EPweightType type)          { fWeightForVZERO = type; }
        // getters 
        TString                 GetJetsName() const                             {return GetJetContainer()->GetArrayName(); }
        TString                 GetTracksName() const                           {return GetTrackContainer()->GetArrayName(); }
        TString                 GetLocalRhoName() const                         {return fLocalRhoName; }
        TArrayD*                GetCentralityClasses() const                    {return fCentralityClasses;}
//        Float_t                 GetCentrality() const;
        TProfile*               GetResolutionParameters(Int_t h, Int_t c) const {return (h==2) ? fProfV2Resolution[c] : fProfV3Resolution[c];}
        TList*                  GetOutputList() const                           {return fOutputList;}
        AliLocalRhoParameter*   GetLocalRhoParameter() const                    {return fLocalRho;}
        Double_t                GetJetRadius() const                            {return GetJetContainer()->GetJetRadius();}
        AliEmcalJet*            GetLeadingJet(AliLocalRhoParameter* localRho = 0x0);
        AliVParticle*           GetLeadingTrack(AliEmcalJet* jet);
        static TH1F*            GetEventPlaneWeights(TH1F* hist, Int_t c);
        static void             PrintTriggerSummary(UInt_t trigger);
        static void             DoSimpleSimulation(Int_t nEvents = 100000, Float_t v2 = 0.02, Float_t v3 = 0.04, Float_t v4 = 0.03);
        void                    ExecMe()                                        {ExecOnce();}
        AliAnalysisTaskJetV3*   ReturnMe()                                      {return this;}
        // local cuts
        void                    SetSoftTrackMinMaxPt(Float_t min, Float_t max)          {fSoftTrackMinPt = min; fSoftTrackMaxPt = max;}
        void                    SetSemiGoodJetMinMaxPhi(Double_t a, Double_t b)         {fSemiGoodJetMinPhi = a; fSemiGoodJetMaxPhi = b;}
        void                    SetSemiGoodTrackMinMaxPhi(Double_t a, Double_t b)       {fSemiGoodTrackMinPhi = a; fSemiGoodTrackMaxPhi = b;}
        // numerical evaluations
        static void             NumericalOverlap(Double_t x1, Double_t x2, Double_t psi2, Double_t &percIn, Double_t &percOut, Double_t &percLost);
        static Int_t            OverlapsWithPlane(Double_t x1, Double_t x2, 
                Double_t a, Double_t b, Double_t c, Double_t d, Double_t e, Double_t phi);
        static Double_t         CalculateEventPlaneChi(Double_t res);
        void                    CalculateEventPlaneVZERO(Double_t vzero[2][2]) const;
        void                    CalculateEventPlaneCombinedVZERO(Double_t* comb) const;
        void                    CalculateEventPlaneTPC(Double_t* tpc);
        void                    CalculateEventPlaneResolution(Double_t vzero[2][2], Double_t* vzeroComb, Double_t* tpc);
        void                    CalculateQvectorVZERO(Double_t Qa2[2], Double_t Qc2[2], Double_t Qa3[2], Double_t Qc3[2]) const;
        void                    CalculateQvectorCombinedVZERO(Double_t Q2[2], Double_t Q3[2]) const;
        void                    CalculateRandomCone(
                Float_t &pt, 
                Float_t &eta, 
                Float_t &phi, 
                AliTrackContainer* tracksCont,
                AliClusterContainer* clusterCont = 0x0,
                AliEmcalJet* jet = 0x0
                ) const;
        Double_t                CalculateQC2(Int_t harm);
        Double_t                CalculateQC4(Int_t harm);
        // helper calculations for the q-cumulant analysis
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
        // event and track selection
        /* inline */    Bool_t PassesCuts(AliVParticle* track) const    { UInt_t rejectionReason = 0; return GetTrackContainer(0)->AcceptParticle(track, rejectionReason); }
        /* inline */    Bool_t PassesCuts(AliEmcalJet* jet)             { 
            if(jet->MaxTrackPt() > fExcludeJetsWithTrackPt) return kFALSE;
            return AcceptJet(jet, 0); 
        }
        /* inline */    Bool_t PassesCuts(AliVCluster* clus) const      { return AcceptCluster(clus, 0); }
        /* inline */    Bool_t PassesSimpleCuts(AliEmcalJet* jet)       {
            Float_t minPhi(GetJetContainer()->GetJetPhiMin()), maxPhi(GetJetContainer()->GetJetPhiMax());
            Float_t minEta(GetJetContainer()->GetJetEtaMin()), maxEta(GetJetContainer()->GetJetEtaMax());
            return (jet/* && jet->Pt() > 1.*/ && jet->Eta() > minEta && jet->Eta() < maxEta && jet->Phi() > minPhi && jet->Phi() < maxPhi && jet->Area() > .557*GetJetRadius()*GetJetRadius()*TMath::Pi());
        }
        Bool_t                  PassesCuts(AliVEvent* event);
        Bool_t                  PassesExperimentalHighLumiCuts(AliAODEvent* event);
        Bool_t                  MultiVertexer(const AliAODEvent* event);
        Double_t                GetWDist(const AliVVertex* v0, const AliVVertex* v1);
        Bool_t                  PassesCuts(const AliVCluster* track) const;
        // filling histograms
        void                    FillHistogramsAfterSubtraction(Double_t psi3, Double_t vzero[2][2], Double_t* vzeroComb, Double_t* tpc);
        void                    FillQAHistograms(AliVTrack* vtrack) const;
        void                    FillQAHistograms(AliVEvent* vevent);
        void                    FillWeightedTrackHistograms() const;
        void                    FillWeightedClusterHistograms() const;
        void                    FillWeightedEventPlaneHistograms(Double_t vzero[2][2], Double_t* vzeroComb, Double_t* tpc) const;
        void                    FillWeightedRhoHistograms();
        void                    FillWeightedDeltaPtHistograms(Double_t psi3) const; 
        void                    FillWeightedJetHistograms(Double_t psi3);
        void                    FillWeightedQAHistograms(AliVTrack* vtrack) const;
        void                    FillWeightedQAHistograms(AliVEvent* vevent);
        void                    FillWeightedTriggerQA(Double_t dPhi, Double_t pt, UInt_t trigger);
        void                    FillAnalysisSummaryHistogram() const;
        virtual void            Terminate(Option_t* option);
        // interface methods for the output file
        void                    SetOutputList(TList* l) {fOutputList = l;}
        TH1F*                   GetResolutionFromOutputFile(detectorType detector, Int_t h = 2, TArrayD* c = 0x0);
        TH1F*                   CorrectForResolutionDiff(TH1F* v, detectorType detector, TArrayD* cen, Int_t c, Int_t h = 2);
        TH1F*                   CorrectForResolutionInt(TH1F* v, detectorType detector, TArrayD* cen, Int_t h = 2);
        TH1F*                   GetDifferentialQC(TProfile* refCumulants, TProfile* diffCumlants, TArrayD* ptBins, Int_t h);
        void                    ReadVZEROCalibration2010h();
        void                    ReadVZEROCalibration2011h();
        void                    ReadVZEROCalibration2015o();
        Int_t                   GetVZEROCentralityBin() const;
        Float_t                 GetCentrality(const char* estimator) const;
    private:
        // analysis flags and settings
        Bool_t                  fRunToyMC;              // run toy mc for fit routine
        Bool_t                  fLocalInit;             //! is the analysis initialized?
        Bool_t                  fAttachToEvent;         // attach local rho to the event
        Bool_t                  fFillHistograms;        // fill histograms
        Bool_t                  fFillQAHistograms;      // fill qa histograms
        Float_t                 fReduceBinsXByFactor;   // reduce the bins on x-axis of histo's by this much
        Float_t                 fReduceBinsYByFactor;   // reduce the bins on y-axis of histo's by this much
        Bool_t                  fNoEventWeightsForQC;   // don't store event weights for qc analysis
        TArrayD*                fCentralityClasses;     //-> centrality classes (maximum 10)
        TArrayI*                fExpectedRuns;          //-> array of expected run numbers, used for QA
        TArrayI*                fExpectedSemiGoodRuns;  //-> array of expected semi-good runs, used for cuts and QA
        TH1F*                   fUserSuppliedV2;        // histo with integrated v2
        TH1F*                   fUserSuppliedV3;        // histo with integrated v3
        TH1F*                   fUserSuppliedR2;        // correct the extracted v2 with this r
        TH1F*                   fUserSuppliedR3;        // correct the extracted v3 with this r
        TH1F*                   fEventPlaneWeights[10]; // weight histos for the event plane (centrality dependent)
        Bool_t                  fAcceptanceWeights;     // store centrality dependent acceptance weights
        Float_t                 fEventPlaneWeight;      //! the actual weight of an event
        AliTrackContainer*   fTracksCont;            //! tracks
        AliClusterContainer*    fClusterCont;           //! cluster container
        AliJetContainer*        fJetsCont;              //! jets
        AliEmcalJet*            fLeadingJet;            //! leading jet
        AliEmcalJet*            fLeadingJetAfterSub;    //! leading jet after background subtraction
        // members
        Int_t                   fNAcceptedTracks;       //! number of accepted tracks
        Int_t                   fNAcceptedTracksQCn;    //! accepted tracks for QCn
        fitModulationType       fFitModulationType;     // fit modulation type
        fitGoodnessTest         fFitGoodnessTest;       // fit goodness test type
        qcRecovery              fQCRecovery;            // recovery type for e-by-e qc method
        Bool_t                  fUsePtWeight;           // use dptdphi instead of dndphi
        Bool_t                  fUsePtWeightErrorPropagation;   // recalculate the bin errors in case of pt weighting 
        Bool_t                  fUse2DIntegration;      // integrate jet background over eta, phi
        detectorType            fDetectorType;          // type of detector used for modulation fit
        analysisType            fAnalysisType;          // analysis type (full or charged jets)
        TString                 fFitModulationOptions;  // fit options for modulation fit
        runModeType             fRunModeType;           // run mode type 
        dataType                fDataType;              // datatype 
        collisionType           fCollisionType;         // collision type
        TRandom3*               fRandom;                //-> dont use gRandom to not interfere with other tasks
        Int_t                   fRunNumber;             //! current runnumber (for QA and jet, track selection)
        Int_t                   fRunNumberCaliInfo;     //! runnumber of the cached calibration info
        Int_t                   fMappedRunNumber;       //! mapped runnumer (for QA)
        Int_t                   fInCentralitySelection; //! centrality bin
        TF1*                    fFitModulation;         //-> modulation fit for rho
        TF1*                    fFitControl;            //-> control fit
        Float_t                 fMinPvalue;             // minimum value of p
        Float_t                 fMaxPvalue;             // maximum value of p
        TString                 fNameSmallRho;          // name of small rho
        AliRhoParameter*        fCachedRho;             //! temp cache for rho pointer
        // additional jet cuts (most are inherited)
        Float_t                 fSoftTrackMinPt;        // min pt for soft tracks
        Float_t                 fSoftTrackMaxPt;        // max pt for soft tracks
        Double_t                fSemiGoodJetMinPhi;     // min phi for semi good tpc runs
        Double_t                fSemiGoodJetMaxPhi;     // max phi for semi good tpc runs
        Double_t                fSemiGoodTrackMinPhi;   // min phi for semi good tpc runs
        Double_t                fSemiGoodTrackMaxPhi;   // max phi for semi good tpc runs
        // general qa histograms
        TH1F*                   fHistCentrality;        //! accepted centrality
        TProfile*               fHistCentralityPercIn;  //! centrality versus perc in
        TProfile*               fHistCentralityPercOut; //! centrality versus perc out
        TProfile*               fHistCentralityPercLost;//! centrality versus perc lost
        TH1F*                   fHistVertexz;           //! accepted verte
        TH2F*                   fHistMultCorAfterCuts;      //! QA profile global and tpc multiplicity after outlier cut
        TH2F*                   fHistMultvsCentr;           //! QA profile of centralty vs multiplicity
        TH2F*                   fHistRunnumbersPhi;     //! run numbers averaged phi
        TH2F*                   fHistRunnumbersEta;     //! run numbers averaged eta
        TH1I*                   fHistRunnumbersCaliInfo;//! calibration info per runnumber
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
        Int_t                   fMaxCones;              // max number of random cones
        Float_t                 fExcludeLeadingJetsFromFit;    // exclude n leading jets from fit
        Float_t                 fExcludeJetsWithTrackPt;// exclude jets with a track with pt higher than this
        Bool_t                  fRebinSwapHistoOnTheFly;       // rebin swap histo on the fly
        Float_t                 fPercentageOfFits;      // save this percentage of fits
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
        TH1F*                   fHistPicoTrackPt[10];   //! pt of all charged tracks
        TH1F*                   fHistPicoTrackMult[10]; //! multiplicity of accepted pico tracks
        TH2F*                   fHistPicoCat1[10];      //! pico tracks spd hit and refit
        TH2F*                   fHistPicoCat2[10];      //! pico tracks wo spd hit w refit, constrained
        TH2F*                   fHistPicoCat3[10];      //! pico tracks wo spd hit wo refit, constrained
        // qa histograms for accepted emcal clusters
        TH1F*                   fHistClusterPt[10];     //! pt emcal clusters
        TH2F*                   fHistClusterEtaPhi[10]; //! eta phi emcal clusters
        TH2F*                   fHistClusterEtaPhiWeighted[10]; //! eta phi emcal clusters, pt weighted
        // qa histograms for triggers
        TH2F*                   fHistTriggerQAIn[10];   //! trigger qa in plane
        TH2F*                   fHistTriggerQAOut[10];  //! trigger qa out of plane
        // qa event planes
        TH2F*                   fHistPsiVZEROAV0M;      //! psi 2 from vzero a
        TH2F*                   fHistPsiVZEROCV0M;      //! psi 2 from vzero c
        TH2F*                   fHistPsiVZEROVV0M;      //! psi 2 from combined vzero
        TH2F*                   fHistPsiTPCV0M;         //! psi 2 from tpc
        TH2F*                   fHistPsiVZEROATRK;      //! psi 2 from vzero a
        TH2F*                   fHistPsiVZEROCTRK;      //! psi 2 from vzero c
        TH2F*                   fHistPsiVZEROTRK;       //! psi 2 from combined vzero
        TH2F*                   fHistPsiTPCTRK;         //! psi 2 from tpc
        TH3F*                   fHistEPCorrelations[10];        //! ep correlations
        TH2F*                   fHistEPCorrAvChi[10];           //! ep corr
        TH2F*                   fHistEPCorrAvSigma[10];         //! ep corr
        TH2F*                   fHistEPCorrChiSigma[10];        //! ep corr
        TH2F*                   fHistIntegralCorrelations[10];  //! correlate polar or local integral
        TProfile*               fProfIntegralCorrelations[10];  //! same qa lot
        TH3F*                   fHistPsiTPCLeadingJet[10];      //! correlation tpc EP, LJ pt
        TH3F*                   fHistPsiVZEROALeadingJet[10];   //! correlation vzeroa EP, LJ pt
        TH3F*                   fHistPsiVZEROCLeadingJet[10];   //! correlation vzeroc EP, LJ pt
        TH3F*                   fHistPsiVZEROCombLeadingJet[10];//! correlation vzerocomb EP, LJ pt
        TH3F*                   fHistPsi3Correlation[10];       //! correlation of event planes
        TH2F*                   fHistLeadingJetBackground[10];  //! geometric correlation of leading jet w/wo bkg subtraction
        // background
        TH1F*                   fHistRhoPackage[10];    //! rho as estimated by emcal jet package
        TH1F*                   fHistRho[10];           //! background
        TH2F*                   fHistRhoVsMult;         //! rho versus multiplicity
        TH2F*                   fHistRhoVsCent;         //! rho veruss centrality
        TH2F*                   fHistRhoAVsMult;        //! rho * A vs multiplicity for all jets
        TH2F*                   fHistRhoAVsCent;        //! rho * A vs centrality for all jets
        TH2F*                   fHistRhoEtaBC[10];      //! rho vs eta before cuts
        // delta pt distributions
        TH2F*                   fHistRCPhiEta[10];              //! random cone eta and phi
        TH2F*                   fHistRhoVsRCPt[10];             //! rho * A vs rcpt
        TH1F*                   fHistRCPt[10];                  //! rcpt
        TH2F*                   fHistDeltaPtDeltaPhi3[10];      //! dpt vs dphi (psi2 - phi)
        TH2F*                   fHistDeltaPtDeltaPhi3Rho0[10];  //! dpt vs dphi, rho_0
        TH2F*                   fHistRCPhiEtaExLJ[10];          //! random cone eta and phi, excl leading jet
        TH2F*                   fHistRhoVsRCPtExLJ[10];         //! rho * A vs rcpt, excl leading jet
        TH1F*                   fHistRCPtExLJ[10];              //! rcpt, excl leading jet
        TH2F*                   fHistDeltaPtDeltaPhi3ExLJ[10];  //! dpt vs dphi, excl leading jet
        TH2F*                   fHistDeltaPtDeltaPhi3ExLJRho0[10];      //! dpt vs dphi, excl leading jet, rho_0
        // jet histograms (after kinematic cuts)
        TH1F*                   fHistJetPtRaw[10];              //! jet pt - no background subtraction
        TH1F*                   fHistJetPt[10];                 //! pt of found jets (background subtracted)
        TH1F*                   fHistJetPtBC[10];               //! jet pt before area cut
        TH2F*                   fHistJetEtaPhi[10];             //! eta and phi correlation
        TH2F*                   fHistJetEtaPhiBC[10];           //! eta and phi correlation before cuts
        TH2F*                   fHistJetPtArea[10];             //! jet pt versus area
        TH2F*                   fHistJetPtAreaBC[10];           //! jet pt versus area before cuts
        TH2F*                   fHistJetPtEta[10];              //! jet pt versus eta (temp control)
        TH2F*                   fHistJetPtConstituents[10];     //! jet pt versus number of constituents
        TH2F*                   fHistJetEtaRho[10];             //! jet eta versus rho
        // in plane, out of plane jet spectra
        TH2F*                   fHistJetPsi3Pt[10];             //! event plane dependence of jet pt
        TH3F*                   fHistJetLJPsi3Pt[10];           //! event plane dependence of jet pt and leading track pt
        TH3F*                   fHistJetLJPsi3PtRatio[10];      //! ratio of leading track v2 to jet v2
        TH2F*                   fHistJetPsi3PtRho0[10];         //! event plane dependence of jet pt vs rho_0
        // vzero event plane calibration cache for 10h data
        Float_t                 fMeanQ[9][2][2];                //! recentering
        Float_t                 fWidthQ[9][2][2];               //! recentering
        Float_t                 fMeanQv3[9][2][2];              //! recentering
        Float_t                 fWidthQv3[9][2][2];             //! recentering
        TH1*                    fMQ[2][2][2];                   //! recentering
        TH1*                    fWQ[2][2][2];                   //! recentering
        TH1*                    fVZEROgainEqualization;         //! equalization histo
        Float_t                 fVZEROApol;                     //! calibration info per disc
        Float_t                 fVZEROCpol;                     //! calibration info per disc
        TArrayD*                fChi2A;                         // chi vs cent for vzero A ep_2
        TArrayD*                fChi2C;                         // chi vs cent for vzero C ep_2
        TArrayD*                fChi3A;                         // chi vs cent for vzero A ep_3
        TArrayD*                fChi3C;                         // chi vs cent for vzero C ep_3
        TArrayD*                fSigma2A;                       // chi vs cent for vzero A ep_2
        TArrayD*                fSigma2C;                       // chi vs cent for vzero C ep_2
        TArrayD*                fSigma3A;                       // chi vs cent for vzero A ep_3
        TArrayD*                fSigma3C;                       // chi vs cent for vzero C ep_3
        EPweightType            fWeightForVZERO;                // use chi weight for vzero
        TFile*                  fOADB;                          //! fOADB
        TH2F*                   fHistQxV0aBC;                   //! qx v0a before cuts
        TH2F*                   fHistQyV0aBC;                   //! qx v0a before cuts
        TH2F*                   fHistQxV0cBC;                   //! qx v0a before cuts
        TH2F*                   fHistQyV0cBC;                   //! qx v0a before cuts
        TH2F*                   fHistQxV0a;                     //! qx v0a before cuts
        TH2F*                   fHistQyV0a;                     //! qx v0a before cuts
        TH2F*                   fHistQxV0c;                     //! qx v0a before cuts
        TH2F*                   fHistQyV0c;                     //! qx v0a before cuts
        TH2F*                   fHistMultVsCellBC;              //! fHistMultVsCellBC
        TH1F*                   fHistEPBC;                      //! fHistEPBC
        TH1F*                   fHistEP;                        //! fHistEP
        TH2F*                   fHistMultVsCell;                //! fHistMultVsCell

        AliAnalysisTaskJetV3(const AliAnalysisTaskJetV3&);                  // not implemented
        AliAnalysisTaskJetV3& operator=(const AliAnalysisTaskJetV3&);       // not implemented

        ClassDef(AliAnalysisTaskJetV3, 1);
};

#endif
