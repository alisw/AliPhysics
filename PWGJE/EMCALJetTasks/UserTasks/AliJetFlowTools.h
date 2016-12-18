/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. * 
 * See cxx source for full Copyright notice */ 
 /* $Id$ */

/* Author: Redmer Alexander Bertens, rbertens@cern.ch, rbertens@nikhef.nl, r.a.bertens@uu.nl
 * see implementation for additional information */

#ifndef AliJetFlowTools_H
#define AliJetFlowTools_H

// root forward declarations
class TH1;
class TF2;
class TH2D;
class TCanvas;
class TString;
class TArrayD;
class TGraph;
class TGraphErrors;
class TObjArray;
// aliroot forward declarations
class AliAnaChargedJetResponseMaker;
class AliUnfolding;
// root includes
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TList.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TProfile.h"
#include "TVirtualPad.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"
// define the following variable to build with debug flags
// #define ALIJETFLOWTOOLS_DEBUG_FLAG

//_____________________________________________________________________________
class AliJetFlowTools {
    public:
        AliJetFlowTools();
    protected:
        ~AliJetFlowTools();     // not implemented (deliberately). object ownership is a bit messy in this class
                                // since most (or all) of the objects are owned by the input and output files
    public:
        // enumerators
        enum unfoldingAlgorithm {       // type of unfolding alrogithm
            kChi2,                      // chi^2 unfolding, implemented in AliUnfolding
            kBayesian,                  // Bayesian unfolding, implemented in RooUnfold
            kBayesianAli,               // Bayesian unfolding, implemented in AliUnfolding
            kSVD,                       // SVD unfolding, implemented in RooUnfold
            kFold,                      // instead of unfolding, fold the input with the response
            kNone };                    // no unfolding
        enum prior {                    // prior that is used for unfolding
            kPriorChi2,                 // prior from chi^2 method
            kPriorMeasured,             // use measured spectrum as prior
            kPriorPythia,               // use pythia spectrum as prior
            kPriorTF1 };                // use properly binned TF1 as prior
        enum histoType {                // histogram identifier, only used internally
            kInPlaneSpectrum,           // default style for spectrum
            kOutPlaneSpectrum,
            kUnfoldedSpectrum,
            kFoldedSpectrum,
            kMeasuredSpectrum,
            kBar,                       // default style for bar histogram
            kRatio,                     // default style for ratio
            kV2,                        // default style for v2
            kDeltaPhi,                  // default for delta phi
            kEmpty };                   // default style
        // setters, interface to the class
        void            SetOffsetStart(Int_t g)         {gOffsetStart           = g;}
        void            SetOffsetStop(Int_t g)          {gOffsetStop            = g;}
        void            SetReductionFactor(Float_t g)   {gReductionFactor       = g;}
        void            SetReductionFactorCorr(Float_t g)       {gReductionFactorCorr   =g;}
        void            SetPwrtTo(Float_t p)            {gPwrtTo                = p;}
        void            SetPwrtToArray(TArrayD* a)      {
            gPwrtToArray = a;   // set the array
            gPwrtTo = -9999;    // tells instance to not use this value
        }
        void            SetPwrtToStatArray(TArrayD* a)  {gPwrtToStatArray       = a;}
        void            SetPivot(Float_t p)             {fPivot                 = p;}
        void            SetConstantUE(Bool_t ue)        {fConstantUE            = ue;}
        void            SetSubdueError(Bool_t b)        {fSubdueError           = b;}
        void            SetSaveFull(Bool_t b)           {fSaveFull              = b;}
        void            SetInputList(TList* list)       {
            fInputList          = list;
            fRefreshInput       = kTRUE;
        }
        void            SetOutputFileName(TString name) {fOutputFileName        = name;}
        void            CreateOutputList(TString name) {
            // create a new output list and add it to the full output
            fListPrefix++;
            if(!fOutputFile) fOutputFile = new TFile(fOutputFileName.Data(), "RECREATE");
            fOutputFile->cd();  // avoid nested dirs
            if(name.EqualTo(fActiveString)) {
                printf(Form(" > Warning, duplicate output list, renaming %s to %s_2 ! < \n", name.Data(), name.Data()));
                name+="_2";
            }
            fActiveString = Form("%i__%s", fListPrefix, name.Data());
            fActiveDir = new TDirectoryFile(fActiveString.Data(), fActiveString.Data());
            fActiveDir->cd();
        }
        void            SetCentralityBin(Int_t bin)             {
            // in case of one centraltiy
            fCentralityArray = new TArrayI(1);
            fCentralityArray->AddAt(bin, 0);
            // for one centrality there's no need for weights
            fCentralityWeights = new TArrayD(1);
            fCentralityWeights->AddAt(1., 0);
        }
        void            SetMergeSpectrumBins(TArrayI* a)        {fMergeBinsArray        =       a;}
        void            SetCentralityBin(TArrayI* bins)         {
            fCentralityArray = bins;
        }
        void            SetCentralityWeight(TArrayD* weights)   {
            fCentralityWeights = weights;
            if(!fCentralityArray) printf(" > Warning: centrality weights set, but bins are not defined! \n");
        }
        void            SetMergeWith(
                TList* l, 
                Int_t c,
                Float_t weight) {
            fMergeWithList      = l; 
            fMergeWithCen       = c;
            fMergeWithWeight    = weight;
        }
        void            SetDetectorResponse(TH2D* dr)           {fDetectorResponse      = dr;}
        void            SetJetFindingEfficiency(TH1D* e)        {fJetFindingEff         = e;}
        void            SetBinsTrue(TArrayD* bins)              {fBinsTrue              = bins;}
        void            SetBinsRec(TArrayD* bins)               {fBinsRec               = bins;}
        void            SetBinsTruePrior(TArrayD* bins)         {fBinsTruePrior         = bins;}
        void            SetBinsRecPrior(TArrayD* bins)          {fBinsRecPrior          = bins;}
        void            SetSVDReg(Int_t r)                      {fSVDRegIn              = r; fSVDRegOut         = r;}
        void            SetSVDReg(Int_t in, Int_t out)          {fSVDRegIn              = in; fSVDRegOut        = out;}
        void            SetSVDToy(Bool_t b, Float_t r)          {fSVDToy                = b; fJetRadius         = r;}
        void            SetBeta(Double_t b)                     {fBetaIn                = b; fBetaOut           = b;}
        void            SetBeta(Double_t i, Double_t o)         {fBetaIn                = i; fBetaOut           = o;}
        void            SetBayesianIter(Int_t i)                {fBayesianIterIn        = i; fBayesianIterOut    = i;}
        void            SetBayesianIter(Int_t i, Int_t o)       {fBayesianIterIn        = i; fBayesianIterOut    = o;}
        void            SetBayesianSmooth(Float_t s)            {fBayesianSmoothIn      = s; fBayesianSmoothOut = s;}
        void            SetBayesianSmooth(Float_t i, Float_t o) {fBayesianSmoothIn      = i; fBayesianSmoothOut = o;}
        void            SetAvoidRoundingError(Bool_t r)         {fAvoidRoundingError    = r;}
        void            SetUnfoldingAlgorithm(unfoldingAlgorithm ua)    {fUnfoldingAlgorithm                    = ua;}
        void            SetPrior(prior p)                       {fPrior                 = p;}
        void            SetPrior(prior p, TF1* function, TArrayD* bins) {
            fPrior = p;
            // set prior to value supplied in TF1
            fPriorUser = new TH1D("prior_user", "prior_user", bins->GetSize()-1, bins->GetArray());
            // loop over bins and fill the histo from the tf1
            for(Int_t i(0); i < fPriorUser->GetNbinsX() + 1; i++) fPriorUser->SetBinContent(i, function->Integral(fPriorUser->GetXaxis()->GetBinLowEdge(i), fPriorUser->GetXaxis()->GetBinUpEdge(i))/fPriorUser->GetXaxis()->GetBinWidth(i));
        }
        void            SetPrior(prior p, TH1D* spectrum)       {fPrior                 = p; fPriorUser         = spectrum;}
        void            SetNormalizeSpectra(Bool_t b)           {fNormalizeSpectra      = b;}
        void            SetNormalizeSpectra(Int_t e)            { // use to normalize to this no of events
            fEventCount         = e;
            fNormalizeSpectra   = kFALSE;
        }
        void            SetSmoothenPrior(Bool_t b, Float_t min = 50., Float_t max = 100., Float_t start= 75., Bool_t counts = kTRUE) {
            fSmoothenPrior      = b;
            fFitMin             = min;
            fFitMax             = max;
            fFitStart           = start;
            fSmoothenCounts     = counts;
        }
        void            SetTestMode(Bool_t t)                   {fTestMode              = t;}
        void            SetEventPlaneResolution(Double_t r)     {fEventPlaneRes         = r;}
        void            SetUseDetectorResponse(Bool_t r)        {fUseDetectorResponse   = r;}
        void            SetUseDptResponse(Bool_t r)             {fUseDptResponse        = r;}
        void            SetTrainPowerFit(Bool_t t)              {fTrainPower            = t;}
        void            SetDphiUnfolding(Bool_t i)              {fDphiUnfolding         = i;}
        void            SetDphiDptUnfolding(Bool_t i)           {fDphiDptUnfolding      = i;}
        void            SetExLJDpt(Bool_t i)                    {fExLJDpt               = i;}
        void            SetWeightFunction(TF1* w)               {fResponseMaker->SetRMMergeWeightFunction(w);}
        void            SetRMS(Bool_t r)                        {fRMS                   = r;}
        void            SetSymmRMS(Bool_t r)                    {fSymmRMS               = r;}
        void            SetRho0(Bool_t r)                       {fRho0                  = r;}
        void            SetBootstrap(Bool_t b, Bool_t r = kTRUE) {
            // note that setting this option will not lead to true resampling
            // but rather to randomly drawing a new distribution from a pdf
            // of the measured distribution
            fBootstrap             = b;
            // by default fully randomize randomizer from system time
            if(r) {
                delete gRandom;
                gRandom = new TRandom3(0);
            }
        }
        void            SetHarmonic(Int_t n)                    {fHarmonic              = n;}
        // main function. buffers about 5mb per call!
        void            Make(TH1* customIn = 0x0, TH1* customOut = 0x0);
        void            MakeAU();       // test function, use with caution (09012014)
        void            Finish() {
            fOutputFile->cd();
            if(fRMSSpectrumIn)  fRMSSpectrumIn->Write();
            if(fRMSSpectrumOut) fRMSSpectrumOut->Write();
            if(fRMSRatio)       fRMSRatio->Write();
            fOutputFile->Close();}
        void            PostProcess(
                TString def,
                Int_t columns = 4,
                Float_t rangeLow = 20,
                Float_t rangeUp = 80,
                TString in = "UnfoldedSpectra.root", 
                TString out = "ProcessedSpectra.root") const;
        void            BootstrapSpectra(
                TString def,
                TString in = "UnfoldedSpectra.root", 
                TString out = "BootstrapSpectra.root") const;
        void            GetNominalValues(
                TH1D*& ratio,
                TGraphErrors*& v2,
                TArrayI* in,
                TArrayI* out,
                TString inFile = "UnfoldedSpectra.root",
                TString outFile = "Nominal.root") const;
        void            GetCorrelatedUncertainty(
                TGraphAsymmErrors*& corrRatio,
                TGraphAsymmErrors*& corrV2,
                TArrayI* variationsIn,
                TArrayI* variationsOut,
                Bool_t sym,
                TArrayI* variantions2ndIn,
                TArrayI* variantions2ndOut,
                Bool_t sym2nd,
                TString type = "",
                TString type2 = "",
                Int_t columns = 4,
                Float_t rangeLow = 20,
                Float_t rangeUp = 80,
                Float_t corr = .5,
                TString in = "UnfoldedSpectra.root", 
                TString out = "CorrelatedUncertainty.root") const;
        void            GetShapeUncertainty(
                TGraphAsymmErrors*& shapeRatio,
                TGraphAsymmErrors*& shapeV2,
                TArrayI* regularizationIn,
                TArrayI* regularizationOut,
                TArrayI* recBinIn = 0x0,
                TArrayI* recBinOut = 0x0,
                TArrayI* methodIn = 0x0,
                TArrayI* methodOut = 0x0,
                Int_t columns = 4,
                Float_t rangeLow = 20,
                Float_t rangeUp = 80,
                Float_t corr = .0,
                TString in = "UnfoldedSpectra.root", 
                TString out = "ShapeUncertainty.root",
                Bool_t regularizationOnV2 = kTRUE     // get uncertainty on yields separately or v2 directly
                ) const;
        Bool_t          SetRawInput (
                TH2D* detectorResponse, // detector response matrix
                TH1D* jetPtIn,          // in plane jet spectrum
                TH1D* jetPtOut,         // out of plane jet spectrum
                TH1D* dptIn,            // in plane delta pt distribution
                TH1D* dptOut,           // out of plane delta pt distribution
                Int_t eventCount = 0);  // event count (optional)
        TH1*            GetUnfoldedSpectrumIn() const   {return fUnfoldedSpectrumIn;}
        TH1*            GetUnfoldedSpectrumOut() const  {return fUnfoldedSpectrumOut;}
        // static const helper functions, mainly histogram manipulation
        static void     RemoveSign(Double_t& d)         {d = TMath::Abs(d);}
        static TH1D*    ResizeXaxisTH1D(TH1D* histo, Int_t low, Int_t up, TString suffix = "");
        static TH2D*    ResizeYaxisTH2D(TH2D* histo, TArrayD* x, TArrayD* y, TString suffix = "");
        static TH2D*    NormalizeTH2D(TH2D* histo, Bool_t noError = kTRUE);
        static TH1*     Bootstrap(TH1* hist, Bool_t kill = kTRUE);
        static TH1D*    RebinTH1D(TH1D* histo, TArrayD* bins, TString suffix = "", Bool_t kill = kTRUE);
        TH2D*           RebinTH2D(TH2D* histo, TArrayD* binsTrue, TArrayD* binsRec, TString suffix = "");
        static TH2D*    MatrixMultiplication(TH2D* a, TH2D* b, TString name = "CombinedResponse");
        static TH1D*    NormalizeTH1D(TH1D* histo, Double_t scale = 1.);
        static TH1D*    MergeSpectrumBins(TArrayI* bins, TH1D* spectrum, TH2D* corr);
        static TGraphErrors*    GetRatio(TH1 *h1 = 0x0, TH1* h2 = 0x0, TString name = "", Bool_t appendFit = kFALSE, Int_t xmax = -1);
        static TGraphErrors*    GetV2(TH1* h1 = 0x0, TH1* h2 = 0x0, Double_t r = 0., TString name = "");
        static TH1D*            GetV2Histo(TH1* h1 = 0x0, TH1* h2 = 0x0, Double_t r = 0., TString name = "");
        static TH1F*            ConvertGraphToHistogram(TGraphErrors* g);
        static TGraphAsymmErrors*       AddHistoErrorsToGraphErrors(TGraphAsymmErrors* g, TH1D* h);
        static Double_t         GetRMSOfTH1(TH1* h, Double_t a, Double_t b);
        static TF1*             GetErrorFromFit(TH1* h1, TH1* h2, Double_t a, Double_t b, 
                Float_t pivot = 50., Bool_t subdueError = kFALSE, 
                TString str = "", Bool_t setContent = kTRUE);
        void     ReplaceBins(TArrayI* array, TGraphAsymmErrors* graph);
        void     ReplaceBins(TArrayI* array, TGraphErrors* graph);
        TGraphAsymmErrors*      GetV2WithSystematicErrors(
                TH1* h1, TH1* h2, Double_t r, TString name, 
                TH1* relativeErrorInUp,
                TH1* relativeErrorInLow,
                TH1* relativeErrorOutUp,
                TH1* relativeErrorOutLow,
                Float_t rho = 0.) const;
        static void     GetSignificance(
                TGraphErrors* n,                // points with stat error
                TGraphAsymmErrors* shape,       // points with shape error
                TGraphAsymmErrors* corr,        // corr with stat error
                Int_t low,                      // pt lower level
                Int_t up                        // pt upper level
        );
        static void     MinimizeChi2nd();
        static Double_t PhenixChi2nd(const Double_t *xx );
        static Double_t ConstructFunctionnd(Double_t *x, Double_t *par);
        static TF2*     ReturnFunctionnd(Double_t &p);
        static void     WriteObject(TObject* object, TString suffix = "", Bool_t kill = kTRUE);
        static TH2D*    ConstructDPtResponseFromTH1D(TH1D* dpt, Bool_t AvoidRoundingError);
        static TH2D*    GetUnityResponse(TArrayD* binsTrue, TArrayD* binsRec, TString suffix = "");
        void            SaveConfiguration(Bool_t convergedIn, Bool_t convergedOut) const;
        static TMatrixD*        CalculatePearsonCoefficients(TMatrixD* covmat);
        static TH1D*    SmoothenPrior(TH1D* spectrum, TF1* function, Double_t min, Double_t max, Double_t start, Bool_t kill = kTRUE, Bool_t counts = kTRUE);
        // set style
        void            SetTitleFontSize(Double_t s)    {fTitleFontSize = s;}
        static void     Style(Bool_t legacy = kFALSE);
        static void     Style(TCanvas* c, TString style = "PEARSON");
        static void     Style(TVirtualPad* c, TString style = "SPECTRUM", Bool_t legacy = kFALSE);
        static void     Style(TLegend* l);
        static void     Style(TH1* h, EColor col = kBlue, histoType = kEmpty, Bool_t legacy = kFALSE);
        static void     Style(TGraph* h, EColor col = kBlue, histoType = kEmpty, Bool_t legacy = kFALSE);
        static TLegend* AddLegend(TVirtualPad* p, Bool_t style = kFALSE) {
            if(!style) return p->BuildLegend(.565, .663, .882, .883);
            else {
                TLegend* l = AddLegend(p, kFALSE);
                Style(l);
                return l;
            }
        }
        static TPaveText*       AddTPaveText(
                // this text appears under the logo
                TString text, 
                Int_t r = 2,
                Double_t a = .587,
                Double_t b = .695,
                Double_t c = .872,
                Double_t d = .801) {
            TPaveText* t(new TPaveText(a, b, c, d, "NDC"));
            t->SetFillColor(0);            
            t->SetBorderSize(0);
            t->AddText(0.,0.,text.Data());
            t->AddText(0., 0., Form("#it{R} = 0.%i anti-#it{k}_{T}, |#eta_{jet}|<%.1f", r, .9-r/10.));
            t->SetTextColor(kBlack);
            t->SetTextFont(42);
            t->SetTextSize(gStyle->GetTextSize()*.8);
            t->Draw("same");
            return t;
        } 
        static TPaveText*       AddText(
                TString text, 
                EColor col,
                Double_t a = .2098,
                Double_t b = .5601,
                Double_t c = .613,
                Double_t d = .6211) {
            TPaveText* t(new TPaveText(a, b, c, d, "NDC"));
            t->SetFillColor(0);            
            t->SetBorderSize(0);
            t->AddText(0.,0.,text.Data());
            t->SetTextColor(col);
            t->SetTextFont(42);
            t->SetTextSize(gStyle->GetTextSize()*.8);
            t->Draw("same");
            return t;
        } 
        static TLatex*          AddLogo(Int_t logo, Double_t xmin = .59, Double_t ymax = .81) {
            if(logo == 0) return AddTLatex(xmin, ymax, "ALICE");
            else if (logo == 1) return AddTLatex(xmin, ymax, "ALICE Preliminary");
            else if (logo == 2) return AddTLatex(xmin, ymax, "ALICE Simulation");
            else if (logo == 3) return AddTLatex(xmin, ymax, "work in progress");
            return 0x0;
        }
        static TLatex*          AddSystem() {
            return AddTLatex(0.55, 87, "Pb-Pb #sqrt{#it{s}}}_{NN} = 2.76 TeV");
        }
        static TLatex*          AddTLatex(Double_t xmin, Double_t ymax, TString string) {
TLatex* tex = new TLatex(xmin, ymax, string.Data());
            tex->SetNDC();
            tex->SetTextFont(42);
            tex->Draw("same");
            return tex;
        }

        static void     SavePadToPDF(TVirtualPad* pad)  {
            if(pad) return;/*pad->SaveAs(Form("%s.pdf", pad->GetName()));*/
            else return;}
        // interface to AliUnfolding, not necessary but nice to have all parameters in one place
        static void     SetMinuitStepSize(Float_t s)    {AliUnfolding::SetMinuitStepSize(s);}
        static void     SetMinuitPrecision(Float_t s)   {AliUnfolding::SetMinuitPrecision(s);}
        static void     SetMinuitPrecision(Int_t i)     {AliUnfolding::SetMinuitMaxIterations(i);}
        static void     SetMinuitStrategy(Double_t s)   {AliUnfolding::SetMinuitStrategy(s);}
        static void     SetDebug(Int_t d)               {AliUnfolding::SetDebug(d);}
    private:
        Bool_t          PrepareForUnfolding(TH1* customIn = 0x0, TH1* customOut = 0x0); 
        Bool_t          PrepareForUnfolding(Int_t low, Int_t up);
        TH1D*           GetPrior(                       const TH1D* measuredJetSpectrum,
                                                        const TH2D* resizedResponse,
                                                        const TH1D* kinematicEfficiency,
                                                        const TH1D* measuredJetSpectrumTrueBins,
                                                        const TString suffix,
                                                        const TH1D* jetFindingEfficiency);
        TH1D*           UnfoldWrapper(                  const TH1D* measuredJetSpectrum, 
                                                        const TH2D* resizedResponse,
                                                        const TH1D* kinematicEfficiency,
                                                        const TH1D* measuredJetSpectrumTrueBins,
                                                        const TString suffix,
                                                        const TH1D* jetFindingEfficiency = 0x0);
        TH1D*           UnfoldSpectrumChi2(             const TH1D* measuredJetSpectrum, 
                                                        const TH2D* resizedResponse,
                                                        const TH1D* kinematicEfficiency,
                                                        const TH1D* measuredJetSpectrumTrueBins,
                                                        const TString suffix,
                                                        const TH1D* jetFindingEfficiency = 0x0);
        TH1D*           UnfoldSpectrumSVD(              const TH1D* measuredJetSpectrum, 
                                                        const TH2D* resizedResponse,
                                                        const TH1D* kinematicEfficiency,
                                                        const TH1D* measuredJetSpectrumTrueBins,
                                                        const TString suffix,
                                                        const TH1D* jetFindingEfficiency = 0x0);
        TH1D*           UnfoldSpectrumBayesianAli(      const TH1D* measuredJetSpectrum, 
                                                        const TH2D* resizedResponse,
                                                        const TH1D* kinematicEfficiency,
                                                        const TH1D* measuredJetSpectrumTrueBins,
                                                        const TString suffix,
                                                        const TH1D* jetFindingEfficiency = 0x0);
        TH1D*           UnfoldSpectrumBayesian(         const TH1D* measuredJetSpectrum, 
                                                        const TH2D* resizedResponse,
                                                        const TH1D* kinematicEfficiency,
                                                        const TH1D* measuredJetSpectrumTrueBins,
                                                        const TString suffix,
                                                        const TH1D* jetFindingEfficiency = 0x0);
        TH1D*           FoldSpectrum(                   const TH1D* measuredJetSpectrum, 
                                                        const TH2D* resizedResponse,
                                                        const TH1D* kinematicEfficiency,
                                                        const TH1D* measuredJetSpectrumTrueBins,
                                                        const TString suffix,
                                                        const TH1D* jetFindingEfficiency = 0x0);
       void            SystematicsWrapper(
                TArrayI* variationsIn,
                TArrayI* variationsOut,
                TH1D*& relativeErrorInUp,
                TH1D*& relativeErrorInLow,
                TH1D*& relativeErrorOutUp,
                TH1D*& relativeErrorOutLow,
                TH1D*& relativeSystematicIn,
                TH1D*& relativeSystematicOut,
                TH1D*& nominal,
                TH1D*& nominalIn,
                TH1D*& nominalOut,
                Int_t columns,
                Float_t rangeLow,
                Float_t rangeUp,
                TFile* readMe, 
                TString source = "",
                Bool_t RMS = kFALSE,
                Bool_t onRatio = kTRUE) const;
        void            DoIntermediateSystematics(
                TArrayI* variationsIn,
                TArrayI* variationsOut,
                TH1D*& relativeErrorInUp,
                TH1D*& relativeErrorInLow,
                TH1D*& relativeErrorOutUp,
                TH1D*& relativeErrorOutLow,
                TH1D*& relativeSystematicIn,
                TH1D*& relativeSystematicOut,
                TH1D*& nominal,
                TH1D*& nominalIn,
                TH1D*& nominalOut,
                Int_t columns,
                Float_t rangeLow,
                Float_t rangeUp,
                TFile* readMe, 
                TString source = "",
                Bool_t RMS = kFALSE) const;
        void            DoIntermediateSystematicsOnV2(
                TArrayI* variationsIn,
                TArrayI* variationsOut,
                TH1D*& relativeErrorInUp,
                TH1D*& relativeErrorInLow,
                TH1D*& relativeErrorOutUp,
                TH1D*& relativeErrorOutLow,
                TH1D*& relativeSystematicIn,
                TH1D*& relativeSystematicOut,
                TH1D*& nominal,
                TH1D*& nominalIn,
                TH1D*& nominalOut,
                Int_t columns,
                Float_t rangeLow,
                Float_t rangeUp,
                TFile* readMe, 
                TString source = "",
                Bool_t RMS = kFALSE) const;

        static void     ResetAliUnfolding();
        static void     SquelchWarning() {
            printf(" >> I squelched a warning, jay, I'm contributing ! << \n");
            return;
        }
        // give object a unique name via the 'protect heap' functions. 
        // may seem redundant, but some internal functions of root (e.g.
        // ProjectionY()) check for existing objects by name and re-use them
        TH1D*           ProtectHeap(TH1D* protect, Bool_t kill = kTRUE, TString suffix = "") const;
        TH2D*           ProtectHeap(TH2D* protect, Bool_t kill = kTRUE, TString suffix = "") const;
        TGraphErrors*   ProtectHeap(TGraphErrors* protect, Bool_t kill = kTRUE, TString suffix = "") const;
        // members, accessible via setters
        Int_t                   fListPrefix;            // prefix for list readability
        AliAnaChargedJetResponseMaker*  fResponseMaker; // utility object
        Bool_t                  fRMS;                   // systematic method
        Bool_t                  fSymmRMS;               // symmetric systematic method
        Bool_t                  fConstantUE;            // assign a constant unfolding error
        Bool_t                  fRho0;                  // use the result obtained with the 'classic' fixed rho
        Bool_t                  fBootstrap;             // use bootstrap resampling of input data
        TF1*                    fPower;                 // smoothening fit
        Bool_t                  fSaveFull;              // save all generated histograms to file
        TString                 fActiveString;          // identifier of active output
        TDirectoryFile*         fActiveDir;             // active directory
        TList*                  fInputList;             // input list
        Bool_t                  fRefreshInput;          // re-read the input (called automatically if input list changes)
        TString                 fOutputFileName;        // output file name
        TFile*                  fOutputFile;            // output file
        TArrayI*                fCentralityArray;       // array of centrality bins that are merged
        TArrayI*                fMergeBinsArray;        // array of pt bins that are merged
        TArrayD*                fCentralityWeights;     // array of centrality weights
        TList*                  fMergeWithList;         // additional input list
        Int_t                   fMergeWithCen;          // centrality bin of additional input list
        Float_t                 fMergeWithWeight;       // weight of additional input list
        TH2D*                   fDetectorResponse;      // detector response
        TH1D*                   fJetFindingEff;         // jet finding efficiency
        Double_t                fBetaIn;                // regularization strength, in plane unfolding
        Double_t                fBetaOut;               // regularization strength, out of plane unfoldign
        Int_t                   fBayesianIterIn;        // bayesian regularization parameter, in plane unfolding
        Int_t                   fBayesianIterOut;       // bayesian regularization parameter, out plane unfolding
        Float_t                 fBayesianSmoothIn;      // bayesian smoothening parameter (AliUnfolding)
        Float_t                 fBayesianSmoothOut;     // bayesian smoothening parameter (AliUnfolding)
        Bool_t                  fAvoidRoundingError;    // set dpt to zero for small values far from the diagonal
        unfoldingAlgorithm      fUnfoldingAlgorithm;    // algorithm used for unfolding
        prior                   fPrior;                 // prior for unfolding
        TH1D*                   fPriorUser;             // user supplied prior (e.g. pythia spectrum)
        TArrayD*                fBinsTrue;              // pt true bins
        TArrayD*                fBinsRec;               // pt rec bins
        TArrayD*                fBinsTruePrior;         // holds true bins for the chi2 prior for SVD. setting this is optional
        TArrayD*                fBinsRecPrior;          // holds rec bins for the chi2 prior for SVD. setting this is optional
        Int_t                   fSVDRegIn;              // svd regularization (a good starting point is half of the number of bins)
        Int_t                   fSVDRegOut;             // svd regularization out of plane
        Bool_t                  fSVDToy;                // use toy to estimate coveriance matrix for SVD method
        Float_t                 fJetRadius;             // jet radius (for SVD toy)
        Int_t                   fEventCount;            // number of events
        Bool_t                  fNormalizeSpectra;      // normalize spectra to event count
        Bool_t                  fSmoothenPrior;         // smoothen the tail of the measured spectrum using a powerlaw fit
        Float_t                 fFitMin;                // lower bound of smoothening fit
        Float_t                 fFitMax;                // upper bound of smoothening fit
        Float_t                 fFitStart;              // from this value, use smoothening
        Bool_t                  fSmoothenCounts;        // fill smoothened spectrum with counts
        Bool_t                  fTestMode;              // unfold with unity response for testing
        Bool_t                  fRawInputProvided;      // input histograms provided, not read from file
        Double_t                fEventPlaneRes;         // event plane resolution for current centrality
        Bool_t                  fUseDetectorResponse;   // add detector response to unfolding
        Bool_t                  fUseDptResponse;        // add dpt response to unfolding
        Bool_t                  fTrainPower;            // don't clear the params of fPower for call to Make
                                                        // might give more stable results, but possibly introduces
                                                        // a bias / dependency on previous iterations
        Bool_t                  fDphiUnfolding;         // do the unfolding in in and out of plane orientation
        Bool_t                  fDphiDptUnfolding;      // do the unfolding in dphi and dpt bins (to fit v2)
        Bool_t                  fExLJDpt;               // exclude randon cones with leading jet
        Double_t                fTitleFontSize;         // title font size
        // members, set internally
        TProfile*               fRMSSpectrumIn;         // rms of in plane spectra of converged unfoldings
        TProfile*               fRMSSpectrumOut;        // rms of out of plane spectra of converged unfoldings
        TProfile*               fRMSRatio;              // rms of ratio of converged unfolded spectra
        TProfile*               fRMSV2;                 // rms of v2 of converged unfolded spectra
        TH2D*                   fDeltaPtDeltaPhi;       // delta pt delta phi distribution
        TH2D*                   fJetPtDeltaPhi;         // jet pt delta phi distribution
        TH1D*                   fSpectrumIn;            // in plane jet pt spectrum
        TH1D*                   fSpectrumOut;           // out of plane jet pt spectrum
        TH1D*                   fDptInDist;             // in plane dpt distribution
        TH1D*                   fDptOutDist;            // out of plane dpt distribution
        TH2D*                   fDptIn;                 // in plane dpt matrix
        TH2D*                   fDptOut;                // out plane dpt matrix
        TH2D*                   fFullResponseIn;        // full response matrix, in plane
        TH2D*                   fFullResponseOut;       // full response matrix, out of plane
        Float_t                 fPivot;                 // pivot of step-function
        Bool_t                  fSubdueError;           // no error > pivot
        TH1*                    fUnfoldedSpectrumIn;    // unfolded spectrum in plane
        TH1*                    fUnfoldedSpectrumOut;   // unfolded spectrum out of plane
        Int_t                   fHarmonic;              // vn harmonic

        static TArrayD*         gV2;                    // internal use only, do not touch these
        static TArrayD*         gStat;                  // internal use only, do not touch these
        static TArrayD*         gShape;                 // internal use only, do not touch these
        static TArrayD*         gCorr;                  // internal use only, do not touch these
        
        static Int_t            gOffsetStart;           // see initialization below
        static Int_t            gOffsetStop;            // see initialization below
        static Float_t          gReductionFactor;       // multiply shape uncertainty by this factor
        static Float_t          gReductionFactorCorr;   // multiply corr uncertainty by this factor
        static Float_t          gPwrtTo;                // p-value will be evaluated wrt y = gPwrtTo
        static TArrayD*         gPwrtToArray;           // p-value will be evaluated wrt y = gPwrtToArray[i]
        static TArrayD*         gPwrtToStatArray;       // stat uncertainty of previous array

        // copy and assignment 
        AliJetFlowTools(const AliJetFlowTools&);             // not implemented
        AliJetFlowTools& operator=(const AliJetFlowTools&);  // not implemented

};
// initialize the static members
TArrayD* AliJetFlowTools::gV2           = new TArrayD(7);       // DO NOT TOUCH - these arrays are filled by the
TArrayD* AliJetFlowTools::gStat         = new TArrayD(7);       // 'GetSignificance' function and
TArrayD* AliJetFlowTools::gShape        = new TArrayD(7);       // then used in the chi2 minimization routine
TArrayD* AliJetFlowTools::gCorr         = new TArrayD(7);       // to calculate the significance of the results
Int_t    AliJetFlowTools::gOffsetStart  =  0;           // start chi2 fit from this bin w.r.t. the binning supplied in the 'GetCorr/GetShape' functions
Int_t    AliJetFlowTools::gOffsetStop   = 0;           // stop chi2 fit at this bin w.r.t. the binning supplied in the 'GetCorr/GetShape' functions
Float_t  AliJetFlowTools::gReductionFactor      = 1.;   // multiply shape uncertainty by this factor
Float_t  AliJetFlowTools::gReductionFactorCorr  = 1.;   // multiply corr uncertainty by this factor
Float_t  AliJetFlowTools::gPwrtTo               = 0.;   // p-value will be evaluated wrt y = gPwrtTo
TArrayD* AliJetFlowTools::gPwrtToArray          = new TArrayD(7);// array of values w.r.t. p value is evaluated
TArrayD* AliJetFlowTools::gPwrtToStatArray      = new TArrayD(7);// array of stat errors on values w.r.t. p values is evaluated
#endif
//_____________________________________________________________________________
