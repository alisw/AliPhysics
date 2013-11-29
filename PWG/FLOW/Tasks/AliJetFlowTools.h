/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. * 
 * See cxx source for full Copyright notice */ 
 /* $Id$ */

/* Author: Redmer Alexander Bertens, rbertens@cern.ch, rbertens@nikhef.nl, r.a.bertens@uu.nl
 * see implementation for additional information */

#ifndef AliJetFlowTools_H
#define AliJetFlowTools_H

// root forward declaration
class TF1;
class TH1D;
class TH2D;
class TCanvas;
class TLegend;
class TString;
class TArrayD;
class TGraphErrors;
class TObjArray;
// aliroot forward declarations
class AliAnaChargedJetResponseMaker;
class AliUnfolding;
// root includes
#include "TMatrixD.h"
#include "TList.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TProfile.h"
#include "TVirtualPad.h"
//_____________________________________________________________________________
class AliJetFlowTools {
    public: 
        AliJetFlowTools();
    protected:
        ~AliJetFlowTools();     // no implemented destructor, but class is not compiled
    public:
        // enumerators
        enum unfoldingAlgorithm {
            kChi2,
            kBayesian,
            kSVD,
            kSVDlegacy,         // first implementation
            kNone };            // type of unfolding algorithm
        enum prior {
            kPriorChi2,
            kPriorMeasured };   // used prior for svd unfolding
        // setters, interface to the class
        void            SetSaveFull(Bool_t b)           {fSaveFull              = b;}
        void            SetInputList(TList* list)       {
            fInputList          = list;
            fRefreshInput       = kTRUE;
        }
        void            SetOutputFileName(TString name) {fOutputFileName        = name;}
        void            CreateOutputList(TString name) {
            // create a new output list and add it to the full output
            if(!fOutputFile) fOutputFile = new TFile(fOutputFileName.Data(), "RECREATE");
            fOutputFile->cd();  // avoid nested dirs
            if(name.EqualTo(fActiveString)) {
                printf(Form(" > Warning, duplicate output list, renaming %s to %s_2 ! < \n", name.Data(), name.Data()));
                name+="_2";
            }
            fActiveString = name;
            fActiveDir = new TDirectoryFile(fActiveString.Data(), fActiveString.Data());
            fActiveDir->cd();
        }
        void            SetCentralityBin(Int_t bin)     {fCentralityBin         = bin;}
        void            SetDetectorResponse(TH2D* dr)   {fDetectorResponse      = dr;}
        void            SetBinsTrue(TArrayD* bins)      {fBinsTrue              = bins;}
        void            SetBinsRec(TArrayD* bins)       {fBinsRec               = bins;}
        void            SetBinsTruePrior(TArrayD* bins) {fBinsTruePrior         = bins;}
        void            SetBinsRecPrior(TArrayD* bins)  {fBinsRecPrior          = bins;}
        void            SetSVDReg(Int_t r)              {fSVDRegIn = r; fSVDRegOut = r;}
        void            SetSVDReg(Int_t in, Int_t out)  {fSVDRegIn = in; fSVDRegOut = out;}
        void            SetSVDToy(Bool_t b, Float_t r)  {fSVDToy = b; fJetRadius = r;}
        void            SetBeta(Double_t b)             {fBetaIn = b; fBetaOut = b;}
        void            SetBeta(Double_t i, Double_t o) {fBetaIn = i; fBetaOut = o;}
        void            SetAvoidRoundingError(Bool_t r) {fAvoidRoundingError    = r;}
        void            SetUnfoldingAlgorithm(unfoldingAlgorithm ua)    {fUnfoldingAlgorithm    = ua;}
        void            SetPrior(prior p)                               {fPrior                 = p;}
        void            SetNormalizeSpectra(Bool_t b)   {fNormalizeSpectra      = b;}
        void            SetNormalizeSpectra(Int_t e)    { // use to normalize to this no of events
            fEventCount = e;
            fNormalizeSpectra = kFALSE;
        }
        void            SetSmoothenSpectrum(Bool_t b, Float_t min = 50., Float_t max = 100., Float_t start= 75.) {
            fSmoothenSpectrum   = b;
            fFitMin             = min;
            fFitMax             = max;
            fFitStart           = start;
        }
        void            SetTestMode(Bool_t t)           {fTestMode                      = t;}
        void            SetNoDphi(Bool_t n)             {fNoDphi                        = n;}
        void            SetEventPlaneResolution(Double_t r)     {fEventPlaneRes         = r;}
        void            SetUseDetectorResponse(Bool_t r)        {fUseDetectorResponse   = r;}
        void            SetTrainPowerFit(Bool_t t)      {fTrainPower                    = t;}
        void            Make();
        void            Finish() {
            fOutputFile->cd();
            fRMSSpectrumIn->Write();
            fRMSSpectrumOut->Write();
            fRMSRatio->Write();
            fOutputFile->Close();}
        void            PostProcess(TString def, TString in = "UnfoldedSpectra.root", TString out = "ProcessedSpectra.root");
        Bool_t          SetRawInput (
                TH2D* detectorResponse, // detector response matrix
                TH1D* jetPtIn,          // in plane jet spectrum
                TH1D* jetPtOut,         // out of plane jet spectrum
                TH1D* dptIn,            // in plane delta pt distribution
                TH1D* dptOut,           // out of plane delta pt distribution
                Int_t eventCount = 0);  // event count (optional)
        // static const helper functions, mainly histogram manipulation
        static TH1D*    ResizeXaxisTH1D(TH1D* histo, Int_t low, Int_t up, TString suffix = "");
        static TH2D*    ResizeYaxisTH2D(TH2D* histo, TArrayD* x, TArrayD* y, TString suffix = "");
        static TH2D*    NormalizeTH2D(TH2D* histo);
        static TH1D*    GetUnfoldingTemplate(TH1D* histo, TArrayD* bins, TString suffix = "");
        TH2D*           RebinTH2D(TH2D* histo, TArrayD* binsTrue, TArrayD* binsRec, TString suffix = "");
        static TH2D*    MatrixMultiplication(TH2D* a, TH2D* b, TString name = "CombinedResponse");
        static TH1D*    NormalizeTH1D(TH1D* histo, Double_t scale = 1.);
        static TGraphErrors*    GetRatio(TH1 *h1 = 0x0, TH1* h2 = 0x0, TString name = "", Bool_t appendFit = kFALSE, Int_t xmax = -1);
        static TGraphErrors*    GetV2(TH1* h1 = 0x0, TH1* h2 = 0x0, Double_t r = .63, TString name = "");
        static void     WriteObject(TObject* object);
        static TH2D*    ConstructDPtResponseFromTH1D(TH1D* dpt, Bool_t AvoidRoundingError);
        static TH2D*    GetUnityResponse(TArrayD* binsTrue, TArrayD* binsRec, TString suffix = "");
        void            SaveConfiguration(Bool_t convergedIn, Bool_t convergedOut);
        static TMatrixD*        CalculatePearsonCoefficients(TMatrixD* covmat);
        static TH1D*    SmoothenSpectrum(TH1D* spectrum, TF1* function, Double_t min, Double_t max, Double_t start, Bool_t kill = kTRUE, Bool_t counts = kTRUE);
        // set styles
        static void     Style(TCanvas* c, TString style = "PEARSON");
        static void     Style(TVirtualPad* c, TString style = "SPECTRUM");
        static TLegend* AddLegend(TVirtualPad* p)       {return p->BuildLegend();}
        // interface to AliUnfolding, not necessary but nice to have all parameters in one place
        static void     SetMinuitStepSize(Float_t s)    {AliUnfolding::SetMinuitStepSize(s);}
        static void     SetMinuitPrecision(Float_t s)   {AliUnfolding::SetMinuitPrecision(s);}
        static void     SetMinuitPrecision(Int_t i)     {AliUnfolding::SetMinuitMaxIterations(i);}
        static void     SetMinuitStrategy(Double_t s)   {AliUnfolding::SetMinuitStrategy(s);}
        static void     SetDebug(Int_t d)               {AliUnfolding::SetDebug(d);}
    private:
        Bool_t          PrepareForUnfolding(); 
        void            GetPrior()                      {; /* FIXME implement */};
        Bool_t          UnfoldSpectrumChi2(             TH1D* resizedJetPt, 
                                                        TH2D* resizedResonse,
                                                        TH1D* kinematicEfficiency,
                                                        TH1D* unfoldingTemplate,
                                                        TH1D *&unfolded,        // careful, pointer reference
                                                        TString suffix);
        Bool_t          UnfoldSpectrumBayesian()        {return kFALSE;}
        Bool_t          UnfoldSpectrumSVDlegacy(        TH1D* resizedJetPt, 
                                                        TH2D* resizedResonse,
                                                        TH1D* kinematicEfficiency,
                                                        TH1D* unfoldingTemplate,
                                                        TH1D *&unfolded,        // careful, pointer reference
                                                        TString suffix);
        Bool_t          UnfoldSpectrumSVD(              TH1D* resizedJetPt, 
                                                        TH2D* resizedResonse,
                                                        TH1D* kinematicEfficiency,
                                                        TH1D* unfoldingTemplate,
                                                        TH1D *&unfolded,        // careful, pointer reference
                                                        TString suffix);
        static void     ResetAliUnfolding();
        // give object a unique name via the 'protect heap' functions. 
        // may seem redundant, but some internal functions of root (e.g.
        // ProjectionY()) check for existing objects by name and re-use them
        TH1D*           ProtectHeap(TH1D* protect, Bool_t kill = kTRUE, TString suffix = "");
        TH2D*           ProtectHeap(TH2D* protect, Bool_t kill = kTRUE, TString suffix = "");
        TGraphErrors*   ProtectHeap(TGraphErrors* protect, Bool_t kill = kTRUE, TString suffix = "");
        // members, accessible via setters
        AliAnaChargedJetResponseMaker*  fResponseMaker; // utility object
        TF1*                    fPower;                 // smoothening fit
        Bool_t                  fSaveFull;              // save all generated histograms to file
        TString                 fActiveString;          // identifier of active output
        TDirectoryFile*         fActiveDir;             // active directory
        TList*                  fInputList;             // input list
        Bool_t                  fRefreshInput;          // re-read the input (called automatically if input list changes)
        TString                 fOutputFileName;        // output file name
        TFile*                  fOutputFile;            // output file
        Int_t                   fCentralityBin;         // centrality bin
        TH2D*                   fDetectorResponse;      // detector response
        Double_t                fBetaIn;                // regularization strength, in plane unfolding
        Double_t                fBetaOut;               // regularization strength, out of plane unfoldign
        Bool_t                  fAvoidRoundingError;    // set dpt to zero for small values far from the diagonal
        unfoldingAlgorithm      fUnfoldingAlgorithm;    // algorithm used for unfolding
        prior                   fPrior;                 // prior for unfolding
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
        Bool_t                  fSmoothenSpectrum;      // smoothen the tail of the measured spectrum using a powerlaw fit
        Float_t                 fFitMin;                // lower bound of smoothening fit
        Float_t                 fFitMax;                // upper bound of smoothening fit
        Float_t                 fFitStart;              // from this value, use smoothening
        Bool_t                  fTestMode;              // unfold with unity response for testing
        Bool_t                  fNoDphi;                // no dphi dependence in unfolding
        Bool_t                  fRawInputProvided;      // input histograms provided, not read from file
        Double_t                fEventPlaneRes;         // event plane resolution for current centrality
        Bool_t                  fUseDetectorResponse;   // add detector response to unfolding
        Bool_t                  fTrainPower;            // don't clear the params of fPower for call to Make
                                                        // might give more stable results, but also introduces
                                                        // a bias / dependency on previous iterations
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
        TH1D*                   fUnfoldedIn;            // unfolded in plane spectrum
        TH1D*                   fUnfoldedOut;           // unfolded out ofplane spectrum
        // copy and assignment 
        AliJetFlowTools(const AliJetFlowTools&);             // not implemented
        AliJetFlowTools& operator=(const AliJetFlowTools&);  // not implemented
};
#endif
//_____________________________________________________________________________
