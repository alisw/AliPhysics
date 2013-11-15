/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. * 
 * See cxx source for full Copyright notice */ 
 /* $Id$ */

/* Author: Redmer Alexander Bertens, rbertens@cern.ch, rbertens@nikhef.nl, r.a.bertens@uu.nl
 * see implementation for additional information */

#ifndef AliJetFlowTools_H
#define AliJetFlowTools_H

// forward declarations
class TH1D;
class TH2D;
class TString;
class TArrayD;
class TGraphErrors;
class TObjArray;

// includes
#include "TMatrixD.h"
#include "TList.h"
#include "TDirectoryFile.h"
#include "TFile.h"

//_____________________________________________________________________________
class AliJetFlowTools {
    public: 
        AliJetFlowTools();
    protected:
        ~AliJetFlowTools();
    public:
        // enumerators
        enum unfoldingAlgorithm {
            kChi2,
            kBayesian,
            kSVD };             // type of unfolding algorithm
        enum prior {
            kPriorChi2,
            kPriorMeasured };   // used prior for svd unfolding

        // setters, setup the procedure
        void            SetInputList(TList* list)       {fInputList             = list;}
        void            SetOutputFileName(TString name) {fOutputFileName        = name;}
        void            CreateOutputList(TString name) {
            // create a new output list and add it to the full output
            if(!fOutputFile) fOutputFile = new TFile(fOutputFileName.Data(), "RECREATE");
            fOutputFile->cd();  // avoid nested dirs
            fActiveString = name;
            fActiveDir = new TDirectoryFile(fActiveString.Data(), fActiveString.Data());
            fActiveDir->cd();
        }
        void            SetCentralityBin(Int_t bin)     {fCentralityBin         = bin;}
        void            SetdPhidPpt(Bool_t inout)       {fdPhidPt               = inout;}
        void            SetDetectorResponse(TH2D* dr)   {fDetectorResponse      = dr;}
        void            SetBinsTrue(TArrayD* bins)      {fBinsTrue              = bins;}
        void            SetBinsRec(TArrayD* bins)       {fBinsRec               = bins;}
        void            SetSVDDraw(Int_t draw)          {fSVDDraw               = draw;}
        void            SetSVDMinMaxIter(Int_t min, Int_t max)          {fIterMin = min; fIterMax = max;}
        void            SetSVDToy(Bool_t b, Float_t r)  {fSVDToy = b; fJetRadius = r;}
        void            SetBeta(Double_t beta)          {fBeta                  = beta;}
        void            SetBetaDOF(Double_t betaDOF)    {fBetaPerDOF            = betaDOF;}
        void            SetAvoidRoundingError(Bool_t r) {fAvoidRoundingError    = r;}
        void            SetUnfoldingAlgorithm(unfoldingAlgorithm ua)    {fUnfoldingAlgorithm    = ua;}
        void            SetPrior(prior p)                               {fPrior                 = p;}
        void            Make();
        void            Finish()                        {fOutputFile->Close();}
        // static const helper functions, mainly histogram manipulation
        static TH1D*    ResizeXaxisTH1D(TH1D* histo, Int_t low, Int_t up, TString suffix = "");
        static TH2D*    ResizeYaxisTH2D(TH2D* histo, TArrayD* x, TArrayD* y, TString suffix = "");
        static TH2D*    NormalizeTH2D(TH2D* histo);
        static TH1D*    GetUnfoldingTemplate(TH1D* histo, TArrayD* bins, TString suffix = "");
        static TH2D*    RebinTH2DX(TH2D* histo, TArrayD* bins, TString suffix = "");
        static TH2D*    RebinTH2DY(TH2D* histo, TArrayD* bins);
        static TH2D*    MatrixMultiplicationTH2D(TH2D* A, TH2D* B, TString name = "CombinedResponse");
        static TH1D*    NormalizeTH1D(TH1D* histo, Double_t scale = 1.);
        static TGraphErrors*    GetRatio(TH1 *h1 = 0x0, TH1* h2 = 0x0, TString name = "", Bool_t appendFit = kFALSE, Double_t low = 0., Double_t up = 100.);
    private:
        Bool_t          PrepareForUnfolding(); 
        Bool_t          UnfoldSpectrumChi2(     TH1D* resizedJetPt, 
                                                TH2D* resizedResonse,
                                                TH1D* kinematicEfficiency,
                                                TH1D* unfoldingTemplate,
                                                TH1D *&unfolded,        // careful, pointer reference
                                                TString suffix);
        Bool_t          UnfoldSpectrumBayesian()        {return kFALSE;}
        Bool_t          UnfoldSpectrumSVD(      TH1D* resizedJetPt, 
                                                TH2D* resizedResonse,
                                                TH1D* kinematicEfficiency,
                                                TH1D* unfoldingTemplate,
                                                TH1D *&unfolded,        // careful, pointer reference
                                                TString suffix);
        TMatrixD*       CalculatePearsonCoefficients(TMatrixD* covmat);
        // members, accessible via setters
        TString                 fActiveString;          // identifier of active output
        TDirectoryFile*         fActiveDir;             // active directory
        TList*                  fInputList;             // input list
        TString                 fOutputFileName;        // output file name
        TFile*                  fOutputFile;            // output file
        Int_t                   fCentralityBin;         // centrality bin
        Bool_t                  fdPhidPt;               // unfold as function of delta phi
        TH2D*                   fDetectorResponse;      // detector response
        Double_t                fBeta;                  // regularization strength
        Double_t                fBetaPerDOF;            // beta per DOF
        Bool_t                  fAvoidRoundingError;    // set dpt to zero for small values far from the diagonal
        unfoldingAlgorithm      fUnfoldingAlgorithm;    // algorithm used for unfolding
        prior                   fPrior;                 // prior for unfolding
        TArrayD*                fBinsTrue;              // pt true bins
        TArrayD*                fBinsRec;               // pt rec bins
        Int_t                   fSVDDraw;               // svd iteration 
        Int_t                   fIterMin;               // min number of iterations (SVD)
        Int_t                   fIterMax;               // max number of iterations (SVD)
        Bool_t                  fSVDToy;                // use toy to estimate coveriance matrix for SVD method
        Float_t                 fJetRadius;             // jet radius (for SVD toy)
        // members, set internally
        TH2D*   fDeltaPtDeltaPhi;               // delta pt delta phi distribution
        TH2D*   fJetPtDeltaPhi;                 // jet pt delta phi distribution
        TH1D*   fSpectrumIn;                    // in plane jet pt spectrum
        TH1D*   fSpectrumOut;                   // out of plane jet pt spectrum
        TH1D*   fDptInDist;                     // in plane dpt distribution
        TH1D*   fDptOutDist;                    // out of plane dpt distribution
        TH2D*   fDptIn;                         // in plane dpt matrix
        TH2D*   fDptOut;                        // out plane dpt matrix
        TH2D*   fFullResponseIn;                // full response matrix, in plane
        TH2D*   fFullResponseOut;               // full response matrix, out of plane
        TH1D*   fUnfoldedIn;                    // unfolded in plane spectrum
        TH1D*   fUnfoldedOut;                   // unfolded out ofplane spectrum
        // copy and assignment 
        AliJetFlowTools(const AliJetFlowTools&);             // not implemented
        AliJetFlowTools& operator=(const AliJetFlowTools&);  // not implemented
};
#endif
//_____________________________________________________________________________
