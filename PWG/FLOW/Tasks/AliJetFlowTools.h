/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. * 
 * See cxx source for full Copyright notice */ 
 /* $Id$ */

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
#include "TFile.h"

//_____________________________________________________________________________
class AliJetFlowTools {
    public: 
        AliJetFlowTools();
    protected:
        ~AliJetFlowTools();
    public:
        // setters, setup the procedure
        void            SetInputList(TList* list)       {fInputList             = list;}
        void            SetOutputFileName(TString name) {fOutputFileName        = name;}
        void            CreateOutputList(TString name) {
            // create a new output list and add it to the full output
            fActiveString = name;
            fOutputList = new TList();
            fOutputList->SetOwner(kTRUE);
            fOutputArray->Add(fOutputList);
        }
        void            SetCentralityBin(Int_t bin)     {fCentralityBin         = bin;}
        void            SetdPhidPpt(Bool_t inout)       {fdPhidPt               = inout;}
        void            SetDetectorResponse(TH2D* dr)   {fDetectorResponse      = dr;}
        void            SetBinsTrue(TArrayD* bins)      {fBinsTrue              = bins;}
        void            SetBinsRec(TArrayD* bins)       {fBinsRec               = bins;}
        void            SetBeta(Double_t beta)          {fBeta                  = beta;}
        void            SetBetaDOF(Double_t betaDOF)    {fBetaPerDOF            = betaDOF;}
        void            Make();
        void            Finish() {
            // write the output lists to a file
            TFile* f = new TFile(fOutputFileName.Data(), "RECREATE");
            for(Int_t i(0); i < fOutputArray->GetSize(); i++) {
                if((TList*)fOutputArray->At(i)) {
                    f->mkdir(Form("slot_%i", i));
                    f->cd(Form("slot_%i", i));
                    ((TList*)fOutputArray->At(i))->Write();
                }
            }
            f->Close();
        }
        // static const helper functions, mainly histogram manipulation
        static TH1D*    ResizeXaxisTH1D(TH1D* histo, Int_t low, Int_t up, TString suffix = "");
        static TH2D*    ResizeYaxisTH2D(TH2D* histo, TArrayD* x, TArrayD* y, TString suffix = "");
        static TH2D*    NormalizeTH2D(TH2D* histo);
        static TH1D*    GetUnfoldingTemplate(TH1D* histo, TArrayD* bins, TString suffix = "");
        static TH2D*    RebinTH2DX(TH2D* histo, TArrayD* bins, TString suffix = "");
        static TH2D*    RebinTH2DY(TH2D* histo, TArrayD* bins);
        static TH2D*    MatrixMultiplicationTH2D(TH2D* A, TH2D* B, TString name = "CombinedResponse");
        static TH1D*    NormalizeTH1D(TH1D* histo, Double_t scale = 1.);
        static TGraphErrors*    divide_histos(  TH1 *h1 = 0x0, 
                                                TH1* h2 = 0x0, 
                                                Double_t xmax=-1.); 
    private:
        Bool_t          PrepareForUnfolding(); 
        void            GetRatios();
        TH1D*           UnfoldSpectrum( TH1D* resizedJetPt, 
                                        TH2D* resizedResonse,
                                        TH1D* kinematicEfficiency,
                                        TH1D* unfoldingTemplate,
                                        TString suffix);
        TMatrixD*       CalculatePearsonCoefficients(TMatrixD* covmat);
        // members, accessible via setters
        TString         fActiveString;          // identifier of active output
        TList*          fInputList;             // input list
        TList*          fOutputList;            // currently active output list
        TObjArray*      fOutputArray;           // output list
        TString         fOutputFileName;        // output file name
        Int_t           fCentralityBin;         // centrality bin
        Bool_t          fdPhidPt;               // unfold as function of delta phi
        TH2D*           fDetectorResponse;      // detector response
        Double_t        fBeta;                  // regularization strength
        Double_t        fBetaPerDOF;            // beta per DOF
        TArrayD*        fBinsTrue;              // pt true bins
        TArrayD*        fBinsRec;               // pt rec bins
        // members, set internally
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
