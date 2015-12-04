#ifndef AliAnalysisTaskJetFlowMC_H
#define AliAnalysisTaskJetFlowMC_H

// root includes
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h" 
#include "TRandom.h"
// aliroot includes
#include "AliAnalysisTaskSE.h"
// forward declarations
class TList;
class TClonesArray;
class TArrayI;
class TVirtualMCDecayer;
class AliPicoTrack;

class AliAnalysisTaskJetFlowMC : public AliAnalysisTaskSE 
{
    public:
        // enumerators
        enum detectorType       {kVZEROA, kVZEROC, kVZEROComb, kFixedEP};  // detector that was used
        // constructors, destructor
        AliAnalysisTaskJetFlowMC();
        AliAnalysisTaskJetFlowMC(const char *name, Bool_t qa = kFALSE, Int_t seed = 0);
        virtual ~AliAnalysisTaskJetFlowMC();
        void    UserCreateOutputObjects();
        TH1F*   BookTH1F(const char* name, const char* x, Int_t bins, Double_t min, Double_t max, Int_t c = -1, Bool_t append = kTRUE);
        TH2F*   BookTH2F(const char* name, const char* x, const char* y, Int_t binsx, Double_t minx, Double_t maxx, Int_t binsy, Double_t miny, Double_t maxy, Int_t c = -1, Bool_t append = kTRUE);

        void    UserExec(Option_t *option);
        void    SetDebugMode(Bool_t d)                          { fDebug = d;}
        void    SetTracksInName(const char *name)               { fTracksInName     = name; }
        void    SetTracksOutName(const char *name)              { fTracksOutName    = name; }
        // additional setters
        void    SetCentralityClasses(TArrayI* c)                { fCentralityClasses = c; }
        void    SetReferenceDetector(detectorType type)         { fDetectorType = type; }
        void    SetDifferentialV2(TF1* v2, Int_t c = 0)         { fFuncDiffV2[c] = v2; }
        void    SetDifferentialV3(TF1* v3, Int_t c = 0)         { fFuncDiffV3[c] = v3; }
        void    SetDifferentialV2(TH1* v2, Int_t c = 0)         { fHistDiffV2[c] = v2; }
        void    SetDifferentialV3(TH1* v3, Int_t c = 0)         { fHistDiffV3[c] = v3; }
        void    SetIntegratedV2(TH1* v2)                        { fHistIntV2 = v2; }
        void    SetIntegratedV3(TH1* v3)                        { fHistIntV3 = v3; }
        void    SetTrackSpectrum(TF1* ts)                       { fTrackSpectrum = ts; }
        void    SetRandomizeEta(Bool_t b)                       { fRandomizeEta = b; }
        void    SetMultiplicity(Int_t m)                        { fMult = m; }
        void    SetReuseTracks(Bool_t r)                        { fReuseTracks = r; }
        void    SetSingleFragmentationJetSpectrum(TF1* js)      { fJetSpectrumSF = js; }
        void    SetNoOfSFJets(Int_t n)                          { fNoOfSFJets = n; }
        void    SetDecayer(TVirtualMCDecayer* d, Int_t c = 1)   { fDecayer = d; fDecayerIterations = c;}
        // additional methods
        void    V2AfterBurner(Double_t& phi, Double_t& eta, Double_t& pt) const;
        void    V3AfterBurner(Double_t& phi, Double_t& eta, Double_t& pt) const;
        void    InjectSingleFragmentationJetSpectrum(Int_t nacc);
        void    CalculateEventPlane();
        // inlined helper calculations
        Double_t        GetTrackPt()       const                { return fTrackSpectrum->GetRandom();}
        Double_t        GetTrackEta()       const               { return gRandom->Uniform(-.9, .9); }
        /* inline */    Double_t GetV2(Double_t pt) const { 
            return (fFuncDiffV2[fCenBin]) ? fFuncDiffV2[fCenBin]->Eval(pt) :
            fHistDiffV2[fCenBin]->GetBinContent(fHistDiffV2[fCenBin]->GetXaxis()->FindBin(pt));
        }
        /* inline */    Double_t GetV3(Double_t pt) const { 
            return (fFuncDiffV3[fCenBin]) ? fFuncDiffV3[fCenBin]->Eval(pt) : 
            fHistDiffV3[fCenBin]->GetBinContent(fHistDiffV3[fCenBin]->GetXaxis()->FindBin(pt));
        }
        /* inline */    void GetFlowFluctuation(Double_t& vn) const {
            vn += TMath::Sqrt(2*(vn*.25)*(vn*.25))*TMath::ErfInverse(2*(gRandom->Uniform(0, fFlowFluctuations))-1); 
        }
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
        /* inline */    void SampleVnFromTF1(Double_t &phi) const {
        phi = (fFuncVn) ? fFuncVn->GetRandom(0., TMath::TwoPi()) : 0; }
        /* inline */    void FillHistogramsOriginalData(Double_t& pt, Double_t& eta, Double_t& phi) const {
            fHistOriginalSpectrum[fCenBin]->Fill(pt);    fHistOriginalEtaPhi[fCenBin]->Fill(eta, phi);
            fHistOriginalDeltaPhi[fCenBin]->Fill(PhaseShift(phi-fPsi2, 2));
        }
        /* inline */    void FillHistogramsToyData(Double_t& pt, Double_t& eta, Double_t& phi, Double_t& vn) const {
            fHistToySpectrum[fCenBin]->Fill(pt);         fHistToyEtaPhi[fCenBin]->Fill(eta, phi);
            fHistToyDeltaPhi[fCenBin]->Fill(PhaseShift(phi-fPsi2, 2));  fHistToyVn[fCenBin]->Fill(pt, vn);
        }
        Int_t   InsertDecayDaughters(AliPicoTrack* mother, TClonesArray* daughters);
        Int_t   InsertDecayDaughters(Double_t pt, Double_t phi, Double_t eta, Double_t mass, Short_t charge, TClonesArray* daughters);
        void    Terminate(Option_t* option);
        void    PrintInfo() const;
    protected:
        Bool_t          fQA;                    //! save QA plots
        TString         fTracksOutName;         // name of output track array
        TString         fTracksInName;          // name of input track array
        TClonesArray   *fTracksIn;              //! track array in
        TClonesArray   *fTracksOut;             //! track array out
        Bool_t          fReuseTracks;           // use original event as template
        Int_t           fMult;                  // multiplicity of new event
        Int_t           fCenBin;                //! centrality bin
        TArrayI*        fCentralityClasses;     // centrality classes (max 10) 
        TF1*            fFuncVn;                //! vn function
        TList*          fOutputList;            // output list
        TF1*            fTrackSpectrum;         // track pt spectrum
        Bool_t          fRandomizeEta;          // randomize eta
        TF1*            fJetSpectrumSF;         // single fragmentation spectrum of jets
        Int_t           fNoOfSFJets;            // number of single fragmentation jets that will be added
        TF1*            fFuncDiffV2[10];        // differential v2 of charged tracks
        TF1*            fFuncDiffV3[10];        // differential v3 of charged tracks
        TH1*            fHistDiffV2[10];        // differential v2 of charged tracks
        TH1*            fHistDiffV3[10];        // differential v3 of charged tracks
        TH1*            fHistIntV2;             // integrated v2 of charged tracks
        TH1*            fHistIntV3;             // integrated v3 of charged tracks
        Float_t         fFlowFluctuations;      // introduce gaussian flow fluctuations of this magnitude
        Int_t           fMaxNumberOfIterations; // max number of iterations for afterburner
        Double_t        fPsi2;                  //! 2nd order event plane orientation
        Double_t        fPsi3;                  //! 3rd order event plane orientation
        Double_t        fPrecisionPhi;          // afterburner precision
        detectorType    fDetectorType;          // type of detector from which the EP is taken
        // output histograms
        TH1F*           fHistOriginalSpectrum[10];      //! original pt spectrum of accepted tracks
        TH2F*           fHistOriginalEtaPhi[10];        //! original eta phi spectrum of accepted tracks
        TH1F*           fHistToySpectrum[10];           //! spectrum of toy (generated) tracks
        TH2F*           fHistToyEtaPhi[10];             //! eta phi spectrum of toy (generated) tracks
        TH1F*           fHistOriginalDeltaPhi[10];      //! original delta phi spectrum
        TH1F*           fHistToyDeltaPhi[10];           //! delta phi spectrum of toy (generated) tracks
        TH2F*           fHistToyVn[10];                 //! generated differential vn values (should equal the differential spectrum)
        TH1F*           fHistSFJetSpectrum;             //! spectrum of generated sf jets
        TH2F*           fHistSFJetEtaPhi;               //! eta phi of generated sf jets
        // decayer
        TVirtualMCDecayer*      fDecayer;               // decayer, needs to be set in macro (avoid dependencies)
        Int_t                   fDecayerIterations;     // max no of possible decay vertices from one track
        TClonesArray*           fDecayerCache[25];      //! cached tparticle's, used as intermediate buffer
        TClonesArray*           fDecayerResults;        //! decayer results
    private:
        AliAnalysisTaskJetFlowMC(const AliAnalysisTaskJetFlowMC&);            // not implemented
        AliAnalysisTaskJetFlowMC &operator=(const AliAnalysisTaskJetFlowMC&); // not implemented

        ClassDef(AliAnalysisTaskJetFlowMC, 4); // Task to generate toy mc PicoTracks based on real events
};
#endif
