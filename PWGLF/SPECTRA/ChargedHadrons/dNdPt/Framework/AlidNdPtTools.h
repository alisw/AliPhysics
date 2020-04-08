PWGLF/SPECTRA/ChargedHadrons/dNdPt/Framework/AliAnalysisTaskBaseWeights.h/// \class AlidNdPtTools
/// \brief Collection of functionality used in dNdPt anlysis
///
/// all methods are static
/// AddAxis,CreateHist,CreateLogHist for easy creation of histograms
///
/// \author Michael Linus Knichel <michael.linus.knichel@cern.ch>, CERN
/// \date Mar 25, 2019

#ifndef AlidNdPtTools_H
#define AlidNdPtTools_H

#include "THnSparse.h"
#include <vector>

class TGraph;
class AliESDtrackCuts;
class AliMCParticle;
class AliMCEvent;

class AlidNdPtTools : public TObject
{
    public:
        virtual                 ~AlidNdPtTools() = 0;

    public:
        enum ParticleType      { kUndefined = -1, kEl = 0, kMu = 1, kPi = 2, kKa = 3, kPr = 4, kOther = 5, kSigmaP = 6, kSigmaM = 7, kXi = 8, kOmega = 9};
        enum ProductionType    { kUnknown = -1, kPrim = 0, kSecDecay = 1, kSecMaterial = 2};
        enum EventType         { kInvalidProcess = -1, kElastic = 0, kND = 1, kDD = 2, kCD = 3, kSD = 4};

        static Long64_t        FillHist(THnBase* s, Double_t x0, Double_t x1=0, Double_t x2=0, Double_t x3=0, Double_t x4=0, Double_t x5=0, Double_t x6=0, Double_t x7 =0, Double_t x8 =0, Double_t x9 =0, Double_t x10 =0, Double_t x11 =0);
        static Long64_t        FillHist(Double_t w, THnBase* s, Double_t x0, Double_t x1=0, Double_t x2=0, Double_t x3=0, Double_t x4=0, Double_t x5=0, Double_t x6=0, Double_t x7 =0, Double_t x8 =0, Double_t x9 =0, Double_t x10 =0, Double_t x11 =0);
        static Long64_t        FillHistWeighted(THnBase* s, std::vector<double> const& val, double weight);
        static Int_t           AddAxis(const char* label, Int_t nbins, Double_t xmin, Double_t xmax, const char* option = 0);                    // options: <none>
        static Int_t           AddAxis(const char* label, const char* title, Int_t nbins, Double_t xmin, Double_t xmax, const char* option = 0); // options: <none>
        static Int_t           AddAxis(const char* label, Int_t nbins, Double_t* xbins, const char* option = 0);        // options: <none>
        static Int_t           AddAxis(const char* label, const char* title, Int_t nbins, Double_t* xbins, const char* option = 0);        // options: <none>
        static Int_t           AddAxis(const char* label, const char* title, const char* option);                                          // options: pt
        static Int_t           AddAxis(const char* label, const char* option);                                          // options: pt
        static Int_t           AddAxis(const char* option);                                                             // options: pt
        static THnSparseD*     CreateHist(const char* name);
        static void            ResetHist() { if (fSparseTmp) { delete fSparseTmp; fSparseTmp=0; } }
        static TH1D*           CreateLogHist(const char* name, const char* title);
        static TH1D*           CreateLogHist(const char* name);
        static void            Log(TH1D* h, const char* name) { if (h) h->Fill(name,1); }
        static Double_t        MCScalingFactor(ProductionType prod, ParticleType part, Double_t pt);          // this is a temp solution, to be replace by MCSpectraWeights
                               //quick and dirty hard coded solution to be replaced
                               //periodindex 0: pp 13 TeV, LHC18b  -- this is the (new) default
                               //periodindex 1: pp 5TeV, LHC17pq   -- previoulsy that was the default! be aware!
        static Double_t        MCScalingFactor(AliMCParticle* part, AliMCEvent* event, Int_t systflag = 0, Int_t periodindex=0);  // access to the scaling

        static ParticleType    ParticleTypeFromPDG(Int_t pdgCode);
        static void            Range2Pi(Double_t &val) {while (val>=2*TMath::Pi()) val-=2*TMath::Pi(); while (val<0) val+=2*TMath::Pi(); } // change range: 0 <= val < 2Pi
        static void            Range1Pi(Double_t &val) {while (val>=TMath::Pi()) val-=2*TMath::Pi(); while (val<-TMath::Pi()) val+=2*TMath::Pi(); }  // change range: -Pi <= val < Pi

        static AliESDtrackCuts* CreateESDtrackCuts(const char* option, int _cutMode=100, bool _SaveHistos=false); // options

    private:
        static THnSparseD*      fSparseTmp;    //! temporary histogram for internal use only
        static TGraph*          fGsscale;       // graph with scaling factors (nominal)
        static TGraph*          fGsscale1;      // graph with scaling factors (syst up)
        static TGraph*          fGsscale2;      // graph with scaling factors (syst down)
        static TGraph*          fGsscaleB;       // graph with scaling factors (nominal)
        static TGraph*          fGsscaleB1;      // graph with scaling factors (syst up)
        static TGraph*          fGsscaleB2;      // graph with scaling factors (syst down)
    /// \cond CLASSIMP
    ClassDef(AlidNdPtTools, 4);
    /// \endcond
};

#endif
