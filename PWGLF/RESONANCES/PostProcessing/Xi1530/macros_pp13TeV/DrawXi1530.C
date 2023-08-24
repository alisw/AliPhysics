#if !defined(__CINT__) || defined(__CLING__)
#include "RooGlobalFunc.h"
#endif
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TH1D.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPRegexp.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVectorD.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include "TMath.h"
#include "Math/MinimizerOptions.h"

#include <TH1D.h>
#include "RooChi2Var.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooMath.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "TAxis.h"
#include "TCanvas.h"

//--------Xi(1530)0 post processing macro ------------//
// 2019
// Name:  Bong-Hwi Lim
// Email: bong-hwi.lim@cern.ch
//
// Based on BSHelper class(Beomsoo Kim)
// Advised by Beomkyu Kim(kimb@cern.ch)

//---------------------------------------
//  typedef, const, global
//---------------------------------------

typedef std::vector<int> Int1D;
typedef std::vector<Int1D> Int2D;
typedef std::vector<Double_t> Double1D;
typedef std::vector<Double1D> Double2D;

const int min_int = std::numeric_limits<int>::min() + 1;

//---------------------------------------
//  Utils
//---------------------------------------

void ErrorExit(TString error, int errnum = 1) {
    cout << "ERROR: " << error << endl;
    gSystem->Exit(errnum);
}
void __DEBUG_BS(int i, const char* fmt, ...) {
    char temp[1024];
    va_list marker;
    va_start(marker, fmt);
    vsprintf(temp, fmt, marker);
    va_end(marker);

    std::cout << Form("DEBUG %d : %s", i, temp) << std::endl;
}

//====================================
//   Manage Root
//====================================

TFile* LoadRoot(TString filename) {
    auto f = gROOT->GetFile(filename);
    if (!f)
        f = TFile::Open(filename);
    if (!f)
        ErrorExit("No File " + filename);
    return f;
}

//---------------------------------------
//  RANGER
//---------------------------------------
class BSRanger {
   public:
    typedef int value_type;
    struct iterator {
        iterator(size_t counter) : counter(counter) {}
        iterator operator++() { return iterator(++counter); }
        bool operator!=(iterator o) { return counter != o.counter; }
        value_type operator*() { return value_type{counter}; }

       private:
        value_type counter;
    };
    BSRanger(int begin, int end) { SetRange(begin, end); }
    void SetRange(int begin, int end) {
        if (begin > end)
            ErrorExit("being>end");
        fBegin = begin;
        fEnd = end;
    }
    iterator begin() { return iterator(fBegin); }
    iterator end() { return iterator(fEnd + 1); }

   private:
    int fBegin;
    int fEnd;
};

template <class T>
BSRanger range(T& t, int begin = 0, int end = -1) {
    return BSRanger(begin, end < 0 ? t.size() - 1 : end);
}
BSRanger range(int end) {
    return BSRanger(0, end - 1);
}
BSRanger range(int begin, int end) {
    return BSRanger(begin, end);
}

BSRanger bin_range(int begin, int end) {
    if (begin < 1)
        ErrorExit("begin for bin_range must larger than 0 ");
    return range(begin, end);
}
BSRanger bin_range(int end) {
    return bin_range(1, end);
}
template <class T>
BSRanger bin_range(const std::vector<T>& t, int begin = 1, int end = min_int) {
    return bin_range(begin,
                     (end < 0 ? std::min(int(t.size() - 1), abs(end)) : end));
}
//---------------------------------------

//====================================
// HistManager
//====================================
class BSTHnSparseHelper {
   public:
    BSTHnSparseHelper();
    BSTHnSparseHelper(THnSparse* h) : fH(h) { LoadBinFromHist(); }
    BSTHnSparseHelper(TObject* o);
    THnSparse* operator->() { return fH; }
    THnSparse* Data() { return fH; }
    THnSparse* Val() { return fH; }
    THnSparse* Hist() { return fH; }

    //================
    //  AXIS
    //================
    int GetAxisID(TString name);
    TAxis* GetAxis(int i);
    TAxis* GetAxis(TString name) { return GetAxis(GetAxisID(name)); }
    int GetNdim() { return fH->GetNdimensions(); };
    Int_t GetNbins(int i) { return GetAxis(i)->GetNbins(); }
    Int_t GetNbins(TString n) { return GetNbins(GetAxisID(n)); }
    void PrintAxis(Option_t* opt = "");
    TH1D* GetTH1(TString name, Int_t xDim, Int1D bin, Option_t* opt = "");
    TH1D* GetTH1(TString name, Int_t xDim, Int2D bin, Option_t* opt = "");
    TH1D* GetTH1(TString name, Int1D bin, Option_t* opt = "") {
        return GetTH1(name, GetNdim() - 1, bin, opt);
    }

    BSRanger BinRange(int i) { return bin_range(GetBin(i)); }
    BSRanger BinRange(TString n) { return bin_range(GetBin(n)); }
    void ClearBin() { LoadBinFromHist(); }
    void ClearBin(int i) { LoadBinFromHist(i); }
    const Double1D& GetBin(int i) { return fCBins[i]; }
    const Double1D& GetBin(TString n) { return GetBin(GetAxisID(n)); }
    void SetBin(Double2D bins);
    void SetBin(int i, Double1D bins);
    void SetBin(TString name, Double1D bins) { SetBin(GetAxisID(name), bins); }
    static BSTHnSparseHelper Load(TString name, TObject* clist = nullptr);
    static BSTHnSparseHelper Load(TString name,
                                  TString fname,
                                  TString lname = "");

   private:
    void LoadBinFromHist();
    void LoadBinFromHist(int iaxis);
    THnSparse* fH;
    Double2D fCBins;
    Double2D fCBinsB;
};

//__________________________________________________________
BSTHnSparseHelper::BSTHnSparseHelper() {
}

//__________________________________________________________
BSTHnSparseHelper::BSTHnSparseHelper(TObject* o) {
    if (!o)
        ErrorExit("No Object");
    fH = dynamic_cast<THnSparse*>(o);
    if (!fH)
        ErrorExit(
            Form("%s is not THnSparse but %s", o->GetName(), o->ClassName()));
    LoadBinFromHist();
}

//__________________________________________________________
TAxis* BSTHnSparseHelper::GetAxis(int i) {
    if (i < 0 || i >= GetNdim())
        ErrorExit(Form("Wrong Axis Index %d of %s", GetNdim(), fH->GetName()));
    return fH->GetAxis(i);
}
//__________________________________________________________
int BSTHnSparseHelper::GetAxisID(TString name) {
    for (int i = 0; i < fH->GetNdimensions(); i++)
        if (name == fH->GetAxis(i)->GetName())
            return i;
    fH->Print();
    PrintAxis();
    ErrorExit("No Axis " + name);
    return -1;
}

//__________________________________________________________
void BSTHnSparseHelper::PrintAxis(Option_t* opt) {
    TString opts = opt;
    for (int i = 0; i < fH->GetNdimensions(); i++) {
        TAxis* ax = fH->GetAxis(i);
        cout << Form("%2d: %10s", i, fH->GetAxis(i)->GetName());
        if (opts == "all") {
            cout << Form(" (nbin/min/max) = %4d %6.2f %8.2f ::", ax->GetNbins(),
                         ax->GetXmin(), ax->GetXmax());
            for (int j = 1; j <= ax->GetNbins() && j < 10; j++) {
                const char* label = ax->GetBinLabel(j);
                if (!TString(label).IsNull())
                    cout << Form(" %s", label);
                else
                    cout << Form(" %.2f", ax->GetBinLowEdge(j));
            }
        }
        cout << endl;
    }
}

//__________________________________________________________
TH1D* BSTHnSparseHelper::GetTH1(TString name,
                                Int_t xDim,
                                vector<int> bin,
                                Option_t* opt) {
    Int2D newbin;
    for (auto b : bin) {
        newbin.push_back({b, b});
    }
    return GetTH1(name, xDim, newbin, opt);
}
//__________________________________________________________
TH1D* BSTHnSparseHelper::GetTH1(TString name,
                                Int_t xDim,
                                Int2D bin,
                                Option_t* opt) {
    //== Histogram Name and Title
    if (name.EndsWith("-"))
        name += Form("%sP%02d", fH->GetName(), xDim);
    TString title = fH->GetTitle();
    for (UInt_t i = 0; i < bin.size(); i++) {
        Int1D b = bin[i];
        Int1D nb = {0, 0};
        int maxnbin = fCBinsB[i].size() - 2;
        auto ax = GetAxis(i);
        if (b[0] < -1 || b[1] < -1 || b[0] > maxnbin || b[1] > maxnbin ||
            b[0] > b[1])
            ErrorExit(Form("Wrong bin : %i : %s : %d %d", i, ax->GetName(),
                           b[0], b[1]));
        if (int(i) == xDim || b.size() == 0 || (b[0] <= 0 && b[1] <= 0)) {
        } else {
            if (b[0] > 0)
                nb[0] = fCBinsB[i][b[0]];
            if (b[1] > 0)
                nb[1] = fCBinsB[i][b[1] + 1] - 1;
            // if( b0 > b1 ) b1 = b1;
            name += Form("%s%03d%03d", GetAxis(i)->GetName(), b[0], b[1]);
            TString label[2];
            label[0] = ax->GetBinLabel(nb[0]);
            label[1] = ax->GetBinLabel(nb[1]);
            if (label[0] == "")
                label[0] = Form("%.6f", ax->GetBinLowEdge(nb[0]));
            else
                label[0] = Form("%d", nb[0]);  // TODO Temp
            if (label[1] == "")
                label[1] = Form("%.6f", ax->GetBinUpEdge(nb[1]));
            else
                label[1] = Form("%d", nb[1]);  // TODO Temp
            TPMERegexp("\\.?0+$").Substitute(label[0], "");
            TPMERegexp("\\.?0+$").Substitute(label[1], "");
            title += Form(" %s:%s", ax->GetName(), label[0].Data());
            if (label[0] != label[1])
                title += Form("-%s", label[1].Data());
        }
        fH->GetAxis(i)->SetRange(nb[0], nb[1]);
    }
    auto h = fH->Projection(xDim, opt);
    h->SetNameTitle(name, title);
    return h;
}

//__________________________________________________________
void BSTHnSparseHelper::SetBin(int iaxis, Double1D bins) {
    auto ax = GetAxis(iaxis);
    if (bins.size() == 0)
        LoadBinFromHist(iaxis);
    else {
        fCBins[iaxis] = bins;
        fCBinsB[iaxis] = {0};
        for (UInt_t ib = 0; ib < bins.size(); ib++) {
            int bin = ax->FindBin(bins[ib]);
            fCBinsB[iaxis].push_back(ax->FindBin(bins[ib]));
        }
    }
}
//__________________________________________________________
void BSTHnSparseHelper::SetBin(Double2D bins) {
    for (UInt_t i = 0; i < bins.size(); i++)
        SetBin(i, bins[i]);
}

//__________________________________________________________
BSTHnSparseHelper BSTHnSparseHelper::Load(TString name, TObject* clist) {
    if (!clist)
        clist = gDirectory;
    auto o = clist->FindObject(name);
    if (!o)
        ErrorExit("No THnSparse " + name);
    return BSTHnSparseHelper(o);
}
//__________________________________________________________
BSTHnSparseHelper BSTHnSparseHelper::Load(TString name,
                                          TString fname,
                                          TString lname) {
    TFile* f = gROOT->GetFile(fname);
    if (!f)
        f = TFile::Open(fname);
    if (!f)
        ErrorExit("No File : " + fname);
    TObject* cl = f;
    if (!lname.IsNull())
        cl = f->Get(lname);
    return Load(name, cl);
}

//__________________________________________________________
void BSTHnSparseHelper::LoadBinFromHist(int iaxis) {
    auto ax = GetAxis(iaxis);
    fCBinsB[iaxis] = {0};
    for (int ib = 1; ib <= ax->GetNbins() + 1; ib++) {
        fCBins[iaxis].push_back(ax->GetBinLowEdge(ib));
        fCBinsB[iaxis].push_back(ib);
    }
}
//__________________________________________________________
void BSTHnSparseHelper::LoadBinFromHist() {
    fCBins.resize(GetNdim());
    fCBinsB.resize(GetNdim());
    for (int i = 0; i < GetNdim(); i++)
        LoadBinFromHist(i);
}
//__________________________________________________________
// Initial values

// constants
const int POLdegree = 1;
const float inf = 1e20;
const float zero = 1e-20;

// Input files
// TString datafile = "ESD1517_1519";  // ESD
// TString datafile = "AOD439_441";  // AOD
//TString datafile = "XeXetest";  // AOD
TString datafile = "PbPb_LHC18qr";  // AOD


// char const* rsnmcfile = "LHC17f2b_cent"; // pPb

TString rsnmcfile = "LHC19h6";        // LHC16 pass2
TString rsnmcfile2 = "LHC18c6b_part1_298";  // LHC16
TString rsnmcfile3 = "LHC18c6b_part2_298";  // LHC16
TString rsnmcfile4 = "LHC18c6a_part1_298";  // LHC17
TString rsnmcfile5 = "LHC18c6a_part2_298";  // LHC17
TString rsnmcfile6 = "LHC18c6a_part3_298";  // LHC17
TString rsnmcfile7 = "LHC18c6c_298";        // LHC15

TString genmcfile = "LHC18l8b";

TString inputDirectory = "Xi1530MB";

// FIT SETUP -----------------
double peakRange[2] = {1.512, 1.552};

char name[500] = "Voigt fit";
TString formula = "pol2(0)";
TString formula_MC = "pol0(0)";
double par[3] = {1.532, 0.003, 0.0091};
double par_2nd[3] = {1.532, 0.0025, 0.0091};
double fix[3] = {0, 0, 1};
double fix_2nd[3] = {1, 0, 1};
double fix_3rd[3] = {1, 1, 1};
double fix_MC[3] = {0, 0, 1};
int vp = 0;
TF1* fb_after;
Double_t bkgintegral = 0;
bool bincount = false;
bool fHM = false;
bool fbkgfit = false;
// ---------------------------

// Mix? LS? Background -------
int bkgtype = 3;  // 2 for LS, 3 for Mix
// ---------------------------

// Analysis options ----------
// pT bins
Double1D ptbin = {0.0, 0.8, 1.2, 1.6, 2.0, 2.4, 3.2, 4.0, 4.8, 5.6, 8.8, 15};
Double1D ptzero = {zero, zero, zero, zero, zero, zero,
                   zero, zero, zero, zero, zero, zero};
Double1D pt_points_e = {};

// Norm region
Int_t normleft = 0;
Double1D NormalizeRange_L = {1.48, 1.50}; // Norminal
Int_t normright = 1;
Double1D NormalizeRange_R = {1.6, 1.8}; // Norminal

// Axis Range
Double1D DrawRange = {1.448, 1.75};
Double1D FitRange = {1.484, 1.675}; // +-rstep
double lstep = 0.001;
double rstep = 0.012;
Double1D IntegralRange = {1.516, 1.548}; // +-
Double1D IntegralRangeMC = {1.4, 1.7};

// Rebin
Int_t rebin = 4;

// Default Canvas
Double_t w = 960;
Double_t h = 720;

TH1D* GetNorBkg(TH1D* hSig, TH1D* hBkg, int isnormleft, int isnormright);
TF1* VoigtFit(TH1D* hSig,
              TString bkgformula,
              double* peakRange,
              double* par,
              double* fix);
TH1D* MakeHistfromArray(char const* name,
                        Double1D dArray,
                        Double1D eArray,
                        Double1D ptbin,
                        const char* foption = "");
double myLevyPt(Double_t* x, Double_t* par);
vector<double> GetNormalisationFactor(double multi_start, double multi_end);
// File handling
TFile* LoadXi1530Results(TString name, TString runnum);
TObject* LoadXi1530ResultList(TFile* fh, TString clistname);
TObject* LoadXi1530ResultList(TString fname,
                              TString clistname,
                              TString runnum = "");
//__________________________________________________________
//__________________________________________________________
//_______________________DrawXi1530_________________________
//__________________________________________________________
//__________________________________________________________
void DrawXi1530(const int sys = 1,
                double multi_start = 0,
                double multi_end = 10,
                char const* inputOptions = "PbPb",
                const int OptionNumber = 1) {
    // Start!
    cout << "== Xi(1530)0 Postprocessing Macro ==" << endl;
    cout << "Cut Systematic option(Default: 1): " << sys << endl;
    cout << "Multiplicity Range: " << multi_start << "-" << multi_end << endl;
    cout << "Input options: " << inputOptions << endl;
    cout << "Option number: " << OptionNumber << endl;

    Double1D centbin = {multi_start, multi_end};
    Double1D fullcentbin = {0, 100};
    Double1D fullcentbin_Effi = {0, 100};
    bool isVertexCutEnd = false;
    bool fpPb = false;
    bool fAA = false;
    bool is5TeV = false;
    bool useFurtherMC = false;
    // Reaction to the option
    TString Options = inputOptions;

    // For fit systematics
    // NOT USING NOW -> Replaced to FitVar, BinCount, NormVar
    if (Options.Contains("5TeV")){
        is5TeV = true;
        useFurtherMC = false;
        datafile = "17qp_wSDD";
        rsnmcfile = "LHC17q_MC";
        genmcfile = "LHC17q_GenMC";
        inputDirectory = "Xi1530MB";
    }
    if (Options.Contains("FitRange"))
        FitRange = {1.49, 1.57};
    if (Options.Contains("NormRange"))
        NormalizeRange_L = {1.48, 1.49};
    if (Options.Contains("BinCount"))
        bincount = true;
    if (Options.Contains("LikeSignBkg")){
        fix_2nd[1] = 1;
        bkgtype = 2;
    }
    if (Options.Contains("NormRight")) {
        normleft = 0;
        normright = 1;
    }
    if (Options.Contains("NormBoth")) {
        normleft = 1;
        normright = 1;
    }
    //---------------------------------------
    if (multi_end < 1) {  // Automatic HM mode
        datafile = "AOD439_441_HM1718";
        inputDirectory = "Xi1530HM";
        fHM = true;
        fAA = true;
        fullcentbin = {0, 0.1};
    }
    // Fit variation check
    if (Options.Contains("FitVarLm")){
        FitRange[0] = FitRange[0]-lstep*rebin*OptionNumber;
    }
    if (Options.Contains("FitVarLp")){
        FitRange[0] = FitRange[0]+rstep*OptionNumber;
    }
    if (Options.Contains("FitVarRm")){
        FitRange[1] = FitRange[1]-lstep*rebin*OptionNumber;
    }
    if (Options.Contains("FitVarRp")){
        FitRange[1] = FitRange[1]+rstep*OptionNumber;
    }
    if (Options.Contains("FitVarBothm")){
        FitRange[0] = FitRange[0]-lstep*rebin*OptionNumber;
        FitRange[1] = FitRange[1]+rstep*OptionNumber;
    }
    if (Options.Contains("FitVarBothp")){
        FitRange[0] = FitRange[0]+lstep*rebin*OptionNumber;
        FitRange[1] = FitRange[1]-rstep*OptionNumber;
    }
    // Norm change
    if (Options.Contains("NormVarm")){
        NormalizeRange_L[0] = NormalizeRange_L[0]-lstep*rebin*OptionNumber;
        NormalizeRange_L[1] = NormalizeRange_L[1]-lstep*rebin*OptionNumber;
        NormalizeRange_R[0] = NormalizeRange_R[0]-rstep*OptionNumber;
        NormalizeRange_R[1] = NormalizeRange_R[1]-rstep*OptionNumber;
    }
    if (Options.Contains("NormVarp")){
        NormalizeRange_L[0] = NormalizeRange_L[0]+lstep*rebin*OptionNumber;
        NormalizeRange_L[1] = NormalizeRange_L[1]+lstep*rebin*OptionNumber;
        NormalizeRange_R[0] = NormalizeRange_R[0]+rstep*OptionNumber;
        NormalizeRange_R[1] = NormalizeRange_R[1]+rstep*OptionNumber;
    }
    if (Options.Contains("NormVarLp")){
        NormalizeRange_L[1] = NormalizeRange_L[1]+lstep*rebin*OptionNumber;
    }
    if (Options.Contains("NormVarLm")){
        NormalizeRange_L[0] = NormalizeRange_L[0]-lstep*rebin*OptionNumber;
    }
    if (Options.Contains("NormVarRp")){
        NormalizeRange_R[1] = NormalizeRange_R[1]+rstep*OptionNumber;
    }
    if (Options.Contains("NormVarRm")){
        NormalizeRange_R[0] = NormalizeRange_R[0]-rstep*OptionNumber;
    }
    //Bin count varation
    if (Options.Contains("BinCountLm")){
        IntegralRange[0] = IntegralRange[0]-0.001*rebin*OptionNumber;
    }
    if (Options.Contains("BinCountLp")){
        IntegralRange[0] = IntegralRange[0]+0.001*rebin*OptionNumber;
    }
    if (Options.Contains("BinCountRm")){
        IntegralRange[1] = IntegralRange[1]-0.001*rebin*OptionNumber;
    }
    if (Options.Contains("BinCountRp")){
        IntegralRange[1] = IntegralRange[1]+0.001*rebin*OptionNumber;
    }
    if (Options.Contains("BinCountBothm")){
        IntegralRange[0] = IntegralRange[0]+0.001*rebin*OptionNumber;
        IntegralRange[1] = IntegralRange[1]-0.001*rebin*OptionNumber;
    }
    if (Options.Contains("BinCountBothp")){
        IntegralRange[0] = IntegralRange[0]-0.001*rebin*OptionNumber;
        IntegralRange[1] = IntegralRange[1]+0.001*rebin*OptionNumber;
    }
    if (Options.Contains("BkgFit")){
        fbkgfit = true;
        //fix_2nd[1] = 1;
        formula = "[0] + [1]*x + [2]*TMath::Sqrt( abs(x - 1.46) )"; // 1.321+0.139
        //formula = "pol2(0)";
        //FitRange={1.496,1.572};
        //FitRange[1] = 1.496;
        FitRange[1] = 1.6;
    }
    if (Options.Contains("BkgFitLm")){
        FitRange[0] = FitRange[0]-lstep*rebin*OptionNumber;
    }
    if (Options.Contains("BkgFitLp")){
        FitRange[0] = FitRange[0]+lstep*rebin*OptionNumber;
    }
    if (Options.Contains("BkgFitRm")){
        FitRange[1] = FitRange[1]-lstep*rebin*OptionNumber;
    }
    if (Options.Contains("BkgFitRp")){
        FitRange[1] = FitRange[1]+lstep*rebin*OptionNumber;
    }
    if (Options.Contains("BkgFitBothm")){
        FitRange[0] = FitRange[0]-lstep*rebin*OptionNumber;
        FitRange[1] = FitRange[1]+lstep*rebin*OptionNumber;
    }
    if (Options.Contains("BkgFitBothp")){
        FitRange[0] = FitRange[0]+lstep*rebin*OptionNumber;
        FitRange[1] = FitRange[1]-lstep*rebin*OptionNumber;
    }

    if (Options.Contains("MCcheck")){
        fullcentbin = centbin;
        fullcentbin_Effi = centbin;
    }
    if (Options.Contains("vertexcut")){
        isVertexCutEnd = true;
    }
    if (Options.Contains("pPb"))
        fpPb = true;
    if (Options.Contains("PbPb"))
        fAA = true;

    // Preparation
    for (int i = 0; i < (int)ptbin.size() - 1; i++)
        pt_points_e.push_back(ptbin.at(i + 1) - ptbin[i]);

    // Latex
    // for memo, small
    TLatex* t = new TLatex();
    t->SetNDC();
    t->SetTextSize(0.04);
    // for warning, big
    TLatex* t2 = new TLatex();
    t2->SetNDC();
    t2->SetTextSize(0.05);
    // for small pad, huge
    TLatex* t3 = new TLatex();
    t3->SetNDC();
    t3->SetTextSize(0.1);

    // Arrays
    Double1D fitmean = {};
    Double1D fitmean_err = {};
    Double1D fitsigma = {};
    Double1D fitsigma_err = {};
    Double1D fitgamma = {};
    Double1D fitgamma_err = {};

    Double1D RawYield = {};
    Double1D RawYield_err = {};

    Double1D Chi2 = {};
    Double1D NDF = {};
    Double1D EventCutRatio = {};
    Double1D EventCutRatio_e = {};
    Double1D Efficiency = {};
    Double1D Efficiency_e = {};

    Double1D fitmean_MC = {};
    Double1D fitmean_MC_err = {};
    Double1D fitsigma_MC = {};
    Double1D fitsigma_MC_err = {};
    Double1D fitgamma_MC = {};
    Double1D fitgamma_MC_err = {};

    //

    double triggereffi = inf;
    double triggereffi_e = zero;

    double vertexeffi = inf;
    double vertexeffi_e = zero;

    // Load DATA
    auto clist = LoadXi1530ResultList(datafile.Data(), inputDirectory);
    auto hInvMass = BSTHnSparseHelper::Load("hInvMass", clist);
    auto clist_MC = LoadXi1530ResultList(
        rsnmcfile.Data(), "Xi1530MB");  // From Resonance Injected MC
    auto hInvMass_MC = BSTHnSparseHelper::Load("hInvMass", clist_MC);
    auto hInvMass_MC_MB = BSTHnSparseHelper::Load("hInvMass", clist_MC);
    auto clist_MC_General = LoadXi1530ResultList(
         genmcfile.Data(), "Xi1530MB");  // From General Purpose MC
         //genmcfile, "Xi1530test");  // From General Purpose MC
    auto hInvMass_MC_General =
        BSTHnSparseHelper::Load("hInvMass", clist_MC_General);
    auto hInvMass_MC_General_MB =
        BSTHnSparseHelper::Load("hInvMass", clist_MC_General);
    
    // pT binning
    hInvMass.SetBin("Pt", ptbin);
    hInvMass_MC.SetBin("Pt", ptbin);
    hInvMass_MC_General.SetBin("Pt", ptbin);
    hInvMass_MC_MB.SetBin("Pt", ptbin);
    hInvMass_MC_General_MB.SetBin("Pt", ptbin);
 

    // Multiplicity percentile binning
    hInvMass.SetBin("Cent", centbin);
    hInvMass_MC.SetBin("Cent", centbin);          
    hInvMass_MC_MB.SetBin("Cent", fullcentbin_Effi);   // for MC reconstruction Efficiency, we use full bin
    hInvMass_MC_General.SetBin("Cent", centbin);  // for Trigger Efficiency
    hInvMass_MC_General_MB.SetBin("Cent", fullcentbin_Effi);  // for MC reconstruction Efficiency, we use full bin
    
    //further RsnMC
    BSTHnSparseHelper hInvMass_MC2;
    BSTHnSparseHelper hInvMass_MC3;
    BSTHnSparseHelper hInvMass_MC4;
    BSTHnSparseHelper hInvMass_MC5;
    BSTHnSparseHelper hInvMass_MC6;
    BSTHnSparseHelper hInvMass_MC7;

    BSTHnSparseHelper hInvMass_MC_MB2;
    BSTHnSparseHelper hInvMass_MC_MB3;
    BSTHnSparseHelper hInvMass_MC_MB4;
    BSTHnSparseHelper hInvMass_MC_MB5;
    BSTHnSparseHelper hInvMass_MC_MB6;
    BSTHnSparseHelper hInvMass_MC_MB7;

    if(useFurtherMC){
        auto clist_MC2 = LoadXi1530ResultList(
            rsnmcfile2.Data(), "Xi1530MB");  // From Resonance Injected MC
        auto clist_MC3 = LoadXi1530ResultList(
            rsnmcfile3.Data(), "Xi1530MB");  // From Resonance Injected MC
        auto clist_MC4 = LoadXi1530ResultList(
            rsnmcfile4.Data(), "Xi1530MB");  // From Resonance Injected MC
        auto clist_MC5 = LoadXi1530ResultList(
            rsnmcfile5.Data(), "Xi1530MB");  // From Resonance Injected MC
        auto clist_MC6 = LoadXi1530ResultList(
            rsnmcfile6.Data(), "Xi1530MB");  // From Resonance Injected MC
        auto clist_MC7 = LoadXi1530ResultList(
            rsnmcfile7.Data(), "Xi1530MB");  // From Resonance Injected MC
        hInvMass_MC2 = BSTHnSparseHelper::Load("hInvMass", clist_MC2);
        hInvMass_MC3 = BSTHnSparseHelper::Load("hInvMass", clist_MC3);
        hInvMass_MC4 = BSTHnSparseHelper::Load("hInvMass", clist_MC4);
        hInvMass_MC5 = BSTHnSparseHelper::Load("hInvMass", clist_MC5);
        hInvMass_MC6 = BSTHnSparseHelper::Load("hInvMass", clist_MC6);
        hInvMass_MC7 = BSTHnSparseHelper::Load("hInvMass", clist_MC7);
        hInvMass_MC_MB2 = BSTHnSparseHelper::Load("hInvMass", clist_MC2);
        hInvMass_MC_MB3 = BSTHnSparseHelper::Load("hInvMass", clist_MC3);
        hInvMass_MC_MB4 = BSTHnSparseHelper::Load("hInvMass", clist_MC4);
        hInvMass_MC_MB5 = BSTHnSparseHelper::Load("hInvMass", clist_MC5);
        hInvMass_MC_MB6 = BSTHnSparseHelper::Load("hInvMass", clist_MC6);
        hInvMass_MC_MB7 = BSTHnSparseHelper::Load("hInvMass", clist_MC7);

        hInvMass_MC2.SetBin("Pt", ptbin);
        hInvMass_MC3.SetBin("Pt", ptbin);
        hInvMass_MC4.SetBin("Pt", ptbin);
        hInvMass_MC5.SetBin("Pt", ptbin);
        hInvMass_MC6.SetBin("Pt", ptbin);
        hInvMass_MC7.SetBin("Pt", ptbin);
        hInvMass_MC_MB2.SetBin("Pt", ptbin);
        hInvMass_MC_MB3.SetBin("Pt", ptbin);
        hInvMass_MC_MB4.SetBin("Pt", ptbin);
        hInvMass_MC_MB5.SetBin("Pt", ptbin);
        hInvMass_MC_MB6.SetBin("Pt", ptbin);
        hInvMass_MC_MB7.SetBin("Pt", ptbin);

        hInvMass_MC2.SetBin("Cent", centbin);
        hInvMass_MC3.SetBin("Cent", centbin);
        hInvMass_MC4.SetBin("Cent", centbin);
        hInvMass_MC5.SetBin("Cent", centbin);
        hInvMass_MC6.SetBin("Cent", centbin);
        hInvMass_MC7.SetBin("Cent", centbin);
        hInvMass_MC_MB2.SetBin(
            "Cent", fullcentbin_Effi);  // for MC reconstruction Efficiency, we
                                        // use full bin
        hInvMass_MC_MB3.SetBin(
            "Cent", fullcentbin_Effi);  // for MC reconstruction Efficiency, we
                                        // use full bin
        hInvMass_MC_MB4.SetBin(
            "Cent", fullcentbin_Effi);  // for MC reconstruction Efficiency, we
                                        // use full bin
        hInvMass_MC_MB5.SetBin(
            "Cent", fullcentbin_Effi);  // for MC reconstruction Efficiency, we
                                        // use full bin
        hInvMass_MC_MB6.SetBin(
            "Cent", fullcentbin_Effi);  // for MC reconstruction Efficiency, we
                                        // use full bin
        hInvMass_MC_MB7.SetBin(
            "Cent", fullcentbin_Effi);  // for MC reconstruction Efficiency, we
                                        // use full bin
    }
    

    cout << "DATA AXES " << endl;
    hInvMass.PrintAxis("all");
    cout << "RSN MC AXES " << endl;
    hInvMass_MC.PrintAxis("all");
    cout << "GEN MC AXES " << endl;
    hInvMass_MC_General.PrintAxis("all");

    // output files
    TString savefile =
        Form("AnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s%i.root", sys,
             multi_start, multi_end, inputOptions, OptionNumber);
    TFile* output = new TFile(savefile.Data(), "RECREATE");

    // Basic QA Plots
    // ------------------------------------------------------------------------------------------
    TH1D* hNumberofEvent = (TH1D*)clist->FindObject(
        "fNormalisationHist");  // N of Event through event cuts
    TH1D* hEventCuts = (TH1D*)clist->FindObject(
        "fCutStats");  // N of Event through event cuts
    TList* QAlist = (TList*)clist->FindObject("EventQA");
    TH1D* hMultQA = (TH1D*)QAlist->FindObject(
        "hMult_QA");  // Multiplicty distribution after all event cuts

    hNumberofEvent->Write("hNumberofEvent");
    hEventCuts->Write("hCutStats");

    // auto hV0M = BSTHnSparseHelper::Load("hV0MSignal", clist);
    // hV0M.SetBin("Cent", fullcentbin_Effi);
    // hV0M.SetBin("SPDNtrk", {0,4000});
    // auto hV0Mbefore = hV0M.GetTH1("V0M", 3, {1, 1, -1, 1});
    // hV0Mbefore->Write("hV0MBefore");
    // auto hV0Mafter = hV0M.GetTH1("V0M", 3, {2, 1, -1, 1});
    // hV0Mafter->Write("hV0MAfter");

    double eventfraction = 0.;
    if (!fHM)
        eventfraction =
            hMultQA->Integral(hMultQA->GetXaxis()->FindBin(multi_start),
                              hMultQA->GetXaxis()->FindBin(multi_end) - 1) /
            hMultQA->Integral(hMultQA->GetXaxis()->FindBin(0.),
                              hMultQA->GetXaxis()->FindBin(100) - 1);
    else
        eventfraction =
            hMultQA->Integral(hMultQA->GetXaxis()->FindBin(multi_start),
                              hMultQA->GetXaxis()->FindBin(multi_end) - 1) /
            hMultQA->Integral(hMultQA->GetXaxis()->FindBin(0.),
                              hMultQA->GetXaxis()->FindBin(0.1) - 1);

    hMultQA->Rebin(10);
    hMultQA->GetXaxis()->SetTitle("Multiplicity Percentile (%)");
    hMultQA->GetYaxis()->SetTitle("# of Event");
    hMultQA->Write("hMultQA");
    // Trigger Efficiency
    // --------------------------------------------------------------------------------------
    auto hInvMass_MC_General_Trig =
        BSTHnSparseHelper::Load("htriggered_CINT7", clist_MC_General);
    auto htrue_cent =
        hInvMass_MC_General_Trig.GetTH1("true", 1, {1, -1, -1}); // MC True INEL>0
    auto hReco_cent =
        hInvMass_MC_General_Trig.GetTH1("Reco", 1, {2, -1, -1}); // True INEL>0 && Triggered
    auto hVtx_cent =
        hInvMass_MC_General_Trig.GetTH1("Vtx", 1, {3, -1, -1}); // True INEL>0 && triggered && good vtx  
    double sumtrue = htrue_cent->Integral(htrue_cent->GetXaxis()->FindBin(multi_start+0.0001),
            htrue_cent->GetXaxis()->FindBin(multi_end-0.0001));
    double sumreco = hReco_cent->Integral(hReco_cent->GetXaxis()->FindBin(multi_start+0.0001),
            hReco_cent->GetXaxis()->FindBin(multi_end-0.0001));
    double sumvtx = hVtx_cent->Integral(hVtx_cent->GetXaxis()->FindBin(multi_start+0.0001),
            hVtx_cent->GetXaxis()->FindBin(multi_end-0.0001));

    triggereffi = sumreco / sumtrue;
    triggereffi_e = sqrt(triggereffi*(1-triggereffi)/sumtrue);
    //triggereffi_e = sqrt(pow(sqrt(sumreco)/sumtrue,2)+pow(sqrt(sumtrue)*sumreco/pow(sumtrue,2),2));
    Double1D tempbin = {multi_start, multi_end};
    if(fAA){
        triggereffi = 1;
        triggereffi_e = zero;
    }
    TH1D* htriggereffi =
        new TH1D("TrigEffi", "Trigger Efficiency", 1, &tempbin[0]);
    htriggereffi->SetBinContent(1, triggereffi);
    htriggereffi->SetBinError(1, triggereffi_e);
    htriggereffi->Write("hTriggerEffi");

    // Vertexer Loss Correction
    // -----------------------------------------------------------------------------------------
    vertexeffi = sumvtx/sumreco;
    vertexeffi_e = sqrt(vertexeffi*(1-vertexeffi)/sumreco);
    //vertexeffi = 1;
    //vertexeffi_e = 1e-10;
    //vertexeffi_e = sqrt(pow(sqrt(sumvtx)/sumreco,2)+pow(sqrt(sumreco)*sumvtx/pow(sumreco,2),2));
    
    TH1D* hvertexeffi =
        new TH1D("hvertexeffi", "Vertex Efficiency", 1, &tempbin[0]);
    hvertexeffi->SetBinContent(1, vertexeffi);
    hvertexeffi->SetBinError(1, vertexeffi_e);
    hvertexeffi->Write("hVertexEffi");
    cout << "number of event after pileupcut: " << sumreco << endl;
    cout << "number of event after vertex selection: " << sumvtx << endl;
    cout << "trigeffi : " << triggereffi << ", vertexeffi: " << vertexeffi << endl;

    if(isVertexCutEnd){
        output->Close();
        gSystem->Exit(1);
    }
    // Signal with Bkg
    // -----------------------------------------------------------------------------------------
    for (auto j : hInvMass.BinRange("Pt")) {
        int n = j - 1;
        // SigBkg
        // ---------------------------------------------------------------------------------------------
        if(n < 1){
            normleft = 1;
            normright = 0;

        }
        else{
            normleft = 1;
            normright = 1;
        }
        auto hSig =
            hInvMass.GetTH1("Signal", 4, {0, 1, 1, j, -1});  // -1 => inv mass
        auto hBkg = hInvMass.GetTH1("Bkg", 4,
                                    {0, bkgtype, 1, j, -1});  // -1 => inv mass
        auto hBkg_norm = GetNorBkg(hSig, hBkg, normleft, normright);
        
        // rebin
        hSig->Rebin(rebin);
        hBkg_norm->Rebin(rebin);
        
        // Axis
        hSig->GetXaxis()->SetTitle("#it{M}_{inv.}(#pi#Xi) (GeV/#it{c}^{2})");
        hSig->GetYaxis()->SetTitle(
            Form("Counts / (%.1f MeV/#it{c}^{2})", 1000 * hSig->GetBinWidth(1)));
        hSig->GetYaxis()->SetLabelSize(0.05);
        hSig->GetYaxis()->SetTitleSize(0.05);
        hSig->GetXaxis()->SetLabelSize(0.05);
        hSig->GetXaxis()->SetTitleSize(0.05);

        // Draw range
        hSig->GetXaxis()->SetRangeUser(DrawRange[0], DrawRange[1]);
        hBkg_norm->GetXaxis()->SetRangeUser(DrawRange[0], DrawRange[1]);

        // Make Norm region plot
        TH1D* temp = (TH1D*)hBkg_norm->Clone();
        for (int k = 0; k < temp->GetNbinsX(); k++) {
            auto check = temp->GetBinCenter(k);
            if (check <= NormalizeRange_L[0])
                temp->SetBinContent(k, 0);
            if (!normleft && (NormalizeRange_L[0] <= check &&
                              check <= NormalizeRange_L[1]))
                temp->SetBinContent(k, 0);
            if (NormalizeRange_L[1] <= check &&
                check <= NormalizeRange_R[0])
                temp->SetBinContent(k, 0);
            if (!normright && (NormalizeRange_R[0] <= check &&
                               check <= NormalizeRange_R[1]))
                temp->SetBinContent(k, 0);
            if (NormalizeRange_R[1] <= check)
                temp->SetBinContent(k, 0);
        }
        temp->SetFillColorAlpha(kRed, 0.05);

        // Draw
        hSig->Write(Form("hSignalOnly_%i", n));
        hBkg_norm->Write(Form("hBkgOnly_%i", n));
        temp->Write(Form("hBkgNorm_%i", n));

        // MC data
        // --------------------------------------------------------------------------------------------
        // After Trigger selection Gen MC
        auto hTrueInput_Gen = // INELg0|vz<10
            hInvMass_MC_General.GetTH1("trueinput", 4, {0, 8, 1, j, -1});
        // After All event cut Gen MC
        auto hInput_Gen = // after all event cut
            hInvMass_MC_General.GetTH1("input", 4, {0, 5, 1, j, -1});
        // True After All event cut
        auto hInput = hInvMass_MC_MB.GetTH1("input", 4, {0, 6, 1, j, -1});
        // My reconstruction in All event cut
        auto hReco = hInvMass_MC_MB.GetTH1("recon", 4, {0, 4, 1, j, -1});

        
        //Further Rsn MC
        if (useFurtherMC) {
            hInput->Add(hInvMass_MC_MB2.GetTH1("input", 4, {1, 6, 1, j, -1}));
            hInput->Add(hInvMass_MC_MB3.GetTH1("input", 4, {1, 6, 1, j, -1}));
            hInput->Add(hInvMass_MC_MB4.GetTH1("input", 4, {1, 6, 1, j, -1}));
            hInput->Add(hInvMass_MC_MB5.GetTH1("input", 4, {1, 6, 1, j, -1}));
            hInput->Add(hInvMass_MC_MB6.GetTH1("input", 4, {1, 6, 1, j, -1}));
            hInput->Add(hInvMass_MC_MB7.GetTH1("input", 4, {1, 6, 1, j, -1}));
            hReco->Add(hInvMass_MC_MB2.GetTH1("recon", 4, {sys, 4, 1, j, -1}));
            hReco->Add(hInvMass_MC_MB3.GetTH1("recon", 4, {sys, 4, 1, j, -1}));
            hReco->Add(hInvMass_MC_MB4.GetTH1("recon", 4, {sys, 4, 1, j, -1}));
            hReco->Add(hInvMass_MC_MB5.GetTH1("recon", 4, {sys, 4, 1, j, -1}));
            hReco->Add(hInvMass_MC_MB6.GetTH1("recon", 4, {sys, 4, 1, j, -1}));
            hReco->Add(hInvMass_MC_MB7.GetTH1("recon", 4, {sys, 4, 1, j, -1}));
        }
        Double_t Input_number_Gen = hInput_Gen->Integral(
            hInput_Gen->GetXaxis()->FindBin(IntegralRangeMC[0]),
            hInput_Gen->GetXaxis()->FindBin(IntegralRangeMC[1]) - 1);

        Double_t True_Input_number_Gen = hTrueInput_Gen->Integral(
            hTrueInput_Gen->GetXaxis()->FindBin(IntegralRangeMC[0]),
            hTrueInput_Gen->GetXaxis()->FindBin(IntegralRangeMC[1]) - 1);
        ;

        Double_t Input_number = hInput->Integral(
            hInput->GetXaxis()->FindBin(IntegralRangeMC[0]),
            hInput->GetXaxis()->FindBin(IntegralRangeMC[1]) - 1);

        Double_t Reco_number =
            hReco->Integral(hReco->GetXaxis()->FindBin(IntegralRangeMC[0]),
                            hReco->GetXaxis()->FindBin(IntegralRangeMC[1]) - 1);

        /*
        // old way
        Double_t Eff_e = sqrt(
            pow(sqrt(Reco_number) / Input_number, 2) +
            pow(sqrt(Input_number) * Reco_number / pow(Input_number, 2), 2));
        */
        Double_t RecEff = Reco_number/Input_number;
        Double_t Eff_e = sqrt(RecEff*(1-RecEff)/Input_number);
        cout << "old rec_err: " << sqrt(
            pow(sqrt(Reco_number) / Input_number, 2) +
            pow(sqrt(Input_number) * Reco_number / pow(Input_number, 2), 2))
            << ", new rec_err: " << Eff_e << endl;
        Double_t ecutratio = True_Input_number_Gen / Input_number_Gen;
        Double_t ecuteffi = Input_number_Gen / True_Input_number_Gen; // same, but efficiency form for the error.
        Double_t ecuteffi_e = sqrt(ecuteffi*(abs(1-ecuteffi))/True_Input_number_Gen); // proper error propagation
        Double_t ecutratio_e = ecutratio * (ecuteffi_e/ecuteffi); // use the percentage of error here.

        cout << "True input: " << True_Input_number_Gen << ", after event cuts: " << Input_number_Gen << endl;
        cout << "SL: " << ecuteffi << ", error: " << ecuteffi_e << ", fraction: " << ecuteffi_e/ecuteffi << endl;
        cout << "old error: " << sqrt(pow(sqrt(True_Input_number_Gen)/Input_number_Gen,2)+pow(sqrt(Input_number_Gen)*True_Input_number_Gen/pow(Input_number_Gen,2),2))
             << ", new error: " << ecutratio_e << endl;
        //Double_t ecutratio_e = sqrt(pow(sqrt(True_Input_number_Gen)/Input_number_Gen,2)+pow(sqrt(Input_number_Gen)*True_Input_number_Gen/pow(Input_number_Gen,2),2));
        if(fHM){
            ecutratio = 1.0;
            ecutratio_e = zero;
        }
        EventCutRatio.push_back(ecutratio);
        EventCutRatio_e.push_back(ecutratio_e);
        Efficiency.push_back(RecEff);
        Efficiency_e.push_back(Eff_e);

        // MC Fit
        hReco->Rebin(rebin);
        hReco->GetXaxis()->SetRangeUser(DrawRange[0], DrawRange[1]);
        hReco->Write(Form("hMCRecon_%i", n));

        TF1* MCfitresult = VoigtFit(hReco, formula_MC, peakRange, par, fix_MC);

        fitmean_MC.push_back(MCfitresult->GetParameter(vp + 1));
        fitmean_MC_err.push_back(MCfitresult->GetParError(vp + 1));
        fitsigma_MC.push_back(MCfitresult->GetParameter(vp + 2));
        fitsigma_MC_err.push_back(MCfitresult->GetParError(vp + 2));
        fitgamma_MC.push_back(MCfitresult->GetParameter(vp + 3));
        fitgamma_MC_err.push_back(MCfitresult->GetParError(vp + 3));

        MCfitresult->Write(Form("fMCFit_%i", n));

        // Fit
        // ------------------------------------------------------------------------------------------------
        TH1D* Sigfit = (TH1D*)hSig->Clone();
        TH1D* Bkgfit = (TH1D*)hBkg_norm->Clone();
        if(!fbkgfit){
            Sigfit->Add(Bkgfit, -1);
        }

        TGaxis::SetMaxDigits(3);
        Sigfit->Write(Form("hSignalBkgSubtraction_%i", n));

        vector<double> par_fit = {0,0,0};
        vector<double> fix_fit = {0,0,0};
        int fitpars;
        if (n == 0){
            for(fitpars = 0; fitpars<fix_fit.size();fitpars++){
                par_fit[fitpars] = par_2nd[fitpars];
                fix_fit[fitpars] = fix_3rd[fitpars];
            }
        }
        else if (n == 1){
            for(fitpars = 0; fitpars<fix_fit.size();fitpars++){
                par_fit[fitpars] = par[fitpars];
                fix_fit[fitpars] = fix_2nd[fitpars];
            }
        }
        else{
            for(fitpars = 0; fitpars<fix_fit.size();fitpars++){
                par_fit[fitpars] = par[fitpars];
                fix_fit[fitpars] = fix[fitpars];
            }
        }
        par_fit[1] = fitsigma_MC[n];
        
        int k;
        TH1D* a = (TH1D*)Sigfit->Clone(Form("%s_nopeak", name));
        for (k = Sigfit->GetXaxis()->FindBin(1.000001 * peakRange[0]);
             k <= Sigfit->GetXaxis()->FindBin(0.999999 * peakRange[1]); k++) {
            a->SetBinContent(k, 0.);
            a->SetBinError(k, 0.);
        }
        ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); // better minimizer!
        TF1* fb = new TF1(Form("%s_back", name), formula.Data(), FitRange[0], FitRange[1]);
        a->Fit(fb, "RQN"); // L
        cout << "First fit: " << fb->GetParameter(0) << ", " << fb->GetParameter(1) << ", " << fb->GetParameter(2) << endl;

        fb->SetLineColor(3);
        fb->SetLineStyle(2);
        //fb->Draw("same");

        vp = fb->GetNpar();

        TF1* fp = new TF1("Voigt+bkg",
                      Form("%s+[%i]*TMath::Voigt(x-[%i],[%i],[%i])", formula.Data(),
                           vp, vp + 1, vp + 2, vp + 3),
                      FitRange[0], FitRange[1]);
        for (k = 0; k < vp; k++) {
            fp->SetParameter(k, fb->GetParameter(k));
            fp->FixParameter(k, fb->GetParameter(k));
        }
        fp->SetParameter(vp, Sigfit->GetBinContent(Sigfit->GetXaxis()->FindBin(
                                 0.5 * (peakRange[0] - peakRange[1]))) -
                                 fb->Eval(0.5 * (peakRange[0] - peakRange[1])));
        for (k = 0; k < vp; k++) {
            fp->SetParameter(vp + k + 1, par_fit[k]);
            fp->FixParameter(vp + k + 1, par_fit[k]);
        }

        fp->SetLineColor(2);
        fp->SetLineStyle(1);

        Sigfit->Fit(fp, "RQN");

        if (!fix_fit[2]) {  // release width
            fp->ReleaseParameter(vp + 3);
            fp->SetParError(vp + 3, 0.1 * fp->GetParameter(vp + 3));
            Sigfit->Fit(fp, "RQN");
        }

        if (!fix_fit[1]) {  // release resolution
            fp->ReleaseParameter(vp + 2);
            fp->SetParLimits(vp+2,0.,0.01);
            fp->SetParError(vp + 2, 0.1 * fp->GetParameter(vp + 2));
            Sigfit->Fit(fp, "RQN");
        }

        if (!fix_fit[0]) {  // release mass
            fp->ReleaseParameter(vp + 1);
            fp->SetParError(vp + 1, 0.1 * fp->GetParameter(vp + 1));
            Sigfit->Fit(fp, "RQN");
        }

        // release background constant parameter
        fp->ReleaseParameter(0);
        fp->SetParError(0, fb->GetParError(0));
        Sigfit->Fit(fp, "RQN");

        // release other background parameters
        for (k = 0; k < vp; k++) {
            fp->ReleaseParameter(k);
            fp->SetParError(k, fb->GetParError(k));
        }
        Sigfit->Fit(fp, "RQN");

        TVirtualFitter::SetMaxIterations(1000000);
        // final fit
        TFitResultPtr fitresult_ptr;
        int fitcounter = 0;
        while(1){
          fitresult_ptr= Sigfit->Fit(fp, "RIS", "", FitRange[0], FitRange[1]);
          
          // fit result
          if(n == 0) break; // 0-0.8 skip
          if(n == 10) break; // 8.8-15 skip

          if(fitresult_ptr) // only converge result
            break;
          fitcounter++;
          if(fitcounter > 10) break;
        }
        fb_after = new TF1(Form("%s_back", name), formula.Data(), 0, 99);
        
        for (k = 0; k < vp; k++) {
            fb_after->SetParameter(k, fp->GetParameter(k));
        }
        fb_after->Draw();
        fb_after->SetRange(IntegralRange[0], IntegralRange[1]);
        fb_after->SetLineColor(3);
        fb_after->SetLineStyle(2);
        fp->Draw("same");

        bkgintegral =
            1000 * fb_after->Integral(IntegralRange[0], IntegralRange[1],1e-6) / rebin;

        fp->Write(Form("fDataFitResult_%i", n));
        TF1* fonly = new TF1("Voigt", "[0]*TMath::Voigt(x-[1],[2],[3])",
                            DrawRange[0], DrawRange[1]);
        fonly->SetParameter(0, fp->GetParameter(vp));
        fonly->SetParError(0, fp->GetParError(vp));
        fonly->SetParameter(1, fp->GetParameter(vp + 1));
        fonly->SetParError(1, fp->GetParError(vp + 1));
        fonly->SetParameter(2, fp->GetParameter(vp + 2));
        fonly->SetParError(2, fp->GetParError(vp + 2));
        fonly->SetParameter(3, fp->GetParameter(vp + 3));
        fonly->SetParError(3, fp->GetParError(vp + 3));
        fonly->Write(Form("fDataFitOnlyResult_%i", n));
        fb_after->Write(Form("fDataFitBkgResult_%i", n));
        
        Double_t integral_result_fullfit = fp->Integral(IntegralRange[0], IntegralRange[1], 1.e-6);
        Double_t integral_e_fullfit = fp->IntegralError(IntegralRange[0], IntegralRange[1],fp->GetParameters(),fitresult_ptr->GetCovarianceMatrix().GetMatrixArray(), 1.e-6);
        Double_t integral_e_ratio = abs(integral_e_fullfit/integral_result_fullfit);

        cout << "full integral: " << 1000*integral_result_fullfit << ", bkg integral: " << bkgintegral << ", sig integral: " << 1000*integral_result_fullfit - bkgintegral << ", e: " << 1000*integral_e_ratio << endl;

        fitmean.push_back(fp->GetParameter(vp + 1));
        fitmean_err.push_back(fp->GetParError(vp + 1));
        fitsigma.push_back(fp->GetParameter(vp + 2));
        fitsigma_err.push_back(fp->GetParError(vp + 2));
        fitgamma.push_back(fp->GetParameter(vp + 3));
        fitgamma_err.push_back(fp->GetParError(vp + 3));
        
        double frawyield = 0.;
        double frawyield_e = 0.;
        double frawyield_bincount = 0.;
        double frawyield_bincount_e = 0.;
        double frawyield_bincount_e_ratio = 0.;
        // Bin count method
        frawyield_bincount =
            Sigfit->Integral(Sigfit->GetXaxis()->FindBin(IntegralRange[0]),
                             Sigfit->GetXaxis()->FindBin(IntegralRange[1] -
                                                         0.001 * rebin));
        for (int bin = Sigfit->GetXaxis()->FindBin(IntegralRange[0]);
             bin <=
             Sigfit->GetXaxis()->FindBin(IntegralRange[1]+0.001*rebin);
             bin++) {
            frawyield_bincount_e += pow(Sigfit->GetBinError(bin), 2);
        }

        frawyield_bincount_e_ratio = abs(sqrt(frawyield_bincount_e)/frawyield_bincount);
        cout << "Bincount error ratio: " << frawyield_bincount_e_ratio <<  endl;

        // lost yield correction (additive)
        double lefttail =
            1000 * fonly->Integral(IntegralRangeMC[0], IntegralRange[0]) /
            rebin;
        double righttail =
            1000 * fonly->Integral(IntegralRange[1], IntegralRangeMC[1]) /
            rebin;
        
        frawyield_bincount = frawyield_bincount + lefttail + righttail - bkgintegral;
        frawyield_bincount_e = frawyield_bincount_e 
            + pow((lefttail + righttail)*(fonly->GetParError(0)/fonly->GetParameter(0)),2);
    
        // Integral from fit function
        frawyield = 1000 * // 1GeV/c2 = 1000 MeV/c2
                    fonly->Integral(IntegralRangeMC[0], IntegralRangeMC[1],1e-6) /
                    rebin;
        
        //frawyield_e = frawyield*(fonly->GetParError(0)/fonly->GetParameter(0)); // old method 
        //frawyield_e = frawyield*integral_e_ratio; // New method
        frawyield_e = frawyield_bincount_e_ratio*frawyield; // bincount error method
            
        if (bincount){
            RawYield.push_back(frawyield_bincount);
            RawYield_err.push_back(sqrt(frawyield_bincount_e));  
        } 
        else {
            RawYield.push_back(frawyield);
            RawYield_err.push_back(frawyield_e);
        }
        
        // Old error: Fit parameter error ratio from normalisation factor of voigt fit [0].
        // New error: Error ratio of full fit(voigt+bkg)
        // Bincount error: error sum of each inv.mass bin ratio to bincount yield
        cout << "Fit status: " << fitresult_ptr << ", yield: " << RawYield[n] << ", old error: " << frawyield*(fonly->GetParError(0)/fonly->GetParameter(0))
                 << ", new error: " << frawyield*integral_e_ratio << ", bincount error: " << RawYield_err[n] << endl;
        
    }
    //-----------------------------------------------------------------------------------------------------
    // MC Efficiency
    auto hEfficiency =
        MakeHistfromArray("MC Efficiency", Efficiency, Efficiency_e, ptbin);
    hEfficiency->SetBinContent(1, 1e-10);
    hEfficiency->SetBinError(1, 1e-10);
    hEfficiency->SetMinimum(4e-4);
    hEfficiency->SetMaximum(4e-1);
    hEfficiency->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hEfficiency->GetYaxis()->SetTitle("Acceptance x Efficiency x BR");
    hEfficiency->Write("hMCReconEffi");

    // ----------------------------------------------------------------------------------------------------
    // Signal loss
    auto hsigloss =
        MakeHistfromArray("Signal loss", EventCutRatio, EventCutRatio_e, ptbin);
    hsigloss->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hsigloss->GetYaxis()->SetTitle("Signal Loss");
    hsigloss->SetMinimum(0);
    hsigloss->Write("hMCSigLoss");

    // Spectra
    // --------------------------------------------------------------------------------------------
    auto hRawYields =
        MakeHistfromArray("Raw Yields", RawYield, RawYield_err, ptbin);
    hRawYields->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRawYields->GetYaxis()->SetTitle("Raw yield");
    hRawYields->Write("hDataRawYield");

    auto hXispectrum = MakeHistfromArray("hXispectrum", ptzero, ptzero, ptbin);
    hXispectrum->SetMarkerStyle(20);
    hXispectrum->SetMarkerSize(1.0);
    hXispectrum->SetTitle("#Xi(1530) spectrum pp #sqrt{s}= 13 TeV");
    hXispectrum->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
    hXispectrum->GetYaxis()->SetTitle(
        "1/N_{event}d^{2}N/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}");

    Double_t value = 0.;
    Double_t base = 0.;
    Double_t err = 0.;
    // pp 13TeV f_norm, which was used in LHC15 period data sets.
    // Not uisng in my analysis. the value itself is comparable.
    //triggereffi = GetNormalisationFactor(multi_start, multi_end)[0];
    //triggereffi_e = GetNormalisationFactor(multi_start, multi_end)[1];
    if(fpPb){
        triggereffi = 1;
        triggereffi_e = 0;
    }
    hXispectrum->SetMinimum(1.0e-3);
    hXispectrum->SetMaximum(3.0);
    //double nOfEventMultibin = hNumberofEvent->GetBinContent(9) * eventfraction;
    double nOfEventMultibin = zero;
    double eventMultiratio = zero;
    if(!fHM) 
        eventMultiratio = (multi_end - multi_start) / 100;
    else
        eventMultiratio = (multi_end - multi_start) / 0.1;
    //nOfEventMultibin = hMultQA->GetBinContent(1);

    nOfEventMultibin = hMultQA->Integral(hMultQA->GetXaxis()->FindBin(multi_start),
                               hMultQA->GetXaxis()->FindBin(multi_end)-1);
    for (int ipTbin = 1; ipTbin < (int)ptbin.size() - 2; ipTbin++) {
        // for debuging
        cout << "RawYield[" << ipTbin << "]: " << RawYield[ipTbin]
             << ", RawYield_err[" << ipTbin << "]: " << RawYield_err[ipTbin]
             << ", Efficiency[" << ipTbin << "]: " << Efficiency[ipTbin]
             << ", EventCutRatio[" << ipTbin << "]: " << EventCutRatio[ipTbin]
             << ", centbin: " << centbin[0]
             << " - " << centbin[1]
             << ", # of Events: " << hNumberofEvent->GetBinContent(5)
             << ", # of Events in mutlbin: " << nOfEventMultibin
             << "(" << eventMultiratio*100
             << "%) , pt_points_e[" << ipTbin << "]: " << pt_points_e[ipTbin] 
             << ", Vertex lost correction: " << vertexeffi
             << ", Trigger Efficiency: " << triggereffi << endl;
        
        // Normalised & Corrected Value -----------------------------
        value = RawYield[ipTbin];
        base = Efficiency[ipTbin] * pt_points_e[ipTbin] * 1 * // Rapidity(0.5 + 0.5)
               nOfEventMultibin;
        if(fpPb) base *= 0.5; // rapidity
        value /= base;

        if(!fAA) value *= triggereffi * vertexeffi;
        value *= EventCutRatio[ipTbin]; // <- Corrected and normalised Xi(1530) + c.c.

        // STATISTICAL ERROR ----------------------------------------
        if(fAA) err = pow(RawYield_err[ipTbin] * EventCutRatio[ipTbin] / base, 2);
        else
            err = pow( 
                (RawYield_err[ipTbin]* triggereffi * vertexeffi * EventCutRatio[ipTbin]) / base
                , 2);
        // MC Efficiency error
        err += pow(
            value * (Efficiency_e[ipTbin] / Efficiency[ipTbin]), 
            2);
        // Signal loss error
        err += pow(value * (EventCutRatio_e[ipTbin] / EventCutRatio[ipTbin]), 2);

        if(!fAA){ // for HM(AA), trigger efficiency is 1, vertex loss is 1.
            // Trigger efficiency error
            err += pow(value * (triggereffi_e/triggereffi), 2);
            // Vertex loss correction error
            err += pow(value * (vertexeffi_e/vertexeffi),2); 
        }
        // Number of event error -> DON'T NEED TO USE!
        // since this error(statistical error) must be included in RawYield_err
        /*
        err += pow(
            value * (sqrt(nOfEventMultibin) / nOfEventMultibin), // Number of event
            2);
        */
        err = sqrt(err);
        // ----------------------------------------------------------

        value /= 2.;  // Xi(1530)0 + c.c.
        err /= 2.;

        cout << "Corrected yield: " << value << endl;
        cout << "Error: " << err << endl;
        hXispectrum->SetBinContent(ipTbin + 1, value);
        hXispectrum->SetBinError(ipTbin + 1, err);
    }
    //hXispectrum->SetBinContent(ptbin.size(), zero);
    //hXispectrum->SetBinError(ptbin.size(), zero);
    hXispectrum->Write("hXiSpectrum");

    TCanvas* cSpectra = new TCanvas("cSpectra", "", w, h);
    cSpectra->Draw();
    cSpectra->SetTickx();
    cSpectra->SetLogy(true);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);

    TH1D* temp_spectra = (TH1D*)hXispectrum->Clone();

    TF1* myLevy = new TF1("myLevy", myLevyPt, 0, 8.5, 3);
    myLevy->SetParName(0, "dN/dy");
    myLevy->SetParName(1, "C");
    myLevy->SetParName(2, "n");
    myLevy->SetParameter(0, .008);
    myLevy->SetParameter(1, .3);
    myLevy->SetParameter(2, 15);
    myLevy->SetParLimits(0, .001, .2);
    myLevy->SetParLimits(1, .1, 1);
    myLevy->SetParLimits(2, 1, 500);

    temp_spectra->Fit(myLevy, "IME", "", 0.8, 8.8);
    auto legend_Spetra_fit = new TLegend(0.57, 0.70, 0.85, 0.83);
    legend_Spetra_fit->SetFillStyle(0);
    legend_Spetra_fit->AddEntry(temp_spectra,
                                "#Xi(1530)^{0} yield 13TeV (0-100)", "PEL");
    legend_Spetra_fit->AddEntry(myLevy, "Levy fit", "L");
    t2->DrawLatex(0.18, 0.92, "#bf{#Xi(1530)^{0} #rightarrow #Xi + #pi}");
    t2->DrawLatex(0.75, 0.92, "#bf{pp 13 TeV}");
    t->DrawLatex(
        0.59, 0.67,
        Form("#bf{dN/dy: %.3f #pm %.3f (x10^{-3})}",
             myLevy->GetParameter(0) * 1e3, myLevy->GetParError(0) * 1e3));
    t->DrawLatex(
        0.59, 0.62,
        Form("#bf{C: %.2f #pm %.2f (x10^{-3})}", myLevy->GetParameter(1) * 1e3,
             myLevy->GetParError(1) * 1e3));
    t->DrawLatex(0.59, 0.57,
                 Form("#bf{n: %.2f #pm %.2f}", myLevy->GetParameter(2),
                      myLevy->GetParError(2)));
    legend_Spetra_fit->Draw();

    TF1* myLevy_full = new TF1("myLevy", myLevyPt, 0, 8.5, 3);
    myLevy_full->SetParName(0, "dN/dy");
    myLevy_full->SetParName(1, "C");
    myLevy_full->SetParName(2, "n");
    myLevy_full->SetParameter(0, myLevy->GetParameter(0));
    myLevy_full->SetParameter(1, myLevy->GetParameter(1));
    myLevy_full->SetParameter(2, myLevy->GetParameter(2));
    myLevy_full->Draw("same");
    myLevy_full->Write("Levy_full");
    cSpectra->Write("cFitSpectra");
    output->Close();
}

TH1D* GetNorBkg(TH1D* hSig, TH1D* hBkg, int isnormleft, int isnormright) {
    Double_t normalization_data = 0.;
    Double_t normalization_mixed = 0.;
    if (isnormleft == 1) {
        normalization_data +=
            hSig->Integral(hSig->GetXaxis()->FindBin(NormalizeRange_L[0]),
                           hSig->GetXaxis()->FindBin(NormalizeRange_L[1]) - 1);
        normalization_mixed +=
            hBkg->Integral(hBkg->GetXaxis()->FindBin(NormalizeRange_L[0]),
                           hBkg->GetXaxis()->FindBin(NormalizeRange_L[1]) - 1);
    }
    if (isnormright == 1) {
        normalization_data +=
            hSig->Integral(hSig->GetXaxis()->FindBin(NormalizeRange_R[0]),
                           hSig->GetXaxis()->FindBin(NormalizeRange_R[1]) - 1);
        normalization_mixed +=
            hBkg->Integral(hBkg->GetXaxis()->FindBin(NormalizeRange_R[0]),
                           hBkg->GetXaxis()->FindBin(NormalizeRange_R[1]) - 1);
    }
    hBkg->Scale(normalization_data / normalization_mixed);
    hBkg->SetLineColor(kRed);
    hBkg->SetMarkerColor(kRed);

    return hBkg;
}
TF1* VoigtFit(TH1D* hSig,
              TString bkgformula,
              double* peakRange,
              double* par,
              double* fix) {
    int j;
    TH1D* a = (TH1D*)hSig->Clone(Form("%s_nopeak", name));
    for (j = hSig->GetXaxis()->FindBin(1.000001 * peakRange[0]);
         j <= hSig->GetXaxis()->FindBin(0.999999 * peakRange[1]); j++) {
        a->SetBinContent(j, 0.);
        a->SetBinError(j, 0.);
    }
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); // better minimizer!
    //ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6); // fit should converge!

    TF1* fb = new TF1(Form("%s_back", name), bkgformula.Data(), FitRange[0], FitRange[1]);
    a->Fit(fb, "RQN"); // L

    fb->SetLineColor(3);
    fb->SetLineStyle(2);
    //fb->Draw("same");

    vp = fb->GetNpar();

    TF1* fp = new TF1("Voigt+bkg",
                      Form("%s+[%i]*TMath::Voigt(x-[%i],[%i],[%i])", bkgformula.Data(),
                           vp, vp + 1, vp + 2, vp + 3),
                      FitRange[0], FitRange[1]);
    for (j = 0; j < vp; j++) {
        fp->SetParameter(j, fb->GetParameter(j));
        fp->FixParameter(j, fb->GetParameter(j));
    }
    fp->SetParameter(vp, hSig->GetBinContent(hSig->GetXaxis()->FindBin(
                             0.5 * (peakRange[0] - peakRange[1]))) -
                             fb->Eval(0.5 * (peakRange[0] - peakRange[1])));
    for (j = 0; j < 3; j++) {
        fp->SetParameter(vp + j + 1, par[j]);
        fp->FixParameter(vp + j + 1, par[j]);
    }

    fp->SetLineColor(2);
    fp->SetLineStyle(1);

    hSig->Fit(fp, "RQN");

    if (!fix[2]) {  // release width
        fp->ReleaseParameter(vp + 3);
        fp->SetParError(vp + 3, 0.1 * fp->GetParameter(vp + 3));
        hSig->Fit(fp, "RQN");
    }

    if (!fix[1]) {  // release resolution
        fp->ReleaseParameter(vp + 2);
        fp->SetParLimits(vp+2,0.,0.01);
        fp->SetParError(vp + 2, 0.1 * fp->GetParameter(vp + 2));
        hSig->Fit(fp, "RQN");
    }

    if (!fix[0]) {  // release mass
        fp->ReleaseParameter(vp + 1);
        fp->SetParError(vp + 1, 0.1 * fp->GetParameter(vp + 1));
        hSig->Fit(fp, "RQN");
    }

    // release background constant parameter
    fp->ReleaseParameter(0);
    fp->SetParError(0, fb->GetParError(0));
    hSig->Fit(fp, "RQN");

    // release other background parameters
    for (j = 1; j < vp; j++) {
        fp->ReleaseParameter(j);
        fp->SetParError(j, fb->GetParError(j));
    }
    hSig->Fit(fp, "RQN");

    // final fit
    hSig->Fit(fp, "RI");

    fb_after = new TF1(Form("%s_back", name), bkgformula.Data(), 0, 99);
    for (j = 0; j < vp; j++) {
        fb_after->SetParameter(j, fp->GetParameter(j));
    }
    fb_after->Draw();
    fb_after->SetRange(FitRange[0], FitRange[1]);
    fb_after->SetLineColor(3);
    fb_after->SetLineStyle(2);
    fp->Draw("same");

    // signal only for integral
    TF1* fonly = new TF1("Voigt", "[0]*TMath::Voigt(x-[1],[2],[3])",
                         0, 99);
    fonly->SetParameter(0, fp->GetParameter(vp));
    fonly->SetParameter(1, fp->GetParameter(vp + 1));
    fonly->SetParameter(2, fp->GetParameter(vp + 2));
    fonly->SetParameter(3, fp->GetParameter(vp + 3));
    fonly->SetLineColor(4);
    fonly->SetLineStyle(2);
    // fonly->Draw("same");

    bkgintegral =
        1000 * fb_after->Integral(IntegralRange[0], IntegralRange[1],1e-6) / rebin;

    return fp;
}

TH1D* MakeHistfromArray(char const* name,
                        Double1D dArray,
                        Double1D eArray,
                        Double1D ptbin,
                        const char* foption) {
    TString option = foption;
    double* ptbin_array = &ptbin[0];

    TH1D* htemp = new TH1D(Form("%s", name), "", ptbin.size() - 1, ptbin_array);
    for (int i = 0; i < (int)dArray.size(); i++) {
        htemp->SetBinContent(i + 1, dArray[i]);
        if (!option.Contains("NOERROR"))
            htemp->SetBinError(i + 1, eArray[i]);
    }
    return htemp;
}
double myLevyPt(Double_t* x, Double_t* par) {
    double lMass = 1.5318;  // Xi mass

    Double_t ldNdy = par[0];          // dN/dy
    Double_t l2pi = 2 * TMath::Pi();  // 2pi
    Double_t lTemp = par[1];          // Temperature
    Double_t lPower = par[2];         // power=n

    Double_t lBigCoef =
        ((lPower - 1) * (lPower - 2)) /
        (l2pi * lPower * lTemp * (lPower * lTemp + lMass * (lPower - 2)));
    Double_t lInPower = 1 + (TMath::Sqrt(x[0] * x[0] + lMass * lMass) - lMass) /
                                (lPower * lTemp);

    return l2pi * ldNdy * x[0] * lBigCoef *
           TMath::Power(lInPower, (-1) * lPower);
}
vector<double> GetNormalisationFactor(double multi_start, double multi_end){
    // Return event normalisation factor with give Multiplicity bin.
    // return {value, err}
    
    vector<double> returnarray;

    //--F_norm--
    // Ref: Centrality-like Dependence of the pseudorapidity distribution 
    //      of charged tracks in pp collisions at sqrt(s) = 13 TeV
    //      https://alice-notes.web.cern.ch/node/510
    // This value has been used for K*, Phi analysis.

    vector<double> factor_multibin = 
    {0,     1,     5,    10,    15,    20,    30,   40,   50,   70, 100};
    vector<double> factor = 
    {0, 0.9825, 0.999, 0.9982, 0.9978, 0.9969, 0.9937, 0.9870,  0.9761, 0.9491, 0.8729};
    vector<double> factor_e = {0,  0.0169,  0.0012,  0.0010,  0.0012,  0.0015, 0.0014, 0.0022, 0.0047, 0.0080, 0.0201};
    
    // input must be in the multiplicity range
    if(std::find(factor_multibin.begin(), factor_multibin.end(), multi_start) == end(factor_multibin))
        return {99,99};
    if(std::find(factor_multibin.begin(), factor_multibin.end(), multi_end) == end(factor_multibin))
        return {99,99};

    // special cases
    if((multi_start == 0) && (multi_end == 100)){
        returnarray = {0.9468, 0.0081};
    }
        else{
        // Common case
        // Value
        vector<double>::iterator itr_left = find(factor_multibin.begin(),
                                factor_multibin.end(),
                                multi_start);
        vector<double>::iterator itr_right = find(factor_multibin.begin(),
                                factor_multibin.end(),
                                multi_end);
        int left = distance(factor_multibin.begin(), itr_left);
        int right = distance(factor_multibin.begin(), itr_right);

        int gap = right - left;

        double result = 0.;
        for(int i = 1; i < gap+1; i++)
            result += factor[i+left]*(factor_multibin[i+left] - factor_multibin[i+left-1]);
            
        result /= (multi_end - multi_start);
        returnarray.push_back(result);

        // Error
        double error = 0.;
        for(int i = 1; i < gap+1; i++)
            error += pow( factor_e[i+left], 2); 
            
        error = sqrt(error);
        returnarray.push_back(error);
    }

    return returnarray;
}
TFile* LoadXi1530Results(TString name, TString runnum) {
    if (!runnum.IsNull())
        runnum = "_" + runnum;
    //auto fname = "/alice/home/blim/postprocessing/data/AnalysisResults_" + name + runnum + ".root";
    auto fname = "AnalysisResults_" + name + runnum + ".root";
    return LoadRoot(fname);
}
//__________________________________________________________
TObject* LoadXi1530ResultList(TFile* fh, TString clistname) {
    auto l = fh->Get(clistname);
    if (!l)
        ErrorExit("No list " + clistname);
    return l;
}
//__________________________________________________________
TObject* LoadXi1530ResultList(TString fname,
                              TString clistname,
                              TString runnum) {
    auto f = LoadXi1530Results(fname, runnum);
    return LoadXi1530ResultList(f, clistname);
}