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
        // cout << "bin: " << b << ", " << b << endl;
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
        // cout << "name: " << name << ", b: " << b[0] << ", " << b[1] << ",
        // maxnbin: " << maxnbin << endl;
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
    //cout << "Axis: " << iaxis << endl;
    if (bins.size() == 0)
        LoadBinFromHist(iaxis);
    else {
        fCBins[iaxis] = bins;
        fCBinsB[iaxis] = {0};
        for (UInt_t ib = 0; ib < bins.size(); ib++) {
            int bin = ax->FindBin(bins[ib]);
            // TODO if( ib>0 && bin <= fCBinsB[i][ib-1] ) ErrorExit("Wrong
            // bin");
            fCBinsB[iaxis].push_back(ax->FindBin(bins[ib]));
            //cout << "bins: " << bins[ib] << ", bin: " << bin << endl;
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
    //std::cout << "axis: " << iaxis << ", nbins: " << ax->GetNbins()+1 <<
    //std::endl;
    for (int ib = 1; ib <= ax->GetNbins() + 1; ib++) {
        fCBins[iaxis].push_back(ax->GetBinLowEdge(ib));
        fCBinsB[iaxis].push_back(ib);
        //std::cout << "fCBinsB[" << iaxis << "].at(" << ib << "): " <<
        //fCBinsB[iaxis].at(ib) << std::endl;
    }
    //std::cout << "fCBinsB[" << iaxis << "].size(): " << fCBinsB[iaxis].size()
    //<< std::endl;
}
//__________________________________________________________
void BSTHnSparseHelper::LoadBinFromHist() {
    fCBins.resize(GetNdim());
    fCBinsB.resize(GetNdim());
    for (int i = 0; i < GetNdim(); i++)
        LoadBinFromHist(i);
}
//__________________________________________________________
//__________________________________________________________
//__________________________________________________________
//__________________________________________________________
//__________________________________________________________
//__________________________________________________________
// Initial values

// constants
const int POLdegree = 1;
const float inf = 1e20;
const float zero = 1e-20;

// Input files
//char const* datafile = "LHC15fi_16deghijklop_17cefgijklor_18bdefhikmnop_SYS_full";
//char const* datafile = "LHC15fi_16deghijklop_17cefgijklmor_18bghi_SYS_HM";
//char const* datafile = "LHC18befghikmnop_SYS_HM";
//char const* datafile = "LHC15fi_16deghlop_17cefijklor_SYS_HM";
// char const* datafile = "Xi1530LHC16klop_17kl_SYS_Mix";
char const* datafile = "LHC15fi_16deghijklop_17cefgijklmor_SYS_light_fixed";
char const* rsnmcfile = "Xi1530LHC18c6b_RsnMC_SYS_fixed";
//char const* rsnmcfile = "LHC18c6bc_RsnMC_SYS_rebin";
//char const* genmcfile = "Xi1530LHC16k_pass2_GenMC_trig_SYS"; // last
char const* genmcfile = "Xi1530LHC16kl_pass2_GenMC_SYS_train";
TString inputDirectory = "Xi1530MB";

// FIT SETUP -----------------
double peakRange[2] = {1.512, 1.552};
//double peakRange[2] = {1.512, 1.56};

char name[500] = "Voigt fit";
TString formula = "pol0(0)";
TString formula_MC = "pol0(0)";
double par[3] = {1.532, 0.003, 0.0091};
double fix[3] = {0, 0, 0};
double MCfix[3] = {0, 0, 1};
double fix_2nd[3] = {1, 0, 1};
double fix_3rd[3] = {1, 1, 0};
double fix_4th[3] = {0, 1, 1};
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
Double1D ptbin = {0.0, 20};
Double1D ptzero = {zero, zero, zero, zero, zero, zero,
                   zero, zero, zero, zero, zero, zero};
Double1D pt_points_e = {};

// Norm region
Int_t normleft = 0;
//Double1D NormalizeRange_L = {1.49, 1.51}; // Norminal
Double1D NormalizeRange_L = {1.48, 1.50}; // Norminal
Int_t normright = 1;
//Double1D NormalizeRange_R = {1.56, 1.58}; // Norminal
Double1D NormalizeRange_R = {1.56, 1.59}; // Norminal

// Axis Range
//Double1D DrawRange = {1.472, 1.592};
Double1D DrawRange = {1.448, 1.62};
//Double1D FitRange = {1.472, 1.592};
Double1D FitRange = {1.484, 1.592};
Double1D IntegralRange = {1.516, 1.548};
//Double1D IntegralRange = {1.4, 1.6};
Double1D IntegralRangeMC = {0, 99};
// Double1D IntegralRange = {0, 99};

// Rebin
Int_t rebin = 4;

// Colors
Int_t colors[] = {633, 810, 800, 830, 840, 840, 870, 864, 890, 617, 619};
Int1D fColors = {kBlack, kRed, kBlue, kGreen + 3, kMagenta + 2,
                 kBlack, kRed, kBlue, kGreen + 3, kMagenta + 2,
                 kBlack, kRed, kBlue, kGreen + 3, kMagenta + 2,
                 kBlack, kRed, kBlue, kGreen + 3, kMagenta + 2};
Int1D fMarkers = {24, 25, 26, 27, 28, 30, 20, 21, 22, 33, 34, 24, 25, 26, 27,
                  28, 30, 20, 21, 22, 33, 34, 24, 25, 26, 27, 28, 30, 20, 21,
                  22, 33, 34, 24, 25, 26, 27, 28, 30, 20, 21, 22, 33, 34};

// Default Canvas
Double_t w = 1920;
Double_t h = 1080;

void SaveCanvas(TCanvas* c, char const* name, TString path, char const* type);
void SavePad(TPad* p, char const* name, TString path, char const* type);
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
TCanvas* MakeRatioCanvas(char const* name, TH1D* hTop, TH1D* hBottom);
double myLevyPt(Double_t* x, Double_t* par);
Double_t linearSqrtBkg(Double_t* x, Double_t* par);
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
void DrawXi1530forXi1820(const int sys = 1,
                double multi_start = 0,
                double multi_end = 100,
                char const* inputOptions = "",
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
    bool isInelastic = false;
    // Reaction to the option
    TString Options = inputOptions;

    // For fit systematics
    // NOT USING NOW -> Replaced to FitVar, BinCount, NormVar
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
        inputDirectory = "Xi1530HM";
        fHM = true;
        fullcentbin = {0, 0.1};
    }
    // Fit variation check
    if (Options.Contains("FitVarLm")){
        FitRange[0] = FitRange[0]-0.001*rebin*OptionNumber;
    }
    if (Options.Contains("FitVarLp")){
        FitRange[0] = FitRange[0]+0.001*rebin*OptionNumber;
    }
    if (Options.Contains("FitVarRm")){
        FitRange[1] = FitRange[1]-0.001*rebin*OptionNumber;
    }
    if (Options.Contains("FitVarRp")){
        FitRange[1] = FitRange[1]+0.001*rebin*OptionNumber;
    }
    if (Options.Contains("FitVarBothm")){
        FitRange[0] = FitRange[0]-0.001*rebin*OptionNumber;
        FitRange[1] = FitRange[1]+0.001*rebin*OptionNumber;
    }
    if (Options.Contains("FitVarBothp")){
        FitRange[0] = FitRange[0]+0.001*rebin*OptionNumber;
        FitRange[1] = FitRange[1]-0.001*rebin*OptionNumber;
    }
    // Norm change
    if (Options.Contains("NormVarm")){
        NormalizeRange_L[0] = NormalizeRange_L[0]-0.001*rebin*OptionNumber;
        NormalizeRange_L[1] = NormalizeRange_L[1]-0.001*rebin*OptionNumber;
        NormalizeRange_R[0] = NormalizeRange_R[0]-0.001*rebin*OptionNumber;
        NormalizeRange_R[1] = NormalizeRange_R[1]-0.001*rebin*OptionNumber;
    }
    if (Options.Contains("NormVarp")){
        NormalizeRange_L[0] = NormalizeRange_L[0]+0.001*rebin*OptionNumber;
        NormalizeRange_L[1] = NormalizeRange_L[1]+0.001*rebin*OptionNumber;
        NormalizeRange_R[0] = NormalizeRange_R[0]+0.001*rebin*OptionNumber;
        NormalizeRange_R[1] = NormalizeRange_R[1]+0.001*rebin*OptionNumber;
    }
    if (Options.Contains("NormVarLp")){
        NormalizeRange_L[1] = NormalizeRange_L[1]+0.001*rebin*OptionNumber;
    }
    if (Options.Contains("NormVarLm")){
        NormalizeRange_L[0] = NormalizeRange_L[0]-0.001*rebin*OptionNumber;
    }
    if (Options.Contains("NormVarRp")){
        normleft = 0;
        normright = 1;
        NormalizeRange_R[1] = NormalizeRange_R[1]+0.001*rebin*OptionNumber;
    }
    if (Options.Contains("NormVarRm")){
        normleft = 0;
        normright = 1;
        NormalizeRange_R[0] = NormalizeRange_R[0]-0.001*rebin*OptionNumber;
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
        //formula = "[0] + [1]*x + [2]*TMath::Sqrt( x - 3.0638 )"; //3.0636 = 2 * 1.5318
        formula = "pol2(0)";
        //FitRange={1.496,1.572};
    }
    if (Options.Contains("MCcheck")){
        fullcentbin = centbin;
    }
    if (Options.Contains("inel")){
        isInelastic = true;
    }
    if (Options.Contains("Corey")){
        formula = "pol0(0)";
        FitRange[0] = 1.5;
        FitRange[1] = 1.56;
    }

    // Preparation
    for (int i = 0; i < (int)ptbin.size() - 1; i++)
        pt_points_e.push_back(ptbin.at(i + 1) - ptbin.at(i));

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

    //Temporal Check!!
    Double1D integralfullfunction = {};
    Double1D integralsigfunction = {};
    Double1D integralbkgfunction = {};
    //

    double triggereffi = inf;
    double triggereffi_e = zero;

    // Load DATA
    auto clist = LoadXi1530ResultList(datafile, inputDirectory);
    auto hInvMass = BSTHnSparseHelper::Load("hInvMass", clist);
    auto clist_MC = LoadXi1530ResultList(
        rsnmcfile, "Xi1530MB");  // From Resonance Injected MC
    auto hInvMass_MC = BSTHnSparseHelper::Load("hInvMass", clist_MC);
    auto hInvMass_MC_MB = BSTHnSparseHelper::Load("hInvMass", clist_MC);
    auto clist_MC_General = LoadXi1530ResultList(
         genmcfile, inputDirectory);  // From General Purpose MC
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
    hInvMass_MC.SetBin("Cent", centbin);          // for Trigger Efficiency
    hInvMass_MC_MB.SetBin("Cent", fullcentbin_Effi);   // for MC reconstruction Efficiency, we use full bin
    hInvMass_MC_General.SetBin("Cent", centbin);  // for Trigger Efficiency
    hInvMass_MC_General_MB.SetBin("Cent", fullcentbin_Effi);  // for MC reconstruction Efficiency, we use full bin

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
        "hEventNumbers");  // N of Event through event cuts
    TH1D* hMultQA = (TH1D*)clist->FindObject(
        "hMult_QA");  // Multiplicty distribution after all event cuts
    TH1D* hMultQA_onlyMult = (TH1D*)clist->FindObject(
        "hMult_QA_onlyMult");  // Multiplicty distribution after only
    // MultiSelection + Trigger selection

    hNumberofEvent->GetXaxis()->SetTitle("Event cuts");
    hNumberofEvent->GetYaxis()->SetTitle("# of Event");
    hNumberofEvent->Write("hNumberofEvent");

    double eventfraction = 0.;
    if (!fHM)
        eventfraction =
            hMultQA->Integral(hMultQA->GetXaxis()->FindBin(multi_start),
                              hMultQA->GetXaxis()->FindBin(multi_end)) /
            hMultQA->Integral(hMultQA->GetXaxis()->FindBin(0.),
                              hMultQA->GetXaxis()->FindBin(100));
    else
        eventfraction =
            hMultQA->Integral(hMultQA->GetXaxis()->FindBin(multi_start),
                              hMultQA->GetXaxis()->FindBin(multi_end)) /
            hMultQA->Integral(hMultQA->GetXaxis()->FindBin(0.),
                              hMultQA->GetXaxis()->FindBin(0.1));

    hMultQA->Rebin(10);
    hMultQA->GetXaxis()->SetTitle("Multiplicity Percentile (%)");
    hMultQA->GetYaxis()->SetTitle("# of Event");
    hMultQA->Write("hMultQA");
    // Trigger Efficiency
    // --------------------------------------------------------------------------------------
    TH1D* hNumberofEvent_Genenral = (TH1D*)clist_MC_General->FindObject(
        "hEventNumbers");  // N of Event through event cuts
    double nEventin_elastic = hNumberofEvent_Genenral->GetBinContent(1); // I think.. it is wrong.
    auto hInvMass_MC_General_Trig_TrueINELg0 =
        BSTHnSparseHelper::Load("htriggered_CINT7", clist_MC_General);
    auto hInvMass_MC_General_Trig_INT7 =
        BSTHnSparseHelper::Load("htriggered_CINT7", clist_MC_General);
    auto htrue_cent =
        hInvMass_MC_General_Trig_TrueINELg0.GetTH1("true", 1, {1, -1, -1});
    auto hReco_cent =
        hInvMass_MC_General_Trig_INT7.GetTH1("Reco", 1, {4, -1, -1});
    double sumtrue = 0;
    double sumreco = 0;
    for (int k = 1; k < (int)htrue_cent->GetNbinsX()+1; k++) {
        double check = htrue_cent->GetBinCenter(k);
        if (((double)centbin[0] < check) && (check < (double)centbin[1])) {
            sumtrue += htrue_cent->GetBinContent(k);
            sumreco += hReco_cent->GetBinContent(k);
        }
    }
    triggereffi = sumreco / sumtrue;
    //triggereffi = sumreco / nEventin_elastic; // for test!
    triggereffi_e = sqrt(pow(sqrt(sumreco)/sumtrue,2)+pow(sqrt(sumtrue)*sumreco/pow(sumtrue,2),2));
    Double1D tempbin = {multi_start, multi_end};
    if(fHM){
        triggereffi = 1;
        triggereffi_e = zero;
    }
    if(isInelastic) {
        // from Normalisation factor study: https://alice-notes.web.cern.ch/node/665
        triggereffi = 0.7448;
        triggereffi_e = 0.0190;
    }
    TH1D* htriggereffi =
        new TH1D("TrigEffi", "Trigger Efficiency", 1, &tempbin[0]);
    htriggereffi->SetBinContent(1, triggereffi);
    htriggereffi->SetBinError(1, triggereffi_e);
    htriggereffi->Write("hTriggerEffi");

    // Signal with Bkg
    // -----------------------------------------------------------------------------------------
    for (auto j : hInvMass.BinRange("Pt")) {
        int n = j - 1;
        // SigBkg
        // ---------------------------------------------------------------------------------------------
        normleft = 1;
        normright = 1;
        auto hSig =
            hInvMass.GetTH1("Signal", 4, {sys, 1, 1, j, -1});  // -1 => inv mass
        auto hBkg = hInvMass.GetTH1("Bkg", 4,
                                    {1, bkgtype, 1, j, -1});  // -1 => inv mass
        auto hBkg_norm = GetNorBkg(hSig, hBkg, normleft, normright);

        // rebin
        hSig->Rebin(rebin);
        hBkg_norm->Rebin(rebin);
        
        // Axis
        hSig->GetXaxis()->SetTitle("Mass(#pi#Xi) (GeV/c^{2})");
        hSig->GetYaxis()->SetTitle(
            Form("# of candidates / %.1f MeV", 1000 * hSig->GetBinWidth(1)));
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
            if (check <= NormalizeRange_L.at(0))
                temp->SetBinContent(k, 0);
            if (!normleft && (NormalizeRange_L.at(0) <= check &&
                              check <= NormalizeRange_L.at(1)))
                temp->SetBinContent(k, 0);
            if (NormalizeRange_L.at(1) <= check &&
                check <= NormalizeRange_R.at(0))
                temp->SetBinContent(k, 0);
            if (!normright && (NormalizeRange_R.at(0) <= check &&
                               check <= NormalizeRange_R.at(1)))
                temp->SetBinContent(k, 0);
            if (NormalizeRange_R.at(1) <= check)
                temp->SetBinContent(k, 0);
        }
        temp->SetFillColorAlpha(kRed, 0.05);

        // Draw
        hSig->Write(Form("hSignalOnly_%i", n));
        hBkg_norm->Write(Form("hBkgOnly_%i", n));
        temp->Write(Form("hBkgNorm_%i", n));

        // Fit
        // ------------------------------------------------------------------------------------------------
        TH1D* Sigfit = (TH1D*)hSig->Clone();
        TH1D* Bkgfit = (TH1D*)hBkg_norm->Clone();
        if(!fbkgfit){
            Sigfit->Add(Bkgfit, -1);
        }
        
        TGaxis::SetMaxDigits(3);
        Sigfit->Write(Form("hSignalBkgSubtraction_%i", n));
        TF1* fitresult = VoigtFit(Sigfit, formula, peakRange, par, fix);
        fitresult->Write(Form("fDataFitResult_%i", n));
        TF1* fonly = new TF1("Voigt", "[0]*TMath::Voigt(x-[1],[2],[3])",
                             DrawRange[0], DrawRange[1]);
        fonly->SetParameter(0, fitresult->GetParameter(vp));
        fonly->SetParError(0, fitresult->GetParError(vp));
        fonly->SetParameter(1, fitresult->GetParameter(vp + 1));
        fonly->SetParError(1, fitresult->GetParError(vp + 1));
        fonly->SetParameter(2, fitresult->GetParameter(vp + 2));
        fonly->SetParError(2, fitresult->GetParError(vp + 2));
        fonly->SetParameter(3, fitresult->GetParameter(vp + 3));
        fonly->SetParError(3, fitresult->GetParError(vp + 3));
        fonly->Write(Form("fDataFitOnlyResult_%i", n));
        fb_after->Write(Form("fDataFitBkgResult_%i", n));

        fitmean.push_back(fitresult->GetParameter(vp + 1));
        fitmean_err.push_back(fitresult->GetParError(vp + 1));
        fitsigma.push_back(fitresult->GetParameter(vp + 2));
        fitsigma_err.push_back(fitresult->GetParError(vp + 2));
        fitgamma.push_back(fitresult->GetParameter(vp + 3));
        fitgamma_err.push_back(fitresult->GetParError(vp + 3));
        integralfullfunction.push_back(1000 * fitresult->Integral(IntegralRange[0], IntegralRange[1]) /
                rebin);
        integralsigfunction.push_back(1000 * fonly->Integral(IntegralRange[0], IntegralRange[1]) /
                rebin);
        integralbkgfunction.push_back(1000 * fb_after->Integral(IntegralRange[0], IntegralRange[1]) /
                rebin);
        
        double frawyield = 0.;
        double frawyield_e = 0.;
        if (bincount) {  // Bin count method
            frawyield =
                Sigfit->Integral(Sigfit->GetXaxis()->FindBin(IntegralRange[0]),
                                 Sigfit->GetXaxis()->FindBin(IntegralRange[1] -
                                                             0.001 * rebin));
            
            for (int bin = Sigfit->GetXaxis()->FindBin(IntegralRange[0]);
                 bin <=
                 Sigfit->GetXaxis()->FindBin(IntegralRange[1]+0.001*rebin);
                 bin++) {
                frawyield_e += pow(Sigfit->GetBinError(bin), 2);
            }
            frawyield += -bkgintegral;

            // lost yield correction (additive)
            double lefttail =
                1000 * fonly->Integral(IntegralRangeMC[0], IntegralRange[0]) /
                rebin;
            double righttail =
                1000 * fonly->Integral(IntegralRange[1], IntegralRangeMC[1]) /
                rebin;
            frawyield = frawyield + lefttail + righttail;
            frawyield_e = frawyield_e 
                + pow(lefttail*fonly->IntegralError(
                    IntegralRangeMC[0], IntegralRange[0]),2)
                + pow(righttail*fonly->IntegralError(
                    IntegralRange[1], IntegralRangeMC[1]),2);
            RawYield_err.push_back(sqrt(frawyield_e));
        } else {  // Integral from fit function
            frawyield = 1000 *
                        fonly->Integral(IntegralRangeMC[0], IntegralRangeMC[1]) /
                        rebin;
            /*
            frawyield_e =
                1000 * fonly->Integral(IntegralRangeMC[0], IntegralRangeMC[1]) *
                fonly->IntegralError(IntegralRangeMC[0], IntegralRangeMC[1]) /
                fonly->Integral(IntegralRangeMC[0], IntegralRangeMC[1]);
            */
            frawyield_e = frawyield*(fonly->GetParError(0)/fonly->GetParameter(0));
            RawYield_err.push_back(frawyield_e);
        }
        RawYield.push_back(frawyield);

        // MC data
        // --------------------------------------------------------------------------------------------
        // After Trigger selection Gen MC
        auto hTrueInput_Gen =
            hInvMass_MC_General_MB.GetTH1("trueinput", 4, {1, 7, 1, j, -1});
        // After All event cut Gen MC
        auto hInput_Gen =
            hInvMass_MC_General_MB.GetTH1("input", 4, {1, 5, 1, j, -1});
        // True After All event cut
        auto hInput = hInvMass_MC_MB.GetTH1("input", 4, {1, 6, 1, j, -1});
        // My reconstruction in All event cut
        auto hReco = hInvMass_MC_MB.GetTH1("recon", 4, {sys, 4, 1, j, -1});

        Double_t Input_number_Gen = 0.;
        Double_t True_Input_number_Gen = 0.;

        Double_t Input_number = 0.;
        Double_t Reco_number = 0.;
        for (int k = 0; k < (int)hInput->GetNbinsX(); k++) {
            auto check = hInput->GetBinCenter(k);
            if ((check >= IntegralRangeMC[0]) &&
                (check <= IntegralRangeMC[1])) {
                True_Input_number_Gen += hTrueInput_Gen->GetBinContent(k);
                Input_number_Gen += hInput_Gen->GetBinContent(k);

                Input_number += hInput->GetBinContent(k);
                Reco_number += hReco->GetBinContent(k);
            }
        }
        Double_t Eff_e = sqrt(
            pow(sqrt(Reco_number) / Input_number, 2) +
            pow(sqrt(Input_number) * Reco_number / pow(Input_number, 2), 2));
        Double_t ecutratio = True_Input_number_Gen / Input_number_Gen;
        Double_t ecutratio_e = sqrt(pow(sqrt(True_Input_number_Gen)/Input_number_Gen,2)+pow(sqrt(Input_number_Gen)*True_Input_number_Gen/pow(Input_number_Gen,2),2));
        EventCutRatio.push_back(ecutratio);
        EventCutRatio_e.push_back(ecutratio_e);
        Efficiency.push_back(Reco_number / Input_number);
        Efficiency_e.push_back(Eff_e);

        // MC Fit
        hReco->Rebin(rebin);
        hReco->GetXaxis()->SetRangeUser(DrawRange[0], DrawRange[1]);
        hReco->Write(Form("hMCRecon_%i", n));

        TF1* MCfitresult = VoigtFit(hReco, formula_MC, peakRange, par, MCfix);

        fitmean_MC.push_back(MCfitresult->GetParameter(vp + 1));
        fitmean_MC_err.push_back(MCfitresult->GetParError(vp + 1));
        fitsigma_MC.push_back(MCfitresult->GetParameter(vp + 2));
        fitsigma_MC_err.push_back(MCfitresult->GetParError(vp + 2));
        fitgamma_MC.push_back(MCfitresult->GetParameter(vp + 3));
        fitgamma_MC_err.push_back(MCfitresult->GetParError(vp + 3));

        MCfitresult->Write(Form("fMCFit_%i", n));
    }
    //-----------------------------------------------------------------------------------------------------
    // MC Efficiency
    auto hEfficiency =
        MakeHistfromArray("MC Efficiency", Efficiency, Efficiency_e, ptbin);
    hEfficiency->SetBinContent(1, 1e-10);
    hEfficiency->SetBinError(1, 1e-10);
    hEfficiency->SetMinimum(4e-4);
    hEfficiency->SetMaximum(4e-1);
    hEfficiency->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hEfficiency->GetYaxis()->SetTitle("Acceptance x Efficiency x BR");
    hEfficiency->Write("hMCReconEffi");

    // ----------------------------------------------------------------------------------------------------
    // Signal loss
    auto hsigloss =
        MakeHistfromArray("Signal loss", EventCutRatio, EventCutRatio_e, ptbin);
    hsigloss->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hsigloss->GetYaxis()->SetTitle("Signal Loss");
    hsigloss->SetMinimum(0);
    hsigloss->Write("hMCSigLoss");

    // Temporal check
    auto hFullIntegral = MakeHistfromArray("hFullIntegral", integralfullfunction, ptzero, ptbin);
    auto hSigIntegral = MakeHistfromArray("hSigIntegral", integralsigfunction, ptzero, ptbin);
    auto hBkgIntegral = MakeHistfromArray("hBkgIntegral", integralbkgfunction, ptzero, ptbin);
    hFullIntegral->Write("hFullIntegral");
    hSigIntegral->Write("hSigIntegral");
    hBkgIntegral->Write("hBkgIntegral");
    
    for (auto j : hInvMass.BinRange("Pt")) {
        int n = j - 1;
        cout << "Full/Sig/Bkg: " << integralfullfunction[n] << "/" << integralsigfunction[n] << "/" << integralbkgfunction[n] << endl;
    }
    

    // Spectra
    // --------------------------------------------------------------------------------------------
    auto hRawYields =
        MakeHistfromArray("Raw Yields", RawYield, RawYield_err, ptbin);
    hRawYields->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hRawYields->GetYaxis()->SetTitle("Raw yield");
    hRawYields->Write("hDataRawYield");

    auto hXispectrum = MakeHistfromArray("hXispectrum", ptzero, ptzero, ptbin);
    hXispectrum->SetMarkerStyle(20);
    hXispectrum->SetMarkerSize(1.0);
    hXispectrum->SetTitle("#Xi(1530) spectrum pp #sqrt{s}= 13 TeV");
    hXispectrum->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hXispectrum->GetYaxis()->SetTitle(
        "1/N_{event}d^{2}N/(dydp_{T}) (GeV/c)^{-1}");

    Double_t value = 0.;
    Double_t base = 0.;
    Double_t err = 0.;
    Double_t f_vtx = 0.931264; 
    //Double_t f_vtx = hNumberofEvent->GetBinContent(6)/hNumberofEvent->GetBinContent(5); 
    //hXispectrum->SetBinContent(1, zero);
    //hXispectrum->SetBinError(1, zero);
    hXispectrum->SetMinimum(1.0e-9);
    hXispectrum->SetMaximum(3.0e-2);
    double nOfEventMultibin = hNumberofEvent->GetBinContent(9) * eventfraction;
    for (int ii = 1; ii < (int)ptbin.size() - 2; ii++) {
        cout << "RawYield.at(ii): " << RawYield.at(ii)
             << ", Efficiency.at(ii): " << Efficiency.at(ii)
             << ", EventCutRatio.at(ii): " << EventCutRatio.at(ii)
             << ", centbin.at(0): " << centbin.at(0)
             << ", centbin.at(1): " << centbin.at(1)
             << ", # of Events: " << hNumberofEvent->GetBinContent(9)
             << ", # of Events in mutlbin: " << nOfEventMultibin
             << ", pt_points_e.at(ii): " << pt_points_e.at(ii) 
             << ", Vertex lost correction: " << f_vtx << endl;
        base = Efficiency.at(ii) * (pt_points_e.at(ii)) *
               nOfEventMultibin / (triggereffi * EventCutRatio.at(ii) * f_vtx);
        err = pow(RawYield_err.at(ii) / base, 2);
        err += pow(
            (RawYield.at(ii) / base) * (Efficiency_e.at(ii) / Efficiency.at(ii)),
            2);
        err += pow(
            (RawYield.at(ii) / base) * (EventCutRatio_e.at(ii) / EventCutRatio.at(ii)),
            2);
        err += pow(
            (RawYield.at(ii) / base) * (triggereffi_e/triggereffi),
            2);
        err += pow(
            (RawYield.at(ii) / base) * (sqrt(nOfEventMultibin) / nOfEventMultibin),
            2);
        err = sqrt(err);
        value = RawYield.at(ii) / base;
        value /= 2.;  // Xi(1530)0 + c.c.

        cout << "Corrected yield: " << value << endl;
        cout << "Error: " << err/2 << endl;
        hXispectrum->SetBinContent(ii + 1, value);
        hXispectrum->SetBinError(ii + 1, err / 2);
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

void SaveCanvas(TCanvas* c,
                char const* name = "temp",
                TString path = "figs/",
                char const* type = "pdf") {
    // Save canvas with path. if path is not there, make a folder.
    if (gSystem->Exec(Form("ls %s", path.Data())) == 256)
        gSystem->Exec(Form("mkdir -p %s", path.Data()));
    c->SaveAs(Form("%s%s.%s", path.Data(), name, type));
}
void SavePad(TPad* p,
             char const* name = "temp",
             TString path = "figs/",
             char const* type = "pdf") {
    // Save Pad with path.
    TCanvas* ctemp = new TCanvas();
    TPad* clone = (TPad*)p->DrawClone();
    clone->SetPad(0, 0, 1, 1);
    SaveCanvas(ctemp, name, path, type);
}

TH1D* GetNorBkg(TH1D* hSig, TH1D* hBkg, int isnormleft, int isnormright) {
    Double_t normalization_data = 0.;
    Double_t normalization_mixed = 0.;

    if (isnormleft == 1) {
        normalization_data +=
            hSig->Integral(hSig->GetXaxis()->FindBin(NormalizeRange_L.at(0)),
                           hSig->GetXaxis()->FindBin(NormalizeRange_L.at(1)));
        normalization_mixed +=
            hBkg->Integral(hBkg->GetXaxis()->FindBin(NormalizeRange_L.at(0)),
                           hBkg->GetXaxis()->FindBin(NormalizeRange_L.at(1)));
    }
    if (isnormright == 1) {
        normalization_data +=
            hSig->Integral(hSig->GetXaxis()->FindBin(NormalizeRange_R.at(0)),
                           hSig->GetXaxis()->FindBin(NormalizeRange_R.at(1)));
        normalization_mixed +=
            hBkg->Integral(hBkg->GetXaxis()->FindBin(NormalizeRange_R.at(0)),
                           hBkg->GetXaxis()->FindBin(NormalizeRange_R.at(1)));
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
              double* fixinput) {
    int j;
    TH1D* a = (TH1D*)hSig->Clone(Form("%s_nopeak", name));
    for (j = hSig->GetXaxis()->FindBin(1.000001 * peakRange[0]);
         j <= hSig->GetXaxis()->FindBin(0.999999 * peakRange[1]); j++) {
        a->SetBinContent(j, 0.);
        a->SetBinError(j, 0.);
    }
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); // better minimizer!
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

    if (!fixinput[2]) {  // release width
        fp->ReleaseParameter(vp + 3);
        fp->SetParError(vp + 3, 0.1 * fp->GetParameter(vp + 3));
        hSig->Fit(fp, "RQN");
    }

    if (!fixinput[1]) {  // release resolution
        fp->ReleaseParameter(vp + 2);
        fp->SetParLimits(vp+2,0.,0.01);
        fp->SetParError(vp + 2, 0.1 * fp->GetParameter(vp + 2));
        hSig->Fit(fp, "RQN");
    }

    if (!fixinput[0]) {  // release mass
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
        1000 * fb_after->Integral(IntegralRange[0], IntegralRange[1]) / rebin;

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
        htemp->SetBinContent(i + 1, dArray.at(i));
        if (!option.Contains("NOERROR"))
            htemp->SetBinError(i + 1, eArray.at(i));
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

Double_t linearSqrtBkg(Double_t* x, Double_t* par) {
    double lMass = 1.5318;  // Xi mass

    return par[0] + x[0]*par[1] + par[2]*TMath::Sqrt(x[0] - 2 * lMass );
}
TFile* LoadXi1530Results(TString name, TString runnum) {
    if (!runnum.IsNull())
        runnum = "_" + runnum;
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