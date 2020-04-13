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
// char const* datafile = "Xi1530LHC16klop_17kl_SYS_Mix";
TString inputDirectory = "Xi1530MB";

// FIT SETUP -----------------
double peakRange[2] = {1.51, 1.55};

char name[500] = "Voigt fit";
char const* formula = "pol2(0)";
char const* formula_MC = "pol0(0)";
double par[3] = {1.532, 0.003, 0.0091};
double par_2nd[3] = {1.532, 0.0025, 0.005};
double fix[3] = {0, 0, 1};
double fix_2nd[3] = {1, 0, 1};
double fix_3rd[3] = {1, 1, 0};
double fix_4th[3] = {0, 1, 1};
int NofParmBkg = 0;
TF1* fb;
Double_t bkgintegral = 0;
bool bincount = false;
bool fHM = false;
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
Int_t normleft = 1;
Double1D NormalizeRange_L = {1.49, 1.51};
Int_t normright = 0;
Double1D NormalizeRange_R = {1.56, 1.58};

// Axis Range
Double1D DrawRange = {1.47, 1.59};
Double1D FitRange = {1.47, 1.59};
Double1D IntegralRange = {1.52, 1.55};
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
              char const* bkgformula,
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

// File handling
TFile* LoadXi1530Results(TString name, TString runnum);
TObject* LoadXi1530ResultList(TFile* fh, TString clistname);
TObject* LoadXi1530ResultList(TString fname,
                              TString clistname,
                              TString runnum = "");
//__________________________________________________________
//__________________________________________________________
//_______________________MCQAXi1530_________________________
//__________________________________________________________
//__________________________________________________________
void MCQAXi1530(TString rsnmcfile = "LHC18c6d_vertexer") {
    // Start!
    const int sys = 1;
    cout << "== Xi(1530)0 MC QA Macro ==" << endl;
    cout << "Cut Systematic option(Default: 1): " << sys << endl;

    Double1D centbin = {0, 0.1};
    Double1D fullcentbin = {0, 0.1};

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
    Double1D Efficiency = {};
    Double1D Efficiency_e = {};

    // Load DATA
    auto clist_MC = LoadXi1530ResultList(
        rsnmcfile.Data(), inputDirectory);  // From Resonance Injected MC
    auto hInvMass_MC = BSTHnSparseHelper::Load("hInvMass", clist_MC);
    auto hInvMass_MC_MB = BSTHnSparseHelper::Load("hInvMass", clist_MC);

    /*
    auto clist_MC_INEL = LoadXi1530ResultList(
        rsnmcfile.Data(), "Xi1530INEL");  // From Resonance Injected MC
    auto hInvMass_MC_INEL = BSTHnSparseHelper::Load("hInvMass", clist_MC);
    */

    // pT binning
    hInvMass_MC_MB.SetBin("Pt", ptbin);
    //hInvMass_MC_INEL.SetBin("Pt", ptbin);

    // Multiplicity percentile binning
    hInvMass_MC.SetBin("Cent", centbin);          // for Trigger Efficiency
    hInvMass_MC_MB.SetBin("Cent", fullcentbin);   // for MC reconstruction Efficiency, we use full bin

    cout << "RSN MC AXES " << endl;
    hInvMass_MC.PrintAxis("all");

    // output files
    TString savefile =
        Form("AnalysisResults_Xi1530MCQA_%s.root", rsnmcfile.Data());
    TFile* output = new TFile(savefile.Data(), "RECREATE");
    
    vector<vector<double>> multibincheck = {
        {0, 100}, {0, 10}, {10, 30}, {30, 50}, {50, 70}, {70, 100}};
    vector<double> trigeffi;
    vector<double> trigeffie;
    vector<double> verteffi;
    vector<double> verteffie;
    //
    // Trigger Efficiency
    // --------------------------------------------------------------------------------------
    TH1D* hNumberofEvent = (TH1D*)clist_MC->FindObject(
        "fNormalisationHist");  // N of Event through event cuts
    auto hInvMass_MC_General_Trig =
        BSTHnSparseHelper::Load("htriggered_CINT7", clist_MC);
    
    /*
    auto hInvMass_MC_General_Trig_INEL =
        BSTHnSparseHelper::Load("htriggered_CINT7", clist_MC_INEL);
    */

    auto htrue_cent =
        hInvMass_MC_General_Trig.GetTH1("true", 1, {1, -1, -1}); // MC True INEL>0
    auto hReco_cent =
        hInvMass_MC_General_Trig.GetTH1("Reco", 1, {2, -1, -1}); // True INEL>0 && Triggered
    auto hVtx_cent =
        hInvMass_MC_General_Trig.GetTH1("Vtx", 1, {3, -1, -1}); // True INEL>0 && triggered && good vtx  

    auto htrue_cent_INEL =  
        hInvMass_MC_General_Trig.GetTH1("true", 1, {6, -1, -1});
    auto hReco_cent_INEL =
        hInvMass_MC_General_Trig.GetTH1("Reco", 1, {7, -1, -1}); // True INEL>0 && Triggered
    auto hVtx_cent_INEL =
        hInvMass_MC_General_Trig.GetTH1("Vtx", 1, {8, -1, -1}); // True INEL>0 && triggered && good vtx  
    for (int imultibin = 0; imultibin < multibincheck.size(); imultibin++) {    
        double sumtrue = htrue_cent->Integral(htrue_cent->GetXaxis()->FindBin(multibincheck[imultibin][0]+0.0001),
                htrue_cent->GetXaxis()->FindBin(multibincheck[imultibin][1]-0.0001));
        double sumreco = hReco_cent->Integral(hReco_cent->GetXaxis()->FindBin(multibincheck[imultibin][0]+0.0001),
                hReco_cent->GetXaxis()->FindBin(multibincheck[imultibin][1]-0.0001));
        double sumvtx = hVtx_cent->Integral(hVtx_cent->GetXaxis()->FindBin(multibincheck[imultibin][0]+0.0001),
                hVtx_cent->GetXaxis()->FindBin(multibincheck[imultibin][1]-0.0001));

        double triggereffi = inf;
        double triggereffi_e = zero;

        double vertexeffi = inf;
        double vertexeffi_e = zero;

        triggereffi = sumreco / sumtrue;
        triggereffi_e = sqrt(triggereffi*(1-triggereffi)/sumtrue);
        //triggereffi_e = sqrt(pow(sqrt(sumreco)/sumtrue,2)+pow(sqrt(sumtrue)*sumreco/pow(sumtrue,2),2));

        // Vertexer Loss Correction
        // -----------------------------------------------------------------------------------------
        vertexeffi = sumvtx/sumreco;
        vertexeffi_e = sqrt(vertexeffi*(1-vertexeffi)/sumreco);
        //vertexeffi_e = sqrt(pow(sqrt(sumvtx)/sumreco,2)+pow(sqrt(sumreco)*sumvtx/pow(sumreco,2),2));
        
        cout << "trigeffi : " << triggereffi << " +- " << triggereffi_e << ", vertexeffi: " << vertexeffi << " +- " << vertexeffi_e << endl;


        if(multibincheck[imultibin][1] - multibincheck[imultibin][0] > 99){
            //MB
            vector<double> MBbin = {0,100};

            TH1D* htriggereffiMB =
                new TH1D("TrigEffiMB","Trigger Efficiency MB", 1, &MBbin[0]);
            htriggereffiMB->SetBinContent(1, triggereffi);
            htriggereffiMB->SetBinError(1, triggereffi_e);
            htriggereffiMB->Write("hTrigEffiMB");

            TH1D* hvertexeffi =
                new TH1D("hVertexEffiMB","Vertex Efficiency MB", 1, &MBbin[0]);
            hvertexeffi->SetBinContent(1, vertexeffi);
            hvertexeffi->SetBinError(1, vertexeffi_e);
            hvertexeffi->Write("hVertexEffiMB");
        }
        else{
            trigeffi.push_back(triggereffi);
            trigeffie.push_back(triggereffi_e);
            verteffi.push_back(vertexeffi);
            verteffie.push_back(vertexeffi_e);
        }
    }
    vector<double> multibin = {0,10,30,50,70,100};

    auto hTrigMult =
        MakeHistfromArray("Trigger Efficiency Multibin", trigeffi, trigeffie, multibin);
    hTrigMult->Write("hTrigEffi");
    auto hVertMult =
        MakeHistfromArray("Vertex Efficiency Multibin", verteffi, verteffie, multibin);
    hVertMult->Write("hVertexEffi");

    double sumtrue_INEL = htrue_cent_INEL->GetEntries();
    double sumreco_INEL = hReco_cent_INEL->GetEntries();
    double sumvtx_INEL = hVtx_cent_INEL->GetEntries();

    double triggereffi = sumreco_INEL / sumtrue_INEL;
    double triggereffi_e = sqrt(triggereffi*(1-triggereffi)/sumtrue_INEL);
    double vertexeffi = sumvtx_INEL/sumreco_INEL;
    double vertexeffi_e = sqrt(vertexeffi*(1-vertexeffi)/sumreco_INEL);

    vector<double> INELbin = {0,100};
    TH1D* htriggereffiINEL =
        new TH1D("TrigEffiINEL","Trigger Efficiency INEL", 1, &INELbin[0]);
    htriggereffiINEL->SetBinContent(1, triggereffi);
    htriggereffiINEL->SetBinError(1, triggereffi_e);
    htriggereffiINEL->Write("hTrigEffiINEL");

    TH1D* hvertexeffiINEL =
        new TH1D("hVertexEffiINEL","Vertex Efficiency MB", 1, &INELbin[0]);
    hvertexeffiINEL->SetBinContent(1, vertexeffi);
    hvertexeffiINEL->SetBinError(1, vertexeffi_e);
    hvertexeffiINEL->Write("hVertexEffiINEL");
    
    //gSystem->Exit(1);

    
    // Signal with Bkg
    // -----------------------------------------------------------------------------------------
    for (auto j : hInvMass_MC_MB.BinRange("Pt")) {
        int n = j - 1;
        // MC data
        // --------------------------------------------------------------------------------------------
        // After PS
        auto hInput = hInvMass_MC_MB.GetTH1("input", 4, {1, 5, 1, j, -1});
        // My reconstruction
        auto hReco = hInvMass_MC_MB.GetTH1("recon", 4, {sys, 4, 1, j, -1});

        double sumtrue = hInput->Integral(hInput->GetXaxis()->FindBin(IntegralRangeMC[0]+0.0001),
            hInput->GetXaxis()->FindBin(IntegralRangeMC[1]-0.0001));
        double sumreco = hReco->Integral(hReco->GetXaxis()->FindBin(IntegralRangeMC[0]+0.0001),
            hReco->GetXaxis()->FindBin(IntegralRangeMC[1]-0.0001));

        Double_t Input_number = 0.;
        Double_t Reco_number = 0.;
        for (int k = 0; k < (int)hInput->GetNbinsX(); k++) {
            auto check = hInput->GetBinCenter(k);
            if ((check >= IntegralRangeMC[0]) &&
                (check <= IntegralRangeMC[1])) {
                Input_number += hInput->GetBinContent(k);
                Reco_number += hReco->GetBinContent(k);
            }
        }
        cout << "Old way: " << Reco_number << " / " << Input_number << " -> " << Reco_number/Input_number <<  endl;
        cout << "New way: " << sumreco << " / " << sumtrue << " -> " << sumreco/sumtrue <<  endl;
        Double_t RecEff = Reco_number/Input_number;
        if(RecEff == 0) RecEff = zero;
        Double_t Eff_e = sqrt(RecEff*(1-RecEff)/Input_number);
        Efficiency.push_back(RecEff);
        Efficiency_e.push_back(Eff_e);
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

    auto hInput = hInvMass_MC.GetTH1("input", 3, {1, 6, 1, -1, -1});
    hInput->SetMinimum(1);
    hInput->SetMaximum(1e5);
    hInput->SetMarkerColor(2);
    hInput->SetLineColor(2);
    hInput->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    auto hReco = hInvMass_MC.GetTH1("recon", 3, {1, 4, 1, -1, -1});
    hReco->SetMinimum(1);
    hReco->SetMaximum(1e5);
    hReco->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");

    TCanvas* cSame = new TCanvas("cSame", "", w, h);
    cSame->Draw();
    cSame->SetTickx();
    cSame->SetLogy(true);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);

    hInput->Draw();
    hInput->Write("hInputMC");
    hReco->Draw("Same");
    hReco->Write("hRecoMC");

    cSame->Write("cMCRecoInput");
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
              char const* bkgformula,
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

    fb = new TF1(Form("%s_back", name), bkgformula, DrawRange[0], DrawRange[1]);
    a->Fit(fb, "RQN");

    fb->SetLineColor(3);
    fb->SetLineStyle(2);
    fb->Draw("same");

    int vp = fb->GetNpar();

    TF1* fp = new TF1("Voigt+bkg",
                      Form("%s+[%i]*TMath::Voigt(x-[%i],[%i],[%i])", bkgformula,
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
    fp->Draw("same");

    // signal only for integral
    TF1* fonly = new TF1("Voigt", "[0]*TMath::Voigt(x-[1],[2],[3])",
                         DrawRange[0], DrawRange[1]);
    fonly->SetParameter(0, fp->GetParameter(vp));
    fonly->SetParameter(1, fp->GetParameter(vp + 1));
    fonly->SetParameter(2, fp->GetParameter(vp + 2));
    fonly->SetParameter(3, fp->GetParameter(vp + 3));
    fonly->SetLineColor(4);
    fonly->SetLineStyle(2);
    // fonly->Draw("same");

    NofParmBkg = vp;

    bkgintegral =
        1000 * fb->Integral(IntegralRange[0], IntegralRange[1]) / rebin;

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