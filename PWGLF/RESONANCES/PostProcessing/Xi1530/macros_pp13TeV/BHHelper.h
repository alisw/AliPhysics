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

// BHHelper
//
// Helper Class for THnSparse analysis
// (Almost)Based on BSHelper(Beomsoo Kim)

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
void __DEBUG_BH(int i, const char* fmt, ...) {
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
class BHRanger {
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
    BHRanger(int begin, int end) { SetRange(begin, end); }
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
BHRanger range(T& t, int begin = 0, int end = -1) {
    return BHRanger(begin, end < 0 ? t.size() - 1 : end);
}
BHRanger range(int end) {
    return BHRanger(0, end - 1);
}
BHRanger range(int begin, int end) {
    return BHRanger(begin, end);
}

BHRanger bin_range(int begin, int end) {
    if (begin < 1)
        ErrorExit("begin for bin_range must larger than 0 ");
    return range(begin, end);
}
BHRanger bin_range(int end) {
    return bin_range(1, end);
}
template <class T>
BHRanger bin_range(const std::vector<T>& t, int begin = 1, int end = min_int) {
    return bin_range(begin,
                     (end < 0 ? std::min(int(t.size() - 1), abs(end)) : end));
}
//---------------------------------------

//====================================
// HistManager
//====================================
class BHTHnSparseHelper {
   public:
    BHTHnSparseHelper();
    BHTHnSparseHelper(THnSparse* h) : fH(h) { LoadBinFromHist(); }
    BHTHnSparseHelper(TObject* o);
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

    BHRanger BinRange(int i) { return bin_range(GetBin(i)); }
    BHRanger BinRange(TString n) { return bin_range(GetBin(n)); }
    void ClearBin() { LoadBinFromHist(); }
    void ClearBin(int i) { LoadBinFromHist(i); }
    const Double1D& GetBin(int i) { return fCBins[i]; }
    const Double1D& GetBin(TString n) { return GetBin(GetAxisID(n)); }
    void SetBin(Double2D bins);
    void SetBin(int i, Double1D bins);
    void SetBin(TString name, Double1D bins) { SetBin(GetAxisID(name), bins); }
    static BHTHnSparseHelper Load(TString name, TObject* clist = nullptr);
    static BHTHnSparseHelper Load(TString name,
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
BHTHnSparseHelper::BHTHnSparseHelper() {
}

//__________________________________________________________
BHTHnSparseHelper::BHTHnSparseHelper(TObject* o) {
    if (!o)
        ErrorExit("No Object");
    fH = dynamic_cast<THnSparse*>(o);
    if (!fH)
        ErrorExit(
            Form("%s is not THnSparse but %s", o->GetName(), o->ClassName()));
    LoadBinFromHist();
}

//__________________________________________________________
TAxis* BHTHnSparseHelper::GetAxis(int i) {
    if (i < 0 || i >= GetNdim())
        ErrorExit(Form("Wrong Axis Index %d of %s", GetNdim(), fH->GetName()));
    return fH->GetAxis(i);
}
//__________________________________________________________
int BHTHnSparseHelper::GetAxisID(TString name) {
    for (int i = 0; i < fH->GetNdimensions(); i++)
        if (name == fH->GetAxis(i)->GetName())
            return i;
    fH->Print();
    PrintAxis();
    ErrorExit("No Axis " + name);
    return -1;
}

//__________________________________________________________
void BHTHnSparseHelper::PrintAxis(Option_t* opt) {
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
TH1D* BHTHnSparseHelper::GetTH1(TString name,
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
TH1D* BHTHnSparseHelper::GetTH1(TString name,
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
void BHTHnSparseHelper::SetBin(int iaxis, Double1D bins) {
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
void BHTHnSparseHelper::SetBin(Double2D bins) {
    for (UInt_t i = 0; i < bins.size(); i++)
        SetBin(i, bins[i]);
}

//__________________________________________________________
BHTHnSparseHelper BHTHnSparseHelper::Load(TString name, TObject* clist) {
    if (!clist)
        clist = gDirectory;
    auto o = clist->FindObject(name);
    if (!o)
        ErrorExit("No THnSparse " + name);
    return BHTHnSparseHelper(o);
}
//__________________________________________________________
BHTHnSparseHelper BHTHnSparseHelper::Load(TString name,
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
void BHTHnSparseHelper::LoadBinFromHist(int iaxis) {
    auto ax = GetAxis(iaxis);
    fCBinsB[iaxis] = {0};
    for (int ib = 1; ib <= ax->GetNbins() + 1; ib++) {
        fCBins[iaxis].push_back(ax->GetBinLowEdge(ib));
        fCBinsB[iaxis].push_back(ib);
    }
}
//__________________________________________________________
void BHTHnSparseHelper::LoadBinFromHist() {
    fCBins.resize(GetNdim());
    fCBinsB.resize(GetNdim());
    for (int i = 0; i < GetNdim(); i++)
        LoadBinFromHist(i);
}
//__________________________________________________________
// Initial values