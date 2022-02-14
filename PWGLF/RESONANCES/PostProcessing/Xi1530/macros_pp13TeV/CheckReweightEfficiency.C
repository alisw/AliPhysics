#include <AliPWGHistoTools.h>
#include "AdditionalFunctions.h"
#include "ReweightEfficiency.C"
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
TFile* LoadXi1530Results(TString name, TString runnum);
TObject* LoadXi1530ResultList(TFile* fh, TString clistname);
TObject* LoadXi1530ResultList(TString fname,
                              TString clistname,
                              TString runnum = "");
void CheckReweightEfficiency(){
    Double1D ptbin = {0.8, 1.2, 1.6, 2.0, 2.4, 3.2, 4.0, 4.8, 5.6, 8.8};
    Double1D centbin = {0.0, 100.0};

    auto clist1 = LoadXi1530ResultList("LHC18c6b4_298", "Xi1530MB");
    auto hInvMass1 = BSTHnSparseHelper::Load("hInvMass", clist1);
    hInvMass1.SetBin("Cent", centbin);
    hInvMass1.SetBin("Pt", ptbin);

    auto hInput = hInvMass1.GetTH1("input", 3, {1, 6, 1, -1, -1});
    auto hReco = hInvMass1.GetTH1("recon", 3, {1, 4, 1, -1, -1});
    
    auto clist2 = LoadXi1530ResultList("LHC18c6c_298", "Xi1530MB");
    auto clist3 = LoadXi1530ResultList("LHC18c6b_part1_298", "Xi1530MB");
    auto clist4 = LoadXi1530ResultList("LHC18c6b_part2_298", "Xi1530MB");
    auto clist5 = LoadXi1530ResultList("LHC18c6a_part1_298", "Xi1530MB");
    auto clist6 = LoadXi1530ResultList("LHC18c6a_part2_298", "Xi1530MB");
    

    auto hInvMass2 = BSTHnSparseHelper::Load("hInvMass", clist2);
    auto hInvMass3 = BSTHnSparseHelper::Load("hInvMass", clist3);
    auto hInvMass4 = BSTHnSparseHelper::Load("hInvMass", clist4);
    auto hInvMass5 = BSTHnSparseHelper::Load("hInvMass", clist5);
    auto hInvMass6 = BSTHnSparseHelper::Load("hInvMass", clist6);
    
    
    
    hInvMass2.SetBin("Cent", centbin);
    hInvMass3.SetBin("Cent", centbin);
    hInvMass4.SetBin("Cent", centbin);
    hInvMass5.SetBin("Cent", centbin);
    hInvMass6.SetBin("Cent", centbin);
    
    hInvMass2.SetBin("Pt", ptbin);
    hInvMass3.SetBin("Pt", ptbin);
    hInvMass4.SetBin("Pt", ptbin);
    hInvMass5.SetBin("Pt", ptbin);
    hInvMass6.SetBin("Pt", ptbin);

    
    hInput->Add(hInvMass2.GetTH1("input", 3, {1, 6, 1, -1, -1}));
    hInput->Add(hInvMass3.GetTH1("input", 3, {1, 6, 1, -1, -1}));
    hInput->Add(hInvMass4.GetTH1("input", 3, {1, 6, 1, -1, -1}));
    hInput->Add(hInvMass5.GetTH1("input", 3, {1, 6, 1, -1, -1}));
    hInput->Add(hInvMass6.GetTH1("input", 3, {1, 6, 1, -1, -1}));
    
    
    
    hReco->Add(hInvMass2.GetTH1("recon", 3, {1, 4, 1, -1, -1}));
    hReco->Add(hInvMass3.GetTH1("recon", 3, {1, 4, 1, -1, -1}));
    hReco->Add(hInvMass4.GetTH1("recon", 3, {1, 4, 1, -1, -1}));
    hReco->Add(hInvMass5.GetTH1("recon", 3, {1, 4, 1, -1, -1}));
    hReco->Add(hInvMass6.GetTH1("recon", 3, {1, 4, 1, -1, -1}));
    

    //-----------------------------------------------------------------
    TString finalfile = "AnalysisResults_Xi1530_systematic_0-100.root";
    TString workdirectory = "/Users/blim/alidock/Postprocessing/data/";
    TFile* inputfile = new TFile(finalfile.Data());
    TH1* hr = (TH1*)inputfile->Get(
        Form("hSpectra_%.2f_%.2f_sys", 0., 100.));
    TH1* hs = (TH1*)inputfile->Get(
        Form("hSpectra_%.2f_%.2f_stat", 0., 100.));

    TH1D* htot = (TH1D*)hr->Clone();
    for (Int_t ibin = 0; ibin < htot->GetNbinsX(); ibin++) {
        htot->SetBinError(
            ibin + 1,
            TMath::Sqrt(
                hr->GetBinError(ibin + 1) * hr->GetBinError(ibin + 1) +
                hs->GetBinError(ibin + 1) * hs->GetBinError(ibin + 1)));
    }

    TF1* myLevy = LevyTsallis("Levy", 1.5318, 15, 0.4, 1.5);
    htot->Fit(myLevy, "0q", "", 0.8, 8.8);

    TFile* f3=new TFile("output_file_new.root","RECREATE","HistoFile");
    

    ReweightEfficiency(htot,myLevy,hInput,hReco,f3,1);
    /*
    hInput->SetMinimum(1);
    hInput->SetMaximum(1e5);
    hInput->SetMarkerColor(2);
    hInput->SetLineColor(2);
    hInput->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hInput->GetXaxis()->SetRangeUser(0,15);
    hReco->SetMinimum(1);
    hReco->SetMaximum(1e5);
    hReco->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hReco->GetXaxis()->SetRangeUser(0,15);
    

    TCanvas* cCanvas = new TCanvas("cCanvas", "cCanvas", 960, 720);
    TGaxis::SetMaxDigits(3);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    cCanvas->SetTickx();
    cCanvas->SetTicky();
    cCanvas->SetTopMargin(0.05);
    cCanvas->SetLeftMargin(0.10);
    //cCanvas->SetBottomMargin(0.01);
    cCanvas->SetRightMargin(0.01);
    cCanvas->SetFillStyle(0);
    cCanvas->SetLogy(true);;
    cCanvas->Draw();

    hInput->Draw();
    hReco->Draw("same");

    htot->Draw();
    */
}