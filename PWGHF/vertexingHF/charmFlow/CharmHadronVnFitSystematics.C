//___________________________________________________________________________________//
// Brief: Macro for the fit systematics of vn analysis using the output of           //
// AliAnalysisTaskSECharmHadronvn                                                    //
// Main Function: CharmHadronVnFitSystematics                                        //
// Author: Fabrizio Grosa, fabrizio.grosa@cern.ch                                    //
//                                                                                   //
// Before running this macro make sure you have installed yaml-cpp on your laptop:   //
// MAC OSX --> brew install yaml-cpp                                                 //
// Ubuntu --> apt install yaml-cpp                                                   //
//___________________________________________________________________________________//

#if !defined (__CINT__) || defined (__CLING__)

#include <string>
#include <vector>

#include "yaml-cpp/yaml.h"

#include <TROOT.h>
#include <Riostream.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TDatabasePDG.h>
#include <TArrayD.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TNtuple.h>
#include <TLatex.h>
#include <TLine.h>

#include "AliHFInvMassFitter.h"
#include "AliHFVnVsMassFitter.h"
#include "AliVertexingHFUtils.h"
#include "AliAnalysisTaskSECharmHadronvn.h"

#endif

using namespace std;

enum fitType{kFreeSigma,kFixSigma,kSimFit,kBinCount};

//___________________________________________________________________________________//
//method prototypes
void CharmHadronVnFitSystematics(string cfgFileName, string refFileName, int refMeasType=kSimFit);
TH1D* ComputeEPresolution(double &resol, double &resolunc, TH3F *hDeltaPsiVsqnVsCentr[3], int harmonic, double qnmin, double qnmax);
TH1D* ComputeSPresolution(double &resol, double &resolunc, TH3F *hQiVsqnVsCentr[3], int harmonic, double qnmin, double qnmax);
void GetInOutOfPlaneInvMassHistos(THnSparseF *sparse, TH1F *&hInvMassInPlane, TH1F *&hInvMassOutOfPlane, int harmonic, double qnmin, double qnmax, double ptmin, double ptmax, bool applyML, double MLmin, double MLmax);
TH1F* GetFuncPhiVsMassHistos(THnSparseF* sparse, TString histoname, int iAxis, double qnmin, double qnmax, double ptmin, double ptmax, double massrebin, double resol, bool applyML, double MLmin, double MLmax);
float GetAveragePtInRange(float &averagePtUnc, THnSparseF *sparse, double qnmin, double qnmax, double ptmin, double ptmax, int bkgfunc, int sgnfunc, bool useRefl, TH1F* hMCRefl,
                          double SoverR, string reflopt, int meson, double massD, bool fixMeanSecP, bool fixSigmaSecP, float sigmaDplus, bool applyML, double MLmin, double MLmax);
double ComputeEPvn(double &vnunc, int harmonic, double nIn, double nInUnc, double nOut, double nOutUnc, double resol, double corr = 0.);
void ApplySelection(THnSparseF *sparse, int axisnum, double min, double max);
void ResetAxes(THnSparseF *sparse, int axisnum = -1);
TList* LoadTListFromTaskOutput(YAML::Node config);
bool LoadD0toKpiReflHistos(string reflFileName, int nPtBins, TH1F* hMCSgn[], TH1F* hMCRefl[]);
bool LoadDsDplusSigma(string sigmaFileName, int nPtBins, float *sigmaDplus);
void SetStyle();
void SetGraphStyle(TGraphAsymmErrors* graph, int color, int markerstyle, float markersize = 1.5, int linewidth = 2);
void SetHistoStyle(TH1* histo, int color, int markerstyle, float markersize = 1.5, int linewidth = 2);
TLatex* BuildTLatex(int color, int font, double fontsize);
void DivideCanvas(TCanvas* c, const unsigned int nPtBins);

//___________________________________________________________________________________//
//main function for vn analysis
void CharmHadronVnFitSystematics(string cfgFileName, string refFileName, int refMeasType) {

    SetStyle();

    //Load configs from yaml file
    YAML::Node config = YAML::LoadFile(cfgFileName);
    if (config.IsNull()) {
        cerr << "Yaml config file not found! Exit" << endl;
        return;
    }

    string mesonname = config["InputFile"]["Meson"].as<string>();
    string flowmethodname = config["InputFile"]["FlowMethod"].as<string>();

    int harmonic = config["AnalysisOptions"]["Harmonic"].as<int>();
    double qnmin = config["AnalysisOptions"]["qnMin"].as<double>();
    double qnmax = config["AnalysisOptions"]["qnMax"].as<double>();
    vector<double> PtMin = config["AnalysisOptions"]["PtMin"].as<vector<double> >();
    vector<double> PtMax = config["AnalysisOptions"]["PtMax"].as<vector<double> >();
    bool useRefl = static_cast<bool>(config["AnalysisOptions"]["IncludeReflections"].as<int>());
    string reflFileName = config["AnalysisOptions"]["ReflFileName"].as<string>();
    string reflopt = config["AnalysisOptions"]["ReflOpt"].as<string>();
    bool fixMeanSecP = static_cast<bool>(config["AnalysisOptions"]["FixMeanSecondPeak"].as<int>());
    bool fixSigmaSecP = static_cast<bool>(config["AnalysisOptions"]["FixSigmaSecondPeak"].as<int>());
    string sigmaFileName = config["AnalysisOptions"]["SigmaFileName"].as<string>();

    vector<double> MassMin = config["FitSystematicsOptions"]["MassMin"].as<vector<double> >();
    vector<double> MassMax = config["FitSystematicsOptions"]["MassMax"].as<vector<double> >();
    vector<int> Rebin = config["FitSystematicsOptions"]["Rebin"].as<vector<int> >();
    vector<string> sBkgFunc = config["FitSystematicsOptions"]["BkgFunc"].as<vector<string> >();
    vector<string> sSgnFunc = config["FitSystematicsOptions"]["SgnFunc"].as<vector<string> >();
    vector<string> sVnBkgFunc = config["FitSystematicsOptions"]["VnBkgFunc"].as<vector<string> >();
    double maxRedChi2 = config["FitSystematicsOptions"]["MaxRedChi2"].as<double>();

    int flowmethod = -1.;
    if(flowmethodname=="EP")
        flowmethod = AliAnalysisTaskSECharmHadronvn::kEP;
    else if(flowmethodname=="EvShapeEP")
        flowmethod = AliAnalysisTaskSECharmHadronvn::kEvShapeEP;
    else if(flowmethodname=="SP")
        flowmethod = AliAnalysisTaskSECharmHadronvn::kSP;
    else if(flowmethodname=="EvShapeSP")
        flowmethod = AliAnalysisTaskSECharmHadronvn::kEvShapeSP;
    else if(flowmethodname=="EPVsMass")
        flowmethod = AliAnalysisTaskSECharmHadronvn::kEPVsMass;
    else if(flowmethodname=="EvShapeEPVsMass")
        flowmethod = AliAnalysisTaskSECharmHadronvn::kEvShapeEPVsMass;

    int meson = -1.;
    double massD = -1., massDplus = -1.;
    int pdgcode = 0;
    if(mesonname=="Dplus") {
        meson = AliAnalysisTaskSECharmHadronvn::kDplustoKpipi;
        pdgcode = 411;
        massD = TDatabasePDG::Instance()->GetParticle(pdgcode)->Mass();
    }
    else if(mesonname=="Dzero") {
        meson = AliAnalysisTaskSECharmHadronvn::kD0toKpi;
        pdgcode = 421;
        massD = TDatabasePDG::Instance()->GetParticle(pdgcode)->Mass();
    }
    else if(mesonname=="Dstar") {
        meson = AliAnalysisTaskSECharmHadronvn::kDstartoKpipi;
        pdgcode = 413;
        massD = (TDatabasePDG::Instance()->GetParticle(413)->Mass() - TDatabasePDG::Instance()->GetParticle(421)->Mass());
    }
    else if(mesonname=="Ds") {
        meson = AliAnalysisTaskSECharmHadronvn::kDstoKKpi;
        pdgcode = 431;
        massD = TDatabasePDG::Instance()->GetParticle(pdgcode)->Mass();
        massDplus = TDatabasePDG::Instance()->GetParticle(411)->Mass();
    }

    if(useRefl && meson!=AliAnalysisTaskSECharmHadronvn::kD0toKpi) useRefl = false;

    const unsigned int nPtBins = PtMin.size();
    double PtLims[nPtBins+1];
    for(unsigned int iPt=0; iPt<nPtBins; iPt++)
        PtLims[iPt] = PtMin[iPt];
    PtLims[nPtBins] = PtMax[nPtBins-1];

    //arguments for ML application
    bool doMLsel = false;
    vector<double> CutValuesMLmin = {};
    vector<double> CutValuesMLmax = {};
    if (config["AnalysisOptions"]["MLSelection"]) {
        doMLsel = static_cast<bool>(config["AnalysisOptions"]["MLSelection"]["ApplyML"].as<int>());
        CutValuesMLmin = config["AnalysisOptions"]["MLSelection"]["CutValuesMin"].as<vector<double> >();
        CutValuesMLmax = config["AnalysisOptions"]["MLSelection"]["CutValuesMax"].as<vector<double> >();
    }

    //Load input file
    TList* list = LoadTListFromTaskOutput(config);
    if(!list) return;
    TH3F* hDeltaPsiVsqnVsCentr[3] = {NULL, NULL, NULL};
    for(int iResHist=0; iResHist<3; iResHist++) {
        if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEPVsMass || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEPVsMass)
            hDeltaPsiVsqnVsCentr[iResHist] = static_cast<TH3F*>(list->FindObject(Form("fHistEvPlaneReso%d",iResHist+1)));
        else if(flowmethod==AliAnalysisTaskSECharmHadronvn::kSP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeSP)
            hDeltaPsiVsqnVsCentr[iResHist] = static_cast<TH3F*>(list->FindObject(Form("hScalProdQnVectors%d",iResHist+1)));
    }
    THnSparseF* sMassVsPtVsPhiVsCentrVsqn = static_cast<THnSparseF*>(list->FindObject("fHistMassPtPhiqnCentr"));

    //Load reflections file
    TH1F* hMCSgn[nPtBins];
    TH1F* hMCRefl[nPtBins];
    if(useRefl)
       useRefl = LoadD0toKpiReflHistos(reflFileName, nPtBins, hMCSgn, hMCRefl);

    //Load D+ peak width
    float DplusSigma[nPtBins];
    if(fixSigmaSecP){
        bool funcOut = LoadDsDplusSigma(sigmaFileName, nPtBins, DplusSigma);
        if(!funcOut) {
            cout<<"Problem in loading of MC D+ sigma, switching to default value: 0.008";
            for(unsigned int iPt = 0; iPt < nPtBins; iPt++)
                DplusSigma[iPt] = 0.008;
        }
    }else{
        for(unsigned int iPt = 0; iPt < nPtBins; iPt++)
            DplusSigma[iPt] = 0.008;
    }

    //Load ref file
    TFile *reffile = TFile::Open(refFileName.data());
    if(!reffile)
        return;
    TString refgraphname = "gvnSimFit";
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
        if(refMeasType==kFreeSigma)
            refgraphname = "gvnFreeSigma";
        else if(refMeasType==kFixSigma)
            refgraphname = "gvnFixSigma";
        else if(refMeasType==kSimFit)
            refgraphname = "gvnSimFit";
        else if(refMeasType==kBinCount)
            refgraphname = "gvnBinCount";
    }
    TGraphAsymmErrors* gVnRef = static_cast<TGraphAsymmErrors*>(reffile->Get(refgraphname.Data()));
    if(!gVnRef) {
        cerr << "Error: ref graph not found! Exit" << endl;
        return;
    }
    reffile->Close();

    //Get EP/SP resolution
    double resol=-1., resolunc=-1.;
    TH1D *hResolVsCentr = NULL;
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEPVsMass || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEPVsMass)
        hResolVsCentr = ComputeEPresolution(resol,resolunc,hDeltaPsiVsqnVsCentr,harmonic,qnmin,qnmax);
    else if(flowmethod==AliAnalysisTaskSECharmHadronvn::kSP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeSP)
        hResolVsCentr = ComputeSPresolution(resol,resolunc,hDeltaPsiVsqnVsCentr,harmonic,qnmin,qnmax);
    TH1D* hResolCentrInt = new TH1D("hResolCentrInt",Form("centrality(%%);;%s",hResolVsCentr->GetYaxis()->GetTitle()),1,hResolVsCentr->GetBinLowEdge(1),hResolVsCentr->GetBinLowEdge(hResolVsCentr->GetNbinsX())+hResolVsCentr->GetBinWidth(1));
    hResolCentrInt->SetBinContent(1,resol);
    hResolCentrInt->SetBinError(1,resolunc);
    SetHistoStyle(hResolCentrInt,kRed,kOpenSquare);

    //TNtuple for multi-trial
    string outputdir = config["OutputDir"]["FitSystematics"].as<string>();
    TFile outFile(Form("%s/FitSystematics_%sv%d_%s_q%d_%0.f-%0.f.root",outputdir.data(),mesonname.data(),harmonic,flowmethodname.data(),harmonic,qnmin,qnmax),"recreate");
    TString vars = "";
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP)
        vars = "PtMin:PtMax:PtCent:PtCentUnc:FitType:BinWidth:MassMin:MassMax:SgnFunc:BkgFunc:RawYieldInPlane:RawYieldOutOfPlane:RawYieldInPlaneUnc:RawYieldOutOfPlaneUnc:MeanInPlane:MeanOutOfPlane:MeanInPlaneUnc:MeanOutOfPlaneUnc:SigmaInPlane:SigmaOutOfPlane:SigmaInPlaneUnc:SigmaOutOfPlaneUnc:RedChi2InPlane:RedChi2OutOfPlane:vn:vnUnc:qnMin:qnMax:TrialNo";
    else
        vars = "PtMin:PtMax:PtCent:PtCentUnc:BinWidth:MassMin:MassMax:SgnFunc:BkgFunc:VnBkgFunc:RawYield:RawYieldUnc:Mean:MeanUnc:Sigma:SigmaUnc:RedChi2:vn:vnUnc:qnMin:qnMax:TrialNo";
    TNtuple* multiTrialNtuple = new TNtuple("multiTrialNtuple","multiTrialNtuple",vars.Data());

    //Histos
    TH1F* hVnResFreeSigma[nPtBins];
    TH1F* hVnResFixSigma[nPtBins];
    TH1F* hVnResBinCounting[nPtBins];
    TH1F* hVnResSimFit[nPtBins];

    TH1F* hVnResMeanVsMethodFreeSigma[nPtBins];
    TH1F* hVnResMeanVsMethodFixSigma[nPtBins];
    TH1F* hVnResMeanVsMethodBinCounting[nPtBins];
    TH1F* hVnResMeanVsMethodSimFit[nPtBins];

    TGraphAsymmErrors* gVnVsTrialFreeSigma[nPtBins];
    TGraphAsymmErrors* gVnVsTrialFixSigma[nPtBins];
    TGraphAsymmErrors* gVnVsTrialBinCounting[nPtBins];
    TGraphAsymmErrors* gVnVsTrialSimFit[nPtBins];

    TLine* lineStatRed[2][nPtBins];

    //Define fitters
    AliHFInvMassFitter* massfitterFreeSigma[2][nPtBins];
    AliHFInvMassFitter* massfitterFixSigma[2][nPtBins];
    AliHFInvMassFitter* massfitterSimFit[2][nPtBins];

    AliHFVnVsMassFitter* vnvsmassfitter[nPtBins];

    //Multi-trial
    unsigned int nVnBkgFunc = (flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) ? 1 : sVnBkgFunc.size();
    int SgnFunc = -1, BkgFunc = -1, VnBkgFunc = -1;
    TH1F *hInvMassInt = NULL, *hVnVsMassToFit = NULL;
    TH1F *hInvMassDeltaPhi[2] = {NULL, NULL}, *hInvMassDeltaPhiToFit[2] = {NULL, NULL};

    int iTrial = 0;
    TString v2measname = (flowmethod==AliAnalysisTaskSECharmHadronvn::kSP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeSP) ? "SP" : "EP";
    for(unsigned int iPt=0; iPt<nPtBins; iPt++) {

        iTrial = 0;

        //Get ref vn
        double vnRef=-1., PtRed=-1.;
        gVnRef->GetPoint(iPt,PtRed,vnRef);
        double statunc = gVnRef->GetErrorYlow(iPt);
        lineStatRed[0][iPt] = new TLine(0.5,statunc,4.5,statunc);
        lineStatRed[1][iPt] = new TLine(0.5,-statunc,4.5,-statunc);
        for(int iLine=0; iLine<2; iLine++) {
            lineStatRed[iLine][iPt]->SetLineWidth(2);
            lineStatRed[iLine][iPt]->SetLineColor(kBlack);
            lineStatRed[iLine][iPt]->SetLineStyle(9);
        }

        //define output histos/graphs
        hVnResFreeSigma[iPt] = new TH1F(Form("hVnResFreeSigma_%d",iPt),Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c};#it{v}_{%d}-#it{v}_{%d}^{ref}{%s};Entries",PtMin[iPt],PtMax[iPt],harmonic,harmonic,v2measname.Data()),100,-0.5,0.5);
        hVnResFixSigma[iPt] = new TH1F(Form("hVnResFixSigma_%d",iPt),Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c};#it{v}_{%d}-#it{v}_{%d}^{ref}{%s};Entries",PtMin[iPt],PtMax[iPt],harmonic,harmonic,v2measname.Data()),100,-0.5,0.5);
        hVnResBinCounting[iPt] = new TH1F(Form("hVnResBinCounting_%d",iPt),Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c};#it{v}_{%d}-#it{v}_{%d}^{ref}{%s};Entries",PtMin[iPt],PtMax[iPt],harmonic,harmonic,v2measname.Data()),100,-0.5,0.5);
        hVnResSimFit[iPt] = new TH1F(Form("hVnResSimFit_%d",iPt),Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c};#it{v}_{%d}-#it{v}_{%d}^{ref}{%s};Entries",PtMin[iPt],PtMax[iPt],harmonic,harmonic,v2measname.Data()),100,-0.5,0.5);
        hVnResFreeSigma[iPt]->SetDirectory(0);
        hVnResFixSigma[iPt]->SetDirectory(0);
        hVnResBinCounting[iPt]->SetDirectory(0);
        hVnResSimFit[iPt]->SetDirectory(0);
        SetHistoStyle(hVnResFreeSigma[iPt],kRed,kFullCircle);
        SetHistoStyle(hVnResFixSigma[iPt],kBlack,kFullCircle);
        SetHistoStyle(hVnResBinCounting[iPt],kGreen+2,kFullCircle);
        SetHistoStyle(hVnResSimFit[iPt],kBlue,kFullCircle);

        hVnResMeanVsMethodFreeSigma[iPt] = new TH1F(Form("hVnResMeanVsMethodFreeSigma_%d",iPt),Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c};;Mean #pm RMS #it{v}_{%d}-#it{v}_{%d}^{ref}{%s}",PtMin[iPt],PtMax[iPt],harmonic,harmonic,v2measname.Data()),4,0.5,4.5);
        hVnResMeanVsMethodFixSigma[iPt] = new TH1F(Form("hVnResMeanVsMethodFreeSigma_%d",iPt),Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c};;Mean #pm RMS #it{v}_{%d}-#it{v}_{%d}^{ref}{%s}",PtMin[iPt],PtMax[iPt],harmonic,harmonic,v2measname.Data()),4,0.5,4.5);
        hVnResMeanVsMethodBinCounting[iPt] = new TH1F(Form("hVnResMeanVsMethodFreeSigma_%d",iPt),Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c};;Mean #pm RMS #it{v}_{%d}-#it{v}_{%d}^{ref}{%s}",PtMin[iPt],PtMax[iPt],harmonic,harmonic,v2measname.Data()),4,0.5,4.5);
        hVnResMeanVsMethodSimFit[iPt] = new TH1F(Form("hVnResMeanVsMethodFreeSigma_%d",iPt),Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c};;Mean #pm RMS #it{v}_{%d}-#it{v}_{%d}^{ref}{%s}",PtMin[iPt],PtMax[iPt],harmonic,harmonic,v2measname.Data()),4,0.5,4.5);
        hVnResMeanVsMethodFreeSigma[iPt]->SetDirectory(0);
        hVnResMeanVsMethodFixSigma[iPt]->SetDirectory(0);
        hVnResMeanVsMethodBinCounting[iPt]->SetDirectory(0);
        hVnResMeanVsMethodSimFit[iPt]->SetDirectory(0);
        SetHistoStyle(hVnResMeanVsMethodFreeSigma[iPt],kRed,kFullCircle);
        SetHistoStyle(hVnResMeanVsMethodFixSigma[iPt],kBlack,kFullTriangleUp);
        SetHistoStyle(hVnResMeanVsMethodBinCounting[iPt],kGreen+2,kFullSquare);
        SetHistoStyle(hVnResMeanVsMethodSimFit[iPt],kBlue,kFullDiamond,2.,2);
        hVnResMeanVsMethodSimFit[iPt]->GetXaxis()->SetBinLabel(1,"Simultaneus Fit");
        hVnResMeanVsMethodFreeSigma[iPt]->GetXaxis()->SetBinLabel(2,"Free Sigma");
        hVnResMeanVsMethodFixSigma[iPt]->GetXaxis()->SetBinLabel(3,"Fix Sigma");
        hVnResMeanVsMethodBinCounting[iPt]->GetXaxis()->SetBinLabel(4,"Bin Counting");

        gVnVsTrialFreeSigma[iPt] = new TGraphAsymmErrors(0);
        gVnVsTrialFixSigma[iPt] = new TGraphAsymmErrors(0);
        gVnVsTrialBinCounting[iPt] = new TGraphAsymmErrors(0);
        gVnVsTrialSimFit[iPt] = new TGraphAsymmErrors(0);
        gVnVsTrialFreeSigma[iPt]->SetName(Form("gVnVsTrialFreeSigma_%d",iPt));
        gVnVsTrialFixSigma[iPt]->SetName(Form("gVnVsTrialFixSigma_%d",iPt));
        gVnVsTrialBinCounting[iPt]->SetName(Form("gVnVsTrialBinCounting_%d",iPt));
        gVnVsTrialSimFit[iPt]->SetName(Form("gVnVsTrialSimFit_%d",iPt));
        SetGraphStyle(gVnVsTrialFreeSigma[iPt],kRed,kFullCircle,0.8,1);
        SetGraphStyle(gVnVsTrialFixSigma[iPt],kGray+2,kFullTriangleUp,0.8,1);
        SetGraphStyle(gVnVsTrialBinCounting[iPt],kGreen+2,kFullSquare,0.8,1);
        SetGraphStyle(gVnVsTrialSimFit[iPt],kBlue,kFullDiamond,1.,1);

        //compute average pT
        double SoverR = 0.;
        if(useRefl)
            SoverR=(hMCRefl[iPt]->Integral(hMCRefl[iPt]->FindBin(1.65),hMCRefl[iPt]->FindBin(2.15)))/(hMCSgn[iPt]->Integral(hMCSgn[iPt]->FindBin(1.65),hMCSgn[iPt]->FindBin(2.15)));

        float averagePtUnc = -1.;
        int bkgfunc = (meson==AliAnalysisTaskSECharmHadronvn::kDstartoKpipi) ? AliHFInvMassFitter::kPowEx : AliHFInvMassFitter::kExpo;
        float averagePt = GetAveragePtInRange(averagePtUnc, sMassVsPtVsPhiVsCentrVsqn, qnmin, qnmax, PtMin[iPt], PtMax[iPt], bkgfunc, SgnFunc, 
                                              useRefl, hMCRefl[iPt], SoverR, reflopt, meson, massD, fixMeanSecP, fixSigmaSecP, DplusSigma[iPt], 
                                              doMLsel, CutValuesMLmin[iPt], CutValuesMLmax[iPt]);

        //get inv-mass histos
        ApplySelection(sMassVsPtVsPhiVsCentrVsqn,1,PtMin[iPt],PtMax[iPt]);
        if(doMLsel)
            ApplySelection(sMassVsPtVsPhiVsCentrVsqn,9,CutValuesMLmin[iPt], CutValuesMLmax[iPt]);

        hInvMassInt = reinterpret_cast<TH1F*>(sMassVsPtVsPhiVsCentrVsqn->Projection(0));
        ResetAxes(sMassVsPtVsPhiVsCentrVsqn,1);
        if(doMLsel)
            ResetAxes(sMassVsPtVsPhiVsCentrVsqn,9);

        if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP)
            GetInOutOfPlaneInvMassHistos(sMassVsPtVsPhiVsCentrVsqn, hInvMassDeltaPhi[0], hInvMassDeltaPhi[1], harmonic, qnmin, qnmax, PtMin[iPt], PtMax[iPt], 
                                         doMLsel, CutValuesMLmin[iPt], CutValuesMLmax[iPt]);

        for(unsigned int iSgnFunc=0; iSgnFunc<sSgnFunc.size(); iSgnFunc++) {

            if(sSgnFunc[iSgnFunc]=="kGaus")
                SgnFunc = AliHFInvMassFitter::kGaus;
            else if(sSgnFunc[iSgnFunc]=="k2Gaus")
                SgnFunc = AliHFInvMassFitter::k2Gaus;
            else if(sSgnFunc[iSgnFunc]=="k2GausSigmaRatioPar")
                SgnFunc = AliHFInvMassFitter::k2GausSigmaRatioPar;

            for(unsigned int iBkgFunc=0; iBkgFunc<sBkgFunc.size(); iBkgFunc++) {

                if(sBkgFunc[iBkgFunc]=="kExpo")
                    BkgFunc = AliHFInvMassFitter::kExpo;
                else if(sBkgFunc[iBkgFunc]=="kLin")
                    BkgFunc = AliHFInvMassFitter::kLin;
                else if(sBkgFunc[iBkgFunc]=="kPol2")
                    BkgFunc = AliHFInvMassFitter::kPol2;
                else if(sBkgFunc[iBkgFunc]=="kNoBk")
                    BkgFunc = AliHFInvMassFitter::kNoBk;
                else if(sBkgFunc[iBkgFunc]=="kPow")
                    BkgFunc = AliHFInvMassFitter::kPow;
                else if(sBkgFunc[iBkgFunc]=="kPowEx")
                    BkgFunc = AliHFInvMassFitter::kPowEx;

                for(unsigned int iVnBkgFunc=0; iVnBkgFunc<sVnBkgFunc.size(); iVnBkgFunc++) {

                    if(sVnBkgFunc[iVnBkgFunc]=="kLin")
                        VnBkgFunc = AliHFVnVsMassFitter::kLin;
                    else if(sVnBkgFunc[iVnBkgFunc]=="kPol2")
                        VnBkgFunc = AliHFVnVsMassFitter::kPol2;
                    else if(sVnBkgFunc[iVnBkgFunc]=="kExpo")
                        VnBkgFunc = AliHFVnVsMassFitter::kExpo;

                    for(unsigned int iRebin=0; iRebin<Rebin.size(); iRebin++) {

                        //rebin histos
                        TH1F* hInvMassIntToFit = static_cast<TH1F*>(hInvMassInt->Clone("hInvMassIntToFit"));
                        hInvMassIntToFit->Rebin(Rebin[iRebin]);

                        if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
                            hInvMassDeltaPhiToFit[0] = static_cast<TH1F*>(hInvMassDeltaPhi[0]->Clone("hInvMassDeltaPhiToFit0"));
                            hInvMassDeltaPhiToFit[0]->Rebin(Rebin[iRebin]);
                            hInvMassDeltaPhiToFit[1] = static_cast<TH1F*>(hInvMassDeltaPhi[1]->Clone("hInvMassDeltaPhiToFit1"));
                            hInvMassDeltaPhiToFit[1]->Rebin(Rebin[iRebin]);
                        }
                        else
                            hVnVsMassToFit = GetFuncPhiVsMassHistos(sMassVsPtVsPhiVsCentrVsqn, Form("hVnVsMassToFit_%d",iPt), 2, qnmin, qnmax, PtMin[iPt], PtMax[iPt], 
                                                                    Rebin[iRebin], resol, doMLsel, CutValuesMLmin[iPt], CutValuesMLmax[iPt]);

                        for(unsigned int iMassMin=0; iMassMin<MassMin.size(); iMassMin++) {
                            for(unsigned int iMassMax=0; iMassMax<MassMax.size(); iMassMax++) {

                                if(useRefl)
                                    SoverR=(hMCRefl[iPt]->Integral(hMCRefl[iPt]->FindBin(MassMin[iMassMin]*1.0001),hMCRefl[iPt]->FindBin(MassMax[iMassMax]*0.9999)))/(hMCSgn[iPt]->Integral(hMCSgn[iPt]->FindBin(MassMin[iMassMin]*1.0001),hMCSgn[iPt]->FindBin(MassMax[iMassMax]*0.9999)));

                                //fit and fill ntuple
                                AliHFInvMassFitter massfitterInt(hInvMassIntToFit,MassMin[iMassMin],MassMax[iMassMax],BkgFunc,SgnFunc);
                                massfitterInt.SetUseLikelihoodFit();
                                massfitterInt.SetInitialGaussianMean(massD);
                                massfitterInt.SetInitialGaussianSigma(0.010);
                                if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi)
                                    massfitterInt.IncludeSecondGausPeak(massDplus,fixMeanSecP,DplusSigma[iPt],fixSigmaSecP);
                                else if(meson==AliAnalysisTaskSECharmHadronvn::kDstartoKpipi)
                                    massfitterInt.SetInitialGaussianSigma(0.001);
                                if(useRefl) {
                                    massfitterInt.SetTemplateReflections(hMCRefl[iPt],reflopt,MassMin[iMassMin],MassMax[iMassMax]);
                                    massfitterInt.SetFixReflOverS(SoverR);
                                }
                                massfitterInt.MassFitter(false);

                                if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {

                                    AliHFInvMassFitter *massfitterFreeSigma[2];
                                    AliHFInvMassFitter *massfitterFixSigma[2];
                                    AliHFInvMassFitter *massfitterSimFit[2];
                                    double rawYieldsFreeSigma[2], rawYieldsBinCount[2], rawYieldsFixSigma[2], rawYieldsSimFit[2], rawYieldsFreeSigmaUnc[2], rawYieldsBinCountUnc[2], rawYieldsFixSigmaUnc[2], rawYieldsSimFitUnc[2];
                                    int fitFreeSigmaStatus[2], fitFixSigmaStatus[2];
                                    ROOT::Fit::FitResult resultSimFit;

                                    for(int iDeltaPhi=0; iDeltaPhi<2; iDeltaPhi++) {
                                        //free sigma
                                        massfitterFreeSigma[iDeltaPhi] = new AliHFInvMassFitter(hInvMassDeltaPhiToFit[iDeltaPhi],MassMin[iMassMin],MassMax[iMassMax],BkgFunc,SgnFunc);
                                        massfitterFreeSigma[iDeltaPhi]->SetUseLikelihoodFit();
                                        massfitterFreeSigma[iDeltaPhi]->SetInitialGaussianMean(massD);
                                        massfitterFreeSigma[iDeltaPhi]->SetInitialGaussianSigma(0.010);
                                        if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi)
                                            massfitterFreeSigma[iDeltaPhi]->IncludeSecondGausPeak(massDplus,fixMeanSecP,DplusSigma[iPt],fixSigmaSecP);
                                        else if(meson==AliAnalysisTaskSECharmHadronvn::kDstartoKpipi)
                                            massfitterInt.SetInitialGaussianSigma(0.001);
                                        if(useRefl) {
                                            massfitterFreeSigma[iDeltaPhi]->SetTemplateReflections(hMCRefl[iPt],reflopt,MassMin[iMassMin],MassMax[iMassMax]);
                                            massfitterFreeSigma[iDeltaPhi]->SetFixReflOverS(SoverR);
                                        }
                                        fitFreeSigmaStatus[iDeltaPhi] = massfitterFreeSigma[iDeltaPhi]->MassFitter(false);
                                        if(fitFreeSigmaStatus[iDeltaPhi] == 1) {
                                            rawYieldsFreeSigma[iDeltaPhi] = massfitterFreeSigma[iDeltaPhi]->GetRawYield();
                                            rawYieldsFreeSigmaUnc[iDeltaPhi] = massfitterFreeSigma[iDeltaPhi]->GetRawYieldError();
                                        }else {
                                            rawYieldsFreeSigma[iDeltaPhi] = 1.;
                                            rawYieldsFreeSigmaUnc[iDeltaPhi] = 1.;
                                        }

                                        //bin counting
                                        if(fitFreeSigmaStatus[iDeltaPhi] == 1)
                                            rawYieldsBinCount[iDeltaPhi] = massfitterFreeSigma[iDeltaPhi]->GetRawYieldBinCounting(rawYieldsBinCountUnc[iDeltaPhi],3,1,pdgcode);
                                        else
                                            rawYieldsBinCount[iDeltaPhi] = 1.;

                                        //fix sigma
                                        massfitterFixSigma[iDeltaPhi] = new AliHFInvMassFitter(hInvMassDeltaPhiToFit[iDeltaPhi],MassMin[iMassMin],MassMax[iMassMax],BkgFunc,SgnFunc);
                                        massfitterFixSigma[iDeltaPhi]->SetUseLikelihoodFit();
                                        massfitterFixSigma[iDeltaPhi]->SetInitialGaussianMean(massD);
                                        massfitterFixSigma[iDeltaPhi]->SetFixGaussianSigma(massfitterInt.GetSigma());
                                        if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi)
                                            massfitterFixSigma[iDeltaPhi]->IncludeSecondGausPeak(massDplus,fixMeanSecP,DplusSigma[iPt],fixSigmaSecP);
                                        if(useRefl) {
                                            massfitterFixSigma[iDeltaPhi]->SetTemplateReflections(hMCRefl[iPt],reflopt,MassMin[iMassMin],MassMax[iMassMax]);
                                            massfitterFixSigma[iDeltaPhi]->SetFixReflOverS(SoverR);
                                        }
                                        fitFixSigmaStatus[iDeltaPhi] = massfitterFixSigma[iDeltaPhi]->MassFitter(false);
                                        if(fitFixSigmaStatus[iDeltaPhi] == 1) {
                                            rawYieldsFixSigma[iDeltaPhi] = massfitterFixSigma[iDeltaPhi]->GetRawYield();
                                            rawYieldsFixSigmaUnc[iDeltaPhi] = massfitterFixSigma[iDeltaPhi]->GetRawYieldError();
                                        }else {
                                            rawYieldsFixSigma[iDeltaPhi] = 1.;
                                            rawYieldsFixSigmaUnc[iDeltaPhi] = 1.;
                                        }
                                        //simultaneus fit
                                        massfitterSimFit[iDeltaPhi] = new AliHFInvMassFitter(hInvMassDeltaPhiToFit[iDeltaPhi],MassMin[iMassMin],MassMax[iMassMax],BkgFunc,SgnFunc);
                                        if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi)
                                            massfitterSimFit[iDeltaPhi]->IncludeSecondGausPeak(massDplus,fixMeanSecP,DplusSigma[iPt],fixSigmaSecP);
                                        if(useRefl) {
                                            massfitterSimFit[iDeltaPhi]->SetTemplateReflections(hMCRefl[iPt],reflopt,MassMin[iMassMin],MassMax[iMassMax]);
                                            massfitterSimFit[iDeltaPhi]->SetFixReflOverS(SoverR);
                                        }
                                    }

                                    if(fitFreeSigmaStatus[0] == 1 && fitFreeSigmaStatus[1] == 1 && fitFixSigmaStatus[0] == 1 && fitFixSigmaStatus[1] == 1) {
                                        vector<unsigned int> commonpars = {static_cast<unsigned int>(massfitterFreeSigma[0]->GetBackgroundFullRangeFunc()->GetNpar())+1,
                                                                           static_cast<unsigned int>(massfitterFreeSigma[0]->GetBackgroundFullRangeFunc()->GetNpar())+2};

                                        if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi){
                                            commonpars.push_back(static_cast<unsigned int>(massfitterFreeSigma[0]->GetBackgroundFullRangeFunc()->GetNpar())+4);
                                            commonpars.push_back(static_cast<unsigned int>(massfitterFreeSigma[0]->GetBackgroundFullRangeFunc()->GetNpar())+5);
                                        }

                                        resultSimFit = AliVertexingHFUtils::DoInPlaneOutOfPlaneSimultaneusFit(massfitterSimFit[0], massfitterSimFit[1], hInvMassDeltaPhiToFit[0], hInvMassDeltaPhiToFit[1], MassMin[iMassMin], MassMax[iMassMax], massD, commonpars);
                                        for(int iDeltaPhi=0; iDeltaPhi<2; iDeltaPhi++) {
                                            rawYieldsSimFit[iDeltaPhi] = massfitterSimFit[iDeltaPhi]->GetSignalFunc()->GetParameter(0) / hInvMassDeltaPhiToFit[0]->GetBinWidth(1);
                                            rawYieldsSimFitUnc[iDeltaPhi] = massfitterSimFit[iDeltaPhi]->GetSignalFunc()->GetParError(0) / hInvMassDeltaPhiToFit[0]->GetBinWidth(1);
                                        }

                                        double vnFreeSigma=0., vnBinCount=0., vnFixSigma=0., vnSimFit=0., vnFreeSigmaUnc=0., vnBinCountUnc=0., vnFixSigmaUnc=0., vnSimFitUnc=0.;
                                        vnFreeSigma = ComputeEPvn(vnFreeSigmaUnc,harmonic,rawYieldsFreeSigma[0],rawYieldsFreeSigmaUnc[0],rawYieldsFreeSigma[1],rawYieldsFreeSigmaUnc[1],resol);
                                        vnBinCount = ComputeEPvn(vnBinCountUnc,harmonic,rawYieldsBinCount[0],rawYieldsBinCountUnc[0],rawYieldsBinCount[1],rawYieldsBinCountUnc[1],resol);
                                        vnFixSigma = ComputeEPvn(vnFixSigmaUnc,harmonic,rawYieldsFixSigma[0],rawYieldsFixSigmaUnc[0],rawYieldsFixSigma[1],rawYieldsFixSigmaUnc[1],resol);
                                        int posRawYieldPar = massfitterFreeSigma[0]->GetBackgroundFullRangeFunc()->GetNpar();
                                        int nTotPars = massfitterFreeSigma[0]->GetMassFunc()->GetNpar();
                                        vnSimFit = ComputeEPvn(vnSimFitUnc,harmonic,rawYieldsSimFit[0],rawYieldsSimFitUnc[0],rawYieldsSimFit[1],rawYieldsSimFitUnc[1],resol,resultSimFit.Correlation(posRawYieldPar,posRawYieldPar+nTotPars));

                                        float array4ntupleFreeSigma[29] = {(float)PtMin[iPt],(float)PtMax[iPt],(float)averagePt,(float)averagePtUnc,(float)kFreeSigma,(float)hInvMassDeltaPhiToFit[0]->GetBinWidth(1),(float)MassMin[iMassMin],(float)MassMax[iMassMax],(float)SgnFunc,(float)BkgFunc,(float)rawYieldsFreeSigma[0],(float)rawYieldsFreeSigma[1],(float)rawYieldsFreeSigmaUnc[0],(float)rawYieldsFreeSigmaUnc[1],(float)massfitterFreeSigma[0]->GetMean(),(float)massfitterFreeSigma[1]->GetMean(),(float)massfitterFreeSigma[0]->GetMeanUncertainty(),(float)massfitterFreeSigma[1]->GetMeanUncertainty(),(float)massfitterFreeSigma[0]->GetSigma(),(float)massfitterFreeSigma[1]->GetSigma(),(float)massfitterFreeSigma[0]->GetSigmaUncertainty(),(float)massfitterFreeSigma[1]->GetSigmaUncertainty(),(float)massfitterFreeSigma[0]->GetReducedChiSquare(),(float)massfitterFreeSigma[1]->GetReducedChiSquare(),(float)vnFreeSigma,(float)vnFreeSigmaUnc,(float)qnmin,(float)qnmax,(float)iTrial};
                                        float array4ntupleFixSigma[29] = {(float)PtMin[iPt],(float)PtMax[iPt],(float)averagePt,(float)averagePtUnc,(float)kFixSigma,(float)hInvMassDeltaPhiToFit[0]->GetBinWidth(1),(float)MassMin[iMassMin],(float)MassMax[iMassMax],(float)SgnFunc,(float)BkgFunc,(float)rawYieldsFixSigma[0],(float)rawYieldsFixSigma[1],(float)rawYieldsFixSigmaUnc[0],(float)rawYieldsFixSigmaUnc[1],(float)massfitterFixSigma[0]->GetMean(),(float)massfitterFixSigma[1]->GetMean(),(float)massfitterFixSigma[0]->GetMeanUncertainty(),(float)massfitterFixSigma[1]->GetMeanUncertainty(),(float)massfitterInt.GetSigma(),(float)massfitterInt.GetSigma(),(float)massfitterInt.GetSigmaUncertainty(),(float)massfitterInt.GetSigmaUncertainty(),(float)massfitterFixSigma[0]->GetReducedChiSquare(),(float)massfitterFixSigma[1]->GetReducedChiSquare(),(float)vnFixSigma,(float)vnFixSigmaUnc,(float)qnmin,(float)qnmax,(float)iTrial};
                                        float array4ntupleSimFit[29] = {(float)PtMin[iPt],(float)PtMax[iPt],(float)averagePt,(float)averagePtUnc,(float)kSimFit,(float)hInvMassDeltaPhiToFit[0]->GetBinWidth(1),(float)MassMin[iMassMin],(float)MassMax[iMassMax],(float)SgnFunc,(float)BkgFunc,(float)rawYieldsSimFit[0],(float)rawYieldsSimFit[1],(float)rawYieldsSimFitUnc[0],(float)rawYieldsSimFitUnc[1],(float)massfitterSimFit[0]->GetSignalFunc()->GetParameter(1),(float)massfitterSimFit[0]->GetSignalFunc()->GetParameter(1),(float)massfitterSimFit[1]->GetSignalFunc()->GetParError(1),(float)massfitterSimFit[1]->GetSignalFunc()->GetParError(1),(float)massfitterSimFit[0]->GetSignalFunc()->GetParameter(2),(float)massfitterSimFit[0]->GetSignalFunc()->GetParameter(2),(float)massfitterSimFit[1]->GetSignalFunc()->GetParError(2),(float)massfitterSimFit[1]->GetSignalFunc()->GetParError(2),(float)resultSimFit.MinFcnValue()/resultSimFit.Ndf(),(float)resultSimFit.MinFcnValue()/resultSimFit.Ndf(),(float)vnSimFit,(float)vnSimFitUnc,(float)qnmin,(float)qnmax,(float)iTrial};
                                        float array4ntupleBinCount[29] = {(float)PtMin[iPt],(float)PtMax[iPt],(float)averagePt,(float)averagePtUnc,(float)kBinCount,(float)hInvMassDeltaPhiToFit[0]->GetBinWidth(1),(float)MassMin[iMassMin],(float)MassMax[iMassMax],(float)SgnFunc,(float)BkgFunc,(float)rawYieldsBinCount[0],(float)rawYieldsBinCount[1],(float)rawYieldsBinCountUnc[0],(float)rawYieldsBinCountUnc[1],-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,(float)vnBinCount,(float)vnBinCountUnc,(float)qnmin,(float)qnmax,(float)iTrial};

                                        multiTrialNtuple->Fill(array4ntupleFreeSigma);
                                        multiTrialNtuple->Fill(array4ntupleFixSigma);
                                        multiTrialNtuple->Fill(array4ntupleSimFit);
                                        multiTrialNtuple->Fill(array4ntupleBinCount);

                                        if(massfitterFreeSigma[0]->GetReducedChiSquare()<maxRedChi2 && massfitterFreeSigma[1]->GetReducedChiSquare()<maxRedChi2) {
                                            hVnResFreeSigma[iPt]->Fill(vnFreeSigma-vnRef);
                                            gVnVsTrialFreeSigma[iPt]->SetPoint(iTrial,iTrial,vnFreeSigma);
                                            gVnVsTrialFreeSigma[iPt]->SetPointError(iTrial,0.,0.,vnFreeSigmaUnc,vnFreeSigmaUnc);
                                        }
                                        if(massfitterFixSigma[0]->GetReducedChiSquare()<maxRedChi2 && massfitterFixSigma[1]->GetReducedChiSquare()<maxRedChi2) {
                                            hVnResFixSigma[iPt]->Fill(vnFixSigma-vnRef);
                                            gVnVsTrialFixSigma[iPt]->SetPoint(iTrial,iTrial,vnFixSigma);
                                            gVnVsTrialFixSigma[iPt]->SetPointError(iTrial,0.,0.,vnFixSigmaUnc,vnFixSigmaUnc);
                                        }
                                        if(massfitterFreeSigma[0]->GetReducedChiSquare()<maxRedChi2 && massfitterFreeSigma[1]->GetReducedChiSquare()<maxRedChi2) {
                                            hVnResBinCounting[iPt]->Fill(vnBinCount-vnRef);
                                            gVnVsTrialBinCounting[iPt]->SetPoint(iTrial,iTrial,vnBinCount);
                                            gVnVsTrialBinCounting[iPt]->SetPointError(iTrial,0.,0.,vnBinCountUnc,vnBinCountUnc);
                                        }
                                        if(resultSimFit.MinFcnValue()/resultSimFit.Ndf()<maxRedChi2) {
                                            hVnResSimFit[iPt]->Fill(vnSimFit-vnRef);
                                            gVnVsTrialSimFit[iPt]->SetPoint(iTrial,iTrial,vnSimFit);
                                            gVnVsTrialSimFit[iPt]->SetPointError(iTrial,0.,0.,vnSimFitUnc,vnSimFitUnc);
                                        }
                                    }

                                    for(int iDeltaPhi=0; iDeltaPhi<2; iDeltaPhi++) {
                                        if(massfitterFreeSigma[iDeltaPhi]) {
                                            delete massfitterFreeSigma[iDeltaPhi];
                                            massfitterFreeSigma[iDeltaPhi] = NULL;
                                        }
                                        if(massfitterFixSigma[iDeltaPhi]) {
                                            delete massfitterFixSigma[iDeltaPhi];
                                            massfitterFixSigma[iDeltaPhi] = NULL;
                                        }
                                        if(massfitterSimFit[iDeltaPhi]) {
                                            delete massfitterSimFit[iDeltaPhi];
                                            massfitterSimFit[iDeltaPhi] = NULL;
                                        }
                                    }
                                }
                                else {
                                    AliHFVnVsMassFitter *vnvsmassfitter = new AliHFVnVsMassFitter(hInvMassIntToFit,hVnVsMassToFit,MassMin[iMassMin],MassMax[iMassMax],BkgFunc,SgnFunc,VnBkgFunc);
                                    vnvsmassfitter->SetHarmonic(harmonic);
                                    vnvsmassfitter->SetInitialGaussianMean(massD,1);
                                    vnvsmassfitter->SetInitialGaussianSigma(massfitterInt.GetSigma(),1);
                                    if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi)
                                        vnvsmassfitter->IncludeSecondGausPeak(massDplus,fixMeanSecP,DplusSigma[iPt],fixSigmaSecP,true);
                                    if(useRefl) {
                                        vnvsmassfitter->SetTemplateReflections(hMCRefl[iPt],reflopt,MassMin[iMassMin],MassMax[iMassMax]);
                                        vnvsmassfitter->SetFixReflOverS(SoverR);
                                        vnvsmassfitter->SetReflVnOption(AliHFVnVsMassFitter::kSameVnSignal);
                                    }
                                    bool isfitok = vnvsmassfitter->SimultaneusFit(false);
                                    if(!isfitok || vnvsmassfitter->GetReducedChiSquare()>10) continue;

                                    float array4ntuple[22] = {(float)PtMin[iPt],(float)PtMax[iPt],(float)averagePt,(float)averagePtUnc,(float)hInvMassIntToFit->GetBinWidth(1),(float)MassMin[iMassMin],(float)MassMax[iMassMax],(float)SgnFunc,(float)BkgFunc,(float)VnBkgFunc,(float)vnvsmassfitter->GetRawYield(),(float)vnvsmassfitter->GetRawYieldUncertainty(),(float)vnvsmassfitter->GetMean(),(float)vnvsmassfitter->GetMeanUncertainty(),(float)vnvsmassfitter->GetSigma(),(float)vnvsmassfitter->GetSigmaUncertainty(),(float)vnvsmassfitter->GetReducedChiSquare(),(float)vnvsmassfitter->GetVn(),(float)vnvsmassfitter->GetVnUncertainty(),(float)qnmin,(float)qnmax,(float)iTrial};

                                    multiTrialNtuple->Fill(array4ntuple);

                                    double vnEst = vnvsmassfitter->GetVn();
                                    if(vnvsmassfitter->GetReducedChiSquare() < maxRedChi2 && TMath::Abs(vnEst - 0.10) > 1e-6) {
                                        hVnResSimFit[iPt]->Fill(vnEst - vnRef);

                                        gVnVsTrialSimFit[iPt]->SetPoint(iTrial,iTrial,vnEst);
                                        gVnVsTrialSimFit[iPt]->SetPointError(iTrial,0.,0.,vnvsmassfitter->GetVnUncertainty(),vnvsmassfitter->GetVnUncertainty());
                                    }
                                }

                                iTrial++;
                            }
                        }

                        if(hVnVsMassToFit) {
                            delete hVnVsMassToFit;
                            hVnVsMassToFit = NULL;
                        }
                        if(hInvMassIntToFit) {
                            delete hInvMassIntToFit;
                            hInvMassIntToFit = NULL;
                        }
                        if(hInvMassDeltaPhiToFit[0]) {
                            delete hInvMassDeltaPhiToFit[0];
                            hInvMassDeltaPhiToFit[0] = NULL;
                        }
                        if(hInvMassDeltaPhiToFit[1]) {
                            delete hInvMassDeltaPhiToFit[1];
                            hInvMassDeltaPhiToFit[1] = NULL;
                        }
                    }
                }
            }
        }

        delete hInvMassInt;
        hInvMassInt = NULL;
        if(hInvMassDeltaPhi[0]) {
            delete hInvMassDeltaPhi[0];
            hInvMassDeltaPhi[0] = NULL;
        }
        if(hInvMassDeltaPhi[1]) {
            delete hInvMassDeltaPhi[1];
            hInvMassDeltaPhi[1] = NULL;
        }
    }

    //plots
    TLegend* legDistr = new TLegend(0.6,0.7,0.8,0.85);
    legDistr->SetTextSize(0.04);
    legDistr->SetFillStyle(0);
    legDistr->AddEntry(hVnResSimFit[0],"Simultaneus fit","l");
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
        legDistr->AddEntry(hVnResFreeSigma[0],"Free sigma","l");
        legDistr->AddEntry(hVnResFixSigma[0],"Fix sigma","l");
        legDistr->AddEntry(hVnResBinCounting[0],"Bin counting","l");
    }

    TLegend* legGraph = new TLegend(0.2,0.15,0.4,0.3);
    legGraph->SetTextSize(0.04);
    legGraph->SetFillStyle(0);
    legGraph->AddEntry(gVnVsTrialSimFit[0],"Simultaneus fit","p");
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
        legGraph->AddEntry(gVnVsTrialFreeSigma[0],"Free sigma","p");
        legGraph->AddEntry(gVnVsTrialFixSigma[0],"Fix sigma","p");
        legGraph->AddEntry(gVnVsTrialBinCounting[0],"Bin counting","p");
    }

    TLegend* legSyst = new TLegend(0.5,0.8,0.8,0.85);
    legSyst->SetTextSize(0.04);
    legSyst->SetFillStyle(0);
    legSyst->AddEntry(lineStatRed[0][0],"Statistical uncertainty","l");

    TLatex* latSimFit[nPtBins];
    TLatex* latFreeSigma[nPtBins];
    TLatex* latFixSigma[nPtBins];
    TLatex* latBinCount[nPtBins];

    TCanvas* cVnVsTrial = new TCanvas("cVnVsTrial","",1920,1080);
    DivideCanvas(cVnVsTrial,nPtBins);

    for(unsigned int iPt=0; iPt<nPtBins; iPt++) {
        cVnVsTrial->cd(iPt+1)->DrawFrame(-0.5,-0.5,iTrial+0.5,0.5,Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c};Trial;#it{v}_{%d}{%s}",PtMin[iPt],PtMax[iPt],harmonic,v2measname.Data()));
        if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
            gVnVsTrialFreeSigma[iPt]->Draw("P");
            gVnVsTrialFixSigma[iPt]->Draw("P");
            gVnVsTrialBinCounting[iPt]->Draw("P");
        }
        gVnVsTrialSimFit[iPt]->Draw("P");
        legGraph->Draw();
    }

    TCanvas* cVnRes = new TCanvas("cVnRes","",1920,1080);
    DivideCanvas(cVnRes,nPtBins);
    for(unsigned int iPt=0; iPt<nPtBins; iPt++) {
        latSimFit[iPt] = BuildTLatex(kBlue,42,0.04);
        latFreeSigma[iPt] = BuildTLatex(kRed,42,0.04);
        latFixSigma[iPt] = BuildTLatex(kBlack,42,0.04);
        latBinCount[iPt] = BuildTLatex(kGreen+2,42,0.04);

        cVnRes->cd(iPt+1)->DrawFrame(-0.5,0.,0.5,iTrial/2,Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c};#it{v}_{%d}-#it{v}_{%d}^{ref}{%s};Entries",PtMin[iPt],PtMax[iPt],harmonic,harmonic,v2measname.Data()));
        hVnResSimFit[iPt]->Draw("same");
        latSimFit[iPt]->DrawLatex(0.2,0.8,Form("RMS = %0.3f",hVnResSimFit[iPt]->GetRMS()));
        latSimFit[iPt]->DrawLatex(0.2,0.75,Form("mean = %0.3f",hVnResSimFit[iPt]->GetMean()));
        if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
            hVnResFreeSigma[iPt]->Draw("same");
            hVnResFixSigma[iPt]->Draw("same");
            hVnResBinCounting[iPt]->Draw("same");
            latFreeSigma[iPt]->DrawLatex(0.2,0.7,Form("RMS = %0.3f",hVnResFreeSigma[iPt]->GetRMS()));
            latFreeSigma[iPt]->DrawLatex(0.2,0.65,Form("mean = %0.3f",hVnResFreeSigma[iPt]->GetMean()));
            latFixSigma[iPt]->DrawLatex(0.2,0.6,Form("RMS = %0.3f",hVnResFixSigma[iPt]->GetRMS()));
            latFixSigma[iPt]->DrawLatex(0.2,0.55,Form("mean = %0.3f",hVnResFixSigma[iPt]->GetMean()));
            latBinCount[iPt]->DrawLatex(0.2,0.5,Form("RMS = %0.3f",hVnResBinCounting[iPt]->GetRMS()));
            latBinCount[iPt]->DrawLatex(0.2,0.45,Form("mean = %0.3f",hVnResBinCounting[iPt]->GetMean()));
        }
        legDistr->Draw();
    }

    TLine* lineatzero = new TLine(0.5,0.,4.5,0.);
    lineatzero->SetLineWidth(2);
    lineatzero->SetLineStyle(9);
    lineatzero->SetLineColor(kBlack);

    TH2F* hFrameSyst[nPtBins];
    TCanvas* cSyst = new TCanvas("cSyst","",1920,1080);
    DivideCanvas(cSyst,nPtBins);
    for(unsigned int iPt=0; iPt<nPtBins; iPt++) {
        double statunc = gVnRef->GetErrorYlow(iPt);
        hFrameSyst[iPt] = new TH2F(Form("hFrameSyst_Pt%d",iPt),Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c};;Mean #pm RMS #it{v}_{%d}-#it{v}_{%d}^{ref}{%s}",PtMin[iPt],PtMax[iPt],harmonic,harmonic,v2measname.Data()),4,0.5,4.5,1,-1.5*statunc,1.5*statunc);
        hFrameSyst[iPt]->SetDirectory(0);
        hFrameSyst[iPt]->GetXaxis()->SetBinLabel(1,"Simultaneus Fit");
        hFrameSyst[iPt]->GetXaxis()->SetBinLabel(2,"Free Sigma");
        hFrameSyst[iPt]->GetXaxis()->SetBinLabel(3,"Fix Sigma");
        hFrameSyst[iPt]->GetXaxis()->SetBinLabel(4,"Bin Counting");
        cSyst->cd(iPt+1);
        hFrameSyst[iPt]->Draw();
        lineatzero->Draw("same");
        lineStatRed[0][iPt]->Draw("same");
        lineStatRed[1][iPt]->Draw("same");
        hVnResMeanVsMethodSimFit[iPt]->SetBinContent(1,hVnResSimFit[iPt]->GetMean());
        hVnResMeanVsMethodSimFit[iPt]->SetBinError(1,hVnResSimFit[iPt]->GetRMS());
        hVnResMeanVsMethodSimFit[iPt]->Draw("same");
        if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
            hVnResMeanVsMethodFreeSigma[iPt]->SetBinContent(2,hVnResFreeSigma[iPt]->GetMean());
            hVnResMeanVsMethodFreeSigma[iPt]->SetBinError(2,hVnResFreeSigma[iPt]->GetRMS());
            hVnResMeanVsMethodFreeSigma[iPt]->Draw("same");
            hVnResMeanVsMethodFixSigma[iPt]->SetBinContent(3,hVnResFixSigma[iPt]->GetMean());
            hVnResMeanVsMethodFixSigma[iPt]->SetBinError(3,hVnResFixSigma[iPt]->GetRMS());
            hVnResMeanVsMethodFixSigma[iPt]->Draw("same");
            hVnResMeanVsMethodBinCounting[iPt]->SetBinContent(4,hVnResBinCounting[iPt]->GetMean());
            hVnResMeanVsMethodBinCounting[iPt]->SetBinError(4,hVnResBinCounting[iPt]->GetRMS());
            hVnResMeanVsMethodBinCounting[iPt]->Draw("same");
        }
        legSyst->Draw();
    }

    //output files
    cSyst->SaveAs(Form("%s/FitSystematics_%sv%d_%s_q%d_%0.f-%0.f.pdf",outputdir.data(),mesonname.data(),harmonic,flowmethodname.data(),harmonic,qnmin,qnmax),"recreate");
    cVnRes->SaveAs(Form("%s/FitSystematics_%sv%d_ResDistr_%s_q%d_%0.f-%0.f.pdf",outputdir.data(),mesonname.data(),harmonic,flowmethodname.data(),harmonic,qnmin,qnmax),"recreate");
    cVnVsTrial->SaveAs(Form("%s/FitSystematics_%sv%d_VsTrial_%s_q%d_%0.f-%0.f.pdf",outputdir.data(),mesonname.data(),harmonic,flowmethodname.data(),harmonic,qnmin,qnmax),"recreate");

    outFile.cd();
    multiTrialNtuple->Write();
    TDirectoryFile dirPlots("Distributions","Distributions");
    dirPlots.cd();
    cVnRes->Write();
    cVnVsTrial->Write();
    cSyst->Write();
    for(unsigned int iPt=0; iPt<nPtBins; iPt++) {
        hVnResSimFit[iPt]->Write();
        gVnVsTrialSimFit[iPt]->Write();
        hVnResMeanVsMethodSimFit[iPt]->Write();
        if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
            hVnResFreeSigma[iPt]->Write();
            gVnVsTrialFreeSigma[iPt]->Write();
            hVnResMeanVsMethodFreeSigma[iPt]->Write();
            hVnResFixSigma[iPt]->Write();
            gVnVsTrialFixSigma[iPt]->Write();
            hVnResMeanVsMethodFixSigma[iPt]->Write();
            hVnResBinCounting[iPt]->Write();
            gVnVsTrialBinCounting[iPt]->Write();
            hVnResMeanVsMethodBinCounting[iPt]->Write();
        }
    }
    outFile.Close();
}

//___________________________________________________________________________________//
//method that returns EP resolution integrated over centrality and the histo of the resolution vs. centrality
TH1D* ComputeEPresolution(double &resol, double &resolunc, TH3F *hDeltaPsiVsqnVsCentr[3], int harmonic, double qnmin, double qnmax) {

    TAxis* centraxis = hDeltaPsiVsqnVsCentr[0]->GetXaxis();
    TAxis* qnaxis = hDeltaPsiVsqnVsCentr[0]->GetYaxis();
    const int ncentrbins = centraxis->GetNbins();
    double centrbinsarray[ncentrbins+1];
    for(int iCentr = 0; iCentr<ncentrbins; iCentr++)
        centrbinsarray[iCentr] = centraxis->GetBinLowEdge(iCentr+1);
    centrbinsarray[ncentrbins] = centraxis->GetBinLowEdge(ncentrbins)+centraxis->GetBinWidth(ncentrbins);
    int qnbinmin = qnaxis->FindBin(qnmin*1.0001);
    int qnbinmax = qnaxis->FindBin(qnmax*0.9999);

    double meanCosnPsi[3], meanUncCosnPsi[3], lowlim[3];
    int nSubEv = 3;
    if(hDeltaPsiVsqnVsCentr[1]->GetEntries()==0) nSubEv = 2;

    TH1D* hResolVsCentr = new TH1D("hResolVsCentr",Form(";centrality (%%);#it{R}_{%d}",harmonic),ncentrbins,centrbinsarray);
    SetHistoStyle(hResolVsCentr,kBlack,kFullCircle);

    for(int iCentr=0; iCentr<ncentrbins+1; iCentr++) {

        int centrbinmin=-1, centrbinmax=-1;
        if(iCentr<ncentrbins) {
            centrbinmin = iCentr+1;
            centrbinmax = iCentr+1;
        }
        else {
            centrbinmin = 0;
            centrbinmax = ncentrbins+1;
        }

        if(nSubEv==3) {
            for(int iResHist=0; iResHist<3; iResHist++) {
                TH1D hDeltaCosnPsi(*(hDeltaPsiVsqnVsCentr[iResHist]->ProjectionZ("hDeltaCosnPsi",centrbinmin,centrbinmax,qnbinmin,qnbinmax)));
                meanCosnPsi[iResHist]    = hDeltaCosnPsi.GetMean();
                meanUncCosnPsi[iResHist] = hDeltaCosnPsi.GetMeanError();
                lowlim[iResHist]         = TMath::Abs(meanCosnPsi[iResHist]-meanUncCosnPsi[iResHist]);
            }
            resol    = TMath::Sqrt(meanCosnPsi[1] * meanCosnPsi[2] / meanCosnPsi[0]);
            resolunc = resol - TMath::Sqrt(lowlim[2] * lowlim[1] / lowlim[0]);
        }
        else {
            TH1F hDeltaCosnPsi(*((TH1F*)hDeltaPsiVsqnVsCentr[0]->ProjectionZ("hDeltaCosnPsi",centrbinmin,centrbinmax,qnbinmin,qnbinmax)));
            resol    = AliVertexingHFUtils::GetFullEvResol(&hDeltaCosnPsi);
            resolunc = TMath::Abs(resol-AliVertexingHFUtils::GetFullEvResolLowLim(&hDeltaCosnPsi));
        }

        if(iCentr<ncentrbins) {
            hResolVsCentr->SetBinContent(iCentr+1,resol);
            hResolVsCentr->SetBinError(iCentr+1,resolunc);
        }
    }

    return hResolVsCentr;
}

//___________________________________________________________________________________//
//method that returns SP resolution integrated over centrality and the histo of the resolution vs. centrality
TH1D* ComputeSPresolution(double &resol, double &resolunc, TH3F *hQiVsqnVsCentr[3], int harmonic, double qnmin, double qnmax) {

    TAxis* centraxis = hQiVsqnVsCentr[0]->GetXaxis();
    TAxis* qnaxis = hQiVsqnVsCentr[0]->GetYaxis();
    const int ncentrbins = centraxis->GetNbins();
    double centrbinsarray[ncentrbins+1];
    for(int iCentr = 0; iCentr<ncentrbins; iCentr++)
        centrbinsarray[iCentr] = centraxis->GetBinLowEdge(iCentr+1);
    centrbinsarray[ncentrbins] = centraxis->GetBinLowEdge(ncentrbins)+centraxis->GetBinWidth(ncentrbins);
    int qnbinmin = qnaxis->FindBin(qnmin*1.0001);
    int qnbinmax = qnaxis->FindBin(qnmax*0.9999);

    double meanScalProd[3], meanUncScalProd[3], lowlim[3];
    int nSubEv = 3;
    if(hQiVsqnVsCentr[1]->GetEntries()==0) nSubEv = 2;

    TH1D* hResolVsCentr = new TH1D("hResolVsCentr",Form(";centrality (%%);#it{R}_{%d}",harmonic),ncentrbins,centrbinsarray);
    SetHistoStyle(hResolVsCentr,kBlack,kFullCircle);

    for(int iCentr=0; iCentr<ncentrbins+1; iCentr++) {

        int centrbinmin=-1, centrbinmax=-1;
        if(iCentr<ncentrbins) {
            centrbinmin = iCentr+1;
            centrbinmax = iCentr+1;
        }
        else {
            centrbinmin = 0;
            centrbinmax = ncentrbins+1;
        }

        if(nSubEv==3) {
            for(int iResHist=0; iResHist<3; iResHist++) {
                TH1D hScalProd(*(hQiVsqnVsCentr[iResHist]->ProjectionZ("hScalProd",centrbinmin,centrbinmax,qnbinmin,qnbinmax)));
                meanScalProd[iResHist]    = hScalProd.GetMean();
                meanUncScalProd[iResHist] = hScalProd.GetMeanError();
                lowlim[iResHist]          = TMath::Abs(meanScalProd[iResHist]-meanUncScalProd[iResHist]);
            }
            resol    = TMath::Sqrt(meanScalProd[1] * meanScalProd[2] / meanScalProd[0]);
            resolunc = resol - TMath::Sqrt(lowlim[2] * lowlim[1] / lowlim[0]);
        }
        else {
            TH1D hScalProd(*((TH1D*)hQiVsqnVsCentr[0]->ProjectionZ("hScalProd",centrbinmin,centrbinmax,qnbinmin,qnbinmax)));
            resol    = TMath::Sqrt(hScalProd.GetMean());
            resolunc = TMath::Sqrt(hScalProd.GetMeanError());
        }

        if(iCentr<ncentrbins) {
            hResolVsCentr->SetBinContent(iCentr+1,resol);
            hResolVsCentr->SetBinError(iCentr+1,resolunc);
        }
    }

    return hResolVsCentr;
}

//___________________________________________________________________________________//
//method that returns in-plane and out-of-plane inv-mass spectra
void GetInOutOfPlaneInvMassHistos(THnSparseF *sparse, TH1F *&hInvMassInPlane, TH1F *&hInvMassOutOfPlane, int harmonic, double qnmin, double qnmax, double ptmin, double ptmax, bool applyML, double MLmin, double MLmax) {

    ApplySelection(sparse,1,ptmin,ptmax);
    ApplySelection(sparse,8,qnmin,qnmax);
    if(applyML)
        ApplySelection(sparse,9,MLmin,MLmax);

    if(harmonic == 2)  {
        //in-plane
        ApplySelection(sparse,2,0.,TMath::Pi()/4);
        hInvMassInPlane = reinterpret_cast<TH1F*>(sparse->Projection(0));
        ResetAxes(sparse,2);
        ApplySelection(sparse,2,3./4*TMath::Pi(),TMath::Pi());
        hInvMassInPlane->Add(reinterpret_cast<TH1F*>(sparse->Projection(0)));
        ResetAxes(sparse,2);
        //out-of-plane
        ApplySelection(sparse,2,TMath::Pi()/4,3./4*TMath::Pi());
        hInvMassOutOfPlane = reinterpret_cast<TH1F*>(sparse->Projection(0));
    }
    else if(harmonic == 3) {
        //in-plane
        ApplySelection(sparse,2,0.,TMath::Pi()/6);
        hInvMassInPlane = reinterpret_cast<TH1F*>(sparse->Projection(0));
        ResetAxes(sparse,2);
        ApplySelection(sparse,2,TMath::Pi()/2,2./3*TMath::Pi());
        hInvMassInPlane->Add(reinterpret_cast<TH1F*>(sparse->Projection(0)));
        ResetAxes(sparse,2);
        //out-of-plane
        ApplySelection(sparse,2,TMath::Pi()/6,TMath::Pi()/2);
        hInvMassOutOfPlane = reinterpret_cast<TH1F*>(sparse->Projection(0));
    }
    ResetAxes(sparse);

    hInvMassInPlane->SetName(Form("hInvMassInPlane_pT_%0.f_%0.f_qn_%0.f_%0.f",ptmin,ptmax,qnmin,qnmax));
    hInvMassOutOfPlane->SetName(Form("hInvMassOutOfPlane_pT_%0.f_%0.f_qn_%0.f_%0.f",ptmin,ptmax,qnmin,qnmax));

    return;
}

//___________________________________________________________________________________//
//method that returns <f(phi)> / R vs. mass histo
TH1F* GetFuncPhiVsMassHistos(THnSparseF* sparse, TString histoname, int iAxis, double qnmin, double qnmax, double ptmin, double ptmax, double massrebin, double resol, bool applyML, double MLmin, double MLmax) {

    ApplySelection(sparse,1,ptmin,ptmax);
    ApplySelection(sparse,8,qnmin,qnmax);
    if(applyML)
        ApplySelection(sparse,9,MLmin,MLmax);

    TH2F* hFuncPhiVsMass2D = reinterpret_cast<TH2F*>(sparse->Projection(0,iAxis));
    TProfile* hProfileFuncPhiVsMass = static_cast<TProfile*>(hFuncPhiVsMass2D->ProfileY("hFuncPhiVsMass2D"));
    //TODO: variable mass bin width
    hProfileFuncPhiVsMass->Rebin(massrebin);
    TH1F* hFuncPhiVsMass = reinterpret_cast<TH1F*>(hProfileFuncPhiVsMass->ProjectionX(histoname.Data()));
    hFuncPhiVsMass->Scale(1./resol);
    hFuncPhiVsMass->GetXaxis()->SetTitleOffset(1.1);

    delete hFuncPhiVsMass2D;
    hFuncPhiVsMass2D = NULL;
    delete hProfileFuncPhiVsMass;
    hProfileFuncPhiVsMass = NULL;

    ResetAxes(sparse);

    return hFuncPhiVsMass;
}

//___________________________________________________________________________________//
//method that returns average pT in [ptmin, ptmax]
float GetAveragePtInRange(float &averagePtUnc, THnSparseF *sparse, double qnmin, double qnmax, double ptmin, double ptmax, int bkgfunc, int sgnfunc, bool useRefl, TH1F* hMCRefl,
                          double SoverR, string reflopt, int meson, double massD, bool fixMeanSecP, bool fixSigmaSecP, float sigmaDplus, bool applyML, double MLmin, double MLmax) {

    ApplySelection(sparse,1,ptmin,ptmax);
    ApplySelection(sparse,8,qnmin,qnmax);
    if(applyML)
        ApplySelection(sparse,9,MLmin,MLmax);

    TH2F* hMassVsPt = reinterpret_cast<TH2F*>(sparse->Projection(0,1));
    TH1F* hMass = reinterpret_cast<TH1F*>(sparse->Projection(0));
    float minfit = hMass->GetBinLowEdge(2);
    float maxfit = hMass->GetBinLowEdge(hMass->GetNbinsX()-1);

    AliHFInvMassFitter massfitter(hMass,minfit,maxfit,bkgfunc,sgnfunc);
    massfitter.SetUseLikelihoodFit();
    massfitter.SetInitialGaussianSigma(0.01);
    massfitter.SetInitialGaussianMean(massD);
    if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi) {
        double massDplus = TDatabasePDG::Instance()->GetParticle(411)->Mass();
        massfitter.IncludeSecondGausPeak(massDplus,fixMeanSecP,sigmaDplus,fixSigmaSecP);
    }
    else if(meson==AliAnalysisTaskSECharmHadronvn::kDstartoKpipi) {
        massfitter.SetInitialGaussianSigma(0.001);
    }
    if(meson==AliAnalysisTaskSECharmHadronvn::kD0toKpi && useRefl) {
        massfitter.SetTemplateReflections(hMCRefl,reflopt,minfit,maxfit);
        massfitter.SetFixReflOverS(SoverR);
    }
    massfitter.MassFitter(false);
    float massFromFit=massfitter.GetMean();
    float sigmaFromFit=massfitter.GetSigma();
    TF1* funcB=massfitter.GetBackgroundRecalcFunc();

    float averagePt = -1;
    AliVertexingHFUtils HFutils;
    HFutils.AveragePt(averagePt,averagePtUnc,ptmin,ptmax,hMassVsPt,massFromFit,sigmaFromFit,funcB,2.5,4.5,0.,3.,1);

    delete hMassVsPt;
    hMassVsPt = NULL;
    delete hMass;
    hMass = NULL;

    ResetAxes(sparse);

    return averagePt;
}

//___________________________________________________________________________________//
//method to compute vn from in-plane and out-of-plane yields
double ComputeEPvn(double &vnunc, int harmonic, double nIn, double nInUnc, double nOut, double nOutUnc, double resol, double corr) {

    double anis = (nIn - nOut) / (nIn + nOut);

    double anisDerivIn  = 2 * nOut / ((nIn + nOut)*(nIn + nOut));
    double anisDerivOut = -2 * nIn / ((nIn + nOut)*(nIn + nOut));

    double anisunc = TMath::Sqrt( anisDerivIn * anisDerivIn * nInUnc * nInUnc + anisDerivOut * anisDerivOut * nOutUnc * nOutUnc + 2 * anisDerivIn * anisDerivOut * nInUnc * nOutUnc * corr);

    double vn = TMath::Pi() / harmonic / harmonic / resol * anis;
    vnunc     = TMath::Pi() / harmonic / harmonic / resol * anisunc;

    return vn;
}

//___________________________________________________________________________________//
//method that applies a selection on sparse axis
void ApplySelection(THnSparseF *sparse, int axisnum, double min, double max) {
    int binmin = sparse->GetAxis(axisnum)->FindBin(min*1.0001);
    int binmax = sparse->GetAxis(axisnum)->FindBin(max*0.9999);
    sparse->GetAxis(axisnum)->SetRange(binmin,binmax);
}

//___________________________________________________________________________________//
//method that resets selections on sparse axes
void ResetAxes(THnSparseF *sparse, int axisnum) {
    if(axisnum >= 0)
        sparse->GetAxis(axisnum)->SetRange(-1,-1);
    else
        for(int iAxis=0; iAxis<sparse->GetNdimensions(); iAxis++) sparse->GetAxis(iAxis)->SetRange(-1,-1);
}

//___________________________________________________________________________________//
//method that returns TList from task output file
TList* LoadTListFromTaskOutput(YAML::Node config) {

    vector<string> filename = config["InputFile"]["FileName"].as<vector<string> >();
    string suffix = config["InputFile"]["Suffix"].as<string>();
    string mesonname = config["InputFile"]["Meson"].as<string>();
    string flowmethodname = config["InputFile"]["FlowMethod"].as<string>();

    TList* list = new TList();
    list->SetOwner();
    TList* listtomerge = new TList();
    for(unsigned int iFile=0; iFile<filename.size(); iFile++) {
        TFile* infile = TFile::Open(filename[iFile].data());
        if(!infile || !infile->IsOpen()) return NULL;
        TDirectoryFile* dir = static_cast<TDirectoryFile*>(infile->Get(Form("PWGHF_D2H_HFvn_%s%s_%s",mesonname.data(),suffix.data(),flowmethodname.data())));
        if(!dir) {
            cerr << Form("TDirectory PWGHF_D2H_HFvn_%s%s_%s not found! Exit",mesonname.data(),suffix.data(),flowmethodname.data()) << endl;
            return NULL;
        }
        TList* listtmp = static_cast<TList*>(dir->Get(Form("coutputvn%s%s_%s",mesonname.data(),suffix.data(),flowmethodname.data())));
        if(!listtmp) {
            cerr << Form("TList coutputvn%s%s_%s not found! Exit",mesonname.data(),suffix.data(),flowmethodname.data()) << endl;
            return NULL;
        }
        if(iFile==0)
            list = listtmp;
        else
            listtomerge->Add(listtmp);
    }
    list->Merge(listtomerge);

    delete listtomerge;
    return list;
}

//__________________________________________________________
//method that loads MC histograms for the reflections of D0
bool LoadD0toKpiReflHistos(string reflFileName, int nPtBins, TH1F* hMCSgn[], TH1F* hMCRefl[]){

    TFile *ReflFile = TFile::Open(reflFileName.data());
    if(!ReflFile){
        cerr << "Error: reflection file "<< reflFileName <<" does not exist! Turning off reflections usage" << endl;
        return false;
    }

    for(int iPt=0; iPt<nPtBins; iPt++){
        hMCSgn[iPt] = NULL;
        hMCSgn[iPt] = static_cast<TH1F*>(ReflFile->Get(Form("histSgn_%d",iPt)));
        if(!hMCSgn[iPt]) {
            cerr << Form("histSgn_%d not found! Turning off reflections usage",iPt) << endl;
            return false;
        }
        hMCSgn[iPt]->SetName(Form("histSgn_%d",iPt));
        hMCSgn[iPt]->SetDirectory(0);
        hMCRefl[iPt] = NULL;
        hMCRefl[iPt] = static_cast<TH1F*>(ReflFile->Get(Form("histRflFittedDoubleGaus_ptBin%d",iPt)));
        if(!hMCRefl[iPt]) {
            cerr << Form("histRflFittedDoubleGaus_ptBin%d not found! Turning off reflections usage",iPt) << endl;
            return false;
        }
        hMCRefl[iPt]->SetName(Form("histRfl_%d",iPt));
        hMCRefl[iPt]->SetDirectory(0);
    }

    ReflFile->Close();

    return true;
}

//__________________________________________________________
//method that loads MC widths for the D+ peak in Ds->KKpi channel
bool LoadDsDplusSigma(string sigmaFileName, int nPtBins, float *sigmaDplus){

    TFile *sigmaFile = TFile::Open(sigmaFileName.data());
    if(!sigmaFile){
        cerr << "Error: sigma file "<< sigmaFileName <<" does not exist!" << endl;
        return false;
    }

    TH1F *hMCSigma = static_cast<TH1F*>(sigmaFile->Get("hRawYieldsSigmaSecondPeak"));
    if(!hMCSigma){
        cerr << "hRawYieldsSigmaSecondPeak not found!" << endl;
        return false;
    }

    for(int iPt=0; iPt<nPtBins; iPt++){
        float sigma = hMCSigma->GetBinContent(iPt+1);
        sigmaDplus[iPt] = sigma;
    }
    sigmaFile->Close();

    return true;
}

//___________________________________________________________________________________//
//method to set plots style
void SetStyle() {
    gROOT->ForceStyle();
    gStyle->SetPadRightMargin(0.035);
    gStyle->SetPadLeftMargin(0.16);
    gStyle->SetPadBottomMargin(0.1);
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetTitleSize(0.045,"xy");
    gStyle->SetLabelSize(0.04,"xy");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetHistLineWidth(2);
    gStyle->SetOptStat(0);
}

//___________________________________________________________________________________//
//method to set plots style
void SetGraphStyle(TGraphAsymmErrors* graph, int color, int markerstyle, float markersize, int linewidth) {
    graph->SetLineColor(color);
    graph->SetMarkerColor(color+1);
    graph->SetMarkerStyle(markerstyle);
    graph->SetMarkerSize(markersize);
    graph->SetLineWidth(linewidth);
}

//___________________________________________________________________________________//
//method to set plots style
void SetHistoStyle(TH1* histo, int color, int markerstyle, float markersize, int linewidth) {
    histo->SetLineColor(color);
    histo->SetMarkerColor(color);
    histo->SetMarkerStyle(markerstyle);
    histo->SetMarkerSize(markersize);
    histo->SetLineWidth(linewidth);
}

//___________________________________________________________________________________//
//method to set plots style
TLatex* BuildTLatex(int color, int font, double fontsize) {
    TLatex* lat = new TLatex();
    lat->SetTextSize(fontsize);
    lat->SetTextColor(color);
    lat->SetTextFont(font);
    lat->SetNDC();
    return lat;
}

//__________________________________________________________________________________________________________________
//helper function to divide canvas
void DivideCanvas(TCanvas* c, const unsigned int nPtBins) {
  if(nPtBins<2)
    return;
  if(nPtBins==2 || nPtBins==3)
    c->Divide(nPtBins,1);
  else if(nPtBins==4 || nPtBins==6 || nPtBins==8)
    c->Divide(nPtBins/2,2);
  else if(nPtBins==5 || nPtBins==7)
    c->Divide((nPtBins+1)/2,2);
  else if(nPtBins==9 || nPtBins==12 || nPtBins==15)
    c->Divide(nPtBins/3,3);
  else if(nPtBins==10 || nPtBins==11)
    c->Divide(4,3);
  else if(nPtBins==13 || nPtBins==14)
    c->Divide(5,3);
  else if(nPtBins>15 && nPtBins<=20 && nPtBins%4==0)
    c->Divide(nPtBins/4,4);
  else if(nPtBins>15 && nPtBins<=20 && nPtBins%4!=0)
    c->Divide(5,4);
  else if(nPtBins==21)
    c->Divide(7,3);
  else if(nPtBins>21 && nPtBins<=25)
    c->Divide(5,5);
  else if(nPtBins>25 && nPtBins%2==0)
    c->Divide(nPtBins/2,2);
  else
    c->Divide((nPtBins+1)/2,2);
}
