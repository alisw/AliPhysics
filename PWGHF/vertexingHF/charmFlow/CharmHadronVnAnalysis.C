//___________________________________________________________________________________//
// Brief: Macro for the analysis of the output of AliAnalysisTaskSECharmHadronvn     //
// Main Function: CharmHadronVnAnalysis                                              //
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

#include "AliHFInvMassFitter.h"
#include "AliHFVnVsMassFitter.h"
#include "AliVertexingHFUtils.h"
#include "AliAnalysisTaskSECharmHadronvn.h"

#endif

using namespace std;

//___________________________________________________________________________________//
//method prototypes
void CharmHadronVnAnalysis(string cfgFileName);
TH1D* ComputeEPresolution(double &resol, double &resolunc, TH3F *hDeltaPsiVsqnVsCentr[3], int harmonic, double qnmin, double qnmax);
TH1D* ComputeSPresolution(double &resol, double &resolunc, TH3F *hQiVsqnVsCentr[3], int harmonic, double qnmin, double qnmax);
void GetInOutOfPlaneInvMassHistos(THnSparseF *sparse, TH1F *&hInvMassInPlane, TH1F *&hInvMassOutOfPlane, int harmonic, double qnmin, double qnmax, double ptmin, double ptmax, bool applyML, double MLmin, double MLmax);
TH1F* GetFuncPhiVsMassHistos(THnSparseF* sparse, TString histoname, int iAxis, double qnmin, double qnmax, double ptmin, double ptmax, double massrebin, bool useVarMassBinning,
                             vector<double> VnVsMassBins, double resol, bool applyML, double MLmin, double MLmax);
float GetAveragePtInRange(float &averagePtUnc, THnSparseF *sparse, double qnmin, double qnmax, double ptmin, double ptmax, int bkgfunc, int sgnfunc, bool useRefl, TH1F* hMCRefl,
                          double SoverR, string reflopt, int meson, double massD, bool fixMeanSecP, bool fixSigmaSecP, float sigmaDplus, bool applyML, double MLmin, double MLmax);
double ComputeEPvn(double &vnunc, int harmonic, double nIn, double nInUnc, double nOut, double nOutUnc, double resol, double corr = 0.);
void ApplySelection(THnSparseF *sparse, int axisnum, double min, double max);
void ResetAxes(THnSparseF *sparse, int axisnum = -1);
TList* LoadTListFromTaskOutput(YAML::Node config);
bool LoadD0toKpiReflHistos(string reflFileName, int nPtBins, TH1F* hMCSgn[], TH1F* hMCRefl[]);
bool LoadDsDplusSigma(string sigmaFileName, int nPtBins, float *sigmaDplus);
double CosnPhi(double *phi, double *pars);
double SinnPhi(double *phi, double *pars);
void SetStyle();
void SetGraphStyle(TGraphAsymmErrors* graph, int color, int markerstyle, float markersize = 1.5, int linewidth = 2);
void SetHistoStyle(TH1* histo, int color, int markerstyle, float markersize = 1.5, int linewidth = 2);

//___________________________________________________________________________________//
//main function for vn analysis
void CharmHadronVnAnalysis(string cfgFileName) {

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
    vector<double> MassMin = config["AnalysisOptions"]["MassMin"].as<vector<double> >();
    vector<double> MassMax = config["AnalysisOptions"]["MassMax"].as<vector<double> >();
    vector<int> Rebin = config["AnalysisOptions"]["Rebin"].as<vector<int> >();
    vector<string> sBkgFunc = config["AnalysisOptions"]["BkgFunc"].as<vector<string> >();
    vector<string> sSgnFunc = config["AnalysisOptions"]["SgnFunc"].as<vector<string> >();
    vector<string> sVnBkgFunc = config["AnalysisOptions"]["VnBkgFunc"].as<vector<string> >();
    int fixMeanVnVsMassFit = config["AnalysisOptions"]["FixMeanVnVsMassFit"].as<int>();
    int fixSigmaVnVsMassFit = config["AnalysisOptions"]["FixSigmaVnVsMassFit"].as<int>();
    bool useVarMassBinning = static_cast<bool>(config["AnalysisOptions"]["UseVarMassBinning"].as<int>());
    vector<double> VnVsMassBins = config["AnalysisOptions"]["VnVsMassBins"].as<vector<double> >();
    bool useRefl = static_cast<bool>(config["AnalysisOptions"]["IncludeReflections"].as<int>());
    string reflFileName = config["AnalysisOptions"]["ReflFileName"].as<string>();
    string reflopt = config["AnalysisOptions"]["ReflOpt"].as<string>();
    bool fixMeanSecP = static_cast<bool>(config["AnalysisOptions"]["FixMeanSecondPeak"].as<int>());
    bool fixSigmaSecP = static_cast<bool>(config["AnalysisOptions"]["FixSigmaSecondPeak"].as<int>());
    string sigmaFileName = config["AnalysisOptions"]["SigmaFileName"].as<string>();

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

    //Check EP residual modulations
    TString detnames[6] = {"TPC","TPCPosEta","TPCNegEta","V0","V0A","V0C"};
    TH3F* hEvPlaneDistrVsqnVsCent[6];
    TH1D* hEvPlaneDistr[6];
    TF1* fCosnPsi[6], *fSinnPsi[6];

    TH1F* hMeanCosnPsi = new TH1F("hMeanCosnPsi",Form(";;<cos(%d#psi_{%d})> / #it{R}_{%d}",harmonic,harmonic,harmonic), 6, 0.5, 6.5);
    TH1F* hMeanSinnPsi = new TH1F("hMeanSinnPsi",Form(";;<sin(%d#psi_{%d})> / #it{R}_{%d}",harmonic,harmonic,harmonic), 6, 0.5, 6.5);
    SetHistoStyle(hMeanCosnPsi,kBlack,kFullSquare);
    SetHistoStyle(hMeanSinnPsi,kBlack,kOpenSquare);

    for(int iDet=0; iDet<6; iDet++) {
        hEvPlaneDistrVsqnVsCent[iDet] = static_cast<TH3F*>(list->FindObject(Form("fHistEvPlaneQncorr%sVsqnVsCent",detnames[iDet].Data())));
        int qnbinmin = hEvPlaneDistrVsqnVsCent[iDet]->GetYaxis()->FindBin(qnmin*1.0001);
        int qnbinmax = hEvPlaneDistrVsqnVsCent[iDet]->GetYaxis()->FindBin(qnmax*0.9999);
        hEvPlaneDistr[iDet] = static_cast<TH1D*>(hEvPlaneDistrVsqnVsCent[iDet]->ProjectionZ(Form("fHistEvPlaneQncorr%s",detnames[iDet].Data()),1,hEvPlaneDistrVsqnVsCent[iDet]->GetXaxis()->GetNbins(),qnbinmin,qnbinmax));
        fCosnPsi[iDet] = new TF1(Form("fCosnPsi%s",detnames[iDet].Data()),CosnPhi,0,TMath::Pi(),3);
        fCosnPsi[iDet]->FixParameter(2,harmonic);
        fCosnPsi[iDet]->SetLineColor(kBlack);
        fSinnPsi[iDet] = new TF1(Form("fSinnPsi%s",detnames[iDet].Data()),SinnPhi,0,TMath::Pi(),3);
        fSinnPsi[iDet]->FixParameter(2,harmonic);
        hEvPlaneDistr[iDet]->Fit(fCosnPsi[iDet],"0");
        hEvPlaneDistr[iDet]->Fit(fSinnPsi[iDet],"0");
        hEvPlaneDistr[iDet]->GetXaxis()->SetTitle(Form("%s %s",detnames[iDet].Data(),hEvPlaneDistr[iDet]->GetXaxis()->GetTitle()));
        hEvPlaneDistr[iDet]->GetYaxis()->SetTitle("Entries");

        hMeanCosnPsi->GetXaxis()->SetBinLabel(iDet+1,detnames[iDet]);
        hMeanSinnPsi->GetXaxis()->SetBinLabel(iDet+1,detnames[iDet]);
        hMeanCosnPsi->SetBinContent(iDet+1,fCosnPsi[iDet]->GetParameter(1));
        hMeanCosnPsi->SetBinError(iDet+1,fCosnPsi[iDet]->GetParError(1));
        hMeanSinnPsi->SetBinContent(iDet+1,fSinnPsi[iDet]->GetParameter(1));
        hMeanSinnPsi->SetBinError(iDet+1,fSinnPsi[iDet]->GetParError(1));
    }

    //Define histos for vn vs. pT
    TH1F*               hInvMassInt[nPtBins];
    TH1F*               hInvMassPhiInt[nPtBins];
    TH1F*               hInvMass[2][nPtBins];
    TH1F*               hVnVsMass[nPtBins];

    AliHFInvMassFitter* massfitterInt[nPtBins];
    AliHFInvMassFitter* massfitterFreeSigma[2][nPtBins];
    AliHFInvMassFitter* massfitterFixSigma[2][nPtBins];
    AliHFInvMassFitter* massfitterSimFit[2][nPtBins];

    AliHFVnVsMassFitter* vnvsmassfitter[nPtBins];

    TGraphAsymmErrors*  gvnFreeSigma = new TGraphAsymmErrors(0);
    TGraphAsymmErrors*  gvnBinCount = new TGraphAsymmErrors(0);
    TGraphAsymmErrors*  gvnFixSigma = new TGraphAsymmErrors(0);
    TGraphAsymmErrors*  gvnSimFit = new TGraphAsymmErrors(0);

    SetGraphStyle(gvnFreeSigma,kRed,kFullCircle);
    SetGraphStyle(gvnBinCount,kGreen+2,kFullTriangleUp);
    SetGraphStyle(gvnFixSigma,kBlack,kFullSquare);
    SetGraphStyle(gvnSimFit,kBlue,kFullDiamond,2.);

    //Phi Modulation histos
    AliHFVnVsMassFitter* sinnphiDvsmassfitter[nPtBins];
    AliHFVnVsMassFitter* cosnphiDvsmassfitter[nPtBins];

    TH1F* hCosnPhiDVsMass[nPtBins];
    TH1F* hSinnPhiDVsMass[nPtBins];

    TH1F* hMeanCosnPhiDVsPt = new TH1F("hMeanCosnPhiDVsPt",Form(";#it{p}_{T} (GeV/#it{c});<cos(%d#varphi_{D})> / #it{R}_{%d}",harmonic,harmonic), nPtBins, PtLims);
    TH1F* hMeanSinnPhiDVsPt = new TH1F("hMeanSinnPhiDVsPt",Form(";#it{p}_{T} (GeV/#it{c});<sin(%d#varphi_{D})> / #it{R}_{%d}",harmonic,harmonic), nPtBins, PtLims);
    SetHistoStyle(hMeanCosnPhiDVsPt,kBlack,kFullCircle);
    SetHistoStyle(hMeanSinnPhiDVsPt,kBlack,kOpenCircle);
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kSP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeSP) {
        hMeanCosnPhiDVsPt->SetName("hMeanUxQyVsPt");
        hMeanCosnPhiDVsPt->SetTitle(Form(";#it{p}_{T} (GeV/#it{c});<u_{D,x}Q_{%d,y}^{A}/M^{A}> / #it{R}_{%d}",harmonic,harmonic));
        hMeanSinnPhiDVsPt->SetName("hMeanUyQxVsPt");
        hMeanSinnPhiDVsPt->SetTitle(Form(";#it{p}_{T} (GeV/#it{c});<u_{D,y}Q_{%d,x}^{A}/M^{A}> / #it{R}_{%d}",harmonic,harmonic));
    }

    //Pars histos
    TH1D* hSigmaInt          = new TH1D("hSigmaInt",";#it{p}_{T} (GeV/#it{c});width (MeV/#it{c}^{2})",nPtBins,PtLims);
    TH1D* hMeanInt           = new TH1D("hMeanInt",";#it{p}_{T} (GeV/#it{c});mean (GeV/#it{c}^{2})",nPtBins,PtLims);
    TH1D* hRedChi2Int        = new TH1D("hRedChi2Int",";#it{p}_{T} (GeV/#it{c});#chi^{2} / ndf",nPtBins,PtLims);
    TH1D* hProbInt           = new TH1D("hProbInt",";#it{p}_{T} (GeV/#it{c});probability",nPtBins,PtLims);
    TH1D* hSigmaSimFit       = new TH1D("hSigmaSimFit",";#it{p}_{T} (GeV/#it{c});width (MeV/#it{c}^{2})",nPtBins,PtLims);
    TH1D* hMeanSimFit        = new TH1D("hMeanSimFit",";#it{p}_{T} (GeV/#it{c});mean (GeV/#it{c}^{2})",nPtBins,PtLims);
    TH1D* hRedChi2SimFit     = new TH1D("hRedChi2SimFit",";#it{p}_{T} (GeV/#it{c});#chi^{2} / ndf",nPtBins,PtLims);
    TH1D* hProbSimFit        = new TH1D("hProbSimFit",";#it{p}_{T} (GeV/#it{c});probability",nPtBins,PtLims);
    TH1D* hRedChi2SBVnPrefit = new TH1D("hRedChi2SBVnPrefit",Form(";#it{p}_{T} (GeV/#it{c});#chi^{2} / ndf SB #it{v}_{%d} prefit",harmonic),nPtBins,PtLims);
    TH1D* hProbSBVnPrefit    = new TH1D("hProbSBVnPrefit",Form(";#it{p}_{T} (GeV/#it{c});probability SB #it{v}_{%d} prefit",harmonic),nPtBins,PtLims);
    SetHistoStyle(hSigmaInt,kGreen+2,kFullTriangleDown);
    SetHistoStyle(hMeanInt,kGreen+2,kFullTriangleDown);
    SetHistoStyle(hRedChi2Int,kGreen+2,kFullTriangleDown);
    SetHistoStyle(hProbInt,kGreen+2,kFullTriangleDown);
    SetHistoStyle(hSigmaSimFit,kBlue,kFullDiamond,2.);
    SetHistoStyle(hMeanSimFit,kBlue,kFullDiamond,2.);
    SetHistoStyle(hRedChi2SimFit,kBlue,kFullDiamond,2.);
    SetHistoStyle(hProbSimFit,kBlue,kFullDiamond,2.);
    SetHistoStyle(hRedChi2SBVnPrefit,kBlue,kFullDiamond,2.);
    SetHistoStyle(hProbSBVnPrefit,kBlue,kFullDiamond,2.);

    TH1D* hSigmaFreeSigma[2];
    TH1D* hMeanFreeSigma[2];
    TH1D* hMeanFixSigma[2];
    TH1D* hRedChi2FreeSigma[2];
    TH1D* hRedChi2FixSigma[2];
    TH1D* hProbFreeSigma[2];
    TH1D* hProbFixSigma[2];
    int markersfreesigma[2] = {kFullCircle,kOpenCircle};
    int markersfixsigma[2] = {kFullSquare,kOpenSquare};
    for(int iDeltaPhi=0; iDeltaPhi<2; iDeltaPhi++) {
        hSigmaFreeSigma[iDeltaPhi]   = new TH1D(Form("hSigmaFreeSigma_phi%d",iDeltaPhi),";#it{p}_{T} (GeV/#it{c});width (MeV/#it{c}^{2})",nPtBins,PtLims);
        hMeanFreeSigma[iDeltaPhi]    = new TH1D(Form("hMeanFreeSigma_phi%d",iDeltaPhi),";#it{p}_{T} (GeV/#it{c});mean (GeV/#it{c}^{2})",nPtBins,PtLims);
        hMeanFixSigma[iDeltaPhi]     = new TH1D(Form("hMeanFixSigma_phi%d",iDeltaPhi),";#it{p}_{T} (GeV/#it{c});mean (GeV/#it{c}^{2})",nPtBins,PtLims);
        hRedChi2FreeSigma[iDeltaPhi] = new TH1D(Form("hRedChi2FreeSigma_phi%d",iDeltaPhi),";#it{p}_{T} (GeV/#it{c});#chi^{2} / ndf",nPtBins,PtLims);
        hRedChi2FixSigma[iDeltaPhi]  = new TH1D(Form("hRedChi2FixSigma_phi%d",iDeltaPhi),";#it{p}_{T} (GeV/#it{c});#chi^{2} / ndf",nPtBins,PtLims);
        hProbFreeSigma[iDeltaPhi]    = new TH1D(Form("hProbFreeSigma_phi%d",iDeltaPhi),";#it{p}_{T} (GeV/#it{c});probability",nPtBins,PtLims);
        hProbFixSigma[iDeltaPhi]     = new TH1D(Form("hProbFixSigma_phi%d",iDeltaPhi),";#it{p}_{T} (GeV/#it{c});probability",nPtBins,PtLims);
        SetHistoStyle(hSigmaFreeSigma[iDeltaPhi],kRed,markersfreesigma[iDeltaPhi]);
        SetHistoStyle(hMeanFreeSigma[iDeltaPhi],kRed,markersfreesigma[iDeltaPhi]);
        SetHistoStyle(hMeanFixSigma[iDeltaPhi],kBlack,markersfixsigma[iDeltaPhi]);
        SetHistoStyle(hRedChi2FreeSigma[iDeltaPhi],kRed,markersfreesigma[iDeltaPhi],2.);
        SetHistoStyle(hRedChi2FixSigma[iDeltaPhi],kBlack,markersfixsigma[iDeltaPhi],2.);
        SetHistoStyle(hProbFreeSigma[iDeltaPhi],kRed,markersfreesigma[iDeltaPhi],2.);
        SetHistoStyle(hProbFixSigma[iDeltaPhi],kBlack,markersfixsigma[iDeltaPhi],2.);
    }

    double rawYieldsFreeSigma[2], rawYieldsBinCount[2], rawYieldsFixSigma[2], rawYieldsSimFit[2], rawYieldsFreeSigmaUnc[2], rawYieldsBinCountUnc[2], rawYieldsFixSigmaUnc[2], rawYieldsSimFitUnc[2];
    for(unsigned int iPt = 0; iPt<nPtBins; iPt++) {

        int SgnFunc = -1, BkgFunc = -1, VnBkgFunc = -1;

        if(sSgnFunc[iPt]=="kGaus")
            SgnFunc = AliHFInvMassFitter::kGaus;
        else if(sSgnFunc[iPt]=="k2Gaus")
            SgnFunc = AliHFInvMassFitter::k2Gaus;
        else if(sSgnFunc[iPt]=="k2GausSigmaRatioPar")
            SgnFunc = AliHFInvMassFitter::k2GausSigmaRatioPar;

        if(sBkgFunc[iPt]=="kExpo")
            BkgFunc = AliHFInvMassFitter::kExpo;
        else if(sBkgFunc[iPt]=="kLin")
            BkgFunc = AliHFInvMassFitter::kLin;
        else if(sBkgFunc[iPt]=="kPol2")
            BkgFunc = AliHFInvMassFitter::kPol2;
        else if(sBkgFunc[iPt]=="kNoBk")
            BkgFunc = AliHFInvMassFitter::kNoBk;
        else if(sBkgFunc[iPt]=="kPow")
            BkgFunc = AliHFInvMassFitter::kPow;
        else if(sBkgFunc[iPt]=="kPowEx")
            BkgFunc = AliHFInvMassFitter::kPowEx;

        if(sVnBkgFunc[iPt]=="kLin")
            VnBkgFunc = AliHFVnVsMassFitter::kLin;
        else if(sVnBkgFunc[iPt]=="kPol2")
            VnBkgFunc = AliHFVnVsMassFitter::kPol2;
        else if(sVnBkgFunc[iPt]=="kExpo")
            VnBkgFunc = AliHFVnVsMassFitter::kExpo;

        double SoverR = 0.;
        if(useRefl)
	        SoverR=(hMCRefl[iPt]->Integral(hMCRefl[iPt]->FindBin(MassMin[iPt]*1.0001),hMCRefl[iPt]->FindBin(MassMax[iPt]*0.9999)))/(hMCSgn[iPt]->Integral(hMCSgn[iPt]->FindBin(MassMin[iPt]*1.0001),hMCSgn[iPt]->FindBin(MassMax[iPt]*0.9999)));

        //compute average pT
        float averagePtUnc = -1.;
        float averagePt = GetAveragePtInRange(averagePtUnc, sMassVsPtVsPhiVsCentrVsqn, qnmin, qnmax, PtMin[iPt], PtMax[iPt], BkgFunc, SgnFunc, 
                                              useRefl, hMCRefl[iPt], SoverR, reflopt, meson, massD, fixMeanSecP, fixSigmaSecP, DplusSigma[iPt], 
                                              doMLsel, CutValuesMLmin[iPt], CutValuesMLmax[iPt]);

        //get histos from sparse
        ApplySelection(sMassVsPtVsPhiVsCentrVsqn,1,PtMin[iPt],PtMax[iPt]);
        if(doMLsel)
            ApplySelection(sMassVsPtVsPhiVsCentrVsqn,9,CutValuesMLmin[iPt], CutValuesMLmax[iPt]);
        hInvMassInt[iPt] = reinterpret_cast<TH1F*>(sMassVsPtVsPhiVsCentrVsqn->Projection(0));

        ApplySelection(sMassVsPtVsPhiVsCentrVsqn,8,qnmin,qnmax);
        hInvMassPhiInt[iPt] = reinterpret_cast<TH1F*>(sMassVsPtVsPhiVsCentrVsqn->Projection(0));

        ResetAxes(sMassVsPtVsPhiVsCentrVsqn,1);
        ResetAxes(sMassVsPtVsPhiVsCentrVsqn,8);
        if(doMLsel)
            ResetAxes(sMassVsPtVsPhiVsCentrVsqn,9);

        if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
            GetInOutOfPlaneInvMassHistos(sMassVsPtVsPhiVsCentrVsqn, hInvMass[0][iPt], hInvMass[1][iPt], harmonic, qnmin, qnmax, PtMin[iPt], PtMax[iPt],
                                         doMLsel, CutValuesMLmin[iPt], CutValuesMLmax[iPt]);
        }
        else if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEPVsMass || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEPVsMass || flowmethod==AliAnalysisTaskSECharmHadronvn::kSP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeSP) {
            hVnVsMass[iPt] = GetFuncPhiVsMassHistos(sMassVsPtVsPhiVsCentrVsqn, Form("hVnVsMass_%d",iPt), 2, qnmin, qnmax, PtMin[iPt], PtMax[iPt], Rebin[iPt], 
                                                    useVarMassBinning, VnVsMassBins, resol, doMLsel, CutValuesMLmin[iPt], CutValuesMLmax[iPt]);
            if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEPVsMass || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEPVsMass)
                hVnVsMass[iPt]->SetTitle(Form(";%s;<cos(%d#Delta#varphi)> / #it{R}_{%d}",hInvMassInt[iPt]->GetXaxis()->GetTitle(),harmonic,harmonic));
            else if(flowmethod==AliAnalysisTaskSECharmHadronvn::kSP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeSP)
                hVnVsMass[iPt]->SetTitle(Form(";%s;<u_{D}Q_{%d,A}/M_{A}> / #it{R}_{%d}",hInvMassInt[iPt]->GetXaxis()->GetTitle(),harmonic,harmonic));
        }

        if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEPVsMass || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEPVsMass || flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
            hCosnPhiDVsMass[iPt] = GetFuncPhiVsMassHistos(sMassVsPtVsPhiVsCentrVsqn, Form("hCosnPhiDVsMass_%d",iPt), 3, qnmin, qnmax, PtMin[iPt], PtMax[iPt], Rebin[iPt],
                                                          useVarMassBinning, VnVsMassBins, resol, doMLsel, CutValuesMLmin[iPt], CutValuesMLmax[iPt]);
            hSinnPhiDVsMass[iPt] = GetFuncPhiVsMassHistos(sMassVsPtVsPhiVsCentrVsqn, Form("hSinnPhiDVsMass_%d",iPt), 4, qnmin, qnmax, PtMin[iPt], PtMax[iPt], Rebin[iPt], 
                                                          useVarMassBinning, VnVsMassBins, resol, doMLsel, CutValuesMLmin[iPt], CutValuesMLmax[iPt]);
            hCosnPhiDVsMass[iPt]->SetTitle(Form(";%s;<cos(%d#varphi_{D})> / #it{R}_{%d}",hInvMassInt[iPt]->GetXaxis()->GetTitle(),harmonic,harmonic));
            hSinnPhiDVsMass[iPt]->SetTitle(Form(";%s;<sin(%d#varphi_{D})> / #it{R}_{%d}",hInvMassInt[iPt]->GetXaxis()->GetTitle(),harmonic,harmonic));
        }
        else {
            hCosnPhiDVsMass[iPt] = GetFuncPhiVsMassHistos(sMassVsPtVsPhiVsCentrVsqn, Form("hUxQyVsMass_%d",iPt), 3, qnmin, qnmax, PtMin[iPt], PtMax[iPt], Rebin[iPt], 
                                                                                          useVarMassBinning, VnVsMassBins, resol, doMLsel, CutValuesMLmin[iPt], CutValuesMLmax[iPt]);
            hSinnPhiDVsMass[iPt] = GetFuncPhiVsMassHistos(sMassVsPtVsPhiVsCentrVsqn, Form("hUyQxVsMass_%d",iPt), 4, qnmin, qnmax, PtMin[iPt], PtMax[iPt], Rebin[iPt], 
                                                                                          useVarMassBinning, VnVsMassBins, resol, doMLsel, CutValuesMLmin[iPt], CutValuesMLmax[iPt]);
            hCosnPhiDVsMass[iPt]->SetTitle(Form(";%s;<u_{D,x}Q_{%dy,A}/M_{A}> / #it{R}_{%d}",hInvMassInt[iPt]->GetXaxis()->GetTitle(),harmonic,harmonic));
            hSinnPhiDVsMass[iPt]->SetTitle(Form(";%s;<u_{D,y}Q_{%dx,A}/M_{A}> / #it{R}_{%d}",hInvMassInt[iPt]->GetXaxis()->GetTitle(),harmonic,harmonic));
        }
        //phi and qn integrated fit
        hInvMassInt[iPt]->Rebin(Rebin[iPt]);
        hInvMassInt[iPt]->GetYaxis()->SetTitle(Form("Counts per %0.f MeV/#it{c}^{2}",hInvMassInt[iPt]->GetBinWidth(1)*1000));
        hInvMassPhiInt[iPt]->Rebin(Rebin[iPt]);
        hInvMassPhiInt[iPt]->GetYaxis()->SetTitle(Form("Counts per %0.f MeV/#it{c}^{2}",hInvMassPhiInt[iPt]->GetBinWidth(1)*1000));
        massfitterInt[iPt] = new AliHFInvMassFitter(hInvMassInt[iPt],MassMin[iPt],MassMax[iPt],BkgFunc,SgnFunc);
        massfitterInt[iPt]->SetUseLikelihoodFit();
        massfitterInt[iPt]->SetInitialGaussianMean(massD);
        massfitterInt[iPt]->SetInitialGaussianSigma(0.010);
        if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi)
            massfitterInt[iPt]->IncludeSecondGausPeak(massDplus,fixMeanSecP,DplusSigma[iPt],fixSigmaSecP);
        else if(meson==AliAnalysisTaskSECharmHadronvn::kDstartoKpipi)
            massfitterInt[iPt]->SetInitialGaussianSigma(0.001);
        if(useRefl) {
            massfitterInt[iPt]->SetTemplateReflections(hMCRefl[iPt],reflopt,MassMin[iPt],MassMax[iPt]);
            massfitterInt[iPt]->SetFixReflOverS(SoverR);
        }
        massfitterInt[iPt]->MassFitter(false);

        hSigmaInt->SetBinContent(iPt+1,massfitterInt[iPt]->GetSigma()*1000);
        hSigmaInt->SetBinError(iPt+1,massfitterInt[iPt]->GetSigmaUncertainty()*1000);
        hMeanInt->SetBinContent(iPt+1,massfitterInt[iPt]->GetMean());
        hMeanInt->SetBinError(iPt+1,massfitterInt[iPt]->GetMeanUncertainty());
        hRedChi2Int->SetBinContent(iPt+1,massfitterInt[iPt]->GetReducedChiSquare());
        hRedChi2Int->SetBinError(iPt+1,1.e-20);
        hProbInt->SetBinContent(iPt+1,massfitterInt[iPt]->GetFitProbability());
        hProbInt->SetBinError(iPt+1,1.e-20);

        ROOT::Fit::FitResult resultSimFit;
        if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
            for(int iDeltaPhi=0; iDeltaPhi<2; iDeltaPhi++) {
                hInvMass[iDeltaPhi][iPt]->Rebin(Rebin[iPt]);
                hInvMass[iDeltaPhi][iPt]->GetYaxis()->SetTitle(Form("Counts per %0.f MeV/#it{c}^{2}",hInvMass[iDeltaPhi][iPt]->GetBinWidth(1)*1000));

                //free sigma
                massfitterFreeSigma[iDeltaPhi][iPt] = new AliHFInvMassFitter(hInvMass[iDeltaPhi][iPt],MassMin[iPt],MassMax[iPt],BkgFunc,SgnFunc);
                massfitterFreeSigma[iDeltaPhi][iPt]->SetUseLikelihoodFit();
                massfitterFreeSigma[iDeltaPhi][iPt]->SetInitialGaussianMean(massD);
                massfitterFreeSigma[iDeltaPhi][iPt]->SetInitialGaussianSigma(0.010);
                if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi)
                    massfitterFreeSigma[iDeltaPhi][iPt]->IncludeSecondGausPeak(massDplus,fixMeanSecP,DplusSigma[iPt],fixSigmaSecP);
                else if(meson==AliAnalysisTaskSECharmHadronvn::kDstartoKpipi)
                    massfitterFreeSigma[iDeltaPhi][iPt]->SetInitialGaussianSigma(0.001);
                if(useRefl) {
                    massfitterFreeSigma[iDeltaPhi][iPt]->SetTemplateReflections(hMCRefl[iPt],reflopt,MassMin[iPt],MassMax[iPt]);
                    massfitterFreeSigma[iDeltaPhi][iPt]->SetFixReflOverS(SoverR);
                }
                massfitterFreeSigma[iDeltaPhi][iPt]->MassFitter(false);

                rawYieldsFreeSigma[iDeltaPhi] = massfitterFreeSigma[iDeltaPhi][iPt]->GetRawYield();
                rawYieldsFreeSigmaUnc[iDeltaPhi] = massfitterFreeSigma[iDeltaPhi][iPt]->GetRawYieldError();

                hSigmaFreeSigma[iDeltaPhi]->SetBinContent(iPt+1,massfitterFreeSigma[iDeltaPhi][iPt]->GetSigma()*1000);
                hSigmaFreeSigma[iDeltaPhi]->SetBinError(iPt+1,massfitterFreeSigma[iDeltaPhi][iPt]->GetSigmaUncertainty()*1000);
                hMeanFreeSigma[iDeltaPhi]->SetBinContent(iPt+1,massfitterFreeSigma[iDeltaPhi][iPt]->GetMean());
                hMeanFreeSigma[iDeltaPhi]->SetBinError(iPt+1,massfitterFreeSigma[iDeltaPhi][iPt]->GetMeanUncertainty());
                hRedChi2FreeSigma[iDeltaPhi]->SetBinContent(iPt+1,massfitterFreeSigma[iDeltaPhi][iPt]->GetReducedChiSquare());
                hRedChi2FreeSigma[iDeltaPhi]->SetBinError(iPt+1,1.e-20);
                hProbFreeSigma[iDeltaPhi]->SetBinContent(iPt+1,massfitterFreeSigma[iDeltaPhi][iPt]->GetFitProbability());
                hProbFreeSigma[iDeltaPhi]->SetBinError(iPt+1,1.e-20);

                //bin counting
                rawYieldsBinCount[iDeltaPhi] = massfitterFreeSigma[iDeltaPhi][iPt]->GetRawYieldBinCounting(rawYieldsBinCountUnc[iDeltaPhi],3.5,0,pdgcode);

                //fix sigma
                massfitterFixSigma[iDeltaPhi][iPt] = new AliHFInvMassFitter(hInvMass[iDeltaPhi][iPt],MassMin[iPt],MassMax[iPt],BkgFunc,SgnFunc);
                massfitterFixSigma[iDeltaPhi][iPt]->SetUseLikelihoodFit();
                massfitterFixSigma[iDeltaPhi][iPt]->SetInitialGaussianMean(massD);
                massfitterFixSigma[iDeltaPhi][iPt]->SetFixGaussianSigma(hSigmaInt->GetBinContent(iPt+1)/1000);
                if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi)
                    massfitterFixSigma[iDeltaPhi][iPt]->IncludeSecondGausPeak(massDplus,fixMeanSecP,DplusSigma[iPt],fixSigmaSecP);
                if(useRefl) {
                    massfitterFixSigma[iDeltaPhi][iPt]->SetTemplateReflections(hMCRefl[iPt],reflopt,MassMin[iPt],MassMax[iPt]);
                    massfitterFixSigma[iDeltaPhi][iPt]->SetFixReflOverS(SoverR);
                }
                massfitterFixSigma[iDeltaPhi][iPt]->MassFitter(false);

                rawYieldsFixSigma[iDeltaPhi] = massfitterFixSigma[iDeltaPhi][iPt]->GetRawYield();
                rawYieldsFixSigmaUnc[iDeltaPhi] = massfitterFixSigma[iDeltaPhi][iPt]->GetRawYieldError();

                hMeanFixSigma[iDeltaPhi]->SetBinContent(iPt+1,massfitterFixSigma[iDeltaPhi][iPt]->GetMean());
                hMeanFixSigma[iDeltaPhi]->SetBinError(iPt+1,massfitterFixSigma[iDeltaPhi][iPt]->GetMeanUncertainty());
                hRedChi2FixSigma[iDeltaPhi]->SetBinContent(iPt+1,massfitterFixSigma[iDeltaPhi][iPt]->GetReducedChiSquare());
                hRedChi2FixSigma[iDeltaPhi]->SetBinError(iPt+1,1.e-20);
                hProbFixSigma[iDeltaPhi]->SetBinContent(iPt+1,massfitterFixSigma[iDeltaPhi][iPt]->GetFitProbability());
                hProbFixSigma[iDeltaPhi]->SetBinError(iPt+1,1.e-20);

                //simultaneus fit
                massfitterSimFit[iDeltaPhi][iPt] = new AliHFInvMassFitter(hInvMass[iDeltaPhi][iPt],MassMin[iPt],MassMax[iPt],BkgFunc,SgnFunc);
                if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi)
                    massfitterSimFit[iDeltaPhi][iPt]->IncludeSecondGausPeak(massDplus,fixMeanSecP,DplusSigma[iPt],fixSigmaSecP);
                if(useRefl) {
                    massfitterSimFit[iDeltaPhi][iPt]->SetTemplateReflections(hMCRefl[iPt],reflopt,MassMin[iPt],MassMax[iPt]);
                    massfitterSimFit[iDeltaPhi][iPt]->SetFixReflOverS(SoverR);
                }
            }
            vector<unsigned int> commonpars = {static_cast<unsigned int>(massfitterFreeSigma[0][iPt]->GetBackgroundFullRangeFunc()->GetNpar())+1,
                                               static_cast<unsigned int>(massfitterFreeSigma[0][iPt]->GetBackgroundFullRangeFunc()->GetNpar())+2};

            if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi){
                commonpars.push_back(static_cast<unsigned int>(massfitterFreeSigma[0][iPt]->GetBackgroundFullRangeFunc()->GetNpar())+4);
                commonpars.push_back(static_cast<unsigned int>(massfitterFreeSigma[0][iPt]->GetBackgroundFullRangeFunc()->GetNpar())+5);
            }

            resultSimFit = AliVertexingHFUtils::DoInPlaneOutOfPlaneSimultaneusFit(massfitterSimFit[0][iPt], massfitterSimFit[1][iPt], hInvMass[0][iPt], hInvMass[1][iPt], MassMin[iPt], MassMax[iPt], massD, commonpars);
            for(int iDeltaPhi=0; iDeltaPhi<2; iDeltaPhi++) {
                rawYieldsSimFit[iDeltaPhi] = massfitterSimFit[iDeltaPhi][iPt]->GetSignalFunc()->GetParameter(0) / hInvMass[0][iPt]->GetBinWidth(1);
                rawYieldsSimFitUnc[iDeltaPhi] = massfitterSimFit[iDeltaPhi][iPt]->GetSignalFunc()->GetParError(0) / hInvMass[0][iPt]->GetBinWidth(1);
            }
            commonpars.clear();

            hSigmaSimFit->SetBinContent(iPt+1,massfitterSimFit[0][iPt]->GetSignalFunc()->GetParameter(2)*1000);
            hSigmaSimFit->SetBinError(iPt+1,massfitterSimFit[0][iPt]->GetSignalFunc()->GetParError(2)*1000);
            hMeanSimFit->SetBinContent(iPt+1,massfitterSimFit[0][iPt]->GetSignalFunc()->GetParameter(1));
            hMeanSimFit->SetBinError(iPt+1,massfitterSimFit[0][iPt]->GetSignalFunc()->GetParError(1));
            hRedChi2SimFit->SetBinContent(iPt+1,resultSimFit.MinFcnValue()/resultSimFit.Ndf());
            hRedChi2SimFit->SetBinError(iPt+1,1.e-20);
            hProbSimFit->SetBinContent(iPt+1,resultSimFit.Prob());
            hProbSimFit->SetBinError(iPt+1,1.e-20);
        }
        else {
            vnvsmassfitter[iPt] = new AliHFVnVsMassFitter(hInvMassPhiInt[iPt],hVnVsMass[iPt],MassMin[iPt],MassMax[iPt],BkgFunc,SgnFunc,VnBkgFunc);
            vnvsmassfitter[iPt]->SetHarmonic(harmonic);
            vnvsmassfitter[iPt]->SetInitialGaussianMean(massD,fixMeanVnVsMassFit);
            vnvsmassfitter[iPt]->SetInitialGaussianSigma(hSigmaInt->GetBinContent(iPt+1)/1000,fixSigmaVnVsMassFit);
            if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi)
                vnvsmassfitter[iPt]->IncludeSecondGausPeak(massDplus,fixMeanSecP,DplusSigma[iPt],fixSigmaSecP,true);
            if(useRefl) {
                vnvsmassfitter[iPt]->SetTemplateReflections(hMCRefl[iPt],reflopt,MassMin[iPt],MassMax[iPt]);
                vnvsmassfitter[iPt]->SetFixReflOverS(SoverR);
                vnvsmassfitter[iPt]->SetReflVnOption(AliHFVnVsMassFitter::kSameVnSignal);
            }
            vnvsmassfitter[iPt]->SimultaneusFit(false);

            hSigmaSimFit->SetBinContent(iPt+1,vnvsmassfitter[iPt]->GetSigma()*1000);
            hSigmaSimFit->SetBinError(iPt+1,vnvsmassfitter[iPt]->GetSigmaUncertainty()*1000);
            hMeanSimFit->SetBinContent(iPt+1,vnvsmassfitter[iPt]->GetMean());
            hMeanSimFit->SetBinError(iPt+1,vnvsmassfitter[iPt]->GetMeanUncertainty());
            hRedChi2SimFit->SetBinContent(iPt+1,vnvsmassfitter[iPt]->GetReducedChiSquare());
            hRedChi2SimFit->SetBinError(iPt+1,1.e-20);
            hProbSimFit->SetBinContent(iPt+1,vnvsmassfitter[iPt]->GetFitProbability());
            hProbSimFit->SetBinError(iPt+1,1.e-20);
            hRedChi2SBVnPrefit->SetBinContent(iPt+1,vnvsmassfitter[iPt]->GetSBVnPrefitReducedChiSquare());
            hRedChi2SBVnPrefit->SetBinError(iPt+1,1.e-20);
            hProbSBVnPrefit->SetBinContent(iPt+1,vnvsmassfitter[iPt]->GetSBVnPrefitProbability());
            hProbSBVnPrefit->SetBinError(iPt+1,1.e-20);
        }

        //compute vn
        double vnFreeSigma=0., vnBinCount=0., vnFixSigma=0., vnSimFit=0., vnFreeSigmaUnc=0., vnBinCountUnc=0., vnFixSigmaUnc=0., vnSimFitUnc=0.;
        if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
            vnFreeSigma = ComputeEPvn(vnFreeSigmaUnc,harmonic,rawYieldsFreeSigma[0],rawYieldsFreeSigmaUnc[0],rawYieldsFreeSigma[1],rawYieldsFreeSigmaUnc[1],resol);
            vnBinCount = ComputeEPvn(vnBinCountUnc,harmonic,rawYieldsBinCount[0],rawYieldsBinCountUnc[0],rawYieldsBinCount[1],rawYieldsBinCountUnc[1],resol);
            vnFixSigma = ComputeEPvn(vnFixSigmaUnc,harmonic,rawYieldsFixSigma[0],rawYieldsFixSigmaUnc[0],rawYieldsFixSigma[1],rawYieldsFixSigmaUnc[1],resol);
            int posRawYieldPar = massfitterFreeSigma[0][iPt]->GetBackgroundFullRangeFunc()->GetNpar();
            int nTotPars = massfitterFreeSigma[0][iPt]->GetMassFunc()->GetNpar();
            vnSimFit = ComputeEPvn(vnSimFitUnc,harmonic,rawYieldsSimFit[0],rawYieldsSimFitUnc[0],rawYieldsSimFit[1],rawYieldsSimFitUnc[1],resol,resultSimFit.Correlation(posRawYieldPar,posRawYieldPar+nTotPars));

            gvnFreeSigma->SetPoint(iPt,averagePt,vnFreeSigma);
            gvnFreeSigma->SetPointError(iPt,averagePt-PtMin[iPt],PtMax[iPt]-averagePt,vnFreeSigmaUnc,vnFreeSigmaUnc);
            gvnBinCount->SetPoint(iPt,averagePt,vnBinCount);
            gvnBinCount->SetPointError(iPt,averagePt-PtMin[iPt],PtMax[iPt]-averagePt,vnBinCountUnc,vnBinCountUnc);
            gvnFixSigma->SetPoint(iPt,averagePt,vnFixSigma);
            gvnFixSigma->SetPointError(iPt,averagePt-PtMin[iPt],PtMax[iPt]-averagePt,vnFixSigmaUnc,vnFixSigmaUnc);
        }
        else if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEPVsMass || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEPVsMass || flowmethod==AliAnalysisTaskSECharmHadronvn::kSP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeSP) {
            vnSimFit = vnvsmassfitter[iPt]->GetVn();
            vnSimFitUnc = vnvsmassfitter[iPt]->GetVnUncertainty();
        }

        gvnSimFit->SetPoint(iPt,averagePt,vnSimFit);
        gvnSimFit->SetPointError(iPt,averagePt-PtMin[iPt],PtMax[iPt]-averagePt,vnSimFitUnc,vnSimFitUnc);

        //Check D-meson phi modulations
        cosnphiDvsmassfitter[iPt] = new AliHFVnVsMassFitter(hInvMassPhiInt[iPt],hCosnPhiDVsMass[iPt],MassMin[iPt],MassMax[iPt],BkgFunc,SgnFunc,VnBkgFunc);
        cosnphiDvsmassfitter[iPt]->SetHarmonic(harmonic);
        cosnphiDvsmassfitter[iPt]->SetInitialGaussianMean(massD,fixMeanVnVsMassFit);
        cosnphiDvsmassfitter[iPt]->SetInitialGaussianSigma(hSigmaInt->GetBinContent(iPt+1)/1000,fixSigmaVnVsMassFit);
        if(useRefl) {
            cosnphiDvsmassfitter[iPt]->SetTemplateReflections(hMCRefl[iPt],reflopt,MassMin[iPt],MassMax[iPt]);
            cosnphiDvsmassfitter[iPt]->SetFixReflOverS(SoverR);
            cosnphiDvsmassfitter[iPt]->SetReflVnOption(AliHFVnVsMassFitter::kSameVnSignal);
        }
        if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi)
            cosnphiDvsmassfitter[iPt]->IncludeSecondGausPeak(massDplus,fixMeanSecP,DplusSigma[iPt],fixSigmaSecP,true);
        cosnphiDvsmassfitter[iPt]->SimultaneusFit(false);

        sinnphiDvsmassfitter[iPt] = new AliHFVnVsMassFitter(hInvMassPhiInt[iPt],hSinnPhiDVsMass[iPt],MassMin[iPt],MassMax[iPt],BkgFunc,SgnFunc,VnBkgFunc);
        sinnphiDvsmassfitter[iPt]->SetHarmonic(harmonic);
        sinnphiDvsmassfitter[iPt]->SetInitialGaussianMean(massD,fixMeanVnVsMassFit);
        sinnphiDvsmassfitter[iPt]->SetInitialGaussianSigma(hSigmaInt->GetBinContent(iPt+1)/1000,fixSigmaVnVsMassFit);
        if(useRefl) {
            sinnphiDvsmassfitter[iPt]->SetTemplateReflections(hMCRefl[iPt],reflopt,MassMin[iPt],MassMax[iPt]);
            sinnphiDvsmassfitter[iPt]->SetFixReflOverS(SoverR);
            sinnphiDvsmassfitter[iPt]->SetReflVnOption(AliHFVnVsMassFitter::kSameVnSignal);
        }
        if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi)
            sinnphiDvsmassfitter[iPt]->IncludeSecondGausPeak(massDplus,fixMeanSecP,DplusSigma[iPt],fixSigmaSecP,true);
        sinnphiDvsmassfitter[iPt]->SimultaneusFit(false);

        hMeanCosnPhiDVsPt->SetBinContent(iPt+1,cosnphiDvsmassfitter[iPt]->GetVn());
        hMeanCosnPhiDVsPt->SetBinError(iPt+1,cosnphiDvsmassfitter[iPt]->GetVnUncertainty());
        hMeanSinnPhiDVsPt->SetBinContent(iPt+1,sinnphiDvsmassfitter[iPt]->GetVn());
        hMeanSinnPhiDVsPt->SetBinError(iPt+1,sinnphiDvsmassfitter[iPt]->GetVnUncertainty());
    }

    //plots
    TLegend* legRes = new TLegend(0.18,0.25,0.4,0.4);
    legRes->SetTextSize(0.04);
    legRes->SetFillStyle(0);
    legRes->AddEntry(hResolVsCentr,"vs. centrality","p");
    legRes->AddEntry(hResolCentrInt,Form("centrality integrated, #it{R}_{%d} = %.4f #pm %.4f",harmonic,resol,resolunc),"p");

    TLegend* legVn = new TLegend(0.6,0.7,0.8,0.9);
    legVn->SetTextSize(0.04);
    legVn->SetFillStyle(0);
    legVn->AddEntry(gvnSimFit,"Simultaneus fit","p");
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
        legVn->AddEntry(gvnFreeSigma,"Free sigma","p");
        legVn->AddEntry(gvnFixSigma,"Fix sigma","p");
        legVn->AddEntry(gvnBinCount,"Bin Counting","p");
    }

    TLegend* legSigma = new TLegend(0.3,0.7,0.8,0.92);
    legSigma->SetTextSize(0.04);
    legSigma->SetFillStyle(0);
    legSigma->AddEntry(hSigmaInt,"#varphi-integrated","p");
    legSigma->AddEntry(hSigmaSimFit,"Simultaneus fit","p");
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
        legSigma->AddEntry(hSigmaFreeSigma[0],"Free sigma - in plane","p");
        legSigma->AddEntry(hSigmaFreeSigma[1],"Free sigma - out of plane","p");
    }

    TLegend* legMean = static_cast<TLegend*>(legSigma->Clone());
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
        legMean->AddEntry(hMeanFixSigma[0],"Fix sigma - in plane","p");
        legMean->AddEntry(hMeanFixSigma[1],"Fix sigma - out of plane","p");
    }

    TLegend* legEPMod = new TLegend(0.3,0.7,0.8,0.92);
    legEPMod->SetTextSize(0.04);
    legEPMod->SetFillStyle(0);
    legEPMod->AddEntry(hMeanCosnPsi,Form("#it{f}(%d#psi_{%d}) = cos(%d#psi_{%d})",harmonic,harmonic,harmonic,harmonic),"p");
    legEPMod->AddEntry(hMeanSinnPsi,Form("#it{f}(%d#psi_{%d}) = sin(%d#psi_{%d})",harmonic,harmonic,harmonic,harmonic),"p");

    TLegend* legPhiDMod = new TLegend(0.3,0.7,0.8,0.92);
    legPhiDMod->SetTextSize(0.04);
    legPhiDMod->SetFillStyle(0);
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kSP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeSP) {
        legPhiDMod->AddEntry(hMeanCosnPhiDVsPt,Form("#it{f}(%d#varphi_{D}) = <u_{D,x}Q_{%d,y}^{A}/M^{A}>",harmonic,harmonic),"p");
        legPhiDMod->AddEntry(hMeanSinnPhiDVsPt,Form("#it{f}(%d#varphi_{D}) = <u_{D,y}Q_{%d,x}^{A}/M^{A}>",harmonic,harmonic),"p");
    }
    else {
        legPhiDMod->AddEntry(hMeanCosnPhiDVsPt,Form("#it{f}(%d#varphi_{D}) = cos(%d#varphi_{D})",harmonic,harmonic),"p");
        legPhiDMod->AddEntry(hMeanSinnPhiDVsPt,Form("#it{f}(%d#varphi_{D}) = sin(%d#varphi_{D})",harmonic,harmonic),"p");
    }

    TCanvas* cResol = new TCanvas("cResol","",800,800);
    double maxRes = 1.;
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kSP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeSP)
        maxRes = 0.05;
    cResol->DrawFrame(hResolVsCentr->GetBinLowEdge(1),0.,hResolVsCentr->GetBinLowEdge(hResolVsCentr->GetNbinsX())+hResolVsCentr->GetBinWidth(1),maxRes,Form(";centrality (%%);#it{R}_{%d}",harmonic));
    hResolVsCentr->Draw("same");
    hResolCentrInt->Draw("same");
    legRes->Draw();

    TCanvas* cEvPlaneDist = new TCanvas("cEvPlaneDist","",1920,1080);
    cEvPlaneDist->Divide(3,2);
    for(int iDet=0; iDet<6; iDet++) {
        cEvPlaneDist->cd(iDet+1);
        hEvPlaneDistr[iDet]->GetYaxis()->SetRangeUser(0.,hEvPlaneDistr[iDet]->GetMaximum()*1.2);
        hEvPlaneDistr[iDet]->Draw("E");
        fCosnPsi[iDet]->Draw("same");
        fSinnPsi[iDet]->Draw("same");
    }

    TCanvas* cEvPlaneModulations = new TCanvas("cEvPlaneModulations","",1920,1080);
    TH2F* hFrameModul = new TH2F("hFrameModul",Form(";;<#it{f}(%d#psi_{%d})>",harmonic,harmonic),6,0.5,6.5,100,-0.1,0.1);
    for(int iDet=0; iDet<6; iDet++)
        hFrameModul->GetXaxis()->SetBinLabel(iDet+1,detnames[iDet]);
    hFrameModul->Draw();
    hMeanCosnPsi->Draw("same");
    hMeanSinnPsi->Draw("same");
    legEPMod->Draw();

    TCanvas *cInvMassFreeSigma = NULL, *cInvMassFixSigma = NULL, *cInvMassSimFit = NULL;
    TCanvas *cVnVsMassSimFit[nPtBins];
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
        cInvMassFreeSigma = new TCanvas("cInvMassFreeSigma","free sigma",1920,1080);
        cInvMassFreeSigma->Divide(nPtBins,2);
        for(unsigned int iPt=0; iPt<nPtBins; iPt++) {
            cInvMassFreeSigma->cd(iPt+1);
            massfitterFreeSigma[0][iPt]->DrawHere(gPad,3,2);
            cInvMassFreeSigma->cd(nPtBins+iPt+1);
            massfitterFreeSigma[1][iPt]->DrawHere(gPad,3,2);
        }

        cInvMassFixSigma = new TCanvas("cInvMassFixSigma","fix sigma",1920,1080);
        cInvMassFixSigma->Divide(nPtBins,2);
        for(unsigned int iPt=0; iPt<nPtBins; iPt++) {
            cInvMassFixSigma->cd(iPt+1);
            massfitterFixSigma[0][iPt]->DrawHere(gPad,3,2);
            cInvMassFixSigma->cd(nPtBins+iPt+1);
            massfitterFixSigma[1][iPt]->DrawHere(gPad,3,2);
        }

        cInvMassSimFit = new TCanvas("cInvMassSimFit","simultaneus fit",1920,1080);
        cInvMassSimFit->Divide(nPtBins,2);
        for(unsigned int iPt=0; iPt<nPtBins; iPt++) {
            cInvMassSimFit->cd(iPt+1);
            massfitterSimFit[0][iPt]->DrawHere(gPad,3,1);
            cInvMassSimFit->cd(nPtBins+iPt+1);
            massfitterSimFit[1][iPt]->DrawHere(gPad,3,1);
        }
    }
    else if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEPVsMass || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEPVsMass || flowmethod==AliAnalysisTaskSECharmHadronvn::kSP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeSP) {
        for(unsigned int iPt=0; iPt<nPtBins; iPt++) {
            cVnVsMassSimFit[iPt] = new TCanvas(Form("cVnVsMassSimFit_pt%d",iPt),"",600,800);
            cVnVsMassSimFit[iPt]->cd();
            vnvsmassfitter[iPt]->DrawHere(gPad);
        }
    }

    TCanvas* cSigma = new TCanvas("cSigma","",800,800);
    cSigma->DrawFrame(PtMin[0],0.,PtMax[nPtBins-1],40.,";#it{p}_{T} (GeV/#it{c});width (MeV/#it{c}^{2})");
    hSigmaInt->Draw("same");
    hSigmaSimFit->Draw("same");
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
        for(int iDeltaPhi=0; iDeltaPhi<2; iDeltaPhi++)
            hSigmaFreeSigma[iDeltaPhi]->Draw("same");
    }
    legSigma->Draw();

    TCanvas* cMean = new TCanvas("cMean","",800,800);
    cMean->DrawFrame(PtMin[0],massD*0.99,PtMax[nPtBins-1],massD*1.01,";#it{p}_{T} (GeV/#it{c});mean (GeV/#it{c}^{2})");
    hMeanInt->Draw("same");
    hMeanSimFit->Draw("same");
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
        for(int iDeltaPhi=0; iDeltaPhi<2; iDeltaPhi++) {
            hMeanFreeSigma[iDeltaPhi]->Draw("same");
            hMeanFixSigma[iDeltaPhi]->Draw("same");
        }
    }
    legMean->Draw();

    TCanvas* cChiSquare = new TCanvas("cChiSquare","",800,800);
    cChiSquare->DrawFrame(PtMin[0],0.,PtMax[nPtBins-1],3.,";#it{p}_{T} (GeV/#it{c});#chi^{2} / ndf");
    hRedChi2Int->Draw("same");
    hRedChi2SimFit->Draw("same");
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
        for(int iDeltaPhi=0; iDeltaPhi<2; iDeltaPhi++) {
            hRedChi2FreeSigma[iDeltaPhi]->Draw("same");
            hRedChi2FixSigma[iDeltaPhi]->Draw("same");
        }
    }
    legMean->Draw();

    TCanvas* cProbability = new TCanvas("cProbability","",800,800);
    cProbability->DrawFrame(PtMin[0],0.,PtMax[nPtBins-1],3.,";#it{p}_{T} (GeV/#it{c});probability");
    hProbInt->Draw("same");
    hProbSimFit->Draw("same");
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
        for(int iDeltaPhi=0; iDeltaPhi<2; iDeltaPhi++) {
            hProbFreeSigma[iDeltaPhi]->Draw("same");
            hProbFixSigma[iDeltaPhi]->Draw("same");
        }
    }
    legMean->Draw();

    TCanvas* cSBVnPrefitChiSquare = NULL;
    TCanvas* cSBVnPrefitProb = NULL;
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kSP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeSP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEPVsMass || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEPVsMass) {
        cSBVnPrefitChiSquare = new TCanvas("cSBVnPrefitChiSquare","",800,800);
        cSBVnPrefitChiSquare->DrawFrame(PtMin[0],0.,PtMax[nPtBins-1],3.,Form(";#it{p}_{T} (GeV/#it{c});#chi^{2} / ndf SB #it{v}_{%d} prefit",harmonic));
        hRedChi2SBVnPrefit->Draw("same");
        cSBVnPrefitProb = new TCanvas("cSBVnPrefitProb","",800,800);
        cSBVnPrefitProb->DrawFrame(PtMin[0],0.,PtMax[nPtBins-1],3.,Form(";#it{p}_{T} (GeV/#it{c});probability SB #it{v}_{%d} prefit",harmonic));
        hProbSBVnPrefit->Draw("same");
    }

    TCanvas *cCosnphiD[nPtBins], *cSinnphiD[nPtBins];
    TCanvas *cPhiDModulationsVsPt = NULL;
    for(unsigned int iPt=0; iPt<nPtBins; iPt++) {
        TString canvasnames[2] = {Form("cCosnphiD_pt%d",iPt),Form("cSinnphiD_pt%d",iPt)};
        if(flowmethod==AliAnalysisTaskSECharmHadronvn::kSP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeSP) {
            canvasnames[0] = Form("cUxQy_pt%d",iPt);
            canvasnames[1] = Form("cUyQx_pt%d",iPt);
        }
        cCosnphiD[iPt] = new TCanvas(canvasnames[0],"",600,800);
        cCosnphiD[iPt]->cd();
        cosnphiDvsmassfitter[iPt]->DrawHere(gPad);
        cSinnphiD[iPt] = new TCanvas(canvasnames[1],"",600,800);
        cSinnphiD[iPt]->cd();
        sinnphiDvsmassfitter[iPt]->DrawHere(gPad);
    }
    cPhiDModulationsVsPt = new TCanvas("cPhiDModulationsVsPt","",800,800);
    cPhiDModulationsVsPt->DrawFrame(PtMin[0],-0.5,PtMax[nPtBins-1],0.5,Form(";#it{p}_{T} (GeV/#it{c});<#it{f}(%d#varphi)> / #it{R}_{%d}",harmonic,harmonic));
    hMeanCosnPhiDVsPt->Draw("same");
    hMeanSinnPhiDVsPt->Draw("same");
    legPhiDMod->Draw();

    TCanvas* cvn = new TCanvas("cvn","",800,800);
    TString v2name = "EP";
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kSP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeSP)
        v2name = "SP";
    cvn->DrawFrame(PtMin[0]-0.5,-0.2,PtMax[nPtBins-1]+0.5,0.4,Form(";#it{p}_{T} (GeV/#it{c});#it{v}_{%d} {%s}",harmonic,v2name.Data()));
    gvnSimFit->Draw("P");
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
        gvnFreeSigma->Draw("P");
        gvnBinCount->Draw("P");
        gvnFixSigma->Draw("P");
    }
    TLine *line = new TLine(PtMin[0], 0., PtMax[nPtBins-1], 0.);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("SAME");
    legVn->Draw();

    //output files
    string outputdir = config["OutputDir"]["Analysis"].as<string>();
    TString resoltypename = "EP";
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kSP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeSP)
        resoltypename = "SP";
    cResol->SaveAs(Form("%s/%sresol_v%d_%s_q%d_%0.f-%0.f.pdf",outputdir.data(),resoltypename.Data(),harmonic,flowmethodname.data(),harmonic,qnmin,qnmax));
    cvn->SaveAs(Form("%s/%sv%d_%s_q%d_%0.f-%0.f.pdf",outputdir.data(),mesonname.data(),harmonic,flowmethodname.data(),harmonic,qnmin,qnmax));

    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
        cInvMassFreeSigma->SaveAs(Form("%s/InvMassFitsFreeSigma%s_%s_q%d_%0.f-%0.f.pdf",outputdir.data(),mesonname.data(),flowmethodname.data(),harmonic,qnmin,qnmax));
        cInvMassFixSigma->SaveAs(Form("%s/InvMassFitsFixSigma%s_%s_q%d_%0.f-%0.f.pdf",outputdir.data(),mesonname.data(),flowmethodname.data(),harmonic,qnmin,qnmax));
        cInvMassSimFit->SaveAs(Form("%s/InvMassSimulFits%s_%s_q%d_%0.f-%0.f.pdf",outputdir.data(),mesonname.data(),flowmethodname.data(),harmonic,qnmin,qnmax));
    }
    else {
        cVnVsMassSimFit[0]->Print(Form("%s/V%dVsMassSimFit%s_%s_q%d_%0.f-%0.f.pdf[",outputdir.data(),harmonic,mesonname.data(),flowmethodname.data(),harmonic,qnmin,qnmax));
        for(unsigned int iPt=0; iPt<nPtBins; iPt++) {
            cVnVsMassSimFit[iPt]->Print(Form("%s/V%dVsMassSimFit%s_%s_q%d_%0.f-%0.f.pdf",outputdir.data(),harmonic,mesonname.data(),flowmethodname.data(),harmonic,qnmin,qnmax));
        }
        cVnVsMassSimFit[nPtBins-1]->Print(Form("%s/V%dVsMassSimFit%s_%s_q%d_%0.f-%0.f.pdf]",outputdir.data(),harmonic,mesonname.data(),flowmethodname.data(),harmonic,qnmin,qnmax));
    }

    TFile outFile(Form("%s/%sv%d_%s_q%d_%0.f-%0.f.root",outputdir.data(),mesonname.data(),harmonic,flowmethodname.data(),harmonic,qnmin,qnmax),"recreate");
    hResolVsCentr->Write();
    hResolCentrInt->Write();
    gvnSimFit->Write("gvnSimFit");
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
        gvnFreeSigma->Write("gvnFreeSigma");
        gvnBinCount->Write("gvnBinCount");
        gvnFixSigma->Write("gvnFixSigma");
    }
    TDirectoryFile dirFits("Fits","Fits");
    dirFits.Write();
    dirFits.cd();
    hSigmaInt->Write();
    hMeanInt->Write();
    hRedChi2Int->Write();
    hProbInt->Write();
    hSigmaSimFit->Write();
    hMeanSimFit->Write();
    hRedChi2SimFit->Write();
    hProbSimFit->Write();
    if(flowmethod==AliAnalysisTaskSECharmHadronvn::kEP || flowmethod==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
        cInvMassFreeSigma->Write();
        cInvMassFixSigma->Write();
        cInvMassSimFit->Write();
        for(int iDeltaPhi=0; iDeltaPhi<2; iDeltaPhi++) {
            hSigmaFreeSigma[iDeltaPhi]->Write();
            hMeanFreeSigma[iDeltaPhi]->Write();
            hRedChi2FreeSigma[iDeltaPhi]->Write();
            hMeanFixSigma[iDeltaPhi]->Write();
            hRedChi2FixSigma[iDeltaPhi]->Write();
        }
    }
    else {
        hRedChi2SBVnPrefit->Write();
        hProbSBVnPrefit->Write();
        for(unsigned int iPt=0; iPt<nPtBins; iPt++) {
            cVnVsMassSimFit[iPt]->Write();
        }
    }
    outFile.cd();
    TDirectoryFile dirEPMod("EPModulations","EPModulations");
    dirEPMod.Write();
    dirEPMod.cd();
    for(int iDet=0; iDet<6; iDet++) {
        hEvPlaneDistr[iDet]->Write();
        fCosnPsi[iDet]->Write();
        fSinnPsi[iDet]->Write();
    }
    hMeanCosnPsi->Write();
    hMeanSinnPsi->Write();
    outFile.cd();
    TDirectoryFile dirPhiMod("PhiDModulations","PhiDModulations");
    dirPhiMod.Write();
    dirPhiMod.cd();
    for(unsigned int iPt=0; iPt<nPtBins; iPt++) {
        cCosnphiD[iPt]->Write();
        cSinnphiD[iPt]->Write();
    }
    hMeanSinnPhiDVsPt->Write();
    hMeanCosnPhiDVsPt->Write();
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
TH1F* GetFuncPhiVsMassHistos(THnSparseF* sparse, TString histoname, int iAxis, double qnmin, double qnmax, double ptmin, double ptmax, double massrebin, bool useVarMassBinning, vector<double> VnVsMassBins, double resol, bool applyML, double MLmin, double MLmax) {

    ApplySelection(sparse,1,ptmin,ptmax);
    ApplySelection(sparse,8,qnmin,qnmax);
    if(applyML)
        ApplySelection(sparse,9,MLmin,MLmax);

    TH2F* hFuncPhiVsMass2D = reinterpret_cast<TH2F*>(sparse->Projection(0,iAxis));
    TProfile* hProfileFuncPhiVsMass = static_cast<TProfile*>(hFuncPhiVsMass2D->ProfileY("hProfileFuncPhiVsMass"));
    TProfile* hProfileFuncPhiVsMassReb = NULL;
    if(!useVarMassBinning) {
        hProfileFuncPhiVsMassReb = static_cast<TProfile*>(hProfileFuncPhiVsMass->Clone("hProfileFuncPhiVsMassReb"));
        hProfileFuncPhiVsMassReb->Rebin(massrebin);
    }
    else {
        hProfileFuncPhiVsMassReb = static_cast<TProfile*>(hProfileFuncPhiVsMass->Rebin(VnVsMassBins.size()-1, "hProfileFuncPhiVsMassReb", VnVsMassBins.data()));
    }
    TH1F* hFuncPhiVsMass = reinterpret_cast<TH1F*>(hProfileFuncPhiVsMassReb->ProjectionX(histoname.Data()));
    hFuncPhiVsMass->Scale(1./resol);
    hFuncPhiVsMass->GetXaxis()->SetTitleOffset(1.1);

    delete hFuncPhiVsMass2D;
    hFuncPhiVsMass2D = NULL;
    delete hProfileFuncPhiVsMass;
    hProfileFuncPhiVsMass = NULL;
    delete hProfileFuncPhiVsMassReb;
    hProfileFuncPhiVsMassReb = NULL;

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
//function to check EP modulations
double CosnPhi(double *phi, double *pars) {

    return  pars[0]*(1+pars[1]*pars[2]*TMath::Cos(pars[2]*phi[0]));
}

//___________________________________________________________________________________//
//function to check EP modulations
double SinnPhi(double *phi, double *pars) {

    return pars[0]*(1+pars[1]*pars[2]*TMath::Sin(pars[2]*phi[0]));
}

//___________________________________________________________________________________//
//method to set plots style
void SetStyle() {
    gROOT->ForceStyle();
    gStyle->SetPadRightMargin(0.035);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.1);
    gStyle->SetPadTopMargin(0.035);
    gStyle->SetTitleSize(0.045,"xy");
    gStyle->SetLabelSize(0.04,"xy");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetHistLineWidth(2);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
}

//___________________________________________________________________________________//
//method to set plots style
void SetGraphStyle(TGraphAsymmErrors* graph, int color, int markerstyle, float markersize, int linewidth) {
    graph->SetLineColor(color);
    graph->SetMarkerColor(color);
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
