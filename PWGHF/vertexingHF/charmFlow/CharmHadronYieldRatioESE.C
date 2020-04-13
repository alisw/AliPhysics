//___________________________________________________________________________________//
// Brief: Macro for the analysis of the output of AliAnalysisTaskSECharmHadronvn     //
// Main Function: CharmHadronYieldRatioESE                                           //
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
#include "AliVertexingHFUtils.h"
#include "AliAnalysisTaskSECharmHadronvn.h"

#endif
using namespace std;

//___________________________________________________________________________________//
//method prototypes
void CharmHadronYieldRatioESE(string cfgFileName);
void ApplySelection(THnSparseF *sparse, int axisnum, double min, double max);
void ResetAxes(THnSparseF *sparse, int axisnum = -1);
TList* LoadTListFromTaskOutput(YAML::Node config);
bool LoadD0toKpiReflHistos(string reflFileName, int nPtBins, TH1F* hMCSgn[], TH1F* hMCRefl[]);
void SetStyle();
void SetGraphStyle(TGraphAsymmErrors* graph, int color, int markerstyle, float markersize = 1.5, int linewidth = 2);
void SetHistoStyle(TH1* histo, int color, int markerstyle, float markersize = 1.5, int linewidth = 2);
void DivideCanvas(TCanvas* c, const int nPtBins);

//___________________________________________________________________________________//
//main function for yield-ratios analysis

void CharmHadronYieldRatioESE(string cfgFileName) {

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
    bool useRefl = static_cast<bool>(config["AnalysisOptions"]["IncludeReflections"].as<int>());
    string reflFileName = config["AnalysisOptions"]["ReflFileName"].as<string>();
    string reflopt = config["AnalysisOptions"]["ReflOpt"].as<string>();

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
    THnSparseF* sMassVsPtVsPhiVsCentrVsqn = static_cast<THnSparseF*>(list->FindObject("fHistMassPtPhiqnCentr"));
    TH3F* hPerqnVsqnVsCentr = static_cast<TH3F*>(list->FindObject("fHistPercqnVsqnVsCentr"));

    //Load reflections file
    TH1F* hMCSgn[nPtBins];
    TH1F* hMCRefl[nPtBins];
    if(useRefl)
       useRefl = LoadD0toKpiReflHistos(reflFileName, nPtBins, hMCSgn, hMCRefl);

    TCanvas* cMassUnb = new TCanvas("cMassUnb","Unbiased",1920,1080);
    TCanvas* cMassFreeSigmaESE = new TCanvas("cMassFreeSigmaESE","ESE - free sigma",1920,1080);
    TCanvas* cMassFixSigmaESE = new TCanvas("cMassFixSigmaESE","ESE - fix sigma",1920,1080);
    TCanvas* cMassSimFitUnb = new TCanvas("cMassSimFitUnb","Unbiased - sim fit",1920,1080);
    TCanvas* cMassSimFitESE = new TCanvas("cMassSimFitESE","ESE - sim fit",1920,1080);
    DivideCanvas(cMassUnb,nPtBins);
    DivideCanvas(cMassFreeSigmaESE,nPtBins);
    DivideCanvas(cMassFixSigmaESE,nPtBins);
    DivideCanvas(cMassSimFitUnb,nPtBins);
    DivideCanvas(cMassSimFitESE,nPtBins);

    TH1D* hRawYieldUnb          = new TH1D("hRawYieldUnb",";#it{p}_{T} (GeV/#it{c}); raw yields",nPtBins,PtLims);
    TH1D* hRawYieldFreeSigmaESE = new TH1D("hRawYieldFreeSigmaESE",";#it{p}_{T} (GeV/#it{c});raw yields",nPtBins,PtLims);
    TH1D* hRawYieldFixSigmaESE  = new TH1D("hRawYieldFixSigmaESE",";#it{p}_{T} (GeV/#it{c});raw yields",nPtBins,PtLims);
    TH1D* hRawYieldSimFitUnb    = new TH1D("hRawYieldSimFitUnb",";#it{p}_{T} (GeV/#it{c});raw yields",nPtBins,PtLims);
    TH1D* hRawYieldSimFitESE    = new TH1D("hRawYieldSimFitESE",";#it{p}_{T} (GeV/#it{c});raw yields",nPtBins,PtLims);
    SetHistoStyle(hRawYieldUnb,kBlack,kFullSquare);
    SetHistoStyle(hRawYieldFreeSigmaESE,kBlue,kFullDiamond,2.);
    SetHistoStyle(hRawYieldFixSigmaESE,kRed,kFullCircle);
    SetHistoStyle(hRawYieldSimFitUnb,kBlack,kOpenSquare);
    SetHistoStyle(hRawYieldSimFitESE,kRed,kOpenCircle);

    TH1D* hMeanUnb          = new TH1D("hMeanUnb",";#it{p}_{T} (GeV/#it{c});mean (GeV/#it{c}^{2}) ",nPtBins,PtLims);
    TH1D* hMeanFreeSigmaESE = new TH1D("hMeanFreeSigmaESE",";#it{p}_{T} (GeV/#it{c});mean (GeV/#it{c}^{2})",nPtBins,PtLims);
    TH1D* hMeanFixSigmaESE  = new TH1D("hMeanFixSigmaESE",";#it{p}_{T} (GeV/#it{c});mean (GeV/#it{c}^{2})",nPtBins,PtLims);
    TH1D* hMeanSimFitUnb    = new TH1D("hMeanSimFitUnb",";#it{p}_{T} (GeV/#it{c});mean (GeV/#it{c}^{2})",nPtBins,PtLims);
    TH1D* hMeanSimFitESE    = new TH1D("hMeanSimFitESE",";#it{p}_{T} (GeV/#it{c});mean (GeV/#it{c}^{2})",nPtBins,PtLims);
    SetHistoStyle(hMeanUnb,kBlack,kFullSquare);
    SetHistoStyle(hMeanFreeSigmaESE,kBlue,kFullDiamond,2.);
    SetHistoStyle(hMeanFixSigmaESE,kRed,kFullCircle);
    SetHistoStyle(hMeanSimFitUnb,kBlack,kOpenSquare);
    SetHistoStyle(hMeanSimFitESE,kRed,kOpenCircle);

    TH1D* hSigmaUnb          = new TH1D("hSigmaUnb",";#it{p}_{T} (GeV/#it{c});width (GeV/#it{c}^{2})",nPtBins,PtLims);
    TH1D* hSigmaFreeSigmaESE = new TH1D("hSigmaFreeSigmaESE",";#it{p}_{T} (GeV/#it{c});width (GeV/#it{c}^{2})",nPtBins,PtLims);
    TH1D* hSigmaFixSigmaESE  = new TH1D("hSigmaFixSigmaESE",";#it{p}_{T} (GeV/#it{c});width (GeV/#it{c}^{2})",nPtBins,PtLims);
    TH1D* hSigmaSimFitUnb    = new TH1D("hSigmaSimFitUnb",";#it{p}_{T} (GeV/#it{c});width (GeV/#it{c}^{2})",nPtBins,PtLims);
    TH1D* hSigmaSimFitESE    = new TH1D("hSigmaSimFitESE",";#it{p}_{T} (GeV/#it{c});width (GeV/#it{c}^{2})",nPtBins,PtLims);
    SetHistoStyle(hSigmaUnb,kBlack,kFullSquare);
    SetHistoStyle(hSigmaFreeSigmaESE,kBlue,kFullDiamond,2.);
    SetHistoStyle(hSigmaFixSigmaESE,kRed,kFullCircle);
    SetHistoStyle(hSigmaSimFitUnb,kBlack,kOpenSquare);
    SetHistoStyle(hSigmaSimFitESE,kRed,kOpenCircle);

    TH1D* hRedChi2Unb          = new TH1D("hRedChi2Unb",";#it{p}_{T} (GeV/#it{c});#chi^{2} / ndf",nPtBins,PtLims);
    TH1D* hRedChi2FreeSigmaESE = new TH1D("hRedChi2FreeSigmaESE",";#it{p}_{T} (GeV/#it{c});#chi^{2} / ndf",nPtBins,PtLims);
    TH1D* hRedChi2FixSigmaESE  = new TH1D("hRedChi2FixSigmaESE",";#it{p}_{T} (GeV/#it{c});#chi^{2} / ndf",nPtBins,PtLims);
    TH1D* hRedChi2SimFitUnb    = new TH1D("hRedChi2SimFitUnb",";#it{p}_{T} (GeV/#it{c});#chi^{2} / ndf",nPtBins,PtLims);
    TH1D* hRedChi2SimFitESE    = new TH1D("hRedChi2SimFitESE",";#it{p}_{T} (GeV/#it{c});#chi^{2} / ndf",nPtBins,PtLims);
    SetHistoStyle(hRedChi2Unb,kBlack,kFullSquare);
    SetHistoStyle(hRedChi2FreeSigmaESE,kBlue,kFullDiamond,2.);
    SetHistoStyle(hRedChi2FixSigmaESE,kRed,kFullCircle);
    SetHistoStyle(hRedChi2SimFitUnb,kBlack,kOpenSquare);
    SetHistoStyle(hRedChi2SimFitESE,kRed,kOpenCircle);

    TH1F *hInvMassUnb[nPtBins], *hInvMassESE[nPtBins];
    AliHFInvMassFitter *massfitterUnb[nPtBins], *massfitterFreeSigmaESE[nPtBins], *massfitterFixSigmaESE[nPtBins], *massfitterSimFitUnb[nPtBins], *massfitterSimFitESE[nPtBins];
    for(unsigned int iPt=0; iPt<nPtBins; iPt++) {

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

        double SoverR = 0.;
        if(useRefl)
	        SoverR=(hMCRefl[iPt]->Integral(hMCRefl[iPt]->FindBin(MassMin[iPt]*1.0001),hMCRefl[iPt]->FindBin(MassMax[iPt]*0.9999)))/(hMCSgn[iPt]->Integral(hMCSgn[iPt]->FindBin(MassMin[iPt]*1.0001),hMCSgn[iPt]->FindBin(MassMax[iPt]*0.9999)));

        ApplySelection(sMassVsPtVsPhiVsCentrVsqn, 1, PtMin[iPt], PtMax[iPt]);
        if(doMLsel)
            ApplySelection(sMassVsPtVsPhiVsCentrVsqn, 9, CutValuesMLmin[iPt], CutValuesMLmax[iPt]);
        hInvMassUnb[iPt] = reinterpret_cast<TH1F*>(sMassVsPtVsPhiVsCentrVsqn->Projection(0));
        hInvMassUnb[iPt]->SetName(Form("hInvMassUnb_Pt%d", iPt));
        ApplySelection(sMassVsPtVsPhiVsCentrVsqn, 8, qnmin, qnmax);
        hInvMassESE[iPt] = reinterpret_cast<TH1F*>(sMassVsPtVsPhiVsCentrVsqn->Projection(0));
        hInvMassESE[iPt]->SetName(Form("hInvMassESE_Pt%d", iPt));

        ResetAxes(sMassVsPtVsPhiVsCentrVsqn);

        //fit unbiased
        TH1F* hMassUnbForFit = reinterpret_cast<TH1F*>(AliVertexingHFUtils::RebinHisto(hInvMassUnb[iPt],Rebin[iPt]));
        hMassUnbForFit->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/#it{c};%s;Counts per %0.f MeV/#it{c}^{2}", PtMin[iPt], PtMax[iPt], hInvMassUnb[iPt]->GetXaxis()->GetTitle(), hMassUnbForFit->GetBinWidth(1)*1000));
        massfitterUnb[iPt] = new AliHFInvMassFitter(hMassUnbForFit,MassMin[iPt],MassMax[iPt],BkgFunc,SgnFunc);
        massfitterUnb[iPt]->SetUseLikelihoodFit();
        massfitterUnb[iPt]->SetInitialGaussianMean(massD);
        massfitterUnb[iPt]->SetInitialGaussianSigma(0.01);

        if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi)
            massfitterUnb[iPt]->IncludeSecondGausPeak(massDplus,false,0.008,false);
        else if(meson==AliAnalysisTaskSECharmHadronvn::kDstartoKpipi)
            massfitterUnb[iPt]->SetInitialGaussianSigma(0.001);
        if(useRefl) {
            massfitterUnb[iPt]->SetTemplateReflections(hMCRefl[iPt],reflopt,MassMin[iPt],MassMax[iPt]);
            massfitterUnb[iPt]->SetFixReflOverS(SoverR);
        }
        massfitterUnb[iPt]->MassFitter(false);
        hRawYieldUnb->SetBinContent(iPt+1,massfitterUnb[iPt]->GetRawYield());
        hRawYieldUnb->SetBinError(iPt+1,massfitterUnb[iPt]->GetRawYieldError());
        hMeanUnb->SetBinContent(iPt+1,massfitterUnb[iPt]->GetMean());
        hMeanUnb->SetBinError(iPt+1,massfitterUnb[iPt]->GetMeanUncertainty());
        hSigmaUnb->SetBinContent(iPt+1,massfitterUnb[iPt]->GetSigma()*1000);
        hSigmaUnb->SetBinError(iPt+1,massfitterUnb[iPt]->GetSigmaUncertainty()*1000);
        hRedChi2Unb->SetBinContent(iPt+1,massfitterUnb[iPt]->GetReducedChiSquare());
        hRedChi2Unb->SetBinError(iPt+1,1.e-20);

        //free sigma
        TH1F* hMassESEForFit = reinterpret_cast<TH1F*>(AliVertexingHFUtils::RebinHisto(hInvMassESE[iPt],Rebin[iPt]));
        hMassESEForFit->SetTitle(Form("%0.f < #it{p}_{T} < %0.f GeV/#it{c};%s;Counts per %0.f MeV/#it{c}^{2}", PtMin[iPt], PtMax[iPt], hInvMassUnb[iPt]->GetXaxis()->GetTitle(), hMassUnbForFit->GetBinWidth(1)*1000));
        massfitterFreeSigmaESE[iPt] = new AliHFInvMassFitter(hMassESEForFit,MassMin[iPt],MassMax[iPt],BkgFunc,SgnFunc);
        massfitterFreeSigmaESE[iPt]->SetUseLikelihoodFit();
        massfitterFreeSigmaESE[iPt]->SetInitialGaussianMean(massD);
        massfitterFreeSigmaESE[iPt]->SetInitialGaussianSigma(0.01);

        if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi)
            massfitterFreeSigmaESE[iPt]->IncludeSecondGausPeak(massDplus,true,0.008,true);
        else if(meson==AliAnalysisTaskSECharmHadronvn::kDstartoKpipi)
            massfitterFreeSigmaESE[iPt]->SetInitialGaussianSigma(0.001);
        if(useRefl) {
            massfitterFreeSigmaESE[iPt]->SetTemplateReflections(hMCRefl[iPt],reflopt,MassMin[iPt],MassMax[iPt]);
            massfitterFreeSigmaESE[iPt]->SetFixReflOverS(SoverR);
        }
        massfitterFreeSigmaESE[iPt]->MassFitter(false);
        hRawYieldFreeSigmaESE->SetBinContent(iPt+1,massfitterFreeSigmaESE[iPt]->GetRawYield());
        hRawYieldFreeSigmaESE->SetBinError(iPt+1,massfitterFreeSigmaESE[iPt]->GetRawYieldError());
        hMeanFreeSigmaESE->SetBinContent(iPt+1,massfitterFreeSigmaESE[iPt]->GetMean());
        hMeanFreeSigmaESE->SetBinError(iPt+1,massfitterFreeSigmaESE[iPt]->GetMeanUncertainty());
        hSigmaFreeSigmaESE->SetBinContent(iPt+1,massfitterFreeSigmaESE[iPt]->GetSigma()*1000);
        hSigmaFreeSigmaESE->SetBinError(iPt+1,massfitterFreeSigmaESE[iPt]->GetSigmaUncertainty()*1000);
        hRedChi2FreeSigmaESE->SetBinContent(iPt+1,massfitterFreeSigmaESE[iPt]->GetReducedChiSquare());
        hRedChi2FreeSigmaESE->SetBinError(iPt+1,1.e-20);

        //fix sigma
        massfitterFixSigmaESE[iPt] = new AliHFInvMassFitter(hMassESEForFit,MassMin[iPt],MassMax[iPt],BkgFunc,SgnFunc);
        massfitterFixSigmaESE[iPt]->SetUseLikelihoodFit();
        massfitterFixSigmaESE[iPt]->SetInitialGaussianMean(massD);
        massfitterFixSigmaESE[iPt]->SetFixGaussianSigma(massfitterUnb[iPt]->GetSigma());

        if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi)
            massfitterFixSigmaESE[iPt]->IncludeSecondGausPeak(massDplus,true,0.008,true);
        if(useRefl) {
            massfitterFixSigmaESE[iPt]->SetTemplateReflections(hMCRefl[iPt],reflopt,MassMin[iPt],MassMax[iPt]);
            massfitterFixSigmaESE[iPt]->SetFixReflOverS(SoverR);
        }
        massfitterFixSigmaESE[iPt]->MassFitter(false);
        hRawYieldFixSigmaESE->SetBinContent(iPt+1,massfitterFixSigmaESE[iPt]->GetRawYield());
        hRawYieldFixSigmaESE->SetBinError(iPt+1,massfitterFixSigmaESE[iPt]->GetRawYieldError());
        hMeanFixSigmaESE->SetBinContent(iPt+1,massfitterFixSigmaESE[iPt]->GetMean());
        hMeanFixSigmaESE->SetBinError(iPt+1,massfitterFixSigmaESE[iPt]->GetMeanUncertainty());
        hSigmaFixSigmaESE->SetBinContent(iPt+1,massfitterFixSigmaESE[iPt]->GetSigma()*1000);
        hSigmaFixSigmaESE->SetBinError(iPt+1,massfitterFixSigmaESE[iPt]->GetSigmaUncertainty()*1000);
        hRedChi2FixSigmaESE->SetBinContent(iPt+1,massfitterFixSigmaESE[iPt]->GetReducedChiSquare());
        hRedChi2FixSigmaESE->SetBinError(iPt+1,1.e-20);

        //simultaneus fit
        massfitterSimFitUnb[iPt] = new AliHFInvMassFitter(hMassUnbForFit,MassMin[iPt],MassMax[iPt],BkgFunc,SgnFunc);
        massfitterSimFitUnb[iPt]->SetUseLikelihoodFit();
        massfitterSimFitUnb[iPt]->SetInitialGaussianMean(massD);
        massfitterSimFitUnb[iPt]->SetInitialGaussianSigma(0.01);

        if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi)
            massfitterSimFitUnb[iPt]->IncludeSecondGausPeak(massDplus,false,0.008,false);
        else if(meson==AliAnalysisTaskSECharmHadronvn::kDstartoKpipi)
            massfitterSimFitUnb[iPt]->SetInitialGaussianSigma(0.001);
        if(useRefl) {
            massfitterSimFitUnb[iPt]->SetTemplateReflections(hMCRefl[iPt],reflopt,MassMin[iPt],MassMax[iPt]);
            massfitterSimFitUnb[iPt]->SetFixReflOverS(SoverR);
        }
        massfitterSimFitUnb[iPt]->MassFitter(false);

        massfitterSimFitESE[iPt] = new AliHFInvMassFitter(hMassESEForFit,MassMin[iPt],MassMax[iPt],BkgFunc,SgnFunc);
        massfitterSimFitESE[iPt]->SetUseLikelihoodFit();
        massfitterSimFitESE[iPt]->SetInitialGaussianMean(massD);
        massfitterSimFitESE[iPt]->SetFixGaussianSigma(0.001);

        if(meson==AliAnalysisTaskSECharmHadronvn::kDstoKKpi)
            massfitterSimFitESE[iPt]->IncludeSecondGausPeak(massDplus,true,0.008,true);
        else if(meson==AliAnalysisTaskSECharmHadronvn::kDstartoKpipi)
            massfitterSimFitESE[iPt]->SetInitialGaussianSigma(0.001);
        if(useRefl) {
            massfitterSimFitESE[iPt]->SetTemplateReflections(hMCRefl[iPt],reflopt,MassMin[iPt],MassMax[iPt]);
            massfitterSimFitESE[iPt]->SetFixReflOverS(SoverR);
        }
        massfitterSimFitESE[iPt]->MassFitter(false);

        vector<unsigned int> commonpars = {static_cast<unsigned int>(massfitterUnb[iPt]->GetBackgroundFullRangeFunc()->GetNpar())+1,static_cast<unsigned int>(massfitterUnb[iPt]->GetBackgroundFullRangeFunc()->GetNpar())+2};

        ROOT::Fit::FitResult resultSimFit = AliVertexingHFUtils::DoInPlaneOutOfPlaneSimultaneusFit(massfitterSimFitUnb[iPt], massfitterSimFitESE[iPt], hMassUnbForFit, hMassESEForFit, MassMin[iPt], MassMax[iPt], massD, commonpars);
        hRawYieldSimFitUnb->SetBinContent(iPt+1,massfitterSimFitUnb[iPt]->GetSignalFunc()->GetParameter(0) / hMassUnbForFit->GetBinWidth(1));
        hRawYieldSimFitUnb->SetBinError(iPt+1,massfitterSimFitUnb[iPt]->GetSignalFunc()->GetParError(0) / hMassUnbForFit->GetBinWidth(1));
        hSigmaSimFitUnb->SetBinContent(iPt+1,massfitterSimFitUnb[iPt]->GetSignalFunc()->GetParameter(2)*1000);
        hSigmaSimFitUnb->SetBinError(iPt+1,massfitterSimFitUnb[iPt]->GetSignalFunc()->GetParError(2)*1000);
        hMeanSimFitUnb->SetBinContent(iPt+1,massfitterSimFitUnb[iPt]->GetSignalFunc()->GetParameter(1));
        hMeanSimFitUnb->SetBinError(iPt+1,massfitterSimFitUnb[iPt]->GetSignalFunc()->GetParError(1));
        hRedChi2SimFitUnb->SetBinContent(iPt+1,resultSimFit.MinFcnValue()/resultSimFit.Ndf());
        hRedChi2SimFitUnb->SetBinError(iPt+1,1.e-20);
        hRawYieldSimFitESE->SetBinContent(iPt+1,massfitterSimFitESE[iPt]->GetSignalFunc()->GetParameter(0) / hMassESEForFit->GetBinWidth(1));
        hRawYieldSimFitESE->SetBinError(iPt+1,massfitterSimFitESE[iPt]->GetSignalFunc()->GetParError(0) / hMassESEForFit->GetBinWidth(1));
        hSigmaSimFitESE->SetBinContent(iPt+1,massfitterSimFitESE[iPt]->GetSignalFunc()->GetParameter(2)*1000);
        hSigmaSimFitESE->SetBinError(iPt+1,massfitterSimFitESE[iPt]->GetSignalFunc()->GetParError(2)*1000);
        hMeanSimFitESE->SetBinContent(iPt+1,massfitterSimFitESE[iPt]->GetSignalFunc()->GetParameter(1));
        hMeanSimFitESE->SetBinError(iPt+1,massfitterSimFitESE[iPt]->GetSignalFunc()->GetParError(1));
        hRedChi2SimFitESE->SetBinContent(iPt+1,resultSimFit.MinFcnValue()/resultSimFit.Ndf());
        hRedChi2SimFitESE->SetBinError(iPt+1,1.e-20);

        commonpars.clear();

        if(nPtBins>1)
            cMassUnb->cd(iPt+1);
        else
            cMassUnb->cd();
        massfitterUnb[iPt]->DrawHere(gPad);

        if(nPtBins>1)
            cMassFreeSigmaESE->cd(iPt+1);
        else
            cMassFreeSigmaESE->cd();
        massfitterFreeSigmaESE[iPt]->DrawHere(gPad);

        if(nPtBins>1)
            cMassFixSigmaESE->cd(iPt+1);
        else
            cMassFixSigmaESE->cd();
        massfitterFixSigmaESE[iPt]->DrawHere(gPad);

        if(nPtBins>1)
            cMassSimFitUnb->cd(iPt+1);
        else
            cMassSimFitUnb->cd();
        massfitterSimFitUnb[iPt]->DrawHere(gPad,3,1);

        if(nPtBins>1)
            cMassSimFitESE->cd(iPt+1);
        else
            cMassSimFitESE->cd();
        massfitterSimFitESE[iPt]->DrawHere(gPad,3,1);
    }

    //Get normalisation
    TH1F* hPerqn = reinterpret_cast<TH1F*>(hPerqnVsqnVsCentr->ProjectionZ());
    int qnbinmin = hPerqn->GetXaxis()->FindBin(qnmin*1.00001);
    int qnbinmax = hPerqn->GetXaxis()->FindBin(qnmax*0.99999);
    double nEvUnb = hPerqn->Integral();
    double nEvESE = hPerqn->Integral(qnbinmin, qnbinmax);
    double rho = TMath::Sqrt(nEvESE/nEvUnb);
    TH1D* hNorm = new TH1D("hNorm",";;number of events",2,0.5,2.5);
    hNorm->GetXaxis()->SetBinLabel(1,"unbiased");
    hNorm->GetXaxis()->SetBinLabel(2,"ESE");
    hNorm->SetBinContent(1,nEvUnb);
    hNorm->SetBinContent(2,nEvESE);

    TH1D* hNormRawYieldUnb          = static_cast<TH1D*>(hRawYieldUnb->Clone("hNormRawYieldUnb"));
    TH1D* hNormRawYieldFreeSigmaESE = static_cast<TH1D*>(hRawYieldFreeSigmaESE->Clone("hNormRawYieldFreeSigmaESE"));
    TH1D* hNormRawYieldFixSigmaESE  = static_cast<TH1D*>(hRawYieldFixSigmaESE->Clone("hNormRawYieldFixSigmaESE"));
    TH1D* hNormRawYieldSimFitUnb    = static_cast<TH1D*>(hRawYieldSimFitUnb->Clone("hNormRawYieldSimFitUnb"));
    TH1D* hNormRawYieldSimFitESE    = static_cast<TH1D*>(hRawYieldSimFitESE->Clone("hNormRawYieldSimFitESE"));
    hNormRawYieldUnb->GetYaxis()->SetTitle("raw yield / N_{events}");
    hNormRawYieldFreeSigmaESE->GetYaxis()->SetTitle("raw yield / N_{events}");
    hNormRawYieldFixSigmaESE->GetYaxis()->SetTitle("raw yield / N_{events}");
    hNormRawYieldSimFitUnb->GetYaxis()->SetTitle("raw yield / N_{events}");
    hNormRawYieldSimFitESE->GetYaxis()->SetTitle("raw yield / N_{events}");
    hNormRawYieldUnb->Scale(1./nEvUnb);
    hNormRawYieldFreeSigmaESE->Scale(1./nEvESE);
    hNormRawYieldFixSigmaESE->Scale(1./nEvESE);
    hNormRawYieldSimFitUnb->Scale(1./nEvUnb);
    hNormRawYieldSimFitESE->Scale(1./nEvESE);

    //Compute ratios
    TH1D* hRatioRawYieldFreeSigmaESE = static_cast<TH1D*>(hNormRawYieldFreeSigmaESE->Clone("hRatioRawYieldFreeSigmaESE"));
    TH1D* hRatioRawYieldFixSigmaESE  = static_cast<TH1D*>(hNormRawYieldFixSigmaESE->Clone("hRatioRawYieldFixSigmaESE"));
    TH1D* hRatioRawYieldSimFitESE    = static_cast<TH1D*>(hNormRawYieldSimFitESE->Clone("hRatioRawYieldSimFitESE"));
    hRatioRawYieldFreeSigmaESE->GetYaxis()->SetTitle("ratio ESE / unbiased");
    hRatioRawYieldFixSigmaESE->GetYaxis()->SetTitle("ratio ESE / unbiased");
    hRatioRawYieldSimFitESE->GetYaxis()->SetTitle("ratio ESE / unbiased");
    hRatioRawYieldFreeSigmaESE->Divide(hNormRawYieldFreeSigmaESE,hNormRawYieldUnb,1.,1.);
    hRatioRawYieldFixSigmaESE->Divide(hRatioRawYieldFixSigmaESE,hNormRawYieldUnb,1.,1.);
    hRatioRawYieldSimFitESE->Divide(hNormRawYieldSimFitESE,hNormRawYieldSimFitUnb,1.,1.);
    for(unsigned int iPt=0; iPt<nPtBins; iPt++) {
        double reluncUnb          = hRawYieldUnb->GetBinError(iPt+1)/hRawYieldUnb->GetBinContent(iPt+1);
        double reluncSimFitUnb    = hRawYieldSimFitUnb->GetBinError(iPt+1)/hRawYieldSimFitUnb->GetBinContent(iPt+1);
        double reluncFreeSigmaESE = hRawYieldFreeSigmaESE->GetBinError(iPt+1)/hRawYieldFreeSigmaESE->GetBinContent(iPt+1);
        double reluncFixSigmaESE  = hRawYieldFixSigmaESE->GetBinError(iPt+1)/hRawYieldFixSigmaESE->GetBinContent(iPt+1);
        double reluncSimFitESE    = hRawYieldSimFitESE->GetBinError(iPt+1)/hRawYieldSimFitESE->GetBinContent(iPt+1);
        double ratioFreeSigmaESE  = hRatioRawYieldFreeSigmaESE->GetBinContent(iPt+1);
        double ratioFixSigmaESE   = hRatioRawYieldFixSigmaESE->GetBinContent(iPt+1);
        double ratioSimFitESE     = hRatioRawYieldSimFitESE->GetBinContent(iPt+1);

        hRatioRawYieldFreeSigmaESE->SetBinError(iPt+1,ratioFreeSigmaESE*TMath::Sqrt(reluncFreeSigmaESE*reluncFreeSigmaESE+reluncUnb*reluncUnb-2*rho*reluncFreeSigmaESE*reluncUnb));
        hRatioRawYieldFixSigmaESE->SetBinError(iPt+1,ratioFixSigmaESE*TMath::Sqrt(reluncFixSigmaESE*reluncFixSigmaESE+reluncUnb*reluncUnb-2*rho*reluncFixSigmaESE*reluncUnb));
        hRatioRawYieldSimFitESE->SetBinError(iPt+1,ratioSimFitESE*TMath::Sqrt(reluncSimFitESE*reluncSimFitESE+reluncUnb*reluncUnb-2*rho*reluncSimFitESE*reluncUnb));
    }

    //plots
    TLine* lineatone = new TLine(PtLims[0],1,PtLims[nPtBins],1);
    lineatone->SetLineColor(kBlack);
    lineatone->SetLineWidth(2);
    lineatone->SetLineStyle(7);

    TLegend* leg = new TLegend(0.6,0.68,0.89,0.88);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(hRawYieldUnb,"Unbiased");
    leg->AddEntry(hRawYieldFreeSigmaESE,"ESE - free sigma");
    leg->AddEntry(hRawYieldFixSigmaESE,"ESE - fix sigma");
    leg->AddEntry(hRawYieldSimFitUnb,"unbiased - sim fit");
    leg->AddEntry(hRawYieldSimFitESE,"ESE - sim fit");

    TLegend* legRatio = new TLegend(0.6,0.73,0.89,0.88);
    legRatio->SetFillStyle(0);
    legRatio->SetBorderSize(0);
    legRatio->SetTextSize(0.04);
    legRatio->AddEntry(hRatioRawYieldFreeSigmaESE,"Free sigma");
    legRatio->AddEntry(hRatioRawYieldFixSigmaESE,"Fix sigma");
    legRatio->AddEntry(hRatioRawYieldSimFitESE,"Simultaneus fit");

    TCanvas* cMean = new TCanvas("cMean","",800,800);
    hMeanUnb->Draw();
    hMeanFreeSigmaESE->Draw("same");
    hMeanFixSigmaESE->Draw("same");
    hMeanSimFitUnb->Draw("same");
    hMeanSimFitESE->Draw("same");
    leg->Draw();

    TCanvas* cSigma = new TCanvas("cSigma","",800,800);
    hSigmaUnb->GetYaxis()->SetRangeUser(0.,hSigmaUnb->GetMaximum()*2.);
    hSigmaUnb->Draw();
    hSigmaFreeSigmaESE->Draw("same");
    hSigmaFixSigmaESE->Draw("same");
    hSigmaSimFitUnb->Draw("same");
    hSigmaSimFitESE->Draw("same");
    leg->Draw();

    TCanvas* cRedChi2 = new TCanvas("cRedChi2","",800,800);
    cRedChi2->DrawFrame(PtLims[0],0.,PtLims[nPtBins],4.,";#it{p}_{T} GeV/#it{c};#chi^{2} / ndf");
    hRedChi2Unb->Draw("same");
    hRedChi2FreeSigmaESE->Draw("same");
    hRedChi2FixSigmaESE->Draw("same");
    hRedChi2SimFitUnb->Draw("same");
    hRedChi2SimFitESE->Draw("same");
    leg->Draw();

    TCanvas* cRawYield = new TCanvas("cRawYield","",800,800);
    hRawYieldUnb->GetYaxis()->SetRangeUser(0.,hRawYieldUnb->GetMaximum()*1.3);
    hRawYieldUnb->Draw();
    hRawYieldFreeSigmaESE->Draw("same");
    hRawYieldFixSigmaESE->Draw("same");
    hRawYieldSimFitUnb->Draw("same");
    hRawYieldSimFitESE->Draw("same");
    leg->Draw();

    TCanvas* cNormRawYield = new TCanvas("cNormRawYield","",800,800);
    hNormRawYieldUnb->GetYaxis()->SetRangeUser(0.,hNormRawYieldUnb->GetMaximum()*1.3);
    hNormRawYieldUnb->Draw();
    hNormRawYieldFreeSigmaESE->Draw("same");
    hNormRawYieldFixSigmaESE->Draw("same");
    hNormRawYieldSimFitUnb->Draw("same");
    hNormRawYieldSimFitESE->Draw("same");
    leg->Draw();

    TCanvas* cRatio = new TCanvas("cRatio","",800,800);
    cRatio->DrawFrame(PtLims[0],0.,PtLims[nPtBins],2.,";#it{p}_{T} GeV/#it{c};ratio ESE / unbiased");
    lineatone->Draw("same");
    hRatioRawYieldFreeSigmaESE->Draw("same");
    hRatioRawYieldFixSigmaESE->Draw("same");
    hRatioRawYieldSimFitESE->Draw("same");
    legRatio->Draw();

    //output files
    string outputdir = config["OutputDir"]["Analysis"].as<string>();
    cMassUnb->SaveAs(Form("%s/InvMassFitsUnbiased%s_q%d_%0.f-%0.f.pdf",outputdir.data(),mesonname.data(),harmonic,qnmin,qnmax));
    cMassFreeSigmaESE->SaveAs(Form("%s/InvMassFitsFreeSigmaESE%s_q%d_%0.f-%0.f.pdf",outputdir.data(),mesonname.data(),harmonic,qnmin,qnmax));
    cMassFixSigmaESE->SaveAs(Form("%s/InvMassFitsFixSigmaESE%s_q%d_%0.f-%0.f.pdf",outputdir.data(),mesonname.data(),harmonic,qnmin,qnmax));
    cMassSimFitUnb->SaveAs(Form("%s/InvMassFitsSimFitUnb%s_q%d_%0.f-%0.f.pdf",outputdir.data(),mesonname.data(),harmonic,qnmin,qnmax));
    cMassSimFitESE->SaveAs(Form("%s/InvMassFitsSimFitESE%s_q%d_%0.f-%0.f.pdf",outputdir.data(),mesonname.data(),harmonic,qnmin,qnmax));
    cNormRawYield->SaveAs(Form("%s/%sNormRawYield_q%d_%0.f-%0.f.pdf",outputdir.data(),mesonname.data(),harmonic,qnmin,qnmax));
    cRatio->SaveAs(Form("%s/%sRawYieldRatio_q%d_%0.f-%0.f.pdf",outputdir.data(),mesonname.data(),harmonic,qnmin,qnmax));

    TFile outfile(Form("%s/%sRawYieldRatio_q%d_%0.f-%0.f.root",outputdir.data(),mesonname.data(),harmonic,qnmin,qnmax),"recreate");
    cRatio->Write();
    cNormRawYield->Write();
    hRatioRawYieldFreeSigmaESE->Write();
    hRatioRawYieldFixSigmaESE->Write();
    hRatioRawYieldSimFitESE->Write();
    hNormRawYieldUnb->Write();
    hNormRawYieldFreeSigmaESE->Write();
    hNormRawYieldFixSigmaESE->Write();
    hNormRawYieldSimFitUnb->Write();
    hNormRawYieldSimFitESE->Write();
    TDirectoryFile dirFits("Fits","Fits");
    dirFits.Write();
    dirFits.cd();
    cMassUnb->Write();
    cMassFreeSigmaESE->Write();
    cMassFixSigmaESE->Write();
    cMassSimFitUnb->Write();
    cMassSimFitESE->Write();
    cRawYield->Write();
    cSigma->Write();
    cMean->Write();
    cRedChi2->Write();
    hRawYieldUnb->Write();
    hRawYieldFreeSigmaESE->Write();
    hRawYieldFixSigmaESE->Write();
    hRawYieldSimFitUnb->Write();
    hRawYieldSimFitESE->Write();
    hMeanUnb->Write();
    hMeanFreeSigmaESE->Write();
    hMeanFixSigmaESE->Write();
    hMeanSimFitUnb->Write();
    hMeanSimFitESE->Write();
    hSigmaUnb->Write();
    hSigmaFreeSigmaESE->Write();
    hSigmaFixSigmaESE->Write();
    hSigmaSimFitUnb->Write();
    hSigmaSimFitESE->Write();
    hRedChi2Unb->Write();
    hRedChi2FreeSigmaESE->Write();
    hRedChi2FixSigmaESE->Write();
    hRedChi2SimFitUnb->Write();
    hRedChi2SimFitESE->Write();
    outfile.cd();
    TDirectoryFile dirMass("InvMassSpectra","InvMassSpectra");
    dirMass.Write();
    dirMass.cd();
    for(unsigned int iPt=0; iPt<nPtBins; iPt++) {
        hInvMassUnb[iPt]->Write();
        hInvMassESE[iPt]->Write();
    }
    outfile.cd();
    TDirectoryFile dirNorm("Norm","Norm");
    dirNorm.Write();
    dirNorm.cd();
    hNorm->Write();
    outfile.Close();
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

//___________________________________________________________________________________//
//method to set plots style
void SetStyle() {
    gROOT->ForceStyle();
    gStyle->SetPadRightMargin(0.035);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.1);
    gStyle->SetPadTopMargin(0.07);
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

//___________________________________________________________________________________//
//method to divide the canvas given a number of pT bins
void DivideCanvas(TCanvas* c, const int nPtBins) {
  if(nPtBins<2)
    c->Divide(1,1);
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
