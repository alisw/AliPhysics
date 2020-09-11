#if !defined(__CINT__) || defined(__CLING__)

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <fstream>

#include "yaml-cpp/yaml.h"

#include <Riostream.h>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>
#include <TList.h>
#include <TMath.h>
#include <TFractionFitter.h>
#include <TObjArray.h>
#include <TDatabasePDG.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TDirectoryFile.h>
#include <TExec.h>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"

#include "AliAnalysisTaskSEHFSystPID.h"

#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// \brief: macro for estimation of PID Systematic uncertainty of the single tracks (pions/kaons/protons)                                            //
//         it needs a config input file such as https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/macros/config_singletrack_pid.yml //
//         it requires ROOT6 for the usage of RDataFrame                                                                                            //
// \main function: AnalysePIDTree(TString cfgFileName)                                                                                              //
// \author: F. Grosa, fabrizio.grosa@cern.ch                                                                                                        //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________
//GLOBAL VARIABLES
enum projVars
{
    kPt,
    kP
};
enum partPID
{
    kEl = 11,
    kMuon = 13,
    kPion = 211,
    kKaon = 321,
    kPr = 2212,
    kAll = -100
};
std::map<int, std::string> pdgNames = {{kEl, "Electron"}, {kMuon, "Muon"}, {kPion, "Pion"}, {kKaon, "Kaon"}, {kPr, "Proton"}, {kAll, "All"}};
std::map<int, int> pdgPosition = {{kEl, 0}, {kMuon, 1}, {kPion, 2}, {kKaon, 3}, {kPr, 4}, {kAll, 5}};
std::map<int, int> pdgColors = {{kEl, kOrange+7}, {kMuon, kGray+1}, {kPion, kRed+1}, {kKaon, kAzure+4}, {kPr, kGreen+2}, {kAll, kBlack}};
std::map<int, int> pdgFillColors = {{kEl, kOrange+7}, {kMuon, kGray+1}, {kPion, kRed+1}, {kKaon, kAzure+4}, {kPr, kGreen+2}, {kAll, kWhite}};
std::map<std::string, int> projVarMap = {{"kP", kP}, {"kPt", kPt}};

const unsigned int nBinsMax = 100;
const unsigned int nEtaBinsMax = 20;
const unsigned int nMaxEff = 20;

//_____________________________________________________
//METHOD PROTOTYPES
void AnalysePIDTree(TString cfgFileName="config_singletrack_pid.yml");
void PerformPIDAnalysis(std::string inFileNameData, std::string inDirNameData, std::string inListNameData,
                        std::string inFileNameMC, std::string inDirNameMC, std::string inListNameMC,
                        std::string outDirName, TString varTitle, int var4proj,
                        unsigned int nBins, std::vector<double> binMins, std::vector<double> binMaxs, double binLims[], std::vector<TString> binLabels,
                        unsigned int nEtaBins, std::vector<double> absEtaBinMins, std::vector<double> absEtaBinMaxs, double binEtaLims[], std::vector<TString> etaBinLabels, YAML::Node config);
void PerformTPCTOFmatchingAnalysis(std::string inFileNameData, std::string inDirNameData, std::string inListNameData,
                                   std::string inFileNameMC, std::string inDirNameMC, std::string inListNameMC,
                                   std::string outDirName, TString varTitle, int var4proj, unsigned int nBins, double binLims[],
                                   unsigned int nEtaBins, std::vector<double> absEtaBinMins, std::vector<double> absEtaBinMaxs, std::vector<TString> etaBinLabels);
void ComputeEfficiency(double num, double den, double &eff, double &effunc);
void GetTOFFractionsFromData(int whichpart, unsigned int iBin, std::map<int, TH1D*> hFracMC, std::map<int, TH1D*> hFracData, std::map<int, TH1D*> hNsigmaMC,
                             TH1D *hNsigmaData, TFractionFitter *&fNsigmaFitter, std::vector<int> &templUsed);
double PDFnsigmaTPCtot(double *nsigma, double *pars);
void PlotQAhistos(TList *listMC, TList *listData, string outDirName);
void DivideCanvas(TCanvas *c, int nBins);
void SetStyle();
void SetTH1Style(TH1 *histo, int markerstyle, int markercolor, float markersize, int linewidth, int linecolor,
                 int fillcolor, float labelsize = -1, float titlesize = -1);

//_____________________________________________________
//METHOD IMPLEMENTATIONS
void AnalysePIDTree(TString cfgFileName)
{

    SetStyle();

    //Load configs from yaml file
    YAML::Node config = YAML::LoadFile(cfgFileName.Data());
    if (config.IsNull())
    {
        std::cerr << "\033[31mERROR: yaml config file not found! Exit\033[0m" << std::endl;
        return;
    }
    else 
    {
        std::cout << "\n\n*******************************************" << std::endl;
        std::cout << Form("\033[32mLoading configuration from file %s\033[0m\n", cfgFileName.Data()) << std::endl;
    }

    // input files and output dir
    std::string inFileNameData = config["inputs"]["data"]["filename"].as<std::string>();
    std::string inDirNameData = config["inputs"]["data"]["dirname"].as<std::string>();
    std::string inListNameData = config["inputs"]["data"]["listname"].as<std::string>();
    std::cout << "\033[32mLoading data tree from\033[0m" << std::endl;
    std::cout << Form("\tfile: %s", inFileNameData.data()) << std::endl;
    std::cout << Form("\t\tdir: %s", inDirNameData.data()) << std::endl;
    std::cout << Form("\t\t\tlist: %s", inListNameData.data()) << std::endl;
    
    std::string inFileNameMC = config["inputs"]["MC"]["filename"].as<std::string>();
    std::string inDirNameMC = config["inputs"]["MC"]["dirname"].as<std::string>();
    std::string inListNameMC = config["inputs"]["MC"]["listname"].as<std::string>();
    std::cout << "\033[32mLoading MC tree from\033[0m" << std::endl;
    std::cout << Form("\tfile: %s", inFileNameMC.data()) << std::endl;
    std::cout << Form("\t\tdir: %s", inDirNameMC.data()) << std::endl;
    std::cout << Form("\t\t\tlist: %s", inListNameMC.data()) << std::endl;

    std::string outDirName = config["output"]["dirname"].as<std::string>();

    // binning related quantities
    int var4proj = projVarMap[config["binning"]["bins"]["var4proj"].as<std::string>()];
    TString varTitle = "", varName = "";
    if (var4proj == kPt)
    {
        varName = "pT";
        varTitle = "#it{p}_{T}";
    }
    else if (var4proj == kP)
    {
        varName = "p";
        varTitle = "#it{p}";
    }
    else
    {
        std::cerr << "\033[31mERROR: you can chose to project the Nsigma distributions vs. pT or p only. Exit\033[0m" << std::endl;
        return;
    }
    std::vector<double> binMins = config["binning"]["bins"]["mins"].as<std::vector<double> >();
    std::vector<double> binMaxs = config["binning"]["bins"]["maxs"].as<std::vector<double> >();
    unsigned int nBins = binMins.size();
    if(nBins > nBinsMax)
    {
        std::cout << Form("\033[33mWARNING: number of bins larger than maximum allowed of %d. Truncating last %d bins\033[0m", nBinsMax, nBinsMax-nBins) << std::endl;
        nBins = nBinsMax;
    }
    std::cout << nBins << Form(" %s bins:", varName.Data()) << std::endl;
    for (unsigned int iBin = 0; iBin < nBins; iBin++)
        std::cout << binMins[iBin] << "-" << binMaxs[iBin] << " GeV/c" << std::endl;
    std::vector<TString> binLabels;
    double binLims[nBins+1];
    for(unsigned int iBin = 0; iBin < nBins; iBin++)
    {
        binLabels.push_back(Form("%s%.0f_%.0f", varName.Data(), binMins[iBin]*10, binMaxs[iBin]*10));
        binLims[iBin] = binMins[iBin];
    }
    binLims[nBins] = binMaxs[nBins-1];
    
    std::vector<double> absEtaBinMins = config["binning"]["absetabins"]["mins"].as<std::vector<double> >();
    std::vector<double> absEtaBinMaxs = config["binning"]["absetabins"]["maxs"].as<std::vector<double> >();
    unsigned int nEtaBins = absEtaBinMins.size();
    std::cout << "\n";
    if(nEtaBins > nEtaBinsMax)
    {
        std::cout << Form("\033[33mWARNING: number of eta bins larger than maximum allowed of %d. Truncating last %d eta bins\033[0m", nBinsMax, nBinsMax-nEtaBinsMax) << std::endl;
        nEtaBins = nEtaBinsMax;
    }
    double binEtaLims[nEtaBins+1];
    for(unsigned int iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++)
        binEtaLims[iEtaBin] = absEtaBinMins[iEtaBin];
    binEtaLims[nEtaBins] = absEtaBinMaxs[nEtaBins-1];

    if(nEtaBins > 1)
    {
        //add also integrated eta bin
        absEtaBinMins.push_back(absEtaBinMins[0]);
        absEtaBinMaxs.push_back(absEtaBinMaxs[nEtaBins-1]);
        nEtaBins++;
        std::cout << "\033[32mAdding integrated eta bin " << absEtaBinMins[nEtaBins-1] << "-" << absEtaBinMaxs[nEtaBins-1] << "\033[0m" <<std::endl;
    }
    std::vector<TString> etaBinLabels;
    if(nEtaBins > 1)
        std::cout << nEtaBins << " eta bins (data only):" << std::endl;
    else
        std::cout << nEtaBins << " eta bins:" << std::endl;
    for (unsigned int iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++)
    {
        if(iEtaBin < nEtaBins-1)
            etaBinLabels.push_back(Form("abseta%.0f_%.0f", absEtaBinMins[iEtaBin]*10, absEtaBinMaxs[iEtaBin]*10));
        else
            etaBinLabels.push_back("etaint");
        std::cout << absEtaBinMins[iEtaBin] << "-" << absEtaBinMaxs[iEtaBin] << std::endl;
    }

    // PID analysis
    bool doPIDstudies = static_cast<bool>(config["PIDstudies"]["enable"].as<int>());
    if(doPIDstudies)
    {
        std::cout << "\n*******************************************\n" << std::endl;
        std::cout << "\e[1m\033[32mStarting PID analysis\033[0m\e[0m\n" << std::endl;

        PerformPIDAnalysis(inFileNameData, inDirNameData, inListNameData, inFileNameMC, inDirNameMC, inListNameMC,
                           outDirName, varTitle, var4proj, nBins, binMins, binMaxs, binLims, binLabels, 
                           nEtaBins, absEtaBinMins, absEtaBinMaxs, binEtaLims, etaBinLabels, config);
    }

    // TPC-TOF matching analysis
    bool doTOFTPCmatchingStudies = static_cast<bool>(config["TPCTOFmatching"]["enable"].as<int>());
    if(doTOFTPCmatchingStudies)
    {
        std::cout << "\n*******************************************\n" << std::endl;
        std::cout << "\e[1m\033[32mStarting TPC-TOF matching analysis\033[0m\e[0m\n" << std::endl;

        PerformTPCTOFmatchingAnalysis(inFileNameData, inDirNameData, inListNameData, inFileNameMC, inDirNameMC, inListNameMC,
                                      outDirName, varTitle, var4proj, nBins, binLims, nEtaBins, absEtaBinMins, absEtaBinMaxs, etaBinLabels);
    }
}

//______________________________________________________
void PerformPIDAnalysis(std::string inFileNameData, std::string inDirNameData, std::string inListNameData,
                        std::string inFileNameMC, std::string inDirNameMC, std::string inListNameMC,
                        std::string outDirName, TString varTitle, int var4proj,
                        unsigned int nBins, std::vector<double> binMins, std::vector<double> binMaxs, double binLims[], std::vector<TString> binLabels,
                        unsigned int nEtaBins, std::vector<double> absEtaBinMins, std::vector<double> absEtaBinMaxs, double binEtaLims[], std::vector<TString> etaBinLabels, YAML::Node config)
{
    // quantities for PID studies
    bool produceQAplots = static_cast<bool>(config["PIDstudies"]["produceQAplots"].as<int>());

    std::vector<int> nSigma4Eff = config["PIDstudies"]["PIDefficiency"]["nSigma"].as<std::vector<int> >();
    std::vector<int> markersEffMC = config["PIDstudies"]["PIDefficiency"]["markersEffMC"].as<std::vector<int> >();
    std::vector<int> markersEffData = config["PIDstudies"]["PIDefficiency"]["markersEffData"].as<std::vector<int> >();
    unsigned int nEff = nSigma4Eff.size();
    if(nEff != markersEffMC.size() || nEff != markersEffData.size())
    {
        std::cerr << "\n\033[31mERROR: number of Nsigma values for efficiencies and number of marker styles not consistent! Exit\033[0m" << std::endl;
        return;
    }
    else
    {
        std::cout << "\nNsigma values to be tested for PID efficiency:" << std::endl;
        for (unsigned int iEff = 0; iEff < nEff; iEff++)
            std::cout << nSigma4Eff[iEff] << std::endl;
    }

    //define histos
    //MC truth
    TH2D *hNsigmaTPCPionVsPtMCTrue, *hNsigmaTPCKaonVsPtMCTrue, *hNsigmaTPCProtonVsPtMCTrue;
    TH2D *hNsigmaTOFPionVsPtMCTrue, *hNsigmaTOFKaonVsPtMCTrue, *hNsigmaTOFProtonVsPtMCTrue;

    // MC truth (one for each p bin)
    std::array<TH1D*, nBinsMax> hNsigmaTPCPionMCTrue, hNsigmaTPCKaonMCTrue, hNsigmaTPCProtonMCTrue;
    std::array<TH1D*, nBinsMax> hNsigmaTOFPionMCTrue, hNsigmaTOFKaonMCTrue, hNsigmaTOFProtonMCTrue;

    //data tagged (one for each p and eta bin)
    std::array<std::array<TH1D*, nBinsMax>, nEtaBinsMax+1> hNsigmaTPCPionDataV0tag, hNsigmaTPCKaonDataKinktag, hNsigmaTPCKaonDataTOFtag, hNsigmaTPCProtonDataV0tag;
    std::array<std::array<TH1D*, nBinsMax>, nEtaBinsMax+1> hNsigmaTOFPionDataV0tag, hNsigmaTOFKaonDataKinktag, hNsigmaTOFKaonDataTPCtag, hNsigmaTOFProtonDataV0tag;

    //MC tagged (one for each p bin, MC no eta differential)
    std::array<std::map<int, TH1D*>, nBinsMax> hNsigmaTPCPionMCV0tag, hNsigmaTPCKaonMCKinktag, hNsigmaTPCKaonMCTOFtag, hNsigmaTPCProtonMCV0tag;
    std::array<std::map<int, TH1D*>, nBinsMax> hNsigmaTOFPionMCV0tag, hNsigmaTOFKaonMCKinktag, hNsigmaTOFKaonMCTPCtag, hNsigmaTOFProtonMCV0tag;

    //load MC inputs
    auto infileMC = TFile::Open(inFileNameMC.data());
    if (!infileMC || !infileMC->IsOpen())
        return;
    auto indirMC = static_cast<TDirectoryFile *>(infileMC->Get(inDirNameMC.data()));
    if (!indirMC)
    {
        std::cerr << Form("\033[31mERROR: TDirectoryFile %s not found in input file for MC! Exit\033[0m", inDirNameMC.data()) << std::endl;
        return;
    }
    auto listMC = static_cast<TList *>(indirMC->Get(inListNameMC.data()));
    if (!listMC)
    {
        std::cerr << Form("\033[31mERROR: TList %s not found in input file for MC! Exit\033[0m", inListNameMC.data()) << std::endl;
        return;
    }
    hNsigmaTPCPionVsPtMCTrue = static_cast<TH2D*>(listMC->FindObject("fHistNsigmaTPCvsPt_Pion"));
    hNsigmaTPCKaonVsPtMCTrue = static_cast<TH2D*>(listMC->FindObject("fHistNsigmaTPCvsPt_Kaon"));
    hNsigmaTPCProtonVsPtMCTrue = static_cast<TH2D*>(listMC->FindObject("fHistNsigmaTPCvsPt_Proton"));
    hNsigmaTOFPionVsPtMCTrue = static_cast<TH2D*>(listMC->FindObject("fHistNsigmaTOFvsPt_Pion"));
    hNsigmaTOFKaonVsPtMCTrue = static_cast<TH2D*>(listMC->FindObject("fHistNsigmaTOFvsPt_Kaon"));
    hNsigmaTOFProtonVsPtMCTrue = static_cast<TH2D*>(listMC->FindObject("fHistNsigmaTOFvsPt_Proton"));
    for (unsigned int iBin = 0; iBin < nBins; iBin++)
    {
        int ptBinMin = hNsigmaTPCPionVsPtMCTrue->GetXaxis()->FindBin(binMins[iBin] * 1.0001);
        int ptBinMax = hNsigmaTPCPionVsPtMCTrue->GetXaxis()->FindBin(binMaxs[iBin] * 0.9999);

        hNsigmaTPCPionMCTrue[iBin] = static_cast<TH1D*>(hNsigmaTPCPionVsPtMCTrue->ProjectionY(Form("hNsigmaTPCPionMCTrue_p%0.f_%0.f", binMins[iBin]*10, binMaxs[iBin]*10), ptBinMin, ptBinMax));
        hNsigmaTPCKaonMCTrue[iBin] = static_cast<TH1D*>(hNsigmaTPCKaonVsPtMCTrue->ProjectionY(Form("hNsigmaTPCKaonMCTrue_p%0.f_%0.f", binMins[iBin]*10, binMaxs[iBin]*10), ptBinMin, ptBinMax));
        hNsigmaTPCProtonMCTrue[iBin] = static_cast<TH1D*>(hNsigmaTPCProtonVsPtMCTrue->ProjectionY(Form("hNsigmaTPCProtonMCTrue_p%0.f_%0.f", binMins[iBin]*10, binMaxs[iBin]*10), ptBinMin, ptBinMax));

        hNsigmaTOFPionMCTrue[iBin] = static_cast<TH1D*>(hNsigmaTOFPionVsPtMCTrue->ProjectionY(Form("hNsigmaTOFPionMCTrue_p%0.f_%0.f", binMins[iBin]*10, binMaxs[iBin]*10), ptBinMin, ptBinMax));
        hNsigmaTOFKaonMCTrue[iBin] = static_cast<TH1D*>(hNsigmaTOFKaonVsPtMCTrue->ProjectionY(Form("hNsigmaTOFKaonMCTrue_p%0.f_%0.f", binMins[iBin]*10, binMaxs[iBin]*10), ptBinMin, ptBinMax));
        hNsigmaTOFProtonMCTrue[iBin] = static_cast<TH1D*>(hNsigmaTOFProtonVsPtMCTrue->ProjectionY(Form("hNsigmaTOFProtonMCTrue_p%0.f_%0.f", binMins[iBin]*10, binMaxs[iBin]*10), ptBinMin, ptBinMax));

        SetTH1Style(hNsigmaTPCPionMCTrue[iBin], kFullCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], pdgFillColors[kAll], 0.055, 0.06);
        SetTH1Style(hNsigmaTPCKaonMCTrue[iBin], kFullCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], pdgFillColors[kAll], 0.055, 0.06);
        SetTH1Style(hNsigmaTPCProtonMCTrue[iBin], kFullCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], pdgFillColors[kAll], 0.055, 0.06);
        SetTH1Style(hNsigmaTOFPionMCTrue[iBin], kFullCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], pdgFillColors[kAll], 0.055, 0.06);
        SetTH1Style(hNsigmaTOFKaonMCTrue[iBin], kFullCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], pdgFillColors[kAll], 0.055, 0.06);
        SetTH1Style(hNsigmaTOFProtonMCTrue[iBin], kFullCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], pdgFillColors[kAll], 0.055, 0.06);
    }

    std::cout << "\n*******************************************\n" << std::endl;
    std::cout << "\033[32mProject MC tree\033[0m\n" << std::endl;
    ROOT::EnableImplicitMT(); //tell ROOT to go parallel
    ROOT::RDataFrame dataFrameMC(Form("%s/%s", inDirNameMC.data(), "fPIDtree"), inFileNameMC);

    TString pSel = "";
    if(var4proj == kP)
        pSel = "pTPC";
    else
        pSel = "pT";

    double etaLims[101];
    for(int iEtaBin = 0; iEtaBin < 101; iEtaBin++)
        etaLims[iEtaBin] = -1. + 0.02*iEtaBin;

    double nSigmaLims[1001];
    for(int iNsigmaBin = 0; iNsigmaBin < 1001; iNsigmaBin++)
        nSigmaLims[iNsigmaBin] = -50. + 0.1*iNsigmaBin;

    double partLims[6];
    for(int iPartBin = 0; iPartBin < 6; iPartBin++)
        partLims[iPartBin] = iPartBin;
    std::map<int, double> partBins = {{kEl, 0.5}, {kMuon, 1.5}, {kPion, 2.5}, {kKaon, 3.5}, {kPr, 4.5}};

    auto dataFrameMCEta = dataFrameMC.Filter(Form("(eta > %f && eta < %f) || (eta > -%f && eta < -%f)",
                                                 absEtaBinMins[nEtaBins-1]*1000, absEtaBinMaxs[nEtaBins-1]*1000,
                                                 absEtaBinMaxs[nEtaBins-1]*1000, absEtaBinMins[nEtaBins-1]*1000));

    std::cout << "Selecting V0 tagged pions" << std::endl;
    TString tagSel = Form("(((tag & %d) > 0) || ((tag & %d) > 0))", AliAnalysisTaskSEHFSystPID::kIsPionFromK0s, AliAnalysisTaskSEHFSystPID::kIsPionFromL);
    auto dataFrameMCSel = dataFrameMCEta.Filter(tagSel.Data());
    auto hNsigmaTPCPionMCV0tagVsPVsPart = dataFrameMCSel.Define("n_sigma_TPC_pi_scaled", "static_cast<float>(n_sigma_TPC_pi)/100")
                                                        .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                        .Define("part", "if(PDGcode == 11) return 0; else if(PDGcode == 13) return 1; else if(PDGcode == 211) return 2; else if(PDGcode == 321) return 3; else if(PDGcode == 2212) return 4; else return -1;")
                                                        .Histo3D({"hNsigmaTPCPionMCV0tagVsPVsPart", "", 5u, partLims, static_cast<int>(nBins), binLims, 1000u, nSigmaLims}, "part", "p_scaled", "n_sigma_TPC_pi_scaled");

    auto hNsigmaTOFPionMCV0tagVsPVsPart = dataFrameMCSel.Define("n_sigma_TOF_pi_scaled", "static_cast<float>(n_sigma_TOF_pi)/100")
                                                        .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                        .Define("part", "if(PDGcode == 11) return 0; else if(PDGcode == 13) return 1; else if(PDGcode == 211) return 2; else if(PDGcode == 321) return 3; else if(PDGcode == 2212) return 4; else return -1;")
                                                        .Histo3D({"hNsigmaTOFPionMCV0tagVsPVsPart", "", 5u, partLims, static_cast<int>(nBins), binLims, 1000u, nSigmaLims}, "part", "p_scaled", "n_sigma_TOF_pi_scaled");

    std::cout << "Selecting kinks tagged kaons" << std::endl;
    tagSel = Form("((tag & %d) > 0)", AliAnalysisTaskSEHFSystPID::kIsKaonFromKinks);
    dataFrameMCSel = dataFrameMCEta.Filter(tagSel.Data());
    auto hNsigmaTPCKaonMCKinktagVsPVsPart = dataFrameMCSel.Define("n_sigma_TPC_K_scaled", "static_cast<float>(n_sigma_TPC_K)/100")
                                                        .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                        .Define("part", "if(PDGcode == 11) return 0; else if(PDGcode == 13) return 1; else if(PDGcode == 211) return 2; else if(PDGcode == 321) return 3; else if(PDGcode == 2212) return 4; else return -1;")
                                                        .Histo3D({"hNsigmaTPCKaonMCKinktagVsPVsPart", "", 5u, partLims, static_cast<int>(nBins), binLims, 1000u, nSigmaLims}, "part", "p_scaled", "n_sigma_TPC_K_scaled");

    auto hNsigmaTOFKaonMCKinktagVsPVsPart = dataFrameMCSel.Define("n_sigma_TOF_K_scaled", "static_cast<float>(n_sigma_TOF_K)/100")
                                                        .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                        .Define("part", "if(PDGcode == 11) return 0; else if(PDGcode == 13) return 1; else if(PDGcode == 211) return 2; else if(PDGcode == 321) return 3; else if(PDGcode == 2212) return 4; else return -1;")
                                                        .Histo3D({"hNsigmaTOFKaonMCKinktagVsPVsPart", "", 5u, partLims, static_cast<int>(nBins), binLims, 1000u, nSigmaLims}, "part", "p_scaled", "n_sigma_TOF_K_scaled");

    std::cout << "Selecting TOF tagged kaons" << std::endl;
    tagSel = Form("((tag & %d) > 0)", AliAnalysisTaskSEHFSystPID::kIsKaonFromTOF);
    dataFrameMCSel = dataFrameMCEta.Filter(tagSel.Data());
    auto hNsigmaTPCKaonMCTOFtagVsPVsPart = dataFrameMCSel.Define("n_sigma_TPC_K_scaled", "static_cast<float>(n_sigma_TPC_K)/100")
                                                        .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                        .Define("part", "if(PDGcode == 11) return 0; else if(PDGcode == 13) return 1; else if(PDGcode == 211) return 2; else if(PDGcode == 321) return 3; else if(PDGcode == 2212) return 4; else return -1;")
                                                        .Histo3D({"hNsigmaTPCKaonMCTOFtagVsPVsPart", "", 5u, partLims, static_cast<int>(nBins), binLims, 1000u, nSigmaLims}, "part", "p_scaled", "n_sigma_TPC_K_scaled");

    std::cout << "Selecting TPC tagged kaons" << std::endl;
    tagSel = Form("((tag & %d) > 0)", AliAnalysisTaskSEHFSystPID::kIsKaonFromTPC);
    dataFrameMCSel = dataFrameMCEta.Filter(tagSel.Data());
    auto hNsigmaTOFKaonMCTPCtagVsPVsPart = dataFrameMCSel.Define("n_sigma_TOF_K_scaled", "static_cast<float>(n_sigma_TOF_K)/100")
                                                        .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                        .Define("part", "if(PDGcode == 11) return 0; else if(PDGcode == 13) return 1; else if(PDGcode == 211) return 2; else if(PDGcode == 321) return 3; else if(PDGcode == 2212) return 4; else return -1;")
                                                        .Histo3D({"hNsigmaTPCKaonMCTPCtagVsPVsPart", "", 5u, partLims, static_cast<int>(nBins), binLims, 1000u, nSigmaLims}, "part", "p_scaled", "n_sigma_TOF_K_scaled");

    std::cout << "Selecting V0 tagged protons" << std::endl;
    tagSel = Form("((tag & %d) > 0)", AliAnalysisTaskSEHFSystPID::kIsProtonFromL);
    dataFrameMCSel = dataFrameMCEta.Filter(tagSel.Data());
    auto hNsigmaTPCProtonMCV0tagVsPVsPart = dataFrameMCSel.Define("n_sigma_TPC_p_scaled", "static_cast<float>(n_sigma_TPC_p)/100")
                                                        .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                        .Define("part", "if(PDGcode == 11) return 0; else if(PDGcode == 13) return 1; else if(PDGcode == 211) return 2; else if(PDGcode == 321) return 3; else if(PDGcode == 2212) return 4; else return -1;")
                                                        .Histo3D({"hNsigmaTPCProtonMCV0tagVsPVsPart", "", 5u, partLims, static_cast<int>(nBins), binLims, 1000u, nSigmaLims}, "part", "p_scaled", "n_sigma_TPC_p_scaled");

    auto hNsigmaTOFProtonMCV0tagVsPVsPart = dataFrameMCSel.Define("n_sigma_TOF_p_scaled", "static_cast<float>(n_sigma_TOF_p)/100")
                                                        .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                        .Define("part", "if(PDGcode == 11) return 0; else if(PDGcode == 13) return 1; else if(PDGcode == 211) return 2; else if(PDGcode == 321) return 3; else if(PDGcode == 2212) return 4; else return -1;")
                                                        .Histo3D({"hNsigmaTOFProtonMCV0tagVsPVsPart", "", 5u, partLims, static_cast<int>(nBins), binLims, 1000u, nSigmaLims}, "part", "p_scaled", "n_sigma_TOF_p_scaled");

    std::cout << "\n\rProcessing pseudorapidity bin 01/01 (\033[32mMC always only eta integrated\033[0m)";

    for(unsigned int iBin = 0; iBin < nBins; iBin++)
    {
        if(iBin == 0)
            std::cout << Form("\n\rProcessing momentum bin %03d/%03d", iBin+1, nBins);
        else
            std::cout << Form("\rProcessing momentum bin %03d/%03d", iBin+1, nBins);
   
        int pMinBin = hNsigmaTPCPionMCV0tagVsPVsPart->GetYaxis()->FindBin(binMins[iBin]*1.0001);
        int pMaxBin = hNsigmaTPCPionMCV0tagVsPVsPart->GetYaxis()->FindBin(binMaxs[iBin]*0.9999);

        for(auto &part : pdgNames)
        {
            int partBinMin = hNsigmaTPCPionMCV0tagVsPVsPart->GetXaxis()->FindBin(partBins[part.first]);
            int partBinMax = hNsigmaTPCPionMCV0tagVsPVsPart->GetXaxis()->FindBin(partBins[part.first]);
            if(part.first == kAll)
            {
                partBinMin = -1;
                partBinMax = -1;
            }

            hNsigmaTPCPionMCV0tagVsPVsPart->GetXaxis()->SetRange(partBinMin, partBinMax);
            hNsigmaTPCPionMCV0tagVsPVsPart->GetYaxis()->SetRange(pMinBin, pMaxBin);
            hNsigmaTPCPionMCV0tag[iBin][part.first] = static_cast<TH1D*>(hNsigmaTPCPionMCV0tagVsPVsPart->Project3D("z"));
            hNsigmaTPCPionMCV0tag[iBin][part.first]->SetNameTitle(Form("hNsigmaTPCPionMCV0tag_%s_%s_%s", part.second.data(), etaBinLabels[nEtaBins-1].Data(), binLabels[iBin].Data()), Form("%.2f < %s < %.2f GeV/#it{c};N_{#sigma}^{TPC}(#pi);Normalised entries", binMins[iBin], varTitle.Data(), binMaxs[iBin]));
            hNsigmaTPCPionMCV0tagVsPVsPart->GetXaxis()->SetRange(-1, -1);
            hNsigmaTPCPionMCV0tagVsPVsPart->GetYaxis()->SetRange(-1, -1);

            hNsigmaTOFPionMCV0tagVsPVsPart->GetXaxis()->SetRange(partBinMin, partBinMax);
            hNsigmaTOFPionMCV0tagVsPVsPart->GetYaxis()->SetRange(pMinBin, pMaxBin);
            hNsigmaTOFPionMCV0tag[iBin][part.first] = static_cast<TH1D*>(hNsigmaTOFPionMCV0tagVsPVsPart->Project3D("z"));
            hNsigmaTOFPionMCV0tag[iBin][part.first]->SetNameTitle(Form("hNsigmaTOFPionMCV0tag_%s_%s_%s", part.second.data(), etaBinLabels[nEtaBins-1].Data(), binLabels[iBin].Data()), Form("%.2f < %s < %.2f GeV/#it{c};N_{#sigma}^{TOF}(#pi);Normalised entries", binMins[iBin], varTitle.Data(), binMaxs[iBin]));
            hNsigmaTOFPionMCV0tagVsPVsPart->GetXaxis()->SetRange(-1, -1);
            hNsigmaTOFPionMCV0tagVsPVsPart->GetYaxis()->SetRange(-1, -1);

            hNsigmaTPCKaonMCKinktagVsPVsPart->GetXaxis()->SetRange(partBinMin, partBinMax);
            hNsigmaTPCKaonMCKinktagVsPVsPart->GetYaxis()->SetRange(pMinBin, pMaxBin);
            hNsigmaTPCKaonMCKinktag[iBin][part.first] = static_cast<TH1D*>(hNsigmaTPCKaonMCKinktagVsPVsPart->Project3D("z"));
            hNsigmaTPCKaonMCKinktag[iBin][part.first]->SetNameTitle(Form("hNsigmaTPCKaonMCKinktag_%s_%s_%s", part.second.data(), etaBinLabels[nEtaBins-1].Data(), binLabels[iBin].Data()), Form("%.2f < %s < %.2f GeV/#it{c};N_{#sigma}^{TPC}(K);Normalised entries", binMins[iBin], varTitle.Data(), binMaxs[iBin]));
            hNsigmaTPCKaonMCKinktagVsPVsPart->GetXaxis()->SetRange(-1, -1);
            hNsigmaTPCKaonMCKinktagVsPVsPart->GetYaxis()->SetRange(-1, -1);

            hNsigmaTOFKaonMCKinktagVsPVsPart->GetXaxis()->SetRange(partBinMin, partBinMax);
            hNsigmaTOFKaonMCKinktagVsPVsPart->GetYaxis()->SetRange(pMinBin, pMaxBin);
            hNsigmaTOFKaonMCKinktag[iBin][part.first] = static_cast<TH1D*>(hNsigmaTOFKaonMCKinktagVsPVsPart->Project3D("z"));
            hNsigmaTOFKaonMCKinktag[iBin][part.first]->SetNameTitle(Form("hNsigmaTOFKaonMCKinktag_%s_%s_%s", part.second.data(), etaBinLabels[nEtaBins-1].Data(), binLabels[iBin].Data()), Form("%.2f < %s < %.2f GeV/#it{c};N_{#sigma}^{TOF}(K);Normalised entries", binMins[iBin], varTitle.Data(), binMaxs[iBin]));
            hNsigmaTOFKaonMCKinktagVsPVsPart->GetXaxis()->SetRange(-1, -1);
            hNsigmaTOFKaonMCKinktagVsPVsPart->GetYaxis()->SetRange(-1, -1);

            hNsigmaTPCKaonMCTOFtagVsPVsPart->GetXaxis()->SetRange(partBinMin, partBinMax);
            hNsigmaTPCKaonMCTOFtagVsPVsPart->GetYaxis()->SetRange(pMinBin, pMaxBin);
            hNsigmaTPCKaonMCTOFtag[iBin][part.first] = static_cast<TH1D*>(hNsigmaTPCKaonMCTOFtagVsPVsPart->Project3D("z"));
            hNsigmaTPCKaonMCTOFtag[iBin][part.first]->SetNameTitle(Form("hNsigmaTPCKaonMCTOFtag_%s_%s_%s", part.second.data(), etaBinLabels[nEtaBins-1].Data(), binLabels[iBin].Data()), Form("%.2f < %s < %.2f GeV/#it{c};N_{#sigma}^{TPC}(K);Normalised entries", binMins[iBin], varTitle.Data(), binMaxs[iBin]));
            hNsigmaTPCKaonMCTOFtagVsPVsPart->GetXaxis()->SetRange(-1, -1);
            hNsigmaTPCKaonMCTOFtagVsPVsPart->GetYaxis()->SetRange(-1, -1);

            hNsigmaTOFKaonMCTPCtagVsPVsPart->GetXaxis()->SetRange(partBinMin, partBinMax);
            hNsigmaTOFKaonMCTPCtagVsPVsPart->GetYaxis()->SetRange(pMinBin, pMaxBin);
            hNsigmaTOFKaonMCTPCtag[iBin][part.first] = static_cast<TH1D*>(hNsigmaTOFKaonMCTPCtagVsPVsPart->Project3D("z"));
            hNsigmaTOFKaonMCTPCtag[iBin][part.first]->SetNameTitle(Form("hNsigmaTOFKaonMCTPCtag_%s_%s_%s", part.second.data(), etaBinLabels[nEtaBins-1].Data(), binLabels[iBin].Data()), Form("%.2f < %s < %.2f GeV/#it{c};N_{#sigma}^{TOF}(K);Normalised entries", binMins[iBin], varTitle.Data(), binMaxs[iBin]));
            hNsigmaTOFKaonMCTPCtagVsPVsPart->GetXaxis()->SetRange(-1, -1);
            hNsigmaTOFKaonMCTPCtagVsPVsPart->GetYaxis()->SetRange(-1, -1);

            hNsigmaTPCProtonMCV0tagVsPVsPart->GetXaxis()->SetRange(partBinMin, partBinMax);
            hNsigmaTPCProtonMCV0tagVsPVsPart->GetYaxis()->SetRange(pMinBin, pMaxBin);
            hNsigmaTPCProtonMCV0tag[iBin][part.first] = static_cast<TH1D*>(hNsigmaTPCProtonMCV0tagVsPVsPart->Project3D("z"));
            hNsigmaTPCProtonMCV0tag[iBin][part.first]->SetNameTitle(Form("hNsigmaTPCProtonMCV0tag_%s_%s_%s", part.second.data(), etaBinLabels[nEtaBins-1].Data(), binLabels[iBin].Data()), Form("%.2f < %s < %.2f GeV/#it{c};N_{#sigma}^{TPC}(p);Normalised entries", binMins[iBin], varTitle.Data(), binMaxs[iBin]));
            hNsigmaTPCProtonMCV0tagVsPVsPart->GetXaxis()->SetRange(-1, -1);
            hNsigmaTPCProtonMCV0tagVsPVsPart->GetYaxis()->SetRange(-1, -1);

            hNsigmaTOFProtonMCV0tagVsPVsPart->GetXaxis()->SetRange(partBinMin, partBinMax);
            hNsigmaTOFProtonMCV0tagVsPVsPart->GetYaxis()->SetRange(pMinBin, pMaxBin);
            hNsigmaTOFProtonMCV0tag[iBin][part.first] = static_cast<TH1D*>(hNsigmaTOFProtonMCV0tagVsPVsPart->Project3D("z"));
            hNsigmaTOFProtonMCV0tag[iBin][part.first]->SetNameTitle(Form("hNsigmaTOFProtonMCV0tag_%s_%s_%s", part.second.data(), etaBinLabels[nEtaBins-1].Data(), binLabels[iBin].Data()), Form("%.2f < %s < %.2f GeV/#it{c};N_{#sigma}^{TOF}(p);Normalised entries", binMins[iBin], varTitle.Data(), binMaxs[iBin]));
            hNsigmaTOFProtonMCV0tagVsPVsPart->GetXaxis()->SetRange(-1, -1);
            hNsigmaTOFProtonMCV0tagVsPVsPart->GetYaxis()->SetRange(-1, -1);

            SetTH1Style(hNsigmaTPCPionMCV0tag[iBin][part.first], kFullCircle, pdgColors[part.first], 0.6, 2, pdgColors[part.first], pdgFillColors[part.first], 0.055, 0.06);
            SetTH1Style(hNsigmaTPCKaonMCKinktag[iBin][part.first], kFullCircle, pdgColors[part.first], 0.6, 2, pdgColors[part.first], pdgFillColors[part.first], 0.055, 0.06);
            SetTH1Style(hNsigmaTPCKaonMCTOFtag[iBin][part.first], kFullCircle, pdgColors[part.first], 0.6, 2, pdgColors[part.first], pdgFillColors[part.first], 0.055, 0.06);
            SetTH1Style(hNsigmaTPCProtonMCV0tag[iBin][part.first], kFullCircle, pdgColors[part.first], 0.6, 2, pdgColors[part.first], pdgFillColors[part.first], 0.055, 0.06);
            SetTH1Style(hNsigmaTOFPionMCV0tag[iBin][part.first], kFullCircle, pdgColors[part.first], 0.6, 2, pdgColors[part.first], pdgFillColors[part.first], 0.055, 0.06);
            SetTH1Style(hNsigmaTOFKaonMCKinktag[iBin][part.first], kFullCircle, pdgColors[part.first], 0.6, 2, pdgColors[part.first], pdgFillColors[part.first], 0.055, 0.06);
            SetTH1Style(hNsigmaTOFKaonMCTPCtag[iBin][part.first], kFullCircle, pdgColors[part.first], 0.6, 2, pdgColors[part.first], pdgFillColors[part.first], 0.055, 0.06);
            SetTH1Style(hNsigmaTOFProtonMCV0tag[iBin][part.first], kFullCircle, pdgColors[part.first], 0.6, 2, pdgColors[part.first], pdgFillColors[part.first], 0.055, 0.06);
        }
    }

    std::cout << "\n\n\033[32mDone\033[0m" << std::endl;

    //load data inputs
    auto infileData = TFile::Open(inFileNameData.data());
    if (!infileData || !infileData->IsOpen())
        return;
    auto indirData = static_cast<TDirectoryFile *>(infileData->Get(inDirNameData.data()));
    if (!indirData)
    {
        std::cerr << Form("TDirectoryFile %s not found in input file for Data! Exit", inDirNameData.data()) << std::endl;
        return;
    }
    auto listData = static_cast<TList *>(indirData->Get(inListNameData.data()));
    if (!listData)
    {
        std::cerr << Form("TList %s not found in input file for Data! Exit", inListNameData.data()) << std::endl;
        return;
    }
    
    std::cout << "\n*******************************************\n" << std::endl;
    std::cout << "\033[32mProject data tree\033[0m\n" << std::endl;
    ROOT::RDataFrame dataFrameData(Form("%s/%s", inDirNameData.data(), "fPIDtree"), inFileNameData);

    // Nsigma vs. eta and p histos    
    auto hNsigmaTPCKaon = dataFrameData.Define("n_sigma_TPC_K_scaled", "static_cast<float>(n_sigma_TPC_K)/100")
                                       .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                       .Histo2D({"hNsigmaTPCKaon", Form("%.2f < %s < %.2f GeV/#it{c};%s (GeV/#it{c});N_{#sigma}^{TOF}(K)", binMins[0], varTitle.Data(), binMaxs[nBins-1], varTitle.Data()), 100u, binMins[0], binMaxs[nBins-1], 1000u, nSigmaLims}, "p_scaled", "n_sigma_TPC_K_scaled");

    auto hNsigmaTOFKaon = dataFrameData.Define("n_sigma_TOF_K_scaled", "static_cast<float>(n_sigma_TOF_K)/100")
                                       .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                       .Histo2D({"hNsigmaTOFKaon", Form("%.2f < %s < %.2f GeV/#it{c};%s (GeV/#it{c});N_{#sigma}^{TOF}(K)", binMins[0], varTitle.Data(), binMaxs[nBins-1], varTitle.Data()), 100u, binMins[0], binMaxs[nBins-1], 1000u, nSigmaLims}, "p_scaled", "n_sigma_TOF_K_scaled");

    std::cout << "Selecting V0 tagged pions" << std::endl;
    tagSel = Form("(((tag & %d) > 0) || ((tag & %d) > 0))", AliAnalysisTaskSEHFSystPID::kIsPionFromK0s, AliAnalysisTaskSEHFSystPID::kIsPionFromL);
    auto dataFrameDataSel = dataFrameData.Filter(tagSel.Data());
    auto hNsigmaTPCPionDataV0tagVsEtaVsP = dataFrameDataSel.Define("n_sigma_TPC_pi_scaled", "static_cast<float>(n_sigma_TPC_pi)/100")
                                                        .Define("eta_scaled", "static_cast<float>(eta)/1000")
                                                        .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                        .Histo3D({"hNsigmaTPCPionDataV0tagVsEtaVsP", "", static_cast<int>(nBins), binLims, 100u, etaLims, 1000u, nSigmaLims}, "p_scaled", "eta_scaled", "n_sigma_TPC_pi_scaled");
    auto hNsigmaTPCPionDataV0tagVsEta = static_cast<TH2D*>(hNsigmaTPCPionDataV0tagVsEtaVsP->Project3D("zy"));
    hNsigmaTPCPionDataV0tagVsEta->SetNameTitle("hNsigmaTPCPionDataV0tagVsEta", Form("%.2f < %s < %.2f GeV/#it{c};#it{#eta};N_{#sigma}^{TPC}(#pi)", binMins[0], varTitle.Data(), binMaxs[nBins-1]));

    auto hNsigmaTOFPionDataV0tagVsEtaVsP = dataFrameDataSel.Define("n_sigma_TOF_pi_scaled", "static_cast<float>(n_sigma_TOF_pi)/100")
                                                           .Define("eta_scaled", "static_cast<float>(eta)/1000")
                                                           .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                           .Histo3D({"hNsigmaTOFPionDataV0tagVsEtaVsP", "", static_cast<int>(nBins), binLims, 100u, etaLims, 1000u, nSigmaLims}, "p_scaled", "eta_scaled", "n_sigma_TOF_pi_scaled");
    auto hNsigmaTOFPionDataV0tagVsEta = static_cast<TH2D*>(hNsigmaTOFPionDataV0tagVsEtaVsP->Project3D("zy"));
    hNsigmaTOFPionDataV0tagVsEta->SetNameTitle("hNsigmaTOFPionDataV0tagVsEta", Form("%.2f < %s < %.2f GeV/#it{c};#it{#eta};N_{#sigma}^{TOF}(#pi)", binMins[0], varTitle.Data(), binMaxs[nBins-1]));

    std::cout << "Selecting kinks tagged kaons" << std::endl;
    tagSel = Form("((tag & %d) > 0)", AliAnalysisTaskSEHFSystPID::kIsKaonFromKinks);
    dataFrameDataSel = dataFrameData.Filter(tagSel.Data());
    auto hNsigmaTPCKaonDataKinktagVsEtaVsP = dataFrameDataSel.Define("n_sigma_TPC_K_scaled", "static_cast<float>(n_sigma_TPC_K)/100")
                                                             .Define("eta_scaled", "static_cast<float>(eta)/1000")
                                                             .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                             .Histo3D({"hNsigmaTPCKaonDataKinktagVsEtaVsP", "", static_cast<int>(nBins), binLims, 100u, etaLims, 1000u, nSigmaLims}, "p_scaled", "eta_scaled", "n_sigma_TPC_K_scaled");
    auto hNsigmaTPCKaonDataKinktagVsEta = static_cast<TH2D*>(hNsigmaTOFPionDataV0tagVsEtaVsP->Project3D("zy"));
    hNsigmaTPCKaonDataKinktagVsEta->SetNameTitle("hNsigmaTPCKaonDataKinktagVsEta", Form("%.2f < %s < %.2f GeV/#it{c};#it{#eta};N_{#sigma}^{TPC}(K)", binMins[0], varTitle.Data(), binMaxs[nBins-1]));

    auto hNsigmaTOFKaonDataKinktagVsEtaVsP = dataFrameDataSel.Define("n_sigma_TOF_K_scaled", "static_cast<float>(n_sigma_TOF_K)/100")
                                                             .Define("eta_scaled", "static_cast<float>(eta)/1000")
                                                             .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                             .Histo3D({"hNsigmaTOFKaonDataKinktagVsEtaVsP", "", static_cast<int>(nBins), binLims, 100u, etaLims, 1000u, nSigmaLims}, "p_scaled", "eta_scaled", "n_sigma_TOF_K_scaled");
    auto hNsigmaTOFKaonDataKinktagVsEta = static_cast<TH2D*>(hNsigmaTOFKaonDataKinktagVsEtaVsP->Project3D("zy"));
    hNsigmaTOFKaonDataKinktagVsEta->SetNameTitle("hNsigmaTOFKaonDataKinktagVsEta", Form("%.2f < %s < %.2f GeV/#it{c};#it{#eta};N_{#sigma}^{TOF}(K)", binMins[0], varTitle.Data(), binMaxs[nBins-1]));

    std::cout << "Selecting TOF tagged kaons" << std::endl;
    tagSel = Form("((tag & %d) > 0)", AliAnalysisTaskSEHFSystPID::kIsKaonFromTOF);
    dataFrameDataSel = dataFrameData.Filter(tagSel.Data());
    auto hNsigmaTPCKaonDataTOFtagVsEtaVsP = dataFrameDataSel.Define("n_sigma_TPC_K_scaled", "static_cast<float>(n_sigma_TPC_K)/100")
                                                            .Define("eta_scaled", "static_cast<float>(eta)/1000")
                                                            .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                            .Histo3D({"hNsigmaTPCKaonDataTOFtagVsEtaVsP", "", static_cast<int>(nBins), binLims, 100u, etaLims, 1000u, nSigmaLims}, "p_scaled", "eta_scaled", "n_sigma_TPC_K_scaled");
    auto hNsigmaTPCKaonDataTOFtagVsEta = static_cast<TH2D*>(hNsigmaTPCKaonDataTOFtagVsEtaVsP->Project3D("zy"));
    hNsigmaTPCKaonDataTOFtagVsEta->SetNameTitle("hNsigmaTPCKaonDataTOFtagVsEta", Form("%.2f < %s < %.2f GeV/#it{c};#it{#eta};N_{#sigma}^{TPC}(K)", binMins[0], varTitle.Data(), binMaxs[nBins-1]));

    auto hNsigmaTOFKaonTOFtagged = dataFrameDataSel.Define("n_sigma_TOF_K_scaled", "static_cast<float>(n_sigma_TOF_K)/100")
                                                   .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                   .Histo2D({"hNsigmaTOFKaonTOFtagged", Form("%.2f < %s < %.2f GeV/#it{c};%s (GeV/#it{c});N_{#sigma}^{TOF}(K)", binMins[0], varTitle.Data(), binMaxs[nBins-1], varTitle.Data()), 100u, binMins[0], binMaxs[nBins-1], 1000u, nSigmaLims}, "p_scaled", "n_sigma_TOF_K_scaled");

    std::cout << "Selecting TPC tagged kaons" << std::endl;
    tagSel = Form("((tag & %d) > 0)", AliAnalysisTaskSEHFSystPID::kIsKaonFromTPC);
    dataFrameDataSel = dataFrameData.Filter(tagSel.Data());
    auto hNsigmaTOFKaonDataTPCtagVsEtaVsP = dataFrameDataSel.Define("n_sigma_TOF_K_scaled", "static_cast<float>(n_sigma_TOF_K)/100")
                                                            .Define("eta_scaled", "static_cast<float>(eta)/1000")
                                                            .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                            .Histo3D({"hNsigmaTOFKaonDataTPCtagVsEtaVsP", "", static_cast<int>(nBins), binLims, 100u, etaLims, 1000u, nSigmaLims}, "p_scaled", "eta_scaled", "n_sigma_TOF_K_scaled");
    auto hNsigmaTOFKaonDataTPCtagVsEta = static_cast<TH2D*>(hNsigmaTOFKaonDataTPCtagVsEtaVsP->Project3D("zy"));
    hNsigmaTOFKaonDataTPCtagVsEta->SetNameTitle("hNsigmaTOFKaonDataTPCtagVsEta", Form("%.2f < %s < %.2f GeV/#it{c};#it{#eta};N_{#sigma}^{TOF}(K)", binMins[0], varTitle.Data(), binMaxs[nBins-1]));

    auto hNsigmaTPCKaonTPCtagged = dataFrameDataSel.Define("n_sigma_TPC_K_scaled", "static_cast<float>(n_sigma_TPC_K)/100")
                                                   .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                   .Histo2D({"hNsigmaTPCKaonTPCtagged", Form("%.2f < %s < %.2f GeV/#it{c};%s (GeV/#it{c});N_{#sigma}^{TOF}(K)", binMins[0], varTitle.Data(), binMaxs[nBins-1], varTitle.Data()), 100u, binMins[0], binMaxs[nBins-1], 1000u, nSigmaLims}, "p_scaled", "n_sigma_TPC_K_scaled");

    std::cout << "Selecting V0 tagged protons\n" << std::endl;
    tagSel = Form("((tag & %d) > 0)", AliAnalysisTaskSEHFSystPID::kIsProtonFromL);
    dataFrameDataSel = dataFrameData.Filter(tagSel.Data());
    auto hNsigmaTPCProtonDataV0tagVsEtaVsP = dataFrameDataSel.Define("n_sigma_TPC_p_scaled", "static_cast<float>(n_sigma_TPC_p)/100")
                                                             .Define("eta_scaled", "static_cast<float>(eta)/1000")
                                                             .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                             .Histo3D({"hNsigmaTPCProtonDataV0tagVsEtaVsP", "", static_cast<int>(nBins), binLims, 100u, etaLims, 1000u, nSigmaLims}, "p_scaled", "eta_scaled", "n_sigma_TPC_p_scaled");
    auto hNsigmaTPCProtonDataV0tagVsEta = static_cast<TH2D*>(hNsigmaTPCProtonDataV0tagVsEtaVsP->Project3D("zy"));
    hNsigmaTPCProtonDataV0tagVsEta->SetNameTitle("hNsigmaTPCProtonDataV0tagVsEta", Form("%.2f < %s < %.2f GeV/#it{c};#it{#eta};N_{#sigma}^{TPC}(p)", binMins[0], varTitle.Data(), binMaxs[nBins-1]));

    auto hNsigmaTOFProtonDataV0tagVsEtaVsP = dataFrameDataSel.Define("n_sigma_TOF_p_scaled", "static_cast<float>(n_sigma_TOF_p)/100")
                                                             .Define("eta_scaled", "static_cast<float>(eta)/1000")
                                                             .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                             .Histo3D({"hNsigmaTOFProtonDataV0tagVsEtaVsP", "", static_cast<int>(nBins), binLims, 100u, etaLims, 1000u, nSigmaLims}, "p_scaled", "eta_scaled", "n_sigma_TOF_p_scaled");
    auto hNsigmaTOFProtonDataV0tagVsEta = static_cast<TH2D*>(hNsigmaTOFProtonDataV0tagVsEtaVsP->Project3D("zy"));
    hNsigmaTOFProtonDataV0tagVsEta->SetNameTitle("hNsigmaTOFProtonDataV0tagVsEta", Form("%.2f < %s < %.2f GeV/#it{c};#it{#eta};N_{#sigma}^{TOF}(p)", binMins[0], varTitle.Data(), binMaxs[nBins-1]));

    for(unsigned int iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++)
    {
        std::cout << Form("\rProcessing pseudorapidity bin %02d/%02d", iEtaBin+1, nEtaBins);
        int etaMinBinPos = hNsigmaTPCPionDataV0tagVsEtaVsP->GetYaxis()->FindBin(absEtaBinMins[iEtaBin]*1.0001);
        int etaMaxBinPos = hNsigmaTPCPionDataV0tagVsEtaVsP->GetYaxis()->FindBin(absEtaBinMaxs[iEtaBin]*0.9999);
        int etaMinBinNeg = hNsigmaTPCPionDataV0tagVsEtaVsP->GetYaxis()->FindBin(-absEtaBinMaxs[iEtaBin]*0.9999);
        int etaMaxBinNeg = hNsigmaTPCPionDataV0tagVsEtaVsP->GetYaxis()->FindBin(-absEtaBinMins[iEtaBin]*1.0001);

        for(unsigned int iBin = 0; iBin < nBins; iBin++)
        {
            if(iBin == 0)
                std::cout << Form("\n\rProcessing momentum bin %03d/%03d", iBin+1, nBins);
            else
                std::cout << Form("\rProcessing momentum bin %03d/%03d", iBin+1, nBins);
            
            int pMinBin = hNsigmaTPCPionDataV0tagVsEtaVsP->GetXaxis()->FindBin(binMins[iBin]*1.0001);
            int pMaxBin = hNsigmaTPCPionDataV0tagVsEtaVsP->GetXaxis()->FindBin(binMaxs[iBin]*0.9999);

            hNsigmaTPCPionDataV0tagVsEtaVsP->GetXaxis()->SetRange(pMinBin, pMaxBin);
            hNsigmaTPCPionDataV0tagVsEtaVsP->GetYaxis()->SetRange(etaMinBinPos, etaMaxBinPos);
            hNsigmaTPCPionDataV0tag[iEtaBin][iBin] = static_cast<TH1D*>(hNsigmaTPCPionDataV0tagVsEtaVsP->Project3D("z"));
            hNsigmaTPCPionDataV0tagVsEtaVsP->GetYaxis()->SetRange(etaMinBinNeg, etaMaxBinNeg);
            hNsigmaTPCPionDataV0tag[iEtaBin][iBin]->Add(static_cast<TH1D*>(hNsigmaTPCPionDataV0tagVsEtaVsP->Project3D("z")));
            hNsigmaTPCPionDataV0tag[iEtaBin][iBin]->SetNameTitle(Form("hNsigmaTPCPionDataV0tag_%s_%s", etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()), Form("%.2f < %s < %.2f GeV/#it{c};N_{#sigma}^{TPC}(#pi);Normalised entries", binMins[iBin], varTitle.Data(), binMaxs[iBin]));
            hNsigmaTPCPionDataV0tagVsEtaVsP->GetXaxis()->SetRange(-1, -1);
            hNsigmaTPCPionDataV0tagVsEtaVsP->GetYaxis()->SetRange(-1, -1);

            hNsigmaTOFPionDataV0tagVsEtaVsP->GetXaxis()->SetRange(pMinBin, pMaxBin);
            hNsigmaTOFPionDataV0tagVsEtaVsP->GetYaxis()->SetRange(etaMinBinPos, etaMaxBinPos);
            hNsigmaTOFPionDataV0tag[iEtaBin][iBin] = static_cast<TH1D*>(hNsigmaTOFPionDataV0tagVsEtaVsP->Project3D("z"));
            hNsigmaTOFPionDataV0tagVsEtaVsP->GetYaxis()->SetRange(etaMinBinNeg, etaMaxBinNeg);
            hNsigmaTOFPionDataV0tag[iEtaBin][iBin]->Add(static_cast<TH1D*>(hNsigmaTOFPionDataV0tagVsEtaVsP->Project3D("z")));
            hNsigmaTOFPionDataV0tag[iEtaBin][iBin]->SetNameTitle(Form("hNsigmaTOFPionDataV0tag_%s_%s", etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()), Form("%.2f < %s < %.2f GeV/#it{c};N_{#sigma}^{TOF}(#pi);Normalised entries", binMins[iBin], varTitle.Data(), binMaxs[iBin]));
            hNsigmaTOFPionDataV0tagVsEtaVsP->GetXaxis()->SetRange(-1, -1);
            hNsigmaTOFPionDataV0tagVsEtaVsP->GetYaxis()->SetRange(-1, -1);

            hNsigmaTPCKaonDataKinktagVsEtaVsP->GetXaxis()->SetRange(pMinBin, pMaxBin);
            hNsigmaTPCKaonDataKinktagVsEtaVsP->GetYaxis()->SetRange(etaMinBinPos, etaMaxBinPos);
            hNsigmaTPCKaonDataKinktag[iEtaBin][iBin] = static_cast<TH1D*>(hNsigmaTPCKaonDataKinktagVsEtaVsP->Project3D("z"));
            hNsigmaTPCKaonDataKinktagVsEtaVsP->GetYaxis()->SetRange(etaMinBinNeg, etaMaxBinNeg);
            hNsigmaTPCKaonDataKinktag[iEtaBin][iBin]->Add(static_cast<TH1D*>(hNsigmaTPCKaonDataKinktagVsEtaVsP->Project3D("z")));
            hNsigmaTPCKaonDataKinktag[iEtaBin][iBin]->SetNameTitle(Form("hNsigmaTPCKaonDataKinktag_%s_%s", etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()), Form("%.2f < %s < %.2f GeV/#it{c};N_{#sigma}^{TPC}(K);Normalised entries", binMins[iBin], varTitle.Data(), binMaxs[iBin]));
            hNsigmaTPCKaonDataKinktagVsEtaVsP->GetXaxis()->SetRange(-1, -1);
            hNsigmaTPCKaonDataKinktagVsEtaVsP->GetYaxis()->SetRange(-1, -1);

            hNsigmaTOFKaonDataKinktagVsEtaVsP->GetXaxis()->SetRange(pMinBin, pMaxBin);
            hNsigmaTOFKaonDataKinktagVsEtaVsP->GetYaxis()->SetRange(etaMinBinPos, etaMaxBinPos);
            hNsigmaTOFKaonDataKinktag[iEtaBin][iBin] = static_cast<TH1D*>(hNsigmaTOFKaonDataKinktagVsEtaVsP->Project3D("z"));
            hNsigmaTOFKaonDataKinktagVsEtaVsP->GetYaxis()->SetRange(etaMinBinNeg, etaMaxBinNeg);
            hNsigmaTOFKaonDataKinktag[iEtaBin][iBin]->Add(static_cast<TH1D*>(hNsigmaTOFKaonDataKinktagVsEtaVsP->Project3D("z")));
            hNsigmaTOFKaonDataKinktag[iEtaBin][iBin]->SetNameTitle(Form("hNsigmaTOFKaonDataKinktag_%s_%s", etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()), Form("%.2f < %s < %.2f GeV/#it{c};N_{#sigma}^{TOF}(K);Normalised entries", binMins[iBin], varTitle.Data(), binMaxs[iBin]));
            hNsigmaTOFKaonDataKinktagVsEtaVsP->GetXaxis()->SetRange(-1, -1);
            hNsigmaTOFKaonDataKinktagVsEtaVsP->GetYaxis()->SetRange(-1, -1);

            hNsigmaTPCKaonDataTOFtagVsEtaVsP->GetXaxis()->SetRange(pMinBin, pMaxBin);
            hNsigmaTPCKaonDataTOFtagVsEtaVsP->GetYaxis()->SetRange(etaMinBinPos, etaMaxBinPos);
            hNsigmaTPCKaonDataTOFtag[iEtaBin][iBin] = static_cast<TH1D*>(hNsigmaTPCKaonDataTOFtagVsEtaVsP->Project3D("z"));
            hNsigmaTPCKaonDataTOFtagVsEtaVsP->GetYaxis()->SetRange(etaMinBinNeg, etaMaxBinNeg);
            hNsigmaTPCKaonDataTOFtag[iEtaBin][iBin]->Add(static_cast<TH1D*>(hNsigmaTPCKaonDataTOFtagVsEtaVsP->Project3D("z")));
            hNsigmaTPCKaonDataTOFtag[iEtaBin][iBin]->SetNameTitle(Form("hNsigmaTPCKaonDataTOFtag_%s_%s", etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()), Form("%.2f < %s < %.2f GeV/#it{c};N_{#sigma}^{TPC}(K);Normalised entries", binMins[iBin], varTitle.Data(), binMaxs[iBin]));
            hNsigmaTPCKaonDataTOFtagVsEtaVsP->GetXaxis()->SetRange(-1, -1);
            hNsigmaTPCKaonDataTOFtagVsEtaVsP->GetYaxis()->SetRange(-1, -1);

            hNsigmaTOFKaonDataTPCtagVsEtaVsP->GetXaxis()->SetRange(pMinBin, pMaxBin);
            hNsigmaTOFKaonDataTPCtagVsEtaVsP->GetYaxis()->SetRange(etaMinBinPos, etaMaxBinPos);
            hNsigmaTOFKaonDataTPCtag[iEtaBin][iBin] = static_cast<TH1D*>(hNsigmaTOFKaonDataTPCtagVsEtaVsP->Project3D("z"));
            hNsigmaTOFKaonDataTPCtagVsEtaVsP->GetYaxis()->SetRange(etaMinBinNeg, etaMaxBinNeg);
            hNsigmaTOFKaonDataTPCtag[iEtaBin][iBin]->Add(static_cast<TH1D*>(hNsigmaTOFKaonDataTPCtagVsEtaVsP->Project3D("z")));
            hNsigmaTOFKaonDataTPCtag[iEtaBin][iBin]->SetNameTitle(Form("hNsigmaTOFKaonDataTPCtag_%s_%s", etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()), Form("%.2f < %s < %.2f GeV/#it{c};N_{#sigma}^{TOF}(K);Normalised entries", binMins[iBin], varTitle.Data(), binMaxs[iBin]));
            hNsigmaTOFKaonDataTPCtagVsEtaVsP->GetXaxis()->SetRange(-1, -1);
            hNsigmaTOFKaonDataTPCtagVsEtaVsP->GetYaxis()->SetRange(-1, -1);

            hNsigmaTPCProtonDataV0tagVsEtaVsP->GetXaxis()->SetRange(pMinBin, pMaxBin);
            hNsigmaTPCProtonDataV0tagVsEtaVsP->GetYaxis()->SetRange(etaMinBinPos, etaMaxBinPos);
            hNsigmaTPCProtonDataV0tag[iEtaBin][iBin] = static_cast<TH1D*>(hNsigmaTPCProtonDataV0tagVsEtaVsP->Project3D("z"));
            hNsigmaTPCProtonDataV0tagVsEtaVsP->GetYaxis()->SetRange(etaMinBinNeg, etaMaxBinNeg);
            hNsigmaTPCProtonDataV0tag[iEtaBin][iBin]->Add(static_cast<TH1D*>(hNsigmaTPCProtonDataV0tagVsEtaVsP->Project3D("z")));
            hNsigmaTPCProtonDataV0tag[iEtaBin][iBin]->SetNameTitle(Form("hNsigmaTPCProtonDataV0tag_%s_%s", etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()), Form("%.2f < %s < %.2f GeV/#it{c};N_{#sigma}^{TPC}(p);Normalised entries", binMins[iBin], varTitle.Data(), binMaxs[iBin]));
            hNsigmaTPCProtonDataV0tagVsEtaVsP->GetXaxis()->SetRange(-1, -1);
            hNsigmaTPCProtonDataV0tagVsEtaVsP->GetYaxis()->SetRange(-1, -1);

            hNsigmaTOFProtonDataV0tagVsEtaVsP->GetXaxis()->SetRange(pMinBin, pMaxBin);
            hNsigmaTOFProtonDataV0tagVsEtaVsP->GetYaxis()->SetRange(etaMinBinPos, etaMaxBinPos);
            hNsigmaTOFProtonDataV0tag[iEtaBin][iBin] = static_cast<TH1D*>(hNsigmaTOFProtonDataV0tagVsEtaVsP->Project3D("z"));
            hNsigmaTOFProtonDataV0tagVsEtaVsP->GetYaxis()->SetRange(etaMinBinNeg, etaMaxBinNeg);
            hNsigmaTOFProtonDataV0tag[iEtaBin][iBin]->Add(static_cast<TH1D*>(hNsigmaTOFProtonDataV0tagVsEtaVsP->Project3D("z")));
            hNsigmaTOFProtonDataV0tag[iEtaBin][iBin]->SetNameTitle(Form("hNsigmaTOFProtonDataV0tag_%s_%s", etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()), Form("%.2f < %s < %.2f GeV/#it{c};N_{#sigma}^{TPC}(p);Normalised entries", binMins[iBin], varTitle.Data(), binMaxs[iBin]));
            hNsigmaTOFProtonDataV0tagVsEtaVsP->GetXaxis()->SetRange(-1, -1);
            hNsigmaTOFProtonDataV0tagVsEtaVsP->GetYaxis()->SetRange(-1, -1);

            SetTH1Style(hNsigmaTPCPionDataV0tag[iEtaBin][iBin], kFullCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], kWhite, 0.055, 0.06);
            SetTH1Style(hNsigmaTPCKaonDataKinktag[iEtaBin][iBin], kFullCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], kWhite, 0.055, 0.06);
            SetTH1Style(hNsigmaTPCKaonDataTOFtag[iEtaBin][iBin], kFullCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], kWhite, 0.055, 0.06);
            SetTH1Style(hNsigmaTPCProtonDataV0tag[iEtaBin][iBin], kFullCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], kWhite, 0.055, 0.06);
            SetTH1Style(hNsigmaTOFPionDataV0tag[iEtaBin][iBin], kFullCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], kWhite, 0.055, 0.06);
            SetTH1Style(hNsigmaTOFKaonDataKinktag[iEtaBin][iBin], kFullCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], kWhite, 0.055, 0.06);
            SetTH1Style(hNsigmaTOFKaonDataTPCtag[iEtaBin][iBin], kFullCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], kWhite, 0.055, 0.06);
            SetTH1Style(hNsigmaTOFProtonDataV0tag[iEtaBin][iBin], kFullCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], kWhite, 0.055, 0.06);
        }
    }

    std::cout << "\n\n\033[32mDone\033[0m" << std::endl;

    TExec *exGrey = new TExec("exGrey","gStyle->SetPalette(kGreyScale);");
    TExec *exBlue = new TExec("exBlue","gStyle->SetPalette(kDeepSea);");
    TExec *exStd = new TExec("exStd","gStyle->SetPalette(kRainBow);");

    if(produceQAplots)
    {
        //compute fractions in data (only TOF)
        std::cout << "\n*******************************************\n" << std::endl;
        std::cout << "\033[32mProduce QA plots\033[0m\n" << std::endl;
        
        PlotQAhistos(listMC, listData, outDirName);

        // Produce plots for 2D histograms

        TLatex *latAll = new TLatex();
        latAll->SetNDC();
        latAll->SetTextFont(42);
        latAll->SetTextColor(kBlack);
        latAll->SetTextSize(0.045);
        TLatex *latTag = new TLatex();
        latTag->SetNDC();
        latTag->SetTextFont(42);
        latTag->SetTextColor(kAzure+4);
        latTag->SetTextSize(0.045);

        TCanvas *cNsigmaTag = new TCanvas("cNsigmaTag", "", 1920, 1080);
        cNsigmaTag->Divide(2, 1);
        cNsigmaTag->cd(1)->SetLogz();
        hNsigmaTPCKaon->GetZaxis()->SetRangeUser(1., hNsigmaTPCKaon->GetMaximum());
        hNsigmaTPCKaon->DrawCopy("axis");
        exGrey->Draw();
        hNsigmaTPCKaon->DrawCopy("col same");
        exBlue->Draw();
        hNsigmaTPCKaonTPCtagged->DrawCopy("col same");
        latAll->DrawLatex(0.65, 0.25, "All tags");
        latTag->DrawLatex(0.65, 0.2, "TPC tag");
        gPad->Update();

        cNsigmaTag->cd(2)->SetLogz();
        hNsigmaTOFKaon->GetZaxis()->SetRangeUser(1., hNsigmaTOFKaon->GetMaximum());
        hNsigmaTOFKaon->DrawCopy("axis");
        exGrey->Draw();
        hNsigmaTOFKaon->DrawCopy("col same");
        exBlue->Draw();
        hNsigmaTOFKaonTOFtagged->DrawCopy("col same");
        latAll->DrawLatex(0.65, 0.25, "All tags");
        latTag->DrawLatex(0.65, 0.2, "TOF tag");
        gPad->Update();

        cNsigmaTag->SaveAs(Form("%s/DistrNsigmaTPC_TOF_tag.pdf", outDirName.data()));

        std::cout << "\033[32mDone\033[0m" << std::endl;
    }

    int nSigmaBinMin = hNsigmaTPCPionDataV0tagVsEta->GetYaxis()->FindBin(-5.*0.9999);
    int nSigmaBinMax = hNsigmaTPCPionDataV0tagVsEta->GetYaxis()->FindBin(5.*0.9999);
    
    TProfile *hProfileNsigmaTPCPionDataV0tagVsEta = hNsigmaTPCPionDataV0tagVsEta->ProfileX("hProfileNsigmaTPCPionDataV0tagVsEta", nSigmaBinMin, nSigmaBinMax);
    SetTH1Style(hProfileNsigmaTPCPionDataV0tagVsEta, 0, 0, 0., 2, kBlack, kWhite);
    TProfile *hProfileNsigmaTPCProtonDataV0tagVsEta = hNsigmaTPCProtonDataV0tagVsEta->ProfileX("hProfileNsigmaTPCProtonDataV0tagVsEta", nSigmaBinMin, nSigmaBinMax);
    SetTH1Style(hProfileNsigmaTPCProtonDataV0tagVsEta, 0, 0, 0., 2, kBlack, kWhite);
    TProfile *hProfileNsigmaTPCKaonDataTOFtagVsEta = hNsigmaTPCKaonDataTOFtagVsEta->ProfileX("hProfileNsigmaTPCKaonDataTOFtagVsEta", nSigmaBinMin, nSigmaBinMax);
    SetTH1Style(hProfileNsigmaTPCKaonDataTOFtagVsEta, 0, 0, 0., 2, kBlack, kWhite);
    TProfile *hProfileNsigmaTPCKaonDataKinktagVsEta = hNsigmaTPCKaonDataKinktagVsEta->ProfileX("hProfileNsigmaTPCKaonDataKinktagVsEta", nSigmaBinMin, nSigmaBinMax);
    SetTH1Style(hProfileNsigmaTPCKaonDataKinktagVsEta, 0, 0, 0., 2, kBlack, kWhite);

    TProfile *hProfileNsigmaTOFPionDataV0tagVsEta = hNsigmaTOFPionDataV0tagVsEta->ProfileX("hProfileNsigmaTOFPionDataV0tagVsEta", nSigmaBinMin, nSigmaBinMax);
    SetTH1Style(hProfileNsigmaTOFPionDataV0tagVsEta, 0, 0, 0., 2, kBlack, kWhite);
    TProfile *hProfileNsigmaTOFProtonDataV0tagVsEta = hNsigmaTOFProtonDataV0tagVsEta->ProfileX("hProfileNsigmaTOFProtonDataV0tagVsEta", nSigmaBinMin, nSigmaBinMax);
    SetTH1Style(hProfileNsigmaTOFProtonDataV0tagVsEta, 0, 0, 0., 2, kBlack, kWhite);
    TProfile *hProfileNsigmaTOFKaonDataTPCtagVsEta = hNsigmaTOFKaonDataTPCtagVsEta->ProfileX("hProfileNsigmaTOFKaonDataTPCtagVsEta", nSigmaBinMin, nSigmaBinMax);
    SetTH1Style(hProfileNsigmaTOFKaonDataTPCtagVsEta, 0, 0, 0., 2, kBlack, kWhite);
    TProfile *hProfileNsigmaTOFKaonDataKinktagVsEta = hNsigmaTOFKaonDataKinktagVsEta->ProfileX("hProfileNsigmaTOFKaonDataKinktagVsEta", nSigmaBinMin, nSigmaBinMax);
    SetTH1Style(hProfileNsigmaTOFKaonDataKinktagVsEta, 0, 0, 0., 2, kBlack, kWhite);

    TCanvas *cNsigmaTPCVsEta = new TCanvas("cNsigmaTPCVsEta", "", 1920, 1080);
    cNsigmaTPCVsEta->Divide(2, 2);
    cNsigmaTPCVsEta->cd(1)->SetLogz();
    hNsigmaTPCPionDataV0tagVsEta->GetYaxis()->SetRangeUser(-5., 5.);
    hNsigmaTPCPionDataV0tagVsEta->DrawCopy("axis");
    exStd->Draw();
    hNsigmaTPCPionDataV0tagVsEta->DrawCopy("colz same");
    hProfileNsigmaTPCPionDataV0tagVsEta->DrawCopy("same");
    cNsigmaTPCVsEta->cd(2)->SetLogz();
    hNsigmaTPCProtonDataV0tagVsEta->GetYaxis()->SetRangeUser(-5., 5.);
    hNsigmaTPCProtonDataV0tagVsEta->DrawCopy("axis");
    exStd->Draw();
    hNsigmaTPCProtonDataV0tagVsEta->DrawCopy("colz same");
    hProfileNsigmaTPCProtonDataV0tagVsEta->DrawCopy("same");
    cNsigmaTPCVsEta->cd(3)->SetLogz();
    hNsigmaTPCKaonDataTOFtagVsEta->GetYaxis()->SetRangeUser(-5., 5.);
    hNsigmaTPCKaonDataTOFtagVsEta->DrawCopy("axis");
    exStd->Draw();
    hNsigmaTPCKaonDataTOFtagVsEta->DrawCopy("colz same");
    hProfileNsigmaTPCKaonDataTOFtagVsEta->DrawCopy("same");
    cNsigmaTPCVsEta->cd(4)->SetLogz();
    hNsigmaTPCKaonDataKinktagVsEta->GetYaxis()->SetRangeUser(-5., 5.);
    hNsigmaTPCKaonDataKinktagVsEta->DrawCopy("axis");
    exStd->Draw();
    hNsigmaTPCKaonDataKinktagVsEta->DrawCopy("colz same");
    hProfileNsigmaTPCKaonDataKinktagVsEta->DrawCopy("same");
    gPad->Update();

    TCanvas *cNsigmaTOFVsEta = new TCanvas("cNsigmaTOFVsEta", "", 1920, 1080);
    cNsigmaTOFVsEta->Divide(2, 2);
    cNsigmaTOFVsEta->cd(1)->SetLogz();
    hNsigmaTOFPionDataV0tagVsEta->GetYaxis()->SetRangeUser(-5., 5.);
    hNsigmaTOFPionDataV0tagVsEta->DrawCopy("axis");
    exStd->Draw();
    hNsigmaTOFPionDataV0tagVsEta->DrawCopy("colz same");
    hProfileNsigmaTOFPionDataV0tagVsEta->DrawCopy("same");
    cNsigmaTOFVsEta->cd(2)->SetLogz();
    hNsigmaTOFProtonDataV0tagVsEta->GetYaxis()->SetRangeUser(-5., 5.);
    hNsigmaTOFProtonDataV0tagVsEta->DrawCopy("axis");
    exStd->Draw();
    hNsigmaTOFProtonDataV0tagVsEta->DrawCopy("colz same");
    hProfileNsigmaTOFProtonDataV0tagVsEta->DrawCopy("same");
    cNsigmaTOFVsEta->cd(3)->SetLogz();
    hNsigmaTOFKaonDataTPCtagVsEta->GetYaxis()->SetRangeUser(-5., 5.);
    hNsigmaTOFKaonDataTPCtagVsEta->DrawCopy("axis");
    exStd->Draw();
    hNsigmaTOFKaonDataTPCtagVsEta->DrawCopy("colz same");
    hProfileNsigmaTOFKaonDataTPCtagVsEta->DrawCopy("same");
    cNsigmaTOFVsEta->cd(4)->SetLogz();
    hNsigmaTOFKaonDataKinktagVsEta->GetYaxis()->SetRangeUser(-5., 5.);
    hNsigmaTOFKaonDataKinktagVsEta->DrawCopy("axis");
    exStd->Draw();
    hNsigmaTOFKaonDataKinktagVsEta->DrawCopy("colz same");
    hProfileNsigmaTOFKaonDataKinktagVsEta->DrawCopy("same");
    gPad->Update();

    cNsigmaTPCVsEta->SaveAs(Form("%s/NsigmaTPCvsEta.pdf", outDirName.data()));
    cNsigmaTOFVsEta->SaveAs(Form("%s/NsigmaTOFvsEta.pdf", outDirName.data()));

    //compute Fraction and contamination in MC
    std::cout << "\n*******************************************\n" << std::endl;
    std::cout << "\033[32mCompute fractions and contaminations in MC\033[0m\n" << std::endl;

    std::array<std::map<int, double>, nBinsMax> intNsTPCPionMCV0tag, intNsTPCKaonMCKinktag, intNsTPCKaonMCTOFtag, intNsTPCProtonMCV0tag;
    std::array<std::map<int, double>, nBinsMax> intNsTOFPionMCV0tag, intNsTOFKaonMCKinktag, intNsTOFKaonMCTPCtag, intNsTOFProtonMCV0tag;

    std::array<std::array<double, nBinsMax>, nEtaBinsMax+1> intNsTPCPionDataV0tag, intNsTPCKaonDataKinktag, intNsTPCKaonDataTOFtag, intNsTPCProtonDataV0tag;
    std::array<std::array<double, nBinsMax>, nEtaBinsMax+1> intNsTOFPionDataV0tag, intNsTOFKaonDataKinktag, intNsTOFKaonDataTPCtag, intNsTOFProtonDataV0tag;

    std::map<int, TH1D*> hFracTPCPionMCV0tag, hFracTPCKaonMCKinktag, hFracTPCKaonMCTOFtag, hFracTPCProtonMCV0tag;
    std::map<int, TH1D*> hFracTOFPionMCV0tag, hFracTOFKaonMCKinktag, hFracTOFKaonMCTPCtag, hFracTOFProtonMCV0tag;

    TLegend *legFracMC = new TLegend(0.2, 0.4, 0.45, 0.7);
    legFracMC->SetTextSize(0.05);

    TCanvas* cFractionTPCMC = new TCanvas(Form("cFractionTPCMC_%s", etaBinLabels[nEtaBins-1].Data()), Form("cFractionTPCMC_%s", etaBinLabels[nEtaBins-1].Data()), 1920, 1080);
    cFractionTPCMC->Divide(2, 2);
    TCanvas* cFractionTOFMC = new TCanvas(Form("cFractionTOFMC_%s", etaBinLabels[nEtaBins-1].Data()), Form("cFractionTOFMC_%s", etaBinLabels[nEtaBins-1].Data()), 1920, 1080);
    cFractionTOFMC->Divide(2, 2);

    for (unsigned int iBin = 0; iBin < nBins; iBin++)
    {
        for (auto &part : pdgNames)
        {
            intNsTPCPionMCV0tag[iBin][part.first] = hNsigmaTPCPionMCV0tag[iBin][part.first]->Integral() / hNsigmaTPCPionMCV0tag[iBin][part.first]->GetBinWidth(1);
            intNsTPCKaonMCKinktag[iBin][part.first] = hNsigmaTPCKaonMCKinktag[iBin][part.first]->Integral() / hNsigmaTPCKaonMCKinktag[iBin][part.first]->GetBinWidth(1);
            intNsTPCKaonMCTOFtag[iBin][part.first] = hNsigmaTPCKaonMCTOFtag[iBin][part.first]->Integral() / hNsigmaTPCKaonMCTOFtag[iBin][part.first]->GetBinWidth(1);
            intNsTPCProtonMCV0tag[iBin][part.first] = hNsigmaTPCProtonMCV0tag[iBin][part.first]->Integral() / hNsigmaTPCProtonMCV0tag[iBin][part.first]->GetBinWidth(1);
            intNsTOFPionMCV0tag[iBin][part.first] = hNsigmaTOFPionMCV0tag[iBin][part.first]->Integral() / hNsigmaTOFPionMCV0tag[iBin][part.first]->GetBinWidth(1);
            intNsTOFKaonMCKinktag[iBin][part.first] = hNsigmaTOFKaonMCKinktag[iBin][part.first]->Integral() / hNsigmaTOFKaonMCKinktag[iBin][part.first]->GetBinWidth(1);
            intNsTOFKaonMCTPCtag[iBin][part.first] = hNsigmaTOFKaonMCTPCtag[iBin][part.first]->Integral() / hNsigmaTOFKaonMCTPCtag[iBin][part.first]->GetBinWidth(1);
            intNsTOFProtonMCV0tag[iBin][part.first] = hNsigmaTOFProtonMCV0tag[iBin][part.first]->Integral() / hNsigmaTOFProtonMCV0tag[iBin][part.first]->GetBinWidth(1);
        }
        for (unsigned int iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++)
        {
            intNsTPCPionDataV0tag[iEtaBin][iBin] = hNsigmaTPCPionDataV0tag[iEtaBin][iBin]->Integral() / hNsigmaTPCPionDataV0tag[iEtaBin][iBin]->GetBinWidth(1);
            intNsTPCKaonDataKinktag[iEtaBin][iBin] = hNsigmaTPCKaonDataKinktag[iEtaBin][iBin]->Integral() / hNsigmaTPCKaonDataKinktag[iEtaBin][iBin]->GetBinWidth(1);
            intNsTPCKaonDataTOFtag[iEtaBin][iBin] = hNsigmaTPCKaonDataTOFtag[iEtaBin][iBin]->Integral() / hNsigmaTPCKaonDataTOFtag[iEtaBin][iBin]->GetBinWidth(1);
            intNsTPCProtonDataV0tag[iEtaBin][iBin] = hNsigmaTPCProtonDataV0tag[iEtaBin][iBin]->Integral() / hNsigmaTPCProtonDataV0tag[iEtaBin][iBin]->GetBinWidth(1);
            intNsTOFPionDataV0tag[iEtaBin][iBin] = hNsigmaTOFPionDataV0tag[iEtaBin][iBin]->Integral() / hNsigmaTOFPionDataV0tag[iEtaBin][iBin]->GetBinWidth(1);
            intNsTOFKaonDataKinktag[iEtaBin][iBin] = hNsigmaTOFKaonDataKinktag[iEtaBin][iBin]->Integral() / hNsigmaTOFKaonDataKinktag[iEtaBin][iBin]->GetBinWidth(1);
            intNsTOFKaonDataTPCtag[iEtaBin][iBin] = hNsigmaTOFKaonDataTPCtag[iEtaBin][iBin]->Integral() / hNsigmaTOFKaonDataTPCtag[iEtaBin][iBin]->GetBinWidth(1);
            intNsTOFProtonDataV0tag[iEtaBin][iBin] = hNsigmaTOFProtonDataV0tag[iEtaBin][iBin]->Integral() / hNsigmaTOFProtonDataV0tag[iEtaBin][iBin]->GetBinWidth(1);
        }
    }

    for (auto &part : pdgNames)
    {
        if(part.first == kAll)
            continue;

        hFracTPCPionMCV0tag[part.first] = new TH1D(Form("hFracTPCPionMCV0tag_%s_%s", part.second.data(), etaBinLabels[nEtaBins-1].Data()),
                                                            Form(";%s (GeV/#it{c});Contamination", varTitle.Data()), static_cast<int>(nBins), binLims);
        hFracTOFPionMCV0tag[part.first] = new TH1D(Form("hFracTOFPionMCV0tag_%s_%s", part.second.data(), etaBinLabels[nEtaBins-1].Data()),
                                                            Form(";%s (GeV/#it{c});Contamination", varTitle.Data()), static_cast<int>(nBins), binLims);
        hFracTPCKaonMCKinktag[part.first] = new TH1D(Form("hFracTPCKaonMCKinktag_%s_%s", part.second.data(), etaBinLabels[nEtaBins-1].Data()),
                                                              Form(";%s (GeV/#it{c});Contamination", varTitle.Data()), static_cast<int>(nBins), binLims);
        hFracTPCKaonMCTOFtag[part.first] = new TH1D(Form("hFracTPCKaonMCTOFtag_%s_%s", part.second.data(), etaBinLabels[nEtaBins-1].Data()),
                                                             Form(";%s (GeV/#it{c});Contamination", varTitle.Data()), static_cast<int>(nBins), binLims);
        hFracTOFKaonMCKinktag[part.first] = new TH1D(Form("hFracTOFKaonMCKinktag_%s_%s", part.second.data(), etaBinLabels[nEtaBins-1].Data()),
                                                              Form(";%s (GeV/#it{c});Contamination", varTitle.Data()), static_cast<int>(nBins), binLims);
        hFracTOFKaonMCTPCtag[part.first] = new TH1D(Form("hFracTOFKaonMCTPCtag_%s_%s", part.second.data(), etaBinLabels[nEtaBins-1].Data()),
                                                             Form(";%s (GeV/#it{c});Contamination", varTitle.Data()), static_cast<int>(nBins), binLims);
        hFracTPCProtonMCV0tag[part.first] = new TH1D(Form("hFracTPCProtonMCV0tag_%s_%s", part.second.data(), etaBinLabels[nEtaBins-1].Data()),
                                                              Form(";%s (GeV/#it{c});Contamination", varTitle.Data()), static_cast<int>(nBins), binLims);
        hFracTOFProtonMCV0tag[part.first] = new TH1D(Form("hFracTOFProtonMCV0tag_%s_%s", part.second.data(), etaBinLabels[nEtaBins-1].Data()),
                                                              Form(";%s (GeV/#it{c});Contamination", varTitle.Data()), static_cast<int>(nBins), binLims);

        SetTH1Style(hFracTPCPionMCV0tag[part.first], kFullCircle, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.045, 0.05);
        SetTH1Style(hFracTOFPionMCV0tag[part.first], kFullCircle, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.045, 0.05);
        SetTH1Style(hFracTPCKaonMCKinktag[part.first], kFullCircle, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.045, 0.05);
        SetTH1Style(hFracTPCKaonMCTOFtag[part.first], kFullCircle, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.045, 0.05);
        SetTH1Style(hFracTOFKaonMCKinktag[part.first], kFullCircle, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.045, 0.05);
        SetTH1Style(hFracTOFKaonMCTPCtag[part.first], kFullCircle, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.045, 0.05);
        SetTH1Style(hFracTPCProtonMCV0tag[part.first], kFullCircle, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.045, 0.05);
        SetTH1Style(hFracTOFProtonMCV0tag[part.first], kFullCircle, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.045, 0.05);
        
        legFracMC->AddEntry(hFracTOFProtonMCV0tag[part.first], part.second.data(), "lp");

        double eff = -1, unc = -1;            
        for (unsigned int iBin = 0; iBin < nBins; iBin++)
        {
            ComputeEfficiency(intNsTPCPionMCV0tag[iBin][part.first], intNsTPCPionMCV0tag[iBin][kAll], eff, unc);
            hFracTPCPionMCV0tag[part.first]->SetBinContent(iBin+1, eff);
            hFracTPCPionMCV0tag[part.first]->SetBinError(iBin+1, unc);

            ComputeEfficiency(intNsTOFPionMCV0tag[iBin][part.first], intNsTOFPionMCV0tag[iBin][kAll], eff, unc);
            hFracTOFPionMCV0tag[part.first]->SetBinContent(iBin+1, eff);
            hFracTOFPionMCV0tag[part.first]->SetBinError(iBin+1, unc);

            ComputeEfficiency(intNsTPCKaonMCKinktag[iBin][part.first], intNsTPCKaonMCKinktag[iBin][kAll], eff, unc);
            hFracTPCKaonMCKinktag[part.first]->SetBinContent(iBin+1, eff);
            hFracTPCKaonMCKinktag[part.first]->SetBinError(iBin+1, unc);

            ComputeEfficiency(intNsTPCKaonMCTOFtag[iBin][part.first], intNsTPCKaonMCTOFtag[iBin][kAll], eff, unc);
            hFracTPCKaonMCTOFtag[part.first]->SetBinContent(iBin+1, eff);
            hFracTPCKaonMCTOFtag[part.first]->SetBinError(iBin+1, unc);

            ComputeEfficiency(intNsTOFKaonMCKinktag[iBin][part.first], intNsTOFKaonMCKinktag[iBin][kAll], eff, unc);
            hFracTOFKaonMCKinktag[part.first]->SetBinContent(iBin+1, eff);
            hFracTOFKaonMCKinktag[part.first]->SetBinError(iBin+1, unc);

            ComputeEfficiency(intNsTOFKaonMCTPCtag[iBin][part.first], intNsTOFKaonMCTPCtag[iBin][kAll], eff, unc);
            hFracTOFKaonMCTPCtag[part.first]->SetBinContent(iBin+1, eff);
            hFracTOFKaonMCTPCtag[part.first]->SetBinError(iBin+1, unc);

            ComputeEfficiency(intNsTPCProtonMCV0tag[iBin][part.first], intNsTPCProtonMCV0tag[iBin][kAll], eff, unc);
            hFracTPCProtonMCV0tag[part.first]->SetBinContent(iBin+1, eff);
            hFracTPCProtonMCV0tag[part.first]->SetBinError(iBin+1, unc);

            ComputeEfficiency(intNsTOFProtonMCV0tag[iBin][part.first], intNsTOFProtonMCV0tag[iBin][kAll], eff, unc);
            hFracTOFProtonMCV0tag[part.first]->SetBinContent(iBin+1, eff);
            hFracTOFProtonMCV0tag[part.first]->SetBinError(iBin+1, unc);
        }
    }

    cFractionTPCMC->cd(1)->DrawFrame(binMins[0], 1.e-5, binMaxs[nBins-1], 10., Form("TPC #pi from K_{s}^{0};%s (GeV/#it{c});Purity / Contamination", varTitle.Data()));
    cFractionTPCMC->cd(1)->SetLogy();
    // cFractionTPCMC->cd(1)->SetLogx();
    for(auto &part : hFracTPCPionMCV0tag)
        hFracTPCPionMCV0tag[part.first]->DrawCopy("same");
    legFracMC->Draw();

    cFractionTPCMC->cd(2)->DrawFrame(binMins[0], 1.e-5, binMaxs[nBins-1], 10., Form("TPC K from kinks;%s (GeV/#it{c});Purity / Contamination", varTitle.Data()));
    cFractionTPCMC->cd(2)->SetLogy();
    // cFractionTPCMC->cd(2)->SetLogx();
    for(auto &part : hFracTPCKaonMCKinktag)
        hFracTPCKaonMCKinktag[part.first]->DrawCopy("same");

    cFractionTPCMC->cd(3)->DrawFrame(binMins[0], 1.e-5, binMaxs[nBins-1], 10., Form("TPC K TOF tagged;%s (GeV/#it{c});Purity / Contamination", varTitle.Data()));
    cFractionTPCMC->cd(3)->SetLogy();
    // cFractionTPCMC->cd(3)->SetLogx();
    for(auto &part : hFracTPCKaonMCTOFtag)
        hFracTPCKaonMCTOFtag[part.first]->DrawCopy("same");

    cFractionTPCMC->cd(4)->DrawFrame(binMins[0], 1.e-5, binMaxs[nBins-1], 10., Form("TPC p from #Lambda;%s (GeV/#it{c});Purity / Contamination", varTitle.Data()));
    cFractionTPCMC->cd(4)->SetLogy();
    // cFractionTPCMC->cd(4)->SetLogx();
    for(auto &part : hFracTPCProtonMCV0tag)
        hFracTPCProtonMCV0tag[part.first]->DrawCopy("same");
    gPad->Update();

    cFractionTOFMC->cd(1)->DrawFrame(binMins[0], 1.e-5, binMaxs[nBins-1], 10., Form("TOF #pi from K_{s}^{0};%s (GeV/#it{c});Purity / Contamination", varTitle.Data()));
    cFractionTOFMC->cd(1)->SetLogy();
    // cFractionTOFMC->cd(1)->SetLogx();
    for(auto &part : hFracTOFPionMCV0tag)
        hFracTOFPionMCV0tag[part.first]->DrawCopy("same");
    legFracMC->Draw();

    cFractionTOFMC->cd(2)->DrawFrame(binMins[0], 1.e-5, binMaxs[nBins-1], 10., Form("TOF K from kinks;%s (GeV/#it{c});Purity / Contamination", varTitle.Data()));
    cFractionTOFMC->cd(2)->SetLogy();
    // cFractionTOFMC->cd(2)->SetLogx();
    for(auto &part : hFracTOFKaonMCKinktag)
        hFracTOFKaonMCKinktag[part.first]->DrawCopy("same");

    cFractionTOFMC->cd(3)->DrawFrame(binMins[0], 1.e-5, binMaxs[nBins-1], 10., Form("TOF K TPC tag;%s (GeV/#it{c});Purity / Contamination", varTitle.Data()));
    cFractionTOFMC->cd(3)->SetLogy();
    // cFractionTOFMC->cd(3)->SetLogx();
    for(auto &part : hFracTOFKaonMCTPCtag)
        hFracTOFKaonMCTPCtag[part.first]->DrawCopy("same");

    cFractionTOFMC->cd(4)->DrawFrame(binMins[0], 1.e-5, binMaxs[nBins-1], 10., Form("TOF p from #Lambda;%s (GeV/#it{c});Purity / Contamination", varTitle.Data()));
    cFractionTOFMC->cd(4)->SetLogy();
    // cFractionTOFMC->cd(45)->SetLogx();
    for(auto &part : hFracTOFProtonMCV0tag)
        hFracTOFProtonMCV0tag[part.first]->DrawCopy("same");
    gPad->Update();

    cFractionTPCMC->SaveAs(Form("%s/PurityTPC_MC.pdf", outDirName.data()));
    cFractionTOFMC->SaveAs(Form("%s/PurityTOF_MC.pdf", outDirName.data()));

    std::cout << "\033[32mDone\033[0m" << std::endl;
    
    //compute fractions in data (TOF only)
    std::cout << "\n*******************************************\n" << std::endl;
    std::cout << "\033[32mCompute fractions and contaminations in data\033[0m\n" << std::endl;

    std::array<std::array<TFractionFitter*, nBinsMax>, nEtaBinsMax+1> fNsTOFKaonTPCtagFitter, fNsTOFProtonV0tagFitter;
    std::array<std::map<int, TH1D*>, nEtaBinsMax+1> hFracTOFKaonDataTPCtag, hFracTOFProtonDataV0tag;
    std::array<std::array<std::map<int, TH1D*>, nBinsMax>, nEtaBinsMax+1> hNsigmaTOFKaonDataTPCtagFit, hNsigmaTOFProtonDataV0tagFit;
    std::array<TCanvas*, nEtaBinsMax+1> cFitResultTOFKaonFromTPCtag, cFitResultTOFProtonFromV0tag, cTOFFractionData;

    TLegend *legTOFFitter = new TLegend(0.14, 0.62, 0.44, 0.86);
    legTOFFitter->SetTextSize(0.04);
    legTOFFitter->AddEntry(hNsigmaTOFKaonDataTPCtag[0][0], "Data", "p");

    TLegend *legFracData = new TLegend(0.3, 0.3, 0.8, 0.8);
    legFracData->SetTextSize(0.05);

    for(unsigned int iEtaBin = 0; iEtaBin<nEtaBins; iEtaBin++)
    {
        cFitResultTOFKaonFromTPCtag[iEtaBin] = new TCanvas(Form("cFitResultTOFKaonFromTPCtag_%s", etaBinLabels[iEtaBin].Data()), Form("cFitResultTOFKaonFromTPCtag_%s", etaBinLabels[iEtaBin].Data()), 1920, 1080);
        DivideCanvas(cFitResultTOFKaonFromTPCtag[iEtaBin], nBins);
        cFitResultTOFProtonFromV0tag[iEtaBin] = new TCanvas(Form("cFitResultTOFProtonFromV0tag_%s", etaBinLabels[iEtaBin].Data()), Form("cFitResultTOFProtonFromV0tag_%s", etaBinLabels[iEtaBin].Data()), 1920, 1080);
        DivideCanvas(cFitResultTOFProtonFromV0tag[iEtaBin], nBins);
        cTOFFractionData[iEtaBin] = new TCanvas(Form("cTOFFractionData_%s", etaBinLabels[iEtaBin].Data()), Form("cTOFFractionData_%s", etaBinLabels[iEtaBin].Data()), 1920, 1080);
        cTOFFractionData[iEtaBin]->Divide(3, 1);

        for (auto &part : pdgNames)
        {
            if(part.first == kAll)
                continue;

            if(iEtaBin == 0)
                legTOFFitter->AddEntry(hNsigmaTOFKaonMCTPCtag[0][part.first], Form("Templ %s", part.second.data()), "l");

            hFracTOFKaonDataTPCtag[iEtaBin][part.first] = new TH1D(Form("hFracTOFKaonDataTPCtag_%s_%s", part.second.data(), etaBinLabels[iEtaBin].Data()), Form(";%s (GeV/#it{c});Fraction", varTitle.Data()), static_cast<int>(nBins), binLims);
            hFracTOFProtonDataV0tag[iEtaBin][part.first] = new TH1D(Form("hFracTOFProtonDataV0tag_%s_%s", part.second.data(), etaBinLabels[iEtaBin].Data()), Form(";%s (GeV/#it{c});Fraction", varTitle.Data()), static_cast<int>(nBins), binLims);
            SetTH1Style(hFracTOFKaonDataTPCtag[iEtaBin][part.first], kOpenSquare, pdgColors[part.first] + 1, 1., 2, pdgColors[part.first], kWhite, 0.045, 0.05);
            SetTH1Style(hFracTOFProtonDataV0tag[iEtaBin][part.first], kOpenSquare, pdgColors[part.first] + 1, 1., 2, pdgColors[part.first], kWhite, 0.045, 0.05);
        }
        for (unsigned int iBin = 0; iBin < nBins; iBin++)
        {
            std::vector<int> templUsedTOFKaonMCTPCtag;
            GetTOFFractionsFromData(kKaon, iBin, hFracTOFKaonMCTPCtag, hFracTOFKaonDataTPCtag[iEtaBin], hNsigmaTOFKaonMCTPCtag[iBin], hNsigmaTOFKaonDataTPCtag[iEtaBin][iBin], fNsTOFKaonTPCtagFitter[iEtaBin][iBin], templUsedTOFKaonMCTPCtag);

            if (templUsedTOFKaonMCTPCtag.size() > 1)
            {
                hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][kAll] = dynamic_cast<TH1D *>(fNsTOFKaonTPCtagFitter[iEtaBin][iBin]->GetPlot());
                hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][kAll]->SetName(Form("hNsigmaTOFKaonDataTPCtagFit_All_%s_%s", etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()));
                SetTH1Style(hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][kAll], 0, kRed, 1., 2, kRed, kWhite);

                for (unsigned int iTempl = 0; iTempl < templUsedTOFKaonMCTPCtag.size(); iTempl++)
                {
                    double frac, err;
                    fNsTOFKaonTPCtagFitter[iEtaBin][iBin]->GetResult(iTempl, frac, err);
                    hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][templUsedTOFKaonMCTPCtag[iTempl]] = dynamic_cast<TH1D *>(fNsTOFKaonTPCtagFitter[iEtaBin][iBin]->GetMCPrediction(iTempl));
                    hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][templUsedTOFKaonMCTPCtag[iTempl]]->SetName(Form("hNsigmaTOFKaonDataTPCtagFit_%s_%s_%s", pdgNames[templUsedTOFKaonMCTPCtag[iTempl]].data(), etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()));
                    hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][templUsedTOFKaonMCTPCtag[iTempl]]->Scale(frac / hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][templUsedTOFKaonMCTPCtag[iTempl]]->Integral() * hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][kAll]->Integral());
                    hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][templUsedTOFKaonMCTPCtag[iTempl]]->SetFillColor(kWhite);
                    hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][templUsedTOFKaonMCTPCtag[iTempl]]->SetFillStyle(0);
                }
            }
            else
            {
                hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][kAll] = dynamic_cast<TH1D *>(hNsigmaTOFKaonDataTPCtag[iEtaBin][iBin]->Clone());
                SetTH1Style(hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][kAll], 0, kRed, 1., 2, kRed, kWhite);
                hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][kKaon] = dynamic_cast<TH1D *>(hNsigmaTOFKaonDataTPCtag[iEtaBin][iBin]->Clone());
                hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][kKaon]->SetFillColor(kWhite);
                hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][kKaon]->SetFillStyle(0);
            }

            if(nBins>1)
                cFitResultTOFKaonFromTPCtag[iEtaBin]->cd(iBin+1)->SetLogy();
            else
                cFitResultTOFKaonFromTPCtag[iEtaBin]->cd()->SetLogy();
            hNsigmaTOFKaonDataTPCtag[iEtaBin][iBin]->DrawCopy("e");
            for(auto &part : pdgNames)
            {
                std::vector<int>::iterator it = find(templUsedTOFKaonMCTPCtag.begin(), templUsedTOFKaonMCTPCtag.end(), part.first);
                if (it != templUsedTOFKaonMCTPCtag.end())
                    hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][part.first]->DrawCopy("hist same");
                else
                    hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][part.first] = nullptr;
            }
            legTOFFitter->Draw();

            std::vector<int> templUsedTOFProtonMCV0tag;
            GetTOFFractionsFromData(kPr, iBin, hFracTOFProtonMCV0tag, hFracTOFProtonDataV0tag[iEtaBin], hNsigmaTOFProtonMCV0tag[iBin], hNsigmaTOFProtonDataV0tag[iEtaBin][iBin], fNsTOFProtonV0tagFitter[iEtaBin][iBin], templUsedTOFProtonMCV0tag);

            if (templUsedTOFProtonMCV0tag.size() > 1)
            {
                hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][kAll] = dynamic_cast<TH1D *>(fNsTOFProtonV0tagFitter[iEtaBin][iBin]->GetPlot());
                hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][kAll]->SetName(Form("hNsigmaTOFProtonDataV0tagFit_All_%s_%s", etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()));
                SetTH1Style(hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][kAll], 0, kRed, 1., 2, kRed, kWhite);

                for (unsigned int iTempl = 0; iTempl < templUsedTOFProtonMCV0tag.size(); iTempl++)
                {
                    double frac, err;
                    fNsTOFProtonV0tagFitter[iEtaBin][iBin]->GetResult(iTempl, frac, err);
                    hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][templUsedTOFProtonMCV0tag[iTempl]] = dynamic_cast<TH1D *>(fNsTOFProtonV0tagFitter[iEtaBin][iBin]->GetMCPrediction(iTempl));
                    hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][templUsedTOFProtonMCV0tag[iTempl]]->SetName(Form("hNsigmaTOFProtonDataV0tagFit_%s_%s_%s", pdgNames[templUsedTOFProtonMCV0tag[iTempl]].data(), etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()));
                    hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][templUsedTOFProtonMCV0tag[iTempl]]->Scale(frac / hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][templUsedTOFProtonMCV0tag[iTempl]]->Integral() * hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][kAll]->Integral());
                    hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][templUsedTOFProtonMCV0tag[iTempl]]->SetFillColor(kWhite);
                    hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][templUsedTOFProtonMCV0tag[iTempl]]->SetFillStyle(0);
                }
            }
            else
            {
                hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][kAll] = dynamic_cast<TH1D *>(hNsigmaTOFProtonDataV0tag[iEtaBin][iBin]->Clone());
                SetTH1Style(hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][kAll], 0, kRed, 1., 2, kRed, kWhite);
                hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][kPr] = dynamic_cast<TH1D *>(hNsigmaTOFProtonDataV0tag[iEtaBin][iBin]->Clone());
                hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][kPr]->SetFillColor(kWhite);
                hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][kPr]->SetFillStyle(0);
            }
            if(nBins>1)
                cFitResultTOFProtonFromV0tag[iEtaBin]->cd(iBin+1)->SetLogy();
            else
                cFitResultTOFProtonFromV0tag[iEtaBin]->cd()->SetLogy();
            hNsigmaTOFProtonDataV0tag[iEtaBin][iBin]->DrawCopy("E");
            for(auto &part : pdgNames)
            {
                std::vector<int>::iterator it = find(templUsedTOFProtonMCV0tag.begin(), templUsedTOFProtonMCV0tag.end(), part.first);
                if (it != templUsedTOFProtonMCV0tag.end())
                    hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][part.first]->DrawCopy("hist same");
                else
                    hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][part.first] = nullptr;
            }
            legTOFFitter->Draw();
        }
        cFitResultTOFKaonFromTPCtag[iEtaBin]->SaveAs(Form("%s/FitNsigmaDistrTOFKaonFromTPCtag_Data_%s.pdf", outDirName.data(), etaBinLabels[iEtaBin].Data()));
        cFitResultTOFProtonFromV0tag[iEtaBin]->SaveAs(Form("%s/FitNsigmaDistrTOFProtonFromV0tag_Data_%s.pdf", outDirName.data(), etaBinLabels[iEtaBin].Data()));
    
        cTOFFractionData[iEtaBin]->cd(1)->DrawFrame(binMins[0], 1.e-5, binMaxs[nBins-1], 10., Form("TOF K TPC tag;%s (GeV/#it{c});Purity / Contamination", varTitle.Data()));
        cTOFFractionData[iEtaBin]->cd(1)->SetLogy();
        // cTOFFractionData[iEtaBin]->cd(1)->SetLogx();
        for (auto &part : hFracTOFKaonDataTPCtag[iEtaBin])
            hFracTOFKaonDataTPCtag[iEtaBin][part.first]->DrawCopy("same");
        cTOFFractionData[iEtaBin]->cd(2)->DrawFrame(binMins[0], 1.e-5, binMaxs[nBins-1], 10., Form("TOF p from #Lambda;%s (GeV/#it{c});Purity / Contamination", varTitle.Data()));
        cTOFFractionData[iEtaBin]->cd(2)->SetLogy();
        // cTOFFractionData[iEtaBin]->cd(2)->SetLogx();
        for (auto &part : hFracTOFProtonDataV0tag[iEtaBin])
            hFracTOFProtonDataV0tag[iEtaBin][part.first]->DrawCopy("same");
        if(iEtaBin == 0)
        {
            for (auto &part : hFracTOFProtonDataV0tag[iEtaBin])
                legFracData->AddEntry(hFracTOFProtonDataV0tag[iEtaBin][part.first], pdgNames[part.first].data(), "lp");
        }
        cTOFFractionData[iEtaBin]->cd(3);
        legFracData->Draw();
        cTOFFractionData[iEtaBin]->SaveAs(Form("%s/PurityTOF_Data_%s.pdf", outDirName.data(), etaBinLabels[iEtaBin].Data()));
    }

    std::cout << "\033[32mDone\033[0m" << std::endl;

    //draw and fit distributions
    std::cout << "\n*******************************************\n" << std::endl;
    std::cout << "\033[32mDraw distributions\033[0m\n" << std::endl;

    //normalise histograms
    for (unsigned int iBin = 0; iBin < nBins; iBin++)
    {
        for (auto part = pdgNames.rbegin(); part != pdgNames.rend(); ++part)
        {
            hNsigmaTPCPionMCV0tag[iBin][part->first]->Scale(1./hNsigmaTPCPionMCV0tag[iBin][kAll]->Integral());
            hNsigmaTPCKaonMCKinktag[iBin][part->first]->Scale(1./hNsigmaTPCKaonMCKinktag[iBin][kAll]->Integral());
            hNsigmaTPCKaonMCTOFtag[iBin][part->first]->Scale(1./hNsigmaTPCKaonMCTOFtag[iBin][kAll]->Integral());
            hNsigmaTPCProtonMCV0tag[iBin][part->first]->Scale(1./hNsigmaTPCProtonMCV0tag[iBin][kAll]->Integral());
            hNsigmaTOFPionMCV0tag[iBin][part->first]->Scale(1./hNsigmaTOFPionMCV0tag[iBin][kAll]->Integral());
            hNsigmaTOFKaonMCKinktag[iBin][part->first]->Scale(1./hNsigmaTOFKaonMCKinktag[iBin][kAll]->Integral());
            hNsigmaTOFKaonMCTPCtag[iBin][part->first]->Scale(1./hNsigmaTOFKaonMCTPCtag[iBin][kAll]->Integral());
            hNsigmaTOFProtonMCV0tag[iBin][part->first]->Scale(1./hNsigmaTOFProtonMCV0tag[iBin][kAll]->Integral());
        }

        for(unsigned int iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++)
        {
            for (auto &part : pdgNames)
            {
                if (part.first != kAll && hFracTOFKaonDataTPCtag[iEtaBin][part.first]->GetBinContent(iBin+1)>0)
                    hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][part.first]->Scale(1./hNsigmaTOFKaonDataTPCtag[iEtaBin][iBin]->Integral());
                if (part.first != kAll && hFracTOFProtonDataV0tag[iEtaBin][part.first]->GetBinContent(iBin+1)>0)
                    hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][part.first]->Scale(1./hNsigmaTOFProtonDataV0tag[iEtaBin][iBin]->Integral());
            }

            hNsigmaTPCPionDataV0tag[iEtaBin][iBin]->Scale(1./hNsigmaTPCPionDataV0tag[iEtaBin][iBin]->Integral());
            hNsigmaTPCKaonDataKinktag[iEtaBin][iBin]->Scale(1./hNsigmaTPCKaonDataKinktag[iEtaBin][iBin]->Integral());
            hNsigmaTPCKaonDataTOFtag[iEtaBin][iBin]->Scale(1./hNsigmaTPCKaonDataTOFtag[iEtaBin][iBin]->Integral());
            hNsigmaTPCProtonDataV0tag[iEtaBin][iBin]->Scale(1./hNsigmaTPCProtonDataV0tag[iEtaBin][iBin]->Integral());
            hNsigmaTOFPionDataV0tag[iEtaBin][iBin]->Scale(1./hNsigmaTOFPionDataV0tag[iEtaBin][iBin]->Integral());
            hNsigmaTOFKaonDataKinktag[iEtaBin][iBin]->Scale(1./hNsigmaTOFKaonDataKinktag[iEtaBin][iBin]->Integral());
            hNsigmaTOFKaonDataTPCtag[iEtaBin][iBin]->Scale(1./hNsigmaTOFKaonDataTPCtag[iEtaBin][iBin]->Integral());
            hNsigmaTOFProtonDataV0tag[iEtaBin][iBin]->Scale(1./hNsigmaTOFProtonDataV0tag[iEtaBin][iBin]->Integral());
        }
    }

    // MC
    TString fitopt = "ML0RQ";
    std::array<std::map<int, TF1*>, nBinsMax> fNsigmaTPCPionMCV0tag, fNsigmaTPCKaonMCKinktag, fNsigmaTPCKaonMCTOFtag, fNsigmaTPCProtonMCV0tag;

    TCanvas *cPionMCV0tagTPC = new TCanvas(Form("cPionMCV0tagTPC_%s", etaBinLabels[nEtaBins-1].Data()), Form("cPionMCV0tagTPC_%s", etaBinLabels[nEtaBins-1].Data()), 1920, 1080);
    DivideCanvas(cPionMCV0tagTPC, nBins);
    TCanvas *cKaonMCKinkstagTPC = new TCanvas(Form("cKaonMCKinkstagTPC_%s", etaBinLabels[nEtaBins-1].Data()), Form("cKaonMCKinkstagTPC_%s", etaBinLabels[nEtaBins-1].Data()), 1920, 1080);
    DivideCanvas(cKaonMCKinkstagTPC, nBins);
    TCanvas *cKaonMCTOFtagTPC = new TCanvas(Form("cKaonMCTOFtagTPC_%s", etaBinLabels[nEtaBins-1].Data()), Form("cKaonMCTOFtagTPC_%s", etaBinLabels[nEtaBins-1].Data()), 1920, 1080);
    DivideCanvas(cKaonMCTOFtagTPC, nBins);
    TCanvas *cProtonMCV0tagTPC = new TCanvas(Form("cProtonMCV0tagTPC_%s", etaBinLabels[nEtaBins-1].Data()), Form("cProtonMCV0tagTPC_%s", etaBinLabels[nEtaBins-1].Data()), 1920, 1080);
    DivideCanvas(cProtonMCV0tagTPC, nBins);
    TCanvas *cPionMCV0tagTOF = new TCanvas(Form("cPionMCV0tagTOF_%s", etaBinLabels[nEtaBins-1].Data()), Form("cPionMCV0tagTOF_%s", etaBinLabels[nEtaBins-1].Data()), 1920, 1080);
    DivideCanvas(cPionMCV0tagTOF, nBins);
    TCanvas *cKaonMCKinkstagTOF = new TCanvas(Form("cKaonMCKinkstagTOF_%s", etaBinLabels[nEtaBins-1].Data()), Form("cKaonMCKinkstagTOF_%s", etaBinLabels[nEtaBins-1].Data()), 1920, 1080);
    DivideCanvas(cKaonMCKinkstagTOF, nBins);
    TCanvas *cKaonMCTPCtagTOF = new TCanvas(Form("cKaonMCTPCtagTOF_%s", etaBinLabels[nEtaBins-1].Data()), Form("cKaonMCTPCtagTOF_%s", etaBinLabels[nEtaBins-1].Data()), 1920, 1080);
    DivideCanvas(cKaonMCTPCtagTOF, nBins);
    TCanvas *cProtonMCV0tagTOF = new TCanvas(Form("cProtonMCV0tagTOF_%s", etaBinLabels[nEtaBins-1].Data()), Form("cProtonMCV0tagTOF_%s", etaBinLabels[nEtaBins-1].Data()), 1920, 1080);
    DivideCanvas(cProtonMCV0tagTOF, nBins);

    TLegend *legMCdist = new TLegend(0.14, 0.62, 0.44, 0.86);
    legMCdist->SetTextSize(0.04);
    legMCdist->AddEntry(hNsigmaTPCPionMCV0tag[0][kAll], "MC All", "p");
    for(auto &part : pdgNames)
    {
        if(part.first == kAll)
            continue;
        legMCdist->AddEntry(hNsigmaTPCPionMCV0tag[0][part.first], Form("MC %s", part.second.data()), "f");
    }
    for (unsigned int iBin = 0; iBin < nBins; iBin++)
    {
        for(auto &part : pdgNames)
        {
            if(part.first == kAll)
                continue;
            fNsigmaTPCPionMCV0tag[iBin][part.first] = new TF1(Form("fNsigmaTPCPionMCV0tag_%s_%s", part.second.data(), binLabels[iBin].Data()), "gaus", -50, 50);
            fNsigmaTPCKaonMCKinktag[iBin][part.first] = new TF1(Form("fNsigmaTPCKaonMCKinktag_%s_%s", part.second.data(), binLabels[iBin].Data()), "gaus", -50, 50);
            fNsigmaTPCKaonMCTOFtag[iBin][part.first] = new TF1(Form("fNsigmaTPCKaonMCTOFtag_%s_%s", part.second.data(), binLabels[iBin].Data()), "gaus", -50, 50);
            fNsigmaTPCProtonMCV0tag[iBin][part.first] = new TF1(Form("fNsigmaTPCProtonMCV0tag_%s_%s", part.second.data(), binLabels[iBin].Data()), "gaus", -50, 50);
        }

        cPionMCV0tagTPC->cd(iBin+1)->SetLogy();
        hNsigmaTPCPionMCV0tag[iBin][kAll]->GetXaxis()->SetNdivisions(505);
        hNsigmaTPCPionMCV0tag[iBin][kAll]->DrawCopy("E");
        for(auto &part : pdgNames)
        {
            if(part.first == kAll)
                continue;
            hNsigmaTPCPionMCV0tag[iBin][part.first]->DrawCopy("hist same");
            if (part.first == kPion || hFracTPCPionMCV0tag[part.first]->GetBinContent(iBin+1) > 1.e-5)
                hNsigmaTPCPionMCV0tag[iBin][part.first]->Fit(fNsigmaTPCPionMCV0tag[iBin][part.first], fitopt.Data());
            else
                fNsigmaTPCPionMCV0tag[iBin][part.first]->SetParameter(0, 0);
        }
        hNsigmaTPCPionMCV0tag[iBin][kAll]->DrawCopy("Esame");
        legMCdist->Draw();

        cKaonMCKinkstagTPC->cd(iBin+1)->SetLogy();
        hNsigmaTPCKaonMCKinktag[iBin][kAll]->GetXaxis()->SetNdivisions(505);
        hNsigmaTPCKaonMCKinktag[iBin][kAll]->DrawCopy("E");
        for(auto &part : pdgNames)
        {
            if(part.first == kAll)
                continue;
            hNsigmaTPCKaonMCKinktag[iBin][part.first]->DrawCopy("hist same");
            if (part.first == kKaon || hFracTPCKaonMCKinktag[part.first]->GetBinContent(iBin+1) > 1.e-5)
                hNsigmaTPCKaonMCKinktag[iBin][part.first]->Fit(fNsigmaTPCKaonMCKinktag[iBin][part.first], fitopt.Data());
            else
                fNsigmaTPCKaonMCKinktag[iBin][part.first]->SetParameter(0, 0);
        }
        hNsigmaTPCKaonMCKinktag[iBin][kAll]->DrawCopy("Esame");
        legMCdist->Draw();

        cKaonMCTOFtagTPC->cd(iBin+1)->SetLogy();
        hNsigmaTPCKaonMCTOFtag[iBin][kAll]->GetXaxis()->SetNdivisions(505);
        hNsigmaTPCKaonMCTOFtag[iBin][kAll]->DrawCopy("E");
        for(auto &part : pdgNames)
        {
            if(part.first == kAll)
                continue;
            hNsigmaTPCKaonMCTOFtag[iBin][part.first]->DrawCopy("hist same");
            if (part.first == kKaon || hFracTPCKaonMCTOFtag[part.first]->GetBinContent(iBin+1) > 1.e-5)
                hNsigmaTPCKaonMCTOFtag[iBin][part.first]->Fit(fNsigmaTPCKaonMCTOFtag[iBin][part.first], fitopt.Data());
            else
                fNsigmaTPCKaonMCTOFtag[iBin][part.first]->SetParameter(0, 0);
        }
        hNsigmaTPCKaonMCTOFtag[iBin][kAll]->DrawCopy("Esame");
        legMCdist->Draw();

        cProtonMCV0tagTPC->cd(iBin+1)->SetLogy();
        hNsigmaTPCProtonMCV0tag[iBin][kAll]->GetXaxis()->SetNdivisions(505);
        hNsigmaTPCProtonMCV0tag[iBin][kAll]->DrawCopy("E");
        for(auto &part : pdgNames)
        {
            if(part.first == kAll)
                continue;
            hNsigmaTPCProtonMCV0tag[iBin][part.first]->DrawCopy("hist same");
            if (part.first == kPr || hFracTPCProtonMCV0tag[part.first]->GetBinContent(iBin+1) > 1.e-5)
                hNsigmaTPCProtonMCV0tag[iBin][part.first]->Fit(fNsigmaTPCProtonMCV0tag[iBin][part.first], fitopt.Data());
            else
                fNsigmaTPCProtonMCV0tag[iBin][part.first]->SetParameter(0, 0);
        }
        hNsigmaTPCProtonMCV0tag[iBin][kAll]->DrawCopy("Esame");
        legMCdist->Draw();

        cPionMCV0tagTOF->cd(iBin+1)->SetLogy();
        hNsigmaTOFPionMCV0tag[iBin][kAll]->GetXaxis()->SetNdivisions(505);
        hNsigmaTOFPionMCV0tag[iBin][kAll]->DrawCopy("E");
        for(auto &part : pdgNames)
        {
            if(part.first == kAll)
                continue;
            hNsigmaTOFPionMCV0tag[iBin][part.first]->DrawCopy("hist same");
        }
        hNsigmaTOFPionMCV0tag[iBin][kAll]->DrawCopy("Esame");
        legMCdist->Draw();

        cKaonMCKinkstagTOF->cd(iBin+1)->SetLogy();
        hNsigmaTOFKaonMCKinktag[iBin][kAll]->GetXaxis()->SetNdivisions(505);
        hNsigmaTOFKaonMCKinktag[iBin][kAll]->DrawCopy("E");
        for(auto &part : pdgNames)
        {
            if(part.first == kAll)
                continue;
            hNsigmaTOFKaonMCKinktag[iBin][part.first]->DrawCopy("hist same");
        }
        hNsigmaTOFKaonMCKinktag[iBin][kAll]->DrawCopy("Esame");
        legMCdist->Draw();

        cKaonMCTPCtagTOF->cd(iBin+1)->SetLogy();
        hNsigmaTOFKaonMCTPCtag[iBin][kAll]->GetXaxis()->SetNdivisions(505);
        hNsigmaTOFKaonMCTPCtag[iBin][kAll]->DrawCopy("E");
        for(auto &part : pdgNames)
        {
            if(part.first == kAll)
                continue;
            hNsigmaTOFKaonMCTPCtag[iBin][part.first]->DrawCopy("hist same");
        }
        hNsigmaTOFKaonMCTPCtag[iBin][kAll]->DrawCopy("Esame");
        legMCdist->Draw();

        cProtonMCV0tagTOF->cd(iBin+1)->SetLogy();
        hNsigmaTOFProtonMCV0tag[iBin][kAll]->GetXaxis()->SetNdivisions(505);
        hNsigmaTOFProtonMCV0tag[iBin][kAll]->DrawCopy("E");
        for(auto &part : pdgNames)
        {
            if(part.first == kAll)
                continue;
            hNsigmaTOFProtonMCV0tag[iBin][part.first]->DrawCopy("hist same");
        }
        hNsigmaTOFProtonMCV0tag[iBin][kAll]->DrawCopy("Esame");
        legMCdist->Draw();
    }

    cPionMCV0tagTPC->SaveAs(Form("%s/DistrNsigmaTPCPionFromV0tag_MC.pdf", outDirName.data()));
    cKaonMCKinkstagTPC->SaveAs(Form("%s/DistrNsigmaTPCKaonFromKinkstag_MC.pdf", outDirName.data()));
    cKaonMCTOFtagTPC->SaveAs(Form("%s/DistrNsigmaTPCKaonFromTOFtag_MC.pdf", outDirName.data()));
    cProtonMCV0tagTPC->SaveAs(Form("%s/DistrNsigmaTPCProtonFromV0tag_MC.pdf", outDirName.data()));
    cPionMCV0tagTOF->SaveAs(Form("%s/DistrNsigmaTOFPionFromV0tag_MC.pdf", outDirName.data()));
    cKaonMCKinkstagTOF->SaveAs(Form("%s/DistrNsigmaTOFKaonFromKinkstag_MC.pdf", outDirName.data()));
    cKaonMCTPCtagTOF->SaveAs(Form("%s/DistrNsigmaTOFKaonFromTPCtag_MC.pdf", outDirName.data()));
    cProtonMCV0tagTOF->SaveAs(Form("%s/DistrNsigmaTOFProtonFromV0tag_MC.pdf", outDirName.data()));

    // Data
    std::array<TCanvas*, nEtaBinsMax+1> cPionDataV0tagTPC, cKaonDataKinkstagTPC, cKaonDataTOFtagTPC, cProtonDataV0tagTPC;
    std::array<TCanvas*, nEtaBinsMax+1> cPionDataV0tagTOF, cKaonDataKinkstagTOF, cKaonDataTPCtagTOF, cProtonDataV0tagTOF;

    std::array<std::array<std::map<int, TF1*>, nBinsMax>, nEtaBinsMax+1> fNsigmaTPCPionDataV0tag, fNsigmaTPCKaonDataKinktag, fNsigmaTPCKaonDataTOFtag, fNsigmaTPCProtonDataV0tag;

    //subtract bkg template for TOF nsigma
    std::array<std::array<TH1D*, nBinsMax>, nEtaBinsMax+1> hNsigmaTOFPionDataV0tagSub, hNsigmaTOFKaonDataKinktagSub, hNsigmaTOFKaonDataTPCtagSub, hNsigmaTOFProtonDataV0tagSub;

    TLegend *legTPCFitter = new TLegend(0.14, 0.62, 0.48, 0.86);
    legTPCFitter->SetTextSize(0.04);
    legTPCFitter->AddEntry(hNsigmaTPCPionDataV0tag[0][0], "Data", "p");

    TLegend *legTOFPionDataV0tag = new TLegend(0.14, 0.7, 0.48, 0.8);
    legTOFPionDataV0tag->SetTextSize(0.05);

    TLegend *legTOFKaonDataKinkstag = new TLegend(0.14, 0.7, 0.48, 0.8);
    legTOFKaonDataKinkstag->SetTextSize(0.05);

    TLegend *legTOFKaonDataTPCtag = new TLegend(0.14, 0.7, 0.48, 0.8);
    legTOFKaonDataTPCtag->SetTextSize(0.05);

    TLegend *legTOFProtonDataV0tag = new TLegend(0.14, 0.7, 0.48, 0.8);
    legTOFProtonDataV0tag->SetTextSize(0.05);

    for(unsigned int iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++)
    {
        cPionDataV0tagTPC[iEtaBin] = new TCanvas(Form("cPionDataV0tagTPC_%s", etaBinLabels[iEtaBin].Data()), Form("cPionDataV0tagTPC_%s", etaBinLabels[iEtaBin].Data()), 1920, 1080);
        DivideCanvas(cPionDataV0tagTPC[iEtaBin], nBins);
        cKaonDataKinkstagTPC[iEtaBin] = new TCanvas(Form("cKaonDataKinkstagTPC_%s", etaBinLabels[iEtaBin].Data()), Form("cKaonDataKinkstagTPC_%s", etaBinLabels[iEtaBin].Data()), 1920, 1080);
        DivideCanvas(cKaonDataKinkstagTPC[iEtaBin], nBins);
        cKaonDataTOFtagTPC[iEtaBin] = new TCanvas(Form("cKaonDataTOFtagTPC_%s", etaBinLabels[iEtaBin].Data()), Form("cKaonDataTOFtagTPC_%s", etaBinLabels[iEtaBin].Data()), 1920, 1080);
        DivideCanvas(cKaonDataTOFtagTPC[iEtaBin], nBins);
        cProtonDataV0tagTPC[iEtaBin] = new TCanvas(Form("cProtonDataV0tagTPC_%s", etaBinLabels[iEtaBin].Data()), Form("cProtonDataV0tagTPC_%s", etaBinLabels[iEtaBin].Data()), 1920, 1080);
        DivideCanvas(cProtonDataV0tagTPC[iEtaBin], nBins);
        cPionDataV0tagTOF[iEtaBin] = new TCanvas(Form("cPionDataV0tagTOF_%s", etaBinLabels[iEtaBin].Data()), Form("cPionDataV0tagTOF_%s", etaBinLabels[iEtaBin].Data()), 1920, 1080);
        DivideCanvas(cPionDataV0tagTOF[iEtaBin], nBins);
        cKaonDataKinkstagTOF[iEtaBin] = new TCanvas(Form("cKaonDataKinkstagTOF_%s", etaBinLabels[iEtaBin].Data()), Form("cKaonDataKinkstagTOF_%s", etaBinLabels[iEtaBin].Data()), 1920, 1080);
        DivideCanvas(cKaonDataKinkstagTOF[iEtaBin], nBins);
        cKaonDataTPCtagTOF[iEtaBin] = new TCanvas(Form("cKaonDataTPCtagTOF_%s", etaBinLabels[iEtaBin].Data()), Form("cKaonDataTPCtagTOF_%s", etaBinLabels[iEtaBin].Data()), 1920, 1080);
        DivideCanvas(cKaonDataTPCtagTOF[iEtaBin], nBins);
        cProtonDataV0tagTOF[iEtaBin] = new TCanvas(Form("cProtonDataV0tagTOF_%s", etaBinLabels[iEtaBin].Data()), Form("cProtonDataV0tagTOF_%s", etaBinLabels[iEtaBin].Data()), 1920, 1080);
        DivideCanvas(cProtonDataV0tagTOF[iEtaBin], nBins);

        for (unsigned int iBin = 0; iBin < nBins; iBin++)
        {
            fNsigmaTPCPionDataV0tag[iEtaBin][iBin][kAll] = new TF1(Form("fNsigmaTPCPionDataV0tag_%s_%s_%s", pdgNames[kAll].data(), etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()), PDFnsigmaTPCtot, -50, 50, pdgPosition[kAll] * 3);
            fNsigmaTPCKaonDataKinktag[iEtaBin][iBin][kAll] = new TF1(Form("fNsigmaTPCKaonDataKinktag_%s_%s_%s", pdgNames[kAll].data(), etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()), PDFnsigmaTPCtot, -50, 50, pdgPosition[kAll] * 3);
            fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][kAll] = new TF1(Form("fNsigmaTPCKaonDataTOFtag_%s_%s_%s", pdgNames[kAll].data(), etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()), PDFnsigmaTPCtot, -50, 50, pdgPosition[kAll] * 3);
            fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][kAll] = new TF1(Form("fNsigmaTPCProtonDataV0tag_%s_%s_%s", pdgNames[kAll].data(), etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()), PDFnsigmaTPCtot, -50, 50, pdgPosition[kAll] * 3);

            for(auto &part : pdgPosition)
            {
                if(part.first == kAll)
                    continue;

                double toll = 100.;
                for (int iPar = 0; iPar < 3; iPar++)
                {
                    if (iPar == 0)
                        toll = 5.;
                    else
                        toll = 0.5;

                    fNsigmaTPCPionDataV0tag[iEtaBin][iBin][kAll]->SetParameter(part.second * 3 + iPar, fNsigmaTPCPionMCV0tag[iBin][part.first]->GetParameter(iPar));
                    fNsigmaTPCKaonDataKinktag[iEtaBin][iBin][kAll]->SetParameter(part.second * 3 + iPar, fNsigmaTPCKaonMCKinktag[iBin][part.first]->GetParameter(iPar));
                    fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][kAll]->SetParameter(part.second * 3 + iPar, fNsigmaTPCKaonMCTOFtag[iBin][part.first]->GetParameter(iPar));
                    fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][kAll]->SetParameter(part.second * 3 + iPar, fNsigmaTPCProtonMCV0tag[iBin][part.first]->GetParameter(iPar));

                    if (part.second * 3 + iPar != pdgPosition[kPion] * 3 + 1)
                    {
                        fNsigmaTPCPionDataV0tag[iEtaBin][iBin][kAll]->SetParLimits(part.second * 3 + iPar, fNsigmaTPCPionMCV0tag[iBin][part.first]->GetParameter(iPar) - TMath::Abs(fNsigmaTPCPionMCV0tag[iBin][part.first]->GetParameter(iPar)) * toll, fNsigmaTPCPionMCV0tag[iBin][part.first]->GetParameter(iPar) + TMath::Abs(fNsigmaTPCPionMCV0tag[iBin][part.first]->GetParameter(iPar)) * toll);
                    }
                    if (part.second * 3 + iPar != pdgPosition[kKaon] * 3 + 1)
                    {
                        fNsigmaTPCKaonDataKinktag[iEtaBin][iBin][kAll]->SetParLimits(part.second * 3 + iPar, fNsigmaTPCKaonMCKinktag[iBin][part.first]->GetParameter(iPar) - TMath::Abs(fNsigmaTPCKaonMCKinktag[iBin][part.first]->GetParameter(iPar)) * toll, fNsigmaTPCKaonMCKinktag[iBin][part.first]->GetParameter(iPar) + TMath::Abs(fNsigmaTPCKaonMCKinktag[iBin][part.first]->GetParameter(iPar)) * toll);

                        fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][kAll]->SetParLimits(part.second * 3 + iPar, fNsigmaTPCKaonMCTOFtag[iBin][part.first]->GetParameter(iPar) - TMath::Abs(fNsigmaTPCKaonMCTOFtag[iBin][part.first]->GetParameter(iPar)) * toll, fNsigmaTPCKaonMCTOFtag[iBin][part.first]->GetParameter(iPar) + TMath::Abs(fNsigmaTPCKaonMCTOFtag[iBin][part.first]->GetParameter(iPar)) * toll);
                    }
                    if (part.second * 3 + iPar != pdgPosition[kPr] * 3 + 1)
                    {
                        fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][kAll]->SetParLimits(part.second * 3 + iPar, fNsigmaTPCProtonMCV0tag[iBin][part.first]->GetParameter(iPar) - TMath::Abs(fNsigmaTPCProtonMCV0tag[iBin][part.first]->GetParameter(iPar)) * toll, fNsigmaTPCProtonMCV0tag[iBin][part.first]->GetParameter(iPar) + TMath::Abs(fNsigmaTPCProtonMCV0tag[iBin][part.first]->GetParameter(iPar)) * toll);
                    }
                }
            }

            hNsigmaTPCPionDataV0tag[iEtaBin][iBin]->Fit(fNsigmaTPCPionDataV0tag[iEtaBin][iBin][kAll], fitopt.Data());
            hNsigmaTPCKaonDataKinktag[iEtaBin][iBin]->Fit(fNsigmaTPCKaonDataKinktag[iEtaBin][iBin][kAll], fitopt.Data());
            hNsigmaTPCKaonDataTOFtag[iEtaBin][iBin]->Fit(fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][kAll], fitopt.Data());
            hNsigmaTPCProtonDataV0tag[iEtaBin][iBin]->Fit(fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][kAll], fitopt.Data());

            for(auto &part : pdgPosition)
            {
                if(part.first == kAll)
                    continue;
                fNsigmaTPCPionDataV0tag[iEtaBin][iBin][part.first] = new TF1(Form("fNsigmaTPCPionDataV0tag_%s_%s_%s", pdgNames[part.first].data(), etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()), "gaus", -50, 50);
                fNsigmaTPCKaonDataKinktag[iEtaBin][iBin][part.first] = new TF1(Form("fNsigmaTPCKaonDataKinktag_%s_%s_%s", pdgNames[part.first].data(), etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()), "gaus", -50, 50);
                fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][part.first] = new TF1(Form("fNsigmaTPCKaonDataTOFtag_%s_%s_%s", pdgNames[part.first].data(), etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()), "gaus", -50, 50);
                fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][part.first] = new TF1(Form("fNsigmaTPCProtonDataV0tag_%s_%s_%s", pdgNames[part.first].data(), etaBinLabels[iEtaBin].Data(), binLabels[iBin].Data()), "gaus", -50, 50);

                for (int iPar = 0; iPar < 3; iPar++)
                {
                    fNsigmaTPCPionDataV0tag[iEtaBin][iBin][part.first]->SetParameter(iPar, fNsigmaTPCPionDataV0tag[iEtaBin][iBin][kAll]->GetParameter(part.second * 3 + iPar));
                    fNsigmaTPCKaonDataKinktag[iEtaBin][iBin][part.first]->SetParameter(iPar, fNsigmaTPCKaonDataKinktag[iEtaBin][iBin][kAll]->GetParameter(part.second * 3 + iPar));
                    fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][part.first]->SetParameter(iPar, fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][kAll]->GetParameter(part.second * 3 + iPar));
                    fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][part.first]->SetParameter(iPar, fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][kAll]->GetParameter(part.second * 3 + iPar));
                    fNsigmaTPCPionDataV0tag[iEtaBin][iBin][part.first]->SetParError(iPar, fNsigmaTPCPionDataV0tag[iEtaBin][iBin][kAll]->GetParError(part.second * 3 + iPar));
                    fNsigmaTPCKaonDataKinktag[iEtaBin][iBin][part.first]->SetParError(iPar, fNsigmaTPCKaonDataKinktag[iEtaBin][iBin][kAll]->GetParError(part.second * 3 + iPar));
                    fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][part.first]->SetParError(iPar, fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][kAll]->GetParError(part.second * 3 + iPar));
                    fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][part.first]->SetParError(iPar, fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][kAll]->GetParError(part.second * 3 + iPar));
                }
                fNsigmaTPCPionDataV0tag[iEtaBin][iBin][part.first]->SetLineColor(pdgColors[part.first]);
                fNsigmaTPCKaonDataKinktag[iEtaBin][iBin][part.first]->SetLineColor(pdgColors[part.first]);
                fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][part.first]->SetLineColor(pdgColors[part.first]);
                fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][part.first]->SetLineColor(pdgColors[part.first]);
                if (iBin == 0 && iEtaBin == 0)
                    legTPCFitter->AddEntry(fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][part.first], Form("Func %s", pdgNames[part.first].data()), "l");
            }

            cPionDataV0tagTPC[iEtaBin]->cd(iBin+1)->SetLogy();
            hNsigmaTPCPionDataV0tag[iEtaBin][iBin]->GetXaxis()->SetNdivisions(505);
            hNsigmaTPCPionDataV0tag[iEtaBin][iBin]->DrawCopy("E");
            for(auto &part : pdgNames)
            {
                if(part.first == kAll)
                    continue;
                fNsigmaTPCPionDataV0tag[iEtaBin][iBin][part.first]->DrawCopy("same");
            }
            legTPCFitter->Draw();

            cKaonDataKinkstagTPC[iEtaBin]->cd(iBin+1)->SetLogy();
            hNsigmaTPCKaonDataKinktag[iEtaBin][iBin]->GetXaxis()->SetNdivisions(505);
            hNsigmaTPCKaonDataKinktag[iEtaBin][iBin]->DrawCopy("E");
            for(auto &part : pdgNames)
            {
                if(part.first == kAll)
                    continue;
                fNsigmaTPCKaonDataKinktag[iEtaBin][iBin][part.first]->DrawCopy("same");
            }
            legTPCFitter->Draw();

            cKaonDataTOFtagTPC[iEtaBin]->cd(iBin+1)->SetLogy();
            hNsigmaTPCKaonDataTOFtag[iEtaBin][iBin]->GetXaxis()->SetNdivisions(505);
            hNsigmaTPCKaonDataTOFtag[iEtaBin][iBin]->DrawCopy("E");
            for(auto &part : pdgNames)
            {
                if(part.first == kAll)
                    continue;
                fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][part.first]->DrawCopy("same");
            }
            legTPCFitter->Draw();

            cProtonDataV0tagTPC[iEtaBin]->cd(iBin+1)->SetLogy();
            hNsigmaTPCProtonDataV0tag[iEtaBin][iBin]->GetXaxis()->SetNdivisions(505);
            hNsigmaTPCProtonDataV0tag[iEtaBin][iBin]->DrawCopy("E");
            for(auto &part : pdgNames)
            {
                if(part.first == kAll)
                    continue;
                fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][part.first]->DrawCopy("same");
            }
            legTPCFitter->Draw();

            //subtract bkg template for TOF nsigma
            TString name = hNsigmaTOFPionDataV0tag[iEtaBin][iBin]->GetName();
            name.ReplaceAll("hNsigmaTOFPionDataV0tag", "hNsigmaTOFPionDataV0tagSub");
            hNsigmaTOFPionDataV0tagSub[iEtaBin][iBin] = static_cast<TH1D*>(hNsigmaTOFPionDataV0tag[iEtaBin][iBin]->Clone(name.Data()));

            name = hNsigmaTOFKaonDataKinktag[iEtaBin][iBin]->GetName();
            name.ReplaceAll("hNsigmaTOFKaonDataKinktag", "hNsigmaTOFKaonDataKinktagSub");
            hNsigmaTOFKaonDataKinktagSub[iEtaBin][iBin] = static_cast<TH1D*>(hNsigmaTOFKaonDataKinktag[iEtaBin][iBin]->Clone(name.Data()));

            name = hNsigmaTOFKaonDataTPCtag[iEtaBin][iBin]->GetName();
            name.ReplaceAll("hNsigmaTOFKaonDataTPCtag", "hNsigmaTOFKaonDataTPCtagSub");
            hNsigmaTOFKaonDataTPCtagSub[iEtaBin][iBin] = static_cast<TH1D*>(hNsigmaTOFKaonDataTPCtag[iEtaBin][iBin]->Clone(name.Data()));

            name = hNsigmaTOFProtonDataV0tag[iEtaBin][iBin]->GetName();
            name.ReplaceAll("hNsigmaTOFProtonDataV0tag", "hNsigmaTOFProtonDataV0tagSub");
            hNsigmaTOFProtonDataV0tagSub[iEtaBin][iBin] = static_cast<TH1D*>(hNsigmaTOFProtonDataV0tag[iEtaBin][iBin]->Clone(name.Data()));

            for(auto &part : pdgNames)
            {
                if(part.first == kAll)
                    continue;

                if (part.first != kPion)
                    hNsigmaTOFPionDataV0tagSub[iEtaBin][iBin]->Add(hNsigmaTOFPionMCV0tag[iBin][part.first], -1.);
                if (part.first != kKaon)
                {
                    hNsigmaTOFKaonDataKinktagSub[iEtaBin][iBin]->Add(hNsigmaTOFKaonMCKinktag[iBin][part.first], -1.);
                    if (hFracTOFKaonDataTPCtag[iEtaBin][part.first]->GetBinContent(iBin+1) > 0)
                        hNsigmaTOFKaonDataTPCtagSub[iEtaBin][iBin]->Add(hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][part.first], -1.);
                    else
                        hNsigmaTOFKaonDataTPCtagSub[iEtaBin][iBin]->Add(hNsigmaTOFKaonMCTPCtag[iBin][part.first], -1.);
                }
                if (part.first != kPr)
                {
                    if (hFracTOFProtonDataV0tag[iEtaBin][part.first]->GetBinContent(iBin+1) > 0)
                        hNsigmaTOFProtonDataV0tagSub[iEtaBin][iBin]->Add(hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][part.first], -1.);
                    else
                        hNsigmaTOFProtonDataV0tagSub[iEtaBin][iBin]->Add(hNsigmaTOFProtonMCV0tag[iBin][part.first], -1.);
                }
            }

            SetTH1Style(hNsigmaTOFPionDataV0tagSub[iEtaBin][iBin], kOpenCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], pdgFillColors[kAll], 0.055, 0.06);
            SetTH1Style(hNsigmaTOFKaonDataKinktagSub[iEtaBin][iBin], kOpenCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], pdgFillColors[kAll], 0.055, 0.06);
            SetTH1Style(hNsigmaTOFKaonDataTPCtagSub[iEtaBin][iBin], kOpenCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], pdgFillColors[kAll], 0.055, 0.06);
            SetTH1Style(hNsigmaTOFProtonDataV0tagSub[iEtaBin][iBin], kOpenCircle, pdgColors[kAll], 0.6, 2, pdgColors[kAll], pdgFillColors[kAll], 0.055, 0.06);
            if (iBin == 0 && iEtaBin == 0)
            {
                legTOFPionDataV0tag->AddEntry(hNsigmaTOFPionDataV0tagSub[iEtaBin][iBin], "Data", "p");
                legTOFPionDataV0tag->AddEntry(hNsigmaTOFPionMCV0tag[iBin][kPion], "MC Pion", "p");
                legTOFKaonDataKinkstag->AddEntry(hNsigmaTOFKaonDataKinktagSub[iEtaBin][iBin], "Data", "p");
                legTOFKaonDataKinkstag->AddEntry(hNsigmaTOFKaonMCKinktag[iBin][kKaon], "MC Kaon", "p");
                legTOFKaonDataTPCtag->AddEntry(hNsigmaTOFKaonDataTPCtagSub[iEtaBin][iBin], "Data", "p");
                legTOFKaonDataTPCtag->AddEntry(hNsigmaTOFKaonMCTPCtag[iBin][kKaon], "MC Kaon", "p");
                legTOFProtonDataV0tag->AddEntry(hNsigmaTOFProtonDataV0tagSub[iEtaBin][iBin], "Data", "p");
                legTOFProtonDataV0tag->AddEntry(hNsigmaTOFKaonMCTPCtag[iBin][kPr], "MC Proton", "p");
            }

            cPionDataV0tagTOF[iEtaBin]->cd(iBin+1)->SetLogy();
            hNsigmaTOFPionDataV0tagSub[iEtaBin][iBin]->GetXaxis()->SetNdivisions(505);
            hNsigmaTOFPionDataV0tagSub[iEtaBin][iBin]->DrawCopy("E");
            hNsigmaTOFPionMCV0tag[iBin][kPion]->DrawCopy("hist same");
            legTOFPionDataV0tag->Draw();
            cKaonDataKinkstagTOF[iEtaBin]->cd(iBin+1)->SetLogy();
            hNsigmaTOFKaonDataKinktagSub[iEtaBin][iBin]->GetXaxis()->SetNdivisions(505);
            hNsigmaTOFKaonDataKinktagSub[iEtaBin][iBin]->DrawCopy("E");
            hNsigmaTOFKaonMCKinktag[iBin][kKaon]->DrawCopy("hist same");
            legTOFKaonDataKinkstag->Draw();
            cKaonDataTPCtagTOF[iEtaBin]->cd(iBin+1)->SetLogy();
            hNsigmaTOFKaonDataTPCtagSub[iEtaBin][iBin]->GetXaxis()->SetNdivisions(505);
            hNsigmaTOFKaonDataTPCtagSub[iEtaBin][iBin]->DrawCopy("E");
            hNsigmaTOFKaonMCTPCtag[iBin][kKaon]->DrawCopy("hist same");
            legTOFKaonDataTPCtag->Draw();
            cProtonDataV0tagTOF[iEtaBin]->cd(iBin+1)->SetLogy();
            hNsigmaTOFProtonDataV0tagSub[iEtaBin][iBin]->GetXaxis()->SetNdivisions(505);
            hNsigmaTOFProtonDataV0tagSub[iEtaBin][iBin]->DrawCopy("E");
            hNsigmaTOFProtonMCV0tag[iBin][kPr]->DrawCopy("hist same");
            legTOFProtonDataV0tag->Draw();
        }
        cPionDataV0tagTPC[iEtaBin]->SaveAs(Form("%s/DistrNsigmaTPCPionFromV0tag_Data_%s.pdf", outDirName.data(), etaBinLabels[iEtaBin].Data()));
        cKaonDataTOFtagTPC[iEtaBin]->SaveAs(Form("%s/DistrNsigmaTPCKaonFromKinkstag_Data_%s.pdf", outDirName.data(), etaBinLabels[iEtaBin].Data()));
        cKaonDataTOFtagTPC[iEtaBin]->SaveAs(Form("%s/DistrNsigmaTPCKaonFromTOFtag_Data_%s.pdf", outDirName.data(), etaBinLabels[iEtaBin].Data()));
        cProtonDataV0tagTPC[iEtaBin]->SaveAs(Form("%s/DistrNsigmaTPCProtonFromV0tag_Data_%s.pdf", outDirName.data(), etaBinLabels[iEtaBin].Data()));
        cPionDataV0tagTOF[iEtaBin]->SaveAs(Form("%s/DistrNsigmaTOFPionFromV0tag_Data_%s.pdf", outDirName.data(), etaBinLabels[iEtaBin].Data()));
        cKaonDataKinkstagTOF[iEtaBin]->SaveAs(Form("%s/DistrNsigmaTOFKaonFromKinkstag_Data_%s.pdf", outDirName.data(), etaBinLabels[iEtaBin].Data()));
        cKaonDataTPCtagTOF[iEtaBin]->SaveAs(Form("%s/DistrNsigmaTOFKaonFromTPCtag_Data_%s.pdf", outDirName.data(), etaBinLabels[iEtaBin].Data()));
        cProtonDataV0tagTOF[iEtaBin]->SaveAs(Form("%s/DistrNsigmaTOFProtonFromV0tag_Data_%s.pdf", outDirName.data(), etaBinLabels[iEtaBin].Data()));
    }

    std::cout << "\033[32mDone\033[0m" << std::endl;

    //compute efficiencies
    std::cout << "\n*******************************************\n" << std::endl;
    std::cout << "\033[32mCompute efficiencies\033[0m\n" << std::endl;

    std::array<TH1D*, nMaxEff> hEffPionTPCMCtrue, hEffPionTOFMCtrue, hEffKaonTPCMCtrue, hEffKaonTOFMCtrue, hEffProtonTPCMCtrue, hEffProtonTOFMCtrue;
    std::array<std::array<TH1D*, nMaxEff>, nEtaBinsMax+1> hEffPionTPCDataV0tag, hEffPionTOFDataV0tag, hEffKaonTPCDataKinktag, hEffKaonTOFDataKinktag, hEffKaonTPCDataTOFtag, hEffKaonTOFDataTPCtag, hEffProtonTPCDataV0tag, hEffProtonTOFDataV0tag;
    std::array<std::array<TH1D*, nMaxEff>, nEtaBinsMax+1> hRatioEffPionTPCDataV0tag, hRatioEffPionTOFDataV0tag, hRatioEffKaonTPCDataKinktag, hRatioEffKaonTOFDataKinktag, hRatioEffKaonTPCDataTOFtag, hRatioEffKaonTOFDataTPCtag, hRatioEffProtonTPCDataV0tag, hRatioEffProtonTOFDataV0tag;

    std::array<TCanvas*, nEtaBinsMax+1> cEffPion, cEffKaon, cEffProton;

    TLegend *legEffPion = new TLegend(0.18, 0.15, 0.925, 0.38);
    legEffPion->SetTextSize(0.045);
    legEffPion->SetNColumns(2);
    legEffPion->AddEntry("", "MC", "");
    legEffPion->AddEntry("", "V0 tag", "");
    TLegend *legEffKaonTPC = new TLegend(0.18, 0.15, 0.925, 0.38);
    legEffKaonTPC->SetTextSize(0.045);
    legEffKaonTPC->SetNColumns(3);
    legEffKaonTPC->AddEntry("", "MC", "");
    legEffKaonTPC->AddEntry("", "TOF tag", "");
    legEffKaonTPC->AddEntry("", "kink tag", "");
    TLegend *legEffKaonTOF = new TLegend(0.18, 0.15, 0.925, 0.38);
    legEffKaonTOF->SetTextSize(0.045);
    legEffKaonTOF->SetNColumns(2);
    legEffKaonTOF->AddEntry("", "MC", "");
    legEffKaonTOF->AddEntry("", "TPC tag", "");
    TLegend *legEffProton = new TLegend(0.18, 0.15, 0.925, 0.38);
    legEffProton->SetTextSize(0.045);
    legEffProton->SetNColumns(2);
    legEffProton->AddEntry("", "MC", "");
    legEffProton->AddEntry("", "V0 tag", "");

    for (unsigned int iEff = 0; iEff < nEff; iEff++)
    {
        int nSigmaBinLow = hNsigmaTOFPionMCTrue[0]->GetXaxis()->FindBin(-nSigma4Eff[iEff] + 1.e-5);
        int nSigmaBinHigh = hNsigmaTOFPionMCTrue[0]->GetXaxis()->FindBin(nSigma4Eff[iEff] - 1.e-5);

        hEffPionTPCMCtrue[iEff] = new TH1D(Form("hEffPionTPCMCtrue_%dsigma", nSigma4Eff[iEff]), Form(";%s (GeV/#it{c});TPC #pi efficiency", varTitle.Data()), static_cast<int>(nBins), binLims);
        hEffPionTOFMCtrue[iEff] = new TH1D(Form("hEffPionTOFMCtrue_%dsigma", nSigma4Eff[iEff]), Form(";%s (GeV/#it{c});TOF #pi efficiency", varTitle.Data()), static_cast<int>(nBins), binLims);

        hEffKaonTPCMCtrue[iEff] = new TH1D(Form("hEffKaonTPCMCtrue_%dsigma", nSigma4Eff[iEff]), Form(";%s (GeV/#it{c});TPC K efficiency", varTitle.Data()), static_cast<int>(nBins), binLims);
        hEffKaonTOFMCtrue[iEff] = new TH1D(Form("hEffKaonTOFMCtrue_%dsigma", nSigma4Eff[iEff]), Form(";%s (GeV/#it{c});TOF K efficiency", varTitle.Data()), static_cast<int>(nBins), binLims);

        hEffProtonTPCMCtrue[iEff] = new TH1D(Form("hEffProtonTPCMCtrue_%dsigma", nSigma4Eff[iEff]), Form(";%s (GeV/#it{c});TPC p efficiency", varTitle.Data()), static_cast<int>(nBins), binLims);
        hEffProtonTOFMCtrue[iEff] = new TH1D(Form("hEffProtonTOFMCtrue_%dsigma", nSigma4Eff[iEff]), Form(";%s (GeV/#it{c});TOF p efficiency", varTitle.Data()), static_cast<int>(nBins), binLims);

        SetTH1Style(hEffPionTPCMCtrue[iEff], markersEffMC[iEff], pdgColors[kPion]-1, 1., 2, pdgColors[kPion]-1, kWhite, 0.045, 0.055);
        SetTH1Style(hEffPionTOFMCtrue[iEff], markersEffMC[iEff], pdgColors[kPion]-1, 1., 2, pdgColors[kPion]-1, kWhite, 0.045, 0.055);
        SetTH1Style(hEffKaonTPCMCtrue[iEff], markersEffMC[iEff], pdgColors[kKaon]-3, 1., 2, pdgColors[kKaon], kWhite, 0.045, 0.055);
        SetTH1Style(hEffKaonTOFMCtrue[iEff], markersEffMC[iEff], pdgColors[kKaon]-3, 1., 2, pdgColors[kKaon], kWhite, 0.045, 0.055);
        SetTH1Style(hEffProtonTPCMCtrue[iEff], markersEffMC[iEff], pdgColors[kPr], 1., 2, pdgColors[kPr], kWhite, 0.045, 0.055);
        SetTH1Style(hEffProtonTOFMCtrue[iEff], markersEffMC[iEff], pdgColors[kPr], 1., 2, pdgColors[kPr], kWhite, 0.045, 0.055);

        double eff = -1, unc = -1;
        for (unsigned int iBin = 0; iBin < nBins; iBin++)
        {
            ComputeEfficiency(hNsigmaTPCPionMCTrue[iBin]->Integral(nSigmaBinLow, nSigmaBinHigh), hNsigmaTPCPionMCTrue[iBin]->Integral(), eff, unc);
            hEffPionTPCMCtrue[iEff]->SetBinContent(iBin+1, eff);
            hEffPionTPCMCtrue[iEff]->SetBinError(iBin+1, unc);

            ComputeEfficiency(hNsigmaTOFPionMCTrue[iBin]->Integral(nSigmaBinLow, nSigmaBinHigh), hNsigmaTOFPionMCTrue[iBin]->Integral(), eff, unc);
            hEffPionTOFMCtrue[iEff]->SetBinContent(iBin+1, eff);
            hEffPionTOFMCtrue[iEff]->SetBinError(iBin+1, unc);

            ComputeEfficiency(hNsigmaTPCKaonMCTrue[iBin]->Integral(nSigmaBinLow, nSigmaBinHigh), hNsigmaTPCKaonMCTrue[iBin]->Integral(), eff, unc);
            hEffKaonTPCMCtrue[iEff]->SetBinContent(iBin+1, eff);
            hEffKaonTPCMCtrue[iEff]->SetBinError(iBin+1, unc);

            ComputeEfficiency(hNsigmaTOFKaonMCTrue[iBin]->Integral(nSigmaBinLow, nSigmaBinHigh), hNsigmaTOFKaonMCTrue[iBin]->Integral(), eff, unc);
            hEffKaonTOFMCtrue[iEff]->SetBinContent(iBin+1, eff);
            hEffKaonTOFMCtrue[iEff]->SetBinError(iBin+1, unc);

            ComputeEfficiency(hNsigmaTPCProtonMCTrue[iBin]->Integral(nSigmaBinLow, nSigmaBinHigh), hNsigmaTPCProtonMCTrue[iBin]->Integral(), eff, unc);
            hEffProtonTPCMCtrue[iEff]->SetBinContent(iBin+1, eff);
            hEffProtonTPCMCtrue[iEff]->SetBinError(iBin+1, unc);

            ComputeEfficiency(hNsigmaTOFProtonMCTrue[iBin]->Integral(nSigmaBinLow, nSigmaBinHigh), hNsigmaTOFProtonMCTrue[iBin]->Integral(), eff, unc);
            hEffProtonTOFMCtrue[iEff]->SetBinContent(iBin+1, eff);
            hEffProtonTOFMCtrue[iEff]->SetBinError(iBin+1, unc);
        }

        for(unsigned int iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++)
        {
            hEffPionTPCDataV0tag[iEtaBin][iEff] = new TH1D(Form("hEffPionTPCDataV0tag_%s_%dsigma", etaBinLabels[iEtaBin].Data(), nSigma4Eff[iEff]), Form(";%s (GeV/#it{c});TPC #pi efficiency", varTitle.Data()), static_cast<int>(nBins), binLims);
            hEffPionTOFDataV0tag[iEtaBin][iEff] = new TH1D(Form("hEffPionTOFDataV0tag_%s_%dsigma", etaBinLabels[iEtaBin].Data(), nSigma4Eff[iEff]), Form(";%s (GeV/#it{c});TOF #pi efficiency", varTitle.Data()), static_cast<int>(nBins), binLims);

            hEffKaonTPCDataKinktag[iEtaBin][iEff] = new TH1D(Form("hEffKaonTPCDataKinktag_%s_%dsigma", etaBinLabels[iEtaBin].Data(), nSigma4Eff[iEff]), Form(";%s (GeV/#it{c});TPC K efficiency", varTitle.Data()), static_cast<int>(nBins), binLims);
            hEffKaonTOFDataKinktag[iEtaBin][iEff] = new TH1D(Form("hEffKaonTOFDataKinktag_%s_%dsigma", etaBinLabels[iEtaBin].Data(), nSigma4Eff[iEff]), Form(";%s (GeV/#it{c});TOF K efficiency", varTitle.Data()), static_cast<int>(nBins), binLims);

            hEffKaonTPCDataTOFtag[iEtaBin][iEff] = new TH1D(Form("hEffKaonTPCDataTOFtag_%s_%dsigma", etaBinLabels[iEtaBin].Data(), nSigma4Eff[iEff]), Form(";%s (GeV/#it{c});TPC K efficiency", varTitle.Data()), static_cast<int>(nBins), binLims);
            hEffKaonTOFDataTPCtag[iEtaBin][iEff] = new TH1D(Form("hEffKaonTOFDataTPCtag_%s_%dsigma", etaBinLabels[iEtaBin].Data(), nSigma4Eff[iEff]), Form(";%s (GeV/#it{c});TOF K efficiency", varTitle.Data()), static_cast<int>(nBins), binLims);

            hEffProtonTPCDataV0tag[iEtaBin][iEff] = new TH1D(Form("hEffProtonTPCDataV0tag_%s_%dsigma", etaBinLabels[iEtaBin].Data(), nSigma4Eff[iEff]), Form(";%s (GeV/#it{c});TPC p efficiency", varTitle.Data()), static_cast<int>(nBins), binLims);
            hEffProtonTOFDataV0tag[iEtaBin][iEff] = new TH1D(Form("hEffProtonTOFDataV0tag_%s_%dsigma", etaBinLabels[iEtaBin].Data(), nSigma4Eff[iEff]), Form(";%s (GeV/#it{c});TOF p efficiency", varTitle.Data()), static_cast<int>(nBins), binLims);

            SetTH1Style(hEffPionTPCDataV0tag[iEtaBin][iEff], markersEffData[iEff], pdgColors[kPion], 1., 2, pdgColors[kPion], kWhite, 0.045, 0.055);
            SetTH1Style(hEffPionTOFDataV0tag[iEtaBin][iEff], markersEffData[iEff], pdgColors[kPion], 1., 2, pdgColors[kPion], kWhite, 0.045, 0.055);
            SetTH1Style(hEffKaonTPCDataKinktag[iEtaBin][iEff], markersEffData[iEff], pdgColors[kKaon]-1, 1., 2, pdgColors[kKaon]-1, kWhite, 0.045, 0.055);
            SetTH1Style(hEffKaonTPCDataTOFtag[iEtaBin][iEff], markersEffData[iEff], pdgColors[kKaon], 1., 2, pdgColors[kKaon], kWhite, 0.045, 0.055);
            SetTH1Style(hEffKaonTOFDataTPCtag[iEtaBin][iEff], markersEffData[iEff], pdgColors[kKaon], 1., 2, pdgColors[kKaon], kWhite, 0.045, 0.055);
            SetTH1Style(hEffProtonTPCDataV0tag[iEtaBin][iEff], markersEffData[iEff], pdgColors[kPr]+1, 1., 2, pdgColors[kPr]+1, kWhite, 0.045, 0.055);
            SetTH1Style(hEffProtonTOFDataV0tag[iEtaBin][iEff], markersEffData[iEff], pdgColors[kPr]+1, 1., 2, pdgColors[kPr]+1, kWhite, 0.045, 0.055);

            hRatioEffPionTPCDataV0tag[iEtaBin][iEff] = static_cast<TH1D*>(hEffPionTPCDataV0tag[iEtaBin][iEff]->Clone(Form("hRatioEffPionTPCDataV0tag_%s_%dsigma", etaBinLabels[iEtaBin].Data(), nSigma4Eff[iEff])));
            hRatioEffPionTPCDataV0tag[iEtaBin][iEff]->GetYaxis()->SetTitle("TPC #pi efficiency ratio data/MC");
            hRatioEffPionTOFDataV0tag[iEtaBin][iEff] = static_cast<TH1D*>(hEffPionTOFDataV0tag[iEtaBin][iEff]->Clone(Form("hRatioEffPionTOFDataV0tag_%s_%dsigma", etaBinLabels[iEtaBin].Data(), nSigma4Eff[iEff])));
            hRatioEffPionTPCDataV0tag[iEtaBin][iEff]->GetYaxis()->SetTitle("TOF #pi efficiency ratio data/MC");
            hRatioEffKaonTPCDataKinktag[iEtaBin][iEff] = static_cast<TH1D*>(hEffKaonTPCDataKinktag[iEtaBin][iEff]->Clone(Form("hRatioEffKaonTPCDataKinktag_%s_%dsigma", etaBinLabels[iEtaBin].Data(), nSigma4Eff[iEff])));
            hRatioEffPionTPCDataV0tag[iEtaBin][iEff]->GetYaxis()->SetTitle("TPC K efficiency ratio data/MC");
            hRatioEffKaonTOFDataKinktag[iEtaBin][iEff] = static_cast<TH1D*>(hEffKaonTOFDataKinktag[iEtaBin][iEff]->Clone(Form("hRatioEffKaonTOFDataKinktag_%s_%dsigma", etaBinLabels[iEtaBin].Data(), nSigma4Eff[iEff])));
            hRatioEffPionTPCDataV0tag[iEtaBin][iEff]->GetYaxis()->SetTitle("TOF K efficiency ratio data/MC");
            hRatioEffKaonTPCDataTOFtag[iEtaBin][iEff] = static_cast<TH1D*>(hEffKaonTPCDataTOFtag[iEtaBin][iEff]->Clone(Form("hRatioEffKaonTPCDataTOFtag_%s_%dsigma", etaBinLabels[iEtaBin].Data(), nSigma4Eff[iEff])));
            hRatioEffPionTPCDataV0tag[iEtaBin][iEff]->GetYaxis()->SetTitle("TPC K efficiency ratio data/MC");
            hRatioEffKaonTOFDataTPCtag[iEtaBin][iEff] = static_cast<TH1D*>(hEffKaonTOFDataTPCtag[iEtaBin][iEff]->Clone(Form("hRatioEffKaonTOFDataTPCtag_%s_%dsigma", etaBinLabels[iEtaBin].Data(), nSigma4Eff[iEff])));
            hRatioEffPionTPCDataV0tag[iEtaBin][iEff]->GetYaxis()->SetTitle("TOF K efficiency ratio data/MC");
            hRatioEffProtonTPCDataV0tag[iEtaBin][iEff] = static_cast<TH1D*>(hEffProtonTPCDataV0tag[iEtaBin][iEff]->Clone(Form("hRatioEffProtonTPCDataV0tag_%s_%dsigma", etaBinLabels[iEtaBin].Data(), nSigma4Eff[iEff])));
            hRatioEffPionTPCDataV0tag[iEtaBin][iEff]->GetYaxis()->SetTitle("TPC p efficiency ratio data/MC");
            hRatioEffProtonTOFDataV0tag[iEtaBin][iEff] = static_cast<TH1D*>(hEffProtonTOFDataV0tag[iEtaBin][iEff]->Clone(Form("hRatioEffProtonTOFDataV0tag_%s_%dsigma", etaBinLabels[iEtaBin].Data(), nSigma4Eff[iEff])));
            hRatioEffPionTPCDataV0tag[iEtaBin][iEff]->GetYaxis()->SetTitle("TOF p efficiency ratio data/MC");

            for (unsigned int iBin = 0; iBin < nBins; iBin++)
            {
                ComputeEfficiency(fNsigmaTPCPionDataV0tag[iEtaBin][iBin][kPion]->Integral(-nSigma4Eff[iEff], nSigma4Eff[iEff]) * intNsTPCPionDataV0tag[iEtaBin][iBin], fNsigmaTPCPionDataV0tag[iEtaBin][iBin][kPion]->Integral(-50, 50) * intNsTPCPionDataV0tag[iEtaBin][iBin], eff, unc);
                hEffPionTPCDataV0tag[iEtaBin][iEff]->SetBinContent(iBin+1, eff);
                hEffPionTPCDataV0tag[iEtaBin][iEff]->SetBinError(iBin+1, unc);

                ComputeEfficiency(fNsigmaTPCKaonDataKinktag[iEtaBin][iBin][kKaon]->Integral(-nSigma4Eff[iEff], nSigma4Eff[iEff]) * intNsTPCKaonDataKinktag[iEtaBin][iBin], fNsigmaTPCKaonDataKinktag[iEtaBin][iBin][kKaon]->Integral(-50, 50) * intNsTPCKaonDataKinktag[iEtaBin][iBin], eff, unc);
                hEffKaonTPCDataKinktag[iEtaBin][iEff]->SetBinContent(iBin+1, eff);
                hEffKaonTPCDataKinktag[iEtaBin][iEff]->SetBinError(iBin+1, unc);

                ComputeEfficiency(fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][kKaon]->Integral(-nSigma4Eff[iEff], nSigma4Eff[iEff]) * intNsTPCKaonDataTOFtag[iEtaBin][iBin], fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][kKaon]->Integral(-50, 50) * intNsTPCKaonDataTOFtag[iEtaBin][iBin], eff, unc);
                hEffKaonTPCDataTOFtag[iEtaBin][iEff]->SetBinContent(iBin+1, eff);
                hEffKaonTPCDataTOFtag[iEtaBin][iEff]->SetBinError(iBin+1, unc);

                ComputeEfficiency(fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][kPr]->Integral(-nSigma4Eff[iEff], nSigma4Eff[iEff]) * intNsTPCProtonDataV0tag[iEtaBin][iBin], fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][kPr]->Integral(-50, 50) * intNsTPCProtonDataV0tag[iEtaBin][iBin], eff, unc);
                hEffProtonTPCDataV0tag[iEtaBin][iEff]->SetBinContent(iBin+1, eff);
                hEffProtonTPCDataV0tag[iEtaBin][iEff]->SetBinError(iBin+1, unc);

                ComputeEfficiency(hNsigmaTOFPionDataV0tagSub[iEtaBin][iBin]->Integral(nSigmaBinLow, nSigmaBinHigh) * intNsTOFPionDataV0tag[iEtaBin][iBin], hNsigmaTOFPionDataV0tagSub[iEtaBin][iBin]->Integral() * intNsTOFPionDataV0tag[iEtaBin][iBin], eff, unc);
                hEffPionTOFDataV0tag[iEtaBin][iEff]->SetBinContent(iBin+1, eff);
                hEffPionTOFDataV0tag[iEtaBin][iEff]->SetBinError(iBin+1, unc);

                ComputeEfficiency(hNsigmaTOFKaonDataTPCtagSub[iEtaBin][iBin]->Integral(nSigmaBinLow, nSigmaBinHigh) * intNsTOFKaonDataTPCtag[iEtaBin][iBin], hNsigmaTOFKaonDataTPCtagSub[iEtaBin][iBin]->Integral() * intNsTOFKaonDataTPCtag[iEtaBin][iBin], eff, unc);
                hEffKaonTOFDataTPCtag[iEtaBin][iEff]->SetBinContent(iBin+1, eff);
                hEffKaonTOFDataTPCtag[iEtaBin][iEff]->SetBinError(iBin+1, unc);

                ComputeEfficiency(hNsigmaTOFProtonDataV0tagSub[iEtaBin][iBin]->Integral(nSigmaBinLow, nSigmaBinHigh) * intNsTOFProtonDataV0tag[iEtaBin][iBin], hNsigmaTOFProtonDataV0tagSub[iEtaBin][iBin]->Integral() * intNsTOFProtonDataV0tag[iEtaBin][iBin], eff, unc);
                hEffProtonTOFDataV0tag[iEtaBin][iEff]->SetBinContent(iBin+1, eff);
                hEffProtonTOFDataV0tag[iEtaBin][iEff]->SetBinError(iBin+1, unc);
            }

            hRatioEffPionTPCDataV0tag[iEtaBin][iEff]->Divide(hEffPionTPCDataV0tag[iEtaBin][iEff], hEffPionTPCMCtrue[iEff], 1., 1., "");
            hRatioEffPionTOFDataV0tag[iEtaBin][iEff]->Divide(hEffPionTOFDataV0tag[iEtaBin][iEff], hEffPionTOFMCtrue[iEff], 1., 1., "");
            hRatioEffKaonTPCDataKinktag[iEtaBin][iEff]->Divide(hEffKaonTPCDataKinktag[iEtaBin][iEff], hEffKaonTPCMCtrue[iEff], 1., 1., "");
            hRatioEffKaonTPCDataTOFtag[iEtaBin][iEff]->Divide(hEffKaonTPCDataTOFtag[iEtaBin][iEff], hEffKaonTPCMCtrue[iEff], 1., 1., "");
            hRatioEffKaonTOFDataTPCtag[iEtaBin][iEff]->Divide(hEffKaonTOFDataTPCtag[iEtaBin][iEff], hEffKaonTOFMCtrue[iEff], 1., 1., "");
            hRatioEffProtonTPCDataV0tag[iEtaBin][iEff]->Divide(hEffProtonTPCDataV0tag[iEtaBin][iEff], hEffProtonTPCMCtrue[iEff], 1., 1., "");
            hRatioEffProtonTOFDataV0tag[iEtaBin][iEff]->Divide(hEffProtonTOFDataV0tag[iEtaBin][iEff], hEffProtonTOFMCtrue[iEff], 1., 1., "");

            if(iEtaBin == 0)
            {
                legEffPion->AddEntry(hEffPionTPCMCtrue[iEff], Form("|N_{#sigma}(#pi)| < %d", nSigma4Eff[iEff]), "p");
                legEffPion->AddEntry(hEffPionTPCDataV0tag[iEtaBin][iEff], Form("|N_{#sigma}(#pi)| < %d", nSigma4Eff[iEff]), "p");
                legEffKaonTPC->AddEntry(hEffKaonTPCMCtrue[iEff], Form("|N_{#sigma}(K)| < %d", nSigma4Eff[iEff]), "p");
                legEffKaonTPC->AddEntry(hEffKaonTPCDataTOFtag[iEtaBin][iEff], Form("|N_{#sigma}(K)| < %d", nSigma4Eff[iEff]), "p");
                legEffKaonTPC->AddEntry(hEffKaonTPCDataKinktag[iEtaBin][iEff], Form("|N_{#sigma}(K)| < %d", nSigma4Eff[iEff]), "p");
                legEffKaonTOF->AddEntry(hEffKaonTPCMCtrue[iEff], Form("|N_{#sigma}(K)| < %d", nSigma4Eff[iEff]), "p");
                legEffKaonTOF->AddEntry(hEffKaonTPCDataTOFtag[iEtaBin][iEff], Form("|N_{#sigma}(K)| < %d", nSigma4Eff[iEff]), "p");
                legEffProton->AddEntry(hEffProtonTPCMCtrue[iEff], Form("|N_{#sigma}(p)| < %d", nSigma4Eff[iEff]), "p");
                legEffProton->AddEntry(hEffProtonTOFDataV0tag[iEtaBin][iEff], Form("|N_{#sigma}(p)| < %d", nSigma4Eff[iEff]), "p");
            }
        }
    }

    for(unsigned int iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++)
    {
        cEffPion[iEtaBin] = new TCanvas(Form("cEffPion_%s", etaBinLabels[iEtaBin].Data()), Form("cEffPion_%s", etaBinLabels[iEtaBin].Data()), 800, 800);
        cEffKaon[iEtaBin] = new TCanvas(Form("cEffKaon_%s", etaBinLabels[iEtaBin].Data()), Form("cEffKaon_%s", etaBinLabels[iEtaBin].Data()), 800, 800);
        cEffProton[iEtaBin] = new TCanvas(Form("cEffProton_%s", etaBinLabels[iEtaBin].Data()), Form("cEffProton_%s", etaBinLabels[iEtaBin].Data()), 800, 800);
        cEffPion[iEtaBin]->Divide(2, 2);
        cEffPion[iEtaBin]->cd(1)->DrawFrame(binLims[0], 0., binLims[nBins], 1., Form(";%s (GeV/#it{c});pion TPC PID efficiency", varTitle.Data()));
        cEffPion[iEtaBin]->cd(2)->DrawFrame(binLims[0], 0., binLims[nBins], 1., Form(";%s (GeV/#it{c});pion TOF PID efficiency", varTitle.Data()));
        cEffPion[iEtaBin]->cd(3)->DrawFrame(binLims[0], 0.5, binLims[nBins], 1.15, Form(";%s (GeV/#it{c}); pion data / MC TPC PID efficiency", varTitle.Data()));
        cEffPion[iEtaBin]->cd(4)->DrawFrame(binLims[0], 0.5, binLims[nBins], 1.15, Form(";%s (GeV/#it{c});pion data / MC TOF PID efficiency", varTitle.Data()));
        cEffKaon[iEtaBin]->Divide(2, 2);
        cEffKaon[iEtaBin]->cd(1)->DrawFrame(binLims[0], 0., binLims[nBins], 1., Form(";%s (GeV/#it{c});kaon TPC PID efficiency", varTitle.Data()));
        cEffKaon[iEtaBin]->cd(2)->DrawFrame(binLims[0], 0., binLims[nBins], 1., Form(";%s (GeV/#it{c});kaon TOF PID efficiency", varTitle.Data()));
        cEffKaon[iEtaBin]->cd(3)->DrawFrame(binLims[0], 0.5, binLims[nBins], 1.15, Form(";%s (GeV/#it{c});kaon data / MC TPC PID efficiency", varTitle.Data()));
        cEffKaon[iEtaBin]->cd(4)->DrawFrame(binLims[0], 0.5, binLims[nBins], 1.15, Form(";%s (GeV/#it{c});kaon data / MC TOF PID efficiency", varTitle.Data()));
        cEffProton[iEtaBin]->Divide(2, 2);
        cEffProton[iEtaBin]->cd(1)->DrawFrame(binLims[0], 0., binLims[nBins], 1., Form(";%s (GeV/#it{c});proton TPC PID efficiency", varTitle.Data()));
        cEffProton[iEtaBin]->cd(2)->DrawFrame(binLims[0], 0., binLims[nBins], 1., Form(";%s (GeV/#it{c});proton TOF PID efficiency", varTitle.Data()));
        cEffProton[iEtaBin]->cd(3)->DrawFrame(binLims[0], 0.5, binLims[nBins], 1.15, Form(";%s (GeV/#it{c});proton data / MC TPC PID efficiency", varTitle.Data()));
        cEffProton[iEtaBin]->cd(4)->DrawFrame(binLims[0], 0.5, binLims[nBins], 1.15, Form(";%s (GeV/#it{c});proton data / MC TOF PID efficiency", varTitle.Data()));

        for (unsigned int iEff = 0; iEff < nEff; iEff++)
        {
            // cEffPion[iEtaBin]->cd(1)->SetLogx();
            cEffPion[iEtaBin]->cd(1)->SetTopMargin(0.035);
            cEffPion[iEtaBin]->cd(1)->SetBottomMargin(0.12);
            cEffPion[iEtaBin]->cd(1)->SetRightMargin(0.035);
            hEffPionTPCMCtrue[iEff]->DrawCopy("same");
            hEffPionTPCDataV0tag[iEtaBin][iEff]->DrawCopy("same");
            legEffPion->Draw();

            // cEffPion[iEtaBin]->cd(2)->SetLogx();
            cEffPion[iEtaBin]->cd(2)->SetTopMargin(0.035);
            cEffPion[iEtaBin]->cd(2)->SetBottomMargin(0.12);
            cEffPion[iEtaBin]->cd(2)->SetRightMargin(0.035);
            hEffPionTOFMCtrue[iEff]->DrawCopy("same");
            hEffPionTOFDataV0tag[iEtaBin][iEff]->DrawCopy("same");
            legEffPion->Draw();

            // cEffPion[iEtaBin]->cd(3)->SetLogx();
            cEffPion[iEtaBin]->cd(3)->SetTopMargin(0.035);
            cEffPion[iEtaBin]->cd(3)->SetBottomMargin(0.12);
            cEffPion[iEtaBin]->cd(3)->SetRightMargin(0.035);
            hRatioEffPionTPCDataV0tag[iEtaBin][iEff]->DrawCopy("same");

            // cEffPion[iEtaBin]->cd(4)->SetLogx();
            cEffPion[iEtaBin]->cd(4)->SetTopMargin(0.035);
            cEffPion[iEtaBin]->cd(4)->SetBottomMargin(0.12);
            cEffPion[iEtaBin]->cd(4)->SetRightMargin(0.035);
            hRatioEffPionTOFDataV0tag[iEtaBin][iEff]->DrawCopy("same");
            cEffPion[iEtaBin]->SaveAs(Form("%s/EfficiencyNsigmaPion_MC_Data_%s.pdf", outDirName.data(), etaBinLabels[iEtaBin].Data()));

            // cEffKaon[iEtaBin]->cd(1)->SetLogx();
            cEffKaon[iEtaBin]->cd(1)->SetTopMargin(0.035);
            cEffKaon[iEtaBin]->cd(1)->SetBottomMargin(0.12);
            cEffKaon[iEtaBin]->cd(1)->SetRightMargin(0.035);
            hEffKaonTPCMCtrue[iEff]->DrawCopy("same");
            hEffKaonTPCDataKinktag[iEtaBin][iEff]->DrawCopy("same");
            hEffKaonTPCDataTOFtag[iEtaBin][iEff]->DrawCopy("same");
            legEffKaonTPC->Draw();

            // cEffKaon[iEtaBin]->cd(2)->SetLogx();
            cEffKaon[iEtaBin]->cd(2)->SetTopMargin(0.035);
            cEffKaon[iEtaBin]->cd(2)->SetBottomMargin(0.12);
            cEffKaon[iEtaBin]->cd(2)->SetRightMargin(0.035);
            hEffKaonTOFMCtrue[iEff]->DrawCopy("same");
            hEffKaonTOFDataTPCtag[iEtaBin][iEff]->DrawCopy("same");
            legEffKaonTOF->Draw();

            // cEffKaon[iEtaBin]->cd(3)->SetLogx();
            cEffKaon[iEtaBin]->cd(3)->SetTopMargin(0.035);
            cEffKaon[iEtaBin]->cd(3)->SetBottomMargin(0.12);
            cEffKaon[iEtaBin]->cd(3)->SetRightMargin(0.035);
            hRatioEffKaonTPCDataKinktag[iEtaBin][iEff]->DrawCopy("same");
            hRatioEffKaonTPCDataTOFtag[iEtaBin][iEff]->DrawCopy("same");

            // cEffKaon[iEtaBin]->cd(4)->SetLogx();
            cEffKaon[iEtaBin]->cd(4)->SetTopMargin(0.035);
            cEffKaon[iEtaBin]->cd(4)->SetBottomMargin(0.12);
            cEffKaon[iEtaBin]->cd(4)->SetRightMargin(0.035);
            hRatioEffKaonTOFDataTPCtag[iEtaBin][iEff]->DrawCopy("same");
            cEffKaon[iEtaBin]->SaveAs(Form("%s/EfficiencyNsigmaKaon_MC_Data_%s.pdf", outDirName.data(), etaBinLabels[iEtaBin].Data()));

            // cEffProton[iEtaBin]->cd(1)->SetLogx();
            cEffProton[iEtaBin]->cd(1)->SetTopMargin(0.035);
            cEffProton[iEtaBin]->cd(1)->SetBottomMargin(0.12);
            cEffProton[iEtaBin]->cd(1)->SetRightMargin(0.035);
            hEffProtonTPCMCtrue[iEff]->DrawCopy("same");
            hEffProtonTPCDataV0tag[iEtaBin][iEff]->DrawCopy("same");
            legEffProton->Draw();

            // cEffProton[iEtaBin]->cd(2)->SetLogx();
            cEffProton[iEtaBin]->cd(2)->SetTopMargin(0.035);
            cEffProton[iEtaBin]->cd(2)->SetBottomMargin(0.12);
            cEffProton[iEtaBin]->cd(2)->SetRightMargin(0.035);
            hEffProtonTOFMCtrue[iEff]->DrawCopy("same");
            hEffProtonTOFDataV0tag[iEtaBin][iEff]->DrawCopy("same");
            legEffProton->Draw();

            // cEffProton[iEtaBin]->cd(3)->SetLogx();
            cEffProton[iEtaBin]->cd(3)->SetTopMargin(0.035);
            cEffProton[iEtaBin]->cd(3)->SetBottomMargin(0.12);
            cEffProton[iEtaBin]->cd(3)->SetRightMargin(0.035);
            hRatioEffProtonTPCDataV0tag[iEtaBin][iEff]->DrawCopy("same");

            // cEffProton[iEtaBin]->cd(4)->SetLogx();
            cEffProton[iEtaBin]->cd(4)->SetTopMargin(0.035);
            cEffProton[iEtaBin]->cd(4)->SetBottomMargin(0.12);
            cEffProton[iEtaBin]->cd(4)->SetRightMargin(0.035);
            hRatioEffProtonTOFDataV0tag[iEtaBin][iEff]->DrawCopy("same");
            cEffProton[iEtaBin]->SaveAs(Form("%s/EfficiencyNsigmaProton_MC_Data_%s.pdf", outDirName.data(), etaBinLabels[iEtaBin].Data()));
        }
    }
    std::cout << "\033[32mDone\033[0m" << std::endl;

    //mean and width of Nsigma distributions
    std::cout << "\n*******************************************\n" << std::endl;
    std::cout << "\033[32mCompute mean and widths of Nsigma TPC distributions for data-driven corrections\033[0m\n" << std::endl;

    TH1D* hMeanPionTPCMCV0tag = new TH1D("hMeanPionTPCMCV0tag", Form(";%s (GeV/#it{c});Mean(#pi)", varTitle.Data()), static_cast<int>(nBins), binLims);
    TH1D* hSigmaPionTPCMCV0tag = new TH1D("hSigmaPionTPCMCV0tag", Form(";%s (GeV/#it{c});Sigma(#pi)", varTitle.Data()), static_cast<int>(nBins), binLims);
    TH1D* hMeanKaonTPCMCTOFtag = new TH1D("hMeanKaonTPCMCTOFtag", Form(";%s (GeV/#it{c});Mean(#pi)", varTitle.Data()), static_cast<int>(nBins), binLims);
    TH1D* hSigmaKaonTPCMCTOFtag = new TH1D("hSigmaKaonTPCMCTOFtag", Form(";%s (GeV/#it{c});Sigma(#pi)", varTitle.Data()), static_cast<int>(nBins), binLims);
    TH1D* hMeanProtonTPCMCV0tag = new TH1D("hMeanProtonTPCMCV0tag", Form(";%s (GeV/#it{c});Mean(#pi)", varTitle.Data()), static_cast<int>(nBins), binLims);
    TH1D* hSigmaProtonTPCMCV0tag = new TH1D("hSigmaProtonTPCMCV0tag", Form(";%s (GeV/#it{c});Sigma(#pi)", varTitle.Data()), static_cast<int>(nBins), binLims);

    SetTH1Style(hMeanPionTPCMCV0tag, markersEffMC[0], pdgColors[kPion]-1, 1., 2, pdgColors[kPion]-1, kWhite, 0.045, 0.055);
    SetTH1Style(hSigmaPionTPCMCV0tag, markersEffMC[0], pdgColors[kPion]-1, 1., 2, pdgColors[kPion]-1, kWhite, 0.045, 0.055);
    SetTH1Style(hMeanKaonTPCMCTOFtag, markersEffMC[0], pdgColors[kKaon]-3, 1., 2, pdgColors[kKaon]-3, kWhite, 0.045, 0.055);
    SetTH1Style(hSigmaKaonTPCMCTOFtag, markersEffMC[0], pdgColors[kKaon]-3, 1., 2, pdgColors[kKaon]-3, kWhite, 0.045, 0.055);
    SetTH1Style(hMeanProtonTPCMCV0tag, markersEffMC[0], pdgColors[kPr], 1., 2, pdgColors[kPr], kWhite, 0.045, 0.055);
    SetTH1Style(hSigmaProtonTPCMCV0tag, markersEffMC[0], pdgColors[kPr], 1., 2, pdgColors[kPr], kWhite, 0.045, 0.055);

    for (unsigned int iBin = 0; iBin < nBins; iBin++)
    {
        hMeanPionTPCMCV0tag->SetBinContent(iBin+1, fNsigmaTPCPionMCV0tag[iBin][kPion]->GetParameter(1));
        hMeanPionTPCMCV0tag->SetBinError(iBin+1, 1.e-20);
        hSigmaPionTPCMCV0tag->SetBinContent(iBin+1, fNsigmaTPCPionMCV0tag[iBin][kPion]->GetParameter(2));
        hSigmaPionTPCMCV0tag->SetBinError(iBin+1, 1.e-20);

        hMeanKaonTPCMCTOFtag->SetBinContent(iBin+1, fNsigmaTPCKaonMCTOFtag[iBin][kKaon]->GetParameter(1));
        hMeanKaonTPCMCTOFtag->SetBinError(iBin+1, 1.e-20);
        hSigmaKaonTPCMCTOFtag->SetBinContent(iBin+1, fNsigmaTPCKaonMCTOFtag[iBin][kKaon]->GetParameter(2));
        hSigmaKaonTPCMCTOFtag->SetBinError(iBin+1, 1.e-20);

        hMeanProtonTPCMCV0tag->SetBinContent(iBin+1, fNsigmaTPCProtonMCV0tag[iBin][kPr]->GetParameter(1));
        hMeanProtonTPCMCV0tag->SetBinError(iBin+1, 1.e-20);
        hSigmaProtonTPCMCV0tag->SetBinContent(iBin+1, fNsigmaTPCProtonMCV0tag[iBin][kPr]->GetParameter(2));
        hSigmaProtonTPCMCV0tag->SetBinError(iBin+1, 1.e-20);
    }

    std::array<TH1D*, nEtaBinsMax+1> hMeanPionTPCDataV0tag, hSigmaPionTPCDataV0tag, hMeanKaonTPCDataTOFtag, hSigmaKaonTPCDataTOFtag, hMeanProtonTPCDataV0tag, hSigmaProtonTPCDataV0tag;

    int nEtaBins4Histos = (nEtaBins == 1) ? nEtaBins : nEtaBins-1;
    TH2D* hMeanPionTPCDataV0tagVsEtaVsP = new TH2D("hMeanPionTPCDataV0tagVsEtaVsP", Form(";%s (GeV/#it{c});|#it{#eta}|;mean", varTitle.Data()), nBins, binLims, nEtaBins4Histos, binEtaLims);
    TH2D* hSigmaPionTPCDataV0tagVsEtaVsP = new TH2D("hSigmaPionTPCDataV0tagVsEtaVsP", Form(";%s (GeV/#it{c});|#it{#eta}|;width", varTitle.Data()), nBins, binLims, nEtaBins4Histos, binEtaLims);
    TH2D* hMeanKaonTPCDataTOFtagVsEtaVsP = new TH2D("hMeanKaonTPCDataTOFtagVsEtaVsP", Form(";%s (GeV/#it{c});|#it{#eta}|;mean", varTitle.Data()), nBins, binLims, nEtaBins4Histos, binEtaLims);
    TH2D* hSigmaKaonTPCDataTOFtagVsEtaVsP = new TH2D("hSigmaKaonTPCDataTOFtagVsEtaVsP", Form(";%s (GeV/#it{c});|#it{#eta}|;width", varTitle.Data()), nBins, binLims, nEtaBins4Histos, binEtaLims);
    TH2D* hMeanProtonTPCDataV0tagVsEtaVsP = new TH2D("hMeanProtonTPCDataV0tagVsEtaVsP", Form(";%s (GeV/#it{c});|#it{#eta}|;mean", varTitle.Data()), nBins, binLims, nEtaBins4Histos, binEtaLims);
    TH2D* hSigmaProtonTPCDataV0tagVsEtaVsP = new TH2D("hSigmaProtonTPCDataV0tagVsEtaVsP", Form(";%s (GeV/#it{c});|#it{#eta}|;width", varTitle.Data()), nBins, binLims, nEtaBins4Histos, binEtaLims);

    TLegend *legPionPars = new TLegend(0.8, 0.7, 0.99, 0.89);
    legPionPars->SetTextSize(0.045);
    TLegend *legKaonPars = new TLegend(0.8, 0.7, 0.99, 0.89);
    legKaonPars->SetTextSize(0.045);
    TLegend *legProtonPars = new TLegend(0.8, 0.7, 0.99, 0.89);
    legProtonPars->SetTextSize(0.045);

    std::array<TCanvas*, nEtaBinsMax+1> cMeanSigma;

    for(unsigned iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++)
    {
        cMeanSigma[iEtaBin] = new TCanvas(Form("cMeanSigma_%s", etaBinLabels[iEtaBin].Data()), Form("cMeanSigma_%s", etaBinLabels[iEtaBin].Data()), 1920, 1080);
        cMeanSigma[iEtaBin]->Divide(3, 2);

        hMeanPionTPCDataV0tag[iEtaBin] = new TH1D(Form("hMeanPionTPCDataV0tag_%s", etaBinLabels[iEtaBin].Data()), Form(";%s (GeV/#it{c});Mean(#pi)", varTitle.Data()), static_cast<int>(nBins), binLims);
        hSigmaPionTPCDataV0tag[iEtaBin] = new TH1D(Form("hSigmaPionTPCDataV0tag_%s", etaBinLabels[iEtaBin].Data()), Form(";%s (GeV/#it{c});Sigma(#pi)", varTitle.Data()), static_cast<int>(nBins), binLims);
        hMeanKaonTPCDataTOFtag[iEtaBin] = new TH1D(Form("hMeanKaonTPCDataTOFtag_%s", etaBinLabels[iEtaBin].Data()), Form(";%s (GeV/#it{c});Mean(#pi)", varTitle.Data()), static_cast<int>(nBins), binLims);
        hSigmaKaonTPCDataTOFtag[iEtaBin] = new TH1D(Form("hSigmaKaonTPCDataTOFtag_%s", etaBinLabels[iEtaBin].Data()), Form(";%s (GeV/#it{c});Sigma(#pi)", varTitle.Data()), static_cast<int>(nBins), binLims);
        hMeanProtonTPCDataV0tag[iEtaBin] = new TH1D(Form("hMeanProtonTPCDataV0tag_%s", etaBinLabels[iEtaBin].Data()), Form(";%s (GeV/#it{c});Mean(#pi)", varTitle.Data()), static_cast<int>(nBins), binLims);
        hSigmaProtonTPCDataV0tag[iEtaBin] = new TH1D(Form("hSigmaProtonTPCDataV0tag_%s", etaBinLabels[iEtaBin].Data()), Form(";%s (GeV/#it{c});Sigma(#pi)", varTitle.Data()), static_cast<int>(nBins), binLims);

        SetTH1Style(hMeanPionTPCDataV0tag[iEtaBin], markersEffData[0], pdgColors[kPion], 1., 2, pdgColors[kPion], kWhite, 0.045, 0.055);
        SetTH1Style(hSigmaPionTPCDataV0tag[iEtaBin], markersEffData[0], pdgColors[kPion], 1., 2, pdgColors[kPion], kWhite, 0.045, 0.055);
        SetTH1Style(hMeanKaonTPCDataTOFtag[iEtaBin], markersEffData[0], pdgColors[kKaon], 1., 2, pdgColors[kKaon], kWhite, 0.045, 0.055);
        SetTH1Style(hSigmaKaonTPCDataTOFtag[iEtaBin], markersEffData[0], pdgColors[kKaon], 1., 2, pdgColors[kKaon], kWhite, 0.045, 0.055);
        SetTH1Style(hMeanProtonTPCDataV0tag[iEtaBin], markersEffData[0], pdgColors[kPr]+1, 1., 2, pdgColors[kPr]+1, kWhite, 0.045, 0.055);
        SetTH1Style(hSigmaProtonTPCDataV0tag[iEtaBin], markersEffData[0], pdgColors[kPr]+1, 1., 2, pdgColors[kPr]+1, kWhite, 0.045, 0.055);

        for (unsigned int iBin = 0; iBin < nBins; iBin++)
        {
            hMeanPionTPCDataV0tag[iEtaBin]->SetBinContent(iBin+1, fNsigmaTPCPionDataV0tag[iEtaBin][iBin][kPion]->GetParameter(1));
            hMeanPionTPCDataV0tag[iEtaBin]->SetBinError(iBin+1, 1.e-20);
            hSigmaPionTPCDataV0tag[iEtaBin]->SetBinContent(iBin+1, fNsigmaTPCPionDataV0tag[iEtaBin][iBin][kPion]->GetParameter(2));
            hSigmaPionTPCDataV0tag[iEtaBin]->SetBinError(iBin+1, 1.e-20);

            hMeanKaonTPCDataTOFtag[iEtaBin]->SetBinContent(iBin+1, fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][kKaon]->GetParameter(1));
            hMeanKaonTPCDataTOFtag[iEtaBin]->SetBinError(iBin+1, 1.e-20);
            hSigmaKaonTPCDataTOFtag[iEtaBin]->SetBinContent(iBin+1, fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][kKaon]->GetParameter(2));
            hSigmaKaonTPCDataTOFtag[iEtaBin]->SetBinError(iBin+1, 1.e-20);

            hMeanProtonTPCDataV0tag[iEtaBin]->SetBinContent(iBin+1, fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][kPr]->GetParameter(1));
            hMeanProtonTPCDataV0tag[iEtaBin]->SetBinError(iBin+1, 1.e-20);
            hSigmaProtonTPCDataV0tag[iEtaBin]->SetBinContent(iBin+1, fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][kPr]->GetParameter(2));
            hSigmaProtonTPCDataV0tag[iEtaBin]->SetBinError(iBin+1, 1.e-20);

            if(nEtaBins == 1 || (iEtaBin < nEtaBins))
            {
                hMeanPionTPCDataV0tagVsEtaVsP->SetBinContent(iBin+1, iEtaBin+1, fNsigmaTPCPionDataV0tag[iEtaBin][iBin][kPion]->GetParameter(1));
                hSigmaPionTPCDataV0tagVsEtaVsP->SetBinContent(iBin+1, iEtaBin+1, fNsigmaTPCPionDataV0tag[iEtaBin][iBin][kPion]->GetParameter(2));
                hMeanKaonTPCDataTOFtagVsEtaVsP->SetBinContent(iBin+1, iEtaBin+1, fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][kKaon]->GetParameter(1));
                hSigmaKaonTPCDataTOFtagVsEtaVsP->SetBinContent(iBin+1, iEtaBin+1, fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][kKaon]->GetParameter(2));
                hMeanProtonTPCDataV0tagVsEtaVsP->SetBinContent(iBin+1, iEtaBin+1, fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][kPr]->GetParameter(1));
                hSigmaProtonTPCDataV0tagVsEtaVsP->SetBinContent(iBin+1, iEtaBin+1, fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][kPr]->GetParameter(2));
            }
        }

        if(iEtaBin == 0)
        {
            legPionPars->AddEntry(hMeanPionTPCMCV0tag, "MC", "p");
            legPionPars->AddEntry(hMeanPionTPCDataV0tag[iEtaBin], "V0 tag", "p");
            legKaonPars->AddEntry(hMeanKaonTPCMCTOFtag, "MC", "p");
            legKaonPars->AddEntry(hMeanKaonTPCDataTOFtag[iEtaBin], "TOF tag", "p");
            legProtonPars->AddEntry(hMeanProtonTPCMCV0tag, "MC", "p");
            legProtonPars->AddEntry(hMeanProtonTPCDataV0tag[iEtaBin], "V0 tag", "p");
        }

        cMeanSigma[iEtaBin]->cd(1)->DrawFrame(binLims[0], -3., binLims[nBins], 3., Form(";%s (GeV/#it{c});mean #it{N}_{#sigma}(#pi)", varTitle.Data()));
        hMeanPionTPCMCV0tag->DrawCopy("same");
        hMeanPionTPCDataV0tag[iEtaBin]->DrawCopy("same");
        legPionPars->Draw();
        cMeanSigma[iEtaBin]->cd(2)->DrawFrame(binLims[0], -3., binLims[nBins], 3., Form(";%s (GeV/#it{c});mean #it{N}_{#sigma}(K)", varTitle.Data()));
        hMeanKaonTPCMCTOFtag->DrawCopy("same");
        hMeanKaonTPCDataTOFtag[iEtaBin]->DrawCopy("same");
        legKaonPars->Draw();
        cMeanSigma[iEtaBin]->cd(3)->DrawFrame(binLims[0], -3., binLims[nBins], 3., Form(";%s (GeV/#it{c});mean #it{N}_{#sigma}(p)", varTitle.Data()));
        hMeanProtonTPCMCV0tag->DrawCopy("same");
        hMeanProtonTPCDataV0tag[iEtaBin]->DrawCopy("same");
        legProtonPars->Draw();
        cMeanSigma[iEtaBin]->cd(4)->DrawFrame(binLims[0], 0., binLims[nBins], 3., Form(";%s (GeV/#it{c});width #it{N}_{#sigma}(#pi)", varTitle.Data()));
        hSigmaPionTPCMCV0tag->DrawCopy("same");
        hSigmaPionTPCDataV0tag[iEtaBin]->DrawCopy("same");
        legPionPars->Draw();
        cMeanSigma[iEtaBin]->cd(5)->DrawFrame(binLims[0], 0., binLims[nBins], 3., Form(";%s (GeV/#it{c});width #it{N}_{#sigma}(K)", varTitle.Data()));
        hSigmaKaonTPCMCTOFtag->DrawCopy("same");
        hSigmaKaonTPCDataTOFtag[iEtaBin]->DrawCopy("same");
        legKaonPars->Draw();
        cMeanSigma[iEtaBin]->cd(6)->DrawFrame(binLims[0], 0., binLims[nBins], 3., Form(";%s (GeV/#it{c});width #it{N}_{#sigma}(p)", varTitle.Data()));
        hSigmaProtonTPCMCV0tag->DrawCopy("same");
        hSigmaProtonTPCDataV0tag[iEtaBin]->DrawCopy("same");
        legProtonPars->Draw();
        cMeanSigma[iEtaBin]->SaveAs(Form("%s/MeanWidthNsigmaTPC_MC_Data_%s.pdf", outDirName.data(), etaBinLabels[iEtaBin].Data()));
    }

    TCanvas* cMeanPionVsEtaVsPt = new TCanvas("cMeanPionVsEtaVsPt", "cMeanPionVsEtaVsPt", 800, 800);
    cMeanPionVsEtaVsPt->SetRightMargin(0.2);
    cMeanPionVsEtaVsPt->SetTopMargin(0.05);
    cMeanPionVsEtaVsPt->SetLogx();
    hMeanPionTPCDataV0tagVsEtaVsP->GetZaxis()->SetRangeUser(-2., 2.);
    hMeanPionTPCDataV0tagVsEtaVsP->GetZaxis()->SetDecimals();
    hMeanPionTPCDataV0tagVsEtaVsP->Draw("colz text");
    cMeanPionVsEtaVsPt->SaveAs(Form("%s/NsigmaTPCPionMeanVsEtaVsPt_Data.pdf", outDirName.data()));
    TCanvas* cSigmaPionVsEtaVsPt = new TCanvas("cSigmaPionVsEtaVsPt", "cSigmaPionVsEtaVsPt", 800, 800);
    cSigmaPionVsEtaVsPt->SetRightMargin(0.2);
    cSigmaPionVsEtaVsPt->SetTopMargin(0.05);
    cSigmaPionVsEtaVsPt->SetLogx();
    hSigmaPionTPCDataV0tagVsEtaVsP->GetZaxis()->SetRangeUser(0., 2.);
    hSigmaPionTPCDataV0tagVsEtaVsP->GetZaxis()->SetDecimals();
    hSigmaPionTPCDataV0tagVsEtaVsP->Draw("colz text");
    cSigmaPionVsEtaVsPt->SaveAs(Form("%s/NsigmaTPCPionWidthVsEtaVsPt_Data.pdf", outDirName.data()));
    TCanvas* cMeanKaonVsEtaVsPt = new TCanvas("cMeanKaonVsEtaVsPt", "cMeanKaonVsEtaVsPt", 800, 800);
    cMeanKaonVsEtaVsPt->SetRightMargin(0.2);
    cMeanKaonVsEtaVsPt->SetTopMargin(0.05);
    cMeanKaonVsEtaVsPt->SetLogx();
    hMeanKaonTPCDataTOFtagVsEtaVsP->GetZaxis()->SetRangeUser(-2., 2.);
    hMeanKaonTPCDataTOFtagVsEtaVsP->GetZaxis()->SetDecimals();
    hMeanKaonTPCDataTOFtagVsEtaVsP->Draw("colz text");
    cMeanKaonVsEtaVsPt->SaveAs(Form("%s/NsigmaTPCKaonWidthVsEtaVsPt_Data.pdf", outDirName.data()));
    TCanvas* cSigmaKaonVsEtaVsPt = new TCanvas("cSigmaKaonVsEtaVsPt", "cSigmaKaonVsEtaVsPt", 800, 800);
    cSigmaKaonVsEtaVsPt->SetRightMargin(0.2);
    cSigmaKaonVsEtaVsPt->SetTopMargin(0.05);
    cSigmaKaonVsEtaVsPt->SetLogx();
    hSigmaKaonTPCDataTOFtagVsEtaVsP->GetZaxis()->SetRangeUser(0., 2.);
    hSigmaKaonTPCDataTOFtagVsEtaVsP->GetZaxis()->SetDecimals();
    hSigmaKaonTPCDataTOFtagVsEtaVsP->Draw("colz text");
    cSigmaKaonVsEtaVsPt->SaveAs(Form("%s/NsigmaTPCKaonWidthVsEtaVsPt_Data.pdf", outDirName.data()));
    TCanvas* cMeanProtonVsEtaVsPt = new TCanvas("cMeanProtonVsEtaVsPt", "cMeanProtonVsEtaVsPt", 800, 800);
    cMeanProtonVsEtaVsPt->SetRightMargin(0.2);
    cMeanProtonVsEtaVsPt->SetTopMargin(0.05);
    cMeanProtonVsEtaVsPt->SetLogx();
    hMeanProtonTPCDataV0tagVsEtaVsP->GetZaxis()->SetRangeUser(-2., 2.);
    hMeanProtonTPCDataV0tagVsEtaVsP->GetZaxis()->SetDecimals();
    hMeanProtonTPCDataV0tagVsEtaVsP->Draw("colz text");
    cMeanProtonVsEtaVsPt->SaveAs(Form("%s/NsigmaTPCProtonMeanVsEtaVsPt_Data.pdf", outDirName.data()));
    TCanvas* cSigmaProtonVsEtaVsPt = new TCanvas("cSigmaProtonVsEtaVsPt", "cSigmaProtonVsEtaVsPt", 800, 800);
    cSigmaProtonVsEtaVsPt->SetRightMargin(0.2);
    cSigmaProtonVsEtaVsPt->SetTopMargin(0.05);
    cSigmaProtonVsEtaVsPt->SetLogx();
    hSigmaProtonTPCDataV0tagVsEtaVsP->GetZaxis()->SetRangeUser(0., 2.);
    hSigmaProtonTPCDataV0tagVsEtaVsP->GetZaxis()->SetDecimals();
    hSigmaProtonTPCDataV0tagVsEtaVsP->Draw("colz text");
    cSigmaProtonVsEtaVsPt->SaveAs(Form("%s/NsigmaTPCProtonWidthVsEtaVsPt_Data.pdf", outDirName.data()));

    std::cout << "\033[32mDone\033[0m" << std::endl;

    //output files
    std::cout << "\n*******************************************\n" << std::endl;
    std::cout << "\033[32mSaving output files\033[0m\n" << std::endl;

    std::array<TDirectoryFile*, nEtaBinsMax+1> dirEtaBinEff;
    TFile outFileEff(Form("%s/PIDEffSystSingleTrack.root", outDirName.data()), "recreate");
    outFileEff.cd();
    for (unsigned int iEff = 0; iEff < nEff; iEff++)
    {
        hEffPionTPCMCtrue[iEff]->Write();
        hEffPionTOFMCtrue[iEff]->Write();
        hEffKaonTPCMCtrue[iEff]->Write();
        hEffKaonTOFMCtrue[iEff]->Write();
        hEffProtonTPCMCtrue[iEff]->Write();
        hEffProtonTOFMCtrue[iEff]->Write();
    }
    for (unsigned int iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++)
    {
        outFileEff.cd();
        dirEtaBinEff[iEtaBin] = new TDirectoryFile(etaBinLabels[iEtaBin], etaBinLabels[iEtaBin]);
        dirEtaBinEff[iEtaBin]->Write();
        dirEtaBinEff[iEtaBin]->cd();
        cEffPion[iEtaBin]->Write();
        cEffKaon[iEtaBin]->Write();
        cEffProton[iEtaBin]->Write();
        for (unsigned int iEff = 0; iEff < nEff; iEff++)
        {
            hEffPionTPCDataV0tag[iEtaBin][iEff]->Write();
            hEffPionTOFDataV0tag[iEtaBin][iEff]->Write();
            hEffKaonTPCDataKinktag[iEtaBin][iEff]->Write();
            hEffKaonTPCDataTOFtag[iEtaBin][iEff]->Write();
            hEffKaonTOFDataTPCtag[iEtaBin][iEff]->Write();
            hEffProtonTPCDataV0tag[iEtaBin][iEff]->Write();
            hEffProtonTOFDataV0tag[iEtaBin][iEff]->Write();
            hRatioEffPionTPCDataV0tag[iEtaBin][iEff]->Write();
            hRatioEffPionTOFDataV0tag[iEtaBin][iEff]->Write();
            hRatioEffKaonTPCDataKinktag[iEtaBin][iEff]->Write();
            hRatioEffKaonTPCDataTOFtag[iEtaBin][iEff]->Write();
            hRatioEffKaonTOFDataTPCtag[iEtaBin][iEff]->Write();
            hRatioEffProtonTPCDataV0tag[iEtaBin][iEff]->Write();
            hRatioEffProtonTOFDataV0tag[iEtaBin][iEff]->Write();
        }
    }
    outFileEff.Close();
    std::cout << Form("File with efficiencies %s/PIDEffSystSingleTrack.root saved", outDirName.data()) << std::endl;    

    std::array<TDirectoryFile*, nEtaBinsMax+1> dirEtaBinDistr;
    TFile outFileDistr(Form("%s/NsigmaPIDdistr.root", outDirName.data()), "recreate");
    outFileDistr.cd();
    TDirectoryFile* dirDataDistr = new TDirectoryFile("Data", "Data");
    dirDataDistr->Write();
    dirDataDistr->cd();
    cNsigmaTPCVsEta->Write();
    cNsigmaTOFVsEta->Write();
    hNsigmaTPCPionDataV0tagVsEta->Write();
    hNsigmaTPCProtonDataV0tagVsEta->Write();
    hNsigmaTPCKaonDataTOFtagVsEta->Write();
    hNsigmaTPCKaonDataKinktagVsEta->Write();
    hNsigmaTOFPionDataV0tagVsEta->Write();
    hNsigmaTOFProtonDataV0tagVsEta->Write();
    hNsigmaTOFKaonDataTPCtagVsEta->Write();
    hNsigmaTOFKaonDataKinktagVsEta->Write();
    hMeanPionTPCDataV0tagVsEtaVsP->Write();
    hSigmaPionTPCDataV0tagVsEtaVsP->Write();
    hMeanKaonTPCDataTOFtagVsEtaVsP->Write();
    hSigmaKaonTPCDataTOFtagVsEtaVsP->Write();
    hMeanProtonTPCDataV0tagVsEtaVsP->Write();
    hSigmaProtonTPCDataV0tagVsEtaVsP->Write();
    for (unsigned int iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++)
    {
        dirDataDistr->cd();
        dirEtaBinDistr[iEtaBin] = new TDirectoryFile(etaBinLabels[iEtaBin], etaBinLabels[iEtaBin]);
        dirEtaBinDistr[iEtaBin]->Write();
        dirEtaBinDistr[iEtaBin]->cd();
        cFitResultTOFKaonFromTPCtag[iEtaBin]->Write();
        cFitResultTOFProtonFromV0tag[iEtaBin]->Write();
        cTOFFractionData[iEtaBin]->Write();
        cPionDataV0tagTPC[iEtaBin]->Write();
        cKaonDataTOFtagTPC[iEtaBin]->Write();
        cKaonDataTOFtagTPC[iEtaBin]->Write();
        cProtonDataV0tagTPC[iEtaBin]->Write();
        cPionDataV0tagTOF[iEtaBin]->Write();
        cKaonDataKinkstagTOF[iEtaBin]->Write();
        cKaonDataTPCtagTOF[iEtaBin]->Write();
        cProtonDataV0tagTOF[iEtaBin]->Write();

        hMeanPionTPCDataV0tag[iEtaBin]->Write();
        hMeanKaonTPCDataTOFtag[iEtaBin]->Write();
        hMeanProtonTPCDataV0tag[iEtaBin]->Write();
        hSigmaPionTPCDataV0tag[iEtaBin]->Write();
        hSigmaKaonTPCDataTOFtag[iEtaBin]->Write();
        hSigmaProtonTPCDataV0tag[iEtaBin]->Write();
        for (unsigned int iBin = 0; iBin < nBins; iBin++)
        {
            hNsigmaTPCPionDataV0tag[iEtaBin][iBin]->Write();
            hNsigmaTPCKaonDataKinktag[iEtaBin][iBin]->Write();
            hNsigmaTPCKaonDataTOFtag[iEtaBin][iBin]->Write();
            hNsigmaTPCProtonDataV0tag[iEtaBin][iBin]->Write();
            hNsigmaTOFPionDataV0tag[iEtaBin][iBin]->Write();
            hNsigmaTOFKaonDataKinktag[iEtaBin][iBin]->Write();
            hNsigmaTOFKaonDataTPCtag[iEtaBin][iBin]->Write();
            hNsigmaTOFProtonDataV0tag[iEtaBin][iBin]->Write();
        }
        for (unsigned int iBin = 0; iBin < nBins; iBin++)
        {
            hNsigmaTOFPionDataV0tagSub[iEtaBin][iBin]->Write();
            hNsigmaTOFKaonDataKinktagSub[iEtaBin][iBin]->Write();
            hNsigmaTOFKaonDataTPCtagSub[iEtaBin][iBin]->Write();
            hNsigmaTOFProtonDataV0tagSub[iEtaBin][iBin]->Write();
        }
        for (unsigned int iBin = 0; iBin < nBins; iBin++)
        {
            for (auto &part : pdgNames)
            {
                if(fNsigmaTPCPionDataV0tag[iEtaBin][iBin][part.first])
                    fNsigmaTPCPionDataV0tag[iEtaBin][iBin][part.first]->Write();
                if(fNsigmaTPCKaonDataKinktag[iEtaBin][iBin][part.first])
                    fNsigmaTPCKaonDataKinktag[iEtaBin][iBin][part.first]->Write();
                if(fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][part.first])
                    fNsigmaTPCKaonDataTOFtag[iEtaBin][iBin][part.first]->Write();
                if(fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][part.first])
                    fNsigmaTPCProtonDataV0tag[iEtaBin][iBin][part.first]->Write();
            }
        }
        for (unsigned int iBin = 0; iBin < nBins; iBin++)
        {
            for (auto &part : pdgNames)
            {
                if (hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][part.first])
                    hNsigmaTOFKaonDataTPCtagFit[iEtaBin][iBin][part.first]->Write();
                if (hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][part.first])
                    hNsigmaTOFProtonDataV0tagFit[iEtaBin][iBin][part.first]->Write();
            }
        }
        for (auto &part : pdgNames)
        {
            if(part.first == kAll)
                continue;
            hFracTOFKaonDataTPCtag[iEtaBin][part.first]->Write();
            hFracTOFProtonDataV0tag[iEtaBin][part.first]->Write();
        }
    }

    outFileDistr.cd();
    TDirectoryFile* dirMCDistr = new TDirectoryFile("MC", "MC");
    dirMCDistr->Write();
    dirMCDistr->cd();
    cPionMCV0tagTPC->Write();
    cKaonMCKinkstagTPC->Write();
    cKaonMCTOFtagTPC->Write();
    cProtonMCV0tagTPC->Write();
    cPionMCV0tagTOF->Write();
    cKaonMCKinkstagTOF->Write();
    cKaonMCTPCtagTOF->Write();
    cProtonMCV0tagTOF->Write();
    hMeanPionTPCMCV0tag->Write();
    hMeanKaonTPCMCTOFtag->Write();
    hMeanProtonTPCMCV0tag->Write();
    hSigmaPionTPCMCV0tag->Write();
    hSigmaKaonTPCMCTOFtag->Write();
    hSigmaProtonTPCMCV0tag->Write();
    cFractionTPCMC->Write();
    cFractionTOFMC->Write();
    for(auto &part : pdgNames)
    {
        if(part.first == kAll)
            continue;
        hFracTPCPionMCV0tag[part.first]->Write();
        hFracTPCKaonMCKinktag[part.first]->Write();
        hFracTPCKaonMCTOFtag[part.first]->Write();
        hFracTPCProtonMCV0tag[part.first]->Write();
        hFracTOFPionMCV0tag[part.first]->Write();
        hFracTOFKaonMCKinktag[part.first]->Write();
        hFracTOFKaonMCTPCtag[part.first]->Write();
        hFracTOFProtonMCV0tag[part.first]->Write();
    }

    for (unsigned int iBin = 0; iBin < nBins; iBin++)
    {
        hNsigmaTPCPionMCTrue[iBin]->Write();
        hNsigmaTPCKaonMCTrue[iBin]->Write();
        hNsigmaTPCProtonMCTrue[iBin]->Write();
        hNsigmaTOFPionMCTrue[iBin]->Write();
        hNsigmaTOFKaonMCTrue[iBin]->Write();
        hNsigmaTOFProtonMCTrue[iBin]->Write();
    }
    for (unsigned int iBin = 0; iBin < nBins; iBin++)
    {
        for (auto &part : pdgNames)
        {
            hNsigmaTPCPionMCV0tag[iBin][part.first]->Write();
            hNsigmaTPCKaonMCKinktag[iBin][part.first]->Write();
            hNsigmaTPCKaonMCTOFtag[iBin][part.first]->Write();
            hNsigmaTPCProtonMCV0tag[iBin][part.first]->Write();
            hNsigmaTOFPionMCV0tag[iBin][part.first]->Write();
            hNsigmaTOFKaonMCKinktag[iBin][part.first]->Write();
            hNsigmaTOFKaonMCTPCtag[iBin][part.first]->Write();
            hNsigmaTOFProtonMCV0tag[iBin][part.first]->Write();
        }
    }
    for (unsigned int iBin = 0; iBin < nBins; iBin++)
    {
        for(auto &part : pdgNames)
        {
            if(part.first == kAll)
                continue;
            fNsigmaTPCPionMCV0tag[iBin][part.first]->Write();
            fNsigmaTPCKaonMCKinktag[iBin][part.first]->Write();
            fNsigmaTPCKaonMCTOFtag[iBin][part.first]->Write();
            fNsigmaTPCProtonMCV0tag[iBin][part.first]->Write();
        }
    }
    outFileDistr.Close();
    
    std::cout << Form("File with distributions %s/NsigmaPIDdistr.root saved\n", outDirName.data()) << std::endl;    
    std::cout << "\033[32mDone\033[0m" << std::endl;
}

//______________________________________________________
void PerformTPCTOFmatchingAnalysis(std::string inFileNameData, std::string inDirNameData, std::string inListNameData,
                                   std::string inFileNameMC, std::string inDirNameMC, std::string inListNameMC,
                                   std::string outDirName, TString varTitle, int var4proj, unsigned int nBins, double binLims[],
                                   unsigned int nEtaBins, std::vector<double> absEtaBinMins, std::vector<double> absEtaBinMaxs, std::vector<TString> etaBinLabels)
{
    //load MC inputs
    auto infileMC = TFile::Open(inFileNameMC.data());
    if (!infileMC || !infileMC->IsOpen())
        return;
    auto indirMC = static_cast<TDirectoryFile *>(infileMC->Get(inDirNameMC.data()));
    if (!indirMC)
    {
        std::cerr << Form("\033[31mERROR: TDirectoryFile %s not found in input file for MC! Exit\033[0m", inDirNameMC.data()) << std::endl;
        return;
    }
    auto listMC = static_cast<TList *>(indirMC->Get(inListNameMC.data()));
    if (!listMC)
    {
        std::cerr << Form("\033[31mERROR: TList %s not found in input file for MC! Exit\033[0m", inListNameMC.data()) << std::endl;
        return;
    }

    std::cout << "\n*******************************************\n" << std::endl;
    std::cout << "\033[32mProject MC tree\033[0m\n" << std::endl;
    ROOT::EnableImplicitMT(); //tell ROOT to go parallel
    ROOT::RDataFrame dataFrameMC(Form("%s/%s", inDirNameMC.data(), "fPIDtree"), inFileNameMC);

    // histos for TPC-TOF marching efficiency vs. p and true particle species    
    std::map<int, TH1D*> hPionMCV0tagWithTOF, hKaonMCTPCtagWithTOF, hProtonMCV0tagWithTOF, hPionMCV0tagAll, hKaonMCTPCtagAll, hProtonMCV0tagAll; 
    TH1D* hTPCTOFMatchEffPionMCV0tag, *hTPCTOFMatchEffKaonMCTPCtag, *hTPCTOFMatchEffProtonMCV0tag, *hTPCTOFMatchEffPionMCtrue, *hTPCTOFMatchEffKaonMCtrue, *hTPCTOFMatchEffProtonMCtrue;
    TH1D* hRatioTPCTOFMatchEffPionMC, *hRatioTPCTOFMatchEffKaonMC, *hRatioTPCTOFMatchEffProtonMC;

    TString pSel = "";
    if(var4proj == kP)
        pSel = "pTPC";
    else
        pSel = "pT";

    double etaLims[101];
    for(int iEtaBin = 0; iEtaBin < 101; iEtaBin++)
        etaLims[iEtaBin] = -1. + 0.02*iEtaBin;

    double hasTOFLims[3] = {-0.5, 0.5, 1.5};

    double partLims[6];
    for(int iPartBin = 0; iPartBin < 6; iPartBin++)
        partLims[iPartBin] = iPartBin;
    std::map<int, double> partBins = {{kEl, 0.5}, {kMuon, 1.5}, {kPion, 2.5}, {kKaon, 3.5}, {kPr, 4.5}};

    auto dataFrameMCEta = dataFrameMC.Filter(Form("(eta > %f && eta < %f) || (eta > -%f && eta < -%f)",
                                                 absEtaBinMins[nEtaBins-1]*1000, absEtaBinMaxs[nEtaBins-1]*1000,
                                                 absEtaBinMaxs[nEtaBins-1]*1000, absEtaBinMins[nEtaBins-1]*1000));

    std::cout << "Selecting V0 tagged pions" << std::endl;
    TString tagSel = Form("(((tag & %d) > 0) || ((tag & %d) > 0))", AliAnalysisTaskSEHFSystPID::kIsPionFromK0s, AliAnalysisTaskSEHFSystPID::kIsPionFromL);
    auto dataFrameMCSel = dataFrameMCEta.Filter(tagSel.Data());
    auto hTOFInfoPionMCV0tagVsPVsPart = dataFrameMCSel.Define("TOF_info", Form("if((trackbits & %d) > 0) return 0; else return 1;", AliAnalysisTaskSEHFSystPID::kHasNoTOF))
                                                      .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                      .Define("part", "if(PDGcode == 11) return 0; else if(PDGcode == 13) return 1; else if(PDGcode == 211) return 2; else if(PDGcode == 321) return 3; else if(PDGcode == 2212) return 4; else return -1;")
                                                      .Histo3D({"hTOFInfoPionMCV0tagVsPVsPart", "", 5u, partLims, static_cast<int>(nBins), binLims, 2u, hasTOFLims}, "part", "p_scaled", "TOF_info");

    std::cout << "Selecting TPC tagged kaons" << std::endl;
    tagSel = Form("((tag & %d) > 0)", AliAnalysisTaskSEHFSystPID::kIsKaonFromTPC);
    dataFrameMCSel = dataFrameMCEta.Filter(tagSel.Data());
    auto hTOFInfoKaonMCTPCtagVsPVsPart = dataFrameMCSel.Define("TOF_info", Form("if((trackbits & %d) > 0) return 0; else return 1;", AliAnalysisTaskSEHFSystPID::kHasNoTOF))
                                                       .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                       .Define("part", "if(PDGcode == 11) return 0; else if(PDGcode == 13) return 1; else if(PDGcode == 211) return 2; else if(PDGcode == 321) return 3; else if(PDGcode == 2212) return 4; else return -1;")
                                                       .Histo3D({"hTOFInfoKaonMCTPCtagVsPVsPart", "", 5u, partLims, static_cast<int>(nBins), binLims, 2u, hasTOFLims}, "part", "p_scaled", "TOF_info");

    std::cout << "Selecting V0 tagged protons" << std::endl;
    tagSel = Form("((tag & %d) > 0)", AliAnalysisTaskSEHFSystPID::kIsProtonFromL);
    dataFrameMCSel = dataFrameMCEta.Filter(tagSel.Data());
    auto hTOFinfoProtonMCV0tagVsPVsPart = dataFrameMCSel.Define("TOF_info", Form("if((trackbits & %d) > 0) return 0; else return 1;", AliAnalysisTaskSEHFSystPID::kHasNoTOF))
                                                        .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                        .Define("part", "if(PDGcode == 11) return 0; else if(PDGcode == 13) return 1; else if(PDGcode == 211) return 2; else if(PDGcode == 321) return 3; else if(PDGcode == 2212) return 4; else return -1;")
                                                        .Histo3D({"hTOFinfoProtonMCV0tagVsPVsPart", "", 5u, partLims, static_cast<int>(nBins), binLims, 2u, hasTOFLims}, "part", "p_scaled", "TOF_info");

    std::cout << "\n\rProcessing pseudorapidity bin 01/01 (\033[32mMC always only eta integrated\033[0m)";
    for(auto &part : pdgNames)
    {
        int partBinMin = hTOFInfoPionMCV0tagVsPVsPart->GetXaxis()->FindBin(partBins[part.first]);
        int partBinMax = hTOFInfoPionMCV0tagVsPVsPart->GetXaxis()->FindBin(partBins[part.first]);
        if(part.first == kAll)
        {
            partBinMin = -1;
            partBinMax = -1;
        }

        hTOFInfoPionMCV0tagVsPVsPart->GetXaxis()->SetRange(partBinMin, partBinMax);
        hTOFInfoPionMCV0tagVsPVsPart->GetZaxis()->SetRange(2, 2);
        hPionMCV0tagWithTOF[part.first] = static_cast<TH1D*>(hTOFInfoPionMCV0tagVsPVsPart->Project3D("y"));
        hPionMCV0tagWithTOF[part.first]->SetNameTitle(Form("hPionMCV0tagWithTOF_%s", part.second.data()),Form(";%s (GeV(#it{c})); Entries", varTitle.Data()));
        hTOFInfoPionMCV0tagVsPVsPart->GetXaxis()->SetRange(-1, -1);
        hTOFInfoPionMCV0tagVsPVsPart->GetZaxis()->SetRange(-1, -1);

        hTOFInfoPionMCV0tagVsPVsPart->GetXaxis()->SetRange(partBinMin, partBinMax);
        hPionMCV0tagAll[part.first] = static_cast<TH1D*>(hTOFInfoPionMCV0tagVsPVsPart->Project3D("y"));
        hPionMCV0tagAll[part.first]->SetNameTitle(Form("hPionMCV0tagAll_%s", part.second.data()), Form(";%s (GeV(#it{c})); Entries", varTitle.Data()));
        hTOFInfoPionMCV0tagVsPVsPart->GetXaxis()->SetRange(-1, -1);

        hTOFInfoKaonMCTPCtagVsPVsPart->GetXaxis()->SetRange(partBinMin, partBinMax);
        hTOFInfoKaonMCTPCtagVsPVsPart->GetZaxis()->SetRange(2, 2);
        hKaonMCTPCtagWithTOF[part.first] = static_cast<TH1D*>(hTOFInfoKaonMCTPCtagVsPVsPart->Project3D("y"));
        hKaonMCTPCtagWithTOF[part.first]->SetNameTitle(Form("hKaonMCTPCtagWithTOF_%s", part.second.data()), Form(";%s (GeV(#it{c})); Entries", varTitle.Data()));
        hTOFInfoKaonMCTPCtagVsPVsPart->GetXaxis()->SetRange(-1, -1);
        hTOFInfoKaonMCTPCtagVsPVsPart->GetZaxis()->SetRange(-1, -1);

        hTOFInfoKaonMCTPCtagVsPVsPart->GetXaxis()->SetRange(partBinMin, partBinMax);
        hKaonMCTPCtagAll[part.first] = static_cast<TH1D*>(hTOFInfoKaonMCTPCtagVsPVsPart->Project3D("y"));
        hKaonMCTPCtagAll[part.first]->SetNameTitle(Form("hKaonMCTPCtagAll_%s", part.second.data()), Form(";%s (GeV(#it{c})); Entries", varTitle.Data()));
        hTOFInfoKaonMCTPCtagVsPVsPart->GetXaxis()->SetRange(-1, -1);

        hTOFinfoProtonMCV0tagVsPVsPart->GetXaxis()->SetRange(partBinMin, partBinMax);
        hTOFinfoProtonMCV0tagVsPVsPart->GetZaxis()->SetRange(2, 2);
        hProtonMCV0tagWithTOF[part.first] = static_cast<TH1D*>(hTOFinfoProtonMCV0tagVsPVsPart->Project3D("y"));
        hProtonMCV0tagWithTOF[part.first]->SetNameTitle(Form("hProtonMCV0tagWithTOF_%s", part.second.data()), Form(";%s (GeV(#it{c})); Entries", varTitle.Data()));
        hTOFinfoProtonMCV0tagVsPVsPart->GetXaxis()->SetRange(-1, -1);
        hTOFinfoProtonMCV0tagVsPVsPart->GetZaxis()->SetRange(-1, -1);

        hTOFinfoProtonMCV0tagVsPVsPart->GetXaxis()->SetRange(partBinMin, partBinMax);
        hProtonMCV0tagAll[part.first] = static_cast<TH1D*>(hTOFinfoProtonMCV0tagVsPVsPart->Project3D("y"));
        hProtonMCV0tagAll[part.first]->SetNameTitle(Form("hProtonMCV0tagAll_%s", part.second.data()), Form(";%s (GeV(#it{c})); Entries", varTitle.Data()));
        hTOFinfoProtonMCV0tagVsPVsPart->GetXaxis()->SetRange(-1, -1);

        if(part.first != kAll)
        {
            SetTH1Style(hPionMCV0tagWithTOF[part.first], kFullCircle, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.055, 0.06);
            SetTH1Style(hPionMCV0tagAll[part.first], kFullCircle, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.055, 0.06);
            SetTH1Style(hKaonMCTPCtagWithTOF[part.first], kFullCircle, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.055, 0.06);
            SetTH1Style(hKaonMCTPCtagAll[part.first], kFullCircle, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.055, 0.06);
            SetTH1Style(hProtonMCV0tagWithTOF[part.first], kFullCircle, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.055, 0.06);
            SetTH1Style(hProtonMCV0tagAll[part.first], kFullCircle, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.055, 0.06);
        }
        else
        {
            SetTH1Style(hPionMCV0tagWithTOF[part.first], kOpenSquare, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.055, 0.06);
            SetTH1Style(hPionMCV0tagAll[part.first], kOpenSquare, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.055, 0.06);
            SetTH1Style(hKaonMCTPCtagWithTOF[part.first], kOpenSquare, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.055, 0.06);
            SetTH1Style(hKaonMCTPCtagAll[part.first], kOpenSquare, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.055, 0.06);
            SetTH1Style(hProtonMCV0tagWithTOF[part.first], kOpenSquare, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.055, 0.06);
            SetTH1Style(hProtonMCV0tagAll[part.first], kOpenSquare, pdgColors[part.first], 1., 2, pdgColors[part.first], kWhite, 0.055, 0.06);
        }
    }

    hTPCTOFMatchEffPionMCV0tag = static_cast<TH1D*>(hPionMCV0tagWithTOF[kAll]->Clone("hTPCTOFMatchEffPionMCV0tag"));
    hTPCTOFMatchEffPionMCV0tag->Divide(hPionMCV0tagWithTOF[kAll], hPionMCV0tagAll[kAll], 1., 1., "B");
    hTPCTOFMatchEffKaonMCTPCtag = static_cast<TH1D*>(hKaonMCTPCtagWithTOF[kAll]->Clone("hTPCTOFMatchEffKaonMCTPCtag"));
    hTPCTOFMatchEffKaonMCTPCtag->Divide(hKaonMCTPCtagWithTOF[kAll], hKaonMCTPCtagAll[kAll], 1., 1., "B");
    hTPCTOFMatchEffProtonMCV0tag = static_cast<TH1D*>(hProtonMCV0tagWithTOF[kAll]->Clone("hTPCTOFMatchEffProtonMCV0tag"));
    hTPCTOFMatchEffProtonMCV0tag->Divide(hProtonMCV0tagWithTOF[kAll], hProtonMCV0tagAll[kAll], 1., 1., "B");

    hTPCTOFMatchEffPionMCtrue = static_cast<TH1D*>(hPionMCV0tagWithTOF[kPion]->Clone("hTPCTOFMatchEffPionMCtrue"));
    hTPCTOFMatchEffPionMCtrue->Divide(hPionMCV0tagWithTOF[kPion], hPionMCV0tagAll[kPion], 1., 1., "B");
    hTPCTOFMatchEffKaonMCtrue = static_cast<TH1D*>(hKaonMCTPCtagWithTOF[kKaon]->Clone("hTPCTOFMatchEffKaonMCtrue"));
    hTPCTOFMatchEffKaonMCtrue->Divide(hKaonMCTPCtagWithTOF[kKaon], hKaonMCTPCtagAll[kKaon], 1., 1., "B");
    hTPCTOFMatchEffProtonMCtrue = static_cast<TH1D*>(hProtonMCV0tagWithTOF[kPr]->Clone("hTPCTOFMatchEffProtonMCtrue"));
    hTPCTOFMatchEffProtonMCtrue->Divide(hProtonMCV0tagWithTOF[kPr], hProtonMCV0tagAll[kPr], 1., 1., "B");

    hRatioTPCTOFMatchEffPionMC = static_cast<TH1D*>(hTPCTOFMatchEffPionMCV0tag->Clone("hRatioTPCTOFMatchEffPionMC"));
    hRatioTPCTOFMatchEffPionMC->Divide(hTPCTOFMatchEffPionMCtrue);
    hRatioTPCTOFMatchEffKaonMC = static_cast<TH1D*>(hTPCTOFMatchEffKaonMCTPCtag->Clone("hRatioTPCTOFMatchEffKaonMC"));
    hRatioTPCTOFMatchEffKaonMC->Divide(hTPCTOFMatchEffKaonMCtrue);
    hRatioTPCTOFMatchEffProtonMC = static_cast<TH1D*>(hTPCTOFMatchEffProtonMCV0tag->Clone("hRatioTPCTOFMatchEffProtonMC"));
    hRatioTPCTOFMatchEffProtonMC->Divide(hTPCTOFMatchEffProtonMCtrue);

    TLegend* letPionMC = new TLegend(0.2, 0.2, 0.5, 0.4);
    letPionMC->SetTextSize(0.045);
    letPionMC->AddEntry(hTPCTOFMatchEffPionMCV0tag, "MC V0 tag", "lp");
    letPionMC->AddEntry(hTPCTOFMatchEffPionMCtrue, "MC true", "lp");

    TLegend* letKaonMC = new TLegend(0.2, 0.2, 0.5, 0.4);
    letKaonMC->SetTextSize(0.045);
    letKaonMC->AddEntry(hTPCTOFMatchEffKaonMCTPCtag, "MC TPC tag", "lp");
    letKaonMC->AddEntry(hTPCTOFMatchEffKaonMCtrue, "MC true", "lp");

    TLegend* letProtonMC = new TLegend(0.2, 0.2, 0.5, 0.4);
    letProtonMC->SetTextSize(0.045);
    letProtonMC->AddEntry(hTPCTOFMatchEffProtonMCV0tag, "MC V0 tag", "lp");
    letProtonMC->AddEntry(hTPCTOFMatchEffProtonMCtrue, "MC true", "lp");

    TCanvas* cTPCTOFMatchEffMCtags = new TCanvas("cTPCTOFMatchEffMCtags", "cTPCTOFMatchEffMCtags", 1920, 1080);
    cTPCTOFMatchEffMCtags->Divide(3, 2);
    cTPCTOFMatchEffMCtags->cd(1)->DrawFrame(binLims[0], 1.e-3, binLims[nBins], 1., Form(";%s (GeV/#it{c});Pion TPC-TOF matching efficiency", varTitle.Data()));
    cTPCTOFMatchEffMCtags->cd(1)->SetLogy();
    hTPCTOFMatchEffPionMCV0tag->DrawCopy("same");
    hTPCTOFMatchEffPionMCtrue->DrawCopy("same");
    letPionMC->Draw();
    cTPCTOFMatchEffMCtags->cd(2)->DrawFrame(binLims[0], 1.e-3, binLims[nBins], 1., Form(";%s (GeV/#it{c});Kaon TPC-TOF matching efficiency", varTitle.Data()));
    cTPCTOFMatchEffMCtags->cd(2)->SetLogy();
    hTPCTOFMatchEffKaonMCTPCtag->DrawCopy("same");
    hTPCTOFMatchEffKaonMCtrue->DrawCopy("same");
    letKaonMC->Draw();
    cTPCTOFMatchEffMCtags->cd(3)->DrawFrame(binLims[0], 1.e-3, binLims[nBins], 1., Form(";%s (GeV/#it{c});Proton TPC-TOF matching efficiency", varTitle.Data()));
    cTPCTOFMatchEffMCtags->cd(3)->SetLogy();
    hTPCTOFMatchEffProtonMCV0tag->DrawCopy("same");
    hTPCTOFMatchEffProtonMCtrue->DrawCopy("same");
    letProtonMC->Draw();
    cTPCTOFMatchEffMCtags->cd(4)->DrawFrame(binLims[0], 0.75, binLims[nBins], 1.25, Form(";%s (GeV/#it{c});Pion match. eff. ratio (MC tag / MC true)", varTitle.Data()));
    hRatioTPCTOFMatchEffPionMC->DrawCopy("same");
    cTPCTOFMatchEffMCtags->cd(5)->DrawFrame(binLims[0], 0.75, binLims[nBins], 1.25, Form(";%s (GeV/#it{c});Kaon match. eff. ratio (MC tag / MC true)", varTitle.Data()));
    hRatioTPCTOFMatchEffKaonMC->DrawCopy("same");
    cTPCTOFMatchEffMCtags->cd(6)->DrawFrame(binLims[0], 0.75, binLims[nBins], 1.25, Form(";%s (GeV/#it{c});Proton match. eff. ratio (MC tag / MC true)", varTitle.Data()));
    hRatioTPCTOFMatchEffProtonMC->DrawCopy("same");

    std::cout << "\n\n\033[32mDone\033[0m" << std::endl;

    //load data inputs
    auto infileData = TFile::Open(inFileNameData.data());
    if (!infileData || !infileData->IsOpen())
        return;
    auto indirData = static_cast<TDirectoryFile *>(infileData->Get(inDirNameData.data()));
    if (!indirData)
    {
        std::cerr << Form("TDirectoryFile %s not found in input file for Data! Exit", inDirNameData.data()) << std::endl;
        return;
    }
    auto listData = static_cast<TList *>(indirData->Get(inListNameData.data()));
    if (!listData)
    {
        std::cerr << Form("TList %s not found in input file for Data! Exit", inListNameData.data()) << std::endl;
        return;
    }

    std::cout << "\n*******************************************\n" << std::endl;
    std::cout << "\033[32mProject data tree\033[0m\n" << std::endl;
    ROOT::RDataFrame dataFrameData(Form("%s/%s", inDirNameData.data(), "fPIDtree"), inFileNameData);

    // histos for TPC-TOF marching efficiency vs. eta and p    
    std::array<TH1D*, nEtaBinsMax+1> hPionDataV0tagWithTOF, hKaonDataTPCtagWithTOF, hProtonDataV0tagWithTOF, hPionDataV0tagAll, hKaonDataTPCtagAll, hProtonDataV0tagAll; 
    std::array<TH1D*, nEtaBinsMax+1> hTPCTOFMatchEffPionDataV0tag, hTPCTOFMatchEffKaonDataTPCtag, hTPCTOFMatchEffProtonDataV0tag;
    std::array<TH1D*, nEtaBinsMax+1> hRatioTPCTOFMatchEffPionDataV0tagMCtrue, hRatioTPCTOFMatchEffKaonDataTPCtagMCtrue, hRatioTPCTOFMatchEffProtonDataV0tagMCtrue;
    std::array<TCanvas*, nEtaBinsMax+1> cTPCTOFMatchEffDatatagsMCtrue;

    TLegend* letPionData = new TLegend(0.2, 0.2, 0.5, 0.4);
    letPionData->SetTextSize(0.045);
    TLegend* letKaonData = new TLegend(0.2, 0.2, 0.5, 0.4);
    letKaonData->SetTextSize(0.045);
    TLegend* letProtonData = new TLegend(0.2, 0.2, 0.5, 0.4);
    letProtonData->SetTextSize(0.045);

    std::cout << "Selecting V0 tagged pions" << std::endl;
    tagSel = Form("(((tag & %d) > 0) || ((tag & %d) > 0))", AliAnalysisTaskSEHFSystPID::kIsPionFromK0s, AliAnalysisTaskSEHFSystPID::kIsPionFromL);
    auto dataFrameDataSel = dataFrameData.Filter(tagSel.Data());
    auto hTOFInfoPionDataV0tagVsEtaVsP = dataFrameDataSel.Define("TOF_info", Form("if((trackbits & %d) > 0) return 0; else return 1;", AliAnalysisTaskSEHFSystPID::kHasNoTOF))
                                                         .Define("eta_scaled", "static_cast<float>(eta)/1000")
                                                         .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                         .Histo3D({"hTOFInfoPionDataV0tagVsEtaVsP", "", static_cast<int>(nBins), binLims, 100u, etaLims, 2u, hasTOFLims}, "p_scaled", "eta_scaled", "TOF_info");

    std::cout << "Selecting TPC tagged kaons" << std::endl;
    tagSel = Form("((tag & %d) > 0)", AliAnalysisTaskSEHFSystPID::kIsKaonFromTPC);
    dataFrameDataSel = dataFrameData.Filter(tagSel.Data());
    auto hTOFInfoKaonDataTPCtagVsEtaVsP = dataFrameDataSel.Define("TOF_info", Form("if((trackbits & %d) > 0) return 0; else return 1;", AliAnalysisTaskSEHFSystPID::kHasNoTOF))
                                                          .Define("eta_scaled", "static_cast<float>(eta)/1000")
                                                          .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                          .Histo3D({"hTOFInfoKaonDataTPCtagVsEtaVsP", "", static_cast<int>(nBins), binLims, 100u, etaLims, 2u, hasTOFLims}, "p_scaled", "eta_scaled", "TOF_info");

    std::cout << "Selecting V0 tagged protons" << std::endl;
    tagSel = Form("((tag & %d) > 0)", AliAnalysisTaskSEHFSystPID::kIsProtonFromL);
    dataFrameDataSel = dataFrameData.Filter(tagSel.Data());
    auto hTOFInfoProtonDataV0tagVsEtaVsP = dataFrameDataSel.Define("TOF_info", Form("if((trackbits & %d) > 0) return 0; else return 1;", AliAnalysisTaskSEHFSystPID::kHasNoTOF))
                                                           .Define("eta_scaled", "static_cast<float>(eta)/1000")
                                                           .Define("p_scaled", Form("static_cast<float>(%s)/1000", pSel.Data()))
                                                           .Histo3D({"hTOFInfoProtonDataV0tagVsEtaVsP", "", static_cast<int>(nBins), binLims, 100u, etaLims, 2u, hasTOFLims}, "p_scaled", "eta_scaled", "TOF_info");

    std::cout << std::endl;
    for(unsigned int iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++)
    {
        std::cout << Form("\rProcessing pseudorapidity bin %02d/%02d", iEtaBin+1, nEtaBins);
        int etaMinBinPos = hTOFInfoPionDataV0tagVsEtaVsP->GetYaxis()->FindBin(absEtaBinMins[iEtaBin]*1.0001);
        int etaMaxBinPos = hTOFInfoPionDataV0tagVsEtaVsP->GetYaxis()->FindBin(absEtaBinMaxs[iEtaBin]*0.9999);
        int etaMinBinNeg = hTOFInfoPionDataV0tagVsEtaVsP->GetYaxis()->FindBin(-absEtaBinMaxs[iEtaBin]*0.9999);
        int etaMaxBinNeg = hTOFInfoPionDataV0tagVsEtaVsP->GetYaxis()->FindBin(-absEtaBinMins[iEtaBin]*1.0001);

        hTOFInfoPionDataV0tagVsEtaVsP->GetYaxis()->SetRange(etaMinBinPos, etaMaxBinPos);
        hTOFInfoPionDataV0tagVsEtaVsP->GetZaxis()->SetRange(2, 2);
        hPionDataV0tagWithTOF[iEtaBin] = static_cast<TH1D*>(hTOFInfoPionDataV0tagVsEtaVsP->Project3D("x"));
        hTOFInfoPionDataV0tagVsEtaVsP->GetYaxis()->SetRange(-1, -1);
        hTOFInfoPionDataV0tagVsEtaVsP->GetYaxis()->SetRange(etaMinBinNeg, etaMaxBinNeg);
        hPionDataV0tagWithTOF[iEtaBin]->Add(static_cast<TH1D*>(hTOFInfoPionDataV0tagVsEtaVsP->Project3D("x")));
        hPionDataV0tagWithTOF[iEtaBin]->SetNameTitle(Form("hPionDataV0tagWithTOF_%s", etaBinLabels[iEtaBin].Data()), Form(";%s (GeV(#it{c})); Entries", varTitle.Data()));
        hTOFInfoPionDataV0tagVsEtaVsP->GetYaxis()->SetRange(-1, -1);
        hTOFInfoPionDataV0tagVsEtaVsP->GetZaxis()->SetRange(-1, -1);

        hTOFInfoPionDataV0tagVsEtaVsP->GetYaxis()->SetRange(etaMinBinPos, etaMaxBinPos);
        hPionDataV0tagAll[iEtaBin] = static_cast<TH1D*>(hTOFInfoPionDataV0tagVsEtaVsP->Project3D("x"));
        hTOFInfoPionDataV0tagVsEtaVsP->GetYaxis()->SetRange(-1, -1);
        hTOFInfoPionDataV0tagVsEtaVsP->GetYaxis()->SetRange(etaMinBinNeg, etaMaxBinNeg);
        hPionDataV0tagAll[iEtaBin]->Add(static_cast<TH1D*>(hTOFInfoPionDataV0tagVsEtaVsP->Project3D("x")));
        hPionDataV0tagAll[iEtaBin]->SetNameTitle(Form("hPionDataV0tagAll_%s", etaBinLabels[iEtaBin].Data()), Form(";%s (GeV(#it{c})); Entries", varTitle.Data()));
        hTOFInfoPionDataV0tagVsEtaVsP->GetYaxis()->SetRange(-1, -1);

        hTOFInfoKaonDataTPCtagVsEtaVsP->GetYaxis()->SetRange(etaMinBinPos, etaMaxBinPos);
        hTOFInfoKaonDataTPCtagVsEtaVsP->GetZaxis()->SetRange(2, 2);
        hKaonDataTPCtagWithTOF[iEtaBin] = static_cast<TH1D*>(hTOFInfoKaonDataTPCtagVsEtaVsP->Project3D("x"));
        hTOFInfoKaonDataTPCtagVsEtaVsP->GetYaxis()->SetRange(-1, -1);
        hTOFInfoKaonDataTPCtagVsEtaVsP->GetYaxis()->SetRange(etaMinBinNeg, etaMaxBinNeg);
        hKaonDataTPCtagWithTOF[iEtaBin]->Add(static_cast<TH1D*>(hTOFInfoKaonDataTPCtagVsEtaVsP->Project3D("x")));
        hKaonDataTPCtagWithTOF[iEtaBin]->SetNameTitle(Form("hKaonDataTPCtagWithTOF_%s", etaBinLabels[iEtaBin].Data()), Form(";%s (GeV(#it{c})); Entries", varTitle.Data()));
        hTOFInfoKaonDataTPCtagVsEtaVsP->GetYaxis()->SetRange(-1, -1);
        hTOFInfoKaonDataTPCtagVsEtaVsP->GetZaxis()->SetRange(-1, -1);

        hTOFInfoKaonDataTPCtagVsEtaVsP->GetYaxis()->SetRange(etaMinBinPos, etaMaxBinPos);
        hKaonDataTPCtagAll[iEtaBin] = static_cast<TH1D*>(hTOFInfoKaonDataTPCtagVsEtaVsP->Project3D("x"));
        hTOFInfoKaonDataTPCtagVsEtaVsP->GetYaxis()->SetRange(-1, -1);
        hTOFInfoKaonDataTPCtagVsEtaVsP->GetYaxis()->SetRange(etaMinBinNeg, etaMaxBinNeg);
        hKaonDataTPCtagAll[iEtaBin]->Add(static_cast<TH1D*>(hTOFInfoKaonDataTPCtagVsEtaVsP->Project3D("x")));
        hKaonDataTPCtagAll[iEtaBin]->SetNameTitle(Form("hKaonDataTPCtagAll_%s", etaBinLabels[iEtaBin].Data()), Form(";%s (GeV(#it{c})); Entries", varTitle.Data()));
        hTOFInfoKaonDataTPCtagVsEtaVsP->GetYaxis()->SetRange(-1, -1);

        hTOFInfoProtonDataV0tagVsEtaVsP->GetYaxis()->SetRange(etaMinBinPos, etaMaxBinPos);
        hTOFInfoProtonDataV0tagVsEtaVsP->GetZaxis()->SetRange(2, 2);
        hProtonDataV0tagWithTOF[iEtaBin] = static_cast<TH1D*>(hTOFInfoProtonDataV0tagVsEtaVsP->Project3D("x"));
        hTOFInfoProtonDataV0tagVsEtaVsP->GetYaxis()->SetRange(-1, -1);
        hTOFInfoProtonDataV0tagVsEtaVsP->GetYaxis()->SetRange(etaMinBinNeg, etaMaxBinNeg);
        hProtonDataV0tagWithTOF[iEtaBin]->Add(static_cast<TH1D*>(hTOFInfoProtonDataV0tagVsEtaVsP->Project3D("x")));
        hProtonDataV0tagWithTOF[iEtaBin]->SetNameTitle(Form("hProtonDataV0tagWithTOF_%s", etaBinLabels[iEtaBin].Data()), Form(";%s (GeV(#it{c})); Entries", varTitle.Data()));
        hTOFInfoProtonDataV0tagVsEtaVsP->GetYaxis()->SetRange(-1, -1);
        hTOFInfoProtonDataV0tagVsEtaVsP->GetZaxis()->SetRange(-1, -1);

        hTOFInfoProtonDataV0tagVsEtaVsP->GetYaxis()->SetRange(etaMinBinPos, etaMaxBinPos);
        hProtonDataV0tagAll[iEtaBin] = static_cast<TH1D*>(hTOFInfoProtonDataV0tagVsEtaVsP->Project3D("x"));
        hTOFInfoProtonDataV0tagVsEtaVsP->GetYaxis()->SetRange(-1, -1);
        hTOFInfoProtonDataV0tagVsEtaVsP->GetYaxis()->SetRange(etaMinBinNeg, etaMaxBinNeg);
        hProtonDataV0tagAll[iEtaBin]->Add(static_cast<TH1D*>(hTOFInfoProtonDataV0tagVsEtaVsP->Project3D("x")));
        hProtonDataV0tagAll[iEtaBin]->SetNameTitle(Form("hProtonDataV0tagAll_%s", etaBinLabels[iEtaBin].Data()), Form(";%s (GeV(#it{c})); Entries", varTitle.Data()));
        hTOFInfoProtonDataV0tagVsEtaVsP->GetYaxis()->SetRange(-1, -1);

        SetTH1Style(hPionDataV0tagWithTOF[iEtaBin], kFullDiamond, pdgColors[kPion]+1, 1., 2, pdgColors[kPion], kWhite, 0.055, 0.06);
        SetTH1Style(hPionDataV0tagAll[iEtaBin], kOpenDiamond, pdgColors[kPion]+1, 1., 2, pdgColors[kPion], kWhite, 0.055, 0.06);
        SetTH1Style(hKaonDataTPCtagWithTOF[iEtaBin], kFullDiamond, pdgColors[kKaon]-1, 1., 2, pdgColors[kKaon], kWhite, 0.055, 0.06);
        SetTH1Style(hKaonDataTPCtagAll[iEtaBin], kOpenDiamond, pdgColors[kKaon]-1, 1., 2, pdgColors[kKaon], kWhite, 0.055, 0.06);
        SetTH1Style(hProtonDataV0tagWithTOF[iEtaBin], kFullDiamond, pdgColors[kPr]+1, 1., 2, pdgColors[kPr], kWhite, 0.055, 0.06);
        SetTH1Style(hProtonDataV0tagAll[iEtaBin], kOpenDiamond, pdgColors[kPr]+1, 1., 2, pdgColors[kPr], kWhite, 0.055, 0.06);

        hTPCTOFMatchEffPionDataV0tag[iEtaBin] = static_cast<TH1D*>(hPionDataV0tagWithTOF[iEtaBin]->Clone("hTPCTOFMatchEffPionDataV0tag"));
        hTPCTOFMatchEffPionDataV0tag[iEtaBin]->Divide(hPionDataV0tagWithTOF[iEtaBin], hPionDataV0tagAll[iEtaBin], 1., 1., "B");
        hTPCTOFMatchEffKaonDataTPCtag[iEtaBin] = static_cast<TH1D*>(hKaonDataTPCtagWithTOF[iEtaBin]->Clone("hTPCTOFMatchEffKaonDataTPCtag"));
        hTPCTOFMatchEffKaonDataTPCtag[iEtaBin]->Divide(hKaonDataTPCtagWithTOF[iEtaBin], hKaonDataTPCtagAll[iEtaBin], 1., 1., "B");
        hTPCTOFMatchEffProtonDataV0tag[iEtaBin] = static_cast<TH1D*>(hProtonDataV0tagWithTOF[iEtaBin]->Clone("hTPCTOFMatchEffProtonDataV0tag"));
        hTPCTOFMatchEffProtonDataV0tag[iEtaBin]->Divide(hProtonDataV0tagWithTOF[iEtaBin], hProtonDataV0tagAll[iEtaBin], 1., 1., "B");

        hRatioTPCTOFMatchEffPionDataV0tagMCtrue[iEtaBin] = static_cast<TH1D*>(hTPCTOFMatchEffPionDataV0tag[iEtaBin]->Clone("hRatioTPCTOFMatchEffPionDataV0tagMCtrue"));
        hRatioTPCTOFMatchEffPionDataV0tagMCtrue[iEtaBin]->Divide(hTPCTOFMatchEffPionMCtrue);
        hRatioTPCTOFMatchEffKaonDataTPCtagMCtrue[iEtaBin] = static_cast<TH1D*>(hTPCTOFMatchEffKaonDataTPCtag[iEtaBin]->Clone("hRatioTPCTOFMatchEffKaonDataTPCtagMCtrue"));
        hRatioTPCTOFMatchEffKaonDataTPCtagMCtrue[iEtaBin]->Divide(hTPCTOFMatchEffKaonMCtrue);
        hRatioTPCTOFMatchEffProtonDataV0tagMCtrue[iEtaBin] = static_cast<TH1D*>(hTPCTOFMatchEffProtonDataV0tag[iEtaBin]->Clone("hRatioTPCTOFMatchEffProtonDataV0tagMCtrue"));
        hRatioTPCTOFMatchEffProtonDataV0tagMCtrue[iEtaBin]->Divide(hTPCTOFMatchEffProtonMCtrue);

        if(iEtaBin == 0)
        {
            letPionData->AddEntry(hTPCTOFMatchEffPionDataV0tag[iEtaBin],"Data V0 tag","lp");
            letPionData->AddEntry(hTPCTOFMatchEffPionMCtrue,"MC true","lp");
            letKaonData->AddEntry(hTPCTOFMatchEffKaonDataTPCtag[iEtaBin],"Data TPC tag","lp");
            letKaonData->AddEntry(hTPCTOFMatchEffKaonMCtrue,"MC true","lp");
            letProtonData->AddEntry(hTPCTOFMatchEffProtonDataV0tag[iEtaBin],"Data V0 tag","lp");
            letProtonData->AddEntry(hTPCTOFMatchEffProtonMCtrue,"MC true","lp");
        }
        cTPCTOFMatchEffDatatagsMCtrue[iEtaBin] = new TCanvas(Form("cTPCTOFMatchEffDatatagsMCtrue_%s", etaBinLabels[iEtaBin].Data()), Form("cTPCTOFMatchEffDatatagsMCtrue_%s", etaBinLabels[iEtaBin].Data()), 1920, 1080);
        cTPCTOFMatchEffDatatagsMCtrue[iEtaBin]->Divide(3, 2);
        cTPCTOFMatchEffDatatagsMCtrue[iEtaBin]->cd(1)->DrawFrame(binLims[0], 1.e-3, binLims[nBins], 1., Form(";%s (GeV/#it{c});Pion TPC-TOF matching efficiency", varTitle.Data()));
        cTPCTOFMatchEffDatatagsMCtrue[iEtaBin]->cd(1)->SetLogy();
        hTPCTOFMatchEffPionMCtrue->DrawCopy("same");
        hTPCTOFMatchEffPionDataV0tag[iEtaBin]->DrawCopy("same");
        letPionData->Draw();
        cTPCTOFMatchEffDatatagsMCtrue[iEtaBin]->cd(2)->DrawFrame(binLims[0], 1.e-3, binLims[nBins], 1., Form(";%s (GeV/#it{c});Kaon TPC-TOF matching efficiency", varTitle.Data()));
        cTPCTOFMatchEffDatatagsMCtrue[iEtaBin]->cd(2)->SetLogy();
        hTPCTOFMatchEffKaonMCtrue->DrawCopy("same");
        hTPCTOFMatchEffKaonDataTPCtag[iEtaBin]->DrawCopy("same");
        letKaonData->Draw();
        cTPCTOFMatchEffDatatagsMCtrue[iEtaBin]->cd(3)->DrawFrame(binLims[0], 1.e-3, binLims[nBins], 1., Form(";%s (GeV/#it{c});Proton TPC-TOF matching efficiency", varTitle.Data()));
        cTPCTOFMatchEffDatatagsMCtrue[iEtaBin]->cd(3)->SetLogy();
        hTPCTOFMatchEffProtonMCtrue->DrawCopy("same");
        hTPCTOFMatchEffProtonDataV0tag[iEtaBin]->DrawCopy("same");
        letProtonData->Draw();
        cTPCTOFMatchEffDatatagsMCtrue[iEtaBin]->cd(4)->DrawFrame(binLims[0], 0.75, binLims[nBins], 1.25, Form(";%s (GeV/#it{c});Pion match. eff. ratio (Data tag / MC true)", varTitle.Data()));
        hRatioTPCTOFMatchEffPionDataV0tagMCtrue[iEtaBin]->DrawCopy("same");
        cTPCTOFMatchEffDatatagsMCtrue[iEtaBin]->cd(5)->DrawFrame(binLims[0], 0.75, binLims[nBins], 1.25, Form(";%s (GeV/#it{c});Kaon match. eff. ratio (Data tag / MC true)", varTitle.Data()));
        hRatioTPCTOFMatchEffKaonDataTPCtagMCtrue[iEtaBin]->DrawCopy("same");
        cTPCTOFMatchEffDatatagsMCtrue[iEtaBin]->cd(6)->DrawFrame(binLims[0], 0.75, binLims[nBins], 1.25, Form(";%s (GeV/#it{c});Proton match. eff. ratio (Data tag / MC true)", varTitle.Data()));
        hRatioTPCTOFMatchEffProtonDataV0tagMCtrue[iEtaBin]->DrawCopy("same");

        cTPCTOFMatchEffDatatagsMCtrue[iEtaBin]->SaveAs(Form("%s/TOFTPCMatchingEfficiency_data_MC_%s.pdf", outDirName.data(), etaBinLabels[iEtaBin].Data()));
    }

    std::cout << "\n\n\033[32mDone\033[0m" << std::endl;

    //output files
    std::cout << "\n*******************************************\n" << std::endl;
    std::cout << "\033[32mSaving output files\033[0m\n" << std::endl;

    std::array<TDirectoryFile*, nEtaBinsMax+1> dirEtaBinEff;
    TFile outFileMatchEff(Form("%s/TPCTOFMatchingEffSystSingleTrack.root", outDirName.data()), "recreate");
    outFileMatchEff.cd();
    cTPCTOFMatchEffMCtags->Write();
    hTPCTOFMatchEffPionMCV0tag->Write();
    hTPCTOFMatchEffPionMCtrue->Write();
    hTPCTOFMatchEffKaonMCTPCtag->Write();
    hTPCTOFMatchEffKaonMCtrue->Write();
    hTPCTOFMatchEffProtonMCV0tag->Write();
    hTPCTOFMatchEffProtonMCtrue->Write();
    hRatioTPCTOFMatchEffPionMC->Write();
    hRatioTPCTOFMatchEffKaonMC->Write();
    hRatioTPCTOFMatchEffProtonMC->Write();
    for (unsigned int iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++)
    {
        outFileMatchEff.cd();
        dirEtaBinEff[iEtaBin] = new TDirectoryFile(etaBinLabels[iEtaBin], etaBinLabels[iEtaBin]);
        dirEtaBinEff[iEtaBin]->Write();
        dirEtaBinEff[iEtaBin]->cd();
        cTPCTOFMatchEffDatatagsMCtrue[iEtaBin]->Write();
        hTPCTOFMatchEffPionDataV0tag[iEtaBin]->Write();
        hTPCTOFMatchEffKaonDataTPCtag[iEtaBin]->Write();
        hTPCTOFMatchEffProtonDataV0tag[iEtaBin]->Write();
        hRatioTPCTOFMatchEffPionDataV0tagMCtrue[iEtaBin]->Write();
        hRatioTPCTOFMatchEffKaonDataTPCtagMCtrue[iEtaBin]->Write();
        hRatioTPCTOFMatchEffProtonDataV0tagMCtrue[iEtaBin]->Write();
    }
    outFileMatchEff.Close();

    std::cout << Form("File with TPC-TOF matching efficiencies %s/TPCTOFMatchingEffSystSingleTrack.root saved", outDirName.data()) << std::endl;
}

//______________________________________________________
void ComputeEfficiency(double num, double den, double &eff, double &effunc)
{

    TH1D htmpnum("htmpnum", "", 1, 0, 1);
    TH1D htmpden("htmpden", "", 1, 0, 1);
    TH1D htmpratio("htmpratio", "", 1, 0, 1);

    htmpnum.SetBinContent(1, num);
    htmpden.SetBinContent(1, den);
    htmpnum.SetBinError(1, TMath::Sqrt(num));
    htmpden.SetBinError(1, TMath::Sqrt(den));
    htmpratio.Divide(&htmpnum, &htmpden, 1., 1., "B");
    eff = htmpratio.GetBinContent(1);
    effunc = htmpratio.GetBinError(1);
}

//______________________________________________________
void GetTOFFractionsFromData(int whichpart, unsigned int iBin, std::map<int, TH1D*> hFracMC, std::map<int, TH1D*> hFracData, std::map<int, TH1D*> hNsigmaMC, TH1D *hNsigmaData, TFractionFitter *&fNsigmaFitter, std::vector<int> &templUsed)
{
    TObjArray *oNsigmaMC = new TObjArray(0);

    for (auto &part : pdgNames)
    {
        if(part.first == kAll)
            continue;

        if (part.first == whichpart || hFracMC[part.first]->GetBinContent(iBin+1) > 1.e-5)
        {
            TH1D *hMC = dynamic_cast<TH1D *>(hNsigmaMC[part.first]->Clone(Form("hMC_%d_%s", iBin, part.second.data())));
            oNsigmaMC->AddLast(hMC);
            templUsed.push_back(part.first);
        }
    }

    if (templUsed.size() == 1)
    {
        for(auto &part : pdgNames)
        {
            if (part.first == whichpart)
                hFracData[part.first]->SetBinContent(iBin+1, 1.);
            else
                hFracData[part.first]->SetBinContent(iBin+1, 0.);
        }
        return;
    }

    if (oNsigmaMC->GetEntries() > 1)
    {

        fNsigmaFitter = new TFractionFitter(hNsigmaData, oNsigmaMC, "Q");
        for (int iEntry = 0; iEntry < oNsigmaMC->GetEntries(); iEntry++)
            fNsigmaFitter->Constrain(iEntry, hFracMC[templUsed[iEntry]]->GetBinContent(iBin+1) * 0.5, hFracMC[templUsed[iEntry]]->GetBinContent(iBin+1) * 2);
        fNsigmaFitter->Fit();

        for(auto &part : pdgNames)
        {
            if(part.first == kAll)
                continue;

            std::vector<int>::iterator it = find(templUsed.begin(), templUsed.end(), part.first);
            if (it != templUsed.end())
            {
                double frac, err;
                int iEntry = static_cast<int>(distance(templUsed.begin(), it));

                fNsigmaFitter->GetResult(iEntry, frac, err);
                hFracData[part.first]->SetBinContent(iBin+1, frac);
                hFracData[part.first]->SetBinError(iBin+1, err);
            }
            else
            {
                hFracData[part.first]->SetBinContent(iBin+1, 0);
                hFracData[part.first]->SetBinError(iBin+1, 0);
            }
        }
    }
    else
    {
        for(auto &part : pdgNames)
        {
            if(part.first == kAll)
                continue;

            if (part.first != whichpart)
            {
                hFracData[part.first]->SetBinContent(iBin+1, 1);
                hFracData[part.first]->SetBinError(iBin+1, 1);
            }
            else
            {
                hFracData[part.first]->SetBinContent(iBin+1, 0);
                hFracData[part.first]->SetBinError(iBin+1, 0);
            }
        }
    }
}

//______________________________________________________
double PDFnsigmaTPCtot(double *nsigma, double *pars)
{

    double partpdf[5];
    double totalpdf = 0;
    for (auto &part : pdgPosition)
    {
        partpdf[part.second] = pars[3 * part.second] * TMath::Gaus(nsigma[0], pars[3 * part.second + 1], pars[3 * part.second + 2]);
        totalpdf += partpdf[part.second];
    }

    return totalpdf;
}

//_____________________________________________________
void PlotQAhistos(TList *listMC, TList *listData, string outDirName)
{

    //V0 QA histos
    TString ArmenterosName[5] = {"All", "K0s", "Lambda", "AntiLambda", "Gamma"};
    TString legV0name[4] = {"K_{s}^{0} #rightarrow #pi^{+}#pi^{-}", "#Lambda #rightarrow p#pi^{-}", "#bar{#Lambda} #rightarrow #bar{p}#pi^{+}", "#gamma #rightarrow e^{+}e^{-}"};
    TExec *exGrey = new TExec("exGrey","gStyle->SetPalette(kGreyScale);");
    TExec *exOrange = new TExec("exBlue","gStyle->SetPalette(kSolar);");
    TExec *exRed = new TExec("exBlue","gStyle->SetPalette(kCherry);");
    TExec *exBlue = new TExec("exBlue","gStyle->SetPalette(kDeepSea);");
    TExec *exGreen = new TExec("exBlue","gStyle->SetPalette(kAvocado);");
    TExec *exStd = new TExec("exStd","gStyle->SetPalette(kRainBow);");
    std::vector<TExec*> ArmenterosPalette = {exGrey, exOrange, exRed, exBlue, exGreen};
    int ArmenterosColor[5] = {kBlack, kOrange+7, kRed+1, kAzure+4, kGreen+2};
    TH2D *hArmenterosMC[5];
    TH2D *hArmenterosData[5];

    //kinks QA histos
    TString KinksName[5] = {"QtVsMassKinks", "PDaughterVsMotherKink", "OpeningAngleVsPMotherKink", "dEdxVsPMotherKink", "NTPCclsVsRadius"};
    TH2D *hKinksMC[5];
    TH2D *hKinksData[5];

    TCanvas *cArmentero = new TCanvas("cArmenteros", "cArmenteros", 1920, 1080);
    TCanvas *cKinks = new TCanvas("cKinks", "cKinks", 1920, 1080);
    cArmentero->Divide(2, 1);
    cKinks->Divide(5, 2);
    TLatex *latV0[5];

    for (int iHisto = 0; iHisto < 5; iHisto++)
    {
        hArmenterosMC[iHisto] = (TH2D *)listMC->FindObject(Form("fHistArmenteroPlot%s", ArmenterosName[iHisto].Data()));
        hArmenterosData[iHisto] = (TH2D *)listData->FindObject(Form("fHistArmenteroPlot%s", ArmenterosName[iHisto].Data()));
        hArmenterosMC[iHisto]->SetDirectory(0);
        hArmenterosData[iHisto]->SetDirectory(0);
        hArmenterosData[iHisto]->GetYaxis()->SetDecimals();
        hArmenterosMC[iHisto]->GetYaxis()->SetDecimals();
        hArmenterosData[iHisto]->SetTitle("Data");
        hArmenterosMC[iHisto]->SetTitle("MC");
        latV0[iHisto] = new TLatex();
        latV0[iHisto]->SetNDC();
        latV0[iHisto]->SetTextFont(42);
        latV0[iHisto]->SetTextColor(ArmenterosColor[iHisto]);
        latV0[iHisto]->SetTextSize(0.045);
        if (iHisto != 0)
        {
            cArmentero->cd(1)->SetLogz();
            if(iHisto == 1)
            {
                hArmenterosMC[iHisto]->DrawCopy("axis");
                ArmenterosPalette[iHisto]->Draw();
                hArmenterosMC[iHisto]->DrawCopy("col same");
            }
            else
            {
                ArmenterosPalette[iHisto]->Draw();
                hArmenterosMC[iHisto]->DrawCopy("col same");
            }
            latV0[iHisto]->DrawLatex(0.6, 0.8 - 0.05 * iHisto, legV0name[iHisto - 1].Data());
            cArmentero->cd(2)->SetLogz();
            if(iHisto == 1)
            {
                hArmenterosData[iHisto]->DrawCopy("axis");
                ArmenterosPalette[iHisto]->Draw();
                hArmenterosData[iHisto]->DrawCopy("col same");
            }
            else
            {
                ArmenterosPalette[iHisto]->Draw();
                hArmenterosData[iHisto]->DrawCopy("col same");
            }
            latV0[iHisto]->DrawLatex(0.6, 0.8 - 0.05 * iHisto, legV0name[iHisto - 1].Data());
        }

        hKinksMC[iHisto] = (TH2D *)listMC->FindObject(Form("fHist%s", KinksName[iHisto].Data()));
        hKinksData[iHisto] = (TH2D *)listData->FindObject(Form("fHist%s", KinksName[iHisto].Data()));
        hKinksMC[iHisto]->SetTitle("MC");
        hKinksData[iHisto]->SetTitle("Data");

        cKinks->cd(1 + iHisto)->SetLogz();
        hKinksMC[iHisto]->DrawCopy("axis");
        exStd->Draw();
        hKinksMC[iHisto]->DrawCopy("colz same");
        cKinks->cd(6 + iHisto)->SetLogz();
        hKinksData[iHisto]->DrawCopy("axis");
        exStd->Draw();
        hKinksData[iHisto]->DrawCopy("colz same");
    }

    TH1D *hKinksMassMC = static_cast<TH1D*>(hKinksMC[0]->ProjectionX("hKinksMassMC"));
    hKinksMassMC->Scale(1./hKinksMassMC->Integral());
    SetTH1Style(hKinksMassMC, kFullSquare, kAzure+4, 1., 2, kAzure+4, kWhite, 0.04, 0.045);
    hKinksMassMC->SetTitle("");
    hKinksMassMC->GetYaxis()->SetTitle("Normalised entries");

    TH1D *hKinksMassData = static_cast<TH1D*>(hKinksData[0]->ProjectionX("hKinksMassData"));
    hKinksMassData->Scale(1./hKinksMassData->Integral());
    SetTH1Style(hKinksMassData, kFullDiamond, kRed+1, 1., 2, kRed+1, kWhite, 0.04, 0.045);
    hKinksMassData->SetTitle("");
    hKinksMassData->GetYaxis()->SetTitle("Normalised entries");

    TLegend *leg = new TLegend(0.2, 0.7, 0.4, 0.8);
    leg->SetTextSize(0.045);
    leg->AddEntry(hKinksMassMC, "MC", "lpe");
    leg->AddEntry(hKinksMassData, "Data", "lpe");
    TCanvas *cKinksMass = new TCanvas("cKinksMass", "", 800, 800);
    cKinksMass->SetTopMargin(0.1);
    hKinksMassMC->DrawCopy();
    hKinksMassData->DrawCopy("same");
    leg->Draw();

    cKinks->SaveAs(Form("%s/KinksQAplots.pdf", outDirName.data()));
    cKinksMass->SaveAs(Form("%s/KinksMass.pdf", outDirName.data()));
    cArmentero->SaveAs(Form("%s/V0QAplots.pdf", outDirName.data()));
}

//_____________________________________________________
void DivideCanvas(TCanvas *c, int nBins)
{
    if (nBins < 2)
        c->cd();
    else if (nBins == 2 || nBins == 3)
        c->Divide(nBins, 1);
    else if (nBins == 4 || nBins == 6 || nBins == 8)
        c->Divide(nBins / 2, 2);
    else if (nBins == 5 || nBins == 7)
        c->Divide((nBins + 1) / 2, 2);
    else if (nBins == 9 || nBins == 12 || nBins == 15)
        c->Divide(nBins / 3, 3);
    else if (nBins == 10 || nBins == 11)
        c->Divide(4, 3);
    else if (nBins == 13 || nBins == 14)
        c->Divide(5, 3);
    else if (nBins > 15 && nBins <= 20 && nBins % 4 == 0)
        c->Divide(nBins / 4, 4);
    else if (nBins > 15 && nBins <= 20 && nBins % 4 != 0)
        c->Divide(5, 4);
    else if (nBins == 21)
        c->Divide(7, 3);
    else if (nBins > 21 && nBins <= 25)
        c->Divide(5, 5);
    else if (nBins > 25 && nBins % 2 == 0)
        c->Divide(nBins / 2, 2);
    else
        c->Divide((nBins + 1) / 2, 2);
}

//_____________________________________________________
void SetStyle()
{
    gStyle->SetLegendBorderSize(1);
    gStyle->SetTitleOffset(1.5, "y");
    gStyle->SetTitleOffset(1.2, "x");
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetTitleSize(0.05, "xyzt");
    gStyle->SetLabelSize(0.045, "xyz");
    gStyle->SetPaintTextFormat(".2f");

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetOptStat(0);

    gROOT->ForceStyle();

    TGaxis::SetMaxDigits(3);
}

//_____________________________________________________
void SetTH1Style(TH1 *histo, int markerstyle, int markercolor, float markersize, int linewidth, int linecolor, int fillcolor, float labelsize, float titlesize)
{

    histo->SetMarkerStyle(markerstyle);
    histo->SetMarkerSize(markersize);
    histo->SetMarkerColor(markercolor);
    histo->SetLineWidth(linewidth);
    histo->SetLineColor(linecolor);
    histo->SetFillColorAlpha(fillcolor, 0.25);

    if (labelsize > 0)
    {
        histo->GetXaxis()->SetLabelSize(labelsize);
        histo->GetYaxis()->SetLabelSize(labelsize);
    }
    if (titlesize > 0)
    {
        histo->GetXaxis()->SetTitleSize(labelsize);
        histo->GetYaxis()->SetTitleSize(labelsize);
        histo->SetTitleSize(labelsize);
    }
}
