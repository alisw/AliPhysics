#if !defined (__CINT__) || defined (__CLING__)

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <regex>

#include "yaml-cpp/yaml.h"

#include "TROOT.h"
#include "Riostream.h"
#include "TH1F.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TArrayD.h"

#include "AliHFSystErr.h"

#endif

//________________________________________________________________________
// Macro to produce pp reference files from the outputs of https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/macros/HFPtSpectrum.C
// Main function: GenppRefFile
// Parameters for the extrapolation are passed to the macro via a yaml config file such as https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/macros/config_ppref.yml
// Possible FONLL extrapolations to be configured:
//    - rapidity
//    - transverse momentum
//
// TODO: add energy extrapolation
//
// Author: F. Grosa, fabrizio.grosa@cern.ch
//________________________________________________________________________

//________________________________________________________________________
// FUNCTION PROTOTYPES
void GenppRefFile(std::string cfgFileName="config_ppref.yml");
void CreateExtrapReferenceHistos(const int nPtBinsMeas, const double* ptLimsMeas, std::vector<double> ptLimsExtr, std::string FONLLFileName, TH1F* hSigmaPP, TGraphAsymmErrors* gSigmaPPSyst, TGraphAsymmErrors* gSigmaPPSystData, TGraphAsymmErrors* gSigmaPPSystFeedDown, TH1F* &hReference, TH1F* &hReferenceSystData, TH1F* &hCombinedReferenceFlag, TGraphAsymmErrors* &gReference, TGraphAsymmErrors* &gReferenceSyst, TGraphAsymmErrors* &gReferenceFdSyst, std::vector<std::string> scalesToExclude, bool useFDUnc=true);
void ApplyYshiftToHisto(TGraphAsymmErrors *gyShift, TH1F* histo);
void ApplyYshiftToGraph(TGraphAsymmErrors *gyShift, TGraphAsymmErrors* graph, bool useUnc=false);
TGraphAsymmErrors* ComputeCorrForYshift(std::string FONLLfileNameNum, std::string FONLLfileNameDen, const int nPtBinsMeas, const double* ptLimsMeas);
std::map<std::string, TH1F*> ReadFONLL(std::string FONLLFileName, std::string suffix="");
void SetStyle();

//________________________________________________________________________
// FUNCTION IMPLEMENTATIONS
void GenppRefFile(std::string cfgFileName) {

    // set drawing style
    SetStyle();

    // load config
    YAML::Node config = YAML::LoadFile(cfgFileName.data());
    std::string mesonName = config["meson"].as<std::string>();
    std::string inFileName = config["infilename"].as<std::string>();
    std::vector<double> ptLimsRebVec = config["ptbins"].as<std::vector<double> >();
    bool doPtExtrap = static_cast<bool>(config["extrap"]["pt"]["doextrap"].as<int>());
    std::string FONLLFileName;
    std::vector<double> ptLimsExtrap;
    std::vector<std::string> scalesToExclude;
    if(doPtExtrap)
    {
        FONLLFileName = config["extrap"]["pt"]["FONLLfilename"].as<std::string>();
        ptLimsExtrap = config["extrap"]["pt"]["ptlims"].as<std::vector<double> >();
        scalesToExclude = config["extrap"]["pt"]["excludescales"].as<std::vector<std::string> >();
    }
    bool doYExtrap = static_cast<bool>(config["extrap"]["y"]["doextrap"].as<int>());
    std::string FONLLFileNameMidY;
    std::string FONLLFileNameShiftY;
    if(doYExtrap)
    {
        FONLLFileNameMidY = config["extrap"]["y"]["FONLLfilename"]["midy"].as<std::string>();
        FONLLFileNameShiftY = config["extrap"]["y"]["FONLLfilename"]["shifty"].as<std::string>();
    }

    // load measured cross section
    TFile* infile = TFile::Open(inFileName.data());
    if(!infile)
        return;

    TH1F* hSigmaPPOrig = static_cast<TH1F*>(infile->Get("histoSigmaCorr"));
    if(!hSigmaPPOrig)
    {
        std::cerr << "ERROR: histo with statistical uncertainties not found! Exit" << std::endl;
        return;
    }
    hSigmaPPOrig->SetDirectory(0);

    TGraphAsymmErrors* gSigmaPPSystDataOrig = nullptr;
    if(mesonName.compare("Lc") != 0)
        gSigmaPPSystDataOrig = static_cast<TGraphAsymmErrors*>(infile->Get("gSigmaCorr"));
    else
        gSigmaPPSystDataOrig = static_cast<TGraphAsymmErrors*>(infile->Get("gSigmaAverageSyst"));
    if(!gSigmaPPSystDataOrig)
    {
        std::cerr << "ERROR: graph with systematic uncertainties not found! Exit" << std::endl;
        return;
    }

    TGraphAsymmErrors* gSigmaPPSystFeedDownOrig = nullptr;
    if(mesonName.compare("DzeroLowPtInclusive") != 0)
    {
        if(mesonName.compare("Lc") != 0)
            gSigmaPPSystFeedDownOrig = static_cast<TGraphAsymmErrors*>(infile->Get("gSigmaCorrConservative"));
        else
        {
            TGraphAsymmErrors* gfPrompt = static_cast<TGraphAsymmErrors*>(infile->Get("gFDLcAverage"));
            gSigmaPPSystFeedDownOrig = static_cast<TGraphAsymmErrors*>(gSigmaPPSystDataOrig->Clone("gSigmaCorrConservative"));
            for(int iPt=0; iPt<gfPrompt->GetN(); iPt++)
            {
                double pt=-1., sigma=-1., fprompt=-1.;
                gfPrompt->GetPoint(iPt, pt, fprompt);
                gSigmaPPSystFeedDownOrig->GetPoint(iPt, pt, sigma);
                double systLow = gfPrompt->GetErrorYlow(iPt) / fprompt * sigma;
                double systHigh = gfPrompt->GetErrorYhigh(iPt) / fprompt * sigma;
                gSigmaPPSystFeedDownOrig->SetPointError(iPt, 0.3, 0.3, systLow, systHigh);
            }
        }

        if(!gSigmaPPSystFeedDownOrig)
        {
            std::cerr << "ERROR: graph with FD systematic uncertainties not found! Exit" << std::endl;
            return;
        }
    }

    AliHFSystErr* syst = static_cast<AliHFSystErr*>(infile->Get("AliHFSystErr"));
    if(!syst)
        syst = static_cast<AliHFSystErr*>(infile->Get("AliHFSystErrTopol"));
    if(!syst)
    {
        std::cerr << "ERROR: AliHFSystErr object not found! Exit" << std::endl;
        return;
    }

    TH1F* hSystTrackingAverage = nullptr;
    TH1F* hSystPtShapeAverage = nullptr;
    TH1F* hSystPidAverage = nullptr;
    TH1F* hSystBRAverage = nullptr;
    TH1F* hSystFDHighAverage = nullptr;
    TH1F* hSystFDLowAverage = nullptr;
    TH1F* hSystRawYieldAverage = nullptr;
    TH1F* hSystCutVarAverage = nullptr;
    if(mesonName.compare("Lc") == 0)
    {
        hSystTrackingAverage = static_cast<TH1F*>(infile->Get("hSystTrackingAverage"));
        hSystPtShapeAverage = static_cast<TH1F*>(infile->Get("hSystPtShapeAverage"));
        hSystPidAverage = static_cast<TH1F*>(infile->Get("hSystPidAverage"));
        hSystBRAverage = static_cast<TH1F*>(infile->Get("hSystBRAverage"));
        hSystFDHighAverage = static_cast<TH1F*>(infile->Get("hSystFDHighAverage"));
        hSystFDLowAverage = static_cast<TH1F*>(infile->Get("hSystFDLowAverage"));
        hSystRawYieldAverage = static_cast<TH1F*>(infile->Get("hSystRawYieldAverage"));
        hSystCutVarAverage = static_cast<TH1F*>(infile->Get("hSystCutVarAverage"));
        hSystTrackingAverage->SetDirectory(0);
        hSystPtShapeAverage->SetDirectory(0);
        hSystPidAverage->SetDirectory(0);
        hSystBRAverage->SetDirectory(0);
        hSystFDHighAverage->SetDirectory(0);
        hSystFDLowAverage->SetDirectory(0);
        hSystRawYieldAverage->SetDirectory(0);
        hSystCutVarAverage->SetDirectory(0);
    }
    infile->Close();

    // rebin if needed
    const int nPtBinsReb = ptLimsRebVec.size()-1;
    double ptLimsReb[nPtBinsReb+1];
    TH1F* hSigmaPPReb = nullptr;
    TGraphAsymmErrors *gSigmaPPSystDataReb = nullptr, *gSigmaPPSystFeedDownReb = nullptr;
    if(nPtBinsReb > 0)
    {
        for(int iPt=0; iPt<=nPtBinsReb; iPt++)
            ptLimsReb[iPt] = ptLimsRebVec[iPt];
        std::cout << hSigmaPPOrig->GetBinError(hSigmaPPOrig->GetNbinsX()) << std::endl;
        for(int iPtOrig=0; iPtOrig<hSigmaPPOrig->GetNbinsX(); iPtOrig++)
        {
            hSigmaPPOrig->SetBinContent(iPtOrig+1, hSigmaPPOrig->GetBinContent(iPtOrig+1)*hSigmaPPOrig->GetBinWidth(iPtOrig+1));
            hSigmaPPOrig->SetBinError(iPtOrig+1, hSigmaPPOrig->GetBinError(iPtOrig+1)*hSigmaPPOrig->GetBinWidth(iPtOrig+1));
        }

        hSigmaPPReb = static_cast<TH1F*>(hSigmaPPOrig->Rebin(nPtBinsReb, Form("%s_reb", hSigmaPPOrig->GetName()), ptLimsReb));
        hSigmaPPReb->Scale(1., "width");
        hSigmaPPOrig->Scale(1., "width");
        std::cout << hSigmaPPOrig->GetBinError(hSigmaPPOrig->GetNbinsX()) << std::endl;
        gSigmaPPSystDataReb = new TGraphAsymmErrors(0);
        gSigmaPPSystFeedDownReb = new TGraphAsymmErrors(0);
        for(int iPtReb=0; iPtReb<nPtBinsReb; iPtReb++)
        {
            double ptCent = hSigmaPPReb->GetBinCenter(iPtReb+1);
            double sigma = hSigmaPPReb->GetBinContent(iPtReb+1);
            double binWidth = hSigmaPPReb->GetBinWidth(iPtReb+1);
            double systUncorr = 0.;
            double systCorr = 0.;
            double systFDLow = 0.;
            double systFDHigh = 0.;

            double ptMin = hSigmaPPReb->GetXaxis()->GetBinLowEdge(iPtReb+1);
            double ptMax = hSigmaPPReb->GetXaxis()->GetBinUpEdge(iPtReb+1);
            int ptBinOrigLow = hSigmaPPOrig->FindBin(ptMin*1.0001);
            int ptBinOrigHigh = hSigmaPPOrig->FindBin(ptMax*0.9999);

            if(mesonName.compare("Lc") != 0)
            {                    
                for(int iPtOrig=ptBinOrigLow; iPtOrig<=ptBinOrigHigh; iPtOrig++)
                {
                    double ptOrig = hSigmaPPOrig->GetBinCenter(iPtOrig);
                    double sigmaOrig  = hSigmaPPOrig->GetBinContent(iPtOrig) * hSigmaPPOrig->GetBinWidth(iPtOrig);
                    double sysRawY    = syst->GetRawYieldErr(ptOrig);
                    double sysCutV    = syst->GetSeleEffErr(ptOrig);
                    double sysPID     = syst->GetPIDEffErr(ptOrig);
                    double sysPtShape = syst->GetMCPtShapeErr(ptOrig);
                    double sysTrack   = syst->GetTrackingEffErr(ptOrig);

                    systUncorr += sysRawY * sysRawY * sigmaOrig * sigmaOrig;
                    systCorr += TMath::Sqrt(sysCutV*sysCutV + sysPID*sysPID + sysPtShape*sysPtShape + sysTrack*sysTrack) * sigmaOrig;
                    for(iPtGraph=0; iPtGraph<gSigmaPPSystFeedDownOrig->GetN(); iPtGraph++)
                    {
                        double ptGraph, y;
                        gSigmaPPSystFeedDownOrig->GetPoint(iPtGraph, ptGraph, y);
                        if(TMath::Abs(ptGraph - ptOrig) < 1.e-3)
                        {
                            systFDLow += gSigmaPPSystFeedDownOrig->GetErrorYlow(iPtGraph) * sigmaOrig;
                            systFDHigh += gSigmaPPSystFeedDownOrig->GetErrorYhigh(iPtGraph) * sigmaOrig;
                            break;
                        }
                    }
                }
            }
            else
            {
                for(int iPtOrig=ptBinOrigLow; iPtOrig<=ptBinOrigHigh; iPtOrig++)
                {
                    double sigmaOrig  = hSigmaPPOrig->GetBinContent(iPtOrig) * hSigmaPPOrig->GetBinWidth(iPtOrig);
                    double sysRawY    = hSystRawYieldAverage->GetBinContent(iPtOrig);
                    double sysCutV    = hSystCutVarAverage->GetBinContent(iPtOrig);
                    double sysBR      = hSystBRAverage->GetBinContent(iPtOrig);
                    double sysPID     = hSystPidAverage->GetBinContent(iPtOrig);
                    double sysPtShape = hSystPtShapeAverage->GetBinContent(iPtOrig);
                    double sysTrack   = hSystTrackingAverage->GetBinContent(iPtOrig);

                    systUncorr += sysRawY * sysRawY * sigmaOrig * sigmaOrig;
                    systCorr += TMath::Sqrt(sysCutV*sysCutV + sysBR*sysBR + sysPID*sysPID + sysPtShape*sysPtShape + sysTrack*sysTrack) * sigmaOrig;
                    systFDLow += hSystFDLowAverage->GetBinContent(iPtOrig) * sigmaOrig;
                    systFDHigh += hSystFDHighAverage->GetBinContent(iPtOrig) * sigmaOrig;
                }
            }
            systUncorr = TMath::Sqrt(systUncorr) / binWidth;
            systCorr /= binWidth;
            double systUncData = TMath::Sqrt(systUncorr*systUncorr + systCorr*systCorr);
            systFDLow /= binWidth;
            systFDHigh /= binWidth;

            gSigmaPPSystDataReb->SetPoint(iPtReb, ptCent, sigma);
            gSigmaPPSystFeedDownReb->SetPoint(iPtReb, ptCent, sigma);
            gSigmaPPSystDataReb->SetPointError(iPtReb, 0.3, 0.3, systUncData, systUncData);
            gSigmaPPSystFeedDownReb->SetPointError(iPtReb, 0.3, 0.3, systFDLow, systFDHigh);
        }
    }

    const int nPtBinsMeas = hSigmaPPOrig->GetNbinsX();
    const TArrayD* ptLimsMeasArray = hSigmaPPOrig->GetXaxis()->GetXbins();
    const double* ptLimsMeas = ptLimsMeasArray->GetArray();

    TH1F* hSigmaPP = nullptr;
    TGraphAsymmErrors *gSigmaPPSystData = nullptr, *gSigmaPPSystFeedDown = nullptr;

    if(hSigmaPPReb)
    {
        hSigmaPP = static_cast<TH1F*>(hSigmaPPReb->Clone("fhScaledData"));
        gSigmaPPSystData = static_cast<TGraphAsymmErrors*>(gSigmaPPSystDataReb->Clone("gScaledDataSystData"));
        if(mesonName.compare("kDzeroLowPtInclusive") != 0)
            gSigmaPPSystFeedDown = static_cast<TGraphAsymmErrors*>(gSigmaPPSystFeedDownReb->Clone("gScaledDataSystFeedDown"));
    }
    else
    {
        hSigmaPP = static_cast<TH1F*>(hSigmaPPOrig->Clone("fhScaledData"));
        gSigmaPPSystData = static_cast<TGraphAsymmErrors*>(gSigmaPPSystDataOrig->Clone("gScaledDataSystData"));
        if(mesonName.compare("kDzeroLowPtInclusive") != 0)
            gSigmaPPSystFeedDown = static_cast<TGraphAsymmErrors*>(gSigmaPPSystFeedDownOrig->Clone("gScaledDataSystFeedDown"));
    }

    if(doPtExtrap)
    {
        if(ptLimsExtrap[0] != hSigmaPP->GetXaxis()->GetBinUpEdge(hSigmaPP->GetNbinsX()))
        {
            std::cerr << ptLimsExtrap[0] << "  " << hSigmaPP->GetXaxis()->GetBinUpEdge(hSigmaPP->GetNbinsX()) << "ERROR: extrapolated and measured pT bins do not match! Please check" << std::endl;
            return;
        }
    }

    TGraphAsymmErrors* gSigmaPPSystTheory = static_cast<TGraphAsymmErrors*>(gSigmaPPSystData->Clone("gScaledDataSystExtrap"));
    TGraphAsymmErrors* gSigmaPPSyst = static_cast<TGraphAsymmErrors*>(gSigmaPPSystData->Clone("gScaledData"));

    for(int iPt=0; iPt<hSigmaPP->GetNbinsX(); iPt++) {
        gSigmaPPSystTheory->SetPointEYlow(iPt, 0.);
        gSigmaPPSystTheory->SetPointEYhigh(iPt, 0.);

        double systdatahigh = gSigmaPPSystData->GetErrorYhigh(iPt);
        double systdatalow = gSigmaPPSystData->GetErrorYlow(iPt);
        double systFDhigh = 0.;
        double systFDlow = 0.;
        if(mesonName.compare("kDzeroLowPtInclusive") != 0)
        {
            systFDhigh = gSigmaPPSystFeedDown->GetErrorYhigh(iPt);
            systFDlow = gSigmaPPSystFeedDown->GetErrorYlow(iPt);
        }

        double systtotlow = TMath::Sqrt(systdatalow*systdatalow+systFDlow*systFDlow);
        double systtothigh = TMath::Sqrt(systdatahigh*systdatahigh+systFDhigh*systFDhigh);
        gSigmaPPSyst->SetPointEYlow(iPt,systtotlow);
        gSigmaPPSyst->SetPointEYhigh(iPt,systtothigh);
    }

    if(doYExtrap)
    {
        TGraphAsymmErrors* gyShift = nullptr;
        if(hSigmaPPReb)
            gyShift = ComputeCorrForYshift(FONLLFileNameMidY, FONLLFileNameShiftY, nPtBinsReb, ptLimsReb);
        else
            gyShift = ComputeCorrForYshift(FONLLFileNameMidY, FONLLFileNameShiftY, nPtBinsMeas, ptLimsMeas);

        ApplyYshiftToHisto(gyShift, hSigmaPP);
        ApplyYshiftToGraph(gyShift, gSigmaPPSyst, true);
        ApplyYshiftToGraph(gyShift, gSigmaPPSystData, false);
        if(mesonName.compare("DzeroLowPtInclusive") != 0)
            ApplyYshiftToGraph(gyShift, gSigmaPPSystFeedDown, false);
        ApplyYshiftToGraph(gyShift, gSigmaPPSystTheory, true);
    }

    TH1F* hReference = nullptr;
    TH1F* hReferenceSystData = nullptr;
    TH1F* hCombinedReferenceFlag = nullptr;
    TGraphAsymmErrors* gReference = nullptr;
    TGraphAsymmErrors* gReferenceSyst = nullptr;
    TGraphAsymmErrors* gReferenceFdSyst = nullptr;
    if(doPtExtrap)
    {
        if(hSigmaPPReb)
            CreateExtrapReferenceHistos(nPtBinsReb, ptLimsReb, ptLimsExtrap, FONLLFileName, hSigmaPP, gSigmaPPSyst, gSigmaPPSystData, gSigmaPPSystFeedDown, hReference, hReferenceSystData,   hCombinedReferenceFlag, gReference, gReferenceSyst, gReferenceFdSyst, scalesToExclude, static_cast<bool>(mesonName.compare("DzeroLowPtInclusive") != 0));
        else
            CreateExtrapReferenceHistos(nPtBinsMeas, ptLimsMeas, ptLimsExtrap, FONLLFileName, hSigmaPP, gSigmaPPSyst, gSigmaPPSystData, gSigmaPPSystFeedDown, hReference, hReferenceSystData,   hCombinedReferenceFlag, gReference, gReferenceSyst, gReferenceFdSyst, scalesToExclude, static_cast<bool>(mesonName.compare("DzeroLowPtInclusive") != 0));
    }

    TH1F* hSigmaPPSyst = static_cast<TH1F*>(hSigmaPP->Clone("fhScaledSystData"));
    for(int iPt=0; iPt<hSigmaPP->GetNbinsX(); iPt++)
    {
        double ptcent = hSigmaPPSyst->GetBinCenter(iPt+1);
        double systunc = -1;
        for(int iPtG=0; iPtG<gSigmaPPSystData->GetN(); iPtG++)
        {
            double ptgraph = -1., ygraph = -1.;
            gSigmaPPSystData->GetPoint(iPtG, ptgraph, ygraph);
            if(TMath::Abs(ptgraph-ptcent)<0.001) {
                systunc = gSigmaPPSystData->GetErrorYlow(iPtG);
                break;
            }
        }
        hSigmaPPSyst->SetBinError(iPt+1, systunc);
    }

    // output file name
    std::string outFileName = mesonName + std::string("_ppreference_pp5TeV");
    if(doYExtrap)
        outFileName.append("_yshift");
    else
        outFileName.append("_noyshift");

    outFileName.append("_pt");
    for(int iPt=0; iPt<nPtBinsMeas+1; iPt++)
    {
        outFileName.append("_");
        if(hSigmaPPReb)
        {
            if(TMath::Abs(static_cast<int>(ptLimsReb[iPt]) - ptLimsReb[iPt]) > 0.1)
                outFileName.append(Form("%0.1f", ptLimsReb[iPt]));
            else
                outFileName.append(Form("%0.f", ptLimsReb[iPt]));
            if(iPt > nPtBinsReb-1)
                break;
        }
        else
        {
            if(TMath::Abs(static_cast<int>(ptLimsMeas[iPt]) - ptLimsMeas[iPt]) > 0.1)
                outFileName.append(Form("%0.1f", ptLimsMeas[iPt]));
            else
                outFileName.append(Form("%0.f", ptLimsMeas[iPt]));
        }
    }

    if(doPtExtrap)
    {
        outFileName.append("_FONLLextrap_pt");
        for(unsigned int iPt=0; iPt<ptLimsExtrap.size(); iPt++)
        {
            outFileName.append("_");
        if(TMath::Abs(static_cast<int>(ptLimsExtrap[iPt]) - ptLimsExtrap[iPt]) > 0.1)
            outFileName.append(Form("%0.1f", ptLimsExtrap[iPt]));
        else
            outFileName.append(Form("%0.f", ptLimsExtrap[iPt]));
        }
    }
    outFileName = std::regex_replace(outFileName, std::regex("\\."), "dot");
    outFileName.append(".root");

    TFile outFile(outFileName.data(), "recreate");
    hSigmaPP->Write();
    hSigmaPPSyst->Write();
    gSigmaPPSyst->Write();
    gSigmaPPSystData->Write();
    if(mesonName.compare("kDzeroLowPtInclusive") != 0)
        gSigmaPPSystFeedDown->Write();
    gSigmaPPSystTheory->Write();
    if(doPtExtrap)
    {
        hReference->Write();
        hReferenceSystData->Write();
        hCombinedReferenceFlag->Write();
        gReference->Write();
        gReferenceSyst->Write();
        gReferenceFdSyst->Write();
    }
    if(syst)
        syst->Write("AliHFSystErr");
    outFile.Close();
}

//________________________________________________________________________
void CreateExtrapReferenceHistos(const int nPtBinsMeas, const double* ptLimsMeas, std::vector<double> ptLimsExtr, std::string FONLLFileName, TH1F* hSigmaPP, TGraphAsymmErrors* gSigmaPPSyst, TGraphAsymmErrors* gSigmaPPSystData, TGraphAsymmErrors* gSigmaPPSystFeedDown, TH1F* &hReference, TH1F* &hReferenceSystData, TH1F* &hCombinedReferenceFlag, TGraphAsymmErrors* &gReference, TGraphAsymmErrors* &gReferenceSyst, TGraphAsymmErrors* &gReferenceFdSyst, std::vector<std::string> scalesToExclude, bool useFDUnc) {

    const int nPtBinsExtr = ptLimsExtr.size() - 1;
    const int nPtBinsNew = nPtBinsMeas + nPtBinsExtr;
    double ptLimsNew[nPtBinsNew+1];
    for(int iPt=0; iPt<nPtBinsMeas; iPt++)
        ptLimsNew[iPt] = ptLimsMeas[iPt];
    for(int iPt=0; iPt<=nPtBinsExtr; iPt++)
        ptLimsNew[nPtBinsMeas+iPt] = ptLimsExtr[iPt];

    hReference = new TH1F("hReference", "", nPtBinsNew, ptLimsNew);
    hReference->SetDirectory(0);
    hReferenceSystData = new TH1F("hReferenceSystData", "", nPtBinsNew, ptLimsNew);
    hReferenceSystData->SetDirectory(0);
    hCombinedReferenceFlag = new TH1F("hCombinedReferenceFlag", "", nPtBinsNew, ptLimsNew);
    hCombinedReferenceFlag->SetDirectory(0);
    gReference = new TGraphAsymmErrors(nPtBinsNew);
    gReference->SetName("gReference");
    gReferenceSyst = new TGraphAsymmErrors(nPtBinsNew);
    gReferenceSyst->SetName("gReferenceSyst");
    if(useFDUnc)
    {
        gReferenceFdSyst = new TGraphAsymmErrors(nPtBinsNew);
        gReferenceFdSyst->SetName("gReferenceFdSyst");
    }

    std::cout << "\nComputing extrapolated pT bins:" << std::endl;
    for(int iPt=0; iPt<nPtBinsExtr; iPt++) {
        std::cout << ptLimsExtr[iPt] << "  " << ptLimsExtr[iPt+1] << std::endl;
    }
    std::cout << "File with FONLL input for extrapolation: " << FONLLFileName << "\n" << std::endl;

    std::map<std::string, int> colors = {{"central", kBlack}, {"min", kGray+1}, {"max", kGray+3}, {"min_sc", kBlue+2}, {"max_sc", kBlue}, {"min_mass", kAzure+4}, {"max_mass", kCyan+3},
                                         {"min_pdf", kTeal+4}, {"max_pdf", kGreen+2}, {"frdot5dot5", kOrange-2}, {"fr22", kOrange+7}, {"fr21", kRed-7}, {"fr12", kRed}, {"fr1dot5", kRed+1}, {"frdot51", kViolet+9}};

    std::map<std::string, TH1F*> hFONLL = ReadFONLL(FONLLFileName);
    std::map<std::string, TH1F*> hFONLLReb, hRatioFONLL;
    std::map<std::string, TF1*> fRatioFONLL;
    double minRatio=1.e9, maxRatio=-1;
    for(auto &el : hFONLL)
    {
        if(scalesToExclude.size() > 0)
        {
            if(std::find(scalesToExclude.begin(), scalesToExclude.end(), el.first) != scalesToExclude.end())
                continue;
        }

        hFONLLReb[el.first] = static_cast<TH1F*>(el.second->Rebin(nPtBinsNew, Form("%s_reb", el.second->GetName()), ptLimsNew));
        hFONLLReb[el.first]->Scale(1., "width");
        hRatioFONLL[el.first] = static_cast<TH1F*>(hSigmaPP->Clone(Form("%s_ratio", el.second->GetName())));
        hRatioFONLL[el.first]->SetStats(0);
        hRatioFONLL[el.first]->SetLineWidth(2);
        hRatioFONLL[el.first]->SetLineColor(colors[el.first]);
        hRatioFONLL[el.first]->SetMarkerColor(colors[el.first]);
        for(int iPt=0; iPt<nPtBinsMeas; iPt++)
        {
            hRatioFONLL[el.first]->SetBinContent(iPt+1, hSigmaPP->GetBinContent(iPt+1)/hFONLLReb[el.first]->GetBinContent(iPt+1));
            hRatioFONLL[el.first]->SetBinError(iPt+1, hSigmaPP->GetBinError(iPt+1)/hFONLLReb[el.first]->GetBinContent(iPt+1));
        }
        fRatioFONLL[el.first] = new TF1(Form("fFONLL%s_ratio", el.first.data()), "pol0", 0., 100.);
        fRatioFONLL[el.first]->SetLineWidth(2);
        fRatioFONLL[el.first]->SetLineColor(colors[el.first]);

        hRatioFONLL[el.first]->Fit(fRatioFONLL[el.first], "Q0", "", 5., ptLimsMeas[nPtBinsMeas]);
        if(fRatioFONLL[el.first]->GetParameter(0) > maxRatio)
            maxRatio = fRatioFONLL[el.first]->GetParameter(0);
        if(fRatioFONLL[el.first]->GetParameter(0) < minRatio)
            minRatio = fRatioFONLL[el.first]->GetParameter(0);
    }

    TAxis* ptAxisMeas = hSigmaPP->GetXaxis();
    TAxis* ptAxisExtr = hFONLLReb["central"]->GetXaxis();
    for(int iPt=0; iPt<nPtBinsNew; iPt++)
    {
        int ptBinMeas = ptAxisMeas->FindBin(hReference->GetBinCenter(iPt+1));
        if(ptBinMeas!=0 && ptBinMeas!=ptAxisMeas->GetNbins()+1)
        {
            int ptBinMeasGraph = -1;
            for(int iPtGraph=0; iPtGraph<gSigmaPPSystData->GetN(); iPtGraph++)
            {
                double ptMeasGraph = -1, sigma = -1;
                gSigmaPPSystData->GetPoint(iPtGraph, ptMeasGraph, sigma);
                if(TMath::Abs(ptMeasGraph-hReference->GetBinCenter(iPt+1)) < 0.01)
                {
                    ptBinMeasGraph = iPtGraph;
                    break;
                }
            }
            double systdatalow = gSigmaPPSystData->GetErrorYlow(ptBinMeasGraph);
            double systdatahigh = gSigmaPPSystData->GetErrorYhigh(ptBinMeasGraph);
            double systFDlow = 0.;
            double systFDhigh = 0.;
            if(useFDUnc) 
            {
                systFDlow = gSigmaPPSystFeedDown->GetErrorYlow(ptBinMeasGraph);
                systFDhigh = gSigmaPPSystFeedDown->GetErrorYhigh(ptBinMeasGraph);
            }

            double systtotlow = gSigmaPPSyst->GetErrorYlow(ptBinMeasGraph);
            double systtothigh = gSigmaPPSyst->GetErrorYhigh(ptBinMeasGraph);

            hReference->SetBinContent(iPt+1, hSigmaPP->GetBinContent(ptBinMeas));
            hReference->SetBinError(iPt+1, hSigmaPP->GetBinError(ptBinMeas));
            hReferenceSystData->SetBinContent(iPt+1, hSigmaPP->GetBinContent(ptBinMeas));
            hReferenceSystData->SetBinError(iPt+1, systdatalow);
            hCombinedReferenceFlag->SetBinContent(iPt+1, 0.);

            gReference->SetPoint(iPt+1, hSigmaPP->GetBinCenter(ptBinMeas), hSigmaPP->GetBinContent(ptBinMeas));
            gReference->SetPointError(iPt+1, 0.1, 0.1, hReference->GetBinError(iPt+1), hReference->GetBinError(iPt+1));
            gReferenceSyst->SetPoint(iPt+1, hSigmaPP->GetBinCenter(ptBinMeas), hSigmaPP->GetBinContent(ptBinMeas));
            gReferenceSyst->SetPointError(iPt+1, 0.1, 0.1, systtotlow, systtothigh);
            if(useFDUnc)
            {
                gReferenceFdSyst->SetPoint(iPt+1, hSigmaPP->GetBinCenter(ptBinMeas), hSigmaPP->GetBinContent(ptBinMeas));
                gReferenceFdSyst->SetPointError(iPt+1, 0.1, 0.1, systFDlow, systFDhigh);
            }
        }
        else
        {
            int ptBinExtr = ptAxisExtr->FindBin(hReference->GetBinCenter(iPt+1));
            double centvalue = hFONLLReb["central"]->GetBinContent(ptBinExtr)*fRatioFONLL["central"]->GetParameter(0);
            double stat = hFONLLReb["central"]->GetBinContent(ptBinExtr)*fRatioFONLL["central"]->GetParError(0);
            double systextrlow = hFONLLReb["central"]->GetBinContent(ptBinExtr)*(fRatioFONLL["central"]->GetParameter(0)-minRatio);
            double systextrhigh = hFONLLReb["central"]->GetBinContent(ptBinExtr)*(maxRatio-fRatioFONLL["central"]->GetParameter(0));

            hReference->SetBinContent(iPt+1, centvalue);
            hReference->SetBinError(iPt+1, stat);
            hReferenceSystData->SetBinContent(iPt+1, centvalue);
            hReferenceSystData->SetBinError(iPt+1, 0.);
            hCombinedReferenceFlag->SetBinContent(iPt+1, 1.);

            gReference->SetPoint(iPt+1, hReference->GetBinCenter(iPt+1), centvalue);
            gReference->SetPointError(iPt+1, 0.1, 0.1, hReference->GetBinError(iPt+1), hReference->GetBinError(iPt+1));
            gReferenceSyst->SetPoint(iPt+1, hReference->GetBinCenter(iPt+1), centvalue);
            gReferenceSyst->SetPointError(iPt+1, 0.1, 0.1, systextrlow, systextrhigh);
            if(useFDUnc) {
                gReferenceFdSyst->SetPoint(iPt+1, hReference->GetBinCenter(iPt+1), centvalue);
                gReferenceFdSyst->SetPointError(iPt+1, 0.1, 0.1, 0., 0.);
            }
        }
    }

    TCanvas* cPtExtrap = new TCanvas("cPtExtrap", "", 800, 800);
    TLegend* legFONLL = new TLegend(0.5, 0.4, 0.8, 0.9);
    legFONLL->SetFillStyle(0);
    legFONLL->SetTextSize(0.04);
    cPtExtrap->DrawFrame(ptLimsMeas[0], 0., ptLimsMeas[nPtBinsMeas], maxRatio*4, ";#it{p}_{T} (GeV/#it{c});meas/FONLL (a.u)");
    for(auto &el : colors)
    {
        if(scalesToExclude.size() > 0)
        {
            if(std::find(scalesToExclude.begin(), scalesToExclude.end(), el.first) != scalesToExclude.end())
                continue;
        }

        hRatioFONLL[el.first]->DrawCopy("same");
        fRatioFONLL[el.first]->Draw("same");
        legFONLL->AddEntry(hRatioFONLL[el.first], el.first.data(), "lp");
    }
    legFONLL->Draw();
}

//________________________________________________________________________
TGraphAsymmErrors* ComputeCorrForYshift(std::string FONLLfileNameNum, std::string FONLLfileNameDen, const int nPtBinsMeas, const double* ptLimsMeas) {

    std::cout << "\nComputing y shift correction:" << std::endl;
    std::cout << "Files with FONLL input:" << std::endl;
    std::cout << "\tNumerator: " << FONLLfileNameNum << std::endl;
    std::cout << "\tDenominator: " << FONLLfileNameDen << "\n" << std::endl;

    std::map<std::string, int> colors = {{"central", kBlack}, {"min", kGray+1}, {"max", kGray+3}, {"min_sc", kBlue+2}, {"max_sc", kBlue}, {"min_mass", kAzure+4}, {"max_mass", kCyan+3},
                                         {"min_pdf", kTeal+4}, {"max_pdf", kGreen+2}, {"frdot5dot5", kOrange-2}, {"fr22", kOrange+7}, {"fr21", kRed-7}, {"fr12", kRed}, {"fr1dot5", kRed+1}, {"frdot51", kViolet+9}};

    std::map<std::string, TH1F*> hFONLLNum = ReadFONLL(FONLLfileNameNum, "num");
    std::map<std::string, TH1F*> hFONLLDen = ReadFONLL(FONLLfileNameDen, "den");
    std::map<std::string, TH1F*> hFONLLDenReb, hFONLLNumReb, hRatioFONLL;
    for(auto &el : hFONLLDen)
    {
        hFONLLNumReb[el.first] = static_cast<TH1F*>(hFONLLNum[el.first]->Rebin(nPtBinsMeas, Form("%s_reb", hFONLLNum[el.first]->GetName()), ptLimsMeas));
        hFONLLDenReb[el.first] = static_cast<TH1F*>(hFONLLDen[el.first]->Rebin(nPtBinsMeas, Form("%s_reb", hFONLLDen[el.first]->GetName()), ptLimsMeas));
        hRatioFONLL[el.first] = static_cast<TH1F*>(hFONLLNumReb[el.first]->Clone(Form("%s_ratio", hFONLLNum[el.first]->GetName())));
        hRatioFONLL[el.first]->Divide(hFONLLDenReb[el.first]);
        hRatioFONLL[el.first]->SetStats(0);
        hRatioFONLL[el.first]->SetLineWidth(2);
        hRatioFONLL[el.first]->SetFillStyle(0);
        hRatioFONLL[el.first]->SetLineColor(colors[el.first]);
    }

    TGraphAsymmErrors* gScal = new TGraphAsymmErrors(0);
    for(int iPt=0; iPt<hRatioFONLL["central"]->GetNbinsX(); iPt++)
    {
        double pt = hRatioFONLL["central"]->GetBinCenter(iPt+1);
        double ptUnc = hRatioFONLL["central"]->GetBinWidth(iPt+1)/2.;
        double ratio = hRatioFONLL["central"]->GetBinContent(iPt+1);
        double ratioMin = ratio;
        double ratioMax = ratio;
        for(auto &el : hRatioFONLL)
        {
            double ratioVar = hRatioFONLL[el.first]->GetBinContent(iPt+1);
            if(ratioVar < ratioMin) ratioMin = ratioVar;
            if(ratioVar > ratioMax) ratioMax = ratioVar;
        }
        double uncRatioLow = ratio - ratioMin;
        double uncRatioHigh = ratioMax - ratio;
        gScal->SetPoint(iPt, pt, ratio);
        gScal->SetPointError(iPt, ptUnc, ptUnc, uncRatioLow, uncRatioHigh);
    }
    gScal->SetFillStyle(1000);
    gScal->SetLineColor(kBlack);
    gScal->SetName("gYcorr");

    TCanvas* cRatioY = new TCanvas("cRatioY", "", 800, 800);
    TLegend* legFONLL = new TLegend(0.2, 0.65, 0.8, 0.9);
    legFONLL->SetFillStyle(0);
    legFONLL->SetNColumns(2);
    legFONLL->SetTextSize(0.04);
    cRatioY->DrawFrame(ptLimsMeas[0], 0.9, ptLimsMeas[nPtBinsMeas], 1.3, ";#it{p}_{T} (GeV/c);FONLL rapidity ratio");
    for(auto &el : hRatioFONLL)
    {
        hRatioFONLL[el.first]->Draw("same");
        legFONLL->AddEntry(hRatioFONLL[el.first], el.first.data(), "lp");
    }
    gScal->Draw("2");
    legFONLL->Draw();

    return gScal;
}

//________________________________________________________________________
void ApplyYshiftToHisto(TGraphAsymmErrors* gyShift, TH1F* histo) {

  for(int iPt=0; iPt<histo->GetNbinsX(); iPt++) {
        double shift = -1.;
        for(int iPtShift=0; iPtShift<gyShift->GetN(); iPtShift++) {
            double ptcent = -1.;
            gyShift->GetPoint(iPtShift, ptcent, shift);
            if(TMath::Abs(histo->GetBinCenter(iPt+1)-ptcent)<0.001)
                break;
        }
        histo->SetBinContent(iPt+1, histo->GetBinContent(iPt+1)*(1./shift));
        histo->SetBinError(iPt+1, histo->GetBinError(iPt+1)*(1./shift));
    }
}

//________________________________________________________________________
void ApplyYshiftToGraph(TGraphAsymmErrors* gyShift, TGraphAsymmErrors* graph, bool useUnc) {

  for(int iPt=0; iPt<graph->GetN(); iPt++) {
        double pt = -1., sigma = -1.;
        graph->GetPoint(iPt, pt, sigma);
        double errlow = graph->GetErrorYlow(iPt);
        double errhigh = graph->GetErrorYhigh(iPt);
        double shift = -1., shiftunclow = 0., shiftunchigh = 0.;

        for(int iPtShift=0; iPtShift<gyShift->GetN(); iPtShift++) {
            double ptcent = -1.;
            gyShift->GetPoint(iPtShift, ptcent, shift);
            shiftunclow = gyShift->GetErrorYlow(iPtShift);
            shiftunchigh = gyShift->GetErrorYhigh(iPtShift);
            if(TMath::Abs(pt-ptcent)<0.001)
                break;
        }

        graph->SetPoint(iPt,pt,sigma*shift);
        if(!useUnc)
            graph->SetPointError(iPt, 0.1, 0.1, errlow*(1./shift), errhigh*(1./shift));
        else 
        {
            double errtotlow = TMath::Sqrt(errlow/sigma*errlow/sigma+shiftunclow/shift*shiftunclow/shift)*sigma*(1./shift);
            double errtothigh = TMath::Sqrt(errhigh/sigma*errhigh/sigma+shiftunchigh/shift*shiftunchigh/shift)*sigma*(1./shift);
            graph->SetPointError(iPt, 0.1, 0.1, errtotlow, errtothigh);
        }
    }
}

//________________________________________________________________________
std::map<std::string, TH1F*> ReadFONLL(std::string FONLLFileName, std::string suffix) {

    std::map<std::string, TH1F*> hFONLL;

    if(FONLLFileName.find("txt") == std::string::npos && FONLLFileName.find("dat") == std::string::npos && FONLLFileName.find("csv") == std::string::npos) {
        std::cerr << "ERROR: Wrong file format! Exit." << std::endl;
        return hFONLL;
    }

    std::ifstream inSet(FONLLFileName.data());
    if(!inSet) {
        std::cerr << "ERROR: Please check if "<< FONLLFileName.data() << " is the right path. Exit." << std::endl;
        return hFONLL;
    }

    std::map<int, std::string> varpos = {{1, "pt"}, {2, "central"}, {3, "min"}, {4, "max"}, {5, "min_sc"}, {6, "max_sc"}, {7, "min_mass"}, {8, "max_mass"},
                                         {9, "min_pdf"}, {10, "max_pdf"}, {11, "frdot5dot5"}, {12, "fr22"}, {13, "fr21"}, {14, "fr12"}, {15, "fr1dot5"}, {16, "frdot51"}};
    std::map<std::string, std::vector<double> > vars;
    for(auto &el : varpos)
        vars[el.second] = {};

    while(!inSet.eof())
    {    
        std::string line;
        std::getline(inSet, line);
        if(line.find("#") != std::string::npos) // do not read commented lines
            continue;

        std::stringstream convert(line);
        std::string sval;
        int iPos = 0;
        while(std::getline(convert, sval, ' '))
        {
            if(std::atof(sval.data()) > 0)
                iPos++;
            if(iPos>0)
                vars[varpos[iPos]].push_back(std::atof(sval.data()));
        }
    }

    unsigned int nPtBins = vars["pt"].size();
    double ptWidth = (vars["pt"][1]-vars["pt"][0])/2;
    double ptMin = vars["pt"][0]-ptWidth;
    double ptMax = vars["pt"][nPtBins-1]+ptWidth;

    for(auto &el : varpos)
    {
        if(el.second.compare("pt") != 0)
            hFONLL[el.second] = new TH1F(Form("hFONLL%s%s", el.second.data(), suffix.data()), "", nPtBins, ptMin, ptMax);
    }

    for(unsigned int iPt=0; iPt<nPtBins; iPt++)
    {
        for(auto &el : varpos)
        {
            if(el.second.compare("pt") != 0)
                hFONLL[el.second]->SetBinContent(iPt+1, vars[el.second][iPt]);
        }
    }

    inSet.close();
    return hFONLL;
}

//________________________________________________________________________
void SetStyle() {
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadRightMargin(0.035);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetTitleSize(0.045, "xy");
    gStyle->SetLabelSize(0.04, "xy");
    gStyle->SetTitleOffset(1.2,"y");
    gStyle->SetTitleOffset(1.2,"x");
    gStyle->SetLegendBorderSize(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    TGaxis::SetMaxDigits(3);
}
