#ifndef ALIT0TIMEAMPLCORR_H
#define ALIT0TIMEAMPLCORR_H

/* Copyright(c) 1998-20016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/**************************************************************************
 * Time - Amplitude correction functions for T0 detector
 * Create correction graphs for PMTs by ESD data trees\
 *
 * 03.02.2016
 * Dmitry.Finogeev@cern.ch
**************************************************************************/

#include<TSystem.h>
#include<TSystemDirectory.h>
#include<TSystemFile.h>
#include<TList.h>
#include<TObjArray.h>
#include<TF1.h>
#include<TH1F.h>
#include<TH2I.h>
#include<TProfile.h>
#include<TString.h>
#include<TChain.h>
#include<TMath.h>
#include<TFile.h>
#include<TSpectrum.h>
#include<TArray.h>
#include<TMath.h>
#include<TGraph.h>
#include<TList.h>
#include<TDirectory.h>
#include<TStopwatch.h>

#include<TGrid.h>
//#include"AliCDBManager.h"

#define NPMT0 24 //number T0 photomultipliers
#define MAX_CHAR_LINE 500 //max char line length

class AliT0TimeAmplCorr
{
public:
    AliT0TimeAmplCorr();
    AliT0TimeAmplCorr(AliT0TimeAmplCorr &source){ Copy(source); }
    AliT0TimeAmplCorr& operator= ( AliT0TimeAmplCorr &source ){ Copy(source); return *this;}

    virtual ~AliT0TimeAmplCorr();

    /*!
     * \brief Set Verbose mode
     * \param verbose 0 - nothing; 1 - errors; 2 - brief status; 3 - full status; 4 - calculated values; 5 - debugging
     */
    void SetVerbose(const Int_t verbose = 1){fVerbose = verbose;}

    /*!
     * \brief Set tree type
     * \param option "RAW", "ESD"
     * \return 1 if successfully else 0
     */
    Int_t SetTreeType(TString option = "RAW");

    /*!
     * \brief Set RUN number
     * \param numrun Run number
     * \return 1 if successfully else 0
     */
    Int_t SetNumRUN(const Int_t numrun);

    /*!
     * \brief Set QTC type
     * \param QTCtype "NEW", "OLD"
     * \return 1 if successfully else 0
     */
    Int_t SetQTCType(TString QTCtype = "OLD");

    void SetRangeVertexSPD(const Int_t rangevertesSPD = 1){fVertexRange = rangevertesSPD;}
    void SetRangeTVDC(const Int_t rangetvdc = 800){frangeTVDC = rangetvdc;}
    void SetRangeCFD(const Int_t rangecfd = 50){frangeCFD = rangecfd;}
    void SetRangeQT1(const Int_t rangeqt1 = 800){frangeQT1 = rangeqt1;}
    void SetRangeNewQTC(const Int_t rangemin = 18000, const Int_t rangemax = 19000){frangeNewQTCmin = rangemin; frangeNewQTCmax = rangemax;}

    /*!
     * \brief Take TTree "treepathname", from file "filename" if tree is in folder
     * treepathname is the path to tree, if tree is in root, "treepathname" is tne tree name
     * tree type must be set with SetTreeType()
     * \param filename File name
     * \param treename TTree name saved in file
     * \return 1 if successfully else 0
     */
    Int_t ConnectSourceFile(TString filename, TString treepathname = "");

    /*!
     * \brief Take TTree "treename" from files
     * file name is configured by mask format as "/tree/treeRUN%i/tree_file_%i.root"\n
     * mask must include two "%i" for numrun and file index respectively\n
     * file index varies from ffirst to flast\n
     * \param filemask example: "/tree/treeRUN%i/tree_file_%i.root"
     * \param numrun
     * \param treename
     * \param ffirst
     * \param flast
     * \return 1 if successfully else 0
     */
    Int_t ConnectSourceFilesByMask(TString filemask, const Int_t numrun, TString treename = "", Int_t ffirst = 0, Int_t flast = 1000);

    /*!
     * \brief Take TTree "treename" from files in directory.
     * File name must contains oblnamepart
     * \param path
     * \param treename
     * \param oblnamepart
     * \return 1 if successfully else 0
     */
    Int_t ConnectSourceFilesFromDir(TString path = "", TString treename = "", TString oblnamepart = "");

    /*!
     * \brief Take from file coefficients for new QTC data
     * RUN number must be set with SetNumRUN(), by default or in case of failure set default values
     * set QTC type as "NEW", so does not need to call SetQTCType("NEW")
     * \param filename
     * \return 1 if successfully else 0
     */
    Int_t SetNewQTCcoefficients(TString filename = "");

    /*!
     * \brief Take QTCCFD plots from root file and merge them to already loaded
     * in case plots have a name "plot_for_pmt1", where 1 is PMT number 1-24,
     * use PlotNameMask = "plot_for_pmt%i"
     * \param filename
     * \return 1 if successfully else 0
     */
    Int_t MergeQTCCFDplots(TString filename, TString PlotNameMask = "");

    /*!
     * \brief MergeQTCCFDplotsFromDir merge all plot files from dirrectory
     * if(PlotNameMask != "") QTCname = Form(PlotNameMask, PMT+1)
     * \param filename
     * \param PlotNameMask
     * \return
     */
    Int_t MergeQTCCFDplotsFromDir(TString path = "");

    /*!
     * \brief Make QTCCFD plots from tree in connected files
     * existing plots will be replaced, also monitoring histograms CFD, QT1, QTC, CFD_for_1_mip will be created
     * for check. Range of fit function in histograms show range of selection for creating QTCCFD plots.
     * Two procedures for mean values may be done: get values from ocdb and find values by fitting histograms
     * \param option
     * "noocdb" - do not use OCD base
     * "nomeanfit" - do not make fitting for mean values
     * "fastmeanscan" - use a part of statistic
     * \return 1 if successfully else 0
     */
    Int_t MakeQTCCFDplots(TString option = "");

    /*!
     * \brief Make projection TGraph function for each created or merged plot
     * To create a correction function a 1D histograms is prepared firstly for plot with CreatePojection():\n
     *  1)QTCCFD plot smoothed if rebinoption contains "s"\n
     *  2)QTCCFD plot rebined with rebinX and rebinY values, and mean number of events in bin AverageEvInBinX is founded.\n
     *  3)Bounds bins are determined with AverageEvInBinX*minbinevedge <= (number of events in bound bin).\n
     *  4)if rebinoption contains "r", QTCCFD plot rebinned along X (amplitude) on condition\n
     *    AverageEvInBinX*minbinevX <= (number of events in new bin)\n
     *  5)1D histogram created by mean and RMS values for Y projection for each new bin\n
     * \n
     * After 1D histogram created, it divided by slices and fitted with CreateCorrectionGraph():\n
     *  1)if fragoption contains "p" maximum and minimum peaks found for 1D histogramm\n
     *  2)ranges between peaks or between first and last bins divided in slices so that the weight of each slice\n
     *    corresponded to (total weight of histogram) / maxnslices. Weight may be varied by fragoption:\n
     *    "e"-by errors; "r"-by ranges; "b"-by bins; all containing option taking in account\n
     *  3)each slice fitted by polynomials with various power, and chosen the one with the best chi-square\n
     *  4)if fragoption contains "nofit", correction function plotted by points in hist, no other procedures are done\n
     *  5)Correction function decremented from QTCCFD plot and chi-square by linear function\n
     *    is calculated for assay the correction function with CreateGraphAssay().\n
     * \n
     * QTCCFD plot are much differ from each other and there are not unique option for all, for this reason\n
     * function MakeATCorrectionForPlots(TString configfilename = "") takes range for main options from file\n
     * <configfilename> and by iterating over all options chooses the best by assay. If file is not defended
     * or is not opened, parameters ranges set by default.\n
     * Configuration file example:\n
     * \n
     * [SLICES]\n
     * 5   80	3\n
     * [REBINX]\n
     * 5   40	3\n
     * [EVINBIN]\n
     * 1    10	3\n
     * # Option for projection algorithm:\n
     * # "n" all projection bins are same and equal rebinX\n
     * # "r" - rebining projection with various bin values with evinbin in each\n
     * # "s" - smoothing plot before making projection\n
     * [REBINOPT]\n
     * r rs s\n
     * [FRAGOPT]\n
     * e r b pe pr bp perb\n
     * [PMTS]\n
     * 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23\n
     * # 0	9 10 16\n
     * \n
     * there # - comments, parameter names are placed in []:\n
     * SLICES - maxnslices; REBINX - rebinX; EVINBIN - minbinevX; REBINOPT - rebinoption; FRAGOPT - fragoption; PMTS - pmts numbers\n
     * after SLICES, REBINX, EVINBIN  parameters follow its minimum, maximum values and number of iterations\n
     * after REBINOPT, FRAGOPT, PMTS follow parameters enumeration. Values divided by space or tab\n
     * \n
     * Resulting correction function may be taken with TObjArray* GetGraphs() or wrote in file with SaveResultsInFile();\n
     * While MakeATCorrectionForPlots as well as correction function, created 1D projection histogram, function assay histograms\n
     * and chi-square trends to monitor the chi-square dependence from iterated options, this information also may be take by\n
     * "get" functions or wrote in file with  common function SaveResultsInFile();
     * \n
     * EXAMPLE:\n
     *
     * T0ampltimeCorrection->MergeQTCCFDplots( "/Plots_mip/ATC_Plots_236349.root" );
     * T0ampltimeCorrection->MergeQTCCFDplots( "/Plots_mip/ATC_Plots_236337.root" );
     * T0ampltimeCorrection->MergeQTCCFDplots( "/Plots_mip/ATC_Plots_236246.root" );
     *
     * if( T0ampltimeCorrection->MakeATCorrectionForPlots("./AT_fast_scan.corr") )
     *   T0ampltimeCorrection->SaveResultsInFile( "comment" );
     *
     * \param configfilename
     * \return
     */
    Int_t MakeATCorrectionForPlots(TString configfilename = "");

    /*!
     * \brief Save all Results In File
     * file name consist of parameters markers and comment
     * return constructed file name
     * rewrites existing information in file
     * \param file name comment
     */
    TString SaveResultsInFile(TString filecomment);

    void ClearPlotsArray(){ ResetArray( &fQTCCFDplot, "recreate" ); }


    //getting results array
    const TObjArray* GetMonHists(){return fMonHists;}
    const TObjArray* GetQTCCFDPlots(){return fQTCCFDplot;}
    const TObjArray* GetProjections(){return fProjection;}
    const TObjArray* GetGraphs(){return fGraphs;}
    const TObjArray* GetAssays(){return fAssay;}
    const TObjArray* GetChiTrends(){return fChiTrends;}


protected:
    //processing parameters
    Int_t fVerbose;
    Bool_t fIsNewQTC;
    Bool_t fIsRAWtree;
    Bool_t fIsPlotsLoaded;

    //timer
    TStopwatch *fTimer;

    //data variables
    Int_t fNumRUN;
    TChain *fTreeChain;

    TFile *fTreeFile; // for ESD
    TTree *fTree;     // for ESD

    TObjArray *fMonHists; //array of histograms for monitoring
    TObjArray *fAssay; //array of histograms for monitoring
    TObjArray *fQTCCFDplot; //QTC-CFD projections for fitting
    TObjArray *fProjection; //QTC-CFD projections for fitting
    TObjArray *fGraphs; //resulting function for projections
    TObjArray *fChiTrends; //array contains trends of bests polynomials chi-square for visual checking

    //values ranges variables
    Int_t fCFDPlotMin;
    Int_t fCFDPlotMax;
    Int_t fQTCPlotMin;
    Int_t fQTCPlotMax;

    Int_t fVertexRange;
    Int_t fMeanTVDC;
    Int_t fMeanCFD[NPMT0];//means values found by histogram fitting
    Int_t fMeanCFDmip[NPMT0];
    Int_t fMeanQT1[NPMT0];
    Int_t fPedestalQTC[NPMT0];

    Double_t fNewQTCa[NPMT0];
    Double_t fNewQTCb[NPMT0];

    Int_t frangeCFD;
    Int_t frangeTVDC;
    Int_t frangeQT1;
    Int_t frangeNewQTCmin;
    Int_t frangeNewQTCmax;

    //service functions:

    //copy this class from source exept results
    void Copy(const AliT0TimeAmplCorr &source);

    //delete all elements in array (recreate - creating void new array, delete - delete array)
    void ResetArray(TObjArray **array, TString option = "delete");

    //set mean values to zero
    void ResetMeanValues(){fMeanTVDC = 0; for(Int_t pmt=0;pmt<NPMT0;pmt++) fMeanCFD[pmt] = fMeanCFDmip[pmt] = fMeanQT1[pmt] = fPedestalQTC[pmt] = 0;  }

    //checking is TTree wrote correctly in file
    Int_t IsTTreeWritedCorrectly(TString filename, TString treename);

    //save all TObjArray elements in file in directory hierarchy replacing all existing element in directory
    Int_t SaveArrayHistsInFile(TObjArray *array, TString filename, TString directoryname);

    //bubble sorting array in increasing
    void ArrayBubbleSortInc(Double_t *const array, const Int_t kElements);

    //fit one histogram by gauss, return mean and sigma
    Int_t FitHistogramm(TH1 *const histogramm, Double_t &mean, Double_t &sigma, Double_t rangeXmin = 0, Double_t rangeXmax = 0);

    //call FitHistogramm in region mean +- RMSrange*RMS
    Int_t FitHistogrammInRMSRange(TH1 *const histogramm, Double_t &mean, Double_t &sigma, const Double_t RMSrange = 1);

    //fit histogramm in first max peak
    Int_t FitHistogrammInFirstMaxPeak(TH1*const historgamm, Double_t &mean, Double_t &sigma, const Double_t peaksRangeMin,
                                      const Double_t peaksRangeMax, const Double_t peakRange, const Double_t peakRMSRange);

    //fit one histogram by gauss in range mean+-n*RMS, return mean and sigma
    // old one, must be replaced
    Int_t FitHistogrammInMaxPeak(TH1 *const histogramm, Double_t &mean, Double_t &sigma);


    //subtract function TF1 from TH1, analog of base ROOT function TH1::Add(TF1*, -1, "")
    Int_t SubtractGraphFromHist(TH1 *const histogramm, const TGraph *graph);

    //Finds mean values by fitting histograms and draw histograms for monitoring, if values were taken from OCDB, compare them
    Int_t GetRAWmeanValuesByFIT(const Int_t statistic);

    //Gets mean values from OCDB
    Int_t GetRAWmeanValuesFromOCDB();

    //make plot histogram CFD vs QTC for fitting by 2D histogram
    TH2* CreateCFDQTCplotFromRAW(const Int_t pmt);
    TH2* CreateCFDQTCplotFromESD(const Int_t pmt);

    //make 1D projection with special binning for future fitting
    //minbinevX minimum event in bin along X; minbinevedge minimum events in edges bins (in parts of mean events in bin);
    //chinbinY rebinning along Y; rebinoption: s-smooth 2Dhist, r-with rebinning
    TH1* CreatePojection(TH2*const sourcehist, TString rebinoption, const Float_t kMinbinevX = 5., const Float_t kMinbinevedge = 0.25,
                         const Int_t kRebinX = 1, const Int_t kRebinY = 1);

    //fit projection by number of polynomials; fragopton p-peak search; e-by errors; r-by ranges; b-by bins;
    TGraph* CreateCorrectionGraph(TH1*const sourceprojection, const Int_t kMaxfuncpower=9, TString fitoption="", TString fragoption = "pe", const Int_t kNslices=3);

    //check correction function, return chi-square of linear fitting corrected 2D histogram (result histogram created and must be deleted if need)
    TH1* CreateGraphAssay(TH2 *const sourceplot, Double_t &chisquare, const TGraph *sourcegraph = 0);

    //function for looking over projection and fitting parameters, pointers show to best option
    Double_t TryProjFitOption(TH2 *sourceplot,  TH1 *sourceprojection, TH1 **bestprojection, TGraph **bestgraph, TH1 **bestassayhist, TH1 **CHItrend, Int_t channel,
                           TString projectionoption, TString fragoption, Int_t nslices, Int_t rebinX, Float_t minbinevX);


};


#endif // ALIT0TIMEAMPLCORR_H

