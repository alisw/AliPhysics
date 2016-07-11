#ifndef ALIANACALOCHANNELANALYSIS_H
#define ALIANACALOCHANNELANALYSIS_H

/// \class AliAnaCaloChannelAnalysis
/// \brief Analyses cell properties and identifies bad cells
///
/// This is used for bad channel identification in EMCal and DCal.
/// The class builds a mean distribution of certain cell observables
/// and compares single cell properties to this mean.
/// That way bad channels (far off the mean) are identified and flagged.
/// A .pdf file with their spectra is created. This should be
/// cross checked by hand.
///
/// \author Eliane Epple <eliane.epple@yale.edu>, Yale Univeristy
/// \author Chiara Bianchin <chiara.bianchin@cern.ch>, Wein University
/// based on the work from
/// \author Alexis Mas <aleximas@if.usp.br> & M. Germain <Marie.Germain@subatech.in2p3.fr>, SUBATECH
/// which is in turn based on getCellsRunQA.C from
/// \author Olga Driga, SUBATECH
/// \date Jun 24, 2016

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#if !defined(__CINT__) || defined(__MAKECINT__) 
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TROOT.h>
#include <TArrayD.h>
#include <Riostream.h>
#include <AliCalorimeterUtils.h>
#endif
using namespace std;
using std::vector;
using std::array;



class AliAnaCaloChannelAnalysis : public TObject {

public:

      AliAnaCaloChannelAnalysis() ;                // default ctor
	  virtual ~AliAnaCaloChannelAnalysis()  { ; }  // virtual dtor
	  AliAnaCaloChannelAnalysis(TString period, TString pass, TString trigger, Int_t runNumber);

	  void Run();

    //Setters
	  void SetExternalMergedFile(TString inputName)     {fExternalFileName = inputName;}
	  void SetInputFileList(TString inputName)          {fRunListFileName  = inputName;}
	  void SetWorkDir(TString inputName)                {fWorkdir          = inputName;}
	  void SetNTrial(Int_t inputNr)                     {fTrial            = inputNr  ;}
  
	  void AddPeriodAnalysis(Int_t criteria, Double_t nsigma, Double_t emin, Double_t emax);


protected:

	  void Init();
	  TString Convert();
	  void BCAnalysis();
	  void PeriodAnalysis(Int_t criterum=7, Double_t nsigma = 4.0, Double_t emin=0.1, Double_t emax=2.0);

	  void FlagAsDead();
	  TH1F* TestCellEandN(Int_t crit, Double_t emin = 0.1, Double_t emax=2., Double_t nsigma = 4.);
	  void TestCellShapes(Int_t crit, Double_t fitemin, Double_t fitemax, Double_t nsigma =4.);
	  void Process(Int_t crit, TH1* inhisto, Double_t nsigma = 4., Int_t dnbins = 200, Double_t dmaxval = -1.);

	  void Draw2(Int_t cell);
	  void SaveBadCellsToPDF(Int_t version, TString pdfName);


	  //Settings for analysed period
	  Int_t   fCurrentRunNumber;            ///< A run number of an analyzed period. This is important for the AliCalorimeterUtils initialization
      TString fPeriod;                      ///< The name of the analyzed period
	  TString fPass;                        ///< Pass of the analyzed data
	  TString fTrigger;                     ///< Selected trigger for the analysis
	  Int_t   fNoOfCells;                   ///< Number of cells in EMCal and DCal
	  Int_t   fGoodCellID;                  ///< ID of a good cell to compare the spectra to

	  //Genergal paths
	  TString fMergeOutput;                 ///< Here the merged files of a period are saved for a later analysis
	  TString fAnalysisOutput;              ///< The list with bad channels and histograms are saved in this folder
	  TString fAnalysisInput;               ///< Here the .root files of each run of the period are saved
	  TString fRunList;                     ///< Thats the full path and name of the file which contains a list of all runs to be merged together

	  //
	  TString fQADirect;                    ///< Dierctory in the QA.root files where the input histograms are stored
	  TString fMergedFileName;              ///< Filename of the .root file containing the merged runs
	  std::vector<TArrayD> fAnalysisVector; ///< Vector of analysis information. Each place is filled with 4 doubles: version, sigma, lower, and upper energy range

	  //Things to be individualized by setters
	  TString fRunListFileName;             ///< This is the name of the file with the run numbers to be merged, by default it's 'runList.txt'
	  TString fWorkdir;                     ///< Directory which contains the folders fMergeOutput, fAnalysisInput and fAnalysisOutput. By default it is './'
	  Int_t   fTrial;                       ///< Number of trial that this specific analyis is. By default '0'
	  TString fExternalFileName;            ///< If you have already a file that contains many runs merged together you can place it in fMergeOutput and set it with SetExternalMergedFile(FileName)

	  //arrays to store information
	  Int_t *fFlag;                         //!<! fFlag[CellID] = 0 (ok),1 (dead),2 (bad by lower),3 (bad by upper)     start at 0 (cellID 0 = histobin 1)

	  //histogram settings
	  Int_t fNMaxCols;                      ///< Maximum No of colums in module (eta direction)
	  Int_t fNMaxRows;                      ///< Maximum No of rows in module   (phi direction)
	  Int_t fNMaxColsAbs;                   ///< Maximum No of colums in Calorimeter
	  Int_t fNMaxRowsAbs;                   ///< Maximum No of rows in Calorimeter

	  //Calorimeter information for the investigated runs
	  AliCalorimeterUtils* fCaloUtils;      //!<! Calorimeter information for the investigated runs

private:
	  AliAnaCaloChannelAnalysis           (const AliAnaCaloChannelAnalysis&); // not implemented
	  AliAnaCaloChannelAnalysis &operator=(const AliAnaCaloChannelAnalysis&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnaCaloChannelAnalysis, 1);
  /// \endcond
};
#endif
