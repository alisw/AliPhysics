#ifndef BADCHANNELANA_H
#define BADCHANNELANA_H

/// \class BadChannelAna
/// \ingroup EMCALPerformanceMacros
/// \brief Analyses cell properties and identifies bad cells
///
/// This is used for bad channel identification in EMCal and DCal.
/// It maks bad channels in several energy intervals and plots output with various
/// information for cross checks.
///
/// The class builds a mean distribution of certain cell observables
/// and compares single cell properties to this mean.
/// That way bad channels (far off the mean) are identified and flagged.
/// A .pdf file with their spectra is created. This should be
/// cross checked by hand.
///
/// \author Eliane Epple <eliane.epple@yale.edu>, Yale University
/// \author Chiara Bianchin <chiara.bianchin@cern.ch>, Wayne State University
/// based on the work from:
/// \author Alexis Mas <aleximas@if.usp.br>, SUBATECH
/// and
/// \author Marie Germain <Marie.Germain@subatech.in2p3.fr>, SUBATECH
/// which is in turn based on getCellsRunQA.C from
/// \author Olga Driga, SUBATECH

/// \date June 29, 2017

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice
*/
 
 
#include <Riostream.h>
#include <TString.h>
#include <TArrayD.h>

class TH1;
class TH1F;
class TH1D;
class TH2F;
class TFile;
class TList;

class AliCalorimeterUtils;


class BadChannelAna : public TObject {
	
public:
      BadChannelAna() ;                // default ctor
	  ~BadChannelAna();                // dtor
	  BadChannelAna(TString period, TString train, TString trigger, Int_t runNumber,Int_t trial, TString workDir, TString listName, Bool_t runByRun=0);

	  void Run(Bool_t mergeOnly=0);

      //Setters
	  void SetExternalMergedFile(TString inputName)        {fExternalFileName   = inputName;}
	  void SetExternalBadMap(TString inputName)            {fExternalBadMapName = inputName;}
	  void SetQAChecks(Bool_t inputBool)                   {fTestRoutine        = inputBool;}
      void SetPrintOutput(Bool_t inputBool)                {fPrint              = inputBool;}
      void SetTrackCellRecord(Bool_t inputBool)            {fTrackCellRecord    = inputBool;}
      void SetStartEndCell(Int_t start, Int_t end)         {fStartCell = start; fNoOfCells = end;}
      void SetLowerBound(Double_t input)                   {fEndLowerBound      = input;}
      void AddManualMasking(std::vector<Int_t> cellVector) {fManualMask.swap(cellVector) ;}
	  void AddMaskSM(Int_t iSM);
	  void RunMaskSM();
	  void AddPeriodAnalysis(Int_t criteria, Double_t nsigma, Double_t emin, Double_t emax);

protected:

	  void Init();
	  TString MergeRuns();
	  void LoadExternalBadMap();
	  void BCAnalysis();
	  void PeriodAnalysis(Int_t criterum=7, Double_t nsigma = 4.0, Double_t emin=0.1, Double_t emax=2.0);

	  TH1F* BuildHitAndEnergyMean(Int_t crit, Double_t emin = 0.1, Double_t emax=2.);
	  TH1F* BuildHitAndEnergyMeanScaled(Int_t crit, Double_t emin = 0.1, Double_t emax=2.);
	  TH1F* BuildTimeMean(Int_t crit, Double_t tmin, Double_t tmax);

	  void FlagAsDead();
	  void FlagAsBad(Int_t crit, TH1F* inhisto, Double_t nsigma = 4., Double_t dnbins = 200);
	  void FlagAsBad_Time(Int_t crit, TH1F* inhisto, Double_t nSig = 3);

	  void SummarizeResultsByFlag();
	  void SummarizeResults();
	  TH1D *BuildMeanFromGood(Int_t warmIn=0);
	  Bool_t CheckDistribution(TH1* ratio, TH1* reference);
	  Bool_t IsCoveredByTRD(Int_t row, Int_t collumn);
	  void SaveBadCellsToPDF(Int_t version, TString pdfName);
	  void PlotFlaggedCells2D(Int_t flagBegin,Int_t flagEnd=-1);
	  void SaveHistoToFile();
	  void PrintCellInfo(Int_t number);

	  //Settings for analysed period
	  Int_t   fCurrentRunNumber;            ///< A run number of an analyzed period. This is important for the AliCalorimeterUtils initialization
      TString fPeriod;                      ///< The name of the analyzed period
	  TString fTrainNo;                     ///< Train number of the analyszed data (can deduce pass & trigger from that etc.)
	  TString fTrigger;                     ///< Selected trigger for the analysis
	  Int_t   fNoOfCells;                   ///< Number of cells in EMCal and DCal
	  Int_t   fCellStartDCal;               ///< ID of the first cell in the DCal
	  Int_t   fStartCell;                   ///< ID of the first cell you want to check
	  Int_t   fStartCellSM[21];             ///< CellIDs of first cell in the  20SMs plus last cell ID
	  Double_t fEndLowerBound;              ///< Lower bound

	  //Genergal paths
	  TString fAnalysisOutput;              ///< The list with bad channels and histograms are saved in this folder
	  TString fAnalysisOutputRbR;           ///< For a compact summary of true run-by-run BC maps
	  TString fAnalysisInput;               ///< Here the .root files of each run of the period are saved
	  TString fRunList;                     ///< Thats the full path and name of the file which contains a list of all runs to be merged together
	  TString fRunListFileName;             ///< This is the name of the file with the run numbers to be merged, by default it's 'runList.txt'
	  TString fWorkdir;                     ///< Directory which contains the folders fMergeOutput, fAnalysisInput and fAnalysisOutput. By default it is './'

	  //
	  TString fQADirect;                    ///< Dierctory in the QA.root files where the input histograms are stored
	  TString fMergedFileName;              ///< Filename of the .root file containing the merged runs
	  std::vector<TArrayD> fAnalysisVector; ///< Vector of analysis information. Each place is filled with 4 doubles: version, sigma, lower, and upper energy range

	  //Things to be individualized by setters
	  Int_t   fTrial;                       ///< Number of trial that this specific analyis is. By default '0' so one can try different settings without overwriting the outputs
	  TString fExternalFileName;            ///< If you have already a file that contains many runs merged together you can place it in fMergeOutput and set it with SetExternalMergedFile(FileName)
	  TString fExternalBadMapName;          ///< Load an external bad map to test the effect on block or a given run
	  Bool_t  fTestRoutine;                 ///< This is a flag, if set true will produce some extra quality check histograms
      Bool_t  fPrint;                       ///< If set true more couts with information of the excluded cells will be printed
	  Bool_t  fTrackCellRecord;             ///< Track the non-zero elements in the flags throughout the routine
	  Bool_t  fRunBRunMap;                  ///< Produce truely run-by-run maps in a separate folder

	  //histogram settings
	  Int_t   fNMaxCols;                    ///< Maximum No of colums in module (eta direction)
	  Int_t   fNMaxRows;                    ///< Maximum No of rows in module   (phi direction)
	  Int_t   fNMaxColsAbs;                 ///< Maximum No of colums in Calorimeter
	  Int_t   fNMaxRowsAbs;                 ///< Maximum No of rows in Calorimeter

	  //arrays to store information
	  Double_t fnEventsInRange;
	  std::vector<Int_t> fFlag;             //!<! fFlag[CellID] = 0 (ok),1 (dead),2 and higher (bad certain criteria) start at 0 (cellID 0 = histobin 1)
	  Int_t fCriterionCounter;              //!<! This value will be written in fflag and updates after each PeriodAnalysis, to distinguish the steps at which cells are marked as bad
	  std::vector<Bool_t> fWarmCell;        //!<! fWarmCell[CellID] = 0 (really bad), fWarmCell[CellID] = 1 (candidate for warm),
	  std::vector<Int_t> fManualMask;       //!<! Is a list of cells that should be addidionally masked by hand.
	  std::vector<Int_t> fSMMask;           //!<! fSMMask is filled with SM numbers that need to be masked by hand.

	  //Calorimeter information for the investigated runs
	  AliCalorimeterUtils* fCaloUtils;      //!<! Calorimeter information for the investigated runs

	  TFile* fRootFile;                     //!<! root file with all histograms from this analysis
	  TH2F* fCellAmplitude;                 //!<! main histogram for the analysis. Cell ID vs. amplitude, read from the input merged file
	  TH2F* fCellTime;                      //!<! possible histogram for the analysis. Cell ID vs. time, read from the input merged file
	  TH1F* fProcessedEvents;               //!<! Stores the number of events in the run
	  TH1F* fhCellFlag;                     //!<! histogram that stores by which flag the cell has been excluded
	  TH1F* fhCellWarm;                     //!<! histogram that stores whether the cell was marked as warm

	  TList* fOutputListBad;                //!<! list with bad channel amplitudes, stored in fRootFile
	  TList* fOutputListBadRatio;           //!<! list with bad channel amplitude ratios, stored in fRootFile
	  TList* fOutputListGood;               //!<! list with good channel amplitudes, stored in fRootFile
	  TList* fOutputListGoodRatio;          //!<! list with good channel amplitude ratios, stored in fRootFile
	  
private:
	  BadChannelAna           (const BadChannelAna&); // not implemented
	  BadChannelAna &operator=(const BadChannelAna&); // not implemented

	  TH1F* fAvgNHitPerEvVsCellId;          //!<! being discussed
	  TH1F* fAvgEngPerHitVsCellId;          //!<! being discussed
	  
	
	/// \cond CLASSIMP
	ClassDef(BadChannelAna, 2);
	/// \endcond
};
#endif
