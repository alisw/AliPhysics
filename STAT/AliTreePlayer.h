#ifndef ALITREEPLAYER_H
#define ALITREEPLAYER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */


/// \ingroup STAT
/// \class AliTreeFormulaF
/// \brief Helper class for AliTreePlayer - formatted string to query tree
/// ### Example usage
/// #### Example: Construct directory path
/// \code
/// AliExternalInfo info;TTree * tree = info.GetTree("QA.TPC","LHC15o","pass1");
/// formula = new AliTreeFormulaF("xxx","xxx/%d{year}/%d{period.GetName()}/%d{pass.GetName()}/%d{run}/xxx",tree);
/// tree->GetEntry(0);
/// formula->PrintValue(0,0,"")
/// Output: "xxx/2015/LHC15o/pass1/245683/xxx"
/// \endcode
///
/// #### Example: Construct directory path
/// \code
///  AliExternalInfo info;TTree * tree = info.GetTree("QA.TPC","LHC15o","pass1");
///  formula = new AliTreeFormulaF("xxx","<a href=\"https://alice-logbook.cern.ch/logbook/date_online.php?p_cont=rund&p_run=%d{run}\"\a>",tree);
///  tree->GetEntry(0);
///  formula->PrintValue(0,0,"")
///  Output: "<a href=\"https://alice-logbook.cern.ch/logbook/date_online.php?p_cont=rund&p_run=245683\"a>"
/// \endcode
///
/// #### Example: Variable formatting
/// * Internally root formatting  are used for the individual data variables. This formatting is different than standard printf formatting
/// * TODO -  Implement own formatting also for variables -
///\code
/// formula = new AliTreeFormulaF("xxx","<a href=\"https://alice-logbook.cern.ch/logbook/date_online.php?p_cont=rund&p_run=%2.1f{run}\"\a>",tree);
/// Ouptut: "<a href=\"https://alice-logbook.cern.ch/logbook/date_online.php?p_cont=rund&p_run=245683.0\"a>"
///\endcode
///\code
/// formula = new AliTreeFormulaF("xxx","<a href=\"https://alice-logbook.cern.ch/logbook/date_online.php?p_cont=rund&p_run=%0.5f{run}\"\a>",tree);
/// Output: "<a href=\"https://alice-logbook.cern.ch/logbook/date_online.php?p_cont=rund&p_run=245683.00000\"a>"
/// \endcode
/// \author Marian Ivanov
///
/// \class AliTreePlayer
/// \brief Set of functions to extend functionality of the TTreePlayer
/// * Input data sources  - TTree Friend trees
/// * Data source can be configured independently using class AliExternalInfo and subdet of fuctionality in the TStatToolkit
///   * friend relation
///   * indices
///   * metatadata for branches,aliases (class, Axis,AxisTitle,Legen, Html, description ...)
/// ### Functionality
///
/// * function to support metadata and collumn annotation
/// * filtering function (branch, aliases metadata)
/// * select function  - export to (json, csv,  html, JIRA) + metadata (not yet implemented)
///   * subset of metadata for columns included in select export outputName.metadata
///         * select function  - for root trees
///   * to be done using TTree::CopyTree  functionality switching ON/OFF selected branches - metadata should be exported automatically (to check)
///
///  See  example usage in the test macro AliTreePlayerTest.C
///  \code
///    AliTreePlayerTest.C::testAll();
///    AliTreePlayerTest.C::testSelectMetadata()
///    AliTreePlayerTest.C::testselectTreeInfo()
///    AliTreePlayerTest.C:testselectWhatWhereOrderByForTRD()
///  \endcode
///
/// \author Marian Ivanov


class TPad;
#include "TTreePlayer.h"
#include "TTreeFormula.h"

class AliTreeFormulaF : public TTreeFormula {
public:
  AliTreeFormulaF();
  ~AliTreeFormulaF();
  AliTreeFormulaF(const char *name, const char *formula, TTree *tree, Int_t debug=0);
  virtual Int_t Compile(const char *expression = "");
  virtual char *PrintValue(Int_t mode = 0) const; // { PrintValue(mode, 0, ""); }
  virtual char *PrintValue(Int_t mode, Int_t instance, const char *decform = "9.9") const;
  virtual void        UpdateFormulaLeaves();
  virtual Int_t       GetNdata(){return 1;} // TODO - support for vectors
public:
  TObjArray *fTextArray;       /// array of text inputs
  TObjArray *fFormatArray;     /// array of format strings to draw
  TObjArray *fFormulaArray;    /// array of TFormulas
  mutable TString    fValue;   /// current cache value of the formula
  Int_t      fDebug;           /// debug level
  ClassDef(AliTreeFormulaF,1)  ///
};


class AliTreePlayer : public TNamed {
public:
  AliTreePlayer() {}
  AliTreePlayer(const char *name, const char *title);
  static TObjArray  * selectMetadata(TTree * tree, TString query, Int_t verbose, TString *idList=NULL);
  static TObjArray  * selectTreeInfo(TTree* tree, TString query,Int_t verbose);
  static Int_t selectWhatWhereOrderBy(TTree * tree, TString what, TString where, TString orderBy,  Int_t firstentry, Int_t nentries, TString outputFormat, TString outputName);
  static TString  printSelectedTreeInfo(TTree*tree, TString infoType,  TString regExpFriend, TString regExpTag, Int_t verbose);
  static TObjArray  * MakeHistograms(TTree * tree, TString hisString, TString defaultCut, Int_t firstEntry, Int_t lastEntry, Int_t chunkSize=-1, Int_t verbose=1);
  static TPad *  DrawHistograms(TPad  * pad, TObjArray * hisArray, TString drawExpression, TObjArray *keepArray=0, Int_t verbose=0);
  static void MakeCacheTree(TTree * tree, TString varList, TString outFile, TString outTree, TCut selection,   Int_t nEntries=-1, Int_t firstEntry=0, const char *fileMode="recreate");
  static Int_t nextPad();

  template <typename T> static Long64_t BinarySearchSmaller(Long64_t n, const T *array, T value);
  enum TStatType {kUndef=-1,kEntries, kSum, kMean, kRMS, kMedian, kLTM, kLTMRMS, kMedianLeft,kMedianRight,kMax,kMin};
  static Int_t GetStatType(const TString &stat);
  static void AddStatInfo(TTree* treeLeft,  TTree * treeRight , const TString refQuery, Double_t deltaT,
		   const TString statString="median:medianLeft:medianRight:RMS:Mean:LTM0.60:LTMRMS0.60:Max:Min",
		   Int_t maxEntries=100000000);
  static TTree *  LoadTrees(const char *inputDataList, const char *  chRegExp, const char * chNotReg,  TString  inputFileSelection, TString axisAlias,  TString axisTitle, Int_t verbose=1);
  THashList *AddMetadata(TTree*, const char *vartagName,const char *varTagValue);
  TNamed *GetMetadata(TTree* tree, const char *vartagName, TString *prefix=0, Bool_t fullMatch=kFALSE);

  /// TODO -
  /// Metadata query from the TStatToolkit (GetMetadata, AddMetadata)
  /// TH1* DrawSorted(const char * expression, const char weights, Bool_t down, Int_t cut);  // draw sorted version of histogram
  /// Friend tree queries for time series using (time in entry or in array)
  ///       - calculate nearest entry
  ///       -- cache adding  branch
  ///       - calculate (index, mean, median, rms ...) entry in interval
  ///       -- cache adding branch
  /// To solve/consider:
  ///       - index branches have to have the same names - indeces to be created in all trees
  ///       - can we change indices? is it needed ?
  ///       - implementations for indices er entry and internal indices within array will be different
  ///       -- are both implementation needed - in which order ...
  ///static TH1* DrawSorted(TTree * tree, const char* varexp, const TCut& selection, Bool_t nbins=20, Long64_t nentries = 1000000000, Long64_t firstentry = 0);
  ///
protected:
  ClassDef(AliTreePlayer,1)  // Extension of the TTreePlayer
};


/// Binary search in an array of n values to locate value.
/// Array is supposed  to be sorted prior to this call.
// function gives nearest element smaller than value (different than in the TMath::BinarySearch
// See TMath::BinarySearch
/// \tparam T     - template (float,double,int ..)
/// \param n      - size of array
/// \param array  - pointer to raw array
/// \param value  - value to find
/// \return nearest element smaller than value.
template <typename T> Long64_t AliTreePlayer::BinarySearchSmaller(Long64_t n, const T *array, T value)
{
  const T* pind;
  pind = std::lower_bound(array, array + n, value);
  return ( pind - array - 1);
}

#endif

