#ifndef ALIMUONANALYSISOUTPUT_H
#define ALIMUONANALYSISOUTPUT_H

/// \class AliMuonAnalysisOutput
/// \brief Output handler for AliVAnalsisMuon
///
/// The class implements functionalities to handle more easily
/// the output of classes based on AliMuonAnalysisOutput
///
/// \author Diego Stocco <dstocco@cern.ch>, Subatech
/// \date Oct 21, 2015


#include "TObject.h"

class TObjArray;
class THashList;
class TAxis;
class AliCounterCollection;
class AliMergeableCollection;

class AliMuonAnalysisOutput : public TObject {
 public:
  AliMuonAnalysisOutput();
  AliMuonAnalysisOutput ( TObjArray* outputList );
  AliMuonAnalysisOutput ( const char *filename, const char *outputName );

  virtual ~AliMuonAnalysisOutput();

  TString GetCentralities ( const TString centralityRange ) const;
  
  // Methods for mergeable object collections
  TObject* GetMergeableObject ( TString physSel, TString trigClassName, TString centrality, TString objectName );
  TObject* GetSum ( TString physSel, TString trigClassNames, TString centrality, TString objectPattern );

  /// Get output list
  TObjArray* GetOutput () { return fOutputList; }
  /// Get counter collection
  AliCounterCollection* GetCounterCollection () { return fCounterCollection; }
  /// Get mergeable collection
  AliMergeableCollection* GetMergeableCollection () { return fMergeableCollection; }

 private:

  AliMuonAnalysisOutput(const AliMuonAnalysisOutput&);
  AliMuonAnalysisOutput& operator=(const AliMuonAnalysisOutput&);
  Bool_t Init();

  AliCounterCollection* fCounterCollection;  //!<! event counters
  AliMergeableCollection* fMergeableCollection; //!<! collection of mergeable objects
  TObjArray* fOutputList;  //!<! List of outputs

  /// \cond CLASSIMP
  ClassDef(AliMuonAnalysisOutput, 0); // Class for better handling of muon analysis output
  /// \endcond
};

#endif
