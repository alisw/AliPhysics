#ifndef ALITRIGCHEFFOUTPUT_H
#define ALITRIGCHEFFOUTPUT_H

/// \class AliTrigChEffOutput
/// \brief Output for Trig chamber effieincy
///
/// The class manipulates the output of AliAnalysisTaskTrigChEff
/// in order to build the trigger chamber efficiency object
/// to be plugged in the OCDB for simulations
///
/// \author Diego Stocco <dstocco@cern.ch>, Subatech
/// \date Oct 21, 2015

#include "AliMuonAnalysisOutput.h"

class TString;
class TList;
class TObjArray;
class TH1;

class AliTrigChEffOutput : public AliMuonAnalysisOutput {
 public:
  AliTrigChEffOutput ( TObjArray* outputList, const char* name = "" );
  AliTrigChEffOutput ( const char *filename, const char *outputName = "testMTRChamberEff" );

  virtual ~AliTrigChEffOutput();

  TList* GetEffHistoList ( TString physSel, TString trigClassNames, TString centrality, Int_t itrackSel, Int_t imatch, Int_t imethod );
  TString GetHistoName ( Int_t itype, Int_t icount, Int_t ichamber, Int_t itrackSel, Int_t imatch, Int_t imethod );
  TH1* GetCountHisto ( Int_t itype, Int_t icount, Int_t ichamber, Int_t itrackSel, Int_t imatch, Int_t imethod );

  enum {
    kBendingEff,     ///< Bending plane fired
    kNonBendingEff,  ///< Non-bending plane fired
    kBothPlanesEff,  ///< Both planes fired
    kAllTracks,      ///< tracks used for calculation
    kNcounts         ///< Number of count type
  };

  enum {
    kHchamberEff,    ///< Counts per cathode histogram index
    kHslatEff,       ///< Counts per slat histogram index
    kHboardEff,      ///< Counts per board histogram index
    kHcheckBoard,    ///< Check rejected tracks per board
    kNhistoTypes     ///< Check rejected tracks per board
  };

  enum {
    kNoMatch,       ///< No match with trigger
    kMatchApt,      ///< Match All Pt
    kMatchLpt,      ///< Match Low Pt
    kMatchHpt,      ///< Match High Pt
    kNtrigMatch     ///< Total number of matched types
  };

  enum {
    kSelectTrack,   ///< Selected track
    kNoSelectTrack, ///< Non selected tracks (includes ghosts)
    kNtrackSel      ///< Total number of track selection
  };

  enum {
    kEffFromTrack,  ///< Hit pattern from tracker track extrapolation
    kEffFromTrig,   ///< Hit pattern from trigger
    kNeffMethods    ///< Total number of efficiency methods
  };

 private:

  /// Dummy
  AliTrigChEffOutput(const AliTrigChEffOutput&);
  /// Dummy
  AliTrigChEffOutput& operator=(const AliTrigChEffOutput&);

  void InitKeys();
 
  TObjArray* fTrackSelKeys;   ///!<! Selection names
  TObjArray* fCountTypeKeys;  ///!<! Count type keys
  TObjArray* fHistoTypeKeys;  ///!<! Base histogram name
  TObjArray* fEffMethodKeys;  ///!<! Efficiency methods keys
  TObjArray* fMatchTrigKeys;  ///!<! Match trigger names

  /// \cond CLASSIMP
  ClassDef(AliTrigChEffOutput, 0); // Trigger chamber efficiencies
  /// \endcond
};

#endif
