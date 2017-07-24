#ifndef ALICALOCALIBSIGNAL_H
#define ALICALOCALIBSIGNAL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//________________________________________________________________________
/// \class AliCaloCalibSignal
/// \ingroup EMCALbase
/// \brief class for signal monitoring and calibration tools
///
///  A help class for monitoring and calibration tools: MOOD, AMORE etc.,
///  that can process events from a standard AliCaloRawStreamV3,
///  most usually from LED/pulser runs. It stores signal info as
///  typical (highest) amplitude vs time in TGraphs (one per channel)
///  or TProfiles if we decide to just store the averages (and not all points) 
///  for the detectors (EMCAL and PHOS).
///
///  Partly based on AliCaloCalibPedestal.
///  It can be created and used a la (ctor):
///  * Create the object for making the histograms
///     * fSignals = new AliCaloCalibSignal( fDetType );
///     * AliCaloCalibSignal knows how many modules we have for PHOS or EMCAL
///     * fNumModules = fSignals->GetModules();
///  * fed an event:
///     * fSignals->ProcessEvent(fCaloRawStream,fRawEventHeaderBase);
///  * get some info:
///     * fSignals->GetXXX..()
/// etc.
///
/// \author: Josh Hamblen (UTenn), original version. 
/// \author D. Silvermyr (ORNL), consultant
//________________________________________________________________________

#include "TString.h"
#include "TTree.h"
#include "AliEMCALGeoParams.h"

class AliCaloRawStreamV3;
class AliCaloAltroMapping;
class AliRawReader;
class AliCaloRawAnalyzer;

class AliCaloCalibSignal : public TObject 
{
  
 public:

  /// The detector types (add common enum as in AliCalibPedestal?)
  enum kDetType {kPhos, kEmCal, kNone}; 
  
  AliCaloCalibSignal(kDetType detectorType = kEmCal); // ctor
  virtual ~AliCaloCalibSignal(); // dtor

private:
  
  // Just declare them, avoid compilation warning
  AliCaloCalibSignal             (const AliCaloCalibSignal & /*sig*/  ); // copy ctor
  AliCaloCalibSignal& operator = (const AliCaloCalibSignal &/*source*/); // assing operator
  
public:
  
  // Event processing methods:
  Bool_t   ProcessEvent(AliRawReader *rawReader);
  Bool_t   ProcessEvent(AliCaloRawStreamV3 *in, UInt_t Timestamp); // added header for time info
 
  Bool_t   CheckFractionAboveAmp(const int *AmpVal, int resultArray[]) const; 
  Bool_t   CheckLEDRefAboveAmp  (const int *AmpVal, int resultArray[]) const;

  // Mapping handling
  AliCaloAltroMapping **GetAltroMapping()                     const { return fMapping; }
  void                  SetAltroMapping(AliCaloAltroMapping **mapp) { fMapping = mapp; }

  // Fitter / Analyzer
  Int_t    GetFittingAlgorithm()        const { return fFittingAlgorithm ; }
  void     SetFittingAlgorithm(Int_t val) ;         
  AliCaloRawAnalyzer *GetRawAnalyzer()  const { return fRawAnalyzer      ; }  

  // Parameter/cut handling
  void     SetParametersFromFile(const char *parameterFile);
  void     WriteParametersToFile(const char *parameterFile);

  ////////////////////////////
  // Simple getters
  // * for TTree
  TTree *  GetTreeAmpVsTime()           const { return fTreeAmpVsTime      ; }
  TTree *  GetTreeAvgAmpVsTime()        const { return fTreeAvgAmpVsTime   ; }
  TTree *  GetTreeLEDAmpVsTime()        const { return fTreeLEDAmpVsTime   ; }
  TTree *  GetTreeLEDAvgAmpVsTime()     const { return fTreeLEDAvgAmpVsTime; }

  // how many points do we have for each tower&gain
  int      GetNHighGain(int imod, int icol, int irow) const 
    { int towId = GetTowerNum(imod, icol, irow); return fNHighGain[towId] ; }	
  int      GetNLowGain (int imod, int icol, int irow) const 
    { int towId = GetTowerNum(imod, icol, irow); return fNLowGain [towId] ; }

  int      GetNHighGain(int towId)      const { return fNHighGain[towId] ; }	
  int      GetNLowGain (int towId)      const { return fNLowGain [towId] ; }	

  // also for LED reference
  int      GetNRef(const int imod, const int istripMod, const int igain) const 
    { int refId = GetRefNum(imod, istripMod, igain); return fNRef[refId] ; }
  int      GetNRef(int refId)           const { return fNRef[refId] ; }

  // Basic info: getters  
  /// \return int, if this is a PHOS or EMCAL object
  kDetType GetDetectorType()            const { return fDetType    ; }
  /// \return string, if this is a PHOS or EMCAL object  
  TString  GetCaloString()              const { return fCaloString ; } 

  int      GetColumns()                 const {return fColumns ; } // The number of columns per module
  int      GetRows()                    const {return fRows    ; } // The number of rows per module
  int      GetLEDRefs()                 const {return fLEDRefs ; } // The number of LED references/monitors per module
  int      GetModules()                 const {return fModules ; } // The number of modules

  int      GetTowerNum(const int imod,  const int icol, const int irow) const 
  { return (imod*fColumns*fRows + icol*fRows + irow) ; } // help index

  /// \return Channel number with gain included
  int      GetChannelNum(const int imod, const int icol, const int irow, const int igain) const 
  { return (igain*fModules*fColumns*fRows + imod*fColumns*fRows + icol*fRows + irow) ; } 
  
  Bool_t   DecodeChannelNum(const int chanId, 
                            int *imod, int *icol, int *irow, int *igain) const; 
  
  // LED reference indexing
  
  /// \return LED channel number with gain included
  int      GetRefNum(const int imod, const int istripMod, const int igain) const 
  { return (igain*fModules*fLEDRefs + imod*fLEDRefs + istripMod) ; } 
  
  Bool_t   DecodeRefNum(const int refId, 
                        int *imod, int *istripMod, int *igain) const;
  
  // Basic Counters
  int      GetNEvents()                 const { return fNEvents         ; }
  int      GetNAcceptedEvents()         const { return fNAcceptedEvents ; }

  ///////////////////////////////
  // Get and Set Cuts
  
  //  * Section for if we should help with the event selection of what is likely LED events
  void     SetAmpCut(double d)                { fAmpCut = d    ; }
  double   GetAmpCut()                  const { return fAmpCut ; }
  
  void     SetReqFractionAboveAmpCutVal(double d) { fReqFractionAboveAmpCutVal = d    ; }
  double   GetReqFractionAboveAmpCutVal()   const { return fReqFractionAboveAmpCutVal ; }
  
  void     SetReqFractionAboveAmp(bool b)     { fReqFractionAboveAmp = b    ; }
  bool     GetReqFractionAboveAmp()     const { return fReqFractionAboveAmp ; }
  
  //  * also for LED Reference/Mon channels
  void     SetAmpCutLEDRef(double d)          { fAmpCutLEDRef = d    ; }
  double   GetAmpCutLEDRef()            const { return fAmpCutLEDRef ; }
 
  void     SetReqLEDRefAboveAmpCutVal(bool b) { fReqLEDRefAboveAmpCutVal = b    ; }
  bool     GetReqLEDRefAboveAmpCutVal() const { return fReqLEDRefAboveAmpCutVal ; }

  // *  We may select to get averaged info
  void     SetUseAverage(bool b)              { fUseAverage = b    ; }
  bool     GetUseAverage()              const { return fUseAverage ; }
  
  /// set length of the interval that should be used for the average calculation (determines number of bins in TProfile)
  void     SetSecInAverage(int secInAverage)  { fSecInAverage = secInAverage ; } 
  int      GetSecInAverage()            const { return fSecInAverage         ; }

  void     SetDownscale(int i)                { fDownscale = i    ; }
  int      GetDownscale()               const { return fDownscale ; }

  // * Info on time since start of run
  //   These times are typically the same, but not necessarily if the events do not come in order 
  //   Could be useful when we know the length of the run (i.e. after it is over), e.g. for PreProcessor
  double   GetHour()                    const { return fHour       ; } // time info for current event
  double   GetCurrentHour()             const { return fHour       ; } // time info for current event (same as GetHour(), just more explicitly named)
  double   GetLatestHour()              const { return fLatestHour ; } // the latest time encountered
  void     SetLatestHour(double d)            { fLatestHour = d    ; } 

  // * RunNumbers : setters and getters
  void     SetRunNumber(int runNo)            { fRunNumber = runNo ; }
  int      GetRunNumber()               const { return fRunNumber  ; }  
  
  // * Start-of-run timestamp : set and get
  void     SetStartTime(int startTime)        { fStartTime = startTime ; }
  int      GetStartTime()               const { return fStartTime      ; }

  /////////////////////////////
  // Analysis functions
  
  void     ResetInfo();
  Bool_t   AddInfo(const AliCaloCalibSignal *sig);  

  // Saving functions
  Bool_t   Save(TString fileName); 
  Bool_t   Analyze();

 private:
 
  void     DeleteTrees(); 
  void     Zero(); 
  void     CreateTrees(); 
  
 private:

  kDetType fDetType;                    ///<  The detector type for this object
  int      fColumns;                    ///<  The number of columns per module
  int      fRows;                       ///<  The number of rows per module
  int      fLEDRefs;	                  ///<  The number of LED references/monitors per module
  int      fModules;	                  ///<  The number of modules
 
  TString  fCaloString;                 ///<  ID for which detector type we have 
  AliCaloAltroMapping **fMapping;       //!<! Altro Mapping object
  
  Int_t    fFittingAlgorithm;           ///<  Select the fitting algorithm
  AliCaloRawAnalyzer *fRawAnalyzer;     //!<! e.g. for sample selection for fits
  
  int      fRunNumber;                  ///<  The run number. Needs to be set by the user.
  int      fStartTime;                  ///<  Time of first event

  double   fAmpCut;                     ///<  Amplitude cut value
  double   fReqFractionAboveAmpCutVal;  ///<  Required fraction that should be above cut
  bool     fReqFractionAboveAmp;        ///<  Flag to select if we should do some event selection based on amplitudes

  double   fAmpCutLEDRef;               ///<  Amplitude cut value for LED reference
  bool     fReqLEDRefAboveAmpCutVal;    ///<  Flag to select if we should require that signal is also seen in LED Reference/Monitoring channel

  double   fHour;                       ///<  Fraction of hour since beginning of run, for amp vs. time graphs, for current event
  double   fLatestHour;                 ///<  Largest fraction of hour since beginning of run, for amp vs. time graphs
  bool     fUseAverage;                 ///<  Flag to average graph points into over a time interval
  int      fSecInAverage;               ///<  Time interval for the graph averaging

  int      fDownscale;                  ///<  To select 1 out every N (fDownscale) events

  // status counters
  int      fNEvents;                    ///<  Number of events processed
  int      fNAcceptedEvents;            ///<  Number of events accepted

  // Constants needed by the class: EMCAL ones are kept in AliEMCALGeoParams.h
  static const int fgkPhosRows    = 64; ///<  Number of rows per module for PHOS
  static const int fgkPhosCols    = 56; ///<  Number of columns per module for PHOS
  static const int fgkPhosLEDRefs =  0; ///<  No LED monitor channels for PHOS
  static const int fgkPhosModules =  5; ///<  Number of modules for PHOS
  
  // From numbers above: EMCal+DCal has more possible towers than PHOS
  static const int fgkMaxTowers   = 23040; ///< AliEMCALGeoParams::fgkEMCALModules * AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;
  
  // for LED references; maximum from EMCAL
  static const int fgkMaxRefs     =   480; ///< AliEMCALGeoParams::fgkEMCALModules * AliEMCALGeoParams::fgkEMCALLEDRefs;

  static const int fgkNumSecInHr  =  3600; ///< Number of seconds in an hour, for the fractional hour conversion on the time graph
  
  // trees
  TTree  * fTreeAmpVsTime;              ///<  Store channel, gain, amp, and time info
  TTree  * fTreeAvgAmpVsTime;           ///<  Store channel, gain, amp, and time info, for averages
  TTree  * fTreeLEDAmpVsTime;           ///<  Store channel, gain, amp, and time info, for LED reference
  TTree  * fTreeLEDAvgAmpVsTime;        ///<  Store channel, gain, amp, and time info, for LED reference - averages

  // counters
  int      fNHighGain[fgkMaxTowers]  ;  ///<  Number of Amp. vs. Time readings per tower
  int      fNLowGain [fgkMaxTowers]  ;  ///<  Number of Amp. vs. Time readings per tower, for low gain
  int      fNRef     [fgkMaxRefs * 2];  ///<  Number of Amp. vs. Time readings per tower, for LED refs; *2 for both gains
  
  /// \cond CLASSIMP
  ClassDef(AliCaloCalibSignal, 9) ;
  /// \endcond
  
};
    
#endif // ALICALOCALIBSIGNAL_H
