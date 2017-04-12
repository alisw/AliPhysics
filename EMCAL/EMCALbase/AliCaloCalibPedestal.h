#ifndef ALICALOCALIBPEDESTAL_H
#define ALICALOCALIBPEDESTAL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//________________________________________________________________________
/// \class AliCaloCalibPedestal
/// \ingroup EMCALbase
/// \brief pedestal/bad map monitoring and calibration tools
///
/// A help class for monitoring and calibration tools: MOOD, AMORE etc.,
/// that can process events from a standard AliCaloRawStreamV3,
/// most usually from LED/pulser runs. It stores signal info as
/// typical (highest) amplitude, and pedestal info in geometrically-binned
/// 2D profiles of the detectors (EMCAL and PHOS).
/// Comparisons (ratios and differences) can be done with references.
//
/// It can be created and used a la (ctor):
///  * Create the object for making the histograms
///    * fPedestals = new AliCaloCalibPedestal( fDetType );
///    * AliCaloCalibPedestal knows how many modules we have for PHOS or EMCAL
///    * fNumModules = fPedestals->GetModules();
/// * fed an event:
///    * fPedestals->ProcessEvent(fCaloRawStream);
/// * asked to draw histograms:
///    * fPedestals->GetDeadMap(i)->Draw("col"); or
///    * fPedestals->GetPeakProfileHighGainRatio((i < fNumModules) ? i : fVisibleModule)->Draw("colz");
/// etc.
/// The pseudo-code examples above were from the first implementation in MOOD (summer 2007).
/// Partly based on AliTPCCalibPedestal.
///
/// \author Timo Alho (Jyvaskyla), original version. 
/// \author D. Silvermyr (ORNL) (consultant)
//________________________________________________________________________

#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2.h"
#include "TObjArray.h"

#include "AliEMCALGeoParams.h"
class AliCaloRawStreamV3;
class AliCaloAltroMapping;
class AliRawReader;

class AliCaloCalibPedestal : public TObject 
{
  
 public:

  /// The detector types
  enum kDetType {kPhos, kEmCal, kNone}; 
  
  /// The entries being put to the deadmap
  enum kDeadMapEntry{kAlive = 0, kDead, kHot, kWarning, kResurrected, kRecentlyDeceased, kNumDeadMapStates}; 
  
  AliCaloCalibPedestal(kDetType detectorType = kEmCal);
  virtual ~AliCaloCalibPedestal();

  // copy ctor, and '=' operator, are not fully tested/debugged yet
  // at least for now; the reference info is not copied from one to the other
  AliCaloCalibPedestal(AliCaloCalibPedestal &ped); 
  AliCaloCalibPedestal& operator = (AliCaloCalibPedestal &source);
  
  // Event processing methods:  
  Bool_t ProcessEvent(AliRawReader *rawReader);
  Bool_t ProcessEvent(AliCaloRawStreamV3  *in);
  
  // Mapping handling
  AliCaloAltroMapping **GetAltroMapping()     const { return fMapping ; }
  void  SetAltroMapping(AliCaloAltroMapping **mapp) { fMapping = mapp ; }

  // Parameter/cut handling
  void SetParametersFromFile(const char *parameterFile);
  void WriteParametersToFile(const char *parameterFile);

  ////////////////////////////
  // Simple getters
  
  // * Main profiles:
  
  /// \return a pointer to the low-gain pedestal profile
  TProfile2D * GetPedProfileLowGain (int i)            { ValidateProfiles() ; return (TProfile2D*)fPedestalLowGain [i] ; }	
  /// \return a pointer to the high-gain pedestal profile
  TProfile2D * GetPedProfileHighGain(int i)            { ValidateProfiles() ; return (TProfile2D*)fPedestalHighGain[i] ; }	
  
  /// \return a pointer to the low-gain LEDRef profile
  TProfile   * GetPedLEDRefProfileLowGain (int i)      { ValidateProfiles() ; return (TProfile*)fPedestalLEDRefLowGain [i] ; }	 
  /// \return a pointer to the high-gain LEDRef profile 
  TProfile   * GetPedLEDRefProfileHighGain(int i)      { ValidateProfiles() ; return (TProfile*)fPedestalLEDRefHighGain[i] ; }	
  
  /// \return a pointer to the low-gain peak-pedestal profile
  TProfile2D * GetPeakProfileLowGain (int i)           { ValidateProfiles() ; return (TProfile2D*)fPeakMinusPedLowGain [i] ; } 
  /// \return a pointer to the high-gain peak-pedestal profile
  TProfile2D * GetPeakProfileHighGain(int i)           { ValidateProfiles() ; return (TProfile2D*)fPeakMinusPedHighGain[i] ; } 	
  
  // * Differences to references:
  
  /// \return a pointer to the low-gain pedestal profile difference
  TProfile2D * GetPedProfileLowGainDiff (int i)        { ValidateComparisonProfiles() ; return (TProfile2D*)fPedestalLowGainDiff [i] ; } 	
   	/// \return a pointer to the high-gain pedestal profile difference
  TProfile2D * GetPedProfileHighGainDiff(int i)        { ValidateComparisonProfiles() ; return (TProfile2D*)fPedestalHighGainDiff[i] ; }
  
  /// \return a pointer to the low-gain LEDRef profile difference
  TProfile   * GetPedLEDRefProfileLowGainDiff (int i)  { ValidateComparisonProfiles() ; return (TProfile*)fPedestalLEDRefLowGainDiff [i] ; } 	
  /// \return a pointer to the high-gain LEDRef profile difference
  TProfile   * GetPedLEDRefProfileHighGainDiff(int i)  { ValidateComparisonProfiles() ; return (TProfile*)fPedestalLEDRefHighGainDiff[i] ; }
  
  /// \return a pointer to the low-gain peak-pedestal profile difference
  TProfile2D * GetPeakProfileLowGainDiff (int i)       { ValidateComparisonProfiles() ; return (TProfile2D*)fPeakMinusPedLowGainDiff [i] ; } 	
  /// \return a pointer to the high-gain peak-pedestal profile difference
  TProfile2D * GetPeakProfileHighGainDiff(int i)       { ValidateComparisonProfiles() ; return (TProfile2D*)fPeakMinusPedHighGainDiff[i] ; } 	
  
  // * Ratio to references:
  
  /// \return a pointer to the low-gain pedestal profile ratio
  TProfile2D * GetPedProfileLowGainRatio (int i)       { ValidateComparisonProfiles() ; return (TProfile2D*)fPedestalLowGainRatio [i] ; } 
  /// \return a pointer to the high-gain pedestal profile ratio
  TProfile2D * GetPedProfileHighGainRatio(int i)       { ValidateComparisonProfiles() ; return (TProfile2D*)fPedestalHighGainRatio[i] ; } 	
  
  /// \return a pointer to the low-gain LEDRef profile ratio
  TProfile   * GetPedLEDRefProfileLowGainRatio (int i) { ValidateComparisonProfiles() ; return (TProfile*)fPedestalLEDRefLowGainRatio [i] ; } 
  /// \return a pointer to the high-gain LEDRef profile ratio 
  TProfile   * GetPedLEDRefProfileHighGainRatio(int i) { ValidateComparisonProfiles() ; return (TProfile*)fPedestalLEDRefHighGainRatio[i] ; } 	
  
  /// \return a pointer to the low-gain peak-pedestal profile ratio
  TProfile2D * GetPeakProfileLowGainRatio (int i)      { ValidateComparisonProfiles() ; return (TProfile2D*)fPeakMinusPedLowGainRatio [i] ; } 
  /// \return a pointer to the high-gain peak-pedestal profile ratio
  TProfile2D * GetPeakProfileHighGainRatio(int i)      { ValidateComparisonProfiles() ; return (TProfile2D*)fPeakMinusPedHighGainRatio[i] ; } 	
  
  /// \return a pointer to the high-gain peak-pedestal histo
  TH2F       * GetPeakHighGainHisto(int i)             { ValidateProfiles() ; return (TH2F*)fPeakMinusPedHighGainHisto[i] ; } 	

  // * Bad channels:
  
  /// \return a pointer to the bad channels map histo
  TH2D       * GetDeadMap(int i)                       { ValidateProfiles() ; return (TH2D*)fDeadMap[i] ; }
//void         SetDeadMap(int i, TH2D *h)              { ((TH2D*)fDeadMap[i])=h ; }

  TObjArray    GetDeadMap()                            { ValidateProfiles() ; return fDeadMap ; }
  void         SetDeadMap(TObjArray map)               { fDeadMap = map ; }

  Bool_t IsBadChannel    (int imod, int icol, int irow) const; 
  void   SetChannelStatus(int imod, int icol, int irow, int status); 
  Int_t  GetChannelStatus(int imod, int icol, int irow) const { return  (Int_t)((TH2D*)fDeadMap[imod])->GetBinContent(icol, irow) ;	}
	
  ////////////////////////////
  // Basic info
  
  // * getters  
  
  /// \return if this is a PHOS or EMCAL object
  kDetType GetDetectorType()         const { return fDetType ; } 
  /// \return if this is a PHOS or EMCAL object
  TString  GetCaloString()           const { return fCaloString ; }  
  
  int      GetColumns()              const { return fColumns ; }  //The number of columns per module
  int      GetRows()                 const { return fRows    ; }  //The number of rows per module
  int      GetLEDRefs()              const { return fLEDRefs ; }  //The number of LED references/monitors per module
  int      GetModules()              const { return fModules ; }  //The number of modules
  int      GetRowMin()               const { return fRowMin  ; }  //for histo def.
  int      GetRowMax()               const { return fRowMax  ; }  //for histo def.
  int      GetRowMultiplier()        const { return fRowMultiplier ; }  //for histo filling

  // * RunNumbers : setters and getters
  void     SetRunNumber(int runNo)         { fRunNumber = runNo ; } 
  int      GetRunNumber()            const { return fRunNumber  ; } 
  int      GetRefRunNumber()         const { if (fReference) return fReference->GetRunNumber(); else return -1 ; } 

  // * Possibility to select only some samples for the pedestal calculation
  
  /// Set select to to use only some range of samples for pedestal calc.
  void     SetSelectPedestalSamples(Bool_t flag = kFALSE) { fSelectPedestalSamples = flag  ; } 
  /// Get select to to use only some range of samples for pedestal calc.
  Bool_t   GetSelectPedestalSamples()               const { return fSelectPedestalSamples  ; } 
  
  void     SetFirstPedestalSample(int i)   { fFirstPedestalSample = i    ; } // first sample to use
  void     SetLastPedestalSample (int i)   { fLastPedestalSample  = i    ; } // last sample to use
  int      GetFirstPedestalSample()  const { return fFirstPedestalSample ; } // first sample to use
  int      GetLastPedestalSample ()  const { return fLastPedestalSample  ; } // last sample to use

  // * Set threshold/event fraction for tower warnings
  void     SetDeadThreshold(int i)         { fDeadThreshold    = i    ; } // peak - pedestal dead threshold
  void     SetWarningThreshold(int i)      { fWarningThreshold = i    ; } // peak - pedestal warning threshold
  void     SetWarningFraction(double d)    { fWarningFraction  = d    ; } // event fraction for warnings
  int      GetDeadThreshold()        const { return fDeadThreshold    ; } // peak - pedestal dead threshold
  int      GetWarningThreshold()     const { return fWarningThreshold ; } // peak - pedestal warning threshold
  double   GetWarningFraction()      const { return fWarningFraction  ; } // event fraction for warnings
  
  // * hot towers
  void     SetHotSigma(double d)           { fHotSigma = d    ; } // rms away from normal
  double   GetHotSigma()             const { return fHotSigma ; }  // rms away from normal

  // * Basic counters
  int      GetNEvents()              const { return fNEvents    ; } 
  int      GetNChanFills()           const { return fNChanFills ; } 
  
  /////////////////////////////
  //Analysis functions
  
  /// \return the number of dead towers, by counting the bins in peak-pedestal smaller than threshold
  void     SetDeadTowerCount(Int_t dead)   { fDeadTowers = dead ; } 
   /// \return the number of dead towers, by counting the bins in peak-pedestal smaller than threshold
  int      GetDeadTowerCount()       const { return fDeadTowers ; }
  /// \return the percentage of dead towers, relative to a full module
  double   GetDeadTowerRatio()       const { return fDeadTowers/(double)(fRows*fColumns) ; } 
  /// \return the new dead towers compared to the reference
  int      GetDeadTowerNew()         const { return fNewDeadTowers ; } 
  /// \return the towers resurrected since the reference run
  int      GetDeadTowerResurrected() const { return fResurrectedTowers ; }  

  void     Reset();

  Bool_t   AddInfo(AliCaloCalibPedestal *ped);  
  
  //////////////////////////////////////////////////////
  //Functions related to comparing this with another (reference) run.
  Bool_t   LoadReferenceCalib(TString fileName, TString objectName); 
  
  ///Get the reference object. Needed for debug, will probably be removed later
  AliCaloCalibPedestal * GetReference() const { return fReference ; }  
  Bool_t                 SetReference(AliCaloCalibPedestal *ref);

  void     ComputeDiffAndRatio();
  void     ComputeDeadTowers(const char * deadMapFile = 0);
  void     ComputeHotAndWarningTowers(const char * hotMapFile = 0);

  // Saving functions
  Bool_t   SaveHistograms(TString fileName, Bool_t saveEmptyHistos = kFALSE); 

  void     Init() { ValidateProfiles() ; } // do basic setup

 private:
  
  void     ValidateProfiles(); 
  void     CompressAndSetOwner();
  void     ValidateComparisonProfiles();
  
  //////////////////////////////////////////////////////
  // The histograms. We use a TObjArray instead of a simple array,because this gives automatic streaming properties for the
  // class. A TClonesArray would be more efficient, but it's a bit more difficult to use and it doesn't matter too much
  // since we have only one object per module in the array anyway.
  
  TObjArray fPedestalLowGain;             ///<  Pedestal info for low gain
  TObjArray fPedestalHighGain;            ///<  Pedestal info for high gain
  TObjArray fPedestalLEDRefLowGain;       ///<  Pedestal LEDRef info for low gain
  TObjArray fPedestalLEDRefHighGain;      ///<  Pedestal LEDRef info for high gain
  TObjArray fPeakMinusPedLowGain;         ///<  (peak-pedestal) info for low gain
  TObjArray fPeakMinusPedHighGain;        ///<  (peak-pedestal) info for high gain

  TObjArray fPeakMinusPedHighGainHisto;   ///<  (peak-pedestal TH2F) info for high gain, used for hot towers eveluation
  
  // The difference of profiles between this and the reference object
  TObjArray fPedestalLowGainDiff;         //!<! Difference of profiles between pedestal low gain and the reference object
  TObjArray fPedestalHighGainDiff;        //!<! Difference of profiles between pedestal high gain and the reference object
  TObjArray fPedestalLEDRefLowGainDiff;   //!<! Difference of profiles between pedestal LED low gain and the reference object 
  TObjArray fPedestalLEDRefHighGainDiff;  //!<! Difference of profiles between pedestal LED high gain and the reference object  
  TObjArray fPeakMinusPedLowGainDiff;     //!<! Difference of profiles between peak low gain and the reference object 
  TObjArray fPeakMinusPedHighGainDiff;    //!<! Difference of profiles between peak high gain and the reference object 
  
  // The ratio of profiles between this and the reference object
  TObjArray fPedestalLowGainRatio;        //!<! The ratio of profiles between pedestal low gain and the reference object
  TObjArray fPedestalHighGainRatio;       //!<! The ratio of profiles between pedestal high gain and the reference object 
  TObjArray fPedestalLEDRefLowGainRatio;  //!<! The ratio of profiles between pedestal LED low gain and the reference object 
  TObjArray fPedestalLEDRefHighGainRatio; //!<! The ratio of profiles between pedestal LED high gain and the reference object  
  TObjArray fPeakMinusPedLowGainRatio;    //!<! The ratio of profiles between peak low gain and the reference object 
  TObjArray fPeakMinusPedHighGainRatio;   //!<! The ratio of profiles between peak high gain and the reference object 
  
  TObjArray fDeadMap;                     ///<  The dead map

  // Status counters
  int       fNEvents;                     ///<  Total number of events processed, 
  int       fNChanFills;                  ///<  Total number of channel fills (NChan * NEvents if not zero-suppressed)

  // The dead tower counts
  int       fDeadTowers;                  ///<  Number of towers found dead.
  int       fNewDeadTowers;               //!<! Towers that have died since the reference run
  int       fResurrectedTowers;           //!<! Towers that have been resurrected from the dead, compared to the reference
  
  AliCaloCalibPedestal * fReference;      //!<! A reference object, for comparing the accumulated results to a previous run
  
  kDetType  fDetType;                     ///<  The detector type for this object
  
  int       fColumns;	                    ///<  The number of columns per module
  int       fRows;	                      ///<  The number of rows per module
  
  int       fLEDRefs;	                    ///<  The number of LED references/monitors per module
  int       fModules;                     ///<  The number of modules
  
  int       fRowMin;                      ///<  Minimum Row number
  int       fRowMax;                      ///<  Maximum Row number
  int       fRowMultiplier;               ///<  Multiplication factor to get proper row range between PHOS and EMCAL
  
  TString   fCaloString;                  ///<  ID for which detector type we have 
  AliCaloAltroMapping **fMapping;         //!<! Altro Mapping object
  
  int       fRunNumber;                   ///<  The run number. Needs to be set by the user.
 
  Bool_t    fSelectPedestalSamples;       ///<  select to to use only some range of samples for pedestal calc.
  int       fFirstPedestalSample;         ///<  first sample to use
  int       fLastPedestalSample;          ///<  Last sample to use

  int       fDeadThreshold;               ///<  Peak - ped threshold used for dead towers evaluation
  int       fWarningThreshold;            ///<  Peak - ped threshold used for warm/warning towers evaluation
  double    fWarningFraction;             ///<  if(Peak - ped) > threshold in more than this fraction of event -> tower is assigned kWarning
  double    fHotSigma;                    ///<  if pedestal rms more than fHotSigma away from normal -> tower is assigned kHot

  // Constants needed by the class: EMCAL ones are kept in AliEMCALGeoParams.h
  static const int fgkPhosRows    = 64;   ///<  Number of rows per module for PHOS
  static const int fgkPhosCols    = 56;   ///<  Number of columns per module for PHOS
  static const int fgkPhosLEDRefs =  1;   ///<  No LED monitor channels for PHOS, set to 1 just to keep code simpler (also create LEDRef histos for PHOS)
  static const int fgkPhosModules =  5;   ///<  Number of modules for PHOS

  /// \cond CLASSIMP
  ClassDef(AliCaloCalibPedestal, 8) ;
  /// \endcond

};
    
#endif // ALICALOCALIBPEDESTAL_H
