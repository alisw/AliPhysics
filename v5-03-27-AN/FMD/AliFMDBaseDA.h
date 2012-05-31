// -*- mode: C++ -*- 
#ifndef ALIFMDBASEDA_H
#define ALIFMDBASEDA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
//
// This class provides a base interface for the Detector Algorithms
// (DA) of the FMD.  At least three implementations are needed:
// AliFMDPedestalDA, AliFMDGainDA and AliFMDPhysicsDA .  These classes
// will provide the calibration data for the AliFMDPreprocessor to be
// used in the shuttle.  The input for this class are raw data
// (AliRawReader) and the output is a comma-separated file
// (std::ofstream) that contains the values defined in the
// implementations of this class.
//
// Author: Hans Hjersing Dalsgaard, hans.dalsgaard@cern.ch
//

#include "TNamed.h"
#include "TObjArray.h"
#include "TString.h"
#include "TArrayS.h"
#include <iosfwd>
#include <fstream>
class AliFMDDigit;
class AliRawReader;
class AliFMDParameters;
class AliFMDRawReader;
class TDirectory;
class TH2;
class TFile;
class TClonesArray;

class AliFMDBaseDA: public TNamed 
{
public:
  /** 
   * Constructor 
   * 
   */
  AliFMDBaseDA() ;
  /** 
   * Copy constructor 
   * 
   * @param baseDA 
   */
  AliFMDBaseDA(const AliFMDBaseDA & baseDA) ;
  //  AliFMDBaseDA& operator = (const AliFMDBaseDA & baseDA) ; 
  /** 
   * Destructor
   * 
   */  
  ~AliFMDBaseDA() ;
  AliFMDBaseDA& operator=(const AliFMDBaseDA&) { return *this; }
  /** 
   * Run this DA
   * 
   * @param fmdReader Raw input reader
   */  
  void Run(AliRawReader* fmdReader);
  /** 
   * Set whether to save diagnostics 
   * 
   * @param save If true, will output diagnostics file
   */
  void SetSaveDiagnostics(Bool_t save) {fSaveHistograms = save;}
  /** 
   * Set whether to make summary histograms to be published to AMORE
   * 
   * @param save If true, will generate summary QA histograms
   */
  void SetMakeSummaries(Bool_t save) {fMakeSummaries = save;}
  /** 
   * Set the number of requried events
   * 
   * @param nEvents Number of event we need
   */
  void SetRequiredEvents(Int_t nEvents) {fRequiredEvents = nEvents;}
  /** 
   * Get the number of required events
   * 
   * 
   * @return number of required events
   */
  Int_t GetRequiredEvents() const {return fRequiredEvents ;}
  /** 
   * Get list of summary histograms 
   *
   * @return Array of summary histograms or null if not defined 
   */
  const TObjArray& GetSummaries() const { return fSummaries; }
  /** 
   * Check if we saw data for detector 
   * 
   * @param det Detector number to check 
   * @return true if the code has seen data from the detector 
   */
  Bool_t HasSeenDetector(UShort_t d) const;
protected:
  /** 
   * Initialize 
   */  
  virtual void Init()  {};
  /** 
   * Fill channels 
   */
  virtual void FillChannels(AliFMDDigit* )  {};
  /** 
   * Analyse a single strip result
   */
  virtual void Analyse(UShort_t, Char_t, UShort_t, UShort_t)  {};
  /** 
   * Write header to output file
   */
  virtual void WriteHeaderToFile()  {};
  /** 
   * Add a strip container 
   */
  virtual void AddChannelContainer(TObjArray*, UShort_t, Char_t, 
				   UShort_t, UShort_t )  {};
  /** 
   * End of event
   */
  virtual void FinishEvent()  {};
  /** 
   * End of run
   */
  virtual void Terminate(TFile* ) {};
  /** 
   * Current event number
   * 
   * 
   * @return 
   */  
  Int_t GetCurrentEvent() const {return fCurrentEvent;}
  /** 
   * Rotate a set of files.   @a base is the basic name of the files.
   * If the file @a base.max exists it is removed. 
   * If the file @a base.n exists (where n < max) it is renamed to @a
   * base.(n-1).  
   * If the file @a base exists, it is renamed to @a base.1 
   * 
   * @param base Base name of the files
   * @param max  Maximum number to keep (minus one for the current).
   */
  void Rotate(const char* base, int max) const;
  /** 
   * Ge the half-ring index 
   * 
   * @param UShort_t 
   * @param Char_t 
   * @param UShort_t 
   * 
   * @return 
   */
  Int_t GetHalfringIndex(UShort_t, Char_t, UShort_t) const;
  /** 
   * Get the pulse size 
   * 
   * @param det   Detector number
   * @param ring  Rin identifier
   * @param board Board number 
   * 
   * @return Pulse step size
   */
  Int_t GetPulseSize(UShort_t det , 
		     Char_t ring, 
		     UShort_t board) 
  {
    return fPulseSize.At(GetHalfringIndex(det,ring,board));
  }
  /** 
   * Get number of events per pulse size 
   * 
   * @param det   Detector number
   * @param ring  Rin identifier
   * @param board Board number 
   * 
   * @return number of events per Pulse size
   */
  Int_t GetPulseLength(UShort_t det, 
		       Char_t ring, 
		       UShort_t board) 
  {
    return fPulseLength.At(GetHalfringIndex(det,ring,board));
  }

  /** 
   * Get the detector path in diagnositcs file
   * 
   * @param det  Detector number
   * @param full If true, return full path
   * 
   * @return Path to detector
   */  
  const char* GetDetectorPath(UShort_t det, Bool_t full=kTRUE) const;
  /** 
   * Get the ring path in diagnositcs file
   * 
   * @param det  Detector number
   * @param ring Ring identifier 
   * @param full If true, return full path
   * 
   * @return Path to ring
   */  
  const char* GetRingPath(UShort_t det, Char_t ring, Bool_t full=kTRUE) const;
  /** 
   * Get the sector path in diagnositcs file
   * 
   * @param det  Detector number
   * @param ring Ring identifier 
   * @param sec  Sector number
   * @param full If true, return full path
   * 
   * @return Path to sector
   */  
  const char* GetSectorPath(UShort_t det, Char_t ring, UShort_t sec, 
			    Bool_t full=kTRUE) const;
  /** 
   * Get the strip path in diagnositcs file
   * 
   * @param det  Detector number
   * @param ring Ring identifier 
   * @param sec  Sector number
   * @param str  Strip number
   * @param full If true, return full path
   * 
   * @return Path to strip
   */  
  const char* GetStripPath(UShort_t det, Char_t ring, UShort_t sec, 
			   UShort_t str, Bool_t full=kTRUE) const;
  /** 
   * Write conditions file 
   * 
   * @param fmdReader Raw input
   */ 
  void WriteConditionsData(AliFMDRawReader* fmdReader);
  /** 
   * Set the current event 
   * 
   * @param currentEvent 
   */
  void SetCurrentEvent(Int_t currentEvent) {fCurrentEvent = currentEvent; }
  /** 
   * Initialize container 
   * 
   * @param dir Directory to make containers in 
   */
  void InitContainer(TDirectory* dir);
  /** 
   * Utility function for defining summary histograms 
   *
   * @param det    Detector 
   * @param ring   Ring identifier 
   * @param prefix Histogram prefix 
   * @param title  Histogram title 
   */
  TH2* MakeSummaryHistogram(const char* prefix, const char* title, 
			    UShort_t det, Char_t ring);
  /** 
   * Make a summary
   * 
   */
  virtual void  MakeSummary(UShort_t, Char_t) { }



  static const UInt_t fgkBaseDDL = 3072;   // base FMD ddl
  //Char_t* fDiagnosticsFilename;
  TString       fDiagnosticsFilename;  // name of diagnostics file
  std::ofstream fOutputFile;           // output file
  std::ofstream fConditionsFile;       // conditions file
  Bool_t        fSaveHistograms;       // save hists or not
  Bool_t        fMakeSummaries;        // save hists or not
  TObjArray     fDetectorArray;        // array indiced by detector
  TArrayS       fPulseSize;            // Pulse size for gain calib
  TArrayS       fPulseLength;          // Pulse length for gain calib
  Bool_t        fSeenDetectors[3];     // Detectors seen so far
  Int_t         fRequiredEvents;       // # events required for this calib
  Int_t         fCurrentEvent;         // the current event       
  UInt_t        fRunno;                // Current run number 
  TObjArray     fSummaries;            // Summary histograms 
  
  ClassDef(AliFMDBaseDA,0) // Base Detector algorithm for all run types

};
//____________________________________________________________________
inline Bool_t
AliFMDBaseDA::HasSeenDetector(UShort_t d) const
{ 
  return (d == 0 || d > 3) ? false : fSeenDetectors[d-1]; 
}

#endif

