/*************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 *************************************************************************
 * $Id$ */
/**
 * @file    AliFMDParameters.cxx
 * @author  Christian Holm Christensen <cholm@nbi.dk>
 * @date    Mon Mar 27 12:44:26 2006
 * @brief   Manager of FMD parameters     
 */
//____________________________________________________________________
//                                                                          
// Forward Multiplicity Detector based on Silicon wafers. 
//
// This class is a singleton that handles various parameters of
// the FMD detectors.  
// The manager normally serves the parameters from the Conditions
// Database (CDB).  These are retrivied by the member function
// `Init'.  Optionally, the class can serve hard-coded constants, if
// no CDB is available. 
//                                                       
#include "AliFMDDebug.h"		   // ALILOG_H
#include "AliFMDParameters.h"	   // ALIFMDPARAMETERS_H
#include "AliFMDGeometry.h"	   // ALIFMDGEOMETRY_H
#include "AliFMDRing.h"	           // ALIFMDRING_H
#include "AliFMDCalibGain.h"       // ALIFMDCALIBGAIN_H
#include "AliFMDCalibPedestal.h"   // ALIFMDCALIBPEDESTAL_H
#include "AliFMDCalibSampleRate.h" // ALIFMDCALIBSAMPLERATE_H
#include "AliFMDCalibStripRange.h" // ALIFMDCALIBSTRIPRANGE_H
#include "AliFMDAltroMapping.h"    // ALIFMDALTROMAPPING_H
#include <AliCDBManager.h>         // ALICDBMANAGER_H
#include <AliCDBEntry.h>           // ALICDBMANAGER_H
#include <AliFMDPreprocessor.h>
#include <AliLog.h>
#include <Riostream.h>
#include <sstream>
#include <TSystem.h>
#include <TArrayF.h>
#include <TH2D.h>

//====================================================================
ClassImp(AliFMDParameters)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDParameters* AliFMDParameters::fgInstance = 0;

//____________________________________________________________________
const char* AliFMDParameters::fgkPulseGain	    = "FMD/Calib/PulseGain";
const char* AliFMDParameters::fgkPedestal	    = "FMD/Calib/Pedestal";
const char* AliFMDParameters::fgkDead	            = "FMD/Calib/Dead";
const char* AliFMDParameters::fgkSampleRate	    = "FMD/Calib/SampleRate";
const char* AliFMDParameters::fgkAltroMap	    = "FMD/Calib/AltroMap";
const char* AliFMDParameters::fgkZeroSuppression    = "FMD/Calib/ZeroSuppression";
const char* AliFMDParameters::fgkStripRange	    = "FMD/Calib/StripRange";
const char* AliFMDParameters::fkPedestalShuttleID   = "pedestals";
const char* AliFMDParameters::fkGainShuttleID       = "gains";
const char* AliFMDParameters::fkConditionsShuttleID = "conditions";

//____________________________________________________________________
AliFMDParameters* 
AliFMDParameters::Instance() 
{
  // 
  // Get static instance 
  //
  if (!fgInstance) fgInstance = new AliFMDParameters;
  return fgInstance;
}

//____________________________________________________________________
AliFMDParameters::AliFMDParameters() 
  : fIsInit(kFALSE),
    fkSiDeDxMip(1.664), 
    fVA1MipRange(0),
    fAltroChannelSize(0),
    fChannelsPerAltro(0),
    fPedestalFactor(0),
    fZSPre(1),
    fZSPost(1),
    fZSPedSubtract(kTRUE),
    fFixedPedestal(100),
    fFixedPedestalWidth(2),
    fFixedZeroSuppression(1),
    fFixedSampleRate(2),
    fFixedThreshold(0),
    fFixedMinStrip(0),
    fFixedMaxStrip(127),
    fFixedPulseGain(2), 
    fEdepMip(0),
    fHasCompleteHeader(kTRUE),
    fZeroSuppression(0), 
    fSampleRate(0), 
    fPedestal(0), 
    fPulseGain(0), 
    fDeadMap(0), 
    fAltroMap(0), 
    fStripRange(0)
{
  //
  // Default constructor 
  //
  SetVA1MipRange();
  SetAltroChannelSize();
  SetChannelsPerAltro();
  SetZeroSuppression();
  SetSampleRate();
  SetPedestal();
  SetPedestalWidth();
  SetPedestalFactor();
  SetThreshold();
  SetStripRange();
  SetGain();
  fAltroMap = new AliFMDAltroMapping;
}

//__________________________________________________________________
void
AliFMDParameters::Init(Bool_t forceReInit, UInt_t what)
{
  // 
  // Initialize the manager.  This tries to read the parameters from
  // CDB.  If that fails, the class uses the hard-coded parameters.
  // 
  // Parameters:
  //    forceReInit Force (re-)initalize flag
  //    what        What to initialize 
  //
  if (forceReInit) fIsInit = kFALSE;
  if (fIsInit) return;
  if (what & kPulseGain)       InitPulseGain();
  if (what & kPedestal)        InitPedestal();
  if (what & kDeadMap)         InitDeadMap();
  if (what & kSampleRate)      InitSampleRate();
  if (what & kZeroSuppression) InitZeroSuppression();
  if (what & kAltroMap)        InitAltroMap();
  if (what & kStripRange)      InitStripRange();
  fIsInit = kTRUE;
}
//__________________________________________________________________
void
AliFMDParameters::Init(AliFMDPreprocessor* pp, Bool_t forceReInit, UInt_t what)
{
  // 
  // Initialize the manager.  This tries to read the parameters from
  // CDB.  If that fails, the class uses the hard-coded parameters.
  // 
  // Parameters:
  //    pp          Preprocessor 
  //    forceReInit Force (re-)initalize flag
  //    what        What to initialize 
  //
  if (forceReInit) fIsInit = kFALSE;
  if (fIsInit) return;
  if (what & kPulseGain)       InitPulseGain(pp);
  if (what & kPedestal)        InitPedestal(pp);
  if (what & kDeadMap)         InitDeadMap(pp);
  if (what & kSampleRate)      InitSampleRate(pp);
  if (what & kZeroSuppression) InitZeroSuppression(pp);
  if (what & kAltroMap)        InitAltroMap(pp);
  if (what & kStripRange)      InitStripRange(pp);
  fIsInit = kTRUE;
}

//__________________________________________________________________
Bool_t
AliFMDParameters::CheckFile(const char* prefix, 
			    const char* path, 
			    int         number, 
			    TString&    f) const
{
  // 
  // Check if the file <i>prefix</i><i>number</i> exists in @a path, 
  // and write the full path to @a f.  
  // 
  // Parameters:
  //    prefix  File prefix (cond, peds, gains, ...)
  //    path    Path to files
  //    number  Detector number (1, 2, or 3)
  //    f       On return full path to file (if found)
  // 
  // Return:
  //    @c true if file exists and is readable, @c false otherwise
  //
  f = (Form("%s%d.csv", prefix, number));
  AliFMDDebug(5, ("Checking if %s exists in %s ...", f.Data(), path));
  f = gSystem->Which(path, f.Data());
  AliFMDDebug(5, ("Got back '%s'", f.Data()));
  return !f.IsNull();
}

//__________________________________________________________________
void
AliFMDParameters::Init(const char* path, Bool_t forceReInit, UInt_t what)
{
  // 
  // Initialize the manager.  This will try to read some calibrations
  // (sample rate, strip range, gains, pedestals) from local comma
  // separated value (CSV) files in the directory pointed at by @a
  // path.  If they are not found, then they will be retrieved from
  // OCDB as appropriately.   Other calibrations are always read from
  // OCDB.  
  // 
  // The CSV files should be named as 
  // 
  // - Pedestals: <tt>peds</tt><i>det_number</i><tt>.csv</tt>
  // - Gains: <tt>gains</tt><i>det_number</i><tt>.csv</tt>
  // - Sample Rate: <tt>conditions</tt><i>det_number</i><tt>.csv</tt>
  // - Strip Range: <tt>conditions</tt><i>det_number</i><tt>.csv</tt>
  //
  // where <i>det_number</i> is the detector number (1, 2, or 3). 
  // 
  // Parameters:
  //    path        Where to look for the CSV files
  //    forceReInit Always reinitialise 
  //    what        What calibrations to load. 
  //  
  if (forceReInit) fIsInit = kFALSE;
  if (fIsInit) return;

  AliFMDCalibStripRange*  range = 0;
  AliFMDCalibSampleRate*  rate  = 0;
  AliFMDCalibPedestal*    peds  = 0;
  AliFMDCalibGain*        gains = 0;

  for (Int_t i = 1; i <= 3; i++) { 
    TString f;
    if (((what & kSampleRate) || (what & kStripRange)) && 
	CheckFile("conditions", path, i, f)) {
      if (!rate  && (what & kSampleRate)) rate  = new AliFMDCalibSampleRate;
      if (!range && (what & kStripRange)) range = new AliFMDCalibStripRange;
      std::ifstream in(f.Data());
      if (range) range->ReadFromFile(in);
      if (rate)  rate->ReadFromFile(in);
      in.close();
    }
    if ((what & kPedestal) && CheckFile("peds", path, i, f)) {
      if (!peds) peds  = new AliFMDCalibPedestal;
      std::ifstream in(f.Data());
      peds->ReadFromFile(in);
      in.close();
    }
    if ((what & kPulseGain) && CheckFile("gains", path, i, f)) { 
      if (!gains) gains = new AliFMDCalibGain;
      std::ifstream in(f.Data());
      gains->ReadFromFile(in);
      in.close();
    }
  }

  if (range) what &= ~kStripRange;
  if (rate)  what &= ~kSampleRate;
  if (peds)  what &= ~kPedestal;
  if (gains) what &= ~kPulseGain;

  Init(kFALSE, what);
  
  if (range) SetStripRange(range);
  if (rate)  SetSampleRate(rate);
  if (peds)  SetPedestal(peds);
  if (gains) SetGain(gains);

  fIsInit = kTRUE;
}

//__________________________________________________________________
void
AliFMDParameters::MakeDeadMap(Float_t maxNoise, 
			      Float_t minGain, 
			      Float_t maxGain)
{
  // 
  // Automatically generate a dead map from the pedestals and gains.
  // A channel is marked as dead of the noise is too high (currently
  // more than 10 ADC counts), or the gain is unreasonable (currently
  // larger than 10, or smaller than 0.1). 
  // 
  // The procedure does not overwrite channels previously marked as
  // dead - e.g., channels marked as dead in the calibration loaded
  // from OCDB will continue to be marked as dead.  That is, this
  // procedure will never make a channel un-dead. 
  // 
  // Parameters:
  //    maxNoise  Maximum noise value before a channel is marked
  // as dead. 
  //    minGain   Minimum value of the calibrated gain before a
  // channel is considered dead. 
  //    maxGain   Maximum value of the calibrated gain before a
  // channel is considered dead. 
  //
  if (fPedestal)  
    fDeadMap = fPedestal->MakeDeadMap(maxNoise, fDeadMap);
  if (fPulseGain) 
    fDeadMap = fPulseGain->MakeDeadMap(minGain, maxGain, fDeadMap);
}
//__________________________________________________________________
#define DET2IDX(det,ring,sec,str) \
  (det * 1000 + (ring == 'I' ? 0 : 512) + str)  
  
//__________________________________________________________________
void
AliFMDParameters::Draw(Option_t* option)
{
  // 
  // Draw parameters. 
  // 
  // Parameters:
  //    option What to draw. Should be one of 
  // - dead	  Dead channels
  // - threshold Threshold
  // - gain	  Gain
  // - pedestal  Pedestal
  // - noise	  Noise (or pedestal width)
  // - zero	  Zero suppression
  // - rate	  Sampling rate (VA1 clock / ALTRO clock)
  // - min	  Minimum strip read out
  // - max 	  Maximum strip read out
  // - map	  hardware address
  //
  TString opt(option);
  enum {
    kLocalPulseGain,       // Path to PulseGain calib object
    kLocalThreshold,       // Path to PulseGain calib object
    kLocalPedestal,        // Path to Pedestal calib object
    kLocalPedestalWidth,   // Path to Pedestal calib object
    kLocalDead,            // Path to Dead calib object
    kLocalSampleRate,      // Path to SampleRate calib object
    kLocalAltroMap,        // Path to AltroMap calib object
    kLocalZeroSuppression, // Path to ZeroSuppression cal object
    kLocalMinStripRange,   // Path to strip range cal object
    kLocalMaxStripRange    // Path to strip range cal object
  } what;
    
  if      (opt.Contains("dead", TString::kIgnoreCase)) 
    what = kLocalDead;
  else if (opt.Contains("threshold",TString::kIgnoreCase)) 
    what = kLocalThreshold;
  else if (opt.Contains("gain",TString::kIgnoreCase)) 
    what = kLocalPulseGain;
  else if (opt.Contains("pedestal",TString::kIgnoreCase)) 
    what = kLocalPedestal;
  else if (opt.Contains("noise",TString::kIgnoreCase)) 
    what = kLocalPedestalWidth;
  else if (opt.Contains("zero",TString::kIgnoreCase)) 
    what = kLocalZeroSuppression;
  else if (opt.Contains("rate",TString::kIgnoreCase)) 
    what = kLocalSampleRate;
  else if (opt.Contains("min",TString::kIgnoreCase)) 
    what = kLocalMinStripRange;
  else if (opt.Contains("max",TString::kIgnoreCase)) 
    what = kLocalMaxStripRange;
  else if (opt.Contains("map",TString::kIgnoreCase)) 
    what = kLocalAltroMap;
  else {
    Warning("Draw", "unknown parameter: %s\n\tShould be one of\n\t"
	    "dead, threshold, gain, pedestal, noise, zero, rate, "
	    "min, max, map",  
	    option); 
    return;
  }

  TArrayD xbins(3 * 512 + 2 * 256 + 5);
  Int_t i = 1;
  Bool_t skip = kTRUE;
  for (UShort_t det = 1; det <= 3; det++) {
    UShort_t nRings = (det == 1 ? 1 : 2);
    for (UShort_t iring = 0; iring < nRings; iring++) {
      UShort_t nStrip  = (iring == 0 ? 512 : 256);
      Char_t   ring    = (iring == 0 ? 'I' : 'O');
      for (UShort_t str = 0; str < nStrip; str++) {
	// UShort_t nSec    = (iring == 0 ? 20  : 40);
	// Char_t   ring    = (iring == 0 ? 'I' : 'O');
	// for (UShort_t sec = 0; sec < nSec; sec++) {
	Int_t idx = DET2IDX(det, ring, 0, str);
	// Int_t idx = DET2IDX(det, ring, sec, 0);
	if (skip) {
	  xbins[i-1] = idx - .5;
	  skip  = kFALSE;
	}
	xbins[i] = idx + .5;
	i++;
      }
      skip = kTRUE;
      i++;
    }
  }
  TArrayD ybins(41);
  for (/*Int_t*/ i = 0; i < ybins.fN; i++) ybins[i] = Float_t(i - .5);
  TH2D* hist = new TH2D("calib", Form("Calibration %s", option), 
			xbins.fN-1, xbins.fArray,  
			ybins.fN-1, ybins.fArray);
  hist->GetXaxis()->SetTitle("1000 #times detector + 512 #times ring + strip");
  hist->GetYaxis()->SetTitle("sector");
  
  // hist->Draw("Lego");
  // return;
  
  for (UShort_t det = 1; det <= 3; det++) {
    UShort_t nRings = (det == 1 ? 1 : 2);
    for (UShort_t iring = 0; iring < nRings; iring++) {
      UShort_t nSector = (iring == 0 ?  20 : 40);
      UShort_t nStrip  = (iring == 0 ? 512 : 256);
      Char_t   ring    = (iring == 0 ? 'I' : 'O');
      for (UShort_t sec = 0; sec < nSector; sec++) {
	for (UShort_t str = 0; str < nStrip; str++) {
	  Int_t idx = DET2IDX(det, ring, sec, str);
	  UShort_t ddl, addr, time, sam=0;
	  Double_t val = 0;
	  switch (what) {
	  case kLocalPulseGain:       // Path to PulseGain calib object
            val = GetPulseGain(det,ring,sec,str); break;
	  case kLocalThreshold:       // Path to PulseGain calib object
            val = GetThreshold(); break;
	  case kLocalPedestal:        // Path to Pedestal calib object
            val = GetPedestal(det,ring,sec,str); break;
	  case kLocalPedestalWidth:   // Path to Pedestal calib object
            val = GetPedestalWidth(det,ring,sec,str); break;
	  case kLocalDead:            // Path to Dead calib object
            val = IsDead(det,ring,sec,str); break;
	  case kLocalSampleRate:      // Path to SampleRate calib object
            val = GetSampleRate(det,ring,sec,str); break;
	  case kLocalAltroMap:        // Path to AltroMap calib object
	    Detector2Hardware(det,ring,sec,str,sam,ddl,addr,time); 
            val = addr; break;
	  case kLocalZeroSuppression: // Path to ZeroSuppression cal object
            val = GetZeroSuppression(det,ring,sec,str); break;
	  case kLocalMinStripRange:   // Path to strip range cal object
            val = GetMinStrip(det,ring,sec,str); break;
	  case kLocalMaxStripRange:    // Path to strip range cal object
            val = GetMaxStrip(det,ring,sec,str); break;
	  }
	  hist->Fill(idx,sec,val);
	  // hist->Fill(idx,str,val);
	}
      }
    }
  }
  hist->Draw("lego");
}

//__________________________________________________________________
void
AliFMDParameters::Print(Option_t* option) const
{
  // Print information. 
  // If option contains an 'A' then everything is printed. 
  // If the option contains the string "FMD" the function will search 
  // for detector, ring, sector, and strip numbers to print, in the
  // format 
  // 
  //    FMD<detector><ring>[<sector>,<string>] 
  // 
  // The wild card '*' means all of <detector>, <ring>, <sector>, or 
  // <strip>. 
  TString opt(option);
  Bool_t showStrips  = opt.Contains("a", TString::kIgnoreCase);
  UShort_t ds[]      = { 1, 2, 3, 0 };
  Char_t   rs[]      = { 'I', 'O', '\0' };
  UShort_t minStrip  = 0;
  UShort_t maxStrip  = 512;
  UShort_t minSector = 0;
  UShort_t maxSector = 40;
  
  
  if (opt.Contains("fmd",TString::kIgnoreCase)) {
    Int_t   i    = opt.Index("fmd",TString::kIgnoreCase);
    Int_t   j    = opt.Index("]",TString::kIgnoreCase);
    if (j != kNPOS)
      showStrips    = kTRUE;
    else 
      j = opt.Length();
    enum {
      kReadDet, 
      kReadRing, 
      kReadLbrack,
      kReadSector,
      kReadComma,
      kReadStrip,
      kReadRbrack, 
      kEnd
    } state = kReadDet;
    std::stringstream s(opt(i+4, j-i-3).Data());
    while (state != kEnd) {
      Char_t tmp = s.peek();
      if (tmp == ' ' || tmp == '\t') {
	s.get();
	continue;
      }
      switch (state) {
      case kReadDet: { // First, try to kRead the detector 
	if (tmp == '*') s.get();
	else { 
	  UShort_t det;
	  s >> det;
	  if (!s.bad()) {
	    ds[0] = det;
	    ds[1] = 0;
	  }
	}
	state = (s.bad() ? kEnd : kReadRing);
      } break;
      case kReadRing: { // Then try to read the ring;
	Char_t ring;
	s >> ring;
	if (ring != '*' && !s.bad()) {
	  rs[0] = ring;
	  rs[1] = '\0';
	}
	state = (s.bad() ? kEnd : kReadLbrack);
      } break;
      case kReadLbrack: { // Try to read a left bracket 
	Char_t lbrack;
	s >> lbrack;
	state = (s.bad() ? kEnd : kReadSector);
      } break;
      case kReadSector: { // Try to read a sector 
	if (tmp == '*') s.get();
	else {
	  UShort_t sec;
	  s >> sec;
	  if (!s.bad()) {
	    minSector = sec;
	    maxSector = sec + 1;
	  }
	}
	state = (s.bad() ? kEnd : kReadComma);
      } break;
      case kReadComma: { // Try to read a left bracket 
	Char_t comma;
	s >> comma;
	state = (s.bad() ? kEnd : kReadStrip);
      } break;
      case kReadStrip: { // Try to read a strip 
	if (tmp == '*') s.get();
	else {
	  UShort_t str;
	  s >> str;
	  if (!s.bad()) {
	    minStrip = str;
	    maxStrip = str + 1;
	  }
	}
	state = (s.bad() ? kEnd : kReadRbrack);
      } break;
      case kReadRbrack: { // Try to read a left bracket 
	Char_t rbrack;
	s >> rbrack;
	state = kEnd;
      } break;
      case kEnd: 
	break;
      }
    }
  }
  UShort_t* dp = ds;
  UShort_t  det;
  while ((det = *(dp++))) {

    Char_t* rp = rs;
    Char_t  ring;
    while ((ring = *(rp++))) {
      if (det == 1 && ring == 'O') continue;
      UShort_t min  = GetMinStrip(det, ring, 0, 0);
      UShort_t max  = GetMaxStrip(det, ring, 0, 0);
      std::cout << "FMD" << det << ring 
		<< "  Strip range: " 
		<< std::setw(3) << min << "," 
		<< std::setw(3) << max << std::endl;

      UShort_t nSec = ( ring == 'I' ? 20  :  40 );
      UShort_t nStr = ( ring == 'I' ? 512 : 256 );
      for (UShort_t sec = minSector; sec < maxSector && sec < nSec; sec++) {

	UShort_t rate = GetSampleRate(det, ring, sec, 0);
	std::cout << "FMD" << det << ring << "[" << std::setw(2) << sec 
		  << "] sample rate: " << rate << std::endl;

	if (!showStrips) continue;
	std::cout 
	  << "  Strip |     Pedestal      |    Gain    | ZS thr. | Address\n" 
	  << "--------+-------------------+------------+---------+---------" 
	  << std::endl;
        for (UShort_t str = minStrip; str < nStr && str < maxStrip; str++) {
	  if (str == minStrip) std::cout << std::setw(3) << sec << ",";
	  else std::cout << "    ";
	  std::cout << std::setw(3) << str << " | ";
	  if (IsDead(det, ring, sec, str)) {
	    std::cout << "dead" << std::endl;
	    continue;
	  }
	  UShort_t ddl, addr, time, sam=0;
	  Detector2Hardware(det, ring, sec, str, sam, ddl, addr, time);
	  std::cout << std::setw(7) << GetPedestal(det, ring, sec, str) 
		    << "+/-" << std::setw(7) 
		    << GetPedestalWidth(det, ring, sec, str) 
		    << " | " << std::setw(10) 
		    << GetPulseGain(det, ring, sec, str) 
		    << " | " << std::setw(7) 
		    << GetZeroSuppression(det, ring, sec, str) 
		    << " | 0x" << std::hex << std::setw(4) 
		    << std::setfill('0') << ddl << ",0x" << std::setw(3) 
		    << addr << std::dec << std::setfill(' ') << std::endl;
        } // for (strip)
      } // for (sector)
      std::cout
	<< "=============================================================" 
	<< std::endl;
    } // while (ring)
  } // while (det)
  
}

//__________________________________________________________________
AliCDBEntry*
AliFMDParameters::GetEntry(const char* path, AliFMDPreprocessor* pp, 
			   Bool_t fatal) const
{
  // 
  // Get an entry from either global AliCDBManager or passed
  // AliFMDPreprocessor. 
  // 
  // Parameters:
  //    path  Path to CDB object. 
  //    pp    AliFMDPreprocessor 
  //    fatal If true, raise a fatal flag if we didn't get the entry.
  // Return:
  //    AliCDBEntry if found 
  // 
  AliCDBEntry* entry = 0;
  if (!pp) {
    AliCDBManager* cdb = AliCDBManager::Instance();
    entry              = cdb->Get(path);
  }
  else {
    const char* third  = gSystem->BaseName(path);
    const char* second = gSystem->BaseName(gSystem->DirName(path));
    entry              = pp->GetFromCDB(second, third);
  }
  if (!entry) { 
    TString msg(Form("No %s found in CDB, perhaps you need to "
		     "use AliFMDCalibFaker?", path));
    if (fatal) { AliFatal(msg.Data()); }
    else       AliLog::Message(AliLog::kWarning, msg.Data(), "FMD", 
			       "AliFMDParameters", "GetEntry", __FILE__, 
			       __LINE__);
    return 0;
  }
  if (entry && AliLog::GetDebugLevel("FMD", "") > 0) { 
    AliInfoF("Got entry %p for %s", entry, path);
    entry->PrintId();
    entry->PrintMetaData();			
    entry->Print();
  }
  return entry;
}

    
//__________________________________________________________________
void
AliFMDParameters::InitPulseGain(AliFMDPreprocessor* pp)
{
  // 
  // Initialize gains.  Try to get them from CDB 
  // 
  // Parameters:
  //    pp Pre-processor if called from shuttle
  //
  AliCDBEntry*   gain     = GetEntry(fgkPulseGain, pp);
  if (!gain) return;
  
  AliFMDDebug(5, ("Got gain from CDB"));
  fPulseGain = dynamic_cast<AliFMDCalibGain*>(gain->GetObject());
  if (!fPulseGain) AliFatal("Invalid pulser gain object from CDB");
  if (!fPulseGain->Values().Ptr()) 
    AliFatal("Empty pulser gain object from CDB");
}
//__________________________________________________________________
void
AliFMDParameters::InitPedestal(AliFMDPreprocessor* pp)
{
  //
  // Initialize pedestals.  Try to get them from CDB
  // 
  // Parameters:
  //    pp Pre-processor if called from shuttle
  //
  AliCDBEntry*   pedestal = GetEntry(fgkPedestal, pp);
  if (!pedestal) return;

  AliFMDDebug(5, ("Got pedestal from CDB"));
  fPedestal = dynamic_cast<AliFMDCalibPedestal*>(pedestal->GetObject());
  if (!fPedestal) AliFatal("Invalid pedestal object from CDB");
  if (!fPedestal->Values().Ptr()) AliFatal("Empty pedestal object from CDB");
}

//__________________________________________________________________
void
AliFMDParameters::InitDeadMap(AliFMDPreprocessor* pp)
{
  //
  // Initialize dead map.  Try to get it from CDB
  // 
  // Parameters:
  //    pp Pre-processor if called from shuttle
  //
  AliCDBEntry*   deadMap  = GetEntry(fgkDead, pp);
  if (!deadMap) return;
  
  AliFMDDebug(5, ("Got dead map from CDB"));
  fDeadMap = dynamic_cast<AliFMDCalibDeadMap*>(deadMap->GetObject());
  if (!fDeadMap) AliFatal("Invalid dead map object from CDB");
  if (!fDeadMap->Ptr()) AliFatal("Empty dead map object from CDB");
}

//__________________________________________________________________
void
AliFMDParameters::InitZeroSuppression(AliFMDPreprocessor* pp)
{
  //
  // Initialize zero suppression thresholds.  Try to get them from CDB
  // 
  // Parameters:
  //    pp Pre-processor if called from shuttle
  //
  AliCDBEntry*   zeroSup  = GetEntry(fgkZeroSuppression, pp);
  if (!zeroSup) return;
  AliFMDDebug(5, ("Got zero suppression from CDB"));
  fZeroSuppression = 
    dynamic_cast<AliFMDCalibZeroSuppression*>(zeroSup->GetObject());
  if (!fZeroSuppression)AliFatal("Invalid zero suppression object from CDB");
  if (!fZeroSuppression->Ptr()) {
    AliWarningF("Empty zero suppression object from CDB, assuming %d",
		fFixedZeroSuppression);
    AliCDBManager* cdbMan = AliCDBManager::Instance();
    if(!cdbMan || !cdbMan->GetCacheFlag())
      delete fZeroSuppression;
    fZeroSuppression = 0;
  }
}

//__________________________________________________________________
void
AliFMDParameters::InitSampleRate(AliFMDPreprocessor* pp)
{
  //
  // Initialize sample rates.  Try to get them from CDB
  // 
  // Parameters:
  //    pp Pre-processor if called from shuttle
  //
  AliCDBEntry*   sampRat  = GetEntry(fgkSampleRate, pp);
  if (!sampRat) return;
  AliFMDDebug(5, ("Got zero suppression from CDB"));
  fSampleRate = dynamic_cast<AliFMDCalibSampleRate*>(sampRat->GetObject());
  if (!fSampleRate) AliFatal("Invalid sample rate object from CDB");
  if (!fSampleRate->Rates().Ptr()) 
    AliFatal("empty sample rate object from CDB");
}

//__________________________________________________________________
void
AliFMDParameters::InitAltroMap(AliFMDPreprocessor* pp)
{
  //
  // Initialize hardware map.  Try to get it from CDB
  // 
  // Parameters:
  //    pp Pre-processor if called from shuttle
  //
  if (fAltroMap) { 
    delete fAltroMap;
    fAltroMap = 0;
  }
  AliCDBEntry*   hwMap    = GetEntry(fgkAltroMap, pp, kFALSE);       
  if (!hwMap) return;

  AliFMDDebug(5, ("Got ALTRO map from CDB"));
  fAltroMap = dynamic_cast<AliFMDAltroMapping*>(hwMap->GetObject());
  if (!fAltroMap) {
    AliFatal("Invalid ALTRO map object from CDB");
    fAltroMap = new AliFMDAltroMapping;
  }
}

//__________________________________________________________________
void
AliFMDParameters::InitStripRange(AliFMDPreprocessor* pp)
{
  //
  // Initialize strip range.  Try to get it from CDB
  // 
  // Parameters:
  //    pp Pre-processor if called from shuttle
  //
  AliCDBEntry*   range    = GetEntry(fgkStripRange, pp);
  if (!range) return;
  AliFMDDebug(5, ("Got strip range from CDB"));
  fStripRange = dynamic_cast<AliFMDCalibStripRange*>(range->GetObject());
  if (!fStripRange) AliFatal("Invalid strip range object from CDB");
  if (!fStripRange->Ranges().Ptr()) 
    AliFatal("Empty strip range object from CDB");
}


//__________________________________________________________________
Float_t
AliFMDParameters::GetThreshold() const
{
  // 
  // Get the threshold in the pulser gain 
  // 
  // 
  // Return:
  //    Threshold from pulser 
  //
  if (!fPulseGain) return fFixedThreshold;
  return fPulseGain->Threshold();
}

//__________________________________________________________________
Float_t
AliFMDParameters::GetPulseGain(UShort_t detector, Char_t ring, 
			       UShort_t sector, UShort_t strip) const
{
  // 
  // Gain of pre-amp. for strip, sector, ring, detector 
  //
  // For simulations this is normally set to 
  //
  // @f[ 
  //  \frac{\mbox{VA1_MIP_Range}{\mbox{ALTRO_channel_size}}\mbox{MIP_Energy_Loss}
  // @f]
  // 
  // 
  // Parameters:
  //    detector Detector # (1-3)
  //    ring     Ring ID ('I' or 'O')
  //    sector   Sector number (0-39)
  //    strip    Strip number (0-511)
  //
  // Return:
  //    Gain of pre-amp.  
  //
  if (!fPulseGain) { 
    if (fFixedPulseGain <= 0)
      fFixedPulseGain = fVA1MipRange * GetEdepMip() / fAltroChannelSize;
    return fFixedPulseGain;
  }  
  AliFMDDebug(50, ("pulse gain for FMD%d%c[%2d,%3d]=%f",
		    detector, ring, sector, strip,
		    fPulseGain->Value(detector, ring, sector, strip)));
  return fPulseGain->Value(detector, ring, sector, strip);
}

//__________________________________________________________________
Bool_t
AliFMDParameters::IsDead(UShort_t detector, Char_t ring, 
			 UShort_t sector, UShort_t strip) const
{
  // 
  // Whether the strip is considered dead
  // 
  // Parameters:
  //    detector Detector # (1-3)
  //    ring     Ring ID ('I' or 'O')
  //    sector   Sector number (0-39)
  //    strip    Strip number (0-511)
  //
  // Return:
  //    @c true if the strip is considered dead, @c false if it's
  // OK.
  //
  if (!fDeadMap) return kFALSE;
  AliFMDDebug(50, ("Dead for FMD%d%c[%2d,%3d]=%s",
		    detector, ring, sector, strip,
		    fDeadMap->operator()(detector, ring, sector, strip) ? 
		    "no" : "yes"));
  return fDeadMap->operator()(detector, ring, sector, strip);
}

//__________________________________________________________________
UShort_t
AliFMDParameters::GetZeroSuppression(UShort_t detector, Char_t ring, 
				     UShort_t sector, UShort_t strip) const
{
  // 
  // zero suppression threshold (in ADC counts)
  // 
  // Parameters:
  //    detector Detector # (1-3)
  //    ring     Ring ID ('I' or 'O')
  //    sector   Sector number (0-39)
  //    strip    Strip number (0-511)
  //
  // Return:
  //    zero suppression threshold (in ADC counts) 
  //
  if (!fZeroSuppression) return fFixedZeroSuppression;

  // In case of empty zero suppression objects. 
  if (!fZeroSuppression->Ptr() || 
      fZeroSuppression->MaxIndex() <= 0) return fFixedZeroSuppression;

  // Need to map strip to ALTRO chip. 
  AliFMDDebug(50, ("zero sup. for FMD%d%c[%2d,%3d]=%d",
		    detector, ring, sector, strip,
		    fZeroSuppression->operator()(detector, ring, 
						 sector, strip)));
  return fZeroSuppression->operator()(detector, ring, sector, strip/128);
}

//__________________________________________________________________
UShort_t
AliFMDParameters::GetSampleRate(UShort_t det, Char_t ring, UShort_t sector, 
				UShort_t str) const
{
  // 
  // Get the sampling rate
  // 
  // Parameters:
  //    detector Detector # (1-3)
  //    ring     Ring ID ('I' or 'O')
  //    sector   Sector number (0-39)
  //    strip    Strip number (0-511)
  //
  // Return:
  //    The sampling rate 
  //
  if (!fSampleRate) return fFixedSampleRate;
  // Need to map sector to digitizier card. 
  UInt_t ret = fSampleRate->Rate(det, ring, sector, str);
  AliFMDDebug(50, ("Sample rate for FMD%d%c[%2d,%3d]=%d", 
		    det, ring, sector, str, ret));
  return ret;
}

//__________________________________________________________________
UShort_t
AliFMDParameters::GetMinStrip(UShort_t det, Char_t ring, UShort_t sector, 
			      UShort_t str) const
{
  // 
  // Get the minimum strip in the read-out range
  // 
  // Parameters:
  //    detector Detector # (1-3)
  //    ring     Ring ID ('I' or 'O')
  //    sector   Sector number (0-39)
  //    strip    Strip number (0-511)
  //
  // Return:
  //    Minimum strip 
  //
  if (!fStripRange) return fFixedMinStrip;
  // Need to map sector to digitizier card. 
  UInt_t ret = fStripRange->Min(det, ring, sector, str);
  AliFMDDebug(50, ("Min strip # for FMD%d%c[%2d,%3d]=%d", 
		    det, ring, sector, str, ret));
  return ret;
}

//__________________________________________________________________
UShort_t
AliFMDParameters::GetMaxStrip(UShort_t det, Char_t ring, UShort_t sector, 
			      UShort_t str) const
{
  // 
  // Get the maximum strip in the read-out range
  // 
  // Parameters:
  //    detector Detector # (1-3)
  //    ring     Ring ID ('I' or 'O')
  //    sector   Sector number (0-39)
  //    strip    Strip number (0-511)
  //
  // Return:
  //    Maximum strip 
  //
  if (!fStripRange) return fFixedMaxStrip;
  // Need to map sector to digitizier card. 
  UInt_t ret = fStripRange->Max(det, ring, sector, str);
  AliFMDDebug(50, ("Max strip # for FMD%d%c[%2d,%3d]=%d", 
		    det, ring, sector, str, ret));
  return ret;
}

//__________________________________________________________________
Float_t
AliFMDParameters::GetPedestal(UShort_t detector, Char_t ring, 
			      UShort_t sector, UShort_t strip) const
{
  // 
  // Get mean of pedestal
  // 
  // Parameters:
  //    detector Detector # (1-3)
  //    ring     Ring ID ('I' or 'O')
  //    sector   Sector number (0-39)
  //    strip    Strip number (0-511)
  //
  // Return:
  //    Mean of pedestal 
  //
  if (!fPedestal) return fFixedPedestal;
  AliFMDDebug(50, ("pedestal for FMD%d%c[%2d,%3d]=%f",
		    detector, ring, sector, strip,
		    fPedestal->Value(detector, ring, sector, strip)));
  return fPedestal->Value(detector, ring, sector, strip);
}

//__________________________________________________________________
Float_t
AliFMDParameters::GetPedestalWidth(UShort_t detector, Char_t ring, 
				   UShort_t sector, UShort_t strip) const
{
  // 
  // Width of pedestal
  // 
  // Parameters:
  //    detector Detector # (1-3)
  //    ring     Ring ID ('I' or 'O')
  //    sector   Sector number (0-39)
  //    strip    Strip number (0-511)
  //
  // Return:
  //    Width of pedestal 
  //
  if (!fPedestal) return fFixedPedestalWidth;
  AliFMDDebug(50, ("pedetal width for FMD%d%c[%2d,%3d]=%f",
		    detector, ring, sector, strip,
		    fPedestal->Width(detector, ring, sector, strip)));
  return fPedestal->Width(detector, ring, sector, strip);
}
  
//__________________________________________________________________
AliFMDAltroMapping*
AliFMDParameters::GetAltroMap() const
{
  // 
  // Get the map that translates hardware to detector coordinates 
  //
  // Return:
  //    Get the map that translates hardware to detector
  // coordinates 
  // 
  return fAltroMap;
}


//____________________________________________________________________
Bool_t 
AliFMDParameters::Hardware2Detector(UShort_t  ddl,       UShort_t addr,
				    UShort_t  timebin,   
				    UShort_t& det,       Char_t&  ring, 
				    UShort_t& sec,       Short_t& str,
				    UShort_t& sam) const
{
  // 
  // Map a hardware address into a detector index. 
  // 
  // Parameters:
  //    ddl        Hardware DDL number 
  //    addr       Hardware address.  
  //    timebin    Timebin 
  //    det        On return, the detector #
  //    ring       On return, the ring ID
  //    sec        On return, the sector #
  //    str        On return, the base of strip #
  //    sam        On return, the sample number for this strip
  //
  // Return:
  //    @c true on success, false otherwise 
  //
  if (!fAltroMap) return kFALSE;
  UShort_t board, chip, chan;
  fAltroMap->ChannelAddress(addr, board, chip, chan);
  return Hardware2Detector(ddl,board,chip,chan,timebin,det,ring,sec,str,sam);
}
//____________________________________________________________________
Bool_t 
AliFMDParameters::Hardware2Detector(UShort_t  ddl,       UShort_t   board,
				    UShort_t  chip,      UShort_t   chan,
				    UShort_t  timebin,   
				    UShort_t& det,       Char_t&   ring, 
				    UShort_t& sec,       Short_t& str,
				    UShort_t& sam) const
{
  // 
  // Map a hardware address into a detector index. 
  // 
  // Parameters:
  //    ddl        Hardware DDL number 
  //    board      FEC number
  //    altro      ALTRO number 
  //    channel    Channel number 
  //    timebin    Timebin 
  //    det        On return, the detector #
  //    ring       On return, the ring ID
  //    sec        On return, the sector #
  //    str        On return, the base of strip #
  //    sam        On return, the sample number for this strip
  //
  // Return:
  //    @c true on success, false otherwise 
  //
  if (!fAltroMap) {
    AliFMDDebug(1, ("No ALTRO map available"));
    return kFALSE;
  }
  if (fAltroMap->DDL2Detector(ddl) < 0) { 
    AliFMDDebug(1, ("Invalid DDL number %d", ddl));
    return kFALSE;
  }
  det = fAltroMap->DDL2Detector(ddl);
  Short_t stripBase = 0;
  if (!fAltroMap->Channel2StripBase(board,chip,chan, ring, sec, stripBase)) {
    AliFMDDebug(1, ("Failed to translate  "
		    "%d/0x%02x/0x%x/0x%x/%04d -> "
		    "FMD%d%c[%2d,%3d] to detector", 
		    ddl, board, chip, chan, timebin, 
		    det, ring, sec, stripBase));
    return kFALSE;
  }
  UShort_t preSamples = GetPreSamples(det, ring, sec, stripBase);
  UShort_t sampleRate = GetSampleRate(det, ring, sec, stripBase);
  Short_t stripOff = 0;
  fAltroMap->Timebin2Strip(sec, timebin, preSamples, sampleRate, stripOff, sam);
  str = stripBase + stripOff;
  AliFMDDebug(50, ("%d/0x%02x/0x%x/0x%x/%04d -> FMD%d%c[%02d,%03d]-%d"
		  " (pre=%2d, rate=%d)", 
		   ddl, board, chip, chan, timebin, 
		   det, ring, sec, str, sam, preSamples, sampleRate));
  return kTRUE;
}


//____________________________________________________________________
Bool_t 
AliFMDParameters::Detector2Hardware(UShort_t  det,        Char_t    ring, 
				    UShort_t  sec,        UShort_t  str,
				    UShort_t  sam, 
				    UShort_t& ddl,        UShort_t& board, 
				    UShort_t& altro,      UShort_t& channel, 
				    UShort_t& timebin) const
{
  // 
  // Map a detector index into a hardware address. 
  // 
  // Parameters:
  //    det         The detector #
  //    ring        The ring ID
  //    sec         The sector #
  //    str         The strip #
  //    sam         The sample number 
  //    ddl         On return, hardware DDL number 
  //    board       On return, the FEC board address (local to DDL)
  //    altro       On return, the ALTRO number (local to FEC)
  //    channel     On return, the channel number (local to ALTRO)
  //    timebin     On return, the timebin number (local to ALTRO)
  //
  // Return:
  //    @c true on success, false otherwise 
  //
  if (!fAltroMap) { 
    AliFMDDebug(1, ("No ALTRO map available"));
    return kFALSE;
  }
  UShort_t preSamples = GetPreSamples(det, ring, sec, str);
  UShort_t sampleRate = GetSampleRate(det, ring, sec, str);
  UShort_t strip      = str - GetMinStrip(det,ring,sec,str);
  return fAltroMap->Detector2Hardware(det, ring, sec, strip, sam,
				      preSamples, sampleRate,
				      ddl, board, altro, channel, timebin);
}

  

//____________________________________________________________________
Bool_t 
AliFMDParameters::Detector2Hardware(UShort_t  det,        Char_t    ring, 
				    UShort_t  sec,        UShort_t  str,
				    UShort_t  sam, 
				    UShort_t&   ddl,        UShort_t&   addr,
				    UShort_t& timebin) const
{
  // 
  // Map a detector index into a hardware address. 
  // 
  // Parameters:
  //    det         The detector #
  //    ring        The ring ID
  //    sec         The sector #
  //    str         The strip #
  //    sam         The sample number 
  //    ddl         On return, hardware DDL number 
  //    addr      On return, hardware address.  
  //    timebin     On return, the timebin number (local to ALTRO)
  //
  // Return:
  //    @c true on success, false otherwise 
  //
  if (!fAltroMap) return kFALSE;
  UShort_t preSamples = GetPreSamples(det, ring, sec, str);
  UShort_t sampleRate = GetSampleRate(det, ring, sec, str);
  UShort_t strip      = str - GetMinStrip(det,ring,sec,str);
  return fAltroMap->Detector2Hardware(det, ring, sec, strip, sam,
				      preSamples, sampleRate,
				      ddl, addr, timebin);
}


//__________________________________________________________________
Float_t
AliFMDParameters::GetEdepMip() const 
{ 
  // 
  // Return:
  //    The average energy deposited by one MIP 
  //
  if (fEdepMip <= 0){
    AliFMDGeometry* fmd = AliFMDGeometry::Instance();
    fEdepMip = (fkSiDeDxMip 
		* fmd->GetRing('I')->GetSiThickness() 
		* fmd->GetSiDensity());
  }
  return fEdepMip;
}
//____________________________________________________________________
Float_t  
AliFMDParameters::GetDACPerMIP() const
{
  // 
  // This is the conversion from Digital-to-Analog-Converter setting
  // to the number of MIPs. The number was measured in the NBI lab during
  // August 2008.
  //
  // Return:
  //    The conversion factor from DAC to ADC 
  //
  return 29.67;
  
}
 
//____________________________________________________________________
//
// EOF
//
