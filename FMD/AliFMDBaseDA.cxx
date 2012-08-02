/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
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
 **************************************************************************/

/** @file    AliFMDBaseDA.cxx
    @author  Hans Hjersing Dalsgaard <canute@nbi.dk>
    @date    Wed Mar 26 11:30:45 2008
    @brief   Base class for detector algorithms.
*/
//
// This is the implementation of the (virtual) base class for the FMD
// detector algorithms(DA). It implements the creation of the relevant
// containers and handles the loop over the raw data. The derived
// classes can control the parameters and action to be taken making
// this the base class for the Pedestal, Gain and Physics DA.
//

#include "AliFMDBaseDA.h"
#include "AliRawReader.h"
#include "AliFMDDigit.h"
#include "AliFMDParameters.h"
#include "AliFMDRawReader.h"
#include "AliFMDCalibSampleRate.h"
#include "AliFMDCalibStripRange.h"
#include "AliLog.h"
#include "AliRawEventHeaderBase.h"
#include "AliFMDDigit.h"
#include <TClonesArray.h>
#include <TFile.h>
#include <TDatime.h>
#include <TSystem.h>
#include <TH2F.h>
#include <TStopwatch.h>
#include <TROOT.h>
#include <TPluginManager.h>
#include <iostream>
#include <fstream>

//_____________________________________________________________________
ClassImp(AliFMDBaseDA)
#if 0 
; // Do not delete  - to let Emacs for mat the code
#endif

//_____________________________________________________________________
const char*
AliFMDBaseDA::GetStripPath(UShort_t det, 
			   Char_t   ring, 
			   UShort_t sec, 
			   UShort_t str, 
			   Bool_t   full) const
{
  // Get the strip path 
  // 
  // Parameters 
  //     det      Detector number
  //     ring     Ring identifier 
  //     sec      Sector number 
  //     str      Strip number
  //     full     If true, return full path 
  // 
  // Return 
  //     The path
  return Form("%s%sFMD%d%c[%02d,%03d]", 
	      (full ? GetSectorPath(det, ring, sec, full) : ""), 
	      (full ? "/" : ""), det, ring, sec, str);
}
//_____________________________________________________________________
const char*
AliFMDBaseDA::GetSectorPath(UShort_t det, 
			    Char_t   ring, 
			    UShort_t sec, 
			    Bool_t   full) const
{
  // Get the strip path 
  // 
  // Parameters 
  //     det      Detector number
  //     ring     Ring identifier 
  //     sec      Sector number 
  //     str      Strip number
  //     full     If true, return full path 
  // 
  // Return 
  //     The path
  return Form("%s%sFMD%d%c[%02d]", 
	      (full ? GetRingPath(det, ring, full) : ""), 
	      (full ? "/" : ""), det, ring, sec);
}
//_____________________________________________________________________
const char*
AliFMDBaseDA::GetRingPath(UShort_t det, 
			  Char_t   ring, 
			  Bool_t   full) const
{
  // Get the strip path 
  // 
  // Parameters 
  //     det      Detector number
  //     ring     Ring identifier 
  //     sec      Sector number 
  //     str      Strip number
  //     full     If true, return full path 
  // 
  // Return 
  //     The path
  return Form("%s%sFMD%d%c", 
	      (full ? GetDetectorPath(det, full) : ""), 
	      (full ? "/" : ""), det, ring);
}
//_____________________________________________________________________
const char*
AliFMDBaseDA::GetDetectorPath(UShort_t det, 
			      Bool_t   full) const
{
  // Get the strip path 
  // 
  // Parameters 
  //     det      Detector number
  //     ring     Ring identifier 
  //     sec      Sector number 
  //     str      Strip number
  //     full     If true, return full path 
  // 
  // Return 
  //     The path
  return Form("%s%sFMD%d", 
	      (full ? fDiagnosticsFilename.Data() : ""), 
	      (full ? ":/" : ""), det);
}

//_____________________________________________________________________
AliFMDBaseDA::AliFMDBaseDA() : 
  TNamed(),
  fDiagnosticsFilename("diagnosticsHistograms.root"),
  fOutputFile(),
  fConditionsFile(),
  fSaveHistograms(kFALSE),
  fMakeSummaries(kFALSE),
  fDetectorArray(),
  fPulseSize(10),
  fPulseLength(10),
  fRequiredEvents(0),
  fCurrentEvent(0), 
  fRunno(0),
  fSummaries(0)
{
  //Constructor
  fSeenDetectors[0] = fSeenDetectors[1] = fSeenDetectors[2] = kFALSE;
  fDetectorArray.SetOwner();
  Rotate("conditions.csv", 3);
  fConditionsFile.open("conditions.csv");
}
//_____________________________________________________________________
AliFMDBaseDA::AliFMDBaseDA(const AliFMDBaseDA & baseDA) : 
  TNamed(baseDA),
  fDiagnosticsFilename(baseDA.fDiagnosticsFilename.Data()),
  fOutputFile(),
  fConditionsFile(),
  fSaveHistograms(baseDA.fSaveHistograms),
  fMakeSummaries(baseDA.fMakeSummaries),
  fDetectorArray(baseDA.fDetectorArray),
  fPulseSize(baseDA.fPulseSize),
  fPulseLength(baseDA.fPulseLength),
  fRequiredEvents(baseDA.fRequiredEvents),
  fCurrentEvent(baseDA.fCurrentEvent),
  fRunno(baseDA.fRunno),
  fSummaries(0)
{
  //Copy constructor
  fSeenDetectors[0] = baseDA.fSeenDetectors[0];
  fSeenDetectors[1] = baseDA.fSeenDetectors[1];
  fSeenDetectors[2] = baseDA.fSeenDetectors[2];

  fDetectorArray.SetOwner();
  
}


//_____________________________________________________________________
AliFMDBaseDA::~AliFMDBaseDA() 
{
  //destructor
}

//_____________________________________________________________________
Bool_t AliFMDBaseDA::HaveEnough(Int_t nEvents) const
{
  return nEvents > GetRequiredEvents();
}
//_____________________________________________________________________
UShort_t AliFMDBaseDA::GetProgress(Int_t nEvents) const
{
  return UShort_t((nEvents *100)/ GetRequiredEvents());
}

//_____________________________________________________________________
void AliFMDBaseDA::Run(AliRawReader* reader) 
{
  //Run the FMD DA
  TFile* diagFile = 0;
  // if (fSaveHistograms)
  //  diagFile = TFile::Open(fDiagnosticsFilename.Data(),"RECREATE");
  
  reader->Reset();
  fRunno = reader->GetRunNumber();

  AliFMDRawReader* fmdReader  = new AliFMDRawReader(reader,0);
  TClonesArray*    digitArray = new TClonesArray("AliFMDDigit",0);
  
  Bool_t sodread = kFALSE;
  
  for(Int_t i=0;i<3;i++) {
    if (!reader->NextEvent()) { 
      // Read Start-of-Run / Start-of-Files event
      AliWarning(Form("Failed to read the %d%s event",
		      i+1, (i == 0 ? "st" : (i == 1 ? "nd" : "rd"))));
      break;
    }
    
    UInt_t eventType = reader->GetType();
    if(eventType == AliRawEventHeaderBase::kStartOfData || 
       eventType == AliRawEventHeaderBase::kFormatError) { 
      
      WriteConditionsData(fmdReader);
      Init();
      sodread = kTRUE;
      break;
    }
  }
  
  InitContainer(diagFile);
  if (AliLog::GetDebugLevel("FMD","") >= 3) { 
    fDetectorArray.ls();
  }
  
  if(!sodread) 
    AliWarning("No SOD event detected!");
  
  int lastProgress = 0;
  
  for(Int_t i = 0; i< 3;i++) fNEventsPerDetector[i] = 0;

  for(Int_t n = 1; !HaveEnough(n); n++) {
    if(!reader->NextEvent()) { n--; continue; }
    SetCurrentEvent(n);
    digitArray->Clear();
    fmdReader->ReadAdcs(digitArray);
    
    Bool_t seen[] = { false, false, false };
    for(Int_t i = 0; i<digitArray->GetEntriesFast();i++) {
      AliFMDDigit* digit = static_cast<AliFMDDigit*>(digitArray->At(i));
      UShort_t det = digit->Detector();
      fSeenDetectors[det-1] = true;
      seen[det-1]           = true;
      FillChannels(digit);
    }
    
    for(Int_t i = 0; i< 3;i++) 
      if (seen[i]) (fNEventsPerDetector[i])++;
      
    FinishEvent();
    
    int progress = GetProgress(n);
    if (progress <= lastProgress) continue;
    lastProgress = progress;
    std::cout << "Progress: " << lastProgress << " / 100 " << std::endl;
    
  }
  
  AliInfoF("Looped over %d events (%d,%d,%d)",GetCurrentEvent(),
	   fNEventsPerDetector[0], 
	   fNEventsPerDetector[1], 
	   fNEventsPerDetector[2]);
  WriteHeaderToFile();
  
  for(UShort_t det=1;det<=3;det++) {
    if (!fSeenDetectors[det-1]) continue;
    std::cout << "FMD" << det << std::endl;
    UShort_t firstRing = (det == 1 ? 1 : 0);
    for (UShort_t ir = firstRing; ir < 2; ir++) {
      Char_t   ring = (ir == 0 ? 'O' : 'I');
      UShort_t nsec = (ir == 0 ? 40  : 20);
      UShort_t nstr = (ir == 0 ? 256 : 512);

      if (fMakeSummaries) MakeSummary(det, ring);

      std::cout << " Ring " << ring << ": " << std::flush;
      for(UShort_t sec =0; sec < nsec;  sec++)  {
  	for(UShort_t strip = 0; strip < nstr; strip++) {
	  Analyse(det,ring,sec,strip);
  	}
	std::cout << '.' << std::flush;
      }
      // if(fSaveHistograms)
      // diagFile->Flush();
      std::cout << "done" << std::endl;
    }
  }
  
  if(fOutputFile.is_open()) {
    fOutputFile.write("# EOF\n",6);
    fOutputFile.close();
  }
  
  Terminate(diagFile);
    
  if(fSaveHistograms ) {
    diagFile = TFile::Open(fDiagnosticsFilename.Data(),"RECREATE");
    fDetectorArray.Write("FMD", TObject::kSingleKey);
    fSummaries.Write();
    AliInfo("Closing diagnostics file - please wait ...");
    // diagFile->Write();
    diagFile->Close();
    AliInfo("done");
    
  }
}
//_____________________________________________________________________

void AliFMDBaseDA::InitContainer(TDirectory* diagFile)
{
  //Prepare container for diagnostics    
  TDirectory* savDir   = gDirectory;

  for(UShort_t det=1;det<=3;det++) {
    TObjArray* detArray = new TObjArray(det == 1 ? 1 : 2);
    detArray->SetOwner();
    detArray->SetName(GetDetectorPath(det, false));
    fDetectorArray.AddAtAndExpand(detArray,det);

    TDirectory* detDir = 0;
    if (diagFile) {
      diagFile->cd();
      detDir = diagFile->mkdir(GetDetectorPath(det, kFALSE));
    }

    UShort_t FirstRing = (det == 1 ? 1 : 0);
    for (UShort_t ir = FirstRing; ir < 2; ir++) {
      Char_t   ring = (ir == 0 ? 'O' : 'I');
      UShort_t nsec = (ir == 0 ? 40  : 20);
      UShort_t nstr = (ir == 0 ? 256 : 512);
      TObjArray* ringArray = new TObjArray(nsec);
      ringArray->SetOwner();
      ringArray->SetName(GetRingPath(det, ring, false));
      detArray->AddAtAndExpand(ringArray,ir);


      TDirectory* ringDir = 0;
      if (detDir) { 
	detDir->cd();
	ringDir = detDir->mkdir(GetRingPath(det,ring, kFALSE));
      }
      

      for(UShort_t sec =0; sec < nsec;  sec++)  {
	TObjArray* sectorArray = new TObjArray(nstr);
	sectorArray->SetOwner();
	sectorArray->SetName(GetSectorPath(det, ring, sec, false));
	ringArray->AddAtAndExpand(sectorArray,sec);


	TDirectory* secDir = 0;
	if (ringDir) { 
	  ringDir->cd();
	  secDir = ringDir->mkdir(GetSectorPath(det, ring, sec, kFALSE));
	}
	
	for(UShort_t strip = 0; strip < nstr; strip++) {
	  if (secDir) { 
	    secDir->cd();
	    secDir->mkdir(GetStripPath(det, ring, sec, strip, kFALSE));
	  }
	  TObjArray* stripArray = new TObjArray(0);
	  stripArray->SetOwner(true);
	  stripArray->SetName(GetStripPath(det, ring, sec, strip, false));
	  sectorArray->AddAtAndExpand(stripArray, strip);
	  AddChannelContainer(stripArray, det, ring, sec, strip);
	}
	AddSectorSummary(sectorArray, det, ring, sec, nstr);
      }
    }
  }
  savDir->cd();
}

//_____________________________________________________________________ 
void AliFMDBaseDA::WriteConditionsData(AliFMDRawReader* fmdReader) 
{
  //Write the conditions data to file
  AliFMDParameters* pars       = AliFMDParameters::Instance();
  fConditionsFile.write(Form("# %s \n",pars->GetConditionsShuttleID()),14);
  TDatime now;
  fConditionsFile << "# This file created from run number " << fRunno 
		  << " at " << now.AsString() << std::endl;
  
  AliFMDCalibSampleRate* sampleRate = new AliFMDCalibSampleRate();
  AliFMDCalibStripRange* stripRange = new AliFMDCalibStripRange();
  
  fmdReader->ReadSODevent(sampleRate,stripRange,fPulseSize,fPulseLength,
			  fSeenDetectors);

  sampleRate->WriteToFile(fConditionsFile, fSeenDetectors);
  stripRange->WriteToFile(fConditionsFile, fSeenDetectors);

 
  // Zero Suppresion
  
  // Strip Range
  
  fConditionsFile.write("# Gain Events \n",15);
  
  for(UShort_t det=1; det<=3;det++) {
    if (!fSeenDetectors[det-1]) { 
      continue;
    }
    UShort_t firstring = (det == 1 ? 1 : 0);
    for(UShort_t iring = firstring; iring <=1;iring++) {
      Char_t ring = (iring == 1 ? 'I' : 'O');
      for(UShort_t board =0 ; board <=1; board++) {
	
	Int_t idx = GetHalfringIndex(det,ring,board);
	
	fConditionsFile << det                     << ','
			<< ring                    << ','
			<< board                   << ','
			<< fPulseLength.At(idx)    << "\n";
	
      }
    }
  }
  
  fConditionsFile.write("# Gain Pulse \n",14);
  
  for(UShort_t det=1; det<=3;det++) {
    if (!fSeenDetectors[det-1]) { 
      continue;
    }
    UShort_t firstring = (det == 1 ? 1 : 0);
    for(UShort_t iring = firstring; iring <=1;iring++) {
      Char_t ring = (iring == 1 ? 'I' : 'O');
      for(UShort_t board =0 ; board <=1; board++) {
	
	Int_t idx = GetHalfringIndex(det,ring,board);
	
	fConditionsFile << det                     << ','
			<< ring                    << ','
			<< board                   << ','
			<< fPulseSize.At(idx)      << "\n";
	
      }
    }
  }
  // sampleRate->WriteToFile(std::cout, fSeenDetectors);
  // stripRange->WriteToFile(std::cout, fSeenDetectors);

  if(fConditionsFile.is_open()) {
    
    fConditionsFile.write("# EOF\n",6);
    fConditionsFile.close();
    
  }
  
}
//_____________________________________________________________________ 
Int_t AliFMDBaseDA::GetHalfringIndex(UShort_t det, Char_t ring, 
				     UShort_t board) const 
{
  // Get the index corresponding to a half-ring 
  // 
  // Parameters: 
  //   det    Detector number 
  //   ring   Ring identifier 
  //   board  Board number 
  //
  // Return 
  //   Internal index of the board 
  UShort_t iring  =  (ring == 'I' ? 1 : 0);
  
  Int_t index = (((det-1) << 2) | (iring << 1) | (board << 0));
  
  return index-2;
  
}
//_____________________________________________________________________ 
void AliFMDBaseDA::Rotate(const char* base, int max) const
{
  // 
  // Rotate a set of files.   base is the basic name of the files.
  // If the file base.max exists it is removed. 
  // If the file base.n exists (where n < max) it is renamed to
  // base.(n-1).  
  // If the file base exists, it is renamed to base.1 
  //
  // Parameters:
  //   base Base name of the files
  //   max  Maximum number to keep (minus one for the current).

  // Note:  TSystem::AccessPathName returns false if the condition is
  // fulfilled! 

  // Check if we have base.max, and if so, remove it. 
  TString testName(Form("%s.%d", base, max));
  if (!gSystem->AccessPathName(testName.Data())) 
    gSystem->Unlink(testName.Data());
    
  // Loop down from max-1 to 1 and move files 
  for (int i = max-1; i >= 1; i--) { 
    testName = Form("%s.%d", base, i);
    if (!gSystem->AccessPathName(testName.Data())) {
      TString newName(Form("%s.%d", base, i+1));
      gSystem->Rename(testName.Data(), newName.Data());
    }
  }

  // If we have the file base, rename it to base.1 
  testName = Form("%s", base);
  if (!gSystem->AccessPathName(testName.Data())){
    TString newName(Form("%s.%d", base, 1));
    gSystem->Rename(testName.Data(), newName.Data());
  }
}

//_____________________________________________________________________ 
TH2*
AliFMDBaseDA::MakeSummaryHistogram(const char* prefix, const char* title, 
				   UShort_t d, Char_t r) 
{
  // 
  // Utility function for defining summary histograms 
  // 
  // Parameters:
  //    det    Detector 
  //    ring   Ring identifier 
  //    prefix Histogram prefix 
  //    title  Histogram title 
  //
  Int_t nX = ((d == 1 || r == 'I' || r == 'i') ?  20 :  40);
  Int_t nY = ((d == 1 || r == 'I' || r == 'i') ? 512 : 256);
  
  TH2* ret = new TH2F(Form("%sFMD%d%c", prefix, d, r), 
		      Form("%s for FMD%d%c", title, d, r), 
		      nX, -0.5, nX-0.5, nY, -0.5, nY-0.5);
  ret->SetXTitle("Sector #");
  ret->SetYTitle("Strip #");
  ret->SetDirectory(0);
  // if (!fSummaries) fSummaries = new TObjArray;
  fSummaries.Add(ret);
  return ret;
}

//_____________________________________________________________________ 
TObjArray*
AliFMDBaseDA::GetDetectorArray(UShort_t det)
{
  if (det < 1 || det > 3) { 
    AliErrorF("Detector index %d out of bounds", det);
    return 0;
  }
  return static_cast<TObjArray*>(fDetectorArray.At(det));
}
//_____________________________________________________________________ 
TObjArray*
AliFMDBaseDA::GetRingArray(UShort_t det, Char_t ring)
{
  Int_t idx = (ring == 'O' || ring == 'o' ? 0 : 1);
  TObjArray* array = GetDetectorArray(det);
  if (!array) return 0;
  array = static_cast<TObjArray*>(array->At(idx));
  if (!array) AliErrorF("No ring array for FMD%d%c (%d)", det, ring, idx);
  return array;
}
//_____________________________________________________________________ 
TObjArray*
AliFMDBaseDA::GetSectorArray(UShort_t det, Char_t ring, UShort_t sector)
{
  TObjArray* array = GetRingArray(det,ring);
  if (!array) return 0;
  array = static_cast<TObjArray*>(array->At(sector));
  if (!array) AliErrorF("No sector array for FMD%d%c[%02d]", det, ring, sector);
  return array;
}
//_____________________________________________________________________ 
TObjArray*
AliFMDBaseDA::GetStripArray(UShort_t det, Char_t ring, 
			    UShort_t sector, UShort_t strip)
{
  TObjArray* array = GetSectorArray(det,ring,sector);
  if (!array) return 0;
  array = static_cast<TObjArray*>(array->At(strip));
  if (!array) AliErrorF("No strip array for FMD%d%c[%02d,%03d]", 
			det, ring, sector, strip);
  return array;
}

//=====================================================================
AliFMDBaseDA::Runner::Runner()
  : fReader(0),
    fDiagFile(""), 
    fDiag(false)
{}

//_____________________________________________________________________ 
void
AliFMDBaseDA::Runner::AddHandlers()
{
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()");
  gROOT->GetPluginManager()->AddHandler("ROOT::Math::Minimizer", "Minuit", 
					"TMinuitMinimizer",
					"Minuit", 
					"TMinuitMinimizer(const char *)");
  gROOT->GetPluginManager()->AddHandler("ROOT::Math::Minimizer", 
					"GSLMultiMin", 
					"ROOT::Math::GSLMinimizer",
					"MathMore", 
					"GSLMinimizer(const char *)");
  gROOT->GetPluginManager()->AddHandler("ROOT::Math::Minimizer", 
					"GSLMultiFit", 
					"ROOT::Math::GSLNLSMinimizer",
					"MathMore", "GSLNLSMinimizer(int)");
  gROOT->GetPluginManager()->AddHandler("ROOT::Math::Minimizer", 
					"GSLSimAn", 
					"ROOT::Math::GSLSimAnMinimizer",
					"MathMore", 
					"GSLSimAnMinimizer(int)");
  gROOT->GetPluginManager()->AddHandler("ROOT::Math::Minimizer", 
					"Linear", 
					"TLinearMinimizer",
					"Minuit", 
					"TLinearMinimizer(const char *)");
  gROOT->GetPluginManager()->AddHandler("ROOT::Math::Minimizer", 
					"Fumili", 
					"TFumiliMinimizer",
					"Fumili", 
					"TFumiliMinimizer(int)");
}
//_____________________________________________________________________ 
void
AliFMDBaseDA::Runner::ShowUsage(std::ostream& o, const char* progname)
{
  o << "Usage: " << progname << " SOURCE [OPTIONS]\n\n"
    << "Options:\n"
    << "   -h,--help                Show this help\n"
    << "   -d,--diagnostics[=FILE]  Write diagnostics to file\n"
    << "   -D,--debug=LEVEL         Set the debug level\n\n"
    << "SOURCE is one of\n"
    << " * FILE.raw                Raw data file\n"
    << " * FILE.root               ROOT'ified raw data file\n"
    << " * collection://FILE.root  Event list in a ROOT file\n"
    << " * collection://FILE       File containing list of ROOT files\n"
    << " * ^FMD                    Monitor source\n"
    << "There are other options too.  Check the AliRawReader docs\n"
    << std::endl;
}
 
//_____________________________________________________________________ 
namespace {
  Bool_t ExtractValue(const TString& arg, TString& val)
  {
    val = "";
    Int_t eq = arg.Index("=");
    if (eq == kNPOS) return false;
    
    val = arg(eq+1, arg.Length()-eq-1);
    return true;
  }
}
      
//_____________________________________________________________________ 
Int_t
AliFMDBaseDA::Runner::Init(int argc, char** argv)
{
  AddHandlers();

  // --- Process the command line ------------------------------------
  TString source;
  Int_t   debugLevel  = 0;
  Bool_t  help        = false;
  
  for (int i = 1; i < argc; i++) { 
    TString arg(argv[i]);
    Bool_t  badOption   = false;
    
    if (arg[0] == '-') { // It's an option 
      if (arg[1] == '-') { // It's a long option 
	TString val;
	if      (arg.EqualTo("--help"))     help = true; 
	else if (arg.BeginsWith("--debug")) {
	  if (ExtractValue(arg, val))
	    debugLevel = val.Atoi();
	}
	else if (arg.BeginsWith("--diagnostics")) { 
	  fDiag = true;
	  if (ExtractValue(arg, val)) 
	    fDiagFile = val;
	}
	else badOption = true;
      }
      else { // Short option 
	TString next(i < argc-1 ? argv[i+1] : "");
	switch (arg[1]) { 
	case 'h': help = true; break;
	case 'd': fDiag = true; 
	  if (!next.IsNull() && next[0] != '-') {
	    fDiagFile = next;
	    i++;
	  }
	  break;
	case 'D': 
	  if (!next.IsNull() && next[0] != '-') {
	    debugLevel = next.Atoi();
	    i++;
	  }
	  break;
	default: badOption = true;
	}
      } // End of options
      if (badOption) { 
	std::cerr << argv[0] << ": Unknown option " << argv[i] 
		  << std::endl;
	return -1;
      }
    }
    else { // source or compatibility debug level 
      if (source.IsNull()) source = arg;
      else                 debugLevel = arg.Atoi();
    }
  }
  
  // --- Check if help was requested ---------------------------------
  if (help) { 
    ShowUsage(std::cout, argv[0]);
    return 1;
  }

  // --- Check if we have a source -----------------------------------
  if (source.IsNull()) { 
    std::cerr << "No source given" << std::endl;
    return -2;
  }

  // --- Initialize various things -----------------------------------
  AliFMDParameters::Instance()->Init(kFALSE,0);

  //This will only work for FDR 1 data. When newer data becomes
  //available the ! must be removed!
  Bool_t old = kTRUE;
  AliFMDParameters::Instance()->UseCompleteHeader(old);
  
  AliLog::EnableDebug(debugLevel > 0);
  AliLog::SetModuleDebugLevel("FMD", debugLevel);

  // --- Make our reader ---------------------------------------------
  fReader = AliRawReader::Create(source);
  if (!fReader) { 
    std::cerr << "Failed to make raw reader for source \"" << source 
	      << "\"" << std::endl;
    return -3;
  }
  return 0;
}

//_____________________________________________________________________ 
Int_t
AliFMDBaseDA::Runner::RunNumber() const
{ 
  if (!fReader) return -1;
  return fReader->GetRunNumber(); 
}

//_____________________________________________________________________ 
void
AliFMDBaseDA::Runner::Exec(AliFMDBaseDA& da)
{
  TStopwatch timer;
  timer.Start();

  da.SetSaveDiagnostics(fDiag || !fDiagFile.IsNull());
  if (!fDiagFile.IsNull()) da.SetDiagnosticsFilename(fDiagFile);

  da.Run(fReader);

  timer.Stop();
  timer.Print();

#ifdef ALI_AMORE
  try { 
    amore::da::AmoreDA myAmore(amore::da::AmoreDA::kSender);

    for (UShort_t det = 1; det <= 3; det++) {
      if (!da.HasSeenDetector(det)) continue;
      TObject* runNo = new TObject;
      runNo->SetUniqueID(fReader->GetRunNumber());
      myAmore.Send(Form("gainRunNoFMD%d", det), runNo);
    }
    
    TIter     next(&da.GetSummaries());
    TObject*  obj = 0;
    while ((obj = next())) 
      myAmore.Send(obj->GetName(), obj);
    
  }
  catch (exception& e) {
    cerr << "Failed to make AMORE instance: " << e.what() << endl;
  }
			       
#endif
}


  
  
  
//_____________________________________________________________________ 

//_____________________________________________________________________ 
//
// EOF
//
