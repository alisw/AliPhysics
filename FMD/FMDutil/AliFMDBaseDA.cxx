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
TString
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
  return TString::Format("%s%sFMD%d%c[%02d,%03d]", 
			 (full ? GetSectorPath(det, ring, sec, full).Data() : ""), 
			 (full ? "/" : ""), det, ring, sec, str);
}
//_____________________________________________________________________
TString
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
  return TString::Format("%s%sFMD%d%c[%02d]", 
			 (full ? GetRingPath(det, ring, full).Data() : ""), 
			 (full ? "/" : ""), det, ring, sec);
}
//_____________________________________________________________________
TString
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
  return TString::Format("%s%sFMD%d%c", 
			 (full ? GetDetectorPath(det, full).Data() : ""), 
			 (full ? "/" : ""), det, ring);
}
//_____________________________________________________________________
TString
AliFMDBaseDA::GetDetectorPath(UShort_t det, 
			      Bool_t   full) const
{
  // Get the strip path 
  // 
  // Parameters 
  //     det      Detector number
  //     ring     Rg identifier 
  //     sec      Sector number 
  //     str      Strip number
  //     full     If true, return full path 
  // 
  // Return 
  //     The path
  return TString::Format("%s%sFMD%d", 
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
  fSummaries(0),
  fAll(false)
{
  //Constructor
  for(Int_t i = 0; i< 3;i++) {
    fSeenDetectors[i] = false;
    fNEventsPerDetector[i] = 0;
  }
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
  fSummaries(0),
  fAll(baseDA.fAll)
{
  //Copy constructor
  for(Int_t i = 0; i< 3;i++) {
    fSeenDetectors[i] = baseDA.fSeenDetectors[0];
    fNEventsPerDetector[i] = baseDA.fNEventsPerDetector[i];
  }

  fDetectorArray.SetOwner();  
  fDetectorArray.ResetBit(TObject::kMustCleanup);
  fSummaries.ResetBit(TObject::kMustCleanup);
  fDetectorArray.ResetBit(TObject::kIsOnHeap);
  fSummaries.ResetBit(TObject::kIsOnHeap);
}


//_____________________________________________________________________
AliFMDBaseDA::~AliFMDBaseDA() 
{
  //destructor
  // fDetectorArray.ls();
  // fSummaries.ls();
  // fDetectorArray.SetOwner(false);
  // fSummaries.SetOwner(false);

}

//_____________________________________________________________________
Bool_t
AliFMDBaseDA::OpenFiles(Bool_t appendRun)
{
  fDetectorArray.SetOwner();
  if (!appendRun || fRunno == 0) {
    Rotate("conditions.csv", 3);
    fConditionsFile.open("conditions.csv");
  }
  else {
    fConditionsFile.open(Form("conditions_%09d.csv", 
			      fRunno));
  }
  if (!fConditionsFile) { 
    Error("OpenFiles", "Failed to open conditions file");
    return false;
  }
  return true;
}

//_____________________________________________________________________
Bool_t AliFMDBaseDA::HaveEnough(Int_t nEvents) const
{
  // if (!fAll) return nEvents > GetRequiredEvents();
  if (nEvents <= 1) return false;

  Bool_t ret = true; // Assume we have it 
  for (Int_t i = 0; i < 3; i++) { 
    if (!fSeenDetectors[i]) continue;
    if (Int_t(fNEventsPerDetector[i]) <= GetRequiredEvents()) ret = false;
  }
  return ret;
}
//_____________________________________________________________________
UShort_t AliFMDBaseDA::GetProgress(Int_t nEvents) const
{
  // if (!fAll) 
  //  return UShort_t((nEvents *100)/ GetRequiredEvents());

  if (nEvents <= 1) return 0;

  Int_t got = 0;
  Int_t total = 0;
  for (Int_t i = 0; i < 3; i++) {
    if (!fSeenDetectors[i]) continue;
    got   += fNEventsPerDetector[i];
    total += GetRequiredEvents();
  }
  return UShort_t(total > 0 ? (got * 100.) / total : 0);
}

//_____________________________________________________________________
Bool_t
AliFMDBaseDA::Run(AliRawReader* reader, Bool_t appendRun, Bool_t isBase) 
{
  //Run the FMD DA
  TFile* diagFile = 0;
  // if (fSaveHistograms)
  //  diagFile = TFile::Open(fDiagnosticsFilename.Data(),"RECREATE");

  reader->Reset();
  fRunno = reader->GetRunNumber();
  if (fRunno <= 0) { 
    TString dateRunNumber(gSystem->Getenv("DATE_RUN_NUMBER"));
    if (!dateRunNumber.IsNull()) fRunno = dateRunNumber.Atoll();
  }

  // Open our output files
  if (!OpenFiles(appendRun)) return false;
 
  // Create our reader 
  AliFMDRawReader* fmdReader  = new AliFMDRawReader(reader,0);
  TClonesArray*    digitArray = new TClonesArray("AliFMDDigit",0);
  
  // Look for SOD  
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
      sodread = kTRUE;
      break;
    }
  }
  if (isBase) 
    // IF this is the base DA, do thing more 
    return true;

  // Initialize everything 
  Init();
  // Initialize our containers 
  InitContainer(diagFile);
  if (AliLog::GetDebugLevel("FMD","") >= 3) { 
    fDetectorArray.ls();
  }
  
  if(!sodread) 
    AliWarning("No SOD event detected!");
  
  int lastProgress = 0;

  // Reset event counters 
  for(Int_t i = 0; i< 3;i++) fNEventsPerDetector[i] = 0;

  // Loop until we have enough 
  for(Int_t n = 1; !HaveEnough(n); n++) {
    AliDebugF(1,"Get the next event %d", n);
    if(!reader->NextEvent()) { n--; continue; }
    UInt_t eventType = reader->GetType();
    AliDebugF(3, "Event type is %d", eventType);
    if(eventType != AliRawEventHeaderBase::kPhysicsEvent) { n--; continue; }

    Bool_t seen[] = { false, false, false };
    SetCurrentEvent(n);
    digitArray->Clear();
    fmdReader->ReadAdcs(digitArray);
  
    for(Int_t i = 0; i<digitArray->GetEntriesFast();i++) {
      AliFMDDigit* digit = static_cast<AliFMDDigit*>(digitArray->At(i));
      UShort_t det = digit->Detector();
      fSeenDetectors[det-1] = true;
      seen[det-1]           = true;

      // Only fill if we do not have enough for this detector 
      if (Int_t(fNEventsPerDetector[det-1]) < GetRequiredEvents()) 
	FillChannels(digit);
    }
    for(Int_t i = 0; i< 3;i++) 
      if (seen[i]) (fNEventsPerDetector[i])++;
      
    FinishEvent();
    
    Int_t nReq = GetRequiredEvents();
    AliDebugF(5, "%9d: %6d/%6d %6d/%6d %6d/%6d", n, 
	      fNEventsPerDetector[0], nReq,
	      fNEventsPerDetector[1], nReq,
	      fNEventsPerDetector[2], nReq);

    int progress = GetProgress(n);
    if (progress <= lastProgress) continue;
    lastProgress = progress;
    std::cout << "Progress: " << lastProgress << " / 100 " << std::endl;
  }

  // Now at end of processing 
  AliInfoF("Looped over %d events (%d,%d,%d)",GetCurrentEvent(),
	   fNEventsPerDetector[0], 
	   fNEventsPerDetector[1], 
	   fNEventsPerDetector[2]);
  WriteHeaderToFile();
  
  // Analyse the data and make summaries 
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
  
  // If we have an output file, close it  
  if(fOutputFile.is_open()) {
    fOutputFile.write("# EOF\n",6);
    fOutputFile.close();
  }
  
  // Do final stuff on the diagnostics file 
  Terminate(diagFile);
    
  // Possibly write histograms to diagnostics file 
  if(fSaveHistograms ) {
    diagFile = TFile::Open(fDiagnosticsFilename.Data(),"RECREATE");
    fDetectorArray.Write("FMD", TObject::kSingleKey);
    fSummaries.Write();
    AliInfo("Closing diagnostics file - please wait ...");
    // diagFile->Write();
    diagFile->Close();
    AliInfo("done");
  }

  // All is good 
  return true;
}
//_____________________________________________________________________

void AliFMDBaseDA::InitContainer(TDirectory* diagFile)
{
  //Prepare container for diagnostics    
  TDirectory* savDir   = gDirectory;
  Bool_t owners = true;

  for(UShort_t det=1;det<=3;det++) {
    Array* detArray = new Array(det == 1 ? 1 : 2);
    detArray->SetOwner(owners);
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
      Array* ringArray = new Array(nsec);
      ringArray->SetOwner(owners);
      ringArray->SetName(GetRingPath(det, ring, false));
      ringArray->ResetBit(TObject::kMustCleanup);
      detArray->AddAtAndExpand(ringArray,ir);


      TDirectory* ringDir = 0;
      if (detDir) { 
	detDir->cd();
	ringDir = detDir->mkdir(GetRingPath(det,ring, kFALSE));
      }
      

      for(UShort_t sec =0; sec < nsec;  sec++)  {
	Array* sectorArray = new Array(nstr);
	sectorArray->SetOwner(owners);
	sectorArray->SetName(GetSectorPath(det, ring, sec, false));
	sectorArray->ResetBit(TObject::kMustCleanup);
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
	  Array* stripArray = new Array(0);
	  stripArray->SetOwner(owners);
	  stripArray->SetName(GetStripPath(det, ring, sec, strip, false));
	  stripArray->ResetBit(TObject::kMustCleanup);
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
  
  TString n = TString::Format("%sFMD%d%c", prefix, d, r);
  TString t = TString::Format("%s for FMD%d%c", title, d, r);
  TH2* ret = new TH2F(n, t, nX, -0.5, nX-0.5, nY, -0.5, nY-0.5);
  ret->SetXTitle("Sector #");
  ret->SetYTitle("Strip #");
  ret->SetDirectory(0);
  // if (!fSummaries) fSummaries = new Array;
  fSummaries.Add(ret);
  return ret;
}

//_____________________________________________________________________ 
AliFMDBaseDA::Array*
AliFMDBaseDA::GetDetectorArray(UShort_t det)
{
  if (det < 1 || det > 3) { 
    AliErrorF("Detector index %d out of bounds", det);
    return 0;
  }
  return static_cast<Array*>(fDetectorArray.At(det));
}
//_____________________________________________________________________ 
AliFMDBaseDA::Array*
AliFMDBaseDA::GetRingArray(UShort_t det, Char_t ring)
{
  Int_t idx = (ring == 'O' || ring == 'o' ? 0 : 1);
  Array* array = GetDetectorArray(det);
  if (!array) return 0;
  array = static_cast<Array*>(array->At(idx));
  if (!array) AliErrorF("No ring array for FMD%d%c (%d)", det, ring, idx);
  return array;
}
//_____________________________________________________________________ 
AliFMDBaseDA::Array*
AliFMDBaseDA::GetSectorArray(UShort_t det, Char_t ring, UShort_t sector)
{
  Array* array = GetRingArray(det,ring);
  if (!array) return 0;
  array = static_cast<Array*>(array->At(sector));
  if (!array) AliErrorF("No sector array for FMD%d%c[%02d]", det, ring, sector);
  return array;
}
//_____________________________________________________________________ 
AliFMDBaseDA::Array*
AliFMDBaseDA::GetStripArray(UShort_t det, Char_t ring, 
			    UShort_t sector, UShort_t strip)
{
  Array* array = GetSectorArray(det,ring,sector);
  if (!array) return 0;
  array = static_cast<Array*>(array->At(strip));
  if (!array) AliErrorF("No strip array for FMD%d%c[%02d,%03d]", 
			det, ring, sector, strip);
  return array;
}

//=====================================================================
AliFMDBaseDA::Runner::Runner()
  : fReader(0),
    fSource(""),
    fDiagFile(""), 
    fDiag(false),
    fAll(false),
    fFast(true),
    fUpload(true),
    fAppendRun(false),
    fOwnUpload(false)
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
    << "   -D,--debug=LEVEL         Set the debug level\n"
    << "   -A,--all                 Try to get data from all detectors (" << fAll << ")\n"
    << "   -F,--fast                Fast exit (" <<fFast <<")\n"
    << "   -U,--upload              Upload to FXS (" <<fUpload <<")\n"
    << "   -R,--append              Append run # to filename (" <<fAppendRun <<")\n"
    << "   -Z,--own                 Use custom upload (" << fOwnUpload << ")\n"
    << "\n"
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
AliFMDBaseDA::Runner::Init(int argc, char** argv, Bool_t reader)
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
	  fDiag = !fDiag;
	  if (ExtractValue(arg, val)) 
	    fDiagFile = val;
	}
	else if (arg.EqualTo("--all"))    fAll  = !fAll;
	else if (arg.EqualTo("--fast"))   fFast = !fFast;
	else if (arg.EqualTo("--upload")) fUpload = !fUpload;
	else if (arg.EqualTo("--append")) fAppendRun = !fAppendRun;
	else if (arg.EqualTo("--own"))    fOwnUpload = !fOwnUpload;
	else badOption = true;
      }
      else { // Short option 
	TString next(i < argc-1 ? argv[i+1] : "");
	switch (arg[1]) { 
	case 'h': help = true; break;
	case 'd': fDiag = !fDiag; 
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
	case 'A': fAll       = !fAll ; break ;
	case 'F': fFast      = !fFast ; break ;
	case 'U': fUpload    = !fUpload ; break ;
	case 'R': fAppendRun = !fAppendRun ; break ;
	case 'Z': fOwnUpload = !fOwnUpload ; break ;
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
      source = arg;
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
  fSource = source;

  // --- Initialize various things -----------------------------------
  AliFMDParameters::Instance()->Init(kFALSE,0);

  //This will only work for FDR 1 data. When newer data becomes
  //available the ! must be removed!
  Bool_t old = kTRUE;
  AliFMDParameters::Instance()->UseCompleteHeader(old);
  
  AliLog::EnableDebug(debugLevel > 0);
  AliLog::SetModuleDebugLevel("FMD", debugLevel);

  // --- Make our reader ---------------------------------------------
  if (!reader) return 0;
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
  Int_t run = -1;
  if (fReader) 
    run = fReader->GetRunNumber(); 
  if (run <= 0) {
    TString dateRunNumber(gSystem->Getenv("DATE_RUN_NUMBER"));
    if (!dateRunNumber.IsNull()) 
      run = dateRunNumber.Atoll();
    else 
      run = -1;
  }
  return run;
}

//_____________________________________________________________________ 
Bool_t
AliFMDBaseDA::Runner::Exec(AliFMDBaseDA& da)
{
  TStopwatch timer;
  timer.Start();

  da.SetSaveDiagnostics(fDiag || !fDiagFile.IsNull());
  da.SetTryAll(fAll);
  if (!fDiagFile.IsNull()) da.SetDiagnosticsFilename(fDiagFile);

  Bool_t ret = da.Run(fReader, fAppendRun);

  timer.Stop();
  timer.Print();

  return ret;
}

#if 0
AliFMDBaseDA::_Array::~_Array()
{
  // Printf("Deleting %s (%p)", this->GetName(), this);
  // SetOwner(false);
  // Clear();
}
#endif
  
  
//_____________________________________________________________________ 

//_____________________________________________________________________ 
//
// EOF
//
 
