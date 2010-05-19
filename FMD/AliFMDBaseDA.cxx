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
#include "iostream"
#include "AliFMDRawReader.h"
#include "AliFMDCalibSampleRate.h"
#include "AliFMDCalibStripRange.h"
#include "AliLog.h"
#include "AliRawEventHeaderBase.h"
#include <TDatime.h>
#include <TSystem.h>
#include <TH2F.h>

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
void AliFMDBaseDA::Run(AliRawReader* reader) 
{
  //Run the FMD DA
  TFile* diagFile = 0;
  if (fSaveHistograms)
    diagFile = TFile::Open(fDiagnosticsFilename.Data(),"RECREATE");

  
  
  
  
  reader->Reset();
  fRunno = reader->GetRunNumber();

  AliFMDRawReader* fmdReader  = new AliFMDRawReader(reader,0);
  TClonesArray*    digitArray = new TClonesArray("AliFMDDigit",0);
  
  Bool_t sodread = kFALSE;
  
  for(Int_t i=0;i<3;i++) {
    reader->NextEvent(); // Read Start-of-Run / Start-of-Files event
    
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
  
  if(!sodread) 
    AliWarning("No SOD event detected!");
  
  int lastProgress = 0;
  
  
  
  for(Int_t n =1;n <= GetRequiredEvents(); n++) {
    if(!reader->NextEvent()) continue;
    SetCurrentEvent(n);
    digitArray->Clear();
    fmdReader->ReadAdcs(digitArray);
    
    for(Int_t i = 0; i<digitArray->GetEntriesFast();i++) {
      AliFMDDigit* digit = static_cast<AliFMDDigit*>(digitArray->At(i));
      fSeenDetectors[digit->Detector()-1] = kTRUE;
      FillChannels(digit);
    }
    
   
    FinishEvent();
    
    int progress = int((n *100)/ GetRequiredEvents()) ;
    if (progress <= lastProgress) continue;
    lastProgress = progress;
    std::cout << "Progress: " << lastProgress << " / 100 " << std::endl;
    
  }
  
  AliInfo(Form("Looped over %d events",GetCurrentEvent()));
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
      if(fSaveHistograms)
      	diagFile->Flush();
      std::cout << "done" << std::endl;
    }
  }
  
  if(fOutputFile.is_open()) {
    fOutputFile.write("# EOF\n",6);
    fOutputFile.close();
  }
  
  Terminate(diagFile);
    
  if(fSaveHistograms ) {
    
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
  TObjArray* detArray;
  TObjArray* ringArray;
  TObjArray* sectorArray;
    
  TDirectory* savDir   = gDirectory;

  for(UShort_t det=1;det<=3;det++) {
    detArray = new TObjArray();
    detArray->SetOwner();
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
      ringArray = new TObjArray();
      ringArray->SetOwner();
      detArray->AddAtAndExpand(ringArray,ir);


      TDirectory* ringDir = 0;
      if (detDir) { 
	detDir->cd();
	ringDir = detDir->mkdir(GetRingPath(det,ring, kFALSE));
      }
      

      for(UShort_t sec =0; sec < nsec;  sec++)  {
	sectorArray = new TObjArray();
	sectorArray->SetOwner();
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
	  AddChannelContainer(sectorArray, det, ring, sec, strip);
	}
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
  Int_t nX = ((d == 1 || r == 'I' || r == 'i') ?  20 :  40);
  Int_t nY = ((d == 1 || r == 'I' || r == 'i') ? 512 : 256);
  
  TH2* ret = new TH2F(Form("%sFMD%d%c", prefix, d, r), 
		      Form("%s for FMD%d%c", title, d, r), 
		      nX, -0.5, nX-0.5, nY, -0.5, nY-0.5);
  ret->SetXTitle("Sector #");
  ret->SetYTitle("Strip #");

  // if (!fSummaries) fSummaries = new TObjArray;
  fSummaries.Add(ret);
  return ret;
}

//_____________________________________________________________________ 
//
// EOF
//
