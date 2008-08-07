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
  return Form("%s%sFMD%d%c", 
	      (full ? GetDetectorPath(det, full) : ""), 
	      (full ? "/" : ""), det, ring);
}
//_____________________________________________________________________
const char*
AliFMDBaseDA::GetDetectorPath(UShort_t det, 
			      Bool_t   full) const
{
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
  fDetectorArray(),
  fPulseSize(10),
  fPulseLength(10),
  fRequiredEvents(0),
  fCurrentEvent(0)
 {
  fDetectorArray.SetOwner();
  fConditionsFile.open("conditions.csv");
}
//_____________________________________________________________________
AliFMDBaseDA::AliFMDBaseDA(const AliFMDBaseDA & baseDA) : 
  TNamed(baseDA),
  fDiagnosticsFilename(baseDA.fDiagnosticsFilename.Data()),
  fOutputFile(),
  fConditionsFile(),
  fSaveHistograms(baseDA.fSaveHistograms),
  fDetectorArray(baseDA.fDetectorArray),
  fPulseSize(baseDA.fPulseSize),
  fPulseLength(baseDA.fPulseLength),
  fRequiredEvents(baseDA.fRequiredEvents),
  fCurrentEvent(baseDA.fCurrentEvent)
{
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
  TFile* diagFile = 0;
  if (fSaveHistograms)
    diagFile = TFile::Open(fDiagnosticsFilename.Data(),"RECREATE");

  
  
  
  InitContainer(diagFile);
  //Init();


  reader->Reset();
  
  AliFMDRawReader* fmdReader  = new AliFMDRawReader(reader,0);
  TClonesArray*    digitArray = new TClonesArray("AliFMDDigit",0);
  
  
  reader->NextEvent(); // Read Start-of-Run event
  reader->NextEvent(); // Read Start-of-Files event
  
  reader->NextEvent(); // Start-of-Data
  UInt_t eventType = reader->GetType();
  Bool_t SOD_read = kFALSE;
  if(eventType == AliRawEventHeaderBase::kStartOfData || 
     eventType == AliRawEventHeaderBase::kFormatError) { 
    
    WriteConditionsData(fmdReader);
    Init();
    SOD_read = kTRUE;
  }
  else {
    AliWarning("No SOD event detected!");
  }
    
  
  int lastProgress = 0;
  
  for(Int_t n =1;n <= GetRequiredEvents(); n++) {
    if(!reader->NextEvent()) continue;
    
    SetCurrentEvent(n);
    digitArray->Clear();
    fmdReader->ReadAdcs(digitArray);
    //std::cout<<"in event "<<*(reader->GetEventId())<<"   "<<n<<std::endl;
    //AliDebug(5, Form("In event # %d with %d entries", 
    //		     *(reader->GetEventId()), digitArray->GetEntriesFast()));
    
    for(Int_t i = 0; i<digitArray->GetEntriesFast();i++) {
      AliFMDDigit* digit = static_cast<AliFMDDigit*>(digitArray->At(i));
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
    std::cout << "FMD" << det << std::endl;
    UShort_t FirstRing = (det == 1 ? 1 : 0);
    for (UShort_t ir = FirstRing; ir < 2; ir++) {
      Char_t   ring = (ir == 0 ? 'O' : 'I');
      UShort_t nsec = (ir == 0 ? 40  : 20);
      UShort_t nstr = (ir == 0 ? 256 : 512);
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
  
  if(fSaveHistograms ) {
    
    Terminate(diagFile);
  
    AliInfo("Closing diagnostics file - please wait ...");
    // diagFile->Write();
    diagFile->Close();
    AliInfo("done");
  }
}
//_____________________________________________________________________

void AliFMDBaseDA::InitContainer(TDirectory* diagFile)
{
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
  AliFMDParameters* pars       = AliFMDParameters::Instance();
  fConditionsFile.write(Form("# %s \n",pars->GetConditionsShuttleID()),14);
  
  AliFMDCalibSampleRate* sampleRate = new AliFMDCalibSampleRate();
  AliFMDCalibStripRange* stripRange = new AliFMDCalibStripRange();
  
  fmdReader->ReadSODevent(sampleRate,stripRange,fPulseSize,fPulseLength);
  
  // Sample Rate
  /*
    UShort_t defSampleRate = 4;
    UShort_t sampleRateFromSOD;
    
  AliFMDCalibSampleRate* sampleRate = new AliFMDCalibSampleRate();

  UShort_t firstStrip = 0;
  UShort_t lastStrip  = 127;
  UShort_t firstStripSOD;
  UShort_t lastStripSOD;
  AliFMDCalibStripRange* stripRange = new AliFMDCalibStripRange();
  
  for(Int_t det=1;det<=3;det++) {
    UShort_t FirstRing = (det == 1 ? 1 : 0);
    for (UShort_t ir = FirstRing; ir < 2; ir++) {
      Char_t   ring = (ir == 0 ? 'O' : 'I');
      UShort_t nsec = (ir == 0 ? 40  : 20);
      for(UShort_t sec =0; sec < nsec;  sec++)  {
	sampleRateFromSOD = defSampleRate;
	sampleRate->Set(det,ring,sec,0,sampleRateFromSOD);
	firstStripSOD = firstStrip;
	lastStripSOD  = lastStrip;
	stripRange->Set(det,ring,sec,0,firstStripSOD,lastStripSOD);
	
      }
    }
  }
  */
  sampleRate->WriteToFile(fConditionsFile);
  stripRange->WriteToFile(fConditionsFile);
  //pars->SetSampleRate(sampleRate);
  //pars->SetStripRange(stripRange);

 
  // Zero Suppresion
  
  // Strip Range
  
  fConditionsFile.write("# Gain Events \n",15);
  
  for(UShort_t det=1; det<=3;det++) {
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
  
  

  // Gain Relevant stuff
  /*
  UShort_t defPulseSize = 32 ; 
  UShort_t defPulseLength = 100 ; 
  UShort_t pulseSizeFromSOD;
  UShort_t pulseLengthFromSOD;  
  
  fPulseSize.Reset(defPulseSize);
  fPulseLength.Reset(defPulseLength);
  
  for(UShort_t det=1;det<=3;det++)
    for(UShort_t iring=0;iring<=1;iring++)
      for(UShort_t board=0;board<=1;board++) {
	pulseSizeFromSOD = defPulseSize;
	pulseLengthFromSOD = defPulseLength;

	fPulseSize.AddAt(pulseSizeFromSOD,GetHalfringIndex(det,iring,board));
	fPulseLength.AddAt(pulseLengthFromSOD,GetHalfringIndex(det,iring,board));
      }
	  
  
  //  fConditionsFile     << defSampleRate   << ',' 
  //		      << timebins     <<"\n";
  */
  if(fConditionsFile.is_open()) {
    
    fConditionsFile.write("# EOF\n",6);
    fConditionsFile.close();
    
  }
  
}
//_____________________________________________________________________ 
Int_t AliFMDBaseDA::GetHalfringIndex(UShort_t det, Char_t ring, UShort_t board) {

  UShort_t iring  =  (ring == 'I' ? 1 : 0);
  
  Int_t index = (((det-1) << 2) | (iring << 1) | (board << 0));
  
  return index-2;
  
}


//_____________________________________________________________________ 
//
// EOF
//
