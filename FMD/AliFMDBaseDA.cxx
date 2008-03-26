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
//This is the implementation of the (virtual) base class for the FMD detector 
//algorithms(DA). It implements the creation of the relevant containers and handles the 
//loop over the raw data. The derived classes can control the parameters and action
//to be taken making this the base class for the Pedestal, Gain and Physics DA.
//

#include "AliFMDBaseDA.h"
#include "iostream"

#include "AliFMDRawReader.h"
#include "AliLog.h"
//_____________________________________________________________________
ClassImp(AliFMDBaseDA)

//_____________________________________________________________________
AliFMDBaseDA::AliFMDBaseDA() : TNamed(),
  fDiagnosticsFilename("diagnosticsHistograms.root"),
  fOutputFile(),
  fSaveHistograms(kFALSE),
  fDetectorArray(),
  fRequiredEvents(0),
  fCurrentEvent(0)
{
  fDetectorArray.SetOwner();
  
}
//_____________________________________________________________________
AliFMDBaseDA::AliFMDBaseDA(const AliFMDBaseDA & baseDA) : 
  TNamed(baseDA),
  fDiagnosticsFilename(baseDA.fDiagnosticsFilename),
  fOutputFile(),
  fSaveHistograms(baseDA.fSaveHistograms),
  fDetectorArray(baseDA.fDetectorArray),
  fRequiredEvents(baseDA.fRequiredEvents),
  fCurrentEvent(baseDA.fCurrentEvent)
{
  fDetectorArray.SetOwner();
  
}

//_____________________________________________________________________
AliFMDBaseDA::~AliFMDBaseDA() {

  //destructor
  
}

//_____________________________________________________________________
void AliFMDBaseDA::Run(AliRawReader* reader) {
  
  
  
  InitContainer();

  Init();

  TFile* diagFile = 0;
  if(fSaveHistograms)
    {
      diagFile = TFile::Open(fDiagnosticsFilename,"RECREATE");
      for(UShort_t det=1;det<=3;det++) {
	UShort_t FirstRing = (det == 1 ? 1 : 0);
	
	for (UShort_t ir = FirstRing; ir < 2; ir++) {
	  Char_t   ring = (ir == 0 ? 'O' : 'I');
	  UShort_t nsec = (ir == 0 ? 40  : 20);
	  UShort_t nstr = (ir == 0 ? 256 : 512);
	  
	  gDirectory->cd(Form("%s:/",fDiagnosticsFilename));
	  gDirectory->mkdir(Form("FMD%d%c",det,ring),Form("FMD%d%c",det,ring));
	  for(UShort_t sec =0; sec < nsec;  sec++)  {
	    gDirectory->cd(Form("%s:/FMD%d%c",fDiagnosticsFilename,det,ring));
	    gDirectory->mkdir(Form("sector_%d",sec));
	    for(UShort_t strip = 0; strip < nstr; strip++) {
	      gDirectory->cd(Form("%s:/FMD%d%c/sector_%d",fDiagnosticsFilename,det,ring,sec));
	      gDirectory->mkdir(Form("strip_%d",strip));
	      
	     }
	  }
	}
      }
      
    }
    
  reader->Reset();
  AliFMDRawReader* fmdReader = new AliFMDRawReader(reader,0);
  TClonesArray* digitArray   = new TClonesArray("AliFMDDigit",0);
  
  reader->NextEvent();
  reader->NextEvent();
  
  
  for(Int_t n =1;n <= GetRequiredEvents(); n++)
    {
      if(!reader->NextEvent()) 
	continue;
      
      SetCurrentEvent(*(reader->GetEventId()));
      
      digitArray->Clear();
      fmdReader->ReadAdcs(digitArray);
      
      // std::cout<<"In event # "<< *(reader->GetEventId()) << "\r"<<std::flush;
            
      for(Int_t i = 0; i<digitArray->GetEntries();i++) {
	AliFMDDigit* digit = static_cast<AliFMDDigit*>(digitArray->At(i));
	FillChannels(digit);
      }
      
      FinishEvent();
      
    }
  AliInfo(Form("Looped over %d events",GetCurrentEvent()));
  WriteHeaderToFile();
  
  for(UShort_t det=1;det<=3;det++) {
    UShort_t FirstRing = (det == 1 ? 1 : 0);
    for (UShort_t ir = FirstRing; ir < 2; ir++) {
      Char_t   ring = (ir == 0 ? 'O' : 'I');
      UShort_t nsec = (ir == 0 ? 40  : 20);
      UShort_t nstr = (ir == 0 ? 256 : 512);
      for(UShort_t sec =0; sec < nsec;  sec++)  {
  	for(UShort_t strip = 0; strip < nstr; strip++) {
	  Analyse(det,ring,sec,strip);
	  
  	}
      }
    }
  }
  if(fOutputFile.is_open()) {
   
    fOutputFile.write("# EOF\n",6);
    fOutputFile.close();
    
  }
  
  if(fSaveHistograms ) {
    AliInfo("Closing diagnostics file...please wait");
       
    delete diagFile;
  }
}
//_____________________________________________________________________

void AliFMDBaseDA::InitContainer(){

  TObjArray* detArray;
  TObjArray* ringArray;
  TObjArray* sectorArray;
    
  for(UShort_t det=1;det<=3;det++) {
    detArray = new TObjArray();
    detArray->SetOwner();
    fDetectorArray.AddAtAndExpand(detArray,det);
    UShort_t FirstRing = (det == 1 ? 1 : 0);
    for (UShort_t ir = FirstRing; ir < 2; ir++) {
      Char_t   ring = (ir == 0 ? 'O' : 'I');
      UShort_t nsec = (ir == 0 ? 40  : 20);
      UShort_t nstr = (ir == 0 ? 256 : 512);
      ringArray = new TObjArray();
      ringArray->SetOwner();
      detArray->AddAtAndExpand(ringArray,ir);
      for(UShort_t sec =0; sec < nsec;  sec++)  {
	sectorArray = new TObjArray();
	sectorArray->SetOwner();
	ringArray->AddAtAndExpand(sectorArray,sec);
	for(UShort_t strip = 0; strip < nstr; strip++) {
	  AddChannelContainer(sectorArray, det, ring, sec, strip);
	}
      }
    }
  }
}
//_____________________________________________________________________ 
//
// EOF
//
