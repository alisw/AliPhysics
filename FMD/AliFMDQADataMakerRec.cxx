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
// --- ROOT system ---
#include <iostream>
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH1I.h> 

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliFMDQADataMakerRec.h"
#include "AliFMDDigit.h"
#include "AliFMDRecPoint.h"
#include "AliQAChecker.h"
#include "AliESDFMD.h"
#include "AliFMDParameters.h"

//_____________________________________________________________________
// This is the class that collects the QA data for the FMD during
// reconstruction.  
//
// The following data types are picked up:
// - digits
// - rec points
// - esd data
// The following data types are not supported (yet):
// - raws
// Author : Hans Hjersing Dalsgaard, hans.dalsgaard@cern.ch
//_____________________________________________________________________

ClassImp(AliFMDQADataMakerRec)
#if 0
; // For Emacs - do not delete!
#endif
           
//_____________________________________________________________________
AliFMDQADataMakerRec::AliFMDQADataMakerRec() : 
  AliQADataMakerRec(AliQA::GetDetName(AliQA::kFMD), 
		    "FMD Quality Assurance Data Maker")
{
  // ctor
}

//_____________________________________________________________________
AliFMDQADataMakerRec::AliFMDQADataMakerRec(const AliFMDQADataMakerRec& /*qadm*/) 
  : AliQADataMakerRec()
{
  // copy ctor 
  // Parameters: 
  //    qadm    Object to copy from
  
}


//_____________________________________________________________________ 

void 
AliFMDQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX task, 
					 TObjArray * list)
{
  // Detector specific actions at end of cycle
  // do the QA checking
  AliLog::Message(5,"FMD: end of detector cycle",
		  "AliFMDQADataMakerRec","AliFMDQADataMakerRec",
		  "AliFMDQADataMakerRec::EndOfDetectorCycle",
		  "AliFMDQADataMakerRec.cxx",95);
  AliQAChecker::Instance()->Run(AliQA::kFMD, task, list);
}

//_____________________________________________________________________ 
void AliFMDQADataMakerRec::InitESDs()
{
  // create Digits histograms in Digits subdir
  TH1F* hEnergyOfESD = new TH1F("hEnergyOfESD","Energy distribution",100,0,3);
  hEnergyOfESD->SetXTitle("Edep/Emip");
  hEnergyOfESD->SetYTitle("Counts");
  Add2ESDsList(hEnergyOfESD, 0);
    
}

//_____________________________________________________________________ 
void AliFMDQADataMakerRec::InitDigits()
{
  // create Digits histograms in Digits subdir
  TH1I* hADCCounts      = new TH1I("hADCCounts","Dist of ADC counts",
				   1024,0,1024);
  hADCCounts->SetXTitle("ADC counts");
  hADCCounts->SetYTitle("");
  Add2DigitsList(hADCCounts, 0);
  
}

//_____________________________________________________________________ 
void AliFMDQADataMakerRec::InitRecPoints()
{
  TH1F* hEnergyOfRecpoints = new TH1F("hEnergyOfRecpoints",
				      "Energy Distribution",100,0,3);
  hEnergyOfRecpoints->SetXTitle("Edep/Emip");
  hEnergyOfRecpoints->SetYTitle("");
  Add2RecPointsList(hEnergyOfRecpoints,0);
}

//_____________________________________________________________________ 
void AliFMDQADataMakerRec::InitRaws()
{
  
}

//_____________________________________________________________________
void AliFMDQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  if(!esd) {
    AliError("FMD ESD object not found!!") ; 
    return;
  }
  AliESDFMD* fmd = esd->GetFMDData();
  if (!fmd) return;
  
  for(UShort_t det=1;det<=3;det++) {
    for (UShort_t ir = 0; ir < 2; ir++) {
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ? 20  : 40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      for(UShort_t sec =0; sec < nsec;  sec++)  {
	for(UShort_t strip = 0; strip < nstr; strip++) {
	  Float_t mult = fmd->Multiplicity(det,ring,sec,strip);
	  if(mult == AliESDFMD::kInvalidMult) continue;
	  
	  GetESDsData(0)->Fill(mult);
	}
      }
    }
  }
}

//_____________________________________________________________________
void AliFMDQADataMakerRec::MakeDigits(TClonesArray * digits)
{
  // makes data from Digits  
  if(!digits)  {
    AliError("FMD Digit object not found!!") ;
    return;
  }
  for(Int_t i=0;i<digits->GetEntries();i++) {
    //Raw ADC counts
    AliFMDDigit* digit = static_cast<AliFMDDigit*>(digits->At(i));
    GetDigitsData(0)->Fill(digit->Counts());
  }
}

//_____________________________________________________________________
void AliFMDQADataMakerRec::MakeDigits(TTree * digitTree)
{
  
  TClonesArray* digits = new TClonesArray("AliFMDDigit", 1000); 
  TBranch*      branch = digitTree->GetBranch("FMD");
  if (!branch) {
    AliWarning("FMD branch in Digit Tree not found") ; 
    return;
  } 
  
  branch->SetAddress(&digits);
  branch->GetEntry(0); 
  MakeDigits(digits); 
}

//_____________________________________________________________________
void AliFMDQADataMakerRec::MakeRaws(AliRawReader* /*rawReader*/)
{
  
}

//_____________________________________________________________________
void AliFMDQADataMakerRec::MakeRecPoints(TTree* clustersTree)
{
  // makes data from RecPoints
  TBranch *fmdbranch = clustersTree->GetBranch("FMD");
  if (!fmdbranch) { 
    AliError("can't get the branch with the FMD recpoints !");
    return;
  }
  
  TClonesArray * fmdrecpoints = new TClonesArray("AliFMDRecPoint", 1000);
  fmdbranch->SetAddress(&fmdrecpoints);
  fmdbranch->GetEntry(0);
    
  TIter next(fmdrecpoints) ; 
  AliFMDRecPoint * rp ; 
  while ((rp = static_cast<AliFMDRecPoint*>(next()))) {
    GetRecPointsData(0)->Fill(rp->Particles()) ;
  }
  fmdrecpoints->Delete();
  delete fmdrecpoints;

}

//_____________________________________________________________________ 
void AliFMDQADataMakerRec::StartOfDetectorCycle()
{
  // What 
  // to 
  // do?
}
//_____________________________________________________________________ 
//
// EOF
//
