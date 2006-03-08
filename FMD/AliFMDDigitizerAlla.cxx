 /**************************************************************************
 * Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
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

 //////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Forward Multiplicity Detector based on Silicon plates                    //
//  This class contains the procedures simulation ADC  signal for            //
//  the Forward Multiplicity detector  : hits -> digits                      //
//  ADC signal consists                                                      //
//   - number of detector;                                                   //
//   - number of ring;                                                       //
//   - number of sector;                                                     //
//   - ADC signal in this channel                                            //
//                                                                           //
 //////////////////////////////////////////////////////////////////////////////

#include <TTree.h> 
#include <TVector.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TRandom.h>
#include "AliLog.h"
#include "AliFMDDigitizerAlla.h"
#include "AliFMD.h"
#include "AliFMDHit.h"
#include "AliFMDDigit.h"
#include "AliRunDigitizer.h"
#include "AliRun.h"
#include "AliLoader.h"
#include "AliRunLoader.h"
#include <stdlib.h>
#include <Riostream.h>

ClassImp(AliFMDDigitizerAlla)
#if 0
  ;
#endif

//___________________________________________
AliFMDDigitizerAlla::AliFMDDigitizerAlla()  
  : AliDigitizer()
{
  // Default ctor - don't use it
}

//___________________________________________
AliFMDDigitizerAlla::AliFMDDigitizerAlla(AliRunDigitizer* manager) 
  : AliDigitizer(manager) 
{
  // ctor which should be used
  //  fDebug =0;
  AliDebug(1," processed");
}

//------------------------------------------------------------------------
AliFMDDigitizerAlla::~AliFMDDigitizerAlla()
{
  // Destructor
}

 //------------------------------------------------------------------------
Bool_t AliFMDDigitizerAlla::Init()
{
  // Initialization
  // cout<<"AliFMDDigitizerAlla::Init"<<endl;
  return kTRUE;
}
 

//---------------------------------------------------------------------

void AliFMDDigitizerAlla::Exec(Option_t * /*option*/)
{

  // Conver hits to digits:
  //  - number of detector;
  //  - number of ring;
  //  - number of sector;
  //  - ADC signal in this channel
  AliDebug(1," start...");
  Float_t de[3][2][50][520];
  for (Int_t i=0; i<3; i++)
    for(Int_t j=0; j<2; j++)
      for(Int_t k=0; k < 50; k++)
	for (Int_t l=0; l < 520; l++)
	  de[i][j][k][l]=0;
  Int_t nstrips[]  = { 512, 256 };
  Int_t nsectors[] = {  20,  40 };
  
  AliRunLoader* outRL = 
    AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  AliLoader* pOutFMD = outRL->GetLoader("FMDLoader");

  AliRunLoader* inRL = 
    AliRunLoader::GetRunLoader(fManager->GetInputFolderName(0));

  if (!inRL) {
    AliError("Can not find Run Loader for input stream 0");
    return;
  }
  if (!inRL->GetAliRun()) inRL->LoadgAlice();

  AliFMD * fFMD = static_cast<AliFMD*>(inRL->GetAliRun()->GetDetector("FMD"));
  if (!fFMD) {
    AliError("Can not get FMD from gAlice");
    return;
  }

  // Loop over files to digitize
  Int_t nFiles = GetManager()->GetNinputs();
  for (Int_t inputFile=0; inputFile < nFiles; inputFile++)  {
    AliDebug(1,Form("Digitizing event number %d",
		    fManager->GetOutputEventNr()));
    if (fFMD) {
       AliRunLoader* inRL = 
	 AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inputFile));
       AliLoader* pInFMD = inRL->GetLoader("FMDLoader");
       pInFMD->LoadHits("READ");
      
      
      TTree* tH = pInFMD->TreeH();
      if (!tH) {
	pInFMD->LoadHits("read");
	tH = pInFMD->TreeH();
      }
      TBranch* brHits = tH->GetBranch("FMD");
      if (brHits)  fFMD->SetHitsAddressBranch(brHits);
      else         AliFatal("Branch FMD hit not found");
      TClonesArray *fFMDhits = fFMD->Hits ();
      
      Int_t ntracks    = tH->GetEntries();
      for (Int_t track = 0; track < ntracks; track++) {
	brHits->GetEntry(track);
	Int_t nhits = fFMDhits->GetEntries ();
	for (Int_t hit = 0; hit < nhits; hit++) {
	  AliFMDHit* fmdHit = 
	    static_cast<AliFMDHit*>(fFMDhits->UncheckedAt(hit));
	  Int_t detector = fmdHit->Detector();
	  Int_t iring    = fmdHit->Ring() == 'I' ? 0 : 1;
	  Int_t sector   = fmdHit->Sector();
	  Int_t strip    = fmdHit->Strip();
	  de[detector-1][iring][sector][strip] += fmdHit->Edep();
	}          //hit loop
      } //track loop
    } //if FMD
  }

 
  // Put noise and make ADC signal
   Float_t mipI = 1.664 * 0.04 * 2.33 / 22400;     // = 6.923e-6;
   for (Int_t detector = 1; detector <= 3; detector++){
     for (Int_t iring = 0; iring < 2; iring++) {
       if (detector == 1 && iring == 1) continue;
       char ring = (iring == 0 ? 'I' : 'O');
       for (Int_t sector = 0; sector < nsectors[iring]; sector++) {
	 for (Int_t strip = 0; strip < nstrips[iring]; strip++) {
	   Int_t signal   = Int_t(de[detector-1][iring][sector][strip] / mipI);
	   Int_t pedestal = Int_t(gRandom->Gaus(500,250));
	   Int_t charge   = signal + pedestal;
	   if(charge <= 500) charge = 500; 
	   //dynamic range from MIP(0.155MeV) to 30MIP(4.65MeV)
	   //1024 ADC channels 
	   Float_t channelWidth = (22400 * 50) / 1024;
	   Int_t   adc          = Int_t(charge / channelWidth);
	   if (adc > 1023)  adc = 1023; 
	   fFMD->AddDigitByFields(detector, ring, sector, strip, adc);
	 } //strip
       } //sector
     } //iring
   } // detector
   
   pOutFMD->LoadDigits("update");
   TTree* treeD = pOutFMD->TreeD();
   if (!treeD) {
     pOutFMD->MakeTree("D");
     treeD = pOutFMD->TreeD();
   }

   treeD->Reset();
   TClonesArray* digits = fFMD->Digits();
   fFMD->MakeBranchInTree(treeD, fFMD->GetName(), &(digits), 4000, 0);
   if (!treeD->GetBranch("FMD")) AliFatal("No branch for FMD digits");
   treeD->Fill();
   pOutFMD->WriteDigits("OVERWRITE");
   pOutFMD->UnloadHits();
   pOutFMD->UnloadDigits();
   fFMD->ResetDigits();
}

//____________________________________________________________________
//
// EOF
//
 




