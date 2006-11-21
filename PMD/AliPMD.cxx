/**************************************************************************
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
 **************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//
//                                                                           //
//  Photon Multiplicity Detector                                             //
//  This class contains the basic functions for the Photon Multiplicity      //
//  Detector. Functions specific to one particular geometry are              //
//  contained in the derived classes                                         //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliPMDClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:sub@vecdec.veccal.ernet.in">Subhasis Chattopadhyay</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TBRIK.h>
#include <TClonesArray.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTree.h>
#include <TVirtualMC.h>

#include "AliLog.h"
#include "AliLoader.h" 
#include "AliPMDLoader.h" 
#include "AliPMD.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliPMDDigitizer.h"
#include "AliPMDhit.h"
#include "AliPMDDDLRawData.h"
#include "AliPMDRawToSDigits.h"
  
ClassImp(AliPMD)
 
//_____________________________________________________________________________
AliPMD::AliPMD()
{
  //
  // Default constructor
  //
  fIshunt = 0;

}
 
//_____________________________________________________________________________
AliPMD::AliPMD(const char *name, const char *title)
  : AliDetector(name,title)
{
  //
  // Default constructor
  //

  // 
  // Allocate the array of hits
  fHits   = new TClonesArray("AliPMDhit",  405);
  gAlice->GetMCApp()->AddHitList(fHits);


  fIshunt =  0;
  
  fPar[0] = 1;
  fPar[1] = 1;
  fPar[2] = 0.8;
  fPar[3] = 0.02;
  fIn[0]  = 6;
  fIn[1]  = 20;
  fIn[2]  = 600;
  fIn[3]  = 27;
  fIn[4]  = 27;
  fGeo[0] = 0;
  fGeo[1] = 0.2;
  fGeo[2] = 4;
  fPadSize[0] = 0.8;
  fPadSize[1] = 1.0;
  fPadSize[2] = 1.2;
  fPadSize[3] = 1.5;
}

AliLoader* AliPMD::MakeLoader(const char* topfoldername)
{
  // Makes PMD Loader
 
  fLoader = new AliPMDLoader(GetName(),topfoldername);
 
  if (fLoader)
    {
      AliDebug(100,"Success");
    }
  else
    {
      AliError("Failure");
    }

  return fLoader;
}

AliPMD::~AliPMD()
{
  //
  // Destructor
  //
}

//_____________________________________________________________________________
void AliPMD::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a PMD hit
  //
  TClonesArray &lhits = *fHits;
  AliPMDhit *newcell, *curcell;
  //  printf("PMD++ Adding energy %f, prim %d, vol %d %d %d %d %d %d %d %d\n",
  // hits[3],gAlice->GetPrimary(track-1),vol[0],vol[1],vol[2],vol[3],
  // vol[4],vol[5],vol[6],vol[7]);

  newcell = new AliPMDhit(fIshunt, track, vol, hits);
  Int_t i;
  for (i=0; i<fNhits; i++) {
    //
    // See if this cell has already been hit
    curcell=(AliPMDhit*) lhits[i];
    if (*curcell==*newcell) {
//        printf("Cell with same numbers found\n") ; curcell->Print();
      *curcell = *curcell+*newcell;
//        printf("Cell after addition\n") ; curcell->Print();
      delete newcell;
      return;
    }
  }
  new(lhits[fNhits++]) AliPMDhit(newcell);
  delete newcell;
}
 
//_____________________________________________________________________________
void AliPMD::BuildGeometry()
{
  //
  // Build simple ROOT TNode geometry for event display
  //

  TNode *node, *top;
  const int kColorPMD  = kRed;

  //
  top=gAlice->GetGeometry()->GetNode("alice");

  // PMD
  new TBRIK("S_PMD","PMD box","void",300,300,5);
  top->cd();
  node = new TNode("PMD","PMD","S_PMD",0,0,-600,"");
  node->SetLineColor(kColorPMD);
  fNodes->Add(node);
}

//_____________________________________________________________________________
void AliPMD::SetPAR(Float_t p1, Float_t p2, Float_t p3,Float_t p4)
{
  //
  // Set PMD parameters
  //
  fPar[0] = p1;
  fPar[1] = p2;
  fPar[2] = p3;
  fPar[3] = p4;
}
 
//_____________________________________________________________________________
void AliPMD::SetIN(Float_t p1, Float_t p2, Float_t p3,Float_t p4,Float_t p5)
{
  //
  // Set PMD parameters
  //
  fIn[0] = p1;
  fIn[1] = p2;
  fIn[2] = p3;
  fIn[3] = p4;
  fIn[4] = p5;
}
 
//_____________________________________________________________________________
void AliPMD::SetGEO(Float_t p1, Float_t p2, Float_t p3)
{
  //
  // Set geometry parameters
  //
  fGeo[0] = p1;
  fGeo[1] = p2;
  fGeo[2] = p3;
}
 
//_____________________________________________________________________________
void AliPMD::SetPadSize(Float_t p1, Float_t p2, Float_t p3,Float_t p4)
{
  //
  // Set pad size
  //
  fPadSize[0] = p1;
  fPadSize[1] = p2;
  fPadSize[2] = p3;
  fPadSize[3] = p4;
}
 
//_____________________________________________________________________________
void AliPMD::StepManager()
{
  //
  // Called at every step in PMD
  //
}

void AliPMD::MakeBranch(Option_t* option)
{
    // Create Tree branches for the PMD
    
    const char *cH = strstr(option,"H");
    if (cH && fLoader->TreeH() && (fHits == 0x0))
      fHits   = new TClonesArray("AliPMDhit",  405);
    
    AliDetector::MakeBranch(option);
}


void AliPMD::SetTreeAddress()
{
  // Set branch address

    if (fLoader->TreeH() && fHits==0x0)
      fHits   = new TClonesArray("AliPMDhit",  405);
      
    AliDetector::SetTreeAddress();
}

//____________________________________________________________________________
void AliPMD::Hits2SDigits()  
{ 
// create summable digits

  AliRunLoader* runLoader = fLoader->GetRunLoader(); 
  AliPMDDigitizer* pmdDigitizer = new AliPMDDigitizer;
  pmdDigitizer->OpengAliceFile(fLoader->GetRunLoader()->GetFileName().Data(),
			       "HS");
  pmdDigitizer->SetZPosition(361.5);

  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    pmdDigitizer->Hits2SDigits(iEvent);
  }
  fLoader->UnloadHits();
  fLoader->UnloadSDigits();
  delete pmdDigitizer;
}
//____________________________________________________________________________
void AliPMD::SDigits2Digits()  
{ 
  // creates sdigits to digits
}
//____________________________________________________________________________
void AliPMD::Hits2Digits()  
{ 
// create digits

  AliRunLoader* runLoader = fLoader->GetRunLoader(); 
  AliPMDDigitizer* pmdDigitizer = new AliPMDDigitizer;
  pmdDigitizer->OpengAliceFile(fLoader->GetRunLoader()->GetFileName().Data(),
			       "HD");
  pmdDigitizer->SetZPosition(361.5);

  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    pmdDigitizer->Hits2Digits(iEvent);
  }
  fLoader->UnloadHits();
  fLoader->UnloadDigits();
  delete pmdDigitizer;

}
// ---------------------------------------------------------------------------
AliDigitizer* AliPMD::CreateDigitizer(AliRunDigitizer* manager) const
{ 
  return new AliPMDDigitizer(manager);
}
// ---------------------------------------------------------------------------
void AliPMD::Digits2Raw()
{ 
// convert digits of the current event to raw data

  fLoader->LoadDigits();
  TTree* digits = fLoader->TreeD();
  if (!digits) {
    AliError("No digits tree");
    return;
  }

  AliPMDDDLRawData rawWriter;
  rawWriter.WritePMDRawData(digits);

  fLoader->UnloadDigits();
}

Bool_t AliPMD::Raw2SDigits(AliRawReader *rawReader)
{
  // converts raw to sdigits
  AliRunLoader* runLoader = fLoader->GetRunLoader(); 
  //runLoader->GetEvent(ievt);

  AliPMDRawToSDigits pmdr2sd;
  pmdr2sd.Raw2SDigits(runLoader, rawReader);
  fLoader->UnloadSDigits();
  return kTRUE;
}

