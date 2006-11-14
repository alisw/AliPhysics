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
//                                                                           //
//  START (T-Zero) Detector                                            //
//  This class contains the base procedures for the START     //
//  detector                                                                 //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliSTARTClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:Alla.Maevskaia@cern.ch">Alla Maevskaia</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>

#include <TFile.h>
#include <TGeometry.h>
#include <TMath.h>
#include <TNode.h>
#include <TParticle.h>
#include <TRandom.h>
#include <TTUBE.h>
#include <TVirtualMC.h>
#include <AliESD.h>

#include "AliLog.h"
#include "AliMC.h"
#include "AliLoader.h"
#include "AliRun.h"
#include "TClonesArray.h"
#include "AliSTART.h"
#include "AliSTARTLoader.h"
#include "AliSTARTdigit.h"
#include "AliSTARThit.h"
#include "AliSTARTDigitizer.h"
#include "AliSTARTRawData.h"
#include "AliSTARTRecPoint.h"
#include "AliLog.h"

ClassImp(AliSTART)

  //static  AliSTARTdigit *digits; 

//_____________________________________________________________________________
AliSTART::AliSTART()
  : AliDetector(), fIdSens(0), fDigits(NULL), fRecPoints(NULL)
{
  //
  // Default constructor for class AliSTART
  //
  fIshunt   = 1;
  fHits     = 0;
  fDigits   = 0;
  fRecPoints = 0;
}
 
//_____________________________________________________________________________
AliSTART::AliSTART(const char *name, const char *title)
  : AliDetector(name,title), fIdSens(0), fDigits(new AliSTARTdigit()), fRecPoints(new AliSTARTRecPoint())
{
  //
  // Standard constructor for START Detector
  //

  
  //
  // Initialise Hit array
  fHits       = new TClonesArray("AliSTARThit",  405);
  gAlice->GetMCApp()->AddHitList(fHits);
  //  fDigits    = new AliSTARTdigit();
  //  fRecPoints = new AliSTARTRecPoint();
  fIshunt     =  1;
  //  fIdSens   =  0;
  //PH  SetMarkerColor(kRed);
}

//_____________________________________________________________________________
AliSTART::~AliSTART() {
  
  if (fHits) {
    fHits->Delete();
    delete fHits;
  }
  /*
  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
    cout<<" delete fDigits; "<<endl;
  }
  if (fRecPoints) {
   fRecPoints ->Delete();
    delete fRecPoints;
    cout<<" delete fRecPoints; "<<endl;
  }
  */ 
}

//_____________________________________________________________________________
void AliSTART::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a START hit
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliSTARThit(fIshunt,track,vol,hits);
}


//_____________________________________________________________________________

void AliSTART::AddDigit(Int_t besttimeright, Int_t besttimeleft, Int_t meantime, 
			Int_t timediff, Int_t sumMult,
			TArrayI *time, TArrayI *adc, TArrayI *timeAmp, TArrayI *adcAmp)
{
  
  //  Add a START digit to the list.
 //
  
  if (!fDigits) {
    fDigits = new AliSTARTdigit();
  }
  fDigits-> SetTimeBestRight(besttimeright);
  fDigits->SetTimeBestLeft(besttimeleft);
  fDigits-> SetMeanTime(meantime);
  fDigits-> SetDiffTime(timediff);
  fDigits-> SetSumMult(sumMult);
  fDigits->SetTime(*time);
  fDigits->SetTimeAmp(*timeAmp);
  fDigits->SetADC(*adc);
  fDigits->SetADCAmp(*adcAmp);
}


//_____________________________________________________________________________
void AliSTART::BuildGeometry()
{
  //
  // Build simple ROOT TNode geometry for event display
  //
  TNode *node, *top;
  const int kColorSTART  = 19;

  top=gAlice->GetGeometry()->GetNode("alice");

  // START define the different volumes
  new TRotMatrix("rotx999","rot999",  90,0,90,90,180,0);

  new TTUBE("S_0ST1","START  volume 1","void",5.,10.7,5.3);
  top->cd();
  node = new TNode("0ST1","0ST01","S_0ST1",0,0,-69.7,"");
  node->SetLineColor(kColorSTART);
  fNodes->Add(node);

  new TTUBE("S_0ST2","START volume 2","void",5.,10.7,5.3);
  top->cd();
  node = new TNode("0ST2","0ST2","S_0ST2",0,0,350,"rotx999");
  node->SetLineColor(kColorSTART);
  fNodes->Add(node);
}
 
//_____________________________________________________________________________
Int_t AliSTART::DistanceToPrimitive(Int_t /*px*/, Int_t /*py*/)
{
  //
  // Calculate the distance from the mouse to the START on the screen
  // Dummy routine
  //
  return 9999;
}
 
//-------------------------------------------------------------------------
void AliSTART::Init()
{
  //
  // Initialis the START after it has been built
  Int_t i;
  //
  if(AliLog::GetGlobalDebugLevel()>0) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" START_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    //
    // Here the START initialisation code (if any!)
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }
}

//---------------------------------------------------------------------------
void AliSTART::MakeBranch(Option_t* option)
{
  //
// Create Tree branches for the START.

 // Options:
  //
  //    H          Make a branch of TClonesArray of AliSTARTHit's
  //    D          Make a branch of TClonesArray of AliSTARTDigit's
  //
  //    R         Make a branch of  AliSTARTRecPoints
  //
  char branchname[20];
  sprintf(branchname,"%s",GetName());

  const char *cH = strstr(option,"H");
  const char *cD = strstr(option,"D");
  const char *cR = strstr(option,"R");

    if (cH && fLoader->TreeH())
  {
     if (fHits == 0x0) fHits  = new TClonesArray("AliSTARThit",  405);
     AliDetector::MakeBranch(option);
  } 
    
    
  if (cD && fLoader->TreeD())
    {
      if (fDigits == 0x0) fDigits  = new AliSTARTdigit();
      //     MakeBranchInTree(fLoader->TreeD(), branchname,
      //		       &fDigits, 405, 0);
      fLoader->TreeD()->Branch(branchname,"AliSTARTdigit",&fDigits,405,1);
      //   fLoader->TreeD()->Print();
    } 
  if (cR && fLoader->TreeR())
    {
      if (fRecPoints == 0x0) fRecPoints  = new AliSTARTRecPoint();
      MakeBranchInTree(fLoader->TreeR(), branchname,
		       &fRecPoints, 405, 0);
    } 
  
}    

//_____________________________________________________________________________
void AliSTART::ResetHits()
{
  AliDetector::ResetHits();
  
}
//____________________________________________________________________
void AliSTART::ResetDigits()
{
  //
  // Reset number of digits and the digits array for this detector
  //
  if (fDigits) fDigits->Clear();
}

//_____________________________________________________________________________
void AliSTART::SetTreeAddress()
{

  TTree    *treeH;
  treeH = TreeH();
  
  if (treeH)
    {
      if (fHits == 0x0) fHits  = new TClonesArray("AliSTARThit",  405);
    }
    
  AliDetector::SetTreeAddress();
  TTree *treeD = fLoader->TreeD();
  if (treeD) {
    if (fDigits == 0x0)  fDigits  = new AliSTARTdigit();
    TBranch* branch = treeD->GetBranch ("START");
    if (branch) branch->SetAddress(&fDigits);
  }

  TTree *treeR = fLoader->TreeR();
  if (treeR) {
    if (fRecPoints == 0x0) fRecPoints  = new  AliSTARTRecPoint()  ;
    TBranch* branch = treeR->GetBranch ("START");
    if (branch) branch->SetAddress(&fRecPoints);
  }
 
}


//_____________________________________________________________________________
void AliSTART::MakeBranchInTreeD(TTree *treeD, const char *file)
{
    //
    // Create TreeD branches for the FMD
    //
    const Int_t kBufferSize = 4000;
    char branchname[20];
    sprintf(branchname,"%s",GetName());
    if(treeD)
     {
       MakeBranchInTree(treeD,  branchname,&fDigits, kBufferSize, file);
     }
}

//_____________________________________________________________________________
AliDigitizer* AliSTART::CreateDigitizer(AliRunDigitizer* manager) const
{
  return new AliSTARTDigitizer(manager);
}
//____________________________________________________________________________
void AliSTART::Digits2Raw()
{
//
// Starting from the START digits, writes the Raw Data objects
//
//  AliSTARTLoader* pStartLoader = (AliSTARTLoader*)fLoader;
  fLoader ->LoadDigits("read");
  TTree* treeD = fLoader->TreeD();
  if (!treeD) {
    AliError("no digits tree");
    return;
  }
  if (fDigits == 0x0)  fDigits  = new AliSTARTdigit();
  
  TBranch *branch = treeD->GetBranch("START");
  if (branch) {
    branch->SetAddress(&fDigits);
  }else{
    AliError("Branch START DIGIT not found");
    exit(111);
  } 
  AliSTARTRawData rawWriter;
  rawWriter.SetVerbose(0);
  
  AliDebug(2,Form(" Formatting raw data for START "));
  branch->GetEntry(0);
  //  rawWriter.RawDataSTART(treeD->GetBranch("START"));
  rawWriter.RawDataSTART(fDigits);
  
  
  fLoader->UnloadDigits();
  
}
