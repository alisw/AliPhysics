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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  FIT ( Fast Interaction Trigger ) Detector                                            //
//  This class contains the base procedures for the FIT    //
//  detector                                                                 //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliFITClass.gif">
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

#include "TClonesArray.h"
#include "TString.h"

#include "AliLoader.h"
#include "AliLog.h"
#include "AliLog.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliFIT.h"
#include "AliFITDigitizer.h"
#include "AliFITDigit.h"
#include "AliFITHits.h"
#include "AliFITRawData.h"
#include "AliFITRawReader.h"

ClassImp(AliFIT)

  //static  AliFITdigit *digits; 

//_____________________________________________________________________________
AliFIT::AliFIT()
  : AliDetector(), 
  fIdSens(0) 
    /* , fDigits(NULL)*/
{
  //
  // Default constructor for class AliFIT
  //
  fIshunt   = 1;
  fHits     = 0;
  fDigits = new TClonesArray("AliFITDigit",100); // from AliDetector
  fNdigits   = 0;
}
 
//_____________________________________________________________________________
AliFIT::AliFIT(const char *name, const char *title)
  : AliDetector(name,title), 
    fIdSens(0) 
    /* , fDigits(new AliFITDigit())*/
{
  //
  // Standard constructor for T0 Detector
  //
  //
  // Initialise Hit array
  AliMC* mc = gAlice->GetMCApp();
  if( mc && mc->GetHitLists() ) {
   fHits = new TClonesArray("AliFITHits",100); // from AliDetector
   mc->AddHitList(fHits);  
  }
 
   fNdigits   = 0;
 fDigits = new TClonesArray("AliFITDigit",100); // from AliDetector
  fIshunt     =  1;
  //  fIdSens   =  0;
  //PH  SetMarkerColor(kRed);
}

//_____________________________________________________________________________
AliFIT::~AliFIT() {
  
  //destructor
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
void AliFIT::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a FIT hit
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliFITHits(fIshunt,track,vol,hits);
}

//____________________________________________________________________________
void AliFIT::AddDigit(Int_t npmt, 
			       Int_t timeCFD, Int_t timeLED, Int_t timeQT0,
			       Int_t timeQT1, Int_t *labels) 
 { 
 
// Adds Digit 
  TClonesArray &ldigits = *fDigits; 
   new(ldigits[fNdigits++]) AliFITDigit( npmt,timeCFD, timeLED, timeQT0, timeQT1, labels);
}


//-------------------------------------------------------------------------
void AliFIT::Init()
{
  //
  // Initialis the T0 after it has been built
  Int_t i;
  //
  //  if(AliLog::GetGlobalDebugLevel()>0) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" FIT_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    //
    // Here the T0 initialisation code (if any!)
    for(i=0;i<80;i++) printf("*");
    printf("\n");
    // }
}

//---------------------------------------------------------------------------
void AliFIT::MakeBranch(Option_t* option)
{
  //
// Create Tree branches for the T0.

 // Options:
  //
  //    H          Make a branch of TClonesArray of AliT0Hit's
  //    D          Make a branch of TClonesArray of AliT0Digit's
  //
  //    R         Make a branch of  AliT0RecPointUps
  //
  //  char branchname[20];
  // sprintf(branchname,"%s",GetName());
  //  strncpy(branchname, GetName(), 20);
  TString branchname = Form("%s", GetName());

  const char *cH = strstr(option,"H");
  const char *cD = strstr(option,"D");
  const char *cS = strstr(option,"S");

    if (cH && fLoader->TreeH())
  {
     if (fHits == 0x0) fHits  = new TClonesArray("AliFITHits",  405);
     AliDetector::MakeBranch(option);
  } 
  if (cD && fLoader->TreeD())
    {
      MakeBranchInTree(fLoader->TreeD(), GetName(),
      		       &fDigits, 100, 0)->SetAddress(&fDigits);
      //      fLoader->TreeD()->Branch(branchname.Data(),"AliFITDigit",fDigits);   
    } 
   if (cS && fLoader->TreeS())
    {
      MakeBranchInTree(fLoader->TreeD(), branchname,
      		       &fDigits, 405, 0);
      // fLoader->TreeS()->Branch(branchname,"AliFITDigit",&fDigits);
    } 
  
}    

//_____________________________________________________________________________
void AliFIT::ResetHits()
{
  //
  //reset hits
  //
  AliDetector::ResetHits();
  
}
//____________________________________________________________________
void AliFIT::ResetDigits()
{
  //
  // Reset number of digits and the digits array for this detector
  //
  if (fDigits) fDigits->Clear();
  fNdigits = 0;
}

//_____________________________________________________________________________
void AliFIT::SetTreeAddress()
{

  TTree    *treeH = fLoader->TreeH();
  
  if (treeH)
    {
      if (fHits == 0x0) fHits  = new TClonesArray("AliFITHits",  405);
    }
    
  AliDetector::SetTreeAddress();
  TTree *treeD = fLoader->TreeD();
  if (treeD) {
    if (fDigits == 0x0)  fDigits  = new TClonesArray("AliFITDigit",100);
    TBranch* branch = treeD->GetBranch ("FIT");
    if (branch) branch->SetAddress(&fDigits);
  }

  // SDigitizer for Federico
  TTree *treeS = fLoader->TreeS();
  if (treeS) {
    //    if (fDigits == 0x0)  fDigits  = new AliFITDigit();
    TBranch* branch = treeS->GetBranch ("FIT");
    if (branch) branch->SetAddress(&fDigits);
  }
 
}


//_____________________________________________________________________________
AliDigitizer* AliFIT::CreateDigitizer(AliDigitizationInput* digInput) const
{

  return new AliFITDigitizer(digInput);
}

//-------------------------------------------------------------------
void AliFIT::Digits2Raw()
{
//
// Starting from the FIT digits, writes the Raw Data objects
//

  fLoader ->LoadDigits("read");
  TTree* treeD = fLoader->TreeD();
  if (!treeD) {
    AliError("no digits tree");
    return;
  }
  TBranch *branch = treeD->GetBranch("FIT");

  AliFITRawData rawWriter;
  rawWriter.SetVerbose(10);
  
  AliDebug(2,Form(" Formatting raw data for FIT "));
  treeD->GetEntry(0);
   rawWriter.RawDataFIT(branch);
  
  
  fLoader->UnloadDigits();
  
}

//____________________________________________________________________________
void AliFIT::Raw2Digits(AliRawReader *rawReader,TTree* digitsTree)
{

 //FIT raw data-> digits conversion
 // reconstruct time information from raw data

  TClonesArray* digits = new TClonesArray ("AliFITDigit", 100);
  digitsTree->Branch("FIT", &digits);


 AliFITRawReader myrawreader(rawReader);
 if (!myrawreader.Next())
   AliDebug(1,Form(" no raw data found!! %i", myrawreader.Next()));

 Int_t allData[500];
  for (Int_t i=0; i<500; i++)  allData[i]=0;
 for (Int_t i=0; i<500; i++) 
   if(myrawreader.GetData(i)>0)  allData[i]=myrawreader.GetData(i);
 
 Int_t timeCFD, timeLED, timeQT1, timeQT0;
 for (Int_t ipmt=0; ipmt<160; ipmt++) {
   if(allData[ipmt]>0) {
     timeCFD = allData[ipmt];
     timeLED = allData[ipmt];
     timeQT0= allData[ipmt+160];
     timeQT1 = allData[ipmt+320];
     AddDigit(ipmt,   timeCFD, timeLED, timeQT0,  timeQT1, 0);
   }
 }
 
 

 digitsTree->Fill(); 
 GetLoader()->WriteDigits("OVERWRITE");//write out digits
 ResetDigits();
 
 }
