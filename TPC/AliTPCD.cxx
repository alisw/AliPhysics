//-----------------------------------------------------------------------------
// 
//  Origin: Marian Ivanov, Uni. of Bratislava, ivanov@fmph.uniba.sk
//
//  Implementation of class AliTPCD
//
//-----------------------------------------------------------------------------

#include <TMath.h>
#include "AliTPCParam.h"
#include "AliTPCRF1D.h"
#include "AliTPCPRF2D.h"
#include "AliTPCD.h"
#include "TClonesArray.h"
#include "TClass.h"
#include "TBranchClones.h"
#include "TTree.h"   
#include "TDirectory.h"



// other include files follow here


ClassImp(AliTPCD)
  //_____________________________________________________________ 
AliTPCD::AliTPCD(Text_t *  name,
		 AliTPCParam * param, AliTPCPRF2D* prf, AliTPCRF1D* prfz) 
{
  //construct new object's or accept objects sent to constructor
  //AliTPCD handle sent object and is repsonsible for
  //deleting it

//Begin_Html
/*
<img src="gif/alitpcd.gif">
*/
//End_Html
  SetName(name);
  if ((param!=0) && 
      ( (param->IsA()->InheritsFrom("AliTPCParam")==kTRUE ) ))
    fParam=param; 
  else
    fParam= new AliTPCParam;
  if ( (prf!=0) && (prf->IsA()->InheritsFrom("AliTPCPRF2D")==kTRUE) )
    fPRF=prf;
  else
    fPRF =  new AliTPCPRF2D;
  if ( (prfz!=0) && (prfz->IsA()->InheritsFrom("AliTPCRF1D")==kTRUE) )
    fRF=prfz;
  else
    fRF =  new AliTPCRF1D(kTRUE);  
  fDigits = new TClonesArray("AliTPCdigit",5000);
  fpthis=this;
}

AliTPCD::~AliTPCD() 
{
  if (fParam!=0) fParam->Delete();
  if (fPRF!=0) fPRF->Delete();
  if (fRF!=0) fRF->Delete();
  if (fDigits!=0) fDigits->Delete();
}


Bool_t  AliTPCD::SetTree(Int_t nevent, TDirectory *dir )
{
  char treeName[100];
  // Get Hits Tree header from file
  sprintf(treeName,"TreeD%d_%s",nevent,GetName());
  fTreeD = (TTree*)dir->Get(treeName);
  if (fTreeD == 0) return kFALSE;
  //set Digit branch 
  TBranch *b = fTreeD->GetBranch("Digits");
  if (b==0) return kFALSE;
  b->SetAddress(&fDigits);
  return kTRUE;
}


Bool_t  AliTPCD::MakeTree(Int_t nevent)
{
  char treeName[100];
  // Get Hits Tree header from file
  sprintf(treeName,"TreeD%d_%s",nevent,GetName());
  fTreeD =  new TTree(treeName,treeName);
  if (fTreeD == 0) return kFALSE;
  //set Digit branch 
  TBranch *b = fTreeD->Branch("Digits",&fDigits,40000);
  if (b==0) return kFALSE;
  b->SetAddress(&fDigits);
 
  return kTRUE;
}
  
void AliTPCD::Fill()
{ 
  if (fTreeD!=0) fTreeD->Fill();    
}





void AliTPCD::Streamer(TBuffer &R__b)
{
 // Stream an object of class AliTPCD. 
  if (R__b.IsReading()) {    
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    TNamed::Streamer(R__b);  
    if (fParam!=0) {
      fParam->Delete();
      fParam = new AliTPCParam;
    }
    if (fPRF!=0) {
      fPRF->Delete();
      fPRF = new AliTPCPRF2D;
    }
    if (fRF!=0) {
      fRF->Delete();
      fRF = new AliTPCRF1D;
    }
    if (fTreeD!=0) {
      fRF->Delete();
      fRF = new AliTPCRF1D;
    }
    R__b >>fParam;
    R__b >>fPRF;
    R__b >>fRF;    
    SetTree();        
  } else {
    R__b.WriteVersion(AliTPCD::IsA());
    TNamed::Streamer(R__b);     
    R__b <<fParam;
    R__b <<fPRF;
    R__b <<fRF;     
    if (fTreeD!=0) fTreeD->Write();         
  } 
}
