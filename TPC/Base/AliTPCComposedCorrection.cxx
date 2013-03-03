
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

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliTPCComposedCorrection class                                             //
//                                                                            //
// This class is creating a correction that is composed out of smaller        //
// corrections.                                                               //
// There are two ways the sub-corrections can be combined into this one:      //
// 1. kParallel: All corrections are applied at the given position x and      //
//    the dx terms are summed up (this commutes).                             //
// 2. kQueue: The corrections are called in order. The first one at the       //
//    given position x resulting in dx1, the second one is called at          //
//    the corrected position (x+dx1) resulting in dx2, the third one          //
//    is then called at position (x+dx1+dx2) and so forth. dx=dx1+dx2+...     //
//    is returned.                                                            //
// For the inverse of the correction this is taken into account by reversing  //
// the order the corrections are applied in the kQueue case (no issue for     //
// kParallel).                                                                //
//                                                                            //
// date: 27/04/2010                                                           //
// Authors: Magnus Mager, Stefan Rossegger, Jim Thomas                       //
//                                                                            //
// Example usage:                                                             //
//                                                                            //
//  AliMagF mag("mag","mag");                                                 //
//  AliTPCExBBShape exb;             // B field shape distortions             //
//  exb.SetBField(&mag);                                                      //
//                                                                            //
//  AliTPCExBTwist twist;            // ExB Twist distortions                 //
//  twist.SetXTwist(0.001);                                                   //
//                                                                            //
//   TObjArray cs;  cs.Add(&exb); cs.Add(&twist);                             //
//                                                                            //
//  AliTPCComposedCorrection cc;                                              //
//  cc.SetCorrections(&cs);                                                   //
//  cc.SetOmegaTauT1T2(wt,T1,T2);                                             //
//  cc.Print("DA");                                                           //               
//  cc.CreateHistoDRPhiinZR(0,100,100)->Draw("surf2");                        //
////////////////////////////////////////////////////////////////////////////////


#include <TCollection.h>
#include <TTimeStamp.h>
#include <TIterator.h>
#include "AliLog.h"

#include "AliTPCComposedCorrection.h"


AliTPCComposedCorrection::AliTPCComposedCorrection() 
  : AliTPCCorrection("composed_correction",
		     "composition of corrections"),
    fCorrections(0),
    fMode(kParallel),
    fWeights(0)  // weights of corrections
{
  //
  // default constructor
  //
}

AliTPCComposedCorrection::AliTPCComposedCorrection(TCollection *corrections,
						   AliTPCComposedCorrection::CompositionType mode)
  : AliTPCCorrection("composed_correction",
		     "composition of corrections"),
    fCorrections(corrections),
    fMode(mode),
    fWeights(0) //weights of correction
{
  //
  // Constructor that defines the set of corrections, this one is composed of.
  //
}

AliTPCComposedCorrection::~AliTPCComposedCorrection() {
  // 
  // destructor
  //
  if (!fCorrections) {
    AliInfo("No Correction-models were set: can not delete them");
  } else {
    TIterator *i=fCorrections->MakeIterator();
    AliTPCCorrection *c;
    while (0!=(c=dynamic_cast<AliTPCCorrection*>(i->Next()))) {
      delete c;
    }
    delete i;
  }
  if (fWeights) delete fWeights;
}

AliTPCCorrection * AliTPCComposedCorrection::GetSubCorrection(Int_t ipos){
  //
  //
  //
  TObjArray *arr = (TObjArray*)fCorrections;
  return (AliTPCCorrection *)arr->At(ipos);
}

AliTPCCorrection * AliTPCComposedCorrection::GetSubCorrection(const char *cname){
  //
  //
  //
  TCollection *arr = fCorrections;
  return (AliTPCCorrection *)arr->FindObject(cname);
}



void AliTPCComposedCorrection::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // This applies all correction and the specified manner (see general
  // class description for details).
  //

  if (!fCorrections) {
    AliInfo("No Corrections-models were set: can not calculate distortions");
    return;
  }
  TIterator *i=fCorrections->MakeIterator();
  AliTPCCorrection *c;
  Int_t weightIndex=0;
  switch (fMode) {
  case kParallel:
    Float_t dxi[3];
    for (int j=0;j<3;++j) dx[j]=0.;
    while (0!=(c=dynamic_cast<AliTPCCorrection*>(i->Next()))) {
      c->GetCorrection(x,roc,dxi);
      Double_t w=1;
      if (fWeights) w=(*fWeights)[weightIndex++];
      for (Int_t j=0;j<3;++j) dx[j]+=w*dxi[j];
    }
    break;
  case kQueue:
    Float_t xi[3];
    for (Int_t j=0;j<3;++j) xi[j]=x[j];
    while (0!=(c=dynamic_cast<AliTPCCorrection*>(i->Next()))) {
      c->GetCorrection(xi,roc,dx);
      Double_t w=1;
      if (fWeights) w=(*fWeights)[weightIndex++];
      for (Int_t j=0;j<3;++j) xi[j]+=w*dx[j];
    }
    for (Int_t j=0;j<3;++j) dx[j]=xi[j]-x[j];
    break;
  }
  delete i;
}

void AliTPCComposedCorrection::GetDistortion(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // This applies all distortions and the specified manner (see general
  // class descxiption for details).
  //

  if (!fCorrections) {
    AliInfo("No Corrections-models were set: can not calculate distortions");
    return;
  }
  TIterator *i=fCorrections->MakeReverseIterator();
  AliTPCCorrection *c;
  Int_t weightIndex=0;
  switch (fMode) {
  case kParallel:
    Float_t dxi[3];
    for (int j=0;j<3;++j) dx[j]=0.;
    while (0!=(c=dynamic_cast<AliTPCCorrection*>(i->Next()))) {
      c->GetDistortion(x,roc,dxi);
      Double_t w=1;
      if (fWeights) w=(*fWeights)[weightIndex++];
      for (Int_t j=0;j<3;++j) dx[j]+=w*dxi[j];
    }
    break;
  case kQueue:
    Float_t xi[3];
    for (Int_t j=0;j<3;++j) xi[j]=x[j];
    while (0!=(c=dynamic_cast<AliTPCCorrection*>(i->Next()))) {
      c->GetDistortion(xi,roc,dx);
      Double_t w=1;
      if (fWeights) w=(*fWeights)[weightIndex++];
      for (Int_t j=0;j<3;++j) xi[j]+=w*dx[j];
    }
    for (Int_t j=0;j<3;++j) dx[j]=xi[j]-x[j];
    break;
  }
  delete i;
}


void AliTPCComposedCorrection::Print(Option_t* option) const {
  //
  // Print function to check which correction classes are used 
  // option=="d" prints details regarding the setted magnitude 
  // option=="a" prints the C0 and C1 coefficents for calibration purposes
  //

  printf("Composed TPC spacepoint correction \"%s\" -- composed of:\n",GetTitle());
  TString opt = option; opt.ToLower();
  Int_t in=1;
  if (!fCorrections) {
    printf("   - composed correction is empty!\n");
    return;
  }
  TIterator *i=fCorrections->MakeIterator();
  AliTPCCorrection *c;
  while (0!=(c=dynamic_cast<AliTPCCorrection*>(i->Next()))) {
    if (opt.Contains("d")) {
      printf("\n");
      printf("%d. %s\t%s\n",in,c->GetTitle(), c->GetName());
      c->Print(option);
    } else {
      printf("%d. %s\t%s\n",in,c->GetTitle(), c->GetName());
    }
    ++in;
  }
  if (in==1) printf("  Info: The correction compound is empty: No corrections set\n");
  delete i;
}


void AliTPCComposedCorrection::Init() {
  //
  // Initialization funtion (not used at the moment)
  //
  if (!fCorrections) {
    AliInfo("No Correction-models were set");
    return;
  }
  TIterator *i=fCorrections->MakeIterator();
  AliTPCCorrection *c;
  while (0!=(c=dynamic_cast<AliTPCCorrection*>(i->Next()))) 
    c->Init();
  delete i;
  
}

void AliTPCComposedCorrection::Update(const TTimeStamp &timeStamp) {
  //
  // Update function 
  //
  if (!fCorrections) {
    AliInfo("No Correction-models were set");
    return;
  }

  TIterator *i=fCorrections->MakeIterator();
  AliTPCCorrection *c;
  while (0!=(c=dynamic_cast<AliTPCCorrection*>(i->Next()))) 
    c->Update(timeStamp);
  delete i;
 
}



void AliTPCComposedCorrection::SetOmegaTauT1T2(Float_t omegaTau,Float_t t1,Float_t t2) {
  //
  // Gives the possibility to set the OmegaTau plus Tensor corrections T1 and T2 (effective omega Tau)
  // to each subcorrection (since they might become event specific due to changing drift velocity)
  //
  // The omegaTau comes idealy from the Database, since it is a function of drift velocity, B and E field 
  // e.g. omegaTau  = -10.0 * Bz * vdrift / Ez ; // with Bz in kG and Ez in V/cm
  //      omegaTau  = -0.325 for Bz=5kG, Ez=400V/cm and vdrift = 2.6cm/muSec
  // The T1 and T2 tensors were measured in a dedicated calibration run
  //
  // Note: overwrites previously set values!
  // 

  if (!fCorrections) {
    AliInfo("No Correction-models were set");
    return;
  }

  TIterator *i=fCorrections->MakeIterator();
  AliTPCCorrection *c;
  while (0!=(c=dynamic_cast<AliTPCCorrection*>(i->Next()))) {
    c->SetOmegaTauT1T2(omegaTau,t1,t2);
  }
  delete i;
}

ClassImp(AliTPCComposedCorrection)
