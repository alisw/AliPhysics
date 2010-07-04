/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-------------------------------------------------------------------------
//     AOD class for  jet backgrounds
//     Author: Christian Klein-Boesing IKP Muenster
//-------------------------------------------------------------------------


#include "AliAODJetEventBackground.h"

using namespace std;

ClassImp(AliAODJetEventBackground)

TString AliAODJetEventBackground::fgkStdBranchName("jeteventbackground");



//______________________________________________________________________________
AliAODJetEventBackground::AliAODJetEventBackground() :
    TNamed()
{
  for(int i = 0;i < kMaxBackground;++i){
    fBackground[i] = 0;
    fSigma[i] = 0;
    fMeanArea[i] = 0; 

  } 
}

//______________________________________________________________________________
AliAODJetEventBackground::~AliAODJetEventBackground() 
{
  //
  // destructor
  //
}

//______________________________________________________________________________
AliAODJetEventBackground::AliAODJetEventBackground(const AliAODJetEventBackground& back) :
    TNamed(back)
{
  //
  // Copy constructor
  //
  for(int i = 0;i < kMaxBackground;++i){
    fBackground[i] = back.fBackground[i];
    fSigma[i] = back.fSigma[i];
    fMeanArea[i] = back.fMeanArea[i]; 
  } 
  
}

//______________________________________________________________________________
AliAODJetEventBackground& AliAODJetEventBackground::operator=(const AliAODJetEventBackground& back)
{
  //
   // Assignment operator
  //

  if(this!=&back) {
    TNamed::operator=(back);
    for(int i = 0;i < kMaxBackground;++i){
      fBackground[i] = back.fBackground[i];
      fSigma[i] = back.fSigma[i];
      fMeanArea[i] = back.fMeanArea[i]; 
    } 
  }

  return *this;
}

void AliAODJetEventBackground::Print(Option_t* /*option*/) const 
{
  //
  // Print information of all data members
  //

  printf("Jet EventBackground :\n");
  for(int i = 0;i < kMaxBackground;++i){
    printf("%d: %3.E GeV Sigma %3.E Mean Area %3.E  \n",i,fBackground[i],fSigma[i],fMeanArea[i]);
  } 
}

void AliAODJetEventBackground::Reset()  
{
  //
  // reset information of all data members
  //
  for(int i = 0;i < kMaxBackground;++i){
    fBackground[i] = 0;
    fSigma[i] = 0;
    fMeanArea[i] = 0; 
  } 
}
