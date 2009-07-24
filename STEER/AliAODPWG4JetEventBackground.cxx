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
//     AOD class for PWG4 jet backgrounds
//     Author: Christian Klein-Boesing IKP Muenster
//-------------------------------------------------------------------------


#include "AliAODPWG4JetEventBackground.h"

using namespace std;

ClassImp(AliAODPWG4JetEventBackground)

TString AliAODPWG4JetEventBackground::fgkStdBranchName("jeteventbackground");



//______________________________________________________________________________
AliAODPWG4JetEventBackground::AliAODPWG4JetEventBackground() :
    TObject()
{
  for(int i = 0;i < kMaxBackground;++i){
    fBackground[i] = 0;
  } 
}

//______________________________________________________________________________
AliAODPWG4JetEventBackground::~AliAODPWG4JetEventBackground() 
{
  //
  // destructor
  //
}

//______________________________________________________________________________
AliAODPWG4JetEventBackground::AliAODPWG4JetEventBackground(const AliAODPWG4JetEventBackground& back) :
    TObject(back)
{
  //
  // Copy constructor
  //
  for(int i = 0;i < kMaxBackground;++i){
    fBackground[i] = back.fBackground[i];
  } 
  
}

//______________________________________________________________________________
AliAODPWG4JetEventBackground& AliAODPWG4JetEventBackground::operator=(const AliAODPWG4JetEventBackground& back)
{
  //
  // Assignment operator
  //

  if(this!=&back) {
    TObject::operator=(back);
    for(int i = 0;i < kMaxBackground;++i){
      fBackground[i] = back.fBackground[i];
    } 
  }

  return *this;
}

void AliAODPWG4JetEventBackground::Print(Option_t* /*option*/) const 
{
  //
  // Print information of all data members
  //

  printf("Jet EventBackground :\n");
  for(int i = 0;i < kMaxBackground;++i){
    printf("%d: %3.E GeV \n",i,fBackground[i]);
  } 
}
