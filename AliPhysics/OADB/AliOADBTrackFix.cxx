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

/* $Id$ */

//-------------------------------------------------------------------------
//     OADB class for run dependent track fixing parameters
//     Convention for phi-dependent data: 0 : 2pi
//     Author: ruben.shahoyan@cern.ch
//-------------------------------------------------------------------------

#include <TGraph.h>
#include <TMath.h>
#include "AliOADBTrackFix.h"
#include "AliLog.h"

ClassImp(AliOADBTrackFix);

//______________________________________________________________________________
AliOADBTrackFix::AliOADBTrackFix()
{
  // Default constructor
  for (int imd=0;imd<kNCorModes;imd++) {
    for (int iside=0;iside<2;iside++) fPtInvCor[imd][iside] = 0;  
    fXIniPtInvCorr[imd] = 0;
  }
  //
}
//______________________________________________________________________________
AliOADBTrackFix::AliOADBTrackFix(const char* name) : 
  TNamed(name, "TrackFix")
{
  // Constructor
  for (int imd=0;imd<kNCorModes;imd++) {
    for (int iside=0;iside<2;iside++) fPtInvCor[imd][iside] = 0;  
    fXIniPtInvCorr[imd] = 0;
  }  
}

//______________________________________________________________________________
AliOADBTrackFix::~AliOADBTrackFix() 
{
  // destructor
  for (int imd=0;imd<kNCorModes;imd++) for (int iside=0;iside<2;iside++) delete fPtInvCor[imd][iside];
  //
}

//______________________________________________________________________________
void AliOADBTrackFix::SetPtInvCorr(int mode,int side, const TGraph* gr)
{
  if (!gr || gr->GetN()<1) {
    AliInfo(Form("Correction for side %d in mode %d is empty",side,mode));
    fPtInvCor[mode][side] = 0;
  }
  fPtInvCor[mode][side] = gr;
}

//______________________________________________________________________________
Double_t AliOADBTrackFix::GetPtInvCorr(int mode, double sideAfrac, double phi) const
{
  // calculate the correction for the track of given model at given phi, provided the sideAfrac fraction of its
  // total lenght in TPC is in side A.
  if (!fPtInvCor[mode][0] || !fPtInvCor[mode][1]) return 0; // no graph 
  while (phi>2*TMath::Pi()) phi -= 2*TMath::Pi();
  while (phi<0) phi += 2*TMath::Pi();  
  int nb = fPtInvCor[mode][0]->GetN();
  int bin = int( phi/(2*TMath::Pi())*nb );
  if (bin==nb) bin = nb-1;
  return sideAfrac*fPtInvCor[mode][0]->GetY()[bin] + (1.-sideAfrac)*fPtInvCor[mode][1]->GetY()[bin];
}


