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

/*
$Log$
*/

///////////////////////////////////////////////////////////////////
//                                                               //
//    Generate the final state of the interaction as the input   //
//    to the MonteCarlo                                          //
//
//Begin_Html
/*
<img src="picts/AliGeneratorClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:andreas.morsch@cern.ch">Andreas Morsch</a>.
</font>
<pre>
*/
//End_Html
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliGenFixed.h"
#include "AliRun.h"
#include "AliConst.h"
#include "AliPDG.h"
#include "TF1.h"
  
ClassImp(AliGenFixed)

//_____________________________________________________________________________
AliGenFixed::AliGenFixed()
  :AliGenerator()
{
  //
  // Default constructor
  //
  fIpart = 0;
}

//_____________________________________________________________________________
AliGenFixed::AliGenFixed(Int_t npart)
  :AliGenerator(npart)
{
  //
  // Standard constructor
  //
  fName="Fixed";
  fTitle="Fixed Particle Generator";
  // Generate Proton by default
  fIpart=kProton;
}

//_____________________________________________________________________________
void AliGenFixed::Generate()
{
  //
  // Generate one trigger
  //
  Float_t polar[3]= {0,0,0};
  Float_t p[3] = {fPMin*TMath::Cos(fPhiMin)*TMath::Sin(fThetaMin),
		  fPMin*TMath::Sin(fPhiMin)*TMath::Sin(fThetaMin),
		  fPMin*TMath::Cos(fThetaMin)};
  Int_t i, nt;
  //
  for(i=0;i<fNpart;i++) {
    gAlice->SetTrack(fTrackIt,-1,fIpart,p,fOrigin.GetArray(),polar,0,"Primary",nt);
  }
}
  
//_____________________________________________________________________________
void AliGenFixed::SetSigma(Float_t sx, Float_t sy, Float_t sz)
{
  //
  // Set the interaction point sigma
  //
  printf("Vertex smearing not implemented for fixed generator\n");
}
