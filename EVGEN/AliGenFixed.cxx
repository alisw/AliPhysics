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
Revision 1.4  2000/12/21 16:24:06  morsch
Coding convention clean-up

Revision 1.3  2000/11/30 07:12:50  alibrary
Introducing new Rndm and QA classes

Revision 1.2  2000/10/02 15:17:54  morsch
Unused includes removed.

Revision 1.1  2000/06/09 20:24:00  morsch
Same class as previously in AliSimpleGen.cxx
All coding rule violations except RS3 corrected (AM)

*/



// Simple particle gun. 
// Momentum, phi and theta of the partice as well as the particle type can be set.
// andreas.morsch@cern.ch
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
#include "AliPDG.h"
  
ClassImp(AliGenFixed)

//_____________________________________________________________________________
AliGenFixed::AliGenFixed()
  :AliGenerator()
{
  //
  // Default constructor
  //
  fIpart = 0;
  fExplicit = kFALSE;
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
  if(!fExplicit) {
    fP[0] = fPMin*TMath::Cos(fPhiMin)*TMath::Sin(fThetaMin);
    fP[1] = fPMin*TMath::Sin(fPhiMin)*TMath::Sin(fThetaMin);
    fP[2] = fPMin*TMath::Cos(fThetaMin);
  }
  Int_t i, nt;
  //
  for(i=0;i<fNpart;i++) 
    gAlice->SetTrack(fTrackIt,-1,fIpart,fP,fOrigin.GetArray(),polar,0,kPPrimary,nt);
  
}
  
//_____________________________________________________________________________
void AliGenFixed::SetSigma(Float_t sx, Float_t sy, Float_t sz)
{
  //
  // Set the interaction point sigma
  //
  printf("Vertex smearing not implemented for fixed generator\n");
}
