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
Revision 1.8  2003/01/14 10:50:18  alibrary
Cleanup of STEER coding conventions

Revision 1.7  2002/02/08 16:50:50  morsch
Add name and title in constructor.

Revision 1.6  2001/07/27 17:09:35  morsch
Use local SetTrack, KeepTrack and SetHighWaterMark methods
to delegate either to local stack or to stack owned by AliRun.
(Piotr Skowronski, A.M.)

Revision 1.5  2000/12/21 16:24:06  morsch
Coding convention clean-up

Revision 1.4  2000/11/30 07:12:49  alibrary
Introducing new Rndm and QA classes

Revision 1.3  2000/10/02 21:28:06  fca
Removal of useless dependecies via forward declarations

Revision 1.2  2000/07/11 18:24:55  fca
Coding convention corrections + few minor bug fixes

Revision 1.1  2000/06/09 20:22:58  morsch
Same class as previously in AliSimpleGen.cxx
All coding rule violations except RS3 corrected (AM)

*/



// Generator for particles in a preset
// kinematic range (flat distribution)
// Note that for a given theta pt and p are not independent 
// Range for only one variable (pt or p) should be given.
//
// Comments and suggestions: andreas.morsch@cern.ch
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

#include "TPDGCode.h"

#include "AliConst.h"
#include "AliGenBox.h"
#include "AliRun.h"

ClassImp(AliGenBox)

//_____________________________________________________________________________
AliGenBox::AliGenBox()
    :AliGenerator()
{
  //
  // Default constructor
  //
  fIpart=0;
}

//_____________________________________________________________________________
AliGenBox::AliGenBox(Int_t npart)
  :AliGenerator(npart)
{
  //
  // Standard constructor
  //
  fName  = "Box";
  fTitle = "Box particle generator";
  // Generate Proton by default
  fIpart=kProton;
}

//_____________________________________________________________________________

void AliGenBox::Generate()
{
  //
  // Generate one trigger
  //
  
    Float_t polar[3]= {0,0,0};
  //
    Float_t origin[3];
    Float_t p[3];
    Int_t i, j, nt;
    Double_t pmom, theta, phi, pt;
    //
    Float_t random[6];
  //
    for (j=0;j<3;j++) origin[j]=fOrigin[j];
    if(fVertexSmear==kPerEvent) {
	Vertex();
	for (j=0;j<3;j++) origin[j]=fVertex[j];
    }

    for(i=0;i<fNpart;i++) {
	Rndm(random,3);
	theta=fThetaMin+random[0]*(fThetaMax-fThetaMin);
	if(TestBit(kMomentumRange)) {
	    pmom=fPMin+random[1]*(fPMax-fPMin);
	    pt=pmom*TMath::Sin(theta);
	} else {

	    pt=fPtMin+random[1]*(fPtMax-fPtMin);
	    pmom=pt/TMath::Sin(theta);
	}
	phi=fPhiMin+random[2]*(fPhiMax-fPhiMin);
	p[0] = pt*TMath::Cos(phi);
	p[1] = pt*TMath::Sin(phi);
	p[2] = pmom*TMath::Cos(theta);

	if(fVertexSmear==kPerTrack) {
	    Rndm(random,6);
	    for (j=0;j<3;j++) {
		origin[j]=fOrigin[j]+fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
		    TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	    }
	}
	SetTrack(fTrackIt,-1,fIpart,p,origin,polar,0,kPPrimary,nt);
    }
}

//_____________________________________________________________________________

void AliGenBox::Init()
{
// Initialisation, check consistency of selected ranges
  if(TestBit(kPtRange)&&TestBit(kMomentumRange)) 
    Fatal("Init","You should not set the momentum range and the pt range!\n");
  if((!TestBit(kPtRange))&&(!TestBit(kMomentumRange))) 
    Fatal("Init","You should set either the momentum or the pt range!\n");
}

