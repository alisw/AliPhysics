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
Revision 1.6  2000/07/11 18:24:59  fca
Coding convention corrections + few minor bug fixes

Revision 1.5  2000/06/08 13:34:50  fca
Better control of momentum range in GenBox

Revision 1.4  1999/09/29 09:24:29  fca
Introduction of the Copyright and cvs Log

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

#include "AliGenerator.h"
#include "AliRun.h"

ClassImp(AliGenerator)

TGenerator* AliGenerator::fgMCEvGen=0;

//____________________________________________________________
AliGenerator::AliGenerator()
{
  //
  // Default constructor
  //
    printf("\n AliGenerator Default Constructor\n\n");
    
    gAlice->SetGenerator(this);
    SetThetaRange(); ResetBit(kThetaRange);
    SetPhiRange(); ResetBit(kPhiRange);
    SetMomentumRange(); ResetBit(kMomentumRange);
    SetPtRange(); ResetBit(kPtRange);
    SetYRange(); ResetBit(kYRange);
    SetNumberParticles();
    SetTrackingFlag();

    fOrigin.Set(3);
    fOsigma.Set(3);
    fOrigin[0]=fOrigin[1]=fOrigin[2]=0;
    fOsigma[0]=fOsigma[1]=fOsigma[2]=0;
    fVMin.Set(3);
    fVMin[0]=fVMin[1]=fVMin[2]=0;
    fVMax.Set(3);
    fVMax[0]=fVMax[1]=fVMax[2]=10000;
}

//____________________________________________________________
AliGenerator::AliGenerator(Int_t npart)
    : TNamed(" "," ")
{
  //
  // Standard constructor
  //
    printf("\n AliGenerator Constructor initializing number of particles \n\n");
    gAlice->SetGenerator(this);
    SetThetaRange(); ResetBit(kThetaRange);
    SetPhiRange(); ResetBit(kPhiRange);
    SetMomentumRange(); ResetBit(kMomentumRange);
    SetPtRange(); ResetBit(kPtRange);
    SetYRange(); ResetBit(kYRange);
    SetTrackingFlag();

    fOrigin.Set(3);
    fOsigma.Set(3);
    fOrigin[0]=fOrigin[1]=fOrigin[2]=0;
    fOsigma[0]=fOsigma[1]=fOsigma[2]=0;
    fVMin.Set(3);
    fVMin[0]=fVMin[1]=fVMin[2]=0;
    fVMax.Set(3);
    fVMax[0]=fVMax[1]=fVMax[2]=10000;

    SetNumberParticles(npart);
}

//____________________________________________________________
AliGenerator::AliGenerator(const AliGenerator &gen) : TNamed(" "," ")
{
  //
  // Copy constructor
  //
  gen.Copy(*this);
}

//____________________________________________________________
AliGenerator & AliGenerator::operator=(const AliGenerator &gen)
{
  //
  // Assignment operator
  //
  gen.Copy(*this);
  return (*this);
}

//____________________________________________________________
void AliGenerator::Copy(AliGenerator &/* gen */) const
{
  //
  // Copy *this onto gen
  //
  Fatal("Copy","Not implemented!\n");
}

//____________________________________________________________
AliGenerator::~AliGenerator()
{
  //
  // Destructor
  //
  fOrigin.Set(0);
  fOsigma.Set(0);
  delete fgMCEvGen;
}

void AliGenerator::Init()
{   
  //
  // Dummy initialisation
  //
}


