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
Revision 1.11  2001/01/26 19:58:48  hristov
Major upgrade of AliRoot code

Revision 1.10  2000/12/21 15:30:18  fca
Correcting coding convention violations

Revision 1.9  2000/10/04 10:08:01  fca
Correction of minor typing mistakes

Revision 1.8  2000/10/02 21:28:14  fca
Removal of useless dependecies via forward declarations

Revision 1.7  2000/07/12 08:56:25  fca
Coding convention correction and warning removal

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
#include "TGenerator.h"
#include "AliRun.h"

ClassImp(AliGenerator)

TGenerator* AliGenerator::fgMCEvGen=0;

//____________________________________________________________
AliGenerator::AliGenerator()
{
  //
  // Default constructor
  //
    if (gAlice->GetDebug()>0)
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
    fVertex.Set(3);

    fOrigin[0]=fOrigin[1]=fOrigin[2]=0;
    fOsigma[0]=fOsigma[1]=fOsigma[2]=0;
    fVertex[0]=fVertex[1]=fVertex[2]=0;

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
    if (gAlice->GetDebug()>0)
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
    fVertex.Set(3);

    fOrigin[0]=fOrigin[1]=fOrigin[2]=0;
    fOsigma[0]=fOsigma[1]=fOsigma[2]=0;
    fVertex[0]=fVertex[1]=fVertex[2]=0;

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

//_______________________________________________________________________
void AliGenerator::SetOrigin(Float_t ox, Float_t oy, Float_t oz)
{
  //
  // Set the vertex for the generated tracks
  //
  fOrigin[0]=ox;
  fOrigin[1]=oy;
  fOrigin[2]=oz;
}

//_______________________________________________________________________
void AliGenerator::SetOrigin(const TLorentzVector &o)
{
  //
  // Set the vertex for the generated tracks
  //
  fOrigin[0]=o[0];
  fOrigin[1]=o[1];
  fOrigin[2]=o[2];
}

//_______________________________________________________________________
void AliGenerator::SetSigma(Float_t sx, Float_t sy, Float_t sz)
{
  //
  // Set the spread of the vertex
  //
  fOsigma[0]=sx;
  fOsigma[1]=sy;
  fOsigma[2]=sz;
}

//_______________________________________________________________________
void AliGenerator::SetMomentumRange(Float_t pmin, Float_t pmax)
{
  //
  // Set the momentum range for the generated particles
  //
  fPMin = pmin;
  fPMax = pmax;
  SetBit(kMomentumRange);
}

//_______________________________________________________________________
void AliGenerator::SetPtRange(Float_t ptmin, Float_t ptmax)
{
  //
  // Set the Pt range for the generated particles
  //
  fPtMin = ptmin;
  fPtMax = ptmax;
  SetBit(kPtRange);
}

//_______________________________________________________________________
void AliGenerator::SetPhiRange(Float_t phimin, Float_t phimax)
{
  //
  // Set the Phi range for the generated particles
  //
  fPhiMin = TMath::Pi()*phimin/180;
  fPhiMax = TMath::Pi()*phimax/180; SetBit(kPhiRange);
}

//_______________________________________________________________________
void AliGenerator::SetYRange(Float_t ymin, Float_t ymax)
{
  //
  // Set the Rapidity range for the generated particles
  //
  fYMin=ymin;
  fYMax=ymax;
  SetBit(kYRange);
}

//_______________________________________________________________________
void AliGenerator::SetVRange(Float_t vxmin, Float_t vxmax,
			     Float_t vymin, Float_t vymax,
			     Float_t vzmin, Float_t vzmax)
{
  //
  // Set the vertex range for the generated particles
  //
  fVMin[0]=vxmin; fVMin[1]=vymin; fVMin[2]=vzmin;
  fVMax[0]=vxmax; fVMax[1]=vymax; fVMax[2]=vzmax;
  SetBit(kVertexRange);
}

//_______________________________________________________________________
void AliGenerator::SetThetaRange(Float_t thetamin, Float_t thetamax)
{
  //
  // Set the theta range for the generated particles
  //
  fThetaMin = TMath::Pi()*thetamin/180;
  fThetaMax = TMath::Pi()*thetamax/180; SetBit(kThetaRange);
}

void AliGenerator::Vertex()
{
  //
  // Obtain vertex for current event from external source or calculated (internal)
  //
    if (fVertexSource == kInternal) {
	VertexInternal();
    } else {
	VertexExternal();
    }
}


void AliGenerator::VertexExternal()
{
    // Dummy !!!!!!
    // Obtain vertex from external source 
    //
    // Should be something like fVertex = gAlice->GetVertex()
    
    fVertex[0]=fVertex[1]=fVertex[2]=0;  
}

void AliGenerator::VertexInternal()
{
    // 
    // Obtain calculated vertex 
    // Default is gaussian smearing
    Float_t random[6];
    Rndm(random,6);
    for (Int_t j = 0; j<3 ; j++) {
	fVertex[j]=
	    fOrigin[j]+fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
	    TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
    }
}


