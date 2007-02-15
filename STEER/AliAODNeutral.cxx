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
//     AOD track base class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include "AliAODNeutral.h"

ClassImp(AliAODNeutral)

//______________________________________________________________________________
AliAODNeutral::AliAODNeutral() : 
  AliVirtualParticle(),
  fEnergy(0),
  fChi2(-999.),
  fID(-999),
  fLabel(-999),
  fCovMatrix(NULL),
  fProdVertex(0x0),
  fPrimTrack(NULL),
  fType(kUndef)
{
  // default constructor

  SetPosition((Float_t*)NULL);
  SetPID((Float_t*)NULL);
}

//______________________________________________________________________________
AliAODNeutral::AliAODNeutral(Int_t id,
			     Int_t label, 
			     Double_t energy,
			     Double_t x[3],
			     Bool_t isDCA,
			     Double_t covMatrix[21],
			     Double_t pid[10],
			     AliAODVertex *prodVertex,
			     AliAODTrack *primTrack,
			     Char_t ttype) :
  AliVirtualParticle(),
  fEnergy(energy),
  fChi2(-999.),
  fID(id),
  fLabel(label),
  fCovMatrix(NULL),
  fProdVertex(prodVertex),
  fPrimTrack(primTrack),
  fType(ttype)
{
  // constructor
 
  SetPosition(x, isDCA);
  if(covMatrix) SetCovMatrix(covMatrix);
  SetPID(pid);

}

//______________________________________________________________________________
AliAODNeutral::AliAODNeutral(Int_t id,
			     Int_t label, 
			     Float_t energy,
			     Float_t x[3],
			     Bool_t isDCA,
			     Float_t covMatrix[21],
			     Float_t pid[10],
			     AliAODVertex *prodVertex,
			     AliAODTrack *primTrack,
			     Char_t ttype) :
  AliVirtualParticle(),
  fEnergy(energy),
  fChi2(-999.),
  fID(id),
  fLabel(label),
  fCovMatrix(NULL),
  fProdVertex(prodVertex),
  fPrimTrack(primTrack),
  fType(ttype)
{
  // constructor
 
  SetPosition(x, isDCA);
  if(covMatrix) SetCovMatrix(covMatrix);
  SetPID(pid);

}


//______________________________________________________________________________
AliAODNeutral::~AliAODNeutral() 
{
  // destructor
  delete fCovMatrix;
}


//______________________________________________________________________________
AliAODNeutral::AliAODNeutral(const AliAODNeutral& trk) :
  AliVirtualParticle(trk),
  fEnergy(trk.fEnergy),
  fChi2(trk.fChi2),
  fID(trk.fID),
  fLabel(trk.fLabel),
  fCovMatrix(NULL),
  fProdVertex(trk.fProdVertex),
  fPrimTrack(trk.fPrimTrack),
  fType(trk.fType)
{
  // Copy constructor

  trk.GetPosition(fPosition);
  if(trk.fCovMatrix) fCovMatrix=new AliAODNeuCov(*trk.fCovMatrix);
  SetPID(trk.fPID);

}

//______________________________________________________________________________
AliAODNeutral& AliAODNeutral::operator=(const AliAODNeutral& trk)
{
  // Assignment operator
  if(this!=&trk) {

    AliVirtualParticle::operator=(trk);

    trk.GetPosition(fPosition);
    trk.GetPID(fPID);

    fChi2 = trk.fEnergy;
    fChi2 = trk.fChi2;

    fID = trk.fID;
    fLabel = trk.fLabel;    
    
    delete fCovMatrix;
    if(trk.fCovMatrix) fCovMatrix=new AliAODNeuCov(*trk.fCovMatrix);
    else fCovMatrix=NULL;
    fProdVertex = trk.fProdVertex;
    fPrimTrack = trk.fPrimTrack;

    fType = trk.fType;
  }

  return *this;
}

//______________________________________________________________________________
template <class T> void AliAODNeutral::SetPosition(const T *x, const Bool_t dca) 
{
  // set the position

  if (x) {
    if (!dca) {
      ResetBit(kIsDCA);

      fPosition[0] = x[0];
      fPosition[1] = x[1];
      fPosition[2] = x[2];
    } else {
      SetBit(kIsDCA);
      // don't know any better yet
      fPosition[0] = -999.;
      fPosition[1] = -999.;
      fPosition[2] = -999.;
    }
  } else {
    ResetBit(kIsDCA);

    fPosition[0] = -999.;
    fPosition[1] = -999.;
    fPosition[2] = -999.;
  }
}

//______________________________________________________________________________
void AliAODNeutral::SetDCA(Double_t d, Double_t z) 
{
  // set the dca
  fPosition[0] = d;
  fPosition[1] = z;
  fPosition[2] = 0.;
  SetBit(kIsDCA);
}

//______________________________________________________________________________
void AliAODNeutral::Print(Option_t* /* option */) const
{
  // prints information about AliAODNeutral

  printf("Object name: %s   Neutral type: %s\n", GetName(), GetTitle()); 
  printf("    energy = %f\n", E());
  printf("      chi2 = %f\n", Chi2());
  printf(" PID object: %p\n", PID());
}

//-------------------------------------------------------------------------
//     AOD track cov matrix base class
//-------------------------------------------------------------------------

ClassImp(AliAODNeutral::AliAODNeuCov)

//______________________________________________________________________________
template <class T> void AliAODNeutral::AliAODNeuCov::GetCovMatrix(T *cmat) const
{
  //
  // Returns the external cov matrix
  //
  cmat[ 0] = fDiag[ 0]*fDiag[ 0];
  cmat[ 2] = fDiag[ 1]*fDiag[ 1];
  cmat[ 5] = fDiag[ 2]*fDiag[ 2];
  cmat[ 9] = fDiag[ 3]*fDiag[ 3];
  //
  cmat[ 1] = fODia[ 0]*fDiag[ 0]*fDiag[ 1];
  cmat[ 3] = fODia[ 1]*fDiag[ 0]*fDiag[ 2];
  cmat[ 4] = fODia[ 2]*fDiag[ 1]*fDiag[ 2];
  cmat[ 6] = fODia[ 3]*fDiag[ 0]*fDiag[ 3];
  cmat[ 7] = fODia[ 4]*fDiag[ 1]*fDiag[ 3];
  cmat[ 8] = fODia[ 5]*fDiag[ 2]*fDiag[ 3];

}


//______________________________________________________________________________
template <class T> void AliAODNeutral::AliAODNeuCov::SetCovMatrix(T *cmat)
{
  //
  // Sets the external cov matrix
  //
  if(cmat) {
    fDiag[ 0] = TMath::Sqrt(cmat[ 0]);
    fDiag[ 1] = TMath::Sqrt(cmat[ 2]);
    fDiag[ 2] = TMath::Sqrt(cmat[ 5]);
    fDiag[ 3] = TMath::Sqrt(cmat[ 9]);

    //
    fODia[ 0] = cmat[ 1]/(fDiag[ 0]*fDiag[ 1]);
    fODia[ 1] = cmat[ 3]/(fDiag[ 0]*fDiag[ 2]);
    fODia[ 2] = cmat[ 4]/(fDiag[ 1]*fDiag[ 2]);
    fODia[ 3] = cmat[ 6]/(fDiag[ 0]*fDiag[ 3]);
    fODia[ 4] = cmat[ 7]/(fDiag[ 1]*fDiag[ 3]);
    fODia[ 5] = cmat[ 8]/(fDiag[ 2]*fDiag[ 3]);
  } else {
    for(Int_t i=0; i< 4; ++i) fDiag[i]=-999.;
    for(Int_t i=0; i<6; ++i) fODia[i]=0.;
  }
}
