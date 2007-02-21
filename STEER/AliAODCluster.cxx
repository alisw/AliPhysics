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
//     AOD cluster base class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include "AliAODCluster.h"

ClassImp(AliAODCluster)

//______________________________________________________________________________
AliAODCluster::AliAODCluster() : 
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
AliAODCluster::AliAODCluster(Int_t id,
			     Int_t label, 
			     Double_t energy,
			     Double_t x[3],
			     Double_t covMatrix[10],
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
 
  SetPosition(x);
  if(covMatrix) SetCovMatrix(covMatrix);
  SetPID(pid);

}

//______________________________________________________________________________
AliAODCluster::AliAODCluster(Int_t id,
			     Int_t label, 
			     Float_t energy,
			     Float_t x[3],
			     Float_t covMatrix[10],
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
 
  SetPosition(x);
  if(covMatrix) SetCovMatrix(covMatrix);
  SetPID(pid);

}


//______________________________________________________________________________
AliAODCluster::~AliAODCluster() 
{
  // destructor
  delete fCovMatrix;
}


//______________________________________________________________________________
AliAODCluster::AliAODCluster(const AliAODCluster& trk) :
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
  if(trk.fCovMatrix) fCovMatrix=new AliAODRedCov<4>(*trk.fCovMatrix);
  SetPID(trk.fPID);

}

//______________________________________________________________________________
AliAODCluster& AliAODCluster::operator=(const AliAODCluster& trk)
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
    if(trk.fCovMatrix) fCovMatrix=new AliAODRedCov<4>(*trk.fCovMatrix);
    else fCovMatrix=NULL;
    fProdVertex = trk.fProdVertex;
    fPrimTrack = trk.fPrimTrack;

    fType = trk.fType;
  }

  return *this;
}

//______________________________________________________________________________
template <class T> void AliAODCluster::SetPosition(const T *x) 
{
  // set the position

  if (x) {
      fPosition[0] = x[0];
      fPosition[1] = x[1];
      fPosition[2] = x[2];
  } else {

    fPosition[0] = -999.;
    fPosition[1] = -999.;
    fPosition[2] = -999.;
  }
}

//______________________________________________________________________________
void AliAODCluster::Print(Option_t* /* option */) const
{
  // prints information about AliAODCluster

  printf("Object name: %s   Cluster type: %s\n", GetName(), GetTitle()); 
  printf("    energy = %f\n", E());
  printf("      chi2 = %f\n", Chi2());
  printf(" PID object: %p\n", PID());
}

