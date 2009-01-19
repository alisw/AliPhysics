/**************************************************************************
 * Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// Base class for AOD reconstructed heavy-flavour decay
//
// Author: A.Dainese, andrea.dainese@lnl.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <TVector3.h>
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"

ClassImp(AliAODRecoDecayHF)

//--------------------------------------------------------------------------
AliAODRecoDecayHF::AliAODRecoDecayHF() :
  AliAODRecoDecay(),
  fOwnPrimaryVtx(0x0),
  fd0err(0x0), 
  fProngID(0x0) 
{
  //
  // Default Constructor
  //
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF::AliAODRecoDecayHF(AliAODVertex *vtx2,Int_t nprongs,Short_t charge,
				     Double_t *px,Double_t *py,Double_t *pz,
				     Double_t *d0,Double_t *d0err) :
  AliAODRecoDecay(vtx2,nprongs,charge,px,py,pz,d0),
  fOwnPrimaryVtx(0x0),
  fd0err(0x0),
  fProngID(0x0) 
{
  //
  // Constructor with AliAODVertex for decay vertex
  //
  fd0err = new Double_t[GetNProngs()];
  for(Int_t i=0; i<GetNProngs(); i++) fd0err[i] = d0err[i];
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF::AliAODRecoDecayHF(AliAODVertex *vtx2,Int_t nprongs,Short_t charge,
				     Double_t *d0,Double_t *d0err) :
  AliAODRecoDecay(vtx2,nprongs,charge,d0),
  fOwnPrimaryVtx(0x0),
  fd0err(0x0),
  fProngID(0x0) 
{
  //
  // Constructor with AliAODVertex for decay vertex and without prongs momenta
  //
  fd0err = new Double_t[GetNProngs()];
  for(Int_t i=0; i<GetNProngs(); i++) fd0err[i] = d0err[i];
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF::AliAODRecoDecayHF(Double_t vtx1[3],Double_t vtx2[3],
				     Int_t nprongs,Short_t charge,
				     Double_t *px,Double_t *py,Double_t *pz,
				     Double_t *d0) :
  AliAODRecoDecay(0x0,nprongs,charge,px,py,pz,d0),
  fOwnPrimaryVtx(0x0),
  fd0err(0x0),
  fProngID(0x0) 
{
  //
  // Constructor that can used for a "MC" object
  //

  fOwnPrimaryVtx = new AliAODVertex(vtx1);

  AliAODVertex *vtx = new AliAODVertex(vtx2);
  SetOwnSecondaryVtx(vtx);

}
//--------------------------------------------------------------------------
AliAODRecoDecayHF::AliAODRecoDecayHF(const AliAODRecoDecayHF &source) :
  AliAODRecoDecay(source),
  fOwnPrimaryVtx(source.fOwnPrimaryVtx),
  fd0err(0x0),
  fProngID(0x0)
{
  //
  // Copy constructor
  //
  if(source.GetNProngs()>0) {
    fd0err = new Double_t[GetNProngs()];
    memcpy(fd0err,source.fd0err,GetNProngs()*sizeof(Double_t));
    if(source.fProngID) {
      fProngID = new UShort_t[GetNProngs()];
      memcpy(fProngID,source.fProngID,GetNProngs()*sizeof(UShort_t));
    }
  }
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF &AliAODRecoDecayHF::operator=(const AliAODRecoDecayHF &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliAODRecoDecay::operator=(source);

  fOwnPrimaryVtx = source.fOwnPrimaryVtx;
  if(source.GetNProngs()>0) {
    fd0err = new Double_t[GetNProngs()];
    memcpy(fd0err,source.fd0err,GetNProngs()*sizeof(Double_t));
    if(source.fProngID) {
      fProngID = new UShort_t[GetNProngs()];
      memcpy(fProngID,source.fProngID,GetNProngs()*sizeof(UShort_t));
    }
  }
  return *this;
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF::~AliAODRecoDecayHF() {
  //  
  // Default Destructor
  //
  if(fOwnPrimaryVtx) delete fOwnPrimaryVtx;
  if(fd0err) delete [] fd0err;
  if(fProngID) delete [] fProngID;
}
//---------------------------------------------------------------------------
