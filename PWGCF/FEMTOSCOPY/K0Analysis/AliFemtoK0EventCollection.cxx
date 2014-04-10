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

////////////////////////////////////////////////////////////////////////////////
//
//  These classes provide storage for event and track information which 
//  are used for same-event and mixed-event analyses in AliFemtoK0Analysis. 
//
//  authors: Dhevan Gangadharan (dhevan.raja.gangadharan@cern.ch)
//           Matthew Steinpreis (matthew.steinpreis@cern.ch)
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoK0EventCollection.h"

AliFemtoK0Particle::AliFemtoK0Particle() :
 fMomentum(),
 fPt(0),
 fMass(0),
 fDaughterID1(0),
 fDaughterID2(0),
 fK0(0),
 fSideLeft(0),
 fSideRight(0),
 fSkipShared(0),
 fV0Dca(0),
 fDDDca(0),
 fDecayLength(0),
 fPosDca(0),
 fNegDca(0),
 fPosPt(0),
 fNegPt(0),
 fPosPhi(0),
 fNegPhi(0),
 fPPos(),
 fPNeg(),
 fPosXYZ(),
 fNegXYZ(),
 fPhi(0),
 fPhiPsi(0),
 fCutPass()
{
  //Default constructor
}
//_____________________________________________________________________________
AliFemtoK0Particle::AliFemtoK0Particle(const AliFemtoK0Particle &obj) :
 fMomentum(),
 fPt(obj.fPt),
 fMass(obj.fMass),
 fDaughterID1(obj.fDaughterID1),
 fDaughterID2(obj.fDaughterID2),
 fK0(obj.fK0),
 fSideLeft(obj.fSideLeft),
 fSideRight(obj.fSideRight),
 fSkipShared(obj.fSkipShared),
 fV0Dca(obj.fV0Dca),
 fDDDca(obj.fDDDca),
 fDecayLength(obj.fDecayLength),
 fPosDca(obj.fPosDca),
 fNegDca(obj.fNegDca),
 fPosPt(obj.fPosPt),
 fNegPt(obj.fNegPt),
 fPosPhi(obj.fPosPhi),
 fNegPhi(obj.fNegPhi),
 fPPos(),
 fPNeg(),
 fPosXYZ(),
 fNegXYZ(),
 fPhi(),
 fPhiPsi(),
 fCutPass()

{
  // copy constructor
}
//_____________________________________________________________________________
AliFemtoK0Particle &AliFemtoK0Particle::operator=(const AliFemtoK0Particle &obj)
{
 //Assignment operator
 if(this == &obj) return *this;

 fMomentum[0] = obj.fMomentum[0];
 fMomentum[1] = obj.fMomentum[1];
 fMomentum[2] = obj.fMomentum[2]; 
 fPt = obj.fPt;
 fMass = obj.fMass;
 fDaughterID1 = obj.fDaughterID1;
 fDaughterID2 = obj.fDaughterID2;
 fK0 = obj.fK0;
 fSideLeft = obj.fSideLeft;
 fSideRight = obj.fSideRight;
 fSkipShared = obj.fSkipShared;
 fV0Dca = obj.fV0Dca;
 fDDDca = obj.fDDDca;
 fDecayLength = obj.fDecayLength;
 fPosDca = obj.fPosDca;
 fNegDca = obj.fNegDca;
 fPosPt = obj.fPosPt;
 fNegPt = obj.fNegPt;
 fPosPhi = obj.fPosPhi;
 fNegPhi = obj.fNegPhi;
 for(int i=0;i<3;i++){
  fPPos[i] = obj.fPPos[i];
  fPNeg[i] = obj.fPNeg[i];
  for(int j=0;j<9;j++){
   fPosXYZ[j][i] = obj.fPosXYZ[j][i];
   fNegXYZ[j][i] = obj.fNegXYZ[j][i];
  }
 }
 for(int i=0;i<4;i++){
  for(int j=0;j<5;j++){
   fCutPass[i][j] = obj.fCutPass[i][j];
 }}
 fPhi = obj.fPhi;
 fPhiPsi = obj.fPhiPsi;
 return (*this);
}
//_____________________________________________________________________________
AliFemtoK0Particle::~AliFemtoK0Particle()
{
  // Destructor
}
//_____________________________________________________________________________

AliFemtoK0Event::AliFemtoK0Event():
 fFillStatus(0),
 fNumV0s(0),
 fK0Particle(0x0)
{
  //Default constructor
}
//_____________________________________________________________________________
AliFemtoK0Event::AliFemtoK0Event(const AliFemtoK0Event &obj) : 
 fFillStatus(obj.fFillStatus),
 fNumV0s(obj.fNumV0s),
 fK0Particle(obj.fK0Particle)
{
  //Copy constructor
}
//_____________________________________________________________________________
AliFemtoK0Event &AliFemtoK0Event::operator=(const AliFemtoK0Event &obj)
{
 //Assignment operator
 if(this == &obj) return *this;

 fFillStatus = obj.fFillStatus;
 fNumV0s = obj.fNumV0s;
 fK0Particle = obj.fK0Particle;

 return (*this);
}
//_____________________________________________________________________________
AliFemtoK0Event::~AliFemtoK0Event()
{
 //Destructor
 if(fK0Particle){
  delete fK0Particle;
  fK0Particle = NULL;
 }
}
//_____________________________________________________________________________
AliFemtoK0EventCollection::AliFemtoK0EventCollection() : 
 fBufferSize(0),
 fLimit(0),
 fEvt(0x0)
{
 //Default constructor
}
//_____________________________________________________________________________
AliFemtoK0EventCollection::AliFemtoK0EventCollection(short a, int lim) :
  fBufferSize(0),
  fLimit(0),
  fEvt(0x0)
{
  //maint constructor
 
  SetBufferSize(a);
  
  fEvt = new AliFemtoK0Event[fBufferSize];  //allocate pointer array of type particle_event
  fLimit = lim;

  for(int ii = 0; ii < fBufferSize; ii++){   //Initialize particle table pointers to NULL
    (fEvt + ii)->fK0Particle = NULL;
    (fEvt + ii)->fNumV0s = 0;
    (fEvt + ii)->fFillStatus = 0;

    (fEvt + ii)->fK0Particle = new AliFemtoK0Particle[fLimit];
  }
}
//_____________________________________________________________________________
AliFemtoK0EventCollection::AliFemtoK0EventCollection(const AliFemtoK0EventCollection &obj)
 : fBufferSize(obj.fBufferSize),
   fLimit(obj.fLimit),
   fEvt(obj.fEvt)
{
 //copy constructor
}
//_____________________________________________________________________________
AliFemtoK0EventCollection &AliFemtoK0EventCollection::operator=(const
AliFemtoK0EventCollection &obj)
{
 //Assignment operator
 if(this == &obj) return *this;
 
 fBufferSize = obj.fBufferSize;
 fLimit = obj.fLimit;
 fEvt = obj.fEvt;

 return (*this);
}
//_____________________________________________________________________________
AliFemtoK0EventCollection::~AliFemtoK0EventCollection()
{
 for(int i =0; i < fBufferSize; i++){
  if((fEvt + i)->fK0Particle != NULL){
   delete [] (fEvt + i)->fK0Particle;
   (fEvt + i)->fK0Particle = NULL;
  }
 }
 
 delete [] fEvt; fEvt = NULL;
}
//_____________________________________________________________________________
void AliFemtoK0EventCollection::FIFOShift(){ //Shift elements in FIFO by one and clear last element in FIFO 
  
  for(unsigned short i=fBufferSize-1 ; i > 0; i--){
    for(int j=0; j<(fEvt + i-1)->fNumV0s; j++) (fEvt + i)->fK0Particle[j] = (fEvt + i-1)->fK0Particle[j];
    (fEvt + i)->fFillStatus = (fEvt + i-1)->fFillStatus;
    (fEvt + i)->fNumV0s = (fEvt + i-1)->fNumV0s;
  }

  (fEvt)->fNumV0s=0;
  (fEvt)->fFillStatus=0;

}
