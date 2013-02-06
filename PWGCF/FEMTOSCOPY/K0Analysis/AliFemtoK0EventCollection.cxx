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
 fXPos(),
 fXNeg(),
 fPPos(),
 fPNeg(),
 fCovPos(),
 fCovNeg()
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
 fXPos(),
 fXNeg(),
 fPPos(),
 fPNeg(),
 fCovPos(),
 fCovNeg()
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
 fXPos[0] = obj.fXPos[0];
 fXPos[1] = obj.fXPos[1];
 fXPos[2] = obj.fXPos[2];
 fXNeg[0] = obj.fXNeg[0];
 fXNeg[1] = obj.fXNeg[1];
 fXNeg[2] = obj.fXNeg[2];
 fPPos[0] = obj.fPPos[0];
 fPPos[1] = obj.fPPos[1];
 fPPos[2] = obj.fPPos[2];
 fPNeg[0] = obj.fPNeg[0];
 fPNeg[1] = obj.fPNeg[1];
 fPNeg[2] = obj.fPNeg[2];
 fCovPos[0] = obj.fCovPos[0];  fCovPos[1] = obj.fCovPos[1];  fCovPos[2] = obj.fCovPos[2];
 fCovPos[3] = obj.fCovPos[3];  fCovPos[4] = obj.fCovPos[4];  fCovPos[5] = obj.fCovPos[5];
 fCovPos[6] = obj.fCovPos[6];  fCovPos[7] = obj.fCovPos[7];  fCovPos[8] = obj.fCovPos[8];
 fCovPos[9] = obj.fCovPos[9];  fCovPos[10] = obj.fCovPos[10];  fCovPos[11] = obj.fCovPos[11];
 fCovPos[12] = obj.fCovPos[12];  fCovPos[13] = obj.fCovPos[13];  fCovPos[14] = obj.fCovPos[14];
 fCovPos[17] = obj.fCovPos[17];  fCovPos[17] = obj.fCovPos[17];  fCovPos[17] = obj.fCovPos[17];
 fCovPos[18] = obj.fCovPos[18];  fCovPos[19] = obj.fCovPos[19];  fCovPos[20] = obj.fCovPos[20];
 fCovNeg[0] = obj.fCovNeg[0];  fCovNeg[1] = obj.fCovNeg[1];  fCovNeg[2] = obj.fCovNeg[2];
 fCovNeg[3] = obj.fCovNeg[3];  fCovNeg[4] = obj.fCovNeg[4];  fCovNeg[5] = obj.fCovNeg[5];
 fCovNeg[6] = obj.fCovNeg[6];  fCovNeg[7] = obj.fCovNeg[7];  fCovNeg[8] = obj.fCovNeg[8];
 fCovNeg[9] = obj.fCovNeg[9];  fCovNeg[10] = obj.fCovNeg[10];  fCovNeg[11] = obj.fCovNeg[11];
 fCovNeg[12] = obj.fCovNeg[12];  fCovNeg[13] = obj.fCovNeg[13];  fCovNeg[14] = obj.fCovNeg[14];
 fCovNeg[17] = obj.fCovNeg[17];  fCovNeg[17] = obj.fCovNeg[17];  fCovNeg[17] = obj.fCovNeg[17];
 fCovNeg[18] = obj.fCovNeg[18];  fCovNeg[19] = obj.fCovNeg[19];  fCovNeg[20] = obj.fCovNeg[20];

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
 if(fK0Particle) delete fK0Particle;
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
  }
 }
 
 delete [] fEvt;
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

