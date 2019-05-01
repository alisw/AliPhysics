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

#include "AliFemtoLambdaEventCollection2.h"

ClassImp(AliFemtoLambdaEventCollection2)

//_____________________________________________________________________________
AliFemtoLambdaEventCollection2::AliFemtoLambdaEventCollection2() : 
 fBufferSize(0),
 fLimit(0),
 fEvt(0x0),
 fNumEvents(0)
{
 //Default constructor
}
//_____________________________________________________________________________
AliFemtoLambdaEventCollection2::AliFemtoLambdaEventCollection2(short a, int lim) :
  fBufferSize(0),
  fLimit(0),
  fEvt(0x0),
  fNumEvents(0)
{
  fBufferSize = a;
  fLimit = lim;
  
  //SetBufferSize(a);
  
  fEvt = new AliFemtoLambdaEvent[fBufferSize];  //allocate pointer array of type particle_event
  

  for(int ii = 0; ii < fBufferSize; ii++)   //Initialize particle table pointers to NULL
    {
      (fEvt + ii)->fLambdaParticle = 0;
      (fEvt + ii)->fAntiLambdaParticle = 0;
      (fEvt + ii)->fProtonParticle = 0;
      (fEvt + ii)->fAntiProtonParticle = 0;
      (fEvt + ii)->fXiParticle = 0;
      
      (fEvt + ii)->fNumV0s = 0;
      (fEvt + ii)->fNumAntiV0s = 0;
      (fEvt + ii)->fNumProtons = 0;
      (fEvt + ii)->fNumAntiProtons = 0;
      (fEvt + ii)->fNumXis = 0;
      
      (fEvt + ii)->fLambdaParticle = new AliFemtoLambdaParticle[fLimit];
      (fEvt + ii)->fAntiLambdaParticle = new AliFemtoLambdaParticle[fLimit];
      (fEvt + ii)->fProtonParticle = new AliFemtoProtonParticle[fLimit];
      (fEvt + ii)->fAntiProtonParticle = new AliFemtoProtonParticle[fLimit];
      (fEvt + ii)->fXiParticle = new AliFemtoXiParticle[fLimit];
    }
}
//_____________________________________________________________________________
AliFemtoLambdaEventCollection2::AliFemtoLambdaEventCollection2(const AliFemtoLambdaEventCollection2 &obj)
 : fBufferSize(obj.fBufferSize),
   fLimit(obj.fLimit),
   fEvt(obj.fEvt),
   fNumEvents(obj.fNumEvents)
{
 //copy constructor
}
//_____________________________________________________________________________
AliFemtoLambdaEventCollection2 &AliFemtoLambdaEventCollection2::operator=(const
AliFemtoLambdaEventCollection2 &obj)
{
 //Assignment operator
 if(this == &obj) return *this;
 
 fBufferSize = obj.fBufferSize;
 fLimit = obj.fLimit;
 fEvt = obj.fEvt;
 fNumEvents = obj.fNumEvents;


 return (*this);
}
//_____________________________________________________________________________
AliFemtoLambdaEventCollection2::~AliFemtoLambdaEventCollection2()
{
 for(int i =0; i < fBufferSize; i++)
   {
     if((fEvt + i)->fLambdaParticle != 0)
       {
	 delete [] (fEvt + i)->fLambdaParticle;
	 (fEvt + i)->fLambdaParticle = 0;
       }
     if((fEvt + i)->fAntiLambdaParticle != 0)
       {
	 delete [] (fEvt + i)->fAntiLambdaParticle;
	 (fEvt + i)->fAntiLambdaParticle = 0;
       }
     if((fEvt + i)->fProtonParticle != 0)
       {
	 delete [] (fEvt + i)->fProtonParticle;
	 (fEvt + i)->fProtonParticle = 0;
       }
     if((fEvt + i)->fAntiProtonParticle != 0)
       {
	 delete [] (fEvt + i)->fAntiProtonParticle;
	 (fEvt + i)->fAntiProtonParticle = 0;
       }
     if((fEvt + i)->fXiParticle != 0)
       {
	 delete [] (fEvt + i)->fXiParticle;
	 (fEvt + i)->fXiParticle = 0;
       }

   }
 
 delete [] fEvt; fEvt = 0;
}
//_____________________________________________________________________________
void AliFemtoLambdaEventCollection2::FIFOShift(){ //Shift elements in FIFO by one and clear last element in FIFO 
  
  for(unsigned short i=fBufferSize-1 ; i > 0; i--)
    {
      for(int j=0; j<(fEvt + i-1)->fNumV0s; j++) (fEvt + i)->fLambdaParticle[j] = (fEvt + i-1)->fLambdaParticle[j];
      for(int j=0; j<(fEvt + i-1)->fNumAntiV0s; j++) (fEvt + i)->fAntiLambdaParticle[j] = (fEvt + i-1)->fAntiLambdaParticle[j];
      for(int j=0; j<(fEvt + i-1)->fNumProtons; j++) (fEvt + i)->fProtonParticle[j] = (fEvt + i-1)->fProtonParticle[j];
      for(int j=0; j<(fEvt + i-1)->fNumAntiProtons; j++) (fEvt + i)->fAntiProtonParticle[j] = (fEvt + i-1)->fAntiProtonParticle[j];
      for(int j=0; j<(fEvt + i-1)->fNumXis; j++) (fEvt + i)->fXiParticle[j] = (fEvt + i-1)->fXiParticle[j];
      
      (fEvt + i)->fNumV0s = (fEvt + i-1)->fNumV0s;
      (fEvt + i)->fNumAntiV0s = (fEvt + i-1)->fNumAntiV0s;
      (fEvt + i)->fNumProtons = (fEvt + i-1)->fNumProtons;
      (fEvt + i)->fNumAntiProtons = (fEvt + i-1)->fNumAntiProtons;
      (fEvt + i)->fNumXis = (fEvt + i-1)->fNumXis;
      (fEvt + i)->fEventNumber = (fEvt + i-1)->fEventNumber;
    }

  (fEvt)->fNumV0s = 0;
  (fEvt)->fNumAntiV0s = 0;
  (fEvt)->fNumProtons = 0;
  (fEvt)->fNumAntiProtons = 0;
  (fEvt)->fNumXis = 0;
  (fEvt)->fEventNumber = 0;
}
