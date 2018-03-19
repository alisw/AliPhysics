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

#include "AliFemtoLambdaEvent.h"
#include "AliFemtoProtonParticle.h"
#include "AliFemtoXiParticle.h"

AliFemtoLambdaEvent::AliFemtoLambdaEvent():
  fEventNumber(0),
  fNumV0s(0),
  fNumAntiV0s(0),
  fNumProtons(0),
  fNumAntiProtons(0),
  fNumXis(0),
  fLambdaParticle(0x0),
  fAntiLambdaParticle(0x0),
  fProtonParticle(0x0),
  fAntiProtonParticle(0x0),
  fXiParticle(0x0)
{
  //Default constructor
}
//_____________________________________________________________________________
AliFemtoLambdaEvent::AliFemtoLambdaEvent(const AliFemtoLambdaEvent &obj) :
  fEventNumber(obj.fEventNumber),
  fNumV0s(obj.fNumV0s),
  fNumAntiV0s(obj.fNumAntiV0s),
  fNumProtons(obj.fNumProtons),
  fNumAntiProtons(obj.fNumAntiProtons),
  fNumXis(obj.fNumXis),
  fLambdaParticle(obj.fLambdaParticle),
  fAntiLambdaParticle(obj.fAntiLambdaParticle),
  fProtonParticle(obj.fProtonParticle),
  fAntiProtonParticle(obj.fProtonParticle),
  fXiParticle(obj.fXiParticle)
{
  //Copy constructor
}
//_____________________________________________________________________________
AliFemtoLambdaEvent &AliFemtoLambdaEvent::operator=(const AliFemtoLambdaEvent &obj)
{
 //Assignment operator
 if(this == &obj) return *this;

 fEventNumber = obj.fEventNumber;
 fNumV0s = obj.fNumV0s;
 fNumAntiV0s = obj.fNumAntiV0s;
 fNumProtons = obj.fNumProtons;
 fNumAntiProtons = obj.fNumAntiProtons;
 fNumXis = obj.fNumXis;
 fLambdaParticle = obj.fLambdaParticle;
 fAntiLambdaParticle = obj.fAntiLambdaParticle;
 fProtonParticle = obj.fProtonParticle;
 fAntiProtonParticle = obj.fAntiProtonParticle;
 fXiParticle = obj.fXiParticle;

 return (*this);
}
//_____________________________________________________________________________
AliFemtoLambdaEvent::~AliFemtoLambdaEvent()
{
 //Destructor
 if(fLambdaParticle)
   {
     delete fLambdaParticle;
     fLambdaParticle = 0;
   }
 if(fAntiLambdaParticle)
    {
     delete fAntiLambdaParticle;
     fAntiLambdaParticle = 0;
    }
 if(fProtonParticle)
   {
     delete fProtonParticle;
     fProtonParticle = 0;
   }
 if(fAntiProtonParticle)
   {
     delete fAntiProtonParticle;
     fAntiProtonParticle = 0;
   }
 if(fXiParticle)
   {
     delete fXiParticle;
     fXiParticle = 0;
   }
}
//_____________________________________________________________________________
