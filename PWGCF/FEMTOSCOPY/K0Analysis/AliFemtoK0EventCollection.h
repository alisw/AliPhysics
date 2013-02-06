#ifndef ALIFEMTOK0EVENTCOLLECTION_H
#define ALIFEMTOK0EVENTCOLLECTION_H
//
//Class AliFemtoK0Particle, AliFemtoK0Event, AliFemtoK0EventCollection
//
//AliFemtoK0Particle, AliFemtoK0Event, AliFemtoK0EventCollection
//authors: 
//        Dhevan Gangadharan (dhevan.raja.gangadharan@cern.ch)
//        Matthew Steinpreis (matthew.steinpreis@cern.ch)
//


#include <iostream>
#include <string>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TBits.h"
#include "TObject.h"
#include "TVector2.h"
#include "AliESDtrack.h"

using namespace std;

class AliFemtoK0Particle  // Reconstructed K0s parameters needed for correlations
{
 public:
  
  AliFemtoK0Particle();
  virtual ~AliFemtoK0Particle();
  AliFemtoK0Particle(const AliFemtoK0Particle &obj);
  AliFemtoK0Particle &operator=(const AliFemtoK0Particle &obj);
 
  double fMomentum[3];  //v0 momentum
  double fPt;           //v0 transverse momentum
  double fMass;         //v0 reconstructed mass
  short fDaughterID1;   //Daughter (pion) AODtrack ID
  short fDaughterID2;   //Daughter (pion) AODtrack ID
  bool fK0;             //if v0 has "good" K0 mass
  bool fSideLeft;
  bool fSideRight;
  bool fSkipShared;
  double fV0Dca;	//used in v0 selection process

  //for single particle histograms
  double fDDDca;	//daughter-daughter DCA
  double fDecayLength; 	//v0 decay length
  double fPosDca;	//positive daughter Dca to prim vert
  double fNegDca;	//negative ""
  double fPosPt; 	//positive daughter pt
  double fNegPt;	//negative daughter pt
  double fPosPhi;	//positive daughter phi
  double fNegPhi;	//negative daughter phi

  //for separation
  double fXPos[3];      //Positive daughter position
  double fXNeg[3];      //Negative daughter position
  double fPPos[3];      //Positive daughter momentum
  double fPNeg[3];      //negative daughter momentum
  double fCovPos[21];   //positive daughter coverity matrix
  double fCovNeg[21];   //negative daughter coverity matrix

  ClassDef(AliFemtoK0Particle, 1);
};

class AliFemtoK0Event // like particle_event
{
 public:

  AliFemtoK0Event();
  virtual ~AliFemtoK0Event();
  AliFemtoK0Event(const AliFemtoK0Event &obj);
  AliFemtoK0Event &operator=(const AliFemtoK0Event &obj);
  
  int fFillStatus;     //tells AliFemtoK0EventCollection to add event
  int fNumV0s;         //number of collected v0s in event
  AliFemtoK0Particle *fK0Particle; //class for K0 parameters needed for CF

  ClassDef(AliFemtoK0Event, 1);
};

class AliFemtoK0EventCollection 
{
  public:
    AliFemtoK0EventCollection();
    AliFemtoK0EventCollection(short,int);
    virtual ~AliFemtoK0EventCollection();
    AliFemtoK0EventCollection(const AliFemtoK0EventCollection &obj);
    AliFemtoK0EventCollection &operator=(const AliFemtoK0EventCollection &obj);

    short fBufferSize; //Size of the Event Storage buffer
    int fLimit;        //Max number of tracks
    AliFemtoK0Event *fEvt; //event class

    void FIFOShift();  //remove/add event (first in, first out)
    void SetBufferSize(short a){fBufferSize = a;} //set size of event buffer
 
    ClassDef(AliFemtoK0EventCollection, 1);
};
#endif

