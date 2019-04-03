#ifndef ALIFEMTOXiParticle_H
#define ALIFEMTOXiParticle_H
//
//Class AliFemtoLambdaParticle, AliFemtoLambdaEvent, AliFemtoLambdaEventCollection
//
//AliFemtoLambdaParticle, AliFemtoLambdaEvent, AliFemtoLambdaEventCollection
//


#include <iostream>
#include <string>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TBits.h"
#include "TObject.h"
#include "AliESDtrack.h"
#include <TVector3.h>

using namespace std;

class AliFemtoXiParticle // Reconstructed Lambdas parameters needed for correlations
{
 public:
  
  AliFemtoXiParticle();
  virtual ~AliFemtoXiParticle();
  //AliFemtoXiParticle(const AliFemtoXiParticle &obj);
  AliFemtoXiParticle &operator=(const AliFemtoXiParticle &obj);
  
  TVector3 fMomentum;  //Xi momentum
  float fPt;           //Xi transverse momentum
  float fMass;         //Xi reconstructed mass
  short fDaughterID1;   //Daughter (proton) AODtrack ID
  short fDaughterID2;   //Daughter (pion) AODtrack ID
  short fBachID; //Bachelor ID of Xi candidate
  bool fXitag;
  float fPointing;

  ClassDef(AliFemtoXiParticle, 3)
};
#endif
