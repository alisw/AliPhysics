#ifndef ALIFEMTOLambdaEVENT_H
#define ALIFEMTOLambdaEVENT_H
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

#include "AliFemtoProtonParticle.h"
#include "AliFemtoLambdaParticle.h"
#include "AliFemtoXiParticle.h"

using namespace std;

class AliFemtoLambdaEvent // like particle_event
{
 public:

  AliFemtoLambdaEvent();
  virtual ~AliFemtoLambdaEvent();
  AliFemtoLambdaEvent(const AliFemtoLambdaEvent &obj);
  AliFemtoLambdaEvent &operator=(const AliFemtoLambdaEvent &obj);

  int fEventNumber;
  int fFillStatus;     //tells AliFemtoLambdaEventCollection to add event
  int fNumV0s;         //number of collected v0s in event
  int fNumAntiV0s;         //number of collected Anti-v0s in event
  int fNumProtons;     //number of primary protons used for correlations
  int fNumAntiProtons;     //number of primary anti-protons used for correlations
  int fNumXis;
  int fNumV0SideBand_left; //number of fake V0s on the left side of Lambda peak
  int fNumV0SideBand_right; //number of fake V0s on the right side of Lambda peak
  int fMultBin;
  int fBoostVal;
  double fSphericity;
  AliFemtoLambdaParticle *fLambdaParticle; //class for Lambda parameters needed for CF
  AliFemtoLambdaParticle *fAntiLambdaParticle; //class for Lambda parameters needed for CF
  AliFemtoProtonParticle *fProtonParticle;
  AliFemtoProtonParticle *fAntiProtonParticle;
  AliFemtoLambdaParticle *fLambdaSideBand_left;
  AliFemtoLambdaParticle *fLambdaSideBand_right;
  AliFemtoXiParticle *fXiParticle;
#ifdef __ROOT__
  ClassDef(AliFemtoLambdaEvent, 1);
#endif
};

#endif
