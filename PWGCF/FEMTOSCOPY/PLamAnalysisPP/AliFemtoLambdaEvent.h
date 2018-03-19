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
  int fNumV0s;         //number of collected v0s in event
  int fNumAntiV0s;         //number of collected Anti-v0s in event
  int fNumProtons;     //number of primary protons used for correlations
  int fNumAntiProtons;     //number of primary anti-protons used for correlations
  int fNumXis;

  AliFemtoLambdaParticle *fLambdaParticle; //class for Lambda parameters needed for CF
  AliFemtoLambdaParticle *fAntiLambdaParticle; //class for Lambda parameters needed for CF
  AliFemtoProtonParticle *fProtonParticle;
  AliFemtoProtonParticle *fAntiProtonParticle;
  AliFemtoXiParticle *fXiParticle;

  ClassDef(AliFemtoLambdaEvent, 2)
};

#endif
