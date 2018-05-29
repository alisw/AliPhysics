#ifndef ALIFEMTOLambdaEVENTCOLLECTION_H
#define ALIFEMTOLambdaEVENTCOLLECTION_H
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

#include "AliFemtoLambdaEvent.h"
#include "AliFemtoLambdaParticle.h"

using namespace std;

class AliFemtoLambdaEventCollection2
{
  public:
    AliFemtoLambdaEventCollection2();
    AliFemtoLambdaEventCollection2(short,int);
    virtual ~AliFemtoLambdaEventCollection2();
    AliFemtoLambdaEventCollection2(const AliFemtoLambdaEventCollection2 &obj);
    AliFemtoLambdaEventCollection2 &operator=(const AliFemtoLambdaEventCollection2 &obj);

    short fBufferSize; //Size of the Event Storage buffer
    int fLimit;        //Max number of tracks
    AliFemtoLambdaEvent *fEvt; //event class
    int fNumEvents;

    void FIFOShift();  //remove/add event (first in, first out)
    void SetBufferSize(short a){fBufferSize = a;} //set size of event buffer

    ClassDef(AliFemtoLambdaEventCollection2, 2)
};
#endif
