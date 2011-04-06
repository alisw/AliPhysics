#ifndef ALIRSNLOOPEFF_H
#define ALIRSNLOOPEFF_H

//
// Class to combine pairs of daughters.
//

#include "AliLog.h"

#include "AliRsnEvent.h"
#include "AliRsnLoop.h"

class AliRsnLoopEff : public AliRsnLoop {
public:

   AliRsnLoopEff(const char *name = "default", Int_t nSteps = 0);
   AliRsnLoopEff(const AliRsnLoopEff &copy);
   AliRsnLoopEff& operator=(const AliRsnLoopEff&);
   ~AliRsnLoopEff();

   AliRsnListOutput* GetOutput()    {return (AliRsnListOutput*)fOutputs[0];}
   void              CreateOutput();

   void              AddStep(TObject *set);
   Int_t             NStepsArray() {return (Int_t)fSteps.GetEntries();}
   Int_t             NStepsAll()   {return fAddSteps + NStepsArray();}
                     
   virtual void      AddOutput(TObject *) { AliWarning("In loops for efficiency it is not allowed to add outputs externally"); }
   virtual Bool_t    Init(const char *prefix, TList *list);

protected:

   Int_t  FindTrack(Int_t label, AliVEvent *event);

   Int_t        fAddSteps; //  number of additional steps
   TObjArray    fSteps;    //  list of cuts for all steps with MC tracks

private:

   ClassDef(AliRsnLoopEff, 1)
};

#endif

