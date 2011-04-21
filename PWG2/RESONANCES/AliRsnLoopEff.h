#ifndef ALIRSNLOOPEFF_H
#define ALIRSNLOOPEFF_H

//
// Class to combine pairs of daughters.
//

#include "AliLog.h"

#include "AliRsnLoop.h"
#include "AliRsnListOutput.h"

class AliRsnLoopEff : public AliRsnLoop {
public:

   AliRsnLoopEff(const char *name = "default", Int_t nSteps = 0, Double_t maxDistPV = 1E-2);
   AliRsnLoopEff(const AliRsnLoopEff &copy);
   AliRsnLoopEff& operator=(const AliRsnLoopEff&);
   ~AliRsnLoopEff();

   AliRsnListOutput* GetOutput()          {return fOutput;}
   void              CreateOutput();
   void              AddStep(TObject *set);                     
   void              SetMaxDistanceFromPV(Double_t value) {fMaxDistPV = value;}
   virtual void      AddOutput(TObject *) {AliWarning("In loops for efficiency it is not allowed to add outputs externally");}
   virtual Bool_t    Init(const char *prefix, TList *list);

protected:

   Int_t             FindTrack(Int_t label, AliVEvent *event);
   Int_t             GetMatchedDaughter(Int_t label, AliRsnEvent *event);
   Double_t          DistanceFromPV(Double_t x, Double_t y, Double_t z);
   Bool_t            CheckDistanceFromPV(Double_t x, Double_t y, Double_t z) {return (DistanceFromPV(x,y,z) <= fMaxDistPV);}

   Int_t              fAddSteps;  //  number of additional steps
   TObjArray          fSteps;     //  list of cuts for all steps with MC tracks
   AliRsnListOutput  *fOutput;    //  unique list output
   Double_t           fVertex[3]; //! primary vertex position
   Double_t           fMaxDistPV; //  maximum allowed distance from primary vertex

private:

   ClassDef(AliRsnLoopEff, 1)
};

#endif

