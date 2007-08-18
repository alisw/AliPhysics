#ifndef ALIACORDETrigger_H
#define ALIACORDETrigger_H

///_________________________________________________________________________
///
///  Class for making  ACORDE Trigger
///_________________________________________________________________________   


#include "AliTriggerDetector.h"
#include "AliTriggerInput.h"

#include "AliACORDELoader.h"
#include "AliACORDEdigit.h"

#include "AliLog.h"


class AliACORDETrigger : public AliTriggerDetector
{
 public:
                   AliACORDETrigger();   // constructor
   virtual        ~AliACORDETrigger(){}  // destructor
   virtual void    CreateInputs();
   virtual void    Trigger();

   virtual Int_t   GetSingleMuon() const {return fSingleMuon;}
   virtual Int_t   GetMultiMuon() const {return fMultiMuon;}
   virtual Bool_t  GetModuleFired(Int_t i) const {return fModuleFired[i-1];}

private:

   Int_t fSingleMuon; // number of module firing the Single Muon trigger
   Int_t fMultiMuon;  // number of modules firing for the Multi Muon trigger
   Bool_t fModuleFired[60]; // modules which have fired

   ClassDef( AliACORDETrigger, 1 )  // ACORDE Trigger Detector class
};

#endif // AliACORDETrigger_H
