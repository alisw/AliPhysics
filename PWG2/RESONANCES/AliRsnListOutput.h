#ifndef ALIRSNLISTOUTPUT_H
#define ALIRSNLISTOUTPUT_H

//
// General class for outputs which can stay into a TList
//

#include <TRef.h>
#include <TNamed.h>
#include <TArrayI.h>
#include <TObjArray.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TNtuple.h>
#include <THnSparse.h>

#include "AliCFContainer.h"

class AliRsnValue;

class AliRsnListOutput : public TNamed {

public:

   enum EOut {
      kHistoDefault,
      kHistoSparse,
      kCFContainer
   };

   AliRsnListOutput(const char *name = "", EOut type = kHistoDefault);
   AliRsnListOutput(const AliRsnListOutput &copy);
   const AliRsnListOutput& operator=(const AliRsnListOutput &copy);
   virtual ~AliRsnListOutput();

   EOut            GetType()               {return  fType;}
   Int_t           GetSteps()              {return  fSteps;}
   TObjArray*      GetValues()             {return &fValues;}
   Int_t           GetNValues()            {return (fNValues = fValues.GetEntries());}
   AliRsnValue*    GetValue(Int_t i)       {return (AliRsnValue*)fValues[i];}
   Int_t           GetIndex()              {return  fIndex;}
   void            SetType(EOut type)      {fType = type;}
   void            SetSteps(Int_t n)       {fSteps = n;}
   void            SetSkipFailed(Bool_t y) {fSkipFailed = y;}

   void            AddValue(AliRsnValue *value);

   virtual void    Reset();
   virtual Bool_t  Init(const char *prefix, TList *list);
   virtual Bool_t  Fill(TObject *target, Int_t step = 0);

private:

   TH1*            CreateHistogram(const char *name);
   THnSparseF*     CreateHistogramSparse(const char *name);
   AliCFContainer* CreateCFContainer(const char *name);

   Bool_t           fSkipFailed;    //  tell to skip fills when one computation fails
   EOut             fType;          //  output format among allowed ones
   Int_t            fSteps;         //  number of steps (only for container)
   TObjArray        fValues;        //  container for all related values
   Int_t            fNValues;       //! number of values (internal use)
   TList           *fList;          //! list containing the output
   Int_t            fIndex;         //  index of object in the list
   
   TArrayD          fArray;         //! temp array of computed values

   ClassDef(AliRsnListOutput, 1)    //  AliRsnListOutput class
};

#endif
