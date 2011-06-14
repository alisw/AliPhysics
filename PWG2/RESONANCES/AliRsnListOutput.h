#ifndef ALIRSNLISTOUTPUT_H
#define ALIRSNLISTOUTPUT_H

//
// Class AliRsnListOutput
//
// This class defines a base classe to implement a Output
// which uses the internal RSN package event format (AliRsnEvent).
// It contains some default flags which turn out to be useful:
//  - a flag to select only the "true" pairs (tracks from same resonance)
//  - a flag to know if the computation is done over two events (mixing)
//
// Any kind of analysis object should be implemented as inheriting from this
// because the AliRsnAnalyzer which executes the analysis will accept a collection
// of such objects, in order to have a unique format of processing method
//
// The user who implements a kind of computation type should inherit from
// this class and override the virtual Outputs defined in it, which
// initialize the final output histogram and define how to process data.
//
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
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

class AliCFContainer;
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

   EOut            GetType() const         {return  fType;}
   Int_t           GetSteps() const        {return  fSteps;}
   TObjArray*      GetValues()             {return &fValues;}
   Int_t           GetNValues()            {return (fNValues = fValues.GetEntries());}
   AliRsnValue*    GetValue(Int_t i) const {return (AliRsnValue*)fValues[i];}
   Int_t           GetIndex() const        {return  fIndex;}
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
