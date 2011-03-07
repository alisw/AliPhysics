//
// Class AliRsn Fcn
//
// This class defines a base classe to implement a typical computation
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
// this class and override the virtual functions defined in it, which
// initialize the final output histogram and define how to process data.
//
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNFUNCTION_H
#define ALIRSNFUNCTION_H

#include <TClonesArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TArrayD.h>

class AliRsnValue;

class AliRsnFunction : public TObject {

public:

   AliRsnFunction(Bool_t useTH1 = kTRUE);
   AliRsnFunction(const AliRsnFunction &copy);
   virtual ~AliRsnFunction() { delete fH1; delete fHSparse; }
   const AliRsnFunction& operator=(const AliRsnFunction &copy);

   virtual const char*  GetName() const;

   Bool_t               IsUsingTH1()  {return fUseTH1;}
   void                 UseTH1()      {fUseTH1 = kTRUE;}
   void                 UseSparse()   {fUseTH1 = kFALSE;}
   
   Bool_t               AddAxis(AliRsnValue* const axis);
   Int_t                GetNumberOfAxes() {return fAxisList.GetEntries();}
   
   TH1*                 CreateHistogram(const char *histoName, const char *histoTitle);
   THnSparseF*          CreateHistogramSparse(const char *histoName, const char *histoTitle);

   Bool_t               Fill(TObject *object);

protected:

   TClonesArray   fAxisList;    //  list of axis
   Bool_t         fUseTH1;      //  use TH1 or not?
   Int_t          fSize;        //  number of dim of output histogram
   TH1           *fH1;          //  output histogram (standard type)
   THnSparseF    *fHSparse;     //  output histogram (sparse type)
   TArrayD        fValues;      //! computed values

   // ROOT dictionary
   ClassDef(AliRsnFunction, 3)
};

#endif
