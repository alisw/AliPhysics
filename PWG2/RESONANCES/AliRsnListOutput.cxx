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

#include <Riostream.h>

#include "AliLog.h"

#include "AliRsnValue.h"
#include "AliRsnLoop.h"

#include "AliRsnListOutput.h"

ClassImp(AliRsnListOutput)

//________________________________________________________________________________________
AliRsnListOutput::AliRsnListOutput(const char *name, AliRsnListOutput::EOut type) :
   TNamed(name, ""),
   fType(type),
   fSteps(0),
   fValues(0),
   fNValues(0),
   fList(0x0),
   fIndex(-1),
   fArray(0)
{
//
// Constructor.
// Requires a name for this object (which will be used to name the output object)
// and the definition of the output type from the built-in enumeration.
//
}

//________________________________________________________________________________________
AliRsnListOutput::AliRsnListOutput(const AliRsnListOutput &copy) :
   TNamed(copy),
   fType(copy.fType),
   fSteps(copy.fSteps),
   fValues(copy.fValues),
   fNValues(copy.fNValues),
   fList(copy.fList),
   fIndex(copy.fIndex),
   fArray(0)
{
//
// Copy constructor.
// Since the pointer objects must be initialized in a second step,
// they are never copied, and then they are initialized to zero.
//
}

//________________________________________________________________________________________
const AliRsnListOutput& AliRsnListOutput::operator=(const AliRsnListOutput& copy)
{
//
// Assignment operator.
// Same consideration as the copiy constructor. In this case, there is
// the possibility to have the output objects alreasy initialized, but
// they are anyway reset.
//

   TNamed::operator=(copy);

   fType = copy.fType;
   fSteps = copy.fSteps;
   fValues = copy.fValues;
   fNValues = copy.fNValues;
   fList = copy.fList;
   fIndex = copy.fIndex;
   fArray = copy.fArray;

   Reset();

   return (*this);
}

//__________________________________________________________________________________________________
AliRsnListOutput::~AliRsnListOutput()
{
//
// Destructor.
// Deletes the output objects.
//

   Reset();
}

//__________________________________________________________________________________________________
void AliRsnListOutput::Reset()
{
//
// Clear all output objects. In general, only one will be initialized at a time.
//

   fList = 0x0;
}

//_____________________________________________________________________________
void AliRsnListOutput::AddValue(AliRsnValue *value)
{
//
// Adds a value computation object to the list.
//

   fValues.AddLast(value);
}


//________________________________________________________________________________________
Bool_t AliRsnListOutput::Init(const char *prefix, TList *list)
{
//
// Initializes the output for this object.
// What kind of output depends on the 'fType' data member,
// and in case it is a CF container, also on the 'fSteps'.
// The object is named with the following criterion:
// <prefix>_<name>_<type>_<varList>
//

   Int_t i;

   // all output objects are cleared
   Reset();

   // count values and set dimension of arrays
   // do also some checks for a good match between type and output
   fNValues = fValues.GetEntries();
   if (fNValues < 1) {
      AliError("Need at least 1 value");
      return kFALSE;
   }
   if (fType == kHistoDefault && fNValues > 3) {
      AliInfo(Form("NValues = %d > 3 --> cannot use a normal histogram, need to use a sparse", fNValues));
      fType = kHistoSparse;
   }
   
   // resize the output array
   fArray.Set(fNValues);

   // create the name
   TString name(Form("%s_%s", prefix, GetName()));
   AliRsnValue *val = 0x0;
   for (i = 0; i < fNValues; i++) {
      val = GetValue(i);
      if (!val) {
         AliError(Form("Slot %d in value list is NULL", i));
         return kFALSE;
      }
      name += '_';
      name += val->GetName();
   }
   
   // allowed objects
   TObject *object = 0x0;

   // initialize appropriate output object
   // and, if successful, insert into the list
   switch (fType) {
      case kHistoDefault:
         name.Append("_hist");
         object = CreateHistogram(name.Data());
         break;
      case kHistoSparse:
         name.Append("_hsparse");
         object = CreateHistogramSparse(name.Data());
         break;
      case kCFContainer:
         name.Append("_cf");
         object = CreateCFContainer(name.Data());
         break;
      default:
         AliWarning("Wrong type output or initialization failure");
   }
   
   if (object) {
      //AliInfo(Form("[%s]: initializing output '%s' (obj name = '%s') with %d values and format %d [%s]", GetName(), name.Data(), object->GetName(), fNValues, fType, object->ClassName()));
      fList = list;
      fList->Add(object);
      fIndex = fList->IndexOf(object);
      return kTRUE;
   }

   return kFALSE;
}

//________________________________________________________________________________________
TH1* AliRsnListOutput::CreateHistogram(const char *name)
{
//
// Initialize the 'default' TH1 output object.
// In case one of the expected axes is NULL, the initialization fails.
//

   // we expect to have maximum 3 axes in this case
   Int_t i, nbins[3] = {0, 0, 0};
   TArrayD  array[3];
   for (i = 0; i < fNValues; i++) {
      AliRsnValue *val = GetValue(i);
      if (!val) {
         AliError(Form("Expected axis %d is NULL", i));
         return 0x0;
      }
      nbins[i] = GetValue(i)->GetArray().GetSize() - 1;
      array[i] = GetValue(i)->GetArray();
   }
   
   TH1 *hist = 0x0;

   // create histogram depending on the number of axes
   switch (fNValues) {
      case 1:
         hist = new TH1F(name, "", nbins[0], array[0].GetArray());
         break;
      case 2:
         hist = new TH2F(name, "", nbins[0], array[0].GetArray(), nbins[1], array[1].GetArray());
         break;
      case 3:
         hist = new TH3F(name, "", nbins[0], array[0].GetArray(), nbins[1], array[1].GetArray(), nbins[2], array[2].GetArray());
         break;
      default:
         AliError(Form("Wrong number of dimensions: %d", fNValues))
         return 0x0;
   }

   if (hist) hist->Sumw2();
   return hist;
}

//________________________________________________________________________________________
THnSparseF* AliRsnListOutput::CreateHistogramSparse(const char *name)
{
//
// Initialize the THnSparse output object.
// In case one of the expected axes is NULL, the initialization fails.
//

   // retrieve binnings and sizes of all axes
   // since the check for null values is done in Init(),
   // we assume that here they must all be well defined
   Int_t i, *nbins = new Int_t[fNValues];
   TArrayD  *array = new TArrayD[fNValues];
   for (i = 0; i < fNValues; i++) {
      nbins[i] = GetValue(i)->GetArray().GetSize() - 1;
      array[i] = GetValue(i)->GetArray();
   }

   // create histogram
   THnSparseF *hist = new THnSparseF(name, "", fNValues, nbins);
   hist->Sumw2();

   // update the various axes using the definitions given in the array of axes here
   for (i = 0; i < fNValues; i++) {
      hist->GetAxis(i)->Set(nbins[i], array[i].GetArray());
   }

   // clear heap
   delete [] nbins;
   delete [] array;

   return hist;
}

//________________________________________________________________________________________
AliCFContainer* AliRsnListOutput::CreateCFContainer(const char *name)
{
//
// Initialize the AliCFContainer output object.
// In case one of the expected axes is NULL, the initialization fails.
//

   // retrieve binnings and sizes of all axes
   // since the check for null values is done in Init(),
   // we assume that here they must all be well defined
   Int_t i, *nbins = new Int_t[fNValues];
   TArrayD  *array = new TArrayD[fNValues];
   for (i = 0; i < fNValues; i++) {
      nbins[i] = GetValue(i)->GetArray().GetSize() - 1;
      array[i] = GetValue(i)->GetArray();
   }

   // create object
   AliCFContainer *cont = new AliCFContainer(name, "", fSteps, fNValues, nbins);

   // set the bin limits for each axis
   for (i = 0; i < fNValues; i++) {
      cont->SetBinLimits(i, array[i].GetArray());
   }

   // clear heap
   delete [] nbins;
   delete [] array;

   return cont;
}

//________________________________________________________________________________________
Bool_t AliRsnListOutput::Fill(TObject *target, Int_t step)
{
//
// Uses the passed argument to compute all values.
// If all computations were successful, fill the output
// Second argument (step) is needed only in case of CF containers.
// Return value is the AND of all computation successes.
//

   Int_t  i;

   // do computations
   Bool_t globalOK = kTRUE;
   AliRsnValue *val = 0x0;
   for (i = 0; i < fNValues; i++) {
      val = GetValue(i);
      if (!val) {
         AliError("NULL value found");
         return kFALSE;
      }
      globalOK = globalOK && val->Eval(target);
      fArray[i] = (Double_t)val->GetComputedValue();
   }
   if (!globalOK) return kFALSE;
   
   // retrieve object
   if (!fList || fIndex < 0) {
      AliError("List not initialized");
      return kFALSE;
   }
   TObject *obj = fList->At(fIndex);
   if (!obj) {
      AliError("Null pointer");
      return kFALSE;
   }
   
   // check
   //AliInfo(Form("[%s] Object index, name, type = %d, %s (%s)", GetName(), fIndex, obj->GetName(), obj->ClassName()));

   // fill output
   if (obj->IsA() == TH1F::Class()) {
      TH1F *h = (TH1F*)obj;
      h->Fill(fArray[0]);
      return kTRUE;
   } else if (obj->IsA() == TH2F::Class()) {
      TH2F *h = (TH2F*)obj;
      h->Fill(fArray[0], fArray[1]);
      return kTRUE;
   } else if (obj->IsA() == TH3F::Class()) {
      TH3F *h = (TH3F*)obj;
      h->Fill(fArray[0], fArray[1], fArray[2]);
      return kTRUE;
   } else if (obj->InheritsFrom(THnSparse::Class())) {
      THnSparseF *h = (THnSparseF*)obj;
      h->Fill(fArray.GetArray());
      return kTRUE;
   } else if (obj->InheritsFrom(AliCFContainer::Class())) {
      AliCFContainer *c = (AliCFContainer*)obj;
      c->Fill(fArray.GetArray(), step);
      return kTRUE;
   } else {
      AliError(Form("Not handled class '%s'", obj->ClassName()));
      return kFALSE;
   }
}
