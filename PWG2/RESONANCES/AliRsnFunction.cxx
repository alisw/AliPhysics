//
// Class AliRsnFunction
//
// This class defines a base classe to implement a function
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

#include <TString.h>
#include <TAxis.h>

#include "AliLog.h"

#include "AliRsnValue.h"

#include "AliRsnFunction.h"

ClassImp(AliRsnFunction)

//________________________________________________________________________________________
AliRsnFunction::AliRsnFunction(Bool_t useTH1) :
   TObject(),
   fAxisList("AliRsnValue", 0),
   fUseTH1(useTH1),
   fSize(0),
   fH1(0x0),
   fHSparse(0x0),
   fValues(0)
{
//
// Constructor.
//
}

//________________________________________________________________________________________
AliRsnFunction::AliRsnFunction(const AliRsnFunction &copy) :
   TObject(copy),
   fAxisList(copy.fAxisList),
   fUseTH1(copy.fUseTH1),
   fSize(copy.fSize),
   fH1(0x0),
   fHSparse(0x0),
   fValues(copy.fValues)
{
//
// Copy constructor.
//
}

//________________________________________________________________________________________
const AliRsnFunction& AliRsnFunction::operator=(const AliRsnFunction& copy)
{
//
// Assignment operator.
//

   fAxisList = copy.fAxisList;
   fUseTH1 = copy.fUseTH1;
   fSize = copy.fSize;
   fValues = copy.fValues;

   if (fH1) delete fH1;
   if (fHSparse) delete fHSparse;

   return (*this);
}

//________________________________________________________________________________________
const char* AliRsnFunction::GetName() const
{
//
// Defines the name of this object according to
// the function type and binning
//

   TString name("");

   TObjArrayIter next(&fAxisList);
   AliRsnValue *axis = 0;

   while ((axis = (AliRsnValue*)next())) {
      if (name.Length() > 1) name += '_';
      name += axis->GetName();
   }

   return name.Data();
}

//________________________________________________________________________________________
Bool_t AliRsnFunction::AddAxis(AliRsnValue *const axis)
{
//
// Try to add a new axis to this function.
// The operation succeeds only if the related value object
// has a target compatible with the function type:
// -- 'single'     functions, for tracks/events: AliRsnDaughter or AliRsnEvent
// -- 'not single' functions, for pairs/events : AliRsnMother or AliRsnEvent
// otherwise the axis is not added.
//
// If more than 3 axes are added, switch to THnSparseF output.
// NOTE: this can cause large files and high memory occupancy.
//

   Int_t size = fAxisList.GetEntries();
   new (fAxisList[size]) AliRsnValue(*axis);

   if (fAxisList.GetEntries() > 3) {
      AliWarning("Adding more than 3 axes will produce a THnSparseD output.");
      fUseTH1 = kFALSE;
   }

   return kTRUE;
}

//________________________________________________________________________________________
TH1* AliRsnFunction::CreateHistogram(const char *histoName, const char *histoTitle)
{
//
// Creates and returns the histogram defined using
// arguments fo name and title, and the first histoDef for binning.
// Variable-sized histogram binning is always called, due to use of histoDef,
// even if the bins are equal, since they are defined in this class.
// Eventually present histoDef's in other slots of array (1, 2) are ignored.
//
// This version produces a THnSparseF.
//

   fSize = fAxisList.GetEntries();
   if (!fSize) {
      AliError("No axes defined");
      return 0x0;
   } else if (fSize > 3) {
      AliError("Too many axes defined for a TH1 object");
      return 0x0;
   }
   
   // initialize the size of values container
   fValues.Set(fSize);

   // retrieve binnings for main and secondary axes
   AliRsnValue *fcnAxis;
   TArrayD      array[3];
   for (Int_t i = 0; i < fSize; i++) {
      fcnAxis = (AliRsnValue*)fAxisList.At(i);
      if (!fcnAxis) {
         AliError("Empty axis");
         array[i].Set(2);
         array[i][0] = -1E5;
         array[i][1] = -1E5;
         continue;
      } else {
         array[i] = fcnAxis->GetArray();
      }
   }

   // create histogram depending on the number of axes
   switch (fSize) {
      case 1:
         fH1 = new TH1F(histoName, histoTitle, array[0].GetSize() - 1, array[0].GetArray());
         break;
      case 2:
         fH1 = new TH2F(histoName, histoTitle, array[0].GetSize() - 1, array[0].GetArray(), array[1].GetSize() - 1, array[1].GetArray());
         break;
      case 3:
         fH1 = new TH3F(histoName, histoTitle, array[0].GetSize() - 1, array[0].GetArray(), array[1].GetSize() - 1, array[1].GetArray(), array[2].GetSize() - 1, array[2].GetArray());
         break;
   }
   fH1->Sumw2();

   return fH1;
}

//________________________________________________________________________________________
THnSparseF* AliRsnFunction::CreateHistogramSparse(const char *histoName, const char *histoTitle)
{
//
// Creates and returns the histogram defined using
// arguments fo name and title, and the first histoDef for binning.
// Variable-sized histogram binning is always called, due to use of histoDef,
// even if the bins are equal, since they are defined in this class.
// Eventually present histoDef's in other slots of array (1, 2) are ignored.
//
// This version produces a THnSparseF.
//

   fSize = fAxisList.GetEntries();
   if (!fSize) {
      AliError("No axes defined");
      return 0x0;
   }
   
   // initialize the size of values container
   fValues.Set(fSize);

   // initialize the array of number of bins for each axis
   // taking it from the stored values, while for the bins
   // they are set as summied and defined later
   Int_t       *nbins   = new Int_t[fSize];
   AliRsnValue *fcnAxis = 0;
   for (Int_t i = 0; i < fSize; i++) {
      fcnAxis = (AliRsnValue*)fAxisList.At(i);
      if (!fcnAxis) {
         nbins[i] = 1;
         AliError("Empty axis");
         continue;
      }
      nbins[i] = fcnAxis->GetArray().GetSize() - 1;
   }

   // create histogram
   fHSparse = new THnSparseF(histoName, histoTitle, fSize, nbins);
   fHSparse->Sumw2();

   // update the various axes using the definitions given in the array of axes here
   for (Int_t i = 0; i < fSize; i++) {
      fcnAxis = (AliRsnValue*)fAxisList.At(i);
      if (!fcnAxis) {
         AliError("Empty axis: doing unique bin betweeen -100000 and 100000");
         continue;
      }
      TAxis* axis = fHSparse->GetAxis(i);
      axis->Set(nbins[i], fcnAxis->GetArray().GetArray());
   }

   delete [] nbins;

   return fHSparse;
}


//________________________________________________________________________________________
Bool_t AliRsnFunction::Fill(TObject *object)
{
//
// Fill function histogram using the passed object
//

   AliDebug(AliLog::kDebug + 2, "->");

   // loop on axes and try to compute values
   // using this object or, as an alternative
   // its reference event
   Int_t  i;
   Bool_t globalOK = kTRUE, computeOK;
   AliRsnValue *fcnAxis = 0;
   for (i = 0; i < fSize; i++) {
      fValues[i] = 0.0;
      computeOK = kFALSE;
      fcnAxis = (AliRsnValue*)fAxisList.At(i);
      if (fcnAxis) {
         computeOK = fcnAxis->Eval(object);
         if (computeOK) fValues[i] = ((Float_t)fcnAxis->GetComputedValue());
      }
      if (!computeOK) globalOK = kFALSE;
   }
   
   // if even one of the computations has failes, the histograms are not filled
   if (!globalOK) return kFALSE;

   // fill histogram
   if (fUseTH1) {
      
      // check presence of output histogram
      if (!fH1) {
         AliError("Required a TH1 which is not initialized");
         return kFALSE;
      }

      // fill according to dimensions
      switch (fSize) {
         case 1: {
            TH1F *h1 = (TH1F*)fH1;
            h1->Fill(fValues[0]);
         }
         break;
         case 2: {
            TH2F *h2 = (TH2F*)fH1;
            h2->Fill(fValues[0], fValues[1]);
         }
         break;
         case 3: {
            TH3F *h3 = (TH3F*)fH1;
            h3->Fill(fValues[0], fValues[1], fValues[2]);
         }
         break;
         default:
            AliError(Form("Wrong size : %d", fSize));
            return kFALSE;
      }
   } else {
      
      // check presence of output histogram
      if (!fHSparse) {
         AliError("Required a THnSparseF which is not initialized");
         return kFALSE;
      }

      fHSparse->Fill(fValues.GetArray());
   }

   AliDebug(AliLog::kDebug + 2, "->");
   return kTRUE;
}
