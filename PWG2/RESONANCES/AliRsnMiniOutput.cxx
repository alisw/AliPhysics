//
// Mini-Output
// All the definitions needed for building a RSN histogram
// including:
// -- properties of resonance (mass, PDG code if needed)
// -- properties of daughters (assigned mass, charges)
// -- definition of output histogram
// 

#include "Riostream.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TList.h"
#include "TMath.h"
#include "THnSparse.h"
#include "TString.h"
#include "TClonesArray.h"

#include "AliRsnMiniParticle.h"
#include "AliRsnMiniPair.h"
#include "AliRsnMiniEvent.h"

#include "AliLog.h"
#include "AliRsnCutSet.h"
#include "AliRsnMiniAxis.h"
#include "AliRsnMiniOutput.h"
#include "AliRsnMiniValue.h"

ClassImp(AliRsnMiniOutput)

//__________________________________________________________________________________________________
AliRsnMiniOutput::AliRsnMiniOutput() :
   TNamed(),
   fOutputType(kTypes),
   fComputation(kComputations),
   fMotherPDG(0),
   fMotherMass(0.0),
   fPairCuts(0x0),
   fOutputID(-1),
   fAxes("AliRsnMiniAxis", 0),
   fComputed(0),
   fPair(),
   fList(0x0)
{
//
// Constructor
//

   fCutID[0] = fCutID[1] = -1;
   fDaughter[0] = fDaughter[1] = AliRsnDaughter::kUnknown;
   fCharge[0] = fCharge[1] = 0;
}

//__________________________________________________________________________________________________
AliRsnMiniOutput::AliRsnMiniOutput(const char *name, EOutputType type, EComputation src) :
   TNamed(name, ""),
   fOutputType(type),
   fComputation(src),
   fMotherPDG(0),
   fMotherMass(0.0),
   fPairCuts(0x0),
   fOutputID(-1),
   fAxes("AliRsnMiniAxis", 0),
   fComputed(0),
   fPair(),
   fList(0x0)
{
//
// Constructor
//

   fCutID[0] = fCutID[1] = -1;
   fDaughter[0] = fDaughter[1] = AliRsnDaughter::kUnknown;
   fCharge[0] = fCharge[1] = 0;
}

//__________________________________________________________________________________________________
AliRsnMiniOutput::AliRsnMiniOutput(const char *name, const char *outType, const char *compType) :
   TNamed(name, ""),
   fOutputType(kTypes),
   fComputation(kComputations),
   fMotherPDG(0),
   fMotherMass(0.0),
   fPairCuts(0x0),
   fOutputID(-1),
   fAxes("AliRsnMiniAxis", 0),
   fComputed(0),
   fPair(),
   fList(0x0)
{
//
// Constructor, with a more user friendly implementation, where
// the user sets the type of output and computations through conventional strings:
//
// Output:
//    -- "HIST"    --> common histogram (up to 3 dimensions)
//    -- "SPARSE"  --> sparse histogram
//                 
// Computation:    
//    -- "EVENT"   --> event-only computations
//    -- "PAIR"    --> track pair computations (default)
//    -- "MIX"     --> event mixing (like track pair, but different events)
//    -- "ROTATE1" --> rotated background (rotate first track)
//    -- "ROTATE2" --> rotated background (rotate second track)
//    -- "TRUE"    --> true pairs (like track pair, but checking that come from same mother)
//    -- "MOTHER"  --> mother (loop on MC directly for mothers --> denominator of efficiency)
//

   TString input;
   
   // understand output type
   input = outType;
   input.ToUpper();
   if (!input.CompareTo("HIST"))
      fOutputType = kHistogram;
   else if (!input.CompareTo("SPARSE"))
      fOutputType = kHistogramSparse;
   else
      AliWarning(Form("String '%s' does not define a meaningful output type", outType));
      
   // understand computation type
   input = compType;
   input.ToUpper();
   if (!input.CompareTo("EVENT"))
      fComputation = kEventOnly;
   else if (!input.CompareTo("PAIR"))
      fComputation = kTrackPair;
   else if (!input.CompareTo("MIX"))
      fComputation = kTrackPairMix;
   else if (!input.CompareTo("ROTATE1"))
      fComputation = kTrackPairRotated1;
   else if (!input.CompareTo("ROTATE2"))
      fComputation = kTrackPairRotated2;
   else if (!input.CompareTo("TRUE"))
      fComputation = kTruePair;
   else if (!input.CompareTo("MOTHER"))
      fComputation = kMother;
   else
      AliWarning(Form("String '%s' does not define a meaningful computation type", compType));
   
   fCutID[0] = fCutID[1] = -1;
   fDaughter[0] = fDaughter[1] = AliRsnDaughter::kUnknown;
   fCharge[0] = fCharge[1] = 0;
}

//__________________________________________________________________________________________________
AliRsnMiniOutput::AliRsnMiniOutput(const AliRsnMiniOutput &copy) :
   TNamed(copy),
   fOutputType(copy.fOutputType),
   fComputation(copy.fComputation),
   fMotherPDG(copy.fMotherPDG),
   fMotherMass(copy.fMotherMass),
   fPairCuts(copy.fPairCuts),
   fOutputID(copy.fOutputID),
   fAxes(copy.fAxes),
   fComputed(copy.fComputed),
   fPair(),
   fList(copy.fList)
{
//
// Copy constructor
//

   Int_t i;
   for (i = 0; i < 2; i++) {
      fCutID[i] = copy.fCutID[i];
      fDaughter[i] = copy.fDaughter[i];
      fCharge[i] = copy.fCharge[i];
   }
}

//__________________________________________________________________________________________________
AliRsnMiniOutput& AliRsnMiniOutput::operator=(const AliRsnMiniOutput &copy)
{
//
// Assignment operator
//

   fOutputType = copy.fOutputType;
   fComputation = copy.fComputation;
   fMotherPDG = copy.fMotherPDG;
   fMotherMass = copy.fMotherMass;
   fPairCuts = copy.fPairCuts;
   fOutputID = copy.fOutputID;
   fAxes = copy.fAxes;
   fComputed = copy.fComputed;
   fList = copy.fList;

   Int_t i;
   for (i = 0; i < 2; i++) {
      fCutID[i] = copy.fCutID[i];
      fDaughter[i] = copy.fDaughter[i];
      fCharge[i] = copy.fCharge[i];
   }
   
   return (*this);
}


//__________________________________________________________________________________________________
void AliRsnMiniOutput::AddAxis(Int_t i, Int_t nbins, Double_t min, Double_t max)
{
//
// Create a new axis reference
//

   Int_t size = fAxes.GetEntries();   
   new (fAxes[size]) AliRsnMiniAxis(i, nbins, min, max);
}

//__________________________________________________________________________________________________
void AliRsnMiniOutput::AddAxis(Int_t i, Double_t min, Double_t max, Double_t step)
{
//
// Create a new axis reference
//

   Int_t size = fAxes.GetEntries();   
   new (fAxes[size]) AliRsnMiniAxis(i, min, max, step);
}

//__________________________________________________________________________________________________
void AliRsnMiniOutput::AddAxis(Int_t i, Int_t nbins, Double_t *values)
{
//
// Create a new axis reference
//

   Int_t size = fAxes.GetEntries();   
   new (fAxes[size]) AliRsnMiniAxis(i, nbins, values);
}

//__________________________________________________________________________________________________
Bool_t AliRsnMiniOutput::Init(const char *prefix, TList *list)
{
//
// Initialize properly the histogram and add it to the argument list
//

   if (!list) {
      AliError("Required an output list");
      return kFALSE;
   }
   
   fList = list;
   Int_t size = fAxes.GetEntries();
   if (size < 1) {
      AliWarning(Form("[%s] Cannot initialize histogram with less than 1 axis", GetName()));
      return kFALSE;
   }

   switch (fOutputType) {
      case kHistogram:
         CreateHistogram(Form("%s_%s", prefix, GetName()));
         return kTRUE;
      case kHistogramSparse:
         CreateHistogramSparse(Form("%s_%s", prefix, GetName()));
         return kTRUE;
      default:
         AliError("Wrong output histogram definition");
         return kFALSE;
   }
}

//__________________________________________________________________________________________________
void AliRsnMiniOutput::CreateHistogram(const char *name)
{
//
// Initialize the 'default' TH1 output object.
// In case one of the expected axes is NULL, the initialization fails.
//

   Int_t size = fAxes.GetEntries();
   AliInfo(Form("Histogram name = '%s', with %d axes", name, size));

   // we expect to have maximum 3 axes in this case
   AliRsnMiniAxis *xAxis = 0x0, *yAxis = 0x0, *zAxis = 0x0;
   if (size >= 1) xAxis = (AliRsnMiniAxis*)fAxes[0];
   if (size >= 2) yAxis = (AliRsnMiniAxis*)fAxes[1];
   if (size >= 3) zAxis = (AliRsnMiniAxis*)fAxes[2];
   
   // create histogram depending on the number of axes
   TH1 *h1 = 0x0;
   if (xAxis && yAxis && zAxis) {
      h1 = new TH3F(name, "", xAxis->NBins(), xAxis->BinArray(), yAxis->NBins(), yAxis->BinArray(), zAxis->NBins(), zAxis->BinArray());
   } else if (xAxis && yAxis) {
      h1 = new TH2F(name, "", xAxis->NBins(), xAxis->BinArray(), yAxis->NBins(), yAxis->BinArray());
   } else if (xAxis) {
      h1 = new TH1F(name, "", xAxis->NBins(), xAxis->BinArray());
   } else {
      AliError("No axis was initialized");
      return;
   }
   
   // switch the correct computation of errors
   if (h1 && fList) {
      h1->Sumw2();
      fList->Add(h1);
      fOutputID = fList->IndexOf(h1);
   }
}

//________________________________________________________________________________________
void AliRsnMiniOutput::CreateHistogramSparse(const char *name)
{
//
// Initialize the THnSparse output object.
// In case one of the expected axes is NULL, the initialization fails.
//

   // retrieve binnings and sizes of all axes
   // since the check for null values is done in Init(),
   // we assume that here they must all be well defined
   Int_t size = fAxes.GetEntries();
   Int_t i, *nbins = new Int_t[size];
   for (i = 0; i < size; i++) {
      AliRsnMiniAxis *axis = (AliRsnMiniAxis*)fAxes[i];
      nbins[i] = axis->NBins();
   }

   // create fHSparseogram
   THnSparseF *h1 = new THnSparseF(name, "", size, nbins);

   // update the various axes using the definitions given in the array of axes here
   for (i = 0; i < size; i++) {
      AliRsnMiniAxis *axis = (AliRsnMiniAxis*)fAxes[i];
      h1->GetAxis(i)->Set(nbins[i], axis->BinArray());
   }

   // clear heap
   delete [] nbins;
   
   // add to list
   if (h1 && fList) {
      h1->Sumw2();
      fList->Add(h1);
      fOutputID = fList->IndexOf(h1);
   }
}

//________________________________________________________________________________________
Bool_t AliRsnMiniOutput::FillEvent(AliRsnMiniEvent *event, TClonesArray *valueList)
{
//
// Compute values for event-based computations (does not use the pair)
//

   // check computation type
   if (fComputation != kEventOnly) {
      AliError("This method can be called only for event-based computations");
      return kFALSE;
   }

   // compute & fill
   ComputeValues(event, valueList);
   FillHistogram();
   return kTRUE;
}

//________________________________________________________________________________________
Bool_t AliRsnMiniOutput::FillMother(const AliRsnMiniPair *pair, AliRsnMiniEvent *event, TClonesArray *valueList)
{
//
// Compute values for mother-based computations
//

   // check computation type
   if (fComputation != kMother) {
      AliError("This method can be called only for mother-based computations");
      return kFALSE;
   }
   
   // copy passed pair info
   fPair = (*pair);
   
   // check pair against cuts
   if (fPairCuts) if (!fPairCuts->IsSelected(&fPair)) return kFALSE;

   // compute & fill
   ComputeValues(event, valueList);
   FillHistogram();
   return kTRUE;
}

//________________________________________________________________________________________
Int_t AliRsnMiniOutput::FillPair(AliRsnMiniEvent *event1, AliRsnMiniEvent *event2, TClonesArray *valueList)
{
//
// Loops on the passed mini-event, and for each pair of particles
// which satisfy the charge and cut requirements defined here, add an entry.
// Returns the number of successful fillings.
//

   // check computation type
   Bool_t okComp = kFALSE;
   if (fComputation == kTrackPair)         okComp = kTRUE;
   if (fComputation == kTrackPairMix)      okComp = kTRUE;
   if (fComputation == kTrackPairRotated1) okComp = kTRUE;
   if (fComputation == kTrackPairRotated2) okComp = kTRUE;
   if (fComputation == kTruePair)          okComp = kTRUE;
   if (!okComp) {
      AliError(Form("[%s] This method can be called only for pair-based computations", GetName()));
      return kFALSE;
   }
   
   // loop variables
   Int_t i1, i2, start, nadded = 0;
   Int_t n1 = event1->Particles().GetEntriesFast();
   Int_t n2 = event2->Particles().GetEntriesFast();
   AliRsnMiniParticle *p1, *p2;
   
   // it is necessary to know if criteria for the two daughters are the same
   // and if the two events are the same or not (mixing)
   Bool_t sameCriteria = ((fCharge[0] == fCharge[1]) && (fCutID[0] == fCutID[1]));
   Bool_t sameEvent = (event1->ID() == event2->ID());
   
   //Int_t nsel1 = event1->CountParticles(fCharge[0], fCutID[0]);
   //Int_t nsel2 = event2->CountParticles(fCharge[1], fCutID[1]);
   //cout << "Charge #1 = " << fCharge[0] << " cut bit #1 = " << fCutID[0] << ", tracks = " << nsel1 << endl;
   //cout << "Charge #2 = " << fCharge[1] << " cut bit #2 = " << fCutID[1] << ", tracks = " << nsel2 << endl;
   //if (!nsel1 || !nsel2) {
   //   cout << "Nothing to mix" << endl;
   //   return 0;
   //}
   
   // external loop
   for (i1 = 0; i1 < n1; i1++) {
      p1 = event1->GetParticle(i1);
      if (p1->Charge() != fCharge[0]) continue;
      if (!p1->HasCutBit(fCutID[0])) continue;
      // define starting point for inner loop
      // if daughter selection criteria (charge, cuts) are the same
      // and the two events coincide, internal loop must start from
      // the first track *after* current one;
      // otherwise it starts from the beginning
      start = ((sameEvent && sameCriteria) ? i1 + 1 : 0);
      // internal loop
      for (i2 = start; i2 < n2; i2++) {
         p2 = event2->GetParticle(i2);
         if (p2->Charge() != fCharge[1]) continue;
         if (!p2->HasCutBit(fCutID[1])) continue;
         //cout << "Matching: " << i1 << ' ' << i2 << endl;
         // avoid to mix a particle with itself
         if (sameEvent && (p1->Index() == p2->Index())) {
            //cout << "-- skipping same index" << endl;
            continue;
         }
         // sum momenta
         fPair.Fill(p1, p2, GetMass(0), GetMass(1), fMotherMass);
         // do rotation if needed
         if (fComputation == kTrackPairRotated1) fPair.InvertP(kTRUE);
         if (fComputation == kTrackPairRotated2) fPair.InvertP(kFALSE);
         // if required, check that this is a true pair
         if (fComputation == kTruePair) {
            if (fPair.Mother() < 0)  {
               continue;
            } else if (TMath::Abs(fPair.MotherPDG()) != fMotherPDG) {
               continue;
            }
            Bool_t decayMatch = kFALSE;
            if (p1->PDGAbs() == AliRsnDaughter::SpeciesPDG(fDaughter[0]) && p2->PDGAbs() == AliRsnDaughter::SpeciesPDG(fDaughter[1]))
               decayMatch = kTRUE;
            if (p2->PDGAbs() == AliRsnDaughter::SpeciesPDG(fDaughter[0]) && p1->PDGAbs() == AliRsnDaughter::SpeciesPDG(fDaughter[1]))
               decayMatch = kTRUE;
            if (!decayMatch) continue;
         }
         // check pair against cuts
         if (fPairCuts) {
            if (!fPairCuts->IsSelected(&fPair)) {
               //cout << "-- cuts not passed: Y = " << TMath::Abs(fPair.Y(0)) << endl;
               continue;
            }
         }
         // get computed values & fill histogram
         ComputeValues(event1, valueList);
         FillHistogram();
         nadded++;
      } // end internal loop
   } // end external loop
   
   return nadded;
}

//________________________________________________________________________________________
void AliRsnMiniOutput::ComputeValues(AliRsnMiniEvent *event, TClonesArray *valueList)
{
//
// Using the arguments and the internal 'fPair' data member,
// compute all values to be stored in the histogram
//

   // check size of computed array
   Int_t size = fAxes.GetEntries();
   if (fComputed.GetSize() != size) fComputed.Set(size);

   Int_t i, ival, nval = valueList->GetEntries();
   
   for (i = 0; i < size; i++) {
      fComputed[i] = 1E20;
      AliRsnMiniAxis *axis = (AliRsnMiniAxis*)fAxes[i];
      if (!axis) {
         AliError("Null axis");
         continue;
      }
      ival = axis->GetValueID();
      if (ival < 0 || ival >= nval) {
         AliError(Form("Required value #%d, while maximum is %d", ival, nval));
         continue;
      }
      AliRsnMiniValue *val = (AliRsnMiniValue*)valueList->At(ival);
      if (!val) {
         AliError(Form("Value in position #%d is NULL", ival));
         continue;
      } 
      // if none of the above exit points is taken, compute value
      fComputed[i] = val->Eval(&fPair, event);
   }
}

//________________________________________________________________________________________
void AliRsnMiniOutput::FillHistogram()
{
//
// Fills the internal histogram using the current values stored in the
// 'fComputed' array, in the order as they are stored, up to the max
// dimension of the initialized histogram itself.
//

   // retrieve object from list
   if (!fList) {
      AliError("List pointer is NULL");
      return;
   }
   TObject *obj = fList->At(fOutputID);

   if (obj->InheritsFrom(TH1F::Class())) {
      ((TH1F*)obj)->Fill(fComputed[0]);
   } else if (obj->InheritsFrom(TH2F::Class())) {
      ((TH2F*)obj)->Fill(fComputed[0], fComputed[1]);
   } else if (obj->InheritsFrom(TH3F::Class())) {
      ((TH3F*)obj)->Fill(fComputed[0], fComputed[1], fComputed[2]);
   } else if (obj->InheritsFrom(THnSparseF::Class())) {
      ((THnSparseF*)obj)->Fill(fComputed.GetArray());
   } else {
      AliError("No output initialized");
   }
}
