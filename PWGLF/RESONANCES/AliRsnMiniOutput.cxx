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
#include "AliAODEvent.h"

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
   fList(0x0),
   fSel1(0),
   fSel2(0),
   fMaxNSisters(-1),
   fCheckP(kFALSE),
   fCheckDecay(kTRUE),
   fCheckFeedDown(kFALSE),   
   fOriginDselection(kFALSE),
   fKeepDfromB(kFALSE),
   fKeepDfromBOnly(kFALSE),
   fRejectIfNoQuark(kFALSE),
   fCheckHistRange(kTRUE)
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
   fList(0x0),
   fSel1(0),
   fSel2(0),
   fMaxNSisters(-1),
   fCheckP(kFALSE),
   fCheckDecay(kTRUE),
   fCheckFeedDown(kFALSE),   
   fOriginDselection(kFALSE),
   fKeepDfromB(kFALSE),
   fKeepDfromBOnly(kFALSE),
   fRejectIfNoQuark(kFALSE),
   fCheckHistRange(kTRUE)
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
   fList(0x0),
   fSel1(0),
   fSel2(0),
   fMaxNSisters(-1),
   fCheckP(kFALSE),
   fCheckDecay(kTRUE),
   fCheckFeedDown(kFALSE),   
   fOriginDselection(kFALSE),
   fKeepDfromB(kFALSE),
   fKeepDfromBOnly(kFALSE),
   fRejectIfNoQuark(kFALSE),
   fCheckHistRange(kTRUE)
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
//    -- "MOTHER_IN_ACC"  --> mother (loop on MC directly for mothers (in a defined acceptance interval)--> needed for efficiency calcutation using  an enriched sample)
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
    else if (!input.CompareTo("MOTHER_IN_ACC"))
      fComputation = kMotherInAcc;   
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
   fList(copy.fList),
   fSel1(0),
   fSel2(0),
   fMaxNSisters(-1),
   fCheckP(kFALSE),
   fCheckDecay(copy.fCheckDecay),
   fCheckFeedDown(kFALSE),   
   fOriginDselection(kFALSE),
   fKeepDfromB(kFALSE),
   fKeepDfromBOnly(kFALSE),
   fRejectIfNoQuark(kFALSE),
   fCheckHistRange(copy.fCheckHistRange)
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
AliRsnMiniOutput &AliRsnMiniOutput::operator=(const AliRsnMiniOutput &copy)
{
//
// Assignment operator
//
   if (this == &copy)
      return *this;
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

   fSel1.Set(0);
   fSel2.Set(0);
   fMaxNSisters = copy.fMaxNSisters;
   fCheckP = copy.fCheckP;
   fCheckDecay = copy.fCheckDecay;
   fCheckFeedDown = copy.fCheckFeedDown;
   fOriginDselection = copy.fOriginDselection;
   fKeepDfromB = copy.fOriginDselection;
   fKeepDfromBOnly = copy.fKeepDfromBOnly;
   fRejectIfNoQuark = copy.fRejectIfNoQuark;
   fCheckHistRange = copy.fCheckHistRange;

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
         if (size <= 3) {
            CreateHistogram(Form("%s_%s", prefix, GetName()));
         } else {
            AliInfo(Form("[%s] Added %d > 3 axes. Creating a sparse histogram", GetName(), size));
            fOutputType = kHistogramSparse;
            CreateHistogramSparse(Form("%s_%s", prefix, GetName()));
         }
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
   if (size >= 1) xAxis = (AliRsnMiniAxis *)fAxes[0];
   if (size >= 2) yAxis = (AliRsnMiniAxis *)fAxes[1];
   if (size >= 3) zAxis = (AliRsnMiniAxis *)fAxes[2];

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

   Int_t size = fAxes.GetEntries();
   AliInfo(Form("Sparse histogram name = '%s', with %d axes", name, size));

   // retrieve binnings and sizes of all axes
   // since the check for null values is done in Init(),
   // we assume that here they must all be well defined
   Int_t i, *nbins = new Int_t[size];
   for (i = 0; i < size; i++) {
      AliRsnMiniAxis *axis = (AliRsnMiniAxis *)fAxes[i];
      nbins[i] = axis->NBins();
   }

   // create fHSparseogram
   THnSparseF *h1 = new THnSparseF(name, "", size, nbins);

   // update the various axes using the definitions given in the array of axes here
   for (i = 0; i < size; i++) {
      AliRsnMiniAxis *axis = (AliRsnMiniAxis *)fAxes[i];
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
Bool_t AliRsnMiniOutput::FillMotherInAcceptance(const AliRsnMiniPair *pair, AliRsnMiniEvent *event, TClonesArray *valueList)
{
//
// Compute values for mother-based computations
//

   // check computation type
   if (fComputation != kMotherInAcc) {
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
Int_t AliRsnMiniOutput::FillPair(AliRsnMiniEvent *event1, AliRsnMiniEvent *event2, TClonesArray *valueList, Bool_t refFirst)
{
//
// Loops on the passed mini-event, and for each pair of particles
// which satisfy the charge and cut requirements defined here, add an entry.
// Returns the number of successful fillings.
// Last argument tells if the reference event for event-based values is the first or the second.
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
   AliRsnMiniParticle *p1, *p2;

   // it is necessary to know if criteria for the two daughters are the same
   // and if the two events are the same or not (mixing)
   //Bool_t sameCriteria = ((fCharge[0] == fCharge[1]) && (fCutID[0] == fCutID[1]));
   Bool_t sameCriteria = ((fCharge[0] == fCharge[1]) && (fDaughter[0] == fDaughter[1]));
   Bool_t sameEvent = (event1->ID() == event2->ID());

   TString selList1  = "";
   TString selList2  = "";
   Int_t   n1 = event1->CountParticles(fSel1, fCharge[0], fCutID[0]);
   Int_t   n2 = event2->CountParticles(fSel2, fCharge[1], fCutID[1]);
   for (i1 = 0; i1 < n1; i1++) selList1.Append(Form("%d ", fSel1[i1]));
   for (i2 = 0; i2 < n2; i2++) selList2.Append(Form("%d ", fSel2[i2]));
   AliDebugClass(1, Form("[%10s] Part #1: [%s] -- evID %6d -- charge = %c -- cut ID = %d --> %4d tracks (%s)", GetName(), (event1 == event2 ? "def" : "mix"), event1->ID(), fCharge[0], fCutID[0], n1, selList1.Data()));
   AliDebugClass(1, Form("[%10s] Part #2: [%s] -- evID %6d -- charge = %c -- cut ID = %d --> %4d tracks (%s)", GetName(), (event1 == event2 ? "def" : "mix"), event2->ID(), fCharge[1], fCutID[1], n2, selList2.Data()));
   if (!n1 || !n2) {
      AliDebugClass(1, "No pairs to mix");
      return 0;
   }

   // external loop
   for (i1 = 0; i1 < n1; i1++) {
      p1 = event1->GetParticle(fSel1[i1]);
      //p1 = event1->GetParticle(i1);
      //if (p1->Charge() != fCharge[0]) continue;
      //if (!p1->HasCutBit(fCutID[0])) continue;
      // define starting point for inner loop
      // if daughter selection criteria (charge, cuts) are the same
      // and the two events coincide, internal loop must start from
      // the first track *after* current one;
      // otherwise it starts from the beginning
      start = ((sameEvent && sameCriteria) ? i1 + 1 : 0);
      AliDebugClass(2, Form("Start point = %d", start));
      // internal loop
      for (i2 = start; i2 < n2; i2++) {
         p2 = event2->GetParticle(fSel2[i2]);
         //p2 = event2->GetParticle(i2);
         //if (p2->Charge() != fCharge[1]) continue;
         //if (!p2->HasCutBit(fCutID[1])) continue;
         // avoid to mix a particle with itself
         if (sameEvent && (p1->Index() == p2->Index())) {
            AliDebugClass(2, "Skipping same index");
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
            } else if (fPair.MotherPDG() != fMotherPDG) {
               continue;
            }
            Bool_t decayMatch = kFALSE;
            if (p1->PDGAbs() == AliRsnDaughter::SpeciesPDG(fDaughter[0]) && p2->PDGAbs() == AliRsnDaughter::SpeciesPDG(fDaughter[1]))
               decayMatch = kTRUE;
            if (p2->PDGAbs() == AliRsnDaughter::SpeciesPDG(fDaughter[0]) && p1->PDGAbs() == AliRsnDaughter::SpeciesPDG(fDaughter[1]))
               decayMatch = kTRUE;
            if (fCheckDecay && !decayMatch) continue;
	    if ( (fMaxNSisters>0) && (p1->NTotSisters()==p2->NTotSisters()) && (p1->NTotSisters()>fMaxNSisters)) continue;
	    if ( fCheckP &&(TMath::Abs(fPair.PmotherX()-(p1->Px(1)+p2->Px(1)))/(TMath::Abs(fPair.PmotherX())+1.e-13)) > 0.00001 && 	  
		          (TMath::Abs(fPair.PmotherY()-(p1->Py(1)+p2->Py(1)))/(TMath::Abs(fPair.PmotherY())+1.e-13)) > 0.00001 &&
     			  (TMath::Abs(fPair.PmotherZ()-(p1->Pz(1)+p2->Pz(1)))/(TMath::Abs(fPair.PmotherZ())+1.e-13)) > 0.00001 ) continue;
	    if ( fCheckFeedDown ){
	    		Int_t pdgGranma = 0;
	  		Bool_t isFromB=kFALSE;
	  		Bool_t isQuarkFound=kFALSE;
			
			if(fPair.IsFromB() == kTRUE) isFromB = kTRUE;
			if(fPair.IsQuarkFound() == kTRUE) isQuarkFound = kTRUE;
	  		if(fRejectIfNoQuark && !isQuarkFound) pdgGranma = -99999;
	  		if(isFromB){
	  		  if (!fKeepDfromB) pdgGranma = -9999; //skip particle if come from a B meson.
	  		}
	  		else{ 
	  		  if (fKeepDfromBOnly) pdgGranma = -999;
			  } 
	  		if (pdgGranma == -99999){
	  			AliDebug(2,"This particle does not have a quark in his genealogy\n");
	  			continue;
	  		}
	  		if (pdgGranma == -9999){
	  			AliDebug(2,"This particle come from a B decay channel but according to the settings of the task, we keep only the prompt charm particles\n");	
	  			continue;
	  		}	
	 
	  		if (pdgGranma == -999){
	  			AliDebug(2,"This particle come from a prompt charm particles but according to the settings of the task, we want only the ones coming from B\n");  
	  			continue;
	  		}	
		    }
         }
         // check pair against cuts
         if (fPairCuts) {
            if (!fPairCuts->IsSelected(&fPair)) continue;
         }
         // get computed values & fill histogram
         nadded++;
         if (refFirst) ComputeValues(event1, valueList); else ComputeValues(event2, valueList);
         FillHistogram();
      } // end internal loop
   } // end external loop

   AliDebugClass(1, Form("Pairs added in total = %4d", nadded));
   return nadded;
}
//___________________________________________________________
void AliRsnMiniOutput::SetDselection(UShort_t originDselection)
{
	// setting the way the D0 will be selected
	// 0 --> only from c quarks
	// 1 --> only from b quarks
	// 2 --> from both c quarks and b quarks
		
	fOriginDselection = originDselection;
	
	if (fOriginDselection == 0) {
		fKeepDfromB = kFALSE;
		fKeepDfromBOnly = kFALSE;
	}
	
	if (fOriginDselection == 1) {
		fKeepDfromB = kTRUE;
		fKeepDfromBOnly = kTRUE;
	}
	
	if (fOriginDselection == 2) {
		fKeepDfromB = kTRUE;
		fKeepDfromBOnly = kFALSE;
	}
	
	return;	
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
      AliRsnMiniAxis *axis = (AliRsnMiniAxis *)fAxes[i];
      if (!axis) {
         AliError("Null axis");
         continue;
      }
      ival = axis->GetValueID();
      if (ival < 0 || ival >= nval) {
         AliError(Form("Required value #%d, while maximum is %d", ival, nval));
         continue;
      }
      AliRsnMiniValue *val = (AliRsnMiniValue *)valueList->At(ival);
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
      ((TH1F *)obj)->Fill(fComputed[0]);
   } else if (obj->InheritsFrom(TH2F::Class())) {
      ((TH2F *)obj)->Fill(fComputed[0], fComputed[1]);
   } else if (obj->InheritsFrom(TH3F::Class())) {
      ((TH3F *)obj)->Fill(fComputed[0], fComputed[1], fComputed[2]);
   } else if (obj->InheritsFrom(THnSparseF::Class())) {
      THnSparseF *h = (THnSparseF *)obj;
      if (fCheckHistRange) {
         for (Int_t iAxis = 0; iAxis<h->GetNdimensions(); iAxis++) {
            if (fComputed.At(iAxis)>h->GetAxis(iAxis)->GetXmax() || fComputed.At(iAxis)<h->GetAxis(iAxis)->GetXmin()) return;
         }
      }
      h->Fill(fComputed.GetArray());
   } else {
      AliError("No output initialized");
   }
}
