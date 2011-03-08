//
// Class AliRsnCutSet
//
// This is the front-end for cut management and checking.
// It must be prepared by adding all required single cuts,
// and then with a logical expression which combines all cuts
// with the "AND", "OR" and "NOT" operators.
//

#include "AliLog.h"

#include "AliRsnExpression.h"

#include "AliRsnCutSet.h"

ClassImp(AliRsnCutSet)

//_____________________________________________________________________________
AliRsnCutSet::AliRsnCutSet() :
   AliRsnTarget(),
   fCuts(0),
   fNumOfCuts(0),
   fCutScheme(""),
   fCutSchemeIndexed(""),
   fBoolValues(0),
   fIsScheme(kFALSE),
   fExpression(0)
{
//
// Constructor without name (not recommended)
//

   fBoolValues = new Bool_t[1];
   AliRsnExpression::fgCutSet = this;
}

//_____________________________________________________________________________
AliRsnCutSet::AliRsnCutSet(const char *name, RSNTARGET target) :
   AliRsnTarget(name, target),
   fCuts(0),
   fNumOfCuts(0),
   fCutScheme(""),
   fCutSchemeIndexed(""),
   fBoolValues(0),
   fIsScheme(kFALSE),
   fExpression(0)
{
//
// Constructor with argument name (recommended)
//

   fBoolValues = new Bool_t[1];
   fExpression = 0;
   AliRsnExpression::fgCutSet = this;
}

//_____________________________________________________________________________
AliRsnCutSet::AliRsnCutSet(const AliRsnCutSet & copy) :
   AliRsnTarget(copy),
   fCuts(copy.fCuts),
   fNumOfCuts(copy.fNumOfCuts),
   fCutScheme(copy.fCutScheme),
   fCutSchemeIndexed(copy.fCutSchemeIndexed),
   fBoolValues(0),
   fIsScheme(copy.fIsScheme),
   fExpression(copy.fExpression)
{
//
// Copy constructor
//

   Int_t i;
   fBoolValues = new Bool_t[fNumOfCuts];
   for (i = 0; i < fNumOfCuts; ++i) {
      fBoolValues[i] = copy.fBoolValues[i];
   }

   AliRsnExpression::fgCutSet = this;
}

//_____________________________________________________________________________
AliRsnCutSet& AliRsnCutSet::operator=(const AliRsnCutSet & copy)
{
//
// Assignment operator.
//

   AliRsnTarget::operator=(copy);

   fCuts = copy.fCuts;
   fNumOfCuts = copy.fNumOfCuts;
   fCutScheme = copy.fCutScheme;
   fCutSchemeIndexed = copy.fCutSchemeIndexed;
   fIsScheme = copy.fIsScheme;
   fExpression = copy.fExpression;

   if (fBoolValues) delete [] fBoolValues;

   Int_t i;
   fBoolValues = new Bool_t[fNumOfCuts];
   for (i = 0; i < fNumOfCuts; ++i) {
      fBoolValues[i] = copy.fBoolValues[i];
   }

   AliRsnExpression::fgCutSet = this;

   return (*this);
}

//_____________________________________________________________________________
AliRsnCutSet::~AliRsnCutSet()
{
//
// Destructor
//

   delete fExpression;
   delete [] fBoolValues;
}

//_____________________________________________________________________________
void AliRsnCutSet::AddCut(AliRsnCut *cut)
{
//
// Add a new cut.
// This must be done for all components of the final expression.
// If the target of the cut does not match the target of this,
// the cut is not added.
//

   if (!cut->IsTarget(GetTargetType())) {
      AliError(Form("Cannot add this cut (cut set name,target = [%s],[%s] --- cut name,target = [%s],[%s]", GetName(), GetTargetTypeName(), cut->GetName(), cut->GetTargetTypeName()));
      return;
   }

   Int_t i;

   AliDebug(AliLog::kDebug, "<-");
   fCuts.Add(cut);
   AliInfo(Form("====> Adding a new cut: [%s]", cut->GetName()));
   //cut->Print();
   fNumOfCuts++;

   if (fBoolValues) delete [] fBoolValues;

   fBoolValues = new Bool_t[fNumOfCuts];
   for (i = 0; i < fNumOfCuts; i++) {
      fBoolValues[i] = kTRUE;
   }

   AliDebug(AliLog::kDebug, Form("%d", fCuts.GetEntriesFast()));
   AliDebug(AliLog::kDebug, "->");
}

//_____________________________________________________________________________
void AliRsnCutSet::ShowCuts() const
{
//
// Prints all cuts
//
   AliRsnCut *cut;

   for (Int_t i = 0; i < fCuts.GetEntriesFast() ; i++) {
      cut = (AliRsnCut*)fCuts.At(i);
      cut->Print();
   }
}

//_____________________________________________________________________________
Bool_t AliRsnCutSet::IsSelected(TObject *object)
{
//
// Checks an object according to the cut expression defined here.
//

   Int_t i;

   if (!fNumOfCuts) return kTRUE;

   Bool_t boolReturn = kTRUE;
   AliRsnCut *cut;
   for (i = 0; i < fNumOfCuts; i++) {
      cut = (AliRsnCut*)fCuts.At(i);
      fBoolValues[i] = cut->IsSelected(object);
   }

   if (fIsScheme) boolReturn = Passed();
   return boolReturn;
}

//_____________________________________________________________________________
void AliRsnCutSet::SetCutScheme(const char *theValue)
{
//
// Define the combination of cuts with logical operators
// and using the names given to all defined cuts.
//

   AliDebug(AliLog::kDebug, "<-");

   fCutScheme = theValue;
   SetCutSchemeIndexed(theValue);
   fIsScheme = kTRUE;
   AliDebug(AliLog::kDebug, "->");
}

//_____________________________________________________________________________
void AliRsnCutSet::SetCutSchemeIndexed(TString theValue)
{
//
// Internal method which indexes all cuts to organize their combo
//

   AliDebug(AliLog::kDebug, "<-");
   theValue.Append(" ");
   // fCutSchemeIndexed = theValue;
   fCutSchemeIndexed = GetCutSchemeIndexed();
   AliDebug(AliLog::kDebug, "->");
}

//_____________________________________________________________________________
Int_t AliRsnCutSet::GetIndexByCutName(TString s)
{
//
// Retrieve the cut index from its name
//

   Int_t i;
   AliRsnCut *cut;

   for (i = 0; i < fCuts.GetEntriesFast(); i++) {
      cut = (AliRsnCut*) fCuts.At(i);
      if (!s.CompareTo(cut->GetName())) return i;
   }

   return -1;
}

//_____________________________________________________________________________
Bool_t AliRsnCutSet::Passed()
{
//
// Combines the cuts according to expression
// and gives a global response to the cut check
//

   AliRsnExpression::fgCutSet = this;
   if (!fExpression) {
      fExpression = new AliRsnExpression(fCutSchemeIndexed);
      AliDebug(AliLog::kDebug, "fExpression was created.");
   }

   if (fCuts.IsEmpty()) return kTRUE;

   return fExpression->Value(*GetCuts());
}

//_____________________________________________________________________________
Bool_t AliRsnCutSet::IsValidScheme()
{
//
// Validity check on cut expression specified by user
//

   TString str(fCutScheme);
   AliRsnCut *cut;
   for (Int_t i = 0; i < fNumOfCuts; i++) {
      cut = (AliRsnCut*)fCuts.At(i);
      str.ReplaceAll(cut->GetName(), "");
   }
   str.ReplaceAll("&", "");
   str.ReplaceAll("!", "");
   str.ReplaceAll("|", "");
   str.ReplaceAll("(", "");
   str.ReplaceAll(")", "");

   if (!str.IsNull()) {
      AliError(Form("Cut scheme '%s' is not valid !!!", fCutScheme.Data()));
      return kFALSE;
   }

   return kTRUE;
//   return (!(ShowCutScheme().Contains("Error")));
}

//_____________________________________________________________________________
TString AliRsnCutSet::ShowCutScheme()
{
//
// Utility method to check validity of expression
//

   return fCutScheme;
//   return fExpression->Unparse();
}

//_____________________________________________________________________________
Int_t AliRsnCutSet::TestExpression(TString opt)
{

//   AliRsnCut *cut1 = new AliRsnCut ("aaa","aaa");
//   cut1->SetCutValues (AliRsnCut::kEsdPt,0.0,1.0);
//   AliRsnCut *cut2 = new AliRsnCut ("bbb","bbb");
//   cut2->SetCutValues (AliRsnCut::kEsdPt,1.,2.0);
//   AliRsnCut *cut3 = new AliRsnCut ("ccc","ccc");
//   cut3->SetCutValues (AliRsnCut::kEsdPt,2.0,3.0);
//
//   AliRsnCutSet* set  = new AliRsnCutSet ("setOne");
//   set->AddCut (cut1);
//   set->AddCut (cut2);
//   set->AddCut (cut3);
//
//   set->SetCutScheme ("(aaa&!(ccc))&(bbb&!(ccc))");
//
//   set->ShowCuts ();

   AliDebug(1, opt.Data());
   return 0;
}

//_____________________________________________________________________________
void AliRsnCutSet::PrintSetInfo()
{
//
// Show data about the cut set
//

   Int_t i;

   AliInfo("========== Rsn Cut Set info ==============");
   AliInfo(Form("Scheme : %s", fCutScheme.Data()));
   AliInfo(Form("Scheme : %s", fCutSchemeIndexed.Data()));
   AliInfo(Form("Num of Cuts: %d", fCuts.GetEntriesFast()));
   AliInfo("====== Cuts ======");
   AliRsnCut *cut;
   for (i = 0; i < fCuts.GetEntriesFast(); i++) {
      cut = (AliRsnCut*) fCuts.At(i);
      if (cut) AliInfo(Form("%d %d", i, fBoolValues[i]));
   }
   AliInfo("========== END Rsn Cut Mgr info ==============");
}

//_____________________________________________________________________________
TString AliRsnCutSet::GetCutSchemeIndexed()
{
//
// Internal method to retrieve the list of cuts with their indexes
// for evaluation of cut expression
//

   AliDebug(AliLog::kDebug, "<-");
   Int_t i;
   TString str(fCutScheme);
   AliDebug(AliLog::kDebug, Form("Num of cuts %d", fCuts.GetEntriesFast()));
   AliRsnCut *cut;
   for (i = 0; i < fCuts.GetEntriesFast(); i++) {
      cut = (AliRsnCut*) fCuts.At(i);
      str.ReplaceAll(cut->GetName(), Form("%d", i));
   }
   AliDebug(AliLog::kDebug, "->");
   return str;
}
