//
// *** Class AliRsnTrainManager ***
//
//  Base class for Action
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Jan Musinsky (jan.musinsky@cern.ch)
//

#include <TError.h>
#include <TObjString.h>
#include <TMap.h>
#include <TClass.h>

#include "AliRsnTrainManager.h"

ClassImp(AliRsnTrainManager)

AliRsnTrainManager *AliRsnTrainManager::fgRsnTrainManager = 0;

//__________________________________________________________________________________________________
AliRsnTrainManager::AliRsnTrainManager(const char *name,const char *title) : TNamed(name,title),
   fGlobals(0)
{
   //
   // Default constructor
   //

   fgRsnTrainManager = this;

   if (TClass::IsCallingNew() != TClass::kDummyNew) {
      fGlobals    = new TMap();
   }

}

//__________________________________________________________________________________________________
AliRsnTrainManager::AliRsnTrainManager(const AliRsnTrainManager &copy) : TNamed(copy),
   fGlobals(copy.fGlobals)
{
   //
   // Copy constructor
   //
}

//__________________________________________________________________________________________________
AliRsnTrainManager &AliRsnTrainManager::operator=(const AliRsnTrainManager &copy)
{
   //
   // Assignment constructor
   //
   TNamed::operator=(copy);
   if (this == &copy)
      return *this;

   fGlobals = copy.fGlobals;
   return (*this);
}

//__________________________________________________________________________________________________
AliRsnTrainManager::~AliRsnTrainManager()
{
   //
   // Destructor
   //
   delete fGlobals;
   fgRsnTrainManager = 0;
}

//______________________________________________________________________________
void AliRsnTrainManager::SetGlobalStr(const char *key, const char *value,Bool_t verbose)
{
   // Define a custom string variable mapped to a global unique name. The variable
   // can be then retrieved by a given analysis macro via GetGlobalStr(key).
   AliRsnTrainManager *mgr = AliRsnTrainManager::GetRsnTrainManager();
   if (!mgr) {
      ::Error("AliRsnTrainManager::SetGlobalStr", "No analysis manager defined");
      return;
   }
   Bool_t valid = kFALSE;
   TString existing = AliRsnTrainManager::GetGlobalStr(key, valid);
   if (valid) {
      if (verbose) ::Error("AliRsnTrainManager::SetGlobalStr", "Global %s = %s already defined.", key, existing.Data());
      return;
   }
   mgr->GetGlobals()->Add(new TObjString(key), new TObjString(value));
}

//______________________________________________________________________________
const char *AliRsnTrainManager::GetGlobalStr(const char *key, Bool_t &valid)
{
   // Static method to retrieve a global variable defined via SetGlobalStr.
   valid = kFALSE;
   AliRsnTrainManager *mgr = AliRsnTrainManager::GetRsnTrainManager();
   if (!mgr) return 0;
   TObject *value = mgr->GetGlobals()->GetValue(key);
   if (!value) return 0;
   valid = kTRUE;
   return value->GetName();
}

//______________________________________________________________________________
void AliRsnTrainManager::SetGlobalInt(const char *key, Int_t value,Bool_t verbose)
{
   // Define a custom integer variable mapped to a global unique name. The variable
   // can be then retrieved by a given analysis macro via GetGlobalInt(key).
   AliRsnTrainManager *mgr = AliRsnTrainManager::GetRsnTrainManager();
   if (!mgr) {
      ::Error("AliRsnTrainManager::SetGlobalStr", "No analysis manager defined");
      return;
   }
   Bool_t valid = kFALSE;
   Int_t existing = AliRsnTrainManager::GetGlobalInt(key, valid);
   if (valid) {
      if (verbose) ::Error("AliRsnTrainManager::SetGlobalInt", "Global %s = %i already defined.", key, existing);
      return;
   }
   mgr->GetGlobals()->Add(new TObjString(key), new TObjString(TString::Format("%i",value)));
}

//______________________________________________________________________________
Int_t AliRsnTrainManager::GetGlobalInt(const char *key, Bool_t &valid)
{
   // Static method to retrieve a global variable defined via SetGlobalInt.
   valid = kFALSE;
   AliRsnTrainManager *mgr = AliRsnTrainManager::GetRsnTrainManager();
   if (!mgr) return 0;
   TObject *value = mgr->GetGlobals()->GetValue(key);
   if (!value) return 0;
   valid = kTRUE;
   TString s = value->GetName();
   return s.Atoi();
}

//______________________________________________________________________________
void AliRsnTrainManager::SetGlobalDbl(const char *key, Double_t value,Bool_t verbose)
{
   // Define a custom double precision variable mapped to a global unique name. The variable
   // can be then retrieved by a given analysis macro via GetGlobalInt(key).
   AliRsnTrainManager *mgr = AliRsnTrainManager::GetRsnTrainManager();
   if (!mgr) {
      ::Error("AliRsnTrainManager::SetGlobalStr", "No analysis manager defined");
      return;
   }
   Bool_t valid = kFALSE;
   Double_t existing = AliRsnTrainManager::GetGlobalDbl(key, valid);
   if (valid) {
      if (verbose) ::Error("AliRsnTrainManager::SetGlobalInt", "Global %s = %g already defined.", key, existing);
      return;
   }
   mgr->GetGlobals()->Add(new TObjString(key), new TObjString(TString::Format("%f",value)));
}

//______________________________________________________________________________
Double_t AliRsnTrainManager::GetGlobalDbl(const char *key, Bool_t &valid)
{
   // Static method to retrieve a global variable defined via SetGlobalDbl.
   valid = kFALSE;
   AliRsnTrainManager *mgr = AliRsnTrainManager::GetRsnTrainManager();
   if (!mgr) return 0;
   TObject *value = mgr->GetGlobals()->GetValue(key);
   if (!value) return 0;
   valid = kTRUE;
   TString s = value->GetName();
   return s.Atof();
}

//______________________________________________________________________________
void AliRsnTrainManager::SetGlobalObj(const char *key, TObject *value,Bool_t verbose)
{
   // Define a custom double precision variable mapped to a global unique name. The variable
   // can be then retrieved by a given analysis macro via GetGlobalInt(key).
   AliRsnTrainManager *mgr = AliRsnTrainManager::GetRsnTrainManager();
   if (!mgr) {
      ::Error("AliRsnTrainManager::SetGlobalStr", "No analysis manager defined");
      return;
   }
   Bool_t valid = kFALSE;
   TObject *existing = AliRsnTrainManager::GetGlobalObj(key, valid);
   if (valid) {
      if (verbose) ::Error("AliRsnTrainManager::SetGlobalInt", "Global %s = %p already defined.", key, existing);
      return;
   }
   mgr->GetGlobals()->Add(new TObjString(key), value->Clone());
}

//______________________________________________________________________________
TObject *AliRsnTrainManager::GetGlobalObj(const char *key, Bool_t &valid)
{
   // Static method to retrieve a global variable defined via SetGlobalDbl.
   valid = kFALSE;
   AliRsnTrainManager *mgr = AliRsnTrainManager::GetRsnTrainManager();
   if (!mgr) return 0;
   TObject *value = mgr->GetGlobals()->GetValue(key);
   return value;
}

//______________________________________________________________________________
void AliRsnTrainManager::Print(Option_t */*option*/) const {



   if (!fGlobals) return;
   Printf("=========================================================");

   Printf("            Rsn Global variables: \n");
   TIter next(fGlobals->GetTable());
   TPair *a;
   TObjString *str;
   while ((a = (TPair *)next())) {
      str = 0;
//    if (a->Key()&& a->Key()->IsOnHeap()){
//    }
      if (a->Value() && a->Value()->IsOnHeap()) {
         str = dynamic_cast<TObjString *>(a->Value());
      }
      if (str) {
         if (str) Printf("  %20s = %-20s",a->Key()->GetName(),str->GetString().Data());
      }
   }
   Printf("\n=========================================================");

   /* Printf("            Rsn Particles and Cuts \n");


      next.Reset();
   // AliRsnCut
      while ((a = (TPair *)next())) {
         str = 0;
   //    if (a->Key()&& a->Key()->IsOnHeap()){
   //    }
         if (a->Value() && a->Value()->IsOnHeap()){
            str = dynamic_cast<TObjString *>(a->Value());
         }
         if (str) {
            if (str) Printf("  %20s = %-20s",a->Key()->GetName(),str->GetString().Data());
         }
      }

      Printf("\n=========================================================");

   */
}
