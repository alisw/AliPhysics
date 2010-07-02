#include "TObjArray.h"
#include "AliLog.h"
#include "AliCDBRunRange.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBStorage.h"
#include "AliTRDcalibDB.h"
#include "AliTRDtrendingManager.h"

ClassImp(AliTRDtrendingManager)

AliTRDtrendingManager* AliTRDtrendingManager::fgInstance=NULL;
Bool_t AliTRDtrendingManager::fgTerminated = kFALSE;

//____________________________________________
AliTRDtrendingManager* AliTRDtrendingManager::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  //
  if (fgTerminated != kFALSE) return NULL;

  if (!fgInstance) {
    fgInstance = new AliTRDtrendingManager();
    AliTRDcalibDB *trd(AliTRDcalibDB::Instance());
    if(!trd){ 
      AliWarningGeneral("AliTRDtrendingManager", "TRD OCDB manager not initialized. No trending DB available.");
    } else {
      const TObjArray *trendMap(NULL/*trd->GetTrendMap()*/);
      if(!trendMap){
        AliWarningGeneral("AliTRDtrendingManager", "No TRD trending DB available  for TRD.");
      } else fgInstance->fEntries=(TObjArray*)trendMap->Clone();
    }
  }

  return fgInstance;
}

//____________________________________________
void AliTRDtrendingManager::Terminate()
{
  //
  // Singleton implementation
  // Deletes the instance of this class and sets the terminated flag,
  // instances cannot be requested anymore
  // This function can be called several times.
  //
  
  fgTerminated = kTRUE;
  
  if (fgInstance != NULL) {
    delete fgInstance;
    fgInstance = NULL;
  }
}

//____________________________________________
AliTRDtrendingManager::AliTRDtrendingManager() 
  : TObject()
  ,fEntries(NULL)
  ,fValue(NULL)
{
  fRunRange[0] = 0; fRunRange[1] = AliCDBRunRange::Infinity();
}

//____________________________________________
AliTRDtrendingManager::~AliTRDtrendingManager()
{
  if(fValue) delete fValue;
  if(fEntries) delete fEntries;
}

//____________________________________________
void AliTRDtrendingManager::AddValue(
  Char_t *class_name
  ,Char_t *name
  ,Char_t *title
  ,Double_t limits[2*(AliTRDtrendValue::kNlevels+1)]
  ,Char_t *messages[AliTRDtrendValue::kNlevels]
  ,Char_t *responsible
  ,Char_t *notifiables
  )
{
// Expert Function !!!
// Add a trend value to the map already loaded
// If no map loaded create a new one from scratch
//
// class_name : name of the performance task 
// name       : name of the value to be trended
// title      : description of the value to be trended
// limits     : array of acceptance limits for this value. The array is organized as follows : 
//               - field 0 and 1 normal limits
//               - field 2 and 3 first level of alarm
//               - ....
//               - field 8 and 9 fourth level of alarm
// messages   : array of alarm messages for each alarm level
// responsible: name and email of the responsible person. Format "name/email"
// notifiables: name and email of the notifiable persons. Format "name1/email1, name2/email2, etc"
//
  AliWarning("*** EXPERT FUNCTION *** This function is adding one trending value to the current DB. Continue if you know what yout do!");

  // create new trending value`
  if(!fValue) fValue = new AliTRDtrendValue(Form("%s_%s", class_name, name), title);
  else new(fValue) AliTRDtrendValue(Form("%s_%s", class_name, name), title);
  fValue->SetLimits(limits);
  for(Int_t ilevel(AliTRDtrendValue::kNlevels); ilevel--;) fValue->SetAlarm(ilevel, messages[ilevel]);
  TString s(responsible);
  TObjArray *r=s.Tokenize("/");
  if(r->GetEntriesFast()!=2){ 
    AliWarning("Responsible name/email incorrectly formated.");
  } else { 
    fValue->SetResponsible(((TObjString*)r->At(0))->String().Data(), ((TObjString*)r->At(1))->String().Data());
  }
  if(notifiables){
    s=notifiables;
    TObjArray *n=s.Tokenize(",");
    for(Int_t in(0); in<TMath::Min(AliTRDtrendValue::kNnotifiable, n->GetEntriesFast()); in++){
      TString ss(((TObjString*)n->At(in))->String());
      r=ss.Tokenize("/");
      if(r->GetEntriesFast()!=2){ 
        AliWarning(Form("Notifiable person name/email incorrectly formated for [%s].", ss.Data()));
      } else { 
        fValue->SetNotifiable(((TObjString*)r->At(0))->String().Data(), ((TObjString*)r->At(1))->String().Data());
      }
    }
  }
  // if no trending map defined create one
  if(!fEntries){
    AliInfo("No trending map loaded. Create one from scratch.");
    fEntries = new TObjArray(50);
    fEntries->SetOwner();
  }

  fEntries->AddLast(new AliTRDtrendValue(*fValue));
}

//____________________________________________
AliTRDtrendValue* AliTRDtrendingManager::GetValue(Char_t *class_name, Char_t *value_name)
{
  if(!fEntries){
    AliError("No trending map defined.");
    return NULL;
  }
  AliTRDtrendValue *val((AliTRDtrendValue*)fEntries->FindObject(Form("%s_%s", class_name, value_name)));
  if(!val){ 
    AliError(Form("Missing trending value %s [%s]", value_name, class_name));
    fEntries->ls();
  }
  return val;
}

//____________________________________________
Bool_t AliTRDtrendingManager::ModifyValue(
  Char_t *class_name
  ,Char_t *name
  ,Char_t *title
  ,Double_t *limits
  ,Char_t **messages
  ,Char_t *responsible
  ,Char_t *notifiables
  )
{
// Expert Function !!!
// Modify a trend value in the map already loaded
// see function AddValue() for explanation of input format. 

  if(!fEntries){
    AliError("No trending map loaded.");
    return kFALSE;
  }
  AliWarning("*** EXPERT FUNCTION *** This function is modifying one trending value to the current DB. Continue if you know what yout do!");

  AliTRDtrendValue *val((AliTRDtrendValue*)fEntries->FindObject(Form("%s_%s", class_name, name)));
  if(!val){ 
    AliError(Form("Missing trending value %s [%s]", name, class_name));
    return kFALSE;
  }  
  
  val->SetTitle(title);
  if(limits) val->SetLimits(limits);
  if(messages){ 
    for(Int_t ilevel(AliTRDtrendValue::kNlevels); ilevel--;) val->SetAlarm(ilevel, messages[ilevel]);
  }
  TString s;
  TObjArray *r(NULL);
  if(responsible){ 
    s=responsible;
    r=s.Tokenize("/");
    if(r->GetEntriesFast()!=2){ 
      AliWarning("Responsible name/email incorrectly formated.");
    } else { 
      val->SetResponsible(((TObjString*)r->At(0))->String().Data(), ((TObjString*)r->At(1))->String().Data());
    }
  }
  if(notifiables){
    s=notifiables;
    TObjArray *n=s.Tokenize(",");
    for(Int_t in(0); in<TMath::Min(AliTRDtrendValue::kNnotifiable, n->GetEntriesFast()); in++){
      TString ss(((TObjString*)n->At(in))->String());
      r=ss.Tokenize("/");
      if(r->GetEntriesFast()!=2){ 
        AliWarning(Form("Notifiable person name/email incorrectly formated for [%s].", ss.Data()));
      } else { 
        val->SetNotifiable(((TObjString*)r->At(0))->String().Data(), ((TObjString*)r->At(1))->String().Data());
      }
    }
  }
  return kTRUE;
}

//____________________________________________
void AliTRDtrendingManager::Print(Option_t *o) const
{
  if(!fEntries){
    AliError("No trending map available.");
    return;
  }

  for(Int_t iv(0); iv<fEntries->GetEntriesFast(); iv++){
    ((AliTRDtrendValue*)fEntries->At(iv))->Print(o);
  }
}

//____________________________________________
void AliTRDtrendingManager::Save()
{
// Saving TRD trending DB to $ALICE_ROOT/OCDB.

  AliWarning("Saving TRD trending DB to $ALICE_ROOT/OCDB.");

  AliCDBMetaData *metaData= new AliCDBMetaData(); 
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Alexander Wilk");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-21-01"); //root version
  metaData->SetComment("TRD trending ");
  
  AliCDBId id("TRD/Calib/Trend", fRunRange[0], fRunRange[1]); 
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBStorage *gStorLoc = man->GetStorage("local://$ALICE_ROOT/OCDB");
  if (!gStorLoc) {
    return;
  }
  gStorLoc->Put(fEntries, id, metaData); 
}
