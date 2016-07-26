#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TObjArray.h>
#include <TFile.h>
#include <TString.h>
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliTPCRecoParam.h"
#include "TMap.h"
#endif

void ModifyObject(TObject* obj);

void applyCorrTPCRecoParam(int useRun=245231,int firstRun=0, int lastRun=-1,
			   TString comment="Cheb.maps used instead of composed correction",
			   const char* destOCDB="local://./",
			   char* cdbPth = "TPC/Calib/RecoParam",
			   TString objIniPath="" // initial TRD alignment, if path, take from specific storage
			   )
{
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  man->SetRun( useRun );
  //
  AliCDBEntry* entry = 0;
  if (objIniPath.EndsWith(".root")) { 
    TFile* flInp = TFile::Open(objIniPath.Data());
    entry = (AliCDBEntry*)flInp->Get("AliCDBEntry");
    flInp->Close();
    delete flInp;
  }
  else {
    if (!objIniPath.IsNull()) man->SetSpecificStorage(cdbPth,objIniPath.Data());
    entry = man->Get(cdbPth);
  }
  TObject *obj = entry->GetObject();
  entry->SetObject(0); // detach
  //
  ModifyObject(obj);
  //
  // store new object
  AliCDBMetaData* mdold = entry->GetMetaData();
  AliCDBMetaData* mdnew = new AliCDBMetaData();
  TString commComb = "";
  if (mdold) {
    mdnew->SetResponsible(mdold->GetResponsible());
    mdnew->SetBeamPeriod(mdold->GetBeamPeriod());
    mdnew->SetAliRootVersion(mdold->GetAliRootVersion());
    commComb += mdold->GetComment();
    commComb += " ";
  }
  commComb += comment;
  mdnew->SetComment(commComb.Data());
  //
  obj->Print();
  //
  man->SetSpecificStorage(cdbPth,destOCDB);
  AliCDBId id(cdbPth,firstRun>=0 ? firstRun : 0,lastRun>=0? lastRun : AliCDBRunRange::Infinity() );
  man->Put(obj,id,mdnew); 
  //
}


void ModifyObject(TObject* obj)
{
  TObjArray* arr = dynamic_cast<TObjArray*>(obj);
  TIter next(arr);
  AliTPCRecoParam* par = 0;
  
  // vector for charging-up masking
  TVectorF zreg(1),zregsigInv(1);
  zreg[0] = -5.;
  zregsigInv[0] = 1./2.0; // 2cm sigma in Z

  while( (par=(AliTPCRecoParam*)next())) {
    int spec = par->GetEventSpecie();
    //    if (spec&AliRecoParam::kLowMult || spec&AliRecoParam::kHighMult) 
    {
      par->SetUseCorrectionMap(kTRUE);
      par->SetUseComposedCorrection(kFALSE);
      par->SetUseAlignmentTime(kFALSE); 
      par->SetAccountDistortions(1);
      //
      // add difference to reference distortions as syst. error
      par->SetUseDistortionFractionAsErrorY(0.3);
      par->SetUseDistortionFractionAsErrorZ(0.3);
      //
      // add difference to reference dispersion as syst. error
      par->SetUseDistDispFractionAsErrorY(3.0);
      par->SetUseDistDispFractionAsErrorZ(3.0);

      // charging-up zone errors
      par->SetSystErrClInnerRegZ( new TVectorF(zreg) );
      par->SetSystErrClInnerRegZSigInv( new TVectorF(zregsigInv) );
      ((double*)par->GetSystematicErrorClusterInner())[0] = 1.0; // 1.5 cm err at maximum
      ((double*)par->GetSystematicErrorClusterInner())[1] = 5.0; // dumped with 5 cm
    }
  }
  //    
  printf("Mod\n"); obj->Print();
}
