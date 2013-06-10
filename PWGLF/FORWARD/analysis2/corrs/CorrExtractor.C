#include <TFile.h>
#include <TList.h>
#include <TParameter.h>
#include <TError.h>
#include "AliCorrectionManagerBase.h"

struct CorrExtractor 
{
  CorrExtractor(AliCorrectionManagerBase* manager)
    : fFile(0), 
      fTop(0), 
      fOut(""),
      fRunNo(0), 
      fSys(0), 
      fSNN(0), 
      fField(0), 
      fMC(false), 
      fSatellite(false),
      fManager(manager)
  {}
  TCollection* GetCollection(TCollection* p, 
			     const TString& name)
  {
    TObject* o = 0;
    if (p == 0) { 
      o = fFile->Get(name);
      if (!o) { 
	Warning("CorrExtractor", "Object %s not found in file", name.Data());
	return 0;
      }
    }
    else {
      o = p->FindObject(name);
      if (!o) { 
	Warning("CorrExtractor", "Object %s not found in %s", 
		name.Data(), p->GetName());
	return 0;
      }
    }      
    if (!o->IsA()->InheritsFrom(TCollection::Class())) {
      Warning("CorrExtractor", "%s in %s is not a collection, but a %s", 
	      name.Data(), (p ? p->GetName() : "file"), o->ClassName());
      return 0;
    }
    return static_cast<TCollection*>(o);
  }
  TCollection* FindCollection(const TString& path)
  {
    if (path.IsNull()) return 0;
    TObjArray*   tokens = path.Tokenize("/");
    TIter        next(tokens);
    TObjString*  token = 0;
    TCollection* p     = 0;
    while ((token = static_cast<TObjString*>(next()))) {
      const TString& t = token->String();
      if (t.IsNull()) continue;
      p = GetCollection(p, t);
      if (!p) break;
    }
    tokens->Delete();
    return p;
    
  }
  
  TObject* FindObject(const TString& path, 
		      const TString& name) 
  {
    if (path.IsNull()) { 
      TObject* o = fFile->Get(name);
      if (!o) { 
	Warning("CorrExtractor", "Object %s not found in file", 
		name.Data());
	return 0;
      }
      return o;
    }
    TCollection* p     = FindCollection(path);
    if (!p) { 
      Warning("CorrExtractor", "Path %s invalid", path.Data());
      return 0;
    }
    return p->FindObject(name);
  }

  Bool_t Init(const TString&        fileName, 
	      const TString&        sumFolder, 
	      const TString&        out)
  {
    fOut  = out;
    Clear();

    fFile = TFile::Open(fileName, "READ");
    if (!fFile) {
      Error("CorrExtractor", "Failed to open \"%s\"", fileName.Data());
      Clear();
      return false;
    }
    TCollection* c = FindCollection(Form("%s/fmdEventInspector", 
					 sumFolder.Data()));
    if (!c) { 
      Error("CorrExtractor", "Couldn't get event inspector list from %s",
	    fileName.Data());
      Clear();
      return false;
    }
    TObject* oSys        = c->FindObject("sys");
    TObject* oSNN        = c->FindObject("sNN");
    TObject* oFld        = c->FindObject("field");
    TObject* oRun        = c->FindObject("runNo");
    TObject* oSat        = c->FindObject("satellite");
    if (oSys) fSys       = oSys->GetUniqueID();
    if (oSNN) fSNN       = oSNN->GetUniqueID();
    if (oFld) fField     = oFld->GetUniqueID();
    if (oRun) fRunNo     = oRun->GetUniqueID();
    if (oSat) fSatellite = oSat->GetUniqueID();

    if (fSys   <= 0 || fSys > 3 ||
	fSNN   <= 0 || 
	fRunNo <= 0) {
      Error("CorrExtractor", "Failed to get settings");
      Clear();
      return false;
    } 
    return true;
  }
  void SetMC(Bool_t mc=true) { fMC = mc; }
  Bool_t Extract(const TClass* cls, const TString& parent)
  {
    return Extract(cls->GetName(), parent);
  }
  Bool_t Extract(const TString& objName, 
		 const TString& parent="") 
  {
    if (!fFile) { 
      Warning("Extract", "No file opened");
      return false;
    }   
    TObject* o = FindObject(parent, objName);
    if (!o) { 
      Warning("Extract", "Object %s not found in collection %s", 
	      objName.Data(), parent.Data());
      return false;
    }
    return fManager->Store(o, 
			   fRunNo, 
			   fSys, 
			   fSNN, 
			   fField, 
			   fMC, 
			   fSatellite, 
			   fOut.Data());
  }

  void Clear()
  {
    if (fFile) fFile->Close();
    fFile      = 0;
    fTop       = 0;
    fRunNo     = 0;
    fSys       = 0;
    fSNN       = 0;
    fField     = 0;
    fMC        = false;
    fSatellite = false;
  }
  TFile*                    fFile;
  TList*                    fTop;
  TString                   fOut;
  ULong_t                   fRunNo;
  UShort_t                  fSys; 
  UShort_t                  fSNN;
  Short_t                   fField;
  Bool_t                    fMC;
  Bool_t                    fSatellite;
  AliCorrectionManagerBase* fManager;
};

//
// EOF
//


