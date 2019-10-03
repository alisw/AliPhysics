#include <TFile.h>
#include <TList.h>
#include <TParameter.h>
#include <TError.h>
#include "AliCorrectionManagerBase.h"

/**
 * Extract corrections from result file 
 * 
 */
struct CorrExtractor 
{
  /** 
   * Constructor 
   * 
   * @param manager Correction manager
   */
  CorrExtractor(AliCorrectionManagerBase* manager)
    : fFile(0), 
      fTop(0), 
      fOut(""),
      fRunNo(0), 
      fSys(0), 
      fSNN(0), 
      fField(999), 
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
  /** 
   * Find a collection in a file 
   * 
   * @param path Path to collection
   * 
   * @return Found collection or null
   */
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
  /** 
   * Find an object 
   * 
   * @param path Path to object 
   * @param name Name of object
   * 
   * @return Found object or null
   */  
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
  /** 
   * Initialize this extactor
   * 
   * @param fileName  File to extract from 
   * @param sumFolder The summed folder 
   * @param out       The result folder
   * 
   * @return true on success
   */
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
      // sumFolder = "forwardQAResults";
      c = FindCollection(Form("%s/fmdEventInspector", 
			      "forwardQAResults"));
      if (!c) {
	Error("CorrExtractor", "Couldn't get event inspector list from %s",
	      fileName.Data());
	Clear();
	return false;
      }
    }
    TObject* oSys        = c->FindObject("sys");
    TObject* oSNN        = c->FindObject("sNN");
    TObject* oFld        = c->FindObject("field");
    TObject* oRun        = c->FindObject("runNo");
    TObject* oSat        = c->FindObject("satellite");
    if (oSys && fSys   <= 0)   fSys       = oSys->GetUniqueID();
    if (oSNN && fSNN   <= 0)   fSNN       = oSNN->GetUniqueID();
    if (oFld && fField >= 999) fField     = oFld->GetUniqueID();
    if (oRun && fRunNo <= 0)   fRunNo     = oRun->GetUniqueID();
    if (oSat)                  fSatellite = oSat->GetUniqueID();

    Bool_t ret = true;
    if (fSys <= 0 || fSys > 5) {
      Error("CorrExtractor", "Invalid collision energy: %d", fSys);
      ret = false;
    }
    if (fSNN <= 0) {
      Error("CorrExtractor", "Invalid collision energy: %d", fSNN);
      ret = false;
    }
    if (fField >= 999) {
      Error("CorrExtractor", "Invalid field value: %d", fField);
      ret = false;
    }
    if (fRunNo <= 0 ){
      Error("CorrExtractor", "Invalid run number: %d", fRunNo);
      ret = false;
    }
    if (!ret) Clear();
    return ret;
  }
  /** 
   * Set whether this is MC or not
   * 
   * @param mc If true, consider this MC 
   */
  void SetMC(Bool_t mc=true) { fMC = mc; }
  /** 
   * Extract the stuff 
   * 
   * @param cls    Class of object
   * @param parent Parent folder 
   * 
   * @return 
   */
  Bool_t Extract(const TClass* cls, const TString& parent)
  {
    return Extract(cls->GetName(), parent);
  }
  /** 
   * Extract the stuff
   * 
   * @param objName  Object name 
   * @param parent   Parent folder 
   * 
   * @return 
   */
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
    if (o->TestBit(1<<15) && !o->TestBit(1<<16)) {
      Warning("Extract", "Object %s is not good", objName.Data());
      TFile* bad = TFile::Open("bad.root", "RECREATE");
      (new TNamed("bad","BadCorrection"))->Write();
      bad->Write();
      bad->Close();
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
  /** 
   * Clear this extractor 
   * 
   */
  void Clear()
  {
    if (!fFile) return;
    fFile->Close();
    fFile      = 0;
    fTop       = 0;
    fRunNo     = 0;
    fSys       = 0;
    fSNN       = 0;
    fField     = 999;
    fMC        = false;
    fSatellite = false;
  }
  TFile*                    fFile;          // Our file
  TList*                    fTop;           // Top list
  TString                   fOut;           // Output 
  ULong_t                   fRunNo;         // Run number
  UShort_t                  fSys;           // System
  UShort_t                  fSNN;           // Collision energy in GeV
  Short_t                   fField;         // L3 field in kG
  Bool_t                    fMC;            // Simulation flag
  Bool_t                    fSatellite;     // Satellite interaction flag
  AliCorrectionManagerBase* fManager;       // Correction manager to use 
};

//
// EOF
//


