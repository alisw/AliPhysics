/**********************************************
 *
 * Class designed to store all multiplicity
 * variables in a TClonesArray
 *
 * Instancing an empty AliMultInput class
 * will allow you to store standard variables
 *
 **********************************************/

#include "TList.h"
#include "TProfile.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliOADBMultSelection.h"
#include <TROOT.h>
#include <TMap.h>

ClassImp(AliMultInput);

const TString AliMultInput::VarName[kNVariables] = {
  "fAmplitude_V0A", "fAmplitude_V0A1", "fAmplitude_V0A2", "fAmplitude_V0A3", "fAmplitude_V0A4",
  "fAmplitude_V0C", "fAmplitude_V0C1", "fAmplitude_V0C2", "fAmplitude_V0C3", "fAmplitude_V0C4",
  "fAmplitude_V0Apartial", "fAmplitude_V0Cpartial", "fAmplitude_V0AEq", "fAmplitude_V0CEq",
  "fAmplitude_OnlineV0A","fAmplitude_OnlineV0C", "fAmplitude_V0AADC", "fAmplitude_V0CADC", "fFlatenicity_V0", 
  "fnSPDClusters", "fnSPDClusters0", "fnSPDClusters1",
  "fMultiplicity_ADA", "fMultiplicity_ADC",
  "fRefMultEta5", "fRefMultEta8", "fnTracklets", "fnTracklets08", "fnTracklets15",
  "fZncEnergy","fZpcEnergy","fZnaEnergy","fZpaEnergy", "fZem1Energy", "fZem2Energy",
  "fZnaTower", "fZncTower", "fZpaTower", "fZpcTower",
  "fZnaFired", "fZncFired", "fZpaFired", "fZpcFired",
  "fNTracks", "fNTracksTPCout", "fNTracksGlobal2015", "fNTracksGlobal2015Trigger", "fNTracksITSsa2010",
  "fNTracksINELgtONE","fNPartINELgtONE","fEvSel_VtxZ",
  "fMC_NPart","fMC_NColl","fMC_NchV0A","fMC_NchV0C",
  "fMC_NchEta05","fMC_NchEta08","fMC_NchEta10","fMC_NchEta14",
  "fMC_b", "fSpherocityMC", "fSpherocityTracksMC"
};

const Bool_t AliMultInput::VarIsInteger[kNVariables] = {
  kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
  kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
  kFALSE, kFALSE, kFALSE, kFALSE,
  kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
  kTRUE, kTRUE, kTRUE,
  kFALSE, kFALSE,
  kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
  kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
  kFALSE, kFALSE, kFALSE, kFALSE,
  kTRUE, kTRUE, kTRUE, kTRUE,
  kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
  kFALSE, kFALSE, kFALSE,
  kTRUE, kTRUE, kTRUE, kTRUE,
  kTRUE, kTRUE, kTRUE, kTRUE,
  kFALSE, kFALSE, kFALSE
};

AliMultInput::AliMultInput() :
TNamed(),
fNVars(0), fVariableList(0x0),
fNVtxZ(0), fVariableVertexZList(0x0),
fMap(0)
{
  // Constructor
  fVariableList = new TList();
  fVariableVertexZList = new TList();
  fVariableVertexZList->SetOwner();
}

AliMultInput::AliMultInput(const char * name, const char * title):
TNamed(name,title),
fNVars(0), fVariableList(0x0),
fNVtxZ(0), fVariableVertexZList(0x0),
fMap(0)
{
  // Constructor
  fVariableList = new TList();
  fVariableVertexZList = new TList();
  fVariableVertexZList->SetOwner();
}

AliMultInput::AliMultInput(const AliMultInput& o)
: TNamed(o),
fNVars(0), fVariableList(0x0),
fNVtxZ(0), fVariableVertexZList(0x0),
fMap(0)
{
  // Constructor
  fVariableList = new TList();
  TIter next(o.fVariableList);
  AliMultVariable* v  = 0;
  while ((v = static_cast<AliMultVariable*>(next())))  AddVariable(v);
  
  // Constructor
  fVariableVertexZList = new TList();
  fVariableVertexZList->SetOwner();
  TIter nextp(o.fVariableVertexZList);
  TProfile* vp  = 0;
  while ((vp = static_cast<TProfile*>(nextp())))  AddVtxZ(vp);
}

AliMultInput& AliMultInput::operator=(const AliMultInput& o)
{
  if (&o == this) return *this;
  SetName(o.GetName());
  SetTitle(o.GetTitle());
  if (!fVariableList) fVariableList = new TList();
  fVariableList->Clear();
  fNVars = 0;
  TIter next(o.fVariableList);
  AliMultVariable* v  = 0;
  while ((v = static_cast<AliMultVariable*>(next())))  AddVariable(v);
  
  if (!fVariableVertexZList) fVariableVertexZList = new TList();
  fVariableVertexZList->SetOwner();
  fVariableVertexZList->Clear();
  fNVtxZ = 0;
  TIter nextp(o.fVariableVertexZList);
  TProfile* vp  = 0;
  while ((vp = static_cast<TProfile*>(nextp())))  AddVtxZ(vp);
  
  return *this;
}

AliMultInput::~AliMultInput(){
  // destructor
  
}
void AliMultInput::AddVariable ( AliMultVariable *lVar )
{
  if (!lVar) return;
  if (!fVariableList) {
    fVariableList = new TList;
    fNVars = 0;
  }
  //Protect against double-declaration
  if ( fVariableList->FindObject( lVar->GetName() ) ){
    Printf("===========================================================================" );
    Printf("                          !!!  WARNING !!!                                 " );
    Printf("Variable named %s already exists, you're doing a double-declaration!",lVar->GetName() );
    Printf("AddVariable call exiting without doing anything... Please check your logic!");
    Printf("===========================================================================" );
    return;
  }
  fVariableList->Add(lVar);
  fNVars++;
}

AliMultVariable* AliMultInput::GetVariable (const TString& lName) const
{
  if (!fVariableList) return 0;
  return static_cast<AliMultVariable*>(fVariableList->FindObject(lName));
}

AliMultVariable* AliMultInput::GetVariable (Long_t iIdx) const
{
  if (!fVariableList) return 0;
  if (iIdx < 0 || iIdx >= fNVars) return 0;
  return static_cast<AliMultVariable*>(fVariableList->At(iIdx));
}

void AliMultInput::AddVtxZ ( TProfile *prof )
{
  //Warning: not protected against naming!
  if (!fVariableVertexZList) {
    fVariableVertexZList = new TList;
    fVariableVertexZList->SetOwner();
    fNVtxZ = 0;
  }
  //Protect against double-declaration
  if ( fVariableVertexZList->FindObject( prof->GetName() ) ){
    Printf("===========================================================================" );
    Printf("                          !!!  WARNING !!!                                 " );
    Printf("A vtx-Z correction named %s already exists!",prof->GetName() );
    Printf("AddVtxZ call exiting without doing anything... Please check your logic!");
    Printf("===========================================================================" );
    return;
  }
  TProfile *lProfClo = (TProfile*) prof->Clone(Form("copy_%s",prof->GetName()));
  fVariableVertexZList->Add(lProfClo);
  fNVtxZ++;
}

TProfile* AliMultInput::GetVtxZProfile (const AliMultVariable *v) const
{
  if (!fMap) return 0;
  TPair* ret = static_cast<TPair*>(fMap->FindObject(v));
  if (!ret) return 0;
  return static_cast<TProfile*>(ret->Value());
}

Double_t AliMultInput::GetVtxZCorrection (const AliMultVariable *v, Double_t lVtxZ) const
{
  Float_t lReturnValue = 1.0;
  if (!fMap) return 0;
  TPair* ret = static_cast<TPair*>(fMap->FindObject(v));
  if (!ret) return 0;
  TProfile *profile = static_cast<TProfile*>(ret->Value());
  Int_t lBin = profile -> FindBin(lVtxZ);
  Float_t lBinContent = profile->GetBinContent(lBin);
  Float_t lBinCenter = profile->GetBinCenter(lBin);
  Float_t lBinWidth = profile->GetBinWidth(lBin);
  if(lBinContent<1e-6) return 1.0; //not dealt with / uncalibrated
  if( lVtxZ >= lBinCenter ){
    //check next bin
    Float_t lFrac = (lVtxZ-lBinCenter)/lBinWidth;
    Float_t lBinContent1 = profile->GetBinContent(lBin+1);
    if(lBinContent1<1e-6) lBinContent1 = lBinContent;
    lReturnValue = (1-lFrac)*lBinContent + lFrac*lBinContent1;
  }else{
    //check preceding bin
    Float_t lFrac = (lBinCenter-lVtxZ)/lBinWidth;
    Float_t lBinContent1 = profile->GetBinContent(lBin-1);
    if(lBinContent1<1e-6) lBinContent1 = lBinContent;
    lReturnValue = (1-lFrac)*lBinContent + lFrac*lBinContent1;
  }
  return lReturnValue;
}

void AliMultInput::ClearVtxZ()
{
  if (fVariableVertexZList)
    fVariableVertexZList->Delete();
  fNVtxZ = fVariableVertexZList->GetEntries(); 
}

void AliMultInput::Clear(Option_t* option)
{
  TIter next(fVariableList);
  AliMultVariable* var = 0;
  while ((var = static_cast<AliMultVariable*>(next()))) {
    var->Clear(option);
  }
}

void AliMultInput::Set(const AliMultInput* other)
{
  TIter next(fVariableList);
  AliMultVariable* var = 0;
  while ((var = static_cast<AliMultVariable*>(next()))) {
    AliMultVariable* ovar = other->GetVariable(var->GetName());
    var->Set(ovar);
  }
}
void AliMultInput::Print(Option_t* option) const
{
  Printf("%s: %s/%s %ld variables", ClassName(),
         GetName(), GetTitle(), fNVars);
  gROOT->IndentLevel();
  Printf("Variables");
  gROOT->IncreaseDirLevel();
  TIter next(fVariableList);
  TObject* o = 0;
  while ((o = next())) {
    gROOT->IndentLevel();
    o->Print(option);
  }
  gROOT->DecreaseDirLevel();
}
//________________________________________________________________
void AliMultInput::SetupAutoVtxZCorrection(){
  Printf("Setting up AliMultVariable <-> TProfile map");
  
  if (fMap) {
      delete fMap;
      fMap = 0;
  }
  fMap = new TMap;
  fMap->SetOwner(false);
  Printf("Variable list size: %i, vertex-Z list size: %i",fVariableList->GetEntries(), fVariableVertexZList->GetEntries());
  
  for(Int_t ii=0 ; ii<kNVariables; ii++){
    //Get var
    AliMultVariable* v = static_cast<AliMultVariable*>(fVariableList->At(ii));
    if (!v) continue;
    
    v->SetUseVertexZCorrection(kFALSE); //assume not found until proven otherwise
    
    //Warning: assumes complete list, usable at calibration time only
    TProfile*   h = static_cast<TProfile*>(fVariableVertexZList->At(ii));
    if (!h) continue;
    
    //Printf("Found vertex-Z calib profile for variable %s, profile name: %s", v->GetName(), h->GetName() );
    v->SetUseVertexZCorrection(kTRUE); //assume not found until proven otherwise
    
    fMap->Add(v, h);
  }
  Printf("Number of automatically calibrated inputs: %i",fMap->GetEntries());
}
//________________________________________________________________
void AliMultInput::SetupAutoVtxZCorrection( const AliOADBMultSelection *oadb){
  Printf("Setting up automatic vertex-Z correction...");
  
  if (fMap) {
      delete fMap;
      fMap = 0;
  }
  fMap = new TMap;
  fMap->SetOwner(false);
  
  for(Int_t ii=0 ; ii<kNVariables; ii++){
    //Get var
    AliMultVariable* v = static_cast<AliMultVariable*>(fVariableList->At(ii));
    if (!v) continue;
    //Printf("Looking for vertex-Z correction for variable: %s", v->GetName() );
    v->SetUseVertexZCorrection(kFALSE); //assume not found until proven otherwise
    
    //Check existence of calib histo
    TString name(Form("hCalibVtx_%s", v->GetName()));
    TProfile*   h = oadb->GetCalibHistoVtx(name);
    if (!h) continue;
    
    Printf("Found vertex-Z calib profile for variable %s, profile name: %s", v->GetName(), h->GetName() );
    v->SetUseVertexZCorrection(kTRUE); //assume not found until proven otherwise
    
    fMap->Add(v, h);
  }
  Printf("Number of automatically calibrated inputs: %i",fMap->GetEntries()); 
}


