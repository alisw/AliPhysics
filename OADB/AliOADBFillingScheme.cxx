//-------------------------------------------------------------------------
//     OADB container for filling scheme information (BX ids, name ...)
//     Author: Michele Floris, CERN
//-------------------------------------------------------------------------

#include "AliOADBFillingScheme.h"
#include "TMap.h"
#include "AliLog.h"
#include "TBrowser.h"
#include "TFolder.h"
#include <iostream>

using namespace std;

ClassImp(AliOADBFillingScheme)


AliOADBFillingScheme::AliOADBFillingScheme() : TNamed("AliOADBFillingScheme", "OADB object storing filling scheme infos"), fFSName(""), fBXIds(0){
  // default ctor

  
}
AliOADBFillingScheme::AliOADBFillingScheme(char* name) : TNamed(name, "OADB object storing filling scheme infos"), fFSName(""), fBXIds(0){
  // ctor
  Init();
}

void AliOADBFillingScheme::Init() {
  // initialize pointers
  fBXIds = new TMap();
  fFSName = "";

}

AliOADBFillingScheme::~AliOADBFillingScheme(){
  // dtor

  if(fBXIds)           delete fBXIds;

}

// AliOADBFillingScheme::AliOADBFillingScheme(const AliOADBFillingScheme& cont) {
//   // Copy ctor
//   AliError("To be implemented");
// }

// AliOADBFillingScheme& AliOADBFillingScheme::operator=(const AliOADBFillingScheme& cont) {
//   //Assignment operator
//   AliError("To be implemented");
// }
  
// Getters

const char *  AliOADBFillingScheme::GetBXIDs(const char * beamSide) const 
{
  //  Returns the bunch crossing numbers for the different beam classes. By default this is empty

  if (!strcmp(beamSide, "AC") && !(TObjString*)fBXIds->GetValue(beamSide)) {

    TString  &bxa =  ((TObjString*)fBXIds->GetValue("A"))->String(); 
    TString  &bxc =  ((TObjString*)fBXIds->GetValue("C"))->String();
    if(bxa.IsNull() && bxc.IsNull()) return "";
    if(bxc.IsNull())         return bxa.Data();
    if(bxa.IsNull())         return bxc.Data();
    static TString bxBoth = bxa.Data();
    bxBoth += bxc.Data();
    return bxBoth.Data();

  } 

  if(!(TObjString*)fBXIds->GetValue(beamSide)) return "";

  TString  &bx =  ((TObjString*)fBXIds->GetValue(beamSide))->String(); 
  if(bx.IsNull()) return "";
  return bx.Data();
  
}

void AliOADBFillingScheme::Browse(TBrowser *b)
{
   // Browse this object.
   // If b=0, there is no Browse call TObject::Browse(0) instead.
   //         This means TObject::Inspect() will be invoked indirectly


  if (b) {
    // Creates a folder for each beam type containing the list of corresponding bx ids
    b->Add(new TObjString(Form("Filling Scheme %s",GetFillingSchemeName())));
    TIterator * mapIter = fBXIds->MakeIterator();
    TObjString * key = 0;
    while ((key = (TObjString*) mapIter->Next())) {
      TFolder * folder = new TFolder(key->String().Data(), "beam side");
      TObjString * value = (TObjString*) fBXIds->GetValue(key);
      TObjArray * tokens = value->String().Tokenize(" ");
      TIterator * tokIter = tokens->MakeIterator();
      TObjString * bxNum = 0;
      while ((bxNum = (TObjString*) tokIter->Next())){
	folder->Add(bxNum);
      }
      b->Add(folder);
      delete tokIter;
    }
    delete mapIter;    
  }     
  else
    TObject::Browse(b);
}

void AliOADBFillingScheme::Print(Option_t* option) const {
  // Print Class contents
  // Option is passed to TMap::Print
  cout << "Filling scheme Name " <<  fFSName.Data() << option << endl;
  fBXIds->Print(option);
}
