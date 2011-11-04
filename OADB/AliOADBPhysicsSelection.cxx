#include "AliOADBPhysicsSelection.h"
#include "TList.h"
#include "TString.h"
#include "TObjString.h"
#include "AliLog.h"
#include "TBrowser.h"
#include "TFolder.h"
#include "TIterator.h"
#include <iostream>

using namespace std;

ClassImp(AliOADBPhysicsSelection)


AliOADBPhysicsSelection::AliOADBPhysicsSelection() : 
TNamed("AliOADBPhysicsSelection", "OADB object for the physics selection"), fNtriggerBits(NTRIGGERBITS), fNtriggerLogics(NTRIGGERLOGICS),fCollTrigClasses(0),fBGTrigClasses(0),fHardwareTrigger(0),fOfflineTrigger(0),fBeamSide(0)
{
  // default ctor
  
}
AliOADBPhysicsSelection::AliOADBPhysicsSelection(const char* name) : 
  TNamed(name, "OADB object for the physics selection"), fNtriggerBits(NTRIGGERBITS), fNtriggerLogics(NTRIGGERLOGICS),fCollTrigClasses(0),fBGTrigClasses(0),fHardwareTrigger(0),fOfflineTrigger(0),fBeamSide(0)
{
  // ctor, better use this one
  Init();
}

void AliOADBPhysicsSelection::Init() {
  //Initialization of pointers

  if(fCollTrigClasses) {
    AliInfo("Already initialized");
    return;
  }

  fNtriggerBits = NTRIGGERBITS;
  fCollTrigClasses = new TList*[fNtriggerBits];
  fBGTrigClasses   = new TList*[fNtriggerBits];
  fHardwareTrigger = new TObjString[fNtriggerLogics];
  fOfflineTrigger  = new TObjString[fNtriggerLogics];  
  fBeamSide = new TMap;
  for(UInt_t ibit = 0; ibit < fNtriggerBits; ibit++){
    fCollTrigClasses[ibit] = new TList;
    fBGTrigClasses  [ibit] = new TList;
  }
  for(UInt_t ilogic = 0; ilogic < fNtriggerLogics; ilogic++){
    fHardwareTrigger[ilogic] = "";
    fOfflineTrigger [ilogic] = "";
  }
  

}

AliOADBPhysicsSelection::~AliOADBPhysicsSelection(){
  // dtor

  if(fHardwareTrigger) delete [] fHardwareTrigger;
  if(fOfflineTrigger)  delete [] fOfflineTrigger ;
  
  for(Int_t ibit = 0; ibit < NTRIGGERBITS; ibit++){
    
    if(fCollTrigClasses)
      if(fCollTrigClasses[ibit]) delete fCollTrigClasses[ibit];
    if(fBGTrigClasses)
      if(fBGTrigClasses  [ibit]) delete fBGTrigClasses  [ibit];
  }

}

AliOADBPhysicsSelection::AliOADBPhysicsSelection(const AliOADBPhysicsSelection& cont) :TNamed("AliOADBPhysicsSelection", "OADB object for the physics selection"), fNtriggerBits(NTRIGGERBITS), fNtriggerLogics(NTRIGGERLOGICS),fCollTrigClasses(0),fBGTrigClasses(0),fHardwareTrigger(0),fOfflineTrigger(0),fBeamSide(0) {
  // Copy ctor
  fCollTrigClasses = cont.fCollTrigClasses;
  fBGTrigClasses = cont.fBGTrigClasses;
  fHardwareTrigger = cont.fHardwareTrigger;
  fOfflineTrigger = cont.fOfflineTrigger;
  fBeamSide = cont.fBeamSide;
  
}

AliOADBPhysicsSelection& AliOADBPhysicsSelection::operator=(const AliOADBPhysicsSelection& other)  {
  //Assignment operator
  TNamed::operator=(other);

  fCollTrigClasses = other.fCollTrigClasses;
  fBGTrigClasses = other.fBGTrigClasses;
  fHardwareTrigger = other.fHardwareTrigger;
  fOfflineTrigger = other.fOfflineTrigger;
  fBeamSide = other.fBeamSide;
  

  return *this;
}
 
void AliOADBPhysicsSelection::AddCollisionTriggerClass(AliVEvent::EOfflineTriggerTypes triggerMask, const char* className,const char * beamSide, UInt_t triggerLogic) {
  // add collision trigger class
  TObjString * tclass = new TObjString(Form("%s &%u *%u",className,triggerMask, triggerLogic));
  fCollTrigClasses[GetActiveBit(triggerMask)]->Add(tclass); 
  SetBeamSide(tclass->String().Data(),beamSide);
}
void AliOADBPhysicsSelection::AddBGTriggerClass       (AliVEvent::EOfflineTriggerTypes triggerMask, const char* className,const char * beamSide, UInt_t triggerLogic) 
{ 
  // add bg trigger class
  TObjString * tclass = new TObjString(Form("%s &%u *%u",className,triggerMask, triggerLogic));
  fBGTrigClasses  [GetActiveBit(triggerMask)]->Add(tclass);
  SetBeamSide(tclass->String().Data(),beamSide);
}

const TString AliOADBPhysicsSelection::GetBeamSide (const char * trigger)  {
  // Associate beam side to trigger class name
  TObjString * cname = new TObjString(trigger);
  CleanKey(cname->String());  
  static TString retValue="";
  retValue = ((TObjString*)fBeamSide->GetValue(cname->String().Data())) ?  
    ((TObjString*)fBeamSide->GetValue(cname->String().Data()))->String() : "" ;
  delete cname;
  return retValue;
}
void AliOADBPhysicsSelection::SetBeamSide (const char * className, const char * side)
{ 
  // return beam side
  TObjString * cname = new TObjString(className);
  CleanKey(cname->String());
  fBeamSide->Add(new TObjString(cname->String().Data()), new TObjString(side));
  delete cname;
}

void AliOADBPhysicsSelection::CleanKey(TString & str) {

  //  Remove all wite spaces and "goodies" of the trigger class string (bx ids, trigger logic...)
  if(str.Index("*")>0)
    str.Remove(str.Index("*")); // keep only the class name (no bx, offline trigger...)   
  if(str.Index("#")>0)
    str.Remove(str.Index("#")); // keep only the class name (no bx, offline trigger...)   
  if(str.Index("&")>0)
    str.Remove(str.Index("&")); // keep only the class name (no bx, offline trigger...)   
  str.ReplaceAll(" ","");

}
 
// Getters

// FIXME: decide about these getters
// TList * AliOADBPhysicsSelection::GetCollTrigClass(AliVEvent::EOfflineTriggerTypes trig) const { 
//   // Returns list of collision trigger classes for the requested bits
//   return GetTrigClass(trig,fCollTrigClasses);
// }

// TList * AliOADBPhysicsSelection::GetBGTrigClass(AliVEvent::EOfflineTriggerTypes trig) const { 
//   // Returns list of background trigger classes for the requested bits
//   return GetTrigClass(trig,fBGTrigClasses);
// }

// TList * AliOADBPhysicsSelection::GetTrigClass(AliVEvent::EOfflineTriggerTypes trig, TList ** listClasses) const { 
//   // Returns list of trigger classes for the requested bits

//   TList * list= new TList;  
//   // Check which bits are on and build the final list
//   for(UInt_t ibit = 0; ibit < fNtriggerBits; ibit++){
//     if ((trig&(0x1<<ibit))) {
//       if(listClasses[ibit]) list->AddAll(listClasses[ibit]);
//       else AliError(Form("List %d not initialized?",ibit));
//     }
//   }
  
//   return list;
// }

void AliOADBPhysicsSelection::Print(Option_t* option) const {

  for(UInt_t ibit = 0; ibit < fNtriggerBits; ibit++){
    cout << "Bit " << ibit << endl;
    
    fCollTrigClasses[ibit]->Print(option);
    fBGTrigClasses[ibit]->Print(option);
    cout << "HW trig " << fHardwareTrigger[ibit].String().Data() << endl;
    cout << "Offline trig " << fHardwareTrigger[ibit].String().Data() << endl;
    
  }
  cout << "Beams: " << endl;
  fBeamSide->Print();

}


void AliOADBPhysicsSelection::Browse(TBrowser *b)
{
   // Browse this object.
   // If b=0, there is no Browse call TObject::Browse(0) instead.
   //         This means TObject::Inspect() will be invoked indirectly
  static TFolder  ** bitFolders = 0;
  if(!bitFolders) {
    bitFolders = new TFolder*[fNtriggerBits];
    for(UInt_t ibit = 0; ibit < fNtriggerBits; ibit++){
      bitFolders[ibit] = 0;
    }
    
  }

  if (b) {
    for(UInt_t ibit = 0; ibit < fNtriggerBits; ibit++){
      if(bitFolders[ibit]) delete bitFolders[ibit];
      bitFolders[ibit] = new TFolder (Form("Bit %2.2d", ibit), "Trigger bit folder");
      bitFolders[ibit]->SetOwner(); // Delete also the TObjString when deleting the folder
      for(UInt_t ilogic = 0; ilogic < fNtriggerLogics; ilogic++){
	if(GetHardwareTrigger(ilogic) != "" || 	GetOfflineTrigger(ilogic) != "")  {
	  bitFolders[ibit]->Add(new TObjString(Form("Hardware Trig    [*%d][%s]",ilogic,GetHardwareTrigger(ilogic).Data())));
	  bitFolders[ibit]->Add(new TObjString(Form("Offline  Trig    [*%d][%s]",ilogic,GetOfflineTrigger(ilogic).Data())));      
	}
      }
      
      TIterator *itColl = fCollTrigClasses[ibit]->MakeIterator();
      TObjString * coll = 0;

      while ((coll = (TObjString*)itColl->Next())) {
	bitFolders[ibit]->Add(new TObjString(Form("Collision Class  [%s] [%s]", coll->String().Data(), 
						  GetBeamSide(coll->String().Data()).Data())));
      }
      delete itColl;

      TIterator *itBG = fBGTrigClasses[ibit]->MakeIterator();
      TObjString * bg = 0;
      while ((bg = (TObjString*)itBG->Next())) {
	bitFolders[ibit]->Add(new TObjString(Form("Background Class [%s] [%s]", bg->String().Data(), 
						  GetBeamSide(bg->String().Data()).Data())));
      }
      delete itBG;

      b->Add(bitFolders[ibit]);

    }    
  }     
   else
      TObject::Browse(b);
}

UInt_t AliOADBPhysicsSelection::GetActiveBit(UInt_t mask) {
  // Returns the active bit index in the mask
  // Assumes only one bit is on.
  // If more than one bit is lit, prints an error and returns the first.
  // If no bit is on, prints an error and returns 0

  const Int_t kNBitsToCheck = 28;

  //  Int_t nbit = sizeof(mask)*8;
  Int_t nbit = kNBitsToCheck;
  Int_t activeBit = -1;
  for(Int_t ibit = 0; ibit < nbit; ibit++){
    if ( mask & (0x1 << ibit) ) {
      if (activeBit == -1) activeBit = ibit;
      else Printf("ERROR (AliTriggerAnalysis::GetActiveBit): More than one bit is on in this mask 0x%x", mask);
    }
  }
  if (activeBit == -1) {
    Printf("ERROR (AliTriggerAnalysis::GetActiveBit): No bit is on");
    activeBit=0;
  }

  return activeBit;

}
