#include "AliHBTReader.h"

#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TClass.h>
#include <Riostream.h>

#include "AliHBTParticleCut.h"


ClassImp(AliHBTReader)
//pure virtual

/*************************************************************************************/

AliHBTReader::AliHBTReader()
{
//constructor
 fCuts = new TObjArray();
 fDirs = 0x0;
}

/*************************************************************************************/
AliHBTReader::AliHBTReader(TObjArray* dirs)
 {
  fCuts = new TObjArray();
  fDirs = dirs;
 }

AliHBTReader::~AliHBTReader()
{
//destructor
 if(fCuts)
  {
   fCuts->SetOwner();
   delete fCuts;
  }
}

/*************************************************************************************/

void AliHBTReader::AddParticleCut(AliHBTParticleCut* cut)
{
 //sets the new cut 
 
  if (!cut) //if cut is NULL return with error
   {
    Error("AddParticleType","NULL pointers are not accepted any more.\nIf You want to accept all particles of this type, set an empty cut ");
    return;
   }
  AliHBTParticleCut *c = (AliHBTParticleCut*)cut->Clone();
  fCuts->Add(c);
}

/*************************************************************************************/

Bool_t AliHBTReader::Pass(AliHBTParticle* p)
 {
 //Method examines whether particle meets all cut and particle type criteria
  
   if(p==0x0)//of corse we not pass NULL pointers
    {
     Warning("Pass()","No Pasaran! We never accept NULL pointers");
     return kTRUE;
    }
   //if no particle is specified, we pass all particles
   //excluding NULL pointers, of course
  if ( fCuts->GetEntriesFast() == 0 ) return kFALSE; //if no cut specified accept all particles
  for(Int_t i=0; i<fCuts->GetEntriesFast(); i++)   
   {
     AliHBTParticleCut &cut = *((AliHBTParticleCut*)fCuts->At(i));
     if(!cut.Pass(p)) return kFALSE;  //accepted
   }
   
  return kTRUE;//not accepted

 }
/*************************************************************************************/

Bool_t  AliHBTReader::Pass(Int_t pid)
{
//this method checks if any of existing cuts accepts this pid particles
//or any cuts accepts all particles

 if(pid == 0)
  return kTRUE;

 if ( fCuts->GetEntriesFast() == 0 ) return kFALSE; //if no cut specified accept all particles
  
 for(Int_t i=0; i<fCuts->GetEntriesFast(); i++)   
   {
     AliHBTParticleCut &cut = *((AliHBTParticleCut*)fCuts->At(i));
     //if some of cuts accepts all particles or some accepts particles of this type, accept
     if ( (cut.GetPID() == 0) || (cut.GetPID() == pid) ) return kFALSE; 
   }
 return kTRUE;
}
/*************************************************************************************/

TString& AliHBTReader::GetDirName(Int_t entry)
 {
   TString* retval;//return value
   if (fDirs ==  0x0)
    {
      retval = new TString(".");
      return *retval;
    }
   
   if ( (entry>fDirs->GetEntries()) || (entry<0))//if out of bounds return empty string
    {                                            //note that entry==0 is accepted even if array is empty (size=0)
      Error("GetDirName","Name out of bounds");
      retval = new TString();
      return *retval;
    }
   
   if (fDirs->GetEntries() == 0)
    { 
      retval = new TString(".");
      return *retval;
    }
   
   TClass *objclass = fDirs->At(entry)->IsA();
   TClass *stringclass = TObjString::Class();
   
   TObjString *dir = (TObjString*)objclass->DynamicCast(stringclass,fDirs->At(entry));
   
   if(dir == 0x0)
    {
      Error("GetDirName","Object in TObjArray is not a TObjString or its descendant");
      retval = new TString();
      return *retval;
    }
   if (gDebug > 0) cout<<"Returned ok "<<dir->String().Data()<<endl;
   return dir->String();
 }

