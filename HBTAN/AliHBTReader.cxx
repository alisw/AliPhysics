#include "AliHBTReader.h"

#include "AliHBTParticleCut.h"


ClassImp(AliHBTReader)
//pure virtual

/*************************************************************************************/

AliHBTReader::AliHBTReader()
{
//constructor
 fCuts = new TObjArray();
}

/*************************************************************************************/

AliHBTReader::~AliHBTReader()
{
//destructor
  fCuts->SetOwner();
  delete fCuts;
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
  
 for(Int_t i=0; i<fCuts->GetEntriesFast(); i++)   
   {
     AliHBTParticleCut &cut = *((AliHBTParticleCut*)fCuts->At(i));
     //if some of cuts accepts all particles or some accepts particles of this type, accept
     if ( (cut.GetPID() == 0) || (cut.GetPID() == pid) ) return kFALSE; 
   }
 return kTRUE;
}
/*************************************************************************************/
