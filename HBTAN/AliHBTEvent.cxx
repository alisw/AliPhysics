//__________________________________________________________
///////////////////////////////////////////////////////////////////
//
// class AliHBTEvent
//
// This class is container for paticles coming from one event
//
// Piotr.Skowronski@cern.ch 
//
///////////////////////////////////////////////////////////////////


#include "AliHBTEvent.h"
#include "AliHBTParticle.h"

ClassImp(AliHBTEvent)

const UInt_t AliHBTEvent::fgkInitEventSize = 100;

/**************************************************************************/ 
 
AliHBTEvent::AliHBTEvent():
 fSize(fgkInitEventSize),
 fParticles(new AliHBTParticle* [fSize]),
 fNParticles(0),
 fOwner(kTRUE),
 fRandomized(kFALSE)
 {
//default constructor   
  if(fgkInitEventSize<1) 
   {
    Fatal("AliHBTEvent::AliHBTEvent()",
          "fgkInitEventSize has a stiupid value (%d). Change it to positive number and recompile",
           fgkInitEventSize);

   }
 }
/**************************************************************************/ 
 
AliHBTEvent::AliHBTEvent(const AliHBTEvent& source):
 TObject(source),
 fSize(source.fSize),
 fParticles(new AliHBTParticle* [fSize]),
 fNParticles(source.fNParticles),
 fOwner(source.fNParticles),
 fRandomized(source.fRandomized)
{
//copy constructor
  for(Int_t i =0; i<fNParticles; i++)
   {
      fParticles[i] = new AliHBTParticle( *(source.fParticles[i]) );
   }
  
}
/**************************************************************************/ 

AliHBTEvent& AliHBTEvent::operator=(const AliHBTEvent source)
{
  // assigment operator
  Reset();
  if (fParticles) delete [] fParticles;
  fSize = source.fSize;
  fParticles = new AliHBTParticle* [fSize];
  fNParticles = source.fNParticles;
  fOwner = source.fNParticles;
  fRandomized = source.fRandomized;
      
  for(Int_t i =0; i<fNParticles; i++)
   {
      fParticles[i] = new AliHBTParticle( *(source.fParticles[i]) );
   }
  return *this;
}
/**************************************************************************/ 
AliHBTEvent::~AliHBTEvent()
 {
//destructor   
  this->Reset();//delete all particles
  if(fParticles)
   { 
    delete [] fParticles; //and delete array itself
   }
  fParticles = 0x0;
 }
/**************************************************************************/ 
void  AliHBTEvent::Reset()
{
  //deletes all particles from the event
  if(fParticles && fOwner)
    {
      for(Int_t i =0; i<fNParticles; i++)
       {
         for (Int_t j = i+1; j<fNParticles; j++)
           if (fParticles[j] == fParticles[i]) fParticles[j] = 0x0;
         delete fParticles[i];
       }
    }
   fNParticles = 0;
} 
/**************************************************************************/ 

AliHBTParticle* AliHBTEvent::GetParticleSafely(Int_t n)
{
  //returns nth particle with range check
  if( (n<0) || (fNParticles<=n) ) return 0x0;
  else return fParticles[n];
}
/**************************************************************************/ 

void  AliHBTEvent:: AddParticle(AliHBTParticle* hbtpart)
{
  //Adds new perticle to the event
  if ( fNParticles+1 >= fSize) Expand(); //if there is no space in array, expand it
  fParticles[fNParticles++] = hbtpart; //add a pointer
}
/**************************************************************************/ 

void  AliHBTEvent::AddParticle(TParticle* part, Int_t idx)
{
  //Adds TParticle to event
  AddParticle( new AliHBTParticle(*part,idx) );
}
/**************************************************************************/ 
void  AliHBTEvent::AddParticle(Int_t pdg, Int_t idx, 
             Double_t px, Double_t py, Double_t pz, Double_t etot,
             Double_t vx, Double_t vy, Double_t vz, Double_t time)
{
  //adds particle to event
  AddParticle(new  AliHBTParticle(pdg,idx,px,py,pz,etot,vx,vy,vz,time) );
}
/**************************************************************************/ 

void AliHBTEvent::Expand()
{
//expands the array with pointers to particles
//about the size defined in fgkInitEventSize

 fSize+=fgkInitEventSize;
 AliHBTParticle** tmpParticles = new AliHBTParticle* [fSize]; //create new array of pointers
 //check if we got memory / if not Abort
 if (!tmpParticles) Fatal("AliHBTEvent::Expand()","No more space in memory");


 for(Int_t i = 0; i<fNParticles ;i++)
  //copy pointers to the new array
  {
    tmpParticles[i] = fParticles[i];
  }
 delete [] fParticles; //delete old array
  fParticles = tmpParticles; //copy new pointer to the array of pointers to particles
}
