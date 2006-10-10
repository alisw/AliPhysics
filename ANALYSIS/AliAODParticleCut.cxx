#include "AliAODParticleCut.h"
//__________________________________________________________________________
////////////////////////////////////////////////////////////////////////////
//                                                                        //
// class AliAODParticleCut                                                //
//                                                                        //
// Classes for single particle cuts.                                      //
// User should use mainly AliAODParticleCut interface methods,            //
// eventually EmptyCut which passes all particles.                        //
//                                                                        //
// There is all interface for setting cuts on all particle properties     //
// The main method is Rejected - which returns                            //
//         True to reject particle                                        //
//         False in case it meets all the criteria of the given cut       //
//                                                                        //
// This class has the list of base particle  cuts that perform check on   //
// single property. Particle  is rejected if any of cuts rejects it.      //
// There are implemented logical base cuts that perform logical           //
// operations on results of two other base cuts. Using them user can      //
// create a tree structure of a base cuts that performs sophisticated     //
// cut.                                                                   //
//                                                                        //
// User can also implement a base cut that performs complicated           //
// calculations, if it is only more convenient and/or efficint.           //
//                                                                        //
// User should delete created cuts  himself                               //
// because when setting a cut, other objects (functions,analyses,         //
// readers, other cuts) make their own copy of a cut.                     //
//                                                                        //
// more info: http://aliweb.cern.ch/people/skowron/analyzer/index.html    //
// responsible: Piotr Skowronski@cern.ch                                  //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>


ClassImp(AliAODParticleCut)
const Int_t AliAODParticleCut::fgkMaxCuts = 50;
/******************************************************************/

AliAODParticleCut::AliAODParticleCut():
 fCuts(new AliAODParticleBaseCut* [fgkMaxCuts]),//last property in the property enum => defines number of properties
 fNCuts(0),
 fPID(0)
{
  //default ctor
}
/******************************************************************/

AliAODParticleCut::AliAODParticleCut(const AliAODParticleCut& in):
  TObject(in),
  fCuts(new AliAODParticleBaseCut* [fgkMaxCuts]),//last property in the property
                                                 //property enum => defines number of properties
  fNCuts(in.fNCuts),
  fPID(in.fPID)
{
  //cpy ctor
  for (Int_t i = 0;i<fNCuts;i++)
   {
     fCuts[i] = (AliAODParticleBaseCut*)in.fCuts[i]->Clone();//create new object (clone) and rember pointer to it
   }
}
/******************************************************************/
AliAODParticleCut& AliAODParticleCut::operator=(const AliAODParticleCut& in)
{
  //assigment operator
  Info("operator=","operator=operator=operator=operator=\noperator=operator=operator=operator=");
  for (Int_t i = 0;i<fNCuts;i++)
   {
     delete fCuts[i];
   }

  fNCuts = in.fNCuts;
  fPID  = in.fPID;
  for (Int_t i = 0;i<fNCuts;i++)
   {
     fCuts[i] = (AliAODParticleBaseCut*)in.fCuts[i]->Clone();//create new object (clone) and rember pointer to it
   }
  return *this;
}    

/******************************************************************/
AliAODParticleCut::~AliAODParticleCut()
{
  //dtor
  for (Int_t i = 0;i<fNCuts;i++)
   {
     delete fCuts[i];
   }
  delete []fCuts;
} 
/******************************************************************/

Bool_t AliAODParticleCut::Rejected(AliVAODParticle* p) const
{
//method checks all the cuts that are set (in the list)
//If any of the baseCuts rejects particle False(rejection) is returned

 if(!p) 
  {
    Warning("Rejected()","No Pasaran! We never accept NULL pointers");
    return kTRUE;
  }
 if( (p->GetPdgCode() != fPID) && ( fPID != 0)) return kTRUE;
 
 for (Int_t i = 0;i<fNCuts;i++)
   {
    if ( (fCuts[i]->Rejected(p)) )
     {
//       fCuts[i]->Print();
       return kTRUE; //if one of the cuts rejects, then reject
     }
   }
  return kFALSE;
}
/******************************************************************/

void AliAODParticleCut::AddBasePartCut(AliAODParticleBaseCut* basecut)
{
  //adds the base pair cut (cut on one value)
 
   if (!basecut) return;
   if( fNCuts == (fgkMaxCuts-1) )
    {
      Warning("AddBasePartCut","Not enough place for another cut");
      return;
    }
   fCuts[fNCuts++]=basecut;
 
}

/******************************************************************/
AliAODParticleBaseCut* AliAODParticleCut::FindCut(AliAODParticleBaseCut::EAODCutProperty property)
{
 //returns pointer to the cut checking the given property
 for (Int_t i = 0;i<fNCuts;i++)
  {
    if (fCuts[i]->GetProperty() == property) 
       return fCuts[i]; //we found the cut we were searching for
  }
 
 return 0x0; //we did not found this cut
 
}
/******************************************************************/

void AliAODParticleCut::SetMomentumRange(Double_t min, Double_t max)
{
  //Sets momentum range
  AliAODMomentumCut* cut= (AliAODMomentumCut*)FindCut(AliAODParticleBaseCut::kAODP);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODMomentumCut(min,max);
}
/******************************************************************/


void AliAODParticleCut::SetPtRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODPtCut* cut= (AliAODPtCut*)FindCut(AliAODParticleBaseCut::kAODPt);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODPtCut(min,max);

}
/******************************************************************/

void AliAODParticleCut::SetEnergyRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODEnergyCut* cut= (AliAODEnergyCut*)FindCut(AliAODParticleBaseCut::kAODE);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODEnergyCut(min,max);
 
}
/******************************************************************/

void AliAODParticleCut::SetRapidityRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(AliAODParticleBaseCut::kAODRapidity);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODRapidityCut(min,max);

}
/******************************************************************/

void AliAODParticleCut::SetPseudoRapidityRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(AliAODParticleBaseCut::kAODPseudoRapidity);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODPseudoRapidityCut(min,max);
 
}
/******************************************************************/

void AliAODParticleCut::SetPxRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(AliAODParticleBaseCut::kAODPx);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODPxCut(min,max);
}
/******************************************************************/

void AliAODParticleCut::SetPyRange(Double_t min, Double_t max)
{  
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(AliAODParticleBaseCut::kAODPy);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODPyCut(min,max);
}
/******************************************************************/

void AliAODParticleCut::SetPzRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(AliAODParticleBaseCut::kAODPz);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODPzCut(min,max);
}
/******************************************************************/

void AliAODParticleCut::SetPhiRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(AliAODParticleBaseCut::kAODPhi);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODPhiCut(min,max);
}
/******************************************************************/

void AliAODParticleCut::SetThetaRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(AliAODParticleBaseCut::kAODTheta);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODThetaCut(min,max);
}
/******************************************************************/

void AliAODParticleCut::SetVxRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(AliAODParticleBaseCut::kAODVx);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODVxCut(min,max);
}
/******************************************************************/

void AliAODParticleCut::SetVyRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(AliAODParticleBaseCut::kAODVy);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODVyCut(min,max);
}
/******************************************************************/

void AliAODParticleCut::SetVzRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(AliAODParticleBaseCut::kAODVz);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODVzCut(min,max);
}

/******************************************************************/
void AliAODParticleCut::Streamer(TBuffer &b)
{
  // Stream all objects in the array to or from the I/O buffer.

   UInt_t R__s, R__c;
   if (b.IsReading()) 
    {
      Int_t i;
      for (i = 0;i<fNCuts;i++) delete fCuts[i];
      b.ReadVersion(&R__s, &R__c);
      TObject::Streamer(b);
      b >> fPID;
      b >> fNCuts;
      for (i = 0;i<fNCuts;i++) 
       {
         b >> fCuts[i];
       }
       b.CheckByteCount(R__s, R__c,AliAODParticleCut::IsA());
    } 
   else 
    {
     R__c = b.WriteVersion(AliAODParticleCut::IsA(), kTRUE);
     TObject::Streamer(b);
     b << fPID;
     b << fNCuts;
     for (Int_t i = 0;i<fNCuts;i++)
      {
       b << fCuts[i];
      }
     b.SetByteCount(R__c, kTRUE);
   }
}
/******************************************************************/

void AliAODParticleCut::Print(const Option_t * /*opt*/) const
{
  //prints all information about the cut to stdout
  cout<<"Printing AliAODParticleCut, this = "<<this<<endl;
  cout<<"fPID  "<<fPID<<endl;
  cout<<"fNCuts  "<<fNCuts <<endl;
  for (Int_t i = 0;i<fNCuts;i++)
      {
       cout<<"  fCuts["<<i<<"]  "<<fCuts[i]<<endl<<"   ";
       fCuts[i]->Print();
      }
}

/******************************************************************/
/******************************************************************/
ClassImp(AliAODParticleEmptyCut)

void AliAODParticleEmptyCut::Streamer(TBuffer &b)
 {
  //stramer
  AliAODParticleCut::Streamer(b);
 }
/******************************************************************/
/******************************************************************/
/******************************************************************/
