#include "AliAODParticleCut.h"
//__________________________________________________________________________
////////////////////////////////////////////////////////////////////////////
//                                                                        //
// class AliAODParticleCut                                                //
//                                                                        //
// Classes for single particle cuts                                       //
// User should use only AliAODParticleCut, eventually                     //
// EmptyCut which passes all particles                                    //
// There is all interface for setting cuts on all particle properties     //
// The main method is Pass - which returns                                //
//         True to reject particle                                        //
//         False in case it meets all the criteria of the given cut       //
//                                                                        //
// User should create (and also destroy) cuts himself                     // 
// and then pass them to the Analysis And Function by a proper method     //
//                                                                        //
//                                                                        //
// more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html   //
// responsible: Piotr Skowronski@cern.ch                                   //
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
 TObject(in)
{
  //cpy ctor
  fCuts = new AliAODParticleBaseCut* [fgkMaxCuts];//last property in the property
                                         //property enum => defines number of properties
  fNCuts = in.fNCuts;
  fPID  = in.fPID;
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

Bool_t AliAODParticleCut::Pass(AliVAODParticle* p) const
{
//method checks all the cuts that are set (in the list)
//If any of the baseCuts rejects particle False(rejection) is returned

 if(!p) 
  {
    Warning("Pass()","No Pasaran! We never accept NULL pointers");
    return kTRUE;
  }
 if( (p->GetPdgCode() != fPID) && ( fPID != 0)) return kTRUE;
 
 for (Int_t i = 0;i<fNCuts;i++)
   {
    if ( (fCuts[i]->Pass(p)) )
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
AliAODParticleBaseCut* AliAODParticleCut::FindCut(EAODCutProperty property)
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
  AliAODMomentumCut* cut= (AliAODMomentumCut*)FindCut(kAODP);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODMomentumCut(min,max);
}
/******************************************************************/


void AliAODParticleCut::SetPtRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODPtCut* cut= (AliAODPtCut*)FindCut(kAODPt);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODPtCut(min,max);

}
/******************************************************************/

void AliAODParticleCut::SetEnergyRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODEnergyCut* cut= (AliAODEnergyCut*)FindCut(kAODE);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODEnergyCut(min,max);
 
}
/******************************************************************/

void AliAODParticleCut::SetRapidityRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(kAODRapidity);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODRapidityCut(min,max);

}
/******************************************************************/

void AliAODParticleCut::SetPseudoRapidityRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(kAODPseudoRapidity);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODPseudoRapidityCut(min,max);
 
}
/******************************************************************/

void AliAODParticleCut::SetPxRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(kAODPx);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODPxCut(min,max);
}
/******************************************************************/

void AliAODParticleCut::SetPyRange(Double_t min, Double_t max)
{  
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(kAODPy);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODPyCut(min,max);
}
/******************************************************************/

void AliAODParticleCut::SetPzRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(kAODPz);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODPzCut(min,max);
}
/******************************************************************/

void AliAODParticleCut::SetPhiRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(kAODPhi);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODPhiCut(min,max);
}
/******************************************************************/

void AliAODParticleCut::SetThetaRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(kAODTheta);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODThetaCut(min,max);
}
/******************************************************************/

void AliAODParticleCut::SetVxRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(kAODVx);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODVxCut(min,max);
}
/******************************************************************/

void AliAODParticleCut::SetVyRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(kAODVy);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODVyCut(min,max);
}
/******************************************************************/

void AliAODParticleCut::SetVzRange(Double_t min, Double_t max)
{
  //name self descriptive
  AliAODParticleBaseCut* cut = FindCut(kAODVz);
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

void AliAODParticleCut::Print(void) const
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

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp(AliAODParticleBaseCut)
void AliAODParticleBaseCut::Print(void) const
{
  // prints the information anout the base cut to stdout
  cout<<"fMin="<<fMin <<", fMax=" <<fMax<<"    ";
  PrintProperty();
}
/******************************************************************/

void AliAODParticleBaseCut::PrintProperty(void) const
{
 //prints the property name 
 switch (fProperty)
  {
   case  kAODP: 
     cout<<"kAODP"; break;
   case  kAODPt: 
     cout<<"kAODPt"; break;
   case  kAODE: 
     cout<<"kAODE"; break;
   case  kAODRapidity: 
     cout<<"kAODRapidity"; break;
   case  kAODPseudoRapidity: 
     cout<<"kAODPseudoRapidity"; break;
   case  kAODPx: 
     cout<<"kAODPx"; break;
   case  kAODPy: 
     cout<<"kAODPy"; break;
   case  kAODPz: 
     cout<<"kAODPz"; break;   
   case  kAODPhi: 
     cout<<"kAODPhi"; break;
   case  kAODTheta: 
     cout<<"kAODTheta"; break;
   case  kAODVx: 
     cout<<"kAODVx"; break;
   case  kAODVy: 
     cout<<"kAODVy"; break;
   case  kAODVz: 
     cout<<"kAODVz"; break;
   case  kAODPid: 
     cout<<"kAODPid"; break;
   case  kAODNone: 
     cout<<"kAODNone"; break;
   default: 
     cout<<"Property Not Found";
  }
 cout<<endl;
}
ClassImp( AliAODMomentumCut )

ClassImp( AliAODPtCut )
ClassImp( AliAODEnergyCut )
ClassImp( AliAODRapidityCut )
ClassImp( AliAODPseudoRapidityCut )
ClassImp( AliAODPxCut )
ClassImp( AliAODPyCut )
ClassImp( AliAODPzCut )
ClassImp( AliAODPhiCut )
ClassImp( AliAODThetaCut )
ClassImp( AliAODVxCut )
ClassImp( AliAODVyCut )
ClassImp( AliAODVzCut )

ClassImp( AliAODPIDCut )

void AliAODPIDCut::Print(void) const
{
  cout<<"PID "<<fPID<<" ";
  AliAODParticleBaseCut::Print();
}

ClassImp( AliAODLogicalOperCut )

AliAODLogicalOperCut::AliAODLogicalOperCut():
 AliAODParticleBaseCut(-10e10,10e10,kAODNone),
 fFirst(new AliAODDummyBaseCut),
 fSecond(new AliAODDummyBaseCut)
{
 //ctor
}
/******************************************************************/

AliAODLogicalOperCut::AliAODLogicalOperCut(AliAODParticleBaseCut* first, AliAODParticleBaseCut* second):
 AliAODParticleBaseCut(-10e10,10e10,kAODNone),
 fFirst((first)?(AliAODParticleBaseCut*)first->Clone():0x0),
 fSecond((second)?(AliAODParticleBaseCut*)second->Clone():0x0)
{
  //ctor
  if ( (fFirst && fSecond) == kFALSE) 
   {
     Fatal("AliAODLogicalOperCut","One of parameters is NULL!");
   }
}
/******************************************************************/

AliAODLogicalOperCut::~AliAODLogicalOperCut()
{
  //destructor
  delete fFirst;
  delete fSecond;
}
/******************************************************************/

Bool_t AliAODLogicalOperCut::AliAODDummyBaseCut::Pass(AliVAODParticle* /*part*/)  const
{
  //checks if particles passes properties defined by this cut
  Warning("Pass","You are using dummy base cut! Probobly some logical cut is not set up properly");
  return kFALSE;//accept
}
/******************************************************************/

void AliAODLogicalOperCut::Streamer(TBuffer &b)
{
  // Stream all objects in the array to or from the I/O buffer.
  UInt_t R__s, R__c;
  if (b.IsReading()) 
   {
     delete fFirst;
     delete fSecond;
     fFirst  = 0x0;
     fSecond = 0x0;

     b.ReadVersion(&R__s, &R__c);
     TObject::Streamer(b);
     b >> fFirst;
     b >> fSecond;
     b.CheckByteCount(R__s, R__c,AliAODLogicalOperCut::IsA());
   } 
  else 
   {
     R__c = b.WriteVersion(AliAODLogicalOperCut::IsA(), kTRUE);
     TObject::Streamer(b);
     b << fFirst;
     b << fSecond;
     b.SetByteCount(R__c, kTRUE);
  }
}

/******************************************************************/
ClassImp(AliAODOrCut)

Bool_t AliAODOrCut::Pass(AliVAODParticle * p) const
{
  //returns true when rejected 
  //AND operation is a little bit misleading but is correct
  //User wants to build logical cuts with natural (positive) logic
  //while AODAN use inernally reverse (returns true when rejected)
  if (fFirst->Pass(p) && fSecond->Pass(p)) return kTRUE;//rejected (both rejected, returned kTRUE)
  return kFALSE;//accepted, at least one accepted (returned kFALSE)
}
/******************************************************************/

ClassImp(AliAODAndCut)

Bool_t AliAODAndCut::Pass(AliVAODParticle * p)  const
{
  //returns true when rejected 
  //OR operation is a little bit misleading but is correct
  //User wants to build logical cuts with natural (positive) logic
  //while AODAN use inernally reverse (returns true when rejected)
  if (fFirst->Pass(p) || fSecond->Pass(p)) return kTRUE;//rejected (any of two rejected(returned kTRUE) )
  return kFALSE;//accepted (both accepted (returned kFALSE))
}
/******************************************************************/
