#include "AliAODParticleBaseCut.h"
//__________________________________________________________________________
////////////////////////////////////////////////////////////////////////////
//                                                                        //
// class AliAODParticleBaseCut                                            //
//                                                                        //
// Set of classes for performing cuts on particle properties of           //
// AliAODParticleBaseCut is a base class for "base                        //
// particle cuts". Further, there are implemented classes that performs   //
// cuts on the most common particle properties like pt, pseudo rapidity,  //
// angles, anergy, etc.                                                   //
//                                                                        //
// There are also implemeted base cuts that perform logical operations    //
// on results of base particle cuts: AliAODOrCut and  AliAODAndCut.       //
//                                                                        //
// Each base cut has a property, thet allows to distinguish them.         //
// This functionality is used by the interface methods of Particle Cut    //
// that allows easy update ranges.                                        //
//                                                                        //
// more info: http://aliweb.cern.ch/people/skowron/analyzer/index.html    //
// responsible: Piotr Skowronski@cern.ch                                  //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include <Riostream.h>

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

Bool_t AliAODLogicalOperCut::AliAODDummyBaseCut::Rejected(AliVAODParticle* /*part*/)  const
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

Bool_t AliAODOrCut::Rejected(AliVAODParticle * p) const
{
  //returns true when rejected 
  //AND operation is a little bit misleading but is correct
  //User wants to build logical cuts with natural (positive) logic
  //while AODAN use inernally reverse (returns true when rejected)
  if (fFirst->Rejected(p) && fSecond->Rejected(p)) return kTRUE;//rejected (both rejected, returned kTRUE)
  return kFALSE;//accepted, at least one accepted (returned kFALSE)
}
/******************************************************************/

ClassImp(AliAODAndCut)

Bool_t AliAODAndCut::Rejected(AliVAODParticle * p)  const
{
  //returns true when rejected 
  //OR operation is a little bit misleading but is correct
  //User wants to build logical cuts with natural (positive) logic
  //while AODAN use inernally reverse (returns true when rejected)
  if (fFirst->Rejected(p) || fSecond->Rejected(p)) return kTRUE;//rejected (any of two rejected(returned kTRUE) )
  return kFALSE;//accepted (both accepted (returned kFALSE))
}
/******************************************************************/
