#include "AliAODPairBaseCut.h"

#include "AliTrackPoints.h"
#include "AliClusterMap.h"


ClassImp(AliAODPairBaseCut)
ClassImp(AliAODQInvCut)
ClassImp(AliAODKtCut)
ClassImp(AliAODQSideLCMSCut)
ClassImp(AliAODQOutLCMSCut)
ClassImp(AliAODQLongLCMSCut)
ClassImp(AliAODDeltaECut)
ClassImp(AliAODDeltaPCut)
ClassImp(AliAODDeltaPvectorCut)
ClassImp(AliAODDeltaPhiCut)
ClassImp(AliAODDeltaThetaCut)

/******************************************************************/
ClassImp(AliAODAvSeparationCut)
    
Double_t AliAODAvSeparationCut::GetValue(AliAODPair* pair) const 
{
  //chacks if avarage distance of two tracks is in given range
  AliTrackPoints* tpts1 = pair->Particle1()->GetTPCTrackPoints();
  if ( tpts1 == 0x0)
   {//it could be simulated pair
//     Warning("GetValue","Track 1 does not have Track Points. Pair NOT Passed.");
     return -1.0;
   }

  AliTrackPoints* tpts2 = pair->Particle2()->GetTPCTrackPoints();
  if ( tpts2 == 0x0)
   {
//     Warning("GetValue","Track 2 does not have Track Points. Pair NOT Passed.");
     return -1.0;
   }
   
  return tpts1->AvarageDistance(*tpts2);
}
/******************************************************************/
ClassImp(AliAODSeparationCut)
    
Double_t AliAODSeparationCut::GetValue(AliAODPair* pair) const 
{
  //chacks if avarage distance of two tracks is in given range
  AliTrackPoints* tpts1 = pair->Particle1()->GetTPCTrackPoints();
  if ( tpts1 == 0x0)
   {//it could be simulated pair
//     Warning("GetValue","Track 1 does not have Track Points. Pair NOT Passed.");
     return -1.0;
   }

  AliTrackPoints* tpts2 = pair->Particle2()->GetTPCTrackPoints();
  if ( tpts2 == 0x0)
   {
//     Warning("GetValue","Track 2 does not have Track Points. Pair NOT Passed.");
     return -1.0;
   }
  Float_t x1=0,y1=0,z1=0; 
  Float_t x2=0,y2=0,z2=0;
  
  tpts1->PositionAt(fPoint,x1,y1,z1);
  tpts2->PositionAt(fPoint,x2,y2,z2);
  Double_t dx1 = x1 - x2;
  Double_t dy1 = y1 - y2;
  Double_t dz1 = z1 - z2;
  Double_t d = TMath::Sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
  return d;
}
/******************************************************************/

ClassImp(AliAODITSSeparationCut)

Bool_t AliAODITSSeparationCut::Rejected(AliAODPair* pair) const
{
 //Checks if two tracks do not cross first pixels too close to each other
 //If two tracks use the same cluster in pixels they are given
 //the same position what skews theta angles (both are the same)
 //These guys create artificial correlation in non-id analyses
 //which is positive for identical polar angles (Qlong=0) 
 //and negative for a little bit different theta angle (Qlong=epsilon)
 //Such tracks "attracks" each other.
 
  AliTrackPoints* tpts1 = pair->Particle1()->GetITSTrackPoints();
  if ( tpts1 == 0x0)
   {//it could be simulated pair
     Warning("Pass","Track 1 does not have ITS Track Points. Pair NOT Passed.");
     return kTRUE;//reject 
   }

  AliTrackPoints* tpts2 = pair->Particle2()->GetITSTrackPoints();
  if ( tpts2 == 0x0)
   {
     Warning("Pass","Track 2 does not have ITS Track Points. Pair NOT Passed.");
     return kTRUE;//reject 
   }
  Float_t  x1=0.0,y1=0.0,z1=0.0,x2=0.0,y2=0.0,z2=0.0;
  tpts1->PositionAt(fLayer,x1,y1,z1);
  tpts2->PositionAt(fLayer,x2,y2,z2);
  
//  Info("Pass","rphi %f z %f",fMin,fMax);
//  Info("Pass","P1: %f %f %f", x1,y1,z1);
//  Info("Pass","P2: %f %f %f", x2,y2,z2);
  
  Double_t dz = TMath::Abs(z1-z2);
  
  //fMax encodes treshold valaue of distance in Z
  if (dz > fMax) return kFALSE;//pair accepted
  
  Double_t drphi = TMath::Hypot(x1-x2,y1-y2);
  
  //fMin encodes treshold valaue of distance in r-phi
  if (drphi > fMin) return kFALSE;
  
  return kTRUE;//they are too close, rejected
}
/******************************************************************/

ClassImp(AliAODCluterOverlapCut)

Double_t  AliAODCluterOverlapCut::GetValue(AliAODPair* pair) const
{
  //Returns Cluter Overlap Factor
  //It ranges between -0.5 (in all padrows both tracks have cluters) 
  // and 1 (in all padrows one track has cluter and second has not)
  // When Overlap Factor is 1 this pair of tracks in highly probable to be
  // splitted track: one particle that is recontructed twise

  AliClusterMap* cm1 = pair->Particle1()->GetClusterMap();
  if ( cm1 == 0x0)
   {
     Warning("GetValue","Track 1 does not have Cluster Map. Returning -0.5.");
     return -.5;
   }

  AliClusterMap* cm2 = pair->Particle2()->GetClusterMap();
  if ( cm2 == 0x0)
   {
     Warning("GetValue","Track 2 does not have Cluster Map. Returning -0.5.");
     return -.5;
   }
  return cm1->GetOverlapFactor(*cm2);
}
/******************************************************************/
ClassImp(AliAODOutSideSameSignCut)

Bool_t AliAODOutSideSameSignCut::Rejected(AliAODPair *p) const
{
  //returns kTRUE if pair DO NOT meet cut criteria
  
  if ( p->GetQOutLCMS()*p->GetQSideLCMS() > 0 ) 
   {
     return kFALSE;//accpeted
   }

  return kTRUE ;//rejected
}
/******************************************************************/
ClassImp(AliAODOutSideDiffSignCut)

Bool_t AliAODOutSideDiffSignCut::Rejected(AliAODPair *p) const
{
  //returns kTRUE if pair DO NOT meet cut criteria
  
  if ( p->GetQOutLCMS()*p->GetQSideLCMS() > 0 ) 
   {
     return kTRUE;//rejected
   }
  
  return kFALSE;//accepted
}
/******************************************************************/
ClassImp( AliAODLogicalOperPairCut )

AliAODLogicalOperPairCut::AliAODLogicalOperPairCut():
 AliAODPairBaseCut(-10e10,10e10,kHbtPairCutPropNone),
 fFirst(new AliAODDummyBasePairCut),
 fSecond(new AliAODDummyBasePairCut)
{
 //ctor
}
/******************************************************************/

AliAODLogicalOperPairCut::AliAODLogicalOperPairCut(AliAODPairBaseCut* first, AliAODPairBaseCut* second):
 AliAODPairBaseCut(-10e10,10e10,kHbtPairCutPropNone),
 fFirst((first)?(AliAODPairBaseCut*)first->Clone():0x0),
 fSecond((second)?(AliAODPairBaseCut*)second->Clone():0x0)
{
  //ctor
  //note that base cuts are copied, not just pointers assigned
  if ( (fFirst && fSecond) == kFALSE) 
   {
     Fatal("AliAODLogicalOperPairCut","One of parameters is NULL!");
   }
}
/******************************************************************/

AliAODLogicalOperPairCut::~AliAODLogicalOperPairCut()
{
  //destructor
  delete fFirst;
  delete fSecond;
}
/******************************************************************/

Bool_t AliAODLogicalOperPairCut::AliAODDummyBasePairCut::Rejected(AliAODPair* /*pair*/)  const
{
  //checks if particles passes properties defined by this cut
  Warning("Pass","You are using dummy base cut! Probobly some logical cut is not set up properly");
  return kFALSE;//accept
}
/******************************************************************/

void AliAODLogicalOperPairCut::Streamer(TBuffer &b)
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
     b.CheckByteCount(R__s, R__c,AliAODLogicalOperPairCut::IsA());
   } 
  else 
   {
     R__c = b.WriteVersion(AliAODLogicalOperPairCut::IsA(), kTRUE);
     TObject::Streamer(b);
     b << fFirst;
     b << fSecond;
     b.SetByteCount(R__c, kTRUE);
  }
}

/******************************************************************/
ClassImp(AliAODOrPairCut)

Bool_t AliAODOrPairCut::Rejected(AliAODPair * p) const
{
  //returns true when rejected 
  //AND operation is a little bit misleading but is correct
  //User wants to build logical cuts with natural (positive) logic
  //while ALIAN use inernally reverse (returns true when rejected)
  if (fFirst->Rejected(p) && fSecond->Rejected(p) ) return kTRUE;//rejected (both rejected, returned kTRUE)
  return kFALSE;//accepted, at least one accepted (returned kFALSE)
}
/******************************************************************/

ClassImp(AliAODAndPairCut)

Bool_t AliAODAndPairCut::Rejected(AliAODPair * p)  const
{
  //returns true when rejected 
  //OR operation is a little bit misleading but is correct
  //User wants to build logical cuts with natural (positive) logic
  //while ALIAN use inernally reverse (returns true when rejected)
  if (fFirst->Rejected(p) || fSecond->Rejected(p)) return kTRUE;//rejected (any of two rejected(returned kTRUE) )
  return kFALSE;//accepted (both accepted (returned kFALSE))
}
/******************************************************************/
