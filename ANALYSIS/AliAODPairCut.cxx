/* $Id$ */
//____________________________________
/////////////////////////////////////////////////////////////////////////
//
// Class AliAODPairCut:
//
// implements cut on the pair of particles
// more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html
// Author: Piotr.Skowronski@cern.ch
//-------------------------------------------------------------------

#include "AliAODPairCut.h"
#include "AliAODPair.h"
#include "AliAODParticleCut.h"
#include "AliTrackPoints.h"
#include "AliClusterMap.h"

ClassImp(AliAODPairCut)
const Int_t AliAODPairCut::fgkMaxCuts = 50;
/**********************************************************/

AliAODPairCut::AliAODPairCut():
  fNCuts(0)
{
  //constructor
  fFirstPartCut = new AliAODEmptyParticleCut(); //empty cuts
  fSecondPartCut= new AliAODEmptyParticleCut(); //empty cuts
    
  fCuts = new AliAODBasePairCut*[fgkMaxCuts];
  for (Int_t i = 0;i<fNCuts;i++)
   {
     fCuts[i] = 0x0;
   }
}
/**********************************************************/

AliAODPairCut::AliAODPairCut(const AliAODPairCut& in):
 TNamed(in)
{
  //copy constructor
  fCuts = new AliAODBasePairCut*[fgkMaxCuts];
  fNCuts = in.fNCuts;

  fFirstPartCut = (AliAODParticleCut*)in.fFirstPartCut->Clone();
  fSecondPartCut = (AliAODParticleCut*)in.fSecondPartCut->Clone();
 
  for (Int_t i = 0;i<fNCuts;i++)
    {
      fCuts[i] = (AliAODBasePairCut*)in.fCuts[i]->Clone();//create new object (clone) and rember pointer to it
    }
}
/**********************************************************/

AliAODPairCut&  AliAODPairCut::operator=(const AliAODPairCut& in)
{
  //assignment operator
  fCuts = new AliAODBasePairCut*[fgkMaxCuts];
  fNCuts = in.fNCuts;

  fFirstPartCut = (AliAODParticleCut*)in.fFirstPartCut->Clone();
  fSecondPartCut = (AliAODParticleCut*)in.fSecondPartCut->Clone();
 
  for (Int_t i = 0;i<fNCuts;i++)
    {
      fCuts[i] = (AliAODBasePairCut*)in.fCuts[i]->Clone();//create new object (clone) and rember pointer to it
    }
  return * this;
}
/**********************************************************/

AliAODPairCut::~AliAODPairCut()
{
  //destructor
  if (fFirstPartCut != fSecondPartCut)
    {
      delete fSecondPartCut;
    }
  delete fFirstPartCut;
  for (Int_t i = 0;i<fNCuts;i++)
    {
      delete fCuts[i];
    }
  delete []fCuts;
} 
/**********************************************************/

/**********************************************************/

void AliAODPairCut::AddBasePairCut(AliAODBasePairCut* basecut)
{
  //adds the base pair cut (cut on one value)
  
  if (!basecut) return;
  if( fNCuts == (fgkMaxCuts-1) )
    {
      Warning("AddBasePairCut","Not enough place for another cut");
      return;
    }
  fCuts[fNCuts++]=basecut;
}
/**********************************************************/

Bool_t AliAODPairCut::Pass(AliAODPair* pair) const
{
  //methods which checks if given pair meets all criteria of the cut
  //if it meets returns FALSE
  //if NOT   returns    TRUE
  if(!pair) 
    {
      Warning("Pass","No Pasaran! We never accept NULL pointers");
      return kTRUE;
    }
  
  //check particle's cuts
  if( ( fFirstPartCut->Pass( pair->Particle1()) ) || 
      ( fSecondPartCut->Pass(pair->Particle2()) )   )
    {  
      return kTRUE;
    }
  return PassPairProp(pair);
}
/**********************************************************/

Bool_t AliAODPairCut::PassPairProp(AliAODPair* pair) const
{
  //methods which checks if given pair meets all criteria of the cut
  //if it meets returns FALSE
  //if NOT   returns    TRUE
  //examine all base pair cuts
  for (Int_t i = 0;i<fNCuts;i++)
    {
      if ( (fCuts[i]->Pass(pair)) ) return kTRUE; //if one of the cuts reject, then reject
    }
  return kFALSE;
}
/**********************************************************/

void AliAODPairCut::Print()
{
 //Prints the cut
  for (Int_t i = 0;i<fNCuts;i++)
    {
      fCuts[i]->Dump();
    }
}
/**********************************************************/

void AliAODPairCut::SetFirstPartCut(AliAODParticleCut* cut)
{
  // set cut for the first particle
  if(!cut) 
    {
      Error("SetFirstPartCut","argument is NULL");
      return;
    }
  delete fFirstPartCut;
  fFirstPartCut = (AliAODParticleCut*)cut->Clone();
  
}
/**********************************************************/

void AliAODPairCut::SetSecondPartCut(AliAODParticleCut* cut)
{
  // set cut for the second particle
  if(!cut) 
    {
      Error("SetSecondPartCut","argument is NULL");
      return;
    }
  delete fSecondPartCut;
  fSecondPartCut = (AliAODParticleCut*)cut->Clone();
}
/**********************************************************/

void AliAODPairCut::SetPartCut(AliAODParticleCut* cut)
{
  //sets the the same cut on both particles
  if(!cut) 
    {
      Error("SetFirstPartCut","argument is NULL");
      return;
    }
  if (fFirstPartCut == fSecondPartCut) fSecondPartCut = 0x0;
  
  delete fFirstPartCut;
  fFirstPartCut = (AliAODParticleCut*)cut->Clone();
  
  delete fSecondPartCut; //even if null should not be harmful
  fSecondPartCut = fFirstPartCut;
}
/**********************************************************/

void AliAODPairCut::SetQInvRange(Double_t min, Double_t max)
{
  // set range of accepted invariant masses
  AliAODQInvCut* cut= (AliAODQInvCut*)FindCut(kHbtPairCutPropQInv);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODQInvCut(min,max);
}
/**********************************************************/
void AliAODPairCut::SetQOutCMSLRange(Double_t min, Double_t max)
{
  // set range of accepted QOut in CMS
  AliAODQOutLCMSCut* cut= (AliAODQOutLCMSCut*)FindCut(kHbtPairCutPropQOutLCMS);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODQOutLCMSCut(min,max);
}

/**********************************************************/
void AliAODPairCut::SetQSideCMSLRange(Double_t min, Double_t max)
{
  // set range of accepted QSide in CMS
  AliAODQSideLCMSCut* cut= (AliAODQSideLCMSCut*)FindCut(kHbtPairCutPropQSideLCMS);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODQSideLCMSCut(min,max);
}

/**********************************************************/
void AliAODPairCut::SetQLongCMSLRange(Double_t min, Double_t max)
{
  // set range of accepted QLong in CMS
  AliAODQLongLCMSCut* cut= (AliAODQLongLCMSCut*)FindCut(kHbtPairCutPropQLongLCMS);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODQLongLCMSCut(min,max);
}

/**********************************************************/

void AliAODPairCut::SetKtRange(Double_t min, Double_t max)
{
  // set range of accepted Kt (?)
  AliAODKtCut* cut= (AliAODKtCut*)FindCut(kHbtPairCutPropKt);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODKtCut(min,max);
}
/**********************************************************/

void AliAODPairCut::SetKStarRange(Double_t min, Double_t max)
{
  // set range of accepted KStar (?)
  AliAODKStarCut* cut= (AliAODKStarCut*)FindCut(kHbtPairCutPropKStar);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODKStarCut(min,max);
}
/**********************************************************/

void AliAODPairCut::SetAvSeparationRange(Double_t min, Double_t max)
{
  //sets avarage separation cut ->Anti-Merging cut
  AliAODBasePairCut* cut= FindCut(kHbtPairCutPropAvSepar);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODAvSeparationCut(min,max);
}
/**********************************************************/

void AliAODPairCut::SetITSSeparation(Int_t layer, Double_t drphi, Double_t dz)
{
  //Anti-Merging Cut for first pixel layer
  AliAODITSSeparationCut* cut= dynamic_cast<AliAODITSSeparationCut*>(FindCut(kHbtPairCutPropPixelSepar));
  if(cut) 
   {
     if (layer == cut->GetLayer())
      {
        cut->SetRange(drphi,dz);//In this cut fMin is drphi, and fMax dz
        return;
      }
   }
  fCuts[fNCuts++] = new AliAODITSSeparationCut(layer,drphi,dz);
//  Info("SetITSSeparation","Added %d at address %#x",fNCuts-1,fCuts[fNCuts-1]);
}
/**********************************************************/

void AliAODPairCut::SetClusterOverlapRange(Double_t min,Double_t max)
{
  //sets cluster overlap factor cut ->Anti-Splitting cut
  //cluster overlap factor ranges between 
  // -0.5 (in all padrows both tracks have cluters) 
  // and 1 (in all padrows one track has cluter and second has not)
  // When Overlap Factor is 1 this pair of tracks in highly probable to be
  // splitted track: one particle that is recontructed twise
  // STAR uses range from -0.5 to 0.6 
  
  AliAODBasePairCut* cut= FindCut(kHbtPairCutPropClOverlap);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODCluterOverlapCut(min,max);
}
/**********************************************************/

AliAODBasePairCut* AliAODPairCut::FindCut(AliAODPairCutProperty property)
{
  // Find the cut corresponding to "property"
  for (Int_t i = 0;i<fNCuts;i++)
    {
      if (fCuts[i]->GetProperty() == property) 
	return fCuts[i]; //we found the cut we were searching for
    }
  
  return 0x0; //we did not found this cut
  
}
/**********************************************************/

void AliAODPairCut::Streamer(TBuffer &b)
{
  // Stream all objects in the array to or from the I/O buffer.
  
  UInt_t R__s, R__c;
  if (b.IsReading()) 
    {
      Version_t v = b.ReadVersion(&R__s, &R__c);
      if (v > -1)
       {
          delete fFirstPartCut;
          delete fSecondPartCut;
          fFirstPartCut = 0x0;
          fSecondPartCut = 0x0;
          TObject::Streamer(b);
          b >> fFirstPartCut;
          b >> fSecondPartCut;
          b >> fNCuts;
          for (Int_t i = 0;i<fNCuts;i++)
           {
             b >> fCuts[i];
           }
        }
      b.CheckByteCount(R__s, R__c,AliAODPairCut::IsA());
    } 
  else 
    {
      R__c = b.WriteVersion(AliAODPairCut::IsA(), kTRUE);
      TObject::Streamer(b);
      
//      printf("Streamer Cut 1 %#x Cut 2 %#x\n",fFirstPartCut,fSecondPartCut);
//      this->Dump();
//      fFirstPartCut->Dump();
      
      b << fFirstPartCut;
      b << fSecondPartCut;
      b << fNCuts;
      for (Int_t i = 0;i<fNCuts;i++)
        {
          b << fCuts[i];
        }
      b.SetByteCount(R__c, kTRUE);
    }
}
/******************************************************************/

ClassImp(AliAODEmptyPairCut)
  
void AliAODEmptyPairCut::Streamer(TBuffer &b)
{
//streamer for empty pair cut
  AliAODPairCut::Streamer(b);
}
/******************************************************************/

ClassImp(AliAODBasePairCut)
ClassImp(AliAODQInvCut)
ClassImp(AliAODKtCut)
ClassImp(AliAODQSideLCMSCut)
ClassImp(AliAODQOutLCMSCut)
ClassImp(AliAODQLongLCMSCut)

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

Bool_t AliAODITSSeparationCut::Pass(AliAODPair* pair) const
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

Bool_t AliAODOutSideSameSignCut::Pass(AliAODPair *p) const
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

Bool_t AliAODOutSideDiffSignCut::Pass(AliAODPair *p) const
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
 AliAODBasePairCut(-10e10,10e10,kHbtPairCutPropNone),
 fFirst(new AliAODDummyBasePairCut),
 fSecond(new AliAODDummyBasePairCut)
{
 //ctor
}
/******************************************************************/

AliAODLogicalOperPairCut::AliAODLogicalOperPairCut(AliAODBasePairCut* first, AliAODBasePairCut* second):
 AliAODBasePairCut(-10e10,10e10,kHbtPairCutPropNone),
 fFirst((first)?(AliAODBasePairCut*)first->Clone():0x0),
 fSecond((second)?(AliAODBasePairCut*)second->Clone():0x0)
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

Bool_t AliAODLogicalOperPairCut::AliAODDummyBasePairCut::Pass(AliAODPair* /*pair*/)  const
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

Bool_t AliAODOrPairCut::Pass(AliAODPair * p) const
{
  //returns true when rejected 
  //AND operation is a little bit misleading but is correct
  //User wants to build logical cuts with natural (positive) logic
  //while ALIAN use inernally reverse (returns true when rejected)
  if (fFirst->Pass(p) && fSecond->Pass(p) ) return kTRUE;//rejected (both rejected, returned kTRUE)
  return kFALSE;//accepted, at least one accepted (returned kFALSE)
}
/******************************************************************/

ClassImp(AliAODAndPairCut)

Bool_t AliAODAndPairCut::Pass(AliAODPair * p)  const
{
  //returns true when rejected 
  //OR operation is a little bit misleading but is correct
  //User wants to build logical cuts with natural (positive) logic
  //while ALIAN use inernally reverse (returns true when rejected)
  if (fFirst->Pass(p) || fSecond->Pass(p)) return kTRUE;//rejected (any of two rejected(returned kTRUE) )
  return kFALSE;//accepted (both accepted (returned kFALSE))
}
/******************************************************************/
