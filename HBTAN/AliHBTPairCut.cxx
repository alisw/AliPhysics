/* $Id$ */
//____________________________________
/////////////////////////////////////////////////////////////////////////
//
// Class AliHBTPairCut:
//
// implements cut on the pair of particles
// more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html
// Author: Piotr.Skowronski@cern.ch
//-------------------------------------------------------------------

#include "AliHBTPairCut.h"
#include "AliHBTPair.h"
#include "AliHBTParticleCut.h"
#include "AliHBTTrackPoints.h"
#include "AliHBTClusterMap.h"

ClassImp(AliHBTPairCut)
const Int_t AliHBTPairCut::fgkMaxCuts = 50;
/**********************************************************/

AliHBTPairCut::AliHBTPairCut():
  fNCuts(0)
{
  //constructor
  fFirstPartCut = new AliHBTEmptyParticleCut(); //empty cuts
  fSecondPartCut= new AliHBTEmptyParticleCut(); //empty cuts
    
  fCuts = new AliHbtBasePairCut*[fgkMaxCuts];
}
/**********************************************************/

AliHBTPairCut::AliHBTPairCut(const AliHBTPairCut& in):
 TNamed(in)
{
  //copy constructor
  fCuts = new AliHbtBasePairCut*[fgkMaxCuts];
  fNCuts = in.fNCuts;

  fFirstPartCut = (AliHBTParticleCut*)in.fFirstPartCut->Clone();
  fSecondPartCut = (AliHBTParticleCut*)in.fSecondPartCut->Clone();
 
  for (Int_t i = 0;i<fNCuts;i++)
    {
      fCuts[i] = (AliHbtBasePairCut*)in.fCuts[i]->Clone();//create new object (clone) and rember pointer to it
    }
}
/**********************************************************/

AliHBTPairCut&  AliHBTPairCut::operator=(const AliHBTPairCut& in)
{
  //assignment operator
  fCuts = new AliHbtBasePairCut*[fgkMaxCuts];
  fNCuts = in.fNCuts;

  fFirstPartCut = (AliHBTParticleCut*)in.fFirstPartCut->Clone();
  fSecondPartCut = (AliHBTParticleCut*)in.fSecondPartCut->Clone();
 
  for (Int_t i = 0;i<fNCuts;i++)
    {
      fCuts[i] = (AliHbtBasePairCut*)in.fCuts[i]->Clone();//create new object (clone) and rember pointer to it
    }
  return * this;
}
/**********************************************************/

AliHBTPairCut::~AliHBTPairCut()
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

void AliHBTPairCut::AddBasePairCut(AliHbtBasePairCut* basecut)
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

Bool_t AliHBTPairCut::Pass(AliHBTPair* pair) const
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

Bool_t AliHBTPairCut::PassPairProp(AliHBTPair* pair) const
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

void AliHBTPairCut::Print()
{
 //Prints the cut
  for (Int_t i = 0;i<fNCuts;i++)
    {
      fCuts[i]->Dump();
    }
}
/**********************************************************/

void AliHBTPairCut::SetFirstPartCut(AliHBTParticleCut* cut)
{
  // set cut for the first particle
  if(!cut) 
    {
      Error("SetFirstPartCut","argument is NULL");
      return;
    }
  delete fFirstPartCut;
  fFirstPartCut = (AliHBTParticleCut*)cut->Clone();
  
}
/**********************************************************/

void AliHBTPairCut::SetSecondPartCut(AliHBTParticleCut* cut)
{
  // set cut for the second particle
  if(!cut) 
    {
      Error("SetSecondPartCut","argument is NULL");
      return;
    }
  delete fSecondPartCut;
  fSecondPartCut = (AliHBTParticleCut*)cut->Clone();
}
/**********************************************************/

void AliHBTPairCut::SetPartCut(AliHBTParticleCut* cut)
{
  //sets the the same cut on both particles
  if(!cut) 
    {
      Error("SetFirstPartCut","argument is NULL");
      return;
    }
  if (fFirstPartCut == fSecondPartCut) fSecondPartCut = 0x0;
  
  delete fFirstPartCut;
  fFirstPartCut = (AliHBTParticleCut*)cut->Clone();
  
  delete fSecondPartCut; //even if null should not be harmful
  fSecondPartCut = fFirstPartCut;
}
/**********************************************************/

void AliHBTPairCut::SetQInvRange(Double_t min, Double_t max)
{
  // set range of accepted invariant masses
  AliHBTQInvCut* cut= (AliHBTQInvCut*)FindCut(kHbtPairCutPropQInv);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliHBTQInvCut(min,max);
}
/**********************************************************/
void AliHBTPairCut::SetQOutCMSLRange(Double_t min, Double_t max)
{
  // set range of accepted QOut in CMS
  AliHBTQOutCMSLCCut* cut= (AliHBTQOutCMSLCCut*)FindCut(kHbtPairCutPropQOutCMSLC);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliHBTQOutCMSLCCut(min,max);
}

/**********************************************************/
void AliHBTPairCut::SetQSideCMSLRange(Double_t min, Double_t max)
{
  // set range of accepted QSide in CMS
  AliHBTQSideCMSLCCut* cut= (AliHBTQSideCMSLCCut*)FindCut(kHbtPairCutPropQSideCMSLC);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliHBTQSideCMSLCCut(min,max);
}

/**********************************************************/
void AliHBTPairCut::SetQLongCMSLRange(Double_t min, Double_t max)
{
  // set range of accepted QLong in CMS
  AliHBTQLongCMSLCCut* cut= (AliHBTQLongCMSLCCut*)FindCut(kHbtPairCutPropQLongCMSLC);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliHBTQLongCMSLCCut(min,max);
}

/**********************************************************/

void AliHBTPairCut::SetKtRange(Double_t min, Double_t max)
{
  // set range of accepted Kt (?)
  AliHBTKtCut* cut= (AliHBTKtCut*)FindCut(kHbtPairCutPropKt);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliHBTKtCut(min,max);
}
/**********************************************************/

void AliHBTPairCut::SetKStarRange(Double_t min, Double_t max)
{
  // set range of accepted KStar (?)
  AliHBTKStarCut* cut= (AliHBTKStarCut*)FindCut(kHbtPairCutPropKStar);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliHBTKStarCut(min,max);
}
/**********************************************************/

void AliHBTPairCut::SetAvSeparationRange(Double_t min, Double_t max)
{
  //sets avarage separation cut ->Anti-Merging cut
  AliHbtBasePairCut* cut= FindCut(kHbtPairCutPropAvSepar);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliHBTAvSeparationCut(min,max);
}
/**********************************************************/

void AliHBTPairCut::SetPixelSeparation(Double_t drphi, Double_t dz)
{
  //Anti-Merging Cut for first pixel layer
  AliHbtBasePairCut* cut= FindCut(kHbtPairCutPropPixelSepar);
  if(cut) cut->SetRange(drphi,dz);//In this cut fMin is drphi, and fMax dz
  else fCuts[fNCuts++] = new AliHBTPixelSeparationCut(drphi,dz);
}
/**********************************************************/

void AliHBTPairCut::SetClusterOverlapRange(Double_t min,Double_t max)
{
  //sets cluster overlap factor cut ->Anti-Splitting cut
  //cluster overlap factor ranges between 
  // -0.5 (in all padrows both tracks have cluters) 
  // and 1 (in all padrows one track has cluter and second has not)
  // When Overlap Factor is 1 this pair of tracks in highly probable to be
  // splitted track: one particle that is recontructed twise
  // STAR uses range from -0.5 to 0.6 
  
  AliHbtBasePairCut* cut= FindCut(kHbtPairCutPropClOverlap);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliHBTCluterOverlapCut(min,max);
}
/**********************************************************/

AliHbtBasePairCut* AliHBTPairCut::FindCut(AliHBTPairCutProperty property)
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

void AliHBTPairCut::Streamer(TBuffer &b)
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
      b.CheckByteCount(R__s, R__c,AliHBTPairCut::IsA());
    } 
  else 
    {
      R__c = b.WriteVersion(AliHBTPairCut::IsA(), kTRUE);
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

ClassImp(AliHBTEmptyPairCut)
  
void AliHBTEmptyPairCut::Streamer(TBuffer &b)
{
//streamer for empty pair cut
  AliHBTPairCut::Streamer(b);
}
/******************************************************************/

ClassImp(AliHbtBasePairCut)
ClassImp(AliHBTQInvCut)
ClassImp(AliHBTKtCut)
ClassImp(AliHBTQSideCMSLCCut)
ClassImp(AliHBTQOutCMSLCCut)
ClassImp(AliHBTQLongCMSLCCut)

/******************************************************************/
ClassImp(AliHBTAvSeparationCut)
    
Double_t AliHBTAvSeparationCut::GetValue(AliHBTPair* pair) const 
{
  //chacks if avarage distance of two tracks is in given range
  AliHBTTrackPoints* tpts1 = pair->Particle1()->GetTrackPoints();
  if ( tpts1 == 0x0)
   {//it could be simulated pair
//     Warning("GetValue","Track 1 does not have Track Points. Pair NOT Passed.");
     return -1.0;
   }

  AliHBTTrackPoints* tpts2 = pair->Particle2()->GetTrackPoints();
  if ( tpts2 == 0x0)
   {
//     Warning("GetValue","Track 2 does not have Track Points. Pair NOT Passed.");
     return -1.0;
   }
   
  return tpts1->AvarageDistance(*tpts2);
}
/******************************************************************/

ClassImp(AliHBTPixelSeparationCut)

Bool_t AliHBTPixelSeparationCut::Pass(AliHBTPair* pair) const
{
 //Checks if two tracks do not cross first pixels too close to each other
 //If two tracks use the same cluster in pixels they are given
 //the same position what skews theta angles (both are the same)
 //These guys create artificial correlation in non-id analyses
 //which is positive for identical polar angles (Qlong=0) 
 //and negative for a little bit different theta angle (Qlong=epsilon)
 //Such tracks "attracks" each other.
 
  AliHBTTrackPoints* tpts1 = pair->Particle1()->GetITSTrackPoints();
  if ( tpts1 == 0x0)
   {//it could be simulated pair
     Warning("Pass","Track 1 does not have ITS Track Points. Pair NOT Passed.");
     return kTRUE;//reject 
   }

  AliHBTTrackPoints* tpts2 = pair->Particle2()->GetITSTrackPoints();
  if ( tpts2 == 0x0)
   {
     Warning("Pass","Track 2 does not have ITS Track Points. Pair NOT Passed.");
     return kTRUE;//reject 
   }
  Float_t  x1=0.0,y1=0.0,z1=0.0,x2=0.0,y2=0.0,z2=0.0;
  tpts1->PositionAt(0,x1,y1,z1);
  tpts2->PositionAt(0,x2,y2,z2);
  
//  Info("Pass","rphi %f z %f",fMin,fMax);
//  Info("Pass","P1: %f %f %f", x1,y1,z1);
//  Info("Pass","P2: %f %f %f", x2,y2,z2);
  
  Double_t dz = TMath::Abs(z1-z2);
  
  //fMax encodes treshold valaue of distance in Z
  if (dz > fMax) return kFALSE;//pair accepted
  
  Double_t drphi = TMath::Hypot(x1-x2,y1-y2);
  
  //fMin encodes treshold valaue of distance in r-phi
  if (drphi > fMin) return kFALSE;
  
  Info("Pass","Rejected !!!!!");
  return kTRUE;//they are too close, rejected
}
/******************************************************************/

ClassImp(AliHBTCluterOverlapCut)

Double_t  AliHBTCluterOverlapCut::GetValue(AliHBTPair* pair) const
{
  //Returns Cluter Overlap Factor
  //It ranges between -0.5 (in all padrows both tracks have cluters) 
  // and 1 (in all padrows one track has cluter and second has not)
  // When Overlap Factor is 1 this pair of tracks in highly probable to be
  // splitted track: one particle that is recontructed twise

  AliHBTClusterMap* cm1 = pair->Particle1()->GetClusterMap();
  if ( cm1 == 0x0)
   {
     Warning("GetValue","Track 1 does not have Cluster Map. Returning -0.5.");
     return -.5;
   }

  AliHBTClusterMap* cm2 = pair->Particle2()->GetClusterMap();
  if ( cm2 == 0x0)
   {
     Warning("GetValue","Track 2 does not have Cluster Map. Returning -0.5.");
     return -.5;
   }
  return cm1->GetOverlapFactor(*cm2);
}
/******************************************************************/
ClassImp(AliHBTOutSideSameSignCut)

Bool_t AliHBTOutSideSameSignCut::Pass(AliHBTPair *p) const
{
  //returns kTRUE if pair DO NOT meet cut criteria
  
  if ( p->GetQOutCMSLC()*p->GetQSideCMSLC() > 0 ) 
   {
     return kFALSE;//accpeted
   }

  return kTRUE ;//rejected
}
/******************************************************************/
ClassImp(AliHBTOutSideDiffSignCut)

Bool_t AliHBTOutSideDiffSignCut::Pass(AliHBTPair *p) const
{
  //returns kTRUE if pair DO NOT meet cut criteria
  
  if ( p->GetQOutCMSLC()*p->GetQSideCMSLC() > 0 ) 
   {
     return kTRUE;//rejected
   }
  
  return kFALSE;//accepted
}
/******************************************************************/
ClassImp( AliHBTLogicalOperPairCut )

AliHBTLogicalOperPairCut::AliHBTLogicalOperPairCut():
 AliHbtBasePairCut(-10e10,10e10,kHbtPairCutPropNone),
 fFirst(new AliHBTDummyBasePairCut),
 fSecond(new AliHBTDummyBasePairCut)
{
 //ctor
}
/******************************************************************/

AliHBTLogicalOperPairCut::AliHBTLogicalOperPairCut(AliHbtBasePairCut* first, AliHbtBasePairCut* second):
 AliHbtBasePairCut(-10e10,10e10,kHbtPairCutPropNone),
 fFirst((first)?(AliHbtBasePairCut*)first->Clone():0x0),
 fSecond((second)?(AliHbtBasePairCut*)second->Clone():0x0)
{
  //ctor
  //note that base cuts are copied, not just pointers assigned
  if ( (fFirst && fSecond) == kFALSE) 
   {
     Fatal("AliHBTLogicalOperPairCut","One of parameters is NULL!");
   }
}
/******************************************************************/

AliHBTLogicalOperPairCut::~AliHBTLogicalOperPairCut()
{
  //destructor
  delete fFirst;
  delete fSecond;
}
/******************************************************************/

Bool_t AliHBTLogicalOperPairCut::AliHBTDummyBasePairCut::Pass(AliHBTPair* /*pair*/)  const
{
  //checks if particles passes properties defined by this cut
  Warning("Pass","You are using dummy base cut! Probobly some logical cut is not set up properly");
  return kFALSE;//accept
}
/******************************************************************/

void AliHBTLogicalOperPairCut::Streamer(TBuffer &b)
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
     b.CheckByteCount(R__s, R__c,AliHBTLogicalOperPairCut::IsA());
   } 
  else 
   {
     R__c = b.WriteVersion(AliHBTLogicalOperPairCut::IsA(), kTRUE);
     TObject::Streamer(b);
     b << fFirst;
     b << fSecond;
     b.SetByteCount(R__c, kTRUE);
  }
}

/******************************************************************/
ClassImp(AliHBTOrPairCut)

Bool_t AliHBTOrPairCut::Pass(AliHBTPair * p) const
{
  //returns true when rejected 
  //AND operation is a little bit misleading but is correct
  //User wants to build logical cuts with natural (positive) logic
  //while HBTAN use inernally reverse (returns true when rejected)
  if (fFirst->Pass(p) && fSecond->Pass(p) ) return kTRUE;//rejected (both rejected, returned kTRUE)
  return kFALSE;//accepted, at least one accepted (returned kFALSE)
}
/******************************************************************/

ClassImp(AliHBTAndPairCut)

Bool_t AliHBTAndPairCut::Pass(AliHBTPair * p)  const
{
  //returns true when rejected 
  //OR operation is a little bit misleading but is correct
  //User wants to build logical cuts with natural (positive) logic
  //while HBTAN use inernally reverse (returns true when rejected)
  if (fFirst->Pass(p) || fSecond->Pass(p)) return kTRUE;//rejected (any of two rejected(returned kTRUE) )
  return kFALSE;//accepted (both accepted (returned kFALSE))
}
/******************************************************************/
