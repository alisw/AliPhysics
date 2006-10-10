#include "AliAODPairCut.h"
/* $Id$ */

/////////////////////////////////////////////////////////////////////////
//
// Class AliAODPairCut:
// implements cut on the pair of particles
// more info: http://aliweb.cern.ch/people/skowron/analyzer/index.html
// Author: Piotr.Skowronski@cern.ch
//-------------------------------------------------------------------

#include "AliAODPair.h"
#include "AliAODParticleCut.h"
//#include "AliTrackPoints.h"
//#include "AliClusterMap.h"

ClassImp(AliAODPairCut)
const Int_t AliAODPairCut::fgkMaxCuts = 50;
/**********************************************************/

AliAODPairCut::AliAODPairCut():
  fFirstPartCut(new AliAODParticleEmptyCut()), //empty cuts
  fSecondPartCut(new AliAODParticleEmptyCut()), //empty cuts
  fCuts(new AliAODPairBaseCut*[fgkMaxCuts]),
  fNCuts(0)
{
  //constructor
    
  for (Int_t i = 0;i<fNCuts;i++)
   {
     fCuts[i] = 0x0;
   }
}
/**********************************************************/

AliAODPairCut::AliAODPairCut(const AliAODPairCut& in):
  TNamed(in),
  fFirstPartCut((AliAODParticleCut*)in.fFirstPartCut->Clone()),
  fSecondPartCut((AliAODParticleCut*)in.fSecondPartCut->Clone()),
  fCuts(new AliAODPairBaseCut*[fgkMaxCuts]),
  fNCuts(in.fNCuts)
{
  //copy constructor

 
  for (Int_t i = 0;i<fNCuts;i++)
    {
      fCuts[i] = (AliAODPairBaseCut*)in.fCuts[i]->Clone();//create new object (clone) and rember pointer to it
    }
}
/**********************************************************/

AliAODPairCut&  AliAODPairCut::operator=(const AliAODPairCut& in)
{
  //assignment operator
  fCuts = new AliAODPairBaseCut*[fgkMaxCuts];
  fNCuts = in.fNCuts;

  fFirstPartCut = (AliAODParticleCut*)in.fFirstPartCut->Clone();
  fSecondPartCut = (AliAODParticleCut*)in.fSecondPartCut->Clone();
 
  for (Int_t i = 0;i<fNCuts;i++)
    {
      fCuts[i] = (AliAODPairBaseCut*)in.fCuts[i]->Clone();//create new object (clone) and rember pointer to it
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

void AliAODPairCut::AddBasePairCut(AliAODPairBaseCut* basecut)
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

Bool_t AliAODPairCut::Rejected(AliAODPair* pair) const
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
  if( ( fFirstPartCut->Rejected( pair->Particle1()) ) || 
      ( fSecondPartCut->Rejected(pair->Particle2()) )   )
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
      if ( (fCuts[i]->Rejected(pair)) ) return kTRUE; //if one of the cuts reject, then reject
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
  AliAODQInvCut* cut= (AliAODQInvCut*)FindCut(AliAODPairBaseCut::kHbtPairCutPropQInv);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODQInvCut(min,max);
}
/**********************************************************/

void AliAODPairCut::SetQOutLCMSRange(Double_t min, Double_t max)
{
  // set range of accepted QOut in CMS
  AliAODQOutLCMSCut* cut= (AliAODQOutLCMSCut*)FindCut(AliAODPairBaseCut::kHbtPairCutPropQOutLCMS);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODQOutLCMSCut(min,max);
}
/**********************************************************/

void AliAODPairCut::SetQSideLCMSRange(Double_t min, Double_t max)
{
  // set range of accepted QSide in CMS
  AliAODQSideLCMSCut* cut= (AliAODQSideLCMSCut*)FindCut(AliAODPairBaseCut::kHbtPairCutPropQSideLCMS);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODQSideLCMSCut(min,max);
}

/**********************************************************/

void AliAODPairCut::SetQLongLCMSRange(Double_t min, Double_t max)
{
  // set range of accepted QLong in CMS
  AliAODQLongLCMSCut* cut= (AliAODQLongLCMSCut*)FindCut(AliAODPairBaseCut::kHbtPairCutPropQLongLCMS);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODQLongLCMSCut(min,max);
}
/**********************************************************/

void AliAODPairCut::SetDeltaERange(Double_t min, Double_t max)
{
  // set range of accepted DeltaE
  AliAODKtCut* cut= (AliAODKtCut*)FindCut(AliAODPairBaseCut::kHbtPairCutPropDeltaE);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODDeltaECut(min,max);
}
/**********************************************************/

void AliAODPairCut::SetDeltaPRange(Double_t min, Double_t max)
{
  // set range of accepted DeltaP
  AliAODKtCut* cut= (AliAODKtCut*)FindCut(AliAODPairBaseCut::kHbtPairCutPropDeltaP);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODDeltaPCut(min,max);
}
/**********************************************************/

void AliAODPairCut::SetKtRange(Double_t min, Double_t max)
{
  // set range of accepted Kt (avarage transverse pair momentum)
  AliAODKtCut* cut= (AliAODKtCut*)FindCut(AliAODPairBaseCut::kHbtPairCutPropKt);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODKtCut(min,max);
}
/**********************************************************/

void AliAODPairCut::SetKStarRange(Double_t min, Double_t max)
{
  // set range of accepted KStar (invariant pair momentum difference (fourvector))
  AliAODKStarCut* cut= (AliAODKStarCut*)FindCut(AliAODPairBaseCut::kHbtPairCutPropKStar);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODKStarCut(min,max);
}
/**********************************************************/

void AliAODPairCut::SetKStarOutRange(Double_t min, Double_t max)
{
  // set range of accepted KStar (invariant pair momentum difference (fourvector))
  AliAODKStarOutCut* cut= (AliAODKStarOutCut*)FindCut(AliAODPairBaseCut::kHbtPairCutPropKStarOut);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODKStarOutCut(min,max);
}
/**********************************************************/

void AliAODPairCut::SetKStarSideRange(Double_t min, Double_t max)
{
  // set range of accepted KStar (invariant pair momentum difference (fourvector))
  AliAODKStarSideCut* cut= (AliAODKStarSideCut*)FindCut(AliAODPairBaseCut::kHbtPairCutPropKStarSide);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODKStarSideCut(min,max);
}
/**********************************************************/

void AliAODPairCut::SetKStarLongRange(Double_t min, Double_t max)
{
  // set range of accepted KStar (invariant pair momentum difference (fourvector))
  AliAODKStarLongCut* cut= (AliAODKStarLongCut*)FindCut(AliAODPairBaseCut::kHbtPairCutPropKStarLong);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODKStarLongCut(min,max);
}
/**********************************************************/

void AliAODPairCut::SetAvSeparationRange(Double_t min, Double_t max)
{
  //sets avarage separation cut ->Anti-Merging cut
  AliAODPairBaseCut* cut= FindCut(AliAODPairBaseCut::kHbtPairCutPropAvSepar);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODAvSeparationCut(min,max);
}
/**********************************************************/

void AliAODPairCut::SetITSSeparation(Int_t layer, Double_t drphi, Double_t dz)
{
  //Anti-Merging Cut for first pixel layer
  AliAODITSSeparationCut* cut= dynamic_cast<AliAODITSSeparationCut*>(FindCut(AliAODPairBaseCut::kHbtPairCutPropPixelSepar));
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
  
  AliAODPairBaseCut* cut= FindCut(AliAODPairBaseCut::kHbtPairCutPropClOverlap);
  if(cut) cut->SetRange(min,max);
  else fCuts[fNCuts++] = new AliAODCluterOverlapCut(min,max);
}
/**********************************************************/

AliAODPairBaseCut* AliAODPairCut::FindCut(AliAODPairBaseCut::EAODPairCutProperty property)
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

ClassImp(AliAODPairEmptyCut)
  
void AliAODPairEmptyCut::Streamer(TBuffer &b)
{
//streamer for empty pair cut
  AliAODPairCut::Streamer(b);
}
/******************************************************************/

