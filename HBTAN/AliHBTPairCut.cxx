/* $Id$ */

//-------------------------------------------------------------------
// Class AliHBTPairCut:
// implements cut on the pair of particles
// more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html
// Author: Piotr.Skowronski@cern.ch
//-------------------------------------------------------------------

#include "AliHBTPairCut.h"
#include "AliHBTPair.h"
#include "AliHBTParticleCut.h"

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

ClassImp(AliHBTEmptyPairCut)
  
void AliHBTEmptyPairCut::Streamer(TBuffer &b)
{
//streamer for empty pair cut
  AliHBTPairCut::Streamer(b);
}

ClassImp(AliHbtBasePairCut)
ClassImp(AliHBTQInvCut)
ClassImp(AliHBTKtCut)
ClassImp(AliHBTQSideCMSLCCut)
ClassImp(AliHBTQOutCMSLCCut)
ClassImp(AliHBTQLongCMSLCCut)
