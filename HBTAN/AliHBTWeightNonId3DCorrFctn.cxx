#include "AliHBTWeightNonId3DCorrFctn.h"
#include <TH1.h>
#include <Riostream.h>

///////////////////////////////////////////////////////
//                                                   //
// AliHBTWeightNonId3DCorrFctn.h                     //
//                                                   //
// Class for calculating 3D non-id correlation       //
// functions using methods of weights                //
//                                                   //
///////////////////////////////////////////////////////

ClassImp(AliHBTWeightNonId3DCorrFctn)

AliHBTWeightNonId3DCorrFctn::AliHBTWeightNonId3DCorrFctn(const char* name, const char* title):
 AliHBTTwoPairFctn1D(name,title),
 fWeightNumOutP(0x0),
 fWeightDenOutP(0x0),
 fWeightRatOutP(0x0),
 fWeightNumOutN(0x0),
 fWeightDenOutN(0x0),
 fWeightRatOutN(0x0),
 fWeightRatOut(0x0),
 fWeightRatOutNOverP(0x0),
 fWeightNumSideP(0x0),
 fWeightDenSideP(0x0),
 fWeightRatSideP(0x0),
 fWeightNumSideN(0x0),
 fWeightDenSideN(0x0),
 fWeightRatSideN(0x0),
 fWeightRatSide(0x0),
 fWeightRatSideNOverP(0x0),
 fWeightNumLongP(0x0),
 fWeightDenLongP(0x0),
 fWeightRatLongP(0x0),
 fWeightNumLongN(0x0),
 fWeightDenLongN(0x0),
 fWeightRatLongN(0x0),
 fWeightRatLong(0x0),
 fWeightRatLongNOverP(0x0)
{
//ctor
}
/******************************************************************/
AliHBTWeightNonId3DCorrFctn::AliHBTWeightNonId3DCorrFctn(const char* name, const char* title, 
							 Int_t nbinsX, Float_t maxXval, Float_t minXval):   
 AliHBTTwoPairFctn1D(name,title,nbinsX,maxXval,minXval),
 fWeightNumOutP(0x0),
 fWeightDenOutP(0x0),
 fWeightRatOutP(0x0),
 fWeightNumOutN(0x0),
 fWeightDenOutN(0x0),
 fWeightRatOutN(0x0),
 fWeightRatOut(0x0),
 fWeightRatOutNOverP(0x0),
 fWeightNumSideP(0x0),
 fWeightDenSideP(0x0),
 fWeightRatSideP(0x0),
 fWeightNumSideN(0x0),
 fWeightDenSideN(0x0),
 fWeightRatSideN(0x0),
 fWeightRatSide(0x0),
 fWeightRatSideNOverP(0x0),
 fWeightNumLongP(0x0),
 fWeightDenLongP(0x0),
 fWeightRatLongP(0x0),
 fWeightNumLongN(0x0),
 fWeightDenLongN(0x0),
 fWeightRatLongN(0x0),
 fWeightRatLong(0x0),
 fWeightRatLongNOverP(0x0)
{
//ctor
}

/******************************************************************/
AliHBTWeightNonId3DCorrFctn::AliHBTWeightNonId3DCorrFctn(const AliHBTWeightNonId3DCorrFctn& in):
 AliHBTTwoPairFctn1D(in),
 fWeightNumOutP(0x0),
 fWeightDenOutP(0x0),
 fWeightRatOutP(0x0),
 fWeightNumOutN(0x0),
 fWeightDenOutN(0x0),
 fWeightRatOutN(0x0),
 fWeightRatOut(0x0),
 fWeightRatOutNOverP(0x0),
 fWeightNumSideP(0x0),
 fWeightDenSideP(0x0),
 fWeightRatSideP(0x0),
 fWeightNumSideN(0x0),
 fWeightDenSideN(0x0),
 fWeightRatSideN(0x0),
 fWeightRatSide(0x0),
 fWeightRatSideNOverP(0x0),
 fWeightNumLongP(0x0),
 fWeightDenLongP(0x0),
 fWeightRatLongP(0x0),
 fWeightNumLongN(0x0),
 fWeightDenLongN(0x0),
 fWeightRatLongN(0x0),
 fWeightRatLong(0x0),
 fWeightRatLongNOverP(0x0)
{
//ctor
}

/******************************************************************/

AliHBTWeightNonId3DCorrFctn::~AliHBTWeightNonId3DCorrFctn()
{
 //dtor
 delete fWeightNumOutP;
 delete fWeightDenOutP;
 delete fWeightRatOutP;
 delete fWeightNumOutN;
 delete fWeightDenOutN;
 delete fWeightRatOutN;
 delete fWeightRatOut;
 delete fWeightRatOutNOverP;
 delete fWeightNumSideP;
 delete fWeightDenSideP;
 delete fWeightRatSideP;
 delete fWeightNumSideN;
 delete fWeightDenSideN;
 delete fWeightRatSideN;
 delete fWeightRatSide;
 delete fWeightRatSideNOverP;
 delete fWeightNumLongP;
 delete fWeightDenLongP;
 delete fWeightRatLongP;
 delete fWeightNumLongN;
 delete fWeightDenLongN;
 delete fWeightRatLongN;
 delete fWeightRatLong;
 delete fWeightRatLongNOverP;
 
}

//-------------------------------------
void AliHBTWeightNonId3DCorrFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
    //Fills the numerator using pair from the same event
    partpair = CheckPair(partpair);
    trackpair = CheckPair(trackpair);
    if(partpair == 0x0) return;
    if(trackpair == 0x0) return;
 
    double tKStar = TMath::Abs(trackpair->GetKStar());
    double tKStarOut = trackpair->GetKStarOut();
    double tKStarSide = trackpair->GetKStarSide();
    double tKStarLong = trackpair->GetKStarLong();

    Double_t weight = 1.0;
    
  if ( trackpair && partpair)
  {
    if ( ( trackpair->Particle1()->GetPdgCode() == partpair->Particle1()->GetPdgCode()) &&
         ( trackpair->Particle2()->GetPdgCode() == partpair->Particle2()->GetPdgCode())    )
    {
	weight=partpair->GetWeight();
    }
    
    if(tKStarOut>0.)
    {
	fWeightNumOutP->Fill(tKStar,weight);
    }
    else
    {
	fWeightNumOutN->Fill(tKStar,weight);
    }

    if(tKStarSide>0.)
    {
	fWeightNumSideP->Fill(tKStar,weight);
    }
    else
    {
	fWeightNumSideN->Fill(tKStar,weight);
    }

    if(tKStarLong>0.)
    {
	fWeightNumLongP->Fill(tKStar,weight);
    }
    else
    {
	fWeightNumLongN->Fill(tKStar,weight);
    }
  }
  
}

/****************************************************************/


void AliHBTWeightNonId3DCorrFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{

  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);

  double tKStar = TMath::Abs(partpair->GetKStar());
  double tKStarOut = trackpair->GetKStarOut();
  double tKStarSide = trackpair->GetKStarSide();
  double tKStarLong = trackpair->GetKStarLong();

  if ( trackpair && partpair)
  {
    if(tKStarOut>0.)
    {
	fWeightDenOutP->Fill(tKStar);
    }
    else
    {
	fWeightDenOutN->Fill(tKStar);
    }

    if(tKStarSide>0.)
    {
	fWeightDenSideP->Fill(tKStar);
    }
    else
    {
	fWeightDenSideN->Fill(tKStar);
    }

    if(tKStarLong>0.)
    {
	fWeightDenLongP->Fill(tKStar);
    }
    else
    {
	fWeightDenLongN->Fill(tKStar);
    }
  }
  
}


/****************************************************************/
void AliHBTWeightNonId3DCorrFctn::Init()
{
    AliHBTTwoPairFctn1D::Init();
 fWeightNumOutP->Reset();
 fWeightDenOutP->Reset();
 fWeightRatOutP->Reset();
 fWeightNumOutN->Reset();
 fWeightDenOutN->Reset();
 fWeightRatOutN->Reset();
 fWeightRatOut->Reset();
 fWeightRatOutNOverP->Reset();

 fWeightNumSideP->Reset();
 fWeightDenSideP->Reset();
 fWeightRatSideP->Reset();
 fWeightNumSideN->Reset();
 fWeightDenSideN->Reset();
 fWeightRatSideN->Reset();
 fWeightRatSide->Reset();
 fWeightRatSideNOverP->Reset();

 fWeightNumLongP->Reset();
 fWeightDenLongP->Reset();
 fWeightRatLongP->Reset();
 fWeightNumLongN->Reset();
 fWeightDenLongN->Reset();
 fWeightRatLongN->Reset();
 fWeightRatLong->Reset();
 fWeightRatLongNOverP->Reset();
}
/******************************************************************/

void AliHBTWeightNonId3DCorrFctn::BuildHistos(Int_t nbins, Float_t max, Float_t min)
{
    
    AliHBTFunction1D::BuildHistos(nbins,max,min);
    
    TString nameNumOutP = "WeightNumOutP";
    TString nameDenOutP = "WeightDenOutP";
    TString nameRatOutP = "WeightRatOutP";
    TString nameNumOutN = "WeightNumOutN";
    TString nameDenOutN = "WeightDenOutN";
    TString nameRatOutN = "WeightRatOutN";
    TString nameRatOut =  "WeightRatOut";
    TString nameRatOutNOverP = "WeightRatOutNOverP";
    TString nameNumSideP = "WeightNumSideP";
    TString nameDenSideP = "WeightDenSideP";
    TString nameRatSideP = "WeightRatSideP";
    TString nameNumSideN = "WeightNumSideN";
    TString nameDenSideN = "WeightDenSideN";
    TString nameRatSideN = "WeightRatSideN";
    TString nameRatSide = "WeightRatSide";
    TString nameRatSideNOverP = "WeightRatSideNOverP";
    TString nameNumLongP = "WeightNumLongP";
    TString nameDenLongP = "WeightDenLongP";
    TString nameRatLongP = "WeightRatLongP";
    TString nameNumLongN = "WeightNumLongN";
    TString nameDenLongN = "WeightDenLongN";
    TString nameRatLongN = "WeightRatLongN";
    TString nameRatLong = "WeightRatLong";
    TString nameRatLongNOverP = "WeightRatLongNOverP";

    fWeightNumOutP = new TH1D(nameNumOutP.Data(),nameNumOutP.Data(),nbins,min,max);
    fWeightDenOutP = new TH1D(nameDenOutP.Data(),nameDenOutP.Data(),nbins,min,max);
    fWeightRatOutP = new TH1D(nameRatOutP.Data(),nameRatOutP.Data(),nbins,min,max);

    fWeightNumOutN = new TH1D(nameNumOutN.Data(),nameNumOutN.Data(),nbins,min,max);
    fWeightDenOutN = new TH1D(nameDenOutN.Data(),nameDenOutN.Data(),nbins,min,max);
    fWeightRatOutN = new TH1D(nameRatOutN.Data(),nameRatOutN.Data(),nbins,min,max);

    fWeightRatOut = new TH1D(nameRatOut.Data(),nameRatOut.Data(),nbins,min,max);
    fWeightRatOutNOverP = new TH1D(nameRatOutNOverP.Data(),nameRatOutNOverP.Data(),nbins,min,max);

    fWeightNumSideP = new TH1D(nameNumSideP.Data(),nameNumSideP.Data(),nbins,min,max);
    fWeightDenSideP = new TH1D(nameDenSideP.Data(),nameDenSideP.Data(),nbins,min,max);
    fWeightRatSideP = new TH1D(nameRatSideP.Data(),nameRatSideP.Data(),nbins,min,max);

    fWeightNumSideN = new TH1D(nameNumSideN.Data(),nameNumSideN.Data(),nbins,min,max);
    fWeightDenSideN = new TH1D(nameDenSideN.Data(),nameDenSideN.Data(),nbins,min,max);
    fWeightRatSideN = new TH1D(nameRatSideN.Data(),nameRatSideN.Data(),nbins,min,max);

    fWeightRatSide = new TH1D(nameRatSide.Data(),nameRatSide.Data(),nbins,min,max);
    fWeightRatSideNOverP = new TH1D(nameRatSideNOverP.Data(),nameRatSideNOverP.Data(),nbins,min,max);

    fWeightNumLongP = new TH1D(nameNumLongP.Data(),nameNumLongP.Data(),nbins,min,max);
    fWeightDenLongP = new TH1D(nameDenLongP.Data(),nameDenLongP.Data(),nbins,min,max);
    fWeightRatLongP = new TH1D(nameRatLongP.Data(),nameRatLongP.Data(),nbins,min,max);

    fWeightNumLongN = new TH1D(nameNumLongN.Data(),nameNumLongN.Data(),nbins,min,max);
    fWeightDenLongN = new TH1D(nameDenLongN.Data(),nameDenLongN.Data(),nbins,min,max);
    fWeightRatLongN = new TH1D(nameRatLongN.Data(),nameRatLongN.Data(),nbins,min,max);

    fWeightRatLong = new TH1D(nameRatLong.Data(),nameRatLong.Data(),nbins,min,max);
    fWeightRatLongNOverP = new TH1D(nameRatLongNOverP.Data(),nameRatLongNOverP.Data(),nbins,min,max);


    fWeightNumOutP->Sumw2();
    fWeightDenOutP->Sumw2();
    fWeightRatOutP->Sumw2();
    fWeightNumOutN->Sumw2();
    fWeightDenOutN->Sumw2();
    fWeightRatOutN->Sumw2();
    fWeightRatOut->Sumw2();
    fWeightRatOutNOverP->Sumw2();
    fWeightNumSideP->Sumw2();
    fWeightDenSideP->Sumw2();
    fWeightRatSideP->Sumw2();
    fWeightNumSideN->Sumw2();
    fWeightDenSideN->Sumw2();
    fWeightRatSideN->Sumw2();
    fWeightRatSide->Sumw2();
    fWeightRatSideNOverP->Sumw2();
    fWeightNumLongP->Sumw2();
    fWeightDenLongP->Sumw2();
    fWeightRatLongP->Sumw2();
    fWeightNumLongN->Sumw2();
    fWeightDenLongN->Sumw2();
    fWeightRatLongN->Sumw2();
    fWeightRatLong->Sumw2();
    fWeightRatLongNOverP->Sumw2();
}

/******************************************************************/

TH1* AliHBTWeightNonId3DCorrFctn::GetResult()
{
    return fWeightRatOutN;
}

/******************************************************************/

void AliHBTWeightNonId3DCorrFctn::WriteFunction()
{

 Double_t outPscale = Scale(fWeightNumOutP,fWeightDenOutP);
 cout<<"WoutPscale = "<<outPscale<<endl;
 fWeightRatOutP->Divide(fWeightNumOutP,fWeightDenOutP,outPscale);

 Double_t outNscale = Scale(fWeightNumOutN,fWeightDenOutN);
 cout<<"WoutNscale = "<<outNscale<<endl;
 fWeightRatOutN->Divide(fWeightNumOutN,fWeightDenOutN,outNscale);

 fWeightRatOutNOverP->Divide(fWeightRatOutN,fWeightRatOutP);


 Double_t sidePscale = Scale(fWeightNumSideP,fWeightDenSideP);
 fWeightRatSideP->Divide(fWeightNumSideP,fWeightDenSideP,sidePscale);

 Double_t sideNscale = Scale(fWeightNumSideN,fWeightDenSideN);
 fWeightRatSideN->Divide(fWeightNumSideN,fWeightDenSideN,sideNscale);

 fWeightRatSideNOverP->Divide(fWeightRatSideN,fWeightRatSideP);


 Double_t longPscale = Scale(fWeightNumLongP,fWeightDenLongP);
 fWeightRatLongP->Divide(fWeightNumLongP,fWeightDenLongP,longPscale);

 Double_t longNscale = Scale(fWeightNumLongN,fWeightDenLongN);
 fWeightRatLongN->Divide(fWeightNumLongN,fWeightDenLongN,longNscale);

 fWeightRatLongNOverP->Divide(fWeightRatLongN,fWeightRatLongP);

 fWeightNumOutP->Write();
 fWeightDenOutP->Write();
 fWeightRatOutP->Write();
 fWeightNumOutN->Write();
 fWeightDenOutN->Write();
 fWeightRatOutN->Write();
 fWeightRatOut->Write();
 fWeightRatOutNOverP->Write();

 fWeightNumSideP->Write();
 fWeightDenSideP->Write();
 fWeightRatSideP->Write();
 fWeightNumSideN->Write();
 fWeightDenSideN->Write();
 fWeightRatSideN->Write();
 fWeightRatSide->Write();
 fWeightRatSideNOverP->Write();

 fWeightNumLongP->Write();
 fWeightDenLongP->Write();
 fWeightRatLongP->Write();
 fWeightNumLongN->Write();
 fWeightDenLongN->Write();
 fWeightRatLongN->Write();
 fWeightRatLong->Write();
 fWeightRatLongNOverP->Write();

}

