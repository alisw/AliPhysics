#include "AliHBTNonId3DCorrFctn.h"
#include <TH1.h>
#include <Riostream.h>

///////////////////////////////////////////////////////
//                                                   //
// AliHBTNonId3DCorrFctn.h                           //
//                                                   //
// Class for calculating 3D non-id correlation       //
// functions                                         //
//                                                   //
///////////////////////////////////////////////////////

ClassImp(AliHBTNonId3DCorrFctn)

AliHBTNonId3DCorrFctn::AliHBTNonId3DCorrFctn(const char* name, const char* title):
 AliHBTOnePairFctn1D(name,title),
 fNumOutP(0x0),
 fDenOutP(0x0),
 fRatOutP(0x0),
 fNumOutN(0x0),
 fDenOutN(0x0),
 fRatOutN(0x0),
 fRatOut(0x0),
 fRatOutNOverP(0x0),
 fNumSideP(0x0),
 fDenSideP(0x0),
 fRatSideP(0x0),
 fNumSideN(0x0),
 fDenSideN(0x0),
 fRatSideN(0x0),
 fRatSide(0x0),
 fRatSideNOverP(0x0),
 fNumLongP(0x0),
 fDenLongP(0x0),
 fRatLongP(0x0),
 fNumLongN(0x0),
 fDenLongN(0x0),
 fRatLongN(0x0),
 fRatLong(0x0),
 fRatLongNOverP(0x0)
{
//ctor
}
/******************************************************************/
AliHBTNonId3DCorrFctn::AliHBTNonId3DCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval):
 AliHBTOnePairFctn1D(name,title,nbins,maxXval,minXval),
 fNumOutP(0x0),
 fDenOutP(0x0),
 fRatOutP(0x0),
 fNumOutN(0x0),
 fDenOutN(0x0),
 fRatOutN(0x0),
 fRatOut(0x0),
 fRatOutNOverP(0x0),
 fNumSideP(0x0),
 fDenSideP(0x0),
 fRatSideP(0x0),
 fNumSideN(0x0),
 fDenSideN(0x0),
 fRatSideN(0x0),
 fRatSide(0x0),
 fRatSideNOverP(0x0),
 fNumLongP(0x0),
 fDenLongP(0x0),
 fRatLongP(0x0),
 fNumLongN(0x0),
 fDenLongN(0x0),
 fRatLongN(0x0),
 fRatLong(0x0),
 fRatLongNOverP(0x0)
{
//ctor
}

/******************************************************************/
AliHBTNonId3DCorrFctn::AliHBTNonId3DCorrFctn(const AliHBTNonId3DCorrFctn& in):
 AliHBTOnePairFctn1D(in),
 fNumOutP(0x0),
 fDenOutP(0x0),
 fRatOutP(0x0),
 fNumOutN(0x0),
 fDenOutN(0x0),
 fRatOutN(0x0),
 fRatOut(0x0),
 fRatOutNOverP(0x0),
 fNumSideP(0x0),
 fDenSideP(0x0),
 fRatSideP(0x0),
 fNumSideN(0x0),
 fDenSideN(0x0),
 fRatSideN(0x0),
 fRatSide(0x0),
 fRatSideNOverP(0x0),
 fNumLongP(0x0),
 fDenLongP(0x0),
 fRatLongP(0x0),
 fNumLongN(0x0),
 fDenLongN(0x0),
 fRatLongN(0x0),
 fRatLong(0x0),
 fRatLongNOverP(0x0)
{
//ctor
}

/******************************************************************/

AliHBTNonId3DCorrFctn::~AliHBTNonId3DCorrFctn()
{
 //dtor
 delete fNumOutP;
 delete fDenOutP;
 delete fRatOutP;
 delete fNumOutN;
 delete fDenOutN;
 delete fRatOutN;
 delete fRatOut;
 delete fRatOutNOverP;
 delete fNumSideP;
 delete fDenSideP;
 delete fRatSideP;
 delete fNumSideN;
 delete fDenSideN;
 delete fRatSideN;
 delete fRatSide;
 delete fRatSideNOverP;
 delete fNumLongP;
 delete fDenLongP;
 delete fRatLongP;
 delete fNumLongN;
 delete fDenLongN;
 delete fRatLongN;
 delete fRatLong;
 delete fRatLongNOverP;
 
}

/******************************************************************/
void AliHBTNonId3DCorrFctn::WriteFunction()
{
    
 Double_t outPscale = Scale(fNumOutP,fDenOutP);
 cout<<"outPscale = "<<outPscale<<endl;
 fRatOutP->Divide(fNumOutP,fDenOutP,outPscale);

 Double_t outNscale = Scale(fNumOutN,fDenOutN);
 cout<<"outNscale = "<<outNscale<<endl;
 fRatOutN->Divide(fNumOutN,fDenOutN,outNscale);

 fRatOutNOverP->Divide(fRatOutN,fRatOutP);


 Double_t sidePscale = Scale(fNumSideP,fDenSideP);
 fRatSideP->Divide(fNumSideP,fDenSideP,sidePscale);

 Double_t sideNscale = Scale(fNumSideN,fDenSideN);
 fRatSideN->Divide(fNumSideN,fDenSideN,sideNscale);

 fRatSideNOverP->Divide(fRatSideN,fRatSideP);


 Double_t longPscale = Scale(fNumLongP,fDenLongP);
 fRatLongP->Divide(fNumLongP,fDenLongP,longPscale);

 Double_t longNscale = Scale(fNumLongN,fDenLongN);
 fRatLongN->Divide(fNumLongN,fDenLongN,longNscale);

 fRatLongNOverP->Divide(fRatLongN,fRatLongP);

 fNumOutP->Write();
 fDenOutP->Write();
 fRatOutP->Write();
 fNumOutN->Write();
 fDenOutN->Write();
 fRatOutN->Write();
 fRatOut->Write();
 fRatOutNOverP->Write();

 fNumSideP->Write();
 fDenSideP->Write();
 fRatSideP->Write();
 fNumSideN->Write();
 fDenSideN->Write();
 fRatSideN->Write();
 fRatSide->Write();
 fRatSideNOverP->Write();

 fNumLongP->Write();
 fDenLongP->Write();
 fRatLongP->Write();
 fNumLongN->Write();
 fDenLongN->Write();
 fRatLongN->Write();
 fRatLong->Write();
 fRatLongNOverP->Write();

}

//-------------------------------------
void AliHBTNonId3DCorrFctn::ProcessSameEventParticles(AliHBTPair* pair)
{
    //Fills the numerator using pair from the same event
    pair = CheckPair(pair);
    if(pair == 0x0) return;
 
    double tKStar = TMath::Abs(pair->GetKStar());
    double tKStarOut = pair->GetKStarOut();
    double tKStarSide = pair->GetKStarSide();
    double tKStarLong = pair->GetKStarLong();

    if(tKStarOut>0.)
    {
	fNumOutP->Fill(tKStar);
    }
    else
    {
	fNumOutN->Fill(tKStar);
    }

    if(tKStarSide>0.)
    {
	fNumSideP->Fill(tKStar);
    }
    else
    {
	fNumSideN->Fill(tKStar);
    }

    if(tKStarLong>0.)
    {
	fNumLongP->Fill(tKStar);
    }
    else
    {
	fNumLongN->Fill(tKStar);
    }
}

/****************************************************************/
void AliHBTNonId3DCorrFctn::Init()
{
    AliHBTOnePairFctn1D::Init();
 fNumOutP->Reset();
 fDenOutP->Reset();
 fRatOutP->Reset();
 fNumOutN->Reset();
 fDenOutN->Reset();
 fRatOutN->Reset();
 fRatOut->Reset();
 fRatOutNOverP->Reset();

 fNumSideP->Reset();
 fDenSideP->Reset();
 fRatSideP->Reset();
 fNumSideN->Reset();
 fDenSideN->Reset();
 fRatSideN->Reset();
 fRatSide->Reset();
 fRatSideNOverP->Reset();

 fNumLongP->Reset();
 fDenLongP->Reset();
 fRatLongP->Reset();
 fNumLongN->Reset();
 fDenLongN->Reset();
 fRatLongN->Reset();
 fRatLong->Reset();
 fRatLongNOverP->Reset();
}
/****************************************************************/


void AliHBTNonId3DCorrFctn::ProcessDiffEventParticles(AliHBTPair* pair)
{

    double tKStar = TMath::Abs(pair->GetKStar());
    double tKStarOut = pair->GetKStarOut();
    double tKStarSide = pair->GetKStarSide();
    double tKStarLong = pair->GetKStarLong();

    if(tKStarOut>0.)
    {
	fDenOutP->Fill(tKStar);
    }
    else
    {
	fDenOutN->Fill(tKStar);
    }

    if(tKStarSide>0.)
    {
	fDenSideP->Fill(tKStar);
    }
    else
    {
	fDenSideN->Fill(tKStar);
    }

    if(tKStarLong>0.)
    {
	fDenLongP->Fill(tKStar);
    }
    else
    {
	fDenLongN->Fill(tKStar);
    }

}


/******************************************************************/

void AliHBTNonId3DCorrFctn::BuildHistos(Int_t nbins, Float_t max, Float_t min)
{
    
    AliHBTFunction1D::BuildHistos(nbins,max,min);
    
    TString nameNumOutP = "NumOutP";
    TString nameDenOutP = "DenOutP";
    TString nameRatOutP = "RatOutP";
    TString nameNumOutN = "NumOutN";
    TString nameDenOutN = "DenOutN";
    TString nameRatOutN = "RatOutN";
    TString nameRatOut =  "RatOut";
    TString nameRatOutNOverP = "RatOutNOverP";
    TString nameNumSideP = "NumSideP";
    TString nameDenSideP = "DenSideP";
    TString nameRatSideP = "RatSideP";
    TString nameNumSideN = "NumSideN";
    TString nameDenSideN = "DenSideN";
    TString nameRatSideN = "RatSideN";
    TString nameRatSide = "RatSide";
    TString nameRatSideNOverP = "RatSideNOverP";
    TString nameNumLongP = "NumLongP";
    TString nameDenLongP = "DenLongP";
    TString nameRatLongP = "RatLongP";
    TString nameNumLongN = "NumLongN";
    TString nameDenLongN = "DenLongN";
    TString nameRatLongN = "RatLongN";
    TString nameRatLong = "RatLong";
    TString nameRatLongNOverP = "RatLongNOverP";

    fNumOutP = new TH1D(nameNumOutP.Data(),nameNumOutP.Data(),nbins,min,max);
    fDenOutP = new TH1D(nameDenOutP.Data(),nameDenOutP.Data(),nbins,min,max);
    fRatOutP = new TH1D(nameRatOutP.Data(),nameRatOutP.Data(),nbins,min,max);

    fNumOutN = new TH1D(nameNumOutN.Data(),nameNumOutN.Data(),nbins,min,max);
    fDenOutN = new TH1D(nameDenOutN.Data(),nameDenOutN.Data(),nbins,min,max);
    fRatOutN = new TH1D(nameRatOutN.Data(),nameRatOutN.Data(),nbins,min,max);

    fRatOut = new TH1D(nameRatOut.Data(),nameRatOut.Data(),nbins,min,max);
    fRatOutNOverP = new TH1D(nameRatOutNOverP.Data(),nameRatOutNOverP.Data(),nbins,min,max);

    fNumSideP = new TH1D(nameNumSideP.Data(),nameNumSideP.Data(),nbins,min,max);
    fDenSideP = new TH1D(nameDenSideP.Data(),nameDenSideP.Data(),nbins,min,max);
    fRatSideP = new TH1D(nameRatSideP.Data(),nameRatSideP.Data(),nbins,min,max);

    fNumSideN = new TH1D(nameNumSideN.Data(),nameNumSideN.Data(),nbins,min,max);
    fDenSideN = new TH1D(nameDenSideN.Data(),nameDenSideN.Data(),nbins,min,max);
    fRatSideN = new TH1D(nameRatSideN.Data(),nameRatSideN.Data(),nbins,min,max);

    fRatSide = new TH1D(nameRatSide.Data(),nameRatSide.Data(),nbins,min,max);
    fRatSideNOverP = new TH1D(nameRatSideNOverP.Data(),nameRatSideNOverP.Data(),nbins,min,max);

    fNumLongP = new TH1D(nameNumLongP.Data(),nameNumLongP.Data(),nbins,min,max);
    fDenLongP = new TH1D(nameDenLongP.Data(),nameDenLongP.Data(),nbins,min,max);
    fRatLongP = new TH1D(nameRatLongP.Data(),nameRatLongP.Data(),nbins,min,max);

    fNumLongN = new TH1D(nameNumLongN.Data(),nameNumLongN.Data(),nbins,min,max);
    fDenLongN = new TH1D(nameDenLongN.Data(),nameDenLongN.Data(),nbins,min,max);
    fRatLongN = new TH1D(nameRatLongN.Data(),nameRatLongN.Data(),nbins,min,max);

    fRatLong = new TH1D(nameRatLong.Data(),nameRatLong.Data(),nbins,min,max);
    fRatLongNOverP = new TH1D(nameRatLongNOverP.Data(),nameRatLongNOverP.Data(),nbins,min,max);


    fNumOutP->Sumw2();
    fDenOutP->Sumw2();
    fRatOutP->Sumw2();
    fNumOutN->Sumw2();
    fDenOutN->Sumw2();
    fRatOutN->Sumw2();
    fRatOut->Sumw2();
    fRatOutNOverP->Sumw2();
    fNumSideP->Sumw2();
    fDenSideP->Sumw2();
    fRatSideP->Sumw2();
    fNumSideN->Sumw2();
    fDenSideN->Sumw2();
    fRatSideN->Sumw2();
    fRatSide->Sumw2();
    fRatSideNOverP->Sumw2();
    fNumLongP->Sumw2();
    fDenLongP->Sumw2();
    fRatLongP->Sumw2();
    fNumLongN->Sumw2();
    fDenLongN->Sumw2();
    fRatLongN->Sumw2();
    fRatLong->Sumw2();
    fRatLongNOverP->Sumw2();
}

 TH1* AliHBTNonId3DCorrFctn::GetResult()
 {
     return fRatOutN;
 }
