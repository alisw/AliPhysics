#include "AliHBTashbtCorrFctn.h"
#include <TH1.h>
#include <Riostream.h>

///////////////////////////////////////////////////////
//                                                   //
// AliHBTashbtCorrFctn.h                             //
//                                                   //
// Class for calculating 3D ashbt correlation        //
// functions                                         //
//                                                   //
///////////////////////////////////////////////////////

ClassImp(AliHBTashbtCorrFctn)

AliHBTashbtCorrFctn::AliHBTashbtCorrFctn(const char* name, const char* title):
 AliHBTOnePairFctn1D(name,title),
    fNumOut1(0x0),
    fNumOut2(0x0),
    fNumOut3(0x0),
    fNumOut4(0x0),
    fNumOut5(0x0),
    fNumOut6(0x0),
    fNumOut7(0x0),
    fNumOut8(0x0),

    fDenOut1(0x0),
    fDenOut2(0x0),
    fDenOut3(0x0),
    fDenOut4(0x0),
    fDenOut5(0x0),
    fDenOut6(0x0),
    fDenOut7(0x0),
    fDenOut8(0x0),
   
    fRatOut1(0x0),
    fRatOut2(0x0),
    fRatOut3(0x0),
    fRatOut4(0x0),
    fRatOut5(0x0),
    fRatOut6(0x0),
    fRatOut7(0x0),
    fRatOut8(0x0),
   

    fNumSide1(0x0),
    fNumSide2(0x0),
    fNumSide3(0x0),
    fNumSide4(0x0),
    fNumSide5(0x0),
    fNumSide6(0x0),
    fNumSide7(0x0),
    fNumSide8(0x0),

    fDenSide1(0x0),
    fDenSide2(0x0),
    fDenSide3(0x0),
    fDenSide4(0x0),
    fDenSide5(0x0),
    fDenSide6(0x0),
    fDenSide7(0x0),
    fDenSide8(0x0),
   
    fRatSide1(0x0),
    fRatSide2(0x0),
    fRatSide3(0x0),
    fRatSide4(0x0),
    fRatSide5(0x0),
    fRatSide6(0x0),
    fRatSide7(0x0),
    fRatSide8(0x0),
   
    fNumLong1(0x0),
    fNumLong2(0x0),
    fNumLong3(0x0),
    fNumLong4(0x0),
    fNumLong5(0x0),
    fNumLong6(0x0),
    fNumLong7(0x0),
    fNumLong8(0x0),

    fDenLong1(0x0),
    fDenLong2(0x0),
    fDenLong3(0x0),
    fDenLong4(0x0),
    fDenLong5(0x0),
    fDenLong6(0x0),
    fDenLong7(0x0),
    fDenLong8(0x0),
   
    fRatLong1(0x0),
    fRatLong2(0x0),
    fRatLong3(0x0),
    fRatLong4(0x0),
    fRatLong5(0x0),
    fRatLong6(0x0),
    fRatLong7(0x0),
    fRatLong8(0x0)


{
//ctor
}
/******************************************************************/
AliHBTashbtCorrFctn::AliHBTashbtCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval):
 AliHBTOnePairFctn1D(name,title,nbins,maxXval,minXval),

    fNumOut1(0x0),
    fNumOut2(0x0),
    fNumOut3(0x0),
    fNumOut4(0x0),
    fNumOut5(0x0),
    fNumOut6(0x0),
    fNumOut7(0x0),
    fNumOut8(0x0),

    fDenOut1(0x0),
    fDenOut2(0x0),
    fDenOut3(0x0),
    fDenOut4(0x0),
    fDenOut5(0x0),
    fDenOut6(0x0),
    fDenOut7(0x0),
    fDenOut8(0x0),
   
    fRatOut1(0x0),
    fRatOut2(0x0),
    fRatOut3(0x0),
    fRatOut4(0x0),
    fRatOut5(0x0),
    fRatOut6(0x0),
    fRatOut7(0x0),
    fRatOut8(0x0),
   

    fNumSide1(0x0),
    fNumSide2(0x0),
    fNumSide3(0x0),
    fNumSide4(0x0),
    fNumSide5(0x0),
    fNumSide6(0x0),
    fNumSide7(0x0),
    fNumSide8(0x0),

    fDenSide1(0x0),
    fDenSide2(0x0),
    fDenSide3(0x0),
    fDenSide4(0x0),
    fDenSide5(0x0),
    fDenSide6(0x0),
    fDenSide7(0x0),
    fDenSide8(0x0),
   
    fRatSide1(0x0),
    fRatSide2(0x0),
    fRatSide3(0x0),
    fRatSide4(0x0),
    fRatSide5(0x0),
    fRatSide6(0x0),
    fRatSide7(0x0),
    fRatSide8(0x0),
   
    fNumLong1(0x0),
    fNumLong2(0x0),
    fNumLong3(0x0),
    fNumLong4(0x0),
    fNumLong5(0x0),
    fNumLong6(0x0),
    fNumLong7(0x0),
    fNumLong8(0x0),

    fDenLong1(0x0),
    fDenLong2(0x0),
    fDenLong3(0x0),
    fDenLong4(0x0),
    fDenLong5(0x0),
    fDenLong6(0x0),
    fDenLong7(0x0),
    fDenLong8(0x0),
   
    fRatLong1(0x0),
    fRatLong2(0x0),
    fRatLong3(0x0),
    fRatLong4(0x0),
    fRatLong5(0x0),
    fRatLong6(0x0),
    fRatLong7(0x0),
    fRatLong8(0x0)


{
BuildHistos(nbins, maxXval, minXval);
//ctor
}

/******************************************************************/
AliHBTashbtCorrFctn::AliHBTashbtCorrFctn(const AliHBTashbtCorrFctn& in):
 AliHBTOnePairFctn1D(in),
    fNumOut1(0x0),
    fNumOut2(0x0),
    fNumOut3(0x0),
    fNumOut4(0x0),
    fNumOut5(0x0),
    fNumOut6(0x0),
    fNumOut7(0x0),
    fNumOut8(0x0),

    fDenOut1(0x0),
    fDenOut2(0x0),
    fDenOut3(0x0),
    fDenOut4(0x0),
    fDenOut5(0x0),
    fDenOut6(0x0),
    fDenOut7(0x0),
    fDenOut8(0x0),
   
    fRatOut1(0x0),
    fRatOut2(0x0),
    fRatOut3(0x0),
    fRatOut4(0x0),
    fRatOut5(0x0),
    fRatOut6(0x0),
    fRatOut7(0x0),
    fRatOut8(0x0),
   

    fNumSide1(0x0),
    fNumSide2(0x0),
    fNumSide3(0x0),
    fNumSide4(0x0),
    fNumSide5(0x0),
    fNumSide6(0x0),
    fNumSide7(0x0),
    fNumSide8(0x0),

    fDenSide1(0x0),
    fDenSide2(0x0),
    fDenSide3(0x0),
    fDenSide4(0x0),
    fDenSide5(0x0),
    fDenSide6(0x0),
    fDenSide7(0x0),
    fDenSide8(0x0),
   
    fRatSide1(0x0),
    fRatSide2(0x0),
    fRatSide3(0x0),
    fRatSide4(0x0),
    fRatSide5(0x0),
    fRatSide6(0x0),
    fRatSide7(0x0),
    fRatSide8(0x0),
   
    fNumLong1(0x0),
    fNumLong2(0x0),
    fNumLong3(0x0),
    fNumLong4(0x0),
    fNumLong5(0x0),
    fNumLong6(0x0),
    fNumLong7(0x0),
    fNumLong8(0x0),

    fDenLong1(0x0),
    fDenLong2(0x0),
    fDenLong3(0x0),
    fDenLong4(0x0),
    fDenLong5(0x0),
    fDenLong6(0x0),
    fDenLong7(0x0),
    fDenLong8(0x0),
   
    fRatLong1(0x0),
    fRatLong2(0x0),
    fRatLong3(0x0),
    fRatLong4(0x0),
    fRatLong5(0x0),
    fRatLong6(0x0),
    fRatLong7(0x0),
    fRatLong8(0x0)
  

 {
//ctor
}

/******************************************************************/

AliHBTashbtCorrFctn::~AliHBTashbtCorrFctn()
{
 //dtor
    delete fNumOut1;
    delete fNumOut2;
    delete fNumOut3;
    delete fNumOut4;
    delete fNumOut5;
    delete fNumOut6;
    delete fNumOut7;
    delete fNumOut8;

    delete fDenOut1;
    delete fDenOut2;
    delete fDenOut3;
    delete fDenOut4;
    delete fDenOut5;
    delete fDenOut6;
    delete fDenOut7;
    delete fDenOut8;
   
    delete fRatOut1;
    delete fRatOut2;
    delete fRatOut3;
    delete fRatOut4;
    delete fRatOut5;
    delete fRatOut6;
    delete fRatOut7;
    delete fRatOut8;
   

    delete fNumSide1;
    delete fNumSide2;
    delete fNumSide3;
    delete fNumSide4;
    delete fNumSide5;
    delete fNumSide6;
    delete fNumSide7;
    delete fNumSide8;

    delete fDenSide1;
    delete fDenSide2;
    delete fDenSide3;
    delete fDenSide4;
    delete fDenSide5;
    delete fDenSide6;
    delete fDenSide7;
    delete fDenSide8;
   
    delete fRatSide1;
    delete fRatSide2;
    delete fRatSide3;
    delete fRatSide4;
    delete fRatSide5;
    delete fRatSide6;
    delete fRatSide7;
    delete fRatSide8;
   
    delete fNumLong1;
    delete fNumLong2;
    delete fNumLong3;
    delete fNumLong4;
    delete fNumLong5;
    delete fNumLong6;
    delete fNumLong7;
    delete fNumLong8;

    delete fDenLong1;
    delete fDenLong2;
    delete fDenLong3;
    delete fDenLong4;
    delete fDenLong5;
    delete fDenLong6;
    delete fDenLong7;
    delete fDenLong8;
   
    delete fRatLong1;
    delete fRatLong2;
    delete fRatLong3;
    delete fRatLong4;
    delete fRatLong5;
    delete fRatLong6;
    delete fRatLong7;
    delete fRatLong8; 

}

/******************************************************************/
void AliHBTashbtCorrFctn::WriteFunction()
{
//out    
   Double_t out1scale = Scale(fNumOut1,fDenOut1);
   cout <<"out1scale = "<<out1scale<<endl;
   fRatOut1->Divide(fNumOut1,fDenOut1,out1scale);
   
   Double_t out2scale = Scale(fNumOut2,fDenOut2);
   cout <<"out2scale = "<<out2scale<<endl;
   fRatOut2->Divide(fNumOut2,fDenOut2,out2scale);
   
   Double_t out3scale = Scale(fNumOut3,fDenOut3);
   cout <<"out3scale = "<<out3scale<<endl;
   fRatOut3->Divide(fNumOut3,fDenOut3,out3scale);
   
   Double_t out4scale = Scale(fNumOut4,fDenOut4);
   cout <<"out4scale = "<<out4scale<<endl;
   fRatOut4->Divide(fNumOut4,fDenOut4,out4scale);
   
   Double_t out5scale = Scale(fNumOut5,fDenOut5);
   cout <<"out5scale = "<<out5scale<<endl;
   fRatOut5->Divide(fNumOut5,fDenOut5,out5scale);
   
   Double_t out6scale = Scale(fNumOut6,fDenOut6);
   cout <<"out6scale = "<<out6scale<<endl;
   fRatOut6->Divide(fNumOut6,fDenOut6,out6scale);
   
   Double_t out7scale = Scale(fNumOut7,fDenOut7);
   cout <<"out7scale = "<<out7scale<<endl;
   fRatOut7->Divide(fNumOut7,fDenOut7,out7scale);
   
   Double_t out8scale = Scale(fNumOut8,fDenOut8);
   cout <<"out8scale = "<<out8scale<<endl;
   fRatOut8->Divide(fNumOut8,fDenOut8,out8scale);

   //side
   Double_t side1scale = Scale(fNumSide1,fDenSide1);
   cout <<"side1scale = "<<side1scale<<endl;
   fRatSide1->Divide(fNumSide1,fDenSide1,side1scale);
   
   Double_t side2scale = Scale(fNumSide2,fDenSide2);
   cout <<"side2scale = "<<side2scale<<endl;
   fRatSide2->Divide(fNumSide2,fDenSide2,side2scale);
   
   Double_t side3scale = Scale(fNumSide3,fDenSide3);
   cout <<"side3scale = "<<side3scale<<endl;
   fRatSide3->Divide(fNumSide3,fDenSide3,side3scale);
   
   Double_t side4scale = Scale(fNumSide4,fDenSide4);
   cout <<"side4scale = "<<side4scale<<endl;
   fRatSide4->Divide(fNumSide4,fDenSide4,side4scale);
   
   Double_t side5scale = Scale(fNumSide5,fDenSide5);
   cout <<"side5scale = "<<side5scale<<endl;
   fRatSide5->Divide(fNumSide5,fDenSide5,side5scale);
   
   Double_t side6scale = Scale(fNumSide6,fDenSide6);
   cout <<"side6scale = "<<side6scale<<endl;
   fRatSide6->Divide(fNumSide6,fDenSide6,side6scale);
   
   Double_t side7scale = Scale(fNumSide7,fDenSide7);
   cout <<"side7scale = "<<side7scale<<endl;
   fRatSide7->Divide(fNumSide7,fDenSide7,side7scale);
   
   Double_t side8scale = Scale(fNumSide8,fDenSide8);
   cout <<"side8scale = "<<side8scale<<endl;
   fRatSide8->Divide(fNumSide8,fDenSide8,side8scale);
   
//long
 
   Double_t long1scale = Scale(fNumLong1,fDenLong1);
   cout <<"long1scale = "<<long1scale<<endl;
   fRatLong1->Divide(fNumLong1,fDenLong1,long1scale);
   
   Double_t long2scale = Scale(fNumLong2,fDenLong2);
   cout <<"long2scale = "<<long2scale<<endl;
   fRatLong2->Divide(fNumLong2,fDenLong2,long2scale);
   
   Double_t long3scale = Scale(fNumLong3,fDenLong3);
   cout <<"long3scale = "<<long3scale<<endl;
   fRatLong3->Divide(fNumLong3,fDenLong3,long3scale);
   
   Double_t long4scale = Scale(fNumLong4,fDenLong4);
   cout <<"long4scale = "<<long4scale<<endl;
   fRatLong4->Divide(fNumLong4,fDenLong4,long4scale);
   
   Double_t long5scale = Scale(fNumLong5,fDenLong5);
   cout <<"long5scale = "<<long5scale<<endl;
   fRatLong5->Divide(fNumLong5,fDenLong5,long5scale);
   
   Double_t long6scale = Scale(fNumLong6,fDenLong6);
   cout <<"long6scale = "<<long6scale<<endl;
   fRatLong6->Divide(fNumLong6,fDenLong6,long6scale);
   
   Double_t long7scale = Scale(fNumLong7,fDenLong7);
   cout <<"long7scale = "<<long7scale<<endl;
   fRatLong7->Divide(fNumLong7,fDenLong7,long7scale);
   
   Double_t long8scale = Scale(fNumLong8,fDenLong8);
   cout <<"long8scale = "<<long8scale<<endl;
   fRatLong8->Divide(fNumLong8,fDenLong8,long8scale);

        fNumOut1->Write();
     fNumOut2->Write();
     fNumOut3->Write();
     fNumOut4->Write();
     fNumOut5->Write();
     fNumOut6->Write();
     fNumOut7->Write();
     fNumOut8->Write();

     fDenOut1->Write();
     fDenOut2->Write();
     fDenOut3->Write();
     fDenOut4->Write();
     fDenOut5->Write();
     fDenOut6->Write();
     fDenOut7->Write();
     fDenOut8->Write();
   
     fRatOut1->Write();
     fRatOut2->Write();
     fRatOut3->Write();
     fRatOut4->Write();
     fRatOut5->Write();
     fRatOut6->Write();
     fRatOut7->Write();
     fRatOut8->Write();
   

     fNumSide1->Write();
     fNumSide2->Write();
     fNumSide3->Write();
     fNumSide4->Write();
     fNumSide5->Write();
     fNumSide6->Write();
     fNumSide7->Write();
     fNumSide8->Write();

     fDenSide1->Write();
     fDenSide2->Write();
     fDenSide3->Write();
     fDenSide4->Write();
     fDenSide5->Write();
     fDenSide6->Write();
     fDenSide7->Write();
     fDenSide8->Write();
   
     fRatSide1->Write();
     fRatSide2->Write();
     fRatSide3->Write();
     fRatSide4->Write();
     fRatSide5->Write();
     fRatSide6->Write();
     fRatSide7->Write();
     fRatSide8->Write();
   
     fNumLong1->Write();
     fNumLong2->Write();
     fNumLong3->Write();
     fNumLong4->Write();
     fNumLong5->Write();
     fNumLong6->Write();
     fNumLong7->Write();
     fNumLong8->Write();

     fDenLong1->Write();
     fDenLong2->Write();
     fDenLong3->Write();
     fDenLong4->Write();
     fDenLong5->Write();
     fDenLong6->Write();
     fDenLong7->Write();
     fDenLong8->Write();
   
     fRatLong1->Write();
     fRatLong2->Write();
     fRatLong3->Write();
     fRatLong4->Write();
     fRatLong5->Write();
     fRatLong6->Write();
     fRatLong7->Write();
     fRatLong8->Write(); 

}

//-------------------------------------
void AliHBTashbtCorrFctn::ProcessSameEventParticles(AliHBTPair* pair)
{
    //Fills the numerator using pair from the same event
    pair = CheckPair(pair);
    if(pair == 0x0) return;
 
   Double_t rplane=0.;   //reaction plane angle - 2 B determined
   Double_t phi=(pair->Particle1()->Phi()+pair->Particle2()->Phi())/2.-rplane; //deltaphi bo nie mam nic innego pod reka
   phi=phi*360/(2*TMath::Pi());
   Double_t qout=pair->GetQOutLCMS();
   Double_t qside=pair->GetQSideLCMS();
   Double_t qlong=pair->GetQLongLCMS();
   
   if((phi>=0.) && (phi <45.))
     {
	fNumOut1->Fill(qout);
	fNumSide1->Fill(qside);
	fNumLong1->Fill(qlong);
     }
   else if((phi>=45.) && (phi <90.))
     {
	fNumOut2->Fill(qout);
	fNumSide2->Fill(qside);
	fNumLong2->Fill(qlong);
     }
   else if((phi>=90.) && (phi <135.))
     {
	fNumOut3->Fill(qout);
	fNumSide3->Fill(qside);
	fNumLong3->Fill(qlong);
     }
   else if((phi>=135.) && (phi <180.))
     {
	fNumOut4->Fill(qout);
	fNumSide4->Fill(qside);
	fNumLong4->Fill(qlong);
     }
   else if((phi>=180.) && (phi <225.))
     {
	fNumOut5->Fill(qout);
	fNumSide5->Fill(qside);
	fNumLong5->Fill(qlong);
     }
   else if((phi>=225.) && (phi <270.))
     {
	fNumOut6->Fill(qout);
	fNumSide6->Fill(qside);
	fNumLong6->Fill(qlong);
     }
   else if((phi>=270.) && (phi <315.))
     {
	fNumOut7->Fill(qout);
	fNumSide7->Fill(qside);
	fNumLong7->Fill(qlong);
     }
   else if((phi>=315.) && (phi <360.))
     {
	fNumOut8->Fill(qout);
	fNumSide8->Fill(qside);
	fNumLong8->Fill(qlong);
     }

   
}

/****************************************************************/
void AliHBTashbtCorrFctn::Init()
{
     fNumOut1->Reset();
     fNumOut2->Reset();
     fNumOut3->Reset();
     fNumOut4->Reset();
     fNumOut5->Reset();
     fNumOut6->Reset();
     fNumOut7->Reset();
     fNumOut8->Reset();

     fDenOut1->Reset();
     fDenOut2->Reset();
     fDenOut3->Reset();
     fDenOut4->Reset();
     fDenOut5->Reset();
     fDenOut6->Reset();
     fDenOut7->Reset();
     fDenOut8->Reset();
   
     fRatOut1->Reset();
     fRatOut2->Reset();
     fRatOut3->Reset();
     fRatOut4->Reset();
     fRatOut5->Reset();
     fRatOut6->Reset();
     fRatOut7->Reset();
     fRatOut8->Reset();
   

     fNumSide1->Reset();
     fNumSide2->Reset();
     fNumSide3->Reset();
     fNumSide4->Reset();
     fNumSide5->Reset();
     fNumSide6->Reset();
     fNumSide7->Reset();
     fNumSide8->Reset();

     fDenSide1->Reset();
     fDenSide2->Reset();
     fDenSide3->Reset();
     fDenSide4->Reset();
     fDenSide5->Reset();
     fDenSide6->Reset();
     fDenSide7->Reset();
     fDenSide8->Reset();
   
     fRatSide1->Reset();
     fRatSide2->Reset();
     fRatSide3->Reset();
     fRatSide4->Reset();
     fRatSide5->Reset();
     fRatSide6->Reset();
     fRatSide7->Reset();
     fRatSide8->Reset();
   
     fNumLong1->Reset();
     fNumLong2->Reset();
     fNumLong3->Reset();
     fNumLong4->Reset();
     fNumLong5->Reset();
     fNumLong6->Reset();
     fNumLong7->Reset();
     fNumLong8->Reset();

     fDenLong1->Reset();
     fDenLong2->Reset();
     fDenLong3->Reset();
     fDenLong4->Reset();
     fDenLong5->Reset();
     fDenLong6->Reset();
     fDenLong7->Reset();
     fDenLong8->Reset();
   
     fRatLong1->Reset();
     fRatLong2->Reset();
     fRatLong3->Reset();
     fRatLong4->Reset();
     fRatLong5->Reset();
     fRatLong6->Reset();
     fRatLong7->Reset();
     fRatLong8->Reset(); 
   

 }
/****************************************************************/


void AliHBTashbtCorrFctn::ProcessDiffEventParticles(AliHBTPair* pair)
{

      Double_t rplane=0.;   //reaction plane angle - 2 B determined
   Double_t phi=(pair->Particle1()->Phi()+pair->Particle2()->Phi())/2.-rplane; //deltaphi bo nie mam nic innego pod reka
   phi=phi*360/(2*TMath::Pi());
   Double_t qout=pair->GetQOutLCMS();
   Double_t qside=pair->GetQSideLCMS();
   Double_t qlong=pair->GetQLongLCMS();
   
   if((phi>=0.) && (phi <45.))
     {
	fDenOut1->Fill(qout);
	fDenSide1->Fill(qside);
	fDenLong1->Fill(qlong);
     }
   else if((phi>=45.) && (phi <90.))
     {
	fDenOut2->Fill(qout);
	fDenSide2->Fill(qside);
	fDenLong2->Fill(qlong);
     }
   else if((phi>=90.) && (phi <135.))
     {
	fDenOut3->Fill(qout);
	fDenSide3->Fill(qside);
	fDenLong3->Fill(qlong);
     }
   else if((phi>=135.) && (phi <180.))
     {
	fDenOut4->Fill(qout);
	fDenSide4->Fill(qside);
	fDenLong4->Fill(qlong);
     }
   else if((phi>=180.) && (phi <225.))
     {
	fDenOut5->Fill(qout);
	fDenSide5->Fill(qside);
	fDenLong5->Fill(qlong);
     }
   else if((phi>=225.) && (phi <270.))
     {
	fDenOut6->Fill(qout);
	fDenSide6->Fill(qside);
	fDenLong6->Fill(qlong);
     }
   else if((phi>=270.) && (phi <315.))
     {
	fDenOut7->Fill(qout);
	fDenSide7->Fill(qside);
	fDenLong7->Fill(qlong);
     }
   else if((phi>=315.) && (phi <360.))
     {
	fDenOut8->Fill(qout);
	fDenSide8->Fill(qside);
	fDenLong8->Fill(qlong);
     }

}


/******************************************************************/

void AliHBTashbtCorrFctn::BuildHistos(Int_t nbins, Float_t max, Float_t min)
{
    
   
   TString nameNumOut1 = "NumOut0-45deg";
   TString nameNumOut2 = "NumOut45-90deg";
   TString nameNumOut3 = "NumOut90-135deg";
   TString nameNumOut4 = "NumOut135-180deg";
   TString nameNumOut5 = "NumOut180-225deg";
   TString nameNumOut6 = "NumOut225-270deg";
   TString nameNumOut7 = "NumOut270-315deg";
   TString nameNumOut8 = "NumOut315-360deg";
   
   
   TString nameDenOut1 = "DenOut0-45deg";
   TString nameDenOut2 = "DenOut45-90deg";
   TString nameDenOut3 = "DenOut90-135deg";
   TString nameDenOut4 = "DenOut135-180deg";
   TString nameDenOut5 = "DenOut180-225deg";
   TString nameDenOut6 = "DenOut225-270deg";
   TString nameDenOut7 = "DenOut270-315deg";
   TString nameDenOut8 = "DenOut315-360deg";
   
   
   TString nameRatOut1 = "RatOut0-45deg";
   TString nameRatOut2 = "RatOut45-90deg";
   TString nameRatOut3 = "RatOut90-135deg";
   TString nameRatOut4 = "RatOut135-180deg";
   TString nameRatOut5 = "RatOut180-225deg";
   TString nameRatOut6 = "RatOut225-270deg";
   TString nameRatOut7 = "RatOut270-315deg";
   TString nameRatOut8 = "RatOut315-360deg";
   
   TString nameNumSide1 = "NumSide0-45deg";
   TString nameNumSide2 = "NumSide45-90deg";
   TString nameNumSide3 = "NumSide90-135deg";
   TString nameNumSide4 = "NumSide135-180deg";
   TString nameNumSide5 = "NumSide180-225deg";
   TString nameNumSide6 = "NumSide225-270deg";
   TString nameNumSide7 = "NumSide270-315deg";
   TString nameNumSide8 = "NumSide315-360deg";
   
   
   TString nameDenSide1 = "DenSide0-45deg";
   TString nameDenSide2 = "DenSide45-90deg";
   TString nameDenSide3 = "DenSide90-135deg";
   TString nameDenSide4 = "DenSide135-180deg";
   TString nameDenSide5 = "DenSide180-225deg";
   TString nameDenSide6 = "DenSide225-270deg";
   TString nameDenSide7 = "DenSide270-315deg";
   TString nameDenSide8 = "DenSide315-360deg";
   
   
   TString nameRatSide1 = "RatSide0-45deg";
   TString nameRatSide2 = "RatSide45-90deg";
   TString nameRatSide3 = "RatSide90-135deg";
   TString nameRatSide4 = "RatSide135-180deg";
   TString nameRatSide5 = "RatSide180-225deg";
   TString nameRatSide6 = "RatSide225-270deg";
   TString nameRatSide7 = "RatSide270-315deg";
   TString nameRatSide8 = "RatSide315-360deg";
   
   TString nameNumLong1 = "NumLong0-45deg";
   TString nameNumLong2 = "NumLong45-90deg";
   TString nameNumLong3 = "NumLong90-135deg";
   TString nameNumLong4 = "NumLong135-180deg";
   TString nameNumLong5 = "NumLong180-225deg";
   TString nameNumLong6 = "NumLong225-270deg";
   TString nameNumLong7 = "NumLong270-315deg";
   TString nameNumLong8 = "NumLong315-360deg";
   
   
   TString nameDenLong1 = "DenLong0-45deg";
   TString nameDenLong2 = "DenLong45-90deg";
   TString nameDenLong3 = "DenLong90-135deg";
   TString nameDenLong4 = "DenLong135-180deg";
   TString nameDenLong5 = "DenLong180-225deg";
   TString nameDenLong6 = "DenLong225-270deg";
   TString nameDenLong7 = "DenLong270-315deg";
   TString nameDenLong8 = "DenLong315-360deg";
   
   
   TString nameRatLong1 = "RatLong0-45deg";
   TString nameRatLong2 = "RatLong45-90deg";
   TString nameRatLong3 = "RatLong90-135deg";
   TString nameRatLong4 = "RatLong135-180deg";
   TString nameRatLong5 = "RatLong180-225deg";
   TString nameRatLong6 = "RatLong225-270deg";
   TString nameRatLong7 = "RatLong270-315deg";
   TString nameRatLong8 = "RatLong315-360deg";

   fNumOut1 = new TH1D(nameNumOut1.Data(),nameNumOut1.Data(),nbins,min,max);
   fNumOut2 = new TH1D(nameNumOut2.Data(),nameNumOut2.Data(),nbins,min,max);
   fNumOut3 = new TH1D(nameNumOut3.Data(),nameNumOut3.Data(),nbins,min,max);
   fNumOut4 = new TH1D(nameNumOut4.Data(),nameNumOut4.Data(),nbins,min,max);
   fNumOut5 = new TH1D(nameNumOut5.Data(),nameNumOut5.Data(),nbins,min,max);
   fNumOut6 = new TH1D(nameNumOut6.Data(),nameNumOut6.Data(),nbins,min,max);
   fNumOut7 = new TH1D(nameNumOut7.Data(),nameNumOut7.Data(),nbins,min,max);
   fNumOut8 = new TH1D(nameNumOut8.Data(),nameNumOut8.Data(),nbins,min,max);
   
   fDenOut1 = new TH1D(nameDenOut1.Data(),nameDenOut1.Data(),nbins,min,max);
   fDenOut2 = new TH1D(nameDenOut2.Data(),nameDenOut2.Data(),nbins,min,max);
   fDenOut3 = new TH1D(nameDenOut3.Data(),nameDenOut3.Data(),nbins,min,max);
   fDenOut4 = new TH1D(nameDenOut4.Data(),nameDenOut4.Data(),nbins,min,max);
   fDenOut5 = new TH1D(nameDenOut5.Data(),nameDenOut5.Data(),nbins,min,max);
   fDenOut6 = new TH1D(nameDenOut6.Data(),nameDenOut6.Data(),nbins,min,max);
   fDenOut7 = new TH1D(nameDenOut7.Data(),nameDenOut7.Data(),nbins,min,max);
   fDenOut8 = new TH1D(nameDenOut8.Data(),nameDenOut8.Data(),nbins,min,max);
   
   fRatOut1 = new TH1D(nameRatOut1.Data(),nameRatOut1.Data(),nbins,min,max);
   fRatOut2 = new TH1D(nameRatOut2.Data(),nameRatOut2.Data(),nbins,min,max);
   fRatOut3 = new TH1D(nameRatOut3.Data(),nameRatOut3.Data(),nbins,min,max);
   fRatOut4 = new TH1D(nameRatOut4.Data(),nameRatOut4.Data(),nbins,min,max);
   fRatOut5 = new TH1D(nameRatOut5.Data(),nameRatOut5.Data(),nbins,min,max);
   fRatOut6 = new TH1D(nameRatOut6.Data(),nameRatOut6.Data(),nbins,min,max);
   fRatOut7 = new TH1D(nameRatOut7.Data(),nameRatOut7.Data(),nbins,min,max);
   fRatOut8 = new TH1D(nameRatOut8.Data(),nameRatOut8.Data(),nbins,min,max);
   
   fNumSide1 = new TH1D(nameNumSide1.Data(),nameNumSide1.Data(),nbins,min,max);
   fNumSide2 = new TH1D(nameNumSide2.Data(),nameNumSide2.Data(),nbins,min,max);
   fNumSide3 = new TH1D(nameNumSide3.Data(),nameNumSide3.Data(),nbins,min,max);
   fNumSide4 = new TH1D(nameNumSide4.Data(),nameNumSide4.Data(),nbins,min,max);
   fNumSide5 = new TH1D(nameNumSide5.Data(),nameNumSide5.Data(),nbins,min,max);
   fNumSide6 = new TH1D(nameNumSide6.Data(),nameNumSide6.Data(),nbins,min,max);
   fNumSide7 = new TH1D(nameNumSide7.Data(),nameNumSide7.Data(),nbins,min,max);
   fNumSide8 = new TH1D(nameNumSide8.Data(),nameNumSide8.Data(),nbins,min,max);
   
   fDenSide1 = new TH1D(nameDenSide1.Data(),nameDenSide1.Data(),nbins,min,max);
   fDenSide2 = new TH1D(nameDenSide2.Data(),nameDenSide2.Data(),nbins,min,max);
   fDenSide3 = new TH1D(nameDenSide3.Data(),nameDenSide3.Data(),nbins,min,max);
   fDenSide4 = new TH1D(nameDenSide4.Data(),nameDenSide4.Data(),nbins,min,max);
   fDenSide5 = new TH1D(nameDenSide5.Data(),nameDenSide5.Data(),nbins,min,max);
   fDenSide6 = new TH1D(nameDenSide6.Data(),nameDenSide6.Data(),nbins,min,max);
   fDenSide7 = new TH1D(nameDenSide7.Data(),nameDenSide7.Data(),nbins,min,max);
   fDenSide8 = new TH1D(nameDenSide8.Data(),nameDenSide8.Data(),nbins,min,max);
   
   fRatSide1 = new TH1D(nameRatSide1.Data(),nameRatSide1.Data(),nbins,min,max);
   fRatSide2 = new TH1D(nameRatSide2.Data(),nameRatSide2.Data(),nbins,min,max);
   fRatSide3 = new TH1D(nameRatSide3.Data(),nameRatSide3.Data(),nbins,min,max);
   fRatSide4 = new TH1D(nameRatSide4.Data(),nameRatSide4.Data(),nbins,min,max);
   fRatSide5 = new TH1D(nameRatSide5.Data(),nameRatSide5.Data(),nbins,min,max);
   fRatSide6 = new TH1D(nameRatSide6.Data(),nameRatSide6.Data(),nbins,min,max);
   fRatSide7 = new TH1D(nameRatSide7.Data(),nameRatSide7.Data(),nbins,min,max);
   fRatSide8 = new TH1D(nameRatSide8.Data(),nameRatSide8.Data(),nbins,min,max);
   
   fNumLong1 = new TH1D(nameNumLong1.Data(),nameNumLong1.Data(),nbins,min,max);
   fNumLong2 = new TH1D(nameNumLong2.Data(),nameNumLong2.Data(),nbins,min,max);
   fNumLong3 = new TH1D(nameNumLong3.Data(),nameNumLong3.Data(),nbins,min,max);
   fNumLong4 = new TH1D(nameNumLong4.Data(),nameNumLong4.Data(),nbins,min,max);
   fNumLong5 = new TH1D(nameNumLong5.Data(),nameNumLong5.Data(),nbins,min,max);
   fNumLong6 = new TH1D(nameNumLong6.Data(),nameNumLong6.Data(),nbins,min,max);
   fNumLong7 = new TH1D(nameNumLong7.Data(),nameNumLong7.Data(),nbins,min,max);
   fNumLong8 = new TH1D(nameNumLong8.Data(),nameNumLong8.Data(),nbins,min,max);
   
   fDenLong1 = new TH1D(nameDenLong1.Data(),nameDenLong1.Data(),nbins,min,max);
   fDenLong2 = new TH1D(nameDenLong2.Data(),nameDenLong2.Data(),nbins,min,max);
   fDenLong3 = new TH1D(nameDenLong3.Data(),nameDenLong3.Data(),nbins,min,max);
   fDenLong4 = new TH1D(nameDenLong4.Data(),nameDenLong4.Data(),nbins,min,max);
   fDenLong5 = new TH1D(nameDenLong5.Data(),nameDenLong5.Data(),nbins,min,max);
   fDenLong6 = new TH1D(nameDenLong6.Data(),nameDenLong6.Data(),nbins,min,max);
   fDenLong7 = new TH1D(nameDenLong7.Data(),nameDenLong7.Data(),nbins,min,max);
   fDenLong8 = new TH1D(nameDenLong8.Data(),nameDenLong8.Data(),nbins,min,max);
   
   fRatLong1 = new TH1D(nameRatLong1.Data(),nameRatLong1.Data(),nbins,min,max);
   fRatLong2 = new TH1D(nameRatLong2.Data(),nameRatLong2.Data(),nbins,min,max);
   fRatLong3 = new TH1D(nameRatLong3.Data(),nameRatLong3.Data(),nbins,min,max);
   fRatLong4 = new TH1D(nameRatLong4.Data(),nameRatLong4.Data(),nbins,min,max);
   fRatLong5 = new TH1D(nameRatLong5.Data(),nameRatLong5.Data(),nbins,min,max);
   fRatLong6 = new TH1D(nameRatLong6.Data(),nameRatLong6.Data(),nbins,min,max);
   fRatLong7 = new TH1D(nameRatLong7.Data(),nameRatLong7.Data(),nbins,min,max);
   fRatLong8 = new TH1D(nameRatLong8.Data(),nameRatLong8.Data(),nbins,min,max);
   
   fNumOut1->Sumw2();
   fNumOut2->Sumw2();
     fNumOut3->Sumw2();
     fNumOut4->Sumw2();
     fNumOut5->Sumw2();
     fNumOut6->Sumw2();
     fNumOut7->Sumw2();
     fNumOut8->Sumw2();

     fDenOut1->Sumw2();
     fDenOut2->Sumw2();
     fDenOut3->Sumw2();
     fDenOut4->Sumw2();
     fDenOut5->Sumw2();
     fDenOut6->Sumw2();
     fDenOut7->Sumw2();
     fDenOut8->Sumw2();
   
     fRatOut1->Sumw2();
     fRatOut2->Sumw2();
     fRatOut3->Sumw2();
     fRatOut4->Sumw2();
     fRatOut5->Sumw2();
     fRatOut6->Sumw2();
     fRatOut7->Sumw2();
     fRatOut8->Sumw2();
   

     fNumSide1->Sumw2();
     fNumSide2->Sumw2();
     fNumSide3->Sumw2();
     fNumSide4->Sumw2();
     fNumSide5->Sumw2();
     fNumSide6->Sumw2();
     fNumSide7->Sumw2();
     fNumSide8->Sumw2();

     fDenSide1->Sumw2();
     fDenSide2->Sumw2();
     fDenSide3->Sumw2();
     fDenSide4->Sumw2();
     fDenSide5->Sumw2();
     fDenSide6->Sumw2();
     fDenSide7->Sumw2();
     fDenSide8->Sumw2();
   
     fRatSide1->Sumw2();
     fRatSide2->Sumw2();
     fRatSide3->Sumw2();
     fRatSide4->Sumw2();
     fRatSide5->Sumw2();
     fRatSide6->Sumw2();
     fRatSide7->Sumw2();
     fRatSide8->Sumw2();
   
     fNumLong1->Sumw2();
     fNumLong2->Sumw2();
     fNumLong3->Sumw2();
     fNumLong4->Sumw2();
     fNumLong5->Sumw2();
     fNumLong6->Sumw2();
     fNumLong7->Sumw2();
     fNumLong8->Sumw2();

     fDenLong1->Sumw2();
     fDenLong2->Sumw2();
     fDenLong3->Sumw2();
     fDenLong4->Sumw2();
     fDenLong5->Sumw2();
     fDenLong6->Sumw2();
     fDenLong7->Sumw2();
     fDenLong8->Sumw2();
   
     fRatLong1->Sumw2();
     fRatLong2->Sumw2();
     fRatLong3->Sumw2();
     fRatLong4->Sumw2();
     fRatLong5->Sumw2();
     fRatLong6->Sumw2();
     fRatLong7->Sumw2();
     fRatLong8->Sumw2();

 }

 TH1* AliHBTashbtCorrFctn::GetResult()
 {
     return fRatOut1;
 }
