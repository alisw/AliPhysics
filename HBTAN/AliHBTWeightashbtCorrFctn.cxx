#include "AliHBTWeightashbtCorrFctn.h"
#include <TH1.h>
#include <Riostream.h>

///////////////////////////////////////////////////////
//                                                   //
// AliHBTWeightashbtCorrFctn.h                             //
//                                                   //
// Class for calculating 3D ashbt correlation        //
// functions                                         //
//                                                   //
///////////////////////////////////////////////////////

ClassImp(AliHBTWeightashbtCorrFctn)

AliHBTWeightashbtCorrFctn::AliHBTWeightashbtCorrFctn(const char* name, const char* title):
 AliHBTOnePairFctn1D(name,title),

    fWeightNumOut1(0x0),
    fWeightNumOut2(0x0),
    fWeightNumOut3(0x0),
    fWeightNumOut4(0x0),
    fWeightNumOut5(0x0),
    fWeightNumOut6(0x0),
    fWeightNumOut7(0x0),
    fWeightNumOut8(0x0),

    fWeightDenOut1(0x0),
    fWeightDenOut2(0x0),
    fWeightDenOut3(0x0),
    fWeightDenOut4(0x0),
    fWeightDenOut5(0x0),
    fWeightDenOut6(0x0),
    fWeightDenOut7(0x0),
    fWeightDenOut8(0x0),
   
    fWeightRatOut1(0x0),
    fWeightRatOut2(0x0),
    fWeightRatOut3(0x0),
    fWeightRatOut4(0x0),
    fWeightRatOut5(0x0),
    fWeightRatOut6(0x0),
    fWeightRatOut7(0x0),
    fWeightRatOut8(0x0),
   

    fWeightNumSide1(0x0),
    fWeightNumSide2(0x0),
    fWeightNumSide3(0x0),
    fWeightNumSide4(0x0),
    fWeightNumSide5(0x0),
    fWeightNumSide6(0x0),
    fWeightNumSide7(0x0),
    fWeightNumSide8(0x0),

    fWeightDenSide1(0x0),
    fWeightDenSide2(0x0),
    fWeightDenSide3(0x0),
    fWeightDenSide4(0x0),
    fWeightDenSide5(0x0),
    fWeightDenSide6(0x0),
    fWeightDenSide7(0x0),
    fWeightDenSide8(0x0),
   
    fWeightRatSide1(0x0),
    fWeightRatSide2(0x0),
    fWeightRatSide3(0x0),
    fWeightRatSide4(0x0),
    fWeightRatSide5(0x0),
    fWeightRatSide6(0x0),
    fWeightRatSide7(0x0),
    fWeightRatSide8(0x0),
   
    fWeightNumLong1(0x0),
    fWeightNumLong2(0x0),
    fWeightNumLong3(0x0),
    fWeightNumLong4(0x0),
    fWeightNumLong5(0x0),
    fWeightNumLong6(0x0),
    fWeightNumLong7(0x0),
    fWeightNumLong8(0x0),

    fWeightDenLong1(0x0),
    fWeightDenLong2(0x0),
    fWeightDenLong3(0x0),
    fWeightDenLong4(0x0),
    fWeightDenLong5(0x0),
    fWeightDenLong6(0x0),
    fWeightDenLong7(0x0),
    fWeightDenLong8(0x0),
   
    fWeightRatLong1(0x0),
    fWeightRatLong2(0x0),
    fWeightRatLong3(0x0),
    fWeightRatLong4(0x0),
    fWeightRatLong5(0x0),
    fWeightRatLong6(0x0),
    fWeightRatLong7(0x0),
    fWeightRatLong8(0x0)


{
//ctor
}
/******************************************************************/
AliHBTWeightashbtCorrFctn::AliHBTWeightashbtCorrFctn(const char* name, const char* title, Int_t nbins, Float_t maxXval, Float_t minXval):
 AliHBTOnePairFctn1D(name,title,nbins,maxXval,minXval),

    fWeightNumOut1(0x0),
    fWeightNumOut2(0x0),
    fWeightNumOut3(0x0),
    fWeightNumOut4(0x0),
    fWeightNumOut5(0x0),
    fWeightNumOut6(0x0),
    fWeightNumOut7(0x0),
    fWeightNumOut8(0x0),

    fWeightDenOut1(0x0),
    fWeightDenOut2(0x0),
    fWeightDenOut3(0x0),
    fWeightDenOut4(0x0),
    fWeightDenOut5(0x0),
    fWeightDenOut6(0x0),
    fWeightDenOut7(0x0),
    fWeightDenOut8(0x0),
   
    fWeightRatOut1(0x0),
    fWeightRatOut2(0x0),
    fWeightRatOut3(0x0),
    fWeightRatOut4(0x0),
    fWeightRatOut5(0x0),
    fWeightRatOut6(0x0),
    fWeightRatOut7(0x0),
    fWeightRatOut8(0x0),
   

    fWeightNumSide1(0x0),
    fWeightNumSide2(0x0),
    fWeightNumSide3(0x0),
    fWeightNumSide4(0x0),
    fWeightNumSide5(0x0),
    fWeightNumSide6(0x0),
    fWeightNumSide7(0x0),
    fWeightNumSide8(0x0),

    fWeightDenSide1(0x0),
    fWeightDenSide2(0x0),
    fWeightDenSide3(0x0),
    fWeightDenSide4(0x0),
    fWeightDenSide5(0x0),
    fWeightDenSide6(0x0),
    fWeightDenSide7(0x0),
    fWeightDenSide8(0x0),
   
    fWeightRatSide1(0x0),
    fWeightRatSide2(0x0),
    fWeightRatSide3(0x0),
    fWeightRatSide4(0x0),
    fWeightRatSide5(0x0),
    fWeightRatSide6(0x0),
    fWeightRatSide7(0x0),
    fWeightRatSide8(0x0),
   
    fWeightNumLong1(0x0),
    fWeightNumLong2(0x0),
    fWeightNumLong3(0x0),
    fWeightNumLong4(0x0),
    fWeightNumLong5(0x0),
    fWeightNumLong6(0x0),
    fWeightNumLong7(0x0),
    fWeightNumLong8(0x0),

    fWeightDenLong1(0x0),
    fWeightDenLong2(0x0),
    fWeightDenLong3(0x0),
    fWeightDenLong4(0x0),
    fWeightDenLong5(0x0),
    fWeightDenLong6(0x0),
    fWeightDenLong7(0x0),
    fWeightDenLong8(0x0),
   
    fWeightRatLong1(0x0),
    fWeightRatLong2(0x0),
    fWeightRatLong3(0x0),
    fWeightRatLong4(0x0),
    fWeightRatLong5(0x0),
    fWeightRatLong6(0x0),
    fWeightRatLong7(0x0),
    fWeightRatLong8(0x0)


{
BuildHistos(nbins, maxXval, minXval);
//ctor
}

/******************************************************************/
AliHBTWeightashbtCorrFctn::AliHBTWeightashbtCorrFctn(const AliHBTWeightashbtCorrFctn& in):
 AliHBTOnePairFctn1D(in),
    fWeightNumOut1(0x0),
    fWeightNumOut2(0x0),
    fWeightNumOut3(0x0),
    fWeightNumOut4(0x0),
    fWeightNumOut5(0x0),
    fWeightNumOut6(0x0),
    fWeightNumOut7(0x0),
    fWeightNumOut8(0x0),

    fWeightDenOut1(0x0),
    fWeightDenOut2(0x0),
    fWeightDenOut3(0x0),
    fWeightDenOut4(0x0),
    fWeightDenOut5(0x0),
    fWeightDenOut6(0x0),
    fWeightDenOut7(0x0),
    fWeightDenOut8(0x0),
   
    fWeightRatOut1(0x0),
    fWeightRatOut2(0x0),
    fWeightRatOut3(0x0),
    fWeightRatOut4(0x0),
    fWeightRatOut5(0x0),
    fWeightRatOut6(0x0),
    fWeightRatOut7(0x0),
    fWeightRatOut8(0x0),
   

    fWeightNumSide1(0x0),
    fWeightNumSide2(0x0),
    fWeightNumSide3(0x0),
    fWeightNumSide4(0x0),
    fWeightNumSide5(0x0),
    fWeightNumSide6(0x0),
    fWeightNumSide7(0x0),
    fWeightNumSide8(0x0),

    fWeightDenSide1(0x0),
    fWeightDenSide2(0x0),
    fWeightDenSide3(0x0),
    fWeightDenSide4(0x0),
    fWeightDenSide5(0x0),
    fWeightDenSide6(0x0),
    fWeightDenSide7(0x0),
    fWeightDenSide8(0x0),
   
    fWeightRatSide1(0x0),
    fWeightRatSide2(0x0),
    fWeightRatSide3(0x0),
    fWeightRatSide4(0x0),
    fWeightRatSide5(0x0),
    fWeightRatSide6(0x0),
    fWeightRatSide7(0x0),
    fWeightRatSide8(0x0),
   
    fWeightNumLong1(0x0),
    fWeightNumLong2(0x0),
    fWeightNumLong3(0x0),
    fWeightNumLong4(0x0),
    fWeightNumLong5(0x0),
    fWeightNumLong6(0x0),
    fWeightNumLong7(0x0),
    fWeightNumLong8(0x0),

    fWeightDenLong1(0x0),
    fWeightDenLong2(0x0),
    fWeightDenLong3(0x0),
    fWeightDenLong4(0x0),
    fWeightDenLong5(0x0),
    fWeightDenLong6(0x0),
    fWeightDenLong7(0x0),
    fWeightDenLong8(0x0),
   
    fWeightRatLong1(0x0),
    fWeightRatLong2(0x0),
    fWeightRatLong3(0x0),
    fWeightRatLong4(0x0),
    fWeightRatLong5(0x0),
    fWeightRatLong6(0x0),
    fWeightRatLong7(0x0),
    fWeightRatLong8(0x0)
  

 {
//ctor
}

/******************************************************************/

AliHBTWeightashbtCorrFctn::~AliHBTWeightashbtCorrFctn()
{
 //dtor
    delete fWeightNumOut1;
    delete fWeightNumOut2;
    delete fWeightNumOut3;
    delete fWeightNumOut4;
    delete fWeightNumOut5;
    delete fWeightNumOut6;
    delete fWeightNumOut7;
    delete fWeightNumOut8;

    delete fWeightDenOut1;
    delete fWeightDenOut2;
    delete fWeightDenOut3;
    delete fWeightDenOut4;
    delete fWeightDenOut5;
    delete fWeightDenOut6;
    delete fWeightDenOut7;
    delete fWeightDenOut8;
   
    delete fWeightRatOut1;
    delete fWeightRatOut2;
    delete fWeightRatOut3;
    delete fWeightRatOut4;
    delete fWeightRatOut5;
    delete fWeightRatOut6;
    delete fWeightRatOut7;
    delete fWeightRatOut8;
   

    delete fWeightNumSide1;
    delete fWeightNumSide2;
    delete fWeightNumSide3;
    delete fWeightNumSide4;
    delete fWeightNumSide5;
    delete fWeightNumSide6;
    delete fWeightNumSide7;
    delete fWeightNumSide8;

    delete fWeightDenSide1;
    delete fWeightDenSide2;
    delete fWeightDenSide3;
    delete fWeightDenSide4;
    delete fWeightDenSide5;
    delete fWeightDenSide6;
    delete fWeightDenSide7;
    delete fWeightDenSide8;
   
    delete fWeightRatSide1;
    delete fWeightRatSide2;
    delete fWeightRatSide3;
    delete fWeightRatSide4;
    delete fWeightRatSide5;
    delete fWeightRatSide6;
    delete fWeightRatSide7;
    delete fWeightRatSide8;
   
    delete fWeightNumLong1;
    delete fWeightNumLong2;
    delete fWeightNumLong3;
    delete fWeightNumLong4;
    delete fWeightNumLong5;
    delete fWeightNumLong6;
    delete fWeightNumLong7;
    delete fWeightNumLong8;

    delete fWeightDenLong1;
    delete fWeightDenLong2;
    delete fWeightDenLong3;
    delete fWeightDenLong4;
    delete fWeightDenLong5;
    delete fWeightDenLong6;
    delete fWeightDenLong7;
    delete fWeightDenLong8;
   
    delete fWeightRatLong1;
    delete fWeightRatLong2;
    delete fWeightRatLong3;
    delete fWeightRatLong4;
    delete fWeightRatLong5;
    delete fWeightRatLong6;
    delete fWeightRatLong7;
    delete fWeightRatLong8; 

}

/******************************************************************/
void AliHBTWeightashbtCorrFctn::WriteFunction()
{
//out    
   Double_t out1scale = Scale(fWeightNumOut1,fWeightDenOut1);
   cout <<"out1scale = "<<out1scale<<endl;
   fWeightRatOut1->Divide(fWeightNumOut1,fWeightDenOut1,out1scale);
   
   Double_t out2scale = Scale(fWeightNumOut2,fWeightDenOut2);
   cout <<"out2scale = "<<out2scale<<endl;
   fWeightRatOut2->Divide(fWeightNumOut2,fWeightDenOut2,out2scale);
   
   Double_t out3scale = Scale(fWeightNumOut3,fWeightDenOut3);
   cout <<"out3scale = "<<out3scale<<endl;
   fWeightRatOut3->Divide(fWeightNumOut3,fWeightDenOut3,out3scale);
   
   Double_t out4scale = Scale(fWeightNumOut4,fWeightDenOut4);
   cout <<"out4scale = "<<out4scale<<endl;
   fWeightRatOut4->Divide(fWeightNumOut4,fWeightDenOut4,out4scale);
   
   Double_t out5scale = Scale(fWeightNumOut5,fWeightDenOut5);
   cout <<"out5scale = "<<out5scale<<endl;
   fWeightRatOut5->Divide(fWeightNumOut5,fWeightDenOut5,out5scale);
   
   Double_t out6scale = Scale(fWeightNumOut6,fWeightDenOut6);
   cout <<"out6scale = "<<out6scale<<endl;
   fWeightRatOut6->Divide(fWeightNumOut6,fWeightDenOut6,out6scale);
   
   Double_t out7scale = Scale(fWeightNumOut7,fWeightDenOut7);
   cout <<"out7scale = "<<out7scale<<endl;
   fWeightRatOut7->Divide(fWeightNumOut7,fWeightDenOut7,out7scale);
   
   Double_t out8scale = Scale(fWeightNumOut8,fWeightDenOut8);
   cout <<"out8scale = "<<out8scale<<endl;
   fWeightRatOut8->Divide(fWeightNumOut8,fWeightDenOut8,out8scale);

   //side
   Double_t side1scale = Scale(fWeightNumSide1,fWeightDenSide1);
   cout <<"side1scale = "<<side1scale<<endl;
   fWeightRatSide1->Divide(fWeightNumSide1,fWeightDenSide1,side1scale);
   
   Double_t side2scale = Scale(fWeightNumSide2,fWeightDenSide2);
   cout <<"side2scale = "<<side2scale<<endl;
   fWeightRatSide2->Divide(fWeightNumSide2,fWeightDenSide2,side2scale);
   
   Double_t side3scale = Scale(fWeightNumSide3,fWeightDenSide3);
   cout <<"side3scale = "<<side3scale<<endl;
   fWeightRatSide3->Divide(fWeightNumSide3,fWeightDenSide3,side3scale);
   
   Double_t side4scale = Scale(fWeightNumSide4,fWeightDenSide4);
   cout <<"side4scale = "<<side4scale<<endl;
   fWeightRatSide4->Divide(fWeightNumSide4,fWeightDenSide4,side4scale);
   
   Double_t side5scale = Scale(fWeightNumSide5,fWeightDenSide5);
   cout <<"side5scale = "<<side5scale<<endl;
   fWeightRatSide5->Divide(fWeightNumSide5,fWeightDenSide5,side5scale);
   
   Double_t side6scale = Scale(fWeightNumSide6,fWeightDenSide6);
   cout <<"side6scale = "<<side6scale<<endl;
   fWeightRatSide6->Divide(fWeightNumSide6,fWeightDenSide6,side6scale);
   
   Double_t side7scale = Scale(fWeightNumSide7,fWeightDenSide7);
   cout <<"side7scale = "<<side7scale<<endl;
   fWeightRatSide7->Divide(fWeightNumSide7,fWeightDenSide7,side7scale);
   
   Double_t side8scale = Scale(fWeightNumSide8,fWeightDenSide8);
   cout <<"side8scale = "<<side8scale<<endl;
   fWeightRatSide8->Divide(fWeightNumSide8,fWeightDenSide8,side8scale);
   
//long
 
   Double_t long1scale = Scale(fWeightNumLong1,fWeightDenLong1);
   cout <<"long1scale = "<<long1scale<<endl;
   fWeightRatLong1->Divide(fWeightNumLong1,fWeightDenLong1,long1scale);
   
   Double_t long2scale = Scale(fWeightNumLong2,fWeightDenLong2);
   cout <<"long2scale = "<<long2scale<<endl;
   fWeightRatLong2->Divide(fWeightNumLong2,fWeightDenLong2,long2scale);
   
   Double_t long3scale = Scale(fWeightNumLong3,fWeightDenLong3);
   cout <<"long3scale = "<<long3scale<<endl;
   fWeightRatLong3->Divide(fWeightNumLong3,fWeightDenLong3,long3scale);
   
   Double_t long4scale = Scale(fWeightNumLong4,fWeightDenLong4);
   cout <<"long4scale = "<<long4scale<<endl;
   fWeightRatLong4->Divide(fWeightNumLong4,fWeightDenLong4,long4scale);
   
   Double_t long5scale = Scale(fWeightNumLong5,fWeightDenLong5);
   cout <<"long5scale = "<<long5scale<<endl;
   fWeightRatLong5->Divide(fWeightNumLong5,fWeightDenLong5,long5scale);
   
   Double_t long6scale = Scale(fWeightNumLong6,fWeightDenLong6);
   cout <<"long6scale = "<<long6scale<<endl;
   fWeightRatLong6->Divide(fWeightNumLong6,fWeightDenLong6,long6scale);
   
   Double_t long7scale = Scale(fWeightNumLong7,fWeightDenLong7);
   cout <<"long7scale = "<<long7scale<<endl;
   fWeightRatLong7->Divide(fWeightNumLong7,fWeightDenLong7,long7scale);
   
   Double_t long8scale = Scale(fWeightNumLong8,fWeightDenLong8);
   cout <<"long8scale = "<<long8scale<<endl;
   fWeightRatLong8->Divide(fWeightNumLong8,fWeightDenLong8,long8scale);

        fWeightNumOut1->Write();
     fWeightNumOut2->Write();
     fWeightNumOut3->Write();
     fWeightNumOut4->Write();
     fWeightNumOut5->Write();
     fWeightNumOut6->Write();
     fWeightNumOut7->Write();
     fWeightNumOut8->Write();

     fWeightDenOut1->Write();
     fWeightDenOut2->Write();
     fWeightDenOut3->Write();
     fWeightDenOut4->Write();
     fWeightDenOut5->Write();
     fWeightDenOut6->Write();
     fWeightDenOut7->Write();
     fWeightDenOut8->Write();
   
     fWeightRatOut1->Write();
     fWeightRatOut2->Write();
     fWeightRatOut3->Write();
     fWeightRatOut4->Write();
     fWeightRatOut5->Write();
     fWeightRatOut6->Write();
     fWeightRatOut7->Write();
     fWeightRatOut8->Write();
   

     fWeightNumSide1->Write();
     fWeightNumSide2->Write();
     fWeightNumSide3->Write();
     fWeightNumSide4->Write();
     fWeightNumSide5->Write();
     fWeightNumSide6->Write();
     fWeightNumSide7->Write();
     fWeightNumSide8->Write();

     fWeightDenSide1->Write();
     fWeightDenSide2->Write();
     fWeightDenSide3->Write();
     fWeightDenSide4->Write();
     fWeightDenSide5->Write();
     fWeightDenSide6->Write();
     fWeightDenSide7->Write();
     fWeightDenSide8->Write();
   
     fWeightRatSide1->Write();
     fWeightRatSide2->Write();
     fWeightRatSide3->Write();
     fWeightRatSide4->Write();
     fWeightRatSide5->Write();
     fWeightRatSide6->Write();
     fWeightRatSide7->Write();
     fWeightRatSide8->Write();
   
     fWeightNumLong1->Write();
     fWeightNumLong2->Write();
     fWeightNumLong3->Write();
     fWeightNumLong4->Write();
     fWeightNumLong5->Write();
     fWeightNumLong6->Write();
     fWeightNumLong7->Write();
     fWeightNumLong8->Write();

     fWeightDenLong1->Write();
     fWeightDenLong2->Write();
     fWeightDenLong3->Write();
     fWeightDenLong4->Write();
     fWeightDenLong5->Write();
     fWeightDenLong6->Write();
     fWeightDenLong7->Write();
     fWeightDenLong8->Write();
   
     fWeightRatLong1->Write();
     fWeightRatLong2->Write();
     fWeightRatLong3->Write();
     fWeightRatLong4->Write();
     fWeightRatLong5->Write();
     fWeightRatLong6->Write();
     fWeightRatLong7->Write();
     fWeightRatLong8->Write(); 

}

//-------------------------------------
void AliHBTWeightashbtCorrFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
    //Fills the numerator using pair from the same event
   
   trackpair = CheckPair(trackpair);
   if(partpair == 0x0) return;
   if(trackpair == 0x0) return;
 
   Double_t weight = 1.0;
   
   Double_t rplane=0.;   //reaction plane angle - 2 B determined
   Double_t phi=(trackpair->Particle1()->Phi()+trackpair->Particle2()->Phi())/2.-rplane; //deltaphi bo nie mam nic innego pod reka
   phi=phi*360/(2*TMath::Pi());
   Double_t qout=trackpair->GetQOutLCMS();
   Double_t qside=trackpair->GetQSideLCMS();
   Double_t qlong=trackpair->GetQLongLCMS();
   
   if ( trackpair && partpair)
     {
	if ( ( trackpair->Particle1()->GetPdgCode() == partpair->Particle1()->GetPdgCode()) &&
	     ( trackpair->Particle2()->GetPdgCode() == partpair->Particle2()->GetPdgCode())    )
	  {
	     weight=partpair->GetWeight();
	  }
	
	
	if((phi>=0.) && (phi <45.))
	  {
	     fWeightNumOut1->Fill(qout,weight);
	     fWeightNumSide1->Fill(qside,weight);
	     fWeightNumLong1->Fill(qlong,weight);
	  }
	else if((phi>=45.) && (phi <90.))
	  {
	     fWeightNumOut2->Fill(qout,weight);
	     fWeightNumSide2->Fill(qside,weight);
	     fWeightNumLong2->Fill(qlong,weight);
	  }
	else if((phi>=90.) && (phi <135.))
	  {
	     fWeightNumOut3->Fill(qout,weight);
	     fWeightNumSide3->Fill(qside,weight);
	     fWeightNumLong3->Fill(qlong,weight);
	  }
	else if((phi>=135.) && (phi <180.))
	  {
	     fWeightNumOut4->Fill(qout,weight);
	     fWeightNumSide4->Fill(qside,weight);
	     fWeightNumLong4->Fill(qlong,weight);
	  }
	else if((phi>=180.) && (phi <225.))
	  {
	     fWeightNumOut5->Fill(qout,weight);
	     fWeightNumSide5->Fill(qside,weight);
	     fWeightNumLong5->Fill(qlong,weight);
	  }
	else if((phi>=225.) && (phi <270.))
	  {
	     fWeightNumOut6->Fill(qout,weight);
	     fWeightNumSide6->Fill(qside,weight);
	     fWeightNumLong6->Fill(qlong,weight);
	  }
	else if((phi>=270.) && (phi <315.))
	  {
	     fWeightNumOut7->Fill(qout,weight);
	     fWeightNumSide7->Fill(qside,weight);
	     fWeightNumLong7->Fill(qlong,weight);
	  }
	else if((phi>=315.) && (phi <360.))
	  {
	     fWeightNumOut8->Fill(qout,weight);
	     fWeightNumSide8->Fill(qside,weight);
	     fWeightNumLong8->Fill(qlong,weight);
	  }
     }
   
   
}

/****************************************************************/
void AliHBTWeightashbtCorrFctn::Init()
{
     fWeightNumOut1->Reset();
     fWeightNumOut2->Reset();
     fWeightNumOut3->Reset();
     fWeightNumOut4->Reset();
     fWeightNumOut5->Reset();
     fWeightNumOut6->Reset();
     fWeightNumOut7->Reset();
     fWeightNumOut8->Reset();

     fWeightDenOut1->Reset();
     fWeightDenOut2->Reset();
     fWeightDenOut3->Reset();
     fWeightDenOut4->Reset();
     fWeightDenOut5->Reset();
     fWeightDenOut6->Reset();
     fWeightDenOut7->Reset();
     fWeightDenOut8->Reset();
   
     fWeightRatOut1->Reset();
     fWeightRatOut2->Reset();
     fWeightRatOut3->Reset();
     fWeightRatOut4->Reset();
     fWeightRatOut5->Reset();
     fWeightRatOut6->Reset();
     fWeightRatOut7->Reset();
     fWeightRatOut8->Reset();
   

     fWeightNumSide1->Reset();
     fWeightNumSide2->Reset();
     fWeightNumSide3->Reset();
     fWeightNumSide4->Reset();
     fWeightNumSide5->Reset();
     fWeightNumSide6->Reset();
     fWeightNumSide7->Reset();
     fWeightNumSide8->Reset();

     fWeightDenSide1->Reset();
     fWeightDenSide2->Reset();
     fWeightDenSide3->Reset();
     fWeightDenSide4->Reset();
     fWeightDenSide5->Reset();
     fWeightDenSide6->Reset();
     fWeightDenSide7->Reset();
     fWeightDenSide8->Reset();
   
     fWeightRatSide1->Reset();
     fWeightRatSide2->Reset();
     fWeightRatSide3->Reset();
     fWeightRatSide4->Reset();
     fWeightRatSide5->Reset();
     fWeightRatSide6->Reset();
     fWeightRatSide7->Reset();
     fWeightRatSide8->Reset();
   
     fWeightNumLong1->Reset();
     fWeightNumLong2->Reset();
     fWeightNumLong3->Reset();
     fWeightNumLong4->Reset();
     fWeightNumLong5->Reset();
     fWeightNumLong6->Reset();
     fWeightNumLong7->Reset();
     fWeightNumLong8->Reset();

     fWeightDenLong1->Reset();
     fWeightDenLong2->Reset();
     fWeightDenLong3->Reset();
     fWeightDenLong4->Reset();
     fWeightDenLong5->Reset();
     fWeightDenLong6->Reset();
     fWeightDenLong7->Reset();
     fWeightDenLong8->Reset();
   
     fWeightRatLong1->Reset();
     fWeightRatLong2->Reset();
     fWeightRatLong3->Reset();
     fWeightRatLong4->Reset();
     fWeightRatLong5->Reset();
     fWeightRatLong6->Reset();
     fWeightRatLong7->Reset();
     fWeightRatLong8->Reset(); 
   

 }
/****************************************************************/


void AliHBTWeightashbtCorrFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{

   trackpair = CheckPair(trackpair);
   Double_t rplane=0.;   //reaction plane angle - 2 B determined
   Double_t phi=(trackpair->Particle1()->Phi()+trackpair->Particle2()->Phi())/2.-rplane; //deltaphi bo nie mam nic innego pod reka
   phi=phi*360/(2*TMath::Pi());
   Double_t qout=trackpair->GetQOutLCMS();
   Double_t qside=trackpair->GetQSideLCMS();
   Double_t qlong=trackpair->GetQLongLCMS();
   

if ( trackpair && partpair)
     {   
	
	if((phi>=0.) && (phi <45.))
	  {
	     fWeightDenOut1->Fill(qout);
	     fWeightDenSide1->Fill(qside);
	     fWeightDenLong1->Fill(qlong);
	  }
	else if((phi>=45.) && (phi <90.))
	  {
	     fWeightDenOut2->Fill(qout);
	     fWeightDenSide2->Fill(qside);
	     fWeightDenLong2->Fill(qlong);
	  }
	else if((phi>=90.) && (phi <135.))
	  {
	     fWeightDenOut3->Fill(qout);
	     fWeightDenSide3->Fill(qside);
	     fWeightDenLong3->Fill(qlong);
	  }
	else if((phi>=135.) && (phi <180.))
	  {
	     fWeightDenOut4->Fill(qout);
	     fWeightDenSide4->Fill(qside);
	     fWeightDenLong4->Fill(qlong);
	  }
	else if((phi>=180.) && (phi <225.))
	  {
	     fWeightDenOut5->Fill(qout);
	     fWeightDenSide5->Fill(qside);
	     fWeightDenLong5->Fill(qlong);
	  }
	else if((phi>=225.) && (phi <270.))
	  {
	     fWeightDenOut6->Fill(qout);
	     fWeightDenSide6->Fill(qside);
	     fWeightDenLong6->Fill(qlong);
	  }
	else if((phi>=270.) && (phi <315.))
	  {
	     fWeightDenOut7->Fill(qout);
	     fWeightDenSide7->Fill(qside);
	     fWeightDenLong7->Fill(qlong);
	  }
	else if((phi>=315.) && (phi <360.))
	  {
	     fWeightDenOut8->Fill(qout);
	     fWeightDenSide8->Fill(qside);
	     fWeightDenLong8->Fill(qlong);
	  }
     }
   
}


/******************************************************************/

void AliHBTWeightashbtCorrFctn::BuildHistos(Int_t nbins, Float_t max, Float_t min)
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

   fWeightNumOut1 = new TH1D(nameNumOut1.Data(),nameNumOut1.Data(),nbins,min,max);
   fWeightNumOut2 = new TH1D(nameNumOut2.Data(),nameNumOut2.Data(),nbins,min,max);
   fWeightNumOut3 = new TH1D(nameNumOut3.Data(),nameNumOut3.Data(),nbins,min,max);
   fWeightNumOut4 = new TH1D(nameNumOut4.Data(),nameNumOut4.Data(),nbins,min,max);
   fWeightNumOut5 = new TH1D(nameNumOut5.Data(),nameNumOut5.Data(),nbins,min,max);
   fWeightNumOut6 = new TH1D(nameNumOut6.Data(),nameNumOut6.Data(),nbins,min,max);
   fWeightNumOut7 = new TH1D(nameNumOut7.Data(),nameNumOut7.Data(),nbins,min,max);
   fWeightNumOut8 = new TH1D(nameNumOut8.Data(),nameNumOut8.Data(),nbins,min,max);
   
   fWeightDenOut1 = new TH1D(nameDenOut1.Data(),nameDenOut1.Data(),nbins,min,max);
   fWeightDenOut2 = new TH1D(nameDenOut2.Data(),nameDenOut2.Data(),nbins,min,max);
   fWeightDenOut3 = new TH1D(nameDenOut3.Data(),nameDenOut3.Data(),nbins,min,max);
   fWeightDenOut4 = new TH1D(nameDenOut4.Data(),nameDenOut4.Data(),nbins,min,max);
   fWeightDenOut5 = new TH1D(nameDenOut5.Data(),nameDenOut5.Data(),nbins,min,max);
   fWeightDenOut6 = new TH1D(nameDenOut6.Data(),nameDenOut6.Data(),nbins,min,max);
   fWeightDenOut7 = new TH1D(nameDenOut7.Data(),nameDenOut7.Data(),nbins,min,max);
   fWeightDenOut8 = new TH1D(nameDenOut8.Data(),nameDenOut8.Data(),nbins,min,max);
   
   fWeightRatOut1 = new TH1D(nameRatOut1.Data(),nameRatOut1.Data(),nbins,min,max);
   fWeightRatOut2 = new TH1D(nameRatOut2.Data(),nameRatOut2.Data(),nbins,min,max);
   fWeightRatOut3 = new TH1D(nameRatOut3.Data(),nameRatOut3.Data(),nbins,min,max);
   fWeightRatOut4 = new TH1D(nameRatOut4.Data(),nameRatOut4.Data(),nbins,min,max);
   fWeightRatOut5 = new TH1D(nameRatOut5.Data(),nameRatOut5.Data(),nbins,min,max);
   fWeightRatOut6 = new TH1D(nameRatOut6.Data(),nameRatOut6.Data(),nbins,min,max);
   fWeightRatOut7 = new TH1D(nameRatOut7.Data(),nameRatOut7.Data(),nbins,min,max);
   fWeightRatOut8 = new TH1D(nameRatOut8.Data(),nameRatOut8.Data(),nbins,min,max);
   
   fWeightNumSide1 = new TH1D(nameNumSide1.Data(),nameNumSide1.Data(),nbins,min,max);
   fWeightNumSide2 = new TH1D(nameNumSide2.Data(),nameNumSide2.Data(),nbins,min,max);
   fWeightNumSide3 = new TH1D(nameNumSide3.Data(),nameNumSide3.Data(),nbins,min,max);
   fWeightNumSide4 = new TH1D(nameNumSide4.Data(),nameNumSide4.Data(),nbins,min,max);
   fWeightNumSide5 = new TH1D(nameNumSide5.Data(),nameNumSide5.Data(),nbins,min,max);
   fWeightNumSide6 = new TH1D(nameNumSide6.Data(),nameNumSide6.Data(),nbins,min,max);
   fWeightNumSide7 = new TH1D(nameNumSide7.Data(),nameNumSide7.Data(),nbins,min,max);
   fWeightNumSide8 = new TH1D(nameNumSide8.Data(),nameNumSide8.Data(),nbins,min,max);
   
   fWeightDenSide1 = new TH1D(nameDenSide1.Data(),nameDenSide1.Data(),nbins,min,max);
   fWeightDenSide2 = new TH1D(nameDenSide2.Data(),nameDenSide2.Data(),nbins,min,max);
   fWeightDenSide3 = new TH1D(nameDenSide3.Data(),nameDenSide3.Data(),nbins,min,max);
   fWeightDenSide4 = new TH1D(nameDenSide4.Data(),nameDenSide4.Data(),nbins,min,max);
   fWeightDenSide5 = new TH1D(nameDenSide5.Data(),nameDenSide5.Data(),nbins,min,max);
   fWeightDenSide6 = new TH1D(nameDenSide6.Data(),nameDenSide6.Data(),nbins,min,max);
   fWeightDenSide7 = new TH1D(nameDenSide7.Data(),nameDenSide7.Data(),nbins,min,max);
   fWeightDenSide8 = new TH1D(nameDenSide8.Data(),nameDenSide8.Data(),nbins,min,max);
   
   fWeightRatSide1 = new TH1D(nameRatSide1.Data(),nameRatSide1.Data(),nbins,min,max);
   fWeightRatSide2 = new TH1D(nameRatSide2.Data(),nameRatSide2.Data(),nbins,min,max);
   fWeightRatSide3 = new TH1D(nameRatSide3.Data(),nameRatSide3.Data(),nbins,min,max);
   fWeightRatSide4 = new TH1D(nameRatSide4.Data(),nameRatSide4.Data(),nbins,min,max);
   fWeightRatSide5 = new TH1D(nameRatSide5.Data(),nameRatSide5.Data(),nbins,min,max);
   fWeightRatSide6 = new TH1D(nameRatSide6.Data(),nameRatSide6.Data(),nbins,min,max);
   fWeightRatSide7 = new TH1D(nameRatSide7.Data(),nameRatSide7.Data(),nbins,min,max);
   fWeightRatSide8 = new TH1D(nameRatSide8.Data(),nameRatSide8.Data(),nbins,min,max);
   
   fWeightNumLong1 = new TH1D(nameNumLong1.Data(),nameNumLong1.Data(),nbins,min,max);
   fWeightNumLong2 = new TH1D(nameNumLong2.Data(),nameNumLong2.Data(),nbins,min,max);
   fWeightNumLong3 = new TH1D(nameNumLong3.Data(),nameNumLong3.Data(),nbins,min,max);
   fWeightNumLong4 = new TH1D(nameNumLong4.Data(),nameNumLong4.Data(),nbins,min,max);
   fWeightNumLong5 = new TH1D(nameNumLong5.Data(),nameNumLong5.Data(),nbins,min,max);
   fWeightNumLong6 = new TH1D(nameNumLong6.Data(),nameNumLong6.Data(),nbins,min,max);
   fWeightNumLong7 = new TH1D(nameNumLong7.Data(),nameNumLong7.Data(),nbins,min,max);
   fWeightNumLong8 = new TH1D(nameNumLong8.Data(),nameNumLong8.Data(),nbins,min,max);
   
   fWeightDenLong1 = new TH1D(nameDenLong1.Data(),nameDenLong1.Data(),nbins,min,max);
   fWeightDenLong2 = new TH1D(nameDenLong2.Data(),nameDenLong2.Data(),nbins,min,max);
   fWeightDenLong3 = new TH1D(nameDenLong3.Data(),nameDenLong3.Data(),nbins,min,max);
   fWeightDenLong4 = new TH1D(nameDenLong4.Data(),nameDenLong4.Data(),nbins,min,max);
   fWeightDenLong5 = new TH1D(nameDenLong5.Data(),nameDenLong5.Data(),nbins,min,max);
   fWeightDenLong6 = new TH1D(nameDenLong6.Data(),nameDenLong6.Data(),nbins,min,max);
   fWeightDenLong7 = new TH1D(nameDenLong7.Data(),nameDenLong7.Data(),nbins,min,max);
   fWeightDenLong8 = new TH1D(nameDenLong8.Data(),nameDenLong8.Data(),nbins,min,max);
   
   fWeightRatLong1 = new TH1D(nameRatLong1.Data(),nameRatLong1.Data(),nbins,min,max);
   fWeightRatLong2 = new TH1D(nameRatLong2.Data(),nameRatLong2.Data(),nbins,min,max);
   fWeightRatLong3 = new TH1D(nameRatLong3.Data(),nameRatLong3.Data(),nbins,min,max);
   fWeightRatLong4 = new TH1D(nameRatLong4.Data(),nameRatLong4.Data(),nbins,min,max);
   fWeightRatLong5 = new TH1D(nameRatLong5.Data(),nameRatLong5.Data(),nbins,min,max);
   fWeightRatLong6 = new TH1D(nameRatLong6.Data(),nameRatLong6.Data(),nbins,min,max);
   fWeightRatLong7 = new TH1D(nameRatLong7.Data(),nameRatLong7.Data(),nbins,min,max);
   fWeightRatLong8 = new TH1D(nameRatLong8.Data(),nameRatLong8.Data(),nbins,min,max);
   
   fWeightNumOut1->Sumw2();
   fWeightNumOut2->Sumw2();
     fWeightNumOut3->Sumw2();
     fWeightNumOut4->Sumw2();
     fWeightNumOut5->Sumw2();
     fWeightNumOut6->Sumw2();
     fWeightNumOut7->Sumw2();
     fWeightNumOut8->Sumw2();

     fWeightDenOut1->Sumw2();
     fWeightDenOut2->Sumw2();
     fWeightDenOut3->Sumw2();
     fWeightDenOut4->Sumw2();
     fWeightDenOut5->Sumw2();
     fWeightDenOut6->Sumw2();
     fWeightDenOut7->Sumw2();
     fWeightDenOut8->Sumw2();
   
     fWeightRatOut1->Sumw2();
     fWeightRatOut2->Sumw2();
     fWeightRatOut3->Sumw2();
     fWeightRatOut4->Sumw2();
     fWeightRatOut5->Sumw2();
     fWeightRatOut6->Sumw2();
     fWeightRatOut7->Sumw2();
     fWeightRatOut8->Sumw2();
   

     fWeightNumSide1->Sumw2();
     fWeightNumSide2->Sumw2();
     fWeightNumSide3->Sumw2();
     fWeightNumSide4->Sumw2();
     fWeightNumSide5->Sumw2();
     fWeightNumSide6->Sumw2();
     fWeightNumSide7->Sumw2();
     fWeightNumSide8->Sumw2();

     fWeightDenSide1->Sumw2();
     fWeightDenSide2->Sumw2();
     fWeightDenSide3->Sumw2();
     fWeightDenSide4->Sumw2();
     fWeightDenSide5->Sumw2();
     fWeightDenSide6->Sumw2();
     fWeightDenSide7->Sumw2();
     fWeightDenSide8->Sumw2();
   
     fWeightRatSide1->Sumw2();
     fWeightRatSide2->Sumw2();
     fWeightRatSide3->Sumw2();
     fWeightRatSide4->Sumw2();
     fWeightRatSide5->Sumw2();
     fWeightRatSide6->Sumw2();
     fWeightRatSide7->Sumw2();
     fWeightRatSide8->Sumw2();
   
     fWeightNumLong1->Sumw2();
     fWeightNumLong2->Sumw2();
     fWeightNumLong3->Sumw2();
     fWeightNumLong4->Sumw2();
     fWeightNumLong5->Sumw2();
     fWeightNumLong6->Sumw2();
     fWeightNumLong7->Sumw2();
     fWeightNumLong8->Sumw2();

     fWeightDenLong1->Sumw2();
     fWeightDenLong2->Sumw2();
     fWeightDenLong3->Sumw2();
     fWeightDenLong4->Sumw2();
     fWeightDenLong5->Sumw2();
     fWeightDenLong6->Sumw2();
     fWeightDenLong7->Sumw2();
     fWeightDenLong8->Sumw2();
   
     fWeightRatLong1->Sumw2();
     fWeightRatLong2->Sumw2();
     fWeightRatLong3->Sumw2();
     fWeightRatLong4->Sumw2();
     fWeightRatLong5->Sumw2();
     fWeightRatLong6->Sumw2();
     fWeightRatLong7->Sumw2();
     fWeightRatLong8->Sumw2();

 }

 TH1* AliHBTWeightashbtCorrFctn::GetResult()
 {
     return fWeightRatOut1;
 }
