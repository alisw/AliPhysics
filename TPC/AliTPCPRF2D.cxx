/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
Revision 1.2  1999/09/29 09:24:34  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//  AliTPCPRF2D -                                                                         //
//  Pad response function object in two dimesions                            //
//  This class contains the basic functions for the                          //
//  calculation of PRF according generic charge distribution                 //
//  In Update function object calculate table of response function           //
//  in discrete x and y position                                             //
// This table is used for interpolation od response function in any position //
// (function GetPRF)                                                          //
//                                                                           // 
//  Origin: Marian Ivanov, Uni. of Bratislava, ivanov@fmph.uniba.sk          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "TMath.h"
#include "AliTPCPRF2D.h"
#include "TF2.h"
#include <iostream.h>
#include <string.h>
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TPaveText.h"
#include "TText.h"

extern TStyle * gStyle;

static const Float_t sqrt12=3.46;
static const Int_t   NPRF = 100;


static Double_t funGauss2D(Double_t *x, Double_t * par)
{ 
  return ( TMath::Exp(-(x[0]*x[0])/(2*par[0]*par[0]))*
	   TMath::Exp(-(x[1]*x[1])/(2*par[1]*par[1])));

}

static Double_t funCosh2D(Double_t *x, Double_t * par)
{
  return ( 1/(TMath::CosH(3.14159*x[0]/(2*par[0]))*
	   TMath::CosH(3.14159*x[1]/(2*par[1]))));
}    

static Double_t funGati2D(Double_t *x, Double_t * par)
{
  //par[1] = is equal to k3X
  //par[0] is equal to pad wire distance
  Float_t K3=par[1];
  Float_t K3R=TMath::Sqrt(K3);
  Float_t K2=(TMath::Pi()/2)*(1-K3R/2.);
  Float_t K1=K2*K3R/(4*TMath::ATan(K3R));
  Float_t l=x[0]/par[0];
  Float_t tan2=TMath::TanH(K2*l);
  tan2*=tan2;
  Float_t res = K1*(1-tan2)/(1+K3*tan2);
 //par[4] = is equal to k3Y
  K3=par[4];
  K3R=TMath::Sqrt(K3);
  K2=(TMath::Pi()/2)*(1-K3R/2.);
  K1=K2*K3R/(4*TMath::ATan(K3R));
  l=x[1]/par[0];
  tan2=TMath::TanH(K2*l);
  tan2*=tan2;
  res = res*K1*(1-tan2)/(1+K3*tan2);  
  return res;  
}   


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

ClassImp(AliTPCPRF2D)

AliTPCPRF2D::AliTPCPRF2D()
{
  ffcharge = 0;
  fNPRF =NPRF ;
  fSigmaX = 0;

  fGRF = 0;
  fkNorm = 1;
  forigsigmaY=0;
  forigsigmaX=0;
  fNdiv = 5;
  //chewron default values   
  SetPad(0.8,0.8);
  SetChevron(0.2,0.0,1.0);
  SetY(-0.2,0.2,2);
  // SetGauss(0.22,0.22,1);  
}

AliTPCPRF2D::~AliTPCPRF2D()
{
  if (ffcharge!=0) delete [] ffcharge;
  if (fGRF !=0 ) fGRF->Delete();
}

void AliTPCPRF2D::SetY(Float_t y1, Float_t y2, Int_t nYdiv)
{
  //
  //set virtual line position
  //first and last line and number of lines
  fNYdiv = nYdiv;
  if (ffcharge!=0) delete [] ffcharge;
  ffcharge = new Float_t[fNPRF*fNYdiv];
  fY1=y1;
  fY2=y2;
}

void AliTPCPRF2D::SetPad(Float_t width, Float_t height)
{
  //set base chevron parameters
 fHeightFull=height;
 fWidth=width;
}
void AliTPCPRF2D::SetChevron(Float_t hstep, 
			Float_t shifty, 
			Float_t fac)
{
  //set shaping of chewron parameters
  fHeightS=hstep;
  fShiftY=shifty;
  fK=fWidth*fac/hstep;
}

void AliTPCPRF2D::SetChParam(Float_t width, Float_t height,
		  Float_t hstep, Float_t shifty, Float_t fac)
{
  SetPad(width,height);
  SetChevron(hstep,shifty,fac);
}


Float_t AliTPCPRF2D::GetPRF(Float_t xin, Float_t yin, Bool_t inter)
{
  if (ffcharge==0) return 0;
  //  Float_t y=Float_t(fNYdiv-1)*(yin-fY1)/(fY2-fY1);
  //transform position to "wire position"
  Float_t y=fDYtoWire*(yin-fY1);
  if (fNYdiv == 1) y=fY1;
  //normaly it find nearest line charge
  if (inter ==kFALSE){   
    Int_t i=Int_t(0.5+y);
    if (y<0) i=Int_t(-0.5+y);
    if ((i<0) || (i>=fNYdiv) ) return 0;
    fcharge   = &(ffcharge[i*fNPRF]);
    return GetPRFActiv(xin);
  }
  else{
    //make interpolation from more fore lines
    Int_t i= Int_t(y);
    if ((i<0) || (i>=fNYdiv) ) return 0;
    Float_t z0=0;
    Float_t z1=0;
    Float_t z2=0;
    Float_t z3=0;
    if (i>0) {
      fcharge =&(ffcharge[(i-1)*fNPRF]);
      z0 = GetPRFActiv(xin);
    }
    fcharge =&(ffcharge[i*fNPRF]);
    z1=GetPRFActiv(xin);
    if ((i+1)<fNYdiv){
      fcharge =&(ffcharge[(i+1)*fNPRF]);
      z2 = GetPRFActiv(xin);
    }
    if ((i+2)<fNYdiv){
      fcharge =&(ffcharge[(i+2)*fNPRF]);
      z3 = GetPRFActiv(xin);
    }
    Float_t a,b,c,d,K,L;
    a=z1;
    b=(z2-z0)/2.;
    K=z2-a-b;
    L=(z3-z1)/2.-b;
    d=L-2*K;
    c=K-d;
    Float_t dy=y-Float_t(i);
        Float_t res = a+b*dy+c*dy*dy+d*dy*dy*dy;  
	//Float_t res = z1*(1-dy)+z2*dy;
    return res;            
  }        
} 


Float_t AliTPCPRF2D::GetPRFActiv(Float_t xin)
{
  //x xin DStep unit
  //return splaine aproximaton 
  Float_t x = (xin*fDStepM1)+fNPRF/2;
  Int_t i = Int_t(x);
  
  if  ( (i>0) && ((i+2)<fNPRF)) {
    Float_t a,b,c,d,K,L;
    a = fcharge[i];
    b = (fcharge[i+1]-fcharge[i-1])*0.5; 
    K = fcharge[i+1]-a-b;
    L = (fcharge[i+2]-fcharge[i])*0.5-b;
    d=L-2.*K;
    c=K-d;
    Float_t dx=x-Float_t(i);
    Float_t res = a+b*dx+c*dx*dx+d*dx*dx*dx;  
    return res;
  }
  else return 0;
}


Float_t  AliTPCPRF2D::GetGRF(Float_t xin, Float_t yin)
{  
  if (fGRF != 0 ) 
    return fkNorm*fGRF->Eval(xin,yin)/fInteg;
      else
    return 0.;
}

   
void AliTPCPRF2D::SetParam( TF2 * GRF,  Float_t kNorm, 
		       Float_t sigmaX, Float_t sigmaY)
{
   if (fGRF !=0 ) fGRF->Delete();
   fGRF = GRF;
   fkNorm = kNorm;
   if (sigmaX ==0) sigmaX=(fWidth+fK*fHeightS)/sqrt12;
   if (sigmaY ==0) sigmaY=(fWidth+fK*fHeightS)/sqrt12;
   forigsigmaX=sigmaX; 
   forigsigmaY=sigmaY; 
   fDStep = TMath::Sqrt(sigmaX*sigmaX+fWidth*fWidth/6.)/10.; 
   //   Update();   
  sprintf(fType,"User");
}
  

void AliTPCPRF2D::SetGauss(Float_t sigmaX, Float_t sigmaY,
		      Float_t kNorm)
{
  fkNorm = kNorm;
  if (fGRF !=0 ) fGRF->Delete();
  fGRF = new TF2("fun",funGauss2D,-5.,5.,-5.,5.,4);
  funParam[0]=sigmaX;
  funParam[1]=sigmaY;  
  funParam[2]=fK;
  funParam[3]=fHeightS;    
  forigsigmaX=sigmaX;
  forigsigmaY=sigmaY;
  fGRF->SetParameters(funParam);
  fDStep = TMath::Sqrt(sigmaX*sigmaX+fWidth*fWidth/6.)/10.; 
  //by default I set the step as one tenth of sigma
  //Update();
  sprintf(fType,"Gauss");
}

void AliTPCPRF2D::SetCosh(Float_t sigmaX, Float_t sigmaY,
		     Float_t kNorm)
{
  fkNorm = kNorm;
  if (fGRF !=0 ) fGRF->Delete();
  fGRF = new TF2("fun",	funCosh2D,-5.,5.,-5.,5.,4);   
  funParam[0]=sigmaX;
  funParam[1]=sigmaY;
  funParam[2]=fK;  
  funParam[3]=fHeightS;
  fGRF->SetParameters(funParam);
  forigsigmaX=sigmaX;
  forigsigmaY=sigmaY;
  fDStep = TMath::Sqrt(sigmaX*sigmaX+fWidth*fWidth/6.)/10.; 
  //by default I set the step as one tenth of sigma
  //Update();
  sprintf(fType,"Cosh");
}

void AliTPCPRF2D::SetGati(Float_t K3X, Float_t K3Y,
		     Float_t padDistance,
		     Float_t kNorm)
{
  fkNorm = kNorm;
  if (fGRF !=0 ) fGRF->Delete();
  fGRF = new TF2("fun",	funGati2D,-5.,5.,-5.,5.,5);  
  fK3X=K3X;
  fK3Y=K3Y;
  fPadDistance=padDistance;
  funParam[0]=padDistance;
  funParam[1]=K3X;
  funParam[2]=fK;  
  funParam[3]=fHeightS;
  funParam[4]=K3Y;
  fGRF->SetParameters(funParam);
  forigsigmaX=padDistance;
  forigsigmaY=padDistance;
  fDStep = TMath::Sqrt(padDistance*padDistance+fWidth*fWidth/6.)/10.; 
  //by default I set the step as one tenth of sigma
  //Update();
  sprintf(fType,"Gati");
}



void AliTPCPRF2D::Update()
{
  for (Int_t i=0; i<fNYdiv; i++){
    if (fNYdiv == 1) fActualY = fY1;
    else
      fActualY = fY1+Float_t(i)*(fY2-fY1)/Float_t(fNYdiv-1);
    fcharge   = &(ffcharge[i*fNPRF]);
    Update1();
  }
}



void AliTPCPRF2D::Update1()
{
  //initialize to 0
  

  Int_t i;
  Float_t x;
  for (i =0; i<fNPRF;i++)  fcharge[i] = 0;
  if ( fGRF == 0 ) return;
  ////////////////////////////////////////////////////////
  //I'm waiting for normal integral
  //in this moment only sum
  Float_t x2=  4*forigsigmaX;
  Float_t y2=  4*forigsigmaY;
  Float_t dx = forigsigmaX/Float_t(fNdiv*6);
  Float_t dy = forigsigmaY/Float_t(fNdiv*6);  
  fInteg  = 0;
  for (x=0.;x<x2;x+=dx)
    for (Float_t y=0;y<y2;y+=dy) fInteg+=fGRF->Eval(x,y)*dx*dy;
  fInteg*=4;
  /////////////////////////////////////////////////////
      
  
  if ( fInteg == 0 ) fInteg = 1; 
  
    //integrate charge over pad for different distance of pad
    for (i =0; i<fNPRF;i++)
      {      //x in cm fWidth in cm
	//calculate integral 
	Float_t xch = fDStep * (Float_t)(i-fNPRF/2);
	Float_t k=1;
	fcharge[i]=0;
	for (Float_t y=-fHeightFull/2.-fShiftY;
	     y<fHeightFull/2.;y+=fHeightS){
	  Float_t y2=TMath::Min((y+fHeightS),Float_t(fHeightFull/2.));
	  Float_t y1=TMath::Max((y),Float_t(-fHeightFull/2.));
	  Float_t x1;
	
	  if (k>0) 
	    x1 = (y2-y1)*fK-(fWidth+fK*fHeightS)/2.;	  
	  else
	    x1 =-(fWidth+fK*fHeightS)/2. ;	  
	  Float_t x2=x1+fWidth;

	  if (y2>y1) {
	    
            if ((x2-x1)*fNdiv<forigsigmaX) dx=(x2-x1);
	    else{
	      dx= forigsigmaX/Float_t(fNdiv);
	      dx = (x2-x1)/Float_t(Int_t(3+(x2-x1)/dx));	  
	    }	    
	    Float_t dy;
	    if ((y2-y1)*fNdiv<forigsigmaY) dy=(y2-y1);
	    else{	      
	      dy= forigsigmaY/Float_t(fNdiv);
	      dy = (y2-y1)/Float_t(Int_t(3+(y2-y1)/dy));
	    }

	    for (x=x1;x<x2;x+=dx)
	      for (Float_t y=y1;y<y2;y+=dy){
		if ( (y>(fActualY-(4.0*forigsigmaY))) &&
		     (y<(fActualY+(4.0*forigsigmaY)))){
		  Float_t xt=x-k*fK*(y-y1); 
		  if ((TMath::Abs(xch-xt)<4*forigsigmaX)){
		    
		    Float_t z0=fGRF->Eval(xch-(xt+dx/2.),fActualY-(y+dy/2.));
		    
		    Float_t z1=fGRF->Eval(xch-(xt+dx/2.),fActualY-y);
		    Float_t z2=fGRF->Eval(xch-xt,fActualY-(y+dy/2.));
		    Float_t z3=fGRF->Eval(xch-(xt-dx/2.),fActualY-y);
		    Float_t z4=fGRF->Eval(xch-xt,fActualY-(y-dy/2.));
		    if (z0<0) z0=0;
		    if (z1<0) z1=0;
		    if (z2<0) z2=0;
		    if (z3<0) z3=0;
		    if (z4<0) z4=0;
		    
		    //	      Float_t a=(z1-z3)/2;
		    //	      Float_t b=(z2-z4)/2;
		    Float_t c= (z3+z1-2*z0)/2.;
		    Float_t d= (z2+z4-2*z0)/2.;
		    Float_t z= (z0+c/12.+d/12.);	      	      		
		    
		    //Float_t z= fGRF->Eval(xch-xt,fActualY-y);
		    if (z>0.)	      fcharge[i]+=z*dx*dy/fInteg;	      
		  }
		}
	      }
	  }
	  k*=-1;
	}
      };   
  
  fSigmaX = 0; 
  Float_t sum =0;
  Float_t mean=0;
  for (x =-fNPRF*fDStep; x<fNPRF*fDStep;x+=fDStep)
    {      //x in cm fWidth in cm
      Float_t weight = GetPRFActiv(x);
      fSigmaX+=x*x*weight; 
      mean+=x*weight;
      sum+=weight;
    };  
  if (sum>0){
    mean/=sum;
    fSigmaX = TMath::Sqrt(fSigmaX/sum-mean*mean);   
  }
  else fSigmaX=0; 
  //calculate conversion coefitient to convert position to virtual wire
  fDYtoWire=Float_t(fNYdiv-1)/(fY2-fY1);
  fDStepM1=1/fDStep;
}

void AliTPCPRF2D::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliTPCPRF2D

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TObject::Streamer(R__b);     
      //read chewron parameters
      R__b >> fSigmaX;
      R__b >> fHeightFull;
      R__b >> fHeightS;
      R__b >> fShiftY;
      R__b >> fWidth;
      R__b >> fK;
      R__b >> fActualY;
      //read charge parameters
      R__b >> fType[0];
      R__b >> fType[1];
      R__b >> fType[2];
      R__b >> fType[3];
      R__b >> fType[4];
      R__b >> forigsigmaX;
      R__b >> forigsigmaY;
      R__b >> fkNorm;
      R__b >> fK3X;
      R__b >> fK3Y;
      R__b >> fPadDistance;
      R__b >> fInteg;
      
      //read functions
      if (fGRF!=0) { 
	delete [] fGRF;  
	fGRF=0;
      }
      if (strncmp(fType,"User",3)==0){
	fGRF= new TF2;
	R__b>>fGRF;   
      }
      if (strncmp(fType,"Gauss",3)==0) 
	fGRF = new TF2("fun",funGauss2D,-5.,5.,-5.,5.,4);
      if (strncmp(fType,"Cosh",3)==0) 
	fGRF = new TF2("fun",funCosh2D,-5.,5.,-5.,5.,4);
       if (strncmp(fType,"Gati",3)==0) 
	fGRF = new TF2("fun",funGati2D,-5.,5.,-5.,5.,5);
      
      //read interpolation parameters
      R__b >>fY1;
      R__b >>fY2;
      R__b >>fNYdiv;  
      R__b >>fDStep;  
      R__b >>fNPRF;
      if (ffcharge!=0) delete [] ffcharge;
      ffcharge = new Float_t[fNPRF*fNYdiv];
      R__b.ReadFastArray(ffcharge,fNPRF*fNYdiv); 
      R__b.ReadFastArray(funParam,5); 
      if (fGRF!=0) fGRF->SetParameters(funParam);
      //calculate conversion coefitient to convert position to virtual wire
      fDYtoWire=Float_t(fNYdiv-1)/(fY2-fY1);
      fDStepM1=1/fDStep;
   } else {
      R__b.WriteVersion(AliTPCPRF2D::IsA());
      TObject::Streamer(R__b);      
      //write chewron parameters
      R__b << fSigmaX;
      R__b << fHeightFull;
      R__b << fHeightS;
      R__b << fShiftY;
      R__b << fWidth;
      R__b << fK;
      R__b << fActualY;
      //write charge parameters
      R__b << fType[0];
      R__b << fType[1];
      R__b << fType[2];
      R__b << fType[3];
      R__b << fType[4];

      R__b << forigsigmaX;
      R__b << forigsigmaY;
      R__b << fkNorm;
      R__b << fK3X;
      R__b << fK3Y;
      R__b << fPadDistance;  
      R__b << fInteg;

      if (strncmp(fType,"User",3)==0)	R__b <<fGRF;         
      //write interpolation parameters
      R__b <<fY1;
      R__b <<fY2;
      R__b <<fNYdiv;   
      R__b <<fDStep;
      R__b <<fNPRF;    
      R__b.WriteFastArray(ffcharge,fNPRF*fNYdiv); 
      R__b.WriteFastArray(funParam,5); 
   }
}




void AliTPCPRF2D::DrawX(Float_t x1 ,Float_t x2,Float_t y, Bool_t inter)
{ 
  if (fGRF==0) return ;
  const Int_t N=100;
  char s[100];
  TCanvas  * c1 = new TCanvas("canPRF","Pad response function",700,900);
  c1->cd();
  TPad * pad1 = new TPad("pad1PRF","",0.05,0.61,0.95,0.97,21);
  pad1->Draw();
  TPad * pad2 = new TPad("pad2PRF","",0.05,0.22,0.95,0.60,21);
  pad2->Draw();

  //  pad1->cd();  
  //pad2->cd();
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0); 
  sprintf(s,"PRF response function for chevron pad");  
  TH1F * hPRFc = new TH1F("hPRFc",s,N+1,x1,x2);
  Float_t x=x1;
  Float_t y1;
  //  Float_t y2;

  for (Float_t i = 0;i<N+1;i++)
    {
      x+=(x2-x1)/Float_t(N);
      y1 = GetPRF(x,y,inter);
      hPRFc->Fill(x,y1);
    };

  pad1->cd();
  fGRF->SetRange(x1,x1,x2,x2); 
  fGRF->SetNpx(25);
  fGRF->SetNpy(25); 
  fGRF->Draw("lego2");
  // hPRFo->Fit("gaus");
  gStyle->SetOptStat(1); 
  pad2->cd();
  hPRFc->Fit("gaus");
  c1->cd(); 
  TPaveText * comment = new TPaveText(0.05,0.02,0.95,0.20,"NDC");
  comment->SetTextAlign(12);
  comment->SetFillColor(42);
  TText *title = comment->AddText("Chevron pad parameters:");
  title->SetTextSize(0.03);
  sprintf(s,"Full height of pad:  %2.2f",fHeightFull);
  comment->AddText(s);
  sprintf(s,"Height of one chevron unit h:  %2.2f cm",2*fHeightS);
  comment->AddText(s);
  sprintf(s,"Width of one chevron unit  w:  %2.2f cm",fWidth);
  comment->AddText(s);
  sprintf(s,"Overlap factor:  %2.2f",fK*fHeightS/fWidth);
  comment->AddText(s);
  sprintf(s,"Y position:  %2.2f ",y);
  comment->AddText(s);
  sprintf(s,"Sigma x of original distribution: %2.2f ",forigsigmaX);
  comment->AddText(s);  
  sprintf(s,"Sigma y of original distribution: %2.2f ",forigsigmaY);
  comment->AddText(s);    
  sprintf(s,"Type of original distribution: %s ",fType);
  comment->AddText(s); 
  comment->Draw();
}



void AliTPCPRF2D::Draw(Float_t x1 ,Float_t x2,Float_t y1, Float_t y2, 
		  Bool_t inter, Int_t Nx, Int_t Ny)
{ 
  char s[100];
  if (fGRF==0) return ;
  TCanvas  * c1 = new TCanvas("canPRF","Pad response function",700,900);
  c1->cd();
  TPad * pad1 = new TPad("pad1PRF","",0.05,0.61,0.95,0.97,21);
  pad1->Draw();
  TPad * pad2 = new TPad("pad2PRF","",0.05,0.22,0.95,0.60,21);
  pad2->Draw();

  //  pad1->cd();  
  //pad2->cd();
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0); 
  sprintf(s,"PRF response function for chevron pad");  
  TH2F * hPRFc = new TH2F("hPRFc",s,Nx+1,x1,x2,Ny+1,y1,y2);
  Float_t dx=(x2-x1)/Float_t(Nx);
  Float_t dy=(y2-y1)/Float_t(Ny) ;
  Float_t x,y,z;
  //  Float_t y2;
  for ( x = x1;x<=x2;x+=dx){
    for(y = y1;y<=y2;y+=dy)
      {
	z = GetPRF(x,y,inter);
	hPRFc->Fill(x,y,z);
      };
  }
  pad1->cd();
  fGRF->SetRange(x1,y1,x2,y2); 
  fGRF->SetNpx(25);
  fGRF->SetNpy(25); 
  fGRF->Draw("lego2");
  // hPRFo->Fit("gaus");
  gStyle->SetOptStat(1); 
  pad2->cd();
  hPRFc->Draw("lego2");
  c1->cd(); 
  TPaveText * comment = new TPaveText(0.05,0.02,0.95,0.20,"NDC");
  comment->SetTextAlign(12);
  comment->SetFillColor(42);
  TText *title = comment->AddText("Chevron pad parameters:");
  title->SetTextSize(0.03);
  sprintf(s,"Full height of pad:  %2.2f",fHeightFull);
  comment->AddText(s);
  sprintf(s,"Height of one chevron unit h:  %2.2f cm",2*fHeightS);
  comment->AddText(s);
  sprintf(s,"Width of one chevron unit  w:  %2.2f cm",fWidth);
  comment->AddText(s);
  sprintf(s,"Overlap factor:  %2.2f",fK*fHeightS/fWidth);
  comment->AddText(s); 
  sprintf(s,"Sigma x of original distribution: %2.2f ",forigsigmaX);
  comment->AddText(s);  
  sprintf(s,"Sigma y of original distribution: %2.2f ",forigsigmaY);
  comment->AddText(s);    
  sprintf(s,"Type of original distribution: %s ",fType);
  comment->AddText(s); 
  comment->Draw();
}

void AliTPCPRF2D::DrawDist(Float_t x1 ,Float_t x2,Float_t y1, Float_t y2, 
		  Bool_t inter, Int_t Nx, Int_t Ny, Float_t thr)
{ 
  const Float_t minth=0.00001;
  if (thr<minth) thr=minth;
  char s[100];
  if (fGRF==0) return ;
  TCanvas  * c1 = new TCanvas("padDistortion","COG distortion",700,900);
  c1->cd();
  TPad * pad1 = new TPad("CHARGE","",0.05,0.61,0.95,0.97,21);
  pad1->Draw();
  TPad * pad2 = new TPad("dist","",0.05,0.22,0.95,0.60,21);
  pad2->Draw();

  //  pad1->cd();  
  //pad2->cd();
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0); 
  sprintf(s,"COG distortion (threshold=%2.2f)",thr);  
  TH2F * hPRFDist = new TH2F("hDistortion",s,Nx+1,x1,x2,Ny+1,y1,y2);
  Float_t dx=(x2-x1)/Float_t(Nx);
  Float_t dy=(y2-y1)/Float_t(Ny) ;
  Float_t x,y,z,ddx;
  //  Float_t y2;
  for ( x = x1;x<(x2+dx/2.);x+=dx)
    for(y = y1;y<=(y2+dx/2.);y+=dy)
      {
	Float_t sumx=0;
	Float_t sum=0;
	for (Float_t padx=-fWidth;padx<(fWidth*1.1);padx+=fWidth)
	  {	    
	    z = GetPRF(x-padx,y,inter);
	    if (z>thr){
	      sum+=z;
	      sumx+=z*padx;
	    }	
	  };	
	if (sum>minth)  
	  {
	    ddx = (x-(sumx/sum));
	  }
	else ddx=-1;
	if (TMath::Abs(ddx)<10) 	hPRFDist->Fill(x,y,ddx);
      }
  pad1->cd();
  fGRF->SetRange(x1,y1,x2,y2); 
  fGRF->SetNpx(25);
  fGRF->SetNpy(25); 
  fGRF->Draw("lego2");
  // hPRFo->Fit("gaus");
  //  gStyle->SetOptStat(1); 
  pad2->cd();
  hPRFDist->Draw("lego2");
  
  c1->cd(); 
  TPaveText * comment = new TPaveText(0.05,0.02,0.95,0.20,"NDC");
  comment->SetTextAlign(12);
  comment->SetFillColor(42);
  //  TText *title = comment->AddText("Distortion of COG method");
  //  title->SetTextSize(0.03);
  TText * title = comment->AddText("Chevron pad parameters:");
  title->SetTextSize(0.03);
  sprintf(s,"Full height of pad:  %2.2f",fHeightFull);
  comment->AddText(s);
  sprintf(s,"Height of one chevron unit h:  %2.2f cm",2*fHeightS);
  comment->AddText(s);
  sprintf(s,"Width of one chevron unit  w:  %2.2f cm",fWidth);
  comment->AddText(s);
  sprintf(s,"Overlap factor:  %2.2f",fK*fHeightS/fWidth);
  comment->AddText(s); 
  sprintf(s,"Sigma x of original distribution: %2.2f ",forigsigmaX);
  comment->AddText(s);  
  sprintf(s,"Sigma y of original distribution: %2.2f ",forigsigmaY);
  comment->AddText(s);    
  sprintf(s,"Type of original distribution: %s ",fType);
  comment->AddText(s); 
  comment->Draw();
  
}

