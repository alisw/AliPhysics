#include "AliHBTFits.h"
//_________________________________________________
///////////////////////////////////////////////////////////////////////////////////
//
// class AliHBTFits
//
// Sets of methods for fittig correlation functions
//
//
//
// 
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////////////////////////////
#include "AliHBTAnalysis.h"

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TError.h>
#include <TDirectory.h>
#include <TMath.h>

ClassImp(AliHBTFits)

TF1*     AliHBTFits::fgF1 = 0x0;
TF1*     AliHBTFits::fgF2 = 0x0;
Double_t *AliHBTFits::fXbins = 0;
Double_t *AliHBTFits::fYbins = 0;
Double_t *AliHBTFits::fZbins = 0;

Double_t *AliHBTFits::fExpT = 0;
Double_t *AliHBTFits::fLongT = 0;
Double_t *AliHBTFits::fIntegrT = 0;

Int_t  AliHBTFits::fNX = 0;
Int_t  AliHBTFits::fNY = 0;
Int_t  AliHBTFits::fNZ = 0;
 Int_t AliHBTFits::fNXY = 0;
Int_t  AliHBTFits::fNT = 0;
Int_t  AliHBTFits::fBinX = 0;
Int_t  AliHBTFits::fBinY = 0;
Int_t  AliHBTFits::fBinZ = 0;
  
Bool_t AliHBTFits::fLongOnly = 1;

/******************************************************************************************/
AliHBTFits::~AliHBTFits()
{
 //dtor
 delete fgF1;
 delete fgF2;
}
/******************************************************************************************/

/******************************************************************************************/
/********    Qout-Qside-Qlong  Cyl Surf   *************************************************/
/******************************************************************************************/

void AliHBTFits::FitQOutQSideQLongCylSurf(const TString& hname,Option_t* fopt, Float_t xmax, Float_t ymax,Float_t zmax)
{
  //Fits 3D histo with QOutQSideQLongCylSurf
 TH3* h = dynamic_cast<TH3*>(gDirectory->Get(hname));
 if (h == 0x0)
  {
    ::Error("FitQOutCylSurf","Can not find histogram %s",hname.Data());
    return;
  }
  
 delete gROOT->GetListOfFunctions()->FindObject("BorisTomasikCylSurfQOutQSideQLong");

 delete  fgF1;
 delete  fgF2;
 // here x is phi
 //par 0: R
 //par 1: Phi
 //par 2: qOut
 //par 3: qSide
 fgF1 = new TF1("ff1","sin([0]*([2]*cos(x)+[3]*sin(x))/0.197)*exp(cos(x)/[1])",-TMath::Pi(),TMath::Pi());
 fgF2 = new TF1("ff2","cos([0]*([2]*cos(x)+[3]*sin(x))/0.197)*exp(cos(x)/[1])",-TMath::Pi(),TMath::Pi());
 fgF1->SetNpx(1000);
 fgF2->SetNpx(1000);
 
 if (xmax <= 0.0) xmax = h->GetXaxis()->GetXmax();
 if (ymax <= 0.0) ymax = h->GetYaxis()->GetXmax();
 if (zmax <= 0.0) zmax = h->GetZaxis()->GetXmax();
 
 TF3 *theory = new TF3("BorisTomasikCylSurfQOutQSideQLong",&AliHBTFits::QOutQSideQLongCylSurf,0.0,xmax,0.0,ymax,0.0,zmax,5);//par :aLambda, aR, aJ,aPhi,
 theory->SetNpx(1000);
 
 theory->SetParameter(0,.75);
 theory->SetParameter(1,6.0);
 theory->SetParameter(2,1.0);
 theory->SetParameter(3,2.0);
 theory->SetParameter(4,0.9);
 
 theory->SetParName(0,"\\lambda");
 theory->SetParName(1,"R");
 theory->SetParName(2,"J");
 theory->SetParName(3,"L");
 theory->SetParName(4,"\\Phi");


/***********************************/
 /***********************************/
  /***********************************/

   Double_t xmin = 0.0;
   Double_t ymin = 0.0;
   Double_t zmin = 0.0;

   Int_t hxfirst = h->GetXaxis()->FindFixBin(xmin); if (hxfirst < 1) hxfirst = 1;
   Int_t hxlast  = h->GetXaxis()->FindFixBin(xmax); if (hxlast > h->GetXaxis()->GetLast()) hxlast = h->GetXaxis()->GetLast();
   Int_t hyfirst = h->GetYaxis()->FindFixBin(ymin); if (hyfirst < 1) hyfirst = 1;
   Int_t hylast  = h->GetYaxis()->FindFixBin(ymax); if (hylast > h->GetYaxis()->GetLast()) hylast = h->GetYaxis()->GetLast();
   Int_t hzfirst = h->GetZaxis()->FindFixBin(zmin); if (hzfirst < 1) hzfirst = 1;
   Int_t hzlast  = h->GetZaxis()->FindFixBin(zmax); if (hzlast > h->GetZaxis()->GetLast()) hzlast = h->GetZaxis()->GetLast();
   
   fNX = hxlast-hxfirst+1;
   fNY = hylast-hyfirst+1;
   fNZ = hzlast-hzfirst+1;
   
   fXbins = new Double_t[fNX];
   fYbins = new Double_t[fNY];
   fZbins = new Double_t[fNZ];
   Int_t c = 0;
   for (Int_t k=hzfirst; k<=hzlast; k++) 
    {
      if (h->GetZaxis()->GetBinCenter(k) > zmax) break;
      fZbins[c++] = TMath::Abs(h->GetZaxis()->GetBinCenter(k));
    }
   fNZ = c;
   c = 0;
   for (Int_t j=hyfirst; j<=hylast; j++)
    {
      if (h->GetYaxis()->GetBinCenter(j) > ymax) break;
      fYbins[c++] = TMath::Abs(h->GetYaxis()->GetBinCenter(j));
    }
   fNY = c; 
   c = 0;
   for (Int_t i=hxfirst; i<=hxlast; i++) 
    {
      
      if (h->GetXaxis()->GetBinCenter(i) > xmax) break;
      fXbins[c++] = TMath::Abs(h->GetXaxis()->GetBinCenter(i));
      ::Info("","%d %d %d %f",i,c-1,fNX,fXbins[c-1]);
    }
   fNX = c;
   
   fNXY=fNX*fNY;
   fNT=fNXY*fNZ;
   
   fExpT = new Double_t[fNT];
   fLongT = new Double_t[fNZ];
   fIntegrT = new Double_t[fNXY];
   
   fBinX = 0;
   fBinY = 0;
   fBinZ = 0;
  /***********************************/
 /***********************************/
/***********************************/
 
   h->Fit(theory,fopt,"E");
   
   delete [] fXbins;
   delete [] fYbins;
   delete [] fZbins;
   
   delete [] fExpT;
   delete [] fLongT;
   delete [] fIntegrT;
   
   fXbins = 0x0;
   fYbins = 0x0;
   fZbins = 0x0;
   fExpT = 0x0;
   fLongT = 0x0;
   fIntegrT = 0x0;

}

/******************************************************************************************/
Double_t AliHBTFits::QOutQSideQLongCylSurf(Double_t *x, Double_t *par)
{
  //QOutQSideQLongCylSurf Function
  
  if (fXbins == 0x0) 
   {
     ::Info("QOutQSideQLongCylSurf","arrays deleted");
     return 1.0;
   }

  
  Double_t qout = x[0];
  Double_t qside = x[1];
  Double_t qlong = x[2];

//  ::Info("QOutQSideQLongCylSurf","out=%f, side=%f, long=%f %d %d %d",qout,qside,qlong,fBinX,fBinY,fBinZ);
  
  if (fXbins[fBinX] != x[0])
   {
     ::Fatal("QOutQSideQLongCylSurf","expected X %f got %f",fXbins[fBinX],qout); 
     return 1.0;
   }
  

  if (fXbins[fBinY] != x[1])
   {
     ::Fatal("QOutQSideQLongCylSurf","expected Y %f got %f",fYbins[fBinY],qside); 
     return 1.0;
   }
  
  if (fXbins[fBinZ] != x[2])
   {
     ::Fatal("QOutQSideQLongCylSurf","expected Z %f got %f",fZbins[fBinZ],qlong); 
     return 1.0;
   }
  
  Double_t aLambda = par[0];
  Double_t aR = par[1];
  Double_t aJ = par[2];
  Double_t aL = par[3];
  Double_t aPhi = par[4];
  
  static Double_t oldLambda = -1;
  static Double_t oldR = -1;
  static Double_t oldJ = -1;
  static Double_t oldL = -1;
  static Double_t oldPhi = -1;
  
  
    
  if ( (oldLambda != aLambda) || (oldR != aR) || (oldL != aL))
   {
    ::Info("QOutQSideQLongCylSurf","out=%f, side=%f, long=%f\n lambda=%f, R=%f, J=%f, L=%f Phi=%f",
                                 qout,qside,qlong,aLambda,aR,aJ,aL,aPhi);
    oldLambda = aLambda;
   }

  Double_t result = aLambda;
  
  Double_t exppart; 
  
  if (aJ == oldJ)
   {
     exppart = fExpT[fBinZ*fNXY + fBinY*fNX + fBinX];
   }
  else
   {
     Int_t cbin = 0;
     Double_t aJ2=aJ*aJ/0.0389273;
     for (Int_t k=0; k<fNZ; k++)
      {
       Double_t qL2 = fZbins[k]*fZbins[k];
       
       for (Int_t j=0; j<fNY; j++)
        {
          Double_t qS2 = fYbins[j]*fYbins[j];
          
          for (Int_t i=0; i<fNX; i++)
           {
             fExpT[cbin++] = TMath::Exp(-(fXbins[i]*fXbins[i]+qS2+qL2)*aJ2);
           }
        }   
      }
     exppart = fExpT[0];
     oldJ = aJ;
   } 
   
  result = aLambda*exppart;
//  ::Info("QOSL","lambda*exp=%f",result);

  Double_t longpart = 1.0;
  
  if (oldL == aL)
   {
     longpart = fLongT[fBinZ];
   }
  else
   {
     if ( aL <= 0.0 )
      {
        for (Int_t i=0; i<fNZ; i++) fLongT[i] = 1.0;
      }
     else
      {
        for (Int_t i=0; i<fNZ; i++) 
         {  
           Double_t qlL = fZbins[i]*aL/(2.0*0.197);
           Double_t sinqlL = TMath::Sin(qlL);
           Double_t sin2qlL = sinqlL*sinqlL;
           fLongT[i] = sin2qlL/(qlL*qlL);
         }
      }
     longpart = fLongT[0];
     oldL = aL;
   }
   
//  ::Info("QOSL","longpart=%f",longpart);
 
  result *=longpart;
  Double_t itgrl;
  
  if ( (oldR != aR) || (oldPhi !=aPhi) )
   {
     if (aPhi <= 0.0) 
      {
       for (Int_t i=0; i<fNXY; i++) 
        {
          fIntegrT[i] = 1.0;
        }
      } 
     else
      {
       Int_t cbin = 0;
       for (Int_t j=0; j<fNY; j++)
        {
          for (Int_t i=0; i<fNX; i++)
           {
             
             fgF1->SetParameter(0,aR);
             fgF2->SetParameter(0,aR);

             fgF1->SetParameter(1,aPhi);
             fgF2->SetParameter(1,aPhi);

             fgF1->SetParameter(2,fXbins[i]);
             fgF2->SetParameter(2,fXbins[i]);

             fgF1->SetParameter(3,fYbins[j]);
             fgF2->SetParameter(3,fYbins[j]);

             Double_t sp = fgF1->Integral(-TMath::Pi(),TMath::Pi(),(Double_t*)0x0,1.e-8);
             Double_t cp = fgF2->Integral(-TMath::Pi(),TMath::Pi(),(Double_t*)0x0,1.e-8);

         //    ::Info("QSideCylSurf","sp=%f, cp=%f",cp,sp);

             fIntegrT[cbin] = cp*cp + sp*sp;
         //    ::Info("QSideCylSurf","2: result=%f",result);

             Double_t oneOverPhi = 1./aPhi;
             Double_t tmp = TMath::BesselI0(oneOverPhi)-TMath::StruveL0(oneOverPhi);
             tmp = tmp*TMath::TwoPi();
             tmp = tmp*tmp;
             if (tmp == 0.0)
              {
                ::Error("QOutCylSurf","Div by 0");
                return 1.0;
              }
             fIntegrT[cbin] /= tmp;
//             ::Info("integr calc","%d (%d %d) %f",cbin,i,j,fIntegrT[cbin]);
             cbin++;
           }
        }   
         
      }
     itgrl = fIntegrT[0];
     oldR = aR;
     oldPhi = aPhi;
   }
  else
   {
     itgrl = fIntegrT[fBinY*fNX+fBinX];
  //   ::Info("integr ret","%d (%d %d) %f",fBinY*fNX+fBinX,fBinX,fBinY,fIntegrT[fBinY*fNX+fBinX]);
   }
  
//  ::Info("QOSL","itgrl/bessel=%f",itgrl);

  result *= itgrl;
   
  fBinX++;
  
  if (fBinX == fNX) 
   {
     fBinX = 0;
     fBinY++;
     if ( fBinY == fNY )
      {
        fBinY = 0;
        fBinZ++;
        if ( fBinZ == fNZ )
         {
           fBinZ = 0;
         }
      }
   }

//  ::Info("QOutQSideQLongCylSurf Return","out=%f, side=%f, long=%f %d %d %d",qout,qside,qlong,fBinX,fBinY,fBinZ);
  
  result+=1.; 
  
//  ::Info("QSideCylSurf","Final: result=%f",result);
//  ::Info("QSideCylSurf","XXX result=%f",XXX(x,par));
//  AliHBTAnalysis::PressAnyKey();
  return result;
}
/******************************************************************************************/
/********    Qout Qside Cyl Surf   ********************************************************/
/******************************************************************************************/

void AliHBTFits::FitQOutQSideCylSurf(const TString& hname,Option_t* fopt,Float_t xmax, Float_t ymax)
{
  //Fits 3D histo with QOutQSideQLongCylSurf
 TH2* h = dynamic_cast<TH2*>(gDirectory->Get(hname));
 if (h == 0x0)
  {
    ::Error("FitQOutCylSurf","Can not find histogram %s",hname.Data());
    return;
  }
  
 delete gROOT->GetListOfFunctions()->FindObject("BorisTomasikCylSurfQOutQSide");

 delete  fgF1;
 delete  fgF2;
 // here x is phi
 //par 0: R
 //par 1: Phi
 //par 2: qOut
 //par 3: qSide
 fgF1 = new TF1("ff1","sin([0]*([2]*cos(x)+[3]*sin(x))/0.197)*exp(cos(x)/[1])",-TMath::Pi(),TMath::Pi());
 fgF2 = new TF1("ff2","cos([0]*([2]*cos(x)+[3]*sin(x))/0.197)*exp(cos(x)/[1])",-TMath::Pi(),TMath::Pi());
 fgF1->SetNpx(1000);
 fgF2->SetNpx(1000);
 
 if (xmax <= 0.0) xmax = h->GetXaxis()->GetXmax();
 if (ymax <= 0.0) ymax = h->GetYaxis()->GetXmax();
 
 TF2 *theory = new TF2("BorisTomasikCylSurfQOutQSide",&AliHBTFits::QOutQSideCylSurf,0.0,xmax,0.0,ymax,4);//par :aLambda, aR, aJ,aPhi,
 theory->SetNpx(1000);
 
 theory->SetParameter(0,.75);
 theory->SetParameter(1,6.0);
 theory->SetParameter(2,1.0);
 theory->SetParameter(3,0.9);
 
 theory->SetParName(0,"\\lambda");
 theory->SetParName(1,"R");
 theory->SetParName(2,"J");
 theory->SetParName(3,"\\Phi");
 
 h->Fit(theory,fopt,"E");
  
}
/******************************************************************************************/

Double_t AliHBTFits::QOutQSideCylSurf(Double_t *x, Double_t *par)
{
  //QOutQSideQLongCylSurf Function
  Double_t qout = x[0];
  Double_t qside = x[1];
  
  Double_t aLambda = par[0];
  Double_t aR = par[1];
  Double_t aJ = par[2];
  Double_t aPhi = par[3];
  
  static Double_t oldlam = aLambda;
  
  if (oldlam != aLambda)
   {
    ::Info("QOutQSideCylSurf","out=%f, side=%f, lambda=%f, R=%f, J=%f, Phi=%f",
                                 qout,qside,aLambda,aR,aJ,aPhi);
    oldlam = aLambda;
   }
  
  Double_t result=aLambda*TMath::Exp(-(qout*qout+qside*qside)*aJ*aJ/0.0389273);
  
  if (aPhi > 0.0)
   {
    fgF1->SetParameter(0,aR);
    fgF2->SetParameter(0,aR);
  
    fgF1->SetParameter(1,aPhi);
    fgF2->SetParameter(1,aPhi);
  
    fgF1->SetParameter(2,qout);
    fgF2->SetParameter(2,qout);

    fgF1->SetParameter(3,qside);
    fgF2->SetParameter(3,qside);
    
    Double_t sp = fgF1->Integral(-TMath::Pi(),TMath::Pi(),(Double_t*)0x0,1.e-8);
    Double_t cp = fgF2->Integral(-TMath::Pi(),TMath::Pi(),(Double_t*)0x0,1.e-8);
    
//    ::Info("QSideCylSurf","sp=%f, cp=%f",cp,sp);
    
    result *= cp*cp + sp*sp;
//    ::Info("QSideCylSurf","2: result=%f",result);

    Double_t oneOverPhi = 1./aPhi;
    Double_t tmp = TMath::BesselI0(oneOverPhi)-TMath::StruveL0(oneOverPhi);
    tmp = tmp*TMath::TwoPi();
    tmp = tmp*tmp;
    if (tmp == 0.0)
     {
       ::Error("QOutCylSurf","Div by 0");
       return 1.0;
     }
    result=result/tmp;
   }
   
  result+=1.; 

//  ::Info("QSideCylSurf","Final: result=%f",result);
  return result;

}
/******************************************************************************************/
/********    Qout  Cyl Surf   *************************************************************/
/******************************************************************************************/

void AliHBTFits::FitQOutCylSurf(const TString& hname, Option_t* fopt, Float_t max)
{
 TH1D* h = dynamic_cast<TH1D*>(gDirectory->Get(hname));
 if (h == 0x0)
  {
    ::Error("FitQOutCylSurf","Can not find histogram %s",hname.Data());
    return;
  }
  
 delete gROOT->GetListOfFunctions()->FindObject("BorisTomasikCylSurfQOut");

 delete  fgF1;
 delete  fgF2;
 fgF1 = new TF1("ff1","sin([0]*[2]*cos(x)/0.197)*exp(cos(x)/[1])",-TMath::Pi(),TMath::Pi());
 fgF2 = new TF1("ff2","cos([0]*[2]*cos(x)/0.197)*exp(cos(x)/[1])",-TMath::Pi(),TMath::Pi());
 fgF1->SetNpx(1000);
 fgF2->SetNpx(1000);
 
 if (max <= 0.0) max = h->GetXaxis()->GetXmax();
 TF1 *theory = new TF1("BorisTomasikCylSurfQOut",&AliHBTFits::QOutCylSurf,0.0,max,3);//par : aR, aPhi
 theory->SetNpx(1000);
 
 theory->SetParameter(0,6.0);
 theory->SetParameter(1,1.0);
 theory->SetParameter(2,0.75);

 theory->FixParameter(0,5.92);
 theory->FixParameter(1,0.864);
 theory->FixParameter(2,0.773);

 theory->SetParName(0,"R_{out}");
 theory->SetParName(1,"\\Phi");
 theory->SetParName(2,"\\lambda");
 
 h->Fit(theory,fopt,"E");
 
}
/******************************************************************************************/

Double_t AliHBTFits::QOutCylSurf(Double_t *x, Double_t *par)
{
  //Qout 
  Double_t qout = x[0];
  Double_t aR = par[0];
  Double_t aPhi = par[1];
  Double_t aLambda = par[2];

//  ::Info("QOutCylSurf","q=%f, lambda=%f, Rout=%f, Phi=%f",qout,aLambda,aR,aPhi);
  static const Double_t j2=1.063*1.063;
  Double_t result=aLambda*TMath::Exp(-qout*qout*j2/0.0389273);
  
//  ::Info("QOutCylSurf","1: result=%f",result);
  
  if (aPhi > 0.0)
   {
    fgF1->SetParameter(0,aR);
    fgF2->SetParameter(0,aR);
  
    fgF1->SetParameter(1,aPhi);
    fgF2->SetParameter(1,aPhi);
  
    fgF1->SetParameter(2,qout);
    fgF2->SetParameter(2,qout);
    
    Double_t sp = fgF1->Integral(-TMath::Pi(),TMath::Pi(),(Double_t*)0x0,1.e-8);
    Double_t cp = fgF2->Integral(-TMath::Pi(),TMath::Pi(),(Double_t*)0x0,1.e-8);
    
//    ::Info("QOutCylSurf","sp=%f, cp=%f",cp,sp);
    
    result = result * (cp*cp + sp*sp);
//    ::Info("QOutCylSurf","2: result=%f",result);

    Double_t oneOverPhi = 1./aPhi;
    Double_t tmp = TMath::BesselI0(oneOverPhi)-TMath::StruveL0(oneOverPhi);
    tmp = tmp*TMath::TwoPi();
    tmp = tmp*tmp;
    if (tmp == 0.0)
     {
       ::Error("QOutCylSurf","Div by 0");
       return 1.0;
     }
    result=result/tmp;
   }
   
  result+=1.; 

//  ::Info("QOutCylSurf","Final: result=%f",result);
  return result;
}
/******************************************************************************************/
/********    Qside  Cyl Surf   ************************************************************/
/******************************************************************************************/

void AliHBTFits::FitQSideCylSurf(const TString& hname, Option_t* fopt, Float_t max)
{
//Fits QSide According to Boris Tomasik formula
 TH1D* h = dynamic_cast<TH1D*>(gDirectory->Get(hname));
 if (h == 0x0)
  {
    ::Error("FitQSideCylSurf","Can not find histogram %s",hname.Data());
    return;
  }
 delete gROOT->GetListOfFunctions()->FindObject("BorisTomasikCylSurfQSide");
 delete fgF1;
 delete fgF2;
 fgF1 = new TF1("ff1","sin([0]*[2]*sin(x)/0.197)*exp(cos(x)/[1])",-TMath::Pi(),TMath::Pi());
 fgF2 = new TF1("ff2","cos([0]*[2]*sin(x)/0.197)*exp(cos(x)/[1])",-TMath::Pi(),TMath::Pi());
 fgF1->SetNpx(1000);
 fgF2->SetNpx(1000);
 
 if (max <= 0.0) max = h->GetXaxis()->GetXmax();
 TF1 *theory = new TF1("BorisTomasikCylSurfQSide",&AliHBTFits::QSideCylSurf,0.0,max,3);//par : aR, aPhi
 theory->SetNpx(1000);
 theory->SetParameter(0,6.0);
 theory->SetParameter(1,1.0);
 theory->SetParameter(2,0.75);
 
 theory->FixParameter(0,5.92);
 theory->FixParameter(1,0.864);
 theory->FixParameter(2,0.773);
 
 theory->SetParName(0,"R_{out}");
 theory->SetParName(1,"\\Phi");
 theory->SetParName(2,"\\lambda");
 
 h->Fit(theory,fopt,"E");
 
}
/******************************************************************************************/

Double_t AliHBTFits::QSideCylSurf(Double_t *x, Double_t *par)
{
  //Qout 
  Double_t qside = x[0];
  Double_t aR = par[0];
  Double_t aPhi = par[1];
  Double_t aLambda = par[2];

//  ::Info("QSideCylSurf","q=%f, lambda=%f, Rside=%f, Phi=%f",qside,aLambda,aR,aPhi);

  static const Double_t j2=1.063*1.063;
  Double_t result=aLambda*TMath::Exp(-qside*qside*j2/0.0389273);
  if (aPhi > 0.0)
   {
    fgF1->SetParameter(0,aR);
    fgF2->SetParameter(0,aR);
  
    fgF1->SetParameter(1,aPhi);
    fgF2->SetParameter(1,aPhi);
  
    fgF1->SetParameter(2,qside);
    fgF2->SetParameter(2,qside);
    
    Double_t sp = fgF1->Integral(-TMath::Pi(),TMath::Pi(),(Double_t*)0x0,1.e-7);
    Double_t cp = fgF2->Integral(-TMath::Pi(),TMath::Pi(),(Double_t*)0x0,1.e-7);

    result *= cp*cp + sp*sp;

    Double_t oneOverPhi = 1./aPhi;
    Double_t tmp = TMath::BesselI0(oneOverPhi)-TMath::StruveL0(oneOverPhi);
    tmp = tmp*TMath::TwoPi();
    tmp = tmp*tmp;
    if (tmp == 0.0)
     {
       ::Error("QOutCylSurf","Div by 0");
       return 1.0;
     }
     
    result/=tmp;
   }
  return result + 1.;
}
/******************************************************************************************/

void AliHBTFits::FitQLongCylSurf(const TString& hname, Option_t* fopt, Float_t max)
{
 //Fits QSide According to Boris Tomasik formula
 TH1D* h = dynamic_cast<TH1D*>(gDirectory->Get(hname));
 if (h == 0x0)
  {
    ::Error("FitQLongCylSurf","Can not find histogram %s",hname.Data());
    return;
  }
  
 delete fgF1;
 delete fgF2;
 fgF1 = new TF1("ff1","sin([0]*[2]*sin(x)/0.197)*exp(cos(x)/[1])",-TMath::Pi(),TMath::Pi());
 fgF2 = new TF1("ff2","cos([0]*[2]*sin(x)/0.197)*exp(cos(x)/[1])",-TMath::Pi(),TMath::Pi());
 fgF1->SetNpx(1000);
 fgF2->SetNpx(1000);
  
 delete gROOT->GetListOfFunctions()->FindObject("FitQLongCylSurf");
 if (max <= 0.0) max = h->GetXaxis()->GetXmax();
 TF1 *theory = new TF1("BorisTomasikCylSurfQLong",&AliHBTFits::QLongCylSurf,0.0,max,3);//par : aL, aPhi,aLambda
 theory->SetNpx(1000);
 theory->SetParameter(0,2.0);//L
 theory->SetParameter(1,1.0);//Phi
 theory->SetParameter(2,0.75);//Lambda
 
 theory->FixParameter(0,1.81);//L
 theory->FixParameter(1,0.864);//Phi
 theory->FixParameter(2,0.773);//Lambda
 
 theory->SetParName(0,"R_{long}");
 theory->SetParName(1,"\\Phi");
 theory->SetParName(2,"\\lambda");

 
 h->Fit(theory,fopt,"E");

}

Double_t AliHBTFits::QLongCylSurf(Double_t *x, Double_t *par)
{
  //Qout 
  Double_t qlong = x[0];
  Double_t aL = par[0];
  Double_t aPhi = par[1];
  Double_t aLambda = par[2];

  static const Double_t j2=1.063*1.063;
  
  fLongOnly = 1;
  Double_t pars[5];
  Double_t xs[5];
  xs[0] = 0.0;
  xs[1] = 0.0;
  xs[2] = qlong;
  pars[0] = aLambda;
  pars[1] = 1.0;
  pars[2] = 1.063;
  pars[3] = aL;
  pars[4] = aPhi;
  fLongOnly = 0;
  return XXX(&(xs[0]),&(pars[0]));
  
  Double_t result=aLambda*TMath::Exp(-qlong*qlong*j2/0.0389273);
  
  
  
//  ::Info("QLongCylSurf","l*exp %f",result);
  Double_t longpart;
  if ( (qlong != 0.0) && (aL != 0.0) )
   {
     Double_t qlL = qlong*aL/(2.0*0.197);
     Double_t sinqlL = TMath::Sin(qlL);
     Double_t sin2qlL = sinqlL*sinqlL;
     longpart = sin2qlL/(qlL*qlL);
     result*= longpart;
   }

//  ::Info("QLongCylSurf","longpart %f",result);
  
  if (aPhi > 0.0)
   {
     Double_t oneOverPhi = 1./aPhi;
     
     Double_t tmp1 = TMath::StruveL0(oneOverPhi)/TMath::BesselI0(oneOverPhi);
     Double_t tmp = 1.0 - tmp1*tmp1;
     result*=tmp;
     ::Info("QLongCylSurf","oneOverPhi %f, tmp1 %f, tmp %f",oneOverPhi ,tmp1, tmp);
   }
//  ::Info("QLongCylSurf","Phipart %f",result);

  result+=1.; 

//  ::Info("QLongCylSurf","Final: result=%f",result);
  return result;
  
}
Double_t AliHBTFits::XXX(Double_t *x, Double_t *par)
{

  //QOutQSideQLongCylSurf Function
  Double_t qout = x[0];
  Double_t qside = x[1];
  Double_t qlong = x[2];
  
  Double_t aLambda = par[0];
  Double_t aR = par[1];
  Double_t aJ = par[2];
  Double_t aL = par[3];
  Double_t aPhi = par[4];
  
  Double_t result=aLambda*TMath::Exp(-(qout*qout+qside*qside+qlong*qlong)*aJ*aJ/0.0389273);
  
//  ::Info("XXX","lambda*exp=%f",result);
  Double_t longpart;
  if ( (qlong != 0.0) && (aL != 0.0) )
   {
     Double_t qlL = qlong*aL/(2.0*0.197);
     Double_t sinqlL = TMath::Sin(qlL);
     Double_t sin2qlL = sinqlL*sinqlL;
     longpart = sin2qlL/(qlL*qlL);
     result*= longpart;
   }
   
//  ::Info("XXX","longpart=%f",longpart);
  
  if (aPhi > 0.0)
   {
    fgF1->SetParameter(0,aR);
    fgF2->SetParameter(0,aR);
  
    fgF1->SetParameter(1,aPhi);
    fgF2->SetParameter(1,aPhi);
  
    fgF1->SetParameter(2,qout);
    fgF2->SetParameter(2,qout);

    fgF1->SetParameter(3,qside);
    fgF2->SetParameter(3,qside);
    
    Double_t sp = fgF1->Integral(-TMath::Pi(),TMath::Pi(),(Double_t*)0x0,1.e-8);
    Double_t cp = fgF2->Integral(-TMath::Pi(),TMath::Pi(),(Double_t*)0x0,1.e-8);
    
    ::Info("QSideCylSurf","sp=%f, cp=%f",cp,sp);
    
    result *= cp*cp + sp*sp;
//    ::Info("QSideCylSurf","2: result=%f",result);

    
    Double_t oneOverPhi = 1./aPhi;
    Double_t tmp = TMath::BesselI0(oneOverPhi)-TMath::StruveL0(oneOverPhi);
    tmp = tmp*TMath::TwoPi();
    tmp = tmp*tmp;
    if (tmp == 0.0)
     {
       ::Error("QOutCylSurf","Div by 0");
       return 1.0;
     }

//    ::Info("XXX","itgrl/bessel=%f",(cp*cp + sp*sp)/tmp);
 
    result=result/tmp;
   }
   
  result+=1.; 

//  ::Info("QSideCylSurf","Final: result=%f",result);
  return result;

}
