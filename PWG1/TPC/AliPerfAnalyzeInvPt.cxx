#include "TNamed.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TF1.h"

//#include "AliPerformancePtCalibMC.h"
//#include "AliPerformancePtCalib.h"
#include "AliPerfAnalyzeInvPt.h"

ClassImp(AliPerfAnalyzeInvPt)


// fit functions
//____________________________________________________________________________________________________________________________________________
   Double_t AliPerfAnalyzeInvPt::Polynomial(Double_t *x, Double_t *par)
{


   if (x[0] > -par[4] && x[0] < par[4]) {
      TF1::RejectPoint();
      return 0;
   }
   return  par[2]+par[0]*pow((x[0]-par[1]),2)+par[3]*pow((x[0]-par[1]),4);
  
}
//____________________________________________________________________________________________________________________________________________
Double_t AliPerfAnalyzeInvPt::PolynomialRejP(Double_t *x, Double_t *par)
{
 
   Double_t pos  = par[5];
   Double_t neg = - par[5];
   pos += par[4];
   neg += par[4];
  
   if (x[0] > neg && x[0] < pos) {
      TF1::RejectPoint();
      return 0;
   }
   return   par[2]+par[0]*pow((x[0]-par[1]),2)+par[3]*pow((x[0]-par[1]),4);
}
//____________________________________________________________________________________________________________________________________________
Double_t AliPerfAnalyzeInvPt::InvGauss(Double_t *x, Double_t *par)
{
   if (x[0] > -par[6] && x[0] < par[6]) {
      TF1::RejectPoint();
      return 0;
   }

   return par[3]+par[0]*TMath::Exp(-0.5*(TMath::Power((x[0]-par[1])/par[2], 2.0)))+par[4]*pow((x[0]-par[1]),2)+par[5]*pow((x[0]-par[1]),4) ;
}
//____________________________________________________________________________________________________________________________________________
Double_t AliPerfAnalyzeInvPt::InvGaussRejP(Double_t *x, Double_t *par)
{
   Double_t pos  = par[7];//0.12;
   Double_t neg = - par[7];//0.12;
   pos += par[6];
   neg += par[6];
  
   if (x[0] > neg && x[0] < pos) {
      TF1::RejectPoint();
      return 0;
   }


   return par[3]+par[0]*TMath::Exp(-0.5*(TMath::Power((x[0]-par[1])/par[2], 2.0)))+par[4]*pow((x[0]-par[1]),2)+par[5]*pow((x[0]-par[1]),4) ;
}


//_____________________________________________________________________________________________________________________________________________
AliPerfAnalyzeInvPt::AliPerfAnalyzeInvPt():
   TNamed("AliPerfAnalyzeInvPt","AliPerfAnalyzeInvPt"),
   fNThetaBins(0), 
   fNPhiBins(0),
   fRange(0),
   fExclRange(0),
   fFitGaus(0) ,
   fHistH2InvPtTheta(0),
   fHistH2InvPtPhi(0), 
   fGrMinPosTheta(0),
   fGrMinPosPhi(0),
   fFitMinPos(0),
   fFitMinPosRejP(0),
   fFitInvGauss(0),
   fFitInvGaussRejP(0)
{
   // Default constructor
  
   fFitGaus = kFALSE;
   fNThetaBins = 0;
   fNPhiBins = 0;
   fRange = 0;
   fExclRange = 0;
   fFitGaus = 0;
   
   // projection histos
   TH1D *fHistFitTheta[100];
   TH1D *fHistFitPhi[100];

   for(Int_t i=0;i<100;i++){
      
      fHistFitTheta[i] = NULL;
      fHistFitPhi[i] = NULL;
   }
  
}
//_____________________________________________________________________________________________________________________________________________
AliPerfAnalyzeInvPt::AliPerfAnalyzeInvPt(Char_t* name="AliAnalyzeInvPt",Char_t* title="AliAnalyzeInvPt"): TNamed(name, title),
   fNThetaBins(0), 
   fNPhiBins(0),
   fRange(0),
   fExclRange(0),
   fFitGaus(0) ,
   fHistH2InvPtTheta(0),
   fHistH2InvPtPhi(0), 
   fGrMinPosTheta(0),
   fGrMinPosPhi(0),
   fFitMinPos(0),
   fFitMinPosRejP(0),
   fFitInvGauss(0),
   fFitInvGaussRejP(0)
   {
      //  Double_t fThetaBins[100] = {0};
      //   Double_t fPhiBins[100] = {0};
  
      fFitGaus = kFALSE;
      fNThetaBins = 0;
      fNPhiBins =0;
      fRange = 0;
      fExclRange = 0;
      fFitGaus = 0;
      // projection histos
      TH1D *fHistFitTheta[100];
      TH1D *fHistFitPhi[100];

      for(Int_t i=0;i<100;i++){
    
	 fHistFitTheta[i] = NULL;
	 fHistFitPhi[i] = NULL;
      }
      fHistH2InvPtTheta = NULL;
      fHistH2InvPtPhi = NULL; 
      fGrMinPosTheta= NULL;
      fGrMinPosPhi= NULL;
   }


   //______________________________________________________________________________________________________________________________________
void AliPerfAnalyzeInvPt::InitHistos(Double_t *binsXTheta,Double_t *fitParamTheta,Double_t *errFitParamTheta,Double_t *binsXPhi,Double_t *fitParamPhi,Double_t *errFitParamPhi){// init Histos


      
  
      fGrMinPosTheta = new TGraphErrors(fNThetaBins,binsXTheta,fitParamTheta,0, errFitParamTheta);  
      fGrMinPosTheta->SetMarkerStyle(20);
      fGrMinPosTheta->SetMarkerColor(2);
      fGrMinPosTheta->SetLineColor(2);
      fGrMinPosTheta->GetYaxis()->SetTitle("min pos (Gev/c)^{-1}");
      fGrMinPosTheta->GetXaxis()->SetTitle("#theta bin no.");
      fGrMinPosTheta->GetYaxis()->SetTitleOffset(1.2);   
      fGrMinPosTheta->SetTitle("#theta bins ");

      fGrMinPosPhi = new TGraphErrors(fNPhiBins,binsXPhi,fitParamPhi,0,errFitParamPhi);  
      fGrMinPosPhi->SetMarkerStyle(20);
      fGrMinPosPhi->SetMarkerColor(4);
      fGrMinPosPhi->SetLineColor(4);
      fGrMinPosPhi->GetYaxis()->SetTitle("min pos (Gev/c)^{-1}");
      fGrMinPosPhi->GetXaxis()->SetTitle("#phi bin no.");
      fGrMinPosPhi->GetYaxis()->SetTitleOffset(1.2);   
      fGrMinPosPhi->SetTitle("#phi bins ");
}

//______________________________________________________________________________________________________________________________________

void AliPerfAnalyzeInvPt::InitFitFcn(){
      // fit functions
      fFitMinPos = new TF1("fFitMinPos", Polynomial,-4.0,4.0,5);
      fFitMinPos->SetLineColor(4);
      fFitMinPos->SetLineWidth(1);
      fFitMinPos->SetParameter(0,1.0);
      fFitMinPos->SetParameter(1,0.0);
      fFitMinPos->SetParameter(2,0.0);

      fFitMinPosRejP = new TF1("fFitMinPosRejP",PolynomialRejP,-4.0,4.0,6);
      fFitMinPosRejP->SetLineColor(2);
      fFitMinPosRejP->SetLineWidth(1);
      fFitMinPosRejP->SetParameter(0,1.0);
      fFitMinPosRejP->SetParameter(1,0.0);
      fFitMinPosRejP->SetParameter(2,0.0);
 
  
      fFitInvGauss = new TF1("fFitInvGauss", InvGauss,-4.0,4.0,6);
      fFitInvGauss->SetLineColor(4);
      fFitInvGauss->SetLineWidth(1);
      fFitInvGauss->SetParameter(2,1.0);
  
      fFitInvGaussRejP = new TF1("fFitInvGaussRejP", InvGaussRejP,-4.0,4.0,7);
      fFitInvGaussRejP->SetLineColor(2);
      fFitInvGaussRejP->SetLineWidth(1);
      fFitInvGaussRejP->SetParameter(2,1.0);
 
   }
//______________________________________________________________________________________________________________________________________
void AliPerfAnalyzeInvPt::StartAnalysis(TH2F *histThetaInvPt, TH2F *histPhiInvPt, TObjArray* aFolderObj){//start Ana

  
   
   if(!histThetaInvPt) {
      Printf("warning: no 1/pt histogram to analyse in theta bins!");
   }
   if(!histPhiInvPt) {
      Printf("warning: no 1/pt histogram to analyse in phit bins!");
   }

   Double_t thetaBins[9] = {0.77,0.97,1.17,1.37,1.57,1.77,1.97,2.17,2.37};                 // theta bins
   Double_t phiBins[13] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.5};           // phi bins
  
   Int_t nThetaBins = 9;
   Int_t nPhiBins = 13;
   Double_t range = 1.0;        // fit range
   Double_t  exclRange  = 0.13; // range of point rejection in fit
   Bool_t fitGaus=kFALSE;       // fit with Gaussian or with polynomial (kFALSE)
   
   if(fNThetaBins == 0){
      fNThetaBins = nThetaBins;
      for(Int_t k = 0;k<fNThetaBins;k++){
	 fThetaBins[k] = thetaBins[k];
      }
      Printf("warning: for theta bins no user intput! bins are set to default values.");
   }

   if(fNPhiBins == 0){
      fNPhiBins = nPhiBins;
      for(Int_t k = 0;k<fNPhiBins;k++){
	 fPhiBins[k] = phiBins[k];
      }
    
      Printf("warning: for phi bins no user intput! bins are set to default values.");
    
   }

   if(fRange==0){
      fRange = range;
      fExclRange = exclRange;
      fFitGaus = fitGaus;
      Printf("warning: no fit range is set. default fitting conditions are used.");
   }

      
   //param arrays
   Double_t fitParamTheta[100],fitParamPhi[100];
   Double_t errFitParamTheta[100], errFitParamPhi[100];
   Double_t binsXTheta[100],binsXPhi[100];

   fHistH2InvPtTheta  = histThetaInvPt;
   fHistH2InvPtPhi    = histPhiInvPt;

   InitFitFcn();
	 
   TCanvas *theCan =new TCanvas("theCan","invPt theta bins",1200,900);
   theCan->Divide((fNThetaBins+2)/2,2);
   TCanvas *phiCan =new TCanvas("phiCan","invPt phi bins",1200,900);
   phiCan->Divide((fNPhiBins+2)/2,2);
   Int_t countPad = 1;
	 
   // analyse 1/pt in bins of theta 
   for(Int_t i=0;i<fNThetaBins;i++){
      TString name = "fit_theta_";
      name +=i;
      Int_t firstBin= fHistH2InvPtTheta->GetYaxis()->FindBin(fThetaBins[i]);
      if(i>0) firstBin +=1;
      Int_t lastBin = fHistH2InvPtTheta->GetYaxis()->FindBin(fThetaBins[i+1]);
      if( i == fNThetaBins-1) {
	 firstBin= fHistH2InvPtTheta->GetYaxis()->FindBin(fThetaBins[0]);
	 lastBin = fHistH2InvPtTheta->GetYaxis()->FindBin(fThetaBins[i]);
      }
    
      fHistH2InvPtTheta->SetName(name.Data());
      fHistFitTheta[i] =  (TH1F*)fHistH2InvPtTheta->ProjectionX("_px",firstBin,lastBin,"e");
      
      Char_t titleTheta[50];
      if(i == fNThetaBins-1) sprintf(titleTheta,"1/pt (GeV/c) integrated over #theta");
      else  sprintf(titleTheta,"1/pt (GeV/c) for #theta range: %1.3f - %1.3f",fThetaBins[i],fThetaBins[i+1]);
      
      fHistFitTheta[i]->SetTitle(titleTheta);
   
      Double_t invPtMinPos  = 0;
      Double_t invPtMinPosErr = 0;
      Double_t invPtMinPosImpr  = 0;
      Double_t invPtMinPosErrImpr = 0;
   

	 
      if(fFitGaus==kFALSE){
	 Printf("making polynomial fit in 1/pt in theta bins");
	 theCan->cd(countPad);
	 MakeFit(fHistFitTheta[i],fFitMinPos, invPtMinPos,invPtMinPosErr, fExclRange,fRange);
	 MakeFitBetter(fHistFitTheta[i],fFitMinPosRejP, invPtMinPosImpr, invPtMinPosErrImpr, invPtMinPos,fExclRange,fRange);

	
	 fHistFitTheta[i]->DrawCopy();
	 fFitMinPos->DrawCopy("L,same");
	 fFitMinPosRejP->DrawCopy("L,same");
      }
      else{
	 Printf("making gauss fit in 1/pt in theta bins");
	 theCan->cd(countPad);
	 MakeFitInvGauss(fHistFitTheta[i],fFitInvGauss, invPtMinPos, invPtMinPosErr,fExclRange,fRange);
	 MakeFitInvGaussBetter(fHistFitTheta[i],fFitInvGaussRejP, invPtMinPosImpr, invPtMinPosErrImpr,invPtMinPos,fExclRange,fRange);
      
	 fHistFitTheta[i]->DrawCopy();
	 fFitInvGauss->DrawCopy("L,same");
	 fFitInvGaussRejP->DrawCopy("L,same");
      }
    
      aFolderObj->Add(fHistFitTheta[i]);
    
      fitParamTheta[i] = invPtMinPosImpr;
      errFitParamTheta[i] = invPtMinPosErrImpr;
    
      binsXTheta[i] = i+1.0;
      countPad++;
   }
      
      
   countPad = 1;

   
   // analyse 1/pt in bins of phi 
  
   for(Int_t i=0;i<fNPhiBins;i++){
      TString name = "fit_phi_";
      name +=i;
    
      fHistH2InvPtPhi->SetName(name.Data());
      Int_t  firstBin = fHistH2InvPtPhi->GetYaxis()->FindBin(fPhiBins[i]);
      if(i>0) firstBin +=1;
      Int_t   lastBin =  fHistH2InvPtPhi->GetYaxis()->FindBin(fPhiBins[i+1]);
      if(i == fNPhiBins-1){
	 firstBin = fHistH2InvPtPhi->GetYaxis()->FindBin(fPhiBins[0]);
	 lastBin =  fHistH2InvPtPhi->GetYaxis()->FindBin(fPhiBins[i]);
      }
      fHistFitPhi[i] =  (TH1F*) fHistH2InvPtPhi->ProjectionX("_px",firstBin,lastBin,"e");
      
      Char_t titlePhi[50];
      if(i == fNPhiBins-1) sprintf(titlePhi,"1/pt (GeV/c) integrated over #phi");
	  else  sprintf(titlePhi,"1/pt (GeV/c) for #phi range: %1.3f - %1.3f",fPhiBins[i],fPhiBins[i+1]);
     
      fHistFitPhi[i]->SetTitle(titlePhi);
  
      Double_t invPtMinPos  = 0;
      Double_t invPtMinPosErr = 0;
      Double_t invPtMinPosImpr  = 0;
      Double_t invPtMinPosErrImpr = 0;
    
      if(fFitGaus==kFALSE){
	 Printf("making polynomial fit in 1/pt in phi bins");
	 phiCan->cd(countPad);
	 MakeFit(fHistFitPhi[i],fFitMinPos, invPtMinPos, invPtMinPosErr,fExclRange,fRange);
	 MakeFitBetter(fHistFitPhi[i],fFitMinPosRejP, invPtMinPosImpr, invPtMinPosErrImpr,invPtMinPos,fExclRange,fRange);
	 
	 fHistFitPhi[i]->DrawCopy();
	 fFitMinPos->DrawCopy("L,same");
	 fFitMinPosRejP->DrawCopy("L,same");

      }
      else {
	 Printf("making gauss fit in 1/pt in phi bins");
	 phiCan->cd(countPad);
	 MakeFitInvGauss(fHistFitPhi[i],fFitInvGauss, invPtMinPos, invPtMinPosErr, exclRange,fRange);
	 MakeFitInvGaussBetter(fHistFitPhi[i],fFitInvGaussRejP, invPtMinPosImpr, invPtMinPosErrImpr, invPtMinPos, fExclRange,fRange);
	 
	 fHistFitPhi[i]->DrawCopy();
	 fFitInvGauss->DrawCopy("L,same");
	 fFitInvGaussRejP->DrawCopy("L,same");
      }
    
      aFolderObj->Add(fHistFitPhi[i]);
    
      fitParamPhi[i] = invPtMinPosImpr;
      errFitParamPhi[i] = invPtMinPosErrImpr;
    
      binsXPhi[i] = i+1.0;
      countPad++;
   }
      
   InitHistos(binsXTheta,fitParamTheta,errFitParamTheta,binsXPhi,fitParamPhi,errFitParamPhi);

   TCanvas *canFitVal = new TCanvas("canFitVal","min pos histos",800,400);
   canFitVal->Divide(2,1);

   canFitVal->cd(1);
   fGrMinPosTheta->Draw("ALP");
   canFitVal->cd(2);
   fGrMinPosPhi->Draw("ALP");

   Printf("AliPerfAnalyzeInvPt: NOTE: last bin is always fit result  of integral over all angle ranges which have been set by user!");
   
   aFolderObj->Add(fGrMinPosTheta);
   aFolderObj->Add(fGrMinPosPhi);
   aFolderObj->Add(fHistH2InvPtTheta);
   aFolderObj->Add(fHistH2InvPtPhi);
  
  
}


//____________________________________________________________________________________________________________________________________________

void AliPerfAnalyzeInvPt::MakeFit(TH1 *hproy, TF1 * fitpb, Double_t &mean, Double_t &errMean, Double_t &excl,Double_t &range)
{
 
   fitpb->SetRange(-range,range);
   fitpb->SetParLimits(1,-0.05,0.05);
   fitpb->SetParameter(0,1.0);
   fitpb->FixParameter(4,excl);
   fitpb->FixParameter(2,0.0);
    
   hproy->Fit(fitpb,"RQM");
   mean = fitpb->GetParameter(1);
   errMean = fitpb->GetParError(1);
   if(mean == 0)  errMean = 0;
}
//____________________________________________________________________________________________________________________________________________
void AliPerfAnalyzeInvPt::MakeFitBetter(TH1 *hproy, TF1 * fitpb2, Double_t &mean, Double_t &errMean, Double_t &f, Double_t &excl, Double_t &range)
{
   fitpb2->FixParameter(5,excl);
   fitpb2->FixParameter(4,f);
   // fitpb2->FixParameter(2,0.0);
   fitpb2->SetRange(-range+f,range+f);// was 0.25 0.45
   fitpb2->SetParLimits(1,-0.05,0.05);
   fitpb2->SetParameter(0,1.0);
   hproy->Fit(fitpb2,"RQM");
   mean = fitpb2->GetParameter(1);
   errMean = fitpb2->GetParError(1);
   if(mean == 0)   errMean = 0;

}
//____________________________________________________________________________________________________________________________________________
void AliPerfAnalyzeInvPt::MakeFitInvGauss(TH1 *hproy, TF1 * fitpb, Double_t &mean, Double_t &errMean, Double_t &excl,Double_t &range)
{
   fitpb->FixParameter(6,excl);
   fitpb->SetRange(-range,range);
   fitpb->SetParameter(0,-1.0);
   fitpb->SetParameter(2,1.0);
   fitpb->SetParameter(3,25000.0);
   fitpb->SetParameter(4,-1.0);
   fitpb->SetParLimits(1,-0.02,0.02);
   hproy->Fit(fitpb,"RQM");
   mean = fitpb->GetParameter(1);
   errMean = fitpb->GetParError(1);
   if(mean == 0)   errMean = 0;
  
}
//____________________________________________________________________________________________________________________________________________
void AliPerfAnalyzeInvPt::MakeFitInvGaussBetter(TH1 *hproy, TF1 * fitpb2, Double_t &mean, Double_t &errMean, Double_t &f, Double_t &excl, Double_t &range)
{
   fitpb2->FixParameter(7,excl);
   fitpb2->FixParameter(6,f);
   fitpb2->SetRange(-range+f,range+f);
   fitpb2->SetParameter(0,-1.0);
   fitpb2->SetParameter(2,1.0);
   fitpb2->SetParameter(3,25000.0);
   fitpb2->SetParameter(4,-1.0);
   fitpb2->SetParLimits(1,-0.02,0.02);
   hproy->Fit(fitpb2,"RQM");
   mean = fitpb2->GetParameter(1);
   errMean = fitpb2->GetParError(1);
   if(mean == 0)   errMean = 0;
   
}


//____________________________________________________________________________________________________________________________________________
// set variables 


void AliPerfAnalyzeInvPt::SetProjBinsPhi(const Double_t *phiBinArray,Int_t nphBins){
   
   fNPhiBins = nphBins;
  
   for(Int_t k = 0;k<fNPhiBins;k++){
      fPhiBins[k] = phiBinArray[k];
   }
   Printf("number of bins in phi set to %i",fNPhiBins);

}
//____________________________________________________________________________________________________________________________________________
void AliPerfAnalyzeInvPt::SetProjBinsTheta(const Double_t *thetaBinArray, Int_t nthBins){
  
   fNThetaBins = nthBins;
   for(Int_t k = 0;k<fNThetaBins;k++){
      fThetaBins[k] = thetaBinArray[k];
   }
   Printf("number of bins in theta set to %i",fNThetaBins);

}
//____________________________________________________________________________________________________________________________________________
void AliPerfAnalyzeInvPt::SetMakeFitOption(const Bool_t setGausFit, const Double_t exclusionR,const Double_t fitR ){

  
   
   fFitGaus = setGausFit;
   fExclRange  = exclusionR;
   fRange = fitR;
  
   if(fFitGaus) Printf("set MakeGausFit with fit range %2.3f and exclusion range in 1/pt: %2.3f",fRange,fExclRange);
   else  Printf("set standard polynomial fit with fit range %2.3f and exclusion range in 1/pt: %2.3f",fRange,fExclRange);
 
}
