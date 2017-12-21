#include "TH1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TStyle.h"

void FitPi0Eff(Int_t cen=6){

  TFile * f1 = new TFile("Eff_MB.root") ;
    
  char key[255] ;
  const Int_t nPID=4 ;
  char cPID[13][25] ;
  sprintf(cPID[0],"Emin3_All") ;
  sprintf(cPID[1],"Emin3_Disp");
  sprintf(cPID[2],"Emin3_CPV") ;
  sprintf(cPID[3],"Emin3_Both"); 
  TH1D * h1[nPID] ;
  TH1D * h1Int[nPID] ;
  for(Int_t iPID=0;iPID<nPID; iPID++){
    printf("iPID=%d \n",iPID) ;
    sprintf(key,"yeild1_GS_%s_cen%d",cPID[iPID],cen) ;
    h1[iPID]=(TH1D*)f1->Get(key) ;
    
    sprintf(key,"yeild1_int_GS_%s_cen%d",cPID[iPID],cen) ;
    printf("key=>%s< \n",key) ;
    h1Int[iPID]=(TH1D*)f1->Get(key) ;
  }
    
  
  Int_t col[13]={kRed,kGray+3,kBlue,kMagenta,kGreen+4,kCyan,kOrange+3,kViolet,kOrange+9,6,5,4,4} ;
  Int_t sym[13]={20,21,22,29,24,25,29,30,31,32,28,20,20} ;

  TF1 * efit = new TF1("efit","exp(-(1.+[1]*x+[2]*x*x+[5]*x*x*x+[7]*x*x*x*x)/([3]*x+[4]*x*x+[6]*x*x*x+[8]*x*x*x*x))",0.8,55.) ;
  
  //  TF1 * fit1 = new TF1("eff_Pi0_Gaus_7TeV","exp(-([0]+[1]*x+[2]*x*x+[3]*x*x*x)/(1.+[4]*x+[5]*x*x+[6]*x*x*x))",0.,30.) ;
  TF1 * fit1[nPID] ;
  TF1 * fit1Int[nPID] ;
  for(Int_t iPID=0;iPID<nPID; iPID++){
    h1[iPID]->SetMarkerStyle(24) ;
    h1[iPID]->SetMarkerColor(col[iPID]) ;
    h1[iPID]->SetLineColor(col[iPID]) ;


    sprintf(key,"eff_Pi0_Gaus_PbPb_%s_cen%d",cPID[iPID],cen) ;
     fit1[iPID] = new TF1(key,"[0]*(1.+[1]*x+[2]*x*x+[5]*x*x*x+[7]*x*x*x*x)/(1+[3]*x+[4]*x*x+[6]*x*x*x+[8]*x*x*x*x)",0.8,55.) ;
//      if(iPID<2)
//        fit1[iPID]->SetParameters(1.42178e-05,2.38111e+01,6.12189e+00,-4.22625e-01,9.19189e-02,1.53728e-01,-5.73258e-03,2.94763e-03,1.46434e-04) ; 
//      else
       fit1[iPID]->SetParameters(4.38276e-06,9.96840e+04,-1.85748e+05,9.56677e+01,-6.56810e+01,4.53750e+04,7.36595e+01,4.54515e+04,2.01751e+01 ) ;
//      fit1[iPID] = new TF1(key,"[0]*TMath::Exp(-(1.+[1]*x+[2]*x*x+[5]*x*x*x)/([3]*x+[4]*x*x+[6]*x*x*x))*(1.-TMath::TanH((x-[7])/[8]))",0.8,55.) ;
//      fit1[iPID]->SetParameters(3.85969e-03,-1.52994e+09,1.34488e+08,7.22354e+08,-1.03081e+09,1.28638e+06,2.29669e+07,3.41922e+01,7.58443e+00) ; 

     
     
    sprintf(key,"eff_int_Pi0_Gaus_PbPb_%s_cen%d",cPID[iPID],cen) ;
     fit1Int[iPID] = new TF1(key,"[0]*(1.+[1]*x+[2]*x*x+[5]*x*x*x+[7]*x*x*x*x)/(1+[3]*x+[4]*x*x+[6]*x*x*x+[8]*x*x*x*x)",0.8,55.) ;
     if(iPID<2)
       fit1Int[iPID]->SetParameters(1.42178e-05,2.38111e+01,6.12189e+00,-4.22625e-01,9.19189e-02,1.53728e-01,-5.73258e-03,2.94763e-03,1.46434e-04) ; 
     else
       fit1Int[iPID]->SetParameters(2.62836e-02,-6.11924e+03,8.07011e+03,3.63094e+04,5.17001e+02,-2.35977e+03,-3.10945e+03,3.23643e+02,8.36184e+02) ;
//     fit1Int[iPID] = new TF1(key,"exp(-(1.+[1]*x+[2]*x*x+[5]*x*x*x+[7]*x*x*x*x)/([3]*x+[4]*x*x+[6]*x*x*x+[8]*x*x*x*x))",0.8,55.) ;
//      fit1Int[iPID]->SetParameters(-9.47206e+05,-2.33561e+04,9.82162e+03,-2.31569e+03,-3.58733e+02,1.06534e+04,1.96734e+03,2.79445e+03,6.48405e+02) ;


    fit1[iPID]->SetLineColor(col[iPID]) ;
    fit1[iPID]->SetLineWidth(2) ;
    fit1Int[iPID]->SetLineColor(col[iPID]) ;
    fit1Int[iPID]->SetLineWidth(2) ;
    fit1Int[iPID]->SetLineStyle(4) ;

    Double_t rangeMin=1.1 ;
    Double_t rangeMax=40. ;
   h1[iPID]->Fit(fit1[iPID],"q","",2.,30.) ;
   h1[iPID]->Fit(fit1[iPID],"M","",rangeMin,rangeMax) ;
 /*  
   TFitResultPtr resultEff =  h1[iPID]->Fit(fit1[iPID],"MSN","",rangeMin,rangeMax) ;
    for(Int_t i=3;i<=h1[iPID]->GetNbinsX();i++){ //start from 3: for some centralities there are problems with eff<1 GeV
       Double_t xmin=h1[iPID]->GetXaxis()->GetBinLowEdge(i) ;
       Double_t xmax=h1[iPID]->GetXaxis()->GetBinUpEdge(i) ;
       Double_t binwidth=h1[iPID]->GetXaxis()->GetBinWidth(i) ;
       h1[iPID]->SetBinContent(i,efit->Integral(xmin,xmax,resultEff->GetParams())/binwidth) ;
       h1[iPID]->SetBinError(i,efit->IntegralError(xmin,xmax,resultEff->GetParams(),resultEff->GetCovarianceMatrix().GetMatrixArray())/binwidth) ;       
       // printf("i=%d, [%f,%f], E=%e, err=%e \n",i,xmin,xmax,efit->Integral(xmin,xmax,resultEffi->GetParams())/binwidth,efit->IntegralError(xmin,xmax,resultEff->GetParams(),resultEff->GetCovarianceMatrix().GetMatrixArray())/binwidth) ;      
    }
//    h1[iPID]->Draw() ; return ;
*/


    h1Int[iPID]->Fit(fit1Int[iPID],"qN","",1.5,25.) ;
    h1Int[iPID]->Fit(fit1Int[iPID],"qN","",rangeMin,rangeMax) ;
/*    
    TFitResultPtr resultEffi = h1Int[iPID]->Fit(fit1Int[iPID],"MSN","",rangeMin,rangeMax) ;
    for(Int_t i=3;i<=h1Int[iPID]->GetNbinsX();i++){
       Double_t xmin=h1Int[iPID]->GetXaxis()->GetBinLowEdge(i) ;
       Double_t xmax=h1Int[iPID]->GetXaxis()->GetBinUpEdge(i) ;
       Double_t binwidth=h1Int[iPID]->GetXaxis()->GetBinWidth(i) ;
       h1Int[iPID]->SetBinContent(i,efit->Integral(xmin,xmax,resultEffi->GetParams())/binwidth) ;
       h1Int[iPID]->SetBinError(i,efit->IntegralError(xmin,xmax,resultEffi->GetParams(),resultEffi->GetCovarianceMatrix().GetMatrixArray())/binwidth) ;
// printf("i=%d, [%f,%f], E=%e, err=%e \n",i,xmin,xmax,efit->Integral(xmin,xmax,resultEffi->GetParams())/binwidth,efit->IntegralError(xmin,xmax,resultEffi->GetParams(),resultEffi->GetCovarianceMatrix().GetMatrixArray())/binwidth) ;
      
    }
    */
//      h1[iPID]->Divide(fit1[iPID]) ;
  }  
  
//    h1[3]->Draw() ; return ;

  sprintf(key,"Efficiency_PID_cen%d",cen) ; 
  TCanvas * cEff = new TCanvas(key,key,10,10,900,900) ;
  cEff->Divide(2,2) ;
  gStyle->SetTitleSize(0.15,"t") ;
  TH1D * hBox = new TH1D("hBox","",250,0.6,40.) ;
   hBox->SetMinimum(1.e-6) ;
   hBox->SetMaximum(0.02) ;
  hBox->SetXTitle("p_{t} (GeV/c)") ;
  hBox->SetYTitle("Eff") ;
  hBox->SetStats(0) ;
  
  for(Int_t iPID=0;iPID<nPID; iPID++){
    cEff->cd(iPID+1) ;
    gPad->SetLogx(1) ;
    gPad->SetLogy(1) ;
    gPad->SetGridx() ;
    gPad->SetGridy() ;
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetFrameBorderMode(0);
    gPad->SetLeftMargin(0.06);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.05);
    gPad->SetRightMargin(0.01);
    hBox->SetTitle(cPID[iPID]) ;
    hBox->DrawClone() ;
//     h1[iPID]->SetMarkerColor(2) ;
    h1[iPID]->SetMarkerStyle(20) ;
    h1[iPID]->Draw("same") ;
    fit1[iPID]->DrawClone("same") ;
    fit1Int[iPID]->DrawClone("same") ;
  }

/*  
  //Evaluate errors
  
  cEff->cd(2) ;  
  TLegend * l = new TLegend(0.6,0.1,0.9,0.5) ;
  l->SetFillColor(0) ;
  l->AddEntry(h1[1],"0.9<p_{t}^{sim}<2.2 GeV","p") ;
  l->AddEntry(fit1[1],"Global fit","l") ;
  l->AddEntry(fit1Int[1],"Global fit, integral","l") ;

  l->Draw() ;

  TCanvas * ratio = new TCanvas("Ratios","",10,10,600,900) ;
  ratio->Divide(2,2) ; 
  hBox->SetMaximum(1.4)  ;
  hBox->SetMinimum(0.6)  ;
  for(Int_t iPID=0;iPID<nPID; iPID++){
    ratio->cd(iPID+1) ;
    gPad->SetLogx(1) ;
    gPad->SetGridx() ;
    gPad->SetGridy() ;
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetFrameBorderMode(0);
    gPad->SetLeftMargin(0.06);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.05);
    gPad->SetRightMargin(0.01);
     hBox->SetTitle(cPID[iPID]) ;
    hBox->DrawClone() ;
    h1[iPID]->Divide(fit1[iPID]) ;
    h1[iPID]->GetXaxis()->SetRangeUser(0.5,5.) ;
    h1[iPID]->DrawClone("same") ;
  }
  
  
*/  

  TFile fout("PHOS_pi0_eff.root","update") ;
  for(Int_t pid=0;pid<nPID; pid++){
    fit1[pid]->Write(0,TObject::kOverwrite) ;
    fit1Int[pid]->Write(0,TObject::kOverwrite) ;
// printf("========pid=%d\n",pid) ;    
//     for(Int_t i=1;i<=h1[pid]->GetNbinsX();i++){ //start from 3: for some centralities there are problems with eff<1 GeV
//       printf("ha(%d)=%f, %f,  %f, %f \n",i,h1[pid]->GetBinContent(i),h1[pid]->GetBinError(i),h1Int[pid]->GetBinContent(i),h1Int[pid]->GetBinError(i)) ;
//     }
    h1[pid]->Write(0,TObject::kOverwrite) ;
    h1Int[pid]->Write(0,TObject::kOverwrite) ;
  }


}
