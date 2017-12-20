/// \ingroup STAT/test
/// \brief  test of AliTMinuiToolkit class - linear fits
/// * For the moment: Test only that code does not crash
/// * Example to run (300 times)
/// \code   .x $AliRoot_SRC/STAT/test/AliTMinuitToolkitTestLinear.C+(300,3) \endcode
/// * Test of linear fit for data with outliers
///   * Compare fits:
///      * 0.) Standard linear pol1 fit
///      * 1.) "robust"  pol1 fit
///      * 2.) Standard pol1 minuit fit
///      * 3.) AliTMinuitToolkit with gaus+cauchy likelihood
///      * 4.) "robust" AliTMinuitToolkit with Huber cost function
///      * 5.) Bootstrap with 20 iteration
///      * 6.) Bootstrap with 100 iteration
///      * 7.) MISAC     with 200 iteration
///   * Input data: pol1 + gaussian noise with tails (random fraction: 20 %, 20 % wider)
///      \code pol1+ ((rmdm>x) ? rndmG :20*rndmG) \endcode
///   * Test example result:
///     * example canvas showing fit quality for few random samples
///     * Canvas with residual histograms for each fit option
///       * misac option - smallest fraction of outlier
///       * robust option - smallest sigma (but some outlier)
///       * standard chi2 minimization (~factor 5-10 times wider than other "robust" method
/// * Example usage in the $AliRoot_SRC/STAT/test/
///   * TODO - save example images in the doxygen/html directory (STEERING part ?) - currently saved by hand - to check with OFFLINE
///   * $AliRoot_SRC/STAT/test/AliTMinuitToolkitTest.C+
///   * $AliRoot_SRC/STAT/test/AliTMinuitToolkitTestLinear.C+
///     *  compare fits example with different likelihood/resp fits strategies for data with outliers
///        \image html AliTMinuiToolkit.TestLinearFitExample.png
///     *  compare performance fits example with different likelihood/resp fits strategies for data with outliers
///        \image html AliTMinuiToolkit.TestLinearFitStatExample.png


#include "AliTMinuitToolkit.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TTreeStream.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TRandom.h"
#include "TStatToolkit.h"
#include "TLegend.h"
#include "TStyle.h"

void AliTMinuitToolkitTestLinear(Int_t nIter, Int_t nDraw){
  //
  // nIter=1000; nDraw=4;
  AliTMinuitToolkit::RegisterDefaultFitters();
  AliTMinuitToolkit::GetPredefinedFitter("pol1R")->SetStreamer("AliTMinuitToolkit_Pol1RStreamer.root");
  ((TFormula*)AliTMinuitToolkit::GetPredefinedFitter("pol1R")->GetLogLikelihoodFunction())->SetParameter(0,0.7);
  //
  TMatrixD &pInit = *((AliTMinuitToolkit::GetPredefinedFitter("pol1R")->GetInitialParam()));
  TMatrixD &pInitH = *((AliTMinuitToolkit::GetPredefinedFitter("pol1H")->GetInitialParam()));
  pInit(0,0)=0; pInit(0,1)=10;
  pInit(1,0)=0; pInit(1,1)=10;
  TTreeSRedirector *pcstream=new TTreeSRedirector("AliTMinuiToolkit.TestLinearFit.root","recreate");
  TCanvas *canvasTestLinearFit=new TCanvas("TestLinearGraphFit","TestLinearGraphFit",1200,1000);
  canvasTestLinearFit->Divide(nDraw,nDraw);
  //
  // 1.) generate random graphs and make fits
  TVectorD fit0(2),fit1(2),fit2(2),fit3(2), fit4(2), fit5(2), fit6(2),fit7(2), fit8(2);
  TVectorD rms5(2), rms6(2),rms7(2),rms8(2);
  TMatrixD mFit2(2,2),mFit3(2,2), mFit4(2,2),mFit5(2,2);
  TF1 *f1= new TF1("fPol1","pol1",-1,1);
  TVectorD vecX(10000);
  TVectorD vecY(10000);
  TVectorD vecE(10000);
  for (Int_t it=0; it<nIter; it++){
    if (it%100==0) printf("%d\n",it);
    Int_t jt=it%(nDraw*nDraw);
    canvasTestLinearFit->cd(jt+1);
    Double_t p0=(gRandom->Rndm()-0.5)*20;
    Double_t p1=(gRandom->Rndm()-0.5)*20;
    Double_t frac=0.2+gRandom->Rndm()*0.2;
    Int_t nPoints=TMath::Nint(10+gRandom->Rndm()*10);
    Double_t width=20;
    for (Int_t iPoint=0; iPoint<nPoints; iPoint++){
      vecX[iPoint]=2*(gRandom->Rndm()-0.5);
      Double_t noise = gRandom->Gaus()*((iPoint<frac*nPoints)?width:1);
      vecY[iPoint]=p0+vecX[iPoint]*p1+noise;
      vecE[iPoint]=1;
    }
    TGraphErrors *gr = new TGraphErrors(nPoints,vecX.GetMatrixArray(),vecY.GetMatrixArray(),0,vecE.GetMatrixArray());
    gr->SetMarkerStyle(21);
    gr->Draw("ap");
    f1->SetLineColor(4);
    gr->Fit(f1,"+","");
    f1->GetParameters(fit0.GetMatrixArray());
    f1->SetLineColor(3);
    gr->Fit(f1,"+N+ROB=0.7","");
    f1->GetParameters(fit1.GetMatrixArray());
    //
    AliTMinuitToolkit*fitter2 = AliTMinuitToolkit::Fit(gr,"pol1","default", "", "funOption(1,3,1)");  // standard minuit
    fit2=*(fitter2->GetParameters());     mFit2=*(fitter2->GetCovarianceMatrix());
    AliTMinuitToolkit*fitter3 = AliTMinuitToolkit::Fit(gr,"pol1R","iter1", "", "funOption(2,3,1)"); // logLike:gaus+cauchy
    fit3=*(fitter3->GetParameters());     mFit3=*(fitter3->GetCovarianceMatrix());
    AliTMinuitToolkit*fitter4 = AliTMinuitToolkit::Fit(gr,"pol1H","", "", "funOption(2,3,2)"); // logLike: huber
    fit4=*(fitter4->GetParameters());     mFit4=*(fitter3->GetCovarianceMatrix());
    AliTMinuitToolkit*fitter5 = AliTMinuitToolkit::Fit(gr,"pol1R","bootstrap20", "","funOption(6,3,4)");  // bootstrap20
    fit5=*(fitter5->GetParameters());     rms5=*(fitter5->GetRMSEstimator());
    AliTMinuitToolkit*fitter6 = AliTMinuitToolkit::Fit(gr,"pol1R","bootstrap50", "", "funOption(3,3,5)" ); // bootstrap50
    fit6=*(fitter6->GetParameters());     rms6=*(fitter6->GetRMSEstimator());
    AliTMinuitToolkit*fitter7 = AliTMinuitToolkit::Fit(gr,"pol1R","misac(4,20)","","funOption(3,3,6)");    // misac 20
    fit7=*(fitter7->GetParameters()); rms7=*(fitter7->GetRMSEstimator());
    AliTMinuitToolkit*fitter8 = AliTMinuitToolkit::Fit(gr,"pol1R","misac(4,50)","","funOption(3,3,6)");    // misac 50
    TMatrixD *misacEstimator=(TMatrixD*)(fitter8->GetMISACEstimators());
    fit8=*(fitter8->GetParameters()); rms8=*(fitter8->GetRMSEstimator());
    gr->Draw("ap");
    TList* farray=gr->GetListOfFunctions();
    TLegend *legend = new TLegend(0.11,0.1,0.89,0.2,"AliTMinuitToolkit::pol1 fit");
    legend->SetBorderSize(0);
    legend->SetNColumns(3);
    for (Int_t iFun=0; iFun<farray->GetEntries(); iFun++){
      legend->AddEntry(farray->At(iFun),farray->At(iFun)->GetName(),"l");
    }
    legend->Draw();
    gr->Draw("p");

    (*pcstream)<<"polFit"<<
      "p0="<<p0<<
      "p1="<<p1<<
      "frac="<<frac<<
      "nPoints="<<nPoints<<
      "gr.="<<gr<<
      "fit0.="<<&fit0<<   // normal
      "fit1.="<<&fit1<<   // robust linear
      "fit2.="<<&fit2<<   // minuit default
      "fit3.="<<&fit3<<   // logLike:gaus+cauchy
      "fit4.="<<&fit4<<   // logLike: huber
      "fit5.="<<&fit5<<   // bootstrap 20
      "fit6.="<<&fit6<<   // bootstrap 50
      "fit7.="<<&fit7<<   // MISAC 20
      "fit8.="<<&fit8<<   // MISAC 50
      "rms5.="<<&rms5<<
      "rms6.="<<&rms6<<
      "rms7.="<<&rms7<<
      "mFit2.="<<&mFit2<<
      "mFit3.="<<&mFit3<<
      "mFit4.="<<&mFit4<<
      "misacE.="<<misacEstimator<<
      "\n";
  }
  canvasTestLinearFit->SaveAs("AliTMinuiToolkit.TestLinearFitExample.png");
  canvasTestLinearFit->SaveAs("AliTMinuiToolkit.TestLinearFitExample.pdf");
  delete pcstream;
  //
  // 2. Draw some results
  //
  TFile * f = TFile::Open("AliTMinuiToolkit.TestLinearFit.root");
  TTree * treeOut=(TTree*)f->Get("polFit");
  TCanvas *canvasStatComparison = new TCanvas("StatComparison","StatComparison",1200,1000);
  canvasStatComparison->Divide(3,3);
  //
  gStyle->SetOptFit();
  const char  *fNames[9]={"ROOT:Default","ROOT:ROB0.70","AliTMinuitToolkit:Default", "LogLike:Gaus(80%)Cauchy(20%)- 1 Iter",\
			  "LogLike:Huber. 1 Iter", "LogLike:Gaus(80%)Cauchy(20%)- Bootstrap(20)", "LogLike:Gaus(80%)Cauchy(20%)- Bootstrap(50)", \
			  "LogLike:Gaus(80%)Cauchy(20%)- MISAC(20)", "LogLike:Gaus(80%)Cauchy(20%)- MISAC(50)"};
  for (Int_t iFit=0; iFit<9; iFit++){
    canvasStatComparison->cd(iFit+1)->SetLogy();
    treeOut->Draw(TString::Format("fit%d.fElements[1]-p1>>hisDelta0_%d(100,-20,20)",iFit,iFit));
    treeOut->GetHistogram()->SetTitle(fNames[iFit]);
    treeOut->GetHistogram()->GetXaxis()->SetTitle("p_{1fit}-p_{1in}");
    treeOut->GetHistogram()->Fit("gaus");
    ::Info("AliTMinuitToolkitTestLinear", Form("Fit%s",fNames[iFit]));
  }
  canvasStatComparison->SaveAs("AliTMinuiToolkit.TestLinearFitStatExample.png");
  canvasStatComparison->SaveAs("AliTMinuiToolkit.TestLinearFitStatExample.pdf");
  //
  // 3. Make more detailed anaysis of the trees
  /*
    TTree * treeBootstrap = ((*(AliTMinuitToolkit::GetPredefinedFitter("pol1R")->GetStreamer()))<<"bootstrap").GetTree();
    TTree * treeTwoFold = ((*(AliTMinuitToolkit::GetPredefinedFitter("pol1R")->GetStreamer()))<<"crossValidation").GetTree();
    TTree * treeMISAC = ((*(AliTMinuitToolkit::GetPredefinedFitter("pol1R")->GetStreamer()))<<"misac").GetTree();
  */

}
