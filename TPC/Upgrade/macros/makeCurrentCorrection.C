/*
.x $HOME/rootlogon.C
.L $ALICE_ROOT/TPC/Upgrade/macros/makeCurrentCorrection.C+

*/

#include "TMath.h"
#include "TRandom.h"
#include "TTreeStream.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "AliTPCParam.h"
#include "AliTPCcalibDB.h"
#include "AliTPCAltroMapping.h"
#include "AliAltroRawStream.h"
#include "AliSysInfo.h"
#include "AliTPCRawStreamV3.h"
#include "AliCDBManager.h"
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"
#include "AliRawReaderRoot.h"
#include "AliRawReader.h"
#include "TH3.h"
#include "TH2.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"
#include "TChain.h"
#include "AliXRDPROOFtoolkit.h"
#include "TLegend.h"
#include "TCut.h"
#include "TGraphErrors.h"
#include "TStatToolkit.h"
#include "TF2.h"

#include "AliDCSSensor.h"
#include "AliCDBEntry.h"
#include "AliDCSSensorArray.h"
#include "TStyle.h"
#include "AliTPCSpaceCharge3D.h"
#include "AliExternalTrackParam.h"
#include "AliTrackerBase.h"
#include "TDatabasePDG.h"
#include "TROOT.h"
#include "AliMathBase.h"
#include "TLatex.h"
//

const Int_t kColors[6]={1,2,3,4,6,7};
const Int_t kMarkers[6]={20,21,24,25,24,25};

void MakeLocalDistortionCurrentTree(Int_t npointsZ=20000, Int_t id=0);

void makeCurrentCorrection(Int_t action=0, Int_t arg0=5000, Int_t arg1=0){
  //
  //
  //
  if (action==0) MakeLocalDistortionCurrentTree(arg0,arg1);
}

void SetGraphTDRStyle(TGraph * graph){
  graph->GetXaxis()->SetLabelSize(0.08);
  graph->GetXaxis()->SetTitleSize(0.08);
  graph->GetYaxis()->SetLabelSize(0.08);
  graph->GetYaxis()->SetTitleSize(0.08);
  graph->GetXaxis()->SetNdivisions(510);
  graph->GetYaxis()->SetNdivisions(505);
}

void MakeLocalDistortionCurrentTree(Int_t npointsZ, Int_t id){
  //
  // Macro to make trees with local distortions 
  // Results are later visualized in the function DrawLocalDistortionPlots()
  //
  /*
    Int_t npointsZ=1000;
    Int_t id=0;
  */
  Int_t nitteration=300;
  TTreeSRedirector *pcstream = new TTreeSRedirector("localCurrent.root","update");
  //
  TFile *fCurrent = TFile::Open("SpaceChargeFluc10_1.root");
  TFile *fRef = TFile::Open("SpaceChargeFluc0_1.root");
  //
  AliTPCSpaceCharge3D* distortion = ( AliTPCSpaceCharge3D*)fCurrent->Get("DistRef"); 
  AliTPCSpaceCharge3D* distortionRef = ( AliTPCSpaceCharge3D*)fRef->Get("DistRef"); 
  //
  Int_t nz= distortion->GetInputSpaceCharge3D()->GetZaxis()->GetNbins();
  TH3 *hisOrig = (TH3*)distortionRef->GetInputSpaceCharge3D();
  printf("Make mean histo\n");
  TStopwatch timer;
  //
  for (Int_t iside=0; iside<2; iside++){
    for (Int_t iphi=1; iphi<distortion->GetInputSpaceCharge3D()->GetXaxis()->GetNbins(); iphi+=3){
      for (Int_t ir=1; ir<distortion->GetInputSpaceCharge3D()->GetYaxis()->GetNbins(); ir+=3){
	Double_t sum=0, sumW=0;
	if (iside==0) for (Int_t iz2=1; iz2<nz/2; iz2++)   {sum+= hisOrig->GetBinContent(iphi,ir,iz2); sumW++;}
	if (iside==1) for (Int_t iz2=nz/2; iz2<nz;  iz2++) {sum+= hisOrig->GetBinContent(iphi,ir,iz2); sumW++;}
	//
	if (iside==0) for (Int_t iz=1; iz<nz/2; iz++) hisOrig->SetBinContent(iphi,ir,iz,sum/sumW);
	if (iside==1) for (Int_t iz=nz/2; iz<=nz; iz++) hisOrig->SetBinContent(iphi,ir,iz,sum/sumW);
      }
    }
  }
  timer.Print();
  printf("Make mean histo\n");
  //
  distortion->InitSpaceCharge3DPoisson(129, 129, 144,nitteration);
  distortionRef->InitSpaceCharge3DPoisson(129, 129, 144,nitteration);
  //
  distortion->AddVisualCorrection(distortion,1);
  distortionRef->AddVisualCorrection(distortionRef,2);
  //
  TVectorD normZR(125), normZRPhi(125), normZZ(125), normZPos(125);
  TVectorD normDZR(125), normDZRPhi(125), normDZZ(125);
  TVectorD normZRChi2(125), normZRPhiChi2(125), normZZChi2(125);
  TVectorD qCurrent(125), qRef(125);
  TVectorD qCurrentInt(125), qRefInt(125);
  
  //
  for (Int_t iz =0; iz<125; iz++){
    for (Int_t iphi=1; iphi<distortion->GetInputSpaceCharge3D()->GetXaxis()->GetNbins(); iphi+=3)
      for (Int_t ir=1; ir<distortion->GetInputSpaceCharge3D()->GetYaxis()->GetNbins(); ir+=3){
	qCurrent[iz]+=	distortion->GetInputSpaceCharge3D()->GetBinContent(iphi,ir,iz+1);
	qRef[iz]+=	distortionRef->GetInputSpaceCharge3D()->GetBinContent(iphi,ir,iz+1);		
      }   
  } 
  //
  for (Int_t iz =0; iz<125; iz++){
    for (Int_t jz =0; jz<125; jz++){
      if (iz<125/2 && jz<=iz){
	qCurrentInt[iz]+=qCurrent[jz];
	qRefInt[iz]+=qRef[jz];
      }
      if (iz>125/2 && jz>=iz){
	qCurrentInt[iz]+=qCurrent[jz];
	qRefInt[iz]+=qRef[jz];
      }
    }
  }
  //    
  //
  for (Int_t iz =0; iz<125; iz++){
    Double_t z0 = -250+iz*4;
    TLinearFitter fitterR(2,"pol1");
    TLinearFitter fitterRPhi(2,"pol1");
    TLinearFitter fitterZ(2,"pol1");
    Double_t xvalue[10]={0};
    for (Int_t ipoint =0; ipoint<npointsZ; ipoint++){      
      Double_t r0   = 95+gRandom->Rndm()*(245-95.);
      Double_t phi0 = gRandom->Rndm()*TMath::TwoPi();      
      if (TMath::Abs(distortion->GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,1,1))>20) continue;
      if (TMath::Abs(distortion->GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,1,2))>20) continue;
      xvalue[0]=distortion->GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,0,2);
      fitterR.AddPoint(xvalue,distortion->GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,0,1));
      xvalue[0]=distortion->GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,1,2);
      fitterRPhi.AddPoint(xvalue,distortion->GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,1,1));
      xvalue[0]=distortion->GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,2,2);
      fitterZ.AddPoint(xvalue,distortion->GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,2,1));
    }    
    fitterR.Eval();
    fitterRPhi.Eval();
    fitterZ.Eval();    
    normZR[iz]=fitterR.GetParameter(1);
    normZRPhi[iz]=fitterRPhi.GetParameter(1);
    normZZ[iz]=fitterZ.GetParameter(1);
    normZRChi2[iz]=TMath::Sqrt(fitterR.GetChisquare()/fitterR.GetNpoints());    
    normZRPhiChi2[iz]=TMath::Sqrt(fitterRPhi.GetChisquare()/fitterRPhi.GetNpoints());    
    normZZChi2[iz]=TMath::Sqrt(fitterZ.GetChisquare()/fitterZ.GetNpoints());    
    //
    if (iz>0){
      normDZR[iz]=(normZR[iz]-normZR[iz-1]);
      normDZRPhi[iz]=(normZRPhi[iz]-normZRPhi[iz-1]);
    }
    normZPos[iz]=z0;    
  }
  {    
    (*pcstream)<<"meanCurrent"<<
      "id="<<id<<                       // lookup ID
      "normZPos.="<<&normZPos<<         // zposition
      //
      "qCurrent.="<<&qCurrent<<         // current measuremn 
      "qRef.="<<&qRef<<                 // current in refenece sample
      "qCurrentInt.="<<&qCurrentInt<<   // integral of current
      "qRefInt.="<<&qRefInt<<           // integral of current 
      //
      //
      //
      "normZR.="<<&normZR<<             // mult. scaling to minimize R distortions
      "normDZR.="<<&normDZR<<           // mult. scaling to minimize R distortions
      "normZRPhi.="<<&normZRPhi<<       // mult. scaling 
      "normDZRPhi.="<<&normDZRPhi<<     // mult. scaling 
      "normZZ.="<<&normZZ<<
      //
      "normZRChi2.="<<&normZRChi2<<            // mult. scaling to minimize R distortions
      "normZRPhiChi2.="<<&normZRPhiChi2<<      // mult. scaling 
      "normZZChi2.="<<&normZZChi2<<
      "\n";
  }
  delete pcstream;
  
  pcstream = new TTreeSRedirector("localCurrent.root","update");
  TTree * treeNormZ= (TTree*)pcstream->GetFile()->Get("meanCurrent");
  TGraphErrors * grZRfit= TStatToolkit::MakeGraphErrors( treeNormZ, "normZR.fElements:normZPos.fElements","",25,2,0.5);
  TGraphErrors * grZRPhifit= TStatToolkit::MakeGraphErrors( treeNormZ, "normZRPhi.fElements:normZPos.fElements","",25,4,0.5);
  grZRfit->Draw("alp");
  grZRPhifit->Draw("lp");
}


void Fit(){
  //
  // Not good  rsponse should be more complicated
  //
  TTreeSRedirector *pcstream = new TTreeSRedirector("localCurrent.root","update");
  TTree * treeCurrent= (TTree*)pcstream->GetFile()->Get("meanCurrent");
  //
  TVectorD* pnormZR=0, *pnormZRPhi=0, *pnormZZ=0, *pnormZPos=0;
  TVectorD* pqCurrent=0, *pqRef=0;
  TVectorD* pqCurrentInt=0, *pqRefInt=0;
  //
  treeCurrent->SetBranchAddress("normZR.",&pnormZR);
  treeCurrent->SetBranchAddress("normZPos.",&pnormZPos);
  treeCurrent->SetBranchAddress("normZRPhi.",&pnormZRPhi);
  treeCurrent->SetBranchAddress("qCurrent.",&pqCurrent);
  treeCurrent->SetBranchAddress("qRef.",&pqRef);
  treeCurrent->SetBranchAddress("qCurrentInt.",&pqCurrentInt);
  treeCurrent->SetBranchAddress("qRefInt.",&pqRefInt);
  Int_t entries = treeCurrent->GetEntries();
  //
  //
  for (Double_t sigma=1; sigma<40; sigma++){    
    for (Int_t entry=0; entry<entries; entry++){
      treeCurrent->GetEntry(entry);
      //
      TVectorD vecCurrentRefInt(125);
      for (Int_t i=0; i<125; i++){
	Double_t sumW=0;
	for (Int_t j=0; j<125; j++){
	  if (((*pnormZPos)[i]*(*pnormZPos)[j])<0) continue;
	  Double_t weight = 0;
	  if ((*pnormZPos)[i]<0) weight=0.5*(TMath::Erf(Double_t(i-j)/Double_t(sigma))+1);
	  if ((*pnormZPos)[i]>0) weight=0.5*(TMath::Erf(Double_t(j-i)/Double_t(sigma))+1);
	  vecCurrentRefInt[i]+=weight*(*pqCurrent)[j]/(*pqRef)[j];	  
	  sumW+=weight;
	}
	vecCurrentRefInt[i]/=sumW;
      }
      (*pcstream)<<"sigmaScan"<<
	"entry="<<entry<<
	"sigma="<<sigma<<
	"normZPos.="<<pnormZPos<<
	"normZR.="<<pnormZR<<
	"normZRPhi.="<<pnormZRPhi<<
	"vecCurrent.="<<&vecCurrentRefInt<<
	"\n";
    }
  }
  delete pcstream; 
  pcstream = new TTreeSRedirector("localCurrent.root","update");
  TTree * treeScan= (TTree*)pcstream->GetFile()->Get("sigmaScan"); 
  treeCurrent= (TTree*)pcstream->GetFile()->Get("meanCurrent");
  treeScan->SetMarkerStyle(25);
  treeScan->SetMarkerSize(0.4); 
  treeScan->Draw("vecCurrent.fElements:normZR.fElements:sigma","abs(normZPos.fElements-100)<20&&sigma>10","colz");

}

void FitLinear(){
  //
  // deltaX_i = A_ij*deltaI_j
  // 
  TTreeSRedirector *pcstream = new TTreeSRedirector("localCurrent.root","update");
  TTree * treeCurrent= (TTree*)pcstream->GetFile()->Get("meanCurrent");
  //
  TVectorD* pnormZR=0, *pnormZRPhi=0, *pnormZZ=0, *pnormZPos=0;
  TVectorD* pqCurrent=0, *pqRef=0;
  TVectorD* pqCurrentInt=0, *pqRefInt=0;
  //
  treeCurrent->SetBranchAddress("normZR.",&pnormZR);
  treeCurrent->SetBranchAddress("normZPos.",&pnormZPos);
  treeCurrent->SetBranchAddress("normZRPhi.",&pnormZRPhi);
  treeCurrent->SetBranchAddress("qCurrent.",&pqCurrent);
  treeCurrent->SetBranchAddress("qRef.",&pqRef);
  treeCurrent->SetBranchAddress("qCurrentInt.",&pqCurrentInt);
  treeCurrent->SetBranchAddress("qRefInt.",&pqRefInt);
  Int_t entries = treeCurrent->GetEntries();
  //
  //
  //
  Int_t nParamSide=62;               // number of z bins on 1 side
  Int_t group=3;                     // grouping of bins
  Int_t nParams=nParamSide/group;    // parameters
  Int_t nParamsFit=nParams*nParams;  //
  TLinearFitter fitter(nParamsFit+1,TString::Format("hyp%d",nParamsFit));
  TVectorD xVector(nParamsFit+1);
  TVectorD pVector(nParamsFit+1);
  //
  for (Int_t ievent=0; ievent<entries; ievent++){
    treeCurrent->GetEntry(ievent);
    for (Int_t iz=0; iz<nParamSide; iz++){      
      Int_t dPar=iz/group;
      if (dPar>nParams-1) dPar=nParams-1;
      // 1.) clear X vectors
      for (Int_t ipar=0; ipar<nParamsFit; ipar++) xVector[ipar]=0;
      // 2.) set   X vector of interest
      for (Int_t cPar=0; cPar<nParams; cPar++){
	xVector[dPar*nParams+cPar] += (*pqCurrent)[cPar*group]/(*pqRef)[cPar*group];
 	xVector[dPar*nParams+cPar] += (*pqCurrent)[TMath::Min(cPar*group+1,nParamSide)]/(*pqRef)[TMath::Min(cPar*group+1,nParamSide)];
 	xVector[dPar*nParams+cPar] += (*pqCurrent)[TMath::Max(cPar*group-1,0)]/(*pqRef)[TMath::Max(cPar*group-1,0)];
      }
      Double_t val = 0;
      val+=(*pnormZR)[dPar*group];
      val+=(*pnormZR)[TMath::Max(dPar*group-1,0)];
      val+=(*pnormZR)[TMath::Min(dPar*group,nParamSide)];
      val/=3.;
      fitter.AddPoint(xVector.GetMatrixArray(),val,0.0035);
    }
  }
  for (Int_t cPar=1; cPar<nParams-1; cPar++){
    //
    for (Int_t ipar=0; ipar<nParamsFit; ipar++) xVector[ipar]=0;
    xVector[(cPar-1)*nParams+cPar-1]=0.5;
    xVector[(cPar)*nParams+cPar]    =-1;
    xVector[(cPar+1)*nParams+cPar+1]=0.5;
    fitter.AddPoint(xVector.GetMatrixArray(),0,0.05);
  }
  fitter.Eval();
  //
  //
  TVectorD fitVector(nParamsFit);
  TVectorD czVector(nParamsFit);
  TVectorD fitVectorErr(nParamsFit);
  fitter.GetParameters(fitVector);
  for (Int_t iPar=0; iPar<nParamsFit; iPar++)  {
    pVector[iPar]=iPar;
    czVector[iPar]=250.*((iPar-1)%nParams)/nParams;
    fitVectorErr[iPar]=fitter.GetParError(iPar);
  }
  Double_t chi2= TMath::Sqrt(fitter.GetChisquare()/fitter.GetNpoints());
  printf("Chi2=%f\n",chi2);
  TGraphErrors * gr  = new TGraphErrors(nParamsFit,  pVector.GetMatrixArray(), fitVector.GetMatrixArray(), 0,fitVectorErr.GetMatrixArray());
  gr->SetMarkerStyle(25); gr->SetMarkerSize(0.3);
  gr->Draw("alp");
  
  TGraphErrors * vgraphs[20]={0};
  for (Int_t ipar=0; ipar<20; ipar++){
    vgraphs[ipar]=new TGraphErrors(nParams, &(czVector.GetMatrixArray()[ipar*nParams+1]), &(fitVector.GetMatrixArray()[ipar*nParams+1]), 0,  &(fitVectorErr.GetMatrixArray()[ipar*nParams+1])); 
    vgraphs[ipar]->GetXaxis()->SetTitle("z_{I} (cm)");
    vgraphs[ipar]->GetYaxis()->SetTitle("#it{A}_{i_{I},j_{#DeltaR}}");
    vgraphs[ipar]->SetMinimum(0);
    vgraphs[ipar]->SetMaximum(0.1);
    SetGraphTDRStyle(vgraphs[ipar]);
  }
  
  TCanvas * canvasFit = new TCanvas("canvasFit","canvasFit",700,700); 
  canvasFit->SetRightMargin(0.05);
  canvasFit->SetLeftMargin(0.15);
  canvasFit->SetBottomMargin(0.18);
  canvasFit->Divide(1,2,0,0);
  TLegend * legend0 = new TLegend(0.4,0.4,0.95,0.95,"Current fluctuation correction matrix. #DeltaR_{zi}=A_{ij}I_{zj}");
  legend0->SetBorderSize(0);
  TLegend * legend1 = new TLegend(0.4,0.5,0.95,0.95,"Current fluctuation correction matrix. #DeltaR_{zi}=A_{ij}I_{zj}");
  legend1->SetBorderSize(0);

  for (Int_t ipar=0; ipar<10; ipar++){
    canvasFit->cd(ipar/5+1)->SetTicks(3,3);
    vgraphs[ipar*2]->SetMarkerStyle(kMarkers[ipar%5]);
    vgraphs[ipar*2]->SetMarkerColor(kColors[ipar%5]);
    vgraphs[ipar*2]->SetLineColor(kColors[ipar%5]);
    if (ipar%5==0) vgraphs[ipar*2]->Draw("alp");
    vgraphs[ipar*2]->Draw("lp");
    if (ipar<5) legend0->AddEntry(vgraphs[ipar*2],TString::Format("z_{#DeltaR}=%1.0f (cm)",250.*ipar/10.),"p");
    if (ipar>=5) legend1->AddEntry(vgraphs[ipar*2],TString::Format("z_{#DeltaR}=%1.0f (cm)",250.*ipar/10.),"p");
  }	
  legend0->SetTextSize(0.05);
  legend1->SetTextSize(0.05);
  //
  canvasFit->cd(1);
  legend0->Draw();
  canvasFit->cd(2);
  legend1->Draw();
  canvasFit->SaveAs("canvasAIJ_CurrenToDR.pdf");
  canvasFit->SaveAs("canvasAIJ_CurrenToDR.png");
}


void FFTexample(){
  //
  TTreeSRedirector *pcstream = new TTreeSRedirector("localCurrent.root","update");
  //
  TFile *fCurrent = TFile::Open("SpaceChargeFluc10_1.root");
  TFile *fRef = TFile::Open("SpaceChargeFluc0_1.root");
  //  
  AliTPCSpaceCharge3D* distortion = ( AliTPCSpaceCharge3D*)fCurrent->Get("DistRef"); 
  AliTPCSpaceCharge3D* distortionRef = ( AliTPCSpaceCharge3D*)fRef->Get("DistRef"); 

  distortion->AddVisualCorrection(distortion,1);
  distortionRef->AddVisualCorrection(distortionRef,2);
  //
  TF2 f2("f2","AliTPCCorrection::GetCorrXYZ(x,0,y,0,1)-1.07*AliTPCCorrection::GetCorrXYZ(x,0,y,0,2)",85,245,0,250);
  TF2 f2c("f2c","AliTPCCorrection::GetCorrXYZ(x,0,y,0,1)-1.07*AliTPCCorrection::GetCorrXYZ(x,0,y,0,2)+sqrt(12.)*(rndm-0.5)*0.04",85,245,0,250);
  f2.SetNpx(16); f2.SetNpy(32);
  f2c.SetNpx(16); f2c.SetNpy(32);
  f2.Draw("colz");
  f2c.Draw("colz");
  
  TH2D * his2D = (TH2D*)f2.GetHistogram();
  TH2D * his2Dc = (TH2D*)f2c.GetHistogram();
  TH2D * his2DSmooth=new TH2D(*his2Dc);
  Double_t sigma = 3;
  Int_t nx=his2D->GetXaxis()->GetNbins();
  Int_t ny=his2D->GetYaxis()->GetNbins();
  //
  for (Int_t ibin=2; ibin<=nx-1;ibin++){
    for (Int_t jbin=2; jbin<ny-1;jbin++){
      TLinearFitter fitter(4,"hyp3");
      for (Int_t di=-3; di<=3; di++)
	for (Int_t dj=-3; dj<=3; dj++){
	  if (ibin+di<=1) continue;
	  if (jbin+dj<=1) continue;
	  if (ibin+di>=nx) continue;
	  if (jbin+dj>=ny) continue;
	  Double_t x[2]={di,dj,dj*dj};
	  Double_t w = 1./(TMath::Gaus(di/sigma)*TMath::Gaus(dj/sigma));	  
	  Double_t y = his2Dc->GetBinContent(ibin+di,jbin+dj);
	  fitter.AddPoint(x,y,w);
	}
      fitter.Eval();
      printf("%d\t%d\t%3.3f\n",ibin,jbin,10*(fitter.GetParameter(0)-his2Dc->GetBinContent(ibin,jbin)));
      his2DSmooth->SetBinContent(ibin,jbin,fitter.GetParameter(0));
    } 
  }
  TH2D hisDiff=*his2DSmooth;
  hisDiff.Add(his2D,-1);
  

  TH1 * hisMag = his2D->FFT(0,"MAG");
  TH1 * hisMagc = his2Dc->FFT(0,"MAG");
  TH1 * hisPhi = his2D->FFT(0,"PH");
  TH1 * hisPhic = his2Dc->FFT(0,"PH");
  hisMag->Draw("surf2");
  hisMagc->Draw("surf2");

    
}

/*
void MakeFFT(){
  //
  //
  //
  Int_t nSize = 125;
  TVirtualFFT *fft_own = TVirtualFFT::FFT(1, &nSize, "R2C ES K");
  fft_own->SetPoints(normZR->GetMatrixArray());
  fft_own->Transform();

  TH1 *hre = 0, *hco=0;  
  hre = TH1::TransformHisto(fft_own, hre, "RE");
  hco = TH1::TransformHisto(fft_own, hco, "C");


  hr->SetTitle("Real part of the 3rd (array) tranfsorm");
  hr->Draw();
  hr->SetStats(kFALSE);
  hr->GetXaxis()->SetLabelSize(0.05);
  hr->GetYaxis()->SetLabelSize(0.05);
  c1_6->cd();
  TH1 *him = 0;
  him = TH1::TransformHisto(fft_own, him, "IM");
  him->SetTitle("Im. part of the 3rd (array) transform");
  him->Draw();
  him->SetStats(kFALSE);
  him->GetXaxis()->SetLabelSize(0.05);
  him->GetYaxis()->SetLabelSize(0.05);
  //
}
*/
