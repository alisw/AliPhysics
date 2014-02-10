/*
.x $ALICE_ROOT/TPC/Upgrade/macros/NimStyle.C

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
#include "AliTPCCorrectionLookupTable.h"
//

const Int_t kColors[6]={1,2,3,4,6,7};
const Int_t kMarkers[6]={20,21,24,25,24,25};
TObjArray garrayFit(3);


void MakeLocalDistortionCurrentTree(Int_t npointsZ=20000, Int_t id=0);
void MakeSmoothKernelStudy(Int_t npoints, Int_t nkernels, Int_t mapID);

void makeCurrentCorrection(Int_t action=0, Int_t arg0=5000, Int_t arg1=0, Int_t arg2=0){
  //
  //
  //
  if (action==0) MakeLocalDistortionCurrentTree(arg0,arg1);
  if (action==1) MakeSmoothKernelStudy(arg0,arg1,arg2);
 
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

void MakeSmoothKernelStudy(Int_t npoints, Int_t nkernels, Int_t mapID){
  //
  // Compare the input and reconstructed distortion maps
  // Input files are expected to have predefined names
  // 
  // Values and smoothed vaules using differnt smoothing algorithms are used
  // Results of given studies will be used to define "optimal" smoothing algorithm
  // and optimal Kernel parameters
  // 
  //
  gRandom->SetSeed(mapID);
  const Double_t cutMaxDist=3;
  const Double_t cutMinDist=0.00001;
  const Int_t kMinPoints=10;
  TFile *foutput = TFile::Open("MeasureResidual.root");
  TFile *finput = TFile::Open("RealResidualScaled.root");
  AliTPCCorrectionLookupTable *mapIn = (AliTPCCorrectionLookupTable *)finput->Get("map");
  AliTPCCorrectionLookupTable *mapOut = (AliTPCCorrectionLookupTable *)foutput->Get("map");
  AliTPCCorrection::AddVisualCorrection(mapIn,1);    // input==1
  AliTPCCorrection::AddVisualCorrection(mapOut,2);   // output (reconstructed==2)
  TF1 * fin = new TF1("fin","AliTPCCorrection::GetCorrXYZ(x,5,10,1,1)",85,245);
  TF1 * fout = new TF1("fout","AliTPCCorrection::GetCorrXYZ(x,5,10,1,2)",85,245);  
  TTreeSRedirector * pcstream = new TTreeSRedirector("smoothLookupStudy.root","recreate");
  //
  //
  // 1. Generate random point of interest
  //
  TVectorD  sigmaKernelR(nkernels);
  TVectorD  sigmaKernelRPhi(nkernels);
  TVectorD  sigmaKernelZ(nkernels);
  //
  TVectorD  fitValuesOut(nkernels);
  TVectorD  fitValuesOutR(nkernels);
  TVectorD  fitValuesIn(nkernels);
  TVectorD  fitValuesInR(nkernels);
  //
  TVectorD  fitValuesOutErr(nkernels);
  TVectorD  fitValuesOutChi2(nkernels);
  TVectorD  fitSumW(nkernels);
  TVectorD  fitValuesOutN(nkernels);
  //
  TLinearFitter *fittersOut[nkernels];
  TLinearFitter *fittersOutR[nkernels];
  TLinearFitter *fittersIn[nkernels];
  TLinearFitter *fittersInR[nkernels];
  for (Int_t ipoint=0; ipoint<npoints; ipoint++){ 
    if (ipoint%10==0) printf("%d\n",ipoint);
    for (Int_t i=0; i<nkernels; i++) {
      fittersOut[i]=new TLinearFitter(7,"hyp6");
      fittersOutR[i]=new TLinearFitter(7,"hyp6");
      fittersIn[i]=new TLinearFitter(7,"hyp6");
      fittersInR[i]=new TLinearFitter(7,"hyp6");
    }
    Double_t phi = gRandom->Rndm()*TMath::TwoPi();
    Double_t r   = 85+gRandom->Rndm()*(245-85);
    Double_t theta   = -0.9+gRandom->Rndm()*1.8;
    Double_t z=r*theta;     
    //
    Double_t x = r*TMath::Cos(phi);
    Double_t y = r*TMath::Sin(phi);
    Double_t drphiInput = AliTPCCorrection::GetCorrXYZ(x,y,z,1,1);   // 
    Double_t drphiOutput = AliTPCCorrection::GetCorrXYZ(x,y,z,1,2);  //
    Double_t drInput = AliTPCCorrection::GetCorrXYZ(x,y,z,0,1);
    Double_t drOutput = AliTPCCorrection::GetCorrXYZ(x,y,z,0,2);
    if (TMath::Abs(drphiOutput)<cutMinDist) continue; // beter condition needed
    if (TMath::Abs(drphiInput)<cutMinDist) continue; // beter condition needed
    if (TMath::Abs(drphiOutput)>cutMaxDist) continue; // beter condition needed
    if (TMath::Abs(drphiInput)>cutMaxDist) continue; // beter condition needed
    //
    for (Int_t i=0; i<nkernels; i++) {
      sigmaKernelR[i]    = 0.5+30.*TMath::Power(gRandom->Rndm(),2); 
      sigmaKernelRPhi[i] = 0.5+30.*TMath::Power(gRandom->Rndm(),2);
      sigmaKernelZ[i]    = 0.5+30.*TMath::Power(gRandom->Rndm(),2);
      fitSumW[i]=0;
    }
    for (Int_t idphi=-12; idphi<=12; idphi++){
      Double_t dphi=idphi*TMath::Pi()/90.;            // we need to know actual segentation of the histograms - lookup should by multiple
      for (Double_t dR=-50; dR<=50.; dR+=5){
	for (Double_t dZ=-50; dZ<=50.; dZ+=5){
	  if (r+dR<85) continue;
	  if (r+dR>245) continue;
	  if (z+dZ<-240) continue;
	  if (z+dZ>240) continue;
	  if (z*(z+dZ)<0) continue;
	  //
	  Double_t x2=(r+dR)*TMath::Cos(phi+dphi);
	  Double_t y2=(r+dR)*TMath::Sin(phi+dphi);
	  Double_t z2=z+dZ;
	  Double_t drphiInput2=AliTPCCorrection::GetCorrXYZ(x2,y2,z2,1,1);
	  Double_t drInput2=AliTPCCorrection::GetCorrXYZ(x2,y2,z2,0,1);
	  Double_t drphiOutput2=AliTPCCorrection::GetCorrXYZ(x2,y2,z2,1,2);
	  Double_t drOutput2=AliTPCCorrection::GetCorrXYZ(x2,y2,z2,0,2);
	  if (TMath::Abs(drphiOutput2)<cutMinDist) continue; // hard cut beter condition needed
	  if (TMath::Abs(drphiOutput2)>cutMaxDist) continue; // beter condition needed
	  if (TMath::Abs(drphiInput2)<cutMinDist) continue; // hard cut beter condition needed
	  if (TMath::Abs(drphiInput2)>cutMaxDist) continue; // beter condition needed

	  Double_t xfit[7]={dphi,dphi*dphi,dR,dR*dR,dZ,dZ*dZ};
	  for (Int_t i=0; i<nkernels; i++) {
	    Double_t weight=1;
	    weight*=TMath::Gaus(dphi*r,0,sigmaKernelRPhi[i]);
	    weight*=TMath::Gaus(dR,0,sigmaKernelR[i]);
	    weight*=TMath::Gaus(dZ,0,sigmaKernelZ[i]);
	    weight+=0.00000001;
	    fitSumW[i]+=weight;
	    fittersOut[i]->AddPoint(xfit,drphiOutput2,0.1/weight);
	    fittersOutR[i]->AddPoint(xfit,drOutput2,0.1/weight);
	    //
	    fittersIn[i]->AddPoint(xfit,drphiInput2,0.1/weight);
	    fittersInR[i]->AddPoint(xfit,drInput2,0.1/weight);
	  }
	}
      }
    }   
    
    for (Int_t ifix=0; ifix<=1; ifix++){
      Bool_t isOK=kTRUE;
      for (Int_t i=0; i<nkernels; i++) {
	if( fittersOut[i]->GetNpoints() < kMinPoints) {
	  isOK=kFALSE;
	  break;
	}
	if (fitSumW[i]<0.01/*kMinWeight*/){
	  isOK=kFALSE;
	  break;
	}
	if (ifix==1){
	  fittersOut[i]->FixParameter(4,0);
	  fittersOut[i]->FixParameter(5,0);
	  fittersOut[i]->FixParameter(6,0);
	  fittersOutR[i]->FixParameter(4,0);
	  fittersOutR[i]->FixParameter(5,0);
	  fittersOutR[i]->FixParameter(6,0);
	  fittersIn[i]->FixParameter(4,0);
	  fittersIn[i]->FixParameter(5,0);
	  fittersIn[i]->FixParameter(6,0);
	  fittersInR[i]->FixParameter(4,0);
	  fittersInR[i]->FixParameter(5,0);
	  fittersInR[i]->FixParameter(6,0);
	}
	fittersOut[i]->Eval();
	fittersOutR[i]->Eval();
	fittersIn[i]->Eval();
	fittersInR[i]->Eval();
	//
	fitValuesOut[i]=fittersOut[i]->GetParameter(0);
	fitValuesOutR[i]=fittersOutR[i]->GetParameter(0);
	fitValuesIn[i]=fittersIn[i]->GetParameter(0);
	fitValuesInR[i]=fittersInR[i]->GetParameter(0);
	//
	fitValuesOutErr[i]=fittersOut[i]->GetParError(0);
	fitValuesOutChi2[i]=TMath::Sqrt(fittersOut[i]->GetChisquare()/fitSumW[i]);
	fitValuesOutN[i]=fittersOut[i]->GetNpoints();
      }
      if (isOK){
      (*pcstream)<<"smoothLookup"<<
	"ipoint"<<ipoint<<
	"mapID="<<mapID<<
	"ifix="<<ifix<<
	// coordinates
	"x="<<x<<
	"y="<<y<<
	"z="<<z<<
	"phi="<<phi<<
	"r="<<r<<
	"z="<<z<<
	// Input output values
	"drphiInput="<<drphiInput<<       // input lookup tables disotrions
	"drphiOutput="<<drphiOutput<<     // reconstructed lookup tables distortion
	"drInput="<<drInput<<       // input lookup tables disotrions
	"drOutput="<<drOutput<<     // reconstructed lookup tables distortion
	// Smoothed values
	"fitValuesIn.="<<&fitValuesIn<<            // smoothed input values - rphi direction
	"fitValuesInR.="<<&fitValuesInR<<          // smoothed input values - r direction
	"fitValuesOut.="<<&fitValuesOut<<          // smoothed measuured values - rphi
	"fitValuesOutR.="<<&fitValuesOutR<<        // smoothed measuured values -r 
	"fitValuesOutErr.="<<&fitValuesOutErr<<
	"fitValuesOutChi2.="<<&fitValuesOutChi2<<
	"fitValuesOutN.="<<&fitValuesOutN<<
	"fitSumW.="<<&fitSumW<<
	// Kernel sigma
	"sigmaKernelR.="<<&sigmaKernelR<<
	"sigmaKernelRPhi.="<<&sigmaKernelRPhi<<
	"sigmaKernelZ.="<<&sigmaKernelZ<<
	"\n";
      }
    }
    for (Int_t i=0; i<nkernels; i++) {
      delete fittersOut[i];
      delete fittersOutR[i];
      delete fittersIn[i];
      delete fittersInR[i];
    }
  }
  delete pcstream;
  //
  //
  /* 
     TFile *ff = TFile::Open("smoothLookupStudy.root")
     TChain * chain = AliXRDPROOFtoolkit::MakeChain("smooth.list", "smoothLookup", 0,100)

   */
}
void MakeSmoothKernelStudyDraw(Int_t npoints){
  //
  // make figure for kernel smoothing Draw
  //
  // npoints=20000
  TFile *finput = TFile::Open("smoothLookupStudy.root");
  TTree * chain = (TTree*)finput->Get("smoothLookup");  
  chain->SetCacheSize(100000000);
  TH2* hisInputsS[10]={0};
  TH3* hisInputs3D[10]={0};
  TH1* hisSigmaS[10]={0};
  //
  // 1.) Resolution not smoothed make nice plots
  //
  TH2 * his2DBase[10]={0};
  TH1 * hisSigmas[10]={0};
  chain->Draw("10*drphiInput:r>>hisIn(30,85,245,100,-2,2)","","colz");
  his2DBase[0]=(TH2*)chain->GetHistogram()->Clone();
  chain->Draw("10*(drphiOutput-drphiInput):r>>hisOutIn(30,85,245,100,-2,2)","","colz");
  his2DBase[1]=(TH2*)chain->GetHistogram()->Clone();
  his2DBase[0]->FitSlicesY(0,0,-1,0,"QNR",&garrayFit);
  hisSigmas[0] = (TH1*)garrayFit.At(2)->Clone();
  his2DBase[1]->FitSlicesY(0,0,-1,0,"QNR",&garrayFit);
  hisSigmas[1] = (TH1*)garrayFit.At(2)->Clone();
  //
  TCanvas *canvasInOut = new TCanvas("deltaInOut","deltaInOut",600,500);
  for (Int_t i=0; i<2; i++) {
    hisSigmas[i]->SetMinimum(0); 
    hisSigmas[i]->SetMaximum(1.5); 
    hisSigmas[i]->SetMarkerStyle(kMarkers[i]);
    hisSigmas[i]->SetMarkerColor(kColors[i]);
    hisSigmas[i]->GetXaxis()->SetTitle("R (cm)");
    hisSigmas[i]->GetYaxis()->SetTitle("#Delta_{R#phi} (mm)");
    if (i==0) hisSigmas[i]->Draw("");
    hisSigmas[i]->Draw("same");
  }
  TLegend * legend = new TLegend(0.4,0.5,0.89,0.89,"Residual #Delta_{r#phi}");
  legend->SetBorderSize(0);
  legend->AddEntry(hisSigmas[0],"#Delta_{current}-k#Delta_{mean}");
  legend->AddEntry(hisSigmas[1],"(#Delta_{current}-k#Delta_{mean})_{align}-(#Delta_{current}-k#Delta_{mean})");
  legend->Draw();
  canvasInOut->SaveAs("canvasDistortionAlignmentInOut.pdf");
  canvasInOut->SaveAs("canvasDistortionAlignmentInOut.png");
  //
  // 1.b) phi modulation plots
  //
  TCanvas * canvasPhi= new TCanvas("canvasPhi","canvasPhi",800,600);
  canvasPhi->Divide(1,3,0,0);
  gStyle->SetOptTitle(1);
  chain->SetAlias("sector","9*phi/pi");
  {
    chain->SetMarkerStyle(25);
    chain->SetMarkerSize(0.3);
    canvasPhi->cd(1);
    chain->Draw("(drphiOutput-drphiInput):sector","abs((drphiOutput-drphiInput))<0.2&&abs(r-100)<20","",npoints);
    canvasPhi->cd(2);
    chain->Draw("(drphiOutput-drphiInput):sector","abs((drphiOutput-drphiInput))<0.2&&abs(r-160)<25","",npoints);
    canvasPhi->cd(3);
    chain->Draw("(drphiOutput-drphiInput):sector","abs((drphiOutput-drphiInput))<0.2&&abs(r-220)<30","",npoints);
  }
  canvasPhi->SaveAs("canvasDistortionPhiSlice.pdf");

  TCanvas * canvasDistortionSliceHisto= new TCanvas("canvasDistortionSliceHisto","canvasDistortionSliceHisto",800,600);
  canvasDistortionSliceHisto->Divide(1,3,0,0);
  gStyle->SetOptTitle(1);
  chain->SetAlias("sector","9*phi/pi");
  {
    chain->SetMarkerStyle(25);
    chain->SetMarkerSize(0.3);
    canvasDistortionSliceHisto->cd(1);
    chain->Draw("(drphiOutput-drphiInput):sector-int(sector)>>hisSec100(40,0,1,50,-0.2,0.2)","abs((drphiOutput-drphiInput))<0.2&&abs(r-110)<30","colz");
    canvasDistortionSliceHisto->cd(2);
    chain->Draw("(drphiOutput-drphiInput):sector-int(sector)>>hisSec160(40,0,1,50,-0.2,0.2)","abs((drphiOutput-drphiInput))<0.2&&abs(r-160)<25","colz");
    canvasDistortionSliceHisto->cd(3);
    chain->Draw("(drphiOutput-drphiInput):sector-int(sector)>>hisSec2200(40,0,1,50,-0.2,0.2)","abs((drphiOutput-drphiInput))<0.2&&abs(r-220)<30","colz");
  }
  canvasDistortionSliceHisto->SaveAs("canvasDistortionSectorSliceHisto.pdf");


  //
  //
  // 2.) Draw plot resolution as fucntion of the kernel smooth sigma
  //      a.) Smoothed input- input
  //      b.) Smoothed
  TObjArray fitResols(1000);
  const char* varSmooth[4]={"Smoothed(Input,Pol1)-Input","Smoothed(Output,Pol1)-Smoothed(Input,Pol1)","Smoothed(Output,Pol2)-Input","Smoothed(Output,Pol1)-Input"};
  //
  chain->Draw("10*(fitValuesIn.fElements-drphiInput):sqrt(sigmaKernelR.fElements**2+sigmaKernelRPhi.fElements**2+sigmaKernelZ.fElements**2):r>>hisSmoothDiff0(5,83,245, 15,5,40,100,-2.2,2.2)","ifix==1&&fitValuesOutErr.fElements<0.2","goff",npoints);
  hisInputs3D[0]= (TH3*)chain->GetHistogram()->Clone();
  chain->Draw("10*(fitValuesOut.fElements-fitValuesIn.fElements):sqrt(sigmaKernelR.fElements**2+sigmaKernelRPhi.fElements**2+sigmaKernelZ.fElements**2):r>>hisSmoothDiff1(5,83,245, 15,5,40,100,-2.2,2.2)","ifix==1&&fitValuesOutErr.fElements<0.2","goff",npoints);
  hisInputs3D[1]= (TH3*)chain->GetHistogram()->Clone();
  chain->Draw("10*(fitValuesOut.fElements-drphiInput):sqrt(sigmaKernelR.fElements**2+sigmaKernelRPhi.fElements**2+sigmaKernelZ.fElements**2):r>>hisSmoothDiff2(5,83,245, 15,5,40,100,-2.2,2.2)","ifix==0&&fitValuesOutErr.fElements<0.2","goff",npoints);
  hisInputs3D[2]= (TH3*)chain->GetHistogram()->Clone();
  chain->Draw("10*(fitValuesOut.fElements-drphiInput):sqrt(sigmaKernelR.fElements**2+sigmaKernelRPhi.fElements**2+sigmaKernelZ.fElements**2):r>>hisSmoothDiff3(5,83,245, 15,5,40,100,-2.2,2.2)","ifix==1&&fitValuesOutErr.fElements<0.2","goff",npoints);
  hisInputs3D[3]= (TH3*)chain->GetHistogram()->Clone();
  //
  //
  //
  TCanvas *canvasSmooth = new TCanvas("DistortionAlignmentSmooth3D","DistortionAlignmentSmooth3D",800,800);
  canvasSmooth->Divide(2,2,0,0);
  gStyle->SetOptStat(0);  
  for (Int_t itype=0; itype<4; itype++){
    canvasSmooth->cd(itype+1);
    TLegend * legend = new TLegend(0.2,0.7,0.9,0.99,varSmooth[itype]);
    legend->SetBorderSize(0);
    for (Int_t ibin=0; ibin<5; ibin++){
      hisInputs3D[itype]->GetXaxis()->SetRange(ibin+1,ibin+1);
      Double_t radius=hisInputs3D[itype]->GetXaxis()->GetBinCenter(ibin+1);
      TH2 * his2D= (TH2*)hisInputs3D[itype]->Project3D("zy");
      his2D->FitSlicesY(0,0,-1,0,"QNR",&garrayFit);
      delete his2D;
      TH1* his1D = (TH1*)garrayFit.At(2)->Clone();
      fitResols.AddLast(his1D);
      his1D->SetMarkerColor(kColors[ibin%6]);
      his1D->SetMarkerStyle(kMarkers[ibin%6]);
      his1D->SetMinimum(0);
      his1D->SetMaximum(0.75);
      his1D->GetXaxis()->SetTitle("#sigma_{kernel} (cm) in 3D (r,r#phi,z)");
      his1D->GetYaxis()->SetTitle("#sigma_{r#phi} (mm)");
      if (ibin==0) his1D->Draw();
      his1D->Draw("same");
      legend->AddEntry(his1D,TString::Format("R=%2.0f (cm)",radius));
    }
    legend->Draw();
  }
  canvasSmooth->SaveAs("DistortionAlignmentSmooth3D.pdf");
  canvasSmooth->SaveAs("DistortionAlignmentSmooth3D.pdf");



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
  f2.SetNpx(32); f2.SetNpy(32);
  f2c.SetNpx(32); f2c.SetNpy(32);
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
	  Double_t x[3]={di,dj,dj*dj};
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
