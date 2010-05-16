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
  TPC Kalman filter implementation for TPC sector alignment.
  Output of the AliTPCcalibAlign is used as a input for TPC global alignment.
  In AliTPCcalibAlign histograms - track parameter matching at sector boundaries are created.
  Each sector is aligned with 5 neighborhoud (sectors)
    1. Up-down    - 1
    2. Left-right - 4

  Sector alignment parameters are obtained finding the alignment parameters, minimizing
  misalignmnet for all piars fo sectors.

  Global minimization-  MakeGlobalAlign
  

  Example usage:
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  //
  Int_t run=117220;
  .x ConfigCalibTrain.C(run)
  calibDB = AliTPCcalibDB::Instance()

  AliTPCkalmanAlign kalmanAlign("TPC align", "TPC align");  // create the object
  kalmanAlign.ReadAlign("CalibObjects.root");               // read the calibration form file         
  kalmanAlign.MakeGlobalAlign();                            // do kalman alignment
  kalmanAlign.DrawDeltaAlign();                             // make QA plot
  //
  

*/
#include "TMath.h"
#include "TTreeStream.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TClonesArray.h"
#include "TCut.h"
#include "TChain.h"
#include "AliXRDPROOFtoolkit.h"
#include "TLegend.h"
#include "TCanvas.h"

#include "AliTPCcalibDB.h"
#include "AliTPCCalROC.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliAlignObjParams.h"
#include "AliTPCROC.h"
#include "AliTracker.h"
#include "TFile.h"
#include "TLinearFitter.h"
#include "AliTPCcalibAlign.h"
#include "TH1.h"
#include "AliTPCCalPad.h"
#include "AliTPCkalmanAlign.h"
#include "TStatToolkit.h"
#include "AliTPCPreprocessorOnline.h"
#include "TPostScript.h"

AliTPCkalmanAlign::AliTPCkalmanAlign():
  TNamed(),
  fCalibAlign(0),     // kalman alignnmnt
  fOriginalAlign(0),   // original alignment 0 read for the OCDB
  fNewAlign(0),
  fPadTime0(0),
  fFitCEGlobal(0),
  fFitCELocal(0)
{
  //
  // Default constructor
  //
  for (Int_t i=0; i<4; i++){
    fDelta1D[i]=0;
    fCovar1D[i]=0;
  }
}

AliTPCkalmanAlign::AliTPCkalmanAlign(const char* name, const char* title): 
  TNamed(name,title),
  fCalibAlign(0),     // kalman alignnmnt
  fOriginalAlign(0),   // original alignment 0 read for the OCDB
  fNewAlign(0),
  fPadTime0(0),
  fFitCEGlobal(0),
  fFitCELocal(0) 
{
  //
  // Default constructor
  //
  for (Int_t i=0; i<4; i++){
    fDelta1D[i]=0;
    fCovar1D[i]=0;
  }
  fFitCEGlobal = new TObjArray(6); 
  fFitCELocal  = new TObjArray(6); 
  for (Int_t ipar=0; ipar<6;ipar++){
    fFitCEGlobal->AddAt(new TVectorD(36),ipar);
    fFitCELocal->AddAt(new TVectorD(36),ipar);
  }
}

void AliTPCkalmanAlign::ReadAlign(const char *fname){
  //
  // Read old alignment used in the reco
  // and the residual histograms
  // WE ASSUME that the OCDB path is set in the same way as done in the calibration
  //
  TFile fcalib(fname);
  TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
  fCalibAlign=0;
  if (array) fCalibAlign=( AliTPCcalibAlign *)array->FindObject("alignTPC");
  fCalibAlign = ( AliTPCcalibAlign *)fcalib.Get("alignTPC");
  //
  // old alignment used
  AliCDBEntry * cdbEntry= AliCDBManager::Instance()->Get("TPC/Align/Data");
  fOriginalAlign =0;
  if (cdbEntry){
    fOriginalAlign = (TClonesArray*)cdbEntry->GetObject();
  }
  
}

void AliTPCkalmanAlign::BookAlign1D(TMatrixD &param, TMatrixD &covar, Double_t mean, Double_t sigma){
  //
  // Book Align 1D parameters and covariance
  //
  param.ResizeTo(72,1);
  covar.ResizeTo(72,72);
  for (Int_t i=0;i<72;i++){
    param(i,0)=mean;
    for (Int_t j=0;j<72;j++) covar(i,j)=0;      
    covar(i,i)=sigma*sigma;
  }
}


void AliTPCkalmanAlign::UpdateAlign1D(Double_t delta, Double_t sigma, Int_t s1, Int_t s2, TMatrixD &vecXk, TMatrixD &covXk){
  //
  // Update 1D kalman filter
  //
  const Int_t knMeas=1;
  const Int_t knElem=72;
  TMatrixD mat1(72,72);            // update covariance matrix
  TMatrixD matHk(1,knElem);        // vector to mesurement
  TMatrixD vecYk(knMeas,1);        // Innovation or measurement residual
  TMatrixD matHkT(knElem,knMeas);  // helper matrix Hk transpose
  TMatrixD matSk(knMeas,knMeas);   // Innovation (or residual) covariance
  TMatrixD matKk(knElem,knMeas);   // Optimal Kalman gain
  TMatrixD covXk2(knElem,knElem);  // helper matrix
  TMatrixD covXk3(knElem,knElem);  // helper matrix
  TMatrixD vecZk(1,1);
  TMatrixD measR(1,1);
  vecZk(0,0)=delta;
  measR(0,0)=sigma*sigma;
  //
  // reset matHk
  for (Int_t iel=0;iel<knElem;iel++) 
    for (Int_t ip=0;ip<knMeas;ip++) matHk(ip,iel)=0; 
  //mat1
  for (Int_t iel=0;iel<knElem;iel++) {
    for (Int_t jel=0;jel<knElem;jel++) mat1(iel,jel)=0;
    mat1(iel,iel)=1;
  }
  //
  matHk(0, s1)=1;
  matHk(0, s2)=-1;
  vecYk = vecZk-matHk*vecXk;               // Innovation or measurement residual
  matHkT=matHk.T(); matHk.T();
  matSk = (matHk*(covXk*matHkT))+measR;    // Innovation (or residual) covariance
  matSk.Invert();
  matKk = (covXk*matHkT)*matSk;            //  Optimal Kalman gain
  vecXk += matKk*vecYk;                    //  updated vector 
  covXk2= (mat1-(matKk*matHk));
  covXk3 =  covXk2*covXk;          
  covXk = covXk3;  
}




void AliTPCkalmanAlign::MakeGlobalAlign(){
  //
  // Combine all pairs of fitters and make global alignemnt
  //
  AliTPCkalmanAlign &kalmanAlign=*this;
  TTreeSRedirector *pcstream = new TTreeSRedirector("AliTPCkalmanAlign.root");
  //
  // get ce info
  //
  AliTPCCalPad *padTime0 = AliTPCcalibDB::Instance()->GetPadTime0();
  TVectorD paramCE[72];
  TMatrixD covarCE[72];
  Int_t  statUpDown=0;   // statistic up down
  Int_t  statLeftRight=0;   // statistic left-right
  Float_t chi2;
  for (Int_t isec=0; isec<72; isec++){
    AliTPCCalROC * roc0  = padTime0->GetCalROC(isec);
    roc0->GlobalFit(0,kFALSE,paramCE[isec],covarCE[isec],chi2,0);
    (*pcstream)<<"ceFit"<<
      "isec="<<isec<<
      "p0.="<<&paramCE[isec]<<
      "\n";
  }

  DumpOldAlignment(pcstream);
  const Int_t kMinEntries=400;
  TMatrixD vec[5];
  TMatrixD cov[5];
  kalmanAlign.BookAlign1D(vec[0],cov[0], 0,0.5);
  kalmanAlign.BookAlign1D(vec[1],cov[1], 0,0.5);
  kalmanAlign.BookAlign1D(vec[2],cov[2], 0,0.01);
  kalmanAlign.BookAlign1D(vec[3],cov[3], 0,0.01);
  //
  TVectorD delta(5);
  TVectorD rms(5);
  TVectorD stat(5);
  TH1 * his=0;
  for (Int_t is0=0;is0<72;is0++)
    for (Int_t is1=0;is1<72;is1++){
      //
      //
      his = kalmanAlign.fCalibAlign->GetHisto(AliTPCcalibAlign::kY,is0,is1);
      if (!his) continue;
      if (his->GetEntries()<kMinEntries) continue;
      delta[0]=his->GetMean();
      rms[0]=his->GetRMS();
      stat[0]=his->GetEntries();
      kalmanAlign.UpdateAlign1D(his->GetMean(),his->GetRMS(),is0,is1, vec[0],cov[0]);
      //     
      his = kalmanAlign.fCalibAlign->GetHisto(AliTPCcalibAlign::kZ,is0,is1);
      if (!his) continue;
      delta[1]=his->GetMean();
      rms[1]=his->GetRMS();
      stat[1]=his->GetEntries();
      kalmanAlign.UpdateAlign1D(his->GetMean(),his->GetRMS(),is0,is1, vec[1],cov[1]);
      //
      his = kalmanAlign.fCalibAlign->GetHisto(AliTPCcalibAlign::kPhi,is0,is1);
      if (!his) continue;
      delta[2] = his->GetMean();
      rms[2]=his->GetRMS();
      stat[2]=his->GetEntries();
      kalmanAlign.UpdateAlign1D(his->GetMean(),his->GetRMS(),is0,is1, vec[2],cov[2]);
      //
      his = kalmanAlign.fCalibAlign->GetHisto(AliTPCcalibAlign::kTheta,is0,is1);
      if (!his) continue;
      delta[3] = his->GetMean();       
      rms[3]=his->GetRMS();
      stat[3]=his->GetEntries();
      kalmanAlign.UpdateAlign1D(his->GetMean(),his->GetRMS(),is0,is1, vec[3],cov[3]);
      if (is1==is0+36) statUpDown+=Int_t(stat[0]);
      if (is1==is0+35) statLeftRight+=Int_t(stat[0]);
    }  

  fDelta1D[0] = (TMatrixD*)vec[0].Clone();
  fDelta1D[1] = (TMatrixD*)vec[1].Clone();
  fDelta1D[2] = (TMatrixD*)vec[2].Clone();
  fDelta1D[3] = (TMatrixD*)vec[3].Clone();
  //
  fCovar1D[0] = (TMatrixD*)cov[0].Clone();
  fCovar1D[1] = (TMatrixD*)cov[1].Clone();
  fCovar1D[2] = (TMatrixD*)cov[2].Clone();
  fCovar1D[3] = (TMatrixD*)cov[3].Clone();
  statUpDown/=36;
  statLeftRight/=36;
  MakeNewAlignment(kTRUE);
  //FitCE();
  for (Int_t is0=0;is0<72;is0++)
    for (Int_t is1=0;is1<72;is1++){
      Bool_t isPair=kFALSE;
      if (TMath::Abs(is0%18-is1%18)<2) isPair=kTRUE;
      if (TMath::Abs(is0%18-is1%18)==17) isPair=kTRUE;
      if (!isPair) continue;
      stat[0]=0; stat[1]=0; stat[2]=0; stat[3]=0; 
      //
      //
      his = kalmanAlign.fCalibAlign->GetHisto(AliTPCcalibAlign::kY,is0,is1);
      if (his){
	delta[0]=his->GetMean();
	rms[0]=his->GetRMS();
	stat[0]=his->GetEntries();
      }
      //     
      his = kalmanAlign.fCalibAlign->GetHisto(AliTPCcalibAlign::kZ,is0,is1);
      if (his) {
	delta[1]=his->GetMean();
	rms[1]=his->GetRMS();
	stat[1]=his->GetEntries();
      }
      //
      his = kalmanAlign.fCalibAlign->GetHisto(AliTPCcalibAlign::kPhi,is0,is1);
      if (his){
	delta[2] = his->GetMean();
	rms[2]=his->GetRMS();
	stat[2]=his->GetEntries();
      }
      //
      his = kalmanAlign.fCalibAlign->GetHisto(AliTPCcalibAlign::kTheta,is0,is1);
      if (his){
	delta[3] = his->GetMean();       
	rms[3]=his->GetRMS();
	stat[3]=his->GetEntries();
      }
      Int_t run = AliCDBManager::Instance()->GetRun();
      Float_t bz = AliTracker::GetBz();
      TVectorD fceG[6],fceL[6];
      for (Int_t ipar=0; ipar<6;ipar++){
	fceG[ipar]=*((TVectorD*)fFitCEGlobal->At(ipar));
	fceL[ipar]=*((TVectorD*)fFitCELocal->At(ipar));
      }
      (*pcstream)<<"kalmanAlignDebug"<<
	"run="<<run<<
	"bz="<<bz<<
	"is0="<<is0<<
	"is1="<<is1<<
	"delta.="<<&delta<<
	"rms.="<<&rms<<
	"stat.="<<&stat<<
	"vec0.="<<&vec[0]<<
	"vec1.="<<&vec[1]<<
	"vec2.="<<&vec[2]<<
	"vec3.="<<&vec[3]<<
	"pceIn0.="<<&paramCE[is0%36]<<
	"pceOut0.="<<&paramCE[is0%36+36]<<
	"pceIn1.="<<&paramCE[is1%36]<<
	"pceOut1.="<<&paramCE[is1%36+36]<<
	"fceG0.="<<&fceG[0]<<  // global fit of CE
	"fceG1.="<<&fceG[1]<<  // global fit of CE
	"fceG2.="<<&fceG[2]<<  // global fit of CE
	"fceG3.="<<&fceG[3]<<  // global fit of CE
	"fceG4.="<<&fceG[4]<<  // global fit of CE
	"fceL5.="<<&fceG[5]<<  // global fit of CE
	"fceL0.="<<&fceL[0]<<  // global fit of CE
	"fceL1.="<<&fceL[1]<<  // global fit of CE
	"fceL2.="<<&fceL[2]<<  // global fit of CE
	"fceL3.="<<&fceL[3]<<  // global fit of CE
	"fceL4.="<<&fceL[4]<<  // global fit of CE
	"fceL5.="<<&fceL[5]<<  // global fit of CE
	"\n";
    }
  
  Int_t run = AliCDBManager::Instance()->GetRun();
  Float_t bz = AliTracker::GetBz();
  (*pcstream)<<"runSummary"<<
    "run="<<run<<                      // run number 
    "bz="<<bz<<                        // bz field
    "statUpDown="<<statUpDown<<        // stat up-down
    "statLeftRight="<<statLeftRight<<  // stat left-right
    "\n";

  delete pcstream;
}






void AliTPCkalmanAlign::UpdateOCDBTime0( AliTPCCalPad  *pad, Int_t ustartRun, Int_t uendRun,  const char* storagePath ){
  //
  // Update OCDB
  // .x ConfigCalibTrain.C(117116)
  // AliTPCcalibDB::Instance()->GetPulserTmean()
  // pad->Add(-pad->GetMean())
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Marian Ivanov");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-25-01"); //root version
  metaData->SetComment("Calibration of the z - time misalignment");
  AliCDBId* id1=NULL;
  id1=new AliCDBId("TPC/Calib/PadTime0", ustartRun, uendRun);
  AliCDBStorage* gStorage = AliCDBManager::Instance()->GetStorage(storagePath);
  gStorage->Put(pad, (*id1), metaData);
}



void AliTPCkalmanAlign::DrawDeltaAlign(){
  //
  // Draw the reuslts of the alignment
  // Residual misalignment in respect with previous alignment shown
  //
  //
  TFile f("AliTPCkalmanAlign.root","update");
  TTree * treeDelta = (TTree*)f.Get("kalmanAlignDebug");
  TH1::AddDirectory(0);
  // 
  treeDelta->SetAlias("sector","is0");
  treeDelta->SetAlias("dYmeas","10*(delta.fElements[0])");
  treeDelta->SetAlias("dZmeas","10*(delta.fElements[1])");
  treeDelta->SetAlias("dPhimeas","1000*(delta.fElements[2])");
  treeDelta->SetAlias("dThetameas","1000*(delta.fElements[3])");
  //
  treeDelta->SetAlias("dYfit","10*(vec0.fElements[is0]-vec0.fElements[is1])");
  treeDelta->SetAlias("dZfit","10*(vec1.fElements[is0]-vec1.fElements[is1])");
  treeDelta->SetAlias("dPhifit","1000*(vec2.fElements[is0]-vec2.fElements[is1])");
  treeDelta->SetAlias("dThetafit","1000*(vec3.fElements[is0]-vec3.fElements[is1])");
  //
  treeDelta->SetMarkerStyle(25);
  treeDelta->SetMarkerColor(4);
  treeDelta->SetLineColor(4);
  const char *type[3]={"up-down","left-right","right-left"};
  const char *gropt[3]={"alp","lp","lp"};
  const char *varTypeY[3]={"dYmeas-dYfit:sector","dYfit:sector","10*vec0.fElements[is0]:sector"};
  const char *varLegendY[3]={"#Delta_{y} fit residual","#Delta_{y} fit value difference","#Delta_{y} sector"};
  const char *varTypeZ[3]={"dZmeas-dZfit:sector","dZfit:sector","10*vec1.fElements[is0]:sector"};
  const char *varLegendZ[3]={"#Delta_{Z} fit residual","#Delta_{Z} fit value difference","#Delta_{Z} sector"};
  const char *varTypeT[3]={"dThetameas-dThetafit:sector","dThetafit:sector","1000*vec3.fElements[is0]:sector"};
  const char *varLegendT[3]={"#Delta_{#theta} fit residual","#Delta_{#theta} fit value difference","#Delta_{#theta} sector"};
  const char *varTypeP[3]={"dPhimeas-dPhifit:sector","dPhifit:sector","1000*vec2.fElements[is0]:sector"};
  const char *varLegendP[3]={"#Delta_{#phi} fit residual","#Delta_{#phi} fit value difference","#Delta_{#phi} sector"};
  TLegend *legend = 0;
  Int_t grcol[3]={2,1,4};
  Int_t entries=0;
  TGraph *grDelta[3]={0,0,0};
  TCanvas * canvasDy=new TCanvas("canvasDy","canvasDy",1200,800);
  TCanvas * canvasDz=new TCanvas("canvasDz","canvasDz",1200,800);
  TCanvas * canvasDphi=new TCanvas("canvasDphi","canvasDphi",1200,800);
  TCanvas * canvasDtheta=new TCanvas("canvasDtheta","canvasDtheta",1200,800);
  canvasDy->Divide(2,2);
  canvasDz->Divide(2,2);
  canvasDtheta->Divide(2,2);
  canvasDphi->Divide(2,2);

  //
  // Dy
  //
  canvasDy->cd(1);
  treeDelta->Draw("dYmeas:dYfit");
  for (Int_t itype=0; itype<3; itype++){
    canvasDy->cd(itype+2);
    entries=treeDelta->Draw(varTypeY[itype],"is1==is0+36","goff");
    grDelta[0]=new TGraph(entries,treeDelta->GetV2(),treeDelta->GetV1());
    entries=treeDelta->Draw(varTypeY[itype],"is1==is0+35","goff");
    grDelta[1]=new TGraph(entries,treeDelta->GetV2(),treeDelta->GetV1());
    entries=treeDelta->Draw(varTypeY[itype],"is1==is0+37","goff");
    grDelta[2]=new TGraph(entries,treeDelta->GetV2(),treeDelta->GetV1());
    legend = new TLegend(0.5,0.7,0.9,0.9, varLegendY[itype]);
    for (Int_t i=0; i<3; i++) {
      grDelta[i]->SetMaximum(1.5);
      grDelta[i]->SetMinimum(-1);
      grDelta[i]->SetTitle(type[i]);
      grDelta[i]->SetMarkerColor(grcol[i]); 
      grDelta[i]->SetLineColor(grcol[i]); 
      grDelta[i]->SetMarkerStyle(25+i); 
      grDelta[i]->GetXaxis()->SetTitle("sector"); 
      grDelta[i]->GetYaxis()->SetTitle("#Delta_{y} (mm)"); 
      if (itype==2 && i>0) continue;
      grDelta[i]->Draw(gropt[i]); 
      legend->AddEntry(grDelta[i]);
    }
    legend->Draw();
  }
  //
  // Dz
  //
  canvasDz->cd();
  canvasDz->cd(1);
  treeDelta->Draw("dZmeas:dZfit");
  for (Int_t itype=0; itype<3; itype++){
    canvasDz->cd(itype+2);
    entries=treeDelta->Draw(varTypeZ[itype],"is1==is0+36","goff");
    grDelta[0]=new TGraph(entries,treeDelta->GetV2(),treeDelta->GetV1());
    entries=treeDelta->Draw(varTypeZ[itype],"is1==is0+35","goff");
    grDelta[1]=new TGraph(entries,treeDelta->GetV2(),treeDelta->GetV1());
    entries=treeDelta->Draw(varTypeZ[itype],"is1==is0+37","goff");
    grDelta[2]=new TGraph(entries,treeDelta->GetV2(),treeDelta->GetV1());
    legend = new TLegend(0.5,0.7,0.9,0.9, varLegendZ[itype]);
    for (Int_t i=0; i<3; i++) {
      grDelta[i]->SetMaximum(1.5);
      grDelta[i]->SetMinimum(-1);
      grDelta[i]->SetTitle(type[i]);
      grDelta[i]->SetMarkerColor(grcol[i]); 
      grDelta[i]->SetLineColor(grcol[i]); 
      grDelta[i]->SetMarkerStyle(25+i); 
      grDelta[i]->GetXaxis()->SetTitle("sector"); 
      grDelta[i]->GetYaxis()->SetTitle("#Delta_{z} (mm)"); 
      if (itype==2 && i>0) continue;
      grDelta[i]->Draw(gropt[i]); 
      legend->AddEntry(grDelta[i]);
    }
    legend->Draw();
  }

  //
  // Dtheta
  //
  canvasDtheta->cd(1);
  treeDelta->Draw("dThetameas:dThetafit");
  for (Int_t itype=0; itype<3; itype++){
    canvasDtheta->cd(itype+2);
    entries=treeDelta->Draw(varTypeT[itype],"is1==is0+36","goff");
    grDelta[0]=new TGraph(entries,treeDelta->GetV2(),treeDelta->GetV1());
    entries=treeDelta->Draw(varTypeT[itype],"is1==is0+35","goff");
    grDelta[1]=new TGraph(entries,treeDelta->GetV2(),treeDelta->GetV1());
    entries=treeDelta->Draw(varTypeT[itype],"is1==is0+37","goff");
    grDelta[2]=new TGraph(entries,treeDelta->GetV2(),treeDelta->GetV1());
    legend = new TLegend(0.5,0.7,0.9,0.9, varLegendT[itype]);
    for (Int_t i=0; i<3; i++) {
      grDelta[i]->SetMaximum(4.);
      grDelta[i]->SetMinimum(-3.);
      grDelta[i]->SetTitle(type[i]);
      grDelta[i]->SetMarkerColor(grcol[i]); 
      grDelta[i]->SetLineColor(grcol[i]); 
      grDelta[i]->SetMarkerStyle(25+i); 
      grDelta[i]->GetXaxis()->SetTitle("sector"); 
      grDelta[i]->GetYaxis()->SetTitle("#Delta_{#theta} (mrad)"); 
      if (itype==2 && i>0) continue;
      grDelta[i]->Draw(gropt[i]); 
      legend->AddEntry(grDelta[i]);
    }
    legend->Draw();
  }

  //
  // Dphi
  //
  canvasDphi->cd();
  canvasDphi->cd(1);
  treeDelta->Draw("dPhimeas:dPhifit");
  for (Int_t itype=0; itype<3; itype++){
    canvasDphi->cd(itype+2);
    entries=treeDelta->Draw(varTypeP[itype],"is1==is0+36","goff");
    grDelta[0]=new TGraph(entries,treeDelta->GetV2(),treeDelta->GetV1());
    entries=treeDelta->Draw(varTypeP[itype],"is1==is0+35","goff");
    grDelta[1]=new TGraph(entries,treeDelta->GetV2(),treeDelta->GetV1());
    entries=treeDelta->Draw(varTypeP[itype],"is1==is0+37","goff");
    grDelta[2]=new TGraph(entries,treeDelta->GetV2(),treeDelta->GetV1());
    legend = new TLegend(0.5,0.7,0.9,0.9, varLegendP[itype]);
    for (Int_t i=0; i<3; i++) {
      grDelta[i]->SetMaximum(4.);
      grDelta[i]->SetMinimum(-3.);
      grDelta[i]->SetTitle(type[i]);
      grDelta[i]->SetMarkerColor(grcol[i]); 
      grDelta[i]->SetLineColor(grcol[i]); 
      grDelta[i]->SetMarkerStyle(25+i); 
      grDelta[i]->GetXaxis()->SetTitle("sector"); 
      grDelta[i]->GetYaxis()->SetTitle("#Delta_{#phi} (mrad)"); 
      if (itype==2 && i>0) continue;
      grDelta[i]->Draw(gropt[i]); 
      legend->AddEntry(grDelta[i]);
    }
    legend->Draw();
  }
  //
  //
  f.cd();
  canvasDy->Write();
  canvasDz->Write();
  canvasDtheta->Write();
  canvasDphi->Write();
  //  
  //
  //
  TPostScript *ps = new TPostScript("alignReport.ps", 112); 
  ps->NewPage();
  canvasDy->Draw();
  ps->NewPage();
  canvasDz->Draw();
  ps->NewPage();
  canvasDtheta->Draw();
  ps->NewPage();
  canvasDphi->Draw();
  ps->Close();
  delete ps;
}



void AliTPCkalmanAlign::DumpOldAlignment(TTreeSRedirector *pcstream){
  // Dump the content of old alignemnt
  // Expected that the old alignmnet is loaded
  //
  if (!fOriginalAlign) return;
  //
  TVectorD localTrans(3);
  TVectorD globalTrans(3);
  TVectorD localRot(3);
  TVectorD globalRot(3);
  AliGeomManager::ELayerID idLayer;
  Int_t idModule=0;
  //
  for (Int_t i=0; i<fOriginalAlign->GetEntries();i++){
    AliAlignObjParams *params = (AliAlignObjParams*)fOriginalAlign->At(i);
    params->GetVolUID(idLayer,idModule);
    params->GetLocalTranslation(localTrans.GetMatrixArray());
    params->GetLocalAngles(localRot.GetMatrixArray());
    params->GetTranslation(globalTrans.GetMatrixArray());
    params->GetAngles(globalRot.GetMatrixArray());
    Int_t sector=idModule;
    if (idLayer>7) sector+=36;
    (*pcstream)<<"oldAlign"<<
      //"idLayer="<<idLayer<<
      "idModule="<<idModule<<
      "sector="<<sector<<
      "lT.="<<&localTrans<<
      "gT.="<<&localTrans<<
      "lR.="<<&localRot<<
      "gR.="<<&globalRot<<
      "\n";
  }
}


void AliTPCkalmanAlign::MakeNewAlignment(Bool_t badd, TTreeSRedirector * pcstream){
  //
  // make a new Alignment entry
  //
  if (!fOriginalAlign) return;
  //
  TVectorD localTrans(3);
  TVectorD globalTrans(3);
  TVectorD localRot(3);
  TVectorD globalRot(3);
  //
  TVectorD localTransNew(3);   // new entries
  TVectorD globalTransNew(3);
  TVectorD localRotNew(3);
  TVectorD globalRotNew(3);
  //
  AliGeomManager::ELayerID idLayer;
  Int_t idModule=0;
  //
  fNewAlign = (TClonesArray*)fOriginalAlign->Clone();
  for (Int_t i=0; i<fOriginalAlign->GetEntries();i++){
    AliAlignObjParams *params = (AliAlignObjParams*)fOriginalAlign->At(i);
    //AliAlignObjParams *paramsNew = (AliAlignObjParams*)fNewAlign->At(i);
    params->GetVolUID(idLayer,idModule);
    Int_t sector=(Int_t)idModule;
    if (idLayer>7) sector+=36;
    params->GetLocalTranslation(localTrans.GetMatrixArray());
    params->GetLocalAngles(localRot.GetMatrixArray());
    params->GetTranslation(globalTrans.GetMatrixArray());
    params->GetAngles(globalRot.GetMatrixArray());
    //
    //
    //
    if (badd){ // addition if 
      localTransNew=localTrans;
      localRotNew=localRot;
    }
    localTransNew[1]=localTransNew[1]-((*fDelta1D[0])(sector,0));
    localRot[0]  =localRot[0]-(*fDelta1D[2])(sector,0);
    //
    if (pcstream) (*pcstream)<<"alignParams"<<
      //"idLayer="<<idLayer<<
      "idModule="<<idModule<<
      "sector="<<sector<<
      "olT.="<<&localTrans<<
      "olR.="<<&localRot<<
      "ogT.="<<&localTrans<<
      "ogR.="<<&globalRot<<
      "nlT.="<<&localTransNew<<
      "nlR.="<<&localRotNew<<
      "ngT.="<<&localTransNew<<
      "ngR.="<<&globalRotNew<<
      "\n";
  }
}



void AliTPCkalmanAlign::DrawAlignmentTrends(){
  //
  // Draw trends of alingment variables
  //
  /*
  */
  AliXRDPROOFtoolkit toolkit;
  TChain * chain = toolkit.MakeChainRandom("align.list.Good","kalmanAlignDebug",0,2000);
  TChain * chainRef = toolkit.MakeChainRandom("alignRef.list","kalmanAlignDebug",0,2000);
  chain->AddFriend(chainRef,"R");
  chainRef->AddFriend(chainRef,"T");
  //cuts
  TCut cutS="stat.fElements[0]>200&&stat.fElements[1]>200&&stat.fElements[3]>200&&stat.fElements[3]>200";   //statistic in the bin
  TCut cutST="T.stat.fElements[0]>200&&T.stat.fElements[1]>200&&T.stat.fElements[3]>200&&T.stat.fElements[3]>200";   //statistic in the bin
  //  TTree *tree  = chain->CopyTree(cutS);
  //TTree *treeR = chainRef->CopyTree(cutST);

  TCanvas * canvasDy= new TCanvas("canvasDy","canvasDy");
  TH1 *his=0;
  TLegend *legend = 0;
  //  Int_t grcol[3]={2,1,4};
  
  legend = new TLegend(0.7,0.6,0.9,0.9, "Alignment #Delta_{y}- Up-Down");
  for (Int_t isec=0; isec<18; isec+=2){
    chain->SetMarkerColor(1+(isec%5));
    chain->SetMarkerStyle(isec+20);
    chain->Draw("10*(delta.fElements[0]-R.delta.fElements[0]):run",cutS+Form("is1==is0+36&&is0==%d",isec),"profgoff");
    his = (TH1*)(chain->GetHistogram()->Clone());
    his->SetName(Form("#Delta_{Y} sector %d",isec));
    his->SetTitle(Form("#Delta_{Y} sector %d",isec));
    his->SetMaximum(1.);
    his->SetMinimum(-1.);
    his->GetYaxis()->SetTitle("#Delta_{y} (mm)");
    his->GetXaxis()->SetTitle("run Number");
    if (isec==0) his->Draw("");
    if (isec>0) his->Draw("same");
    legend->AddEntry(his);
  }
  legend->Draw();
  canvasDy->Draw();
}






void AliTPCkalmanAlign::FitCE(){
  //
  // fit CE 
  // 1. Global fit - gy and gx
  // 2. Local X fit common 
  // 3. Sector fit 
  //
  AliTPCPreprocessorOnline * preprocesor = new AliTPCPreprocessorOnline;
  //
  AliTPCCalPad *padTime0 = AliTPCcalibDB::Instance()->GetPadTime0();
  AliTPCCalPad *padNoise = AliTPCcalibDB::Instance()->GetPadNoise();
  AliTPCCalPad * ceTmean = AliTPCcalibDB::Instance()->GetCETmean();  // CE information
  AliTPCCalPad * ceTrms  = AliTPCcalibDB::Instance()->GetCETrms();
  AliTPCCalPad * ceQmean = AliTPCcalibDB::Instance()->GetCEQmean();  
  AliTPCCalPad * pulserTmean = AliTPCcalibDB::Instance()->GetPulserTmean(); //
  AliTPCCalPad * pulserTrms = AliTPCcalibDB::Instance()->GetPulserTrms();
  AliTPCCalPad * pulserQmean = AliTPCcalibDB::Instance()->GetPulserQmean();
  AliTPCCalPad * dmap0  = AliTPCcalibDB::Instance()->GetDistortionMap(0);   // distortion maps
  AliTPCCalPad * dmap1  = AliTPCcalibDB::Instance()->GetDistortionMap(1);
  AliTPCCalPad * dmap2  = AliTPCcalibDB::Instance()->GetDistortionMap(2);
  pulserTmean->Add(-pulserTmean->GetMean());
  //
  preprocesor->AddComponent(padTime0->Clone());
  preprocesor->AddComponent(padNoise->Clone());
  preprocesor->AddComponent(pulserTmean->Clone());
  preprocesor->AddComponent(pulserQmean->Clone());
  preprocesor->AddComponent(pulserTrms->Clone());
  preprocesor->AddComponent(ceTmean->Clone());
  preprocesor->AddComponent(ceQmean->Clone());
  preprocesor->AddComponent(ceTrms->Clone());
  preprocesor->AddComponent(dmap0->Clone());
  preprocesor->AddComponent(dmap1->Clone());
  preprocesor->AddComponent(dmap2->Clone());
  preprocesor->DumpToFile("cetmean.root");

  TCut cutNoise="abs(PadNoise.fElements/PadNoise_Median-1)<0.3";
  TCut cutPulserT="abs(PulserTrms.fElements/PulserTrms_Median-1)<0.2";
  TCut cutPulserQ="abs(PulserQmean.fElements/PulserQmean_Median-1)<0.2";
  TCut cutCEQ="CEQmean.fElements>50";
  TCut cutCET="abs(CETmean.fElements)<2";
  TCut cutAll=cutNoise+cutPulserT+cutPulserQ+cutCEQ+cutCET;
  //
  //
  TFile * f = new TFile("cetmean.root");
  TTree * chain = (TTree*) f->Get("calPads");
  Int_t entries = chain->Draw("1",cutAll,"goff");
  if (entries<200000) return;  // no calibration available - pulser or CE or noise

  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD param;
  TMatrixD covar;
  //
  // make a aliases
  AliTPCkalmanAlign::MakeAliasCE(chain);
  TString  fstringG="";              // global part
  //
  fstringG+="Gy++";                  // par 1 - global y
  fstringG+="Gx++";                  // par 2 - global x
  // 
  fstringG+="isin++";                // delta IROC-OROC offset
  fstringG+="Lx++";                  // common slope 
  fstringG+="Lx*isin++";             // delta slope 
  fstringG+="Ly++";                  // common slope 
  fstringG+="Ly*isin++";             // delta slope 
  TVectorD vecG[2];
  TString * strFitG=0;
  TString * strFitLX=0;
  //
  strFitG = TStatToolkit::FitPlane(chain,"deltaT", fstringG.Data(),"sideA"+cutAll, chi2,npoints,vecG[0],covar,-1,0, 10000000, kFALSE);
  chain->SetAlias("tfitGA",strFitG->Data());
  strFitG->Tokenize("++")->Print();
  printf("chi2=%f\n",TMath::Sqrt(chi2/npoints));
  //
  strFitG = TStatToolkit::FitPlane(chain,"deltaT", fstringG.Data(),"sideC"+cutAll, chi2,npoints,vecG[1],covar,-1,0, 10000000, kFALSE);
  chain->SetAlias("tfitGC",strFitG->Data());
  strFitG->Tokenize("++")->Print();
  printf("chi2=%f\n",TMath::Sqrt(chi2/npoints));
  //
  AliTPCCalPad *padFitG =AliTPCCalPad::CreateCalPadFit("1++gy/500.++gx/500.++0+++0++0++0++0",vecG[0],vecG[1]);
  AliTPCCalPad *padFitLX=AliTPCCalPad::CreateCalPadFit("0++0++0++(sector<36)++(lx-133)/100++(sector<36)*(lx-133)/100.++(ly)/100++(sector<36)*(ly)/100.",vecG[0],vecG[1]);
  // swap a side and c side
  AliTPCCalPad *padFitGSwap =AliTPCCalPad::CreateCalPadFit("1++gy/500.++gx/500.++0+++0++0++0++0",vecG[1],vecG[0]);
  AliTPCCalPad *padFitLXSwap=AliTPCCalPad::CreateCalPadFit("0++0++0++(sector<36)++(lx-133)/100++(sector<36)*(lx-133)/100.++(ly)/100++(sector<36)*(ly)/100.",vecG[1],vecG[0]);
  padFitG->SetName("CEG");
  padFitLX->SetName("CELX");
  padFitGSwap->SetName("CEGS");
  padFitLXSwap->SetName("CELXS");
  preprocesor->AddComponent(padFitG->Clone());
  preprocesor->AddComponent(padFitLX->Clone());
  preprocesor->AddComponent(padFitGSwap->Clone());
  preprocesor->AddComponent(padFitLXSwap->Clone());
  preprocesor->DumpToFile("cetmean.root");   // add it to the file
  //
  // make local fits
  //
  f = new TFile("cetmean.root");
  chain = (TTree*) f->Get("calPads");
  AliTPCkalmanAlign::MakeAliasCE(chain);
  TString  fstringL="";              // local fit
  //                                 // 0. delta common
  fstringL+="isin++";                // 1. delta IROC-OROC offset
  fstringL+="Lx++";                  // 2. common slope 
  fstringL+="Lx*isin++";             // 3. delta slope 
  fstringL+="Ly++";                  // 2. common slope 
  fstringL+="Ly*isin++";             // 3. delta slope 
  TVectorD vecL[36];
  TVectorD dummy(6);
  AliTPCCalPad *padFitLCE = new AliTPCCalPad("LocalCE","LocalCE");
  AliTPCCalPad *padFitTmpCE;
  for (Int_t isec=0; isec<36; isec++){
    TCut cutSector=Form("(sector%36)==%d",isec);
    strFitLX = TStatToolkit::FitPlane(chain,"deltaT-CEG.fElements-CELX.fElements", fstringL.Data(),cutSector+cutAll+"abs(deltaT-CEG.fElements-CELX.fElements)<0.4", chi2,npoints,vecL[isec],covar,-1,0, 10000000, kFALSE);
    printf("sec=%d\tchi2=%f\n",isec,TMath::Sqrt(chi2/npoints));
    //
    TString fitL=Form("((sector%36)==%d)++((sector%36)==%d)*(sector<36)++((sector%36)==%d)*(lx-133)/100.++((sector%36)==%d)*(sector<36)*(lx-133)/100.++((sector%36)==%d)*(ly)/100.++((sector%36)==%d)*(sector<36)*(ly)/100.",isec,isec,isec,isec);
    if (isec<18) padFitTmpCE=AliTPCCalPad::CreateCalPadFit(fitL.Data(),vecL[isec],dummy);
    if (isec>=18) padFitTmpCE=AliTPCCalPad::CreateCalPadFit(fitL.Data(),dummy,vecL[isec]);
    padFitLCE->Add(padFitTmpCE);
  }
  //
  padFitLCE->SetName("CELocal");
  preprocesor->AddComponent(padFitLCE->Clone());
  preprocesor->DumpToFile("cetmean.root");   // add it to the file
  //
  // write data to array
  //
  fFitCEGlobal = new TObjArray(6); 
  fFitCELocal  = new TObjArray(6); 
  for (Int_t ipar=0; ipar<6;ipar++){
    fFitCEGlobal->AddAt(new TVectorD(36),ipar);
    fFitCELocal->AddAt(new TVectorD(36),ipar);
    //
    TVectorD &fvecG = *((TVectorD*)fFitCEGlobal->At(ipar));
    TVectorD &fvecL = *((TVectorD*)fFitCELocal->At(ipar));
    //
    for (Int_t isec=0; isec<36;isec++){      
      fvecL[isec]=vecL[isec][ipar];
      if (ipar>0){
	if (isec<18)  fvecG[isec]=vecG[0][ipar+2];
	if (isec>=18) fvecG[isec]=vecG[1][ipar+2];
      }
    }
  }
  //
  // 
  //
}

void AliTPCkalmanAlign::MakeAliasCE(TTree * chain){
  //
  // make a aliases of pad variables
  //
  chain->SetAlias("side","(-1+(sector%36<18)*2)");
  chain->SetAlias("sideA","(sector%36<18)");
  chain->SetAlias("sideC","(sector%36>=18)");
  chain->SetAlias("isin","(sector<36)");
  chain->SetAlias("deltaT","CETmean.fElements-PulserTmean.fElements");
  chain->SetAlias("timeP","PulserTmean.fElements");
  chain->SetAlias("Gy","(gy.fElements/500.)");
  chain->SetAlias("Gx","(gx.fElements/500.)");
  chain->SetAlias("Lx","(lx.fElements-133)/100.");   // lx in meters
  chain->SetAlias("Ly","(ly.fElements)/100.");
  chain->SetAlias("La","(ly.fElements/lx.fElements/0.155)");
  chain->SetAlias("deltaT","(CETmean.fElements-PulserTmean.fElements)");
}
