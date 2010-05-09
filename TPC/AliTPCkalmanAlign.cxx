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
  Int_t run=117092;
  .x ConfigCalibTrain.C(run)

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


#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliAlignObjParams.h"
#include "AliTPCROC.h"
#include "TFile.h"
#include "TLinearFitter.h"
#include "AliTPCcalibAlign.h"
#include "TH1.h"
#include "AliTPCCalPad.h"
#include "AliTPCkalmanAlign.h"
#include "TLegend.h"
#include "TCanvas.h"

AliTPCkalmanAlign::AliTPCkalmanAlign():
  TNamed(),
  fCalibAlign(0),     // kalman alignnmnt
  fOriginalAlign(0),   // original alignment 0 read for the OCDB
  fNewAlign(0)
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
  fNewAlign(0)     // kalman alignnmnt
{
  //
  // Default constructor
  //
  for (Int_t i=0; i<4; i++){
    fDelta1D[i]=0;
    fCovar1D[i]=0;
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
  TH1 * his=0;
  for (Int_t is0=0;is0<72;is0++)
    for (Int_t is1=0;is1<72;is1++){
      //
      //
      his = kalmanAlign.fCalibAlign->GetHisto(AliTPCcalibAlign::kY,is0,is1);
      if (!his) continue;
      if (his->GetEntries()<kMinEntries) continue;
      delta[0]=his->GetMean();
      kalmanAlign.UpdateAlign1D(his->GetMean(),his->GetRMS(),is0,is1, vec[0],cov[0]);
      //     
      his = kalmanAlign.fCalibAlign->GetHisto(AliTPCcalibAlign::kZ,is0,is1);
      if (!his) continue;
      delta[1]=his->GetMean();
      kalmanAlign.UpdateAlign1D(his->GetMean(),his->GetRMS(),is0,is1, vec[1],cov[1]);
      //
      his = kalmanAlign.fCalibAlign->GetHisto(AliTPCcalibAlign::kPhi,is0,is1);
      if (!his) continue;
      delta[2] = his->GetMean();
      kalmanAlign.UpdateAlign1D(his->GetMean(),his->GetRMS(),is0,is1, vec[2],cov[2]);
      //
      his = kalmanAlign.fCalibAlign->GetHisto(AliTPCcalibAlign::kTheta,is0,is1);
      if (!his) continue;
      delta[3] = his->GetMean();       
      kalmanAlign.UpdateAlign1D(his->GetMean(),his->GetRMS(),is0,is1, vec[3],cov[3]);
    }  
  for (Int_t is0=0;is0<72;is0++)
    for (Int_t is1=0;is1<72;is1++){
      //
      //
      his = kalmanAlign.fCalibAlign->GetHisto(AliTPCcalibAlign::kY,is0,is1);
      if (!his) continue;
      if (his->GetEntries()<kMinEntries) continue;
      delta[0]=his->GetMean();
      kalmanAlign.UpdateAlign1D(his->GetMean(),his->GetRMS(),is0,is1, vec[0],cov[0]);
      //     
      his = kalmanAlign.fCalibAlign->GetHisto(AliTPCcalibAlign::kZ,is0,is1);
      if (!his) continue;
      delta[1]=his->GetMean();
      kalmanAlign.UpdateAlign1D(his->GetMean(),his->GetRMS(),is0,is1, vec[1],cov[1]);
      //
      his = kalmanAlign.fCalibAlign->GetHisto(AliTPCcalibAlign::kPhi,is0,is1);
      if (!his) continue;
      delta[2] = his->GetMean();
      kalmanAlign.UpdateAlign1D(his->GetMean(),his->GetRMS(),is0,is1, vec[2],cov[2]);
      //
      his = kalmanAlign.fCalibAlign->GetHisto(AliTPCcalibAlign::kTheta,is0,is1);
      if (!his) continue;
      delta[3] = his->GetMean();       
      kalmanAlign.UpdateAlign1D(his->GetMean(),his->GetRMS(),is0,is1, vec[3],cov[3]);
      (*pcstream)<<"kalmanAlignDebug"<<
	"is0="<<is0<<
	"is1="<<is1<<
	"delta.="<<&delta<<
	"vec0.="<<&vec[0]<<
	"vec1.="<<&vec[1]<<
	"vec2.="<<&vec[2]<<
	"vec3.="<<&vec[3]<<
	"\n";
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
  canvasDy->SaveAs("alignDy.pdf");
  canvasDz->SaveAs("alignDz.pdf");
  canvasDtheta->SaveAs("alignDtheta.pdf");
  canvasDphi->SaveAs("alignDphi.pdf");
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
    AliAlignObjParams *paramsNew = (AliAlignObjParams*)fNewAlign->At(i);
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
//     localTrans[1]=localTrans[1]-(*fDelta1D[0])[sector];
//     localRot[0]  =localRot[0]-(*fDelta1D[0])[sector];
    //
    if (pcstream) (*pcstream)<<"newAlign"<<
      //"idLayer="<<idLayer<<
      "idModule="<<idModule<<
      "sector="<<sector<<
      "olT.="<<&localTrans<<
      "ogT.="<<&localTrans<<
      "olR.="<<&localRot<<
      "ogR.="<<&globalRot<<
      "\n";
  }
  

}
