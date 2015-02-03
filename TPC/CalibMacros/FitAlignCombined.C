/// \file FitAlignCombined.C
///
/// \author marian.ivanov@cern.ch
///
/// Macro to create  alignment/distortion maps
/// As a input output of AliTPCcalibAlign and AliTPCcalibTime is used.
/// distortion lookup tables are used.
/// 
/// Input file mean.root with distortion tree expected to be in directory:
/// ../mergeField0/mean.root
/// 
/// The ouput file fitAlignCombined.root contains:
/// 1. Resulting (residual) AliTPCCalibMisalignment 
/// 2. QA fit plots
/// 
/// Functions documented inside:
/// 
/// RegisterAlignFunction();
/// MakeAlignFunctionGlobal();
/// MakeAlignFunctionSector();
/// 
/// Example usage:
///
/// ~~~ 
/// .x ~/NimStyle.C
/// gROOT->Macro("~/rootlogon.C");
/// gSystem->Load("libANALYSIS");
/// gSystem->Load("libSTAT");
/// gSystem->Load("libTPCcalib");
/// gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros -I$ALICE_ROOT/TPC/TPC -I$ALICE_ROOT/STAT");
/// .L $ALICE_ROOT/TPC/CalibMacros/FitAlignCombined.C+
/// .x ConfigCalibTrain.C(119047)
///
/// FitAlignCombinedCorr();
/// ~~~

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TH1D.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TVectorD.h"
#include "TTreeStream.h"
#include "TFile.h"
#include "TChain.h"
#include "AliTPCcalibAlign.h"
#include "AliTPCROC.h"
#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TROOT.h"
#include "TProfile.h"
#include "AliTPCPreprocessorOnline.h"
#include "AliTPCcalibDB.h"
#include "AliTPCkalmanAlign.h"
#include "TPostScript.h"
#include "TLegend.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "AliLog.h"
#include "AliTPCExBEffectiveSector.h"
#include "AliTPCComposedCorrection.h"
#include "AliTPCCalibGlobalMisalignment.h"
//
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "TStatToolkit.h"
#include "AliAlignObjParams.h"
#include "AliTPCParam.h"
#endif


void DrawFitQA();

AliTPCROC *proc = AliTPCROC::Instance();
AliTPCCalibGlobalMisalignment *combAlignOCDBOld=0;  // old misalignment from OCDB - used for reconstruction
AliTPCCalibGlobalMisalignment *combAlignOCDBNew=0;  // new misalignment
AliTPCCalibGlobalMisalignment *combAlignGlobal=0;   // delta misalignment - global part
AliTPCCalibGlobalMisalignment *combAlignLocal=0;    // delta misalignment - sector/local part
//

TTreeSRedirector *pcstream= 0;    // Workspace to store the fit parameters and QA histograms
//
TChain * chain=0;
TChain * chainZ=0;                // trees with z measurement
TTree * tree=0;
TTree *itsdy=0;
TTree *itsdyP=0;
TTree *itsdyM=0;
TTree *itsdp=0;
TTree *itsdphiP=0;
TTree *itsdphiM=0;

TTree *itsdz=0;
TTree *itsdt=0;
TTree *tofdy=0;
TTree *trddy=0;
//
TTree *vdy=0;
TTree *vdyP=0;
TTree *vdyM=0;
TTree *vds=0;
TTree *vdz=0;
TTree *vdt=0;

TCut cutS="entries>1000&&abs(snp)<0.2&&abs(theta)<1.";


void RegisterAlignFunction(){
  /// Register primitive alignment components.
  /// Linear conbination of primitev forulas used for fit
  /// The nominal delta 1 mm in shift and 1 mrad in rotation
  /// Primitive formulas registeren in AliTPCCoreection::AddvisualCorrection
  /// 0 - deltaX
  /// 1 - deltaY
  /// 2 - deltaZ
  /// 3 - rot0 (phi)
  /// 4 - rot1 (theta)
  /// 5 - rot2

  TGeoHMatrix matrixX;
  TGeoHMatrix matrixY;
  TGeoHMatrix matrixZ;
  TGeoRotation rot0;
  TGeoRotation rot1;
  TGeoRotation rot2;  //transformation matrices
  TGeoRotation rot90;  //transformation matrices
  matrixX.SetDx(0.1); matrixY.SetDy(0.1); matrixZ.SetDz(0.1); //1 mm translation
  rot0.SetAngles(0.001*TMath::RadToDeg(),0,0);
  rot1.SetAngles(0,0.001*TMath::RadToDeg(),0);
  rot2.SetAngles(0,0,0.001*TMath::RadToDeg());
  //how to get rot02 ?
  rot90.SetAngles(0,90,0);
  rot2.MultiplyBy(&rot90,kTRUE);
  rot90.SetAngles(0,-90,0);
  rot2.MultiplyBy(&rot90,kFALSE);
  AliTPCCalibGlobalMisalignment *alignRot0  =new  AliTPCCalibGlobalMisalignment;
  alignRot0->SetAlignGlobal(&rot0);
  AliTPCCalibGlobalMisalignment *alignRot1=new  AliTPCCalibGlobalMisalignment;
  alignRot1->SetAlignGlobal(&rot1);
  AliTPCCalibGlobalMisalignment *alignRot2=new  AliTPCCalibGlobalMisalignment;
  alignRot2->SetAlignGlobal(&rot2);
  //
  AliTPCCalibGlobalMisalignment *alignTrans0  =new  AliTPCCalibGlobalMisalignment;
  alignTrans0->SetAlignGlobal(&matrixX);
  AliTPCCalibGlobalMisalignment *alignTrans1=new  AliTPCCalibGlobalMisalignment;
  alignTrans1->SetAlignGlobal(&matrixY);
  AliTPCCalibGlobalMisalignment *alignTrans2=new  AliTPCCalibGlobalMisalignment;
  alignTrans2->SetAlignGlobal(&matrixZ);
  AliTPCCorrection::AddVisualCorrection(alignTrans0  ,0);
  AliTPCCorrection::AddVisualCorrection(alignTrans1  ,1);
  AliTPCCorrection::AddVisualCorrection(alignTrans2  ,2);
  AliTPCCorrection::AddVisualCorrection(alignRot0    ,3);
  AliTPCCorrection::AddVisualCorrection(alignRot1    ,4);
  AliTPCCorrection::AddVisualCorrection(alignRot2    ,5);
  TObjArray arrAlign(6);
  arrAlign.AddAt(alignTrans0->Clone(),0);
  arrAlign.AddAt(alignTrans1->Clone(),1);
  arrAlign.AddAt(alignTrans2->Clone(),2);
  arrAlign.AddAt(alignRot0->Clone(),3);
  arrAlign.AddAt(alignRot1->Clone(),4);
  arrAlign.AddAt(alignRot2->Clone(),5);
  //combAlign.SetCorrections((TObjArray*)arrAlign.Clone());
}

AliTPCCalibGlobalMisalignment * MakeAlignFunctionGlobal(TVectorD paramYGlobal){
  /// Take a fit parameters and make a combined correction
  /// 1. Take the common part
  /// 3. Make combined AliTPCCalibGlobalMisalignment - register it
  /// 4. Compare the aliases with fit values - IT is OK

  AliTPCCalibGlobalMisalignment *alignGlobal  =new  AliTPCCalibGlobalMisalignment;
  TGeoHMatrix matGlobal; // global parameters
  TGeoHMatrix matDelta;  // delta A side - C side
  //
  TGeoHMatrix matrixX;
  TGeoHMatrix matrixY;
  TGeoHMatrix matrixZ;
  TGeoRotation rot0;
  TGeoRotation rot1;
  TGeoRotation rot2;  //transformation matrices
  TGeoRotation rot90;  //transformation matrices
  //
  // 
  matrixX.SetDx(0.1*paramYGlobal[1]); 
  matrixY.SetDy(0.1*paramYGlobal[2]); 
  rot0.SetAngles(0.001*TMath::RadToDeg()*paramYGlobal[3],0,0);
  rot1.SetAngles(0,0.001*TMath::RadToDeg()*paramYGlobal[4],0);
  rot2.SetAngles(0,0,0.001*TMath::RadToDeg()*paramYGlobal[5]);
  rot90.SetAngles(0,90,0);
  rot2.MultiplyBy(&rot90,kTRUE);
  rot90.SetAngles(0,-90,0);
  rot2.MultiplyBy(&rot90,kFALSE);
  matGlobal.Multiply(&matrixX);
  matGlobal.Multiply(&matrixY);
  matGlobal.Multiply(&rot0);
  matGlobal.Multiply(&rot1);
  matGlobal.Multiply(&rot2);
  alignGlobal->SetAlignGlobal((TGeoMatrix*)matGlobal.Clone());
  AliTPCCorrection::AddVisualCorrection(alignGlobal    ,100);  
  //
  /*
    chain->Draw("mean-deltaG:int(sector)",cutS+"type==0&&refX==0&&theta>0","profsame");
    chain->Draw("mean-deltaG:int(sector+18)",cutS+"type==0&&refX==0&&theta<0","profsame");

    chain->Draw("deltaG:fitYGlobal",cutS+"type==0");
    chain->Draw("deltaG:fitYGlobal",cutS+"type==2");

  */

  return alignGlobal;
}


AliTPCCalibGlobalMisalignment * MakeAlignFunctionSector(TVectorD paramYLocal){
  /// Take a fit parameters and make a combined correction:
  /// Only delta local Y and delta phi are fitted - not sensityvity for other parameters
  /// Algorithm:
  ///   1. Loop over sectors
  ///   2. Make combined AliTPCCalibGlobalMisalignment - register it
  ///   3. Compare the aliases with fit values - IT is OK

  AliTPCCalibGlobalMisalignment *alignLocal  =new  AliTPCCalibGlobalMisalignment;
  TGeoHMatrix matrixX;
  TGeoHMatrix matrixY;
  TGeoHMatrix matrixZ;
  TGeoRotation rot0;
  TGeoRotation rot1;
  TGeoRotation rot2;  //transformation matrices
  TGeoRotation rot90;  //transformation matrices
  TObjArray * array = new TObjArray(72);
  TVectorD vecLY(72);
  TVectorD vecGY(72);
  TVectorD vecGX(72);
  TVectorD vecPhi(72);
  TVectorD vecSec(72);
  //
  // 
  {for (Int_t isec=0; isec<72; isec++){
      Int_t sector=isec%18;
      Int_t offset=sector*4+1;
      if ((isec%36)>=18) offset+=2;
      Double_t phi= (Double_t(sector)+0.5)*TMath::Pi()/9.;
      //
      Double_t dly = paramYLocal[offset];
      Double_t dphi= paramYLocal[offset+1];
      Double_t dgx =   TMath::Cos(phi)*0+TMath::Sin(phi)*dly;
      Double_t dgy =   TMath::Sin(phi)*0-TMath::Cos(phi)*dly;
      vecSec[isec]=isec;
      vecLY[isec]=dly;
      vecGX[isec]=dgx;
      vecGY[isec]=dgy;
      vecPhi[isec]=dphi;
      //
      matrixX.SetDx(dgx); 
      matrixY.SetDy(dgy); 
      rot0.SetAngles(0.001*TMath::RadToDeg()*dphi,0,0);
      TGeoHMatrix matrixSector; // global parameters
      matrixSector.Multiply(&matrixX);
      matrixSector.Multiply(&matrixY);
      matrixSector.Multiply(&rot0);
      array->AddAt(matrixSector.Clone(),isec);
    }}
  alignLocal->SetAlignSectors(array);
  AliTPCCorrection::AddVisualCorrection(alignLocal,101);  
  //
  //
  TGraph * graphLY = new TGraph(72,vecSec.GetMatrixArray(), vecLY.GetMatrixArray());
  TGraph * graphGX = new TGraph(72,vecSec.GetMatrixArray(), vecGX.GetMatrixArray());
  TGraph * graphGY = new TGraph(72,vecSec.GetMatrixArray(), vecGY.GetMatrixArray());
  TGraph * graphPhi = new TGraph(72,vecSec.GetMatrixArray(), vecPhi.GetMatrixArray());
  graphLY->SetMarkerStyle(25); graphGX->SetMarkerStyle(25);  graphGY->SetMarkerStyle(25);   graphPhi->SetMarkerStyle(25);
  graphLY->SetName("DeltaLY"); graphLY->SetTitle("#Delta_{ly} (mm)");
  graphGX->SetName("DeltaGX"); graphGX->SetTitle("#Delta_{gx} (mm)");
  graphGY->SetName("DeltaGY"); graphGY->SetTitle("#Delta_{gy} (mm)");
  graphPhi->SetName("DeltaPhi"); graphPhi->SetTitle("#Delta_{phi} (mrad)");
  //
  graphLY->Write("grDeltaLY");
  graphGX->Write("grDeltaGX");
  graphGY->Write("grDeltaGY");
  graphPhi->Write("grDeltaPhi");
  // Check:
  /*
    graphLY.Draw("alp")
    chain->Draw("mean-deltaG:int(sector)",cutS+"type==0&&refX==0&&theta>0","profsame");
    chain->Draw("mean-deltaG:int(sector+18)",cutS+"type==0&&refX==0&&theta<0","profsame");

    graphPhi.Draw("alp")
    chain->Draw("1000*(mean-deltaG):int(sector)",cutS+"type==2&&refX==0&&theta>0","profsame");
    chain->Draw("1000*(mean-deltaG):int(sector+18)",cutS+"type==2&&refX==0&&theta<0","profsame");

    // 
    chain->Draw("deltaL:fitYLocal",cutS+"type==2&&refX==0","");  // phi shift - OK
    chain->Draw("deltaL:fitYLocal",cutS+"type==0&&refX==0","");  // OK
    chain->Draw("deltaL:fitYLocal",cutS+"type==0",""); // OK

  */
  
  return alignLocal;
}





void LoadTrees(){
  /// make  sector alignment - using Kalman filter method -AliTPCkalmanAlign
  /// Combined information is used, mean residuals are minimized:
  ///
  /// 1. TPC-ITS alignment
  /// 2. TPC vertex alignment

  TFile *f0 = new TFile("../mergeField0/mean.root");
  TFile *fP= new TFile("../mergePlus/mean.root");
  TFile *fM= new TFile("../mergeMinus/mean.root");
  //
  itsdy=(TTree*)f0->Get("ITSdy");
  itsdyP=(TTree*)fP->Get("ITSdy");
  itsdyM=(TTree*)fM->Get("ITSdy");
  itsdp=(TTree*)f0->Get("ITSdsnp");
  itsdphiP=(TTree*)fP->Get("ITSdsnp");
  itsdphiM=(TTree*)fM->Get("ITSdsnp");

  itsdz=(TTree*)f0->Get("ITSdz");
  itsdt=(TTree*)f0->Get("ITSdtheta");
  tofdy=(TTree*)f0->Get("TOFdy");
  trddy=(TTree*)f0->Get("TRDdy");
  //
  vdy=(TTree*)f0->Get("Vertexdy");
  vdyP=(TTree*)fP->Get("Vertexdy");
  vdyM=(TTree*)fM->Get("Vertexdy");
  vds=(TTree*)f0->Get("Vertexdsnp");
  vdz=(TTree*)f0->Get("Vertexdz");
  vdt=(TTree*)f0->Get("Vertexdtheta");
  chain    = new TChain("mean","mean");
  chainZ    = new TChain("mean","mean");
  chain->AddFile("../mergeField0/mean.root",10000000,"ITSdy");
  chain->AddFile("../mergeField0/mean.root",10000000,"ITSdsnp");
  chain->AddFile("../mergeField0/mean.root",10000000,"Vertexdy");
  chain->AddFile("../mergeField0/mean.root",10000000,"Vertexdsnp");
  chain->AddFile("../mergeField0/mean.root",10000000,"TOFdy");
  chain->AddFile("../mergeField0/mean.root",10000000,"TRDdy");
  //
  chainZ->AddFile("../mergeField0/mean.root",10000000,"Vertexdz");
  chainZ->AddFile("../mergeField0/mean.root",10000000,"Vertexdtheta");
  chainZ->AddFile("../mergeField0/mean.root",10000000,"ITSdz");
  chainZ->AddFile("../mergeField0/mean.root",10000000,"ITSdtheta");

  //
  itsdy->AddFriend(itsdp,"Phi");
  itsdy->AddFriend(itsdphiP,"PhiP");
  itsdy->AddFriend(itsdphiM,"PhiM");
  itsdy->AddFriend(itsdz,"Z");
  itsdy->AddFriend(itsdt,"T");
  itsdy->AddFriend(itsdyP,"YP");
  itsdy->AddFriend(itsdyM,"YM");
  //
  itsdy->AddFriend(vdy,"V");
  itsdy->AddFriend(vdyP,"VP");
  itsdy->AddFriend(vdyP,"VM");
  itsdy->AddFriend(tofdy,"TOF");
  itsdy->AddFriend(trddy,"TRD");
  itsdy->AddFriend(vds,"VPhi");
  itsdy->AddFriend(vdz,"VZ");
  itsdy->AddFriend(vdt,"VT");
  itsdy->AddFriend(tofdy,"TOF.");
  itsdy->SetMarkerStyle(25);
  itsdy->SetMarkerSize(0.4);
  tree=itsdy;
  tree->SetAlias("side","(-1+2*(theta>0))");
  chain->SetAlias("side","(-1+2*(theta>0))");
  chain->SetMarkerStyle(25);
  chain->SetMarkerSize(0.5);
}


void FitAlignCombinedCorr(){
  /// Fit Global X and globalY shift at vertex and at ITS

  RegisterAlignFunction();
  LoadTrees();
  combAlignOCDBOld = AliTPCCalibGlobalMisalignment::CreateOCDBAlign();
  //

  pcstream= new TTreeSRedirector("fitAlignCombined.root"); 

  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  // Fit vectors
  // global fits  
  TVectorD vecSec(37);
  TVectorD paramYGlobal;
  TVectorD paramYLocal;
  TMatrixD covar;
  // make Aliases
  {for (Int_t isec=0; isec<18; isec++){ // sectors
      tree->SetAlias(Form("sec%d",isec), Form("(abs(sector-%f)<0.5)",isec+0.5));
      chain->SetAlias(Form("sec%d",isec), Form("(abs(sector-%f)<0.5)",isec+0.5));
    }}
  for (Int_t isec=-1;isec<36; isec++){
    vecSec[isec+1]=isec;
  }
  chain->SetAlias("err","((type==0)*0.1+(type==2)*0.001)");
  // delta phi
  chain->SetAlias("dpgx","(0+0)");
  chain->SetAlias("dpgy","(0+0)");
  chain->SetAlias("dpgr0","((AliTPCCorrection::GetCorrSector(sector,245,theta,1,3)-AliTPCCorrection::GetCorrSector(sector,85,theta,1,3))/160)");
  chain->SetAlias("dpgr1","((AliTPCCorrection::GetCorrSector(sector,245,theta,1,4)-AliTPCCorrection::GetCorrSector(sector,85,theta,1,4))/160)");
  chain->SetAlias("dpgr2","((AliTPCCorrection::GetCorrSector(sector,245,theta,1,5)-AliTPCCorrection::GetCorrSector(sector,85,theta,1,5))/160)");
  chain->SetAlias("deltaPG","((AliTPCCorrection::GetCorrSector(sector,245,theta,1,100)-AliTPCCorrection::GetCorrSector(sector,85,theta,1,100))/160)");
  chain->SetAlias("deltaPL","((AliTPCCorrection::GetCorrSector(sector,245,theta,1,101)-AliTPCCorrection::GetCorrSector(sector,85,theta,1,101))/160)");
  // delta y at 85 cm
  chain->SetAlias("dygxT","(AliTPCCorrection::GetCorrSector(sector,85,theta,1,0)+0)"); 
  chain->SetAlias("dygyT","(AliTPCCorrection::GetCorrSector(sector,85,theta,1,1)+0)");
  chain->SetAlias("dyr0T","(AliTPCCorrection::GetCorrSector(sector,85,theta,1,3)+0)");
  chain->SetAlias("dyr1T","(AliTPCCorrection::GetCorrSector(sector,85,theta,1,4)+0)");
  chain->SetAlias("dyr2T","(AliTPCCorrection::GetCorrSector(sector,85,theta,1,5)+0)");
  chain->SetAlias("deltaYGT","(AliTPCCorrection::GetCorrSector(sector,85,theta,1,100)+0)");
  chain->SetAlias("deltaYLT","(AliTPCCorrection::GetCorrSector(sector,85,theta,1,101)+0)");
  //
  // delta y at reference X
  chain->SetAlias("dygx","(dygxT+0)");   // due global x shift
  chain->SetAlias("dygy","(dygyT+0)");   // due global y shift
  chain->SetAlias("dyr0","(dyr0T+dpgr0*(refX-85.))");  // due rotation 0
  chain->SetAlias("dyr1","(dyr1T+dpgr1*(refX-85.))");  // due rotation 1
  chain->SetAlias("dyr2","(dyr2T+dpgr2*(refX-85.))");  // due rotation 2
  chain->SetAlias("deltaYG","(deltaYGT+deltaPG*(refX-85.))");  // due global distortion
  chain->SetAlias("deltaYL","(deltaYLT+deltaPL*(refX-85.))");  // due local distortion

  chain->SetAlias("deltaG","(type==0)*deltaYG+(type==2)*deltaPG");  //alias to global function
  chain->SetAlias("deltaL","(type==0)*deltaYL+(type==2)*deltaPL");  //alias to local function
  //
  TString  fstringGlobal="";		 
  // global part
  fstringGlobal+="((type==0)*dygx+0)++";
  fstringGlobal+="((type==0)*dygy+0)++";
  fstringGlobal+="((type==0)*dyr0+(type==2)*dpgr0)++";
  fstringGlobal+="((type==0)*dyr1+(type==2)*dpgr1)++";
  fstringGlobal+="((type==0)*dyr2+(type==2)*dpgr2)++";
  //
  // Make global fits
  //
  TString *strFitYGlobal = TStatToolkit::FitPlane(chain,"mean:err", fstringGlobal.Data(),cutS+"", chi2,npoints,paramYGlobal,covar,-1,0, 10000000, kTRUE);
  combAlignGlobal =MakeAlignFunctionGlobal(paramYGlobal);
  printf("chi2=%f\n",TMath::Sqrt(chi2/npoints));
  chain->SetAlias("fitYGlobal",strFitYGlobal->Data());
  paramYGlobal.Print();
  paramYGlobal.Write("paramYGlobal");
  //
  //
  //
  //
  {for (Int_t isec=0; isec<18; isec++){
      //A side
      chain->SetAlias(Form("glyASec%d",isec),Form("((type==0)*sec%d*(theta>0))",isec)); // ly shift   
      chain->SetAlias(Form("gxASec%d",isec),Form("((type==0)*dygx*sec%d*(theta>0))",isec)); // gx shift
      chain->SetAlias(Form("gyASec%d",isec),Form("((type==0)*dygy*sec%d*(theta>0))",isec)); // gy shift
      chain->SetAlias(Form("gpASec%d",isec),Form("(((type==0)*dyr0+(type==2)*dpgr0)*sec%d)*(theta>0)",isec)); // phi rotation
      //C side
      chain->SetAlias(Form("glyCSec%d",isec),Form("((type==0)*sec%d*(theta<0))",isec)); // ly shift   
      chain->SetAlias(Form("gxCSec%d",isec),Form("((type==0)*dygx*sec%d*(theta<0))",isec)); 
      chain->SetAlias(Form("gyCSec%d",isec),Form("((type==0)*dygy*sec%d*(theta<0))",isec));
      chain->SetAlias(Form("gpCSec%d",isec),Form("(((type==0)*dyr0+(type==2)*dpgr0)*sec%d)*(theta<0)",isec)); // phi rotation
    }}

  TString  fstringLocal="";		 
  {for (Int_t isec=0; isec<18; isec++){   
      //fstringLocal+=Form("((type==0)*dygx*sec%d*(theta>0))++",isec);
      //      fstringLocal+=Form("(gxASec%d)++",isec);
      // fstringLocal+=Form("(gyASec%d)++",isec);
      fstringLocal+=Form("(glyASec%d)++",isec);
      fstringLocal+=Form("(gpASec%d)++",isec);
      fstringLocal+=Form("(glyCSec%d)++",isec);
      //fstringLocal+=Form("(gxCSec%d)++",isec);
      //fstringLocal+=Form("(gyCSec%d)++",isec);
      fstringLocal+=Form("(gpCSec%d)++",isec);
    }}
  		 
  TString *strFitYLocal = TStatToolkit::FitPlane(chain,"(mean-deltaG):err", fstringLocal.Data(),cutS+"", chi2,npoints,paramYLocal,covar,-1,0, 10000000, kTRUE);
  printf("chi2=%f\n",TMath::Sqrt(chi2/npoints));
  chain->SetAlias("fitYLocal",strFitYLocal->Data());
  combAlignLocal =MakeAlignFunctionSector(paramYLocal);
  //paramYLocal.Print();
  paramYLocal.Write("paramYLocal");
  //
  // scaling  - why do we need it?
  TString  fstringScale="";		 
  fstringScale+="((type==1)*refX+(type==2))*dsec++";
  fstringScale+="((type==1)*refX+(type==2))*dsec*cos(phi)++";
  fstringScale+="((type==1)*refX+(type==2))*dsec*sin(phi)++";
  fstringScale+="((type==1)*refX+(type==2))*dsec*side++";
  fstringScale+="((type==1)*refX+(type==2))*dsec*side*cos(phi)++";
  fstringScale+="((type==1)*refX+(type==2))*dsec*side*sin(phi)++";
  //
  //
  //
  pcstream->GetFile()->cd();
  DrawFitQA();  
  //
  AliTPCCalibGlobalMisalignment * combAlignOCDBNew0 =( AliTPCCalibGlobalMisalignment *)combAlignOCDBOld->Clone();  
  combAlignOCDBNew0->AddAlign(*combAlignGlobal);
  combAlignOCDBNew0->AddAlign(*combAlignLocal);
  combAlignOCDBNew= AliTPCCalibGlobalMisalignment::CreateMeanAlign(combAlignOCDBNew0);
  //
  combAlignGlobal->SetName("fcombAlignGlobal");
  combAlignLocal->SetName("fcombAlignLocal");
  combAlignOCDBOld->SetName("fcombAlignOCDBOld");
  combAlignOCDBNew->SetName("fcombAlignOCDBNew");
  combAlignOCDBNew0->SetName("fcombAlignOCDBNew0");
  //
  combAlignGlobal->Write("fcombAlignGlobal");
  combAlignLocal->Write("fcombAlignLocal");
  combAlignOCDBOld->Write("fcombAlignOCDBOld");
  combAlignOCDBNew->Write("fcombAlignOCDBNew");
  combAlignOCDBNew0->Write("fcombAlignOCDBNew0");
  combAlignOCDBNew->Print();
  delete pcstream;
}

void DrawFitQA(){
 /// MakeQA plot 1D

  TCanvas c;
  c.SetLeftMargin(0.15);
  chain->Draw("1000*(mean-deltaG)>>his(100,-1.5,1.5)",cutS+"type==2&&refX==0","");
  chain->GetHistogram()->SetName("TPC_VertexDphiGlobal1D");
  chain->GetHistogram()->SetTitle("TPC-Vertex #Delta_{#phi} - (Global Fit)");
  chain->GetHistogram()->GetXaxis()->SetTitle("#Delta_{#phi} (mrad)");
  chain->GetHistogram()->Fit("gaus","qnr");
  chain->GetHistogram()->Write();
  //
  chain->Draw("1000*(mean-deltaG)>>his(100,-1.5,1.5)",cutS+"type==2&&abs(refX-85)<2","");
  chain->GetHistogram()->SetName("TPC_ITSDphiGlobal1D");
  chain->GetHistogram()->SetTitle("TPC-ITS #Delta_{#phi} - (Global Fit)");
  chain->GetHistogram()->GetXaxis()->SetTitle("#Delta_{#phi} (mrad)");
  chain->GetHistogram()->Fit("gaus","qnr");
  chain->GetHistogram()->Write();
  //
  chain->Draw("10*(mean-deltaG)>>his(100,-1,1)",cutS+"type==0&&abs(refX-85)<2","");
  chain->GetHistogram()->SetName("TPC_ITSDRPhiGlobal1D");
  chain->GetHistogram()->SetTitle("TPC-ITS #Delta_{r#phi} - (Global Fit)");
  chain->GetHistogram()->GetXaxis()->SetTitle("#Delta_{r#phi} (mm)");
  chain->GetHistogram()->Fit("gaus","qnr");
  chain->GetHistogram()->Write();
  //
  chain->Draw("10*(mean-deltaG)>>his(100,-1,1)",cutS+"type==0&&abs(refX-0)<2","");
  chain->GetHistogram()->SetName("TPC-VertexDRPhiGlobal1D");
  chain->GetHistogram()->SetTitle("TPC-Vertex #Delta_{r#phi} - (Global Fit)");
  chain->GetHistogram()->GetXaxis()->SetTitle("#Delta_{r#phi} (mm)");
  chain->GetHistogram()->Fit("gaus","qnr");
  chain->GetHistogram()->Write();
  //
  // Make QA plot 3D
  //
  chain->Draw("1000*(mean-deltaG):sector:abs(theta)",cutS+"theta>0&&type==2&&refX==0","colz");
  chain->GetHistogram()->SetName("TPC_VertexDphi3DGlobal3D");
  chain->GetHistogram()->SetTitle("TPC-Vertex #Delta_{#phi} - (Global Fit)");
  chain->GetHistogram()->GetYaxis()->SetTitle("#Delta_{#phi} (mrad)");
  chain->GetHistogram()->GetXaxis()->SetTitle("sector");
  chain->GetHistogram()->GetZaxis()->SetTitle("tan(#Theta)");
  chain->GetHistogram()->Draw("colz");
  chain->GetHistogram()->Write();
  //
  chain->Draw("1000*(mean-deltaG):sector:abs(theta)",cutS+"theta>0&&type==2&&abs(refX-85)<2","colz");
  chain->GetHistogram()->SetName("TPC_ITSDphi3DGlobal3D");
  chain->GetHistogram()->SetTitle("TPC-ITS #Delta_{#phi} - (Global Fit)");
  chain->GetHistogram()->GetYaxis()->SetTitle("#Delta_{#phi} (mrad)");
  chain->GetHistogram()->GetXaxis()->SetTitle("sector");
  chain->GetHistogram()->GetZaxis()->SetTitle("tan(#Theta)");
  chain->GetHistogram()->Draw("colz");
  chain->GetHistogram()->Write();
  //
  chain->Draw("10*(mean-deltaG):sector:abs(theta)",cutS+"theta>0&&type==0&&abs(refX-85)<2","colz");
  chain->GetHistogram()->SetName("TPC_ITSDRphi3DGlobal3D");
  chain->GetHistogram()->SetTitle("TPC-ITS #Delta_{r#phi} - (Global Fit)");
  chain->GetHistogram()->GetYaxis()->SetTitle("#Delta_{r#phi} (mm)");
  chain->GetHistogram()->GetXaxis()->SetTitle("sector");
  chain->GetHistogram()->GetZaxis()->SetTitle("tan(#Theta)");
  chain->GetHistogram()->Draw("colz");
  chain->GetHistogram()->Write();
  //
  chain->Draw("10*(mean-deltaG):sector:abs(theta)",cutS+"theta>0&&type==0&&abs(refX-0)<2","colz");
  chain->GetHistogram()->SetName("TPC_VertexDRphi3DGlobal3D");
  chain->GetHistogram()->SetTitle("TPC-Vertex #Delta_{r#phi} - (Global Fit)");
  chain->GetHistogram()->GetYaxis()->SetTitle("#Delta_{r#phi} (mm)");
  chain->GetHistogram()->GetXaxis()->SetTitle("sector");
  chain->GetHistogram()->GetZaxis()->SetTitle("tan(#Theta)");
  chain->GetHistogram()->Draw("colz");
  chain->GetHistogram()->Write();
  //
  //
  //
  chain->Draw("1000*(mean-deltaG-deltaL):sector:abs(theta)",cutS+"theta>0&&type==2&&refX==0","colz");
  chain->GetHistogram()->SetName("TPC_VertexDphi3DLocal3D");
  chain->GetHistogram()->SetTitle("TPC-Vertex #Delta_{#phi} - (Local Fit)");
  chain->GetHistogram()->GetYaxis()->SetTitle("#Delta_{#phi} (mrad)");
  chain->GetHistogram()->GetXaxis()->SetTitle("sector");
  chain->GetHistogram()->GetZaxis()->SetTitle("tan(#Theta)");
  chain->GetHistogram()->Draw("colz");
  chain->GetHistogram()->Write();
  //
  chain->Draw("1000*(mean-deltaG-deltaL):sector:abs(theta)",cutS+"theta>0&&type==2&&abs(refX-85)<2","colz");
  chain->GetHistogram()->SetName("TPC_ITSDphi3DLocal3D");
  chain->GetHistogram()->SetTitle("TPC-ITS #Delta_{#phi} - (Local Fit)");
  chain->GetHistogram()->GetYaxis()->SetTitle("#Delta_{#phi} (mrad)");
  chain->GetHistogram()->GetXaxis()->SetTitle("sector");
  chain->GetHistogram()->GetZaxis()->SetTitle("tan(#Theta)");
  chain->GetHistogram()->Draw("colz");
  chain->GetHistogram()->Write();
  //
  chain->Draw("10*(mean-deltaG-deltaL):sector:abs(theta)",cutS+"theta>0&&type==0&&abs(refX-85)<2","colz");
  chain->GetHistogram()->SetName("TPC_ITSDRphi3DLocal3D");
  chain->GetHistogram()->SetTitle("TPC-ITS #Delta_{r#phi} - (Local Fit)");
  chain->GetHistogram()->GetYaxis()->SetTitle("#Delta_{r#phi} (mm)");
  chain->GetHistogram()->GetXaxis()->SetTitle("sector");
  chain->GetHistogram()->GetZaxis()->SetTitle("tan(#Theta)");
  chain->GetHistogram()->Draw("colz");
  chain->GetHistogram()->Write();
  //
  chain->Draw("10*(mean-deltaG-deltaL):sector:abs(theta)",cutS+"theta>0&&type==0&&abs(refX-0)<2","colz");
  chain->GetHistogram()->SetName("TPC_VertexDRphiLocal3D");
  chain->GetHistogram()->SetTitle("TPC-Vertex #Delta_{r#phi} - (Local Fit)");
  chain->GetHistogram()->GetYaxis()->SetTitle("#Delta_{r#phi} (mm)");
  chain->GetHistogram()->GetXaxis()->SetTitle("sector");
  chain->GetHistogram()->GetZaxis()->SetTitle("tan(#Theta)");
  chain->GetHistogram()->Draw("colz");
  chain->GetHistogram()->Write();
}




void FitAlignCombined0(){
  /// Fit Global X and globalY shift at vertex and at ITS

  TTreeSRedirector *pcstream= new TTreeSRedirector("fitAlignCombined.root"); 

  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  // Fit vectors
  // global fits  
  TVectorD vecSec(37);
  TVectorD paramYVGlobal;
  TVectorD paramYITSGlobal;
  TVectorD paramPhiVGlobal;
  TVectorD paramPhiITSGlobal;
  // local fits
  TVectorD paramYVLocal;
  TVectorD paramPhiVLocal;
  TVectorD paramYITSLocal;
  TVectorD paramPhiITSLocal;
  TMatrixD covar;
  // make Aliases
  {for (Int_t isec=0; isec<18; isec++){ // sectors
      tree->SetAlias(Form("sec%d",isec), Form("(abs(sector-%f)<0.5)",isec+0.5));
    }}
  for (Int_t isec=-1;isec<36; isec++){
    vecSec[isec+1]=isec;
  }
  //
  TString  fstringGlobal="";		 
  fstringGlobal+="side++";
  fstringGlobal+="theta++";
  fstringGlobal+="cos(phi)*(theta>0)++";  // GX - A side
  fstringGlobal+="sin(phi)*(theta>0)++";  // GY - A side
  fstringGlobal+="cos(phi)*(theta<0)++";  // GX - C side
  fstringGlobal+="sin(phi)*(theta<0)++";  // GY - C side
  fstringGlobal+="cos(phi)*(theta)++";    // theta - get rid of rotation 
  fstringGlobal+="sin(phi)*(theta)++";    // theta - get rid of rotation
  //
  // Local Fit
  //  
  TString  fstringLocal="";		 
  {for (Int_t sector=0; sector<18; sector++){   
      fstringLocal+=Form("((sec%d)*theta>0)++",sector);
      fstringLocal+=Form("((sec%d)*theta<0)++",sector);
    }}
  //
  // Make global fits
  //
  TString *strFitYVGlobal = TStatToolkit::FitPlane(tree,"V.mean", fstringGlobal.Data(),cutS+"", chi2,npoints,paramYVGlobal,covar,-1,0, 10000000, kFALSE);
  printf("chi2=%f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitYVGlobal",strFitYVGlobal->Data());
  tree->Draw("V.mean:fitYVGlobal",cutS);
  //
  TString *strFitPhiVGlobal = TStatToolkit::FitPlane(tree,"VPhi.mean", fstringGlobal.Data(),cutS+"", chi2,npoints,paramPhiVGlobal,covar,-1,0, 10000000, kFALSE);
  printf("chi2=%f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitPhiVGlobal",strFitPhiVGlobal->Data());
  tree->Draw("VPhi.mean:fitPhiVGlobal",cutS);
  //
  TString *strFitYITSGlobal = TStatToolkit::FitPlane(tree,"mean", fstringGlobal.Data(),cutS+"", chi2,npoints,paramYITSGlobal,covar,-1,0, 10000000, kFALSE);
  printf("chi2=%f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitYITSGlobal",strFitYITSGlobal->Data());
  tree->Draw("mean:fitYITSGlobal",cutS);
  TString *strFitPhiITSGlobal = TStatToolkit::FitPlane(tree,"Phi.mean", fstringGlobal.Data(),cutS+"", chi2,npoints,paramPhiITSGlobal,covar,-1,0, 10000000, kFALSE);
  printf("chi2=%f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitPhiITSGlobal",strFitPhiITSGlobal->Data());
  tree->Draw("Phi.mean:fitPhiITSGlobal",cutS);
  //
  // Residuals to the global fit
  //
  tree->SetAlias("dyVG", "V.mean-fitYVGlobal");      
  tree->SetAlias("dphiVG","(VPhi.mean-fitPhiVGlobal)");
  tree->SetAlias("dyITSG","mean-fitYITSGlobal");
  tree->SetAlias("dphiITSG","(Phi.mean-fitPhiITSGlobal)");

  TString *strFitYVLocal = TStatToolkit::FitPlane(tree,"dyVG", fstringLocal.Data(),cutS+"", chi2,npoints,paramYVLocal,covar,-1,0, 10000000, kFALSE);
  printf("chi2=%f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitYVLocal",strFitYVLocal->Data());
  tree->Draw("dyVG-fitYVLocal:sector:abs(theta)",cutS+"theta>0","colz");  
  //
  TString *strFitPhiVLocal = TStatToolkit::FitPlane(tree,"dphiVG", fstringLocal.Data(),cutS+"", chi2,npoints,paramPhiVLocal,covar,-1,0, 10000000, kFALSE);
  printf("chi2=%f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitPhiVLocal",strFitPhiVLocal->Data());
  tree->Draw("dphiVG-fitPhiVLocal:sector:abs(theta)",cutS+"theta>0","colz");  
  //
  //
  //
  TString *strFitYITSLocal = TStatToolkit::FitPlane(tree,"dyITSG", fstringLocal.Data(),cutS+"", chi2,npoints,paramYITSLocal,covar,-1,0, 10000000, kFALSE);
  printf("chi2=%f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitYITSLocal",strFitYITSLocal->Data());
  tree->Draw("dyITSG-fitYITSLocal:sector:abs(theta)",cutS+"theta>0","colz");
  //
  TString *strFitPhiITSLocal = TStatToolkit::FitPlane(tree,"dphiITSG", fstringLocal.Data(),cutS+"", chi2,npoints,paramPhiITSLocal,covar,-1,0, 10000000, kFALSE);
  printf("chi2=%f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitPhiITSLocal",strFitPhiITSLocal->Data());
  tree->Draw("dphiITSG-fitPhiITSLocal:sector:abs(theta)",cutS+"theta>0","colz");

  //
  TH1D *hisA = 0;
  TH1D *hisC = 0;
  TVectorD vSec(36);
  TVectorD vecDPhi(36);
  TVectorD vecDLY(36);
  TVectorD vecDGX(36);
  TVectorD vecDGY(36);

  //
  tree->SetAlias("phiMean","(fitPhiVLocal+fitPhiITSLocal+fitPhiVGlobal+fitPhiITSGlobal)*0.5");
  tree->SetAlias("yMeanV","((fitYITSLocal+fitYITSGlobal-85*phiMean)+(fitYVLocal+fitYVGlobal))*0.5");
  tree->SetAlias("yMeanITS","((fitYITSLocal+fitYITSGlobal)+(fitYVLocal+fitYVGlobal+85*phiMean))*0.5");
  

  tree->Draw("phiMean:int(sector)>>hhhA(18,0.18)",cutS+"abs(theta)<0.15&&theta>0","prof"); 
  hisA = (TH1D*) tree->GetHistogram()->Clone();
  tree->Draw("phiMean:int(sector)>>hhhC(18,0.18)",cutS+"abs(theta)<0.15&&theta<0","prof"); 
  hisC = (TH1D*) tree->GetHistogram()->Clone();
  
  for (Int_t isec=0; isec<36; isec++){
    vSec[isec]=isec;
    if (isec<18)  vecDPhi[isec]=hisA->GetBinContent(isec%18+1);
    if (isec>=18) vecDPhi[isec]=hisC->GetBinContent(isec%18+1);    
  }
  tree->Draw("yMeanV:int(sector)>>hhhA(18,0.18)",cutS+"abs(theta)<0.15&&theta>0","prof"); 
  hisA = (TH1D*) tree->GetHistogram()->Clone();
  tree->Draw("yMeanV:int(sector)>>hhhC(18,0.18)",cutS+"abs(theta)<0.15&&theta<0","prof"); 
  hisC = (TH1D*) tree->GetHistogram()->Clone();
  //
  for (Int_t isec=0; isec<36; isec++){
    if (isec<18)  vecDLY[isec]=hisA->GetBinContent(isec%18+1);
    if (isec>=18) vecDLY[isec]=hisC->GetBinContent(isec%18+1);    
    vecDGX[isec]=-vecDLY[isec]*TMath::Sin(TMath::Pi()*(isec+0.5)/9.);
    vecDGY[isec]= vecDLY[isec]*TMath::Cos(TMath::Pi()*(isec+0.5)/9.);
  }

 


  //
  // Store results 
  //
  {(*pcstream)<<"fitLocal"<<
      //global fits
      "sec.="<<&vSec<<
      "pYVGlobal.="<<&paramYVGlobal<<
      "pYITSGlobal.="<<&paramYITSGlobal<<
      "pPhiVGlobal.="<<&paramPhiVGlobal<<
      "pPhiITSGlobal.="<<&paramPhiITSGlobal<<
      // local fits
      "pYVLocal.="<<&paramYVLocal<<
      "pPhiVLocal.="<<&paramPhiVLocal<<
      "pYITSLocal.="<<&paramYITSLocal<<
      "pPhiITSLocal.="<<&paramPhiITSLocal<<
      //
      // combined
      "vDPhi.="<<&vecDPhi<<
      "vDLY.="<<&vecDLY<<
      "vDGX.="<<&vecDGX<<
      "vDGY.="<<&vecDGY<<
      "\n";
  }
  paramYVGlobal.Write("paramYVGlobal");
  paramYITSGlobal.Write("paramYITSGlobal");
  paramPhiVGlobal.Write("paramPhiVGlobal");
  paramPhiITSGlobal.Write("paramPhiITSGlobal");
  // local fits
  paramYVLocal.Write("paramYVLocal");
  paramPhiVLocal.Write("paramPhiVLocal");
  paramYITSLocal.Write("paramYITSLocal");
  paramPhiITSLocal.Write("paramPhiITSLocal");
  vecDPhi.Write("vecDPhi");
  vecDLY.Write("vecDLY");
  vecDGX.Write("vecDGX");
  vecDGY.Write("vecDGY");



  tree->Draw("dyVG:sector:abs(theta)",cutS+"theta>0","colz");
  tree->GetHistogram()->Write("fitYVGlobal");
  tree->Draw("dphiVG:sector:abs(theta)",cutS+"theta>0","colz");
  tree->GetHistogram()->Write("fitPhiVGlobal");
  tree->Draw("dyITSG:sector:abs(theta)",cutS+"theta>0","colz");
  tree->GetHistogram()->Write("fitYITSGlobal");
  tree->Draw("dphiITSG:sector:abs(theta)",cutS+"theta>0","colz");
  tree->GetHistogram()->Write("fitPhiITSGlobal");
  //
  tree->Draw("dyVG-fitYVLocal:sector:abs(theta)",cutS+"theta>0","colz");
  tree->GetHistogram()->Write("fitYVLocal");
  tree->Draw("dphiVG-fitPhiVLocal:sector:abs(theta)",cutS+"theta>0","colz");
  tree->GetHistogram()->Write("fitPhiVLocal");
  tree->Draw("dyITSG-fitYITSLocal:sector:abs(theta)",cutS+"theta>0","colz");
  tree->GetHistogram()->Write("fitYITSLocal");
  tree->Draw("dphiITSG-fitPhiITSLocal:sector:abs(theta)",cutS+"theta>0","colz");
  tree->GetHistogram()->Write("fitPhiITSLocal");
  delete pcstream;
}

void FitAlignCombined(){
  /// make  sector alignment - using Kalman filter method -AliTPCkalmanAlign
  /// Combined information is used, mean residuals are minimized:
  ///
  /// 1. TPC-TPC sector alignment
  /// 2. TPC-ITS alignment
  /// 3. TPC vertex alignment

  TFile fcalib("../mergeField0/TPCAlignObjects.root");
  AliTPCcalibAlign * align = ( AliTPCcalibAlign *)fcalib.Get("alignTPC");

  TCut cutQ="entries>1000&&abs(theta)<0.5&&abs(snp)<0.2&&YP.entries>50&&YM.entries>50";
  TH1F his1("hdeltaY1","hdeltaY1",100,-0.5,0.5);
  TMatrixD vecAlign(72,1);
  TMatrixD covAlign(72,72);
  TMatrixD vecAlignY(72,1);
  TMatrixD covAlignY(72,72);
  TMatrixD vecAlignTheta(72,1);
  TMatrixD covAlignTheta(72,72);
  TMatrixD vecAlignZ(72,1);
  TMatrixD covAlignZ(72,72);
  AliTPCkalmanAlign::BookAlign1D(vecAlign,covAlign,0,0.002);
  AliTPCkalmanAlign::BookAlign1D(vecAlignY,covAlignY,0,0.1);
  AliTPCkalmanAlign::BookAlign1D(vecAlignTheta,covAlignTheta,0,0.002);
  AliTPCkalmanAlign::BookAlign1D(vecAlignZ,covAlignZ,0,0.1);
  TVectorD vecITSY(72);
  TVectorD vecITSYPM(72);
  TVectorD vecITSPhi(72);
  TVectorD vecITSPhiPM(72);
  TVectorD vecVY(72);
  TVectorD vecVS(72);
  TVectorD vecITSTheta(72);
  TVectorD vecVTheta(72);
  {for (Int_t isec0=0; isec0<36; isec0++){
      Double_t phi0=(isec0%18+0.5)*TMath::Pi()/9.;
      if (phi0>TMath::Pi()) phi0-=TMath::TwoPi();
      Int_t iside0=(isec0%36<18)? 0:1;
      TCut cutSector=Form("abs(%f-phi)<0.14",phi0);
      TCut cutSide = (iside0==0)? "theta>0":"theta<0";
      //
      itsdy->Draw("mean",cutQ+cutSector+cutSide);
      Double_t meanITSY=itsdy->GetHistogram()->GetMean();
      vecITSY[isec0]=meanITSY;
      vecITSY[isec0+36]=meanITSY;
      //
      itsdy->Draw("(YP.mean+YM.mean)*0.5",cutQ+cutSector+cutSide);
      Double_t meanITSYPM=itsdy->GetHistogram()->GetMean();
      vecITSYPM[isec0]=meanITSYPM;
      vecITSYPM[isec0+36]=meanITSYPM;
      //
      itsdy->Draw("Phi.mean",cutQ+cutSector+cutSide);
      Double_t meanITSPhi=itsdy->GetHistogram()->GetMean();
      vecITSPhi[isec0]=meanITSPhi;
      vecITSPhi[isec0+36]=meanITSPhi;
      //
      itsdy->Draw("(PhiP.mean+PhiM.mean)*0.5",cutQ+cutSector+cutSide);
      Double_t meanITSPhiPM=itsdy->GetHistogram()->GetMean();
      vecITSPhiPM[isec0]=meanITSPhiPM;
      vecITSPhiPM[isec0+36]=meanITSPhiPM;
      //
      //
      itsdy->Draw("VPhi.mean",cutQ+cutSector+cutSide);
      Double_t meanVS=itsdy->GetHistogram()->GetMean();
      vecVS[isec0]=meanVS;
      vecVS[isec0+36]=meanVS;
      //
      itsdy->Draw("V.mean",cutQ+cutSector+cutSide);
      Double_t meanVY=itsdy->GetHistogram()->GetMean();
      vecVY[isec0]=meanVY;
      vecVY[isec0+36]=meanVY;
      //
      itsdy->Draw("T.mean",cutQ+cutSector+cutSide);
      Double_t meanITST=itsdy->GetHistogram()->GetMean();
      vecITSTheta[isec0]=meanITST;
      vecITSTheta[isec0+36]=meanITST;
      // 
      itsdy->Draw("VT.mean",cutQ+cutSector+cutSide);
      Double_t meanVT=itsdy->GetHistogram()->GetMean();
      vecVTheta[isec0]=meanVT;
      vecVTheta[isec0+36]=meanVT;
    }
  }
  AliTPCROC * roc = AliTPCROC::Instance();
  Double_t fXIROC = (roc->GetPadRowRadii(0,0)+roc->GetPadRowRadii(0,roc->GetNRows(0)-1))*0.5;
  Double_t fXOROC = (roc->GetPadRowRadii(36,0)+roc->GetPadRowRadii(36,roc->GetNRows(36)-1))*0.5;
  Double_t fXIO       = ( roc->GetPadRowRadii(0,roc->GetNRows(0)-1)+roc->GetPadRowRadii(36,0))*0.5;

  TTreeSRedirector *pcstream=new TTreeSRedirector("combAlign.root");
  
  {
    for (Int_t iter=0; iter<5; iter++){
    for (Int_t isec0=0; isec0<72; isec0++){
    for (Int_t isec1=0; isec1<72; isec1++){
      //if (isec0%36!=isec0%36) continue;
      if (iter==0 && isec0%36!=isec1%36) continue;
      if (iter==1 && isec0%18!=isec1%18) continue;
      TH1 * his = align->GetHisto(AliTPCcalibAlign::kY,isec0,isec1);
      TH1 * hisPhi = align->GetHisto(AliTPCcalibAlign::kPhi,isec0,isec1);
      TH1 * hisTheta = align->GetHisto(AliTPCcalibAlign::kTheta,isec0,isec1);
      TH1 * hisZ = align->GetHisto(AliTPCcalibAlign::kZ,isec0,isec1);
      Int_t side0=(isec0%36)<18 ? 1:-1; 
      Int_t side1=(isec1%36)<18 ? 1:-1;       
      Double_t weightTPC= 0.002;
      if (isec0%18==isec1%18) weightTPC=0.0005;
      if (isec0%36==isec1%36) weightTPC=0.00005;
      if (side0!=side1) weightTPC*=2.;
      Double_t weightTPCT= 0.001;
      if (isec0%18==isec1%18) weightTPC=0.0005;
      if (isec0%36==isec1%36) weightTPC=0.00005;
      if (TMath::Abs(isec0%18-isec1%18)==1) weightTPC = 0.0005;
      if (TMath::Abs(isec0%36-isec1%36)==1) weightTPC = 0.0001;

      //
      Double_t meanITS0=vecITSY[isec0];
      Double_t meanITSPM0=vecITSYPM[isec0];
      Double_t meanITSPhi0=vecITSPhi[isec0];
      Double_t meanITSPhiPM0=vecITSPhiPM[isec0];
      Double_t meanVY0=vecVY[isec0];
      Double_t meanVPhi0=vecVS[isec0];
      Double_t meanITSTheta0=vecITSTheta[isec0];
      Double_t meanVTheta0=vecVTheta[isec0];
      //
      Double_t meanITS1=vecITSY[isec1];
      Double_t meanITSPM1=vecITSYPM[isec1];
      Double_t meanITSPhi1=vecITSPhi[isec1];
      Double_t meanITSPhiPM1=vecITSPhiPM[isec1];
      Double_t meanVY1=vecVY[isec1];
      Double_t meanVPhi1=vecVS[isec1];
      Double_t meanITSTheta1=vecITSTheta[isec1];
      Double_t meanVTheta1=vecVTheta[isec1];
      //
      Double_t meanPhi0 = (meanITSPhi0+meanVPhi0+2*meanITS0/83.)/4.;
      Double_t meanPhi1 = (meanITSPhi1+meanVPhi1+2*meanITS1/83.)/4.;
      //
      //
      if (iter>2 &&isec0==isec1){
	AliTPCkalmanAlign::UpdateAlign1D(-meanITS0/83, 0.001, isec0,  vecAlign,covAlign);
	AliTPCkalmanAlign::UpdateAlign1D(-meanITSPhi0, 0.001, isec0,  vecAlign,covAlign);
	AliTPCkalmanAlign::UpdateAlign1D(-meanVPhi0, 0.001, isec0,  vecAlign,covAlign);
	//
	AliTPCkalmanAlign::UpdateAlign1D(-meanPhi0, 0.001, isec0,  vecAlign,covAlign);
	AliTPCkalmanAlign::UpdateAlign1D(-meanVTheta0, 0.001, isec0,  vecAlignTheta,covAlignTheta);
	AliTPCkalmanAlign::UpdateAlign1D(-meanITSTheta0, 0.001, isec0,  vecAlignTheta,covAlignTheta);
      }
      if (iter>2){
	AliTPCkalmanAlign::UpdateAlign1D(meanPhi1-meanPhi0, 0.001, isec0,isec1,  vecAlign,covAlign);
	AliTPCkalmanAlign::UpdateAlign1D(meanITSPhiPM1-meanITSPhiPM0, 0.001, isec0,isec1,  vecAlign,covAlign);
	AliTPCkalmanAlign::UpdateAlign1D((meanITSPM1-meanITSPM0)/83., 0.001, isec0,isec1,  vecAlign,covAlign);
	AliTPCkalmanAlign::UpdateAlign1D(meanVTheta1-meanVTheta0, 0.001, isec0,isec1,  vecAlignTheta,covAlignTheta);
	AliTPCkalmanAlign::UpdateAlign1D(meanITSTheta1-meanITSTheta0, 0.001, isec0,isec1,  vecAlignTheta,covAlignTheta);
      }
      //
      if (!his) continue;
      if (!hisPhi) continue;
      if (!hisTheta) continue;
      if (his->GetEntries()<100) continue;
      Double_t xref=fXIO;
      if (isec0<36&&isec1<36) xref=fXIROC;
      if (isec0>=36&&isec1>=36) xref=fXOROC;
      Double_t meanTPC=his->GetMean();
      Double_t meanTPCPhi=hisPhi->GetMean();
      Double_t meanTPCTheta=hisTheta->GetMean();
      Double_t meanTPCZ=hisZ->GetMean();
      //
      //
      Double_t kalman0= vecAlign(isec0,0);
      Double_t kalman1= vecAlign(isec1,0);
      Double_t kalmanY0= vecAlignY(isec0,0);
      Double_t kalmanY1= vecAlignY(isec1,0);
      Double_t kalmanTheta0= vecAlignTheta(isec0,0);
      Double_t kalmanTheta1= vecAlignTheta(isec1,0);
      Double_t kalmanZ0= vecAlignZ(isec0,0);
      Double_t kalmanZ1= vecAlignZ(isec1,0);
      //
      //
      AliTPCkalmanAlign::UpdateAlign1D((meanTPC)/xref, weightTPC, isec0,isec1,  vecAlign,covAlign);
      AliTPCkalmanAlign::UpdateAlign1D(meanTPCPhi, weightTPC*10,isec0,isec1,  vecAlign,covAlign);
      //
      AliTPCkalmanAlign::UpdateAlign1D(meanTPCTheta, weightTPCT, isec0,isec1,  vecAlignTheta,covAlignTheta);
      if (side0==side1) AliTPCkalmanAlign::UpdateAlign1D(meanTPCZ, weightTPCT*100., isec0,isec1,  vecAlignZ,covAlignZ);      
      //printf("isec0\t%d\tisec1\t%d\t%f\t%f\t%f\n",isec0,isec1, meanTPC, meanITS0,meanITS1);

      if (iter>=0) (*pcstream)<<"align"<<
	"iter="<<iter<<           
	"xref="<<xref<<                 // reference X
	"isec0="<<isec0<<               // sector number
	"isec1="<<isec1<<
	"side0="<<side0<<
	"side1="<<side1<<
	//TPC
	"mTPC="<<meanTPC<<              // delta Y  / xref
	"mTPCPhi="<<meanTPCPhi<<        // delta Phi
	"mTPCZ="<<meanTPCZ<<        // delta Z
	"mTPCTheta="<<meanTPCTheta<<    // delta Theta
	//ITS
	"mITS0="<<meanITS0<<  
	"mITS1="<<meanITS1<<
	"mITSPhi0="<<meanITSPhi0<<
	"mITSPhi1="<<meanITSPhi1<<
	//
	"mITSPM0="<<meanITSPM0<<  
	"mITSPM1="<<meanITSPM1<<
	"mITSPhiPM0="<<meanITSPhiPM0<<
	"mITSPhiPM1="<<meanITSPhiPM1<<
	//
	"mITSTheta0="<<meanITSTheta0<<
	"mITSTheta1="<<meanITSTheta1<<
	//Vertex
	"mVY0="<<meanVY0<<
	"mVY1="<<meanVY1<<
	"mVPhi0="<<meanVPhi0<<
	"mVPhi1="<<meanVPhi1<<
	"mVTheta0="<<meanVTheta0<<
	"mVTheta1="<<meanVTheta1<<
	// Vertex+ITS mean
	"mPhi0="<<meanPhi0<<
	"mPhi1="<<meanPhi1<<
	//Kalman
	"kY0="<<kalmanY0<<
	"kY1="<<kalmanY1<<
	"kPhi0="<<kalman0<<
	"kPhi1="<<kalman1<<
	"kTheta0="<<kalmanTheta0<<
	"kTheta1="<<kalmanTheta1<<
	"kZ0="<<kalmanZ0<<
	"kZ1="<<kalmanZ1<<
	"\n";
    }          
    }
    }
  }
  pcstream->GetFile()->cd();
  vecAlign.Write("alignPhiMean");
  covAlign.Write("alingPhiCovar");
  vecAlignTheta.Write("alignThetaMean");
  covAlignTheta.Write("alingThetaCovar");
  vecAlignZ.Write("alignZMean");
  covAlignZ.Write("alingZCovar");
  delete pcstream;
  TFile f("combAlign.root");
  TTree * treeA = (TTree*)f.Get("align"); 
  treeA->SetMarkerStyle(25);
  treeA->SetMarkerSize(0.5);
}





void UpdateOCDBAlign(){
  /// Store resulting OCDB entry
  /// 0. Setup OCDB to get necccessary old entries - not done here
  /// 1. Get old OCDB entry
  /// 2. Get delta alignment
  /// 3. Add delta alignment
  /// 4. Store new alignment in

  AliCDBEntry * entry = AliCDBManager::Instance()->Get("TPC/Align/Data");
  TClonesArray * array = (TClonesArray*)entry->GetObject();
  Int_t entries = array->GetEntries();
  TClonesArray * arrayNew=(TClonesArray*)array->Clone();
  TFile f("combAlign.root");
  TMatrixD *matPhi=(TMatrixD*)f.Get("alignPhiMean");  
  //
  //
  { for (Int_t i=0;i<entries; i++){
      //
      //
      AliAlignObjParams *alignP = (AliAlignObjParams*)array->UncheckedAt(i);
      AliAlignObjParams *alignPNew = (AliAlignObjParams*)arrayNew->UncheckedAt(i);
      Int_t imod;
      AliGeomManager::ELayerID ilayer;
      alignP->GetVolUID(ilayer, imod);
      if (ilayer==AliGeomManager::kInvalidLayer) continue;
      Int_t sector=imod;
      if (ilayer==AliGeomManager::kTPC2) sector+=36;
      Double_t transOld[3], rotOld[3];
      alignP->GetTranslation(transOld);  // in cm 
      alignP->GetAngles(rotOld);         // in degrees
      printf("%d\t%d\t%d\t",ilayer, imod,sector);
      printf("%f mrad \t%f mrad\n",1000*rotOld[2]*TMath::DegToRad(), 1000*((*matPhi)(sector,0)));
      alignPNew->SetRotation(rotOld[0],rotOld[1], rotOld[2]+(*matPhi)(sector,0)*TMath::RadToDeg());
      alignPNew->Print("kokot");
      alignP->Print("kokot");	
    }
  }

  

  TString ocdbStorage="local://"+gSystem->GetFromPipe("pwd")+"/OCDBUpdate";
  TString userName=gSystem->GetFromPipe("echo $USER");
  TString date=gSystem->GetFromPipe("date");

  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Marian Ivanov");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-26-02"); //root version  
  metaData->SetComment(Form("Correction calibration. User: %s\n Data%s",userName.Data(),date.Data()));

  AliCDBId* id1=NULL;
  id1=new AliCDBId("TPC/Align/Data", 0, AliCDBRunRange::Infinity());
  AliCDBStorage* gStorage = AliCDBManager::Instance()->GetStorage(ocdbStorage);
  gStorage->Put(arrayNew, (*id1), metaData);

  TFile falign("falign.root","recreate");
  arrayNew->Write("new");
  array->Write("old");
  falign.Close();
}



void UpdateOCDBAlign0(){
  /// Store resulting OCDB entry
  /// 0. Setup OCDB to get necccessary old entries - not done here
  /// 1. Get old OCDB entry
  /// 2. Get delta alignment
  /// 3. Add delta alignment
  /// 4. Store new alignment in

  AliCDBEntry * entry = AliCDBManager::Instance()->Get("TPC/Align/Data");
  TClonesArray * array = (TClonesArray*)entry->GetObject();
  Int_t entries = array->GetEntries();
  TClonesArray * arrayNew=(TClonesArray*)array->Clone();
  TFile f("fitAlignCombined.root");
  TVectorD *vecDPhi = (TVectorD*) f.Get("vecDPhi");
  TVectorD *vecDGX = (TVectorD*) f.Get("vecDGX");
  TVectorD *vecDGY = (TVectorD*) f.Get("vecDGY");
  
  //
  //
  { for (Int_t i=0;i<entries; i++){
      //
      //
      AliAlignObjParams *alignP = (AliAlignObjParams*)array->UncheckedAt(i);
      AliAlignObjParams *alignPNew = (AliAlignObjParams*)arrayNew->UncheckedAt(i);
      Int_t imod;
      AliGeomManager::ELayerID ilayer;
      alignP->GetVolUID(ilayer, imod);
      if (ilayer==AliGeomManager::kInvalidLayer) continue;
      Int_t sector=imod;
      if (ilayer==AliGeomManager::kTPC2) sector+=36;
      Double_t transOld[3], rotOld[3];
      alignP->GetTranslation(transOld);  // in cm 
      alignP->GetAngles(rotOld);         // in degrees
      printf("%d\t%d\t%d\t",ilayer, imod,sector);
      printf("%f mrad \t%f mrad\n",1000*rotOld[2]*TMath::DegToRad(), 1000*((*vecDPhi)[sector%36]));
      alignPNew->SetRotation(rotOld[0],rotOld[1], rotOld[2]-(*vecDPhi)[sector%36]*TMath::RadToDeg());
      alignPNew->SetTranslation(transOld[0]-(*vecDGX)[sector%36],transOld[1]-(*vecDGY)[sector%36], transOld[2]);
      alignPNew->Print("kokot");
      alignP->Print("kokot");	
    }
  }

  

  TString ocdbStorage="local://"+gSystem->GetFromPipe("pwd")+"/OCDBUpdate";
  TString userName=gSystem->GetFromPipe("echo $USER");
  TString date=gSystem->GetFromPipe("date");

  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Marian Ivanov");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-26-02"); //root version  
  metaData->SetComment(Form("Correction calibration. User: %s\n Data%s",userName.Data(),date.Data()));

  AliCDBId* id1=NULL;
  id1=new AliCDBId("TPC/Align/Data", 0, AliCDBRunRange::Infinity());
  AliCDBStorage* gStorage = AliCDBManager::Instance()->GetStorage(ocdbStorage);
  gStorage->Put(arrayNew, (*id1), metaData);

  TFile falign("falign.root","recreate");
  arrayNew->Write("new");
  array->Write("old");
  falign.Close();
}
