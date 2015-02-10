/// \file FitRodShift.C
///
/// \author marian.ivanov@cern.ch, Stefan.Rossegger@cern.ch
///
/// Macro to fit  alignment/E field distortion maps
/// As a input output cluster distortion maps (produced by  AliTPCcalibAling class)
/// are used. Cluster distrotion maps granularity (180 *phi x44 *theta x 53 R)
///
/// In total 440 parameters to fit using global fit:
///
/// 1. Rotation and translation for each sector (72x2)
/// 2. Rotation and translation for each quadrant of OROC (36x4x2)
/// 3. Rod/strip shifts in IFC and OFC (18 sectors x 2 sides x 2 ) 
/// 4. Rotated clips in IFC and OFC 
///
///
/// Input file mean.root with distortion tree expected to be in directory:
/// ../mergeField0/clusterDY.root
/// The ouput file fitRodShift.root contains:
/// 1. Resulting (residual) AliTPCCalibMisalignment  and AliTPCFCVoltError3D classes
/// 2. All important temporary results are stored in the workspace (TTreeSRedirector)  associated
///    with the file - fitrodShift
///    2.a  Fit parameters with errors
///    2.b  QA default graphs
///    2.c  QA defualt histograms
///
///
///
/// Functions:
/// 1. LoadModels()                   - load models to fit - Rod shift (10,20)
/// 2. RegisterAlignFunction()        - register align functions (10-16)
/// 3. LoadTrees                      - load trees and make aliases
/// 4. PrintFit                       - helper function to print substring of the fit 
/// 5. DeltaLookup                    - function to calulate the distortion for given distortion cunction
///                                   -
/// 6. MakeAliases                    - Make tree aliases shortcuts for the fitting function
///
/// 7. MakeAlignCorrection            - Crete alignment entry to be stored in the OCDB
///                                   - Dump fit parameters and errors
/// 8. MakeQuadrantCorrection         - 
///                                   - Dump fit parameters and errors
///
/// 9. FitRodShifts                   - main minimization function
///
/// ~~~
/// .x ~/rootlogon.C
/// .x ~/NimStyle.C
/// .L $ALICE_ROOT/TPC/CalibMacros/FitRodShift.C+
/// FitRodShift(kFALSE);
/// ~~~

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"
#include "TTreeStream.h"
#include "TMath.h" 
#include "TGraph.h" 
#include "TRandom.h"
#include "TTree.h"
#include "TF1.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TAxis.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "AliTPCParamSR.h"
#include "TDatabasePDG.h"
#include "AliTPCFCVoltError3D.h"
#include "AliTPCROCVoltError3D.h"
#include "AliTPCComposedCorrection.h"
#include "AliTPCBoundaryVoltError.h"
#include "AliCDBEntry.h"
#include "AliTPCROC.h"
#include <TStatToolkit.h>
#include "TCut.h"
#include "TGraphErrors.h"
#include "AliTPCCalibGlobalMisalignment.h"
#endif

TTreeSRedirector *pcWorkspace= 0;    // Workspace to store the fit parameters and QA histograms
TTree * treeDY= 0;                       //distortion tree - Dy
TTree * treeDZ= 0;                       //distortion tree - Dz
// fit models ...
AliTPCFCVoltError3D *rotOFC  =0;          // fit models
AliTPCFCVoltError3D *rodOFC1 =0;
AliTPCFCVoltError3D *rodOFC2 =0;
AliTPCFCVoltError3D *rotIFC  =0;
AliTPCFCVoltError3D *rodIFC1 =0;
AliTPCFCVoltError3D *rodIFC2 =0;
AliTPCFCVoltError3D *rodAll  =0;        //all lookup initialized
//
AliTPCComposedCorrection *corFit0 = 0;  // composed correction
void FitFunctionQA();
void MakeQA();
void DrawAlignParam();

void PrintFit(TString fitString, TString filter){
  /// helper function to print substring of the fit
  /// TString fitString =*strFitGA;
  /// TString filter="rot";

  TObjArray *arr = fitString.Tokenize("++");
  Int_t entries=arr->GetEntries();
  for (Int_t i=0; i<entries; i++){
    TString aaa = arr->At(i)->GetName();
    if (aaa.Contains(filter)==1){
      printf("%d\t%s\n",i,aaa.Data());
    }    
  }
}


void LoadModels(){
  /// Load the models from the file
  /// Or create it

  Int_t volt = 1;
  TFile f("OCDB-RodShifts.root");
  AliCDBEntry *entry = (AliCDBEntry*) f.Get("AliCDBEntry");
  if (entry) { // load file    
    AliTPCComposedCorrection *cor = (AliTPCComposedCorrection*) entry->GetObject();    
    cor->Print();
    TCollection *iter = cor->GetCorrections();    
    rotOFC = (AliTPCFCVoltError3D*)iter->FindObject("rotOFC");
    rodOFC1 = (AliTPCFCVoltError3D*)iter->FindObject("rodOFC1");
    rodOFC2 = (AliTPCFCVoltError3D*)iter->FindObject("rodOFC2");
    rotIFC = (AliTPCFCVoltError3D*)iter->FindObject("rotIFC");
    rodIFC1 = (AliTPCFCVoltError3D*)iter->FindObject("rodIFC1");
    rodIFC2 = (AliTPCFCVoltError3D*)iter->FindObject("rodIFC2");
    rodAll = (AliTPCFCVoltError3D*)iter->FindObject("rodAll");
    // rotOFC->CreateHistoDZinXY(1,500,500)->Draw("surf2");
  } else {    
    // OFC 
    rotOFC = new AliTPCFCVoltError3D();
    rotOFC->SetOmegaTauT1T2(0,1,1);
    rotOFC->SetRotatedClipVoltA(1,volt);
    rotOFC->SetRotatedClipVoltC(1,volt);
    //
    rodOFC1 = new AliTPCFCVoltError3D();
    rodOFC1->SetOmegaTauT1T2(0,1,1);
    rodOFC1->SetRodVoltShiftA(18,volt);
    rodOFC1->SetRodVoltShiftC(18,volt);
    //
    rodOFC2 = new AliTPCFCVoltError3D();
    rodOFC2->SetOmegaTauT1T2(0,1,1);
    rodOFC2->SetCopperRodShiftA(18,volt);
    rodOFC2->SetCopperRodShiftC(18,volt);    
    // IFC     
    rotIFC = new AliTPCFCVoltError3D();
    rotIFC->SetOmegaTauT1T2(0,1,1);
    rotIFC->SetRotatedClipVoltA(0,volt);
    rotIFC->SetRotatedClipVoltC(0,volt);
    //
    rodIFC1 = new AliTPCFCVoltError3D();
    rodIFC1->SetOmegaTauT1T2(0,1,1);
    rodIFC1->SetRodVoltShiftA(0,volt);
    rodIFC1->SetRodVoltShiftC(0,volt);
    //
    rodIFC2 = new AliTPCFCVoltError3D();
    rodIFC2->SetOmegaTauT1T2(0,1,1);
    rodIFC2->SetCopperRodShiftA(0,volt);
    rodIFC2->SetCopperRodShiftC(0,volt);
    // dummy object with all inits
    rodAll = new AliTPCFCVoltError3D();
    Double_t volt0=0.0000000001;
    rodAll->SetRotatedClipVoltA(1,volt0);
    rodAll->SetRotatedClipVoltC(1,volt0);
    rodAll->SetRodVoltShiftA(18,volt0);
    rodAll->SetRodVoltShiftC(18,volt0);
    rodAll->SetCopperRodShiftA(18,volt0);
    rodAll->SetCopperRodShiftC(18,volt0);    
    rodAll->SetRotatedClipVoltA(0,volt0);
    rodAll->SetRotatedClipVoltC(0,volt0);
    rodAll->SetRodVoltShiftA(0,volt0);
    rodAll->SetRodVoltShiftC(0,volt0);
    rodAll->SetCopperRodShiftA(0,volt0);
    rodAll->SetCopperRodShiftC(0,volt0);
    rodAll->SetOmegaTauT1T2(0,1,1);
    //
    //
    // Initialization of the lookup tables
    //
    printf(" ------- OFC rotated clip:\n"); rotOFC->InitFCVoltError3D();
    printf(" ------- OFC rod & strip:\n");  rodOFC1->InitFCVoltError3D();
    printf(" ------- OFC copper rod:\n");   rodOFC2->InitFCVoltError3D();
    printf(" ------- IFC rotated clip:\n"); rotIFC->InitFCVoltError3D();
    printf(" ------- IFC rod & strip:\n");  rodIFC1->InitFCVoltError3D();
    printf(" ------- IFC copper rod:\n");   rodIFC2->InitFCVoltError3D();
    printf(" ------- Dummy all:\n");        rodAll->InitFCVoltError3D();
    // give names
    rotOFC->SetName("rotOFC");rotOFC->SetTitle("rotOFC");
    rodOFC1->SetName("rodOFC1");rodOFC1->SetTitle("rodOFC1");
    rodOFC2->SetName("rodOFC2");rodOFC2->SetTitle("rodOFC2");
    rotIFC->SetName("rotIFC");rotIFC->SetTitle("rotIFC");
    rodIFC1->SetName("rodIFC1");rodIFC1->SetTitle("rodIFC1");
    rodIFC2->SetName("rodIFC2");rodIFC2->SetTitle("rodIFC2");
    rodAll->SetName("rodAll");rodAll->SetTitle("rodAll");
    //
    // save in file
    AliTPCComposedCorrection *corFit0 = new AliTPCComposedCorrection();
    TObjArray *cs = new TObjArray();
    cs->Add(rotIFC); cs->Add(rotOFC);
    cs->Add(rodIFC1); cs->Add(rodOFC1);
    cs->Add(rodIFC2); cs->Add(rodOFC2);
    cs->Add(rodAll);
    corFit0->SetCorrections(cs);
    corFit0->SetOmegaTauT1T2(0,1.,1.);
    corFit0->Print();    
    corFit0->StoreInOCDB(0,1000000,"testForFit");    
  }
  //
  AliTPCCorrection::AddVisualCorrection(rotOFC,0); 
  AliTPCCorrection::AddVisualCorrection(rodOFC1,1); 
  AliTPCCorrection::AddVisualCorrection(rodOFC2,2); 
  AliTPCCorrection::AddVisualCorrection(rotIFC,3); 
  AliTPCCorrection::AddVisualCorrection(rodIFC1,4); 
  AliTPCCorrection::AddVisualCorrection(rodIFC2,5); 
  AliTPCCorrection::AddVisualCorrection(rodAll,6); 
}


void RegisterAlignFunction(){
  /// Register primitive alignment components.
  /// Linear conbination of primitev forulas used for fit
  /// The nominal delta 1 mm in shift and 1 mrad in rotation
  /// Primitive formulas registeren in AliTPCCoreection::AddvisualCorrection
  /// 10 - deltaX
  /// 11 - deltaY
  /// 12 - deltaZ
  /// 13 - rot0 (phi)
  /// 14 - rot1 (theta)
  /// 15 - rot2

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
  AliTPCCorrection::AddVisualCorrection(alignTrans0  ,10);
  AliTPCCorrection::AddVisualCorrection(alignTrans1  ,11);
  AliTPCCorrection::AddVisualCorrection(alignTrans2  ,12);
  AliTPCCorrection::AddVisualCorrection(alignRot0    ,13);
  AliTPCCorrection::AddVisualCorrection(alignRot1    ,14);
  AliTPCCorrection::AddVisualCorrection(alignRot2    ,15);
  TObjArray arrAlign(6);
  arrAlign.AddAt(alignTrans0->Clone(),0);
  arrAlign.AddAt(alignTrans1->Clone(),1);
  arrAlign.AddAt(alignTrans2->Clone(),2);
  arrAlign.AddAt(alignRot0->Clone(),3);
  arrAlign.AddAt(alignRot1->Clone(),4);
  arrAlign.AddAt(alignRot2->Clone(),5);
  //combAlign.SetCorrections((TObjArray*)arrAlign.Clone());
}


void LoadTrees(){
  /// 1. Load trees
  /// 2. Make standard cuts - filter the tree
  /// 3. makeStandard aliases

  TFile *fy = new TFile("clusterDY.root");
  TFile *fz = new TFile("clusterDZ.root");
  treeDY= (TTree*)fy->Get("delta");
  treeDZ= (TTree*)fz->Get("delta");
  TCut cutAll = "entries>3000&&abs(kZ)<1&&localX>80&&localX<246&&abs(sector-int(sector)-0.5)<0.41&&abs(localX-134)>2&&rmsG<0.5&&abs(localX-189)>3";

  treeDY->Draw(">>outListY",cutAll,"entryList");
  TEntryList *elistOutY = (TEntryList*)gDirectory->Get("outListY");
  treeDY->SetEntryList(elistOutY);
  treeDZ->Draw(">>outListZ",cutAll,"entryList");
  TEntryList *elistOutZ = (TEntryList*)gDirectory->Get("outListZ");
  treeDZ->SetEntryList(elistOutZ);
  //
  // Make aliases
  //
  treeDY->SetAlias("dsec","(sector-int(sector)-0.5)");
  treeDY->SetAlias("dsec0","(sector-int(sector))");
  treeDY->SetAlias("signy","(-1.+2*(sector-int(sector)>0.5))");
  treeDY->SetAlias("dx","(localX-165.5)");   // distance X to ref. plane
  treeDY->SetAlias("dxm","0.001*(localX-165.5)");
  treeDY->SetAlias("rx","(localX/166.5)");   //ratio distrnace to reference plane
  treeDY->SetAlias("dq1","(((q1==0)*(-rx)+q1*(1-rx))*signy)");
  //
  {for (Int_t isec=0; isec<18; isec++){ // sectors
      treeDY->SetAlias(Form("sec%d",isec), Form("(abs(sector-%3.1lf)<0.5)",isec+0.5));
      treeDZ->SetAlias(Form("sec%d",isec), Form("(abs(sector-%3.1lf)<0.5)",isec+0.5));
    }}
  treeDY->SetAlias("gy","localX*sin(pi*sector/9)");
  treeDY->SetAlias("gx","localX*cos(pi*sector/9)");
  treeDY->SetAlias("side","(-1.+2.*(kZ>0))");  
  treeDY->SetAlias("drphi","mean");
  treeDY->SetMarkerStyle(25);
}

Double_t DeltaLookup(Double_t sector, Double_t localX, Double_t kZ, Double_t xref, Int_t value, Int_t corr){
  /// Distortion maps are calculated relative to the reference X plane
  /// The same procedure applied for given correction
  /// fd(sector, localX, kZ) ==>  fd(sector, localX, kZ)-fd(sector, xref, kZ)*localX/xref

  Double_t distortion    = AliTPCCorrection::GetCorrSector(sector,localX,kZ,value,corr);
  Double_t distortionRef = AliTPCCorrection::GetCorrSector(sector,xref,kZ,value,corr)*localX/xref;
  return distortion-distortionRef;
}

void MakeAliases(){
  /// make alias names

  AliTPCROC * roc = AliTPCROC::Instance();
  Double_t xref  = ( roc->GetPadRowRadii(0,0)+roc->GetPadRowRadii(36,roc->GetNRows(36)-1))*0.5;
  treeDY->SetAlias("iroc","(localX<134)");  // IROC
  treeDY->SetAlias("oroc","(localX>134)");  // OROC
  treeDY->SetAlias("q1","abs(localX-161)<28"); // OROC first haf
  treeDY->SetAlias("q2","(localX>189)");       // OROC second half
  //
  treeDY->SetAlias("dIFC","abs(localX-83)");    // distance to IFC
  treeDY->SetAlias("dOFC","abs(localX-250)");   // distance to OFC
  treeDY->SetAlias("errY","(1.+1/(1+(dIFC/2.))+1/(1+dOFC/2.))");

  //
  // rotated clip - IFC --------------------
  treeDY->SetAlias("rotClipI","DeltaLookup(sector,localX,kZ,165.6,1,3+0)");
  // rotated clip - OFC --------------------
  treeDY->SetAlias("rotClipO","DeltaLookup(sector,localX,kZ,165.6,1,0+0)");
  for (Int_t isec=0;isec<18; isec++){    
    //
    // alignment IROC - shift cm - rotation mrad
    treeDY->SetAlias(Form("dyI_%d",isec),Form("(iroc*sec%d)",isec));
    treeDY->SetAlias(Form("drotI_%d",isec),Form("(iroc*sec%d*(localX-%3.1lf)/1000.)",isec,xref));
    // alignment OROC
    treeDY->SetAlias(Form("dyO_%d",isec),Form("(sec%d*oroc-sec%d*localX/165.6)",isec,isec)); // fix local Y shift - do not use it
    treeDY->SetAlias(Form("drotO_%d",isec),Form("(sec%d*oroc*(localX-%3.1lf)/1000.)",isec,xref));
    //
    // quadrant shift OROC 
    treeDY->SetAlias(Form("dqO0_%d",isec),Form("(sec%d*sign(dsec)*(q1-localX/165.6))",isec));
    // quadrant delta rotation OFC inner
    treeDY->SetAlias(Form("dqOR0_%d",isec),Form("(sec%d*sign(dsec)*q1*(localX-165.6)/1000.)",isec));
    //
    treeDY->SetAlias(Form("dqO1_%d",isec),Form("(q2*sec%d*sign(dsec))",isec));  // delta
    treeDY->SetAlias(Form("dqOR1_%d",isec),Form("(q2*sec%d*sign(dsec)*(localX-165.6)/1000.)",isec));
    treeDY->SetAlias(Form("dqO2_%d",isec),Form("(q2*sec%d)",isec)); //common
    treeDY->SetAlias(Form("dqOR2_%d",isec),Form("(q2*sec%d*(localX-165.6)/1000.)",isec));
    //
    //
    // shifted rod -  OFC
    treeDY->SetAlias(Form("rodStripO_%d",isec),Form("DeltaLookup(sector-%d,localX,kZ,165.6,1,1+0)",isec));
    // Shifted cooper rod OFC 
    treeDY->SetAlias(Form("coppRodO_%d",isec),Form("DeltaLookup(sector-%d,localX,kZ,165.6,1,2+0)",isec)); 
    // shifted rod - IFC 
    treeDY->SetAlias(Form("rodStripI_%d",isec),Form("DeltaLookup(sector-%d,localX,kZ,165.6,1,4+0)",isec)); 
    // Shifted cooper rod IFC 
    treeDY->SetAlias(Form("coppRodI_%d",isec),Form("DeltaLookup(sector-%d,localX,kZ,165.6,1,5+0)",isec)); 
  }
  //
  // fitted correction - Using the saved coposed corrections
  //
  treeDY->SetAlias("dAlign","DeltaLookup(sector,localX,kZ,165.6,1,1000+0)");
  treeDY->SetAlias("dQuadrant","DeltaLookup(sector,localX,kZ,165.6,1,1100+0)");
  treeDY->SetAlias("dField","DeltaLookup(sector,localX,kZ,165.6,1,100+0)");
  treeDY->SetAlias("dAll","(dAlign+dQuadrant+dField)");
}

AliTPCCalibGlobalMisalignment *MakeAlignCorrection(TVectorD paramA, TVectorD paramC,  TMatrixD covar, Double_t chi2){
  /// Make a global alignmnet
  /// Take a fit parameters and make a combined correction:
  /// Only delta local Y and delta phi are fitted - not sensitivity for other parameters
  /// GX and GY shift extracted per side.
  ///
  /// Algorithm:
  ///   1. Loop over sectors
  ///   2. Make combined AliTPCCalibGlobalMisalignment

  AliTPCCalibGlobalMisalignment *alignLocal  =new  AliTPCCalibGlobalMisalignment;  
  Int_t offset=3;
  AliTPCROC * roc = AliTPCROC::Instance();
  Double_t xref  = ( roc->GetPadRowRadii(0,0)+roc->GetPadRowRadii(36,roc->GetNRows(36)-1))*0.5;
  // 
  pcWorkspace->GetFile()->cd();
  TObjArray * array = new TObjArray(72);
  
  for (Int_t side=-1; side<2; side+=2){
    TVectorD &param = (side==1) ? paramA : paramC;
    TGeoHMatrix matrixGX;
    TGeoHMatrix matrixGY;
    matrixGX.SetDx(0.1*param[1]); 
    matrixGY.SetDy(0.1*param[2]); 
    //  
    for (Int_t isec=0; isec<18; isec++){ 
      TGeoHMatrix matrixOROC;
      TGeoHMatrix matrixIROC; 
      TGeoRotation rotIROC;
      TGeoRotation rotOROC;
      Double_t phi= (Double_t(isec)+0.5)*TMath::Pi()/9.;
      //
      Double_t drotIROC = -param[offset+isec+18]*0.001;
      Double_t drotOROC = -param[offset+isec+36]*0.001;
      Double_t dlyIROC  = param[offset+isec];
      Double_t dlyIROC0 = param[offset+isec]+drotIROC*xref;     // shift at x=0
      Double_t dlyOROC  = 0;
      Double_t dlyOROC0 = 0+drotOROC*xref;  // shift at x=0
      //
      Double_t dgxIROC = TMath::Cos(phi)*0+TMath::Sin(phi)*dlyIROC0;
      Double_t dgyIROC = TMath::Sin(phi)*0-TMath::Cos(phi)*dlyIROC0;
      Double_t dgxOROC = TMath::Cos(phi)*0+TMath::Sin(phi)*dlyOROC0;
      Double_t dgyOROC = TMath::Sin(phi)*0-TMath::Cos(phi)*dlyOROC0;    
      Double_t errYIROC=TMath::Sqrt(covar(offset+isec, offset+isec)*chi2);
      Double_t errPhiIROC=TMath::Sqrt(covar(offset+isec+18, offset+isec+18)*chi2);
      Double_t errPhiOROC=TMath::Sqrt(covar(offset+isec+36, offset+isec+36)*chi2);
      matrixIROC.SetDx(dgxIROC);
      matrixIROC.SetDy(dgyIROC);
      matrixOROC.SetDx(dgxOROC);
      matrixOROC.SetDy(dgyOROC);
      rotIROC.SetAngles(drotIROC*TMath::RadToDeg(),0,0);
      matrixIROC.Multiply(&matrixGX);
      matrixIROC.Multiply(&matrixGY);
      matrixIROC.Multiply(&rotIROC);
      rotOROC.SetAngles(drotOROC*TMath::RadToDeg(),0,0);
      matrixOROC.Multiply(&matrixGX);
      matrixOROC.Multiply(&matrixGY);
      matrixOROC.Multiply(&rotOROC);
      if (side>0){
	array->AddAt(matrixIROC.Clone(),isec);
	array->AddAt(matrixOROC.Clone(),isec+36);
      }
      if (side<0){
	array->AddAt(matrixIROC.Clone(),isec+18);
	array->AddAt(matrixOROC.Clone(),isec+36+18);
      }
      Double_t rms=TMath::Sqrt(chi2); //chi2 normalized 1/npoints before
      (*pcWorkspace)<<"align"<<
	"isec="<<isec<<           // sector
	"side="<<side<<           // side
	"phi="<<phi<<             // phi
	"rms="<<rms<<             // rms in cm
	"cov.="<<&covar<<
	// errors
	"errYIROC="<<errYIROC<<      //error in local y
	"errPhiIROC="<<errPhiIROC<<  //error in phi IROC
	"errPhiOROC="<<errPhiOROC<<  //error in phi OROC
	//
	"dlyIROC="<<dlyIROC<<        //delta Local Y at refX=165 -IROC
	"dlyIROC0="<<dlyIROC0<<      //delta Local Y at refX=0   -IROC
	"dgxIROC="<<dgxIROC<<        //
	"dgyIROC="<<dgyIROC<<
	"drotIROC="<<drotIROC<<
	//
	"dlyOROC="<<dlyOROC<<       // assuming 0 at 
	"dlyOROC0="<<dlyOROC0<<
	"dgxOROC="<<dgxOROC<<
	"dgyOROC="<<dgyOROC<<
	"drotOROC="<<drotOROC<<
	"\n";
    } 
  }
  alignLocal->SetAlignSectors(array);
  AliTPCCorrection::AddVisualCorrection(alignLocal,1000);
  /*
    Check:
    alignLocal->CreateHistoDRPhiinXY(10,100,100)->Draw("colz");  // OK
    treeDY->Draw("DeltaLookup(sector,localX,kZ,165.6,1,1001):localX:int(sector)","kZ>0","colz");
    //
    (( (*pcWorkspace)<<"align"))->GetTree()->Draw("dlyOROC-dlyOROC:isec","","*")
    treeDY->Draw("tfitA:DeltaLookup(sector,localX,kZ,165.6,1,1001):sector","kZ>0","colz");

  */
  return alignLocal;
}

AliTPCCalibGlobalMisalignment *MakeQuadrantCorrection(TVectorD paramA, TVectorD paramC, TMatrixD covar, Double_t chi2){
  /// Make a global alignmnet
  /// side= 1 - A side
  /// side=-1 - C side
  /// Take a fit parameters and make a combined correction:
  /// Only delta local Y and delta phi are fitted - not sensitivity for other parameters
  /// GX and GY shift extracted per side.
  ///
  /// Algorithm:
  ///   1. Loop over sectors
  ///   2. Make combined AliTPCCalibGlobalMisalignment

  AliTPCCalibGlobalMisalignment *alignLocalQuadrant  =new  AliTPCCalibGlobalMisalignment;  
  //
  Int_t offset=3+3*18;
  TVectorD quadrantQ0(36);   //dly+=sign*(*fQuadrantQ0)[isec];  // shift in cm
  TVectorD quadrantRQ0(36);  //dly+=sign*(*fQuadrantRQ0)[isec]*(pos[0]-xref);      
  TVectorD quadrantQ1(36);   //dly+=sign*(*fQuadrantQ1)[isec];  // shift in cm
  TVectorD quadrantRQ1(36);  //dly+=sign*(*fQuadrantRQ1)[isec]*(pos[0]-xref);      
  TVectorD quadrantQ2(36);   //dly+=(*fQuadrantQ2)[isec];  // shift in cm
  TVectorD quadrantRQ2(36);  //dly+=(*fQuadrantRQ2)[isec]*(pos[0]-xref);      

  for (Int_t side=-1; side<2; side+=2){
    Int_t offsetSec= (side==1) ? 0:18;
    TVectorD &param = (side==1) ? paramA : paramC;
    for (Int_t isec=0; isec<18; isec++){ 
      Double_t  q0=-param[offset+isec+0*18];
      Double_t rq0=-param[offset+isec+1*18]*0.001;
      Double_t  q1=-param[offset+isec+2*18];
      Double_t rq1=-param[offset+isec+3*18]*0.001;
      Double_t  q2=-param[offset+isec+4*18];
      Double_t rq2=-param[offset+isec+5*18]*0.001;
      //
      Double_t  sq0=TMath::Sqrt(covar(offset+isec+0*18,offset+isec+0*18)*chi2);
      Double_t srq0=TMath::Sqrt(covar(offset+isec+1*18,offset+isec+1*18)*chi2)*0.001;
      Double_t  sq1=TMath::Sqrt(covar(offset+isec+2*18,offset+isec+2*18)*chi2);
      Double_t srq1=TMath::Sqrt(covar(offset+isec+3*18,offset+isec+3*18)*chi2)*0.001;
      Double_t  sq2=TMath::Sqrt(covar(offset+isec+4*18,offset+isec+4*18)*chi2);
      Double_t srq2=TMath::Sqrt(covar(offset+isec+5*18,offset+isec+5*18)*chi2)*0.001;
      //
      quadrantQ0[offsetSec+isec]=q0;
      quadrantRQ0[offsetSec+isec]=rq0;
      quadrantQ1[offsetSec+isec]=q1;
      quadrantRQ1[offsetSec+isec]=rq1;
      quadrantQ2[offsetSec+isec]=q2;
      quadrantRQ2[offsetSec+isec]=rq2;
      Double_t rms=TMath::Sqrt(chi2); //chi2 normalized 1/npoints before
      (*pcWorkspace)<<"quadrant"<<
	"isec="<<isec<<           // sector
	"side="<<side<<           // side
	"rms="<<rms<<             // rms in cm
	"cov.="<<&covar<<
	//
	"q0="<<q0<<               // quadrant alignment
	"rq0="<<rq0<<
	"q1="<<q1<<
	"rq1="<<rq1<<
	"q2="<<q2<<
	"rq2="<<rq2<<
	"sq0="<<sq0<<             //rms of quadrant parameters 
	"srq0="<<srq0<<
	"sq1="<<sq1<<
	"srq1="<<srq1<<
	"sq2="<<sq2<<
	"srq2="<<srq2<<
	"\n";    
    }
  }
  alignLocalQuadrant->SetQuadranAlign(&quadrantQ0, &quadrantRQ0, &quadrantQ1, &quadrantRQ1, &quadrantQ2, &quadrantRQ2);
  AliTPCCorrection::AddVisualCorrection(alignLocalQuadrant,1100);  
  /*
    (( (*pcWorkspace)<<"quadrant"))->GetTree()->Draw("q0:isec","","*")

    alignLocalQuadrant->CreateHistoDRPhiinXY(10,100,100)->Draw("colz");  //OK
    
    treeDY->Draw("tfitA-DeltaLookup(sector,localX,kZ,165.6,1,1000):DeltaLookup(sector,localX,kZ,165.6,1,1100):int(sector)","kZ>0&&oroc","colz");
  */
  return alignLocalQuadrant;
}

AliTPCFCVoltError3D* MakeEfieldCorrection(TVectorD paramA, TVectorD paramC, TMatrixD covar, Double_t chi2){
  /// Make a global  AliTPCFCVoltError3D object
  ///
  /// Take a fit parameters and make a combined correction:
  /// Only delta local Y and delta phi are fitted - not sensitivity for other parameters
  /// GX and GY shift extracted per side.
  ///
  /// Algorithm:
  ///   1. Loop over sectors
  ///   2. Make combined AliTPCCalibGlobalMisalignment

  Int_t offset=3+9*18;
  //
  AliTPCFCVoltError3D* corrField = new AliTPCFCVoltError3D;
  //
  Double_t rotIROCA=paramA[offset+0];
  Double_t rotIROCC=paramC[offset+0];
  Double_t rotOROCA=paramA[offset+1];
  Double_t rotOROCC=paramC[offset+1];
  //
  corrField->SetRotatedClipVoltA(0,paramA[offset+0]);
  corrField->SetRotatedClipVoltC(0,paramC[offset+0]);
  corrField->SetRotatedClipVoltA(1,paramA[offset+1]);
  corrField->SetRotatedClipVoltC(1,paramC[offset+1]);
  {for (Int_t isec=0; isec<18; isec++){
      corrField->SetRodVoltShiftA(isec,paramA[offset+2+36+isec]);     // rod shift IFC
      corrField->SetRodVoltShiftA(18+isec,paramA[offset+2+isec]);     // rod shift OFC
      corrField->SetRodVoltShiftC(isec,paramC[offset+2+36+isec]);    // rod shift IFC
      corrField->SetRodVoltShiftC(18+isec,paramC[offset+2+isec]);    // rod shift OFC
      Double_t rodIROCA=paramA[offset+2+36+isec];
      Double_t rodIROCC=paramC[offset+2+36+isec];
      Double_t rodOROCA=paramA[offset+2+isec];
      Double_t rodOROCC=paramC[offset+2+isec];
      //
      Double_t srodIROCA=TMath::Sqrt(covar(offset+2+36+isec,offset+2+36+isec)*chi2);
      Double_t srodOROCA=TMath::Sqrt(covar(offset+2+isec,offset+2+isec)*chi2);
      Double_t phi= (Double_t(isec)+0.5)*TMath::Pi()/9.;
      Double_t rms=TMath::Sqrt(chi2); //chi2 normalized 1/npoints before
      (*pcWorkspace)<<"field"<<
	"isec="<<isec<<           // sector
	//"side="<<side<<           // side
	"phi="<<phi<<             // phi
	"rms="<<rms<<             // rms in cm
	"cov.="<<&covar<<
	//
	"rotIROCA="<<rotIROCA<<  
	"rotIROCC="<<rotIROCC<<
	"rotOROCA="<<rotOROCA<<
	"rotOROCC="<<rotOROCC<<
	//
	"rodIROCA="<<rodIROCA<<
	"rodIROCC="<<rodIROCC<<
	"rodOROCA="<<rodOROCA<<
	"rodOROCC="<<rodOROCC<<
	//
	"srodIROCA="<<srodIROCA<<
	"srodOROCA="<<srodOROCA<<
	"srodIROCC="<<srodIROCA<<
	"srodOROCC="<<srodOROCA<<
	"\n";
    }    
  }
  corrField->SetOmegaTauT1T2(0,1.,1.);
  corrField->InitFCVoltError3D();
  AliTPCCorrection::AddVisualCorrection(corrField,100);  
  return corrField;
}


void FitRodShift(Bool_t flagIFCcopper = kTRUE) {
  /// Main fit function
  /// In total 440 parameters to fit using global fit:
  ///   1. Rotation and translation for each sector (72x2)
  ///   2. Rotation and translation for each quadrant of OROC (36x4x2)
  ///   3. Rod/strip shifts in IFC and OFC (18 sectors x 2 sides x 2 )
  ///   4. Rotated clips in IFC and OFC

  LoadTrees();
  LoadModels();
  RegisterAlignFunction();
  MakeAliases();

  pcWorkspace= new TTreeSRedirector("fitAlignLookup.root");
  TFormula::SetMaxima(10000);
  //
  // 1. fit Global
  //
  Double_t chi2G=0;   Int_t npointsG=0;  TVectorD paramG;  TMatrixD covarG;
  TString fstringGlobal="";  
  fstringGlobal+="DeltaLookup(sector,localX,kZ,165.6,1,10-0)++"; // deltaGX
  fstringGlobal+="DeltaLookup(sector,localX,kZ,165.6,1,11-0)++"; // deltaGY
  fstringGlobal+="DeltaLookup(sector,localX,kZ,165.6,1,10-0)*side++"; // deltaGX - sides
  fstringGlobal+="DeltaLookup(sector,localX,kZ,165.6,1,11-0)*side++"; // deltaGY -sides
  //
  TString *strFitGlobal = TStatToolkit::FitPlane(treeDY,"drphi", fstringGlobal.Data(),"1", chi2G,npointsG,paramG,covarG,-1,0, 10000000, kTRUE);
  treeDY->SetAlias("fitYGlobal",strFitGlobal->Data());
  strFitGlobal->Tokenize("++")->Print();
  printf("chi2=%f\n",TMath::Sqrt(chi2G/npointsG));
  // testplots - if necessary - e.g.
  // treeDY->Draw("AliTPCCorrection::GetCorrSector(sector,localX,kZ,1,2):sector:localX","Cut","colz")

  // FIT PREPARATION +++++++++++++++++++++++++++++++++
  // define the fit function

  TString  fstring="";         

  //
  // Alignment
  //
  {
    fstring+="DeltaLookup(sector,localX,kZ,165.6,1,10-0)++"; // deltaGX
    fstring+="DeltaLookup(sector,localX,kZ,165.6,1,11-0)++"; // deltaGY
    // alignment IROC
    for (Int_t i=0;i<18;i++){
      fstring+=Form("dyI_%d++",i);  // alignment - linear shift (in cm)
    }
    for (Int_t i=0;i<18;i++){
      fstring+=Form("drotI_%d++",i);  // alignment - rotation (in mrad)
    }
    // alignment OROC
    for (Int_t i=0;i<18;i++){ // shift of OROC is not allowed
      //fstring+=Form("dyO_%d++",i);  // alignment - linear shift (in cm)
    }
    for (Int_t i=0;i<18;i++){
      fstring+=Form("drotO_%d++",i);  // alignment - rotation (in mrad)
    }
    //
    for (Int_t i=0;i<18;i++){
      fstring+=Form("dqO0_%d++",i);  // alignment - quadrant shift OROC in
    }
    for (Int_t i=0;i<18;i++){
      fstring+=Form("dqOR0_%d++",i);  // alignment - quadrant rotation OROC in
    }
    for (Int_t i=0;i<18;i++){
      fstring+=Form("dqO1_%d++",i);  // alignment - quadrant shift OROC out
    }
    for (Int_t i=0;i<18;i++){
      fstring+=Form("dqOR1_%d++",i);  // alignment - quadrant rotation OROC out
    }
    for (Int_t i=0;i<18;i++){
      fstring+=Form("dqO2_%d++",i);  // alignment - quadrant shift OROC out
    }
    for (Int_t i=0;i<18;i++){
      fstring+=Form("dqOR2_%d++",i);  // alignment - quadrant rotation OROC out
    }
  }
  //
  // Field distortion
  //
  {
    // rotated clip - IFC --------------------
    treeDY->SetAlias("rotClipI","DeltaLookup(sector,localX,kZ,165.6,1,3+0)");
    fstring+="rotClipI++";
    // rotated clip - OFC --------------------
    treeDY->SetAlias("rotClipO","DeltaLookup(sector,localX,kZ,165.6,1,0+0)");
    fstring+="rotClipO++";
    
    // shifted rod -  OFC
    for (Int_t i=0;i<18;i++){
      fstring+=Form("rodStripO_%d++",i);
    }    
    // Shifted cooper rod OFC     
    for (Int_t i=0;i<18;i++){
      fstring+=Form("coppRodO_%d++",i);
    }    
    // shifted rod - IFC 
    for (Int_t i=0;i<18;i++){
      fstring+=Form("rodStripI_%d++",i);
    }    
    // Shifted cooper rod IFC 
    if (flagIFCcopper) for (Int_t i=0;i<18;i++){
	fstring+=Form("coppRodI_%d++",i);
    }  
  }
  //fstring+=fstringGlobal;
 
  // FIT ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  Double_t chi2A=0;   Int_t npointsA=0;  TVectorD paramA;  TMatrixD covarA;
  Double_t chi2C=0;   Int_t npointsC=0;  TVectorD paramC;  TMatrixD covarC;
  
  
  printf("Fitting A side\n");
  TCut cutAllA = "Cut&&kZ>0";
  TString *strFitGA = TStatToolkit::FitPlane(treeDY,"drphi:errY", fstring.Data(),"kZ>0", chi2A,npointsA,paramA,covarA,-1,0, 10000000, kTRUE);
  treeDY->SetAlias("tfitA",strFitGA->Data());
  // strFitGA->Tokenize("++")->Print();
  printf("chi2=%f\n",TMath::Sqrt(chi2A/npointsA));
  printf("Sigma Y:\t%f (mm)\n",10.*TMath::Sqrt(covarA(3,3)*chi2A/npointsA));
  printf("IROC Sigma Angle:\t%f (mrad)\n",TMath::Sqrt(covarA(3+18,3+18)*chi2A/npointsA));
  printf("OROC Sigma Angle:\t%f (mrad)\n",TMath::Sqrt(covarA(3+36,3+36)*chi2A/npointsA));

  printf("Fitting C side\n");
  TCut cutAllC = "Cut&&kZ<0";
  TString *strFitGC = TStatToolkit::FitPlane(treeDY,"drphi:errY", fstring.Data(),"kZ<0", chi2C,npointsC,paramC,covarC,-1,0, 10000000, kTRUE);
  treeDY->SetAlias("tfitC",strFitGC->Data());
  //  strFitGC->Tokenize("++")->Print();
  printf("chi2=%f\n",TMath::Sqrt(chi2C/npointsC));
  printf("Sigma Y:\t%f (mm)\n",10.*TMath::Sqrt(covarC(3,3)*chi2C/npointsC));
  printf("IROC Sigma Angle:\t%f (mrad)\n",TMath::Sqrt(covarC(3+18,3+18)*chi2C/npointsC));
  printf("OROC Sigma Angle:\t%f (mrad)\n",TMath::Sqrt(covarC(3+36,3+36)*chi2C/npointsC));

  //
  // make combined correction
  //
  TH2 *hdist =0;
  //
  AliTPCCalibGlobalMisalignment * alignLocal = MakeAlignCorrection(paramA, paramC, covarA, chi2A/npointsA);
  alignLocal->SetName("alignLocal");
  alignLocal->Write("alignLocal");
  hdist = alignLocal->CreateHistoDRPhiinXY(10,250,250);
  hdist->SetName("AlignLocalAside"); hdist->SetTitle("Alignment map (A side)");
  hdist->Write("AlignLocalAside");
  hdist = alignLocal->CreateHistoDRPhiinXY(-10,250,250);
  hdist->SetName("AlignLocalCside"); hdist->SetTitle("Alignment map (C side)");
  hdist->Write("AlignLocalCside");
  //
  AliTPCCalibGlobalMisalignment *alignQuadrant = MakeQuadrantCorrection(paramA, paramC,  covarA, chi2A/npointsA);
  alignQuadrant->SetName("alignQuadrant");
  alignQuadrant->Write("alignQuadrant");
  hdist = alignQuadrant->CreateHistoDRPhiinXY(10,250,250);
  hdist->SetName("AlignQuadrantAside"); hdist->SetTitle("Quadrant Alignment map (A side)");
  hdist->Write("AlignQuadrantAside");
  hdist = alignQuadrant->CreateHistoDRPhiinXY(-10,250,250);
  hdist->SetName("AlignQuadrantCside"); hdist->SetTitle("Quadrant Alignment map (C side)");
  hdist->Write("AlignQuadrantCside");
  //
  AliTPCFCVoltError3D* corrField = MakeEfieldCorrection(paramA, paramC, covarA, chi2A/npointsA);
  corrField->SetName("corrField");
  corrField->Write("corrField");

  hdist = corrField->CreateHistoDRPhiinXY(10,250,250);
  hdist->SetName("AlignEfieldAside"); hdist->SetTitle("Efield Alignment map (A side)");
  hdist->Write("AlignEfieldAside");
  hdist = corrField->CreateHistoDRPhiinXY(-10,250,250);
  hdist->SetName("AlignEfieldCside"); hdist->SetTitle("Efield Alignment map (C side)");
  hdist->Write("AlignEfieldCside");


  //
  delete pcWorkspace;
  MakeQA();
  return;  
}


void MakeQA(){
  LoadTrees();
  LoadModels();
  RegisterAlignFunction();
  MakeAliases();
  TFile f("fitAlignLookup.root","update");  
  AliTPCCalibGlobalMisalignment*     alignLocal =  (AliTPCCalibGlobalMisalignment*)f.Get("alignLocal");
  AliTPCCalibGlobalMisalignment*  alignQuadrant =  (AliTPCCalibGlobalMisalignment*)f.Get("alignQuadrant");
  AliTPCFCVoltError3D      *corrField= (AliTPCFCVoltError3D*)f.Get("corrField");
  AliTPCCorrection::AddVisualCorrection(alignLocal,1000);
  AliTPCCorrection::AddVisualCorrection(alignQuadrant,1100);
  AliTPCCorrection::AddVisualCorrection(corrField,100);
  FitFunctionQA();
  DrawAlignParam();
  f.Close();
}

void FitFunctionQA(){
  ///

  TH1 *his=0;
  TCanvas *canvasDist= new TCanvas("FitQA","fitQA",1200,800);
  canvasDist->Divide(2,2);
  //
  canvasDist->cd(1)->SetLogy(kFALSE); canvasDist->cd(1)->SetRightMargin(0.15);
  treeDY->Draw("10*mean:sector:localX","kZ<0&&localX<134","colz"); 
  his= (TH1*)treeDY->GetHistogram()->Clone();   his->SetName("DeltaRPhi1");
  his->GetXaxis()->SetTitle("Sector"); 
  his->GetZaxis()->SetTitle("R (cm)"); 
  his->GetYaxis()->SetTitle("#Delta_{r#phi} (mm)");
  his->Draw("cont2"); 
  his->Draw("colz");
  //
  canvasDist->cd(2)->SetLogy(kFALSE); canvasDist->cd(2)->SetRightMargin(0.15);
  treeDY->Draw("10*mean:sector:localX","kZ<0&&localX>160","colz"); 
  his= (TH1*)treeDY->GetHistogram()->Clone();   his->SetName("DeltaRPhi2");
  his->GetXaxis()->SetTitle("Sector"); 
  his->GetZaxis()->SetTitle("R (cm)"); 
  his->GetYaxis()->SetTitle("#Delta_{r#phi} (mm)");
  his->Draw("cont2"); 
  his->Draw("colz");
  //
  canvasDist->cd(3)->SetLogy(kFALSE); canvasDist->cd(3)->SetRightMargin(0.15);
  treeDY->Draw("10*(mean-dAll):sector:localX","kZ<0&&localX<134","colz"); 
  his= (TH1*)treeDY->GetHistogram()->Clone();   his->SetName("DeltaRPhi3"); his->SetTitle("Delta #RPhi");
  his->GetXaxis()->SetTitle("Sector"); 
  his->GetZaxis()->SetTitle("R (cm)"); 
  his->GetYaxis()->SetTitle("#Delta_{r#phifit}-#Delta_{r#phi} (mm)"); 
  his->Draw("cont2"); 
  his->Draw("colz");
  //
  canvasDist->cd(4)->SetLogy(kFALSE); canvasDist->cd(4)->SetRightMargin(0.15);
  treeDY->SetMarkerColor(1);
  treeDY->Draw("10*mean:10*dAll:localX","","colz"); 
  his= (TH1*)treeDY->GetHistogram()->Clone();   his->SetName("DeltaRPhi4");
  his->GetXaxis()->SetTitle("Fit value #Delta_{r#phi} (mm)"); 
  his->GetYaxis()->SetTitle("#Delta_{r#phi} (mm)"); 
  his->GetZaxis()->SetTitle("R (cm)"); 
  his->Draw("cont2");
  his->Draw("colz");
  //
  canvasDist->Write("FitQA");
}

void DrawAlignParam(){
  ///

  TFile f("fitAlignLookup.root");
  TTree * treeAlign=(TTree*)f.Get("align");
  TTree * treeQuadrant=(TTree*)f.Get("quadrant");
  TCanvas *canvasAlign=new TCanvas("align","align",1000,800); 
  TH1 *his=0;
  canvasAlign->Divide(2,2);
  treeAlign->SetMarkerStyle(25);
  treeQuadrant->SetMarkerStyle(25);
  gStyle->SetOptStat(kTRUE);
  //
  canvasAlign->cd(1); canvasAlign->cd(1)->SetRightMargin(0.15);
  treeAlign->Draw("10*dlyIROC*side:isec:1+side","","colz");
  his= (TH1*)treeAlign->GetHistogram()->Clone();   his->SetName("dlyIROC");his->SetTitle("IROC Alignment #Delta_{r#phi}");
  his->GetXaxis()->SetTitle("sector"); 
  his->GetYaxis()->SetTitle("#Delta_{r#phi} (mm)"); 
  his->GetZaxis()->SetTitle("side");
  his->Draw("cont2"); 
  his->Draw("colz");
  //
  canvasAlign->cd(2); canvasAlign->cd(2)->SetRightMargin(0.15);
  treeAlign->Draw("1000*drotIROC*side:isec:1+side","","colz");
  his= (TH1*)treeAlign->GetHistogram()->Clone();   his->SetName("drotIROC");his->SetTitle("IROC Angular Alignment #Delta_{r#phi}");
  his->GetXaxis()->SetTitle("sector"); 
  his->GetYaxis()->SetTitle("#Delta_{r#phi} (mrad)"); 
  his->GetZaxis()->SetTitle("side");
  his->Draw("cont2");
  his->Draw("colz");
  //
  canvasAlign->cd(4);canvasAlign->cd(4)->SetRightMargin(0.15);
  treeAlign->Draw("1000*drotOROC*side:isec:1+side","","colz");
  his= (TH1*)treeAlign->GetHistogram()->Clone();   his->SetName("drotOROC");his->SetTitle("OROC Angular Alignment #Delta_{r#phi}");
  his->GetXaxis()->SetTitle("sector"); 
  his->GetYaxis()->SetTitle("#Delta_{r#phi} (mrad)"); 
  his->GetZaxis()->SetTitle("side");
  his->Draw("cont2");
  his->Draw("colz");
  //
  canvasAlign->cd(3);canvasAlign->cd(3)->SetRightMargin(0.15);
  treeQuadrant->Draw("10*q2*side:isec:1+side","","colz");
  his= (TH1*)treeQuadrant->GetHistogram()->Clone();   his->SetName("drphiOROC");his->SetTitle("OROC Alignment  Outer Quadrant #Delta_{r#phi}");
  his->GetXaxis()->SetTitle("sector"); 
  his->GetYaxis()->SetTitle("#Delta_{r#phi} (mm)"); 
  his->GetZaxis()->SetTitle("side");
  his->Draw("cont2");
  his->Draw("colz");


}
