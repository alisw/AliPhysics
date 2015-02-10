/// \file CalibPID.C
///
/// 1. dump information to the tree
///
/// ~~~{.cpp}
/// gSystem->Load("libANALYSIS");
/// gSystem->Load("libTPCcalib");
/// gSystem->Load("libSTAT");
///
/// .L $ALICE_ROOT/TPC/CalibMacros/CalibPID.C+
/// .x ../ConfigOCDB.C
/// paramCl = AliTPCcalibDB::Instance()->GetClusterParam();
/// Init("calibPID06");
/// LookupHisto() // change SetRange in LookupHisto if needed !, check with pid->GetHistQtot()->Projection(0,1)->Draw("colz")
/// ~~~
///
/// exit aliroot
///
/// 2. update the OCDB
///
/// ~~~{.cpp}
/// gSystem->Load("libANALYSIS");
/// gSystem->Load("libTPCcalib");
/// gSystem->Load("libSTAT");
///
/// .L $ALICE_ROOT/TPC/CalibMacros/CalibPID.C+
/// .x ../ConfigOCDB.C
/// paramCl = AliTPCcalibDB::Instance()->GetClusterParam();
///
/// TFile fff("lookupdEdx.root")
/// TTree * treeDump =0;
/// TObjArray fitArr;
/// treeDump = (TTree*)fff.Get("dumpdEdx");
///
/// TCut cutAll = "meangTot>0.0&&sumMax>150&&sumTot>150&&rmsgMax/rmsMax<1.5&&abs(p3)<1&&isOK";
/// treeDump->Draw("meanTotP:ipad","meangTot>0&&isOK"+cutAll,"*")
///
/// FitFit(kTRUE)
/// StoreParam("local:///lustre/alice/akalweit/OCDBforMC") // specify corresponding location before !!
/// ~~~

#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TCut.h"
#include "TMatrixD.h"
#include "TVectorD.h"
//
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBRunRange.h"
#include "AliCDBStorage.h"
#include "AliTPCClusterParam.h"
//
#include "AliTPCcalibPID.h"
#include "AliTPCcalibDB.h"
#include "TStatToolkit.h"

AliTPCClusterParam * paramCl=0;
AliTPCcalibPID * pid =0;
TObjArray fitArr;
TTree * treeDump =0;
TTree * treeQtot;
TTree * treeQmax;
TTree * treeRatioQmax;
TTree * treeRatioQtot;
Int_t kmicolors[10]={1,2,3,4,6,7,8,9,10,11};
Int_t kmimarkers[10]={21,22,23,24,25,26,27,28,29,30};

//                     0      1   2    3        4   5    6    7
//                    dE/dx,  z, phi, theta,    p,  bg, ncls  type
//Int_t binsQA[7]    = {150, 10,  10,    10,   50, 50,  8};

void Init(char* name="calibPID06"){
  ///

  TFile fcalib("CalibObjectsTrain2.root");
  //TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib"); // old interface
  pid = ( AliTPCcalibPID *) fcalib.Get(name);
  TString axisName[9];
  axisName[0]  ="dE/dx"; axisName[1]  ="z (cm)";
  axisName[2]  ="sin(#phi)"; axisName[3]  ="tan(#theta)";
  axisName[4]  ="p (GeV)"; axisName[5]  ="#beta#gamma";
  axisName[6]  ="N_{cl}";

  pid->GetHistQtot()->SetTitle("Q_{tot};(z,sin(#phi),tan(#theta),p,betaGamma,ncls); TPC signal Q_{tot} (a.u.)");
  pid->GetHistQmax()->SetTitle("Q_{max};(z,sin(#phi),tan(#theta),p,betaGamma,ncls); TPC signal Q_{max} (a.u.)");

  for (Int_t ivar2=0;ivar2<7;ivar2++){      
    pid->GetHistQmax()->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
    pid->GetHistQtot()->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
    pid->GetHistRatioMaxTot()->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
    pid->GetHistRatioQmax()->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
    pid->GetHistRatioQtot()->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
    //
    pid->GetHistQmax()->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
    pid->GetHistQtot()->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());

    pid->GetHistRatioMaxTot()->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
    pid->GetHistRatioQmax()->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
    pid->GetHistRatioQtot()->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
  }
}


void StoreParam(char* localStorage = "local:///lustre/alice/akalweit/OCDB"){
  Int_t runNumber = 0;
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("AliTPCClusterParam");
  metaData->SetResponsible("Alexander Kalweit");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-24-00"); 
  metaData->SetComment("October runs calibration");
  AliCDBId id1("TPC/Calib/ClusterParam", runNumber, AliCDBRunRange::Infinity());
  AliCDBStorage *gStorage = AliCDBManager::Instance()->GetStorage(localStorage);
  gStorage->Put(paramCl, id1, metaData);

}




void SetRange(Int_t index, Float_t min, Float_t max){  
  pid->GetHistQmax()->GetAxis(index)->SetRangeUser(min,max);
  pid->GetHistQtot()->GetAxis(index)->SetRangeUser(min,max);

  pid->GetHistRatioMaxTot()->GetAxis(index)->SetRangeUser(min,max);
  pid->GetHistRatioQtot()->GetAxis(index)->SetRangeUser(min,max);
  pid->GetHistRatioQmax()->GetAxis(index)->SetRangeUser(min,max);
  pid->GetHistRatioTruncQtot()->GetAxis(index)->SetRangeUser(min,max);
  pid->GetHistRatioTruncQmax()->GetAxis(index)->SetRangeUser(min,max);
}


void SetType(Int_t type){
  pid->GetHistQmax()->GetAxis(7)->SetRange(type,type);
  pid->GetHistQtot()->GetAxis(7)->SetRange(type,type);
  //
  //
  pid->GetHistRatioMaxTot()->GetAxis(7)->SetRange(type,type);
  pid->GetHistRatioQtot()->GetAxis(7)->SetRange(type,type);
  pid->GetHistRatioQmax()->GetAxis(7)->SetRange(type,type);
  pid->GetHistRatioTruncQtot()->GetAxis(7)->SetRange(type,type);
  pid->GetHistRatioTruncQmax()->GetAxis(7)->SetRange(type,type);

}




void ReadTrees(){
  TFile f0("dumpQtot.root");
  treeQtot = (TTree*)f0.Get("Dump");
  TFile f1("dumpQmax.root");
  treeQmax = (TTree*)f1.Get("Dump");
  TFile f2("dumpRatioQtot.root");
  treeRatioQtot = (TTree*)f2.Get("Dump");
  TFile f3("dumpRatioQmax.root");
  treeRatioQmax = (TTree*)f3.Get("Dump");
}



void Fit(Bool_t updateParam=kFALSE){
  /// align pads

  TStatToolkit toolkit;
  Double_t chi2;
  TVectorD paramTot[5], paramMax[5];
  TString *strQT[4], *strQM[4];
  TVectorD paramTotRatio[5], paramMaxRatio[5];
  TString *strQTRatio[5], *strQMRatio[5];

  TMatrixD covar;
  Int_t npoints;
  treeQmax->SetAlias("cdr","(1-1/(1+dr))");
  treeQmax->SetAlias("cty","(1-1/sqrt(1+ty^2))");
  treeQmax->SetAlias("ctz","(1-1/sqrt(1+p3^2))");
  treeQtot->SetAlias("cdr","(1-1/(1+dr))");
  treeQtot->SetAlias("cty","(1-1/sqrt(1+ty^2))");
  treeQtot->SetAlias("ctz","(1-1/sqrt(1+p3^2))");
  //
  treeRatioQmax->SetAlias("cdr","(1-1/(1+dr))");
  treeRatioQmax->SetAlias("cty","(1-1/sqrt(1+ty^2))");
  treeRatioQmax->SetAlias("ctz","(1-1/sqrt(1+p3^2))");
  treeRatioQtot->SetAlias("cdr","(1-1/(1+dr))");
  treeRatioQtot->SetAlias("cty","(1-1/sqrt(1+ty^2))");
  treeRatioQtot->SetAlias("ctz","(1-1/sqrt(1+p3^2))");


  TString strFit="";
  //
  strFit+="dr++";
  strFit+="abs(p3)++";
  //
  strFit+="(1-(sy*sz)/0.5)++";
  strFit+="(1-(sy/0.5))++";
  strFit+="(1-sz)++";

  //
  TCut cutAll="(p>1&&dr>0.1&&val<3)"; 
  TCut cutPads[4];
  TCut cutPadsW[4];
  for (Int_t ipad=0;ipad<4;ipad++){
    printf("fitting pad\t%d\n",ipad);
    cutPads[ipad]=Form("abs(type-%f)<0.2&&(p>1&&dr>0.1&&val<3)",Double_t(ipad+1.5));
    cutPadsW[ipad]=Form("(abs(type-%f)<0.2&&(p>1&&dr>0.1&&val<3))*bincont",Double_t(ipad+1.5));
    //
    strQT[ipad] = toolkit.FitPlane(treeQtot,"val:1/sqrt(bincont)",strFit.Data(),cutAll+cutPads[ipad],chi2,npoints,paramTot[ipad],covar);
    strQM[ipad] = toolkit.FitPlane(treeQmax,"val:1/sqrt(bincont)",strFit.Data(),cutAll+cutPads[ipad],chi2,npoints,paramMax[ipad],covar);

    strQTRatio[ipad] = toolkit.FitPlane(treeRatioQtot,"val:1/sqrt(bincont)",strFit.Data(),cutAll+cutPads[ipad],chi2,npoints,paramTotRatio[ipad],covar);
    strQMRatio[ipad] = toolkit.FitPlane(treeRatioQmax,"val:1/sqrt(bincont)",strFit.Data(),cutAll+cutPads[ipad],chi2,npoints,paramMaxRatio[ipad],covar);

    treeRatioQmax->SetAlias(Form("fitQM%d",ipad),strQM[ipad]->Data()); 
    treeRatioQmax->SetAlias(Form("fitQT%d",ipad),strQT[ipad]->Data()); 
    treeRatioQtot->SetAlias(Form("fitQM%d",ipad),strQM[ipad]->Data()); 
    treeRatioQtot->SetAlias(Form("fitQT%d",ipad),strQT[ipad]->Data()); 
    treeQmax->SetAlias(Form("fitQM%d",ipad),strQM[ipad]->Data()); 
    treeQmax->SetAlias(Form("fitQT%d",ipad),strQT[ipad]->Data()); 
    treeQtot->SetAlias(Form("fitQM%d",ipad),strQM[ipad]->Data()); 
    treeQtot->SetAlias(Form("fitQT%d",ipad),strQT[ipad]->Data()); 
    
    treeQmax->SetAlias(Form("fitQMR%d",ipad),strQMRatio[ipad]->Data()); 
    treeQmax->SetAlias(Form("fitQTR%d",ipad),strQTRatio[ipad]->Data()); 
    treeQtot->SetAlias(Form("fitQMR%d",ipad),strQMRatio[ipad]->Data()); 
    treeQtot->SetAlias(Form("fitQTR%d",ipad),strQTRatio[ipad]->Data()); 
    treeRatioQmax->SetAlias(Form("fitQMR%d",ipad),strQMRatio[ipad]->Data()); 
    treeRatioQmax->SetAlias(Form("fitQTR%d",ipad),strQTRatio[ipad]->Data()); 
    treeRatioQtot->SetAlias(Form("fitQMR%d",ipad),strQMRatio[ipad]->Data()); 
    treeRatioQtot->SetAlias(Form("fitQTR%d",ipad),strQTRatio[ipad]->Data()); 
  }
  //
  //
  //
  for (Int_t ipad=0;ipad<4;ipad++){
    for (Int_t icorr=1;icorr<4;icorr++) {
      paramMax[ipad][icorr]/=paramMax[ipad][0];
      paramTot[ipad][icorr]/=paramTot[ipad][0];
      paramMaxRatio[ipad][icorr]/=paramMaxRatio[ipad][0];
      paramTotRatio[ipad][icorr]/=paramTotRatio[ipad][0];
    }
  }

  for (Int_t ipad=0;ipad<3;ipad++){
    for (Int_t icorr=1;icorr<4;icorr++) {
      paramMax[ipad][icorr]+=(paramMax[3][icorr]+paramMax[ipad][icorr])*0.5;
      paramTot[ipad][icorr]+=(paramTot[3][icorr]+paramTot[ipad][icorr])*0.5;
      if (updateParam){
	(*paramCl->fQNormCorr)(ipad+6,icorr) = paramTot[ipad][icorr+1];
	(*paramCl->fQNormCorr)(ipad+9,icorr) = paramMax[ipad][icorr+1];
      }

    }
  }


//  for (Int_t ipad=0;ipad<3;ipad++){
//    TVectorD *vecMax =  (TVectorD*)(paramCl->fQNormGauss->At(3*1+ipad));
//    TVectorD *vecTot =  (TVectorD*)(paramCl->fQNormGauss->At(3*0+ipad));
//    for (Int_t icorr=1;icorr<4;icorr++){
//      (*vecMax)[icorr]+=paramMax[ipad][icorr]/(*vecMax)[0];
//      (*vecTot)[icorr]+=paramTot[ipad][icorr]/(*vecTot)[0];
//    }
//  }


}

void fitQdep(){ 
  SetRange(6,100,160);
  SetRange(4,0.5,2);
  SetRange(0,0.5,1.5);
  TF1 f1("f1","([0]+[1]/x^2)",0.5,2);
  SetType(2);
  (pid->GetHistRatioQtot()->Projection(0,4))->ProfileX()->Fit(&f1);
  SetType(3);
  (pid->GetHistRatioQtot()->Projection(0,4))->ProfileX()->Fit(&f1);
  SetType(4);
  (pid->GetHistRatioQtot()->Projection(0,4))->ProfileX()->Fit(&f1);
  //
  SetType(2);
  (pid->GetHistRatioQmax()->Projection(0,4))->ProfileX()->Fit(&f1);
  SetType(3);
  (pid->GetHistRatioQmax()->Projection(0,4))->ProfileX()->Fit(&f1);
  SetType(4);
  (pid->GetHistRatioQmax()->Projection(0,4))->ProfileX()->Fit(&f1);

  
}

void LookupHisto(Int_t minTracks=200, Float_t minp=20, Float_t maxp=10000){
  TF1 f1("myg","gaus",0,10);
  Int_t dim[4]={0,1,2,3};
  pid->GetHistQtot()->GetAxis(6)->SetRangeUser(100,160);
  pid->GetHistQtot()->GetAxis(0)->SetRangeUser(0.1,6.); // important adaption to be done here ...
  pid->GetHistQmax()->GetAxis(6)->SetRangeUser(100,160);
  pid->GetHistQmax()->GetAxis(0)->SetRangeUser(0.1,6);  // important adaption to be done here ...

  pid->GetHistQmax()->GetAxis(4)->SetRangeUser(minp,maxp);
  pid->GetHistQtot()->GetAxis(4)->SetRangeUser(minp,maxp);
  
  
  THnSparse *hisTot[4];
  TH3F      *hisTotMean[4];
  Float_t    meanTotPad[4];
  Float_t    rmsTotPad[4];
  //
  THnSparse *hisMax[4];
  TH3F      *hisMaxMean[4];
  Float_t    meanMaxPad[4];
  Float_t    rmsMaxPad[4];

  TObjArray *array = new TObjArray(8);
    
  for (Int_t ipad=0;ipad<4;ipad++){
    pid->GetHistQtot()->GetAxis(7)->SetRange(2+ipad,2+ipad);
    pid->GetHistQmax()->GetAxis(7)->SetRange(2+ipad,2+ipad);
    hisTot[ipad] = pid->GetHistQtot()->Projection(4,dim);
    hisMax[ipad] = pid->GetHistQmax()->Projection(4,dim);
    Float_t drmin = hisTot[ipad]->GetAxis(1)->GetXmin();
    Float_t drmax = hisTot[ipad]->GetAxis(1)->GetXmax();
    Float_t p2min = hisTot[ipad]->GetAxis(2)->GetXmin();
    Float_t p2max = hisTot[ipad]->GetAxis(2)->GetXmax();
    Float_t p3min = hisTot[ipad]->GetAxis(3)->GetXmin();
    Float_t p3max = hisTot[ipad]->GetAxis(3)->GetXmax();
    meanTotPad[ipad] =hisTot[ipad]->Projection(0)->GetMean(); 
    meanMaxPad[ipad] =hisMax[ipad]->Projection(0)->GetMean(); 
    rmsTotPad[ipad] =hisTot[ipad]->Projection(0)->GetRMS(); 
    rmsMaxPad[ipad] =hisMax[ipad]->Projection(0)->GetRMS(); 
    
    hisTotMean[ipad]=new TH3F(Form("htot%d",ipad),Form("htot%d",ipad),
			      hisTot[ipad]->GetAxis(1)->GetNbins(), drmin,drmax,
			      hisTot[ipad]->GetAxis(2)->GetNbins(), p2min,p2max,
			      hisTot[ipad]->GetAxis(3)->GetNbins(), p3min,p3max);
    hisMaxMean[ipad]=new TH3F(Form("hmax%d",ipad),Form("hmax%d",ipad),
			      hisMax[ipad]->GetAxis(1)->GetNbins(), drmin,drmax,
			      hisMax[ipad]->GetAxis(2)->GetNbins(), p2min,p2max,
			      hisMax[ipad]->GetAxis(3)->GetNbins(), p3min,p3max);
    array->AddAt(hisTotMean[ipad],ipad);
    array->AddAt(hisMaxMean[ipad],ipad+4);
  }
  
  for (Int_t ibindr=1;ibindr<hisTot[0]->GetAxis(1)->GetNbins()+1; ibindr++)
    for (Int_t ibinty=1;ibinty<hisTot[0]->GetAxis(2)->GetNbins()+1; ibinty++)
      for (Int_t ibintz=1;ibintz<hisTot[0]->GetAxis(3)->GetNbins()+1; ibintz++)
	for (Int_t ipad=0;ipad<4; ipad++){
	  hisTotMean[ipad]->SetBinContent(ibindr,ibinty,ibintz,meanTotPad[ipad]);
	  hisMaxMean[ipad]->SetBinContent(ibindr,ibinty,ibintz,meanMaxPad[ipad]);
	}
  
  
  TTreeSRedirector * cstream = new TTreeSRedirector("lookupdEdx.root");
  
  for (Int_t ibindr=1;ibindr<hisTot[0]->GetAxis(1)->GetNbins()+1; ibindr+=1)
    for (Int_t ibinty=1;ibinty<hisTot[0]->GetAxis(2)->GetNbins()+1; ibinty+=1)
      for (Int_t ibintz=1;ibintz<hisTot[0]->GetAxis(3)->GetNbins()+1; ibintz+=1)
	for (Int_t ipad=0;ipad<4; ipad++){
	  //
	  Float_t dr =  hisTot[ipad]->GetAxis(1)->GetBinCenter(ibindr);
	  Float_t p2 =  hisTot[ipad]->GetAxis(2)->GetBinCenter(ibinty);
	  Float_t p3 =  hisTot[ipad]->GetAxis(3)->GetBinCenter(ibintz);
	  Int_t sumTot=0, sumMax=0;
	  TH1D * hisproTot =0; 
	  TH1D * hisproMax =0;
	  Float_t delta=0.1;
	  while(1){
	    hisTot[ipad]->GetAxis(1)->SetRangeUser(dr-delta,dr+delta);
	    hisTot[ipad]->GetAxis(2)->SetRangeUser(p2-delta,p2+delta);
	    hisTot[ipad]->GetAxis(3)->SetRangeUser(p3-delta,p3+delta);
	    //
	    hisMax[ipad]->GetAxis(1)->SetRangeUser(dr-delta,dr+delta);
	    hisMax[ipad]->GetAxis(2)->SetRangeUser(p2-delta,p2+delta);
	    hisMax[ipad]->GetAxis(3)->SetRangeUser(p3-delta,p3+delta);
	    // 4 sigma cut
	    hisTot[ipad]->GetAxis(0)->SetRangeUser(meanTotPad[ipad]-4.*rmsTotPad[ipad],meanTotPad[ipad]+4.*rmsTotPad[ipad]);
	    hisMax[ipad]->GetAxis(0)->SetRangeUser(meanMaxPad[ipad]-4.*rmsMaxPad[ipad],meanMaxPad[ipad]+4.*rmsMaxPad[ipad]);
	    hisproTot = hisTot[ipad]->Projection(0);
	    hisproMax = hisMax[ipad]->Projection(0);
	    sumTot  = hisproTot->GetSum();
	    sumMax  = hisproMax->GetSum();
	    if (sumTot>minTracks) break;
	    if (delta>0.35) break;
	    delete hisproTot;
	    delete hisproMax;
	    delta+=0.1;
	  }
	  //
	  Float_t meanTot = hisproTot->GetMean();
	  Float_t rmsTot  = hisproTot->GetRMS();
	  Float_t meanMax = hisproMax->GetMean();
	  Float_t rmsMax  = hisproMax->GetRMS();
	  //

	  Float_t meangTot = -1;
	  Float_t rmsgTot  = -1;
	  Float_t meangMax = -1;
	  Float_t rmsgMax  = -1;
	  if(sumTot>minTracks/2 && rmsTot>0.02&&rmsMax>0.02){
	    f1.SetParameters(hisproTot->GetMaximum(),meanTot,rmsTot);
	    hisproTot->Fit(&f1,"QNR","",0.1,10);
	    meangTot = f1.GetParameter(1);
	    rmsgTot  = f1.GetParameter(2);
	    f1.SetParameters(hisproMax->GetMaximum(),meanMax,rmsMax);
	    hisproMax->Fit(&f1,"QNR","",0.1,10);
	    meangMax = f1.GetParameter(1);
	    rmsgMax  = f1.GetParameter(2);
	  }

	  TH1 *htemp = 0;
	  htemp = hisTot[ipad]->Projection(1);
	  Float_t drc =  htemp->GetMean();
	  delete htemp;
	  htemp = hisTot[ipad]->Projection(2);
	  Float_t p2c =  htemp->GetMean();
	  delete htemp;
	  htemp = hisTot[ipad]->Projection(3);
	  Float_t p3c =  htemp->GetMean();
	  delete htemp;
	  //
	  Bool_t isOK=kTRUE;
	  if (meanTot==0)   {meanTot=meanTotPad[ipad]; isOK=kFALSE;}
	  if (meanMax==0)   {meanMax=meanMaxPad[ipad]; isOK=kFALSE;}
	  //
	  if (TMath::Abs(rmsTot/meanTot-0.12)>0.04) {meanTot=meanTotPad[ipad]; isOK=kFALSE;}
	  if (TMath::Abs(rmsMax/meanMax-0.12)>0.04) {meanMax=meanMaxPad[ipad]; isOK=kFALSE;}
	  //
	  hisTotMean[ipad]->SetBinContent(ibindr,ibinty,ibintz,meanTot);
	  hisMaxMean[ipad]->SetBinContent(ibindr,ibinty,ibintz,meanMax);
	  //
	  printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\n", ipad, dr,p2,p3, meanTot, rmsTot, meangTot);
	  printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\n", ipad, dr,p2,p3, meanMax, rmsMax, meangMax);
	  delete hisproTot;
	  delete hisproMax;
	  (*cstream)<<"dumpdEdx"<<
	    "ipad="<<ipad<<
	    "isOK="<<isOK<<
	    "sumTot="<<sumTot<<
	    "sumMax="<<sumMax<<
	    "meanTotP="<<meanTotPad[ipad]<<
	    "meanMaxP="<<meanMaxPad[ipad]<<
	    "meanTot="<<meanTot<<
	    "rmsTot="<<rmsTot<<
	    "meanMax="<<meanMax<<
	    "rmsMax="<<rmsMax<<
	    //
	    "meangTot="<<meangTot<<
	    "rmsgTot="<<rmsgTot<<
	    "meangMax="<<meangMax<<
	    "rmsgMax="<<rmsgMax<<
	    //
	    "dr="<<dr<<
	    "p2="<<p2<<
	    "p3="<<p3<<
	    "drc="<<drc<<
	    "p2c="<<p2c<<
	    "p3c="<<p3c<<
	    "\n";
	}
  TFile *fstream = cstream->GetFile();
  fstream->cd();
  array->Write("histos",TObject::kSingleKey);
  delete cstream;  
}

void FitFit(Bool_t updateParam=kFALSE){
  
  TStatToolkit toolkit;
  Double_t chi2Tot[5],chi2Max[5];
  TVectorD paramTot[5], paramMax[5];
  TString *strQT[4], *strQM[4];
  TVectorD paramTotRatio[5], paramMaxRatio[5];
  TString *strQTRatio[5], *strQMRatio[5];
  TMatrixD covar;
  Int_t npoints;
  TCut cutPads[5];
  treeDump->SetAlias("drm","(0.5-abs(drc))");
  treeDump->SetAlias("dri","(1-abs(drc))");
  treeDump->SetAlias("ty","tan(asin(p2c))**2");
  treeDump->SetAlias("tz","(abs(p3c)**2*sqrt(1+ty))");
  //
  TString strFit="";
  strFit+="drm++";
  strFit+="ty++";
  strFit+="tz++";
  //
  //strFit+="sqrt(dri*ty)++";
  //strFit+="sqrt(dri*tz)++";
  //strFit+="sqrt(ty*tz)++";
  
  TCut cutAll = "meangTot>0.0&&sumMax>150&&sumTot>150&&rmsgMax/rmsMax<1.5&&abs(p3)<1&&isOK"; // MY MODIFICATION ! isOK added
  for (Int_t ipad=0;ipad<3;ipad++){ // MY MODIFICATION ! <3 instead of 4
    cutPads[ipad]=Form("ipad==%d",ipad);
    //
    strQT[ipad] = toolkit.FitPlane(treeDump,"meangTot/meanTotP",strFit.Data(),cutAll+cutPads[ipad],chi2Tot[ipad],npoints,paramTot[ipad],covar);
    chi2Tot[ipad]=TMath::Sqrt(chi2Tot[ipad]/npoints);
    printf("Tot%d\t%f\t%s\t\n",ipad,chi2Tot[ipad],strQT[ipad]->Data());
    //
    strQM[ipad] = toolkit.FitPlane(treeDump,"meangMax/meanMaxP",strFit.Data(),cutAll+cutPads[ipad],chi2Max[ipad],npoints,paramMax[ipad],covar);
    chi2Max[ipad]=TMath::Sqrt(chi2Max[ipad]/npoints);
    printf("Max%d\t%f\t%s\t\n",ipad,chi2Max[ipad],strQM[ipad]->Data()); 
    //
    strQTRatio[ipad] = toolkit.FitPlane(treeDump,"meangTot/meangMax",strFit.Data(),cutAll+cutPads[ipad],chi2Max[ipad],npoints,paramTotRatio[ipad],covar);
    chi2Max[ipad]=TMath::Sqrt(chi2Max[ipad]/npoints);
    printf("Ratio%d\t%f\t%s\t\n\n",ipad,chi2Max[ipad],strQTRatio[ipad]->Data()); 
    //
    treeDump->SetAlias(Form("fitQT%d",ipad),strQT[ipad]->Data());
    treeDump->SetAlias(Form("fitQM%d",ipad),strQM[ipad]->Data());
    treeDump->SetAlias(Form("fitQTM%d",ipad),strQM[ipad]->Data());
  }


  TMatrixD mat(6,6);
  for (Int_t ipad=0; ipad<3; ipad++){
    for (Int_t icorr=0; icorr<3; icorr++){
      if (updateParam){
	(*paramCl->fQNormCorr)(ipad+6,icorr) += paramTot[ipad][icorr+1];
	(*paramCl->fQNormCorr)(ipad+9,icorr) += paramMax[ipad][icorr+1];
      }
      mat(ipad+0,icorr) = paramTot[ipad][icorr+1];
      mat(ipad+3,icorr) = paramMax[ipad][icorr+1];
    }
  }

  Float_t normShortTot;
  Float_t normMedTot;
  Float_t normLongTot;
  Float_t normTotMean;
  //
  Float_t normShortMax;
  Float_t normMedMax;
  Float_t normLongMax;
  Float_t normMaxMean;
  TH1D * dummyHist = new TH1D("dummyHist","absolute gain alignment of pads",100,0,10);
  //
  treeDump->Draw("meanTotP>>dummyHist","meangTot>0&&isOK&&ipad==0"+cutAll,"");
  normShortTot = dummyHist->GetMean();
  treeDump->Draw("meanTotP>>dummyHist","meangTot>0&&isOK&&ipad==1"+cutAll,"");
  normMedTot = dummyHist->GetMean();
  treeDump->Draw("meanTotP>>dummyHist","meangTot>0&&isOK&&ipad==2"+cutAll,"");
  normLongTot = dummyHist->GetMean();
  normTotMean = (normShortTot+normMedTot+normLongTot)/3.;
  //
  treeDump->Draw("meanMaxP>>dummyHist","meangTot>0&&isOK&&ipad==0"+cutAll,"");
  normShortMax = dummyHist->GetMean();
  treeDump->Draw("meanMaxP>>dummyHist","meangTot>0&&isOK&&ipad==1"+cutAll,"");
  normMedMax = dummyHist->GetMean();
  treeDump->Draw("meanMaxP>>dummyHist","meangTot>0&&isOK&&ipad==2"+cutAll,"");
  normLongMax = dummyHist->GetMean();
  normMaxMean = (normShortMax+normMedMax+normLongMax)/3.;
  
  cout << "Tot: " << normShortTot << " " << normMedTot << " " << normLongTot << endl;
  cout << "Max: " << normShortMax << " " << normMedMax << " " << normLongMax << endl;

  // hand alignment of pads to be improved:
  (*paramCl->fQNormCorr)(0,5) *= (normShortTot/normTotMean);
  (*paramCl->fQNormCorr)(1,5) *= (normMedTot/normTotMean);
  (*paramCl->fQNormCorr)(2,5) *= (normLongTot/normTotMean);
  //
  (*paramCl->fQNormCorr)(3,5) *= (normShortMax/normMaxMean);
  (*paramCl->fQNormCorr)(4,5) *= (normMedMax/normMaxMean);
  (*paramCl->fQNormCorr)(5,5) *= (normLongMax/normMaxMean);


}



