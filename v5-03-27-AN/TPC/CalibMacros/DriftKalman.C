
// See http://en.wikipedia.org/wiki/Kalman_filter
//
/*
.x ~/UliStyle.C
gSystem->AddIncludePath("-I$ALICE_ROOT/STAT");
gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
gSystem->AddIncludePath("-I$ALICE_ROOT/TPC");
gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
gSystem->Load("libSTAT.so");

.L $ALICE_ROOT/TPC/CalibMacros/DriftKalman.C++ 
  //
  //
  //
  TFile fgoofie("goofieTree.root");
  TTree *chainTime = (TTree*)fgoofie.Get("drift");
  chainTime->SetAlias("mcut","abs(areaF/areaN-1.117)<0.1");
  chainTime->SetAlias("atpcor","(1000*press/(3.329*(tempF)))-1");
  chainTime->SetAlias("dvdrift","(((vdrift/2.636)-1))");
  chainTime->SetAlias("trCorr","0");
 
  Int_t npoints= chainTime->Draw("time:(dvdrift)-trCorr:(atpcor)","mcut");
  Double_t *timeV = new Double_t[npoints];
  Double_t *vdriftV = new Double_t[npoints];
  Double_t *cPTV = new Double_t[npoints];
  memcpy(timeV,chainTime->GetV1(),npoints*sizeof(Double_t));
  memcpy(vdriftV,chainTime->GetV2(),npoints*sizeof(Double_t));
  memcpy(cPTV,chainTime->GetV3(),npoints*sizeof(Double_t));

  Double_t rms=0.006, corr;
  covXk(0,0)=0.01; covXk(1,1)=0.1; covXk(2,2)=0.1; covXk(1,0)=0; covXk(0,1)=0;
  vecXk(0,0)=0; vecXk(1,0)=1; vecXk(2,0)=0;          //initial stat

  recVDriftSort(npoints,timeV,vdriftV,cPTV,rms,10/(250*24*60*60.),kFALSE);
  recVDriftSort(npoints,timeV,vdriftV,cPTV,rms,10/(250*24*60*60.),kTRUE); 
  TFile f("realvdrift.root");
  drift->Draw("vecXk.fElements[1]:time","isOK","")
 

  //
  // Real data test usage
  //

  AliXRDPROOFtoolkit tool;
  TChain * chainTime = tool.MakeChain("time.txt","timeInfo",0,10200);
  chainTime->Lookup();
  //
  
  Double_t rms, corr;
  chainTime->SetAlias("atpcorP","press0/(3.329*(temp0+273.15))-1");
  chainTime->SetAlias("mcutP","press0>0");
  chainTime->SetAlias("ncutP","press0<1");
  chainTime->SetAlias("atpcorG","temp1.fElements[14]/(0.003378*(temp0+273.15))-1");
  chainTime->SetAlias("mcutG","temp1.fElements[14]>0");

  chainTime->SetAlias("atpcor","(mcutP*atpcorP+(ncutP)*atpcorG)");
  chainTime->SetAlias("mcut","temp1.fElements[14]>0");
  chainTime->SetAlias("dvdrift","(fDz/500)");
  
  fitLinear(chainTime, rms,corr);
 

  Int_t npoints= chainTime->Draw("time:(dvdrift)-trCorr:(atpcor)","mcut&&abs(fDz-500*driftG)<1");
  Double_t *timeV = new Double_t[npoints];
  Double_t *vdriftV = new Double_t[npoints];
  Double_t *cPTV = new Double_t[npoints];
  memcpy(timeV,chainTime->GetV1(),npoints*sizeof(Double_t));
  memcpy(vdriftV,chainTime->GetV2(),npoints*sizeof(Double_t));
  memcpy(cPTV,chainTime->GetV3(),npoints*sizeof(Double_t));
  

  kalmanTime.Init(0,0,-1,0.01,0.01);

  recVDriftSort(npoints,timeV,vdriftV,cPTV,rms,1/(250*24*60*60.),kFALSE);
  recVDriftSort(npoints,timeV,vdriftV,cPTV,rms,1/(250*24*60*60.),kTRUE); 

  TFile f("realvdrift.root");
  drift->Draw("k.fState.fElements[0]:time","isOK","");
*/

#include "TMatrixD.h"
#include "TRandom.h"
#include "TTreeStream.h"
#include "TMath.h"
#include "TString.h"
#include "TChain.h"
#include "TStatToolkit.h"
#include "AliTPCkalmanTime.h"


Double_t cshift  = 0;
Double_t cPT     = 0;
Double_t floatS  = 0.00005;   // 
Double_t floatTP = 0.0005;    //
Double_t noiseV  = 0.3/250.;  // cosmic drift precission
Double_t vdrift    = 0;

AliTPCkalmanTime kalmanTime;


void recVDrift(Int_t npoints, Double_t *timeV, Double_t *vdriftV, Double_t * cPTV, Double_t noise, Double_t sshift){
  //
  // chainTime->Draw("time:(dvdrift):(atpcor)","mcut&&abs(((dvdrift))-(atpcor))<0.004");
  //
  TTreeSRedirector *pcstream = new TTreeSRedirector("realvdrift.root");
  //
  Int_t nused=0;
  for (Int_t ipoint=0; ipoint<npoints;ipoint++){    
    //
    if (ipoint==0) kalmanTime.fTime=timeV[ipoint];
    TMatrixD &vecXk = *(kalmanTime.fState);
    TMatrixD &covXk = *(kalmanTime.fCovariance);
    Double_t fdrift = vecXk(0,0)+ cPTV[ipoint]*vecXk(1,0);
    Double_t ddrift = fdrift-vdriftV[ipoint];
    Double_t sdrift = TMath::Sqrt(covXk(0,0)+noise*noise);
    Bool_t isOK = nused<10||TMath::Abs(ddrift)<5*sdrift;
    kalmanTime.Propagate(timeV[ipoint],sshift);
    kalmanTime.Update(vdriftV[ipoint],noise,cPTV[ipoint],pcstream);
    nused++;
    (*pcstream)<<"drift"<<
      "ipoint="<<ipoint<<      
      "time="<<timeV[ipoint]<<
      "isOK="<<isOK<<
      "noise="<<noise<<
      "cPT="<<cPTV[ipoint]<<
      "vdrift="<<vdriftV[ipoint]<<
      "fdrift="<<fdrift<<
      "k.="<<&kalmanTime<<
      "\n";
  }
  delete pcstream;
}


void recVDriftSort(Int_t npoints, Double_t *timeV, Double_t *vdriftV, Double_t * cPTV, Double_t noise, Double_t sshift, Bool_t sort){
  //
  // sort entries before fit
  //
  Int_t *index = new Int_t[npoints];
  Double_t *tmpsort = new Double_t[npoints];
  TMath::Sort(npoints, timeV, index,sort);
  //
  for (Int_t i=0; i<npoints;i++){
    tmpsort[i]= timeV[index[i]];
  }
  for (Int_t i=0; i<npoints;i++){
    timeV[i] =tmpsort[i];
  }
  //
  for (Int_t i=0; i<npoints;i++){
    tmpsort[i]= vdriftV[index[i]];
  }
  for (Int_t i=0; i<npoints;i++){
    vdriftV[i] =tmpsort[i];
  }
  //
  for (Int_t i=0; i<npoints;i++){
    tmpsort[i]= cPTV[index[i]];
  }
  for (Int_t i=0; i<npoints;i++){
    cPTV[i] =tmpsort[i];
  }
  delete [] index;
  delete [] tmpsort;
  recVDrift(npoints,timeV,vdriftV,cPTV,noise,sshift);
}


void  fitLinear(TChain * chainTime, Double_t& rms, Double_t & corr){
  //
  // Fit globaly - make aliases;
  //
  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD fitParam;
  TMatrixD covMatrix;
  TString  fitStr="";  
  fitStr+="(atpcor)++";
  fitStr+="(trigger==1)++";
  fitStr+="(trigger==2)++";
  fitStr+="(trigger==4)++";
  fitStr+="(trigger==8)++";

  TString *strvdpt1 = 0;
  strvdpt1 = toolkit.FitPlane(chainTime,"(dvdrift)",fitStr.Data(),"mcut" , chi2,npoints,fitParam,covMatrix); 
  chainTime->SetAlias("driftG",strvdpt1->Data());

  strvdpt1 =toolkit.FitPlane(chainTime,"(dvdrift)",fitStr.Data(),"mcut&&abs(fDz-500*driftG)<2" , chi2,npoints,fitParam,covMatrix);    
  chainTime->SetAlias("driftG",strvdpt1->Data());
  
 
  rms = TMath::Sqrt(chi2/npoints);
  corr = fitParam[1];  // cPT parameter
  //
  //
  fitStr="";
  fitStr+="(trigger==1)++";
  fitStr+="(trigger==2)++";
  fitStr+="(trigger==4)++";
  fitStr+="(trigger==8)++";
  TString fitExpr = "(dvdrift)-";
  fitExpr+=fitParam[1];
  fitExpr+="* (atpcor)";
  TString *strCorr = toolkit.FitPlane(chainTime,fitExpr.Data(),fitStr.Data(),"abs(fDz-500*driftG)<1" , chi2,npoints,fitParam,covMatrix);  
  chainTime->SetAlias("trCorr",strCorr->Data());


  fitStr="";  
  fitStr+="(atpcor)++";
  fitStr+="(atpcor)^2++";
  fitStr+="(trigger==1)++";
  fitStr+="(trigger==2)++";
  fitStr+="(trigger==4)++";
  fitStr+="(trigger==8)++";
  TString *strvdpt2 = 0;
  strvdpt2 =toolkit.FitPlane(chainTime,"(dvdrift)",fitStr.Data(),"mcut&&abs(fDz-500*driftG)<2" , chi2,npoints,fitParam,covMatrix);    
  chainTime->SetAlias("driftG2",strvdpt2->Data());


}
