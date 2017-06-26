/*
Class for fitting of distortion maps using phsyical models

*/

/*
  .x $NOTES/aux/NimStyle.C(2)
  gSystem->AddIncludePath("-I$AliRoot_SRC/TPC/TPCcalib/");
  .L $AliRoot_SRC/TPC/TPCcalib/AliTPCDistortionFit.cxx+
  AliTPCDistortionFit::RegisterFitters();

  AliTPCDistortionFit::LoadDistortionMaps(245683, "local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/");
  AliTPCDistortionFit::LoadDistortionMaps(246855, "local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/");

  AliTPCDistortionFit::LoadDistortionTree("/data/alien/alice/data/2015/LHC15o/000245683/cpass0_pass1/ResidualMerge/TPCSPCalibration/1448920523_1448922567_000245683/voxelResTree.root");
  AliTPCDistortionFit::LoadDistortionTree("/data/alien/alice/data/2015/LHC15o/000245683/cpass0_pass1/ResidualMerge/TPCSPCalibration/1448918479_1448920523_000245683/voxelResTree.root");
  AliTPCDistortionFit::LoadDistortionTree("/data/alien/alice/data/2015/LHC15o/000245683/cpass0_pass1/ResidualMerge/TPCSPCalibration/1448922567_1448924611_000245683/voxelResTree.root");


*/
#include <fstream>
#include <iostream>
#include <map>
#include <sys/types.h>
#include <sys/stat.h>
#include <string>
#include "TMath.h"
#include "TTree.h"
#include "TKey.h"
#include "TStatToolkit.h"
#include "TPRegexp.h"
#include "AliTPCChebCorr.h"
#include "AliTPCChebCorr.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGRPObject.h"
#include "AliTMinuitToolkit.h"
#include "TStyle.h"
#include "AliTPCROC.h"
#include "AliLumiTools.h"
#include "AliTPCRecoParam.h"
#include "TLine.h"
#include "AliExternalInfo.h"
#include "AliLHCData.h"
#include "AliTPCDistortionFit.h"
using namespace std;
 
using std::cout;


ClassImp(AliTPCDistortionFit)

map<string,int>    AliTPCDistortionFit::fgkMapNameHash;              // map  name->index
map<int,string>  AliTPCDistortionFit::fgkMapHashName;                // map  index->name
map<int,const AliTPCChebCorr *> AliTPCDistortionFit::fgkMapCheb;     // map  index->Cheb param
TTree *  AliTPCDistortionFit::fgkDistortionTree;     


Double_t AliTPCDistortionFit::LineFieldLocal(const Double_t *x, const Double_t *param){
  // Simple E filed
  // x        - input vector
  //   x[0]   - return vector description
  //     0    - Er
  //     1    - Erphi
  //     2    - Dr
  //     3    - Drphi
  //     4    - Dr-bckg
  //     5    - Drphi-bckg
  //
  //   x[1]   - r
  //   x[2]   - rphi
  //   x[3]   - z
  // param    - description of charge
  //      [0] - line Q
  //      [1] - r0
  //      [2] - rphi0
  //      [3] - scale
  //      [4] - wt
  //      [5] - symmetry plane for mirror charge - not used
  //      [6] - E nominal
  //
  if (x[0]==4 || x[0]==5 ){
    Double_t drphi=TMath::Pi()*0.3*x[1]/9.;    
    Double_t values[4], distance2[4];
    Double_t xxx[4];
    for (Int_t i=0; i<4; i++){   // calculate si
      xxx[0]=x[0]-2;   xxx[1]=x[1]; xxx[2]=0; xxx[3]=x[3];
      xxx[1]=(i<2)?    x[1]-10:x[1]+10;
      xxx[2]=(i%2==0)? -drphi:drphi;
      values[i]=AliTPCDistortionFit::LineFieldLocal(xxx,param);
      distance2[i]=(x[1]-xxx[1])*(x[1]-xxx[1])+  (x[2]-xxx[2])*(x[2]-xxx[2]); 
    }
    Double_t sum=0, sumw=0;
    for (Int_t i=0; i<4; i++){
      sum+=values[i]/(1.+distance2[i]);
      sumw+=1/(1.+distance2[i]);
    }
    Double_t bckg=sum/sumw;
    xxx[0]=x[0]-2;  xxx[1]=x[1]; xxx[2]=x[2]; xxx[3]=x[3];
    Double_t val=AliTPCDistortionFit::LineFieldLocal(xxx,param);
    return val-bckg;
  }

  Double_t rdist2= (x[1]-param[1])*(x[1]-param[1])+(x[2]-param[2])*(x[2]-param[2]);
  Double_t eField=param[0]/TMath::Sqrt(rdist2+param[3]*param[3]);
  Double_t eR=eField*(x[1]-param[1])/TMath::Sqrt(rdist2);
  Double_t eRPhi=eField*(x[2]-param[2])/TMath::Sqrt(rdist2);
  //
  Double_t xmirror=param[5]-(param[1]-param[5]);
  Double_t rdist2M= (x[1]-xmirror)*(x[1]-xmirror)+(x[2]-param[2])*(x[2]-param[2]);
  Double_t eFieldM=param[0]/TMath::Sqrt(rdist2M+param[3]*param[3]);
  Double_t eRM=eFieldM*(x[1]-xmirror)/TMath::Sqrt(rdist2M);
  Double_t eRPhiM=eFieldM*(x[2]-param[2])/TMath::Sqrt(rdist2M);
  //
  eR+=eRM;
  eRPhi+=eRPhiM;
  //
  if (x[0]==0) return eR;
  if (x[0]==1) return eRPhi;
  Double_t drift=250-TMath::Abs(x[3]);
  Double_t dR=    drift*(eR-param[4]*eRPhi)/param[6];
  Double_t dRPhi= drift*(eR*param[4]+eRPhi)/param[6];
  if (x[0]==2) return dR;
  if (x[0]==3) return dRPhi;
}

Double_t AliTPCDistortionFit::NLinesFieldLocal(const Double_t *x, const Double_t *param){
  //  Mulitiline response
  //   x       - input vector
  //   x[0]   - return vector description
  //     0    - Er
  //     1    - Erphi
  //     2    - Dr
  //     3    - Drphi
  //     4    - Dr-bckg
  //     5    - Drphi-bckg
  //
  //   x[1]   - r
  //   x[2]   - rphi
  //   x[3]   - z
  // param    - description of charge
  //
  //      [0] - scale
  //      [1] - wt
  //      [2] - symmetry plane for mirror charge - not used
  //      [3] - E nominal
  //      [4] - nlines
  //      [5] - line Q
  //      [6] - r0
  //      [7] - rphi0 ...

  Double_t paramLocal[7];
  paramLocal[3]=param[0]; // scale 
  paramLocal[4]=param[1]; // wt
  paramLocal[5]=param[2]; // symetry axis
  paramLocal[6]=param[3]; // E nominal
  Int_t nLines=param[4];
  if (nLines<=0) return 0;
  Double_t sum=0;
  for (Int_t iLine=0; iLine<nLines; iLine++){
    paramLocal[0]=param[5+iLine*3];
    paramLocal[1]=param[6+iLine*3];
    paramLocal[2]=param[7+iLine*3];
    sum+=AliTPCDistortionFit::LineFieldLocal(x,paramLocal);
  }
  return sum;
}


Int_t AliTPCDistortionFit::LoadDistortionMaps(Int_t run, const char *storage){
  //
  // Loading and registering reference correction maps for given run
  //
  /*
    AliTPCDistortionFit::LoadDistortionMaps(245683, "local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/");
    AliTPCDistortionFit::LoadDistortionMaps(246855, "local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/");
    AliTPCDistortionFit::PrintMap();
  */
  AliCDBManager * man  = AliCDBManager::Instance();
  if (storage) man->SetDefaultStorage(storage);
  man->SetRun(run);
  AliCDBEntry* entry=0;
  entry = man->Get("GRP/GRP/Data");  
  AliGRPObject * grp = (AliGRPObject *) entry->GetObject();
  Double_t polarity = (grp->GetL3Polarity()<=0)? 1:-1;
  //
  TString mapRefName=TString::Format("R%d.refMap",run);    
  AliCDBEntry* entryC = man->Get("TPC/Calib/CorrectionMapsRef");  
  TObjArray *referenceCheb=(TObjArray *)entryC->GetObject();
  AliTPCChebCorr *cheb=(grp->GetL3Polarity()<=0) ?(AliTPCChebCorr *)referenceCheb->At(0) : (AliTPCChebCorr *)referenceCheb->At(1);
  if (cheb==NULL) cheb=(AliTPCChebCorr *)referenceCheb->At(0); // what is the logic- why sometimes 2 maps somtimes 1 map???? to aks Ruben
  RegisterMap(mapRefName,cheb);
  cheb->Init();
  //
  AliCDBEntry* entryM = man->Get("TPC/Calib/DistortionMaps");  
  TObjArray *distortionMaps=(TObjArray *)entryM->GetObject();
  Int_t entries= distortionMaps->GetEntries();
  for (Int_t i=0; i<entries; i++){
    AliTPCChebCorr *chebMap=(AliTPCChebCorr *)distortionMaps->At(i);
    RegisterMap(chebMap->GetName(), cheb);
    chebMap->Init();
  }
  //delete entry;    //  
  //  delete entryC;   //  memory leak here
}


TString AliTPCDistortionFit::LoadDistortionTree(const char *chinput){
  //
  // Load distortion tree and set metadata
  /*
    chinput="/data/alien/alice/data/2015/LHC15o/000245683/cpass0_pass1/ResidualMerge/TPCSPCalibration/1448920523_1448922567_000245683/voxelResTree.root";
    chiinput="/data/alien/alice/data/2015/LHC15o/000246855/cpass0_pass1/ResidualMerge/TPCSPCalibration/1449886884_1449888523_000246855/voxelResTree.root";
    AliTPCDistortionFit::LoadDistortionTree("/data/alien/alice/data/2015/LHC15o/000245683/cpass0_pass1/ResidualMerge/TPCSPCalibration/1448920523_1448922567_000245683/voxelResTree.root");
  */  
  TFile * finput=TFile::Open(chinput);
  TTree * treeLoad =(TTree*)finput->Get("voxRes");
  if (treeLoad==NULL){
    ::Error(" AliTPCDistortionFit::LoadDistortionTree","File or tree not available %s",chinput);
    return "";
  }
  AliTPCChebCorr * cheb=0;    
  for (Int_t i=0; i<finput->GetListOfKeys()->GetEntries(); i++){             // name of Cheb can be different - but should be only one
    TString clname=((TKey*)finput->GetListOfKeys()->At(i))->GetClassName();
    if (clname.Contains("AliTPCChebCor")) {
      cheb = (AliTPCChebCorr*)finput->Get(finput->GetListOfKeys()->At(i)->GetName());
      cheb->Init();
    }
  }
  if (cheb==NULL){
    ::Error(" AliTPCDistortionFit::LoadDistortionTree","Cheb map not avalible in file  %s",chinput);
    return "";
  }
  Int_t run=cheb->GetRun();
  RegisterMap(cheb->GetName(),cheb);
  treeLoad->SetMarkerStyle(21);
  treeLoad->SetMarkerSize(0.6);
  if (fgkDistortionTree==NULL) {
    TFile * finput2=TFile::Open(chinput);
    fgkDistortionTree=(TTree*)finput2->Get("voxRes");
    SetMetadata(fgkDistortionTree,"",run);
    fgkDistortionTree->SetMarkerStyle(25);
    fgkDistortionTree->SetMarkerSize(0.4);
  }  
  TString refMap=TString::Format("R%d.refMap",run);
  if (GetCheb(refMap.Data())==NULL) LoadDistortionMaps(run,NULL);
  fgkDistortionTree->AddFriend(treeLoad,cheb->GetName());
  SetMetadata(fgkDistortionTree,cheb->GetName(),run); 
  return cheb->GetName();
}


void  AliTPCDistortionFit::SetMetadata(TTree * tree, TString friendName, Int_t run){
  //
  //
  //
  TString  refMap=TString::Format("R%d.refMap",run);
  Int_t hashRef=refMap.Hash();
  Int_t    friendRef=friendName.Hash();

  TString prefix="";
  if (friendName.Length()>0) {
    prefix=friendName+".";
    TTree *ftree=tree->GetFriend(friendName);
    TList * alist= tree->GetListOfAliases();
    TPRegexp matchF("R.*.");
    for (Int_t i=0; i<alist->GetEntries(); i++){
      if (matchF.Match(alist->At(i)->GetName())) continue;
      ftree->SetAlias((alist->At(i)->GetName()),alist->At(i)->GetTitle());
    }
    ftree->SetAlias("ddX",TString::Format("%s.D[0]-AliTPCDistortionFit::Eval(%d,sector,x,y2xAV,z2xAV,0+0)",friendName.Data(),hashRef).Data());
    ftree->SetAlias("ddY",TString::Format("%s.D[1]-AliTPCDistortionFit::Eval(%d,sector,x,y2xAV,z2xAV,1+0)",friendName.Data(),hashRef).Data());
    ftree->SetAlias("ddZ",TString::Format("%s.D[2]-AliTPCDistortionFit::Eval(%d,sector,x,y2xAV,z2xAV,2+0)",friendName.Data(),hashRef).Data());
    ftree->SetAlias("ddXS",TString::Format("%s.DS[0]-AliTPCDistortionFit::Eval(%d,sector,x,y2xAV,z2xAV,0+0)",friendName.Data(),hashRef).Data());
    ftree->SetAlias("ddYS",TString::Format("%s.DS[1]-AliTPCDistortionFit::Eval(%d,sector,x,y2xAV,z2xAV,1+0)",friendName.Data(),hashRef).Data());
    ftree->SetAlias("ddZS",TString::Format("%s.DS[2]-AliTPCDistortionFit::Eval(%d,sector,x,y2xAV,z2xAV,2+0)",friendName.Data(),hashRef).Data());
    //
    ftree->SetAlias("dXBckg",TString::Format("AliTPCDistortionFit::EvalSectorBckg(%d,%d,fsector,x,z2xAV,0+0,10,0.3)",friendRef,hashRef).Data());
    ftree->SetAlias("dYBckg",TString::Format("AliTPCDistortionFit::EvalSectorBckg(%d,%d,fsector,x,z2xAV,1+0,10,0.3)",friendRef,hashRef).Data());
    ftree->SetAlias("dXCorr",TString::Format("AliTPCDistortionFit::EvalSector(%d,%d,fsector,x,z2xAV,0+0)-dXBckg",friendRef,hashRef).Data());
    ftree->SetAlias("dYCorr",TString::Format("AliTPCDistortionFit::EvalSector(%d,%d,fsector,x,z2xAV,1+0)-dYBckg",friendRef,hashRef).Data());
    //
  }else{
    tree->SetAlias("epsilonfCCM","(8.85418781762e+1)");  //epsilon fC/(V.cm)  https://en.wikipedia.org/wiki/Vacuum_permittivity
    tree->SetAlias("nDiv","(30.+0)");                    // rphi division to extract parameterization
    tree->SetAlias("omegatau","(-0.35+0)");              // omega tau
    tree->SetAlias("EZnorm","(400+0)");                  // EZnom
    tree->SetAlias("padrow","x");                        // pad-row
    tree->SetAlias("sector","bsec");                     // sector
    tree->SetAlias("fsector","bsec+0.5+9.*(y2xAV)/pi");  //"float" sector position using av position in bin
    //tree->SetAlias("csector",TString::Format("bsec+((y2x+0.5)/%d)",ndiv));  //"float" sector using center position of bin
    tree->SetAlias("lx","xAV");                          // local X
    tree->SetAlias("lylx","y2xAV");                      // local y/localX
    tree->SetAlias("lz","z2xAV*xAV");                    // local Z
    tree->SetAlias("r","sqrt(lx**2+(y2xAV*xAV)**2)");    // radius
    tree->SetAlias("gx","r*cos(pi*fsector/9)");          // global x
    tree->SetAlias("gy","r*sin(pi*fsector/9)");          // global y
    tree->SetAlias("side","(-1+2*(bsec<18))");           // side
    TStatToolkit::AddMetadata(tree,"gx.AxisTitle","x_{G} (cm)");
    TStatToolkit::AddMetadata(tree,"gy.AxisTitle","y_{G} (cm)");
    TStatToolkit::AddMetadata(tree,"r.AxisTitle","r (cm)");
    TStatToolkit::AddMetadata(tree,"lz.AxisTitle","z (cm)");
    TStatToolkit::AddMetadata(tree,"lx.AxisTitle","x_{L} (cm)");
    TStatToolkit::AddMetadata(tree,"fsector.AxisTitle","sector position");
    TStatToolkit::AddMetadata(tree,"padrow.AxisTitle","pad-row");
  }
  //  
  TStatToolkit::AddMetadata(tree,prefix+"dXdX.AxisTitle","#Delta_{x}(row+2)-#Delta_{x}(row-2) (cm)");
  TStatToolkit::AddMetadata(tree,prefix+"dYdY.AxisTitle","#Delta_{y}(fsec+0.05)-#Delta_{y}(fsec-0.05) (cm)");

  TStatToolkit::AddMetadata(tree,prefix+"ddX.AxisTitle","#Delta_{x}-#Delta_{xref} (cm)");
  TStatToolkit::AddMetadata(tree,prefix+"ddY.AxisTitle","#Delta_{y}-#Delta_{yref} (cm)");
  TStatToolkit::AddMetadata(tree,prefix+"ddZ.AxisTitle","#Delta_{z}-#Delta_{zref} (cm)");
  TStatToolkit::AddMetadata(tree,prefix+"ddXS.AxisTitle","#Delta_{x}-#Delta_{xref} (cm)");
  TStatToolkit::AddMetadata(tree,prefix+"ddYS.AxisTitle","#Delta_{y}-#Delta_{yref} (cm)");
  TStatToolkit::AddMetadata(tree,prefix+"ddZS.AxisTitle","#Delta_{z}-#Delta_{zref} (cm)");
  TStatToolkit::AddMetadata(tree,prefix+"dX.AxisTitle","#Delta_{x} (cm)");
  TStatToolkit::AddMetadata(tree,prefix+"dXS.AxisTitle","smoothed #Delta_{x} (cm)");
  TStatToolkit::AddMetadata(tree,prefix+"dXC.AxisTitle","cheb. #Delta_{x} (cm)");
  TStatToolkit::AddMetadata(tree,prefix+"dY.AxisTitle","#Delta_{y} (cm)");
  TStatToolkit::AddMetadata(tree,prefix+"dYS.AxisTitle","smoothed #Delta_{y} (cm)");
  TStatToolkit::AddMetadata(tree,prefix+"dYC.AxisTitle","cheb. #Delta_{y} (cm)");
  TStatToolkit::AddMetadata(tree,prefix+"dZ.AxisTitle","#Delta_{z} (cm)");
  TStatToolkit::AddMetadata(tree,prefix+"dZS.AxisTitle","smoothed #Delta_{z} (cm)");
  TStatToolkit::AddMetadata(tree,prefix+"dZC.AxisTitle","cheb. #Delta_{z} (cm)");
  
}









void AliTPCDistortionFit::RegisterMap(TString mapName,  AliTPCChebCorr *map){
  //
  // register map 
  // 
  Int_t hashIndex= mapName.Hash();
  fgkMapNameHash[mapName.Data()]=hashIndex;
  fgkMapHashName[hashIndex]=mapName;
  fgkMapCheb[hashIndex]=map;  
}

void AliTPCDistortionFit::PrintMap(TPRegexp *filter){
  // Loop through the map (Would be much easier in c++11)
  // looping over map with const_iterator
  typedef std::map<string,int>::const_iterator it_type;
  for(it_type iterator = fgkMapNameHash.begin(); iterator != fgkMapNameHash.end(); ++iterator) {
    if (filter){
      if (filter->Match(iterator->first)==0) continue;
    }
    std::cout << iterator->first << " " << iterator->second << "\n";
  }
  return;
}



Double_t  AliTPCDistortionFit::Eval(Int_t hashIndex, int sector, int row, float y2x, float z2x, int dimOut){
  // static interface to Eval function
  const AliTPCChebCorr *cheb=fgkMapCheb[hashIndex];
  if (cheb==0) return 0;
  return cheb->Eval(sector,row, y2x, z2x,dimOut);
}



Double_t AliTPCDistortionFit::EvalEfield(Int_t hashIndex, Int_t hashref, int sector, int row, float y2x, float z2x, int dimOut, Int_t dir, Double_t wt, Double_t dz, Double_t Ez){
  //
  // calculate E field:
  //
  // dimOut  : 0-Er, 1-Erphi,2- Ez
  // dz      : bin width to calculate numberical derivative
  // return value:
  //        numerical derivative return V/cm 
  if (sector>=18) Ez*=-1;
  const AliTPCChebCorr *cheb=fgkMapCheb[hashIndex];
  Double_t localX=cheb->GetPadRowX()[row];
  Double_t dtheta=dz/localX;
  Double_t drdz=0.5*(Eval(hashIndex, sector,row, y2x,z2x+dtheta,0)-Eval(hashIndex, sector,row, y2x,  z2x-dtheta,0))/dz;
  Double_t drdzRef=0.5*(Eval(hashref, sector,row, y2x, z2x+dtheta,0)-Eval(hashref, sector,row, y2x,  z2x-dtheta,0))/dz;
  Double_t drphidz=0.5*(Eval(hashIndex, sector,row, y2x,  z2x+dtheta,1)-Eval(hashIndex, sector,row, y2x, z2x-dtheta,1))/dz;
  Double_t drphidzRef=0.5*(Eval(hashref, sector,row, y2x,  z2x+dtheta,1)-Eval(hashref, sector,row, y2x,  z2x-dtheta,1))/dz;
  drdz-=drdzRef;
  drphidz-=drphidzRef;
  //
  Double_t c0=1./(1+wt*wt); // ignoring "tensor terms"
  Double_t c1=dir*wt/c0;            // ignoring "tensor terms"
  //
  Double_t erez=   (drdz-wt*drphidz)/(c0+c1*c1/c0);
  Double_t erphiez=(drphidz+wt*drdz)/(c0+c1*c1/c0);
  if (dimOut==0) return Ez*(drdz-wt*drphidz)/(c0+c1*c1/c0);
  if (dimOut==1) return Ez*(drphidz+wt*drdz)/(c0+c1*c1/c0);
  // 
  // tree.Draw("AliTPCDistortionFit::EvalEfield(hashIndex,hashRefIndex,  sector, padrow, y2xAV, z2xAV, 0, polarity, omegatau,2.5,EZnorm):padrow:fsector>>his(540,0,36,160,0,160)","fitOK&&z2x==1","profcolz")
  // tree.Draw("AliTPCDistortionFit::EvalEfield(hashIndex,hashRefIndex,  sector, padrow, y2xAV, z2xAV, 1, polarity, omegatau,2.5,EZnorm):padrow:fsector>>his(540,0,36,160,0,160)","fitOK&&z2x==1","profcolz")
  // tree.Draw("AliTPCDistortionFit::EvalEfield(hashIndex,hashRefIndex,  sector, padrow, y2xAV, z2xAV, 1, polarity, omegatau,2.5,EZnorm):fsector:padrow","fitOK&&z2x==1","colz");
  return 0;
}
Double_t AliTPCDistortionFit::EvalRho(Int_t hashIndex, Int_t hashref, int sector, int row, float y2x, float z, Int_t dir, Double_t wt, Double_t dz,  Double_t dr, Double_t drphi,  Double_t Ez){
  // retun numerical derivative V/(cm^2)
  // to get density  multiply by epsilon ( C/(Vm) )
  // second order correction needed to transform dr, drphi 
  const Double_t sectorBin=1/15.;
  const Double_t phiBin=sectorBin*TMath::Pi()/9.; 
  Int_t row0=TMath::Max(Int_t(row-dr),0);
  Int_t row1=TMath::Min(Int_t(row+dr),158);
  const AliTPCChebCorr *cheb=fgkMapCheb[hashIndex];
  Double_t drused=cheb->GetPadRowX()[row1]-cheb->GetPadRowX()[row0];
  Double_t r=cheb->GetPadRowX()[row];
  Double_t fsector=(sector+0.5+y2x*9/TMath::Pi());   // float sector
  Float_t  drphiUsed=drphi;
  if (drphiUsed<=0) drphiUsed=r*phiBin;
  Double_t fsector0=fsector;
  Double_t fsector1=fsector+(9*drphiUsed/r)/TMath::Pi();

  Int_t    sector0=TMath::Nint(fsector0-0.5);
  Int_t    sector1=TMath::Nint(fsector1-0.5)%36;
  Double_t y2x0= TMath::Pi()*(fsector0-(sector0+0.5))/9.;
  Double_t y2x1= TMath::Pi()*(fsector1-(sector1+0.5))/9.;
  Double_t derdr= (EvalEfield(hashIndex,hashref, sector,row1, y2x,z,0,dir,wt,dz,Ez)-EvalEfield(hashIndex,hashref, sector,row0, y2x,z,0,dir,wt,dz,Ez))/drused;  
  Double_t derphidr= (EvalEfield(hashIndex,hashref, sector1,row, y2x1,z,1,dir,wt,dz,Ez)-EvalEfield(hashIndex,hashref, sector0,row, y2x0,z,1,dir,wt,dz,Ez))/drphiUsed;
  // (V/cm)/cm
  //   tree.Draw("EvalRho(hashIndex,hashRefIndex,  sector, padrow, y2xAV, z2xAV, polarity, omegatau,5,5,5,EZnorm):fsector:padrow","fitOK&&z2x==1","colz");
  return derdr+derphidr;
}

Double_t AliTPCDistortionFit::EvalSectorBckg(Int_t hashIndex, Int_t hashRef, Double_t fsector, int row, float z2x, int dimOut, Double_t dRow, Double_t dSec){
  //
  // evaluate sector backround caluclating weighted mean for 4 points around point of interest
  //
  AliTPCROC *roc = AliTPCROC::Instance();
  Int_t nRowsAll=roc->GetNRows(0)+roc->GetNRows(36);
  Double_t fsectorM=TMath::Nint(fsector);
  Double_t r0=(row<Int_t(roc->GetNRows(0))) ?  roc->GetPadRowRadii(0, row): roc->GetPadRowRadii(36, row-roc->GetNRows(0));
  Double_t rphi0=(fsector-fsectorM)*r0;  
  Double_t values[4], distance2[4];
  for (Int_t i=0; i<4; i++){
    Double_t rowC=(i<2)? row-dRow:row+dRow;
    Double_t secC=(i%2==0)? fsectorM-dSec:fsectorM+dSec;
    if (rowC<0) rowC=0;
    if (rowC>=nRowsAll) rowC=nRowsAll-1;
    Double_t rC=(rowC<Int_t(roc->GetNRows(0))) ?  roc->GetPadRowRadii(0, rowC): roc->GetPadRowRadii(36, rowC-roc->GetNRows(0));
    Double_t rphiC=(secC-fsectorM)*r0;
    values[i]=AliTPCDistortionFit::EvalSector(hashIndex, hashRef, secC, rowC, z2x,dimOut);
    distance2[i]=((rC-r0)*(rC-r0)+(rphi0-rphiC)*(rphi0-rphiC));    
  }
  Double_t sum=0, sumw=0;
  for (Int_t i=0; i<4; i++){
    sum+=values[i]/(1.+distance2[i]);
    sumw+=1/(1.+distance2[i]);
  }
  return sum/sumw;
}

Double_t AliTPCDistortionFit::EvalSector(Int_t hashIndex, Int_t hashRef, Double_t fsector, int row, float z2x, int dimOut){
  return AliTPCDistortionFit::EvalSector(hashIndex, fsector, row, z2x,dimOut,0,0)-AliTPCDistortionFit::EvalSector(hashRef, fsector, row, z2x,dimOut,0,0);
}

Double_t AliTPCDistortionFit::EvalSector(Int_t hashIndex, Double_t fsector, int row, float z2x, int dimOut, int interpol, Double_t mndiv){  
  //
  //const Double_t mndiv=1./30.;
  if (fsector<0) fsector+=18;
  if (fsector>36) fsector-=18;
  if (row<0) row=0;
  if (row>158) row=158;

  const Double_t dphi=TMath::Pi()/9.;
  Int_t sector=Int_t(fsector);
  Double_t dsector=fsector-sector;
  if (interpol==0 || (dsector>mndiv && dsector<1-mndiv)){
    if (sector<0) sector=0;
    if (sector>35) sector=35;
    Double_t y2x=dphi*(fsector-(sector+0.5));
    return Eval(hashIndex,sector,row, y2x, z2x,dimOut);
  }
  
  Double_t sector0,sector1, y2x0,y2x1;
  Double_t dbin=0;
  if (dsector<mndiv){
    sector0=(sector+35)%36;
    sector1=sector;
    y2x0=dphi*((1.-mndiv)-0.5);
    y2x1=dphi*(mndiv-0.5);
    dbin=(fsector-sector);
  }else{
    sector0=sector;
    sector1=(sector+1)%36;
    y2x0=dphi*((1.-mndiv)-0.5);
    y2x1=dphi*(mndiv-0.5);    
    dbin=(fsector-sector-1);
  }
  if (interpol==1){
    Double_t val0=Eval(hashIndex,sector0,row, y2x0, z2x,dimOut);  // Eval(-mndiv)
    Double_t val1=Eval(hashIndex,sector1,row, y2x1, z2x,dimOut);  // Eval(mndiv)    
    Double_t linear=((dbin+mndiv)*val1+(mndiv-dbin)*val0)/(2*mndiv);
    return linear;
  }
  if (interpol==2){
    Double_t val1=Eval(hashIndex,sector0,row, y2x0, z2x,dimOut);  // Eval(-1*mndiv)
    Double_t val2=Eval(hashIndex,sector1,row, y2x1, z2x,dimOut);  // Eval(-1*mndiv)
    //
    Double_t dval1=(Eval(hashIndex,sector0,row, y2x0+0.05*dphi*mndiv, z2x,dimOut)-Eval(hashIndex,sector0,row, y2x0-0.05*dphi*mndiv, z2x,dimOut))/(0.1*mndiv);
    Double_t dval2=(Eval(hashIndex,sector1,row, y2x1+0.05*dphi*mndiv, z2x,dimOut)-Eval(hashIndex,sector1,row, y2x1-0.05*dphi*mndiv, z2x,dimOut))/(0.1*mndiv);
    //   q=(1-t)y_{1}+ty_{2}+t(1-t)(a(1-t)+bt)}    https://en.wikipedia.org/wiki/Spline_interpolation
    //
    Double_t t= (dbin+mndiv)/(2.*mndiv);  // working in sector unit
    Double_t a= dval1*(2.*mndiv)-(val2-val1);
    Double_t b= -dval2*(2.*mndiv)+(val2-val1);
    Double_t value=(1-t)*val1+t*val2+t*(1-t)*(a*(1-t)+b*t);
    return value;     
  }
  return 0;
}




Double_t AliTPCDistortionFit::EvalEfieldSector(Int_t hashIndex, Int_t hashref, Float_t fsector, int row, float z2x, int dimOut, Int_t polarity, Double_t wt, Double_t dz, Double_t Ez, Int_t interpolationType, Double_t mBinSize){
  //
  // hash          : map hash  
  // hashref       : reference map hash
  // 
  // fsector       : relative sector position 
  // z2x           : z/x coordinate
  // dimOut        : 0-Er, 1-Erphi,2- Ez
  // polarity           : direction of magnetic field
  // wt            : omega tau
  // dz            : bin width to calculate numberical derivative
  // interpolationType : 0-constant 1.)  linear 2.) cubic spline
  // mBinSize      : 1/binSize used for map extraction
  // return value:
  //        numerical derivative return V/cm 
  //        use the same formula as in  AliTPCSpaceCharge3D::GetCorrection
  //
  if (fsector>=18) {
    Ez*=-1;
  }
  const AliTPCChebCorr *cheb=fgkMapCheb[hashIndex];
  Double_t localX=cheb->GetPadRowX()[row];
  Double_t dtheta=dz/localX;  // convert dz->dTheta used in parametrization
  //
  Double_t drdz=0.5*(EvalSector(hashIndex, fsector,row, z2x+dtheta,0,interpolationType, mBinSize)- \
		     EvalSector(hashIndex, fsector,row,   z2x-dtheta,0,interpolationType, mBinSize))/dz;
  Double_t drdzRef=0.5*(EvalSector(hashref, fsector,row,  z2x+dtheta,0,interpolationType, mBinSize)- \
			EvalSector(hashref, fsector,row,   z2x-dtheta,0,interpolationType, mBinSize))/dz;
  Double_t drphidz=0.5*(EvalSector(hashIndex, fsector,row,   z2x+dtheta,1,interpolationType, mBinSize)-\
			EvalSector(hashIndex, fsector,row,  z2x-dtheta,1,interpolationType, mBinSize))/dz;
  Double_t drphidzRef=0.5*(EvalSector(hashref, fsector,row,   z2x+dtheta,1,interpolationType, mBinSize)-\
			   EvalSector(hashref, fsector,row, z2x-dtheta,1,interpolationType, mBinSize))/dz;
  //
  Double_t c0=1./(1+wt*wt);         // Decompose E and ExB component
  Double_t c1=polarity*wt/c0;       // 
  Double_t erez=  drdz-drphidz*c1/c0;
  Double_t erphiez=drphidz+drdz*c1/c0;
  if (dimOut==0) return Ez*(erez);
  if (dimOut==1) return Ez*(erphiez);
  // 
  // tree.Draw("AliTPCDistortionFit::EvalEfieldSector(hashIndex,hashRefIndex,  fsector, padrow, z2xAV, 1, polarity, omegatau,5,EZnorm,1,1./nDiv):fsector:padrow","fitOK&&z2x==1","colz");

}


Double_t AliTPCDistortionFit::EvalRhoSector(Int_t hashIndex, Int_t hashref, Float_t fsector, int row, float z2x, Int_t polarity, Double_t wt, Double_t dz,  Double_t dr, Double_t drphi,  Double_t Ez, Int_t interpolationType, Double_t mBinSize){
  // retun numerical derivative V/(cm^2)
  // to get density  multiply by epsilon ( C/(Vm) )
  // second order correction needed to transform dr, drphi 
  Int_t row0=TMath::Max(Int_t(row-dr),0);
  Int_t row1=TMath::Min(Int_t(row+dr),158);
  const AliTPCChebCorr *cheb=fgkMapCheb[hashIndex];
  Double_t drused=cheb->GetPadRowX()[row1]-cheb->GetPadRowX()[row0];
  Double_t r=cheb->GetPadRowX()[row];
  Double_t dsector=9.*(drphi/r)/TMath::Pi();
  Double_t fsector0=fsector-dsector;
  Double_t fsector1=fsector+dsector;
  if (fsector0<0)   fsector0+=36;  //!!!!! should treat A side C side 
  if (fsector1>=36) fsector1-=36;
  //
  Double_t derdr= (EvalEfieldSector(hashIndex,hashref, fsector,row1,z2x,0,polarity,wt,dz,Ez,interpolationType,mBinSize)-EvalEfieldSector(hashIndex,hashref, fsector,row0,z2x,0,polarity,wt,dz,Ez,interpolationType,mBinSize))/drused;  
  Double_t derphidr= (EvalEfieldSector(hashIndex,hashref, fsector1,row,z2x,1,polarity,wt,dz,Ez,interpolationType,mBinSize)-EvalEfieldSector(hashIndex,hashref, fsector0,row,z2x,1,polarity,wt,dz,Ez,interpolationType,mBinSize))/drphi;

  // (V/cm)/cm
  //   tree.Draw("AliTPCDistortionFit::EvalRhoSector(hashIndex,hashRefIndex,  fsector, padrow, z2xAV, polarity, omegatau,5,5,5,EZnorm,2,1./nDiv):fsector:padrow","fitOK&&z2x==1","colz");
  return derdr+derphidr;

}

Int_t AliTPCDistortionFit::RegisterFitters(){
  //
  //
  //  
  TF1 *likeGausCachy = new TF1("likeGausCachy", AliTMinuitToolkit::GaussCachyLogLike,-10,10,2);
  likeGausCachy->SetParameters(0.8,1);
  //  1 line function
  TF1 *funLine1 = new TF1("funLine1",AliTPCDistortionFit::LineFieldLocal,-20,20,7);
  AliTMinuitToolkit * fitterLine1 = new AliTMinuitToolkit();
  fitterLine1->SetVerbose(0x1); fitterLine1->SetFitFunction(funLine1,kTRUE);
  //fitterLine1->SetLogLikelihoodFunction(likeGausCachy);
  TMatrixD initPar1(7,4),   param1(7,1);       // initial parameters of 1D fits  - see functionAliTPCDistortionFit::LineFieldLocal
  initPar1(0,0)=1;   initPar1(0,1)=1;    initPar1(0,2)=0;    initPar1(0,3)=0;            // q
  initPar1(1,0)=120; initPar1(1,1)=10;   initPar1(1,2)=80;   initPar1(1,3)=140;          // r
  initPar1(2,0)=0;   initPar1(2,1)=1;    initPar1(2,2)=0;    initPar1(2,3)=0;            // rphi
  initPar1(3,0)=1;   initPar1(3,1)=0;    initPar1(3,2)=0.01; initPar1(3,3)=5;            // scale distance
  initPar1(4,0)=0.35;initPar1(4,1)=0.1;  initPar1(4,2)=-0.5; initPar1(4,3)=0.5;          // wt   
  initPar1(5,0)=40;  initPar1(5,1)=0;    initPar1(5,2)=0;    initPar1(5,3)=0.0;          // symetry plane
  initPar1(6,0)=400; initPar1(6,1)=0;    initPar1(6,2)=0;    initPar1(6,3)=0.0;          // Ez 
  for (Int_t ipar=0; ipar<7; ipar++) {funLine1->SetParameter(ipar, initPar1(ipar,0)); param1(ipar,0)=initPar1(ipar,0);}
  fitterLine1->SetInitialParam(&initPar1);
  AliTMinuitToolkit::SetPredefinedFitter("fitterLine1",fitterLine1);  
  //
  //  N line function
  TF1 *funLineN = new TF1("funLineN",AliTPCDistortionFit::NLinesFieldLocal,-20,20,21);
  AliTMinuitToolkit * fitterLineN = new AliTMinuitToolkit();
  fitterLineN->SetVerbose(0x1); fitterLineN->SetFitFunction(funLineN,kTRUE);
  TMatrixD initParN(21,4),   paramN(21,1);            // initial parameters of 1D fits  - see functionAliTPCDistortionFit::LineFieldLocal
  initParN(0,0)=0.1;   initParN(0,1)=0;    initParN(0,2)=0.01; initParN(0,3)=5;            // scale distance
  initParN(1,0)=0.35;   initParN(1,1)=0.0;  initParN(1,2)=-0.5; initParN(1,3)=0.5;          // wt   
  initParN(2,0)=400;   initParN(2,1)=0;    initParN(2,2)=0;    initParN(2,3)=0.0;          // Ez 
  initParN(3,0)=40;    initParN(3,1)=0;    initParN(3,2)=0;    initParN(3,3)=0.0;          /// symetry plane
  initParN(4,0)=1;     initParN(4,1)=0;    initParN(4,2)=0;    initParN(4,3)=0.0;          // nlines
  for (Int_t ipar=4; ipar<21; ipar++) {initParN(ipar,1)=1;}
  for (Int_t ipar=0; ipar<21; ipar++) {funLineN->SetParameter(ipar, initParN(ipar,0));}
  fitterLineN->SetInitialParam(&initParN);
  AliTMinuitToolkit::SetPredefinedFitter("fitterLineN",fitterLineN);  
  //

}

void AliTPCDistortionFit::MakeFitExample1(Int_t run, const char * chinput, const char * ocdbPath){
  /*
    Int_t run=245683; const char * chinput="/data/alien/alice/data/2015/LHC15o/000245683/cpass0_pass1/ResidualMerge/TPCSPCalibration/1448920523_1448922567_000245683/voxelResTree.root", ocdbPath="local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/"
    AliTPCDistortionFit::MakeFitExample1(run,chinput);
  */
  TTreeSRedirector *pcstream = new TTreeSRedirector("makeFit1D.root","recreate");
  // 1 line fitter  - working for sectors with one hotspots
  //
  AliTPCDistortionFit::RegisterFitters();
  AliTPCDistortionFit::LoadDistortionMaps(run, ocdbPath);
  TString fName =AliTPCDistortionFit::LoadDistortionTree(chinput);
  
  TObjString oName(fName.Data());
  Int_t hashID=fName.Hash();
  AliTMinuitToolkit * fitter = AliTMinuitToolkit::GetPredefinedFitter("fitterLine1");
  TGraph * lumiGraphMap = AliLumiTools::GetLumiGraph(AliTPCRecoParam::kCorrMapNoScaling,run);
  const AliTPCChebCorr *cheb = AliTPCDistortionFit::GetCheb(fName.Data());
  AliLHCData* fLHCDataBckg = (AliLHCData*)((AliCDBEntry*)(AliCDBManager::Instance()->Get("GRP/GRP/LHCData")))->GetObject();
  Double_t polarity=0;
  if (cheb->GetFieldType()== AliTPCChebCorr::kFieldPos) polarity=1;
  if (cheb->GetFieldType()== AliTPCChebCorr::kFieldNeg) polarity=-1;

  TMatrixD initPar1(9,4),   param1(9,1);       // initial parameters of 1D fits  - see functionAliTPCDistortionFit::LineFieldLocal
  initPar1(0,0)=10;   initPar1(0,1)=1;    initPar1(0,2)=0;    initPar1(0,3)=0;            // q
  initPar1(1,0)=110; initPar1(1,1)=5;   initPar1(1,2)=80;   initPar1(1,3)=140;          // r
  initPar1(2,0)=0;   initPar1(2,1)=1;    initPar1(2,2)=-4;   initPar1(2,3)=4;            // rphi
  initPar1(3,0)=0.1; initPar1(3,1)=0;    initPar1(3,2)=0.01; initPar1(3,3)=5;            // scale distance
  initPar1(4,0)=polarity*0.35; initPar1(4,1)=0;    initPar1(4,2)=0;    initPar1(4,3)=0.0;          // wt   
  initPar1(5,0)=0;   initPar1(5,1)=0;    initPar1(5,2)=0;    initPar1(5,3)=80;            // symetry plane
  initPar1(6,0)=400; initPar1(6,1)=0;    initPar1(6,2)=0;    initPar1(6,3)=0.0;          // Ez 
  fitter->SetInitialParam(&initPar1);
  //
  fitter->SetVerbose(AliTMinuitToolkit::kPrintAll); 
  AliTPCDistortionFit::fgkDistortionTree->SetAlias("fx","(((Entry$%2)==0)+0)");
  AliTPCDistortionFit::fgkDistortionTree->SetAlias("fy","(((Entry$%2)!=0)+0)");
  AliTPCDistortionFit::fgkDistortionTree->SetAlias("fitCut",TString::Format("lx<150&&%s.fitOK&&%s.dispOK",fName.Data(),fName.Data()));  
  Int_t fsectors[16]={2,   4,  6,   7,    9,   10,  11,  13,   16,  19,  20,   29, 30,   31,   35};
  Int_t frows[16]   ={116, 104,104, 107,  114, 122, 130, 110,  95,  110, 117,  100, 110,  111,  122};
  Double_t chi2=0;
  Int_t flag= Int_t(AliTMinuitToolkit::kStreamFcnPoint);
  Int_t npar=4;
  TLatex latex;
  latex.SetTextSize(0.035);

  for (Int_t isec=0; isec<15; isec++){
    Int_t jsec=fsectors[isec];
    initPar1(1,0)=frows[isec];
    AliTPCDistortionFit::fgkDistortionTree->SetAlias("sectorCut",TString::Format("abs(fsector-%d)<0.25",jsec).Data());
    AliTPCDistortionFit::fgkDistortionTree->SetAlias("dYSector",TString::Format("lx*pi*(fsector-%d)/9+%s.dYS",jsec,fName.Data()).Data());
    AliTPCDistortionFit::fgkDistortionTree->SetAlias("lXCorr",TString::Format("lx+%s.dXS",fName.Data()));
    AliTPCDistortionFit::fgkDistortionTree->SetAlias("dXRun",TString::Format("%s.dXCorr",fName.Data()));
    AliTPCDistortionFit::fgkDistortionTree->SetAlias("dYRun",TString::Format("%s.dYCorr",fName.Data()));
    {
      //fitter->FillFitter(AliTPCDistortionFit::fgkDistortionTree,"dXRun:0.0025", "4:lXCorr:dYSector:lz","sectorCut&&fitCut",0,10000000, kTRUE);
      //fitter->FillFitter(AliTPCDistortionFit::fgkDistortionTree,"dYRun:1.", "5:lXCorr:dYSector:lz","sectorCut&&fitCut",0,10000000, kFALSE);
      fitter->FillFitter(AliTPCDistortionFit::fgkDistortionTree,"dYRun:1.", "5:lXCorr:dYSector:lz","sectorCut&&fitCut",0,10000000, kTRUE);
      fitter->Fit(); printf("Fit++\t%f\n",sqrt(fitter->GetChisquare()/fitter->GetPoints()->GetNrows()));
      fitter->SetVerbose(0);
      fitter->Bootstrap(10,0,0);
      if ((*(fitter->GetRMSEstimator()))[1]>1.) {
	fitter->Fit(); printf("Fit++\t%f\n",sqrt(fitter->GetChisquare()/fitter->GetPoints()->GetNrows()));
      }
    }  
    fitter->SetStreamer("testFitLine1.root");
    fitter->SetVerbose( AliTMinuitToolkit::kPrintAll| AliTMinuitToolkit::kStreamFcn);
    fitter->FillFitter(AliTPCDistortionFit::fgkDistortionTree,"dXRun:1", "4:lXCorr:dYSector:lz","sectorCut&&fitCut",0,10000000);
    AliTMinuitToolkit::FitterFCN(npar,0,chi2,(Double_t*)(fitter->GetParameters()->GetMatrixArray()), flag);
    fitter->FillFitter(AliTPCDistortionFit::fgkDistortionTree,"dYRun:1", "5:lXCorr:dYSector:lz","sectorCut&&fitCut",0,10000000);
    AliTMinuitToolkit::FitterFCN(npar,0,chi2,(Double_t*)(fitter->GetParameters()->GetMatrixArray()), flag);
    ((*(fitter->GetStreamer())).GetFile())->Flush();
    TTree * t = ((*(fitter->GetStreamer()))<<"fcnDebugPoint").GetTree();
    t->ResetBranchAddresses();
    t->SetMarkerStyle(25); t->SetMarkerSize(0.5);    
    t->SetAlias("xlocal","x.fElements[1]");
    t->SetAlias("yedge","100*x.fElements[2]/xlocal");
    TCanvas * canvasFit = new TCanvas("canvasFit","canvasFit",1300,900);
    canvasFit->Divide(3,2);
    Double_t maxY,minY;
    // Retrieve the stat box
    gStyle->SetTitleOffset(0.92,"X");
    gStyle->SetTitleOffset(1.2,"Y");
    gStyle->SetTitleOffset(1,"Z");
    for (Int_t ipad=0; ipad<6;ipad++){
      Int_t icol=ipad%3;
      Int_t irow=ipad/3;
      TVirtualPad *vpad=canvasFit->cd(ipad+1);
      vpad->SetPad(icol/3., irow/2., (icol+1)/3., (irow+1)/2.);
      if (irow==1) {
	vpad->SetBottomMargin(0);
	vpad->SetTopMargin(0);
      }else{
	vpad->SetTopMargin(0);
	vpad->SetBottomMargin(0.15);
      }
      if (icol<2) {
	vpad->SetRightMargin(0);
      }else {
	vpad->SetRightMargin(0.2);
      }
      if (icol>0) {
	vpad->SetLeftMargin(0); 
      }else{
	vpad->SetLeftMargin(0.2); 
      }

      vpad->Draw();
      canvasFit->cd();
    }
    TVectorD rmsEstimator(*(fitter->GetRMSEstimator()));
    {       
      canvasFit->cd(4);
      t->Draw("val:xlocal:yedge","x.fElements[0]==4&&abs(abs(x.fElements[3]/xlocal)-0.1)<0.1","colz");
      maxY=TMath::Nint(TMath::MaxElement(t->GetSelectedRows(), t->GetV1())+0.5)+0.5;
      minY=TMath::Nint(TMath::MinElement(t->GetSelectedRows(), t->GetV1())-0.5)-0.5;
      t->GetHistogram()->GetYaxis()->SetLimits(minY,maxY);
      t->GetHistogram()->GetYaxis()->SetTitle("#Delta_{R} (cm)"); 
      latex.DrawLatexNDC(0.25,0.85,"Data");
      canvasFit->cd(5);
      t->Draw("fun:xlocal:yedge","x.fElements[0]==4&&abs(abs(x.fElements[3]/xlocal)-0.1)<0.1","colz");
      t->GetHistogram()->GetYaxis()->SetLimits(minY,maxY);
      latex.DrawLatexNDC(0.15,0.85,"Line charge fit");
      canvasFit->cd(6);
      t->Draw("val-fun:xlocal:yedge","x.fElements[0]==4&&abs(abs(x.fElements[3]/xlocal)-0.1)<0.1","colz");
      t->GetHistogram()->GetYaxis()->SetLimits(minY,maxY);
      t->GetHistogram()->GetZaxis()->SetTitle("100 #Delta_{edge}/x_{L} (unit)");
      latex.DrawLatexNDC(0.15,0.85,"Data-Fit");
      canvasFit->cd(1);
      t->Draw("val:xlocal:yedge","x.fElements[0]==5&&abs(abs(x.fElements[3]/xlocal)-0.1)<0.1","colz");
      maxY=TMath::Nint(TMath::MaxElement(t->GetSelectedRows(), t->GetV1())+0.5)+0.5;
      minY=TMath::Nint(TMath::MinElement(t->GetSelectedRows(), t->GetV1())-0.5)-0.5;
      t->GetHistogram()->GetYaxis()->SetLimits(minY,maxY);
      t->GetHistogram()->GetYaxis()->SetTitle("#Delta_{R#phi} (cm)");
      t->GetHistogram()->GetXaxis()->SetTitle("local X (cm)");
      canvasFit->cd(2);
      t->Draw("fun:xlocal:yedge","x.fElements[0]==5&&abs(abs(x.fElements[3]/xlocal)-0.1)<0.1","colz");
      t->GetHistogram()->GetYaxis()->SetLimits(minY,maxY);
      t->GetHistogram()->GetXaxis()->SetTitle("local X (cm)");
      canvasFit->cd(3);
      t->Draw("val-fun:xlocal:yedge","x.fElements[0]==5&&abs(abs(x.fElements[3]/xlocal)-0.1)<0.1","colz");
      t->GetHistogram()->GetYaxis()->SetLimits(minY,maxY);
      t->GetHistogram()->GetZaxis()->SetTitle("100 #Delta_{edge}/x_{L} (unit)");
      t->GetHistogram()->GetXaxis()->SetTitle("local X (cm)");
      //t->GetHistogram()->GetListOfFunctions()->Add(fitter->GetFormula());
      //t->GetHistogram()->Draw("colz");
      latex.DrawLatexNDC(0.05,0.9,fName.Data());
      latex.DrawLatexNDC(0.05,0.85,TString::Format("Sector=%d",jsec).Data());
      latex.DrawLatexNDC(0.05,0.8,TString::Format("Q=%2.2f#pm%2.2f",(*(fitter->GetParameters()))[0],rmsEstimator[0]).Data());
      latex.DrawLatexNDC(0.05,0.75,TString::Format("R=%2.2f#pm%2.2f (cm)",(*(fitter->GetParameters()))[1],rmsEstimator[1]).Data());
      latex.DrawLatexNDC(0.05,0.70,TString::Format("R#phi=%2.2f#pm%2.2f (cm)",(*(fitter->GetParameters()))[2],rmsEstimator[2]).Data());
      //      latex.DrawLatexNDC(0.05,0.65,TString::Format("scale=%2.2f#pm%2.2f (cm)",(*(fitter->GetParameters()))[3],rmsEstimator[3]).Data());
      latex.DrawLatexNDC(0.05,0.65,TString::Format("#omega#tau=%2.2f#pm%2.2f",(*(fitter->GetParameters()))[4],rmsEstimator[4]).Data());
    }
    canvasFit->SaveAs(TString::Format("canvasDistortionFitLine1_Sec%d.png",fsectors[isec]).Data());
    //const AliTPCChebCorr *cheb = AliTPCDistortionFit::GetCheb(fName.Data());
    TVectorD *param= (TVectorD *)fitter->GetParameters();
    TMatrixD *covar=(TMatrixD*)fitter->GetCovarianceMatrix();
    TVectorD *rms= (TVectorD *)fitter->GetRMSEstimator();
    Double_t chi2=fitter->GetChisquare();
    Double_t tC=cheb->GetTimeStampCenter();
    Double_t tB=cheb->GetTimeStampStart();
    Double_t tE= cheb->GetTimeStampEnd();
    Double_t lumi = cheb->GetLuminosityCOG(lumiGraphMap,tB,tE);
    //
    const Int_t nBGs = AliLHCData::kNBGs;
    Double_t bckg[5]={0};
    for (Int_t ibg=0;ibg<nBGs;ibg++){
      Int_t counter=0;
      for (Int_t it=tB;it<tE; it++){ // get mean bckg
	bckg[ibg] += fLHCDataBckg->GetBckgAlice(ibg,Double_t(it));  // function returns differnt values in case of calling with interger and case of double invocation
	counter++;
      }
      bckg[ibg]/=counter;
    }
 
    //
    pcstream->GetFile()->cd();
    {
      (*pcstream)<<"fit1D"<<
	"run="<<run<<
	"name.="<<&oName<<
	"hashID="<<hashID<<
	"sec="<<jsec<<
	"par.="<<param<<
	"cov.="<<covar<<
	"rms.="<<rms<<
	"chi2="<<chi2<<
	"tC="<<tC<<
	"tB="<<tB<<
	"tE="<<tE<<
	"lumi="<<lumi;
      for (Int_t ibg=0;ibg<5;ibg++){
	(*pcstream)<<"fit1D"<<
	  TString::Format("bckg%d=",ibg)<<bckg[ibg]<<
	  "\n";      
      }
    }
  } 
  (((*pcstream)<<"fit1D").GetTree())->Write();
  pcstream->GetFile()->Flush();
  delete pcstream;
}



void MakeFitExample4(Int_t run=245683, const char * chinput="/data/alien/alice/data/2015/LHC15o/000245683/cpass0_pass1/ResidualMerge/TPCSPCalibration/1448920523_1448922567_000245683/voxelResTree.root"){
  //
  // N(4) line fitter  - working for sectors with 0-N  hotspots
  //
  AliTPCDistortionFit::RegisterFitters();
  AliTPCDistortionFit::LoadDistortionMaps(run, "local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/");
  TString fName =AliTPCDistortionFit::LoadDistortionTree(chinput);
  AliTMinuitToolkit * fitterN = AliTMinuitToolkit::GetPredefinedFitter("fitterLineN");
   
  TMatrixD initParN(21,4),   paramN(21,1);       // initial parameters of 1D fits  - see functionAliTPCDistortionFit::LineFieldLocal
  initParN(0,0)=1;   initParN(0,1)=0;    initParN(0,2)=0.01; initParN(0,3)=5;            // scale distance
  initParN(1,0)=0.3; initParN(1,1)=0;    initParN(1,2)=-0.5; initParN(1,3)=0.5;          // wt   
  initParN(2,0)=400; initParN(2,1)=0;    initParN(2,2)=0;    initParN(2,3)=0.0;          // Ez 
  initParN(3,0)=4;   initParN(3,1)=0;    initParN(3,2)=0;    initParN(3,3)=0.0;          // Ez 
  Double_t deltaR=(134-85)/8.;
  Double_t r0=85+deltaR;
  for (Int_t ir=0; ir<4; ir++){
    Int_t offset=4+ir*3;   
    Double_t R=r0+deltaR*ir*2;
    initParN(offset+0,0)=3;     initParN(offset+0,1)=1;    initParN(offset+0,2)=0;            initParN(offset+0,3)=1000.0;     // q
    initParN(offset+1,0)=R;     initParN(offset+1,1)=5;    initParN(offset+1,2)=R-deltaR+1;   initParN(offset+1,3)=R+deltaR-1;    // q
    initParN(offset+2,0)=0;     initParN(offset+2,1)=0;    initParN(offset+2,2)=-2;           initParN(offset+2,3)=2;            // q
  }
  fitterN->SetInitialParam(&initParN);
  //
  fitterN->SetVerbose(AliTMinuitToolkit::kPrintAll); 
  AliTPCDistortionFit::fgkDistortionTree->SetAlias("fx","(((Entry$%2)==0)+0)");
  AliTPCDistortionFit::fgkDistortionTree->SetAlias("fy","(((Entry$%2)!=0)+0)");
  AliTPCDistortionFit::fgkDistortionTree->SetAlias("fitCut",TString::Format("lx<150&&%s.fitOK&&%s.dispOK",fName.Data(),fName.Data()));  
  Int_t fsectors[12]={2,4,6,7,9,10, 16,20,29, 30,31,35};
  Double_t chi2=0;
  Int_t flag= Int_t(AliTMinuitToolkit::kStreamFcnPoint);
  Int_t npar=4;
  TLatex latex;
  latex.SetTextSize(0.035);

  for (Int_t isec=0; isec<12; isec++){
    Int_t jsec=fsectors[isec];
    AliTPCDistortionFit::fgkDistortionTree->SetAlias("sectorCut",TString::Format("abs(fsector-%d)<0.3",jsec).Data());
    AliTPCDistortionFit::fgkDistortionTree->SetAlias("dYSector",TString::Format("lx*pi*(fsector-%d)/9+%s.dYS",jsec,fName.Data()).Data());
    AliTPCDistortionFit::fgkDistortionTree->SetAlias("lXCorr",TString::Format("lx+%s.dXS",fName.Data()));
    AliTPCDistortionFit::fgkDistortionTree->SetAlias("dXRun",TString::Format("%s.ddXS",fName.Data()));
    AliTPCDistortionFit::fgkDistortionTree->SetAlias("dYRun",TString::Format("%s.ddYS",fName.Data()));
    {
      //      fitterN->FillFitter(AliTPCDistortionFit::fgkDistortionTree,"fx*dXRun+fy*dYRun:1+fy", "2+(fy>0):lXCorr:dYSector:lz","sectorCut&&fitCut",0,10000000, kFALSE);
      //fitterN->FillFitter(AliTPCDistortionFit::fgkDistortionTree,"dXRun:0.5", "2:lXCorr:dYSector:lz","sectorCut&&fitCut",0,10000000, kTRUE);
      //      fitterN->FillFitter(AliTPCDistortionFit::fgkDistortionTree,"dYRun:1.", "3:lXCorr:dYSector:lz","sectorCut&&fitCut",0,10000000, kFALSE);
      fitterN->FillFitter(AliTPCDistortionFit::fgkDistortionTree,"dYRun:1.", "3:lXCorr:dYSector:lz","sectorCut&&fitCut",0,10000000, kTRUE);
      fitterN->Fit(); printf("Fit++\t%f\n",sqrt(fitterN->GetChisquare()/fitterN->GetPoints()->GetNrows()));
      fitterN->SetVerbose(0);
      fitterN->Bootstrap(10,0,0);
    }  
    fitterN->SetStreamer("testFitLine1.root");
    fitterN->SetVerbose( AliTMinuitToolkit::kPrintAll| AliTMinuitToolkit::kStreamFcn);
    fitterN->FillFitter(AliTPCDistortionFit::fgkDistortionTree,"dXRun:1", "2:lXCorr:dYSector:lz","sectorCut&&fitCut",0,10000000);
    AliTMinuitToolkit::FitterFCN(npar,0,chi2,(Double_t*)(fitterN->GetParameters()->GetMatrixArray()), flag);
    fitterN->FillFitter(AliTPCDistortionFit::fgkDistortionTree,"dYRun:1", "3:lXCorr:dYSector:lz","sectorCut&&fitCut",0,10000000);
    AliTMinuitToolkit::FitterFCN(npar,0,chi2,(Double_t*)(fitterN->GetParameters()->GetMatrixArray()), flag);
    ((*(fitterN->GetStreamer())).GetFile())->Flush();
    TTree * t = ((*(fitterN->GetStreamer()))<<"fcnDebugPoint").GetTree();
    t->ResetBranchAddresses();
    t->SetMarkerStyle(25); t->SetMarkerSize(0.5);    
    t->SetAlias("xlocal","x.fElements[1]");
    t->SetAlias("yedge","100*x.fElements[2]/xlocal");
    TCanvas * canvasFit = new TCanvas("canvasFit","canvasFit",1400,900);
    canvasFit->Divide(3,2);
    Double_t maxY,minY;
    // Retrieve the stat box
    gStyle->SetTitleOffset(0.92,"X");
    gStyle->SetTitleOffset(1.2,"Y");
    gStyle->SetTitleOffset(1,"Z");
    for (Int_t ipad=0; ipad<6;ipad++){
      Int_t icol=ipad%3;
      Int_t irow=ipad/3;
      TVirtualPad *vpad=canvasFit->cd(ipad+1);
      vpad->SetPad(icol/3., irow/2., (icol+1)/3., (irow+1)/2.);
      if (irow==1) {
	vpad->SetBottomMargin(0);
	vpad->SetTopMargin(0);
      }else{
	vpad->SetTopMargin(0);
	vpad->SetBottomMargin(0.15);
      }
      if (icol<2) {
	vpad->SetRightMargin(0);
      }else {
	vpad->SetRightMargin(0.2);
      }
      if (icol>0) {
	vpad->SetLeftMargin(0); 
      }else{
	vpad->SetLeftMargin(0.2); 
      }

      vpad->Draw();
      canvasFit->cd();
    }
    TVectorD rmsEstimator(*(fitterN->GetRMSEstimator()));
    {       
      canvasFit->cd(4);
      t->Draw("val:xlocal:yedge","x.fElements[0]==2&&abs(abs(x.fElements[3]/xlocal)-0.1)<0.1","colz");
      maxY=TMath::Nint(TMath::MaxElement(t->GetSelectedRows(), t->GetV1())+0.5)+0.5;
      minY=TMath::Nint(TMath::MinElement(t->GetSelectedRows(), t->GetV1())-0.5)-0.5;
      t->GetHistogram()->GetYaxis()->SetLimits(minY,maxY);
      t->GetHistogram()->GetYaxis()->SetTitle("#Delta_{R} (cm)"); 
      latex.DrawLatexNDC(0.25,0.85,"Data");
      canvasFit->cd(5);
      t->Draw("fun:xlocal:yedge","x.fElements[0]==2&&abs(abs(x.fElements[3]/xlocal)-0.1)<0.1","colz");
      t->GetHistogram()->GetYaxis()->SetLimits(minY,maxY);
      latex.DrawLatexNDC(0.15,0.85,"Line charge fit");
      canvasFit->cd(6);
      t->Draw("val-fun:xlocal:yedge","x.fElements[0]==2&&abs(abs(x.fElements[3]/xlocal)-0.1)<0.1","colz");
      t->GetHistogram()->GetYaxis()->SetLimits(minY,maxY);
      t->GetHistogram()->GetZaxis()->SetTitle("100 #Delta_{edge}/x_{L} (unit)");
      latex.DrawLatexNDC(0.15,0.85,"Data-Fit");
      canvasFit->cd(1);
      t->Draw("val:xlocal:yedge","x.fElements[0]==3&&abs(abs(x.fElements[3]/xlocal)-0.1)<0.1","colz");
      maxY=TMath::Nint(TMath::MaxElement(t->GetSelectedRows(), t->GetV1())+0.5)+0.5;
      minY=TMath::Nint(TMath::MinElement(t->GetSelectedRows(), t->GetV1())-0.5)-0.5;
      t->GetHistogram()->GetYaxis()->SetLimits(minY,maxY);
      t->GetHistogram()->GetYaxis()->SetTitle("#Delta_{R#phi} (cm)");
      t->GetHistogram()->GetXaxis()->SetTitle("local X (cm)");
      canvasFit->cd(2);
      t->Draw("fun:xlocal:yedge","x.fElements[0]==3&&abs(abs(x.fElements[3]/xlocal)-0.1)<0.1","colz");
      t->GetHistogram()->GetYaxis()->SetLimits(minY,maxY);
      t->GetHistogram()->GetXaxis()->SetTitle("local X (cm)");
      canvasFit->cd(3);
      t->Draw("val-fun:xlocal:yedge","x.fElements[0]==3&&abs(abs(x.fElements[3]/xlocal)-0.1)<0.1","colz");
      t->GetHistogram()->GetYaxis()->SetLimits(minY,maxY);
      t->GetHistogram()->GetZaxis()->SetTitle("100 #Delta_{edge}/x_{L} (unit)");
      t->GetHistogram()->GetXaxis()->SetTitle("local X (cm)");
      //t->GetHistogram()->GetListOfFunctions()->Add(fitterN->GetFormula());
      //t->GetHistogram()->Draw("colz");
      latex.DrawLatexNDC(0.05,0.9,fName.Data());
      latex.DrawLatexNDC(0.05,0.85,TString::Format("Sector=%d",jsec).Data());
      latex.DrawLatexNDC(0.05,0.8,TString::Format("Q=%2.2f#pm%2.2f",(*(fitterN->GetParameters()))[0],rmsEstimator[0]).Data());
      latex.DrawLatexNDC(0.05,0.75,TString::Format("R=%2.2f#pm%2.2f (cm)",(*(fitterN->GetParameters()))[1],rmsEstimator[1]).Data());
      latex.DrawLatexNDC(0.05,0.70,TString::Format("R#phi=%2.2f#pm%2.2f (cm)",(*(fitterN->GetParameters()))[2],rmsEstimator[2]).Data());
      //      latex.DrawLatexNDC(0.05,0.65,TString::Format("scale=%2.2f#pm%2.2f (cm)",(*(fitterN->GetParameters()))[3],rmsEstimator[3]).Data());
      latex.DrawLatexNDC(0.05,0.65,TString::Format("#omega#tau=%2.2f#pm%2.2f",(*(fitterN->GetParameters()))[4],rmsEstimator[4]).Data());
    }
    canvasFit->SaveAs(TString::Format("canvasDistortionFitLineN_Sec%d.png",fsectors[isec]).Data());
  }


}


void DrawSummary(const char * period, const char *pass){
  //
  // const char * period="LHC15o"; 
  // 
  TFile *ff = TFile::Open("makeFitExample1.root");
  TTree * tree = (TTree*)ff->Get("fit1D");
  AliExternalInfo info;
  TTree * treeLogbook0= info.GetTree("Logbook",period,"");
  TTree * treeLogbook= info.GetChain("Logbook","LHC***","");
  treeLogbook->BuildIndex("run");
  tree->AddFriend(treeLogbook,"Logbook");
  TLine lineGap, lineROC, lineMedian;
  lineGap.SetLineStyle(2); lineGap.SetLineWidth(2);  lineGap.SetLineColor(3);
  lineROC.SetLineStyle(2); lineROC.SetLineWidth(2);   lineROC.SetLineColor(4);    
  lineMedian.SetLineStyle(1); lineMedian.SetLineWidth(3);   lineMedian.SetLineColor(2);    
  TLatex latex, latexMed;
  latex.SetTextSize(0.15);
  latexMed.SetTextSize(0.15);  latexMed.SetTextColor(2);
  Int_t fsectors[16]={2,   4,  6,   7,    9,   10, 19,  20,    29, 30,   31,   35};
  TCanvas *canvasInfo  = new TCanvas("canvasInfo","canvasInfo",1400,1000);
  canvasInfo->SetLeftMargin(0.12);
  canvasInfo->SetBottomMargin(0.2);
  canvasInfo->Divide(2,6,0,0);  
  gStyle->SetOptTitle(0);
  //
  gStyle->SetTitleSize(0.12,"Y");
  gStyle->SetLabelSize(0.12,"Y");
  gStyle->SetTitleYOffset(0.4);
  Int_t npoints=0;
  for (Int_t iSec=0;iSec<12; iSec++){
    canvasInfo->cd(iSec+1);
    TCut cut = TString::Format("sec==%d",fsectors[iSec]).Data();    
    TGraph *gr = TStatToolkit::MakeGraphSparse(tree,"par.fElements[2]:run:rms.fElements[2]",cut,25,1,0);
    gr->SetMaximum(2.5);
    gr->SetMinimum(-2.5);
    gr->GetXaxis()->SetTitle("run");
    gr->GetYaxis()->SetTitle("hot spot r#phi (cm)");
    Double_t median=TMath::Mean(gr->GetN(),gr->GetY());    
    gr->Draw("ap");
    npoints=gr->GetN();
    lineMedian.DrawLine(0.0,median,npoints,median);
    lineGap.DrawLine(0.0,0.2,npoints,0.2);
    lineGap.DrawLine(0.0,-0.2,npoints,-0.2);
    lineROC.DrawLine(0.0,1.3,npoints,1.3);
    lineROC.DrawLine(0.0,-1.3,npoints,-1.3);
    latex.DrawLatexNDC(0.2,0.85 , TString::Format("sector %d",fsectors[iSec]).Data());
    latexMed.DrawLatexNDC(0.2,0.7 , TString::Format("#Delta_{R#phi}=%2.2f (cm)",median).Data());
  }
  canvasInfo->SaveAs("hotSpotsRPhiPos.png");
  
  for (Int_t iSec=0;iSec<12; iSec++){
    canvasInfo->cd(iSec+1);
    TCut cut = TString::Format("sec==%d",fsectors[iSec]).Data();    
    TGraph *gr = TStatToolkit::MakeGraphSparse(tree,"par.fElements[1]:run:rms.fElements[2]",cut,25,1,0);
    gr->GetXaxis()->SetTitle("run");
    gr->GetYaxis()->SetTitle("hot spot r(cm)");
    Double_t median=TMath::Mean(gr->GetN(),gr->GetY());
    Double_t medianI=TMath::Nint(median);
    gr->SetMinimum(medianI-5);
    gr->SetMaximum(medianI+5);
    gr->Draw("ap");
    lineMedian.DrawLine(0.0,median,npoints,median);
    latex.DrawLatexNDC(0.2,0.85 , TString::Format("Sector %d",fsectors[iSec]).Data());
    latexMed.DrawLatexNDC(0.2,0.7 , TString::Format("R=%2.2f (cm)",median).Data());
  }
  canvasInfo->SaveAs("hotSpotsRPos.png");

  gStyle->SetOptFit(1);
  for (Int_t iSec=0;iSec<12; iSec++){
    canvasInfo->cd(iSec+1);
    TCut cutP = TString::Format("sec==%d&&L3_magnetCurrent>0",fsectors[iSec]).Data();    
    TCut cutN = TString::Format("sec==%d&&L3_magnetCurrent<0",fsectors[iSec]).Data();    
    TGraph *grP = TStatToolkit::MakeGraphErrors(tree,"par.fElements[0]:lumi:rms.fElements[0]",cutP,25,2,0.6);
    TGraph *grN = TStatToolkit::MakeGraphErrors(tree,"par.fElements[0]:lumi:rms.fElements[0]",cutN,21,4,0.6);
    grP->GetXaxis()->SetTitle("lumi");
    grP->GetYaxis()->SetTitle("charge density (V)");
    grP->Draw("ap");
    grN->Draw("p");
    //grP->Fit("pol1");
    //grN->Fit("pol1");
    grP->SetMaximum(TMath::Max( grP->GetMaximum(), grN->GetMaximum()));
    grP->Draw("ap");
    grN->Draw("p");    
    latex.DrawLatexNDC(0.3,0.7 , TString::Format("sector %d",fsectors[iSec]).Data());
  }
  canvasInfo->SaveAs("hotSpotsLumi.png");


}





