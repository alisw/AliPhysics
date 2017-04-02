/*
Class for fitting of distortion maps using phsyical models

*/

/*
  gSystem->AddIncludePath("-I$AliRoot_SRC/TPC/TPCcalib/");
  .L $AliRoot_SRC/TPC/TPCcalib/AliTPCDistortionFit.cxx+
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

using namespace std;

#include "AliTPCDistortionFit.h"
using std::cout;

ClassImp(AliTPCDistortionFit)

map<string,int>    AliTPCDistortionFit::fgkMapNameHash;              // map  name->index
map<int,string>  AliTPCDistortionFit::fgkMapHashName;                // map  index->name
map<int,const AliTPCChebCorr *> AliTPCDistortionFit::fgkMapCheb;
TTree *  AliTPCDistortionFit::fgkDistortionTree;     



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


void AliTPCDistortionFit::LoadDistortionTree(const char *chinput){
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
    return;
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
    return;
  }
  Int_t run=cheb->GetRun();
  RegisterMap(cheb->GetName(),cheb);
  treeLoad->SetMarkerStyle(21);
  treeLoad->SetMarkerSize(0.6);
  if (fgkDistortionTree==NULL) {
    TFile * finput2=TFile::Open(chinput);
    fgkDistortionTree=(TTree*)finput2->Get("voxRes");
    SetMetadata(fgkDistortionTree,"",run);
  }  
  TString refMap=TString::Format("R%d.refMap",run);
  if (GetCheb(refMap.Data())==NULL) LoadDistortionMaps(run,NULL);
  fgkDistortionTree->AddFriend(treeLoad,cheb->GetName());
  SetMetadata(fgkDistortionTree,cheb->GetName(),run); 
}


void  AliTPCDistortionFit::SetMetadata(TTree * tree, TString friendName, Int_t run){
  //
  //
  //
  TString  refMap=TString::Format("R%d.refMap",run);
  Int_t hashRef=refMap.Hash();
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
    ftree->SetAlias("dX",TString::Format("%s.D[1]-AliTPCDistortionFit::Eval(%d,sector,x,y2xAV,z2xAV,0+0)",friendName.Data(),hashRef).Data());
    ftree->SetAlias("dY",TString::Format("%s.D[1]-AliTPCDistortionFit::Eval(%d,sector,x,y2xAV,z2xAV,1+0)",friendName.Data(),hashRef).Data());
    ftree->SetAlias("dZ",TString::Format("%s.D[1]-AliTPCDistortionFit::Eval(%d,sector,x,y2xAV,z2xAV,2+0)",friendName.Data(),hashRef).Data());
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

void AliTPCDistortionFit::PrintMap(){
  // Loop through the map (Would be much easier in c++11)
  // looping over map with const_iterator
  typedef std::map<string,int>::const_iterator it_type;
  for(it_type iterator = fgkMapNameHash.begin(); iterator != fgkMapNameHash.end(); ++iterator) {
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


