/*
  
gSystem->AddIncludePath("-I$ALICE_ROOT/TPC");
gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");

// gROOT->LoadMacro("$ALICE_ROOT/TPC/CalibMacros/CalibAlignKalman.C+");
// AliTPCTransformation::BuildBasicFormulas();

// AliXRDPROOFtoolkit tool;
// chainPoints = tool.MakeChainRandom("align.txt","trackPoints",0,50000);
// chainPoints->Lookup();

// CalibAlignKalman(40000);
// kalmanFit0->DumpCorelation(0.8);
// TFile f("kalmanfitTPC.root");


*/
#include <fstream>

#include "TSystem.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TMath.h"
#include "TBits.h"
#include "TFormula.h"
#include "TF1.h"
#include "TLinearFitter.h"
#include "TFile.h"
#include "TChain.h"
#include "TCut.h"
#include "TEntryList.h"
#include "TH1F.h"


#include "TTreeStream.h"
#include "AliTrackPointArray.h"
#include "AliLog.h"
#include "AliTPCTransformation.h"
#include "AliTPCkalmanFit.h"
#include "AliXRDPROOFtoolkit.h"


//
TChain *chainPoints=0; 
TEntryList *elist=0;
AliTPCkalmanFit * kalmanFit0=0;
AliTPCkalmanFit * kalmanFitIdeal=0;

AliTPCkalmanFit *  CalibAlignKalmanFit(Int_t maxTracks);
void FilterTracks();
AliTPCkalmanFit * SetupFit();
void              AddFitFieldCage(AliTPCkalmanFit *kalmanFit);
void              AddPhiScaling(AliTPCkalmanFit *kalmanFit);
void              AddZShift(AliTPCkalmanFit *kalmanFit);
void              AddZRotation(AliTPCkalmanFit *kalmanFit);
void              AddLocalXYMisalignment(AliTPCkalmanFit *kalmanFit);
AliTrackPointArray *FilterPoints(AliTrackPointArray &points, Int_t dir, TTreeSRedirector *pcstream);
AliTPCkalmanFit *   FitPointsLinear(Int_t maxTracks);

void CalibAlignKalman(Int_t npoints, Int_t maxFiles, Int_t startFile){
  //
  //
  //
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
  
  AliTPCTransformation::BuildBasicFormulas();
  AliXRDPROOFtoolkit tool;
  chainPoints = tool.MakeChainRandom("align.txt","trackPoints",0, maxFiles, startFile);
  chainPoints->Lookup();
  CalibAlignKalmanFit(npoints);
}



//
//
// track point cuts
//
const Float_t krmsYcut    = 0.3;
const Float_t krmsZcut    = 0.3;
const Float_t kSigmaCut   = 10.;
const Int_t   knclCut     =  80;
const Double_t kArmCut    = 50.;
//
// Track cuts
//
TCut cutsP3="abs(p0Out.fP[3]-p0In.fP[3])<0.0002&&abs(p1Out.fP[3]-p1In.fP[3])<0.0002";
TCut cutsP4="abs(p0Out.fP[4]-p0In.fP[4])<0.005&&abs(p1Out.fP[4]-p1In.fP[4])<0.005";
TCut cutP3 ="abs(p1In.fP[3]+p0In.fP[3]-0.0015)<0.0025";
TCut cutP4 ="abs(p1In.fP[4])<0.03&&abs(p0In.fP[4])<0.03";
TCut cutSide="p1Out.fP[1]*p0Out.fP[1]>0&&p1In.fP[1]*p0In.fP[1]>0";
TCut cutAll=cutSide+cutsP4+cutsP3+cutP3+cutP4+"abs(mag)<0.01&&ncont>0&&p.fNPoints>120";




AliTPCkalmanFit *  CalibAlignKalmanFit(Int_t maxTracks){
  //
  //
  AliTPCTransformation::BuildBasicFormulas();
  FilterTracks();
  kalmanFit0     = SetupFit();  
  kalmanFitIdeal = SetupFit();
  kalmanFit0->Init();
  kalmanFitIdeal->Init();
  return FitPointsLinear(maxTracks);
}



AliTPCkalmanFit * SetupFit(){
  //
  AliTPCkalmanFit *kalmanFit =  new AliTPCkalmanFit;
  AddFitFieldCage(kalmanFit); 
  AddPhiScaling(kalmanFit);
  AddZShift(kalmanFit); 
  AddZRotation(kalmanFit); 
  AddLocalXYMisalignment(kalmanFit);  
  return kalmanFit;
}


void FilterTracks(){
  //
  //
  //
  chainPoints->Draw(">>listEL",cutAll,"entryList");
  elist = (TEntryList*)gDirectory->Get("listEL");
  chainPoints->SetEntryList(elist);
  elist->SetDirectory(0);
}



AliTPCkalmanFit * FitPointsLinear(Int_t maxTracks){
  //
  //
  //
  // create debug streeamers
  TTreeSRedirector *pcstream      = new TTreeSRedirector("kalmanfitTPC.root");  
  TTreeSRedirector *pcstreamIdeal = new TTreeSRedirector("kalmanfitTPCOrig.root");  
  //
  //
  AliTrackPointArray *points=0;
  Float_t mag=0;
  Int_t   time=0;
  chainPoints->SetBranchAddress("p.",&points);
  chainPoints->SetBranchAddress("mag",&mag);
  chainPoints->SetBranchAddress("time",&time);
  Int_t accepted=0;
  //
  for (Int_t itrack=0;itrack<elist->GetN(); itrack++){
    Int_t entry=chainPoints->GetEntryNumber(itrack);
    chainPoints->GetEntry(entry);
    if (TMath::Abs(mag)>0.01) {printf("mag- Cut not accpeted\n"); continue;}
    if (points->GetNPoints()<knclCut) {printf("ncl - Cut not accpeted\n"); continue;}
    if (accepted>maxTracks) break;
    for (Int_t idir=-1; idir<=1; idir++){
      AliTrackPointArray *pointsF = FilterPoints(*points,idir, pcstream);
      if (!pointsF) continue;
      if (idir==0)  accepted++;
      //
      if (accepted%50==0) {
	kalmanFit0->FitTrackLinear(*pointsF, 10, pcstream);
      }else{
	kalmanFit0->FitTrackLinear(*pointsF, 1, 0);
      }    
      if (idir==0) kalmanFit0->DumpTrackLinear(*pointsF,pcstream);
      if (idir==0) kalmanFitIdeal->DumpTrackLinear(*pointsF,pcstreamIdeal);
      if (accepted%25==0) printf("%d\n", accepted);
      delete pointsF;
    }
  }
  pcstream->GetFile()->cd();
  kalmanFit0->Write("kalmanFit");
  pcstreamIdeal->GetFile()->cd();
  kalmanFitIdeal->Write("kalmanFitIdeal");
  delete pcstream;  
  delete pcstreamIdeal;  
  return kalmanFit0;   
}




AliTrackPointArray *FilterPoints(AliTrackPointArray &points, Int_t dir, TTreeSRedirector */*pcstream*/){
  //
  //
  //
  TLinearFitter lfitY(2,"pol1");
  TLinearFitter lfitZ(2,"pol1");
  TVectorD vecZ(2);
  TVectorD vecY(2);
  //
  lfitY.StoreData(kTRUE);
  lfitZ.StoreData(kTRUE);
  Int_t npoints = points.GetNPoints();
  if (npoints<2) return 0;
  Double_t currentAlpha = TMath::ATan2(points.GetY()[npoints-1]-points.GetY()[0], points.GetX()[npoints-1]-points.GetX()[0]);  
  Double_t ca = TMath::Cos(currentAlpha);
  Double_t sa = TMath::Sin(currentAlpha);
  //
  // 1.b Fit the track in the rotated frame - MakeSeed 
  //
  Double_t maxX =-10000, minX=10000;
  for (Int_t ipoint=0; ipoint<npoints-1; ipoint++){
    Double_t rx =   ca*points.GetX()[ipoint]+sa*points.GetY()[ipoint];
    Double_t ry =  -sa*points.GetX()[ipoint]+ca*points.GetY()[ipoint];
    Double_t rz =  points.GetZ()[ipoint];
    if (dir== 1 && rx<0) continue;
    if (dir==-1 && rx>0) continue;
    if (maxX<rx) maxX=rx;
    if (minX>rx) minX=rx;
    lfitY.AddPoint(&rx,ry,1);
    lfitZ.AddPoint(&rx,rz,1);
  }
  if (TMath::Abs(maxX-minX)<kArmCut) return 0;
  if (lfitY.GetNpoints()<knclCut) return 0;
  //
  lfitY.Eval();
  lfitZ.Eval();
  lfitY.GetParameters(vecY);
  lfitZ.GetParameters(vecZ);
  //
  Double_t chi2Y = lfitY.GetChisquare()/lfitY.GetNpoints();
  Double_t chi2Z = lfitZ.GetChisquare()/lfitZ.GetNpoints();
  if (TMath::Sqrt(chi2Y)>krmsYcut) return 0;
  if (TMath::Sqrt(chi2Z)>krmsZcut) return 0;
  //
  //
  Int_t accepted=0;
  AliTrackPoint point;
  AliTrackPointArray *pointsF=0;
  for (Int_t iter=0; iter<2;iter++){
    for (Int_t ipoint=0; ipoint<npoints-1; ipoint++){
      //
      if (!points.GetPoint(point,ipoint)) continue;
      Double_t rx =   ca*points.GetX()[ipoint]+sa*points.GetY()[ipoint];
      Double_t ry =  -sa*points.GetX()[ipoint]+ca*points.GetY()[ipoint];
      Double_t rz =  points.GetZ()[ipoint];
      if (dir== 1 && rx<0) continue;
      if (dir==-1 && rx>0) continue;
      Double_t erry = TMath::Sqrt(chi2Y);
      Double_t errz = TMath::Sqrt(chi2Z);
      Double_t fy = vecY[0]+vecY[1]*rx;
      Double_t fz = vecZ[0]+vecZ[1]*rx;
      if (TMath::Abs(fy-ry)>erry*kSigmaCut) continue;
      if (TMath::Abs(fz-rz)>errz*kSigmaCut) continue;
      if (pointsF) pointsF->AddPoint(accepted,&point);
      accepted++;
    }
    if (accepted<knclCut) break;
    if (iter==0) pointsF = new AliTrackPointArray(accepted);
    accepted=0;
  }
  return pointsF;
}




void  AddFitFieldCage(AliTPCkalmanFit *kalmanFit){
  //
  // Add radial scaling due field cage
  //
  TVectorD fpar(10);
  AliTPCTransformation * transformation=0;
  char tname[100];
  //
  // linear R scaling and shift
  //  
  for (Int_t iside=0; iside<=1; iside++)
    for (Int_t ipolR=0; ipolR<2; ipolR++){
      for (Int_t ipolZ=0; ipolZ<3; ipolZ++){
	fpar[0]=ipolR;
	fpar[1]=ipolZ;
	if (ipolR+ipolZ==0) continue;
	sprintf(tname,"tTPCscalingRPolR%dDr%dSide%d",ipolR,ipolZ,iside);
	transformation = new AliTPCTransformation(tname,AliTPCTransformation::BitsSide(iside),"TPCscalingRPol",0,0,  1);
	transformation->SetParams(0,0.2,0,&fpar);
	kalmanFit->AddCalibration(transformation);      
	//      
      }
    }
  //
  //
  //Inner field cage  
  for (Int_t iside=0; iside<=1; iside++)
    for (Int_t ipol=0; ipol<3; ipol++){
      fpar[0]=ipol; 
      sprintf(tname,"tTPCscalingRIFC%dSide%d",ipol,iside);
      transformation = new AliTPCTransformation(tname,AliTPCTransformation::BitsSide(iside),"TPCscalingRIFC",0,0,   1);
      transformation->SetParams(0,0.2,0,&fpar);
      kalmanFit->AddCalibration(transformation);
    }
  //
  //
  //Outer field cage  
  for (Int_t iside=0; iside<=1; iside++)
    for (Int_t ipol=0; ipol<3; ipol++){
      fpar[0]=ipol;
      //Outer field cage
      sprintf(tname,"tTPCscalingROFC%dSide%d",ipol,iside);
      transformation = new AliTPCTransformation(tname,AliTPCTransformation::BitsSide(iside),"TPCscalingROFC",0,0,  1);
      transformation->SetParams(0,0.2,0,&fpar);
      kalmanFit->AddCalibration(transformation);
    }
}


void AddPhiScaling(AliTPCkalmanFit *kalmanFit){
  //
  // Add local phi scaling
  //
  TVectorD fpar(10);
  AliTPCTransformation * transformation=0;
  fpar[0]=1;
  transformation = new AliTPCTransformation("tscalingLocalPhi", AliTPCTransformation::BitsAll(), 0,"TPCscalingPhiLocal",0,  1);
  transformation->SetParams(0,0.2,0,&fpar);
  kalmanFit->AddCalibration(transformation);  
  //
}

void  AddZShift(AliTPCkalmanFit *kalmanFit){
  TVectorD fpar(10);
  AliTPCTransformation * transformation=0;
  TBits maskInnerA(72);
  TBits maskInnerC(72);
  TBits maskOuter(72);
  for (Int_t i=0;i<72;i++){
    if (i<36){
      if (i%36<18)  maskInnerA[i]=kTRUE;
      if (i%36>=18) maskInnerC[i]=kTRUE;
    }
    if (i>=36){
      maskOuter[i]=kTRUE;
    }
  }
  //
  transformation = new AliTPCTransformation("tTPCDeltaZIROCA", new TBits(maskInnerA), 0,0, "TPCDeltaZ",  0);
  transformation->SetParams(0,0.2,0,&fpar);
  kalmanFit->AddCalibration(transformation);
  transformation = new AliTPCTransformation("tTPCDeltaZIROCC", new TBits(maskInnerC), 0,0, "TPCDeltaZ",  0);
  transformation->SetParams(0,0.2,0,&fpar);
  kalmanFit->AddCalibration(transformation);
  //
  transformation = new AliTPCTransformation("tTPCDeltaZOROCML", new TBits(maskOuter), 0,0, "TPCDeltaZMediumLong",  0);
  transformation->SetParams(0,0.2,0,&fpar);
  kalmanFit->AddCalibration(transformation);
  //
}

void AddZRotation(AliTPCkalmanFit *kalmanFit){
  //
  //
  //
  TVectorD fpar(10);
  //  AliTPCTransformation::BitsSide(iside)
  fpar[0]=0; fpar[1]=0;
  AliTPCTransformation * transformation=0;
  transformation = new AliTPCTransformation("tTPCTiltingZAside00",AliTPCTransformation::BitsSide(0) , 0,0, "TPCTiltingZ",  0);
  transformation->SetParams(0,0.4,0,&fpar);
  kalmanFit->AddCalibration(transformation);
  transformation = new AliTPCTransformation("tTPCTiltingZCside00",AliTPCTransformation::BitsSide(1) , 0,0, "TPCTiltingZ",  0);
  transformation->SetParams(0,0.4,0,&fpar);
  kalmanFit->AddCalibration(transformation);
  //
  fpar[0]=1; fpar[1]=0;
  transformation = new AliTPCTransformation("tTPCTiltingZAside10",AliTPCTransformation::BitsSide(0) , 0,0, "TPCTiltingZ",  0);
  transformation->SetParams(0,0.4,0,&fpar);
  kalmanFit->AddCalibration(transformation);
  transformation = new AliTPCTransformation("tTPCTiltingZCside10",AliTPCTransformation::BitsSide(1) , 0,0, "TPCTiltingZ",  0);
  transformation->SetParams(0,0.4,0,&fpar);
  kalmanFit->AddCalibration(transformation);
  //
  //
  fpar[0]=0; fpar[1]=1;
  transformation = new AliTPCTransformation("tTPCTiltingZAside01",AliTPCTransformation::BitsSide(0) , 0,0, "TPCTiltingZ",  0);
  transformation->SetParams(0,0.4,0,&fpar);
  kalmanFit->AddCalibration(transformation);
  transformation = new AliTPCTransformation("tTPCTiltingZCside01",AliTPCTransformation::BitsSide(1) , 0,0, "TPCTiltingZ",  0);
  transformation->SetParams(0,0.4,0,&fpar);
  kalmanFit->AddCalibration(transformation);  
}



void  AddLocalXYMisalignment(AliTPCkalmanFit *kalmanFit){
  //
  //
  //
  TVectorD fpar(10);
  AliTPCTransformation * transformation=0;
  TBits maskInnerA(72);
  TBits maskInnerC(72);
  for (Int_t i=0;i<72;i++){
    if (i<36){
      if (i%36<18)  maskInnerA[i]=kTRUE;
      if (i%36>=18) maskInnerC[i]=kTRUE;
    }
  }
  //
  //
  transformation = new AliTPCTransformation("tTPCDeltaLxIROCA", new TBits(maskInnerA), "TPClocaldLxdGX","TPClocaldLxdGY",0,  0);
  transformation->SetParams(0,0.2,0,&fpar);
  kalmanFit->AddCalibration(transformation);
  transformation = new AliTPCTransformation("tTPCDeltaLxIROCC", new TBits(maskInnerC), "TPClocaldLxdGX","TPClocaldLxdGY",0,  0);
  transformation->SetParams(0,0.2,0,&fpar);
  kalmanFit->AddCalibration(transformation);
  //
  transformation = new AliTPCTransformation("tTPCDeltaLyIROCA", new TBits(maskInnerA), "TPClocaldLydGX","TPClocaldLydGY",0,  0);
  transformation->SetParams(0,0.2,0,&fpar);
  kalmanFit->AddCalibration(transformation);
  transformation = new AliTPCTransformation("tTPCDeltaLyIROCC", new TBits(maskInnerC), "TPClocaldLydGX","TPClocaldLydGY",0,  0);
  transformation->SetParams(0,0.2,0,&fpar);
  kalmanFit->AddCalibration(transformation);
}

void MergeKalman(const char * list = "kalmanFit.list"){
  //
  //
  //
  ifstream in;
  in.open(list);
  TString currentFile;
  kalmanFit0= 0;
  Int_t counter=0;
  while(in.good()) {
    in >> currentFile;
    printf("%d\t%d\t%s\n", counter,currentFile.Length(),currentFile.Data());
    if (currentFile.Length()==0) continue;
    TFile * ffit = TFile::Open(currentFile.Data());
    AliTPCkalmanFit * fit = ( AliTPCkalmanFit *)ffit->Get("kalmanFit");
    if (!fit) continue;
    if (!kalmanFit0) {kalmanFit0= fit; continue;};
    kalmanFit0->Add(fit);
    delete fit;
    delete ffit;
    counter++;
  }
}


/*
  myvar=0; 
  ntracks=10000000
  bDir=`pwd`
  while [ $myvar -ne 190 ] ; do mkdir kalmanDir$myvar; cd kalmanDir$myvar; cp $bDir/align.txt .;  bsub -q proof command aliroot  -q -b  "CalibAlignKalman.C($ntracks,5,$myvar)" ; myvar=$(( $myvar + 5 )) ; echo $myvar ; cd $bDir; echo $bDir; done

*/












