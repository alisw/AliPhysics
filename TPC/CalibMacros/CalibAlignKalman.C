/// \file CalibAlignKalman.C
///
/// ~~~{.cpp}
/// gSystem->AddIncludePath("-I$ALICE_ROOT/TPC");
/// gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
/// gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
///
/// gROOT->LoadMacro("$ALICE_ROOT/TPC/CalibMacros/CalibAlignKalman.C+");
///
/// AliTPCTransformation::BuildBasicFormulas();
///
/// AliXRDPROOFtoolkit tool;
/// chainPoints = tool.MakeChainRandom("align.txt","trackPoints",0,50000);
/// chainPoints->Lookup();
///
/// chainMS = tool.MakeChainRandom("kalmanFit.list","kf",0,50000);
/// chainMS->Lookup();
///
/// chainFP = tool.MakeChainRandom("kalmanFit.list","filter",0,50000);
/// chainFP->Lookup();
///
/// CalibAlignKalmanFit(40000,1);
/// kalmanFit0->DumpCorelation(0.8);
/// TFile f("kalmanfitTPC.root");
/// ~~~


#ifdef __CINT__
#else
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
#include "THnSparse.h"


#include "AliSysInfo.h"
#include "AliExternalTrackParam.h"
#include "TTreeStream.h"
#include "AliTrackPointArray.h"
#include "AliLog.h"
#include "AliTPCTransformation.h"
#include "AliTPCkalmanFit.h"
#include "AliMathBase.h"
#include "AliXRDPROOFtoolkit.h"
#endif


//
//
// track point cuts
//
const Float_t krmsYcut      = 0.2;   // cluster RMS cut in Y - residuals form local  fit
const Float_t krmsZcut      = 0.2;   // cluster RMS cut in Z - residuals from local  fit 
const Float_t krmsYcutGlobal= 0.2;   // cluster RMS cut in Y - residuals form global fit
const Float_t krmsZcutGlobal= 2.0;   // cluster RMS cut in Z - residuals from global fit 
const Float_t kSigmaCut     = 5.;    // clusters to be removed 
const Int_t   knclCut       =  80;   // minimal number of clusters
const Double_t kArmCut      = 50.;   // minimal level arm
const Double_t kNsigma      = 5.;    // cut on track level
//
// mult. scattering cuts    
/*
TCut cutClY("rms09.fElements[0]<0.15");
TCut cuttY("rms09.fElements[1]/rms09.fElements[0]<0.9");
TCut cutkY("rms09.fElements[2]<0.015");
TCut cutClZ("rms09.fElements[3]<0.2");
TCut cuttZ("rms09.fElements[4]/rms09.fElements[3]<0.9");
TCut cutkZ("rms09.fElements[5]<0.015");
TCut cutMS=cutClY+cuttY+cutkY+cutClZ+cuttZ+cutkZ
*/
TMatrixD cutMatrix(4*7,2);
const Double_t rmsCut09[6]={0.15,0.9,0.015, 0.2, 0.9, 0.015};



//
//
Int_t          toSkip       = 2;     // 
Int_t          toSkipOffset = 0;
Int_t          toSkipTrack       = 2;
Int_t          toSkipTrackOffset = 0;
Int_t          isFilterTest = 0;
//
// Track cuts
//
TCut * cSide[4]={0,0,0,0};
TCut *cP3[4]={0,0,0,0};
TCut *cSP3[4]={0,0,0,0};
TCut *cP4[4]={0,0,0,0};
TCut *cM4[4]={0,0,0,0};
TCut *cA[4]={0,0,0,0};

TCut cutAll="";


//
//
//
TChain *chainPoints=0; 
TEntryList *elist=0;
AliTPCkalmanFit * kalmanFitNew=0;
AliTPCkalmanFit * kalmanFitOrig=0;
AliTPCkalmanFit * kalmanFitApply=0;

AliTPCkalmanFit *  CalibAlignKalmanFit(Int_t maxTracks, Int_t trackDump);
void FilterTracks();
AliTPCkalmanFit * SetupFit();
//
void              AddFitFieldCage(AliTPCkalmanFit *kalmanFit);
void              AddPhiScaling(AliTPCkalmanFit *kalmanFit);
void              AddZShift(AliTPCkalmanFit *kalmanFit, Int_t ncos, Int_t nsin );
void              AddZTilting(AliTPCkalmanFit *kalmanFit, Int_t ncos, Int_t nsin);
void              AddLocalXYMisalignment(AliTPCkalmanFit *kalmanFit);
void              AddLocalXYMisalignmentSector(AliTPCkalmanFit *kalmanFit);
void              AddAlignOROCIROCFourier(AliTPCkalmanFit *kalmanFit, Int_t ncos, Int_t nsin);
void              AddAlignSectorFourier(AliTPCkalmanFit *kalmanFit, Int_t ncos, Int_t nsin);
void              AddDrift(AliTPCkalmanFit *kalmanFit);


AliTrackPointArray *FilterPoints(AliTrackPointArray &points, Int_t dir, TTreeSRedirector *pcstream);

TVectorD * EstimateScatteringKalmanLinear(AliTrackPointArray &points, AliExternalTrackParam &p0, AliExternalTrackParam &p1 , TTreeSRedirector *pcstream);
AliTrackPointArray *SkipPoints(AliTrackPointArray &points, Int_t nskip, Int_t nskipOffset);

AliTPCkalmanFit *   FitPointsLinear(Int_t maxTracks, Int_t trackDump);

void CalibAlignKalman(Int_t npoints, Int_t maxFiles, Int_t startFile, Int_t trackDump, Int_t nSkipTrack, Int_t nSkipTrackOffset, Int_t nSkip, Int_t nSkipOffset, Int_t bfilterTest){
  ///

  AliTPCTransformation::BuildBasicFormulas();
  toSkip=nSkip;
  toSkipOffset= nSkipOffset;
  toSkipTrack = nSkipTrack;
  toSkipTrackOffset = nSkipTrackOffset;
  isFilterTest = bfilterTest;
  //
  // read the transformation to be applied
  TFile ftrafo("kalmanFitApply.root");
  kalmanFitApply = (AliTPCkalmanFit *)ftrafo.Get("kalmanFitNew");
  if (kalmanFitApply) {
    printf("Loaded transforamtion\n");
    kalmanFitApply->DumpCalib("IROCOROC");
    kalmanFitApply->InitTransformation();
  }else{
    printf("Not trnasformation specified\n");
  }
  //
  //
  //
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
  //
  AliXRDPROOFtoolkit tool;
  chainPoints = tool.MakeChainRandom("align.txt","trackPoints",0, maxFiles, startFile);
  CalibAlignKalmanFit(npoints, trackDump);
}




AliTPCkalmanFit *  CalibAlignKalmanFit(Int_t maxTracks, Int_t trackDump){
  /// Fitting procedure

  AliTPCTransformation::BuildBasicFormulas();
  FilterTracks();
  kalmanFitNew     = SetupFit();  
  kalmanFitOrig  = SetupFit();
  kalmanFitNew->Init();
  kalmanFitOrig->Init();
  return FitPointsLinear(maxTracks,trackDump);
}



AliTPCkalmanFit * SetupFit(){
  ///

  AliTPCkalmanFit *kalmanFit =  new AliTPCkalmanFit;
  AddFitFieldCage(kalmanFit); 
  AddPhiScaling(kalmanFit);
  AddDrift(kalmanFit);
  AddZShift(kalmanFit,3,3); 
  AddZTilting(kalmanFit,3,3); 
  //  AddLocalXYMisalignment(kalmanFit);  
  //  AddLocalXYMisalignmentSector(kalmanFit);  
  AddAlignSectorFourier(kalmanFit,4,4); 
  AddAlignOROCIROCFourier(kalmanFit,5,5);
  kalmanFit->Init();
  return kalmanFit;
}


void FilterTracks(){
  ///

  cSide[0] = new TCut("cutAA","p0In.fP[1]>0&&p1In.fP[1]>0");
  cSide[1] = new TCut("cutCC","p0In.fP[1]<0&&p1In.fP[1]<0");
  cSide[2] = new TCut("cutAC","p0In.fP[1]>0&&p1In.fP[1]<0");
  cSide[3] = new TCut("cutCA","p0In.fP[1]<0&&p1In.fP[1]>0");
  //  
  TH1F * phisP3 = new TH1F("hhisP3","hhisP3",100,-0.01,0.01);
  TH1F * phisSP3 = new TH1F("hhisSP3","hhisSP3",100,-0.001,0.001);
  TH1F * phisP4 = new TH1F("hhisP4","hhisP4",100,-0.1,0.1);
  TH1F * phisM4 = new TH1F("hhisM4","hhisM4",100,-0.01,0.01);
  
  //

  TF1 *fg = new TF1("fg","gaus");
  for (Int_t iter=0; iter<2; iter++){
    for (Int_t ic=0;ic<4;ic++){
      if (!cA[ic]){
	cA[ic]=new TCut;
	*cA[ic]= *cSide[ic];
      }
      //
      // cutP3
      //
      chainPoints->Draw("p0In.fP[3]+p1In.fP[3]>>hhisP3",*cA[ic],"goff");
      phisP3->Fit(fg,"QNR","QNR",phisP3->GetBinCenter(phisP3->GetMaximumBin())-0.003,phisP3->GetBinCenter(phisP3->GetMaximumBin())+0.003);
      cP3[ic]=new TCut(Form("abs(p0In.fP[3]+p1In.fP[3]-%f)<%f",fg->GetParameter(1), fg->GetParameter(2)*kNsigma));
      cutMatrix(7*ic+0,0) = fg->GetParameter(1);
      cutMatrix(7*ic+0,1) = fg->GetParameter(2);
      //
      // cutSP3
      //
      chainPoints->Draw("p0Out.fP[3]-p0In.fP[3]>>hhisSP3",*cA[ic],"goff");
      phisSP3->Fit(fg,"QNR","QNR",phisP3->GetBinCenter(phisSP3->GetMaximumBin())-0.0015,phisP3->GetBinCenter(phisSP3->GetMaximumBin())+0.0015);
      cSP3[ic]=new TCut(Form("abs(p0Out.fP[3]-p0In.fP[3]-%f)<%f",fg->GetParameter(1), fg->GetParameter(2)*kNsigma));
      cutMatrix(7*ic+1,0) = fg->GetParameter(1);
      cutMatrix(7*ic+1,1) = fg->GetParameter(2);
      //
      chainPoints->Draw("p1Out.fP[3]-p1In.fP[3]>>hhisSP3",*cA[ic],"goff");
      phisSP3->Fit(fg,"QNR","QNR",phisP3->GetBinCenter(phisSP3->GetMaximumBin())-0.0015,phisP3->GetBinCenter(phisSP3->GetMaximumBin())+0.0015);
      *cSP3[ic]+=Form("abs(p1Out.fP[3]-p1In.fP[3]-%f)<%f",fg->GetParameter(1), fg->GetParameter(2)*kNsigma);  
      cutMatrix(7*ic+2,0) = fg->GetParameter(1);
      cutMatrix(7*ic+2,1) = fg->GetParameter(2);    
      //
      // cutP4
      //
      chainPoints->Draw("p0Out.fP[4]>>hhisP4",*cA[ic],"goff");
      phisP4->Fit(fg,"QNR","QNR",phisP4->GetBinCenter(phisP4->GetMaximumBin())-0.03,phisP4->GetBinCenter(phisP4->GetMaximumBin())+0.03);
      cP4[ic]=new TCut(Form("abs(p0Out.fP[4]-%f)<%f",fg->GetParameter(1), fg->GetParameter(2)*kNsigma));
      cutMatrix(7*ic+3,0) = fg->GetParameter(1);
      cutMatrix(7*ic+3,1) = fg->GetParameter(2);    
      chainPoints->Draw("p1Out.fP[4]>>hhisP4",*cA[ic],"goff");
      phisP4->Fit(fg,"QNR","QNR",phisP4->GetBinCenter(phisP4->GetMaximumBin())-0.03,phisP4->GetBinCenter(phisP4->GetMaximumBin())+0.03);
      *cP4[ic]+=Form("abs(p1Out.fP[4]-%f)<%f",fg->GetParameter(1), fg->GetParameter(2)*kNsigma);
      cutMatrix(7*ic+4,0) = fg->GetParameter(1);
      cutMatrix(7*ic+4,1) = fg->GetParameter(2);    

    //
      // cutM4
      //
      chainPoints->Draw("p0Out.fP[4]-p0In.fP[4]>>hhisM4",*cA[ic],"goff");
      phisM4->Fit(fg,"QNR","QNR",phisM4->GetBinCenter(phisM4->GetMaximumBin())-0.03,phisM4->GetBinCenter(phisM4->GetMaximumBin())+0.03);
      cM4[ic]=new TCut(Form("abs(p0Out.fP[4]-p0In.fP[4]-%f)<%f",fg->GetParameter(1), fg->GetParameter(2)*kNsigma));
      cutMatrix(7*ic+5,0) = fg->GetParameter(1);
      cutMatrix(7*ic+5,1) = fg->GetParameter(2);    

      chainPoints->Draw("p1Out.fP[4]-p1In.fP[4]>>hhisM4",*cA[ic],"goff");
      phisM4->Fit(fg,"QNR","QNR",phisM4->GetBinCenter(phisM4->GetMaximumBin())-0.03,phisM4->GetBinCenter(phisM4->GetMaximumBin())+0.03);
      *cM4[ic]+=Form("abs(p1Out.fP[4]-p1In.fP[4]-%f)<%f",fg->GetParameter(1), fg->GetParameter(2)*kNsigma);
      cutMatrix(7*ic+6,0) = fg->GetParameter(1);
      cutMatrix(7*ic+6,1) = fg->GetParameter(2);    
      //
      //
      //
      cA[ic]=new TCut;
      *cA[ic]= *cSide[ic]+*cP3[ic]+*cSP3[ic]+*cP4[ic]+*cM4[ic];
    }
  }
  cutMatrix.Print();

  cutAll = (*cA[0])||(*cA[1])||(*cA[2])||(*cA[3])+"abs(mag)<0.01&&ncont>0&&p.fNPoints>120";

  delete  phisP3; // = new TH1F("hhisP3","hhisP3",100,-0.01,0.01);
  delete  phisSP3; // = new TH1F("hhisSP3","hhisSP3",100,-0.001,0.001);
  delete  phisP4;// = new TH1F("hhisP4","hhisP4",100,-0.1,0.1);
  delete  phisM4;// = new TH1F("hhisM4","hhisM4",100,-0.01,0.01);



  chainPoints->Draw(">>listEL",cutAll,"entryList");
  elist = (TEntryList*)gDirectory->Get("listEL");
  chainPoints->SetEntryList(elist);
  elist->SetDirectory(0);
}



AliTPCkalmanFit * FitPointsLinear(Int_t maxTracks, Int_t trackDump){
  // create debug streamers

  TTreeSRedirector *pcstream      = new TTreeSRedirector("kalmanfitTPC.root");  
  TTreeSRedirector *pcstreamOrig = new TTreeSRedirector("kalmanfitTPCOrig.root");  
  pcstream->GetFile()->cd();
  cutMatrix.Write("cutMarix");
  elist->Write("eventList");
  //
  //
  AliTrackPointArray *pointsNS=0;
  Float_t mag=0;
  Int_t   time=0;
  AliExternalTrackParam *param0=0;
  AliExternalTrackParam *param1=0;
  chainPoints->SetBranchAddress("p.",&pointsNS);
  chainPoints->SetBranchAddress("p0In.",&param0);
  chainPoints->SetBranchAddress("p1In.",&param1);
  chainPoints->SetBranchAddress("mag",&mag);
  chainPoints->SetBranchAddress("time",&time);
  Int_t accepted=0;
  printf("\n*\n*\n*Selected entries = %d\n*\n*\n*",Int_t(elist->GetN()));

  //
  for (Int_t itrack=0;itrack<elist->GetN(); itrack++){
    if (itrack%toSkipTrack!=toSkipTrackOffset) continue;   
    Int_t entry=chainPoints->GetEntryNumber(itrack);
    chainPoints->GetEntry(entry);
    if (accepted>maxTracks) break;
    //
    AliTrackPointArray *points = AliTPCkalmanFit::SortPoints(*pointsNS);    
    if (kalmanFitApply)  kalmanFitApply->ApplyCalibration(points,-1.);
    //
    // estimate and filter scattering
    //
    TVectorD *vecRMS09 = EstimateScatteringKalmanLinear(*points,*param0,*param1,pcstream);
    if (!vecRMS09) continue;
    Bool_t isOK=kTRUE;
    if ((*vecRMS09)[0] >rmsCut09[0]) isOK=kFALSE;
    if ((*vecRMS09)[1]/(*vecRMS09)[0] >rmsCut09[1]) isOK=kFALSE;;
    if ((*vecRMS09)[2] >rmsCut09[2]) isOK=kFALSE;
    if ((*vecRMS09)[3] >rmsCut09[3]) isOK=kFALSE;
    if ((*vecRMS09)[4]/(*vecRMS09)[0] >rmsCut09[4]) isOK=kFALSE;
    if ((*vecRMS09)[5] >rmsCut09[5]) isOK=kFALSE;
    if (!isOK || isFilterTest) {
      delete points;
      continue;
    }
    kalmanFitNew->PropagateTime(time);
    //
    //
    for (Int_t idir=-1; idir<=1; idir++){
      AliTrackPointArray *pointsF = FilterPoints(*points,idir, pcstream);      
      if (!pointsF) continue;
      AliTrackPointArray *spointsF = 0;     
      // we skip points for alignemnt  but not for QA
      if (idir==0)  spointsF = SkipPoints(*pointsF, toSkip*2, toSkipOffset);
      if (idir!=0)  spointsF = SkipPoints(*pointsF, toSkip,   toSkipOffset);
      if (idir==0)  accepted++;
      //
      if (accepted%50==0) {
	kalmanFitNew->FitTrackLinear(*pointsF, pcstream);
      }else{
	if (idir==0) kalmanFitNew->FitTrackLinear(*spointsF, 0);
	if (idir!=0) kalmanFitNew->FitTrackLinear(*spointsF, 0);
      }    
      if (idir==0) kalmanFitNew->DumpTrackLinear(*pointsF,pcstream);
      if (idir==0) kalmanFitOrig->DumpTrackLinear(*pointsF,pcstreamOrig);
      if (accepted%trackDump==0) {
	printf("%d\n", accepted);
      }
      AliSysInfo::AddStamp("trackFit", accepted,itrack);
      delete pointsF;
      delete spointsF;
    }
    delete points;
  }
  pcstream->GetFile()->cd();
  kalmanFitNew->Write("kalmanFit");
  pcstreamOrig->GetFile()->cd();
  kalmanFitOrig->Write("kalmanFitOrig");
  pcstreamOrig->GetFile()->cd();
  if (kalmanFitApply) kalmanFitApply->Write("kalmanFitApply");
  
  delete pcstream;  
  delete pcstreamOrig;  
  return kalmanFitNew;   
}

void  QAPointsLinear(Int_t maxTracks, Int_t trackDump){
  /// check  the consistency of kalman fit
  /// Apply transformation

  // create debug streeamers
  TTreeSRedirector *pcstreamNonCalib      = new TTreeSRedirector("kalmanfitTPCQANonCalib.root");
  TTreeSRedirector *pcstreamCalib         = new TTreeSRedirector("kalmanfitTPCQACalib.root");
  TTreeSRedirector *pcstream=0;
  AliTPCkalmanFit  *kalmanFitters[6]={0,0,0,0,0,0};
  for (Int_t i=0;i<6;i++){
    kalmanFitters[i]=SetupFit();
  }
  //
  AliTrackPointArray *points=0;
  AliExternalTrackParam *param0=0;
  AliExternalTrackParam *param1=0;
  Float_t mag=0;
  Int_t   time=0;
  chainPoints->SetBranchAddress("p.",&points);
  chainPoints->SetBranchAddress("mag",&mag);
  chainPoints->SetBranchAddress("time",&time);
  chainPoints->SetBranchAddress("p0In.",&param0);
  chainPoints->SetBranchAddress("p1In.",&param1);

  Int_t accepted=0;
  //


  for (Int_t itrack=0;itrack<elist->GetN(); itrack++){
    if (itrack%toSkipTrack!=toSkipTrackOffset) continue;   
    Int_t entry=chainPoints->GetEntryNumber(itrack);
    chainPoints->GetEntry(entry);
    if (accepted>maxTracks) break;
    //

    AliTrackPointArray pointsCalib(*points); 
    for (Int_t iscalib=0; iscalib<1;iscalib++){
      if (iscalib>0)  kalmanFitNew->ApplyCalibration(&pointsCalib,-1.);
      if (iscalib==0) pcstream=pcstreamNonCalib;
      if (iscalib>0)  pcstream=pcstreamCalib;
      for (Int_t idir=-1; idir<=1; idir++){
	AliTrackPointArray *pointsF = FilterPoints(pointsCalib,idir, pcstream);
	if (!pointsF) continue;
	//
	if (idir==0) accepted++;
	kalmanFitters[iscalib*3+idir+1]->DumpTrackLinear(*pointsF,0);
	EstimateScatteringKalmanLinear(*pointsF,*param0,*param1,pcstream);
	delete pointsF;
      }
    }
    if (accepted%trackDump==0) {
      printf("%d\n", accepted);
    }
  }
  pcstreamCalib->GetFile()->cd();
  kalmanFitters[0]->Write("fitUpNonCalib");
  kalmanFitters[1]->Write("fitUpDownNonCalib");
  kalmanFitters[2]->Write("fitDownNonCalib");
  kalmanFitters[3]->Write("fitUpCalib");
  kalmanFitters[4]->Write("fitUpDownCalib");
  kalmanFitters[5]->Write("fitDownCalib");
  delete pcstreamCalib;  
  delete pcstreamNonCalib;  
}


void  TestScattering(Int_t maxTracks, Int_t trackDump){
  /// test Multiple scattering algorithm
  /// Apply transformation
  ///
  /// create debug streeamers

  TTreeSRedirector *pcstream      = new TTreeSRedirector("kalmanfitTPCMS.root");
  //
  //
  AliTrackPointArray *points=0;
  AliExternalTrackParam *param0=0;
  AliExternalTrackParam *param1=0;
  Float_t mag=0;
  Int_t   time=0;
  chainPoints->SetBranchAddress("p.",&points);
  chainPoints->SetBranchAddress("mag",&mag);
  chainPoints->SetBranchAddress("time",&time);
  chainPoints->SetBranchAddress("p0In.",&param0);
  chainPoints->SetBranchAddress("p1In.",&param1);
  Int_t accepted=0;
  //
  for (Int_t itrack=0;itrack<elist->GetN(); itrack++){
    Int_t entry=chainPoints->GetEntryNumber(itrack);
    chainPoints->GetEntry(entry);
    if (accepted>maxTracks) break;
    //
    AliTrackPointArray *pointsSorted = AliTPCkalmanFit::SortPoints(*points);
    EstimateScatteringKalmanLinear(*pointsSorted,*param0,*param1,pcstream);
    accepted++;
    if (accepted%trackDump==0) {
      printf("%d\n", accepted);
    }
    delete pointsSorted;
  }
  delete pcstream;  
}





AliTrackPointArray *SkipPoints(AliTrackPointArray &points, Int_t nskip, Int_t nskipOffset){
  /// create new array with skipped points

  Int_t npoints = points.GetNPoints();
  Int_t npointsF = (npoints-nskipOffset-1)/nskip;
  AliTrackPoint point;
  AliTrackPointArray *pointsF= new AliTrackPointArray(npointsF);
  Int_t used=0;
  for (Int_t ipoint=nskipOffset; ipoint<npoints; ipoint+=nskip){
    //
    if (!points.GetPoint(point,ipoint)) continue;
    pointsF->AddPoint(used,&point);
    used++;
    if (used==npointsF) break;
  }
  return pointsF;
}



AliTrackPointArray *FilterPoints(AliTrackPointArray &points, Int_t dir, TTreeSRedirector *pcstream){
  ///  Filter points - input points for KalmanFilter

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
  if (TMath::Sqrt(chi2Y)>krmsYcutGlobal) return 0;
  if (TMath::Sqrt(chi2Y)>krmsYcutGlobal) return 0;
  //
  //
  Int_t accepted=0;
  Int_t toBeUsed    =0;
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
      accepted++;
      if (pointsF) pointsF->AddPoint(toBeUsed,&point);
      toBeUsed++;
    }
    if (pcstream){
      (*pcstream)<<"filter"<<
	"iter="<<iter<<
	"accepted="<<accepted<<
	"minX="<<minX<<
	"maxX="<<maxX<<
	"vY.="<<&vecY<<
	"vZ.="<<&vecZ<<
	"chi2Y="<<chi2Y<<
	"chi2Z="<<chi2Z<<
	"\n";
    }
    if (accepted<knclCut) break;
    if (iter==0) pointsF = new AliTrackPointArray(toBeUsed);
    accepted=0;
    toBeUsed=0;
  }
  return pointsF;
}


AliTrackPointArray * SortPoints(AliTrackPointArray &points){
  /// Creates the array  - points sorted according radius - neccessay for kalman fit
  ///
  /// 0. choose the frame - rotation angle

  Int_t npoints = points.GetNPoints();
  if (npoints<1) return 0;
  Double_t currentAlpha = TMath::ATan2(points.GetY()[npoints-1]-points.GetY()[0], points.GetX()[npoints-1]-points.GetX()[0]);  
  Double_t ca = TMath::Cos(currentAlpha);
  Double_t sa = TMath::Sin(currentAlpha);
  //
  // 1. sort the points
  //
  Double_t *rxvector = new Double_t[npoints];
  Int_t    *indexes  = new Int_t[npoints];
  for (Int_t ipoint=0; ipoint<npoints-1; ipoint++){
    rxvector[ipoint]=ca*points.GetX()[ipoint]+sa*points.GetY()[ipoint];
  }
  TMath::Sort(npoints, rxvector,indexes,kFALSE);
  AliTrackPoint point;
  AliTrackPointArray *pointsSorted= new AliTrackPointArray(npoints);
  for (Int_t ipoint=0; ipoint<npoints; ipoint++){
    if (!points.GetPoint(point,indexes[ipoint])) continue;
    pointsSorted->AddPoint(ipoint,&point);
  }
  delete [] rxvector;
  delete [] indexes;
  return pointsSorted;
}

TVectorD *  EstimateScatteringKalmanLinear(AliTrackPointArray &points, AliExternalTrackParam &p0, AliExternalTrackParam &p1 , TTreeSRedirector *pcstream){
  /// Algorithm - 0. Fit the track forward and backward
  ///           - 1. Store the current parameters in each point

  const Int_t kMinPoints= 70;
  const Double_t kResY  = 0.1;
  const Double_t kResZ  = 0.1;
  const Double_t kMisY  = 0.02;
  const Double_t kMisZ  = 0.02;
  const Double_t kLArm  = 120.;
  const Double_t kACsideFac = 10.;
  const Double_t kMS0   = kResY/20.;
    
  Int_t npoints = points.GetNPoints();
  if (npoints<kMinPoints) return 0;

  TVectorD *vecPos[11]={0,0,0,0,0,0,0,0,0,0,0};
  for (Int_t i=0;i<11;i++){
    vecPos[i]=new TVectorD(npoints);
  }

  //
  Double_t currentAlpha = TMath::ATan2(points.GetY()[npoints-1]-points.GetY()[0], points.GetX()[npoints-1]-points.GetX()[0]);  
  Double_t ca  = TMath::Cos(currentAlpha);
  Double_t sa  = TMath::Sin(currentAlpha);
  //
  //
  TMatrixD trParamY(2,1),trCovarY(2,2);
  TMatrixD trParamZ(2,1),trCovarZ(2,2);
  Int_t kNmeas  = 1; 
  Int_t nelem   = 2;
  //
  TMatrixD matHk(kNmeas,nelem);     // vector to mesurement
  TMatrixD matHkT(nelem,kNmeas);    // helper matrix Hk transpose
  TMatrixD matFk(nelem,nelem);      
  TMatrixD matFkT(nelem,nelem);      
  TMatrixD vecYk(kNmeas,1);         // Innovation or measurement residual
  TMatrixD vecZk(kNmeas,1);         // Innovation or measurement residual
  TMatrixD measR(kNmeas,kNmeas);
  TMatrixD matSk(kNmeas,kNmeas);    // Innovation (or residual) covariance
  TMatrixD matKk(nelem,kNmeas);     // Optimal Kalman gain
  TMatrixD covXk2(nelem,nelem);     // helper matrix
  TMatrixD covXk3(nelem,nelem);     // helper matrix
  TMatrixD mat1(nelem,nelem);
  mat1(0,0)=1.; mat1(0,1)=0.;
  mat1(1,0)=0.; mat1(1,1)=1.;
  matHk(0,0)=1;
  matHk(0,1)=0;

  Double_t lastX    = 0;
  Int_t lastVolId=-1;
  for (Int_t idir=0; idir<2;idir++){
    // fit direction
    //
    for (Int_t ip=0; ip<npoints; ip++){
      Int_t ipoint= ip;
      if (idir>0) ipoint = npoints-ip-1;
      Double_t rx =   ca*points.GetX()[ipoint]+sa*points.GetY()[ipoint];
      Double_t ry =  -sa*points.GetX()[ipoint]+ca*points.GetY()[ipoint];
      Double_t rz =  points.GetZ()[ipoint];
      Int_t  volId= points.GetVolumeID()[ipoint];
      //
      if (ip==0){
	// set initital parameters and covariance - use first and middle point
	Double_t rxm =   ca*points.GetX()[npoints/2]+sa*points.GetY()[npoints/2];
	Double_t rym =  -sa*points.GetX()[npoints/2]+ca*points.GetY()[npoints/2];
	Double_t rzm =  points.GetZ()[npoints/2];
	trParamY(0,0) = ry; 
	trParamY(1,0) = (rym-ry)/(rxm-rx);
	trParamZ(0,0) = rz;
	trParamZ(1,0) = (rzm-rz)/(rxm-rx);
	//
	trCovarY(0,0) = kResY*kResY;
	trCovarY(1,1) = (kResY*kResY)/((rxm-rx)*(rxm-rx));
	trCovarZ(0,0) = kResZ*kResZ;
	trCovarZ(1,1) = (kResZ*kResZ)/((rxm-rx)*(rxm-rx));
	lastX     = rx;
	lastVolId = volId;
      }
      //
      // Propagate
      //
      if ((volId%36)<18 && (lastVolId%36)>=18){
	// A - C side cross
	trCovarY(0,0)+=kMisY*kMisY*kACsideFac;
	trCovarZ(0,0)+=kMisZ*kMisZ*kACsideFac;
	trCovarY(1,1)+=kACsideFac*(kMisY*kMisY)/(kLArm*kLArm);;
	trCovarZ(1,1)+=kACsideFac*(kMisZ*kMisZ)/(kLArm*kLArm);
	lastVolId=volId;
      }
      if (volId!=lastVolId){
	// volumeID change
	trCovarY(0,0)+=kMisY*kMisY;
	trCovarZ(0,0)+=kMisZ*kMisZ;
	trCovarY(1,1)+=(kMisY*kMisY)/(kLArm*kLArm);;
	trCovarZ(1,1)+=(kMisZ*kMisZ)/(kLArm*kLArm);
	lastVolId=volId;
      }
      //
      // Propagate
      //
      Double_t deltaX=rx-lastX;
      trParamY(0,0)+=deltaX*trParamY(1,0);
      trParamZ(0,0)+=deltaX*trParamZ(1,0);
      matFk(0,0)=1; matFk(0,1)=deltaX;
      matFk(1,0)=0; matFk(1,1)=1;
      matFkT=matFk.T(); matFk.T();
      covXk2=matFk*trCovarY*matFkT;
      trCovarY=covXk2;
      covXk2=matFk*trCovarZ*matFkT;
      trCovarZ=covXk2;

      // multiple scattering
      trCovarY(1,1)+=TMath::Abs(deltaX)*kMS0*kMS0;
      trCovarZ(1,1)+=TMath::Abs(deltaX)*kMS0*kMS0;
      lastX=rx;
      //
      // Update
      //
      for (Int_t coord=0; coord<2;coord++){
	TMatrixD* pvecXk = (coord==0)? &trParamY: &trParamZ;
	TMatrixD* pcovXk = (coord==0)? &trCovarY: &trCovarZ;
	//
	TMatrixD& vecXk = *pvecXk;
	TMatrixD& covXk = *pcovXk;
	measR(0,0) = (coord==0) ? kResY:kResZ;
	vecZk(0,0) = (coord==0) ? ry:rz;
	//
	vecYk = vecZk-matHk*vecXk;               // Innovation or measurement residual
	matHkT=matHk.T(); matHk.T();
	matSk = (matHk*(covXk*matHkT))+measR;    // Innovation (or residual) covariance
	matSk.Invert();
	matKk = (covXk*matHkT)*matSk;            //  Optimal Kalman gain
	//
	covXk2= (mat1-(matKk*matHk));
	covXk3 =  covXk2*covXk;          
	vecXk += matKk*vecYk;                    //  updated vector 
	covXk = covXk3; 
      }
      //
      // store parameters
      //
      (*vecPos[0])[ipoint]=rx;
      (*vecPos[1])[ipoint]=ry;
      (*vecPos[2])[ipoint]=rz;

      (*vecPos[4*idir+0+3])[ipoint]=trParamY(0,0);
      (*vecPos[4*idir+1+3])[ipoint]=trParamY(1,0);
      //
      (*vecPos[4*idir+2+3])[ipoint]=trParamZ(0,0);
      (*vecPos[4*idir+3+3])[ipoint]=trParamZ(1,0);
    }
  }
  //
  //
  //
  TVectorD vec(npoints);
  TVectorD rms(6);
  TVectorD rms09(6); // robust RMS - fraction 0.9
  TVectorD mean09(6); // robust RMS - fraction 0.9
  Double_t meanR,rmsR;
  Int_t npoints09 = Int_t(npoints*0.9);
  //
  vec=(*(vecPos[3])-*(vecPos[1]));
  rms[0]=TMath::RMS(npoints, vec.GetMatrixArray());   // rms cluster y
  AliMathBase::EvaluateUni(npoints, vec.GetMatrixArray(),meanR,rmsR,npoints09);
  rms09[0]=rmsR;
  mean09[0]=meanR;
  vec=(*(vecPos[7])-*(vecPos[3]));
  rms[1]=TMath::RMS(npoints, vec.GetMatrixArray());   // rms track y
  AliMathBase::EvaluateUni(npoints, vec.GetMatrixArray(),meanR,rmsR,npoints09);
  rms09[1]=rmsR;
  mean09[1]=meanR;
  vec=(*(vecPos[8])-*(vecPos[4]));
  rms[2]=TMath::RMS(npoints, vec.GetMatrixArray());   // rms track ky
  AliMathBase::EvaluateUni(npoints, vec.GetMatrixArray(),meanR,rmsR,npoints09);
  rms09[2]=rmsR;
  mean09[2]=meanR;
  //
  vec=(*(vecPos[5])-*(vecPos[2]));
  rms[3]=TMath::RMS(npoints, vec.GetMatrixArray());   // rms cluster z
  AliMathBase::EvaluateUni(npoints, vec.GetMatrixArray(),meanR,rmsR,npoints09);
  rms09[3]=rmsR;
  mean09[3]=meanR;
  vec=(*(vecPos[9])-*(vecPos[5]));
  rms[4]=TMath::RMS(npoints, vec.GetMatrixArray());   // rms track z
  AliMathBase::EvaluateUni(npoints, vec.GetMatrixArray(),meanR,rmsR,npoints09);
  rms09[4]=rmsR;
  mean09[4]=meanR;
  vec=(*(vecPos[10])-*(vecPos[6]));
  rms[5]=TMath::RMS(npoints, vec.GetMatrixArray());   // rms track kz
  AliMathBase::EvaluateUni(npoints, vec.GetMatrixArray(),meanR,rmsR,npoints09);
  rms09[5]=rmsR;
  mean09[5]=meanR;


  (*pcstream)<<"kf"<<
    "p.="<<&points<<    
    "p0.="<<&p0<<
    "p1.="<<&p1<<
    "rms.="<<&rms<<
    "rms09.="<<&rms09<<
    "mean09.="<<&mean09<<
    "py.="<<&trParamY<<
    "pz.="<<&trParamZ<<
    "cy.="<<&trCovarY<<
    "cz.="<<&trCovarY<<
    "\n";
  for (Int_t i=0;i<11;i++){
    delete vecPos[i];
  }
  return new TVectorD(rms09);
}









void  AddFitFieldCage(AliTPCkalmanFit *kalmanFit){
  /// Add radial scaling due field cage

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
  /// Add linear local phi scaling -
  /// separate IROC/OROC  - A side/C side
  ///    "tscalingLocalPhiIROC"
  ///    "tscalingLocalPhiOROC"

  TBits maskInner(72);
  TBits maskOuter(72);
  for (Int_t i=0;i<72;i++){
    if (i<36){
      maskInner[i]=kTRUE;
    }
    if (i>=36){
      maskOuter[i]=kTRUE;
    }
  }
  TVectorD fpar(10);
  AliTPCTransformation * transformation=0;
  fpar[0]=1; 
  transformation = new AliTPCTransformation("tscalingLocalPhiIROC", new TBits(maskInner), 0,"TPCscalingPhiLocal",0,  1);
  transformation->SetParams(0,0.02,0,&fpar);
  kalmanFit->AddCalibration(transformation);  
  transformation = new AliTPCTransformation("tscalingLocalPhiOROC", new TBits(maskOuter), 0,"TPCscalingPhiLocal",0,  1);
  transformation->SetParams(0,0.02,0,&fpar);
  kalmanFit->AddCalibration(transformation);  
  //
}

void AddDrift(AliTPCkalmanFit *kalmanFit){
  /// Add drift velocity transformation

  TVectorD fpar(10);
  AliTPCTransformation * transformation=0;
  fpar[0]=1;
  transformation = new AliTPCTransformation("tTPCscalingZDrift", AliTPCTransformation::BitsAll(), 0,0,"TPCscalingZDr",  0);
  transformation->SetParams(0,4.,1.,&fpar);
  kalmanFit->AddCalibration(transformation);  
  //
  transformation = new AliTPCTransformation("tTPCscalingZDriftGy", AliTPCTransformation::BitsAll(), 0,0,"TPCscalingZDrGy",  0);
  transformation->SetParams(0,0.2,0.0,&fpar);
  kalmanFit->AddCalibration(transformation);  
}




void  AddZShift(AliTPCkalmanFit *kalmanFit, Int_t ncos, Int_t nsin){
  ///

  TVectorD fpar(10);
  fpar[0]=0; fpar[1]=0; fpar[2]=0;
  char tname[1000];
  AliTPCTransformation * transformation=0;
  //
  //
  //
  for (Int_t i=0; i<=ncos;i++){
    fpar[0]=i;  // cos frequency
    fpar[1]=0;  // no sinus
    fpar[2]=1;  // relative misalignment
    //
    sprintf(tname,"tTPCDeltaZIROCOROCSideA_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), 0,0, "TPCDeltaZ" ,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    sprintf(tname,"tTPCDeltaZIROCOROCSideC_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), 0,0, "TPCDeltaZ" ,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    //
    //
    fpar[2]=0;  // absolute  misalignment
    sprintf(tname,"tTPCDeltaZSectorSideA_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), 0,0, "TPCDeltaZ" ,  0);
    transformation->SetParams(0,0.1,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    sprintf(tname,"tTPCDeltaZSectorSideC_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), 0,0, "TPCDeltaZ" ,  0);
    transformation->SetParams(0,0.1,0,&fpar);
    kalmanFit->AddCalibration(transformation);
  }

  for (Int_t i=1; i<=nsin;i++){
    fpar[0]=0;  // cos frequency
    fpar[1]=i;  // sinus frequncy
    fpar[2]=1;  // relative misalignment
    //
    sprintf(tname,"tTPCDeltaZIROCOROCSideA_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), 0,0, "TPCDeltaZ" ,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    sprintf(tname,"tTPCDeltaZIROCOROCSideC_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), 0,0, "TPCDeltaZ" ,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    //
    //
    fpar[2]=0;  // absolute  misalignment
    sprintf(tname,"tTPCDeltaZSectorSideA_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), 0,0, "TPCDeltaZ" ,  0);
    transformation->SetParams(0,0.1,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    sprintf(tname,"tTPCDeltaZSectorSideC_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), 0,0, "TPCDeltaZ" ,  0);
    transformation->SetParams(0,0.1,0,&fpar);
    kalmanFit->AddCalibration(transformation);
  }

}




void AddZTilting(AliTPCkalmanFit *kalmanFit, Int_t ncos, Int_t nsin){
  /// z tilting absolute (sector) and relative (IROC-OROC)

  TVectorD fpar(10);
  fpar[0]=0; fpar[1]=0; fpar[2]=0;
  char tname[1000];
  AliTPCTransformation * transformation=0;
  //
  //
  //
  for (Int_t i=0; i<=ncos;i++){
    fpar[0]=i;  // cos frequency
    fpar[1]=0;  // sinus frequency
    fpar[2]=1;  // relative misalignment
    //
    sprintf(tname,"tTPCTiltingZIROCOROCSideA_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), 0,0, "TPCTiltingZ" ,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    sprintf(tname,"tTPCTiltingZIROCOROCSideC_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), 0,0, "TPCTiltingZ" ,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    //
    //
    fpar[2]=0;  // absolute  misalignment
    sprintf(tname,"tTPCTiltingZSectorSideA_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), 0,0, "TPCTiltingZ" ,  0);
    transformation->SetParams(0,0.1,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    sprintf(tname,"tTPCTiltingZSectorSideC_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), 0,0, "TPCTiltingZ" ,  0);
    transformation->SetParams(0,0.1,0,&fpar);
    kalmanFit->AddCalibration(transformation);
  }

  for (Int_t i=1; i<=nsin;i++){
    fpar[0]=0;  // cos frequency
    fpar[1]=i;  // sinus frequncy
    fpar[2]=1;  // relative misalignment
    //
    sprintf(tname,"tTPCTiltingZIROCOROCSideA_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), 0,0, "TPCTiltingZ" ,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    sprintf(tname,"tTPCTiltingZIROCOROCSideC_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), 0,0, "TPCTiltingZ" ,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    //
    //
    fpar[2]=0;  // absolute  misalignment
    sprintf(tname,"tTPCTiltingZSectorSideA_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), 0,0, "TPCTiltingZ" ,  0);
    transformation->SetParams(0,0.1,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    sprintf(tname,"tTPCTiltingZSectorSideC_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), 0,0, "TPCTiltingZ" ,  0);
    transformation->SetParams(0,0.1,0,&fpar);
    kalmanFit->AddCalibration(transformation);
  }
}



void  AddLocalXYMisalignment(AliTPCkalmanFit *kalmanFit){
  ///

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

void  AddLocalXYMisalignmentSector(AliTPCkalmanFit *kalmanFit){
  ///

  TVectorD fpar(10);
  AliTPCTransformation * transformation=0;
  Int_t fixSector =4;
  //
  for (Int_t isec=0; isec<36;isec++){
    TBits mask(72);
    mask[isec]=kTRUE;
    char tname[1000];
    //
    sprintf(tname,"tTPClocalLxIROC%d",isec);
    transformation = new AliTPCTransformation(tname, new TBits(mask), "TPClocaldLydGX","TPClocaldLydGY",0,  0);
    transformation->SetParams(0,0.2,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    sprintf(tname,"tTPClocalLyIROC%d",isec);
    transformation = new AliTPCTransformation(tname, new TBits(mask), "TPClocaldLydGX","TPClocaldLydGY",0,  0);
    transformation->SetParams(0,0.2,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    
    sprintf(tname,"tTPClocalRzIROC%d",isec);
    transformation = new AliTPCTransformation(tname, new TBits(mask), "TPClocaldLydGX","TPClocaldRzdGY",0,  0);
    transformation->SetParams(0,0.2,0,&fpar);
    kalmanFit->AddCalibration(transformation);

  }
  //
  for (Int_t isec=0; isec<36;isec++){
    if (isec%18==fixSector) continue;
    TBits mask(72);
    mask[isec]   =kTRUE;
    mask[isec+36]=kTRUE;    
    char tname[1000];
    //
    sprintf(tname,"tTPClocalLxSector%d",isec);
    transformation = new AliTPCTransformation(tname, new TBits(mask), "TPClocaldLydGX","TPClocaldLydGY",0,  0);
    transformation->SetParams(0,0.2,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    sprintf(tname,"tTPClocalLySector%d",isec);
    transformation = new AliTPCTransformation(tname, new TBits(mask), "TPClocaldLydGX","TPClocaldLydGY",0,  0);
    transformation->SetParams(0,0.2,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    
    sprintf(tname,"tTPClocalRzSector%d",isec);
    transformation = new AliTPCTransformation(tname, new TBits(mask), "TPClocaldLydGX","TPClocaldRzdGY",0,  0);
    transformation->SetParams(0,0.2,0,&fpar);
    kalmanFit->AddCalibration(transformation);
  }
  

}


void  AddAlignOROCIROCFourier(AliTPCkalmanFit *kalmanFit, Int_t ncos, Int_t nsin){
  ///

  TVectorD fpar(10);
  AliTPCTransformation * transformation=0;

  for (Int_t i=0; i<=ncos;i++){
    char tname[1000];
    fpar[0]=i;  // cos frequency
    fpar[1]=0;  // no sinus
    fpar[2]=1;  // relative misalignment
    //
    // Local X shift
    //
    sprintf(tname,"tTPClocalLxIROCOROCSideA_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), "TPClocaldLxdGX","TPClocaldLxdGY",0,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    sprintf(tname,"tTPClocalLxIROCOROCSideC_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), "TPClocaldLxdGX","TPClocaldLxdGY",0,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    //
    // Local Y shift
    //
    sprintf(tname,"tTPClocalLyIROCOROCSideA_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), "TPClocaldLydGX","TPClocaldLydGY",0,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    sprintf(tname,"tTPClocalLyIROCOROCSideC_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), "TPClocaldLydGX","TPClocaldLydGY",0,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    //
    //
    // Z rotation
    //
    sprintf(tname,"tTPClocalRzIROCOROCSideA_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), "TPClocaldRzdGX","TPClocaldRzdGY",0,  0);
    transformation->SetParams(0,0.3,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    sprintf(tname,"tTPClocalRzIROCOROCSideC_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), "TPClocaldRzdGX","TPClocaldRzdGY",0,  0);
    transformation->SetParams(0,0.3,0,&fpar);
    kalmanFit->AddCalibration(transformation);
  }
  //
  for (Int_t i=1; i<=nsin;i++){
    char tname[1000];
    fpar[0]=0;  // cos frequency
    fpar[1]=i;  // sinus frequency
    fpar[2]=1;  // relative misalignment
    //
    // Local X shift
    //
    sprintf(tname,"tTPClocalLxIROCOROCSideA_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), "TPClocaldLxdGX","TPClocaldLxdGY",0,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    sprintf(tname,"tTPClocalLxIROCOROCSideC_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), "TPClocaldLxdGX","TPClocaldLxdGY",0,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    //
    // Local Y shift
    //
    sprintf(tname,"tTPClocalLyIROCOROCSideA_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), "TPClocaldLydGX","TPClocaldLydGY",0,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    sprintf(tname,"tTPClocalLyIROCOROCSideC_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), "TPClocaldLydGX","TPClocaldLydGY",0,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    //
    //
    // Z rotation
    //
    sprintf(tname,"tTPClocalRzIROCOROCSideA_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), "TPClocaldRzdGX","TPClocaldRzdGY",0,  0);
    transformation->SetParams(0,0.3,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    sprintf(tname,"tTPClocalRzIROCOROCSideC_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), "TPClocaldRzdGX","TPClocaldRzdGY",0,  0);
    transformation->SetParams(0,0.3,0,&fpar);
    kalmanFit->AddCalibration(transformation);
  }
}

void  AddAlignSectorFourier(AliTPCkalmanFit *kalmanFit, Int_t ncos, Int_t nsin){
  ///

  TVectorD fpar(10);
  AliTPCTransformation * transformation=0;

  for (Int_t i=0; i<=ncos;i++){
    char tname[1000];
    fpar[0]=i;  // cos frequency
    fpar[1]=0;  // no sinus
    fpar[2]=0;  // absolute misalignment
    if (i>0){
      //
      // A side is reference
      // local x 
      sprintf(tname,"tTPClocalLxSectorSideA_Cos%d",i);
      transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), "TPClocaldLxdGX","TPClocaldLxdGY",0,  0);
      transformation->SetParams(0,0.03,0,&fpar);
      kalmanFit->AddCalibration(transformation);
      if (i>1){
	//
	// Local Y shift
	//
	sprintf(tname,"tTPClocalLySectorSideA_Cos%d",i);
	transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), "TPClocaldLydGX","TPClocaldLydGY",0,  0);
	transformation->SetParams(0,0.03,0,&fpar);
	kalmanFit->AddCalibration(transformation);
      }
      //
      //
    }
    //
    // C side to vary
    // local x 
    sprintf(tname,"tTPClocalLxSectorSideC_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), "TPClocaldLxdGX","TPClocaldLxdGY",0,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    // Local Y shift
    //
    sprintf(tname,"tTPClocalLySectorSideC_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), "TPClocaldLydGX","TPClocaldLydGY",0,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //    
    //
    // Z rotation - independent
    //
    sprintf(tname,"tTPClocalRzSectorSideA_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), "TPClocaldRzdGX","TPClocaldRzdGY",0,  0);
    transformation->SetParams(0,0.3,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    sprintf(tname,"tTPClocalRzSectorSideC_Cos%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), "TPClocaldRzdGX","TPClocaldRzdGY",0,  0);
    transformation->SetParams(0,0.3,0,&fpar);
    kalmanFit->AddCalibration(transformation);
  }




  //
  //
  //
  for (Int_t i=1; i<=nsin;i++){
    char tname[1000];
    fpar[0]=0;  // non cos frequency
    fpar[1]=i;  // sinus frequency
    fpar[2]=0;  // absolute misalignment
    //
    // Local X shift
    //
    sprintf(tname,"tTPClocalLxSectorSideA_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), "TPClocaldLxdGX","TPClocaldLxdGY",0,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    sprintf(tname,"tTPClocalLxSectorSideC_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), "TPClocaldLxdGX","TPClocaldLxdGY",0,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    //
    // Local Y shift
    //
    sprintf(tname,"tTPClocalLySectorSideA_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), "TPClocaldLydGX","TPClocaldLydGY",0,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    sprintf(tname,"tTPClocalLySectorSideC_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), "TPClocaldLydGX","TPClocaldLydGY",0,  0);
    transformation->SetParams(0,0.03,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    //
    //
    //
    // Z rotation
    //
    sprintf(tname,"tTPClocalRzSectorSideA_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(0), "TPClocaldRzdGX","TPClocaldRzdGY",0,  0);
    transformation->SetParams(0,0.3,0,&fpar);
    kalmanFit->AddCalibration(transformation);
    sprintf(tname,"tTPClocalRzSectorSideC_Sin%d",i);
    transformation = new AliTPCTransformation(tname, AliTPCTransformation::BitsSide(1), "TPClocaldRzdGX","TPClocaldRzdGY",0,  0);
    transformation->SetParams(0,0.3,0,&fpar);
    kalmanFit->AddCalibration(transformation);
  }
}

void SelectPixel(){
  for (Int_t i=0;i<8;i++){
    kalmanFitNew->fLinearTrackDelta[i]->GetAxis(2)->SetRangeUser(-10,10);
    kalmanFitOrig->fLinearTrackDelta[i]->GetAxis(2)->SetRangeUser(-10,10);
    kalmanFitNew->fLinearTrackDelta[i]->GetAxis(3)->SetRangeUser(-15,15);
    kalmanFitOrig->fLinearTrackDelta[i]->GetAxis(3)->SetRangeUser(-15,15);
  }
}

void SelectNonPixelA(){
  for (Int_t i=0;i<8;i++){
    kalmanFitNew->fLinearTrackDelta[i]->GetAxis(2)->SetRangeUser(10,250);
    kalmanFitOrig->fLinearTrackDelta[i]->GetAxis(2)->SetRangeUser(10,250);
    kalmanFitNew->fLinearTrackDelta[i]->GetAxis(3)->SetRangeUser(0,250);
    kalmanFitOrig->fLinearTrackDelta[i]->GetAxis(3)->SetRangeUser(0,250);
  }
}

void SelectNonPixelC(){
  for (Int_t i=0;i<8;i++){
    kalmanFitNew->fLinearTrackDelta[i]->GetAxis(2)->SetRangeUser(10,250);
    kalmanFitOrig->fLinearTrackDelta[i]->GetAxis(2)->SetRangeUser(10,250);
    kalmanFitNew->fLinearTrackDelta[i]->GetAxis(3)->SetRangeUser(-250,0);
    kalmanFitOrig->fLinearTrackDelta[i]->GetAxis(3)->SetRangeUser(-250,0);
  }
}


void DumpQA1D(  TObjArray &arrayOut){
  ///

  TF1 fg("fg","gaus");
  TMatrixD sideARMS(8,2);
  TMatrixD sideCRMS(8,2);
  TMatrixD sideACRMS(8,2);
  TH1 *his=0;
  //
  // A side
  //
  SelectNonPixelA();
  for (Int_t i=0; i<8;i++){
    his = kalmanFitOrig->fLinearTrackDelta[i]->Projection(0);
    his->Fit(&fg,"","", -0.15,0.15);
    sideARMS(i,0) = fg.GetParameter(2);
    his->SetDirectory(0);
    his->SetName(Form("Original SideA_%s",his->GetName()));
    his->SetTitle(Form("Original SideA_%s",his->GetTitle()));
    arrayOut.AddLast(his);
    his = kalmanFitNew->fLinearTrackDelta[i]->Projection(0);
    his->Fit(&fg,"","", -0.15,0.15);
    sideARMS(i,1) = fg.GetParameter(2);
    his->SetDirectory(0);
    his->SetName(Form("Aligned SideA_%s",his->GetName()));
    his->SetTitle(Form("Aligned SideA_%s",his->GetTitle()));
    arrayOut.AddLast(his);
  }
  //
  // C side
  //
  SelectNonPixelC();
  for (Int_t i=0; i<8;i++){
    his = kalmanFitOrig->fLinearTrackDelta[i]->Projection(0);
    his->Fit(&fg,"","", -0.15,0.15);
    sideCRMS(i,0) = fg.GetParameter(2);
    his->SetDirectory(0);
    his->SetName(Form("Original SideC_%s",his->GetName()));
    his->SetTitle(Form("Original SideC_%s",his->GetTitle()));
    arrayOut.AddLast(his);
    his = kalmanFitNew->fLinearTrackDelta[i]->Projection(0);
    his->Fit(&fg,"","", -0.15,0.15);
    sideCRMS(i,1) = fg.GetParameter(2);
    his->SetDirectory(0);
    his->SetName(Form("Aligned SideC_%s",his->GetName()));
    his->SetTitle(Form("Aligned SideC_%s",his->GetTitle()));
    arrayOut.AddLast(his);
  }
  //
  // AC side
  //
  SelectPixel();
  for (Int_t i=0; i<8;i++){
    his = kalmanFitOrig->fLinearTrackDelta[i]->Projection(0);
    his->Fit(&fg,"","", -0.15,0.15);
    sideACRMS(i,0) = fg.GetParameter(2);
    his->SetDirectory(0);
    his->SetName(Form("Original SideAC_%s",his->GetName()));
    his->SetTitle(Form("Original SideAC_%s",his->GetTitle()));
    arrayOut.AddLast(his);
    his = kalmanFitNew->fLinearTrackDelta[i]->Projection(0);
    his->Fit(&fg,"","", -0.15,0.15);
    sideACRMS(i,1) = fg.GetParameter(2);
    his->SetDirectory(0);
    his->SetName(Form("Aligned SideC_%s",his->GetName()));
    his->SetTitle(Form("Aligned SideC_%s",his->GetTitle()));
    arrayOut.AddLast(his);
  }
  printf("DumQA\n");
  sideARMS.Print();
  sideCRMS.Print();
  sideACRMS.Print();
}

void MakeFolders(TObjArray * arrayOut){
  ///

  TFolder *folderBase = new TFolder("TPC align","TPC align");
  //
  //
  TString maskDelta[8];
  TString maskAlign[2]={"Orig","Alig");
  for (Int_t i=0; i<8;i++){
    maskDelta[i] = kalmanFitNew->fLinearTrackDelta[i]->GetName();
  }
  
  Int_t entries = arrayOut->GetEntries();
  //
}




void MergeKalman(const char * list = "kalmanFit.list"){
  ///

  ifstream in;
  in.open(list);
  TString currentFile;
  kalmanFitNew= 0;
  Int_t counter=0;
  while(in.good()) {
    //
    // calibrated
    //
    in >> currentFile;    
    printf("%d\t%d\t%s\n", counter,currentFile.Length(),currentFile.Data());
    if (currentFile.Length()==0) continue;
    TFile * ffit = TFile::Open(currentFile.Data());
    TEntryList * tlist = (TEntryList*) ffit->Get("eventList");
    TMatrixD * cMatrix = (TMatrixD*) ffit->Get("cutMarix");
    if (tlist&&cMatrix){
      printf("Track entries=%d\n",tlist->GetN());
      if (cMatrix->Sum()<=0.00000001) {
	printf("Problem with track selection\n");
	continue;
      }
    }else{
      printf("Problem with track selection\n");
      continue;
    }
    //

    AliTPCkalmanFit * fit = ( AliTPCkalmanFit *)ffit->Get("kalmanFit");
    if (!fit) continue;
    if (!kalmanFitNew) {kalmanFitNew= fit; continue;};
    kalmanFitNew->Add(fit);
    printf("Selected entries=%f\n",fit->fLinearTrackDelta[0]->GetEntries());
    //delete tlist;
    //delete cMatrix;
    delete fit;
    delete ffit;
    //
    // original
    //
    currentFile.ReplaceAll("TPC","TPCOrig");
    printf("%d\t%d\t%s\n", counter,currentFile.Length(),currentFile.Data());
    if (currentFile.Length()==0) continue;
    ffit = TFile::Open(currentFile.Data());
    fit = ( AliTPCkalmanFit *)ffit->Get("kalmanFitOrig");
    if (!fit) continue;
    if (!kalmanFitOrig) {kalmanFitOrig= fit; continue;};
    kalmanFitOrig->Add(fit);
    delete fit;
    delete ffit;
    counter++;
  }
  //
  // save merged results
  //
  TFile f("mergeKalmanFit.root","recreate");
  kalmanFitNew->Write("kalmanFitNew");
  kalmanFitOrig->Write("kalmanFitOrig");
  f.Close();
}



/*
  Example - submit alignment as batch jobs
  ifile=0; 
  ntracksSkip=200
  ntracksSkipOffset=0
  nclustersSkip=2
  nclustersSkipOffset=0
  ntracks=1000000
  ndump=5
  isTest=0
  bDir=`pwd`
  ver=aaa;

  while [ $ntracksSkipOffset -lt $ntracksSkip ] ; do
    nclustersSkipOffset=0;
    while [ $nclustersSkipOffset -lt $nclustersSkip ] ; do
      echo Tr $ntracksSkipOffset  Cl $nclustersSkipOffset;     
      ver=kalmanDirTrack$ntracksSkipOffset$nclustersSkipOffset
      echo New Directory $ver
      mkdir $ver; 
      cd $ver; 
      cp $bDir/align.txt . ;  
      ln -sf $bDir/kalmanFitApply.root . ;
      echo  aliroot  -q -b  "$ALICE_ROOT/TPC/CalibMacros/CalibAlignKalman.C($ntracks,10000,0,$ndump,$ntracksSkip,$ntracksSkipOffset, $nclustersSkip,$nclustersSkipOffset,$isTest)" ; 
      bsub -q alice-t3_8h -o `pwd`/output.log command aliroot  -q -b  "$ALICE_ROOT/TPC/CalibMacros/CalibAlignKalman.C($ntracks,10000,0,$ndump,$ntracksSkip,$ntracksSkipOffset, $nclustersSkip,$nclustersSkipOffset,$isTest)"  ; 
      nclustersSkipOffset=$(( $nclustersSkipOffset + 1 ))
      cd $bDir; 
      echo $bDir; 
    done;
    cd $bDir; 
    echo $bDir; 
    ntracksSkipOffset=$(( $ntracksSkipOffset + 1 ))
    echo Tr $ntracksSkipOffset;   
done



*/











