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

//-------------------------------------------------------------------
// Base class to test the extrapolation performance from TPC to outer 
// detectors. Several member functions AddTracks() with different
// arguments can be called to execute extrapolation and several 
// THnSparse are filled with residuls and basic track information
//
// Anthor: R.Ma, M.Ivanov 02/04/2011
//--------------------------------------------------------------------
/* $Id:$ */

#include "AliTrackComparison.h"
#include "AliExternalTrackParam.h"
#include "AliTrackReference.h"
#include "AliTracker.h"
#include "AliTrackPointArray.h"
#include "THnSparse.h"
#include "TParticle.h"
#include "TMatrix.h"
#include "TParticlePDG.h"
#include "TParticle.h"
#include "TTreeStream.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TTreeStream.h"

ClassImp(AliTrackComparison)

//________________________________________________________________________
AliTrackComparison::AliTrackComparison() 
  :TNamed()
  , fStep(1)
  , fLowBinDY(-10)
  , fUpBinDY(10)
  , fLowBinDZ(-10)
  , fUpBinDZ(10)
  , fLowBinDSnp(-0.5)
  , fUpBinDSnp(0.5)
  , fLowBinDTheta(-0.5)
  , fUpBinDTheta(0.5)
  , fLowBin1Pt(-3)
  , fUpBin1Pt(3)
  , fLowBin1PtLoss(-2)
  , fUpBin1PtLoss(2)
  , fNBinsDY(50)
  , fNBinsDZ(50)
  , fNBinsDSnp(50)
  , fNBinsDTheta(50)
  , fNBins1Pt(50)
  , fNBins1PtLoss(100)
  , fLayerID(-1)
  , fFillAll(kTRUE)
  , fNCombineBin(1)
{  
  //
  // Default constructor
  //
  for (Int_t i=0;i<6; i++) fResolHisto[i]=0;
}


//________________________________________________________________________
AliTrackComparison::AliTrackComparison(const Text_t *name, const Text_t *title)
  :TNamed(name,title)
  , fStep(1)
  , fLowBinDY(-10)
  , fUpBinDY(10)
  , fLowBinDZ(-10)
  , fUpBinDZ(10)
  , fLowBinDSnp(-0.5)
  , fUpBinDSnp(0.5)
  , fLowBinDTheta(-0.5)
  , fUpBinDTheta(0.5)
  , fLowBin1Pt(-3)
  , fUpBin1Pt(3)
  , fLowBin1PtLoss(-2)
  , fUpBin1PtLoss(2)
  , fNBinsDY(50)
  , fNBinsDZ(50)
  , fNBinsDSnp(50)
  , fNBinsDTheta(50)
  , fNBins1Pt(50)
  , fNBins1PtLoss(100)
  , fLayerID(-1)
  , fFillAll(kTRUE)
  , fNCombineBin(1)
{
  //
  // Non default cosntructor
  //
  for (Int_t i=0;i<6; i++) fResolHisto[i]=0;
}

//________________________________________________________________________
AliTrackComparison::AliTrackComparison(const AliTrackComparison& comp)
  :TNamed(comp)
  , fStep(comp.fStep)
  , fLowBinDY(comp.fLowBinDY)
  , fUpBinDY(comp.fUpBinDY)
  , fLowBinDZ(comp.fLowBinDZ)
  , fUpBinDZ(comp.fUpBinDZ)
  , fLowBinDSnp(comp.fLowBinDSnp)
  , fUpBinDSnp(comp.fUpBinDSnp)
  , fLowBinDTheta(comp.fLowBinDTheta)
  , fUpBinDTheta(comp.fUpBinDTheta)
  , fLowBin1Pt(comp.fLowBin1Pt)
  , fUpBin1Pt(comp.fUpBin1Pt)
  , fLowBin1PtLoss(comp.fLowBin1PtLoss)
  , fUpBin1PtLoss(comp.fUpBin1PtLoss)
  , fNBinsDY(comp.fNBinsDY)
  , fNBinsDZ(comp.fNBinsDZ)
  , fNBinsDSnp(comp.fNBinsDSnp)
  , fNBinsDTheta(comp.fNBinsDTheta)
  , fNBins1Pt(comp.fNBins1Pt)
  , fNBins1PtLoss(comp.fNBins1PtLoss)
  , fLayerID(comp.fLayerID)
  , fFillAll(comp.fFillAll)
  , fNCombineBin(comp.fNCombineBin)
{
  //
  // copy constructor
  //

  for (Int_t i=0;i<6; i++) fResolHisto[i]=comp.fResolHisto[i];
}

//________________________________________________________________________
AliTrackComparison& AliTrackComparison::operator=(const AliTrackComparison& comp)
{
  //
  //
  //
  if(this != &comp) {
    TNamed::operator=(comp);
  }
  return *this;
}

//________________________________________________________________________
void AliTrackComparison::Init(){
  //
  //Initilized the output histograms
  //
  MakeHistos();

}

//________________________________________________________________________
AliTrackComparison::~AliTrackComparison(){
  //
  //
  //
}
 
//________________________________________________________________________
void AliTrackComparison::Analyze() {
  //
  //
  //
}

//________________________________________________________________________
Int_t AliTrackComparison::AddTracks(AliTrackReference *ref0,  AliTrackReference *ref1, Double_t mass, Int_t charge)
{
  // Make track param out of track reference
  // Test propagation from ref0 to ref1
  // Fill the THnSparse with ref0 and ref1

  AliExternalTrackParam *param0 = 0;
  AliExternalTrackParam *param1 = 0;
   
  param0=MakeTrack(ref0,charge);
  param1=MakeTrack(ref1,charge);

  if (!param0 || !param1) return 0;
  
  Double_t tr1Pt = param0->GetSigned1Pt();

  if(!PropagateToPoint(param0,param1,mass)) return 0;

  FillHistos(param0,param1,tr1Pt);
  return 1;
}

//________________________________________________________________________
Int_t AliTrackComparison::AddTracks(AliExternalTrackParam *param0,  AliExternalTrackParam *param1, Double_t mass)
{
  //Test propagation from param0 to param1
  //Fill the THnSparse with param0 and param1
  //
  Double_t tr1Pt=param0->GetSigned1Pt();
  if( !PropagateToPoint(param0,param1,mass)) return 0;
  FillHistos(param0,param1,tr1Pt);
  return 1;
}

//________________________________________________________________________
Int_t AliTrackComparison::AddTracks(AliExternalTrackParam *param0,  AliTrackPoint *point1, Double_t mass, Double_t energy, Double_t *vxyz)
{
  //Test propagation from param0 to point1
  //This function is usually called in real data analysis when only the
  //position of the track point is known. In this case, the angles 
  //in the track param are not the angles of the track momentum, but position.
  //Only the first two and last two THnSparse are meaninglful and should be 
  //filled. Set this via SetFillAll(kFALSE)

  Double_t gPos[3]= {point1->GetX(),point1->GetY(),point1->GetZ()};

  Double_t pos[3], pxyz[3];
  for(Int_t i=0; i<3; i++) 
    pos[i] = gPos[i]-vxyz[i];
  Double_t R = TMath::Sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[3]);
  for(Int_t i=0; i<3; i++) pxyz[i]= energy*pos[i]/R;

  Double_t cv[21];
  for (Int_t i=0; i<21;i++) cv[i]=0;
  AliExternalTrackParam * param1 = new AliExternalTrackParam(gPos,pxyz,cv,param0->Charge());

  if(!param1) return 0;
  Double_t tr1Pt = param0->GetSigned1Pt();
  if(!PropagateToPoint(param0,param1,mass)) return 0;

  FillHistos(param0,param1,tr1Pt);
  return 1;
}

//________________________________________________________________________
Bool_t AliTrackComparison::PropagateToPoint(AliExternalTrackParam *param0, AliExternalTrackParam *param1,  Double_t mass)
{
  //
  //Extrapolate is performed 
  //
  Double_t radius = param1->GetX();
  param0->Rotate(param1->GetAlpha());
  AliTracker::PropagateTrackToBxByBz(param0, radius+fStep, mass, fStep, kFALSE,0.99,-1);
  Bool_t isOK = param0->PropagateTo(radius,AliTracker::GetBz());
  return isOK;
}

//________________________________________________________________________
AliExternalTrackParam * AliTrackComparison::MakeTrack(const AliTrackReference* ref, Int_t charge)
{
  //
  // Make track out of the track ref
  // the covariance matrix - equal 0 - starting from ideal MC position
  Double_t xyz[3]={ref->X(),ref->Y(),ref->Z()};
  Double_t pxyz[3]={ref->Px(),ref->Py(),ref->Pz()};
  Double_t cv[21];
  for (Int_t i=0; i<21;i++) cv[i]=0;
  AliExternalTrackParam * param = new AliExternalTrackParam(xyz,pxyz,cv,charge);
  return param;
}

//________________________________________________________________________
Long64_t AliTrackComparison::Merge(TCollection *const li) {
  //
  //Merge the comparison objects
  //   
  TIterator* iter = li->MakeIterator();
  AliTrackComparison* comp = 0;
  TString strName(GetName());
  while ((comp = (AliTrackComparison*)iter->Next())) {
    if (!comp->InheritsFrom(AliTrackComparison::Class())) {
      return -1;
    }
    if (strName.CompareTo(comp->GetName())!=0) return -1;
    // add histograms here...
    Add(comp);
  }
  return 0;  
}

//________________________________________________________________________
void  AliTrackComparison::Add(AliTrackComparison *const comp){
  //
  // Add THnSparse
  //
  if (!fResolHisto) return;
  for (Int_t i=0;i<6;i++){
    THnSparse * h0 = (THnSparse*)fResolHisto[i];
    THnSparse * h1 = (THnSparse*)comp->GetHnSparse(i);
    if (h0&&h1) h0->Add(h1);
  }
}


//________________________________________________________________________
void AliTrackComparison::MakeHistos()
{
  //
  //Called in Init() to initialize histograms
  // 
  Double_t xminTrack[5], xmaxTrack[5];
  Int_t   binsTrack[5];
  TString axisName[5];
  TString axisTitle[5];
  TString hisNames[6]={"DeltaY","DeltaZ","DeltaSnp","DeltaTheta","Delta1Pt","1PtLoss"};
  Double_t lowBins[6]={fLowBinDY, fLowBinDZ, fLowBinDSnp, fLowBinDTheta, fLowBin1Pt, fLowBin1PtLoss};
  Double_t upBins[6]={fUpBinDY, fUpBinDZ, fUpBinDSnp, fUpBinDTheta, fUpBin1Pt, fUpBin1PtLoss};
  Int_t nBins[6] = {fNBinsDY, fNBinsDZ, fNBinsDSnp, fNBinsDTheta, fNBins1Pt, fNBins1PtLoss};
  //
  axisName[0]   ="#Delta";
  axisTitle[0]  ="#Delta";
  //
  binsTrack[1] =40;
  xminTrack[1] =-1.0; xmaxTrack[1]=1.0;
  axisName[1]  ="tanTheta";
  axisTitle[1]  ="tan(#Theta)";
  //
  binsTrack[2] =180;
  xminTrack[2] =-TMath::Pi(); xmaxTrack[2]=TMath::Pi(); 
  axisName[2]  ="phi";
  axisTitle[2]  ="#phi";
  //
  binsTrack[3] =120;
  xminTrack[3] =-3.; xmaxTrack[3]=3.;   // 0.33 GeV cut 
  axisName[3]  ="1pt";
  axisTitle[3]  ="1/pt";
  //
  binsTrack[4] =22;
  xminTrack[4] =0; xmaxTrack[4]=22;  
  axisName[4]  ="LayerID";
  axisTitle[4] ="LayerID";

  for (Int_t ihis=0; ihis<6; ihis++){
    // modify ranges for bin 0
    //
    binsTrack[0]=nBins[ihis];
    xminTrack[0]=lowBins[ihis];
    xmaxTrack[0]=upBins[ihis];
    fResolHisto[ihis] = new THnSparseF(hisNames[ihis],hisNames[ihis],  5, binsTrack,xminTrack, xmaxTrack);
    for (Int_t ivar=0;ivar<5;ivar++){
      fResolHisto[ihis]->GetAxis(ivar)->SetName(axisName[ivar].Data());
      fResolHisto[ihis]->GetAxis(ivar)->SetTitle(axisTitle[ivar].Data());
    }
  }
}

//________________________________________________________________________
void AliTrackComparison::FillHistos(AliExternalTrackParam *param0, AliExternalTrackParam *param1, Double_t tr1Pt)
{
  //Fill the THnSparse. 
  //In case of not filling all, only the 4 out of 6 THnSparse are filed
  //
  Double_t dY     = param1->GetY()-param0->GetY();        //Residual in Y
  Double_t dZ     = param1->GetZ()-param0->GetZ();        //Residual in Z
  Double_t dSnp   = param1->GetSnp()-param0->GetSnp();    //Residual in sin(phi)
  Double_t dTheta = param1->Theta()-param0->Theta();      //Residual in Theta
  Double_t d1Pt   = param1->GetSigned1Pt()-param0->GetSigned1Pt(); //Residual in 1/pT
  Double_t d1PtLoss = param0->GetSigned1Pt()-tr1Pt;       //Corrected energy loss

  Double_t globalPos[3];
  param1->GetXYZ(globalPos);
  //printf("param1: Atan(y,x)=%5.3f | phi=%5.3f\n",TMath::ATan2(globalPos[1],globalPos[0]),param1->Phi());
  //printf("trkP=%5.3f | param0=%5.3f | param1=%5.3f\n",trP,param0->GetP(),param1->GetP());

  Double_t residual[6] = {dY,dZ,dSnp,dTheta,d1Pt,d1PtLoss};
  Double_t argu[6][5];
  for(Int_t j=0; j<6; j++)
    {
      if(!fFillAll && j<4 && j>1) continue;
      argu[j][0] = residual[j];
      argu[j][1] = param1->GetTgl();
      argu[j][2] = TMath::ATan2(globalPos[1],globalPos[0]);
      argu[j][3] = param1->GetSigned1Pt();
      argu[j][4] = fLayerID;
      fResolHisto[j]->Fill(argu[j]);
    } 
}

//________________________________________________________________________
void   AliTrackComparison::MakeDistortionMap(Double_t refX, Int_t type){
  //
  // make a distortion map out of the residual histogram
  // Results are written to the debug streamer - pcstream
  // Parameters:
  //   pcstream   - file to write the tree
  //   run        - run number
  //   refX       - track matching reference X
  //   type       - 0- y 1-z,2 -snp, 3-theta, 4-1/pt, 5-corrected 1/pT
  // THnSparse axes:
  // OBJ: TAxis     #Delta  #Delta
  // OBJ: TAxis     tanTheta        tan(#Theta)
  // OBJ: TAxis     phi     #phi
  // OBJ: TAxis     1/pt    1/pt
  // OBJ: TAxis     LayerID LayerID  
  // marian.ivanov@cern.ch
  
  //Double_t refX=365;   //temporary
  Int_t    run=0;   //temporary
  TString hname = Form("%s_%d",GetName(),type);
  THnSparse * his0= GetHnSparse(type);
  
  TString dsName;
  dsName=GetName();
  dsName+=Form("Debug_%d.root",type);
  printf(" Create debug streamer \n");
  dsName.ReplaceAll(" ","");
  TTreeSRedirector *pcstream = new TTreeSRedirector(dsName.Data());


  const Int_t kMinEntries=10;
  Double_t bz=AliTrackerBase::GetBz();
  Int_t idim[4]={0,1,2,3};
  //
  //
  //
  Int_t nbins3=his0->GetAxis(3)->GetNbins();
  Int_t first3=his0->GetAxis(3)->GetFirst();
  Int_t last3 =his0->GetAxis(3)->GetLast();
  //
  for (Int_t ibin3=first3; ibin3<last3; ibin3+=1){   // axis 3 - 1/pt
    his0->GetAxis(3)->SetRange(TMath::Max(ibin3-fNCombineBin,1),TMath::Min(ibin3+fNCombineBin,nbins3));
    Double_t      x3= his0->GetAxis(3)->GetBinCenter(ibin3);
    THnSparse * his3= his0->Projection(3,idim);         //projected histogram according selection 3
    //
    Int_t nbins2    = his3->GetAxis(2)->GetNbins();
    Int_t first2    = his3->GetAxis(2)->GetFirst();
    Int_t last2     = his3->GetAxis(2)->GetLast();
    //
    for (Int_t ibin2=first2; ibin2<last2; ibin2+=1){   // axis 2 - phi
      his3->GetAxis(2)->SetRange(TMath::Max(ibin2-fNCombineBin,1),TMath::Min(ibin2+fNCombineBin,nbins2));
      Double_t      x2= his3->GetAxis(2)->GetBinCenter(ibin2);
      THnSparse * his2= his3->Projection(2,idim);         //projected histogram according selection 2
      Int_t nbins1     = his2->GetAxis(1)->GetNbins();
      Int_t first1     = his2->GetAxis(1)->GetFirst();
      Int_t last1      = his2->GetAxis(1)->GetLast();
      for (Int_t ibin1=first1; ibin1<last1; ibin1++){   //axis 1 - tan(theta)
	//
	Double_t       x1= his2->GetAxis(1)->GetBinCenter(ibin1);
	Double_t binWidth= his2->GetAxis(1)->GetBinWidth(ibin1);

	his2->GetAxis(1)->SetRange(TMath::Max(ibin1-fNCombineBin,1),TMath::Min(ibin1+fNCombineBin,nbins1));
	if (TMath::Abs(x1)<(fNCombineBin+0.5)*binWidth){
	  if (x1<0) his2->GetAxis(1)->SetRange(TMath::Max(ibin1-fNCombineBin,1),TMath::Min(ibin1,nbins1));
	  if (x1>0) his2->GetAxis(1)->SetRange(TMath::Max(ibin1,1),TMath::Min(ibin1+fNCombineBin,nbins1));
	}
	TH1 * hisDelta = his2->Projection(0);
	//
	Double_t entries = hisDelta->GetEntries();
	Double_t mean=0, rms=0;
	if (entries>kMinEntries){
	  mean    = hisDelta->GetMean(); 
	  rms = hisDelta->GetRMS(); 
	}
	Double_t sector = 9.*x2/TMath::Pi();
	if (sector<0) sector+=18;
	Double_t dsec = sector-Int_t(sector)-0.5;
	Double_t quantiles[5]={0};
	Double_t mean75=0, rms75=0;
	Double_t prob[5]={0, 0.25,0.5,0.75,1};
	if (entries>kMinEntries){
	  hisDelta->GetQuantiles(5,quantiles,prob);
	  if(x3>0)hisDelta->GetXaxis()->SetRangeUser(quantiles[0], quantiles[3]);
	  else hisDelta->GetXaxis()->SetRangeUser(quantiles[1], quantiles[4]);
	  rms75=hisDelta->GetRMS();
	  mean75=hisDelta->GetMean();
	}
	
	Double_t pt = TMath::Abs(x3);
	Double_t z=refX*x1;
	(*pcstream)<<"distortion"<<
	  "run="<<run<<
	  "bz="<<bz<<
	  "theta="<<x1<<          // tan(theta)
	  "phi="<<x2<<            // phi
	  "z="<<z<<               // dummy z
	  "mpt="<<x3<<            // signed 1/pt
	  "1Pt="<<pt<<
	  "entries="<<entries<<   // entries
	  "mean="<<mean<<         // normal mean
	  //
	  "q25="<<quantiles[1]<<  // quantiles of distribution - importnat for asymetric distibutions
	  "q50="<<quantiles[2]<<  // quantiles of distribution - importnat for asymetric distibutions
	  "q75="<<quantiles[3]<<  // quantiles of distribution - importnat for asymetric distibutions
	  "mean75="<<mean75<<     // mean of the truncated distribution 0-75%
	  "rms75="<<rms75<<       // rms of the truncated distibution up to 0-75%
	  //
	  "rms="<<rms<<           // normal rms
	  "refX="<<refX<<         // track matching refernce plane
	  "type="<<type<<         // tpye of residuals
	  "sector="<<sector<<     // sector number according to TPC standard
	  "dsec="<<dsec<<
	  "\n";
	delete hisDelta;
	//printf("%f\t%f\t%f\t%f\t%f\n",x3,x2,x1, entries,mean);
      }
      delete his2;
    }
    delete his3;
  }

  printf(" Finished! Close debug streamer \n");
  delete pcstream;
}

/*
  .x ~/rootlogon.C
   //TFile f("outFile.root");
   //TList * list = f.Get("jthaeder_OutterTracking");
  TFile f("/u/alice-st/marr/EnergyCorrection/TrainAnalysis/trunk/marr_TrackingTest.root");
  TList * list = f.Get("marr_TrackingTest");

  AliTrackComparison * compTOF = list->FindObject("TPCOutToTOFIn");
  AliTrackComparison * compEMCAL = list->FindObject("TPCOutToEMCalInPion");
  AliTrackComparison * compEMCALEl = list->FindObject("TPCOutToEMCalInElec");
  AliTrackComparison * compHMPID = list->FindObject("TPCOutToHMPIDIn");
  Double_t refX=438;

  compEMCAL-> MakeDistortionMap(refX,0);
  compEMCAL-> MakeDistortionMap(refX,1);
  compEMCAL-> MakeDistortionMap(refX,4);
  compEMCALEl-> MakeDistortionMap(refX,0);
  compEMCALEl-> MakeDistortionMap(refX,1);
  compEMCALEl-> MakeDistortionMap(refX,4);
  
  compTOF->SetNCombineBin(3)
  compTOF->MakeDistortionMap(365,4);
  compTOF->MakeDistortionMap(365,0);
  
  TFile f("emcalDistortion.root");
  

  // ideal case - no mislaignment dy only due energy loss correction
  TPCOutToEMCalInElec_0->Draw("mean:mpt","entries>10");
  //
  TPCOutToEMCalInElec_4->Draw("mean/abs(mpt):mpt","entries>10");
  // delta 1/pt
  TPCOutToEMCalInPion_4->Draw("mean:mpt","entries>10");
  // (delta pt)/pt
  TPCOutToEMCalInPion_4->Draw("mean/abs(mpt):mpt","entries>10");  // pions  - bias +-1 %
  //
  TPCOutToEMCalInElec_4->Draw("mean/abs(mpt):mpt","entries>10");  // electrons  - bias +-20 %
  //




*/
