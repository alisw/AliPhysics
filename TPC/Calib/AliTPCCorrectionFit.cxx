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
  Responsible: marian.ivanov@cern.ch 
  Tools for fitting of the space point distortion parameters.
  Functionality  
  
 
1. Creation of the distortion maps from the residual histograms

2. Making fit trees

3. Space point distortion not directly observable. Instead a derived variables
   like DCA at vertex, local y distortion in the TPC-*TOF,TRD,ITS) matching
   in all 5 tracking parameters are obsereved.
   In the AliTPCcorrection fir code we calculate the derivative of given variables
   dO_{i}/dp_{i}

4. Global fit - later
   d0 = sum{ki*dO_{i}/dp_{i}}  - linear fitting of the amplitudes ki

 Some functions, for the moment function present in the AliTPCPreprocesorOffline, some will be 
 extracted from the old macros 


*/

#include "Riostream.h"
#include <fstream>
#include "TMap.h"
#include "TGraphErrors.h"
#include "AliExternalTrackParam.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "THnBase.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TPad.h"
#include "TH2D.h"
#include "TH3D.h"
#include "AliTPCROC.h"
#include "AliTPCCalROC.h"
#include "AliESDfriend.h"
#include "AliTPCcalibTime.h"
#include "AliSplineFit.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliTPCcalibBase.h"
#include "AliTPCcalibDB.h"
#include "AliTPCcalibDButil.h"
#include "AliRelAlignerKalman.h"
#include "AliTPCParamSR.h"
#include "AliTPCcalibTimeGain.h"
#include "AliTPCcalibGainMult.h"
#include "AliSplineFit.h"
#include "AliTPCComposedCorrection.h"
#include "AliTPCExBTwist.h"
#include "AliTPCCalibGlobalMisalignment.h"
#include "TStatToolkit.h"
#include "TChain.h"
#include "TCut.h"
#include "AliTrackerBase.h"
#include "AliTPCCorrectionFit.h"
#include "AliTPCLaserTrack.h"
#include "TDatabasePDG.h"
#include "AliTPCcalibAlign.h"
#include "AliLog.h"
#include "AliRieman.h"

ClassImp(AliTPCCorrectionFit)

AliTPCCorrectionFit::AliTPCCorrectionFit():
  TNamed("TPCCorrectionFit","TPCCorrectionFit")
{
  //
  // default constructor
  //
}

AliTPCCorrectionFit::~AliTPCCorrectionFit() {
  //
  // Destructor
  //
}


Double_t AliTPCCorrectionFit::EvalAt(Double_t phi, Double_t refX, Double_t theta, Int_t corr, Int_t ptype){
  //
  // Evalution at point using the lienar approximation
  //
  Double_t sector = 9*phi/TMath::Pi();
  if (sector<0) sector+=18;
  Double_t y85=AliTPCCorrection::GetCorrSector(sector,85,theta,1,corr);
  Double_t y245=AliTPCCorrection::GetCorrSector(sector,245,theta,1,corr);
  if (ptype==0) return y85+(y245-y85)*(refX-85.)/(245.-85.);
  if (ptype==2) return (y245-y85)/(245.-85.);
  return 0;
}




Double_t AliTPCCorrectionFit::EvalAtPar(Double_t phi0, Double_t snp, Double_t refX, Double_t theta, Int_t corr, Int_t ptype, Int_t nsteps){
  //
  // Fit the distortion along the line with the parabolic model
  // We assume that the track are primaries  - where the vertex is at (0,0,0)
  //
  // Parameters:
  //  phi0  -  phi at the entrance of the TPC
  //  snp   -  local inclination angle at the entrance of the TPC
  //  refX  -  ref X where the distortion is evanluated
  //  theta -  inclination angle
  //  corr  -  internal number of the correction and distortion 
  //  ptype -  0 - local y distortion
  //           //1 - local z distortion
  //           2 - local derivative distortion
  //           //3
  //           4 - distortion in the curvature
  //  nsteps - number of fit points
  //
  // return value -  distortion at point refX with type ptype
  //
  static TLinearFitter fitter(3,"pol2"); 
  fitter.ClearPoints();
  if (nsteps<3) nsteps=3;
  Double_t deltaX=(245-85)/(nsteps);
  for (Int_t istep=0; istep<(nsteps+1); istep++){
    //
    Double_t localX =85.+deltaX*istep;
    Double_t localPhi=phi0+deltaX*snp*istep;
    Double_t sector = 9*localPhi/TMath::Pi();
    if (sector<0) sector+=18;
    Double_t dy=AliTPCCorrection::GetCorrSector(sector,localX,theta,1,corr);
    Double_t dlocalX=AliTPCCorrection::GetCorrSector(sector,localX,theta,0,corr);
    Double_t x[1]={localX-dlocalX};
    fitter.AddPoint(x,dy);
  }
  fitter.Eval();
  Double_t par[3];
  par[0]=fitter.GetParameter(0);
  par[1]=fitter.GetParameter(1);
  par[2]=fitter.GetParameter(2);

  if (ptype==0) return par[0]+par[1]*refX+par[2]*refX*refX;
  if (ptype==2) return par[1]+2*par[2]*refX;
  if (ptype==4) return par[2];
  return 0;
}


Double_t AliTPCCorrectionFit::EvalAtHelix(Double_t phi0, Double_t snp, Double_t refX, Double_t theta, Int_t corr, Int_t ptype, Int_t nsteps){
  //
  // Fit the distortion along the line with the helix model
  // FIXME - original trajectory to be changed - AliHelix to be used
  // We assume that the track are primaries  - where the vertex is at (0,0,0)
  //
  // Parameters:
  //  phi0  -  phi at the entrance of the TPC
  //  snp   -  local inclination angle at the entrance of the TPC
  //  refX  -  ref X where the distortion is evanluated
  //  theta -  inclination angle
  //  corr  -  internal number of the correction and distortion 
  //  ptype -  0 - local y distortion
  //           //1 - local z distortion
  //           2 - local derivative distortion
  //           //3
  //           4 - distortion in the curvature
  //  nsteps - number of fit points
  //
  // return value -  distortion at point refX with type ptype
  //
  if (nsteps<3) nsteps=3;
  Double_t deltaX=(245-85)/(nsteps);
  AliRieman rieman(nsteps);

  for (Int_t istep=0; istep<(nsteps+1); istep++){
    //
    Double_t localX =85.+deltaX*istep;
    Double_t localPhi=phi0+deltaX*snp*istep;
    Double_t sector = 9*localPhi/TMath::Pi();
    if (sector<0) sector+=18;
    Double_t dy=AliTPCCorrection::GetCorrSector(sector,localX,theta,1,corr);
    Double_t dlocalX=AliTPCCorrection::GetCorrSector(sector,localX,theta,0,corr);
    Double_t x[1]={localX-dlocalX};
    Double_t z=theta*x[0];
    rieman.AddPoint(x[0],dy,z,0.1,0.1);
  }
  rieman.Update();
  //
 
  if (ptype==0) return rieman.GetYat(refX);
  if (ptype==2) return rieman.GetDYat(refX);
  if (ptype==4) return rieman.GetC();
  return 0;
}




void AliTPCCorrectionFit::CreateAlignMaps(Double_t bz, Int_t run){
  //
  // Create cluster distortion map
  //
  TFile *falign = TFile::Open("CalibObjects.root");
  TObjArray * arrayAlign = (TObjArray *)falign->Get("TPCAlign");
  if (!arrayAlign) {
    AliWarningGeneral("AliTPCCorrectionFit::CreateAlignMaps","Alignment was not included in the calibration task");
    return;
  }
  AliTPCcalibAlign * align =  (AliTPCcalibAlign *)arrayAlign->FindObject("alignTPC");
  if (!align) {
      AliWarningGeneral("AliTPCCorrectionFit::CreateAlignMaps","Alignment was not included in the calibration task");
    return;
  }
  TTreeSRedirector * pcstream = new TTreeSRedirector("TPCAlign.root");

  THnBase * hdY = (THnBase*)align->GetClusterDelta(0);
  //THnBase * hdZ = (THnBase*)align->GetClusterDelta(1);
  AliTPCCorrectionFit::MakeClusterDistortionMap(hdY,pcstream,0, bz);
  //  AliTPCCorrectionFit::MakeClusterDistortionMap(hdZ,pcstream,1, bz);

  const char * hname[5]={"dy","dz","dsnp","dtheta","d1pt"};
  for (Int_t ihis=0; ihis<4; ihis++){
    THnSparse * hisAlign =align->GetTrackletDelta(ihis);
    AliTPCCorrectionFit::MakeDistortionMapSector(hisAlign, pcstream, hname[ihis], run, ihis,bz);
  }
  delete pcstream;
  delete falign;
}





void  AliTPCCorrectionFit::MakeClusterDistortionMap(THnBase * hisInput,TTreeSRedirector *pcstream , Int_t ptype, Int_t dtype){
  //
  // Make cluster residual map from the n-dimensional histogram
  // hisInput supposed to have given format:
  //     4 Dim:
  //      delta, 
  //      sector 
  //      localX
  //      kZ
  // Vertex position assumed to be at (0,0,0)          
  //
  //TTreeSRedirector *pcstream=new TTreeSRedirector(sname);
  //
  Int_t nbins1=hisInput->GetAxis(1)->GetNbins();
  Int_t nbins2=hisInput->GetAxis(2)->GetNbins();
  Int_t nbins3=hisInput->GetAxis(3)->GetNbins();
  TF1 *fgaus=0;
  TH3F * hisResMap3D = 
    new TH3F("his3D","his3D",
	     nbins1,hisInput->GetAxis(1)->GetXmin(), hisInput->GetAxis(1)->GetXmax(),
	     nbins2,hisInput->GetAxis(2)->GetXmin(), hisInput->GetAxis(2)->GetXmax(),
	     nbins3,hisInput->GetAxis(3)->GetXmin(), hisInput->GetAxis(3)->GetXmax());
  hisResMap3D->GetXaxis()->SetTitle("sector");
  hisResMap3D->GetYaxis()->SetTitle("localX");
  hisResMap3D->GetZaxis()->SetTitle("kZ");

  TH2F * hisResMap2D[4] ={0,0,0,0};
  for (Int_t i=0; i<4; i++){
    hisResMap2D[i]=
      new TH2F(Form("his2D_0%d",i),Form("his2D_0%d",i),
	       nbins1,hisInput->GetAxis(1)->GetXmin(), hisInput->GetAxis(1)->GetXmax(),
	       nbins2,hisInput->GetAxis(2)->GetXmin(), hisInput->GetAxis(2)->GetXmax());
    hisResMap2D[i]->GetXaxis()->SetTitle("sector");
    hisResMap2D[i]->GetYaxis()->SetTitle("localX");
  }
  //
  //
  //
  TF1 * f1= 0;
  Int_t axis0[4]={0,1,2,3};
  Int_t axis1[4]={0,1,2,3};
  Int_t counter=0;
  for (Int_t ibin1=1; ibin1<nbins1; ibin1+=1){
    // phi- sector  range
    hisInput->GetAxis(1)->SetRange(ibin1-1,ibin1+1);
    THnBase *his1=(THnBase *)hisInput->ProjectionND(4,axis0); 
    Double_t sector=hisInput->GetAxis(1)->GetBinCenter(ibin1);
    //
    for (Int_t ibin2=1; ibin2<nbins2; ibin2+=1){
      // local x range
      // kz fits
      his1->GetAxis(2)->SetRange(ibin2-1,ibin2+1);
      THnBase *his2=(THnBase *)his1->ProjectionND(4,axis1); 
      Double_t localX=hisInput->GetAxis(2)->GetBinCenter(ibin2);
      //      
      //A side
      his2->GetAxis(3)->SetRangeUser(0.01,0.3);
      TH1 * hisA = his2->Projection(0);
      Double_t meanA= hisA->GetMean();
      Double_t rmsA= hisA->GetRMS();
      Double_t entriesA= hisA->GetEntries();
      delete hisA;
      //C side
      his2->GetAxis(3)->SetRangeUser(0.01,0.3);
      TH1 * hisC = his2->Projection(0);
      Double_t meanC= hisC->GetMean();
      Double_t rmsC= hisC->GetRMS();
      Double_t entriesC= hisC->GetEntries();
      delete hisC;
      his2->GetAxis(3)->SetRangeUser(-1.2,1.2);      
      TH2 * hisAC       = his2->Projection(0,3);
      TProfile *profAC  = hisAC->ProfileX(); 
      delete hisAC;
      profAC->Fit("pol1","QNR","QNR",0.05,1);
      if (!f1) f1=(TF1*)gROOT->FindObject("pol1");
      Double_t offsetA=f1->GetParameter(0);
      Double_t slopeA=f1->GetParameter(1);
      Double_t offsetAE=f1->GetParError(0);
      Double_t slopeAE=f1->GetParError(1); 
      Double_t chi2A=f1->GetChisquare()/f1->GetNumberFreeParameters();
      profAC->Fit("pol1","QNR","QNR",-1.1,-0.1);
      f1=(TF1*)gROOT->FindObject("pol1");
      Double_t offsetC=f1->GetParameter(0);
      Double_t slopeC=f1->GetParameter(1); 
      Double_t offsetCE=f1->GetParError(0);
      Double_t slopeCE=f1->GetParError(1); 
      Double_t chi2C=f1->GetChisquare()/f1->GetNumberFreeParameters();
      if (counter%50==0) printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\n", sector,localX, entriesA+entriesC, slopeA,slopeC, chi2A, chi2C);
      counter++;
      (*pcstream)<<"deltaFit"<<
	"sector="<<sector<<
	"localX="<<localX<<
	"meanA="<<meanA<<
	"rmsA="<<rmsA<<
	"entriesA="<<entriesA<<
	"meanC="<<meanC<<
	"rmsC="<<rmsC<<
	"entriesC="<<entriesC<<
	"offsetA="<<offsetA<<
	"slopeA="<<slopeA<<
	"offsetAE="<<offsetAE<<
	"slopeAE="<<slopeAE<<
	"chi2A="<<chi2A<<
	"offsetC="<<offsetC<<
	"slopeC="<<slopeC<<
	"offsetCE="<<offsetCE<<
	"slopeCE="<<slopeCE<<
	"chi2C="<<chi2C<<
	"\n";
      //
      hisResMap2D[0]->SetBinContent(ibin1,ibin2, offsetA);
      hisResMap2D[1]->SetBinContent(ibin1,ibin2, slopeA);
      hisResMap2D[2]->SetBinContent(ibin1,ibin2, offsetC);
      hisResMap2D[3]->SetBinContent(ibin1,ibin2, slopeC);
      
      for (Int_t ibin3=1; ibin3<nbins3; ibin3++){
	Double_t kZ=hisInput->GetAxis(3)->GetBinCenter(ibin3);
	if (TMath::Abs(kZ)<0.05) continue;  // crossing 
	his2->GetAxis(3)->SetRange(ibin3,ibin3);
	if (TMath::Abs(kZ)>0.15){
	  his2->GetAxis(3)->SetRange(ibin3,ibin3);
	}
	TH1 * his = his2->Projection(0);
	Double_t mean= his->GetMean();
	Double_t rms= his->GetRMS();
	Double_t entries= his->GetEntries();
	//printf("%f\t%f\t%f\t%f\t%f\t%f\n", sector,localX,kZ, entries, mean,rms);
	hisResMap3D->SetBinContent(ibin1,ibin2,ibin3, mean);
	Double_t phi=TMath::Pi()*sector/9;
	if (phi>TMath::Pi()) phi+=TMath::Pi();
	Double_t meanG=0;
	Double_t rmsG=0;
	if (entries>50){
	  if (!fgaus) {	    
	    his->Fit("gaus","Q","goff");
	    fgaus= (TF1*)((his->GetListOfFunctions()->FindObject("gaus"))->Clone());
	  }
	  if (fgaus) {
	    his->Fit(fgaus,"Q","goff");
	    meanG=fgaus->GetParameter(1);
	    rmsG=fgaus->GetParameter(2);
	  }
	}
	Double_t dsec=sector-Int_t(sector)-0.5;
	Double_t snp=dsec*TMath::Pi()/9.;
	(*pcstream)<<"delta"<<
	  "ptype="<<ptype<<
	  "dtype="<<dtype<<
	  "sector="<<sector<<
	  "dsec="<<dsec<<
	  "snp="<<snp<<
	  "phi="<<phi<<
	  "localX="<<localX<<
	  "kZ="<<kZ<<
	  "theta="<<kZ<<
	  "mean="<<mean<<
	  "rms="<<rms<<
	  "meanG="<<meanG<<
	  "rmsG="<<rmsG<<
	  "entries="<<entries<<
	  "meanA="<<meanA<<
	  "rmsA="<<rmsA<<
	  "entriesA="<<entriesA<<
	  "meanC="<<meanC<<
	  "rmsC="<<rmsC<<
	  "entriesC="<<entriesC<<
	  "offsetA="<<offsetA<<
	  "slopeA="<<slopeA<<
	  "chi2A="<<chi2A<<
	  "offsetC="<<offsetC<<
	  "slopeC="<<slopeC<<
	  "chi2C="<<chi2C<<
	  "\n";
	delete his;
      }
      delete his2;
    }
    delete his1;
  }
  hisResMap3D->Write();
  hisResMap2D[0]->Write();
  hisResMap2D[1]->Write();
  hisResMap2D[2]->Write();
  hisResMap2D[3]->Write();
  //  delete pcstream;
}



void   AliTPCCorrectionFit::MakeDistortionMap(THnSparse * his0, TTreeSRedirector * const pcstream, const char* hname, Int_t run, Float_t refX, Int_t type, Int_t integ, Double_t bz){
  //
  // make a distortion map out ou fthe residual histogram
  // Results are written to the debug streamer - pcstream
  // Parameters:
  //   his0       - input (4D) residual histogram
  //   pcstream   - file to write the tree
  //   run        - run number
  //   refX       - track matching reference X
  //   type       - 0- y 1-z,2 -snp, 3-theta, 4=1/pt
  // THnSparse axes:
  // OBJ: TAxis     #Delta  #Delta
  // OBJ: TAxis     tanTheta        tan(#Theta)
  // OBJ: TAxis     phi     #phi
  // OBJ: TAxis     snp     snp

  // marian.ivanov@cern.ch
  const Int_t kMinEntries=10;
  Int_t idim[4]={0,1,2,3};
  //
  //
  //
  Int_t nbins3=his0->GetAxis(3)->GetNbins();
  Int_t first3=his0->GetAxis(3)->GetFirst();
  Int_t last3 =his0->GetAxis(3)->GetLast();
  //
  for (Int_t ibin3=first3; ibin3<last3; ibin3+=1){   // axis 3 - local angle
    his0->GetAxis(3)->SetRange(TMath::Max(ibin3-integ,1),TMath::Min(ibin3+integ,nbins3));
    Double_t      x3= his0->GetAxis(3)->GetBinCenter(ibin3);
    THnSparse * his3= his0->Projection(3,idim);         //projected histogram according selection 3
    //
    Int_t nbins2    = his3->GetAxis(2)->GetNbins();
    Int_t first2    = his3->GetAxis(2)->GetFirst();
    Int_t last2     = his3->GetAxis(2)->GetLast();
    //
    for (Int_t ibin2=first2; ibin2<last2; ibin2+=1){   // axis 2 - phi
      his3->GetAxis(2)->SetRange(TMath::Max(ibin2-integ,1),TMath::Min(ibin2+integ,nbins2));
      Double_t      x2= his3->GetAxis(2)->GetBinCenter(ibin2);
      THnSparse * his2= his3->Projection(2,idim);         //projected histogram according selection 2
      Int_t nbins1     = his2->GetAxis(1)->GetNbins();
      Int_t first1     = his2->GetAxis(1)->GetFirst();
      Int_t last1      = his2->GetAxis(1)->GetLast();
      for (Int_t ibin1=first1; ibin1<last1; ibin1++){   //axis 1 - theta
	//
	Double_t       x1= his2->GetAxis(1)->GetBinCenter(ibin1);
	his2->GetAxis(1)->SetRange(TMath::Max(ibin1-1,1),TMath::Min(ibin1+1,nbins1));
	if (TMath::Abs(x1)<0.1){
	  if (x1<0) his2->GetAxis(1)->SetRange(TMath::Max(ibin1-1,1),TMath::Min(ibin1,nbins1));
	  if (x1>0) his2->GetAxis(1)->SetRange(TMath::Max(ibin1,1),TMath::Min(ibin1+1,nbins1));
	}
	if (TMath::Abs(x1)<0.06){
	  his2->GetAxis(1)->SetRange(TMath::Max(ibin1,1),TMath::Min(ibin1,nbins1));
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
	Double_t z=refX*x1;
	(*pcstream)<<hname<<
	  "run="<<run<<
	  "bz="<<bz<<
	  "theta="<<x1<<
	  "phi="<<x2<<
	  "z="<<z<<            // dummy z
	  "snp="<<x3<<
	  "entries="<<entries<<
	  "mean="<<mean<<
	  "rms="<<rms<<
	  "refX="<<refX<<   // track matching refernce plane
	  "type="<<type<<   //
	  "sector="<<sector<<
	  "dsec="<<dsec<<
	  "\n";
	delete hisDelta;
	//printf("%f\t%f\t%f\t%f\t%f\n",x3,x2,x1, entries,mean);
      }
      delete his2;
    }
    delete his3;
  }
}



void   AliTPCCorrectionFit::MakeDistortionMapCosmic(THnSparse * hisInput, TTreeSRedirector * const pcstream, const char* hname, Int_t run, Float_t refX, Int_t type){
  //
  // make a distortion map out ou fthe residual histogram
  // Results are written to the debug streamer - pcstream
  // Parameters:
  //   his0       - input (4D) residual histogram
  //   pcstream   - file to write the tree
  //   run        - run number
  //   refX       - track matching reference X
  //   type       - 0- y 1-z,2 -snp, 3-theta, 4=1/pt
  // marian.ivanov@cern.ch
  //
  //  Histo axeses
  //   Collection name='TObjArray', class='TObjArray', size=16
  //  0. OBJ: TAxis     #Delta  #Delta
  //  1. OBJ: TAxis     N_{cl}  N_{cl}
  //  2. OBJ: TAxis     dca_{r} (cm)    dca_{r} (cm)
  //  3. OBJ: TAxis     z (cm)  z (cm)
  //  4. OBJ: TAxis     sin(#phi)       sin(#phi)
  //  5. OBJ: TAxis     tan(#theta)     tan(#theta)
  //  6. OBJ: TAxis     1/pt (1/GeV)    1/pt (1/GeV)
  //  7. OBJ: TAxis     pt (GeV)        pt (GeV)
  //  8. OBJ: TAxis     alpha   alpha
  const Int_t kMinEntries=10;
  //
  //  1. make default selections
  //
  TH1 * hisDelta=0;
  Int_t idim0[4]={0 , 5, 8,  3};   // delta, theta, alpha, z
  hisInput->GetAxis(1)->SetRangeUser(110,190);   //long tracks
  hisInput->GetAxis(2)->SetRangeUser(-10,35);    //tracks close to beam pipe
  hisInput->GetAxis(4)->SetRangeUser(-0.3,0.3); //small snp at TPC entrance
  hisInput->GetAxis(7)->SetRangeUser(3,100); //"high pt tracks"
  hisDelta= hisInput->Projection(0);
  hisInput->GetAxis(0)->SetRangeUser(-6.*hisDelta->GetRMS(), +6.*hisDelta->GetRMS());
  delete hisDelta;
  THnSparse *his0=  hisInput->Projection(4,idim0);
  //
  // 2. Get mean in diferent bins
  //
  Int_t nbins1=his0->GetAxis(1)->GetNbins();
  Int_t first1=his0->GetAxis(1)->GetFirst();
  Int_t last1 =his0->GetAxis(1)->GetLast();
  //
  Double_t bz=AliTrackerBase::GetBz();
  Int_t idim[4]={0,1, 2,  3};  // delta, theta,alpha,z
  //
  for (Int_t ibin1=first1; ibin1<=last1; ibin1++){   //axis 1 - theta
    //
    Double_t       x1= his0->GetAxis(1)->GetBinCenter(ibin1);  
    his0->GetAxis(1)->SetRange(TMath::Max(ibin1-1,1),TMath::Min(ibin1+1,nbins1));
    //
    THnSparse * his1 = his0->Projection(4,idim);  // projected histogram according range1
    Int_t nbins3     = his1->GetAxis(3)->GetNbins();
    Int_t first3     = his1->GetAxis(3)->GetFirst();
    Int_t last3      = his1->GetAxis(3)->GetLast();
    //
    for (Int_t ibin3=first3-1; ibin3<=last3; ibin3+=1){   // axis 3 - z at "vertex"
      his1->GetAxis(3)->SetRange(TMath::Max(ibin3-1,1),TMath::Min(ibin3+1,nbins3));
      Double_t      x3= his1->GetAxis(3)->GetBinCenter(ibin3);
      if (ibin3<first3) {
	his1->GetAxis(3)->SetRangeUser(-1,1);
	x3=0;
      }
      THnSparse * his3= his1->Projection(4,idim);         //projected histogram according selection 3
      Int_t nbins2    = his3->GetAxis(2)->GetNbins();
      Int_t first2    = his3->GetAxis(2)->GetFirst();
      Int_t last2     = his3->GetAxis(2)->GetLast();
      //
      for (Int_t ibin2=first2; ibin2<=last2; ibin2+=1){
	his3->GetAxis(2)->SetRange(TMath::Max(ibin2-1,1),TMath::Min(ibin2+1,nbins2));
	Double_t x2= his3->GetAxis(2)->GetBinCenter(ibin2);
	hisDelta = his3->Projection(0);
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
	Double_t snp=0;  // dummy snp - equal 0
	(*pcstream)<<hname<<
	  "run="<<run<<
	  "bz="<<bz<<            // magnetic field
	  "theta="<<x1<<         // theta
	  "phi="<<x2<<           // phi (alpha)
	  "z="<<x3<<             // z at "vertex"
	  "snp="<<snp<<          // dummy snp
	  "entries="<<entries<<  // entries in bin
	  "mean="<<mean<<        // mean
	  "rms="<<rms<<
	  "refX="<<refX<<        // track matching refernce plane
	  "type="<<type<<        // parameter type
	  "sector="<<sector<<    // sector
	  "dsec="<<dsec<<        // dummy delta sector
	  "\n";
	delete hisDelta;
	printf("%f\t%f\t%f\t%f\t%f\n",x1,x3,x2, entries,mean);
      }
      delete his3;
    }
    delete his1;
  }
  delete his0;
}



void   AliTPCCorrectionFit::MakeDistortionMapSector(THnSparse * hisInput, TTreeSRedirector * const pcstream, const char* hname, Int_t run, Int_t type, Double_t bz){
  //
  // make a distortion map out of the residual histogram
  // Results are written to the debug streamer - pcstream
  // Parameters:
  //   his0       - input (4D) residual histogram
  //   pcstream   - file to write the tree
  //   run        - run number
  //   type       - 0- y 1-z,2 -snp, 3-theta
  //   bz         - magnetic field
  // marian.ivanov@cern.ch

  //Collection name='TObjArray', class='TObjArray', size=16
  //0  OBJ: TAxis     delta   delta
  //1  OBJ: TAxis     phi     phi
  //2  OBJ: TAxis     localX  localX
  //3  OBJ: TAxis     kY      kY
  //4  OBJ: TAxis     kZ      kZ
  //5  OBJ: TAxis     is1     is1
  //6  OBJ: TAxis     is0     is0
  //7. OBJ: TAxis     z       z
  //8. OBJ: TAxis     IsPrimary       IsPrimary

  const Int_t kMinEntries=10;
  THnSparse * hisSector0=0;
  TH1 * htemp=0;    // histogram to calculate mean value of parameter
  //  Double_t bz=AliTrackerBase::GetBz();

  //
  // Loop over pair of sector:
  // isPrim         - 8  ==> 8
  // isec0          - 6  ==> 7
  //   isec1        - 5  ==> 6
  //     refX       - 2  ==> 5
  //
  //     phi        - 1  ==> 4
  //       z        - 7  ==> 3
  //         snp    - 3  ==> 2
  //           theta- 4  ==> 1
  //                  0  ==> 0;           
  for (Int_t isec0=0; isec0<72; isec0++){
    Int_t index0[9]={0, 4, 3, 7, 1, 2, 5, 6,8}; //regroup indeces
    //
    //hisInput->GetAxis(8)->SetRangeUser(-0.1,0.4);  // select secondaries only ? - get out later ?
    hisInput->GetAxis(6)->SetRangeUser(isec0-0.1,isec0+0.1);
    hisSector0=hisInput->Projection(7,index0);
    //
    //
    for (Int_t isec1=isec0+1; isec1<72; isec1++){    
      //if (isec1!=isec0+36) continue;
      if ( TMath::Abs((isec0%18)-(isec1%18))>1.5 && TMath::Abs((isec0%18)-(isec1%18))<16.5) continue;
      printf("Sectors %d\t%d\n",isec1,isec0);
      hisSector0->GetAxis(6)->SetRangeUser(isec1-0.1,isec1+0.1);      
      TH1 * hisX=hisSector0->Projection(5);
      Double_t refX= hisX->GetMean();
      delete hisX;
      TH1 *hisDelta=hisSector0->Projection(0);
      Double_t dmean = hisDelta->GetMean();
      Double_t drms = hisDelta->GetRMS();
      hisSector0->GetAxis(0)->SetRangeUser(dmean-5.*drms, dmean+5.*drms);
      delete hisDelta;
      //
      //  1. make default selections
      //
      Int_t idim0[5]={0 , 1, 2, 3, 4}; // {delta, theta, snp, z, phi }
      THnBase *hisSector1=  hisSector0->ProjectionND(5,idim0);
      //
      // 2. Get mean in diferent bins
      //
      Int_t idim[5]={0, 1, 2,  3, 4};  // {delta, theta-1,snp-2 ,z-3, phi-4}
      //
      //      Int_t nbinsPhi=hisSector1->GetAxis(4)->GetNbins();
      Int_t firstPhi=hisSector1->GetAxis(4)->GetFirst();
      Int_t lastPhi =hisSector1->GetAxis(4)->GetLast();
      //
      for (Int_t ibinPhi=firstPhi; ibinPhi<=lastPhi; ibinPhi+=2){   //axis 4 - phi
	//
	// Phi loop
	//
	Double_t       xPhi= hisSector1->GetAxis(4)->GetBinCenter(ibinPhi);         
	Double_t psec    = (9*xPhi/TMath::Pi());
	if (psec<0) psec+=18;
	Bool_t isOK0=kFALSE;
	Bool_t isOK1=kFALSE;
	if (TMath::Abs(psec-isec0%18-0.5)<1. || TMath::Abs(psec-isec0%18-17.5)<1.)  isOK0=kTRUE;
	if (TMath::Abs(psec-isec1%18-0.5)<1. || TMath::Abs(psec-isec1%18-17.5)<1.)  isOK1=kTRUE;
	if (!isOK0) continue;
	if (!isOK1) continue;
	//
	hisSector1->GetAxis(4)->SetRange(TMath::Max(ibinPhi-2,firstPhi),TMath::Min(ibinPhi+2,lastPhi));
	if (isec1!=isec0+36) {
	  hisSector1->GetAxis(4)->SetRange(TMath::Max(ibinPhi-3,firstPhi),TMath::Min(ibinPhi+3,lastPhi));
	}
	//
	htemp = hisSector1->Projection(4);
	xPhi=htemp->GetMean();
	delete htemp;
	THnBase * hisPhi = hisSector1->ProjectionND(4,idim);
	//Int_t nbinsZ     = hisPhi->GetAxis(3)->GetNbins();
	Int_t firstZ     = hisPhi->GetAxis(3)->GetFirst();
	Int_t lastZ      = hisPhi->GetAxis(3)->GetLast();
	//
	for (Int_t ibinZ=firstZ; ibinZ<=lastZ; ibinZ+=2){   // axis 3 - z
	  //
	  // Z loop
	  //
	  hisPhi->GetAxis(3)->SetRange(TMath::Max(ibinZ,firstZ),TMath::Min(ibinZ,lastZ));
	  if (isec1!=isec0+36) {
	    hisPhi->GetAxis(3)->SetRange(TMath::Max(ibinZ-1,firstZ),TMath::Min(ibinZ-1,lastZ));	    
	  }
	  htemp = hisPhi->Projection(3);
	  Double_t      xZ= htemp->GetMean();
	  delete htemp;
	  THnBase * hisZ= hisPhi->ProjectionND(3,idim);         
	  //projected histogram according selection 3 -z
	  //
	  //
	  //Int_t nbinsSnp    = hisZ->GetAxis(2)->GetNbins();
	  Int_t firstSnp    = hisZ->GetAxis(2)->GetFirst();
	  Int_t lastSnp     = hisZ->GetAxis(2)->GetLast();
	  for (Int_t ibinSnp=firstSnp; ibinSnp<=lastSnp; ibinSnp+=2){   // axis 2 - snp
	    //
	    // Snp loop
	    //
	    hisZ->GetAxis(2)->SetRange(TMath::Max(ibinSnp-1,firstSnp),TMath::Min(ibinSnp+1,lastSnp));
	    if (isec1!=isec0+36) {
	      hisZ->GetAxis(2)->SetRange(TMath::Max(ibinSnp-2,firstSnp),TMath::Min(ibinSnp+2,lastSnp));
	    }
	    htemp = hisZ->Projection(2);
	    Double_t      xSnp= htemp->GetMean();
	    delete htemp;
	    THnBase * hisSnp= hisZ->ProjectionND(2,idim);         
	    //projected histogram according selection 2 - snp
	    
	    //Int_t nbinsTheta    = hisSnp->GetAxis(1)->GetNbins();
	    Int_t firstTheta    = hisSnp->GetAxis(1)->GetFirst();
	    Int_t lastTheta     = hisSnp->GetAxis(1)->GetLast();
	    //
	    for (Int_t ibinTheta=firstTheta; ibinTheta<=lastTheta; ibinTheta+=2){  // axis1 theta
	      
	      
	      hisSnp->GetAxis(1)->SetRange(TMath::Max(ibinTheta-2,firstTheta),TMath::Min(ibinTheta+2,lastTheta));
	      if (isec1!=isec0+36) {
		 hisSnp->GetAxis(1)->SetRange(TMath::Max(ibinTheta-3,firstTheta),TMath::Min(ibinTheta+3,lastTheta));		 
	      }
	      htemp = hisSnp->Projection(1);	      
	      Double_t xTheta=htemp->GetMean();
	      delete htemp;
	      hisDelta = hisSnp->Projection(0);
	      //
	      Double_t entries = hisDelta->GetEntries();
	      Double_t mean=0, rms=0;
	      if (entries>kMinEntries){
		mean    = hisDelta->GetMean(); 
		rms = hisDelta->GetRMS(); 
	      }
	      Double_t sector = 9.*xPhi/TMath::Pi();
	      if (sector<0) sector+=18;
	      Double_t dsec = sector-Int_t(sector)-0.5;
	      Int_t dtype=1;  // TPC alignment type
	      (*pcstream)<<hname<<
		"run="<<run<<
		"bz="<<bz<<             // magnetic field
		"ptype="<<type<<         // parameter type
		"dtype="<<dtype<<         // parameter type
		"isec0="<<isec0<<       // sector 0 
		"isec1="<<isec1<<       // sector 1		
		"sector="<<sector<<     // sector as float
		"dsec="<<dsec<<         // delta sector
		//
		"theta="<<xTheta<<      // theta
		"phi="<<xPhi<<          // phi (alpha)	      
		"z="<<xZ<<              // z
		"snp="<<xSnp<<          // snp
		
		//
		"entries="<<entries<<  // entries in bin
		"mean="<<mean<<        // mean
		"rms="<<rms<<          // rms 
		"refX="<<refX<<        // track matching reference plane
		"\n";
	      delete hisDelta;
	      //printf("%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",isec0, isec1, xPhi,xZ,xSnp, xTheta, entries,mean);
	      //
	    }//ibinTheta
	    delete hisSnp;
	  } //ibinSnp
	  delete hisZ;
	}//ibinZ
	delete hisPhi;
      }//ibinPhi
      delete hisSector1;      
    }//isec1
    delete hisSector0;
  }//isec0
}





// void AliTPCCorrectionFit::MakeLaserDistortionTree(TTree* tree, TObjArray */*corrArray*/, Int_t /*itype*/){
//   //
//   // Make a laser fit tree for global minimization
//   //  
//   AliTPCcalibDB*  calib=AliTPCcalibDB::Instance();  
//   AliTPCCorrection * correction = calib->GetTPCComposedCorrection();  
//   if (!correction) correction = calib->GetTPCComposedCorrection(AliTrackerBase::GetBz());  
//   correction->AddVisualCorrection(correction,0);  //register correction

//   //  AliTPCTransform *transform = AliTPCcalibDB::Instance()->GetTransform() ;
//   //AliTPCParam     *param     = AliTPCcalibDB::Instance()->GetParameters();
//   //
//   const Double_t cutErrY=0.05;
//   const Double_t kSigmaCut=4;
//   //  const Double_t cutErrZ=0.03;
//   const Double_t kEpsilon=0.00000001;
//   //  const Double_t kMaxDist=1.;  // max distance - space correction
//   TVectorD *vecdY=0;
//   TVectorD *vecdZ=0;
//   TVectorD *veceY=0;
//   TVectorD *veceZ=0;
//   AliTPCLaserTrack *ltr=0;
//   AliTPCLaserTrack::LoadTracks();
//   tree->SetBranchAddress("dY.",&vecdY);
//   tree->SetBranchAddress("dZ.",&vecdZ);
//   tree->SetBranchAddress("eY.",&veceY);
//   tree->SetBranchAddress("eZ.",&veceZ);
//   tree->SetBranchAddress("LTr.",&ltr);
//   Int_t entries= tree->GetEntries();
//   TTreeSRedirector *pcstream= new TTreeSRedirector("distortionLaser_0.root");
//   Double_t bz=AliTrackerBase::GetBz();
//   // 
//   //  Double_t globalXYZ[3];
//   //Double_t globalXYZCorr[3];
//   for (Int_t ientry=0; ientry<entries; ientry++){
//     tree->GetEntry(ientry);
//     if (!ltr->GetVecGX()){
//       ltr->UpdatePoints();
//     }
//     //
//     TVectorD fit10(5);
//     TVectorD fit5(5);
//     printf("Entry\t%d\n",ientry);
//     for (Int_t irow0=0; irow0<158; irow0+=1){
//       //       
//       TLinearFitter fitter10(4,"hyp3");
//       TLinearFitter fitter5(2,"hyp1");
//       Int_t sector= (Int_t)(*ltr->GetVecSec())[irow0];
//       if (sector<0) continue;
//       //if (TMath::Abs(vecdY->GetMatrixArray()[irow0])<kEpsilon) continue;

//       Double_t refX= (*ltr->GetVecLX())[irow0];
//       Int_t firstRow1 = TMath::Max(irow0-10,0);
//       Int_t lastRow1  = TMath::Min(irow0+10,158);
//       Double_t padWidth=(irow0<64)?0.4:0.6;
//       // make long range fit
//       for (Int_t irow1=firstRow1; irow1<=lastRow1; irow1++){
// 	if (TMath::Abs((*ltr->GetVecSec())[irow1]-sector)>kEpsilon) continue;
// 	if (veceY->GetMatrixArray()[irow1]>cutErrY) continue;
// 	if (TMath::Abs(vecdY->GetMatrixArray()[irow1])<kEpsilon) continue;
// 	Double_t idealX= (*ltr->GetVecLX())[irow1];
// 	Double_t idealY= (*ltr->GetVecLY())[irow1];
// 	//	Double_t idealZ= (*ltr->GetVecLZ())[irow1];
// 	Double_t gx= (*ltr->GetVecGX())[irow1];
// 	Double_t gy= (*ltr->GetVecGY())[irow1];
// 	Double_t gz= (*ltr->GetVecGZ())[irow1];
// 	Double_t measY=(*vecdY)[irow1]+idealY;
// 	Double_t deltaR = GetCorrXYZ(gx, gy, gz, 0,0);
// 	// deltaR = R distorted -R ideal
// 	Double_t xxx[4]={idealX+deltaR-refX,TMath::Cos(idealY/padWidth), TMath::Sin(idealY/padWidth)};
// 	fitter10.AddPoint(xxx,measY,1);
//       }
//       Bool_t isOK=kTRUE;
//       Double_t rms10=0;//TMath::Sqrt(fitter10.GetChisquare()/(fitter10.GetNpoints()-4));
//       Double_t mean10  =0;//   fitter10.GetParameter(0);
//       Double_t slope10  =0;//   fitter10.GetParameter(0);
//       Double_t cosPart10  = 0;//  fitter10.GetParameter(2);
//       Double_t sinPart10   =0;//  fitter10.GetParameter(3); 

//       if (fitter10.GetNpoints()>10){
// 	fitter10.Eval();
// 	rms10=TMath::Sqrt(fitter10.GetChisquare()/(fitter10.GetNpoints()-4));
// 	mean10      =   fitter10.GetParameter(0);
// 	slope10     =   fitter10.GetParameter(1);
// 	cosPart10   =   fitter10.GetParameter(2);
// 	sinPart10   =  fitter10.GetParameter(3); 
// 	//
// 	// make short range fit
// 	//
// 	for (Int_t irow1=firstRow1+5; irow1<=lastRow1-5; irow1++){
// 	  if (TMath::Abs((*ltr->GetVecSec())[irow1]-sector)>kEpsilon) continue;
// 	  if (veceY->GetMatrixArray()[irow1]>cutErrY) continue;
// 	  if (TMath::Abs(vecdY->GetMatrixArray()[irow1])<kEpsilon) continue;
// 	  Double_t idealX= (*ltr->GetVecLX())[irow1];
// 	  Double_t idealY= (*ltr->GetVecLY())[irow1];
// 	  //	  Double_t idealZ= (*ltr->GetVecLZ())[irow1];
// 	  Double_t gx= (*ltr->GetVecGX())[irow1];
// 	  Double_t gy= (*ltr->GetVecGY())[irow1];
// 	  Double_t gz= (*ltr->GetVecGZ())[irow1];
// 	  Double_t measY=(*vecdY)[irow1]+idealY;
// 	  Double_t deltaR = GetCorrXYZ(gx, gy, gz, 0,0);
// 	  // deltaR = R distorted -R ideal 
// 	  Double_t expY= mean10+slope10*(idealX+deltaR-refX);
// 	  if (TMath::Abs(measY-expY)>kSigmaCut*rms10) continue;
// 	  //
// 	  Double_t corr=cosPart10*TMath::Cos(idealY/padWidth)+sinPart10*TMath::Sin(idealY/padWidth);
// 	  Double_t xxx[4]={idealX+deltaR-refX,TMath::Cos(idealY/padWidth), TMath::Sin(idealY/padWidth)};
// 	  fitter5.AddPoint(xxx,measY-corr,1);
// 	}     
//       }else{
// 	isOK=kFALSE;
//       }
//       if (fitter5.GetNpoints()<8) isOK=kFALSE;

//       Double_t rms5=0;//TMath::Sqrt(fitter5.GetChisquare()/(fitter5.GetNpoints()-4));
//       Double_t offset5  =0;//  fitter5.GetParameter(0);
//       Double_t slope5   =0;//  fitter5.GetParameter(0); 
//       if (isOK){
// 	fitter5.Eval();
// 	rms5=TMath::Sqrt(fitter5.GetChisquare()/(fitter5.GetNpoints()-4));
// 	offset5  =  fitter5.GetParameter(0);
// 	slope5   =  fitter5.GetParameter(0); 
//       }
//       //
//       Double_t dtype=5;
//       Double_t ptype=0;
//       Double_t phi   =(*ltr->GetVecPhi())[irow0];
//       Double_t theta =ltr->GetTgl();
//       Double_t mean=(vecdY)->GetMatrixArray()[irow0];
//       Double_t gx=0,gy=0,gz=0;
//       Double_t snp = (*ltr->GetVecP2())[irow0];
//       Int_t bundle= ltr->GetBundle();
//       Int_t id= ltr->GetId();
//       //      Double_t rms = err->GetMatrixArray()[irow];
//       //
//       gx = (*ltr->GetVecGX())[irow0];
//       gy = (*ltr->GetVecGY())[irow0];
//       gz = (*ltr->GetVecGZ())[irow0];
//       Double_t dRrec = GetCorrXYZ(gx, gy, gz, 0,0);
//       fitter10.GetParameters(fit10);
//       fitter5.GetParameters(fit5);      
//       Double_t idealY= (*ltr->GetVecLY())[irow0];
//       Double_t measY=(*vecdY)[irow0]+idealY;
//       Double_t corr=cosPart10*TMath::Cos(idealY/padWidth)+sinPart10*TMath::Sin(idealY/padWidth);
//       if (TMath::Max(rms5,rms10)>0.06) isOK=kFALSE;
//       //
//       (*pcstream)<<"fitFull"<<  // dumpe also intermediate results
// 	"bz="<<bz<<         // magnetic filed used
// 	"dtype="<<dtype<<   // detector match type
// 	"ptype="<<ptype<<   // parameter type
// 	"theta="<<theta<<   // theta
// 	"phi="<<phi<<       // phi 
// 	"snp="<<snp<<       // snp
// 	"sector="<<sector<<
// 	"bundle="<<bundle<<
// // 	//	"dsec="<<dsec<<
// 	"refX="<<refX<<      // reference radius
// 	"gx="<<gx<<         // global position
// 	"gy="<<gy<<         // global position
// 	"gz="<<gz<<         // global position
// 	"dRrec="<<dRrec<<      // delta Radius in reconstruction
//  	"id="<<id<<     //bundle
// 	"rms10="<<rms10<<
// 	"rms5="<<rms5<<
// 	"fit10.="<<&fit10<<
// 	"fit5.="<<&fit5<<
// 	"measY="<<measY<<
// 	"mean="<<mean<<
// 	"idealY="<<idealY<<
// 	"corr="<<corr<<
// 	"isOK="<<isOK<<
// 	"\n";
//     }
//   }
//   delete pcstream;
// }


void AliTPCCorrectionFit::MakeTrackDistortionTree(TTree *tinput, Int_t dtype, Int_t ptype, const TObjArray * corrArray, Int_t step, Int_t offset, Bool_t debug ){
  //
  // Make a fit tree:
  // For each partial correction (specified in array) and given track topology (phi, theta, snp, refX)
  // calculates partial distortions
  // Partial distortion is stored in the resulting tree
  // Output is storred in the file distortion_<dettype>_<partype>.root
  // Partial  distortion is stored with the name given by correction name
  //
  //
  // Parameters of function:
  // input     - input tree
  // dtype     - distortion type 0 - ITSTPC,  1 -TPCTRD, 2 - TPCvertex , 3 - TPC-TOF,  4 - TPCTPC track crossing 
  // ppype     - parameter type
  // corrArray - array with partial corrections
  // step      - skipe entries  - if 1 all entries processed - it is slow
  // debug     0 if debug on also space points dumped - it is slow

  const Double_t kMaxSnp = 0.85;  
  const Double_t kcutSnp=0.25;
  const Double_t kcutTheta=1.;
  const Double_t kRadiusTPC=85;
  //  AliTPCROC *tpcRoc =AliTPCROC::Instance();  
  //
  const Double_t kMass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  //  const Double_t kB2C=-0.299792458e-3;
  const Int_t kMinEntries=20; 
  Double_t phi,theta, snp, mean,rms, entries,sector,dsec;
  Float_t refX;  
  Int_t run;
  tinput->SetBranchAddress("run",&run);
  tinput->SetBranchAddress("theta",&theta);
  tinput->SetBranchAddress("phi", &phi);
  tinput->SetBranchAddress("snp",&snp);
  tinput->SetBranchAddress("mean",&mean);
  tinput->SetBranchAddress("rms",&rms);
  tinput->SetBranchAddress("entries",&entries);
  tinput->SetBranchAddress("sector",&sector);
  tinput->SetBranchAddress("dsec",&dsec);
  tinput->SetBranchAddress("refX",&refX);
  TTreeSRedirector *pcstream = new TTreeSRedirector(Form("distortion%d_%d_%d.root",dtype,ptype,offset));
  //
  Int_t nentries=tinput->GetEntries();
  Int_t ncorr=corrArray->GetEntries();
  Double_t corrections[100]={0}; //
  Double_t tPar[5];
  Double_t cov[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Int_t dir=0;
  if (dtype==5 || dtype==6) dtype=4;
  if (dtype==0) { dir=-1;}
  if (dtype==1) { dir=1;}
  if (dtype==2) { dir=-1;}
  if (dtype==3) { dir=1;}
  if (dtype==4) { dir=-1;}
  //
  for (Int_t ientry=offset; ientry<nentries; ientry+=step){
    tinput->GetEntry(ientry);
    if (TMath::Abs(snp)>kMaxSnp) continue;
    tPar[0]=0;
    tPar[1]=theta*refX;
    if (dtype==2)  tPar[1]=theta*kRadiusTPC;
    tPar[2]=snp;
    tPar[3]=theta;
    tPar[4]=(gRandom->Rndm()-0.5)*0.02;  // should be calculated - non equal to 0
    if (dtype==4){
      // tracks crossing CE
      tPar[1]=0;   // track at the CE
      //if (TMath::Abs(theta) <0.05) continue;  // deep cross
    }

    if (TMath::Abs(snp) >kcutSnp) continue;
    if (TMath::Abs(theta) >kcutTheta) continue;
    printf("%f\t%f\t%f\t%f\t%f\t%f\n",entries, sector,theta,snp, mean,rms);
    Double_t bz=AliTrackerBase::GetBz();
    if (dtype !=4) { //exclude TPC  - for TPC mainly non primary tracks
      if (dtype!=2 && TMath::Abs(bz)>0.1 )  tPar[4]=snp/(refX*bz*kB2C*2);
      
      if (dtype==2 && TMath::Abs(bz)>0.1 )  {
	tPar[4]=snp/(kRadiusTPC*bz*kB2C*2);//
	// snp at the TPC inner radius in case the vertex match used
      }
    }
    //
    tPar[4]+=(gRandom->Rndm()-0.5)*0.02;
    AliExternalTrackParam track(refX,phi,tPar,cov);
    Double_t xyz[3];
    track.GetXYZ(xyz);
    Int_t id=0;
    Double_t pt=1./tPar[4];
    Double_t dRrec=0; // dummy value - needed for points - e.g for laser
    //if (ptype==4 &&bz<0) mean*=-1;  // interpret as curvature -- COMMENTED out - in lookup signed 1/pt used
    Double_t refXD=refX;
    (*pcstream)<<"fit"<<
      "run="<<run<<       // run number
      "bz="<<bz<<         // magnetic filed used
      "dtype="<<dtype<<   // detector match type
      "ptype="<<ptype<<   // parameter type
      "theta="<<theta<<   // theta
      "phi="<<phi<<       // phi 
      "snp="<<snp<<       // snp
      "mean="<<mean<<     // mean dist value
      "rms="<<rms<<       // rms
      "sector="<<sector<<
      "dsec="<<dsec<<
      "refX="<<refXD<<         // referece X as double
      "gx="<<xyz[0]<<         // global position at reference
      "gy="<<xyz[1]<<         // global position at reference
      "gz="<<xyz[2]<<         // global position at reference	
      "dRrec="<<dRrec<<      // delta Radius in reconstruction
      "pt="<<pt<<            // pt
      "id="<<id<<             // track id
      "entries="<<entries;// number of entries in bin
    //
    Bool_t isOK=kTRUE;
    if (entries<kMinEntries) isOK=kFALSE;
    //
    if (dtype!=4) for (Int_t icorr=0; icorr<ncorr; icorr++) {
      AliTPCCorrection *corr = (AliTPCCorrection*)corrArray->At(icorr);
      corrections[icorr]=0;
      if (entries>kMinEntries){
	AliExternalTrackParam trackIn(refX,phi,tPar,cov);
	AliExternalTrackParam *trackOut = 0;
	if (debug) trackOut=corr->FitDistortedTrack(trackIn, refX, dir,pcstream);
	if (!debug) trackOut=corr->FitDistortedTrack(trackIn, refX, dir,0);
	if (dtype==0) {dir= -1;}
	if (dtype==1) {dir=  1;}
	if (dtype==2) {dir= -1;}
	if (dtype==3) {dir=  1;}
	//
	if (trackOut){
	  if (!AliTrackerBase::PropagateTrackTo(&trackIn,refX,kMass,5,kTRUE,kMaxSnp)) isOK=kFALSE;
	  if (!trackOut->Rotate(trackIn.GetAlpha())) isOK=kFALSE;
	  if (!AliTrackerBase::PropagateTrackTo(trackOut,trackIn.GetX(),kMass,5,kFALSE,kMaxSnp)) isOK=kFALSE;
	  //	  trackOut->PropagateTo(trackIn.GetX(),AliTrackerBase::GetBz());
	  //	  
	  corrections[icorr]= trackOut->GetParameter()[ptype]-trackIn.GetParameter()[ptype];
	  delete trackOut;      
	}else{
	  corrections[icorr]=0;
	  isOK=kFALSE;
	}
	//if (ptype==4 &&bz<0) corrections[icorr]*=-1;  // interpret as curvature - commented out
      }      
      (*pcstream)<<"fit"<<
	Form("%s=",corr->GetName())<<corrections[icorr];   // dump correction value
    }
  
    if (dtype==4) for (Int_t icorr=0; icorr<ncorr; icorr++) {
      //
      // special case of the TPC tracks crossing the CE
      //
      AliTPCCorrection *corr = (AliTPCCorrection*)corrArray->At(icorr);
      corrections[icorr]=0;
      if (entries>kMinEntries){
	AliExternalTrackParam trackIn0(refX,phi,tPar,cov); //Outer - direction to vertex
	AliExternalTrackParam trackIn1(refX,phi,tPar,cov); //Inner - direction magnet 
	AliExternalTrackParam *trackOut0 = 0;
	AliExternalTrackParam *trackOut1 = 0;
	//
	if (debug)  trackOut0=corr->FitDistortedTrack(trackIn0, refX, dir,pcstream);
	if (!debug) trackOut0=corr->FitDistortedTrack(trackIn0, refX, dir,0);
	if (debug)  trackOut1=corr->FitDistortedTrack(trackIn1, refX, -dir,pcstream);
	if (!debug) trackOut1=corr->FitDistortedTrack(trackIn1, refX, -dir,0);
	//
	if (trackOut0 && trackOut1){
	  if (!AliTrackerBase::PropagateTrackTo(&trackIn0,refX,kMass,5,kTRUE,kMaxSnp))  isOK=kFALSE;
	  if (!AliTrackerBase::PropagateTrackTo(&trackIn0,refX,kMass,1,kFALSE,kMaxSnp)) isOK=kFALSE;
	  if (!trackOut0->Rotate(trackIn0.GetAlpha())) isOK=kFALSE;
	  if (!AliTrackerBase::PropagateTrackTo(trackOut0,trackIn0.GetX(),kMass,5,kFALSE,kMaxSnp)) isOK=kFALSE;
	  //
	  if (!AliTrackerBase::PropagateTrackTo(&trackIn1,refX,kMass,5,kTRUE,kMaxSnp)) isOK=kFALSE;
	  if (!trackIn1.Rotate(trackIn0.GetAlpha()))  isOK=kFALSE;
	  if (!AliTrackerBase::PropagateTrackTo(&trackIn1,trackIn0.GetX(),kMass,1,kFALSE,kMaxSnp)) isOK=kFALSE;
	  if (!trackOut1->Rotate(trackIn1.GetAlpha())) isOK=kFALSE;	  
	  if (!AliTrackerBase::PropagateTrackTo(trackOut1,trackIn1.GetX(),kMass,5,kFALSE,kMaxSnp)) isOK=kFALSE;
	  //
	  corrections[icorr] = (trackOut0->GetParameter()[ptype]-trackIn0.GetParameter()[ptype]);
	  corrections[icorr]-= (trackOut1->GetParameter()[ptype]-trackIn1.GetParameter()[ptype]);
	  if (isOK)
	    if ((TMath::Abs(trackOut0->GetX()-trackOut1->GetX())>0.1)||
		(TMath::Abs(trackOut0->GetX()-trackIn1.GetX())>0.1)||
		(TMath::Abs(trackOut0->GetAlpha()-trackOut1->GetAlpha())>0.00001)||
		(TMath::Abs(trackOut0->GetAlpha()-trackIn1.GetAlpha())>0.00001)||
		(TMath::Abs(trackIn0.GetTgl()-trackIn1.GetTgl())>0.0001)||
		(TMath::Abs(trackIn0.GetSnp()-trackIn1.GetSnp())>0.0001)
		){
	      isOK=kFALSE;
	    }	  	  
	  delete trackOut0;      
	  delete trackOut1;    	  
	}else{
	  corrections[icorr]=0;
	  isOK=kFALSE;
	}
	//
	//if (ptype==4 &&bz<0) corrections[icorr]*=-1;  // interpret as curvature - commented out no in lookup
      }      
      (*pcstream)<<"fit"<<
	Form("%s=",corr->GetName())<<corrections[icorr];   // dump correction value
    }
    //
    (*pcstream)<<"fit"<<"isOK="<<isOK<<"\n";
  }


  delete pcstream;
}



void AliTPCCorrectionFit::MakeSectorDistortionTree(TTree *tinput, Int_t dtype, Int_t ptype, const TObjArray * corrArray, Int_t step, Int_t offset, Bool_t debug ){
  //
  // Make a fit tree:
  // For each partial correction (specified in array) and given track topology (phi, theta, snp, refX)
  // calculates partial distortions
  // Partial distortion is stored in the resulting tree
  // Output is storred in the file distortion_<dettype>_<partype>.root
  // Partial  distortion is stored with the name given by correction name
  //
  //
  // Parameters of function:
  // input     - input tree
  // dtype     - distortion type 10 - IROC-OROC 
  // ppype     - parameter type
  // corrArray - array with partial corrections
  // step      - skipe entries  - if 1 all entries processed - it is slow
  // debug     0 if debug on also space points dumped - it is slow

  const Double_t kMaxSnp = 0.8;  
  const Int_t kMinEntries=200; 
  //  AliTPCROC *tpcRoc =AliTPCROC::Instance();  
  //
  const Double_t kMass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  //  const Double_t kB2C=-0.299792458e-3;
  Double_t phi,theta, snp, mean,rms, entries,sector,dsec,globalZ;
  Int_t isec1, isec0;
  Double_t refXD;
  Float_t refX;
  Int_t run;
  tinput->SetBranchAddress("run",&run);
  tinput->SetBranchAddress("theta",&theta);
  tinput->SetBranchAddress("phi", &phi);
  tinput->SetBranchAddress("snp",&snp);
  tinput->SetBranchAddress("mean",&mean);
  tinput->SetBranchAddress("rms",&rms);
  tinput->SetBranchAddress("entries",&entries);
  tinput->SetBranchAddress("sector",&sector);
  tinput->SetBranchAddress("dsec",&dsec);
  tinput->SetBranchAddress("refX",&refXD);
  tinput->SetBranchAddress("z",&globalZ);
  tinput->SetBranchAddress("isec0",&isec0);
  tinput->SetBranchAddress("isec1",&isec1);
  TTreeSRedirector *pcstream = new TTreeSRedirector(Form("distortionSector%d_%d_%d.root",dtype,ptype,offset));
  //
  Int_t nentries=tinput->GetEntries();
  Int_t ncorr=corrArray->GetEntries();
  Double_t corrections[100]={0}; //
  Double_t tPar[5];
  Double_t cov[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Int_t dir=0;
  //
  for (Int_t ientry=offset; ientry<nentries; ientry+=step){
    tinput->GetEntry(ientry);
    refX=refXD;
    Int_t id=-1;
    if (TMath::Abs(TMath::Abs(isec0%18)-TMath::Abs(isec1%18))==0) id=1;  // IROC-OROC - opposite side
    if (TMath::Abs(TMath::Abs(isec0%36)-TMath::Abs(isec1%36))==0) id=2;  // IROC-OROC - same side
    if (dtype==10  && id==-1) continue;
    //
    dir=-1;
    tPar[0]=0;
    tPar[1]=globalZ;
    tPar[2]=snp;
    tPar[3]=theta;
    tPar[4]=(gRandom->Rndm()-0.1)*0.2;  //
    Double_t pt=1./tPar[4];
    //
    printf("%f\t%f\t%f\t%f\t%f\t%f\n",entries, sector,theta,snp, mean,rms);
    Double_t bz=AliTrackerBase::GetBz();
    AliExternalTrackParam track(refX,phi,tPar,cov);    
    Double_t xyz[3],xyzIn[3],xyzOut[3];
    track.GetXYZ(xyz);
    track.GetXYZAt(85,bz,xyzIn);    
    track.GetXYZAt(245,bz,xyzOut);    
    Double_t phiIn  = TMath::ATan2(xyzIn[1],xyzIn[0]);
    Double_t phiOut = TMath::ATan2(xyzOut[1],xyzOut[0]);
    Double_t phiRef = TMath::ATan2(xyz[1],xyz[0]);
    Int_t sectorRef = TMath::Nint(9.*phiRef/TMath::Pi()-0.5);
    Int_t sectorIn  = TMath::Nint(9.*phiIn/TMath::Pi()-0.5);
    Int_t sectorOut = TMath::Nint(9.*phiOut/TMath::Pi()-0.5);
    //
    Bool_t isOK=kTRUE; 
    if (sectorIn!=sectorOut) isOK=kFALSE;  // requironment - cluster in the same sector
    if (sectorIn!=sectorRef) isOK=kFALSE;  // requironment - cluster in the same sector
    if (entries<kMinEntries/(1+TMath::Abs(globalZ/100.))) isOK=kFALSE;  // requironment - minimal amount of tracks in bin
    // Do downscale
    if (TMath::Abs(theta)>1) isOK=kFALSE;
    //
    Double_t dRrec=0; // dummy value - needed for points - e.g for laser
    //
    (*pcstream)<<"fit"<<
      "run="<<run<<       //run
      "bz="<<bz<<         // magnetic filed used
      "dtype="<<dtype<<   // detector match type
      "ptype="<<ptype<<   // parameter type
      "theta="<<theta<<   // theta
      "phi="<<phi<<       // phi 
      "snp="<<snp<<       // snp
      "mean="<<mean<<     // mean dist value
      "rms="<<rms<<       // rms
      "sector="<<sector<<
      "dsec="<<dsec<<
      "refX="<<refXD<<         // referece X
      "gx="<<xyz[0]<<         // global position at reference
      "gy="<<xyz[1]<<         // global position at reference
      "gz="<<xyz[2]<<         // global position at reference	
      "dRrec="<<dRrec<<      // delta Radius in reconstruction
      "pt="<<pt<<      //pt
      "id="<<id<<             // track id
      "entries="<<entries;// number of entries in bin
    //
    AliExternalTrackParam *trackOut0 = 0;
    AliExternalTrackParam *trackOut1 = 0;
    AliExternalTrackParam *ptrackIn0 = 0;
    AliExternalTrackParam *ptrackIn1 = 0;

    for (Int_t icorr=0; icorr<ncorr; icorr++) {
      //
      // special case of the TPC tracks crossing the CE
      //
      AliTPCCorrection *corr = (AliTPCCorrection*)corrArray->At(icorr);
      corrections[icorr]=0;
      if (entries>kMinEntries &&isOK){
	AliExternalTrackParam trackIn0(refX,phi,tPar,cov);
	AliExternalTrackParam trackIn1(refX,phi,tPar,cov);
	ptrackIn1=&trackIn0;
	ptrackIn0=&trackIn1;
	//
	if (debug)  trackOut0=corr->FitDistortedTrack(trackIn0, refX, dir,pcstream);
	if (!debug) trackOut0=corr->FitDistortedTrack(trackIn0, refX, dir,0);
	if (debug)  trackOut1=corr->FitDistortedTrack(trackIn1, refX, -dir,pcstream);
	if (!debug) trackOut1=corr->FitDistortedTrack(trackIn1, refX, -dir,0);
	//
	if (trackOut0 && trackOut1){
	  //
	  if (!AliTrackerBase::PropagateTrackTo(&trackIn0,refX,kMass,1,kTRUE,kMaxSnp))  isOK=kFALSE;
	  if (!AliTrackerBase::PropagateTrackTo(&trackIn0,refX,kMass,1,kFALSE,kMaxSnp)) isOK=kFALSE;
	  // rotate all tracks to the same frame
	  if (!trackOut0->Rotate(trackIn0.GetAlpha())) isOK=kFALSE;
	  if (!trackIn1.Rotate(trackIn0.GetAlpha()))  isOK=kFALSE;
	  if (!trackOut1->Rotate(trackIn0.GetAlpha())) isOK=kFALSE;	  
	  //
	  if (!AliTrackerBase::PropagateTrackTo(trackOut0,refX,kMass,1,kFALSE,kMaxSnp)) isOK=kFALSE;
	  if (!AliTrackerBase::PropagateTrackTo(&trackIn1,refX,kMass,1,kFALSE,kMaxSnp)) isOK=kFALSE;
	  if (!AliTrackerBase::PropagateTrackTo(trackOut1,refX,kMass,1,kFALSE,kMaxSnp)) isOK=kFALSE;
	  //
	  corrections[icorr] = (trackOut0->GetParameter()[ptype]-trackIn0.GetParameter()[ptype]);
	  corrections[icorr]-= (trackOut1->GetParameter()[ptype]-trackIn1.GetParameter()[ptype]);
	  (*pcstream)<<"fitDebug"<< // just to debug the correction
	    "mean="<<mean<<
	    "pIn0.="<<ptrackIn0<<
	    "pIn1.="<<ptrackIn1<<
	    "pOut0.="<<trackOut0<<
	    "pOut1.="<<trackOut1<<
	    "refX="<<refXD<<
	    "\n";
	  delete trackOut0;      
	  delete trackOut1;      
	}else{
	  corrections[icorr]=0;
	  isOK=kFALSE;
	}
      }      
      (*pcstream)<<"fit"<<
	Form("%s=",corr->GetName())<<corrections[icorr];   // dump correction value
    }
    //
    (*pcstream)<<"fit"<<"isOK="<<isOK<<"\n";
  }
  delete pcstream;
}


