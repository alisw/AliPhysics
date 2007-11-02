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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//     Class to analyse tracks for calibration                               //
//     to be used as a component in AliTPCSelectorTracks                     //
//     In the constructor you have to specify name and title                 //
//     to get the Object out of a file.                                      //
//     The parameter 'clusterParam', a AliTPCClusterParam object             //
//      (needed for TPC cluster error and shape parameterization)            //
//     Normally you get this object out of the file 'TPCClusterParam.root'   //
//     In the parameter 'cuts' the cuts are specified, that decide           //
//     weather a track will be accepted for calibration or not.              //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//
// ROOT includes 
//
#include <iostream>
#include <fstream>
using namespace std;

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH3F.h>
//
#include <TPDGCode.h>
#include <TStyle.h>
#include "TLinearFitter.h"
#include "TMatrixD.h"
#include "TTreeStream.h"
#include "TF1.h"
#include <TCanvas.h>
#include <TGraph2DErrors.h>
#include "TPostScript.h"
#include "TCint.h"

//
// AliROOT includes 
//
#include "AliMagF.h"
#include "AliTracker.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h" 
#include "AliTPCseed.h"
#include "AliTPCclusterMI.h"
#include "AliTPCROC.h"

#include "AliTPCParamSR.h"
#include "AliTrackPointArray.h"
#include "AliTPCcalibTracks.h"
#include "AliTPCClusterParam.h"


ClassImp(AliTPCcalibTracks)

AliTPCParam  param;


AliTPCcalibTracks::AliTPCcalibTracks():TNamed() {
   // 
   // AliTPCcalibTracks default constructor
   //    
   cout << "AliTPCcalibTracks' default constructor called" << endl;
   fClusterParam = 0;
   fArrayAmpRow  = 0;
   fArrayAmp     = 0; 
   fArrayQDY     = 0; 
   fArrayQDZ     = 0; 
   fArrayQRMSY   = 0;
   fArrayQRMSZ   = 0;
   fDeltaY       = 0;
   fDeltaZ       = 0;
   fResolY       = 0;
   fResolZ       = 0;
   fRMSY         = 0;
   fRMSZ         = 0;
   fHclus        = 0;
   fROC          = 0;
   fCuts         = 0;
}   


AliTPCcalibTracks::AliTPCcalibTracks(AliTPCcalibTracks* ct){
   // 
   // AliTPCcalibTracks copy constructor
   // 
   
   cout << " ***** this is AliTPCcalibTracks' copy constructor ***** " << endl;
   
   Int_t length = ct->fArrayAmpRow->GetEntriesFast();
   fArrayAmpRow = new TObjArray(length);
   fArrayAmp = new TObjArray(length);
   for (Int_t i = 0; i < length; i++) {
      fArrayAmpRow->AddAt( (TProfile*)ct->fArrayAmpRow->At(i)->Clone(), i);
      fArrayAmp->AddAt( ((TProfile*)ct->fArrayAmp->At(i)->Clone()), i);
   }
   
   length = ct->fArrayQDY->GetEntriesFast();
   fArrayQDY= new TObjArray(length);
   fArrayQDZ= new TObjArray(length);
   fArrayQRMSY= new TObjArray(length);
   fArrayQRMSZ= new TObjArray(length);
   for (Int_t i = 0; i < length; i++) {
      fArrayQDY->AddAt( ((TH1F*)ct->fArrayQDY->At(i)->Clone()), i);
      fArrayQDZ->AddAt( ((TH1F*)ct->fArrayQDZ->At(i)->Clone()), i);
      fArrayQRMSY->AddAt( ((TH1F*)ct->fArrayQRMSY->At(i)->Clone()), i);
      fArrayQRMSZ->AddAt( ((TH1F*)ct->fArrayQRMSZ->At(i)->Clone()), i);
   }
   
   length = ct->fResolY->GetEntriesFast();
   fResolY = new TObjArray(length);
   fResolZ = new TObjArray(length);
   fRMSY = new TObjArray(length);
   fRMSZ = new TObjArray(length);
   for (Int_t i = 0; i < length; i++) {
      fResolY->AddAt( ((TH1F*)ct->fResolY->At(i)->Clone()), i);
      fResolZ->AddAt( ((TH1F*)ct->fResolZ->At(i)->Clone()), i);
      fRMSY->AddAt( ((TH1F*)ct->fRMSY->At(i)->Clone()), i);
      fRMSZ->AddAt( ((TH1F*)ct->fRMSZ->At(i)->Clone()), i);
   } 
   
   fDeltaY =  (TH1F*)ct->fDeltaY->Clone();
   fDeltaZ =  (TH1F*)ct->fDeltaZ->Clone();
   
   fHclus = (TH1I*)ct->fHclus->Clone();
   
   fCuts = new AliTPCcalibTracksCuts(ct->fCuts->GetMinClusters(), ct->fCuts->GetMinRatio(), 
      ct->fCuts->GetMax1pt(), ct->fCuts->GetEdgeYXCutNoise(), ct->fCuts->GetEdgeThetaCutNoise());
}


AliTPCcalibTracks::~AliTPCcalibTracks() {
   // 
   // AliTPCcalibTracks destructor
   // 
   cout << "AliTPCcalibTracks' destuctor called." << endl;
   Int_t length = 0;
   if (fArrayAmpRow) length = fArrayAmpRow->GetEntriesFast();
   for (Int_t i = 0; i < length; i++){
      delete fArrayAmpRow->At(i);  
      delete fArrayAmp->At(i);  
   }
   delete fArrayAmpRow;
   delete fArrayAmp;
   
   delete fDeltaY;
   delete fDeltaZ;
   
   if (fResolY) length = fResolY->GetEntriesFast();
   for (Int_t i = 0; i < length; i++){
      delete fResolY->At(i);
      delete fResolZ->At(i);
      delete fRMSY->At(i);
      delete fRMSZ->At(i);
   }
   delete fResolY;
   delete fResolZ;
   delete fRMSY;
   delete fRMSZ;
   
   if (fArrayQDY) length = fArrayQDY->GetEntriesFast();
   for (Int_t i = 0; i < length; i++){
      delete fArrayQDY->At(i);
      delete fArrayQDZ->At(i);
      delete fArrayQRMSY->At(i);
      delete fArrayQRMSZ->At(i);
   }
   delete fArrayQDY;
   delete fArrayQDZ;
   delete fArrayQRMSY;
   delete fArrayQRMSZ;   
   delete fHclus;
   //delete fDebugStream;
   //fDebugStream->Close();;
}


AliTPCcalibTracks::AliTPCcalibTracks(const Text_t *name, const Text_t *title, AliTPCClusterParam *clusterParam,  AliTPCcalibTracksCuts* cuts) : 
   TNamed(name, title),
   fHclus(0)
   // 
   // AliTPCcalibTracks constructor
   // specify 'name' and 'title' of your object
   // specify 'clusterParam', (needed for TPC cluster error and shape parameterization)
   // In the parameter 'cuts' the cuts are specified, that decide           
   // weather a track will be accepted for calibration or not.              
   // 
   // All histograms are instatiated in this constructor.
   // 
 {
   G__SetCatchException(0);     
   param.Update();
   
   fClusterParam = clusterParam;
   if (fClusterParam){
     fClusterParam->SetInstance(fClusterParam);
   }
   else {
     printf("Cluster Param not found\n");
   } 
   fCuts = cuts;
   fDebugStream = new TTreeSRedirector("TPCSelectorDebug.root");     // needs investigation !!!!!
   
   TH1::AddDirectory(kFALSE);
   
   char chname[1000];
   TProfile * prof1=0;
   TH1F     * his1 =0;
   fHclus = new TH1I("hclus","Number of clusters",100,0,200);     // valgrind 3
   
   // Amplitude  - sector - row histograms 
   fArrayAmpRow = new TObjArray(72);
   fArrayAmp    = new TObjArray(72);
   for (Int_t i = 0; i < 36; i++){   
      sprintf(chname,"Amp_row_Sector%d",i);
      prof1 = new TProfile(chname,chname,63,0,64);          // valgrind 3   193,536 bytes in 354 blocks are still reachable 
      prof1->SetXTitle("Pad row");
      prof1->SetYTitle("Mean Max amplitude");
      fArrayAmpRow->AddAt(prof1,i);
      sprintf(chname,"Amp_row_Sector%d",i+36);
      prof1 = new TProfile(chname,chname,96,0,97);       // valgrind 3   3,912 bytes in 6 blocks are possibly lost
      prof1->SetXTitle("Pad row");
      prof1->SetYTitle("Mean Max  amplitude");
      fArrayAmpRow->AddAt(prof1,i+36);
      
      // amplitude
      sprintf(chname,"Amp_Sector%d",i);
      his1 = new TH1F(chname,chname,250,0,500);         // valgrind 
      his1->SetXTitle("Max Amplitude (ADC)");
      fArrayAmp->AddAt(his1,i);
      sprintf(chname,"Amp_Sector%d",i+36);
      his1 = new TH1F(chname,chname,200,0,600);         // valgrind 3   13,408,208 bytes in 229 blocks are still reachable
      his1->SetXTitle("Max Amplitude (ADC)");
      fArrayAmp->AddAt(his1,i+36);
   }
   
   TH1::AddDirectory(kFALSE);
   
   fDeltaY = new TH1F("DeltaY","DeltaY",100,-1,1);
   fDeltaZ = new TH1F("DeltaZ","DeltaZ",100,-1,1);
   
   fResolY = new TObjArray(3);
   fResolZ = new TObjArray(3);
   fRMSY   = new TObjArray(3);
   fRMSZ   = new TObjArray(3);
   TH3F * his3D;
   //
   his3D = new TH3F("Resol Y0","Resol Y0", 5,20,250, 4, 0,1., 50, -1,1);
   fResolY->AddAt(his3D,0);	
   his3D = new TH3F("Resol Y1","Resol Y1", 5,20,250, 4, 0,1., 50, -1,1);
   fResolY->AddAt(his3D,1);
   his3D = new TH3F("Resol Y2","Resol Y2", 5,20,250, 4, 0,0.8, 50, -1,1);
   fResolY->AddAt(his3D,2);
   //
   his3D = new TH3F("Resol Z0","Resol Z0", 5,20,250, 4, 0,1, 50, -1,1);
   fResolZ->AddAt(his3D,0);
   his3D = new TH3F("Resol Z1","Resol Z1", 5,20,250, 4, 0,1, 50, -1,1);
   fResolZ->AddAt(his3D,1);
   his3D = new TH3F("Resol Z2","Resol Z2", 5,20,250, 4, 0,1, 50, -1,1);
   fResolZ->AddAt(his3D,2);
   //
   his3D = new TH3F("RMS Y0","RMS Y0", 5,20,250, 4, 0,1., 50, 0,0.8);
   fRMSY->AddAt(his3D,0);
   his3D = new TH3F("RMS Y1","RMS Y1", 5,20,250, 4, 0,1., 50, 0,0.8);
   fRMSY->AddAt(his3D,1);
   his3D = new TH3F("RMS Y2","RMS Y2", 5,20,250, 4, 0,0.8, 50, 0,0.8);
   fRMSY->AddAt(his3D,2);
   //
   his3D = new TH3F("RMS Z0","RMS Z0", 5,20,250, 4, 0,1, 50, 0,0.8);
   fRMSZ->AddAt(his3D,0);
   his3D = new TH3F("RMS Z1","RMS Z1", 5,20,250, 4, 0,1, 50, 0,0.8);
   fRMSZ->AddAt(his3D,1);
   his3D = new TH3F("RMS Z2","RMS Z2", 5,20,250, 4, 0,1, 50, 0,0.8);
   fRMSZ->AddAt(his3D,2);
   //
      
   TH1::AddDirectory(kFALSE);
   
   fArrayQDY = new TObjArray(300);
   fArrayQDZ = new TObjArray(300);
   fArrayQRMSY = new TObjArray(300);
   fArrayQRMSZ = new TObjArray(300);
   for (Int_t iq = 0; iq <= 10; iq++){
      for (Int_t ipad = 0; ipad < 3; ipad++){
         Int_t   bin   = GetBin(iq, ipad);
         Float_t qmean = GetQ(bin);
         char name[200];
         sprintf(name,"ResolY Pad%d Qmiddle%f",ipad, qmean);
         his3D = new TH3F(name, name, 20,10,250, 20, 0,1.5, 50, -1,1);
         fArrayQDY->AddAt(his3D, bin);
         sprintf(name,"ResolZ Pad%d Qmiddle%f",ipad, qmean);
         his3D = new TH3F(name, name, 20,10,250, 20, 0,1.5, 50, -1,1);
         fArrayQDZ->AddAt(his3D, bin);
         sprintf(name,"RMSY Pad%d Qmiddle%f",ipad, qmean);
         his3D = new TH3F(name, name, 20,10,250, 20, 0,1.5, 50, 0,1);
         fArrayQRMSY->AddAt(his3D, bin);
         sprintf(name,"RMSZ Pad%d Qmiddle%f",ipad, qmean);
         his3D = new TH3F(name, name, 20,10,250, 20, 0,1.5, 50, 0,1);
         fArrayQRMSZ->AddAt(his3D, bin);
      }
   } 
}    
   
  
void AliTPCcalibTracks::AddInfo(TChain * chain, char* fileName){
   // 
   // Add the neccessary information for process to the chain 
   // (cluster parametrization)
   // 
   TFile clusterParamFile(fileName);
   AliTPCClusterParam *clusterParam  =  (AliTPCClusterParam *) clusterParamFile.Get("Param");
   chain->GetUserInfo()->AddLast((TObject*)clusterParam);
}

void AliTPCcalibTracks::AddCuts(TChain * chain, char* ctype){
   // 
   // add predefined cuts to the chain for processing
   // (creates AliTPCcalibTracksCuts object)
   // the cuts are set in the following order:
   // fMinClusters (number of clusters)
   // fMinRatio 
   // fMax1pt   1  over p_t
   // fEdgeYXCutNoise
   // fEdgeThetaCutNoise
   // 
   // The following predefined sets of cuts can be selected:
   // laser:      20, 0.4, 0.5, 0.13, 0.018
   // cosmic:     20, 0.4, 0.5, 0.13, 0.01
   // lowflux:    20, 0.4, 5, 0.2, 0.0001
   // highflux:   20, 0.4, 5, 0.2, 0.0001
   // 
   
   TString cutType(ctype);
   cutType.ToUpper();
   AliTPCcalibTracksCuts *cuts = 0;
   if (cutType == "LASER")
      cuts = new AliTPCcalibTracksCuts(20, 0.4, 0.5, 0.13, 0.018);
   else if (cutType == "COSMIC")
      cuts = new AliTPCcalibTracksCuts(20, 0.4, 0.5, 0.13, 0.018);
   else if (cutType == "LOWFLUX")
      cuts = new AliTPCcalibTracksCuts(20, 0.4, 5, 0.2, 0.0001);
   else if (cutType == "HIGHFLUX")
      cuts = new AliTPCcalibTracksCuts(20, 0.4, 5, 0.2, 0.0001);
   else {
      cuts = new AliTPCcalibTracksCuts(20, 0.4, 5, 0.2, 0.0001);
      cerr << "WARNING! unknown type '" << ctype << "', cuts set to default values for cosmics." << endl;
   }
   chain->GetUserInfo()->AddLast(cuts);
}

   
void AliTPCcalibTracks::Process(AliTPCseed *track, AliESDtrack *esd){
   // 
   // To be called in the selector
   // first AcceptTrack is evaluated, then calls all the following analyse functions: 
   // FillResolutionHistoLocal(track)
   // AlignUpDown(track, esd)
   // 
   if (AcceptTrack(track)) {
      FillResolutionHistoLocal(track);
      AlignUpDown(track, esd);
   }
}



Int_t AliTPCcalibTracks::GetBin(Float_t q, Int_t pad){
  //
  // calculate bins for given q and pad type 
  // used in TObjArray
  //
  Int_t res = TMath::Max( TMath::Nint((TMath::Sqrt(q) - 3.)), 0 );  
  res *= 3;
  res += pad;
  return res;
}


Int_t AliTPCcalibTracks::GetBin(Int_t iq, Int_t pad){
  //
  // calculate bins for given iq and pad type 
  // used in TObjArray
  //
  return iq*3+pad;;
}


Float_t AliTPCcalibTracks::GetQ(Int_t bin){
   // 
   // (bin / 3 + 3)^2
   // 
   Int_t bin0 = bin / 3;
   bin0 += 3;
   return bin0 * bin0;
}


Float_t AliTPCcalibTracks::TPCBetheBloch(Float_t bg){
   //
   // Bethe-Bloch energy loss formula
   //
   const Double_t kp1=0.76176e-1;
   const Double_t kp2=10.632;
   const Double_t kp3=0.13279e-4;
   const Double_t kp4=1.8631;
   const Double_t kp5=1.9479;
   Double_t dbg = (Double_t) bg;
   Double_t beta = dbg/TMath::Sqrt(1.+dbg*dbg);
   Double_t aa = TMath::Power(beta,kp4);
   Double_t bb = TMath::Power(1./dbg,kp5);
   bb=TMath::Log(kp3+bb);
   return ((Float_t)((kp2-aa-bb)*kp1/aa));
}


Bool_t AliTPCcalibTracks::AcceptTrack(AliTPCseed * track){
  //
  // Function, that decides wheather a given track is accepted for 
  // the analysis or not. 
  // The cuts are specified in the AliTPCcalibTracksCuts object 'fCuts'
  //
  const Int_t   kMinClusters  = fCuts->GetMinClusters();
  const Float_t kMinRatio     = fCuts->GetMinRatio();
  const Float_t kMax1pt       = fCuts->GetMax1pt();
  const Float_t kEdgeYXCutNoise    = fCuts->GetEdgeYXCutNoise();
  const Float_t kEdgeThetaCutNoise = fCuts->GetEdgeThetaCutNoise();
  //
  // edge induced noise tracks - NEXT RELEASE will be removed during tracking
  if ( TMath::Abs(track->GetY() / track->GetX()) > kEdgeYXCutNoise )
    if ( TMath::Abs(track->GetTgl()) < kEdgeThetaCutNoise ) return kFALSE;
  if (track->GetNumberOfClusters() < kMinClusters) return kFALSE;
  Float_t ratio = track->GetNumberOfClusters() / (track->GetNFoundable() + 1.);
  if (ratio < kMinRatio) return kFALSE;
  Float_t mpt = track->GetSigned1Pt();
  if (TMath::Abs(mpt) > kMax1pt) return kFALSE;
  //if (TMath::Abs(track->GetZ())>240.) return kFALSE;
  //if (TMath::Abs(track->GetZ())<10.) return kFALSE;
  //if (TMath::Abs(track->GetTgl())>0.03) return kFALSE;
  
  return kTRUE;
}


void AliTPCcalibTracks::FillHistoCluster(AliTPCseed * track){
  // 
  // fill fArrayAmpRow
  // 72 TProfiles, one for each ROC with amplitudes vs. row
  // Is this function used somewhere???
  // 
  const Int_t kFirstLargePad = 127;
  const Float_t kLargePadSize = 1.5;
  for (Int_t irow = 0; irow < 159; irow++){
    AliTPCclusterMI * cluster = track->GetClusterPointer(irow);
    if (!cluster) continue;
    Int_t sector = cluster->GetDetector();
    if (cluster->GetQ() <= 0) continue;
    Float_t max = cluster->GetMax();
    printf ("irow, kFirstLargePad = %d, %d \n", irow, kFirstLargePad);
    if ( irow >= kFirstLargePad) {
      max /= kLargePadSize;
    }
    TProfile *profAmpRow = (TProfile*)fArrayAmpRow->At(sector);
    profAmpRow->Fill(cluster->GetRow(), max);
  }  
}


void  AliTPCcalibTracks::FillResolutionHistoLocal(AliTPCseed * track){
   //
   // fill resolution histograms - localy - tracklet in the neighborhood
   // write debug information to 'TPCSelectorDebug.root'
   // 
   // 
   // loop over all padrows along the track
   // fit tracklets (length: 13 rows) calculate mean chi^2 for this track-fit in Y and Z direction
   // 
   // loop again over all padrows along the track
   // fit tracklet (clusters in given padrow +- kDelta padrows) 
   // with polynom of 2nd order and two polynoms of 1st order
   // take both polynoms of 1st order, calculate difference of their parameters
   // add covariance matrixes and calculate chi2 of this difference
   // if this chi2 is bigger than a given threshold, assume that the current cluster is
   // a kink an goto next padrow
   // if not:
   // fill fArrayAmpRow, array with amplitudes vs. row for given sector
   // fill fArrayAmp, array with amplitude histograms for give sector
   // fill fRMSY, fRMSZ, fArrayQRMSY and fArrayQRMSZ, fDeltaY, fDeltaZ, fResolY, fResolZ, fArrayQDY, fArrayQDY
   // 
   // write debug information to 'TPCSelectorDebug.root'
   // 

  
  const Int_t   kDelta    = 10;          // delta rows to fit
  const Float_t kMinRatio = 0.75;        // minimal ratio
  const Float_t kCutChi2  = 6.;          // cut chi2 - left right  - kink removal
  const Float_t kErrorFraction = 0.5;    // use only clusters with small interpolation error - for error param
  const Int_t   kFirstLargePad = 127;    // medium pads -> long pads
  const Float_t kLargePadSize  = 1.5;    // factor between medium and long pads' area
  static TLinearFitter fitterY0(2,"pol1");
  static TLinearFitter fitterZ0(2,"pol1");
  static TLinearFitter fitterY1(2,"pol1");
  static TLinearFitter fitterZ1(2,"pol1");   // valgrind 3   20,484 bytes in 435 blocks are indirectly lost
  static TLinearFitter fitterY2(3,"pol2");
  static TLinearFitter fitterZ2(3,"pol2");
  TVectorD paramY0(2);
  TVectorD paramZ0(2);
  TVectorD paramY1(2);
  TVectorD paramZ1(2);
  TVectorD paramY2(3);
  TVectorD paramZ2(3);
  TMatrixD matrixY0(2,2);
  TMatrixD matrixZ0(2,2);
  TMatrixD matrixY1(2,2);
  TMatrixD matrixZ1(2,2);
  
  // estimate mean error
  Int_t nTrackletsAll = 0;
  Int_t nClusters     = 0;
  Float_t csigmaY     = 0;
  Float_t csigmaZ     = 0;
  Int_t sectorG       = -1;
  
  fHclus->Fill(track->GetNumberOfClusters());
  
  for (Int_t irow = 0; irow < 159; irow++){
    // loop over all rows along the track
    // fit tracklets (length: 13 rows) with pol2 in Y and Z direction
    // calculate mean chi^2 for this track-fit in Y and Z direction
    AliTPCclusterMI * cluster0 = track->GetClusterPointer(irow);
    if (!cluster0) continue;  // no cluster found
    Int_t sector = cluster0->GetDetector();
    if (sector != sectorG){
      // track leaves sector before it crossed enough rows to fit / initialization
      nClusters = 0;
      fitterY2.ClearPoints();
      fitterZ2.ClearPoints();
      sectorG = sector;
    }
    else {
      nClusters++;
      Double_t x = cluster0->GetX();
      fitterY2.AddPoint(&x, cluster0->GetY(), 1);
      fitterZ2.AddPoint(&x, cluster0->GetZ(), 1);
      //
      if ( nClusters >= kDelta + 3 ){  
        // if more than 13 (kDelta+3) rows / clusters were added to the fitters
        // fit the tracklet, increase trackletCounter
	fitterY2.Eval();
	fitterZ2.Eval();
	nTrackletsAll++;
	csigmaY += fitterY2.GetChisquare() / (nClusters - 3.);
	csigmaZ += fitterZ2.GetChisquare() / (nClusters - 3.);
	nClusters = -1;
	fitterY2.ClearPoints();
	fitterZ2.ClearPoints();
      }
    }
  }      // for (Int_t irow = 0; irow < 159; irow++)
  // mean chi^2 for all tracklet fits in Y and in Z direction
  csigmaY = TMath::Sqrt(csigmaY / nTrackletsAll);
  csigmaZ = TMath::Sqrt(csigmaZ / nTrackletsAll);
  
  //
  //
  //
  for (Int_t irow = 0; irow < 159; irow++){
    // loop again over all rows along the track
    // do analysis
    // 
    Int_t nclFound     = 0;
//     Int_t nclFoundable = 0;
    AliTPCclusterMI * cluster0 = track->GetClusterPointer(irow);
    if (!cluster0) continue;
    Int_t sector = cluster0->GetDetector();
    Float_t xref = cluster0->GetX();
    
    
    // what is the following loop good for?????
/*    
    // check the neighborhood occupancy - (Delta ray - noise removal)
    for (Int_t idelta = -kDelta; idelta <= kDelta; idelta++){
      // loop over irow +- kDelta rows (neighboured rows)
      // increase nclFoundable and nclFound
      // nclFoundable == nclFound !!!!!!!!!
      if (idelta == 0) continue;    // check neighbourhood rows, not the row itself
      if (idelta + irow < 0 || idelta + irow > 159) continue;   // don't go out of ROC
      AliTPCclusterMI * clusterD = track->GetClusterPointer(irow);
      if (!clusterD) continue;      // no cluster found in row
      if ( clusterD->GetDetector() != sector) continue;     // track leaves ROC
      if (clusterD->GetType() < 0) continue;      
      nclFoundable++;                                             // ???????????????????????????????????
      nclFound++;                                                 // nclFoundable == nclFound !!!!!!!!!
    }    // neighbourhood-loop
    
    if (nclFound < kDelta * kMinRatio) continue;    // if not enough clusters found in neighbourhood
    if ( Float_t(nclFound) / Float_t(nclFoundable) < kMinRatio ) continue;  // if ratio between foundable and found clusters is too bad
  
*/    
    // Make Fit
    fitterY2.ClearPoints();
    fitterZ2.ClearPoints();
    fitterY0.ClearPoints();
    fitterZ0.ClearPoints();
    fitterY1.ClearPoints();
    fitterZ1.ClearPoints();
    nclFound   = 0;
    Int_t ncl0 = 0;
    Int_t ncl1 = 0;
    
   // fit tracklet (clusters in given padrow +- kDelta padrows) 
   // with polynom of 2nd order and two polynoms of 1st order
   // take both polynoms of 1st order, calculate difference of their parameters
   // add covariance matrixes and calculate chi2 of this difference
   // if this chi2 is bigger than a given threshold, assume that the current cluster is
   // a kink an goto next padrow
    
    
    for (Int_t idelta = -kDelta; idelta <= kDelta; idelta++){
      // loop over irow +- kDelta rows (neighboured rows)
      // 
      // 
      if (idelta == 0) continue;
      if (idelta + irow < 0 || idelta + irow > 159) continue;   // don't go out of ROC
      AliTPCclusterMI * cluster = track->GetClusterPointer(irow + idelta);
      if (!cluster) continue;
      if (cluster->GetType() < 0) continue;
      if (cluster->GetDetector() != sector) continue;
      Double_t x = cluster->GetX() - xref;  // x = differece: current cluster - cluster @ irow
      nclFound++;
      if (idelta < 0){
	ncl0++;
	fitterY0.AddPoint(&x, cluster->GetY(), csigmaY);
	fitterZ0.AddPoint(&x, cluster->GetZ(), csigmaZ);
      }
      if (idelta > 0){
	ncl1++;
	fitterY1.AddPoint(&x, cluster->GetY(), csigmaY);
	fitterZ1.AddPoint(&x, cluster->GetZ(), csigmaZ);
      }
      fitterY2.AddPoint(&x, cluster->GetY(), csigmaY);  
      fitterZ2.AddPoint(&x, cluster->GetZ(), csigmaZ);  
    }  // loop over neighbourhood for fitter filling 
                                                    
    if (nclFound < kDelta * kMinRatio) continue;    // if not enough clusters found in neighbourhood goto next padrow
    fitterY2.Eval();
    fitterZ2.Eval();
    Double_t chi2 = (fitterY2.GetChisquare() + fitterZ2.GetChisquare()) / (2. * nclFound - 6.);
    if (chi2 > kCutChi2) continue;   // if chi^2 is too big goto next padrow
    
    // REMOVE KINK
    // only when there are enough clusters (4) in each direction
    if (ncl0 > 4){
      fitterY0.Eval();
      fitterZ0.Eval();
    }
    if (ncl1 > 4){
      fitterY1.Eval();
      fitterZ1.Eval();
    }
    
    if (ncl0 > 4 && ncl1 > 4){
      fitterY0.GetCovarianceMatrix(matrixY0);
      fitterY1.GetCovarianceMatrix(matrixY1);
      fitterZ0.GetCovarianceMatrix(matrixZ0);
      fitterZ1.GetCovarianceMatrix(matrixZ1);
      fitterY1.GetParameters(paramY1);
      fitterZ1.GetParameters(paramZ1);
      fitterY0.GetParameters(paramY0);
      fitterZ0.GetParameters(paramZ0);
      paramY0 -= paramY1;
      paramZ0 -= paramZ1;
      matrixY0 += matrixY1;
      matrixZ0 += matrixZ1;
      Double_t chi2 = 0;
      
      TMatrixD difY(2, 1, paramY0.GetMatrixArray());
      TMatrixD difYT(1, 2, paramY0.GetMatrixArray());
      matrixY0.Invert();
      TMatrixD mulY(matrixY0, TMatrixD::kMult, difY);
      TMatrixD chi2Y(difYT, TMatrixD::kMult, mulY);
      chi2 += chi2Y(0, 0);
      
      TMatrixD difZ(2, 1, paramZ0.GetMatrixArray());
      TMatrixD difZT(1, 2, paramZ0.GetMatrixArray());
      matrixZ0.Invert();
      TMatrixD mulZ(matrixZ0, TMatrixD::kMult, difZ);
      TMatrixD chi2Z(difZT, TMatrixD::kMult, mulZ);
      chi2 += chi2Z(0, 0);      
      
      // REMOVE KINK
      if (chi2 * 0.25 > kCutChi2) continue;   // if chi2 is too big goto next padrow
      // fit tracklet with polynom of 2nd order and two polynoms of 1st order
      // take both polynoms of 1st order, calculate difference of their parameters
      // add covariance matrixes and calculate chi2 of this difference
      // if this chi2 is bigger than a given threshold, assume that the current cluster is
      // a kink an goto next padrow
    }
    // current padrow has no kink
    
    // get fit parameters from pol2 fit: 
    Double_t paramY[4], paramZ[4];
    paramY[0] = fitterY2.GetParameter(0);
    paramY[1] = fitterY2.GetParameter(1);
    paramY[2] = fitterY2.GetParameter(2);
    paramZ[0] = fitterZ2.GetParameter(0);
    paramZ[1] = fitterZ2.GetParameter(1);
    paramZ[2] = fitterZ2.GetParameter(2);    
    
    Double_t tracky = paramY[0];
    Double_t trackz = paramZ[0];
    Float_t  deltay = tracky - cluster0->GetY();
    Float_t  deltaz = trackz - cluster0->GetZ();
    Float_t  angley = paramY[1] - paramY[0] / xref;
    Float_t  anglez = paramZ[1];
    
    Float_t max = cluster0->GetMax();
    TProfile *profAmpRow =  (TProfile*)fArrayAmpRow->At(sector);
    if ( irow >= kFirstLargePad) max /= kLargePadSize;
    profAmpRow->Fill( (Double_t)cluster0->GetRow(), max );
    TH1F *hisAmp =  (TH1F*)fArrayAmp->At(sector);
    hisAmp->Fill(max);
    
    Int_t ipad = 0;
    if (cluster0->GetDetector() >= 36) {
      ipad = 1;
      if (cluster0->GetRow() > 63) ipad = 2;
    }
    
    TH3F * his3 = 0;
    his3 = (TH3F*)fRMSY->At(ipad);
    if (his3) his3->Fill(250 - TMath::Abs(cluster0->GetZ()), TMath::Abs(angley), TMath::Sqrt(cluster0->GetSigmaY2()) );
    his3 = (TH3F*)fRMSZ->At(ipad);
    if (his3) his3->Fill( 250 - TMath::Abs(cluster0->GetZ()), TMath::Abs(anglez), TMath::Sqrt(cluster0->GetSigmaZ2()) );
    
    his3 = (TH3F*)fArrayQRMSY->At(GetBin(cluster0->GetMax(), ipad));
    if (his3) his3->Fill( 250 - TMath::Abs(cluster0->GetZ()), TMath::Abs(angley), TMath::Sqrt(cluster0->GetSigmaY2()) );
    his3 = (TH3F*)fArrayQRMSZ->At(GetBin(cluster0->GetMax(), ipad));
    if (his3) his3->Fill( 250 - TMath::Abs(cluster0->GetZ()), TMath::Abs(anglez), TMath::Sqrt(cluster0->GetSigmaZ2()) );

    // Fill resolution histograms
    Bool_t useForResol = kTRUE;
    if (fitterY2.GetParError(0) > kErrorFraction * csigmaY) useForResol = kFALSE;

    if (useForResol){
      fDeltaY->Fill(deltay);
      fDeltaZ->Fill(deltaz);
      his3 = (TH3F*)fResolY->At(ipad);
      if (his3) his3->Fill( 250 - TMath::Abs(cluster0->GetZ()), TMath::Abs(angley), deltay );
      his3 = (TH3F*)fResolZ->At(ipad);
      if (his3) his3->Fill( 250 - TMath::Abs(cluster0->GetZ()), TMath::Abs(anglez), deltaz );
      his3 = (TH3F*)fArrayQDY->At(GetBin(cluster0->GetMax(), ipad));
      if (his3) his3->Fill( 250 - TMath::Abs(cluster0->GetZ()),TMath::Abs(angley), deltay );
      his3 = (TH3F*)fArrayQDZ->At(GetBin(cluster0->GetMax(), ipad));
      if (his3) his3->Fill( 250 - TMath::Abs(cluster0->GetZ()),TMath::Abs(anglez), deltaz );
    }
    
    
    if (useForResol && nclFound > 2 * kMinRatio * kDelta){
      // 
      // fill resolution trees
      // 
      static TLinearFitter fitY0(3, "pol2");
      static TLinearFitter fitZ0(3, "pol2");
      static TLinearFitter fitY2(5, "hyp4");
      static TLinearFitter fitZ2(5, "hyp4");
      static TLinearFitter fitY2Q(5, "hyp4");
      static TLinearFitter fitZ2Q(5, "hyp4");
      static TLinearFitter fitY2S(5, "hyp4");
      static TLinearFitter fitZ2S(5, "hyp4");
      fitY0.ClearPoints();
      fitZ0.ClearPoints();
      fitY2.ClearPoints();
      fitZ2.ClearPoints();
      fitY2Q.ClearPoints();
      fitZ2Q.ClearPoints();
      fitY2S.ClearPoints();
      fitZ2S.ClearPoints();
      
      for (Int_t idelta = -kDelta; idelta <= kDelta; idelta++){
         // loop over irow +- kDelta rows (neighboured rows)
         // 
         // 
	if (idelta == 0) continue;
	if (idelta + irow < 0 || idelta + irow > 159) continue;   // don't go out of ROC
//	if (idelta + irow > 159) continue;
	AliTPCclusterMI * cluster = track->GetClusterPointer(irow + idelta);
	if (!cluster) continue;
	if (cluster->GetType() < 0) continue;
	if (cluster->GetDetector() != sector) continue;
	Double_t x = cluster->GetX() - xref;
	Double_t sigmaY0 = fClusterParam->GetError0Par( 0, ipad, (250.0 - TMath::Abs(cluster->GetZ())), TMath::Abs(angley) );
	Double_t sigmaZ0 = fClusterParam->GetError0Par( 1, ipad, (250.0 - TMath::Abs(cluster->GetZ())), TMath::Abs(anglez) );
	//
	Double_t sigmaYQ = fClusterParam->GetErrorQPar( 0, ipad, (250.0 - TMath::Abs(cluster->GetZ())), TMath::Abs(angley), TMath::Abs(cluster->GetMax()) );
	Double_t sigmaZQ = fClusterParam->GetErrorQPar( 1, ipad, (250.0 - TMath::Abs(cluster->GetZ())), TMath::Abs(anglez), TMath::Abs(cluster->GetMax()) );
	Double_t sigmaYS = fClusterParam->GetErrorQParScaled( 0, ipad, (250.0 - TMath::Abs(cluster->GetZ())), 
                                                           TMath::Abs(angley), TMath::Abs(cluster->GetMax()) );
	Double_t sigmaZS = fClusterParam->GetErrorQParScaled( 1, ipad, (250.0 - TMath::Abs(cluster->GetZ())), 
                                                           TMath::Abs(anglez), TMath::Abs(cluster->GetMax()) );
	Float_t rmsYFactor = fClusterParam->GetShapeFactor( 0, ipad,(250.0 - TMath::Abs(cluster->GetZ())),
							   TMath::Abs(anglez), TMath::Abs(cluster->GetMax()),
							   TMath::Sqrt(cluster0->GetSigmaY2()), 0 );
	Float_t rmsZFactor = fClusterParam->GetShapeFactor(0, ipad,(250.0 - TMath::Abs(cluster->GetZ())),
							   TMath::Abs(anglez), TMath::Abs(cluster->GetMax()),
							   TMath::Sqrt(cluster0->GetSigmaZ2()),0 );
	sigmaYS  = TMath::Sqrt(sigmaYS * sigmaYS + rmsYFactor * rmsYFactor / 12.);
	sigmaZS  = TMath::Sqrt(sigmaZS * sigmaZS + rmsZFactor * rmsZFactor / 12. + rmsYFactor * rmsYFactor / 24.);
	//
	if (kDelta != 0){
	  fitY0.AddPoint(&x, cluster->GetY(), sigmaY0);
	  fitZ0.AddPoint(&x, cluster->GetZ(), sigmaZ0);
	}
	Double_t xxx[4];
	xxx[0] = ( (idelta+irow) % 2 == 0 ) ? 1 : 0;
	xxx[1] = x;
	xxx[2] = ( (idelta+irow) % 2 == 0 ) ? x : 0;
	xxx[3] = x * x;	
	fitY2.AddPoint(xxx, cluster->GetY(), sigmaY0);
	fitY2Q.AddPoint(xxx, cluster->GetY(), sigmaYQ);
	fitY2S.AddPoint(xxx, cluster->GetY(), sigmaYS);
	fitZ2.AddPoint(xxx, cluster->GetZ(), sigmaZ0);
	fitZ2Q.AddPoint(xxx, cluster->GetZ(), sigmaZQ);
	fitZ2S.AddPoint(xxx, cluster->GetZ(), sigmaZS);
	//
      }  // neigbouhood-loop
      //
      fitY0.Eval();
      fitZ0.Eval();
      fitY2.Eval();
      fitZ2.Eval();
      fitY2Q.Eval();
      fitZ2Q.Eval();
      fitY2S.Eval();
      fitZ2S.Eval();
      Float_t chi2Y0 = fitY0.GetChisquare() / (nclFound-3.);
      Float_t chi2Z0 = fitZ0.GetChisquare() / (nclFound-3.);
      Float_t chi2Y2 = fitY2.GetChisquare() / (nclFound-5.);
      Float_t chi2Z2 = fitZ2.GetChisquare() / (nclFound-5.);
      Float_t chi2Y2Q = fitY2Q.GetChisquare() / (nclFound-5.);
      Float_t chi2Z2Q = fitZ2Q.GetChisquare() / (nclFound-5.);
      Float_t chi2Y2S = fitY2S.GetChisquare() / (nclFound-5.);
      Float_t chi2Z2S = fitZ2S.GetChisquare() / (nclFound-5.);
      //
      static  TVectorD    parY0(3);
      static  TMatrixD    matY0(3, 3);
      static  TVectorD    parZ0(3);
      static  TMatrixD    matZ0(3, 3);
      fitY0.GetParameters(parY0);
      fitY0.GetCovarianceMatrix(matY0);
      fitZ0.GetParameters(parZ0);
      fitZ0.GetCovarianceMatrix(matZ0);
      //
      static  TVectorD    parY2(5);
      static  TMatrixD    matY2(5,5);
      static  TVectorD    parZ2(5);
      static  TMatrixD    matZ2(5,5);
      fitY2.GetParameters(parY2);
      fitY2.GetCovarianceMatrix(matY2);
      fitZ2.GetParameters(parZ2);
      fitZ2.GetCovarianceMatrix(matZ2);
      //
      static  TVectorD    parY2Q(5);
      static  TMatrixD    matY2Q(5,5);
      static  TVectorD    parZ2Q(5);
      static  TMatrixD    matZ2Q(5,5);
      fitY2Q.GetParameters(parY2Q);
      fitY2Q.GetCovarianceMatrix(matY2Q);
      fitZ2Q.GetParameters(parZ2Q);
      fitZ2Q.GetCovarianceMatrix(matZ2Q);
      static  TVectorD    parY2S(5);
      static  TMatrixD    matY2S(5,5);
      static  TVectorD    parZ2S(5);
      static  TMatrixD    matZ2S(5,5);
      fitY2S.GetParameters(parY2S);
      fitY2S.GetCovarianceMatrix(matY2S);
      fitZ2S.GetParameters(parZ2S);
      fitZ2S.GetCovarianceMatrix(matZ2S);
      Float_t sigmaY0   = TMath::Sqrt(matY0(0,0));
      Float_t sigmaZ0   = TMath::Sqrt(matZ0(0,0));
      Float_t sigmaDY0  = TMath::Sqrt(matY0(1,1));
      Float_t sigmaDZ0  = TMath::Sqrt(matZ0(1,1));
      Float_t sigmaY2   = TMath::Sqrt(matY2(1,1));
      Float_t sigmaZ2   = TMath::Sqrt(matZ2(1,1));
      Float_t sigmaDY2  = TMath::Sqrt(matY2(3,3));
      Float_t sigmaDZ2  = TMath::Sqrt(matZ2(3,3));
      Float_t sigmaY2Q  = TMath::Sqrt(matY2Q(1,1));
      Float_t sigmaZ2Q  = TMath::Sqrt(matZ2Q(1,1));
      Float_t sigmaDY2Q = TMath::Sqrt(matY2Q(3,3));
      Float_t sigmaDZ2Q = TMath::Sqrt(matZ2Q(3,3));
      Float_t sigmaY2S  = TMath::Sqrt(matY2S(1,1));
      Float_t sigmaZ2S  = TMath::Sqrt(matZ2S(1,1));
      Float_t sigmaDY2S = TMath::Sqrt(matY2S(3,3));
      Float_t sigmaDZ2S = TMath::Sqrt(matZ2S(3,3));
      
      // Error parameters
      Float_t csigmaY0 = fClusterParam->GetError0Par(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),TMath::Abs(angley));
      Float_t csigmaZ0 = fClusterParam->GetError0Par(1,ipad,(250.0-TMath::Abs(cluster0->GetZ())),TMath::Abs(anglez));
      //
      Float_t csigmaYQ = fClusterParam->GetErrorQPar(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						     TMath::Abs(angley), TMath::Abs(cluster0->GetMax()));
      Float_t csigmaZQ = fClusterParam->GetErrorQPar(1,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						       TMath::Abs(anglez),TMath::Abs(cluster0->GetMax()));
      Float_t csigmaYS = fClusterParam->GetErrorQParScaled(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						     TMath::Abs(angley), TMath::Abs(cluster0->GetMax()));
      Float_t csigmaZS = fClusterParam->GetErrorQParScaled(1,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						       TMath::Abs(anglez),TMath::Abs(cluster0->GetMax()));
      
      // RMS parameters
      Float_t meanRMSY = 0;
      Float_t meanRMSZ = 0;
      Int_t   nclRMS = 0;
      for (Int_t idelta = -2; idelta <= 2; idelta++){
        // loop over neighbourhood
	if (idelta+irow < 0 || idelta+irow > 159) continue;
// 	if (idelta+irow>159) continue;
	AliTPCclusterMI * cluster = track->GetClusterPointer(irow+idelta);
	if (!cluster) continue;
	meanRMSY += TMath::Sqrt(cluster->GetSigmaY2());
	meanRMSZ += TMath::Sqrt(cluster->GetSigmaZ2());
	nclRMS++;
      }
      meanRMSY /= nclRMS; 
      meanRMSZ /= nclRMS; 

      Float_t rmsY      = TMath::Sqrt(cluster0->GetSigmaY2());  
      Float_t rmsZ      = TMath::Sqrt(cluster0->GetSigmaZ2());
      Float_t rmsYT     = fClusterParam->GetRMSQ(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						TMath::Abs(angley), TMath::Abs(cluster0->GetMax()));
      Float_t rmsZT     = fClusterParam->GetRMSQ(1,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						TMath::Abs(anglez), TMath::Abs(cluster0->GetMax()));
      Float_t rmsYT0    = fClusterParam->GetRMS0(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						 TMath::Abs(angley));
      Float_t rmsZT0    = fClusterParam->GetRMS0(1,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						 TMath::Abs(anglez));
      Float_t rmsYSigma = fClusterParam->GetRMSSigma(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						     TMath::Abs(anglez), TMath::Abs(cluster0->GetMax()));
      Float_t rmsZSigma = fClusterParam->GetRMSSigma(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						     TMath::Abs(anglez), TMath::Abs(cluster0->GetMax()));
      Float_t rmsYFactor = fClusterParam->GetShapeFactor(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
							 TMath::Abs(anglez), TMath::Abs(cluster0->GetMax()),
							 rmsY,meanRMSY);
      Float_t rmsZFactor = fClusterParam->GetShapeFactor(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
							 TMath::Abs(anglez), TMath::Abs(cluster0->GetMax()),
							 rmsZ,meanRMSZ);
      
      // cluster debug
      (*fDebugStream)<<"ResolCl"<<	// valgrind 3   40,000 bytes in 1 blocks are possibly 
      "Sector="<<sector<<
      "Cl.="<<cluster0<<
      "CSigmaY0="<<csigmaY0<<   // cluster errorY
      "CSigmaYQ="<<csigmaYQ<<   // cluster errorY - q scaled
      "CSigmaYS="<<csigmaYS<<   // cluster errorY - q scaled
      "CSigmaZ0="<<csigmaZ0<<   // 
      "CSigmaZQ="<<csigmaZQ<<
      "CSigmaZS="<<csigmaZS<<
      "shapeYF="<<rmsYFactor<<
      "shapeZF="<<rmsZFactor<<
      "rmsY="<<rmsY<<
      "rmsZ="<<rmsZ<<
      "rmsYM="<<meanRMSY<<
      "rmsZM="<<meanRMSZ<<
      "rmsYT="<<rmsYT<<
      "rmsZT="<<rmsZT<<
      "rmsYT0="<<rmsYT0<<
      "rmsZT0="<<rmsZT0<<
      "rmsYS="<<rmsYSigma<<  
      "rmsZS="<<rmsZSigma<<
      "IPad="<<ipad<<
      "Ncl="<<nclFound<<	
      "PY0.="<<&parY0<<
      "PZ0.="<<&parZ0<<
      "SigmaY0="<<sigmaY0<< 
      "SigmaZ0="<<sigmaZ0<< 
      "angley="<<angley<<
      "anglez="<<anglez<<
      

      "\n";

//       tracklet dubug
      (*fDebugStream)<<"ResolTr"<<	
      "IPad="<<ipad<<
      "Sector="<<sector<<
      "Ncl="<<nclFound<<	
      "chi2Y0="<<chi2Y0<<
      "chi2Z0="<<chi2Z0<<
      "chi2Y2="<<chi2Y2<<
      "chi2Z2="<<chi2Z2<<
      "chi2Y2Q="<<chi2Y2Q<<
      "chi2Z2Q="<<chi2Z2Q<<
      "chi2Y2S="<<chi2Y2S<<
      "chi2Z2S="<<chi2Z2S<<
      "PY0.="<<&parY0<<
      "PZ0.="<<&parZ0<<
      "PY2.="<<&parY2<<
      "PZ2.="<<&parZ2<<
      "PY2Q.="<<&parY2Q<<
      "PZ2Q.="<<&parZ2Q<<
      "PY2S.="<<&parY2S<<
      "PZ2S.="<<&parZ2S<<
      "SigmaY0="<<sigmaY0<< 
      "SigmaZ0="<<sigmaZ0<< 
      "SigmaDY0="<<sigmaDY0<< 
      "SigmaDZ0="<<sigmaDZ0<< 
      "SigmaY2="<<sigmaY2<< 
      "SigmaZ2="<<sigmaZ2<< 
      "SigmaDY2="<<sigmaDY2<< 
      "SigmaDZ2="<<sigmaDZ2<< 
      "SigmaY2Q="<<sigmaY2Q<< 
      "SigmaZ2Q="<<sigmaZ2Q<< 
      "SigmaDY2Q="<<sigmaDY2Q<< 
      "SigmaDZ2Q="<<sigmaDZ2Q<< 
      "SigmaY2S="<<sigmaY2S<< 
      "SigmaZ2S="<<sigmaZ2S<< 
      "SigmaDY2S="<<sigmaDY2S<< 
      "SigmaDZ2S="<<sigmaDZ2S<< 
	"angley="<<angley<<
	"anglez="<<anglez<<
      

      "\n";
    }  // if (useForResol && nclFound > 2 * kMinRatio * kDelta)
  }    // loop over all padrows along the track: for (Int_t irow = 0; irow < 159; irow++)
}  // FillResolutionHistoLocal(...)


void  AliTPCcalibTracks::AlignUpDown(AliTPCseed * track, AliESDtrack * esdTrack){
  //
  // Make simple parabolic fit
  // Write debug information to 'TPCSelectorDebug.root'
  //
  const Int_t kMinClusters = 60;
  const Int_t kMinClustersSector = 15;
  const Float_t kSigmaCut = 6;
  const Float_t kMaxTan = TMath::Tan(TMath::Pi() * 10. / 180.);
  const Float_t kDeadZone = 6.;
  const Float_t kMinZ     = 15;
  if (track->GetNumberOfClusters() < kMinClusters) return;
  if (TMath::Abs(track->GetZ()) < kMinZ) return;
  //
  Int_t nclUp   = 0;
  Int_t nclDown = 0;
  Int_t rSector =-1;
  Float_t refX  = (param.GetInnerRadiusUp() + param.GetOuterRadiusLow()) * 0.5;
  
  for (Int_t irow = 0; irow < 159; irow++){
    AliTPCclusterMI * cluster0 = track->GetClusterPointer(irow);
    if (!cluster0) continue;
    Int_t sector = cluster0->GetDetector();
    if (rSector < 0) rSector = sector % 36;
    if (sector % 36 != rSector) continue;
    if ( ((TMath::Abs(cluster0->GetY()) - kDeadZone) / cluster0->GetX()) > kMaxTan ) continue;  //remove edge clusters
    if (sector > 35) nclUp++;
    if (sector < 36) nclDown++;
  }  // loop over padrows
  if (nclUp < kMinClustersSector) return;
  if (nclDown < kMinClustersSector) return;
  
  TLinearFitter fitterY(5,"hyp4");  //fitter with common 2 nd derivation
  TLinearFitter fitterZ(5,"hyp4");
  TLinearFitter fitterY0(3,"pol2");   // valgrind 3  58,142 bytes in 2,117 blocks are indirectly lost
                                      // valgrind    608,151 (198,860 direct, 409,291 indirect) bytes in 1,108 blocks are definitely lost
  TLinearFitter fitterZ0(3,"pol2");
  TLinearFitter fitterY1(3,"pol2");   // valgrind   4,284 bytes in 21 blocks are possibly lost 
  TLinearFitter fitterZ1(3,"pol2");   // valgrind   6 blocks possibly lost   57,956 bytes in 844 blocks are indirectly lost
  
  Float_t msigmay = 1;
  Float_t msigmaz = 1;
  Float_t param0[3];
  Float_t param1[3];
  Float_t angley = 0;
  Float_t anglez = 0;
  
  for (Int_t iter = 0; iter < 3; iter++){
    nclUp  = 0;
    nclDown= 0;
    for (Int_t irow = 0; irow < 159; irow++){
      AliTPCclusterMI * cluster0 = track->GetClusterPointer(irow);
      if (!cluster0) continue;
      Int_t sector = cluster0->GetDetector();
      if (sector % 36 != rSector) continue;
      Double_t y = cluster0->GetY();
      Double_t z = cluster0->GetZ();
      //remove edge clusters
      if ( (iter == 0) && ((TMath::Abs(cluster0->GetY()) - kDeadZone) / cluster0->GetX()) > kMaxTan ) continue;  
      if (iter > 0){
	Float_t tx = cluster0->GetX() - refX;
	Float_t ty = 0;
	if (sector < 36){
	  ty = param0[0] + param0[1] * tx + param0[2] * tx * tx;
	}else{
	  ty = param1[0] + param1[1] * tx + param1[2] * tx * tx;	  
	}
	if (((TMath::Abs(ty)-kDeadZone)/cluster0->GetX())>kMaxTan) continue;
	if (TMath::Abs(ty-y)>kSigmaCut*(msigmay+0.2)) continue;
      }
      Int_t  ipad = 0;
      if (cluster0->GetDetector() >= 36) {
	ipad = 1;
	if (cluster0->GetRow()>63) ipad=2;
      }
      //
      Float_t sigmaY = msigmay;
      Float_t sigmaZ = msigmay;      
      if (iter == 2){
	sigmaY = fClusterParam->GetErrorQParScaled(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
							   TMath::Abs(angley), TMath::Abs(cluster0->GetMax()));
	sigmaZ = fClusterParam->GetErrorQParScaled(1,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
							   TMath::Abs(anglez),TMath::Abs(cluster0->GetMax()));
      }  // iter == 2
      Double_t deltaX = cluster0->GetX() - refX;
      Double_t x[5];
      x[0] = (ipad==0) ? 0:1;
      x[1] = deltaX;
      x[2] = (ipad==0) ? 0:deltaX;
      x[3] = deltaX*deltaX;
      if (ipad < 2){
	fitterY.AddPoint(x,y,sigmaY);
	fitterZ.AddPoint(x,z,sigmaZ);
      }
      if (ipad == 0){
	nclDown++;
	fitterY0.AddPoint(&deltaX,y,sigmaY);
	fitterZ0.AddPoint(&deltaX,z,sigmaZ);
      }
      if (ipad == 1){
	nclUp++;
	fitterY1.AddPoint(&deltaX,y,sigmaY);
	fitterZ1.AddPoint(&deltaX,z,sigmaZ);
      }
    }    // loop over padrows
    if (nclUp < kMinClustersSector) continue;
    if (nclDown < kMinClustersSector) continue;
    fitterY.Eval();
    fitterZ.Eval();
    fitterY0.Eval();
    fitterZ0.Eval();
    fitterY1.Eval();
    fitterZ1.Eval();
    param0[0] = fitterY0.GetParameter(0);
    param0[1] = fitterY0.GetParameter(1);
    param0[2] = fitterY0.GetParameter(2);
    param1[0] = fitterY1.GetParameter(0);
    param1[1] = fitterY1.GetParameter(1);
    param1[2] = fitterY1.GetParameter(2);
    //
    angley = fitterY.GetParameter(2);
    anglez = fitterZ.GetParameter(2);
    //
    TVectorD    parY(5);
    TMatrixD    matY(5,5);
    TVectorD    parZ(5);
    TMatrixD    matZ(5,5);
    Double_t    chi2Y = fitterY.GetChisquare() / (nclUp+nclDown); 
    Double_t    chi2Z = fitterZ.GetChisquare() / (nclUp+nclDown); 
    fitterY.GetParameters(parY);
    fitterY.GetCovarianceMatrix(matY);
    fitterZ.GetParameters(parZ);
    fitterZ.GetCovarianceMatrix(matZ); 
    if (iter == 0) {
      msigmay = msigmay*TMath::Sqrt(chi2Y);
      msigmaz = msigmaz*TMath::Sqrt(chi2Z);
    }
    Float_t sigmaY  = TMath::Sqrt(matY(1,1)*chi2Y);
    Float_t sigmaDY = TMath::Sqrt(matY(3,3)*chi2Y);
    Float_t sigmaDDY = TMath::Sqrt(matY(4,4)*chi2Y);
    Float_t sigmaZ  = TMath::Sqrt(matZ(1,1)*chi2Z);
    Float_t sigmaDZ = TMath::Sqrt(matZ(3,3)*chi2Z);
    Float_t sigmaDDZ = TMath::Sqrt(matZ(4,4)*chi2Z);
    //
    TVectorD    parY0(3);
    TMatrixD    matY0(3,3);
    TVectorD    parZ0(3);
    TMatrixD    matZ0(3,3);
    Double_t    chi2Y0= fitterY0.GetChisquare()/(nclDown); 
    Double_t    chi2Z0= fitterZ0.GetChisquare()/(nclDown); 
    fitterY0.GetParameters(parY0);
    fitterY0.GetCovarianceMatrix(matY0);
    fitterZ0.GetParameters(parZ0);
    fitterZ0.GetCovarianceMatrix(matZ0); 
    Float_t sigmaY0  = TMath::Sqrt(matY0(0,0)*chi2Y0);
    Float_t sigmaDY0 = TMath::Sqrt(matY0(1,1)*chi2Y0);
    Float_t sigmaDDY0 = TMath::Sqrt(matY0(2,2)*chi2Y0);
    Float_t sigmaZ0  = TMath::Sqrt(matZ0(0,0)*chi2Z0);
    Float_t sigmaDZ0 = TMath::Sqrt(matZ0(1,1)*chi2Z0);
    Float_t sigmaDDZ0 = TMath::Sqrt(matZ0(2,2)*chi2Z0);
    //
    TVectorD    parY1(3);
    TMatrixD    matY1(3,3);
    TVectorD    parZ1(3);
    TMatrixD    matZ1(3,3);
    Double_t    chi2Y1= fitterY1.GetChisquare()/(nclUp); 
    Double_t    chi2Z1= fitterZ1.GetChisquare()/(nclUp); 
    fitterY1.GetParameters(parY1);
    fitterY1.GetCovarianceMatrix(matY1);
    fitterZ1.GetParameters(parZ1);
    fitterZ1.GetCovarianceMatrix(matZ1); 
    Float_t sigmaY1  = TMath::Sqrt(matY1(0,0)*chi2Y1);
    Float_t sigmaDY1 = TMath::Sqrt(matY1(1,1)*chi2Y1);
    Float_t sigmaDDY1 = TMath::Sqrt(matY1(2,2)*chi2Y1);
    Float_t sigmaZ1  = TMath::Sqrt(matZ1(0,0)*chi2Z1);
    Float_t sigmaDZ1 = TMath::Sqrt(matZ1(1,1)*chi2Z1);
    Float_t sigmaDDZ1 = TMath::Sqrt(matZ1(2,2)*chi2Z1);
    const AliESDfriendTrack * ftrack = esdTrack->GetFriendTrack();
    AliTrackPointArray *points = (AliTrackPointArray*)ftrack->GetTrackPointArray();
    
    if (iter>0) (*fDebugStream)<<"Align"<<   // valgrind    85,932 bytes in 543 blocks are still reachable
      "track.="<<track<<
      "Iter="<<iter<<
      "xref="<<refX<<
      "Points="<<points<<
      "Sector="<<rSector<<
      "nclUp="<<nclUp<<
      "nclDown="<<nclDown<<
      "angley="<<angley<<
      "anglez="<<anglez<<
      
      "chi2Y="<<chi2Y<<
      "chi2Z="<<chi2Z<<
      "parY.="<<&parY<<
      "parZ.="<<&parZ<<
      "matY.="<<&matY<<
      "matZ.="<<&matZ<<
      "sigmaY="<<sigmaY<<
      "sigmaZ="<<sigmaZ<<
      "sigmaDY="<<sigmaDY<<
      "sigmaDZ="<<sigmaDZ<<
      "sigmaDDY="<<sigmaDDY<<
      "sigmaDDZ="<<sigmaDDZ<<
      
      "chi2Y0="<<chi2Y0<<
      "chi2Z0="<<chi2Z0<<
      "parY0.="<<&parY0<<
      "parZ0.="<<&parZ0<<
      "matY0.="<<&matY0<<
      "matZ0.="<<&matZ0<<
      "sigmaY0="<<sigmaY0<<
      "sigmaZ0="<<sigmaZ0<<
      "sigmaDY0="<<sigmaDY0<<
      "sigmaDZ0="<<sigmaDZ0<<
      "sigmaDDY0="<<sigmaDDY0<<
      "sigmaDDZ0="<<sigmaDDZ0<<
      
      "chi2Y1="<<chi2Y1<<
      "chi2Z1="<<chi2Z1<<
      "parY1.="<<&parY1<<
      "parZ1.="<<&parZ1<<
      "matY1.="<<&matY1<<
      "matZ1.="<<&matZ1<<
      "sigmaY1="<<sigmaY1<<
      "sigmaZ1="<<sigmaZ1<<
      "sigmaDY1="<<sigmaDY1<<
      "sigmaDZ1="<<sigmaDZ1<<
      "sigmaDDY1="<<sigmaDDY1<<
      "sigmaDDZ1="<<sigmaDDZ1<<
      "\n";
  }      // for (Int_t iter = 0; iter < 3; iter++)
}


TH2D * AliTPCcalibTracks::MakeDiff(TH2D * hfit, TF2 * func){
   // 
   // creates a new histogram which contains the difference between
   // the histogram hfit and the function func
   // 
  TH2D * result = (TH2D*)hfit->Clone();      // valgrind 3   40,139 bytes in 11 blocks are still reachable
  result->SetTitle(Form("%s fit residuals",result->GetTitle()));
  result->SetName(Form("%s fit residuals",result->GetName()));
  TAxis *xaxis = hfit->GetXaxis();
  TAxis *yaxis = hfit->GetYaxis();
  Double_t x[2];
  for (Int_t biny = 0; biny <= yaxis->GetNbins(); biny++) {
    x[1]  = yaxis->GetBinCenter(biny);
    for (Int_t binx = 0; binx <= xaxis->GetNbins(); binx++) {
      x[0]  = xaxis->GetBinCenter(binx);
      Int_t bin = hfit->GetBin(binx, biny);
      Double_t val = hfit->GetBinContent(bin);
//      result->SetBinContent( bin, (val - func->Eval(x[0], x[1])) / func->Eval(x[0], x[1]) );
      result->SetBinContent( bin, (val / func->Eval(x[0], x[1])) - 1 );
    }    
  }
  return result;
}



void  AliTPCcalibTracks::SetStyle(){
   // 
   // set style, can be called by all draw functions
   // 
   gROOT->SetStyle("Plain");
   gStyle->SetFillColor(10);
   gStyle->SetPadColor(10);
   gStyle->SetCanvasColor(10);
   gStyle->SetStatColor(10);
   gStyle->SetPalette(1,0);
   gStyle->SetNumberContours(60);
}


void AliTPCcalibTracks::Draw(Option_t* opt){
   // 
   // draw-function of AliTPCcalibTracks
   // will draws some exemplaric pictures
   // 
   
   cout << "***** not yet implemented *****" << endl;
   fDeltaY->Draw(opt);

}


void AliTPCcalibTracks::MakeReport(Int_t stat, char* pathName){ 
   // 
   // all functions are called, that produce pictures
   // the histograms are written to the directory 'pathName'
   // 'stat' is a threshhold: only histograms with more than 'stat' entries are wirtten to file
   // 'stat' is also the number of minEntries for MakeResPlotsQTree
   // 


   MakeAmpPlots(stat, pathName);
   MakeDeltaPlots(pathName);
   FitResolutionNew(pathName);
   FitRMSNew(pathName);
//    cout << "now comes MakeResPlotsQ" << endl;
//    MakeResPlotsQ(1, 1); 
//    cout << "now comes MakeResPlotsQTree" << endl;
   MakeResPlotsQTree(stat, pathName);

}
   

void AliTPCcalibTracks::MakeAmpPlots(Int_t stat, char* pathName){ 
   // 
   // creates several plots:
   // fArrayAmp.ps, fArrayAmpRow.ps and DeltaYZ.ps
   // fArrayAmp.ps: one histogram per sector, the histogram shows the charge per cluster
   // fArrayAmpRow.ps: one histogram per sector, mean max. amplitude vs. pad row with landau fit
   // DeltaYZ.ps: DeltaY and DeltaZ histogram with gaus fit
   // Empty histograms (sectors without data) are not written to file
   // the ps-files are written to the directory 'pathName', that is created if it does not exist
   // 'stat': only histograms with more than 'stat' entries are written to file.
   // 
   
   SetStyle();
   gSystem->MakeDirectory(pathName);
   gSystem->ChangeDirectory(pathName);
   
   TCanvas* c1 = new TCanvas();     // valgrind 3 ???  634 bytes in 28 blocks are still reachable
   TPostScript *ps; 
   // histograms with accumulated amplitude for all IROCs and OROCs
   TH1F *allAmpHisIROC = ((TH1F*)(fArrayAmp->At(0))->Clone());
   allAmpHisIROC->SetName("Amp all IROCs");
   allAmpHisIROC->SetTitle("Amp all IROCs");
   TH1F *allAmpHisOROC = ((TH1F*)(fArrayAmp->At(36))->Clone());
   allAmpHisOROC->SetName("Amp all OROCs");
   allAmpHisOROC->SetTitle("Amp all OROCs");
   
   
   ps = new TPostScript("fArrayAmp.ps", 112);
   cout << "creating fArrayAmp.ps..." << endl;
   for (Int_t i = 0; i < fArrayAmp->GetEntriesFast(); i++){
      if ( ((TH1F*)fArrayAmp->At(i))->GetEntries() < stat  ) continue;
      ps->NewPage();
      ((TH1F*)fArrayAmp->At(i))->Draw();
      c1->Update();              // valgrind 3
      if (i > 0 && i < 36) {
         allAmpHisIROC->Add(((TH1F*)fArrayAmp->At(i)));
         allAmpHisOROC->Add(((TH1F*)fArrayAmp->At(i+36)));
      }
   }
   ps->NewPage();
   allAmpHisIROC->Draw();
   c1->Update();              // valgrind
   ps->NewPage();
   allAmpHisOROC->Draw();
   c1->Update();
   ps->Close();
   delete ps;
   
   TH1F *his = 0;
   Double_t min = 0;
   Double_t max = 0;
   ps = new TPostScript("fArrayAmpRow.ps", 112);
   cout << "creating fArrayAmpRow.ps..." << endl;
   for (Int_t i = 0; i < fArrayAmpRow->GetEntriesFast(); i++){
      his = (TH1F*)fArrayAmpRow->At(i);
      if (his->GetEntries() < stat) continue;
      ps->NewPage();
      min = TMath::Max( his->GetBinCenter(his->GetMaximumBin() )-100., 0.);
      max = his->GetBinCenter(5*his->GetMaximumBin()) + 100;
      his->SetAxisRange(min, max);
      his->Fit("pol3", "q", "", min, max);
      // his->Draw("error");    // don't use this line when you don't want to have empty pages in the ps-file
      c1->Update();
   }
   ps->Close();
   delete ps;

   delete c1;
   gSystem->ChangeDirectory("..");
}


void AliTPCcalibTracks::MakeDeltaPlots(char* pathName){
   // 
   // creates several plots:
   // DeltaYZ.ps: DeltaY and DeltaZ histogram with gaus fit
   // the ps-files are written to the directory 'pathName', that is created if it does not exist
   // 
   
   SetStyle();
   gSystem->MakeDirectory(pathName);
   gSystem->ChangeDirectory(pathName);
   
   TCanvas* c1 = new TCanvas();     // valgrind 3 ???  634 bytes in 28 blocks are still reachable
   TPostScript *ps; 
   Double_t min = 0;
   Double_t max = 0;
   
   ps = new TPostScript("DeltaYZ.ps", 112);
   cout << "creating DeltaYZ.ps..." << endl;
   min = fDeltaY->GetBinCenter(fDeltaY->GetMaximumBin())-20;
   max = fDeltaY->GetBinCenter(fDeltaY->GetMaximumBin())+20;
   fDeltaY->SetAxisRange(min, max);
   ps->NewPage();
//    fDeltaY->Draw();
   fDeltaY->Fit("gaus","q","",min, max);        // valgrind 3  7 block possibly lost   2,400 bytes in 1 blocks are still reachable
   c1->Update();
   ps->NewPage();
   max = fDeltaZ->GetBinCenter(fDeltaZ->GetMaximumBin())+20;
   min = fDeltaZ->GetBinCenter(fDeltaZ->GetMaximumBin())-20;
   fDeltaZ->SetAxisRange(min, max);
//    fDeltaZ->Draw();
   fDeltaZ->Fit("gaus","q","",min, max);
   c1->Update();
   ps->Close();
   delete ps;
   delete c1;

   gSystem->ChangeDirectory("..");
}


void AliTPCcalibTracks::FitResolutionNew(char* pathName){
   // 
   // calculates different resulution fits in Y and Z direction
   // the histograms are written to 'ResolutionYZ.ps'
   // writes calculated resolution to 'resol.txt'
   // all files are stored in the directory pathName
   // 
   
   SetStyle();
   gSystem->MakeDirectory(pathName);
   gSystem->ChangeDirectory(pathName);
   
   TCanvas c;
   c.Divide(2,1); 
   cout << "creating ResolutionYZ.ps..." << endl;
   TPostScript *ps = new TPostScript("ResolutionYZ.ps", 112); 
   TF2 *fres = new TF2("fres","TMath::Sqrt([0]*[0]+[1]*[1]*x+[2]*[2]*y*y)",0,250,0,1);
   fres->SetParameter(0,0.02);
   fres->SetParameter(1,0.0054);
   fres->SetParameter(2,0.13);  
   
   TH1::AddDirectory(kTRUE);  // TH3F::FitSlicesZ() writes histograms into the current directory
   
   // create histogramw for Y-resolution
   TH3F * hisResY0 = (TH3F*)fResolY->At(0);
   hisResY0->FitSlicesZ();
   TH2D * hisResY0_2 = (TH2D*)gDirectory->Get("Resol Y0_2");
   TH3F * hisResY1 = (TH3F*)fResolY->At(1); 
   hisResY1->FitSlicesZ();
   TH2D * hisResY1_2 = (TH2D*)gDirectory->Get("Resol Y1_2");
   TH3F * hisResY2 = (TH3F*)fResolY->At(2);
   hisResY2->FitSlicesZ();
   TH2D * hisResY2_2 = (TH2D*)gDirectory->Get("Resol Y2_2");
    //
   ps->NewPage();
   c.cd(1);
   hisResY0_2->Fit(fres, "q");      // valgrind    132,072 bytes in 6 blocks are indirectly lost
   hisResY0_2->Draw("surf1"); 
   c.cd(2);
   MakeDiff(hisResY0_2,fres)->Draw("surf1");
   c.Update();
   //   c.SaveAs("ResolutionYPad0.eps");
   ps->NewPage();
   c.cd(1);
   hisResY1_2->Fit(fres, "q");
   hisResY1_2->Draw("surf1");
   c.cd(2);
   MakeDiff(hisResY1_2,fres)->Draw("surf1");
   c.Update();
   //   c.SaveAs("ResolutionYPad1.eps");
   ps->NewPage();
   c.cd(1);
   hisResY2_2->Fit(fres, "q");
   hisResY2_2->Draw("surf1"); 
   c.cd(2);
   MakeDiff(hisResY2_2,fres)->Draw("surf1");
   c.Update();
//    c.SaveAs("ResolutionYPad2.eps");
   
   // create histogramw for Z-resolution
   TH3F * hisResZ0 = (TH3F*)fResolZ->At(0);
   hisResZ0->FitSlicesZ();
   TH2D * hisResZ0_2 = (TH2D*)gDirectory->Get("Resol Z0_2");
   TH3F * hisResZ1 = (TH3F*)fResolZ->At(1);
   hisResZ1->FitSlicesZ();
   TH2D * hisResZ1_2 = (TH2D*)gDirectory->Get("Resol Z1_2");
   TH3F * hisResZ2 = (TH3F*)fResolZ->At(2);
   hisResZ2->FitSlicesZ();
   TH2D * hisResZ2_2 = (TH2D*)gDirectory->Get("Resol Z2_2");
   
   ps->NewPage();
   c.cd(1);
   hisResZ0_2->Fit(fres, "q");
   hisResZ0_2->Draw("surf1");
   c.cd(2);
   MakeDiff(hisResZ0_2,fres)->Draw("surf1");
   c.Update();
//    c.SaveAs("ResolutionZPad0.eps");
   ps->NewPage();
   c.cd(1);
   hisResZ1_2->Fit(fres, "q");
   hisResZ1_2->Draw("surf1"); 
   c.cd(2);
   MakeDiff(hisResZ1_2,fres)->Draw("surf1");
   c.Update();
//    c.SaveAs("ResolutionZPad1.eps");
   ps->NewPage();
   c.cd(1);
   hisResZ2_2->Fit(fres, "q");
   hisResZ2_2->Draw("surf1");  
   c.cd(2);
   MakeDiff(hisResZ2_2,fres)->Draw("surf1");
   c.Update();
//    c.SaveAs("ResolutionZPad2.eps");
   ps->Close();
   delete ps;
   
   // write calculated resoltuions to 'resol.txt'
   ofstream fresol("resol.txt");
   fresol<<"Pad 0.75 cm"<<"\n";
   hisResY0_2->Fit(fres, "q");                     // valgrind
   fresol<<"Y\t"<<fres->GetParameter(0)<<"\t"<<fres->GetParameter(1)<<"\t"<<fres->GetParameter(2)<<"\n";
   hisResZ0_2->Fit(fres, "q");
   fresol<<"Z\t"<<fres->GetParameter(0)<<"\t"<<fres->GetParameter(1)<<"\t"<<fres->GetParameter(2)<<"\n";
   //
   fresol<<"Pad 1.00 cm"<<1<<"\n";
   hisResY1_2->Fit(fres, "q");                     // valgrind
   fresol<<"Y\t"<<fres->GetParameter(0)<<"\t"<<fres->GetParameter(1)<<"\t"<<fres->GetParameter(2)<<"\n";
   hisResZ1_2->Fit(fres, "q");
   fresol<<"Z\t"<<fres->GetParameter(0)<<"\t"<<fres->GetParameter(1)<<"\t"<<fres->GetParameter(2)<<"\n";
   //
   fresol<<"Pad 1.50 cm"<<0<<"\n";
   hisResY2_2->Fit(fres, "q");
   fresol<<"Y\t"<<fres->GetParameter(0)<<"\t"<<fres->GetParameter(1)<<"\t"<<fres->GetParameter(2)<<"\n";
   hisResZ2_2->Fit(fres, "q");
   fresol<<"Z\t"<<fres->GetParameter(0)<<"\t"<<fres->GetParameter(1)<<"\t"<<fres->GetParameter(2)<<"\n";
   
   TH1::AddDirectory(kFALSE);
   gSystem->ChangeDirectory("..");
   delete fres;
}


void AliTPCcalibTracks::FitRMSNew(char* pathName){
   // 
   // calculates different resulution-rms fits in Y and Z direction
   // the histograms are written to 'RMS_YZ.ps'
   // writes calculated resolution rms to 'rms.txt'
   // all files are stored in the directory pathName
   // 
   
   SetStyle();
   gSystem->MakeDirectory(pathName);
   gSystem->ChangeDirectory(pathName);
   
   TCanvas c;        // valgrind 3   42,120 bytes in 405 blocks are still reachable   23,816 bytes in 229 blocks are still reachable
   c.Divide(2,1); 
   cout << "creating RMS_YZ.ps..." << endl;
   TPostScript *ps = new TPostScript("RMS_YZ.ps", 112); 
   TF2 *frms = new TF2("fres","TMath::Sqrt([0]*[0]+[1]*[1]*x+[2]*[2]*y*y)",0,250,0,1);
   frms->SetParameter(0,0.02);
   frms->SetParameter(1,0.0054);
   frms->SetParameter(2,0.13);  
   
   TH1::AddDirectory(kTRUE);  // TH3F::FitSlicesZ() writes histograms into the current directory
   
   // create histogramw for Y-RMS   
   TH3F * hisResY0 = (TH3F*)fRMSY->At(0);
   hisResY0->FitSlicesZ();
   TH2D * hisResY0_2 = (TH2D*)gDirectory->Get("RMS Y0_1");
   TH3F * hisResY1 = (TH3F*)fRMSY->At(1);
   hisResY1->FitSlicesZ();
   TH2D * hisResY1_2 = (TH2D*)gDirectory->Get("RMS Y1_1");
   TH3F * hisResY2 = (TH3F*)fRMSY->At(2);
   hisResY2->FitSlicesZ();
   TH2D * hisResY2_2 = (TH2D*)gDirectory->Get("RMS Y2_1");
   //
   ps->NewPage();
   c.cd(1);
   hisResY0_2->Fit(frms, "qn0"); 
   hisResY0_2->Draw("surf1"); 
   c.cd(2);
   MakeDiff(hisResY0_2,frms)->Draw("surf1");
   c.Update();
//    c.SaveAs("RMSYPad0.eps");
   ps->NewPage();
   c.cd(1);
   hisResY1_2->Fit(frms, "qn0");               // valgrind   several blocks possibly lost
   hisResY1_2->Draw("surf1");
   c.cd(2);
   MakeDiff(hisResY1_2,frms)->Draw("surf1");
   c.Update();
//    c.SaveAs("RMSYPad1.eps");
   ps->NewPage();
   c.cd(1);
   hisResY2_2->Fit(frms, "qn0");
   hisResY2_2->Draw("surf1"); 
   c.cd(2);
   MakeDiff(hisResY2_2,frms)->Draw("surf1");
   c.Update();
//    c.SaveAs("RMSYPad2.eps");
   
   // create histogramw for Z-RMS   
   TH3F * hisResZ0 = (TH3F*)fRMSZ->At(0);
   hisResZ0->FitSlicesZ();
   TH2D * hisResZ0_2 = (TH2D*)gDirectory->Get("RMS Z0_1");
   TH3F * hisResZ1 = (TH3F*)fRMSZ->At(1); 
   hisResZ1->FitSlicesZ();
   TH2D * hisResZ1_2 = (TH2D*)gDirectory->Get("RMS Z1_1");
   TH3F * hisResZ2 = (TH3F*)fRMSZ->At(2); 
   hisResZ2->FitSlicesZ();
   TH2D * hisResZ2_2 = (TH2D*)gDirectory->Get("RMS Z2_1");
   //
   ps->NewPage();
   c.cd(1);
   hisResZ0_2->Fit(frms, "qn0");         // valgrind
   hisResZ0_2->Draw("surf1");
   c.cd(2);
   MakeDiff(hisResZ0_2,frms)->Draw("surf1");
   c.Update();
//    c.SaveAs("RMSZPad0.eps");
   ps->NewPage();
   c.cd(1);
   hisResZ1_2->Fit(frms, "qn0");
   hisResZ1_2->Draw("surf1"); 
   c.cd(2);
   MakeDiff(hisResZ1_2,frms)->Draw("surf1");
   c.Update();
//    c.SaveAs("RMSZPad1.eps");
   ps->NewPage();
   c.cd(1);
   hisResZ2_2->Fit(frms, "qn0");         // valgrind  1 block possibly lost
   hisResZ2_2->Draw("surf1");  
   c.cd(2);
   MakeDiff(hisResZ2_2,frms)->Draw("surf1");
   c.Update();
//    c.SaveAs("RMSZPad2.eps");
   
   // write calculated resoltuion rms to 'rms.txt'
   ofstream filerms("rms.txt");
   filerms<<"Pad 0.75 cm"<<"\n";
   hisResY0_2->Fit(frms, "qn0");         // valgrind   23 blocks indirectly lost
   filerms<<"Y\t"<<frms->GetParameter(0)<<"\t"<<frms->GetParameter(1)<<"\t"<<frms->GetParameter(2)<<"\n";
   hisResZ0_2->Fit(frms, "qn0");         // valgrind   23 blocks indirectly lost
   filerms<<"Z\t"<<frms->GetParameter(0)<<"\t"<<frms->GetParameter(1)<<"\t"<<frms->GetParameter(2)<<"\n";
   //
   filerms<<"Pad 1.00 cm"<<1<<"\n";
   hisResY1_2->Fit(frms, "qn0");         // valgrind      3,256 bytes in 22 blocks are indirectly lost 
   filerms<<"Y\t"<<frms->GetParameter(0)<<"\t"<<frms->GetParameter(1)<<"\t"<<frms->GetParameter(2)<<"\n";
   hisResZ1_2->Fit(frms, "qn0");         // valgrind    66,036 bytes in 3 blocks are still reachable
   filerms<<"Z\t"<<frms->GetParameter(0)<<"\t"<<frms->GetParameter(1)<<"\t"<<frms->GetParameter(2)<<"\n";
   //
   filerms<<"Pad 1.50 cm"<<0<<"\n";
   hisResY2_2->Fit(frms, "qn0");      // valgrind   40,139 bytes in 11 blocks are still reachable   330,180 bytes in 15 blocks are possibly lost
   filerms<<"Y\t"<<frms->GetParameter(0)<<"\t"<<frms->GetParameter(1)<<"\t"<<frms->GetParameter(2)<<"\n";
   hisResZ2_2->Fit(frms, "qn0");
   filerms<<"Z\t"<<frms->GetParameter(0)<<"\t"<<frms->GetParameter(1)<<"\t"<<frms->GetParameter(2)<<"\n";

   TH1::AddDirectory(kFALSE);
   gSystem->ChangeDirectory("..");
   ps->Close();
   delete ps;
   delete frms;
}



void AliTPCcalibTracks::MakeResPlotsQ(Int_t minEntries,  Bool_t bDraw, char* pathName){
  //
  // make resolution Plots
  // not yet finished function
  //
  cout << " not yet finished function MakeResPlotsQ" << endl;
  gSystem->MakeDirectory(pathName);
  gSystem->ChangeDirectory(pathName);
  
  TLinearFitter fitter(3,"hyp2");
  TLinearFitter fitterQ(3,"hyp2");
  Int_t npointsFit =0;
  TF2 *fres = new TF2("fres","TMath::Sqrt([0]*[0]+[1]*[1]*x+[2]*[2]*y*y)",0,250,0,1);
  fres->SetParameter(0,0.02);
  fres->SetParameter(1,0.0054);
  fres->SetParameter(2,0.13);  
  fres->SetParLimits(0,0,0.1);
  TGraph2DErrors *graph = new TGraph2DErrors(1000);

  for (Int_t idim = 0; idim < 2; idim++){
    const char* dim = (idim==0)? "Y":"Z";
    for (Int_t ipad = 0; ipad < 3; ipad++){
      //
      printf("Direction %s\t",dim);
      printf("Pad %d\n",ipad);    
      npointsFit = 0;
      fitter.ClearPoints();
      fitterQ.ClearPoints();
      for (Int_t iq = 0; iq < 4; iq++){
	Int_t   bin   = GetBin(iq, ipad);
	Float_t qmean = GetQ(bin);
	char name[200];
	sprintf(name, "Resol%s Pad%d Qmiddle%f",dim, ipad, qmean);
        TH3F *hisin = 0;
        if (idim == 0) hisin = (TH3F*)fArrayQDY->At(bin);
        if (idim == 1) hisin = (TH3F*)fArrayQDZ->At(bin);
	if (!hisin) continue;
	//printf("Q\t%f\t",qmean);
        cout << "calling FitProjections" << endl;
	TObjArray * array = FitProjections(hisin, 2, minEntries, bDraw);
	//
	// Fit resolution
	//
        cout << "fit resolution" << endl;
         cout << "array->GetEntriesFast(): " << array->GetEntriesFast() << endl;
         //array->Print();
	for (Int_t ipoint = 0; ipoint < array->GetEntriesFast(); ipoint++){
	  TVectorF * vector = (TVectorF*)array->At(ipoint);
	  if (!vector) continue;
	  //vector->Print();
	  Double_t val   = (*vector)[4];	  
	  Double_t error = (*vector)[5];
	  error *= 2 * val;
	  //error += val * val * 0.05;  // 5% error addition
	  val *= val;
	  //
	  // mean fitter
	  //
	  Double_t zmean = (*vector)[0];
	  Double_t angle = (*vector)[2]* (*vector)[2];
	  Double_t x[2];
	  x[0] = zmean;
	  x[1] = angle;
	  fitter.AddPoint(x,val,error);
	  //
	  // with q fitter
	  //
	  Double_t zmeanq = (*vector)[0]/qmean;
	  Double_t angleq = (*vector)[2]* (*vector)[2];
	  Double_t xq[2];
	  xq[0] = zmeanq;
	  xq[1] = angleq;
	  fitterQ.AddPoint(xq,val,error);
	  graph->SetPoint(npointsFit,(*vector)[0],(*vector)[2],(*vector)[4]);
	  graph->SetPointError(npointsFit,(*vector)[1],(*vector)[3],(*vector)[5]);
	  npointsFit++;	  
	}
      }
      printf("NPoints = %d \n", npointsFit);
      cout << "evaluating fitters" << endl;
      fitter.Eval();
      fitterQ.Eval();
      TGraph2DErrors *graphVal = new TGraph2DErrors(npointsFit);
      TGraph2DErrors *graphDif = new TGraph2DErrors(npointsFit);
/*      graphVal->SetDirectory(0);
      graphDif->SetDirectory(0);*/
      for (Int_t ipoint=0; ipoint<npointsFit; ipoint++){
	Double_t x[3], sigma[3];
	x[0] = graph->GetX()[ipoint]+0.01*Float_t(ipad);
	x[1] = graph->GetY()[ipoint]+0.01*ipad;
	x[2] = graph->GetZ()[ipoint];
	sigma[0] = graph->GetErrorX(ipoint);
	sigma[1] = graph->GetErrorY(ipoint);
	sigma[3] = graph->GetErrorZ(ipoint);
	graphVal->SetPoint(ipoint,x[0],x[1],x[2]);
	graphVal->SetPointError(ipoint,sigma[0],sigma[1],sigma[2]);	 
      }
      fres->SetParameter(0,0.02);
      fres->SetParameter(1,0.0054);
      fres->SetParameter(2,0.13);  
      
      graphVal->Fit(fres,"q");
      for (Int_t ipoint=0; ipoint<npointsFit; ipoint++){
	Double_t x[3], sigma[3];
	x[0] = graph->GetX()[ipoint]+0.01*Float_t(ipad);
	x[1] = graph->GetY()[ipoint]+0.01*ipad;
	x[2] = graph->GetZ()[ipoint]-fres->Eval(x[0],x[1]);
	sigma[0] = graph->GetErrorX(ipoint);
	sigma[1] = graph->GetErrorY(ipoint);
	sigma[3] = graph->GetErrorZ(ipoint);
	graphDif->SetPoint(ipoint,x[0],x[1],x[2]);
	graphDif->SetPointError(ipoint,sigma[0],sigma[1],sigma[2]);	 
      }
      Double_t p0 = TMath::Sqrt(TMath::Abs(fitter.GetParameter(0)));
      Double_t p1 = TMath::Sqrt(TMath::Abs(fitter.GetParameter(1)));
      Double_t p2 = TMath::Sqrt(TMath::Abs(fitter.GetParameter(2)));
      Double_t chi2= fitter.GetChisquare()/npointsFit;
      printf("Linear fit  - chi2 %f\t%f\t%f\t%f\n",chi2,p0,p1,p2);
      Double_t p0q = TMath::Sqrt(TMath::Abs(fitterQ.GetParameter(0)));
      Double_t p1q = TMath::Sqrt(TMath::Abs(fitterQ.GetParameter(1)));
      Double_t p2q = TMath::Sqrt(TMath::Abs(fitterQ.GetParameter(2)));
      Double_t chi2q= fitterQ.GetChisquare()/npointsFit;
      printf("Linear fitQ - chi2 %f\t%f\t%f\t%f\n",chi2q,p0q,p1q,p2q);
 
      printf("Graph fit   - chi2 %f\t%f\t%f\t%f\n",fres->GetChisquare(),fres->GetParameter(0),fres->GetParameter(1),fres->GetParameter(2));     
    }
  }
  gSystem->ChangeDirectory("..");
}


TObjArray* AliTPCcalibTracks::FitProjections(TH3F * hfit, Int_t val, Int_t minEntry, Bool_t bDraw){
  //
  // ???
  // needed by MakeResPlotsQ
  //
  
  cout << "AliTPCcalibTracks::FitProjections started" << endl;
  TObjArray * array = new TObjArray(100);
  TAxis *xaxis  = hfit->GetXaxis();
  TAxis *yaxis  = hfit->GetYaxis();
  Double_t x[2];
  char name[200];
  sprintf(name,"Histos%s_%d", hfit->GetName(),val);
  TCanvas *canvas = 0;
  if (bDraw){
    canvas = new TCanvas(name,name);
    canvas->Divide(yaxis->GetNbins(),xaxis->GetNbins());
  }
  Int_t count =1;
  Int_t zaehler = 0;
  for (Int_t biny = 1; biny <= yaxis->GetNbins(); biny++) {
    x[1]  = yaxis->GetBinCenter(biny);
    for (Int_t binx = 1; binx <= xaxis->GetNbins(); binx++) {
      x[0]  = xaxis->GetBinCenter(binx);
      char name[200];
      sprintf(name, "%s x %f y %f", hfit->GetName(),x[0],x[1]);
      TH1D * projection = (TH1D*)( hfit->ProjectionZ(name, binx, binx, biny, biny) );
      projection->SetDirectory(0);
      //
      //
      if (projection->GetEntries() < minEntry){
	Float_t meanz = 0;
	Float_t meanangle = 0;
	Double_t entries = 0;
	projection->Clear();
	for (Int_t dbin = 0; dbin <= 4; dbin++)
	  for (Int_t dbiny2 = -3; dbiny2 <= 3; dbiny2++) {
	    for (Int_t dbinx2 = -3; dbinx2 <= 3; dbinx2++){
              zaehler++; 
	      if (TMath::Abs(dbinx2) + TMath::Abs(dbiny2) != dbin) continue;
	      Int_t binx2 = binx + dbinx2;
	      Int_t biny2 = biny + dbiny2;
// 	      if (binx2 < 1 || biny2 < 1 || binx2 > xaxis->GetNbins() || biny2 > yaxis->GetNbins()) continue;
	      if (binx2 < 1) continue;
	      if (biny2 < 1) continue;
	      if (binx2 > xaxis->GetNbins()) continue;
	      if (biny2 > yaxis->GetNbins()) continue;
	      TH1D * projection2 = (TH1D*)(hfit->ProjectionZ("Temp", binx2, binx2, biny2, biny2));
	      //projection2->SetDirectory(0);
	      projection->Add(projection2);
	      entries += projection2->GetEntries();
	      meanz     += projection2->GetEntries() * xaxis->GetBinCenter(binx2);
	      meanangle += projection2->GetEntries() * yaxis->GetBinCenter(biny2);
	      if (entries > minEntry) break;
	    }
	    if (entries > minEntry) break;
	  }
	if ( projection->GetEntries()<minEntry) continue;
	meanz/=entries;
	meanangle/=entries;
	x[0] = meanz;
	x[1] = meanangle;
      }
      
      if (projection->GetEntries() < minEntry) continue;
      //
      Float_t xmin = projection->GetMean()-2.*projection->GetRMS()-0.02;
      Float_t xmax = projection->GetMean()+2.*projection->GetRMS()+0.02;
      //      printf("%f\t%f\n",xmin,xmax);
      if (bDraw) canvas->cd(count);
      projection->Fit("gaus","q","",xmin,xmax);
      TVectorF *pvector = new TVectorF(6);
      TVectorF & vector = *pvector;
      vector[0] = x[0];
      vector[1] = xaxis->GetBinWidth(binx);
      vector[2] = x[1];
      vector[3] = yaxis->GetBinWidth(biny);
      vector[4] = projection->GetFunction("gaus")->GetParameter(2);
      vector[5] = projection->GetFunction("gaus")->GetParError(2);
      array->AddLast(pvector);
      count++;
    }    
  }
  
//   cout << "the inner loop was executed "<< zaehler << " times!" << endl;
  //
  //
  //
  
  TH1::AddDirectory(kTRUE);  // TH3F::FitSlicesZ() writes histograms into the current directory
 
//  cout << "critical point in FitProjections" << endl;
  
  TF2 *fres = new TF2("fres","TMath::Sqrt([0]*[0]+[1]*[1]*x+[2]*[2]*y*y)",0,250,0,1);
  fres->SetParameter(0,0.02);
  fres->SetParameter(1,0.0054);
  fres->SetParameter(2,0.13);  
  //
  hfit->FitSlicesZ();
  sprintf(name,"%s_%d", hfit->GetName(),val);
  TH2D * his_2 = (TH2D*)gDirectory->Get(name);
  his_2->Fit(fres,"q");
  his_2->SetDirectory(0);
  if (bDraw){
    TCanvas *canvas2 = new TCanvas( hfit->GetName(), hfit->GetName());
    canvas2->Divide(1,2);
    canvas2->cd(1);
    his_2->Draw("surf1");
    canvas2->cd(2);
    MakeDiff(his_2,fres)->Draw("surf1");
  }
  //
  printf("fres_Parameter[0]: %f\tfres_Parameter[1]: %f\tfres_Parameter[2]: %f\n",fres->GetParameter(0),fres->GetParameter(1),fres->GetParameter(2));
  
  TH1::AddDirectory(kFALSE);  
  cout << "end of FitProjections, array->GetEntriesFast():" << array->GetEntriesFast() << endl;
  return array;
}



void AliTPCcalibTracks::MakeResPlotsQTree(Int_t minEntries, char* pathName){
  //
  //  Make tree with resolution parameters
  //  the result is written to 'resol.root' in directory 'pathname'
  //
  
   cout << " ***** this is MakeResPlotsQTree *****" << endl;
   cout << "    relax, the calculation will take a while..." << endl;
  
   gSystem->MakeDirectory(pathName);
   gSystem->ChangeDirectory(pathName);
   TTreeSRedirector fTreeResol("resol.root");
   
   TH3F *resArray[2][3][11];
   TH3F *rmsArray[2][3][11];
  
   // load histograms into resArraz and rmsArray
   for (Int_t idim = 0; idim < 2; idim++){
      for (Int_t ipad = 0; ipad < 3; ipad++){
         for (Int_t iq = 0; iq <= 10; iq++){
            rmsArray[idim][ipad][iq]=0;
            resArray[idim][ipad][iq]=0;
            Int_t bin = GetBin(iq,ipad); 
//             Double_t qCenter = GetQ(bin);          // unused variable !!!
            
            TH3F *hresl = 0;
            if (idim == 0) hresl = (TH3F*)fArrayQDY->At(bin);
            if (idim == 1) hresl = (TH3F*)fArrayQDZ->At(bin);
            if (!hresl) continue;
            resArray[idim][ipad][iq] = (TH3F*) hresl->Clone();
            resArray[idim][ipad][iq]->SetDirectory(0);
            
            TH3F * hreslRMS = 0;
            if (idim == 0) hreslRMS = (TH3F*)fArrayQRMSY->At(bin);
            if (idim == 1) hreslRMS = (TH3F*)fArrayQRMSZ->At(bin);
            if (!hreslRMS) continue;
            rmsArray[idim][ipad][iq] = (TH3F*) hreslRMS->Clone();
            rmsArray[idim][ipad][iq]->SetDirectory(0);
         }
      }
   }
    
   cout << "Histograms loaded, starting to proces..." << endl;
   
   //--------------------------------------------------------------------------------------------
  
   //Double_t qCenter; 
   Double_t qMean;
   Double_t zMean, angleMean, zCenter, angleCenter;
   Double_t zSigma, angleSigma;
   Int_t loopCounter = 1;
  
   for (Int_t idim = 0; idim < 2; idim++){
      // Loop y-z corrdinate
      for (Int_t ipad = 0; ipad < 3; ipad++){
         // loop pad type
         
        // printf("%d\t%d\n", idim, ipad);
         for (Int_t iq = -1; iq < 10; iq++){
            // LOOP Q
            cout << "Loop-counter, this is loop " << loopCounter << " of 66, (" 
                 << (Int_t)((loopCounter)/66.*100) << "% done), " 
                 << "idim = " << idim << ", ipad = " << ipad << ", iq = " << iq << "  \r" << std::flush;
            loopCounter++;
            
            TH3F *hres = 0;
            TH3F *hrms = 0;
            qMean = 0;
            if (iq == -1){
               // integrated spectra
               Float_t entriesQ = 0;
               for (Int_t iql = 0; iql < 10; iql++){    
                  Int_t bin = GetBin(iql,ipad); 
                  TH3F *hresl = resArray[idim][ipad][iql];
                  TH3F *hrmsl = rmsArray[idim][ipad][iql];
                  if (!hresl) continue;
                  if (!hrmsl) continue;	    
                  entriesQ += hresl->GetEntries();
                  qMean += hresl->GetEntries() * GetQ(bin);      
                  if (!hres) {
                     hres = (TH3F*)hresl->Clone();
                     hrms = (TH3F*)hrmsl->Clone();
                  }
                  else{
                     hres->Add(hresl);
                     hrms->Add(hrmsl);
                  }
               }
               qMean /= entriesQ;
               qMean *= -1.;  // integral mean charge
            }
            else {
               Float_t entriesQ = 0;
               for (Int_t iql = iq - 1; iql <= iq + 1; iql++){		    
                  if (iql < 0) continue;
                  Int_t bin = GetBin(iql,ipad);
                  // qCenter   = GetQ(bin);  
                  TH3F * hresl = resArray[idim][ipad][iql];
                  TH3F * hrmsl = rmsArray[idim][ipad][iql];
                  if (!hresl) continue;
                  if (!hrmsl) continue;
                  entriesQ += hresl->GetEntries(); 
                  qMean += hresl->GetEntries() * GetQ(bin);      
                  if (!hres) {
                     hres = (TH3F*) hresl->Clone();
                     hrms = (TH3F*) hrmsl->Clone();
                  }
                  else{
                     hres->Add(hresl);
                     hrms->Add(hrmsl);
                  }
               }
               qMean/=entriesQ;
            }
      
            TAxis *xaxis  = hres->GetXaxis();
            TAxis *yaxis  = hres->GetYaxis();
            TAxis *zaxis  = hres->GetZaxis();
            TAxis *zaxisrms  = hrms->GetZaxis();
            for (Int_t biny = 1; biny <= yaxis->GetNbins(); biny++) {
               // angle loop
               angleCenter = yaxis->GetBinCenter(biny);
               for (Int_t binx = 1; binx <= xaxis->GetNbins(); binx++) {
                  // z - loop
                  zCenter    = xaxis->GetBinCenter(binx);
                  zMean      = zCenter;
                  angleMean  = angleCenter;
                  zSigma     = xaxis->GetBinWidth(binx);
                  angleSigma = yaxis->GetBinWidth(biny); 
                  
                  
                  // create 2 1D-Histograms, projectionRes and projectionRms
                  // here it is possible to speed up the program by using the loop and the other 
                  // TH1D *projectionRes... and TH1D *projectionRms... statements
                  // but there is a bug somewhere....
                  char name[200];
                  sprintf(name,"%s x %f y %f", hres->GetName(),zCenter,angleCenter);
//                   TH1D *projectionRes = new TH1D(name, name, zaxis->GetNbins(), zaxis->GetXmin(), zaxis->GetXmax());
//                  TH1D * projectionRes = (TH1D*)(hres->ProjectionZ(name,binx,binx, biny,biny));
                  TH1D * projectionRes = new TH1D(name,name,zaxis->GetNbins(),zaxis->GetXmin(), zaxis->GetXmax());
                  

                  sprintf(name,"%s x %f y %f", hrms->GetName(),zCenter,angleCenter);
//                   TH1D *projectionRms = new TH1D(name, name, zaxis->GetNbins(), zaxis->GetXmin(), zaxis->GetXmax());
//                  TH1D * projectionRms = (TH1D*)(hrms->ProjectionZ(name,binx,binx, biny,biny));
                  TH1D * projectionRms =  new TH1D(name,name,zaxis->GetNbins(),zaxisrms->GetXmin(), zaxisrms->GetXmax());
                  
/*                  
                  for (Int_t ibin3 = 1; ibin3 < zaxis->GetNbins(); ibin3++) {
//                   fill 1D-Histograms
                     projectionRes->Fill(ibin3, hres->GetBinContent(binx, biny, ibin3));
                     projectionRms->Fill(ibin3, hrms->GetBinContent(binx, biny, ibin3));
                  }
*/                  
                  projectionRes->SetDirectory(0);
                  projectionRms->SetDirectory(0);
                  Double_t entries = projectionRes->GetEntries();
		  Int_t    nbins   =0;
                  if (projectionRes->GetEntries() < minEntries){
                     // if not enough statistic
                     zMean = 0;
                     angleMean = 0;
                     entries =0;
                     //projectionRes->Clear();
                     projectionRms->Clear();	      
                     for (Int_t dbin = 0; dbin <= 8; dbin++){
                        for (Int_t dbiny2 = -1; dbiny2 <= 1; dbiny2++) {
                           for (Int_t dbinx2 = -3; dbinx2 <= 3; dbinx2++){
                              if (TMath::Abs(dbinx2) + TMath::Abs(dbiny2) != dbin) continue;
                              Int_t binx2 = binx + dbinx2;
                              Int_t biny2 = biny + dbiny2;
                              if (binx2 < 1) continue;
                              if (biny2 < 1) continue;
                              if (binx2 >= xaxis->GetNbins()) continue;
                              if (biny2 >= yaxis->GetNbins()) continue;
			      nbins++;
			      //
			      //
			      // Fill resolution histo
                              for (Int_t ibin3 = 1; ibin3 < zaxis->GetNbins(); ibin3++) {
				Int_t content = hres->GetBinContent(binx2, biny2, ibin3);
				
				projectionRes->Fill(zaxis->GetBinCenter(ibin3), hres->GetBinContent(binx2, biny2, ibin3));
                                 entries   += hres->GetBinContent(binx2, biny2, ibin3);
                                 zMean     += hres->GetBinContent(binx2, biny2, ibin3) * xaxis->GetBinCenter(binx2);
                                 angleMean += hres->GetBinContent(binx2, biny2, ibin3) * yaxis->GetBinCenter(biny2);
                              }  // ibin3 loop
			      // fill RMS histo
			      for (Int_t ibin3 = 1; ibin3 < zaxisrms->GetNbins(); ibin3++) {
				projectionRms->Fill(zaxisrms->GetBinCenter(ibin3), hrms->GetBinContent(binx2, biny2, ibin3));
			      }

                           }  //dbinx2 loop
                        }  // dbiny2 loop
			if (entries > minEntries) break;
		     }
                     if ( entries< minEntries) continue;
                     zMean /= entries;
                     angleMean /= entries;
                  }     // if (projectionRes->GetEntries() < minEntries)
                  
                  if (projectionRes->GetSum() > minEntries) {
                     //  when enough statistic is accumulated
                     Float_t entries2 = projectionRes->GetSum();
                     Float_t xmin    = projectionRes->GetMean() - 2. * projectionRes->GetRMS() - 0.2;
                     Float_t xmax    = projectionRes->GetMean() + 2. * projectionRes->GetRMS() + 0.2;
                     projectionRes->Fit("gaus","q","",xmin,xmax);
                     Float_t resol   = projectionRes->GetFunction("gaus")->GetParameter(2);
                     Float_t sigma   = projectionRes->GetFunction("gaus")->GetParError(2);
		     Float_t meanR   = projectionRes->GetMean();
		     Float_t sigmaR  = projectionRes->GetRMS();
		     
                     //
                     xmin = projectionRms->GetMean() - 2. * projectionRes->GetRMS() - 0.2;
                     xmax = projectionRms->GetMean() + 2. * projectionRes->GetRMS() + 0.2;
                     projectionRms->Fit("gaus","q","",xmin,xmax);
                     Float_t rmsMean    = projectionRms->GetFunction("gaus")->GetParameter(1);
                     Float_t errorRMS   = projectionRms->GetFunction("gaus")->GetParError(1);
                     Float_t rmsSigma   = projectionRms->GetFunction("gaus")->GetParameter(2);
                     Float_t errorSigma = projectionRms->GetFunction("gaus")->GetParError(2);
                     //
                     Float_t length = 0.75;
                     if (ipad == 1) length = 1;
                     if (ipad == 2) length = 1.5;
                     
                     fTreeResol<<"Resol"<<
                        "Entries="<<entries<<
		       "Entries2="<<entries2<<
		       "nbins="<<nbins<<
                        "Dim="<<idim<<
                        "Pad="<<ipad<<
                        "Length="<<length<<
                        "QMean="<<qMean<<
                        "Zc="<<zCenter<<
                        "Zm="<<zMean<<
                        "Zs="<<zSigma<<
                        "AngleC="<<angleCenter<<
                        "AngleM="<<angleMean<<
                        "AngleS="<<angleSigma<<
                        "Resol="<<resol<<
                        "Sigma="<<sigma<<
                        "MeanR="<<meanR<<
                        "SigmaR="<<sigmaR<<
		        //
                        "RMSm="<<rmsMean<<
                        "RMSs="<<rmsSigma<<
                        "RMSe0="<<errorRMS<<
                        "RMSe1="<<errorSigma<<
                        "\n";
		     /*
		     projectionRes->SetDirectory(fTreeResol.GetFile());
		     projectionRes->Write(projectionRes->GetName());
		     projectionRes->SetDirectory(0);
		     projectionRms->SetDirectory(fTreeResol.GetFile());
		     projectionRms->Write(projectionRms->GetName());
		     projectionRes->SetDirectory(0);
		     */
                  }
                  delete projectionRes;
                  delete projectionRms;
               }
            }
         }
      }
   }
   cout << endl;
   cout << "MakeResPlotsQTree done, results are in 'resol.root'." << endl;
   gSystem->ChangeDirectory("..");

}



Long64_t AliTPCcalibTracks::Merge(TCollection *collectionList) {
   // 
   // function to merge several AliTPCcalibTracks objects after PROOF calculation
   // The object's histograms are merged via their merge functions
   // 
   
   cout << " *****  this is AliTPCcalibTracks::Merge(TCollection *collectionList)  *****"<< endl;  
   if (!collectionList) return 0;
   if (collectionList->IsEmpty()) return -1;
   
   cout << "the collectionList contains " << collectionList->GetEntries() << " entries." << endl;     //    REMOVE THIS LINE!!!!!!!!!!!!!!!!!1
   
   // create a list for each data member
   TList* deltaYList = new TList;
   TList* deltaZList = new TList;
   TList* arrayAmpRowList = new TList;
   TList* arrayAmpList = new TList;
   TList* arrayQDYList = new TList;
   TList* arrayQDZList = new TList;
   TList* arrayQRMSYList = new TList;
   TList* arrayQRMSZList = new TList;
   TList* resolYList = new TList;
   TList* resolZList = new TList;
   TList* rMSYList = new TList;
   TList* rMSZList = new TList;
   
   TList* nRowsList = new TList;
   TList* nSectList = new TList;
   TList* fileNoList = new TList;
   
   TIterator *listIterator = collectionList->MakeIterator();
   AliTPCcalibTracks *calibTracks = 0;
   cout << "start to iterate, filling lists" << endl;                      //    REMOVE THIS LINE!!!!!!!!!!!!!!!!!1
   Int_t counter = 0;
   while ( (calibTracks = (AliTPCcalibTracks*)listIterator->Next()) ){
      // loop over all entries in the collectionList and get dataMembers into lists
      if (!calibTracks) continue;
      deltaYList->Add( calibTracks->GetfDeltaY() );
      deltaZList->Add( calibTracks->GetfDeltaZ() );
      arrayAmpRowList->Add(calibTracks->GetfArrayAmpRow());
      arrayAmpList->Add(calibTracks->GetfArrayAmp());
      arrayQDYList->Add(calibTracks->GetfArrayQDY());
      arrayQDZList->Add(calibTracks->GetfArrayQDZ());
      arrayQRMSYList->Add(calibTracks->GetfArrayQRMSY());
      arrayQRMSZList->Add(calibTracks->GetfArrayQRMSZ());
      resolYList->Add(calibTracks->GetfResolY());
      resolZList->Add(calibTracks->GetfResolZ());
      rMSYList->Add(calibTracks->GetfRMSY());
      rMSZList->Add(calibTracks->GetfRMSZ());
      counter++;
   }
   
   // reset data members
   cout << "histogram's reset-functins are called... " << endl; //    REMOVE THIS LINE!!!!!!!!!!!!!!!!!1
   fDeltaY->Reset();
   fDeltaZ->Reset();
   for (Int_t i = 0; i < fArrayAmpRow->GetEntriesFast(); i++ ) 
      ((TProfile*)(fArrayAmpRow->At(i)))->Reset();
   for (Int_t i = 0; i < fArrayAmp->GetEntriesFast(); i++ ) 
      ((TProfile*)(fArrayAmp->At(i)))->Reset();
   for (Int_t i = 0; i < fArrayQDY->GetEntriesFast(); i++)
      ((TH3F*)(fArrayQDY->At(i)))->Reset();
   for (Int_t i = 0; i < fArrayQDZ->GetEntriesFast(); i++)
      ((TH3F*)(fArrayQDZ->At(i)))->Reset();
   for (Int_t i = 0; i < fArrayQRMSY->GetEntriesFast(); i++)
      ((TH3F*)(fArrayQRMSY->At(i)))->Reset();
   for (Int_t i = 0; i < fArrayQRMSZ->GetEntriesFast(); i++)
      ((TH3F*)(fArrayQRMSZ->At(i)))->Reset();
   for (Int_t i = 0; i < fResolY->GetEntriesFast(); i++) {
      ((TH3F*)(fResolY->At(i)))->Reset();
      ((TH3F*)(fResolZ->At(i)))->Reset();
      ((TH3F*)(fRMSY->At(i)))->Reset();
      ((TH3F*)(fRMSZ->At(i)))->Reset();
   }
               
   // merge data members
   cout << "histogram's merge-functins are called... " << endl;    //    REMOVE THIS LINE!!!!!!!!!!!!!!!!!1
   fDeltaY->Merge(deltaYList);
   fDeltaZ->Merge(deltaZList);
   
   TObjArray* objarray = 0;
   TH1* hist = 0;
   TList* histList = 0;
   TIterator *objListIterator = 0;
   
   cout << "merging fArrayAmpRows..." << endl;
   // merge fArrayAmpRows
   for (Int_t i = 0; i < fArrayAmpRow->GetEntriesFast(); i++ ) {  // loop over data member, i<72
      objListIterator = arrayAmpRowList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayAmpRowList, get TObjArray, get object at position i, cast it into TProfile
         if (!objarray) continue;
         hist = (TProfile*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TProfile*)(fArrayAmpRow->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }
   
   cout << "merging fArrayAmps..." << endl;
   // merge fArrayAmps
   for (Int_t i = 0; i < fArrayAmp->GetEntriesFast(); i++ ) {  // loop over data member, i<72
      TIterator *objListIterator = arrayAmpList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayAmpList, get TObjArray, get object at position i, cast it into TH1F
         if (!objarray) continue;
         hist = (TH1F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH1F*)(fArrayAmp->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }
   
   cout << "merging fArrayQDY..." << endl;
   // merge fArrayQDY
   for (Int_t i = 0; i < fArrayQDY->GetEntriesFast(); i++) { // loop over data member, i < 300
      objListIterator = arrayQDYList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayQDYList, get TObjArray, get object at position i, cast it into TH3F
         if (!objarray) continue;
         hist = (TH3F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH3F*)(fArrayQDY->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }

   cout << "merging fArrayQDZ..." << endl;
   // merge fArrayQDZ
   for (Int_t i = 0; i < fArrayQDZ->GetEntriesFast(); i++) { // loop over data member, i < 300
      objListIterator = arrayQDZList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayQDZList, get TObjArray, get object at position i, cast it into TH3F
         if (!objarray) continue;
         hist = (TH3F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH3F*)(fArrayQDZ->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }

   cout << "merging fArrayQRMSY..." << endl;
   // merge fArrayQRMSY
   for (Int_t i = 0; i < fArrayQRMSY->GetEntriesFast(); i++) { // loop over data member, i < 300
      objListIterator = arrayQRMSYList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayQDZList, get TObjArray, get object at position i, cast it into TH3F
         if (!objarray) continue;
         hist = (TH3F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH3F*)(fArrayQRMSY->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }   

   cout << "merging fArrayQRMSZ..." << endl;
   // merge fArrayQRMSZ
   for (Int_t i = 0; i < fArrayQRMSZ->GetEntriesFast(); i++) { // loop over data member, i < 300
      objListIterator = arrayQRMSZList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayQDZList, get TObjArray, get object at position i, cast it into TH3F
         if (!objarray) continue;
         hist = (TH3F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH3F*)(fArrayQRMSZ->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }      
   
   cout << "starting to merge the rest: fResolY, fResolZ , fRMSY, fRMSZ..." << endl;
   // merge fResolY
   for (Int_t i = 0; i < fResolY->GetEntriesFast(); i++) { // loop over data member, i < 3
      objListIterator = resolYList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayQDZList, get TObjArray, get object at position i, cast it into TH3F
         if (!objarray) continue;
         hist = (TH3F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH3F*)(fResolY->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }
   
   // merge fResolZ
   for (Int_t i = 0; i < fResolZ->GetEntriesFast(); i++) { // loop over data member, i < 3
      objListIterator = resolZList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayQDZList, get TObjArray, get object at position i, cast it into TH3F
         if (!objarray) continue;
         hist = (TH3F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH3F*)(fResolZ->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }

   // merge fRMSY
   for (Int_t i = 0; i < fRMSY->GetEntriesFast(); i++) { // loop over data member, i < 3
      objListIterator = rMSYList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayQDZList, get TObjArray, get object at position i, cast it into TH3F
         if (!objarray) continue;
         hist = (TH3F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH3F*)(fRMSY->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }
         
   // merge fRMSZ
   for (Int_t i = 0; i < fRMSZ->GetEntriesFast(); i++) { // loop over data member, i < 3
      objListIterator = rMSZList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayQDZList, get TObjArray, get object at position i, cast it into TH3F
         if (!objarray) continue;
         hist = (TH3F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH3F*)(fRMSZ->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }
   
   delete deltaYList;
   delete deltaZList;
   delete arrayAmpRowList;
   delete arrayAmpList;
   delete arrayQDYList;
   delete arrayQDZList;
   delete arrayQRMSYList;
   delete arrayQRMSZList;
   delete resolYList;
   delete resolZList;
   delete rMSYList;
   delete rMSZList;
   delete nRowsList;
   delete nSectList;
   delete fileNoList;
   delete listIterator;
   
   cout << "merging done!" << endl;
   
   return 1;
}


AliTPCcalibTracks* AliTPCcalibTracks::TestMerge(AliTPCcalibTracks *ct, AliTPCClusterParam *clusterParam, Int_t nCalTracks){
   // 
   // function to test AliTPCcalibTrack::Merge:
   // in the file 'f' is a AliTPCcalibTrack object with name "calibTracks"
   // this object is appended 'nCalTracks' times to a TList
   // A new AliTPCcalibTrack object is created which merges the list
   // this object is returned
   // 
   /*
   .L AliTPCcalibTracks.cxx+g
   TFile f("Output.root");
   AliTPCcalibTracks* calTracks = (AliTPCcalibTracks*)f.Get("calibTracks");
   //f.Close();
   TFile clusterParamFile("/u/lbozyk/calibration/workdir/calibTracks/TPCClusterParam.root");
   AliTPCClusterParam *clusterParam  =  (AliTPCClusterParam *) clusterParamFile.Get("Param"); 
   clusterParamFile.Close();

   AliTPCcalibTracks::TestMerge(calTracks, clusterParam);
   */
   TList *list = new TList();
   if (ct == 0 || clusterParam == 0) return 0;
   cout << "making list with " << nCalTracks << " AliTPCcalibTrack objects" << endl;
   for (Int_t i = 0; i < nCalTracks; i++) {
      list->Add(new AliTPCcalibTracks(ct));
      if (i%10==0) cout << "Adding element " << i << " of " << nCalTracks << endl;
   }
   
   // only for check at the end
   AliTPCcalibTracks* cal1 = new AliTPCcalibTracks(ct);
   Double_t cal1Entries = ((TH1F*)cal1->GetfArrayAmpRow()->At(5))->GetEntries();
//    Double_t cal1Entries = 5; //((TH1F*)ct->GetfArrayAmpRow()->At(5))->GetEntries();

   cout  << "The list contains " << list->GetEntries() << " entries. " << endl;
   
   
   AliTPCcalibTracksCuts *cuts = new AliTPCcalibTracksCuts(20, 0.4, 0.5, 0.13, 0.018);
   AliTPCcalibTracks* cal = new AliTPCcalibTracks("calTracksMerged", "calTracksMerged", clusterParam, cuts);
   cal->Merge(list);
   
   cout << "cal->GetfArrayAmpRow()->At(5)->Print():" << endl;
   cal->GetfArrayAmpRow()->At(5)->Print();
   Double_t calEntries = ((TH1F*)cal->GetfArrayAmpRow()->At(5))->GetEntries();
   
   cout << "cal1->GetfArrayAmpRow()->At(5))->GetEntries() = " << cal1Entries << endl;
   cout << " cal->GetfArrayAmpRow()->At(5))->GetEntries() = " <<  calEntries << endl;
   printf("That means there were %f / %f = %f AliTPCcalibTracks-Objects merged. \n", 
      calEntries, cal1Entries, ((Double_t)calEntries/cal1Entries));
   
   return cal;

}



