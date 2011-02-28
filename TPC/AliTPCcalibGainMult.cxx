
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


Basic calibration and QA class for the TPC gain calibration based on tracks from BEAM events.


Send comments etc. to: A.Kalweit@gsi.de, marian.ivanov@cern.ch
*/


#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TProfile.h"

#include "AliTPCcalibDB.h"
#include "AliTPCclusterMI.h"
#include "AliTPCClusterParam.h"
#include "AliTPCseed.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliTPCParam.h"

#include "AliComplexCluster.h"
#include "AliTPCclusterMI.h"

#include "AliLog.h"

#include "AliTPCcalibGainMult.h"

#include "TTreeStream.h"


ClassImp(AliTPCcalibGainMult)


AliTPCcalibGainMult::AliTPCcalibGainMult() 
  :AliTPCcalibBase(),
   fMIP(0),
   fLowerTrunc(0),
   fUpperTrunc(0),
   fUseMax(kFALSE),
   fHistNTracks(0),
   fHistClusterShape(0),
   fHistQA(0),
   fHistGainSector(0),
   fHistPadEqual(0),
   fHistGainMult(0)
{  
  //
  // Empty default cosntructor
  //
  AliInfo("Default Constructor");  
}


AliTPCcalibGainMult::AliTPCcalibGainMult(const Text_t *name, const Text_t *title) 
  :AliTPCcalibBase(),
   fMIP(0),
   fLowerTrunc(0),
   fUpperTrunc(0),
   fUseMax(kFALSE),
   fHistNTracks(0),
   fHistClusterShape(0),
   fHistQA(0),
   fHistGainSector(0),
   fHistPadEqual(0),
   fHistGainMult(0)
{
  //
  //
  //  
  SetName(name);
  SetTitle(title);
  //
  fMIP = 50.;
  fLowerTrunc = 0.02; // IMPORTANT CHANGE --> REMOVE HARDWIRED TRUNCATION FROM TRACKER
  fUpperTrunc = 0.6;
  fUseMax = kTRUE; // IMPORTANT CHANGE FOR PbPb; standard: kFALSE;
  //
  fHistNTracks = new TH1F("ntracks","Number of Tracks per Event; number of tracks per event; number of tracks",1001,-0.5,1000.5);
  fHistClusterShape = new TH1F("fHistClusterShape","cluster shape; rms meas. / rms exp.;",300,0,3);
  fHistQA = new TH3F("fHistQA","dEdx; momentum p (GeV); TPC signal (a.u.); pad",500,0.1,20.,500,0.,500,6,-0.5,5.5);
  BinLogX(fHistQA);
  //
  //
  //                          MIP, sect,  pad (short,med,long,full,oroc),   run,      ncl
  Int_t binsGainSec[5]    = { 145,   72,    4,  10000000,   65};
  Double_t xminGainSec[5] = { 10., -0.5, -0.5,      -0.5, -0.5}; 
  Double_t xmaxGainSec[5] = {300., 71.5,  3.5, 9999999.5, 64.5};
  TString axisNameSec[5]={"Q","sector","pad type","run", "ncl"};
  TString axisTitleSec[5]={"Q (a.u)","sector","pad type","run","ncl"};
  //
  fHistGainSector = new THnSparseF("fHistGainSector","0:MIP, 1:sect, 2:pad, 3:run, 4:ncl", 5, binsGainSec, xminGainSec, xmaxGainSec);
  for (Int_t iaxis=0; iaxis<5;iaxis++){
    fHistGainSector->GetAxis(iaxis)->SetName(axisNameSec[iaxis]);
    fHistGainSector->GetAxis(iaxis)->SetTitle(axisTitleSec[iaxis]);
  }
  //
  //
  //
  Int_t binsPadEqual[6]    = { 200, 200,    4,   20,   50, 100};
  Double_t xminPadEqual[6] = { 0.5, 0.5, -0.5,    0, -250,   0}; 
  Double_t xmaxPadEqual[6] = { 1.5, 1.5,  3.5, 13000,  250,   3};
  TString axisNamePadEqual[6]   = {"dEdxRatioMax","dEdxRatioTot","padType","mult","driftlength", "1_pt"};
  TString axisTitlePadEqual[6]  = {"dEdx_padRegion/mean_dEdx Qmax", "dEdx_padRegion/mean_dEdx Qtot","padType","mult","driftlength", "1/pt"};
  //
  fHistPadEqual = new THnSparseF("fHistPadEqual","0:dEdx_pad/dEdx_mean, 1:pad, 2:mult, 3:drift, 4:1/pt", 6, binsPadEqual, xminPadEqual, xmaxPadEqual);
  for (Int_t iaxis=0; iaxis<6;iaxis++){
    fHistPadEqual->GetAxis(iaxis)->SetName(axisNamePadEqual[iaxis]);
    fHistPadEqual->GetAxis(iaxis)->SetTitle(axisTitlePadEqual[iaxis]);
  }
  //
  //
  //                    MIP Qmax, MIP Qtot,  z,  pad, vtx. contribut., ncl
  Int_t binsGainMult[6]    = { 145,  145,   25,    4,  100,  80};
  Double_t xminGainMult[6] = { 10.,  10.,    0, -0.5,    0, -0.5}; 
  Double_t xmaxGainMult[6] = {300., 300.,  250,  3.5, 13000, 159.5};
  TString axisNameMult[6]={"Qmax","Qtot","drift","padtype""multiplicity","ncl"};
  TString axisTitleMult[6]={"Qmax (a.u)","Qtot (a.u.)","driftlenght l (cm)","Pad Type","multiplicity","ncl"};
  //
  fHistGainMult = new THnSparseF("fHistGainMult","MIP Qmax, MIP Qtot, z, type, vtx. contribut.", 6, binsGainMult, xminGainMult, xmaxGainMult); 
  for (Int_t iaxis=0; iaxis<6;iaxis++){
    fHistGainMult->GetAxis(iaxis)->SetName(axisNameMult[iaxis]);
    fHistGainMult->GetAxis(iaxis)->SetTitle(axisTitleMult[iaxis]);
  }
  //
  AliInfo("Non Default Constructor");  
  //
}


AliTPCcalibGainMult::~AliTPCcalibGainMult(){
  //
  //
  //
  delete fHistNTracks;            //  histogram showing number of ESD tracks per event
  delete fHistClusterShape;       //  histogram to check the cluster shape
  delete fHistQA;                 //  dE/dx histogram showing the final spectrum
  //
  delete fHistGainSector;   //  histogram which shows MIP peak for each of the 3x36 sectors (pad region)
  delete fHistPadEqual;     //  histogram for the equalization of the gain in the different pad regions -> pass0
  delete fHistGainMult;     //  histogram which shows decrease of MIP signal as a function


}



void AliTPCcalibGainMult::Process(AliESDEvent *event) {
  //
  //
  //
  if (!event) {
    Printf("ERROR: ESD not available");
    return;
  }  
  Int_t ntracks=event->GetNumberOfTracks(); 
  fHistNTracks->Fill(ntracks);
  
  AliESDfriend *esdFriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!esdFriend) {
   Printf("ERROR: esdFriend not available");
   return;
  }
  UInt_t runNumber = event->GetRunNumber();
  Int_t nContributors = event->GetNumberOfTracks();
  //
  // track loop
  //
  for (Int_t i=0;i<ntracks;++i) {
    //
    AliESDtrack *track = event->GetTrack(i);
    if (!track) continue;
    //   
    AliExternalTrackParam * trackIn  = (AliExternalTrackParam *)track->GetInnerParam();
    if (!trackIn) continue;
  
    // calculate necessary track parameters
    Double_t meanP = trackIn->GetP();
    Int_t ncls = track->GetTPCNcls();

    if (ncls < 80) continue;     
    
    // exclude tracks which do not look like primaries or are simply too short or on wrong sectors

    if (TMath::Abs(trackIn->Eta()) > 0.8) continue;
    UInt_t status = track->GetStatus();
    if ((status&AliESDtrack::kTPCrefit)==0) continue;
    //if (track->GetNcls(0) < 3) continue; // ITS clusters
    Float_t dca[2], cov[3];
    track->GetImpactParameters(dca,cov);
    Float_t primVtxDCA = TMath::Sqrt(dca[0]*dca[0]);
    if (primVtxDCA > 10 || primVtxDCA < 0.00001) continue;
    if (TMath::Abs(dca[1]) > 5) continue;
    //
    // require that the track does not cross any dead area
    //
    //if (track->GetTPCNclsF() < 158) continue;
    //
    //if (seed->CookShape(1) > 1) continue;
    //if (TMath::Abs(trackIn->GetY()) > 20) continue;
    //if (TMath::Abs(d)>20) continue;   // distance to the 0,0; select only tracks which cross chambers under proper angle
    //if (TMath::Abs(trackIn->GetSnp()) > 0.6) continue;
    if (primVtxDCA < 3 && track->GetNcls(0) > 3 && track->GetKinkIndex(0) == 0 && ncls > 100) fHistQA->Fill(meanP, track->GetTPCsignal(), 5);

    // Get seeds
    AliESDfriendTrack *friendTrack = esdFriend->GetTrack(i);
    if (!friendTrack) continue;
    TObject *calibObject;
    AliTPCseed *seed = 0;
    for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
      if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }    

    if (seed) {
      //
      const AliExternalTrackParam * trackOut = friendTrack->GetTPCOut();
      if (!trackIn) continue;
      if (!trackOut) continue;
      Double_t meanDrift = 250 - 0.5*TMath::Abs(trackIn->GetZ() + trackOut->GetZ());
      //
      for (Int_t irow =0; irow<160;irow++)    {
	AliTPCTrackerPoint * point = seed->GetTrackPoint(irow);
	if (point==0) continue;
	AliTPCclusterMI * cl = seed->GetClusterPointer(irow);
	if (cl==0) continue;	
	//
	Float_t rsigmay =  TMath::Sqrt(point->GetSigmaY());
	fHistClusterShape->Fill(rsigmay);
      }
      //
      Int_t row0 = 0;
      Int_t row1 = 160;
      //
      Double_t signalShortMax = seed->CookdEdxAnalytical(fLowerTrunc,fUpperTrunc,1,0,62);
      Double_t signalMedMax   = seed->CookdEdxAnalytical(fLowerTrunc,fUpperTrunc,1,63,126);
      Double_t signalLongMax  = seed->CookdEdxAnalytical(fLowerTrunc,fUpperTrunc,1,127,159);
      Double_t signalMax      = seed->CookdEdxAnalytical(fLowerTrunc,fUpperTrunc,1,row0,row1);
      Double_t signalArrayMax[4] = {signalShortMax, signalMedMax, signalLongMax, signalMax};
      //
      Double_t signalShortTot = seed->CookdEdxAnalytical(fLowerTrunc,fUpperTrunc,0,0,62);
      Double_t signalMedTot   = seed->CookdEdxAnalytical(fLowerTrunc,fUpperTrunc,0,63,126);
      Double_t signalLongTot  = seed->CookdEdxAnalytical(fLowerTrunc,fUpperTrunc,0,127,159);
      Double_t signalTot      = seed->CookdEdxAnalytical(fLowerTrunc,fUpperTrunc,0,row0,row1);
      Double_t signalArrayTot[4] = {signalShortTot, signalMedTot, signalLongTot, signalTot};
      //
      Double_t mipSignalShort = fUseMax ? signalShortMax : signalShortTot;
      Double_t mipSignalMed   = fUseMax ? signalMedMax   : signalMedTot;
      Double_t mipSignalLong  = fUseMax ? signalLongMax  : signalLongTot;
      Double_t mipSignalOroc  = seed->CookdEdxAnalytical(fLowerTrunc,fUpperTrunc,fUseMax,63,159);
      Double_t signal =  fUseMax ? signalMax  : signalTot;
      //
      fHistQA->Fill(meanP, mipSignalShort, 0);
      fHistQA->Fill(meanP, mipSignalMed, 1);
      fHistQA->Fill(meanP, mipSignalLong, 2);
      fHistQA->Fill(meanP, signal, 3);
      fHistQA->Fill(meanP, mipSignalOroc, 4);
      //
      // "dEdxRatioMax","dEdxRatioTot","padType","mult","driftlength", "1_pt"
      Float_t meanMax = (1/3.)*(signalArrayMax[0] + signalArrayMax[1] + signalArrayMax[2]);
      Float_t meanTot = (1/3.)*(signalArrayTot[0] + signalArrayTot[1] + signalArrayTot[2]); 
      if (meanMax < 1e-5 || meanTot < 1e-5) continue;
      for(Int_t ipad = 0; ipad < 4; ipad ++) {
	Double_t vecPadEqual[6] = {signalArrayMax[ipad]/meanMax, signalArrayTot[ipad]/meanTot, ipad, nContributors, meanDrift, track->OneOverPt()};
	fHistPadEqual->Fill(vecPadEqual);
      }
      //
      if (meanP > 0.4 && meanP < 0.55) {
	Double_t vecMult[6] = {seed->CookdEdxAnalytical(fLowerTrunc,fUpperTrunc,1,row0,row1),
			       seed->CookdEdxAnalytical(fLowerTrunc,fUpperTrunc,0,row0,row1),
			       meanDrift,
			       3,
			       nContributors,
			       ncls};
	//
	fHistGainMult->Fill(vecMult);
	vecMult[0]=mipSignalShort; vecMult[1]=mipSignalShort; vecMult[3]=0;
	fHistGainMult->Fill(vecMult);
	vecMult[0]=mipSignalMed; vecMult[1]=mipSignalMed; vecMult[3]=1;
	fHistGainMult->Fill(vecMult);
	vecMult[0]=mipSignalLong; vecMult[1]=mipSignalLong; vecMult[3]=2;
	fHistGainMult->Fill(vecMult);
	//
      }
      //
      //
      if (meanP > 0.5 || meanP < 0.4) continue; // only MIP pions
      //
      // for each track, we look at the three different pad regions, split it into tracklets, check if the sector does not change and fill the histogram
      //
      Bool_t isNotSplit[3] = {kTRUE, kTRUE, kTRUE}; //  short, medium, long (true if the track is not split between two chambers)
      //
      Double_t sector[4] = {-1, -1, -1, -1}; // sector number short, medium, long, all
      Int_t ncl[3] = {0,0,0};
      //
      for (Int_t irow=0; irow < 159; irow++){
	Int_t padRegion = 0;
	if (irow > 62) padRegion = 1;
	if (irow > 126) padRegion = 2;
	//
	AliTPCclusterMI* cluster = seed->GetClusterPointer(irow);
	if (!cluster) continue;
	if (sector[padRegion] == -1) {
	  sector[padRegion] = cluster->GetDetector();
	  continue;
	}
	if (sector[padRegion] != -1 && sector[padRegion] != cluster->GetDetector()) isNotSplit[padRegion] = kFALSE;
	ncl[padRegion]++;
      }
      //
      //                        MIP, sect,  pad,   run
      //
      Double_t vecMip[5] = {mipSignalShort, mipSignalMed, mipSignalLong, signal, mipSignalOroc};
      //
      for(Int_t ipad = 0; ipad < 3; ipad++) {
	//
	Double_t vecGainSec[5] = {vecMip[ipad], sector[ipad], ipad, runNumber, ncl[ipad]};
	if (isNotSplit[ipad]) fHistGainSector->Fill(vecGainSec);
      }
    }
   
  }    
}  


void AliTPCcalibGainMult::MakeLookup(THnSparse * /*hist*/, Char_t * /*outputFile*/) {
  //
  // Not  yet implemented
  //
}


void AliTPCcalibGainMult::Analyze() {


  return;

}


Long64_t AliTPCcalibGainMult::Merge(TCollection *li) {

  TIterator* iter = li->MakeIterator();
  AliTPCcalibGainMult* cal = 0;

  while ((cal = (AliTPCcalibGainMult*)iter->Next())) {
    if (!cal->InheritsFrom(AliTPCcalibGainMult::Class())) {
      Error("Merge","Attempt to add object of class %s to a %s", cal->ClassName(), this->ClassName());
      return -1;
    }
    
    if (cal->GetHistNTracks()) fHistNTracks->Add(cal->GetHistNTracks());
    if (cal->GetHistClusterShape()) fHistClusterShape->Add(cal->GetHistClusterShape());
    if (cal->GetHistQA()) fHistQA->Add(cal->GetHistQA());
    if (cal->GetHistGainSector()) fHistGainSector->Add(cal->GetHistGainSector());
    if (cal->GetHistPadEqual()) fHistPadEqual->Add(cal->GetHistPadEqual());
    if (cal->GetHistGainMult()) fHistGainMult->Add(cal->GetHistGainMult());
 
  }
  
  return 0;
  
}



void AliTPCcalibGainMult::BinLogX(const TH1 *h) {

  // Method for the correct logarithmic binning of histograms

  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
   
  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (int i = 1; i <= bins; i++) {
   newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;
  
  
}


void AliTPCcalibGainMult::UpdateGainMap() {
  //
  // read in the old gain map and scale it appropriately...
  //
  /*
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  //
  TFile jj("Run0_999999999_v1_s0.root");
  AliTPCCalPad * pad = AliCDBEntry->GetObject()->Clone();
  TFile hh("output.root");
  AliTPCcalibGainMult * gain = calibTracksGain;
  TH2D * histGainSec = gain->GetHistGainSector()->Projection(0,1);
  //
  TObjArray arr;
  histGainSec->FitSlicesY(0, 0, -1, 0, "QNR", &arr);
  TH1D * meanGainSec = arr->At(1);
  Double_t gainsIROC[36];
  Double_t gainsOROC[36];
  Double_t gains[72];
  //
  for(Int_t isec = 1; isec < meanGainSec->GetNbinsX() + 1; isec++) {
    cout << isec << " " << meanGainSec->GetXaxis()->GetBinCenter(isec) << " " <<meanGainSec->GetBinContent(isec) << endl;
    gains[isec-1] = meanGainSec->GetBinContent(isec);
    if (isec < 37) {
      gainsIROC[isec-1] = meanGainSec->GetBinContent(isec);
    } else {
      gainsOROC[isec - 36 -1] = meanGainSec->GetBinContent(isec);
    }
  }
  Double_t meanIroc = TMath::Mean(36, gainsIROC);
  Double_t meanOroc = TMath::Mean(36, gainsIROC);
  for(Int_t i = 0; i < 36; i++) gains[i] /= meanIroc;
  for(Int_t i = 36; i < 72; i++) gains[i] /= meanOroc;
  //
  for(Int_t i = 0; i < 72; i++) {
    AliTPCCalROC * chamber = pad->GetCalROC(i);
    chamber->Multiply(gains[i]);
    cout << i << " "<< chamber->GetMean() << endl;
  }
  //
  // update the OCDB
  //
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("AliTPCCalPad");
  metaData->SetResponsible("Alexander Kalweit");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("04-19-05"); //root version
  metaData->SetComment("New gain map for 1600V OROC gain increase and equalization. Valid for runs starting after technical stop beginning of September.");
  AliCDBId id1("TPC/Calib/GainFactorDedx", 131541, AliCDBRunRange::Infinity()); // important: new gain runs here..
  AliCDBStorage * gStorage = AliCDBManager::Instance()->GetStorage("local:///d/alice05/akalweit/projects/OCDBupdate/HighGain_2010-09-03/OCDB/");
  gStorage->Put(pad, id1, metaData);
  */
  
}

void AliTPCcalibGainMult::UpdateClusterParam() {
  //
  //
  //
  /*
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  //
  TFile ff("OldClsParam.root");
  AliTPCClusterParam * param = AliCDBEntry->GetObject()->Clone();
 
  TFile hh("output.root");
  AliTPCcalibGainMult * gain = calibGainMult;
  TH2D * histGainSec = gain->GetHistGainSector()->Projection(0,2);
  TObjArray arr;
  histGainSec->FitSlicesY(0, 0, -1, 0, "QNR", &arr);
  histGainSec->Draw("colz");
  TH1D * fitVal = arr.At(1);
  fitVal->Draw("same");
  param->GetQnormCorrMatrix()->Print();
  param->GetQnormCorrMatrix()(0,5) *= fitVal->GetBinContent(1)/fitVal->GetBinContent(1); // short pads Qtot
  param->GetQnormCorrMatrix()(1,5) *= fitVal->GetBinContent(2)/fitVal->GetBinContent(1); // med pads Qtot
  param->GetQnormCorrMatrix()(2,5) *= fitVal->GetBinContent(3)/fitVal->GetBinContent(1); // long pads Qtot
  //
  param->GetQnormCorrMatrix()(3,5) *= fitVal->GetBinContent(1)/fitVal->GetBinContent(1); // short pads Qmax -> scaling assumed
  param->GetQnormCorrMatrix()(4,5) *= fitVal->GetBinContent(2)/fitVal->GetBinContent(1); // med pads Qmax -> scaling assumed
  param->GetQnormCorrMatrix()(5,5) *= fitVal->GetBinContent(3)/fitVal->GetBinContent(1); // long pads Qmax -> scaling assumed
  //
  TFile jj("newClusterParam.root","RECREATE");
  param->Write();
  param->GetQnormCorrMatrix()->Print();
  //
  // update the OCDB
  // 
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("AliTPCClusterParam");
  metaData->SetResponsible("Alexander Kalweit");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("04-19-04"); //root version
  metaData->SetComment("1600V OROC / hard thres. / new algorithm");
  AliCDBId id1("TPC/Calib/ClusterParam", 0, AliCDBRunRange::Infinity()); // important: new gain runs here..
  AliCDBStorage * gStorage = AliCDBManager::Instance()->GetStorage("local:///lustre/alice/akalweit/baseline/CalibrationEntries/OldThres_NewAlgo_PP");
  gStorage->Put(param, id1, metaData);
  */
  

}

