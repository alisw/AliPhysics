
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
#include "TVectorF.h"
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
#include "TDatabasePDG.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliESDv0.h"
#include "AliESDkink.h"
#include "AliRecoParam.h"
#include "AliTracker.h"
#include "AliTPCTransform.h"
#include "AliTPCROC.h"
#include "TROOT.h"

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
   fHistGainMult(0),
   fPIDMatrix(0),
   fHistdEdxMap(0),
   fHistdEdxMax(0),
   fHistdEdxTot(0),
   fdEdxTree(0),
   fBBParam(0)
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
   fHistGainMult(0),
   fPIDMatrix(0),
   fHistdEdxMap(0),
   fHistdEdxMax(0),
   fHistdEdxTot(0),
   fdEdxTree(0),
   fBBParam(0)
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
  AliTPCcalibBase::BinLogX(fHistQA);
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
  Int_t binsPadEqual[6]    = { 400, 400,    4,   20,   50, 100};
  Double_t xminPadEqual[6] = { 0.0, 0.0, -0.5,    0, -250,   0}; 
  Double_t xmaxPadEqual[6] = { 2.0, 2.0,  3.5, 13000,  250,   3};
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
  //
  //                    dedx maps - bigger granulatity in phi -
  //                                to construct the dedx sector/phi map
  Int_t    binsGainMap[4]  = { 100,  90,             10,   6};
  Double_t xminGainMap[4]  = { 0.3,  -TMath::Pi(),    0,   0}; 
  Double_t xmaxGainMap[4]  = {   2,  TMath::Pi(),     1,   6};
  TString  axisNameMap[4]  = {"Q_Qexp","phi",      "1/Qexp","Pad Type"};
  TString  axisTitleMap[4] = {"Q/Q_{exp}","#phi (a.u.)","1/Q_{exp}","Pad Type"};
  //
  fHistdEdxMap = new THnSparseF("fHistdEdxMap","fHistdEdxMap", 4, binsGainMap, xminGainMap, xmaxGainMap); 
  for (Int_t iaxis=0; iaxis<4;iaxis++){
    fHistdEdxMap->GetAxis(iaxis)->SetName(axisNameMap[iaxis]);
    fHistdEdxMap->GetAxis(iaxis)->SetTitle(axisTitleMap[iaxis]);
  }
  //
  //
  //
  //                    dedx maps
  Int_t    binsGainMax[6]  = { 100,  10,  10,   10, 5,     3};
  Double_t xminGainMax[6]  = { 0.5,   0,   0,    0, 0,     0}; 
  Double_t xmaxGainMax[6]  = { 1.5,   1, 1.0,  1.0, 3000,  3};
  TString  axisNameMax[6]  = {"Q_Qexp","1/Qexp",  "phi","theta","mult", "Pad Type"};
  TString  axisTitleMax[6] = {"Q/Q_{exp}","1/Qexp", "#phi","#theta","mult","Pad Type"};
  //
  fHistdEdxMax = new THnSparseF("fHistdEdxMax","fHistdEdxMax", 6, binsGainMax, xminGainMax, xmaxGainMax); 
  fHistdEdxTot = new THnSparseF("fHistdEdxTot","fHistdEdxTot", 6, binsGainMax, xminGainMax, xmaxGainMax); 
  for (Int_t iaxis=0; iaxis<6;iaxis++){
    fHistdEdxMax->GetAxis(iaxis)->SetName(axisNameMax[iaxis]);
    fHistdEdxMax->GetAxis(iaxis)->SetTitle(axisTitleMax[iaxis]);
    fHistdEdxTot->GetAxis(iaxis)->SetName(axisNameMax[iaxis]);
    fHistdEdxTot->GetAxis(iaxis)->SetTitle(axisTitleMax[iaxis]);
  }
  //
  AliInfo("Non Default Constructor");  
}


AliTPCcalibGainMult::~AliTPCcalibGainMult(){
  //
  // Destructor
  //
  delete fHistNTracks;            //  histogram showing number of ESD tracks per event
  delete fHistClusterShape;       //  histogram to check the cluster shape
  delete fHistQA;                 //  dE/dx histogram showing the final spectrum
  //
  delete fHistGainSector;   //  histogram which shows MIP peak for each of the 3x36 sectors (pad region)
  delete fHistPadEqual;     //  histogram for the equalization of the gain in the different pad regions -> pass0
  delete fHistGainMult;     //  histogram which shows decrease of MIP signal as a function
  //
  delete fHistdEdxMap;
  delete fHistdEdxMax;
  delete fHistdEdxTot;
  delete fdEdxTree;
}



void AliTPCcalibGainMult::Process(AliESDEvent *event) {
  //
  // Main function of the class
  // 1. Select Identified  particles - for identified particles the flag in the PID matrix is stored
  //    1.a) ProcessV0s  - select Electron (gamma coversion) and pion canditates (from K0s) 
  //    1.b) ProcessTOF  - select - Proton, kaon and pions candidates
  //                       AS THE TOF not calibrated yet in Pass0 - we are calibrating the TOF T0 in this function    
  //    1.c) ProcessCosmic - select cosmic mumn candidates   - too few entries - not significant for the calibration
  //    1.d) ProcessKinks - select Kaon and pion candidates. From our experience (see Kink debug streamer), the angular cut for kink daughter is not sufficient - big contamination - delta rays, hadronic  interaction (proton)
  //          - NOT USED for the 
  //  
  // 2. Loop over tracks   
  //     2.a DumpTrack() - for identified particles dump the track and dEdx information into the tree (for later fitting)
  // 3. Actual fitting for the moment macro

  //
  // Criteria for the track selection
  //
  const Int_t kMinNCL=80;     // minimal number of cluster  - tracks accepted for the dedx calibration
  const Double_t kMaxEta=0.8; // maximal eta fo the track to be accepted
  const Double_t kMaxDCAR=10; // maximal DCA R of the track
  const Double_t kMaxDCAZ=5;  // maximal DCA Z of the track
  const Double_t kMIPPt=0.45; // MIP pt
  
  if (!event) {
    Printf("ERROR: ESD not available");
    return;
  }  
  fCurrentEvent=event  ;
  fMagF = event->GetMagneticField();
  Int_t ntracks=event->GetNumberOfTracks();  
  AliESDfriend *esdFriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!esdFriend) {
    //Printf("ERROR: esdFriend not available");
    delete fPIDMatrix;
    return;
  }
  if (!(esdFriend->TestSkipBit())) fPIDMatrix= new TMatrixD(ntracks,5);
  fHistNTracks->Fill(ntracks);
  ProcessCosmic(event);  // usually not enogh statistic

  if (esdFriend->TestSkipBit()) {
    return;
  }
  //
  ProcessV0s(event);   // 
  ProcessTOF(event);   //
  ProcessKinks(event); // not relyable
  DumpHPT(event);      // 
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

    if (ncls < kMinNCL) continue;     
    // exclude tracks which do not look like primaries or are simply too short or on wrong sectors
    if (TMath::Abs(trackIn->Eta()) > kMaxEta) continue;

    UInt_t status = track->GetStatus();
    if ((status&AliESDtrack::kTPCrefit)==0) continue;
    //if (track->GetNcls(0) < 3) continue; // ITS clusters
    Float_t dca[2], cov[3];
    track->GetImpactParameters(dca,cov);
    Float_t primVtxDCA = TMath::Sqrt(dca[0]*dca[0]);
    if (primVtxDCA > kMaxDCAR || primVtxDCA < 0.00001) continue;
    if (TMath::Abs(dca[1]) > kMaxDCAZ) continue;
    //
    //
    // fill Alexander QA histogram
    //
    if (primVtxDCA < 3 && track->GetNcls(0) > 3 && track->GetKinkIndex(0) == 0 && ncls > 100) fHistQA->Fill(meanP, track->GetTPCsignal(), 5);

    // Get seeds
    AliESDfriendTrack *friendTrack = esdFriend->GetTrack(i);
    if (!friendTrack) continue;
    TObject *calibObject;
    AliTPCseed *seed = 0;
    for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
      if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }    
    if (seed) DumpTrack(track, friendTrack,seed,i); // MI implementation for the identified particles
    //
    if (seed) { // seed the container with track parameters and the clusters
      // 
      const AliExternalTrackParam * trackOut = friendTrack->GetTPCOut();  // tack at the outer radius of TPC
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
      //      if (meanP > 0.4 && meanP < 0.55) {
      if ( TMath::Abs(meanP-kMIPPt)<0.05 ) {
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
      if ( TMath::Abs(meanP-kMIPPt)>0.05 ) continue;  // only MIP pions
      //if (meanP > 0.5 || meanP < 0.4) continue; // only MIP pions
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
	// AK. -  run Number To be removed - not needed 
	Double_t vecGainSec[5] = {vecMip[ipad], sector[ipad], ipad, runNumber, ncl[ipad]};
	if (isNotSplit[ipad]) fHistGainSector->Fill(vecGainSec);
      }
    }
   
  }    

  delete fPIDMatrix;
}  


void AliTPCcalibGainMult::MakeLookup(THnSparse * /*hist*/, Char_t * /*outputFile*/) {
  //
  // Not  yet implemented
  //
}


void AliTPCcalibGainMult::Analyze() {
  //
  // Not implemented
  //

  return;

}


Long64_t AliTPCcalibGainMult::Merge(TCollection *li) {
  //
  // merging of the component
  //

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
    if (cal->fHistdEdxMap){
      if (fHistdEdxMap) fHistdEdxMap->Add(cal->fHistdEdxMap);
    }
    if (cal->fHistdEdxMax){
      if (fHistdEdxMax) fHistdEdxMax->Add(cal->fHistdEdxMax);
    }
    if (cal->fHistdEdxTot){
      if (fHistdEdxTot) fHistdEdxTot->Add(cal->fHistdEdxTot);
    }
    // 
    // Originally we tireied to write the tree to the same file as other calibration components
    // We failed in merging => therefore this optio  was disabled
    //
    //    if (cal->fdEdxTree && cal->fdEdxTree->GetEntries()>0) {
    //       if (fdEdxTree) {
    // 	const Int_t kMax=100000;
    // 	Int_t entriesSum = (Int_t)fdEdxTree->GetEntries();
    // 	Int_t entriesCurrent = (Int_t)cal->fdEdxTree->GetEntries();
    // 	Int_t entriesCp=TMath::Min((Int_t) entriesCurrent*(kMax*entriesSum),entriesCurrent);
    // // 	cal->fdEdxTree->SetBranchStatus("track*",0);
    // // 	cal->fdEdxTree->SetBranchStatus("vertex*",0);
    // // 	cal->fdEdxTree->SetBranchStatus("tpcOut*",0);
    // // 	cal->fdEdxTree->SetBranchStatus("vec*",0);
    // // 	fdEdxTree->SetBranchStatus("track*",0);
    // // 	fdEdxTree->SetBranchStatus("vertex*",0);
    // // 	fdEdxTree->SetBranchStatus("tpcOut*",0);
    // // 	fdEdxTree->SetBranchStatus("vec*",0);
    // 	fdEdxTree->Print();
    // 	fdEdxTree->Dump();
    // 	fdEdxTree->GetEntry(0);
    // 	for (Int_t i=0; i<entriesCurrent; i++){
    // 	  cal->fdEdxTree->CopyAddresses(fdEdxTree);
    // 	  cal->fdEdxTree->GetEntry(i);
    // 	  fdEdxTree->Fill();
    // 	}		     
    // 	TObjArray *brarray =  cal->fdEdxTree->GetListOfBranches(); 
    // 	for (Int_t i=0; i<brarray->GetEntries(); i++) {TBranch * br = (TBranch *)brarray->At(i); br->SetAddress(0); }      
    //       }
    //       else{
    // 	fdEdxTree = (TTree*)(cal->fdEdxTree->Clone());
    // 	TObjArray *brarray =  fdEdxTree->GetListOfBranches(); 
    // 	for (Int_t i=0; i<brarray->GetEntries(); i++) {TBranch * br = (TBranch *)brarray->At(i); br->SetAddress(0);}      	
    //       }
    //}
    
  }
  
  return 0;
  
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


void AliTPCcalibGainMult::DumpTrack(AliESDtrack * track, AliESDfriendTrack *ftrack, AliTPCseed * seed, Int_t index){
  //
  // dump interesting tracks
  // 1. track at MIP region
  // 2. highly ionizing protons
  // pidType: 0 - unselected 
  //          1 - TOF
  //          2 - V0
  //          4 - Cosmic
  //          or of value
  //
  const Int_t    kMax=10000;
  const Int_t    kMinRows=80;
  const Double_t kDCAcut=30;
  //
  // Bethe Bloch paramerization
  //
  Double_t kp1= 0.0851148;
  Double_t kp2= 9.25771;
  Double_t kp3= 2.6558e-05;
  Double_t kp4= 2.32742;
  Double_t kp5= 1.83039;
  if (fBBParam){
    kp1=(*fBBParam)[0];
    kp2=(*fBBParam)[1];
    kp3=(*fBBParam)[2];
    kp4=(*fBBParam)[3];
    kp5=(*fBBParam)[4];
  }
  //
  //AliTPCROC *roc = AliTPCROC::Instance();
  TDatabasePDG *pdg = TDatabasePDG::Instance();

  Int_t nclITS   = track->GetNcls(0);
  Int_t ncl   = track->GetTPCncls();
  Double_t ncl21 = track->GetTPCClusterInfo(3,1);
  Double_t ncl20 = track->GetTPCClusterInfo(3,0);
  //
  if (!seed) return;
  if (ncl21<kMinRows) return;  
  static Int_t counter=0;
  static Int_t counterHPT=0;
  //
  static TH1F     *hisBB=new TH1F("hisBB","hisBB",20,0.1,1.00);  // bethe bloch histogram  = 
  //                                                                 used to cover more homegenously differnt dEdx regions
  static Double_t massPi = pdg->GetParticle("pi-")->Mass();      // 
  static Double_t massK  = pdg->GetParticle("K-")->Mass();
  static Double_t massP  = pdg->GetParticle("proton")->Mass();
  static Double_t massE  = pdg->GetParticle("e-")->Mass();
  static Double_t massMuon  = pdg->GetParticle("mu-")->Mass();
  static Double_t radius0= roc->GetPadRowRadiiLow(roc->GetNRows(0)/2);
  static Double_t radius1= roc->GetPadRowRadiiUp(30);
  static Double_t radius2= roc->GetPadRowRadiiUp(roc->GetNRows(36)-15);

  AliESDVertex *vertex= (AliESDVertex *)fCurrentEvent->GetPrimaryVertex();
  //
  // Estimate current MIP position - 
  //
  static Double_t mipArray[kMax];               // mip array
  static Int_t    counterMIP0=0;          
  static Double_t    medianMIP0=100000;         // current MIP median position - estimated after some mimnimum number of MIP tracks

  if (TMath::Abs(track->GetP()-0.5)<0.1&&track->GetTPCsignal()/medianMIP0<1.5){
    mipArray[counterMIP0%kMax]= track->GetTPCsignal();
    counterMIP0++;
    if (counterMIP0>10) medianMIP0=TMath::Median(counterMIP0%kMax, mipArray);
  }
  // the PID as defiend from the external sources
  //
  Int_t isElectron   =  TMath::Nint((*fPIDMatrix)(index,0));
  Int_t isMuon       =  TMath::Nint((*fPIDMatrix)(index,1));
  Int_t isPion       =  TMath::Nint((*fPIDMatrix)(index,2));
  Int_t isKaon       =  TMath::Nint((*fPIDMatrix)(index,3));
  Int_t isProton     =  TMath::Nint((*fPIDMatrix)(index,4));
  Float_t dca[2];
  track->GetImpactParameters(dca[0],dca[1]);
  //
  if ( (isMuon==0 && isElectron==0)  && (TMath::Sqrt(dca[0]*dca[0]+dca[1]*dca[1])>kDCAcut) ) return;
  Double_t normdEdx= track->GetTPCsignal()/(medianMIP0); // TPC signal normalized to the MIP
  //
  AliExternalTrackParam * trackIn  = (AliExternalTrackParam *)track->GetInnerParam();
  AliExternalTrackParam * trackOut = (AliExternalTrackParam *)track->GetOuterParam();
  AliExternalTrackParam * tpcOut   = (AliExternalTrackParam *)ftrack->GetTPCOut();
  if (!trackIn) return;
  if (!trackOut) return;
  if (!tpcOut) return;
  if (trackIn->GetZ()*trackOut->GetZ()<0) return;  // remove crossing tracks
  //
  // calculate local and global angle
  Int_t side = (trackIn->GetZ()>0)? 1:-1;
  Double_t tgl=trackIn->GetTgl();
  Double_t gangle[3]={0,0,0};
  Double_t langle[3]={0,0,0};
  Double_t length[3]={0,0,0};
  Double_t pxpypz[3]={0,0,0};
  Double_t bz=fMagF;
  trackIn->GetXYZAt(radius0,bz,pxpypz);            // get the global position  at the middle of the IROC
  gangle[0]=TMath::ATan2(pxpypz[1],pxpypz[0]);     // global angle IROC 
  trackIn->GetXYZAt(radius1,bz,pxpypz);            // get the global position at the middle of the OROC - medium pads      
  gangle[1]=TMath::ATan2(pxpypz[1],pxpypz[0]);     // global angle OROC
  trackOut->GetXYZAt(radius2,bz,pxpypz);           // get the global position at the middle of OROC - long pads
  gangle[2]=TMath::ATan2(pxpypz[1],pxpypz[0]);
  //
  trackIn->GetPxPyPzAt(radius0,bz,pxpypz);               //get momentum vector 
  langle[0]=TMath::ATan2(pxpypz[1],pxpypz[0])-gangle[0];  //local angle between padrow and track IROC  
  trackIn->GetPxPyPzAt(radius1,bz,pxpypz); 
  langle[1]=TMath::ATan2(pxpypz[1],pxpypz[0])-gangle[1];                                           
  trackOut->GetPxPyPzAt(radius2,bz,pxpypz);               //                                     OROC medium    
  langle[2]=TMath::ATan2(pxpypz[1],pxpypz[0])-gangle[2];
  for (Int_t i=0; i<3; i++){
    if (langle[i]>TMath::Pi())  langle[i]-=TMath::TwoPi();
    if (langle[i]<-TMath::Pi()) langle[i]+=TMath::TwoPi();
    length[i]=TMath::Sqrt(1+langle[i]*langle[i]+tgl*tgl);  // the tracklet length
  }
  //
  // Select the kaons and Protons which are "isolated" in TPC dedx curve
  // 
  //
  Double_t dedxP = AliExternalTrackParam::BetheBlochAleph(track->GetInnerParam()->GetP()/massP,kp1,kp2,kp3,kp4,kp5);
  Double_t dedxK = AliExternalTrackParam::BetheBlochAleph(track->GetInnerParam()->GetP()/massK,kp1,kp2,kp3,kp4,kp5);
  if (dedxP>2 || dedxK>2){
    if (track->GetP()<1.2 && normdEdx>1.8&&counterMIP0>10){ // not enough from TOF and V0s triggered by high dedx
      // signing the Proton  and kaon - using the "bitmask" bit 1 and 2 is dedicated for V0s and TOF selected       
      if ( TMath::Abs(normdEdx/dedxP-1)<0.3)  isProton+=4;
      if ( TMath::Abs(normdEdx/dedxK-1)<0.3)  isKaon+=4;
      if (normdEdx/dedxK>1.3) isProton+=8;
      if (normdEdx/dedxP<0.7) isKaon+=8;
    }
  }
  //
  //
  //
  Double_t mass = 0;  
  Bool_t isHighPt = ((TMath::Power(1/track->Pt(),4)*gRandom->Rndm())<0.005);  // rnadomly selected HPT tracks
  // there are selected for the QA of the obtained calibration
  Bool_t isMIP    =  TMath::Abs(track->GetInnerParam()->P()-0.4)<0.005&&(counter<kMax); //
  // REMINDER - it is not exactly MIP - we select the regtion where the Kaon and Electrons are well separated

  if (isElectron>0) mass = massE;
  if (isProton>0)   mass = massP;
  if (isKaon>0)     mass = massK;
  if (isMuon>0)     mass = massMuon;
  if (isPion>0)     mass = massPi;
  if (isHighPt)     mass = massPi;  //assign mass of pions
  if (isMIP&&track->GetTPCsignal()/medianMIP0<1.5)   mass = massPi;  //assign mass of pions
  if (mass==0)      return;
  //
  // calculate expected dEdx
  Double_t dedxDef= 0;
  Double_t dedxDefPion= 0,dedxDefProton=0, dedxDefKaon=0;
  Double_t pin=trackIn->GetP();
  Double_t pout=trackOut->GetP();
  Double_t p=(pin+pout)*0.5;  // momenta as the mean between tpc momenta at inner and outer wall of the TPC
  if (mass>0) dedxDef = AliExternalTrackParam::BetheBlochAleph(p/mass,kp1,kp2,kp3,kp4,kp5); 
  dedxDefPion = AliExternalTrackParam::BetheBlochAleph(p/massPi,kp1,kp2,kp3,kp4,kp5); 
  dedxDefProton = AliExternalTrackParam::BetheBlochAleph(p/massP,kp1,kp2,kp3,kp4,kp5); 
  dedxDefKaon = AliExternalTrackParam::BetheBlochAleph(p/massK,kp1,kp2,kp3,kp4,kp5); 
  //
  // dEdx Truncated mean vectros with differnt tuncation 
  // 11 truncations array -  0-10  - 0~50%  11=100%
  // 3 Regions            -  IROC,OROC0, OROC1
  // 2 Q                  -  total charge and max charge
  // Log                  -  Logarithmic mean used
  // Up/Dwon              -  Upper half or lower half of truncation used
  // RMS                  -  rms of the distribction (otherwise truncated mean)
  // M2 suffix            -  second moment ( truncated) 
  TVectorF truncUp(11);
  TVectorF truncDown(11);
  TVectorF vecAllMax(11);
  TVectorF vecIROCMax(11);
  TVectorF vecOROCMax(11);
  TVectorF vecOROC0Max(11);
  TVectorF vecOROC1Max(11);
  //
  TVectorF vecAllTot(11);
  TVectorF vecIROCTot(11);
  TVectorF vecOROCTot(11);
  TVectorF vecOROC0Tot(11);
  TVectorF vecOROC1Tot(11);
  //
  TVectorF vecAllTotLog(11);
  TVectorF vecIROCTotLog(11);
  TVectorF vecOROCTotLog(11);
  TVectorF vecOROC0TotLog(11);
  TVectorF vecOROC1TotLog(11);
  //
  TVectorF vecAllTotUp(11);
  TVectorF vecIROCTotUp(11);
  TVectorF vecOROCTotUp(11);
  TVectorF vecOROC0TotUp(11);
  TVectorF vecOROC1TotUp(11);
  //
  TVectorF vecAllTotDown(11);
  TVectorF vecIROCTotDown(11);
  TVectorF vecOROCTotDown(11);
  TVectorF vecOROC0TotDown(11);
  TVectorF vecOROC1TotDown(11);

  TVectorF vecAllTotRMS(11);
  TVectorF vecIROCTotRMS(11);
  TVectorF vecOROCTotRMS(11);
  TVectorF vecOROC0TotRMS(11);
  TVectorF vecOROC1TotRMS(11);
  //
  TVectorF vecAllTotM2(11);
  TVectorF vecIROCTotM2(11);
  TVectorF vecOROCTotM2(11);
  TVectorF vecOROC0TotM2(11);
  TVectorF vecOROC1TotM2(11);
  //
  TVectorF vecAllTotMS(11);
  TVectorF vecIROCTotMS(11);
  TVectorF vecOROCTotMS(11);
  TVectorF vecOROC0TotMS(11);
  TVectorF vecOROC1TotMS(11);
  //
  // Differnt number of clusters definitions - in separate regions of the TPC
  // 20  -    ratio - found/findabel
  // 21  -    number of clusters used for given dEdx calculation
  //
  // suffix - 3 or 4 -  number of padrows before and after given row to define findable row
  //
  Double_t ncl20All  = seed->CookdEdxAnalytical(0.0,1, 1 ,0,159,3);
  Double_t ncl20IROC = seed->CookdEdxAnalytical(0.,1, 1 ,0,63,3);
  Double_t ncl20OROC = seed->CookdEdxAnalytical(0.,1, 1 ,64,159,3);
  Double_t ncl20OROC0= seed->CookdEdxAnalytical(0.,1, 1 ,64,128,3);
  Double_t ncl20OROC1= seed->CookdEdxAnalytical(0.,1, 1 ,129,159,3);
  //
  Double_t ncl20All4  = seed->CookdEdxAnalytical(0.0,1, 1 ,0,159,3,4);
  Double_t ncl20IROC4 = seed->CookdEdxAnalytical(0.,1, 1 ,0,63,3,4);
  Double_t ncl20OROC4 = seed->CookdEdxAnalytical(0.,1, 1 ,64,159,3,4);
  Double_t ncl20OROC04= seed->CookdEdxAnalytical(0.,1, 1 ,64,128,3,4);
  Double_t ncl20OROC14= seed->CookdEdxAnalytical(0.,1, 1 ,129,159,3,4);
  //
  Double_t ncl20All3  = seed->CookdEdxAnalytical(0.0,1, 1 ,0,159,3,3);
  Double_t ncl20IROC3 = seed->CookdEdxAnalytical(0.,1, 1 ,0,63,3,3);
  Double_t ncl20OROC3 = seed->CookdEdxAnalytical(0.,1, 1 ,64,159,3,3);
  Double_t ncl20OROC03= seed->CookdEdxAnalytical(0.,1, 1 ,64,128,3,3);
  Double_t ncl20OROC13= seed->CookdEdxAnalytical(0.,1, 1 ,129,159,3,3);
  //
  Double_t ncl21All  = seed->CookdEdxAnalytical(0.0,1, 1 ,0,159,2);
  Double_t ncl21IROC = seed->CookdEdxAnalytical(0.,1, 1 ,0,63,2);
  Double_t ncl21OROC = seed->CookdEdxAnalytical(0.,1, 1 ,64,159,2);
  Double_t ncl21OROC0= seed->CookdEdxAnalytical(0.,1, 1 ,64,128,2);
  Double_t ncl21OROC1= seed->CookdEdxAnalytical(0.,1, 1 ,129,159,2);
  // calculate truncated dEdx - mean rms M2 ... 
  Int_t ifrac=0;
  for (Int_t ifracDown=0; ifracDown<1; ifracDown++){
    for (Int_t ifracUp=0; ifracUp<11; ifracUp++){
      Double_t fracDown = 0.0+Double_t(ifracDown)*0.05;
      Double_t fracUp = 0.5+Double_t(ifracUp)*0.05;
      vecAllMax[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 1 ,0,159,0);
      vecIROCMax[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 1 ,0,63,0);
      vecOROCMax[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 1 ,64,159,0);
      vecOROC0Max[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 1 ,64,128,0);
      vecOROC1Max[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 1 ,129,159,0);
      //
      vecAllTot[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,0,159,0);
      vecIROCTot[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,0,63,0);
      vecOROCTot[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,64,159,0);
      vecOROC0Tot[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,64,128,0);
      vecOROC1Tot[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,129,159,0);
      //
      vecAllTotLog[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,0,159,0,2,1);
      vecIROCTotLog[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,0,63,0,2,1);
      vecOROCTotLog[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,64,159,0,2,1);
      vecOROC0TotLog[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,64,128,0,2,1);
      vecOROC1TotLog[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,129,159,0,2,1);
      //
      vecAllTotUp[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,0,159,4,2,1);
      vecIROCTotUp[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,0,63,4,2,1);
      vecOROCTotUp[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,64,159,4,2,1);
      vecOROC0TotUp[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,64,128,4,2,1);
      vecOROC1TotUp[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,129,159,4,2,1);
      //
      vecAllTotDown[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,0,159,5,2,1);
      vecIROCTotDown[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,0,63,5,2,1);
      vecOROCTotDown[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,64,159,5,2,1);
      vecOROC0TotDown[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,64,128,5,2,1);
      vecOROC1TotDown[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,129,159,5,2,1);
      //
      vecAllTotRMS[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,0,159,1,2,0);
      vecIROCTotRMS[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,0,63,1,2,0);
      vecOROCTotRMS[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,64,159,1,2,0);
      vecOROC0TotRMS[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,64,128,1,2,0);
      vecOROC1TotRMS[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,129,159,1,2,0);
      //
      vecAllTotM2[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,0,159,6,2,1);
      vecIROCTotM2[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,0,63,6,2,1);
      vecOROCTotM2[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,64,159,6,2,1);
      vecOROC0TotM2[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,64,128,6,2,1);
      vecOROC1TotM2[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,129,159,6,2,1);
      //
      vecAllTotMS[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,0,159,8,2,1);
      vecIROCTotMS[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,0,63,8,2,1);
      vecOROCTotMS[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,64,159,8,2,1);
      vecOROC0TotMS[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,64,128,8,2,1);
      vecOROC1TotMS[ifrac]= seed->CookdEdxAnalytical(fracDown,fracUp, 0 ,129,159,8,2,1);
      truncUp[ifrac]=fracUp;
      truncDown[ifrac]=fracDown;
      ifrac++;
    }
  }
  //
  // Fill histograms
  //
  if ((isKaon||isProton||isPion||isElectron||isMIP||isMuon&&(!isHighPt)) && dedxDef>0) {
    //
    Int_t ncont = vertex->GetNContributors();
    for (Int_t ipad=0; ipad<3; ipad++){
      // histogram with enahanced phi granularity - to make gain phi maps
      Double_t xxx[4]={0,gangle[ipad],1./dedxDef,ipad*2+((side>0)?0:1)};
      Double_t nclR=0;
      if (ipad==0)  {xxx[0]=vecIROCTot[4]/medianMIP0; nclR=ncl21IROC/63.;}
      if (ipad==1)  {xxx[0]=vecOROC0Tot[4]/medianMIP0;nclR=ncl21OROC0/63.;}
      if (ipad==2)  {xxx[0]=vecOROC1Tot[4]/medianMIP0;nclR=ncl21OROC1/32.;}
      xxx[0]/=dedxDef;
      if (xxx[0]>0) xxx[0]=1./xxx[0];
      if (TMath::Abs(langle[ipad])<0.25&&nclR>0.4)  fHistdEdxMap->Fill(xxx);
    }
    for (Int_t ipad=0; ipad<3; ipad++){
      //
      // this are histogram to define  overall main gain correction
      // Maybe dead end - we can not put all info which can be used into the THnSparse
      // It is keeped there for educational point of view
      //
      Double_t xxx[6]={0,1./dedxDef, TMath::Abs(langle[ipad]), TMath::Abs(tgl), ncont, ipad };
      if (ipad==0)  {xxx[0]=vecIROCTot[4]/medianMIP0;}
      if (ipad==1)  {xxx[0]=vecOROC0Tot[4]/medianMIP0;}
      if (ipad==2)  {xxx[0]=vecOROC1Tot[4]/medianMIP0;}
      xxx[0]/=dedxDef;
      if (xxx[0]>0) xxx[0]=1./xxx[0];
      if (xxx[0]>0) fHistdEdxTot->Fill(xxx);
      if (ipad==0)  {xxx[0]=vecIROCMax[4]/medianMIP0;}
      if (ipad==1)  {xxx[0]=vecOROC0Max[4]/medianMIP0;}
      if (ipad==2)  {xxx[0]=vecOROC1Max[4]/medianMIP0;}
      xxx[0]=dedxDef;
      if (xxx[0]>0) xxx[0]=1./xxx[0];
      fHistdEdxMax->Fill(xxx);
    }
  }  
  //
  // Downscale  selected tracks before filling the tree
  //
  Bool_t isSelected = kFALSE;  
  //
  if (isKaon||isProton||isPion||isElectron||isMIP||isMuon) isSelected=kTRUE;
  isHighPt = kFALSE;
  if (!isSelected) isHighPt = ((TMath::Power(1/track->Pt(),4)*gRandom->Rndm())<0.005);  
  if (counter>kMax && ((1/track->Pt()*gRandom->Rndm())>kMax/counter)) return; 
  isSelected|=isHighPt;
  //
  //
  //
  // Equalize statistic in BB bins - special care about pions
  Int_t entriesBB = (Int_t)hisBB->GetEntries();
  if ((isElectron==0 &&isMuon==0 && p<2.) && entriesBB>20 &&dedxDef>0.01){
    Int_t bin = hisBB->GetXaxis()->FindBin(1./dedxDef);
    Double_t cont = hisBB->GetBinContent(bin);
    Double_t mean =(entriesBB)/20.;
    if ((isPion>0)  && gRandom->Rndm()*cont > 0.1*mean) return;
    if ((isPion==0) && gRandom->Rndm()*cont > 0.25*mean) return;
  }  
  if (!isSelected) return;
  if (dedxDef>0.01) hisBB->Fill(1./dedxDef);  
  //
  if (isHighPt) counterHPT++;
  counter++;  
  //
  TTreeSRedirector * pcstream =  GetDebugStreamer();
  Double_t ptrel0 = AliTPCcalibDB::GetPTRelative(fTime,fRun,0);
  Double_t ptrel1 = AliTPCcalibDB::GetPTRelative(fTime,fRun,1);
  Int_t sectorIn   = Int_t(18+9*(trackIn->GetAlpha()/TMath::Pi()))%18;
  Int_t sectorOut  = Int_t(18+9*(trackOut->GetAlpha()/TMath::Pi()))%18;
  //
  if (pcstream){
    (*pcstream)<<"dump"<<
      "vertex.="<<vertex<<
      "bz="<<fMagF<<
      "ptrel0="<<ptrel0<<
      "ptrel1="<<ptrel1<<
      "sectorIn="<<sectorIn<<
      "sectorOut="<<sectorOut<<
      "side="<<side<<
      // pid type
      "isMuon="<<isMuon<<
      "isProton="<<isProton<<
      "isKaon="<<isKaon<<
      "isPion="<<isPion<<
      "isElectron="<<isElectron<<
      "isMIP="<<isMIP<<
      "isHighPt="<<isHighPt<<
      "mass="<<mass<<
      "dedxDef="<<dedxDef<<
      "dedxDefPion="<<dedxDefPion<<
      "dedxDefKaon="<<dedxDefKaon<<
      "dedxDefProton="<<dedxDefProton<<
      //
      "nclITS="<<nclITS<<
      "ncl="<<ncl<<
      "ncl21="<<ncl21<<
      "ncl20="<<ncl20<<
      //
      "ncl20All="<<ncl20All<<
      "ncl20OROC="<<ncl20OROC<<
      "ncl20IROC="<<ncl20IROC<<
      "ncl20OROC0="<<ncl20OROC0<<
      "ncl20OROC1="<<ncl20OROC1<<
      //
      "ncl20All4="<<ncl20All4<<
      "ncl20OROC4="<<ncl20OROC4<<
      "ncl20IROC4="<<ncl20IROC4<<
      "ncl20OROC04="<<ncl20OROC04<<
      "ncl20OROC14="<<ncl20OROC14<<
      //
      "ncl20All3="<<ncl20All3<<
      "ncl20OROC3="<<ncl20OROC3<<
      "ncl20IROC3="<<ncl20IROC3<<
      "ncl20OROC03="<<ncl20OROC03<<
      "ncl20OROC13="<<ncl20OROC13<<
      //
      "ncl21All="<<ncl21All<<
      "ncl21OROC="<<ncl21OROC<<
      "ncl21IROC="<<ncl21IROC<<
      "ncl21OROC0="<<ncl21OROC0<<
      "ncl21OROC1="<<ncl21OROC1<<  
      //track properties
      "langle0="<<langle[0]<<
      "langle1="<<langle[1]<<
      "langle2="<<langle[2]<<
      "gangle0="<<gangle[0]<<   //global angle phi IROC 
      "gangle1="<<gangle[1]<<   //                 OROC medium 
      "gangle2="<<gangle[2]<<   //                 OROC long
      "L0="<<length[0]<<
      "L1="<<length[1]<<
      "L2="<<length[2]<<
      "p="<<p<<
      "pin="<<pin<<
      "pout="<<pout<<
      "tgl="<<tgl<<
      "track.="<<track<<
      "trackIn.="<<trackIn<<
      "trackOut.="<<trackOut<<
      "tpcOut.="<<tpcOut<<
      "medianMIP0="<<medianMIP0<<    // median MIP position as estimated from the array of (not only) "MIPS"
      //dedx 
      "truncUp.="<<&truncUp<<
      "truncDown.="<<&truncDown<<
      "vecAllMax.="<<&vecAllMax<<
      "vecIROCMax.="<<&vecIROCMax<<
      "vecOROCMax.="<<&vecOROCMax<<
      "vecOROC0Max.="<<&vecOROC0Max<<
      "vecOROC1Max.="<<&vecOROC1Max<<
      //
      "vecAllTot.="<<&vecAllTot<<
      "vecIROCTot.="<<&vecIROCTot<<
      "vecOROCTot.="<<&vecOROCTot<<
      "vecOROC0Tot.="<<&vecOROC0Tot<<
      "vecOROC1Tot.="<<&vecOROC1Tot<<
      //
      "vecAllTotLog.="<<&vecAllTotLog<<
      "vecIROCTotLog.="<<&vecIROCTotLog<<
      "vecOROCTotLog.="<<&vecOROCTotLog<<
      "vecOROC0TotLog.="<<&vecOROC0TotLog<<
      "vecOROC1TotLog.="<<&vecOROC1TotLog<<
      //
      "vecAllTotUp.="<<&vecAllTotUp<<
      "vecIROCTotUp.="<<&vecIROCTotUp<<
      "vecOROCTotUp.="<<&vecOROCTotUp<<
      "vecOROC0TotUp.="<<&vecOROC0TotUp<<
      "vecOROC1TotUp.="<<&vecOROC1TotUp<<
      //
      "vecAllTotDown.="<<&vecAllTotDown<<
      "vecIROCTotDown.="<<&vecIROCTotDown<<
      "vecOROCTotDown.="<<&vecOROCTotDown<<
      "vecOROC0TotDown.="<<&vecOROC0TotDown<<
      "vecOROC1TotDown.="<<&vecOROC1TotDown<<
      //
      "vecAllTotRMS.="<<&vecAllTotRMS<<
      "vecIROCTotRMS.="<<&vecIROCTotRMS<<
      "vecOROCTotRMS.="<<&vecOROCTotRMS<<
      "vecOROC0TotRMS.="<<&vecOROC0TotRMS<<
      "vecOROC1TotRMS.="<<&vecOROC1TotRMS<<
      //
      "vecAllTotM2.="<<&vecAllTotM2<<
      "vecIROCTotM2.="<<&vecIROCTotM2<<
      "vecOROCTotM2.="<<&vecOROCTotM2<<
      "vecOROC0TotM2.="<<&vecOROC0TotM2<<
      "vecOROC1TotM2.="<<&vecOROC1TotM2<<
      //
      "vecAllTotMS.="<<&vecAllTotMS<<
      "vecIROCTotMS.="<<&vecIROCTotMS<<
      "vecOROCTotMS.="<<&vecOROCTotMS<<
      "vecOROC0TotMS.="<<&vecOROC0TotMS<<
      "vecOROC1TotMS.="<<&vecOROC1TotMS<<
      "\n";
  }
}




void AliTPCcalibGainMult::ProcessV0s(AliESDEvent * event){
  //
  // Select the K0s and gamma  - and sign daughter products 
  //  
  TTreeSRedirector * pcstream =  GetDebugStreamer();
  AliKFParticle::SetField(event->GetMagneticField()); 
  AliESDfriend *esdFriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!esdFriend) {
    //Printf("ERROR: esdFriend not available");
   return;
  }
  if (esdFriend->TestSkipBit()) return;
  //
  // 
  TDatabasePDG *pdg = TDatabasePDG::Instance();  
  const Double_t kChi2Cut=5;
  const Double_t kMinR=2;
  const Int_t    kMinNcl=110;
  const Double_t kMinREl=5;
  const Double_t kMaxREl=70;
  //
  Int_t nv0 = event->GetNumberOfV0s(); 
  AliESDVertex *vertex= (AliESDVertex *)event->GetPrimaryVertex();
  AliKFVertex kfvertex=*vertex;
  //
  for (Int_t iv0=0;iv0<nv0;iv0++){
    AliESDv0 *v0 = event->GetV0(iv0);
    if (!v0) continue;
    if (v0->GetOnFlyStatus()<0.5) continue;
    if (v0->GetPindex()<0) continue;
    if (v0->GetNindex()<0) continue;
    if (TMath::Max(v0->GetPindex(), v0->GetNindex())>event->GetNumberOfTracks()) continue;
    //
    //   
    AliExternalTrackParam pp=(v0->GetParamP()->GetSign()>0) ? (*(v0->GetParamP())):(*(v0->GetParamN()));
    AliExternalTrackParam pn=(v0->GetParamP()->GetSign()>0) ? (*(v0->GetParamN())):(*(v0->GetParamP()));
    AliKFParticle kfp1( pp, 211 );
    AliKFParticle kfp2( pn, -211 );
    //
    AliKFParticle *v0KFK0 = new AliKFParticle(kfp1,kfp2);
    AliKFParticle *v0KFK0CV = new AliKFParticle(*v0KFK0);
    v0KFK0CV->SetProductionVertex(kfvertex);
    v0KFK0CV->TransportToProductionVertex();
    Double_t chi2K0 = v0KFK0CV->GetChi2();
    if (chi2K0>kChi2Cut) continue;
    if (v0->GetRr()<kMinR) continue;
    Bool_t isOKC=TMath::Max(v0->GetCausalityP()[0],v0->GetCausalityP()[1])<0.7&&TMath::Min(v0->GetCausalityP()[2],v0->GetCausalityP()[3])>0.2;
    //
    Double_t effMass22=v0->GetEffMass(2,2);
    Double_t effMass42=v0->GetEffMass(4,2);
    Double_t effMass24=v0->GetEffMass(2,4);
    Double_t effMass00=v0->GetEffMass(0,0);
    AliKFParticle *v0KFK0CVM = new AliKFParticle(*v0KFK0CV);
    v0KFK0CVM->SetMassConstraint(pdg->GetParticle("K_S0")->Mass());
    Bool_t isV0= kFALSE;
    //    
    Double_t d22   = TMath::Abs(effMass22-pdg->GetParticle("K_S0")->Mass());
    Double_t d42   = TMath::Abs(effMass42-pdg->GetParticle("Lambda0")->Mass());
    Double_t d24   = TMath::Abs(effMass24-pdg->GetParticle("Lambda0")->Mass());
    Double_t d00   = TMath::Abs(effMass00);
    //
    Bool_t isKaon      = d22<0.01 && d22< 0.3 * TMath::Min(TMath::Min(d42,d24),d00);
    Bool_t isLambda    = d42<0.01 && d42< 0.3 * TMath::Min(TMath::Min(d22,d24),d00);
    Bool_t isAntiLambda= d24<0.01 && d24< 0.3 * TMath::Min(TMath::Min(d22,d42),d00);
    Bool_t isGamma     = d00<0.02 && d00< 0.3 * TMath::Min(TMath::Min(d42,d24),d22);
    //
    if (isGamma  &&  (isKaon||isLambda||isAntiLambda)) continue;
    if (isLambda &&  (isKaon||isGamma||isAntiLambda)) continue;
    if (isKaon   &&  (isLambda||isGamma||isAntiLambda)) continue;    
    Double_t sign= v0->GetParamP()->GetSign()* v0->GetParamN()->GetSign();
    if (sign>0) continue;
    isV0=isKaon||isLambda||isAntiLambda||isGamma;
    if (!(isV0)) continue;
    if (isGamma&&v0->GetRr()<kMinREl) continue;
    if (isGamma&&v0->GetRr()>kMaxREl) continue;
    if (!isOKC) continue;
    //
    Int_t pindex = (v0->GetParamP()->GetSign()>0) ? v0->GetPindex() : v0->GetNindex();
    Int_t nindex = (v0->GetParamP()->GetSign()>0) ? v0->GetNindex() : v0->GetPindex();
    AliESDtrack * trackP = event->GetTrack(pindex);
    AliESDtrack * trackN = event->GetTrack(nindex);
    if (!trackN) continue;
    if (!trackP) continue;
    Int_t nclP= (Int_t)trackP->GetTPCClusterInfo(2,1);
    Int_t nclN= (Int_t)trackN->GetTPCClusterInfo(2,1);
    if (TMath::Min(nclP,nclN)<kMinNcl) continue;
    Double_t eta = TMath::Max(TMath::Abs(trackP->Eta()), TMath::Abs(trackN->Eta()));
    if (TMath::Abs(eta)>1) continue;
    //
    //
    AliESDfriendTrack *friendTrackP = esdFriend->GetTrack(pindex);
    AliESDfriendTrack *friendTrackN = esdFriend->GetTrack(nindex);
    if (!friendTrackP) continue;
    if (!friendTrackN) continue;
    TObject *calibObject;
    AliTPCseed *seedP = 0;
    AliTPCseed *seedN = 0;
    for (Int_t l=0;(calibObject=friendTrackP->GetCalibObject(l));++l) {
      if ((seedP=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }    
    for (Int_t l=0;(calibObject=friendTrackN->GetCalibObject(l));++l) {
      if ((seedN=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }   
    if (isGamma){
      if ( TMath::Abs((trackP->GetTPCsignal()/(trackN->GetTPCsignal()+0.0001)-1)>0.3)) continue;
    }
    if (isGamma)   (*fPIDMatrix)(pindex, 0)+=2;
    if (isGamma)   (*fPIDMatrix)(nindex, 0)+=2;
    //
    if (isKaon)    (*fPIDMatrix)(pindex, 2)+=2;
    if (isKaon)    (*fPIDMatrix)(nindex, 2)+=2;
    //
    //
    if (pcstream){
      (*pcstream)<<"v0s"<<
	"isGamma="<<isGamma<<
	"isKaon="<<isKaon<<
	"isLambda="<<isLambda<<
	"isAntiLambda="<<isAntiLambda<<
	"chi2="<<chi2K0<<
	"trackP.="<<trackP<<
	"trackN.="<<trackN<<
	"v0.="<<v0<<
	"\n";
    }
  }
}




void AliTPCcalibGainMult::ProcessCosmic(const AliESDEvent * event) {
  //
  // Find cosmic pairs trigger by random trigger
  // 
  // 
  AliTPCTransform *transform = AliTPCcalibDB::Instance()->GetTransform() ;
  AliTPCParam     *param     = AliTPCcalibDB::Instance()->GetParameters();

  AliESDVertex *vertexSPD =  (AliESDVertex *)event->GetPrimaryVertexSPD();
  AliESDVertex *vertexTPC =  (AliESDVertex *)event->GetPrimaryVertexTPC(); 
  AliESDfriend *esdFriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  const Double_t kMinPt=4;
  const Double_t kMinPtMax=0.8;
  const Double_t kMinNcl=159*0.5;
  const Double_t kMaxDelta[5]={2,600,0.02,0.02,0.1};
  Int_t ntracks=event->GetNumberOfTracks(); 
  const Double_t kMaxImpact=80;
  //  Float_t dcaTPC[2]={0,0};
  // Float_t covTPC[3]={0,0,0};

  UInt_t specie = event->GetEventSpecie();  // skip laser events
  if (specie==AliRecoParam::kCalib) return;
  

  for (Int_t itrack0=0;itrack0<ntracks;itrack0++) {
    AliESDtrack *track0 = event->GetTrack(itrack0);
    if (!track0) continue;
    if (!track0->IsOn(AliESDtrack::kTPCrefit)) continue;

    if (TMath::Abs(AliTracker::GetBz())>1&&track0->Pt()<kMinPt) continue;
    if (track0->GetTPCncls()<kMinNcl) continue;
    if (TMath::Abs(track0->GetY())<2*kMaxDelta[0]) continue; 
    if (TMath::Abs(track0->GetY())>kMaxImpact) continue; 
    if (track0->GetKinkIndex(0)>0) continue;
    const Double_t * par0=track0->GetParameter(); //track param at rhe DCA
    //rm primaries
    //
    for (Int_t itrack1=itrack0+1;itrack1<ntracks;itrack1++) {
      AliESDtrack *track1 = event->GetTrack(itrack1);
      if (!track1) continue;  
      if (!track1->IsOn(AliESDtrack::kTPCrefit)) continue;
      if (track1->GetKinkIndex(0)>0) continue;
      if (TMath::Abs(AliTracker::GetBz())>1&&track1->Pt()<kMinPt) continue;
      if (track1->GetTPCncls()<kMinNcl) continue;
      if (TMath::Abs(AliTracker::GetBz())>1&&TMath::Max(track1->Pt(), track0->Pt())<kMinPtMax) continue;
      if (TMath::Abs(track1->GetY())<2*kMaxDelta[0]) continue;
      if (TMath::Abs(track1->GetY())>kMaxImpact) continue; 
      //
      const Double_t* par1=track1->GetParameter(); //track param at rhe DCA
      //
      Bool_t isPair=kTRUE;
      for (Int_t ipar=0; ipar<5; ipar++){
	if (ipar==4&&TMath::Abs(AliTracker::GetBz())<1) continue; // 1/pt not defined for B field off
	if (TMath::Abs(TMath::Abs(par0[ipar])-TMath::Abs(par1[ipar]))>kMaxDelta[ipar]) isPair=kFALSE;
      }
      if (!isPair) continue;
      if (TMath::Abs(TMath::Abs(track0->GetAlpha()-track1->GetAlpha())-TMath::Pi())>kMaxDelta[2]) isPair=kFALSE;
      //delta with correct sign
      if  (TMath::Abs(par0[0]+par1[0])>kMaxDelta[0]) isPair=kFALSE; //delta y   opposite sign
      if  (TMath::Abs(par0[3]+par1[3])>kMaxDelta[3]) isPair=kFALSE; //delta tgl opposite sign
      if  (TMath::Abs(AliTracker::GetBz())>1 && TMath::Abs(par0[4]+par1[4])>kMaxDelta[4]) isPair=kFALSE; //delta 1/pt opposite sign
      if (!isPair) continue;
      TString filename(AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile()->GetName());
      Int_t eventNumber = event->GetEventNumberInFile(); 
      Bool_t hasFriend=(esdFriend) ? (esdFriend->GetTrack(itrack0)!=0):0; 
      Bool_t hasITS=(track0->GetNcls(0)+track1->GetNcls(0)>4);
      printf("DUMPHPTCosmic:%s|%f|%d|%d|%d\n",filename.Data(),(TMath::Min(track0->Pt(),track1->Pt())), eventNumber,hasFriend,hasITS);
      //
      //       
      TTreeSRedirector * pcstream =  GetDebugStreamer();
      Int_t ntracksSPD = vertexSPD->GetNContributors();
      Int_t ntracksTPC = vertexTPC->GetNContributors();
      //
      if (pcstream){
	(*pcstream)<<"cosmicPairsAll"<<
	  "run="<<fRun<<              //  run number
	  "event="<<fEvent<<          //  event number
	  "time="<<fTime<<            //  time stamp of event
	  "trigger="<<fTrigger<<      //  trigger
	  "triggerClass="<<&fTriggerClass<<      //  trigger
	  "bz="<<fMagF<<             //  magnetic field
	  //
	  "nSPD="<<ntracksSPD<<
	  "nTPC="<<ntracksTPC<<
	  "vSPD.="<<vertexSPD<<         //primary vertex -SPD
	  "vTPC.="<<vertexTPC<<         //primary vertex -TPC
	  "t0.="<<track0<<              //track0
	  "t1.="<<track1<<              //track1
	  "\n";      
      }
      //
      AliESDfriendTrack *friendTrack0 = esdFriend->GetTrack(itrack0);
      if (!friendTrack0) continue;
      AliESDfriendTrack *friendTrack1 = esdFriend->GetTrack(itrack1);
      if (!friendTrack1) continue;
      TObject *calibObject;
      AliTPCseed *seed0 = 0;   
      AliTPCseed *seed1 = 0;
      //
      for (Int_t l=0;(calibObject=friendTrack0->GetCalibObject(l));++l) {
	if ((seed0=dynamic_cast<AliTPCseed*>(calibObject))) break;
      }
      for (Int_t l=0;(calibObject=friendTrack1->GetCalibObject(l));++l) {
	if ((seed1=dynamic_cast<AliTPCseed*>(calibObject))) break;
      }
      //
      if (pcstream){
	(*pcstream)<<"cosmicPairs"<<
	  "run="<<fRun<<              //  run number
	  "event="<<fEvent<<          //  event number
	  "time="<<fTime<<            //  time stamp of event
	  "trigger="<<fTrigger<<      //  trigger
	  "triggerClass="<<&fTriggerClass<<      //  trigger
	  "bz="<<fMagF<<             //  magnetic field
	  //
	  "nSPD="<<ntracksSPD<<
	  "nTPC="<<ntracksTPC<<
	  "vSPD.="<<vertexSPD<<         //primary vertex -SPD
	  "vTPC.="<<vertexTPC<<         //primary vertex -TPC
	  "t0.="<<track0<<              //track0
	  "t1.="<<track1<<              //track1
 	  "ft0.="<<friendTrack0<<       //track0
 	  "ft1.="<<friendTrack1<<       //track1
 	  "s0.="<<seed0<<               //track0
 	  "s1.="<<seed1<<               //track1
	  "\n";      
      }
      if (!seed0) continue;
      if (!seed1) continue;
      Int_t nclA0=0, nclC0=0;     // number of clusters
      Int_t nclA1=0, nclC1=0;     // number of clusters
      //      
      for (Int_t irow=0; irow<159; irow++){
	AliTPCclusterMI *cluster0=seed0->GetClusterPointer(irow);
	AliTPCclusterMI *cluster1=seed1->GetClusterPointer(irow);
	if (cluster0){
	  if (cluster0->GetQ()>0 &&  cluster0->GetDetector()%36<18)  nclA0++;
	  if (cluster0->GetQ()>0 &&  cluster0->GetDetector()%36>=18) nclC0++;
	}
	if (cluster1){
	  if (cluster1->GetQ()>0 &&  cluster1->GetDetector()%36<18)  nclA1++;
	  if (cluster1->GetQ()>0 &&  cluster1->GetDetector()%36>=18) nclC1++;
	}
      }
      Int_t cosmicType=0;  // types of cosmic topology
      if ((nclA0>nclC0) && (nclA1>nclC1)) cosmicType=0; // AA side
      if ((nclA0<nclC0) && (nclA1<nclC1)) cosmicType=1; // CC side
      if ((nclA0>nclC0) && (nclA1<nclC1)) cosmicType=2; // AC side
      if ((nclA0<nclC0) && (nclA1>nclC1)) cosmicType=3; // CA side
      if (cosmicType<2) continue; // use only crossing tracks
      //
      Double_t deltaTimeCluster=0;
      deltaTimeCluster=0.5*(track1->GetZ()-track0->GetZ())/param->GetZWidth();
      if (nclA0>nclC0) deltaTimeCluster*=-1; // if A side track
      //
      for (Int_t irow=0; irow<159; irow++){
	AliTPCclusterMI *cluster0=seed0->GetClusterPointer(irow);
	if (cluster0 &&cluster0->GetX()>10){
	  Double_t x0[3]={cluster0->GetRow(),cluster0->GetPad(),cluster0->GetTimeBin()+deltaTimeCluster};
	  Int_t index0[1]={cluster0->GetDetector()};
	  transform->Transform(x0,index0,0,1);  
	  cluster0->SetX(x0[0]);
	  cluster0->SetY(x0[1]);
	  cluster0->SetZ(x0[2]);
	  //
	}
	AliTPCclusterMI *cluster1=seed1->GetClusterPointer(irow);
	if (cluster1&&cluster1->GetX()>10){
	  Double_t x1[3]={cluster1->GetRow(),cluster1->GetPad(),cluster1->GetTimeBin()+deltaTimeCluster};
	  Int_t index1[1]={cluster1->GetDetector()};
	  transform->Transform(x1,index1,0,1);  
	  cluster1->SetX(x1[0]);
	  cluster1->SetY(x1[1]);
	  cluster1->SetZ(x1[2]);
	}	
      }
      //
      //
      if (fPIDMatrix){
	(*fPIDMatrix)(itrack0,1)+=4;  //
	(*fPIDMatrix)(itrack1,1)+=4;  //
      }
    }
  }
}



void AliTPCcalibGainMult::ProcessKinks(const AliESDEvent * event){
  //
  //
  //
  AliKFParticle::SetField(event->GetMagneticField()); 
  AliESDfriend *esdFriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!esdFriend) {
    //Printf("ERROR: esdFriend not available");
    return;
  }
  //  if (esdFriend->TestSkipBit()) return;
  //
  // 
  const Double_t kChi2Cut=10;
  const Double_t kMinR=100;
  const Double_t kMaxR=230;
  const Int_t    kMinNcl=110;
  //
  Int_t nkinks = event->GetNumberOfKinks(); 
  AliESDVertex *vertex= (AliESDVertex *)event->GetPrimaryVertex();
  AliKFVertex kfvertex=*vertex;
  TTreeSRedirector * pcstream =  GetDebugStreamer();
  //
  for (Int_t ikink=0;ikink<nkinks;ikink++){
    AliESDkink *kink = event->GetKink(ikink);
    if (!kink) continue;
    if (kink->GetIndex(0)<0) continue;
    if (kink->GetIndex(1)<0) continue;
    if (TMath::Max(kink->GetIndex(1), kink->GetIndex(0))>event->GetNumberOfTracks()) continue;
    //
    //   
    AliExternalTrackParam pd=kink->RefParamDaughter();
    AliExternalTrackParam pm=kink->RefParamMother();
    AliKFParticle kfpd( pd, 211 );
    AliKFParticle kfpm( pm, -13 );
    //
    AliKFParticle *v0KF = new AliKFParticle(kfpm,kfpd); 
    v0KF->SetVtxGuess(kink->GetPosition()[0],kink->GetPosition()[1],kink->GetPosition()[2]);
    Double_t chi2 = v0KF->GetChi2();
    AliESDtrack * trackM = event->GetTrack(kink->GetIndex(0));
    AliESDtrack * trackD = event->GetTrack(kink->GetIndex(1));
    if (!trackM) continue;
    if (!trackD) continue;
    Int_t nclM= (Int_t)trackM->GetTPCClusterInfo(2,1);
    Int_t nclD= (Int_t)trackD->GetTPCClusterInfo(2,1);
    Double_t eta = TMath::Max(TMath::Abs(trackM->Eta()), TMath::Abs(trackD->Eta()));
    Double_t kx= v0KF->GetX();
    Double_t ky= v0KF->GetY();
    Double_t kz= v0KF->GetZ();
    Double_t ex= v0KF->GetErrX();
    Double_t ey= v0KF->GetErrY();
    Double_t ez= v0KF->GetErrZ();
    //
    Double_t radius=TMath::Sqrt(kx*kx+ky*ky);
    Double_t alpha=TMath::ATan2(ky,kx);
    if (!pd.Rotate(alpha)) continue;
    if (!pm.Rotate(alpha)) continue;
    if (!pd.PropagateTo(radius,event->GetMagneticField())) continue;
    if (!pm.PropagateTo(radius,event->GetMagneticField())) continue;
    Double_t pos[2]={0,kz};
    Double_t cov[3]={ex*ex+ey*ey,0,ez*ez};
    pd.Update(pos,cov);
    pm.Update(pos,cov);
    //
    if (pcstream){
      (*pcstream)<<"kinks"<<
	"eta="<<eta<<
	"nclM="<<nclM<<
	"nclD="<<nclD<<
	"kink.="<<kink<<
	"trackM.="<<trackM<<
	"trackD.="<<trackD<<
	"pm.="<<&pm<<             //updated parameters
	"pd.="<<&pd<<             // updated parameters
	"v0KF.="<<v0KF<<
	"chi2="<<chi2<<
	"\n";
    }
    /*
      TCut cutQ="chi2<10&&kink.fRr>90&&kink.fRr<220";
      TCut cutRD="20*sqrt(pd.fC[14])<abs(pd.fP[4])&&trackD.fTPCsignal>10&&trackD.fTPCsignalN>50";
      
    */
    //
    if (chi2>kChi2Cut) continue;
    if (kink->GetR()<kMinR) continue;
    if (kink->GetR()>kMaxR) continue;
    if ((nclM+nclD)<kMinNcl) continue;
    if (TMath::Abs(eta)>1) continue;
    //
    //
    AliESDfriendTrack *friendTrackM = esdFriend->GetTrack(kink->GetIndex(0));
    AliESDfriendTrack *friendTrackD = esdFriend->GetTrack(kink->GetIndex(1));
    if (!friendTrackM) continue;
    if (!friendTrackD) continue;
    TObject *calibObject;
    AliTPCseed *seedM = 0;
    AliTPCseed *seedD = 0;
    for (Int_t l=0;(calibObject=friendTrackM->GetCalibObject(l));++l) {
      if ((seedM=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }    
    for (Int_t l=0;(calibObject=friendTrackD->GetCalibObject(l));++l) {
      if ((seedD=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }    
  }
}

void AliTPCcalibGainMult::DumpHPT(const AliESDEvent * event){
  //
  // Function to select the HPT tracks and events
  // It is used to select event with HPT - list used later for the raw data downloading
  //                                     - and reconstruction
  // Not actualy used for the calibration of the data

  TTreeSRedirector * pcstream =  GetDebugStreamer();
  AliKFParticle::SetField(event->GetMagneticField()); 
  AliESDfriend *esdFriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!esdFriend) {
    //Printf("ERROR: esdFriend not available");
   return;
  }
  if (esdFriend->TestSkipBit()) return;

  Int_t ntracks=event->GetNumberOfTracks(); 
  //
  for (Int_t i=0;i<ntracks;++i) {
    //
    AliESDtrack *track = event->GetTrack(i);
    if (!track) continue;
    if (track->Pt()<4) continue; 
    UInt_t status = track->GetStatus();
    //   
    AliExternalTrackParam * trackIn  = (AliExternalTrackParam *)track->GetInnerParam();
    if (!trackIn) continue;
    if ((status&AliESDtrack::kTPCrefit)==0) continue;
    if ((status&AliESDtrack::kITSrefit)==0) continue;
    AliESDfriendTrack *friendTrack = esdFriend->GetTrack(i);
    if (!friendTrack) continue;
    AliExternalTrackParam * itsOut = (AliExternalTrackParam *)(friendTrack->GetITSOut());
    if (!itsOut) continue;
    AliExternalTrackParam * itsOut2 = (AliExternalTrackParam *)(friendTrack->GetITSOut()->Clone());
    AliExternalTrackParam * tpcIn2 = (AliExternalTrackParam *)(trackIn->Clone());
    if (!itsOut2->Rotate(trackIn->GetAlpha())) continue;
    //Double_t xmiddle=0.5*(itsOut2->GetX()+tpcIn2->GetX());
    Double_t xmiddle=(itsOut2->GetX());
    if (!itsOut2->PropagateTo(xmiddle,event->GetMagneticField())) continue;
    if (!tpcIn2->PropagateTo(xmiddle,event->GetMagneticField())) continue;
    //
    AliExternalTrackParam * tpcInner = (AliExternalTrackParam *)(track->GetTPCInnerParam());
    if (!tpcInner) continue;
    tpcInner->Rotate(track->GetAlpha());
    tpcInner->PropagateTo(track->GetX(),event->GetMagneticField());
    //
    // tpc constrained
    //
    AliExternalTrackParam * tpcInnerC = (AliExternalTrackParam *)(track->GetTPCInnerParam()->Clone());
    if (!tpcInnerC) continue;
    tpcInnerC->Rotate(track->GetAlpha());
    tpcInnerC->PropagateTo(track->GetX(),event->GetMagneticField());
    Double_t dz[2],cov[3];
    AliESDVertex *vtx= (AliESDVertex *)event->GetPrimaryVertex();
  
    if (!tpcInnerC->PropagateToDCA(vtx, event->GetMagneticField(), 3, dz, cov)) continue;
    Double_t covar[6]; vtx->GetCovMatrix(covar);
    Double_t p[2]={tpcInnerC->GetParameter()[0]-dz[0],tpcInnerC->GetParameter()[1]-dz[1]};
    Double_t c[3]={covar[2],0.,covar[5]};
    //
    Double_t chi2C=tpcInnerC->GetPredictedChi2(p,c);
    tpcInnerC->Update(p,c);

    if (pcstream){
      (*pcstream)<<"hpt"<<
	"run="<<fRun<<
	"time="<<fTime<<
	"vertex="<<vtx<<
	"bz="<<fMagF<<
	"track.="<<track<<
	"tpcInner.="<<tpcInner<<
	"tpcInnerC.="<<tpcInnerC<<
	"chi2C="<<chi2C<<
	//
	"its.="<<itsOut<<
	"its2.="<<itsOut2<<
	"tpc.="<<trackIn<<
	"tpc2.="<<tpcIn2<<
	"\n";
    }
  }
}



void AliTPCcalibGainMult::ProcessTOF(const AliESDEvent * event){
  //
  // 1. Loop over tracks
  // 2. Fit T0
  // 3. Sign positivelly identified tracks
  // 
  const Double_t kMaxDelta=1000;
  const Double_t kOrbit=50000; // distance in the time beween two orbits in the TOF time units  - 50000=50 ns
  const Double_t kMaxD=20;
  const Double_t kRMS0=200; 
  const Double_t kMaxDCAZ=10;
  AliESDVertex *vtx= (AliESDVertex *)event->GetPrimaryVertex();
  //
  TTreeSRedirector * pcstream =  GetDebugStreamer();
  AliKFParticle::SetField(event->GetMagneticField()); 
  AliESDfriend *esdFriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  if (!esdFriend) {
    //Printf("ERROR: esdFriend not available");
   return;
  }
  //if (esdFriend->TestSkipBit()) return;

  Int_t ntracks=event->GetNumberOfTracks(); 
  //
  Double_t deltaTPion[10000];
  Double_t medianT0=0;
  Double_t meanT0=0;
  Double_t rms=10000;
  Int_t counter=0;
  //
  // Get Median time for pion hypothesy
  //
  for (Int_t iter=0; iter<3; iter++){
    counter=0;
    for (Int_t i=0;i<ntracks;++i) {
      //
      AliESDtrack *track = event->GetTrack(i);
      if (!track) continue;
      if (!track->IsOn(AliESDtrack::kTIME)) continue;
      if (TMath::Abs(track->GetZ())>kMaxDCAZ) continue;         // remove overlaped events
      if (TMath::Abs(track->GetTOFsignalDz())>kMaxD) continue;
      Double_t times[1000];
      track->GetIntegratedTimes(times);
      Int_t norbit=TMath::Nint((track->GetTOFsignal()-times[2])/kOrbit);
      Double_t torbit=norbit*kOrbit; 
      if (iter==1 &&TMath::Abs(times[2]-times[3])<3*rms) continue;  // skip umbigous points - kaon pion
      //
      Int_t indexBest=2;
      if (iter>1){
	for (Int_t j=3; j<5; j++) 
	  if (TMath::Abs(track->GetTOFsignal()-times[j]-torbit-medianT0)<TMath::Abs(track->GetTOFsignal()-times[j]-torbit-medianT0)) indexBest=j; 
      }
      //
      if (iter>0) if (TMath::Abs(track->GetTOFsignal()-times[indexBest]-torbit-medianT0)>3*(kRMS0+rms)) continue;
      if (iter>0) if (TMath::Abs(track->GetTOFsignal()-times[indexBest]-torbit-medianT0)>kMaxDelta) continue;
      deltaTPion[counter]=track->GetTOFsignal()-times[indexBest]-torbit;
      counter++;
    }    
    if (counter<2) return;
    medianT0=TMath::Median(counter,deltaTPion);    
    meanT0=TMath::Median(counter,deltaTPion);    
    rms=TMath::RMS(counter,deltaTPion);    
  }
  if (counter<3) return;
  //
  // Dump
  //
  for (Int_t i=0;i<ntracks;++i) {
    //
    AliESDtrack *track = event->GetTrack(i);
    if (!track) continue;
    if (!track->IsOn(AliESDtrack::kTIME)) continue;
    if (TMath::Abs(track->GetZ())>kMaxDCAZ) continue;          //remove overlapped events
    if (TMath::Abs(track->GetTOFsignalDz())>kMaxD) continue;
    Double_t times[1000];
    track->GetIntegratedTimes(times);  
    Int_t norbit=TMath::Nint((track->GetTOFsignal()-times[2])/kOrbit);
    Double_t torbit=norbit*kOrbit;
    if (rms<=0) continue;
    //
    Double_t tPion  = (track->GetTOFsignal()-times[2]-medianT0-torbit);
    Double_t tKaon  = (track->GetTOFsignal()-times[3]-medianT0-torbit);
    Double_t tProton= (track->GetTOFsignal()-times[4]-medianT0-torbit);
    Double_t tElectron= (track->GetTOFsignal()-times[0]-medianT0-torbit);
    //
    Bool_t isPion   = (TMath::Abs(tPion/rms)<6)   && TMath::Abs(tPion)<(TMath::Min(TMath::Abs(tKaon), TMath::Abs(tProton))-rms);
    Bool_t isKaon   = (TMath::Abs(tKaon/rms)<3)   && TMath::Abs(tKaon)<0.2*(TMath::Min(TMath::Abs(tPion), TMath::Abs(tProton))-3*rms);
    Bool_t isProton = (TMath::Abs(tProton/rms)<6) && TMath::Abs(tProton)<0.5*(TMath::Min(TMath::Abs(tKaon), TMath::Abs(tPion))-rms);
    Bool_t isElectron = (TMath::Abs(tElectron/rms)<6) && TMath::Abs(tElectron)<0.2*(TMath::Min(TMath::Abs(tKaon), TMath::Abs(tPion))-rms) &&TMath::Abs(tElectron)<0.5*(TMath::Abs(tPion)-rms);

    if (isPion)   (*fPIDMatrix)(i,2)+=1;
    if (isKaon)   (*fPIDMatrix)(i,3)+=1;
    if (isProton) (*fPIDMatrix)(i,4)+=1;
    //    if (isElectron) (*fPIDMatrix)(i,0)+=1;
    //
    if (pcstream){
      // debug streamer to dump the information 
    (*pcstream)<<"tof"<<
      "isPion="<<isPion<<
      "isKaon="<<isKaon<<
      "isProton="<<isProton<<
      "isElectron="<<isElectron<<
      //
      "counter="<<counter<<
      "torbit="<<torbit<<
      "norbit="<<norbit<<
      "medianT0="<<medianT0<<
      "meanT0="<<meanT0<<
      "rmsT0="<<rms<<
      "track.="<<track<<
      "vtx.="<<vtx<<
      "\n";
    }

  }
  /*
    tof->SetAlias("isProton","(abs(track.fTOFsignal-track.fTrackTime[4]-medianT0-torbit)<(0.5*abs(track.fTOFsignal-track.fTrackTime[3]-medianT0-torbit)-rmsT0))");
    tof->SetAlias("isPion","(abs(track.fTOFsignal-track.fTrackTime[2]-medianT0-torbit)<(0.5*abs(track.fTOFsignal-track.fTrackTime[3]-medianT0-torbit)-rmsT0))");
    tof->SetAlias("isKaon","(abs(track.fTOFsignal-track.fTrackTime[3]-medianT0-torbit)<(0.5*abs(track.fTOFsignal-track.fTrackTime[2]-medianT0-torbit)-rmsT0))&&(abs(track.fTOFsignal-track.fTrackTime[3]-medianT0-torbit)<(0.5*abs(track.fTOFsignal-track.fTrackTime[4]-medianT0-torbit)-rmsT0))");
    
   */

}


// void AliTPCcalibGainMult::Terminate(){
//   //
//   // Terminate function
//   // call base terminate + Eval of fitters
//   //
//    Info("AliTPCcalibGainMult","Terminate");
//    TTreeSRedirector *pcstream = GetDebugStreamer();
//    if (pcstream){
//      TTreeStream &stream = (*pcstream)<<"dump";
//      TTree* tree = stream.GetTree();
//      if (tree) if ( tree->GetEntries()>0){
//        TObjArray *array =  tree->GetListOfBranches(); 
//        for (Int_t i=0; i<array->GetEntries(); i++) {TBranch * br = (TBranch *)array->At(i); br->SetAddress(0);}      
//        gDirectory=gROOT;
//        fdEdxTree=tree->CloneTree(10000);
//      }
//    }
//    AliTPCcalibBase::Terminate();
// }

