/*

  //
  // Macro to compare  tracks used 2 different reconstruction algorithms (Offline and HLT)
  // (In general the approach  can be used for any pair of the esds - track containers)
  //   
  // Main functions - Dump  content of HLT and OFFLINE esd files  - full esd tracks 
  //  1. Dump to the tree the pair of close tracks ( n sigma distance at the TPC entrance)
  //  2. Dump to the tree the track and corrsponding counter of tracks from other container 
  //    Compare two events - track by track comparison
  //       0. Filter "Good" tpc tracks - see function IsSelected(track)
  //       1. Dump OFFline -> HLT tracks into trees  
  //           a.) tree "offhlt"  - counters+ Ofline track + corresponding  Hlt track
  //           b.) tree "offhlt0" - offline track + hlt track counter 
  //       2. HLT and offline
  //     
  //  3. To save the CPU and disk space the input tracks to compare are pt downscaled
  //      see function IsDownscaled();
  //
  //  Visualization part is currently work in progress
  //  
  // Responsible:
  // dumping part:         marian.ivanov@cern.ch
  // visualization part:   
  // 
  Example usage:
  aliroot -b -q $ALICE_PHYSICS/TPC/macros/compareTracks.C+
  //
  //
  .x $HOME/rootlogon.C
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");
  gSystem->Load("libPWG0selectors");
  .L $ALICE_PHYSICS/PWGPP/TPC/macros/compareTracks.C+
  //
  compareTracks("compare.list");
  //
  TFile f("dump.root");
*/
#include <exception>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "AliESDEvent.h"
#include "TTreeStream.h"
#include "TRandom.h"
#include "AliXRDPROOFtoolkit.h"
#include "AliSysInfo.h"

//
//
//
void CompareFile(TTreeSRedirector *pcstream,  TFile *fOff, TFile *fHLT);
void CompareEvents(AliESDEvent *evoff, AliESDEvent *evhlt, TTreeSRedirector *pcstream);
Bool_t AreTracksCloseFast(AliESDtrack *track1, AliESDtrack *track2, AliESDEvent * event, AliExternalTrackParam &inner2);
Bool_t IsSelected(AliESDtrack *track);
Bool_t IsDownscaled(AliESDtrack *track);
//
Double_t gptdownscale=100;  //pt downscale parameter
//
TTree * chainOFFHLT=0;
TTree * chainOFFHLT0=0;
TTree * chainHLTOFF=0;
TTree * chainHLTOFF0=0;

void compareTracks(const char * flistOFF="compareOFF.list", const char * flistHLT="compareHLT.list", Double_t downscale=10000){
    //
    // Compare HLT and OFFLINE esd files 
    // Input: 
    //   flist - the ascii files with the filename - ESD form the offline
    //         - hlt path constructed  from the offline ESD file repalcing the OFFLINE with HLT
    //
    gptdownscale=downscale;
    TTreeSRedirector *pcstream = new TTreeSRedirector("dump.root");
    ifstream inOFF;
    inOFF.open(flistOFF);
    ifstream inHLT;
    inHLT.open(flistHLT);


    // Read the input list of files and add them to the chain
    TString currentFileOff;
    TString currentFileHlt;
    {  while( (inOFF.good()) && (inHLT.good())) {
		        inOFF >> currentFileOff;
			inHLT >> currentFileHlt;
	                //currentFileHlt=currentFileOff;
	                //currentFileHlt.ReplaceAll("/OFFLINE/","/HLT/");
	                printf("%s\n%s\n\n",currentFileOff.Data(), currentFileHlt.Data());
	                TFile *fOff = TFile::Open(currentFileOff.Data(),"READ");
		        TFile *fHlt = TFile::Open(currentFileHlt.Data(),"READ");
	                if ( fOff && fHlt){
		           CompareFile(pcstream, fOff, fHlt);
		        }
		      delete fOff;
		      delete fHlt;
    }}
    delete pcstream;
		        /* 
			 *      TTreeSRedirector *pcstream = new TTreeSRedirector("dump.root");
			 *           TFile *fOff = TFile::Open("/lustre/alice/mknichel/cpass1/2011-12-13_0154/OFFLINE/167693/11000167693001.11/AliESDs.root","READ");
			 *                TFile *fHLT = TFile::Open("/lustre/alice/mknichel/cpass1/2011-12-13_0154/HLT/167693/11000167693001.11/AliESDs.root","READ");
			 *                     CompareFile(pcstream, fOff, fHLT);
			 *                          delete pcstream;
			 *                            */
}

void CompareFile(TTreeSRedirector *pcstream,  TFile *fOff, TFile *fHLT){
  //
  // Compare 2 esd files - fOff and fHLT
  // Pairs of tracks dumpt in the trees - stored in the filee connected to the pcstream
  // 
  TTree *tOff = (TTree*)fOff->Get("esdTree");
  TTree *tHLT = (TTree*)fHLT->Get("esdTree");
  //    
  Int_t nevents=tOff->GetEntries();
  //
  AliESDEvent *evoff = new AliESDEvent();
  AliESDEvent *evhlt = new AliESDEvent();  
  evoff->ReadFromTree(tOff);
  evhlt->ReadFromTree(tHLT);
  //
  for(Long64_t iev = 0; iev < nevents; iev++){
    tOff->GetEvent(iev);
    tHLT->GetEvent(iev);
    CompareEvents(evoff,evhlt,pcstream);    
    AliSysInfo::AddStamp(fOff->GetName(),iev);
  }
  delete evoff;
  delete evhlt;
}



void CompareEvents(AliESDEvent *evoff, AliESDEvent *evhlt, TTreeSRedirector *pcstream){
  //
  // Main function to dump events:
  //
  // compare two events - track by track comparison
  // 0. Filter "Good" tpc tracks - see function IsSelected(track)
  // 1. Dump OFFline -> HLT tracks into trees  
  //    a.) tree "offhlt" - counters+ Ofline track + corresponding  Hlt track
  //    b.) tree "offhlt" - counters +Offline track
  // 2. HLT and offline
  //
  // 3. To save the CPU and disk space the input tracks to compare are pt downscaled
  //    see function Isdownscaled();

  Int_t nOffTracks = evoff->GetNumberOfTracks();
  Int_t nHLTTracks = evhlt->GetNumberOfTracks();
  Int_t iev= evoff->GetEventNumberInFile();
  if (nOffTracks==0) return;
  AliESDtrack *offlinetrack=0, *hlttrack=0;
  Double_t predchi2 = 0;  
  //for(Long64_t iev = 120; iev < 125; iev++)
  TObjArray arrayOFF(20000);
  TObjArray arrayHLT(20000);
  AliExternalTrackParam  inner2; //working track ref
  //
  //
  //
  Int_t nselOFF = 0;
  Int_t nselHLT = 0;
  // filter offline
  for(Int_t ioff = 0; ioff < nOffTracks; ioff++){
    offlinetrack = evoff->GetTrack(ioff);
    if (!IsSelected(offlinetrack)) continue;
    arrayOFF.AddAt(offlinetrack,nselOFF);	
    nselOFF++;
  }
  //filter hlt
  for(Int_t ihlt = 0; ihlt< nHLTTracks; ihlt++){
    hlttrack = evhlt->GetTrack(ihlt);
    if (!IsSelected(hlttrack)) continue;
    arrayHLT.AddAt(hlttrack,nselHLT);	
    nselHLT++;
  }
  printf("%d\t%d\t%d\t%d\n", nOffTracks, nselOFF,  nHLTTracks, nselHLT);
  //
  // - compare offline - hlt
  TStopwatch timer;
  //
  for (Int_t ioff = 0; ioff < nselOFF; ioff++){
    offlinetrack = (AliESDtrack*)arrayOFF.At(ioff);
    Float_t ncl21off= offlinetrack->GetTPCClusterInfo(2,1);
    Float_t ncl20off= offlinetrack->GetTPCClusterInfo(2,0);
    if (IsDownscaled(offlinetrack)) continue;
    Int_t counter=0;
    //
    for(Int_t ihlt = 0; ihlt < nselHLT; ihlt++){
      hlttrack = (AliESDtrack*)arrayHLT.At(ihlt);
      Bool_t close = AreTracksCloseFast(offlinetrack,hlttrack,evhlt,inner2);	   
      if (close){
	Float_t ncl21hlt= hlttrack->GetTPCClusterInfo(2,1);
	Float_t ncl20hlt= hlttrack->GetTPCClusterInfo(2,0);
	counter++;
	predchi2=offlinetrack->GetInnerParam()->GetPredictedChi2(&inner2);
	(*pcstream)<<"offhlt"<<
	  //multiplicity
	  "nHLT="<<nselHLT<<
	  "nOFF="<<nselOFF<<
	  "ncl20off="<<ncl20off<<
	  "ncl21off="<<ncl21off<<	
	  "ncl20hlt="<<ncl20hlt<<
	  "ncl21hlt="<<ncl21hlt<<	       
	  "counter="<<counter<<
	  //
	  "chi2="<<predchi2<<
	  "track1.="<<offlinetrack<<
	  "track2.="<<hlttrack<<
	  "inner2.="<<&inner2<<
	  "\n";
      }
    }
    Bool_t isDownscaled=IsDownscaled(offlinetrack);	
    (*pcstream)<<"offhlt0"<<
      "iev="<<iev<<	       
      "nHLT="<<nselHLT<<
      "nOFF="<<nselOFF<<
      //
      "ncl20="<<ncl20off<<
      "ncl21="<<ncl21off<<
      "track.="<<offlinetrack<<
      "counter="<<counter<<
      "isDownscaled="<<isDownscaled<<
      "\n";
  }
  timer.Print();
  //
  for(Int_t ihlt = 0; ihlt < nselHLT; ihlt++){
    Int_t counter=0;      
    hlttrack = (AliESDtrack*)arrayHLT.At(ihlt);	
    Float_t ncl21hlt= hlttrack->GetTPCClusterInfo(2,1);
    Float_t ncl20hlt= hlttrack->GetTPCClusterInfo(2,0);
    if (IsDownscaled(hlttrack)) continue;      
    for (Int_t ioff = 0; ioff < nselOFF; ioff++){
      offlinetrack = (AliESDtrack*)arrayOFF.At(ioff);	
      Bool_t close = AreTracksCloseFast(hlttrack,offlinetrack,evhlt,inner2);	   
      if (close){
	Float_t ncl21off= offlinetrack->GetTPCClusterInfo(2,1);
	Float_t ncl20off= offlinetrack->GetTPCClusterInfo(2,0);
	predchi2=hlttrack->GetInnerParam()->GetPredictedChi2(&inner2);
	counter++;
	(*pcstream)<<"hltoff"<<
	  //multiplicity
	  "nHLT="<<nselHLT<<
	  "nOFF="<<nselOFF<<
	  "ncl20off="<<ncl20off<<
	  "ncl21off="<<ncl21off<<
	  "ncl20hlt="<<ncl20hlt<<
	  "ncl21hlt="<<ncl21hlt<<	       
	  "counter="<<counter<<
	  //
	  "chi2="<<predchi2<<
	  "track2.="<<offlinetrack<<
	  "track1.="<<hlttrack<<
	  "inner2.="<<&inner2<<
	  "\n";
      }
    }
    Bool_t isDownscaled=IsDownscaled(hlttrack);	
    (*pcstream)<<"hltoff0"<<
      "iev="<<iev<<	       
      "nHLT="<<nselHLT<<
      "nOFF="<<nselOFF<<
      //
      "ncl20="<<ncl20hlt<<
      "ncl21="<<ncl21hlt<<
      "track.="<<hlttrack<<
      "counter="<<counter<<
      "isDownscaled="<<isDownscaled<<
      "\n";
  } 
  timer.Print();
}

Bool_t IsSelected(AliESDtrack *track){
  //
  // modified by:
  // philipp.luettig@cern.ch
  // 
  //
  if (track->IsOn(0x40)==0) return kFALSE;                // Refit
  if (TMath::Abs(track->GetTgl())>1.1)  return kFALSE;    // tangent lambda
  if (track->GetTPCClusterInfo(2,1)<50) return kFALSE;    // cluster information number of crossed rows >50
  return kTRUE;
}

Bool_t IsDownscaled(AliESDtrack *track){
  //
  // Downscale randomly low pt tracks
  //
  //return kFALSE;
  Double_t scalempt= TMath::Min(1./TMath::Abs(track->GetParameter()[4]),10.);
  if (TMath::Exp(2*scalempt)<gptdownscale*gRandom->Rndm()) return kTRUE;
  return kFALSE;
}


Bool_t AreTracksCloseFast(AliESDtrack *track1, AliESDtrack *track2, AliESDEvent * event, AliExternalTrackParam &inner2){
  //
  // 
  // Fast comparison uning track param close to the prim vertex and at the inner wall of the TPC
  //
  // 1. Fast cut on the invariant (under rotation) variable P3 and P4 (means pz/pt and 1/pt)
  // 2. Slower cuts - parameters at the entrance of the TPC (tracks to be propagated)
  //
  // In case the tracks are relativelaly close -the inner2 parameters are created
  //                                           -track 2 propagated and rotated to the same position as the track1 
  // 
  const Double_t absCut[5] ={5,10, 0.02,0.02,1};   // abs cut values  
  const Double_t pullCut[5]={6,100,6,   100,6};    // pull cut values
  //   *              "External" track parametrisation class                       *
  //   *                                                                           *
  //   *      external param0:   local Y-coordinate of a track (cm)                *
  //   *      external param1:   local Z-coordinate of a track (cm)                *
  //   *      external param2:   local sine of the track momentum azimuthal angle  *
  //   *      external param3:   tangent of the track momentum dip angle           *
  //   *      external param4:   1/pt (1/(GeV/c))                                  *

  //
  // 
  const Double_t kTglCut=0.1;    
  const Double_t k1PtCut=0.5;
  //  const Double_t kAlphaCut=0.2;
  const Double_t kTglCutSigma=10;
  const Double_t k1PtCutSigma=10;
  //
  //
  if(!track1) return kFALSE;
  if(!track2) return kFALSE;  
  const Double_t *param1 = track1->GetParameter();
  const Double_t *param2 = track2->GetParameter();
  const Double_t *param1I = track1->GetInnerParam()->GetParameter();
  const Double_t *param2I = track2->GetInnerParam()->GetParameter();
  const Double_t *covar1 = track1->GetCovariance();
  const Double_t *covar2 = track2->GetCovariance();
  const Double_t *covar1I = track1->GetInnerParam()->GetCovariance();
  const Double_t *covar2I = track2->GetInnerParam()->GetCovariance();
  //
  if (TMath::Abs(param1[3]-param2[3])>kTglCut) return kFALSE;	
  //if (TMath::Abs(param1[4]-param2[4])>k1PtCut) return kFALSE;	
  if (TMath::Abs(param1I[3]-param2I[3])>kTglCut) return kFALSE;	
  if (TMath::Abs(param1I[4]-param2I[4])>k1PtCut) return kFALSE;	
  Double_t dalpha = TMath::Abs(track1->GetAlpha()-track2->GetAlpha());
  if (dalpha>TMath::Pi()) dalpha-=TMath::Abs(dalpha-TMath::TwoPi());
  //if (dalpha > kAlphaCut) return kFALSE;

  //
  Int_t index22=track1->GetIndex(2,2);	       
  Int_t index33=track1->GetIndex(3,3);	       
  Int_t index44=track1->GetIndex(4,4);	       
  if (TMath::Abs(param1[3]-param2[3])/TMath::Sqrt(TMath::Max(covar1[index33],covar2[index33]))>kTglCutSigma) return kFALSE;
  //if (TMath::Abs(param1[4]-param2[4])/TMath::Sqrt(TMath::Max(covar1[index44],covar2[index44]))>k1PtCutSigma) return kFALSE;  
  if (TMath::Abs(param1I[3]-param2I[3])/TMath::Sqrt(TMath::Max(covar1I[index33],covar2I[index33]))>kTglCutSigma) return kFALSE;
  if (TMath::Abs(param1I[4]-param2I[4])/TMath::Sqrt(TMath::Max(covar1I[index44],covar2I[index44]))>k1PtCutSigma) return kFALSE;  
  if (TMath::Abs(dalpha)/TMath::Sqrt(TMath::Max(covar1[index22],covar2[index22]))>k1PtCutSigma) return kFALSE;  
  //
  // 2. Slow cuts on the paramters at the entrance of the TPC
  //
  inner2=*(track2->GetInnerParam());
  inner2.Rotate(track1->GetInnerParam()->GetAlpha());
  inner2.PropagateTo(track1->GetInnerParam()->GetX(),  event->GetMagneticField());
  const Double_t *pinner2  = inner2.GetParameter();
  const Double_t *pcovar2  = inner2.GetCovariance();
  //
  Bool_t isOK = kTRUE;
  for (Int_t ipar=0; ipar<5; ipar++){
    // 
    Int_t index=track1->GetIndex(ipar,ipar);	      
    if (TMath::Abs(pinner2[ipar]-param1I[ipar])>absCut[ipar]) isOK=kFALSE;
    if (TMath::Abs(pinner2[ipar]-param1I[ipar])>pullCut[ipar]*TMath::Sqrt(TMath::Max(covar1I[index],pcovar2[index]))) isOK=kFALSE;
  }
  if (!isOK) return kFALSE; 
  return kTRUE;
}


void MakeChain(){
  //
  //
  AliXRDPROOFtoolkit toolkit;
  chainOFFHLT= toolkit.MakeChainRandom("dumpHLTOFFLINE.list","offhlt",0,100);
  chainOFFHLT0=toolkit.MakeChainRandom("dumpHLTOFFLINE.list","offhlt0",0,100);
  chainHLTOFF=toolkit.MakeChainRandom("dumpHLTOFFLINE.list","hltoff",0,100);
  chainHLTOFF0=toolkit.MakeChainRandom("dumpHLTOFFLINE.list","hltoff0",0,100);

 
}

void DrawDiffPt(){
  //
  // Draw difference between the HLT and offline tracks
  //
  TCut cut="sqrt(chi2)<10&&ncl21off>120";
  TCut cutNoiseEvent = "abs(nHLT/nOFF-1)<0.2";   //mainly laser events
  //
  // 1. check the edge effect 1/pt resolution TPC only pull
  // ...
  chainOFFHLT->Draw("(track1.fIp.fP[4]-track2.fIp.fP[4])/sqrt(max(track1.fIp.fC[14],track2.fIp.fC[14])):sign(inner2.fP[4])*inner2.fP[0]/inner2.fX>>hisTPCEdge(50,-0.18,0.18,100,-6,6)",cut+"abs(track1.fP[4])<0.25","colz",200000);
  /*
    hisTPCEdge->FitSlicesY();
    hisTPCEdge_2->GetXaxis()->SetTitle("q*ly/lx");
    hisTPCEdge_2->GetYaxis()->SetTitle("#Delta_{1/pt}/#sigma_{1/pt}");
    hisTPCEdge_2->Draw();
  */
  // 2. check the edge effect 1/pt resolution combined 
  chainOFFHLT->Draw("(track1.fP[4]-track2.fP[4])/sqrt(max(track1.fC[14],track2.fC[14])):sign(inner2.fP[4])*inner2.fP[0]/inner2.fX>>hisCombEdge(50,-0.18,0.18,100,-6,6)",cut+"abs(track1.fP[4])<0.25","colz",200000);
  /*
    hisCombEdge->FitSlicesY();
    hisCombEdge_2->Draw();
  */
  // 3. Combined momentum resolution as function of the inverse moment
  chainOFFHLT->Draw("(track1.fIp.fP[4]-track2.fIp.fP[4])/sqrt(max(track1.fIp.fC[14],track2.fIp.fC[14])):abs(track1.fP[4])>>hisTPCP4(20,-0.0,1,100,-6,6)",cut+"abs(track1.fP[4])<1","colz",200000);
  /*
    hisTPCP4->FitSlicesY();  
    hisTPCP4_2->GetXaxis()->SetTitle("1/p_{t} (1/GeV))");
    hisTPCP4_2->GetYaxis()->SetTitle("#Delta_{1/pt}/#sigma_{1/pt}");
    hisTPCP4_2->Draw();
  */
  // 4. Combined momentum resolution as function of the inverse moment
  chainOFFHLT->Draw("(track1.fP[4]-track2.fP[4])/sqrt(max(track1.fC[14],track2.fC[14])):abs(track1.fP[4])>>hisCombP4(20,-0.0,1,100,-6,6)",cut+"abs(track1.fP[4])<1","colz",200000);
  /*
    hisCombP4->FitSlicesY();
    hisCombP4_2->Draw();
  */

}
//
void DrawDiffEff(){
  //
  //
  //
  TCut cutEff = "abs(track.fIp.fP[4])<5&&abs(track.fIp.fP[1])<90";
  TCut cutNoiseEvent = "abs(nHLT/nOFF-1)<0.2"; // HLT cluster finder more sensitive to the noise
  //
  // 
  //
  chainOFFHLT0->Draw("counter==0:(nOFF+nHLT)/2.>>effOccuOFFHLT(20,0,8000)",cutNoiseEvent+cutEff+"abs(track.fP[4])<1&&ncl21>120","prof",50000);
  chainHLTOFF0->Draw("counter==0:(nOFF+nHLT)/2.>>effOccuHLTOFF(20,0,8000)",cutNoiseEvent+cutEff+"abs(track.fP[4])<1&&ncl21>120","prof",50000);

 chainOFFHLT0->Draw("counter==0:track.fTPCncls>>effNCLOFFHLT(40,0,160)",cutNoiseEvent+cutEff+"abs(track.fP[4])<1","prof",50000);
 chainHLTOFF0->Draw("counter==0:track.fTPCncls>>effNCLHLTOFF(40,0,160)",cutNoiseEvent+cutEff+"abs(track.fP[4])<1","prof",50000);
 /*
   effNCLOFFHLT->SetMarkerStyle(25);
   effNCLHLTOFF->SetMarkerStyle(25);
   effNCLOFFHLT->SetMarkerColor(2);
   effNCLHLTOFF->SetMarkerColor(4);
   effNCLOFFHLT->Draw();
   effNCLHLTOFF->Draw("same");
  */


  chainHLTOFF0->Draw("counter==0:sign(track.fIp.fP[4])*track.fIp.fP[0]/track.fIp.fX>>profTPCEdge(50,-0.18,0.18)",cutNoiseEvent+cutEff+"abs(track.fP[4])<0.25","prof",50000);

  chainOFFHLT0->Draw("counter==0:sign(track.fIp.fP[4])*track.fIp.fP[0]/track.fIp.fX>>profTPCEdge(50,-0.18,0.18)",cutNoiseEvent+cutEff+"abs(track.fP[4])<1","prof",50000);

  chainOFFHLT0->Draw("track.fTPCncls:sign(track.fIp.fP[4])*track.fIp.fP[0]/track.fIp.fX>>profTPCEdge(50,-0.18,0.18)",cutNoiseEvent+cutEff+"abs(track.fP[4])<1","prof",50000);


}

/*
  This is  shell script real example  to submit jobs for the track comparison: 
  //
  //
  rm -rf list*
  rm -rf dirlist*
  wdir=`pwd`
  split offline.list --lines=50 list -d
  for a in `ls list*`; do
    mkdir dir$a
    cd dir$a
    mv ../$a compare.list
    bsub -q proof -oo outcompare.log aliroot -b -q $ALICE_PHYSICS/PWGPP/TPC/macroscompareTracks.C+
    cd $wdir
  done;

  find  `pwd`/  | grep .root > dumpHLTOFFLINE.list

*/
