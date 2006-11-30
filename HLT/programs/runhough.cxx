// $Id$

// Author: Constantin Loizides <loizides@ikf.uni-frankfurt.de
//*-- Copyright &copy ALICE HLT Group

#include <AliHLTStandardIncludes.h>

#include <AliHLTRootTypes.h>
#include <AliHLTLogging.h>
#include <AliHLTLogger.h>
#include <AliHLTTransform.h>
#include <AliHLTTrack.h>
#include <AliHLTTrackArray.h>
#include <AliHLTHoughTrack.h>
#include <AliHLTClustFinderNew.h>
#include <AliHLTMemHandler.h>
#include <AliHLTSpacePointData.h>
#include <AliHLTHoughBaseTransformer.h>
#include <AliHLTHoughTransformer.h>
#include <AliHLTHoughTransformerLUT.h>
#include <AliHLTHoughTransformerVhdl.h>
#include <AliHLTHoughMaxFinder.h>
#include <AliHLTHough.h>

#ifndef no_root
#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#endif

#if __GNUC__ == 3
using namespace std;
#endif

int main(Int_t argc,Char_t **argv)
{
  Int_t sl=0;
  Int_t sh=0;
  Int_t segs=100;

  AliHLTLogger l;
  l.Set(AliHLTLogger::kAll);
  l.UseStderr();
  //l.UseStdout();
  //l.UseStream();

  Char_t path[1024];
  if(argc<2){
    cout<<"Usage: runhough path [slow] [shigh] [segs]"<<endl;
    exit(1);
  }
  strcpy(path,argv[1]);
  if (argc>2) {
    sl=atoi(argv[2]);
  }
  if (argc>3) {
    sh=atoi(argv[3]);
  }
  if (argc>4) {
    segs=atoi(argv[4]);
  }

  //AliHLTFFloat::SetParams(10000);
  AliHLTTransform::Init(path);

#if 0
  runhough(sl,sh,path,segs);
#else //do some comparison tests

  AliHLTHoughBaseTransformer *fh1 = new AliHLTHoughTransformerVhdl(0,0,segs);
  AliHLTHoughBaseTransformer *fh2 = new AliHLTHoughTransformerLUT(0,0,segs);

  fh1->CreateHistograms(64,0.1,64,-30.,30.);
  fh2->CreateHistograms(64,0.1,64,-30.,30.);

  fh1->Print();

#endif

  //AliHLTFFloat::PrintStat();
  exit(0);
}


//----------------------------------------------------------------------------------
// dont look beyond...
//----------------------------------------------------------------------------------

#if 0
void runhough(Int_t sl,Int_t sh,Char_t *path,Int_t n_eta_segments, Int_t show_seg=-1)
{

  Bool_t binary = kTRUE;
  Bool_t bit8 = kTRUE;
  Int_t tv=1;
  Int_t th=14000;

  AliHLTHough *hough = new AliHLTHough();
  hough->Init(path,binary,n_eta_segments,bit8,tv);
  hough->GetMaxFinder()->SetThreshold(th);
  Int_t ntracks=0;

  for(Int_t slice=sl;slice<=sh;slice++){
    hough->ReadData(slice);
    hough->Transform();
    hough->AddAllHistograms();
    hough->FindTrackCandidates();
    //hough->Evaluate(5);

    AliHLTTrackArray *tracks = (AliHLTTrackArray*)hough->GetTracks(0);
    ntracks=tracks->GetNTracks();
    for(int i=0; i<ntracks; i++)
      {
	AliHLTHoughTrack *track = (AliHLTHoughTrack*)tracks->GetCheckedTrack(i);
	if(!track) continue;
	if(sl==sh) cout<<"pt "<<track->GetPt()<<" psi "<<track->GetPsi()<<" eta "<<track->GetEta()<<" etaindex "<<track->GetEtaIndex()<<" weight "<<track->GetWeight()<<endl;
	if(show_seg<0) show_seg=track->GetEtaIndex();
      }

    cout<<"Found total "<<tracks->GetNTracks()<<" tracks"<<endl;
    hough->WriteTracks(slice);
  }

  //if((ntracks>0)&&(sl==sh)) display(hough,show_seg);
}
#endif

#if 0
void display(AliHLTHough *hough,Int_t eta_index)
{
  //Display the data/tracks in eta_index
  
  hough->InitEvaluate();
  AliHLTHistogram *digitd = new AliHLTHistogram("Digits display","",250,0,250,250,-125,125);
  AliHLTHistogram *trackd = new AliHLTHistogram("Found tracks display","",250,0,250,250,-125,125);
  for(int i=0; i<6; i++)
    hough->GetEval(i)->DisplayEtaSlice(eta_index,digitd);
  
  float xyz[3];
  tracks = (AliHLTTrackArray*)hough->GetTracks(0);
  for(int i=0; i<tracks->GetNTracks(); i++)
    {
      AliHLTHoughTrack *track = (AliHLTHoughTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      if(track->GetEtaIndex() != eta_index) continue;

      for(int j=0; j<176; j++)
	{
	  track->GetCrossingPoint(j,xyz);
	  trackd->Fill(xyz[0],xyz[1],1);
	}
    }
  
  //Draw the parameter space
  TCanvas *c1 = new TCanvas("c1","",2);
  hough->GetTransformer(0)->GetHistogram(eta_index)->Draw("box");
  
  //Draw the tracks
  TCanvas *c2 = new TCanvas("c2","",2);
  digitd->Draw();
  trackd->Draw("same");
  ((TH1F*)trackd->GetRootHisto())->SetMarkerColor(2);
}
#endif

#if 0
struct GoodTrack
{
  Int_t event;
  Int_t label;
  Double_t eta;
  Int_t code;
  Double_t px,py,pz;
  Double_t pt;
  Int_t nhits;
};

void geteff(char *fname)
{
  GoodTrack gt[15000];
  int counter=0;
  ifstream in(fname);
  if(!in)
    {
      cerr<<"Could not open "<<fname<<endl;
      return;
    }
  while(in>>gt[counter].event>>gt[counter].label>>gt[counter].code
	>>gt[counter].px>>gt[counter].py>>gt[counter].pz>>gt[counter].pt>>gt[counter].eta>>gt[counter].nhits)
    counter++;
  
  char filename[100];
  file = new AliHLTMemHandler();
}
#endif
