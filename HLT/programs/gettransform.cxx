//#define DOTESTS /* undefine this for getting the parameters */

#ifdef DOTESTS
#define use_root
#endif

#include <stream.h>
#include <libgen.h>
#include <math.h>
#include <AliL3RootTypes.h>
#include <AliL3Transform.h>
#include <AliL3Logging.h>
#include <AliL3Logger.h>
#include <AliL3MemHandler.h>

#ifdef use_root
#include <TApplication.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1.h>
#include <TMath.h>
#endif

/**
   This program extracts parameters and lookup tables needed for the
   vhdl implementation of the Hough transform. 
*/

int main(int argc,char **argv)
{
  Int_t patch=0;
  Int_t slice=0;
  Char_t path[1000];

  AliL3Logger l;
  l.Set(AliL3Logger::kAll);
  l.UseStderr();
  //l.UseStdout();
  //l.UseStream();
    
  if (argc>1) {
    slice=atoi(argv[1]);
  }
  if (argc>2) {
    patch=atoi(argv[2]);
  }  
  if (argc>3) {
    strcpy(path,argv[3]);
  } else strcpy(path,"/tmp/data/RawData/slice0");
  if(argc>4){
    cout<<"Usage: transform [slice] [patch] [path]"<<endl;
    exit(1);
  }

  AliL3Transform::Init(path);
  //cerr << "Transform version: " << AliL3Transform::GetVersion() << endl;

  Int_t npads=0;
  Int_t sector=0;
  Int_t sector_row=0;

  Float_t Xt[200];
  Float_t Yt[200];


#ifdef DOTESTS /* do some tests in root histos */
#ifdef use_root
  Float_t xbins[126];
  for(int t=0;t<126;t++){
    Float_t eta=exp(2.*1./125*t);
    xbins[t]=(eta-1)/(1+eta);
    //cout << eta << " " << xbins[t] << endl;
  }

  TApplication theApp("App", &argc, argv);
  
  TCanvas *c = new TCanvas("c", "The Histogram Canvas", 800, 800);
  TPad *pad1 = new TPad("pad1","Pad for Eta",0.05,0.50,0.95,0.95,21);
  TPad *pad2 = new TPad("pad2","Pad for Tests",0.05,0.05,0.95,0.45,21);
  //  TH1F *h1f = new TH1F("h1f","Eta Distribution",125,0,1.);
  TH1F *h1f = new TH1F("h1f","Test Error",100,90,110);
  TH1F *h2f = new TH1F("h2f","Test Distribution",1000,0,10);
  //  TH1F *h2f = new TH1F("h2f","Test Distribution",125,xbins);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
#endif

  Int_t ntime=AliL3Transform::GetNTimeBins();
  Float_t xyz[3]={0,0,0};
  Float_t rpe[3]={0,0,0};
  int rebinned=0;
  for(Int_t rr=AliL3Transform::GetFirstRow(patch);rr<=AliL3Transform::GetLastRow(patch);rr++){
    npads=AliL3Transform::GetNPads(rr);
    AliL3Transform::Slice2Sector(slice,rr,sector,sector_row);
    for(Int_t pp=0;pp<npads;pp++){
      for (Int_t tt=0;tt<ntime;tt++){
	AliL3Transform::Raw2Local(xyz,sector,sector_row,pp,tt);
	AliL3Transform::XYZtoRPhiEta(rpe,xyz);
#ifdef use_root
	if(!rebinned) {
	  h2f->SetBins(500,0,252.*252/xyz[0]/xyz[0]);
	  rebinned=1;
	}
	//h1f->Fill(rpe[2]);
	//h2f->Fill(xyz[2]/rpe[0]);
	Float_t A2=xyz[2]*xyz[2]/(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
	Float_t Sroot=1./TMath::Sqrt(1.+A2);
	Int_t bin=h2f->FindBin(A2);
	Double_t binCenter = h2f->GetBinCenter(bin);
	Double_t binCenterval = 1./TMath::Sqrt(1.+binCenter);
	Float_t err=(Sroot*100./binCenterval);
	h1f->Fill(err);
	h2f->Fill(A2);
	//printf("%.3f %.3f %.3f\n",Sroot,binCenter,err);
#endif
	//printf("%d %d %d %f %f %f %f %f %f\n",rr,pp,tt,xyz[0],xyz[1],xyz[2],rpe[0],rpe[1],rpe[2]);
      }
    }
  }
#else /* do the extraction of the parameters */
  Int_t firstrow=AliL3Transform::GetFirstRow(patch);
  Int_t lastrow=AliL3Transform::GetLastRow(patch);
  AliL3Transform::Slice2Sector(slice,firstrow,sector,sector_row);

  Float_t ytabval=0;
  //Int_t maxrow=NRows[patch][1]-NRows[patch][0];
  Int_t maxrow=AliL3Transform::GetNRows(patch);  
  Float_t padpitch=0;
  if(sector<AliL3Transform::GetNSectorLow())
    padpitch=AliL3Transform::GetPadPitchWidthLow();
  else
    padpitch=AliL3Transform::GetPadPitchWidthUp();  

  printf("MinRow: %d\nMaxRow: %d\n",0,maxrow);
  printf("YPadWidth: %.2f\n",padpitch);
  printf("ZSign: %d\n",slice < 18 ? 1:-1);
  //printf("ZWidth: %.2f\n",AliL3Transform::GetZWidth());
  printf("ZWidth: %.2f\n",AliL3Transform::GetZLength()+AliL3Transform::GetZOffset());
  printf("TimeWidth: %.2f\n\n",AliL3Transform::GetZWidth());

  //calculating lookup tables for slice and patch!
  for(Int_t rr=firstrow;rr<=lastrow;rr++){
    npads=AliL3Transform::GetNPads(rr);
    
    //Y(row,pad)=pad*padpitch-ytabval(row);
    ytabval=0.5*(npads-1)*padpitch;

    Xt[rr-firstrow]=AliL3Transform::Row2X(rr);
    Yt[rr-firstrow]=ytabval;
    //row in patch: X(row) Y_part(row)
    printf("Row: %d X: %.2f Y: %.2f\n",rr-firstrow,AliL3Transform::Row2X(rr),ytabval);
  }
  printf("\n\nVHDL-Output for LUT:\nX_table := (");
  for(int i=0;i<maxrow-1;i++){
    printf("%.2f, ",Xt[i]);
  }
  printf("%.2f);\n",Xt[maxrow-1]);
  printf("Y_table := (");
  for(int i=0;i<maxrow-1;i++){
    printf("%.2f, ",Yt[i]);
  }
  printf("%.2f);\n",Yt[maxrow-1]);
#endif

#ifdef DOTESTS 
#ifdef use_root
  h1f->Draw();
  pad2->cd();
  h2f->Draw();
  c->Update();

  theApp.Run();
#endif
#endif

  exit(0);
}
  
//#############################################################

#if 0
  AliL3MemHandler file; //Does all the file/data handling

  //Open the data file:
  if(!file.SetBinaryInput(digitfile))
    {
      cerr<<"Error opening file "<< digitfile <<endl;
      return -1;
    }
  
  //Store the data in memory, and get the pointer to it:
  unsigned int ndigits=0;
  AliL3DigitRowData *digits=0;
  digits=(AliL3DigitRowData*)file.CompBinary2Memory(ndigits);
  file.CloseBinaryInput();
  
  /*  
  for(Int_t r=0; r<31; r++) {
    UInt_t padrow=digits->fRow;
    AliL3DigitData *dPt = (AliL3DigitData*)digits->fDigitData;
    cout<<"padrow "<<padrow<<" ndigits "<<digits->fNDigit<<endl;
	  
    for(Int_t d=0; d<digits->fNDigit; d++) {
      cout<<" padrow "<<padrow<<" pad "<<(int)dPt[d]->fPad<<" time " <<(int)dPt[d]->fTime<<" charge "<<(int)dPt[d]->fCharge<<endl;
    }
    file->UpdateRowPointer(digits);
  }
  return 100;
  */
    
  digits=(AliL3DigitRowData*)file->GetDataPointer(ndigits);
  AliL3Transform transform; //Storing all detector-spesific quantities, needed by the clusterfinder.
  AliL3ClustFinderNew cf(&transform); //The cluster finder itself.

  //Switch off deconvolution:
  cf.SetDeconv(false);
  cf.SetXYError(0.2);
  cf.SetZError(0.3);

  //Init cluster finder
  cf.InitSlice(0,0,0,20,10000);
  
  //Give the data pointer to the cluster finder
  cf.Read(ndigits,digits);
  cout << digits << endl;
  //Start processing:
  cf.ProcessDigits();
  
/* 
 int slice=1;
  a = new AliL3Hough(digitfile,kTRUE,n_eta_segments);
  
  a->ReadData(slice);
  
  a->Transform();
  a->AddAllHistograms();
    
  c1 = new TCanvas("c1","",2);
  a->GetTransformer(0)->GetHistogram(0)->Draw("box");
  
  
  a->FindTrackCandidates();
  //a->Evaluate(2);
 */ 
  /*
  //Display the track:
  dhist = new AliL3Histogram("dhist","",250,0,250,250,-125,125);
  for(int i=0; i<6; i++)
    a->GetEval(i)->DisplayEtaSlice(0,dhist);
  c2 = new TCanvas("c2","",2);
  dhist->Draw();
  */
  
  //a->WriteTracks();//"../compression/");
    
/*
  tracks = (AliL3TrackArray*)a->GetTracks(0);
   
  float xyz[3];
  int row=0;
  cout<<"Found "<<tracks->GetNTracks()<<" tracks"<<endl;
  for(int i=0; i<tracks->GetNTracks(); i++)
    {
      AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      //track->GetCrossingPoint(row,xyz);
      //transform->Local2Raw(xyz,1,row);
      //cout<<"Recon. pad "<<(int)xyz[1]<<" time "<<(int)xyz[2]<<endl;
      cout<<"Pt "<<track->GetPt()<<" phi "<<track->GetPhi0()<<" charge "<<track->GetCharge()<<" weigth "<<track->GetWeight()<<endl;
    }
*/
  
  return;
}

#endif  
