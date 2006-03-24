#ifndef __CINT__
  #include "alles.h"
  #include "AliTPCtracker.h"
  #include "TView.h"
  #include "TPolyMarker3D.h"
  #include "AliSimDigits.h"
  #include "AliRunLoader.h"
  #include "AliLoader.h"
  #include "AliTPCParamSR.h"
#endif

//  
//  display 3D digits
//  input parameter is event number  - threshol to the noise
//  if sdigits=kTRUE it will display Sdigits - otherwise it display digits
//  signal event is displayed with yellow color

Int_t AliTPCDisplayDigits3Dnew(Int_t eventn=0, Int_t noiseth=15, Bool_t sdigits=kFALSE) {
   cerr<<"Displaying digits...\n";

//   TFile *file=TFile::Open("galice.root");
//   if (!file->IsOpen()) {cerr<<"Can't open galice.root !\n"; return 1;}

//   TFile *cf=TFile::Open("galice.root");
//    if (!cf->IsOpen()){cerr<<"Can't open AliTPCclusters.root !\n"; return 3;}

   AliRunLoader* rl = AliRunLoader::Open();
   rl->GetEvent(eventn);

   AliLoader* tpcl = (AliLoader*)rl->GetLoader("TPCLoader");
   if (tpcl == 0x0)
    {
      cerr<<"Can not get TPC Loader"<<endl;
      return 1;
    }

   rl->CdGAFile();
   AliTPCParam *param=(AliTPCParam *)gDirectory->Get("75x40_100x60");
   if(param){
     cerr<<"2 pad-length geom hits with 3 pad-lengths geom parameters\n";
     delete param;
     param = new AliTPCParamSR();
   }
   else
   {
     param=(AliTPCParamSR *)gDirectory->Get("75x40_100x60_150x60");
   }

   if (!param) {cerr<<"TPC parameters have not been found !\n"; return 2;}

   TCanvas *c1=new TCanvas("ddisplay", "Digits display",0,0,700,730);
   TView *v=new TView(1);
   v->SetRange(-430,-560,-430,430,560,1710);
   c1->Clear();
   c1->SetFillColor(1);
   c1->SetTheta(90.);
   c1->SetPhi(0.);

   AliTPCDigitsArray *digarr=new AliTPCDigitsArray;
   digarr->Setup(param);
   
  
  TTree* tree;
   
   char  cname[100];
   if (!sdigits)
     {
       tpcl->LoadDigits();
       tree = tpcl->TreeD();
     }
   else
     {
       tpcl->LoadSDigits();
       tree = tpcl->TreeS();
     }

// some "constants"
   Int_t markerColorSignal = 5;
   Int_t markerColorBgr = 2;
   Int_t MASK = 10000000;

   Int_t imarBgr=0;
   Int_t imarSignal=0;
   Int_t ncl=0;
   
   digarr->ConnectTree(tree);
   
   Int_t nrows=Int_t(digarr->GetTree()->GetEntries());
   Int_t all0=0;
   for (Int_t n=0; n<nrows; n++) {
       AliSimDigits *s=(AliSimDigits*)digarr->LoadEntry(n);
       Int_t sec,row;
       param->AdjustSectorRow(s->GetID(),sec,row);
       Int_t npads, sign;
       {
       const Int_t kNIS=param->GetNInnerSector(), kNOS=param->GetNOuterSector();
       if (sec < kNIS) {
         npads = param->GetNPadsLow(row);
         sign = (sec < kNIS/2) ? 1 : -1;
       } else {
         npads = param->GetNPadsUp(row);
         sign = ((sec-kNIS) < kNOS/2) ? 1 : -1;
       }
       }

       AliSimDigits *digrow = (AliSimDigits*)digarr->GetRow(sec,row);
       ncl=0;
       if (digrow->First()){
       while(digrow->Next()) ncl++;
       ncl++;
       }
       TPolyMarker3D *pm=new TPolyMarker3D(ncl);
       TPolyMarker3D *pmSignal=new TPolyMarker3D(ncl); // polymarker for signal
       imarBgr=0;
       imarSignal=0;
       if (digrow->First()) do {
       Short_t dig=digrow->CurrentDigit();
       Double_t y = (digrow->CurrentColumn()- 0.5 - 0.5*npads)*param->GetPadPitchWidth(sec);
       Double_t z = sign*(param->GetZLength()-param->GetZWidth()*digrow->CurrentRow());       
       Double_t x=param->GetPadRowRadii(sec,row);
//       cout<<dig<<endl;
       if (dig<noiseth) continue;
//       cout<<"\nAbove noise Threshold";
       Int_t trackId = digrow->GetTrackID(digrow->CurrentRow(),digrow->CurrentColumn(),0);      
       Float_t cs, sn, tmp;
       param->AdjustCosSin(sec,cs,sn);
       tmp = x*cs-y*sn; y= x*sn+y*cs; x=tmp;
       if (trackId<0) continue;  //if signal from croostalk - we don't have track ID information
//       cout<<"Track ID > 0";
             if (trackId < MASK-1) {
         pmSignal->SetPoint(imarSignal,x,y,z);
         imarSignal++;
       } else {
         pm->SetPoint(imarBgr,x,y,z);
         imarBgr++;
       }
       } while (digrow->Next());

// change color for signal
       digarr->ClearRow(sec,row);
       pm->SetMarkerSize(1); 
       pm->SetMarkerColor(markerColorBgr);
       pm->SetMarkerStyle(1);
       pm->Draw();

       pmSignal->SetMarkerSize(1); 
       pmSignal->SetMarkerColor(markerColorSignal);
       pmSignal->SetMarkerStyle(1);
       pmSignal->Draw();
       all0+=imarSignal;
//       cout<<"imarSignal ="<<imarSignal<<"   imarBgr ="<<imarBgr<<"   ncl ="<<ncl<<endl;
   }
   printf("%d\n",all0);
   
   

   delete digarr;

   //draw TPC skeleton
   rl->CdGAFile();
   TGeometry *geom=(TGeometry*)gDirectory->Get("AliceGeom");
   TNode * main = (TNode*)((geom->GetListOfNodes())->First());
   TIter next(main->GetListOfNodes());
   TNode  *module=0;
   while((module = (TNode*)next())) {
     char ch[100];
     sprintf(ch,"%s\n",module->GetTitle());
     //printf("%s\n",module->GetTitle());
     if (ch[0]=='T'&&ch[1]=='P' && ch[2]=='C')  //if TPC draw
       module->SetVisibility(3);
     else
       module->SetVisibility(-1);
   }
     
   
   geom->Draw("same");
   //v->Draw();
   c1->Modified(); c1->Update(); 
   
   delete rl;
   return 0;
}

