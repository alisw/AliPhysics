#ifndef __CINT__
  #include "alles.h"
  #include "AliTPCtracker.h"
  #include "TView.h"
  #include "TPolyMarker3D.h"

#endif
Int_t AliTPCDisplayDigits3D(Int_t eventn=0, Int_t noiseth=15) {
   cerr<<"Displaying digits...\n";

   TFile *file=TFile::Open("galice.root");
   if (!file->IsOpen()) {cerr<<"Can't open galice.root !\n"; return 1;}

   TFile *cf=TFile::Open("galice.root");
   // if (!cf->IsOpen()){cerr<<"Can't open AliTPCclusters.root !\n"; return 3;}

   AliTPCParam *param=(AliTPCParam *)cf->Get("75x40_100x60");
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
   char  cname[100];
   sprintf(cname,"TreeD_75x40_100x60_%d",eventn);

   digarr->ConnectTree(cname);
   Int_t nrows=Int_t(digarr->GetTree()->GetEntries());
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
       Int_t ncl=0;
       if (digrow->First()){
	 while(digrow->Next()) ncl++;
	 ncl++;
       }
       TPolyMarker3D *pm=new TPolyMarker3D(ncl);
       digrow->First();
       Int_t imar=0;
       do {
	 Short_t dig=digrow->CurrentDigit();
	 Double_t y = (digrow->CurrentColumn()- 0.5 - 0.5*npads)*param->GetPadPitchWidth(sec);
	 Double_t z = sign*(param->GetZLength()-param->GetZWidth()*digrow->CurrentRow());	 
	 Double_t x=param->GetPadRowRadii(sec,row);
	 if (dig<noiseth) continue;
	 Float_t cs, sn, tmp;
	 param->AdjustCosSin(sec,cs,sn);
	 tmp = x*cs-y*sn; y= x*sn+y*cs; x=tmp;
	 pm->SetPoint(imar,x,y,z);
	 imar++;
       } while (digrow->Next());
       digarr->ClearRow(sec,row);
       pm->SetMarkerSize(1); pm->SetMarkerColor(2); pm->SetMarkerStyle(1);
       pm->Draw();
   }
   delete digarr;
   cf->Close();
   //draw TPC
   TGeometry *geom=(TGeometry*)file->Get("AliceGeom");
   TList *list = geom->GetListOfNodes();
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

   file->Close();
   return 0;
}

