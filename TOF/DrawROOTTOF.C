#ifndef __CINT__
  #include "TView.h"
  #include "TPolyMarker3D.h"

#endif

Int_t DrawROOTTOF() {
  //
  // author: F. Pierella
  // report bug to pierella@bo.infn.it
  //
  // use case
  // generate an event with TOF included
  // start aliroot
  // .x DrawROOTTOF.C
  cerr<<"ROOT TOF Geometry...\n";
  cerr<<"(TRD and TPC are included)\n";
  TFile *file=TFile::Open("galice.root");
  if (!file->IsOpen()) {cerr<<"Can't open galice.root !\n"; return 1;}
  
  TCanvas *c1=new TCanvas("ddisplay", "TOF display",0,0,700,730);
  TView *v=new TView(1);
  v->SetRange(-430,-560,-430,430,560,1710);
  c1->Clear();
  c1->SetFillColor(10);
  c1->SetTheta(90.);
  c1->SetPhi(0.);
  
  
  //draw TOF with TRD and TPC included
  TGeometry *geom=(TGeometry*)file->Get("AliceGeom");
  TList *list = geom->GetListOfNodes();
  TNode * main = (TNode*)((geom->GetListOfNodes())->First());
  TIter next(main->GetListOfNodes());
  TNode  *module=0;
  while((module = (TNode*)next())) {
    char ch[100];
    sprintf(ch,"%s\n",module->GetTitle());
    //printf("%s\n",module->GetTitle());
    if ((ch[0]=='F'&&ch[1]=='T' && ch[2]=='O') || (ch[0]=='T'&&ch[1]=='R' && ch[2]=='D') || (ch[0]=='T'&&ch[1]=='P' && ch[2]=='C')){  //if TOF or TPC or TRD draw
      module->SetVisibility(3);
    }else{
      module->SetVisibility(-1);
    }
  }
  
  c1->cd();
  geom->Draw("same");
  //v->Draw();
  c1->Modified(); c1->Update(); 
  
  //
  //  Draw the geometry using the x3d viewver.

  c1->x3d();
  //
  // once in x3d viewer, type m to see the menu.
  // For example typing r will show a solid model of this geometry.
  
  file->Close();
  return 0;
}
