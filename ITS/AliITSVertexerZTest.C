#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <Riostream.h>
#include <AliRun.h>
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliITSVertexerPPZ.h>

#endif

void AliITSVertexerZTest(Float_t delphi=0.05,Float_t window=3.,Float_t initx=0., Float_t inity=0.,TString FileWithKine="galice.root", TString FileWithRecP="galice.root"){
  // delphi ---> azimuthal range to accept tracklets
  // window ---> window in Z around the peak of tracklets proj. in mm
  Int_t kDebug = 0;
  TH2F *diff2 = new TH2F("diff2","Zfound vs Ztrue",100,-6,6,100,-6,6); 
  TH1F *diff1 = new TH1F("diff1","Zfound - Ztrue(#mu m)",100,-500,500);
  delete gAlice;
  gAlice = 0;
  TFile *in = new TFile("galice.root");
  gAlice = (AliRun*)in->Get("gAlice");
  if(FileWithKine != FileWithRecP)gAlice->SetTreeRFileName(FileWithRecP);
  TFile *fo = new TFile("vertici.root","recreate");
  AliITSVertexerPPZ *dovert = new AliITSVertexerPPZ(in,fo,initx,inity);
  dovert->SetDebug(0);
  dovert->SetDiffPhiMax(delphi);
  dovert->SetWindow(window);
  dovert->PrintStatus();
  Int_t meno100=0;
  Int_t meno200=0;
  Int_t meno110=0;
  Int_t sigmazero=0;
  AliITSVertex *vert = 0;
  for(Int_t i=0; i<gAlice->TreeE()->GetEntries(); i++){
    gAlice->GetEvent(i);
    // The true Z coord. is fetched for comparison
    AliHeader *header = gAlice->GetHeader();
    AliGenEventHeader* genEventHeader = header->GenEventHeader();
    TArrayF primaryVertex(3);
    genEventHeader->PrimaryVertex(primaryVertex);
    vert = dovert->FindVertexForCurrentEvent(i);
    if(kDebug>0){
      // Prints the results
      cout <<"========================================================\n";
      cout << "Event number: "<<i<<")  Z Vertex:"<<endl;
      if(vert){
	cout<<"FOUND: "<<vert->GetZv()<<"; ";
	cout<<vert->GetZRes()<<"; "<<vert->GetNContributors()<<endl;
      }
      else {
	cout<<"NOT FOUND - fZFound= "<<dovert->GetZFound();
	cout<<" fZsig= "<<dovert->GetZsig()<<endl;
	if(dovert->GetZFound() == -100) meno100++;
	if(dovert->GetZFound() == -200) meno200++;
	if(dovert->GetZFound() == -110) meno110++;
      }
      cout <<" True Z position "<<primaryVertex[2]<<", diff= ";
      cout<<(primaryVertex[2]-dovert->GetZFound())*10000.<<endl;
    }
    if(vert){
      Double_t pos[3];
      for(Int_t kk=0;kk<3;kk++)pos[kk]=(Double_t)primaryVertex[kk];
      vert->SetTruePos(pos);
      Float_t found = vert->GetZv();
      diff2->Fill(primaryVertex[2],found);
      found = 10000.*(found-primaryVertex[2]);
      if(vert->GetZRes()!=0){
	diff1->Fill(found);
      }
      else {
	sigmazero++;
      }
      dovert->WriteCurrentVertex();
    }
  }
  in->Close();
  fo->Close();
  delete in;
  delete fo;
  if(kDebug>0){
    cout<<"Number of bad vertices with code  -100: "<<meno100<<endl;
    cout<<"Number of bad vertices with code  -110: "<<meno110<<endl;
    cout<<"Number of bad vertices with code  -200: "<<meno200<<endl;
    cout<<"Only one tracklet (sigma = 0) "<<sigmazero<<endl;
  }

}
