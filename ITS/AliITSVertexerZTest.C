#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TGeoManager.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <Riostream.h>
#include <AliRun.h>
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliGeomManager.h>
#include <AliITSVertexerZ.h>
#include <AliRunLoader.h>
#include <AliITSLoader.h>

#endif

void AliITSVertexerZTest(Float_t delphi=0.05,Float_t window=3.,Float_t initx=0., Float_t inity=0.,TString FileIn="galice.root"){
  // delphi ---> azimuthal range to accept tracklets
  // window ---> window in Z around the peak of tracklets proj. in mm
  Int_t kDebug = 50;
  TH2F *diff2 = new TH2F("diff2","Zfound vs Ztrue",100,-6,6,100,-6,6); 
  TH1F *diff1 = new TH1F("diff1","Zfound - Ztrue(#mu m)",100,-500,500);
  delete gAlice;
  gAlice = 0;
  AliRunLoader *rl = AliRunLoader::Open(FileIn.Data());
  Int_t retval = rl->LoadgAlice();
  if(retval){
    cerr<<"AliITSVertexerZTest.C: AliRun object not found"<<endl;
    return;
  }
  retval = rl->LoadHeader();
  if (retval){
    cerr<<"AliITSVertexerZTest.C : LoadHeader returned error"<<endl;
    return;
  }
  retval = rl->LoadKinematics();
  if (retval){
    cerr<<"AliITSVertexerZTest.C : LoadKinematics returned error"<<endl;
    return;
  }
  
  AliGeomManager::LoadGeometry("geometry.root");

  AliITSVertexerZ *dovert = new AliITSVertexerZ(initx,inity);
  dovert->Init("default");
  //dovert->SetDebug(0);
  //  dovert->SetDiffPhiMax(delphi);
  //  dovert->SetWindow(window);
  dovert->PrintStatus();
  Int_t sigmazero=0;
  AliESDVertex *vert = 0;
  AliITSLoader* itsloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  itsloader->LoadRecPoints("read");
  for(Int_t i=0; i<rl->TreeE()->GetEntries(); i++){
    rl->GetEvent(i);
    // The true Z coord. is fetched for comparison
    AliHeader *header = rl->GetHeader();
    AliGenEventHeader* genEventHeader = header->GenEventHeader();
    TArrayF primaryVertex(3);
    genEventHeader->PrimaryVertex(primaryVertex);
    TTree* cltree = itsloader->TreeR();
    vert = dovert->FindVertexForCurrentEvent(cltree);
    if(kDebug>0){
      // Prints the results
      cout <<"========================================================\n";
      cout << "Event number: "<<i<<")  Z Vertex:"<<endl;
      if(vert){
	cout<<"FOUND: "<<vert->GetZv()<<"; ";
	cout<<vert->GetZRes()<<"; "<<vert->GetNContributors()<<endl;
	cout <<" True Z position "<<primaryVertex[2]<<", diff= ";
	cout<<(primaryVertex[2]-vert->GetZv())*10000.<<endl;
      } else {
	cout<<"NOT FOUND"<<endl;
      }
    }
    if(vert){
      Float_t found = vert->GetZv();
      diff2->Fill(primaryVertex[2],found);
      found = 10000.*(found-primaryVertex[2]);
      if(vert->GetZRes()!=0){
	diff1->Fill(found);
      } else {
	sigmazero++;
      }
      dovert->WriteCurrentVertex();
    }
  }
  if(kDebug>0){
    cout<<"Only one tracklet (sigma = 0) "<<sigmazero<<endl;
  }
  
}
