#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TClassTable.h>
#include <TFile.h>
#include <TTree.h>
#include <Riostream.h>
#include <AliRun.h>
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliITSVertexerIons.h>
#include <AliRunLoader.h>
#include <AliITSVertexerIons.h>
#include <AliITSLoader.h>
#include <unistd.h>

#endif

void AliITSFindPrimaryVertex(Int_t evNumber1=0,Int_t NumbofEv=1, const char *filename="galice.root") {

  Int_t evNumber2 = evNumber1+NumbofEv;
  
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->Macro("loadlibs.C");
  } else {
    if(gAlice){
      delete gAlice->GetRunLoader();
      delete gAlice;
      gAlice=0;
    }
  }
  
  AliRunLoader* rl = AliRunLoader::Open("galice.root");
  if (rl == 0x0){
    ::Error("AliITSFindPrimaryVertex.C","Can not open session RL=NULL");
    return;
  }
  
  // Open output file for vertices (default name: ITS.Vertex.root 
  // and Create vertexer

  AliITSVertexerIons *vertexer = new AliITSVertexerIons("default");
  //vertexer->SetDebug(1);
  
  AliESDVertex *V;
  //   Loop over events 
   
  for (int nev=evNumber1; nev< evNumber2; nev++) {
    cout<<"=============================================================\n";
    cout<<" Processing event "<<nev<<endl;
    cout<<"=============================================================\n";
    rl->GetEvent(nev);
    AliHeader *header = rl->GetHeader();
    AliGenEventHeader* genEventHeader = header->GenEventHeader();
    TArrayF primaryVertex(3);
    genEventHeader->PrimaryVertex(primaryVertex);
    
    AliGenHijingEventHeader* hijingHeader = (AliGenHijingEventHeader*)  genEventHeader;
    Float_t b = hijingHeader->ImpactParameter();   
    cout << "Impact parameter = " << b << " fm" << endl;

    TStopwatch timer;
    timer.Start();

    V=vertexer->FindVertexForCurrentEvent(nev);

    TVector3 vtrue(primaryVertex[0],primaryVertex[1],primaryVertex[2]);
    TVector3 vfound(V->GetXv(),V->GetYv(),V->GetZv());
    TVector3 dif=vtrue-vfound;
    cout << "True vertex coordinates (cm) = " << vtrue.X() << " " << vtrue.Y() << " " << vtrue.Z() << endl;
    cout << "Found vertex coordinates  (cm) = " << vfound.X() << " " << vfound.Y() << " " << vfound.Z() << endl;    cout << "Difference true - found (cm) = " << dif.Mag() << " " << dif.X() << " " << dif.Y() << " " << dif.Z() << endl;
    
    if(V){
      Double_t pos[3];
      for(Int_t kk=0;kk<3;kk++)pos[kk]=(Double_t)primaryVertex[kk];
      V->SetTruePos(pos);
    }
    timer.Stop();
    timer.Print();
    
    vertexer->WriteCurrentVertex(); 
  } 
  

}

