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
#include <AliITSLoader.h>

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
 
  Int_t retval = rl->LoadHeader();
  if (retval){
    cerr<<"AliITSFindPrimaryVertex.C : LoadHeader returned error"<<endl;
    return;
  }
  retval = rl->LoadKinematics();
  if (retval){
    cerr<<"AliITSFindPrimaryVertex.C : LoadKinematics returned error"<<endl;
    return;
  }

  // Open output file for vertices (default name: ITS.Vertex.root 
  // and Create vertexer
  AliITSVertexerIons *vertexer = new AliITSVertexerIons("default");
  AliITSVertex *V;
  //   Loop over events 
  //
 
  for (int nev=evNumber1; nev< evNumber2; nev++) {
    cout<<"=============================================================\n";
    cout<<" Processing event "<<nev<<endl;
    cout<<"=============================================================\n";
    rl->GetEvent(nev);
    cout << "nev         " << nev <<endl;
    // The true Z coord. is fetched for comparison
    AliHeader *header = rl->GetHeader();
    AliGenEventHeader* genEventHeader = header->GenEventHeader();
    TArrayF primaryVertex(3);
    genEventHeader->PrimaryVertex(primaryVertex);
  
    TStopwatch timer;
    timer.Start();

    V=vertexer->FindVertexForCurrentEvent(nev);
    if(V){
      Double_t pos[3];
      for(Int_t kk=0;kk<3;kk++)pos[kk]=(Double_t)primaryVertex[kk];
      V->SetTruePos(pos);
    }
    timer.Stop();
    timer.Print();
    if(!V)continue;
    cout << endl << "Xv = " << V->GetXv() << " cm" << endl;
    cout << "X resolution = " << V->GetXRes()*10000 << " microns"  << endl;
    cout << "Signal/Noise for X = " << V->GetXSNR() << endl;
    cout << endl << "Yv = " << V->GetYv() << " cm"  << endl;
    cout << "Y resolution = " << V->GetYRes()*10000 << " microns"  << endl;
    cout << "Signal/Noise for Y = " << V->GetYSNR() << endl;
    cout << endl << "Zv = " << V->GetZv() << " cm" << endl;
    cout << "Z Resolution = " << V->GetZRes()*10000 << " microns" << endl;
    cout << "Signal/Noise for Z = " << V->GetZSNR() <<endl;
	 	 
    vertexer->WriteCurrentVertex(); 
  } 
  

}

