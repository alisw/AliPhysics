#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"

#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODCluster.h"

#endif

void ReadAOD(const char *fileName = "AliAOD.root") {

  // open input file and get the TTree
  TFile inFile(fileName, "READ");
  if (!inFile.IsOpen()) return;

  TTree *aodTree = (TTree*)inFile.Get("aodTree");

  AliAODEvent *ev = new AliAODEvent();
  ev->ReadFromTree(aodTree);

  // loop over events
  Int_t nEvents = aodTree->GetEntries();
  for (Int_t nEv = 0; nEv < nEvents; nEv++) {
    cout << "Event: " << nEv+1 << "/" << nEvents << endl;

    // read events
    aodTree->GetEvent(nEv);

    //print event info
    ev->GetHeader()->Print();

    // loop over tracks
    Int_t nTracks = ev->GetNTracks();
    for (Int_t nTr = 0; nTr < nTracks; nTr++) {
      
      AliAODTrack *tr = ev->GetTrack(nTr);

      // print track info
      cout << nTr+1 << "/" << nTracks << ": track pt: " << tr->Pt();
      if (tr->GetProdVertex()) {
	cout << ", vertex z of this track: " << tr->GetProdVertex()->GetZ();
      }
      cout << endl;
    }

    // loop over vertices
    Int_t nVtxs = ev->GetNVertices();
    for (Int_t nVtx = 0; nVtx < nVtxs; nVtx++) {
      
      // print track info
      cout << nVtx+1 << "/" << nVtxs << ": vertex z position: " << ev->GetVertex(nVtx)->GetZ() << endl;
    }
  }
  
  return;
}
