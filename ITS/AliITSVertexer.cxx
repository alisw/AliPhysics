#include <AliESDVertex.h>
#include <AliITSVertexer.h>
#include <AliRunLoader.h>
#include <AliITSLoader.h>

ClassImp(AliITSVertexer)

//////////////////////////////////////////////////////////////////////
// Base class for primary vertex reconstruction                     //
// AliESDVertexer is a class for full 3D primary vertex finding     //
// derived classes: AliITSVertexerIons AliITSvertexerPPZ            //
//                  AliITSVertexerTracks                            //
//////////////////////////////////////////////////////////////////////

//______________________________________________________________________
AliITSVertexer::AliITSVertexer():AliVertexer() {
  // Default Constructor
  //SetUseV2Clusters(kTRUE);
}

AliITSVertexer::AliITSVertexer(TString filename) {
  // Standard constructor
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  if(!rl){
    Fatal("AliITSVertexer","Run Loader not found");
  }
  if(rl->LoadgAlice()){
    Fatal("AliITSVertexer","The AliRun object is not available - nothing done");
  }
  fCurrentVertex  = 0;   
  SetDebug();
  SetFirstEvent(0);
  SetLastEvent(0);
  rl->LoadHeader();
  AliITSLoader* itsLoader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  if(!filename.Contains("default"))itsLoader->SetVerticesFileName(filename);
  if(!filename.Contains("null"))itsLoader->LoadVertices("recreate");
  itsLoader->LoadRecPoints();
  Int_t lst;
  if(rl->TreeE()){
    lst = static_cast<Int_t>(rl->TreeE()->GetEntries());
    SetLastEvent(lst-1);
  }
  //SetUseV2Clusters(kTRUE);
}

//______________________________________________________________________
AliITSVertexer::AliITSVertexer(const AliITSVertexer &vtxr) : AliVertexer(vtxr) {
  // Copy constructor
  // Copies are not allowed. The method is protected to avoid misuse.
  Error("AliITSVertexer","Copy constructor not allowed\n");
}

//______________________________________________________________________
AliITSVertexer& AliITSVertexer::operator=(const AliITSVertexer& /* vtxr */){
  // Assignment operator
  // Assignment is not allowed. The method is protected to avoid misuse.
  Error("= operator","Assignment operator not allowed\n");
  return *this;
}


//______________________________________________________________________
void AliITSVertexer::WriteCurrentVertex(){
  // Write the current AliVertex object to file fOutFile
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliITSLoader* itsLoader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  fCurrentVertex->SetName("Vertex");
  //  const char * name = fCurrentVertex->GetName();
  //  itsLoader->SetVerticesContName(name);
  Int_t rc = itsLoader->PostVertex(fCurrentVertex);
  rc = itsLoader->WriteVertices();
}
/*
//______________________________________________________________________
void AliITSVertexer::Clusters2RecPoints
(const TClonesArray *clusters, Int_t idx, TClonesArray *points) {
  //------------------------------------------------------------
  // Conversion AliITSRecPoint -> AliITSRecPoints for the ITS
  // module "idx" (entry in the tree with the clusters).
  // Simplified version, supposed to work with the pixels only !
  //------------------------------------------------------------
  const Int_t klastSPD1=79; //let's hope the number of the SPDs will not change
  const Int_t klastSPD2=239;//let's hope the number of the SPDs will not change

  Float_t yshift = 0; //see AliITSclustererV2.cxx about these shifts
  Float_t zshift[4] = {-10.708000, -3.536000, 3.536000, 10.708000}; //let's hope the positioning of the SPDs will not change

  if (idx<=klastSPD1) {
    yshift=0.248499;  //let's hope the positioning of the SPDs will not change
  } else if (idx<=klastSPD2) {
    yshift=3.096207;  //let's hope the positioning of the SPDs will not change
  } else {
    Fatal("Clusters2RecPoints","This is not an SPD module ! %d",idx);
  }

  TClonesArray &pn=*points;
  Int_t ncl=clusters->GetEntriesFast();
  for (Int_t i=0; i<ncl; i++) {
    AliITSRecPoint p;
    AliITSRecPoint *c = (AliITSRecPoint *)clusters->UncheckedAt(i);

    Float_t x=c->GetY();  if (idx<=klastSPD1) x=-x;
    x+=yshift;

    Float_t z=c->GetZ();
    z=-z; z+=zshift[idx%4];

    p.SetDetLocalX(x);
    p.SetDetLocalZ(z);
    p.SetQ(c->GetQ());
    p.SetSigmaDetLocX2(c->GetSigmaY2());
    p.SetSigmaZ2(c->GetSigmaZ2());
    p.SetLabel(c->GetLabel(0),0);
    p.SetLabel(c->GetLabel(1),1);
    p.SetLabel(c->GetLabel(2),2);

    new (pn[i]) AliITSRecPoint(p);
  }

}

*/
