#include <AliESDVertex.h>
#include <AliITSVertexer.h>
#include <AliRunLoader.h>
#include <AliITSLoader.h>
#include <AliMultiplicity.h>
#include <AliITSMultReconstructor.h>

ClassImp(AliITSVertexer)

//////////////////////////////////////////////////////////////////////
// Base class for primary vertex reconstruction                     //
// AliESDVertexer is a class for full 3D primary vertex finding     //
// derived classes: AliITSVertexerIons AliITSvertexerPPZ            //
//                  AliITSVertexer3D                                //
//////////////////////////////////////////////////////////////////////

//______________________________________________________________________
AliITSVertexer::AliITSVertexer():AliVertexer() {
  // Default Constructor
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
  fDebug = 0;
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
void AliITSVertexer::FindMultiplicity(Int_t evnumber){
  // Invokes AliITSMultReconstructor to determine the
  // charged multiplicity in the pixel layers
  if(fMult){delete fMult; fMult = 0;}
  Bool_t success=kTRUE;
  if(!fCurrentVertex)success=kFALSE;
  if(fCurrentVertex && fCurrentVertex->GetNContributors()<1)success=kFALSE;
  if(!success){
    AliWarning("Tracklets multiplicity not determined because the primary vertex was not found");
    return;
  }
  AliITSMultReconstructor* multReco = new AliITSMultReconstructor();
  AliRunLoader *rl =AliRunLoader::GetRunLoader();
  AliITSLoader* itsLoader = (AliITSLoader*)rl->GetLoader("ITSLoader");
  multReco->SetGeometry(itsLoader->GetITSgeom());
  itsLoader->LoadRecPoints();
  rl->GetEvent(evnumber);
  TTree* itsClusterTree = itsLoader->TreeR();
  if (!itsClusterTree) {
    AliError(" Can't get the ITS cluster tree !\n");
    return;
  }
  Double_t vtx[3];
  fCurrentVertex->GetXYZ(vtx);
  Float_t vtxf[3];
  for(Int_t i=0;i<3;i++)vtxf[i]=vtx[i];
  multReco->SetHistOn(kFALSE);
  multReco->Reconstruct(itsClusterTree,vtxf,vtxf);
  Int_t notracks=multReco->GetNTracklets();
  Float_t *tht = new Float_t [notracks];
  Float_t *phi = new Float_t [notracks];
  Float_t *dphi = new Float_t [notracks];
  for(Int_t i=0;i<multReco->GetNTracklets();i++){
    tht[i] = multReco->GetTracklet(i)[0];
    phi[i] =  multReco->GetTracklet(i)[1];
    dphi[i] = multReco->GetTracklet(i)[2];
  }
  Int_t nosingleclus=multReco->GetNSingleClusters();
  Float_t *ths = new Float_t [nosingleclus];
  Float_t *phs = new Float_t [nosingleclus];
  for(Int_t i=0;i<nosingleclus;i++){
    ths[i] = multReco->GetCluster(i)[0];
    phs[i] =  multReco->GetCluster(i)[1];
  }
  fMult = new AliMultiplicity(notracks,tht,phi,dphi,nosingleclus,ths,phs);
  delete [] tht;
  delete [] phi;
  delete [] dphi;
  delete [] ths;
  delete [] phs;
  itsLoader->UnloadRecPoints();
  delete multReco;
  return;
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
