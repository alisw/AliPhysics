#include <Riostream.h>
#include <AliITSVertex.h>
#include <AliITSVertexer.h>
#include <AliRunLoader.h>
#include <AliITSLoader.h>

ClassImp(AliITSVertexer)

//////////////////////////////////////////////////////////////////////
// Base class for primary vertex reconstruction                     //
// AliITSVertex is a class for full 3D primary vertex finding       //
//////////////////////////////////////////////////////////////////////

//______________________________________________________________________
AliITSVertexer::AliITSVertexer() {
  // Default Constructor

    fCurrentVertex  = 0;
    SetDebug();
    SetFirstEvent(0);
    SetLastEvent(0);
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
  AliITSLoader* ITSloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  if(filename.Data()!="default")ITSloader->SetVerticesFileName(filename);
  ITSloader->LoadVertices("recreate");
  ITSloader->LoadRecPoints();
  Int_t lst;
  if(rl->TreeE()){
    lst = static_cast<Int_t>(rl->TreeE()->GetEntries());
    SetLastEvent(lst-1);
  }
}

//______________________________________________________________________
AliITSVertexer::~AliITSVertexer() {
  // Default Destructor
  // The objects poited by the following pointers are not owned
  // by this class and are not deleted

    fCurrentVertex  = 0;
}

//______________________________________________________________________
void AliITSVertexer::WriteCurrentVertex(){
  // Write the current AliVertex object to file fOutFile
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliITSLoader* ITSloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  fCurrentVertex->SetName("Vertex");
  //  const char * name = fCurrentVertex->GetName();
  //  ITSloader->SetVerticesContName(name);
  Int_t rc = ITSloader->PostVertex(fCurrentVertex);
  rc = ITSloader->WriteVertices();
}
