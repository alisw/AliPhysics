#include <AliITSVertex.h>
#include <AliITSVertexer.h>

ClassImp(AliITSVertexer)

//////////////////////////////////////////////////////////////////////
// Base class for primary vertex reconstruction                     //
// AliITSVertex is a class for full 3D primary vertex finding       //
//////////////////////////////////////////////////////////////////////

//______________________________________________________________________
AliITSVertexer::AliITSVertexer() {
  // Default Constructor

    fCurrentVertex  = 0;
    fInFile         = 0;
    fOutFile        = 0;
    SetDebug();
    SetFirstEvent(0);
    SetLastEvent(0);
}

AliITSVertexer::AliITSVertexer(TFile *infile, TFile *outfile) {
  // Standard constructor
  fCurrentVertex  = 0;   
  SetInputFile(infile);
  SetOutputFile(outfile);
  SetDebug();
  SetFirstEvent(0);
  SetLastEvent(0);
  if(gAlice){
    Int_t lst;
    if(gAlice->TreeE()){
      lst = static_cast<Int_t>(gAlice->TreeE()->GetEntries());
      SetLastEvent(lst);
    }
  }
}

//______________________________________________________________________
AliITSVertexer::~AliITSVertexer() {
  // Default Destructor
  // The objects poited by the following pointers are not owned
  // by this class and are not deleted

    fCurrentVertex  = 0;
    fInFile         = 0;
    fOutFile        = 0;
}

//______________________________________________________________________
void AliITSVertexer::WriteCurrentVertex(){
  // Write the current AliVertex object to file fOutFile
  if(!fOutFile){
    Error("WriteCurrentEvent","The output file is not defined");
    return;
  }
  TDirectory *curdir = gDirectory;
  fOutFile->cd();
  fCurrentVertex->Write();
  curdir->cd();
  fCurrentVertex = 0;
}
