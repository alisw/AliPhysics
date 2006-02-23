#if !defined(__CINT__) || defined(__MAKECINT__)
//-- --- standard headers------------- 
#include <Riostream.h>
//--------Root headers ---------------
#include <TSystem.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TObject.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TKey.h>
//----- AliRoot headers ---------------
#include "alles.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliMagF.h"
#include "AliESDVertex.h"
#include "AliITSVertexer.h"
#include "AliITSVertexerTracks.h"
#include "AliTracker.h"
//-------------------------------------
#endif
void AliITSVertexerTracksTest(Int_t evFirst=0,Int_t evLast=0,Bool_t esd=kTRUE,
			      const Char_t *galiceName="galice.root",
			      const Char_t *inName="AliESDs.root",
                              const Char_t *outName="AliESDs.root") {
  /*******************************************************************
   *                                                                 *
   * Test macro for vertexing in pp using tracks.                    *
   * If(esd) {                                                       *
   *   Input file must contain ESD events with tracks reco in ITS    *
   *   if output file = input file, vertices will be stored in the   *
   *   ESD events and these rewritten to the same file; otherwise    *
   *   they will be written to a new file                            *
   * } else {                                                        *
   *   input file must contain a tree with AliITStracksV2 and        *
   *   can be written to the same file (output file = input file)    *
   *   or to a new one                                               *
   * }                                                               *
   * If file galice.root is present field will be read from it       *
   * otherwise it can be set here "by hand".                         *
   *                                                                 *
   * Origin: A.Dainese, Padova  andrea.dainese@pd.infn.it            *
   *******************************************************************/     

  TStopwatch timer;
  timer.Start();

  // Look for field value in galice.root
  Double_t field = 0.4;
  if(!gSystem->AccessPathName(galiceName,kFileExists)) {
    AliRunLoader *rl = AliRunLoader::Open(galiceName);
    rl->LoadgAlice();
    AliMagF *fiel = (AliMagF*)gAlice->Field();
    field=(Double_t)fiel->SolenoidField()/10.;
    printf(" B = %3.1f read from gAlice and set\n",field);

    AliTracker::SetFieldMap(gAlice->Field(),1); // 1 means uniform magnetic field
  } else {
    printf(" File galice.root not found: default 0.4 T being used!\n");
  }

  // Open input and output files
  TFile *inFile  = TFile::Open(inName);
  TFile *outFile = TFile::Open(outName,"update");

  // Create vertexer
  AliITSVertexerTracks *vertexer = 
    new AliITSVertexerTracks(inFile,outFile,evFirst,evLast);
  // Find vertices
  if(esd) {
    vertexer->FindVerticesESD();
  } else {
    vertexer->FindVertices();
  }

  timer.Stop(); timer.Print();

  delete vertexer;

  inFile->Close();
  outFile->Close();
  delete inFile;
  delete outFile;

  return;
}
//----------------------------------------------------------------------------
void VertexRecoInESDChain() {

  TStopwatch timer;
  timer.Start();

  // Look for field value in galice.root
  Double_t field = 0.4;
  if(!gSystem->AccessPathName("galice.root",kFileExists)) {
    AliRunLoader *rl = AliRunLoader::Open("galice.root");
    rl->LoadgAlice();
    AliMagF *fiel = (AliMagF*)gAlice->Field();
    field=(Double_t)fiel->SolenoidField()/10.;
    printf(" B = %3.1f read from gAlice and set\n",field);

    AliTracker::SetFieldMap(gAlice->Field(),1); // 1 means uniform magnetic field


  } else {
    printf(" File galice.root not found: default 0.4 T being used!\n");
  }

  // open file with ESD events (this would not be necessary in the chain...)
  if(gSystem->AccessPathName("AliESDs.root",kFileExists))
    { printf(" ESD file not found!\n"); return; }
  TFile *esdFile = new TFile("AliESDs.root","update");
 
  AliITSVertexerTracks *vertexer = new AliITSVertexerTracks();
  vertexer->SetField(field);

  Int_t n=0;
  TKey *key=0;
  TIter next(esdFile->GetListOfKeys());

  //******* The loop over events
  while ((key=(TKey*)next())!=0) {
    printf("--- Processing event number : %d",n++);

    AliESD *event=(AliESD*)key->ReadObj();

    // find vertex and store it in the ESD (only pos and cov matrix)
    vertexer->FindVertexForCurrentEvent(event);

    // write event on file
    event->Write();

  }

  timer.Stop(); timer.Print();

  delete vertexer;
  esdFile->Close();
  delete esdFile;

  return;
}







