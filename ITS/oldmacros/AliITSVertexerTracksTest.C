#if !defined(__CINT__) || defined(__MAKECINT__)
//-- --- standard headers------------- 
#include <Riostream.h>
//--------Root headers ---------------
#include <TSystem.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TObject.h>
#include <TTree.h>
//----- AliRoot headers ---------------
#include "alles.h"
#include "AliRun.h"
#include "AliMagF.h"
#include "AliKalmanTrack.h"
#include "AliITSVertex.h"
#include "AliITSVertexer.h"
#include "AliITSVertexerTracks.h"
//-------------------------------------
#endif
void AliITSVertexerTracksTest(Int_t evFirst=0,Int_t evLast=0,
			      const Char_t *galiceName="galice.root",
			      const Char_t *trksName="AliITStracksV2.root",
			      const Char_t *vtxName="AliITSVertices.root") {
  /*******************************************************************
   *                                                                 *
   * Test macro for vertexing in pp using tracks.                    *
   * Input file must contain trees with AliITStrackV2 objects.       *
   * Output file can be the same file with the tracks                *
   * or another file.                                                *
   * If the file galice.root is available, B is taken from there,    *
   * otherwise is can be set here "by hand".                         *
   *                                                                 *
   * Origin: A.Dainese, Padova  andrea.dainese@pd.infn.it            *
   *******************************************************************/     

  // Look for field value in galice.root
  Double_t field = 0.4;
  if(!gSystem->AccessPathName(galiceName,kFileExists)) {
    TFile *galice = new TFile(galiceName);
    gAlice = (AliRun*)galice->Get("gAlice");
    AliMagF *fiel = (AliMagF*)gAlice->Field();
    field=(Double_t)fiel->SolenoidField()/10.;
    AliKalmanTrack::SetConvConst(100/0.299792458/field);
    printf(" B = %3.1f read from gAlice and set\n",field);
    delete gAlice;
    gAlice = 0;
    galice->Close();
    delete galice;
  } else {
    printf(" File galice.root not found: default 0.4 T being used!\n");
  }

  // Open input and output files
  TFile *inFile = TFile::Open(trksName);
  TFile *outFile = TFile::Open(vtxName,"recreate");

  // Create vertexer
  AliITSVertexerTracks *vertexer = 
    new AliITSVertexerTracks(inFile,outFile,field);
  vertexer->SetFirstEvent(evFirst);
  vertexer->SetLastEvent(evLast);
  vertexer->SetDebug(0);
  vertexer->SetUseThrustFrame(0);
  vertexer->PrintStatus();
  // Find vertices
  vertexer->FindVertices();


  delete vertexer;



  inFile->Close();
  outFile->Close();
  delete inFile;
  delete outFile;

  return;
}
