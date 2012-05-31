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
#include <AliHeader.h>
#include <AliGenEventHeader.h>
//-------------------------------------
#endif
void AliITSVertexerTracksTest2(Int_t evFirst=0,Int_t evLast=0,
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
  Int_t kDebug = 0;
  TFile *galice = 0;
  if(!gSystem->AccessPathName(galiceName,kFileExists)) {
    galice = new TFile(galiceName);
    gAlice = (AliRun*)galice->Get("gAlice");
    AliMagF *fiel = TGeoGlobalMagField::Instance()->GetField();
    field=(Double_t)fiel->SolenoidField()/10.;
    AliKalmanTrack::SetConvConst(100/0.299792458/field);
    printf(" B = %3.1f read from gAlice and set\n",field);
  } else {
    printf(" File galice.root not found: default 0.4 T being used!\n");
  }

  // Open input and output files
  TFile *inFile = TFile::Open(trksName);
  TFile *outFile = TFile::Open(vtxName,"recreate");

  // Create vertexer
  AliITSVertexerTracks *vertexer = 
    new AliITSVertexerTracks(inFile,outFile,field,0.,0.);
  vertexer->SetDebug(0);
  vertexer->SetUseThrustFrame(0);
  vertexer->PrintStatus();

  AliITSVertex *vert = 0;
  // Find vertices

  for(Int_t i=evFirst; i<=evLast; i++){
    if(i%100==0)cout<<"processing event "<<i<<endl;
    gAlice->GetEvent(i);
    // The true Z coord. is fetched for comparison
    AliHeader *header = gAlice->GetHeader();
    AliGenEventHeader* genEventHeader = header->GenEventHeader();
    TArrayF primaryVertex(3);
    genEventHeader->PrimaryVertex(primaryVertex);
    vert = vertexer->FindVertexForCurrentEvent(i);
    if(kDebug>0){
      // Prints the results
      cout <<"========================================================\n";
      cout << "Event number: "<<i<<")  Z Vertex:"<<endl;
      if(vert){
	cout<<"FOUND: "<<vert->GetZv()<<"; "<<vert->GetZRes()<<endl;
	cout <<" True Z position "<<primaryVertex[2]<<endl;
	cout<<", diff= "<<(primaryVertex[2]-vert->GetZv())*10000.<<endl;
      }
      else {
	cout<<"NOT FOUND "<<endl;
      }
    }
    if(vert){
      Double_t pos[3];
      for(Int_t kk=0;kk<3;kk++)pos[kk]=(Double_t)primaryVertex[kk];
      vert->SetTruePos(pos);
      vertexer->WriteCurrentVertex();
    }
  }

  delete vertexer;



  inFile->Close();
  outFile->Close();
  delete inFile;
  delete outFile;
  galice->Close();
  delete galice;
  return;
}
