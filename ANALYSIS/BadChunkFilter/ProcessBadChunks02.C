#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace std;

int ProcessBadChunks02(){
  // Open the file
  TFile* file = TFile::Open("BadChunkTree.root", "READ");
  if (!file || !file->IsOpen()) {
    cout<<"File not found!"<<endl;
    return 1;
  }
  // Get the tree
  TTree* ftree = (TTree*)file->FindObjectAny("fTree");
  if (!ftree) {
    cout<<"File doesn't contain fTree!"<<endl;
    return 2;
  }
  //Print out "I'm alive"
  cout<<"----------------------------------------------------"<<endl;
  cout<<" Bad Chunk Analysis Macro"<<endl;
  cout<<"----------------------------------------------------"<<endl;  
  cout<<" Entries in Event Tree......"<<ftree->GetEntries()<<endl;
  cout<<"----------------------------------------------------"<<endl;
  Int_t llRunNumber = 0;
  TString * llFileName = 0x0;
  Int_t llNTracks = 0;
  Int_t llNGlobalTracks = 0;
  
  ftree->SetBranchAddress("fRunNumber"      ,&llRunNumber     );
  ftree->SetBranchAddress("fFileName"       , &llFileName   );
  ftree->SetBranchAddress("fNTracks"        ,&llNTracks       );
  ftree->SetBranchAddress("fNGlobalTracks"  ,&llNGlobalTracks );

  //Output ostreams: text files with chunk location, good and bad
  
  filebuf fbgood;
  fbgood.open ("goodguys.txt",ios::out);
  ostream osgood(&fbgood);

  filebuf fbbad;
  fbbad.open ("badguys.txt",ios::out);
  ostream osbad(&fbbad);
  
  //Main Tree Loop
  Long_t lProcessedChunks = 0;
  Long_t lProcessedChunksGood = 0;
  Long_t lProcessedChunksBad = 0;
  Long_t lTracks = 0;
  Long_t lGlobalTracks = 0;
  TString *lChunkName = new TString();
  
  //First Chunk entry
  ftree->GetEntry(0);
  *lChunkName = *llFileName;
  
  cout<<"TEST  "<<lChunkName->Data() <<", track = "<<lTracks<<", globals = "<<lGlobalTracks<<" " << llFileName->Data()  << endl;
  
  for(Long_t iEvent = 1; iEvent<ftree->GetEntries(); iEvent++){
    ftree->GetEntry(iEvent);
    if( !llFileName->EqualTo(*lChunkName) ){
      lProcessedChunks++;
      //Change Chunk
      if( lTracks > 0 && lGlobalTracks ==0){
        //Candidate bad chunk found!
        cout<<"BAD CHUNK at "<<lChunkName->Data() <<", track = "<<lTracks<<", globals = "<<lGlobalTracks<< endl;
        osbad<<lChunkName->Data()<<endl;
        lProcessedChunksBad++;
      }else{
        //This looks OK...
        //cout<<"GOOD CHUNK at "<<lChunkName<<endl;
        osgood<<lChunkName->Data()<<endl;
        lProcessedChunksGood++;
      }
      //Get ready to loop over new chunk
      *lChunkName = *llFileName;
      lTracks = llNTracks;
      lGlobalTracks = llNGlobalTracks;
      if(lProcessedChunks%5000==0) cout<<"---> Processed "<<lProcessedChunks<<"..."<<endl;
    }else{
      lTracks        += llNTracks;
      lGlobalTracks  += llNGlobalTracks;
    }
  }
  //CLOSE the processing: one extra pass...
  lProcessedChunks++;
  if( lTracks > 0 && lGlobalTracks ==0){
    //Candidate bad chunk found!
    cout<<"BAD CHUNK at "<<lChunkName->Data()<<endl;
    osbad<<lChunkName->Data() <<endl;
    lProcessedChunksBad++;
  }else{
    //This looks OK...
    //cout<<"GOOD CHUNK at "<<lChunkName<<endl;
    osgood<<lChunkName->Data()<<endl;
    lProcessedChunksGood++;
  }
  fbgood.close();
  fbbad.close();
  cout<<"----------------------------------------------------"<<endl;
  cout<<"Processed chunks, total..: "<<lProcessedChunks<<endl;
  cout<<"Processed chunks, good...: "<<lProcessedChunksGood<<endl;
  cout<<"Processed chunks, bad....: "<<lProcessedChunksBad<<endl;
  cout<<"----------------------------------------------------"<<endl;
  cout<<endl;
  cout<<"---> Good chunks saved to \"goodguys.txt\""<<endl;
  cout<<"---> Bad chunks saved to \"badguys.txt\""<<endl;
  cout<<endl;
  cout<<"DONE!"<<endl;
  return 0;
   
}
