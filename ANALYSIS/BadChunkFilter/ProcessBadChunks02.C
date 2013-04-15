/********************************************************************
 
 Bad Chunks Checking code, 15th April 2013
 
 --- This version is a bit more "automatized" in that you give it a
 dataset string as a parameter and it spits out appropriately named 
 text files. Some small customization for each usage case (output 
 directory, in the first lines of the function below) may still 
 be needed when being used in general. 
 
 --- Also note: the code expects to find an "output" directory to 
 store text files with summaries. If it does not exist, please create
 it or also chang the path!
 
********************************************************************/
 

#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace std;

int ProcessBadChunks02(TString lDataset ){
  
  //Path to output files: change as needed
  TString lFileName = "/Volumes/MyPassport/work/download/badchunk/";
  lFileName.Append(lDataset.Data());
  lFileName.Append(".root");
  
  //Print out "I'm alive"
  cout<<"----------------------------------------------------"<<endl;
  cout<<" Bad Chunk Analysis Macro"<<endl;
  cout<<"----------------------------------------------------"<<endl;
  cout<<" Dataset identified as......: "<<lDataset<<endl;
  cout<<" Filename to open...........: "<<lFileName<<endl;
  
  //Open...
  // Open the file
  TFile* file = TFile::Open(lFileName, "READ");
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
  
  cout<<" Entries in Event Tree......: "<<ftree->GetEntries()<<endl;
  cout<<"----------------------------------------------------"<<endl;
  //return 0;
  Int_t llRunNumber = 0;
  TString * llFileName = 0x0;
  Int_t llNTracks = 0;
  Int_t llNGlobalTracks = 0;
  
  ftree->SetBranchAddress("fRunNumber"      ,&llRunNumber     );
  ftree->SetBranchAddress("fFileName"       , &llFileName   );
  ftree->SetBranchAddress("fNTracks"        ,&llNTracks       );
  ftree->SetBranchAddress("fNGlobalTracks"  ,&llNGlobalTracks );

  //Output ostreams: text files with chunk location, good and bad
  
  //Save by default to "output" directory
  TString lGoodName = "output/GoodChunks-";
  lGoodName.Append(lDataset.Data());
  lGoodName.Append(".txt");
  TString lBadName = "output/BadChunks-";
  lBadName.Append(lDataset.Data());
  lBadName.Append(".txt");

  TString lReport = "output/Datasetreport-";
  lReport.Append(lDataset.Data());
  lReport.Append(".txt");
  
  filebuf fbgood;
  fbgood.open (lGoodName,ios::out);
  ostream osgood(&fbgood);

  filebuf fbbad;
  fbbad.open (lBadName,ios::out);
  ostream osbad(&fbbad);

  filebuf fbreport;
  fbreport.open (lReport,ios::out);
  ostream osreport(&fbreport);
  
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
  cout<<"Corruption rate..........: "<<((double)lProcessedChunksBad)/((double)lProcessedChunks)<<endl;
  cout<<"Corruption rate, percent.: "<<100.*((double)lProcessedChunksBad)/((double)lProcessedChunks)<<" percent"<<endl;
  cout<<"----------------------------------------------------"<<endl;

  //Write report file 
  osreport<<"----------------------------------------------------"<<endl;
  osreport<<" "<<lDataset<<" Processed chunks, total..: "<<lProcessedChunks<<endl;
  osreport<<" "<<lDataset<<" Processed chunks, good...: "<<lProcessedChunksGood<<endl;
  osreport<<" "<<lDataset<<" Processed chunks, bad....: "<<lProcessedChunksBad<<endl;
  osreport<<" "<<lDataset<<" Corruption rate..........: "<<((double)lProcessedChunksBad)/((double)lProcessedChunks)<<endl;
  osreport<<" "<<lDataset<<" Corruption rate, percent.: "<<100.*((double)lProcessedChunksBad)/((double)lProcessedChunks)<<" percent"<<endl;
  osreport<<"----------------------------------------------------"<<endl;
  
  
  fbreport.close();
  cout<<endl;
  cout<<"---> Good chunks saved to \"goodguys.txt\""<<endl;
  cout<<"---> Bad chunks saved to \"badguys.txt\""<<endl;
  cout<<endl;
  cout<<"DONE!"<<endl;
  return 0;
   
}
