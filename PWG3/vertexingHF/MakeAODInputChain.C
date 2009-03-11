TChain *MakeAODInputChain(const char* collectionfileAOD,
			  const char* collectionfileAODfriend) {
  //
  // Check one-to-one correspondence of the two collections
  // Create AOD chain with friend AODVertexingHF chain
  // Origin: A.Rossi, andrea.rossi@ts.infn.it
  //

  TAlienCollection *collectionAOD       = TAlienCollection::Open(collectionfileAOD);
  TAlienCollection *collectionAODfriend = TAlienCollection::Open(collectionfileAODfriend);
    

  TGridResult *tagResultAOD = collectionAOD->GetGridResult("",0,0);
  TGridResult *tagResultAODFriend = collectionAODfriend->GetGridResult("",0,0);
  //  tagResultAOD->Print();
  TMap *mappa;
  TObjString *filename;
  Int_t event,dir,pos;
  string str;
  Int_t lastAOD[1000],lastAODFriend[1000];
  TArrayI ***arrAOD=new TArrayI**[1000];
  TArrayI ***arrAODFriend=new TArrayI**[1000];

  for(Int_t j=0;j<1000;j++) {
    arrAOD[j]=new TArrayI*[2];
    arrAOD[j][0]=new TArrayI(100);
    arrAOD[j][0]->Reset(0);
    arrAOD[j][1]=new TArrayI(100);
    arrAOD[j][1]->Reset(-1);

    arrAODFriend[j]=new TArrayI*[2];
    arrAODFriend[j][0]=new TArrayI(100);
    arrAODFriend[j][0]->Reset(0);
    arrAODFriend[j][1]=new TArrayI(100);
    arrAODFriend[j][1]->Reset(-1);

    lastAOD[j]=0;
    lastAODFriend[j]=0;
  }

  for(Int_t j=0;j<tagResultAOD->GetEntries();j++) {  
    mappa=(TMap*)tagResultAOD->At(j);
    filename=(TObjString*)mappa->GetValue("turl");
    str=filename->GetString();
    pos=str.find_last_of("/");
    str.string::replace(pos,50,"");
    pos=str.find_last_of("/");
    str.string::replace(pos,1," ");
    sscanf(str.data(),"%*s %d",&event);
    pos=str.find_last_of("/");
    str.string::replace(pos,1," ");
    sscanf(str.data(),"%*s %d",&dir);
    arrAOD[event][0]->AddAt(dir,lastAOD[event]);
    arrAOD[event][1]->AddAt(j,lastAOD[event]);
    // printf("Adding AOD, event:%d,dir:%d,position:%d\n",event,dir,j);
    lastAOD[event]++;
  }
  
  for(Int_t j=0;j<tagResultAODFriend->GetEntries();j++) {  
    mappa=(TMap*)tagResultAODFriend->At(j);
    filename=(TObjString*)mappa->GetValue("turl");
    str=filename->GetString();
    pos=str.find_last_of("/");
    str.string::replace(pos,50,"");
    pos=str.find_last_of("/");
    str.string::replace(pos,1," ");
    sscanf(str.data(),"%*s %d",&event);
    pos=str.find_last_of("/");
    str.string::replace(pos,1," ");
    sscanf(str.data(),"%*s %d",&dir);
    arrAODFriend[event][0]->AddAt(dir,lastAODFriend[event]);
    arrAODFriend[event][1]->AddAt(j,lastAODFriend[event]);
    //printf("Adding AODFriend, event:%d,dir:%d,position:%d\n",event,dir,j);
    lastAODFriend[event]++;
  }
 
  TChain *chainAOD = new TChain("aodTree");
  TChain *chainAODfriend = new TChain("aodTree");
  
  for(Int_t ev=0;ev<1000;ev++) {
    for(Int_t j=0;j<lastAOD[ev];j++) {
      dir= arrAOD[ev][0]->At(j);
      for(Int_t k=0;k<lastAODFriend[ev];k++) {
	if(arrAODFriend[ev][0]->At(k)==dir){
	  chainAOD->Add((((TObjString*)((TMap*)tagResultAOD->At(arrAOD[ev][1]->At(j)))->GetValue("turl"))->GetString()).Data());
	  chainAODfriend->Add((((TObjString*)((TMap*)tagResultAODFriend->At(arrAODFriend[ev][1]->At(k)))->GetValue("turl"))->GetString()).Data());
	  printf("Events: %d, adding AOD at position %d,dir:%d posarray:%d \n and Friend at %d, dir:%d , posarray: %d \n \n",ev,arrAOD[ev][1]->At(j),arrAOD[ev][0]->At(j),j,arrAODFriend[ev][1]->At(k),arrAODFriend[ev][0]->At(k),k);
	  break;
	}
      }
    }
  }
  

  chainAOD->AddFriend(chainAODfriend);
 
  return chainAOD;
}
//----------------------------------------------------------------------------
TChain *MakeAODInputChain(const char* pathname="",
			  Int_t firstdir=1,Int_t lastdir=-1) {
  //
  // Create AOD chain with friend AODVertexingHF chain
  // Example path: "alien:///alice/cern.ch/user/r/rbala/analysis/out_lhcw/290001/"
  // Origin: A.Rossi, andrea.rossi@ts.infn.it
  //

 
  TChain *chainAOD = new TChain("aodTree");
  TChain *chainAODfriend = new TChain("aodTree");
 
  if(lastdir==-1) { // only one pair of files
    chainAOD->Add("AliAOD.root");
    chainAODfriend->Add("AliAOD.VertexingHF.root"); 
  } else {
    // set the path to the files (can be local or on alien)
    for(Int_t idir=firstdir; idir<=lastdir; idir++) {
      TString aodname=pathname;
      TString aodHFname=pathname;
      aodname+=idir;
      aodHFname+=idir;
      aodname.Append("/AliAOD.root");
      aodHFname.Append("/AliAOD.VertexingHF.root");
      chainAOD->Add(aodname.Data());
      chainAODfriend->Add(aodHFname.Data());
    }
  }

  chainAOD->AddFriend(chainAODfriend);
 
  return chainAOD;
}
