void AliITSDumpVertices(Int_t firstEv=0, Int_t noev=1,
                        TString fileimp="galice.root"){
  // This is a simple example on how to access the vertex objects
  Int_t evmax = firstEv+noev;
  AliRunLoader *rl = AliRunLoader::Open(fileimp.Data());
  AliITSLoader* ITSloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  ITSloader->LoadVertices();
  for(Int_t i=firstEv; i<evmax; i++){
    rl->GetEvent(i);
    AliESDVertex *vert = ITSloader->GetVertex();
    if(vert){
      cout <<"===============================================\n";
      cout <<" Event n. "<<i<<endl;
      vert->PrintStatus();
    }
  }
}
