void AliITSDumpVertices(Int_t firstEv=0, Int_t noev=1,
       TString fileimp="AliITSVertices.root",
       TString objbasename="VertexTracks_"){
  // This is a simple example on how to access the vertex objects
  // The default object base name is VertexTracks_  for AliITSVertexerTracks
  //                                 Vertex_        for AliITSVertexerPPZ
  //                                 Vertex_        for AliITSVertexerIons
  Int_t evmax = firstEv+noev;
  TFile *file = new TFile(fileimp);
  for(Int_t i=firstEv; i<evmax; i++){
    TString name = objbasename;
    name += i;
    vert = (AliITSVertex*)file->Get(name.Data());
    if(vert){
      cout <<"===============================================\n";
      cout <<" Event n. "<<i<<endl;
      vert->PrintStatus();
    }
  }
}
