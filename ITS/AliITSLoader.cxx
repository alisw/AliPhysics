#include "AliITSLoader.h"




//////////////////////////////////////////////////////////////////////////////////////////
// Loader for ITS
// it manages the I/O for:
// raw clusters, primary vertices
// V0 and cascade
// and tracks propagated to the origin
//////////////////////////////////////////////////////////////////////////////////////////
const TString AliITSLoader::fgkDefaultRawClustersContainerName = "TreeC";
const TString AliITSLoader::fgkDefaultBackTracksContainerName = "TreeB";
const TString AliITSLoader::fgkDefaultVerticesContainerName = "Vertex";
const TString AliITSLoader::fgkDefaultV0ContainerName = "V0";
const TString AliITSLoader::fgkDefaultCascadeContainerName = "Cascade";

ClassImp(AliITSLoader)

/*****************************************************************************/ 

  AliITSLoader::AliITSLoader():AliLoader(){
  // Default constructor

}

/*****************************************************************************/ 
AliITSLoader::AliITSLoader(const Char_t *name,const Char_t *topfoldername): AliLoader(name,topfoldername){
  //ctor   
  AliDataLoader* rawClustersDataLoader = new AliDataLoader(fDetectorName + ".RawCl.root",fgkDefaultRawClustersContainerName,"Raw Clusters");
  fDataLoaders->Add(rawClustersDataLoader);
  rawClustersDataLoader->SetEventFolder(fEventFolder);
  rawClustersDataLoader->SetFolder(GetDetectorDataFolder());

  AliDataLoader* backTracksDataLoader = 
                 new AliDataLoader(fDetectorName + ".BackTracks.root",fgkDefaultBackTracksContainerName,"Back Propagated Tracks");
  fDataLoaders->Add(backTracksDataLoader);
  backTracksDataLoader->SetEventFolder(fEventFolder);
  backTracksDataLoader->SetFolder(GetDetectorDataFolder());

  AliDataLoader* vertexDataLoader = new AliDataLoader(fDetectorName + ".Vertex.root",fgkDefaultVerticesContainerName,"Primary Vertices","O");
  fDataLoaders->Add(vertexDataLoader);
  vertexDataLoader->SetEventFolder(fEventFolder);
  vertexDataLoader->SetFolder(GetDetectorDataFolder());

  AliDataLoader* v0DataLoader = new AliDataLoader(fDetectorName + ".V0s.root",fgkDefaultV0ContainerName,"V0 Vertices");
  fDataLoaders->Add(v0DataLoader);
  v0DataLoader->SetEventFolder(fEventFolder);
  v0DataLoader->SetFolder(GetDetectorDataFolder());
   
  AliDataLoader* cascadeDataLoader = new AliDataLoader(fDetectorName + ".Cascades.root",fgkDefaultCascadeContainerName,"Cascades");
  fDataLoaders->Add(cascadeDataLoader);
  cascadeDataLoader->SetEventFolder(fEventFolder);
  cascadeDataLoader->SetFolder(GetDetectorDataFolder());
   
}
/*****************************************************************************/ 

AliITSLoader::AliITSLoader(const Char_t *name,TFolder *topfolder): AliLoader(name,topfolder) {
  //ctor  
  AliDataLoader*  rawClustersDataLoader = new AliDataLoader(fDetectorName + ".RawCl.root",fgkDefaultRawClustersContainerName,"Raw Clusters"); 
  fDataLoaders->Add(rawClustersDataLoader);
  rawClustersDataLoader->SetEventFolder(fEventFolder);
  rawClustersDataLoader->SetFolder(GetDetectorDataFolder());

  AliDataLoader*  backTracksDataLoader = 
                  new AliDataLoader(fDetectorName + ".BackTracks.root",fgkDefaultBackTracksContainerName,"Back Propagated Tracks");
  fDataLoaders->Add(backTracksDataLoader);
  backTracksDataLoader->SetEventFolder(fEventFolder);
  backTracksDataLoader->SetFolder(GetDetectorDataFolder());

  AliDataLoader* vertexDataLoader = new AliDataLoader(fDetectorName + ".Vertex.root",fgkDefaultVerticesContainerName,"Primary Vertices","O");
  fDataLoaders->Add(vertexDataLoader);
  vertexDataLoader->SetEventFolder(fEventFolder);
  vertexDataLoader->SetFolder(GetDetectorDataFolder());

  AliDataLoader* v0DataLoader = new AliDataLoader(fDetectorName + ".V0.root",fgkDefaultV0ContainerName,"V0 Vertices");
  fDataLoaders->Add(v0DataLoader);
  v0DataLoader->SetEventFolder(fEventFolder);
  v0DataLoader->SetFolder(GetDetectorDataFolder());
   
  AliDataLoader* cascadeDataLoader = new AliDataLoader(fDetectorName + ".Cascade.root",fgkDefaultCascadeContainerName,"Cascades");
  fDataLoaders->Add(cascadeDataLoader);
  cascadeDataLoader->SetEventFolder(fEventFolder);
  cascadeDataLoader->SetFolder(GetDetectorDataFolder());
   
}
/*****************************************************************************/ 
AliITSLoader::~AliITSLoader()
{
 //destructor
  AliDataLoader* dl = 0;
  UnloadRawClusters();
  dl = GetRawClLoader();
  fDataLoaders->Remove(dl);

  UnloadBackTracks();
  dl = GetBackTracksDataLoader();
  fDataLoaders->Remove(dl);

  UnloadVertices();
  dl = GetVertexDataLoader();
  fDataLoaders->Remove(dl);

  UnloadV0s();
  dl = GetV0DataLoader();
  fDataLoaders->Remove(dl);

  UnloadCascades();
  dl = GetCascadeDataLoader();
  fDataLoaders->Remove(dl);

}

void AliITSLoader::MakeTree(Option_t *opt){
  // invokes AliLoader::MakeTree + specific ITS tree(s)
  // Valid options: H,S,D,R,T and C (C=raw clusters)
  AliLoader::MakeTree(opt);
  const char *oC = strstr(opt,"C");
  if (oC) MakeRawClustersContainer();

  const char *oB = strstr(opt,"B");
  if (oB) MakeBackTracksContainer();
  
  const char *oV0 = strstr(opt,"V0");
  if (oV0) MakeV0Container();
  
  const char *oX = strstr(opt,"X");
  if (oX) MakeCascadeContainer();
  
}
