#include "AliITSLoader.h"
#include <AliRunLoader.h>

#include <TTree.h>
#include <TFile.h>

const TString AliITSLoader::fgkDefaultRawClustersContainerName = "TreeC";
const TString AliITSLoader::fgkDefaultBackTracksContainerName = "TreeB";
const TString AliITSLoader::fgkDefaultVerticesContainerName = "Vertex";

ClassImp(AliITSLoader)

/*****************************************************************************/ 
AliITSLoader::AliITSLoader(const Char_t *name,const Char_t *topfoldername):
 AliLoader(name,topfoldername),
 fRawClustersDataLoader(fDetectorName + ".RawCl.root",fgkDefaultRawClustersContainerName,"Raw Clusters"),
 fBackTracksDataLoader(fDetectorName + ".BackTracks.root",fgkDefaultBackTracksContainerName,"Back Propagated Tracks"),
  fVertexDataLoader(fDetectorName + ".Vertex.root",fgkDefaultVerticesContainerName,"Primary Vertices","O")
{
//ctor   
   fDataLoaders->Add(&fRawClustersDataLoader);
   fRawClustersDataLoader.SetEventFolder(fEventFolder);
   fRawClustersDataLoader.SetFolder(GetDetectorDataFolder());

   fDataLoaders->Add(&fBackTracksDataLoader);
   fBackTracksDataLoader.SetEventFolder(fEventFolder);
   fBackTracksDataLoader.SetFolder(GetDetectorDataFolder());

   fDataLoaders->Add(&fVertexDataLoader);
   fVertexDataLoader.SetEventFolder(fEventFolder);
   fVertexDataLoader.SetFolder(GetDetectorDataFolder());
   
}
/*****************************************************************************/ 

AliITSLoader::AliITSLoader(const Char_t *name,TFolder *topfolder):
 AliLoader(name,topfolder),
 fRawClustersDataLoader(fDetectorName + ".RawCl.root",fgkDefaultRawClustersContainerName,"Raw Clusters"),
 fBackTracksDataLoader(fDetectorName + ".BackTracks.root",fgkDefaultBackTracksContainerName,"Back Propagated Tracks"),
  fVertexDataLoader(fDetectorName + ".Vertex.root",fgkDefaultVerticesContainerName,"Primary Vertices","O")
{
//ctor   
   fDataLoaders->Add(&fRawClustersDataLoader);
   fRawClustersDataLoader.SetEventFolder(fEventFolder);
   fRawClustersDataLoader.SetFolder(GetDetectorDataFolder());

   fDataLoaders->Add(&fBackTracksDataLoader);
   fBackTracksDataLoader.SetEventFolder(fEventFolder);
   fBackTracksDataLoader.SetFolder(GetDetectorDataFolder());

   fDataLoaders->Add(&fVertexDataLoader);
   fVertexDataLoader.SetEventFolder(fEventFolder);
   fVertexDataLoader.SetFolder(GetDetectorDataFolder());

}
/*****************************************************************************/ 
AliITSLoader::~AliITSLoader()
{
 //destructor
  UnloadRawClusters();
  fDataLoaders->Remove(&fRawClustersDataLoader);

  UnloadBackTracks();
  fDataLoaders->Remove(&fBackTracksDataLoader);

  UnloadVertices();
  fDataLoaders->Remove(&fVertexDataLoader);

}

void AliITSLoader::MakeTree(Option_t *opt){
  // invokes AliLoader::MakeTree + specific ITS tree(s)
  // Valid options: H,S,D,R,T and C (C=raw clusters)
  AliLoader::MakeTree(opt);
  const char *oC = strstr(opt,"C");
  if (oC) MakeRawClustersContainer();

  const char *oB = strstr(opt,"B");
  if (oB) MakeBackTracksContainer();
  
}
