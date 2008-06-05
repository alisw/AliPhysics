#include <TFile.h>
#include "AliGeomManager.h"
#include "AliHeader.h"
#include "AliITSDetTypeRec.h"
#include "AliITSInitGeometry.h"
#include "AliITSMeanVertexer.h"
#include "AliITSLoader.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include "AliRunLoader.h"
#include "AliITSVertexer3D.h"
#include "AliESDVertex.h"
#include "AliMeanVertex.h"
#include "AliMultiplicity.h"

const TString AliITSMeanVertexer::fgkMVFileNameDefault = "MeanVertex.root";


ClassImp(AliITSMeanVertexer)

///////////////////////////////////////////////////////////////////////
//                                                                   //
// Class to compute vertex position using SPD local reconstruction   //
// An average vertex position using all the events                   //
// is built and saved                                                //
// Individual vertices can be optionally stored                      //
// Origin: M.Masera  (masera@to.infn.it)                             //
// Usage:
// AliITSMeanVertexer mv("RawDataFileName");
// mv.SetGeometryFileName("FileWithGeometry.root"); // def. geometry.root
// ....  optional setters ....
// mv.Reconstruct();  // SPD local reconstruction
// mv.DoVertices(); 
// Resulting AliMeanVertex object is stored in a file named fMVFileName
///////////////////////////////////////////////////////////////////////

/* $Id$ */

//______________________________________________________________________
AliITSMeanVertexer::AliITSMeanVertexer():TObject(),
fLoaderFileName(),
fGeometryFileName(),
fMVFileName(),
fRawReader(),
fRunLoader(),
fNoEventsContr(0),
fTotTracklets(0.),
fAverTracklets(0.),
fSigmaOnAverTracks(0.), 
fFilterOnContributors(0),
fFilterOnTracklets(0),
fWriteVertices(kFALSE)
{
  // Default Constructor
  for(Int_t i=0;i<3;i++){
    fWeighPos[i] = 0.;
    fWeighSig[i] = 0.;
    fAverPos[i] = 0.;
    for(Int_t j=0; j<3;j++)fAverPosSq[i][j] = 0.;
  }
  
}

//______________________________________________________________________
AliITSMeanVertexer::AliITSMeanVertexer(TString &filename):TObject(),
fLoaderFileName(),
fGeometryFileName(),
fMVFileName(fgkMVFileNameDefault),
fRawReader(),
fRunLoader(),
fNoEventsContr(0),
fTotTracklets(0.),
fAverTracklets(0.),
fSigmaOnAverTracks(0.), 
fFilterOnContributors(0),
fFilterOnTracklets(0),
fWriteVertices(kTRUE)
{
  // Standard constructor

  for(Int_t i=0;i<3;i++){
    fWeighPos[i] = 0.;
    fWeighSig[i] = 0.;
    fAverPos[i] = 0.;
    for(Int_t j=0; j<3;j++)fAverPosSq[i][j] = 0.;
  }
  SetLoaderFileName();
  SetGeometryFileName();

  Init(filename);
}
//______________________________________________________________________
AliITSMeanVertexer::AliITSMeanVertexer(TString &filename, 
                                       TString &loaderfilename, 
				       TString &geometryfilename):TObject(),
fLoaderFileName(),
fGeometryFileName(),
fMVFileName(fgkMVFileNameDefault),
fRawReader(),
fRunLoader(),
fNoEventsContr(0), 
fTotTracklets(0.),
fAverTracklets(0.),
fSigmaOnAverTracks(0.), 
fFilterOnContributors(0),
fFilterOnTracklets(0),
fWriteVertices(kTRUE)
{
  // Standard constructor with explicit geometry file name assignment
  for(Int_t i=0;i<3;i++){
    fWeighPos[i] = 0.;
    fWeighSig[i] = 0.;
    fAverPos[i] = 0.;
    for(Int_t j=0; j<3;j++)fAverPosSq[i][j] = 0.;
  }
  SetLoaderFileName(loaderfilename);
  SetGeometryFileName(geometryfilename);
  Init(filename);
}

//______________________________________________________________________
void AliITSMeanVertexer::Init(TString &filename){
  // Initialization part common to different constructors
  if(filename.IsNull()){
    AliFatal("Please, provide a valid file name for raw data file\n");
  }
  // if file name ends with root a raw reader ROOT is assumed
  if(filename.EndsWith(".root")){
    fRawReader = new AliRawReaderRoot(filename);
  }
  else {  // DATE raw reader is assumed
    fRawReader = new AliRawReaderDate(filename);
    fRawReader->SelectEvents(7);
  }
  fRunLoader = AliRunLoader::Open(fLoaderFileName.Data(),AliConfig::GetDefaultEventFolderName(),"recreate");
  fRunLoader->MakeTree("E");
  Int_t iEvent = 0;
  while (fRawReader->NextEvent()) {
    fRunLoader->SetEventNumber(iEvent);
    fRunLoader->GetHeader()->Reset(fRawReader->GetRunNumber(), 
				   iEvent, iEvent);
    fRunLoader->MakeTree("H");
    fRunLoader->TreeE()->Fill();
    iEvent++;
  }
  fRawReader->RewindEvents();
  fRunLoader->SetNumberOfEventsPerFile(iEvent);
  fRunLoader->WriteHeader("OVERWRITE");
 Int_t retval = AliConfig::Instance()->AddDetector(fRunLoader->GetEventFolder(),"ITS","ITS");
 if(retval != 0)AliFatal("Not able to add ITS detector");
  AliITSLoader *loader = new AliITSLoader("ITS",fRunLoader->GetEventFolder()->GetName());
  fRunLoader->AddLoader(loader);
  fRunLoader->CdGAFile();
  fRunLoader->Write(0, TObject::kOverwrite);
  // Initialize geometry
 
  AliGeomManager::LoadGeometry(fGeometryFileName.Data());

  AliITSInitGeometry initgeom;
  AliITSgeom *geom = initgeom.CreateAliITSgeom();
  printf("Geometry name: %s \n",(initgeom.GetGeometryName()).Data());
  loader->SetITSgeom(geom);
  // Initialize filter values to their defaults
  SetFilterOnContributors();
  SetFilterOnTracklets();
}

//______________________________________________________________________
AliITSMeanVertexer::AliITSMeanVertexer(const AliITSMeanVertexer &vtxr) : TObject(vtxr),
fLoaderFileName(vtxr.fLoaderFileName),
fGeometryFileName(vtxr.fGeometryFileName),
fMVFileName(vtxr.fMVFileName),
fRawReader(vtxr.fRawReader),
fRunLoader(vtxr.fRunLoader),
fNoEventsContr(vtxr.fNoEventsContr),
fTotTracklets(vtxr.fTotTracklets),
fAverTracklets(vtxr.fAverTracklets),
fSigmaOnAverTracks(vtxr.fSigmaOnAverTracks),  
fFilterOnContributors(vtxr.fFilterOnContributors),
fFilterOnTracklets(vtxr.fFilterOnTracklets),
fWriteVertices(vtxr.fWriteVertices)
{
  // Copy constructor
  // Copies are not allowed. The method is protected to avoid misuse.
  AliFatal("Copy constructor not allowed\n");

}

//______________________________________________________________________
AliITSMeanVertexer& AliITSMeanVertexer::operator=(const AliITSMeanVertexer&  /* vtxr */ ){
  // Assignment operator
  AliError("Assignment operator not allowed\n");
  return *this;
}

//______________________________________________________________________
AliITSMeanVertexer::~AliITSMeanVertexer() {
  // Destructor
  delete fRawReader;
  delete fRunLoader;

}

//______________________________________________________________________
void AliITSMeanVertexer::Reconstruct(){
  // Performs SPD local reconstruction
  AliITSLoader* loader = static_cast<AliITSLoader*>(fRunLoader->GetLoader("ITSLoader"));
  if (!loader) {
    AliFatal("ITS loader not found");
    return;
  }
  AliITSDetTypeRec* rec = new AliITSDetTypeRec();
  rec->SetITSgeom(loader->GetITSgeom());
  rec->SetDefaults();

  rec->SetDefaultClusterFindersV2(kTRUE);

  Int_t nEvents = fRunLoader->GetNumberOfEvents();
  fRawReader->RewindEvents();
  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    fRawReader->NextEvent();
    fRunLoader->GetEvent(iEvent);
    AliDebug(1,Form(">>>>>>>   Processing event number: %d",iEvent));
    loader->LoadRecPoints("update");
    loader->CleanRecPoints();
    loader->MakeRecPointsContainer();
    TTree *tR = loader->TreeR();
    if(!tR){
      AliFatal("Tree R pointer not found - Abort \n");
      break;
    }
    rec->DigitsToRecPoints(fRawReader,tR,"SPD");
    rec->ResetRecPoints();
    rec->ResetClusters();    
    loader->WriteRecPoints("OVERWRITE");
    loader->UnloadRecPoints();
  }

}

//______________________________________________________________________
void AliITSMeanVertexer::DoVertices(){
  // Loop on all events and compute 3D vertices
  AliITSLoader* loader = static_cast<AliITSLoader*>(fRunLoader->GetLoader("ITSLoader"));
  AliITSVertexer3D *dovert = new AliITSVertexer3D();
  AliESDVertex *vert = 0;
  Int_t nevents = fRunLoader->TreeE()->GetEntries();
  for(Int_t i=0; i<nevents; i++){
    fRunLoader->GetEvent(i);
    TTree* cltree = loader->TreeR();
    vert = dovert->FindVertexForCurrentEvent(cltree);
    AliMultiplicity *mult = dovert->GetMultiplicity();
    if(Filter(vert,mult)){
      AddToMean(vert); 
      if(fWriteVertices){
	loader->PostVertex(vert);
	loader->WriteVertices();
      }
    }
    else {
      if(vert)delete vert;
    }
  }
  if(ComputeMean()){
    Double_t cov[6];
    cov[0] =  fAverPosSq[0][0];  // variance x
    cov[1] =  fAverPosSq[0][1];  // cov xy
    cov[2] =  fAverPosSq[1][1];  // variance y
    cov[3] =  fAverPosSq[0][2];  // cov xz
    cov[4] =  fAverPosSq[1][2];  // cov yz
    cov[5] =  fAverPosSq[2][2];  // variance z
    AliMeanVertex mv(fWeighPos,fWeighSig,cov,nevents,fTotTracklets,fAverTracklets,fSigmaOnAverTracks);
    mv.SetTitle("Mean Vertex");
    mv.SetName("Meanvertex");
    AliDebug(1,Form("Contrib av. trk = %10.2f ",mv.GetAverageNumbOfTracklets()));
    AliDebug(1,Form("Sigma %10.4f ",mv.GetSigmaOnAvNumbOfTracks()));
    TFile fmv(fMVFileName.Data(),"recreate");
    mv.Write();
    fmv.Close();
  }
  else {
    AliError(Form("Evaluation of mean vertex not possible. Number of used events = %d",fNoEventsContr));
  }
  delete dovert;
}

//______________________________________________________________________
Bool_t AliITSMeanVertexer::Filter(AliESDVertex *vert,AliMultiplicity *mult){
  // Apply selection criteria to events
  Bool_t status = kFALSE;
  if(!vert || !mult)return status;
  Int_t ncontr = vert->GetNContributors();
  Int_t ntracklets = mult->GetNumberOfTracklets();
  AliDebug(1,Form("Number of contributors = %d",ncontr));
  AliDebug(1,Form("Number of tracklets = %d",ntracklets));
  if(ncontr>fFilterOnContributors && ntracklets > fFilterOnTracklets) status = kTRUE;
  fAverTracklets += ntracklets;
  fSigmaOnAverTracks += ntracklets*ntracklets;
  return status;
}

//______________________________________________________________________
void AliITSMeanVertexer::AddToMean(AliESDVertex *vert){
  // update mean vertex
  Double_t currentPos[3],currentSigma[3];
  vert->GetXYZ(currentPos);
  vert->GetSigmaXYZ(currentSigma);
  Bool_t goon = kTRUE;
  for(Int_t i=0;i<3;i++)if(currentSigma[i] == 0.)goon = kFALSE;
  if(!goon)return;
  for(Int_t i=0;i<3;i++){
    fWeighPos[i]+=currentPos[i]/currentSigma[i]/currentSigma[i];
    fWeighSig[i]+=1./currentSigma[i]/currentSigma[i];
    fAverPos[i]+=currentPos[i];
  }
  for(Int_t i=0;i<3;i++){
    for(Int_t j=i;j<3;j++){
      fAverPosSq[i][j] += currentPos[i] * currentPos[j];
    }
  }
  fNoEventsContr++;
}

//______________________________________________________________________
Bool_t AliITSMeanVertexer::ComputeMean(){
  // compute mean vertex
  if(fNoEventsContr < 2) return kFALSE;
  Double_t nevents = fNoEventsContr;
  for(Int_t i=0;i<3;i++){
    fWeighPos[i] /= fWeighSig[i]; 
    fWeighSig[i] = 1./TMath::Sqrt(fWeighSig[i]);
    fAverPos[i] /= nevents;
  } 
  for(Int_t i=0;i<3;i++){
    for(Int_t j=i;j<3;j++){
      fAverPosSq[i][j] /= (nevents -1.);
      fAverPosSq[i][j] -= nevents/(nevents -1.)*fAverPos[i]*fAverPos[j];
    }
  } 
  fTotTracklets = fAverTracklets;  //  total number of tracklets used 
  fAverTracklets /= nevents;
  fSigmaOnAverTracks /= (nevents - 1);
  Double_t tmp = nevents/(nevents -1.)*fAverTracklets*fAverTracklets;
  fSigmaOnAverTracks -= tmp;
  fSigmaOnAverTracks = TMath::Sqrt(fSigmaOnAverTracks);
  return kTRUE;
}
