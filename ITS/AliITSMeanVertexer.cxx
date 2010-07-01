#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include "AliGeomManager.h"
#include "AliITSDetTypeRec.h"
#include "AliITSInitGeometry.h"
#include "AliITSMeanVertexer.h"
#include "AliITSRecPointContainer.h"
#include "AliITSLoader.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include "AliRunLoader.h"
#include "AliITSVertexer3D.h"
#include "AliITSVertexer3DTapan.h"
#include "AliESDVertex.h"
#include "AliMeanVertex.h"
#include "AliMultiplicity.h"

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
AliITSMeanVertexer::AliITSMeanVertexer(Bool_t mode):TObject(),
fDetTypeRec(NULL),
fVertexXY(NULL),
fVertexZ(NULL),
fNoEventsContr(0),
fTotTracklets(0.),
fAverTracklets(0.),
fTotTrackletsSq(0.),
fSigmaOnAverTracks(0.), 
fFilterOnContributors(0),
fFilterOnTracklets(0),
fMode(mode),
fVertexer(NULL)
{
  // Default Constructor
  for(Int_t i=0;i<3;i++){
    fWeighPosSum[i] = 0.;
    fWeighSigSum[i] = 0.;
    fAverPosSum[i] = 0.;
    fWeighPos[i] = 0.;
    fWeighSig[i] = 0.;
    fAverPos[i] = 0.;
    for(Int_t j=0; j<3;j++)fAverPosSq[i][j] = 0.;
    for(Int_t j=0; j<3;j++)fAverPosSqSum[i][j] = 0.;
  }

  // Histograms initialization
  const Float_t xLimit = 5.0, yLimit = 5.0, zLimit = 50.0;
  const Float_t xDelta = 0.02, yDelta = 0.02, zDelta = 0.2;
  fVertexXY = new TH2F("VertexXY","Vertex Diamond (Y vs X)",
		       2*(Int_t)(xLimit/xDelta),-xLimit,xLimit,
		       2*(Int_t)(yLimit/yDelta),-yLimit,yLimit);
  fVertexXY->SetXTitle("X , cm");
  fVertexXY->SetYTitle("Y , cm");
  fVertexXY->SetOption("colz");
  fVertexZ  = new TH1F("VertexZ"," Longitudinal Vertex Profile",
		       2*(Int_t)(zLimit/zDelta),-zLimit,zLimit);
  fVertexZ->SetXTitle("Z , cm");
}

//______________________________________________________________________
Bool_t AliITSMeanVertexer::Init() {
  // Initialize filters
  // Initialize geometry
  // Initialize ITS classes
 
  AliGeomManager::LoadGeometry();
  if (!AliGeomManager::ApplyAlignObjsFromCDB("ITS")) return kFALSE;

  AliITSInitGeometry initgeom;
  AliITSgeom *geom = initgeom.CreateAliITSgeom();
  if (!geom) return kFALSE;
  printf("Geometry name: %s \n",(initgeom.GetGeometryName()).Data());

  fDetTypeRec = new AliITSDetTypeRec();
  fDetTypeRec->SetLoadOnlySPDCalib(kTRUE);
  fDetTypeRec->SetITSgeom(geom);
  fDetTypeRec->SetDefaults();
  fDetTypeRec->SetDefaultClusterFindersV2(kTRUE);

  // Initialize filter values to their defaults
  SetFilterOnContributors();
  SetFilterOnTracklets();

  // Instatiate vertexer
  if (!fMode) {
    fVertexer = new AliITSVertexer3DTapan(1000);
  }
  else {
    fVertexer = new AliITSVertexer3D();
    fVertexer->SetDetTypeRec(fDetTypeRec);
    fVertexer->SetComputeMultiplicity(kTRUE);
  }
  return kTRUE;
}

//______________________________________________________________________
AliITSMeanVertexer::~AliITSMeanVertexer() {
  // Destructor
  delete fDetTypeRec;
  delete fVertexXY;
  delete fVertexZ;
  delete fVertexer;
}

//______________________________________________________________________
Bool_t AliITSMeanVertexer::Reconstruct(AliRawReader *rawReader){
  // Performs SPD local reconstruction
  // and vertex finding
  // returns true in case a vertex is found

  // Run SPD cluster finder
  AliITSRecPointContainer::Instance()->PrepareToRead();
  TTree* clustersTree = new TTree("TreeR", "Reconstructed Points Container"); //make a tree
  fDetTypeRec->DigitsToRecPoints(rawReader,clustersTree,"SPD");

  Bool_t vtxOK = kFALSE;
  AliESDVertex *vtx = fVertexer->FindVertexForCurrentEvent(clustersTree);
  if (!fMode) {
    if (TMath::Abs(vtx->GetChi2()) < 0.1) vtxOK = kTRUE;
  }
  else {
    AliMultiplicity *mult = fVertexer->GetMultiplicity();
    if(Filter(vtx,mult)) vtxOK = kTRUE;
  }
  delete clustersTree;
  if (vtxOK) AddToMean(vtx);
  if (vtx) delete vtx;

  return vtxOK;
}

//______________________________________________________________________
void AliITSMeanVertexer::WriteVertices(const char *filename){
  // Compute mean vertex and
  // store it along with the histograms
  // in a file
  
  TFile fmv(filename,"update");

  if(ComputeMean()){
    Double_t cov[6];
    cov[0] =  fAverPosSq[0][0];  // variance x
    cov[1] =  fAverPosSq[0][1];  // cov xy
    cov[2] =  fAverPosSq[1][1];  // variance y
    cov[3] =  fAverPosSq[0][2];  // cov xz
    cov[4] =  fAverPosSq[1][2];  // cov yz
    cov[5] =  fAverPosSq[2][2];  // variance z
    AliMeanVertex mv(fWeighPos,fWeighSig,cov,fNoEventsContr,fTotTracklets,fAverTracklets,fSigmaOnAverTracks);
    mv.SetTitle("Mean Vertex");
    mv.SetName("MeanVertex");
    AliDebug(1,Form("Contrib av. trk = %10.2f ",mv.GetAverageNumbOfTracklets()));
    AliDebug(1,Form("Sigma %10.4f ",mv.GetSigmaOnAvNumbOfTracks()));
    // we have to add chi2 here
    AliESDVertex vtx(fWeighPos,cov,0,TMath::Nint(fAverTracklets),"MeanVertexPos");

    mv.Write(mv.GetName(),TObject::kOverwrite);
    vtx.Write(vtx.GetName(),TObject::kOverwrite);
  }
  else {
    AliError(Form("Evaluation of mean vertex not possible. Number of used events = %d",fNoEventsContr));
  }

  fVertexXY->Write(fVertexXY->GetName(),TObject::kOverwrite);
  fVertexZ->Write(fVertexZ->GetName(),TObject::kOverwrite);
  fmv.Close();
}

//______________________________________________________________________
Bool_t AliITSMeanVertexer::Filter(AliESDVertex *vert,AliMultiplicity *mult){
  // Apply selection criteria to events
  Bool_t status = kFALSE;
  if(!vert || !mult)return status;
  // Remove vertices reconstructed with vertexerZ
  if(strcmp(vert->GetName(),"SPDVertexZ") == 0) return status;
  Int_t ncontr = vert->GetNContributors();
  Int_t ntracklets = mult->GetNumberOfTracklets();
  AliDebug(1,Form("Number of contributors = %d",ncontr));
  AliDebug(1,Form("Number of tracklets = %d",ntracklets));
  if(ncontr>fFilterOnContributors && ntracklets > fFilterOnTracklets) status = kTRUE;
  fTotTracklets += ntracklets;
  fTotTrackletsSq += ntracklets*ntracklets;
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
    fWeighPosSum[i]+=currentPos[i]/currentSigma[i]/currentSigma[i];
    fWeighSigSum[i]+=1./currentSigma[i]/currentSigma[i];
    fAverPosSum[i]+=currentPos[i];
  }
  for(Int_t i=0;i<3;i++){
    for(Int_t j=i;j<3;j++){
      fAverPosSqSum[i][j] += currentPos[i] * currentPos[j];
    }
  }

  fVertexXY->Fill(currentPos[0],currentPos[1]);
  fVertexZ->Fill(currentPos[2]);

  fNoEventsContr++;
}

//______________________________________________________________________
Bool_t AliITSMeanVertexer::ComputeMean(){
  // compute mean vertex
  if(fNoEventsContr < 2) return kFALSE;
  Double_t nevents = fNoEventsContr;
  for(Int_t i=0;i<3;i++){
    fWeighPos[i] = fWeighPosSum[i]/fWeighSigSum[i]; 
    fWeighSig[i] = 1./TMath::Sqrt(fWeighSigSum[i]);
    fAverPos[i] = fAverPosSum[i]/nevents;
  } 
  for(Int_t i=0;i<3;i++){
    for(Int_t j=i;j<3;j++){
      fAverPosSq[i][j] = fAverPosSqSum[i][j]/(nevents -1.);
      fAverPosSq[i][j] -= nevents/(nevents -1.)*fAverPos[i]*fAverPos[j];
    }
  } 

  fAverTracklets = fTotTracklets/nevents;
  fSigmaOnAverTracks = fTotTrackletsSq/(nevents - 1);
  fSigmaOnAverTracks -= nevents/(nevents -1.)*fAverTracklets*fAverTracklets;
  fSigmaOnAverTracks = TMath::Sqrt(fSigmaOnAverTracks);
  return kTRUE;
}
