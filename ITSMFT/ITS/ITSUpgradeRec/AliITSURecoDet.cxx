#include <TGeoVolume.h>
#include <TGeoTube.h>
#include <TGeoManager.h>
#include <TTree.h>
#include "AliITSURecoDet.h"
#include "AliITSUGeomTGeo.h"
#include "AliITSUClusterPix.h"
#include "AliITSUReconstructor.h"



ClassImp(AliITSURecoDet)


const Char_t* AliITSURecoDet::fgkBeamPipeVolName = "IP_PIPE";

//______________________________________________________
AliITSURecoDet::AliITSURecoDet()
:  fNLayers(0)
  ,fNLayersActive(0)
  ,fRMax(-1)
  ,fRMin(-1)
  ,fRITSTPCRef(-1)
  ,fLayers(0)
  ,fLayersActive(0)
  ,fGeom(0)
{
  // def. c-tor
}

//______________________________________________________
AliITSURecoDet::AliITSURecoDet(AliITSUGeomTGeo* geom, const char* name)
:  fNLayers(0)
  ,fNLayersActive(0)
  ,fRMax(-1)
  ,fRMin(-1)
  ,fRITSTPCRef(-1)
  ,fLayers(0)
  ,fLayersActive(0)
  ,fGeom(geom)
{
  // def. c-tor
  SetNameTitle(name,name);
  fLayers.SetOwner(kTRUE);        // layers belong to this array
  fLayersActive.SetOwner(kFALSE); // this one just points on active layers in fLayers
  Build();
}

//______________________________________________________
AliITSURecoDet::~AliITSURecoDet()
{
  // def. d-tor
  fLayersActive.Clear(); 
  fLayers.Clear();         // owned!
}

//______________________________________________________
void AliITSURecoDet::Print(Option_t* opt) const			      
{
  //print 
  printf("Detector %s, %d layers, %d active layers\n",GetName(),GetNLayers(),GetNLayersActive());
  TString opts = opt; opts.ToLower();
  if (opts.Contains("lr")) {
    for (int i=0;i<GetNLayers();i++) GetLayer(i)->Print(opt);
    printf("ITS-TPC matching reference R: %.3f\n",fRITSTPCRef);
  }
}

//______________________________________________________
void AliITSURecoDet::AddLayer(const AliITSURecoLayer* lr)
{
  //add new layer
  fLayers.AddLast((TObject*)lr);
  fNLayers++;
  if (lr->IsActive()) {
    fLayersActive.AddLast((TObject*)lr);
    fNLayersActive++;
  }
}

//______________________________________________________
Bool_t AliITSURecoDet::Build()
{
  // build detector from TGeo
  //
  if (!fGeom) AliFatal("Geometry interface is not set");
  int nlr = fGeom->GetNLayers();
  if (!nlr) AliFatal("No geometry loaded");
  //
  // build active ITS layers
  for (int ilr=0;ilr<nlr;ilr++) {
    int lrTyp = fGeom->GetLayerChipTypeID(ilr);
    // name layer according its active id, detector type and segmentation tyoe
    AliITSURecoLayer* lra = new AliITSURecoLayer(Form("Lr%d%s%d",ilr,fGeom->GetChipTypeName(lrTyp),
						      lrTyp%AliITSMFTAux::kMaxSegmPerChipType),
						 ilr,fGeom);
    lra->SetPassive(kFALSE);
    AddLayer(lra);
  }
  //
  // build passive ITS layers
  //
  double rMin,rMax,zMin,zMax;
  // beam pipe
  TGeoVolume *v = gGeoManager->GetVolume(fgkBeamPipeVolName);
  AliITSURecoLayer* lrp = 0;
  if (!v) AliWarning("No beam pipe found in geometry");
  else {
    TGeoTube *t=(TGeoTube*)v->GetShape();
    rMin = t->GetRmin();
    rMax = t->GetRmax();
    zMin =-t->GetDz();
    zMax = t->GetDz();
    lrp = new AliITSURecoLayer("BeamPipe");
    lrp->SetRMin(rMin);
    lrp->SetRMax(rMax);
    lrp->SetR(0.5*(rMin+rMax));
    lrp->SetZMin(zMin);
    lrp->SetZMax(zMax);
    lrp->SetPassive(kTRUE);
    AddLayer(lrp);
    //
  }
  //
  // TPC-ITS wall
  const AliITSURecoParam* recopar = AliITSUReconstructor::GetRecoParam();
  if (recopar) {
    lrp = new AliITSURecoLayer("TPC-ITSwall");
    lrp->SetRMin(AliITSUReconstructor::GetRecoParam()->GetTPCITSWallRMin());
    lrp->SetRMax(AliITSUReconstructor::GetRecoParam()->GetTPCITSWallRMax());
    lrp->SetR(0.5*(lrp->GetRMin()+lrp->GetRMax()));
    lrp->SetZMin(-AliITSUReconstructor::GetRecoParam()->GetTPCITSWallZSpanH());
    lrp->SetZMax( AliITSUReconstructor::GetRecoParam()->GetTPCITSWallZSpanH());
    lrp->SetMaxStep( AliITSUReconstructor::GetRecoParam()->GetTPCITSWallMaxStep());
    lrp->SetPassive(kTRUE);
    AddLayer(lrp);
  }
  else {
    AliWarning("RecoParam is not available, TPC-ITS wall is not set");
  }
  //
  IndexLayers();
  Print("lr");
  return kTRUE;
}

//______________________________________________________
void AliITSURecoDet::IndexLayers()
{
  // sort and index layers
  const Double_t kRMargin = 1e-2; // 100 micron margin
  fLayersActive.Sort();
  for (int i=0;i<fNLayersActive;i++) GetLayerActive(i)->SetActiveID(i);
  fLayers.Sort();
  for (int i=0;i<fNLayers;i++) GetLayer(i)->SetID(i);
  if (fNLayers>0) {
    SetRMin(GetLayer(0)->GetRMin()-kRMargin);
    SetRMax(GetLayer(fNLayers-1)->GetRMax()+kRMargin);
  }
  //
  // define the reference R for ITS/TPC matching: outside of last layer but before TPC materials
  const double kOffsLastActR = 5.; // offset of reference layer wrt last active R
  int lastActive = GetLrIDActive(fNLayersActive-1);
  AliITSURecoLayer* lrA = GetLayer(lastActive); // last active
  double rref = lrA->GetRMax() + kOffsLastActR;
  //
  if (lastActive <  fNLayers-1) { // there are material layers outside ...
    AliITSURecoLayer* lrL = GetLayer(lastActive+1);
    if (lrL->GetRMin()<=rref) rref = lrL->GetRMin();
    if (rref - lrA->GetRMax()<kOffsLastActR) {
      AliError(Form("The ITS-TPC matching reference R=%.2f is too close to last active R=%.3f",rref,lrA->GetRMax()));
    }
  }
  SetRITSTPCRef(rref);
  //
}

//______________________________________________________
Int_t AliITSURecoDet::FindLastLayerID(Double_t r, int dir) const
{
  // find the last layer which the particle moving in direction dir (1:outward,-1:inward) 
  // will traverse on its way to radius r 
  int ilr;
  //
  if (dir>0) {
    for (ilr=0;ilr<fNLayers;ilr++) {
      AliITSURecoLayer* lr = GetLayer(ilr);
      if ( r<lr->GetR(-dir) ) break;  // this layer at least entered
    }
    return --ilr;  // -1 will correspond to point below the smalles layer
  }
  else {
    for (ilr=fNLayers;ilr--;) {
      AliITSURecoLayer* lr = GetLayer(ilr);
      if ( r>lr->GetR(-dir) ) break; // this layer at least entered
    }
    ilr++;
    return ilr<fNLayers ? ilr:-1; // -1 will correspond to point above outer layer
  }
  //
}

//______________________________________________________
Int_t AliITSURecoDet::FindFirstLayerID(Double_t r, int dir) const
{
  // find the first layer which the particle moving in direction dir (1:outward,-1:inward) 
  // will traverse starting from radius r 
  int ilr;
  //
  if (dir>0) {
    for (ilr=0;ilr<fNLayers;ilr++) {
      AliITSURecoLayer* lr = GetLayer(ilr);
      if ( r<lr->GetR(dir) ) break;  // this layer at least entered
    }
    return ilr<fNLayers ? ilr:-1;  // -1 will correspond to point above outer leayer
  }
  else {
    for (ilr=fNLayers;ilr--;) {
      AliITSURecoLayer* lr = GetLayer(ilr);
      if ( r>lr->GetR(dir) ) break; // this layer at least entered
    }
    return ilr; // -1 will correspond to point below inner layer
  }
  //
}

//______________________________________________________
void AliITSURecoDet::CreateClusterArrays()
{
  // create cluster arrays for active layers
  for (int ilr=0;ilr<fNLayersActive;ilr++) {
    AliITSURecoLayer*  lr = GetLayerActive(ilr);
    lr->SetOwnsClusterArray(kTRUE);
    int tpDet = fGeom->GetLayerChipTypeID(ilr)/AliITSMFTAux::kMaxSegmPerChipType;
    //
    if (tpDet == AliITSMFTAux::kChipTypePix) {
      lr->SetClusters(new TClonesArray(AliITSUClusterPix::Class()));
    }
    else {
      AliFatal(Form("Unknown detector type %d",tpDet));
    }
    //
  }
  //
}

//_____________________________________________________________________________
Int_t AliITSURecoDet::LoadClusters(TTree* treeRP) 
{
  // read clusters from the tree, if it is provided
  if (!treeRP) return 0;
  for (int ilr=fNLayersActive;ilr--;) {
    TBranch* br = treeRP->GetBranch(Form("ITSRecPoints%d",ilr));
    if (!br) AliFatal(Form("Provided cluster tree does not contain branch for layer %d",ilr));
    br->SetAddress( GetLayerActive(ilr)->GetClustersAddress() );
  }
  return treeRP->GetEntry(0); // we are still in 1 ev/tree mode...
}

Int_t FindIndex(AliITSUClusterPix **parr, Int_t i, AliITSUClusterPix *c) {
  Int_t b=0;
  AliITSUClusterPix *cc=parr[b];
  if (c->Compare(cc) <= 0) return b;

  Int_t e=b+i-1;
  cc=parr[e];
  if (c->Compare(cc) >  0) return e+1;

  Int_t m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    cc=parr[m];
    if (c->Compare(cc) > 0) b=m+1;
    else e=m; 
  }
  return m;   
}

//______________________________________________________
void AliITSURecoDet::SortClusters(AliITSUClusterPix::SortMode_t mode)
{
  // process clsuters according to requested mode
  AliITSUClusterPix::SetSortMode( mode );
  for (int ilr=fNLayersActive;ilr--;) {
      TClonesArray *arr=GetLayerActive(ilr)->GetClusters();

      const Int_t kNcl=arr->GetEntriesFast();
      AliITSUClusterPix **parr=new AliITSUClusterPix*[kNcl];
      parr[0]=(AliITSUClusterPix *)arr->UncheckedAt(0);
      for (Int_t j=1; j<kNcl; j++) {
	  AliITSUClusterPix *c=(AliITSUClusterPix*)arr->UncheckedAt(j);
          Int_t i=FindIndex(parr,j,c);
          Int_t k=j-i;
          memmove(parr+i+1 ,parr+i,k*sizeof(AliITSUClusterPix*));
          parr[i]=c;
      }

      TObject **cont=arr->GetObjectRef();
      for (Int_t i=0; i<kNcl; i++) cont[i]=parr[i];
      delete[] parr;

  }
}

