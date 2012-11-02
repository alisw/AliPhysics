#include "AliITSURecoDet.h"
#include "AliITSUGeomTGeo.h"
#include "AliITSsegmentation.h"
#include "AliITSUSegmentationPix.h"

ClassImp(AliITSURecoDet)

//______________________________________________________
AliITSURecoDet::AliITSURecoDet(const char* name)
:  fNLayers(0)
  ,fNLayersActive(0)
  ,fRMax(-1)
  ,fRMin(-1)
  ,fLayers(0)
  ,fLayersActive(0)
  ,fSegmentations(0)
  ,fITSGeom(0)
{
  // def. c-tor
  SetNameTitle(name,name);
  fLayers.SetOwner(kTRUE);        // layers belong to this array
  fLayersActive.SetOwner(kFALSE); // this one just points on active layers in fLayers
  fSegmentations.SetOwner(kTRUE); // segmentations are owned by the detector
}

//______________________________________________________
AliITSURecoDet::~AliITSURecoDet()
{
  // def. d-tor
  fLayersActive.Clear(); 
  fLayers.Clear();         // owned!
  fSegmentations.Clear();  // owned!
  delete fITSGeom;
}

//______________________________________________________
void AliITSURecoDet::Print(Option_t* opt) const			      
{
  //print 
  printf("Detector %s, %d layers, %d active layers\n",GetName(),GetNLayers(),GetNLayersActive());
  TString opts = opt; opts.ToLower();
  if (opts.Contains("lr")) for (int i=0;i<GetNLayers();i++) GetLayer(i)->Print(opt);
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
  fITSGeom = new AliITSUGeomTGeo(kTRUE);
  int nlr = fITSGeom->GetNLayers();
  if (!nlr) AliFatal("No geometry loaded");
  AliITSUSegmentationPix::LoadSegmentations(&fSegmentations, fITSGeom->GetITSsegmentationFileName());
  if (!fSegmentations.GetEntriesFast()) AliFatal(Form("Segmentations from %s are not loaded",fITSGeom->GetITSsegmentationFileName()));
  //
  // build active ITS layers
  for (int ilr=0;ilr<nlr;ilr++) {
    int lrTyp = fITSGeom->GetLayerDetTypeID(ilr);
    int nLad = fITSGeom->GetNLadders(ilr);
    int nDet = fITSGeom->GetNDetectors(ilr);
    // name layer according its active id, detector type and segmentation tyoe
    AliITSUSegmentationPix* segm = (AliITSUSegmentationPix*)fSegmentations.At(lrTyp);
    if (!segm) { AliFatal(Form("Did not find segmentation type %d",lrTyp)); continue; }
    AliITSURecoLayer* lra = new AliITSURecoLayer( Form("Lr%d%s%d",ilr,fITSGeom->GetDetTypeName(lrTyp),
						       lrTyp%AliITSUGeomTGeo::kMaxSegmPerDetType),
						  ilr,nLad*nDet,fITSGeom,segm);
    lra->Build();
    AddLayer(lra);
  }
  return kTRUE;
}
