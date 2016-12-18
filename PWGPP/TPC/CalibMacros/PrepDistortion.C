#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliTPCDcalibRes.h"
#include <TSting.h>
#endif

// This macro reads aldeady processed AliTPCDcalibRes object from 
// alitpcdcalibres.root file together with corresponding voxel tree from pathToOrig, 
// adds to AliTPCDcalibRes object Cheb. parameterization for the distortions and
// saves it locally. On demand the control tree is produced.

void PrepDistortion(const char* pathToOrig, Bool_t testTree=kFALSE)
{
  TString path = pathToOrig;
  if (path.EndsWith("alitpcdcalibres.root")) path.ReplaceAll("alitpcdcalibres.root","");
  if (path.BeginsWith("alien://") && !gGrid && !TGrid::Connect("alien://")) {
    ::Error("PrepDistortion","kill: Failed to open alien connection");
    exit(1);    
  }
  path = path.Strip(TString::kTrailing,'/');
  AliTPCDcalibRes* clb = AliTPCDcalibRes::Load(path.Data());
  if (!clb) {
    ::Error("PrepDistortion","kill: Failed to extract original calib. object from %s",path.Data());
    exit(1);
  }
  
  // this is to recreate corrections object
  // clb->ReProcessFromResVoxTree(Form("%s/voxelResTree.root",path.Data()), kFALSE);
  //
  // this is to prepare data for distortion object, while keeping existing correction object
  if (! clb->LoadDataFromResTree(Form("%s/voxelResTree.root",path.Data())) ) return; 
  clb->ReProcessResiduals();
  //
  clb->CreateDistortionObject();
  if (testTree) clb->WriteDistCorTestTree();
  clb->Save();
}
