/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercialf purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

///////////////////////////////////////////////
// TRD alignment task
//
// Some documentation
//
// Authors :
//
////////////////////////////////////////////////

#include "TTreeStream.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"

#include "AliTrackPointArray.h"
#include "AliLog.h"

#include "AliTRDgeometry.h"
#include "AliTRDcluster.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"

#include "AliTRDalignmentTask.h"

ClassImp(AliTRDalignmentTask)

//________________________________________________________
AliTRDalignmentTask::AliTRDalignmentTask()
  :AliTRDrecoTask()
  ,fTree(0x0)
  ,fArray(0x0)
{
  InitFunctorList();
}

AliTRDalignmentTask::AliTRDalignmentTask(char* name)
  :AliTRDrecoTask(name, "TRD alignment")
  ,fTree(0x0)
  ,fArray(0x0)
{
  InitFunctorList();
  DefineOutput(2, TTree::Class());
}

//________________________________________________________
AliTRDalignmentTask::~AliTRDalignmentTask()
{
  if (fArray) delete fArray;
}


//________________________________________________________
void AliTRDalignmentTask::UserCreateOutputObjects()
{
  // spatial resolution

  if(!HasFunctorList()) InitFunctorList();
  OpenFile(2, "RECREATE");

  fTree = new TTree("spTree", "Tree with track space point arrays");
  fTree->Branch("SP","AliTrackPointArray", &fArray);
}


//________________________________________________________
void AliTRDalignmentTask::UserExec(Option_t *opt)
{
// Documentation to come

  AliTRDrecoTask::UserExec(opt);
  PostData(2, fTree);
}


//________________________________________________________
TH1* AliTRDalignmentTask::PlotTrackPoints(const AliTRDtrackV1 *track)
{
// Documentation to come 

  if(track) fkTrack = track;
  if(!fkTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }

  if (fArray) delete fArray;
  fArray = new AliTrackPointArray(fkTrack->GetNumberOfTracklets());

  // Filling the track points array
  Float_t x, y, z;
  AliTrackPoint p; Int_t ip = 0;
  AliTRDseedV1 *tracklet = 0x0;
  for(Int_t il=0; il<AliTRDgeometry::kNlayer; il++){
    if(!(tracklet = fkTrack->GetTracklet(il))) continue;
    if(!tracklet->IsOK()) continue;
    tracklet->Fit(kTRUE);

    x = tracklet->GetX0();
    y = tracklet->GetYfit(0)-tracklet->GetYfit(1)*(tracklet->GetX0()-x);
    z = tracklet->GetZfit(0);
    p.SetXYZ(x,y,z);
    fArray->AddPoint(ip++, &p);
  }
  fTree->Fill();

/*  if(fDebugLevel>=1){
    Float_t yt = fRim.GetYat(x[il]);
    (*fDebugStream) << "TrkltResiduals"
      << "layer="  << il
      << "x="      <<x[il]
      << "y="      <<y[il]
      << "yt="     <<yt
      << "dydx="   << dydx[il]
      << "dy="     << dy
      << "\n";
  }*/
  return 0x0;
}



//________________________________________________________
Bool_t AliTRDalignmentTask::IsIdenticalWithOneOf(AliTrackPoint * const p, AliTrackPointArray *parray, int nmax) {

  // Is the point p identical with one of the points on the list parray?
  // This is a fix for aliroot 4-16-Rev-01 (and before) writing some 
  // spurious unitialized points. 
 
  for (int i=0; i<parray->GetNPoints() && i<nmax; i++) {
    AliTrackPoint pa;
    parray->GetPoint(pa,i);
    //printf("comparing %7.3f with %7.3f\n",p->GetY(),pa.GetY());
    if (p->GetResidual(pa,0)<1e-8) return kTRUE;
    //printf("different\n");
  }
  return kFALSE;
}
