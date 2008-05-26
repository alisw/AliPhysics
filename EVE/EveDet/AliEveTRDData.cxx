// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "TEveTrans.h"

#include "AliEveTRDData.h"
#include "AliEveTRDModuleImp.h"

#include "AliLog.h"

#include "AliTRDhit.h"
#include "AliTRDcluster.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "AliTRDpadPlane.h"
#include "AliTRDgeometry.h"
#include "AliTRDtransform.h"
#include "AliTRDdigitsManager.h"

ClassImp(AliEveTRDHits)
ClassImp(AliEveTRDDigits)
ClassImp(AliEveTRDClusters)
ClassImp(AliEveTRDTracklet)
ClassImp(AliEveTRDTrack)

///////////////////////////////////////////////////////////
/////////////   AliEveTRDDigits       /////////////////////
///////////////////////////////////////////////////////////

//______________________________________________________________________________
AliEveTRDDigits::AliEveTRDDigits(AliEveTRDChamber *p) :
  TEveQuadSet("digits", ""), fParent(p), fBoxes(), fData()
{
  // Constructor.
}

//______________________________________________________________________________
AliEveTRDDigits::~AliEveTRDDigits()
{
//  AliInfo(GetTitle());
}

//______________________________________________________________________________
void AliEveTRDDigits::ComputeRepresentation()
{
  // Calculate digits representation according to user settings. The
  // user can set the following parameters:
  // - digits scale (log/lin)
  // - digits threshold
  // - digits apparence (quads/boxes)

  TEveQuadSet::Reset(TEveQuadSet::kQT_RectangleYZ, kTRUE, 64);

  Double_t scale, dy, dz;
  Int_t q, color;
  Int_t nrows = fParent->fNrows,
        ncols = fParent->fNcols,
        ntbs  = fParent->fNtime,
        det   = fParent->GetID();
  Float_t threshold = fParent->GetDigitsThreshold();

  AliTRDtransform transform(det);
  AliTRDgeometry *geo = fParent->fGeo;
  AliTRDpadPlane *pp = geo->GetPadPlane(geo->GetPlane(det), geo->GetChamber(det));

  // express position in tracking coordinates
  fData.Expand();
  for (Int_t ir = 0; ir < nrows; ir++) {
    dz = pp->GetRowSize(ir);
    for (Int_t ic = 0; ic < ncols; ic++) {
      dy = pp->GetColSize(ic);
      for (Int_t it = 0; it < ntbs; it++) {
        q = fData.GetDataUnchecked(ir, ic, it);
        if (q < threshold) continue;

        Double_t x[6] = {0., 0., Double_t(q), 0., 0., 0.}; 
        Int_t  roc[3] = {ir, ic, 0}; 
        Bool_t    out = kTRUE;
        transform.Transform(&x[0], &roc[0], UInt_t(it), out, 0);

        scale = q < 512 ? q/512. : 1.;
        color  = 50+int(scale*50.);
    
        AddQuad(x[1]-.45*dy, x[2]-.5*dz*scale, x[0], .9*dy, dz*scale);
        QuadValue(q);
        QuadColor(Color_t(color));
        QuadId(new TNamed(Form("Charge%d", q), "dummy title"));
      }  // end time loop
    }  // end col loop
  }  // end row loop
  fData.Compress(1);
  
  // rotate to global coordinates
  //RefitPlex();
  TEveTrans& t = RefMainTrans();
  t.SetRotByAngles((geo->GetSector(det)+.5)*AliTRDgeometry::GetAlpha(), 0.,0.);
}

//______________________________________________________________________________
void AliEveTRDDigits::SetData(AliTRDdigitsManager *digits)
{
  // Set data source.

  fData.Allocate(fParent->fNrows, fParent->fNcols, fParent->fNtime);
  //	digits->Expand();
  for (Int_t  row = 0;  row <  fParent->fNrows;  row++)
    for (Int_t  col = 0;  col <  fParent->fNcols;  col++)
      for (Int_t time = 0; time < fParent->fNtime; time++) {
        if(digits->GetDigitAmp(row, col, time, fParent->GetID()) < 0) continue;
        fData.SetDataUnchecked(row, col, time, digits->GetDigitAmp(row, col, time, fParent->GetID()));
      }
}


//______________________________________________________________________________
void AliEveTRDDigits::Paint(Option_t *option)
{
  // Paint the object.

  if(fParent->GetDigitsBox()) fBoxes.Paint(option);
  else TEveQuadSet::Paint(option);
}

//______________________________________________________________________________
void AliEveTRDDigits::Reset()
{
  // Reset raw and visual data.

  TEveQuadSet::Reset(TEveQuadSet::kQT_RectangleYZ, kTRUE, 64);
  // MT fBoxes.fBoxes.clear();
  fData.Reset();
}

///////////////////////////////////////////////////////////
/////////////   AliEveTRDHits         /////////////////////
///////////////////////////////////////////////////////////

//______________________________________________________________________________
AliEveTRDHits::AliEveTRDHits(AliEveTRDChamber *p) :
  TEvePointSet("hits", 20), fParent(p)
{
  // Constructor.
  SetTitle(Form("Hits for Det %d", p->GetID()));
}

//______________________________________________________________________________
AliEveTRDHits::~AliEveTRDHits()
{
  //AliInfo(GetTitle());
}

//______________________________________________________________________________
void AliEveTRDHits::PointSelected(Int_t n)
{
  // Handle an individual point selection from GL.

  fParent->SpawnEditor();
  AliTRDhit *h = dynamic_cast<AliTRDhit*>(GetPointId(n));
  printf("\nDetector             : %d\n", h->GetDetector());
  printf("Region of production : %c\n", h->FromAmplification() ? 'A' : 'D');
  printf("TR photon            : %s\n", h->FromTRphoton() ? "Yes" : "No");
  printf("Charge               : %d\n", h->GetCharge());
  printf("MC track label       : %d\n", h->GetTrack());
  printf("Time from collision  : %f\n", h->GetTime());
}


///////////////////////////////////////////////////////////
/////////////   AliEveTRDClusters         /////////////////////
///////////////////////////////////////////////////////////

//______________________________________________________________________________
AliEveTRDClusters::AliEveTRDClusters(AliEveTRDChamber *p):AliEveTRDHits(p)
{
  // Constructor.
  SetName("clusters");
  SetTitle(Form("Clusters for Det %d", p->GetID()));
}

//______________________________________________________________________________
void AliEveTRDClusters::PointSelected(Int_t n)
{
  // Handle an individual point selection from GL.

  fParent->SpawnEditor();
  AliTRDcluster *c = dynamic_cast<AliTRDcluster*>(GetPointId(n));
  printf("\nDetector             : %d\n", c->GetDetector());
  printf("Charge               : %f\n", c->GetQ());
  printf("Sum S                : %4.0f\n", c->GetSumS());
  printf("Time bin             : %d\n", c->GetLocalTimeBin());
  printf("Signals              : ");
  Short_t *cSignals = c->GetSignals();
  for(Int_t ipad=0; ipad<7; ipad++) printf("%d ", cSignals[ipad]); printf("\n");
  printf("Central pad          : %d\n", c->GetPadCol());
  printf("MC track labels      : ");
  for(Int_t itrk=0; itrk<3; itrk++) printf("%d ", c->GetLabel(itrk)); printf("\n");
  // Bool_t	AliCluster::GetGlobalCov(Float_t* cov) const
  // Bool_t	AliCluster::GetGlobalXYZ(Float_t* xyz) const
  // Float_t	AliCluster::GetSigmaY2() const
  // Float_t	AliCluster::GetSigmaYZ() const
  // Float_t	AliCluster::GetSigmaZ2() const
}

///////////////////////////////////////////////////////////
/////////////   AliEveTRDTracklet         /////////////////////
///////////////////////////////////////////////////////////

//______________________________________________________________________________
AliEveTRDTracklet::AliEveTRDTracklet():TEveLine(), AliTRDseedV1()
{
  // Constructor.
  SetName("tracklet");
  //SetTitle(Form("Clusters for Det %d", p->GetID()));
}

///////////////////////////////////////////////////////////
/////////////   AliEveTRDTrack         /////////////////////
///////////////////////////////////////////////////////////

//______________________________________________________________________________
AliEveTRDTrack::AliEveTRDTrack():TEveLine(), AliTRDtrackV1()
{
  // Constructor.
  SetName("track");
  //SetTitle(Form("Clusters for Det %d", p->GetID()));
}

