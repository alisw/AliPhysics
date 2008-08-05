// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "TVector.h"
#include "TLinearFitter.h"
#include "TEveTrans.h"

#include "AliEveTRDData.h"
#include "AliEveTRDModuleImp.h"

#include "AliLog.h"
#include "AliPID.h"
#include "AliTrackPointArray.h"

#include "AliTRDhit.h"
#include "AliTRDcluster.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "AliTRDtrackerV1.h"
#include "AliTRDpadPlane.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDdataArrayDigits.h"
#include "AliTRDSignalIndex.h"
#include "AliTRDgeometry.h"
#include "AliTRDtransform.h"
#include "AliTRDReconstructor.h"
#include "AliTRDrecoParam.h"

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

  if(!fData.HasData()){
    return;
  }

  TEveQuadSet::Reset(TEveQuadSet::kQT_RectangleYZ, kTRUE, 64);

  Double_t scale, dy, dz;
  Int_t q, color;
  Int_t nrows = fData.GetNrow(),
        ncols = fData.GetNcol(),
        ntbs  = fData.GetNtime(),
        det   = fParent->GetID();
  Float_t threshold = fParent->GetDigitsThreshold();

  AliTRDtransform transform(det);
  AliTRDgeometry *geo = fParent->fGeo;
  AliTRDpadPlane *pp = geo->GetPadPlane(geo->GetLayer(det), geo->GetStack(det));

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
        QuadColor(color);
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

  Int_t det = fParent->GetID();
  AliTRDdataArrayDigits *data = digits->GetDigits(det);
  if(!data->HasData()) return;
  data->Expand();

  AliTRDSignalIndex *indexes = digits->GetIndexes(det);
  if(!indexes->IsAllocated()) digits->BuildIndexes(det);

  if(!fData.HasData()) fData.Allocate(data->GetNrow(), data->GetNcol(), data->GetNtime());
  fData.Expand();

  Int_t row, col, time, adc;
  indexes->ResetCounters();
  while (indexes->NextRCIndex(row, col)){
    indexes->ResetTbinCounter();
    while (indexes->NextTbinIndex(time)){
      if(data->IsPadCorrupted(row, col, time)){
        // we should mark this position
        break;
      }
      adc = data->GetData(row, col, time);
      if(adc <= 1) continue;
      fData.SetDataUnchecked(row, col, time, adc);
      //fIndex->AddIndexTBin(row,col,time);
      //printf("\tr[%d] c[%d] t[%d] ADC[%d]\n", row, col, time, adc);
    } 
  }
  fData.Compress(1);
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
AliEveTRDHits::AliEveTRDHits() : TEvePointSet("hits", 20)
{
  // Constructor.
  SetMarkerSize(.1);
  SetMarkerColor(2);
  SetOwnIds(kTRUE);
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
AliEveTRDClusters::AliEveTRDClusters():AliEveTRDHits()
{
  // Constructor.
  SetName("clusters");

  SetMarkerSize(.4);
  SetMarkerStyle(24);
  SetMarkerColor(kGray);
  SetOwnIds(kTRUE);
}

//______________________________________________________________________________
void AliEveTRDClusters::PointSelected(Int_t n)
{
  // Handle an individual point selection from GL.

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
AliEveTRDTracklet::AliEveTRDTracklet(AliTRDseedV1 *trklt):TEveLine()
  ,fClusters(0x0)
{
  // Constructor.
  SetName("tracklet");
  
  SetUserData(trklt);
  Int_t det = -1, sec;
  Float_t g[3];
  AliTRDcluster *c = 0x0;
  AddElement(fClusters = new AliEveTRDClusters());
  for(Int_t ic=0; ic<35; ic++){
    if(!(c = trklt->GetClusters(ic))) continue;
    det = c->GetDetector();
    c->GetGlobalXYZ(g); 
    Int_t id = fClusters->SetNextPoint(g[0], g[1], g[2]);    
    fClusters->SetPointId(id, new AliTRDcluster(*c));
  } 

  SetTitle(Form("Det[%d] Plane[%d] P[%7.3f]", det, trklt->GetPlane(), trklt->GetMomentum()));
  SetLineColor(kRed);
  //SetOwnIds(kTRUE);
  
  sec = det/30;
  Double_t alpha = AliTRDgeometry::GetAlpha() * (sec<9 ? sec + .5 : sec - 17.5); 
  Double_t x0 = trklt->GetX0(), 
    y0f = trklt->GetYfit(0), 
    ysf = trklt->GetYfit(1),
    z0r = trklt->GetZref(0), 
    zsr = trklt->GetZref(1);
  Double_t xg =  x0 * TMath::Cos(alpha) - y0f * TMath::Sin(alpha); 
  Double_t yg = x0 * TMath::Sin(alpha) + y0f * TMath::Cos(alpha);
  SetPoint(0, xg, yg, z0r);
  //SetPointId(0, new AliTRDseedV1(*trackletObj));
  Double_t x1 = x0-3.5, 
    y1f = y0f - ysf*3.5,
    z1r = z0r - zsr*3.5; 
  xg =  x1 * TMath::Cos(alpha) - y1f * TMath::Sin(alpha); 
  yg = x1 * TMath::Sin(alpha) + y1f * TMath::Cos(alpha);
  SetPoint(1, xg, yg, z1r);
}

//______________________________________________________________________________
void AliEveTRDTracklet::ProcessData()
{
  AliTRDseedV1 *tracklet = (AliTRDseedV1*)GetUserData();
  tracklet->Print();
}

///////////////////////////////////////////////////////////
/////////////   AliEveTRDTrack         /////////////////////
///////////////////////////////////////////////////////////

//______________________________________________________________________________
AliEveTRDTrack::AliEveTRDTrack(AliTRDtrackV1 *trk) : TEveLine()
{
  // Constructor.
  SetName("track");

  SetUserData(trk);
  
  AliTRDtrackerV1::SetNTimeBins(24);
  AliTrackPoint points[AliTRDtrackV1::kMAXCLUSTERSPERTRACK];
  Int_t nc = 0, sec = -1; Float_t alpha = 0.;
  AliTRDcluster *c = 0x0;
  AliTRDseedV1 *tracklet = 0x0;
  for(Int_t il=0; il<AliTRDgeometry::kNlayer; il++){
    if(!(tracklet = trk->GetTracklet(il))) continue;
    if(!tracklet->IsOK()) continue;
    AddElement(fTracklet[il] = new AliEveTRDTracklet(tracklet));

    for(Int_t ic=34; ic>=0; ic--){
      if(!(c = tracklet->GetClusters(ic))) continue;
      if(sec<0){
        sec = c->GetDetector()/30;
        alpha = AliTRDgeometry::GetAlpha() * (sec<9 ? sec + .5 : sec - 17.5); 
      }
      points[nc].SetXYZ(c->GetX(),0.,0.);
      nc++;
    }
  }

  //AliTRDtrackerV1::FitRiemanTilt(trk, 0x0, kTRUE, nc, points);
  //AliTRDtrackerV1::FitKalman(trk, 0x0, kFALSE, nc, points);
  AliTRDtrackerV1::FitLine(trk, 0x0, kFALSE, nc, points);

  Float_t global[3];
  for(Int_t ip=0; ip<nc; ip++){
    points[ip].Rotate(-alpha).GetXYZ(global);
    SetNextPoint(global[0], global[1], global[2]);
  }

  SetMarkerColor(kBlack);
  SetLineColor(kBlack);
  SetSmooth(kTRUE);
}

//______________________________________________________________________________
void AliEveTRDTrack::ProcessData()
{
  AliTRDtrackV1 *track = (AliTRDtrackV1*)GetUserData();
  AliInfo(Form("Clusters[%d]", track->GetNumberOfClusters()));
  
  track->CookPID();
  printf("PIDLQ : "); for(int is=0; is<AliPID::kSPECIES; is++) printf("%s[%5.2f] ", AliPID::ParticleName(is), 1.E2*track->GetPID(is)); printf("\n");
}

