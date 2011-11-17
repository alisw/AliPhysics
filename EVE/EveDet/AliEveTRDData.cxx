// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "TROOT.h"
#include "TStyle.h"
#include "TVector.h"
#include "TLinearFitter.h"
#include "TCanvas.h"
#include "TGeoMatrix.h"
#include "TGeoManager.h"

#include "TEveTrans.h"
#include "TEveManager.h"

#include "EveBase/AliEveEventManager.h"

#include "AliEveTRDData.h"
#include "AliEveTRDModuleImp.h"
#include "AliEveTRDLoader.h"
#include "AliEveTRDLoaderImp.h"

#include "AliGeomManager.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliTrackPointArray.h"
#include "AliRieman.h"

#include "AliTRDhit.h"
#include "AliTRDcluster.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackletMCM.h"
#include "AliTRDtrackletWord.h"
#include "AliTRDmcmSim.h"
#include "AliTRDtrackV1.h"
#include "AliTRDtrackerV1.h"
#include "AliTRDpadPlane.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDmcmSim.h"
#include "AliTRDarrayADC.h"
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
ClassImp(AliEveTRDTrackletOnline)
ClassImp(AliEveTRDmcm)

///////////////////////////////////////////////////////////
/////////////   AliEveTRDDigits       /////////////////////
///////////////////////////////////////////////////////////

//______________________________________________________________________________
AliEveTRDDigits::AliEveTRDDigits(AliEveTRDChamber *p) 
  :TEveQuadSet("digits", "")
  ,fParent(p)
{
  // Constructor.
  SetOwnIds(kTRUE);
  gStyle->SetPalette(1, 0);
  SetPalette(new TEveRGBAPalette(0, 512));
  Reset(TEveQuadSet::kQT_RectangleYZ, kFALSE, 32);
}

//______________________________________________________________________________
AliEveTRDDigits::~AliEveTRDDigits()
{
//  AliInfo(GetTitle());
}

//______________________________________________________________________________
void AliEveTRDDigits::SetData(AliTRDdigitsManager *digits)
{
  // Set data source.

  Int_t det(fParent->GetID());
  AliTRDarrayADC *data = digits->GetDigits(det);
  if(!data->GetDim()) return;
  data->Expand();

  AliTRDSignalIndex *indexes = digits->GetIndexes(det);
  if(!indexes->IsAllocated()) digits->BuildIndexes(det);

  Double_t scale, dy, dz;
  Int_t ly    = AliTRDgeometry::GetLayer(det),
        stk   = AliTRDgeometry::GetStack(det),
        sec   = AliTRDgeometry::GetSector(det),
        vid   = AliGeomManager::LayerToVolUID(AliGeomManager::kTRD1 + ly, stk + AliTRDgeometry::Nstack() * sec);
  SetNameTitle(Form("digits%03d", det), Form("D-%03d [%02d_%d_%d]", det, sec, stk, ly));
  Short_t sig[7]={0,0,0,10,0,0,0};

  AliTRDtransform transform(det);
  AliTRDpadPlane *pp(fParent->fGeo->GetPadPlane(ly, stk));

  Int_t row, col;
  AliTRDcluster c;
  indexes->ResetCounters();
  while (indexes->NextRCIndex(row, col)){
    dz = pp->GetRowSize(row);
    dy = pp->GetColSize(col);
    Short_t *const adc = data->GetDataAddress(row,col);
    for (Int_t time(0); time<30; time++){     
      if(data->IsPadCorrupted(row, col, time)) break;
      if(adc[time] <= 1) continue;
      new (&c) AliTRDcluster(det, col, row, time, sig, vid);
      transform.Transform(&c);

      scale = adc[time] < 512 ? adc[time]/512. : 1.;
      AddQuad(c.GetY()-0.5*dy, c.GetZ()-0.5*dz*scale, c.GetX(), dy*0.95, dz*scale);
      QuadValue(Int_t(adc[time]));
      QuadId(new TNamed(Form("ADC %d", adc[time]), Form("det[%3d(%02d_%d_%d)] col[%3d] row[%2d] tb[%2d]", det, sec, stk, ly, col, row, time)));
    } 
  }

  // rotate to global coordinates
  RefitPlex();
  TEveTrans& t = RefMainTrans();
  t.SetRotByAngles((sec+.5)*AliTRDgeometry::GetAlpha(), 0.,0.);
}


// //______________________________________________________________________________
// void AliEveTRDDigits::Paint(Option_t *option)
// {
//   // Paint the object.
// 
//   if(fParent->GetDigitsBox()) fBoxes.Paint(option);
//   else TEveQuadSet::Paint(option);
// }

// //______________________________________________________________________________
// void AliEveTRDDigits::Reset()
// {
//   // Reset raw and visual data.
// 
//   TEveQuadSet::Reset(TEveQuadSet::kQT_RectangleYZ, kTRUE, 64);
//   // MT fBoxes.fBoxes.clear();
//   fData.Reset();
// }

///////////////////////////////////////////////////////////
/////////////   AliEveTRDHits         /////////////////////
///////////////////////////////////////////////////////////

//______________________________________________________________________________
AliEveTRDHits::AliEveTRDHits() : TEvePointSet("hits", 20)
{
  // Constructor.
  SetMarkerSize(.1);
  SetMarkerColor(kGreen);
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

  AliTRDhit *h = NULL;
  if(!(h = dynamic_cast<AliTRDhit*>(GetPointId(n)))) return;
  printf("Id[%3d] Det[%3d] Reg[%c] TR[%c] Q[%3d] MC[%d] t[%f]\n", 
    n, h->GetDetector(), 
    h->FromAmplification() ? 'A' : 'D', 
    h->FromTRphoton() ? 'y' : 'n', 
    h->GetCharge(), h->GetTrack(), h->GetTime());
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
  if(!c) return;
  c->Print();
  Emit("PointSelected(Int_t)", n);
  // Bool_t	AliCluster::GetGlobalCov(Float_t* cov) const
  // Bool_t	AliCluster::GetGlobalXYZ(Float_t* xyz) const
  // Float_t	AliCluster::GetSigmaY2() const
  // Float_t	AliCluster::GetSigmaYZ() const
  // Float_t	AliCluster::GetSigmaZ2() const
}

//______________________________________________________________________________
void AliEveTRDClusters::Print(Option_t *o) const
{
  AliTRDcluster *c = NULL;

  for(Int_t n = GetN(); n--;){
    if(!(c = dynamic_cast<AliTRDcluster*>(GetPointId(n)))) continue;
    c->Print(o);
  }
}

//______________________________________________________________________________
void AliEveTRDClusters::Load(const Char_t *w) const
{
  Int_t typ = -1;
  if(strcmp(w, "hit")==0) typ = 0;
  else if(strcmp(w, "dig")==0) typ = 1;
  else if(strcmp(w, "cls")==0) typ = 2;
  else if(strcmp(w, "all")==0) typ = 3;
  else{
    AliInfo("The following arguments are accepted:");
    AliInfo("   \"hit\" : loading of MC hits");
    AliInfo("   \"dig\" : loading of digits");
    AliInfo("   \"cls\" : loading of reconstructed clusters");
    AliInfo("   \"all\" : loading of MC hits+digits+clusters");
    return;
  }

  AliTRDcluster *c = NULL;
  Int_t n = 0;
  while((n = GetN() && !(c = dynamic_cast<AliTRDcluster*>(GetPointId(n))))) n++;
  if(!c) return;

  Int_t det = c->GetDetector();
  AliEveTRDLoader *loader = NULL;
  switch(typ){
  case 0:  
    loader = new AliEveTRDLoader("Hits");
    if(!loader->Open("TRD.Hits.root")){ 
      delete loader;
      return;
    }
    loader->SetDataType(AliEveTRDLoader::kTRDHits);
    break;
  case 1:
    loader = new AliEveTRDLoader("Digits");
    if(!loader->Open("TRD.Digits.root")){ 
      delete loader;
      return;
    }
    loader->SetDataType(AliEveTRDLoader::kTRDDigits);
    break;
  case 2:
    loader = new AliEveTRDLoader("Clusters");
    if(!loader->Open("TRD.RecPoints.root")){ 
      delete loader;
      return;
    }
    loader->SetDataType(AliEveTRDLoader::kTRDClusters);
    break;
  case 3:
    loader = new AliEveTRDLoaderSim("MC");
    if(!loader->Open("galice.root")){ 
      delete loader;
      return;
    }
    loader->SetDataType(AliEveTRDLoader::kTRDHits | AliEveTRDLoader::kTRDDigits | AliEveTRDLoader::kTRDClusters);
    break;
  default: return;
  }

  loader->AddChambers(AliTRDgeometry::GetSector(det),AliTRDgeometry::GetStack(det), AliTRDgeometry::GetLayer(det));
  // load first event
  loader->GoToEvent(AliEveEventManager::GetCurrent()->GetEventId());
  
  // register loader with alieve
  gEve->AddElement(loader->GetChamber(det), *(BeginParents()));
  //loader->SpawnEditor();
  gEve->Redraw3D();
}

///////////////////////////////////////////////////////////
/////////////   AliEveTRDTracklet         /////////////////////
///////////////////////////////////////////////////////////

//______________________________________________________________________________
AliEveTRDTracklet::AliEveTRDTracklet(AliTRDseedV1 *trklt):TEveLine()
  ,fClusters(NULL)
{
  // Constructor.
  SetName("tracklet");
  
  if(!gGeoManager){ 
    AliEveEventManager::AssertGeometry();
    AliInfo(Form("gGeo[%p] Closed[%c]", (void*)gGeoManager, gGeoManager->IsClosed()?'y':'n'));
  }
  SetUserData(trklt);
  Float_t dx;
  Float_t x0   = trklt->GetX0();
  Float_t y0   = trklt->GetYref(0);
  Float_t z0   = trklt->GetZref(0);
  Float_t dydx = trklt->GetYref(1);
  Float_t dzdx = trklt->GetZref(1);
  Float_t  tilt = trklt->GetTilt();
  Float_t g[3];
  AliTRDcluster *c = NULL;
  for(Int_t ic=0; ic<AliTRDseedV1::kNclusters; ic++){
    if(!(c = trklt->GetClusters(ic))) continue;
    if(!fClusters) AddElement(fClusters = new AliEveTRDClusters());
    dx = x0 - c->GetX();
    //Float_t yt = y0 - dx*dydx;
    Float_t zt = z0 - dx*dzdx;
    // backup yc - for testing purposes
    Float_t yc = c->GetY(); 
    c->SetY(yc-tilt*(c->GetZ()-zt));
    c->GetGlobalXYZ(g); 
    Int_t id = fClusters->SetNextPoint(g[0], g[1], g[2]);
    c->SetY(yc);
    fClusters->SetPointId(id, new AliTRDcluster(*c));
  } 
  if(fClusters){
    fClusters->SetName("TC clusters");
    fClusters->SetTitle(Form("N[%d]", trklt->GetN()));
    fClusters->SetMarkerColor(kMagenta);
  }

  SetTitle(Form("Det[%d] RC[%c] Layer[%d] P[%7.3f]", trklt->GetDetector(), trklt->IsRowCross()?'y':'n', trklt->GetPlane(), trklt->GetMomentum()));
  SetLineColor(kRed);
  //SetOwnIds(kTRUE);
  
  // init tracklet line
  Int_t sec = AliTRDgeometry::GetSector(trklt->GetDetector());
  Double_t alpha = AliTRDgeometry::GetAlpha() * (sec<9 ? sec + .5 : sec - 17.5); 

  //trklt->Fit(kTRUE);
  y0   = trklt->GetYfit(0);
  dydx = trklt->GetYfit(1);
  Double_t xg =  x0 * TMath::Cos(alpha) - y0 * TMath::Sin(alpha); 
  Double_t yg = x0 * TMath::Sin(alpha) + y0 * TMath::Cos(alpha);
  SetPoint(0, xg, yg, z0); 
  //SetPoint(0, x0, y0, z0);


  dx = .5*AliTRDgeometry::CamHght()+AliTRDgeometry::CdrHght();
  x0 -= dx; 
  y0 -= dydx*dx,
  z0 -= dzdx*dx; 
  xg = x0 * TMath::Cos(alpha) - y0 * TMath::Sin(alpha); 
  yg = x0 * TMath::Sin(alpha) + y0 * TMath::Cos(alpha);
  SetPoint(1, xg, yg, z0);
  //SetPoint(1, x0, y0, z0);
}


AliEveTRDTracklet::~AliEveTRDTracklet() 
{

}

//______________________________________________________________________________
void AliEveTRDTracklet::Print(Option_t *o) const
{
  AliTRDseedV1 *tracklet = (AliTRDseedV1*)GetUserData();
  if(!tracklet) return;
  tracklet->Print(o);
}

///////////////////////////////////////////////////////////
/////////////   AliEveTRDTrack         /////////////////////
///////////////////////////////////////////////////////////

//______________________________________________________________________________
AliEveTRDTrack::AliEveTRDTrack(AliTRDtrackV1 *trk) 
  :TEveLine()
  ,fTrackState(0)
  ,fESDStatus(0)
  ,fAlpha(0.)
  ,fPoints(NULL)
  ,fRim(NULL)
{
  // Constructor.
  SetUserData(trk);
  SetName("");

//  AliTRDtrackerV1::SetNTimeBins(24);

  fRim = new AliRieman(trk->GetNumberOfClusters());
  AliTRDseedV1 *tracklet = NULL;
  for(Int_t il=0; il<AliTRDgeometry::kNlayer; il++){
    if(!(tracklet = trk->GetTracklet(il))) continue;
    if(!tracklet->IsOK()) continue;
    AddElement(new AliEveTRDTracklet(tracklet));

//     tracklet->ResetClusterIter(kFALSE);
//     while((c = tracklet->PrevCluster())){
//    AliTRDcluster *c(NULL);
/*    for(Int_t ic=AliTRDseedV1::kNtb; ic--;){
      if(!(c=tracklet->GetClusters(ic))) continue;
      Float_t xc = c->GetX();
      Float_t yc = c->GetY();
      Float_t zc = c->GetZ();
      Float_t zt = tracklet->GetZref(0) - (tracklet->GetX0()-xc)*tracklet->GetZref(1); 
      yc -= tracklet->GetTilt()*(zc-zt);
      fRim->AddPoint(xc, yc, zc, .05, 2.3);
    }*/
  }
  if(trk->GetNumberOfTracklets()>1) fRim->Update();
  SetStatus(fTrackState);
}

//______________________________________________________________________________
AliEveTRDTrack::~AliEveTRDTrack()
{
  if(fPoints) delete [] fPoints; fPoints = NULL;
  //delete dynamic_cast<AliTRDtrackV1*>(GetUserData());
}

//______________________________________________________________________________
void AliEveTRDTrack::Print(Option_t *o) const
{
  AliTRDtrackV1 *track = (AliTRDtrackV1*)GetUserData();
  if(!track) return;
  track->Print(o);
}

//______________________________________________________________________________
void AliEveTRDTrack::SetStatus(UChar_t s)
{
  // nothing to be done
  if(fPoints && fTrackState == s) return;
  //return;
  const Int_t nc = AliTRDtrackV1::kMAXCLUSTERSPERTRACK;
  AliTRDtrackV1 *trk(NULL);
  if(!(trk=static_cast<AliTRDtrackV1*>(GetUserData()))) {
    AliError("Failed casting data to TRD track.");
    return;
  }

  Bool_t BUILD = kFALSE;
  if(!fPoints){ 
    fPoints = new AliTrackPoint[nc];

    // define the radial span of the track in the TRD
    Double_t xmin = -1., xmax = -1.;
    Int_t det = 0;
    AliTRDseedV1 *trklt = NULL;
    for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
      if(!(trklt = trk->GetTracklet(ily))) continue;
      if(xmin<0.) xmin = trklt->GetX0() - AliTRDgeometry::CamHght() - AliTRDgeometry::CdrHght();
      if(trklt->GetX0()>xmax) xmax = trklt->GetX0();
      det = trklt->GetDetector();
    }
    Int_t sec = det/AliTRDgeometry::kNdets;
    fAlpha = AliTRDgeometry::GetAlpha() * (sec<9 ? sec + .5 : sec - 17.5); //trk->GetAlpha()

    Double_t dx =(xmax - xmin)/nc;
    for(Int_t ip=0; ip<nc; ip++){
      fPoints[ip].SetXYZ(xmin, 0., 0.);
      xmin+=dx;
    }
    BUILD = kTRUE;
  }

  // select track model
  if(BUILD || ((s&12) != (fTrackState&12))){
    if(TESTBIT(s, kTrackCosmics)){
      //printf("Straight track\n");
      AliTRDtrackerV1::FitLine(trk, NULL, kFALSE, nc, fPoints);
    } else {
      if(TESTBIT(s, kTrackModel)){
        //printf("Kalman track\n");
        if(trk->GetNumberOfTracklets() >=4) AliTRDtrackerV1::FitKalman(trk, NULL, kFALSE, nc, fPoints);
      } else { 
        //printf("Rieman track rim[%p] nc[%d]\n", (void*)fRim, nc);
        //if(trk->GetNumberOfTracklets() >=4) AliTRDtrackerV1::FitRiemanTilt(trk, NULL, kTRUE, nc, fPoints);
        Float_t x = 0.;
        for(Int_t ip = nc; ip--;){
          x = fPoints[ip].GetX();
          //printf("%2d x[%f] y[%f] z[%f]\n", ip, x, fRim->GetYat(x), fRim->GetZat(x));
          fPoints[ip].SetXYZ(x, fRim->GetYat(x), fRim->GetZat(x));
        }
      }
    }
  
    Float_t global[3];
    for(Int_t ip=0; ip<nc; ip++){
      fPoints[ip].Rotate(-fAlpha).GetXYZ(global);
      SetPoint(ip, global[0], global[1], global[2]);
      //printf("*** %2d x[%f] y[%f] z[%f]\n", ip, global[0], global[1], global[2]);
    
    }
    SetSmooth(kTRUE);
  }

  // set color
  if(BUILD || ((s&3) != (fTrackState&3))){
    if(TESTBIT(s, kSource)){
      //printf("Source color\n");
      if(fESDStatus&AliESDtrack::kTRDin){
        SetMarkerColor(kGreen);
        SetLineColor(kGreen);
      } else {
        SetMarkerColor(kMagenta);
        SetLineColor(kMagenta);
      }
    } else {
      if(TESTBIT(s, kPID) == AliTRDpidUtil::kLQ){
        //printf("PID color kLQPID\n");
        //trk->GetReconstructor()->SetOption("!nn");
      } else {
        //printf("PID color kNNPID\n");
        //trk->GetReconstructor()->SetOption("nn");
      }
      //trk->CookPID();
  
      Int_t species = 0; Float_t pid = 0.;
      for(Int_t is=0; is<AliPID::kSPECIES; ++is) 
        if(trk->GetPID(is) > pid){
          pid = trk->GetPID(is);
		  species = (AliPID::EParticleType) is;
        }
      switch(species){
      case AliPID::kElectron:
        SetMarkerColor(kRed);
        SetLineColor(kRed);
        break;
      default:
        SetMarkerColor(kBlue);
        SetLineColor(kBlue);
        break;
      }
    }
    SetLineWidth(2);
  }
  
  const Char_t *model = "line";
  if(!TESTBIT(s, kTrackCosmics)){
    if(TESTBIT(s, kTrackModel)) model = "kalman";
    else model = "rieman";
  }
  Int_t species = 0; Float_t pid = 0.;
  for(Int_t is=0; is<AliPID::kSPECIES; is++) 
    if(trk->GetPID(is) > pid){
      pid = trk->GetPID(is);
      species = is;
    }

  SetTitle(Form(
    "Tracklets[%d] Clusters[%d]\n"
    "Reconstruction Source[%s]\n"
    "PID[%4.1f %4.1f %4.1f %4.1f %4.1f]\n"
    "MC[%d]", trk->GetNumberOfTracklets(), trk->GetNumberOfClusters(), fESDStatus&AliESDtrack::kTRDin ? "barrel" : "sa",
    1.E2*trk->GetPID(0), 1.E2*trk->GetPID(1),
    1.E2*trk->GetPID(2), 1.E2*trk->GetPID(3), 1.E2*trk->GetPID(4), trk->GetLabel()));

  if(GetName()){
    char id[6]; snprintf(id, 6, "%s", GetName());
    SetName(Form("%s %s", id, AliPID::ParticleName(species)));
  }

  // save track status
  fTrackState = s;
}

//______________________________________________________________________________
void AliEveTRDTrack::Load(const Char_t *what) const
{
// Spread downwards to tracklets the command "what"

  const AliEveTRDTracklet* trklt(NULL);
  TEveElement::List_ci itrklt=BeginChildren();
  while(itrklt!=EndChildren()){
    if((trklt = dynamic_cast<const AliEveTRDTracklet*>(*itrklt))) trklt->Load(what);
    itrklt++;
  }
}

//______________________________________________________________________________
AliEveTRDTrackletOnline::AliEveTRDTrackletOnline(AliTRDtrackletMCM *tracklet) :
  TEveLine(),
  fDetector(-1),
  fROB(-1),
  fMCM(-1)
{
  AliTRDtrackletMCM *trkl = new AliTRDtrackletMCM(*tracklet);
  SetUserData(trkl);

  fDetector = trkl->GetDetector();
  fROB = trkl->GetROB();
  fMCM = trkl->GetMCM();

  SetName("sim. tracklet");
  SetTitle(Form("Det: %i, ROB: %i, MCM: %i, Label: %i\n0x%08x", 
                trkl->GetDetector(), trkl->GetROB(), trkl->GetMCM(), trkl->GetLabel(),
                trkl->GetTrackletWord()));
  SetLineColor(kGreen);
  SetLineWidth(3);

  AliTRDgeometry geo;
  TGeoHMatrix *matrix = geo.GetClusterMatrix(trkl->GetDetector());

  fDetector = trkl->GetDetector();
  fROB = trkl->GetROB();
  fMCM = trkl->GetMCM();
  
  Float_t length = 3.;
  Double_t x[3];
  Double_t p[3];
  Double_t p2[3];
  x[0] = AliTRDgeometry::AnodePos(); 
  x[1] = trkl->GetY();
  x[2] = trkl->GetLocalZ();

  matrix->LocalToMaster(x, p);
  geo.RotateBack(trkl->GetDetector(), p, p2);
  SetPoint(0, p2[0], p2[1], p2[2]);

  x[0] -= length;
  x[1] -= length * trkl->GetdYdX();
  matrix->LocalToMaster(x, p);
  p[2] *= p[0] / (p[0] + length);
  geo.RotateBack(trkl->GetDetector(), p, p2);
  SetPoint(1, p2[0], p2[1], p2[2]);
}

AliEveTRDTrackletOnline::AliEveTRDTrackletOnline(AliTRDtrackletWord *tracklet) :
  TEveLine(),
  fDetector(-1),
  fROB(-1),
  fMCM(-1)
{
  AliTRDtrackletWord *trkl = new AliTRDtrackletWord(*tracklet);
  SetUserData(trkl);

  fDetector = trkl->GetDetector();
  fROB = trkl->GetROB(); 
  fMCM = trkl->GetMCM(); 

  SetName("raw tracklet");
  SetTitle(Form("Det: %i, ROB: %i, MCM: %i, Label: %i\n0x%08x", 
                trkl->GetDetector(), fROB, fMCM, -1,
                trkl->GetTrackletWord()));
  SetLineColor(kRed);
  SetLineWidth(3);

  AliTRDgeometry geo;
  TGeoHMatrix *matrix = geo.GetClusterMatrix(trkl->GetDetector());

  Float_t length = 3.;
  Double_t x[3];
  Double_t p[3];
  Double_t p2[3];
  x[0] = AliTRDgeometry::AnodePos();
  x[1] = trkl->GetY();
  x[2] = trkl->GetLocalZ();
  
  matrix->LocalToMaster(x, p);
  geo.RotateBack(trkl->GetDetector(), p, p2);
  SetPoint(0, p2[0], p2[1], p2[2]);

  x[0] -= length;
  x[1] -= length * trkl->GetdYdX();
  matrix->LocalToMaster(x, p);
  p[2] *= p[0] / (p[0] + length);
  geo.RotateBack(trkl->GetDetector(), p, p2);
  SetPoint(1, p2[0], p2[1], p2[2]);
}

AliEveTRDTrackletOnline::~AliEveTRDTrackletOnline() 
{
  delete ((AliTRDtrackletBase*) GetUserData());
  SetUserData(0x0);
}

void AliEveTRDTrackletOnline::ShowMCM(Option_t *opt) const
{
  if (fDetector < 0 || fROB < 0 || fMCM < 0)
    return;

  AliEveTRDmcm *evemcm = new AliEveTRDmcm();
  evemcm->Init(fDetector, fROB, fMCM);
  evemcm->LoadDigits();
  evemcm->Draw(opt);

  TEveElementList *mcmlist = NULL;
  if (gEve->GetCurrentEvent()) 
    mcmlist = (TEveElementList*) gEve->GetCurrentEvent()->FindChild("TRD MCMs");
  if (!mcmlist) {
    mcmlist = new TEveElementList("TRD MCMs");
    gEve->AddElement(mcmlist);
  }
  gEve->AddElement(evemcm, mcmlist);
}


AliEveTRDmcm::AliEveTRDmcm() :
  TEveElement(),
  TNamed(),
  fMCM(new AliTRDmcmSim())
{
  SetName("MCM");
  SetTitle("Unknown MCM");
}

AliEveTRDmcm::~AliEveTRDmcm()
{
  delete fMCM;
}

Bool_t AliEveTRDmcm::Init(Int_t det, Int_t rob, Int_t mcm)
{
  SetName(Form("MCM: Det. %i, ROB %i, MCM %i", det, rob, mcm));
  SetTitle(Form("MCM: Det. %i, ROB %i, MCM %i", det, rob, mcm));
  fMCM->Init(det, rob, mcm);
  fMCM->Reset();
  return kTRUE;
}

Bool_t AliEveTRDmcm::LoadDigits()
{
  AliRunLoader *rl = AliEveEventManager::AssertRunLoader();
  return fMCM->LoadMCM(rl, fMCM->GetDetector(), 
                       fMCM->GetRobPos(), fMCM->GetMcmPos());
}

Bool_t AliEveTRDmcm::Filter()
{
  fMCM->Filter();
  return kTRUE;
}

Bool_t AliEveTRDmcm::Tracklet()
{
  fMCM->Tracklet();
  return kTRUE;
}

void AliEveTRDmcm::Draw(Option_t* option)
{
  const char *mcmname = Form("mcm_%i_%i_%i", fMCM->GetDetector(),
                             fMCM->GetRobPos(), fMCM->GetMcmPos());

  TCanvas *c = dynamic_cast<TCanvas*> (gROOT->FindObject(mcmname));
  if (!c)
    c = gEve->AddCanvasTab("TRD MCM");
  c->SetTitle(Form("MCM %i on ROB %i of det. %i", 
                   fMCM->GetMcmPos(), fMCM->GetRobPos(), fMCM->GetDetector()));
  c->SetName(mcmname);
  c->cd();
  fMCM->Draw(option);
}

Bool_t AliEveTRDmcm::AssignPointer(const char* ptrname)
{
  gROOT->ProcessLine(Form("AliTRDmcmSim* %s = (AliTRDmcmSim *)%p", ptrname, (void*)fMCM));
  return kTRUE;
}
