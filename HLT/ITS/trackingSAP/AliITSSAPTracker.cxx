#include "AliITSSAPTracker.h"
#include "AliITSSAPLayer.h"
#include "AliITSRecPoint.h"
#include "AliGeomManager.h"
#include "AliVParticle.h"
#include "AliSymMatrix.h"
//
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliStack.h"
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TFile.h>
#include "AliLog.h"

#include "AliHLTITSTrackPoint.h"
#include "Riostream.h"

ClassImp(AliITSSAPTracker)

const Float_t AliITSSAPTracker::fgkZSpanITS[AliITSSAPTracker::kMaxLrITS] = 
{ 36. ,14.1,14.1,  38., 22.2,29.7, 51.   ,43.1,48.9};

const Float_t AliITSSAPTracker::fgkRLayITS[AliITSSAPTracker::kMaxLrITS] = 
  { 2.94, 3.9,7.6, 11.04, 15.0,23.9, 29.44 ,38.0,43.0};

const Float_t AliITSSAPTracker::fgkRSpanITS[AliITSSAPTracker::kMaxLrITS] = // half span in R
  { 0.04, 0.5,0.5, 0.5, 0.8, 0.8, 0.5 ,0.6,0.6};

const Int_t    AliITSSAPTracker::fgkPassivLrITS[AliITSSAPTracker::kNLrPassive] = 
  {AliITSSAPTracker::kLrBeamPime,AliITSSAPTracker::kLrShield1,AliITSSAPTracker::kLrShield2};

const Int_t    AliITSSAPTracker::fgkActiveLrITS[AliITSSAPTracker::kNLrActive] = 
  {AliITSSAPTracker::kLrSPD1,AliITSSAPTracker::kLrSPD2,
   AliITSSAPTracker::kLrSDD1,AliITSSAPTracker::kLrSDD2,
   AliITSSAPTracker::kLrSSD1,AliITSSAPTracker::kLrSSD2};

const Int_t    AliITSSAPTracker::fgkLr2Active[AliITSSAPTracker::kMaxLrITS] = // conversion to active lr.
  {-1, 0, 1, -1, 2, 3, -1, 4, 5};

const Float_t AliITSSAPTracker::fgkRhoLITS[AliITSSAPTracker::kMaxLrITS] = {
  0.162802, 0.321960,0.354588, 0.274995, 0.193789,0.198168, 0.435372, 0.195828,0.226940};

const Float_t AliITSSAPTracker::fgkX2X0ITS[AliITSSAPTracker::kMaxLrITS] = {
  0.002757, 0.011660,0.012614, 0.006488, 0.007714,0.007916, 0.012689, 0.007849,0.009128};


const Double_t AliITSSAPTracker::fgkClSystYErr2[AliITSSAPTracker::kNLrActive] = 
  {0.0010*0.0010, 0.0030*0.0030, 0.0500*0.0500, 0.0500*0.0500, 0.0020*0.0020, 0.0020*0.0020};

const Double_t AliITSSAPTracker::fgkClSystZErr2[AliITSSAPTracker::kNLrActive] = 
  {0.0050*0.0050, 0.0050*0.0050, 0.0050*0.0050, 0.0050*0.0050, 0.1000*0.1000, 0.1000*0.1000};


const Int_t    AliITSSAPTracker::fgkLrDefBins[AliITSSAPTracker::kNLrActive][2] = // n bins in z, phi
  { {20,20}, {20,20}, {20,20}, {20,20}, {20,20}, {20,20} };

const Float_t AliITSSAPTracker::fgkDefMass = 0.14;
const Int_t   AliITSSAPTracker::fgkDummyLabel = -3141593;

#ifdef _TIMING_
const char* AliITSSAPTracker::fgkSWNames[AliITSSAPTracker::kNSW] = {
  "Total"
  ,"Tracklets"
  ,"Tracks"
  ,"Vertex"
};
#endif


//______________________________________________
AliITSSAPTracker::AliITSSAPTracker() :
  fSPD2Discard()
  ,fTracklets()
  ,fSPD1Tracklet()
  ,fBlacklist(0)
  ,fPhiShift(0.0045)
  ,fSigThetaTracklet(0.025)
  ,fSigPhiTracklet(0.08)
  ,fChi2CutTracklet(1.5)
  ,fPhiShiftSc(0.)
  ,fDThetaTrackletSc(0)
  ,fDPhiTrackletSc(0)
  ,fBz(5.0)
  ,fMaxRSPDVtx(1.5)
  ,fDPhiTol(0.)
  ,fDThSig2Inv(0.)
  ,fDPhSig2Inv(0.)
  //
  ,fMinPt(0.3)
  ,fCurvMax(0)
  ,fZSPD2CutMin(1e9)
  ,fZSPD2CutMax(-1e9)
  ,fMaxChi2Tr2Cl(40)
  ,fAddErr2YspdVtx(0.02*0.02)
  ,fAddErr2ZspdVtx(0.04*0.04)
  //
  ,fMaxDRPhi(1.0)
  ,fMaxDZ(1.0)
  //
  ,fMissChi2Penalty(3)
  ,fMaxMissedLayers(1)
  ,fNTracks(0)
  ,fMaxTrackletsToRunTracking(99999)
  ,fMaxVtxIter(5)
  ,fStopScaleChange(0.8)
  ,fTracks()
  ,fTrackVertex()
  ,fFitVertex(kTRUE)
  //
  ,fSPDVertex(0)
#ifdef _CONTROLH_
  ,fHTrackletMC(0),fHTrackletAll(0),fHTrackletFake(0),fHTrackMC(0),fHTrackAll(0),fHTrackFake(0)
  ,fHVtxDiffXY(0)
  ,fHVtxDiffXMlt(0),fHVtxDiffYMlt(0),fHVtxDiffZMlt(0)
  ,fHVtxPullXMlt(0),fHVtxPullYMlt(0),fHVtxPullZMlt(0)
  ,fHVtxMCSPDDiffXY(0)
  ,fHVtxMCSPDDiffXMlt(0),fHVtxMCSPDDiffYMlt(0),fHVtxMCSPDDiffZMlt(0)
  ,fHVtxMCSPDPullXMlt(0),fHVtxMCSPDPullYMlt(0),fHVtxMCSPDPullZMlt(0)
  ,fHChi2NDFvsPT(0),fHChi2vsNC(0)
  ,fHVtxMltRef(0),fHVtxOKMlt(0),fHVtxDiffZ(0),fHVtxMCSPDDiffZ(0)
#endif
{
  // def. c-tor
  for (int i=kNLrActive;i--;) fLayers[i] = 0;
}

//______________________________________________
AliITSSAPTracker::~AliITSSAPTracker()
{
  // d-tor
  for (int i=0;i<kNLrActive;i++) delete fLayers[i];
}

//______________________________________________
void AliITSSAPTracker::Init()
{
  // init tracker
  //
  if (!AliGeomManager::GetGeometry()) {
    AliGeomManager::LoadGeometry("geometry.root");
    AliGeomManager::ApplyAlignObjsFromCDB("ITS");
  }
  //
  for (int i=0;i<kNLrActive;i++) {
    int iAct = fgkActiveLrITS[i];
    fLayers[i] = new AliITSSAPLayer(i,fgkZSpanITS[iAct]+1,fgkLrDefBins[i][0],fgkLrDefBins[i][1]);
    fSkipLayer[i] = kFALSE;
    fNSigma2[i] = 7*7;
    fYToler2[i] = 0.2*0.2;
    fZToler2[i] = 0.2*0.2;
    fChi2TotCut[i] = 0;
  }  
  fMaxDRPhi = 1.;
  fMaxDZ    = 1.;
  fChi2TotCut[1] = 40; // 2 cl+vtx -> NDF=1
  fChi2TotCut[2] = 40; 
  fChi2TotCut[3] = 30; 
  fChi2TotCut[4] = 35; 
  fChi2TotCut[5] = 40; 
  //
  fMissChi2Penalty = 3;
  fMaxMissedLayers = 1;
  //
  // auxialary precalculated variables
  if (fChi2CutTracklet<0.1) fChi2CutTracklet = 0.1;
  double scl = TMath::Sqrt(fChi2CutTracklet);
  fDThetaTrackletSc = fSigThetaTracklet*scl;
  fDPhiTrackletSc   = fSigPhiTracklet*scl;
  //
  fDThSig2Inv = 1./(fSigThetaTracklet*fSigThetaTracklet);
  fDPhSig2Inv = 1./(fSigPhiTracklet*fSigPhiTracklet);
  //
  fBlacklist = new TBits(100*100);
  //
  SetMaxVtxIter();
  SetStopScaleChange();
  //
#ifdef _TIMING_
  for (int i=kNSW;i--;) {
    fSW[i].Stop();
    fSW[i].Reset();
  }
#endif
  //
#ifdef _CONTROLH_
  BookHistos();
#endif
}

//______________________________________________
void AliITSSAPTracker::ProcessEvent()
{
  // do full reconstruction
#ifdef _TIMING_
  fSW[kSWTotal].Start(0);
  fSW[kSWTracklets].Start(0);
#endif
  //
  fNTracks = 0;  
  FindTracklets();
  //
  if (GetNTracklets()>fMaxTrackletsToRunTracking) return;
  //
#ifdef _TIMING_
  fSW[kSWTracklets].Stop();
  fSW[kSWTracks].Start(0);
#endif
  //
  Tracklets2Tracks();
  RefitInward();
#ifdef _TIMING_
  fSW[kSWTracks].Stop();
  fSW[kSWVertex].Start(0);
#endif
  if (fFitVertex) {
    if (FitTrackVertex()) {
#ifdef _DEBUG_
      printf("FittedVertex: "); fTrackVertex.Print(();
      printf("SPD   Vertex: "); fSPDVertex->Print();
#endif
    }
  }
  //
  //
#ifdef _TIMING_
  fSW[kSWVertex].Stop();
  fSW[kSWTotal].Stop();
  PrintTiming();
#endif
  //
#ifdef _CONTROLH_
  FillRecoStat();
#endif
  /*
  PrintTracklets();
  PrintTracks();  
  if (fSPDVertex) {printf("SPDvtx: "); fSPDVertex->Print();}
  printf("TRKVtx: "); fTrackVertex.Print();
  */
}


//______________________________________________
void AliITSSAPTracker::Clear(Option_t*)
{
  // reset event info
  ClearTracklets();
  ClearTracks();
  for (int i=kNLrActive;i--;) {
    fNClusters[i] = 0;
    if (fLayers[i]) fLayers[i]->Clear();
  }
}

//______________________________________________
void AliITSSAPTracker::ClearTracklets()
{
  // reset tracklets info
  fSPD2Discard.clear();
  fTracklets.clear();
  fSPD1Tracklet.clear();
  if (fBlacklist) fBlacklist->ResetAllBits();
}


//______________________________________________
void AliITSSAPTracker::AddCluster(AliITSRecPoint* cl)
{
  // add cluster to corresponding layer
  if (!cl->Misalign()) AliWarning("Can't misalign this cluster !"); 
  fLayers[cl->GetLayer()]->AddCluster(cl); 
}

//______________________________________________
Bool_t AliITSSAPTracker::FindTracklets()
{
  // find SPD tracklets
  //
  if (!fSPDVertex) {
    //    AliInfo("No SPD vertex set");
    return kFALSE;
  }
  float rv2 = fSPDVertex->GetX()*fSPDVertex->GetX()+fSPDVertex->GetY()*fSPDVertex->GetY();
  if (rv2>fMaxRSPDVtx*fMaxRSPDVtx) {
#ifdef _DEBUG_
    AliInfo("SPD vertex is too far from beam line");
    fSPDVertex->Print();
#endif
    return kFALSE;    
  } 
  fPhiShiftSc = fPhiShift*TMath::Abs(fBz/5.0);
  fDPhiTol = fDPhiTrackletSc + fPhiShiftSc;
  //
  AliITSSAPLayer &spdL1 = *fLayers[kALrSPD1];
  AliITSSAPLayer &spdL2 = *fLayers[kALrSPD2];
  spdL1.SortClusters(fSPDVertex);
  spdL2.SortClusters(fSPDVertex);
  fNClusters[0] = spdL1.GetNClusters();
  fNClusters[1] = spdL2.GetNClusters();
  //
  if (fNClusters[0]<1 || fNClusters[1]<1) return kFALSE;
  //
  fSPD2Discard.resize(fNClusters[1]);
  fSPD1Tracklet.resize(fNClusters[0]);
  //
  fBlacklist->SetBitNumber(TMath::Max(fNClusters[0]*fNClusters[1],10000),kFALSE); // to reserve the space
  //
  int nfound;
  do {
    nfound = 0;
    for (int icl2=fNClusters[1];icl2--;) if (!fSPD2Discard[icl2]) nfound += AssociateClusterOfL2(icl2);
  } while(nfound);
  //
  for (int itr=GetNTracklets();itr--;) CookLabel(fTracklets[itr]); 
  //
  return kTRUE;
}

//______________________________________________
Int_t AliITSSAPTracker::AssociateClusterOfL2(int icl2)
{
  // find SPD1 cluster matching to SPD2 cluster icl2
  AliITSSAPLayer &spdL1 = *fLayers[kALrSPD1];
  AliITSSAPLayer &spdL2 = *fLayers[kALrSPD2];
  AliITSSAPLayer::ClsInfo* cli2 = spdL2.GetClusterInfo(icl2);
  // expected z at SPD1
  float zV = fSPDVertex->GetZ();
  float z2 = cli2->z - zV;
  float tg2Inv = z2/cli2->r;
  float dzt = (1.+tg2Inv*tg2Inv)*fDThetaTrackletSc;
  float dz = dzt*fgkRLayITS[kLrSPD1] + TMath::Abs(tg2Inv)*fgkRSpanITS[kLrSPD1]; // uncertainty from dTheta and from Layer1 R spread
  float zL1 = zV + tg2Inv*fgkRLayITS[kLrSPD1]; // center of expected Z1
  int nsel1 = spdL1.SelectClusters(zL1-dz,zL1+dz, cli2->phi-fDPhiTol,cli2->phi+fDPhiTol);
  if (!nsel1) {
    fSPD2Discard[icl2] = true;
    return 0; // no candidates
  }
  float chiBest = 9999;
  SPDtracklet_t trk;
  trk.id1 = -1;
  int icl1,nCand=0;
  while ( (icl1=spdL1.GetNextClusterInfoID())!=-1) {  // loop over matching clusters of lr1
    if (IsBlacklisted(icl1,icl2)) continue;
    AliITSSAPLayer::ClsInfo* cli1 = spdL1.GetClusterInfo(icl1);
    float z1 = cli1->z - zV;
    float tg1Inv = z1/cli1->r;
    //
    float dTheta = (tg2Inv-tg1Inv)/(1.+tg1Inv*tg1Inv);        // fast check on theta
    if (TMath::Abs(dTheta)>fDThetaTrackletSc) continue;
    //
    float dPhi = cli1->phi - cli2->phi;                       // fast check on phi
    if (dPhi>TMath::Pi()) dPhi = TMath::TwoPi()-dPhi;
    else if (dPhi<-TMath::Pi()) dPhi += TMath::TwoPi();
    double dPhiS = TMath::Abs(dPhi)-fPhiShiftSc;
    if (TMath::Abs(dPhiS)>fDPhiTrackletSc) continue;
    //
    float chi2 = dTheta*dTheta*fDThSig2Inv + dPhiS*dPhiS*fDPhSig2Inv; // check final chi2
    if (chi2>1.) {
      Blacklist(icl1,icl2);
      continue;
    }
    nCand++;
    if (chi2>chiBest) continue;
    // check if cl1 is already associated with better 
    trk.id1 = icl1;
    trk.id2 = icl2;
    trk.dtht = dTheta;
    trk.dphi = dPhi;
    trk.chi2 = chiBest = chi2;
  }
  //
  if (trk.id1!=-1) { // check if there is no better icl1 candidate for icl2
    int oldId = fSPD1Tracklet[trk.id1];
    if (!oldId) { // store new tracklet
      fTracklets.push_back(trk);
      fSPD1Tracklet[trk.id1] = fTracklets.size(); // refer from clusters to tracklet (id+1)
      fSPD2Discard[icl2] = true; // mark as used
      Blacklist(trk.id1,trk.id2);
      return 1;
    }
    SPDtracklet_t& oldTrk = (SPDtracklet_t&)fTracklets[--oldId];
    if (oldTrk.chi2 < trk.chi2) { // previous is better 
      Blacklist(trk.id1,trk.id2);  // shall we blacklist new combination?
      if (nCand==1)  fSPD2Discard[icl2] = true; // there was just 1 candidate and it is discarded
      return 0;
    }
    // new combination is better, overwrite the old one with new one, marking old L2 cluster free
    fSPD2Discard[oldTrk.id2] = false; // mark as free
    fSPD2Discard[icl2] = true; // mark as used
    oldTrk = trk;         // new combination is better, overwrite it with new one
    Blacklist(trk.id1,trk.id2);
    return 1;
  }
  //
  fSPD2Discard[icl2] = true; // no chance to find partner for this cluster
  return 0;
  //
}


//______________________________________________
void AliITSSAPTracker::Tracklets2Tracks()
{
  // try to extend tracklets to outer layers
  int nTrk = GetNTracklets();
  if (!nTrk) return;
  //
  CalcAuxTracking(); // RS??? do we need to repeat this?
  //
  for (int ila=kALrSDD1;ila<kNLrActive;ila++) {
    if (fSkipLayer[ila]) continue;
    fLayers[ila]->SortClusters(0);
    fNClusters[ila] = fLayers[ila]->GetNClusters();
  }
  //
  fTracks.resize(nTrk);

  //
  for (int itr=0;itr<nTrk;itr++) {
    SPDtracklet_t& trlet = fTracklets[itr];
    //
#ifdef _DEBUG_
    printf("TestTracklet %d\t|",itr);
    int stat = GetTrackletMCTruth(trlet);
    //
    int nmiss=0;
    for (int i=2;i<kNLrActive;i++) {
      printf("%c", (stat&(0x1<<i)) ? '*':'-'); 
      if (!(stat&(0x1<<i))) nmiss++;
    }
    printf("|\n");
    PrintTracklet(itr);
#endif
    //
    float zspd2 = fLayers[kALrSPD2]->GetClusterInfo(trlet.id2)->z;
    if (zspd2<fZSPD2CutMin || zspd2>fZSPD2CutMax) continue;
    ITStrack_t &track = fTracks[fNTracks];
    if (!CreateTrack(track, trlet)) continue;
    track.trackletID = itr;
    Bool_t res;
#ifdef _DEBUG_
    double xyz[3];
    track.paramOut.GetXYZAt(0,fBz,xyz);
    printf("process track pt:%f XYZ: %+.4f %+.4f %+.4f\n",track.paramOut.Pt(),xyz[0],xyz[1],xyz[2]);
#endif
    for (int lrID=kLrShield1;lrID<kMaxLrITS;lrID++) {
      res = FollowToLayer(track,lrID) && IsAcceptableTrack(track);
      if (!res) break;
    }
#ifdef _DEBUG_
    printf("%s:%d\n",res ? "OK" : "Fail",nmiss<=fMaxMissedLayers);
#endif
    if (!res) continue;
    track.paramOut.ResetBit(kInvalidBit); // flag that outward fit succeeded
    CookLabel(track);
    fNTracks++;
    //
  }  
}

//______________________________________________
Bool_t AliITSSAPTracker::IsAcceptableTrack(const AliITSSAPTracker::ITStrack_t& /*track*/) const
{
  // check if the track is acceptable
  return kTRUE;
}

//______________________________________________
void AliITSSAPTracker::PrintTrack(const AliITSSAPTracker::ITStrack_t& track) const
{
  // print track info
  printf("Chi2 = %f for %d clusters. Tracklet %d\n",track.chi2,track.ncl,track.trackletID);
  //  
  for (int ilr=0;ilr<kNLrActive;ilr++) {
    if (track.clID[ilr]<0) continue;
    AliITSRecPoint* cl = fLayers[ilr]->GetClusterSorted(track.clID[ilr]);
    printf("L%d #%4d ",ilr,track.clID[ilr]);
    for (int i=0;i<3;i++) printf("%d ",cl->GetLabel(i)); printf("\n");
  }
  track.paramOut.Print();
  track.paramInw.Print();  
}

//______________________________________________
void AliITSSAPTracker::PrintTracklets() const
{
  // print traklets info
  int ntr = fTracklets.size();
  printf("NTracklets: %d\n",ntr);
  printf("Nspd1: %4d Nspd2: %4d, Ntracklets: %d\n",fNClusters[0],fNClusters[1],ntr);
  for (int itr=0;itr<ntr;itr++) PrintTracklet(itr);
  //
}

//______________________________________________
void AliITSSAPTracker::PrintTracklet(Int_t itr) const
{
  // print single tracklet
  const SPDtracklet_t* trk = &fTracklets[itr];
  AliITSRecPoint* cl1 = fLayers[kALrSPD1]->GetClusterSorted(trk->id1);
  AliITSRecPoint* cl2 = fLayers[kALrSPD2]->GetClusterSorted(trk->id2);
  AliITSSAPLayer::ClsInfo_t* cli0 = fLayers[kALrSPD1]->GetClusterInfo(trk->id1);
  printf("#%3d Phi:%+.3f Eta:%+.3f Dphi:%+.3f Dtht:%+.3f Chi2:%.3f | Lbl:",
	 itr,cli0->phi,
	 -TMath::Log(TMath::Tan(TMath::ATan2(cli0->r,cli0->z-fSPDVertex->GetZ())/2.)),
	 trk->dphi,trk->dtht,trk->chi2);
  int lab=-1,lb = -1;
  for (int i=0;i<3;i++) if ( (lb=cl1->GetLabel(i))>=0 ) {if (lab<0)lab=lb; printf(" %5d",lb);} printf("|");
  for (int i=0;i<3;i++) if ( (lb=cl2->GetLabel(i))>=0 ) printf(" %5d",lb); 
  printf("| ->%d\n",trk->label);
  lab = TMath::Abs(trk->label);
  //
  AliStack* stack = 0;
  AliRunLoader* rl = AliRunLoader::Instance();
  if (lab>=0 && rl && (stack=rl->Stack())) {
    TParticle* mctr = stack->Particle(lab);
    if (mctr) {
      TParticlePDG* mctrPDG = mctr->GetPDG();
      if (mctrPDG) {
	double qpt = mctrPDG->Charge()>0 ? mctr->Pt() : -mctr->Pt();
	printf("MCTrack: Prim:%d Vxyz: {%+.4f %+.4f %+.4f} 1/pt: %.3f tgl: %.3f\n",
	       stack->IsPhysicalPrimary(lab),
	       mctr->Vx(),mctr->Vy(),mctr->Vz(),
	       TMath::Abs(qpt)>0 ? 1./qpt : 9999., TMath::Tan(TMath::Pi()/2. - mctr->Theta()));
      }
    }
  }
}


//______________________________________________
void AliITSSAPTracker::PrintTracks() const
{
  // print tracks info
  printf("NTracks: %d\n",fNTracks);
  for (int itr=0;itr<fNTracks;itr++) PrintTrack(fTracks[itr]);
  //
}


//______________________________________________
void AliITSSAPTracker::CalcAuxTracking()
{
  // precalculate auxilarry variables for tracking
  //
  // largest track curvature to search
  const double ztolerEdge = 1.0;
  fCurvMax = TMath::Abs(fBz*kB2C/fMinPt);
  double thMin =-1e9;
  double thMax = 1e9;
  for (int ilA=kNLrActive-1;ilA>kALrSPD2;ilA--) {
    if (GetSkipLayer(ilA)) continue;
    int ilr=fgkActiveLrITS[ilA];
    double r   = fgkRLayITS[ilr] - fgkRSpanITS[ilr];
    double dz = fgkZSpanITS[ilr]+ztolerEdge+fDThetaTrackletSc*r;
    double ri  = 1./r;
    double tmin= (-dz-fSPDVertex->GetZ())*ri;
    double tmax= ( dz-fSPDVertex->GetZ())*ri;
    if (tmin>thMin) thMin = tmin;
    if (tmax<thMax) thMax = tmax;
  }
  double r = fgkRLayITS[kLrSPD2] + fgkRSpanITS[kLrSPD2];
  fZSPD2CutMin = fSPDVertex->GetZ()+thMin*r; // min Z of SPD2 in tracklet to consider tracking
  fZSPD2CutMax = fSPDVertex->GetZ()+thMax*r; // max Z of SPD2 in tracklet to consider tracking
  //
}

//______________________________________________
Bool_t AliITSSAPTracker::CreateTrack(AliITSSAPTracker::ITStrack_t& track, 
				     AliITSSAPTracker::SPDtracklet_t& trlet)
{
  // create track seed from tracklet
  // init track
  track.label = trlet.label;
  //
  AliITSSAPLayer::ClsInfo_t *cli1=fLayers[kALrSPD1]->GetClusterInfo(trlet.id1);
  AliITSSAPLayer::ClsInfo_t *cli2=fLayers[kALrSPD2]->GetClusterInfo(trlet.id2);
  AliITSRecPoint *cl1=fLayers[kALrSPD1]->GetClusterUnSorted(cli1->index);
  AliITSRecPoint *cl2=fLayers[kALrSPD2]->GetClusterUnSorted(cli2->index);
  int det1 = cl1->GetVolumeId()-fLayers[kALrSPD1]->GetVIDOffset();
  int det2 = cl2->GetVolumeId()-fLayers[kALrSPD2]->GetVIDOffset();
  AliITSSAPLayer::ITSDetInfo_t& detInfo1 = fLayers[kALrSPD1]->GetDetInfo(det1);
  AliITSSAPLayer::ITSDetInfo_t& detInfo2 = fLayers[kALrSPD2]->GetDetInfo(det2);
  //
  // crude momentun estimate
  float dx=cli1->x-cli2->x,dy=cli1->y-cli2->y,d=TMath::Sqrt(dx*dx+dy*dy);
  float qptInv = fBz ? 2*TMath::Sin(cli2->phi-cli1->phi)/d/fBz/kB2C : 0; // positive particle goes anticlockwise in B+
  //
  // we initialize the seed in the tracking frame of 1st detector
  float xv= fSPDVertex->GetX()*detInfo1.cosTF + fSPDVertex->GetY()*detInfo1.sinTF;
  float yv=-fSPDVertex->GetX()*detInfo1.sinTF + fSPDVertex->GetY()*detInfo1.cosTF;
  float zv= fSPDVertex->GetZ();
  float par[5] = {yv, zv, (float)TMath::Sin(cli1->phi-detInfo1.phiTF), (cli1->z-zv)/cli1->r, qptInv};
  double covVtx[6]; 
  fSPDVertex->GetCovarianceMatrix(covVtx);
  float cov[15] = {float(covVtx[0]+covVtx[2] + fAddErr2YspdVtx),
		   0, float(covVtx[5] + fAddErr2ZspdVtx),
		   0,0,1,
		   0,0,0,1,
		   0,0,0,0,100*100};
  AliExternalTrackParam& param = track.paramOut;
  param.Set(xv, detInfo1.phiTF, par, cov);
  track.chi2 = 0;   // chi2 at 1st two point is 0
  // go to 1st layer, ignoring the MS (errors are anyway not defined)
  if (!param.PropagateTo(detInfo1.xTF+cl1->GetX(), fBz)) return kFALSE;
  Double_t cpar0[2]={ cl1->GetY(), cl1->GetZ()};
  Double_t ccov0[3]={ cl1->GetSigmaY2() + GetClSystYErr2(kALrSPD1), 0., cl1->GetSigmaZ2() + GetClSystZErr2(kALrSPD1)};
  if (!param.Update(cpar0,ccov0)) return kFALSE;
  if (!param.CorrectForMeanMaterial(fgkX2X0ITS[kLrSPD1],-fgkRhoLITS[kLrSPD1],fgkDefMass)) return kFALSE;
  // go to 2nd layer
  if (!param.Rotate(detInfo2.phiTF)) return kFALSE;
  if (!param.PropagateTo(detInfo2.xTF+cl2->GetX(), fBz)) return kFALSE;
  Double_t cpar1[2]={ cl2->GetY(), cl2->GetZ()};
  Double_t ccov1[3]={ cl2->GetSigmaY2() + GetClSystYErr2(kALrSPD2), 0., cl2->GetSigmaZ2() + GetClSystZErr2(kALrSPD2)};
  track.chi2 += param.GetPredictedChi2(cpar1,ccov1);
  if (!param.Update(cpar1,ccov1)) return kFALSE;
#ifdef _CONTROLH_
  FillTrackingControlHistos(1,track.label,&param,cpar1,ccov1,cl2);
#endif  
  //
  track.clID[0] = trlet.id1;
  track.clID[1] = trlet.id2;
  track.clID[2] = track.clID[3] = track.clID[4] = track.clID[5] = -1;
  track.ncl = 2;
  track.nmiss=0;
  //
  param.SetBit(kInvalidBit); // flag that track is not yer refitted outward 
  track.paramOut.SetBit(kInvalidBit); // flag that track was not refitter inward
  return kTRUE;
}

//______________________________________________
Bool_t AliITSSAPTracker::CrossPassiveLayer(AliExternalTrackParam& param, Int_t lrID)
{
  // cross the layer, applying mat. corrections
  double xStart=param.GetX();
  double xToGo = GetXatLabRLin(param,fgkRLayITS[lrID]);
  if (xToGo<0 || !param.PropagateTo(xToGo,fBz)) return kFALSE;
  double x2x0=fgkX2X0ITS[lrID],xrho=fgkRhoLITS[lrID];
  if (xStart<xToGo) xrho = -xrho; // inward propagation
  return param.CorrectForMeanMaterial(x2x0,xrho,fgkDefMass,kFALSE);
//
}

//______________________________________________
Bool_t AliITSSAPTracker::FollowToLayer(AliITSSAPTracker::ITStrack_t& track, Int_t lrID)
{
  // take track to given layer, searching hits if needed and applying mat. corrections
  int lrIDA = fgkLr2Active[lrID]; // active layer ID
  if (lrIDA<0 || fSkipLayer[lrIDA]) return CrossPassiveLayer(track.paramOut,lrID);
  //
  AliExternalTrackParam trCopy(track.paramOut);
  double xToGo = GetXatLabRLin(trCopy,fgkRLayITS[lrID]); // aproximate X at lrID
  if (!trCopy.PropagateTo(xToGo,fBz)) return kFALSE;
  double xyz[3];
  trCopy.GetXYZ(xyz);
  double phi=TMath::ATan2(xyz[1],xyz[0]),z=trCopy.GetZ();
  // we need track errors in the plane nearly tangential to crossing point
  if (!trCopy.Rotate(phi)) return kFALSE;
  double drphi = TMath::Sqrt(trCopy.GetSigmaY2()*fNSigma2[lrIDA]+fYToler2[lrIDA]);
  if (drphi>fMaxDRPhi) drphi = fMaxDRPhi;
  double dphi = drphi/fgkRLayITS[lrID];
  double dz   = TMath::Sqrt(trCopy.GetSigmaZ2()*fNSigma2[lrIDA]+fZToler2[lrIDA]);
  if (dz>fMaxDZ) dz = fMaxDZ;
  AliITSSAPLayer* lrA = fLayers[lrIDA];
  int nCl = lrA->SelectClusters(z-dz,z+dz,phi-dphi,phi+dphi);
  Bool_t updDone = kFALSE;
  //
#ifdef _DEBUG_
  printf("at Lr%d, Ncl:%d ",lrIDA,nCl);
  trCopy.Print();
#endif
  //
  if (nCl) {
    int icl,iclBest=-1;
    double chi2Best = fMaxChi2Tr2Cl;
    AliITSRecPoint* bestCl = 0;
    AliExternalTrackParam bestTr;
    //
#ifdef _DEBUG_
    int iclt=0;
#endif
    while ( (icl=lrA->GetNextClusterInfoID())!=-1) {
      AliITSSAPLayer::ClsInfo_t *cli = lrA->GetClusterInfo(icl);
      AliITSRecPoint *cl=lrA->GetClusterUnSorted(cli->index);
      int detId = cl->GetVolumeId()-lrA->GetVIDOffset();
      AliITSSAPLayer::ITSDetInfo_t& detInfo = lrA->GetDetInfo(detId);
      trCopy = track.paramOut;
      if (!trCopy.Propagate(detInfo.phiTF, detInfo.xTF+cl->GetX(), fBz)) continue;
      double cpar[2]={ cl->GetY(), cl->GetZ()};
      double ccov[3]={ cl->GetSigmaY2() + GetClSystYErr2(lrIDA) , 0., cl->GetSigmaZ2() + GetClSystZErr2(lrIDA)};
      double chi2cl = trCopy.GetPredictedChi2(cpar,ccov);
      //
#ifdef _DEBUG_      
      float clXYZ[3]; cl->GetGlobalXYZ(clXYZ);
      double trXYZ[3]; trCopy.GetXYZ(trXYZ);
      Float_t xCl, alphaCl; 
      cl->GetXAlphaRefPlane(xCl,alphaCl);
      //
      printf("cl%d Chi2:%.2f Dyz: %+e %+e Err: %e %e %e |Lb:",iclt++,chi2cl, 
	     cl->GetY()-trCopy.GetY(),cl->GetZ()-trCopy.GetZ(),
	     TMath::Sqrt(ccov[0]),ccov[1],TMath::Sqrt(ccov[2])); //TMP
      for (int j=0;j<3;j++) if (cl->GetLabel(j)>=0) printf(" %d",cl->GetLabel(j)); printf("\n");
      printf("CL: X:%.4f Alp:%+.4f XYZ: %+.4f %+.4f %+.4f\n",xCl,alphaCl,clXYZ[0],clXYZ[1],clXYZ[2]);
      printf("TR: X:%.4f Alp:%+.4f XYZ: %+.4f %+.4f %+.4f\n",detInfo.xTF,detInfo.phiTF,trXYZ[0],trXYZ[1],trXYZ[2]);
      trCopy.Print();
#endif
      //
#ifdef _CONTROLH_
      FillTrackingControlHistos(lrIDA,track.label,&trCopy,cpar,ccov,cl);
#endif	
      //
      if (chi2cl>fMaxChi2Tr2Cl) continue;
      //    SaveCandidate(lrIDA,trCopy,chi2cl,icl);  // RS: do we need this?
      if (chi2cl>chi2Best) continue;
      chi2Best = chi2cl;
      iclBest = icl;
      bestCl = cl;
      bestTr = trCopy;
      if (nCl==1) { // in absence of competitors, do the fit on spot
	if (!bestTr.Update(cpar,ccov)) return kFALSE;
	updDone = kTRUE;
      }
    }
#ifdef _DEBUG_
    printf("Lr%d -> %f\n",lrIDA,chi2Best);
#endif
    //
    if (bestCl) {
      if (!updDone) {
	double cpar[2]={ bestCl->GetY(), bestCl->GetZ()};
	double ccov[3]={ bestCl->GetSigmaY2(), 0., bestCl->GetSigmaZ2()}; // RS: add syst errors    
	if (!bestTr.Update(cpar,ccov)) return kFALSE;
	updDone = kTRUE;
      }
      track.paramOut = bestTr;
      track.clID[lrIDA] = iclBest;      
      track.ncl++;
      track.chi2 += chi2Best;      
    }
  }
  //
  if (!updDone) {
    if (++track.nmiss > fMaxMissedLayers)  return kFALSE;
    track.paramOut = trCopy;
    track.chi2 += fMissChi2Penalty;
  }
  //
#ifdef _CONTROLH_
  int ndf = 2*track.ncl-5;
  if (ndf>0) {
    fHChi2vsNC->Fill(track.ncl,track.chi2); 
    if (lrID==kNLrActive-1) fHChi2NDFvsPT->Fill(track.paramOut.Pt(),track.chi2/ndf);
  }
#endif
  if (track.chi2 > GetChi2TotCut(track.ncl+1)) return kFALSE;
  //
  return track.paramOut.CorrectForMeanMaterial(fgkX2X0ITS[lrID],-fgkRhoLITS[lrID],fgkDefMass,kFALSE);
  //
}

//______________________________________________
void AliITSSAPTracker::CookLabel(AliITSSAPTracker::ITStrack_t& track)
{
  // cook mc label for the track
  track.label = fgkDummyLabel;
  if (!track.ncl) return;
  const int kMaxLbPerCl = 3;
  int lbID[kNLrActive*6],lbStat[kNLrActive*6];
  Int_t nLab=0;
  for (int i=kNLrActive;i--;) {
    int clid = track.clID[i];
    if (clid<0) continue;
    AliITSRecPoint* cl = fLayers[i]->GetClusterSorted(clid);
    for (int imc=0;imc<kMaxLbPerCl;imc++) { // labels within single cluster
      int trLb = cl->GetLabel(imc);
      if (trLb<0) break;
      // search this mc track in already accounted ones
      int iLab;
      for (iLab=0;iLab<nLab;iLab++) if (lbID[iLab]==trLb) break;
      if (iLab<nLab) lbStat[iLab]++;
      else {
	lbID[nLab] = trLb;
	lbStat[nLab++] = 1;
      }
    } // loop over given cluster's labels
  } // loop over all clusters
  //
  if (nLab) {
    int maxLab=0;
    for (int ilb=nLab;ilb--;) if (lbStat[maxLab]<lbStat[ilb]) maxLab=ilb;
    track.label = lbStat[maxLab]==track.ncl ? lbID[maxLab] : -lbID[maxLab];
  }
  //
}

//______________________________________________
void AliITSSAPTracker::CookLabel(AliITSSAPTracker::SPDtracklet_t& tracklet)
{
  // cook mc label for the tracklet
  tracklet.label = fgkDummyLabel;
  const int kMaxLbPerCl = 3;
  int lbID[kNLrActive*6],lbStat[kNLrActive*6];
  Int_t nLab=0;
  for (int i=2;i--;) {
    int clid = i ? tracklet.id2 : tracklet.id1;
    AliITSRecPoint* cl = fLayers[i]->GetClusterSorted(clid);
    for (int imc=0;imc<kMaxLbPerCl;imc++) { // labels within single cluster
      int trLb = cl->GetLabel(imc);
      if (trLb<0) break;
      // search this mc track in already accounted ones
      int iLab;
      for (iLab=0;iLab<nLab;iLab++) if (lbID[iLab]==trLb) break;
      if (iLab<nLab) lbStat[iLab]++;
      else {
	lbID[nLab] = trLb;
	lbStat[nLab++] = 1;
      }
    } // loop over given cluster's labels
  } // loop over all clusters
  //
  if (nLab) {
    int maxLab=0;
    for (int ilb=nLab;ilb--;) if (lbStat[maxLab]<lbStat[ilb]) maxLab=ilb;
    tracklet.label = lbStat[maxLab]==2 ? lbID[maxLab] : -lbID[maxLab];
  }
  //
}

//______________________________________________
Double_t AliITSSAPTracker::GetXatLabRLin(AliExternalTrackParam& track, double r)
{
  // X of track circle intersection in current tracking frame, neglecting the curvature
  // Solution of equation (x+d)^2+(y+b*d)^2 - r^2, where x,y are current coordinates of 
  // track and d=X-x0. b = tg(phi)
  //double sn=tr.GetSnp();
  double sn=track.GetSnp();
  if (TMath::Abs(sn)>kAlmost1) return -999;
  double x=track.GetX(), y=track.GetY();
  double cs2=(1.-sn)*(1.+sn), tg=sn/TMath::Sqrt(cs2);
  double t0=x+tg*y, t1=x*x+y*y-r*r, det=t0*t0-t1/cs2;
  if (det<0) return -999; // does not touch circle
  det = TMath::Sqrt(det);
  return x+(det-t0)*cs2;
  //
}

//______________________________________________
Int_t AliITSSAPTracker::GetTrackletMCTruth(AliITSSAPTracker::SPDtracklet_t& trlet) const
{
  int status = 0;
  AliITSSAPLayer::ClsInfo_t *cli1=fLayers[kALrSPD1]->GetClusterInfo(trlet.id1);
  AliITSSAPLayer::ClsInfo_t *cli2=fLayers[kALrSPD2]->GetClusterInfo(trlet.id2);
  AliITSRecPoint *cl1=fLayers[kALrSPD1]->GetClusterUnSorted(cli1->index);
  AliITSRecPoint *cl2=fLayers[kALrSPD2]->GetClusterUnSorted(cli2->index);
  //
  int lab = -1;
  //
  for (int i=0;i<3;i++) {
    int lb1 = cl1->GetLabel(i); 
    if (lb1<0) continue;
    for (int j=0;j<3;j++) {
      int lb2 = cl2->GetLabel(i); 
      if (lb2<0) break;
      if (lb1==lb2) {lab = lb1; break;}
    }
    if (lab>=0) break;
  }
  if (lab<0) return 0;
  //
  for (int ila=kALrSDD1;ila<kNLrActive;ila++) {
    for (int icl=fNClusters[ila];icl--;) {
      AliITSRecPoint *cl=fLayers[ila]->GetClusterUnSorted(icl);
      for (int i=0;i<3;i++) {
	if (cl->GetLabel(i)<0) break;
	if (cl->GetLabel(i)==lab) {status |= 0x1<<ila; break;}
      }
      if (status & (0x1<<ila)) break;
    }
  }
  return status;
}

//______________________________________________
Bool_t AliITSSAPTracker::RefitInward(int itr)
{
  // refit track inward with material correction
  ITStrack_t &track = fTracks[itr];
  AliExternalTrackParam &trout = track.paramOut;
  if (trout.TestBit(kInvalidBit)) return kFALSE;
  AliExternalTrackParam &trin  = track.paramInw;
  trin = trout;
  int ilA = kNLrActive;
  for (;ilA--;) {                    // find outermost layer with cluster
    if (track.clID[ilA]<0) continue;
    break;
  }
  int ilStart = fgkActiveLrITS[ilA]; // corresponding total lr id 
  AliITSSAPLayer* lrA = fLayers[ilA];
  AliITSRecPoint *cl=lrA->GetClusterSorted(track.clID[ilA]);
  AliITSSAPLayer::ITSDetInfo_t& detInfo = lrA->GetDetInfo(cl->GetVolumeId()-lrA->GetVIDOffset());
  if (!trin.RotateParamOnly(detInfo.phiTF)) return kFALSE;
  if (!trin.PropagateParamOnlyTo(detInfo.xTF+cl->GetX(), fBz)) return kFALSE;
  // init with outer cluster y,z and slopes, q/pt of outward track
  double par[5] = {cl->GetY(), cl->GetZ(), trin.GetSnp(), trin.GetTgl(), trin.GetSigned1Pt()}; 
  double cov[15] = {cl->GetSigmaY2() + GetClSystYErr2(kALrSPD1), 
		   0., cl->GetSigmaZ2() + GetClSystZErr2(kALrSPD1),
		   0,0,1,
		   0,0,0,1,
		   0,0,0,0,100*100};
  trin.Set(double(detInfo.xTF+cl->GetX()),double(detInfo.phiTF), par, cov);
  // !!! no material correction is needed: errors are not defined yer
  //
  for (int ilr=ilStart;ilr--;) {
    //
    if ( (ilA=fgkLr2Active[ilr])<0 || track.clID[ilA]<0) { // either passive layer or no cluster
      if (CrossPassiveLayer(trin,ilr)) continue;
      else return kFALSE;
    }
    // there is a cluster, need to update
    lrA = fLayers[ilA];
    cl = lrA->GetClusterSorted(track.clID[ilA]);
    AliITSSAPLayer::ITSDetInfo_t& detInfo1 = lrA->GetDetInfo(cl->GetVolumeId()-lrA->GetVIDOffset());
    if (!trin.Propagate(detInfo1.phiTF, detInfo1.xTF+cl->GetX(), fBz)) return kFALSE;
    double cpar[2]={ cl->GetY(), cl->GetZ()};
    double ccov[3]={ cl->GetSigmaY2() + GetClSystYErr2(ilA) , 0., cl->GetSigmaZ2() + GetClSystZErr2(ilA)};
    if (!trin.Update(cpar,ccov)) return kFALSE;
    //
    // correct for layer materials
    if (!trin.CorrectForMeanMaterial(fgkX2X0ITS[ilr],fgkRhoLITS[ilr],fgkDefMass,kFALSE)) return kFALSE;
    //
  }
  //
  // now go to PCA to vertex
  //double dca[2],dcaCov[3];
  if (!trin.PropagateToDCA(fSPDVertex,fBz,fgkRLayITS[kLrBeamPime])) return kFALSE; //,dca,dcaCov);
  //
  trin.ResetBit(kInvalidBit); // flag that inward fit succeeded
  return kTRUE;
  //
}

//______________________________________________
void AliITSSAPTracker::RefitInward()
{
  // refit tracks inward with material correction
  for (int itr=fNTracks;itr--;) {
    if (!RefitInward(itr)) {
#ifdef _DEBUG_
      printf("RefitInward failed for track %d\n",itr);
      PrintTrack(fTracks[itr]);
#endif
    }
  }
  //
}


//______________________________________________
Bool_t AliITSSAPTracker::FitTrackVertex()
{
  // Fit the vertexTracks. The inner tracks must be already propagated to the SPD vertex.
  // In this case straight line extrapolation can be used
  //
  const double kTiny = 1e-9;
  const double kTukey2 = 6;
  fTrackVertex.SetNContributors(0); // invalidate
  if (fNTracks<3) return kFALSE;
  fTrackVertex.SetXv(fSPDVertex->GetX());
  fTrackVertex.SetYv(fSPDVertex->GetY());
  fTrackVertex.SetZv(fSPDVertex->GetZ());
  //
  double vtxXYZ[3];
  fSPDVertex->GetXYZ(vtxXYZ); // initial vertex
  double scaleSigma2=9;       // initial sigma scaling
  double dz[2],covdum[3],*covt;
  //

#ifdef _DEBUG_
  AliRunLoader* rl = AliRunLoader::Instance();
  AliHeader* hd = 0;
  AliGenEventHeader* hdmc=0;
  TArrayF vtxMC(3);
  if (rl && (hd=rl->GetHeader()) && (hdmc=hd->GenEventHeader())) {
    hdmc->PrimaryVertex(vtxMC);
  }
#endif
  //
  int nIter = 0;
  while(nIter++<fMaxVtxIter) {
    int ntAcc = 0;
    double wghSum=0,wghChi2=0; 
    double cxx=0,cxy=0,cxz=0,cx0=0,cyy=0,cyz=0,cy0=0,czz=0,cz0=0;
    //
    for (int itr=fNTracks;itr--;) {
      //
      AliExternalTrackParam& trc = fTracks[itr].paramInw;
      if (trc.TestBit(kInvalidBit)) continue; // the track is invalidated, skip
      trc.ResetBit(kVtUsedBit);
      //
      double *param = (double*)trc.GetParameter();
      double *covar = (double*)trc.GetCovariance();
      //
      double  x0=trc.GetX();
      double &y0=param[0];
      double &z0=param[1];
      double sn=param[2];
      double cs2=(1.-sn)*(1.+sn);
      if (cs2<kAlmost0) continue;
      double cs=TMath::Sqrt(cs2), tgp=sn/cs, tgl=param[3]/cs;
      // assume straight track equation Y=y0+tgp*X, Z=z0+tgl*X in tracking frame
      //
      double alp = trc.GetAlpha();
      sn = TMath::Sin(alp); // parameters for rotation of vertex to
      cs = TMath::Cos(alp); // tracking frame
      //
      double &syy=covar[0], &syz=covar[1], &szz=covar[2];
      double detI = syy*szz - syz*syz;
      if (TMath::Abs(detI)<kAlmost0) return kFALSE;
      detI = 1./detI;
      double syyI = szz*detI;
      double szzI = syy*detI;
      double syzI =-syz*detI;
      //
      // determine weight of the track
      double vlocX = vtxXYZ[0]*cs+vtxXYZ[1]*sn;
      double vlocY =-vtxXYZ[0]*sn+vtxXYZ[1]*cs;
      double vlocZ = vtxXYZ[2];
      double dy    = y0 + tgp*(vlocX-x0) - vlocY;
      double dz    = z0 + tgl*(vlocX-x0) - vlocZ;
      double chi2T = 0.5*(dy*dy*syyI + dz*dz*szzI) + dy*dz*syzI; 
      double wghT = (1-chi2T/kTukey2/scaleSigma2);
      if (wghT<kTiny)  continue;
      wghSum  += wghT;
      wghChi2 += wghT*chi2T;
      //
      syyI *= wghT;
      syzI *= wghT;
      szzI *= wghT;
      //
      trc.SetBit(kVtUsedBit);
      //
      //      printf("VTXFIT Bef %d X0= %+.4f Z= %+.4f Y=%+.4f\n",itr, x0, z0, y0);
      //
      double tmpSP = sn*tgp;
      double tmpCP = cs*tgp;
      double tmpSC = sn+tmpCP;
      double tmpCS =-cs+tmpSP;
      double tmpCL = cs*tgl;
      double tmpSL = sn*tgl;
      double tmpYXP = y0-tgp*x0;
      double tmpZXL = z0-tgl*x0;
      //
      double tmpCLzz = tmpCL*szzI;
      double tmpSLzz = tmpSL*szzI;
      double tmpSCyz = tmpSC*syzI;
      double tmpCSyz = tmpCS*syzI;
      double tmpCSyy = tmpCS*syyI;
      double tmpSCyy = tmpSC*syyI;
      double tmpSLyz = tmpSL*syzI;
      double tmpCLyz = tmpCL*syzI;
      //
      cxx += tmpCL*(tmpCLzz+tmpSCyz+tmpSCyz)+tmpSC*tmpSCyy;          // dchi^2/dx/dx
      cxy += tmpCL*(tmpSLzz+tmpCSyz)+tmpSL*tmpSCyz+tmpSC*tmpCSyy;    // dchi^2/dx/dy
      cxz += -sn*syzI-tmpCLzz-tmpCP*syzI;                            // dchi^2/dx/dz
      cx0 += -(tmpCLyz+tmpSCyy)*tmpYXP-(tmpCLzz+tmpSCyz)*tmpZXL;     // RHS 
      //
      //double cyx
      cyy += tmpSL*(tmpSLzz+tmpCSyz+tmpCSyz)+tmpCS*tmpCSyy;          // dchi^2/dy/dy
      cyz += -(tmpCSyz+tmpSLzz);                                     // dchi^2/dy/dz
      cy0 += -tmpYXP*(tmpCSyy+tmpSLyz)-tmpZXL*(tmpCSyz+tmpSLzz);     // RHS
      //
      //double czx
      //double czy
      czz += szzI;                                                    // dchi^2/dz/dz
      cz0 += tmpZXL*szzI+tmpYXP*syzI;                                 // RHS
      //
      ntAcc++;
    }
    //
    if (ntAcc<2) break;   // failed
    //
    double vec[3] = {cx0,cy0,cz0};
    AliSymMatrix mat(3);
    mat(0,0) = cxx;
    mat(0,1) = cxy;
    mat(0,2) = cxz;
    mat(1,1) = cyy;
    mat(1,2) = cyz;
    mat(2,2) = czz;
    // 
#ifdef _DEBUG_
    printf("MatBefore: \n"); mat.Print("d");
#endif
    if (!mat.SolveChol(vec,kTRUE)) return kFALSE;
#ifdef _DEBUG_
    printf("MatAfter : \n"); mat.Print("d");
#endif
    //
    double scaleSigma2New = wghChi2/wghSum;
    //
#ifdef _DEBUG_
    double dVtX = vec[0] - vtxXYZ[0];
    double dVtY = vec[1] - vtxXYZ[1];
    double dst2 = dVtX*dVtX+dVtY*dVtY;
    double dVtZ = vec[2] - vtxXYZ[2];
    printf("VTIter%d %d %d  %+e %+e %e %+e %.3f %.3f  %e %e  %e %e %e\n",
	   nIter,ntAcc,fNTracks,dVtX,dVtY,dVtZ,dst2,
	   scaleSigma2,scaleSigma2New, wghChi2,wghSum,
	   vec[0]-vtxMC[0],
	   vec[1]-vtxMC[1],
	   vec[2]-vtxMC[2]
	   );
#endif
    //
    double vtCov[6] = {mat(0,0),mat(0,1),mat(1,1),mat(0,2),mat(1,2),mat(2,2)};
    fTrackVertex.SetXYZ(vec);
    fTrackVertex.SetCovarianceMatrix(vtCov);
    fTrackVertex.SetNContributors(ntAcc);
    //
    if (scaleSigma2<1. && 
	scaleSigma2New/scaleSigma2>fStopScaleChange) break;
    scaleSigma2 = scaleSigma2New;
    for (int i=3;i--;) vtxXYZ[i] = vec[i];
    //
  }  
  // calculate explicitly chi2
  double chiTRC = 0;
  double chiSPD = 0;
  //
  for (int itr=fNTracks;itr--;) {
    AliExternalTrackParam& trc = fTracks[itr].paramInw;
    if (trc.TestBit(kInvalidBit)) continue; // the track is invalidated, skip
    AliExternalTrackParam trT(trc);
    AliExternalTrackParam trS(trc);
    trT.PropagateToDCA(&fTrackVertex,fBz,10,dz,covdum);
    covt = (double*)trT.GetCovariance();
    double detI = covt[0]*covt[2] - covt[1]*covt[1];
    detI = 1./detI;
    double syyI = covt[2]*detI;
    double szzI = covt[0]*detI;
    double syzI =-covt[1]*detI;
    chiTRC += dz[0]*dz[0]*syyI + dz[1]*dz[1]*szzI + 2*dz[0]*dz[1]*syzI;
    //
    trS.PropagateToDCA(fSPDVertex,fBz,10,dz,covdum);
    covt = (double*)trT.GetCovariance();
    detI = covt[0]*covt[2] - covt[1]*covt[1];
    detI = 1./detI;
    syyI = covt[2]*detI;
    szzI = covt[0]*detI;
    syzI =-covt[1]*detI;
    chiSPD += dz[0]*dz[0]*syyI + dz[1]*dz[1]*szzI + 2*dz[0]*dz[1]*syzI;
    //    printf("VTXFITChi2 Aft %d X0= %+.4f Z= %+.4f Y=:%+.4f SPD: X:%+.4f Z:%+.4f Y:%+.4f\n",itr, 
    //	   trT.GetX(), trT.GetZ(), trT.GetY(),
    //	   trS.GetX(), trS.GetZ(), trS.GetY());
  }
#ifdef _DEBUG_    
  /*
  //-------------------------TMP>>>
  AliRunLoader* rl = AliRunLoader::Instance();
  AliHeader* hd = 0;
  AliGenEventHeader* hdmc=0;
  TArrayF vtxMC(3);
  if (rl && (hd=rl->GetHeader()) && (hdmc=hd->GenEventHeader())) {
    hdmc->PrimaryVertex(vtxMC);
  }
  printf("VTFIT %f %f %f %d %8.2f %8.2f   %.4f %.4f %.4f   %.4f %.4f %.4f\n",
	 vtxMC[0],vtxMC[1],vtxMC[2],
	 ntAcc,chiTRC,chiSPD,
	 fTrackVertex.GetX(),fTrackVertex.GetY(),fTrackVertex.GetZ(),
	 fSPDVertex->GetX(),fSPDVertex->GetY(),fSPDVertex->GetZ());
  //-------------------------TMP<<<
  */
#endif
  //
  return kTRUE;
}

#ifdef _CONTROLH_
//______________________________________________
void AliITSSAPTracker::FillRecoStat()
{
  // fill data for preformance study
  //
  AliStack* stack = 0;
  AliRunLoader* rl = AliRunLoader::Instance();
  if (!rl || !(stack=rl->Stack())) return;
  //
  TBits patternMC;
  enum {kIsPrim=kNLrActive,kValidTracklet,kValidTrack,kRecDone,kBitPerTrack};
  int nTrkMC = stack->GetNtrack();
  patternMC.SetBitNumber(nTrkMC*kBitPerTrack,0);
  //
  // fill MC track patterns
  for (int ilr=kNLrActive;ilr--;) {
    AliITSSAPLayer *lr = fLayers[ilr];
    int ncl = lr->GetNClusters();
    for (int icl=ncl;icl--;) {
      AliITSRecPoint* cl = lr->GetClusterUnSorted(icl);
      for (int j=0;j<3;j++) {
	int lb = cl->GetLabel(j);
	if (lb<0 || lb>=nTrkMC) break;
	patternMC.SetBitNumber(lb*kBitPerTrack+ilr,kTRUE);
      }
    }
  }
  // set reconstructability
  for (int itr=nTrkMC;itr--;) {
    int bitoffs = itr*kBitPerTrack;
    Bool_t isPrim = stack->IsPhysicalPrimary(itr);
    patternMC.SetBitNumber(bitoffs+kIsPrim,isPrim);
    if (patternMC.TestBitNumber(bitoffs+kALrSPD1) && patternMC.TestBitNumber(bitoffs+kALrSPD2)) {
      patternMC.SetBitNumber(bitoffs+kValidTracklet,kTRUE);
      //
      TParticle* mctr = stack->Particle(itr);
      fHTrackletMC->Fill(mctr->Pt(),isPrim);
      // check outer layers reconstructability
      int nmiss = 0;
      for (int il=kALrSDD1;il<=kALrSSD2;il++) if (!fSkipLayer[il] && !patternMC.TestBitNumber(bitoffs+il)) nmiss++;
      if (nmiss<=fMaxMissedLayers) {
	patternMC.SetBitNumber(bitoffs+kValidTrack);
	fHTrackMC->Fill(mctr->Pt(),isPrim);
      }
    }
  }
  //
  int nTrk = GetNTracklets();
  if (!nTrk) return;
  for (int itr=0;itr<nTrk;itr++) {
    SPDtracklet_t& trlet = fTracklets[itr];
    //    PrintTracklet(itr);
    //
    int lbl = trlet.label;
    if (lbl==fgkDummyLabel) continue;
    int lblA = TMath::Abs(lbl);
    int bitoffs = lblA*kBitPerTrack;
    Bool_t isPrim = patternMC.TestBitNumber(bitoffs+kIsPrim);
    TParticle* mctr = stack->Particle(lblA);
    double pt = mctr->Pt();
    fHTrackletAll->Fill(pt,isPrim);
    if (lbl<0) fHTrackletFake->Fill(pt,isPrim);
  }
  //
  nTrk = GetNTracks();
  for (int itr=0;itr<nTrk;itr++) {
    ITStrack_t &track = fTracks[itr];
    CookLabel(track);
    //
    int lbl = track.label;
    if (lbl==fgkDummyLabel) continue;
    int lblA = TMath::Abs(lbl);
    int bitoffs = lblA*kBitPerTrack;
    Bool_t isPrim = patternMC.TestBitNumber(bitoffs+kIsPrim);
    TParticle* mctr = stack->Particle(lblA);
    double pt = mctr->Pt();
    Bool_t clone = patternMC.TestBitNumber(bitoffs+kRecDone); // was the track already reconstructed?
    float bn = isPrim ? (clone ? 2:1):(clone ? -1:0);  // fill clones in over/underflow
    fHTrackAll->Fill(pt,bn);
    patternMC.SetBitNumber(bitoffs+kRecDone);
    if (lbl<0) fHTrackFake->Fill(pt,bn);
  }
  //
  AliHeader* hd = rl->GetHeader();
  AliGenEventHeader* hdmc;
  TArrayF vtxMC;
  if (hd && (hdmc=hd->GenEventHeader())) hdmc->PrimaryVertex(vtxMC);
  //
  nTrk = GetNTracklets();
  fHVtxMltRef->Fill(nTrk);
  if (fTrackVertex.GetStatus()==1) {
    if (hdmc) {
      double dx = vtxMC[0]-fTrackVertex.GetX();
      double dy = vtxMC[1]-fTrackVertex.GetY();
      double dz = vtxMC[2]-fTrackVertex.GetZ();
      fHVtxDiffXY->Fill(dx,dy);
      fHVtxDiffZ->Fill(dz);
      fHVtxDiffXMlt->Fill(nTrk, dx);
      fHVtxDiffYMlt->Fill(nTrk, dy);
      fHVtxDiffZMlt->Fill(nTrk, dz);
      //
      double sig[3];
      fTrackVertex.GetSigmaXYZ(sig);    
      if (sig[0]>0) fHVtxPullXMlt->Fill(nTrk, dx/sig[0]);
      if (sig[1]>0) fHVtxPullYMlt->Fill(nTrk, dy/sig[1]);
      if (sig[2]>0) fHVtxPullZMlt->Fill(nTrk, dz/sig[2]);
    }
    fHVtxOKMlt->Fill(nTrk);
  } 
  //
  if (fSPDVertex->GetStatus()==1 && hdmc) {
    double dx = vtxMC[0]-fSPDVertex->GetX();
    double dy = vtxMC[1]-fSPDVertex->GetY();
    double dz = vtxMC[2]-fSPDVertex->GetZ();
    fHVtxMCSPDDiffXY->Fill(dx,dy);
    fHVtxMCSPDDiffZ->Fill(dz);
    fHVtxMCSPDDiffXMlt->Fill(nTrk, dx);
    fHVtxMCSPDDiffYMlt->Fill(nTrk, dy);
    fHVtxMCSPDDiffZMlt->Fill(nTrk, dz);
    //
    double sig[3];
    fSPDVertex->GetSigmaXYZ(sig);    
    if (sig[0]>0) fHVtxMCSPDPullXMlt->Fill(nTrk, dx/sig[0]);
    if (sig[1]>0) fHVtxMCSPDPullYMlt->Fill(nTrk, dy/sig[1]);
    if (sig[2]>0) fHVtxMCSPDPullZMlt->Fill(nTrk, dz/sig[2]);
    //
  }
  //
}

//______________________________________________
void AliITSSAPTracker::BookHistos()
{
  // book control histos
  const int kNBinMlt=20, kNBPt=15, kNBDiffVtx=50, kNResBins=250,kNPullBins=50,kNChiClBins=50,kNBPullVtx=50;
  const double kMinMlt=1,kMaxMlt=5000,kMinPt=0.01,kMaxPt=3, kMaxDiffVtx=0.05, kMaxResidYZ=2.5,kMaxPullYZ=10,kChiClMax=100,kMaxPullVtx=10;
  //
  double* axLogPt  = DefLogAx(kMinPt,kMaxPt,kNBPt);
  double* axLogMlt = DefLogAx(kMinMlt,kMaxMlt,kNBinMlt);

  for (int ilr=0;ilr<kNLrActive;ilr++) {
    //
    // ----------------- These are histos to be filled during tracking
    // PropagateBack and RefitInward will be stored among the histos of 1st pass
    //
    int ilrS = ilr*10;
    TString ttl = Form("residY%d",ilr);
    TH2F* hdy = new TH2F(ttl.Data(),ttl.Data(),kNBPt,axLogPt,kNResBins,-kMaxResidYZ,kMaxResidYZ);
    fArrHisto.AddAtAndExpand(hdy,ilrS+kHResidY);
    hdy->SetDirectory(0);
    //
    ttl = Form("residYPull%d",ilr);	
    TH2F* hdyp = new TH2F(ttl.Data(),ttl.Data(),kNBPt,axLogPt,kNPullBins,-kMaxPullYZ,kMaxPullYZ);
    fArrHisto.AddAtAndExpand(hdyp,ilrS+kHPullY);
    hdyp->SetDirectory(0);
    //
    ttl = Form("residZ%d",ilr);	
    TH2F* hdz = new TH2F(ttl.Data(),ttl.Data(),kNBPt,axLogPt,kNResBins,-kMaxResidYZ,kMaxResidYZ);
    fArrHisto.AddAtAndExpand(hdz,ilrS+kHResidZ);
    hdz->SetDirectory(0);
    //
    ttl = Form("residZPull%d",ilr);		
    TH2F* hdzp = new TH2F(ttl.Data(),ttl.Data(),kNBPt,axLogPt,kNPullBins,-kMaxPullYZ,kMaxPullYZ);
    hdzp->SetDirectory(0);
    fArrHisto.AddAtAndExpand(hdzp,ilrS+kHPullZ);
    //
    ttl = Form("chi2Cl%d",ilr);		
    TH2F* hchi = new TH2F(ttl.Data(),ttl.Data(),kNBPt,axLogPt, kNChiClBins,0.,kChiClMax);
    hchi->SetDirectory(0);
    fArrHisto.AddAtAndExpand(hchi,ilrS+kHChi2Cl);
  } // loop over layers
  //
  fHChi2NDFvsPT = new TH2F("chi2ndfPT","chi2/ndf total vs pt",kNBPt,axLogPt, kNChiClBins,0.,kChiClMax);
  fArrHisto.AddLast(fHChi2NDFvsPT);
  fHChi2NDFvsPT->SetDirectory(0);
  //
  fHChi2vsNC = new TH2F("chi2NC","chi2 total vs NCl",kNLrActive-2,2.5,kNLrActive+0.5, kNChiClBins,0.,kChiClMax);
  fArrHisto.AddLast(fHChi2vsNC);
  fHChi2vsNC->SetDirectory(0);

  // SPDvertex vs MC
  fHVtxMCSPDDiffXY = new TH2F("vtxMCSPDDiffXY","vtxMC-vtxSPD XY",kNBDiffVtx,-kMaxDiffVtx,kMaxDiffVtx,
			    kNBDiffVtx,-kMaxDiffVtx,kMaxDiffVtx);
  fArrHisto.AddLast(fHVtxMCSPDDiffXY);
  fHVtxMCSPDDiffXY->SetDirectory(0);
  //
  fHVtxMCSPDDiffZ = new TH1F("vtxMCSPDDiffZ","vtxMC-vtxSPD Z",kNBDiffVtx,-kMaxDiffVtx,kMaxDiffVtx);
  fArrHisto.AddLast(fHVtxMCSPDDiffZ);
  fHVtxMCSPDDiffZ->SetDirectory(0);
  //
  fHVtxMCSPDDiffXMlt = new TH2F("VtxMCSPDDiffXMlt","vX_{MC}-vX_{SPD} vs mlt",kNBinMlt,axLogMlt, 
			   kNBDiffVtx,-kMaxDiffVtx,kMaxDiffVtx);
  fArrHisto.AddLast(fHVtxMCSPDDiffXMlt);
  fHVtxMCSPDDiffXMlt->SetDirectory(0);
  //
  fHVtxMCSPDDiffYMlt = new TH2F("VtxMCSPDDiffYMlt","vY_{MC}-vY_{rec} vs mlt",kNBinMlt,axLogMlt, 
			   kNBDiffVtx,-kMaxDiffVtx,kMaxDiffVtx);
  fArrHisto.AddLast(fHVtxMCSPDDiffYMlt);
  fHVtxMCSPDDiffYMlt->SetDirectory(0);
  //
  fHVtxMCSPDDiffZMlt = new TH2F("VtxMCSPDDiffZMlt","vZ_{MC}-vZ_{rec} vs mlt",kNBinMlt,axLogMlt, 
			   kNBDiffVtx,-kMaxDiffVtx,kMaxDiffVtx);
  fArrHisto.AddLast(fHVtxMCSPDDiffZMlt);
  fHVtxMCSPDDiffZMlt->SetDirectory(0);
  //
  //
  fHVtxMCSPDPullXMlt = new TH2F("VtxMCSPDPullXMlt","Pull vX_{MC}-vX_{SPD} vs mlt",kNBinMlt,axLogMlt, 
			   kNBPullVtx,-kMaxPullVtx,kMaxPullVtx);
  fArrHisto.AddLast(fHVtxMCSPDPullXMlt);
  fHVtxMCSPDPullXMlt->SetDirectory(0);
  //
  fHVtxMCSPDPullYMlt = new TH2F("VtxMCSPDPullYMlt","Pull vY_{MC}-vY_{rec} vs mlt",kNBinMlt,axLogMlt, 
			   kNBPullVtx,-kMaxPullVtx,kMaxPullVtx);
  fArrHisto.AddLast(fHVtxMCSPDPullYMlt);
  fHVtxMCSPDPullYMlt->SetDirectory(0);
  //
  fHVtxMCSPDPullZMlt = new TH2F("VtxMCSPDPullZMlt","Pull vZ_{MC}-vZ_{rec} vs mlt",kNBinMlt,axLogMlt, 
			   kNBPullVtx,-kMaxPullVtx,kMaxPullVtx);
  fArrHisto.AddLast(fHVtxMCSPDPullZMlt);
  fHVtxMCSPDPullZMlt->SetDirectory(0);
   //
  fHTrackletMC = new TH2F("MCRefTracklet","MCRef Tracklet",kNBPt,axLogPt, 2, -0.5, 1.5);
  fHTrackletMC->SetXTitle("p_{T}");
  fHTrackletMC->GetYaxis()->SetBinLabel(1,"sec");
  fHTrackletMC->GetYaxis()->SetBinLabel(2,"prim");
  fArrHisto.AddLast(fHTrackletMC);
  fHTrackletMC->SetDirectory(0);
  //
  fHTrackletAll = new TH2F("TrackletAll","Tracklet All rec",kNBPt,axLogPt, 2, -0.5, 1.5);
  fHTrackletAll->SetXTitle("p_{T}");
  fHTrackletAll->GetYaxis()->SetBinLabel(1,"sec");
  fHTrackletAll->GetYaxis()->SetBinLabel(2,"prim");
  fArrHisto.AddLast(fHTrackletAll);
  fHTrackletAll->SetDirectory(0);
  //
  fHTrackletFake = new TH2F("TrackletFake","Tracklet Fake rec",kNBPt,axLogPt, 2, -0.5, 1.5);
  fHTrackletFake->SetXTitle("p_{T}");
  fHTrackletFake->GetYaxis()->SetBinLabel(1,"sec");
  fHTrackletFake->GetYaxis()->SetBinLabel(2,"prim");
  fArrHisto.AddLast(fHTrackletFake);
  fHTrackletFake->SetDirectory(0);
  //
  fHTrackMC = new TH2F("MCRefTrack","MCRef Track",kNBPt,axLogPt, 2, -0.5, 1.5);
  fHTrackMC->SetXTitle("p_{T}");
  fHTrackMC->GetYaxis()->SetBinLabel(1,"sec");
  fHTrackMC->GetYaxis()->SetBinLabel(2,"prim");
  fArrHisto.AddLast(fHTrackMC);
  fHTrackMC->SetDirectory(0);
  //
  fHTrackAll = new TH2F("TrackAll","Track All rec",kNBPt,axLogPt, 2, -0.5, 1.5);
  fHTrackAll->SetXTitle("p_{T}");
  fHTrackAll->GetYaxis()->SetBinLabel(1,"sec");
  fHTrackAll->GetYaxis()->SetBinLabel(2,"prim");
  fArrHisto.AddLast(fHTrackAll);
  fHTrackAll->SetDirectory(0);
  //
  fHTrackFake = new TH2F("TrackFake","Track Fake rec",kNBPt,axLogPt, 2, -0.5, 1.5);
  fHTrackFake->SetXTitle("p_{T}");
  fHTrackFake->GetYaxis()->SetBinLabel(1,"sec");
  fHTrackFake->GetYaxis()->SetBinLabel(2,"prim");
  fArrHisto.AddLast(fHTrackFake);
  fHTrackFake->SetDirectory(0);
  //
  fHVtxMltRef = new TH1F("vtxRef","vtxRef",kNBinMlt,axLogMlt);
  fHVtxMltRef->SetXTitle("sqrt(Ntracklets)");
  fArrHisto.AddLast(fHVtxMltRef);
  fHVtxMltRef->SetDirectory(0);
  //
  fHVtxOKMlt = new TH1F("vtxOK","vtxOK",kNBinMlt,axLogMlt);
  fHVtxOKMlt->SetXTitle("sqrt(Ntracklets)");
  fArrHisto.AddLast(fHVtxOKMlt);
  fHVtxOKMlt->SetDirectory(0);
  //
  fHVtxDiffXY = new TH2F("vtxMCDiffXY","vtxMC-vtxRec XY",kNBDiffVtx,-kMaxDiffVtx,kMaxDiffVtx,
			 kNBDiffVtx,-kMaxDiffVtx,kMaxDiffVtx);
  fArrHisto.AddLast(fHVtxDiffXY);
  fHVtxDiffXY->SetDirectory(0);
  //
  fHVtxDiffZ = new TH1F("vtxMCDiffZ","vtxMC-vtxRec Z",kNBDiffVtx,-kMaxDiffVtx,kMaxDiffVtx);
  fArrHisto.AddLast(fHVtxDiffZ);
  fHVtxDiffZ->SetDirectory(0);
  //
  fHVtxDiffXMlt = new TH2F("VtxDiffXMlt","vX_{MC}-vX_{rec} vs mlt",kNBinMlt,axLogMlt, 
			   kNBDiffVtx,-kMaxDiffVtx,kMaxDiffVtx);
  fArrHisto.AddLast(fHVtxDiffXMlt);
  fHVtxDiffXMlt->SetDirectory(0);
  //
  fHVtxDiffYMlt = new TH2F("VtxDiffYMlt","vY_{MC}-vY_{rec} vs mlt",kNBinMlt,axLogMlt, 
			   kNBDiffVtx,-kMaxDiffVtx,kMaxDiffVtx);
  fArrHisto.AddLast(fHVtxDiffYMlt);
  fHVtxDiffYMlt->SetDirectory(0);
  //
  fHVtxDiffZMlt = new TH2F("VtxDiffZMlt","vZ_{MC}-vZ_{rec} vs mlt",kNBinMlt,axLogMlt, 
			   kNBDiffVtx,-kMaxDiffVtx,kMaxDiffVtx);
  fArrHisto.AddLast(fHVtxDiffZMlt);
  fHVtxDiffZMlt->SetDirectory(0);
  //
  fHVtxPullXMlt = new TH2F("VtxPullXMlt","Pull vX_{MC}-vX_{rec} vs mlt",kNBinMlt,axLogMlt, 
			   kNBPullVtx,-kMaxPullVtx,kMaxPullVtx);
  fArrHisto.AddLast(fHVtxPullXMlt);
  fHVtxPullXMlt->SetDirectory(0);
  //
  fHVtxPullYMlt = new TH2F("VtxPullYMlt","Pull vY_{MC}-vY_{rec} vs mlt",kNBinMlt,axLogMlt, 
			   kNBPullVtx,-kMaxPullVtx,kMaxPullVtx);
  fArrHisto.AddLast(fHVtxPullYMlt);
  fHVtxPullYMlt->SetDirectory(0);
  //
  fHVtxPullZMlt = new TH2F("VtxPullZMlt","Pull vZ_{MC}-vZ_{rec} vs mlt",kNBinMlt,axLogMlt, 
			   kNBPullVtx,-kMaxPullVtx,kMaxPullVtx);
  fArrHisto.AddLast(fHVtxPullZMlt);
  fHVtxPullZMlt->SetDirectory(0);
  //
}

//______________________________________________
void AliITSSAPTracker::SaveHistos(const char* outFName)
{
  // save control histos
  TString fnms = outFName;
  if (fnms.IsNull()) fnms = "XXXITSTrackerControlH.root";
  TFile* fl = TFile::Open(fnms.Data(),"recreate");
  if (!fl) {
    printf("Failed to open output file %s\n",fnms.Data());
    return;
  }
  fArrHisto.Write();
  fl->Close();
  delete fl;
  printf("Stored control histos in %s\n",fnms.Data());
  //
}

//______________________________________________
void AliITSSAPTracker::FillTrackingControlHistos(int lrID,int lbl,const AliExternalTrackParam* track,
						 const double cpar[2],const double ccov[3], 
						 const AliITSRecPoint* bestCl)
{
  // fill control histos for tracking for correct matches
  Bool_t corr = kFALSE;
  for (int i=0;i<3;i++) if (bestCl->GetLabel(i)==lbl) {corr=kTRUE; break;}
  if (!corr) return;
  double pt = track->Pt();
  double dy = cpar[0]-track->GetY();
  double dz = cpar[1]-track->GetZ();
  double sgy = TMath::Sqrt(ccov[0]+track->GetSigmaY2());
  double sgz = TMath::Sqrt(ccov[2]+track->GetSigmaZ2());
  int lrIDS = lrID*10;
  ((TH2F*)fArrHisto[lrIDS+kHResidY])->Fill(pt,dy);
  ((TH2F*)fArrHisto[lrIDS+kHPullY])->Fill(pt,dy/sgy);
  ((TH2F*)fArrHisto[lrIDS+kHResidZ])->Fill(pt,dz);
  ((TH2F*)fArrHisto[lrIDS+kHPullZ])->Fill(pt,dz/sgz);
  ((TH2F*)fArrHisto[lrIDS+kHChi2Cl])->Fill(pt,track->GetPredictedChi2(cpar,ccov));
  //
}

//______________________________________________
Double_t* AliITSSAPTracker::DefLogAx(double xMn,double xMx, int nbin)
{
  // get array for log axis
  if (xMn<=0 || xMx<=xMn || nbin<2) {
    printf("Wrong axis request: %f %f %d\n",xMn,xMx,nbin);
    return 0;
  }
  double dx = log(xMx/xMn)/nbin;
  double *xax = new Double_t[nbin+1];
  for (int i=0;i<=nbin;i++) xax[i]= xMn*exp(dx*i);
  return xax;
}

#endif
//
#ifdef _TIMING_
//______________________________________________
void AliITSSAPTracker::PrintTiming()
{
  // print timing info
  for (int i=0;i<kNSW;i++) {printf("%-10s:\t",fgkSWNames[i]); fSW[i].Print();}
}
#endif


Int_t AliITSSAPTracker::GetTrackPoint( Int_t iLayer, Int_t clusterIndex, AliHLTITSTrackPoint& p ) const 
{
  //--------------------------------------------------------------------
  // Get track space point with index clusterIndex
  //--------------------------------------------------------------------

  const AliITSRecPoint* cl = fLayers[iLayer]->GetClusterSorted( clusterIndex );

  if( !cl ){
    AliError( "AliITSSAPTracker: wrong cluster pointer" );
    p.Reset();
    return -1;
  }

  Int_t idet = cl->GetDetectorIndex();

  bool ok = cl->GetGlobalXYZ(p.fXYZ) && cl->GetGlobalCov(p.fCov);
  if( !ok ){
    AliError( "AliITSSAPTracker: can not get global coordinates of a cluster" );
    p.Reset();
    return -1;
  }

  p.fCharge = cl->GetQ();
  p.fChargeRatio = cl->GetChargeRatio();
  p.fClusterType = cl->GetClusterType();
  p.fDriftTime = cl->GetDriftTime();

  AliGeomManager::ELayerID layerID = AliGeomManager::kInvalidLayer; 
  switch (iLayer) {
  case 0:
    layerID = AliGeomManager::kSPD1;
    break;
  case 1:
    layerID = AliGeomManager::kSPD2;
    break;
  case 2:
    layerID = AliGeomManager::kSDD1;
    break;
  case 3:
    layerID = AliGeomManager::kSDD2;
    break;
  case 4:
    layerID = AliGeomManager::kSSD1;
    break;
  case 5:
    layerID = AliGeomManager::kSSD2;
    break;
  default:
    AliWarning(Form("AliITSSAPTracker: Wrong layer index in ITS (%d) !",iLayer));
    break;
  };
  p.fVolumeID = AliGeomManager::LayerToVolUID(layerID,idet);
  return 0;
}
