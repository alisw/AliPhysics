#include "KMCDetector.h"
#include <TMath.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TF1.h>
#include <TMatrixD.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TFormula.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TText.h>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TStopwatch.h>
#include <TH2.h>
#include <TH1.h>
#include <TArrayI.h>
#include <AliLog.h>

/***********************************************************

Fast Simulation tool for Inner Tracker Systems

***********************************************************/


#define RIDICULOUS 999999 // A ridiculously large resolution (cm) to flag a dead detector

#define Luminosity    1.e27       // Luminosity of the beam (LHC HI == 1.e27, RHIC II == 8.e27 )
#define SigmaD        6.0         // Size of the interaction diamond (cm) (LHC = 6.0 cm)
#define dNdEtaMinB    1//950//660//950           // Multiplicity per unit Eta  (AuAu MinBias = 170, Central = 700)
// #define dNdEtaCent    2300//15000 //1600//2300        // Multiplicity per unit Eta  (LHC at 5.5 TeV not known)

#define CrossSectionMinB         8    // minB Cross section for event under study (PbPb MinBias ~ 8 Barns)
#define AcceptanceOfTpcAndSi     1 //1//0.60 //0.35  // Assumed geometric acceptance (efficiency) of the TPC and Si detectors
#define UPCBackgroundMultiplier  1.0   // Increase multiplicity in detector (0.0 to 1.0 * UPCRate ) (eg 1.0)
#define OtherBackground          0.0   // Increase multiplicity in detector (0.0 to 1.0 * minBias)  (eg 0.0)
#define EfficiencySearchFlag     2     // Define search method:
                                       // -> ChiSquarePlusConfLevel = 2, ChiSquare = 1, Simple = 0.  

#define PionMass                 0.139  // Mass of the Pion
#define KaonMass                 0.498  // Mass of the Kaon
#define D0Mass                   1.865  // Mass of the D0

//TMatrixD *probKomb; // table for efficiency kombinatorics

ClassImp(KMCProbe)

Int_t    KMCProbe::fgNITSLayers = 0;
Double_t KMCProbe::fgMissingHitPenalty = 2.;
//__________________________________________________________________________
KMCProbe& KMCProbe::operator=(const KMCProbe& src) 
{
  if (this!=&src) {
    AliExternalTrackParam::operator=(src);
    fMass = src.fMass;
    fChi2 = src.fChi2;
    fHits = src.fHits;
    fFakes = src.fFakes;
    fNHits = src.fNHits;
    fNHitsITS = src.fNHitsITS;
    fNHitsITSFake = src.fNHitsITSFake;
    fInnLrCheck     = src.fInnLrCheck;
    for (int i=kMaxITSLr;i--;) fClID[i] = src.fClID[i];
  }
  return *this;
}

//__________________________________________________________________________
KMCProbe::KMCProbe(KMCProbe& src) 
  : AliExternalTrackParam(src),
    fMass(src.fMass),
    fChi2(src.fChi2),
    fHits(src.fHits),
    fFakes(src.fFakes),
    fNHits(src.fNHits),
    fNHitsITS(src.fNHitsITS),
    fNHitsITSFake(src.fNHitsITSFake),
    fInnLrCheck(src.fInnLrCheck)
{
  for (int i=kMaxITSLr;i--;) fClID[i] = src.fClID[i];
}

//__________________________________________________________________________
void KMCProbe::ResetCovMat()
{
  // reset errors
  double *trCov  = (double*)GetCovariance();
  double *trPars = (double*)GetParameter();
  const double kLargeErr2Coord = 50*50;
  const double kLargeErr2Dir = 0.6*0.6;
  const double kLargeErr2PtI = 0.5*0.5;
  for (int ic=15;ic--;) trCov[ic] = 0.;
  trCov[kY2]   = trCov[kZ2]   = kLargeErr2Coord; 
  trCov[kSnp2] = trCov[kTgl2] = kLargeErr2Dir;
  trCov[kPtI2] = kLargeErr2PtI*trPars[kPtI]*trPars[kPtI];
  //
}

//__________________________________________________________________________
void KMCProbe::Print(Option_t* option) const
{
  printf("M=%.3f Chi2=%7.2f (Norm:%6.2f|%d) Hits: Total:%d ITS:%d ITSFakes:%d | Y:%+8.4f Z: %+8.4f |", 
	 fMass,fChi2,GetNormChi2(kTRUE),fInnLrCheck,fNHits,fNHitsITS,fNHitsITSFake, GetY(),GetZ());
  for (int i=0;i<fgNITSLayers;i++) {
    if (!IsHit(i)) printf(".");
    else printf("%c",IsHitFake(i) ? '-':'+');
  }
  printf("|%s\n",IsKilled() ? " KILLED":"");
  TString opt = option;
  if (!opt.IsNull()) AliExternalTrackParam::Print(option);
}

//__________________________________________________________________________
Int_t KMCProbe::Compare(const TObject* obj) const
{
  // compare to tracks
  const KMCProbe* trc = (KMCProbe*) obj;
  if (trc->IsKilled()) {
    if (IsKilled()) return 0;
    return -1;
  }
  else if (IsKilled()) return 1;
  double chi2a = GetNormChi2(kTRUE);
  double chi2b = trc->GetNormChi2(kTRUE);
  if (chi2a<chi2b) return  -1;
  if (chi2a>chi2b) return   1;
  return 0;
}

//__________________________________________________________________________
Bool_t KMCProbe::GetXatLabR(Double_t r,Double_t &x, Double_t bz, Int_t dir) const
{
  // Get local X of the track position estimated at the radius lab radius r. 
  // The track curvature is accounted exactly
  //
  // The flag "dir" can be used to remove the ambiguity of which intersection to take (out of 2 possible)
  // 0  - take the intersection closest to the current track position
  // >0 - go along the track (increasing fX)
  // <0 - go backward (decreasing fX)
  //
  // special case of R=0
  if (r<kAlmost0) {x=0; return kTRUE;}

  const double* pars = GetParameter();
  const Double_t &fy=pars[0], &sn = pars[2];
  //
  double fx = GetX();
  double crv = GetC(bz);
  if (TMath::Abs(crv)<=kAlmost0) { // this is a straight track
    if (TMath::Abs(sn)>=kAlmost1) { // || to Y axis
      double det = (r-fx)*(r+fx);
      if (det<0) return kFALSE;     // does not reach raduis r
      x = fx;
      if (dir==0) return kTRUE;
      det = TMath::Sqrt(det);
      if (dir>0) {                       // along the track direction
	if (sn>0) {if (fy>det)  return kFALSE;} // track is along Y axis and above the circle
	else      {if (fy<-det) return kFALSE;} // track is against Y axis amd belo the circle
      }
      else if(dir>0) {                                    // agains track direction
	if (sn>0) {if (fy<-det) return kFALSE;} // track is along Y axis
        else if (fy>det)  return kFALSE;        // track is against Y axis
      }
    }
    else if (TMath::Abs(sn)<=kAlmost0) { // || to X axis
      double det = (r-fy)*(r+fy);
      if (det<0) return kFALSE;     // does not reach raduis r
      det = TMath::Sqrt(det);
      if (!dir) {
	x = fx>0  ? det : -det;    // choose the solution requiring the smalest step
	return kTRUE;
      }
      else if (dir>0) {                    // along the track direction
	if      (fx > det) return kFALSE;  // current point is in on the right from the circle
	else if (fx <-det) x = -det;       // on the left
	else               x =  det;       // within the circle
      }
      else {                               // against the track direction
	if      (fx <-det) return kFALSE;  
	else if (fx > det) x =  det;
	else               x = -det;
      }
    }
    else {                                 // general case of straight line
      double cs = TMath::Sqrt((1-sn)*(1+sn));
      double xsyc = fx*sn-fy*cs;
      double det = (r-xsyc)*(r+xsyc);
      if (det<0) return kFALSE;    // does not reach raduis r
      det = TMath::Sqrt(det);
      double xcys = fx*cs+fy*sn;
      double t = -xcys;
      if (dir==0) t += t>0 ? -det:det;  // chose the solution requiring the smalest step
      else if (dir>0) {                 // go in increasing fX direction. ( t+-det > 0)
	if (t>=-det) t += -det;         // take minimal step giving t>0
	else return kFALSE;             // both solutions have negative t
      }
      else {                            // go in increasing fx direction. (t+-det < 0)
	if (t<det) t -= det;            // take minimal step giving t<0
	else return kFALSE;             // both solutions have positive t
      }
      x = fx + cs*t;
    }
  }
  else {                                 // helix
    // get center of the track circle
    double tR = 1./crv;   // track radius (for the moment signed)
    double cs = TMath::Sqrt((1-sn)*(1+sn));
    double x0 = fx - sn*tR;
    double y0 = fy + cs*tR;
    double r0 = TMath::Sqrt(x0*x0+y0*y0);
    //    printf("Xc:%+e Yc:%+e\n",x0,y0);
    //
    if (r0<=kAlmost0) {
      AliDebug(2,Form("r0 = %f",r0));
      return kFALSE;
    }            // the track is concentric to circle
    tR = TMath::Abs(tR);
    double tR2r0 = tR/r0;
    double g = 0.5*(r*r/(r0*tR) - tR2r0 - 1./tR2r0);
    double det = (1.-g)*(1.+g);
    if (det<0) {
      AliDebug(2,Form("g=%f tR=%f r0=%f\n",g,tR, r0));
      return kFALSE;
    }         // does not reach raduis r
    det = TMath::Sqrt(det);
    //
    // the intersection happens in 2 points: {x0+tR*C,y0+tR*S} 
    // with C=f*c0+-|s0|*det and S=f*s0-+c0 sign(s0)*det
    // where s0 and c0 make direction for the circle center (=x0/r0 and y0/r0)
    //
    double tmp = 1.+g*tR2r0;
    x = x0*tmp; 
    double y = y0*tmp;
    if (TMath::Abs(y0)>kAlmost0) { // when y0==0 the x,y is unique
      double dfx = tR2r0*TMath::Abs(y0)*det;
      double dfy = tR2r0*x0*TMath::Sign(det,y0);
      if (dir==0) {                    // chose the one which corresponds to smallest step 
	double delta = (x-fx)*dfx-(y-fy)*dfy; // the choice of + in C will lead to smaller step if delta<0
	if (delta<0) x += dfx;
	else         x -= dfx;
      }
      else if (dir>0) {  // along track direction: x must be > fx
	x -= dfx; // try the smallest step (dfx is positive)
	if (x<fx && (x+=dfx+dfx)<fx) return kFALSE;
      }
      else { // backward: x must be < fx
	x += dfx; // try the smallest step (dfx is positive)
	if (x>fx && (x-=dfx+dfx)>fx) return kFALSE;
      }
    }
    else { // special case: track touching the circle just in 1 point
      if ( (dir>0&&x<fx) || (dir<0&&x>fx) ) return kFALSE; 
    }
  }
  //
  return kTRUE;
}

//____________________________________
Bool_t KMCProbe::PropagateToR(double r, double b, int dir) 
{
  // go to radius R
  //
  double xR = 0;
  double rr = r*r;
  int iter = 0;
  const double kTiny = 1e-4;
  while(1) {
    if (!GetXatLabR(r ,xR, b, dir)) {
      //      printf("Track with pt=%f cannot reach radius %f\n",Pt(),r);
      //      Print("l");
      return kFALSE;
    }
    
    if (!PropagateTo(xR, b)) {
      if (AliLog::GetGlobalDebugLevel()>2) {
	printf("Failed to propagate to X=%f for R=%f\n",xR,r); 
	Print("l"); 
      }
      return kFALSE;
    }
    double rcurr2 = xR*xR + GetY()*GetY();
    if (TMath::Abs(rcurr2-rr)<kTiny || rr<kAlmost0) return kTRUE;
    //
    // two radii correspond to this X...
    double pos[3]; GetXYZ(pos);
    double phi = TMath::ATan2(pos[1],pos[0]);
    if (!Rotate(phi)) {
      if (AliLog::GetGlobalDebugLevel()>2) {
	printf("Failed to rotate to %f to propagate to R=%f\n",phi,r); 
	Print("l"); 
      }
      return kFALSE;
    }
    if (++iter>8) {
      if (AliLog::GetGlobalDebugLevel()>2) {
	printf("Failed to propagate to R=%f after %d steps\n",r,iter); 
	Print("l"); 
      }
      return kFALSE;
    }
  } 
  return kTRUE;
}


//__________________________________________________________________________
Bool_t KMCProbe::CorrectForMeanMaterial(const KMCLayer* lr, Bool_t inward)
{
  //  printf("before at r=%.1f p=%.4f\n",lr->fR, P());
  if (AliExternalTrackParam::CorrectForMeanMaterial(lr->fx2X0, inward ? lr->fXRho : -lr->fXRho, GetMass() , kTRUE)) {
    //  printf("after  at r=%.1f p=%.4f\n",lr->fR, P());
    return kTRUE;
  }
  AliDebug(2,Form("Failed to apply material correction, X/X0=%.4f", lr->fx2X0));
  if (AliLog::GetGlobalDebugLevel()>1) Print();
  return kFALSE;
}

/////////////////////////////////////////////////////////////////////////////
ClassImp(KMCCluster)

//_________________________________________________________________________
KMCCluster::KMCCluster(KMCCluster &src) 
: TObject(src),
  fY(src.fY),fZ(src.fZ),fX(src.fX),fPhi(src.fPhi)
{}

//__________________________________________________________________________
KMCCluster& KMCCluster::operator=(const KMCCluster& src) 
{
  if (this!=&src) {
    TObject::operator=(src);
    fY = src.fY;
    fZ = src.fZ;
    fX = src.fX;
    fPhi = src.fPhi;
  }
  return *this;
}

//_________________________________________________________________________
void KMCCluster::Print(Option_t *) const 
{
  printf(" Local YZ = (%3.4lf,%3.4lf) | X=%3.4lf  phi: %+.3f %s\n",fY,fZ,fX,fPhi,IsKilled()?"Killed":""); 
}

/////////////////////////////////////////////////////////////////////////////
ClassImp(KMCLayer)

Double_t KMCLayer::fgDefEff = 1.0;
//__________________________________________________________________________
KMCLayer::KMCLayer(char *name) : 
  TNamed(name,name),fR(0),fx2X0(0),fPhiRes(0),fZRes(0),fEff(0),fIsDead(kFALSE),fType(-1),fActiveID(-1),fSig2EstD(999),fSig2EstZ(999),
  fClCorr(),fClMC(),fClBg("KMCCluster",5), fTrCorr(), fTrMC("KMCProbe",5)
{
  Reset();
}

//__________________________________________________________________________
void KMCLayer::Reset() 
{
  fTrCorr.Reset();
  fClCorr.Reset();
  ResetMC();
  fSig2EstD = fSig2EstZ = 999;
  //
}

//__________________________________________________________________________
KMCProbe* KMCLayer::AddMCTrack(KMCProbe* src) 
{
  int ntr = GetNMCTracks(); 
  KMCProbe* prb = 0;
  if (src) prb = new(fTrMC[ntr]) KMCProbe(*src);
  else     prb = new(fTrMC[ntr]) KMCProbe();
  if (!IsDead()) prb->ResetHit(GetActiveID());
  return prb;
}

//__________________________________________________________________________
void KMCLayer::Print(Option_t *opt) const
{
  printf("Lr%3d(A%3d) %10s R=%5.1f X2X0=%.3f XRho=%.3f SigY=%.4f SigZ=%.4f Eff:%4.2f\n",GetUniqueID(),fActiveID,GetName(), fR, fx2X0,fXRho,fPhiRes,fZRes,fEff);
  TString opts = opt; opts.ToLower();
  if (opts.Contains("c")) {
    printf("Clusters: MC: %+7.4f:%+7.4f Ideal: %+7.4f:%+7.4f  NBgCl: %3d NTrMC: %4d\n",fClMC.fY,fClMC.fZ, fClCorr.fY,fClCorr.fZ, GetNBgClusters(),GetNMCTracks());
  }
}

/////////////////////////////////////////////////////////////////////////////
Double_t KMCDetector::fgVtxConstraint[2]={-1,-1};

ClassImp(KMCDetector)
KMCDetector::KMCDetector() :
TNamed("test_detector","detector"),
  fNLayers(0),
  fNActiveLayers(0),
  fNActiveITSLayers(0),
  fLastActiveLayer(-1),
  fLastActiveITSLayer(-1),
  fLastActiveLayerTracked(-1),
  fBFieldG(5.),
  fLhcUPCscale(1.0),    
  fIntegrationTime(0.02), // in ms
  fConfLevel(0.0027),      // 0.27 % -> 3 sigma confidence
  fAvgRapidity(0.45),      // Avg rapidity, MCS calc is a function of crossing angle
  fParticleMass(0.140),    // Standard: pion mass 
  fMaxRadiusSlowDet(10.),
  fAtLeastCorr(-1),     // if -1, then correct hit on all ITS layers
  fAtLeastFake(1),       // if at least x fakes, track is considered fake ...
  fdNdEtaCent(2300),
  fMaxChi2Cl(25.),
  fMaxNormChi2NDF(5.),
  fMinITSHits(4),
//
  fMaxChi2ClSQ(4.), // precalulated internally
  fMaxSeedToPropagate(300),
//
  fUseBackground(kFALSE),
  fBgYMin(1e6),fBgYMax(-1e6),fBgZMin(1e6),fBgZMax(-1e6),
  fBgYMinTr(100),fBgYMaxTr(100),fBgZMinTr(100),fBgZMaxTr(100),fNBgLimits(0),
  fDensFactorEta(1.),
//
  fUpdCalls(0),
  fHMCLrResidRPhi(0),
  fHMCLrResidZ(0),
  fHMCLrChi2(0),
  fPattITS(0)
{
  //
  // default constructor
  //
  //  fLayers = new TObjArray();
  RequireMaxChi2Cl(fMaxChi2Cl); // just to precalulate default square
}

KMCDetector::KMCDetector(char *name, char *title)
  : TNamed(name,title),
    fNLayers(0),
    fNActiveLayers(0),
    fNActiveITSLayers(0),
    fLastActiveLayer(-1),
    fLastActiveITSLayer(-1),    
    fLastActiveLayerTracked(-1),
    fBFieldG(5.0),
    fLhcUPCscale(1.0),
    fIntegrationTime(0.02),  // in ms
    fConfLevel(0.0027),      // 0.27 % -> 3 sigma confidence
    fAvgRapidity(0.45),      // Avg rapidity, MCS calc is a function of crossing angle
    fParticleMass(0.140),     // Standard: pion mass
    fMaxRadiusSlowDet(10.),
    fAtLeastCorr(-1),     // if -1, then correct hit on all ITS layers
    fAtLeastFake(1),       // if at least x fakes, track is considered fake ...
    fdNdEtaCent(2300),
    fMaxChi2Cl(9.),
    fMaxNormChi2NDF(5.),
    fMinITSHits(4),
    //
    fMaxChi2ClSQ(3.), // precalulated internally
    fMaxSeedToPropagate(50),
    //
    fUseBackground(kFALSE),
    fBgYMin(1e6),fBgYMax(-1e6),fBgZMin(1e6),fBgZMax(-1e6),
    fBgYMinTr(100),fBgYMaxTr(100),fBgZMinTr(100),fBgZMaxTr(100),fNBgLimits(0),
    fDensFactorEta(1.),
    //
    fUpdCalls(0),
    fHMCLrResidRPhi(0),
    fHMCLrResidZ(0),
    fHMCLrChi2(0),
    fPattITS(0)
{
  //
  // default constructor, that set the name and title
  //
  //  fLayers = new TObjArray();
}
KMCDetector::~KMCDetector() { // 
  // virtual destructor
  //
  //  delete fLayers;
}

void KMCDetector::InitMCWatchHistos()
{
  // init utility histos used for MC tuning
  enum {kNBinRes=1000};
  const double kMaxResidRPhi=1.0,kMaxResidZ=1.0,kMaxChi2=50;
  int nlr = fNActiveITSLayers;
  TString nam = "mc_residrhi";
  TString tit = "MC $Delta Cl-Tr R#phi";
  fHMCLrResidRPhi = new TH2F(nam.Data(),tit.Data(),kNBinRes,-kMaxResidRPhi,kMaxResidRPhi,nlr,-0.5,nlr-0.5);
  fHMCLrResidRPhi->GetXaxis()->SetTitle("cl-tr #Delta r#phi");
  fHMCLrResidRPhi->Sumw2();
    //
  nam = "mc_residz";
  tit = "MC $Delta Cl-Tr Z";
  fHMCLrResidZ = new TH2F(nam.Data(),tit.Data(),kNBinRes,-kMaxResidRPhi,kMaxResidZ,nlr,-0.5,nlr-0.5);
  fHMCLrResidZ->GetXaxis()->SetTitle("cl-tr #Delta z");
  fHMCLrResidZ->Sumw2();
    //
  nam = "mc_chi2";
  tit = "MC #chi^{2} Cl-Tr Z";
  fHMCLrChi2 = new TH2F(nam.Data(),tit.Data(),kNBinRes,-kMaxResidRPhi,kMaxChi2,nlr,-0.5,nlr-0.5);
  fHMCLrChi2->GetXaxis()->SetTitle("cl-tr #chi^{2}");
  fHMCLrChi2->Sumw2();
  //
  SetBit(kUtilHisto);
}


void KMCDetector::AddLayer(char *name, Float_t radius, Float_t x2X0, Float_t xrho, Float_t phiRes, Float_t zRes, Float_t eff) {
  //
  // Add additional layer to the list of layers (ordered by radius)
  // 

  KMCLayer *newLayer = (KMCLayer*) fLayers.FindObject(name);

  if (!newLayer) {
    newLayer = new KMCLayer(name);
    newLayer->fR = radius;
    newLayer->fx2X0 = x2X0;
    newLayer->fXRho  = xrho;
    newLayer->fPhiRes = phiRes;
    newLayer->fZRes = zRes;
    eff = TMath::Min(1.f,eff);
    newLayer->fEff = eff <0 ? KMCLayer::GetDefEff() : eff;
    newLayer->fActiveID = -2;
    TString lname = name;
    newLayer->fType = KMCLayer::kTypeNA;
    if      (lname.Contains("tpc")) newLayer->fType = KMCLayer::kTPC;
    else if (lname.Contains("its")) newLayer->fType = KMCLayer::kITS;
    if (lname.Contains("vertex")) newLayer->SetBit(KMCLayer::kBitVertex);
    //
    if (newLayer->fType==KMCLayer::kTypeNA) printf("Attention: the layer %s has undefined type\n",name);
    //
    newLayer->fIsDead =  (newLayer->fPhiRes==RIDICULOUS && newLayer->fZRes==RIDICULOUS);
    //
    if (fLayers.GetEntries()==0) 
      fLayers.Add(newLayer);
    else {
      //
      for (Int_t i = 0; i<fLayers.GetEntries(); i++) {
	KMCLayer *l = (KMCLayer*)fLayers.At(i);
	if (radius<l->fR) { fLayers.AddBefore(l,newLayer); break; }
	  if (radius>l->fR && (i+1)==fLayers.GetEntries() ) fLayers.Add(newLayer); // even bigger then last one
      }
      //
    }
    //
    ClassifyLayers();
    //
  } else {
    printf("Layer with the name %s does already exist\n",name);
  }
}

//____________________________________________________________
void KMCDetector::ClassifyLayers()
{
  // assign active Id's, etc
  fLastActiveLayer = -1;
  fLastActiveITSLayer = -1;
  fNActiveLayers = 0;
  fNActiveITSLayers = 0;
  //
  int nl = GetNLayers();
  for (int il=0;il<nl;il++) {
    KMCLayer* lr = GetLayer(il);
    lr->SetUniqueID(il);
    if (!lr->IsDead()) {
      fLastActiveLayer = il; 
      lr->fActiveID = fNActiveLayers++;
      if (lr->IsITS()) {
	fLastActiveITSLayer = il;
	fNActiveITSLayers++;
      }
    }
  }
  //
  KMCProbe::SetNITSLayers(fNActiveITSLayers);
}

//____________________________________________________________
void KMCDetector::KillLayer(char *name) {
  //
  // Marks layer as dead. Contribution only by Material Budget
  //

  KMCLayer *tmp = (KMCLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot mark as dead\n",name);
  else {
     tmp->fPhiRes = 999999;
     tmp->fZRes = 999999;
     tmp->fIsDead = kTRUE;
     ClassifyLayers();
  }
}

void KMCDetector::SetRadius(char *name, Float_t radius) {
  //
  // Set layer radius [cm]
  //
  KMCLayer *tmp = (KMCLayer*) fLayers.FindObject(name);
  //
  if (!tmp) {
    printf("Layer %s not found - cannot set radius\n",name);
  } else {
      
    Float_t tmpRadL  = tmp->fx2X0;
    Float_t tmpPhiRes = tmp->fPhiRes;
    Float_t tmpZRes = tmp->fZRes;
    Float_t tmpXRho = tmp->fXRho;
    RemoveLayer(name); // so that the ordering is correct
    AddLayer(name,radius,tmpRadL,tmpXRho,tmpPhiRes,tmpZRes);
  }
}

Float_t KMCDetector::GetRadius(char *name) {
  //
  // Return layer radius [cm]
  //

  KMCLayer *tmp = (KMCLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot get radius\n",name);
  else 
    return tmp->fR;

  return 0;
}

void KMCDetector::SetRadiationLength(char *name, Float_t x2X0) {
  //
  // Set layer material [cm]
  //

  KMCLayer *tmp = (KMCLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot set layer material\n",name);
  else {
    tmp->fx2X0 = x2X0;
  }
}

Float_t KMCDetector::GetRadiationLength(char *name) {
  //
  // Return layer radius [cm]
  //

  KMCLayer *tmp = (KMCLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot get layer material\n",name);
  else 
    return tmp->fx2X0;
    
  return 0;
  
}

void KMCDetector::SetResolution(char *name, Float_t phiRes, Float_t zRes) {
  //
  // Set layer resolution in [cm]
  //

  KMCLayer *tmp = (KMCLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot set resolution\n",name);
  else {
    tmp->fPhiRes = phiRes;
    tmp->fZRes = zRes;
    tmp->fIsDead = (zRes==RIDICULOUS && phiRes==RIDICULOUS);
    ClassifyLayers();
  }
}

Float_t KMCDetector::GetResolution(char *name, Int_t axis) {
  //
  // Return layer resolution in [cm]
  // axis = 0: resolution in rphi
  // axis = 1: resolution in z
  //

  KMCLayer *tmp = GetLayer(name);
  if (!tmp) 
    printf("Layer %s not found - cannot get resolution\n",name);
  else {
    if (axis==0) return tmp->fPhiRes;
    if (axis==1) return tmp->fZRes;
    printf("error: axis must be either 0 or 1 (rphi or z axis)\n");
  }
  return 0;
}

void KMCDetector::SetLayerEfficiency(char *name, Float_t eff) {
  //
  // Set layer efficnecy (prop that his is missed within this layer) 
  //

  KMCLayer *tmp = (KMCLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot set layer efficiency\n",name);
  else {
    tmp->fEff = eff;
  }
}

Float_t KMCDetector::GetLayerEfficiency(char *name) {
  //
  // Get layer efficnecy (prop that his is missed within this layer) 
  //

  KMCLayer *tmp = (KMCLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot get layer efficneicy\n",name);
  else 
    return tmp->fEff;
    
  return 0;
  
}

void KMCDetector::RemoveLayer(char *name) {
  //
  // Removes a layer from the list
  //

  KMCLayer *tmp = (KMCLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot remove it\n",name);
  else {
    fLayers.Remove(tmp);
    ClassifyLayers();
  }
}


void KMCDetector::PrintLayout() {
  //
  // Prints the detector layout
  //

  printf("Detector %s: \"%s\"\n",GetName(),GetTitle());
  
  if (fLayers.GetEntries()>0) 
    printf("  Name \t\t r [cm] \t  X0 \t  phi & z res [um]\n");

  KMCLayer *tmp = 0;
  for (Int_t i = 0; i<fLayers.GetEntries(); i++) {
    tmp = (KMCLayer*)fLayers.At(i);
  
    // don't print all the tpc layers
    TString name(tmp->GetName());
    if (name.Contains("tpc") && (!name.Contains("tpc_0")) ) continue;

    printf("%d. %s \t %03.2f   \t%1.4f\t  ",i,
 	   tmp->GetName(), tmp->fR, tmp->fx2X0);
    if (tmp->fPhiRes==RIDICULOUS) 
      printf("  -  ");
    else
      printf("%3.0f   ",tmp->fPhiRes*10000);
    if (tmp->fZRes==RIDICULOUS) 
      printf("  -\n");
    else
      printf("%3.0f\n",tmp->fZRes*10000);
  }
}

void KMCDetector::PlotLayout(Int_t plotDead) {
  //
  // Plots the detector layout in Front view
  //

  Double_t x0=0, y0=0;

  TGraphErrors *gr = new TGraphErrors();
  gr->SetPoint(0,0,0);
  KMCLayer *lastLayer = (KMCLayer*)fLayers.At(fLayers.GetEntries()-1);  Double_t maxRad = lastLayer->fR;
  gr->SetPointError(0,maxRad,maxRad);
  gr->Draw("APE");
  

  KMCLayer *tmp = 0;
  for (Int_t i = fLayers.GetEntries()-1; i>=0; i--) {
    tmp = (KMCLayer*)fLayers.At(i);
  

    Double_t txtpos = tmp->fR;
    if ((tmp->IsDead())) txtpos*=-1; //
    TText *txt = new TText(x0,txtpos,tmp->GetName());
    txt->SetTextSizePixels(5); txt->SetTextAlign(21);
    if (!tmp->IsDead() || plotDead) txt->Draw();

    TEllipse *layEl = new TEllipse(x0,y0,tmp->fR);
    //  layEl->SetFillColor(5);
    layEl->SetFillStyle(5001);
    layEl->SetLineStyle(tmp->IsDead()+1); // dashed if not active
    layEl->SetLineColor(4);
    TString name(tmp->GetName());
    if (!tmp->IsDead()) layEl->SetLineWidth(2);
    if (name.Contains("tpc") )  layEl->SetLineColor(29);

    if (!tmp->IsDead() || plotDead) layEl->Draw();
  
  }

}



void KMCDetector::AddTPC(Float_t phiResMean, Float_t zResMean, Int_t skip) {
  //
  // Emulates the TPC
  // 
  // skip=1: Use every padrow, skip=2: Signal in every 2nd padrow 


  AddLayer((char*)"IFCtpc",   77.8,0.01367, 0); // Inner Field cage (RS: set correct xrho for eloss)
  
  // % Radiation Lengths ... Average per TPC row  (i.e. total/159 )
  Float_t x2X0PerRow = 0.000036;
  
  Float_t tpcInnerRadialPitch  =    0.75 ;    // cm
  Float_t tpcMiddleRadialPitch =    1.0  ;    // cm
  Float_t tpcOuterRadialPitch  =    1.5  ;    // cm
  //  Float_t tpcInnerPadWidth     =    0.4  ;    // cm
  //  Float_t tpcMiddlePadWidth    =    0.6   ;   // cm
  //  Float_t tpcOuterPadWidth     =    0.6   ;   // cm
  Float_t innerRows            =   63 ;
  Float_t middleRows           =   64  ;
  Float_t outerRows            =   32  ;
  Float_t tpcRows            =   (innerRows + middleRows + outerRows) ;
  Float_t rowOneRadius         =   85.2  ;    // cm
  Float_t row64Radius          =  135.1  ;    // cm
  Float_t row128Radius         =  199.2  ;    // cm                       
 
  for ( Int_t k = 0 ; k < tpcRows ; k++ ) {
    
    Float_t rowRadius =0;
    if (k<innerRows) 
      rowRadius =  rowOneRadius + k*tpcInnerRadialPitch ;
    else if ( k>=innerRows && k<(innerRows+middleRows) )
      rowRadius =  row64Radius + (k-innerRows+1)*tpcMiddleRadialPitch ;
    else if (k>=(innerRows+middleRows) && k<tpcRows )
      rowRadius = row128Radius + (k-innerRows-middleRows+1)*tpcOuterRadialPitch ;

    if ( k%skip == 0 )
      AddLayer(Form("tpc_%d",k),rowRadius,x2X0PerRow,0, phiResMean,zResMean);    
    else 
      AddLayer(Form("tpc_%d",k),rowRadius,x2X0PerRow,0); // non "active" row
    
  
  }
 
}

void KMCDetector::RemoveTPC() {

  // flag as dead, although resolution is ok ... makes live easier in the prints ... ;-)
  KMCLayer *tmp = 0;
  for (Int_t i = 0; i<fLayers.GetEntries(); i++) {
    tmp = (KMCLayer*)fLayers.At(i);  
    TString name(tmp->GetName());
    if (name.Contains("tpc")) { RemoveLayer((char*)name.Data()); i--; }
  }
  RemoveLayer((char*)"IFC");
  
}


Double_t KMCDetector::ThetaMCS ( Double_t mass, Double_t x2X0, Double_t momentum ) const
{
  //
  // returns the Multiple Couloumb scattering angle (compare PDG boolet, 2010, equ. 27.14)
  //

  Double_t beta  =  momentum / TMath::Sqrt(momentum*momentum+mass*mass)  ;
  Double_t theta =  0.0 ;    // Momentum and mass in GeV
  // if ( RadLength > 0 ) theta  =  0.0136 * TMath::Sqrt(RadLength) / ( beta * momentum );
  if ( x2X0 > 0 ) theta  =  0.0136 * TMath::Sqrt(x2X0) / ( beta * momentum ) * (1+0.038*TMath::Log(x2X0)) ;
  return (theta) ;
}


Double_t KMCDetector::ProbGoodHit ( Double_t radius, Double_t searchRadiusRPhi, Double_t searchRadiusZ ) 
{
  // Based on work by Howard Wieman: http://rnc.lbl.gov/~wieman/GhostTracks.htm 
  // and http://rnc.lbl.gov/~wieman/HitFinding2D.htm
  // This is the probability of getting a good hit using 2D Gaussian distribution function and infinite search radius
  Double_t sx, sy, goodHit ;
  sx = 2 * TMath::Pi() *  searchRadiusRPhi * searchRadiusRPhi * HitDensity(radius) ;
  sy = 2 * TMath::Pi() *  searchRadiusZ    * searchRadiusZ    * HitDensity(radius) ;
  goodHit =  TMath::Sqrt(1./((1+sx)*(1+sy)))  ;
  return ( goodHit ) ;
}


Double_t KMCDetector::ProbGoodChiSqHit ( Double_t radius, Double_t searchRadiusRPhi, Double_t searchRadiusZ ) 
{
  // Based on work by Victor Perevoztchikov and Howard Wieman: http://rnc.lbl.gov/~wieman/HitFinding2DXsq.htm
  // This is the probability of getting a good hit using a Chi**2 search on a 2D Gaussian distribution function
  Double_t sx, goodHit ;
  sx = 2 * TMath::Pi() *  searchRadiusRPhi * searchRadiusZ * HitDensity(radius) ;
  goodHit =  1./(1+sx) ;
  return ( goodHit ) ;  
}

Double_t KMCDetector::ProbGoodChiSqPlusConfHit ( Double_t radius, Double_t leff, Double_t searchRadiusRPhi, Double_t searchRadiusZ ) 
{
  // Based on work by Ruben Shahoyen 
  // This is the probability of getting a good hit using a Chi**2 search on a 2D Gaussian distribution function
  // Plus, in addition, taking a "confidence level" and the "layer efficiency" into account 
  // Following is correct for 2 DOF

  Double_t c = -2 *TMath::Log(fConfLevel); // quantile at cut of confidence level
  Double_t alpha = (1 + 2 * TMath::Pi() * HitDensity(radius) * searchRadiusRPhi * searchRadiusZ)/2; 
  Double_t goodHit = leff/(2*alpha) * (1 - TMath::Exp(-alpha*c));
  return ( goodHit ) ;  
}

Double_t KMCDetector::ProbNullChiSqPlusConfHit ( Double_t radius, Double_t leff, Double_t searchRadiusRPhi, Double_t searchRadiusZ ) 
{
  // Based on work by Ruben Shahoyen 
  // This is the probability to not have any match to the track (see also :ProbGoodChiSqPlusConfHit:)

  Double_t c = -2 *TMath::Log(fConfLevel); // quantile at cut of confidence level
  Double_t alpha = (1 + 2 * TMath::Pi() * HitDensity(radius) * searchRadiusRPhi * searchRadiusZ)/2; 
  Double_t nullHit = (1-leff+fConfLevel*leff)*TMath::Exp(-c*(alpha-1./2));
  return ( nullHit ) ;  
}

Double_t KMCDetector::HitDensity ( Double_t radius ) 
{
  // Background (0-1) is included via 'OtherBackground' which multiplies the minBias rate by a scale factor.
  // UPC electrons is a temporary kludge that is based on Kai Schweda's summary of Kai Hainken's MC results
  // See K. Hencken et al. PRC 69, 054902 (2004) and PPT slides by Kai Schweda.
  // Note that this function assumes we are working in CM and CM**2 [not meters].
  // Based on work by Yan Lu 12/20/2006, all radii and densities in centimeters or cm**2.

  //  Double_t MaxRadiusSlowDet = 0.1; //?   // Maximum radius for slow detectors.  Fast detectors 
                                        // and only fast detectors reside outside this radius.
  Double_t arealDensity = 0 ;
  if (radius<0.01) return 0;

  if ( radius >= fMaxRadiusSlowDet ) 
    {
      arealDensity  = OneEventHitDensity(fdNdEtaCent,radius)  ; // Fast detectors see central collision density (only)
      arealDensity += OtherBackground*OneEventHitDensity(dNdEtaMinB,radius)  ;  // Increase density due to background 
    }

  if (radius < fMaxRadiusSlowDet )
    { // Note that IntegratedHitDensity will always be minB one event, or more, even if integration time => zero.
      arealDensity  = OneEventHitDensity(fdNdEtaCent,radius) 
	            + IntegratedHitDensity(dNdEtaMinB,radius) 
	            + UpcHitDensity(radius) ;
      arealDensity += OtherBackground*IntegratedHitDensity(dNdEtaMinB,radius) ;  
      // Increase density due to background 
    } 

  return ( arealDensity ) ;  
}


double KMCDetector::OneEventHitDensity( Double_t multiplicity, Double_t radius ) const
{
  // This is for one event at the vertex.  No smearing.
  double den   = multiplicity / (2.*TMath::Pi()*radius*radius) * fDensFactorEta ; // 2 eta ?
  // note: surface of sphere is  '4*pi*r^2'
  //       surface of cylinder is '2*pi*r* h' 
  return den ;
} 


double KMCDetector::IntegratedHitDensity(Double_t multiplicity, Double_t radius)
{ 
  // The integral of minBias events smeared over a gaussian vertex distribution.
  // Based on work by Yan Lu 12/20/2006, all radii in centimeters.

  Double_t zdcHz = Luminosity * 1.e-24 * CrossSectionMinB ;
  Double_t den   = zdcHz * fIntegrationTime/1000. * multiplicity * Dist(0., radius) / (2.*TMath::Pi()*radius) ;

  // Note that we do not allow the rate*time calculation to fall below one minB event at the vertex.
  double dens1 = OneEventHitDensity(multiplicity,radius);
  if ( den < dens1 )  den = dens1;

  return den ;
} 


double KMCDetector::UpcHitDensity(Double_t radius)
{ 
  // QED electrons ...

  Double_t mUPCelectrons ;                                 ;  
  //  mUPCelectrons =  fLhcUPCscale * (1.23 - radius/6.5)      ;  // Fit to Kai Schweda summary tables at RHIC * 'scale' for LHC
  mUPCelectrons = fLhcUPCscale*5456/(radius*radius)/dNdEtaMinB;      // Fit to 'Rossegger,Sadovsky'-Alice simulation
  if ( mUPCelectrons < 0 ) mUPCelectrons =  0.0             ;  // UPC electrons fall off quickly and don't go to large R
  mUPCelectrons *= IntegratedHitDensity(dNdEtaMinB,radius) ;  // UPCs increase Mulitiplicty ~ proportional to MinBias rate
  mUPCelectrons *= UPCBackgroundMultiplier                 ;  // Allow for an external multiplier (eg 0-1) to turn off UPC

  return mUPCelectrons ;
} 


double KMCDetector::Dist(double z, double r)
{
  // Convolute dEta/dZ  distribution with assumed Gaussian of vertex z distribution
  // Based on work by Howard Wieman http://rnc.lbl.gov/~wieman/HitDensityMeasuredLuminosity7.htm
  // Based on work by Yan Lu 12/20/2006, all radii and Z location in centimeters.
  Int_t    index  =  1     ;     // Start weight at 1 for Simpsons rule integration
  Int_t    nsteps =  301   ;     // NSteps must be odd for Simpson's rule to work
  double   dist   =  0.0   ;
  double   dz0    =  ( 4*SigmaD - (-4)*SigmaD ) / (nsteps-1)  ;  //cm
  double    z0    =  0.0   ;     //cm
  for(int i=0; i<nsteps; i++){
    if ( i == nsteps-1 ) index = 1 ;
    z0 = -4*SigmaD + i*dz0 ;
    dist += index * (dz0/3.) * (1/sqrt(2.*TMath::Pi())/SigmaD) * exp(-z0*z0/2./SigmaD/SigmaD) * 
      (1/sqrt((z-z0)*(z-z0) + r*r)) ;
    if ( index != 4 ) index = 4; else index = 2 ;
  }
  return dist; 
}

#define  PZero   0.861  // Momentum of back to back decay particles in the CM frame
#define  EPiZero 0.872  // Energy of the pion from a D0 decay at rest
#define  EKZero  0.993  // Energy of the Kaon from a D0 decay at rest

Double_t KMCDetector::D0IntegratedEfficiency( Double_t pt, Double_t corrEfficiency[][20] ) const {
  // Math from Ron Longacre.  Note hardwired energy to bin conversion for PtK and PtPi.

  Double_t const1  =  pt / D0Mass ;
  Double_t const2  =  TMath::Sqrt(pt*pt+D0Mass*D0Mass) / D0Mass ;
  Double_t sum, ptPi, ptK ;
  Double_t effp, effk ;

  sum = 0.0 ;
  for ( Int_t k = 0 ; k < 360 ; k++ )   {
    
    Double_t theta = k * TMath::Pi() / 180. ;
    
    ptPi = TMath::Sqrt( 
		       PZero*PZero*TMath::Cos(theta)*TMath::Cos(theta)*const2*const2 +
		       const1*const1*EPiZero*EPiZero -
		       2*PZero*TMath::Cos(theta)*const2*const1*EPiZero +
		       PZero*PZero*TMath::Sin(theta)*TMath::Sin(theta)
		       ) ;
    
    ptK = TMath::Sqrt( 
		      PZero*PZero*TMath::Cos(theta)*TMath::Cos(theta)*const2*const2 +
		      const1*const1*EKZero*EKZero +
		      2*PZero*TMath::Cos(theta)*const2*const1*EKZero +
		      PZero*PZero*TMath::Sin(theta)*TMath::Sin(theta)
		      ) ;

    // JT Test Remove 100 MeV/c in pt to simulate eta!=0 decays
    Int_t pionindex = (int)((ptPi-0.1)*100.0 - 65.0*TMath::Abs(fBFieldG)) ; 
    Int_t kaonindex = (int)((ptK -0.1)*100.0 - 65.0*TMath::Abs(fBFieldG)) ; 
      
    if ( pionindex >= 20 ) pionindex = 399 ;
    if ( pionindex >= 0 )   effp = corrEfficiency[0][pionindex] ;
    if ( pionindex <  0 )   effp = (corrEfficiency[0][1]-corrEfficiency[0][0])*pionindex + corrEfficiency[0][0] ; // Extrapolate if reqd
    if ( effp < 0 )         effp = 0 ;

    if ( kaonindex >= 20 ) kaonindex = 399 ;
    if ( kaonindex >= 0 )   effk = corrEfficiency[1][kaonindex] ;
    if ( kaonindex <  0 )   effk = (corrEfficiency[1][1]-corrEfficiency[1][0])*kaonindex + corrEfficiency[1][0] ; // Extrapolate if reqd
    if ( effk < 0 )         effk = 0 ;

    // Note that we assume that the Kaon Decay efficiency has already been inlcuded in the kaon efficiency used here.
      
    sum += effp * effk ;
 
  }    
  
  Double_t mean =sum/360; 
  return mean ;
  
}

KMCProbe* KMCDetector::PrepareKalmanTrack(double pt, double lambda, double mass, int charge, double phi, double x,double y, double z)
{
  // Prepare trackable Kalman track at the farthest position
  //
  // Set track parameters
  // Assume track started at (0,0,0) and shoots out on the X axis, and B field is on the Z axis
  fProbe.Reset();
  fProbe.SetMass(mass);
  KMCProbe* probe = new KMCProbe(fProbe);
  double *trPars = (double*)probe->GetParameter();
  double *trCov  = (double*)probe->GetCovariance();
  double xyz[3] = {x,y,z};
  probe->Global2LocalPosition(xyz,phi);
  probe->Set(xyz[0],phi,trPars,trCov);
  trPars[KMCProbe::kY] = xyz[1];
  trPars[KMCProbe::kZ] = xyz[2];
  trPars[KMCProbe::kSnp] = 0;                       //            track along X axis at the vertex
  trPars[KMCProbe::kTgl] = TMath::Tan(lambda);                // dip
  trPars[KMCProbe::kPtI] = charge/pt;               //            q/pt      
  //
  // put tiny errors to propagate to the outer-most radius
  trCov[KMCProbe::kY2] = trCov[KMCProbe::kZ2] = trCov[KMCProbe::kSnp2] = trCov[KMCProbe::kTgl2] = trCov[KMCProbe::kPtI2] = 1e-20;
  fProbe = *probe;  // store original track
  //
  // propagate to last layer
  fLastActiveLayerTracked = 0;
  for (Int_t j=0; j<=fLastActiveLayer; j++) {
    KMCLayer* lr = GetLayer(j);
    lr->Reset();
    //
    if (!PropagateToLayer(probe,lr,1)) break;
    if (!probe->CorrectForMeanMaterial(lr, kFALSE)) break;
    //
    lr->fClCorr.Set(probe->GetY(),probe->GetZ(), probe->GetX(), probe->GetAlpha());
    if (!lr->IsDead()) fLastActiveLayerTracked = j;
  }
  probe->ResetCovMat();// reset cov.matrix
  printf("Last active layer trracked: %d (out of %d)\n",fLastActiveLayerTracked,fLastActiveLayer);
  //
  return probe;
}



TGraph * KMCDetector::GetGraphMomentumResolution(Int_t color, Int_t linewidth) {
  //
  // returns the momentum resolution 
  //
  
  TGraph *graph = new TGraph(20, fTransMomenta, fMomentumRes);
  graph->SetTitle("Momentum Resolution .vs. Pt" ) ;
  //  graph->GetXaxis()->SetRangeUser(0.,5.0) ;
  graph->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->SetTitle("Momentum Resolution (%)") ;
  graph->GetYaxis()->CenterTitle();

  graph->SetMaximum(20) ;
  graph->SetMinimum(0.1) ;
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetLineWidth(linewidth);

  return graph;

}

TGraph * KMCDetector::GetGraphPointingResolution(Int_t axis, Int_t color, Int_t linewidth) {
 
  // Returns the pointing resolution
  // axis = 0 ... rphi pointing resolution
  // axis = 1 ... z pointing resolution
  //

  TGraph * graph =  0;

  if (axis==0) {
    graph = new TGraph ( 20, fTransMomenta, fResolutionRPhi ) ;
    graph->SetTitle("R-#phi Pointing Resolution .vs. Pt" ) ;
    graph->GetYaxis()->SetTitle("R-#phi Pointing Resolution (#mum)") ;
  } else {
    graph =  new TGraph ( 20, fTransMomenta, fResolutionZ ) ;
    graph->SetTitle("Z Pointing Resolution .vs. Pt" ) ;
    graph->GetYaxis()->SetTitle("Z Pointing Resolution (#mum)") ;
  }
  
  graph->SetMinimum(1) ;
  graph->SetMaximum(300.1) ;
  graph->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->CenterTitle();
  
  graph->SetLineWidth(linewidth);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  
  return graph;

}


TGraph * KMCDetector::GetGraphPointingResolutionTeleEqu(Int_t axis,Int_t color, Int_t linewidth) {
  //
  // returns the Pointing resolution (accoring to Telescope equation)
  // axis =0 ... in rphi
  // axis =1 ... in z
  //
  
  Double_t resolution[20];

  Double_t layerResolution[2];
  Double_t layerRadius[2];
  Double_t layerThickness[2];

  Int_t count =0; // search two first active layers
  printf("Telescope equation for layers:  ");
  for (Int_t i = 0; i<fLayers.GetEntries(); i++) {
    KMCLayer *l = (KMCLayer*)fLayers.At(i);
    if (!l->IsDead() && l->fR>0) {
      layerRadius[count]     = l->fR;
      layerThickness[count]  = l->fx2X0;
      if (axis==0) {
	layerResolution[count] = l->fPhiRes;
      } else {
	layerResolution[count] = l->fZRes;
      }
      printf("%s, ",l->GetName());
      count++;
    }
    if (count>=2) break;	
  }
  printf("\n");

  Double_t pt, momentum, thickness,aMCS ;
  Double_t lambda = TMath::Pi()/2.0 - 2.0*TMath::ATan(TMath::Exp(-1*fAvgRapidity)); 

  for ( Int_t i = 0 ; i < 20 ; i++ ) { 
    // Reference data as if first two layers were acting all alone 
    pt  =  fTransMomenta[i]  ;
    momentum = pt / TMath::Cos(lambda)   ;  // Total momentum
    resolution[i] =  layerResolution[0]*layerResolution[0]*layerRadius[1]*layerRadius[1] 
      +  layerResolution[1]*layerResolution[1]*layerRadius[0]*layerRadius[0] ;
    resolution[i] /= ( layerRadius[1] - layerRadius[0] ) * ( layerRadius[1] - layerRadius[0] ) ;
    thickness = layerThickness[0] / TMath::Sin(TMath::Pi()/2 - lambda) ;
    aMCS = ThetaMCS(fParticleMass, thickness, momentum) ;
    resolution[i] += layerRadius[0]*layerRadius[0]*aMCS*aMCS ;
    resolution[i] =  TMath::Sqrt(resolution[i]) * 10000.0 ;  // result in microns
  }



  TGraph* graph = new TGraph ( 20, fTransMomenta, resolution ) ;
   
  if (axis==0) {
    graph->SetTitle("RPhi Pointing Resolution .vs. Pt" ) ;
    graph->GetYaxis()->SetTitle("RPhi Pointing Resolution (#mum) ") ;
  } else {
    graph->SetTitle("Z Pointing Resolution .vs. Pt" ) ;
    graph->GetYaxis()->SetTitle("Z Pointing Resolution (#mum) ") ;
  }
  graph->SetMinimum(1) ;
  graph->SetMaximum(300.1) ;
  graph->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->CenterTitle();
  
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetLineStyle(kDashed);
  graph->SetLineWidth(linewidth);

  return graph;

}

TGraph * KMCDetector::GetGraphRecoEfficiency(Int_t particle,Int_t color, Int_t linewidth) {
  //
  // particle = 0 ... choosen particle (setted particleMass)
  // particle = 1 ... Pion
  // particle = 2 ... Kaon
  // particle = 3 ... D0
  //
  Double_t lambda = TMath::Pi()/2.0 - 2.0*TMath::ATan(TMath::Exp(-1*fAvgRapidity)); 
  
  Double_t particleEfficiency[20]; // with chosen particle mass
  Double_t kaonEfficiency[20], pionEfficiency[20], d0efficiency[20]; 
  Double_t partEfficiency[2][20];
  
  if (particle != 0) {
    // resulting Pion and Kaon efficiency scaled with overall efficiency
    Double_t doNotDecayFactor;
    for ( Int_t massloop = 0 ; massloop < 2 ; massloop++) { //0-pion, 1-kaon
      
      for ( Int_t j = 0 ; j < 20 ; j++ ) { 
	// JT Test Let the kaon decay.  If it decays inside the TPC ... then it is gone; for all decays < 130 cm.
	Double_t momentum = fTransMomenta[j] / TMath::Cos(lambda)           ;  // Total momentum at average rapidity
	if ( massloop == 1 ) { // KAON
	  doNotDecayFactor  = TMath::Exp(-130/(371*momentum/KaonMass)) ;  // Decay length for kaon is 371 cm.
	  kaonEfficiency[j] = fEfficiency[1][j] * AcceptanceOfTpcAndSi*doNotDecayFactor ;
	} else { // PION
	  doNotDecayFactor = 1.0 ;
	  pionEfficiency[j] = fEfficiency[0][j] * AcceptanceOfTpcAndSi*doNotDecayFactor ;	
	}
	partEfficiency[0][j] = pionEfficiency[j];
	partEfficiency[1][j] = kaonEfficiency[j];
      }      
    }
    
    // resulting estimate of the D0 efficiency
    for ( Int_t j = 0 ; j < 20 ; j++ ) {
      d0efficiency[j] = D0IntegratedEfficiency(fTransMomenta[j],partEfficiency);
    }
  } else { 
    for ( Int_t j = 0 ; j < 20 ; j++ ) { 
      particleEfficiency[j] = fEfficiency[2][j]* AcceptanceOfTpcAndSi;
      // NOTE: Decay factor (see kaon) should be included to be realiable
    }
  }

  for ( Int_t j = 0 ; j < 20 ; j++ ) { 
    pionEfficiency[j]     *= 100;
    kaonEfficiency[j]     *= 100;
    d0efficiency[j]       *= 100;
    particleEfficiency[j] *= 100;
  }
 
  TGraph * graph =  0;
  if (particle==0) {
    graph = new TGraph ( 20, fTransMomenta, particleEfficiency ) ; // choosen mass
    graph->SetLineWidth(1);
  }  else if (particle==1) {
    graph = new TGraph ( 20, fTransMomenta, pionEfficiency ) ;
    graph->SetLineWidth(1);
  }  else if (particle ==2) {
    graph = new TGraph ( 20, fTransMomenta, kaonEfficiency ) ;
    graph->SetLineWidth(1);
  }  else if (particle ==3) {
    graph = new TGraph ( 20, fTransMomenta, d0efficiency ) ;
    graph->SetLineStyle(kDashed);
  } else 
    return 0;

  graph->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->SetTitle("Efficiency (%)") ;
  graph->GetYaxis()->CenterTitle();
	  
  graph->SetMinimum(0.01) ; 
  graph->SetMaximum(100)  ; 

  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetLineWidth(linewidth);

  return graph;
}

TGraph * KMCDetector::GetGraphRecoFakes(Int_t particle,Int_t color, Int_t linewidth) {
  //
  // particle = 0 ... choosen particle (setted particleMass)
  // particle = 1 ... Pion
  // particle = 2 ... Kaon
  //

  Double_t lambda = TMath::Pi()/2.0 - 2.0*TMath::ATan(TMath::Exp(-1*fAvgRapidity)); 
  
  Double_t particleFake[20]; // with chosen particle mass
  Double_t kaonFake[20], pionFake[20];
  Double_t partFake[2][20];
  
  if (particle != 0) {
    // resulting Pion and Kaon efficiency scaled with overall efficiency
    Double_t doNotDecayFactor;
    for ( Int_t massloop = 0 ; massloop < 2 ; massloop++) { //0-pion, 1-kaon
      
      for ( Int_t j = 0 ; j < 20 ; j++ ) { 
	// JT Test Let the kaon decay.  If it decays inside the TPC ... then it is gone; for all decays < 130 cm.
	Double_t momentum = fTransMomenta[j] / TMath::Cos(lambda)           ;  // Total momentum at average rapidity
	if ( massloop == 1 ) { // KAON
	  doNotDecayFactor  = TMath::Exp(-130/(371*momentum/KaonMass)) ;  // Decay length for kaon is 371 cm.
	  kaonFake[j] = fFake[1][j] /( doNotDecayFactor) ;
	} else { // PION
	  pionFake[j] = fFake[0][j] ;	
	}
	partFake[0][j] = pionFake[j];
	partFake[1][j] = kaonFake[j];
      }      
    }
    
  } else { 
    for ( Int_t j = 0 ; j < 20 ; j++ ) { 
      particleFake[j] = fFake[2][j];
      // NOTE: Decay factor (see kaon) should be included to be realiable
    }
  }

  for ( Int_t j = 0 ; j < 20 ; j++ ) { 
    pionFake[j]     *= 100;
    kaonFake[j]     *= 100;
    particleFake[j] *= 100;
  }
 
  TGraph * graph =  0;
  if (particle==0) {
    graph = new TGraph ( 20, fTransMomenta, particleFake ) ; // choosen mass
    graph->SetLineWidth(1);
  }  else if (particle==1) {
    graph = new TGraph ( 20, fTransMomenta, pionFake ) ;
    graph->SetLineWidth(1);
  }  else if (particle ==2) {
    graph = new TGraph ( 20, fTransMomenta, kaonFake ) ;
    graph->SetLineWidth(1);
  } 
  
  graph->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->SetTitle("Fake (%)") ;
  graph->GetYaxis()->CenterTitle();
	  
  graph->SetMinimum(0.01) ; 
  graph->SetMaximum(100)  ; 

  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetLineWidth(linewidth);

  return graph;
}


TGraph* KMCDetector::GetGraphImpactParam(Int_t mode, Int_t axis, Int_t color, Int_t linewidth) {
  //
  // returns the Impact Parameter d0 (convolution of pointing resolution and vtx resolution)
  // mode 0: impact parameter (convolution of pointing and vertex resolution)
  // mode 1: pointing resolution
  // mode 2: vtx resolution 
  
  
  TGraph *graph = new TGraph();

  //  TFormula vtxResRPhi("vtxRes","50-2*x"); // 50 microns at pt=0, 15 microns at pt =20 ?
  TFormula vtxResRPhi("vtxRes","35/(x+1)+10"); // 
  TFormula vtxResZ("vtxResZ","600/(x+6)+10"); // 
    
  TGraph *trackRes = GetGraphPointingResolution(axis,1);
  Double_t *pt = trackRes->GetX();
  Double_t *trRes = trackRes->GetY();
  for (Int_t ip =0; ip<trackRes->GetN(); ip++) {
    Double_t vtxRes = 0;
    if (axis==0) 
      vtxRes = vtxResRPhi.Eval(pt[ip]);
    else 
      vtxRes = vtxResZ.Eval(pt[ip]);
    
    if (mode==0)
      graph->SetPoint(ip,pt[ip],TMath::Sqrt(vtxRes*vtxRes+trRes[ip]*trRes[ip]));
    else if (mode ==1)
      graph->SetPoint(ip,pt[ip],trRes[ip]);
    else
      graph->SetPoint(ip,pt[ip],vtxRes);
  }
  
  graph->SetTitle("d_{0} r#phi resolution .vs. Pt" ) ;
  graph->GetYaxis()->SetTitle("d_{0} r#phi resolution (#mum)") ;
  
  graph->SetMinimum(1) ;
  graph->SetMaximum(300.1) ;
  graph->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->CenterTitle();
  
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetLineWidth(linewidth);

  return graph;

}

TGraph* KMCDetector::GetGraph(Int_t number, Int_t color, Int_t linewidth) {
  // 
  // returns graph according to the number
  //
  switch(number) {
  case 1:
    return GetGraphPointingResolution(0,color, linewidth); // dr
  case 2:
    return GetGraphPointingResolution(1,color, linewidth); // dz
  case 3:
    return GetGraphPointingResolutionTeleEqu(0,color, linewidth); // dr - tele
  case 4:
    return GetGraphPointingResolutionTeleEqu(1,color, linewidth); // dz - tele
  case 5:
    return GetGraphMomentumResolution(color, linewidth); // pt resolution
  case 10:
    return GetGraphRecoEfficiency(0, color, linewidth);  // tracked particle
  case 11:
    return GetGraphRecoEfficiency(1, color, linewidth);  // eff. pion
  case 12:
    return GetGraphRecoEfficiency(2, color, linewidth);  // eff. kaon
  case 13: 
    return GetGraphRecoEfficiency(3, color, linewidth);  // eff. D0
  case 15:
    return GetGraphRecoFakes(0, color, linewidth);  // Fake tracked particle
  case 16:
    return GetGraphRecoFakes(1, color, linewidth);  // Fake pion
  case 17:
    return GetGraphRecoFakes(2, color, linewidth);  // Fake kaon
  default:
    printf(" Error: chosen graph number not valid\n");
  }
  return 0;

}

void KMCDetector::MakeAliceAllNew(Bool_t flagTPC,Bool_t flagMon, int setVer) {
  
  // All New configuration with X0 = 0.3 and resolution = 4 microns
  
  AddLayer((char*)"bpipe_its",2.0,0.0022, 0.092); // beam pipe, 0.5 mm Be
  AddLayer((char*)"vertex_its",     0,     0); // dummy vertex for matrix calculation
  if (fgVtxConstraint[0]>0 && fgVtxConstraint[1]>0) {
    printf("vertex will work as constraint: %.4f %.4f\n",fgVtxConstraint[0],fgVtxConstraint[1]);
    SetResolution((char*)"vertex_its",fgVtxConstraint[0],fgVtxConstraint[1]);
  }
  //
  // new ideal Pixel properties?
  Double_t x0     = 0.0050;
  Double_t resRPhi = 0.0006;
  Double_t resZ   = 0.0006;
  Double_t xrho = 0.0116;  // assume 0.5mm of silicon for eloss
  //
  if (flagMon) {
    x0 *= 3./5.;
    xrho *=3./5.;
    resRPhi = 0.0004;
    resZ   = 0.0004;
  }

  double sclD = 1;
  double sclZ = 1;

  // default all new  
  if (setVer<=0) {
    AddLayer((char*)"its1",  2.2 ,  x0, xrho, resRPhi, resZ); 
    AddLayer((char*)"its2",  3.8 ,  x0, xrho, resRPhi, resZ); 
    AddLayer((char*)"its3",  6.8 ,  x0, xrho, resRPhi, resZ); 
    AddLayer((char*)"its4", 12.4 ,  x0, xrho, resRPhi, resZ); 
    AddLayer((char*)"its5", 23.5 ,  x0, xrho, resRPhi, resZ); 
    AddLayer((char*)"its6", 39.6 ,  x0, xrho, resRPhi, resZ); 
    AddLayer((char*)"its7", 43.0 ,  x0, xrho, resRPhi, resZ); 
  }
  else if (setVer==1) {
    AddLayer((char*)"its1",  2.2 ,  x0, xrho, resRPhi, resZ); 
    AddLayer((char*)"its2",  2.8 ,  x0, xrho, resRPhi, resZ); 
    AddLayer((char*)"its3",  3.6 ,  x0, xrho, resRPhi, resZ); 
    AddLayer((char*)"its4", 20.0 ,  x0, xrho, resRPhi*sclD, resZ*sclZ); 
    AddLayer((char*)"its5", 22.0 ,  x0, xrho, resRPhi*sclD, resZ*sclZ); 
    AddLayer((char*)"its6", 43.0 ,  x0, xrho, resRPhi*sclD, resZ*sclZ); 
    AddLayer((char*)"its7", 43.6 ,  x0, xrho, resRPhi*sclD, resZ*sclZ); 
    //
    /*
    UInt_t patt[] = {
      KMCTrackSummary::Bits(1,1,1),
      KMCTrackSummary::Bits(0,0,0,1,1),
      KMCTrackSummary::Bits(0,0,0,0,0,1,1)
    };
    RequirePattern( patt, sizeof(patt)/sizeof(UInt_t));
    */
  }
  else if (setVer==2) {
    AddLayer((char*)"its1",  2.2 ,  x0, xrho, resRPhi, resZ); 
    AddLayer((char*)"its1a", 2.8 ,  x0, xrho, resRPhi, resZ); 
    AddLayer((char*)"its2",  3.6 ,  x0, xrho, resRPhi, resZ); 
    AddLayer((char*)"its2a", 4.2 ,  x0, xrho, resRPhi, resZ); 
    AddLayer((char*)"its3", 20.0 ,  x0, xrho, resRPhi*sclD, resZ*sclZ); 
    AddLayer((char*)"its4", 22.0 ,  x0, xrho, resRPhi*sclD, resZ*sclZ); 
    AddLayer((char*)"its5", 33.0 ,  x0, xrho, resRPhi*sclD, resZ*sclZ); 
    AddLayer((char*)"its6", 43.0 ,  x0, xrho, resRPhi*sclD, resZ*sclZ); 
    AddLayer((char*)"its7", 43.6 ,  x0, xrho, resRPhi*sclD, resZ*sclZ); 
    //
    /*
    UInt_t patt[] = {
      KMCTrackSummary::Bits(1,1,1,1),
      KMCTrackSummary::Bits(0,0,0,0,1,1),
      KMCTrackSummary::Bits(0,0,0,0,0,0,1,1,1)
    };
    RequirePattern( patt, sizeof(patt)/sizeof(UInt_t));
    */
  }
   /*
  // last 2 layers strips
  double resRPStr=0.0020, resZStr = 0.0830;
  AddLayer((char*)"its1",  2.2 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its2",  3.8 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its3",  6.8 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its4", 12.4 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its5", 23.5 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its6", 39.6 ,  x0, xrho, resRPStr, resZStr); 
  AddLayer((char*)"its7", 43.0 ,  x0, xrho, resRPStr, resZStr); 
   */
  //*/
  /*
  // resolution scaled with R as res(2.2) * R/2.2
  AddLayer((char*)"its1",  2.2 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its2",  3.8 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its3",  6.8 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its4", 12.4 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its5", 23.5 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its6", 39.6 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its7", 43.0 ,  x0, xrho, resRPhi, resZ); 
  KMCLayer* lr0 = 0;
  for (int i=0;i<GetNActiveITSLayers();i++) {
    KMCLayer *lr = GetActiveLayer(i);
    if (lr->IsVertex()) continue;
    if (lr0==0) {
      lr0=lr; 
      //  continue;
    }
    double scl = 5*lr->GetRadius()/lr0->GetRadius();
    SetResolution((char*)lr->GetName(), resRPhi*scl, resZ*scl*4);
  }
  */
  /*
  // 1st 2 layers double
  AddLayer((char*)"its1",  2.2 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its1a",  2.8 ,  x0, xrho, resRPhi, resZ); 

  AddLayer((char*)"its2",  3.8 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its2a",  4.2 ,  x0, xrho, resRPhi, resZ); 

  AddLayer((char*)"its3a",  6.4 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its3",  6.8 ,  x0, xrho, resRPhi, resZ); 

  AddLayer((char*)"its4", 12.4 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its5", 23.5 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its6", 39.6 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its7", 43.0 ,  x0, xrho, resRPhi, resZ); 
  */
  /*
  // last 4 layers doubled
  AddLayer((char*)"its1",  2.2 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its2",  3.8 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its3",  6.8 ,  x0, xrho, resRPhi, resZ); 

  AddLayer((char*)"its4", 12.4 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its4a", 12.8 ,  x0, xrho, resRPhi, resZ); 

  AddLayer((char*)"its5", 23.5 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its5a", 23.9 ,  x0, xrho, resRPhi, resZ); 

  AddLayer((char*)"its6", 39.6 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its6a", 40.0 ,  x0, xrho, resRPhi, resZ); 

  AddLayer((char*)"its7", 43.0 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its7a", 43.4 ,  x0, xrho, resRPhi, resZ); 
  */

  /* //last 3 lr together 
  AddLayer((char*)"its1",  2.2 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its2",  3.8 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its3",  6.8 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its4", 12.4 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its5", 42.2 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its6", 42.6 ,  x0, xrho, resRPhi, resZ); 
  AddLayer((char*)"its7", 43.0 ,  x0, xrho, resRPhi, resZ); 
  */
  if (flagTPC) {
    AddTPC(0.1,0.1);                        // TPC
  }
  PrintLayout();
}

void KMCDetector::MakeAliceCurrent(Bool_t flagTPC, Int_t AlignResiduals) {

  // Numbers taken from 
  // 2010 JINST 5 P03003 - Alignment of the ALICE Inner Tracking System with cosmic-ray tracks
  // number for misalingment: private communication with Andrea Dainese

  //  1.48e-01, 2.48e-01,2.57e-01, 1.34e-01, 3.34e-01,3.50e-01, 2.22e-01, 2.38e-01,2.25e-01};

  AddLayer((char*)"vertex_its",     0,     0); // dummy vertex for matrix calculation
  if (fgVtxConstraint[0]>0 && fgVtxConstraint[1]>0) {
    printf("vertex will work as constraint: %.4f %.4f",fgVtxConstraint[0],fgVtxConstraint[1]);
    SetResolution((char*)"vertex_its",fgVtxConstraint[0],fgVtxConstraint[1]);
  }
  //
  AddLayer((char*)"bpipe_its",2.94,0.0022, 1.48e-01); // beam pipe
  AddLayer((char*)"tshld1_its",11.5,0.0065, 1.34e-01); // Thermal shield  // 1.3% /2
  AddLayer((char*)"tshld2_its",31.0,0.0108, 2.22e-01); // Thermal shield  // 1.3% /2
  //
  if (flagTPC) {
    AddTPC(0.1,0.1);                        // TPC
  }
  // Adding the ITS - current configuration
  
  if (AlignResiduals==0) {
    AddLayer((char*)"spd1_its", 3.9, 0.0114, 2.48e-01, 0.0012, 0.0130);
    AddLayer((char*)"spd2_its", 7.6, 0.0114, 2.57e-01, 0.0012, 0.0130);
    AddLayer((char*)"sdd1_its",15.0, 0.0113, 3.34e-01, 0.0035, 0.0025);
    AddLayer((char*)"sdd2_its",23.9, 0.0126, 3.50e-01, 0.0035, 0.0025);
    AddLayer((char*)"ssd1_its",38.0, 0.0083, 2.38e-01, 0.0020, 0.0830);
    AddLayer((char*)"ssd2_its",43.0, 0.0086, 2.25e-01, 0.0020, 0.0830);
    /*
    UInt_t patt[] = {
      KMCTrackSummary::Bits(1,1),
      KMCTrackSummary::Bits(0,0,1,1),
      KMCTrackSummary::Bits(0,0,0,0,1,1)
    };
    RequirePattern( patt, sizeof(patt)/sizeof(UInt_t));
    */
  } else if (AlignResiduals==1) {

    // tracking errors ...
    // (Additional systematic errors due to misalignments) ... 
    // itsRecoParam->SetClusterMisalErrorYBOn(0.0010,0.0030,0.0500,0.0500,0.0020,0.0020);  // [cm]
    // itsRecoParam->SetClusterMisalErrorZBOn(0.0050,0.0050,0.0050,0.0050,0.1000,0.1000);

    AddLayer((char*)"spd1_its", 3.9, 0.0114, 2.48e-01, TMath::Sqrt(0.0012*0.0012+0.0010*0.0010), 
	     TMath::Sqrt(0.0130*0.0130+0.0050*0.0050));
    AddLayer((char*)"spd2_its", 7.6, 0.0114, 2.57e-01, TMath::Sqrt(0.0012*0.0012+0.0030*0.0030),
	     TMath::Sqrt(0.0130*0.0130+0.0050*0.0050));
    AddLayer((char*)"sdd1_its",15.0, 0.0113, 3.34e-01, TMath::Sqrt(0.0035*0.0035+0.0100*0.0100),
	     TMath::Sqrt(0.0025*0.0025+0.0050*0.0050));
    AddLayer((char*)"sdd2_its",23.9, 0.0126, 3.50e-01, TMath::Sqrt(0.0035*0.0035+0.0100*0.0100),
	     TMath::Sqrt(0.0025*0.0025+0.0050*0.0050));
    AddLayer((char*)"ssd1_its",38.0, 0.0083, 2.38e-01, TMath::Sqrt(0.0020*0.0020+0.0020*0.0020), 
	     TMath::Sqrt(0.0830*0.0830+0.1000*0.1000));
    AddLayer((char*)"ssd2_its",43.0, 0.0086, 2.25e-01, TMath::Sqrt(0.0020*0.0020+0.0020*0.0020),
	     TMath::Sqrt(0.0830*0.0830+0.1000*0.1000));   
    
  } else if (AlignResiduals==2) {

    
    // tracking errors ... PLUS ... module misalignment
    
    // itsRecoParam->SetClusterMisalErrorYBOn(0.0010,0.0030,0.0500,0.0500,0.0020,0.0020);  // [cm]
    // itsRecoParam->SetClusterMisalErrorZBOn(0.0050,0.0050,0.0050,0.0050,0.1000,0.1000);
    
    //  the ITS modules are misalignment with small gaussian smearings with
    //  sigmarphi ~ 8, 10, 10 micron in SPD, SDD, SSD
    
    AddLayer((char*)"spd1_its", 3.9, 0.0114, 2.48e-01, TMath::Sqrt(0.0012*0.0012+0.0010*0.0010+0.0008*0.0008), 
	     TMath::Sqrt(0.0130*0.0130+0.0050*0.0050));
    AddLayer((char*)"spd2_ots", 7.6, 0.0114, 2.57e-01, TMath::Sqrt(0.0012*0.0012+0.0030*0.0030+0.0008*0.0008),
	     TMath::Sqrt(0.0130*0.0130+0.0050*0.0050));
    AddLayer((char*)"sdd1_ots",15.0, 0.0113, 3.34e-01, TMath::Sqrt(0.0035*0.0035+0.0500*0.0500+0.0010*0.0010),
	     TMath::Sqrt(0.0025*0.0025+0.0050*0.0050));
    AddLayer((char*)"sdd2_its",23.9, 0.0126, 3.50e-01, TMath::Sqrt(0.0035*0.0035+0.0500*0.0500+0.0010*0.0010),
	     TMath::Sqrt(0.0025*0.0025+0.0050*0.0050));
    AddLayer((char*)"ssd1_its",38.0, 0.0083, 2.38e-01, TMath::Sqrt(0.0020*0.0020+0.0020*0.0020+0.0010*0.0010), 
	     TMath::Sqrt(0.0830*0.0830+0.1000*0.1000));
    AddLayer((char*)"ssd2_its",43.0, 0.0086, 2.25e-01, TMath::Sqrt(0.0020*0.0020+0.0020*0.0020+0.0010*0.0010),
	     TMath::Sqrt(0.0830*0.0830+0.1000*0.1000)); 

  } else {
      
      //  the ITS modules are misalignment with small gaussian smearings with
      //  sigmarphi ~ 8, 10, 10 micron in SPD, SDD, SSD
      //  unknown in Z ????

    AddLayer((char*)"spd1_its", 3.9, 0.0114, 2.48e-01, TMath::Sqrt(0.0012*0.0012+0.0008*0.0008), 
	     TMath::Sqrt(0.0130*0.0130+0.000*0.000));
    AddLayer((char*)"spd2_its", 7.6, 0.0114, 2.57e-01, TMath::Sqrt(0.0012*0.0012+0.0008*0.0008),
	     TMath::Sqrt(0.0130*0.0130+0.000*0.000));
    AddLayer((char*)"sdd1_its",15.0, 0.0113, 3.34e-01, TMath::Sqrt(0.0035*0.0035+0.0010*0.0010),
	     TMath::Sqrt(0.0025*0.0025+0.000*0.000));
    AddLayer((char*)"sdd2_its",23.9, 0.0126, 3.50e-01, TMath::Sqrt(0.0035*0.0035+0.0010*0.0010),
	     TMath::Sqrt(0.0025*0.0025+0.000*0.000));
    AddLayer((char*)"ssd1_its",38.0, 0.0083, 2.38e-01, TMath::Sqrt(0.0020*0.0020+0.0010*0.0010), 
	     TMath::Sqrt(0.0830*0.0830+0.000*0.000));
    AddLayer((char*)"ssd2_its",43.0, 0.0086, 2.25e-01, TMath::Sqrt(0.0020*0.0020+0.0010*0.0010),
	     TMath::Sqrt(0.0830*0.0830+0.000*0.000));       
  }
  
}


void KMCDetector::MakeStandardPlots(Bool_t add, Int_t color, Int_t linewidth,Bool_t onlyPionEff) {
  //
  // Produces the standard performace plots
  //
 
  if (!add) {

    TCanvas *c1 = new TCanvas("c1","c1");//,100,100,500,500);  
    c1->Divide(2,2);
    
    c1->cd(1);  gPad->SetGridx();   gPad->SetGridy(); 
    gPad->SetLogx(); 
    TGraph *eff = GetGraphRecoEfficiency(1,color,linewidth);
    eff->SetTitle("Efficiencies");
    eff->Draw("AL");
    if (!onlyPionEff) {
      GetGraphRecoEfficiency(2,color,linewidth)->Draw("L");
      GetGraphRecoEfficiency(3,color,linewidth)->Draw("L");
    }
    c1->cd(2); gPad->SetGridx();   gPad->SetGridy(); 
    gPad->SetLogy();  gPad->SetLogx(); 
    GetGraphMomentumResolution(color,linewidth)->Draw("AL");
    
    c1->cd(3); gPad->SetGridx();   gPad->SetGridy(); 
    gPad->SetLogx(); 
    GetGraphPointingResolution(0,color,linewidth)->Draw("AL");
    
    c1->cd(4); gPad->SetGridx();   gPad->SetGridy(); 
    gPad->SetLogx(); 
    GetGraphPointingResolution(1,color,linewidth)->Draw("AL");

  } else {

    TVirtualPad *c1 = gPad->GetMother();

    c1->cd(1);
    GetGraphRecoEfficiency(1,color,linewidth)->Draw("L");
    if (!onlyPionEff) {
      GetGraphRecoEfficiency(2,color,linewidth)->Draw("L");
      GetGraphRecoEfficiency(3,color,linewidth)->Draw("L");
    }
    c1->cd(2); GetGraphMomentumResolution(color,linewidth)->Draw("L");
    
    c1->cd(3); GetGraphPointingResolution(0,color,linewidth)->Draw("L");
    
    c1->cd(4); GetGraphPointingResolution(1,color,linewidth)->Draw("L");
    
  }

}

void KMCDetector::ApplyMS(KMCProbe* trc, double x2X0) const
{
  // simulate random modification of track params due to the MS
  if (x2X0<=0) return;
  double alpha = trc->GetAlpha(); // store original alpha
  double mass = trc->GetMass();
  //
  double snp = trc->GetSnp();
  double dip = trc->GetTgl();
  Double_t angle=TMath::Sqrt((1.+ dip*dip)/((1-snp)*(1.+snp)));
  x2X0 *= angle;
  //
  static double covCorr[15],covDum[21]={0};
  static double mom[3],pos[3];
  double *cov = (double*) trc->GetCovariance();
  memcpy(covCorr,cov,15*sizeof(double));
  trc->GetXYZ(pos);
  trc->GetPxPyPz(mom);
  double pt2 = mom[0]*mom[0]+mom[1]*mom[1];
  double pt = TMath::Sqrt(pt2);
  double ptot2 = pt2 + mom[2]*mom[2];
  double ptot  = TMath::Sqrt(ptot2);
  double beta = ptot/TMath::Sqrt(ptot2 + mass*mass);
  double sigth = TMath::Sqrt(x2X0)*0.014/(ptot*beta);
  //
  // a la geant
  double phiSC = gRandom->Rndm()*TMath::Pi();
  double thtSC = gRandom->Gaus(0,1.4142*sigth);
  //  printf("MS phi: %+.5f tht: %+.5f\n",phiSC,thtSC);
  double sn = TMath::Sin(thtSC);
  double dx = sn*TMath::Sin(phiSC);
  double dy = sn*TMath::Cos(phiSC);  
  double dz = TMath::Cos(thtSC);
  double v[3];
  //  printf("Before: %+.3e %+.3e %+.3e | MS: %+.3e %+.3e\n",mom[0],mom[1],mom[2],thtSC,phiSC);
  for (int i=3;i--;) mom[i] /= ptot;
  double vmm = TMath::Sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
  if (!IsZero(pt)) {
    double pd1 = mom[0]/vmm;
    double pd2 = mom[1]/vmm;
    v[0] = pd1*mom[2]*dx - pd2*dy + mom[0]*dz;
    v[1] = pd2*mom[2]*dx + pd1*dy + mom[1]*dz;
    v[2] = -vmm*dx                + mom[2]*dz;
  }
  else {
    v[0] = dx;
    v[1] = dy;
    v[2] = dz*TMath::Sign(1.,mom[2]);
  }
  double nrm = TMath::Sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  //  printf("before :%+e %+e %+e  || %+e %+e %+e %+e\n",mom[0],mom[1],mom[2],  sigth, x2X0, pt, beta);
  //  trc->Print();
  // direction cosines -> p
  for (int i=3;i--;) mom[i] = ptot*v[i]/nrm;
  //  printf("After : %+.3e %+.3e %+.3e\n",mom[0],mom[1],mom[2]);
  trc->Set(pos,mom,covDum,trc->Charge());
  //
  trc->Rotate(alpha);
  memcpy(cov,covCorr,15*sizeof(double));
  //
}

//________________________________________________________________________________
Bool_t KMCDetector::TransportKalmanTrackWithMS(KMCProbe *probTr, int maxLr) const
{
  // Transport track till layer maxLr, applying random MS
  //
  for (Int_t j=0; j<maxLr; j++) {
    KMCLayer* lr0 = (KMCLayer*)fLayers.At(j);
    KMCLayer* lr  = (KMCLayer*)fLayers.At(j+1);
    if (!lr0->IsVertex()) {
      ApplyMS(probTr,lr0->fx2X0); // apply MS
      if (!probTr->CorrectForMeanMaterial(lr0,kFALSE)) return kFALSE; 
    }
    //
    if (!PropagateToLayer(probTr,lr,1)) return kFALSE;
    // store randomized cluster local coordinates and phi
    lr->ResetBgClusters();
    double rz,ry;
    gRandom->Rannor(rz,ry);
    lr->GetMCCluster()->Set(probTr->GetY()+ry*lr->GetPhiRes(),probTr->GetZ()+rz*lr->GetZRes(), 
			       probTr->GetX(), probTr->GetAlpha() );    
    //
  }
  //
  return kTRUE;
}

//________________________________________________________________________________
Bool_t KMCDetector::SolveSingleTrack(Double_t mass, Double_t pt, Double_t eta, TObjArray* sumArr,int nMC, int offset)
{
  // analityc and fullMC (nMC trials) evaluaion of tracks with given kinematics.
  // the results are filled in KMCTrackSummary objects provided via summArr array
  //
  int progressP = 10;//0; // print progress in percents
  //
  progressP = int(nMC*0.01*progressP);

  if (!SolveSingleTrackViaKalman(mass,pt,eta)) return kFALSE;
  //
  // Store non-updated track errors of inward propagated seed >>>>>>>>
  int maxLr = fLastActiveITSLayer + offset;
  if (maxLr >= fLastActiveLayerTracked-1) maxLr = fLastActiveLayerTracked;
  KMCProbe probeTmp = fProbe; // original probe at vertex
  KMCLayer* lr = 0;
  for (Int_t j=1; j<=maxLr; j++) {
    lr = GetLayer(j);
    //    printf("Here0: %d\n",j);
    if (!PropagateToLayer(&probeTmp,lr,1)) return 0;
    if (j!=maxLr) if (!probeTmp.CorrectForMeanMaterial(lr, kFALSE)) return 0;
    //    printf("Prelim. Err at lr:%8s | %7.3f %7.3f\n",lr->GetName(),TMath::Sqrt(probeTmp.GetSigmaY2()),TMath::Sqrt(probeTmp.GetSigmaZ2()));
  }
  for (Int_t j=maxLr; j>0; j--) {
    lr = GetLayer(j);
    //    printf("Here1: %d\n",j);
    if (j!=maxLr) if (!PropagateToLayer(&probeTmp,lr,-1)) return 0;
    lr->fSig2EstD = probeTmp.GetSigmaY2();
    lr->fSig2EstZ = probeTmp.GetSigmaZ2();
    //    probeTmp.Print("l");
    printf("Natural Err at lr:%8s | %7.3f %7.3f\n",lr->GetName(),TMath::Sqrt(lr->fSig2EstD),TMath::Sqrt(lr->fSig2EstZ));
    if (!probeTmp.CorrectForMeanMaterial(lr, kTRUE)) return 0;
  }
  // Store non-updated track errors of inward propagated seed <<<<<<<<
  //
  int nsm = sumArr ? sumArr->GetEntriesFast() : 0;
  KMCLayer* vtx = GetLayer(0);
  //
  for (int i=0;i<nsm;i++) {
    KMCTrackSummary* tsm = (KMCTrackSummary*)sumArr->At(i);
    if (!tsm) continue;
    tsm->SetRefProbe( GetProbeTrack() ); // attach reference track (generated)
    tsm->SetAnProbe( vtx->GetAnProbe() ); // attach analitycal solution
  }
  //
  TStopwatch sw;
  sw.Start();
  for (int it=0;it<nMC;it++) {
    printf("ev: %d\n",it);
    SolveSingleTrackViaKalmanMC(offset);
    KMCProbe* trc = vtx->GetWinnerMCTrack();
    vtx->GetMCTracks()->Print();
    if (progressP==1 || (progressP>0 &&  (it%progressP)==0)) {
      printf("%d%% done |",it*100/nMC); 
      sw.Stop(); sw.Print(); sw.Start(kFALSE);
    }
    for (int ism=nsm;ism--;) { // account the track in each of summaries
      KMCTrackSummary* tsm = (KMCTrackSummary*)sumArr->At(ism);
      if (!tsm) continue;
      tsm->AddUpdCalls(GetUpdCalls());
      tsm->AddTrack(trc); 
    }
  }
  //
  sw.Stop();
  printf("Total time: "); sw.Print();
  return kTRUE;
}

//________________________________________________________________________________
Int_t KMCDetector::GetLayerID(int actID) const
{
  // find physical layer id from active id
  if (actID<0 || actID>fNActiveLayers) return -1;
  int start = actID<fNActiveITSLayers ? fLastActiveITSLayer : fLastActiveLayer;
  for (int i=start+1;i--;) {
    if (GetLayer(i)->GetActiveID()==actID) return i;   
  }
  return -1;
}

//________________________________________________________________________________
KMCProbe* KMCDetector::KalmanSmooth(int actLr, int actMin,int actMax) const
{
  // estimate kalman smoothed track params at given active lr 
  // from fit at layers actMin:actMax (excluding actLr)
  // SolveSingleTrackViaKalman must have been called before
  //
  if (actMin>actMax) swap(actMin,actMax);
  if (actMax>=fNActiveLayers) actMax = fNActiveLayers-1;
  int nlrfit = actMax-actMin;
  if (actLr>=actMin && actLr<=actMax) nlrfit-=1;
  if (nlrfit<2) {AliInfo("Need a least 2 active layers in the fit"); return 0;}
  static KMCProbe iwd,owd;
  //
  // find phisical layer id's
  int pLr  = GetLayerID(actLr);
  int pMin = GetLayerID(actMin);
  int pMax = GetLayerID(actMax);
  //
  //  printf(">>> %d %d %d\n",pLr, pMin,pMax);
  Bool_t useIwd=kFALSE, useOwd=kFALSE;
  if (pLr<pMax) { // need inward piece
    iwd = GetLayer(pMax)->fTrCorr;
    iwd.ResetCovMat();
    iwd.GetHitsPatt() = 0;
    for (int i=pMax;i>=pLr;i--) {
      KMCLayer* lr = GetLayer(i);
      //      printf("IWD %d\n",i);
      if (!lr->IsDead() && i!=pLr && i>=pMin) if (!UpdateTrack(&iwd,lr,&lr->fClCorr))  return 0;
      if (i!=pLr) {
	if (!iwd.CorrectForMeanMaterial(lr,kTRUE)) return 0; // correct for materials of this layer
	if (!PropagateToLayer(&iwd,GetLayer(i-1),-1)) return 0;      // propagate to next layer
      }
      //  printf("IWD%d:  ",i); iwd.Print("l");
    }
    useIwd = kTRUE;
  }
  if (pLr>pMin) { // need outward piece
    owd = GetLayer(pMin)->fTrCorr;
    owd.ResetCovMat();
    owd.GetHitsPatt() = 0;
    for (int i=pMin;i<=pLr;i++) {
      KMCLayer* lr = GetLayer(i);
      //      printf("OWD %d\n",i);
      if (!lr->IsDead() && i!=pLr && i<=pMax) if (!UpdateTrack(&owd,lr,&lr->fClCorr))  return 0;
      if (i!=pLr) {
	if (!owd.CorrectForMeanMaterial(lr,0)) return 0; // correct for materials of this layer
	if (!PropagateToLayer(&owd,GetLayer(i+1), 1)) return 0;      // propagate to next layer
      }
      //      printf("OWD%d:  ",i); owd.Print("l");
    }
    useOwd = kTRUE;
  }
  //
  // was this extrapolation outside the fit range?
  if (!useIwd) return (KMCProbe*)&owd; 
  if (!useOwd) return (KMCProbe*)&iwd;
  //
  // weight both tracks
  if (!iwd.Propagate(owd.GetAlpha(),owd.GetX(),fBFieldG)) return 0;
  double meas[2] = {owd.GetY(),owd.GetZ()};
  double measErr2[3] = {owd.GetSigmaY2(), owd.GetSigmaZY(), owd.GetSigmaZ2()};
  //  printf("Weighting\n");
  //  owd.Print("l");
  //  iwd.Print("l");
  if (!iwd.Update(meas,measErr2)) return 0;
  iwd.GetHitsPatt() |= owd.GetHitsPatt();

  //  printf("->\n");
  //  iwd.Print("l");

  return (KMCProbe*)&iwd;
  //
}

//________________________________________________________________________________
KMCProbe* KMCDetector::KalmanSmoothFull(int actLr, int actMin,int actMax) const
{
  // estimate kalman smoothed track params at given active lr 
  // from fit at layers actMin:actMax (excluding actLr)
  // SolveSingleTrackViaKalman must have been called before
  //
  static TClonesArray prediction("KMCProbe",10);
  static TClonesArray update("KMCProbe",10);
  static KMCProbe res;
  //
  if (actMin>actMax) swap(actMin,actMax);
  int nlrfit = actMax-actMin;
  if (actLr>=actMin && actLr<=actMax) nlrfit-=1;
  if (nlrfit<2) {AliInfo("Need a least 2 active layers in the fit"); return 0;}
  //
  // find phisical layer id's
  int pLr  = GetLayerID(actLr);
  int pMin = GetLayerID(actMin);
  int pMax = GetLayerID(actMax);
  //
  int dir=0,dirInt=0;
  if      (pLr<=pMin) dir=-1; // inward extrapolation
  else if (pLr>=pMax) dir= 1; // outward extrapolation
  else if (actMax-actLr >= actLr-actMin) dirInt = -1; // inward  interpolation (the test point is closer to inner layer)
  else    dirInt = 1;                                 // outward interpolation (the test point is closer to outer layer)
  //
  if (dir!=0) { // no sens to do smoothing: simple Kalman filtering extrapolation
    int start = dir<0 ? pMax : pMin;
    res = GetLayer(start)->fTrCorr;
    res.ResetCovMat();
    KMCLayer* lr = 0;
    for (int i=(dir<0?pMax:pMin); i!=pLr; i+=dir) { // track till nearest layer to pLr
      lr = GetLayer(i);
      if (!lr->IsDead() && !(i<pMin ||i>pMax)) if (!UpdateTrack(&res,lr,&lr->fClCorr))  return 0; // update only with layers in fit range
      if (!res.CorrectForMeanMaterial(lr,dir<0 ? kTRUE:kFALSE))   return 0; // correct for materials of this layer
      if (!PropagateToLayer(&res,GetLayer(i+dir),dir))            return 0; // propagate to next layer     
    }
    if (!res.CorrectForMeanMaterial(lr,dir<0 ? kTRUE:kFALSE))   return 0; // correct for materials of this nearest layer
    if (!PropagateToLayer(&res,GetLayer(pLr), dir)) return 0; // propagate to test layer
    return (KMCProbe*)&res;
  }
  //
  // too bad, need to do real filtering
  //
  int start = dirInt<0 ? pMax : pMin;
  int stop  = dirInt<0 ? pMin-1 : pMax+1;
  res = GetLayer(start)->fTrCorr;
  res.ResetCovMat();
  KMCLayer* lr = 0;
  int count = 0;
  for (int i=start; i!=stop; i+=dirInt) { // track in full range, storing updates and predictions
    new(prediction[count]) KMCProbe(res);
    lr = GetLayer(i);
    if (!lr->IsDead() && i!=pLr) if (!UpdateTrack(&res,lr,&lr->fClCorr))  return 0; // update only with layers in fit range
    new(update[count]) KMCProbe(res);
    if (!res.CorrectForMeanMaterial(lr,dir<0 ? kTRUE:kFALSE))   return 0; // correct for materials of this layer
    if (!PropagateToLayer(&res,GetLayer(i+dir),dir))            return 0; // propagate to next layer     
    count++;
  }
  return (KMCProbe*)&res;
  //
}

//________________________________________________________________________________
Bool_t KMCDetector::SolveSingleTrackViaKalman(Double_t mass, Double_t pt, Double_t eta)
{
  // analytical estimate of tracking resolutions
  //  fProbe.SetUseLogTermMS(kTRUE);
  //
  if (fMinITSHits>fNActiveITSLayers) {fMinITSHits = fNActiveITSLayers; printf("Redefined request of min N ITS hits to %d\n",fMinITSHits);}
  if (TMath::Abs(eta)<1e-3) fDensFactorEta = 1.;
  else {
    fDensFactorEta = TMath::Tan( 2.*TMath::ATan(TMath::Exp(-TMath::Abs(eta))) );
    fDensFactorEta = 1./TMath::Sqrt( 1. + 1./fDensFactorEta/fDensFactorEta);
  }
  double lambda = TMath::Pi()/2.0 - 2.0*TMath::ATan(TMath::Exp(-eta)); 
  KMCProbe* probe = PrepareKalmanTrack(pt,lambda,mass,-1);
  if (!probe) return kFALSE;
  //
  KMCLayer *lr = 0;
  //
  //
  // Start the track fitting --------------------------------------------------------
  //
  // Back-propagate the covariance matrix along the track. 
  // Kalman loop over the layers
  //
  KMCProbe* currTr = 0;
  lr = (KMCLayer*)fLayers.At(fLastActiveLayerTracked);
  lr->fTrCorr = *probe;
  delete probe; // rethink...
  //
  for (Int_t j=fLastActiveLayerTracked; j--; ) {  // Layer loop
    //
    KMCLayer *lrP = lr;
    lr = (KMCLayer*)fLayers.At(j);
    //
    lr->fTrCorr = lrP->fTrCorr;
    currTr = &lr->fTrCorr;
    currTr->ResetHit(lrP->GetActiveID());
    //
    // if there was a measurement on prev layer, update the track
    if (!lrP->IsDead()) { // include measurement
      KMCCluster cl(currTr->GetY(),currTr->GetZ(), currTr->GetX(), currTr->GetAlpha());
      if (!UpdateTrack(currTr,lrP,&cl))  return kFALSE;
    }
    if (!currTr->CorrectForMeanMaterial(lrP,kTRUE)) return kFALSE; // correct for materials of this layer
    if (!PropagateToLayer(currTr,lr,-1)) return kFALSE;      // propagate to current layer
    //
  } // end loop over layers
  //
  return kTRUE;
}

//____________________________________________________________
Bool_t KMCDetector::SolveSingleTrackViaKalmanMC(int offset)
{
  // MC estimate of tracking resolutions/effiencies. Requires that the SolveSingleTrackViaKalman
  // was called before, since it uses data filled by this method
  //
  // The MC tracking will be done starting from fLastActiveITSLayer + offset (before analytical estimate will be used)
  //
  // At this point, the fProbe contains the track params generated at vertex.
  // Clone it and propagate to target layer to generate hit positions affected by MS
  //
  fUpdCalls = 0.;
  KMCProbe *currTrP=0,*currTr=0;
  int maxLr = fLastActiveITSLayer + offset;
  if (maxLr >= fLastActiveLayerTracked-1) maxLr = fLastActiveLayerTracked;
  ResetMCTracks(maxLr);
  KMCLayer* lr = (KMCLayer*)fLayers.At(maxLr);
  currTr = lr->AddMCTrack(&fProbe); // start with original track at vertex
  //
  if (!TransportKalmanTrackWithMS(currTr, maxLr)) return kFALSE; // transport it to outermost layer where full MC is done
  //
  if (fLastActiveITSLayer<fLastActiveLayerTracked) { // prolongation from TPC
    // start from correct track propagated from above till maxLr
    double *covMS = (double*)currTr->GetCovariance();
    const double *covIdeal =lr->fTrCorr.GetCovariance();
    for (int i=15;i--;) covMS[i] = covIdeal[i];
  }
  else { // ITS SA: randomize the starting point
    //    double *pars = (double*)currTr->GetParameter();
    //    pars[0] += gRandom->Gaus(0,TMath::Sqrt(currTr->GetSigmaY2()));
    //    pars[1] += gRandom->Gaus(0,TMath::Sqrt(currTr->GetSigmaZ2()));
    //
    currTr->ResetCovMat();
    /*
    double *trCov  = (double*)currTr->GetCovariance();
    double *trPars = (double*)currTr->GetParameter();
    const double kLargeErr2PtI = 0.3*0.3;
    trCov[14] = TMath::Max(trCov[14],kLargeErr2PtI*trPars[4]*trPars[4]);
    */
  }
  //
  for (Int_t j=maxLr; j--; ) {  // Layer loop
    //
    KMCLayer *lrP = lr;
    lr = (KMCLayer*)fLayers.At(j);
    int ntPrev = lrP->GetNMCTracks();
    //
    if (lrP->IsDead()) { // for passive layer just propagate the copy of all tracks of prev layer >>>
      for (int itrP=ntPrev;itrP--;) { // loop over all tracks from previous layer
	currTrP = lrP->GetMCTrack(itrP); if (currTrP->IsKilled()) continue;
	currTr = lr->AddMCTrack( currTrP );
	if (!currTr->CorrectForMeanMaterial(lrP,kTRUE)) {currTr->Kill(); continue;} // correct for materials of prev. layer
	if (!PropagateToLayer(currTr,lr,-1))      {currTr->Kill(); continue;} // propagate to current layer
      }
      continue;
    } // treatment of dead layer <<<
    //
    if (lrP->IsTPC()) { // we don't consider bg hits in TPC, just update with MC cluster
      for (int itrP=ntPrev;itrP--;) { // loop over all tracks from previous layer
	currTrP = lrP->GetMCTrack(itrP); if (currTrP->IsKilled()) continue;
	currTr = lr->AddMCTrack( currTrP );
	if (!UpdateTrack(currTr, lrP, lrP->GetMCCluster(), kTRUE)) {currTr->Kill(); continue;} // update with correct MC cl.
	if (!currTr->CorrectForMeanMaterial(lrP,kTRUE)) {currTr->Kill(); continue;} // correct for materials of prev. layer
	if (!PropagateToLayer(currTr,lr,-1))      {currTr->Kill(); continue;} // propagate to current layer
      }
      continue;
    } // treatment of ideal (TPC?) layer <<<
    //
    // active layer under eff. study (ITS?): propagate copy of every track to MC cluster frame (to have them all in the same frame)
    // and calculate the limits of bg generation
    KMCCluster* clMC = lrP->GetMCCluster();
    if (lrP->GetLayerEff()<gRandom->Rndm()) clMC->Kill(); // simulate inefficiency
    ResetSearchLimits();
    int nseeds = 0;
    for (int itrP=ntPrev;itrP--;) { // loop over all tracks from previous layer
      currTrP = lrP->GetMCTrack(itrP); if (currTrP->IsKilled()) continue;
      currTr = lr->AddMCTrack( currTrP );
      if (!currTr->PropagateToCluster(clMC,fBFieldG)) {currTr->Kill(); continue;} // go to MC cluster
      if ( !(currTr->GetNITSHits()>0 && currTr->GetNITSHits()==currTr->GetNFakeITSHits()) ) UpdateSearchLimits(currTr, lrP); // RS
      nseeds++;
    }
    //
    //    printf("%3d seeds\n",nseeds);
    if (fUseBackground && lrP->IsITS()) GenBgClusters(lrP); //  generate background hits
    //
    ntPrev = lr->GetNMCTracks();
    for (int itr=ntPrev;itr--;) { // loop over all tracks PROPAGATED from previous layer to clusters frame on previous layer
      currTrP = lr->GetMCTrack(itr); // this is a seed from prev layer. The new clusters are attached to its copies, the seed itself
                                     // will be propagated w/o cluster update if it does not violate requested "reconstructed" track settings
      if (currTrP->IsKilled()) continue;
      //printf("Check    %d %p %d\n",itr,currTrP,currTrP->GetUniqueID()); currTrP->Print();
      CheckTrackProlongations(currTrP, lr, lrP);
      if (NeedToKill(currTrP)) currTrP->Kill(); // kill track which was not updated at lrP
      //currTrP->Kill(); // kill track which was not updated at lrP
    }
    //  
    lr->GetMCTracks()->Sort();
    int ntTot = lr->GetNMCTracks(); // propagate max amount of allowed tracks to current layer
    if (ntTot>fMaxSeedToPropagate && fMaxSeedToPropagate>0) {
      for (int itr=ntTot;itr>=fMaxSeedToPropagate;itr--)  lr->GetMCTracks()->RemoveAt(itr);
      ntTot = fMaxSeedToPropagate;
    }
    //
    for (int itr=ntTot;itr--;) {
      currTr = lr->GetMCTrack(itr);
      if (!currTr->CorrectForMeanMaterial(lrP,kTRUE)) {currTr->Kill();continue;} // correct for materials of prev. layer
      if (!PropagateToLayer(currTr,lr,-1))      {currTr->Kill();continue;} // propagate to current layer
    }
    AliDebug(1,Form("Got %d tracks on layer %s",ntTot,lr->GetName()));
    //    lr->GetMCTracks()->Print();
    //
  } // end loop over layers    
  //
  // do we use vertex constraint?
  KMCLayer *vtx = GetLayer(0);
  if (!vtx->IsDead() && vtx->IsITS()) {
    int ntr = vtx->GetNMCTracks();
    for (int itr=0;itr<ntr;itr++) {
      currTr = vtx->GetMCTrack(itr);
      if (currTr->IsKilled()) continue;
      KMCCluster* clv = vtx->GetMCCluster();
      double meas[2] = {clv->GetY(),clv->GetZ()};
      double measErr2[3] = {vtx->fPhiRes*vtx->fPhiRes,0,vtx->fZRes*vtx->fZRes};
      double chi2v = currTr->GetPredictedChi2(meas,measErr2);
      currTr->AddHit(vtx->GetActiveID(), chi2v, -1);
      currTr->SetInnerLrChecked(vtx->GetActiveID());
      if (NeedToKill(currTr)) currTr->Kill();
      // if (vtx->IsITS()) {if (!UpdateTrack(currTr, vtx, vtx->GetMCCluster(), kFALSE)) {currTr->Kill();continue;}}
    }
  }
  EliminateUnrelated();
  
  return kTRUE;
}

//____________________________________________________________________________
Bool_t KMCDetector::PropagateToLayer(KMCProbe* trc, KMCLayer* lr, int dir) const
{
  // bring the track to layer and rotat to frame normal to its surface
  if (!trc->PropagateToR(lr->fR,fBFieldG, dir)) return kFALSE;
  //
  // rotate to frame with X axis normal to the surface (defined by ideal track)
  if (!lr->IsVertex()) {
    double pos[3];
    trc->GetXYZ(pos);  // lab position
    double phi = TMath::ATan2(pos[1],pos[0]);
    if ( TMath::Abs(TMath::Abs(phi)-TMath::Pi()/2)<1e-3) phi = 0;
    if (!trc->Rotate(phi)) {
      AliDebug(2,Form("Failed to rotate to the frame (phi:%+.3f) of layer at %.2f at XYZ: %+.3f %+.3f %+.3f",
		      phi,lr->fR,pos[0],pos[1],pos[2]));
      if (AliLog::GetGlobalDebugLevel()>1) trc->Print("l");
      return kFALSE;
    }
  }
  //
  return kTRUE;
}

//____________________________________________________________________________
Bool_t KMCDetector::UpdateTrack(KMCProbe* trc, KMCLayer* lr, KMCCluster* cl, Bool_t goToCluster) const
{
  // update track with measured cluster
  // propagate to cluster
  double meas[2] = {cl->GetY(),cl->GetZ()}; // ideal cluster coordinate
  double measErr2[3] = {lr->fPhiRes*lr->fPhiRes,0,lr->fZRes*lr->fZRes};
  //
  if (goToCluster) if (!trc->PropagateToCluster(cl,fBFieldG)) return kFALSE; // track was not propagated to cluster frame
  //
  double chi2 = trc->GetPredictedChi2(meas,measErr2);
  //  if (chi2>fMaxChi2Cl) return kTRUE; // chi2 is too large
  //  
  if (!trc->Update(meas,measErr2)) {
    AliDebug(2,Form("layer %s: Failed to update the track by measurement {%.3f,%3f} err {%.3e %.3e %.3e}",
		    lr->GetName(),meas[0],meas[1], measErr2[0],measErr2[1],measErr2[2]));
    if (AliLog::GetGlobalDebugLevel()>1) trc->Print("l");
    return kFALSE;
  }
  trc->AddHit(lr->GetActiveID(), chi2);
  //
  return kTRUE;
}

//____________________________________________________________________________
Int_t KMCDetector::GenBgClusters(KMCLayer* lr)
{
  // Generate fake clusters in precalculated RPhi,Z range
  if (fNBgLimits<1) return 0; // limits were not set - no track was prolongated
  //
  // Fix search limits to avoid seeds which will anyway point very far from the vertex
  double tolY = TMath::Sqrt(lr->fSig2EstD)*fMaxChi2ClSQ;
  double tolZ = TMath::Sqrt(lr->fSig2EstZ)*fMaxChi2ClSQ;
  
  //  printf("Before: Y: %+6.3f : %+6.3f tolY: %6.3f || Z: %+6.3f : %+6.3f tolZ: %6.3f\n",fBgYMin,fBgYMax,tolY, fBgZMin,fBgZMax,tolZ);
  if (fBgYMin < lr->fClCorr.fY-tolY) fBgYMin = lr->fClCorr.fY-tolY;
  if (fBgYMax > lr->fClCorr.fY+tolY) fBgYMax = lr->fClCorr.fY+tolY;
  if (fBgZMin < lr->fClCorr.fZ-tolZ) fBgZMin = lr->fClCorr.fZ-tolZ;
  if (fBgZMax > lr->fClCorr.fZ+tolZ) fBgZMax = lr->fClCorr.fZ+tolZ;
  //printf("After: Y: %+6.3f : %+6.3f tolY: %6.3f || Z: %+6.3f : %+6.3f tolZ: %6.3f\n",fBgYMin,fBgYMax,tolY, fBgZMin,fBgZMax,tolZ);
  //
  double dy = fBgYMax - fBgYMin;
  double dz = fBgZMax - fBgZMin;
  double surf = dy*dz;               // surface of generation
  if (surf<0) return 0;
  double poissProb = surf*HitDensity(lr->fR)*lr->GetLayerEff();
  AliDebug(2,Form("Bg for Lr %s (r=%.2f) : Density %.2f on surface %.2e [%+.4f : %+.4f][%+.4f %+.4f]",
		  lr->GetName(),lr->fR,HitDensity(lr->fR),surf,fBgYMin,fBgYMax,fBgZMin,fBgZMax));
  int nFakesGen = gRandom->Poisson( poissProb ); // preliminary number of extra clusters to test
  KMCCluster *refCl = lr->GetMCCluster();
  double sig2y = lr->GetPhiRes()*lr->GetPhiRes();
  double sig2z = lr->GetZRes()*lr->GetZRes();
  for (int ic=nFakesGen;ic--;) {
    double y = fBgYMin+dy*gRandom->Rndm();
    double z = fBgZMin+dz*gRandom->Rndm();
    double dfy = y-refCl->GetY();
    double dfz = z-refCl->GetZ();
    double dist = (dfy*dfy)/sig2y + (dfz*dfz)/sig2z;
    if (dist<4) continue; // avoid overlap with MC cluster
    lr->AddBgCluster(y, z, refCl->GetX(), refCl->GetPhi());
  }
  AliDebug(2,Form("Added %6d noise clusters on lr %s (poisson Prob=%8.2f for surface %.2e) DY:%7.4f DZ: %7.4f",
		   lr->GetNBgClusters(),lr->GetName(),poissProb,surf,dy,dz));
  return nFakesGen;
  //
}

//____________________________________________________________________________
void KMCDetector::UpdateSearchLimits(KMCProbe* probe, KMCLayer* lr)
{
  // define the search window for track on layer (where the bg hist will be generated)
  static double *currYMin = fBgYMinTr.GetArray();
  static double *currYMax = fBgYMaxTr.GetArray();
  static double *currZMin = fBgZMinTr.GetArray();
  static double *currZMax = fBgZMaxTr.GetArray();
  //
  double sizeY = probe->GetSigmaY2(), sizeZ = probe->GetSigmaZ2();
  //  /*
  if (probe->GetNITSHits()<3) sizeY = 10*lr->fSig2EstD;
  if (probe->GetNITSHits()<2) sizeZ = 10*lr->fSig2EstZ;
  sizeY = fMaxChi2ClSQ*TMath::Sqrt(sizeY+lr->fPhiRes*lr->fPhiRes); // max deviation in rphi to accept
  sizeZ = fMaxChi2ClSQ*TMath::Sqrt(sizeZ+lr->fZRes*lr->fZRes); // max deviation in dz to accept
  //  */
  //
  /*
  if (probe->GetNITSHits()<3) sizeY = 1./(1./lr->fSig2EstD + 1./sizeY);
  if (probe->GetNITSHits()<2) sizeZ = 1./(1./lr->fSig2EstZ + 1./sizeZ);
  sizeY = fMaxChi2ClSQ*TMath::Sqrt(sizeY); // max deviation in rphi to accept
  sizeZ = fMaxChi2ClSQ*TMath::Sqrt(sizeZ); // max deviation in dz to accept
  */
  //
  //  if (sizeY>2) sizeY=2;
  //  if (sizeZ>2) sizeZ=2;
  //  printf("Sizes at %s: %.5f %.5f\n",lr->GetName(), sizeY,sizeZ);
  //
  if (fNBgLimits>=fBgYMinTr.GetSize()) { // expand arrays, update pointers
    fBgYMinTr.Set(2*(fNBgLimits+1));
    fBgYMaxTr.Set(2*(fNBgLimits+1));
    fBgZMinTr.Set(2*(fNBgLimits+1));
    fBgZMaxTr.Set(2*(fNBgLimits+1));
    currYMin = fBgYMinTr.GetArray();
    currYMax = fBgYMaxTr.GetArray();
    currZMin = fBgZMinTr.GetArray();
    currZMax = fBgZMaxTr.GetArray();
  }
  if (fBgYMin > (currYMin[fNBgLimits]=probe->GetY()-sizeY) ) fBgYMin = currYMin[fNBgLimits];
  if (fBgYMax < (currYMax[fNBgLimits]=probe->GetY()+sizeY) ) fBgYMax = currYMax[fNBgLimits];
  if (fBgZMin > (currZMin[fNBgLimits]=probe->GetZ()-sizeZ) ) fBgZMin = currZMin[fNBgLimits];
  if (fBgZMax < (currZMax[fNBgLimits]=probe->GetZ()+sizeZ) ) fBgZMax = currZMax[fNBgLimits];
  if (AliLog::GetGlobalDebugLevel()>=2) {
    probe->Print("l");
    AliInfo(Form("Seed%3d Lr %s limits for y:%+8.4f z:%+8.4f [%+.4f : %+.4f][%+.4f %+.4f]",fNBgLimits,lr->GetName(),probe->GetY(),probe->GetZ(),currYMin[fNBgLimits],currYMax[fNBgLimits],currZMin[fNBgLimits],currZMax[fNBgLimits]));
    AliInfo(Form("Global Limits Lr %s                            [%+.4f : %+.4f][%+.4f %+.4f]",lr->GetName(),fBgYMin,fBgYMax,fBgZMin,fBgZMax));
    AliInfo(Form("MC Cluster: %+.4f : %+.4f",lr->fClMC.fY, lr->fClMC.fZ));
  }
  probe->SetUniqueID(fNBgLimits++);
  //
  if (lr->IsITS() && probe->GetNFakeITSHits()==0) {
    if (fHMCLrResidRPhi) fHMCLrResidRPhi->Fill(probe->GetY() - lr->GetMCCluster()->GetY(), lr->GetActiveID());
    if (fHMCLrResidZ)    fHMCLrResidZ->Fill(probe->GetZ() - lr->GetMCCluster()->GetZ(),lr->GetActiveID());
  }
  //
}

//____________________________________________________________________________
void KMCDetector::CheckTrackProlongations(KMCProbe *probe, KMCLayer* lr, KMCLayer* lrP)
{
  // explore prolongation of probe from lrP to lr with all possible clusters of lrP
  // the probe is already brought to clusters frame
  int nCl = lrP->GetNBgClusters();
  double measErr2[3] = {lrP->fPhiRes*lrP->fPhiRes,0,lrP->fZRes*lrP->fZRes};
  double meas[2] = {0,0};
  UInt_t tmpID = probe->GetUniqueID();
  double yMin = fBgYMinTr[tmpID];
  double yMax = fBgYMaxTr[tmpID];
  double zMin = fBgZMinTr[tmpID];
  double zMax = fBgZMaxTr[tmpID];
  //
  probe->SetInnerLrChecked(lrP->GetActiveID());
  for (int icl=-1;icl<nCl;icl++) {
    KMCCluster* cl = icl<0 ? lrP->GetMCCluster() : lrP->GetBgCluster(icl);  // -1 is for true MC cluster
    if (cl->IsKilled()) {
      if (AliLog::GetGlobalDebugLevel()>1) {printf("Skip cluster %d ",icl); cl->Print();}
      continue;
    }
    double y = cl->GetY();
    double z = cl->GetZ();
    AliDebug(2,Form("Check seed%d against cl#%d out of %d at layer %s | y:%+8.4f z:%+8.4f [%+.4f:%+.4f]  [%+.4f:%+.4f]",tmpID,icl,nCl,lrP->GetName(),y,z,yMin,yMax,zMin,zMax));
    if (AliLog::GetGlobalDebugLevel()>0) {
      if (icl==-1 && probe->GetNFakeITSHits()==0) {
	meas[0] = y; meas[1] = z;
	double chi2a = probe->GetPredictedChi2(meas,measErr2);
	if (chi2a>fMaxChi2Cl || (y<yMin || y>yMax) || (z<zMin || z>zMax)) {
	  probe->Print();
	  printf("Loosing good point (y:%+8.4f z:%+8.4f) on lr %s: chi2: %.2f  | dy:%+8.4f dz:%+8.4f [%+.4f:%+.4f]  [%+.4f:%+.4f] |x: %.2f %.2f | phi: %.2f %.2f\n",
		 y,z,lrP->GetName(),chi2a,y-probe->GetY(),z-probe->GetZ(),yMin,yMax,zMin,zMax, probe->GetX(), cl->GetX(), probe->GetAlpha(), cl->GetPhi());
	}
      }
    }
    if (y<yMin || y>yMax) continue; // preliminary check on Y
    if (z<zMin || z>zMax) continue; // preliminary check on Z
    meas[0] = y; meas[1] = z;
    double chi2 = probe->GetPredictedChi2(meas,measErr2);
    if (fHMCLrChi2 && probe->GetNFakeITSHits()==0 && icl==-1) fHMCLrChi2->Fill(chi2,lrP->GetActiveID());
    AliDebug(2,Form("Seed-to-cluster chi2 = Chi2=%.2f",chi2));
    if (chi2>fMaxChi2Cl) continue;
    // 
    // update track copy
    KMCProbe* newTr = lr->AddMCTrack( probe );
    fUpdCalls++;
    if (!newTr->Update(meas,measErr2)) {
      AliDebug(2,Form("Layer %s: Failed to update the track by measurement {%.3f,%3f} err {%.3e %.3e %.3e}",
		      lrP->GetName(),meas[0],meas[1], measErr2[0],measErr2[1],measErr2[2]));
      if (AliLog::GetGlobalDebugLevel()>1) newTr->Print("l");
      newTr->Kill();
      continue;
    }
    newTr->AddHit(lrP->GetActiveID(), chi2, icl);
    if (AliLog::GetGlobalDebugLevel()>1) {
      AliInfo("Cloned updated track is:");
      newTr->Print();
    }
    if (NeedToKill(newTr)) newTr->Kill();
  }
  //
}

//____________________________________________________________________________
void KMCDetector::ResetMCTracks(Int_t maxLr)
{
  int nl = GetNLayers();
  if (maxLr<0 || maxLr>=nl) maxLr = nl-1;
  for (int i=maxLr+1;i--;) GetLayer(i)->ResetMCTracks();
}

//____________________________________________________________________________
Bool_t KMCDetector::NeedToKill(KMCProbe* probe) const
{
  // check if the seed at given layer (last one where update was tried) 
  // still has chances to be reconstructed
  const Bool_t kModeKillMiss = kFALSE;
  //
  Bool_t kill = kFALSE;
  while (1) {
    int il = probe->GetInnerLayerChecked();
    int nITS = probe->GetNITSHits();
    int nITSMax = nITS + il; // maximum it can have
    if (nITSMax<fMinITSHits) {
      kill = kTRUE; 
      break;
    }    // has no chance to collect enough ITS hits
    //
    int ngr = fPattITS.GetSize();
    if (ngr>0) { // check pattern
      UInt_t patt = probe->GetHitsPatt();
      // complete the layers not checked yet
      for (int i=il;i--;) patt |= (0x1<<i);
      for (int ig=ngr;ig--;) 
	if (!(((UInt_t)fPattITS[ig]) & patt)) {
	  kill = kTRUE; 
	  break;
	}
      //
    }
    //
    if (nITS>2) {  // check if smallest possible norm chi2/ndf is acceptable
      double chi2min = probe->GetChi2();
      if (kModeKillMiss) {
	int nMiss = fNActiveITSLayers - probe->GetInnerLayerChecked() - nITS; // layers already missed
	chi2min = nMiss*probe->GetMissingHitPenalty();
      }
      chi2min /= ((nITSMax<<1)-KMCProbe::kNDOF);
      if (chi2min>fMaxNormChi2NDF) {
	kill = kTRUE; 
	break;
      }
    }
    //
    // loose vertex constraint
    double dst;
    if (nITS>=2) {
      probe->GetZAt(0,fBFieldG,dst);
      //printf("Zd (F%d): %f\n",probe->GetNFakeITSHits(),dst);
      if (TMath::Abs(dst)>10.) {
	kill = kTRUE; 
	break;
      }
    }
    if (nITS>=3) {
      probe->GetYAt(0,fBFieldG,dst);
      //printf("Dd (F%d): %f\n",probe->GetNFakeITSHits(),dst);
      if (TMath::Abs(dst)>10.) {
	kill = kTRUE; 
	break;
      }
    }
    //
    break;
  }
  if (kill && AliLog::GetGlobalDebugLevel()>1 && probe->GetNFakeITSHits()==0) {
    printf("Killing good seed, last upd layer was %d\n",probe->GetInnerLayerChecked());
    probe->Print("l");
  }
  return kill;
}

//_____________________________________________________________________
void KMCDetector::EliminateUnrelated()
{
  // kill useless tracks
  KMCLayer* lr = GetLayer(0);
  int ntr = lr->GetNMCTracks();
  int nval = 0;
  for (int itr=0;itr<ntr;itr++) {
    KMCProbe* probe = lr->GetMCTrack(itr);
    if (probe->IsKilled()) continue;
    if (probe->GetNITSHits()-probe->GetNFakeITSHits()<1) {probe->Kill(); continue;}
    nval++;
  }
  lr->GetMCTracks()->Sort();
  const int kDump = 0;
  if (kDump>0) {
    printf("Valid %d out of %d\n",nval, ntr);
    ntr = ntr>kDump ? kDump:0;
    for (int itr=0;itr<ntr;itr++) {
      lr->GetMCTrack(itr)->Print();
    }
  }
}

//_____________________________________________________________________
void KMCDetector::RequirePattern(UInt_t *patt, int groups)
{
  if (groups<1) {fPattITS.Set(0); return;}
  fPattITS.Set(groups);
  for (int i=0;i<groups;i++) fPattITS[i] = patt[i];
}

//_____________________________________________________________________
void KMCDetector::CalcHardSearchLimits(double dzv)
{
  //
  TArrayD zlims;
  zlims.Set(fNActiveITSLayers);
  for (int il0=0;il0<fNActiveITSLayers;il0++) {
    KMCLayer* lr0 = GetActiveLayer(il0);
    double angZTol = dzv/lr0->GetRadius();
    for (int il1=0;il1<fNActiveITSLayers;il1++) {
      if (il1==il0) continue;
      KMCLayer* lr1 = GetActiveLayer(il1);
      double ztol = angZTol*TMath::Abs(lr0->GetRadius() - lr1->GetRadius());
      if (ztol>zlims[il1]) zlims[il1] = ztol;
    }
  }
  //
  for (int il=0;il<fNActiveITSLayers;il++) printf("ZTol%d: %8.4f\n",il,zlims[il]);
}

//_______________________________________________________
double KMCDetector::PropagateBack(KMCProbe* trc) 
{
  static KMCProbe bwd;
  bwd = *trc;
  bwd.ResetCovMat();
  static double measErr2[3] = {0,0,0};
  static double meas[2] = {0,0};
  int icl = 0;
  double chi2Tot = 0;
  for (int il=1;il<=fLastActiveITSLayer;il++) {
    KMCLayer* lr = GetLayer(il);
    if (!PropagateToLayer(&bwd,lr,1)) return -1;
    int aID = lr->GetActiveID();
    if (aID>-1 && (icl=bwd.fClID[aID])>=-1) {
      KMCCluster* clMC =  icl<0 ? lr->GetMCCluster() : lr->GetBgCluster(icl);
      if (!bwd.PropagateToCluster(clMC,fBFieldG)) return -1;
      meas[0] = clMC->GetY(); meas[1] = clMC->GetZ();
      measErr2[0] = lr->fPhiRes*lr->fPhiRes;
      measErr2[2] = lr->fZRes*lr->fZRes;
      double chi2a = bwd.GetPredictedChi2(meas,measErr2);
      chi2Tot += chi2a;
      printf("Chis %d (cl%+3d): t2c: %6.3f tot: %6.3f\n",aID,icl,chi2a, chi2Tot);
      bwd.Update(meas,measErr2);
      bwd.AddHit(aID, chi2a, icl);	
    }
    if (!bwd.CorrectForMeanMaterial(lr,kFALSE)) return -1;
  }
  return chi2Tot;
}


//////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------------------------
ClassImp(KMCTrackSummary)

Int_t KMCTrackSummary::fgSumCounter = 0;

//____________________________________________________________________________
KMCTrackSummary::KMCTrackSummary(const char* name,const char* title, int nlr) :
  TNamed(name,name), fNITSLayers(nlr),
  fPattITS(0),fPattITSFakeExcl(0),fPattITSCorr(0),
  fHMCChi2(0),fHMCSigDCARPhi(0),fHMCSigDCAZ(0),fHMCSigPt(0),fHMCNCl(0),fHMCFakePatt(0),
  fCountAcc(0),fCountTot(0),fUpdCalls(0),
  fRefProbe(),fAnProbe()
{
  // create summary structure for specific track conditions
  //
  SetMinMaxClITS();
  SetMinMaxClITSFake();
  SetMinMaxClITSCorr();
  //
  if (name==0 && title==0 && nlr==0) return;  // default dummy contructor
  //
  enum {kNBinRes=1000};
  const double kMaxChi2(10),kMaxResRPhi=0.1, kMaxResZ=0.1,kMaxResPt=0.3;
  //
  TString strN = name, strT = title;
  if (strN.IsNull()) strN = Form("TrSum%d",fgSumCounter);
  if (strT.IsNull()) strT = strN;
  fgSumCounter++;
  //
  if (fNITSLayers<1) {
    fNITSLayers=KMCProbe::GetNITSLayers(); 
    if (fNITSLayers<1) {AliError("N ITS layer is not provided and not available from KMCProbe::GetNITSLayers()"); exit(1);}
    AliInfo(Form("%s No N Layers provided, set %d",strN.Data(),fNITSLayers));
  }
  nlr = fNITSLayers;
  //
  TString nam,tit;
  //
  nam = Form("%s_mc_chi2",strN.Data());
  tit = Form("%s MC #chi^{2}",strT.Data());
  fHMCChi2 = new TH1F(nam.Data(),tit.Data(),kNBinRes,0,kMaxChi2);
  fHMCChi2->GetXaxis()->SetTitle("#chi^{2}");
  fHMCChi2->Sumw2();
  //
  nam = Form("%s_mc_DCArphi",strN.Data());
  tit = Form("%s MC #DeltaDCA r#phi",strT.Data());
  fHMCSigDCARPhi = new TH1F(nam.Data(),tit.Data(),kNBinRes,-kMaxResRPhi,kMaxResRPhi);
  fHMCSigDCARPhi->GetXaxis()->SetTitle("#Delta r#phi");
  fHMCSigDCARPhi->Sumw2();
  //
  nam = Form("%s_mc_DCAz",strN.Data());
  tit = Form("%s MC #DeltaDCA Z",strT.Data());
  fHMCSigDCAZ = new TH1F(nam.Data(),tit.Data(),kNBinRes,-kMaxResZ,kMaxResZ);
  fHMCSigDCAZ->GetXaxis()->SetTitle("#Delta Z");
  fHMCSigDCAZ->Sumw2();
  //
  nam = Form("%s_mc_pt",strN.Data());
  tit = Form("%s MC $Deltap_{T}/p_{T}",strT.Data());
  fHMCSigPt = new TH1F(nam.Data(),tit.Data(),2*kNBinRes,-kMaxResPt,kMaxResPt);
  fHMCSigPt->GetXaxis()->SetTitle("#Deltap_{T}/p_{T}");
  fHMCSigPt->Sumw2();
  //
  nam = Form("%s_n_cl",strN.Data());
  tit = Form("%s N Clusters",strT.Data());
  fHMCNCl = new TH2F(nam.Data(),tit.Data(),nlr,1,nlr+1,nlr,1,nlr+1);
  fHMCNCl->GetXaxis()->SetTitle("N Clusters Tot");
  fHMCNCl->GetYaxis()->SetTitle("N Clusters Fake");
  fHMCNCl->Sumw2();
  //
  nam = Form("%s_fake_cl",strN.Data());
  tit = Form("%s Fake Clusters",strT.Data());
  fHMCFakePatt = new TH1F(nam.Data(),tit.Data(),nlr,-0.5,nlr-0.5);
  fHMCFakePatt->GetXaxis()->SetTitle("Fake clusters pattern");
  fHMCFakePatt->Sumw2();
  //
}

//____________________________________________________________________________
KMCTrackSummary::~KMCTrackSummary()
{
  if (fHMCChi2)         delete fHMCChi2;
  if (fHMCSigDCARPhi)   delete fHMCSigDCARPhi;
  if (fHMCSigDCAZ)      delete fHMCSigDCAZ;
  if (fHMCSigPt)        delete fHMCSigPt;
  if (fHMCNCl)          delete fHMCNCl;
  if (fHMCFakePatt)     delete fHMCFakePatt;
}


//____________________________________________________________________________
void KMCTrackSummary::Add(const KMCTrackSummary* src, double scl)
{
  if (fHMCChi2)         fHMCChi2->Add(src->fHMCChi2,scl);
  if (fHMCSigDCARPhi)   fHMCSigDCARPhi->Add(src->fHMCSigDCARPhi,scl);
  if (fHMCSigDCAZ)      fHMCSigDCAZ->Add(src->fHMCSigDCAZ,scl);
  if (fHMCSigPt)        fHMCSigPt->Add(src->fHMCSigPt,scl);
  if (fHMCNCl)          fHMCNCl->Add(src->fHMCNCl,scl);
  if (fHMCFakePatt)     fHMCFakePatt->Add(src->fHMCFakePatt,scl);
  fCountAcc += src->fCountAcc;
  fUpdCalls += src->fUpdCalls;
}

//____________________________________________________________________________
void KMCTrackSummary::PrependTNamed(TNamed* obj, const char* nm, const char* tit)
{
  // prepend  name, title by this prefix
  TString name = obj->GetName(),title = obj->GetTitle();
  name.Prepend(nm); title.Prepend(tit);
  obj->SetNameTitle(name.Data(),title.Data());
  //
}
  
//____________________________________________________________________________
void KMCTrackSummary::SetNamePrefix(const char* pref)
{
  // prepend all names by this prefix
  TString nmT,nmP = pref;
  if (nmP.IsNull()) return;
  nmP = Form("%s_",pref);
  nmT = Form("%s ",pref);
  PrependTNamed(this, nmP.Data(), nmT.Data());
  PrependTNamed(fHMCChi2,       nmP.Data(), nmT.Data());
  PrependTNamed(fHMCSigDCARPhi, nmP.Data(), nmT.Data());
  PrependTNamed(fHMCSigDCAZ,    nmP.Data(), nmT.Data());
  PrependTNamed(fHMCSigPt,      nmP.Data(), nmT.Data());
  PrependTNamed(fHMCNCl,        nmP.Data(), nmT.Data());
  PrependTNamed(fHMCFakePatt,   nmP.Data(), nmT.Data());
  //  
}

//____________________________________________________________________________
Bool_t KMCTrackSummary::CheckTrack(KMCProbe* trc)
{
  // check if the track satisfies to selection conditions
  fCountTot++;
  if (!trc) return kFALSE;
  if (OutOfRange(trc->GetNITSHits(), fMinMaxClITS)) return kFALSE;
  if (OutOfRange(trc->GetNFakeITSHits(), fMinMaxClITSFake)) return kFALSE;
  if (OutOfRange(trc->GetNITSHits()-trc->GetNFakeITSHits(), fMinMaxClITSCorr)) return kFALSE;
  //
  // check layer patterns if requested
  UInt_t patt  = trc->GetHitsPatt();
  UInt_t pattF = trc->GetFakesPatt();
  UInt_t pattC = patt&(~pattF);    // correct hits
  // is there at least one hit in each of requested layer groups?
  for (int ip=fPattITS.GetSize();ip--;)     if (!CheckPattern(patt,fPattITS[ip])) return kFALSE;
  // is there at least one fake in any of requested layer groups?
  for (int ip=fPattITSFakeExcl.GetSize();ip--;) if (CheckPattern(pattF,fPattITSFakeExcl[ip])) return kFALSE;
  // is there at least one correct hit in each of requested layer groups?
  for (int ip=fPattITSCorr.GetSize();ip--;) if (!CheckPattern(pattC,fPattITSCorr[ip])) return kFALSE;
  //
  fCountAcc++;
  return kTRUE;
}

//____________________________________________________________________________
void KMCTrackSummary::AddTrack(KMCProbe* trc)
{
  // fill track info
  if (!CheckTrack(trc)) return;
  //
  fHMCChi2->Fill(trc->GetNormChi2(kFALSE));
  fHMCSigDCARPhi->Fill(trc->GetY());
  fHMCSigDCAZ->Fill(trc->GetZ());
  if (fRefProbe.Pt()>0) fHMCSigPt->Fill( trc->Pt()/fRefProbe.Pt()-1.);
  //  printf("Pts: %.3f %.3f -> %.3f\n",trc->Pt(),fRefProbe.Pt(),trc->Pt()/fRefProbe.Pt()-1.);
  //  trc->Print("l");
  //  fRefProbe.Print("l");
  fHMCNCl->Fill(trc->GetNITSHits(),trc->GetNFakeITSHits());
  for (int i=trc->GetNITSLayers();i--;) if (trc->IsHitFake(i)) fHMCFakePatt->Fill(i);
  //
}

//____________________________________________________________________________
UInt_t KMCTrackSummary::Bits(Bool_t l0, Bool_t l1, Bool_t l2, Bool_t l3, Bool_t l4, Bool_t l5, Bool_t l6, Bool_t l7,Bool_t  l8, Bool_t l9,
			  Bool_t l10,Bool_t l11,Bool_t l12,Bool_t l13,Bool_t l14,Bool_t l15,Bool_t l16,Bool_t l17,Bool_t l18,Bool_t l19,
			  Bool_t l20,Bool_t l21,Bool_t l22,Bool_t l23,Bool_t l24,Bool_t l25,Bool_t l26,Bool_t l27,Bool_t l28,Bool_t l29,
			  Bool_t l30,Bool_t l31)
{
  // create corresponding bit pattern
  UInt_t patt = 0;
  if (l0 ) patt |= (0x1<<0);    if (l1 ) patt |= (0x1<<1);   if (l2 ) patt |= (0x1<<2);
  if (l3 ) patt |= (0x1<<3);    if (l4 ) patt |= (0x1<<4);   if (l5 ) patt |= (0x1<<5);
  if (l6 ) patt |= (0x1<<6);    if (l7 ) patt |= (0x1<<7);   if (l8 ) patt |= (0x1<<8);
  if (l9 ) patt |= (0x1<<9);    if (l10) patt |= (0x1<<10);  if (l11) patt |= (0x1<<11);
  if (l12) patt |= (0x1<<12);   if (l13) patt |= (0x1<<13);  if (l14) patt |= (0x1<<14);
  if (l15) patt |= (0x1<<15);   if (l16) patt |= (0x1<<16);  if (l17) patt |= (0x1<<17);
  if (l18) patt |= (0x1<<18);   if (l19) patt |= (0x1<<19);  if (l20) patt |= (0x1<<20);
  if (l21) patt |= (0x1<<21);   if (l22) patt |= (0x1<<22);  if (l23) patt |= (0x1<<23);
  if (l24) patt |= (0x1<<24);   if (l25) patt |= (0x1<<25);  if (l26) patt |= (0x1<<26);
  if (l27) patt |= (0x1<<27);   if (l28) patt |= (0x1<<28);  if (l29) patt |= (0x1<<29);
  if (l30) patt |= (0x1<<30);   if (l31) patt |= (0x1<<31);   
  return patt;
}

//__________________________________________________________________________
void KMCTrackSummary::Print(Option_t* ) const
{
  printf("%s: summary for track M=%5.3f pT: %6.3f eta: %.2f\n", GetName(),
	 fRefProbe.GetMass(),fRefProbe.Pt(), fRefProbe.Eta());
  printf("Cuts on NCl ITS: Tot: %2d - %2d Fake: %2d - %2d Corr: %2d - %2d\n",
	 fMinMaxClITS[0],fMinMaxClITS[1], 
	 fMinMaxClITSFake[0],fMinMaxClITSFake[1],
	 fMinMaxClITSCorr[0],fMinMaxClITSCorr[1]);
  //
  int nlr = fNITSLayers;
  if (fPattITS.GetSize()) {
    printf("Require at least 1 hit in groups: ");
    printf("Hits obligatory in groups: ");
    for (int i=fPattITS.GetSize();i--;) {
      UInt_t pat = (UInt_t)fPattITS[i];
      printf("[");
      for (int j=0;j<nlr;j++) if (pat&(0x1<<j)) printf("%d ",j);
      printf("] ");
    }
    printf("\n");
  }
  //
  if (fPattITSFakeExcl.GetSize()) {
    printf("Fake hit forbidden in groups    : ");
    for (int i=fPattITSFakeExcl.GetSize();i--;) {
      UInt_t pat = (UInt_t)fPattITSFakeExcl[i];
      printf("[");
      for (int j=0;j<nlr;j++) if (pat&(0x1<<j)) printf("%d ",j);
      printf("] ");
    }
    printf("\n");
  }
  //
  if (fPattITSCorr.GetSize()) {
    printf("Correct hit obligatory in groups: ");
    for (int i=fPattITSCorr.GetSize();i--;) {
      UInt_t pat = (UInt_t)fPattITSCorr[i];
      printf("[");
      for (int j=0;j<nlr;j++) if (pat&(0x1<<j)) printf("%d ",j);
      printf("] ");
    }
    printf("\n");
  }
  //  
  printf("Entries: Tot: %.4e Acc: %.4e -> Eff: %.4f+-%.4f %.2e KMC updates\n",fCountTot,fCountAcc,GetEff(),GetEffErr(),GetUpdCalls());
}

