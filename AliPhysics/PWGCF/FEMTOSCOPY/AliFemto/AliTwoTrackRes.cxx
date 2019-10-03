////////////////////////////////////////////////////////////////////////////////
//345678901234567890123456789012345678901234567890123456789012345678901234567890
//       1         2         3         4         5         6         7         8
//
// Tool to study two-track effects in ALICE for femtoscopic analyses
// J. Mercado <mercado@physi.uni-heidelberg.de> Last modified: 20.01.2011
//
////////////////////////////////////////////////////////////////////////////////

#include "AliTwoTrackRes.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TLeaf.h"
#include "TNtuple.h"
#include "TRandom2.h"

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliTwoTrackRes);
  /// \endcond
#endif

//______________________________________________________________________________
// Constructor(s)

AliTwoTrackRes::AliTwoTrackRes(const char *name) :
  AliAnalysisTask(name,""), fChain(0), fESDEvent(0),
  fOutContainer(0), fTrackCuts(0), fNTuple1(0),
  fNTuple2(0), fP1(), fP2(), fPb1(), fPb2(), fP(), fQ(), fTpcEnt1(), fTpcEnt2(),
  fTpcDist(), fOutFilename()
{
  DefineInput(0, TChain::Class());     // Slot input 0 reads from a TChain
  DefineOutput(0, TObjArray::Class()); // Slot output 0 writes into a TObjArray
}

AliTwoTrackRes::AliTwoTrackRes(const AliTwoTrackRes& aTwoTrackRes) :
  AliAnalysisTask(aTwoTrackRes), fChain(0), fESDEvent(0), fOutContainer(0),
  fTrackCuts(0), fNTuple1(0), fNTuple2(0), fP1(), fP2(), fPb1(), fPb2(), fP(),
  fQ(), fTpcEnt1(), fTpcEnt2(), fTpcDist(), fOutFilename()
{
  //Copy constructor
  fChain = aTwoTrackRes.fChain;
  fESDEvent = aTwoTrackRes.fESDEvent;
  fOutContainer = aTwoTrackRes.fOutContainer;
  fTrackCuts = aTwoTrackRes.fTrackCuts;
  fNTuple1 = aTwoTrackRes.fNTuple1;
  fNTuple2 = aTwoTrackRes.fNTuple2;
  fP1 = aTwoTrackRes.fP1;
  fP2 = aTwoTrackRes.fP2;
  fPb1 = aTwoTrackRes.fPb1;
  fPb2 = aTwoTrackRes.fPb2;
  fP = aTwoTrackRes.fP;
  fQ = aTwoTrackRes.fQ;
  fTpcEnt1 = aTwoTrackRes.fTpcEnt1;
  fTpcEnt2 = aTwoTrackRes.fTpcEnt2;
  fTpcDist = aTwoTrackRes.fTpcDist;
  fOutFilename = aTwoTrackRes.fOutFilename;
}

AliTwoTrackRes& AliTwoTrackRes::operator=(const AliTwoTrackRes& aTwoTrackRes)
{
  // Assignment operator
  if (this == &aTwoTrackRes)
    return *this;
  fChain = aTwoTrackRes.fChain;
  fESDEvent = aTwoTrackRes.fESDEvent;
  fOutContainer = aTwoTrackRes.fOutContainer;
  fTrackCuts = aTwoTrackRes.fTrackCuts;
  fNTuple1 = aTwoTrackRes.fNTuple1;
  fNTuple2 = aTwoTrackRes.fNTuple2;
  fP1 = aTwoTrackRes.fP1;
  fP2 = aTwoTrackRes.fP2;
  fPb1 = aTwoTrackRes.fPb1;
  fPb2 = aTwoTrackRes.fPb2;
  fP = aTwoTrackRes.fP;
  fQ = aTwoTrackRes.fQ;
  fTpcEnt1 = aTwoTrackRes.fTpcEnt1;
  fTpcEnt2 = aTwoTrackRes.fTpcEnt2;
  fTpcDist = aTwoTrackRes.fTpcDist;
  fOutFilename = aTwoTrackRes.fOutFilename;
  return *this;
}

AliTwoTrackRes::~AliTwoTrackRes() {printf("AliTwoTrackRes destroyed\n");}

void AliTwoTrackRes::ConnectInputData(Option_t *) {
//______________________________________________________________________________
// Connect input data and initialize track cuts

  fChain = (TChain*)GetInputData(0);
  fESDEvent = new AliESDEvent();
  fESDEvent->ReadFromTree(fChain);

  // Cuts to select primary tracks (ITS+TPC)
  fTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
  Double_t cov1, cov2, cov3, cov4, cov5; // diagonal cov. matrix elements
  Double_t nSigma;                       // max. DCA to primary vertex
  Int_t minNClustersTPC;                 // min. number of clusters per TPC tracks
  Double_t maxChi2PerClusterTPC;         // max. chi2 per cluster per TPC track
  Int_t cutMode = 1;                     // select cut mode
  if (cutMode == 1) {
  cov1 = 2; cov2 = 2; cov3 = 0.5; cov4 = 0.5; cov5 = 2;
  nSigma = 3; minNClustersTPC = 75; maxChi2PerClusterTPC = 3.5;
  fTrackCuts->SetMaxCovDiagonalElements(cov1, cov2, cov3, cov4, cov5);
  fTrackCuts->SetMaxNsigmaToVertex(nSigma);
  fTrackCuts->SetRequireSigmaToVertex(kTRUE);
  fTrackCuts->SetRequireTPCRefit(kTRUE);
  fTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts->SetMinNClustersTPC(minNClustersTPC);
  fTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
  TString tag("Global tracking");}
}

void AliTwoTrackRes::CreateOutputObjects() {
//______________________________________________________________________________
// Create output objects

  fNTuple1 = new TNtuple("nt1","True pairs",
  "pt1:eta1:phi1:nsh1:pt2:eta2:phi2:nsh2:qinv:mindist:dist:corr:qfac");
  fNTuple2 = new TNtuple("nt2","Mixed pairs",
  "pt1:eta1:phi1:nsh1:pt2:eta2:phi2:nsh2:qinv:mindist:dist:corr:qfac");
  Int_t c = 0;
  fOutContainer = new TObjArray(2);
  fOutContainer->AddAt(fNTuple1, c++);
  fOutContainer->AddAt(fNTuple2, c++);
}

void AliTwoTrackRes::Exec(Option_t *) {
//______________________________________________________________________________
// Create true and mixed pairs keeping some track parameters

  double bfield = 5.0;
  static int nr=0;
  if (nr == 0) printf("\tStarting event loop...\n");
  printf("\rProcessing event %8d", nr);
  Double_t mpi = 0.13957; // [GeV/c^2]
  Double_t pidTrk1[AliPID::kSPECIES], pidTrk2[AliPID::kSPECIES];
  Int_t tpcIn = 80;   // [cm]
  Int_t tpcOut = 250; // [cm]
  Double_t tdist[170];
  Double_t tdistrot[170];
  Double_t tpcEnt1[3], tpcEnt2[3], pos1[3];
  TVector3 x1, x2, diff;
  TBits clu1, clu2, sha1, sha2;
  TRandom2 rnd;
  Int_t  ntracks = fESDEvent->GetNumberOfTracks();
  for(Int_t itrack = 0; itrack < ntracks; itrack++) {
    AliESDtrack *track1 = fESDEvent->GetTrack(itrack);
    AliExternalTrackParam *trp1 = const_cast<AliExternalTrackParam*>
      (track1->GetTPCInnerParam());
    if (!trp1) continue;
    if (!track1->IsOn(AliESDtrack::kTPCpid)) continue;
    track1->GetTPCpid(pidTrk1);
    Int_t q1 = trp1->Charge();
    if (!((fTrackCuts->AcceptTrack(track1)) && (q1 == 1) &&
	  (pidTrk1[AliPID::kPion]+pidTrk1[AliPID::kMuon] > 0.5))) continue;
    if (!track1->GetInnerXYZ(tpcEnt1)) continue;
    clu1 = track1->GetTPCClusterMap();
    sha1 = track1->GetTPCSharedMap();
    SetTr1(track1->Pt(), track1->Eta(), track1->Phi(), mpi);
    SetTpcEnt1(tpcEnt1[0], tpcEnt1[1], tpcEnt1[2]);
    for(Int_t jtrack = 0; jtrack < itrack; jtrack++) {
      AliESDtrack *track2 = fESDEvent->GetTrack(jtrack);
      AliExternalTrackParam *trp2 = const_cast<AliExternalTrackParam*>
	(track2->GetTPCInnerParam());
      if (!trp2) continue;
      if (!track2->IsOn(AliESDtrack::kTPCpid)) continue;
      track2->GetTPCpid(pidTrk2);
      Int_t q2 = trp2->Charge();
      if (!((fTrackCuts->AcceptTrack(track2)) && (q2 == 1) &&
	    (pidTrk2[AliPID::kPion]+pidTrk2[AliPID::kMuon] > 0.5))) continue;
      if (!track2->GetInnerXYZ(tpcEnt2)) continue;
      clu2 = track2->GetTPCClusterMap();
      sha2 = track2->GetTPCSharedMap();
      SetTr2(track2->Pt(), track2->Eta(), track2->Phi(), mpi);
      SetTpcEnt2(tpcEnt2[0], tpcEnt2[1], tpcEnt2[2]);
      for (Int_t i = tpcIn; i < tpcOut; i++) { // Minimum distance
	trp1->GetDistance(trp2, (double) i, pos1, bfield);
	x1.SetXYZ(pos1[0], pos1[1], pos1[2]);
	tdist[i-tpcIn] = x1.Mag();
	x1.SetXYZ(-pos1[0], -pos1[1], pos1[2]);
	tdistrot[i-tpcIn] = x1.Mag();
      }
      Double_t mindist = 100000;
      Int_t jmin=0;
      for (Int_t j = 0; j < tpcOut-tpcIn; j++) {
	if (tdist[j] < mindist) {jmin=j;  mindist = tdist[j]; }
      }
      //      Double_t mindist = MinDist(track1, track2);
      Double_t dist = Dist();
      Double_t dphi = DPhi();
      Double_t deta = DEta();
      Int_t    nsh1 = GetNSha(clu1, sha1);
      Int_t    nsh2 = GetNSha(clu2, sha2);
      Double_t corr = Corr(clu1, clu2, sha1, sha2);
      Double_t qfac = Qfac(clu1, clu2, sha1, sha2);
      if ((TMath::Abs(track1->Eta())>0.8)&&(TMath::Abs(track2->Eta())>0.8)) continue;
      if ((TMath::Abs(dphi)<0.35)&&(deta<0.35)) {
      FillNTuple1(mindist,dist,corr,qfac,nsh1,nsh2);}    // True
      Double_t tr2rot = RotTr2Phi();                // Rotate trck2
      SetTr2(track2->Pt(), track2->Eta(), tr2rot, mpi);
      tpcEnt2[0] = -tpcEnt2[0];
      tpcEnt2[1] = -tpcEnt2[1];
      Double_t distrot = Dist();
      Double_t dphirot = DPhi();
      Double_t mindistrot = 100000;
      jmin=0;
      for (Int_t j = 0; j < tpcOut-tpcIn; j++) {
	if (tdistrot[j] < mindistrot) {jmin=j;  mindistrot = tdistrot[j]; }
      }
      if ((TMath::Abs(dphirot)<0.35)&&(deta<0.35)) {
	if (rnd.Rndm() < 0.5) NoSwap();
	else Swap();
	FillNTuple2(mindistrot,distrot,corr,qfac,nsh1,nsh2);} // Mixed
    }
  }
  PostData(0, fOutContainer);
  nr++;
}

void AliTwoTrackRes::Terminate(Option_t *) {
//______________________________________________________________________________
// Write output and clean up

  fOutContainer = (TObjArray*)GetOutputData(0);
  TFile *f1  = new TFile( fOutFilename, "RECREATE" );
  fOutContainer->Write();
  f1->Flush();
  f1->Close();
  delete f1;
  delete fChain;
  delete fNTuple1;
  delete fNTuple2;
  printf("\n");
}

//______________________________________________________________________________
// Miscellaneous methods

// Set tracks
void AliTwoTrackRes::SetTr1(double pt1, double eta1, double phi1, double m) {
  fP1.SetPtEtaPhiM(pt1, eta1, phi1, m);}
void AliTwoTrackRes::SetTr2(double pt2, double eta2, double phi2, double m) {
  fP2.SetPtEtaPhiM(pt2, eta2, phi2, m);}

// Set nominal TPC entrance coordinates
void AliTwoTrackRes::SetTpcEnt1(double x1, double y1, double z1) {
  fTpcEnt1.SetX(x1); fTpcEnt1.SetY(y1); fTpcEnt1.SetZ(z1);}
void AliTwoTrackRes::SetTpcEnt2(double x2, double y2, double z2) {
  fTpcEnt2.SetX(x2); fTpcEnt2.SetY(y2); fTpcEnt2.SetZ(z2);}

double AliTwoTrackRes::MinDist(AliExternalTrackParam *trk1,
			       AliExternalTrackParam *trk2) {
// Calculate minimum track separation within the TPC

  int tpcIn = 0;   // [cm]
  int tpcOut = 170; // [cm]
  double tdist[170], pos[3];
  TVector3 x;
  for (int i = tpcIn; i < tpcOut; i++) {
    trk1->GetDistance(trk2, i, pos, 5000);
    x.SetXYZ(pos[0], pos[1], pos[2]);
    tdist[i-tpcIn] = x.Mag();
  }
  double maxdist = 0.0;
  for (int j = 0; j < tpcOut-tpcIn; j++) {
    if (tdist[j] > maxdist) { maxdist = tdist[j]; }
  }
  double mindist = maxdist;
  for (int j = 0; j < tpcOut-tpcIn; j++) {
    if (tdist[j] < mindist) { mindist = tdist[j]; }
  }
  return mindist;}

int AliTwoTrackRes::GetNSha(TBits cl, TBits sh) {
// Get number of shared clusters

  int ncl = cl.GetNbits();
  int sum = 0;
  for(int i = 0; i < ncl; i++) {
    if (!cl.TestBitNumber(i)) continue;
    int n = sh.TestBitNumber(i);
    sum += n;}
  return sum;}

double AliTwoTrackRes::Corr(TBits cl1,  TBits cl2, TBits sh1, TBits sh2) {
// Calculate correlation coefficient

  int ncl1 = cl1.GetNbits();
  int ncl2 = cl2.GetNbits();
  double sumN = 0;  double sumX = 0;  double sumY = 0;
  double sumXX = 0; double sumYY = 0; double sumXY = 0; double corr = -2.0;
  for(int i = 0; i < ncl1 && i < ncl2; i++) {
    if (!(cl1.TestBitNumber(i)&&cl2.TestBitNumber(i))) continue;
    int x = sh1.TestBitNumber(i);
    int y = sh2.TestBitNumber(i);
    sumN += 1.0;
    sumX += x;
    sumY += y;
    sumXX += x*x;
    sumYY += y*y;
    sumXY += x*y;
  }
  double meanX = sumX/sumN;
  double meanY = sumY/sumN;
  double meanXX = sumXX/sumN;
  double meanYY = sumYY/sumN;
  double meanXY = sumXY/sumN;
  double sX = TMath::Sqrt(TMath::Abs(meanXX-meanX*meanX));
  double sY = TMath::Sqrt(TMath::Abs(meanYY-meanY*meanY));
  if (sX*sY!=0) corr = (meanXY-meanX*meanY)/(sX*sY);
  return corr;}

double AliTwoTrackRes::Qfac(TBits cl1,  TBits cl2, TBits sh1, TBits sh2) {
// Quality factor from AliFemto

  int ncl1 = cl1.GetNbits();
  int ncl2 = cl2.GetNbits();
  int sumCls = 0; int sumSha = 0; int sumQ = 0;
  double shfrac = 0; double qfactor = 0;
  for(int i = 0; i < ncl1 && i < ncl2; i++) {
    if (cl1.TestBitNumber(i) && cl2.TestBitNumber(i)) { // Both clusters
      if (sh1.TestBitNumber(i) && sh2.TestBitNumber(i)) { // Shared
	sumQ++;
	sumCls+=2;
	sumSha+=2;}
      else {sumQ--; sumCls+=2;}
    }
    else if (cl1.TestBitNumber(i) || cl2.TestBitNumber(i)) { // Non shared
      sumQ++;
      sumCls++;}
  }
  if (sumCls>0) {
    qfactor = sumQ*1.0/sumCls;
    shfrac = sumSha*1.0/sumCls;
  }
  return qfactor;
}

// Rotate second track for mixed pairs
double AliTwoTrackRes::RotTr2Phi() {
  double rot = TVector2::Phi_mpi_pi(fP2.Phi()+TMath::Pi());
  fTpcEnt2.SetPhi(TVector2::Phi_mpi_pi(fTpcEnt2.Phi()+TMath::Pi()));
  return rot;}

// Fill NTuples
void AliTwoTrackRes::FillNTuple1(double minsep, double sep, double corr,
				 double qf, int ns1, int ns2) {
  fNTuple1->Fill(fP1.Pt(),fP1.Eta(),fP1.Phi(),ns1,fP2.Pt(),fP2.Eta(),
		 fP2.Phi(),ns2,Qinv(),minsep,sep,corr,qf);}
void AliTwoTrackRes::FillNTuple2(double minsep, double sep, double corr,
				 double qf, int ns1, int ns2) {
  fNTuple2->Fill(fPb1.Pt(),fPb1.Eta(),fPb1.Phi(),ns1,fPb2.Pt(),fPb2.Eta(),
		 fPb2.Phi(),ns2,Qinv(),minsep,sep,corr,qf);}

//______________________________________________________________________________
// EOF
