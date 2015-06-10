#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include <TRandom.h>
#include <TObjArray.h>
#include "DetectorK.h"
#include "AliExternalTrackParam.h"
#endif

const double PTMIN = 0.03;
const double PTMAX = 30.0;
const Bool_t VERBOSE = kFALSE;
const double kMaxChi2 = 10;
const int    kNbChi2 = 100;
const int    kNDF = 5;

const double kFactdNdEta = 2;

const UInt_t kFailedBit = 0x1<<14;
// You have to load the class before ... ;-)
// .L DetectorK.cxx++


const double restR[] = {85.20, 100.20, 120.45,   136.10, 160.10, 220.20, 290};

TrackSol* ts=0;
DetectorK* its=0;

TH2F* chiComp=0;

TH1F* solPtHisto=0;
TH1F* effITSSAHisto = 0;
TH1F* hChiSig = 0;
TObjArray  arrMatchHisto;
TObjArray  arrFakeHisto;
TObjArray  chiBgMatch;
TObjArray* sigPoolITS=0,*sigPoolGlo=0;
TObjArray* bgPoolITS=0;
Double_t wghEtaRange=1;

TObjArray* GenerateSigPool(double ptMin=0.1, double ptMax=20, int nPt=100);
TObjArray* GenerateBgPool(int nFactor, double ptMin, double ptMax, double refR, double maxDZ, TObjArray* sigPool);
Double_t   GetRandomPt(double ptMin,double ptMax);
TH1F*      TestMatching(TrackSol* trGlo, double rLastUpd, double refR, int nPhiTest);
TObjArray* TestMatchingN(TrackSol* trGlo, const double* rLastUpd, int nRLastUpd, double refR, int nPhiTest);
TH1F*      PrepChi2Ndf(int ndf=kNDF, int ngen=1000000);
Double_t   CalcCorrProb(TH1* chiSig, TH1* chiBg, double &smBgRet, double maxChi=-1);
Bool_t     PassFastCheck(AliExternalTrackParam &tr1, AliExternalTrackParam &tr2, double cut=10);

void testDetRS(double ptMin=0.1, double ptMax=20, int nPt=100, int nFactor=200, double maxDZ=5, double refR=45) 
{
  //
  if (ptMin<PTMIN || ptMin>PTMAX || ptMax<PTMIN || ptMax>PTMAX || ptMin>=ptMax) {
    printf("Pt range should be within %f %f\n",PTMIN,PTMAX);
    return;
  }
  its = new DetectorK("ALICE","ITS");
  its->SetAvgRapidity(0.2);
  //
  its->AddLayer((char*)"vertex",   0,      0); // dummy vertex for matrix calculation
  its->AddLayer((char*)"bpipe",  2.0, 0.0022);
  // new ideal Pixel properties?
  Double_t x0IB     = 0.003;
  Double_t x0OB     = 0.008;
  Double_t resRPhiIB     = 0.0004;
  Double_t resZIB        = 0.0004;
  Double_t resRPhiOB     = 0.0004;
  Double_t resZOB        = 0.0004;
  Double_t eff           = 0.95;
  //
  its->AddLayer((char*)"ddd1",  2.32 ,  x0IB, resRPhiIB, resZIB,eff); 
  its->AddLayer((char*)"ddd2",  3.13 ,  x0IB, resRPhiIB, resZIB,eff); 
  its->AddLayer((char*)"ddd3",  3.91 ,  x0IB, resRPhiIB, resZIB,eff); 
  its->AddLayer((char*)"ddd4",  19.41,  x0OB, resRPhiOB, resZOB,eff); 
  its->AddLayer((char*)"ddd5",  24.71 ,  x0OB, resRPhiOB, resZOB,eff); 
  its->AddLayer((char*)"ddd6",  35.33 ,  x0OB, resRPhiOB, resZOB,eff); 
  its->AddLayer((char*)"ddd7",  40.53 ,  x0OB, resRPhiOB, resZOB,eff); 
  its->AddLayer((char*)"shell",   43.00 ,  0.01);
  //
  printf("Reference radius: %f\n",refR);
  //
  its->SetAtLeastHits(5);
  its->SetAtLeastCorr(5);
  its->SetAtLeastFake(1);
  //
  its->PrintLayout();
  //
  sigPoolITS = GenerateSigPool(ptMin,ptMax,nPt);
  bgPoolITS  = GenerateBgPool(nFactor, ptMin,ptMax, refR, maxDZ, sigPoolITS);
  //
  its->AddTPC(0.1,0.1);

  its->AddTRD();

  //
  sigPoolGlo = GenerateSigPool(ptMin,ptMax,nPt);
  //
  //  its->SolveTrack(*ts);
  //  its->CalcITSEff(*ts,kTRUE);
  //
  effITSSAHisto = (TH1F*)solPtHisto->Clone("itsSAeff");
  effITSSAHisto->SetTitle("ITS SA eff");
  effITSSAHisto->Reset();
  hChiSig = PrepChi2Ndf();
  int nrTest = sizeof(restR)/sizeof(double);
  for (int ir=0;ir<nrTest;ir++) {
    TH1F* hmt = (TH1F*)solPtHisto->Clone(Form("its_tpc_match_rupd%.1f_rref%.1f",restR[ir],refR));
    hmt->SetTitle(Form("ITS-TPC match eff, R TPCmin = %.1f, Rmatch = f%.1f",restR[ir],refR));
    hmt->Reset();
    arrMatchHisto.AddLast(hmt);
    TH1F* hfm = (TH1F*)solPtHisto->Clone(Form("its_tpc_fake_rupd%.1f_rref%.1f",restR[ir],refR));
    hfm->SetTitle(Form("ITS-TPC fake match eff, R TPCmin = %.1f, Rmatch = f%.1f",restR[ir],refR));
    hfm->Reset();
    arrFakeHisto.AddLast(hfm);
  }
  //
  for (int ipt=0;ipt<nPt;ipt++) {
    ts = (TrackSol*)sigPoolITS->At(ipt);
    if (ts && !ts->TestBit(kFailedBit)) effITSSAHisto->SetBinContent(ipt+1, ts->fProb[2][0]); // combined SA eff
    //
    ts = (TrackSol*)sigPoolGlo->At(ipt);
    if (ts->TestBit(kFailedBit)) continue;
    TObjArray* harr = TestMatchingN(ts,restR,nrTest,refR,50);
    printf("ipt = %d, pt=%.3f : %d histos vs %d tested\n",ipt,ts->fPt, harr ? harr->GetEntries() : 0, nrTest);
    if (!harr) continue;
    chiBgMatch.AddAtAndExpand(harr,ipt);
    //
    for (int ir=0;ir<nrTest;ir++) {
      TH1* hbg = (TH1*)harr->At(ir);
      if (!hbg) continue;
      double bg = 0;
      double vl = CalcCorrProb(hChiSig,hbg,bg); // matching prob
      double effITSSA = effITSSAHisto->GetBinContent(ipt+1);
      ((TH1*)arrMatchHisto.At(ir))->SetBinContent(ipt+1, vl*effITSSA);
      ((TH1*)arrFakeHisto.At(ir))->SetBinContent(ipt+1, (1.-vl)*effITSSA + bg*(1.-effITSSA));
    }
  }
  //
}

TObjArray* GenerateSigPool(double ptMin, double ptMax, int nPt)
{
  // evaluate signal tracks
  TObjArray* arr = new TObjArray(nPt);
  //
  if (!solPtHisto) {
    double dx = log(ptMax/ptMin)/nPt;
    double *xax = new Double_t[nPt+1];
    for (int i=0;i<=nPt;i++) xax[i]= ptMin*exp(dx*i);
    solPtHisto = new TH1F("solPtHisto","solved pt",nPt, xax);
    delete[] xax;
  }
  //
  for (int ip=0;ip<nPt;ip++) {
    double pt = solPtHisto->GetBinCenter(ip+1);
    ts = new TrackSol(its->GetNumberOfLayers(), pt, its->GetAvgRapidity(), -1);
    if (!its->SolveTrack(*ts) || !its->CalcITSEff(*ts,VERBOSE)) ts->SetBit(kFailedBit);
    arr->AddAtAndExpand(ts,ip);
  }
  //
  return arr;
}

TObjArray* GenerateBgPool(int nFactor, double ptMin, double ptMax, double refR, double maxDZ, TObjArray* sigPool)
{
  // Generate bg ITS tracks, which may fall within maxDZ of the average 
  // rapidity track extrapolation at refR
  //
  enum {kY,kZ,kSnp,kTgl,kPtI};              // track parameter aliases
  enum {kY2,kYZ,kZ2,kYSnp,kZSnp,kSnp2,kYTgl,kZTgl,kSnpTgl,kTgl2,kYPtI,kZPtI,kSnpPtI,kTglPtI,kPtI2,kNCovPar}; // cov.matrix aliases
  //
  // estimate eta range of generation
  double thtAv = 2*TMath::ATan(TMath::Exp(-its->GetAvgRapidity()));
  double zAv = refR*TMath::Tan(thtAv);
  double zMin = zAv - maxDZ;
  double zMax = zAv + maxDZ;
  double thtMin = TMath::ATan(zMin/refR);
  double thtMax = TMath::ATan(zMax/refR);
  double etaMax =-TMath::Log( TMath::Tan(thtMin/2.));
  double etaMin =-TMath::Log( TMath::Tan(thtMax/2.));
  //
  wghEtaRange = (etaMax-etaMin); // weight wrt to dn/deta due to the restricted eta range
  AliExternalTrackParam probTr;
  double bGauss = its->GetBField()*10;
  //
  int ndNdEta = its->GetdNdEtaCent()*kFactdNdEta;
  int nGen = ndNdEta*nFactor;
  //
  int lrRef = its->FindLayerID(refR,0); // compare at this layer
  //
  TObjArray* genArr = new TObjArray();
  //
  for (int it=0;it<nGen;it++) {
    double pt = GetRandomPt(ptMin,ptMax);
    //
    int ptBin = solPtHisto->FindBin(pt)-1;
    TrackSol* trSol = (TrackSol*)sigPool->At(ptBin); // fetch the solved track for this pt
    //
    if (trSol->TestBit(kFailedBit)) continue;
    //
    double eta = etaMin+(etaMax-etaMin)*gRandom->Rndm();
    int q = gRandom->Rndm()>0.5 ? -1:1;
    double lambda = TMath::Pi()/2.0 - 2.0*TMath::ATan(TMath::Exp(-eta));
    probTr.Reset();
    double *trPars = (double*)probTr.GetParameter();
    double *trCov  = (double*)probTr.GetCovariance();
    trPars[kY] = 0;                         // start from Y = 0
    trPars[kZ] = 0;                         //            Z = 0 
    trPars[kSnp] = 0;                       //            track along X axis at the vertex
    trPars[kTgl] = TMath::Tan(lambda);      //            dip
    trPars[kPtI] = q/pt;                    //            q/pt      
    trCov[kY2] = trCov[kZ2] = trCov[kSnp2] = trCov[kTgl2] = trCov[kPtI2] = 1e-9;
    if (!its->PropagateToR(&probTr,refR,bGauss,1)) {
      printf("FAILED to propagate to r=%f\n",refR);
      probTr.Print();
      continue;
    }
    double pos[3];
    probTr.GetXYZ(pos);  // lab position
    double phi = TMath::ATan2(pos[1],pos[0]);
    if ( TMath::Abs(TMath::Abs(phi)-TMath::Pi()/2)<1e-3) phi = 0;//TMath::Sign(TMath::Pi()/2 - 1e-3,phi);
    if (!probTr.Rotate(phi))  {
      printf("FAILED to rotate to phi=%f\n",phi);
      probTr.Print();
      continue;
    }
    //
    // assign the error matrix of the ITS outward track with closest pt, estimated with FT
    AliExternalTrackParam* tr = (AliExternalTrackParam*)trSol->fTrackOutA[lrRef]; // and track pars at requested layer
    //printf("Pick at lr %d   %p\n",lrRef,tr);
    //if (tr) tr->Print();
    memcpy(trCov, tr->GetCovariance(), kNCovPar*sizeof(double));
    //
    genArr->AddLast(new AliExternalTrackParam(probTr));
    solPtHisto->Fill(pt);
  }
  //
  return genArr;
}

Double_t GetRandomPt(double ptMin,double ptMax)
{
  // generate random thermal pt
  static TF1* dndpt = 0;
  if (!dndpt) {
    dndpt = new TF1("dndpt","pow(sqrt(x*x+0.14*0.14),1.5)*exp(-sqrt(x*x+0.14*0.14)/0.17)",PTMIN,PTMAX);
    dndpt->SetNpx(1000);
  }
  return dndpt->GetRandom(ptMin,ptMax);
  //
}

TH1F* TestMatching(TrackSol* trGlo, double rLastUpd, double refR, int nPhiTest)
{
  //
  int maxLrAcc = trGlo->fTrackInw.GetEntries()-1-3; // need at least 3 points fitted
  Int_t lrLastUpd = its->FindLayerID(rLastUpd, -1); // find the layer below last update R
  if (lrLastUpd<0 || trGlo->TestBit(kFailedBit)) return 0;
  //
  if (lrLastUpd >= maxLrAcc) return 0;
  // pick the extrapolation to the layer below last requested update layer
  AliExternalTrackParam* trSignal = (AliExternalTrackParam*)trGlo->fTrackInw[lrLastUpd];
  if (!trSignal) return 0; 
  AliExternalTrackParam trSignalC = *trSignal; // use its copy
  //
  if (!its->ExtrapolateToR(&trSignalC, refR, trGlo->fMass)) return 0;
  // bring to the tracking frame of the matching layer
  //
  double pos[3];
  trSignalC.GetXYZ(pos);  // lab position
  double phi = TMath::ATan2(pos[1],pos[0]);
  if ( TMath::Abs(TMath::Abs(phi)-TMath::Pi()/2)<1e-3) phi = 0;//TMath::Sign(TMath::Pi()/2 - 1e-3,phi);
  if (!trSignalC.Rotate(phi))  {
    printf("FAILED to rotate to phi=%f\n",phi);
    trSignalC.Print();
    return 0;
  }
  //
  int nbg = bgPoolITS->GetEntriesFast();
  double posBg[3];
  double xBg,alpBg,parBgOr[5],covBgOr[15],*parBg,*covBg;
  AliExternalTrackParam trBg;
  //
  TString tts = Form("hchi_pt%.1f_r%.1f",trSignalC.Pt()*1000,rLastUpd);
  TH1F* hchi = new TH1F(tts.Data(),tts.Data(),kNbChi2,0,kMaxChi2);
  hchi->Sumw2();
  double bGauss = its->GetBField()*10;
  //
  for (int ib=0;ib<nbg;ib++) {
    trBg = *((AliExternalTrackParam*)bgPoolITS->At(ib));
    trBg.GetXYZ(posBg);  // lab position
    //
    xBg     = trBg.GetX();
    alpBg   = trBg.GetAlpha();
    parBg = (double*)trBg.GetParameter();
    covBg = (double*)trBg.GetCovariance();
    //
    memcpy(&parBgOr,parBg,5*sizeof(double));
    memcpy(&covBgOr,covBg,15*sizeof(double));
    //
    // simulate random phi
    for (int iphi=0;iphi<nPhiTest;iphi++) {
      double alpBgR = (gRandom->Rndm() - 0.5)*TMath::Pi();
      trBg.Set(xBg,alpBgR,parBg,covBg);
      //
      if (trBg.Rotate(trSignalC.GetAlpha()) && trBg.PropagateTo( trSignalC.GetX(), bGauss)) {
	if (!PassFastCheck(trSignalC,trBg)) continue;
	double chi2 = trSignalC.GetPredictedChi2(&trBg);
	//	trSignalC.Print();
	//	trBg.Print();
	//	printf("chi2: %f\n\n\n",chi2);
	hchi->Fill(chi2/kNDF);
      }
      //
      trBg.Set(xBg,alpBg,parBgOr, covBgOr);  // restore original track
    }
    //
  }
  //
  double nTest = nbg*nPhiTest;
  if (nTest<1) return 0;
  double nCorr = wghEtaRange*its->GetdNdEtaCent()*kFactdNdEta;
  printf("Expected: %.0f Tested: %.0f\n",nCorr, nTest);
  hchi->Scale(nCorr/nTest);
  return hchi;
}

//______________________________________________________________________________
TObjArray* TestMatchingN(TrackSol* trGlo, const double* rLastUpd, int nRLastUpd, double refR, int nPhiTest)
{
  //
  if (trGlo->TestBit(kFailedBit)) return 0;
  AliExternalTrackParam trSignalC[nRLastUpd];
  //
  TObjArray *arrH = new TObjArray(nRLastUpd);
  int maxLrAcc = trGlo->fTrackInw.GetEntries()-1-3; // need at least 3 points fitted
  //
  int sgRef = -1;
  for (int ir=0;ir<nRLastUpd;ir++) {
    trSignalC[ir].SetBit(kFailedBit);
    double rl = rLastUpd[ir];
    Int_t lrLastUpd = its->FindLayerID(rl, -1); // find the layer below last update R
    if (lrLastUpd<0) continue;
    if (lrLastUpd >= maxLrAcc) continue;
    //
    // pick the extrapolation to the layer below last requested update layer
    AliExternalTrackParam* trSignal = (AliExternalTrackParam*)trGlo->fTrackInw[lrLastUpd];
    if (!trSignal) continue; 
    trSignalC[ir] = *trSignal; // use its copy
    //
    //    printf("\n\n\nPicked at lr=%d r=%f  ",lrLastUpd, rl); trSignal->Print();
    if (!its->ExtrapolateToR(&trSignalC[ir], refR, trGlo->fMass)) return 0;
    // bring to the tracking frame of the matching layer
    //
    double pos[3];
    trSignalC[ir].GetXYZ(pos);  // lab position
    double phi = TMath::ATan2(pos[1],pos[0]);
    if (!trSignalC[ir].Rotate(phi))  {
      printf("FAILED to rotate to phi=%f\n",phi);
      trSignalC[ir].Print();
      continue;
    }
    trSignalC[ir].ResetBit(kFailedBit); // validate the track
    //
    //    printf("Taken to r=%f  ",refR); trSignalC[ir].Print();
    //
    TString tts = Form("hchi_pt%.1f_r%.1f",trGlo->fPt*1000,rl);
    TH1F* hchi = new TH1F(tts.Data(),tts.Data(),kNbChi2,0,kMaxChi2);
    hchi->Sumw2();
    arrH->AddAt(hchi,ir);
    if (sgRef<0) sgRef = ir;
  }
  //
  if (sgRef<0) return 0;
  //
  int nbg = bgPoolITS->GetEntriesFast();
  double posBg[3];
  double xBg,alpBg,parBgOr[5],covBgOr[15],*parBg,*covBg;
  AliExternalTrackParam trBg;
  //
  double bGauss = its->GetBField()*10;
  //
  for (int ib=0;ib<nbg;ib++) {
    trBg = *((AliExternalTrackParam*)bgPoolITS->At(ib));
    trBg.GetXYZ(posBg);  // lab position
    //
    xBg     = trBg.GetX();
    alpBg   = trBg.GetAlpha();
    parBg = (double*)trBg.GetParameter();
    covBg = (double*)trBg.GetCovariance();
    //
    memcpy(&parBgOr,parBg,5*sizeof(double));
    memcpy(&covBgOr,covBg,15*sizeof(double));
    //
    // simulate random phi
    for (int iphi=0;iphi<nPhiTest;iphi++) {
      double alpBgR = (gRandom->Rndm() - 0.5)*TMath::Pi();
      trBg.Set(xBg,alpBgR,parBg,covBg);
      //
      if (trBg.Rotate(trSignalC[sgRef].GetAlpha()) && trBg.PropagateTo( trSignalC[sgRef].GetX(), bGauss)) {
	for (int ir=0;ir<nRLastUpd;ir++) {
	  if (trSignalC[ir].TestBit(kFailedBit)) continue;
	  if (!PassFastCheck(trSignalC[ir],trBg)) continue;
	  double chi2 = trSignalC[ir].GetPredictedChi2(&trBg);
	  ((TH1*)arrH->UncheckedAt(ir))->Fill(chi2/kNDF);
	}
	//
      }
      trBg.Set(xBg,alpBg,parBgOr, covBgOr);  // restore original track
      //
    }
  }
  //
  double nTest = nbg*nPhiTest;
  if (nTest<1) return 0;
  double nCorr = wghEtaRange*its->GetdNdEtaCent()*kFactdNdEta;
  printf("Expected: %.0f Tested: %.0f\n",nCorr, nTest);
  for (int ir=0;ir<nRLastUpd;ir++) {
    TH1* hchi = (TH1*)arrH->UncheckedAt(ir);
    if (!hchi) continue;
    hchi->Scale(nCorr/nTest);
  }
  return arrH;
}

//_________________________________________________________________________________________
TH1F* PrepChi2Ndf(int ndf, int ngen)
{
  // generate distribution for chi2 with given ndf
  TString ttl = Form("Chi2_NDF%d",ndf);
  TH1F* hh = new TH1F(ttl.Data(),ttl.Data(),kNbChi2,0,kMaxChi2);
  for (int i=0;i<ngen;i++) hh->Fill(TMath::ChisquareQuantile(gRandom->Rndm(),ndf)/ndf);
  hh->Scale(1./ngen);
  return hh;
}

//____________________________________________________
Double_t CalcCorrProb(TH1* chiSig, TH1* chiBg, double &smBgRet, double maxChi) 
{
  // Calculate probability of all fake matches being worse than correct one
  // The chi2 histos should be normalized
  int nb = chiSig->GetNbinsX();
  if (maxChi>0) {
    int nt = chiSig->FindBin(maxChi);
    if (nt<nb) nb = nt;
  }
  //
  double bgchiC = maxChi;
  if (bgchiC<0) bgchiC = 5.;
  int maxBgBinCut = chiSig->FindBin(bgchiC);
  double smBg = 0;
  smBgRet = 0;
  double probTotCorr = 0;
  for (int i=1;i<=nb;i++) {
    double wSig = chiSig->GetBinContent(i);
    probTotCorr += wSig*TMath::Exp(-smBg);  // Poisson prob of not having a bg with better chi2
    smBg += chiBg->GetBinContent(i);
    if (i<=maxBgBinCut) smBgRet = smBg;
  }
  //
  smBg = 1.-TMath::Exp(smBg); // prob. to find bg within a cut in absence of correct match
  return probTotCorr;
}

//__________________________________________________
Bool_t PassFastCheck(AliExternalTrackParam &tr1, AliExternalTrackParam &tr2, double cut)
{
  enum {kY,kZ,kSnp,kTgl,kPtI};              // track parameter aliases
  enum {kY2,kYZ,kZ2,kYSnp,kZSnp,kSnp2,kYTgl,kZTgl,kSnpTgl,kTgl2,kYPtI,kZPtI,kSnpPtI,kTglPtI,kPtI2}; // cov.matrix aliases
  const int erId[5] = {kY2,kZ2,kSnp2,kTgl2,kPtI2};
  //
  // fast check ignoring non-diagonal elements
  double *cov1 = (double*)tr1.GetCovariance();
  double *cov2 = (double*)tr2.GetCovariance();
  double *par1 = (double*)tr1.GetParameter();
  double *par2 = (double*)tr2.GetParameter();
  //
  double chi2 = 0;
  for (int i=0;i<5;i++) {
    double df = par1[i]-par2[i];
    double err2 = cov1[erId[i]]+cov2[erId[i]];
    chi2 += df*df/err2;
  }
  return (chi2>cut*kNDF) ? kFALSE : kTRUE;
  //
}
