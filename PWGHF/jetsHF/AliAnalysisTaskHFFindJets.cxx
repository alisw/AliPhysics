#include "AliAnalysisTaskHFFindJets.h"

#include <TChain.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TSystem.h>
#include <TH1F.h>
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliVertexerTracks.h"
#include "AliFJWrapper.h"
#include "FJ_includes.h"
#include "AliAnalysisTaskSE.h"


ClassImp(AliAnalysisTaskHFFindJets)
//______________________________________________________________________________
AliAnalysisTaskHFFindJets::AliAnalysisTaskHFFindJets():
  AliAnalysisTaskSE("HFFindJets"),
	fOutput{nullptr},
	hpt_nocuts{nullptr},
	htgl_nocuts{nullptr},
	hpt_cuts{nullptr},
	hdcatoprimxy_cuts{nullptr},
	htgl_cuts{nullptr},
	hvx{nullptr},
	hvy{nullptr},
	hvz{nullptr},
	hvx3{nullptr},
	hvy3{nullptr},
	hvz3{nullptr},
	hitsmap{nullptr},

	hvertexx{nullptr},
	hvertexy{nullptr},
	hvertexz{nullptr},

	hdecayxyz{nullptr},
	hdecayxy{nullptr},
	hmass0{nullptr},
	hmassP{nullptr},
	hptD0{nullptr},
	hptprong0{nullptr},
	hptprong1{nullptr},
	hd0{nullptr},
	hd0d0{nullptr},
	hImpParErr{nullptr},
	hDecLenErr{nullptr},
	hDecLenXYErr{nullptr},
	hCovPVXX{nullptr},
	hCovSVXX{nullptr},

	hjetpt{nullptr},
	hjetzg{nullptr},
	hjetrg{nullptr},
	hjetnsd{nullptr}
{
  InitDefault();

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}
//_______________________________________________________________________________________
void AliAnalysisTaskHFFindJets::InitDefault()
{
   kbitDplus = 0;
   kbitDs = 1;
   kbitLc = 2;
   
   fMassDzero = TDatabasePDG::Instance()->GetParticle(421)->Mass();
   fMassDplus = TDatabasePDG::Instance()->GetParticle(411)->Mass();
   fMassDs = TDatabasePDG::Instance()->GetParticle(431)->Mass();
   fMassLambdaC = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
                                                        	
}
//_______________________________________________________________________________________
char* AliAnalysisTaskHFFindJets::GetJsonString(const char* jsonFileName, const char* key)
{
  FILE* fj = fopen(jsonFileName, "r");
  char line[500];
  char* value = 0x0;
  while (!feof(fj)) {
    fgets(line, 500, fj);
    if (strstr(line, key)) {
      value = strtok(line, ":");
      value = strtok(NULL, ":");
      break;
    }
  }
  fclose(fj);
  return value;
}

//_______________________________________________________________________________________
int AliAnalysisTaskHFFindJets::GetJsonInteger(const char* jsonFileName, const char* key)
{
  FILE* fj = fopen(jsonFileName, "r");
  char line[500];
  int value = -999;
  while (!feof(fj)) {
    fgets(line, 500, fj);
    if (strstr(line, key)) {
      char* token = strtok(line, ":");
      token = strtok(NULL, ":");
      TString temp = token;
      temp.ReplaceAll("\"", "");
      temp.ReplaceAll(",", "");
      value = temp.Atoi();
      break;
    }
  }
  fclose(fj);
  return value;
}

//_______________________________________________________________________________________
bool AliAnalysisTaskHFFindJets::GetJsonBool(const char* jsonFileName, const char* key)
{
  FILE* fj = fopen(jsonFileName, "r");
  char line[500];
  bool value = false;
  while (!feof(fj)) {
    fgets(line, 500, fj);
    if (strstr(line, key)) {
      char* token = strtok(line, ":");
      token = strtok(NULL, ":");
      TString temp = token;
      temp.ReplaceAll("\"", "");
      temp.ReplaceAll(",", "");
      if (temp.Contains("true"))
        value = true;
      break;
    }
  }
  fclose(fj);
  return value;
}

//_______________________________________________________________________________________
float AliAnalysisTaskHFFindJets::GetJsonFloat(const char* jsonFileName, const char* key)
{
  FILE* fj = fopen(jsonFileName, "r");
  char line[500];
  float value = -999.;
  while (!feof(fj)) {
    fgets(line, 500, fj);
    if (strstr(line, key)) {
      char* token = strtok(line, ":");
      token = strtok(NULL, ":");
      TString temp = token;
      temp.ReplaceAll("\"", "");
      temp.ReplaceAll(",", "");
      value = temp.Atof();
      break;
    }
  }
  fclose(fj);
  return value;
}

//_______________________________________________________________________________________
void AliAnalysisTaskHFFindJets::ReadJson()
{
  char* value = GetJsonString("dpl-config_std.json", "aod-file");
  printf("%s\n", value);
  int i3p = GetJsonInteger("dpl-config_std.json", "do3prong");
  printf("%d\n", i3p);
  int tin = GetJsonInteger("dpl-config_std.json", "triggerindex");
  printf("%d\n", tin);
  float minpt = GetJsonFloat("dpl-config_std.json", "ptmintrack");
  printf("%f\n", minpt);
  float dcatoprimxymin = GetJsonFloat("dpl-config_std.json", "dcatoprimxymin");
  printf("%f\n", dcatoprimxymin);
  bool doit = GetJsonBool("dpl-config_std.json", "b_propdca");
  printf("%d\n", doit);
  float chi2 = GetJsonFloat("dpl-config_std.json", "d_minrelchi2change");
  printf("%f\n", chi2);
}

//_______________________________________________________________________________________
void AliAnalysisTaskHFFindJets::UserCreateOutputObjects() 
{  
	// create output histos
	fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("OutputHistos1");
    
    hpt_nocuts = new TH1F("hpt_nocuts", " ; pt tracks (#GeV) ; Entries", 100, 0, 10.);
    fOutput->Add(hpt_nocuts);
	htgl_nocuts = new TH1F("htgl_nocuts", "tgl tracks (#GeV)", 100, -5., 5.);
	fOutput->Add(htgl_nocuts);
	hpt_cuts = new TH1F("hpt_cuts", " ; pt tracks (#GeV) ; Entries", 100, 0, 10.);
	fOutput->Add(hpt_cuts);
	hdcatoprimxy_cuts = new TH1F("hdcatoprimxy_cuts", "dca xy to prim. vtx (cm)", 100, -1.0, 1.0);
	fOutput->Add(hdcatoprimxy_cuts);
	htgl_cuts = new TH1F("htgl_cuts", "tgl tracks (#GeV)", 100, -5., 5.);
	fOutput->Add(htgl_cuts);
	hvx = new TH1F("hvx", " Secondary vertex ; X vertex (cm) ; Entries", 1000, -2.0, 2.0);
	fOutput->Add(hvx);
	hvy = new TH1F("hvy", " Secondary vertex ; Y vertex (cm) ; Entries", 1000, -2.0, 2.0);
	fOutput->Add(hvy);
	hvz = new TH1F("hvz", " Secondary vertex ; Z vertex (cm) ; Entries", 1000, -20.0, 20.0);
	fOutput->Add(hvz);
	hvx3 = new TH1F("hvx3", " Secondary vertex 3prong ; X vertex (cm) ; Entries", 1000, -2.0, 2.0);
	fOutput->Add(hvx3);
	hvy3 = new TH1F("hvy3", " Secondary vertex 3prong ; Y vertex (cm) ; Entries", 1000, -2.0, 2.0);
	fOutput->Add(hvy3);
	hvz3 = new TH1F("hvz3", " Secondary vertex 3prong ; Z vertex (cm) ; Entries", 1000, -20.0, 20.0);
	fOutput->Add(hvz3);
	hitsmap = new TH1F("hitsmap", "hitsmap_cuts", 64, -0.5, 63.5);
	fOutput->Add(hitsmap);

	hvertexx = new TH1F("hvertexx", " Primary vertex ; X vertex (cm) ; Entries", 100, -0.5, 0.5);
	fOutput->Add(hvertexx);
	hvertexy = new TH1F("hvertexy", " Primary vertex ; Y vertex (cm) ; Entries", 100, -0.5, 0.5);
	fOutput->Add(hvertexy);
	hvertexz = new TH1F("hvertexz", " Primary vertex ; Z vertex (cm) ; Entries", 100, -20.0, 20.0);
	fOutput->Add(hvertexz);

	hdecayxyz = new TH1F("hdecayxyz", "hdecayxyz", 200, 0., 2.0);
	fOutput->Add(hdecayxyz);
	hdecayxy = new TH1F("hdecayxy", "hdecayxy", 200, 0., 2.0);
	fOutput->Add(hdecayxy);
	hmass0 = new TH1F("hmass0", "; Inv Mass (GeV/c^{2})", 500, 0, 5.0);
	fOutput->Add(hmass0);
	hmassP = new TH1F("hmassP", "; Inv Mass (GeV/c^{2})", 500, 1.6, 2.1);
	fOutput->Add(hmassP);
	hptD0 = new TH1F("hptD0", " ; pt D0 (#GeV) ; Entries", 100, 0, 10.);
	fOutput->Add(hptD0);
	hptprong0 = new TH1F("hptprong0", " ; pt prong0 (#GeV) ; Entries", 100, 0, 10.);
	fOutput->Add(hptprong0);
	hptprong1 = new TH1F("hptprong1", " ; pt prong1 (#GeV) ; Entries", 100, 0, 10.);
	fOutput->Add(hptprong1);
	hd0 = new TH1F("hd0", "dca xy to prim. vertex (cm)", 100, -1.0, 1.0);
	fOutput->Add(hd0);
	hd0d0 = new TH1F("hd0d0", "product of dca xy to prim. vertex (cm^{2})", 500, -1.0, 1.0);
	fOutput->Add(hd0d0);
	hImpParErr = new TH1F("hImpParErr", "impact parameter error", 100, -1.0, 1.0);
	fOutput->Add(hImpParErr);
	hDecLenErr = new TH1F("hDecLenErr", "decay length error", 100, 0., 1.0);
	fOutput->Add(hDecLenErr);
	hDecLenXYErr = new TH1F("hDecLenXYErr", "decay length XY error", 100, 0., 1.0);
	fOutput->Add(hDecLenXYErr);
	hCovPVXX = new TH1F("hCovPVXX", "XX element of PV cov. matrix", 100, 0., 1.0e-4);
	fOutput->Add(hCovPVXX);
	hCovSVXX = new TH1F("hCovSVXX", "XX element of SV cov. matrix", 100, 0., 0.2);
	fOutput->Add(hCovSVXX);

	hjetpt = new TH1F("hjetpt", " ; pt jet (#GeV) ; Entries", 100, 0., 100.);
	fOutput->Add(hjetpt);
	hjetzg = new TH1F("hjetzg", " ; jet zg ; Entries", 100, 0., 100.);
	fOutput->Add(hjetzg);
	hjetrg = new TH1F("hjetrrg", " ; jet rg ; Entries", 100, 0., 100.);
	fOutput->Add(hjetrg);
	hjetnsd = new TH1F("hjetnsd", " ; jet nsd ; Entries", 100, 0., 100.);
	fOutput->Add(hjetnsd);

	PostData(1,fOutput);
}

//_______________________________________________________________________________________
Int_t AliAnalysisTaskHFFindJets::GetpTBin(Double_t candpT)
{

  Double_t pTBins[npTBins + 1] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 12.0, 16.0, 20.0, 24.0, 36.0, 50.0, 100.0};
  for (Int_t i = 0; i < npTBins; i++) {
    if (candpT >= pTBins[i] && candpT < pTBins[i + 1]) {
      return i;
    }
  }
  return -1;
}

//_______________________________________________________________________________________
Bool_t AliAnalysisTaskHFFindJets::GetTrackMomentumAtSecVert(AliESDtrack* tr, AliAODVertex* secVert, Double_t momentum[3], float fBzkG)
{
  /// fast calculation (no covariance matrix treatment) of track momentum at secondary vertex

  Double_t alpha = tr->GetAlpha();
  Double_t sn = TMath::Sin(alpha), cs = TMath::Cos(alpha);
  Double_t x = tr->GetX(), y = tr->GetParameter()[0], snp = tr->GetParameter()[2];
  Double_t xv = secVert->GetX() * cs + secVert->GetY() * sn;
  Double_t yv = -secVert->GetX() * sn + secVert->GetY() * cs;
  x -= xv;
  y -= yv;
  Double_t crv = tr->GetC(fBzkG);
  if (TMath::Abs(fBzkG) < 0.000001)
    crv = 0.;
  double csp = TMath::Sqrt((1. - snp) * (1. + snp));

  Double_t tgfv = -(crv * x - snp) / (crv * y + csp);
  cs = 1. / TMath::Sqrt(1 + tgfv * tgfv);
  sn = cs < 1. ? tgfv * cs : 0.;

  x = xv * cs + yv * sn;
  Double_t alpNew = alpha + TMath::ASin(sn);
  Double_t ca = TMath::Cos(alpNew - alpha), sa = TMath::Sin(alpNew - alpha);
  Double_t p2 = tr->GetSnp();
  Double_t xNew = tr->GetX() * ca + tr->GetY() * sa;
  Double_t p2New = p2 * ca - TMath::Sqrt((1. - p2) * (1. + p2)) * sa;
  momentum[0] = tr->GetSigned1Pt();
  momentum[1] = p2New * (x - xNew) * tr->GetC(fBzkG);
  momentum[2] = tr->GetTgl();
  Bool_t retCode = tr->Local2GlobalMomentum(momentum, alpNew);
  return retCode;
}

//_______________________________________________________________________________________
Bool_t AliAnalysisTaskHFFindJets::SingleTrkCuts(AliESDtrack* trk, AliESDtrackCuts* esdTrackCuts, AliESDVertex* fV1, Double_t fBzkG)
{
  if (!trk->PropagateToDCA(fV1, fBzkG, kVeryBig))
    return kFALSE;
  trk->RelateToVertex(fV1, fBzkG, kVeryBig);
  return esdTrackCuts->AcceptTrack(trk);
}

//_______________________________________________________________________________________
Bool_t AliAnalysisTaskHFFindJets::SingleTrkCutsSimple(AliESDtrack* trk, Int_t minclutpc, int ptmintrack, double dcatoprimxymin, AliESDVertex* fV1, Double_t fBzkG)
{
  Int_t status = trk->GetStatus();
  bool sel_track = status & AliESDtrack::kITSrefit && (trk->HasPointOnITSLayer(0) || trk->HasPointOnITSLayer(1));
  sel_track = sel_track && trk->GetNcls(1) >= minclutpc;
  sel_track = sel_track && trk->Pt() > ptmintrack;
  AliExternalTrackParam* track = (AliExternalTrackParam*)trk;
  double b[2];
  double bCov[3];
  track->PropagateToDCA(fV1, fBzkG, 100., b, bCov);
  sel_track = sel_track && abs(b[0]) > dcatoprimxymin;
  return sel_track;
}

//_______________________________________________________________________________________
Int_t AliAnalysisTaskHFFindJets::TwoProngSelectionCuts(AliAODRecoDecayHF2Prong* cand, Double_t candpTMin, Double_t candpTMax)
{

  bool isD0 = true;
  bool isD0bar = true;
  Double_t candpT = cand->Pt();
  if (candpT < candpTMin || candpT >= candpTMax) return 0;
  Int_t pTBin = GetpTBin(candpT);
  if (pTBin==-1) return 0;
  if (cand->Prodd0d0() > fCuts[pTBin][7]) return 0;
  if (cand->CosPointingAngle() < fCuts[pTBin][8]) return 0;
  if (cand->CosPointingAngleXY() < fCuts[pTBin][9]) return 0;
  if (cand->NormalizedDecayLengthXY() < fCuts[pTBin][10]) return 0;
  Double_t decayLengthCut = TMath::Min((cand->P() * 0.0066) + 0.01, 0.06);
  if (TMath::Abs(cand->Normalizedd0Prong(0)) < 0.5 || TMath::Abs(cand->Normalizedd0Prong(1)) < 0.5) return 0;
  if (cand->DecayLength() * cand->DecayLength() < decayLengthCut * decayLengthCut) return 0;
  // if (cand->NormalizedDecayLength() * cand->NormalizedDecayLength() < 1.0) return 0;
  if (TMath::Abs(cand->InvMassD0()-fMassDzero) > fCuts[pTBin][0] ) isD0=false;
  if (TMath::Abs(cand->InvMassD0bar()-fMassDzero) > fCuts[pTBin][0] ) isD0bar=false;
  if (!isD0 && !isD0bar) return 0;

  if (cand->Pt2Prong(0) < fCuts[pTBin][4]*fCuts[pTBin][4] || cand->Pt2Prong(1) < fCuts[pTBin][3]*fCuts[pTBin][3] ) isD0=false;
  if (cand->Pt2Prong(0) < fCuts[pTBin][3]*fCuts[pTBin][3] || cand->Pt2Prong(1) < fCuts[pTBin][4]*fCuts[pTBin][4] ) isD0bar=false;
  if (!isD0 && !isD0bar) return 0;

  if (TMath::Abs(cand->Getd0Prong(0)) > fCuts[pTBin][6] || TMath::Abs(cand->Getd0Prong(1)) > fCuts[pTBin][5] ) isD0=false;
  if (TMath::Abs(cand->Getd0Prong(0)) > fCuts[pTBin][5] || TMath::Abs(cand->Getd0Prong(1)) > fCuts[pTBin][6] ) isD0bar=false;
  if (!isD0 && !isD0bar) return 0;

  Double_t cosThetaStarD0,cosThetaStarD0bar;
  cand->CosThetaStarD0(cosThetaStarD0,cosThetaStarD0bar);
  if (TMath::Abs(cosThetaStarD0) > fCuts[pTBin][2] ) isD0=false;
  if (TMath::Abs(cosThetaStarD0bar) > fCuts[pTBin][2] ) isD0bar=false;
  if (!isD0 && !isD0bar) return 0;


  Int_t returnValue=0;
  if(isD0) returnValue+=1;
  if(isD0bar) returnValue+=2;
  return returnValue;
}

//_______________________________________________________________________________________
AliESDVertex* AliAnalysisTaskHFFindJets::ReconstructSecondaryVertex(AliVertexerTracks* vt, TObjArray* trkArray, AliESDVertex* primvtx, double rmax)
{

  vt->SetVtxStart(primvtx);

  AliESDVertex* trkv = (AliESDVertex*)vt->VertexForSelectedESDTracks(trkArray);
  if (trkv->GetNContributors() != trkArray->GetEntriesFast())
    return 0x0;
  Double_t vertRadius2 = trkv->GetX() * trkv->GetX() + trkv->GetY() * trkv->GetY();
  if(vertRadius2>rmax*rmax) return 0x0;
  return trkv;
}

//_______________________________________________________________________________________
AliAODVertex* AliAnalysisTaskHFFindJets::ConvertToAODVertex(AliESDVertex* trkv)
{
  Double_t pos_[3], cov_[6], chi2perNDF_;
  trkv->GetXYZ(pos_);       // position
  trkv->GetCovMatrix(cov_); // covariance matrix
  chi2perNDF_ = trkv->GetChi2toNDF();
  double dispersion_ = trkv->GetDispersion();
  //  printf(" pos_ %f %f %f \n", pos_[0], pos_[1], pos_[2]);
  AliAODVertex* vertexAOD = new AliAODVertex(pos_, cov_, chi2perNDF_, 0x0, -1, AliAODVertex::kUndef, 2);
  return vertexAOD;
}

//_______________________________________________________________________________________
Int_t AliAnalysisTaskHFFindJets::SelectInvMassAndPt3prong(TObjArray* trkArray, AliAODRecoDecay* rd4massCalc3)
{

  Int_t retval = 0;
  Double_t momentum[3];
  Double_t px[3], py[3], pz[3];
  for (Int_t iTrack = 0; iTrack < 3; iTrack++) {
    AliESDtrack* track = (AliESDtrack*)trkArray->UncheckedAt(iTrack);
    track->GetPxPyPz(momentum);
    px[iTrack] = momentum[0];
    py[iTrack] = momentum[1];
    pz[iTrack] = momentum[2];
  }
  UInt_t pdg3[3];
  Int_t nprongs = 3;
  rd4massCalc3->SetPxPyPzProngs(nprongs, px, py, pz);
  Double_t minv2, mrange;
  Double_t lolim, hilim;
  mrange = 0.1;
  lolim = fMassDplus - mrange;
  hilim = fMassDplus + mrange;
  pdg3[0] = 211;
  pdg3[1] = 321;
  pdg3[2] = 211;
  minv2 = rd4massCalc3->InvMass2(nprongs, pdg3);
  if (minv2 > lolim * lolim && minv2 < hilim * hilim)
    retval += (1 << kbitDplus);
  lolim = fMassDs - mrange;
  hilim = fMassDs + mrange;
  for (Int_t ih = 0; ih < 2; ih++) {
    Int_t k = ih * 2;
    pdg3[k] = 321;
    pdg3[1] = 321;
    pdg3[2 - k] = 211;
    minv2 = rd4massCalc3->InvMass2(nprongs, pdg3);
    if (minv2 > lolim * lolim && minv2 < hilim * hilim && (retval & (1 << kbitDs)) == 0)
      retval += (1 << kbitDs);
  }
  lolim = fMassLambdaC - mrange;
  hilim = fMassLambdaC + mrange;
  pdg3[0] = 2212;
  pdg3[1] = 321;
  pdg3[2] = 211;
  minv2 = rd4massCalc3->InvMass2(nprongs, pdg3);
  if (minv2 > lolim * lolim && minv2 < hilim * hilim && (retval & (1 << kbitLc)) == 0)
    retval += (1 << kbitLc);
  pdg3[0] = 211;
  pdg3[1] = 321;
  pdg3[2] = 2212;
  minv2 = rd4massCalc3->InvMass2(nprongs, pdg3);
  if (minv2 > lolim * lolim && minv2 < hilim * hilim && (retval & (1 << kbitLc)) == 0)
    retval += (1 << kbitLc);

  return retval;
}

//_______________________________________________________________________________________
AliAODRecoDecayHF2Prong* AliAnalysisTaskHFFindJets::Make2Prong(TObjArray* twoTrackArray, AliAODVertex* secVert, Double_t fBzkG)
{

  AliESDtrack* track_0 = (AliESDtrack*)twoTrackArray->UncheckedAt(0);
  AliESDtrack* track_1 = (AliESDtrack*)twoTrackArray->UncheckedAt(1);

  Double_t px[2], py[2], pz[2], d0[2], d0err[2];
  Double_t momentum[3];
  GetTrackMomentumAtSecVert(track_0, secVert, momentum, fBzkG);
  px[0] = momentum[0];
  py[0] = momentum[1];
  pz[0] = momentum[2];
  GetTrackMomentumAtSecVert(track_1, secVert, momentum, fBzkG);
  px[1] = momentum[0];
  py[1] = momentum[1];
  pz[1] = momentum[2];

  Float_t d0z0f[2], covd0z0f[3];
  track_0->GetImpactParameters(d0z0f, covd0z0f);
  d0[0] = d0z0f[0];
  d0err[0] = TMath::Sqrt(covd0z0f[0]);
  track_1->GetImpactParameters(d0z0f, covd0z0f);
  d0[1] = d0z0f[0];
  d0err[1] = TMath::Sqrt(covd0z0f[0]);

  Double_t xdummy, ydummy;
  float dcap1n1 = track_0->GetDCA(track_1, fBzkG, xdummy, ydummy);

  AliAODRecoDecayHF2Prong* the2Prong = new AliAODRecoDecayHF2Prong(0x0, px, py, pz, d0, d0err, dcap1n1);
  AliAODVertex* ownsecv=secVert->CloneWithoutRefs();
  the2Prong->SetOwnSecondaryVtx(ownsecv);
  return the2Prong;
}

//_______________________________________________________________________________________
AliAODRecoDecayHF3Prong* AliAnalysisTaskHFFindJets::Make3Prong(TObjArray* threeTrackArray, AliAODVertex* secVert, Double_t fBzkG)
{

  AliESDtrack* track_0 = (AliESDtrack*)threeTrackArray->UncheckedAt(0);
  AliESDtrack* track_1 = (AliESDtrack*)threeTrackArray->UncheckedAt(1);
  AliESDtrack* track_2 = (AliESDtrack*)threeTrackArray->UncheckedAt(2);

  Double_t px[3], py[3], pz[3], d0[3], d0err[3];
  Double_t momentum[3];
  GetTrackMomentumAtSecVert(track_0, secVert, momentum, fBzkG);
  px[0] = momentum[0];
  py[0] = momentum[1];
  pz[0] = momentum[2];
  GetTrackMomentumAtSecVert(track_1, secVert, momentum, fBzkG);
  px[1] = momentum[0];
  py[1] = momentum[1];
  pz[1] = momentum[2];
  GetTrackMomentumAtSecVert(track_2, secVert, momentum, fBzkG);
  px[2] = momentum[0];
  py[2] = momentum[1];
  pz[2] = momentum[2];
  Float_t d0z0f[2], covd0z0f[3];
  track_0->GetImpactParameters(d0z0f, covd0z0f);
  d0[0] = d0z0f[0];
  d0err[0] = TMath::Sqrt(covd0z0f[0]);
  track_1->GetImpactParameters(d0z0f, covd0z0f);
  d0[1] = d0z0f[0];
  d0err[1] = TMath::Sqrt(covd0z0f[0]);
  track_2->GetImpactParameters(d0z0f, covd0z0f);
  d0[2] = d0z0f[0];
  d0err[2] = TMath::Sqrt(covd0z0f[0]);

  Double_t xdummy, ydummy;
  float dcap1n1 = track_0->GetDCA(track_1, fBzkG, xdummy, ydummy);
  float dcap2n1 = track_2->GetDCA(track_1, fBzkG, xdummy, ydummy);
  float dcap1p2 = track_0->GetDCA(track_2, fBzkG, xdummy, ydummy);
  Double_t dca[3] = {dcap1n1, dcap2n1, dcap1p2};
  Double_t dispersion = 0;
  Double_t dist12 = 0.;
  Double_t dist23 = 0.;
  Short_t charge = (Short_t)(track_0->Charge() + track_1->Charge() + track_2->Charge());

  // construct the candidate passing a NULL pointer for the secondary vertex to avoid creation of TRef
  AliAODRecoDecayHF3Prong* the3Prong = new AliAODRecoDecayHF3Prong(0x0, px, py, pz, d0, d0err, dca, dispersion, dist12, dist23, charge);
  // add a pointer to the secondary vertex via SetOwnSecondaryVtx (no TRef created)
  // AliAODVertex* ownsecv=secVert->CloneWithoutRefs();
  // the3Prong->SetOwnSecondaryVtx(ownsecv);
  return the3Prong;
}

//_______________________________________________________________________________________
void AliAnalysisTaskHFFindJets::InitFromJson(TString esdFileName)
{
  // read configuration from json file
  if (esdFileName != "" && gSystem->Exec(Form("ls %s > /dev/null", esdFileName.Data())) == 0) {
	printf("------Read configuration from JSON file------\n");
		
	minncluTPC = 50;
	dcatoprimxymin = 0.;
	
    TFile* esdfile = TFile::Open(esdFileName.Data()); 
	tree = (TTree*)esdfile->Get("esdTree");
		
	printf("Read configuration from JSON file\n");
	selectD0 = GetJsonInteger(esdFileName.Data(), "d_selectionFlagD0");
	selectD0bar = GetJsonInteger(esdFileName.Data(), "d_selectionFlagD0bar");
	printf("D0 cuts: %d, D0bar cuts: %d\n", selectD0, selectD0bar);
	ptmintrack = GetJsonFloat(esdFileName.Data(), "ptmintrack");
	printf("Min pt track = %f\n", ptmintrack);
	do3Prongs = GetJsonInteger(esdFileName.Data(), "do3prong");
	printf("do3prong     = %d\n", do3Prongs);
	minncluTPC = GetJsonInteger(esdFileName.Data(), "d_tpcnclsfound");
	printf("minncluTPC   = %d\n", minncluTPC);
	dcatoprimxymin = GetJsonFloat(esdFileName.Data(), "dcatoprimxymin");
	printf("dcatoprimxymin   = %f\n", dcatoprimxymin);
	printf("Read configuration from JSON file\n");
	candpTMin = GetJsonFloat(esdFileName.Data(), "d_pTCandMin");
	printf("Min pt 2prong cand = %f\n", candpTMin);
	candpTMax = GetJsonInteger(esdFileName.Data(), "d_pTCandMax");
	printf("Max pt 2prong cand = %f\n", candpTMax);
	d_maxr = GetJsonFloat(esdFileName.Data(), ", d_maxr");
	printf("Max DCA radius = %f\n", d_maxr);
		
	esdTrackCuts->SetPtRange(ptmintrack, 1.e10);
	esdTrackCuts->SetEtaRange(-0.8, +0.8);
	esdTrackCuts->SetMinNClustersTPC(minncluTPC);
	esdTrackCuts->SetRequireITSRefit(kTRUE);
	esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                         AliESDtrackCuts::kAny);
	esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
	esdTrackCuts->SetMaxDCAToVertexZ(3.2);
	esdTrackCuts->SetMaxDCAToVertexXY(2.4);
	esdTrackCuts->SetDCAToVertex2D(kTRUE);
				
   }else{
   AliError(Form("Json configuration file %s not found\n",esdFileName.Data()));
  }	
}

//_______________________________________________________________________________________
void AliAnalysisTaskHFFindJets::UserExec(Option_t *)
{
  
  AliESDEvent *esd = (AliESDEvent*) (InputEvent());
  if(!esd) {
    printf("AliAnalysisTaskHFFindJets:::UserExec(): bad ESD\n");
    return;
  }
  
  Double_t d03[3] = {0., 0., 0.};
  AliAODRecoDecay* rd4massCalc3 = new AliAODRecoDecay(0x0, 3, 1, d03);
  Double_t covMatrix[6];
  
  for (Int_t iEvent = 0; iEvent < tree->GetEntries(); iEvent++) {
    tree->GetEvent(iEvent);
    
    TString trClass = esd->GetFiredTriggerClasses();
    if (triggerstring != "" && !trClass.Contains(triggerstring))
      continue;
    printf("\n------------ Run %d Event: %d  Tracks %d ------------------\n",
           esd->GetRunNumber(), iEvent, esd->GetNumberOfTracks());
    printf("      Fired Trigger Classes %s\n", trClass.Data());

    Int_t maxTracksToProcess = 9999999; /// temporary to limit the time duration of tests
    Int_t totTracks = TMath::Min(maxTracksToProcess, esd->GetNumberOfTracks());

    AliESDVertex* primvtx = (AliESDVertex*)esd->GetPrimaryVertex();
    TString title = primvtx->GetTitle();
    if (primvtx->IsFromVertexer3D() || primvtx->IsFromVertexerZ())
      continue;
    if (primvtx->GetNContributors() < 2)
      continue;
    hvertexx->Fill(primvtx->GetX());
    hvertexy->Fill(primvtx->GetY());
    hvertexz->Fill(primvtx->GetZ());
    AliAODVertex *vertexAODp = ConvertToAODVertex(primvtx);
    if (triggerstring != "" && !trClass.Contains(triggerstring))
      continue;

    Double_t fBzkG = (Double_t)esd->GetMagneticField();

    
    // Apply single track cuts and flag them
    UChar_t* status = new UChar_t[totTracks];
    for (Int_t iTrack = 0; iTrack < totTracks; iTrack++) {
      status[iTrack] = 0;
      AliESDtrack* track = esd->GetTrack(iTrack);
      hpt_nocuts->Fill(track->Pt());
      htgl_nocuts->Fill(track->GetTgl());
      if (SingleTrkCutsSimple(track, minncluTPC, ptmintrack, dcatoprimxymin, primvtx, fBzkG))
        status[iTrack] = 1; //FIXME
    }

    TObjArray* twoTrackArray = new TObjArray(2);
    TObjArray* threeTrackArray = new TObjArray(3);

    AliVertexerTracks* vt = new AliVertexerTracks(fBzkG);

    Double_t mom0[3], mom1[3], mom2[3];
    for (Int_t iPosTrack_0 = 0; iPosTrack_0 < totTracks; iPosTrack_0++) {
      AliESDtrack* track_p0 = esd->GetTrack(iPosTrack_0);
      track_p0->GetPxPyPz(mom0);
      if (status[iPosTrack_0] == 0)
        continue;
      AliExternalTrackParam* trackext = (AliExternalTrackParam*)track_p0;
      double b[2];
      double bCov[3];
      trackext->PropagateToDCA(primvtx, fBzkG, 100., b, bCov);
      hpt_cuts->Fill(track_p0->Pt());
      hdcatoprimxy_cuts->Fill(b[0]);
      htgl_cuts->Fill(track_p0->GetTgl());
      hitsmap->Fill(track_p0->GetITSClusterMap());
      if (track_p0->Charge() < 0)
        continue;

      for (Int_t iNegTrack_0 = 0; iNegTrack_0 < totTracks; iNegTrack_0++) {
        AliESDtrack* track_n0 = esd->GetTrack(iNegTrack_0);
        track_n0->GetPxPyPz(mom1);
        if (track_n0->Charge() > 0)
          continue;
        if (status[iNegTrack_0] == 0)
          continue;

        twoTrackArray->AddAt(track_p0, 0);
        twoTrackArray->AddAt(track_n0, 1);
        AliESDVertex* trkv = ReconstructSecondaryVertex(vt, twoTrackArray, primvtx, d_maxr);
        if (trkv == 0x0) {
          twoTrackArray->Clear();
          continue;
        }

        //  printf(" px track_0 %.4f, track_1 %.4f \n", TMath::Max(track_0->Px(),track_1->Px()),TMath::Min(track_0->Px(),track_1->Px()));
        hvx->Fill(trkv->GetX());
        hvy->Fill(trkv->GetY());
        hvz->Fill(trkv->GetZ());
        double deltax = trkv->GetX() - primvtx->GetX();
        double deltay = trkv->GetY() - primvtx->GetY();
        double deltaz = trkv->GetZ() - primvtx->GetZ();
        double decaylength = TMath::Sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);
        double decaylengthxy = TMath::Sqrt(deltax * deltax + deltay * deltay);

        AliAODVertex* vertexAOD = ConvertToAODVertex(trkv);
        delete trkv;
        AliAODRecoDecayHF2Prong* the2Prong = Make2Prong(twoTrackArray, vertexAOD, fBzkG);
        the2Prong->SetOwnPrimaryVtx(vertexAODp);

        Int_t twoProngSelection = 3;
        if (selectD0 + selectD0bar > 0)
          twoProngSelection = TwoProngSelectionCuts(the2Prong, candpTMin, candpTMax);
        Double_t m0 = the2Prong->InvMassD0();
        Double_t m0b = the2Prong->InvMassD0bar();
        if (twoProngSelection > 0) {
          if (selectD0 == 0 || twoProngSelection == 1 || twoProngSelection == 3) hmass0->Fill(m0);
          if (selectD0bar == 0 || twoProngSelection == 2 || twoProngSelection == 3) hmass0->Fill(m0b);
          hdecayxyz->Fill(decaylength);
          hdecayxy->Fill(decaylengthxy);
          hptD0->Fill(the2Prong->Pt());
          hptprong0->Fill(the2Prong->PtProng(0));
          hptprong1->Fill(the2Prong->PtProng(1));
          hd0->Fill(the2Prong->Getd0Prong(0));
          hd0->Fill(the2Prong->Getd0Prong(1));
          hd0d0->Fill(the2Prong->Prodd0d0());
          hImpParErr->Fill(the2Prong->Getd0errProng(0));
          hImpParErr->Fill(the2Prong->Getd0errProng(1));
          hDecLenErr->Fill(the2Prong->DecayLengthError());
          hDecLenXYErr->Fill(the2Prong->DecayLengthXYError());
          the2Prong->GetPrimaryVtx()->GetCovMatrix(covMatrix);
          hCovPVXX->Fill(covMatrix[0]);
          the2Prong->GetSecondaryVtx()->GetCovMatrix(covMatrix);
          hCovSVXX->Fill(covMatrix[0]);

	  AliFJWrapper *fFastJetWrapper;
	  fFastJetWrapper = new AliFJWrapper("fFastJetWrapper","fFastJetWrapper");
	  fFastJetWrapper->Clear();
	  fFastJetWrapper->SetR(0.4); 
	  fFastJetWrapper->SetAlgorithm(fastjet::JetAlgorithm::antikt_algorithm);
	  fFastJetWrapper->SetRecombScheme(fastjet::RecombinationScheme::E_scheme);
	  fFastJetWrapper->SetStrategy(fastjet::Strategy::Best);
	  fFastJetWrapper->SetGhostArea(0.005); 
	  fFastJetWrapper->SetAreaType(fastjet::AreaType::passive_area);


	  bool isHFJet=false;
	  fFastJetWrapper->Clear();
	  for (Int_t iTrack = 0; iTrack < totTracks; iTrack++) {
	    AliESDtrack* track = esd->GetTrack(iTrack);
	    if (track->Pt() >= 0.15 && TMath::Abs(track->Eta()) < 0.9){
	      if (iTrack==iNegTrack_0 || iTrack==iPosTrack_0) continue;
	      fFastJetWrapper->AddInputVector(track->Px(), track->Py(), track->Pz(), TMath::Sqrt(track->P()*track->P()+0.13957*0.13957),iTrack+2);
	    }
	  }
	  fFastJetWrapper->AddInputVector(the2Prong->Px(), the2Prong->Py(), the2Prong->Pz(), the2Prong->ED0(),1);

	  fFastJetWrapper->Run();
	  std::vector<fastjet::PseudoJet> jets = fFastJetWrapper->GetInclusiveJets();
	  for (Int_t ijet=0; ijet<jets.size(); ijet++){
	    isHFJet=false;
	    fastjet::PseudoJet jet = jets[ijet];
	    if (jet.pt() < 0.15 || jet.perp() >= 1000.0 || TMath::Abs(jet.eta()) >= 0.5) continue;
	    std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(ijet));
	    for (Int_t iconstituent=0; iconstituent<constituents.size(); iconstituent++){
	      if (constituents[iconstituent].user_index()==1) isHFJet=true;
	      break;
	    }
	    if(!isHFJet) continue;
	    hjetpt->Fill(jet.pt());
	    fastjet::JetDefinition subJetDef(fastjet::JetAlgorithm::cambridge_algorithm , 0.4*2.5,fastjet::RecombinationScheme::E_scheme, fastjet::Best);
	    try{
	      fastjet::ClusterSequence reclusterSeq(constituents, subJetDef);
	      std::vector<fastjet::PseudoJet> reclusteredJet =  reclusterSeq.inclusive_jets(0.0);
	      reclusteredJet = sorted_by_pt(reclusteredJet);
         
	      fastjet::PseudoJet daughterSubJet = reclusteredJet[0];
	      fastjet::PseudoJet parentSubJet1; 
	      fastjet::PseudoJet parentSubJet2;

	      Float_t zg=-1.0,rg=-1.0;
	      Int_t nsd=0;
	      bool softDropped=false;
	      bool isHFSubJet=false;
	      std::vector<fastjet::PseudoJet> constituentsSubJet;
	      while(daughterSubJet.has_parents(parentSubJet1,parentSubJet2)){
		isHFSubJet=false;
		constituentsSubJet=parentSubJet1.constituents();
		for (Int_t iconstituent=0; iconstituent<constituentsSubJet.size(); iconstituent++){
		  if (constituentsSubJet[iconstituent].user_index()==1) isHFSubJet=true;
		}
		if (!isHFSubJet) std::swap(parentSubJet1,parentSubJet2);
		zg=parentSubJet2.perp()/(parentSubJet1.perp()+parentSubJet2.perp());
		rg=parentSubJet1.delta_R(parentSubJet2);

		if (zg >= 0.1*TMath::Power(rg/0.4,0.0)){
		  if(!softDropped){
		    hjetzg->Fill(zg);
		    hjetrg->Fill(rg);
		    softDropped=true;
		  }
		  nsd++;
		}
		daughterSubJet=parentSubJet1;
	      }
	      hjetnsd->Fill(nsd);
	    }catch (fastjet::Error) {}
	  }
	  delete fFastJetWrapper;
        }
        delete the2Prong;
        delete vertexAOD;
        //  printf(" masses = %f %f\n",TMath::Max(m0,m0b),TMath::Min(m0,m0b));
        if (do3Prongs) {
          for (Int_t iPosTrack_1 = iPosTrack_0 + 1; iPosTrack_1 < totTracks; iPosTrack_1++) {
            AliESDtrack* track_p1 = esd->GetTrack(iPosTrack_1);
            if (!track_p1)
              continue;
            if (track_p1->Charge() < 0)
              continue;
            track_p1->GetPxPyPz(mom2);
            if (status[iPosTrack_1] == 0)
              continue;
            // order tracks according to charge: +-+
            threeTrackArray->AddAt(track_p0, 0);
            threeTrackArray->AddAt(track_n0, 1);
            threeTrackArray->AddAt(track_p1, 2);
            Int_t massSel = SelectInvMassAndPt3prong(threeTrackArray, rd4massCalc3);
            if (massSel == 0) {
              threeTrackArray->Clear();
              continue;
            }
            AliESDVertex* trkv3 = ReconstructSecondaryVertex(vt, threeTrackArray, primvtx, d_maxr);
            if (trkv3 == 0x0) {
              threeTrackArray->Clear();
              continue;
            }
            AliAODVertex* vertexAOD3 = ConvertToAODVertex(trkv3);
            AliAODRecoDecayHF3Prong* the3Prong = Make3Prong(threeTrackArray, vertexAOD3, fBzkG);
            //  the3Prong->SetOwnPrimaryVtx(vertexAODp);
            if (massSel & (1 << kbitDplus)) {
              Double_t mp = the3Prong->InvMassDplus();
              hmassP->Fill(mp);
              hvx3->Fill(trkv3->GetX());
              hvy3->Fill(trkv3->GetY());
              hvz3->Fill(trkv3->GetZ());
            }
            delete trkv3;
            delete the3Prong;
            delete vertexAOD3;
            threeTrackArray->Clear();
          }
          for (Int_t iNegTrack_1 = iNegTrack_0 + 1; iNegTrack_1 < totTracks; iNegTrack_1++) {
            AliESDtrack* track_n1 = esd->GetTrack(iNegTrack_1);
            if (!track_n1)
              continue;
            if (track_n1->Charge() > 0)
              continue;
            track_n1->GetPxPyPz(mom2);
            if (status[iNegTrack_1] == 0)
              continue;
            // order tracks according to charge: -+-
            threeTrackArray->AddAt(track_n0, 0);
            threeTrackArray->AddAt(track_p0, 1);
            threeTrackArray->AddAt(track_n1, 2);
            Int_t massSel = SelectInvMassAndPt3prong(threeTrackArray, rd4massCalc3);
            if (massSel == 0) {
              threeTrackArray->Clear();
              continue;
            }
            AliESDVertex* trkv3 = ReconstructSecondaryVertex(vt, threeTrackArray, primvtx, d_maxr);
            if (trkv3 == 0x0) {
              threeTrackArray->Clear();
              continue;
            }
            AliAODVertex* vertexAOD3 = ConvertToAODVertex(trkv3);
            AliAODRecoDecayHF3Prong* the3Prong = Make3Prong(threeTrackArray, vertexAOD3, fBzkG);
            //  the3Prong->SetOwnPrimaryVtx(vertexAODp);
            if (massSel & (1 << kbitDplus)) {
              Double_t mp = the3Prong->InvMassDplus();
              hmassP->Fill(mp);
            }
            delete trkv3;
            delete the3Prong;
            delete vertexAOD3;
            threeTrackArray->Clear();
          }
        }
        twoTrackArray->Clear();
      }
      //  delete vertexAODp;
    }
    delete[] status;
    delete vt;
    delete twoTrackArray;
    delete threeTrackArray;
  }
  delete esdTrackCuts;
  delete rd4massCalc3;
  
  PostData(1,fOutput);
}

//_______________________________________________________________________________________
void AliAnalysisTaskHFFindJets::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }
  return;
}
