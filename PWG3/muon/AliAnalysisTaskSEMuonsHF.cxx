#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TClonesArray.h>

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliAODMCParticle.h"
#include "AliVVertex.h"
#include "AliVHeader.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDMuonTrack.h"
#include "AliAODTrack.h"
#include "AliMuonsHFHeader.h"
#include "AliAODMuonTrack.h"
#include "AliAODMuonPair.h"
#include "AliMCMuonTrack.h"
#include "AliMCMuonPair.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEMuonsHF.h"

ClassImp(AliAnalysisTaskSEMuonsHF)

//_____________________________________________________________________________
AliAnalysisTaskSEMuonsHF::AliAnalysisTaskSEMuonsHF() :
AliAnalysisTaskSE(),
fAnaMode(0),
fIsOutputTree(kFALSE),
fIsUseMC(kFALSE),
fHeader(0),
fMuTrkClArr(0),
fMuPairClArr(0),
fListHisAtEvLevel(0),
fListHisSingleMuon(0),
fListHisDimuon(0)
{
  //
  // Default constructor
  //
  fSingleMuonCuts[0] = -1.;
  fSingleMuonCuts[1] = 999999.;
  fSingleMuonCuts[2] = -1.;
  fSingleMuonCuts[3] = 999999.;
  fSingleMuonCuts[4] = -1.;
  fSingleMuonCuts[5] = 999999.;
  fSingleMuonCuts[6] = -1.;
  fSingleMuonCuts[7] = 999999.;
  fSingleMuonCuts[8] = -1.;
  fSingleMuonCuts[9] = 4.;
}

//_____________________________________________________________________________
AliAnalysisTaskSEMuonsHF::AliAnalysisTaskSEMuonsHF(const char *name) :
AliAnalysisTaskSE(name),
fAnaMode(0),
fIsOutputTree(kFALSE),
fIsUseMC(kFALSE),
fHeader(0),
fMuTrkClArr(0),
fMuPairClArr(0),
fListHisAtEvLevel(0),
fListHisSingleMuon(0),
fListHisDimuon(0)
{
  //
  // Constructor
  //
  fSingleMuonCuts[0] = -1.;
  fSingleMuonCuts[1] = 999999.;
  fSingleMuonCuts[2] = -1.;
  fSingleMuonCuts[3] = 999999.;
  fSingleMuonCuts[4] = -1.;
  fSingleMuonCuts[5] = 999999.;
  fSingleMuonCuts[6] = -1.;
  fSingleMuonCuts[7] = 999999.;
  fSingleMuonCuts[8] = -1.;
  fSingleMuonCuts[9] = 4.;

  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskSEMuonsHF::~AliAnalysisTaskSEMuonsHF()
{
  //
  // Default destructor
  //
  if (fHeader) { delete fHeader; fHeader=NULL; }
  if (fMuTrkClArr)  { delete fMuTrkClArr;  fMuTrkClArr =NULL; }
  if (fMuPairClArr) { delete fMuPairClArr; fMuPairClArr=NULL; }

  if (fListHisSingleMuon) { delete fListHisSingleMuon; fListHisSingleMuon=NULL; }
  if (fListHisDimuon)     { delete fListHisDimuon;     fListHisDimuon=NULL;     }
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::UserCreateOutputObjects()
{
  CreateOutoutHistosAtEvnetLevel();
  if (fAnaMode!=2) CreateOutputHistosSingleMuon();
  if (fAnaMode!=1) CreateOutputHistosDimuon();

  if (!fHeader) {
    fHeader = new AliMuonsHFHeader();
    fHeader->SetName("MuonsHFHeader");
  }
  if (!fMuTrkClArr) {
    if (fIsUseMC) fMuTrkClArr = new TClonesArray("AliMCMuonTrack", 0);
    else fMuTrkClArr = new TClonesArray("AliAODMuonTrack", 0);
    fMuTrkClArr->SetName("MuonTrack");
  }

  if (fAnaMode!=1 && !fMuPairClArr) {
    if (fIsUseMC) fMuPairClArr = new TClonesArray("AliMCMuonPair", 0);
    else fMuPairClArr = new TClonesArray("AliAODMuonPair", 0);
    fMuPairClArr->SetName("MuonPair");
  }

  if (fIsOutputTree) {
    AddAODBranch("AliMuonsHFHeader", &fHeader);
    AddAODBranch("TClonesArray", &fMuTrkClArr);
    if (fAnaMode!=1) AddAODBranch("TClonesArray", &fMuPairClArr);
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::UserExec(Option_t *)
{
  AliVEvent *event = dynamic_cast<AliVEvent*>(InputEvent());
  //if (AODEvent() && IsStandardAOD()) event = dynamic_cast<AliVEvent*>(AODEvent());
  TString evName = event->IsA()->GetName();
  Bool_t isAOD = ((evName=="AliAODEvent") ? kTRUE : kFALSE);
  fHeader->SetHeader(((AliVHeader*)event->GetHeader()));
  fHeader->SetVertex(((AliVVertex*)event->GetPrimaryVertex()));
  event = 0x0;

  Int_t ntrks = 0;
  AliAODEvent *aod = 0;
  AliESDEvent *esd = 0;
  TClonesArray *mcClArr = 0;
  AliMCEventHandler *mcH = 0;
  if (isAOD) {
    aod = dynamic_cast<AliAODEvent*>(InputEvent());
    //if (AODEvent() && IsStandardAOD()) aod = dynamic_cast<AliAODEvent*>(AODEvent());
    if (!aod) {
      AliError("AOD event not found. Nothing done!");
      return;
    }
    if (fIsUseMC) {
      mcClArr = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
      if (!mcClArr) {
        AliError("MC Array not found. Nothing done!");
        return;
      }
    }
    ntrks = aod->GetNTracks();
  } else {
    esd = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!esd) {
      AliError("ESD event not found. Nothing done!");
      return;
    }
    if (fIsUseMC) {
      mcH = (AliMCEventHandler*)
              ((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
      if (!mcH) {
        AliError("MC Handler not found. Nothing done!");
        return;
      }
    }
    ntrks = esd->GetNumberOfMuonTracks();
  }

  fMuTrkClArr->Delete();
  TClonesArray &muTrkRef = *fMuTrkClArr;
  Int_t countN = fMuTrkClArr->GetEntriesFast();

  Int_t theMult = 0;
  AliAODTrack *trkAOD = 0;
  AliESDMuonTrack *trkESD = 0;
  AliAODMuonTrack *track = 0;
  AliMCMuonTrack *trkMC = 0;
  for (Int_t itrk=0; itrk<ntrks; itrk++) {  // loop over all tracks
    track=0; trkMC=0;
    if (isAOD) {
      trkAOD = (AliAODTrack*)aod->GetTrack(itrk);
      if (!trkAOD->IsMuonTrack()) { trkAOD=0; continue; }
      if (fIsUseMC) trkMC = new AliMCMuonTrack(trkAOD, mcClArr);
      else track = new AliAODMuonTrack(trkAOD);
    } else {
      trkESD = (AliESDMuonTrack*)esd->GetMuonTrack(itrk);
      if (!trkESD->ContainTrackerData()) { trkESD=0; continue; }
      if (fIsUseMC) trkMC = new AliMCMuonTrack(trkESD, esd, mcH);
      else track = new AliAODMuonTrack(trkESD);
    }

    if (fAnaMode!=2) {
      if (track) FillDistributionsSingleMuon(track);
      else FillDistributionsSingleMuon((AliAODMuonTrack*)trkMC, trkMC->GetSource());
    }

    if (fIsOutputTree || fAnaMode!=1) {
      if (track) new(muTrkRef[countN++]) AliAODMuonTrack(*track);
      else new(muTrkRef[countN++]) AliMCMuonTrack(*trkMC);
    }

    theMult++;
    trkESD = 0;
    trkAOD = 0;
    if (track) { delete track; track=0; }
    if (trkMC) { delete trkMC; trkMC=0; }
  }  // end loop of all tracks
  fHeader->SetMultSingleMuon(theMult);
  FillDistributionsAtEventLeavel();
  PostData(1, fListHisAtEvLevel);
  if (fAnaMode!=2) PostData(2, fListHisSingleMuon);
  if (fAnaMode==1) return;

  theMult = 0;
  fMuPairClArr->Delete();
  TClonesArray &muPairRef = *fMuPairClArr;
  countN = fMuPairClArr->GetEntriesFast();
  AliAODMuonPair *pair = 0;
  AliMCMuonPair *pairMC = 0;
  ntrks = fMuTrkClArr->GetEntriesFast();
  for (Int_t itrk=0; itrk<ntrks-1; itrk++) {  // 1st loop over muon tracks
    for (Int_t jtrk=itrk+1; jtrk<ntrks; jtrk++) {  // 2nd loop over muon tracks
      pair=0; pairMC=0;
      if (fIsUseMC)
        pairMC = new AliMCMuonPair((AliMCMuonTrack*)fMuTrkClArr->At(itrk), (AliMCMuonTrack*)fMuTrkClArr->At(jtrk));
      else
        pair = new AliAODMuonPair((AliAODMuonTrack*)fMuTrkClArr->At(itrk), (AliAODMuonTrack*)fMuTrkClArr->At(jtrk));

      if (pair) FillDistributionsDimuon(pair);
      else FillDistributionsDimuon(pairMC, pairMC->GetSource());

      if (fIsOutputTree) {
        if (pair) new(muPairRef[countN++]) AliAODMuonPair(*pair);
        else new(muPairRef[countN++]) AliMCMuonPair(*pairMC);
      }

      theMult++;
      if (pair)   { delete pair;   pair  =0; }
      if (pairMC) { delete pairMC; pairMC=0; }
    }  // end 2nd loop of muon traks
  }  // end 1st loop of muon tracks
  fHeader->SetMultDimuon(theMult);

  PostData(3, fListHisDimuon);
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::Terminate(Option_t *)
{
  printf("This is Terminate!\n");
  // adding the fitting
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::CreateOutoutHistosAtEvnetLevel()
{
  if (!fListHisAtEvLevel) fListHisAtEvLevel = new TList();
  Int_t nbins[kNHistEv];
  Double_t xlow[kNHistEv], xup[kNHistEv];
  TString name[kNHistEv], axis[kNHistEv], unit[kNHistEv];
  nbins[kHistVx]=200; xlow[kHistVx]  =-50.; xup[kHistVx]  =50.; name[kHistVx]  =  "vx"; axis[kHistVx] ="v_{x}"; unit[kHistVx]="cm";
  nbins[kHistVy]=200; xlow[kHistVy]  =-50.; xup[kHistVy]  =50.; name[kHistVy]  =  "vy"; axis[kHistVy] ="v_{y}"; unit[kHistVy]="cm";
  nbins[kHistVz]=200; xlow[kHistVz]  =-50.; xup[kHistVz]  =50.; name[kHistVz]  =  "vz"; axis[kHistVz] ="v_{z}"; unit[kHistVz]="cm";
  nbins[kHistMult]=9; xlow[kHistMult]=-0.5; xup[kHistMult]=8.5; name[kHistMult]="mult"; axis[kHistMult]="mult"; unit[kHistMult]="";

  TH1F *histo = 0;
  TString histName, histTitle;
  for (Int_t ihist=0; ihist<kNHistEv; ihist++) {  // loop over histos
    histName  = Form("hsMu%s", name[ihist].Data());
    histTitle = Form("%s of event", axis[ihist].Data());
    histo = new TH1F(histName.Data(), histTitle.Data(), nbins[ihist], xlow[ihist], xup[ihist]);
    histo->GetXaxis()->SetTitle(Form("%s %s", axis[ihist].Data(), unit[ihist].Data()));
    histo->Sumw2();
    fListHisAtEvLevel->AddAt(histo, ihist);
  }  // end loop over histos

  return;
}


//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::CreateOutputHistosSingleMuon()
{
  if (!fListHisSingleMuon) fListHisSingleMuon = new TList();

  Int_t nbins[kNHistMu];
  Double_t xlow[kNHistMu], xup[kNHistMu];
  TString name[kNHistMu], axis[kNHistMu], unit[kNHistMu];
  nbins[kHistP]  =300 ; xlow[kHistP]=   0. ; xup[kHistP]  =150.; name[kHistP]  ="P"  ; axis[kHistP]  ="P"      ; unit[kHistP]  ="[GeV/c]";
  nbins[kHistPt] =120 ; xlow[kHistPt]=  0. ; xup[kHistPt] = 30.; name[kHistPt] ="Pt" ; axis[kHistPt] ="p_{t}"  ; unit[kHistPt] ="[GeV/c]";
  nbins[kHistEta]=50  ; xlow[kHistEta]=-4. ; xup[kHistEta]=-2.5; name[kHistEta]="Eta"; axis[kHistEta]="#eta"   ; unit[kHistEta]="";
  nbins[kHistDca]=1000; xlow[kHistDca]= 0. ; xup[kHistDca]=500.; name[kHistDca]="Dca"; axis[kHistDca]="DCA"    ; unit[kHistDca]="cm";
  nbins[kHistTrg]=   4; xlow[kHistTrg]=-0.5; xup[kHistTrg]= 3.5; name[kHistTrg]="trg"; axis[kHistTrg]="Trigger"; unit[kHistTrg]="";

  TH1F *histo = 0;
  TString histName, histTitle;
  Int_t index=-1, lowBound=0;
  for (Int_t ihist=0; ihist<kNHistMu; ihist++) {  // loop over histos
    index = lowBound + ihist;
    histName  = Form("hsMu%s", name[ihist].Data());
    histTitle = Form("%s of muon candidates", axis[ihist].Data());
    histo = new TH1F(histName.Data(), histTitle.Data(), nbins[ihist], xlow[ihist], xup[ihist]);
    histo->GetXaxis()->SetTitle(Form("%s %s", axis[ihist].Data(), unit[ihist].Data()));
    histo->Sumw2();
    fListHisSingleMuon->AddAt(histo, index);
  }  // end loop over histos

  if (fIsUseMC) CreateOutputHistosSingleMuonMC(nbins, xlow, xup, name, axis, unit);
  TH2F *hVzPt = new TH2F("hVzPt", "v_{z} vs. p_{t}", 300, -30., 30., 120, 0., 30);
  hVzPt->GetXaxis()->SetTitle("v_{z} [cm]");
  hVzPt->GetYaxis()->SetTitle("p_{t} [GeV/c]");
  fListHisSingleMuon->AddLast(hVzPt);
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::CreateOutputHistosDimuon()
{
  if (!fListHisDimuon) fListHisDimuon = new TList();

  TString dimuName[kNMuMus], dimuTitle[kNMuMus];
  dimuName[kMuNMuP]="MuNMuP"; dimuTitle[kMuNMuP]="#mu^{-}#mu^{+}";
  dimuName[kMuNMuN]="MuNMuN"; dimuTitle[kMuNMuN]="#mu^{-}#mu^{-}";
  dimuName[kMuPMuP]="MuPMuP"; dimuTitle[kMuPMuP]="#mu^{+}#mu^{+}";

  Int_t nbins[kNHistDimu];
  Double_t xlow[kNHistDimu], xup[kNHistDimu];
  TString name[kNHistDimu], axis[kNHistDimu], unit[kNHistDimu];
  nbins[kInvM]  =300; xlow[kInvM]  =0.; xup[kInvM]  =30.; name[kInvM]  ="InvM"; axis[kInvM]  ="M";     unit[kInvM]  ="[GeV/c^{2}]";
  nbins[kPtPair]=120; xlow[kPtPair]=0.; xup[kPtPair]=60.; name[kPtPair]="Pt";   axis[kPtPair]="p_{t}"; unit[kPtPair]="[GeV/c]";

  TH1F *histo = 0;
  TString histName, histTitle;
  Int_t index=-1, lowBound=0;
  for (Int_t idmu=0; idmu<kNMuMus; idmu++) {  // loop over different kinds of dimuon
    for (Int_t ihist=0; ihist<kNHistDimu; ihist++) {  // loop over histos
      index = lowBound + ihist + idmu*kNHistDimu;
      histName  = Form("h%s%s", dimuName[idmu].Data(), name[ihist].Data());
      histTitle = Form("%s of %s candidates", axis[ihist].Data(), dimuTitle[idmu].Data());
      histo = new TH1F(histName.Data(), histTitle.Data(), nbins[ihist], xlow[ihist], xup[ihist]);
      histo->GetXaxis()->SetTitle(Form("%s(%s) %s", axis[ihist].Data(), dimuTitle[idmu].Data(), unit[ihist].Data()));
      histo->Sumw2();
      fListHisDimuon->AddAt(histo, index);
    }  // end loop over histos
  }  // end loop of different kinds of dimuon

  if (fIsUseMC) CreateOutputHistosDimuonMC(nbins, xlow, xup, name, axis, unit, dimuName, dimuTitle);
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::CreateOutputHistosSingleMuonMC(Int_t *nbins, Double_t *xlow, Double_t *xup,
                                                              TString *name, TString *axis, TString *unit)
{
  // define the histos of single muons
  TString srcName[kNSingleMuSrcs];
  srcName[kBeautyMu]    = "Bottom";
  srcName[kCharmMu]     = "Charm";
  srcName[kPrimaryMu]   = "Primary";
  srcName[kSecondaryMu] = "Secondary";
  srcName[kNotMu]       = "NotMu";

  TH1F *histo = 0;
  TString histName, histTitle;
  Int_t index=-1, lowBound=kNHistMu;
  for (Int_t src=0; src<kNSingleMuSrcs; src++) {  // loop over single muon src
    for (Int_t ihist=0; ihist<kNHistMu; ihist++) {  // loop over histos
      index = lowBound + ihist + src*kNHistMu;
      histName  = Form("hsMu%s_%s", srcName[src].Data(), name[ihist].Data());
      histTitle = Form("%s of muon candidates<-%s", axis[ihist].Data(), srcName[src].Data());
      histo = new TH1F(histName.Data(), histTitle.Data(), nbins[ihist], xlow[ihist], xup[ihist]);
      histo->GetXaxis()->SetTitle(Form("%s %s", axis[ihist].Data(), unit[ihist].Data()));
      histo->Sumw2();
      fListHisSingleMuon->AddAt(histo, index);
    }  // end loop over histos
  }  // // end loop of single muon src

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::CreateOutputHistosDimuonMC(Int_t *nbins, Double_t *xlow, Double_t *xup,
                                                          TString *name, TString *axis, TString *unit,
                                                          TString *dimuName, TString *dimuTitle)
{
  TString srcName[kNDimuSrcs];
  srcName[kBBdiff]    = "BBdiff";
  srcName[kBchain]    = "Bchain";
  srcName[kDDdiff]    = "DDdiff";
  srcName[kDchain]    = "Dchain";
  srcName[kResonance] = "Resonance";
  srcName[kUncorr]    = "Uncorr";

  TH1F *histo = 0;
  TString histName, histTitle;
  Int_t index=-1, lowBound=kNMuMus*kNHistDimu;
  for (Int_t idmu=0; idmu<kNMuMus; idmu++) {  // loop over different kinds of dimuon
    for (Int_t ihist=0; ihist<kNHistDimu; ihist++) {  // loop over histos
      for (Int_t src=0; src<kNDimuSrcs; src++) {  // loop over dimuon src
        index = lowBound + src + ihist*kNDimuSrcs + idmu*kNHistDimu*kNDimuSrcs;
        histName  = Form("hs%s%s_%s", dimuName[idmu].Data(), srcName[src].Data(), name[ihist].Data());
        histTitle = Form("%s of %s<-%s", name[ihist].Data(), dimuName[idmu].Data(), srcName[src].Data());
        histo = new TH1F(histName.Data(), histTitle.Data(), nbins[ihist], xlow[ihist], xup[ihist]);
        histo->GetXaxis()->SetTitle(Form("%s(%s) %s", axis[ihist].Data(), dimuTitle[idmu].Data(), unit[ihist].Data()));
        fListHisDimuon->AddAt(histo, index);
      }  // end loop of dimuon src
    }  // end loop over histos
  }  // end loop of different kinds of dimuon

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::FillDistributionsAtEventLeavel()
{
  Double_t disEv[kNHistEv];
  disEv[kHistVx]   = fHeader->GetXv();
  disEv[kHistVy]   = fHeader->GetYv();
  disEv[kHistVz]   = fHeader->GetZv();
  disEv[kHistMult] = (Double_t)fHeader->GetMultSingleMuon();

  Int_t index=-1, lowBound=0;
  for (Int_t ihist=0; ihist<kNHistEv; ihist++) {
    index = lowBound + ihist;
    ((TH1F*)fListHisAtEvLevel->At(index))->Fill(disEv[ihist]);
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::FillDistributionsSingleMuon(AliAODMuonTrack *track, Int_t src)
{
  if (fHeader->IsUnrecoVertex()) return;
  if (!(track->SelectSingleMuon(fSingleMuonCuts))) return;

  TLorentzVector lorentzP = track->GetP();
  Double_t disMu[kNHistMu];
  disMu[kHistP] = lorentzP.P();
  disMu[kHistPt] = lorentzP.Pt();
  disMu[kHistEta] = lorentzP.Eta();
  disMu[kHistDca] = track->GetDCA();
  disMu[kHistTrg] = (Double_t)track->GetTrigger();

  Int_t index=-1, lowBound=0;
  for (Int_t ihist=0; ihist<kNHistMu; ihist++) {
    index = lowBound + ihist;
    ((TH1F*)fListHisSingleMuon->At(index))->Fill(disMu[ihist]);
  } 

  if (fIsUseMC && src>=0) FillDistributionsSingleMuonMC(disMu, src);

  Double_t vtx[3];
  fHeader->GetXYZ(vtx);
  ((TH2F*)fListHisSingleMuon->FindObject("hVzPt"))->Fill(vtx[2], disMu[kHistPt]);
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::FillDistributionsSingleMuonMC(Double_t *disMu, Int_t src)
{
  if (src<0) return;
  Int_t index=-1, lowBound=kNHistMu;
  for (Int_t ihist=0; ihist<kNHistMu; ihist++) {
    index = lowBound + ihist + src*kNHistMu;
    ((TH1F*)fListHisSingleMuon->At(index))->Fill(disMu[ihist]);
  }
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::FillDistributionsDimuon(AliAODMuonPair *pair, Int_t src)
{
  if (!(pair->GetTrack(0)->SelectSingleMuon(fSingleMuonCuts)) ||
      !(pair->GetTrack(1)->SelectSingleMuon(fSingleMuonCuts))) return;

  TLorentzVector lorentzP = pair->GetP();
  Double_t disDimu[kNHistDimu];
  disDimu[kInvM]   = lorentzP.Mag();
  disDimu[kPtPair] = lorentzP.Pt();

  Int_t dimuK=-1, charge=pair->GetCharge();
  if (charge==0) dimuK = kMuNMuP;
  if (charge<0)  dimuK = kMuNMuN;
  if (charge>0)  dimuK = kMuPMuP;

  Int_t index=-1, lowBound=0;
  for (Int_t ihist=0; ihist<kNHistDimu; ihist++) {
    index = lowBound + ihist + dimuK*kNHistDimu;
    ((TH1F*)fListHisDimuon->At(index))->Fill(disDimu[ihist]);
  }

  if (fIsUseMC && src>=0) FillDistributionsDimuonMC(disDimu, dimuK, src);
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEMuonsHF::FillDistributionsDimuonMC(Double_t *disDimu, Int_t dimuK, Int_t src)
{
  if (src<0) return;
  Int_t index=-1, lowBound=kNMuMus*kNHistDimu;
  for (Int_t ihist=0; ihist<kNHistDimu; ihist++) {
    index = lowBound + src + ihist*kNDimuSrcs + dimuK*kNHistDimu*kNDimuSrcs;
    ((TH1F*)fListHisDimuon->At(index))->Fill(disDimu[ihist]);
  }
  return;
}
