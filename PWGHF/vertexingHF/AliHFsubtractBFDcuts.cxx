/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////////////////////////////////
/// AliHFsubtractBFDcuts                                                  //
///                                                                       //
///  Class for storing and handling D0 meson candidates properties        //
///  for estimating the feed-down fraction using several sets of cuts     //
///     Andrea Rossi <andrea.rossi@cern.ch>                               //
///     Felix Reidt  <felix.reidt@cern.ch>                                //
///                                                                       //
////////////////////////////////////////////////////////////////////////////
#include "AliHFsubtractBFDcuts.h"

#include <vector>

#include "TNamed.h"
#include "THnSparse.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TObjString.h"
#include "TMath.h"

#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliVertexerTracks.h"
#include "AliExternalTrackParam.h"
#include "AliNeutralTrackParam.h"

/// \cond CLASSIMP
ClassImp(AliHFsubtractBFDcuts);
/// \endcond

AliHFsubtractBFDcuts::AliHFsubtractBFDcuts()
  : TNamed()
  , fIsMC(kTRUE)
  , fCheckAcceptance(kTRUE)
  , fResolveResonances(kTRUE)
  , fTHnGenStep(0x0)
  , fTHnData(0x0)
  , fTHnMC(0x0)
  , fMCarray(0x0)
  , fAODtracks(0x0)
  , fLabCand(-1)
  , fPtCand(0.)
  , fNprongs((UInt_t)-1)
  , fNprongsInAcc((UInt_t)-1)
  , fFoundElectron(kFALSE)
  , fMotherPt(-1.)
  , fGenerateDecayList(kFALSE) // DON'T ACTIVATE, DOESN'T MERGE
  , fDecayProngs()
  , fDecayStrList(0x0)
  , fQAhists(0x0)
{
  /// default constructor
}

AliHFsubtractBFDcuts::AliHFsubtractBFDcuts(const char* name, const char* title)
  : TNamed(name,title)
  , fIsMC(kTRUE)
  , fCheckAcceptance(kTRUE)
  , fResolveResonances(kTRUE)
  , fTHnGenStep(0x0)
  , fTHnData(0x0)
  , fTHnMC(0x0)
  , fMCarray(0x0)
  , fAODtracks(0x0)
  , fLabCand(-1)
  , fPtCand(0.)
  , fNprongs((UInt_t)-1)
  , fNprongsInAcc((UInt_t)-1)
  , fFoundElectron(kFALSE)
  , fMotherPt(-1.)
  , fGenerateDecayList(kFALSE) // DON'T ACTIVATE, DOESN'T MERGE
  , fDecayProngs()
  , fDecayStrList(0x0)
  , fQAhists(0x0)
{
  /// default constructor with name
  fDecayStrList = new TList();
  fDecayStrList->SetOwner();
  fQAhists = new TList();
  fQAhists->SetOwner();
}

AliHFsubtractBFDcuts::AliHFsubtractBFDcuts(const AliHFsubtractBFDcuts& c)
  : TNamed(c.GetName(),c.GetTitle())
  , fIsMC(c.fIsMC)
  , fCheckAcceptance(c.fCheckAcceptance)
  , fResolveResonances(c.fResolveResonances)
  , fTHnGenStep(c.fTHnGenStep)
  , fTHnData(c.fTHnData)
  , fTHnMC(c.fTHnMC)
  , fMCarray(c.fMCarray)
  , fAODtracks(c.fAODtracks)
  , fLabCand(c.fLabCand)
  , fPtCand(c.fPtCand)
  , fNprongs(c.fNprongs)
  , fNprongsInAcc(c.fNprongsInAcc)
  , fFoundElectron(c.fFoundElectron)
  , fMotherPt(c.fMotherPt)
  , fGenerateDecayList(c.fGenerateDecayList)
  , fDecayProngs(c.fDecayProngs)
  , fDecayStrList(c.fDecayStrList) // FIXME: TList copy contructor not implemented
  , fQAhists(c.fQAhists)
{
  /// copy constructor
}

AliHFsubtractBFDcuts AliHFsubtractBFDcuts::operator=(const AliHFsubtractBFDcuts& c)
{
  /// assignment operator
  return c;
}

AliHFsubtractBFDcuts::~AliHFsubtractBFDcuts() {
  AliDebug(3, "TODO: implement me!");
}

void AliHFsubtractBFDcuts::InitHistos(){
  // mass, pt, normLXY, cosPointXY, normL, cosPoint, LXY, L
  Int_t dimAxes[8]={500  ,96 ,200 ,100   , 200,100   , 100 ,100  };
  Double_t  min[8]={  1.7, 0.,  0.,  0.99,  0.,  0.99,  0. ,  0. };
  Double_t  max[8]={  2.2,24.,100.,  1.  ,100.,  1.  ,  1.0,  1.0};
  fTHnData=new THnSparseF("fCutsDataFD","fCutsDataFD",8,dimAxes,min,max);
  fTHnData->GetAxis(0)->SetName("mass");
  fTHnData->GetAxis(0)->SetTitle("Mass (K,#pi) (GeV/#it{c^{2}})");
  fTHnData->GetAxis(1)->SetName("pt");
  fTHnData->GetAxis(1)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fTHnData->GetAxis(2)->SetName("NormDecLengthXY");
  fTHnData->GetAxis(2)->SetTitle("Normalized XY decay length");
  fTHnData->GetAxis(3)->SetName("CosPointXY");
  fTHnData->GetAxis(3)->SetTitle("Cos#theta_{point}^{XY}");
  fTHnData->GetAxis(4)->SetName("NormDecLength");
  fTHnData->GetAxis(4)->SetTitle("Normalized decay length");
  fTHnData->GetAxis(5)->SetName("CosPoint");
  fTHnData->GetAxis(5)->SetTitle("Cos#theta_{point}");
  fTHnData->GetAxis(6)->SetName("DecLengthXY");
  fTHnData->GetAxis(6)->SetTitle("XY decay length (cm)");
  fTHnData->GetAxis(7)->SetName("DecLength");
  fTHnData->GetAxis(7)->SetTitle("Decay length (cm)");

  // pt, normLXY, cosPointXY, #prongs, mother pt, normL, cosPoint, LXY, L
  Int_t dimAxesMC[10]={96 ,200 ,100   ,20 ,24 ,200 ,100   ,100  ,100  ,2 };
  Double_t  minMC[10]={ 0.,  0.,  0.99, 0., 0.,  0.,  0.99,  0. ,  0. ,0.};
  Double_t  maxMC[10]={24.,100.,  1.  ,20.,24.,100.,  1.  ,  1.0,  1.0,2.};
  fTHnMC=new THnSparseF("fCutsMCFD","fCutsMCFD",10,dimAxesMC,minMC,maxMC);
  fTHnMC->GetAxis(0)->SetName("pt");
  fTHnMC->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fTHnMC->GetAxis(1)->SetName("NormDecLengthXY");
  fTHnMC->GetAxis(1)->SetTitle("Normalized XY decay length");
  fTHnMC->GetAxis(2)->SetName("CosPointXY");
  fTHnMC->GetAxis(2)->SetTitle("Cos#theta_{point}^{XY}");
  fTHnMC->GetAxis(3)->SetName("nProngs");
  fTHnMC->GetAxis(3)->SetTitle("Number of Decay Prongs");
  fTHnMC->GetAxis(4)->SetName("Mother_pt");
  fTHnMC->GetAxis(4)->SetTitle("Mother #it{p}_{T} (GeV/#it{c})");
  fTHnMC->GetAxis(5)->SetName("NormDecLength");
  fTHnMC->GetAxis(5)->SetTitle("Normalized decay length");
  fTHnMC->GetAxis(6)->SetName("CosPoint");
  fTHnMC->GetAxis(6)->SetTitle("Cos#theta_{point}");
  fTHnMC->GetAxis(7)->SetName("DecLengthXY");
  fTHnMC->GetAxis(7)->SetTitle("XY decay length (cm)");
  fTHnMC->GetAxis(8)->SetName("DecLength");
  fTHnMC->GetAxis(8)->SetTitle("Decay length (cm)");
  fTHnMC->GetAxis(9)->SetName("ContainsElectron");
  fTHnMC->GetAxis(9)->SetTitle("ContainsElectron");

  Int_t dimAxesGen[6]={96 ,20 ,24 ,100 ,100 ,2 };
  Double_t  minGen[6]={ 0., 0., 0.,  0.,  0.,0.};
  Double_t  maxGen[6]={24.,20.,24.,  1.,  1.,2.};
  fTHnGenStep=new THnSparseF("fPtMCGenStep","fPtMCGenStep",6,dimAxesGen,minGen,maxGen);
  fTHnGenStep->GetAxis(0)->SetName("pt");
  fTHnGenStep->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fTHnGenStep->GetAxis(1)->SetName("nProngs");
  fTHnGenStep->GetAxis(1)->SetTitle("Number of Decay Prongs");
  fTHnGenStep->GetAxis(2)->SetName("Mother_pt");
  fTHnGenStep->GetAxis(2)->SetTitle("Mother #it{p}_{T} (GeV/#it{c})");
  fTHnGenStep->GetAxis(3)->SetName("DecLengthXY");
  fTHnGenStep->GetAxis(3)->SetTitle("XY decay length (cm)");
  fTHnGenStep->GetAxis(4)->SetName("DecLength");
  fTHnGenStep->GetAxis(4)->SetTitle("Decay length (cm)");
  fTHnGenStep->GetAxis(5)->SetName("ContainsElectron");
  fTHnGenStep->GetAxis(5)->SetTitle("ContainsElectron");

  fQAhists->Add(new TH1F("hRapidityDist"           , "All Particles;y;Counts (a.u.)"                                                                                                ,800,-4., 4.                            )); //  0
  fQAhists->Add(new TH1F("hRapidityDistStable"     , "All Stable Particles;y;Counts (a.u.)"                                                                                         ,800,-4., 4.                            )); //  1
  fQAhists->Add(new TH2F("hDecayLengthB"           , ";D meson #it{p}_{T};Decay Length B Meson (cm);counts (a.u.)"                                                                  , 24, 0.,24.  , 200,0.,2.0              )); //  2
  fQAhists->Add(new TH2F("hDecayLengthBxy"         , ";D meson #it{p}_{T};Decay Length B Meson in XY direction (cm);counts (a.u.)"                                                  , 24, 0.,24.  , 200,0.,2.0              )); //  3
  fQAhists->Add(new TH3F("hDCABhypothesis"         , "All Particles;D meson #it{p}_{T};Distance of Closest Approach (cm);Impact parameter of D0 x track (cm^{2})"                   , 24, 0.,24.  , 100,0.,1.0 ,100,0.,0.001)); //  4
  fQAhists->Add(new TH3F("hDCABprongs"             , "Only Real B Decay Prongs;D meson #it{p}_{T};Distance of Closest Approach (cm);Impact parameter of D0 x track (cm^{2})"        , 24, 0.,24.  , 100,0.,1.0 ,100,0.,0.001)); //  5
  fQAhists->Add(new TH3F("hVtxDistBprongs"         , "Only Real B Decay Prongs;D meson #it{p}_{T};Distance to Primary Vertex (cm);Error Distance to Primary Vertex (cm)"            , 24, 0.,24.  , 100,0.,1.0 ,100,0.,1.0  )); //  6
  fQAhists->Add(new TH3F("hVtxDistXYBprongs"       , "Only Real B Decay Prongs;D meson #it{p}_{T};Distance to Primary Vertex in XY (cm);Error Distance to Primary Vertex in XY (cm)", 24, 0.,24.  , 100,0.,1.0 ,100,0.,1.0  )); //  7
  fQAhists->Add(new TH3F("hVtxDistBhypothesis"     , "All Particles;D meson #it{p}_{T};Distance to Primary Vertex (cm);Error Distance to Primary Vertex (cm)"                       , 24, 0.,24.  , 100,0.,1.0 ,100,0.,1.0  )); //  8
  fQAhists->Add(new TH3F("hVtxDistXYBhypthesis"    , "All Particles;D meson #it{p}_{T};Distance to Primary Vertex in XY (cm);Error Distance to Primary Vertex in XY (cm)"           , 24, 0.,24.  , 100,0.,1.0 ,100,0.,1.0  )); //  9
  fQAhists->Add(new TH2F("hNormVtxDistBprongs"     , "Only Real B Decay Prongs;D meson #it{p}_{T};Distance to Primary Vertex (cm);counts (a.u.)"                                    , 24, 0.,24.  , 100,0.,1.0              )); // 10
  fQAhists->Add(new TH2F("hNormVtxDistXYBprongs"   , "Only Real B Decay Prongs;D meson #it{p}_{T};Distance to Primary Vertex in XY (cm);counts (a.u.)"                              , 24, 0.,24.  , 100,0.,1.0              )); // 11
  fQAhists->Add(new TH2F("hNormVtxDistBhypothesis" , "All Particles;D meson #it{p}_{T};Distance to Primary Vertex (cm);counts (a.u.)"                                               , 24, 0.,24.  , 100,0.,1.0              )); // 12
  fQAhists->Add(new TH2F("hNormVtxDistXYBhypthesis", "All Particles;D meson #it{p}_{T};Distance to Primary Vertex in XY (cm);counts (a.u.)"                                         , 24, 0.,24.  , 100,0.,1.0              )); // 13
  fQAhists->Add(new TH2F("hVtxPrecision"           , ";#it{d}_{z}^{Vertex} (cm);#it{d}_{xyz}^{Vertex} (cm);counts (a.u.)"                                                           ,100, 0., 0.05,1000,0.,0.05             )); // 14
  fQAhists->Add(new TH2F("hDecayLengthPrecision"   , ";D meson #it{p}_{T};Decay Length Residual (cm);counts (a.u.)"                                                                 , 24, 0.,24.  ,1000,0.,0.05             )); // 15
  fQAhists->Add(new TH2F("hDecayLengthXYPrecision" , ";D meson #it{p}_{T};XY decay Length Residual(cm);counts (a.u.)"                                                               , 24, 0.,24.  ,1000,0.,0.05             )); // 16
  fQAhists->ls();

  return;
}

void AliHFsubtractBFDcuts::FillGenStep(AliAODMCParticle* dzeropart,Double_t pt/*=-1.*/,Double_t weight/*=1.*/,TClonesArray* mcArray/*=0x0*/,AliAODMCHeader* mcHeader/*=0x0*/){
  fMCarray=mcArray;
  if (pt<0) {
    pt=dzeropart->Pt();
  }
  if (fIsMC && fMCarray) {
    fNprongs=0;
    fNprongsInAcc=0;
    fFoundElectron=kFALSE;
    fDecayChain=kFALSE; // TODO: use this value
    fMotherPt=pt;
    fLabCand=dzeropart->GetLabel();
    fPtCand=dzeropart->Pt();
    Double_t vtxDist[]   = { mcHeader->GetVtxX()-((AliAODMCParticle*)mcArray->At(dzeropart->GetDaughterFirst()))->Xv(),
                             mcHeader->GetVtxY()-((AliAODMCParticle*)mcArray->At(dzeropart->GetDaughterFirst()))->Yv(),
                             mcHeader->GetVtxZ()-((AliAODMCParticle*)mcArray->At(dzeropart->GetDaughterFirst()))->Zv() };
    Double_t decayLength   = TMath::Sqrt(vtxDist[0]*vtxDist[0]+vtxDist[1]*vtxDist[1]+vtxDist[2]*vtxDist[2]);
    Double_t decayLengthXY = TMath::Sqrt(vtxDist[0]*vtxDist[0]+vtxDist[1]*vtxDist[1]);
    if (!AnalyseDecay(fGenerateDecayList, kTRUE)) {
      AliDebug(3, "Error during the decay type determination!");
    }
    Double_t entry[] = {pt,(Double_t)fNprongsInAcc,fMotherPt,decayLength,decayLengthXY,(Double_t)fFoundElectron};
    fTHnGenStep->Fill(entry,weight);
    // y distribution of all particles
    for (Int_t i=0; i<fMCarray->GetEntriesFast(); ++i) {
      Double_t y=((AliAODMCParticle*)fMCarray->UncheckedAt(i))->Y();
      ((TH1F*)fQAhists->At(0))->Fill(y); // all particles
      if (IsStable(i)) ((TH1F*)fQAhists->At(1))->Fill(y); // all stable particles
    }
  }
  else {
    Double_t entry[] = {pt,0.,-1., 0., 0.};
    fTHnGenStep->Fill(entry,weight);
  }
  fMCarray=0x0;
  return;
}

void AliHFsubtractBFDcuts::FillSparses(AliAODRecoDecayHF2Prong* dzerocand,Int_t isSelected,Double_t pt,Double_t massD0,Double_t massD0bar,Double_t weight,TClonesArray* mcArray,AliAODEvent* aodEvent, AliAODMCHeader* mcHeader){
  fMCarray=mcArray;
  if (aodEvent) {
    fAODtracks=aodEvent->GetTracks();
    fBkG      =aodEvent->GetMagneticField();
    fPriVtx   =aodEvent->GetPrimaryVertex();
  }
  if(isSelected<=0||isSelected>3){
    Printf("isSelected = %d", isSelected);
    return;
  }
  if(massD0<0){
    dzerocand->InvMassD0(massD0,massD0bar);
    Printf("mass D0 =%f, massD0bar = %f ", massD0, massD0bar);
  }
  if(pt<0){
    pt=dzerocand->Pt();
  }
  Double_t normalDecayLengXY=dzerocand->NormalizedDecayLengthXY();
  Double_t cptangXY=dzerocand->CosPointingAngleXY();
  Double_t normalDecayLeng=dzerocand->NormalizedDecayLength();
  Double_t cptang=dzerocand->CosPointingAngle();
  Double_t decayLengXY=dzerocand->DecayLengthXY();
  Double_t decayLeng=dzerocand->DecayLength();
  // mass, pt, normLXY, cosPointXY, normL, cosPoint, LXY, L
  Double_t pointData[]={massD0,pt,normalDecayLengXY,cptangXY,normalDecayLeng,cptang,decayLengXY,decayLeng};

  if(isSelected==1||isSelected==3){
    fTHnData->Fill(pointData,weight);
  }
  if(isSelected==2||isSelected==3){
    pointData[0]=massD0bar;
    fTHnData->Fill(pointData,weight);
  }

  if(fIsMC && fMCarray) {
    fNprongs=0;
    fNprongsInAcc=0;
    fDecayChain=kFALSE; // TODO: use this value
    fFoundElectron=kFALSE;
    fMotherPt=pt;
    fD0Cand=dzerocand;
    if (fD0Cand) fD0CandParam=new AliNeutralTrackParam(fD0Cand);
    if (!GetCandidateLabel() || !AnalyseDecay(kFALSE, kFALSE)) {
      AliDebug(3, "Error during the decay type determination!");
    }
    // pt, normLXY, cosPointXY, #prongs, mother pt, normL, cosPoint, LXY, L
    Double_t pointMC[]={pt,normalDecayLengXY,cptangXY,(Double_t)fNprongsInAcc,fMotherPt,normalDecayLeng,cptang,decayLengXY,decayLeng,(Double_t)fFoundElectron};
    fTHnMC->Fill(pointMC, weight);

    if (mcHeader) {
      Double_t dist[] = { fPriVtx->GetX()-mcHeader->GetVtxX(),
                          fPriVtx->GetY()-mcHeader->GetVtxY(),
                          fPriVtx->GetZ()-mcHeader->GetVtxZ() };
      Int_t labDau0=((AliAODTrack*)fD0Cand->GetDaughter(0))->GetLabel();
      if (labDau0<0) {
        labDau0 *= -1.;
        AliDebug(1, "Negative Daughter label");
      }
      Double_t decLengMC3D[]   = { mcHeader->GetVtxX()-((AliAODMCParticle*)mcArray->At(labDau0))->Xv(),
                                   mcHeader->GetVtxY()-((AliAODMCParticle*)mcArray->At(labDau0))->Yv(),
                                   mcHeader->GetVtxZ()-((AliAODMCParticle*)mcArray->At(labDau0))->Zv() };
      Double_t decayLengMC   = TMath::Sqrt(decLengMC3D[0]*decLengMC3D[0]+decLengMC3D[1]*decLengMC3D[1]+decLengMC3D[2]*decLengMC3D[2]);
      Double_t decayLengXYMC = TMath::Sqrt(decLengMC3D[0]*decLengMC3D[0]+decLengMC3D[1]*decLengMC3D[1]);


      ((TH2F*)fQAhists->At(14))->Fill(TMath::Sqrt(dist[2]*dist[2]),
                                      TMath::Sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]));
      ((TH2F*)fQAhists->At(15))->Fill(pt, decayLengMC-decayLeng);
      ((TH2F*)fQAhists->At(16))->Fill(pt, decayLengXYMC-decayLengXY);
    }
  }
  fMCarray=0x0;
  fAODtracks=0x0;
  fBkG=0.;
  fD0Cand=0x0;
  return;
}

Bool_t AliHFsubtractBFDcuts::GetCandidateLabel(){
  if (fD0Cand) {
    Int_t labDau0=((AliAODTrack*)fD0Cand->GetDaughter(0))->GetLabel();
    if (labDau0<0) {
      labDau0 *= -1.;
      AliDebug(1, "Negative Daughter label");
    }
    AliAODMCParticle* firstDau=(AliAODMCParticle*)fMCarray->UncheckedAt(labDau0);
    fLabCand = firstDau->GetMother();
    return kTRUE;
  }
  else {
    AliDebug(3, "Could not obtain the label of the candidate");
    return kFALSE;
  }
  return kFALSE;
}

Bool_t AliHFsubtractBFDcuts::AnalyseDecay(Bool_t generateString, Bool_t mcOnly) {
  AliAODMCParticle* cand=(AliAODMCParticle*)fMCarray->UncheckedAt(fLabCand);
  fLabMother = cand->GetMother();
  AliAODMCParticle* mother = (AliAODMCParticle*)fMCarray->UncheckedAt(fLabMother);
  Int_t pdgMother = mother->GetPdgCode();
  if (pdgMother<0) pdgMother*=-1; // treat particles and anti-particles the same way
  if (pdgMother==4 || pdgMother==2212 || pdgMother==2112) { // prompt production
    fNprongs=1;
    fNprongsInAcc=1;
    return kTRUE;
  }
  if ((pdgMother%1000)/100==4 || (pdgMother%10000)/1000==4) {
    // chained decay of charmed hadrons, using recursion to resolve it
    fDecayChain=kTRUE;
    fLabCand=fLabMother;
    return AnalyseDecay(generateString, mcOnly);
  }
  if ((pdgMother%1000)/100!=5 && (pdgMother%10000)/1000!=5) {
    AliDebug(3, "Found strange decay, expected the mother to be a beauty hadron!");
    fNprongs=0;
    fNprongsInAcc=0;
    return kFALSE;
  }
  CountProngs(fLabMother, fLabCand, generateString, mcOnly); // count the prongs

  if (generateString) {
    // Store the decay
    std::sort(fDecayProngs.begin(), fDecayProngs.end());
    TString decayStr = "";
    for (ULong64_t i=0; i<fDecayProngs.size(); ++i) {
      decayStr = (i==0) ? Form("%d", fDecayProngs[i]) : Form("%s_%d", decayStr.Data(), fDecayProngs[i]);
    }
    //decayStr = Form("%s__%d_%d_%d", decayStr.Data(), fNprongs, fNprongsInAcc, fDecayProngs.size());
    TObjString* str = new TObjString(decayStr);
    if (!fDecayStrList->FindObject(str)) fDecayStrList->Add(str); // only allow unique entries
    fDecayProngs.clear();
  }

  if (mcOnly) {
    Double_t decayLengthB;   // Decay length B meson
    Double_t decayLengthBxy; // Decay length B meson (xy-plane)
    Double_t originB[3] = {0.,0.,0.};
    if(!mother->XvYvZv(originB)) AliDebug(3, "Couldn't determine MC origin of the beauty hadron");
    Double_t originD[3] = {0.,0.,0.};
    if(!cand->XvYvZv(originD)) AliDebug(3, "Couldn't determine MC origin of the charmed hadron");
    decayLengthBxy = TMath::Sqrt((originB[0]-originD[0])*(originB[0]-originD[0])+
                                 (originB[1]-originD[1])*(originB[1]-originD[1]));
    decayLengthB   = TMath::Sqrt(decayLengthBxy*decayLengthBxy+(originB[2]-originD[2])*(originB[2]-originD[2]));
    ((TH2F*)fQAhists->At(2))->Fill(fPtCand,decayLengthB);
    ((TH2F*)fQAhists->At(3))->Fill(fPtCand,decayLengthBxy);
  }
  else {
    for (Int_t iAODtrack=0; iAODtrack<fAODtracks->GetEntriesFast(); ++iAODtrack) {
      CheckBhypothesis(iAODtrack, kFALSE);
    }
  }

  fMotherPt=mother->Pt();
  return kTRUE;
}

void AliHFsubtractBFDcuts::CountProngs(Int_t labCurrMother, Int_t labCurrExcl,
                                       Bool_t generateString, Bool_t mcOnly) {
  for (Int_t iMCParticle=0; iMCParticle<fMCarray->GetEntriesFast(); ++iMCParticle) {
    if (iMCParticle!=labCurrExcl) {
      if (((AliAODMCParticle*)fMCarray->UncheckedAt(iMCParticle))->GetMother()==labCurrMother) {
        if (!fResolveResonances || IsStable(iMCParticle)) {
          if (generateString) fDecayProngs.push_back(((AliAODMCParticle*)fMCarray->UncheckedAt(iMCParticle))->GetPdgCode());
          ++fNprongs;
          if (!mcOnly) {
            for (Int_t iAODtrack=0; iAODtrack<fAODtracks->GetEntriesFast(); ++iAODtrack) {
              AliAODTrack* aodTrack=(AliAODTrack*)fAODtracks->UncheckedAt(iAODtrack);
              if (aodTrack->GetLabel()==iMCParticle) CheckBhypothesis(iAODtrack, kTRUE);
            }
          }
          if (!fCheckAcceptance || IsInAcceptance(iMCParticle)) {
            ++fNprongsInAcc;
          }
        }
        else CountProngs(iMCParticle, -1, generateString, mcOnly);
      }
    }
    else {
      ++fNprongs; // candidate is only counted as a single prong
      ++fNprongsInAcc;
      if (generateString) fDecayProngs.push_back(((AliAODMCParticle*)fMCarray->UncheckedAt(iMCParticle))->GetPdgCode());
    }
  }
}

Bool_t AliHFsubtractBFDcuts::IsStable(Int_t labProng) {
  const Int_t stablePartPdgs[] = { 11, 13, 211, 321, 2212, 12, 14, 22, 111, 130 };
  const Int_t nStablePartPdgs  = sizeof(stablePartPdgs)/sizeof(Int_t);
  AliAODMCParticle* prong = (AliAODMCParticle*)fMCarray->UncheckedAt(labProng);
  Int_t pdgProng = prong->GetPdgCode();
  if (pdgProng<0) pdgProng*=-1; // treat particles and anti-particles the same way
  for (Int_t iPdg=0; iPdg<nStablePartPdgs; ++iPdg) {
    if (pdgProng == 11) fFoundElectron = kTRUE;
    if (stablePartPdgs[iPdg] == pdgProng) return kTRUE;
  }
  return kFALSE;
}

Bool_t AliHFsubtractBFDcuts::IsInAcceptance(Int_t labProng) const {
  AliDebug(1, "AliHFsubtractBFDcuts::IsInAcceptance(...) hasn't been implemented yet, prong");
  Double_t eta    = ((AliAODMCParticle*)fMCarray->UncheckedAt(labProng))->Eta();
  Double_t pt     = ((AliAODMCParticle*)fMCarray->UncheckedAt(labProng))->Pt();
  Short_t charge = ((AliAODMCParticle*)fMCarray->UncheckedAt(labProng))->Charge();
  if ((pt>0.15) && (eta>-0.9) && (eta<0.9) && (charge!=0)) return kTRUE;
  return kFALSE;
}

Bool_t AliHFsubtractBFDcuts::CheckBhypothesis(Int_t iAODtrack, Bool_t Bprong) {
  AliExternalTrackParam* t = new AliExternalTrackParam();
  t->CopyFromVTrack((AliVTrack*)fAODtracks->UncheckedAt(iAODtrack));

  TObjArray* tracks = new TObjArray(2);
  tracks->AddAt(t,           0);
  tracks->AddAt(fD0CandParam,1);

  AliAODVertex* bVtx = RecBvtx(tracks);
  if(!bVtx) {
    AliDebug(3, "Couldn't reconstruct B meson vertex!");
    delete t; t=0x0;
    return kFALSE;
  }

  Double_t vtxDist      = (fPriVtx) ? bVtx->DistanceToVertex(fPriVtx)        : -1.0  ;
  Double_t vtxDistXY    = (fPriVtx) ? bVtx->DistanceXYToVertex(fPriVtx)      : -2.e10;
  Double_t vtxDistErr   = (fPriVtx) ? bVtx->ErrorDistanceToVertex(fPriVtx)   : -1.0  ;
  Double_t vtxDistErrXY = (fPriVtx) ? bVtx->ErrorDistanceXYToVertex(fPriVtx) : -2.e10;

  const Double_t maxD = 1.;

  // Propagate candidates to secondary vertex
  Double_t dz[2],cov[3];
  t->PropagateToDCA(bVtx,fBkG,maxD,dz,cov);
  fD0CandParam->PropagateToDCA(bVtx,fBkG,maxD,dz,cov);

  // Impact parameters
  t->PropagateToDCA(fPriVtx,fBkG,maxD,dz,cov);
  Double_t d0Track=dz[0];
  Double_t d0TrackErr=TMath::Sqrt(cov[0]);
  fD0CandParam->PropagateToDCA(fPriVtx,fBkG,maxD,dz,cov);
  Double_t d0D0Cand=dz[0];
  Double_t d0D0CandErr=TMath::Sqrt(cov[0]);

  // distance of closest approach of the D0 and the track
  Double_t xDCAtrack, xDCAD0;
  Double_t dcaB=t->GetDCA(fD0CandParam,fBkG,xDCAtrack,xDCAD0);
  if (!Bprong) {
    ((TH3F*)fQAhists->At(4))->Fill(fPtCand,dcaB,d0D0Cand*d0Track);
    ((TH3F*)fQAhists->At(6))->Fill(fPtCand,vtxDist,vtxDistErr);
    ((TH3F*)fQAhists->At(7))->Fill(fPtCand,vtxDistXY,vtxDistErrXY);
    ((TH2F*)fQAhists->At(10))->Fill(fPtCand,vtxDist/vtxDistErr);
    ((TH2F*)fQAhists->At(11))->Fill(fPtCand,vtxDistXY/vtxDistErrXY);
  }
  else {
    ((TH3F*)fQAhists->At(5))->Fill(fPtCand,dcaB,d0D0Cand*d0Track);
    ((TH3F*)fQAhists->At(8))->Fill(fPtCand,vtxDist,vtxDistErr);
    ((TH3F*)fQAhists->At(9))->Fill(fPtCand,vtxDistXY,vtxDistErrXY);
    ((TH2F*)fQAhists->At(12))->Fill(fPtCand,vtxDist/vtxDistErr);
    ((TH2F*)fQAhists->At(13))->Fill(fPtCand,vtxDistXY/vtxDistErrXY);
  }

  delete tracks; tracks=0x0;
  delete bVtx  ; bVtx  =0x0;
  delete t     ; t     =0x0;
  return kTRUE;
}

AliAODVertex* AliHFsubtractBFDcuts::RecBvtx(TObjArray* tracks) const {
  AliESDVertex* vtxESD=0x0;
  AliVertexerTracks* vertexer = new AliVertexerTracks(fBkG);
  vertexer->SetVtxStart((AliESDVertex*)fPriVtx);
  vtxESD = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(tracks);
  delete vertexer; vertexer=0x0;
  if(!vtxESD || vtxESD->GetNContributors()!=tracks->GetEntriesFast()) {
    AliDebug(3, "Couldn't reconstruct B Meson vertex");
    delete vtxESD; vtxESD=0x0;
    return 0x0;
  }
  Double_t rVtxSq=vtxESD->GetX()*vtxESD->GetX()+vtxESD->GetY()*vtxESD->GetY();
  if(rVtxSq>8.){
    // vertex outside beam pipe, reject candidate to avoid propagation through material
    delete vtxESD; vtxESD=0x0;
    return 0x0;
  }
  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  vtxESD->GetXYZ(pos);
  vtxESD->GetCovMatrix(cov);
  chi2perNDF=vtxESD->GetChi2toNDF();
  delete vtxESD; vtxESD=0x0;
  return new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,0);
}
