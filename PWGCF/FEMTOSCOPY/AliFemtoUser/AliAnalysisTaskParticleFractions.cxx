#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TObjArray.h"
#include "TString.h"
#include "TParticle.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"

#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODpidUtil.h"
#include "AliAnalysisUtils.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliInputEventHandler.h"

#include "AliAnalysisTaskParticleFractions.h"

ClassImp(AliAnalysisTaskParticleFractions)

using std::cout;
using std::endl;

//_______________________________________________________

AliAnalysisTaskParticleFractions::AliAnalysisTaskParticleFractions(const Char_t *partName) :
  AliAnalysisTaskSE(partName), centrality(0), fpidResponse(0), fHistoList(0),  fHistEv(0)
{
  for(Int_t i = 0; i < multbins*parttypes; i++) {
    fParticleOriginMC[i] = NULL;
    fParticleOriginRec[i] = NULL;
  }

  /* create output */
  fHistoList = new TList();
  fHistoList->SetOwner(kTRUE);

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//_______________________________________________________

AliAnalysisTaskParticleFractions::~AliAnalysisTaskParticleFractions()
{
}

//_______________________________________________________

void AliAnalysisTaskParticleFractions::UserCreateOutputObjects()
{
  TString hname1, hname2;
  TString htitle1, htitle2;
  TString hname1M, hname2M, hname;
  TString htitle1M, htitle2M, htitle;
  TString parttypename = "None";

  for (Int_t j = 0; j < parttypes; j++)  {
    if (j==0) parttypename="All";
    else if (j==1) parttypename="Pion";
    else if (j==2) parttypename="Kaon";
    else if (j==3) parttypename="Proton";
    else if (j==4) parttypename="Lambda";

    for (Int_t i = 0; i < multbins; i++)  {
      hname1  = "hParticleOriginMC"; hname1+=i; hname1+=parttypename;
      htitle1 = "Particle Origin MC"; htitle1+=i; htitle1+=parttypename;
      fParticleOriginMC[i*parttypes+j] = new TH1F(hname1.Data(),htitle1.Data(),6000, 0, 6000.);;

      hname2  = "hParticleOriginREC"; hname2+=i; hname2+=parttypename;
      htitle2 = "Particle Origin Rec"; htitle2+=i; htitle2+=parttypename;
      fParticleOriginRec[i*parttypes+j] = new TH1F(hname2.Data(),htitle2.Data(),6000, 0, 6000.);;
    }
  }

  fHistEv = new TH1F("fHistEv", "Multiplicity", 100, 0, 100);
  fHistoList->Add(fHistEv);

  for (Int_t i = 0; i < multbins*parttypes; i++) {
    fHistoList->Add(fParticleOriginMC[i]);
    fHistoList->Add(fParticleOriginRec[i]);
  }

  //********** PID ****************
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fpidResponse = inputHandler->GetPIDResponse();
  cout<<"*******"<< fpidResponse<<endl;
  // ************************

  PostData(1, fHistoList);
}

//_____________________________________________________________________

bool AliAnalysisTaskParticleFractions::IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
{
  if (mom > 0.5) {
    if (TMath::Hypot( nsigmaTOFPi, nsigmaTPCPi ) < 3)
      return true;
	}
  else {
    if (TMath::Abs(nsigmaTPCPi) < 3)
      return true;
  }
  return false;
}

bool AliAnalysisTaskParticleFractions::IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
{
  if (mom > 0.5) {
    if (TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < 3)
      return true;
	}
  else {
    if (TMath::Abs(nsigmaTPCK) < 3)
      return true;
  }
  return false;
}

bool AliAnalysisTaskParticleFractions::IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP)
{
  if (mom > 0.5) {
    if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < 3)
      return true;
	}
  else {
    if (TMath::Abs(nsigmaTPCP) < 3)
      return true;
  }
  return false;
}


bool AliAnalysisTaskParticleFractions::IsElectron(float nsigmaTPCE, float nsigmaTPCPi,float nsigmaTPCK, float nsigmaTPCP)
{
  if (TMath::Abs(nsigmaTPCE)<3 && TMath::Abs(nsigmaTPCPi)>3 && TMath::Abs(nsigmaTPCK)>3 && TMath::Abs(nsigmaTPCP)>3)
    return true;
  else
    return false;
}

//_______________________________________________________

void AliAnalysisTaskParticleFractions::UserExec(Option_t *)
{
  /***Get Event****/
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!aodEvent) return;
  AliAODHeader *fAODheader = dynamic_cast<AliAODHeader*>(aodEvent->GetHeader());
  if(!fAODheader) AliFatal("Not a standard AOD");
  AliCentrality* alicent= aodEvent->GetCentrality(); //in PbPb and pPb
  Double_t mult = alicent->GetCentralityPercentile("V0M"); //in PbPb

  // EVENT SELECTION ********************
  //collision candidate
  if (!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB)) return;

  //****** Multiplicity selection *********
  Int_t fcent = 0;
  // Int_t fcent = -999;
  // if (mult >= 0 && mult <=20)  fcent = 0;
  // else if (mult >= 20 && mult <=39) fcent = 1;
  // else if (mult >= 40 && mult <=59) fcent = 2;
  // else if (mult >= 60 && mult <=90) fcent = 3;
  // else  return;

  const AliAODVertex* vertex =(AliAODVertex*) aodEvent->GetPrimaryVertex();
  if (!vertex || vertex->GetNContributors()<=0) return;
  vertex->GetPosition(fV1);

  AliAnalysisUtils *anaUtil=new AliAnalysisUtils();

  Bool_t fpA2013 = kFALSE;
  Bool_t fMVPlp = kFALSE;
  Bool_t fisPileUp = kTRUE;
  Int_t fMinPlpContribMV = 0;
  Int_t fMinPlpContribSPD = 3;

  if(fpA2013)
    if (anaUtil->IsVertexSelected2013pA(aodEvent)==kFALSE) return;

  if (fMVPlp) anaUtil->SetUseMVPlpSelection(kTRUE);
  else anaUtil->SetUseMVPlpSelection(kFALSE);

  if (fMinPlpContribMV) anaUtil->SetMinPlpContribMV(fMinPlpContribMV);
  if (fMinPlpContribSPD) anaUtil->SetMinPlpContribSPD(fMinPlpContribSPD);

  if (fisPileUp)
    if (anaUtil->IsPileUpEvent(aodEvent)) return;

  delete anaUtil;
  Float_t zvtx = vertex->GetZ();
  if (TMath::Abs(zvtx) > 10) return;

  fHistEv->Fill(mult);
  cout << mult<< endl;
  //**** getting MC array ******
  TClonesArray  *arrayMC;
  arrayMC = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));

//copying pid information for FB 128
  int labels[20000];
  for (int il=0; il<20000; il++) labels[il] = -1;

  cout << aodEvent->GetNumberOfTracks() << endl;
  // looking for global tracks and saving their numbers to copy from them PID information to TPC-only tracks in the main loop over tracks
  for (int i=0;i<aodEvent->GetNumberOfTracks();i++) {
    const AliAODTrack *aodtrack=dynamic_cast<const AliAODTrack*>(aodEvent->GetTrack(i));
    if (!aodtrack) AliFatal("Not a standard AOD");
    if (!aodtrack->TestFilterBit(128)) {
      if (aodtrack->GetID() < 0) continue;
      labels[aodtrack->GetID()] = i;
    }
  }

  //RECONSTRUCTED TRACKS
  TObjArray recoParticleArray[parttypes];

  //loop over AOD tracks
  for (Int_t iTracks = 0; iTracks < aodEvent->GetNumberOfTracks(); iTracks++) {
    //get track
    AliAODTrack *track = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTracks));
    if (!track) AliFatal("Not a standard AOD");
    if (!track) continue;

    UInt_t filterBit = 128;
    if(!track->TestFilterBit(filterBit))continue;

    if (track->Eta() < -0.8 || track->Eta() > 0.8) continue;
    if (track->Pt() < 0.2 || track->Pt() > 20) continue;

    // //DCA
    // Double_t DCAXY;
    // Double_t DCAZ;
    // //  if(filterBit==(1 << (7))){
    // DCAXY = -TMath::Abs(track->DCA());
    // DCAZ = -TMath::Abs(track->ZAtDCA());

    // if(!(DCAXY==-999 || DCAZ==-999)){
    //   //if(TMath::Abs(DCAXY) > 0.0182 + 0.035*TMath::Power(track->Pt(), -1.01)) continue; //XY, Pt dep
    //   //no DCA cut
    //   //if(TMath::Abs(DCAXY) > 1000.0) {continue;} //XY
    //   //if(TMath::Abs(DCAZ) > 1000.0) {continue;} //Z
    // }
    // else {
    //   // code from Michael and Prabhat from AliAnalysisTaskDptDptCorrelations
    //   // const AliAODVertex* vertex = (AliAODVertex*) aodEvent->GetPrimaryVertex(); (already defined above)
    //   float vertexX  = -999.;
    //   float vertexY  = -999.;
    //   float vertexZ  = -999.;

    //   if(vertex) {
    //     Double32_t fCov[6];
    //     vertex->GetCovarianceMatrix(fCov);
    //     if(vertex->GetNContributors() > 0) {
    //       if(fCov[5] != 0) {
    //         vertexX = vertex->GetX();
    //         vertexY = vertex->GetY();
    //         vertexZ = vertex->GetZ();
    //       }
    //     }
    //   }

    //   Double_t pos[3];
    //   track->GetXYZ(pos);

    //   Double_t DCAX = pos[0] - vertexX;
    //   Double_t DCAY = pos[1] - vertexY;
    //   DCAZ = pos[2] - vertexZ;
    //   DCAXY = TMath::Sqrt((DCAX*DCAX) + (DCAY*DCAY));
    // }

    AliAODTrack* aodtrackpid;

    //for FB 128 - tpc only tracks
    if(filterBit==(1 << (7))) {
      aodtrackpid = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(labels[-1-aodEvent->GetTrack(iTracks)->GetID()]));
      if(!aodtrackpid) AliFatal("Not a standard AOD");
    }
    else
      aodtrackpid = track;

    AliMCEvent* mcEvent = MCEvent();
    if (!mcEvent) {
      AliError("ERROR: Could not retrieve MC event");
      return;
    }

    Int_t nMCParticles = mcEvent->GetNumberOfTracks();
    Int_t labelp = TMath::Abs(track->GetLabel());
    if (labelp > nMCParticles) continue;

    AliAODMCParticle *AODmcTrack = (AliAODMCParticle*) mcEvent->GetTrack(labelp);

    //Electron rejection
    double nSigmaTPCPi = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kPion);
    double nSigmaTPCK = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kKaon);
    double nSigmaTPCP = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kProton);
    double nSigmaTPCe = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kElectron);
    if(IsElectron(nSigmaTPCe,nSigmaTPCPi,nSigmaTPCK,nSigmaTPCP)) continue;

    float px=track->Px(); float py=track->Py();  float ph=atan2(py,px); //track->Phi()
    float tPt = track->Pt();

    //PID monitors
    double nSigmaTOFPi = fpidResponse->NumberOfSigmasTOF(aodtrackpid,AliPID::kPion);
    double nSigmaTOFK = fpidResponse->NumberOfSigmasTOF(aodtrackpid,AliPID::kKaon);
    double nSigmaTOFP = fpidResponse->NumberOfSigmasTOF(aodtrackpid,AliPID::kProton);

    float tdEdx = aodtrackpid->GetTPCsignal();
    float tTofSig = aodtrackpid->GetTOFsignal();
    double pidTime[5]; aodtrackpid->GetIntegratedTimes(pidTime);

    if (!arrayMC) continue;
    //get coresponding MC particle
    Int_t label = TMath::Abs(track->GetLabel());
    AliAODMCParticle *MCtrk = (AliAODMCParticle*)arrayMC->At(label);

    //********* PID - pions ********
    if (IsPionNSigma(track->P(),nSigmaTPCPi, nSigmaTOFPi)){
      //if (!MCtrk)
      continue;
      //recoParticleArray[1].Add(MCtrk);
    }
    //********* PID - kaons ********
    if (IsKaonNSigma(track->P(),nSigmaTPCK, nSigmaTOFK)){
      //if (!MCtrk)
      continue;
      //recoParticleArray[2].Add(MCtrk);
    }
    //********* PID - protons ********
    if (IsProtonNSigma(track->P(),nSigmaTPCP, nSigmaTOFP)){
      if (!MCtrk) continue;
      recoParticleArray[3].Add(MCtrk);
    }

    double LambdaMass = 1.115683;
    // loop over V0s
    for (Int_t i = 0; i < aodEvent->GetNumberOfV0s(); i++) {
      AliAODv0* aodv0 = aodEvent->GetV0(i);
      if (!aodv0) continue;
      if (aodv0->GetNDaughters()>2) continue;
      if (aodv0->GetNProngs()>2) continue;
      if (aodv0->GetCharge()!=0) continue;
      if (aodv0->ChargeProng(0)==aodv0->ChargeProng(1)) continue;
      if (aodv0->CosPointingAngle(fV1)<0.9993) continue;
      if (aodv0->DecayLength(fV1)>60) continue;

      AliAODTrack* daughterTrackPos = (AliAODTrack*)aodv0->GetDaughter(0); //getting positive daughter track
      AliAODTrack* daughterTrackNeg = (AliAODTrack*)aodv0->GetDaughter(1); //getting negative daughter track
      if (!daughterTrackPos) continue; //daughter tracks must exist
      if (!daughterTrackNeg) continue;
      if (daughterTrackNeg->Charge() == daughterTrackPos->Charge() ) continue; //and have different charge

      if (aodv0->Eta() < -0.8 || aodv0->Eta() > 0.8) continue;
      if (daughterTrackPos->Eta() < -0.8 || daughterTrackPos->Eta() > 0.8) continue;
      if (daughterTrackNeg->Eta() < -0.8 || daughterTrackNeg->Eta() > 0.8) continue;
      if (aodv0->Pt() < 0.5 || aodv0->Pt() > 5) continue;

      if (daughterTrackPos->Pt() < 0.5 || daughterTrackPos->Pt() > 4) continue;
      if (daughterTrackNeg->Pt() < 0.16 || daughterTrackNeg->Pt() > 4) continue;
      if (daughterTrackPos->GetTPCNcls() < 80) continue;
      if (daughterTrackNeg->GetTPCNcls() < 80) continue;
      if (daughterTrackPos->Chi2perNDF() > 4) continue;
      if (daughterTrackNeg->Chi2perNDF() > 4) continue;

      if (!(daughterTrackPos->GetStatus() & AliESDtrack::kTPCrefit)) continue;
      if (!(daughterTrackNeg->GetStatus() & AliESDtrack::kTPCrefit)) continue;

      if (TMath::Abs(aodv0->DcaV0Daughters())>0.4) continue;
      if (TMath::Abs(aodv0->DcaPosToPrimVertex())<0.1 || TMath::Abs(aodv0->DcaNegToPrimVertex())<0.3) continue;
      if (TMath::Abs(aodv0->DcaV0ToPrimVertex())>1.0) continue;

      double nSigmaTPCP = fpidResponse->NumberOfSigmasTPC(daughterTrackPos,AliPID::kProton);
      double nSigmaTOFP = fpidResponse->NumberOfSigmasTOF(daughterTrackPos,AliPID::kProton);
      double nSigmaTPCPi = fpidResponse->NumberOfSigmasTPC(daughterTrackNeg,AliPID::kPion);
      double nSigmaTOFPi = fpidResponse->NumberOfSigmasTOF(daughterTrackNeg,AliPID::kPion);

      if (!IsProtonNSigma(daughterTrackPos->Pt(), nSigmaTPCP, nSigmaTOFP)) continue;
      if (!IsPionNSigma(daughterTrackNeg->Pt(), nSigmaTPCPi, nSigmaTOFPi)) continue;
      if (aodv0->MassLambda()<LambdaMass-0.0041 || aodv0->MassLambda()>LambdaMass+0.0041) continue;

      if (daughterTrackPos->GetLabel() > 0 && daughterTrackNeg->GetLabel() > 0 ) {

        AliAODMCParticle* mcParticlePos = (AliAODMCParticle*)arrayMC->At(daughterTrackPos->GetLabel());
        AliAODMCParticle* mcParticleNeg = (AliAODMCParticle*)arrayMC->At(daughterTrackNeg->GetLabel());

        if ((mcParticlePos!=NULL) && (mcParticleNeg!=NULL)){
          int motherOfPosID = mcParticlePos->GetMother();
          int motherOfNegID = mcParticleNeg->GetMother();
          if ((motherOfPosID > -1) && (motherOfPosID == motherOfNegID)){
            AliAODMCParticle *MCv0 = (AliAODMCParticle*)arrayMC->At(motherOfPosID); //our V0 particle
            if (!MCv0) continue;
            recoParticleArray[4].Add(MCv0);
          }
        }
      }
    }

    if (!MCtrk) continue;
    if(MCtrk->Charge()==0){cout<<"!!!"<<endl; continue;}
    recoParticleArray[0].Add(MCtrk);
  }

  // MONTECARLO PARTICLES
  if(!arrayMC){
    AliError("Array of MC particles not found");
    return;
  }
  cout << arrayMC->GetEntries() << endl;
  // loop over MC stack
  for (Int_t ipart = 0; ipart < arrayMC->GetEntries(); ipart++) {
    AliAODMCParticle *MCtrk = (AliAODMCParticle*)arrayMC->At(ipart);
    if (!MCtrk) continue;
    Int_t PDGcode = MCtrk->GetPdgCode();
    if (! (PDGcode==3122 || PDGcode==2212)) continue;
    Int_t MCTrkMotherID = MCtrk->GetMother();
    Int_t PDGcodeMother = 0;
    if (MCTrkMotherID > -1) { // particle has a mother
      AliAODMCParticle* MCTrkMother = (AliAODMCParticle*)arrayMC->At(MCTrkMotherID);
      PDGcodeMother = MCTrkMother->GetPdgCode();
    }

    if (MCtrk->Eta() < -0.8 || MCtrk->Eta() > 0.8) continue;
    if (MCtrk->Pt() < 0.2 || MCtrk->Pt() > 5) continue;

    // // filling histograms for MC truth particles
    // if(PDGcode==211){
    //   continue;
    //   //fParticleOriginMC[fcent*parttypes+1]->Fill(PDGcodeMother);
    // }
    // else if(PDGcode==321) {
    //   continue;
    //   //fParticleOriginMC[fcent*parttypes+2]->Fill(PDGcodeMother);
    // }
    // else
    if(PDGcode==2212) {
      fParticleOriginMC[fcent*parttypes+3]->Fill(PDGcodeMother);
    }
    else if(PDGcode==3122) {
      fParticleOriginMC[fcent*parttypes+4]->Fill(PDGcodeMother);
    }

    // //Filling data from MC truth particles only for particles that were reconstruced
    // if (recoParticleArray[1].Contains(MCtrk)){ //Pions
    //   if(PDGcode==211)
    //   fParticleOriginRec[fcent*parttypes+1]->Fill(PDGcodeMother);
    // }
    // if (recoParticleArray[2].Contains(MCtrk)){ //Kaons
    //   if(PDGcode==321)
    //   fParticleOriginRec[fcent*parttypes+2]->Fill(PDGcodeMother);
    // }
    if (recoParticleArray[3].Contains(MCtrk)){ //Protons
      if(PDGcode==2212){
        fParticleOriginRec[fcent*parttypes+3]->Fill(PDGcodeMother);
      }
    }
    if (recoParticleArray[4].Contains(MCtrk)){ //Lambdas
      if(PDGcode==3122){
        fParticleOriginRec[fcent*parttypes+4]->Fill(PDGcodeMother);
      }
    }
  }
  PostData(1, fHistoList);
}
