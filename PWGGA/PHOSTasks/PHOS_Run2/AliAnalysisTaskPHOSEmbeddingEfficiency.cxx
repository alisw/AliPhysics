#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "THnSparse.h"
#include "THashList.h"
#include "TMath.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliLog.h"

#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"

#include "AliESDtrackCuts.h"
#include "AliESDHeader.h"
#include "AliESDEvent.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloCluster.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"

#include "AliVEvent.h"
#include "AliVHeader.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"

#include "AliMultSelection.h"
#include "AliEventplane.h"
#include "AliQnCorrectionsManager.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliQnCorrectionsQnVector.h"

#include "AliOADBContainer.h"

#include "AliPID.h"
#include "AliPIDResponse.h"

#include "AliAODMCHeader.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODInputHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"


#include "AliPHOSEventCuts.h"
#include "AliPHOSClusterCuts.h"
#include "AliPHOSJetJetMC.h"
#include "AliPHOSTriggerHelper.h"
#include "AliAnalysisTaskPHOSPi0EtaToGammaGamma.h"

#include "AliAnalysisTaskPHOSEmbeddingEfficiency.h"

// Author: Daiki Sekihata (Hiroshima University)

using namespace std;

ClassImp(AliAnalysisTaskPHOSEmbeddingEfficiency)

//________________________________________________________________________
AliAnalysisTaskPHOSEmbeddingEfficiency::AliAnalysisTaskPHOSEmbeddingEfficiency(const char *name):
  AliAnalysisTaskPHOSPi0EtaToGammaGamma(name),
  fParticleName(""),
  fMCArray(0x0),
  fWeightCen0005(0x0)
{
  // Constructor


  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, THashList::Class());

}
//________________________________________________________________________
AliAnalysisTaskPHOSEmbeddingEfficiency::~AliAnalysisTaskPHOSEmbeddingEfficiency()
{
  delete fWeightCen0005;
  fWeightCen0005 = 0x0;


}
//________________________________________________________________________
void AliAnalysisTaskPHOSEmbeddingEfficiency::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  fWeightCen0005 = new TF1("fWeightCen0005" ,"64.41*TMath::Power(x,-(5.88 + -92.9/(TMath::Power(x,4.12) + 54.1)))" ,0.,100.);
  fWeightCen0005->SetNpx(1000);

  AliAnalysisTaskPHOSPi0EtaToGammaGamma::UserCreateOutputObjects();

  const Int_t NpTgg = 101;
  Double_t pTgg[NpTgg]={};
  for(Int_t i=0;i<50;i++)     pTgg[i] = 0.1 * i;            //every 0.1 GeV/c, up to 5 GeV/c
  for(Int_t i=50;i<60;i++)    pTgg[i] = 0.5 * (i-50) + 5.0; //every 0.5 GeV/c, up to 10 GeV/c
  for(Int_t i=60;i<NpTgg;i++) pTgg[i] = 1.0 * (i-60) + 10.0;//every 1.0 GeV/c, up to 50 GeV/c

  const Int_t Npar = 3;
  const TString parname[Npar] = {"Pi0","Eta","Gamma"};
  for(Int_t ipar=0;ipar<Npar;ipar++){
    TH1F *h1Pt = new TH1F(Form("hGenEmbedded%sPt",parname[ipar].Data()        ),Form("generated %s pT",parname[ipar].Data()        ),NpTgg-1,pTgg);
    h1Pt->Sumw2();
    fOutputContainer->Add(h1Pt);

    TH2F *h2EtaPhi = new TH2F(Form("hGenEmbedded%sEtaPhi",parname[ipar].Data()),Form("generated %s eta vs phi",parname[ipar].Data()),200,-1,1,60,0,TMath::TwoPi());
    h2EtaPhi->Sumw2();
    fOutputContainer->Add(h2EtaPhi);

    TH2F *h2EtaPt = new TH2F(Form("hGenEmbedded%sEtaPt",parname[ipar].Data()  ),Form("generated %s eta vs pT",parname[ipar].Data() ),200,-1,1,NpTgg-1,pTgg);
    h2EtaPt->Sumw2();
    fOutputContainer->Add(h2EtaPt);

  }//end of particle loop


  PostData(1,fOutputContainer);

}
//________________________________________________________________________
void AliAnalysisTaskPHOSEmbeddingEfficiency::UserExec(Option_t *option) 
{
  // Main loop
  // Called for each event

  if(!fPHOSEventCuts){
    AliError("fPHOSEventCuts is not set! return");
    return;
  }
  if(!fPHOSClusterCuts){
    AliError("fPHOSClusterCuts is not set! return");
    return;
  }

  fEvent = dynamic_cast<AliVEvent*>(InputEvent());
  if(!fEvent){
    AliError("event is not available.");
    return;
  }

  fESDEvent = dynamic_cast<AliESDEvent*>(fEvent);
  fAODEvent = dynamic_cast<AliAODEvent*>(fEvent);

  FillHistogramTH1(fOutputContainer,"hEventSummary",1);//all

  const AliVVertex *vVertex = fEvent->GetPrimaryVertex();
  fVertex[0] = vVertex->GetX();
  fVertex[1] = vVertex->GetY();
  fVertex[2] = vVertex->GetZ();
  FillHistogramTH1(fOutputContainer,"hVertexZ" ,fVertex[2]);
  fZvtx = (Int_t)((fVertex[2]+10.)/2.);//it should be 0-9.
  if(fZvtx < 0) fZvtx = 0;//protection to avoid fZvtx = -1.
  if(fZvtx > 9) fZvtx = 9;//protection to avoid fZvtx = 10.

  Int_t Ncontributor  = vVertex->GetNContributors();

  //centrality estimation

  Float_t fCentralityV0M = -1.;
  Float_t fCentralityCL0 = -1.;
  Float_t fCentralityCL1 = -1.;
  Float_t fCentralityV0A = -1.;
  Float_t fCentralityV0C = -1.;
  Float_t fCentralityZNA = -1.;
  Float_t fCentralityZNC = -1.;

  if(fEstimator.Contains("V0") || fEstimator.Contains("ZN") || fEstimator.Contains("CL")){

    //Get Centrality
    fMultSelection = (AliMultSelection*)fEvent->FindListObject("MultSelection");
    if(!fMultSelection){
      //If you get this warning (and fCentralityV0M 300) please check that the AliMultSelectionTask actually ran (before your task)
      AliWarning("AliMultSelection object not found!");
      return;
    }
    else{
      fCentralityV0M  = fMultSelection->GetMultiplicityPercentile("V0M");
      fCentralityCL0  = fMultSelection->GetMultiplicityPercentile("CL0");
      fCentralityCL1  = fMultSelection->GetMultiplicityPercentile("CL1");
      fCentralityV0A  = fMultSelection->GetMultiplicityPercentile("V0A");
      fCentralityV0C  = fMultSelection->GetMultiplicityPercentile("V0C");
      fCentralityZNA  = fMultSelection->GetMultiplicityPercentile("ZNA");
      fCentralityZNC  = fMultSelection->GetMultiplicityPercentile("ZNC");
      fCentralityMain = fMultSelection->GetMultiplicityPercentile(fEstimator);
    }

  }
  else if(fEstimator.Contains("HybridTrack")){
    //hybrid track multiplicity 
    //done manually in this task.
    Int_t NHybrid = 0; 
    const Int_t trackMult = fEvent->GetNumberOfTracks();
    if(fESDEvent){
      for(Int_t itrack=0;itrack<trackMult;itrack++){
        AliESDtrack *esdtrack = (AliESDtrack*)fEvent->GetTrack(itrack);
        if(TMath::Abs(esdtrack->Eta()) > 0.8) continue;

        if(fESDtrackCutsGlobal           ->AcceptTrack(esdtrack)//select global track
        || fESDtrackCutsGlobalConstrained->AcceptTrack(esdtrack)){//select complementary track
          NHybrid++;
        }

      }//end of track loop
    }//end of ESD
    else if(fAODEvent){
      for(Int_t itrack=0;itrack<trackMult;itrack++){
        AliAODTrack *aodtrack = (AliAODTrack*)fEvent->GetTrack(itrack);
        if(TMath::Abs(aodtrack->Eta()) > 0.8) continue;

        if(aodtrack->IsHybridGlobalConstrainedGlobal()){//hybrid track
          NHybrid++;
        }

      }//end of track loop
    }//end of AOD
    fCentralityMain = (Float_t)NHybrid;
  }
  else if(fEstimator.Contains("SPDTracklet")){
    //hybrid track multiplicity 
    fCentralityMain = (Float_t)(fEvent->GetMultiplicity()->GetNumberOfTracklets());
  }
  else{
    AliInfo(Form("%s is not supported. return",fEstimator.Data()));
    return;
  }

  if(fCentralityMain < fCentralityMin || fCentralityMax < fCentralityMain){
    AliInfo(Form("Reject this event because centrality %s %f %% is out of the configuration of this task.", fEstimator.Data(),fCentralityMain));
    return;
  }

  FillHistogramTH2(fOutputContainer,Form("hCentrality%svsNContributor",fEstimator.Data()),fCentralityMain,Ncontributor);

  UInt_t fSelectMask = fInputHandler->IsEventSelected();
  Bool_t isINT7selected = fSelectMask & AliVEvent::kINT7;

  if(!fIsPHOSTriggerAnalysis && !isINT7selected){
    AliInfo("INT7 Event is rejected by IsEventSelected()");
    return;
  }

  //event selection
  if(!(fPHOSEventCuts->AcceptEvent(fEvent))){
    AliInfo("event is rejected.");
    return;
  }

  //<- end of event selection
  //-> start physics analysis

  fPHOSClusterArray = (TClonesArray*)fEvent->FindListObject(Form("PHOSEmbeddedDiffClusterArray_%s",fParticleName.Data()));

  if(!fPHOSClusterArray){
    AliWarning("fPHOSClusterArray object not found!");
    return;
  }

  AliInfo(Form("Particle : %s , Ncluster produced by embedding = %d.",fParticleName.Data(),fPHOSClusterArray->GetEntriesFast()));

  FillHistogramTH2(fOutputContainer,"hCentralityV0MvsCL0",fCentralityV0M,fCentralityCL0);
  FillHistogramTH2(fOutputContainer,"hCentralityV0MvsCL1",fCentralityV0M,fCentralityCL1);
  FillHistogramTH2(fOutputContainer,"hCentralityCL0vsCL1",fCentralityCL0,fCentralityCL1);
  FillHistogramTH2(fOutputContainer,"hCentralityV0AvsV0C",fCentralityV0A,fCentralityV0C);
  FillHistogramTH2(fOutputContainer,"hCentralityZNAvsZNC",fCentralityZNA,fCentralityZNC);

  FillHistogramTH1(fOutputContainer,"hVertexZSelectEvent" ,fVertex[2]);
  FillHistogramTH1(fOutputContainer,"hEventSummary",2);//selected event

  if(fIsFlowTask){
    //fFlowQnVectorMgr->GetQnVectorList()->Print("",-1);

    if(fHarmonics < 0){
      AliError(Form("Qn Flow vector correction flag is ON, but fHarmonics is not set. (it is %d now).",fHarmonics));
      return;
    }

    TList* qnlist = fFlowQnVectorMgr->GetQnVectorList();

    const AliQnCorrectionsQnVector *QnVectorTPCDet[3];
    Double_t TPCEP[3] = {};
    for(Int_t i=0;i<3;i++){
      QnVectorTPCDet[i] = GetQnVectorFromList(qnlist,Form("%s%s",fTPCEPName[i].Data(),fQnEstimator.Data()),"","");
      if(QnVectorTPCDet[i]) TPCEP[i] = QnVectorTPCDet[i]->EventPlane(fHarmonics);
      if(TPCEP[i] < 0) TPCEP[i] += 2./(Double_t) fHarmonics * TMath::Pi();
      FillHistogramTH2(fOutputContainer,Form("hCentrality%svsEventPlane%s%s",fEstimator.Data(),fTPCEPName[i].Data(),fQnEstimator.Data()),fCentralityMain,TPCEP[i]);
      AliInfo(Form("harmonics %d | TPC sub detector name %s : event plane = %f (rad).",fHarmonics,fTPCEPName[i].Data(),TPCEP[i]));
    }

    const AliQnCorrectionsQnVector *QnVectorV0Det[3];
    Double_t V0EP[3]  = {};
    for(Int_t i=0;i<3;i++){
      QnVectorV0Det[i]  = GetQnVectorFromList(qnlist,Form("%s%s",fV0EPName[i].Data(),fQnEstimator.Data()),"","");
      if(QnVectorV0Det[i]) V0EP[i] = QnVectorV0Det[i]->EventPlane(fHarmonics);
      if(V0EP[i] < 0)  V0EP[i]  += 2./(Double_t) fHarmonics * TMath::Pi();
      FillHistogramTH2(fOutputContainer,Form("hCentrality%svsEventPlane%s%s",fEstimator.Data(),fV0EPName[i].Data(),fQnEstimator.Data()),fCentralityMain,V0EP[i]);
      AliInfo(Form("harmonics %d | V0  sub detector name %s : event plane = %f (rad).",fHarmonics,fV0EPName[i].Data() ,V0EP[i]));
    }

    //0 < event plane < 2*pi/fHarmonics.
    //fEventPlane = TPCEP[0];//full acceptance of TPC
    fEventPlane = V0EP[0];//full V0

    Double_t Q1[2] = {};
    Double_t Q2[2] = {};
    Double_t Q3[2] = {};

    Q1[0] = QnVectorV0Det[1]->Qx(fHarmonics);//V0A
    Q1[1] = QnVectorV0Det[1]->Qy(fHarmonics);//V0A
    Q2[0] = QnVectorV0Det[2]->Qx(fHarmonics);//V0C
    Q2[1] = QnVectorV0Det[2]->Qy(fHarmonics);//V0C

    Q3[0] = QnVectorTPCDet[0]->Qx(fHarmonics);//full acceptance of TPC
    Q3[1] = QnVectorTPCDet[0]->Qy(fHarmonics);//full acceptance of TPC

    fQVector1.Set(Q1[0],Q1[1]);
    fQVector2.Set(Q2[0],Q2[1]);

    TVector2 QVector3(Q3[0],Q3[1]);//full acceptance of TPC

    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsQ1x",fEstimator.Data()),fCentralityMain,Q1[0]);//V0A
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsQ1y",fEstimator.Data()),fCentralityMain,Q1[1]);//V0A
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsQ2x",fEstimator.Data()),fCentralityMain,Q2[0]);//V0C
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsQ2y",fEstimator.Data()),fCentralityMain,Q2[1]);//V0C
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsQ3x",fEstimator.Data()),fCentralityMain,Q3[0]);//full acceptance of TPC
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsQ3y",fEstimator.Data()),fCentralityMain,Q3[1]);//full acceptance of TPC

    //Double_t sp = fQ1[0] * fQ2[0] + fQ1[1] * fQ2[1];//scalar product between Q1 vector and Q2 vector
    Double_t sp12 = fQVector1 * fQVector2;//scalar product between Q1 vector and Q2 vector
    Double_t sp23 = fQVector2 *  QVector3;//scalar product between Q2 vector and Q2 vector
    Double_t sp31 =  QVector3 * fQVector1;//scalar product between Q3 vector and Q1 vector

    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsSPQ1Q2",fEstimator.Data()),fCentralityMain,sp12);//mean value of sp is denominator of vn{SP}
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsSPQ2Q3",fEstimator.Data()),fCentralityMain,sp23);//mean value of sp is denominator of vn{SP}
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsSPQ3Q1",fEstimator.Data()),fCentralityMain,sp31);//mean value of sp is denominator of vn{SP}
    AliInfo(Form("Q1x = %e , Q1y = %e , Q2x = %e , Q2y = %e , Q3x = %e , Q3y = %e ,  SP12 = %e ,  SP23 = %e ,  SP31 = %e",Q1[0],Q1[1],Q2[0],Q2[1],Q3[0],Q3[1],sp12,sp23,sp31));

    const Double_t delta = 2. * TMath::Pi() / Double_t(fHarmonics) / 12.;
    fEPBin = (Int_t)((fEventPlane) / delta);//it should be 0-11.
    if(fEPBin < 0)  fEPBin =  0;//protection to avoid fEPBin = -1.
    if(fEPBin > 11) fEPBin = 11;//protection to avoid fEPBin = 12.

    //for event plane resolution
    //cos V0(main)-TPCA-TPCC
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsCosDeltaEventPlane%s%s",fEstimator.Data(),fV0EPName[0].Data() ,fTPCEPName[1].Data()),fCentralityMain,TMath::Cos(fHarmonics * (V0EP[0]  - TPCEP[1])));
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsCosDeltaEventPlane%s%s",fEstimator.Data(),fV0EPName[0].Data() ,fTPCEPName[2].Data()),fCentralityMain,TMath::Cos(fHarmonics * (V0EP[0]  - TPCEP[2])));
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsCosDeltaEventPlane%s%s",fEstimator.Data(),fTPCEPName[1].Data(),fTPCEPName[2].Data()),fCentralityMain,TMath::Cos(fHarmonics * (TPCEP[1] - TPCEP[2])));

    //cos TPC(main)-V0A-V0C
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsCosDeltaEventPlane%s%s",fEstimator.Data(),fTPCEPName[0].Data(),fV0EPName[1].Data()),fCentralityMain,TMath::Cos(fHarmonics * (TPCEP[0] - V0EP[1])));
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsCosDeltaEventPlane%s%s",fEstimator.Data(),fTPCEPName[0].Data(),fV0EPName[2].Data()),fCentralityMain,TMath::Cos(fHarmonics * (TPCEP[0] - V0EP[2])));
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsCosDeltaEventPlane%s%s",fEstimator.Data(),fV0EPName[1].Data() ,fV0EPName[2].Data()),fCentralityMain,TMath::Cos(fHarmonics * (V0EP[1]  - V0EP[2])));
  }
  else fEPBin = 0;
 
  AliInfo(Form("Collision system = %d | fCentralityMain estimated by %s = %f %% | Zvtx = %f cm , fZvtx = %d | Harmonics = %d , fEventPlane = %f (rad.) , fEPBin = %d |",fCollisionSystem,fEstimator.Data(),fCentralityMain,fVertex[2],fZvtx,fHarmonics,fEventPlane,fEPBin));

  if(fRunNumber != fEvent->GetRunNumber()){ // Check run number
    fRunNumber = fEvent->GetRunNumber();
    fPHOSGeo = GetPHOSGeometry();
  }

  fPIDResponse = fInputHandler->GetPIDResponse();
  if(!fPIDResponse){
    AliWarning("fPIDResponse does not exist! This is not crucial in photon analysis in calorimeters.");
  }

  //track QA
  TrackQA();

  if(fIsMC){//fill cluster occupancy in M.C. to obtain total(i.e. HIJNG/HIJING+JJ) Ncluster.
    Int_t multPHOSClustAll = 0;
    for(Int_t i1=0;i1<fPHOSClusterArray->GetEntriesFast();i1++){
      AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
      if(!fPHOSClusterCuts->AcceptPhoton(ph)) continue;
      multPHOSClustAll++;
    }
    FillHistogramTH2(fOutputContainer,Form("hCentrality%svsPHOSClusterMultiplicityMC",fEstimator.Data()),fCentralityMain,multPHOSClustAll);
  }

  if(!fPHOSEvents[fZvtx][fEPBin]) fPHOSEvents[fZvtx][fEPBin] = new TList();
  TList *prevPHOS = fPHOSEvents[fZvtx][fEPBin];

  GetEmbeddedMCInfo();
  ProcessMC();
  SetWeightToClusters();

  AliAnalysisTaskPHOSPi0EtaToGammaGamma::ClusterQA();
  FillPhoton();
  if(fParticleName.Contains("Pi0") || fParticleName.Contains("Eta")){
    FillMgg();
    FillMixMgg();
  }
  AliAnalysisTaskPHOSPi0EtaToGammaGamma::EstimatePIDCutEfficiency();

  //Now we either add current events to stack or remove
  //If no photons in current event - no need to add it to mixed
  if(fPHOSClusterArray->GetEntriesFast() > 0){
    TClonesArray *clone = new TClonesArray(*fPHOSClusterArray);

    //prevPHOS->AddFirst(fPHOSClusterArray);
    prevPHOS->AddFirst(clone);
    //fPHOSClusterArray=0;
    //clone = 0;

    if(prevPHOS->GetSize() > fNMixed){//Remove redundant events
      TClonesArray * tmp = static_cast<TClonesArray*>(prevPHOS->Last());
      prevPHOS->RemoveLast();
      delete tmp;
      tmp = NULL;
    }
  }

  PostData(1, fOutputContainer);
}
//________________________________________________________________________
void AliAnalysisTaskPHOSEmbeddingEfficiency::Terminate(Option_t *option) 
{
  //Called once at the end of the query
  //In principle, this function is not needed...

  AliInfo(Form("%s is done.",GetName()));

}
//________________________________________________________________________
void AliAnalysisTaskPHOSEmbeddingEfficiency::ProcessMC()
{
  //This is for analyzing general purpose MC such as pure PYTHIA, HIJING, DPMJET, PHOJET and so on.
  //get MC information
  TF1 *f1weight = 0x0;
  if(fParticleName.Contains("Pi0"))        f1weight = GetAdditionalPi0PtWeightFunction(fCentralityMain);
  else if(fParticleName.Contains("Eta"))   f1weight = GetAdditionalEtaPtWeightFunction(fCentralityMain);
  else if(fParticleName.Contains("Gamma")) f1weight = GetAdditionalGammaPtWeightFunction(fCentralityMain);

  AliAODMCParticle *p_origin = (AliAODMCParticle*)fMCArray->At(0);//0 is always generated particle by AliGenBox.
  Double_t pT_origin = p_origin->Pt();
  Double_t weight = f1weight->Eval(pT_origin) * pT_origin;

  Int_t genID = -1;
  Double_t pT=0, rapidity=0, phi=0;
  Int_t pdg = 0;
  TString parname = "";
  TString genname = "";

  const Int_t Ntrack = fMCArray->GetEntriesFast();

  for(Int_t i=0;i<Ntrack;i++){
    AliAODMCParticle *p = (AliAODMCParticle*)fMCArray->At(i);
    genID = p->GetGeneratorIndex();
    pT = p->Pt();
    rapidity = p->Y();
    phi = p->Phi();
    pdg = p->PdgCode();

    //rapidity is Y(), but, pseudo-rapidity is Eta();

    if(pT < 1e-3) continue;//reject below 1 MeV
    if(TMath::Abs(rapidity) > 0.5) continue;

    //if(R(p) > 1.0) continue;

    if(pdg==111){//pi0
      parname = "Pi0";
    }
    else if(pdg==221){//eta
      parname = "Eta";

    }
    else if(pdg==22){//gamma
      parname = "Gamma";
    }
    else{
      continue;
    }

    //Double_t eta = p->Eta();
    //Printf("particle %d is generated at eta = %f , phi = %f and pT = %f.",pdg,eta,phi,pT);


    FillHistogramTH1(fOutputContainer,Form("hGenEmbedded%sPt"    ,parname.Data()),pT          ,weight);
    FillHistogramTH2(fOutputContainer,Form("hGenEmbedded%sEtaPhi",parname.Data()),rapidity,phi,weight);
    FillHistogramTH2(fOutputContainer,Form("hGenEmbedded%sEtaPt" ,parname.Data()),rapidity,pT ,weight);

  }//end of generated particle loop

}
//________________________________________________________________________
void AliAnalysisTaskPHOSEmbeddingEfficiency::SetWeightToClusters()
{
  TF1 *f1weight = 0x0;
  if(fParticleName.Contains("Pi0"))        f1weight = GetAdditionalPi0PtWeightFunction(fCentralityMain);
  else if(fParticleName.Contains("Eta"))   f1weight = GetAdditionalEtaPtWeightFunction(fCentralityMain);
  else if(fParticleName.Contains("Gamma")) f1weight = GetAdditionalGammaPtWeightFunction(fCentralityMain);

  AliAODMCParticle *p_origin = (AliAODMCParticle*)fMCArray->At(0);//0 is always generated particle by AliGenBox.
  Double_t pT_origin = p_origin->Pt();
  Double_t weight = f1weight->Eval(pT_origin) * pT_origin;

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();
  for(Int_t i=0;i<multClust;i++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i);
    ph->SetWeight(weight);

  }

}
//________________________________________________________________________
Int_t AliAnalysisTaskPHOSEmbeddingEfficiency::FindCommonParent(Int_t iPart, Int_t jPart)
{
  //check if there is a common parent for particles i and j
  // -1: no common parent or wrong iPart/jPart

  Int_t ntrack = fMCArray->GetEntriesFast();
  if(iPart==-1 || iPart>=ntrack || jPart==-1 || jPart>=ntrack) return -1;

  Int_t iprim1 = iPart;

  while(iprim1>-1){
    Int_t iprim2=jPart;

    while(iprim2>-1){
      if(iprim1==iprim2) return iprim1;
      //iprim2 = GetParticle(iprim2)->GetMother();
      iprim2 = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iprim2))->GetMother();
    }

    //iprim1 = GetParticle(iprim1)->GetMother();
    //iprim1 = (AliAODMCParticle*)(fMCArray->At(iprim1))->GetMother();
    iprim1 = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iprim1))->GetMother();
  }

  return -1;
}
//________________________________________________________________________
void AliAnalysisTaskPHOSEmbeddingEfficiency::GetEmbeddedMCInfo()
{
  fMCArray = 0x0;
  //if(fParticleName.Contains("Pi0"))        fMCArray = (TClonesArray*)fEvent->FindListObject(Form("%s_pi0"  ,AliAODMCParticle::StdBranchName()));
  //else if(fParticleName.Contains("Eta"))   fMCArray = (TClonesArray*)fEvent->FindListObject(Form("%s_eta"  ,AliAODMCParticle::StdBranchName()));
  //else if(fParticleName.Contains("Gamma")) fMCArray = (TClonesArray*)fEvent->FindListObject(Form("%s_gamma",AliAODMCParticle::StdBranchName()));

  fMCArray = (TClonesArray*)fEvent->FindListObject(Form("%s_%s",AliAODMCParticle::StdBranchName(),fParticleName.Data()));

  if(!fMCArray){
    AliError("Could not retrieve AOD event!");
    return;
  }
}
//________________________________________________________________________
void AliAnalysisTaskPHOSEmbeddingEfficiency::FillPhoton() 
{
  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  Double_t pT=0,energy=0;
  Double_t phi = -999, dphi = -999.;
  Double_t eff=1;
  TF1 *f1tof = GetTOFCutEfficiencyFunction();
  Double_t trgeff=1;
  TF1 *f1trg = GetTriggerEfficiencyFunction();
  Double_t value[4] = {};
  Double_t sp1 = -999;
  Double_t sp2 = -999;

  Double_t weight = 1.;
  Int_t primary = -1;

  for(Int_t iph=0;iph<multClust;iph++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(iph);
    if(!fPHOSClusterCuts->AcceptPhoton(ph)) continue;

     if(fIsPHOSTriggerAnalysis){
      if(!fPHOSTriggerHelper->IsOnActiveTRUChannel(ph)) continue;
      if(!fIsMC && !ph->IsTrig()) continue;//it is meaningless to focus on photon without fired trigger in PHOS triggered data.
    }

    weight = 1.;
    if(fIsMC){
      primary = ph->GetPrimary();
      weight = ph->GetWeight();
    }

    pT = ph->Pt();
    energy = ph->Energy();
    phi = ph->Phi();

    if(fUseCoreEnergy){
      pT = (ph->GetMomV2())->Pt();
      energy = (ph->GetMomV2())->Energy();
      phi = (ph->GetMomV2())->Phi();
    }

    eff = f1tof->Eval(energy);

    //0 < photon phi < 2pi
    if(phi < 0) phi += TMath::TwoPi();
    TVector2 vg(TMath::Cos(fHarmonics * phi),TMath::Sin(fHarmonics * phi));

    if(fIsFlowTask){
      dphi = DeltaPhiIn0Pi(phi - fEventPlane);
      sp1 = vg * fQVector1;
      sp2 = vg * fQVector2;
    }
    else{
      dphi = phi;
      sp1 = 0;
      sp2 = 0;
    }

    value[0] = pT;
    value[1] = TMath::Cos(fHarmonics * dphi);
    value[2] = sp1;//reserved by sp
    value[3] = sp2;//reserved by sp

    //FillHistogramTH2(fOutputContainer,"hPhotonPt",pT,TMath::Cos(fHarmonics * dphi),weight);
    FillSparse(fOutputContainer,"hSparsePhoton",value,weight * 1/trgeff);

    if(ph->IsTOFOK()){
      //FillHistogramTH2(fOutputContainer,"hPhotonPt_TOF",pT,TMath::Cos(fHarmonics * dphi),1/eff * weight);
      FillSparse(fOutputContainer,"hSparsePhoton_TOF",value,1/eff * weight * 1/trgeff);
    }

  }//end of cluster loop

}
//________________________________________________________________________
void AliAnalysisTaskPHOSEmbeddingEfficiency::FillMgg() 
{
  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();
  TLorentzVector p12, p12core;

  Double_t m12=0,pt12=0,asym=0;
  Double_t e1=0,e2=0;
  Double_t phi = -999, dphi = -999.;

  Double_t eff1=1, eff2=1, eff12=1;
  TF1 *f1tof = GetTOFCutEfficiencyFunction();

  Double_t trgeff1=1;
  Double_t trgeff2=1;
  Double_t trgeff12=1;
  TF1 *f1trg = GetTriggerEfficiencyFunction();

  Double_t value[6] = {};
  Double_t sp1 = -999;
  Double_t sp2 = -999;

  Double_t weight = 1., w1 = 1., w2 = 1.;

  Int_t primary1 = -1;
  Int_t primary2 = -1;
  Double_t TrueK0SPt = 0;
  Double_t TrueL0Pt = 0;

  Int_t commonID = -1;

  for(Int_t i1=0;i1<multClust-1;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(fIsPHOSTriggerAnalysis && !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph1)) continue;

    for(Int_t i2=i1+1;i2<multClust;i2++){
      AliCaloPhoton *ph2 = (AliCaloPhoton*)fPHOSClusterArray->At(i2);
      if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;
      if(fIsPHOSTriggerAnalysis && !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph2)) continue;

      if(!fIsMC && fIsPHOSTriggerAnalysis && (!ph1->IsTrig() && !ph2->IsTrig())) continue;//it is meaningless to reconstruct invariant mass with FALSE-FALSE combination in PHOS triggered data.

      e1 = ph1->Energy();
      e2 = ph2->Energy();

      p12  = *ph1 + *ph2;
      m12  = p12.M();
      pt12 = p12.Pt();
      phi  = p12.Phi();
      asym = TMath::Abs((ph1->Energy()-ph2->Energy())/(ph1->Energy()+ph2->Energy()));//always full energy

      if(fUseCoreEnergy){
        p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
        m12     = p12core.M();
        pt12    = p12core.Pt();
        phi     = p12core.Phi();

        e1 = (ph1->GetMomV2())->Energy();
        e2 = (ph2->GetMomV2())->Energy();
      }

      eff1 = f1tof->Eval(e1);
      eff2 = f1tof->Eval(e2);
      eff12 = eff1 * eff2;

      if(!fIsMC && fIsPHOSTriggerAnalysis){
        trgeff1  = f1trg->Eval(e1);
        trgeff2  = f1trg->Eval(e2);
        trgeff12 = trgeff1 + trgeff2 - (trgeff1 * trgeff2);//logical OR
      }

      weight = 1.;
      if(fIsMC){
        w1= ph1->GetWeight();
        primary1 = ph1->GetPrimary();

        w2 = ph2->GetWeight();
        primary2 = ph2->GetPrimary();

        weight = w1;//common weighting to all generated particles in embedding.

        //commonID = FindCommonParent(primary1,primary2);
        //if(commonID > -1) weight = w1;
        //else weight = w1*w2;

      }//end of if fIsMC

      if(phi < 0) phi += TMath::TwoPi();

      TVector2 vgg(TMath::Cos(fHarmonics * phi),TMath::Sin(fHarmonics * phi));

      if(fIsFlowTask){
        dphi = DeltaPhiIn0Pi(phi - fEventPlane);
        sp1 = vgg * fQVector1;
        sp2 = vgg * fQVector2;
      }
      else{
        dphi = phi;
        sp1 = 0;
        sp2 = 0;
      }

      value[0] = m12;
      value[1] = pt12;
      value[2] = asym;
      value[3] = TMath::Cos(fHarmonics * dphi);
      value[4] = sp1;//reserved by sp
      value[5] = sp2;//reserved by sp

      if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMgg_M%d%d",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12,weight * 1/trgeff12);
      //FillHistogramTH3(fOutputContainer,"hMgg",m12,pt12,TMath::Cos(fHarmonics * dphi),weight);
      FillSparse(fOutputContainer,"hSparseMgg",value,weight * 1/trgeff12);

      //FillHistogramTH2(fOutputContainer,"hAsymvsMgg",asym,m12,weight);
      //if(0.12 < m12 && m12 < 0.15) FillHistogramTH2(fOutputContainer,"hAsymvsPt",asym,pt12,weight);
      //if(asym < 0.8) FillHistogramTH3(fOutputContainer,"hMgg_asym08",m12,pt12,TMath::Cos(fHarmonics * dphi),weight);

      if(ph1->IsTOFOK() && ph2->IsTOFOK()){
        //FillHistogramTH3(fOutputContainer,"hMgg_TOF",m12,pt12,TMath::Cos(fHarmonics * dphi),1/eff12 * weight);

        FillSparse(fOutputContainer,"hSparseMgg_TOF",value,1/eff12 * weight * 1/trgeff12);

        if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMgg_M%d%d_TOF",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12,1/eff12 * weight * 1/trgeff12);

        //if(asym < 0.8) FillHistogramTH3(fOutputContainer,"hMgg_TOF_asym08",m12,pt12,TMath::Cos(fHarmonics * dphi),1/eff12 * weight);

      }//end of TOF cut

    }//end of ph2

  }//end of ph1

}
//________________________________________________________________________
void AliAnalysisTaskPHOSEmbeddingEfficiency::FillMixMgg() 
{
  TList *prevPHOS = fPHOSEvents[fZvtx][fEPBin];

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  TLorentzVector p12, p12core;
  Double_t m12=0,pt12=0,asym=0;
  Double_t phi = -999, dphi = -999.;
  Double_t weight = 1., w1 = 1., w2 = 1.;

  Double_t eff1=1, eff2=1, eff12=1;
  Double_t e1=0,e2=0;
  TF1 *f1tof = GetTOFCutEfficiencyFunction();
  Double_t trgeff1=1;
  Double_t trgeff2=1;
  Double_t trgeff12=1;
  TF1 *f1trg = GetTriggerEfficiencyFunction();

  Double_t value[4] = {};
  Double_t sp1 = -999;
  Double_t sp2 = -999;

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;

    for(Int_t ev=0;ev<prevPHOS->GetSize();ev++){
      TClonesArray *mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev));

      for(Int_t i2=0;i2<mixPHOS->GetEntriesFast();i2++){
        AliCaloPhoton *ph2 = (AliCaloPhoton*)mixPHOS->At(i2);
        if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;

        if(!fIsMC && fIsPHOSTriggerAnalysis && (!ph1->IsTrig() && !ph2->IsTrig())) continue;//it is meaningless to reconstruct invariant mass with FALSE-FALSE combination in PHOS triggered data.

        e1 = ph1->Energy();
        e2 = ph2->Energy();

        p12  = *ph1 + *ph2;
        m12  = p12.M();
        pt12 = p12.Pt();
        phi  = p12.Phi();
        asym = TMath::Abs((ph1->Energy()-ph2->Energy())/(ph1->Energy()+ph2->Energy()));

        if(fUseCoreEnergy){
          p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
          m12     = p12core.M();
          pt12    = p12core.Pt();
          phi     = p12core.Phi();

          e1 = (ph1->GetMomV2())->Energy();
          e2 = (ph2->GetMomV2())->Energy();
        }

        eff1 = f1tof->Eval(e1);
        eff2 = f1tof->Eval(e2);
        eff12 = eff1 * eff2;
        weight = 1.;

        if(!fIsMC && fIsPHOSTriggerAnalysis){
          trgeff1  = f1trg->Eval(e1);
          trgeff2  = f1trg->Eval(e2);
          trgeff12 = trgeff1 + trgeff2 - (trgeff1 * trgeff2);//logical OR
        }

        if(fIsMC){
          w1= ph1->GetWeight();
          w2 = ph2->GetWeight();

          weight = w1*w2;

        }//end of if fIsMC

        if(phi < 0) phi += TMath::TwoPi();
        TVector2 vgg(TMath::Cos(fHarmonics * phi),TMath::Sin(fHarmonics * phi));

        if(fIsFlowTask){
          dphi = DeltaPhiIn0Pi(phi - fEventPlane);
          sp1 = vgg * fQVector1;
          sp2 = vgg * fQVector2;
        }
        else{
          dphi = phi;
          sp1 = 0;
          sp2 = 0;
        }

        value[0] = m12;
        value[1] = pt12;
        value[2] = asym;
        value[3] = TMath::Cos(fHarmonics * dphi);

        if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMixMgg_M%d%d",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12);
        FillSparse(fOutputContainer,"hSparseMixMgg",value,weight * 1/trgeff12);

        //FillHistogramTH3(fOutputContainer,"hMixMgg",m12,pt12,TMath::Cos(fHarmonics * dphi),weight);
        //if(asym < 0.8) FillHistogramTH3(fOutputContainer,"hMixMgg_asym08",m12,pt12,TMath::Cos(fHarmonics * dphi),weight);

        if(ph1->IsTOFOK() && ph2->IsTOFOK()){
          FillSparse(fOutputContainer,"hSparseMixMgg_TOF",value,1/eff12 * weight * 1/trgeff12);
          if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMixMgg_M%d%d_TOF",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12,1/eff12);

          //FillHistogramTH3(fOutputContainer,"hMixMgg_TOF",m12,pt12,TMath::Cos(fHarmonics * dphi),1/eff12 * weight);
          //if(asym < 0.8) FillHistogramTH3(fOutputContainer,"hMixMgg_TOF_asym08",m12,pt12,TMath::Cos(fHarmonics * dphi),1/eff12 * weight);

        }//end of TOF cut

      }//end of ph2

    }//end of mix

  }//end of ph1

}

//________________________________________________________________________
//________________________________________________________________________
