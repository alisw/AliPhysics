#include "TF1.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TMath.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "THashList.h"
#include "TChain.h"

#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliCaloPhoton.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliPHOSGeometry.h"
#include "AliOADBContainer.h"
#include "AliPHOSCalibData.h"
#include "AliPHOSAodCluster.h"
#include "AliAODMCParticle.h"

#include "AliAnalysisTaskPHOSObjectCreator.h"
#include "AliAnalysisTaskPHOSEmbeddedDiffObjectCreator.h"

ClassImp(AliAnalysisTaskPHOSEmbeddedDiffObjectCreator)

//________________________________________________________________________
AliAnalysisTaskPHOSEmbeddedDiffObjectCreator::AliAnalysisTaskPHOSEmbeddedDiffObjectCreator(const char *name):
  AliAnalysisTaskPHOSObjectCreator(name),
  fParticleName(""),
  fPHOSCalibData(0x0)
{
  // Constructor

  fPHOSCalibData = new AliPHOSCalibData();
  for(Int_t module=1; module<=5; module++) {
    for(Int_t column=1; column<=56; column++) {
      for(Int_t row=1; row<=64; row++) {
        fPHOSCalibData->SetADCchannelEmc(module,column,row,1.);
      }
    }
  } 

 
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, THashList::Class());

}
//________________________________________________________________________
AliAnalysisTaskPHOSEmbeddedDiffObjectCreator::~AliAnalysisTaskPHOSEmbeddedDiffObjectCreator()
{
  //destructor

  delete fPHOSCalibData;


}
//________________________________________________________________________
void AliAnalysisTaskPHOSEmbeddedDiffObjectCreator::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  Init();//initialize PHOS object array

  //fOutputContainer = new THashList();
  //fOutputContainer->SetOwner(kTRUE);
  AliAnalysisTaskPHOSObjectCreator::UserCreateOutputObjects();
  PostData(1,fOutputContainer);
}
//________________________________________________________________________
void AliAnalysisTaskPHOSEmbeddedDiffObjectCreator::UserExec(Option_t *option) 
{
  // Main loop
  // Called for each event

  fEvent = dynamic_cast<AliVEvent*>(InputEvent());
  if(!fEvent){
    AliError("event is not available.");
    return;
  }

  fAODEvent = dynamic_cast<AliAODEvent*>(fEvent);
  //fESDEvent = dynamic_cast<AliESDEvent*>(fEvent);
  //ESD is not supported at all.

  fMCArrayAOD = 0x0;
  fMCArrayAOD = (TClonesArray*)fEvent->FindListObject(Form("%s_%s",AliAODMCParticle::StdBranchName(),fParticleName.Data()));
  if(!fMCArrayAOD){
    AliError("Could not retrieve MC array!");
    return;
  }

  if(fRunNumber != fEvent->GetRunNumber()) { // Check run number
    fRunNumber = fEvent->GetRunNumber();
    SetGeometry();
//    InitBadMap();
  }

  const AliVVertex *vVertex = fEvent->GetPrimaryVertex();
  fVertex[0] = vVertex->GetX();
  fVertex[1] = vVertex->GetY();
  fVertex[2] = vVertex->GetZ();

  fPHOSObjectArray->Clear();

  TLorentzVector p1,p1core;

  //Double_t energy=0, tof=-999, M20=0, M02=0, R2=999, coreE=0, coreR2=999;
  Double_t energy=0, coreE=0;
  Double_t tof=-999, M20=0, M02=0, R2=999, coreR2=999;
  Double_t r=999;

  Int_t relId[4] = {};
  Int_t module=-1, cellx=-1, cellz=-1;
  Int_t digMult=0;
  Float_t position[3] = {}; 
  Double_t distance=0;
  Int_t inPHOS=0;
  //const Double_t vtx000[3] = {0.,0.,0.};//this is because of generated pi0 from 000.

  TClonesArray * clustersEmb = (TClonesArray*)   fEvent->FindListObject(Form("Embedded%sCaloClusters",fParticleName.Data())) ;
  AliAODCaloCells * cellsEmb = (AliAODCaloCells*)fEvent->FindListObject(Form("Embedded%sPHOScells"   ,fParticleName.Data())) ;
  const Int_t multClustEmb = clustersEmb->GetEntriesFast();
  AliInfo(Form("%d clusters (PHOS CPV) are detected after embedding.",multClustEmb));

  TClonesArray * clustersUE = (TClonesArray*)fAODEvent->GetCaloClusters();//clusters in underlying event.
  Int_t multClustUE = fAODEvent->GetNumberOfCaloClusters();//real data is treated as underlying event.
  //AliVCaloCells *cellsUE = dynamic_cast<AliVCaloCells*>(fEvent->GetPHOSCells());

  if(!fUsePHOSTender) CalibrateEmbeddedClusters(clustersEmb, cellsEmb);

  for(Int_t iclu=0; iclu<multClustEmb; iclu++){//embedded cluster loop
    AliAODCaloCluster *cluster = (AliAODCaloCluster*)clustersEmb->At(iclu);

    if(cluster->GetType() != AliVCluster::kPHOSNeutral
        || cluster->E() < 0.1 // Emin cut
      ) continue;

    if(cluster->GetLabel() < 0) continue;

    Bool_t IsSameFound = kFALSE;

    for(Int_t icluUE=0;icluUE<multClustUE;icluUE++){
      AliAODCaloCluster *clusterUE = (AliAODCaloCluster*)clustersUE->At(icluUE);

      if(clusterUE->GetType() != AliVCluster::kPHOSNeutral
          || clusterUE->E() < 0.1 // noise cut
        ) continue;


//      cout << "UE E = " << clusterUE->E() << " , Ncell = " << clusterUE->GetNCells() << " , label = " << clusterUE->GetLabel() << endl;

      Bool_t issame = IsSameCluster(cluster,clusterUE);

      if(issame){
        IsSameFound = kTRUE;
        //cout << "same cluster is found!" << endl;
        break;
      }

    }//end of cluster loop

    if(IsSameFound) continue;
    AliInfo(Form("Particle : %s , embedded cluster is found. E = %f GeV , Ncell = %d , MC label = %d.",fParticleName.Data(),cluster->E(), cluster->GetNCells(), cluster->GetLabel()));

    distance = cluster->GetDistanceToBadChannel();//in cm

    Int_t Ncell = cluster->GetNCells();
    for(Int_t i=0;i<cluster->GetNCells();i++){
      Int_t absId = cluster->GetCellAbsId(i);
      Double_t amp = cellsEmb->GetCellAmplitude(absId);
      //cout << "cell amplitude = " << amp << endl;

      if(amp < 2e-6) //less than 2 keV
        Ncell--;
    }
    //cluster->SetNCells(Ncell);

    //if(distance < fMinDistBC) continue;

    cluster->GetPosition(position);
    TVector3 global1(position);
    fPHOSGeo->GlobalPos2RelId(global1,relId);

    module = relId[0];
    cellx  = relId[2];
    cellz  = relId[3];

    if(fIsM4Excluded && module==4) continue;

    if(module < 1 || 4 < module){
      AliError(Form("Wrong module number %d",module));
      return;
    }

    if(!fUsePHOSTender && !IsGoodChannel("PHOS",module,cellx,cellz)) continue;

    energy = cluster->E();
    digMult = cluster->GetNCells();
    tof = cluster->GetTOF();
    Int_t primary = cluster->GetLabel();

    ////if simulated photon energy is less than a half of a cluster energy, reject such cluster.
    //AliAODMCParticle *p = (AliAODMCParticle*)fMCArrayAOD->At(primary);
    //Double_t Etrue = p->E();//energy of photon
    //if(Etrue < 0.5 * energy){
    //  AliInfo(Form("energy of photon = %e GeV , reconstructed cluster energy = %e GeV. Eture < 0.5 * Erec. reject this cluster.",Etrue,energy));
    //  continue;
    //}

    cluster->GetMomentum(p1,fVertex);
    cluster->GetMomentum(p1core,fVertex);

    new((*fPHOSObjectArray)[inPHOS]) AliCaloPhoton(p1.Px(),p1.Py(),p1.Pz(),p1.E());
    AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSObjectArray->At(inPHOS); 
    ph->SetCluster(cluster);
    ph->SetModule(module);
    ph->SetNCells(digMult); 
    ph->SetTime(tof);//unit of second
    ph->SetTOFBit(TMath::Abs(tof*1e+9) < fBunchSpace/2.);
    ph->SetDistToBadfp(distance/2.2);//in unit of cells with floating point. 2.2 cm is crystal size
    ph->SetEMCx((Double_t)position[0]);
    ph->SetEMCy((Double_t)position[1]);
    ph->SetEMCz((Double_t)position[2]);

    //Bool_t sure = kFALSE;
    //Int_t primary = FindPrimary(ph,sure);

    ph->SetPrimary(primary);//MC label is meaningless in embedding.
    ph->SetWeight(1.);

    //When PHOSTender is applied, GetM20(), GetM02() returns M20,M02 of core of cluster respectively.
    M20 = cluster->GetM20();//M20 is short axis of elliptic shower shape.
    M02 = cluster->GetM02();//M02 is long  axis of elliptic shower shape.
    R2 = cluster->GetDispersion();//full dispersion
    coreE = cluster->GetCoreEnergy();
    coreR2 = cluster->Chi2();//core dispersion
    r = cluster->GetEmcCpvDistance();

    ph->SetLambdas(M20,M02);
    ph->SetNsigmaCPV(r);
    ph->SetNsigmaFullDisp(TMath::Sqrt(R2));
    ph->SetNsigmaCoreDisp(TMath::Sqrt(coreR2));

    p1core *= coreE/energy; //use core energy in PbPb.
    ph->SetMomV2(&p1core);//core energy

    //printf("energy = %e GeV, coreE = %e GeV\n",ph->Energy(),ph->GetMomV2()->Energy());

    inPHOS++;
  }//end of cluster loop


  EstimateSTDCutEfficiency(fPHOSObjectArray);

  const Int_t Nph = fPHOSObjectArray->GetEntries();

  for(Int_t iph=0;iph<Nph;iph++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSObjectArray->At(iph);
    AliVCluster *cluster = (AliVCluster*)ph->GetCluster();
    if(!PassSTDCut(cluster)){
      //AliInfo("cluster is rejected by the standard cluster cut.");
      fPHOSObjectArray->Remove(ph);
    }
  }
  fPHOSObjectArray->Compress();

  AliVEvent* input = InputEvent();
  TObject* outO = input->FindListObject(Form("PHOSEmbeddedDiffClusterArray_%s",fParticleName.Data()));

  if(!outO){
    fPHOSObjectArray->SetName(Form("PHOSEmbeddedDiffClusterArray_%s",fParticleName.Data()));
    input->AddObject(fPHOSObjectArray);
  }


}
//________________________________________________________________________
void AliAnalysisTaskPHOSEmbeddedDiffObjectCreator::Terminate(Option_t *option) 
{
  //Called once at the end of the query
  


}
//________________________________________________________________________
void AliAnalysisTaskPHOSEmbeddedDiffObjectCreator::CalibrateEmbeddedClusters(TClonesArray *clustersEmb, AliAODCaloCells *cellsEmb) 
{
  //instead of PHOSTender, embedded clusters should be de/calibrated here.
  //Embedding study supports only AOD.
  //Thus, this function supports only AOD clusters.
  //this is copied from AliPHOSTenderSupply.cxx/h

  const Double_t logWeight = 4.5;
  TVector3 vertex;
  vertex.SetXYZ(fVertex[0],fVertex[1],fVertex[2]);

  const Int_t multClust = clustersEmb->GetEntriesFast();

  for (Int_t i=0; i<multClust; i++) {
    AliAODCaloCluster *clu = (AliAODCaloCluster*)clustersEmb->At(i);

    if(clu->GetType() != AliVCluster::kPHOSNeutral) continue;
    if(clu->GetLabel() < 0) continue;

    //cout << "energy before calibration = " << clu->E() << endl;

    Float_t  positionOld[3]={};
    clu->GetPosition(positionOld);
    //cout << "before recalibrate x = " << positionOld[0] << " , y = " << positionOld[1] << " , z = " << positionOld[2] << endl;
    TVector3 globalOld(positionOld) ;

    //Apply re-Calibreation
    AliPHOSAodCluster cluPHOS(*clu);
    cluPHOS.Recalibrate(fPHOSCalibData,cellsEmb); // modify the cell energies
    cluPHOS.EvalAll(logWeight,vertex);         // recalculate the cluster parameters
    cluPHOS.SetE(fNonLinCorr->Eval(cluPHOS.E()));// Users's nonlinearity

    Float_t  position[3] = {};
    cluPHOS.GetPosition(position);

    //cout << "after recalibrate x = " << position[0] << " , y = " << position[1] << " , z = " << position[2] << endl;

    clu->SetPosition(position);                       //rec.point position in MARS
    TVector3 global(position) ;
    Int_t relId[4] = {};
    fPHOSGeo->GlobalPos2RelId(global,relId) ;
    Int_t mod  = relId[0] ;
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3] ;
    if(!IsGoodChannel("PHOS",mod,cellX,cellZ) ) {
      clu->SetE(0.) ;
      continue ;
    }
    TVector3 locPosOld; //Use it to re-calculate distance to track
    fPHOSGeo->Global2Local(locPosOld,globalOld,mod) ;

    //Double_t ecore=CoreEnergy(&cluPHOS,cellsEmb);
    Double_t ecore=fNonLinCorr->Eval(CoreEnergy(&cluPHOS));

    clu->SetE(cluPHOS.E());                           //total particle energy
    clu->SetCoreEnergy(ecore);                  //core particle energy

    //Eval FullDispersion
    clu->SetDispersion(TestLambda(clu->E(),cluPHOS.GetM20(),cluPHOS.GetM02())) ;

    //Eval CoreDispersion
    Double_t m02=0.,m20=0.;
    EvalLambdas(&cluPHOS,m02, m20);
    clu->SetChi2(TestCoreLambda(clu->E(),m20,m02));                     //not yet implemented
    clu->SetM02(m02) ;               //second moment M2x
    clu->SetM20(m20) ;               //second moment M2z

    //correct distance to track
    Double_t pttrack=0.;
    Int_t charge=0;
    Double_t dx=999.,dz=999. ;
    TVector3 locPos;
    fPHOSGeo->Global2Local(locPos,global,mod) ;
    Int_t itr=FindTrackMatching(mod,&locPos,dx,dz,pttrack,charge) ;
    clu->SetTrackDistance(dx,dz);
    Double_t r = 999. ; //Big distance
    if(itr > 0) {
      r=TestCPVRun2(dx, dz, pttrack,charge) ;
    }
    clu->SetEmcCpvDistance(r); //Distance in sigmas
    if(itr >= 0){ //there is a track
      //Remove existing
      AliAODTrack * tr = (AliAODTrack*)fAODEvent->GetTrack(itr) ;
      Int_t ntrM = clu->GetNTracksMatched();
      if(ntrM>0){
        AliAODTrack * trStored = (AliAODTrack*)clu->GetTrackMatched(0);
        if(trStored != tr){
          clu->RemoveTrackMatched(trStored);
          clu->AddTrackMatched(tr) ;
        }
      }
      else{
        clu->AddTrackMatched(tr) ;
      }
    }

    Double_t minDist=clu->GetDistanceToBadChannel() ;//Already calculated
    DistanceToBadChannel(mod,&locPos,minDist);
    clu->SetDistanceToBadChannel(minDist) ;

  }

}
//________________________________________________________________________
Bool_t AliAnalysisTaskPHOSEmbeddedDiffObjectCreator::IsSameCluster(AliVCluster * c1, AliVCluster * c2) const
{
  //Compare clusters before and after embedding
  //clusters are the same if
  // - Energy changed less than 1%  (numerical accuracy in reconstruction)
  // - lists of digits are the same

  //cout << "e1 = " << c1->E() << " , e2 = " << c2->E() << endl;
  //cout << "N1 = " << c1->GetNCells() << " , N2 = " << c2->GetNCells() << endl;

  if(c1->GetNCells() != c2->GetNCells()) return kFALSE;//this is different cluster

  //c1 is embedded cluster, c2 is UE.

  //Double_t dE = c1->E() - c2->E();
  //if(TMath::Abs(dE) > 0.01*c1->E()) return kFALSE;//this is different cluster

  Double_t ratioE = c1->E() / c2->E();
  if(TMath::Abs(1. - ratioE) > 0.01) return kFALSE;//this is different cluster
  

  UShort_t *list1 = c1->GetCellsAbsId();
  UShort_t *list2 = c2->GetCellsAbsId();

  for(Int_t i=0; i< c1->GetNCells(); i++){
    //printf("list1[%d] = %d , list2[%d] = %d\n",i,list1[i],i,list2[i]);
    if(list1[i] != list2[i]) return kFALSE;//this is different cluster
  }

  return kTRUE;//they are same clusters

}
//________________________________________________________________________
