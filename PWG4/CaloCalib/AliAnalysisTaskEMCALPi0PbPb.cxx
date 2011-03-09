// $Id$

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "AliESDVertex.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDInputHandler.h"

#include "AliAnalysisTaskEMCALPi0PbPb.h"
#include "AliCentrality.h"

#include "AliEMCALGeoUtils.h"
// example of an analysis task creating a Neutral detector QA Histogram 

ClassImp(AliAnalysisTaskEMCALPi0PbPb)
//________________________________________________________________________
AliAnalysisTaskEMCALPi0PbPb::AliAnalysisTaskEMCALPi0PbPb(const char *name) 
  : AliAnalysisTaskSE(name),
  fOutputDataESD(0),
  fhNEvents(0),
  fhCentVsClusMultEMC(0x0),
  fhCentVsCellMultEMC(0x0),
  fhCentVsClusMultEMCAll(0),
  fhCentVsCellMultEMCAll(0),
  fhCentVsInMVsPtEMC(0),
  fhCentVsInMVsPhiEMC(0),
  fhAsyVsPt(0),
  fhCentVsColuRowEMC(0x0),
  fhCentVsColuRowEnerEMC(0x0),
  fhMggPtEr(0x0),
  fhPtEnSg(0x0)   
{
  // Constructor

  nModuleEMC = 4,
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
  // Initialize the PHOS geometry
}
//________________________________________________________________________
AliAnalysisTaskEMCALPi0PbPb::~AliAnalysisTaskEMCALPi0PbPb() {

  if(fOutputDataESD) delete fOutputDataESD;
  fOutputDataESD = 0;

}
//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
 
  OpenFile(1); 
  fOutputDataESD = new TList();
  fOutputDataESD->SetOwner();
  
  fhNEvents   = new TH1I ("fhNEvents", "Number of events analyzed",1,0.,1.);
  fOutputDataESD-> Add(fhNEvents);

  fhCentVsClusMultEMC = new TH2F *[nModuleEMC]; 
  fhCentVsCellMultEMC = new TH2F *[nModuleEMC];
  fhCentVsColuRowEMC      = new TH3F *[nModuleEMC];
  fhCentVsColuRowEnerEMC  = new TH3F *[nModuleEMC];

  char key[256];
  char title[256];

    fhCentVsClusMultEMCAll = new TH2F("fhCentVsClusMultEMCAll","EMCal Cluster Multiplicity Dis. at centrality",500, 0, 500, 100, 0., 100.);
   fhCentVsCellMultEMCAll = new TH2F ("fhCentVsCellMultEMCAll","EMCal Cell Multiplicity Dis. at centrality",2000, 0, 2000, 100, 0., 100.);
    fOutputDataESD->Add(fhCentVsClusMultEMCAll);
    fOutputDataESD->Add(fhCentVsCellMultEMCAll);

  fhCentVsInMVsPtEMC =  new TH3F("fhCentVsInMVsPtEMC","EMCal Inva Mass Vs Pt Vs central 2 cluster",500, 0, 1, 400, 0, 40, 100, 0, 100);
  fhCentVsInMVsPhiEMC = new TH3F("fhCentVsInMVsPhiEMC","EMCal Inva Mass Vs Phi Vs central 2 cluster", 500, 0, 1, 360, 0, 7, 100, 0, 100); 
  fhAsyVsPt = new TH2F("fhAsyVsPt", "Asymmetry Vs Pt 2 cluster",100, 0, 1, 400, 0, 40);
  
  fhMggPtEr = new TH3F("fhMggPtEr","Invariant mass vs pT vs Energy ratio",100,0,1,100,0,20,100,0,1);
  fhMggPtEr->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}_{}]");
  fhMggPtEr->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  fhMggPtEr->GetZaxis()->SetTitle("E_{r}");
  fhPtEnSg  = new TH3F("fhPtEnSg","photon Pt vs Energy vs Large axis",100,0,20,100,0,50,100,0,50); 
  fhPtEnSg->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fhPtEnSg->GetYaxis()->SetTitle("E [GeV]");
  fhPtEnSg->GetZaxis()->SetTitle("E#sigma_{max} [GeV]");

 fOutputDataESD->Add(fhCentVsInMVsPtEMC);
 fOutputDataESD->Add(fhCentVsInMVsPhiEMC);
 fOutputDataESD->Add(fhAsyVsPt);
 fOutputDataESD->Add(fhMggPtEr);
 fOutputDataESD->Add(fhPtEnSg);

  for(Int_t i=0; i<nModuleEMC; i++) { // for  EMCal
    // for cluster
    sprintf(key,"fhCentVsClusMultEMC_Modu%d",i);
    sprintf(title,"EMCal Module %d Cluster Multiplicity Dis. at centrality",i);
    fhCentVsClusMultEMC[i] = new TH2F (key,title, 500, 0,500, 100, 0.,100.);
   // for Cell
   sprintf(key,"fhCentVsCellMultEMC_Modu%d",i);
   sprintf(title,"EMCal Module %d Cell Multiplicity Dis. at centrality",i);
   fhCentVsCellMultEMC[i] = new TH2F (key,title, 2000, 0, 2000, 100, 0., 100.);
  // for occupancy
   sprintf(key, "fhCentVsColuRowEMC_Modu%d",i); 
   sprintf(title,"entries of cell in Module %d", i); 
   fhCentVsColuRowEMC[i] = new TH3F (key,title, 49, 0, 49, 25, 0, 25, 100, 0., 100.);
  //for energy 
   sprintf(key, "fhCentVsColuRowEMC_Modu%d",i);
   sprintf(title,"energy of cell in Module %d", i);
   fhCentVsColuRowEnerEMC[i] = new TH3F (key,title, 49, 0, 49, 25, 0, 25, 100, 0., 100.);
 
  fOutputDataESD->Add(fhCentVsClusMultEMC[i]);
  fOutputDataESD->Add(fhCentVsCellMultEMC[i]);
  fOutputDataESD->Add(fhCentVsColuRowEMC[i]);
  fOutputDataESD->Add(fhCentVsColuRowEnerEMC[i]);
  } 


 PostData(1, fOutputDataESD); 
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
      AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
      TRefArray * clusterListESD = new TRefArray();
      event->GetEMCALClusters(clusterListESD); 
      TClonesArray * clusterListAOD = new TClonesArray();
      if(AODEvent()){clusterListAOD = dynamic_cast<TClonesArray*> (AODEvent()->FindListObject("newEMCALClusters"));}

/*
  const AliESDVertex *vertex = event->GetPrimaryVertex();
    if (!vertex)   return;  
    if(! vertex->GetStatus()) return;
    if(TMath::Abs(vertex->GetZv())>10) return;

  Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);// | AliVEvent::kMBNoTRD);

  TString trigClasses = event->GetFiredTriggerClasses();

  if(trigClasses.Contains("CINT1") || trigClasses.Contains("CSH")  || trigClasses.Contains("CMU") || trigClasses.Contains("CBEAM") ) {
     printf("You are analyzing the REAL DATA!!! \n");
     if(!isSelected) return;

  }
*/  

  fhNEvents->Fill(0);

  // centrality
  AliCentrality *centrality =(AliCentrality *) event->GetCentrality();
   Float_t opt = (Float_t)centrality->GetCentralityPercentile("V0M");




  AliEMCALGeoUtils * fGeom   = new AliEMCALGeoUtils("EMCAL_FIRSTYEARV1","EMCAL");
    fGeom->GetEMCGeometry();

  //_____________ for EMCal  cell loop________________________________________________
  AliESDCaloCells *CellEMCAL   = event->GetEMCALCells();
   Int_t multCells = CellEMCAL ->GetNumberOfCells();
   fhCentVsCellMultEMCAll->Fill(multCells,opt);
   Int_t CellModuCount[4]={0,0,0,0};
  for(Int_t icell=0; icell<multCells; icell++){
    Int_t AbsID = CellEMCAL->GetCellNumber(icell); 
//    cout<<"EMCalAbsId =   "<<AbsID<<endl;  
    Int_t iSM1=-1, iTower1=-1, Iphi1=-1, Ieta1=-1;
    fGeom->GetCellIndex(AbsID, iSM1, iTower1, Iphi1, Ieta1 ); 
    Int_t iPhi1=-1, iEta1=-1;
    fGeom->GetCellPhiEtaIndexInSModule(iSM1, iTower1, Iphi1, Ieta1, iPhi1,iEta1);
    for(Int_t i=0; i<nModuleEMC; i++){
      if(i==iSM1){
        CellModuCount[i]++; 
        fhCentVsColuRowEMC[i]->Fill(iEta1,iPhi1,opt,1);// if hit on or not
        fhCentVsColuRowEnerEMC[i]->Fill(iEta1,iPhi1,opt,CellEMCAL->GetCellAmplitude(AbsID));// if hit on or not
      }
    }
  }//cell
    for(Int_t i=0; i<nModuleEMC; i++){
   fhCentVsCellMultEMC[i]->Fill(CellModuCount[i],opt); 
    }  

//__________for EMCal  cluster loop______________________________________________
  Int_t multCaloClust = 0;
  if(!AODEvent()) { multCaloClust = clusterListESD->GetEntriesFast(); 
  }else{            multCaloClust = clusterListAOD->GetEntriesFast(); }
  Double_t vtx0[3] = {0,0,0};
//  Double_t ClustEnergy = 0.;
  Int_t multEMCaClust =0;
  Int_t ClustModuCount[4] ={0,0,0,0};
  
  if(multCaloClust>1){  
  for(Int_t i1=0; i1<multCaloClust-1; i1++){
    
      multEMCaClust++ ; 
      TLorentzVector pemc1;
   if(!AODEvent()){ ((AliESDCaloCluster*) clusterListESD->At(i1))->GetMomentum(pemc1,vtx0);
   }else{           ((AliAODCaloCluster*) clusterListAOD->At(i1))->GetMomentum(pemc1,vtx0);}
// cluster showershape
    Double_t qua1;
    if(!AODEvent()){
      AliESDCaloCluster * ESDcluster = (AliESDCaloCluster*) clusterListESD->At(i1);
      Double_t sigmaMax1 = AliAnalysisTaskEMCALPi0PbPb::GetSigmaMax(ESDcluster, event);
      fhPtEnSg->Fill(pemc1.Pt(), pemc1.E(), sigmaMax1*pemc1.E());
      for(Int_t j=0; j<ESDcluster->GetNCells();j++){
       qua1<event->GetEMCALCells()->GetCellAmplitude(ESDcluster->GetCellAbsId(j))/ESDcluster->E() ? qua1 = event->GetEMCALCells()->GetCellAmplitude(ESDcluster->GetCellAbsId(j))/ESDcluster->E() : qua1 = qua1;
      }
    }else{
      AliAODCaloCluster * AODcluster = (AliAODCaloCluster*) clusterListAOD->At(i1);
      Double_t sigmaMax1 = GetSigmaMax(AODcluster, AODEvent());
      fhPtEnSg->Fill(pemc1.Pt(), pemc1.E(), sigmaMax1*pemc1.E());
      for(Int_t j=0; j<AODcluster->GetNCells();j++){
       qua1<AODEvent()->GetEMCALCells()->GetCellAmplitude(AODcluster->GetCellAbsId(j))/AODcluster->E() ? qua1 = AODEvent()->GetEMCALCells()->GetCellAmplitude(AODcluster->GetCellAbsId(j))/AODcluster->E() : qua1 = qua1;
      }
    }

    Int_t cellAbsId1=-1;
    fGeom->GetAbsCellIdFromEtaPhi(pemc1.Eta(), pemc1.Phi(), cellAbsId1);
    Int_t iSM1=-1, iTower1=-1, Iphi1=-1, Ieta1=-1;
    fGeom->GetCellIndex(cellAbsId1, iSM1, iTower1, Iphi1, Ieta1 );
      for(Int_t i=0; i<nModuleEMC; i++){
        if(i==iSM1){
          ClustModuCount[i]++;
        }
      }
  for(Int_t i2=i1; i2<multCaloClust; i2++){
   TLorentzVector pemc2;
   if(!AODEvent()){ ((AliESDCaloCluster*)clusterListESD->At(i2))->GetMomentum(pemc2,vtx0);
   }else{           ((AliAODCaloCluster*)clusterListAOD->At(i2))->GetMomentum(pemc2,vtx0);}
// cluster showershape
    Double_t qua2;
    if(!AODEvent()){
      AliESDCaloCluster * ESDcluster = (AliESDCaloCluster*) clusterListESD->At(i2);
      for(Int_t j=0; j<ESDcluster->GetNCells();j++){
       qua2<event->GetEMCALCells()->GetCellAmplitude(ESDcluster->GetCellAbsId(j))/ESDcluster->E() ? qua2 = event->GetEMCALCells()->GetCellAmplitude(ESDcluster->GetCellAbsId(j))/ESDcluster->E() : qua2 = qua2;
      }
    }else{
      AliAODCaloCluster * AODcluster = (AliAODCaloCluster*) clusterListAOD->At(i2);
      for(Int_t j=0; j<AODcluster->GetNCells();j++){
       qua2<AODEvent()->GetEMCALCells()->GetCellAmplitude(AODcluster->GetCellAbsId(j))/AODcluster->E() ? qua2 = AODEvent()->GetEMCALCells()->GetCellAmplitude(AODcluster->GetCellAbsId(j))/AODcluster->E() : qua2 = qua2;
      }
    }

   TLorentzVector pemc12  = pemc1  + pemc2;
   Double_t asy  = TMath::Abs((pemc1.E()-pemc2.E())/(pemc1.E()+pemc2.E()));
   Double_t qua = TMath::Max(qua1,qua2);
//   fhEMCALAsymPt->Fill(asy, p12.Pt());
   fhCentVsInMVsPtEMC->Fill(pemc12.M(), pemc12.Pt(), opt);
   fhCentVsInMVsPhiEMC->Fill(pemc12.M(), pemc12.Phi(), opt); 
   fhAsyVsPt->Fill(asy, pemc12.Pt());
   fhMggPtEr->Fill(pemc12.M(), pemc12.Pt(), qua);

  }//loop second cluster
   }// EMCal first cluster loop
  } //multCaloCluster>1

     for(Int_t i=0; i<nModuleEMC; i++){ // fill Module
          fhCentVsClusMultEMC[i]->Fill(ClustModuCount[i],opt);
      }
    fhCentVsClusMultEMCAll->Fill(multEMCaClust,opt);

   
// Post output data.
  PostData(1, fOutputDataESD);
}      

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
}
//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPi0PbPb::GetSigmaMax(AliESDCaloCluster * cluster, AliESDEvent * event){
 Double_t sigmaMax;
 Double_t Ec = cluster->E();
 Double_t Xc = 0 ;
 Double_t Yc = 0 ;
 Double_t Sxx = 0 ;
 Double_t Sxy = 0 ;
 Double_t Syy = 0 ;
 AliEMCALGeoUtils * fGeom   = new AliEMCALGeoUtils("EMCAL_FIRSTYEARV1","EMCAL");
 fGeom->GetEMCGeometry();
 for(Int_t j=0; j<cluster->GetNCells();j++){
   TVector3 pos;
   fGeom->GetGlobal(cluster->GetCellAbsId(j),pos);
   Xc = Xc + event->GetEMCALCells()->GetCellAmplitude(cluster->GetCellAbsId(j))*pos.X()/Ec; 
   Yc = Yc + event->GetEMCALCells()->GetCellAmplitude(cluster->GetCellAbsId(j))*pos.Y()/Ec;    
 }
 for(Int_t j=0; j<cluster->GetNCells();j++){
   TVector3 pos;
   fGeom->GetGlobal(cluster->GetCellAbsId(j),pos);
   Sxx = Sxx + event->GetEMCALCells()->GetCellAmplitude(cluster->GetCellAbsId(j))*(pos.X()-Xc)*(pos.X()-Xc)/Ec; 
   Syy = Syy + event->GetEMCALCells()->GetCellAmplitude(cluster->GetCellAbsId(j))*(pos.Y()-Yc)*(pos.Y()-Yc)/Ec; 
   Sxy = Sxy + event->GetEMCALCells()->GetCellAmplitude(cluster->GetCellAbsId(j))*(pos.X()-Xc)*(pos.Y()-Yc)/Ec; 
 } 
 sigmaMax = (Sxx + Syy + TMath::Sqrt((Sxx-Syy)*(Sxx-Syy)+4.0*Sxy*Sxy))/2.0;
 sigmaMax = TMath::Sqrt(sigmaMax); 
 return  sigmaMax;
}
//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPi0PbPb::GetSigmaMax(AliAODCaloCluster * cluster, AliAODEvent * event){
 Double_t sigmaMax;
 Double_t Ec = cluster->E();
 Double_t Xc = 0 ;
 Double_t Yc = 0 ;
 Double_t Sxx = 0 ;
 Double_t Sxy = 0 ;
 Double_t Syy = 0 ;
 AliEMCALGeoUtils * fGeom   = new AliEMCALGeoUtils("EMCAL_FIRSTYEARV1","EMCAL");
 fGeom->GetEMCGeometry();
 for(Int_t j=0; j<cluster->GetNCells();j++){
   TVector3 pos;
   fGeom->GetGlobal(cluster->GetCellAbsId(j),pos);
   Xc = Xc + event->GetEMCALCells()->GetCellAmplitude(cluster->GetCellAbsId(j))*pos.X()/Ec; 
   Yc = Yc + event->GetEMCALCells()->GetCellAmplitude(cluster->GetCellAbsId(j))*pos.Y()/Ec;    
 }
 for(Int_t j=0; j<cluster->GetNCells();j++){
   TVector3 pos;
   fGeom->GetGlobal(cluster->GetCellAbsId(j),pos);
   Sxx = Sxx + event->GetEMCALCells()->GetCellAmplitude(cluster->GetCellAbsId(j))*(pos.X()-Xc)*(pos.X()-Xc)/Ec; 
   Syy = Syy + event->GetEMCALCells()->GetCellAmplitude(cluster->GetCellAbsId(j))*(pos.Y()-Yc)*(pos.Y()-Yc)/Ec; 
   Sxy = Sxy + event->GetEMCALCells()->GetCellAmplitude(cluster->GetCellAbsId(j))*(pos.X()-Xc)*(pos.Y()-Yc)/Ec; 
 } 
 sigmaMax = (Sxx + Syy + TMath::Sqrt((Sxx-Syy)*(Sxx-Syy)+4.0*Sxy*Sxy))/2.0;
 sigmaMax = TMath::Sqrt(sigmaMax); 
 return  sigmaMax;

}
