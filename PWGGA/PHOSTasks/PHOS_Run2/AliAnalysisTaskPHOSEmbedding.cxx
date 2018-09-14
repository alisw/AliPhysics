#include <stdio.h>
#include <iostream>

#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TRandom3.h"
#include "TGrid.h"
#include "TROOT.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TList.h"
#include "THashList.h"

#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliLog.h"

#include "AliESDtrackCuts.h"
#include "AliESDHeader.h"
#include "AliESDEvent.h"
#include "AliESDVZERO.h"
#include "AliESDTZERO.h"
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
#include "AliPDG.h"

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

#include "TGeoManager.h"
#include "TGeoGlobalMagField.h"
#include "AliGeomManager.h"
#include "AliGRPObject.h"
#include "AliMagF.h"

#include "AliPHOSGeometry.h"
#include "AliPHOSDigit.h"
#include "AliPHOSReconstructor.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSCalibData.h"
#include "AliCDBManager.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "AliCDBMetaData.h"
#include "AliPHOSAodCluster.h"
#include "AliPHOSEsdCluster.h"
#include "AliOADBContainer.h"

#include "AliAnalysisTaskPHOSEmbedding.h"

using namespace std;

// Author: Daiki Sekihata (Hiroshima University)

ClassImp(AliAnalysisTaskPHOSEmbedding)

//____________________________________________________________________________________________________________________________________
AliAnalysisTaskPHOSEmbedding::AliAnalysisTaskPHOSEmbedding(const char *name):
  AliAnalysisTaskSE(name),
  //AliAnalysisTaskESDfilter(name),
  fOutputContainer(0x0),
  fHistoFileID(0x0),
  fHistoEventID(0x0),
  fHistoStatus(0x0),
  fHistoPt(0x0),
  fHistoEtaPhi(0x0),
  fHistoEtaPt(0x0),
  fParticle(""),
  fEvent(0x0),
  fRandom3(0x0),
  fAODPathArray(0x0),
  fAODPath(""),
  fAODInput(0x0),
  fAODTree(0x0),
  fAODEvent(0x0),
  fEventCounter(0),
  fStartID(0),
  fNeventMC(0),
  fCellsPHOS(0x0),
  fDigitsTree(0x0),
  fClustersTree(0x0),
  fDigitsArr(0x0),
  fPHOSReconstructor(0x0),
  fClusterizer(0x0),
  fInitialized(kFALSE),
  fRunNumber(0),
  fESDEvent(0x0),
  fMCArray(0x0),
  fEmbeddedClusterArray(0x0),
  fEmbeddedCells(0x0),
  fSignalECorrection(1.),
  fSignalCalibData(0x0),
  fUserNonLin(0x0),
  fApplyZS(kFALSE),
  fZSThreshold(0.020)
{
  // Constructor

  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");
  cdb->SetRun(246982);//dummy run number//necessary run number is set in Init()
  fUserNonLin = new TF1("fUserNonLin","1.",0,100);

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, THashList::Class());

}
//____________________________________________________________________________________________________________________________________
AliAnalysisTaskPHOSEmbedding::~AliAnalysisTaskPHOSEmbedding()
{
  if(fRandom3){
    delete fRandom3;
    fRandom3 = 0x0;
  }

  if(fUserNonLin){
    delete fUserNonLin;
    fUserNonLin = 0x0;
  }

  if(fPHOSReconstructor){
    delete fPHOSReconstructor;
    fPHOSReconstructor = 0x0;
  }

  if(fMCArray){
    delete fMCArray;
    fMCArray = 0x0;
  }

  if(fEmbeddedClusterArray){
    delete fEmbeddedClusterArray;
    fEmbeddedClusterArray = 0x0;
  }

  if(fEmbeddedCells){
    delete fEmbeddedCells;
    fEmbeddedCells = 0x0;
  }

  if(fAODEvent){
    //remove previous input MC event
    delete fAODEvent;
    fAODEvent = 0x0;
  }

  if(fHistoFileID){
    delete fHistoFileID;
    fHistoFileID = 0x0;
  }

  if(fHistoEventID){
    delete fHistoEventID;
    fHistoEventID = 0x0;
  }

  if(fHistoStatus){
    delete fHistoStatus;
    fHistoStatus = 0x0;
  }

  if(fHistoPt){
    delete fHistoPt;
    fHistoPt = 0x0;
  }

  if(fHistoEtaPhi){
    delete fHistoEtaPhi;
    fHistoEtaPhi = 0x0;
  }

  if(fHistoEtaPt){
    delete fHistoEtaPt;
    fHistoEtaPt = 0x0;
  }

}
//____________________________________________________________________________________________________________________________________
void AliAnalysisTaskPHOSEmbedding::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);

  const Int_t Nfile = 500;
  const Int_t Nev = 30e+3;
  fHistoFileID  = new TH1F(Form("hEventFileID_%s" ,fParticle.Data()),Form("file index in MC str array %s;file ID;Number of events",fParticle.Data()),Nfile+1,-0.5,Nfile+0.5);
  fHistoEventID = new TH1F(Form("hEventEventID_%s",fParticle.Data()),Form("event index in MC AOD %s;event ID;Number of events"    ,fParticle.Data()),Nev+1  ,-0.5,Nev+0.5);
  fOutputContainer->Add(fHistoFileID);
  fOutputContainer->Add(fHistoEventID);

  fHistoStatus = new TH1F(Form("hIOStatus_%s",fParticle.Data()),Form("byte count read from MC AOD %s;kbyte;Number of events",fParticle.Data()),502,-2,500);//unit of kilo byte
  fOutputContainer->Add(fHistoStatus);

  fHistoPt     = new TH1F(Form("hMotherPt_%s"    ,fParticle.Data()),Form("p_{T} of %s;p_{T} (GeV/c);Number of counts",fParticle.Data()),500,0,50);
  fHistoEtaPhi = new TH2F(Form("hMotherEtaPhi_%s",fParticle.Data()),Form("y vs. #phi of %s;#phi (radian);rapidity"   ,fParticle.Data()),60,0,TMath::TwoPi(),200,-1,1);
  fHistoEtaPt  = new TH2F(Form("hMotherEtaPt_%s" ,fParticle.Data()),Form("y vs p_{T} of %s;rapidity;p_{T} (GeV/c)"   ,fParticle.Data()),200,-1,1,500,0,50);

  fHistoPt    ->Sumw2();
  fHistoEtaPhi->Sumw2();
  fHistoEtaPt ->Sumw2();

  fOutputContainer->Add(fHistoPt    );
  fOutputContainer->Add(fHistoEtaPhi);
  fOutputContainer->Add(fHistoEtaPt );


  fMCArray = new TClonesArray("AliAODMCParticle",50);

  fEmbeddedClusterArray = new TClonesArray("AliAODCaloCluster",200);
  //fEmbeddedClusterArray->SetName(Form("Embedded%sCaloClusters",fParticle.Data()));

  fEmbeddedCells = new AliAODCaloCells();
  //fEmbeddedCells->SetName(Form("Embedded%sPHOScells",fParticle.Data()));
  
  PostData(1,fOutputContainer);

}
//____________________________________________________________________________________________________________________________________
void AliAnalysisTaskPHOSEmbedding::UserExec(Option_t *option) 
{
  // Main loop
  // Called for each event
  fEvent = dynamic_cast<AliVEvent*>(InputEvent());
  if(!fEvent){
    AliError("event is not available.");
    return;
  }

  if(fRunNumber != fEvent->GetRunNumber()){
    fRunNumber = fEvent->GetRunNumber();
    Init();
  }

  if(fAODInput){
    AliInfo("fAODInput exists. No need to open external MC file.");
  }
  else{
    AliInfo("fAODInput does not exist. Get fAODInput from external MC file.");
    Bool_t IsFileOK = OpenAODFile();//fEventCounter is reset here.
    if(!IsFileOK){
      AliError(Form("external file (%s) could not be opened. return.",fAODPath.Data()));
      return;
    }
  }

  Long64_t evID = (fStartID + (Long64_t)fEventCounter) % (Long64_t)fNeventMC;
  fHistoEventID->Fill(evID);
  //Long64_t evID = Long64_t(fNeventMC * fRandom3->Rndm());
  //fHistoEventID->Fill(evID);
  Int_t status = fAODTree->GetEvent(evID,0);//second argument is for "getall" in TBranch. //default is 0 (no getall).

  fHistoStatus->Fill((Double_t)status/1.e+3);
  if(status <= 0){
    //status is data size in unit of byte read from branch.
    //if this is less than or equal to 0, I/O problem.
    //https://root.cern.ch/root/html534/TTree.html#TTree:GetEntry
    AliInfo("There is I/O problem in external MC AOD file.");
    return;
  }

  AliInfo(Form("fStartID = %lld , fEventCounter = %d , Extract event ID %lld , %d byte read-out from %s.",fStartID,fEventCounter,evID,status,fAODPath.Data()));
  ConvertAODtoESD();

  //start embedding
  CopyRecalibrateDigits();

  if(fDigitsArr){
    delete fDigitsArr;
    fDigitsArr = 0x0;
  }
  fDigitsArr = new TClonesArray("AliPHOSDigit",4*56*64);
  fDigitsArr->Clear();

  gROOT->cd(); //make sure that the digits and RecPoints Trees are memory resident
  if(fDigitsTree){
    delete fDigitsTree;
    fDigitsTree = 0x0;
  }
  fDigitsTree = new TTree("digitstree","digitstree");
  fDigitsTree->Branch("PHOS","TClonesArray", &fDigitsArr, 32000);

  if(fClustersTree){
    delete fClustersTree;
    fClustersTree = 0x0;
  }
  fClustersTree = new TTree("clustertree","clustertree");

  //Remember number of Clusters before we added new ones...
  Int_t fNCaloClustersOld = fEvent->GetNumberOfCaloClusters();

  Int_t nPHOSBefore=0;
  Int_t nCPVBefore=0;
  Int_t nEMCALBefore=0;

  for(Int_t iClust=0; iClust<fEvent->GetNumberOfCaloClusters(); iClust++){
    AliVCluster *cluster = (AliVCluster*)fEvent->GetCaloCluster(iClust);

    if(cluster->GetType()      == AliVCluster::kPHOSNeutral) nPHOSBefore++;
    else if(cluster->GetType() == AliVCluster::kPHOSCharged) nCPVBefore++;
    else                                                     nEMCALBefore++;
  }

  AliInfo(Form("Before embedding: Nall = %d , nPHOS = %d , nCPV = %d , nEMCAL = %d.",fNCaloClustersOld, nPHOSBefore, nCPVBefore, nEMCALBefore));

  AliInfo(Form("Ncluster in sigle simulation = %d",fAODEvent->GetNumberOfCaloClusters()));
  AliAODCaloCells *cellsEmb = (AliAODCaloCells*)fAODEvent->GetPHOSCells();
  Int_t NcellEmb = cellsEmb->GetNumberOfCells();
  AliInfo(Form("Before embedding: N PHOS cells in simulation = %d",NcellEmb));

  //-------------------------------------------------------------------------------------
  //Transform CaloCells into Digits which can be used for standard reconstruction
  //Add signal digits to the event
  //-------------------------------------------------------------------------------------

  //First copy data digits
  Int_t ndigit=0 ;
  for (Short_t icell = 0; icell < fCellsPHOS->GetNumberOfCells(); icell++) {
    Short_t id=0;
    Double_t time=0., amp=0. ;
    Int_t mclabel;
    Double_t efrac =0. ;
    if(fCellsPHOS->GetCell(icell, id, amp, time,mclabel, efrac) != kTRUE) break;

    Int_t idLong=id ;
    if(id<0)idLong= -id+56*64*5 ; //CPV digits
    new((*fDigitsArr)[ndigit]) AliPHOSDigit(-1,idLong,float(amp),float(time),ndigit);
    ndigit++;
  }

  //Add Digits from Signal
  TClonesArray sdigits("AliPHOSDigit",56*64*5) ;
  Int_t isdigit=0 ;
  if(fAODEvent){
    AliAODCaloCells* cellsS = fAODEvent->GetPHOSCells();
    Int_t cellLabels[1000]={0} ;       //1000 should be enough for simulated
    Int_t cellSecondLabels[1000]={0} ; //low-statistics event.
    for(Int_t i=0;i<cellsS->GetNumberOfCells();i++){
      cellLabels[i]=-1 ;
      cellSecondLabels[i]=-1;
    }
    //------------------------------------------------------------------------------------
    //Ancestry information
    //Celect digits contributing to fAODEvent clusters and add primary information
    //(it is not stored in CaloCells)
    //------------------------------------------------------------------------------------
    //sdigits.Expand(cellsS->GetNumberOfCells());

    for(Int_t i=0; i<fAODEvent->GetNumberOfCaloClusters(); i++) {
      //cluster from embedded signal
      AliAODCluster *clus = (AliAODCaloCluster*)fAODEvent->GetCaloCluster(i);

      if(!clus->IsPHOS()) continue;

      AliInfo(Form("Ecluster in M.C. = %e GeV",clus->E()));

      Int_t label = clus->GetLabel();
      Int_t label2 = -1 ;
      if (clus->GetNLabels()>=2) label2 = clus->GetLabelAt(1) ;

      UShort_t * index    = clus->GetCellsAbsId() ;
      for(Int_t ic=0; ic < clus->GetNCells(); ic++ ){
        for (Int_t icell = 0; icell < cellsS->GetNumberOfCells(); icell++){
          Short_t cellNumber;
          Double_t cellAmplitude=0., cellTime=0. ;
          Int_t mclabel;
          Double_t efrac =0. ;
          cellsS->GetCell(icell, cellNumber, cellAmplitude, cellTime,mclabel,efrac) ;
          Int_t longCellNumber=cellNumber ;
          if(cellNumber<0)longCellNumber= -cellNumber+56*64*5 ; //CPV digits
          if(longCellNumber==index[ic]){
            cellLabels[icell]=label;
            cellSecondLabels[icell]=label2;
            break ;
          }
        }
      }
    }

    for (Int_t icell = 0; icell < cellsS->GetNumberOfCells(); icell++) {
      Short_t cellNumber;
      Double_t cellAmplitude=0., cellTime=0. ;
      Int_t mclabel;
      Double_t efrac =0. ;
      if (cellsS->GetCell(icell, cellNumber, cellAmplitude, cellTime,mclabel,efrac) != kTRUE) break;
      //Add only digits related to the cluster, no noisy digits...
      if(cellLabels[icell]==-1){
        continue ;
      }

      cellAmplitude = DecalibrateSignal(cellAmplitude,cellNumber);//decalibration from Dmitri, Daiki added this line on 21.December.2016

      if(fApplyZS && (cellAmplitude < fZSThreshold)){
        AliInfo(Form("cellAmplitude %e GeV is less than ZS threshold (%e GeV). This cell will not be merged.",cellAmplitude,fZSThreshold));
        continue;
      }

      Int_t longCellNumber=cellNumber ;
      if(cellNumber<0)longCellNumber= -cellNumber+56*64*5 ; //CPV digits
      new(sdigits[isdigit]) AliPHOSDigit(cellLabels[icell],longCellNumber,float(cellAmplitude),float(cellTime),isdigit);
      isdigit++;
    }
  }


  //Merge digits
  Int_t icurrent = 0 ; //index of the last used digit in underlying event
  fDigitsArr->Expand(ndigit+isdigit) ;
  for(Int_t i=0; i<isdigit;i++){
    AliPHOSDigit * sdigit=static_cast<AliPHOSDigit*>(sdigits.At(i)) ;
    Bool_t added=kFALSE ;
    for(Int_t id=icurrent;id<ndigit;id++){
      AliPHOSDigit * digit=static_cast<AliPHOSDigit*>(fDigitsArr->At(id)) ;
      if(sdigit->GetId() ==  digit->GetId() ){
        *digit += *sdigit ;  //add energies
        icurrent=id+1 ;
        added=kTRUE ;
        break ; //no more digits with same ID in the list
      }
      if(sdigit->GetId() < digit->GetId() ){
        icurrent=id ;
        break ; //no more digits with same ID in the list
      }
    }
    if(!added){
      new((*fDigitsArr)[ndigit]) AliPHOSDigit(*sdigit) ;
      ndigit++ ;
    }
  }

  sdigits.Clear();

  //Change Amp back from Energy to ADC counts
  //Note that Reconstructor uses best ("new") calibration
  for(Int_t i=0; i<ndigit;i++){
    AliPHOSDigit * digit=static_cast<AliPHOSDigit*>(fDigitsArr->At(i)) ;
    Float_t calib =fPHOSReconstructor->Calibrate(1.,digit->GetId()) ;
    if(calib>0.)
      digit->SetEnergy(digit->GetEnergy()/calib) ;
  }

  fDigitsArr->Sort() ;
  for (Int_t i = 0 ; i < ndigit ; i++) {
    AliPHOSDigit *digit = static_cast<AliPHOSDigit*>( fDigitsArr->At(i) ) ;
    digit->SetIndexInList(i) ;
  }
  fDigitsTree->Fill();

  //clusterize and make tracking//no track matching here!
  fPHOSReconstructor->Reconstruct(fDigitsTree,fClustersTree) ;

  //Note that the current ESDEvent will be modified!
  fPHOSReconstructor->FillESD(fDigitsTree, fClustersTree, fESDEvent);
  //this modified event should go to ConvertEmbeddedClusters

  Int_t nPHOSAfter=0;
  Int_t nCPVAfter=0;
  Int_t nEMCALAfter=0;
  
  for(Int_t iClust=0; iClust<fESDEvent->GetNumberOfCaloClusters(); iClust++) {
    AliESDCaloCluster * cluster = (AliESDCaloCluster*)fESDEvent->GetCaloCluster(iClust);

    if(cluster->IsPHOS()){
      if(cluster->GetType() == AliVCluster::kPHOSNeutral) nPHOSAfter++;
      if(cluster->GetType() == AliVCluster::kPHOSCharged) nCPVAfter++;
    }
    else{
      nEMCALAfter++;
    }

  }

  AliInfo(Form("After embedding: Nall = %d , nPHOSAfter = %d , nCPVAfter = %d , nEMCALAfter = %d.", fESDEvent->GetNumberOfCaloClusters(), nPHOSAfter, nCPVAfter, nEMCALAfter));
  const Int_t Ncell = fESDEvent->GetPHOSCells()->GetNumberOfCells();
  AliInfo(Form("%d PHOS cells after embedding.",Ncell));

  ConvertESDtoAOD();

  //convert AliAODMCParticle
  TClonesArray *mcarray_org = dynamic_cast<TClonesArray*>(fAODEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  const Int_t Npart = mcarray_org->GetEntries();

  fMCArray->Clear();
  //fMCArray->Expand(Npart);

  for(Int_t i=0;i<Npart;i++){
    AliAODMCParticle* aodpart = (AliAODMCParticle*)mcarray_org->At(i);
    //printf("PDG = %d , X = %e , Y = %e , Z = %e\n",aodpart->GetPdgCode(), aodpart->Xv(),aodpart->Yv(),aodpart->Zv());
    new ((*fMCArray)[i]) AliAODMCParticle(*aodpart);
  }

  const Int_t Np = fMCArray->GetEntries();
  for(Int_t i=0;i<Np;i++){
    AliAODMCParticle *p = (AliAODMCParticle*)fMCArray->At(i);
    Double_t pT  = p->Pt();
    Double_t y   = p->Y();
    Double_t phi = p->Phi();

    Double_t xv = p->Xv();
    Double_t yv = p->Yv();
    Double_t R  = TMath::Sqrt(xv*xv + yv*yv);
    Int_t pdg   = p->GetPdgCode();

    if(R > 1.) continue;

    if(fParticle.Contains("Pi0",TString::kIgnoreCase) && pdg == 111){
      fHistoPt->Fill(pT);
      fHistoEtaPhi->Fill(phi,y);
      fHistoEtaPt->Fill(y,pT);
    }
    else if(fParticle.Contains("Eta",TString::kIgnoreCase) && pdg == 221){
      fHistoPt->Fill(pT);
      fHistoEtaPhi->Fill(phi,y);
      fHistoEtaPt->Fill(y,pT);
    }
    else if(fParticle.Contains("Gamma",TString::kIgnoreCase) && pdg == 11){
      fHistoPt->Fill(pT);
      fHistoEtaPhi->Fill(phi,y);
      fHistoEtaPt->Fill(y,pT);
    }

  }//end of MC particles loop



  AliVEvent* input = InputEvent();

  const TString brname = Form("%s_%s",AliAODMCParticle::StdBranchName(),fParticle.Data());
  TObject* outMC = input->FindListObject(brname);
  AliInfo(Form("MC array name = %s , entry = %d.",brname.Data(),Npart));
  if(!outMC){
    AliInfo(Form("%s object is not found in input event. Add it.",brname.Data()));
    fMCArray->SetName(brname);
    input->AddObject(fMCArray);
  }

  const TString brname_cluster = Form("Embedded%sCaloClusters",fParticle.Data());
  AliInfo(Form("name of embedded cluster array = %s",brname_cluster.Data()));
  TObject* outClusters = input->FindListObject(brname_cluster);
  if(!outClusters){
    AliInfo(Form("%s object is not found in input event. Add it.",brname_cluster.Data()));
    fEmbeddedClusterArray->SetName(brname_cluster);
    input->AddObject(fEmbeddedClusterArray);
  }
 
  const TString brname_cell = Form("Embedded%sPHOScells",fParticle.Data());
  AliInfo(Form("name of embedded cells = %s",brname_cell.Data()));
  TObject* outCells = input->FindListObject(brname_cell);
  if(!outCells){
    AliInfo(Form("%s object is not found in input event. Add it.",brname_cell.Data()));
    fEmbeddedCells->SetName(brname_cell);
    input->AddObject(fEmbeddedCells);
  }

  fEventCounter++;

  //PostData(1, fOutputContainer);
}
//____________________________________________________________________________________________________________________________________
void AliAnalysisTaskPHOSEmbedding::Terminate(Option_t *option) 
{
  //Called once at the end of the query
  //In principle, this function is not needed...
  if(fAODInput && fAODInput->IsOpen()) fAODInput->Close();

  AliInfo(Form("%s is done.",GetName()));

}
//____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskPHOSEmbedding::UserNotify() 
{
  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  TTree* tr = aodH->GetTree();
  TDirectory* dr = tr->GetDirectory();
  TString DataPath = dr->GetPath();

  TObjArray *tx = DataPath.Tokenize("/");
  //tx->Print();
  Int_t Nt = tx->GetEntries();
  //for (Int_t i = 0; i < Nt; i++) std::cout << "i = " << i << " , " << ((TObjString *)(tx->At(i)))->String() << std::endl;

  TString filestr = ((TObjString *)(tx->At(Nt-2)))->String();
  Int_t fileID = filestr.Atoi();
  TString runstr = ((TObjString *)(tx->At(5)))->String();
  Int_t runNum = runstr.Atoi();

  if(!fRandom3){
    UInt_t seed = UInt_t(runNum * 1e+4) + UInt_t(fileID);
    AliInfo(Form("seed is %u",seed));
    fRandom3 = new TRandom3(seed);
  }

  return kTRUE;
}
//____________________________________________________________________________________________________________________________________
Int_t AliAnalysisTaskPHOSEmbedding::SelectAODFile()
{
  Int_t Nfile = fAODPathArray->GetEntries();
  Int_t ientry = (Int_t)(Nfile * fRandom3->Rndm());

  TObjString *objStr = (TObjString*) fAODPathArray->At(ientry);
  if(!objStr){
    AliError("Could not get path of AOD MC file.");
    return -1;
  }

  fAODPath = objStr->GetString();
  AliInfo(Form("External MC file is %s.",fAODPath.Data()));
  fHistoFileID->Fill(ientry);

  return Nfile;
}
//____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskPHOSEmbedding::OpenAODFile()
{
  AliInfo("Open a new file.");

  if(!gGrid){
    AliInfo("Trying to connect to AliEn ...");
    TGrid::Connect("alien://");
  } 
  if(fAODInput && fAODInput->IsOpen()) fAODInput->Close();

  Int_t Nfile = SelectAODFile();//1st trial
  fAODInput = TFile::Open(fAODPath,"TIMEOUT=180");

  if(!fAODInput){
    AliInfo(Form("fAODInput (%s) is not properly opened. return kFALSE",fAODPath.Data()));
    return kFALSE;
  }

  AliInfo(Form("%d files are stored in fAODPathArray.",Nfile));

  fAODTree = (TTree*)fAODInput->Get("aodTree");
  if(!fAODTree){
    AliError(Form("aodTree from an external AOD MC input (%s) is not found.",fAODPath.Data()));
    return kFALSE;
  }

  if(fAODEvent){
    //remove previous input MC event
    delete fAODEvent;
    fAODEvent = 0x0;
  }
  fAODEvent = new AliAODEvent();
  fAODEvent->ReadFromTree(fAODTree);
  fNeventMC = fAODTree->GetEntries();
  fStartID = Long64_t(fNeventMC * fRandom3->Rndm());

  AliInfo(Form("%d events are stored in %s.",fNeventMC,fAODPath.Data()));
  AliInfo(Form("fStartID = %lld",fStartID));

  return kTRUE;
}
//____________________________________________________________________________________________________________________________________
void AliAnalysisTaskPHOSEmbedding::CopyRecalibrateDigits()
{
  //Recalibrate digits if there is better calibration ("new")
  //exists than one used in reconstruction ("ESD")

  AliAODEvent *event = dynamic_cast<AliAODEvent*>(fEvent);
  if(fCellsPHOS){
    delete fCellsPHOS;
    fCellsPHOS = 0x0;
  }
  fCellsPHOS = new AliAODCaloCells();
  fCellsPHOS->CreateContainer(event->GetPHOSCells()->GetNumberOfCells());
  AliInfo(Form("%d PHOS cells in real data.",event->GetPHOSCells()->GetNumberOfCells()));

  for (Short_t icell=0;icell<event->GetPHOSCells()->GetNumberOfCells();icell++){
    Short_t id=0;
    Double_t time=0., amp=0.;
    Int_t mclabel = -1;
    Double_t efrac =0.;
    if(event->GetPHOSCells()->GetCell(icell, id, amp, time, mclabel, efrac) != kTRUE) break;
    fCellsPHOS->SetCell(icell, id, amp, time, mclabel, efrac);
  }

}
//____________________________________________________________________________________________________________________________________
void AliAnalysisTaskPHOSEmbedding::Init()
{

  AliAODEvent *event = dynamic_cast<AliAODEvent*>(InputEvent());

  if(!event){
    AliError("Can not obtain InputEvent!") ;
    return;
  }

  Int_t runNum = event->GetRunNumber();
  AliCDBManager::Instance()->SetRun(runNum);
  //AliCDBManager::Instance()->SetDefaultStorage("raw://");

  if(fPHOSReconstructor){
    delete fPHOSReconstructor;
    fPHOSReconstructor = 0x0;
  }
  fPHOSReconstructor = new AliPHOSReconstructor("Run2");

  AliCDBPath path("PHOS","Calib","RecoParam");
  AliCDBEntry *entry = AliCDBManager::Instance()->Get(path.GetPath());
  if(!entry){
    AliError(Form("Can not get OCDB entry %s",path.GetPath().Data())) ;
    return ;
  }

  TObjArray* recoParamArray = (TObjArray*)entry->GetObject();
  AliPHOSRecoParam* recoParam = (AliPHOSRecoParam*)recoParamArray->At(2);
  fPHOSReconstructor->SetRecoParam(recoParam) ;

  InitMF() ;
  InitGeometry() ;

  //Get OADB calibration for this run for de-calibration of signal
  AliOADBContainer calibContainer("phosRecalibration");
  calibContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSCalibrations.root","phosRecalibration");
  TObjArray *recalib = (TObjArray*)calibContainer.GetObject(runNum,"PHOSRecalibration");
  if(!recalib){
    AliWarning(Form("Can not read calibrations for run %d, do not apply OADB de-calibration\n",runNum)) ;
  }
  else{
    const Int_t recoPass=1;
    fSignalCalibData = (AliPHOSCalibData*)recalib->At(recoPass-1) ;
  }

}
//____________________________________________________________________________________________________________________________________
void AliAnalysisTaskPHOSEmbedding::InitMF()
{
  //------------------------------------
  // Initialization of the Mag.Fiels from GRP entry
  // Copied from AliReconstruction::InitGRP()
  //------------------------------------
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  AliGRPObject * aGRPData=0 ;
  if (entry) {

    TMap* m = dynamic_cast<TMap*>(entry->GetObject());  // old GRP entry

    if (m) {
       AliInfo("Found a TMap in GRP/GRP/Data, converting it into an AliGRPObject");
       m->Print();
       aGRPData = new AliGRPObject();
       aGRPData->ReadValuesFromMap(m);
    }

    else {
       AliInfo("Found an AliGRPObject in GRP/GRP/Data, reading it");
       aGRPData = dynamic_cast<AliGRPObject*>(entry->GetObject());  // new GRP entry
       entry->SetOwner(0);
    }
  }

  if (!aGRPData) {
     AliError("No GRP entry found in OCDB!");
     return ;
  }
  //*** Dealing with the magnetic field map

 TString lhcState = aGRPData->GetLHCState();
  if (lhcState==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the LHC state ! Using UNKNOWN");
    lhcState = "UNKNOWN";
  }

  TString beamType = aGRPData->GetBeamType();
  if (beamType==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam type ! Using UNKNOWN");
    beamType = "UNKNOWN";
  }

  Float_t beamEnergy = aGRPData->GetBeamEnergy();
  if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
    beamEnergy = 0;
  }

  TString runType = aGRPData->GetRunType();
  if (runType==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the run type ! Using UNKNOWN");
    runType = "UNKNOWN";
  }

  if ( TGeoGlobalMagField::Instance()->IsLocked() ) {
    if (TGeoGlobalMagField::Instance()->GetField()->TestBit(AliMagF::kOverrideGRP)) {
      AliInfo("PHOSEmbedding: MF is locked - GRP information will be ignored !");
      AliInfo("Running with the externally locked B field !");
    }
    else {
      AliInfo("Destroying existing B field instance!");
      delete TGeoGlobalMagField::Instance();
    }
  }
  if ( !TGeoGlobalMagField::Instance()->IsLocked() ) {
    // Construct the field map out of the information retrieved from GRP.
    Bool_t ok = kTRUE;
    // L3
    Float_t l3Current = aGRPData->GetL3Current((AliGRPObject::Stats)0);
    if (l3Current == AliGRPObject::GetInvalidFloat()) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 current !");
      ok = kFALSE;
    }

    Char_t l3Polarity = aGRPData->GetL3Polarity();
    if (l3Polarity == AliGRPObject::GetInvalidChar()) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 polarity !");
      ok = kFALSE;
    }

    // Dipole
    Float_t diCurrent = aGRPData->GetDipoleCurrent((AliGRPObject::Stats)0);
    if (diCurrent == AliGRPObject::GetInvalidFloat()) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole current !");
      ok = kFALSE;
    }

    Char_t diPolarity = aGRPData->GetDipolePolarity();
    if (diPolarity == AliGRPObject::GetInvalidChar()) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole polarity !");
      ok = kFALSE;
    }

    // read special bits for the polarity convention and map type
    Int_t  polConvention = aGRPData->IsPolarityConventionLHC() ? AliMagF::kConvLHC : AliMagF::kConvDCS2008;
    Bool_t uniformB = aGRPData->IsUniformBMap();

    if (ok) {
      AliMagF* fld = AliMagF::CreateFieldMap(TMath::Abs(l3Current) * (l3Polarity ? -1:1),
                                             TMath::Abs(diCurrent) * (diPolarity ? -1:1),
                                             polConvention,uniformB,beamEnergy, beamType.Data());
      if (fld) {
        TGeoGlobalMagField::Instance()->SetField( fld );
        TGeoGlobalMagField::Instance()->Lock();
        AliInfo("Running with the B field constructed out of GRP !");
      }
      else AliFatal("Failed to create a B field map !");
    }
    else AliFatal("B field is neither set nor constructed from GRP ! Exitig...");
  }

}
//____________________________________________________________________________________________________________________________________
void AliAnalysisTaskPHOSEmbedding::InitGeometry()
{
  // Import ideal TGeo geometry and apply misalignment
  if (!gGeoManager) {
    AliGeomManager::LoadGeometry("geometry.root");
    if (!gGeoManager) {
      AliFatal("Can not load geometry");
    }
    if(!AliGeomManager::CheckSymNamesLUT("PHOS")) {
      AliFatal("CheckSymNamesLUT");
    }
  }
 
  TString detStr = "PHOS";
  TString loadAlObjsListOfDets = "PHOS";

  if(AliGeomManager::GetNalignable("GRP") != 0) loadAlObjsListOfDets.Prepend("GRP "); //add alignment objects for non-sensitive modules

  AliGeomManager::ApplyAlignObjsFromCDB(loadAlObjsListOfDets.Data());

  AliCDBManager::Instance()->UnloadFromCache("*/Align/*");
  AliCDBManager::Instance()->UnloadFromCache("GRP/Geometry/Data");

}
//____________________________________________________________________________________________________________________________________
void AliAnalysisTaskPHOSEmbedding::ConvertAODtoESD()
{
  //this ESD event is necessary for AliPHOSReconstructor::FillESD.

  if(fESDEvent){
    delete fESDEvent;
    fESDEvent = 0x0;
  }

  const AliVVertex *vVertex = fEvent->GetPrimaryVertex();
  Double_t vtx_data[3] = {-999,-999,-999};
  vtx_data[0] = vVertex->GetX();
  vtx_data[1] = vVertex->GetY();
  vtx_data[2] = vVertex->GetZ();
  AliInfo(Form("vetex position in AOD | X = %f cm , Y = %f cm , Z = %f cm.",vtx_data[0],vtx_data[1],vtx_data[2]));

  //create AliESDEvent and set vertex position from AOD.
  //only vertex position is necessary for AliPHOSPID/AliPHOSTrackSegmentsMaker
  //think again AliPHOSReconstructor::FillESD should be called or I should write them here manually.
  //at least, track maching is re-done in PHOSTender.

  fESDEvent = new AliESDEvent();
  fESDEvent->CreateStdContent();//this is important to set PHOSCells in AliESDEvent.

  AliESDVertex *vertex = new AliESDVertex();
  vertex->SetXYZ(vtx_data);
  fESDEvent->SetPrimaryVertexTracks(vertex);
  fESDEvent->SetPrimaryVertexSPD(vertex);
  fESDEvent->SetPrimaryVertexTPC(vertex);

  delete vertex;
  vertex = 0x0;

  const Int_t Ncell = fEvent->GetPHOSCells()->GetNumberOfCells();
  AliESDCaloCells *phoscells = dynamic_cast<AliESDCaloCells*>(fESDEvent->GetPHOSCells());
  phoscells->CreateContainer(Ncell);
  AliInfo(Form("%d PHOS cells in real data before embedding.",Ncell));

  for (Short_t icell=0;icell<Ncell;icell++){
    Short_t id=0;
    Double_t time=0., amp=0.;
    Int_t mclabel = -1;
    Double_t efrac =0.;
    if(fEvent->GetPHOSCells()->GetCell(icell, id, amp, time, mclabel, efrac) != kTRUE) break;
    phoscells->SetCell(icell, id, amp, time, mclabel, efrac);
  }

  const AliESDVertex *esdVertexBest = fESDEvent->GetPrimaryVertex();
  Double_t vtxBest[3] = {};
  vtxBest[0] = esdVertexBest->GetX();
  vtxBest[1] = esdVertexBest->GetY();
  vtxBest[2] = esdVertexBest->GetZ();
  AliInfo(Form("vetex position in ESD | X = %f cm , Y = %f cm , Z = %f cm.",vtxBest[0],vtxBest[1],vtxBest[2]));

}
//____________________________________________________________________________________________________________________________________
void AliAnalysisTaskPHOSEmbedding::ConvertESDtoAOD()
{
  fEmbeddedClusterArray->Clear();
  fEmbeddedCells->Clear("");

  // Access to the AOD container of clusters
  Int_t jClusters = 0;

  const Int_t Ncluster = fESDEvent->GetNumberOfCaloClusters();
  //fEmbeddedClusterArray->Expand(Ncluster);

  //PHOS + CPV clusters are stored.
  for (Int_t iClust=0; iClust<Ncluster; iClust++) {
    AliESDCaloCluster * cluster = (AliESDCaloCluster*)fESDEvent->GetCaloCluster(iClust);

    Int_t  id        = cluster->GetID();
    Int_t  nLabel    = cluster->GetNLabels();
    Int_t *labels    = cluster->GetLabels();

    Float_t energy = cluster->E();
    Float_t posF[3] = {0.};
    cluster->GetPosition(posF);

    AliAODCaloCluster *caloCluster = new((*fEmbeddedClusterArray)[jClusters++]) AliAODCaloCluster(id,
        nLabel,
        labels,
        energy,
        posF,
        NULL,
        cluster->GetType(),0);

    //Double_t dx=cluster->GetTrackDx() ;
    //Double_t dz=cluster->GetTrackDz() ;
    Float_t cpv=99. ; //No track matched by default

//no track matching
//    TArrayI * itracks = cluster->GetTracksMatched() ;
//    if(itracks->GetSize()>0){
//      Int_t iTr = itracks->At(0);
//      if(iTr>=0 && iTr<fESDEvent->GetNumberOfTracks()){
//        AliESDtrack *track = fESDEvent->GetTrack(iTr) ;
//        Double_t pt = track->Pt() ;
//        Short_t charge = track->Charge() ;
//        cpv=TestCPVRun2(dx, dz, pt,charge) ;
//      }
//    }

    caloCluster->SetCaloCluster(cluster->GetDistanceToBadChannel(),
        cluster->GetDispersion(),
        cluster->GetM20(), cluster->GetM02(),
        cpv,
        cluster->GetNExMax(),cluster->GetTOF()) ;

    caloCluster->SetPIDFromESD(cluster->GetPID());
    caloCluster->SetNCells(cluster->GetNCells());
    caloCluster->SetCellsAbsId(cluster->GetCellsAbsId());
    caloCluster->SetCellsAmplitudeFraction(cluster->GetCellsAmplitudeFraction());
  }

  //for(Int_t i=0;i<fEmbeddedClusterArray->GetEntriesFast();i++){
  //  AliAODCaloCluster *caloCluster =(AliAODCaloCluster *)fEmbeddedClusterArray->At(i) ;
  //  caloCluster->GetID();
  //}

  //Now cells
  if(fESDEvent->GetPHOSCells()){ // protection against missing ESD information
    AliESDCaloCells &fESDEventPHcells = *(fESDEvent->GetPHOSCells());
    Int_t nPHcell = fESDEventPHcells.GetNumberOfCells() ;
    fEmbeddedCells->CreateContainer(nPHcell);
    fEmbeddedCells->SetType(AliAODCaloCells::kPHOSCell);
    for (Int_t iCell = 0; iCell < nPHcell; iCell++) {
      //fEmbeddedCells->SetCell(iCell,fESDEventPHcells.GetCellNumber(iCell),fESDEventPHcells.GetAmplitude(iCell),fESDEventPHcells.GetTime(iCell));
      fEmbeddedCells->SetCell(iCell,fESDEventPHcells.GetCellNumber(iCell),fESDEventPHcells.GetAmplitude(iCell),fESDEventPHcells.GetTime(iCell),fESDEventPHcells.GetMCLabel(iCell),fESDEventPHcells.GetEFraction(iCell),fESDEventPHcells.GetHighGain(iCell));

    }
    fEmbeddedCells->Sort();
  }

  AliInfo(Form("entry of fEmbeddedClusterArray = %d (PHOS+CPV) and Ncell = %d (PHOS+CPV) after embedding.",fEmbeddedClusterArray->GetEntriesFast(),fEmbeddedCells->GetNumberOfCells()));

}
//____________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskPHOSEmbedding::DecalibrateSignal(Double_t cellAmplitude,Int_t cellNumber)
{
  //Apply de-calibration inverse to the calibration, stored in OADB

  AliInfo(Form("before correction : cellAmplitude = %e (GeV), fSignalECorrection = %e , nonlin = %e",cellAmplitude,fSignalECorrection,fUserNonLin->Eval(cellAmplitude)));
  //Apply overall energy correction and nonlin
  cellAmplitude *= fUserNonLin->Eval(cellAmplitude);
  cellAmplitude *= fSignalECorrection;

  if(!fSignalCalibData){
    AliInfo("fSignalCalibData is not applied.");
    return cellAmplitude;
  }

  Int_t relId[4];
  AliPHOSGeometry * phosgeom = AliPHOSGeometry::GetInstance();
  phosgeom->AbsToRelNumbering(cellNumber,relId);
  if(relId[1]!=0) //CPV
    return  cellAmplitude ;
  Int_t   module = relId[0];
  Int_t   column = relId[3];
  Int_t   row    = relId[2];
  Double_t c = fSignalCalibData->GetADCchannelEmc(module,column,row);
  AliInfo(Form("calibration co-efficient at M%d : X%d : Z%d in OADB is %e",module,row,column,c));

  if(c>0){
    cellAmplitude *= 1./c;
    AliInfo(Form("after correction : cellAmplitude = %e (GeV)",cellAmplitude));
    return  cellAmplitude;
  }
  else{
    //cellAmplitude *= 1./TMath::Abs(c);
    //AliInfo(Form("after correction : cellAmplitude = %e (GeV)",cellAmplitude));
    return  cellAmplitude;
  }
}
//____________________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________________

