#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TList.h>
#include <THashList.h>
#include <TClonesArray.h>  
#include <TGeoGlobalMagField.h>
#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TClonesArray.h"
#include "TEfficiency.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskGammaPHOSPP.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "AliAODEvent.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAODTrack.h"
#include "AliPHOSAodCluster.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliTriggerAnalysis.h"
#include "AliVVZERO.h"
#include "AliPHOSTriggerUtils.h"
#include "AliAODMCParticle.h"
#include "AliPHOSTenderSupply.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliInputEventHandler.h"
#include "AliTriggerAnalysis.h"
#include "AliOADBContainer.h"

#include "AliAnalysisManager.h"

// Analysis task to fill histograms with PHOS ESD clusters and cells
// Authors: Yuri Kharlov
// Date   : 28.05.2009

ClassImp(AliAnalysisTaskGammaPHOSPP)

//________________________________________________________________________
AliAnalysisTaskGammaPHOSPP::AliAnalysisTaskGammaPHOSPP(const char *name) : AliAnalysisTaskSE(name),
  fOutputContainer(0), fOutputContainer2(0),  
  fESDtrackCuts(0),
  fEvent(0), 
  fPHOSEvent(0), 
  fnCINT1B(0), fnCINT1A(0), fnCINT1C(0), fnCINT1E(0), 
  fEventVtxExist(0), fEventVtxZ10cm(0), fEventPileup(0), fEventV0AND(0), 
  fPHOSGeo(0),
  fInPHOS(0), 
  fMCArray(0),
  fPIDResponse(0x0), 
  fBCgap(525.e-9),
  fTOFcut(12.5e-09),
  fEventCounter(0), 
  fWeightFunction(new TF1("fWeightFunction", "((([0]+([1]*x))+([2]*(x*x)))/((1.+([3]*x))+([4]*(x*x))))+([5]*x)", 0., 100.)),
  fCurrFileName(0), 
  fCheckMCCrossSection(kFALSE),
  fh1Xsec(0), fh1Trials(0), fAvgTrials(-1), 
  fTriggerAnalysis(new AliTriggerAnalysis),
  NsigmaCPV(2.5),
  NsigmaDisp(2.5)

{

  fPidCuts.emplace_back("all", "no cuts");
  fPidCuts.emplace_back("cpv", "cpv");
  fPidCuts.emplace_back("disp", "disp");
  fPidCuts.emplace_back("both", "both");
// Constructor
  Int_t nBin=10 ;
  for (Int_t i=0; i<nBin; i++) {
    for (Int_t j=0; j<2; j++) {
      fPHOSEvents[i][j] = 0 ;
    }	 
   }
  
  // Absolute recalibration for LHC11a. Use SetRecalib(mod,recalib) to change it
  fRecalib[0] = 0.9942;
  fRecalib[1] = 0.9822;
  fRecalib[2] = 1.0072;

  // Output slots #0 write into a TH1 container
  DefineOutput(1,THashList::Class());
  DefineOutput(2,THashList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::UserCreateOutputObjects()
{
  // Create histograms, called once

  // AOD histograms
  //
  if (fOutputContainer != NULL) {
    delete fOutputContainer;
  }

  if (fOutputContainer2 != NULL) {
    delete fOutputContainer2;
  }


  fOutputContainer  = new THashList();
  fOutputContainer2 = new THashList();

  fOutputContainer ->SetOwner(kTRUE);
  fOutputContainer2->SetOwner(kTRUE);
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if (!man) {
      AliFatal("Could not find manager");
  }	

  AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*> (man->GetInputEventHandler());
  if (!inputHandler) {
      AliFatal("No input event handler");
  }	 

  fPIDResponse = dynamic_cast<AliPIDResponse *>(inputHandler->GetPIDResponse());
  if (!fPIDResponse) {
      AliFatal("PIDResponse object was not created"); // Escalated to fatal. This task is unusable without PID response.
  }	

  fPIDResponse -> SetUseTPCMultiplicityCorrection(kFALSE);
  fPIDResponse->SetCurrentMCEvent(MCEvent()); //!!
 
  fh1Xsec = new TH1F("hXsec","xsec from pyxsec.root", 1, 0, 1);
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
  fOutputContainer->Add(fh1Xsec);
  
  fh1Trials = new TH1F("hTrials","trials root file", 1, 0, 1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
  fOutputContainer->Add(fh1Trials);

  AddQAHistograms();
  AddMassHistograms();
  AddClusterHistograms();
  AddTrackHistograms();
  AddMCHistograms();

  PostData(1, fOutputContainer);
  PostData(2, fOutputContainer2);
}

/*
================================================================================
  // Main loop, called for each event

================================================================================
*/
void AliAnalysisTaskGammaPHOSPP::UserExec(Option_t *) 
{
  // Get event
  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fEvent) {
     Printf("ERROR: Could not retrieve event");
     return;
  }
  fLHCRunN = fEvent->GetRunNumber() < 224994 ? 1 : 2;

  // Initialize geometry
  if (fEventCounter == 0) {
    fPHOSGeo = AliPHOSGeometry::GetInstance() ; // use tender
  
    if (!fPHOSGeo) {
        AliInfo("PHOS geometry not initialized, initializing it for you");
        if(fLHCRunN == 1)
          fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ; // Run1 geometry
        else
          fPHOSGeo = AliPHOSGeometry::GetInstance("Run2") ;
    }
  }

  // Vertex 
  const AliAODVertex *aodVertex5 =    fEvent->GetPrimaryVertex();
  // const AliAODVertex *aodVertexSPD  = fEvent->GetPrimaryVertexSPD();
  fVtx5[0] = aodVertex5->GetX();
  fVtx5[1] = aodVertex5->GetY();
  fVtx5[2] = aodVertex5->GetZ();

  // Filter events
  Bool_t acceptEvent = kFALSE;
  acceptEvent = AcceptEvent(fEvent);
  if (!acceptEvent) return;
  
  //PHOS event
  if (fPHOSEvent) {
    fPHOSEvent->Clear() ;
  }  
  else
    fPHOSEvent = new TClonesArray("AliCaloPhoton", 50) ;

  // PID
  fPIDResponse -> SetUseTPCMultiplicityCorrection(kFALSE);
  AliPIDCombined *pidcomb=new AliPIDCombined();
  pidcomb->SetDefaultTPCPriors();
  pidcomb->SetSelectedSpecies(AliPID::kSPECIESC);
  pidcomb->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF|AliPIDResponse::kDetITS|AliPIDResponse::kDetTRD);

  // Event centrality 
  fEventCentrality = GetEventCentrality(fEvent);

  // Process MC
  fMCArray = (TClonesArray*)fEvent->FindListObject(AliAODMCParticle::StdBranchName());
  ProcessMC() ;
 
  // Notify 
  Notify();

  // Check PHOS and EMCAL clusters
  PHOSvsEMCALClusters();
  
  // PHOS cells 
  AnalyzeCells();

  // PHOS clusters *
  fInPHOS = 0 ;
  for (Int_t ic = 0; ic < fEvent->GetNumberOfCaloClusters(); ic++) {
    AliAODCaloCluster *clu1 = fEvent->GetCaloCluster(ic);
    SelectCluster(clu1);
  }

  // Photons 
  for (Int_t iph = 0; iph < fInPHOS; iph++) {
    AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSEvent->At(iph) ;
    FillOnePhotonHistograms(ph);
  }
 
  // Two photon histograms
  FillTwoPhotonHistograms();
  MixPhotons();

  // Events 
  FillHistogram("hEventCounter", 0.5);
  FillHistogram("hEventCounterCentrality", fEventCentrality + 0.5);
  fEventCounter++;
}

/*============================================================================*/
// Terminate
void AliAnalysisTaskGammaPHOSPP::FinishTaskOutput()
{
  GammaEfficiencies();
  CutEfficiencies();
}

//________________________________________________________________________
Bool_t AliAnalysisTaskGammaPHOSPP::AcceptEvent(AliAODEvent *aodEvent)
{
  FillHistogram("hSelEvents",0) ; // All events accepted by Physics Selection

  TString trigClasses = aodEvent->GetFiredTriggerClasses();

  if (trigClasses.Contains("FAST")  && !trigClasses.Contains("ALL")) {
    AliWarning(Form("Skip event with triggers %s", trigClasses.Data()));
    return kFALSE;
  }

  if (trigClasses.Contains("CINT1B")) fnCINT1B++;
  if (trigClasses.Contains("CINT1A")) fnCINT1A++;
  if (trigClasses.Contains("CINT1C")) fnCINT1C++;
  if (trigClasses.Contains("CINT1-E")) fnCINT1E++;

  // Event selection flags

  fEventVtxExist    = kFALSE;
  fEventVtxZ10cm    = kFALSE;
  fEventPileup      = kFALSE;
  fEventV0AND       = kFALSE;

  Int_t eventNumberInFile = aodEvent->GetEventNumberInFile();
  if (eventNumberInFile == 0) return kFALSE;

  if (aodEvent->GetPrimaryVertex()->GetNContributors() < 1 && !fMCArray) {
      fEventVtxExist    = kFALSE;
  }    
  else     
    fEventVtxExist    = kTRUE; 

  if (aodEvent->IsPileupFromSPD()) {
    fEventPileup = kTRUE;
  }  

  const AliAODVertex *aodVertex5 =    aodEvent->GetPrimaryVertex();
  const AliAODVertex *aodVertexSPD  = aodEvent->GetPrimaryVertexSPD();

  FillHistogram("hNvertexTracks", aodVertex5->GetNContributors());

  if (aodEvent->GetPrimaryVertex() && aodEvent->GetPrimaryVertex()->GetNContributors() > 0) {
    FillHistogram("hZvertex", aodVertex5->GetZ());
    if (TMath::Abs(aodVertex5->GetZ()) < 10.) {
      fEventVtxZ10cm = kTRUE;
    }
  }

  if (aodEvent->IsPileupFromSPD()) {
    fEventPileup = kTRUE;
    FillHistogram("hNPileupVtx", aodEvent->GetNumberOfPileupVerticesSPD());
    for (Int_t puVtx = 0;  puVtx < fEvent->GetNumberOfPileupVerticesSPD(); puVtx++) {
      Double_t dZpileup = aodVertexSPD->GetZ() - fEvent->GetPileupVertexSPD(puVtx)->GetZ();
      FillHistogram("hZPileupVtx", dZpileup);
    }
  }

  ULong64_t trigmask = aodEvent->GetTriggerMask();
  
  fEventV0AND = (trigmask & (1ull << (AliTriggerAnalysis::kV0AND-1)));

  // Fill event statistics for different selection criteria

  FillHistogram("hSelEvents",1) ;
  if (fEventVtxExist) 
    FillHistogram("hSelEvents",2) ;
  if (fEventVtxExist && fEventVtxZ10cm)
    FillHistogram("hSelEvents",3) ;
  if (fEventVtxExist && fEventVtxZ10cm && fEventV0AND)
    FillHistogram("hSelEvents",4) ;
  if (fEventVtxExist && fEventVtxZ10cm && fEventV0AND && fEventPileup)
    FillHistogram("hSelEvents",5) ;
  if (fEventPileup)
    FillHistogram("hSelEvents",6) ;
  if(fEventV0AND)
    FillHistogram("hSelEvents",7) ;
  if(fEventVtxZ10cm)
    FillHistogram("hSelEvents",8) ;  
        
  if (aodVertex5->GetNContributors() < 1  && !fMCArray)
     return kFALSE;
  if (aodVertexSPD->GetNContributors() < 1 && !fMCArray) 
     return kFALSE;
  if (!fEventVtxExist) 
     return kFALSE;
  if (!fEventVtxZ10cm) 
     return kFALSE;
  if (fEventPileup)
     return kFALSE;
  
  return kTRUE;
}
/*============================================================================*/
Int_t AliAnalysisTaskGammaPHOSPP::GetEventCentrality(AliAODEvent *event)
{
  //Calculate charged multiplicity   Bool_t   IsCPVOK(void)const {return fCpv;}
  
  for (Int_t i=0; i < event->GetNumberOfTracks(); i++) {
      AliAODTrack* track = dynamic_cast<AliAODTrack*>(event->GetTrack(i)) ;
      FillHistogram("hTrackCharge", track->Charge());
  }

  Int_t trackMult = event->GetNumberOfTracks() ;

  Float_t tV0A = event->GetVZEROData()->GetV0ATime();
  Float_t tV0C = event->GetVZEROData()->GetV0CTime();
  FillHistogram("hV0Atime", tV0A);
  FillHistogram("hV0Ctime", tV0C);
  FillHistogram("hV0AV0Ctime", tV0A, tV0C);  
  FillHistogram("hTrackMult", trackMult+0.5) ;

  Int_t centr = 0;
  if (trackMult <= 2)
    centr = 0 ;
  else 
    if (trackMult <= 5)
      centr = 1 ;
    else
      if (trackMult <= 9)
        centr = 2 ;
      else
        if (trackMult <= 14)
          centr = 3 ;
        else
          if (trackMult <= 22)
            centr = 4 ;
          else
            if (trackMult <= 35)
              centr = 5 ;
            else
              if (trackMult <= 50)
                centr = 6 ;
              else
                centr = 7 ;

  return centr;
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::FillHistogram(const char * key,Double_t x)const
{
  //FillHistogram
  TH1I * tmpI = dynamic_cast<TH1I*>((fOutputContainer->FindObject(key))?(fOutputContainer->FindObject(key)):(fOutputContainer2->FindObject(key)) ) ;
  if(tmpI) {
    tmpI->Fill(x) ;
    return ;
  }
  TH1F * tmpF = dynamic_cast<TH1F*>((fOutputContainer->FindObject(key))?(fOutputContainer->FindObject(key)):(fOutputContainer2->FindObject(key)) ) ;
  if(tmpF) {
    tmpF->Fill(x) ;
    return ;
  }
  TH1D * tmpD = dynamic_cast<TH1D*>((fOutputContainer->FindObject(key))?(fOutputContainer->FindObject(key)):(fOutputContainer2->FindObject(key)) ) ;
  if(tmpD) {
    tmpD->Fill(x) ;
    return ;
  }
  AliInfo(Form("can not find histogram <%s> (Fill I)",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::FillHistogram(const char * key,Double_t x,Double_t y)const
{
  //FillHistogram
  TObject * tmp = (fOutputContainer->FindObject(key))?(fOutputContainer->FindObject(key)):(fOutputContainer2->FindObject(key)) ;

  if (!tmp) {
    AliInfo(Form("can not find histogram <%s> (Fill II) ",key)) ;
    return ;
  }
  if (tmp->IsA() == TClass::GetClass("TH1F")) {
    ((TH1F*)tmp)->Fill(x,y) ;
    return ;
  }
  if (tmp->IsA() == TClass::GetClass("TH2F")) {
    ((TH2F*)tmp)->Fill(x,y) ;
    return ;
  }
  AliError(Form("Calling FillHistogram with 2 parameters for histo <%s> of type %s",key,tmp->IsA()->GetName())) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const
{
  //Fills 1D histograms with key
  TObject * tmp = (fOutputContainer->FindObject(key))?(fOutputContainer->FindObject(key)):(fOutputContainer2->FindObject(key)) ;

  if (!tmp) {
    AliInfo(Form("can not find histogram <%s> (Fill III) ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")) {
    ((TH2F*)tmp)->Fill(x,y,z) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH3F")) {
    ((TH3F*)tmp)->Fill(x,y,z) ;
    return ;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::ProcessMC()
{
  if (!fMCArray) return; 

  FillHistogram("hEventCounterMC",0.5);
  
  for (Int_t i = 0; i < fMCArray->GetEntriesFast(); i++) {
    AliAODMCParticle* particle =  (AliAODMCParticle*) fMCArray->At(i);
    
    FillHistogram("hMCTrackCharge", particle->Charge());

    if (TMath::Abs(particle->GetPdgCode()) == 111)
      FillHistogram("hPi0MC", particle->Pt(), particle->Y(), Weight(particle));
    else if (particle->GetPdgCode() == 211)
      FillHistogram("hPiPlusMC", particle->Pt(), particle->Y(), Weight(particle)); 
    else if (particle->GetPdgCode() == -211)
      FillHistogram("hPiMinusMC", particle->Pt(), particle->Y(), Weight(particle));                  
    else if (TMath::Abs(particle->GetPdgCode()) == 221)
         FillHistogram("hEtaMC", particle->Pt(), particle->Y(), Weight(particle)); 
    else if (TMath::Abs(particle->GetPdgCode()) == 331)
         FillHistogram("hEtaPrimeMC", particle->Pt(), particle->Y(), Weight(particle)); 
    else if (TMath::Abs(particle->GetPdgCode()) == 223)
         FillHistogram("hOmegaMC", particle->Pt(), particle->Y(), Weight(particle)); 
    else if (TMath::Abs(particle->GetPdgCode()) == 130)
         FillHistogram("hK0LMC", particle->Pt(), particle->Y(), Weight(particle));        
    else if (TMath::Abs(particle->GetPdgCode()) == 310)
         FillHistogram("hK0SMC", particle->Pt(), particle->Y(),  Weight(particle));         
    else if (TMath::Abs(particle->GetPdgCode() == 2212)) 
         FillHistogram("hProtonMC",particle->Pt(), particle->Y(), Weight(particle));
    else if (TMath::Abs(particle->GetPdgCode()==2112)) 
         FillHistogram("hNeutrMC",particle->Pt(), particle->Y(), Weight(particle));
    else if (TMath::Abs(particle->GetPdgCode()==321)) 
         FillHistogram("hKchMC",particle->Pt(), particle->Y(), Weight(particle));
    else if (TMath::Abs(particle->GetPdgCode()) == 22 || TMath::Abs(particle->GetPdgCode()) == 11) {           
        if (TMath::Abs(particle->GetPdgCode()) == 22)
           FillHistogram("hGammaMC_all", particle->Pt(), particle->Y());

     if (particle->IsSecondaryFromMaterial()) {
           FillHistogram(Form("h%sMC_FromMaterial",  particle->GetPdgCode() == 22 ? "Gamma" : "Beta"), 
                                                      particle->Pt(), particle->Y(), Weight(particle));
     }						       
     else {
        Int_t iMother = particle->GetMother();
        AliAODMCParticle* mparticle = (AliAODMCParticle*) fMCArray->At(iMother);
        if (mparticle->GetPdgCode() == 111 && mparticle->GetLabel() > -1) {
          Int_t iMother2 = mparticle->GetMother();
          if (((AliAODMCParticle*) fMCArray->At(iMother2))->GetPdgCode() == 310 
              || ((AliAODMCParticle*) fMCArray->At(iMother2))->GetPdgCode() == 130) { 
            mparticle = (AliAODMCParticle*) fMCArray->At(iMother2);
          }  
     } 
     FillHistogram("hGammaMCSources", particle->Pt(), mparticle->GetPdgCode(), Weight(mparticle));
     if (TMath::Hypot(particle->Xv(), particle->Yv()) < 1.0 
                             && mparticle->GetPdgCode() != 130
                             && mparticle->GetPdgCode() != 310
                             && particle->GetPdgCode() == 22) {              
       FillHistogram("hGammaMC_true", particle->Pt(), particle->Y(), Weight(mparticle));
      }	    
     }     
    }
  }
}
//------------------------------------------------------------------------------
void AliAnalysisTaskGammaPHOSPP::AnalyzeCells()
{
  AliAODCaloCells *cells = fEvent->GetPHOSCells();

  Int_t multCells = cells->GetNumberOfCells();
  FillHistogram("hCellMultEvent",multCells);

  Float_t  energy;
  Int_t    mod1, relId[4], cellAbsId, cellX, cellZ;

  // Single loop over cells

  Int_t nCellModule[4] = {0, 0, 0, 0};
  for (Int_t iCell=0; iCell<multCells; iCell++) {
    cellAbsId = cells->GetCellNumber(iCell);
    fPHOSGeo->AbsToRelNumbering(cellAbsId,relId);
    mod1  = relId[0];
    cellX = relId[2];
    cellZ = relId[3] ;
    energy = cells->GetAmplitude(iCell);
    FillHistogram("hCellEnergy",energy);
    if (mod1==1) {
      nCellModule[0]++;
      FillHistogram("hCellEnergyM1",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM1",cellX,cellZ,1.);
      FillHistogram("hCellEXZM1",cellX,cellZ,energy);
     } else if (mod1==2) {
      nCellModule[1]++;
      FillHistogram("hCellEnergyM2",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM2",cellX,cellZ,1.);
      FillHistogram("hCellEXZM2",cellX,cellZ,energy);
    } else if (mod1==3) {
      nCellModule[2]++;
      FillHistogram("hCellEnergyM3",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM3",cellX,cellZ,1.);
      FillHistogram("hCellEXZM3",cellX,cellZ,energy);
    } else if (mod1==4) {
      nCellModule[3]++;
      FillHistogram("hCellEnergyM4",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM4",cellX,cellZ,1.);
      FillHistogram("hCellEXZM4",cellX,cellZ,energy);
    }
  }
  FillHistogram("hCellMultEventM1",nCellModule[0]);
  FillHistogram("hCellMultEventM2",nCellModule[1]);
  FillHistogram("hCellMultEventM3",nCellModule[2]);
  FillHistogram("hCellMultEventM4",nCellModule[3]);
}
//------------------------------------------------------------------------------
void AliAnalysisTaskGammaPHOSPP::SelectCluster(AliAODCaloCluster *clu1)
{
  //Analyze clusters and select photons for analysis
  Int_t  multPHOSClust[5]  = {0, 0, 0, 0, 0} ; 
  
  TLorentzVector p1, p11 ;
  
  Int_t    mod1, relId[4], cellAbsId, cellX, cellZ ;
  Float_t  position[3] ;  
  Int_t digMult ;
  Double_t energy, weight ;

  if (!clu1->IsPHOS()) return ; 
       
  clu1->GetPosition(position) ;
  TVector3 global1(position)  ;
  cellAbsId = clu1->GetCellAbsId(0) ;
  fPHOSGeo->AbsToRelNumbering(cellAbsId, relId) ;
  mod1  = relId[0] ;
  cellX = relId[2] ;
  cellZ = relId[3] ;
  energy = clu1->E() ;
  digMult = clu1->GetNCells() ; 

  if (mod1 < 1 || mod1 > 4) {
    Printf("Wrong module number %d", mod1) ;
    return ;
  }

  FillHistogram("hDistanceToBadChannel", clu1->GetDistanceToBadChannel());
  
  if (clu1->GetType() != AliVCluster::kPHOSNeutral) return ;

  if (fLHCRunN == 2 && !fMCArray && TMath::Abs(clu1->GetTOF()) > fTOFcut) return ; // TOF cut for real data only!
  
  multPHOSClust[0]++ ;
  FillHistogram("hClusterEnergy",energy) ;
  FillHistogram("hCellMultClu_all",digMult) ;
  FillHistogram("hClusterEvsN_all",energy,digMult) ;
  FillHistogram(Form("hClusterEnergyM%d", mod1), energy) ;
  FillHistogram(Form("hClusterEvsN_all_M%d", mod1), energy, digMult) ;
  FillHistogram(Form("hCellMultClu_all_M%d", mod1), digMult) ;
  FillHistogram(Form("hCluNXZM%d", mod1), cellX, cellZ, 1.) ;
  FillHistogram(Form("hCluEXZM%d", mod1), cellX, cellZ, clu1->E());

  if (clu1->GetEmcCpvDistance() > 2.5) {
    FillHistogram("hClusterEvsN_cpv",energy,digMult) ;
    FillHistogram("hCellMultClu_cpv",digMult) ;
    FillHistogram(Form("hClusterEvsN_cpv_M%d", mod1), energy, digMult) ;
    FillHistogram(Form("hCellMultClu_cpv_M%d", mod1), digMult) ;
  }
  if (clu1->Chi2() < 2.5) {
    FillHistogram("hClusterEvsN_disp",energy,digMult);
    FillHistogram("hCellMultClu_disp",digMult);
    FillHistogram(Form("hClusterEvsN_disp_M%d", mod1), energy, digMult);
    FillHistogram(Form("hCellMultClu_disp_M%d", mod1), digMult);
  }
  if (clu1->GetEmcCpvDistance() > 2.5 && clu1->Chi2() < 2.5) {
    FillHistogram("hClusterEvsN_both",energy,digMult);
    FillHistogram("hCellMultClu_both",digMult);
    FillHistogram(Form("hClusterEvsN_both_M%d", mod1), energy,digMult);
    FillHistogram(Form("hCellMultClu_both_M%d", mod1), digMult);
  }

  if (mod1 == 1) multPHOSClust[1]++ ;
  else if (mod1 == 2) multPHOSClust[2]++;
  else if (mod1 == 3) multPHOSClust[3]++;
  else if (mod1 == 4) multPHOSClust[4]++;    
   
  AliPHOSAodCluster cluPHOS1(*clu1);

  Double_t en = clu1->E();

  cluPHOS1.SetE(NonlinearCorrection(en));

  cluPHOS1.GetMomentum(p1 , fVtx0);
  cluPHOS1.GetMomentum(p11, fVtx5);
   
  Double_t pAbs = p1.P();
  Double_t pT   = p1.Pt();
  Double_t pX   = p1.Px();
  Double_t pY   = p1.Py();

  if (pAbs<1.e-10) return ;

  Double_t kappa = pAbs - TMath::Power(0.135, 2)/4./pAbs;

  FillHistogram("hPhotonKappa", kappa);
  FillHistogram("hPhotonPt", pT);
  FillHistogram("hPhotonPx", pX);
  FillHistogram("hPhotonPy", pY);

  if (clu1->E() < 0.3) return;
  if (fLHCRunN == 1) {
     if (clu1->GetNCells() < 3)  return ; //no cuts on cell multiplicity only in MC for Run 2 
  } else {
     if (!fMCArray) {
       if (clu1->E() < 1. && clu1->GetNCells() < 1) return;
       if (clu1->E() > 1. && clu1->GetNCells() < 2) return;
     }
  }

  if (clu1->GetM02() < 0.1) return ;    

  FillHistogram("hEmcCPVDistance", clu1->GetEmcCpvDistance());
  Int_t iPrimaryAtVertex = GetPrimaryLabelAtVertex(clu1);
  
  if (!fMCArray) {
    weight = 1.0;
  }  
  else {
     weight = Weight((AliAODMCParticle*)fMCArray->At(iPrimaryAtVertex));       
  }

  if (!fMCArray && fLHCRunN == 1) {
    p11 *= fRecalib[mod1-1] ; 
    p1  *= fRecalib[mod1-1];
  }

  new((*fPHOSEvent)[fInPHOS]) AliCaloPhoton(p1.X(),p1.Y(),p1.Z(),p1.E());
  
  AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSEvent->At(fInPHOS) ;
  ph->SetModule(mod1) ;
  ph->SetMomV2(&p11) ;
  ph->SetNCells(clu1->GetNCells());
  ph->SetEMCx(global1.X()); 
  ph->SetEMCy(global1.Y());
  ph->SetEMCz(global1.Z());
  ph->SetNsigmaCPV(clu1->GetEmcCpvDistance());
  ph->SetNsigmaFullDisp(TMath::Sqrt(clu1->Chi2()));
  ph->SetCPVBit(ph->GetNsigmaCPV() < NsigmaCPV);
  ph->SetDispBit(ph->GetNsigmaFullDisp() < NsigmaDisp);
  ph->SetBC(TestBC(clu1->GetTOF()));
  ph->SetPrimary(GetPrimaryLabel(clu1));
  ph->SetPrimaryAtVertex(GetPrimaryLabelAtVertex(clu1));
  ph->SetWeight(weight); 
  
  TestMatchingTrackPID(ph, p11.Pt());

  FillHistogram("hvt0vsvt5", p11.Pt()- p1.Pt());
  FillHistogram("hBC", TestBC(clu1->GetTOF()) + 0.5);
  FillHistogram("hTOF", clu1->GetTOF());
  FillHistogram("hWeights", ph->GetWeight());      

  fInPHOS++ ;
}

//===========================================================================
void AliAnalysisTaskGammaPHOSPP::FillOnePhotonHistograms(AliCaloPhoton *ph)
{
   TLorentzVector p11 = *(ph->GetMomV2());
   
   if (fLHCRunN == 1 && ph->GetBC() != 0 ) return; //Run 1

   Double_t weight = ph->GetWeight();
   Int_t pdg = 0;
   Int_t pdg_naive = 0;
   if (fMCArray) {
     pdg = ((AliAODMCParticle*)fMCArray->At(ph->GetPrimaryAtVertex())) -> GetPdgCode();
     pdg_naive = ((AliAODMCParticle*)fMCArray->At(ph->GetPrimary())) -> GetPdgCode(); //
   }
   
   Int_t sm1 = ph->Module();

   Double_t pt   = p11.Pt();
   Double_t ptMC = TestGammaPt(ph);

   if (pdg == 22) {
      FillHistogram("hMatrixEff_all", pt, ptMC);   
      if (pdg_naive == 22)     
        FillHistogram("hMatrixEff_gamma_all", pt, ptMC); 
      if (TMath::Abs(pdg_naive) == 11)     
        FillHistogram("hMatrixEff_beta_all", pt, ptMC); 
   } 

   FillHistogram("hCaloPhotonPt_all", pt, weight);
   FillHistogram("hCaloPhotonPdgvsPt_all", pt, pdg, weight);
   FillHistogram("hCaloPhotonPdgvsPt_all_naive",  pt, pdg_naive, weight);

   FillHistogram("hCaloPhotonPtvsNcl_all", pt, ph->GetNCells() + 0.5);
   FillHistogram("hCentralityvsClustPt_all", pt, fEventCentrality + 0.5);	
   FillHistogram(Form("hCaloPhotonPt_all_M%d", sm1),  pt, weight);

   if (ph->IsCPVOK()) {
     FillHistogram("hCaloPhotonPt_cpv", pt, weight );
     FillHistogram("hCaloPhotonPdgvsPt_cpv", pt, pdg, weight); 
     FillHistogram("hCaloPhotonPdgvsPt_cpv_naive",  pt, pdg_naive, weight);  

     FillHistogram("hCaloPhotonPtvsNcl_cpv", pt, ph->GetNCells()  + 0.5);
     FillHistogram("hCentralityvsClustPt_cpv", pt, fEventCentrality + 0.5);	
     FillHistogram(Form("hCaloPhotonPt_cpv_M%d", sm1),  pt, weight);

     if (pdg == 22) {
       FillHistogram("hMatrixEff_cpv", pt, ptMC, weight);
       if (pdg_naive == 22)
         FillHistogram("hMatrixEff_gamma_cpv", pt, ptMC, weight);
       if (TMath::Abs(pdg_naive) == 11)
         FillHistogram("hMatrixEff_beta_cpv", pt, ptMC, weight);  
     }        
   }

   if (ph->IsDispOK()) {
     FillHistogram("hCaloPhotonPt_disp", pt, weight);
     FillHistogram("hCaloPhotonPdgvsPt_disp", pt, pdg, weight);
     FillHistogram("hCaloPhotonPdgvsPt_disp_naive",  pt, pdg_naive, weight);

     FillHistogram("hCaloPhotonPtvsNcl_disp", pt, ph->GetNCells()  + 0.5);
     FillHistogram("hCentralityvsClustPt_disp", pt, fEventCentrality + 0.5);	
     FillHistogram(Form("hCaloPhotonPt_disp_M%d", sm1),  pt, weight);

     if (pdg == 22) {
       FillHistogram("hMatrixEff_disp", pt, ptMC, weight);
       if (pdg_naive == 22)
         FillHistogram("hMatrixEff_gamma_disp", pt, ptMC, weight);
       if (TMath::Abs(pdg_naive) == 11)
         FillHistogram("hMatrixEff_beta_disp", pt, ptMC, weight);  
     }                    
   }
   if (ph->IsCPVOK() && ph->IsDispOK()) { 
     FillHistogram("hCaloPhotonPt_both", pt, weight);
     FillHistogram("hCaloPhotonPdgvsPt_both", pt, pdg, weight);
     FillHistogram("hCaloPhotonPdgvsPt_both_naive",  pt, pdg_naive, weight);

     FillHistogram("hCaloPhotonPtvsNcl_both", pt, ph->GetNCells() + 0.5, weight);
     FillHistogram("hCentralityvsClustPt_both", pt, fEventCentrality + 0.5, weight);	
     FillHistogram(Form("hCaloPhotonPt_both_M%d", sm1),  pt, weight);

     if (pdg == 22) {
       FillHistogram("hMatrixEff_both", pt, ptMC, weight);
       if (pdg_naive == 22)
         FillHistogram("hMatrixEff_gamma_both", pt, ptMC, weight);
       if (TMath::Abs(pdg_naive) == 11)
         FillHistogram("hMatrixEff_beta_both", pt, ptMC, weight);  
     }      
   }
}

//============================================================================
void AliAnalysisTaskGammaPHOSPP::FillTwoPhotonHistograms()
{
  TLorentzVector p1, p2, p12, pv1, pv2, pv12, p11;

  for (Int_t i1 = 0; i1 < fInPHOS-1; i1 ++ ) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    if (fLHCRunN == 1 && ph1->GetBC() != 0 ) continue; //Run 1
    for (Int_t i2 = i1+1; i2 < fInPHOS; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent->At(i2) ;
      if (fLHCRunN == 1 && ph2->GetBC() != 0 ) continue; //Run 1
      pv1 = *(ph1->GetMomV2());
      pv2 = *(ph2->GetMomV2());
      Double_t P1 = pv1.P();
      Double_t P2 = pv2.P();
      Double_t asym = TMath::Abs((P1-P2)/(P1+P2));
      FillHistogram("hAsym", asym);
      p12  = *ph1  + *ph2;
      pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
      Double_t ma12 = pv12.M();
      Double_t pt12 = pv12.Pt();
      Int_t sm1 = ph1->Module();
      Int_t sm2 = ph2->Module();

      FillHistogram("hMassPt_all", ma12 , pt12 );
      FillHistogram("hMassPt_asym_all", ma12, pt12, asym);
      if (sm1 == sm2) FillHistogram(Form("hMassPt_all_M%d", sm1), ma12, pt12);

      if (ph1->IsCPVOK() && ph2->IsCPVOK()) {
         FillHistogram("hMassPt_cpv", ma12, pt12);
         FillHistogram("hMassPt_asym_cpv", ma12, pt12, asym);
         if (sm1 ==sm2) FillHistogram(Form("hMassPt_cpv_M%d", sm1), ma12, pt12);
      }

      if (ph1->IsDispOK() && ph2->IsDispOK() ) {
         FillHistogram("hMassPt_disp", ma12, pt12);
         FillHistogram("hMassPt_asym_disp", ma12, pt12, asym);
         if (sm1 ==sm2) FillHistogram(Form("hMassPt_disp_M%d", sm1), ma12, pt12);
         if (ph1->IsCPVOK() && ph2->IsCPVOK() ) {
            FillHistogram("hMassPt_both", ma12, pt12);
            FillHistogram("hMassPt_asym_both", ma12, pt12, asym);
            if (sm1 ==sm2) FillHistogram(Form("hMassPt_both_M%d", sm1), ma12, pt12);
            if (sm1 ==sm2 && asym < 0.1) FillHistogram(Form("hMassPt_both_M%d", sm1), ma12, pt12);
         }
      }

      FillHistogram("hMassSingle_all", ma12, ph1->Pt()) ;
      FillHistogram("hMassSingle_all", ma12, ph2->Pt()) ;
      FillHistogram("hMassSingle_asym_all", ma12, ph1->Pt(), asym);
      FillHistogram("hMassSingle_asym_all", ma12, ph2->Pt(), asym);

      FillHistogram(Form("hMassSingle_all_M%d", sm1), ma12, ph1->Pt());
      FillHistogram(Form("hMassSingle_all_M%d", sm2), ma12, ph2->Pt());

      if (ph1->IsCPVOK()) {
        FillHistogram("hMassSingle_cpv",ma12,ph1->Pt()) ;
        FillHistogram(Form("hMassSingle_cpv_M%d", sm1), ma12, ph1->Pt());
        FillHistogram("hMassSingle_asym_cpv", ma12, ph1->Pt(), asym);
      }

      if (ph2->IsCPVOK()) {
        FillHistogram("hMassSingle_cpv",ma12,ph2->Pt()) ;
        FillHistogram(Form("hMassSingle_cpv_M%d", sm2), ma12, ph2->Pt());
        FillHistogram("hMassSingle_asym_cpv", ma12, ph2->Pt(), asym);
      }

      if (ph1->IsDispOK()) {
        FillHistogram("hMassSingle_disp",ma12,ph1->Pt()) ;   
        FillHistogram(Form("hMassSingle_disp_M%d", sm1), ma12, ph1->Pt());
        FillHistogram("hMassSingle_asym_disp", ma12, ph1->Pt(), asym);
      }
      if (ph2->IsDispOK()) {
        FillHistogram("hMassSingle_disp", ma12, ph2->Pt()) ;
        FillHistogram(Form("hMassSingle_disp_M%d", sm2), ma12, ph2->Pt());
        FillHistogram("hMassSingle_asym_disp", ma12, ph2->Pt(), asym);
      }

      if (ph1->IsCPVOK() && ph1->IsDispOK()) {
        FillHistogram("hMassSingle_both", ma12, ph1->Pt()) ;
        FillHistogram(Form("hMassSingle_both_M%d", sm1), ma12, ph1->Pt());
        FillHistogram("hMassSingle_asym_both", ma12, ph1->Pt(), asym);
      }   
      if (ph2->IsCPVOK() && ph2->IsDispOK()) {
        FillHistogram("hMassSingle_both",ma12,ph2->Pt()) ;
        FillHistogram(Form("hMassSingle_both_M%d", sm2), ma12, ph2->Pt());
        FillHistogram("hMassSingle_asym_both", ma12, ph2->Pt(), asym);
      }
  } //end of loop  i2 
 } //end of loop   i1
}
//==============================================================================
void AliAnalysisTaskGammaPHOSPP::MixPhotons()
{
  TLorentzVector p1, p2, p12, pv1, pv2, pv12, p11;

  Int_t zvtx = (Int_t)((fVtx5[2] + 10.)/2.) ;
  if(zvtx < 0) zvtx = 0 ;
  if(zvtx > 9) zvtx = 9 ;
  
  Int_t centr = 0;
  
  if(!fPHOSEvents[zvtx][centr]) fPHOSEvents[zvtx][centr]=new TList() ;
  
  TList * prevPHOS = fPHOSEvents[zvtx][centr] ;

  for (Int_t i1 = 0; i1 < fInPHOS; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    if (fLHCRunN == 1 && ph1->GetBC() != 0 ) continue; //Run 1
    for (Int_t ev = 0; ev < prevPHOS->GetSize(); ev ++) {
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev)) ;
      for (Int_t i2=0; i2 < mixPHOS->GetEntriesFast();i2++) {
        AliCaloPhoton * ph2=(AliCaloPhoton*)mixPHOS->At(i2) ;
        if (fLHCRunN == 1 && ph2->GetBC() != 0 ) continue; //Run 1
        pv1 = *(ph1->GetMomV2());
        pv2 = *(ph2->GetMomV2());
        Double_t P1 = pv1.P();
        Double_t P2 = pv2.P();
        Double_t asym = TMath::Abs((P1-P2)/(P1+P2));
        FillHistogram("hAsym_mix", asym);
        p12  = *ph1  + *ph2;
        pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
        Double_t ma12 = pv12.M();
        Double_t pt12 = pv12.Pt();
        Int_t sm1 = ph1->Module();
        Int_t sm2 = ph2->Module();

        if (fLHCRunN == 1 || (fLHCRunN == 2 && !fMCArray)) {
           if (ph1->GetNCells() < 3 || ph2->GetNCells() < 3) return;      
        }

      FillHistogram("hMiMassPt_all", ma12 , pt12 );
      if (sm1 == sm2) FillHistogram(Form("hMiMassPt_all_M%d", sm1), ma12, pt12);

      if (ph1->IsCPVOK() && ph2->IsCPVOK()) {
         FillHistogram("hMiMassPt_cpv", ma12, pt12);
         if (sm1 ==sm2) FillHistogram(Form("hMiMassPt_cpv_M%d", sm1), ma12, pt12);
      }

      if (ph1->IsDispOK() && ph2->IsDispOK() ) {
         FillHistogram("hMiMassPt_disp", ma12, pt12);
         if (sm1 ==sm2) FillHistogram(Form("hMiMassPt_disp_M%d", sm1), ma12, pt12);
         if (ph1->IsCPVOK() && ph2->IsCPVOK() ) {
            FillHistogram("hMiMassPt_both", ma12, pt12);
            if (sm1 ==sm2) FillHistogram(Form("hMiMassPt_both_M%d", sm1), ma12, pt12);
            if (sm1 ==sm2 && asym < 0.1) FillHistogram(Form("hMiMassPt_both_M%d", sm1), ma12, pt12);
         }
      }

      FillHistogram("hMiMassSingle_all", ma12, ph1->Pt()) ;
      FillHistogram("hMiMassSingle_all", ma12, ph2->Pt()) ;
      FillHistogram("hMiMassSingle_asym_all", ma12, ph1->Pt(), asym);
      FillHistogram("hMiMassSingle_asym_all", ma12, ph2->Pt(), asym);

      FillHistogram(Form("hMiMassSingle_all_M%d", sm1), ma12, ph1->Pt());
      FillHistogram(Form("hMiMassSingle_all_M%d", sm2), ma12, ph2->Pt());

      if (ph1->IsCPVOK()) {
        FillHistogram("hMiMassSingle_cpv",ma12,ph1->Pt()) ;
        FillHistogram(Form("hMiMassSingle_cpv_M%d", sm1), ma12, ph1->Pt());
        FillHistogram("hMiMassSingle_asym_cpv", ma12, ph1->Pt(), asym);
      }

      if (ph2->IsCPVOK()) {
        FillHistogram("hMiMassSingle_cpv",ma12,ph2->Pt()) ;
        FillHistogram(Form("hMiMassSingle_cpv_M%d", sm2), ma12, ph2->Pt());
        FillHistogram("hMiMassSingle_asym_cpv", ma12, ph2->Pt(), asym);
      }

      if (ph1->IsDispOK()) {
        FillHistogram("hMiMassSingle_disp",ma12,ph1->Pt()) ;   
        FillHistogram(Form("hMiMassSingle_disp_M%d", sm1), ma12, ph1->Pt());
        FillHistogram("hMiMassSingle_asym_disp", ma12, ph1->Pt(), asym);
      }
      if (ph2->IsDispOK()) {
        FillHistogram("hMiMassSingle_disp", ma12, ph2->Pt()) ;
        FillHistogram(Form("hMiMassSingle_disp_M%d", sm2), ma12, ph2->Pt());
        FillHistogram("hMiMassSingle_asym_disp", ma12, ph2->Pt(), asym);
      }

      if (ph1->IsCPVOK() && ph1->IsDispOK()) {
        FillHistogram("hMiMassSingle_both", ma12, ph1->Pt()) ;
        FillHistogram(Form("hMiMassSingle_both_M%d", sm1), ma12, ph1->Pt());
        FillHistogram("hMiMassSingle_asym_both", ma12, ph1->Pt(), asym);
      }   
      if (ph2->IsCPVOK() && ph2->IsDispOK()) {
        FillHistogram("hMiMassSingle_both",ma12,ph2->Pt()) ;
        FillHistogram(Form("hMiMassSingle_both_M%d", sm2), ma12, ph2->Pt());
        FillHistogram("hMiMassSingle_asym_both", ma12, ph2->Pt(), asym);
      }
    } // end of loop i2
  }
 } // end of loop i1
 
  if (fPHOSEvent->GetEntriesFast() > 0) {
    prevPHOS->AddFirst(fPHOSEvent) ;
    fPHOSEvent=0;
    if (prevPHOS->GetSize() > 100) {//Remove redundant events
      TClonesArray * tmp = static_cast<TClonesArray*>(prevPHOS->Last()) ;
      prevPHOS->RemoveLast() ;
      delete tmp ;
    }
  }
}

//=================================== Returns label at vertex ==================
Int_t AliAnalysisTaskGammaPHOSPP::GetPrimaryLabelAtVertex(AliVCluster *clu)
{
   if (!fMCArray) 
     return 0;
      
   Int_t iPrimaryAtVertex = clu->GetLabel();
   AliAODMCParticle *particle0 =  (AliAODMCParticle*) fMCArray->At(iPrimaryAtVertex);
   
   if (particle0-> IsSecondaryFromMaterial()) {
    // Printf("Secondary from the material, Epart = %f, Eclust = %f" , particle0 ->E(), clu->E());
   //  return 0;
     FillHistogram("hCaloPhotonPdgvsPt_FromMaterial", clu->E(), particle0->GetPdgCode());
   }

   while (TMath::Hypot(particle0->Xv(), particle0->Yv()) > 1.0) {
      iPrimaryAtVertex = particle0->GetMother();
      particle0 = (AliAODMCParticle*) fMCArray->At(particle0->GetMother());
   }

   Int_t nn = iPrimaryAtVertex;

   if (particle0->GetPdgCode() == 22 || particle0->GetPdgCode() == 11) {
     //iPrimaryAtVertex = particle0->GetMother();
     for (Int_t i = 0; i < nn ; i++) {
       AliAODMCParticle* particle =  (AliAODMCParticle*) fMCArray->At(i);
       if(particle->GetPdgCode() != 310 && particle->GetPdgCode() != 130) 
         continue;
       Int_t iSecondDaughter = particle->GetDaughterLabel(1); 
       if(iSecondDaughter != iPrimaryAtVertex) 
         continue;
       else
         iPrimaryAtVertex = i;
     }
   } else 
       iPrimaryAtVertex = nn;
  
   return iPrimaryAtVertex;   
}
//=========================== Returns label of the impinging particle =========
Int_t AliAnalysisTaskGammaPHOSPP::GetPrimaryLabel(AliVCluster *clu)
{
   if (!fMCArray) 
     return 0;
      
   return clu->GetLabel();
}
//=================================== TestGamma returns momentum==============
Double_t AliAnalysisTaskGammaPHOSPP::TestGammaPt(AliCaloPhoton *ph)
{

   if (!fMCArray) 
     return 0;

   Int_t iPrimaryAtVertex = ph->GetPrimaryAtVertex();
   AliAODMCParticle *particle0 = (AliAODMCParticle*) fMCArray->At(iPrimaryAtVertex);

   return (particle0->Pt());
}

//=======================================
Int_t AliAnalysisTaskGammaPHOSPP::TestTrack(AliAODTrack *track)
{
   if (!fMCArray) 
     return 0;

   Int_t TrackLabel=track->GetLabel();

   //if(TrackLabel < -2) return 0;

   AliAODMCParticle *TrackParticle =  (AliAODMCParticle*) fMCArray->At(TrackLabel);
   if (!TrackParticle) 
     return 0;

 //  if (((AliAODMCParticle*) TrackParticle)->IsSecondaryFromWeakDecay()) return 0;
 //  if (!TrackParticle->IsPhysicalPrimary()) return 0;

   Int_t TrackPDG = TrackParticle->GetPdgCode();

   return TrackPDG;
}
//=======================================
Double_t AliAnalysisTaskGammaPHOSPP::Weight(AliAODMCParticle *particleAtVertex)
{
   //mt scaling for Run1 MC production
   if (!fMCArray) 
     return 1.0; // Data

   if (fLHCRunN == 2)
     return 1.0; // Run 2
   
   Double_t pdg = TMath::Abs(particleAtVertex->GetPdgCode());

   if (pdg == 11 || pdg == 22) {
           particleAtVertex = (AliAODMCParticle*)fMCArray->At(particleAtVertex->GetMother());
           pdg = TMath::Abs(particleAtVertex->GetPdgCode());
   }

   if (pdg == 111 || pdg  == 211 ) {
          // fWeightFunction->SetParameters(0.611073, -0.0222529, 0.190541, -0.416579, 0.396059, 0.611073);
	  Double_t parsPi0[6] = {0.54556579, -0.16503941, 0.24846410, -0.65202206, 0.51017892, 0.0033008629};
          fWeightFunction->SetParameters(parsPi0);
   }   else if (pdg == 221 || pdg == 331 || pdg == 223 ) {   
         //fWeightFunction->SetParameters(0.0601459, 0, 4.11665, 0, 6.46838, -0.00319589);
	 Double_t parsEta[6]={0.692057, -1.60313, 1.12976, -2.15741, 1.46253, 0.00414204};
	 fWeightFunction->SetParameters(parsEta);
       } else if(pdg == 130 || pdg == 310 || pdg == 311 || pdg == 321 ) {    
             //fWeightFunction->SetParameters(0.708656, 0.355564, -0.00468263, 0.0570132, 0.076876, 0.0382327);
	     Double_t parsKa[6] = {0.718241, 0.173219, 0.0504723, -0.111518, 0.118757, 0.0134929};
             fWeightFunction->SetParameters(parsKa);
         } else if (pdg > 1000) {
             //fWeightFunction->SetParameters(0.215726, 0.292934, 0.163074, -0.460113, 0.219988, -0.0903996);
	     Double_t parsPN[6] = {0.209599, 0.298912, 0.139915, -0.450081, 0.207574, -0.0861399};
	     fWeightFunction->SetParameters(parsPN);
           } else fWeightFunction->SetParameters(1.0, 0., 0., 0., 0., 0.);

   return fWeightFunction->Eval(particleAtVertex->Pt());
}

//=======================================
Int_t AliAnalysisTaskGammaPHOSPP::TestBC(Double_t tof)
{
  Int_t bc = (Int_t)(TMath::Ceil((tof + fBCgap/2)/fBCgap) - 1);
  return bc;
}

//===============================================
Bool_t AliAnalysisTaskGammaPHOSPP::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  //
  
  //if(!fCheckMCCrossSection) return kTRUE;

  // Fetch the aod also from the input in,
  // have todo it in notify
  
  Float_t xsection = 0;
  Float_t trials   = 1;
  fAvgTrials = -1;
  
  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if(!tree) return kFALSE;
  
  TFile *curfile = tree->GetCurrentFile();
  
  if(!curfile) return kFALSE;
  
  if(fCurrFileName == curfile->GetName()) return kFALSE;
  
  fCurrFileName = TString(curfile->GetName());
  
  if (!fh1Xsec||!fh1Trials) {
  //  Printf("%s%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__);
    return kFALSE;
  }
  
  Bool_t ok = PythiaInfoFromFile(fCurrFileName, xsection, trials);
  
  if (!ok) 
    return kFALSE;
  
  fh1Xsec->Fill("<#sigma>",xsection);
  
  // construct a poor man average trials
  Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
  
  if(trials >= nEntries && nEntries > 0.) fAvgTrials = trials/nEntries;
  
  fh1Trials->Fill("#sum{ntrials}",trials);
  
  //printf("AliAnalysisTaskGammaPHOSPP::Notify() - xs %f, trial %f, avg trials %f\n",xsection,trials, fAvgTrials);
  
  if (fDebug) 
    Printf("Reading File %s",fInputHandler->GetTree()->GetCurrentFile()->GetName());
  
  return kTRUE;
}

//_____________________________________________________________________________________________________
Bool_t AliAnalysisTaskGammaPHOSPP::PythiaInfoFromFile(TString file,Float_t & xsec,Float_t & trials)
{
  //
  // get the cross section and the trails either from pyxsec.root or from pysec_hists.root
  // This is to called in Notify and should provide the path to the AOD/ESD file
    
  xsec   = 0;
  trials = 1;
  
  if (file.Contains("root_archive.zip#")) {
    Ssiz_t pos1 = file.Index("root_archive",12,0,TString::kExact);
    Ssiz_t pos  = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  } else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  
  //Printf("%s",file.Data());
  
  TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); // problem that we cannot really test the existance of a file in a archive so we have to lvie with open error message from root
  if (!fxsec) {
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
    if (!fxsec) {
      // not a severe condition but inciate that we have no information
      return kFALSE;
    } else {
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0);
      if (!key) {
        fxsec->Close();
        return kFALSE;
      }
      
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if (!list) {
        fxsec->Close();
        return kFALSE;
      }
      
      xsec    = ((TProfile*)list->FindObject("h1Xsec"))  ->GetBinContent(1);
      trials  = ((TH1F*)    list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } else {
      TTree *xtree = (TTree*)fxsec->Get("Xsection");
      if(!xtree) {
        fxsec->Close();
        return kFALSE;
      }
      
      UInt_t   ntrials  = 0;
      Double_t  xsection  = 0;
      xtree->SetBranchAddress("xsection",&xsection);
      xtree->SetBranchAddress("ntrials",&ntrials);
      xtree->GetEntry(0);
      trials = ntrials;
      xsec = xsection;
      fxsec->Close();
  }
  return kTRUE;
}

//=============================================================================

Bool_t AliAnalysisTaskGammaPHOSPP::PhotonWithinPeak(Double_t Minv, Double_t pt)
{
  const Int_t Nbins = 33;

  Double_t xbins[Nbins+1] = { 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2,
                              2.2, 2.4, 2.6, 2.8, 3, 3.2, 
                              3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 7, 8, 10, 12, 14, 16, 
                              18, 20, 25};

  Double_t xmass  [Nbins] = { 0.13531, 0.13531,  0.13531,  0.13531,
                              0.13531, 0.135298, 0.135323, 0.135263, 0.135279, 0.135261, 
                              0.135328,0.135201, 0.135349, 0.135275, 0.1353, 0.1354, 
                              0.135014,0.13545,  0.135336, 0.134832, 0.135384, 0.135042, 
                              0.13505, 0.134911, 0.135474, 0.134483, 0.13551, 0.135829, 
                              0.130291,0.134933, 0.145022, 0.139125, 0.135036};

  Double_t xwidth [Nbins] = { 5.47384, 5.47384, 5.47384, 5.47384,
                              5.47384, 5.19915, 5.16869, 5.00068, 5.03612, 4.92791, 
		                      4.99442, 4.98119, 5.14243, 5.01174, 5.09855, 4.82253, 
                              5.10046, 5.09045, 5.22204, 5.20583, 5.30179, 5.00537, 
                              5.09301, 5.57206, 5.64171, 4.92138, 5.49863, 7.24342, 
                              8.35561, 5.34015, 8.93808, 2.65933, 1.80287};

  Int_t k = 0;

  while (pt > xbins[k] && k < Nbins+1) 
     k = k + 1;
  return(Minv < xmass[k] + xwidth[k] && Minv > xmass[k] - xwidth[k]);
}

//=================
void AliAnalysisTaskGammaPHOSPP::TestMatchingTrackPID(AliCaloPhoton *ph, Double_t pt)
{
    AliVCluster *clu1 = ph->GetCluster();

    const Bool_t CPVBit  = ph->IsCPVOK();
    const Bool_t DispBit = ph->IsDispOK();

    const Int_t NTracksMatched = clu1->GetNTracksMatched();

    const Double_t dx = clu1->GetTrackDx(); 
    const Double_t dz = clu1->GetTrackDz(); 

    const Double_t dist = TMath::Hypot(dx,dz);
   
    Double_t pBayesMatched[AliPID::kSPECIESC];

    AliPIDCombined *pidcomb=new AliPIDCombined();

    pidcomb->SetDefaultTPCPriors();
    //pidcomb->SetEnablePriors(kFALSE);
    pidcomb->SetSelectedSpecies(AliPID::kSPECIESC);
    pidcomb->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF|AliPIDResponse::kDetITS|AliPIDResponse::kDetTRD);
      
    if (NTracksMatched == 0) return;

    FillHistogram("hTracks_matched", 0.5);
     FillHistogram("hDistance", dist);
     AliAODTrack* trackMatched= dynamic_cast<AliAODTrack*>(clu1->GetTrackMatched(0));
 
     if (trackMatched->TestFilterBit(32) &&  trackMatched->GetTPCsignal() > 0.) {
   
       //UInt_t oo = pidcomb->ComputeProbabilities(trackMatched, fPIDResponse, pBayesMatched);
       
       Bool_t pidPion3 = kFALSE , pidKaon3 = kFALSE , pidProton3 = kFALSE , pidElectron3 = kFALSE, pidUndef3 = kFALSE ;
       Bool_t pidPion1 = kFALSE , pidKaon1 = kFALSE , pidProton1 = kFALSE,  pidElectron1 = kFALSE, pidUndef1 = kFALSE;
      
       Float_t  nsigmaElectron =   TMath::Abs( fPIDResponse->NumberOfSigmasTPC( trackMatched, AliPID::kElectron )) ;
       Float_t  nsigmaPion =       TMath::Abs( fPIDResponse->NumberOfSigmasTPC( trackMatched, AliPID::kPion )) ;
       Float_t  nsigmaKaon =       TMath::Abs( fPIDResponse->NumberOfSigmasTPC( trackMatched, AliPID::kKaon ) ) ;
       Float_t  nsigmaProton =     TMath::Abs( fPIDResponse->NumberOfSigmasTPC( trackMatched, AliPID::kProton ) );
   
      // smallest sigma
       if (( nsigmaPion < nsigmaKaon) && (nsigmaPion < nsigmaProton) && (nsigmaPion < nsigmaElectron) && (nsigmaPion < 3.)) 
         pidPion3= kTRUE;
       if (( nsigmaProton < nsigmaKaon) && (nsigmaProton < nsigmaPion) && (nsigmaProton < nsigmaElectron) &&  (nsigmaProton < 3.)) 
         pidProton3= kTRUE;
       if (( nsigmaKaon < nsigmaPion) && (nsigmaKaon < nsigmaProton) && (nsigmaKaon < nsigmaElectron) &&  (nsigmaKaon < 3.))
         pidKaon3= kTRUE;
       if (( nsigmaElectron < nsigmaKaon) && (nsigmaElectron < nsigmaPion) && (nsigmaElectron < nsigmaProton) && (nsigmaElectron < 3.))
         pidElectron3= kTRUE;
       if (!pidPion3 && !pidProton3 && !pidKaon3 && !pidElectron3)
         pidUndef3=kTRUE;
       if (pidPion3)    
         FillHistogram("hpid3", 0.5);     
       if (pidProton3)   
         FillHistogram("hpid3", 1.5);     
       if (pidKaon3)     
         FillHistogram("hpid3", 2.5);     
       if (pidElectron3) 
         FillHistogram("hpid3", 3.5);     
       if (pidUndef3)    
         FillHistogram("hpid3", 4.5);     
   
       if (( nsigmaPion < nsigmaKaon) && (nsigmaPion < nsigmaProton) && (nsigmaPion < nsigmaElectron) && (nsigmaPion < 1.))
         pidPion1= kTRUE;
       if (( nsigmaProton < nsigmaKaon) && (nsigmaProton < nsigmaPion) && (nsigmaProton < nsigmaElectron) &&  (nsigmaProton < 1.)) 
         pidProton1= kTRUE;
       if (( nsigmaKaon < nsigmaPion) && (nsigmaKaon < nsigmaProton) && (nsigmaKaon < nsigmaElectron) &&  (nsigmaKaon < 1.)) 
         pidKaon1= kTRUE;
       if (( nsigmaElectron < nsigmaKaon) && (nsigmaElectron < nsigmaPion) && (nsigmaElectron < nsigmaProton) && (nsigmaElectron < 1.)) pidElectron1= kTRUE;
       if (!pidPion1 && !pidProton1 && !pidKaon1 && !pidElectron1) pidUndef1=kTRUE;
   
       if (pidPion1)     FillHistogram("hpid1", 0.5);     
       if (pidProton1)   FillHistogram("hpid1", 1.5);     
       if (pidKaon1)     FillHistogram("hpid1", 2.5);     
       if (pidElectron1) FillHistogram("hpid1", 3.5);     
       if (pidUndef1)    FillHistogram("hpid1", 4.5);     
       
       const Int_t nmaxMatched = TMath::LocMax(AliPID::kSPECIESC, pBayesMatched);
       const Double_t distEmcCpv = clu1->GetEmcCpvDistance();
       const Int_t  trackPdg = TestTrack(trackMatched);  
   
       FillHistogram("hTracksofClusts", pt, nmaxMatched+0.5);
       FillHistogram("hEnTrackvsClust", clu1->E(), trackMatched->E());
        
       if (pidPion3 && trackMatched->Charge() > 0) {
         FillHistogram("hTracksOfPi_ThreeSigma", pt, distEmcCpv);
         FillHistogram("hTracksOfPi_ThreeSigma_label", pt, trackPdg );
         if (DispBit) {
           FillHistogram("hTracksOfPi_ThreeSigma_disp",       pt, distEmcCpv);
           FillHistogram("hTracksOfPi_ThreeSigma_disp_label", pt, trackPdg );
         }
       } 
       if (pidPion1 && trackMatched->Charge() > 0) {
         FillHistogram("hTracksOfPi_OneSigma",       pt,distEmcCpv);
         FillHistogram("hTracksOfPi_OneSigma_label", pt, trackPdg );
         if (DispBit) {
           FillHistogram("hTracksOfPi_OneSigma_disp",       pt, distEmcCpv);
           FillHistogram("hTracksOfPi_OneSigma_disp_label", pt, trackPdg );
         }
       } 
       if (pidPion3 && trackMatched->Charge() < 0) {
         FillHistogram("hTracksOfAntiPi_ThreeSigma", pt,distEmcCpv);
         FillHistogram("hTracksOfAntiPi_ThreeSigma_label", pt, trackPdg );
         if (DispBit) {
           FillHistogram("hTracksOfAntiPi_ThreeSigma_disp",       pt, distEmcCpv);
           FillHistogram("hTracksOfAntiPi_ThreeSigma_disp_label", pt, trackPdg );
         }
       } 
       if (pidPion1 && trackMatched->Charge() < 0) {
         FillHistogram("hTracksOfAntiPi_OneSigma", pt,distEmcCpv);
         FillHistogram("hTracksOfAntiPi_OneSigma_label", pt, trackPdg );
         if (DispBit) {
           FillHistogram("hTracksOfAntiPi_OneSigma_disp",       pt, distEmcCpv);
           FillHistogram("hTracksOfAntiPi_OneSigma_disp_label", pt, trackPdg );
         }
       } 
       if (pidProton3 && trackMatched->Charge() > 0) { 
         FillHistogram("hTracksOfPr_ThreeSigma",pt,distEmcCpv);
         FillHistogram("hTracksOfPr_ThreeSigma_label", pt, trackPdg );
         if (DispBit) {
           FillHistogram("hTracksOfPr_ThreeSigma_disp_label", pt, trackPdg );
           FillHistogram("hTracksOfPr_ThreeSigma_disp",pt,distEmcCpv);
         }
       } 
       if (pidProton1 && trackMatched->Charge() > 0) { 
         FillHistogram("hTracksOfPr_OneSigma",pt,distEmcCpv);
         FillHistogram("hTracksOfPr_OneSigma_label", pt, trackPdg );
         if (DispBit) {
           FillHistogram("hTracksOfPr_OneSigma_disp_label", pt, trackPdg );
           FillHistogram("hTracksOfPr_OneSigma_disp",pt,distEmcCpv);
         }
       }
       if (pidProton3 && trackMatched->Charge() < 0) { 
         FillHistogram("hTracksOfAntiPr_ThreeSigma",pt,distEmcCpv);
         FillHistogram("hTracksOfAntiPr_ThreeSigma_label",pt,trackPdg);
         if (DispBit) {
           FillHistogram("hTracksOfAntiPr_ThreeSigma_disp",       pt, distEmcCpv);
           FillHistogram("hTracksOfAntiPr_ThreeSigma_disp_label", pt, trackPdg);
         }
       }
       if (pidProton1 && trackMatched->Charge() < 0) { 
         FillHistogram("hTracksOfAntiPr_OneSigma",pt,distEmcCpv);
         FillHistogram("hTracksOfAntiPr_OneSigma_label",pt,trackPdg);
         if (DispBit) {
           FillHistogram("hTracksOfAntiPr_OneSigma_disp",       pt, distEmcCpv);
           FillHistogram("hTracksOfAntiPr_OneSigma_disp_label", pt, trackPdg);
         }
       }
       if (pidKaon3 && trackMatched->Charge() > 0) { 
         FillHistogram("hTracksOfKa_ThreeSigma",       pt,distEmcCpv);
         FillHistogram("hTracksOfKa_ThreeSigma_label", pt, trackPdg );
         if (DispBit) {
           FillHistogram("hTracksOfKa_ThreeSigma_disp",       pt, distEmcCpv);
           FillHistogram("hTracksOfKa_ThreeSigma_disp_label", pt, trackPdg );
         }
       } 
       if (pidKaon1 && trackMatched->Charge() > 0) { 
         FillHistogram("hTracksOfKa_OneSigma",       pt,distEmcCpv);
         FillHistogram("hTracksOfKa_OneSigma_label", pt, trackPdg );
         if (DispBit) {
           FillHistogram("hTracksOfKa_OneSigma_disp",       pt, distEmcCpv);
           FillHistogram("hTracksOfKa_OneSigma_disp_label", pt, trackPdg );
         }
       } 
       if (pidKaon3 && trackMatched->Charge() < 0) { 
         FillHistogram("hTracksOfAntiKa_ThreeSigma",       pt,distEmcCpv);
         FillHistogram("hTracksOfAntiKa_ThreeSigma_label", pt, trackPdg );
         if (DispBit) { 
           FillHistogram("hTracksOfAntiKa_ThreeSigma_disp",   pt,distEmcCpv);
           FillHistogram("hTracksOfAntiKa_ThreeSigma_disp_label", pt, trackPdg );
         }
       } 
       if (pidKaon1 && trackMatched->Charge() < 0) { 
         FillHistogram("hTracksOfAntiKa_OneSigma",       pt,distEmcCpv);
         FillHistogram("hTracksOfAntiKa_OneSigma_label", pt, trackPdg );
         if (DispBit) { 
           FillHistogram("hTracksOfAntiKa_OneSigma_disp",       pt,distEmcCpv);
           FillHistogram("hTracksOfAntiKa_OneSigma_disp_label", pt, trackPdg );
         }
       } 
       if (pidElectron3 && trackMatched->Charge() > 0) {
         FillHistogram("hTracksOfAntiBeta_ThreeSigma",	  pt,distEmcCpv);
         FillHistogram("hTracksOfAntiBeta_ThreeSigma_label", pt, trackPdg );
         if (DispBit) {
           FillHistogram("hTracksOfAntiBeta_ThreeSigma_disp",	 pt, distEmcCpv);
           FillHistogram("hTracksOfAntiBeta_ThreeSigma_disp_label", pt, trackPdg );
         }
  	}
  	if (pidElectron1 && trackMatched->Charge() > 0) {
         FillHistogram("hTracksOfAntiBeta_OneSigma",	pt,distEmcCpv);
         FillHistogram("hTracksOfAntiBeta_OneSigma_label", pt, trackPdg );
         if (DispBit) {
           FillHistogram("hTracksOfAntiBeta_OneSigma_disp",       pt, distEmcCpv);
           FillHistogram("hTracksOfAntiBeta_OneSigma_disp_label", pt, trackPdg );
         }
  	}
       if (pidElectron3 && trackMatched->Charge() < 0) {
         FillHistogram("hTracksOfBeta_ThreeSigma",	  pt,distEmcCpv);
         FillHistogram("hTracksOfBeta_ThreeSigma_label", pt, trackPdg );
         if (DispBit) {
           FillHistogram("hTracksOfBeta_ThreeSigma_disp",	 pt, distEmcCpv);
           FillHistogram("hTracksOfBeta_ThreeSigma_disp_label", pt, trackPdg );
         }
  	}
  	if (pidElectron1 && trackMatched->Charge() < 0) {
         FillHistogram("hTracksOfBeta_OneSigma",	pt,distEmcCpv);
         FillHistogram("hTracksOfBeta_OneSigma_label", pt, trackPdg );
         if (DispBit) {
           FillHistogram("hTracksOfBeta_OneSigma_disp",       pt, distEmcCpv);
           FillHistogram("hTracksOfBeta_OneSigma_disp_label", pt, trackPdg );
         }
  	} else {
           FillHistogram("hTracksOfOthers", pt, distEmcCpv);
           if (DispBit) 
             FillHistogram("hTracksOfOthers_disp", pt, distEmcCpv);
       }
  
       if (CPVBit) {
         if (nmaxMatched == 2) { 
           FillHistogram("hTracksOfPiClose", pt);
           if (DispBit) 
             FillHistogram("hTracksOfPiCloseDispOK", pt);
         } else if(nmaxMatched==4) { 
             FillHistogram("hTracksOfPrClose", pt);
             if (DispBit) 
               FillHistogram("hTracksOfPrCloseDispOK", pt);
         } else if (nmaxMatched==3) { 
             FillHistogram("hTracksOfKaClose", pt);
             if (DispBit) 
                FillHistogram("hTracksOfKaCloseDispOK", pt);
         } else {
             FillHistogram("hTracksOfOthersClose",pt);
             FillHistogram("hTracksOfOthersCloseDispOK",pt);
         }               
       }
     }
}

//===========================================================================//
void AliAnalysisTaskGammaPHOSPP::PHOSvsEMCALClusters()
{
   Int_t multPHOS = 0 , multEMCAL = 0;

   AliAODCaloCluster *clu2; 
   for (Int_t ic = 0; ic < fEvent->GetNumberOfCaloClusters(); ic++) {
      clu2 = fEvent->GetCaloCluster(ic);
      if (clu2->IsPHOS() ) 
        multPHOS  =  multPHOS  + 1;
      if (clu2->IsEMCAL()) 
        multEMCAL =  multEMCAL + 1;
   }

   //Printf("There are %d caloclusters in this event: %d in PHOS, %d in EMCAL", 
   //fEvent->GetNumberOfCaloClusters(), multPHOS, multEMCAL);

   FillHistogram("hPHOSvsEMCAL", 0.5, multPHOS );
   FillHistogram("hPHOSvsEMCAL", 1.5, multEMCAL);
}

//===========================================================================//
Double_t  AliAnalysisTaskGammaPHOSPP::NonlinearCorrection(Double_t en)
{
//   Double_t en = clu->E();
   if (!fMCArray) 
      return (en);
   else
/*
   if (fEvent->GetRunNumber() > 110000 && fEvent->GetRunNumber() < 150000) { // 7 TeV
      Double_t calib=1.008;
      Double_t ParA=0.015;
      Double_t ParB=0.4;
      return (0.0241+1.0504*en+0.000249*en*en)*calib*(1+ParA/(1.+en*en/ParB/ParB)) ;
      //return clu->E() * calib * (1+ParA/(1+pow(clu->E()/ParB, 2)));
   } else
*/
   if (fLHCRunN == 2) { // Run 2
     Double_t calib = 0.9; //!!! 1st iteration!!!
     Double_t ParA = 0.132081;
     Double_t ParB = -4.91409e-03;

     return en*calib*(1+ParA*TMath::Exp(ParB*en)); 
   } else if (fLHCRunN == 1) { // Run 1
       Double_t calib=1.008;
       Double_t ParA=0.015;
       Double_t ParB=0.4;
       return (1+ParA/(1+pow(en/ParB,2)))*en*calib;
   } else   
   return (en);
}
//===========================================================================//
void AliAnalysisTaskGammaPHOSPP::GammaEfficiencies()
{
  if (!fMCArray) {
    return;
  }
   
  printf("Calculating photon detection efficiencies!!!\n");
  
  TH2F *h2 = (TH2F*)fOutputContainer2->FindObject("hGammaMC_true");
  TH1D *htot = (TH1D*)h2->ProjectionX("hgamma", 70, 170);
  htot->SetTitle("total #gamma-s");
  
  std::vector<std::pair<TString, TString>> cuts;
  cuts.emplace_back("all",  "no cuts");
  cuts.emplace_back("cpv",  "cpv");
  cuts.emplace_back("disp", "disp");
  cuts.emplace_back("both", "cpv+disp");

  TEfficiency *pEff;

  for (auto cut : fPidCuts) {
    TH2F *h2_1 = (TH2F*)fOutputContainer2->FindObject(Form("hCaloPhotonPdgvsPt_%s", cut.first.Data()));
    TH1D *hpass = (TH1D*)h2_1->ProjectionX(Form("hpass_%s", cut.first.Data()), 4000 + 22 + 1, 4000 + 22 + 1);
    hpass->SetTitle(Form("Passed #gamma-s, %s", cut.second.Data()));
    pEff = new TEfficiency(*hpass, *htot);
    pEff->SetTitle(Form("#gamma efficienty, %s", cut.second.Data()));
    pEff->SetName(Form("gamma_efficiency_%s", cut.first.Data()));
    fOutputContainer2->Add(pEff);
    fOutputContainer2->Add(hpass);
  }
  fOutputContainer2->Add(htot);
}

//===========================================================================//
void AliAnalysisTaskGammaPHOSPP::CutEfficiencies() 
{

  printf("Calculating cut efficiencies for clusters...\n");

  std::vector<std::pair<TString, TString>> cuts;
//  cuts.emplace_back("all",  "no cuts");
  cuts.emplace_back("cpv",  "cpv/all");
  cuts.emplace_back("disp", "disp/all");
  cuts.emplace_back("both", "(cpv + disp)/all");

  TEfficiency *pEff;
  
  TH1F *htot = (TH1F*)fOutputContainer->FindObject("hCaloPhotonPt_all");
   
  for (auto cut : fPidCuts) {
    if (cut.first.CompareTo("all")==0) continue;
    TH1F *hpass = (TH1F*)fOutputContainer->FindObject(Form("hCaloPhotonPt_%s", cut.first.Data()));
    pEff = new TEfficiency(*hpass, *htot);
    pEff->SetName(Form("cut_efficiency_%s", cut.first.Data())); 
    fOutputContainer2->Add(pEff);
  }
}

//===========================================================================//
void AliAnalysisTaskGammaPHOSPP::AddClusterHistograms()
{
  const Int_t nPt      = 400;
  const Double_t ptMin = 0;
  const Double_t ptMax = 40;

  for (auto cut : fPidCuts) {
    fOutputContainer->Add(new TH1F(Form("hCaloPhotonPt_%s", cut.first.Data()),"Cluster P_{t} spectrum", nPt, ptMin, ptMax));
    for (Int_t imod = 1; imod <5; imod++) {
      fOutputContainer->Add(new TH1F(Form("hCaloPhotonPt_%s_M%d", cut.first.Data(), imod), Form("Cluster P_{t} spectrum, M%d", imod), nPt, ptMin, ptMax));
    }
    fOutputContainer->Add(new TH2F(Form("hCaloPhotonPtvsNcl_%s", cut.first.Data()),"Minumum number of cells in a cluster vs cluster P_{t}", nPt, ptMin, ptMax, 5, 3, 8));
    fOutputContainer->Add(new TH2F(Form("hCentralityvsClustPt_%s", cut.first.Data()),"Centrality vs ClustPt", nPt, ptMin, ptMax, 8, 0., 8.));
  }

  for (auto cut : fPidCuts) {
    fOutputContainer2->Add(new TH2F(Form("hCaloPhotonPdgvsPt_%s", cut.first.Data()),"Cluster pdg vs p_{T}." , nPt, ptMin, ptMax, 8000, -4000, 4000));
    fOutputContainer2->Add(new TH2F(Form("hCaloPhotonPdgvsPt_%s_naive", cut.first.Data()),"Cluster pdg vs p_{T} (naive)." , nPt, ptMin, ptMax, 8000, -4000, 4000));
    fOutputContainer2->Add(new TH2F(Form("hMatrixEff_%s", cut.first.Data()), "Efficiency matrix", nPt, ptMin, ptMax, nPt, ptMin, ptMax));
    fOutputContainer2->Add(new TH2F(Form("hMatrixEff_gamma_%s", cut.first.Data()), "Efficiency matrix for #gamma", nPt, ptMin, ptMax, nPt, ptMin, ptMax));     
    fOutputContainer2->Add(new TH2F(Form("hMatrixEff_beta_%s", cut.first.Data()), "Efficiency matrix for #beta^{#pm}", nPt, ptMin, ptMax, nPt, ptMin, ptMax));     
  }

}

//===========================================================================//
void AliAnalysisTaskGammaPHOSPP::AddMassHistograms() 
{
  const Int_t nM       = 750;
  const Double_t mMin  = 0.0;
  const Double_t mMax  = 1.5;
  const Int_t nPt      = 400;
  const Double_t ptMin = 0;
  const Double_t ptMax = 40;

  fOutputContainer->Add(new TH1F("hAsym",     "Asymmetry, abs((p1-p2)/(p1+p2))", 100, 0, 1));
  fOutputContainer->Add(new TH1F("hAsym_mix", "Asymmetry, abs((p1-p2)/(p1+p2))", 100, 0, 1));

  for (auto cut : fPidCuts) {
    fOutputContainer->Add(new TH2F(Form("hMassPt_%s", cut.first.Data()) ,"(M,p_{T})_{#gamma#gamma}"   ,nM, mMin, mMax, nPt, ptMin, ptMax));
    fOutputContainer->Add(new TH2F(Form("hMiMassPt_%s", cut.first.Data()) ,   "(M,p_{T})_{#gamma#gamma}"   ,nM, mMin, mMax, nPt, ptMin, ptMax));
    fOutputContainer->Add(new TH3F(Form("hMassPt_asym_%s", cut.first.Data()), "(M,p_{T})_{#gamma#gamma}" ,nM, mMin, mMax, nPt, ptMin, ptMax, 10, 0., 1.));
    fOutputContainer->Add(new TH3F(Form("hMiMassPt_asym_%s", cut.first.Data()), "(M,p_{T})_{#gamma#gamma}" ,nM, mMin, mMax, nPt, ptMin, ptMax, 10, 0., 1.));

    fOutputContainer->Add(new TH2F(Form("hMassSingle_%s", cut.first.Data()),  "(M,p_{T})_{#gamma#gamma}" ,nM, mMin, mMax, nPt, ptMin, ptMax));
    fOutputContainer->Add(new TH2F(Form("hMiMassSingle_%s", cut.first.Data()),  "(M,p_{T})_{#gamma#gamma}" ,nM, mMin, mMax, nPt, ptMin, ptMax));
    fOutputContainer->Add(new TH3F(Form("hMassSingle_asym_%s", cut.first.Data()), "(M,p_{T})_{#gamma#gamma}" ,nM, mMin, mMax, nPt, ptMin, ptMax, 10, 0., 1.));
    fOutputContainer->Add(new TH3F(Form("hMiMassSingle_asym_%s", cut.first.Data()), "(M,p_{T})_{#gamma#gamma}" ,nM, mMin, mMax, nPt, ptMin, ptMax, 10, 0., 1.));

    for (Int_t iMod = 1; iMod < 5; iMod ++) {  
      fOutputContainer->Add(new TH2F(Form("hMassPt_%s_M%d",       cut.first.Data(), iMod), Form("(M,p_{T})_{#gamma#gamma}, M%d",  iMod), nM, mMin, mMax, nPt, ptMin, ptMax));
      fOutputContainer->Add(new TH2F(Form("hMiMassPt_%s_M%d",     cut.first.Data(), iMod), Form("(M,p_{T})_{#gamma#gamma}, M%d",  iMod), nM, mMin, mMax, nPt, ptMin, ptMax));
      fOutputContainer->Add(new TH2F(Form("hMassSingle_%s_M%d",   cut.first.Data(), iMod), Form("(M,p_{T})_{#gamma#gamma}, M%d",  iMod), nM, mMin, mMax, nPt, ptMin, ptMax));
      fOutputContainer->Add(new TH2F(Form("hMiMassSingle_%s_M%d", cut.first.Data(), iMod), Form("(M,p_{T})_{#gamma#gamma}, M%d",  iMod), nM, mMin, mMax, nPt, ptMin, ptMax));
    }
  }
}

//===========================================================================//
void AliAnalysisTaskGammaPHOSPP::AddMCHistograms() 
{
  const Int_t nPt      = 400;
  const Double_t ptMin = 0;
  const Double_t ptMax = 40;

  fOutputContainer2->Add(new TH1F("hEventCounterMC","Count events", 1, 0, 1));
  
  fOutputContainer2->Add(new TH1F("hWeights","Particle weights", 2000, 0., 2.));

  fOutputContainer2->Add(new TH2F("hMCPdgvsPt","Pdg code vs particle p_{T}", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH1F("hMCPt","MC particles p_{T}", nPt, ptMin, ptMax));

  fOutputContainer2->Add(new TH2F("hPi0MC", "MC distribution of #pi^{0}-s ", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new TH2F("hPiPlusMC", "MC distribution of #pi^{+}-s ", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new TH2F("hPiMinusMC", "MC distribution of #pi^{-}-s ", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hEtaMC", "MC distribution of #eta^{0}-s ", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hEtaPrimeMC", "MC distribution of #eta'", nPt, ptMin, ptMax,240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hOmegaMC", "MC distribution of #omega", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hK0SMC", "MC distribution of #K_{0}^{S}-s", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hK0LMC", "MC distribution of #K_{0}^{L}-s", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hProtonMC", "MC distribution of protons", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hNeutrMC", "MC distribution of neutrons", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hKchMC", "MC distribution of neutrons", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hPdgvsPt_MCTracks","MC Particle PDG code vs Pt", nPt, ptMin, ptMax , 8000, -4000, 4000));

  fOutputContainer2->Add(new  TH2F("hGammaMC_all", "MC distribution of #gamma-s", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hGammaMC_true", "MC distribution of #gamma-s at IP", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  
  fOutputContainer2->Add(new  TH2F("hGammaMCSources", "MC distribution of #gamma-s", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new  TH2F("hCaloPhotonPdgvsPt_FromMaterial", "Cluster pdg vs p_{T} from material" , nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new  TH2F("hGammaMC_FromMaterial", "MC distribution of #gamma-s", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hBetaMC_FromMaterial", "MC distribution of #gamma-s", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  
  fOutputContainer2->Add(new TH1I("hpid3", "Pid #sigma<3, #pi, p, K, #beta, Undef", 5, 0., 5.));
  fOutputContainer2->Add(new TH1I("hpid1", "Pid #sigma<1, #pi, p, K, #beta, Undef", 5, 0., 5.));
}

//============================================================
void AliAnalysisTaskGammaPHOSPP::AddTrackHistograms() 
{
  
  Int_t nPt      = 400;
  Double_t ptMin = 0;
  Double_t ptMax = 40;
  
  
  fOutputContainer->Add(new TH1F("hDistanceToBadChannel", "Distance to bad channel in cm", 1000, 0, 100));

  fOutputContainer ->Add(new TH1F("hTrackCharge","Charge of track", 9, -4.5, 4.5));
  fOutputContainer2->Add(new TH1F("hMCTrackCharge","Charge of MC track", 9, -4.5, 4.5));
 
  fOutputContainer->Add(new TH2F("hTracksPt2D","Tracks ID vs p_{T}", nPt, ptMin, ptMax, 9, 0.,9.));
  fOutputContainer->Add(new TH1F("hTracksPt_beta","Electron tracks p_{T}", nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH1F("hTracksPt_mu","#mu tracks p_{T}", nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH1F("hTracksPt_pi","#pi^{#pm} tracks p_{T}", nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH1F("hTracksPt_ka","K^{#pm} tracks p_{T}", nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH1F("hTracksPt_p","p, #bar{p} tracks vs p_{T}", nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH1F("hTracksPt_de","^{2}_{1}H tracks vs p_{T}", nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH1F("hTracksPt_tri","^{3}_{1}H tracks vs p_{T}", nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH1F("hTracksPt_he3","^{3}_{2}He tracks vs p_{T}", nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH1F("hTracksPt_alpha","^{4}_{2} He tracks vs p_{T}", nPt, ptMin, ptMax));

  fOutputContainer->Add(new TH2F("hTracksofClusts","Tracks vc p_{T}.", nPt, ptMin, ptMax, 9, 0.,9.));
  fOutputContainer->Add(new TH2F("hTracksofClusts_cpv","Tracks vc p_{T}, CPV cut.", nPt, ptMin, ptMax, 9, 0.,9.));
  fOutputContainer->Add(new TH2F("hTracksofClusts_disp","Tracks vc p_{T}, shape cut.", nPt, ptMin, ptMax, 9, 0.,9.));
  fOutputContainer->Add(new TH2F("hTracksofClusts_both","Tracks vc p_{T}, CPV + shape cuts." , nPt, ptMin, ptMax, 9, 0.,9.));

  fOutputContainer->Add(new TH2F("hTracksOfPi_OneSigma", "#pi^{+} clusters", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfPi_OneSigma_disp", "#pi^{+} clusters, Disp cut", nPt, ptMin, ptMax, 1000., 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiPi_OneSigma", "#pi^{-} clusters", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiPi_OneSigma_disp", "#pi^{-} clusters, Disp cut", nPt, ptMin, ptMax, 1000., 0., 100.));

  fOutputContainer->Add(new TH2F("hTracksOfPi_ThreeSigma", "#pi^{+} clusters", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfPi_ThreeSigma_disp", "#pi^{+} clusters, Disp cut", nPt, ptMin, ptMax, 1000., 0., 100.));  
  fOutputContainer->Add(new TH2F("hTracksOfAntiPi_ThreeSigma", "#pi^{-} clusters", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiPi_ThreeSigma_disp", "#pi^{-} clusters, Disp cut", nPt, ptMin, ptMax, 1000., 0., 100.));

  fOutputContainer->Add(new TH2F("hTracksOfPr_OneSigma", "p clusters", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfPr_OneSigma_disp", "p clusters, Disp cut", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiPr_OneSigma", "#bar{p} clusters", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiPr_OneSigma_disp", "#bar{p} clusters, Disp cut", nPt, ptMin, ptMax, 1000, 0., 100.));

  fOutputContainer->Add(new TH2F("hTracksOfPr_ThreeSigma", "p clusters", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfPr_ThreeSigma_disp", "p clusters, Disp cut", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiPr_ThreeSigma", "#bar{p} clusters", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiPr_ThreeSigma_disp", "#bar{p} clusters, Disp cut", nPt, ptMin, ptMax, 1000, 0., 100.));
 
  fOutputContainer->Add(new TH2F("hTracksOfKa_OneSigma", "K^{+} clusters", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfKa_OneSigma_disp", "K^{+} clusters, Disp cut", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiKa_OneSigma", "K^{-} clusters", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiKa_OneSigma_disp", "K^{-} clusters, Disp cut", nPt, ptMin, ptMax, 1000, 0., 100.));

  fOutputContainer->Add(new TH2F("hTracksOfKa_ThreeSigma", "K^{+} clusters", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfKa_ThreeSigma_disp", "K^{+} clusters, Disp cut", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiKa_ThreeSigma", "K^{-} clusters", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiKa_ThreeSigma_disp", "K^{-} clusters, Disp cut", nPt, ptMin, ptMax, 1000, 0., 100.));

  fOutputContainer->Add(new TH2F("hTracksOfBeta_OneSigma", "#beta^{-} clusters", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfBeta_OneSigma_disp", "#beta^{-} clusters, Disp cut", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiBeta_OneSigma", "#beta^{+} clusters", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiBeta_OneSigma_disp", "#beta^{+} clusters, Disp cut", nPt, ptMin, ptMax, 1000, 0., 100.));

  fOutputContainer->Add(new TH2F("hTracksOfBeta_ThreeSigma", "#beta^{-} clusters", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfBeta_ThreeSigma_disp", "#beta^{-} clusters, Disp cut", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiBeta_ThreeSigma", "#beta^{+} clusters", nPt, ptMin, ptMax, 1000, 0., 100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiBeta_ThreeSigma_disp", "#beta^{+} clusters, Disp cut", nPt, ptMin, ptMax, 1000, 0., 100.));

  fOutputContainer->Add(new TH1F("hTracksOfOthers", "Other clusters", nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hTracksOfOthers_disp", "Other clusters, Disp cut", nPt, ptMin, ptMax, 1000, 0., 100.));

  //==== MC ====

  fOutputContainer2->Add(new TH2F("hTracksOfPi_OneSigma_label", "PDG #pi^{+}-s vs p_{T}", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfPi_OneSigma_disp_label", "PDG #pi^{+}-s vs p_{T}, dispersion cut", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiPi_OneSigma_label", "PDG #pi^{-}-s vs p_{T}", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiPi_OneSigma_disp_label", "PDG #pi^{-}-s vs p_{T}, dispersion cut", nPt, ptMin, ptMax, 8000, -4000, 4000));

  fOutputContainer2->Add(new TH2F("hTracksOfPi_ThreeSigma_label", "PDG #pi^{+}-s vs p_{T}", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfPi_ThreeSigma_disp_label", "PDG #pi^{+}-s vs p_{T}, dispersion cut", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiPi_ThreeSigma_label", "PDG #pi^{-}-s vs p_{T}", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiPi_ThreeSigma_disp_label", "PDG #pi^{-}-s vs p_{T}, dispersion cut", nPt, ptMin, ptMax, 8000, -4000, 4000));

  fOutputContainer2->Add(new TH2F("hTracksOfPr_OneSigma_label", "PDG p-s vs p_{T}", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfPr_OneSigma_disp_label", "PDG p-s vs p_{T}, dispersion cut", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiPr_OneSigma_label", "PDG #bar{p}-s vs p_{T}", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiPr_OneSigma_disp_label", "PDG #bar{p}-s vs p_{T}, dispersion cut", nPt, ptMin, ptMax, 8000, -4000, 4000));

  fOutputContainer2->Add(new TH2F("hTracksOfPr_ThreeSigma_label", "PDG p-s vs p_{T}", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfPr_ThreeSigma_disp_label", "PDG p-s vs p_{T}, dispersion cut", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiPr_ThreeSigma_label", "PDG #bar{p}-s vs p_{T}", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiPr_ThreeSigma_disp_label", "PDG #bar{p}-s vs p_{T}, dispersion cut", nPt, ptMin, ptMax, 8000, -4000, 4000));

  fOutputContainer2->Add(new TH2F("hTracksOfKa_OneSigma_label", "PDG K^{+}-s vs p_{T}", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfKa_OneSigma_disp_label", "PDG K^{+}-s vs p_{T}, dispersion cut", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiKa_OneSigma_label", "PDG K^{-}-s vs p_{T}", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiKa_OneSigma_disp_label", "PDG K^{-}-s vs p_{T}, dispersion cut", nPt, ptMin, ptMax, 8000, -4000, 4000));

  fOutputContainer2->Add(new TH2F("hTracksOfKa_ThreeSigma_label", "PDG K^{+}-s vs p_{T}", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfKa_ThreeSigma_disp_label", "PDG Ka^{+}-s vs p_{T}, dispersion cut", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiKa_ThreeSigma_label", "PDG K^{-}-s vs p_{T}", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiKa_ThreeSigma_disp_label", "PDG K^{-}-s vs p_{T}, dispersion cut", nPt, ptMin, ptMax, 8000, -4000, 4000));

  fOutputContainer2->Add(new TH2F("hTracksOfBeta_OneSigma_label", "PDG #beta^{-}-s vs p_{T}", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfBeta_OneSigma_disp_label", "PDG #beta^{-}-s vs p_{T}, dispersion cut", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiBeta_OneSigma_label", "PDG #beta^{+}-s vs p_{T}", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiBeta_OneSigma_disp_label", "PDG #beta^{+}-s vs p_{T}, dispersion cut", nPt, ptMin, ptMax, 8000, -4000, 4000));

  fOutputContainer2->Add(new TH2F("hTracksOfBeta_ThreeSigma_label", "PDG #beta^{-}-s vs p_{T}", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfBeta_ThreeSigma_disp_label", "PDG #beta^{-}-s vs p_{T}, dispersion cut", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiBeta_ThreeSigma_label", "PDG #beta^{+}-s vs p_{T}", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiBeta_ThreeSigma_disp_label", "PDG #beta^{+}-s vs p_{T}, dispersion cut", nPt, ptMin, ptMax, 8000, -4000, 400));

//=================
  fOutputContainer2->Add(new TH1F("hTracksOfPiClose","Tracks of #pi^{#pm} close to cluster", nPt, ptMin, ptMax));
  fOutputContainer2->Add(new TH1F("hTracksOfPiCloseDispOK","Tracks of #pi^{#pm} close to cluster", nPt, ptMin, ptMax));
  fOutputContainer2->Add(new TH1F("hTracksOfPrClose","Tracks of p#bar{p} close to cluster", nPt, ptMin, ptMax));
  fOutputContainer2->Add(new TH1F("hTracksOfPrCloseDispOK","Tracks of p#bar{p} close to cluster", nPt, ptMin, ptMax));
  fOutputContainer2->Add(new TH1F("hTracksOfKaClose","Tracks of K^{#pm} close to cluster", nPt, ptMin, ptMax));
  fOutputContainer2->Add(new TH1F("hTracksOfKaCloseDispOK","Tracks of Ka^{#pm} close to cluster", nPt, ptMin, ptMax));
  fOutputContainer2->Add(new TH1F("hTracksOfOthersClose","Tracks of K^{#pm} close to cluster", nPt, ptMin, ptMax));
  fOutputContainer2->Add(new TH1F("hTracksOfOthersCloseDispOK","Tracks of Ka^{#pm} close to cluster", nPt, ptMin, ptMax));
        
    
  fOutputContainer2->Add(new TH1F("hTracks_matched","Nr of matched tracks per cluster",20, 0,20));
  fOutputContainer2->Add(new TH2F("hEnTrackvsClust","Energies of track vs cluster", nPt, ptMin, ptMax, nPt, ptMin, ptMax));
  fOutputContainer2->Add(new TH1F("hEmcCPVDistance","EMC to CPV distance", 100, 0., 100.));
  fOutputContainer2->Add(new TH1F("hDistance","Distance from cluster to associated track, cm", 100, 0., 100.));
}
//================================================================================================
//
void AliAnalysisTaskGammaPHOSPP::AddQAHistograms()
{
  const Int_t nPt      = 400;
  const Double_t ptMin = 0;
  const Double_t ptMax = 40;

  fOutputContainer ->Add(new TH1F("hEventCounter","Count accepted events", 1, 0, 1));
  fOutputContainer ->Add(new TH1F("hEventCounterCentrality","Count events in centrality bins", 8, 0. , 8.));
  fOutputContainer ->Add(new TH1F("hPHOSCounter","Count accepted photons", 2, 0, 2));
 
  fOutputContainer->Add(new TH1I("hCellMultEvent"  ,"PHOS cell multiplicity per event"    ,2000, 0,2000));
  fOutputContainer->Add(new TH1I("hClusterMult"      ,"CaloCluster multiplicity"     , 100, 0, 100));
  fOutputContainer->Add(new TH1I("hPHOSClusterMult"  ,"PHOS cluster multiplicity"    , 100, 0, 100));
  fOutputContainer->Add(new TH1F("hCellEnergy"  ,"Cell energy"            ,5000, 0.,50.));
  fOutputContainer->Add(new TH1F("hClusterEnergy"  ,"Cluster energy"      ,5000, 0.,50.));

  for (Int_t imod = 1; imod < 5; imod ++) {
    fOutputContainer->Add(new TH1I(Form("hCellMultEventM%d", imod),Form("PHOS cell multiplicity per event, M%id", imod), 2000, 0, 2000));
    fOutputContainer->Add(new TH1I(Form("hPHOSClusterMultM%d",imod),Form("PHOS cluster multiplicity, M%d",imod), 100, 0, 100));
    fOutputContainer->Add(new TH1F(Form("hCellEnergyM%d",imod),Form("Cell energy in module %d", imod),5000, 0.,50.));
    fOutputContainer->Add(new TH1F(Form("hClusterEnergyM%d", imod),Form("Cluster energy, M%d",imod)  ,5000, 0.,50.));
    fOutputContainer->Add(new TH2F(Form("hCellNXZM%d", imod), Form("Cell (X,Z), M%d" , imod), 64, 0.5, 64.5, 56, 0.5,56.5));
    fOutputContainer->Add(new TH2F(Form("hCellEXZM%d", imod), Form("Cell E(X,Z), M%d", imod), 64, 0.5, 64.5, 56, 0.5,56.5));
    fOutputContainer->Add(new TH2F(Form("hCluNXZM%d",  imod), Form("Clu (X,Z), M%d",   imod), 64, 0.5, 64.5, 56, 0.5,56.5));
    fOutputContainer->Add(new TH2F(Form("hCluEXZM%d",  imod), Form("Clu E(X,Z), M%d",  imod), 64, 0.5, 64.5, 56, 0.5,56.5));
  }

  for (auto cut : fPidCuts) {
    fOutputContainer->Add(new TH2F(Form("hClusterEvsN_%s", cut.first.Data())  ,"Cluster energy vs digit multiplicity"    ,5000, 0.,50., 40, 0., 40.));
    for (Int_t imod = 1; imod < 5; imod ++) {
      fOutputContainer->Add(new TH2F(Form("hClusterEvsN_%s_M%d", cut.first.Data(), imod),Form("Cluster energy vs digit multiplicity, M%d",imod),5000, 0.,50., 40, 0., 40.));
  } 
    fOutputContainer->Add(new TH1I(Form("hCellMultClu_%s", cut.first.Data())  ,"Cell multiplicity per cluster"    ,200, 0,200));
    for (Int_t imod = 1; imod < 5; imod ++) {
      fOutputContainer->Add(new TH1I(Form("hCellMultClu_%s_M%d", cut.first.Data(), imod),Form("Cell multiplicity per cluster, M%d", imod),200, 0,200));
    }
  }

  fOutputContainer->Add(new TH1I("hModule","Module events", 5, 0., 5.));
  fOutputContainer->Add(new TH1F("hSelEvents","Selected events",9, -0.5, 8.5));
  fOutputContainer->Add(new TH1F("hSelEventspi0","Selected events for pi0", 9, -0.5, 8.5)); 

  fOutputContainer->Add(new TH1F("hPhotonKappa","#kappa(#gamma)", 400, 0., 40.));
  fOutputContainer->Add(new TH1F("hPhotonPt","p_{T}(#gamma)", 400, 0., 40.));
  fOutputContainer->Add(new TH1F("hPhotonPx","p_{x}(#gamma)", 400, 0., 40.));
  fOutputContainer->Add(new TH1F("hPhotonPy","p_{y}(#gamma)", 400, 0., 40.));

  fOutputContainer->Add(new TH1F("hTrigClass","Trigger class",5, 0.5,5.5));
  fOutputContainer->Add(new TH1F("hNPileupVtx","Number of SPD pileup vertices", 10, 0., 10.));
  fOutputContainer->Add(new TH1F("hZPileupVtx", "Z pileup", 200, -50., 50.));

  fOutputContainer->Add(new TH1F("hZvertex","Z vertex",200, -50.,+50.));
  fOutputContainer->Add(new TH1F("hNvertexTracks","N of primary tracks from the primary vertex", 150, 0., 150.));
  fOutputContainer->Add(new TH1F("hTrackMult","Charged track multiplicity", 150, 0., 150.));

  fOutputContainer->Add(new TH1F("hV0Atime","V0A time", 1200, -6.e-6,+6.e-6));
  fOutputContainer->Add(new TH1F("hV0Ctime","V0C time", 1200, -6.e-6,+6.e-6));
  fOutputContainer->Add(new TH2F("hV0AV0Ctime","V0A time vs V0C time", 120, -6.e-6,+6.e-6 , 120, -6.e-6,+6.e-6));

  fOutputContainer->Add(new TH1D("hTOF","Time of Flight", 1000, 0, 1e-6));
  fOutputContainer->Add(new TH1I("hBC", "BC", 1000, 0, 1000));

  fOutputContainer->Add(new TH1F("hvt0vsvt5", "Difference of momenta for different vertices", 200, -1., 1.));

  fOutputContainer->Add(new TH2F("hClustM02","Cluster M02 vs p_{T}", nPt, ptMin, ptMax, 100, 0, 10));
  fOutputContainer->Add(new TH2F("hClustM20","Cluster M20 vs p_{T}", nPt, ptMin, ptMax, 100, 0, 10));
  fOutputContainer ->Add(new TH1F("hPHOSvsEMCAL", "PHOS and EMCAL clusters", 2, 0, 2));
}
