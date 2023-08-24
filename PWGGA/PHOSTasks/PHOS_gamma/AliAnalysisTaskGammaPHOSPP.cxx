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
  fEventVtxExists(0), fEventVtxZ10cm(0), fEventPileup(0), fEventV0AND(0), 
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
  fNsigmaCPV(2.5),
  fNsigmaDisp(2.5)

{
  // define cust
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

  AddOnePhotonHistograms();
  AddTwoPhotonHistograms();
  AddQAHistograms();
  AddTrackHistograms();
  AddMCHistograms();

  PostData(1, fOutputContainer);
  PostData(2, fOutputContainer2);
}

//________________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::UserExec(Option_t *)  // Main loop, called for each event
{
  // Get event
  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fEvent) {
     Printf("ERROR: Could not retrieve evenat!");
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

  //Get MC
  fMCArray = (TClonesArray*)fEvent->FindListObject(AliAODMCParticle::StdBranchName());

  // Filter events
  Bool_t acceptEvent = kFALSE;
  acceptEvent = AcceptEvent(fEvent);
  if (!acceptEvent) return;

  //PHOS event
  if (fPHOSEvent) {
    fPHOSEvent->Clear() ;
  }  
    
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
  ProcessMC() ;
 
  // Notify 
  Notify();

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

//_______________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::FinishTaskOutput() // Terminate
{
  GammaAcceptance();
  GammaEfficiencies();
  CutEfficiencies();
}

//________________________________________________________________________
Bool_t AliAnalysisTaskGammaPHOSPP::AcceptEvent(AliAODEvent *event)
{

  // Event selection flags
  fEventVtxExists    = kFALSE;
  fEventVtxZ10cm     = kFALSE;
  fEventPileup       = kFALSE;
  fEventV0AND         = kFALSE;
  
  FillHistogram("hSelEvents", 0) ; // Events accepted by the physics selection

  TString trigClasses = event->GetFiredTriggerClasses();

  if (trigClasses.Contains("FAST")  && !trigClasses.Contains("ALL")) {
    AliWarning(Form("Skip event with triggers %s", trigClasses.Data()));
    return kFALSE;
  }
  else FillHistogram("hSelEvents", 1) ; // After filtering by trigger classes

  const AliAODVertex *aodVertex5 =    event->GetPrimaryVertex();
  const AliAODVertex *aodVertexSPD  = event->GetPrimaryVertexSPD();
  
  FillHistogram("hNContributors", aodVertex5->GetNContributors());
  
  if (!aodVertex5) {
     fEventVtxExists    = kFALSE;
  }
  else fEventVtxExists = kTRUE;

  if (!fEventVtxExists) return kFALSE;
    else FillHistogram("hSelEvents", 2); //vtx exists
  
  if (aodVertex5 && aodVertex5->GetNContributors() < 1)  return kFALSE;
    else FillHistogram("hSelEvents", 3);  // aod vertex exists and has more than one contributors

  if (aodVertexSPD && aodVertexSPD->GetNContributors() < 1)  return kFALSE;
    else FillHistogram("hSelEvents", 4); // SPDvertex exists and has more than one contributors

  fVtx5[0] = aodVertex5->GetX();
  fVtx5[1] = aodVertex5->GetY();
  fVtx5[2] = aodVertex5->GetZ();

  FillHistogram("hZvertex", aodVertex5->GetZ());
  FillHistogram("hSPDvertex", aodVertexSPD->GetZ());
  FillHistogram("hDistSPDTrackvertex", TMath::Abs(aodVertex5->GetZ() - aodVertexSPD->GetZ()));
  FillHistogram("hTrackvsSPDVertices", aodVertexSPD->GetZ(), aodVertex5->GetZ());
  if (TMath::Abs(aodVertex5->GetZ()) < 10.) fEventVtxZ10cm = kTRUE; 
  if (!fEventVtxZ10cm) return kFALSE;
    else FillHistogram("hSelEvents", 5); //interaction vertex within 10 cm from the center

  if (event->IsPileupFromSPD()) {
    fEventPileup = kTRUE;
    FillHistogram("hNPileupVtx", event->GetNumberOfPileupVerticesSPD());
    for (Int_t puVtx = 0;  puVtx < event->GetNumberOfPileupVerticesSPD(); puVtx++) {
      Double_t dZpileup = aodVertexSPD->GetZ() - event->GetPileupVertexSPD(puVtx)->GetZ();
      FillHistogram("hZPileupVtx", dZpileup);
    }
  }

  if (fEventPileup) return kFALSE;
    else FillHistogram("hSelEvents", 6); //no pileup from SPD

  if (trigClasses.Contains("CINT1B")) fnCINT1B++;
  if (trigClasses.Contains("CINT1A")) fnCINT1A++;
  if (trigClasses.Contains("CINT1C")) fnCINT1C++;
  if (trigClasses.Contains("CINT1-E")) fnCINT1E++;

  // ULong64_t trigmask = event->GetTriggerMask();
  fEventV0AND = fTriggerAnalysis->IsOfflineTriggerFired(event, AliTriggerAnalysis::kV0AND);

  if (!fEventV0AND) return kFALSE;
    else FillHistogram("hSelEvents", 7); // V0AND trigger

  Float_t tV0A = event->GetVZEROData()->GetV0ATime();
  Float_t tV0C = event->GetVZEROData()->GetV0CTime();
  FillHistogram("hV0Atime",    tV0A);
  FillHistogram("hV0Ctime",    tV0C);
  FillHistogram("hV0AV0Ctime", tV0A, tV0C);  

  return kTRUE;
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskGammaPHOSPP::GetEventCentrality(AliAODEvent *event)
{
  //Calculate charged multiplicity   Bool_t   IsCPVOK(void)const {return fCpv;}
  const Int_t trackMult = event->GetNumberOfTracks() ;

  FillHistogram("hTrackMult", trackMult + 0.5) ;

  if (trackMult <= 2)
    return 0 ;
  else  if (trackMult <= 5)
    return 1 ;
  else if (trackMult <= 9)
     return  2 ;
  else if (trackMult <= 14)
      return 3 ;
  else if (trackMult <= 22)
     return 4 ;
  else if (trackMult <= 35)
     return 5 ;
  else if (trackMult <= 50)
     return 6 ;
  else return 7 ;
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

  cluPHOS1.SetE(NonlinearMCCorrection(en));

  cluPHOS1.GetMomentum(p1 , fVtx0);
  cluPHOS1.GetMomentum(p11, fVtx5);
   
  Double_t pAbs = p1.P();
  Double_t pT   = p1.Pt();
  Double_t pX   = p1.Px();
  Double_t pY   = p1.Py();

  if (pAbs < 1.e-10) return ;

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
 
  TestMatchingTrackPID(clu1, p11.Pt());
   
  if (!fMCArray) {
    weight = 1.0;
  }  else {
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
  ph->SetCPVBit(ph->GetNsigmaCPV() > fNsigmaCPV);
  ph->SetDispBit(ph->GetNsigmaFullDisp() < fNsigmaDisp);
  ph->SetBC(TestBC(clu1->GetTOF()));
  ph->SetPrimary(GetPrimaryLabel(clu1));
  ph->SetPrimaryAtVertex(GetPrimaryLabelAtVertex(clu1));
  ph->SetWeight(weight); 
  

  FillHistogram("hvt0vsvt5", p11.Pt()- p1.Pt());
  FillHistogram("hBC", TestBC(clu1->GetTOF()) + 0.5);
  FillHistogram("hTOF", clu1->GetTOF());
  FillHistogram("hWeights", ph->GetWeight());      

  fInPHOS++ ;
}

//___________________________________________________________________________
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

     Double_t enMC   = ((AliAODMCParticle*)fMCArray->At(ph->GetPrimary()))->E();
     Double_t enMeas = ph->E();
     if (pdg_naive == 22 || pdg_naive == 11) 
       FillHistogram("hEnergyResolution_EM", enMC, enMeas - enMC);
     else
       FillHistogram("hEnergyResolution_other", enMC, enMeas - enMC);
  }
   
   Int_t sm1 = ph->Module();

   Double_t pt   = p11.Pt();
   Double_t ptMC = TestGammaPt(ph);

   std::vector<TString> passed_cuts = {"all"};
   Bool_t CPVBit  =  ph->IsCPVOK();
   Bool_t DispBit =  ph->IsDispOK();
   if (CPVBit)  passed_cuts.emplace_back("cpv");
   if (DispBit) passed_cuts.emplace_back("disp");
   if (CPVBit && DispBit) passed_cuts.emplace_back("both");

   for (auto cut : passed_cuts) {
     if (pdg == 22) {
        FillHistogram(Form("hMatrixEff_%s", cut.Data()), pt, ptMC);   
        if (pdg_naive == 22)     
          FillHistogram(Form("hMatrixEff_gamma_%s", cut.Data()), pt, ptMC); 
        if (TMath::Abs(pdg_naive) == 11)     
          FillHistogram(Form("hMatrixEff_beta_%s", cut.Data()), pt, ptMC); 
     } 

     FillHistogram(Form("hCaloPhotonPt_%s", cut.Data()), pt, weight);
     FillHistogram(Form("hCaloPhotonPdgvsPt_%s", cut.Data()), pt, pdg, weight);
     FillHistogram(Form("hCaloPhotonPdgvsPt_%s_naive", cut.Data()),  pt, pdg_naive, weight);

     FillHistogram(Form("hCaloPhotonPtvsNcl_%s", cut.Data()), pt, ph->GetNCells() + 0.5);
     FillHistogram(Form("hCentralityvsClustPt_%s", cut.Data()), pt, fEventCentrality + 0.5);	
     FillHistogram(Form("hCaloPhotonPt_%s_M%d", cut.Data(), sm1),  pt, weight);
   }
/*
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
  */ 
}

//____________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::FillTwoPhotonHistograms()
{
  TLorentzVector p1, p2, p12, pv1, pv2, pv12, p11;

  for (Int_t i1 = 0; i1 < fInPHOS-1; i1 ++ ) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    if (fLHCRunN == 1 && ph1->GetBC() != 0 ) continue; //Run 1
    pv1 = *(ph1->GetMomV2());
    Double_t P1 = pv1.P();
    Double_t pt1 = ph1->Pt();
    Int_t sm1 = ph1->Module();

    for (Int_t i2 = i1+1; i2 < fInPHOS; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent->At(i2) ;
      if (fLHCRunN == 1 && ph2->GetBC() != 0 ) continue; //Run 1
      pv2 = *(ph2->GetMomV2());
      Double_t P2 = pv2.P();
      Double_t pt2 = ph2->Pt();
      Int_t sm2 = ph2->Module();

      Double_t asym = TMath::Abs((P1-P2)/(P1+P2));
      FillHistogram("hAsym", asym);
      p12  = *ph1  + *ph2;
      pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
      Double_t ma12 = pv12.M();
      Double_t pt12 = pv12.Pt();

      std::vector<TString> passed_cuts = {"all"};
      Bool_t CPVBit = ph1->IsCPVOK() && ph2->IsCPVOK();
      Bool_t DispBit = ph1->IsDispOK() && ph2->IsDispOK();
      if (CPVBit) passed_cuts.emplace_back("cpv");
      if (DispBit) passed_cuts.emplace_back("disp");
      if (CPVBit && DispBit) passed_cuts.emplace_back("both");

      for (auto cut : passed_cuts) {
        FillHistogram(Form("hMassPt_%s", cut.Data()), ma12 , pt12 );
        FillHistogram(Form("hMassPt_asym_%s", cut.Data()), ma12, pt12, asym);
        if (sm1 == sm2) FillHistogram(Form("hMassPt_%s_M%d", cut.Data(), sm1), ma12, pt12);
        FillHistogram(Form("hMassSingle_%s", cut.Data()), ma12, pt1) ;
        FillHistogram(Form("hMassSingle_%s", cut.Data()), ma12, pt2) ;
        FillHistogram(Form("hMassSingle_asym_%s", cut.Data()), ma12, pt1, asym);
        FillHistogram(Form("hMassSingle_asym_%s", cut.Data()), ma12, pt2, asym);
        FillHistogram(Form("hMassSingle_%s_M%d", cut.Data(), sm1), ma12, pt1);
        FillHistogram(Form("hMassSingle_%s_M%d", cut.Data(), sm2), ma12, pt2);
      }
  } //end of loop  i2 
 } //end of loop   i1
}
//______________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::MixPhotons()
{
  TLorentzVector p1, p2, p12, pv1, pv2, pv12, p11;

  Int_t zvtx = (Int_t)((fVtx5[2] + 10.)/2.) ;
  if(zvtx < 0) zvtx = 0 ;
  if(zvtx > 9) zvtx = 9 ;
  
  Int_t centr = 0;
  
  if(!fPHOSEvents[zvtx][centr]) fPHOSEvents[zvtx][centr]=new TList() ;
  
  TList * prevPHOS = fPHOSEvents[zvtx][centr] ;

  for (Int_t i1 = 0; i1 < fInPHOS-1; i1 ++ ) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    if (fLHCRunN == 1 && ph1->GetBC() != 0 ) continue; //Run 1
    pv1 = *(ph1->GetMomV2());
    Double_t P1 = pv1.P();
    Double_t pt1 = ph1->Pt();
    Int_t sm1 = ph1->Module();

    for (Int_t ev = 0; ev < prevPHOS->GetSize(); ev ++) {
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev)) ;

      for (Int_t i2 = 0; i2 < mixPHOS->GetEntriesFast(); i2 ++) {
        AliCaloPhoton * ph2=(AliCaloPhoton*)mixPHOS->At(i2) ;
        if (fLHCRunN == 1 && ph2->GetBC() != 0 ) continue; //Run 1
        pv2 = *(ph2->GetMomV2());
        Double_t P2 = pv2.P();
        Double_t pt2 = ph2->Pt();
        Int_t sm2 = ph2->Module();

        Double_t asym = TMath::Abs((P1-P2)/(P1+P2));
        FillHistogram("hAsym", asym);
        p12  = *ph1  + *ph2;
        pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
        Double_t ma12 = pv12.M();
        Double_t pt12 = pv12.Pt();

        std::vector<TString> passed_cuts = {"all"};
        Bool_t cpvBit1 = ph1->IsCPVOK(), cpvBit2 = ph1->IsCPVOK();
        Bool_t dispBit1  = ph1->IsDispOK(), dispBit2 = ph2->IsDispOK();
        if (cpvBit1 && cpvBit2) passed_cuts.emplace_back("cpv");
        if (dispBit1 && dispBit2) passed_cuts.emplace_back("disp");
        if ((cpvBit1 && cpvBit2) && (dispBit1 && dispBit2)) passed_cuts.emplace_back("both");
        
        if (cpvBit1 || cpvBit2) passed_cuts.emplace_back("all_cpv");
        if (dispBit1 || dispBit2) passed_cuts.emplace_back("all_disp");

        if (!cpvBit1 || !cpvBit2) passed_cuts.emplace_back("all_anticpv");
        if (!dispBit1 || !dispBit2) passed_cuts.emplace_back("all_antidisp");

        if ((cpvBit1 && !dispBit2) || (!dispBit1 && cpvBit2)) passed_cuts.emplace_back("cpv_antidisp");
        if ((!cpvBit1 && dispBit2) || (dispBit1 && !cpvBit2)) passed_cuts.emplace_back("anticpv_disp");

        for (auto cut : passed_cuts) {
          FillHistogram(Form("hMiMassPt_%s", cut.Data()), ma12 , pt12 );
          FillHistogram(Form("hMiMassPt_asym_%s", cut.Data()), ma12, pt12, asym);
          if (sm1 == sm2) FillHistogram(Form("hMiMassPt_%s_M%d", cut.Data(), sm1), ma12, pt12);
          FillHistogram(Form("hMiMassSingle_%s", cut.Data()), ma12, pt1) ;
          FillHistogram(Form("hMiMassSingle_%s", cut.Data()), ma12, pt2) ;
          FillHistogram(Form("hMiMassSingle_asym_%s", cut.Data()), ma12, pt1, asym);
          FillHistogram(Form("hMiMassSingle_asym_%s", cut.Data()), ma12, pt2, asym);
          FillHistogram(Form("hMiMassSingle_%s_M%d", cut.Data(), sm1), ma12, pt1);
          FillHistogram(Form("hMiMassSingle_%s_M%d", cut.Data(), sm2), ma12, pt2);
        }
      }
    }
  }  

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
//_____________________________________________________________________________
Int_t AliAnalysisTaskGammaPHOSPP::GetPrimaryLabelAtVertex(AliVCluster *clu) //Returns label at vertex
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

//_____________________________________________________________________________
Int_t AliAnalysisTaskGammaPHOSPP::GetPrimaryLabel(AliVCluster *clu) //Returns label of the impinging particle 
{
   if (!fMCArray) 
     return 0;
      
   return clu->GetLabel();
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskGammaPHOSPP::TestGammaPt(AliCaloPhoton *ph) //returns momentum of MC particle
{

   if (!fMCArray) 
     return 0;

   Int_t iPrimaryAtVertex = ph->GetPrimaryAtVertex();
   AliAODMCParticle *particle0 = (AliAODMCParticle*) fMCArray->At(iPrimaryAtVertex);

   return (particle0->Pt());
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskGammaPHOSPP::TestTrack(AliAODTrack *track)
{
   if (!fMCArray) 
     return 0;

   Int_t TrackLabel=track->GetLabel();

   //if(TrackLabel < -2) return 0;

   AliAODMCParticle *TrackParticle =  (AliAODMCParticle*) fMCArray->At(TrackLabel);
   if (!TrackParticle) 
     return 0;

   Int_t TrackPDG = TrackParticle->GetPdgCode();

   return TrackPDG;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskGammaPHOSPP::Weight(AliAODMCParticle *particleAtVertex)
{
   //weighting for Run1 MC production
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

//_____________________________________________________________________________
Int_t AliAnalysisTaskGammaPHOSPP::TestBC(Double_t tof)
{
  Int_t bc = (Int_t)(TMath::Ceil((tof + fBCgap/2)/fBCgap) - 1);
  return bc;
}

//_____________________________________________________________________________
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

//_____________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::TestMatchingTrackPID(AliAODCaloCluster *clu1, Double_t pt)
{
    const Bool_t CPVBit = clu1->GetEmcCpvDistance() > fNsigmaCPV;
    const Bool_t DispBit = clu1->Chi2() < fNsigmaDisp;

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
    if (( nsigmaElectron < nsigmaKaon) && (nsigmaElectron < nsigmaPion) && (nsigmaElectron < nsigmaProton) && (nsigmaElectron < 1.)) 
      pidElectron1= kTRUE;
    if (!pidPion1 && !pidProton1 && !pidKaon1 && !pidElectron1) 
      pidUndef1=kTRUE;
   
    if (pidPion1)     FillHistogram("hpid1", 0.5);     
    if (pidProton1)   FillHistogram("hpid1", 1.5);     
    if (pidKaon1)     FillHistogram("hpid1", 2.5);     
    if (pidElectron1) FillHistogram("hpid1", 3.5);     
    if (pidUndef1)    FillHistogram("hpid1", 4.5);     
      
    const Int_t nmaxMatched = TMath::LocMax(AliPID::kSPECIESC, pBayesMatched);
    const Double_t distEmcCpv = clu1->GetEmcCpvDistance();
    const Int_t  trackPdg = TestTrack(trackMatched);  
 
    FillHistogram("hTracksOfClusts", pt, nmaxMatched+0.5);
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

//___________________________________________________________________________
Double_t  AliAnalysisTaskGammaPHOSPP::NonlinearMCCorrection(Double_t en)
{
   if (!fMCArray) 
      return (en);
   else

   if (fLHCRunN == 2) { 
     Double_t calib = 0.9; 
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

//____________________________________________________________________________

void AliAnalysisTaskGammaPHOSPP::GammaAcceptance()
{

  if (!fMCArray) {
    return;
  }

  auto h2 = (TH2F*)fOutputContainer2->FindObject("hGammaMC_true");
  auto h1 = (TH1D*)h2->ProjectionX("h1", 71, 170);
  
  auto h2_1 = (TH2F*)fOutputContainer2->FindObject("hCaloPhotonPdgvsPt_all");
  auto h1_1 = (TH1D*)h2_1->ProjectionX("h1_1", 4000 + 22 + 1, 4000 + 22 + 1);

  Double_t xbins[2] = {0.3, 20};

  auto htot  = (TH1D*)h1  ->Rebin(1, "htot",  xbins);
  auto hacc  = (TH1D*)h1_1->Rebin(1, "hacc", xbins); 

  hacc->Divide(hacc, htot,1, 1, "b");

  hacc->SetTitle("#gamma acceptance;;");

  hacc->SetTitle("#gamma acceptance");

  fOutputContainer2->Add(hacc);
 
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::GammaEfficiencies()
{

  if (!fMCArray) {
    return;
  }
   
  TH2F *h2 = (TH2F*)fOutputContainer2->FindObject("hGammaMC_true");
  TH1D *htot = (TH1D*)h2->ProjectionX("htotal", 71, 170);
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

//_____________________________________________________________________________
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

//___________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::AddOnePhotonHistograms()
{
  const Int_t nPt      = 400;
  const Double_t ptMin = 0;
  const Double_t ptMax = 40;

  for (auto cut : fPidCuts) {
    fOutputContainer->Add(new TH1F(Form("hCaloPhotonPt_%s", cut.first.Data()),"Photon p_{T};p_{T}, GeV/c", nPt, ptMin, ptMax));
    for (Int_t imod = 1; imod <5; imod++) {
      fOutputContainer->Add(new TH1F(Form("hCaloPhotonPt_%s_M%d", cut.first.Data(), imod), Form("Photon p_{T}, M%d;p_{T}, GeV/c", imod), nPt, ptMin, ptMax));
    }
    fOutputContainer->Add(new TH2F(Form("hCaloPhotonPtvsNcl_%s", cut.first.Data()),"Number of cells in a cluster vs photon p_{t};p_{T}, GeV/c;N_{cells}", nPt, ptMin, ptMax, 5, 3, 8));
    fOutputContainer->Add(new TH2F(Form("hCentralityvsClustPt_%s", cut.first.Data()),"Centrality vs photon p_{T};p_{T}, GeV/c;centrality", nPt, ptMin, ptMax, 8, 0., 8.));
  }

  for (auto cut : fPidCuts) {
    fOutputContainer2->Add(new TH2F(Form("hCaloPhotonPdgvsPt_%s", cut.first.Data()),"Photon pdg vs p_{T},p_{T}, GeV/c;pdg" , nPt, ptMin, ptMax, 8000, -4000, 4000));
    fOutputContainer2->Add(new TH2F(Form("hCaloPhotonPdgvsPt_%s_naive", cut.first.Data()),"Photon pdg vs p_{T} (naive);p_{T}, GeV/c;pdg (naive)" , nPt, ptMin, ptMax, 8000, -4000, 4000));
    fOutputContainer2->Add(new TH2F(Form("hMatrixEff_%s", cut.first.Data()), "Efficiency matrix;p_{T}^{meas}, GeV/c;p_{T}^{gen}, GeV/c", nPt, ptMin, ptMax, nPt, ptMin, ptMax));
    fOutputContainer2->Add(new TH2F(Form("hMatrixEff_gamma_%s", cut.first.Data()), "Efficiency matrix for #gamma;p_{T}^{meas}, GeV/c;p_{T}^{gen}, GeV/c", nPt, ptMin, ptMax, nPt, ptMin, ptMax));     
    fOutputContainer2->Add(new TH2F(Form("hMatrixEff_beta_%s", cut.first.Data()), "Efficiency matrix for #beta^{#pm};p_{T}^{meas}, GeV/c;p_{T}^{gen}, GeV/c", nPt, ptMin, ptMax, nPt, ptMin, ptMax));     
  }

}

//___________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::AddTwoPhotonHistograms() 
{
  const Int_t nM       = 750;
  const Double_t mMin  = 0.0;
  const Double_t mMax  = 1.5;
  const Int_t nPt      = 400;
  const Double_t ptMin = 0;
  const Double_t ptMax = 40;

  fOutputContainer->Add(new TH1F("hAsym",     "Asymmetry, abs((p1-p2)/(p1+p2));asym", 100, 0, 1));
  fOutputContainer->Add(new TH1F("hAsym_mix", "Asymmetry, abs((p1-p2)/(p1+p2));asym", 100, 0, 1));

  for (auto cut : fPidCuts) {
    TString histTitle = "(M,p_{T})_{#gamma#gamma};M^{#gamma#gamma}_{inv};p_{1T}+p_{2T}, GeV/c";
    fOutputContainer->Add(new TH2F(Form("hMassPt_%s",   cut.first.Data()),  histTitle  ,nM, mMin, mMax, nPt, ptMin, ptMax));
    fOutputContainer->Add(new TH2F(Form("hMiMassPt_%s", cut.first.Data()),  histTitle  ,nM, mMin, mMax, nPt, ptMin, ptMax));
    histTitle = "(M,p_{T})_{#gamma#gamma};M^{#gamma#gamma}_{inv};p_{1T}+p_{2T}, GeV/c;asym";
    fOutputContainer->Add(new TH3F(Form("hMassPt_asym_%s",  cut.first.Data()),  histTitle ,nM, mMin, mMax, nPt, ptMin, ptMax, 10, 0., 1.));
    fOutputContainer->Add(new TH3F(Form("hMiMassPt_asym_%s", cut.first.Data()), histTitle ,nM, mMin, mMax, nPt, ptMin, ptMax, 10, 0., 1.));
    histTitle = "(M,p_{T})_{#gamma#gamma};M^{#gamma#gamma}_{inv};p_{T}, GeV/c";
    fOutputContainer->Add(new TH2F(Form("hMassSingle_%s", cut.first.Data()),    histTitle , nM, mMin, mMax, nPt, ptMin, ptMax));
    fOutputContainer->Add(new TH2F(Form("hMiMassSingle_%s", cut.first.Data()),  histTitle,  nM, mMin, mMax, nPt, ptMin, ptMax));
    histTitle = "(M,p_{T})_{#gamma#gamma};M^{#gamma#gamma}_{inv};p_{T}, GeV/c;asym";
    fOutputContainer->Add(new TH3F(Form("hMassSingle_asym_%s",   cut.first.Data()), histTitle, nM, mMin, mMax, nPt, ptMin, ptMax, 10, 0., 1.));
    fOutputContainer->Add(new TH3F(Form("hMiMassSingle_asym_%s", cut.first.Data()), histTitle, nM, mMin, mMax, nPt, ptMin, ptMax, 10, 0., 1.));

    for (Int_t iMod = 1; iMod < 5; iMod ++) {  
      fOutputContainer->Add(new TH2F(Form("hMassPt_%s_M%d",       cut.first.Data(), iMod), Form("(M,p_{T})_{#gamma#gamma}, M%d",  iMod), nM, mMin, mMax, nPt, ptMin, ptMax));
      fOutputContainer->Add(new TH2F(Form("hMiMassPt_%s_M%d",     cut.first.Data(), iMod), Form("(M,p_{T})_{#gamma#gamma}, M%d",  iMod), nM, mMin, mMax, nPt, ptMin, ptMax));
      fOutputContainer->Add(new TH2F(Form("hMassSingle_%s_M%d",   cut.first.Data(), iMod), Form("(M,p_{T})_{#gamma#gamma}, M%d",  iMod), nM, mMin, mMax, nPt, ptMin, ptMax));
      fOutputContainer->Add(new TH2F(Form("hMiMassSingle_%s_M%d", cut.first.Data(), iMod), Form("(M,p_{T})_{#gamma#gamma}, M%d",  iMod), nM, mMin, mMax, nPt, ptMin, ptMax));
    }
  }
}

//___________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::AddMCHistograms() 
{
  const Int_t nPt      = 400;
  const Double_t ptMin = 0;
  const Double_t ptMax = 40;

  fOutputContainer2->Add(new TH1F("hEventCounterMC", "Count events;;events accepted", 1, 0, 1));
  
  fOutputContainer2->Add(new TH1F("hWeights","Particle weights;weight", 2000, 0., 2.));

  fOutputContainer2->Add(new TH2F("hMCPdgvsPt","Pdg code vs particle p_{T};p_{T}, GeV/c;pdg", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH1F("hMCPt","MC particles p_{T}", nPt, ptMin, ptMax));

  fOutputContainer2->Add(new TH2F("hPi0MC", "#pi^{0}-s;p_{T};y", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new TH2F("hPiPlusMC",  "#pi^{+};p_{T};y", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new TH2F("hPiMinusMC", "#pi^{-};p_{T};y", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hEtaMC", "#eta;p_{T};y", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hEtaPrimeMC", "#eta';p_{T};y", nPt, ptMin, ptMax,240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hOmegaMC", "#omega;p_{T};y", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hK0SMC", "K_{0}^{S};p_{T};y", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hK0LMC", "K_{0}^{L};p_{T};y", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hProtonMC", "p+#bar{p};p_{T};y", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hNeutrMC", "n+#bar{n};p_{T};y", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hKchMC", "K^{#pm};p_{T};y", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hPdgvsPt_MCTracks","pdg code vs p_{T};p_{T};pdg", nPt, ptMin, ptMax , 8000, -4000, 4000));

  fOutputContainer2->Add(new  TH2F("hGammaMC_all", "all #gamma;p_{T};y", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hGammaMC_true", "#gamma at IP;p_{T};y", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  
  fOutputContainer2->Add(new  TH2F("hGammaMCSources", "sources of #gamma;p_{T};y", nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new  TH2F("hCaloPhotonPdgvsPt_FromMaterial", "photon pdg vs p_{T} from material;p_{T};pdg" , nPt, ptMin, ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new  TH2F("hGammaMC_FromMaterial", "MC distribution of #gamma-s", nPt, ptMin, ptMax, 240, -1.2, 1.2));
  fOutputContainer2->Add(new  TH2F("hBetaMC_FromMaterial", "MC distribution of #gamma-s", nPt, ptMin, ptMax, 240, -1.2, 1.2));

  fOutputContainer2->Add(new TH1I("hpid3", "Pid #sigma<3, #pi, p, K, #beta, Undef", 5, 0., 5.));
  fOutputContainer2->Add(new TH1I("hpid1", "Pid #sigma<1, #pi, p, K, #beta, Undef", 5, 0., 5.));
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::AddTrackHistograms() 
{
  
  const Int_t nPt      = 400;
  const Double_t ptMin = 0;
  const Double_t ptMax = 40;
  
  fOutputContainer->Add(new TH1F("hDistanceToBadChannel", "Distance to bad channel in cm", 1000, 0, 100));

  fOutputContainer2->Add(new TH1F("hMCTrackCharge","Charge of MC track", 9, -4.5, 4.5));
 
  fOutputContainer->Add(new TH2F("hTracksOfClusts","Tracks vs p_{T}.", nPt, ptMin, ptMax, 9, 0.,9.));
  fOutputContainer->Add(new TH2F("hTracksOfClusts_cpv","Tracks vs p_{T}, CPV cut.", nPt, ptMin, ptMax, 9, 0.,9.));
  fOutputContainer->Add(new TH2F("hTracksOfClusts_disp","Tracks vs p_{T}, shape cut.", nPt, ptMin, ptMax, 9, 0.,9.));
  fOutputContainer->Add(new TH2F("hTracksOfClusts_both","Tracks vs p_{T}, CPV + shape cuts." , nPt, ptMin, ptMax, 9, 0.,9.));

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

  // MC 
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

//________________________________________________________________________________________________
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

  fOutputContainer->Add(new TH2F("hEnergyResolution_EM", "Energy resolution;E_{MC};E_{meas}-E_{MC}, EM probes", 400, 0, 40, 320, -1.6, 1.6));
  fOutputContainer->Add(new TH2F("hEnergyResolution_other", "Energy resolution;E_{MC};E_{meas}-E_{MC}, other probes", 400, 0, 40, 320, -1.6, 1.6));

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

  fOutputContainer->Add(new TH1F("hPhotonKappa","#kappa(#gamma)", 400, 0., 40.));
  fOutputContainer->Add(new TH1F("hPhotonPt","p_{T}(#gamma)", 400, 0., 40.));
  fOutputContainer->Add(new TH1F("hPhotonPx","p_{x}(#gamma)", 400, 0., 40.));
  fOutputContainer->Add(new TH1F("hPhotonPy","p_{y}(#gamma)", 400, 0., 40.));

  fOutputContainer->Add(new TH1F("hTrigClass","Trigger class",5, 0.5,5.5));
  fOutputContainer->Add(new TH1F("hNPileupVtx","Number of SPD pileup vertices", 10, 0., 10.));
  fOutputContainer->Add(new TH1F("hZPileupVtx", "Z pileup", 200, -50., 50.));

  fOutputContainer->Add(new TH1F("hZvertex","Z vertex",400, -20.,+20.));
  fOutputContainer->Add(new TH1F("hSPDvertex","SPD vertex",400, -20.,+20.));
  fOutputContainer->Add(new TH1F("hDistSPDTrackvertex","z-distance between track and SPD vertex", 200, 0., 20.));
  fOutputContainer->Add(new TH2F("hTrackvsSPDVertices", "Track vs SPD vertices;SPDV, cm;TrackV, cm", 400, -20, 20, 400, -20, 20));
  fOutputContainer->Add(new TH1F("hNContributors","N of primary tracks from the primary vertex", 150, 0., 150.));
  fOutputContainer->Add(new TH1F("hTrackMult","Charged track multiplicity", 150, 0., 150.));

  fOutputContainer->Add(new TH1F("hV0Atime","V0A time", 1200, -6.e-6,+6.e-6));
  fOutputContainer->Add(new TH1F("hV0Ctime","V0C time", 1200, -6.e-6,+6.e-6));
  fOutputContainer->Add(new TH2F("hV0AV0Ctime","V0A time vs V0C time", 120, -6.e-6,+6.e-6 , 120, -6.e-6,+6.e-6));

  fOutputContainer->Add(new TH1D("hTOF","Time of Flight", 1000, 0, 1e-6));
  fOutputContainer->Add(new TH1I("hBC", "BC", 1000, 0, 1000));

  fOutputContainer->Add(new TH1F("hvt0vsvt5", "Difference of momenta for different vertices", 200, -1., 1.));

  fOutputContainer->Add(new TH2F("hClustM02","Cluster M02 vs p_{T}", nPt, ptMin, ptMax, 100, 0, 10));
  fOutputContainer->Add(new TH2F("hClustM20","Cluster M20 vs p_{T}", nPt, ptMin, ptMax, 100, 0, 10));
}
