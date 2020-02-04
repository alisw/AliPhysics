/**************************************************************************
 * Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//*************************************************************************
// Class AliAnalysisTaskSEDvsRT
// AliAnalysisTaskSE for charmed hadrons vs. transverse activity analysis
// Authors: Jeremy Wilkinson,
/////////////////////////////////////////////////////////////


#include <TClonesArray.h>
#include <TCanvas.h>
#include <TList.h>
#include <TString.h>
#include <TDatabasePDG.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TProfile.h>
#include "AliAnalysisManager.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsLctoV0.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEDvsRT.h"
#include "AliNormalizationCounter.h"
#include "AliVertexingHFUtils.h"
#include "AliAODVZERO.h"
#include "AliESDUtils.h"


//________________________________________________________________________
AliAnalysisTaskSEDvsRT::AliAnalysisTaskSEDvsRT():
   AliAnalysisTaskSE(),
   fOutput(0),
   fListCuts(0),
   fHistNEvents(0),
   fRDCutsAnalysis(cuts)
   {
      ///default constructor
   }

AliAnalysisTaskSEDvsRT::AliAnalysisTaskSEDvsRT(const char *name, Int_t pdgMeson, AliRDHFCuts *cuts):
   AliAnalysisTaskSE(name),
   fOutput(0),
   fListCuts(0),
   fHistNEvents(0),
   fRDCutsAnalysis(cuts)
   {
      ///Default constructor
      DefineOutput(1,TList::Class()); //Output slot 1: fOutput TList container
      DefineOutput(2,TList::Class()); //Output slot 2: fListCuts
      DefineOutput(3,TList::Class()); //Output slot 3: normalisation counters
   
   
//________________________________________________________________________
AliAnalysisTaskSEDvsRT::~AliAnalysisTaskSEDvsRT()
{
   //
   /// Standard destructor
   //
   delete fOutput;
}

//________________________________________________________________________
void AliAnalysisTaskSEDvsRT::Init(){
   //
   /// Initialisation
   //
   
   printf("AliAnalysisTaskSEDvsRT_0::Init() \n");
   
   fListCuts = new TList();
   fListCuts->SetOwner();
   fListCuts->SetName("CutsList");
   
   if(fPdgSpecies==411){
    AliRDHFCutsDplustoKpipi* copycut=new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi*>(fRDCutsAnalysis)));
    copycut->SetName("AnalysisCutsDplus");
    fListCuts->Add(copycut);
  }else if(fPdgSpecies==421){
    AliRDHFCutsD0toKpi* copycut=new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fRDCutsAnalysis)));
    copycut->SetName("AnalysisCutsDzero");
    fListCuts->Add(copycut);
  }else if(fPdgSpecies==413){
    AliRDHFCutsDStartoKpipi* copycut=new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fRDCutsAnalysis)));
    copycut->SetName("AnalysisCutsDStar");
    fListCuts->Add(copycut);
  }else if(fPdgSpecies==431){
    AliRDHFCutsDstoKKpi* copycut=new AliRDHFCutsDstoKKpi(*(static_cast<AliRDHFCutsDstoKKpi*>(fRDCutsAnalysis)));
    copycut->SetName("AnalysisCutsDs");
    fListCuts->Add(copycut);
  }else if(fPdgSpecies==4122){
    if(fLctoV0){
      AliRDHFCutsLctoV0* copycut=new AliRDHFCutsLctoV0(*(static_cast<AliRDHFCutsLctoV0*>(fRDCutsAnalysis)));
      copycut->SetName("AnalysisCutsLc2pK0S");
      fListCuts->Add(copycut);
      }else{
      AliRDHFCutsLctopKpi *copycut=new AliRDHFCutsLctopKpi(*(static_cast<AliRDHFCutsLctopKpi*>(fRDCutsAnalysis)));
      copycut->SetName("LctopKpiProdCuts");
      fListCuts->Add(copycut);
      }
  }
   PostData(2,fListCuts);
}

//________________________________________________________________________
void AliAnalysisTaskSEDvsRT::UserCreateOutputObjects()
{
   /// Create output containers
   
   if (fDebug > 1) printf("AliAnalysisTaskSEDvsRT::UserCreateOutputObjects() \n");
   
   // TList for output
   fOutput = new TList();
   fOutput->SetOwner();
   fOutput->SetName("OutputHistos");
   
   fOutputCounters = new TList();
   fOutputCounters->SetOwner();
   fOutputCounters->SetName("OutputCounters");
   fOutputCounters->Add(fCounter);
   
   
   fPtvsMassvsRTToward = new TH3F("hPtvsMassvsRTToward", "D candidates in toward region: p_{T} vs mass vs R_{T};R_{T};Mass [GeV/c^{2}];p_{T} [GeV/c]",nRTBins,firstRTbin,lastRTbin, fNMassBins,fLowmasslimit,fUpmasslimit,42,0.,24.);
fPtvsMassvsRTAway = new TH3F("hPtvsMassvsRTAway", "D candidates in away region: p_{T} vs mass vs R_{T};R_{T};Mass [GeV/c^{2}];p_{T} [GeV/c]",nRTBins,firstRTbin,lastRTbin, fNMassBins,fLowmasslimit,fUpmasslimit,42,0.,24.);
   
   
   
   
   PostData(1,fOutput);
   PostData(3,fOutputCounters);
   
   return;   
}

//________________________________________________________________________
void AliAnalysisTaskSEDvsRT::UserExec(Option_t */*option*/)
{
   /// Execute analysis for current event:
   /// heavy flavour candidates as function of RT
   
   AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
   
   //  AliAODTracklets* tracklets = aod->GetTracklets();
   //Int_t ntracklets = tracklets->GetNumberOfTracklets();
   if(fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      return;
    }
  }
   
    
  TClonesArray *arrayCand = 0;
  TString arrayName="";
  UInt_t pdgDau[3];
  Int_t nDau=0;
  Int_t selbit=0;
  if(fPdgSpecies==411){
    arrayName="Charm3Prong";
    pdgDau[0]=211; pdgDau[1]=321; pdgDau[2]=211; 
    nDau=3;
    selbit=AliRDHFCuts::kDplusCuts;
  }else if(fPdgSpecies==421){
    arrayName="D0toKpi";
    pdgDau[0]=211; pdgDau[1]=321; pdgDau[2]=0;
    nDau=2;
    selbit=AliRDHFCuts::kD0toKpiCuts;
  }else if(fPdgSpecies==413){
    arrayName="Dstar";
    pdgDau[0]=321; pdgDau[1]=211; pdgDau[2]=0; // Quoting here D0 daughters (D* ones on another variable later)
    nDau=2;
    selbit=AliRDHFCuts::kDstarCuts;
  }else if(fPdgSpecies==431){
    arrayName="Charm3Prong";
    pdgDau[0]=321; pdgDau[1]=321; pdgDau[2]=211;
    nDau=3;
    selbit=AliRDHFCuts::kDsCuts;
  }else if(fPdgSpecies==4122){
    if(fLctoV0){
    arrayName="CascadesHF";
    pdgDau[0]=211; pdgDau[1]=211; pdgDau[2]=0; // Quoting here K0S daughters (Lc ones on another variable later)
    nDau=2;
    selbit=AliRDHFCuts::kLctoV0Cuts;
    }else{
    arrayName="Charm3Prong";
    pdgDau[0]=2212; pdgDau[1]=321; pdgDau[2]=211;
    nDau=3;
    selbit=AliRDHFCuts::kLcCuts;
    }
  }

//Arrays of daughters for inv mass calculation
    UInt_t pdgdaughtersD0[2] = {211,321}; //pi, K
    UInt_t pdgdaughtersD0bar[2] = {321,211}; //K, pi
    UInt_t pdgDsKKpi[3] = {321,321,211};
    UInt_t pdgDspiKK[3] = {311, 321,321};  
  
  
  
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      arrayCand=(TClonesArray*)aodFromExt->GetList()->FindObject(arrayName.Data());
    }
  } else if(aod) {
    arrayCand=(TClonesArray*)aod->GetList()->FindObject(arrayName.Data());
  }

  if(!aod || !arrayCand) {
    printf("AliAnalysisTaskSEDvsMultiplicity::UserExec: Charm3Prong branch not found!\n");
    return;
  }

  if(fisPPbData && fReadMC){
    Int_t runnumber = aod->GetRunNumber();
    if(aod->GetTriggerMask()==0 && 
       (runnumber>=195344 && runnumber<=195677)){
      AliDebug(3,"Event rejected because of null trigger mask");
      return;
    }
  }


  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex()||TMath::Abs(aod->GetMagneticField())<0.001) return;
  
  
  ///!----- l.797-833 for multiplicity estimation; RT determination goes here
  
  Double_t rtval= 0;

  fCounter->StoreEvent(aod,fRDCutsAnalysis,fReadMC,RT);
  fHistNEvents->Fill(0);
  
  //!----l.839-888: multiplicity correction
  
  Bool_t isEvSel = fRDCutsAnalysis->IsEventSelected(aod); //?
  
   switch (fRDCutsAnalysis->GetWhyRejection()) {
      case 5: fHistNEvents->Fill(3); break;
      case 7: fHistNEvents->Fill(4); break;
      case 6: fHistNEvents->Fill(5); break;
      case 1: fHistNEvents->Fill(6); break;
      
   }
  
  
  Bool_t isEvPSRejected = fRDCutsAnalysis->IsEventRejectedDuePhysicsSelection();
  Bool_t isEvTrigNameRejected = fRDCutsAnalysis->IsEventRejectedDueToTrigger();
  Bool_t isEvPileUpRejected = fRDCutsAnalysis->IsEventRejectedDueToPileup();
  Bool_t isEvNoVtxRejected = fRDCutsAnalysis->IsEventRejectedDueToNotRecoVertex();
  Bool_t isEvVtxContribRejected = fRDCutsAnalysis->IsEventRejectedDueToVertexContributors();
  Bool_t isEvVtxRangeRejected= fRDCutsAnalysis->IsEventRejectedDueToZVertexOutsideFiducialRegion();
  Bool_t isEvCentralityRejected = fRDCutsAnalysis->IsEventRejectedDueToCentrality();
  if(!isEvPSRejected){
    fHistNtrUnCorrPSSel->Fill(countMult);
    fHistNtrCorrPSSel->Fill(countCorr);
    if(!isEvTrigNameRejected){
      fHistNtrUnCorrPSTrigSel->Fill(countMult);
      if(!isEvPileUpRejected){
	fHistNtrUnCorrPSTrigPileUpSel->Fill(countMult);
	if(!isEvNoVtxRejected){
	  fHistNtrUnCorrPSTrigPileUpVtxSel->Fill(countMult);
	  if(!isEvVtxContribRejected){
	    fHistNtrUnCorrPSTrigPileUpVtxContSel->Fill(countMult);
	    if(!isEvVtxRangeRejected){
	      fHistNtrUnCorrPSTrigPileUpVtxRangeSel->Fill(countMult);
	      if(!isEvCentralityRejected){
		fHistNtrUnCorrPSTrigPileUpVtxRangeCentrSel->Fill(countMult);
	      }
	    }
	  }
	}
      }
    }
  }
  
  if(!isEvSel)return;
  fHistNEvents->Fill(2);

//l.948-1094: MC multiplicity counting/reweighting  
  
  Int_t nCand = arrayCand->GetEntriesFast();
  Int_t nSelectedNoPID=0,nSelectedPID=0,nSelectedInMassPeak=0;
  Double_t mD0PDG    = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t mDplusPDG = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  Double_t mDstarPDG = TDatabasePDG::Instance()->GetParticle(413)->Mass();
  Double_t mDsPDG    = TDatabasePDG::Instance()->GetParticle(431)->Mass();
  Double_t mLcPDG    = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  
  // PDG of daughters for D*
  UInt_t pdgDgDStartoD0pi[2] = {421, 211};
  
  // PDG of daughters for Lc2pK0
  UInt_t pdgDgLctopK0S[2] = {2212, 310};
  
   // omitting "aveMult" l.1110
  
  //TODO:  Loop on candidates, perform selection, determine phi of candidate wrt leading particle, fill corresponding histo (toward/away)
  
  for (Int_t iCand = 0; iCand < nCand; iCand++) {  //Loop over candidates
     
     AliAODRecoDecayHF *d =(AliAODRecoDecayHF*)arrayCand->UncheckedAt(iCand);
     AliAODRecoCascadeHF *dCascade = NULL;
     if (fPdgSpecies == 413 || (fPdgSpecies==4122 && fLctoV0)) dCascade = (AliAODRecoCascadeHF*)d;
     
     fHistNEvents->Fill(7);
     
     if(fPdgSpecies == 4122) { //V0 selection for LctopK0
        if(fLctoV0) {
           AliAODv0 *v0part = (AliAODv0*)dCascade->Getv0();
           Bool_t onFlyV0 = v0part->GetOnFlyStatus();
           if (onFlyV0)  {fHistNEvents->Fill(8); continue;}
        } else if (!fLctoV0 && fUseBit && !d->HasSelectionBit(selbit)) { // check selection bit for LcpKpi
           fHistNEvents->Fill(8); continue;
        }
     }
    
     if (fPdgSpecies != 4122 && fUseBit && !d->HasSelectionBit(selbit)) { //check selection bit for all others
         fHistNEvents->Fill(8); continue;
     }

     Double_t ptCant = d->Pt();
     Double_t rapid = d->Y(fPdgSpecies);
     Bool_t isInFidAcc = fDRCutsAnalysis->IsInFiducialAcceptance(ptCand,rapid);
     if (!isInFidAcc) continue;
     
     Int_t labD = -1;
     if(fReadMC) {
      if(fPdgMeson==413){
         labD = dCascade->MatchToMC(fPdgMeson,421,(Int_t*)pdgDgDStartoD0pi,(Int_t*)pdgDau,arrayMC);
      } else if(fPdgMeson==4122 && fLctoV0){
         labD = dCascade->MatchToMC(fPdgMeson,pdgDgLctopK0S[1],(Int_t*)pdgDgLctopK0S,(Int_t*)pdgDau,arrayMC,kTRUE);
      } else {
         labD = d->MatchToMC(fPdgMeson,arrayMC,nDau,(Int_t*)pdgDau);
      }
      FillMCMassHistos(arrayMC,labD, countMult,nchWeight);
     }
     
     Int_t passAllCuts = 0, passTopolCuts = 0;
     
     if (fPdgSpecies == 4122) {  //special cuts cases for Lc
        if (fLctoV0) { //LctopK0
           passAllCuts = (((fRDCutsAnalysis->IsSelected(d, AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLctoK0Spr)) == (AliRDHFCutsLctoV0::kLctoK0Spr));
           passTopolCuts = (((fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr));
        }
        else { //LctopKpi
           passAllCuts = fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kAll,aod);
           passTopolCuts= fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kCandidate,aod);
        }
     }
     else { //everything that isn't Lc
         passAllCuts = fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kAll,aod);
         passTopolCuts = fRDCutsAnalysis->GetIsSelectedCuts();
     }
     
     if (fPdgSpecies != 431 && passTopolCuts == 0) continue;
     nSelectedNoPID++;
     fHistNEvents->Fill(9);
     if (fPdgSpecies == 431 && passAllCuts == 0) continue;
     if (passAllCuts) {
        nSelectedPID++;
        fHistNEvents->Fill(10);
     }
     //daughter subtraction from multiplicity, l.1169-1192
     
     Bool_t isPrimary = kTRUE;
     Double_t trueImpParXY = 9999.;
     Double_t impparXY = d->ImpParXY()*10000.;
     Double_t dlen = 0.1;
     Double_t RTforCand = -1;
     Double_t candDeltaPhi = -1; //delta-phi wrt leading particle
     Double_t mass[2];  //Calculate inv mass and check if in peak region
     if (fPdgSpecies == 411) {   //D+
         mass[0] = d->InvMass(nDau,pdgDau);
         mass[1] = -1.;
         if (TMath::Abs(mass[0]-mDplusPDG) < 0.02) nSelectedInMassPeak++;
     } else if (fPdgSpecies == 421) { //D0

         mass[0] = d->InvMass(2,pdgdaughtersD0);
         mass[1] = d->InvMass(2,pdgdaughtersD0bar);
         if (TMath::Abs(mass[0]-mD0PDG)<0.02 || TMath::Abs((mass[1]-mD0PDG) < 0.02) nSelectedInMassPeak++;
     } else if (fPdgSpecies == 413) { //D*
         mass[0] = dCascade->DeltaInvMass();
         mass[1] = -1.;
         if (TMath::Abs(mass[0] - (mDstarPDG-mD0PDG)) < 0.0015 ) nSelectedInMassPeak++;
     } else if (fPdgSpecies == 431) { //Ds
         mass[0] = d->InvMass(nDau,pdgDsKKpi);
         mass[1] = d->InvMass(nDau,pdgDspiKK);
         if (TMath::Abs(mass[0]-mDsPDG) < 0.02 || TMath::Abs(mass[1]-mDsPDG) < 0.02) nSelectedInMassPeak++;
     } else if (fPdgSpecies == 4122) { //Lc
         if (fLctoV0) { //pK0 inv mass
            mass[0] = d->InvMass(2,pdgDgLctopK0S);
            mass[1] = -1.;
            if (TMath::Abs(mass[0]-mLcPDG)< 0.02) nSelectedInMassPeak++;
         } else { //pKpi inv mass
            UInt_t pdgpKpi[3] = {2212,321,211};
            UInt_t pdgpiKp[3] = {211, 321,2212};
            if (passTopolCuts == 3 || passTopolCuts == 1) mass[0] = d->InvMass(3,pdgpKpi);
            if (passTopolCuts >= 2 ) mass[1] = d->InvMass(3,pdgpiKp);
            if (TMath::Abs(mass[0] - mLcPDG) < 0.02 ) nSelectedInMassPeak++;
         }
     }
            
     for (Int_t iHyp = 0; iHyp < 2; iHyp++) {
        if (mass[iHyp] < 0.) continue; //for D+ , D*, Lc2pK0, only one mass hypothesis
        Double_t invMass = mass[iHyp];
        Double_t arrayForSparse[5] = {invMass, ptCand, impparXY, dlen, RTForCand, candDeltaPhi};  //! RT value to be refined here
     
     
     
     }
     
     
  
  
  
  
  
  PostData(1, fOutput);
  PostData(2, fListCuts);
  PostData(3, fOutputCounters);
   
}

void AliAnalysisTaskSEDvsRT::Terminate(Option_t */*option*/)
{
   /// Terminate analysis
   ///
   
   if (fDebug > 1) printf("AliAnalysisTaskSEDvsRT: Terminate() \n");
   fOutput = dynamic_cast<TList*> (GetOutputData(1));
   if (!fOutput) {
      printf("ERROR: fOutput not found\n");
      return;
   }
   
   fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));
   if (!fHistNEvents){
      printf("ERROR: fHistNEvents not found\n");
      return;
   }
   printf("Number of events analysed = %d\n",(Int_t)fHistNEvents->GetBinContent(3));
   return; 
}
