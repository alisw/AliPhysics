/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

//_________________________________________________________________________
// Gerenate a random trigger, input for other analysis
// Set flat energy distribution over acceptance of EMCAL, PHOS or CTS
// Be careful, correlate only with Min Bias events this trigger
//
//
//-- Author: Gustavo Conesa (LPSC-Grenoble)
//_________________________________________________________________________


// --- ROOT system ---
#include <TH2F.h>
#include <TClonesArray.h>

//---- AliRoot system ----
#include "AliAnaRandomTrigger.h"
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliEMCALGeometry.h"

ClassImp(AliAnaRandomTrigger)
  
//__________________________________________
AliAnaRandomTrigger::AliAnaRandomTrigger() : 
    AliAnaCaloTrackCorrBaseClass(),
    fDetector("EMCAL"), fRandom(0), fNRandom(0),
    fhE(0),             fhPt(0),
    fhPhi(0),           fhEta(0), 
    fhEtaPhi(0) 
{
  //Default Ctor

  //Initialize parameters
  InitParameters();

}

//_____________________________________________________________________________________
Bool_t AliAnaRandomTrigger::ExcludeDeadBadRegions(const Float_t eta, const Float_t phi)
{
  // Check if there is a dead or bad region in a detector
  // Now only EMCAL
  
  if(fDetector!="EMCAL") return kFALSE;
  
  //-------------------------------------
  // Get the corresponding cell in EMCAL, check if it exists in acceptance (phi gaps, borders)
  //-------------------------------------

  Int_t absId = -1;
  if(!GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(eta,phi, absId)) return kTRUE; // remove if out of EMCAL acceptance, phi gaps
  
  Int_t icol = -1, irow = -1, iRCU = -1;
  Int_t sm = GetCaloUtils()->GetModuleNumberCellIndexes(absId,"EMCAL", icol, irow, iRCU);
  
  //printf("eta %f, phi %f, ieta %d, iphi %d, sm %d\n",eta,phi,icol,irow,sm);
  
  //-------------------------------------
  // Remove in case of close to border, by default always 1 but check what was set in reco utils
  //-------------------------------------

  Bool_t okrow = kFALSE;
	Bool_t okcol = kFALSE;
  Int_t nborder = GetCaloUtils()->GetEMCALRecoUtils()->GetNumberOfCellsFromEMCALBorder();
  if (nborder<1) nborder = 1;
  
  // Rows
  if(sm < 10)
  {
    if(irow >= nborder && irow < 24-nborder) okrow =kTRUE; 
  }
  else
  {
    if((GetCaloUtils()->EMCALGeometryName()).Contains("12SM")) // 1/3 SM
    {
      if(irow >= nborder && irow < 8-nborder) okrow =kTRUE; 
    }
    else // 1/2 SM
    {
      if(irow >= nborder && irow <12-nborder) okrow =kTRUE; 
    }
  }
  
  // Columns
  if(sm%2==0)
  {
    if(icol >= nborder)     okcol = kTRUE;	
  }
  else 
  {
    if(icol <  48-nborder)  okcol = kTRUE;	
  }
  
  //printf("okcol %d, okrow %d\n",okcol,okrow);
  if (!okcol || !okrow) return kTRUE; 
  
  //-------------------------------------
  // Check if the cell or those around are bad
  //-------------------------------------

  if(GetCaloUtils()->GetEMCALChannelStatus(sm,icol, irow) > 0) return kTRUE ; // trigger falls into a bad channel

  // Check if close there was a bad channel
//  for(Int_t i = -1; i <= 1; i++)
//  {
//    for(Int_t j = -1; j <= 1; j++)
//    {
//      //printf("\t check icol %d, irow %d \n",icol+i, irow+j);
//      if(GetCaloUtils()->GetEMCALChannelStatus(sm,icol+i, irow+j) > 0) return kTRUE ; // trigger falls into a bad channel
//      //printf("\t ok\n");
//    }
//  }

   //printf("\t OK\n");
  
  return kFALSE;
  
}


//__________________________________________________
TObjString *  AliAnaRandomTrigger::GetAnalysisCuts()
{  	
  //Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaRandomTrigger ---\n") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Detector: %s\n"    , fDetector.Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"N per event = %d\n", fNRandom       ) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Min E   = %3.2f - Max E   = %3.2f\n", GetMinPt(), GetMaxPt()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Min Eta = %3.2f - Max Eta = %3.2f\n", fEtaCut[0], fEtaCut[1]) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Min Phi = %3.2f - Max Phi = %3.2f\n", fPhiCut[0], fPhiCut[1]) ;
  parList+=onePar ;
   
  return new TObjString(parList) ;
}

//_______________________________________________________
TList *  AliAnaRandomTrigger::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in fOutputContainer
  
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("RandomTrigger") ; 
  
  Int_t nptbins  = GetHistogramRanges()->GetHistoPtBins(); Int_t nphibins = GetHistogramRanges()->GetHistoPhiBins(); Int_t netabins = GetHistogramRanges()->GetHistoEtaBins();
  Float_t ptmax  = GetHistogramRanges()->GetHistoPtMax();  Float_t phimax = GetHistogramRanges()->GetHistoPhiMax();  Float_t etamax = GetHistogramRanges()->GetHistoEtaMax();
  Float_t ptmin  = GetHistogramRanges()->GetHistoPtMin();  Float_t phimin = GetHistogramRanges()->GetHistoPhiMin();  Float_t etamin = GetHistogramRanges()->GetHistoEtaMin();	

  fhE  = new TH1F ("hE","Random E distribution", nptbins,ptmin,ptmax); 
  fhE->SetXTitle("E (GeV)");
  outputContainer->Add(fhE);
  
  fhPt  = new TH1F ("hPt","Random p_{T} distribution", nptbins,ptmin,ptmax); 
  fhPt->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPt);
  
  fhPhi  = new TH2F ("hPhi","Random #phi distribution",
                        nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
  fhPhi->SetYTitle("#phi (rad)");
  fhPhi->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPhi);
  
  fhEta  = new TH2F ("hEta","Random #eta distribution",
                        nptbins,ptmin,ptmax, netabins,etamin,etamax); 
  fhEta->SetYTitle("#eta ");
  fhEta->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhEta);
  
  fhEtaPhi  = new TH2F ("hEtaPhi","Random #eta vs #phi ",netabins,etamin,etamax, nphibins,phimin,phimax); 
  fhEtaPhi->SetXTitle("#eta ");
  fhEtaPhi->SetYTitle("#phi (rad)");  
  outputContainer->Add(fhEtaPhi);
    
  return outputContainer;

}

//___________________________________________
void AliAnaRandomTrigger::InitParameters()
{ 
  //Initialize the parameters of the analysis.
  SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  SetOutputAODName("RandomTrigger");

  AddToHistogramsName("AnaRandomTrigger_");
  
  fNRandom   = 1    ;
  fPhiCut[0] = 0.   ;
  fPhiCut[1] = TMath::TwoPi() ;
  fEtaCut[0] =-1.   ;
  fEtaCut[1] = 1.   ;
  
}

//____________________________________________________________
void AliAnaRandomTrigger::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");	

  printf("Detector = %s\n",  fDetector.Data());
  printf("Min E   = %3.2f - Max E   = %3.2f\n", GetMinPt(), GetMaxPt());
  printf("Min Eta = %3.2f - Max Eta = %3.2f\n", fEtaCut[0], fEtaCut[1]);
  printf("Min Phi = %3.2f - Max Phi = %3.2f\n", fPhiCut[0], fPhiCut[1]);

} 

//______________________________________________
void  AliAnaRandomTrigger::MakeAnalysisFillAOD() 
{
  // Do analysis and fill aods
  // Generate particle randomly
  // fNRandom particles per event
  
  for(Int_t irandom = 0; irandom < fNRandom; irandom++)
  {
    // Get the random variables of the trigger
    Float_t pt  = fRandom.Uniform(GetMinPt(), GetMaxPt());
    Float_t eta = fRandom.Uniform(fEtaCut[0], fEtaCut[1]);
    Float_t phi = fRandom.Uniform(fPhiCut[0], fPhiCut[1]);
    
    // Check if particle falls into a dead region, if inside, get new
    Bool_t excluded =  ExcludeDeadBadRegions(eta,phi);
    
    // if excluded, generate a new trigger until accepted
    while (excluded)
    {
      pt  = fRandom.Uniform(GetMinPt(), GetMaxPt());
      eta = fRandom.Uniform(fEtaCut[0], fEtaCut[1]);
      phi = fRandom.Uniform(fPhiCut[0], fPhiCut[1]);
      
      excluded = ExcludeDeadBadRegions(eta,phi);
    }
    
    // Create the AOD trigger object
    TLorentzVector mom;
    mom.SetPtEtaPhiM(pt,eta,phi,0);
    
    AliAODPWG4Particle trigger = AliAODPWG4Particle(mom);
    trigger.SetDetector(fDetector);
    
    if(GetDebug() > 1) 
      printf("AliAnaRandomTrigger::MakeAnalysisFillAOD() -  iRandom %d, Trigger e %2.2f pt %2.2f, phi %2.2f, eta %2.2f \n",
             irandom, trigger.E(), trigger.Pt(), trigger.Phi(), trigger.Eta());
    
    AddAODParticle(trigger);
  }
  
  if(GetDebug() > 0) 	
    printf("AliAnaRandomTrigger::MakeAnalysisFillAOD() - Final aod branch entries %d\n", GetOutputAODBranch()->GetEntriesFast());   
} 

//_____________________________________________________
void  AliAnaRandomTrigger::MakeAnalysisFillHistograms() 
{
  // Fill histograms
  
  //Loop on stored AODParticles
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  
  if(GetDebug() > 0) 
    printf("AliAnaRandomTrigger::MakeAnalysisFillHistograms() - aod branch entries %d, fNRandom %d\n", naod, fNRandom);
  
  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    AliAODPWG4Particle* trigger =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    
    fhPt    ->Fill(trigger->Pt());
    fhE     ->Fill(trigger->E());
    fhPhi   ->Fill(trigger->Pt(), trigger->Phi());
    fhEta   ->Fill(trigger->Pt(), trigger->Eta());
    fhEtaPhi->Fill(trigger->Eta(),trigger->Phi());
    
  }// aod branch loop
  
}
