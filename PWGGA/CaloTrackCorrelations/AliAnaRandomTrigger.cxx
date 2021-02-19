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

// --- ROOT system ---
#include <TObjString.h>
#include <TH2F.h>
#include <TClonesArray.h>

//---- AliRoot system ----
#include "AliAnaRandomTrigger.h"
#include "AliCaloTrackParticleCorrelation.h"
#include "AliEMCALGeometry.h"

/// \cond CLASSIMP
ClassImp(AliAnaRandomTrigger) ;
/// \endcond

//__________________________________________
/// Default Constructor. Initialize parameters.
//__________________________________________
AliAnaRandomTrigger::AliAnaRandomTrigger() :
    AliAnaCaloTrackCorrBaseClass(),
    fTriggerDetector(kEMCAL),
    fTriggerDetectorString("EMCAL"),
    fRandom(0),         fNRandom(0),
    fMomentum(),
    fhPt(0),
    fhPhi(0),           fhEta(0), 
    fhEtaPhi(0) 
{
  InitParameters();
}

//_________________________________________________________________________
/// Check if there is a dead or bad region in a detector.
/// Only EMCAL for now.
//_________________________________________________________________________
Bool_t AliAnaRandomTrigger::ExcludeDeadBadRegions(Float_t eta, Float_t phi)
{
  if ( fTriggerDetector != kEMCAL && fTriggerDetector != kDCAL ) return kFALSE;
  
  //-------------------------------------
  // Get the corresponding cell in EMCAL, check if it exists in acceptance (phi gaps, borders)
  //-------------------------------------

  Int_t absId = -1;
  if ( !GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(eta,phi, absId) ) 
    return kTRUE; // remove if out of EMCAL acceptance, phi gaps
  
  Int_t icol = -1, irow = -1, iRCU = -1;
  Int_t sm = GetCaloUtils()->GetModuleNumberCellIndexes(absId,kEMCAL, icol, irow, iRCU);
  
  //printf("Exclude? eta %f, phi %f, ieta %d, iphi %d, sm %d\n",eta,phi,icol,irow,sm);
  
  //-------------------------------------
  // Remove in case of close to border, by default always 1 but check what was set in reco utils
  //-------------------------------------

  Bool_t okrow = kFALSE;
	Bool_t okcol = kFALSE;
  Int_t nborder = GetCaloUtils()->GetEMCALRecoUtils()->GetNumberOfCellsFromEMCALBorder();
  if ( nborder < 1 ) nborder = 1;
  
  // Rows
  if ( (sm < 10 && sm >=  0) || // EMCal full
       (sm < 18 && sm >= 12)   )// DCal full
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
  if ( sm%2==0 )
  {
    if ( icol >= nborder )     okcol = kTRUE;	
  }
  else 
  {
    if ( icol <  48-nborder )  okcol = kTRUE;	
  }
  
  //printf("Exclude: okcol %d, okrow %d\n",okcol,okrow);
  if (!okcol || !okrow) return kTRUE; 
  
  //-------------------------------------
  // Check if the cell or those around are bad
  //-------------------------------------

  Int_t status = 0;
  if ( GetCaloUtils()->GetEMCALChannelStatus(sm,icol, irow,status) ) 
  {
    //printf("Exclude: bad channel\n");
    return kTRUE ; // trigger falls into a bad channel
  }
  
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
/// Save parameters used for analysis.
//__________________________________________________
TObjString *  AliAnaRandomTrigger::GetAnalysisCuts()
{  	
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaRandomTrigger ---:") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Detector: %s;"    , fTriggerDetectorString.Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"N per event = %d;", fNRandom       ) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Min E   = %3.2f - Max E   = %3.2f;", GetMinPt(), GetMaxPt()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Min Eta = %3.2f - Max Eta = %3.2f;", fEtaCut[0], fEtaCut[1]) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Min Phi = %3.2f - Max Phi = %3.2f;", fPhiCut[0], fPhiCut[1]) ;
  parList+=onePar ;
   
  return new TObjString(parList) ;
}

//____________________________________________________
/// Create histograms to be saved in output file and
/// store them in fOutputContainer
//____________________________________________________
TList *  AliAnaRandomTrigger::GetCreateOutputObjects()
{
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("RandomTrigger") ; 
  
  Int_t nptbins  = GetHistogramRanges()->GetHistoPtBins(); Int_t nphibins = GetHistogramRanges()->GetHistoPhiBins(); Int_t netabins = GetHistogramRanges()->GetHistoEtaBins();
  Float_t ptmax  = GetHistogramRanges()->GetHistoPtMax();  Float_t phimax = GetHistogramRanges()->GetHistoPhiMax();  Float_t etamax = GetHistogramRanges()->GetHistoEtaMax();
  Float_t ptmin  = GetHistogramRanges()->GetHistoPtMin();  Float_t phimin = GetHistogramRanges()->GetHistoPhiMin();  Float_t etamin = GetHistogramRanges()->GetHistoEtaMin();	
  
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

//________________________________________
/// Initialize the parameters of the analysis.
//________________________________________
void AliAnaRandomTrigger::InitParameters()
{ 
  SetOutputAODClassName("AliCaloTrackParticleCorrelation");
  SetOutputAODName("RandomTrigger");

  AddToHistogramsName("AnaRandomTrigger_");
  
  fNRandom   = 1    ;
  fPhiCut[0] = 0.   ;
  fPhiCut[1] = TMath::TwoPi() ;
  fEtaCut[0] =-1.   ;
  fEtaCut[1] = 1.   ;
}

//_________________________________________________________
/// Print some relevant parameters set for the analysis.
//_________________________________________________________
void AliAnaRandomTrigger::Print(const Option_t * opt) const
{
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");	

  printf("Detector = %s\n",  fTriggerDetectorString.Data());
  printf("Min E   = %3.2f - Max E   = %3.2f\n", GetMinPt(), GetMaxPt());
  printf("Min Eta = %3.2f - Max Eta = %3.2f\n", fEtaCut[0], fEtaCut[1]);
  printf("Min Phi = %3.2f - Max Phi = %3.2f\n", fPhiCut[0], fPhiCut[1]);
} 

//______________________________________________
/// Do analysis and fill aods.
/// Generate particle randomly.
/// fNRandom particles per event.
//______________________________________________
void  AliAnaRandomTrigger::MakeAnalysisFillAOD()
{
  for(Int_t irandom = 0; irandom < fNRandom; irandom++)
  {
    // Get the random variables of the trigger
    Float_t pt  = -1000;
    Float_t eta = -1000;
    Float_t phi = -1000;
    Bool_t  in  = kFALSE;
    Bool_t  exc = kTRUE; 
    
//    printf("Random generation ranges: pt [%2.2f,%2.2f], eta [%2.2f,%2.2f], phi [%2.2f,%2.2f]\n",
//           GetMinPt(), GetMaxPt(),
//           fEtaCut[0], fEtaCut[1],
//           fPhiCut[0]*TMath::RadToDeg(), fPhiCut[1]*TMath::RadToDeg());
    
    // If excluded, generate a new trigger until accepted
    while ( exc || !in )
    {
      pt  = fRandom.Uniform(GetMinPt(), GetMaxPt());
      eta = fRandom.Uniform(fEtaCut[0], fEtaCut[1]);
      phi = fRandom.Uniform(fPhiCut[0], fPhiCut[1]);
      in  = kFALSE;
      exc = kTRUE; 
      
      //printf("Generated pt %2.2f, eta %2.2f, phi %2.2f\n",pt,eta,phi*TMath::RadToDeg());
      
      // Check if particle falls into a predefined/standard detector region
      if ( IsFiducialCutOn() )
      {
        in  = GetFiducialCut()->IsInFiducialCut(eta,phi,fTriggerDetector) ;
        //printf("\t in %s acceptance? %d\n",fTriggerDetectorString.Data(), in);
      }
      else
      {
        in = kTRUE;
      }
      
      // Check if particle falls into a dead region, if inside, get new
      if ( in )
      {
        exc = ExcludeDeadBadRegions(eta,phi);
        //printf("\t dead area? %d\n",exc);
      }
    }
    
    // Create the AOD trigger object
    fMomentum.SetPtEtaPhiM(pt,eta,phi,0);
    
    AliCaloTrackParticle trigger = AliCaloTrackParticle(fMomentum);
    trigger.SetDetectorTag(fTriggerDetector);
    trigger.SetSModNumber(GetModuleNumber(&trigger));
    
    AliDebug(1,Form("iRandom %d, Trigger e %2.2f pt %2.2f, phi %2.2f, eta %2.2f, SM %d",
                    irandom, trigger.E(), trigger.Pt(), trigger.Phi(), trigger.Eta(), trigger.GetSModNumber()));
    
    AddAODParticle(trigger);
  }
  
  AliDebug(1,Form("Final aod branch entries %d", GetOutputAODBranch()->GetEntriesFast()));
} 

//_____________________________________________________
// Fill control histograms with generated trigger kinematics.
//_____________________________________________________
void  AliAnaRandomTrigger::MakeAnalysisFillHistograms()
{
  // Loop on stored AODParticles
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  
  AliDebug(1,Form("AOD branch entries %d, fNRandom %d", naod, fNRandom));
  
  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    AliCaloTrackParticle* trigger =  (AliCaloTrackParticle*) (GetOutputAODBranch()->At(iaod));
    
    fhPt    ->Fill(trigger->Pt (),                 GetEventWeight());
    fhPhi   ->Fill(trigger->Pt (), trigger->Phi(), GetEventWeight());
    fhEta   ->Fill(trigger->Pt (), trigger->Eta(), GetEventWeight());
    fhEtaPhi->Fill(trigger->Eta(), trigger->Phi(), GetEventWeight());
  }// aod branch loop
}

//_________________________________________________________
/// Set the detrimeter for the analysis.
//_________________________________________________________
void AliAnaRandomTrigger::SetTriggerDetector(TString det)
{
  fTriggerDetectorString = det;
  
  if     (det=="EMCAL") fTriggerDetector = kEMCAL;
  else if(det=="PHOS" ) fTriggerDetector = kPHOS;
  else if(det=="CTS")   fTriggerDetector = kCTS;
  else if(det=="DCAL")  fTriggerDetector = kDCAL;
  else if(det.Contains("DCAL") && det.Contains("PHOS")) fTriggerDetector = kDCALPHOS;
  else AliFatal(Form("Detector < %s > not known!", det.Data()));
}

//______________________________________________________
// Set the calorimeter for the analysis.
//______________________________________________________
void AliAnaRandomTrigger::SetTriggerDetector(Int_t det)
{
  fTriggerDetector = det;
  
  if     ( det == kEMCAL   ) fTriggerDetectorString = "EMCAL";
  else if( det == kPHOS    ) fTriggerDetectorString = "PHOS";
  else if( det == kCTS     ) fTriggerDetectorString = "CTS";
  else if( det == kDCAL    ) fTriggerDetectorString = "DCAL";
  else if( det == kDCALPHOS) fTriggerDetectorString = "DCAL_PHOS";
  else AliFatal(Form("Detector < %d > not known!", det));
}


