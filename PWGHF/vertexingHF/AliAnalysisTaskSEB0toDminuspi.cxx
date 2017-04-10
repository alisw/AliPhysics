/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appeuear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//
//
//                  Base class for B0 analysis
//
//
//  Cuts arew centralized in AliRDHFCutsDplustoKpipi
//

#include <TSystem.h>
#include <TParticle.h>
#include <TH1I.h>
#include <TObjArray.h>
#include "TROOT.h"
#include <TDatabasePDG.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
//#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliVertex.h"
#include "AliVVertex.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliExternalTrackParam.h"
#include "AliNeutralTrackParam.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliNormalizationCounter.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskSEB0toDminuspi.h"
#include "AliVertexingHFUtils.h"

ClassImp(AliAnalysisTaskSEB0toDminuspi)

//__________________________________________________________________________
AliAnalysisTaskSEB0toDminuspi::AliAnalysisTaskSEB0toDminuspi():  
  AliAnalysisTaskSE(),
  fEvents(0),
  fOutput(0),
  fOutputBins(0),
  //fOutputPID(0),
  fArrayB0(0),
  //HISTO
  fCEvents(0),
  fB0InvMass(0), 
  fDplusInvMass(0),  
  fB0InvMass_PassCut1(0),   
  fB0InvMassMinusD_PassCut1(0),   
  fDplusInvMass_PassCut1(0), 
  fDplusInvMass_DplusMCmatch(0),
  fDplusInvMass_MCmatch(0),
  fBInvMass_MCmatch(0),
  fBInvMass_MCmatch_PassCut1(0),
  fB0InvMassMinusD_MCmatch_PassCut1(0),
  
  fpiPt(0),
  fpiEta(0),
  fpiPhi(0),
  fDplusPt(0),
  fDplusEta(0),
  fDplusPhi(0),
  fDplusPointing(0),
  fDplusDecayLength(0),
  fDplusNormalizedDecayLength(0),
  fDplusSigmaVert(0),
  fBPt(0),
  fBEta(0),
  fBPhi(0),
  fBPointing(0),
  fBPointingXY(0),
  fBDecayLength(0),
  fBNormalizedDecayLength(0),
  fDplusd0(0),
  fpid0(0),
  fproductd0(0),
  fpiPt_MCmatch(0),
  fpiEta_MCmatch(0),
  fpiPhi_MCmatch(0),
  fDplusPt_MCmatch(0),
  fDplusEta_MCmatch(0),
  fDplusPhi_MCmatch(0),
  fDplusPointing_MCmatch(0),
  fDplusDecayLength_MCmatch(0),
  fDplusNormalizedDecayLength_MCmatch(0),
  fDplusSigmaVert_MCmatch(0),
  fBPt_MCmatch(0),
  fBEta_MCmatch(0),
  fBPhi_MCmatch(0),
  fBPointing_MCmatch(0),
  fBPointingXY_MCmatch(0),
  fBDecayLength_MCmatch(0),
  fBNormalizedDecayLength_MCmatch(0),
  fDplusd0_MCmatch(0),
  fpid0_MCmatch(0),
  fproductd0_MCmatch(0),

  
  fd0MMExpD(0),
  fd0MMExpDMCmatch(0),
  fd0MMExppi(0),
  fd0MMExppiMCmatch(0),
  fd0MMExpDDaughters(0),
  fd0MMExpDDaughters_MCmatch(0),
  fcosoa(0),
  fcosoa_MCmatch(0),

  fBDplusPt(0),
  fBDplusPt_MCmatch(0),
  fBpiPt(0),
  fBpiPt_MCmatch(0),
  
  
  fCounter(0),
  fCutsDistr(kFALSE),
  fListCuts(0),
  fEtaSelection(0),
  fReadMC(kFALSE),
  fUseBit(kTRUE),
  fUseQuarkTagInKine(kTRUE),
  fRDCutsAnalysis(0),
  kbins(kTRUE),
  fInvertedTopomaticCutOnDDaughters(0),
  fpTD(0),
  fpTpi(0), 
  fprodd0(999), 
  fcosB(0),
  fcosBXY(0),
  fdlB(0),
  fNdlBXY(0),
  fTopomaticD(999),
  fTopomaticpi(999),
  fcosoaDpi(-999),
  
  fDplusMassLowLimit(1.845),
  fDplusMassUpLimit(1.895),
  fBMassLowLimit(4.),
  fBMassUpLimit(7),
  fBinWidth(0.002)

  
  

  //
  // Default ctor
  //
  
{
    
  for(Int_t i=0; i<kMaxPtBins; i++){
    fB0InvMassBin[i]=0;
    fDplusInvMassBin[i]=0;
    fB0InvMassBin_MCmatch[i]=0;
    fB0InvMassBin_MCmatch_PassCut1[i]=0;
    fDplusInvMassBin_PassCut1[i]=0;
    fB0InvMassBin_MCmatch_PassCut1[i]=0;
    fDplusInvMassBin_MCmatch[i]=0;
    fB0InvMassBin_PassCut1[i]=0;
    fDplusInvMassBin_PassCut1[i]=0;
    fB0InvMassBinMinusD_PassCut1[i]=0;
    fB0InvMassBinMinusD_MCmatch_PassCut1[i]=0;
    
    fBDecayLengthBin[i]=0;
    fBNormalizedDecayLengthBin[i]=0;
    fBPointingBin[i]=0;
    fBPointingXYBin[i]=0;
    fDplusd0Bin[i]=0;
    fBpTBin[i]=0;
    fDpluspTBin[i]=0;
    fproductd0Bin[i]=0;
    fDplusPointingBin[i]=0;
    fDplusDecayLengthBin[i]=0;
    fDplusNormalizedDecayLengthBin[i]=0;
    fDplusSigmaVertBin[i]=0;
      
    fBDecayLengthBin_MCmatch[i]=0;
    fBNormalizedDecayLengthBin_MCmatch[i]=0;
    fBPointingBin_MCmatch[i]=0;
    fBPointingXYBin_MCmatch[i]=0;
    fDplusd0Bin_MCmatch[i]=0;
    fBpTBin_MCmatch[i]=0;
    fDpluspTBin_MCmatch[i]=0;
    fproductd0Bin_MCmatch[i]=0;
    fDplusPointingBin_MCmatch[i]=0;
    fDplusDecayLengthBin_MCmatch[i]=0;
    fDplusNormalizedDecayLengthBin_MCmatch[i]=0;
    fDplusSigmaVertBin_MCmatch[i]=0;
    
    fd0MMExpDBin[i]=0;               
    fd0MMExpDBin_MCmatch[i]=0;              
    fd0MMExppiBin[i]=0;            
    fd0MMExppiBin_MCmatch[i]=0;
    fd0MMExpDDaughtersBin[i]=0;
    fd0MMExpDDaughtersBin_MCmatch[i]=0;



    fBDplusPtBin[i]=0;
    fBDplusPtBin_MCmatch[i]=0;
    fBpiPtBin[i]=0;
    fBpiPtBin_MCmatch[i]=0;

 
  
  }
  
}
//___________________________________________________________________________
AliAnalysisTaskSEB0toDminuspi::AliAnalysisTaskSEB0toDminuspi(const Char_t* name, AliRDHFCutsDplustoKpipi* dpluscutsana) :
  AliAnalysisTaskSE(name),
  fEvents(0),
  fOutput(0),
  fOutputBins(0),
  //fOutputPID(0),
  //fCuts(0),
  fArrayB0(0),
  //HISTO
  fCEvents(0),
  fB0InvMass(0),    
  fDplusInvMass(0),
  fB0InvMass_PassCut1(0),
  fB0InvMassMinusD_PassCut1(0),   
  fDplusInvMass_PassCut1(0),  
  fDplusInvMass_DplusMCmatch(0),
  fBInvMass_MCmatch(0),
  fBInvMass_MCmatch_PassCut1(0),
  fB0InvMassMinusD_MCmatch_PassCut1(0),
 
  fpiPt(0),
  fpiEta(0),
  fpiPhi(0),
  fDplusPt(0),
  fDplusEta(0),
  fDplusPhi(0),
  fDplusPointing(0),
  fDplusDecayLength(0),
  fDplusNormalizedDecayLength(0),
  fDplusSigmaVert(0),
  fBPt(0),
  fBEta(0),
  fBPhi(0),
  fBPointing(0),
  fBPointingXY(0),
  fBDecayLength(0),
  fBNormalizedDecayLength(0),
  fDplusd0(0),
  fpid0(0),
  fproductd0(0),
  fpiPt_MCmatch(0),
  fpiEta_MCmatch(0),
  fpiPhi_MCmatch(0),
  fDplusPt_MCmatch(0),
  fDplusEta_MCmatch(0),
  fDplusPhi_MCmatch(0),
  fDplusPointing_MCmatch(0),
  fDplusDecayLength_MCmatch(0),
  fDplusNormalizedDecayLength_MCmatch(0),
  fDplusSigmaVert_MCmatch(0),
  fBPt_MCmatch(0),
  fBEta_MCmatch(0),
  fBPhi_MCmatch(0),
  fBPointing_MCmatch(0),
  fBPointingXY_MCmatch(0),
  fBDecayLength_MCmatch(0),
  fBNormalizedDecayLength_MCmatch(0),
  fDplusd0_MCmatch(0),
  fpid0_MCmatch(0),
  fproductd0_MCmatch(0),
  fDplusInvMass_MCmatch(0),


  fd0MMExpD(0),
  fd0MMExpDMCmatch(0),
  fd0MMExppi(0),
  fd0MMExppiMCmatch(0),
  fd0MMExpDDaughters(0),
  fd0MMExpDDaughters_MCmatch(0),
  fcosoa(0),
  fcosoa_MCmatch(0),

  fBDplusPt(0),
  fBDplusPt_MCmatch(0),
  fBpiPt(0),
  fBpiPt_MCmatch(0),

  fCounter(0),
  fCutsDistr(kFALSE),
  fListCuts(0),
  fEtaSelection(0),
  fReadMC(kFALSE),
  fUseBit(kTRUE),
  fUseQuarkTagInKine(kTRUE),
  fRDCutsAnalysis(dpluscutsana),
  kbins(kTRUE),
  fInvertedTopomaticCutOnDDaughters(0),
  fpTD(0),
  fpTpi(0), 
  fprodd0(999), 
  fcosB(0),
  fcosBXY(0),
  fdlB(0),
  fNdlBXY(0),
  fTopomaticD(999),
  fTopomaticpi(999),
  fcosoaDpi(-999),
 
  fDplusMassLowLimit(1.845),
  fDplusMassUpLimit(1.895),
  fBMassLowLimit(4),
  fBMassUpLimit(7),
  fBinWidth(0.002)
  


{  

  //standard contructor
  for(Int_t i=0; i<kMaxPtBins; i++){
    fB0InvMassBin[i]=0;
    fDplusInvMassBin[i]=0;
    fB0InvMassBin_MCmatch[i]=0;
    fB0InvMassBin_MCmatch_PassCut1[i]=0;
    fDplusInvMassBin_MCmatch[i]=0;
    fB0InvMassBin_PassCut1[i]=0;
    fDplusInvMassBin_PassCut1[i]=0;
    fB0InvMassBinMinusD_PassCut1[i]=0;
    fB0InvMassBinMinusD_MCmatch_PassCut1[i]=0;
    
    fBDecayLengthBin[i]=0;
    fBNormalizedDecayLengthBin[i]=0;
    fBPointingBin[i]=0;
    fBPointingXYBin[i]=0;
    fDplusd0Bin[i]=0;
    fBpTBin[i]=0;
    fDpluspTBin[i]=0;
    fproductd0Bin[i]=0;
    fDplusPointingBin[i]=0;
    fDplusDecayLengthBin[i]=0;
    fDplusNormalizedDecayLengthBin[i]=0;
    fDplusSigmaVertBin[i]=0;
      
    fBDecayLengthBin_MCmatch[i]=0;
    fBNormalizedDecayLengthBin_MCmatch[i]=0;
    fBPointingBin_MCmatch[i]=0;
    fBPointingXYBin_MCmatch[i]=0;
    fDplusd0Bin_MCmatch[i]=0;
    fBpTBin_MCmatch[i]=0;
    fDpluspTBin_MCmatch[i]=0;
    fproductd0Bin_MCmatch[i]=0;
    fDplusPointingBin_MCmatch[i]=0;
    fDplusDecayLengthBin_MCmatch[i]=0;
    fDplusNormalizedDecayLengthBin_MCmatch[i]=0;
    fDplusSigmaVertBin_MCmatch[i]=0;
    
    fd0MMExpDBin[i]=0;               
    fd0MMExpDBin_MCmatch[i]=0;              
    fd0MMExppiBin[i]=0;            
    fd0MMExppiBin_MCmatch[i]=0;
    fd0MMExpDDaughtersBin[i]=0;
    fd0MMExpDDaughtersBin_MCmatch[i]=0;

 
    fBDplusPtBin[i]=0;
    fBDplusPtBin_MCmatch[i]=0;
    fBpiPtBin[i]=0;
    fBpiPtBin_MCmatch[i]=0;
    
  
  }  
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSEB0toDminuspi","Calling Constructor");

  fRDCutsAnalysis=dpluscutsana;

  DefineOutput(1,TList::Class());  //conters
  DefineOutput(2,TList::Class());  //All Entries output
  DefineOutput(3,AliNormalizationCounter::Class());   // normalization
  DefineOutput(4,AliRDHFCutsDplustoKpipi::Class());   //My private output
}
//_________________________________________________________________
void  AliAnalysisTaskSEB0toDminuspi::SetBMassLimits(Float_t lowlimit, Float_t uplimit){
  /// set invariant mass limits
  if(uplimit>lowlimit)
    {
      Float_t bw=GetBinWidth();
      fBMassUpLimit=uplimit;
      fBMassLowLimit= lowlimit;
      SetBinWidth(bw);
    }
}
//________________________________________________________________
void AliAnalysisTaskSEB0toDminuspi::SetBinWidth(Float_t w){
  Float_t width=w;
  Int_t nbins=(Int_t)((fBMassUpLimit-fBMassLowLimit)/width+0.5);
  Int_t missingbins=4-nbins%4;
  nbins=nbins+missingbins;
  width=(fBMassUpLimit-fBMassLowLimit)/nbins;
  if(missingbins!=0){
    printf("AliAnalysisTaskSEDplus::SetBinWidth: W-bin width of %f will produce histograms not rebinnable by 4. New width set to %f\n",w,width);
  }
  else{
    if(fDebug>1) printf("AliAnalysisTaskSEDplus::SetBinWidth: width set to %f\n",width);
  }
  fBinWidth=width;
}
//_________________________________________________________________
Int_t AliAnalysisTaskSEB0toDminuspi::GetNBinsHistos(){
  return (Int_t)((fBMassUpLimit-fBMassLowLimit)/fBinWidth+0.5);
}
//___________________________________________________________________________
AliAnalysisTaskSEB0toDminuspi::~AliAnalysisTaskSEB0toDminuspi() {
  //
  // destructor
  //
  Info("~AliAnalysisTaskSEB0toDminuspi","Calling Destructor");
  
  delete fOutput;
  delete fOutputBins;
  //delete fOutputPID;
  //delete fCuts;
  delete fListCuts;
  delete fRDCutsAnalysis;
  //~fArrayB0();
  //HISTOs
  delete fCEvents;
  delete fB0InvMass; 
  delete fDplusInvMass;
  delete fB0InvMass_PassCut1;
  delete fB0InvMassMinusD_PassCut1;   
  delete fDplusInvMass_PassCut1; 
  delete fDplusInvMass_DplusMCmatch;
  delete fDplusInvMass_MCmatch;  
  delete fBInvMass_MCmatch;  
  delete fBInvMass_MCmatch_PassCut1;  
  delete fB0InvMassMinusD_MCmatch_PassCut1;
  delete fpiPt;
  delete fpiEta;
  delete fpiPhi;
  delete fDplusPt;
  delete fDplusEta;
  delete fDplusPhi;
  delete fDplusPointing;
  delete fDplusDecayLength;
  delete fDplusNormalizedDecayLength;
  delete fDplusSigmaVert;
  delete fBPt;
  delete fBEta;
  delete fBPhi;
  delete fBPointing;
  delete fBPointingXY;
  delete fBDecayLength;
  delete fBNormalizedDecayLength;
  delete fDplusd0;
  delete fpid0;
  delete fproductd0;
  delete fpiPt_MCmatch;
  delete fpiEta_MCmatch;
  delete fpiPhi_MCmatch;
  delete fDplusPt_MCmatch;
  delete fDplusEta_MCmatch;
  delete fDplusPhi_MCmatch;
  delete fDplusPointing_MCmatch;
  delete fDplusDecayLength_MCmatch;
  delete fDplusNormalizedDecayLength_MCmatch;
  delete fDplusSigmaVert_MCmatch;  
  delete fBPt_MCmatch;
  delete fBEta_MCmatch;
  delete fBPhi_MCmatch;
  delete fBPointing_MCmatch;
  delete fBPointingXY_MCmatch;
  delete fBDecayLength_MCmatch;
  delete fBNormalizedDecayLength_MCmatch;
  delete fDplusd0_MCmatch;
  delete fpid0_MCmatch;  
  delete fproductd0_MCmatch;
   

  delete fd0MMExpD;
  delete fd0MMExpDMCmatch;
  delete fd0MMExppi;
  delete fd0MMExppiMCmatch;
  delete fd0MMExpDDaughters;
  delete fd0MMExpDDaughters_MCmatch;
  delete fcosoa;
  delete fcosoa_MCmatch;
  delete fBDplusPt;
  delete fBDplusPt_MCmatch;
  delete fBpiPt;
  delete fBpiPt_MCmatch;

  
  for(Int_t i=0; i<kMaxPtBins; i++){
    delete fB0InvMassBin[i];
    delete fDplusInvMassBin[i];
    delete fB0InvMassBin_MCmatch[i];
    delete fDplusInvMassBin_MCmatch[i];
    delete fB0InvMassBin_PassCut1[i];
    delete fB0InvMassBin_MCmatch_PassCut1[i];
    delete fDplusInvMassBin_PassCut1[i]; 
    delete fB0InvMassBinMinusD_PassCut1[i];
    delete fB0InvMassBinMinusD_MCmatch_PassCut1[i];
 
   
    delete fBDecayLengthBin[i];
    delete fBNormalizedDecayLengthBin[i];
    delete fBPointingBin[i];
    delete fBPointingXYBin[i];   
    delete fDplusd0Bin[i];
    delete fBpTBin[i];
    delete fDpluspTBin[i];
    delete fproductd0Bin[i];
    delete fDplusPointingBin[i];
    delete fDplusDecayLengthBin[i];
    delete fDplusNormalizedDecayLengthBin[i];
    delete fDplusSigmaVertBin[i];

    delete fBDecayLengthBin_MCmatch[i];
    delete fBNormalizedDecayLengthBin_MCmatch[i];
    delete fBPointingBin_MCmatch[i];
    delete fBPointingXYBin_MCmatch[i];
    delete fDplusd0Bin_MCmatch[i];
    delete fBpTBin_MCmatch[i];
    delete fDpluspTBin_MCmatch[i];
    delete fproductd0Bin_MCmatch[i];
    delete fDplusPointingBin_MCmatch[i];
    delete fDplusDecayLengthBin_MCmatch[i];
    delete fDplusNormalizedDecayLengthBin_MCmatch[i];
    delete fDplusSigmaVertBin_MCmatch[i];


    delete fd0MMExpDBin[i];
    delete fd0MMExpDBin_MCmatch[i];
    delete fd0MMExppiBin[i];
    delete fd0MMExppiBin_MCmatch[i];
    delete fd0MMExpDDaughtersBin[i];
    delete fd0MMExpDDaughtersBin_MCmatch[i];

    delete fBDplusPtBin[i];
    delete fBDplusPtBin_MCmatch[i];
    delete fBpiPtBin[i];
    delete fBpiPtBin_MCmatch[i];

   
  }
}
//_________________________________________________
void AliAnalysisTaskSEB0toDminuspi::Init(){
  //
  // Initialization
  //

  if(fDebug > 1) printf("AnalysisTaskSEB0toDminuspi::Init() \n");
  fListCuts=new TList();
  fListCuts->SetOwner();
  AliRDHFCutsDplustoKpipi *analysis = new AliRDHFCutsDplustoKpipi(*fRDCutsAnalysis);
  analysis->SetName("AnalysisCuts");
  
  fListCuts->Add(analysis);
  
  PostData(4,fListCuts);
  

  return;
}

//_________________________________________________
void AliAnalysisTaskSEB0toDminuspi::UserExec(Option_t *)
{
  // user exec
  if (!fInputEvent) {
    Error("UserExec","NO EVENT FOUND!");
    return;
  }
  
  fEvents++;
  
  TClonesArray *array3Prong = 0;
  
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  
  Double_t bz= (Double_t)aodEvent->GetMagneticField(); 
  
  fCEvents->Fill(1);

  //Extract the Charm3Prong object from the AOD
  
  if(!aodEvent && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      array3Prong=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
    }
  } else {
    array3Prong=(TClonesArray*)aodEvent->GetList()->FindObject("Charm3Prong");
  }
  
  
  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return;
   fCEvents->Fill(2);
 
  
  if (array3Prong->GetEntriesFast()==0) return;
  
  Bool_t isEvSel=fRDCutsAnalysis->IsEventSelected(aodEvent);


 TClonesArray *mcArray = 0; 
  AliAODMCHeader *mcHeader=0;
  
  if(fReadMC){
   mcArray = (TClonesArray*)aodEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName());
   if(!mcArray) {
     printf("AliAnalysisTaskSEDplus::UserExec: MC particles branch not found!\n");
     return;
   }
   // load MC header
    mcHeader =  (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSEDplus::UserExec: MC header branch not found!\n");
      return;
    }
  }

  //  if(fReadMC && fStepMCAcc)FillMCAcceptanceHistos(arrayMC, mcHeader);


  
  if(!isEvSel) return;
  fCounter->StoreEvent(aodEvent,fRDCutsAnalysis,fReadMC);
   
  fCEvents->Fill(3);
  // counters for efficiencies
  //Int_t icountReco = 0;
  Int_t nOS=0;
  
  //B0 and D+ prongs needed to MatchToMC and MatchToMCB method
  Int_t pdgDgDplustoKpipi[3]={321,211,211};
  Int_t pdgDgB0toDminuspi[2]={411,211};
  
  // AOD primary vertex
  AliAODVertex *primaryVertex = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if(!primaryVertex) return;
  if(primaryVertex->GetNContributors()<1) return;
  
  if (!array3Prong){
    AliInfo("Could not find array of Charm3Prong, skipping the event");
    return;
  }
  else AliDebug(2, Form("Found %d vertices",array3Prong->GetEntriesFast())); 
  
  Int_t nSelectedAna =0;
  Int_t nSelectedProd =0;
	
  // loop over the D+ candidates and reconstruct the secondary vertex //
  //                                                                  //
  //------------------------------------------------------------------//

 
  Double_t xdummy=0.,ydummy=0.,dca,dcaB, e[2];
  Double_t d0z0B[2],covd0z0B[3],d0B[2],d0errB[2],d0z0[2],covd0z0[3],d0[2],d0err[2];
  
  AliAODRecoDecayHF2Prong *B0Dminuspi =0x0;
  
  Int_t n3Prong = array3Prong->GetEntriesFast();  
  AliAODRecoDecayHF3Prong *theDminusparticle =0x0;
  Int_t labDp=-1;		//if Dminus is matched with MC
    
  Int_t B0evtcount =0;  
  
  //Loop over the Dplus candidates
  
  for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
  // D+ candidate
    theDminusparticle = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);
    fCEvents->Fill(4);
  
    if(fUseBit && !theDminusparticle->HasSelectionBit(AliRDHFCuts::kDplusCuts)){
	   fCEvents->Fill(5);
  	   continue;
	}
    Int_t passTopolAndPIDCuts=fRDCutsAnalysis->IsSelected(theDminusparticle,AliRDHFCuts::kAll,aodEvent);

      if(!fRDCutsAnalysis->GetIsSelectedCuts()) {
	 fCEvents->Fill(6);
      continue;
      }
      if(!fRDCutsAnalysis->GetIsSelectedPID()) {
	   fCEvents->Fill(6);
        continue;
      }

   
      Double_t etaD=theDminusparticle->Eta();
      Double_t phiD=theDminusparticle->Phi();
      if(fEtaSelection!=0){
	if(fEtaSelection==1 && etaD<0)	 continue;
	if(fEtaSelection==-1 && etaD>0)  continue;
      }
      
      

      Bool_t unsetvtx=kFALSE;
      if(!theDminusparticle->GetOwnPrimaryVtx()){
	    theDminusparticle->SetOwnPrimaryVtx(primaryVertex);
	    unsetvtx=kTRUE;
      }

      Double_t ptCand = theDminusparticle->Pt();
    
      Bool_t recVtx=kFALSE;
      AliAODVertex *origownvtx=0x0;
      if(fRDCutsAnalysis->GetIsPrimaryWithoutDaughters()){
	    if(theDminusparticle->GetOwnPrimaryVtx()) origownvtx = new AliAODVertex(*theDminusparticle->GetOwnPrimaryVtx());	
	    if(fRDCutsAnalysis->RecalcOwnPrimaryVtx(theDminusparticle,aodEvent)) recVtx=kTRUE;
	    else fRDCutsAnalysis->CleanOwnPrimaryVtx(theDminusparticle,aodEvent,origownvtx);
      }
    

      Bool_t isPrimary=kFALSE;
      Bool_t isFeeddown=kFALSE;
      Float_t pdgCode=-2;
      Float_t trueImpParXY=0.;
      Double_t ptB=-1.5;
      Int_t orig=0;
      
      if(TMath::Abs(theDminusparticle->Charge())!=1){
        fCEvents->Fill(7);
        continue;
      }
      
      if(fReadMC){
  	    labDp = theDminusparticle->MatchToMC(411,mcArray,3,pdgDgDplustoKpipi);
	    if(labDp>=0){
	      AliAODMCParticle *partDp = (AliAODMCParticle*)mcArray->At(labDp);
	      if(!partDp)continue;
	      orig=AliVertexingHFUtils::CheckOrigin(mcArray,partDp,fUseQuarkTagInKine);//Prompt = 4, FeedDown = 5
	      pdgCode=TMath::Abs(partDp->GetPdgCode());
	      if(orig==4){
	        isPrimary=kTRUE;
	        isFeeddown=kFALSE;
	        fCEvents->Fill(8);
	        AliAODMCParticle *partMomDp = (AliAODMCParticle*)mcArray->At(partDp->GetMother());
	      }else if(orig==5){
                isPrimary=kFALSE;
	        isFeeddown=kTRUE;
	        AliAODMCParticle *partMomDp = (AliAODMCParticle*)mcArray->At(partDp->GetMother());
		if(!partMomDp)continue;
		fCEvents->Fill(9); 
		// AliAODMCParticle *partGranMomDp = (AliAODMCParticle*)mcArray->At(partMomDp->GetMother());	        
	        if(TMath::Abs(partMomDp->GetPdgCode())==511 && ((partMomDp->GetLastDaughter()-partMomDp->GetFirstDaughter())==1)){
		  AliAODMCParticle* B0Dau2=0;
		  if(labDp== (partMomDp->GetFirstDaughter()-1)) B0Dau2 = (AliAODMCParticle*)mcArray->At(partMomDp->GetLastDaughter()-1);
	          else if(labDp== (partMomDp->GetLastDaughter()-1)) B0Dau2 = (AliAODMCParticle*)mcArray->At(partMomDp->GetFirstDaughter()-1);
	          if (B0Dau2) {if(TMath::Abs(B0Dau2->GetPdgCode())==211) fCEvents->Fill(10);}
		
		}
	        ptB=AliVertexingHFUtils::GetBeautyMotherPt(mcArray,partDp);
	      }else{
	        pdgCode=-3;
	      }
	    }
	    else{
	      pdgCode=-1;
	    }
      }



      


      
      Double_t invMass=theDminusparticle->InvMassDplus();
      Double_t rapid=theDminusparticle->YDplus();

      Bool_t isFidAcc=fRDCutsAnalysis->IsInFiducialAcceptance(ptCand,rapid);
      if(!isFidAcc){fCEvents->Fill(19);}//continue;}


   
      
      

      fCEvents->Fill(11);
  

      if(recVtx)fRDCutsAnalysis->CleanOwnPrimaryVtx(theDminusparticle,aodEvent,origownvtx);
      
      if(unsetvtx) theDminusparticle->UnsetOwnPrimaryVtx();

    

      //D particle treated as an external track to create the pairs D-pion
      AliExternalTrackParam *trackDminus = new AliExternalTrackParam();
      trackDminus->CopyFromVTrack(theDminusparticle);
      fDplusInvMass->Fill(theDminusparticle->InvMassDplus());
      if(labDp>=0)fDplusInvMass_DplusMCmatch->Fill(theDminusparticle->InvMassDplus());
      
    
      // get the ID of the D+ daughters
      Double_t Id_prong0 = theDminusparticle->GetProngID(0);
      Double_t Id_prong1 = theDminusparticle->GetProngID(1);
      Double_t Id_prong2 = theDminusparticle->GetProngID(2);
    
      //TObjArray b012;
    

    
    // ----------------------------------------------- ------------------------------------------- //
    //                                   B0 reconstruction                                         //
    //---------------------------------------------------------------------------------------------//
    
    Int_t iB0count=0;
    fArrayB0.clear();
    Int_t arrayMCB0label[50];    

    // select a region around Dplus mass (1.869), ~ 3 sigma of the integrated D invariant mass distribution
    //  if(theDminusparticle->InvMassDplus()<(1.895) && theDminusparticle->InvMassDplus()>(1.845)){
    if( theDminusparticle->InvMassDplus()>(fDplusMassLowLimit) && theDminusparticle->InvMassDplus()<(fDplusMassUpLimit)){
    
      // loop on all the track of the events to be paired with the D 
      for (Int_t i=0; i<aodEvent->GetNumberOfTracks(); i++){
           
	
	AliAODTrack* aodTrack = (AliAODTrack*)aodEvent->GetTrack(i);
	// remove the Dminus prongs
	if(aodTrack->GetID() == Id_prong0 ||aodTrack->GetID() == Id_prong1 || aodTrack->GetID() == Id_prong2) continue;
	
	fCEvents->Fill(21);
	

	if(!(aodTrack->TestFilterMask(BIT(4)))){fCEvents->Fill(20); continue;}
	
	// remove very low pt tracks
	if(aodTrack->Pt()<0.5) continue;

	//salections on pion tracks



		       
	if (TMath::Abs(aodTrack->Eta())>0.8)continue; // cut in eta pi 

	//check of the track PID, 3 sigma on TPC and TOF
	Double_t nSigmaTPCpi=fRDCutsAnalysis->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(aodTrack,AliPID::kPion);
	Double_t nSigmaTOFpi=fRDCutsAnalysis->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(aodTrack,AliPID::kPion);
	if(nSigmaTPCpi<-3 || nSigmaTPCpi>3 ){fCEvents->Fill(12);
	  continue;}
	if(nSigmaTOFpi<-3 || nSigmaTOFpi>3){
	  fCEvents->Fill(12);
	  continue;}
	//remove pair with D and pi with the same charge
	if(aodTrack->Charge()==theDminusparticle->Charge()){
	  fCEvents->Fill(13);
	  continue;
	}
	
	fCEvents->Fill(14);
	
	// reconstruct the new vertex of the 2 prongs
	AliExternalTrackParam et4; 
	et4.CopyFromVTrack(static_cast<AliAODTrack*>(aodEvent->GetTrack(i)));
	
	// here dca, d0 and d0rr
	trackDminus->PropagateToDCA(primaryVertex,bz,100.,d0z0B,covd0z0B);
	d0B[0]=d0z0B[0];
	d0errB[0] = TMath::Sqrt(covd0z0B[0]);
	et4.PropagateToDCA(primaryVertex,bz,100.,d0z0B,covd0z0B);
	d0B[1]=d0z0B[0];
	d0errB[1] = TMath::Sqrt(covd0z0B[0]);   
	
	dcaB=et4.GetDCA(trackDminus,bz,xdummy,ydummy);
	
	TObjArray b012_B;	      	     
	b012_B.Add(trackDminus); b012_B.Add(&et4); 
	
	// B0 vertex
	AliAODVertex *vtxNewB0 = RecalculateVertex(primaryVertex,&b012_B,bz);
	if(!vtxNewB0) {
	 continue;}
	
	
	// fill the momenta of the prongs
	Double_t pxB[2],pyB[2],pzB[2];
	pxB[0] = trackDminus-> Px();
	pyB[0] = trackDminus-> Py();
	pzB[0] = trackDminus-> Pz();
	pxB[1] = aodTrack-> Px();
	pyB[1] = aodTrack-> Py();
	pzB[1] = aodTrack-> Pz();
	
	UShort_t idB[2];
	//Int_t prongs_ID[2];
	
	// make the vertex of the B0 
	
	idB[0]= 0; 
	idB[1]= aodTrack->GetID();
	
	// thid is the D* vertex
	B0Dminuspi = new AliAODRecoDecayHF2Prong(vtxNewB0,pxB,pyB,pzB,d0B,d0errB,dcaB);
	if(!B0Dminuspi) continue;
	

	//add the daughter tracks to the B vertex 
	B0Dminuspi->GetSecondaryVtx()->AddDaughter(theDminusparticle);
	B0Dminuspi->GetSecondaryVtx()->AddDaughter(aodTrack);
	B0Dminuspi->SetPrimaryVtxRef((AliAODVertex*)aodEvent->GetPrimaryVertex());
	B0Dminuspi->SetProngIDs(2,idB);
    
	
	
	fCEvents->Fill(15);
	
	fArrayB0.resize(iB0count+1);
	fArrayB0[iB0count]=B0Dminuspi;
	
	
	
	Int_t isB0 = 0; 
	arrayMCB0label[iB0count]=-1;
	
	if(fReadMC){

	  // find associated MC particle for B0 -> pi DminustoKpipi	  
	  arrayMCB0label[iB0count] = B0Dminuspi->MatchToMCB3Prong(511,411, pdgDgB0toDminuspi, pdgDgDplustoKpipi, mcArray);
	  if(arrayMCB0label[iB0count]>(-1)){	    
	    isB0 = 1;

	   
	    fCEvents->Fill(16);
	    //	   fd0MMExpDDaughters_MCmatch->Fill(dd0min);

	  }      
	}

    
    
    
	
	nSelectedProd++;
	nSelectedAna++;
       	

	iB0count++;

       	
      

      }
      
   

    }
  
    if(iB0count>0) FillSpectrum(iB0count, arrayMCB0label,bz);
    
  
  fArrayB0.clear();
  
  B0evtcount+=iB0count;
  
  
  }
  
  fCounter->StoreCandidates(aodEvent,nSelectedProd,kTRUE);  
  fCounter->StoreCandidates(aodEvent,nSelectedAna,kFALSE); 

  
  AliDebug(2, Form("Found %i Reco particles that are D+!!",nOS));
  AliDebug(2, Form("Found %i Reco particles that are B0!!",B0evtcount));
  
  PostData(1,fOutput);
  PostData(2,fOutputBins);
  PostData(3,fCounter);

}
//________________________________________ terminate ___________________________
void AliAnalysisTaskSEB0toDminuspi::Terminate(Option_t*)
{    
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  //Info("Terminate","");
  AliAnalysisTaskSE::Terminate();
  
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  
  fCEvents        = dynamic_cast<TH1F*>(fOutput->FindObject("fCEvents"));
 

   

  return;
}
//___________________________________________________________________________
void AliAnalysisTaskSEB0toDminuspi::UserCreateOutputObjects() { 
 // output
  Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
  
  //slot #1  
  //OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("PtIntegrated");

  fOutputBins = new TList();
  fOutputBins->SetOwner();
  fOutputBins->SetName("PtBins");


    
  // define histograms
  DefineHistograms();

 //Counter for Normalization
 fCounter = new AliNormalizationCounter(Form("%s",GetOutputSlot(3)->GetContainer()->GetName()));
 fCounter->Init();


  PostData(1,fOutput);
  PostData(2,fOutputBins);


  return;
}
//___________________________________ hiostograms _______________________________________
void  AliAnalysisTaskSEB0toDminuspi::DefineHistograms(){
  // Create histograms

  fCEvents = new TH1F("fCEvents","conter",30,0.5,30.5);
  fCEvents->GetXaxis()->SetBinLabel(1,"n Events");
  fCEvents->GetXaxis()->SetBinLabel(2,"w primary vtx");
  fCEvents->GetXaxis()->SetBinLabel(3,"selected/ w array3prong");
  fCEvents->GetXaxis()->SetBinLabel(4,"D+ candidates");
  fCEvents->GetXaxis()->SetBinLabel(5,"D+ rejected for DplusCuts");
  fCEvents->GetXaxis()->SetBinLabel(6,"D+ rejected for topological + PID selection ");
  fCEvents->GetXaxis()->SetBinLabel(7,"D+ rejected for charge");
  fCEvents->GetXaxis()->SetBinLabel(8,"MC D+ primary orig=4");
  fCEvents->GetXaxis()->SetBinLabel(9,"MC D+ feed-down  orig=5");
  fCEvents->GetXaxis()->SetBinLabel(10,"MC B->D pi");
  
  fCEvents->GetXaxis()->SetBinLabel(11,"D+ candidates to be associated with a pion");
  fCEvents->GetXaxis()->SetBinLabel(12,"pion cadidate rejected for PID");
  fCEvents->GetXaxis()->SetBinLabel(13,"pion candidate rejected because of charge");
  fCEvents->GetXaxis()->SetBinLabel(14,"pion candidate passing selections");

  fCEvents->GetXaxis()->SetBinLabel(15,"B0 candidate");
  fCEvents->GetXaxis()->SetBinLabel(16,"B0 candidate matched MC");
  fCEvents->GetXaxis()->SetBinLabel(17,"B0 passing set cut 1");
  fCEvents->GetXaxis()->SetBinLabel(18,"B0 passing set cut 1 matched MC");
  fCEvents->GetXaxis()->SetBinLabel(19,"D + not in fiducial accaptance");
  fCEvents->GetXaxis()->SetBinLabel(20,"aodtrack rejected for filterbit 4 ");
  fCEvents->GetXaxis()->SetBinLabel(21,"pion tracks before cuts");

  fCEvents->SetStats(kTRUE);
  fCEvents->GetXaxis()->SetTitle("1");
  fCEvents->GetYaxis()->SetTitle("counts");



  Int_t nbins=GetNBinsHistos();

  fB0InvMass = new TH1F("fB0InvMass", "B0 invariant mass;GeV", nbins, fBMassLowLimit, fBMassUpLimit);
  fB0InvMass->Sumw2();
  fB0InvMassMinusD_PassCut1 = new TH1F("fB0InvMassMinusD_PassCut1", "B0-D invariant mass;GeV", nbins, fBMassLowLimit-1.869, fBMassUpLimit-1.869);
  fB0InvMassMinusD_PassCut1->Sumw2();
  fDplusInvMass = new TH1F("fDplusInvMass", "Dplus invariant mass;GeV", 500, 1.75, 1.95 );
  fDplusInvMass->Sumw2();
  fB0InvMass_PassCut1 = new TH1F("fB0InvMass_PassCut1", "B0 invariant mass with Cut1 implemented;GeV", nbins, fBMassLowLimit, fBMassUpLimit);
  fB0InvMass_PassCut1->Sumw2();
  fB0InvMassMinusD_MCmatch_PassCut1 = new TH1F("fB0InvMassMinusD_MCmatch_PassCut1", "B0-D invariant mass;GeV", nbins, fBMassLowLimit-1.869, fBMassUpLimit-1.869);
  fB0InvMassMinusD_MCmatch_PassCut1->Sumw2();
  fDplusInvMass_PassCut1 = new TH1F("fDplusInvMass_PassCut1", "D+ invariant mass with Cut1 implemented;GeV", 500, 1.75, 1.95);
  fDplusInvMass_PassCut1->Sumw2();
  fDplusInvMass_DplusMCmatch = new TH1F("fDplusInvMass_DplusMCmatch", "Dplus invariant mass (Dplus matched);GeV", 500, 1.75, 1.95);
  fDplusInvMass_DplusMCmatch->Sumw2();



  fpiPt = new TH1F("fpiPt", "Pt of pion (from B0 selected)", 200, 0, 20);
  fpiEta = new TH1F("fpiEta", "Eta of pion (from B0 selected)", 200, -1, 1);
  fpiPhi = new TH1F("fpiPhi", "Phi of pion (from B0 selected)", 200, 0, 7);
  fDplusPt = new TH1F("fDplusPt", "Pt of Dplus (from B0 selected)",100, 0, 50);
  fDplusEta = new TH1F("fDplusEta", "Eta of Dplus (from B0 selected)", 200, -1, 1);
  fDplusPhi = new TH1F("fDplusPhi", "Phi of Dplus (from B0 selected)", 200, 0, 7);
  fDplusPointing = new TH1F("fDplusPointing", "CosPointingAngle of Dplus (from B0 selected)", 200, 0., 1.);
  fDplusDecayLength = new TH1F("fDplusDecayLength", "Decay Length of D+ (B0 selected)", 200, 0, 0.5); 
  fDplusNormalizedDecayLength = new TH1F("fDplusNormalizedDecayLength", "Normalized Decay Length of D+ (B0 selected)", 200, 0, 20);
  fDplusSigmaVert = new TH1F("fDplusSigmaVert", "Sigma vertex of D+ (B0 selected)", 200, 0, 0.5); 
  fBPt = new TH1F("fBPt", "Pt of B0 (B0 selected)", 100, 0, 50);
  fBEta = new TH1F("fBEta", "Eta of B0 (B0 selected)", 200, -1, 1);
   fBPhi = new TH1F("fBPhi", "Phi of B0 (B0 selected)", 200, 0, 7);
  fBPointing = new TH1F("fBPointing", "CosPointingAngle of B0 (B0 selected)", 200, 0., 1.);
  fBPointingXY = new TH1F("fBPointingXY", "CosPointingAngleXY of B0 (B0 selected)", 200, 0., 1.);
  fBDecayLength = new TH1F("fBDecayLength", "Decay Length of B0 (B0 selected)", 200, 0, 0.5);
  fBNormalizedDecayLength = new TH1F("fBNormalizedDecayLength", "Normalized Decay Length of B0 (B0 selected)", 200, 0, 20);
  fDplusd0 = new TH1F("fDplusd0", "|d0| Dplus to prim vtx(B0 selected)", 200, 0, 0.25);
  fpid0 = new TH1F("fpid0", "|d0| pion to prim vtx(B0 selected)", 200, 0, 0.25);
  fproductd0 = new TH1F("fproductd0", "d0 pion x d0 Dplus (B0 selected)", 200, -0.001, 0.001);
  fpiPt_MCmatch = new TH1F("fpiPt_MCmatch", "Pt of pion (MCmatched B0)", 200, 0, 20);
  fpiEta_MCmatch = new TH1F("fpiEta_MCmatch", "Eta of pion (MCmatched B0)", 200, -1, 1);
  fpiPhi_MCmatch = new TH1F("fpiPhi_MCmatch", "Phi of pion (MCmatched B0)", 200, 0, 7);
  fDplusPt_MCmatch = new TH1F("fDplusPt_MCmatch", "Pt of Dplus (MCmatched B0)", 100, 0, 50);
  fDplusEta_MCmatch = new TH1F("fDplusEta_MCmatch", "Eta of Dplus (MCmatched B0)", 200, -1, 1);
  fDplusPhi_MCmatch = new TH1F("fDplusPhi_MCmatch", "Phi of Dplus (MCmatched B0)", 200, 0, 7);
  fDplusPointing_MCmatch = new TH1F("fDplusPointing_MCmatch", "CosPointingAngle of Dplus (MCmatched B0)", 200, 0., 1.);
  fDplusDecayLength_MCmatch = new TH1F("fDplusDecayLength_MCmatch", "Decay Length of D+ (MCmatched B0)", 200, 0, 0.5);
  fDplusNormalizedDecayLength_MCmatch = new TH1F("fDplusNormalizedDecayLength_MCmatch", "Normalized Decay Length of D+ (MCmatched B0)", 200, 0, 20);
  fDplusSigmaVert_MCmatch = new TH1F("fDplusSigmaVert_MCmatch", "Sigma vertex of D+ (MC matched B0)", 200, 0, 0.5); 
  fBPt_MCmatch = new TH1F("fBPt_MCmatch", "Pt of B0 (MCmatched B0)", 100, 0, 50);
  fBEta_MCmatch = new TH1F("fBEta_MCmatch", "Eta of B0 (MCmatched B0)", 200, -1, 1);
  fBPhi_MCmatch = new TH1F("fBPhi_MCmatch", "Phi of B0 (MCmatched B0)", 200, 0, 7);
  fBPointing_MCmatch = new TH1F("fBPointing_MCmatch", "CosPointingAngle of B0 (MCmatched B0)", 200, 0., 1.);
  fBPointingXY_MCmatch = new TH1F("fBPointingXY_MCmatch", "CosPointingAngleXY of B0 (MC matched B0)", 200, 0., 1.);
  fBDecayLength_MCmatch = new TH1F("fBDecayLength_MCmatch", "Decay Length of B0 (MCmatched B0)", 200, 0, 0.5);
  fBNormalizedDecayLength_MCmatch = new TH1F("fBNormalizedDecayLength_MCmatch", "Normalized Decay Length of B0 (MCmatched B0)", 200, 0, 20);
  fDplusd0_MCmatch = new TH1F("fDplusd0_MCmatch", "|d0| Dplus to prim vtx(MCmatched B0)", 200, 0, 0.25);
  fpid0_MCmatch = new TH1F("fpid0_MCmatch", "|d0| pi to prim vtx(MCmatched B0)", 200, 0, 0.25);
  fproductd0_MCmatch = new TH1F("fproductd0_MCmatch", "d0 pion x d0 Dplus (MCmatched B0)", 200, -0.001, 0.001);
  fDplusInvMass_MCmatch = new TH1F("fDplusInvMass_MCmatch", "D+ invariant mass (MCmatched B0);GeV", 500, 1.75, 1.95);
  fBInvMass_MCmatch = new TH1F("fBInvMass_MCmatch", "B0 invariant mass (MCmatched B0);GeV", nbins, fBMassLowLimit, fBMassUpLimit);
  fBInvMass_MCmatch_PassCut1 = new TH1F("fBInvMass_MCmatch_PassCut1", "B0 invariant mass (MCmatched B0);GeV", nbins, fBMassLowLimit, fBMassUpLimit);
 
  fd0MMExpD = new TH1F("fd0MMExpD", "Topomatic variable for D", 200, -10.0,10.);
  fd0MMExpDMCmatch = new TH1F("fd0MMExpDMCmatch", "Topomatic variable for D from B matched to MC", 200, -10,10.); 
  fd0MMExppi = new TH1F("fd0MMExppi", "Topomatic variable for pi ", 200, -10.0,10.);
  fd0MMExppiMCmatch = new TH1F("fd0MMExppiMCmatch", "Topomatic variable for pi from B  matched to MC", 200, -10,10.);
  fd0MMExpDDaughters = new TH1F("fd0MMExpDDaughters", "Topomatic variable for D daughters", 200, -10.0,10.);
  fd0MMExpDDaughters_MCmatch = new TH1F("fd0MMExpDDaughters_MCmatch", "Topomatic variable for D daughters with B matched to MC", 200, -10,10.); 

  fcosoa= new TH1F("fcosoa", "cos opening angle", 200, -1.10,1.10);
  fcosoa_MCmatch= new TH1F("fcosoa_MCmatch", "cos opening angle MC", 200, -1.10,1.1);
 
  
  fBDplusPt = new TH2F("fBDplusPt", ";B pT;D pT", 100, 0, 50,100,0,50);
  fBpiPt = new TH2F("fBpiPt", ";B pT;pi pT", 100, 0, 50,100,0,50);
  fBDplusPt_MCmatch = new TH2F("fBDplusPt_MCmatch", " ;B pT;D pT", 100, 0, 50,100,0,50);
  fBpiPt_MCmatch = new TH2F("fBpiPt_MCmatch", " ;B pT;pi pT", 100, 0, 50,100,0,50);
 

  
  fOutput->Add(fCEvents);
  fOutput->Add(fB0InvMass);
  fOutput->Add(fDplusInvMass);
  fOutput->Add(fB0InvMass_PassCut1);
  fOutput->Add(fDplusInvMass_PassCut1);
  fOutput->Add(fB0InvMassMinusD_PassCut1);
  
  fOutput->Add(fpiPt);
  fOutput->Add(fpiEta);
  fOutput->Add(fpiPhi);
  fOutput->Add(fDplusPt);
  fOutput->Add(fDplusEta);
  fOutput->Add(fDplusPhi);
  fOutput->Add(fDplusPointing);
  fOutput->Add(fDplusDecayLength);
  fOutput->Add(fDplusNormalizedDecayLength);
  fOutput->Add(fDplusSigmaVert);
  fOutput->Add(fBPt);
  fOutput->Add(fBEta);
  fOutput->Add(fBPhi);
  fOutput->Add(fBPointing);
  fOutput->Add(fBPointingXY);
  fOutput->Add(fBDecayLength);
  fOutput->Add(fBNormalizedDecayLength);
  fOutput->Add(fDplusd0);
  fOutput->Add(fpid0);
  fOutput->Add(fproductd0);
  fOutput->Add(fd0MMExpD);
  fOutput->Add(fd0MMExppi);
  fOutput->Add(fd0MMExpDDaughters);
  fOutput->Add(fBDplusPt);
  fOutput->Add(fBpiPt);
  fOutput->Add(fcosoa);


  if(fReadMC){
    fOutput->Add(fDplusInvMass_MCmatch);
    fOutput->Add(fBInvMass_MCmatch);    
    fOutput->Add(fBInvMass_MCmatch_PassCut1);    
    fOutput->Add(fB0InvMassMinusD_MCmatch_PassCut1);    
    fOutput->Add(fDplusInvMass_DplusMCmatch);

    fOutput->Add(fpiPt_MCmatch);
    fOutput->Add(fpiEta_MCmatch);
    fOutput->Add(fpiPhi_MCmatch);
    fOutput->Add(fDplusPt_MCmatch);
    fOutput->Add(fDplusEta_MCmatch);
    fOutput->Add(fDplusPhi_MCmatch);
    fOutput->Add(fDplusPointing_MCmatch);
    fOutput->Add(fDplusDecayLength_MCmatch);
    fOutput->Add(fDplusNormalizedDecayLength_MCmatch);
    fOutput->Add(fDplusSigmaVert_MCmatch);  
    fOutput->Add(fBPt_MCmatch);
    fOutput->Add(fBEta_MCmatch);
    fOutput->Add(fBPhi_MCmatch);
    fOutput->Add(fBPointing_MCmatch);
    fOutput->Add(fBPointingXY_MCmatch);
    fOutput->Add(fBDecayLength_MCmatch);
    fOutput->Add(fBNormalizedDecayLength_MCmatch);
    fOutput->Add(fDplusd0_MCmatch);
    fOutput->Add(fpid0_MCmatch);
    fOutput->Add(fproductd0_MCmatch);
    fOutput->Add(fd0MMExpDMCmatch);
    fOutput->Add(fd0MMExppiMCmatch);
    fOutput->Add(fd0MMExpDDaughters_MCmatch);
    fOutput->Add(fBDplusPt_MCmatch);
    fOutput->Add(fBpiPt_MCmatch);
     fOutput->Add(fcosoa_MCmatch);

  }
  

  

  TString HistName;
  Int_t fNPtBins=fRDCutsAnalysis->GetNPtBins();
  
  
   for(Int_t i=0; i<fNPtBins; i++){
    HistName.Form("fB0InvMass_bin%d",i);
    fB0InvMassBin[i]= new TH1F(HistName.Data(), HistName.Data(), nbins, fBMassLowLimit, fBMassUpLimit);
    fB0InvMassBin[i]->Sumw2();
    HistName.Form("fDplusInvMass_bin%d",i);
    fDplusInvMassBin[i]= new TH1F(HistName.Data(), HistName.Data(), 500, 1.819, 1.919 );
    fDplusInvMassBin[i]->Sumw2();
    HistName.Form("fB0InvMass_bin%d_MCmatch",i);
    fB0InvMassBin_MCmatch[i]= new TH1F(HistName.Data(), HistName.Data(), nbins, fBMassLowLimit, fBMassUpLimit);
    fB0InvMassBin_MCmatch[i]->Sumw2();
     HistName.Form("fB0InvMass_bin%d_MCmatch_PassCut1",i);
    fB0InvMassBin_MCmatch_PassCut1[i]= new TH1F(HistName.Data(), HistName.Data(), nbins, fBMassLowLimit, fBMassUpLimit);
    fB0InvMassBin_MCmatch_PassCut1[i]->Sumw2();
 
    HistName.Form("fB0InvMassMinusD_bin%d_PassCut1",i);
    fB0InvMassBinMinusD_PassCut1[i]= new TH1F(HistName.Data(), HistName.Data(), nbins, fBMassLowLimit-1.869, fBMassUpLimit-1.869);
    fB0InvMassBinMinusD_PassCut1[i]->Sumw2();
    HistName.Form("fDplusInvMass_bin%d_MCmatch",i);
    fDplusInvMassBin_MCmatch[i]= new TH1F(HistName.Data(), HistName.Data(), 500, 1.819, 1.919 );
    fDplusInvMassBin_MCmatch[i]->Sumw2();
    HistName.Form("fB0InvMass_bin%d_PassCut1",i);
    fB0InvMassBin_PassCut1[i]= new TH1F(HistName.Data(), HistName.Data(), nbins, fBMassLowLimit, fBMassUpLimit);
    fB0InvMassBin_PassCut1[i]->Sumw2();
    HistName.Form("fB0InvMassMinusD_bin%d_MCmatch_PassCut1",i);
    fB0InvMassBinMinusD_MCmatch_PassCut1[i]= new TH1F(HistName.Data(), HistName.Data(), nbins, fBMassLowLimit-1.869, fBMassUpLimit-1.869);
    fB0InvMassBinMinusD_MCmatch_PassCut1[i]->Sumw2();
    HistName.Form("fDplusInvMass_bin%d_PassCut1",i);
    fDplusInvMassBin_PassCut1[i]= new TH1F(HistName.Data(), HistName.Data(), 500, 1.819, 1.919 );
   
    
   
    HistName.Form("fBDecayLengthBin_bin%d",i);
    fBDecayLengthBin[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0, 0.5);
    HistName.Form("fBNormalizedDecayLengthBin_bin%d",i);
    fBNormalizedDecayLengthBin[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0, 20);
    HistName.Form("fBPointingBin_bin%d",i);
    fBPointingBin[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0, 1.);
    HistName.Form("fBPointingXYBin_bin%d",i);
    fBPointingXYBin[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0, 1.);
    HistName.Form("fDplusd0Bin_bin%d",i);
    fDplusd0Bin[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0, 0.25);
    HistName.Form("fBpTBin_bin%d",i);
    fBpTBin[i]= new TH1F(HistName.Data(), HistName.Data(), 100, 0, 50);
    HistName.Form("fDpluspTBin_bin%d",i);
    fDpluspTBin[i]= new TH1F(HistName.Data(), HistName.Data(), 100, 0, 50);
    HistName.Form("fproductd0Bin_bin%d",i);
    fproductd0Bin[i]= new TH1F(HistName.Data(), HistName.Data(), 200, -0.001, 0.001);
    HistName.Form("fDplusPointingBin_bin%d",i);
    fDplusPointingBin[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0., 1.);
    HistName.Form("fDplusDecayLengthBin_bin%d",i);
    fDplusDecayLengthBin[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0, 0.5);
    HistName.Form("fDplusNormalizedDecayLengthBin_bin%d",i);
    fDplusNormalizedDecayLengthBin[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0, 20);
    HistName.Form("fDplusSigmaVertBin_bin%d",i);
    fDplusSigmaVertBin[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0, 0.5);
    
    HistName.Form("fBDecayLengthBin_MCmatch_bin%d",i);
    fBDecayLengthBin_MCmatch[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0, 0.5);
    HistName.Form("fBNormalizedDecayLengthBin_MCmatch_bin%d",i);
    fBNormalizedDecayLengthBin_MCmatch[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0, 20);
    HistName.Form("fBPointingBin_MCmatch_bin%d",i);
    fBPointingBin_MCmatch[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0, 1.);
    HistName.Form("fBPointingXYBin_MCmatch_bin%d",i);
    fBPointingXYBin_MCmatch[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0, 1.);
    HistName.Form("fDplusd0Bin_MCmatch_bin%d",i);
    fDplusd0Bin_MCmatch[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0, 0.25);
    HistName.Form("fBpTBin_MCmatch_bin%d",i);
    fBpTBin_MCmatch[i]= new TH1F(HistName.Data(), HistName.Data(), 100, 0, 50);
    HistName.Form("fDpluspTBin_MCmatch_bin%d",i);
    fDpluspTBin_MCmatch[i]= new TH1F(HistName.Data(), HistName.Data(), 100, 0, 50);
    HistName.Form("fproductd0Bin_MCmatch_bin%d",i);
    fproductd0Bin_MCmatch[i]= new TH1F(HistName.Data(), HistName.Data(), 200, -0.001, 0.001);
    HistName.Form("fDplusPointingBin_MCmatch_bin%d",i);
    fDplusPointingBin_MCmatch[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0., 1.);
    HistName.Form("fDplusDecayLengthBin_MCmatch_bin%d",i);
    fDplusDecayLengthBin_MCmatch[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0, 0.5);
    HistName.Form("fDplusNormalizedDecayLengthBin_MCmatch_bin%d",i);
    fDplusNormalizedDecayLengthBin_MCmatch[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0, 20);
    HistName.Form("fDplusSigmaVertBin_MCmatch_bin%d",i);
    fDplusSigmaVertBin_MCmatch[i]= new TH1F(HistName.Data(), HistName.Data(), 200, 0, 0.5);
   
    HistName.Form("fd0MMExpD_bin%d",i);
    fd0MMExpDBin[i] = new TH1F(HistName.Data(), HistName.Data(), 200, -10.0,10.);
    HistName.Form("fd0MMExpD_MCmatch_bin%d",i);
    fd0MMExpDBin_MCmatch[i] = new TH1F(HistName.Data(), HistName.Data(), 200, -10.0,10.);
    HistName.Form("fd0MMExppi_bin%d",i);
    fd0MMExppiBin[i] = new TH1F(HistName.Data(), HistName.Data(), 200, -10.0,10.);
    HistName.Form("fd0MMExppi_MCmatch_bin%d",i);
    fd0MMExppiBin_MCmatch[i] = new TH1F(HistName.Data(), HistName.Data(), 200, -10.0,10.);
    HistName.Form("fd0MMExpDDaughtersBin%d",i);
    fd0MMExpDDaughtersBin[i] = new TH1F(HistName.Data(), HistName.Data(), 200, -10.0,10.);
    HistName.Form(" fd0MMExpDDaughters_MCmatch_bin%d",i);
    fd0MMExpDDaughtersBin_MCmatch[i] = new TH1F(HistName.Data(), HistName.Data(), 200, -10.0,10.);
   
    HistName.Form("fBDplusPtBin_bin%d",i);
    fBDplusPtBin[i]= new TH2F(HistName.Data(), HistName.Data(), 100, 0, 50,100,0,50);
    HistName.Form("fBpiPtBin_bin%d",i);
    fBpiPtBin[i]= new TH2F(HistName.Data(), HistName.Data(), 100, 0, 50,100,0,50);
    HistName.Form("fBDplusPtBin_bin%d_MCmatch",i);
    fBDplusPtBin_MCmatch[i]= new TH2F(HistName.Data(), HistName.Data(), 100, 0, 50,100,0,50);
    HistName.Form("fBpiPtBin_bin%d_MCmatch",i);
    fBpiPtBin_MCmatch[i]= new TH2F(HistName.Data(), HistName.Data(), 100, 0, 50,100,0,50);
    

   

    if(kbins){
    fOutputBins->Add(fB0InvMassBin[i]);
    fOutputBins->Add(fB0InvMassBinMinusD_PassCut1[i]);
    fOutputBins->Add(fDplusInvMassBin[i]); 
    fOutputBins->Add(fB0InvMassBin_PassCut1[i]);
    fOutputBins->Add(fDplusInvMassBin_PassCut1[i]);
    
    
    fOutputBins->Add(fBDecayLengthBin[i]);
    fOutputBins->Add(fBNormalizedDecayLengthBin[i]);
    fOutputBins->Add(fBPointingBin[i]);
    fOutputBins->Add(fBPointingXYBin[i]);
    fOutputBins->Add(fDplusd0Bin[i]);
    fOutputBins->Add(fBpTBin[i]);
    fOutputBins->Add(fDpluspTBin[i]);
    fOutputBins->Add(fproductd0Bin[i]);
    fOutputBins->Add(fDplusPointingBin[i]);
    fOutputBins->Add(fDplusDecayLengthBin[i]);
    fOutputBins->Add(fDplusNormalizedDecayLengthBin[i]);
    fOutputBins->Add(fDplusSigmaVertBin[i]);

    fOutputBins->Add(fd0MMExpDBin[i]);
    fOutputBins->Add(fd0MMExppiBin[i]);
    fOutputBins->Add(fd0MMExpDDaughtersBin[i]);
    fOutputBins->Add(fBDplusPtBin[i]);
    fOutputBins->Add(fBpiPtBin[i]);
    if(fReadMC){
      fOutputBins->Add(fB0InvMassBin_MCmatch[i]);
      fOutputBins->Add(fB0InvMassBin_MCmatch_PassCut1[i]); 
      fOutputBins->Add(fB0InvMassBinMinusD_MCmatch_PassCut1[i]);
      fOutputBins->Add(fDplusInvMassBin_MCmatch[i]);
      fOutputBins->Add(fBDecayLengthBin_MCmatch[i]);
      fOutputBins->Add(fBNormalizedDecayLengthBin_MCmatch[i]);
      fOutputBins->Add(fBPointingBin_MCmatch[i]);
      fOutputBins->Add(fBPointingXYBin_MCmatch[i]);
      fOutputBins->Add(fDplusd0Bin_MCmatch[i]);
      fOutputBins->Add(fBpTBin_MCmatch[i]);
      fOutputBins->Add(fDpluspTBin_MCmatch[i]);
      fOutputBins->Add(fproductd0Bin_MCmatch[i]);
      fOutputBins->Add(fDplusPointingBin_MCmatch[i]);
      fOutputBins->Add(fDplusDecayLengthBin_MCmatch[i]);
      fOutputBins->Add(fDplusNormalizedDecayLengthBin_MCmatch[i]);
      fOutputBins->Add(fDplusSigmaVertBin_MCmatch[i]);
      fOutputBins->Add(fd0MMExpDBin_MCmatch[i]);
      fOutputBins->Add(fd0MMExppiBin_MCmatch[i]);
      fOutputBins->Add(fd0MMExpDDaughtersBin_MCmatch[i]);
      fOutputBins->Add(fBDplusPtBin_MCmatch[i]);
      fOutputBins->Add(fBpiPtBin_MCmatch[i]);


    
    }
    }
   }
   

  

  
  return;
}
//________________________________________________________________________
void AliAnalysisTaskSEB0toDminuspi::FillSpectrum(Int_t nB0, Int_t* arrayMClabel,Double_t Bfield){
  //
  // Fill histos for B0 spectrum
  //
  if (nB0 ==0) return;
  
  //Note: all this B0 candidates are "made" from the same Dplus, that's why fArrayB0[0]
  AliAODRecoDecayHF3Prong* Dmeson = (AliAODRecoDecayHF3Prong*)fArrayB0[0]->GetSecondaryVtx()->GetDaughter(0);	
  UInt_t prongs[2]={411,211};
  
  Double_t d0MM[3];Double_t dd0min;
  for(Int_t i=0; i<3; i++) {
    Double_t diffIP, errdiffIP;
   Dmeson->Getd0MeasMinusExpProng(i,Bfield,diffIP,errdiffIP);
    Double_t normdd0= diffIP/errdiffIP;
    d0MM[i]=normdd0;
    
    if(i==0) dd0min=normdd0;
    else if(TMath::Abs(normdd0)<TMath::Abs(dd0min)) dd0min=normdd0;
    
      }
  
  //if(TMath::Abs(dd0min)<fInvertedTopomaticCutOnDDaughters)continue;
      
   

  
  Double_t Dmass= Dmeson->InvMassDplus();
  Double_t B0InvMass, d01, d02, productd0;
  
  for (int iB0 =0; iB0<nB0; iB0++){
  

	
    B0InvMass = fArrayB0[iB0]->InvMass2Prongs(0 , 1 ,prongs[0], prongs[1]);
    AliAODTrack* aodTrack = (AliAODTrack*)fArrayB0[iB0]->GetSecondaryVtx()->GetDaughter(1);
    d01 = TMath::Abs(fArrayB0[iB0]->Getd0Prong(0));
    d02 = TMath::Abs(fArrayB0[iB0]->Getd0Prong(1));
    productd0 = fArrayB0[iB0]->Getd0Prong(0) * fArrayB0[iB0]->Getd0Prong(1);

    Double_t diffIPMCD, errdiffIPMCD;
    fArrayB0[iB0]->Getd0MeasMinusExpProng(0,Bfield,diffIPMCD,errdiffIPMCD);
    Double_t normdd0MCD= diffIPMCD/errdiffIPMCD;
    
    Double_t diffIPMCpi, errdiffIPMCpi;
    fArrayB0[iB0]->Getd0MeasMinusExpProng(1,Bfield,diffIPMCpi,errdiffIPMCpi);
    Double_t normdd0MCpi= diffIPMCpi/errdiffIPMCpi;
    
    Double_t cosoa=(fArrayB0[iB0]->PxProng(0)*fArrayB0[iB0]->PxProng(1)+fArrayB0[iB0]->PyProng(0)*fArrayB0[iB0]->PyProng(1)+fArrayB0[iB0]->PzProng(0)*fArrayB0[iB0]->PzProng(1))/(TMath::Sqrt(fArrayB0[iB0]->P2Prong(0))*TMath::Sqrt(fArrayB0[iB0]->P2Prong(1)));
    

     fB0InvMass->Fill(B0InvMass);
     if(arrayMClabel[iB0] > (-1)){
       fBInvMass_MCmatch->Fill(B0InvMass);}
	 
       if(kbins){
       Int_t index = fRDCutsAnalysis->PtBin(fArrayB0[iB0]->Pt());
       fB0InvMassBin[index]->Fill(B0InvMass);
       fDplusInvMassBin[index]->Fill(Dmeson->InvMassDplus());
       if(arrayMClabel[iB0] > (-1))   {
	 	 fB0InvMassBin_MCmatch[index]->Fill(B0InvMass);
       }
     }
     
       
     Bool_t PassCut1=0; 
    //Cuts array: 0= D pT,1= pi pT, 2= prod d0, 3= cos_point B, 4= cos_pointXY B, 5= decl B, 6=NdlXY B,7= Topomatic D, 8= topmatico pi, 9=cos opening angle 
     
     if(Dmeson->Pt()>fpTD &&  aodTrack->Pt()>fpTpi &&productd0 < fprodd0 && fArrayB0[iB0]->CosPointingAngle()>fcosB && fArrayB0[iB0]->CosPointingAngleXY()>fcosBXY &&
	 fArrayB0[iB0]->DecayLength()>fdlB && fArrayB0[iB0]->NormalizedDecayLengthXY()> fNdlBXY  && TMath::Abs(normdd0MCpi)<fTopomaticpi && TMath::Abs(normdd0MCD)<fTopomaticD && cosoa>fcosoaDpi  && TMath::Abs(dd0min)>fInvertedTopomaticCutOnDDaughters) PassCut1=1;
   
	
	
    if(PassCut1){

         fCEvents->Fill(17);
      fB0InvMass_PassCut1->Fill(B0InvMass);
      fB0InvMassMinusD_PassCut1->Fill(B0InvMass-Dmass);  	 
      fDplusInvMass_PassCut1->Fill(Dmeson->InvMassDplus());
      //histo for the D+ (ex. if two B are found with the same D it will be counted twice)
      fBDplusPt->Fill(fArrayB0[iB0]->Pt(),Dmeson->Pt());
      fBpiPt->Fill(fArrayB0[iB0]->Pt(),aodTrack->Pt());
      fDplusPt->Fill(Dmeson->Pt());
      fDplusEta->Fill(Dmeson->Eta());
      fDplusPhi->Fill(Dmeson->Phi());
      fDplusPointing->Fill(Dmeson->CosPointingAngle());
      fDplusDecayLength->Fill(Dmeson->DecayLength());
      fDplusNormalizedDecayLength->Fill(Dmeson->NormalizedDecayLengthXY());
      fDplusSigmaVert->Fill(Dmeson->GetSigmaVert());
      fDplusd0->Fill(d01);
   
      
      
      fpiPt->Fill(aodTrack->Pt());
      fpiEta->Fill(aodTrack->Eta());
      fpiPhi->Fill(aodTrack->Phi());
      fBPt->Fill(fArrayB0[iB0]->Pt());
      fBEta->Fill(fArrayB0[iB0]->Eta());
      fBPhi->Fill(fArrayB0[iB0]->Phi());
      fBDecayLength->Fill(fArrayB0[iB0]->DecayLength());
      fBNormalizedDecayLength->Fill(fArrayB0[iB0]->NormalizedDecayLengthXY());
      fBPointing->Fill(fArrayB0[iB0]->CosPointingAngle());
      fBPointingXY->Fill(fArrayB0[iB0]->CosPointingAngleXY());
      fpid0->Fill(d02);
      fproductd0->Fill(productd0);
      fd0MMExpD->Fill(normdd0MCD);
      fd0MMExppi->Fill(normdd0MCpi);
      fd0MMExpDDaughters->Fill(dd0min);

      fcosoa->Fill(cosoa); 
    
      
      
    if(arrayMClabel[iB0] > (-1)){
            	
      fCEvents->Fill(18);
      fBInvMass_MCmatch_PassCut1->Fill(B0InvMass); 
      fB0InvMassMinusD_MCmatch_PassCut1->Fill(B0InvMass-Dmass); 
      fBDplusPt_MCmatch->Fill(fArrayB0[iB0]->Pt(),Dmeson->Pt());
      fBpiPt_MCmatch->Fill(fArrayB0[iB0]->Pt(),aodTrack->Pt());
     
      fpiPt_MCmatch->Fill(aodTrack->Pt());
      fpiEta_MCmatch->Fill(aodTrack->Eta());
      fpiPhi_MCmatch->Fill(aodTrack->Phi());
      fBPt_MCmatch->Fill(fArrayB0[iB0]->Pt());
      fBEta_MCmatch->Fill(fArrayB0[iB0]->Eta());
      fBPhi_MCmatch->Fill(fArrayB0[iB0]->Phi());
      fBDecayLength_MCmatch->Fill(fArrayB0[iB0]->DecayLength());
      fBNormalizedDecayLength_MCmatch->Fill(fArrayB0[iB0]->NormalizedDecayLengthXY());
      fBPointing_MCmatch->Fill(fArrayB0[iB0]->CosPointingAngle());
      fBPointingXY_MCmatch->Fill(fArrayB0[iB0]->CosPointingAngleXY());
      fDplusd0_MCmatch->Fill(d01);
      fpid0_MCmatch->Fill(d02);
      fproductd0_MCmatch->Fill(productd0);
      fDplusPt_MCmatch->Fill(Dmeson->Pt());
      fDplusEta_MCmatch->Fill(Dmeson->Eta());
      fDplusPhi_MCmatch->Fill(Dmeson->Phi());
      fDplusPointing_MCmatch->Fill(Dmeson->CosPointingAngle());
      fDplusDecayLength_MCmatch->Fill(Dmeson->DecayLength());
      fDplusNormalizedDecayLength_MCmatch->Fill(Dmeson->NormalizedDecayLengthXY());
      fDplusSigmaVert_MCmatch->Fill(Dmeson->GetSigmaVert());
      fDplusInvMass_MCmatch->Fill(Dmeson->InvMassDplus());
      fd0MMExpDMCmatch->Fill(normdd0MCD);
      fd0MMExppiMCmatch->Fill(normdd0MCpi);
      fd0MMExpDDaughters_MCmatch->Fill(dd0min);
      fcosoa_MCmatch->Fill(cosoa); 
   
    }
    
 
    
    if(kbins){
      Int_t index = fRDCutsAnalysis->PtBin(fArrayB0[iB0]->Pt());

      fB0InvMassBin_PassCut1[index]->Fill(B0InvMass);
      fB0InvMassBinMinusD_PassCut1[index]->Fill(B0InvMass-Dmass);
      fDplusInvMassBin_PassCut1[index]->Fill(Dmeson->InvMassDplus());
      fBDecayLengthBin[index]->Fill(fArrayB0[iB0]->DecayLength());
      fBNormalizedDecayLengthBin[index]->Fill(fArrayB0[iB0]->NormalizedDecayLengthXY());
      fBPointingBin[index]->Fill(fArrayB0[iB0]->CosPointingAngle());
      fBPointingXYBin[index]->Fill(fArrayB0[iB0]->CosPointingAngleXY());
      fDplusd0Bin[index]->Fill(d01);
      fBpTBin[index]->Fill(fArrayB0[iB0]->Pt());
      fDpluspTBin[index]->Fill(Dmeson->Pt());
      fproductd0Bin[index]->Fill(productd0);
      fDplusPointingBin[index]->Fill(Dmeson->CosPointingAngle());
      fDplusDecayLengthBin[index]->Fill(Dmeson->DecayLength());
      fDplusNormalizedDecayLengthBin[index]->Fill(Dmeson->NormalizedDecayLengthXY());    
      fDplusSigmaVertBin[index]->Fill(Dmeson->GetSigmaVert());
      fd0MMExpDBin[index]->Fill(normdd0MCD);
      fd0MMExppiBin[index]->Fill(normdd0MCpi);
      fd0MMExpDDaughtersBin[index]->Fill(dd0min);

      
      fBDplusPtBin[index]->Fill(fArrayB0[iB0]->Pt(),Dmeson->Pt());
      fBpiPtBin[index]->Fill(fArrayB0[iB0]->Pt(),aodTrack->Pt());
	   
      if(arrayMClabel[iB0] > (-1)){
	
	fB0InvMassBinMinusD_MCmatch_PassCut1[index]->Fill(B0InvMass-Dmass);
	fB0InvMassBin_MCmatch_PassCut1[index]->Fill(B0InvMass);
	fDplusInvMassBin_MCmatch[index]->Fill(Dmeson->InvMassDplus());
	fBDecayLengthBin_MCmatch[index]->Fill(fArrayB0[iB0]->DecayLength());
	fBNormalizedDecayLengthBin_MCmatch[index]->Fill(fArrayB0[iB0]->NormalizedDecayLengthXY());
	fBPointingBin_MCmatch[index]->Fill(fArrayB0[iB0]->CosPointingAngle());
	fBPointingXYBin_MCmatch[index]->Fill(fArrayB0[iB0]->CosPointingAngleXY());
	fDplusd0Bin_MCmatch[index]->Fill(d01);
	fBpTBin_MCmatch[index]->Fill(fArrayB0[iB0]->Pt());
	fDpluspTBin_MCmatch[index]->Fill(Dmeson->Pt());
	fproductd0Bin_MCmatch[index]->Fill(productd0);
	fDplusPointingBin_MCmatch[index]->Fill(Dmeson->CosPointingAngle());
	fDplusDecayLengthBin_MCmatch[index]->Fill(Dmeson->DecayLength());
	fDplusNormalizedDecayLengthBin_MCmatch[index]->Fill(Dmeson->NormalizedDecayLengthXY());
	fDplusSigmaVertBin_MCmatch[index]->Fill(Dmeson->GetSigmaVert());
    	fd0MMExpDBin_MCmatch[index]->Fill(normdd0MCD);
	fd0MMExppiBin_MCmatch[index]->Fill(normdd0MCpi);
	fd0MMExpDDaughtersBin_MCmatch[index]->Fill(dd0min);
	fBDplusPtBin_MCmatch[index]->Fill(fArrayB0[iB0]->Pt(),Dmeson->Pt());
	fBpiPtBin_MCmatch[index]->Fill(fArrayB0[iB0]->Pt(),aodTrack->Pt());
	    
	    
      }   
    }
   }	
  }
  

  return;
}


//___________________________________________________________________________-

AliAODVertex* AliAnalysisTaskSEB0toDminuspi::RecalculateVertex(const AliVVertex *primary,TObjArray *tracks,Double_t bField) {
  //
  // Helper function to recalculate a vertex.
  //

  AliESDVertex *vertexESD = 0;
  AliAODVertex *vertexAOD = 0;
  //Double_t covmatrix[6];
  // AliVertexerTracks
  AliVertexerTracks *vertexer = new AliVertexerTracks(bField);
  vertexer->SetVtxStart((AliESDVertex*)primary);//primary vertex
  vertexESD = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(tracks);
  delete vertexer; vertexer=NULL;

  if(!vertexESD) return vertexAOD;


  if(vertexESD->GetNContributors()!=tracks->GetEntriesFast()) {

    delete vertexESD; vertexESD=NULL;
    return vertexAOD;
  }

  // Double_t pos[3],cov[6],chi2perNDF;
  // vertexESD->GetXYZ(pos);
  // convert to AliAODVertex
  //
  Double_t dispersion;
  Double_t pos[3],cov[6],chi2perNDF;
  for(Int_t a=0;a<3;a++)pos[a]=0.;
  for(Int_t b=0;b<6;b++)cov[b]=0.;
  chi2perNDF=0;
  //
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix


  Double_t vertRadius2=pos[0]*pos[0]+pos[1]*pos[1];
  if(vertRadius2>8.){//(2.82)^2 radius beam pipe
    //  vertex outside beam pipe, reject candidate to avoid propagation through material
    delete vertexESD; vertexESD=NULL;
    return vertexAOD;
  }
  
  chi2perNDF = vertexESD->GetChi2toNDF();
  dispersion = vertexESD->GetDispersion();
  delete vertexESD; vertexESD=NULL;
  Int_t nprongs= tracks->GetEntriesFast();
  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,nprongs);
  return vertexAOD;

}



//----------------------------------------------------------------------------
