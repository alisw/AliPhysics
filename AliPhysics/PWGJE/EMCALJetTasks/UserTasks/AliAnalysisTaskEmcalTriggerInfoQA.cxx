#include "AliEmcalTriggerSetupInfo.h"
#include "AliAnalysisTaskEmcalTriggerInfoQA.h"

#include <Riostream.h>
#include <ctime>
#include <TString.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TRandom.h>
#include <TRandom3.h>

#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliEmcalJet.h"
#include "AliEMCALGeometry.h"

#include "AliEMCALTriggerPatchInfo.h"

ClassImp(AliAnalysisTaskEmcalTriggerInfoQA)

//________________________________________________________________________
AliAnalysisTaskEmcalTriggerInfoQA::AliAnalysisTaskEmcalTriggerInfoQA() : 
    AliAnalysisTaskEmcal(),
    fOutput(0),
    fHistos(0),
    fTriggersInfo(0),
    fTriggerSetup(0),
    fIsInitialized(kFALSE),
    fCaloTriggerPatchInfoName("EmcalTriggers"),
    fCaloTriggerSetupInfoName("EmcalTriggerSetup")
{
}

//________________________________________________________________________
AliAnalysisTaskEmcalTriggerInfoQA::AliAnalysisTaskEmcalTriggerInfoQA(const char *name) :
    AliAnalysisTaskEmcal(name),
    fOutput(0),
    fHistos(0),
    fTriggersInfo(0),
    fTriggerSetup(0),
    fIsInitialized(kFALSE),
    fCaloTriggerPatchInfoName("EmcalTriggers"),
    fCaloTriggerSetupInfoName("EmcalTriggerSetup")
{
    DefineOutput(1,TList::Class());  // for output list
}

//________________________________________________________________________
AliAnalysisTaskEmcalTriggerInfoQA::~AliAnalysisTaskEmcalTriggerInfoQA()
{
//   // Destructor. Clean-up the output list, but not the histograms that are put inside
//   // (the list is owner and will clean-up these histograms). Protect in PROOF case.
//     if (fHistos && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())
//     {
//         delete fOutput;
//     }
//     delete fTrackCuts;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalTriggerInfoQA::UserCreateOutputObjects()
{
  Int_t i;
  
  // Create histograms
  // Called once (on the worker node)
  fIsInitialized=kFALSE;
  fOutput = new TList();
  fOutput->SetOwner();  // IMPORTANT!

  char buf[100];
  char strId[3][10];
  sprintf( strId[0], "JetLow" );
  sprintf( strId[1], "JetHigh" );
  sprintf( strId[2], "JetMain" );

  fHistos = new TH1*[5000];

  for( i = 0; i < 100; i++ )
    fHistos[i] = 0;
  
  sprintf( buf, "fhPatchNTotal" );
  fHistos[0] = new TH1D( buf, buf, 200, -0.5, 199.5 );
  fHistos[0]->SetXTitle( "N_{patch} in event" );
  fHistos[0]->SetYTitle( "N" );

  // common histos
  for( i = 0; i < 3; i++ ){
    
    sprintf( buf, "fhPatchN%s", strId[i] );
    fHistos[1000*(i+1)+0] = new TH1D( buf, buf, 200, -0.5, 199.5 );
    fHistos[1000*(i+1)+0]->SetXTitle( "N_{patch} in event" );
    fHistos[1000*(i+1)+0]->SetYTitle( "N" );

    sprintf( buf, "fhThreshold%s", strId[i] );
    fHistos[1000*(i+1)+1] = new TH1D( buf, buf, 500, -0.5, 2999.5 );
    fHistos[1000*(i+1)+1]->SetXTitle( "thresholds [ADC]" );
    fHistos[1000*(i+1)+1]->SetYTitle( "N" );

    sprintf( buf, "fhThresholdGeV%s", strId[i] );
    fHistos[1000*(i+1)+2] = new TH1D( buf, buf, 500, 0, 250 );
    fHistos[1000*(i+1)+2]->SetXTitle( "thresholds [GeV]" );
    fHistos[1000*(i+1)+2]->SetYTitle( "N" );
  
    sprintf( buf, "fhADCPatch%s", strId[i] );
    fHistos[1000*(i+1)+3] = new TH1D( buf, buf, 500, -0.5, 2999.5 );
    fHistos[1000*(i+1)+3]->SetXTitle( "E [ADC]" );
    fHistos[1000*(i+1)+3]->SetYTitle( "N" );

    sprintf( buf, "fhADCPatchGeV%s", strId[i] );
    fHistos[1000*(i+1)+4] = new TH1D( buf, buf, 500, 0, 250 );
    fHistos[1000*(i+1)+4]->SetXTitle( "E [GeV]" );
    fHistos[1000*(i+1)+4]->SetYTitle( "N" );
  
    sprintf( buf, "fhECells%s", strId[i] );
    fHistos[1000*(i+1)+5] = new TH1D( buf, buf, 500, 0, 250 );
    fHistos[1000*(i+1)+5]->SetXTitle( "E cells [GeV]" );
    fHistos[1000*(i+1)+5]->SetYTitle( "N" );
  
    sprintf( buf, "fhPatchSpacial%s", strId[i] );
    fHistos[1000*(i+1)+6] = new TH2D( buf, buf, 9, -0.5, 8.5,
                                                13, -0.5, 12.5 );
    fHistos[1000*(i+1)+6]->SetXTitle( "cells in -eta [8 cells]" );
    fHistos[1000*(i+1)+6]->SetYTitle( "cells in phi [8 cells]" );
    fHistos[1000*(i+1)+6]->SetOption( "COLZ" );
    fHistos[1000*(i+1)+6]->SetStats( kFALSE );
    
    sprintf( buf, "fhPatchSpacialADC%s", strId[i] );
    fHistos[1000*(i+1)+7] = new TH2D( buf, buf, 9, -0.5, 8.5,
                                                13, -0.5, 12.5 );
    fHistos[1000*(i+1)+7]->SetXTitle( "cells in -eta [8 cells]" );
    fHistos[1000*(i+1)+7]->SetYTitle( "cells in phi [8 cells]" );
    fHistos[1000*(i+1)+7]->SetOption( "COLZ" );
    fHistos[1000*(i+1)+7]->SetStats( kFALSE );

    sprintf( buf, "fhPatchSpacialECells%s", strId[i] );
    fHistos[1000*(i+1)+8] = new TH2D( buf, buf, 9, -0.5, 8.5,
                                                13, -0.5, 12.5 );
    fHistos[1000*(i+1)+8]->SetXTitle( "cells in -eta [8 cells]" );
    fHistos[1000*(i+1)+8]->SetYTitle( "cells in phi [8 cells]" );
    fHistos[1000*(i+1)+8]->SetOption( "COLZ" );
    fHistos[1000*(i+1)+8]->SetStats( kFALSE );

    sprintf( buf, "fhPatchCenter%s", strId[i] );
    fHistos[1000*(i+1)+9] = new TH2D( buf, buf, 100, -0.7, 0.7,
                                                100, 0, TMath::Pi()*2 );
    fHistos[1000*(i+1)+9]->SetXTitle( "#eta" );
    fHistos[1000*(i+1)+9]->SetYTitle( "#phi [rad]" );
    fHistos[1000*(i+1)+9]->SetOption( "COLZ" );
    fHistos[1000*(i+1)+9]->SetStats( kFALSE );
    
    sprintf( buf, "fhPatchCenterMass%s", strId[i] );
    fHistos[1000*(i+1)+10] = new TH2D( buf, buf, 100, -0.7, 0.7,
                                                100, 0, TMath::Pi()*2 );
    fHistos[1000*(i+1)+10]->SetXTitle( "#eta" );
    fHistos[1000*(i+1)+10]->SetYTitle( "#phi [rad]" );
    fHistos[1000*(i+1)+10]->SetOption( "COLZ" );
    fHistos[1000*(i+1)+10]->SetStats( kFALSE );
  
    sprintf( buf, "fhEInPatchSpacial%s", strId[i] );
    fHistos[1000*(i+1)+11] = new TH2D( buf, buf, 32, -0.5, 31.5,
                                                32, -0.5, 31.5 );
    fHistos[1000*(i+1)+11]->SetXTitle( "cells in -eta" );
    fHistos[1000*(i+1)+11]->SetYTitle( "cells in phi" );
    fHistos[1000*(i+1)+11]->SetOption( "COLZ" );
    fHistos[1000*(i+1)+11]->SetStats( kFALSE );

    sprintf( buf, "fhADCvsEcells%s", strId[i] );
    fHistos[1000*(i+1)+12] = new TH2D( buf, buf,  500, -0.5, 2999.5,
                                                500, 0, 250 );
    fHistos[1000*(i+1)+12]->SetXTitle( "E [ADC]" );
    fHistos[1000*(i+1)+12]->SetYTitle( "E cells [GeV]" );
    fHistos[1000*(i+1)+12]->SetOption( "COLZ" );
    fHistos[1000*(i+1)+12]->SetStats( kFALSE );

  }

  for( i = 0; i < 5000; i++ )
    if( fHistos[i] != 0 )
      fOutput->Add( fHistos[i] );
    
 
  // Post data for ALL output slots >0 here,
  // To get at least an empty histogram
  // 1 is the outputnumber of a certain weg of task 1
  PostData(1, fOutput);
}

void AliAnalysisTaskEmcalTriggerInfoQA::UserExecOnce()
{
  AliAnalysisTaskEmcal::ExecOnce();

  // Get the event tracks from PicoTracks
  fTriggersInfo =dynamic_cast <TClonesArray*>(InputEvent()->FindListObject( fCaloTriggerPatchInfoName.Data() ));
  
  fTriggerSetup = dynamic_cast <AliEmcalTriggerSetupInfo*>(InputEvent()->FindListObject( fCaloTriggerSetupInfoName.Data() ));
    
  fIsInitialized=kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskEmcalTriggerInfoQA::UserExec(Option_t *) 
{
  Int_t nPatch, iPatch, iMainPatch, nJetLow, nJetHigh, type;
  AliEMCALTriggerPatchInfo *patch;
  
  if (fIsInitialized==kFALSE)
    UserExecOnce();

    
  if( !fTriggersInfo ){
    AliError(Form("Calo triggers info container %s not available.", fCaloTriggerPatchInfoName.Data()));
    return;
  }
  if( !fTriggerSetup ){
    AliError(Form("Calo trigger setup info container %s not available.", fCaloTriggerSetupInfoName.Data()));
    return;
  }
  if( !fCaloCells ){
    AliError(Form("Calo cells container %s not available.", fCaloCellsName.Data()));
    return;
  }
  
  // first get the thresholds
  fHistos[1000+1]->Fill( fTriggerSetup->GetThresholdJetLow() );
  fHistos[2000+1]->Fill( fTriggerSetup->GetThresholdJetHigh() );
  fHistos[1000+2]->Fill( fTriggerSetup->GetThresholdGeVRoughJetLow() );
  fHistos[2000+2]->Fill( fTriggerSetup->GetThresholdGeVRoughJetHigh() );
  
  // now go through patchs
  nPatch = fTriggersInfo->GetEntries();
  
  fHistos[0]->Fill( nPatch );
  
  nJetLow = 0;
  nJetHigh = 0;
  iMainPatch = -1;
  
  for( iPatch = 0; iPatch < nPatch; iPatch++ ){
    
    patch = (AliEMCALTriggerPatchInfo*)fTriggersInfo->At( iPatch );
    
    // check if high/low threshold
    // high overrides the low, to avoid double counting
    type = -1;
    if( patch->IsJetHigh() ){
      type = 1;
      nJetHigh++;
    }
    else if( patch->IsJetLow() ){
      type = 0;
      nJetLow++;
    }
    
    if( type == -1 )
      continue;
    
    FillPatch( patch, type );
    
    // save main pach position
    if( patch->IsMainTrigger() )
      iMainPatch = iPatch;

  } // patches
  
  // fill jet counts per event
  fHistos[1000]->Fill( nJetLow );
  fHistos[2000]->Fill( nJetHigh );
  
  //  main trigger patch check -----------------------------
  if( iMainPatch == -1 ){
      // count of patches
      fHistos[3000]->Fill( 0 );
  }
  else{
    patch = (AliEMCALTriggerPatchInfo*)fTriggersInfo->At( iMainPatch );
    
    fHistos[3000]->Fill( 1 );
    
    // what threshold was the main trigger taken with
    if( patch->IsJetHigh() ){
      fHistos[3000+1]->Fill( fTriggerSetup->GetThresholdJetHigh() );
      fHistos[3000+2]->Fill( fTriggerSetup->GetThresholdGeVRoughJetHigh() );
    }
    else{
      fHistos[3000+1]->Fill( fTriggerSetup->GetThresholdJetLow() );
      fHistos[3000+2]->Fill( fTriggerSetup->GetThresholdGeVRoughJetLow() );
    }
    
    FillPatch( patch, 2 );
  }

  PostData(1, fOutput);
}


//________________________________________________________________________
void AliAnalysisTaskEmcalTriggerInfoQA::FillPatch( AliEMCALTriggerPatchInfo *patch, Int_t type ){
  
  // fills the patch parameters
  
  Int_t globCol, globRow, i, j, k, absId, cellAbsId[4];;
  
  // patch energies
  fHistos[1000*(type+1)+3]->Fill( patch->GetADCAmp() );
  fHistos[1000*(type+1)+4]->Fill( patch->GetADCAmpGeVRough() );
  fHistos[1000*(type+1)+5]->Fill( patch->GetPatchE() );
  
  //cout << "amp: " << patch->GetADCAmp() << endl;
  
  // get corner, convert from cells to trigger channels
  globCol = patch->GetEdgeCellX() / 2;
  globRow = patch->GetEdgeCellY() / 2;
  
  // fill in patch steps (8 cells = 4 trigger channels)
  ((TH2D*)fHistos[1000*(type+1)+6])->Fill( globCol/4, globRow/4 );
  ((TH2D*)fHistos[1000*(type+1)+7])->Fill( globCol/4, globRow/4, patch->GetADCAmp() );
  ((TH2D*)fHistos[1000*(type+1)+8])->Fill( globCol/4, globRow/4, patch->GetPatchE() );
  
  ((TH2D*)fHistos[1000*(type+1)+12])->Fill( patch->GetADCAmp(), patch->GetPatchE() );

  // phi/eta
  ((TH2D*)fHistos[1000*(type+1)+9])->Fill( patch->GetEtaGeo(), patch->GetPhiGeo() );
  ((TH2D*)fHistos[1000*(type+1)+10])->Fill( patch->GetEtaCM(), patch->GetPhiCM() );

  // E inside patch distribution
  // get the absolute trigger ID
  fGeom->GetAbsFastORIndexFromPositionInEMCAL( globCol, globRow, absId );
  // convert to the 4 absId of the cells composing the trigger channel
  fGeom->GetCellIndexFromFastORIndex( absId, cellAbsId );
  
  // sum the available energy in the 32/32 window of cells
  // step over trigger channels and get all the corresponding cells
  for( i = 0; i < 16; i++ ){
    for( j = 0; j < 16; j++ ){
      // get the 4 cells composing the trigger channel
      fGeom->GetAbsFastORIndexFromPositionInEMCAL( globCol+i, globRow+j, absId );
      fGeom->GetCellIndexFromFastORIndex( absId, cellAbsId );
      // add amplitudes and find patch edges
      for( k = 0; k < 4; k++ ){
        ((TH2D*)fHistos[1000*(type+1)+11])->Fill( i*2+k%2, j*2+k/2,
                          fCaloCells->GetCellAmplitude( cellAbsId[k] ));
      }
    }
  } // 32x32 cell window
 
}

//________________________________________________________________________
void AliAnalysisTaskEmcalTriggerInfoQA::Terminate(Option_t *) //specify what you want to have done
{
    // Called once at the end of the query. Done nothing here.
}
