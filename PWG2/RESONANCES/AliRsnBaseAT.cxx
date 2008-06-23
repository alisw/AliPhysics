#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "AliLog.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliRsnEvent.h"


#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODInputHandler.h"

#include "AliRsnBaseAT.h"


ClassImp ( AliRsnBaseAT );

//________________________________________________________________________
AliRsnBaseAT::AliRsnBaseAT ( const char *name /*,Bool_t isSecondChain*/ )
    : AliAnalysisTask ( name, "" ) /*,fIsSecondChain ( isSecondChain )*/
{
  InitIOVars ();
  DefineInput ( 0, TChain::Class() );
//   if ( IsSecondChain() )
//     DefineInput ( 1, TChain::Class() );
}

void AliRsnBaseAT::InitIOVars ()
{
  AliDebug ( AliLog::kDebug, "<-" );

  fNumOfEvents=0;
  for ( Int_t i=0;i<3;i++ )
  {
    fChain[i]=0;
    fRSN[i] = 0;
    fESD[i] = 0;
    fAOD[i] = 0;
    fInputType[i] = kRSN;
  }
  
  fAnalysisMgr=0;
  
  AliDebug ( AliLog::kDebug, "->" );
}

Bool_t AliRsnBaseAT::Notify()
{
  return AliAnalysisTask::Notify();
}

//________________________________________________________________________
void AliRsnBaseAT::ConnectInputData ( Option_t * )
{
  ConnectInputDataByInputType ( fInputType[0],0 );
//   if ( IsSecondChain() )
//     ConnectInputDataByInputType ( fInputType[1],1 );

}

void AliRsnBaseAT::ConnectInputDataByInputType ( EInputType type ,Short_t inputIndex )
{
  AliDebug ( AliLog::kDebug, "<-" );

  switch ( type )
  {
    case kAOD:
    {
      ConnectAOD ( inputIndex );
      break;
    }
    case kESD:
    {
      ConnectESD ( inputIndex );
      break;
    }
    case kESDMC:
    {
      ConnectESDMC ( inputIndex );
      break;
    }
    case kMC:
      AliError ( "Not Implemented Yet ..." );
      break;
    case kRSN:
    {
      ConnectRSN ( inputIndex );
      break;
    }
    default:
      AliError ( "Type not supported ..." );
      break;
  }
  AliDebug ( AliLog::kDebug, "->" );
}

void AliRsnBaseAT::ConnectRSN ( Short_t inputIndex )
{
  AliDebug ( AliLog::kDebug, "<-" );
  char ** address = ( char ** ) GetBranchAddress ( inputIndex, "RsnEvents" );
  if ( address )
  {
    fRSN[inputIndex] = ( AliRsnEvent* ) ( *address );
  }
  else
  {
//     fRSN[inputIndex] = new AliRsnEvent();
    fRSN[inputIndex] = 0;
    SetBranchAddress ( inputIndex, "RsnEvents", &fRSN[inputIndex] );
  }
  AliDebug ( AliLog::kDebug, "->" );
}

void AliRsnBaseAT::ConnectESD(Short_t inputIndex)
{
  AliDebug ( AliLog::kDebug, "<-" );

//   fAnalysisMgr->SetInputEventHandler ( new AliESDInputHandler() );

  TTree* tree = dynamic_cast<TTree*> ( GetInputData ( inputIndex ) );
  if ( !tree ) { AliError ( "Could not read chain from input slot 0" ); }
  else
  {
    // Disable all branches, we want to process only MC
    tree->SetBranchStatus ( "*", kFALSE );
    tree->SetBranchStatus ( "fTracks.*", kTRUE );

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> ( AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler() );

    if ( !esdH ) { AliError ( "Could not get ESDInputHandler" ); }
    else
      fESD[inputIndex] = esdH->GetEvent();
  }
  AliDebug ( AliLog::kDebug, "->" );

}

void AliRsnBaseAT::ConnectESDMC(Short_t inputIndex)
{
  AliDebug ( AliLog::kDebug, "<-" );

//   fAnalysisMgr->SetInputEventHandler ( new AliESDInputHandler() );
//   fAnalysisMgr->SetMCtruthEventHandler ( new AliMCEventHandler() );


  TTree* tree = dynamic_cast<TTree*> ( GetInputData ( inputIndex ) );
  if ( !tree ) { AliError ( "Could not read chain from input slot 0" ); }
  else
  {
    // Disable all branches, we want to process only MC
    tree->SetBranchStatus ( "*", kFALSE );
    tree->SetBranchStatus ( "fTracks.*", kTRUE );

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> ( AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler() );

    if ( !esdH ) { AliError ( "Could not get ESDInputHandler" ); }
    else
      fESD[inputIndex] = esdH->GetEvent();
  }
  AliDebug ( AliLog::kDebug, "->" );

}

void AliRsnBaseAT::ConnectAOD(Short_t inputIndex)
{
  AliDebug ( AliLog::kDebug, "<-" );
  
  //   fAnalysisMgr->SetInputEventHandler ( new AliAODInputHandler());
  
  TTree* tree = dynamic_cast<TTree*> ( GetInputData ( inputIndex ) );
  if ( !tree ) { AliError ( "Could not read chain from input slot 0" );}
  else
  {
    AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> ( AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler() );

    if ( !aodH ) { AliError ( "Could not get AODInputHandler" ); }
    else
    {
      fAOD[inputIndex] = aodH->GetEvent();
    }
  }
  AliDebug ( AliLog::kDebug, "->" );
}
