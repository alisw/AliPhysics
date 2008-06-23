#include <TSystem.h>
#include <TFile.h>

#include "AliLog.h"

#include "AliAnalysisManager.h"
#include "AliRsnPairMgr.h"
#include "AliRsnEventBuffer.h"

#include "AliMCEventHandler.h"

#include "AliRsnAnalysisAT.h"

ClassImp ( AliRsnAnalysisAT )
AliRsnAnalysisAT::AliRsnAnalysisAT ( const char * name )
    : AliRsnBaseAT ( name )
{
  InitIOVars ();
  DefineInput ( 1, AliRsnPairMgr::Class() );

  DefineOutput ( 0, TList::Class() );
}

AliRsnAnalysisAT::~AliRsnAnalysisAT()
{
}

void AliRsnAnalysisAT::InitIOVars()
{
  AliDebug ( AliLog::kDebug, "<-" );
  AliRsnBaseAT::InitIOVars();

  fPairMgr = 0;
  fRsnMVEventBuffer = 0;
  fOutList = 0;

  for ( Int_t i=0;i<100 ;i++ )
    for ( Int_t j=0;j<100 ;j++ )
    {
      fHist[i][j] = 0;
    }

  AliDebug ( AliLog::kDebug, "->" );
}

void AliRsnAnalysisAT::LocalInit()
{
  AliDebug ( AliLog::kDebug, "<-" );
  fPairMgr = dynamic_cast<AliRsnPairMgr*> ( GetInputData ( 1 ) );
  fPairMgr->PrintPairs();
  AliDebug ( AliLog::kDebug, "->" );
}

Bool_t AliRsnAnalysisAT::Notify()
{
  AliDebug ( AliLog::kDebug, "<-" );
  fChain[0] = ( TChain* ) GetInputData ( 0 );
  if ( fChain[0] )
  {
    TFile *f = fChain[0]->GetCurrentFile();
    if ( f ) { AliInfo ( Form ( "Processing file %s",  f->GetName() ) );}
    else AliError ( "fTree->GetCurrentFile() is 0" );
    AliInfo ( Form ( "NumOfEvents %d",  fChain[0]->GetTree()->GetEntries() ) );
  }
  else
  {
    AliError ( "fChain[0] not available" );
  }

  AliDebug ( AliLog::kDebug, "->" );
  return   AliRsnBaseAT::Notify();
}

void AliRsnAnalysisAT::CreateOutputObjects()
{
  AliDebug ( AliLog::kDebug, "<-" );

  fPairMgr = dynamic_cast<AliRsnPairMgr*> ( GetInputData ( 1 ) );
  OpenFile ( 0 );
  fOutList = new TList();

  AliRsnPair *def=0;
  for ( Int_t i=0;i< fPairMgr->GetPairs()->GetEntriesFast();i++ )
  {
    def = ( AliRsnPair * ) fPairMgr->GetPairs()->At ( i );
    for ( Int_t j=0;j<def->GetCutMgr()->GetEntriesFast() ;j++ )
    {
      fHist[i][j] = def->GenerateEffMassHist ( j );
      fOutList->Add ( fHist[i][j] );
    }
  }

  fRsnMVEventBuffer = new AliRsnEventBuffer ( 1000 );
//   fRsnMVEventBuffer = new AliRsnEventBuffer ( 100 ,kFALSE );
  AliDebug ( AliLog::kDebug, "->" );

}

void AliRsnAnalysisAT::Exec ( Option_t * option )
{
  TTree *tree = ( ( TChain* ) GetInputData ( 0 ) )->GetTree();
  Long64_t ientry = ( Long64_t ) tree->GetReadEntry();

  if ( ientry%100==0 )
    AliInfo ( Form ( "Event #%d",ientry ) );

//   AliRsnEvent *curEvent = GetRsnMVEventFromInputType();
//   if ( !curEvent ) { AliError ( "Could not get AliRsnEvent from GetRsnMVEventFromInputType(). Skipping..." ); return; }

//   ProcessEventAnalysis ( curEvent );
//   PostEventProcess();

//   if (ientry%10000==0)
//   AliInfo(Form("Event #%d",ientry));
  PostData ( 0, fOutList );
}

void AliRsnAnalysisAT::Terminate ( Option_t * )
{
  AliDebug ( AliLog::kDebug, "<-" );
  fOutList = dynamic_cast<TList*> ( GetOutputData ( 0 ) );
  if ( !fOutList ) { AliError ( " fOutList not available" ); return; }
  fOutList->Print();

  AliDebug ( AliLog::kDebug, "->" );
}

void AliRsnAnalysisAT::Cleanup()
{
  AliInfo ( Form ( "Cleaning up in worker %s ...",gSystem->HostName() ) );
//   fRsnMVEventBuffer->ClearBuffer();
//   AliRsnPair *def = ( AliRsnPair * ) fPairMgr->GetPairs()->At ( 0 );
//   def->DoCleanUpAfterOneEvent();
}

void AliRsnAnalysisAT::ProcessEventAnalysis ( AliRsnEvent *curEvent )
{



  fRsnMVEventBuffer->AddEvent ( curEvent );
  AliRsnPair *def=0;
  AliRsnEvent *event=0;
  Int_t numOfTracks;
  for ( Int_t i=0;i< fPairMgr->GetPairs()->GetEntriesFast();i++ )
  {
    def = ( AliRsnPair * ) fPairMgr->GetPairs()->At ( i );
    def->SetRsnMVEventBuffer ( fRsnMVEventBuffer );
    for ( Int_t j=0;j<def->GetCutMgr()->GetEntriesFast() ;j++ )
    {
      event = fRsnMVEventBuffer->GetCurrentEvent();
      numOfTracks = event->GetMultiplicity();
//       AliInfo ( Form ( "%d",event->GetMultiplicity() ) );

      if ( numOfTracks>0 )
        def->ProcessPair ( event ,fHist[i][j] ,j );
    }
  }

}



AliRsnEvent * AliRsnAnalysisAT::GetRsnMVEventFromInputType ( const Short_t & index )
{
  switch ( fInputType[index] )
  {
    case kAOD:
    {
      return GetRsnMVFromAOD ( index );
      break;
    }
    case kESD:
    {
      AliWarning ( "Not Implemented Yet ..." );
      return GetRsnMVFromESD ( index );
      break;
    }
    case kESDMC:
    {
      AliWarning ( "Not Implemented Yet ..." );
      return GetRsnMVFromESDMC ( index );
      break;
    }
    case kMC:
      AliWarning ( "Not Implemented Yet ..." );
      return ( AliRsnEvent* ) 0x0;
      break;
    case kRSN:
    {
      return GetRsnMVFromRSN();
      break;
    }
    default:
      AliError ( "Type not supported ..." );
      return ( AliRsnEvent* ) 0x0;
      break;
  }
  return ( AliRsnEvent* ) 0x0;
}

void AliRsnAnalysisAT::PostEventProcess ( const Short_t & index )
{
  switch ( fInputType[index] )
  {
    case kAOD:
      break;
    case kESD:
      break;
    case kESDMC:
      break;
    case kMC:
      break;
    case kRSN:
    {
      if ( fRsnMVEventBuffer->GetDeleteBufferWhenReset() == kFALSE )
      {
        fRSN[index] = ( AliRsnEvent* ) fRsnMVEventBuffer->GetNextEvent();
        SetBranchAddress ( 0 , "RsnEvents", &fRSN[index] );
      }
      break;
    }
    default:
      break;
  }

}

AliRsnEvent * AliRsnAnalysisAT::GetRsnMVFromAOD ( const Short_t & index )
{

  if ( !fAOD[index] ) { AliError ( "fAOD not available." ); return ( AliRsnEvent * ) 0x0; }


//   fRSN[0] = new AliRsnEvent();
//   fRSN[0]->Init();
//   fRSN[0]->BuildEvent ( fAOD[index] );
//   return fRSN[0];

  return ( AliRsnEvent* ) 0x0;

}

AliRsnEvent * AliRsnAnalysisAT::GetRsnMVFromESD ( const Short_t & index )
{
  if ( !fESD[index] ) { AliError ( "fESD not available." ); return ( AliRsnEvent * ) 0x0; }

//   fRSN[0] = new AliRsnEvent();
//   fRSN[0]->Init();
//   fRSN[0]->BuildEvent ( fESD[index] );
//   return fRSN[0];

  return ( AliRsnEvent* ) 0x0;
}

AliRsnEvent * AliRsnAnalysisAT::GetRsnMVFromESDMC ( const Short_t & index )
{

  if ( !fESD[index] ) { AliError ( "fESD not available." ); return ( AliRsnEvent * ) 0x0; }
  AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> ( AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() );
  if ( !mcHandler ) { AliError ( "Could not retrieve MC event handler" ); return ( AliRsnEvent * ) 0x0; }

//   fRSN[0] = new AliRsnEvent();
//   fRSN[0]->Init();
//   fRSN[0]->BuildEvent ( fESD[index] ,mcHandler );
//   return fRSN[0];

  return ( AliRsnEvent* ) 0x0;
}

AliRsnEvent * AliRsnAnalysisAT::GetRsnMVFromRSN ( const Short_t & index )
{
  AliRsnEvent *event = fRSN[index];
  if ( fRsnMVEventBuffer->GetDeleteBufferWhenReset() == kTRUE )
  {
    event = ( AliRsnEvent * ) fRSN[index]->Clone();
  }
//   AliInfo ( Form ( "%p %p",event,fRSN[index] ) );
  return event;
}
