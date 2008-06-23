#include "AliLog.h"

#include "AliRsnEvent.h"
#include "AliRsnPairParticle.h"
#include "AliRsnPairDef.h"

#include "AliRsnPair.h"

ClassImp ( AliRsnPair )

AliRsnPair::AliRsnPair()
    : TObject() ,fPairDef(),
    fPairType ( kESDNoPID )
{
  fCurrentCutMgr=0;
  fMass[0]=-1.0;
  fMass[1]=-1.0;
  fRsnMVEventBuffer=0;
  fNumOfMixEvent=0;
  fIsSignSame = kFALSE;
  fIsFilledOnlyInHistRange = kTRUE;

}

AliRsnPair::AliRsnPair ( AliRsnPair::EPairType type, AliRsnPairDef * pairDef,Int_t numOfMix) : TObject(),
    fPairDef ( *pairDef ),
    fPairType ( type )
{
  fCurrentCutMgr=0;
  fMass[0]=-1.0;
  fMass[1]=-1.0;
  fRsnMVEventBuffer=0;
  fNumOfMixEvent=numOfMix;
  fIsSignSame = kFALSE;
  fIsFilledOnlyInHistRange = kTRUE;
}

AliRsnPair::~AliRsnPair()
{
}

TString AliRsnPair::GetEffMassHistName ( Int_t index )
{

  fCurrentCutMgr = ( AliRsnCutMgr* ) fCutMgrs.UncheckedAt ( index );

  TString sName;
  sName += GetPairTypeName ( fPairType );
  sName += GetESDParticleName ( fPairDef.GetType( 0 ) );
  sName += fPairDef.GetCharge ( 0 );
  sName += GetESDParticleName ( fPairDef.GetType( 1 ) );
  sName += fPairDef.GetCharge ( 1 );
  sName += "_";
  if ( fCurrentCutMgr )
    sName += fCurrentCutMgr->GetName();
  else
    sName += "NoCut";
  sName += "_";
  sName += "[";
  sName += Form ( "%.2f", fPairDef.GetMin() );
  sName += "-";
  sName += Form ( "%.2f",fPairDef.GetMax() );
  sName += "]";

  return sName;
}

TString AliRsnPair::GetEffMassHistTitle ( Int_t index )
{

  fCurrentCutMgr= ( AliRsnCutMgr* ) fCutMgrs.UncheckedAt ( index );

  TString sTitle;
  sTitle += GetPairTypeName ( fPairType );
  sTitle += GetESDParticleName ( fPairDef.GetType ( 0 ) );
  sTitle += fPairDef.GetCharge ( 0 );
  sTitle += GetESDParticleName ( fPairDef.GetType ( 1 ) );
  sTitle += fPairDef.GetCharge ( 1 );
  sTitle += " ";
  if ( fCurrentCutMgr )
    sTitle += fCurrentCutMgr->GetTitle();
  else
    sTitle += "NoCut";
  return sTitle;
}

TH1F * AliRsnPair::GenerateEffMassHist ( Int_t index )
{
  return new TH1F ( GetEffMassHistName ( index ).Data(),GetEffMassHistTitle ( index ).Data(),fPairDef.GetNBins(),fPairDef.GetMin(),fPairDef.GetMax() );

}

TString AliRsnPair::GetESDParticleName (  AliRsnPID::EType type )
{

  switch ( type )
  {
    case AliRsnPID::kElectron : return ( "e" );break;
    case AliRsnPID::kMuon     : return ( "mu" );break;
    case AliRsnPID::kPion     : return ( "pi" );break;
    case AliRsnPID::kKaon     : return ( "K" );break;
    case AliRsnPID::kProton   : return ( "p" );break;
    case AliRsnPID::kUnknown  : return ( "unknown" );
    default:
      AliWarning ( "Unrecognized value of EParticle argument" );
      break;
  }

  return "";

}

TString AliRsnPair::GetPairTypeName ( EPairType type )
{
  switch ( type )
  {
    case kESDNoPID : return ( "ESDNOPID_" );break;
    case kESDNoPIDMix : return ( "ESDNOPIDMIX_" );break;
    case kESDNormal : return ( "ESDNORMAL_" );break;
    case kESDMix : return ( "ESDMIX_" );break;
    case kMCNoPID : return ( "MCNOPID_" );break;
    case kMCNormal : return ( "MCNORMAL_" );break;
    case kMCMix : return ( "MCMIX_" );break;
    case kMCSignalOnly : return ( "MCSIGNAL_" );break;
    case kMCBackgroundOnly : return ( "MCBKGONLY_" );break;
    default:
      AliWarning ( "Unrecognized value of EPairTypeName argument" );
      break;
  }

  return "NOTYPE";
}

void AliRsnPair::AddCutMgr ( AliRsnCutMgr * theValue )
{
  fCutMgrs.Add ( theValue );
}

void AliRsnPair::ProcessPair ( AliRsnEvent * event,TH1F*hist ,Int_t index )
{
  AliDebug ( AliLog::kDebug+2,"<-" );

  switch ( fPairType )
  {
    case kESDNoPID :
      DoESDNoPID ( event,hist,index );
      break;
    case kESDNoPIDMix :
      DoESDNoPIDMix ( event,hist,index );
      break;
    case kESDNormal :
      DoESDNormal ( event,hist,index );
      break;
    case kESDMix :
      DoESDMix ( event,hist,index );
      break;
    case kMCNoPID :
      DoMCNoPID ( event,hist,index );
      break;
    case kMCNormal :
      DoMCNormal ( event,hist,index );
      break;
    case kMCSignalOnly :
      DoMCNormal ( event,hist,index );
      break;
    case kMCBackgroundOnly :
      DoMCNormal ( event,hist,index );
      break;
    default:
      AliWarning ( "Wrong fPaitType Skipping pair..." );
      break;
  }

  AliDebug ( AliLog::kDebug+2,"->" );
}

void AliRsnPair::DoCleanUpAfterOneEvent()
{
}

void AliRsnPair::DoLoopPairESD ( AliRsnEvent * event1, TArrayI * array1, AliRsnEvent * event2, TArrayI * array2, TH1F * hist, Int_t index )
{
  AliDebug ( AliLog::kDebug+2,"<-" );

  AliDebug ( AliLog::kDebug+2,Form ( "NumArray1 = %d\tNumArray2 = %d",array1->GetSize(),array2->GetSize() ) );

  fCurrentCutMgr= ( AliRsnCutMgr* ) fCutMgrs.UncheckedAt ( index );
  Int_t startj=0;
  Double_t effMass=0;
  Double_t histMin = fPairDef.GetMin();
  Double_t histMax = fPairDef.GetMax();
  AliRsnDaughter *daughter1=0;
  AliRsnDaughter *daughter2=0;
  Int_t  howManuFilledHist=0;
  for ( Int_t i=0;i<array1->GetSize();i++ )
  {


    daughter1 = ( AliRsnDaughter * ) event1->GetTrack(array1->At ( i ) );
    if ( !daughter1 ) continue;

    if ( fCurrentCutMgr )
      if ( ! ( fCurrentCutMgr->IsSelected ( AliRsnCut::kParticle,daughter1 ) ) ) continue;

    daughter2 = 0;
    if ( fIsSignSame ) startj=i+1;
    for ( Int_t j=startj;j<array2->GetSize();j++ )
    {
      daughter2 = ( AliRsnDaughter * ) event2->GetTrack(array2->At ( j ) );
      if ( !daughter2 ) continue;

      if ( fCurrentCutMgr )
        if ( ! ( fCurrentCutMgr->IsSelected ( AliRsnCut::kParticle,daughter2 ) ) ) continue;

//       AliRsnPairParticle effMassPart;
      fEffMassParticle.FillPairParticle( daughter1 ,daughter2 );
      if ( fCurrentCutMgr )
        if ( ! ( fCurrentCutMgr->IsSelected ( AliRsnCut::kPair,&fEffMassParticle ) ) ) continue;


      if ( fMass[0]<0.0&&fMass[1]<0.0 )
        effMass =  fEffMassParticle.GetESDEffectiveMass() ;
      else
        effMass =  fEffMassParticle.GetESDEffectiveMass ( fMass[0],fMass[1] ) ;

//       if ( fIsSignSame )
//         if ( effMass < 0.988 )
//       {
      //
//         effMassPart.PrintInfo("pt");
//         AliInfo(Form("effMass = %f",effMass));
//       }


      if ( fIsFilledOnlyInHistRange )
        if ( ! ( ( effMass>=histMin ) && ( effMass<=histMax ) ) ) continue;

      hist->Fill ( effMass );

      howManuFilledHist++;
      AliDebug ( AliLog::kDebug+2,Form ( "i=%d j=%d",i,j ) );
    }
  }

//   AliInfo (Form ( "%d",tmpNum));

  AliDebug ( AliLog::kDebug+2,Form ( "NumOfFilledHist = %d",howManuFilledHist ) );

  AliDebug ( AliLog::kDebug+2,"->" );
}

void AliRsnPair::DoLoopPairMC ( AliRsnEvent * event1, TArrayI * array1, AliRsnEvent * event2, TArrayI * array2, TH1F * hist, Int_t index )
{
  AliDebug ( AliLog::kDebug+2,"<-" );

  AliDebug ( AliLog::kDebug+2,Form ( "NumArray1 = %d\tNumArray2 = %d",array1->GetSize(),array2->GetSize() ) );

  fCurrentCutMgr = ( AliRsnCutMgr* ) fCutMgrs.UncheckedAt ( index );
  Int_t startj=0;
  Double_t effMass=0;
  Double_t histMin = fPairDef.GetMin();
  Double_t histMax = fPairDef.GetMax();
  AliRsnDaughter *daughter1=0;
  AliRsnDaughter *daughter2=0;
  Int_t  howManuFilledHist=0;
  Bool_t isSignal=kFALSE;
  for ( Int_t i=0;i<array1->GetSize();i++ )
  {

    daughter1 = ( AliRsnDaughter * ) event1->GetTrack( array1->At ( i ) );
    if ( !daughter1 ) continue;

    if ( fCurrentCutMgr )
      if ( ! ( fCurrentCutMgr->IsSelected ( AliRsnCut::kParticle,daughter1 ) ) ) continue;

    daughter2 = 0;
    if ( fIsSignSame ) startj=i+1;
    for ( Int_t j=startj;j<array2->GetSize();j++ )
    {
      daughter2 = ( AliRsnDaughter * ) event2->GetTrack(array2->At ( j ) );
      if ( !daughter2 ) continue;

      if ( fCurrentCutMgr )
        if ( ! ( fCurrentCutMgr->IsSelected ( AliRsnCut::kParticle,daughter2 ) ) ) continue;

//       AliRsnPairParticle effMassPart;
      fEffMassParticle.FillPairParticle( daughter1 ,daughter2 );
      if ( fCurrentCutMgr )
        if ( ! ( fCurrentCutMgr->IsSelected ( AliRsnCut::kPair,&fEffMassParticle ) ) ) continue;

      
      if ( ( fPairType == kMCSignalOnly ) || ( fPairType == kMCBackgroundOnly ) )
      {
        isSignal = ( daughter1->GetMCInfo()->MotherPDG() == fPairDef.GetMotherPDG() ) && ( daughter2->GetMCInfo()->MotherPDG() == fPairDef.GetMotherPDG() );
        if ( ( fPairType == kMCSignalOnly ) && ( !isSignal ) ) continue;
        if ( ( fPairType == kMCBackgroundOnly ) && ( isSignal ) ) continue;
      }
      
      if ( fMass[0]<0.0&&fMass[1]<0.0 )
        effMass =  fEffMassParticle.GetMCEffectiveMass() ;
      else
        effMass =  fEffMassParticle.GetMCEffectiveMass ( fMass[0],fMass[1] ) ;

//       if ( fIsSignSame )
//         if ( effMass < 0.988 )
//       {
      //
//         effMassPart.PrintInfo("pt");
//         AliInfo(Form("effMass = %f",effMass));
//       }


      if ( fIsFilledOnlyInHistRange )
        if ( ! ( ( effMass>=histMin ) && ( effMass<=histMax ) ) ) continue;

      hist->Fill ( effMass );

      howManuFilledHist++;
      AliDebug ( AliLog::kDebug+2,Form ( "i=%d j=%d",i,j ) );
    }
  }

//   AliInfo (Form ( "%d",tmpNum));

  AliDebug ( AliLog::kDebug+2,Form ( "NumOfFilledHist = %d",howManuFilledHist ) );

  AliDebug ( AliLog::kDebug+2,"->" );
}

void AliRsnPair::DoESDNoPID ( AliRsnEvent * event, TH1F * hist ,Int_t index )
{

  Char_t chargeChar1 = fPairDef.GetCharge ( 0 );
  Char_t chargeChar2 = fPairDef.GetCharge ( 1 );
  
  Int_t chargeIndex1 = fPairDef.GetCharge ( 0 ) =='+'? 0 : 1;
  Int_t chargeIndex2 = fPairDef.GetCharge ( 1 ) =='+'? 0 : 1;
  
  fIsSignSame = ( chargeIndex1 == chargeIndex2 );
  Int_t numOfTracks = event->GetMultiplicity();
//   if (numOfTracks<=0) return;
  
//   AliInfo(Form("%d", numOfTracks));
  Int_t counter1=0,counter2=0;
  TArrayI array1 ( numOfTracks ),array2 ( numOfTracks );
  TArrayI *arraytmp=0;
//   AliRsnDaughter *dtmp=0;
  Int_t j;
  for ( Int_t i=0; i < AliRsnPID::kSpecies ;i++ )
  {
    arraytmp = ( TArrayI* ) event->GetTracksArray ( chargeChar1, (AliRsnPID::EType )i );
    for ( j=0;j< arraytmp->GetSize();j++ )
    {
      array1.AddAt ( arraytmp->At ( j ),counter1++ );
    }
    arraytmp = ( TArrayI* ) event->GetTracksArray ( chargeChar2, (AliRsnPID::EType )i );
    for ( j=0;j< arraytmp->GetSize();j++ )
    {
      array2.AddAt ( arraytmp->At ( j ),counter2++ );
    }
  }

  array1.Set ( counter1 );
  array2.Set ( counter2 );

//   AliInfo(Form("%d %d",array1.GetSize(),array2.GetSize()));
  DoLoopPairESD ( event,&array1,event,&array2,hist,index );

}

void AliRsnPair::DoESDNoPIDMix ( AliRsnEvent * event, TH1F * hist, Int_t index )
{

  Char_t chargeChar1 = fPairDef.GetCharge ( 0 );
  Char_t chargeChar2 = fPairDef.GetCharge ( 1 );
  
  Int_t chargeIndex1 = fPairDef.GetCharge ( 0 ) =='+'? 0 : 1;
  Int_t chargeIndex2 = fPairDef.GetCharge ( 1 ) =='+'? 0 : 1;
  
  fIsSignSame = ( chargeIndex1 == chargeIndex2 );

  Long64_t currBuffIndex = GetRsnMVEventBuffer()->GetEventsBufferIndex();
  if ( currBuffIndex < 2 ) return;

  
  AliRsnEvent *evCurrEvent = ( AliRsnEvent * ) GetRsnMVEventBuffer()->GetEvent ( currBuffIndex-1 );
  if ( ! evCurrEvent ) {AliWarning ( Form ( "Event Not found" ) ); return;}
  
  Int_t numOfTracks = evCurrEvent->GetMultiplicity();
  Int_t counter1=0;
  TArrayI arrayCurrEvent ( numOfTracks );
  TArrayI *arraytmp=0;
  Int_t j;
  for ( Int_t i=0; i < AliRsnPID::kSpecies ;i++ )
  {
    arraytmp = ( TArrayI* ) evCurrEvent->GetTracksArray ( chargeChar1,(AliRsnPID::EType ) i );
    for ( j=0;j< arraytmp->GetSize();j++ )
    {
      arrayCurrEvent.AddAt ( arraytmp->At ( j ),counter1++ );
    }
  }
  
  arrayCurrEvent.Set ( counter1 );
  
  Int_t numMix = 0;
  for ( Int_t i=currBuffIndex-2;i>=0;i-- )
  {
    if ( ++numMix>fNumOfMixEvent ) break;


    AliRsnEvent *evMix = ( AliRsnEvent * ) GetRsnMVEventBuffer()->GetEvent ( i ) ;
//     if ( ! evMix ) {AliWarning ( Form ( "Event evMix Not found" ) ); continue;}
    numOfTracks = evMix->GetMultiplicity();
    Int_t counter2=0;
    TArrayI arrayMix ( numOfTracks );
    arraytmp=0;
    for ( Int_t i=0; i < AliRsnPID::kSpecies ;i++ )
    {
      arraytmp = ( TArrayI* ) evMix->GetTracksArray ( chargeChar2,(AliRsnPID::EType ) i );
      if (arraytmp->GetSize() > numOfTracks) {AliError(Form("%d %d",arraytmp->GetSize(),numOfTracks));continue;}
      for ( j=0;j< arraytmp->GetSize();j++ )
      {
        arrayMix.AddAt ( arraytmp->At ( j ),counter2++ );
      }
    }
    
    arrayMix.Set ( counter2  );
    DoLoopPairESD ( evCurrEvent,&arrayCurrEvent,evMix,&arrayMix,hist,index );
  }
}

void AliRsnPair::DoESDNormal ( AliRsnEvent * event, TH1F * hist ,Int_t index )
{

  Int_t chargeIndex1 = fPairDef.GetCharge ( 0 ) =='+'? 0 : 1;
  Int_t chargeIndex2 = fPairDef.GetCharge ( 1 ) =='+'? 0 : 1;
  fIsSignSame = ( chargeIndex1 == chargeIndex2 );

  TArrayI *array1 = ( TArrayI* ) event->GetTracksArray ( chargeIndex1, fPairDef.GetType( 0 ) );
  TArrayI *array2 = ( TArrayI* ) event->GetTracksArray ( chargeIndex2, fPairDef.GetType( 1 ) );
  
  DoLoopPairESD ( event, array1,event,array2,hist,index );
}

void AliRsnPair::PrepareMixForPair ( AliRsnEvent * event,TTree *tree )
{
}

void AliRsnPair::DoESDMix ( AliRsnEvent * event, TH1F * hist, Int_t index )
{
  Int_t chargeIndex1 = fPairDef.GetCharge ( 0 ) =='+'? 0 : 1;
  Int_t chargeIndex2 = fPairDef.GetCharge ( 1 ) =='+'? 0 : 1;
  fIsSignSame = ( chargeIndex1 == chargeIndex2 );

  Long64_t currBuffIndex = GetRsnMVEventBuffer()->GetEventsBufferIndex();
  if ( currBuffIndex < 2 ) return;

  AliRsnEvent *evCurrEvent = ( AliRsnEvent * ) GetRsnMVEventBuffer()->GetEvent ( currBuffIndex-1 );
  if ( ! evCurrEvent ) {AliWarning ( Form ( "Event Not found" ) ); return;}
  TArrayI  *arrayCurrEvent = ( TArrayI* ) evCurrEvent->GetTracksArray ( chargeIndex1, fPairDef.GetType ( 0 ) );

  Int_t numMix = 0;
  TArrayI* arrayMix=0;
  for ( Int_t i=currBuffIndex-2;i>=0;i-- )
  {
    if ( ++numMix>fNumOfMixEvent ) break;


    AliRsnEvent *evMix = ( AliRsnEvent * ) GetRsnMVEventBuffer()->GetEvent ( i ) ;
    if ( ! evMix ) {AliWarning ( Form ( "Event Not found" ) ); continue;}
    arrayMix = ( TArrayI* ) evMix->GetTracksArray ( chargeIndex2,fPairDef.GetType ( 1 ) );
    DoLoopPairESD ( evCurrEvent,arrayCurrEvent,evMix,arrayMix,hist,index );
  }

}

void AliRsnPair::DoMCNoPID ( AliRsnEvent * event, TH1F * hist ,Int_t index )
{

/*
  Int_t chargeIndex1 = fPairDef.GetCharge ( 0 ) =='+'? 0 : 1;
  Int_t chargeIndex2 = fPairDef.GetCharge ( 1 ) =='+'? 0 : 1;
  fIsSignSame = ( chargeIndex1 == chargeIndex2 );

//   TRefArray array1,array2;
//   TRefArray *arraytmp=0;
//   AliRsnDaughter *dtmp=0;
//   for ( Int_t i=0; i < AliRsnPID::kSpecies ;i++ )
//   {
//     arraytmp = ( TRefArray* ) event->GetMCTrackRefs ( chargeIndex1, i );
//     for ( Int_t j=0;j< arraytmp->GetEntriesFast();j++ )
//     {
//       dtmp = ( AliRsnDaughter * ) arraytmp->At ( j );
//       array1.Add ( dtmp );
//     }
//     arraytmp = ( TRefArray* ) event->GetMCTrackRefs ( chargeIndex2, i );
//     for ( Int_t j=0;j< arraytmp->GetEntriesFast();j++ )
//     {
//       dtmp = ( AliRsnDaughter * ) arraytmp->At ( j );
//       array2.Add ( dtmp );
//     }
//   }
//   DoLoopPairMC ( &array1,&array2,hist,index );

  Int_t numOfTracks = event->GetMultiplicity();
//   AliInfo(Form("%d",numOfTracks));
  Int_t counter1=0,counter2=0;
  TArrayI array1 ( numOfTracks ),array2 ( numOfTracks );
  TArrayI *arraytmp=0;
  Int_t j;
  for ( Int_t i=0; i < AliRsnPID::kSpecies ;i++ )
  {
    arraytmp = ( TArrayI* ) event->GetMCTrackArray ( chargeIndex1, i );
    for ( j=0;j< arraytmp->GetSize();j++ )
    {
      array1.AddAt ( arraytmp->At ( j ),counter1++ );
    }
    arraytmp = ( TArrayI* ) event->GetMCTrackArray ( chargeIndex2, i );
    for ( j=0;j< arraytmp->GetSize();j++ )
    {
      array2.AddAt ( arraytmp->At ( j ),counter2++ );
    }
  }

  array1.Set ( counter1 );
  array2.Set ( counter2 );

//   AliInfo(Form("%d %d",array1.GetSize(),array2.GetSize()));
  DoLoopPairMC ( event,&array1,event,&array2,hist,index );*/

}

void AliRsnPair::DoMCNormal ( AliRsnEvent * event, TH1F * hist ,Int_t index )
{
//   Int_t chargeIndex1 = fPairDef.GetCharge ( 0 ) =='+'? 0 : 1;
//   Int_t chargeIndex2 = fPairDef.GetCharge ( 1 ) =='+'? 0 : 1;
//   fIsSignSame = ( chargeIndex1 == chargeIndex2 );
// //   TRefArray *array1 = ( TRefArray* ) event->GetMCTrackRefs ( chargeIndex1, ( Int_t ) fPairDef.GetESDParticle ( 0 ) );
// //   TRefArray *array2 = ( TRefArray* ) event->GetMCTrackRefs ( chargeIndex2, ( Int_t ) fPairDef.GetESDParticle ( 1 ) );
// //
// //   DoLoopPairMC ( array1,array2,hist,index );
// 
//   TArrayI *array1 = ( TArrayI* ) event->GetMCTrackArray ( chargeIndex1, ( Int_t ) fPairDef.GetESDParticle ( 0 ) );
//   TArrayI *array2 = ( TArrayI* ) event->GetMCTrackArray ( chargeIndex2, ( Int_t ) fPairDef.GetESDParticle ( 1 ) );
//   DoLoopPairMC ( event, array1, event,array2,hist,index );
// 

}


