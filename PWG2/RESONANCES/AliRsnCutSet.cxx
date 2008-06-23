#include "AliLog.h"

#include "AliRsnCut.h"
#include "AliRsnExpression.h"

#include "AliRsnCutSet.h"

ClassImp ( AliRsnCutSet )

AliRsnCutSet::AliRsnCutSet()
    : TNamed(),fNumOfCuts ( 0 ),
    fCutScheme ( "" ),
    fCutSchemeIndexed ( "" ),
    fBoolValues ( 0 ),fIsScheme ( kFALSE )
{
  fBoolValues = new Bool_t[1];
//   fExpression = new AliRsnExpression ( fCutSchemeIndexed );
  fExpression = 0;
  AliRsnExpression::sCutSet = this;
}

AliRsnCutSet::AliRsnCutSet ( TString name )
    : TNamed ( name,name ),fNumOfCuts ( 0 ),
    fCutScheme ( "" ),
    fCutSchemeIndexed ( "" ),
    fBoolValues ( 0 ),fIsScheme ( kFALSE )
{
  fBoolValues = new Bool_t[1];
  fExpression = 0;
  AliRsnExpression::sCutSet = this;
}

AliRsnCutSet::AliRsnCutSet ( const AliRsnCutSet & copy )
    :    TNamed ( ( TNamed ) copy ),fCuts ( copy.fCuts ),fNumOfCuts ( copy.fNumOfCuts ),
    fCutScheme ( copy.fCutScheme ),
    fCutSchemeIndexed ( copy.fCutSchemeIndexed ),
    fIsScheme ( copy.fIsScheme ),
    fExpression ( copy.fExpression )
{
  AliRsnExpression::sCutSet = this;
}

AliRsnCutSet::~AliRsnCutSet()
{
  delete fExpression;
  delete [] fBoolValues;
}

void AliRsnCutSet::AddCut ( AliRsnCut * cut )
{
  AliDebug ( AliLog::kDebug,"<-" );
  fCuts.Add ( cut );
  fNumOfCuts++;
  if ( fBoolValues )
    delete fBoolValues;
  fBoolValues = new Bool_t[fNumOfCuts];
  for ( Int_t i=0;i<fNumOfCuts ; i++ )
  {
    fBoolValues[i] = kTRUE;
  }

  AliDebug ( AliLog::kDebug,Form ( "%d",fCuts.GetEntriesFast() ) );
  AliDebug ( AliLog::kDebug,"->" );
}

void AliRsnCutSet::ShowCuts ( )
{
//   AliRsnCut *cut;
//   for ( Int_t i=0; i<fCuts.GetEntriesFast() ;i++ )
//   {
//     cut = ( AliRsnCut* ) fCuts.At ( i );
//     AliInfo ( Form ( "%s (\"%s\") [%.2f - %.2f]",cut->GetName(),cut->GetTitle(),cut->GetMin(),cut->GetMax() ) );
//   }
}

Bool_t AliRsnCutSet::IsSelected ( AliRsnCut::ECutSetType type, AliRsnDaughter *daughter )
{
  if ( !fNumOfCuts )
    return kTRUE;

  Bool_t boolReturn = kTRUE;
  AliRsnCut *cut;
  for ( Int_t i=0; i<fNumOfCuts ;i++ )
  {
    cut = ( AliRsnCut* ) fCuts.At ( i );
    fBoolValues[i] = cut->IsSelected ( type,daughter );
//     AliInfo(Form("%s %d",cut->GetName(),fBoolValues[i]));
  }

  if ( fIsScheme )
    boolReturn = Passed();

  return boolReturn;
}

Bool_t AliRsnCutSet::IsSelected ( AliRsnCut::ECutSetType type, AliRsnPairParticle * pair )
{
  if ( !fNumOfCuts )
    return kTRUE;

  Bool_t boolReturn = kTRUE;
  AliRsnCut *cut;
  for ( Int_t i=0; i<fNumOfCuts ;i++ )
  {
    cut = ( AliRsnCut* ) fCuts.At ( i );
    fBoolValues[i] = cut->IsSelected (type,pair );
//     AliInfo(Form("%s %d",cut->GetName(),fBoolValues[i]));
  }

  if ( fIsScheme )
    boolReturn = Passed();

  return boolReturn;
}


void AliRsnCutSet::SetCutScheme ( const TString & theValue )
{
  AliDebug ( AliLog::kDebug,"<-" );
  fCutScheme = theValue;
  SetCutSchemeIndexed ( theValue );
  fIsScheme = kTRUE;
  AliDebug ( AliLog::kDebug,"->" );
}


void AliRsnCutSet::SetCutSchemeIndexed ( TString theValue )
{
  AliDebug ( AliLog::kDebug,"<-" );
//   fCutSchemeIndexed = theValue;
  fCutSchemeIndexed = GetCutSchemeIndexed();
  AliDebug ( AliLog::kDebug,"->" );
}

Int_t AliRsnCutSet::GetIndexByCutName ( TString s )
{
  AliRsnCut *cut;
  for ( Int_t i=0 ;i< fCuts.GetEntriesFast() ;i++ )
  {
    cut = ( AliRsnCut* ) fCuts.At ( i );
//     if ( !cut->GetName().CompareTo ( s ) )
    if ( !s.CompareTo ( cut->GetName() ) )
    {
      return i;
    }
  }

  return -1;
}

Bool_t AliRsnCutSet::Passed()
{
  AliRsnExpression::sCutSet = this;
  if ( !fExpression )
  {
    fExpression = new AliRsnExpression ( fCutSchemeIndexed );
    AliDebug ( AliLog::kDebug,"fExpression was created." );
  }
  return fExpression->Value ( *GetCuts() );
}

Bool_t AliRsnCutSet::IsValidScheme()
{
  return ( ! ( ShowCutScheme().Contains ( "Error" ) ) );
}

TString AliRsnCutSet::ShowCutScheme()
{
  return fExpression->Unparse();
}

Int_t AliRsnCutSet::TestExpression ( TString opt )
{
//   AliRsnCut *cut1 = new AliRsnCut ( "aaa","aaa" );
//   cut1->SetCutValues ( AliRsnCut::kEsdPt,0.0,1.0 );
//   AliRsnCut *cut2 = new AliRsnCut ( "bbb","bbb" );
//   cut2->SetCutValues ( AliRsnCut::kEsdPt,1.,2.0 );
//   AliRsnCut *cut3 = new AliRsnCut ( "ccc","ccc" );
//   cut3->SetCutValues ( AliRsnCut::kEsdPt,2.0,3.0 );
//
//   AliRsnCutSet* set  = new AliRsnCutSet ( "setOne" );
//   set->AddCut ( cut1 );
//   set->AddCut ( cut2 );
//   set->AddCut ( cut3 );
//
//   set->SetCutScheme ( "(aaa&!(ccc))&(bbb&!(ccc))" );
//
//   set->ShowCuts ();

  return 0;
}

void AliRsnCutSet::PrintSetInfo()
{


  AliInfo ( "========== Rsn Cut Set info ==============" );
  AliInfo ( Form ( "Sheme : %s",fCutScheme.Data() ) );
  AliInfo ( Form ( "Sheme : %s",fCutSchemeIndexed.Data() ) );
  AliInfo ( Form ( "Num of Cuts: %d", fCuts.GetEntriesFast() ) );
  AliInfo ( "====== Cuts ======" );
  AliRsnCut *cut;
  for ( Int_t i=0 ;i< fCuts.GetEntriesFast() ;i++ )
  {
    cut = ( AliRsnCut* ) fCuts.At ( i );
    if ( cut )
      AliInfo ( Form ( "%d %d",i,fBoolValues[i] ) );

  }
  AliInfo ( "========== END Rsn Cut Mgr info ==============" );

}

TString AliRsnCutSet::GetCutSchemeIndexed()
{
  AliDebug ( AliLog::kDebug,"<-" );
  TString str ( fCutScheme );
  AliDebug ( AliLog::kDebug,Form ( "Num of cuts %d",fCuts.GetEntriesFast() ) );
  AliRsnCut *cut;
  for ( Int_t i=0; i<fCuts.GetEntriesFast();i++ )
  {
    cut = ( AliRsnCut* ) fCuts.At ( i );
    str.ReplaceAll ( cut->GetName(),Form ( "%d",i ) );
  }
  AliDebug ( AliLog::kDebug,"->" );
  return str;
}

