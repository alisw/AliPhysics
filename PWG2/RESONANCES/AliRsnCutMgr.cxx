#include "AliLog.h"

#include "AliRsnDaughter.h"
#include "AliRsnCutSet.h"
#include "AliRsnCutMgr.h"

ClassImp ( AliRsnCutMgr )

AliRsnCutMgr::AliRsnCutMgr()
    : TNamed("defaultName","default Tilte")
{
  for ( Int_t i=0 ;i<AliRsnCut::kLastCutSetIndex ;i++ )
  {
    fCutSets[i] = 0;
  }

}

AliRsnCutMgr::AliRsnCutMgr ( const char * name, const char * title )
    : TNamed ( name,title )
{

  for ( Int_t i=0 ;i<AliRsnCut::kLastCutSetIndex ;i++ )
  {
    fCutSets[i] = 0 ;
  }

}

AliRsnCutMgr::~AliRsnCutMgr()
{
  for ( Int_t i=0 ;i<AliRsnCut::kLastCutSetIndex ;i++ )
  {
    delete fCutSets[i];
  }

}

void AliRsnCutMgr::SetCutSet ( AliRsnCut::ECutSetType type, AliRsnCutSet* cutset )
{
  if ( !fCutSets[type] )
    fCutSets[type] = ( AliRsnCutSet* ) cutset->Clone() ;
  AliDebug (AliLog::kDebug ,Form ( "DatasetName %s",fCutSets[type]->GetName() ) );
}


Bool_t AliRsnCutMgr::IsSelected ( AliRsnCut::ECutSetType type,TObject *obj )
{
  AliDebug ( AliLog::kDebug,"<-" );

  if (!fCutSets[type]) return kTRUE;

  switch ( type )
  {
    case AliRsnCut::kParticle:
//       AliDebug ( AliLog::kDebug,"kParticle" );
      return fCutSets[type]->IsSelected (type, ( AliRsnDaughter * ) obj );
      break;
    case AliRsnCut::kPair:
//       AliInfo ( Form("%p" ,fCutSets[type]));
      return fCutSets[type]->IsSelected (type, ( AliRsnPairParticle * ) obj );
      break;
    default:
      AliWarning ( "Wrong ECutSetType selected." );
      return kTRUE;
      break;

  }
  return kTRUE;
}
