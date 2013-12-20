#include "AliAnalysisMuMuEventCutter.h"

/**
 *
 * \ingroup pwg-muon-mumu
 *
 * \class AliAnalysisMuMuEventCutter
 *
 * We're grouping here various event cut methods that can be used
 * as AliAnalysisMuMuCutElement. For instance :
 *
 * - \ref IsPhysicsSelected this is the normal physics selection check
 * - \ref IsPhysicsSelectedVDM version of the physics selection used for VdM analysis
 * - \ref IsSPDzVertexInRange selects event with a SPDvertex in a given range
 * - \ref IsAbsZBelowValue selects event with |Zvertex| (being SPD or not) below a given value
 * - \ref IsAbsZBelowValue as above but requesting explicitely the SPD vertex
 * - \ref IsSPDzQA whether the SPD vertex is a good one
 */

#include "TObjString.h"
#include "AliLog.h"
#include "AliMuonEventCuts.h"
#include "TList.h"
#include "Riostream.h"
#include "AliVVertex.h"
#include "AliAODVertex.h"
#include "AliVVZERO.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "TMath.h"
#include "AliESDTZERO.h"
#include "AliAODTZERO.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"

ClassImp(AliAnalysisMuMuEventCutter)

//______________________________________________________________________________
AliAnalysisMuMuEventCutter::AliAnalysisMuMuEventCutter(TList* triggerClasses)
: TObject(), fMuonEventCuts(0x0)
{
  /// ctor
  TString tclasses;
  
  if ( !triggerClasses )
  {
    tclasses = "ANY";
  }
  else
  {
    TObjString* tname;
    TIter next(triggerClasses);
    
    while ( ( tname = static_cast<TObjString*>(next()) ) )
    {
      if (tclasses.Length()>0)
      {
        tclasses += ",";
      }
      
      tclasses += tname->String();
    }
  }
  
  MuonEventCuts()->SetTrigClassPatterns(tclasses);
}

//______________________________________________________________________________
AliAnalysisMuMuEventCutter::~AliAnalysisMuMuEventCutter()
{
  /// dtor
  delete fMuonEventCuts;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuEventCutter::SelectTriggerClass(const TString& firedTriggerClasses,
                                                      TString& acceptedClasses,
                                                      UInt_t L0, UInt_t L1, UInt_t L2) const
{
  /// Forward the trigger class selection to MuonEventCuts::GetSelectedTrigClassesInEvent
  acceptedClasses = "";
  
  TIter next(MuonEventCuts()->GetSelectedTrigClassesInEvent(firedTriggerClasses,L0,L1,L2));
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    acceptedClasses += str->String();
    acceptedClasses += " ";
  }
  return (acceptedClasses.Length()>0);
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuEventCutter::IsPhysicsSelected(const AliInputEventHandler& eventHandler) const
{
  /// Whether or not the event is physics selected
  return const_cast<AliInputEventHandler&>(eventHandler).IsEventSelected() & AliVEvent::kAny;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuEventCutter::IsPhysicsSelectedVDM(const AliVEvent& event) const
{
  // cut used in vdM scans
  
  AliVVZERO* vzero = event.GetVZEROData();
  
  if (vzero)
  {
    Float_t v0a = vzero->GetV0ATime();
    Float_t v0c = vzero->GetV0CTime();
    
    Float_t v0diff = v0a-v0c;
    Float_t v0sum = v0a+v0c;
    
    if ( ( v0sum > 10.5 && v0sum < 18 ) && ( v0diff > 4 && v0diff < 12 ) )
    {
      return kTRUE;
    }
  }
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuEventCutter::IsAbsZBelowValue(const AliVEvent& event, const Double_t& z) const
{
  // Checks if the absolute value of the Z component of the primary vertex is below a certain value
  
  const AliVVertex* vertex = event.GetPrimaryVertex();
  return (TMath::Abs(vertex->GetZ())<=z);
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuEventCutter::IsAbsZSPDBelowValue(const AliVEvent& event, const Double_t& z) const
{
  // Checks if the absolute value of the SPD Z component of the given vertex is below a certain value
  AliAODVertex* SPDVertex = static_cast<const AliAODEvent*>(&event)->GetPrimaryVertexSPD();
  
  if ( !SPDVertex )
  {
    AliError("SPD |z| cut requested and no SPD vertex found in the event");
    return kFALSE;
  }
  
  else return (TMath::Abs(SPDVertex->GetZ())<=z);
}


//_____________________________________________________________________________
Bool_t AliAnalysisMuMuEventCutter::IsSPDzVertexInRange(AliVEvent& event, const Double_t& zMin, const Double_t& zMax) const
{
  /// Whether or not the SPD Z vertex is in the range [zMin,zMax[
  
  AliAODVertex* SPDVertex = static_cast<AliAODEvent*>(&event)->GetPrimaryVertexSPD();
  
  if ( !SPDVertex )
  {
    AliError("Cut on SPD z Vertex requested for an event with no SPD vertex info");
    return kFALSE;
  }
  Double_t zV = SPDVertex->GetZ();
  
  if ( zV >= zMin && zV < zMax ) return kTRUE;
  else return kFALSE;
}


//_____________________________________________________________________________
Bool_t AliAnalysisMuMuEventCutter::HasSPDVertex(AliVEvent& event) const
{
  /// Does the event have a SPD vertex ?
  AliAODVertex* SPDVertex = static_cast<AliAODEvent*>(&event)->GetPrimaryVertexSPD();
  if ( SPDVertex && SPDVertex->GetNContributors() > 0) return kTRUE;
  else return kFALSE;
}


//_____________________________________________________________________________
Bool_t AliAnalysisMuMuEventCutter::IsSPDzQA(const AliVEvent& event, /*const AliVVertex& vertex2Test,*/ const Double_t& zResCut, const Double_t& zDifCut) const
{
  // Checks if the value of the Z component of the given vertex fullfills the quality assurance condition
  
  const AliVVertex* vertex = event.GetPrimaryVertex();
  const AliAODVertex* vertex2Test = static_cast<const AliAODEvent*>(&event)->GetPrimaryVertexSPD();
  
  Double_t cov[6]={0};
  vertex2Test->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  
  if ( (zRes < zResCut) && TMath::Abs(vertex2Test->GetZ() - vertex->GetZ()) <= zDifCut )
  {
    return kTRUE;
  }
  else return kFALSE;
}


//_____________________________________________________________________________
Bool_t AliAnalysisMuMuEventCutter::IsMeandNchdEtaInRange(AliVEvent& event, const Double_t& dNchdEtaMin, const Double_t& dNchdEtaMax) const
{
  TList* nchList = static_cast<TList*>(event.FindListObject("NCH"));
  
  if (!nchList || nchList->IsEmpty())
  {
    AliFatal("No NCH information found in event. Nch analysis MUST be executed to apply a NCH cut");
    return kFALSE;
  }
  
  Int_t i(0);
  
  while ( nchList->At(i)->IsA() != TObjString::Class() ) //Asign a name to find it by name
  {
    i++;
  }
  
  TObjString* value = static_cast<TObjString*>(nchList->At(i));
  Double_t meandNchdEta = value->String().Atof();
  
  if ( meandNchdEta >= dNchdEtaMin && meandNchdEta < dNchdEtaMax ) return kTRUE;
  else return kFALSE;
}


//_____________________________________________________________________________
AliMuonEventCuts*
AliAnalysisMuMuEventCutter::MuonEventCuts() const
{
  /// Return the single instance of AliMuonEventCuts object we're using
  
  if (!fMuonEventCuts)
  {
    fMuonEventCuts = new AliMuonEventCuts("EventCut","");
  }
  return fMuonEventCuts;
}


//_____________________________________________________________________________
void AliAnalysisMuMuEventCutter::NameOfIsSPDzVertexInRange(TString& name, const Double_t& zMin, const Double_t& zMax) const
{
  name.Form("SPDZBTW%3.2f_%3.2f",zMin,zMax);
}
//_____________________________________________________________________________
void AliAnalysisMuMuEventCutter::NameOfIsAbsZBelowValue(TString& name, const Double_t& z) const
{
  name.Form("ABSZLT%3.2f",z);
}

//_____________________________________________________________________________
void AliAnalysisMuMuEventCutter::NameOfIsAbsZSPDBelowValue(TString& name, const Double_t& z) const
{
  name.Form("SPDABSZLT%3.2f",z);
}


//_____________________________________________________________________________
void AliAnalysisMuMuEventCutter::NameOfIsSPDzQA(TString& name, const Double_t& zResCut, const Double_t& zDifCut) const
{
  name.Form("SPDZQA_RES%3.2f_ZDIF%3.2f",zResCut,zDifCut);
}

//_____________________________________________________________________________
void AliAnalysisMuMuEventCutter::NameOfIsMeandNchdEtaInRange(TString& name, const Double_t& dNchdEtaMin, const Double_t& dNchdEtaMax) const
{
  name.Form("MEANDNDETABTW%3.2f_%3.2f",dNchdEtaMin,dNchdEtaMax);
}


//_____________________________________________________________________________
Bool_t AliAnalysisMuMuEventCutter::IsTZEROPileUp(const AliVEvent& event) const
{
  Bool_t pileupFlag(kFALSE);
  
  if ( event.IsA() == AliESDEvent::Class() )
  {
    const AliESDTZERO* tzero = static_cast<AliESDEvent&>(const_cast<AliVEvent&>(event)).GetESDTZERO();
    if ( tzero )
    {
      pileupFlag = tzero->GetPileupFlag();
    }
  }
  else if ( event.IsA() == AliAODEvent::Class() )
  {
    AliAODTZERO* tzero = static_cast<const AliAODEvent&>(event).GetTZEROData();
    if ( tzero )
    {
      pileupFlag = tzero->GetPileupFlag();
    }
  }
  return pileupFlag;
}
