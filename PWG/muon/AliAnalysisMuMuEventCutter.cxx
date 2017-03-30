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
#include "AliAnalysisMuonUtility.h"
#include "TList.h"
#include "TTree.h"
#include "Riostream.h"
#include "AliVVertex.h"
#include "AliAODVertex.h"
#include "AliVVZERO.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "TMath.h"
#include "TParameter.h"
#include "AliESDTZERO.h"
#include "AliAODTZERO.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenCocktailEventHeader.h"

ClassImp(AliAnalysisMuMuEventCutter)

//______________________________________________________________________________
AliAnalysisMuMuEventCutter::AliAnalysisMuMuEventCutter(TRootIOCtor* /*ioCtor*/)
: TObject(), fMuonEventCuts(0x0)
{
  /// default io ctor
}

//______________________________________________________________________________
AliAnalysisMuMuEventCutter::AliAnalysisMuMuEventCutter(const char* triggerClasses, const char* triggerInputsMap)
: TObject(), fMuonEventCuts(0x0)
{
  /// ctor
  TString tclasses(triggerClasses);

  if ( !triggerClasses )
  {
    tclasses = "ANY";
  }

  TString tinputs(triggerInputsMap);

  if ( !triggerInputsMap )
  {
    tinputs = "";
  }

  MuonEventCuts()->SetTrigClassPatterns(tclasses,tinputs);
}

//______________________________________________________________________________
AliAnalysisMuMuEventCutter::AliAnalysisMuMuEventCutter(TList* triggerClasses, TList* triggerInputsMap)
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

  TString tinputs;

  if ( !triggerInputsMap )
  {
    tinputs = "";
  }
  else
  {
    TObjString* tinputsname;
    TIter next(triggerInputsMap);

    while ( ( tinputsname = static_cast<TObjString*>(next()) ) )
    {
      if (tinputs.Length()>0)
      {
        tinputs += ",";
      }

      tinputs += tinputsname->String();
    }
  }

  MuonEventCuts()->SetTrigClassPatterns(tclasses,tinputs);

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

// //_____________________________________________________________________________
// Bool_t AliAnalysisMuMuEventCutter::SelectTriggerClassWithInputHandler(const AliInputEventHandler& eventHandler,
//                                                       TString& acceptedClasses) const
// {
//   /// Forward the trigger class selection to MuonEventCuts::GetSelectedTrigClassesInEvent
//   acceptedClasses = "";
//
//   TIter next(MuonEventCuts()->GetSelectedTrigClassesInEvent(&eventHandler));
//   TObjString* str;
//
//   while ( ( str = static_cast<TObjString*>(next()) ) )
//   {
//     acceptedClasses += str->String();
//     acceptedClasses += " ";
//   }
//   return (acceptedClasses.Length()>0);
// }

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuEventCutter::IsPhysicsSelectedANY(const AliInputEventHandler& eventHandler) const
{
  /// Whether or not the event is physics selected
  return const_cast<AliInputEventHandler&>(eventHandler).IsEventSelected() & AliVEvent::kAny;
}

Bool_t AliAnalysisMuMuEventCutter::IsPhysicsSelectedINT7(const AliInputEventHandler& eventHandler) const
{
  /// Whether or not the event is physics selected
  return const_cast<AliInputEventHandler&>(eventHandler).IsEventSelected() & AliVEvent::kINT7;
}

Bool_t AliAnalysisMuMuEventCutter::IsPhysicsSelectedINT8(const AliInputEventHandler& eventHandler) const
{
  /// Whether or not the event is physics selected
  return const_cast<AliInputEventHandler&>(eventHandler).IsEventSelected() & AliVEvent::kINT8;
}

Bool_t AliAnalysisMuMuEventCutter::IsPhysicsSelectedMUL(const AliInputEventHandler& eventHandler) const
{
  /// Whether or not the event is physics selected
  return const_cast<AliInputEventHandler&>(eventHandler).IsEventSelected() & AliVEvent::kMuonUnlikeLowPt7;
}

Bool_t AliAnalysisMuMuEventCutter::IsPhysicsSelectedMULORMLL(const AliInputEventHandler& eventHandler) const
{
  /// Whether or not the event is physics selected
  return const_cast<AliInputEventHandler&>(eventHandler).IsEventSelected() & ( AliVEvent::kMuonUnlikeLowPt7 | AliVEvent::kMuonLikeLowPt7 );
}

Bool_t AliAnalysisMuMuEventCutter::IsPhysicsSelectedINT7inMUON(const AliInputEventHandler& eventHandler) const
{
  /// Whether or not the event is physics selected
  return const_cast<AliInputEventHandler&>(eventHandler).IsEventSelected() & AliVEvent::kINT7inMUON;
}

Bool_t AliAnalysisMuMuEventCutter::IsPhysicsSelectedMSL(const AliInputEventHandler& eventHandler) const
{
  /// Whether or not the event is physics selected
  return const_cast<AliInputEventHandler&>(eventHandler).IsEventSelected() & AliVEvent::kMuonSingleLowPt7;
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

    if ( ( v0sum > 11.5 && v0sum < 17.5 ) && ( v0diff > 5.5 && v0diff < 11.5 ) )
    {
      return kTRUE;
    }
  }
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuEventCutter::IsMCEventNSD(const AliVEvent& event) const
{
  // Check for headers
  //FIXME: For now this is only valid for DPMJET
  //FIXME: working only on AODs

  if ( static_cast<const AliVEvent*>(&event)->IsA() == AliESDEvent::Class() )
  {
    AliWarning("Not implemented for ESDs yet");
    return kFALSE;
  }

  const AliAODEvent* eventAOD = static_cast<const AliAODEvent*>(&event);

  AliAODMCHeader* mcHeader = static_cast<AliAODMCHeader*>(eventAOD->FindListObject(AliAODMCHeader::StdBranchName()));

  if(mcHeader)
  {
    TList* lheaders = mcHeader->GetCocktailHeaders();

    if ( lheaders->GetEntries() > 1 ) AliWarning("There is more than one header: The simulation is a cocktail");

    AliGenEventHeader* mcGenH(0x0);
    AliGenHijingEventHeader* hHijing(0x0);
    AliGenDPMjetEventHeader* hDpmJet(0x0);
    TIter next(lheaders); // Get the iterator on the list of cocktail headers

//    lheaders->Print();

    while ( (mcGenH = static_cast<AliGenEventHeader*>(next())) ) // Loop over the cocktail headers
    {
//      std::cout << mcGenH->GetName() << std::endl;
      if (mcGenH->InheritsFrom(AliGenHijingEventHeader::Class()))
      {
        hHijing = static_cast<AliGenHijingEventHeader*>(mcGenH);
      }
      if (mcGenH->InheritsFrom(AliGenDPMjetEventHeader::Class()))
      {
        hDpmJet = static_cast<AliGenDPMjetEventHeader*>(mcGenH);
      }
    } // End of loop over cocktail headers

    if ( !hDpmJet && !hHijing )
    {
      AliError("No GenHeader found");
      return kFALSE;
    }
    else if ( hDpmJet )
    {
      Int_t nsd1,nsd2,ndd;
      Int_t npProj = hDpmJet->ProjectileParticipants();
      Int_t npTgt = hDpmJet->TargetParticipants();

      // In p-Pb collisions: npProj >=1 (Pb) and npTgt = 1 (p) ->We have to use npProj
      // In Pb-p collisions: npProj =1 (p) and npTgt >= 1 (Pb) ->We have to use npTgt
      if ( npTgt >= npProj ) npProj = npTgt; // In order to use the correct value to compare to nsd1 and nsd2 (If in Pb-p we use npProj(=1), as soon as we have 1 nucleon difracted the event is flagged as SD, which is not correct)

//      Int_t pType = hDpmJet->ProcessType();
//      std::cout << "Nof Proj part = " << npProj << " ; Nof Tgt part = " << npTgt << std::endl;
      hDpmJet->GetNDiffractive(nsd1,nsd2,ndd);

      if (ndd==0 && (npProj==nsd1 || npProj==nsd2))
      {
        return kFALSE; // reject SD
      }
      else return kTRUE; //Is NSD event

    }
    else
    {
      AliWarning("Implement the Hijing section");
      return kFALSE;
    }
  }
  else
  {
    AliError("No MCheader in MCEvent");
    return kFALSE;
  }

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

  Double_t SPDzv(0.);
  Bool_t vertexFound(kFALSE);

  const AliVVertex* SPDVertex = event.GetPrimaryVertexSPD();
  if ( SPDVertex )
  {
    vertexFound = kTRUE;
    SPDzv = SPDVertex->GetZ();
  }

  if ( !vertexFound )
  {
    AliError("SPD |z| cut requested and no SPD vertex found in the event");
    return kFALSE;
  }

  else return (TMath::Abs(SPDzv)<z);
}


//_____________________________________________________________________________
Bool_t AliAnalysisMuMuEventCutter::IsSPDzVertexInRange(AliVEvent& event, const Double_t& zMin, const Double_t& zMax) const
{
  /// Whether or not the SPD Z vertex is in the range [zMin,zMax[

  const AliVVertex* SPDVertex = event.GetPrimaryVertexSPD();

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
  const AliVVertex* SPDVertex = event.GetPrimaryVertexSPD();
  if ( SPDVertex && SPDVertex->GetNContributors() > 0) return kTRUE;
  else return kFALSE;
}


//_____________________________________________________________________________
Bool_t AliAnalysisMuMuEventCutter::IsSPDzQA(const AliVEvent& event, /*const AliVVertex& vertex2Test,*/ const Double_t& zResCut, const Double_t& zDifCut) const
{
  // Checks if the value of the Z component of the given vertex fullfills the quality assurance condition
  Double_t zRes,zvertex;

  const AliVVertex* vertex = event.GetPrimaryVertex();

  if ( vertex )
  {
    const AliVVertex* SPDVertex = event.GetPrimaryVertexSPD();
    if ( SPDVertex )
    {
     if ( SPDVertex->GetNContributors() > 0 && !SPDVertex->IsFromVertexerZ() )
     {
       Double_t cov[6]={0};
       SPDVertex->GetCovarianceMatrix(cov);
       zRes = TMath::Sqrt(cov[5]);
       zvertex = SPDVertex->GetZ();

         if ( (zRes <= zResCut) && TMath::Abs(zvertex - vertex->GetZ()) <= zDifCut )
         {
           return kTRUE;
         }
         else return kFALSE;
       }
       else return kFALSE;
    }
    else
    {
      AliError("Cut on SPD z Vertex requested for an event with no SPD vertex info");
      return kFALSE;
    }
  }
  else
  {
    AliError("Cut on SPD z Vertex requested for an event with no vertex info");
    return kFALSE;
  }
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
  Bool_t parFound(kFALSE);
  TParameter<Double_t>* eventdNchdEta;

  while ( !parFound )
  {
    while ( nchList->At(i)->IsA() != TParameter<Double_t>::Class() )
    {
      i++;
    }

    eventdNchdEta = static_cast<TParameter<Double_t>*>(nchList->At(i));

    if ( TString(eventdNchdEta->GetName()).Contains("MeandNchdEta") ) parFound = kTRUE;
  }

  Double_t meandNchdEta = eventdNchdEta->GetVal();

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
