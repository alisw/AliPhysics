#include "AliAnalysisNanoAODCuts.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliNanoAODHeader.h"
#include "AliNanoAODTrack.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODv0.h"
#include "AliEventCuts.h"

ClassImp(AliAnalysisNanoAODTrackCuts)
ClassImp(AliAnalysisNanoAODV0Cuts)
ClassImp(AliAnalysisNanoAODEventCuts)
ClassImp(AliNanoAODSimpleSetter)


AliAnalysisNanoAODTrackCuts::AliAnalysisNanoAODTrackCuts():
AliAnalysisCuts(), fBitMask(1), fMinPt(0), fMaxEta(10)
{
  // default ctor 
}

Bool_t AliAnalysisNanoAODTrackCuts::IsSelected(TObject* obj)
{
  // Returns true if the track is good!
  AliAODTrack* track = dynamic_cast<AliAODTrack*>(obj);
  
  
  if (!track->TestFilterBit(fBitMask))    return kFALSE;
  if (track->Pt() < fMinPt)               return kFALSE;
  if (TMath::Abs(track->Eta()) > fMaxEta) return kFALSE; 
  
  return kTRUE;  
}

Bool_t AliAnalysisNanoAODV0Cuts::IsSelected(TObject* obj)
{
  // Returns true if the track is good!
  // TODO this is only an example implementation.
  float v0pTmin = 0.2;
  float v0EtaMax = 0.9;
  float TransverseRadiusMin = 0.2;
  float TransverseRadiusMax = 100;
  float CPAmin = 0.98;
  float DCADaugv0VtxMax = 1.5;
  float DCADaugPrimVtxMin = 0.05;
  AliAODv0* v0 = static_cast<AliAODv0*>(obj);
  AliAODEvent* evt = static_cast<AliAODEvent*>(v0->GetEvent());

  if(v0->GetNProngs() != 2) {
    return false;
  }
  if(v0->GetNDaughters() != 2) {
    return false;
  }
  if (v0->GetOnFlyStatus()) {
    return false;
  }
  if (v0->Pt() > v0pTmin) {
    return false;
  }
  float xvP = evt->GetPrimaryVertex()->GetX();
  float yvP = evt->GetPrimaryVertex()->GetY();
  float zvP = evt->GetPrimaryVertex()->GetZ();
  double vecTarget[3] = { xvP, yvP, zvP };
  float TransverseRadius = v0->DecayLengthXY(vecTarget);
  if (TransverseRadius > TransverseRadiusMax
      || TransverseRadius < TransverseRadiusMin) {
    return false;
  }
  if (v0->CosPointingAngle(vecTarget) > CPAmin) {
    return false;
  }
  if (v0->Eta() < v0EtaMax) {
    return false;
  }
  if (v0->DcaV0Daughters() < DCADaugv0VtxMax) {
    return false;
  }
  if (v0->DcaPosToPrimVertex() < DCADaugPrimVtxMin
      || v0->DcaNegToPrimVertex() < DCADaugPrimVtxMin) {
    return false;
  }
  return true;
}

AliAnalysisNanoAODEventCuts::AliAnalysisNanoAODEventCuts():
  AliAnalysisCuts(), 
  fTrackCut(0),
  fMinMultiplicity(-1),
  fMaxMultiplicity(-1),
  fEventCuts()
{
  // default ctor   
  
}

Bool_t AliAnalysisNanoAODEventCuts::IsSelected(TObject* obj)
{
  // Returns true if object accepted on the event level
  
  AliVEvent* evt = dynamic_cast<AliVEvent*>(obj);
  
  if (!fEventCuts.AcceptEvent(evt))
    return kFALSE;
  
  if (fTrackCut != 0)
  {
    Int_t trackCount = 0;
    for (Int_t j=0; j<evt->GetNumberOfTracks(); j++)
      if (fTrackCut->IsSelected(evt->GetTrack(j)))
        trackCount++;
      
    if (fMinMultiplicity > 0 && trackCount < fMinMultiplicity)
      return kFALSE;
    if (fMaxMultiplicity > 0 && trackCount > fMaxMultiplicity)
      return kFALSE;
  }
      
  return kTRUE;
}

void AliNanoAODSimpleSetter::Init(AliNanoAODHeader* head, TString varListHeader)
{
  // Initialize header
  
  TObjArray *vars = varListHeader.Tokenize(",");
  TIter it(vars);
  TObjString *token  = 0;
  Int_t index=0;
  Int_t indexInt=0;
  
  std::map<TString,int> cstMap = head->GetMapCstVar();

  while ((token = (TObjString*) it.Next())) {
    TString var = token->GetString().Strip(TString::kBoth, ' ');

    if(var == "FiredTriggerClasses"    ){ head->SetFiredTriggerClassesIndex(indexInt); indexInt++; continue;}
    else if(var == "BunchCrossNumber"  ){ head->SetBunchCrossNumberIndex(indexInt); indexInt++; continue;}
    else if(var == "OrbitNumber"       ){ head->SetOrbitNumberIndex(indexInt); indexInt++; continue;}
    else if(var == "PeriodNumber"      ){ head->SetPeriodNumberIndex(indexInt); indexInt++; continue;}
    else if(var == "Centr"      ) head->SetCentrIndex      (index);
    else if(var == "CentrTRK"   ) head->SetCentrTRKIndex   (index);
    else if(var == "CentrCL0"   ) head->SetCentrCL0Index   (index);
    else if(var == "CentrCL1"   ) head->SetCentrCL1Index   (index);
    else if(var == "MagField"   ) head->SetMagFieldIndex   (index);
    else if(var == "OfflineTrigger"   ) head->SetOfflineTriggerIndex   (index);
    else if(var == "RunNumber"  ) head->SetRunNumberIndex  (index);
    else if(var.BeginsWith("cst") || var.BeginsWith("MultSelection.")) {
      cstMap[var] = index;
      fMultMap[var(14, var.Length())] = index;
      std::cout << "ADDING " << index << " " << cstMap[var] << " " << var.Data() << std::endl;
    }
    else 
      AliFatal(Form("Invalid var [%s]", var.Data()));

    index++;
  }

  vars->Delete();
  delete vars;

  head->SetMapCstVar(cstMap);
  
  fInitialized = kTRUE; 
}

void AliNanoAODSimpleSetter::SetNanoAODHeader(const AliAODEvent* event, AliNanoAODHeader* head, TString varListHeader) 
{
  if (!fInitialized)
    Init(head, varListHeader);
  
  AliAODHeader * header = dynamic_cast<AliAODHeader*>(event->GetHeader());
  if (!header) AliFatal("Not a standard AOD");

  // Set custom nano aod vars
  Double_t centrV0M=-1;
  Double_t centrTRK=-1;
  Double_t centrCL1=-1;
  Double_t centrCL0=-1;

  //2015 multiplicity selection
  AliMultSelection* MultSelection = dynamic_cast<AliMultSelection*> (event->FindListObject("MultSelection"));

  if(MultSelection){
    centrV0M = MultSelection->GetMultiplicityPercentile("V0M");
    centrCL1 = MultSelection->GetMultiplicityPercentile("CL1");
    centrCL0 = MultSelection->GetMultiplicityPercentile("CL0");
    centrTRK = MultSelection->GetMultiplicityPercentile("TRK");
    
    for (std::map<TString,int>::iterator it = fMultMap.begin(); it != fMultMap.end(); it++)
      head->SetVar(it->second, MultSelection->GetMultiplicityPercentile(it->first));
  }else{
    //2011 
    AliCentrality * centralityObj = header->GetCentralityP();
    centrV0M = centralityObj->GetCentralityPercentile("V0M");
    centrTRK = centralityObj->GetCentralityPercentile("TRK");
    centrCL1 = centralityObj->GetCentralityPercentile("CL1");
    centrCL0 = centralityObj->GetCentralityPercentile("CL0");
  }

  Double_t magfield = header->GetMagneticField();
  UInt_t offlineTrigger = header->GetOfflineTrigger();
  Int_t runNumber = event->GetRunNumber();
  TString firedTriggerClasses = header->GetFiredTriggerClasses();
  UShort_t bunchCrossNumber = header->GetBunchCrossNumber();
  UInt_t orbitNumber = header->GetOrbitNumber();
  UInt_t periodNumber = header->GetPeriodNumber();

  if ((head->GetFiredTriggerClassesIndex())!=-1)  head->SetFiredTriggerClasses(firedTriggerClasses);
  if ((head->GetBunchCrossNumberIndex())!=-1)     head->SetBunchCrossNumber(bunchCrossNumber);
  if ((head->GetOrbitNumberIndex())!=-1)          head->SetOrbitNumber(orbitNumber);
  if ((head->GetPeriodNumberIndex())!=-1)         head->SetPeriodNumber(periodNumber);
  if ((head->GetCentrIndex())!=-1)     head->SetVar(head->GetCentrIndex()    ,           centrV0M );
  if ((head->GetCentrTRKIndex())!=-1)  head->SetVar(head->GetCentrTRKIndex() ,           centrTRK );
  if ((head->GetCentrCL1Index())!=-1)  head->SetVar(head->GetCentrCL1Index() ,           centrCL1 );
  if ((head->GetCentrCL0Index())!=-1)  head->SetVar(head->GetCentrCL0Index() ,           centrCL0 );
  if ((head->GetMagFieldIndex())!=-1)  head->SetVar(head->GetMagFieldIndex() ,           magfield );
  if ((head->GetOfflineTriggerIndex())!=-1)  head->SetVar(head->GetOfflineTriggerIndex() , Double_t(offlineTrigger));
  if ((head->GetRunNumberIndex())!=-1) head->SetVar(head->GetRunNumberIndex(), Double_t(runNumber));
}

void AliNanoAODSimpleSetter::SetNanoAODTrack (const AliAODTrack * aodTrack, AliNanoAODTrack * nanoTrack) 
{
  // Set custom variables in the special track

  // PID
  
  // Initialize PID track mapping
  static const Bool_t bPIDFilled = AliNanoAODTrack::InitPIDIndex();
  
  // PID variables
  if (bPIDFilled)
  {
    // Get the PID info
    static AliPIDResponse* pidResponse = 0;
    if (!pidResponse) {
      AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
      AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
      pidResponse = inputHandler->GetPIDResponse();
    }
    
    if (!pidResponse)
      AliFatal("PID response not available but fields requested");

    for (Int_t r = AliNanoAODTrack::kSigmaTPC; r<AliNanoAODTrack::kLAST; r++) {
      for (Int_t p = 0; p<AliPID::kSPECIESC; p++) {
        Int_t index = AliNanoAODTrack::GetPIDIndex((AliNanoAODTrack::ENanoPIDResponse) r, (AliPID::EParticleType) p);
        if (index == -1)
          continue;
        Double_t value = 0;
        if (r == AliNanoAODTrack::kSigmaTPC)
          value = pidResponse->NumberOfSigmasTPC(aodTrack, (AliPID::EParticleType) p);
        else if (r == AliNanoAODTrack::kSigmaTOF)
          value = pidResponse->NumberOfSigmasTOF(aodTrack, (AliPID::EParticleType) p);

        nanoTrack->SetVar(index, value);
      }
    }
  }

  // additional fields which are to be moved to AliNanoAODTrack
  static const Int_t cstTOFBunchCrossing = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstTOFBunchCrossing");
  if (cstTOFBunchCrossing != -1)
    nanoTrack->SetVar(cstTOFBunchCrossing, aodTrack->GetTOFBunchCrossing());
}
