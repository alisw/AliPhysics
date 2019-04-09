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

AliAnalysisNanoAODV0Cuts::AliAnalysisNanoAODV0Cuts()
    : AliAnalysisCuts(),
      fSelectOnFly(false),
      fOnFlyStatus(false),
      fv0pTMin(-1),
      fv0EtaMax(-1),
      fTransverseRadiusMin(-1),
      fTransverseRadiusMax(-1),
      fCPAMin(-1),
      fDCADaugv0VtxMax(-1),
      fDCADaugPrimVtxMin(-1),
      fDaughEtaMax(-1),
      fDaugMinClsTPC(-1),
      fLambdaDaugnSigTPCMax(-1),
      fCheckDaughterPileup(kFALSE)
{
}

Bool_t AliAnalysisNanoAODV0Cuts::IsSelected(TObject* obj) 
{
  // V0 selection
  
  AliAODv0* v0 = dynamic_cast<AliAODv0*>(obj);
  if (!v0)
    AliFatal("Did not pass a V0");

  if (v0->GetNProngs() != 2)
    return false;

  if (v0->GetNDaughters() != 2)
    return false;

  if (fSelectOnFly && v0->GetOnFlyStatus() == fOnFlyStatus)
    return false;

  if (fv0pTMin > 0 && v0->Pt() < fv0pTMin)
    return false;

  if (fv0EtaMax > 0 && TMath::Abs(v0->Eta()) > fv0EtaMax)
    return false;

  if (fTransverseRadiusMin > 0 || fCPAMin > 0) {
    const AliAODEvent* evt = static_cast<const AliAODEvent*>(static_cast<AliAODTrack *>(v0->GetDaughter(0))->GetEvent());
    if (!evt)
      AliFatal("No event but cut on transverse radius or/and pointing angle requested");

    Float_t xvP = evt->GetPrimaryVertex()->GetX();
    Float_t yvP = evt->GetPrimaryVertex()->GetY();
    Float_t zvP = evt->GetPrimaryVertex()->GetZ();
    Double_t vecTarget[3] = { xvP, yvP, zvP };
    Float_t TransverseRadius = v0->DecayLengthXY(vecTarget);
    if (TransverseRadius > fTransverseRadiusMax || TransverseRadius < fTransverseRadiusMin)
      return false;
    if (fCPAMin > 0 && v0->CosPointingAngle(vecTarget) > fCPAMin)
      return false;
  }
  
  if (fDCADaugv0VtxMax > 0 && v0->DcaV0Daughters() > fDCADaugv0VtxMax)
    return false;
  if (fDCADaugPrimVtxMin > 0 && (v0->DcaPosToPrimVertex() < fDCADaugPrimVtxMin || v0->DcaNegToPrimVertex() < fDCADaugPrimVtxMin))
    return false;

  //TODO: For the preselection it does not matter if the assignment positive and
  //negative is correct, during creation of the Nano AOD this should be checked
  AliAODTrack *pTrack = static_cast<AliAODTrack *>(v0->GetDaughter(0));
  AliAODTrack *nTrack = static_cast<AliAODTrack *>(v0->GetDaughter(1));
  if (fDaughEtaMax > 0 && (TMath::Abs(pTrack->Eta()) > fDaughEtaMax || TMath::Abs(nTrack->Eta()) > fDaughEtaMax))
    return false;
    
  if (fDaugMinClsTPC > 0 && (pTrack->GetTPCNcls() < fDaugMinClsTPC || nTrack->GetTPCNcls() < fDaugMinClsTPC))
    return false;

  if (fLambdaDaugnSigTPCMax > 0) {
    static AliPIDResponse* pidResponse = 0;
    if (!pidResponse) {
      AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
      AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
      pidResponse = inputHandler->GetPIDResponse();
      if (!pidResponse)
        AliFatal("No PID response but PID selection requested");
    }
    AliPIDResponse::EDetPidStatus statusPosTPC = pidResponse->CheckPIDStatus(AliPIDResponse::kTPC, pTrack);
    AliPIDResponse::EDetPidStatus statusNegTPC = pidResponse->CheckPIDStatus(AliPIDResponse::kTPC, nTrack);
    if (!(statusPosTPC == AliPIDResponse::kDetPidOk && statusNegTPC == AliPIDResponse::kDetPidOk))
      return false;

    Float_t nSigPosPion = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC, pTrack, AliPID::kPion);
    Float_t nSigPosProton = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC, pTrack, AliPID::kProton);
    Float_t nSigNegPion = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC, nTrack, AliPID::kPion);
    Float_t nSigNegProton = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC, nTrack, AliPID::kProton);

    // if the daughter tracks are not a proton or a pion within loose cuts, the candidate can be rejected
    if (!((nSigPosPion < fLambdaDaugnSigTPCMax && nSigNegProton < fLambdaDaugnSigTPCMax) ||
          (nSigNegPion < fLambdaDaugnSigTPCMax && nSigPosProton < fLambdaDaugnSigTPCMax)))
      return false;
  }
  
  if (fCheckDaughterPileup) {
    if (!((pTrack->GetTOFBunchCrossing() == 0) || (nTrack->GetTOFBunchCrossing() == 0))) {
      Bool_t PileUpPass = false;
      for (Int_t iLay = 0; iLay < 6; ++iLay) {
        //do not use SDD for Pile Up Rejection since detector is too slow
        if (iLay == 2 || iLay == 3)
          continue;

        if (pTrack->HasPointOnITSLayer(iLay) || nTrack->HasPointOnITSLayer(iLay)) {
          PileUpPass = true;
          break;
        }
      }
      if (!PileUpPass)
        return false;
    }
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

  // TPC clusters
  static const Int_t cstTPCClusterInfo21 = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstTPCClusterInfo21");
  if (cstTPCClusterInfo21 != -1)
    nanoTrack->SetVar(cstTPCClusterInfo21, aodTrack->GetTPCClusterInfo(2,1));

  static const Int_t cstIntID = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstIntID");
  if (cstIntID != -1)
    nanoTrack->SetVar(cstIntID, aodTrack->GetID());
}
