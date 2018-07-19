#include "AliAnalysisNanoAODCuts.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliNanoAODHeader.h"
#include "AliNanoAODTrack.h"
#include "AliMultSelection.h"
#include <iomanip>

ClassImp(AliAnalysisNanoAODTrackCuts)
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

AliAnalysisNanoAODEventCuts::AliAnalysisNanoAODEventCuts():
  AliAnalysisCuts(), 
  fVertexRange(-1),
  fTrackCut(0),
  fMinMultiplicity(-1),
  fMaxMultiplicity(-1)
{
  // default ctor   
}

Bool_t AliAnalysisNanoAODEventCuts::IsSelected(TObject* obj)
{
  // Returns true if object accepted on the event level
  
  AliAODEvent * evt = dynamic_cast<AliAODEvent*>(obj);
  
  if (fVertexRange > 0)
  {
    AliAODVertex * vertex = evt->GetPrimaryVertex();
    if (!vertex) 
      return kFALSE;
    
    if (vertex->GetNContributors() < 1) 
    {
      // SPD vertex cut
      vertex = evt->GetPrimaryVertexSPD();    
      if (!vertex || vertex->GetNContributors() < 1) 
        return kFALSE;
    }    
    
    if (TMath::Abs(vertex->GetZ()) > fVertexRange) 
      return kFALSE;
  }
  
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


void AliNanoAODSimpleSetter::SetNanoAODHeader(const AliAODEvent * event   , AliNanoAODHeader * head , TString varListHeader  ) {

  AliAODHeader * header = dynamic_cast<AliAODHeader*>(event->GetHeader());
  if (!header) AliFatal("Not a standard AOD");

  // Set custom nano aod vars
  Double_t centrV0M=-1;
  Double_t centrTRK=-1;
  Double_t centrCL1=-1;
  Double_t centrCL0=-1;

      //2015 multiplicity selection
      AliMultSelection *MultSelection = 0x0; 
      MultSelection = (AliMultSelection *) event->FindListObject("MultSelection");

      if(MultSelection){
        centrV0M = MultSelection->GetMultiplicityPercentile("V0M");
        centrCL1 = MultSelection->GetMultiplicityPercentile("CL1");
        centrCL0 = MultSelection->GetMultiplicityPercentile("CL0");
        centrTRK = MultSelection->GetMultiplicityPercentile("TRK");

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

  static const char * validatorString[] = {"Centr","CentrTRK","CentrCL0","CentrCL1", "MagField", "OfflineTrigger", "RunNumber", "FiredTriggerClasses", "BunchCrossNumber", "OrbitNumber", "PeriodNumber", 0};
  TObjArray * vars = varListHeader.Tokenize(",");
  //Int_t size = vars->GetSize();
  TIter it(vars);
  TObjString *token  = 0;
  Int_t index=0;
  Int_t indexInt=0;

  std::map<TString,int> cstMap = head->GetMapCstVar();

  while ((token = (TObjString*) it.Next())) {
    TString var = token->GetString().Strip(TString::kBoth, ' ');

    // Check if string is in the allowed list
    Bool_t isValid = kFALSE;
    Int_t ivalidator = 0;
    while (validatorString[ivalidator]) {
      if(var == validatorString[ivalidator++]) isValid = kTRUE;
    }

    if (!( isValid || var.BeginsWith("cst"))) AliFatal(Form("Invalid var [%s]", var.Data()));

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
    else {
      cstMap[var] = index;
      std::cout << "ADDING " << index << " " << cstMap[var] << " " << var.Data() << std::endl;
      
    }

    index++;
  }
  //size = index;
  if(vars){
    vars->Delete();
    delete vars;
  }
  head->SetMapCstVar(cstMap);

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
