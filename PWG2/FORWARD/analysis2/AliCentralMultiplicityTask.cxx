//====================================================================
// 
// Base class for classes that calculate the multiplicity in the
// central region event-by-event
// 
// Inputs: 
//   - AliESDEvent 
//
// Outputs: 
//   - AliAODCentralMult 
// 
// Histograms 
//   
// Corrections used 
#include "AliCentralMultiplicityTask.h"
#include "AliForwardCorrectionManager.h"
#include "AliForwardUtil.h"
#include "AliLog.h"
#include "AliAODHandler.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliMultiplicity.h"
#include <TROOT.h>
#include <TFile.h>
#include <TError.h>
#include <iostream>
#include <iomanip>

//====================================================================
AliCentralMultiplicityTask::AliCentralMultiplicityTask(const char* name) 
  : AliAnalysisTaskSE(name),
    fData(0),
    fList(0),
    fAODCentral(kFALSE),
    fManager()
{
  
  DefineOutput(1, TList::Class());
}
//____________________________________________________________________
AliCentralMultiplicityTask::AliCentralMultiplicityTask() 
  : AliAnalysisTaskSE(),
    fData(0),
    fList(0),
    fAODCentral(),
    fManager()
{
}
//____________________________________________________________________
AliCentralMultiplicityTask::AliCentralMultiplicityTask(const AliCentralMultiplicityTask& o)
  : AliAnalysisTaskSE(o),
    fData(o.fData),
    fList(o.fList),
    fAODCentral(o.fAODCentral),
    fManager(o.fManager)
{
}
//____________________________________________________________________
AliCentralMultiplicityTask&
AliCentralMultiplicityTask::operator=(const AliCentralMultiplicityTask& o)
{
  fData       = o.fData;
  fList       = o.fList;
  fAODCentral = o.fAODCentral;
  fManager    = o.fManager;
  return *this;
}
//____________________________________________________________________
void AliCentralMultiplicityTask::UserCreateOutputObjects() 
{

  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah) AliFatal("No AOD output handler set in analysis manager");
  
  
  TObject* obj = &fAODCentral;
  ah->AddBranch("AliAODCentralMult", &obj);
  
  fList = new TList();
  PostData(1,fList);
  
}
//____________________________________________________________________
void AliCentralMultiplicityTask::UserExec(Option_t* /*option*/) 
{
  
  AliESDInputHandler* eventHandler = 
    dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()
				       ->GetInputEventHandler());
  if (!eventHandler) {
    AliWarning("No inputhandler found for this event!");
    return;}
  
  AliESDEvent* esd = eventHandler->GetEvent();
  
  //Selecting only events with |valid vertex| < 10 cm
  const AliESDVertex* vertex = esd->GetPrimaryVertexSPD();
  if (!vertex) return;
  if(!(vertex->GetStatus())) return;
  if(vertex->GetNContributors() <= 0) return ;
  if(vertex->GetZRes() > 0.1 ) return;
  Double_t vertexXYZ[3]={0,0,0};
  vertex->GetXYZ(vertexXYZ);
  if(TMath::Abs(vertexXYZ[2]) > 10) return;
  
  Double_t delta           = 2 ;
  Double_t vertexBinDouble = (vertexXYZ[2] + 10) / delta;
  //HHD: The vtxbins are 1-10, not 0-9
  Int_t    vtxbin          = Int_t(vertexBinDouble + 1) ; 
  
  // Make sure AOD is filled
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah)  
    AliFatal("No AOD output handler set in analysis manager");
  
  ah->SetFillAOD(kTRUE);
  
  //Doing analysis
  fAODCentral.Clear("");
  TH2D *aodHist = &(fAODCentral.GetHistogram());
  
  const AliMultiplicity* spdmult = esd->GetMultiplicity();
  //Filling clusters in layer 1 used for tracklets...
  for(Int_t j = 0; j< spdmult->GetNumberOfTracklets();j++)
    aodHist->Fill(spdmult->GetEta(j),spdmult->GetPhi(j));

  //...and then the unused ones in layer 1 
  for(Int_t j = 0; j< spdmult->GetNumberOfSingleClusters();j++) 
    aodHist->Fill(-TMath::Log(TMath::Tan(spdmult->GetThetaSingle(j)/2.)),
		  spdmult->GetPhiSingle(j));
  
  
  // Corrections
  TH2D* hSecMap     = fManager.GetSecMapCorrection(vtxbin);
  TH1D* hAcceptance = fManager.GetAcceptanceCorrection(vtxbin);
  if (!hSecMap)     AliFatal("No secondary map!");
  if (!hAcceptance) AliFatal("No acceptance!");
    
  aodHist->Divide(hSecMap);
  
  for(Int_t nx = 1; nx <= aodHist->GetNbinsX(); nx++) {
    Float_t acccor = hAcceptance->GetBinContent(nx);
    
    Bool_t etabinSeen = kFALSE;  
    for(Int_t ny = 1; ny <= aodHist->GetNbinsY(); ny++) {
      Float_t aodvalue = aodHist->GetBinContent(nx,ny);
      Float_t seccor = hSecMap->GetBinContent(nx,ny);
      if(seccor > 0.5) etabinSeen = kTRUE;
      if(aodvalue < 0.000001) { aodHist->SetBinContent(nx,ny, 0); continue; }
      
      Float_t aodnew   = aodvalue / acccor ;
      aodHist->SetBinContent(nx,ny, aodnew);
      Float_t aodErr   = aodHist->GetBinError(nx,ny);
      Float_t accErr   = hAcceptance->GetBinError(nx);
      Float_t error    = aodnew *TMath::Sqrt(TMath::Power(aodErr/aodvalue,2) +
					     TMath::Power(accErr/acccor,2) );
      aodHist->SetBinError(nx,ny,error);
      
    }
    //Filling underflow bin if we eta bin is in range
    if(etabinSeen) aodHist->SetBinContent(nx,0, 1.);
  }  

  PostData(1,fList);
}
//____________________________________________________________________
void AliCentralMultiplicityTask::Terminate(Option_t* /*option*/) 
{
}
//____________________________________________________________________
void
AliCentralMultiplicityTask::Print(Option_t* /*option*/) const
{
}
//====================================================================
AliCentralMultiplicityTask::Manager::Manager() :
  fAcceptancePath("$ALICE_ROOT/PWG2/FORWARD/corrections/CentralAcceptance"),
  fSecMapPath("$ALICE_ROOT/PWG2/FORWARD/corrections/CentralSecMap"),
  fAcceptance(),
  fSecmap(),
  fAcceptanceName("centralacceptance"),
  fSecMapName("centralsecmap")
{


}
//____________________________________________________________________
AliCentralMultiplicityTask::Manager::Manager(const Manager& o) 
  :fAcceptancePath(o.fAcceptancePath),
   fSecMapPath(o.fSecMapPath),
   fAcceptance(o.fAcceptance),
   fSecmap(o.fSecmap),
   fAcceptanceName(o.fAcceptanceName),
   fSecMapName(o.fSecMapName) 
{}
//____________________________________________________________________
AliCentralMultiplicityTask::Manager&
AliCentralMultiplicityTask::Manager::operator=(const Manager& o)
{
  fAcceptancePath = o.fAcceptancePath;
  fSecMapPath     = o.fSecMapPath;
  fAcceptance     = o.fAcceptance;
  fSecmap         = o.fSecmap;
  fAcceptanceName = o.fAcceptanceName;
  fSecMapName     = o.fSecMapName;
  return *this;
}

//____________________________________________________________________
const char* 
AliCentralMultiplicityTask::Manager::GetFullFileName(UShort_t what, 
						     UShort_t sys, 
						     UShort_t sNN,  
						     Short_t  field) const
{
  return Form("%s/%s",
	      what == 0 ? GetSecMapPath() : GetAcceptancePath(), 
	      GetFileName(what, sys, sNN, field));
}

//____________________________________________________________________
const char* 
AliCentralMultiplicityTask::Manager::GetFileName(UShort_t  what ,
						 UShort_t  sys, 
						 UShort_t  sNN,
						 Short_t   field) const
{
  // Must be static - otherwise the data may disappear on return from
  // this member function
  static TString fname = "";
  
  fname = "";
  switch(what) {
  case 0:  fname.Append(fSecMapName.Data());     break;
  case 1:  fname.Append(fAcceptanceName.Data()); break;
  default:
    ::Error("GetFileName", 
	    "Invalid indentifier %d for central object, must be 0 or 1!", what);
    break;
  }
  fname.Append(Form("_%s_%04dGeV_%c%1dkG.root", 
		    AliForwardUtil::CollisionSystemString(sys), 
		    sNN, (field < 0 ? 'm' : 'p'), TMath::Abs(field)));
  
  return fname.Data();
}

//____________________________________________________________________
TH2D* 
AliCentralMultiplicityTask::Manager::GetSecMapCorrection(UShort_t vtxbin) const
{
  if (!fSecmap) { 
    ::Warning("GetSecMapCorrection","No secondary map defined");
    return 0;
  }
  return fSecmap->GetCorrection(vtxbin);
}
//____________________________________________________________________
TH1D* 
AliCentralMultiplicityTask::Manager::GetAcceptanceCorrection(UShort_t vtxbin) 
  const 
{
  if (!fAcceptance) { 
    ::Warning("GetAcceptanceCorrection","No acceptance map defined");
    return 0;
  }
  return fAcceptance->GetCorrection(vtxbin);
}

//____________________________________________________________________
void 
AliCentralMultiplicityTask::Manager::Init(UShort_t  sys, 
					  UShort_t  sNN,
					  Short_t   field) 
{
  TFile fsec(GetFullFileName(0,sys,sNN,field));
  fSecmap = 
    dynamic_cast<AliCentralCorrSecondaryMap*>(fsec.Get(fSecMapName.Data()));  
  if(!fSecmap) {
    ::Error("Init", "no central Secondary Map found!") ;
    return;
  }
  TFile facc(GetFullFileName(1,sys,sNN,field));
  fAcceptance = 
    dynamic_cast<AliCentralCorrAcceptance*>(facc.Get(fAcceptanceName.Data()));
  if(!fAcceptance) {
    ::Error("Init", "no central Acceptance found!") ;
    return;
  }

}
//
// EOF
//
