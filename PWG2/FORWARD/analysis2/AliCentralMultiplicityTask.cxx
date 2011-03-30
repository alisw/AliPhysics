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
#include "AliFMDEventInspector.h"
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
    fManager(),
    fUseSecondary(true),
    fUseAcceptance(true),
    fFirstEventSeen(false)
{
  // 
  // Constructor 
  //   
  DefineOutput(1, TList::Class());
  fBranchNames = 
    "ESD:AliESDRun.,AliESDHeader.,AliMultiplicity.,"
    "SPDVertex.,PrimaryVertex.";
}
//____________________________________________________________________
AliCentralMultiplicityTask::AliCentralMultiplicityTask() 
  : AliAnalysisTaskSE(),
    fData(0),
    fList(0),
    fAODCentral(),
    fManager(),
    fUseSecondary(true),
    fUseAcceptance(true),
    fFirstEventSeen(false)
{
  // 
  // Constructor 
  // 
}
//____________________________________________________________________
AliCentralMultiplicityTask::AliCentralMultiplicityTask(const AliCentralMultiplicityTask& o)
  : AliAnalysisTaskSE(o),
    fData(o.fData),
    fList(o.fList),
    fAODCentral(o.fAODCentral),
    fManager(o.fManager),
    fUseSecondary(o.fUseSecondary),
    fUseAcceptance(o.fUseAcceptance),
    fFirstEventSeen(o.fFirstEventSeen)
{
  //
  // Copy constructor 
  // 
}
//____________________________________________________________________
AliCentralMultiplicityTask&
AliCentralMultiplicityTask::operator=(const AliCentralMultiplicityTask& o)
{
  // 
  // Assignment operator 
  //
  fData           = o.fData;
  fList           = o.fList;
  fAODCentral     = o.fAODCentral;
  fManager        = o.fManager;
  fUseSecondary   = o.fUseSecondary;
  fUseAcceptance  = o.fUseAcceptance;
  fFirstEventSeen = o.fFirstEventSeen;
  return *this;
}
//____________________________________________________________________
void AliCentralMultiplicityTask::UserCreateOutputObjects() 
{
  // 
  // Create output objects 
  // 
  //

  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah) AliFatal("No AOD output handler set in analysis manager");
  
  
  TObject* obj = &fAODCentral;
  ah->AddBranch("AliAODCentralMult", &obj);
  
  fList = new TList();
  fList->SetOwner();
  PostData(1,fList);
  
}
//____________________________________________________________________
void AliCentralMultiplicityTask::UserExec(Option_t* /*option*/) 
{
  // 
  // Process each event 
  // 
  // Parameters:
  //    option Not used
  //  
  fAODCentral.Clear("");
  
  AliESDInputHandler* eventHandler = 
    dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()
				       ->GetInputEventHandler());
  if (!eventHandler) {
    AliWarning("No inputhandler found for this event!");
    return;}
  
  AliESDEvent* esd = eventHandler->GetEvent();
  
  if(!GetManager().IsInit() && !fFirstEventSeen) {
    AliFMDEventInspector inspector;
    inspector.ReadRunDetails(esd);
    GetManager().Init(inspector.GetCollisionSystem(),
		      inspector.GetEnergy(),
		      inspector.GetField());
    
    AliInfo("Manager of corrections in AliCentralMultiplicityTask init");
    fFirstEventSeen = kTRUE;
  }
    
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
  TH1D* hAcceptance = fManager.GetAcceptanceCorrection(vtxbin);
  TH2D* hSecMap     = fManager.GetSecMapCorrection(vtxbin);
  
  if (!hSecMap)     AliFatal("No secondary map!");
  if (!hAcceptance) AliFatal("No acceptance!");
    
  if (fUseSecondary && hSecMap) aodHist->Divide(hSecMap);
  
  for(Int_t nx = 1; nx <= aodHist->GetNbinsX(); nx++) {
    Float_t accCor = hAcceptance->GetBinContent(nx);
    Float_t accErr = hAcceptance->GetBinError(nx);
    
    Bool_t etabinSeen = kFALSE;  
    for(Int_t ny = 1; ny <= aodHist->GetNbinsY(); ny++) {
      // Get currrent value 
      Float_t aodValue = aodHist->GetBinContent(nx,ny);
      Float_t aodErr   = aodHist->GetBinError(nx,ny);

      // Set underflow bin
      Float_t secCor   = 0;
      if(hSecMap) secCor   = hSecMap->GetBinContent(nx,ny);
      if (secCor > 0.5) etabinSeen = kTRUE;
      if (aodValue < 0.000001) { aodHist->SetBinContent(nx,ny, 0); continue; }

      if (!fUseAcceptance) continue; 

      // Acceptance correction 
      if (accCor   < 0.000001) accCor = 1;
      Float_t aodNew   = aodValue / accCor ;
      Float_t error    = aodNew*TMath::Sqrt(TMath::Power(aodErr/aodValue,2) +
					    TMath::Power(accErr/accCor,2) );
      aodHist->SetBinContent(nx,ny, aodNew);
      //test
      aodHist->SetBinError(nx,ny,error);
      aodHist->SetBinError(nx,ny,aodErr);
      
    }
    //Filling underflow bin if we eta bin is in range
    if(etabinSeen) aodHist->SetBinContent(nx,0, 1.);
  }  

  PostData(1,fList);
}
//____________________________________________________________________
void AliCentralMultiplicityTask::Terminate(Option_t* /*option*/) 
{
  // 
  // End of job
  // 
  // Parameters:
  //    option Not used 
  //
}
//____________________________________________________________________
void
AliCentralMultiplicityTask::Print(Option_t* /*option*/) const
{
  // 
  // Print information 
  // 
  // Parameters:
  //    option Not used
  //
}
//====================================================================
AliCentralMultiplicityTask::Manager::Manager() :
  fAcceptancePath("$ALICE_ROOT/PWG2/FORWARD/corrections/CentralAcceptance"),
  fSecMapPath("$ALICE_ROOT/PWG2/FORWARD/corrections/CentralSecMap"),
  fAcceptance(),
  fSecmap(),
  fAcceptanceName("centralacceptance"),
  fSecMapName("centralsecmap"),
  fIsInit(kFALSE)
{
  //
  // Constructor 
  // 
}
//____________________________________________________________________
AliCentralMultiplicityTask::Manager::Manager(const Manager& o) 
  :fAcceptancePath(o.fAcceptancePath),
   fSecMapPath(o.fSecMapPath),
   fAcceptance(o.fAcceptance),
   fSecmap(o.fSecmap),
   fAcceptanceName(o.fAcceptanceName),
   fSecMapName(o.fSecMapName),
   fIsInit(o.fIsInit)
{
  //
  // Copy Constructor 
  // 
}
//____________________________________________________________________
AliCentralMultiplicityTask::Manager&
AliCentralMultiplicityTask::Manager::operator=(const Manager& o)
{
  //
  // Assignment operator  
  // 
  fAcceptancePath = o.fAcceptancePath;
  fSecMapPath     = o.fSecMapPath;
  fAcceptance     = o.fAcceptance;
  fSecmap         = o.fSecmap;
  fAcceptanceName = o.fAcceptanceName;
  fSecMapName     = o.fSecMapName;
  fIsInit         = o.fIsInit;
  return *this;
}

//____________________________________________________________________
const char* 
AliCentralMultiplicityTask::Manager::GetFullFileName(UShort_t what, 
						     UShort_t sys, 
						     UShort_t sNN,  
						     Short_t  field) const
{
  // 
  // Get full path name to object file 
  // 
  // Parameters:
  //    what   What to get 
  //    sys    Collision system
  //    sNN    Center of mass energy 
  //    field  Magnetic field 
  // 
  // Return:
  //    
  //
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
  // 
  // Get the full path name 
  // 
  // Parameters:
  //    what   What to get
  //    sys    Collision system
  //    sNN    Center of mass energy 
  //    field  Magnetic field 
  // 
  // Return:
  //    
  //
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
  // 
  // Get the secondary map
  // 
  // Parameters:
  //    vtxbin 
  // 
  // Return:
  //    
  //
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
  // 
  // Get the acceptance correction 
  // 
  // Parameters:
  //    vtxbin 
  // 
  // Return:
  //    
  //
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
  // 
  // Initialize 
  // 
  // Parameters:
  //    sys    Collision system (1: pp, 2: PbPb)
  //    sNN    Center of mass energy per nucleon pair [GeV]
  //    field  Magnetic field [kG]
  //
  if(fIsInit) ::Warning("Init","Already initialised - overriding...");
  
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
  
  if(fSecmap && fAcceptance) {
    fIsInit = kTRUE;
    ::Info("Init", 
	   "Central Manager initialised for sys %d, energy %d, field %d",sys,sNN,field);

  }
  
}
//
// EOF
//
