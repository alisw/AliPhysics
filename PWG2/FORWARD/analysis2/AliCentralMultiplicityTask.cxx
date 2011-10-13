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
#include "AliAODForwardMult.h"
#include "AliForwardUtil.h"
#include "AliLog.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliMultiplicity.h"
#include <TROOT.h>
#include <TFile.h>
#include <TError.h>
#include <TSystem.h>
#include <iostream>
#include <iomanip>

//====================================================================
AliCentralMultiplicityTask::AliCentralMultiplicityTask(const char* name) 
  : AliAnalysisTaskSE(name),
    fInspector("centralEventInspector"),
    fData(0),
    fList(0),
    fAODCentral(kFALSE),
    fManager(),
    fUseSecondary(true),
    fUseAcceptance(true),
    fFirstEventSeen(false), 
    fIvz(0),
    fNClusterTracklet(0),
    fClusterPerTracklet(0),
    fNCluster(0),
    fNTracklet(0),
    fEtaMin(0),
    fEtaMax(0)
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
    fInspector(),
    fData(0),
    fList(0),
    fAODCentral(),
    fManager(),
    fUseSecondary(true),
    fUseAcceptance(true),
    fFirstEventSeen(false), 
    fIvz(0),
    fNClusterTracklet(0),
    fClusterPerTracklet(0),
    fNCluster(0),
    fNTracklet(0),
    fEtaMin(0),
    fEtaMax(0)
{
  // 
  // Constructor 
  // 
}
//____________________________________________________________________
AliCentralMultiplicityTask::AliCentralMultiplicityTask(const AliCentralMultiplicityTask& o)
  : AliAnalysisTaskSE(o),
    fInspector(o.fInspector),
    fData(o.fData),
    fList(o.fList),
    fAODCentral(o.fAODCentral),
    fManager(o.fManager),
    fUseSecondary(o.fUseSecondary),
    fUseAcceptance(o.fUseAcceptance),
    fFirstEventSeen(o.fFirstEventSeen), 
    fIvz(0),
    fNClusterTracklet(o.fNClusterTracklet),
    fClusterPerTracklet(o.fClusterPerTracklet),
    fNCluster(o.fNCluster),
    fNTracklet(o.fNTracklet),
    fEtaMin(o.fEtaMin),
    fEtaMax(o.fEtaMax)      
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
  fInspector         = o.fInspector;
  fData              = o.fData;
  fList              = o.fList;
  fAODCentral        = o.fAODCentral;
  fManager           = o.fManager;
  fUseSecondary      = o.fUseSecondary;
  fUseAcceptance     = o.fUseAcceptance;
  fFirstEventSeen    = o.fFirstEventSeen;
  fIvz               = 0; 
  fNClusterTracklet  = o.fNClusterTracklet;
  fClusterPerTracklet= o.fClusterPerTracklet;
  fNCluster          = o.fNCluster;
  fNTracklet         = o.fNTracklet;
  fEtaMin            = o.fEtaMin;
  fEtaMax            = o.fEtaMax;
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

  fInspector.DefineOutput(fList);

  PostData(1,fList);  
}

//____________________________________________________________________
AliESDEvent*
AliCentralMultiplicityTask::GetESDEvent()
{
  //
  // Get the ESD event. IF this is the first event, initialise
  //
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) {
    AliWarning("No ESD event found for input event");
    return 0;
  }
  
  // IF we've read the first event already, just return the event 
  if (fFirstEventSeen) return esd;
  
  // Read the details of the rung 
  fInspector.ReadRunDetails(esd);

  // If we weren't initialised before (i.e., in the setup), do so now. 
  if (!GetManager().IsInit()) {
    GetManager().Init(fInspector.GetCollisionSystem(),
		      fInspector.GetEnergy(),
		      fInspector.GetField());
    AliInfo("Manager of corrections in AliCentralMultiplicityTask init");
  }

  // Check for existence and get secondary map 
  AliCentralCorrSecondaryMap* secMap = GetManager().GetSecMap();
  if (!secMap) AliFatal("No secondary map defined!");
  const TAxis& vaxis = secMap->GetVertexAxis();

  FindEtaLimits();

  fNClusterTracklet = new TH2D("nClusterVsnTracklet", 
			       "Total number of cluster vs number of tracklets",
			       100, 0, 100, 100, 0, 100);
  fNClusterTracklet->SetDirectory(0);
  fNClusterTracklet->SetXTitle("# of free clusters");
  fNClusterTracklet->SetYTitle("# of tracklets");
  fNClusterTracklet->SetStats(0);
  fList->Add(fNClusterTracklet);

  Int_t    nEta = 80;
  Double_t lEta = 2;
  fClusterPerTracklet = new TH2D("clusterPerTracklet", 
				 "N_{free cluster}/N_{tracklet} vs. #eta", 
				 nEta,-lEta,lEta, 101, -.05, 10.05);
  fClusterPerTracklet->SetDirectory(0);
  fClusterPerTracklet->SetXTitle("#eta");
  fClusterPerTracklet->SetYTitle("N_{free cluster}/N_{tracklet}");
  fClusterPerTracklet->SetStats(0);
  fList->Add(fClusterPerTracklet);

  // Cache histograms 
  fNCluster = new TH1D("cacheCluster", "", nEta,-lEta,lEta);
  fNCluster->SetDirectory(0);
  fNCluster->Sumw2();
		       
  fNTracklet = new TH1D("cacheTracklet", "", nEta,-lEta,lEta);
  fNTracklet->SetDirectory(0);
  fNTracklet->Sumw2();

  // Initialize the inspecto 
  fInspector.Init(vaxis);
  fFirstEventSeen = kTRUE;

  // Print some information 
  Print();

  return esd;
}
//____________________________________________________________________
void
AliCentralMultiplicityTask::MarkEventForStore() const
{
  // Make sure the AOD tree is filled 
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah)  
    AliFatal("No AOD output handler set in analysis manager");

  ah->SetFillAOD(kTRUE);
}

//____________________________________________________________________
void AliCentralMultiplicityTask::FindEtaLimits()
{
  AliCentralCorrSecondaryMap* secMap = GetManager().GetSecMap();
  
  const TAxis& vaxis = secMap->GetVertexAxis();
  
  fEtaMin.Set(vaxis.GetNbins());
  fEtaMax.Set(vaxis.GetNbins());
  
  TList* hits = new TList;
  hits->SetOwner();
  hits->SetName("hitMaps");
  fList->Add(hits);
  
  TList* secs = new TList;
  secs->SetOwner();
  secs->SetName("secondaryMaps");
  fList->Add(secs);
  TH2D* hCoverage = new TH2D("coverage", "#eta coverage per v_{z}", 
			     secMap->GetCorrection(UShort_t(1))->GetXaxis()->GetNbins(),
			     secMap->GetCorrection(UShort_t(1))->GetXaxis()->GetXmin(),
			     secMap->GetCorrection(UShort_t(1))->GetXaxis()->GetXmax(),
			     vaxis.GetNbins(),vaxis.GetXmin(),vaxis.GetXmax());
  hCoverage->SetDirectory(0);
  hCoverage->SetXTitle("#eta");
  hCoverage->SetYTitle("v_{z} [cm]");
  hCoverage->SetZTitle("n_{bins}");
  fList->Add(hCoverage);
  
  for (Int_t v = 1; v <= vaxis.GetNbins(); v++) { 
    TH2D* corr = secMap->GetCorrection(UShort_t(v));
    TH1D* proj = corr->ProjectionX(Form("secCor%02d", v));
    proj->Scale(1. / corr->GetNbinsY());
    proj->SetTitle(Form("Projection of secondary correction "
			"for %+5.1f<v_{z}<%+5.1f",
			vaxis.GetBinLowEdge(v), vaxis.GetBinUpEdge(v)));
    proj->SetYTitle("#LT 2^{nd} correction#GT");
    proj->SetDirectory(0);
    proj->SetMarkerStyle(20);
    proj->SetMarkerColor(kBlue+1);
    secs->Add(proj);
    
    TH2D* obg = static_cast<TH2D*>(corr->Clone(Form("secCor2DFiducial%02d",v)));
    obg->SetDirectory(0);
    secs->Add(obg);
    
    TH1D* after = static_cast<TH1D*>(proj->Clone(Form("secCorFiducial%02d",v)));
    after->SetDirectory(0);
    after->SetMarkerColor(kRed+1);
    secs->Add(after);
    
    TH2D* data = static_cast<TH2D*>(corr->Clone(Form("hitMap%02d",v)));
    //d->SetTitle(Form("hitMap%02d",v));
    data->SetTitle(Form("d^{2}N/d#eta d#phi "
			"for %+5.1f<v_{z}<%+5.1f",
			vaxis.GetBinLowEdge(v), vaxis.GetBinUpEdge(v)));
    data->GetZaxis()->SetTitle("");
    data->SetMarkerColor(kBlack);
    data->SetMarkerStyle(1);
    hits->Add(data);
    
    TH1D* hAcceptance = fManager.GetAcceptanceCorrection(v);
    TH1D* accClone   = static_cast<TH1D*>(hAcceptance->Clone(Form("acceptance%02d",v)));
    secs->Add(accClone);
    
    // Double_t prev = 0;
    for (Int_t e = 1; e <= proj->GetNbinsX(); e++) { 
      Double_t c = proj->GetBinContent(e);
      if (c > .5 /*&& TMath::Abs(c - prev) < .1*c*/) {
	fEtaMin[v-1] = e;
	break;
      }
      // prev = c;
      after->SetBinContent(e, 0);
      after->SetBinError(e, 0);
      for(Int_t nn =1; nn <=obg->GetNbinsY();nn++)
	obg->SetBinContent(e,nn,0);
      
    }
    for (Int_t e = proj->GetNbinsX(); e >= 1; e--) { 
      Double_t c = proj->GetBinContent(e);
      if (c > .5 /*&& TMath::Abs(c - prev) < .1*c*/) {
	fEtaMax[v-1] = e;
	break;
      }
      // prev = c;
      after->SetBinContent(e, 0);
      after->SetBinError(e, 0);
      for(Int_t nn =1; nn <=obg->GetNbinsY();nn++)
	obg->SetBinContent(e,nn,0);
      
    }
    
    for (Int_t nn = fEtaMin[v-1]; nn<=fEtaMax[v-1]; nn++) { 
      hCoverage->SetBinContent(nn,v,1);
    }
    
  }
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
  fIvz = 0;

  AliESDEvent* esd = GetESDEvent();
  
  Bool_t   lowFlux   = kFALSE;
  UInt_t   triggers  = 0;
  UShort_t ivz       = 0;
  Double_t vz        = 0;
  Double_t cent      = -1;
  UShort_t nClusters = 0;
  UInt_t   found     = fInspector.Process(esd, triggers, lowFlux, 
					  ivz, vz, cent, nClusters);

  // No event or no trigger 
  if (found & AliFMDEventInspector::kNoEvent)    return;
  if (found & AliFMDEventInspector::kNoTriggers) return;
  
  // Make sure AOD is filled
  MarkEventForStore();

  if (found == AliFMDEventInspector::kNoSPD)      return;
  if (found == AliFMDEventInspector::kNoVertex)   return;
  if (triggers & AliAODForwardMult::kPileUp)      return;
  if (found == AliFMDEventInspector::kBadVertex)  return; // Out of range
  
  //Doing analysis
  fIvz = ivz;
  const AliMultiplicity* spdmult = esd->GetMultiplicity();

  TH2D& aodHist = fAODCentral.GetHistogram();

  ProcessESD(aodHist, spdmult);
  CorrectData(aodHist, ivz);
  //Producing hit maps
  TList* hitList = static_cast<TList*>(fList->FindObject("hitMaps"));
  TH2D* data  = 0;
  if(hitList)
    data = static_cast<TH2D*>(hitList->At(ivz-1));
  if(data)
    data->Add(&aodHist);
  
  PostData(1,fList);
}
//____________________________________________________________________
void 
AliCentralMultiplicityTask::ProcessESD(TH2D& aodHist, 
				       const AliMultiplicity* spdmult) const
{
  fNTracklet->Reset();
  fNCluster->Reset();

  //Filling clusters in layer 1 used for tracklets...
  for(Int_t j = 0; j< spdmult->GetNumberOfTracklets();j++) {
    Double_t eta = spdmult->GetEta(j);
    fNTracklet->Fill(eta);
    aodHist.Fill(eta,spdmult->GetPhi(j));
  }

  //...and then the unused ones in layer 1 
  for(Int_t j = 0; j< spdmult->GetNumberOfSingleClusters();j++) {
    Double_t eta = -TMath::Log(TMath::Tan(spdmult->GetThetaSingle(j)/2.));
    fNCluster->Fill(eta);
    aodHist.Fill(eta, spdmult->GetPhiSingle(j));
  }
  fNClusterTracklet->Fill(fNCluster->GetEntries(), 
			  fNTracklet->GetEntries());
  
  fNCluster->Divide(fNTracklet);
  for (Int_t j = 1; j <= fNCluster->GetNbinsX(); j++)  
    fClusterPerTracklet->Fill(fNCluster->GetXaxis()->GetBinCenter(j), 
			      fNCluster->GetBinContent(j));

}

//____________________________________________________________________
void 
AliCentralMultiplicityTask::CorrectData(TH2D& aodHist, UShort_t vtxbin) const
{  
  // Corrections
  TH1D* hAcceptance = fManager.GetAcceptanceCorrection(vtxbin);
  TH2D* hSecMap     = fManager.GetSecMapCorrection(vtxbin);
  
  if (!hSecMap)     AliFatal("No secondary map!");
  if (!hAcceptance) AliFatal("No acceptance!");
    
  if (fUseSecondary && hSecMap) aodHist.Divide(hSecMap);
  
  for(Int_t nx = 1; nx <= aodHist.GetNbinsX(); nx++) {
    Float_t accCor = hAcceptance->GetBinContent(nx);
    Float_t accErr = hAcceptance->GetBinError(nx);

    Bool_t fiducial = true;
    if (nx < fEtaMin[vtxbin-1] || nx > fEtaMax[vtxbin-1]) 
      fiducial = false;
    //  Bool_t etabinSeen = kFALSE;  
    for(Int_t ny = 1; ny <= aodHist.GetNbinsY(); ny++) {
#if 1
      if (!fiducial) { 
	aodHist.SetBinContent(nx, ny, 0);
	aodHist.SetBinError(nx, ny, 0);
	continue;
      }
#endif	
      // Get currrent value 
      Float_t aodValue = aodHist.GetBinContent(nx,ny);
      Float_t aodErr   = aodHist.GetBinError(nx,ny);

#if 0 // This is done once in the FindEtaBins function
      // Set underflow bin
      Float_t secCor   = 0;
      if(hSecMap)       secCor     = hSecMap->GetBinContent(nx,ny);
      if (secCor > 0.5) etabinSeen = kTRUE;
#endif
      if (aodValue < 0.000001) { 
	aodHist.SetBinContent(nx,ny, 0); 
	aodHist.SetBinError(nx,ny, 0); 
	continue; 
      }
      if (!fUseAcceptance) continue; 

      // Acceptance correction 
      if (accCor   < 0.000001) accCor = 1;
      Float_t aodNew   = aodValue / accCor ;
      Float_t error    = aodNew*TMath::Sqrt(TMath::Power(aodErr/aodValue,2) +
					    TMath::Power(accErr/accCor,2) );
      aodHist.SetBinContent(nx,ny, aodNew);
      //test
      aodHist.SetBinError(nx,ny,error);
      aodHist.SetBinError(nx,ny,aodErr);
    }
    //Filling underflow bin if we eta bin is in range
    if (fiducial) aodHist.SetBinContent(nx,0, 1.);
    // if (etabinSeen) aodHist.SetBinContent(nx,0, 1.);
  }  
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
AliCentralMultiplicityTask::Print(Option_t* option) const
{
  // 
  // Print information 
  // 
  // Parameters:
  //    option Not used
  //
  std::cout << ClassName() << ": " << GetName() << "\n" 
	    << std::boolalpha
	    << "  Use secondary correction:  " << fUseSecondary << '\n'
	    << "  Use acceptance correction: " << fUseAcceptance << '\n' 
	    << "  Off-line trigger mask:  0x" 
	    << std::hex     << std::setfill('0') 
	    << std::setw (8) << fOfflineTriggerMask 
	    << std::dec     << std::setfill (' ') 
	    << std::noboolalpha << std::endl;
  AliCentralCorrSecondaryMap* secMap = GetManager().GetSecMap();
  if (secMap) {
    const TAxis& vaxis = secMap->GetVertexAxis();
    std::cout << "  Eta ranges:\n"
	    << "     Vertex        | Eta bins\n"
	      << "   bin     range   | \n"
	      << "   ----------------+-----------" << std::endl;
    for (Int_t v = 1; v <= vaxis.GetNbins(); v++) { 
      std::cout << "   " << std::setw(2) << v << "  " 
		<< std::setw(5) << vaxis.GetBinLowEdge(v) << "-"
		<< std::setw(5) << vaxis.GetBinUpEdge(v) << " | "
		<< std::setw(3) << fEtaMin[v-1] << "-" 
		<< std::setw(3) << fEtaMax[v-1] << std::endl;
    }
  }

  gROOT->IncreaseDirLevel();
  fManager.Print(option);
  fInspector.Print(option);
  gROOT->DecreaseDirLevel();
  
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
  
  switch(what) {
  case 0:  fname = fSecMapName;     break;
  case 1:  fname = fAcceptanceName; break;
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
	   "Central Manager initialised for %s, energy %dGeV, field %dkG",
	   sys == 1 ? "pp" : sys == 2 ? "PbPb" : "unknown", sNN,field);
  }  
}
//____________________________________________________________________
Bool_t
AliCentralMultiplicityTask::Manager::WriteFile(UShort_t what, 
					       UShort_t sys, 
					       UShort_t sNN, 
					       Short_t  fld, 
					       TObject* obj, 
					       Bool_t   full) const
{
  // 
  // Write correction output to (a temporary) file 
  // 
  // Parameters: 
  //   What     What to write 
  //   sys      Collision system (1: pp, 2: PbPb)
  //   sNN      Center of mass energy per nucleon (GeV)
  //   fld      Field (kG)
  //   obj      Object to write 
  //   full     if true, write to full path, otherwise locally
  // 
  // Return: 
  //   true on success. 
  TString ofName;
  if (!full)
    ofName = GetFileName(what, sys, sNN, fld);
  else 
    ofName = GetFullFileName(what, sys, sNN, fld);
  if (ofName.IsNull()) { 
    AliErrorGeneral("Manager",Form("Unknown object type %d", what));
    return false;
  }
  TFile* output = TFile::Open(ofName, "RECREATE");
  if (!output) { 
    AliErrorGeneral("Manager",Form("Failed to open file %s", ofName.Data()));
    return false;
  }
  
  TString oName(GetObjectName(what));
  Int_t ret = obj->Write(oName);
  if (ret <= 0) { 
    AliErrorGeneral("Manager",Form("Failed to write %p to %s/%s (%d)", 
				   obj, ofName.Data(), oName.Data(), ret));
    return false;
  }

  ret = output->Write();
  if (ret < 0) { 
    AliErrorGeneral("Manager",
		    Form("Failed to write %s to disk (%d)", ofName.Data(),ret));
    return false;
  }
  output->ls();
  output->Close();
  
  TString cName(obj->IsA()->GetName());
  AliInfoGeneral("Manager",
		 Form("Wrote %s object %s to %s\n",
		      cName.Data(),oName.Data(), ofName.Data()));
  if (!full) { 
    TString dName(GetFileDir(what));
    AliInfoGeneral("Manager",
		   Form("%s should be copied to %s\n"
			"Do for example\n\t"
			"aliroot $ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/"
			"MoveCorrections.C\\(%d\\)\nor\n\t"
			"cp %s %s/", 
			ofName.Data(),dName.Data(), 
			what, ofName.Data(), 
			gSystem->ExpandPathName(dName.Data())));


  }
  return true;
}

//____________________________________________________________________
void 
AliCentralMultiplicityTask::Manager::Print(Option_t* option) const
{
  // 
  // Print information to standard output 
  //
  std::cout << " AliCentralMultiplicityTask::Manager\n" 
	    << std::boolalpha 
	    << "  Initialized:     " << fIsInit << '\n'
	    << "  Acceptance path: " << fAcceptancePath << '\n'
	    << "  Acceptance name: " << fAcceptanceName << '\n'
	    << "  Acceptance:      " << fAcceptance << '\n'
	    << "  Secondary path:  " << fSecMapPath << '\n'
	    << "  Secondary name:  " << fSecMapName << '\n'
	    << "  Secondary map:   " << fSecmap 
	    << std::noboolalpha << std::endl;
  if (fAcceptance) fAcceptance->Print(option);
  if (fSecmap)     fSecmap->Print(option);
}

//
// EOF
//
