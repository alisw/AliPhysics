#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliESDVZERO.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliAnalysisTaskHFEIPDistribution.h"
#include "AliHFEpidTOF.h"
#include "AliHFEpid.h"
#include "AliHFEtools.h"
//#include "AliESDTrdTracklet.h"
//#include "AliTRDgeometry.h"
//#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"
#include "AliHFEVZEROEventPlane.h"
#include "AliEventplane.h"
#include <AliOADBContainer.h>


ClassImp(AliAnalysisTaskHFEIPDistribution)

//________________________________________________________________________
AliAnalysisTaskHFEIPDistribution::AliAnalysisTaskHFEIPDistribution()
  : AliAnalysisTaskSE(), fAOD(0), fOutputContainer(0), EP2040(0), EP2040Corrected(0), EP2040V0A(0), EP2040V0C(0), TPCnSigma(0), EPCent(0), EPCentCorrected(0), EPCentV0A(0), EPCentV0C(0), DeltaPhi(0), fExtraCuts(0), fIPData(0), fpTIP2040IP(0), fpTIP2040OOP(0), fpTIP3050IP(0), fpTIP3050OOP(0), EventSelectionSteps(0)
{
  // default Constructor
  // Define input and output slots here
  
  // HFE cuts
    /*hfetrackCuts = new AliHFEcuts("V0trackCuts", "Track Cuts for tagged track Analysis");
    hfetrackCuts->CreateStandardCuts();

    hfetrackCuts->SetMinNClustersTPC(110);
    hfetrackCuts->SetMinNClustersTPCPID(80);
    hfetrackCuts->SetFractionOfSharedTPCClusters(1.1);
    hfetrackCuts->SetMinRatioTPCclusters(0.6);
    hfetrackCuts->SetTPCmodes(3,4);
    hfetrackCuts->SetMinNClustersITS(4);
    hfetrackCuts->SetCutITSpixel(AliHFEextraCuts::kBoth);
    hfetrackCuts->SetCheckITSLayerStatus(kFALSE);
    //hfetrackCuts->UnsetVertexRequirement();
    hfetrackCuts->SetTOFPIDStep(kTRUE);
    hfetrackCuts->Initialize();*/
    
    //fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");

    // end HFE cuts
  
}

//________________________________________________________________________
AliAnalysisTaskHFEIPDistribution::AliAnalysisTaskHFEIPDistribution(const char *name)
  : AliAnalysisTaskSE(name), fAOD(0), fOutputContainer(0),
  EP2040(0), EP2040Corrected(0), EP2040V0A(0), EP2040V0C(0), TPCnSigma(0), EPCent(0), EPCentCorrected(0), EPCentV0A(0), EPCentV0C(0), DeltaPhi(0), fExtraCuts(0), fIPData(0), fpTIP2040IP(0), fpTIP2040OOP(0), fpTIP3050IP(0), fpTIP3050OOP(0), EventSelectionSteps(0)
{
  // HFE cuts
    /*hfetrackCuts = new AliHFEcuts("V0trackCuts", "Track Cuts for tagged track Analysis");
    hfetrackCuts->CreateStandardCuts();

    hfetrackCuts->SetMinNClustersTPC(110);
    hfetrackCuts->SetMinNClustersTPCPID(80);
    hfetrackCuts->SetFractionOfSharedTPCClusters(1.1);
    hfetrackCuts->SetMinRatioTPCclusters(0.6);
    hfetrackCuts->SetTPCmodes(3,4);
    hfetrackCuts->SetMinNClustersITS(4);
    hfetrackCuts->SetCutITSpixel(AliHFEextraCuts::kBoth);
    hfetrackCuts->SetCheckITSLayerStatus(kFALSE);
    //hfetrackCuts->UnsetVertexRequirement();
    hfetrackCuts->SetTOFPIDStep(kTRUE);
    hfetrackCuts->Initialize();
    
    //fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
*/
    // end HFE cuts
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
   DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  //DefineOutput(1, TObjArray::Class());
    DefineOutput(1, TList::Class());
    //    DefineOutput(2, TNtuple::Class());
    

}

AliAnalysisTaskHFEIPDistribution::~AliAnalysisTaskHFEIPDistribution()
{
}


//________________________________________________________________________
void AliAnalysisTaskHFEIPDistribution::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
    Double_t ptbinningX[19] = {0., 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 2., 2.5, 3., 4., 6., 8., 10., 12., 16., 20.}; // Apr 2018 binning
    
    EPCent = new TH2D(Form("EPCent"),Form("EPCent"), 20, 0., TMath::Pi(), 10, 0., 100.);
    EPCentV0A = new TH2D(Form("EPCentV0A"),Form("EPCentV0A"), 20, 0., TMath::Pi(), 10, 0., 100.);
    EPCentV0C = new TH2D(Form("EPCentV0C"),Form("EPCentV0C"), 20, 0., TMath::Pi(), 10, 0., 100.);
    EPCentCorrected = new TH2D(Form("EPCentCorrected"),Form("EPCentCorrected"), 20, 0., TMath::Pi(), 10, 0., 100.);
    
    EP2040 = new TH1D(Form("EP2040"),Form("EP2040"), 100, 0., TMath::Pi());
    EP2040Corrected = new TH1D(Form("EP2040Corrected"),Form("EP2040Corrected"), 100, 0., TMath::Pi());
    EP2040V0A = new TH1D(Form("EP2040V0A"),Form("EP2040V0A"), 100, 0., TMath::Pi());
    EP2040V0C = new TH1D(Form("EP2040V0C"),Form("EP2040V0C"), 100, 0., TMath::Pi());
    DeltaPhi = new TH1D(Form("DeltaPhi"),Form("DeltaPhi"), 40, 0., TMath::Pi());
    
    TPCnSigma = new TH2D("TPCnSigma", "", 18, ptbinningX, 100, -10., 5.);
    fIPData =  new TH2D("fIPData", "", 18, ptbinningX, 400, -0.2, 0.2);

    fpTIP2040IP =  new TH2D("pTIP2040IP", "", 18, ptbinningX, 400, -0.2, 0.2);
    fpTIP2040OOP =  new TH2D("pTIP2040OOP", "", 18, ptbinningX, 400, -0.2, 0.2);
    fpTIP3050IP =  new TH2D("pTIP3050IP", "", 18, ptbinningX, 400, -0.2, 0.2);
    fpTIP3050OOP =  new TH2D("pTIP3050OOP", "", 18, ptbinningX, 400, -0.2, 0.2);

    EventSelectionSteps = new TH1D(Form("EventSelectionSteps"),Form("EventSelectionSteps"), 10, -0.5, 9.5);

    
    fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
    
    fOutputContainer = new TObjArray(1);
    fOutputContainer->SetName(GetName());
    fOutputContainer->SetOwner();
    fOutputContainer->Add(EPCent);
    fOutputContainer->Add(EPCentV0A);
    fOutputContainer->Add(EPCentV0C);
    fOutputContainer->Add(EPCentCorrected);
    fOutputContainer->Add(EP2040);
    fOutputContainer->Add(EP2040Corrected);
    fOutputContainer->Add(EP2040V0A);
    fOutputContainer->Add(EP2040V0C);
    fOutputContainer->Add(DeltaPhi);
    fOutputContainer->Add(TPCnSigma);
    fOutputContainer->Add(fIPData);
    fOutputContainer->Add(fpTIP2040IP);
    fOutputContainer->Add(fpTIP2040OOP);
    fOutputContainer->Add(fpTIP3050IP);
    fOutputContainer->Add(fpTIP3050OOP);
    fOutputContainer->Add(EventSelectionSteps);

    PostData(1, fOutputContainer);
    
}

//_____________________________________________________________________________
void AliAnalysisTaskHFEIPDistribution::UserExec(Option_t *)
{
  //
  // Called for each event
    //
    

    //AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    //AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());

    if (!fAOD)
    {
	printf("ERROR: Could not get AODInputHandler \n");
    }
    else fAOD = dynamic_cast<AliAODEvent*>(InputEvent());


    Process(fAOD);


    // Post output data.
    //PostData(1, fOutput);
    //PostData(1, fOutputContainer);
    //    PostData(2, rec);
    //    PostData(2, treeoutput);
}

Bool_t AliAnalysisTaskHFEIPDistribution::PassesTrackCuts(AliAODTrack *track)
{
    if(TMath::Abs(track->Eta())>0.8) return kFALSE;
    ULong_t status = track->GetStatus();
    // Basic tracking
    if(!((status & AliVTrack::kITSrefit) && (status & AliVTrack::kTPCrefit))) return kFALSE;
    // ITS tracking
    if(!(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1) && track->GetITSNcls()>=4)) return kFALSE;
    Int_t nclustersITS(track->GetITSclusters(NULL)),
          nclustersTPC(track->GetTPCNcls()),
          nclustersTPCall(track->GetTPCClusterMap().CountBits()),
          nclustersTPCshared(0);
    UChar_t nfindableTPC = track->GetTPCNclsF();
    Double_t FoundOverFindable = (nfindableTPC ? static_cast<Float_t>(nclustersTPC)/static_cast<Float_t>(nfindableTPC) : 0);
    // TPC Quality
    if(!(track->GetTPCNcls()>=110 && track->GetTPCsignalN()>=80 && FoundOverFindable>0.6 && track->GetTPCchi2()/Double_t(nclustersTPC)<4.)) return kFALSE;
    // ITS Quality
    if(track->GetITSchi2()/double(track->GetITSNcls()) > 10.) return kFALSE; 
      UChar_t ITSShared = track->GetITSSharedClusterMap();
      Int_t nSharedITS = 0;
      for(int i=0;i<6;i++)
          if((ITSShared >> i) & 1)
              nSharedITS++;
    if(nSharedITS > 3) return kFALSE; 
    // TOF quality
    if(!((status & AliVTrack::kTOFpid))) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTaskHFEIPDistribution::PassesElectronPID(AliAODTrack *track, AliPIDResponse *pid)
{
    if(pid->NumberOfSigmasTPC(track, AliPID::kElectron) >3. || pid->NumberOfSigmasTPC(track, AliPID::kElectron) <-0.5) return kFALSE;
    if(TMath::Abs(pid->NumberOfSigmasTOF(track, AliPID::kElectron)>3.))return kFALSE;
    return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskHFEIPDistribution::Process(AliAODEvent *const aodEvent)
{
  // Main loop
  // Called for each event
  EventSelectionSteps->Fill(4);
  if(!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB))) return;
  EventSelectionSteps->Fill(5);

  if (!aodEvent) {
    Printf("ERROR: aodEvent not available");
    //return;
  }
  AliInputEventHandler* handler = dynamic_cast<AliInputEventHandler*>( AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!(handler)) printf("AOD inputhandler not available \n");


  AliPIDResponse *pid = NULL;
  if(handler){
    pid = handler->GetPIDResponse();
  } else {
    AliError("No Handler");
  }
  if(!pid){
    AliError("No PID response");
    return;
  }
  
    if(!fExtraCuts){
    fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
  }
  fExtraCuts->SetRecEventInfo(aodEvent);

  /*AliMCEventHandler* MCEventHandler = (AliMCEventHandler*)AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
  //AliMCEvent* fMCEvent=0;
  if(MCEventHandler)
  fMCEvent = MCEventHandler->MCEvent();*/

  Float_t centrality = -1.;
  AliCentrality *hicent = aodEvent->GetCentrality();
  centrality = hicent->GetCentralityPercentile("V0M");
  
  const AliAODVertex *vertex = aodEvent->GetPrimaryVertex();
  const AliAODVertex *vertexSPD = aodEvent->GetPrimaryVertexSPD();
  Double_t vcov[6];
  vertex->GetCovMatrix(vcov);
  Double_t vtx[3];
  Double_t vtxSPD[3];
  vertex->GetXYZ(vtx);
  vertexSPD->GetXYZ(vtxSPD);

  bool analyzeEvent=(centrality>=20.0 && centrality<=50.0 && TMath::Sqrt(vcov[5]) < 0.25 && TMath::Abs(vtx[2])<10. && TMath::Abs(vtx[2] - vtxSPD[2]) < 0.5);
  //hfetrackCuts->SetRecEvent(aodEvent);

   if(centrality>=20.0 && centrality<=40.0)
   {
      EventSelectionSteps->Fill(0);
      if(TMath::Abs(vtx[2])<10.)
      {
        EventSelectionSteps->Fill(1);
        if(TMath::Abs(vtx[2] - vtxSPD[2]) < 0.5)
        {
          EventSelectionSteps->Fill(2);
          if(TMath::Sqrt(vcov[5])) EventSelectionSteps->Fill(3);
        }
      }
   }

  ULong_t status;
  AliEventplane* vEPa;
  Double_t qVx, qVy;
  double epCorrArray[2];
  double epCorr;
  if(analyzeEvent)
  {

      vEPa = fInputEvent->GetEventplane();
      Float_t V0PlanePhi = TVector2::Phi_0_2pi(vEPa->CalculateVZEROEventPlane(fInputEvent,10,2,qVx,qVy));
      if(V0PlanePhi>TMath::Pi()) V0PlanePhi -=TMath::Pi();
      Float_t V0APlanePhi = TVector2::Phi_0_2pi(vEPa->CalculateVZEROEventPlane(fInputEvent,8,2,qVx,qVy));
      if(V0APlanePhi > TMath::Pi()) V0APlanePhi = V0APlanePhi - TMath::Pi();
      Float_t V0CPlanePhi = TVector2::Phi_0_2pi(vEPa->CalculateVZEROEventPlane(fInputEvent,9,2,qVx,qVy));
      if(V0CPlanePhi > TMath::Pi()) V0CPlanePhi = V0CPlanePhi - TMath::Pi();
      
      
      EPCent->Fill(V0PlanePhi, centrality);
      EPCentV0A->Fill(V0APlanePhi, centrality);
      EPCentV0C->Fill(V0CPlanePhi, centrality);
      //EPCentCorrected->Fill(epCorr, centrality);
      EP2040->Fill(V0PlanePhi);
      //EP2040Corrected->Fill(epCorr);
      EP2040V0A->Fill(V0APlanePhi);
      EP2040V0C->Fill(V0CPlanePhi);
      
      double DPhi =0.;
      Float_t dcaxy = -999., dcaz = -999.;
      Double_t dcaErr, dcaxyD;
      AliAODTrack *track = 0x0;
      for(Int_t itrack = 0; itrack < fInputEvent->GetNumberOfTracks(); itrack++)
      {
          // Run track loop
          track = dynamic_cast<AliAODTrack *>(fInputEvent->GetTrack(itrack));
          if(!track) continue;
          if(!PassesTrackCuts(track)) continue;
          if(!PassesElectronPID(track, pid)) continue;
          TPCnSigma->Fill(track->Pt(),pid->NumberOfSigmasTPC(track, AliPID::kElectron));
          DPhi = track->Phi() - V0PlanePhi - TMath::Pi();
          while(DPhi<0.) DPhi += TMath::Pi();
          
          fExtraCuts->GetImpactParameters((AliVTrack *)track,dcaxy,dcaz);
          fExtraCuts->GetHFEImpactParameters((AliVTrack *)track,dcaxyD,dcaErr);
          fIPData->Fill(track->Pt(), dcaxyD);
          if(track->Pt() > 0.5 && centrality>=20.0 && centrality<=50.0)
          {
            DeltaPhi->Fill(DPhi);
            if(DPhi/2. < TMath::Pi()/4.) // IP
            {
              if(centrality>=20.0 && centrality<=40.0) fpTIP2040IP->Fill(track->Pt(), dcaxyD*track->Charge()*TMath::Sign(1.,aodEvent->GetMagneticField()));
              if(centrality>=30.0 && centrality<=50.0) fpTIP3050IP->Fill(track->Pt(), dcaxyD*track->Charge()*TMath::Sign(1.,aodEvent->GetMagneticField()));
            }
            else
            {
              if(centrality>=20.0 && centrality<=40.0) fpTIP2040OOP->Fill(track->Pt(), dcaxyD*track->Charge()*TMath::Sign(1.,aodEvent->GetMagneticField()));
              if(centrality>=30.0 && centrality<=50.0) fpTIP3050OOP->Fill(track->Pt(), dcaxyD*track->Charge()*TMath::Sign(1.,aodEvent->GetMagneticField()));
            }
          }

          
      }
      
  }
  

  PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskHFEIPDistribution::Terminate(const Option_t *)
{
//   treeoutput->Write();

}
