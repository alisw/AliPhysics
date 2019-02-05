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
#include "AliAODv0KineCuts.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliESDVZERO.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliAnalysisTaskHFEIPCorrection.h"
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
#include "AliMultSelection.h"
#include "AliVertexerTracks.h"

// This analysis builds histograms to enable and check the impact parameter
// correction in phi, z, and pT for the 15o PbPb data set
// Author: Martin Voelkl


ClassImp(AliAnalysisTaskHFEIPCorrection)

//________________________________________________________________________
AliAnalysisTaskHFEIPCorrection::AliAnalysisTaskHFEIPCorrection()
  : AliAnalysisTaskSE(), fAOD(0), fOutputContainer(0), EP2040(0), EP2040Corrected(0), EP2040V0A(0), EP2040V0C(0), TPCnSigma(0), EPCent(0), EPCentCorrected(0), EPCentV0A(0), EPCentV0C(0), DeltaPhi(0), fExtraCuts(0), fIPData(0), fpTIP2040IP(0), fpTIP2040OOP(0), fpTIP3050IP(0), fpTIP3050OOP(0), EventSelectionSteps(0), fPionV0pTRNoCuts(0), fPionV0pTRWithCuts(0), fAODV0Cuts(0), fPionV0pTTPC(0), fPionV0pTTPCWithCuts(0), fDCARegionRun(0), fDCAPhiZHadrons(0), fDCAPhiZHadronsC(0), fDCAPhiZHadronsC2(0), fDCAPhiZHadrons2nd(0), fpTPhiZHadrons(0), fDCAPhipTHadrons(0), fDCAPhipTHadronsC(0), fDCAPhipTHadronsC2(0), fDCAWErrHadrons(0), fDCAHadronsFineBins(0), fDCAKaons(0), fDCAWErrKaons(0), fDCAKaonsFineBins(0), fDCAvsCorrected(0)
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
AliAnalysisTaskHFEIPCorrection::AliAnalysisTaskHFEIPCorrection(const char *name)
  : AliAnalysisTaskSE(name), fAOD(0), fOutputContainer(0), EP2040(0), EP2040Corrected(0), EP2040V0A(0), EP2040V0C(0), TPCnSigma(0), EPCent(0), EPCentCorrected(0), EPCentV0A(0), EPCentV0C(0), DeltaPhi(0), fExtraCuts(0), fIPData(0), fpTIP2040IP(0), fpTIP2040OOP(0), fpTIP3050IP(0), fpTIP3050OOP(0), EventSelectionSteps(0), fPionV0pTRNoCuts(0), fPionV0pTRWithCuts(0), fAODV0Cuts(0), fPionV0pTTPC(0), fPionV0pTTPCWithCuts(0), fDCARegionRun(0), fDCAPhiZHadrons(0), fDCAPhiZHadronsC(0), fDCAPhiZHadronsC2(0), fDCAPhiZHadrons2nd(0), fpTPhiZHadrons(0), fDCAPhipTHadrons(0), fDCAPhipTHadronsC(0), fDCAPhipTHadronsC2(0), fDCAWErrHadrons(0), fDCAHadronsFineBins(0), fDCAKaons(0), fDCAWErrKaons(0), fDCAKaonsFineBins(0), fDCAvsCorrected(0)
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

AliAnalysisTaskHFEIPCorrection::~AliAnalysisTaskHFEIPCorrection()
{
}


//________________________________________________________________________
void AliAnalysisTaskHFEIPCorrection::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
    Double_t ptbinningX[19] = {0., 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 2., 2.5, 3., 4., 6., 8., 10., 12., 16., 20.}; // Apr 2018 binning
    //Double_t IPBinning[401];
    
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

    fPionV0pTRNoCuts = new TH3D(Form("fPionV0pTRNoCuts"),Form("fPionV0pTRNoCuts"), 40, 0., 10., 80, 0., 20., 10, 0., 100.);
    fPionV0pTRWithCuts = new TH3D(Form("fPionV0pTRWithCuts"),Form("fPionV0pTRWithCuts"), 40, 0., 10., 80, 0., 20., 10, 0., 100.);
    fPionV0pTTPC = new TH2D(Form("fPionV0pTTPC"),Form("fPionV0pTTPC"), 18, ptbinningX, 200, -10., 10.);
    fPionV0pTTPCWithCuts = new TH2D(Form("fPionV0pTTPCWithCuts"),Form("fPionV0pTTPCWithCuts"), 18, ptbinningX, 200, -10., 10.);

    EventSelectionSteps = new TH1D(Form("EventSelectionSteps"),Form("EventSelectionSteps"), 10, -0.5, 9.5);

    fDCARegionRun = new TH3D(Form("fDCARegionRun"),Form("fDCARegionRun"), 400, -0.05, 0.05, 5, 0.5, 5.5, 183, 0.5, 183.5);
    fDCAPhiZHadrons = new TH3D(Form("fDCAPhiZHadrons"),Form("fDCAPhiZHadrons"), 40, 0., 2.*3.14159, 400, -0.05, 0.05, 12, -12., 12.);
    fDCAPhiZHadronsC = new TH3D(Form("fDCAPhiZHadronsC"),Form("fDCAPhiZHadronsC"), 40, 0., 2.*3.14159, 400, -0.05, 0.05, 12, -12., 12.);
    fDCAPhiZHadronsC2 = new TH3D(Form("fDCAPhiZHadronsC2"),Form("fDCAPhiZHadronsC2"), 40, 0., 2.*3.14159, 400, -0.05, 0.05, 12, -12., 12.);
    fDCAPhiZHadrons2nd = new TH3D(Form("fDCAPhiZHadrons2nd"),Form("fDCAPhiZHadrons2nd"), 40, 0., 2.*3.14159, 400, -0.05, 0.05, 12, -12., 12.);
    fpTPhiZHadrons = new TH3D(Form("fpTPhiZHadrons"),Form("fpTPhiZHadrons"), 40, 0., 2.*3.14159, 40, 0., 10., 12, -12., 12.);
    fDCAPhipTHadrons = new TH3D(Form("fDCAPhipTHadrons"),Form("fDCAPhipTHadrons"), 400, -0.05, 0.05, 40, 0., 2.*3.14159, 40, 0., 10.);
    fDCAPhipTHadronsC = new TH3D(Form("fDCAPhipTHadronsC"),Form("fDCAPhipTHadronsC"), 400, -0.05, 0.05, 40, 0., 2.*3.14159, 40, 0., 10.);
    fDCAPhipTHadronsC2 = new TH3D(Form("fDCAPhipTHadronsC2"),Form("fDCAPhipTHadronsC2"), 400, -0.05, 0.05, 40, 0., 2.*3.14159, 40, 0., 10.);
    fDCAWErrHadrons = new TH3D(Form("fDCAWErrHadrons"),Form("fDCAWErrHadrons"), 80, 0., 10., 400, -0.2, 0.2, 100, 0.0, 0.01);
    fDCAHadronsFineBins = new TH3D(Form("fDCAHadronsFineBins"),Form("fDCAHadronsFineBins"), 80, 0., 10., 400, -0.2, 0.2, 10, 0., 100.);
    fDCAKaons = new TH2D(Form("fDCAKaons"),Form("fDCAKaons"), 18, ptbinningX, 400, -0.2, 0.2);
    fDCAWErrKaons = new TH3D(Form("fDCAWErrKaons"),Form("fDCAWErrKaons"), 80, 0., 10., 400, -0.2, 0.2, 100, 0.0, 0.01);
    fDCAKaonsFineBins = new TH3D(Form("fDCAKaonsFineBins"),Form("fDCAKaonsFineBins"), 80, 0., 10., 400, -0.2, 0.2, 10, 0., 100.);
    fDCAvsCorrected = new TH2D(Form("fDCAvsCorrected"),Form("fDCAvsCorrected"), 200, -0.03, 0.03, 200, -0.03, 0.03);
    
    fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
    fAODV0Cuts = new AliAODv0KineCuts();
    
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
    fOutputContainer->Add(fPionV0pTRNoCuts);
    fOutputContainer->Add(fPionV0pTRWithCuts);
    fOutputContainer->Add(fPionV0pTTPC);
    fOutputContainer->Add(fPionV0pTTPCWithCuts);
    fOutputContainer->Add(EventSelectionSteps);
    fOutputContainer->Add(fDCARegionRun);
    fOutputContainer->Add(fDCAPhiZHadrons);
    fOutputContainer->Add(fDCAPhiZHadronsC);
    fOutputContainer->Add(fDCAPhiZHadronsC2);
    fOutputContainer->Add(fDCAPhiZHadrons2nd);
    fOutputContainer->Add(fpTPhiZHadrons);
    fOutputContainer->Add(fDCAPhipTHadrons);
    fOutputContainer->Add(fDCAPhipTHadronsC);
    fOutputContainer->Add(fDCAPhipTHadronsC2);
    fOutputContainer->Add(fDCAWErrHadrons);
    fOutputContainer->Add(fDCAHadronsFineBins);
    fOutputContainer->Add(fDCAKaons);
    fOutputContainer->Add(fDCAWErrKaons);
    fOutputContainer->Add(fDCAKaonsFineBins);
    fOutputContainer->Add(fDCAvsCorrected);

    PostData(1, fOutputContainer);
    
}

//_____________________________________________________________________________
void AliAnalysisTaskHFEIPCorrection::UserExec(Option_t *)
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

Bool_t AliAnalysisTaskHFEIPCorrection::PassesTrackCuts(AliAODTrack *track)
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
    //if(track->GetTOFsignal()>=99999) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTaskHFEIPCorrection::PassesElectronPID(AliAODTrack *track, AliPIDResponse *pid)
{
    if(pid->NumberOfSigmasTPC(track, AliPID::kElectron) >3. || pid->NumberOfSigmasTPC(track, AliPID::kElectron) <-0.5) return kFALSE;
    if(TMath::Abs(pid->NumberOfSigmasTOF(track, AliPID::kElectron))>3.) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTaskHFEIPCorrection::PassesPionPID(AliAODTrack *track, AliPIDResponse *pid)
{
    if(pid->NumberOfSigmasTPC(track, AliPID::kPion) > 3. || pid->NumberOfSigmasTPC(track, AliPID::kPion) < -1.) return kFALSE;
    if(TMath::Abs(pid->NumberOfSigmasTOF(track, AliPID::kPion))>3.) return kFALSE; // Should be basically same as electron
    return kTRUE;
}

Bool_t AliAnalysisTaskHFEIPCorrection::PassesKaonPID(AliAODTrack *track, AliPIDResponse *pid)
{
    if(pid->NumberOfSigmasTPC(track, AliPID::kKaon) >3. || pid->NumberOfSigmasTPC(track, AliPID::kKaon) <-3.) return kFALSE;
    if(TMath::Abs(pid->NumberOfSigmasTOF(track, AliPID::kKaon))>2.) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTaskHFEIPCorrection::PassesMinimalTrackCuts(AliAODTrack *track)
{
    if(TMath::Abs(track->Eta())>0.8) return kFALSE;
    ULong_t status = track->GetStatus();
    // Basic tracking
    if(!((status & AliVTrack::kITSrefit) && (status & AliVTrack::kTPCrefit))) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTaskHFEIPCorrection::PassesITSTrackCuts(AliAODTrack *track)
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
    // ITS Quality
    if(track->GetITSchi2()/double(track->GetITSNcls()) > 10.) return kFALSE; 
      UChar_t ITSShared = track->GetITSSharedClusterMap();
      Int_t nSharedITS = 0;
      for(int i=0;i<6;i++)
          if((ITSShared >> i) & 1)
              nSharedITS++;
    if(nSharedITS > 3) return kFALSE; 
    return kTRUE;
}

Int_t AliAnalysisTaskHFEIPCorrection::IsInMisalignedRegion(AliAODTrack *track, double vtxz)
{
  double halfstavephi = 2. * TMath::Pi() / 40.;
  double phi = track->Phi();
  double zSPD = vtxz+4.5*TMath::Tan(TMath::Pi()/2.-2.*TMath::ATan(TMath::Exp(-track->Eta()))); // z position of track at first layer
  if(phi>6.*halfstavephi && phi<=9.*halfstavephi && zSPD>0.) return 1;
  if(phi>17.*halfstavephi && phi<=19.*halfstavephi && zSPD>0.) return 2;
  if(phi>21.*halfstavephi && phi<=23.*halfstavephi) return 3;
  if(phi>25.*halfstavephi && phi<=27.*halfstavephi && zSPD<0.) return 4;
  return 0;
}

AliAODVertex * AliAnalysisTaskHFEIPCorrection::CorrectVertex(AliAODEvent *aodEvent, double vtxz) // Vertex without using excluded regions
{
  AliAODVertex * vtxcorr;
  AliAODTrack *track = 0x0;
  Int_t * skippedTracks = new Int_t[fInputEvent->GetNumberOfTracks()];
  Int_t nskipped = 0;
  // First fill array of tracks to be disregarded in primary vertex estimation
  for(Int_t itrack = 0; itrack < fInputEvent->GetNumberOfTracks(); itrack++)
      {
        // Run track loop
        track = dynamic_cast<AliAODTrack *>(fInputEvent->GetTrack(itrack));
        if(track->Pt()>0.5 && TMath::Abs(track->Eta())<1.)
        {
          if(IsInMisalignedRegion(track, vtxz)>0)
          {
            if(((Int_t) track->GetID())>0)
            {
              skippedTracks[nskipped] = (Int_t) track->GetID();
              nskipped++;
            }
          }
        }
      }
  // Now reconstruct the primary vertex without these tracks
  AliVertexerTracks vertexer(aodEvent->GetMagneticField());
  vertexer.SetITSMode();
  vertexer.SetMinClusters(3);
  vertexer.SetConstraintOff();
  vertexer.SetSkipTracks(nskipped,skippedTracks);
  AliESDVertex *vtxESDNew = vertexer.FindPrimaryVertex(aodEvent);

  // convert to AliAODVertex -- copied from hfe task
  Double_t pos[3],cov[6],chi2perNDF;
  vtxESDNew->GetXYZ(pos); // position
  vtxESDNew->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vtxESDNew->GetChi2toNDF();
  delete vtxESDNew; vtxESDNew=NULL;

  vtxcorr = new AliAODVertex(pos,cov,chi2perNDF); // Makes AOD vertex

  delete[] skippedTracks;
  return vtxcorr;
}

void AliAnalysisTaskHFEIPCorrection::GetTrackImpactParameter(AliAODEvent *aodEvent, AliAODTrack *track, AliAODVertex * pvtx, Double_t &dcaxy)
{
  Double_t dcaD[3]; // 2 should be sufficient, 3 to check if reason for crash
  Double_t covD[3];
  Double_t fMagField;
  const Double_t kBeampiperadius=3.;
  AliAODTrack aodtrack(*track);
  AliExternalTrackParam etp; etp.CopyFromVTrack(&aodtrack);
  fMagField = aodEvent->GetMagneticField();
  etp.PropagateToDCA(pvtx, fMagField, kBeampiperadius, dcaD, covD);
  dcaxy = dcaD[0];
}

void AliAnalysisTaskHFEIPCorrection::GetCorrectedImpactParameter(AliAODTrack *track, Double_t primVertexZ, Double_t &dcaxy)
{
  // primVertexZ is the z Position of the primary vertex
  // Using 40 phi bins, 12 z bins (-12 to +12 cm) and smooth pT correction
  int PhiBin = int(track->Phi()/(2.*TMath::Pi())*40); // phi bin
  double z_SPD1 = primVertexZ+4.5*TMath::Tan(TMath::Pi()/2.-2.*TMath::ATan(TMath::Exp(-track->Eta()))); // z position of track in inner SPD layer
  int zBin = int((z_SPD1+12.)/24.*12.); // z bin
  int pt = track->Pt();
  double oldIP, dcaErr;
  fExtraCuts->GetHFEImpactParameters((AliVTrack *)track,oldIP,dcaErr);
    
  double correctionMatrix[40][12] = {
  { 0.000282 , -0.000227 , 0.000402 , -6.42e-05 , -0.000242 , -0.000356 , -0.000262 , -0.00021 , -0.000715 , -0.000565 , -0.000648 , -0.000761 },
  { 0.000421 , 0.000477 , -0.000102 , -0.000113 , -0.000104 , -0.000293 , -0.000179 , 0.000147 , -1.83e-05 , -5.21e-05 , -0.000121 , -0.000278 },
  { 0.000863 , 0.000673 , 0.000415 , 0.000274 , 0.000308 , -6.39e-05 , 0.000325 , 0.000224 , 0.000298 , 0.00082 , 0.000559 , 0.000102 },
  { 0.000879 , 0.000633 , 0.000367 , 0.000184 , 0.000104 , -2.96e-07 , 0.000388 , 0.000921 , 0.000851 , 0.000652 , 0.000803 , 0.000146 },
  { -0.000574 , -0.000638 , -0.000881 , -0.000655 , -0.000316 , -0.000386 , 0.00055 , 0.00042 , 0.000497 , 0.000866 , 0.000479 , 0.000237 },
  { -0.000856 , -0.000353 , -0.000599 , -0.000615 , -0.000297 , -0.000407 , 0.000165 , 0.000986 , 0.00086 , 0.000697 , 0.00057 , 8.85e-05 },
  { 0.000385 , 0.000558 , 0.000624 , 0.000267 , 0.000137 , -8.53e-05 , -0.00304 , -0.00281 , -0.00237 , -0.00267 , -0.00212 , -0.00201 },
  { 0.000115 , 0.000256 , -0.000154 , 0.000274 , 0.000809 , 3.29e-06 , -0.00459 , -0.00541 , -0.00498 , -0.00471 , -0.00384 , -0.00318 },
  { -0.000358 , -5.75e-05 , -0.0011 , -0.000182 , 0.00044 , -0.000504 , -0.00215 , -0.00243 , -0.00262 , -0.00289 , -0.00203 , -0.00149 },
  { 0 , 0 , 0 , 0 , 0 , 0.00161 , 0.00131 , 0.0014 , 0.00131 , 0.000964 , 0.000752 , 0.000883 },
  { 0 , 0 , 0 , 0 , 0 , 0.000754 , 0.000305 , 0.000778 , 0.000645 , 0.000356 , 0.000827 , 0.000366 },
  { 0 , 0 , 0 , 0 , 0 , 0.000307 , -0.000151 , 0.000491 , 0.000437 , 0.00028 , 0.000803 , 0.000594 },
  { 0 , 0 , 0 , 0 , 0 , -0.000132 , -0.000395 , 0.000656 , 0.000668 , 3.55e-05 , 0.000759 , 0.000789 },
  { 0 , 0 , 0 , 0 , 0 , -0.00119 , -0.00106 , -0.00049 , -0.000364 , -0.000216 , 0.000302 , 0.000324 },
  { 0 , 0.0126 , 0.0115 , 0.0123 , 0.0136 , 0.00968 , -0.00106 , -0.000168 , -0.000155 , -0.000705 , 2.18e-05 , -0.000317 },
  { 0.000674 , 0.000879 , 0.000719 , 0.000805 , 0.00116 , 0.00102 , 0.000206 , 0.000692 , 0.000293 , -2.6e-05 , 0.000829 , 0.000428 },
  { 0.000889 , 0.00113 , 0.000817 , 0.000894 , 0.000999 , 0.000636 , -0.000117 , 0.000528 , 0.000339 , -0.000119 , 0.000381 , 0.000298 },
  { 0.000872 , 0.000899 , 0.000515 , 0.000672 , 0.000904 , 0.000564 , 0.00125 , 0.00283 , 0.0028 , 0.00287 , 0.00268 , 0.00162 },
  { 0.00029 , 0.000905 , -0.000312 , 0.000636 , 0.000714 , 0.000135 , 2.45e-05 , 0.00222 , 0.00226 , 0.00172 , 0.00181 , 0.00161 },
  { 0.000591 , 0.00127 , 0.000893 , 0.000979 , 0.000803 , 0.000464 , 0.000324 , 0.000737 , 0.000505 , 0.000277 , 0.000441 , 0.000405 },
  { 0.000941 , 0.00101 , 0.000617 , 0.00085 , 0.000619 , 0.000551 , 0.000364 , 0.000229 , 0.000328 , -2.29e-05 , 0.000284 , 0.000372 },
  { -0.0017 , -0.00181 , -0.00224 , -0.00248 , -0.00283 , -0.00306 , -0.00364 , -0.00398 , -0.00349 , -0.00282 , -0.0025 , -0.0026 },
  { -0.00175 , -0.00173 , -0.00232 , -0.00253 , -0.00326 , -0.00328 , -0.00404 , -0.00388 , -0.00418 , -0.00391 , -0.00352 , -0.00311 },
  { -0.000551 , -0.0005 , -0.000489 , -0.000495 , -0.00059 , -0.000369 , 6.83e-05 , 0 , 0 , 0 , 0 , 0 },
  { -0.000858 , -0.000756 , -0.0012 , -0.001 , -0.000262 , -0.00121 , -0.00133 , -0.00186 , -0.00222 , -0.00169 , -0.00115 , -0.000857 },
  { -0.00221 , -0.00167 , -0.00176 , -0.00186 , -0.00186 , -0.00198 , -0.00192 , -0.00146 , -0.00157 , -0.00164 , -0.00152 , -0.000787 },
  { -0.00146 , -0.00169 , -0.00248 , -0.0024 , -0.00253 , -0.00249 , -0.00111 , -0.00153 , -0.00175 , -0.00126 , -0.00143 , -0.000942 },
  { 0 , 0 , 0 , 0 , 0 , 0 , -0.00113 , -0.00132 , -0.00143 , -0.00147 , -0.0013 , -0.00125 },
  { 0.0016 , 0.000813 , 0.00113 , 0.00141 , 0.000816 , 0.000618 , 0.000258 , -0.000228 , -0.000502 , 3.33e-05 , -0.00053 , -0.000589 },
  { -0.000203 , 0.001 , 0.000475 , 0.000464 , 0.00106 , 0.000523 , -0.000281 , 0.000132 , -2.44e-05 , -0.000762 , -0.000489 , -6.63e-05 },
  { 0.000707 , 0.00201 , 0.00207 , 0.00147 , 0.00181 , 0.00265 , -9.27e-05 , -0.000863 , -0.000497 , -0.00028 , -0.000821 , 7.88e-05 },
  { 0.00133 , 0.00155 , 0.00136 , 0.00136 , 0.00155 , 0.00139 , -0.000212 , -0.000295 , -0.00067 , -0.000596 , 0.000263 , -0.000287 },
  { 0.000444 , 0.000197 , 0.000446 , 5.43e-05 , -0.000232 , -0.000543 , -0.00094 , -0.000797 , -0.000624 , -0.000737 , -0.000781 , -0.000657 },
  { 0.000157 , 0.000202 , -0.000169 , -0.000444 , -0.000434 , -0.000902 , -0.0015 , -0.00141 , -0.00154 , -0.00147 , -0.00141 , -0.00116 },
  { 0.000945 , 0.000666 , 0.000954 , 0.000526 , 0.000714 , 0.000461 , 0.00105 , 0 , 0 , 0 , 0 , 0 },
  { 0.000543 , 0.00043 , 0.000527 , 0.00033 , 0.000586 , 0.000737 , 0.00065 , 0.000351 , -0.000478 , -0.00115 , 2.38e-05 , -0.000237 },
  { -0.000644 , -0.000751 , -0.000288 , -0.000726 , -0.000879 , -0.000496 , -0.000765 , -0.000472 , -0.000507 , -0.000369 , -0.000369 , -0.000609 },
  { -0.00107 , -0.00082 , -0.000995 , -0.000987 , -0.000978 , -0.00113 , -0.000399 , -8.63e-05 , -8.52e-06 , 0.000562 , 0.000521 , 0.000368 },
  { -0.0009 , -0.000544 , -0.000525 , -0.000643 , -0.000497 , -0.000681 , 0.000283 , 0.000284 , 0.000252 , 0.000792 , 0.00084 , 0.000589 },
  { -0.00091 , -0.000822 , -0.000716 , -0.000655 , -0.000486 , -0.000996 , -0.000367 , -1.48e-05 , -0.0009 , -0.000892 , -0.000219 , -0.000516 }
  }; // Hardcoded, 15o correction matrix for the amplitude parameter of the correction function, bins are [phibin][zbin]
  Double_t correctedDCA = 1.;
  if(pt > 0.) correctedDCA = oldIP - correctionMatrix[PhiBin][zBin]*1.22/(1.+1./(0.57+9.8/(pt*pt)));
  dcaxy = correctedDCA;
}



Int_t AliAnalysisTaskHFEIPCorrection::ReturnRunBin(Int_t RunNr)
{
  Int_t runList[183]={
  246994, 246991, 246989, 246984, 246982, 246980, 246949, 246948, 246945, 246942, 246937, 246930, 246928, 246871, 246870, 246867, 246865, 246864, 246859, 246858, 246855, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246806, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246755, 246751, 246750, 246676, 246675, 246671, 246648, 246639, 246583, 246575, 246568, 246567, 246553, 246543, 246540, 246495, 246493, 246488, 246487, 246434, 246433, 246431, 246428, 246424, 246392, 246391, 246390, 246276, 246275, 246272, 246271, 246225, 246222, 246220, 246217, 246187, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246148, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245996, 245963, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245793, 245785, 245775, 245766, 245759, 245752, 245738, 245731, 245729, 245705, 245702, 245700, 245692, 245683, 245554, 245545, 245544, 245543, 245542, 245540, 245535, 245507, 245505, 245504, 245501, 245497, 245496, 245454, 245453, 245452, 245450, 245446, 245441, 245439, 245411, 245410, 245409, 245407, 245401, 245397, 245396, 245353, 245349, 245347, 245346, 245345, 245343, 245341, 245259, 245256, 245253, 245233, 245232, 245231, 245230, 245152, 245151, 245148, 245146, 245145, 245068, 245066, 245064, 245061, 244983, 244982, 244980, 244975, 244972, 244918, 244917, 244911, 244889, 244827, 244824};
  for(int i=0;i<183;i++)
    if(runList[i]==RunNr) return i+1;
  return 0;
}


//________________________________________________________________________
void AliAnalysisTaskHFEIPCorrection::Process(AliAODEvent *const aodEvent)
{
  // Main loop
  // Called for each event
  EventSelectionSteps->Fill(4);
  //if(!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB))) //return;
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

  //Float_t centrality = -1.;  // In Run2, centrality is estimated differently
  //AliCentrality *hicent = aodEvent->GetCentrality();
  //centrality = hicent->GetCentralityPercentile("V0M");

Float_t centrality = -1;
Float_t fV0Cent = -999.;
Float_t fV0CentCalib = -999.;
 
AliMultSelection *MultSelection = 0x0;
MultSelection = (AliMultSelection*)aodEvent->FindListObject("MultSelection");
if(!MultSelection){
  AliWarning("AliMultSelection object not found!");
}else{
  fV0Cent = MultSelection->GetMultiplicityPercentile("V0M", false);
  fV0CentCalib = MultSelection->GetMultiplicityPercentile("V0M", true);
  centrality = fV0CentCalib;
}

  const AliAODVertex *vertex = aodEvent->GetPrimaryVertex();
  const AliAODVertex *vertexSPD = aodEvent->GetPrimaryVertexSPD();
  Double_t vcov[6];
  vertex->GetCovMatrix(vcov);
  Double_t vtx[3];
  Double_t vtxSPD[3];
  vertex->GetXYZ(vtx);
  vertexSPD->GetXYZ(vtxSPD);

  AliAODVertex *correctedVertex;
  Bool_t madecorvtx = kFALSE;

  bool analyzeEvent=(TMath::Sqrt(vcov[5]) < 0.25 && TMath::Abs(vtx[2])<10. && TMath::Abs(vtx[2] - vtxSPD[2]) < 0.5 && centrality <= 100.);
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
  Double_t IP=0.;

  AliAODv0 * v0;
  AliAODTrack* V0Daughter[2];
  Int_t V0MotherPdg, V0Daughter1Pdg, V0Daughter2Pdg;
  double recoRadius=0.;
  Double_t IPCorrected = 0.;
  Int_t RunBin = ReturnRunBin(aodEvent->GetRunNumber());

  if(analyzeEvent)
  {
      correctedVertex = CorrectVertex(aodEvent, vtx[2]); // created without problematic ITS regions
      madecorvtx = true;
      Double_t vtxcorr[3];
      correctedVertex->GetXYZ(vtxcorr);
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
        
    if(centrality>=20.0 && centrality<=40.0)
    {
      EP2040->Fill(V0PlanePhi);
      //EP2040Corrected->Fill(epCorr);
      EP2040V0A->Fill(V0APlanePhi);
      EP2040V0C->Fill(V0CPlanePhi);
    }
      double DPhi =0.;
      Float_t dcaxy = -999., dcaz = -999.;
      Double_t dcaErr, dcaxyD;
      Double_t dcaxyC;
      AliAODTrack *track = 0x0;
      for(Int_t itrack = 0; itrack < fInputEvent->GetNumberOfTracks(); itrack++)
      {
        // Run track loop
        track = dynamic_cast<AliAODTrack *>(fInputEvent->GetTrack(itrack));
        if(!track) continue;
        if(!PassesTrackCuts(track)) continue;
        if(PassesElectronPID(track, pid))
        {
          if(centrality>=20.0 && centrality<=50.0) TPCnSigma->Fill(track->Pt(),pid->NumberOfSigmasTPC(track, AliPID::kElectron));
          DPhi = track->Phi() - V0PlanePhi - TMath::Pi();
          while(DPhi<0.) DPhi += TMath::Pi();
          
          fExtraCuts->GetImpactParameters((AliVTrack *)track,dcaxy,dcaz);
          fExtraCuts->GetHFEImpactParameters((AliVTrack *)track,dcaxyD,dcaErr);
          IP = dcaxyD*track->Charge()*TMath::Sign(1.,aodEvent->GetMagneticField());
          if(centrality>=20.0 && centrality<=50.0) fIPData->Fill(track->Pt(), dcaxyD);
          if(track->Pt() > 0.5)
          {
            DeltaPhi->Fill(DPhi);
            if(DPhi/2. < TMath::Pi()/4.) // IP
            {
              if(centrality>=20.0 && centrality<=40.0) fpTIP2040IP->Fill(track->Pt(), IP);
              if(centrality>=30.0 && centrality<=50.0) fpTIP3050IP->Fill(track->Pt(), IP);
            }
            else
            {
              if(centrality>=20.0 && centrality<=40.0) fpTIP2040OOP->Fill(track->Pt(), IP);
              if(centrality>=30.0 && centrality<=50.0) fpTIP3050OOP->Fill(track->Pt(), IP);
            }
          }
        }
        if((PassesPionPID(track, pid) || PassesKaonPID(track, pid)) && track->Pt() > 0.5) // To avoid calling the IP calculation for all tracks
        {
          fExtraCuts->GetHFEImpactParameters((AliVTrack *)track,dcaxyD,dcaErr);
          GetTrackImpactParameter(aodEvent, track, correctedVertex, dcaxyC);
          IP = dcaxyD*track->Charge()*TMath::Sign(1.,aodEvent->GetMagneticField());
          GetCorrectedImpactParameter(track, vtx[2], IPCorrected);
          if(PassesPionPID(track, pid))
          {
            if(centrality>=0.0 && centrality<=80.0 && track->Pt()>1.)
            {
              if(track->Pt()>3. && IsInMisalignedRegion(track, vtx[2])>0) fDCAvsCorrected->Fill(dcaxyD, dcaxyC);
              fDCAPhiZHadrons->Fill(track->Phi(), dcaxyD, vtx[2]+4.5*TMath::Tan(TMath::Pi()/2.-2.*TMath::ATan(TMath::Exp(-track->Eta()))));
              fDCAPhiZHadronsC->Fill(track->Phi(), dcaxyC, vtx[2]+4.5*TMath::Tan(TMath::Pi()/2.-2.*TMath::ATan(TMath::Exp(-track->Eta()))));
              fDCAPhiZHadronsC2->Fill(track->Phi(), IPCorrected, vtx[2]+4.5*TMath::Tan(TMath::Pi()/2.-2.*TMath::ATan(TMath::Exp(-track->Eta()))));
              fDCAPhiZHadrons2nd->Fill(track->Phi(), dcaxyD, vtx[2]+7.*TMath::Tan(TMath::Pi()/2.-2.*TMath::ATan(TMath::Exp(-track->Eta()))));
              if(IsInMisalignedRegion(track, vtx[2])>0)fDCARegionRun->Fill(dcaxyD, IsInMisalignedRegion(track, vtx[2]) ,RunBin);
            }
            if(centrality>=0.0 && centrality<=80.0 && track->Pt()>0.5) fpTPhiZHadrons->Fill(track->Phi(), track->Pt(), vtx[2]+7.*TMath::Tan(TMath::Pi()/2.-2.*TMath::ATan(TMath::Exp(-track->Eta()))));
            if(centrality>=0.0 && centrality<=80.0 && track->Pt()>0.5) fDCAPhipTHadrons->Fill(dcaxyD, track->Phi(), track->Pt());
            if(centrality>=0.0 && centrality<=80.0 && track->Pt()>0.5) fDCAPhipTHadronsC->Fill(dcaxyC, track->Phi(), track->Pt());
            if(centrality>=0.0 && centrality<=80.0 && track->Pt()>0.5) fDCAPhipTHadronsC2->Fill(IPCorrected, track->Phi(), track->Pt());
            if(centrality>=20.0 && centrality<=50.0) fDCAWErrHadrons->Fill(track->Pt(), IP, dcaxyD/dcaErr);
            fDCAHadronsFineBins->Fill(track->Pt(), IP, centrality);
          }
          if(PassesKaonPID(track, pid))
          {
            if(centrality>=20.0 && centrality<=50.0) fDCAKaons->Fill(track->Pt(), IP);
            if(centrality>=20.0 && centrality<=50.0) fDCAWErrKaons->Fill(track->Pt(), IP, dcaxyD/dcaErr);
            fDCAKaonsFineBins->Fill(track->Pt(), IP, centrality);
          }
        }
      }
    fAODV0Cuts->SetEvent(aodEvent);
    for(int i= 0;i<aodEvent->GetNumberOfV0s();i++)
    {
      //
      v0 = aodEvent->GetV0(i);
      //if(!PassesTrackCuts(track));
      if(fAODV0Cuts->ProcessV0(v0, V0MotherPdg, V0Daughter1Pdg, V0Daughter2Pdg))
      {
        recoRadius = v0->RadiusSecVtx();
        if(TMath::Abs(V0MotherPdg)==310 && recoRadius>0.5)
        {
          V0Daughter[0] = dynamic_cast<AliAODTrack *> (v0->GetSecondaryVtx()->GetDaughter(0)); // This is how to get the daughter particles in AODs, apparently
          V0Daughter[1] = dynamic_cast<AliAODTrack *> (v0->GetSecondaryVtx()->GetDaughter(1));
          if(PassesMinimalTrackCuts(V0Daughter[0])) fPionV0pTRNoCuts->Fill(V0Daughter[0]->Pt(), recoRadius, centrality);
          if(PassesMinimalTrackCuts(V0Daughter[1])) fPionV0pTRNoCuts->Fill(V0Daughter[1]->Pt(), recoRadius, centrality);
          if(PassesITSTrackCuts(V0Daughter[0])) fPionV0pTRWithCuts->Fill(V0Daughter[0]->Pt(), recoRadius, centrality);
          if(PassesITSTrackCuts(V0Daughter[1])) fPionV0pTRWithCuts->Fill(V0Daughter[1]->Pt(), recoRadius, centrality);
          if(centrality>=20.0 && centrality<=50.0){
          if(PassesMinimalTrackCuts(V0Daughter[0])) fPionV0pTTPC->Fill(V0Daughter[0]->Pt(), pid->NumberOfSigmasTPC(V0Daughter[0], AliPID::kPion));
          if(PassesMinimalTrackCuts(V0Daughter[1])) fPionV0pTTPC->Fill(V0Daughter[1]->Pt(), pid->NumberOfSigmasTPC(V0Daughter[1], AliPID::kPion));
          if(PassesMinimalTrackCuts(V0Daughter[0]) && V0Daughter[0]->GetTPCNcls()>=110 && V0Daughter[0]->GetTPCsignalN()>=80) fPionV0pTTPCWithCuts->Fill(V0Daughter[0]->Pt(), pid->NumberOfSigmasTPC(V0Daughter[0], AliPID::kPion));
          if(PassesMinimalTrackCuts(V0Daughter[1]) && V0Daughter[1]->GetTPCNcls()>=110 && V0Daughter[1]->GetTPCsignalN()>=80) fPionV0pTTPCWithCuts->Fill(V0Daughter[1]->Pt(), pid->NumberOfSigmasTPC(V0Daughter[1], AliPID::kPion));}
        }
      }
    }
  }
  
  if(madecorvtx) delete correctedVertex;
  PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskHFEIPCorrection::Terminate(const Option_t *)
{
//   treeoutput->Write();

}
