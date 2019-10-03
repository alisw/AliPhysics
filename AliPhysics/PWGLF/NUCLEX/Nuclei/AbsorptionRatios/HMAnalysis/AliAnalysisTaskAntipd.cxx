#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TProfile.h"

#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>
#include <AliInputEventHandler.h>
#include <AliAODTrack.h>
#include "AliPIDResponse.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"

#include "AliAnalysisTaskAntipd.h"

class AliAnalysisTaskAntipd;
using namespace std;

ClassImp(AliAnalysisTaskAntipd)

//________________________________________________________________________
AliAnalysisTaskAntipd::AliAnalysisTaskAntipd()
: AliAnalysisTaskSE()
,fVEvent(0)
,fAODEvent(0)
,fOutputList(0)
,fOutputEvent(0)
,fOutputProtons(0)
,fOutputAProtons(0)
,fOutputDeuterons(0)
,fOutputADeuterons(0)
,fOutputHe3(0)
,fOutputAHe3(0)
,fHistTrackCuts(0)
,fEventStat(0)
,fVertexZ(0)
,fMultPercentileV0M(0)
,fMultPercentileV0MZoomed(0)
,fHistsProton()
,fHistsAProton()
,fHistsDeuteron()
,fHistsADeuteron()
,fHistsHe3()
,fHistsAHe3()
,fPIDResponse(0)
,fFilterBit(256)   // TPC and ITS tracking requirements correspond to AOD Filter Bit 256
,fLowPCut(0.0)
,fHighPCut(1e30)
,fEtaCut(1.0)
,fMinClIts(0)
,fMaxDCAxy(10.)
,fMaxDCAz(10.)
,fMaxTPCnSigma(10.)
,fMaxTOFnSigma(10.)
,fUseTOFPidCut(kFALSE)
,fMinITSnSigma(-10.)
,fMaxITSnSigma(10.)
,fMomTOFProt(10.)
,fMomTOFDeut(10.)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}

//________________________________________________________________________
AliAnalysisTaskAntipd::AliAnalysisTaskAntipd(const char* name)
: AliAnalysisTaskSE(name)
,fVEvent(0)
,fAODEvent(0)
,fOutputList(0)
,fOutputEvent(0)
,fOutputProtons(0)
,fOutputAProtons(0)
,fOutputDeuterons(0)
,fOutputADeuterons(0)
,fOutputHe3(0)
,fOutputAHe3(0)
,fHistTrackCuts(0)
,fEventStat(0)
,fVertexZ(0)
,fMultPercentileV0M(0)
,fMultPercentileV0MZoomed(0)
,fHistsProton()
,fHistsAProton()
,fHistsDeuteron()
,fHistsADeuteron()
,fHistsHe3()
,fHistsAHe3()
,fPIDResponse(0)
,fFilterBit(256)
,fLowPCut(0.0)
,fHighPCut(1e30)
,fEtaCut(1.0)
,fMinClIts(0)
,fMaxDCAxy(10.)
,fMaxDCAz(10.)
,fMaxTPCnSigma(10.)
,fMaxTOFnSigma(10.)
,fUseTOFPidCut(kFALSE)
,fMinITSnSigma(-10.)
,fMaxITSnSigma(10.)
,fMomTOFProt(10.)
,fMomTOFDeut(10.)
{
    // constructor

    fHistsProton.clear();
    fHistsAProton.clear();
    fHistsDeuteron.clear();
    fHistsADeuteron.clear();
    fHistsHe3.clear();
    fHistsAHe3.clear();


    DefineInput(0, TChain::Class());    // a 'chain' of events created by the analysis manager automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
    // one can add more output objects by calling DefineOutput(2, classname::Class())
    // make sure to connect them properly in AddTask and to call PostData() for all of them!

}

//________________________________________________________________________
AliAnalysisTaskAntipd::~AliAnalysisTaskAntipd()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;
    }
}

//________________________________________________________________________
void AliAnalysisTaskAntipd::UserCreateOutputObjects()
{
    // Called once at the start of the analysis

    // retrieve PID object from the input handler
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    // create output lists
    fOutputList = new TList();
    fOutputList->SetName("Output");
    fOutputList->SetOwner(kTRUE); // the list is owner of all objects it contains and will delete them

    fOutputEvent = new TList();
    fOutputEvent->SetName("Event");
    fOutputEvent->SetOwner(kTRUE);

    fOutputProtons = new TList();
    fOutputProtons->SetName("Protons");
    fOutputProtons->SetOwner(kTRUE);

    fOutputAProtons = new TList();
    fOutputAProtons->SetName("AProtons");
    fOutputAProtons->SetOwner(kTRUE);

    fOutputDeuterons = new TList();
    fOutputDeuterons->SetName("Deuterons");
    fOutputDeuterons->SetOwner(kTRUE);

    fOutputADeuterons = new TList();
    fOutputADeuterons->SetName("ADeuterons");
    fOutputADeuterons->SetOwner(kTRUE);

    fOutputHe3 = new TList(); // helium3  //
    fOutputHe3->SetName("Helium3");
    fOutputHe3->SetOwner(kTRUE);

    fOutputAHe3 = new TList(); // Anti-helium3  //
    fOutputAHe3->SetName("AHelium3");
    fOutputAHe3->SetOwner(kTRUE);

    fOutputList->Add(fOutputEvent);
    fOutputList->Add(fOutputProtons);
    fOutputList->Add(fOutputAProtons);
    fOutputList->Add(fOutputDeuterons);
    fOutputList->Add(fOutputADeuterons);
    fOutputList->Add(fOutputHe3);
    fOutputList->Add(fOutputAHe3);


    // create event histograms
    fEventStat = new TH1I("fEventStat", "Event statistics", kNbinsEvent, 0, kNbinsEvent);
    fEventStat->GetXaxis()->SetBinLabel(1,"After Phys. Sel. and trigger");
    fEventStat->GetXaxis()->SetBinLabel(2,"kINT7 selected (cross-check)");
    fEventStat->GetXaxis()->SetBinLabel(3,"DAQ imcomplete (cross-check)");
    fEventStat->GetXaxis()->SetBinLabel(4,"V0 timing (cross-check)");
    fEventStat->GetXaxis()->SetBinLabel(5,"Cluster-tracklet cut (cross-check)");
    fEventStat->GetXaxis()->SetBinLabel(6,"Vertex Z position");
    fEventStat->GetXaxis()->SetBinLabel(7,"N of vtx contrib");
    fEventStat->GetXaxis()->SetBinLabel(8,"SPD pile-up in mult. bins");
    fEventStat->GetXaxis()->SetBinLabel(9,"Vertex displacement");
    fEventStat->GetXaxis()->SetBinLabel(10,"Vertex Z resolution");
    fOutputEvent->Add(fEventStat);

    fVertexZ = new TH1F("fVertexZ", "Vertex Z distribution", 300, -15., 15.);
    fVertexZ->GetXaxis()->SetTitle("vertex Z, cm");
    fVertexZ->Sumw2();
    fOutputEvent->Add(fVertexZ);

  // percentile --- Measuring the Multiplicity
    fMultPercentileV0M = new TH1F("fMultPercentileV0M","", 250., 0., 250.);
    fVertexZ->GetXaxis()->SetTitle("Multiplicity Percentile");
    fVertexZ->Sumw2();
    fOutputEvent->Add(fMultPercentileV0M);

    fMultPercentileV0MZoomed = new TH1F("fMultPercentileV0MZoomed","", 100., 0., 1.);
    fVertexZ->GetXaxis()->SetTitle("Multiplicity Percentile");
    fVertexZ->Sumw2();
    fOutputEvent->Add(fMultPercentileV0MZoomed);

    // track cuts config
    fHistTrackCuts = new TProfile("fHistTrackCuts", "TrackCuts config", 9, 0, 9);
    fHistTrackCuts->GetXaxis()->SetBinLabel(1,"Filter Bit");
    fHistTrackCuts->GetXaxis()->SetBinLabel(2,"low #it{p} cut");
    fHistTrackCuts->GetXaxis()->SetBinLabel(3,"#eta cut");
    fHistTrackCuts->GetXaxis()->SetBinLabel(4,"Min. N ITS cl");
    fHistTrackCuts->GetXaxis()->SetBinLabel(5,"Max. DCAxy");
    fHistTrackCuts->GetXaxis()->SetBinLabel(6,"Max. DCAz");
    fHistTrackCuts->GetXaxis()->SetBinLabel(7,"TPCnSigma cut");
    fHistTrackCuts->GetXaxis()->SetBinLabel(8,"TOFnSigma cut");
    fHistTrackCuts->GetXaxis()->SetBinLabel(9,"TOFmomentum prot");
    fHistTrackCuts->GetXaxis()->SetBinLabel(10,"TOFmomentum deut");
    fOutputList->Add(fHistTrackCuts);

    fHistTrackCuts->Fill(0.0,fFilterBit);
    fHistTrackCuts->Fill(1,fLowPCut);
    fHistTrackCuts->Fill(2,fEtaCut);
    fHistTrackCuts->Fill(3,fMinClIts);
    fHistTrackCuts->Fill(4,fMaxDCAxy);
    fHistTrackCuts->Fill(5,fMaxDCAz);
    fHistTrackCuts->Fill(6,fMaxTPCnSigma);
    fHistTrackCuts->Fill(7,fMaxTOFnSigma);
    fHistTrackCuts->Fill(8,fMomTOFProt);
    fHistTrackCuts->Fill(9,fMomTOFDeut);

    // (anti)proton histograms
    CreateHistosTrack(fHistsProton);
    for (Int_t i=0;i<(int)(fHistsProton.size()); i++){
        fOutputProtons->Add(fHistsProton.at(i));
    }

    CreateHistosTrack(fHistsAProton);
    for (Int_t i=0;i<(int)(fHistsProton.size()); i++){
        fOutputAProtons->Add(fHistsAProton.at(i));
    }

    // (anti)deuteron histograms
    CreateHistosTrack(fHistsDeuteron);
    for (Int_t i=0;i<(int)(fHistsDeuteron.size()); i++){
        fOutputDeuterons->Add(fHistsDeuteron.at(i));
    }

    CreateHistosTrack(fHistsADeuteron);
    for (Int_t i=0;i<(int)(fHistsADeuteron.size()); i++){
        fOutputADeuterons->Add(fHistsADeuteron.at(i));
    }

  // (anti)helium histograms
     CreateHistosTrack(fHistsHe3);
     for (Int_t i=0;i< (int)(fHistsHe3.size()); i++){
     fOutputHe3->Add(fHistsHe3.at(i));
    }

    CreateHistosTrack(fHistsAHe3);
    for (Int_t i=0;i< (int)(fHistsAHe3.size()); i++){
    fOutputAHe3->Add(fHistsAHe3.at(i));
   }

    PostData(1, fOutputList);
    // postdata will notify the analysis manager of changes / updates to the
    // fOutputList object. the manager will in the end take care of writing output to file
    // so it needs to know what's in the output
}

//________________________________________________________________________
void AliAnalysisTaskAntipd::UserExec(Option_t *)
{
    // Main loop
    // Called for each event

    fVEvent = dynamic_cast<AliVEvent*>(InputEvent());
    if (!fVEvent) {
        //printf("fVEvent not available\n");
        return;
    }

    fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAODEvent) {
        printf("fAODEvent not available\n");
        return;
    }

    fEventStat->Fill(kSelectedEvents); // all events after trigger and physics selection
//    // only for cross-check: is kINT7 selected
//    UInt_t fSelectMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
//    Bool_t isINT7selected = fSelectMask & AliVEvent::kINT7;
//    if (!isINT7selected){
//        PostData(1, fOutputList);
//        return;
//    }
    fEventStat->Fill(kINT7selected);

//    // only for cross-check: Incomplete events from DAQ
//    if (fVEvent->IsIncompleteDAQ()) { // probably not implemented for AOD?
//        PostData(1, fOutputList);
//        return;
//    }
    fEventStat->Fill(kDAQincomplete);

//    //only for cross-check: VZERO timing desicion
//    Int_t fV0ADecision;
//    Int_t fV0CDecision;
//    AliVVZERO* vzeroData = fVEvent->GetVZEROData();
//    fV0ADecision = vzeroData->GetV0ADecision();
//    fV0CDecision = vzeroData->GetV0CDecision();
//    if (!fV0ADecision || !fV0CDecision) {
//        PostData(1, fOutputList);
//        return;
//    }
    fEventStat->Fill(kV0timing);

//    //only for cross-check: SPD clusters vs tracklets cut
//    AliAnalysisUtils *utils = new AliAnalysisUtils();
//    //utils->SetASPDCvsTCut(100.); // 65 by default
//    //utils->SetBSPDCvsTCut(5.); // 4 by default
//    Bool_t ClustersVsTrackletBG = utils->IsSPDClusterVsTrackletBG(fVEvent);
//    if (ClustersVsTrackletBG) {
//        PostData(1, fOutputList);
//        return;
//    }
    fEventStat->Fill(kClusterTrackletCut);

    // vertex filter
    const AliVVertex *vertex = fVEvent->GetPrimaryVertex();
    Bool_t lVertexAcceptable = (TMath::Abs(vertex->GetZ()) <= 10.);
    Bool_t lVertexNcontrib   = (vertex->GetNContributors() > 0);

    if (!lVertexAcceptable){
        PostData(1, fOutputList);
        return;
    }
    fEventStat->Fill(kVertexZ);

    if (!lVertexNcontrib){
        PostData(1, fOutputList);
        return;
    }
    fEventStat->Fill(kVtxNcontrib);

    // SPD pile-up in mult. bins
    if (fVEvent->IsPileupFromSPDInMultBins()){
        PostData(1, fOutputList);
        return;
    }
    fEventStat->Fill(kSPDPileUp);

//    // spd vertex and its resolution
//    const AliVVertex* vtxSPD = fVEvent->GetPrimaryVertexSPD();
//    Double_t cov[6]={0};
//    vtxSPD->GetCovarianceMatrix(cov);
//    Double_t zRes = TMath::Sqrt(cov[5]);
//    if (TMath::Abs(vtxSPD->GetZ() - vertex->GetZ()) > 0.3){
//        PostData(1, fOutputList);
//        return;
//    }
    fEventStat->Fill(kVtxDisplace);

//    if (vtxSPD->IsFromVertexerZ() && zRes > 0.3 ) {
//        PostData(1, fOutputList);
//        return;
//    }
    fEventStat->Fill(kVtxRes);

    // fill event histograms
    fVertexZ->Fill(vertex->GetZ());

    // Centrality percentile
    Float_t lPercentile = 300;
    AliMultSelection *MultSelection = 0x0;
    MultSelection = (AliMultSelection * ) fVEvent->FindListObject("MultSelection");
    if( !MultSelection) {
   //If you get this warning (and lPercentiles 300) please check that the   AliMultSelectionTask actually ran (before your task)
   AliWarning("AliMultSelection object not found!");
}else{
   lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
}

   fMultPercentileV0M->Fill(lPercentile);
   fMultPercentileV0MZoomed->Fill(lPercentile);

    // track loop
    for (Int_t iTrack = 0; iTrack < fAODEvent->GetNumberOfTracks(); iTrack++){

        AliAODTrack* track = (AliAODTrack*)fAODEvent->GetTrack(iTrack);
        if (!track) { Printf("ERROR: Could not receive track %d", iTrack); continue; }
        // Int_t label = track->GetLabel();
        Double_t trackP   = track->P();
        Double_t trackEta = track->Eta();

        // ===================== track cuts =====================
        // filter bit, pt and eta cuts, N of ITS clusters
        if (!(track->TestFilterBit(fFilterBit)))          continue;
        if (trackP < fLowPCut || trackP > fHighPCut)      continue;
        if (TMath::Abs(trackEta) > fEtaCut)               continue;
        if (track->GetITSNcls() < fMinClIts)              continue;

        // DCA information
        double DCAxy = -99.;
        double DCAz  = -99.;
        double dcaVals[2] = {-99., -99.};
        double covar[3]={0.,0.,0.};
        AliAODTrack copy(*track);
        const AliVVertex *AODeventVtx = copy.GetAODEvent()->GetPrimaryVertex();
        Double_t AODeventMagneticF = copy.GetAODEvent()->GetMagneticField();
        if (copy.PropagateToDCA(AODeventVtx, AODeventMagneticF,10, dcaVals, covar)){
            DCAxy = dcaVals[0];
            DCAz  = dcaVals[1];
        } else {
            DCAxy = -99.; // track->DCA();
            DCAz  = -99.; // track->ZAtDCA();
        }

        // cut on DCA
        if (TMath::Abs(DCAz ) > fMaxDCAz)     continue;
        if (TMath::Abs(DCAxy) > fMaxDCAxy)    continue;

        // ===================== PID selection =====================
        // TPC
        Float_t nSigmaTPCprot = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        Float_t nSigmaTPCdeut = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
        Float_t nSigmaTPCHe3 = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kHe3); //

        // ITS
        Float_t nSigmaITSprot = fPIDResponse->NumberOfSigmasITS(track, AliPID::kProton);
        Float_t nSigmaITSdeut = fPIDResponse->NumberOfSigmasITS(track, AliPID::kDeuteron);

        // TOF
        Bool_t isTOF = (track->GetStatus() & AliVTrack::kTOFout)
                    && (track->GetStatus() & AliVTrack::kTIME);


        // optinally cut on TOFnSigma at high momentum
        Float_t nSigmaTOFprot = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
        Float_t nSigmaTOFdeut = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron);

        Bool_t TOFpid_prot = isTOF;
        Bool_t TOFpid_deut = isTOF;

        if (fUseTOFPidCut){
            TOFpid_prot = (isTOF && TMath::Abs(nSigmaTOFprot) < fMaxTOFnSigma);
            TOFpid_deut = (isTOF && TMath::Abs(nSigmaTOFdeut) < fMaxTOFnSigma);
        }

        // is (anti-)proton selected?
        Bool_t isProtonSelected   = (TMath::Abs(nSigmaTPCprot) < fMaxTPCnSigma) &&
                                    (!(trackP > fMomTOFProt && !TOFpid_prot));

        // is (anti-)deuteron selected?
        Bool_t isDeuteronSelected = (TMath::Abs(nSigmaTPCdeut) < fMaxTPCnSigma) &&
                                    (!( trackP < fMomTOFDeut &&
                                    (nSigmaITSdeut < fMinITSnSigma ||
                                     nSigmaITSdeut > fMaxITSnSigma))) &&
                                    (!(trackP > fMomTOFDeut && !TOFpid_deut));

        // is (anti-)Helium-3 selected?
      Bool_t isHe3Selected   = (TMath::Abs(nSigmaTPCHe3) < fMaxTPCnSigma);

        // ===================== fill track histos =====================
        if (isProtonSelected){
            if (track->Charge() > 0){
                FillHistosTrack(fHistsProton,track);

            } else if (track->Charge() < 0){
                FillHistosTrack(fHistsAProton,track);

            }
        }

        if (isDeuteronSelected){
            if (track->Charge() > 0){
                FillHistosTrack(fHistsDeuteron,track);

            } else if (track->Charge() < 0){
                FillHistosTrack(fHistsADeuteron,track);

            }
        }

        if (isHe3Selected){
           if (track->Charge() > 0){
              FillHistosTrack(fHistsHe3,track);

             } else if (track->Charge() < 0){
                 FillHistosTrack(fHistsAHe3,track);

             }
         }

    } // track loop

    PostData(1, fOutputList);


} // end of UserExec

//________________________________________________________________________
void AliAnalysisTaskAntipd::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}

//________________________________________________________________________
void AliAnalysisTaskAntipd::CreateHistosTrack(vector<TH1*> &histos)
{

    TH1F* fHistP = new TH1F("fHistP","p distribution;#it{p}, GeV/#it{c};counts",100, 0., 10.);
    fHistP->Sumw2();
    histos.push_back(fHistP);  //0

    TH1F* fHistPt = new TH1F("fHistPt","pt distribution;#it{p}_{T}, GeV/#it{c};counts",100, 0., 10.);
    fHistPt->Sumw2();
    histos.push_back(fHistPt);  //1

    TH2F* fHistEtaPhi = new TH2F("fHistEtaPhi","eta-phi distribution;#eta;#varphi",100, -1.0, 1.0, 360, 0, 2*TMath::Pi());
    histos.push_back(fHistEtaPhi);  //2

    // tracking
    TH1F* fHistTPCCrossedRows = new TH1F("fHistTPCCrossedRows","TPC crossed rows;N TPC crossed rows;counts",200,0,200);
    histos.push_back(fHistTPCCrossedRows);  //3

    TH1F* fHistTPCClusters = new TH1F("fHistTPCClusters","TPC clusters;N TPC clusters;counts",200,0,200);
    histos.push_back(fHistTPCClusters); //4

    TH1F* fHistTPCCRoverFind = new TH1F("fHistTPCCRoverFind","TPC crossed rows / findable;TPC crossed rows / findable;counts",150,0,1.5);
    histos.push_back(fHistTPCCRoverFind);  //5

    TH1F* fHistTPCFracShared = new TH1F("fHistTPCFracShared","TPC shared clusters;TPC shared clusters;counts",250,-1.0,1.5);
    histos.push_back(fHistTPCFracShared); //6

    TH1F* fHistTPCSignalN = new TH1F("fHistTPCSignalN","Number of PID clusters TPC;N of TPC PID clusters;counts",200,0.,200.);
    histos.push_back(fHistTPCSignalN); //7

    TH1F* fHistTPCSignalNfrac = new TH1F("fHistTPCSignalNfrac","Fraction of PID clusters TPC;Fraction of TPC PID clusters;counts",120,0.,1.2);
    histos.push_back(fHistTPCSignalNfrac);  //8

    TH1F* fHistITSnCls = new TH1F("fHistITSnCls","Number of ITS clusters;N of ITS clusters;counts",10,0.,10.);
    histos.push_back(fHistITSnCls);  //9

    TH1F* fHistChi2 = new TH1F("fHistChi2","track chi2;#chi^{2};counts",100,0,10);
    histos.push_back(fHistChi2);  //10

    TH1F* fHistDCAxy = new TH1F("fHistDCAxy","DCA xy;DCA_{xy};counts",400,-2.0,2.0);
    histos.push_back(fHistDCAxy); //11

    TH1F* fHistDCAz = new TH1F("fHistDCAz","DCA z;DCA_{z};counts",400,-2.0,2.0);
    histos.push_back(fHistDCAz);  //12

    TH2F* fHistDCAxyDCAz = new TH2F("fHistDCAxyDCAz","DCAxy vs DCAz;DCA_{xy};DCA_{z}",400,-2.0,2.0, 400, -2.0, 2.0);
    histos.push_back(fHistDCAxyDCAz);  //13

    // PID histos
    TH2F* fHistTPCSignal = new TH2F("fHistTPCSignal","TPC dE/dx;#it{p}, GeV/#it{c};TPC dE/dx",500,0.,10.0,750,0.,750.);
    histos.push_back(fHistTPCSignal);  //14

    TH2F* fHistTPCnSigmaProt = new TH2F("fHistTPCnSigmaProt","TPC nSigma prot;#it{p}, GeV/#it{c};TPCn#sigma_{prot}",100,0.,10.0,100,-5.,5.);
    histos.push_back(fHistTPCnSigmaProt);  //15

    TH2F* fHistTPCnSigmaDeut = new TH2F("fHistTPCnSigmaDeut","TPC nSigma deut;#it{p}, GeV/#it{c};TPCn#sigma_{deut}",100,0.,10.0,100,-5.,5.);
    histos.push_back(fHistTPCnSigmaDeut); //16

    TH2F* fHistTPCnSigmaHe3 = new TH2F("fHistTPCnSigmaHe3","TPC nSigma He3;#it{p}, GeV/#it{c};TPCn#sigma_{he3}",100,0.,10.0,100,-5.,5.);
    histos.push_back(fHistTPCnSigmaHe3);  //17

   //TH2F* fHistTPCnSigmaAHe3 = new TH2F("fHistTPCnSigmaAHe3","TPC nSigma He3;#it{p}, GeV/#it{c};TPCn#sigma_{he3}",100,0.,10.0,100,-5.,5.);
    //histos.push_back(fHistTPCnSigmaAHe3);  // fHistsAHe3

    TH2F* fHistTOFSignal = new TH2F("fHistTOFSignal","TOF signal;#it{p}, GeV/#it{c};TOF signal",500,0.,10.0,200,0.,2.);
    histos.push_back(fHistTOFSignal); //18

    TH2F* fHistTOFnSigmaProt = new TH2F("fHistTOFnSigmaProt","TOF nSigma prot;#it{p}, GeV/#it{c};TOFn#sigma_{prot}",100,0.,10.0,100,-5.,5.);
    histos.push_back(fHistTOFnSigmaProt); //19

    TH2F* fHistTOFnSigmaDeut = new TH2F("fHistTOFnSigmaDeut","TOF nSigma deut;#it{p}, GeV/#it{c};TOFn#sigma_{deut}",100,0.,10.0,100,-5.,5.);
    histos.push_back(fHistTOFnSigmaDeut); //20

    TH2F* fHistTOFnSigmaPion = new TH2F("fHistTOFnSigmaPion","TOF nSigma pion;#it{p}, GeV/#it{c};TOFn#sigma_{pion}",100,0.,10.0,100,-5.,5.);
    histos.push_back(fHistTOFnSigmaPion); //21

    TH3F* fHistTOFmass2_DCAxy_p = new TH3F("fHistTOFmass2_DCAxy_p","TOF mass2 vs DCA_{xy} vs #it{p};#it{p}, GeV/#it{c};TOF mass^{2}, (GeV/#it{c}^{2})^{2};DCA_{xy}, cm;",50,0.,5.,500,0.0,5.0,200,-1.0,1.0);
    histos.push_back(fHistTOFmass2_DCAxy_p); //22

    TH2F* fHistITSSignal = new TH2F("fHistITSSignal","ITS signal;#it{p}, GeV/#it{c};ITS signal",500,0.,10.0,200,0.,200.);
    histos.push_back(fHistITSSignal); //23

    TH2F* fHistITSnSigmaProt = new TH2F("fHistITSnSigmaProt","ITS nSigma prot;#it{p}, GeV/#it{c};ITSn#sigma_{prot}",100,0.,10.0,200,-10.,10.);
    histos.push_back(fHistITSnSigmaProt); //24

    TH2F* fHistITSnSigmaDeut = new TH2F("fHistITSnSigmaDeut","ITS nSigma deut;#it{p}, GeV/#it{c};ITSn#sigma_{deut}",100,0.,10.0,200,-10.,10.);
    histos.push_back(fHistITSnSigmaDeut); //25

}

//________________________________________________________________________
void AliAnalysisTaskAntipd::FillHistosTrack(vector<TH1*> &histos, AliAODTrack *track)
{
    Double_t trackP   = track->P();
    Double_t trackPt  = track->Pt();
    Double_t trackEta = track->Eta();
    Double_t trackPhi = track->Phi();

    // TOF information
    Float_t fbetaTOF  = GetTOFBeta(track);
    Float_t fmass2TOF = GetMass2TOF(fbetaTOF,track);

    // PID information
    Float_t nSigmaTPCprot = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    Float_t nSigmaITSprot = fPIDResponse->NumberOfSigmasITS(track, AliPID::kProton);
    Float_t nSigmaTOFprot = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);

    Float_t nSigmaTPCdeut = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
    Float_t nSigmaITSdeut = fPIDResponse->NumberOfSigmasITS(track, AliPID::kDeuteron);
    Float_t nSigmaTOFdeut = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron);

    Float_t nSigmaTPCHe3 = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kHe3); //

    Float_t nSigmaTOFpion = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);

    // DCA information
    double DCAxy = -99.;
    double DCAz  = -99.;
    double dcaVals[2] = {-99., -99.};
    double covar[3]={0.,0.,0.};
    AliAODTrack copy(*track);
    const AliVVertex *AODeventVtx = copy.GetAODEvent()->GetPrimaryVertex();
    Double_t AODeventMagneticF = copy.GetAODEvent()->GetMagneticField();
    if (copy.PropagateToDCA(AODeventVtx, AODeventMagneticF,10, dcaVals, covar)){
        DCAxy = dcaVals[0];
        DCAz  = dcaVals[1];
    } else {
        DCAxy = -99.; // track->DCA();
        DCAz  = -99.; // track->ZAtDCA();
    }

    // if DCAxy > cut for final analysis -> fill only 3d histogram
    if (TMath::Abs(DCAxy) > 0.1){

        ((TH3F*)(histos.at(22)))->Fill(trackP,fmass2TOF,DCAxy);  // changed from 21 histos to 22 when added He3
        return;

    } else {

        histos.at(0)->Fill(trackP);
        histos.at(1)->Fill(trackPt);
        histos.at(2)->Fill(trackEta,trackPhi);
        histos.at(3)->Fill(track->GetTPCClusterInfo(2, 1));
        histos.at(4)->Fill(track->GetTPCNcls());
        if (!track->GetTPCNclsF()) {
            histos.at(5)->Fill(0);
        } else {
            histos.at(5)->Fill(track->GetTPCClusterInfo(2, 1)/
                               double(track->GetTPCNclsF()));
        }
        if (!track->GetTPCNcls()){
            histos.at(6)->Fill(-1);
        } else {
            histos.at(6)->Fill(track->GetTPCnclsS()/double(track->GetTPCNcls()));
        }
        histos.at(7)->Fill(track->GetTPCsignalN());
        if (!track->GetTPCNcls()){
            histos.at(8)->Fill(-1);
        } else {
            histos.at(8)->Fill(track->GetTPCsignalN()/double(track->GetTPCNcls()));
        }
        histos.at(9)->Fill(track->GetITSNcls());
        histos.at(10)->Fill(track->Chi2perNDF());
        histos.at(11)->Fill(DCAxy);
        histos.at(12)->Fill(DCAz);
        histos.at(13)->Fill(DCAxy,DCAz);

        // fill PID histos
        histos.at(14)->Fill(trackP,track->GetTPCsignal());
        histos.at(15)->Fill(trackP,nSigmaTPCprot);
        histos.at(16)->Fill(trackP,nSigmaTPCdeut);
        histos.at(17)->Fill(trackP,nSigmaTPCHe3);
        histos.at(18)->Fill(trackP,fbetaTOF);
        histos.at(19)->Fill(trackP,nSigmaTOFprot);
        histos.at(20)->Fill(trackP,nSigmaTOFdeut);
        histos.at(21)->Fill(trackP,nSigmaTOFpion);
        ((TH3F*)(histos.at(22)))->Fill(trackP,fmass2TOF,DCAxy);
        histos.at(23)->Fill(trackP,track->GetITSsignal());
        histos.at(24)->Fill(trackP,nSigmaITSprot);
        histos.at(25)->Fill(trackP,nSigmaITSdeut);

    }

}

//________________________________________________________________________
Float_t AliAnalysisTaskAntipd::GetTOFBeta(AliAODTrack *track)
{
    float beta = -999;
    double integratedTimes[9] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};

    track->GetIntegratedTimes(integratedTimes);

    const float c = 2.99792457999999984e-02;
    float p = track->P();
    float l = integratedTimes[0]*c;

    float trackT0 = fPIDResponse->GetTOFResponse().GetStartTime(p);

    float timeTOF = track->GetTOFsignal() - trackT0;
    if(timeTOF > 0){
        beta  = l/timeTOF/c;
    }
    return beta;
}

//________________________________________________________________________
Float_t AliAnalysisTaskAntipd::GetMass2TOF(Float_t beta, AliAODTrack *track)
{
    Float_t p = track->P();
    Float_t mass2sq = -999;
    if(!(beta==0)){
        mass2sq = ((1/(beta*beta))-1)*(p*p);
    }
    return mass2sq;
}
