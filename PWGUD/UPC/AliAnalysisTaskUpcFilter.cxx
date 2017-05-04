
//_____________________________________________________________________________
//    Class for semiforward UPC filter
//    Author: Jaroslav Adam
//
//    Fill structure of AliUPCEvent
//_____________________________________________________________________________

// c++ headers
#include <iostream>
#include <string.h>

// root headers
#include "TList.h"
#include "TH1I.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TObjString.h"
#include "TFile.h"
#include "TArrayF.h"

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliAODVertex.h"
#include "AliESDVZERO.h"
#include "AliESDAD.h"
#include "AliAODAD.h"
#include "AliVAD.h"
#include "AliESDZDC.h"
#include "AliESDVertex.h"
#include "AliAODTrack.h"
#include "AliMultiplicity.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliTriggerAnalysis.h"
#include "AliMuonTrackCuts.h"
#include "AliESDtrackCuts.h"
#include "AliAODPid.h"
#include "AliGenEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliPIDResponse.h"
#include "AliESDHLTtrack.h"

// my headers
#include "AliUPCTrack.h"
#include "AliUPCMuonTrack.h"
#include "AliUPCEvent.h"
#include "AliAnalysisTaskUpcFilter.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskUpcFilter);

// task for upc semiforward filter
// jaroslav.adam@cern.ch

//_____________________________________________________________________________
AliAnalysisTaskUpcFilter::AliAnalysisTaskUpcFilter(const char *name)
 :AliAnalysisTaskSE(name),
  fIsESD(0), fIsMC(0), fFillSPD(0), fMuonCuts(0x0), fTriggerAna(0x0), fCutsList(0x0), fPIDResponse(0x0), fMuonCutsPassName(0x0),
  fHistList(0x0), fCounter(0x0), fTriggerCounter(0x0), fMuonCounter(0x0),
  fUPCEvent(0x0), fUPCTree(0x0)
{

  // Constructor
  for(Int_t itrg=0; itrg<fgkNtrg; itrg++) fTrgMask[itrg] = kTRUE;

  fMuonCutsPassName = new TObjString();
  fMuonCutsPassName->SetString("muon_calo_pass1");

  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());

}//AliAnalysisTaskUpcFilter

//_____________________________________________________________________________
AliAnalysisTaskUpcFilter::~AliAnalysisTaskUpcFilter()
{
  // destructor

  if(fHistList) {delete fHistList; fHistList = 0x0;}
  if(fCounter) {delete fCounter; fCounter = 0x0;}
  if(fTriggerCounter) {delete fTriggerCounter; fTriggerCounter = 0x0;}
  if(fMuonCounter) {delete fMuonCounter; fMuonCounter = 0x0;}
  if(fUPCEvent) {delete fUPCEvent; fUPCEvent = 0x0;}
  if(fUPCTree) {delete fUPCTree; fUPCTree = 0x0;}
  if(fMuonCuts) {delete fMuonCuts; fMuonCuts = 0x0;}
  if(fTriggerAna) {delete fTriggerAna; fTriggerAna = 0x0;}
  if(fCutsList) {delete[] fCutsList; fCutsList=0x0;}
  if(fPIDResponse) {delete fPIDResponse; fPIDResponse=0x0;}
  if(fMuonCutsPassName) {delete fMuonCutsPassName; fMuonCutsPassName=0x0;}

}//~AliAnalysisTaskUpcFilter

//_____________________________________________________________________________
void AliAnalysisTaskUpcFilter::SetAllTrg(Bool_t set)
{
  //activate or deactivate all trigger classes by the mask 'fTrgMask'
  //
  //put kFALSE to mask all triggers and kTRUE to unmask
  //
  //by defalut all triggers are active

  for(Int_t itrg=0; itrg<fgkNtrg; itrg++) fTrgMask[itrg] = set;

}//SetAllTrg

//_____________________________________________________________________________
void AliAnalysisTaskUpcFilter::SetTrgClass(Int_t idx, Bool_t set)
{
  //activate or deactivate the trigger classe at index 'idx' by putting the value of 'set'
  //
  //put kFALSE to mask the trigger and kTRUE to unmask

  fTrgMask[idx] = set;

}//SetTrgClass

//_____________________________________________________________________________
void AliAnalysisTaskUpcFilter::UserCreateOutputObjects()
{
  //muon track cuts
  fMuonCuts = new AliMuonTrackCuts("StdMuonCuts", "StdMuonCuts");
  fMuonCuts->SetFilterMask ( AliMuonTrackCuts::kMuPdca );
  fMuonCuts->Print("mask");
  fMuonCuts->SetAllowDefaultParams(kTRUE);
  fMuonCuts->SetPassName( (fMuonCutsPassName->GetString()).Data() );

  //trigger analysis for SPD FO fired chips in MC
  fTriggerAna = new AliTriggerAnalysis();
  fTriggerAna->SetAnalyzeMC( fIsMC );

  //PID response, input handler from AliAnalysisTaskSE
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) ( AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler() );
  fPIDResponse = inputHandler->GetPIDResponse();

  //ESD track cuts
  if( fIsESD ) {
    fCutsList = new AliESDtrackCuts*[32];
    for(UInt_t i=0; i<32; i++) fCutsList[i] = 0x0;

    // bits 0 - 10 set in $ALICE_ROOT/ANALYSIS/ESDfilter/macros/AddTaskESDFilter.C

    // Cuts on primary tracks
    AliESDtrackCuts* esdTrackCutsL = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fCutsList[0] = esdTrackCutsL;

    // ITS stand-alone tracks
    AliESDtrackCuts* esdTrackCutsITSsa = new AliESDtrackCuts("ITS stand-alone Track Cuts", "ESD Track Cuts");
    esdTrackCutsITSsa->SetRequireITSStandAlone(kTRUE);
    fCutsList[1] = esdTrackCutsITSsa;

    // standard cuts with very loose DCA
    AliESDtrackCuts* esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); 
    esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
    esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
    esdTrackCutsH->SetDCAToVertex2D(kTRUE);
    fCutsList[4] = esdTrackCutsH;

    // standard cuts with tight DCA cut
    AliESDtrackCuts *esdTrackCutsH2 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
    fCutsList[5] = esdTrackCutsH2;

    // bits 11 - 31 free

    // selection by Daniele, UPC meeting, 15 July 2014
    // https://indico.cern.ch/event/330483/contribution/3/material/slides/0.pdf

    AliESDtrackCuts *cuts11 = new AliESDtrackCuts();
    cuts11->SetMinNClustersTPC(70);
    cuts11->SetRequireTPCRefit(kTRUE);
    cuts11->SetRequireITSRefit(kTRUE);
    cuts11->SetMaxChi2PerClusterTPC(4);
    cuts11->SetAcceptKinkDaughters(kFALSE);
    cuts11->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    cuts11->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.01");
    cuts11->SetMaxDCAToVertexZ(2);
    cuts11->SetDCAToVertex2D(kFALSE);
    cuts11->SetRequireSigmaToVertex(kFALSE);
    cuts11->SetMaxChi2PerClusterITS(5.);
    fCutsList[11] = cuts11;

    AliESDtrackCuts *cuts12 = new AliESDtrackCuts();
    cuts12->SetMinNClustersTPC(70);
    cuts12->SetRequireTPCRefit(kTRUE);
    cuts12->SetRequireITSRefit(kTRUE);
    cuts12->SetMaxChi2PerClusterTPC(4);
    cuts12->SetAcceptKinkDaughters(kFALSE);
    cuts12->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    cuts12->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.01");
    cuts12->SetMaxDCAToVertexZ(2);
    cuts12->SetDCAToVertex2D(kFALSE);
    cuts12->SetRequireSigmaToVertex(kFALSE);
    fCutsList[12] = cuts12;

    // setting in Evgeny's trees
    AliESDtrackCuts *cuts13 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
    cuts13->SetMaxChi2TPCConstrainedGlobal(1e+10);
    fCutsList[13] = cuts13;

    // ITS stand-alone tracks
    AliESDtrackCuts* esdTrackCutsPureITSsa = new AliESDtrackCuts("ITS pure stand-alone Track Cuts", "ESD Track Cuts");
    esdTrackCutsPureITSsa->SetRequireITSPureStandAlone(kTRUE);
    fCutsList[14] = esdTrackCutsPureITSsa;

  }

  //output histograms
  fHistList = new TList();
  fHistList->SetOwner();

  fCounter = new TH1I("fCounter", "fCounter", 100, 1, 101);
  fHistList->Add(fCounter);

  fTriggerCounter = new TH2I("fTriggerCounter", "fTriggerCounter", 126000, 154000, 280000, fgkNtrg+1, 0, fgkNtrg+1);
  fHistList->Add(fTriggerCounter);

  fMuonCounter = new TH1I("fMuonCounter", "fMuonCounter", 4, 1, 5);
  fHistList->Add(fMuonCounter);

  //errors in counter start at 51
  if( !fPIDResponse ) fCounter->Fill( kPidErr ); // no PID response

  //output tree
  fUPCEvent  = new AliUPCEvent();
  fUPCEvent->SetIsMC( fIsMC );
  fUPCEvent->SetIsESD( fIsESD );

  TDirectory *pwd = gDirectory;
  OpenFile(1);
  fUPCTree = new TTree("fUPCTree", "fUPCTree");
  pwd->cd();
  fUPCTree->Branch("fUPCEvent", &fUPCEvent);

  PostData(1, fUPCTree);
  PostData(2, fHistList);

}//UserCreateOutputObjects

//_____________________________________________________________________________
void AliAnalysisTaskUpcFilter::NotifyRun()
{

  fMuonCuts->SetRun(fInputHandler); // use input handler from AliAnalysisTaskSE

  if( fIsESD ) { ((AliESDEvent*) InputEvent())->InitMagneticField(); }

}//NotifyRun

//_____________________________________________________________________________
void AliAnalysisTaskUpcFilter::UserExec(Option_t *) 
{

  // input event
  AliVEvent *vEvent = dynamic_cast<AliVEvent*>(InputEvent());
  if(!vEvent) return;

  fUPCEvent->ClearEvent();

  fCounter->Fill( kAna ); // analyzed events

  // trigger
  TString trigger = vEvent->GetFiredTriggerClasses();
  Bool_t trgClasses[fgkNtrg]; // array of fired trigger classes
  for(Int_t itrg=0; itrg<fgkNtrg; itrg++) trgClasses[itrg] = kFALSE;

  //list of trigger classes

  trgClasses[ 1] = trigger.Contains("CMUP6-B"); // p-Pb FW
  trgClasses[ 2] = trigger.Contains("CMUP3-B"); // Pb-p FW
  trgClasses[ 3] = trigger.Contains("CMUP8-B"); // Pb-p FW

  trgClasses[ 4] = trigger.Contains("CMUP7-B"); // p-Pb SFW
  trgClasses[ 5] = trigger.Contains("CMUP5-B"); // Pb-p SFW
  trgClasses[ 6] = trigger.Contains("CMUP9-B"); // Pb-p SFW

  trgClasses[ 7] = trigger.Contains("CCUP7-B"); // p-Pb Pb-p CEN

  trgClasses[ 8]= trigger.Contains("CMUP1-B"); // PbPb FW  !0VBA & 0VBC & 0MSL

  trgClasses[ 9]= trigger.Contains("CTRUE-B"); // p-Pb control trigger

  trgClasses[10] = trigger.Contains("CCUP2-B");      // !0VBA & !0VBC & 0SH1 & 0OM2
  trgClasses[11] = trigger.Contains("CCUP4-B");      // !0VBA & !0VBC & 0SH1 & 0OMU
  trgClasses[12] = trigger.Contains("CCUP8-B");     // *0VBA *0VBC *0UBA *0UBC 0STP 0OMU              (=CTEST57-B)
  trgClasses[13] = trigger.Contains("CCUP9-B");     // *0VBA *0VBC *0UBA *0UBC 0STP                   (=CTEST59-B)
  trgClasses[14] = trigger.Contains("CCUP10-B");    // *0VBA *0VBC *0UBA *0UBC 0SH1                   (=CTEST58-B)

  trgClasses[15] = trigger.Contains("CMUP10-B");    // *0VBA *0UBA *0UBC 0MSL                         (=CTEST63-B)
  trgClasses[16] = trigger.Contains("CMUP11-B");    // !0VBA & !0UBA & !0UBC & 0MUL                   (=CTEST64-B)
  trgClasses[17] = trigger.Contains("CMUP12-B");    // !0VBA & !0UBA & !0UBC & 0MSL & 0SMB

  trgClasses[18] = trigger.Contains("CTEST62-B");    // !0VBA & !0UBA & !0UBC & 0VBC & 0MSL
  trgClasses[19] = trigger.Contains("CTEST63-B");    // !0VBA & !0UBA & !0UBC & 0MSL
  trgClasses[20] = trigger.Contains("CTEST64-B");    // !0VBA & !0UBA & !0UBC & 0MUL

  trgClasses[21] = trigger.Contains("CTEST57-B");    // !0VBA & !0VBC & !0UBA & !0UBC & 0STP & 0OMU
  trgClasses[22] = trigger.Contains("CTEST58-B");    // !0VBA & !0VBC & !0UBA & !0UBC & 0SH1
  trgClasses[23] = trigger.Contains("CTEST59-B");   // !0VBA & !0VBC & !0UBA & !0UBC & 0STP
  trgClasses[24] = trigger.Contains("CTEST60-B");   // !0VBA & !0VBC & !0UBA & !0UBC & 0OM2
  trgClasses[25] = trigger.Contains("CTEST61-B");   // !0VBA & !0VBC & !0UBA & !0UBC & 0OMU
  
  trgClasses[26] = trigger.Contains("CMUP14-B");   // 0MSL & !0VBA & !0UBA
  trgClasses[27] = trigger.Contains("CMUP15-B");   // *0VBA *0UBA *0VC5 0SMB *0SH2 0MSL
  trgClasses[28] = trigger.Contains("CMUP16-B");   // 0MSL *0VBA *0UBA *0UGC *0VGA
  trgClasses[29] = trigger.Contains("CMUP17-B");   // *0VBA *0UBA *0VC5 0SMB *0SH2 0MSL *0UGC *0VGA
  trgClasses[30] = trigger.Contains("CMUP21-B");   // *0VBA *0UBA *0VBC 0SH1 *0SH2 *0UGC *0VGA
  trgClasses[31] = trigger.Contains("CMUP22-B");   // *0UBC *0UGC *0VBA *0VGA *0SH2 *0VC5 0MSL 0SMB
  trgClasses[32] = trigger.Contains("CMUP23-B");   // *0UBC *0UGC *0VBA *0VGA *0SH2 *0VC5 0MUL
  
  trgClasses[33] = trigger.Contains("CCUP14-B");   // *0VBA *0UBA *0VC5 0OMU 0STG
  trgClasses[34] = trigger.Contains("CCUP15-B");   // *0VBA *0UBA *0VC5 0SH1
  trgClasses[35] = trigger.Contains("CCUP16-B");   // *0VBA *0UBA *0VC5 0OMU 0STG *0UGC *0VGA
  trgClasses[36] = trigger.Contains("CCUP17-B");   // *0VBA *0UBA *0VC5 0SH1 *0UGC *0VGA
  trgClasses[37] = trigger.Contains("CCUP20-B");   // *0VBA *0UBA *0VBC 0OMU 0STG *0SH2 *0UGC *0VG
  trgClasses[38] = trigger.Contains("CCUP21-B");   // *0VBA *0UBA *0VBC 0SH1 *0SH2 *0UGC *0VGA
  trgClasses[39] = trigger.Contains("CCUP22-B");   // *0UBC *0UGC *0VBA *0VGA *0SH2 *0VBC 0STG 0OMU
  trgClasses[40] = trigger.Contains("CCUP23-B");   // *0UBC *0UGC *0VBA *0VGA *0SH2 *0VBC 0SH1

  //end of list of trigger classes

  Bool_t isTrg = kFALSE;
  for(Int_t itrg=1; itrg<fgkNtrg; itrg++) {
    if( !trgClasses[itrg] || !fTrgMask[itrg] ) continue;

    //unmasked trigger at itrg is fired
    fUPCEvent->SetTriggerClass( itrg , kTRUE );
    fTriggerCounter->Fill( vEvent->GetRunNumber() , itrg );
    isTrg = kTRUE;
  }
  if(!isTrg && !fIsMC) {PostData(2, fHistList); return;}
  //event passed the trigger

  fCounter->Fill( kTrg ); //  events after trigger (ESD and AOD)

  // ESD / AOD specific tasks: MC, SPD FO, L0 inputs, SPD tracklets, tracks, ZDC tdc
  Bool_t stat;
  if( fIsESD ) stat = RunESD();
  if(!fIsESD ) stat = RunAOD();
  if( !stat && !fIsMC ) {PostData(2, fHistList); return;}
  //event passed ESD / AOD specific selection

  fCounter->Fill( kSpecific ); // events after ESD / AOD specific part

  //L1 trigger inputs for 1ZED, idx = 14
  UInt_t inputsL1 = vEvent->GetHeader()->GetL1TriggerInputs();
  Bool_t inp1ZED = inputsL1 & (1<<14);
  //use UPC event internal flag bits to save 1ZED
  fUPCEvent->ResetFlagBits();     //first put all flag bits to 0
  fUPCEvent->SetIsESD( fIsESD );  //set bits 0 and 1 for data type and mc
  fUPCEvent->SetIsMC( fIsMC );
  if(inp1ZED) fUPCEvent->SetFlagBit( (UChar_t) 2 ); // use bit 2 for 1ZED, bit is set when 1ZED is fired, remains cleared otherwise

  // input data
  const char *filnam = ((TTree*) GetInputData(0))->GetCurrentFile()->GetName();
  // reconstruction pass: -1 = unknown, 1 = pass1, 2 = pass2
  fUPCEvent->SetRecoPass( -1 );
  if( strstr(filnam,"/muon_calo_pass1/") ) {fUPCEvent->SetRecoPass( 1 ); fCounter->Fill( kPass1 );}
  if( strstr(filnam,"/pass2/") ) {fUPCEvent->SetRecoPass( 2 ); fCounter->Fill( kPass2 );}
  if( fUPCEvent->GetRecoPass() < 0 ) fCounter->Fill( kPassX );
  fUPCEvent->SetInputFileName( filnam );
  fUPCEvent->SetEventNumber( ((TTree*) GetInputData(0))->GetTree()->GetReadEntry() );
  fUPCEvent->SetRunNumber( vEvent->GetRunNumber() );

  //VZERO
  AliVVZERO *dataVZERO = dynamic_cast<AliVVZERO*>(vEvent->GetVZEROData());
  if(!dataVZERO) {PostData(2, fHistList); return;}

  fUPCEvent->SetV0ADecision( dataVZERO->GetV0ADecision() );
  fUPCEvent->SetV0CDecision( dataVZERO->GetV0CDecision() );
  for (UInt_t iv=0; iv<32; iv++) {
    if( dataVZERO->BBTriggerV0C((Int_t)iv) ) fUPCEvent->SetBBtriggerV0Cmask(iv);
    if( dataVZERO->GetBBFlag((Int_t)iv) ) fUPCEvent->SetBBFlagV0Cmask(iv);
    if( dataVZERO->BBTriggerV0A((Int_t)iv) ) fUPCEvent->SetBBtriggerV0Amask(iv);
    if( dataVZERO->GetBBFlag((Int_t)iv+32) ) fUPCEvent->SetBBFlagV0Amask(iv);
  }

  //AD
  AliVAD *dataAD = dynamic_cast<AliVAD*>(vEvent->GetADData());
  if(dataAD) {
    fUPCEvent->SetADADecision( dataAD->GetADADecision() );
    fUPCEvent->SetADCDecision( dataAD->GetADCDecision() );
    for (UInt_t iv=0; iv<8; iv++) {
      if( dataAD->BBTriggerADC((Int_t)iv) ) fUPCEvent->SetBBtriggerADCmask(iv);
      if( dataAD->GetBBFlag((Int_t)iv) ) fUPCEvent->SetBBFlagADCmask(iv);
      if( dataAD->BBTriggerADA((Int_t)iv) ) fUPCEvent->SetBBtriggerADAmask(iv);
      if( dataAD->GetBBFlag((Int_t)iv+8) ) fUPCEvent->SetBBFlagADAmask(iv);
  }
  } else {
    fUPCEvent->SetADADecision( -999 );
    fUPCEvent->SetADCDecision( -999 );
  }

  //ZDC
  AliVZDC *dataZDC = dynamic_cast<AliVZDC*>(vEvent->GetZDCData());
  if(!dataZDC) {PostData(2, fHistList); return;}

  //energy in ZDC
  Double_t eZnc=0., eZpc=0., eZna=0., eZpa=0.;
  for(Int_t i=0; i<5; i++) {
    eZnc += dataZDC->GetZNCTowerEnergy()[i];
    eZpc += dataZDC->GetZPCTowerEnergy()[i];
    eZna += dataZDC->GetZNATowerEnergy()[i];
    eZpa += dataZDC->GetZPATowerEnergy()[i];
  }
  fUPCEvent->SetZNCEnergy( eZnc );
  fUPCEvent->SetZPCEnergy( eZpc );
  fUPCEvent->SetZNAEnergy( eZna );
  fUPCEvent->SetZPAEnergy( eZpa );

  //default primary vertex
  const AliVVertex *vtx = vEvent->GetPrimaryVertex();
  if(!vtx) {PostData(2, fHistList); return;}
  Double_t pos[3]; vtx->GetXYZ(pos);
  Double_t cov[6]; vtx->GetCovarianceMatrix(cov);
  fUPCEvent->SetPrimaryVertex(pos, vtx->GetChi2perNDF(), cov, vtx->GetNContributors());
  const char *vtxtitle = vtx->GetTitle();
  fUPCEvent->SetPrimaryVertexTitle( vtxtitle );

  //fUPCEvent->MakeArrayD(1);
  //fUPCEvent->GetArrayD()->SetAt( eZnc, 0 ); // example of arrayD usage

  fCounter->Fill( kWritten ); // events written to the tree (ESD and AOD)

  fUPCTree ->Fill();

  PostData(1, fUPCTree);
  PostData(2, fHistList);

}//UserExec

//_____________________________________________________________________________
Bool_t AliAnalysisTaskUpcFilter::RunAOD()
{

  //input AOD event
  AliAODEvent *aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!aodEvent) return kFALSE;

  fCounter->Fill( kAOD ); // AOD analyzed events

  if(fIsMC) RunAODMC( (TClonesArray*) aodEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName()),
                      (AliAODMCHeader*) aodEvent->FindListObject(AliAODMCHeader::StdBranchName())
                    );

  //nominal interaction point
  static AliAODVertex *vtx0 = 0x0;
  if( !vtx0 ) {
    // make a vertex at coordinates 0, 0, 0
    vtx0 = new AliAODVertex();
    vtx0->SetX(0.); vtx0->SetY(0.); vtx0->SetZ(0.);
  }

  Double_t pxpypz[3];
  UChar_t maskMan;
  Double_t xyzDca[2], cov[3];
  Float_t b[2], covF[3];
  Int_t nmun=0, ncen=0;
  Bool_t pdca;
  // AOD tracks loop
  for(Int_t itr=0; itr<aodEvent->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(itr));
    if( !trk ) continue;

    //muon track
    if( trk->IsMuonTrack() ) {

      fMuonCounter->Fill( kMunAll );

      //select only tracks with good eta and Rabs
      if( trk->GetRAtAbsorberEnd() < 17.5 || trk->GetRAtAbsorberEnd() > 89.5 ) continue;
      fMuonCounter->Fill( kMunRabs );
      if( trk->Eta() < -4.0 || trk->Eta() > -2.5 ) continue;
      fMuonCounter->Fill( kMunEta );

      pdca = fMuonCuts->IsSelected(trk);
      if(!pdca) continue;
      fMuonCounter->Fill( kMunPDCA );

      //muon track is accepted to put to the output

      AliUPCMuonTrack *upcMuon = fUPCEvent->AddMuonTrack();
      fCounter->Fill( kMunTrack ); // added muon track
      upcMuon->SetPtEtaPhi( trk->Pt(), trk->Eta(), trk->Phi() );
      upcMuon->SetCharge( trk->Charge() );
      upcMuon->SetMatchTrigger( trk->GetMatchTrigger() );
      upcMuon->SetRAtAbsorberEnd( trk->GetRAtAbsorberEnd() );
      upcMuon->SetChi2perNDF( trk->Chi2perNDF() );
      upcMuon->SetDCA( trk->DCA() );
      upcMuon->SetPxDCA( pdca );

      nmun++;
    }
    //central track
    else {

      maskMan = 0;
      if( trk->HasPointOnITSLayer(0) )                 { maskMan |= 1 << 0; }
      if( trk->HasPointOnITSLayer(1) )                 { maskMan |= 1 << 1; }
      if( trk->GetStatus() & AliESDtrack::kITSrefit )  { maskMan |= 1 << 2; }
      if( trk->GetStatus() & AliESDtrack::kTPCrefit )  { maskMan |= 1 << 3; }
      if( trk->GetTOFsignal() < 99999. && trk->GetTOFsignal() > 0. )  { maskMan |= 1 << 4; }

      //selection for at least one point in SPD
      //if( !(maskMan & (1 << 0)) && !(maskMan & (1 << 1)) ) continue;
	/*/
      trk->GetPxPyPz(pxpypz);
      AliAODPid *apid = trk->GetDetPid();

      AliUPCTrack *upcTrack = fUPCEvent->AddTrack();
      fCounter->Fill( kCenTrack ); // added central track
      upcTrack->SetPxPyPz( pxpypz );
      upcTrack->SetMaskMan( maskMan );
      upcTrack->SetFilterMap( trk->GetFilterMap() );
      upcTrack->SetCharge( trk->Charge() );
      upcTrack->SetChi2perNDF( trk->Chi2perNDF() );
      if(apid) {
        upcTrack->SetTPCmomentum( apid->GetTPCmomentum() );
        upcTrack->SetTPCsignal( apid->GetTPCsignal() );
      }
      upcTrack->SetTPCNcls( trk->GetTPCNcls() );
      upcTrack->SetTPCCrossedRows( (Float_t) trk->GetTPCNCrossedRows() );
      upcTrack->SetTPCNclsF( trk->GetTPCNclsF() );
      upcTrack->SetTPCNclsS( trk->GetTPCnclsS() );
      upcTrack->SetITSClusterMap( trk->GetITSClusterMap() );
      upcTrack->SetTOFsignal( trk->GetTOFsignal() );

      //DCA to default primary vertex
      AliAODTrack *track1 = (AliAODTrack*) trk->Clone("track1");
      track1->PropagateToDCA(aodEvent->GetPrimaryVertex(), aodEvent->GetMagneticField(), 9999., xyzDca, cov);
      for(Int_t i=0; i<2; i++) {b[i] = (Float_t) xyzDca[i]; covF[i] = (Float_t) cov[i];}
      upcTrack->SetImpactParameters(b, covF);
      delete track1;

      //DCA to SPD vertex
      track1 = (AliAODTrack*) trk->Clone("track1");
      track1->PropagateToDCA(aodEvent->GetPrimaryVertexSPD(), aodEvent->GetMagneticField(), 9999., xyzDca, cov);
      for(Int_t i=0; i<2; i++) {b[i] = (Float_t) xyzDca[i]; covF[i] = (Float_t) cov[i];}
      upcTrack->SetImpactParametersSPD(b, covF);
      delete track1;

      //DCA to nominal interaction point
      track1 = (AliAODTrack*) trk->Clone("track1");
      track1->PropagateToDCA(vtx0, aodEvent->GetMagneticField(), 9999., xyzDca, cov);
      for(Int_t i=0; i<2; i++) {b[i] = (Float_t) xyzDca[i]; covF[i] = (Float_t) cov[i];}
      upcTrack->SetImpactParametersIP(b, covF);
      delete track1;
      /*/

      ncen++;
    }
  }// AOD tracks loop

  //selection for at least one muon track
  if( nmun < 1 ) return kFALSE;

  //selection for at least one muon and at least one central track
  //if( nmun < 1 || ncen < 1 ) return kFALSE;


  // L0 trigger inputs
  fUPCEvent->SetL0inputs( aodEvent->GetHeader()->GetL0TriggerInputs() );

  // Tracklets
  fUPCEvent->SetNumberOfTracklets( aodEvent->GetTracklets()->GetNumberOfTracklets() );

  // AOD ZDC for timing
  AliAODZDC *dataZDCAOD = aodEvent->GetZDCData();
  if(!dataZDCAOD) {PostData(2, fHistList); return kFALSE;}

  fUPCEvent->SetZNCTime( dataZDCAOD->GetZNCTime() );
  fUPCEvent->SetZNATime( dataZDCAOD->GetZNATime() );

  //SPD primary vertex in AOD
  AliAODVertex *vtx = aodEvent->GetPrimaryVertexSPD();
  if(!vtx) {PostData(2, fHistList); return kFALSE;}
  Double_t posVtx[3]; vtx->GetXYZ( posVtx );
  Double_t covVtx[6]; vtx->GetCovarianceMatrix( covVtx );
  fUPCEvent->SetPrimaryVertexSPD(posVtx, vtx->GetChi2perNDF(), covVtx, vtx->GetNContributors());
  const char *vtxtitle = vtx->GetTitle();
  fUPCEvent->SetPrimaryVertexSPDtitle( vtxtitle );

  return kTRUE;

}//RunAOD

//_____________________________________________________________________________
void AliAnalysisTaskUpcFilter::RunAODMC(TClonesArray *arrayMC, AliAODMCHeader *headerMC)
{
  // run over AOD mc particles

  if( !arrayMC || !headerMC ) return;

  //generated vertex in AOD MC
  Double_t vtxD[3]; headerMC->GetVertex( vtxD );
  Float_t vtxF[3]; for(Int_t i=0; i<3; i++) vtxF[i] = (Float_t) vtxD[i];
  fUPCEvent->SetPrimaryVertexMC( vtxF );

  //loop over mc particles
  for(Int_t imc=0; imc<arrayMC->GetEntriesFast(); imc++) {
    AliAODMCParticle *aodmc = dynamic_cast<AliAODMCParticle*>(arrayMC->At(imc));
    if(!aodmc) continue;

    //if(aodmc->GetMother() >= 0) continue;

    TParticle *part = fUPCEvent->AddMCParticle();
    part->SetMomentum(aodmc->Px(), aodmc->Py(), aodmc->Pz(), aodmc->E());
    part->SetProductionVertex(aodmc->Xv(), aodmc->Yv(), aodmc->Zv(), aodmc->T());
    part->SetFirstMother(aodmc->GetMother());
    part->SetLastDaughter(aodmc->GetNDaughters());
    part->SetPdgCode(aodmc->GetPdgCode());
    part->SetUniqueID(imc);
  }//loop over mc particles

}//RunAODMC

//_____________________________________________________________________________
Bool_t AliAnalysisTaskUpcFilter::RunESD()
{

  // Input ESD event
  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!esdEvent) return kFALSE;

  fCounter->Fill( kESD ); // ESD analyzed events

  if( fFillSPD || fIsMC ) {
    //SPD FO fired chips
    fUPCEvent->SetNSPDfiredInner( fTriggerAna->SPDFiredChips(esdEvent,1,kFALSE,1) );
    fUPCEvent->SetNSPDfiredOuter( fTriggerAna->SPDFiredChips(esdEvent,1,kFALSE,2) );
    fUPCEvent->SetFastOrFiredChips( (TBits*) &esdEvent->GetMultiplicity()->GetFastOrFiredChips() );
  }

  if(fIsMC) RunESDMC();

  static AliESDVertex *vtx0 = 0x0;
  if( !vtx0 ) {
    // make a vertex at coordinates 0, 0, 0
    vtx0 = new AliESDVertex();
    vtx0->SetXv(0.); vtx0->SetYv(0.); vtx0->SetZv(0.);
  }

  // ESD central tracks
  Double_t pxpypz[3];
  Float_t b[2]; Float_t cov[3];
  UChar_t maskMan;
  UInt_t filterMap;
  Int_t nmun=0, ncen=0;

  //ESD central tracks loop
  for(Int_t itr=0; itr<esdEvent->GetNumberOfTracks(); itr++) {
    AliESDtrack *eTrack = esdEvent->GetTrack(itr);
    if( !eTrack ) continue;

    // manual filter
    maskMan = 0;
    if( eTrack->HasPointOnITSLayer(0) )                 { maskMan |= 1 << 0; }
    if( eTrack->HasPointOnITSLayer(1) )                 { maskMan |= 1 << 1; }
    if( eTrack->GetStatus() & AliESDtrack::kITSrefit )  { maskMan |= 1 << 2; }
    if( eTrack->GetStatus() & AliESDtrack::kTPCrefit )  { maskMan |= 1 << 3; }
    if( eTrack->GetTOFsignal() < 99999. && eTrack->GetTOFsignal() > 0. )  { maskMan |= 1 << 4; }
    if( eTrack->GetKinkIndex(0) > 0 )                   { maskMan |= 1 << 5; } // bit set if track is kink candidate

    //selection for at least one point in SPD
    //if( !(maskMan & (1 << 0)) && !(maskMan & (1 << 1)) ) continue;

    //selection for ITS refit
    //if( !(maskMan & (1 << 2)) ) continue;

    //take out the selections for reconstruction without the TPC, write all the tracks

    //central track accepted to write to the UPC event

    // various configurations of AliESDtrackCuts
    filterMap = 0;
    for(UInt_t i=0; i<32; i++) {
      if( !fCutsList[i] ) continue;

      if( fCutsList[i]->AcceptTrack(eTrack) ) { filterMap |= 1 << i; }
    }

    eTrack->GetPxPyPz(pxpypz);

    //fill the UPC track
    AliUPCTrack *upcTrack = fUPCEvent->AddTrack();
    fCounter->Fill( kCenTrack ); // added central track
    upcTrack->SetPxPyPz( pxpypz );
    upcTrack->SetMaskMan( maskMan );
    upcTrack->SetFilterMap( filterMap );
    upcTrack->SetCharge( eTrack->Charge() );
    upcTrack->SetChi2perNDF( eTrack->GetTPCchi2()/((Double_t) eTrack->GetTPCNcls()) ); // TPC chi2 used to fill this data member in esd
    upcTrack->SetTPCmomentum( eTrack->GetTPCmomentum() );
    upcTrack->SetTPCsignal( eTrack->GetTPCsignal() );
    upcTrack->SetTPCNcls( eTrack->GetTPCNcls() );
    upcTrack->SetTPCCrossedRows( eTrack->GetTPCCrossedRows() );
    upcTrack->SetTPCNclsF( eTrack->GetTPCNclsF() );
    upcTrack->SetTPCNclsS( eTrack->GetTPCnclsS() );
    upcTrack->SetITSchi2perNDF( eTrack->GetITSchi2()/((Double_t) eTrack->GetNcls(0)) );
    upcTrack->SetITSClusterMap( eTrack->GetITSClusterMap() );
    upcTrack->SetTOFsignal( eTrack->GetTOFsignal() );
    upcTrack->MakeArrayD(1);
    upcTrack->GetArrayD()->SetAt( eTrack->GetITSsignal() , 0 ); // put ITS signal at slot 0
    upcTrack->MakeArrayInt(2);
    upcTrack->GetArrayInt()->SetAt( (Int_t) eTrack->GetTPCNcls(), 0 ); // example of arrayI usage
    upcTrack->GetArrayInt()->SetAt( (Int_t) eTrack->GetTPCNclsF(), 1 );

    //TPC PID
    // AliPID::EParticleType = { kElectron = 0,  kMuon = 1,  kPion = 2,  kKaon = 3,  kProton = 4 .... }
    if( fPIDResponse && (maskMan & (1 << 3)) ) {
      for(Int_t itype=0; itype<=4; itype++) {
        upcTrack->SetNSigmasTPC( itype, fPIDResponse->NumberOfSigmasTPC( eTrack, (AliPID::EParticleType) itype ) );
      }
    }
    //TOF PID
    if( fPIDResponse && (maskMan & (1 << 4)) ) {
      for(Int_t itype=0; itype<=4; itype++) {
        upcTrack->SetNSigmasTOF( itype, fPIDResponse->NumberOfSigmasTOF( eTrack, (AliPID::EParticleType) itype ) );
      }
    }

    //DCA to default primary vertex
    eTrack->GetImpactParameters(b, cov);
    upcTrack->SetImpactParameters(b, cov);

    //DCA to SPD vertex
    AliESDtrack *track1 = (AliESDtrack*) eTrack->Clone("track1");
    track1->RelateToVertex( esdEvent->GetPrimaryVertexSPD(), esdEvent->GetMagneticField(), 9999. );
    track1->GetImpactParameters(b, cov);
    upcTrack->SetImpactParametersSPD(b, cov);
    delete track1;

    //DCA to nominal interaction point
    track1 = (AliESDtrack*) eTrack->Clone("track1");
    track1->RelateToVertex( vtx0, esdEvent->GetMagneticField(), 9999. );
    track1->GetImpactParameters(b, cov);
    upcTrack->SetImpactParametersIP(b, cov);
    delete track1;

    ncen++;
  } //ESD central tracks loop

  // ESD muon tracks
  //muon tracks loop
  for(Int_t itr=0; itr<esdEvent->GetNumberOfMuonTracks(); itr++) {
    AliESDMuonTrack *mTrack = esdEvent->GetMuonTrack(itr);
    if( !mTrack ) continue;

    //select only tracks with good eta and Rabs
    if( mTrack->Eta() < -4.0 || mTrack->Eta() > -2.5 ) continue;
    if( mTrack->GetRAtAbsorberEnd() < 17.5 || mTrack->GetRAtAbsorberEnd() > 89.5 ) continue;

    //muon track is accepted to put to the output

    AliUPCMuonTrack *upcMuon = fUPCEvent->AddMuonTrack();
    fCounter->Fill( kMunTrack ); // added muon track
    upcMuon->Clear();
    upcMuon->SetPtEtaPhi( mTrack->Pt(), mTrack->Eta(), mTrack->Phi() );
    upcMuon->SetCharge( mTrack->Charge() );
    upcMuon->SetMatchTrigger( mTrack->GetMatchTrigger() );
    upcMuon->SetRAtAbsorberEnd( mTrack->GetRAtAbsorberEnd() );
    upcMuon->SetChi2perNDF( mTrack->GetChi2()/((Double_t) mTrack->GetNDF()) );
    upcMuon->SetDCA( mTrack->GetDCA() );
    upcMuon->SetPxDCA( fMuonCuts->IsSelected(mTrack) );
    upcMuon->MakeArrayD(1);
    upcMuon->GetArrayD()->SetAt( mTrack->Eta(), 0 ); // example of arrayD usage

    nmun++;
  } //muon tracks loop

  //selection for at least one muon or central track
  //if( nmun + ncen < 1 ) return kFALSE;

  //selection for at least one muon and at least one central track
  //if( nmun < 1 || ncen < 1 ) return kFALSE;

  // L0 trigger inputs
  fUPCEvent->SetL0inputs( esdEvent->GetHeader()->GetL0TriggerInputs() );

  //Tracklets
  fUPCEvent->SetNumberOfTracklets( esdEvent->GetMultiplicity()->GetNumberOfTracklets() );

  // ESD ZDC for TDC
  AliESDZDC *dataZDCESD = esdEvent->GetESDZDC();
  if(!dataZDCESD) {PostData(2, fHistList); return kFALSE;}

  Bool_t znctdc = kFALSE, zpctdc = kFALSE, znatdc = kFALSE, zpatdc = kFALSE;
  Int_t tdcsum[4] = {0,0,0,0}, channeltdc[4];
  channeltdc[0] = dataZDCESD->GetZNCTDCChannel();
  channeltdc[1] = dataZDCESD->GetZPCTDCChannel();
  channeltdc[2] = dataZDCESD->GetZNATDCChannel();
  channeltdc[3] = dataZDCESD->GetZPATDCChannel();
  Float_t znctime = 0., znatime = 0.;
  for(Int_t iz=0;iz<4;iz++) {
    for(Int_t i=0; i<4; i++) tdcsum[i] += dataZDCESD->GetZDCTDCData(channeltdc[i],iz);
    znctime += dataZDCESD->GetZDCTDCCorrected(channeltdc[0],iz);
    znatime += dataZDCESD->GetZDCTDCCorrected(channeltdc[2],iz);
  }
  if( tdcsum[0] != 0 ) znctdc = kTRUE;
  if( tdcsum[1] != 0 ) zpctdc = kTRUE;
  if( tdcsum[2] != 0 ) znatdc = kTRUE;
  if( tdcsum[3] != 0 ) zpatdc = kTRUE;

  fUPCEvent->SetZNCtdc( znctdc );
  fUPCEvent->SetZPCtdc( zpctdc );
  fUPCEvent->SetZNAtdc( znatdc );
  fUPCEvent->SetZPAtdc( zpatdc );
  fUPCEvent->SetZNCtdcData( tdcsum[0] );
  fUPCEvent->SetZPCtdcData( tdcsum[1] );
  fUPCEvent->SetZNAtdcData( tdcsum[2] );
  fUPCEvent->SetZPAtdcData( tdcsum[3] );
  fUPCEvent->SetZNCTime( znctime );
  fUPCEvent->SetZNATime( znatime );

  //SPD primary vertex in ESD
  const AliESDVertex *vtx = esdEvent->GetPrimaryVertexSPD();
  if(!vtx) {PostData(2, fHistList); return kFALSE;}
  Double_t posVtx[3]; vtx->GetXYZ( posVtx );
  Double_t covVtx[6]; vtx->GetCovarianceMatrix( covVtx );
  fUPCEvent->SetPrimaryVertexSPD(posVtx, vtx->GetChi2perNDF(), covVtx, vtx->GetNContributors());
  const char *vtxtitle = vtx->GetTitle();
  fUPCEvent->SetPrimaryVertexSPDtitle( vtxtitle );

  return kTRUE;

}//RunESD

//_____________________________________________________________________________
void AliAnalysisTaskUpcFilter::RunESDMC()
{
  // ESD MC particles

  AliMCEvent *mcEvent = MCEvent();
  if(!mcEvent) {cout<<"Error: no ESD MC found"<<endl; return;}

  //generated vertex
  TArrayF vtx(3);
  mcEvent->GenEventHeader()->PrimaryVertex(vtx);
  fUPCEvent->SetPrimaryVertexMC( vtx.GetArray() );

  //loop over mc particles
  for(Int_t imc=0; imc<mcEvent->GetNumberOfTracks(); imc++) {
    AliMCParticle *esdmc = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(imc));
    if(!esdmc) continue;

    if(esdmc->GetMother() >= 0) continue;

    TParticle *part = fUPCEvent->AddMCParticle();
    part->SetMomentum(esdmc->Px(), esdmc->Py(), esdmc->Pz(), esdmc->E());
    part->SetProductionVertex(esdmc->Xv(), esdmc->Yv(), esdmc->Zv(), 0.);
    part->SetFirstMother(esdmc->GetMother());
    part->SetLastDaughter(esdmc->GetLastDaughter()-esdmc->GetFirstDaughter()+1);
    part->SetPdgCode(esdmc->PdgCode());
    part->SetUniqueID(imc);
  }//loop over mc particles

}//RunESDMC

//_____________________________________________________________________________
void AliAnalysisTaskUpcFilter::Terminate(Option_t *) 
{

  cout<<"Analysis complete."<<endl;
}//Terminate






























