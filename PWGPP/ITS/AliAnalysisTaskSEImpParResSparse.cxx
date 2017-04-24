/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for the study of the impact parameter resolution
//
// Authors:A.Dainese,    andrea.dainese@pd.infn.it
//     and Xianbao Yuan, yuanxb@iopp.ccnu.edu.cn; xianbao.yuan@pd.infn.it
/////////////////////////////////////////////////////////

#include <TList.h>
#include <TH1F.h>
#include <THnSparse.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "TParticle.h"
#include "AliGeomManager.h"
#include "AliMultiplicity.h"
#include "AliTriggerClass.h"
#include "AliTriggerCluster.h"
#include "AliTriggerConfiguration.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTrackPointArray.h"
#include "AliMCEventHandler.h"
#include "AliGenEventHeader.h"
#include "AliMCEvent.h"
#include "AliAODHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"   
#include "AliVertexerTracks.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliPID.h"
#include "AliITSgeomTGeo.h"
#include "AliAnalysisTaskSEImpParResSparse.h"

ClassImp(AliAnalysisTaskSEImpParResSparse)

//________________________________________________________________________
AliAnalysisTaskSEImpParResSparse::AliAnalysisTaskSEImpParResSparse():
AliAnalysisTaskSE(),
  fIsAOD(kFALSE),
  fReadMC(kFALSE),
  fSelectedPdg(-1),
  fUseDiamond(kFALSE),
  fSkipTrack(kTRUE),
  fMinMult(0),
  fMaxMult(1000000),
  fCheckSDDIsIn(0),
  fUseTriggerSelection(kFALSE),
  fTriggerClass(""),
  fTriggerMask(AliVEvent::kAny),
  fTrigConfig(0),
  fOCDBPath(""),
  fESDtrackCuts(0),
  fNentries(0),
  fMultiplicity(0),
  fUseCutGeoNcrNcl(kFALSE),
  fDeadZoneWidth(3.),
  fCutGeoNcrNclLength(130.),
  fCutGeoNcrNclGeom1Pt(1.5),
  fCutGeoNcrNclFractionNcr(0.85),
  fCutGeoNcrNclFractionNcl(0.7),
  fTrackType(0),
  fFillSparseForExpert(kFALSE),
  fImpParrphiSparsePtBchargePhi(0),
  fImpParrphiSparsePtzVtxEtaPhi(0),
  fImpParrphiSparsePtEtaPhi(0),
  fImpParPullrphiSparsePtEtaPhi(0),
  fImpParzSparsePtBchargePhi(0),
  fImpParzSparsePtzVtxEtaPhi(0),
  fImpParzSparsePtEtaPhi(0),
  fImpParPullzSparsePtEtaPhi(0),
  fPtDistrib(0),
  fhPtWeights(0x0),
  fUseptWeights(0),
  fScalingFactPtWeight(1.0),
  fOutput(0)
{
  //
  // Default constructor
  //
}

//________________________________________________________________________
AliAnalysisTaskSEImpParResSparse::AliAnalysisTaskSEImpParResSparse(const char *name):
  AliAnalysisTaskSE(name),
  fIsAOD(kFALSE),
  fReadMC(kFALSE),
  fSelectedPdg(-1),
  fUseDiamond(kFALSE),
  fSkipTrack(kTRUE),
  fMinMult(0),
  fMaxMult(1000000),
  fCheckSDDIsIn(0),
  fUseTriggerSelection(kFALSE),
  fTriggerClass(""),
  fTriggerMask(AliVEvent::kAny),
  fTrigConfig(0),
  fOCDBPath(""),
  fESDtrackCuts(0),
  fNentries(0),
  fMultiplicity(0),
  fUseCutGeoNcrNcl(kFALSE),
  fDeadZoneWidth(3.),
  fCutGeoNcrNclLength(130.),
  fCutGeoNcrNclGeom1Pt(1.5),
  fCutGeoNcrNclFractionNcr(0.85),
  fCutGeoNcrNclFractionNcl(0.7),
  fTrackType(0),
  fFillSparseForExpert(kFALSE),
  fImpParrphiSparsePtBchargePhi(0),
  fImpParrphiSparsePtzVtxEtaPhi(0),
  fImpParrphiSparsePtEtaPhi(0),
  fImpParPullrphiSparsePtEtaPhi(0),
  fImpParzSparsePtBchargePhi(0),
  fImpParzSparsePtzVtxEtaPhi(0),
  fImpParzSparsePtEtaPhi(0),
  fImpParPullzSparsePtEtaPhi(0),
  fPtDistrib(0),
  fhPtWeights(0x0),
  fUseptWeights(0),
  fScalingFactPtWeight(1.0),
  fOutput(0)
{
  //
  // Default constructor
  //

  DefineOutput(1, TList::Class());  //My private output
  
}

//________________________________________________________________________
AliAnalysisTaskSEImpParResSparse::~AliAnalysisTaskSEImpParResSparse()
{
  //
  // default distructor  
  // 
  //
  if (fESDtrackCuts)                                       { delete fESDtrackCuts;  fESDtrackCuts = 0;  }
  if (fNentries)                                           { delete fNentries;     fNentries    =0x0; }
  if (fMultiplicity)                                       { delete fMultiplicity;     fMultiplicity    =0x0; }

  delete fImpParrphiSparsePtBchargePhi;
  delete fImpParrphiSparsePtzVtxEtaPhi;
  delete fImpParrphiSparsePtEtaPhi;
  delete fImpParPullrphiSparsePtEtaPhi;
  delete fImpParzSparsePtBchargePhi;
  delete fImpParzSparsePtzVtxEtaPhi;
  delete fImpParzSparsePtEtaPhi;
  delete fImpParPullzSparsePtEtaPhi;
  delete fPtDistrib;
  delete fhPtWeights;

} 
//______________________________________________________________________________________________________
void AliAnalysisTaskSEImpParResSparse::UserCreateOutputObjects()
{
  // 
  // Create the output container
  //
  
  if(fDebug>1) printf("AnalysisTaskSEImpParRes::UserCreateOutputObjects() \n");
  
  // Several histograms are more conveniently managed in a TList
  if (!fOutput) {
    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("SparseList");
  }

  //THnSparses for experts
  if(fFillSparseForExpert){
    Int_t nbinsImpParSparse1[5] =       {3000, 50, 4, 2, 2};
    Double_t limitLowImpParSparse1[5] = {-1500., 0.1, 0., 0., 0.};
    Double_t limitUpImpParSparse1[5] =  {1500., 25., 4., 2., 2.};
    TString axTitle1[5]={"imp. par. (#mum)",
			 "#it{p}_{T} (GeV/c)",
			 "#phi",
			 "mag. field",
			 "charge"};
    fImpParrphiSparsePtBchargePhi=new THnSparseF("fImpParrphiSparsePtBchargePhi","fImpParrphiSparsePtBchargePhi",5,nbinsImpParSparse1,limitLowImpParSparse1,limitUpImpParSparse1);
    for(Int_t iax=0; iax<5; iax++) fImpParrphiSparsePtBchargePhi->GetAxis(iax)->SetTitle(axTitle1[iax].Data());
    BinLogAxis(fImpParrphiSparsePtBchargePhi, 1);
    fOutput->Add(fImpParrphiSparsePtBchargePhi);


    Int_t nbinsImpParSparse2[5] =       {3000, 50, 4, 2, 10};
    Double_t limitLowImpParSparse2[5] = {-1500., 0.1, 0., 0., -10.};
    Double_t limitUpImpParSparse2[5] =  {1500., 25., 4., 2., +10.};
    TString axTitle2[7]={"imp. par. (#mum)",
			 "#it{p}_{T} (GeV/c)",
			 "#phi",
			 "#eta",
			 "{z}_{vtx} (cm)"};
    fImpParrphiSparsePtzVtxEtaPhi=new THnSparseF("fImpParrphiSparsePtzVtxEtaPhi","fImpParrphiSparsePtzVtxEtaPhi",5,nbinsImpParSparse2,limitLowImpParSparse2,limitUpImpParSparse2);
    for(Int_t iax=0; iax<5; iax++) fImpParrphiSparsePtzVtxEtaPhi->GetAxis(iax)->SetTitle(axTitle2[iax].Data());
    BinLogAxis(fImpParrphiSparsePtzVtxEtaPhi, 1);
    fOutput->Add(fImpParrphiSparsePtzVtxEtaPhi);


    Int_t nbinsImpParSparse3[5] =       {1000, 50, 4, 2, 2};
    Double_t limitLowImpParSparse3[5] = {-1500., 0.1, 0., 0., 0.};
    Double_t limitUpImpParSparse3[5] =  {1500., 25., 4., 2., 2.};
    TString axTitle3[5]={"imp. par. (#mum)",
			 "#it{p}_{T} (GeV/c)",
			 "#phi",
			 "mag. field",
			 "charge"};
    fImpParzSparsePtBchargePhi=new THnSparseF("fImpParzSparsePtBchargePhi","fImpParzSparsePtBchargePhi",5,nbinsImpParSparse3,limitLowImpParSparse3,limitUpImpParSparse3);
    for(Int_t iax=0; iax<5; iax++) fImpParzSparsePtBchargePhi->GetAxis(iax)->SetTitle(axTitle3[iax].Data());
    BinLogAxis(fImpParzSparsePtBchargePhi, 1);
    fOutput->Add(fImpParzSparsePtBchargePhi);

    Int_t nbinsImpParSparse4[5] =       {1000, 50, 4, 2, 10};
    Double_t limitLowImpParSparse4[5] = {-1500., 0.1, 0., 0., -10.};
    Double_t limitUpImpParSparse4[5] =  {1500., 25., 4., 2., +10.};
    TString axTitle4[5]={"imp. par. (#mum)",
			 "#it{p}_{T} (GeV/c)",
			 "#phi",
			 "#eta",
			 "{z}_{vtx} (cm)"};
    fImpParzSparsePtzVtxEtaPhi=new THnSparseF("fImpParzSparsePtzVtxEtaPhi","fImpParzSparsePtzVtxEtaPhi",5,nbinsImpParSparse4,limitLowImpParSparse4,limitUpImpParSparse4);
    for(Int_t iax=0; iax<5; iax++) fImpParzSparsePtzVtxEtaPhi->GetAxis(iax)->SetTitle(axTitle4[iax].Data());
    BinLogAxis(fImpParzSparsePtzVtxEtaPhi, 1);
    fOutput->Add(fImpParzSparsePtzVtxEtaPhi);
  }

  //default THnSparses
  Int_t nbinsImpParSparse_rphi[4] =       {3000, 50, 4, 2};
  Double_t limitLowImpParSparse_rphi[4] = {-1500., 0.1, 0., 0.};
  Double_t limitUpImpParSparse_rphi[4] =  {1500., 25., 4., 2.};
  TString axTitle_rphi[4]={"rphi imp. par. (#mum)",
			   "#it{p}_{T} (GeV/c)",
			   "#phi",
			   "#eta"};
  fImpParrphiSparsePtEtaPhi=new THnSparseF("fImpParrphiSparsePtEtaPhi","fImpParrphiSparsePtEtaPhi",4,nbinsImpParSparse_rphi,limitLowImpParSparse_rphi,limitUpImpParSparse_rphi);
  for(Int_t iax=0; iax<4; iax++) fImpParrphiSparsePtEtaPhi->GetAxis(iax)->SetTitle(axTitle_rphi[iax].Data());
  BinLogAxis(fImpParrphiSparsePtEtaPhi, 1);
  fOutput->Add(fImpParrphiSparsePtEtaPhi);
  Int_t nbinsImpParSparse_z[4] =       {1000, 50, 4, 2};
  Double_t limitLowImpParSparse_z[4] = {-1500., 0.1, 0., 0.};
  Double_t limitUpImpParSparse_z[4] =  {1500., 25., 4., 2.};
  TString axTitle_z[4]={"z imp. par. (#mum)",
		       "#it{p}_{T} (GeV/c)",
		       "#phi",
		       "#eta"};
  fImpParzSparsePtEtaPhi=new THnSparseF("fImpParzSparsePtEtaPhi","fImpParzSparsePtEtaPhi",4,nbinsImpParSparse_z,limitLowImpParSparse_z,limitUpImpParSparse_z);
  for(Int_t iax=0; iax<4; iax++) fImpParzSparsePtEtaPhi->GetAxis(iax)->SetTitle(axTitle_z[iax].Data());
  BinLogAxis(fImpParzSparsePtEtaPhi, 1);
  fOutput->Add(fImpParzSparsePtEtaPhi);

  //pulls
  Int_t nbinsImpParSparse_pull[4] =       {400, 50, 4, 2};
  Double_t limitLowImpParSparse_pull[4] = {-10., 0.1, 0., 0.};
  Double_t limitUpImpParSparse_pull[4] =  {10., 25., 4., 2.};
  TString axTitle_rphi_pull[4]={"rphi pull",
				"#it{p}_{T} (GeV/c)",
				"#phi",
				"#eta"};
  fImpParPullrphiSparsePtEtaPhi=new THnSparseF("fImpParPullrphiSparsePtEtaPhi","fImpParPullrphiSparsePtEtaPhi",4,nbinsImpParSparse_pull,limitLowImpParSparse_pull,limitUpImpParSparse_pull);
  for(Int_t iax=0; iax<4; iax++) fImpParPullrphiSparsePtEtaPhi->GetAxis(iax)->SetTitle(axTitle_rphi_pull[iax].Data());
  BinLogAxis(fImpParPullrphiSparsePtEtaPhi, 1);
  fOutput->Add(fImpParPullrphiSparsePtEtaPhi);
  TString axTitle_z_pull[4]={"z pull",
			     "#it{p}_{T} (GeV/c)",
			     "#phi",
			     "#eta"};
  fImpParPullzSparsePtEtaPhi=new THnSparseF("fImpParPullzSparsePtEtaPhi","fImpParPullzSparsePtEtaPhi",4,nbinsImpParSparse_pull,limitLowImpParSparse_pull,limitUpImpParSparse_pull);
  for(Int_t iax=0; iax<4; iax++) fImpParPullzSparsePtEtaPhi->GetAxis(iax)->SetTitle(axTitle_z_pull[iax].Data());
  BinLogAxis(fImpParPullzSparsePtEtaPhi, 1);
  fOutput->Add(fImpParPullzSparsePtEtaPhi);
  

  fPtDistrib = new TH1F("fhpTdistr",";p_{T} (GeV/c)",500,0.1,25.);
  BinLogPtAxis(fPtDistrib);
  fOutput->Add(fPtDistrib);

  fhPtWeights = new TH1F("hptw",";p_{T} (GeV/c);weight",50,0.1,25.);
  BinLogPtAxis(fhPtWeights);
  ConfigurePtWeights();
  fOutput->Add(fhPtWeights);
  
  if(!fNentries) fNentries = new TH1F("hNentries", "number of entries", 26, 0., 40.);
  if(!fMultiplicity) fMultiplicity = new TH1F("fMultiplicity", "number of hits in SPD layer 1", 100, 0., 10000.);
 
  fOutput->Add(fNentries);
  fOutput->Add(fMultiplicity);

  PostData(1, fOutput);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEImpParResSparse::UserExec(Option_t */*option*/)
{
  if(fDebug>1) printf("AnalysisTaskSEImpParRes::UserExec() \n");
  //
  // Track selection and filling of d0 histograms
  //
  AliVEvent* event = dynamic_cast<AliVEvent*>(InputEvent());
  if (!event) {
    AliError("event not found. Nothing done!");
    return;
  }

  // only events in the requested multiplicity range
  TString firedTriggerClasses="";
  Int_t runNumber=0;
  if(fIsAOD){
    Int_t nclsITS = 0;
    runNumber=((AliAODEvent*)event)->GetRunNumber();
    nclsITS = ((AliVAODHeader*)((AliAODEvent*)event)->GetHeader())->GetNumberOfITSClusters(1);
    if(nclsITS<fMinMult || nclsITS>fMaxMult) return;
    firedTriggerClasses=((AliAODEvent*)event)->GetFiredTriggerClasses();
    if(!firedTriggerClasses.Contains(fTriggerClass.Data())) {/*Printf("%s",fTriggerClass.Data());*/ return;}
    if(fUseTriggerSelection){
      Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
      if(!isSelected) {/*Printf("EVENT NOT SELECTED");*/ return;}
    }
  }
  else{
    runNumber=((AliESDEvent*)event)->GetRunNumber();
    if(!IsSelectedCentrality(((AliESDEvent*)event))) return;
    firedTriggerClasses=((AliESDEvent*)event)->GetFiredTriggerClasses();
    if(!firedTriggerClasses.Contains(fTriggerClass.Data())) return;
  }

  Bool_t sddIsIn=kTRUE;
  if(fCheckSDDIsIn) {

    if(!fTrigConfig) {    
      AliCDBManager* man = AliCDBManager::Instance();
      if(fOCDBPath.Contains("OCDB")) { // when running in the QAtrain this is not called (OCBD is already set)
        man->SetDefaultStorage(fOCDBPath.Data());
        man->SetRun(runNumber);
      }
      if(!man) {      
        AliFatal("CDB not set but needed by AliAnalysisTaskITSTrackingCheck");
        return;    
      }  
      AliCDBEntry* eT=(AliCDBEntry*)man->Get("GRP/CTP/Config");    
      if(eT) {      
        fTrigConfig=(AliTriggerConfiguration*)eT->GetObject();    
      }    
      if(!eT || !fTrigConfig) {      
        AliError("Cannot retrieve CDB entry for GRP/CTP/Config");      
        return;    
      }
    }

    if(fIsAOD){
      const TObjArray& classesArray=fTrigConfig->GetClasses();
      ULong64_t trigMask=((AliAODEvent*)event)->GetTriggerMask();
      Int_t nclasses = classesArray.GetEntriesFast();
      for(Int_t iclass=0; iclass < nclasses; iclass++ ) 
	{
	  AliTriggerClass* trclass = (AliTriggerClass*)classesArray.At(iclass);
	  ULong64_t classMask=trclass->GetMask();
	  if(trigMask & classMask)
	    {
	      TString detList=trclass->GetCluster()->GetDetectorsInCluster();
	      if(detList.Contains("ITSSDD")) sddIsIn = kTRUE;
	      else sddIsIn = kFALSE;
	    }
	}
      //sddIsIn = kFALSE;
    }
    else sddIsIn=((AliESDEvent*)event)->IsDetectorInTriggerCluster("ITSSDD",fTrigConfig);
    if(fCheckSDDIsIn==1 && !sddIsIn) return;
    if(fCheckSDDIsIn==-1 && sddIsIn) return;
  }

  fNentries->Fill(1);
  if(fIsAOD)fMultiplicity->Fill(((AliVAODHeader*)((AliAODEvent*)event)->GetHeader())->GetNumberOfITSClusters(1));
  else{
    const AliMultiplicity *alimult =((AliESDEvent*)event)->GetMultiplicity();
    Int_t nclsSPDouter=0;
    if(alimult) {
      nclsSPDouter = alimult->GetNumberOfITSClusters(1);
    }
    fMultiplicity->Fill(nclsSPDouter);
  }
  
  Int_t nTrks = event->GetNumberOfTracks();
  Bool_t highMult=(nTrks>500 ? kTRUE : kFALSE);
  //Printf("%d",nTrks);

  Double_t vtxTrue[3];
  AliMCEvent* mcEvent=0x0;
  TClonesArray *mcArray=0;
  AliESDVertex *vtxESDTrue=0;
  AliVVertex *vtxVSkip=0;
  AliVVertex *vtxVRec=0;  
  AliVVertex* primaryVtx=0;


  // event primary vertex
  AliVertexerTracks vertexer0(event->GetMagneticField());
  vertexer0.SetITSMode();
  vertexer0.SetMinClusters(3);
  //if (fTrackType==1) vertexer0.SetITSpureSA(kTRUE);
  if(highMult) vertexer0.SetITSMode(0.1,0.1,0.5,5,1,3.,100.,1000.,3.,30.,1,1);
  if(fUseDiamond){
    // diamond constraint
    Float_t diamondcovxy[3];
    event->GetDiamondCovXY(diamondcovxy);
    Double_t pos[3]={event->GetDiamondX(),event->GetDiamondY(),0.};
    Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.};
    AliESDVertex diamond(pos,cov,1.,1);
    vertexer0.SetVtxStart(&diamond);
  }
  vtxVRec=(AliVVertex*)vertexer0.FindPrimaryVertex(event);
  if(!vtxVRec) {/*Printf("VERTEX REC NOT FOUND");*/ return;}
  if(vtxVRec->GetNContributors()<1){
    delete vtxVRec; vtxVRec=NULL;
    {/*Printf("VERTEX REC NUMBER OF CONTRIBUTORS < 1");*/return;}
  }
  
  if (fReadMC) {
    if (fIsAOD){
      mcArray = dynamic_cast<TClonesArray*>(((AliAODEvent*)event)->FindListObject(AliAODMCParticle::StdBranchName()));
      if(!mcArray){
	AliError("Clould not find Monte-Carlo in AOD");
	return;
      }
      AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(((AliAODEvent*)(event))->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
      if (!mcHeader) {
	AliError("Could not find MC Header in AOD");
	return;
      }
      
      Double_t mcVertex[3]={9999.,9999.,9999.};
      mcHeader->GetVertex(mcVertex);
      vtxTrue[0]=mcVertex[0];vtxTrue[1]=mcVertex[1];vtxTrue[2]=mcVertex[2];
      Double_t sigmaTrue[3]={0., 0., 0.,};
      vtxESDTrue = new AliESDVertex(vtxTrue,sigmaTrue);
    }//end if isAOD
    else{
      AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      if (!eventHandler) {
	Printf("ERROR: Could not retrieve MC event handler");
	return;
      }
      
      mcEvent = eventHandler->MCEvent();
      if (!mcEvent) {
	Printf("ERROR: Could not retrieve MC event");
	return;
      }
      
      //load MC header for ESD;//see $ALICE_PHYSICS/PWGPP/global/AliAnalysisTaskSEVertexESD.cxx
      AliHeader *mcHeader = eventHandler->MCEvent()->Header();
      if (!mcHeader) {
	AliDebug(AliLog::kError, "Header not available");
	return;
      }

      AliGenEventHeader* genHeader = mcHeader->GenEventHeader();  
      TArrayF mcVertex(3);
      mcVertex[0]=9999.; mcVertex[1]=9999.; mcVertex[2]=9999.;
      genHeader->PrimaryVertex(mcVertex);
      vtxTrue[0]=mcVertex[0];vtxTrue[1]=mcVertex[1];vtxTrue[2]=mcVertex[2];
      Double_t sigmaTrue[3]={0., 0., 0.,};
      //mcHeader->GetVertex(vtxTrue);//note the vtxTrue is void here,so must need the next line.
      //AliESDVertex *vtxESDTrue = new AliESDVertex(vtxTrue,sigmaTrue);
      vtxESDTrue = new AliESDVertex(vtxTrue,sigmaTrue);
      
    }//end else (!isAOD)
  }
  
  Double_t beampiperadius=3.;
  AliVTrack *vtrack = 0;
  Int_t pdgCode=0;
  Int_t trkLabel;
  TParticle  *part =0;
  AliAODMCParticle *AODpart=0;
  Int_t npointsITS=0,npointsSPD=0;
  Int_t skipped[2];
  Double_t dzRec[2], covdzRec[3], dzRecSkip[2], covdzRecSkip[3], dzTrue[2], covdzTrue[3];
  Double_t dz[2], covdz[3];
  Double_t pt;
  Int_t bin;
  Int_t nClsTotTPC=0;
  Bool_t haskITSrefit=kFALSE;
  Bool_t haskTPCrefit=kFALSE;
  Int_t charge=0;
  Double_t phi=0.;
  Double_t theta=0.;
  Double_t eta=0.;
  Double_t pointrphi[4];
  Double_t pullrphi[4];
  Double_t pointrphi1[5];
  Double_t pointrphi2[5];
  Double_t pointz[4];
  Double_t pullz[4];
  Double_t pointz1[5];
  Double_t pointz2[5];

  //Printf("\nSTART LOOP OVER TRACKS\n");

  for (Int_t it=0; it<nTrks; it++){ //start loop over tracks
    vtrack = (AliVTrack*)event->GetTrack(it);
    if(!vtrack) {/*Printf("TRACK NOT FOUND");*/ continue;}

    //Printf("ETA %f",eta);
    
    npointsITS=0; npointsSPD=0;
    if(fIsAOD){
      haskITSrefit=(((AliAODTrack*)vtrack)->GetStatus()&AliESDtrack::kITSrefit);
      haskTPCrefit=(((AliAODTrack*)vtrack)->GetStatus()&AliESDtrack::kTPCrefit);
      nClsTotTPC=((AliAODTrack*)vtrack)->GetTPCNcls();
      //Printf("nClsTotTPC %d",nClsTotTPC);
      if(!haskITSrefit) continue;
      for(Int_t ilayer=0; ilayer<6; ilayer++){
	if (ilayer<2 && ((AliAODTrack*)vtrack)->HasPointOnITSLayer(ilayer)) npointsSPD++;
	if (((AliAODTrack*)vtrack)->HasPointOnITSLayer(ilayer)) npointsITS++;  
      }
      //Printf("npointsITS %d",npointsITS);
    }
    else {
      haskITSrefit=(((AliESDtrack*)vtrack)->GetStatus()&AliESDtrack::kITSrefit);
      haskTPCrefit=(((AliESDtrack*)vtrack)->GetStatus()&AliESDtrack::kTPCrefit);
      nClsTotTPC=((AliESDtrack*)vtrack)->GetTPCNcls();
      if(!haskITSrefit) continue;
      for (Int_t ilayer=0; ilayer<6; ilayer++){ 
	if (ilayer<2 && ((AliESDtrack*)vtrack)->HasPointOnITSLayer(ilayer)) npointsSPD++;
	if (((AliESDtrack*)vtrack)->HasPointOnITSLayer(ilayer)) npointsITS++;  
      }
    }

    pt = vtrack->Pt();
    Double_t weight=pt<fhPtWeights->GetBinLowEdge(fhPtWeights->GetNbinsX()+1) ? fhPtWeights->GetBinContent(fhPtWeights->FindBin(pt)) : 1.;
    if (pt > 1000.) continue;
    if( ((Double_t)pt*10000.)-((Long_t)(pt*10000.))>weight) continue;
    
    pullrphi[1]=pt;
    pointrphi[1]=pt;
    pointrphi1[1]=pt;
    pointrphi2[1]=pt;
    pullz[1]=pt;
    pointz[1]=pt;
    pointz1[1]=pt;
    pointz2[1]=pt;

    eta = vtrack->Eta();
    if(eta<-0.8 || eta>0.8) continue;
    if(eta<0.) {
      pullrphi[3]=0.;
      pullz[3]=0.;
      pointrphi[3]=0.;
      pointz[3]=0.;
      pointrphi2[3]=0.;
      pointz2[3]=0.;
    }
    else if(eta>=0.) {
      pullrphi[3]=1.;
      pullz[3]=1.;
      pointrphi[3]=1.;
      pointz[3]=1.;
      pointrphi2[3]=1.;
      pointz2[3]=1.;
    }

    charge=vtrack->Charge();
    if(charge<0.) {pointrphi1[4]=0.; pointz1[4]=0.;}
    else if(charge>0.) {pointrphi1[4]=1.; pointz1[4]=1.;}
    //Printf("charge %d",charge);

    phi=vtrack->Phi();
    Int_t phibin=PhiBin(phi);
    if(phibin<0) continue;
    pullrphi[2]=phibin;
    pullz[2]=phibin;
    pointrphi[2]=phibin;
    pointz[2]=phibin;
    pointrphi1[2]=phibin;
    pointz1[2]=phibin;
    pointrphi2[2]=phibin;
    pointz2[2]=phibin;
    //Printf("phibin %d",phibin);


    Float_t magField=event->GetMagneticField();
    if(magField<0.) {pointrphi1[3]=0.; pointz1[3]=0.;}
    else if(magField>0.) {pointrphi1[3]=1.; pointz1[3]=1.;}
    

    //MC
    if (fReadMC){
      trkLabel = vtrack->GetLabel();
      if(trkLabel<0) continue;
      if(fIsAOD && mcArray){
	AODpart = (AliAODMCParticle*)mcArray->At(trkLabel);
	if(!AODpart) printf("NOPART\n");
	pdgCode = TMath::Abs(AODpart->GetPdgCode());	
      }
      if(!fIsAOD && mcEvent) {
	part = ((AliMCParticle*)mcEvent->GetTrack(trkLabel))->Particle();
	pdgCode = TMath::Abs(part->GetPdgCode());
      }
      //pdgCode = TMath::Abs(part->GetPdgCode());
      //printf("pdgCode===%d\n", pdgCode);
      if(fSelectedPdg>0 && pdgCode!=fSelectedPdg) continue;
    }
    
      
    //Get specific primary vertex--Reconstructed primary vertex do not include the track considering.
    AliVertexerTracks vertexer(event->GetMagneticField());
    vertexer.SetITSMode();
    //if (fTrackType==1) vertexer.SetITSpureSA(kTRUE);
    vertexer.SetMinClusters(3);
    if(fUseDiamond){
      Float_t diamondcovxy[3];
      event->GetDiamondCovXY(diamondcovxy);
      Double_t pos[3]={event->GetDiamondX(),event->GetDiamondY(),0.};
      Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.};
      AliESDVertex diamond(pos,cov,1.,1);
      vertexer.SetVtxStart(&diamond);
    }
    skipped[0] = (Int_t)vtrack->GetID();
    vertexer.SetSkipTracks(1,skipped);      
    // create vertex with new!
    if(!highMult && fSkipTrack) {
      vtxVSkip = (AliVVertex*)vertexer.FindPrimaryVertex(event);
      if(!vtxVSkip) {/*Printf("VERTEX SKIP NOT FOUND");*/ continue;}
      if(vtxVSkip->GetNContributors()<1) {
	delete vtxVSkip; vtxVSkip=NULL;
	continue;
      }
    } // else {
      // vtxVSkip = new AliVVertex(); produce error!!!
      // } 
 

    // Select primary particle if MC event (for ESD event), Rprod < 1 micron
    if(fReadMC){
      if(fIsAOD){
	if((AODpart->Xv()-vtxTrue[0])*(AODpart->Xv()-vtxTrue[0])+
	   (AODpart->Yv()-vtxTrue[1])*(AODpart->Yv()-vtxTrue[1])
	   > 0.0001*0.0001) {
	  delete vtxVSkip; vtxVSkip=NULL;
	  continue;
	}
      }
      else{
	if((part->Vx()-vtxTrue[0])*(part->Vx()-vtxTrue[0])+
	   (part->Vy()-vtxTrue[1])*(part->Vy()-vtxTrue[1])
	   > 0.0001*0.0001) {
	  delete vtxVSkip; vtxVSkip=NULL;
	  continue;
	}
      }
    }
    
    
    // compute impact parameters
    // wrt event vertex
    vtrack->PropagateToDCA(vtxVRec, event->GetMagneticField(), beampiperadius, dzRec, covdzRec);
    dz[0]=dzRec[0]; 
    dz[1]=dzRec[1];
    covdz[0]=covdzRec[0];
    covdz[1]=covdzRec[1];
    covdz[2]=covdzRec[2];
    // wrt event vertex without this track
    if(!highMult && fSkipTrack) {
      vtrack->PropagateToDCA(vtxVSkip, event->GetMagneticField(), beampiperadius, dzRecSkip, covdzRecSkip);
      dz[0]=dzRecSkip[0]; 
      dz[1]=dzRecSkip[1];
      covdz[0]=covdzRecSkip[0];
      covdz[1]=covdzRecSkip[1];
      covdz[2]=covdzRecSkip[2];
    } else if(!fSkipTrack) {
      dz[0]=dzRec[0]; 
      dz[1]=dzRec[1];
      covdz[0]=covdzRec[0];
      covdz[1]=covdzRec[1];
      covdz[2]=covdzRec[2];
    } else {
      dz[0]=0; 
      dz[1]=0;
      covdz[0]=0;
      covdz[1]=0;
      covdz[2]=0;
    }
    //delete vtxVSkip; vtxVSkip=NULL; // not needed anymore
 
    if(fReadMC) {
      vtrack->PropagateToDCA(vtxESDTrue, event->GetMagneticField(), beampiperadius, dzTrue, covdzTrue);
      dz[0]=dzTrue[0]; 
      dz[1]=dzTrue[1];
      covdz[0]=covdzTrue[0];
      covdz[1]=covdzTrue[1];
      covdz[2]=covdzTrue[2];
    }
    if(covdz[0]<1.e-13 || covdz[2]<1.e-13) continue;
    pointrphi[0]=10000.*dz[0];
    pointrphi1[0]=10000.*dz[0];
    pointrphi2[0]=10000.*dz[0];
    pointz[0]=10000.*dz[1];
    pointz1[0]=10000.*dz[1];
    pointz2[0]=10000.*dz[1];
    pullrphi[0]=dz[0]/TMath::Sqrt(covdz[0]);
    pullz[0]=dz[1]/TMath::Sqrt(covdz[2]);
    //    Printf("point %f",pointrphi[0]);

    if(fReadMC) primaryVtx=vtxESDTrue;
    else if(fSkipTrack) primaryVtx=vtxVSkip;
    else primaryVtx=vtxVRec;

    Double_t zvtx=primaryVtx->GetZ();
    pointrphi2[4]=zvtx;
    pointz2[4]=zvtx;

    // 6 points in ITS
    if(fTrackType==0 && (haskTPCrefit && nClsTotTPC>=70 && (npointsITS==6 || (npointsITS==4 && !sddIsIn)))) {
      if(fFillSparseForExpert){
	fImpParrphiSparsePtBchargePhi->Fill(pointrphi1);
	fImpParrphiSparsePtzVtxEtaPhi->Fill(pointrphi2);
	fImpParzSparsePtBchargePhi->Fill(pointz1);
	fImpParzSparsePtzVtxEtaPhi->Fill(pointz2);
      }
      fPtDistrib->Fill(pt);
      fImpParrphiSparsePtEtaPhi->Fill(pointrphi);
      fImpParzSparsePtEtaPhi->Fill(pointz);
      fImpParPullrphiSparsePtEtaPhi->Fill(pullrphi);
      fImpParPullzSparsePtEtaPhi->Fill(pullz);
    }
    // ITS standalone
    else if(fTrackType==1 && nClsTotTPC==0 && haskITSrefit && npointsSPD>0 && npointsITS>=4) {
      if(fFillSparseForExpert){
	fPtDistrib->Fill(pt);
	fImpParrphiSparsePtBchargePhi->Fill(pointrphi1);
	fImpParrphiSparsePtzVtxEtaPhi->Fill(pointrphi2);
	fImpParzSparsePtBchargePhi->Fill(pointz1);
	fImpParzSparsePtzVtxEtaPhi->Fill(pointz2);
      }
      fPtDistrib->Fill(pt);
      fImpParrphiSparsePtEtaPhi->Fill(pointrphi);
      fImpParzSparsePtEtaPhi->Fill(pointz);
      fImpParPullrphiSparsePtEtaPhi->Fill(pullrphi);
      fImpParPullzSparsePtEtaPhi->Fill(pullz);
    }
    // ESD TRACK CUTS
    else if(fTrackType==2 && IsTrackSelected(vtrack,primaryVtx,fESDtrackCuts,event)){
       if(fFillSparseForExpert){
	fPtDistrib->Fill(pt);
	fImpParrphiSparsePtBchargePhi->Fill(pointrphi1);
	fImpParrphiSparsePtzVtxEtaPhi->Fill(pointrphi2);
	fImpParzSparsePtBchargePhi->Fill(pointz1);
	fImpParzSparsePtzVtxEtaPhi->Fill(pointz2);
      }
      fPtDistrib->Fill(pt);
      fImpParrphiSparsePtEtaPhi->Fill(pointrphi);
      fImpParzSparsePtEtaPhi->Fill(pointz);
      fImpParPullrphiSparsePtEtaPhi->Fill(pullrphi);
      fImpParPullzSparsePtEtaPhi->Fill(pullz);
    }  
    
    

  }//end loop over tracks

  
  delete vtxVSkip; vtxVSkip=NULL;
  delete vtxVRec;  vtxVRec=NULL;
  delete vtxESDTrue; vtxESDTrue=NULL;
  PostData(1, fOutput);
 

  return;
}
//________________________________________________________________________
void AliAnalysisTaskSEImpParResSparse::BinLogAxis(const THnSparseF *h, Int_t axisNumber) {
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis *axis = h->GetAxis(axisNumber);
  int bins = axis->GetNbins();
  
  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
     
  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
    
  for (int i = 1; i <= bins; i++) {
    newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;
    
}
//________________________________________________________________________
void AliAnalysisTaskSEImpParResSparse::BinLogPtAxis(TH1F *h) {
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis *axis = ((TH1F*)h)->GetXaxis();
  int bins = axis->GetNbins();
  
  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
     
  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
    
  for (int i = 1; i <= bins; i++) {
    newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;
    
} 
//________________________________________________________________________

Int_t AliAnalysisTaskSEImpParResSparse::PhiBin(Double_t phi) const { 
  Double_t pi=TMath::Pi();
  if(phi>2.*pi || phi<0.) return -1;
  if((phi<=(pi/4.)) || (phi>7.*(pi/4.))) return 0;
  if((phi>(pi/4.)) && (phi<=3.*(pi/4.))) return 1;
  if((phi>3.*(pi/4.)) && (phi<=5.*(pi/4.))) return 2;
  if((phi>(5.*pi/4.)) && (phi<=7.*(pi/4.))) return 3;
  return -1;
}
//___________________________________________________________________________
void AliAnalysisTaskSEImpParResSparse::Terminate(Option_t */*option*/) {
  //
  // Terminate analysis
  //

  if (fDebug>1) printf("AnalysisTaskSEImpParRes: Terminate() \n");

  return;
}
//__________________________________________________________________________
Int_t AliAnalysisTaskSEImpParResSparse::ClusterTypeOnITSLayer(AliESDtrack *track,
							    Int_t layer) const {
  //
  // Returns cluster type on ITS layer. Returns -1 if no cluster on this layer
  //
  Int_t ctype=-1;

  if(layer<0 || layer>5) return ctype;
  if(!track->HasPointOnITSLayer(layer)) return ctype;
  
  const AliTrackPointArray *array = track->GetTrackPointArray();
  if(!array) {
    //    printf("No tracks points avaialble: check ESDfriends\n");
    return ctype;
  }
  AliTrackPoint point;
  Int_t ipt,volId,modId,layerId;
  for(ipt=0; ipt<array->GetNPoints(); ipt++) {
    array->GetPoint(point,ipt);
    volId = point.GetVolumeID();
    if(volId<=0) continue;
    layerId = AliGeomManager::VolUIDToLayer(volId,modId);
    if(layerId==layer+1 && !point.IsExtra()) {
      ctype = point.GetClusterType();
      break;
    }
  }
  return ctype;
}
//---------------------------------------------------------------------------
Bool_t AliAnalysisTaskSEImpParResSparse::IsSelectedCentrality(AliESDEvent *esd) const
{
  //
  // check if events is in the required multiplicity range
  //

  const AliMultiplicity *alimult = esd->GetMultiplicity();
  Int_t ntrklets=1;
  Int_t nclsSPDouter=0;
  if(alimult) {
    ntrklets = alimult->GetNumberOfTracklets();
    nclsSPDouter = alimult->GetNumberOfITSClusters(1);
  }

  if(nclsSPDouter<fMinMult || nclsSPDouter>fMaxMult) return kFALSE;


  return kTRUE;
}

//----------------------------------------------------------------------------------
Bool_t AliAnalysisTaskSEImpParResSparse::IsTrackSelected(AliVTrack *track, AliVVertex *primary, AliESDtrackCuts *cuts, const AliVEvent* aod) const{

  if(!cuts) return kTRUE;
  Bool_t retval = kTRUE;
  if(fIsAOD) {
    AliESDtrack esdTrack(track);
    esdTrack.SetTPCClusterMap(((AliAODTrack*)track)->GetTPCClusterMap());
    esdTrack.SetTPCSharedMap(((AliAODTrack*)track)->GetTPCSharedMap());
    esdTrack.SetTPCPointsF(((AliAODTrack*)track)->GetTPCNclsF());
    esdTrack.RelateToVertex((AliESDVertex*)primary,0.,3.);
    if(!cuts->IsSelected(&esdTrack)) retval = kFALSE;
    //Printf("NOT Checking geo cut for AOD: %d",retval);

    if(fUseCutGeoNcrNcl){
      Float_t nCrossedRowsTPC = esdTrack.GetTPCCrossedRows();
      Float_t lengthInActiveZoneTPC=esdTrack.GetLengthInActiveZone(0,fDeadZoneWidth,220.,aod->GetMagneticField());
      Double_t cutGeoNcrNclLength=fCutGeoNcrNclLength-TMath::Power(TMath::Abs(esdTrack.GetSigned1Pt()),fCutGeoNcrNclGeom1Pt);
      Bool_t isOK=kTRUE;
      if (lengthInActiveZoneTPC<cutGeoNcrNclLength) isOK=kFALSE;
      if (nCrossedRowsTPC<fCutGeoNcrNclFractionNcr*cutGeoNcrNclLength) isOK=kFALSE;
      if (esdTrack.GetTPCncls()<fCutGeoNcrNclFractionNcl*cutGeoNcrNclLength) isOK=kFALSE;
      if(!(isOK && retval)) retval = kFALSE;
      //Printf("Checking geo cut for AOD: %f %f isOK=%d retval=%d",lengthInActiveZoneTPC,cutGeoNcrNclLength,isOK,retval);
    }
    
  }
  else {
    AliESDtrack *esdTrack = (AliESDtrack*)track;
    if(!cuts->IsSelected(esdTrack)) retval = kFALSE;
  }
  return retval;
}
//----------------------------------------------------------------------------------
void AliAnalysisTaskSEImpParResSparse::ConfigurePtWeights(){

  if(fUseptWeights==0){ // no weights
    for(Int_t i=0; i<fhPtWeights->GetNbinsX(); i++){
      fhPtWeights->SetBinContent(i+1,1.);
    }
  }
  else if(fUseptWeights==1){ //pp weights
    for(Int_t i=0; i<fhPtWeights->GetXaxis()->FindBin(7.); i++){
      fhPtWeights->SetBinContent(i+1,0.01);
    }
    for(Int_t i=fhPtWeights->GetXaxis()->FindBin(7.); i<fhPtWeights->GetNbinsX(); i++){
      fhPtWeights->SetBinContent(i+1,1.);
    }
  }
  else if(fUseptWeights==2){ //pPb weights
    for(Int_t i=0; i<fhPtWeights->GetXaxis()->FindBin(7.); i++){
      fhPtWeights->SetBinContent(i+1,0.001);
    }
    for(Int_t i=fhPtWeights->GetXaxis()->FindBin(7.); i<fhPtWeights->GetNbinsX(); i++){
      fhPtWeights->SetBinContent(i+1,0.5);
    }
  }
  else if(fUseptWeights==3){ //PbPb weights
    for(Int_t i=0; i<fhPtWeights->GetXaxis()->FindBin(7.); i++){
      fhPtWeights->SetBinContent(i+1,0.0001);
    }
    for(Int_t i=fhPtWeights->GetXaxis()->FindBin(7.); i<fhPtWeights->GetNbinsX(); i++){
      fhPtWeights->SetBinContent(i+1,0.05);
    }
  }

  fhPtWeights->Scale(1./fScalingFactPtWeight);
  
  
}
