
/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//*************************************************************************
// Class AliAnalysisTaskSEDplusCorrelations
// AliAnalysisTaskSE for the D+ candidates Invariant Mass Histogram and 
//comparison of heavy-flavour decay candidates
// to MC truth (kinematics stored in the AOD)
// Authors: Sadhana Dash (correlation)
/////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TList.h>
#include <TString.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODPidHF.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliVertexingHFUtils.h"
#include "AliAnalysisTaskSEDplusCorrelations.h"
#include "AliNormalizationCounter.h"
#include "AliVParticle.h"
#include "AliHFAssociatedTrackCuts.h"
#include "AliReducedParticle.h"
#include "AliHFCorrelator.h"


using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSEDplusCorrelations)


//________________________________________________________________________
AliAnalysisTaskSEDplusCorrelations::AliAnalysisTaskSEDplusCorrelations():
AliAnalysisTaskSE(),
  fOutput(0),
  fCorrelator(0),
  fSelect(0),
  fDisplacement(0),
  fHistNEvents(0),
  fTrig(kTRUE),
  fEventMix(0),
  fUpmasslimit(1.965),
  fLowmasslimit(1.765),
  fNPtBins(0),
  fBinWidth(0.002),
  fListCuts(0),
//fListCutsAsso(0), 
  fRDCutsAnalysis(0),
  fCuts(0),
  fCounter(0),
  fReadMC(kFALSE),
  fUseBit(kTRUE),
  fMixing(kFALSE),
  fSystem(kFALSE),
  fReco(kTRUE)
{
  // Default constructor
  
  for(Int_t i=0;i<3*kMaxPtBins;i++){
   
    if(fMassHistK0S[i])fMassHistK0S[i]=0;
    if(fLeadPt[i])fLeadPt[i]=0;
    if(fMassHist[i])fMassHist[i]=0;
    if(fMassHistTC[i])fMassHistTC[i]=0;
    if(fMassHistOrigC[i])fMassHistOrigC[i]=0;
    if(fMassHistOrigB[i])fMassHistOrigB[i]=0;
    if(fMassHistMC[i])fMassHistMC[i]=0;
  }
  
  for(Int_t i=0;i<kMaxPtBins+1;i++){
    fArrayBinLimits[i]=0;
  }

}

//________________________________________________________________________
AliAnalysisTaskSEDplusCorrelations::AliAnalysisTaskSEDplusCorrelations(const char *name,AliRDHFCutsDplustoKpipi *dpluscutsana, AliHFAssociatedTrackCuts *asstrkcuts):
  AliAnalysisTaskSE(name),
  fOutput(0),
  fCorrelator(0),
  fSelect(0),
  fDisplacement(0),
  fHistNEvents(0),
  fTrig(kTRUE),
  fEventMix(0),
  fUpmasslimit(1.965),
  fLowmasslimit(1.765),
  fNPtBins(0),
  fBinWidth(0.002),
  fListCuts(0),
  //fListCutsAsso(0),
  fRDCutsAnalysis(dpluscutsana),
  fCuts(asstrkcuts),
  fCounter(0),
  fReadMC(kFALSE),
  fUseBit(kTRUE),
  fMixing(kFALSE),
  fSystem(kFALSE),
  fReco(kTRUE)
{
  // 
  // Standrd constructor
  //
  fNPtBins=fRDCutsAnalysis->GetNPtBins();
    
    
  for(Int_t i=0;i<3*kMaxPtBins;i++){
   
    if(fMassHistK0S[i])fMassHistK0S[i]=0;
    if(fLeadPt[i])fLeadPt[i]=0;
    if(fMassHist[i])fMassHist[i]=0;
    if(fMassHistTC[i])fMassHistTC[i]=0;
    if(fMassHistOrigC[i])fMassHistOrigC[i]=0;
    if(fMassHistOrigB[i])fMassHistOrigB[i]=0;
    if(fMassHistMC[i])fMassHistMC[i]=0;
   
    
    
  }
  
  for(Int_t i=0;i<kMaxPtBins+1;i++){
    fArrayBinLimits[i]=0;
  }
  
  
  // Default constructor
  // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
  // Output slot #2 writes cut to private output
  DefineOutput(2,TList::Class());
  // Output slot #3 writes cut to private output
  DefineOutput(3,AliNormalizationCounter::Class());
  //DefineOutput(4,AliHFAssociatedTrackCuts::Class());
  
 
}

//________________________________________________________________________
AliAnalysisTaskSEDplusCorrelations::~AliAnalysisTaskSEDplusCorrelations()
{
  //
  // Destructor
  //
  if(fOutput) {delete fOutput; fOutput = 0;}
  if(fListCuts) {delete fListCuts; fListCuts = 0;}
  if(fRDCutsAnalysis) {delete fRDCutsAnalysis; fRDCutsAnalysis = 0;}
  if(fCuts) {delete fCuts; fCuts = 0;}
  if(fCounter) {delete fCounter; fCounter = 0;}
  if(fCorrelator) {delete fCorrelator; fCorrelator=0;}



}  
//_________________________________________________________________
void  AliAnalysisTaskSEDplusCorrelations::SetMassLimits(Float_t range){
  // set invariant mass limits
  Float_t bw=GetBinWidth();
  fUpmasslimit = 1.865+range;
  fLowmasslimit = 1.865-range;
  SetBinWidth(bw);
}
//_________________________________________________________________
void  AliAnalysisTaskSEDplusCorrelations::SetMassLimits(Float_t lowlimit, Float_t uplimit){
  // set invariant mass limits
  if(uplimit>lowlimit)
    {
      Float_t bw=GetBinWidth();
      fUpmasslimit = lowlimit;
      fLowmasslimit = uplimit;
      SetBinWidth(bw);
    }
}
//________________________________________________________________
void AliAnalysisTaskSEDplusCorrelations::SetBinWidth(Float_t w){
  // Define width of mass bins
  Float_t width=w;
  Int_t nbins=(Int_t)((fUpmasslimit-fLowmasslimit)/width+0.5);
  Int_t missingbins=4-nbins%4;
  nbins=nbins+missingbins;
  width=(fUpmasslimit-fLowmasslimit)/nbins;
  if(missingbins!=0){
    printf("AliAnalysisTaskSEDplusCorrelations::SetBinWidth: W-bin width of %f will produce histograms not rebinnable by 4. New width set to %f\n",w,width);
  }
  else{
    if(fDebug>1) printf("AliAnalysisTaskSEDplusCorrelations::SetBinWidth: width set to %f\n",width);
  }
  fBinWidth=width;
}
//_________________________________________________________________
Int_t AliAnalysisTaskSEDplusCorrelations::GetNBinsHistos(){
  // Compute number of mass bins
  return (Int_t)((fUpmasslimit-fLowmasslimit)/fBinWidth+0.5);
}


//__________________________________________
void AliAnalysisTaskSEDplusCorrelations::Init(){
  //
  // Initialization
  //
  if(fDebug > 1) printf("AnalysisTaskSEDplusCorrelations::Init() \n");
  

  fListCuts=new TList();

  // fListCutsAsso=new TList();
  
  AliRDHFCutsDplustoKpipi *analysis = new AliRDHFCutsDplustoKpipi(*fRDCutsAnalysis);
  analysis->SetName("AnalysisCuts");

  // AliHFAssociatedTrackCuts *trkcuts = new AliHFAssociatedTrackCuts(*fCuts);
  //trkcuts->SetName("Assotrkcuts");
  
  fListCuts->Add(analysis);
  //fListCutsAsso->Add(trkcuts);

   

  PostData(2,fListCuts);
  //≈ß PostData(4,fListCutsAsso);
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDplusCorrelations::UserCreateOutputObjects()
{
  // Create the output container
  //

  if(fDebug > 1) printf("AnalysisTaskSEDplusCorrelations::UserCreateOutputObjects() \n");
  // Correlator creation and definition

  Double_t Pi = TMath::Pi();

  // Set up the correlator


  fCorrelator = new AliHFCorrelator("Correlator",fCuts,fSystem); // fCuts is the hadron cut object, fSystem to switch between pp or PbPb

  fCorrelator->SetDeltaPhiInterval((-0.5-1./32)*Pi,(1.5-1./32)*Pi); // set correct phi interval

  fCorrelator->SetEventMixing(fMixing); //set kFALSE/kTRUE for mixing Off/On

  fCorrelator->SetAssociatedParticleType(fSelect); // set 1/2/3 for hadron/kaons/kzeros

  fCorrelator->SetApplyDisplacementCut(fDisplacement); //set kFALSE/kTRUE for using the displacement cut
  
  fCorrelator->SetUseMC(fReadMC);
  fCorrelator->SetPIDmode(2);
  fCorrelator->SetUseReco(fReco); // to analyse reco/kine
 
  // For Event Mixing
 
  Bool_t pooldef = fCorrelator->DefineEventPool();
  
  if(!pooldef) AliInfo("Warning:: Event pool not defined properly");


  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  TString hisname;
  Int_t index=0;
  Int_t nbins=GetNBinsHistos();


  
  for(Int_t i=0;i<fNPtBins;i++){

    index=GetHistoIndex(i);
    
    hisname.Form("hMassPtK0S%d",i);
    fMassHistK0S[index]=new TH1F(hisname.Data(),hisname.Data(),500,0.3,0.8);
    fMassHistK0S[index]->Sumw2();

    hisname.Form("hLeadPt%d",i);
    fLeadPt[index]=new TH1F(hisname.Data(),hisname.Data(),500,0.0,50);
    fLeadPt[index]->Sumw2();

          
    hisname.Form("hMassPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();

      
    hisname.Form("hMassPt%dTC",i);
    fMassHistTC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistTC[index]->Sumw2();

       
    index=GetSignalHistoIndex(i); 
    
    hisname.Form("hSigPtK0S%d",i);
    fMassHistK0S[index]=new TH1F(hisname.Data(),hisname.Data(),500,0.3,0.8);
    fMassHistK0S[index]->Sumw2();

    hisname.Form("hSigLeadPt%d",i);
    fLeadPt[index]=new TH1F(hisname.Data(),hisname.Data(),500,0.0,50);
    fLeadPt[index]->Sumw2();
    
    hisname.Form("hSigPt%d",i);
    fMassHist[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHist[index]->Sumw2();

    hisname.Form("hSigOrigCPt%d",i);
    fMassHistOrigC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistOrigC[index]->Sumw2();

    hisname.Form("hSigOrigBPt%d",i);
    fMassHistOrigB[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistOrigB[index]->Sumw2();

      
    hisname.Form("hSigMCPt%d",i);
    fMassHistMC[index]=new TH1F(hisname.Data(),hisname.Data(),nbins,fLowmasslimit,fUpmasslimit);
    fMassHistMC[index]->Sumw2();
   
     }
    

  for(Int_t i=0; i<3*fNPtBins; i++){
   
    fOutput->Add(fMassHistK0S[i]);
    fOutput->Add(fLeadPt[i]);
    fOutput->Add(fMassHist[i]);
    fOutput->Add(fMassHistTC[i]);
    fOutput->Add(fMassHistOrigC[i]);
    fOutput->Add(fMassHistOrigB[i]);
    fOutput->Add(fMassHistMC[i]);
   
  }
  
  
  TH1D * EventTypeMC = new TH1D("EventTypeMC","EventTypeMC",100,-0.5,99.5);
  if(fReadMC) fOutput->Add(EventTypeMC);

    	
  //Event Mixing histos (for control plots)

  Int_t NumberEvents = fCuts->GetMaxNEventsInPool();
  Int_t NofTracks = fCuts->GetMinNTracksInPool();
  Int_t NofCentBins = fCuts->GetNCentPoolBins();
  Double_t * CentBins = fCuts->GetCentPoolBins();
  Int_t NumberZVtxBins = fCuts->GetNZvtxPoolBins();
  Double_t *ZVtxBins = fCuts->GetZvtxPoolBins();
  
  
  
  fEventMix = new TH2D("EventMixingCheck","EventMixingCheck",100,0,600,100,-15,15);
  
  fEventMix->GetXaxis()->SetTitle("Multiplicity ");
  fEventMix->GetYaxis()->SetTitle("Z vertex [cm]");
  
  if(fMixing)fOutput->Add(fEventMix);
  
  TH3D * EventsInPool = new TH3D("EventsInPool","Number of events in  pool",NofCentBins,0,600,NumberZVtxBins,-15,15,1000000,0,1000000);
  
  EventsInPool->GetXaxis()->SetTitle("Multiplicity ");
  EventsInPool->GetYaxis()->SetTitle("Z vertex [cm]");
  EventsInPool->GetZaxis()->SetTitle("Number of events in pool");
  if(fMixing) fOutput->Add(EventsInPool);
  
  Int_t MaxNofTracks = (NumberEvents+1)*NofTracks;
  
  
  TH3D * NofTracksInPool = new TH3D("NofTracksInPool","Number of tracks in  pool",NofCentBins,0,500,NumberZVtxBins,-15,15,MaxNofTracks,0,MaxNofTracks);
  NofTracksInPool->GetXaxis()->SetTitle("Multiplicity ");
  NofTracksInPool->GetYaxis()->SetTitle("Z vertex [cm]");
  NofTracksInPool->GetZaxis()->SetTitle("Number of tracks");
  
  if(fMixing) fOutput->Add(NofTracksInPool);
  
  TH2D * NofPoolBinCalls = new TH2D("NofPoolBinCalls","Calls per pool bin",NofCentBins,CentBins,NumberZVtxBins,ZVtxBins);
  NofPoolBinCalls->GetXaxis()->SetTitle("Multiplicity ");
  NofPoolBinCalls->GetYaxis()->SetTitle("Z vertex [cm]");
  if(fMixing) fOutput->Add(NofPoolBinCalls);
  
    	
  
  fHistNEvents = new TH1F("fHistNEvents", "number of events ",10,-0.5,10.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"nEvents accepted");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"nEvent with good vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"Total number of candidates");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"Number without bitmask"); 
  fHistNEvents->GetXaxis()->SetBinLabel(6,"Number after Physics Selection");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"Total number in Fiducial Accpt");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"Total no. of good candidates");

 
  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);  
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);
  
  // Counter for Normalization
  TString normName="NormalizationCounter";
  AliAnalysisDataContainer *cont = GetOutputSlot(3)->GetContainer();
  
  if(cont)normName=(TString)cont->GetName();
  fCounter = new AliNormalizationCounter(normName.Data());
  fCounter->Init();

  CreateCorrelationObjs(); 

  PostData(1,fOutput);      
  PostData(3,fCounter);      
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDplusCorrelations::UserExec(Option_t */*option*/)
{
  // Do the analysis
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  TClonesArray *array3Prong = 0;

    if(fReco) std::cout << "USING RECONSTRUCTION" << std::endl;
   if(!fReco) std::cout << "USING MC TRUTH" << std::endl;
  
  if(!fMixing){
    if(fSelect==1) cout << "TASK::Correlation with hadrons on SE "<< endl;
    if(fSelect==2) cout << "TASK::Correlation with kaons on SE "<< endl;
    if(fSelect==3) cout << "TASK::Correlation with kzeros on SE "<< endl;
  }  
  if(fMixing){
    if(fSelect==1) cout << "TASK::Correlation with hadrons on ME "<< endl;
    if(fSelect==2) cout << "TASK::Correlation with kaons on ME "<< endl;
    if(fSelect==3) cout << "TASK::Correlation with kzeros on ME "<< endl;
  }
  
  
  if(!aod && AODEvent() && IsStandardAOD()) {

    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      array3Prong=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
    }
  } else if(aod) {
    array3Prong=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
  }

  if(!aod || !array3Prong) {
    printf("AliAnalysisTaskSEDplusCorrelations::UserExec: Charm3Prong branch not found!\n");
    return;
  }

  
  // the AODs with null vertex pointer didn't pass the PhysSel
  
  if(!aod->GetPrimaryVertex()||TMath::Abs(aod->GetMagneticField())<0.001) return;
  fCounter->StoreEvent(aod,fRDCutsAnalysis,fReadMC);
  
  fHistNEvents->Fill(0); // count event
  

  Bool_t isEvSel=fRDCutsAnalysis->IsEventSelected(aod);  
   

  
  PostData(1,fOutput);
  
  if(!isEvSel)return;
  
  fHistNEvents->Fill(1);

  // set PIDResponse for associated tracks
  
  fCorrelator->SetPidAssociated(); 

  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;
  
  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  
  
  TString primTitle = vtx1->GetTitle();
  
  if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0){
  fHistNEvents->Fill(2);
  }
  
  
  
 
  // load MC particles
  if(fReadMC){
    
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      printf("AliAnalysisTaskSEDplusCorrelations::UserExec: MC particles branch not found!\n");
      return;
    }
    
    // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSEDplusCorrelations::UserExec: MC header branch not found!\n");
      return;
    }

       Bool_t isMCeventgood = kFALSE;
       
		
       Int_t eventType = mcHeader->GetEventType();
       Int_t NMCevents = fCuts->GetNofMCEventType();
	       
	
       	for(Int_t k=0; k<NMCevents; k++){
       	Int_t * MCEventType = fCuts->GetMCEventType();

       	if(eventType == MCEventType[k]) isMCeventgood= kTRUE;

	}
	((TH1D*)fOutput->FindObject("EventTypeMC"))->Fill(eventType);
		
		if(NMCevents && !isMCeventgood){
		  
		std::cout << "The MC event " << eventType << " not interesting for this analysis: skipping" << std::endl;
		return;	
		}
    
  }

  //HFCorrelators initialization (for this event)

  fCorrelator->SetAODEvent(aod); // set the event to be processed

  Bool_t correlatorON = fCorrelator->Initialize(); //define the pool for mixing

  if(!correlatorON) {
    AliInfo("AliHFCorrelator didn't initialize the pool correctly or processed a bad event");
    return;
  }
  if(fReadMC) fCorrelator->SetMCArray(arrayMC);

  Int_t n3Prong = 0;
 
  if(fReco) n3Prong = array3Prong->GetEntriesFast();

  
  if(!fReco) n3Prong = arrayMC->GetEntriesFast();  //for MC kine
  
  
   printf("Number of D+->Kpipi: %d and of tracks: %d\n",n3Prong,aod->GetNumberOfTracks());  
  
  Int_t index;
  Int_t pdgDgDplustoKpipi[3]={321,211,211};
  
  Int_t nSelectedloose=0,nSelectedtight=0;

  Double_t ptCand;
  Double_t phiCand;
  Double_t etaCand;
  Double_t invMass = -1;
  
  Int_t nIDs[3] = {-9999999};
     
  Bool_t isDplus = kFALSE;
  
  Int_t labDp = -1;
  
  for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {

    AliAODMCParticle* DplusMC;

    if(fReco){
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);
  
    fHistNEvents->Fill(3);


    
    if(d->Pt()<2.0) continue;  // Dplus below 2.0 not considered.

    

    //      cout<<multipli<<"    multi"<<endl;
    
    

    if(fUseBit && !d->HasSelectionBit(AliRDHFCuts::kDplusCuts)){
      fHistNEvents->Fill(4);
      continue;
    }

    Int_t passTightCuts=fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kAll,aod);
    
      if(!fRDCutsAnalysis->GetIsSelectedCuts()) continue;
    
      fHistNEvents->Fill(5);
    
    
    if(fReadMC){
      labDp = d->MatchToMC(411,arrayMC,3,pdgDgDplustoKpipi);
      
      //  cout<<labDp<<"***"<<endl;
      if(labDp>=0){
	isDplus = kTRUE;

      }
    }    
    
    invMass=d->InvMassDplus();
    Double_t rapid=d->YDplus();
    etaCand = d->Eta();
    phiCand = d->Phi();
    ptCand = d->Pt();
    
    Bool_t isFidAcc=fRDCutsAnalysis->IsInFiducialAcceptance(ptCand,rapid);

    if(isFidAcc){
      fHistNEvents->Fill(6);
      nSelectedloose++;}
      if(!isFidAcc) continue;         
      if(!passTightCuts)continue; 
     
    fHistNEvents->Fill(7);
    
    nSelectedtight++;
    
    
    nIDs[0] = ((AliAODTrack*)d->GetDaughter(0))->GetID();
    nIDs[1] = ((AliAODTrack*)d->GetDaughter(1))->GetID();
    nIDs[2] = ((AliAODTrack*)d->GetDaughter(2))->GetID();
    
    }

    Int_t labDplus=-1;			 	
    
    if(!fReco){
           
      DplusMC = dynamic_cast<AliAODMCParticle*>(arrayMC->At(i3Prong));
    
      if (!DplusMC) {
	AliWarning("Dplus MC Particle not found "); 
	
       continue;
      }
      labDplus = DplusMC->GetLabel();
      
      Int_t PDG =TMath::Abs(DplusMC->PdgCode()); 
      
      if(PDG !=411) continue; // not a Dplus
      ptCand = DplusMC->Pt();
      phiCand = DplusMC->Phi();
      etaCand = DplusMC->Eta();

      Bool_t isFidAccMC =fRDCutsAnalysis->IsInFiducialAcceptance(ptCand,DplusMC->Y());
      
      if(!isFidAccMC)continue;
           
    }				
    
    //cout << "PHI D+ = " << phiCand << endl;	
    // cout << "ETA D+ = " << etaCand << endl;
    // cout << "PT D+ = " << ptCand << endl;
    
   
    
    
        Int_t iPtBin = fRDCutsAnalysis->PtBin(ptCand);				    
    
    if(iPtBin>=0){
      
      index=GetHistoIndex(iPtBin);
      
     cout<<"*****"<<iPtBin<<endl;
      Double_t effTrig;
      if (fTrig){
	
	Double_t multipli  = AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1,1.); 
        
	effTrig = fCuts->GetTrigWeight(ptCand,multipli);

	//	cout<<"*****"<<effTrig<<"  "<<multipli<<"   "<<iPtBin<<endl;
	
      }
      else
	{
	  
	 effTrig = 1.0;
	
	}
     
      
      Double_t trigweight = 1/effTrig ;

      // cout<<"****"<<trigweight<<"  ***"<<endl;
      
      if(fReco && !fMixing){
	
	fMassHist[index]->Fill(invMass); //without trig weight
	
	fMassHistTC[index]->Fill(invMass,trigweight); //with trig wt
      }
     
      if( fReco && fReadMC && isDplus && !fMixing) {
	
	index=GetSignalHistoIndex(iPtBin);
	
	 fMassHistMC[index]->Fill(invMass,trigweight);

         if(labDp>=0){
     
	 AliAODMCParticle *partDp = (AliAODMCParticle*)arrayMC->At(labDp);
	 Int_t DplusSource = CheckOrigin(arrayMC,partDp);
	 
	  if(DplusSource == 4){ // is from charm 
	    
	     fMassHistOrigC[index]->Fill(invMass,trigweight);
	  }

	  if(DplusSource == 5){ // is from beauty
	    
	    fMassHistOrigB[index]->Fill(invMass,trigweight);


	  }	
	 }
      }
    
      
      if(!fReco && fReadMC && !fMixing){

	
	if(labDplus>=0){

	  fMassHist[index]->Fill(1.869);
	
	  AliAODMCParticle *kineDp = (AliAODMCParticle*)arrayMC->At(labDplus);
	  
	  Int_t DplusSource = CheckOrigin(arrayMC,kineDp);
	  
	  
	  if(DplusSource==4){ // is from charm
	    
	     ((TH1F*)fOutput->FindObject(Form("histOrig_Dplus_Bin%d",iPtBin)))->Fill(0.);
	  }

	  if(DplusSource==5){ // is from beauty

	  
	    ((TH1F*)fOutput->FindObject(Form("histOrig_Dplus_Bin%d",iPtBin)))->Fill(1.); 
	  }	
	}
      }

          
	  //Dplus info

      Double_t phiDplus = fCorrelator->SetCorrectPhiRange(phiCand);
      fCorrelator->SetTriggerParticleProperties(ptCand,phiDplus,etaCand);
	    
	  
      
      Double_t ptlead = 0;
      Double_t philead = 0;
      Double_t etalead = 0;
      Double_t refpt = 0;
     
      Int_t LeadSource = 0;
      
      
      
      Bool_t execPool = fCorrelator->ProcessEventPool();
      
      //printf("*************: %d\n",execPool);
      
      if(fMixing && !execPool) {
	AliInfo("Mixed event analysis: pool is not ready");
	continue;
      }

      Int_t NofEventsinPool = 1;
      if(fMixing) {
	NofEventsinPool = fCorrelator->GetNofEventsInPool();

	//cout<<"*********"<<NofEventsinPool<<endl;
      }
      
      
      
      for (Int_t jMix =0; jMix < NofEventsinPool; jMix++){// loop on events in the pool; if it is SE analysis, it stops at 1
	
	Bool_t analyzetracks = fCorrelator->ProcessAssociatedTracks(jMix);
	if(!analyzetracks) {
	  AliInfo("AliHFCorrelator::Cannot process the track array");
	  continue;
	}
	
		
	//Int_t NofTracks = fCorrelator->GetNofTracks();

	//cout<<"***number of tracks****"<<fCorrelator->GetNofTracks()<<endl;
	
	for (Int_t iTrack = 0;iTrack<fCorrelator->GetNofTracks();iTrack++) {	               
	  Bool_t runcorrelation = fCorrelator->Correlate(iTrack);
	  
	  if(!runcorrelation) continue;
	  
	  Double_t DeltaPhi = fCorrelator->GetDeltaPhi();
	 
	  
	  Double_t DeltaEta = fCorrelator->GetDeltaEta();
	  
	  AliReducedParticle* redpart = fCorrelator->GetAssociatedParticle();
	  Double_t phiHad=redpart->Phi();
	  Double_t ptHad=redpart->Pt();
	  Double_t etaHad=redpart->Eta();
	  Int_t label = redpart->GetLabel();
	  Int_t trackid = redpart->GetID();

	  Double_t effi = redpart->GetWeight();
	  Double_t weight = (1/effi)*(trigweight);
	  
	  phiHad = fCorrelator->SetCorrectPhiRange(phiHad);
	  //  cout<<effi<<"******"<<endl;



	  if (!fMixing && fReco){
	    if( trackid == nIDs[0] || trackid == nIDs[1] || trackid == nIDs[2]) continue;  //discard the Dplus daughters
	  }

	  if(fReco && !fReadMC)
	    {
	      FillCorrelations(ptHad,invMass,DeltaPhi,DeltaEta,iPtBin,fSelect,weight);
	    }

	 
	    if(fReco && fReadMC && isDplus){
	   	    
	     if(label<0) continue;
	     
	     AliAODMCParticle* mchad = (AliAODMCParticle*)arrayMC->At(label);
	     
	     if (!mchad) continue;
	     
	     if (!mchad->IsPhysicalPrimary()) continue; //reject the Reco track if correspondent Kine track is not primary
	     
	     if(labDp < 0) continue;
	     
	     AliAODMCParticle *reDp = (AliAODMCParticle*)arrayMC->At(labDp);
	     
	     Int_t RDSource = CheckOrigin(arrayMC,reDp);
	     
	     Int_t MCSource = CheckTrackOrigin(arrayMC,(AliAODMCParticle*)arrayMC->At(label)) ;// check source of associated particle (hadron/kaon/K0)

	     
	    FillMCRCorrelations(ptHad,invMass,DeltaPhi,DeltaEta,iPtBin,MCSource,RDSource,weight);   
	     
	    }
	     	
	      // montecarlo kine
	      if( !fReco && fReadMC){
		
		if(label<0) continue;
		
		if(TMath::Abs(etaHad) > 0.8) continue;
		
		AliAODMCParticle *hadMC = (AliAODMCParticle*)arrayMC->At(label);
		if(!hadMC) continue;
		if (!hadMC->IsPhysicalPrimary()) continue;
		if(labDplus<0) continue;
		  
		AliAODMCParticle *partDp = (AliAODMCParticle*)arrayMC->At(labDplus);
		if(IsDDaughter(partDp,hadMC,arrayMC)) continue;
		  
		
		AliAODMCParticle *trueDp = (AliAODMCParticle*)arrayMC->At(labDplus);
		Int_t DSource = CheckOrigin(arrayMC,trueDp);
	      	
		Int_t PartSource = CheckTrackOrigin(arrayMC,(AliAODMCParticle*)arrayMC->At(label)) ;// check source of associated particle (hadron/kaon/K0)
		
		
		FillMCTruthCorrelations(ptHad,DeltaPhi,DeltaEta,iPtBin,PartSource,DSource,fSelect);   
	      }
	      
	      
	      
	      // For leading particle
	      
	          if (ptHad > refpt) {
		//refpt = ptHad; ptlead = ptHad;

		if(fReadMC){
		 if(label<0) continue;
		 if(!(AliAODMCParticle*)arrayMC->At(label))continue;
		 
		 LeadSource = CheckTrackOrigin(arrayMC,(AliAODMCParticle*)arrayMC->At(label)) ;// check source of associated particle (hadron/kaon/K0)
		}

		philead = phiCand - phiHad;
		etalead = etaCand - etaHad;
		if (philead < (-1)*TMath::Pi()/2) philead += 2*TMath::Pi();
		
		if (philead > 3*TMath::Pi()/2) philead -= 2*TMath::Pi();
		
		refpt = ptHad; ptlead = ptHad;
	      }
	}

	/*if (fReco){
	 Double_t fillSparseLeadDplus[3] = {philead,invMass,etalead};
	}
        else{
	  Double_t fillSparseLeadDplus[3] = {philead,1.87,etalead};
	  }


		if(fReco && !fReadMC && !fMixing){

	((THnSparseF*)fOutput->FindObject(Form("hPhi_Lead_Bin%d",iPtBin)))->Fill(fillSparseLeadDplus,eweight);

	}

	if(fReco && fReadMC && isDplus && !fMixing){

	((THnSparseF*)fOutput->FindObject(Form("hPhi_Lead_Bin%d",iPtBin)))->Fill(fillSparseLeadDplus,eweight);
       
	

	if(LeadSource>=1&&LeadSource<=4){ // is from charm 
	
	  ((THnSparseF*)fOutput->FindObject(Form("hPhi_Lead_From_c_Bin%d",iPtBin)))->Fill(fillSparseLeadDplus,eweight);

	  
	}	
	if(LeadSource>=4&&LeadSource<=8){ // is from beauty
	

	  ((THnSparseF*)fOutput->FindObject(Form("hPhi_Lead_From_b_Bin%d",iPtBin)))->Fill(fillSparseLeadDplus,eweight);


	  
	}	
	
    // from NHF
	if(LeadSource==0){ // is from charm ->D
	  

	  ((THnSparseF*)fOutput->FindObject(Form("hPhi_Lead_From_nhf_Bin%d",iPtBin)))->Fill(fillSparseLeadDplus,eweight);

	}	
	
	fLeadPt[index]->Fill(ptlead);
	}
	
	if(!fReco && fReadMC && !fMixing){

	  index=GetSignalHistoIndex(iPtBin);
	  
	  
	  	((THnSparseF*)fOutput->FindObject(Form("hPhi_Lead_Bin%d",iPtBin)))->Fill(fillSparseLeadDplus,eweight);

	  
	  if(LeadSource>=1&&LeadSource<=4){ // is from charm
 
 ((THnSparseF*)fOutput->FindObject(Form("hPhi_Lead_From_c_Bin%d",iPtBin)))->Fill(fillSparseLeadDplus,eweight);


	    
	  }	
	  if(LeadSource>=4&&LeadSource<=8){ // is from beauty

((THnSparseF*)fOutput->FindObject(Form("hPhi_Lead_From_b_Bin%d",iPtBin)))->Fill(fillSparseLeadDplus,eweight);

	  }	
	  
	  // from NHF
	  if(LeadSource==0){ // is from charm ->D
	    
((THnSparseF*)fOutput->FindObject(Form("hPhi_Lead_From_nhf_Bin%d",iPtBin)))->Fill(fillSparseLeadDplus,eweight);
	    

		}
	  fLeadPt[index]->Fill(ptlead);
	  
	  }*/
	
	  } //jmix
      
	
    } //pt bin 
    
    
    
  }

  if(fMixing){
    Bool_t updated = fCorrelator->PoolUpdate();
	
    if(!updated) AliInfo("Pool was not updated");
  }
  fCounter->StoreCandidates(aod,nSelectedloose,kTRUE);
  fCounter->StoreCandidates(aod,nSelectedtight,kFALSE);
  	
      
  PostData(1,fOutput); 
  PostData(2,fListCuts);
  PostData(3,fCounter);    
  //  PostData(4,fListCutsAsso);    
  return;
}

//________________________________________________________________________

void AliAnalysisTaskSEDplusCorrelations::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEDplusCorrelations: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
 
  return;
}


//________________________________________________________________________
void AliAnalysisTaskSEDplusCorrelations::FillCorrelations(Double_t ptTrack, Double_t mass, Double_t deltaPhi, Double_t deltaEta, Int_t ind, Int_t sel, Double_t eweight) const{
  

  //Filling THnSparse for correlations , only for charged tracks
  
 
  if(sel==1){	  	
      
    if(!fMixing){
      
    Double_t ptLim_Sparse = ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_Bin%d",ind)))->GetAxis(2)->GetXmax(); //all plots have same axes...
    Double_t EtaLim_Sparse = ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_Bin%d",ind)))->GetAxis(3)->GetXmax();
    
    if(ptTrack > ptLim_Sparse) ptTrack = ptLim_Sparse-0.01;   if(deltaEta > EtaLim_Sparse) deltaEta = EtaLim_Sparse-0.01;
    if(deltaEta < -EtaLim_Sparse) deltaEta = -EtaLim_Sparse+0.01;
  
    //variables for filling histos
  
    Double_t fillSparseDplus[4] = {deltaPhi,mass,ptTrack,deltaEta};
   
   

   
      //sparse fill for data (tracks) + weighted
    
    ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_Bin%d",ind)))->Fill(fillSparseDplus,eweight);
          
    }
    
    if(fMixing) {
    

    //variables for filling histos
    
      Double_t fillSparseDplusMix[3] = {deltaPhi,mass,deltaEta};


      ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_Bin%d_EvMix",ind)))->Fill(fillSparseDplusMix,eweight);
    }

  }

  if(sel==2){

  }
  if(sel==3){

  }
  
  return;
}



//________________________________________________________________________



void AliAnalysisTaskSEDplusCorrelations::FillMCTruthCorrelations( Double_t ptTrack,Double_t deltaPhi, Double_t deltaEta, Int_t ind, Int_t mcSource, Int_t origDplus, Int_t sel) const{

   // Filling histos with MC Kine information

     Double_t invm = 1.876;
     
   Double_t ptLim_Sparse = ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_Bin%d",ind)))->GetAxis(2)->GetXmax(); //all plots have same axes...
   Double_t EtaLim_Sparse = ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_Bin%d",ind)))->GetAxis(3)->GetXmax();
    if(ptTrack > ptLim_Sparse) ptTrack = ptLim_Sparse-0.01;
    if(deltaEta > EtaLim_Sparse) deltaEta = EtaLim_Sparse-0.01;
    if(deltaEta < -EtaLim_Sparse) deltaEta = -EtaLim_Sparse+0.01;

    
  
    //variables for filling histos
     
    Double_t fillSparseDplus[4] = {deltaPhi,invm,ptTrack,deltaEta};
    

  if(sel==1){
    
    if(!fMixing){
    
    ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_Bin%d",ind)))->Fill(fillSparseDplus);
    
    if(origDplus==4&&mcSource>=1&&mcSource<=3){ // is from charm 
      
      ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_HFC_Bin%d",ind)))->Fill(fillSparseDplus);
    }
    if(origDplus==5&&mcSource>=4&&mcSource<=8){ // is from beauty 
      
      ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_HFB_Bin%d",ind)))->Fill(fillSparseDplus);
    }
    
    if(mcSource==0){ // is from  NHF
      
      ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_NHF_Bin%d",ind)))->Fill(fillSparseDplus);
    }
    
    }
    
    if(fMixing) {


    //variables for filling histos
      Double_t fillSparseDplusMix[3] = {deltaPhi,invm,deltaEta};
      
      ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_Bin%d_EvMix",ind)))->Fill(fillSparseDplusMix);
    }  
    
  }
  
  if(sel==2){
    
  }
  
  if(sel==3){
    
  }
  
  
  return;
}


void AliAnalysisTaskSEDplusCorrelations::FillMCRCorrelations(Double_t ptTrack, Double_t mass, Double_t deltaPhi, Double_t deltaEta, Int_t ind, Int_t mcSource, Int_t origDplus, Double_t eweight) const{
  
  // Filling histos with MC information

    
   Double_t ptLim_Sparse = ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_Bin%d",ind)))->GetAxis(2)->GetXmax(); //all plots have same axes...
  Double_t EtaLim_Sparse = ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_Bin%d",ind)))->GetAxis(3)->GetXmax();
    if(ptTrack > ptLim_Sparse) ptTrack = ptLim_Sparse-0.01;
     if(deltaEta > EtaLim_Sparse) deltaEta = EtaLim_Sparse-0.01;
    if(deltaEta < -EtaLim_Sparse) deltaEta = -EtaLim_Sparse+0.01;


  if(!fMixing){
 
    Double_t fillSparseDplus[4] = {deltaPhi,mass,ptTrack,deltaEta};
    

    
    ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_Bin%d",ind)))->Fill(fillSparseDplus,eweight);
  
     if(origDplus==4&&mcSource>=1&&mcSource<=3){ // is from charm 

    ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_HFC_Bin%d",ind)))->Fill(fillSparseDplus,eweight);
     }
     if(origDplus==5&&mcSource>=4&&mcSource<=8){ // is from beauty 

      ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_HFB_Bin%d",ind)))->Fill(fillSparseDplus,eweight);
     }

     if(mcSource==0){ // is from  NHF

      ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_NHF_Bin%d",ind)))->Fill(fillSparseDplus,eweight);
     }
  }
         
     if(fMixing) {
    
    //variables for filling histos
    Double_t fillSparseDplusMix[3] = {deltaPhi,mass,deltaEta};

    ((THnSparseF*)fOutput->FindObject(Form("hPhi_Charg_Bin%d_EvMix",ind)))->Fill(fillSparseDplusMix,eweight);
     }


  return;
}

Bool_t AliAnalysisTaskSEDplusCorrelations::IsDDaughter(AliAODMCParticle* d, AliAODMCParticle* track, TClonesArray* mcArray) const {

  //Daughter removal in MCKine case
  Bool_t isDaughter = kFALSE;
  Int_t labelDplus = d->GetLabel();
  
  Int_t mother = track->GetMother();

  
  while (mother > 0){
    AliAODMCParticle* mcMoth = dynamic_cast<AliAODMCParticle*>(mcArray->At(mother)); 
    if (mcMoth){
      if (mcMoth->GetLabel() == labelDplus) isDaughter = kTRUE;
      mother = mcMoth->GetMother(); //goes back to last generation
    } else{
      AliError(" mother particle not found!");
      break;
    }
  }

  return isDaughter;
}


Int_t AliAnalysisTaskSEDplusCorrelations::CheckTrackOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const {		
  //
  // checks on particle (#) origin:
  // 0) Not HF
  // 1) D->#
  // 2) D->X->#
  // 3) c hadronization
  // 4) B->#
  // 5) B->X-># (X!=D)
  // 6) B->D->#
  // 7) B->D->X->#
  // 8) b hadronization
  //
  if(fDebug>2) printf("AliAnalysisTaskSEDplusCorrelations::CheckTrkOrigin() \n");
	
  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPartCandidate->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  Bool_t isDdaugh=kFALSE;
  Bool_t isDchaindaugh=kFALSE;
  Bool_t isBdaugh=kFALSE;
  Bool_t isBchaindaugh=kFALSE;
  Bool_t isQuarkFound=kFALSE;

  while (mother > 0){
    istep++;
    AliAODMCParticle* mcMoth = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcMoth){
      pdgGranma = mcMoth->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isBchaindaugh=kTRUE;
        if(istep==1) isBdaugh=kTRUE;
      }
      if ((abspdgGranma > 400 && abspdgGranma < 500) || (abspdgGranma > 4000 && abspdgGranma < 5000)){
	isDchaindaugh=kTRUE;
        if(istep==1) isDdaugh=kTRUE;
      }
      if(abspdgGranma==4 || abspdgGranma==5) {isQuarkFound=kTRUE; if(abspdgGranma==5) isFromB = kTRUE;}
      mother = mcMoth->GetMother();
    }else{
      AliError("Failed casting the mother particle!");
      break;
    }
  }

  //decides what to return based on the flag status
  if(isQuarkFound) {
    if(!isFromB) {  //charm
      if(isDdaugh) return 1; //charm immediate
      else if(isDchaindaugh) return 2; //charm chain
      else return 3; //charm hadronization
    }
    else { //beauty
      if(isBdaugh) return 4; //b immediate
      else if(isBchaindaugh) { //b chain
        if(isDchaindaugh) {
          if(isDdaugh) return 6; //d immediate
          return 7; //d chain
          }
        else return 5; //b, not d
      }
      else return 8; //b hadronization
    }
  }
  else return 0; //no HF quark
}

Int_t AliAnalysisTaskSEDplusCorrelations::CheckOrigin(TClonesArray* arrayMC,  AliAODMCParticle *mcPartCandidate) const {		
  //
  // checking whether the mother of the particles come from a charm or a bottom quark
  //
	
  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPartCandidate->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  Bool_t isQuarkFound=kFALSE;
  while (mother >0 ){
    istep++;
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma){
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isFromB=kTRUE;
      }
      if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
      mother = mcGranma->GetMother();
    }else{
      AliError("Failed casting the mother particle!");
      break;
    }
  }
  
  if(isFromB) return 5;
  else return 4;
}


void AliAnalysisTaskSEDplusCorrelations::CreateCorrelationObjs() {
//

  TString namePlot = "";
  Int_t nbinsmass=GetNBinsHistos();
  //These for limits in THnSparse (one per bin, same limits). 
  //Vars: DeltaPhi, InvMass, PtTrack, Displacement, DeltaEta
  Int_t nBinsPhi[4] = {32,nbinsmass,30,16};
  Double_t binMinPhi[4] = {-TMath::Pi()/2.+TMath::Pi()/32.,fLowmasslimit,0.,-1.6};  //is the minimum for all the bins
  Double_t binMaxPhi[4] = {3*TMath::Pi()/2.+TMath::Pi()/32.,fUpmasslimit,30.0,1.6};  //is the maximum for all the bins

  //Vars: DeltaPhi, InvMass, DeltaEta
  Int_t nBinsMix[3] = {32,nbinsmass,16};
  Double_t binMinMix[3] = {-TMath::Pi()/2.+TMath::Pi()/32.,fLowmasslimit,-1.6};  //is the minimum for all the bins
  Double_t binMaxMix[3] = {3*TMath::Pi()/2.+TMath::Pi()/32.,fUpmasslimit,1.6};  //is the maximum for all the bins

  for(Int_t i=0;i<fNPtBins;i++){

    if(!fMixing) {
      //THnSparse plots: correlations for various invariant mass (MC and data)
      namePlot="hPhi_Charg_Bin";
      namePlot+=i;

      THnSparseF *hPhi = new THnSparseF(namePlot.Data(), "Azimuthal correlation; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
      hPhi->Sumw2();
      fOutput->Add(hPhi);

      namePlot="hPhi_Kaon_Bin";
      namePlot+=i;

      THnSparseF *hPhiK = new THnSparseF(namePlot.Data(), "Azimuthal correlation; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
      hPhiK->Sumw2();
      fOutput->Add(hPhiK);

      namePlot="hPhi_K0_Bin";
      namePlot+=i;

      THnSparseF *hPhiK0 = new THnSparseF(namePlot.Data(), "Azimuthal correlation; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
      hPhiK0->Sumw2();
      fOutput->Add(hPhiK0);
  
      //histos for c/b origin for D+ (MC only)
      if (fReadMC) {

        //generic origin for tracks
        namePlot="hPhi_Charg_HFC_Bin";
        namePlot+=i;

        THnSparseF *hPhiH_c = new THnSparseF(namePlot.Data(), "Azimuthal correlation - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiH_c->Sumw2();
        fOutput->Add(hPhiH_c);

        namePlot="hPhi_Kaon_HFC_Bin";
        namePlot+=i;

        THnSparseF *hPhiK_c = new THnSparseF(namePlot.Data(), "Azimuthal correlation - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiK_c->Sumw2();
        fOutput->Add(hPhiK_c);

        namePlot="hPhi_K0_HFC_Bin";
        namePlot+=i;

        THnSparseF *hPhiK0_c = new THnSparseF(namePlot.Data(), "Azimuthal correlation - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiK0_c->Sumw2();
        fOutput->Add(hPhiK0_c);
  
        namePlot="hPhi_Charg_HFB_Bin";
        namePlot+=i;

        THnSparseF *hPhiH_b = new THnSparseF(namePlot.Data(), "Azimuthal correlation - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiH_b->Sumw2();
        fOutput->Add(hPhiH_b);

        namePlot="hPhi_Kaon_HFB_Bin";
        namePlot+=i;

        THnSparseF *hPhiK_b = new THnSparseF(namePlot.Data(), "Azimuthal correlation - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiK_b->Sumw2();
        fOutput->Add(hPhiK_b);

        namePlot="hPhi_K0_HFB_Bin";
        namePlot+=i;

        THnSparseF *hPhiK0_b = new THnSparseF(namePlot.Data(), "Azimuthal correlation - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiK0_b->Sumw2();
        fOutput->Add(hPhiK0_b);

	namePlot="hPhi_Charg_NHF_Bin";
        namePlot+=i;

        THnSparseF *hPhiH_NHF = new THnSparseF(namePlot.Data(), "Azimuthal correlation - Non HF; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiH_NHF->Sumw2();
        fOutput->Add(hPhiH_NHF);

        namePlot="hPhi_Kaon_NHF_Bin";
        namePlot+=i;

        THnSparseF *hPhiK_Non = new THnSparseF(namePlot.Data(), "Azimuthal correlation - Non HF; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiK_Non->Sumw2();
        fOutput->Add(hPhiK_Non);

        namePlot="hPhi_K0_NHF_Bin";
        namePlot+=i;

        THnSparseF *hPhiK0_Non = new THnSparseF(namePlot.Data(), "Azimuthal correlation - Non HF; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiK0_Non->Sumw2();
        fOutput->Add(hPhiK0_Non);
      }

      //leading hadron correlations
      namePlot="hPhi_Lead_Bin";
      namePlot+=i;

      THnSparseF *hCorrLead = new THnSparseF(namePlot.Data(), "Leading particle correlations; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
      hCorrLead->Sumw2();
      fOutput->Add(hCorrLead);

      if (fReadMC) {
        namePlot="hPhi_Lead_From_c_Bin";
        namePlot+=i;

        THnSparseF *hCorrLead_c = new THnSparseF(namePlot.Data(), "Leading particle correlations - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
        hCorrLead_c->Sumw2();
        fOutput->Add(hCorrLead_c);
  
        namePlot="hPhi_Lead_From_b_Bin";
        namePlot+=i;
  
        THnSparseF *hCorrLead_b = new THnSparseF(namePlot.Data(), "Leading particle correlations - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
        hCorrLead_b->Sumw2();
        fOutput->Add(hCorrLead_b);
  


        namePlot="hPhi_Lead_From_nhf_Bin";
        namePlot+=i;
  
        THnSparseF *hCorrLead_Non = new THnSparseF(namePlot.Data(), "Leading particle correlations - Non HF; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
        hCorrLead_Non->Sumw2();
        fOutput->Add(hCorrLead_Non);
      }
      
           

    //pT distribution histos
    namePlot = "hist_Pt_Charg_Bin"; namePlot+=i;
    TH1F *hPtHad = new TH1F(namePlot.Data(), "Charged track pT (in D+ events); p_{T} (GeV/c)",240,0.,20.);
    hPtHad->SetMinimum(0);
    fOutput->Add(hPtHad);

    namePlot = "hist_Pt_Kcharg_Bin"; namePlot+=i;
    TH1F *hPtH = new TH1F(namePlot.Data(), "Hadrons pT (in D+ events); p_{T} (GeV/c)",240,0.,12.);
    hPtH->SetMinimum(0);
    fOutput->Add(hPtH);

    namePlot = "hist_Pt_K0_Bin"; namePlot+=i;
    TH1F *hPtK = new TH1F(namePlot.Data(), "Kaons pT (in D+ events); p_{T} (GeV/c)",240,0.,12.);
    hPtK->SetMinimum(0);
    fOutput->Add(hPtK);

     
    }

    if(fMixing) {
      //THnSparse plots for event mixing!


      namePlot="hPhi_K0_Bin";
      namePlot+=i;namePlot+="_EvMix";

      THnSparseF *hPhiK_EvMix = new THnSparseF(namePlot.Data(), "Az. corr. EvMix; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
      hPhiK_EvMix->Sumw2();
      fOutput->Add(hPhiK_EvMix);

      namePlot="hPhi_Kcharg_Bin";
      namePlot+=i;namePlot+="_EvMix";
  
      THnSparseF *hPhiH_EvMix = new THnSparseF(namePlot.Data(), "Az. corr. EvMix; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
      hPhiK_EvMix->Sumw2();
      fOutput->Add(hPhiH_EvMix);

      namePlot="hPhi_Charg_Bin";
      namePlot+=i;namePlot+="_EvMix";

      THnSparseF *hPhiC_EvMix = new THnSparseF(namePlot.Data(), "Az. corr. EvMix; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
      hPhiC_EvMix->Sumw2();
      fOutput->Add(hPhiC_EvMix);
    }
  }

  //out of bin loop
  if(!fMixing) {
    TH1F *hCountC = new TH1F("hist_Count_Charg", "Charged track counter; # Tracks",100,0.,100.);
    hCountC->SetMinimum(0);
    fOutput->Add(hCountC);

    TH1F *hCountH = new TH1F("hist_Count_Kcharg", "Hadrons counter; # Tracks",20,0.,20.);
    hCountH->SetMinimum(0);
    fOutput->Add(hCountH);

    TH1F *hCountK = new TH1F("hist_Count_K0", "Kaons counter; # Tracks",20,0.,20.);
    hCountK->SetMinimum(0);
    fOutput->Add(hCountK);
  }

  if (fReadMC) {
    TH1D *hEventTypeMC = new TH1D("EventTypeMC","EventTypeMC",100,-0.5,99.5);
    fOutput->Add(hEventTypeMC); 
  }

  // if (fFillGlobal) { //all-events plots
    //pt distributions

   if(!fMixing) {
    TH1F *hPtHAll = new TH1F("hist_Pt_Charg_AllEv", "Charged track pT (All); p_{T} (GeV/c)",240,0.,12.);
    hPtHAll->SetMinimum(0);
    fOutput->Add(hPtHAll);

    TH1F *hPtKAll = new TH1F("hist_Pt_Kaons_AllEv", "Kaons pT (All); p_{T} (GeV/c)",240,0.,12.);
    hPtKAll->SetMinimum(0);
    fOutput->Add(hPtKAll);

    TH1F *hPtK0All = new TH1F("hist_Pt_K0_AllEv", "K0 pT (All); p_{T} (GeV/c)",240,0.,12.);
    hPtK0All->SetMinimum(0);
    fOutput->Add(hPtK0All);

   
      //phi distributions
      TH1F *hPhiDistHAll = new TH1F("hist_PhiDistr_Charg", "Charged track phi distr. (All); p_{T} (GeV/c)",64,0,6.283);
      hPhiDistHAll->SetMinimum(0);
      fOutput->Add(hPhiDistHAll);

      TH1F *hPhiDistKAll = new TH1F("hist_PhiDistr_Kcharg", "Kaons phi distr. (All); p_{T} (GeV/c)",64,0,6.283);
      hPhiDistKAll->SetMinimum(0);
      fOutput->Add(hPhiDistKAll);

      TH1F *hPhiDistK0All = new TH1F("hist_PhiDistr_K0", "K0 phi distr. (All); p_{T} (GeV/c)",64,0,6.283);
      hPhiDistK0All->SetMinimum(0);
      fOutput->Add(hPhiDistK0All);

      TH1F *hPhiDistDAll = new TH1F("hist_PhiDistr_Dplus", "D^{+} phi distr. (All); p_{T} (GeV/c)",64,0,6.283);
      hPhiDistDAll->SetMinimum(0);
      fOutput->Add(hPhiDistDAll);
   }
    //}

  //for MC analysis only
  for(Int_t i=0;i<fNPtBins;i++) {

    if (fReadMC && !fMixing) {

      
      //origin of tracks histos
      namePlot="histOrig_Had_Bin";  namePlot+=i;
      TH1F *hOrigin_had = new TH1F(namePlot.Data(), "Origin of associated tracksin MC",10,-0.5,9.5);
      hOrigin_had->SetMinimum(0);
      
      hOrigin_had->GetXaxis()->SetBinLabel(1,"All ");
      hOrigin_had->GetXaxis()->SetBinLabel(2," from Heavy flavour");
      hOrigin_had->GetXaxis()->SetBinLabel(3," from c->D");
      hOrigin_had->GetXaxis()->SetBinLabel(4," from b->D");
      hOrigin_had->GetXaxis()->SetBinLabel(5," from b->B");
      hOrigin_had->GetXaxis()->SetBinLabel(6," from NHF");
      if(fReadMC) fOutput->Add(hOrigin_had);

      namePlot="histOrig_Kaon_Bin";  namePlot+=i;
      TH1F *hOrigin_kaon = new TH1F(namePlot.Data(), "Origin of kaons in MC",10,-0.5,9.5);
      hOrigin_kaon->SetMinimum(0);
      hOrigin_kaon->GetXaxis()->SetBinLabel(1,"All Kaons");
      hOrigin_kaon->GetXaxis()->SetBinLabel(2,"Kaons from Heavy flavour");
      hOrigin_kaon->GetXaxis()->SetBinLabel(3,"Kaons from Charm");
      hOrigin_kaon->GetXaxis()->SetBinLabel(4,"Kaons from Beauty");
      hOrigin_kaon->GetXaxis()->SetBinLabel(5,"Kaons from NHF");
      if(fReadMC) fOutput->Add(hOrigin_kaon);
      

      namePlot="histOrig_K0_Bin";  namePlot+=i;
      TH1F *hOrigin_K0 = new TH1F(namePlot.Data(), "Origin of Kshorts in MC",10,-0.5,9.5);
      hOrigin_K0->SetMinimum(0);

      hOrigin_K0->GetXaxis()->SetBinLabel(1,"All K0s");
      hOrigin_K0->GetXaxis()->SetBinLabel(2,"K0s from Heavy flavour");
      hOrigin_K0->GetXaxis()->SetBinLabel(3,"K0s from Charm");
      hOrigin_K0->GetXaxis()->SetBinLabel(4,"K0s from Beauty");
      hOrigin_K0->GetXaxis()->SetBinLabel(5,"K0s from NHF");
      if(fReadMC) fOutput->Add(hOrigin_K0);
    }

    if (fReadMC) {
      //origin of Dplus histos
      namePlot="histOrig_Dplus_Bin";  namePlot+=i;
      TH1F *hOrigin_Dplus = new TH1F(namePlot.Data(),"Origin of D+ in MC",2,0.,2.);
      hOrigin_Dplus->SetMinimum(0);
      hOrigin_Dplus->GetXaxis()->SetBinLabel(1,"From c");
      hOrigin_Dplus->GetXaxis()->SetBinLabel(2,"From b");
      fOutput->Add(hOrigin_Dplus);
    }
  }
}









