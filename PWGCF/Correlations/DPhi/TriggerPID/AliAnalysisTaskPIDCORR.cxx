/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

//Modified.....

#include "AliAnalysisTaskPIDCORR.h"

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliPIDResponse.h"
#include "AliPID.h"
#include "AliAODInputHandler.h"
#include "AliVEvent.h"
#include "AliVParticle.h"

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "TParticle.h"
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"

#include "AliEventPoolManager.h"


ClassImp(AliAnalysisTaskPIDCORR)
ClassImp(AliPIDCorrParticle)

//________________________________________________________________________
AliAnalysisTaskPIDCORR::AliAnalysisTaskPIDCORR() // All data members should be initialised here
:AliAnalysisTaskSE(),
  fAOD(0),
  fAODVertex(0),
  fPIDResponse(0),
  fOutputList(0),
  fHistPt(0),
  fHistdEdx(0),
  fTriggerPhiAll(0),
  fTriggerPhiPion(0),
  fTriggerPhiKaonProton(0),
  fTriggerEtaAll(0),
  fTriggerEtaPion(0),
  fTriggerEtaKaonProton(0),
  fAssoPhi(0),
  fAssoEta(0),
  fGlobalTracks(0),
  fArrayMC(0),
  fPoolMgr(0),
  EffEtaTrigPr(1),
  EffEtaTrigPi(1),
  EffEtaTrigAll(1)
  
{
    
  
 

  	for(int i=0;i<10;i++)
	{
	fDihadronCorrelation[i]=NULL;
	fDihadronCorrelationPion[i]=NULL;
	fDihadronCorrelationProtonKaon[i]=NULL;
	fHistNSigmaTPCPion[i]=NULL;
	fMixedEvent[i]=NULL;
	fMixedPion[i]=NULL;
	fMixedProtonKaon[i]=NULL;
	}

	for(int j=0;j<3;j++)
	{
	fHistNSAll[j]=NULL;
	fHistNSPion[j]=NULL;
	fHistNSProton[j]=NULL;
	fHistASAll[j]=NULL;
	fHistASPion[j]=NULL;
	fHistASProton[j]=NULL;
	fHistBgAll[j]=NULL;
	fHistBgPion[j]=NULL;
	fHistBgProton[j]=NULL;
	fHistBulkAll[j]=NULL;
	fHistBulkPion[j]=NULL;
	fHistBulkProton[j]=NULL;
	
	
	}
  
  
}

//________________________________________________________________________
AliAnalysisTaskPIDCORR::AliAnalysisTaskPIDCORR(const char *name) // All data members should be initialised here
  :AliAnalysisTaskSE(name),
   fAOD(0),
  fAODVertex(0),
  fPIDResponse(0),
  fOutputList(0),
  fHistPt(0),
  fHistdEdx(0),
  fTriggerPhiAll(0),
  fTriggerPhiPion(0),
  fTriggerPhiKaonProton(0),
  fTriggerEtaAll(0),
  fTriggerEtaPion(0),
  fTriggerEtaKaonProton(0),
  fAssoPhi(0),
  fAssoEta(0),
  fGlobalTracks(0),
  fArrayMC(0),
  fPoolMgr(0),
  EffEtaTrigPr(1),
  EffEtaTrigPi(1),
  EffEtaTrigAll(1)
  
{

   

	for(int i=0;i<10;i++)
	{
	fDihadronCorrelation[i]=NULL;
	fDihadronCorrelationPion[i]=NULL;
	fDihadronCorrelationProtonKaon[i]=NULL;
	fHistNSigmaTPCPion[i]=NULL;
	fMixedEvent[i]=NULL;
	fMixedPion[i]=NULL;
	fMixedProtonKaon[i]=NULL;
	}

	for(int j=0;j<3;j++)
	{
	fHistNSAll[j]=NULL;
	fHistNSPion[j]=NULL;
	fHistNSProton[j]=NULL;
	fHistASAll[j]=NULL;
	fHistASPion[j]=NULL;
	fHistASProton[j]=NULL;
	fHistBgAll[j]=NULL;
	fHistBgPion[j]=NULL;
	fHistBgProton[j]=NULL;
	fHistBulkAll[j]=NULL;
	fHistBulkPion[j]=NULL;
	fHistBulkProton[j]=NULL;
	
	
	}
    
  
  // The last in the above list should not have a comma after it
  
  // Constructor
  // Define input and output slots here (never in the dummy constructor)
  // Input slot #0 works with a TChain - it is connected to the default input container
  // Output slot #1 writes into a TH1 container
  
  DefineOutput(1, TList::Class());                                        // for output list
}

//________________________________________________________________________
AliAnalysisTaskPIDCORR::~AliAnalysisTaskPIDCORR()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if (fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutputList;
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskPIDCORR::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (on the worker node)

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager(); 
  AliInputEventHandler *inputHandler = (AliInputEventHandler*)man->GetInputEventHandler(); 
  fPIDResponse = inputHandler->GetPIDResponse();
  
  fOutputList = new TList();
  fOutputList->SetOwner();  // IMPORTANT!

 
  
  
  

  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 15, 0.1, 3.1);
  fHistdEdx = new TH2F("fHistdEdx","TPC signal distribution",1000,0.,50.,1000,0,500);

  
  TString HistName;
  
	for(int i=0;i<10;i++)
	{
	HistName="fDihadronCorrelation";HistName+=i;
	fDihadronCorrelation[i]= new TH2F(HistName.Data(),"#delta#eta-#delta#phi correlation",36,-TMath::Pi()/2,3.*TMath::Pi()/2,32,-1.6,1.6);
	fDihadronCorrelation[i]->Sumw2();
	fOutputList->Add(fDihadronCorrelation[i]);

	HistName="fDihadronCorrelationPion";HistName+=i;
	fDihadronCorrelationPion[i]= new TH2F(HistName.Data(),"#delta#eta-#delta#phi correlation",36,-TMath::Pi()/2,3.*TMath::Pi()/2,32,-1.6,1.6);
	fDihadronCorrelationPion[i]->Sumw2();
	fOutputList->Add(fDihadronCorrelationPion[i]);

	HistName="fDihadronCorrelationProtonKaon";HistName+=i;
	fDihadronCorrelationProtonKaon[i]= new TH2F(HistName.Data(),"#delta#eta-#delta#phi correlation",36,-TMath::Pi()/2,3.*TMath::Pi()/2,32,-1.6,1.6);
	fDihadronCorrelationProtonKaon[i]->Sumw2();
	fOutputList->Add(fDihadronCorrelationProtonKaon[i]);

	HistName="fMixedEvent";HistName+=i;
	fMixedEvent[i]= new TH2F(HistName.Data(),"#delta#eta-#delta#phi correlation",36,-TMath::Pi()/2,3.*TMath::Pi()/2,32,-1.6,1.6);
	fMixedEvent[i]->Sumw2();
	fOutputList->Add(fMixedEvent[i]);

	HistName="fMixedPion";HistName+=i;
	fMixedPion[i]= new TH2F(HistName.Data(),"#delta#eta-#delta#phi correlation",36,-TMath::Pi()/2,3.*TMath::Pi()/2,32,-1.6,1.6);
	fMixedPion[i]->Sumw2();
	fOutputList->Add(fMixedPion[i]);

	HistName="fMixedProtonKaon";HistName+=i;
	fMixedProtonKaon[i]= new TH2F(HistName.Data(),"#delta#eta-#delta#phi correlation",36,-TMath::Pi()/2,3.*TMath::Pi()/2,32,-1.6,1.6);
	fMixedProtonKaon[i]->Sumw2();
	fOutputList->Add(fMixedProtonKaon[i]);
	
	HistName="fHistNSigmaTPCPion";HistName+=i;
	fHistNSigmaTPCPion[i]= new TH1F(HistName.Data(),"NSigma distribution for all particles",200,-10,10);
	fHistNSigmaTPCPion[i]->Sumw2();
	fOutputList->Add(fHistNSigmaTPCPion[i]);
	}

	for(Int_t l=0;l<3;l++)
	{
	HistName="fHistNSPtSpectraAll";HistName+=l;
	fHistNSAll[l]= new TH1F(HistName.Data(),"Pt_spectra of Associated",15,0.5,3.5);
	fHistNSAll[l]->Sumw2();
	fOutputList->Add(fHistNSAll[l]);
	
	HistName="fHistNSPtSpectraPion";HistName+=l;
	fHistNSPion[l]= new TH1F(HistName.Data(),"Pt_spectra of Associated",15,0.5,3.5);
	fHistNSPion[l]->Sumw2();
	fOutputList->Add(fHistNSPion[l]);
	
	HistName="fHistNSPtSpectraProton";HistName+=l;
	fHistNSProton[l]= new TH1F(HistName.Data(),"Pt_spectra of Associated",15,0.5,3.5);
	fHistNSProton[l]->Sumw2();
	fOutputList->Add(fHistNSProton[l]);

	HistName="fHistASPtSpectraAll";HistName+=l;
	fHistASAll[l]= new TH1F(HistName.Data(),"Pt_spectra of Associated",15,0.5,3.5);
	fHistASAll[l]->Sumw2();
	fOutputList->Add(fHistASAll[l]);

	HistName="fHistASPtSpectraPion";HistName+=l;
	fHistASPion[l]= new TH1F(HistName.Data(),"Pt_spectra of Associated",15,0.5,3.5);
	fHistASPion[l]->Sumw2();
	fOutputList->Add(fHistASPion[l]);

	HistName="fHistASPtSpectraProton";HistName+=l;
	fHistASProton[l]= new TH1F(HistName.Data(),"Pt_spectra of Associated",15,0.5,3.5);
	fHistASProton[l]->Sumw2();
	fOutputList->Add(fHistASProton[l]);

	HistName="fHistBgPtSpectraAll";HistName+=l;
	fHistBgAll[l]= new TH1F(HistName.Data(),"Pt_spectra of Associated",15,0.5,3.5);
	fHistBgAll[l]->Sumw2();
	fOutputList->Add(fHistBgAll[l]);

	HistName="fHistBgPtSpectraPion";HistName+=l;
	fHistBgPion[l]= new TH1F(HistName.Data(),"Pt_spectra of Associated",15,0.5,3.5);
	fHistBgPion[l]->Sumw2();
	fOutputList->Add(fHistBgPion[l]);

	HistName="fHistBgPtSpectraProton";HistName+=l;
	fHistBgProton[l]= new TH1F(HistName.Data(),"Pt_spectra of Associated",15,0.5,3.5);
	fHistBgProton[l]->Sumw2();
	fOutputList->Add(fHistBgProton[l]);

	

	HistName="fHistBulkPtSpectraAll";HistName+=l;
	fHistBulkAll[l]= new TH1F(HistName.Data(),"Pt_spectra of Associated",15,0.5,3.5);
	fHistBulkAll[l]->Sumw2();
	fOutputList->Add(fHistBulkAll[l]);

	HistName="fHistBulkPtSpectraPion";HistName+=l;
	fHistBulkPion[l]= new TH1F(HistName.Data(),"Pt_spectra of Associated",15,0.5,3.5);
	fHistBulkPion[l]->Sumw2();
	fOutputList->Add(fHistBulkPion[l]);

	HistName="fHistBulkPtSpectraProton";HistName+=l;
	fHistBulkProton[l]= new TH1F(HistName.Data(),"Pt_spectra of Associated",15,0.5,3.5);
	fHistBulkProton[l]->Sumw2();
	fOutputList->Add(fHistBulkProton[l]);

	}

	fTriggerPhiAll = new TH1F("fTriggerPhiAll","Phi distribution of  All trigger particles",340,0.0,2*TMath::Pi());
	fTriggerPhiPion = new TH1F("fTriggerPhiPion","Phi distribution of Pion trigger particles",340,0.0,2*TMath::Pi());
	fTriggerPhiKaonProton = new TH1F("fTriggerPhiKaonProton","Phi distribution of Kaon & Proton trigger particles",340,0.0,2*TMath::Pi());

	fTriggerEtaAll = new TH1F("fTriggerEtaAll","Eta distribution of  All trigger particles",16,-0.8,0.8);
	fTriggerEtaPion = new TH1F("fTriggerEtaPion","Eta distribution of  Pion trigger particles",16,-0.8,0.8);
	fTriggerEtaKaonProton = new TH1F("fTriggerEtaKaonProton","Eta distribution of  Kaon & Proton trigger particles",16,-0.8,0.8);

	fAssoPhi = new TH1F("fAssoPhi","Phi distribution of  Associated particles",340,0.0,2*TMath::Pi());
	fAssoEta = new TH1F("fAssoEta","Eta distribution of  Associated particles",16,-0.8,0.8); 



  fOutputList->Add(fTriggerPhiAll);
  fOutputList->Add(fTriggerPhiPion);
  fOutputList->Add(fTriggerPhiKaonProton);
  fOutputList->Add(fTriggerEtaAll);
  fOutputList->Add(fTriggerEtaPion);
  fOutputList->Add(fTriggerEtaKaonProton);
  fOutputList->Add(fAssoPhi);
  fOutputList->Add(fAssoEta);
  fOutputList->Add(fHistPt);
  fOutputList->Add(fHistdEdx);

	
	
  

  
          
     
//Mixing
DefineEventPool();
	
  
  PostData(1, fOutputList);              // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
Bool_t AliAnalysisTaskPIDCORR :: SelectEvent(AliAODVertex* vertex)

{

	//
	// Event Selection.
	//
	
	
	Double_t primVtx[3];
	vertex->GetXYZ(primVtx);
	if (TMath::Sqrt(primVtx[0]*primVtx[0] + primVtx[1]*primVtx[1])>1. || TMath::Abs(primVtx[2])>10.) {
		
		return kFALSE;
	}
	
	return kTRUE;

}
//_______________________________________________________________________
Int_t AliAnalysisTaskPIDCORR::ClassifyTrack(AliAODTrack* track)
{
//Int_t Classification=999;
Double_t pt = track->Pt();
Double_t eta = track->Eta();
//Double_t phi = track->Phi();

if (!(track->TestFilterMask(1<<7))) {
		
		return 0;
	}

if(pt< 0.2)
	{
	 return 0;
}
	
if (TMath::Abs(eta)>=0.8) {
		
		return 0;
    }

//DCA cut

/*if(PtDependentCut)
{


}*/

if(track->Charge()==0) return 0;

AliAODTrack *globaltrack=GetGlobalTrack(track->GetID());
if(!globaltrack) return  0;


Float_t dxy, dz;
		  
dxy = track->DCA();
dz  = track->ZAtDCA();
//cout<<dxy<<"   "<<dz<<endl;
if(TMath ::Abs(dxy) >2.4 || TMath ::Abs(dz)>3.2) return 0;

Double_t chi2ndf = track->Chi2perNDF();
if(chi2ndf>4.0) return 0;

Double_t nclus = track->GetTPCClusterInfo(2,1);
if(nclus<80) return 0;

//select primary:
if(track->GetType() != AliAODTrack::kPrimary) return 0;

return 1;
}
//________________________________________________________________
void AliAnalysisTaskPIDCORR::FillGlobalTracksArray() {

	
	
	if (!fAOD) {
        //if (fVerbose>0) cout << "AliAnalysisTaskDiHadronPID::FillGlobalTracksArray -> ERROR: fAODEvent not set." << endl;
		return;
	}

	
	
	fGlobalTracks = new TObjArray();
	AliAODTrack* track = 0x0;
		
	for (Int_t iTrack = 0; iTrack < fAOD->GetNumberOfTracks(); iTrack++) {
		
		track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
		if(!track) AliFatal("Not a standard AOD");
		
		// I.e., if it does NOT pass the filtermask.
		if (!(track->TestFilterMask(1<<7))) {
            if (track->GetID()>-1) fGlobalTracks->AddAtAndExpand(track,track->GetID());
            //cout<<"Track ID: "<<track->GetID()<<" Partner ID: "<<(-track->GetID()-1)<<endl;
		}
        
	}
	}
//_______________________________________________________________________

AliAODTrack* AliAnalysisTaskPIDCORR::GetGlobalTrack(Int_t trackID) {
	
	
	
	AliAODTrack* partner = 0x0;
	    
    partner = (AliAODTrack*)(fGlobalTracks->At(-trackID-1));
	  
	
	return partner;
	
}


//________________________________________________________________________
void AliAnalysisTaskPIDCORR::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  // Post output data.
	
  AliVEvent *event = InputEvent();
  if (!event) {  return; }

 
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
   // printf("ERROR: fAOD not available\n");
    return;
  }




Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
    if(!isSelected) return;

fAODVertex = fAOD->GetPrimaryVertex();
	if (!fAODVertex) {
		//if (fVerbose>0) cout << "AliAnalysisTaskDiHadronPID::UserExec -> ERROR: No AliAODVertex pointer could be created." << endl;
		return;
	}

//select vertex
if(!SelectEvent(fAODVertex)) return;

/*fHistVx->Fill(fAODVertex->GetX());
fHistVy->Fill(fAODVertex->GetY());
fHistVz->Fill(fAODVertex->GetZ());*/


//Float_t bSign = (fAOD->GetMagneticField() > 0) ? 1 : -1;

/*fArrayMC = (TClonesArray*)fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if (!fArrayMC) {
   // AliFatal("Error: MC particles branch not found!\n");
return;
  }
*/
/*AliAODMCHeader *mcHdr=NULL;
  mcHdr=(AliAODMCHeader*)fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());  
  if(!mcHdr) {
    //printf("MC header branch not found!\n");
    return;
	}
*/
// Fill the TObjArray which holds Global tracks.
	FillGlobalTracksArray();

// Create object arrays for associateds.
	TObjArray *associateds	= new TObjArray();


Int_t TriggerIndx=-99999;
Double_t TriggerPtMin=4.00;
Double_t TriggerPtMax=8.00;
Double_t TriggerPhi=1e-10;
Double_t TriggerEta=1e-10;
Double_t TriggerPt=TriggerPtMin;





  // Track loop to fill a pT spectrum
  for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++)
	 {
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));
    if(!track) AliFatal("Not a standard AOD");
    if (!track) {
      //printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
   	 }
	
	Int_t tracktype=ClassifyTrack(track);
	
	if(tracktype==0) continue;

	if(tracktype==1)
	{

	AliAODTrack *globaltrack=GetGlobalTrack(track->GetID());
	Double_t pT=globaltrack->Pt();
	Double_t Eta=globaltrack->Eta();
	
	
	Double_t dEdx=globaltrack->GetTPCsignal();
	if(dEdx==0) continue;
	

	fHistdEdx->Fill(globaltrack->P(),dEdx);
	fHistPt->Fill(pT);
	

        if(pT>TriggerPtMin && pT<=TriggerPtMax)
	{
	TriggerPt=pT;
	TriggerPhi=globaltrack->Phi();
	TriggerEta=globaltrack->Eta();
	TriggerIndx=track->GetID();
	TriggerPtMin=TriggerPt;
	}
	
	if(pT>0.2 && TMath::Abs(Eta)<0.8)
	associateds->AddLast(globaltrack);
	
	
	}

    
  } //track loop 

	
	Bool_t IsPion=kFALSE;
	Bool_t IsKaonProton=kFALSE;
	//cout<<"***************************************************************************"<<TriggerEta<<endl;
	if(TMath::Abs(TriggerEta)>0.8) return;
	
	if(TriggerIndx!=-99999)
	{
	
	//Float_t bSign = (fAOD->GetMagneticField() > 0) ? 1 : -1;
	
	AliAODTrack *currentAssociateds=0x0;
	
	AliAODTrack *TriggerTrck=GetGlobalTrack(TriggerIndx);
	
	if(!TriggerTrck) return;
	
	
		
	
	
	
	
	Int_t TrgrPtBin=GetTriggerPtBin(TriggerTrck);
	
	if(TrgrPtBin> 7) return;
	
	Bool_t CheckTrackPIDqa=GetTrackStatus(TriggerTrck);
	if(!CheckTrackPIDqa) return;
	
	Double_t nSigma_TPC=(fPIDResponse->NumberOfSigmasTPC(TriggerTrck, (AliPID ::EParticleType)2));

	
	fHistNSigmaTPCPion[TrgrPtBin]->Fill(nSigma_TPC);
	//if(nSigma_TPC >= 4.0 || nSigma_TPC< -7.5) return;
	if(nSigma_TPC > 0.0 && nSigma_TPC< 4.0) IsPion=kTRUE;
	if(nSigma_TPC > -6.0 && nSigma_TPC< -3.0) IsKaonProton=kTRUE;
	//cout<<"NSigma_TPC********************************************************************************"<<nSigma_TPC<<"   "<<TriggerTrck->GetTPCsignalN()<<endl;
	
	EffEtaTrigAll=GetEtaCorrectionFactorTrigAll(TriggerEta) ;
	fTriggerPhiAll->Fill(TriggerPhi);
	fTriggerEtaAll->Fill(TriggerEta,1./(EffEtaTrigAll));
	if(IsPion)
	{
	EffEtaTrigPi=GetEtaCorrectionFactorTrigPion(TriggerEta) ;
	fTriggerPhiPion->Fill(TriggerPhi);
	fTriggerEtaPion->Fill(TriggerEta,1./(EffEtaTrigPi));
	}
	if(IsKaonProton)
	{
	EffEtaTrigPr=GetEtaCorrectionFactorTrigProton(TriggerEta) ;//Eff value of triggers
	fTriggerPhiKaonProton->Fill(TriggerPhi);
	fTriggerEtaKaonProton->Fill(TriggerEta,1./(EffEtaTrigPr));
	
	//cout<<"PROTON TRIGGER***************************************************************"<<"   "<<EffEtaTrigP<<endl;
	}

	

	for(Int_t iasso=0;iasso<associateds->GetEntriesFast();iasso++ )
	{
	currentAssociateds = (AliAODTrack*)(associateds->At(iasso));
	Double_t assoPt=currentAssociateds->Pt();

	Int_t AssoPtBin=(Int_t)(2*assoPt);

	if(assoPt<4.0)
	{
	//Bool_t TwoTrackEff=TwoTrackEfficiency(TriggerTrck,currentAssociateds,0.02,bSign);

	//if(!TwoTrackEff) continue;
	
//Identify Associated**********************************************************************

IdentifyAssociated(TriggerEta,TriggerPhi,IsPion,IsKaonProton,currentAssociateds);

//Identify Associated***********************************************************************

	fAssoPhi->Fill(currentAssociateds->Phi());
	fAssoEta->Fill(currentAssociateds->Eta());

	Double_t delphi= PhiRange(TriggerPhi - currentAssociateds->Phi());
	Double_t deleta=TriggerEta - currentAssociateds->Eta();
	
	Float_t EffEtaAsso=GetEtaCorrectionFactorAsso(currentAssociateds->Eta());//Eff value of associateds 

	Float_t WeightSAll=(1./(EffEtaAsso*EffEtaTrigAll));
	fDihadronCorrelation[AssoPtBin]->Fill(delphi,deleta,WeightSAll);
	if(IsPion)
	{
	Float_t WeightSPi=(1./(EffEtaAsso*EffEtaTrigPi));
	fDihadronCorrelationPion[AssoPtBin]->Fill(delphi,deleta,WeightSPi);

	}
	if(IsKaonProton)
	{
	
	Float_t WeightSPr=(1./(EffEtaAsso*EffEtaTrigPr));
	//cout<<WeightS<<endl;
	fDihadronCorrelationProtonKaon[AssoPtBin]->Fill(delphi,deleta,WeightSPr);//Applied Efficiency Correction
	}

	}

	}
	
	}

//Mixing starts here----------------------------------------------------------------------------------------------------------------
	
	Double_t multiplicity=fAOD->GetNumberOfTracks();
	Double_t MultipORcent=multiplicity;
	AliAODVertex *vtx =fAOD->GetPrimaryVertex();
	Double_t zvertex =vtx->GetZ();
	Double_t poolmin=0;
	Double_t poolmax=150;
	if(TMath:: Abs(zvertex)>=10 || MultipORcent>poolmax || MultipORcent<poolmin) return ;
	AliEventPool* pool =NULL; 
	pool=fPoolMgr->GetEventPool(MultipORcent, zvertex);
	if(pool==NULL)
	{
	//AliInfo(Form("No pool found for multiplicity = %f, zVtx = %f cm", MultipORcent, zvertex));
	 return;
	}
	if (pool->IsReady() || pool->NTracksInPool() > 50000/ 10 || pool->GetCurrentNEvents() >=5)
	{
	const Int_t nMix = pool->GetCurrentNEvents();
	
	for (Int_t jMix=0; jMix<nMix; jMix++)
	{
	TObjArray* bgTracks = NULL;
	bgTracks=pool->GetEvent(jMix);
	if(!bgTracks) continue;
	for(Int_t itrack=0;itrack<bgTracks->GetEntriesFast();++itrack)
	{
	AliPIDCorrParticle *PIDCorrParticle =NULL;
	PIDCorrParticle=(AliPIDCorrParticle*)(bgTracks->At(itrack));
	Double_t assoPt=PIDCorrParticle->Pt();
	Int_t Ptbinbg=(Int_t)(assoPt*2);
	if(assoPt<4.0)
	{
	
	Float_t EffEtaAssoMixed=GetEtaCorrectionFactorAsso(PIDCorrParticle->Eta());
	
	Double_t mixedDelPhi=PhiRange(TriggerPhi-PIDCorrParticle->Phi());
	Double_t mixedDelEta=TriggerEta-PIDCorrParticle->Eta();

	
	if(TriggerIndx!=-99999)
	{
	Float_t WeightMAll=(1./(EffEtaAssoMixed*EffEtaTrigAll));
	fMixedEvent[Ptbinbg]->Fill(mixedDelPhi,mixedDelEta,WeightMAll);
	if(IsPion)
	{ 
	Float_t WeightMPi=(1./(EffEtaAssoMixed*EffEtaTrigPi));
	fMixedPion[Ptbinbg]->Fill(mixedDelPhi,mixedDelEta,WeightMPi);

	}
	if(IsKaonProton)
	{
	Float_t WeightMPr=(1./(EffEtaAssoMixed*EffEtaTrigPr));
	//cout<<WeightM<<endl;
	fMixedProtonKaon[Ptbinbg]->Fill(mixedDelPhi,mixedDelEta,WeightMPr);
	}
	}
	}
	}
	}
	}
	
	
	
	//Update pool
	TObjArray* objArray = NULL;
	objArray = (TObjArray*)AcceptTracksforMixing(fAOD);
	if(!objArray) return;
	if(objArray->GetEntriesFast()>0) {
	pool->UpdatePool(objArray);
	} 

	

  PostData(1, fOutputList);
} 
//_______________________________________________________________________________________

Double_t AliAnalysisTaskPIDCORR::PhiRange(Double_t DPhi)

{
	//
	// Puts the argument in the range [-pi/2,3 pi/2].
	//
	
	if (DPhi < -TMath::Pi()/2) DPhi += 2*TMath::Pi();
	if (DPhi > 3*TMath::Pi()/2) DPhi -= 2*TMath::Pi();	

	return DPhi;
	
} 
 
//________________________________________________________________________
Int_t AliAnalysisTaskPIDCORR:: GetTriggerPtBin(AliAODTrack *track)
{
Int_t Bin;
Double_t Pt_Max=track->Pt();
if(Pt_Max<4.0) Bin=999;
if(Pt_Max>=4.0 && Pt_Max<4.5) Bin=0;
if(Pt_Max>=4.5 && Pt_Max<5.0) Bin=1;
if(Pt_Max>=5.0 && Pt_Max<5.5) Bin=2;
if(Pt_Max>=5.5 && Pt_Max<6.0) Bin=3;
if(Pt_Max>=6.0 && Pt_Max<6.5) Bin=4;
if(Pt_Max>=6.5 && Pt_Max<7.0) Bin=5;
if(Pt_Max>=7.0 && Pt_Max<7.5) Bin=6;
if(Pt_Max>=7.5 && Pt_Max<8.0) Bin=7;
if(Pt_Max>=8.0) Bin=999;
return Bin;


}

//______________________________________________________________________
Bool_t AliAnalysisTaskPIDCORR :: GetTrackStatus(AliAODTrack *track)
{
  if ((track->GetStatus() & AliAODTrack::kTPCin   ) == 0) return kFALSE;
  if ((track->GetStatus() & AliAODTrack::kTPCrefit) == 0) return kFALSE;
  if ((track->GetStatus() & AliAODTrack::kITSrefit) == 0) return kFALSE;
 return kTRUE;
}
//______________________________________________________________________

void AliAnalysisTaskPIDCORR ::DefineEventPool()
{
const Int_t MaxNofEvents=1000;
const Int_t MinNofTracks=50000;
const Int_t NofCentBins=3;
Double_t CentralityBins[NofCentBins+1]={2,20,50,150};
const Int_t NofVrtxBins=10;
Double_t ZvrtxBins[NofVrtxBins+1]={-10,-8,-6,-4,-2,0,2,4,6,8,10};


fPoolMgr = new AliEventPoolManager(MaxNofEvents,MinNofTracks,NofCentBins,CentralityBins,NofVrtxBins,ZvrtxBins);
}
//_______________________________________________________________________
TObjArray *AliAnalysisTaskPIDCORR::AcceptTracksforMixing(AliAODEvent *inputEvent)
{
	Int_t nTracks = inputEvent->GetNumberOfTracks();
	
	

	TObjArray* tracksClone = new TObjArray;
	tracksClone->SetOwner(kTRUE);

	for (Int_t iTrack=0; iTrack<nTracks; ++iTrack)
{
	AliAODTrack* track = dynamic_cast<AliAODTrack*>(inputEvent->GetTrack(iTrack));
	if(!track) AliFatal("Not a standard AOD");
	
	Int_t trackclass=ClassifyTrack(track);
	if (trackclass==1)
	{
	//AliAODTrack *globaltrack=GetGlobalTrack(track->GetID());
	//Double_t MpT=globaltrack->Pt();
	//Double_t MEta=globaltrack->Eta();
	//if(MpT>0.2 && TMath::Abs(MEta)<=0.8)
	tracksClone->Add(new AliPIDCorrParticle(track->Eta(), track->Phi(), track->Pt(),track->Charge()));
	}
}

return tracksClone;
}
//_________________________________________________________________________________________________________________________


Bool_t AliAnalysisTaskPIDCORR::TwoTrackEfficiency(AliAODTrack *trig,AliAODTrack *asso,Float_t ftwoTrackEfficiencyCutValue,Float_t bSign) 
{
 if(!trig || !asso) return kFALSE;
 //if (twoTrackEfficiencyCut)
 //{
	  // the variables & cuthave been developed by the HBT group 
	  // see e.g. https://indico.cern.ch/materialDisplay.py?contribId=36&sessionId=6&materialId=slides&confId=142700

	  Float_t phi1 = trig->Phi();
	  Float_t pt1 = trig->Pt();
	  Float_t charge1 = trig->Charge();
	  Float_t eta1=trig->Eta();
  
	  Float_t phi2 = asso->Phi();
	  Float_t pt2 = asso->Pt();
	  Float_t charge2 = asso->Charge();
	  Float_t eta2=asso->Eta();

	  Float_t deta=eta1-eta2;     
	  // optimization
	  if (TMath::Abs(deta) < ftwoTrackEfficiencyCutValue * 2.5 * 3)
	  {
	    // check first boundaries to see if is worth to loop and find the minimum
	    Float_t dphistar1 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 0.8, bSign);
	    Float_t dphistar2 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 2.5, bSign);
	    
	    const Float_t kLimit = ftwoTrackEfficiencyCutValue * 3;

	    Float_t dphistarminabs = 1e5;
	    //Float_t dphistarmin = 1e5;
	    if (TMath::Abs(dphistar1) < kLimit || TMath::Abs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0)
	    {
	      for (Double_t rad=0.8; rad<2.51; rad+=0.01) 
	      {
		Float_t dphistar = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, rad, bSign);

		Float_t dphistarabs = TMath::Abs(dphistar);
		
		if (dphistarabs < dphistarminabs)
		{
		  //dphistarmin = dphistar;
		  dphistarminabs = dphistarabs;
		}
	      }
	      
	      //fTwoTrackDistancePt[0]->Fill(deta, dphistarmin, TMath::Abs(pt1 - pt2));
	      
	      if (dphistarminabs < ftwoTrackEfficiencyCutValue && TMath::Abs(deta) < ftwoTrackEfficiencyCutValue)
	      {
// 		Printf("Removed track pair %d %d with %f %f %f %f %f %f %f %f %f", i, j, deta, dphistarminabs, phi1, pt1, charge1, phi2, pt2, charge2, bSign);
		return kFALSE;
	      }

    	      //fTwoTrackDistancePt[1]->Fill(deta, dphistarmin, TMath::Abs(pt1 - pt2));
	    }
	  }
	  //}
	  return kTRUE;
}
//_________________________________________________________________________________________________________________________________________

Bool_t AliAnalysisTaskPIDCORR::TwoTrackEfficiencyBg(AliAODTrack *trig,AliPIDCorrParticle *asso,Float_t ftwoTrackEfficiencyCutValue,Float_t bSign) 
{
 if(!trig || !asso) return kFALSE;
 //if (twoTrackEfficiencyCut)
 //{
	  // the variables & cuthave been developed by the HBT group 
	  // see e.g. https://indico.cern.ch/materialDisplay.py?contribId=36&sessionId=6&materialId=slides&confId=142700

	  Float_t phi1 = trig->Phi();
	  Float_t pt1 = trig->Pt();
	  Float_t charge1 = trig->Charge();
	  Float_t eta1=trig->Eta();
  
	  Float_t phi2 = asso->Phi();
	  Float_t pt2 = asso->Pt();
	  Float_t charge2 = asso->Charge();
	  Float_t eta2=asso->Eta();

	  Float_t deta=eta1-eta2;     
	  // optimization
	  if (TMath::Abs(deta) < ftwoTrackEfficiencyCutValue * 2.5 * 3)
	  {
	    // check first boundaries to see if is worth to loop and find the minimum
	    Float_t dphistar1 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 0.8, bSign);
	    Float_t dphistar2 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 2.5, bSign);
	    
	    const Float_t kLimit = ftwoTrackEfficiencyCutValue * 3;

	    Float_t dphistarminabs = 1e5;
	    //Float_t dphistarmin = 1e5;
	    if (TMath::Abs(dphistar1) < kLimit || TMath::Abs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0)
	    {
	      for (Double_t rad=0.8; rad<2.51; rad+=0.01) 
	      {
		Float_t dphistar = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, rad, bSign);

		Float_t dphistarabs = TMath::Abs(dphistar);
		
		if (dphistarabs < dphistarminabs)
		{
		  //dphistarmin = dphistar;
		  dphistarminabs = dphistarabs;
		}
	      }
	      
	      //fTwoTrackDistancePt[0]->Fill(deta, dphistarmin, TMath::Abs(pt1 - pt2));
	      
	      if (dphistarminabs < ftwoTrackEfficiencyCutValue && TMath::Abs(deta) < ftwoTrackEfficiencyCutValue)
	      {
// 		Printf("Removed track pair %d %d with %f %f %f %f %f %f %f %f %f", i, j, deta, dphistarminabs, phi1, pt1, charge1, phi2, pt2, charge2, bSign);
		return kFALSE;
	      }

    	      //fTwoTrackDistancePt[1]->Fill(deta, dphistarmin, TMath::Abs(pt1 - pt2));
	    }
	  }
	  //}
	  return kTRUE;
}

//__________________________________________________________________________________________________________________________________________

Float_t  AliAnalysisTaskPIDCORR::GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign)
{ 
  //
  // calculates dphistar
  //
  
  Float_t dphistar = phi1 - phi2 - charge1 * bSign * TMath::ASin(0.075 * radius / pt1) + charge2 * bSign * TMath::ASin(0.075 * radius / pt2);
  
  static const Double_t kPi = TMath::Pi();
  
  // circularity
//   if (dphistar > 2 * kPi)
//     dphistar -= 2 * kPi;
//   if (dphistar < -2 * kPi)
//     dphistar += 2 * kPi;
  
  if (dphistar > kPi)
    dphistar = kPi * 2 - dphistar;
  if (dphistar < -kPi)
    dphistar = -kPi * 2 - dphistar;
  if (dphistar > kPi) // might look funny but is needed
    dphistar = kPi * 2 - dphistar;
  
  return dphistar;
}
//_______________________________________________________________________
Float_t AliAnalysisTaskPIDCORR::GetEtaCorrectionFactorAsso(Double_t Eta)
{


Int_t Etabin=(16/1.6)*(-0.8-Eta);
Int_t BIN=TMath::Abs(Etabin);

Float_t eff[16]={0.773274,0.783099,0.774943,0.769163,0.759582,0.74767,0.732028,0.724823,0.756575,0.780094,0.788648,0.79306,0.796114,0.801123,0.799477,0.779948};

return eff[BIN];



}
//____________________________________________________________________________________
Float_t AliAnalysisTaskPIDCORR::GetEtaCorrectionFactorTrigAll(Double_t Eta)
{


Int_t Etabin=(16/1.6)*(-0.8-Eta);
Int_t BIN=TMath::Abs(Etabin);



Float_t EffAll[16]={0.705833,0.715773,0.72209,0.718895,0.703701,0.691354,0.679242,0.65723,0.698546,0.72464,0.730565,0.732914,0.733499,0.746637,0.738133,0.720495};

return EffAll[BIN];

}

//_____________________________________________________________________________________________
Float_t AliAnalysisTaskPIDCORR::GetEtaCorrectionFactorTrigPion(Double_t Eta)
{


Int_t Etabin=(16/1.6)*(-0.8-Eta);
Int_t BIN=TMath::Abs(Etabin);


Float_t EffPion[16]={0.36819,0.369563,0.378747,0.378856,0.361426,0.372088,0.355179,0.344626,0.376348,0.387998,0.381553,0.373258,0.397368,0.380326,0.381536,0.369155};

return EffPion[BIN];


}

//______________________________________________________________________________________________________
Float_t AliAnalysisTaskPIDCORR::GetEtaCorrectionFactorTrigProton(Double_t Eta)
{


Int_t Etabin=(16/1.6)*(-0.8-Eta);
Int_t BIN=TMath::Abs(Etabin);

Float_t EffProton[16]={0.698916,0.708791,0.686047,0.663609,0.621142,0.57571,0.556255,0.518182,0.54031,0.58828,0.62885,0.675739,0.686916,0.707395,0.714397,0.693769};
return EffProton[BIN];




}
//_______________________________________________________________________
Bool_t AliAnalysisTaskPIDCORR::CheckTOF(AliVTrack * trk)
{
	//in addition to KTOFout and kTIME we look at the pt
  if(trk->Pt()<0.6)return kFALSE;

//check if the particle has TOF Matching
  UInt_t status;
  status=trk->GetStatus();
  if((status&AliVTrack::kTOFout)==0 || (status&AliVTrack::kTIME)==0)return kFALSE;
  else return kTRUE;
  


}
//________________________________________________________________________
void AliAnalysisTaskPIDCORR::IdentifyAssociated(Double_t Eta_trig,Double_t Phi_trig,Bool_t PION,Bool_t PROTON,AliAODTrack *TRK)
{
Bool_t HasTOFPID=CheckTOF(TRK);
//Bool_t HasTPCPID=GetTrackStatus(TRK);

if(HasTOFPID)
{
Double_t Eta_asso=TRK->Eta();
Double_t Phi_asso=TRK->Phi();

Int_t Ptype=GetTOFPID(TRK);
if(Ptype>2) return;

//cout<<"TOFPID"<<HasTOFPID<<"  "<<"ParticleType"<<Ptype<<endl;
//Int_t Ptype=GetTPCTOFPID(TRK);

if(TMath::Abs(Eta_trig-Eta_asso)<=0.5)
	{

	if(TMath::Abs(Phi_trig-Phi_asso)<=0.5)
	{
	fHistNSAll[Ptype]->Fill(TRK->Pt());
	if(PION)
	fHistNSPion[Ptype]->Fill(TRK->Pt());
	if(PROTON)
	fHistNSProton[Ptype]->Fill(TRK->Pt());
	}

	if((Phi_trig-Phi_asso)>=TMath::Pi()-0.5 && (Phi_trig-Phi_asso)<=TMath::Pi()+0.5)
	{
	fHistASAll[Ptype]->Fill(TRK->Pt());
	if(PION)
	fHistASPion[Ptype]->Fill(TRK->Pt());
	if(PROTON)
	fHistASProton[Ptype]->Fill(TRK->Pt());
	}

	if((Phi_trig-Phi_asso)>=0.5*TMath::Pi()-0.4 && (Phi_trig-Phi_asso)<=0.5*TMath::Pi()+0.2)
	{
	fHistBgAll[Ptype]->Fill(TRK->Pt());
	if(PION)
	fHistBgPion[Ptype]->Fill(TRK->Pt());
	if(PROTON)
	fHistBgProton[Ptype]->Fill(TRK->Pt());
	}


	}

if(TMath::Abs(Eta_trig-Eta_asso)>=1.0)
	{
	fHistBulkAll[Ptype]->Fill(TRK->Pt());
	if(PION)
	fHistBulkPion[Ptype]->Fill(TRK->Pt());
	if(PROTON)
	fHistBulkProton[Ptype]->Fill(TRK->Pt());
	}





}


}

//_______________________________________________________________________
Int_t AliAnalysisTaskPIDCORR::GetTOFPID(AliAODTrack *track)
{

Double_t nsigmaTOFkProton=999.,nsigmaTOFkKaon=999.,nsigmaTOFkPion=999.;

    nsigmaTOFkProton = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
    nsigmaTOFkKaon   = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon); 
    nsigmaTOFkPion   = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);

//cout<<nsigmaTOFkProton<<"   "<<nsigmaTOFkKaon<<"    "<< nsigmaTOFkPion<<endl;

	if(TMath::Abs(nsigmaTOFkPion)<= 2.0) return 0;
	if(TMath::Abs(nsigmaTOFkKaon)<= 2.0) return 1;
	if(TMath::Abs(nsigmaTOFkProton)<= 2.0) return 2;
else return 999;


}
//_______________________________________________________________________
Int_t AliAnalysisTaskPIDCORR::GetTPCTOFPID(AliAODTrack *track)
{

//___TPC
  Double_t nsigmaTPCkProton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
  Double_t nsigmaTPCkKaon   = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon); 
  Double_t nsigmaTPCkPion   = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion); 
//___TOF

Double_t nsigmaTOFkProton=999.,nsigmaTOFkKaon=999.,nsigmaTOFkPion=999.;

    nsigmaTOFkProton = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
    nsigmaTOFkKaon   = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon); 
    nsigmaTOFkPion   = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);

//__TPCTOF

 Double_t nsigmaTPCTOFkProton=999.,nsigmaTPCTOFkKaon=999.,nsigmaTPCTOFkPion=999.;
    nsigmaTPCTOFkProton = TMath::Sqrt((nsigmaTPCkProton*nsigmaTPCkProton+nsigmaTOFkProton*nsigmaTOFkProton));
    nsigmaTPCTOFkKaon   = TMath::Sqrt((nsigmaTPCkKaon*nsigmaTPCkKaon+nsigmaTOFkKaon*nsigmaTOFkKaon));
    nsigmaTPCTOFkPion   = TMath::Sqrt((nsigmaTPCkPion*nsigmaTPCkPion+nsigmaTOFkPion*nsigmaTOFkPion));

    if(TMath::Abs(nsigmaTPCTOFkPion)<= 3.0) return 0;
    if(TMath::Abs(nsigmaTPCTOFkKaon)<= 3.0) return 1;
    if(TMath::Abs(nsigmaTPCTOFkProton)<= 3.0) return 2;

    return -1;
}



//________________________________________________________________________
void AliAnalysisTaskPIDCORR::Terminate(Option_t *) 
{
  // Draw result to screen, or perform fitting, normalizations
  // Called once at the end of the query
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if(!fOutputList) { Printf("ERROR: could not retrieve TList fOutput"); return; }
  
  
}
 
		  
