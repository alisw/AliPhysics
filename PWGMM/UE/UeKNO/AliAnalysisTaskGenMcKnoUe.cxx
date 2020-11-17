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
 * Author: Ahsan Mehmood Khan(ahsan.mehmood.khan@cern.ch)                 * 
 *         Feng Fan (Feng.Fan@cern.ch)				                                  *
 *         Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)                *
 **************************************************************************/

/* AliAnaysisTaskMcKnoUe source code
 * The analysis task produce all the histos needed for MC closure test studies
 * Results include both KNO properties and traditional UE analysis
 */

class TTree;

class AliPPVsMultUtils;
class AliESDtrackCuts;


#include <Riostream.h>
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TList.h"

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVTrack.h"
//#include "AliV0vertexer.h"


#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDVZERO.h"

#include "AliESDInputHandler.h"
#include "AliMultiplicity.h"
#include "AliPPVsMultUtils.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskESDfilter.h"

#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"

#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include <AliHeader.h>

#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultInput.h"
#include "AliMultVariable.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"

#include "AliESDUtils.h"

#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
//#include "AliGenPythiaEventHeader.h"

#include "AliEventCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliPPVsMultUtils.h"
#include <AliESDVertex.h>
#include <AliMultiplicity.h>
#include <TTree.h>
#include <TMath.h>
#include <TRandom.h>
#include <TDirectory.h>
#include <TBits.h>
#include <AliAnalysisFilter.h>

using std::cout;
using std::endl;

#include "AliAnalysisTaskGenMcKnoUe.h"

 
const Char_t * NameOfRegion[3]={"NS","AS","TS"};
const Int_t NchNBins = 100;
Double_t NchBins[NchNBins+1]={-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5,29.5,30.5,31.5,32.5,33.5,34.5,35.5,36.5,37.5,38.5,39.5,40.5,41.5,42.5,43.5,44.5,45.5,46.5,47.5,48.5,49.5,50.5,51.5,52.5,53.5,54.5,55.5,56.5,57.5,58.5,59.5,60.5,61.5,62.5,63.5,64.5,65.5,66.5,67.5,68.5,69.5,70.5,71.5,72.5,73.5,74.5,75.5,76.5,77.5,78.5,79.5,80.5,81.5,82.5,83.5,84.5,85.5,86.5,87.5,88.5,89.5,90.5,91.5,92.5,93.5,94.5,95.5,96.5,97.5,98.5,99.5};

const Int_t pTNBins = 36;
Double_t pTNBins1[pTNBins+1] = {
	0.0,  0.1,  0.15,  0.2,  0.25,  0.3,  0.35,  0.4,  0.45,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5, 5.0, 6.0,    7.0,  8.0,  9.0,  10.0,  12.0,  14.0,  16.0,  18.0,  20.0,  25.0,  30.0,  40.0,  50.0
};

const Int_t pTNBinsL = 24;
Double_t pTNBins1L[pTNBinsL+1] = {
	0.15, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50,
	5.00, 6.00, 7.00, 8.00, 9.00, 10.0, 12.0, 14.0, 16.0, 18.0,
	20.0, 25.0, 30.0, 40.0, 50.0
};



const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;
class AliAnalysisTaskGenMcKnoUe;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskGenMcKnoUe) // classimp: necessary for root

AliAnalysisTaskGenMcKnoUe::AliAnalysisTaskGenMcKnoUe() : AliAnalysisTaskSE(),
    
    fMC(0x0),fMcHandler(0x0),fMCStack(0),fEtaCut(0.8),fPtMin(0.5),fOutputList(0),fGenLeadPhi(0),fGenLeadPt(0),fGenLeadIn(0),  fHistEvt(0x0),hCounter(0),hPtLeadingGenAll(0),hPtLeadingTrue(0)

{
	for(Int_t i=0;i<3;++i){ 
		
		hPtVsPtLeadingTrue[i]=0;
		pNumDenTrueAll[i]=0;
		pSumPtTrueAll[i]=0;
        
        hNumDen[i]=0;
        hSumPt[i]=0;
        
		pNumDenTrue[i]=0;
		pSumPtTrue[i]=0;
	}
	
	// default constructor, don't allocate memory here!  this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskGenMcKnoUe::AliAnalysisTaskGenMcKnoUe(const char* name) : AliAnalysisTaskSE(name),
	
    fMC(0x0),fMcHandler(0x0),fMCStack(0),fEtaCut(0.8),fPtMin(0.5),fOutputList(0),fGenLeadPhi(0),fGenLeadPt(0),fGenLeadIn(0),fHistEvt(0x0),hCounter(0),hPtLeadingGenAll(0),hPtLeadingTrue(0)

{
	for(Int_t i=0;i<3;++i){

		hPtVsPtLeadingTrue[i]=0;
		pNumDenTrueAll[i]=0;
		pSumPtTrueAll[i]=0;
        
        hNumDen[i]=0;
        hSumPt[i]=0;
    
		pNumDenTrue[i]=0;
		pSumPtTrue[i]=0;

	}
	
	// constructor
	DefineInput(0, TChain::Class());    // define the input of the analysis: in this case you take a 'chain' of events
	// this chain is created by the analysis manager, so no need to worry about it, does its work automatically
	DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms

}
//_____________________________________________________________________________
AliAnalysisTaskGenMcKnoUe::~AliAnalysisTaskGenMcKnoUe()
{
	// destructor
	if(fOutputList) {
		delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
		fOutputList = 0x0;
	}

}
//_____________________________________________________________________________
void AliAnalysisTaskGenMcKnoUe::UserCreateOutputObjects()
{

	// create output objects

	OpenFile(1);
	fOutputList = new TList();          // this is a list which will contain all of your histograms
	// at the end of the analysis, the contents of this list are written  to the output file
	fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects and will delete them if requested

    hCounter = new TH1D("hCounter","Counter; sel; Nev",3,0,3);
    fOutputList->Add(hCounter);
    
    fHistEvt = new TH1I("fHistEvt","fHistEvt",2,0,2) ;
    fHistEvt->GetXaxis()->SetBinLabel(1,"All events");
    fHistEvt->GetXaxis()->SetBinLabel(2,"All particles");
    fHistEvt->Sumw2();
    fOutputList->Add(fHistEvt);
    
    hPtLeadingGenAll = new TH1D("hPtLeadingGenAll","gen pTleading before any selection",pTNBinsL,pTNBins1L);
    fOutputList->Add(hPtLeadingGenAll);
    
    for(Int_t i=0;i<3;++i){
    pNumDenTrueAll[i] = new TProfile(Form("pNumDenTrueAll_%s",NameOfRegion[i]),"",pTNBinsL,pTNBins1L);
               fOutputList->Add(pNumDenTrueAll[i]);

               pSumPtTrueAll[i] = new TProfile(Form("pSumPtTrueAll_%s",NameOfRegion[i]),"",pTNBinsL,pTNBins1L);
               fOutputList->Add(pSumPtTrueAll[i]);
               
    }
               hPtLeadingTrue = new TH1D("hPtLeadingTrue","",pTNBinsL,pTNBins1L);
               fOutputList->Add(hPtLeadingTrue);

		// UE analysis
		for(Int_t i=0;i<3;++i){
			
            
            hNumDen[i]= new TH2D(Form("hNumDen_%s",NameOfRegion[i]),"",pTNBinsL,pTNBins1L,NchNBins,NchBins);
            fOutputList->Add(hNumDen[i]);
            hSumPt[i]= new TH2D(Form("hSumPt_%s",NameOfRegion[i]),"",pTNBinsL,pTNBins1L,pTNBins,pTNBins1);
            fOutputList->Add(hSumPt[i]);
            
			hPtVsPtLeadingTrue[i] = new TH2D(Form("hPtVsPtLeadingTrue_%s",NameOfRegion[i]),"",pTNBinsL,pTNBins1L,pTNBins,pTNBins1);
			fOutputList->Add(hPtVsPtLeadingTrue[i]);

			pNumDenTrue[i] = new TProfile(Form("pNumDenTrue_%s",NameOfRegion[i]),"",pTNBinsL,pTNBins1L);
			fOutputList->Add(pNumDenTrue[i]);

			pSumPtTrue[i] = new TProfile(Form("pSumPtTrue_%s",NameOfRegion[i]),"",pTNBinsL,pTNBins1L);
			fOutputList->Add(pSumPtTrue[i]);
		}

   
	//fEventCuts.AddQAplotsToList(fOutputList);
	PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the

}
//_____________________________________________________________________________
void AliAnalysisTaskGenMcKnoUe::UserExec(Option_t *)
{

	// ### Initialize
    
    fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    
    // ### MC handler
	
        if(fMcHandler)
		fMC = fMcHandler->MCEvent();
    
        else { if(fDebug > 1) printf("AliAnalysisTaskGenUeNchTS::Handler() fMcHandler = NULL\n"); return; }
    
    // ### MC event
    
		if(!fMC){
			Printf("%s:%d MCEvent not found in Input Manager",(char*)__FILE__,__LINE__);
			this->Dump();
			return;
		}
		fMCStack = ((AliMCEvent*)fMC)->Stack();
    
    if (!fMCStack) {
        Printf("ERROR: Could not retrieve MC stack \n");
        cout << "Name of the file with pb :" << endl; //fInputHandler->GetTree()->GetCurrentFile()->GetName() << 
        return;
    }
	
	hCounter->Fill(0);

    // ### MC event selection
    
    Bool_t isEventMCSelected = IsMCEventSelected(fMC);
    if( !isEventMCSelected )
       // cout << " nnn 888888 888  88  Name of the file with pb :" << endl;
        return;

    AliHeader* headerMC = fMC->Header();
    
	Bool_t isGoodVtxPosMC = kFALSE;
	
		AliGenEventHeader* genHeader = headerMC->GenEventHeader();
		TArrayF vtxMC(3); // primary vertex  MC 
		vtxMC[0]=9999; vtxMC[1]=9999;  vtxMC[2]=9999; //initialize with dummy
		if (genHeader) {
			genHeader->PrimaryVertex(vtxMC);
		}
		if(TMath::Abs(vtxMC[2])<=10)
			isGoodVtxPosMC = kTRUE;

		// Before trigger selection
	
	GetGenLeadingObject();// leading particle at gen level
	if(isGoodVtxPosMC)
		hPtLeadingGenAll->Fill(fGenLeadPt);

    // GENERATOR LEVEL
    
    vector<Double_t> ue_gen;// 0: nch_near, 1: nch_away, 2: nch_trans, 3: sumpt_near, ...
	GetMeanGenUEObservables(ue_gen);

	for(Int_t i=0;i<3;++i){
		
		if(isGoodVtxPosMC){
			pNumDenTrueAll[i]->Fill(fGenLeadPt,ue_gen[i]);
			pSumPtTrueAll[i]->Fill(fGenLeadPt,ue_gen[i+3]);
            
		}
	}
   
    
	hCounter->Fill(1);

	
			if(isGoodVtxPosMC){
				// UE analysis
				if(fGenLeadPt>=fPtMin){
				GetGenUEObservables();
			}

		}
	
	ue_gen.clear();
	
	PostData(1, fOutputList); // stream the result of this event to the output manager which will write it to a file

}

//______________________________________________________________________________
void AliAnalysisTaskGenMcKnoUe::Terminate(Option_t *)
{

}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskGenMcKnoUe::IsMCEventSelected(TObject* obj){

    Bool_t isSelected = kTRUE;

    AliMCEvent *event = 0x0;
    event = dynamic_cast<AliMCEvent*>(obj);
    if( !event )
        isSelected = kFALSE;

    if( isSelected )
        FillHisto("fHistEvt",0.5);

    return isSelected;
    
}
//---------------------------------------
inline void AliAnalysisTaskGenMcKnoUe::FillHisto(const char* objkey, Double_t x)
{
    TH1* hTmp = 0;
    hTmp = static_cast<TH1*>(fOutputList->FindObject(objkey));
    if(!hTmp){
        AliError(Form("Cannot find histogram: %s",objkey)) ;
        return;
    }
    hTmp->Fill(x);
}
//------------------------------------------

void AliAnalysisTaskGenMcKnoUe::GetMeanGenUEObservables(vector<Double_t> &genArray){

    genArray.clear();
    

    Int_t nch_top[3];
    Double_t sumpt_top[3];
    for(Int_t i=0;i<3;++i){
        nch_top[i]=0;
        sumpt_top[i]=0;
    }


        for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {

            if(i==fGenLeadIn)
                continue;

            AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
            if (!particle) continue;
            
            if (!fMC->IsPhysicalPrimary(i)) continue;
            if (particle->Charge() == 0) continue;
            if ( TMath::Abs(particle->Eta()) > fEtaCut )continue;
            if( particle->Pt() < fPtMin)continue;

            Double_t DPhi = DeltaPhi(particle->Phi(), fGenLeadPhi);

            // definition of the topological regions
            if(TMath::Abs(DPhi)<pi/3.0){// near side
                nch_top[0]++; sumpt_top[0]+=particle->Pt();
            }
            else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
                nch_top[1]++; sumpt_top[1]+=particle->Pt();
            }
            else{// transverse side
                nch_top[2]++; sumpt_top[2]+=particle->Pt();
            }

        }
    
    
        for(Int_t i=0;i<3;++i)
            genArray.push_back(1.0*nch_top[i]);
        for(Int_t i=0;i<3;++i)
            genArray.push_back(sumpt_top[i]);
    
}
//-------------------------------------------------
void AliAnalysisTaskGenMcKnoUe::GetGenLeadingObject() {

	Double_t flPt = 0;// leading pT
	Double_t flPhi = 0;
	Int_t flIndex = 0;

	
		for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {

			AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
			if (!particle) continue;
            
            FillHisto("fHistEvt",1.5);
            
            
			if (!fMC->IsPhysicalPrimary(i)) continue;  //  Particles produced including products of strong and electromagnetic decays but excluding feed-down from weak decays of strange particles like Ks,Lambda etc)
			if (particle->Charge() == 0) continue;
			if ( TMath::Abs(particle->Eta()) > fEtaCut )continue;
			if( particle->Pt() < fPtMin)continue;

			if (flPt<particle->Pt()){
				flPt = particle->Pt();
				flPhi = particle->Phi();
				flIndex = i;
                
			}
		}

		fGenLeadPhi = flPhi;
		fGenLeadPt  = flPt;
		fGenLeadIn  = flIndex;
        cout <<"Lead pt = %f  \n"<< fGenLeadPt;
}
//----------------------
void AliAnalysisTaskGenMcKnoUe::GetGenUEObservables(){

    Int_t nch_top[3];
    Double_t sumpt_top[3];
    for(Int_t i=0;i<3;++i){
        nch_top[i]=0;
        sumpt_top[i]=0;
    }


    for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {

        if(i==fGenLeadIn)
            continue;

        AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
        if (!particle) continue;

        if (!fMC->IsPhysicalPrimary(i)) continue;
        if (particle->Charge() == 0) continue;
        if ( TMath::Abs(particle->Eta()) > fEtaCut )continue;
        if( particle->Pt() < fPtMin)continue;

        Double_t DPhi = DeltaPhi(particle->Phi(), fGenLeadPhi);

        // definition of the topological regions
        if(TMath::Abs(DPhi)<pi/3.0){// near side
            nch_top[0]++; sumpt_top[0]+=particle->Pt();
            hPtVsPtLeadingTrue[0]->Fill(fGenLeadPt,particle->Pt());
        }
        else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
            nch_top[1]++; sumpt_top[1]+=particle->Pt();
            hPtVsPtLeadingTrue[1]->Fill(fGenLeadPt,particle->Pt());
        }
        else{// transverse side
            nch_top[2]++; sumpt_top[2]+=particle->Pt();
            hPtVsPtLeadingTrue[2]->Fill(fGenLeadPt,particle->Pt());
        }



    }

    for(Int_t i=0;i<3;++i){
        
        pNumDenTrue[i]->Fill(fGenLeadPt,nch_top[i]);
        pSumPtTrue[i]->Fill(fGenLeadPt,sumpt_top[i]);
        hNumDen[i]->Fill(fGenLeadPt,nch_top[i]);
        hSumPt[i]->Fill(fGenLeadPt,sumpt_top[i]);
    }
    hPtLeadingTrue->Fill(fGenLeadPt);

}

//____________________________________________________________

Double_t AliAnalysisTaskGenMcKnoUe::DeltaPhi(Double_t phia, Double_t phib,
		Double_t rangeMin, Double_t rangeMax)
{
	Double_t dphi = -999;
	Double_t pi = TMath::Pi();

	if (phia < 0)         phia += 2*pi;
	else if (phia > 2*pi) phia -= 2*pi;
	if (phib < 0)         phib += 2*pi;
	else if (phib > 2*pi) phib -= 2*pi;
	dphi = phib - phia;
	if (dphi < rangeMin)      dphi += 2*pi;
	else if (dphi > rangeMax) dphi -= 2*pi;

	return dphi;
}

