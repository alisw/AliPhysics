
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *				       					  *
 * Authors: Svein Lindal, Daniel Lohner					  *
 * Version 1.0								  *
 *									  *
 * Permission to use, copy, modify and distribute this software and its	  *
 * documentation strictly for non-commercial purposes is hereby granted	  *
 * without fee, provided that the above copyright notice appears in all	  *
 * copies and that both the copyright notice and this permission notice	  *
 * appear in the supporting documentation. The authors make no claims	  *
 * about the suitability of this software for any purpose. It is	  *
 * provided "as is" without express or implied warranty.		  *
 **************************************************************************/

////////////////////////////////////////////////
//--------------------------------------------- 
// Class handling all kinds of selection cuts for
// Gamma Conversion analysis
//---------------------------------------------
////////////////////////////////////////////////


#include "AliDalitzElectronCuts.h"
#include "AliAODConversionPhoton.h"
#include "AliKFVertex.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "TH1.h"
#include "TH2.h"
#include "AliMCEvent.h"
#include "TObjString.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "TList.h"
class iostream;

using namespace std;

ClassImp(AliDalitzElectronCuts)


const char* AliDalitzElectronCuts::fgkCutNames[AliDalitzElectronCuts::kNCuts] = {
"MaxChi2TPCConstrainedGlobal",
"ededxSigmaITSCut",
"ededxSigmaTPCCut",
"pidedxSigmaTPCCut",
"piMinMomdedxSigmaTPCCut",
"piMaxMomdedxSigmaTPCCut",
"LowPRejectionSigmaCut",
"kTOFelectronPID",
"clsITSCut",
"clsTPCCut",
"EtaCut",
"PsiPair",
"RejectSharedElecGamma",
"MaxChi2PerClusterTPC",
"MaxChi2PerClusterITS",
"PtCut",
"DCAcut",
"MassCut",
"Weights",
"VPhotonMCPSmearing"
};

//________________________________________________________________________
AliDalitzElectronCuts::AliDalitzElectronCuts(const char *name,const char *title) : AliAnalysisCuts(name,title),
    fHistograms(NULL),
    fPIDResponse(NULL),
    fesdTrackCuts(NULL),
    fEtaCut(0.9),
    fDoEtaCut(kFALSE),
    fPtMinCut(0.0),
    fPtMaxCut(9999),
    fRadiusCut(1000.0),
    fPsiPairCut(0.45),
    fDeltaPhiCutMin(0.),
    fDeltaPhiCutMax(0.12),
    fMinClsTPC(0), // minimum clusters in the TPC
    fMinClsTPCToF(0), // minimum clusters to findable clusters
    fDodEdxSigmaITSCut(kFALSE),
    fDodEdxSigmaTPCCut(kTRUE),
    fDoTOFsigmaCut(kFALSE), // RRnewTOF
    fDoRejectSharedElecGamma(kFALSE),
    fDoPsiPairCut(kFALSE),
    fPIDnSigmaAboveElectronLineITS(100),
    fPIDnSigmaBelowElectronLineITS(-100),
    fPIDnSigmaAboveElectronLineTPC(100),
    fPIDnSigmaBelowElectronLineTPC(-100),
    fPIDnSigmaAbovePionLineTPC(0),
    fPIDnSigmaAbovePionLineTPCHighPt(-100),
    fTofPIDnSigmaAboveElectronLine(100), // RRnewTOF
    fTofPIDnSigmaBelowElectronLine(-100), // RRnewTOF
    fPIDMinPnSigmaAbovePionLineTPC(0),
    fPIDMaxPnSigmaAbovePionLineTPC(0),
    fDoKaonRejectionLowP(kFALSE),
    fDoProtonRejectionLowP(kFALSE),
    fDoPionRejectionLowP(kFALSE),
    fPIDnSigmaAtLowPAroundKaonLine(0),
    fPIDnSigmaAtLowPAroundProtonLine(0),
    fPIDnSigmaAtLowPAroundPionLine(0),
    fPIDMinPKaonRejectionLowP(1.5),
    fPIDMinPProtonRejectionLowP(2.0),
    fPIDMinPPionRejectionLowP(0.5),
    fUseCorrectedTPCClsInfo(kFALSE),
    fUseCrossedRows(kFALSE),
    fUseTOFpid(kFALSE),
    fRequireTOF(kFALSE),
    fDoMassCut(kFALSE),
    fDoMassMinCut(kFALSE),
    fMassCutLowPt(999.),
    fMassCutHighPt(999.),
    fMassCutPtMin(-100.0),
    fMassMinCut(-999.),
    fDoWeights(kFALSE),
    fUseVPhotonMCPSmearing(kFALSE),
    fUseElectronMCPSmearing(kFALSE),
    fCutString(NULL),
    fCutStringRead(""),
    hCutIndex(NULL),
    hdEdxCuts(NULL),
    hITSdEdxbefore(NULL),
    hITSdEdxafter(NULL),
    hTPCdEdxbefore(NULL),
    hTPCdEdxafter(NULL),
    hTPCdEdxSignalbefore(NULL),
    hTPCdEdxSignalafter(NULL),
    hTOFbefore(NULL),
    hTOFafter(NULL),
    hTrackDCAxyPtbefore(NULL),
    hTrackDCAxyPtafter(NULL),
    hTrackDCAzPtbefore(NULL),
    hTrackDCAzPtafter(NULL),
    hTrackNFindClsPtTPCbefore(NULL),
    hTrackNFindClsPtTPCafter(NULL),
    hTrackPosEtabeforeDedx(NULL),
    hTrackNegEtabeforeDedx(NULL),
    hTrackPosEtaafterDedx(NULL),
    hTrackNegEtaafterDedx(NULL)
   {
    InitPIDResponse();
    for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=0;}
    fCutString=new TObjString((GetCutNumber()).Data());

   //fesdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
   // Using standard function for setting Cuts
    Bool_t selectPrimaries=kFALSE;
    fesdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
}

//________________________________________________________________________
AliDalitzElectronCuts::~AliDalitzElectronCuts() {
    // Destructor
  //Deleting fHistograms leads to seg fault it it's added to output collection of a task
  // if(fHistograms)
  // 	delete fHistograms;
  // fHistograms = NULL;

   if(fCutString != NULL){
      delete fCutString;
      fCutString = NULL;
   }


}

//________________________________________________________________________
void AliDalitzElectronCuts::InitCutHistograms(TString name, Bool_t preCut,TString cutNumber){

    // Initialize Cut Histograms for QA (only initialized and filled if function is called)

     TH1::AddDirectory(kFALSE);

     TString cutName = "";
    
     if( cutNumber==""){
         cutName = GetCutNumber().Data();
     }
     else {
            cutName = cutNumber.Data();
     } 

    if(fHistograms != NULL){
	delete fHistograms;
	fHistograms=NULL;
    }
    if(fHistograms==NULL){
	fHistograms=new TList();
	if(name=="")fHistograms->SetName(Form("ElectronCuts_%s",cutName.Data()));
	else fHistograms->SetName(Form("%s_%s",name.Data(),cutName.Data()));
    }
    
    
    Int_t kDedxSignalbins = 200;
    
     const Int_t kDCABins=62;
    
     Double_t binsDCADummy[63]={-3.0,-2.7,-2.4,-2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,-0.25,-0.2,-0.19,-0.18,-0.17,-0.16,-0.15,-0.14,-0.13,-0.12,-0.11,-0.10,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.25,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3.0};

     const Int_t kPtBins=110;
     Double_t binsPtDummy[kPtBins+1];
     const Int_t kPBins = 109;
     Double_t binsPDummy[kPBins+1];
     binsPtDummy[0]=0.0;
     binsPDummy[0]=0.05;
     
        for(Int_t i=1;i<kPtBins+1;i++)
        {
                if(binsPtDummy[i-1]+0.05<1.01)
                        binsPtDummy[i]=binsPtDummy[i-1]+0.05;
                else
                        binsPtDummy[i]=binsPtDummy[i-1]+0.1;
		
        }
        for(Int_t i=1; i <kPBins+1;i++){
		  
		  if( binsPDummy[i-1]+0.05<1.01)
		        binsPDummy[i] = binsPDummy[i-1]+0.05;
		  else
			binsPDummy[i] = binsPDummy[i-1]+0.1;
		
	}
       


    hCutIndex=new TH1F(Form("IsElectronSelected %s",cutName.Data()),"IsElectronSelected",10,-0.5,9.5);
    hCutIndex->GetXaxis()->SetBinLabel(kElectronIn+1,"in");
    hCutIndex->GetXaxis()->SetBinLabel(kNoTracks+1,"no tracks");
    hCutIndex->GetXaxis()->SetBinLabel(kTrackCuts+1,"Track cuts");
    hCutIndex->GetXaxis()->SetBinLabel(kdEdxCuts+1,"dEdx");
    hCutIndex->GetXaxis()->SetBinLabel(kElectronOut+1,"out");
    fHistograms->Add(hCutIndex);



    // dEdx Cuts
    hdEdxCuts=new TH1F(Form("dEdxCuts %s",cutName.Data()),"dEdxCuts",10,-0.5,9.5);
    hdEdxCuts->GetXaxis()->SetBinLabel(1,"in");
    hdEdxCuts->GetXaxis()->SetBinLabel(2,"ITSelectron");
    hdEdxCuts->GetXaxis()->SetBinLabel(3,"TPCelectron");
    hdEdxCuts->GetXaxis()->SetBinLabel(4,"TPCpion");
    hdEdxCuts->GetXaxis()->SetBinLabel(5,"TPCpionhighp");
    hdEdxCuts->GetXaxis()->SetBinLabel(6,"TPCkaonlowprej");
    hdEdxCuts->GetXaxis()->SetBinLabel(7,"TPCprotonlowprej");
    hdEdxCuts->GetXaxis()->SetBinLabel(8,"TPCpionlowprej");
    hdEdxCuts->GetXaxis()->SetBinLabel(9,"TOFelectron");
    hdEdxCuts->GetXaxis()->SetBinLabel(10,"out");
    fHistograms->Add(hdEdxCuts);
    


    TAxis *AxisBeforeITS  = NULL;
    TAxis *AxisBeforedEdx = NULL;
    TAxis *AxisBeforeTOF  = NULL;
    TAxis *AxisBeforedEdxSignal = NULL;

    if(preCut){


       hITSdEdxbefore=new TH2F(Form("Electron_ITS_before %s",cutName.Data()),"ITS dEdx electron before" ,kPBins,binsPDummy,200,-10,10);
       fHistograms->Add(hITSdEdxbefore);
       AxisBeforeITS = hITSdEdxbefore->GetXaxis();

       hTPCdEdxbefore=new TH2F(Form("Electron_dEdx_before %s",cutName.Data()),"dEdx electron before" ,kPBins,binsPDummy,200,-10,10);
       fHistograms->Add(hTPCdEdxbefore);
       AxisBeforedEdx = hTPCdEdxbefore->GetXaxis();

       hTPCdEdxSignalbefore=new TH2F(Form("Electron_dEdxSignal_before %s",cutName.Data()),"dEdx electron signal before" ,kPBins,binsPDummy,kDedxSignalbins,0.0,200);
       fHistograms->Add(hTPCdEdxSignalbefore);
       AxisBeforedEdxSignal = hTPCdEdxSignalbefore->GetXaxis();

       hTOFbefore=new TH2F(Form("Electron_TOF_before %s",cutName.Data()),"TOF electron before" ,kPBins,binsPDummy,200,-10,10);
       fHistograms->Add(hTOFbefore);
       AxisBeforeTOF = hTOFbefore->GetXaxis();
       
       hTrackDCAxyPtbefore = new TH2F(Form("hTrack_DCAxy_Pt_before %s",cutName.Data()),"DCAxy Vs Pt of tracks before",kDCABins,binsDCADummy,kPtBins,binsPtDummy);
       fHistograms->Add(hTrackDCAxyPtbefore); 	
       
       hTrackDCAzPtbefore  = new TH2F(Form("hTrack_DCAz_Pt_before %s",cutName.Data()), "DCAz  Vs Pt of tracks before",kDCABins,binsDCADummy,kPtBins,binsPtDummy);
       fHistograms->Add(hTrackDCAzPtbefore); 
       
       hTrackNFindClsPtTPCbefore = new TH2F(Form("hTrack_NFindCls_Pt_TPC_before %s",cutName.Data()),"Track: N Findable Cls TPC Vs Pt before",60,0,1.5,kPtBins,binsPtDummy);
       fHistograms->Add(hTrackNFindClsPtTPCbefore); 
	
       

    }


    hITSdEdxafter=new TH2F(Form("Electron_ITS_after %s",cutName.Data()),"ITS dEdx electron after" ,kPBins,binsPDummy,200, -10,10);
    fHistograms->Add(hITSdEdxafter);

    hTPCdEdxafter=new TH2F(Form("Electron_dEdx_after %s",cutName.Data()),"dEdx electron after" ,kPBins,binsPDummy,200, -10,10);
    fHistograms->Add(hTPCdEdxafter);

    hTPCdEdxSignalafter=new TH2F(Form("Electron_dEdxSignal_after %s",cutName.Data()),"dEdx electron signal after" ,kPBins,binsPDummy,kDedxSignalbins,0.0,200);
    fHistograms->Add(hTPCdEdxSignalafter);

    hTOFafter=new TH2F(Form("Electron_TOF_after %s",cutName.Data()),"TOF electron after" ,kPBins,binsPDummy,200,-6,10);
    fHistograms->Add(hTOFafter);
      
    hTrackDCAxyPtafter  = new TH2F(Form("hTrack_DCAxy_Pt_after %s",cutName.Data()),"DCAxy Vs Pt of tracks after",kDCABins,binsDCADummy,kPtBins,binsPtDummy);
    fHistograms->Add(hTrackDCAxyPtafter); 
    
    hTrackDCAzPtafter  = new TH2F(Form("hTrack_DCAz_Pt_after %s",cutName.Data()), "DCAz Vs Pt of tracks  after",kDCABins,binsDCADummy,kPtBins,binsPtDummy);
    fHistograms->Add(hTrackDCAzPtafter); 
    
    hTrackNFindClsPtTPCafter = new TH2F(Form("hTrack_NFindCls_Pt_TPC_after %s",cutName.Data()),"Track: N Findable Cls TPC Vs Pt after",60,0,1.5,kPtBins,binsPtDummy);
    fHistograms->Add(hTrackNFindClsPtTPCafter); 
    
    hTrackPosEtabeforeDedx = new TH1F(Form("hTrack_Pos_Eta_before_Dedx %s",cutName.Data()),"hTrack_Pos_Eta_before_Dedx",600,-1.5,1.5);
    fHistograms->Add(hTrackPosEtabeforeDedx);
    
    hTrackNegEtabeforeDedx = new TH1F(Form("hTrack_Neg_Eta_before_Dedx %s",cutName.Data()),"hTrack_Neg_Eta_before_Dedx",600,-1.5,1.5);
    fHistograms->Add(hTrackNegEtabeforeDedx);
    
    hTrackPosEtaafterDedx  = new TH1F(Form("hTrack_Pos_Eta_after_Dedx %s",cutName.Data()),"hTrack_Pos_Eta_after_Dedx",600,-1.5,1.5);
    fHistograms->Add(hTrackPosEtaafterDedx);
    
    hTrackNegEtaafterDedx  = new TH1F(Form("hTrack_Neg_Eta_afterDedx %s",cutName.Data()),"hTrack_Neg_Eta_after_Dedx",600,-1.5,1.5);
    fHistograms->Add(hTrackNegEtaafterDedx);
    
    

    TAxis *AxisAfter = hTPCdEdxafter->GetXaxis(); 
    Int_t bins = AxisAfter->GetNbins();
    Double_t from = AxisAfter->GetXmin();
    Double_t to = AxisAfter->GetXmax();
    Double_t *newBins = new Double_t[bins+1];
    newBins[0] = from;
    Double_t factor = TMath::Power(to/from, 1./bins);
    for(Int_t i=1; i<=bins; ++i) newBins[i] = factor * newBins[i-1];
    AxisAfter->Set(bins, newBins);
    AxisAfter = hTOFafter->GetXaxis(); 
    AxisAfter->Set(bins, newBins);
    AxisAfter = hITSdEdxafter->GetXaxis();
    AxisAfter->Set(bins,newBins); 
    AxisAfter = hTPCdEdxSignalafter->GetXaxis();
    AxisAfter->Set(bins,newBins);
    
    if(preCut){
       AxisBeforeITS->Set(bins, newBins);
       AxisBeforedEdx->Set(bins, newBins);
       AxisBeforedEdxSignal->Set(bins,newBins);
       AxisBeforeTOF->Set(bins, newBins);
       
    }
    delete [] newBins;

    TH1::AddDirectory(kTRUE);        

    // Event Cuts and Info
}


//________________________________________________________________________
Bool_t AliDalitzElectronCuts::InitPIDResponse(){

// Set Pointer to AliPIDResponse

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();

  if(man) {

    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
    if(fPIDResponse)return kTRUE;

  }

  return kFALSE;
}
///________________________________________________________________________
Bool_t AliDalitzElectronCuts::ElectronIsSelectedMC(Int_t labelParticle,AliMCEvent *mcEvent)
{   
        if( labelParticle < 0 || labelParticle >= mcEvent->GetNumberOfTracks() ) return kFALSE;
        //if( mcEvent->IsPhysicalPrimary(labelParticle) == kFALSE ) return kFALSE; //Ask Ana

        TParticle* particle = mcEvent->Particle(labelParticle);

        if( TMath::Abs( particle->GetPdgCode() ) != 11 )  return kFALSE;
        
        if( fDoEtaCut ){
	  if( particle->Eta() > fEtaCut  || particle->Eta() < -fEtaCut  )
	  return kFALSE;
	}
        

return kTRUE;
}


///________________________________________________________________________
Bool_t AliDalitzElectronCuts::ElectronIsSelected(AliESDtrack* lTrack)
{
    //Selection of Reconstructed electrons
    
    
    Float_t b[2];
    Float_t bCov[3];
    lTrack->GetImpactParameters(b,bCov);
   
    if (bCov[0]<=0 || bCov[2]<=0) {
	AliDebug(1, "Estimated b resolution lower or equal zero!");
	bCov[0]=0; bCov[2]=0;
    }
    
    

    Float_t dcaToVertexXY = b[0];
    Float_t dcaToVertexZ  = b[1];
    Double_t clsToF = GetNFindableClustersTPC(lTrack);
   
   if( hTrackDCAxyPtbefore) hTrackDCAxyPtbefore->Fill(dcaToVertexXY,lTrack->Pt());
   if( hTrackDCAzPtbefore ) hTrackDCAzPtbefore->Fill( dcaToVertexZ, lTrack->Pt());
   if( hTrackNFindClsPtTPCbefore ) hTrackNFindClsPtTPCbefore->Fill( clsToF, lTrack->Pt());
   
   

    if(hCutIndex)hCutIndex->Fill(kElectronIn);

    if (lTrack == NULL){
      if(hCutIndex)hCutIndex->Fill(kNoTracks);
         return kFALSE;  
    }   
       
    if ( ! lTrack->GetConstrainedParam() ){
        return kFALSE;
    }
    AliVTrack * track = dynamic_cast<AliVTrack*>(lTrack);


    // Track Cuts
    if( !TrackIsSelected(lTrack) ){
         if(hCutIndex)hCutIndex->Fill(kTrackCuts);
         return kFALSE;
    }

    if( lTrack->GetSign() > 0.0 ){
       
     if (hTrackPosEtabeforeDedx) hTrackPosEtabeforeDedx->Fill(lTrack->Eta());
      
    } else{
      
      if(hTrackNegEtabeforeDedx) hTrackNegEtabeforeDedx->Fill(lTrack->Eta());
      
    }
    
    
    // dEdx Cuts
    if( ! dEdxCuts( track ) ) {
         if(hCutIndex)hCutIndex->Fill(kdEdxCuts);
         return kFALSE;

    }
    
    if( lTrack->GetSign() > 0.0 ){
       
      if( hTrackPosEtaafterDedx) hTrackPosEtaafterDedx->Fill(lTrack->Eta());
      
    } else{
      
      if( hTrackNegEtaafterDedx) hTrackNegEtaafterDedx->Fill(lTrack->Eta());
      
    }
    
    

    //Electron passed the cuts
    if(hCutIndex)hCutIndex->Fill(kElectronOut);
    
    if( hTrackDCAxyPtafter) 	   hTrackDCAxyPtafter->Fill(dcaToVertexXY,lTrack->Pt());
    if( hTrackDCAzPtafter ) 	   hTrackDCAzPtafter->Fill(dcaToVertexZ,lTrack->Pt());
    if( hTrackNFindClsPtTPCafter ) hTrackNFindClsPtTPCafter->Fill( clsToF, lTrack->Pt());


    return kTRUE;
}

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::TrackIsSelected(AliESDtrack* lTrack) {
    // Track Selection for Photon Reconstruction
    
         
    Double_t clsToF = GetNFindableClustersTPC(lTrack);
    

    if( ! fesdTrackCuts->AcceptTrack(lTrack) ){

        return kFALSE;
    }
    
    if( fDoEtaCut ) {
      if(  lTrack->Eta() > fEtaCut  || lTrack->Eta() < -fEtaCut ) {
        return kFALSE;
      }
   }
   
   
   if( lTrack->Pt() < fPtMinCut || lTrack->Pt() > fPtMaxCut ) {
     
	return kFALSE;
	
   }

   

    if( clsToF < fMinClsTPCToF){
    return kFALSE;
    }

    

   return kTRUE;
}
///________________________________________________________________________
Bool_t AliDalitzElectronCuts::dEdxCuts(AliVTrack *fCurrentTrack){

    // Electron Identification Cuts for Photon reconstruction

    if(!fPIDResponse){  InitPIDResponse();  }// Try to reinitialize PID Response
    if(!fPIDResponse){  AliError("No PID Response"); return kFALSE;}// if still missing fatal error



    //cout<<"dEdxCuts: //////////////////////////////////////////////////////////////////////////"<<endl;



    Int_t cutIndex=0;

    if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
    if(hITSdEdxbefore)hITSdEdxbefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasITS(fCurrentTrack, AliPID::kElectron));
    if(hTPCdEdxbefore)hTPCdEdxbefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTPC(fCurrentTrack, AliPID::kElectron));
    if(hTPCdEdxSignalbefore)hTPCdEdxSignalbefore->Fill(fCurrentTrack->P(),TMath::Abs(fCurrentTrack->GetTPCsignal()));


    cutIndex++;


		if( fDodEdxSigmaITSCut == kTRUE ){


			if( fPIDResponse->NumberOfSigmasITS(fCurrentTrack,AliPID::kElectron)<fPIDnSigmaBelowElectronLineITS ||
					fPIDResponse->NumberOfSigmasITS(fCurrentTrack,AliPID::kElectron)> fPIDnSigmaAboveElectronLineITS ){
				
				if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
				return kFALSE;
      }
			
		}
		
		if(hITSdEdxafter)hITSdEdxafter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasITS(fCurrentTrack, AliPID::kElectron));
		
		
		cutIndex++;
		
		
		if(fDodEdxSigmaTPCCut == kTRUE){
			
			
      // TPC Electron Line
      if( fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)<fPIDnSigmaBelowElectronLineTPC ||
					fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)>fPIDnSigmaAboveElectronLineTPC){
				
				if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
				return kFALSE;
      }
      cutIndex++;
			
      // TPC Pion Line
			if( fCurrentTrack->P()>fPIDMinPnSigmaAbovePionLineTPC && fCurrentTrack->P()<fPIDMaxPnSigmaAbovePionLineTPC ){
				if(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLineTPC     &&
					 fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLineTPC &&
					 fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)<fPIDnSigmaAbovePionLineTPC){
					
					if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
					return kFALSE;
				}
			}
			cutIndex++;
			
			// High Pt Pion rej
			if( fCurrentTrack->P()>fPIDMaxPnSigmaAbovePionLineTPC ){
				if(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLineTPC &&
					 fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLineTPC&&
					 fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)<fPIDnSigmaAbovePionLineTPCHighPt){
					
					if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
					return kFALSE;
				}
			}
			
			cutIndex++;
		}

		else{ cutIndex+=3; }


		if(   fDoKaonRejectionLowP == kTRUE   ){

			if( fCurrentTrack->P() < fPIDMinPKaonRejectionLowP ){
				
				if( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kKaon))<fPIDnSigmaAtLowPAroundKaonLine){
					
					if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
					
					return kFALSE;
				}
			}
		}
		cutIndex++;
		
		if(   fDoProtonRejectionLowP == kTRUE    ){
			
			if( fCurrentTrack->P()  < fPIDMinPProtonRejectionLowP ){
				if( TMath::Abs(   fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kProton))<fPIDnSigmaAtLowPAroundProtonLine){
					
					if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
					return kFALSE;
				}
			}
		}
		cutIndex++;
		
		if(fDoPionRejectionLowP == kTRUE){
			if( fCurrentTrack->P() < fPIDMinPPionRejectionLowP ){
				if( TMath::Abs( fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)) < fPIDnSigmaAtLowPAroundPionLine ){
					
					if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
					return kFALSE;
				}
			}
		}
		cutIndex++;
		
		
		if( ( fCurrentTrack->GetStatus() & AliESDtrack::kTOFpid ) && ( !( fCurrentTrack->GetStatus() & AliESDtrack::kTOFmismatch) ) ){
			if(hTOFbefore) hTOFbefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron));
			if(fUseTOFpid){
        if(fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron)>fTofPIDnSigmaAboveElectronLine ||
           fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron)<fTofPIDnSigmaBelowElectronLine ){
					if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
					return kFALSE;
        }
			}
			if(hTOFafter)hTOFafter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron));
		}
		else if ( fRequireTOF == kTRUE ) {
			
			if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
			return kFALSE;
		}
		cutIndex++;
		
		if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
		if(hTPCdEdxafter)hTPCdEdxafter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTPC(fCurrentTrack, AliPID::kElectron));
		if(hTPCdEdxSignalafter)hTPCdEdxSignalafter->Fill(fCurrentTrack->P(),TMath::Abs(fCurrentTrack->GetTPCsignal()));
		
		return kTRUE;
}
///________________________________________________________________________


AliVTrack *AliDalitzElectronCuts::GetTrack(AliVEvent * event, Int_t label){
    //Returns pointer to the track with given ESD label
    //(Important for AOD implementation, since Track array in AOD data is different
    //from ESD array, but ESD tracklabels are stored in AOD Tracks)

  AliESDEvent * esdEvent = dynamic_cast<AliESDEvent*>(event);
  if(esdEvent) {
  	if(label > event->GetNumberOfTracks() ) return NULL;
  	AliESDtrack * track = esdEvent->GetTrack(label);
  	return track;
	
  } else { 
	for(Int_t ii=0; ii<event->GetNumberOfTracks(); ii++) {
	  AliVTrack * track = dynamic_cast<AliVTrack*>(event->GetTrack(ii));
	  
	  if(track) { 
		if(track->GetID() == label) {
		  return track;
		}
	  }
	}
    }
  
  cout << "track not found " << label << " " << event->GetNumberOfTracks() << endl;
  return NULL;
}
///________________________________________________________________________
Bool_t AliDalitzElectronCuts::RejectSharedElecGamma(TList *photons, Int_t indexEle){


     for(Int_t i = 0;i<photons->GetEntries();i++){

      AliAODConversionPhoton *photonComp = (AliAODConversionPhoton*) photons->At(i);

      Int_t posLabel = photonComp->GetTrackLabelPositive();
      Int_t negLabel = photonComp->GetTrackLabelNegative();

      if( (photonComp->GetConversionRadius() < fRadiusCut) && (posLabel == indexEle || negLabel == indexEle) ){
        return kFALSE;
      }
     }

   return kTRUE;
}
Bool_t AliDalitzElectronCuts::MassCut(Double_t pi0CandidatePt , Double_t vphotonCandidateMass){
  
	if( pi0CandidatePt < fMassCutPtMin ){
	  
	      if( vphotonCandidateMass < fMassCutLowPt ){
		    return kTRUE;
	      }
		
	}
	else{
	  
	       if( vphotonCandidateMass < fMassCutHighPt ){
		    return kTRUE;
	       }
	      
	}
	
	return kFALSE;
	  
}

Double_t AliDalitzElectronCuts::GetNFindableClustersTPC(AliESDtrack* lTrack){
  
  
  Double_t clsToF=0;
  
  if( fUseCrossedRows == kFALSE ) {

    if ( !fUseCorrectedTPCClsInfo ){
        if(lTrack->GetTPCNclsF()!=0){

              clsToF = (Double_t)lTrack->GetNcls(1)/(Double_t)lTrack->GetTPCNclsF();
        }// Ncluster/Nfindablecluster
    }
    else {

              //clsToF = lTrack->GetTPCClusterInfo(2,0,GetFirstTPCRow(photon->GetConversionRadius()));
              clsToF = lTrack->GetTPCClusterInfo(2,0); //NOTE ask friederike
                
    }
  } else  {
   
	 Float_t nCrossedRowsTPC = lTrack->GetTPCCrossedRows();
	 clsToF = 1.0;
	  if ( lTrack->GetTPCNclsF()>0 ) {
	      clsToF = nCrossedRowsTPC / lTrack->GetTPCNclsF();
	  }	
    }
  
  return clsToF;
   
}

/*
Double_t AliDalitzElectronCuts::GetPsiPair( const AliESDtrack *trackPos, const AliESDtrack *trackNeg )
{
//
// This angle is a measure for the contribution of the opening in polar
// direction ??0 to the opening angle ?? Pair
//
// Ref. Measurement of photons via conversion pairs with the PHENIX experiment at RHIC
//      Master Thesis. Thorsten Dahms. 2005
// https://twiki.cern.ch/twiki/pub/ALICE/GammaPhysicsPublications/tdahms_thesis.pdf
//
        Double_t momPos[3];
        Double_t momNeg[3];
        if( trackPos->GetConstrainedPxPyPz(momPos) == 0 ) trackPos->GetPxPyPz( momPos );
        if( trackNeg->GetConstrainedPxPyPz(momNeg) == 0 ) trackNeg->GetPxPyPz( momNeg );

        TVector3 posDaughter;
        TVector3 negDaughter;

        posDaughter.SetXYZ( momPos[0], momPos[1], momPos[2] );
        negDaughter.SetXYZ( momNeg[0], momNeg[1], momNeg[2] );

        Double_t deltaTheta = negDaughter.Theta() - posDaughter.Theta();
        Double_t openingAngle =  posDaughter.Angle( negDaughter );  //TMath::ACos( posDaughter.Dot(negDaughter)/(negDaughter.Mag()*posDaughter.Mag()) );

        if( openingAngle < 1e-20 ) return 0.;

        Double_t psiAngle = TMath::ASin( deltaTheta/openingAngle );

        return psiAngle;
}*/

Bool_t AliDalitzElectronCuts::IsFromGammaConversion( Double_t psiPair, Double_t deltaPhi )
{
//
// Returns true if it is a gamma conversion according to psi pair value
//
        return ( (deltaPhi > fDeltaPhiCutMin  &&  deltaPhi < fDeltaPhiCutMax) &&
        TMath::Abs(psiPair) < ( fPsiPairCut - fPsiPairCut/fDeltaPhiCutMax * deltaPhi ) );
}

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::UpdateCutString(cutIds cutID, Int_t value) {
///Update the cut string (if it has been created yet)

  if(fCutString && fCutString->GetString().Length() == kNCuts) {
         cout << "Updating cut id in spot number " << cutID << " to " << value << endl;
	fCutString->SetString(GetCutNumber());
  } else {
         cout << "fCutString not yet initialized, will not be updated" << endl;
	return kFALSE;
  }
 // cout << fCutString->GetString().Data() << endl;
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::InitializeCutsFromCutString(const TString analysisCutSelection ) {
  fCutStringRead = Form("%s",analysisCutSelection.Data());
  
   // Initialize Cuts from a given Cut string

//   out<<"Set Cut Number: "<<analysisCutSelection.Data()<<endl;
  AliInfo(Form("Set ElectronCuts Number: %s",analysisCutSelection.Data()));
  
  if(analysisCutSelection.Length()!=kNCuts) {
	AliError(Form("Cut selection has the wrong length! size is %d, number of cuts is %d", analysisCutSelection.Length(), kNCuts));
	return kFALSE;
  }
  if(!analysisCutSelection.IsAlnum()){
	AliError("Cut selection is not alphanumeric");
	return kFALSE;
  }
  
  TString analysisCutSelectionLowerCase = Form("%s",analysisCutSelection.Data());
  analysisCutSelectionLowerCase.ToLower();
  const char *cutSelection = analysisCutSelectionLowerCase.Data();
  #define ASSIGNARRAY(i)  fCuts[i] = ((int)cutSelection[i]>=(int)'a') ? cutSelection[i]-'a'+10 : cutSelection[i]-'0'
  for(Int_t ii=0;ii<kNCuts;ii++){
      ASSIGNARRAY(ii);
  }

 
  // Set Individual Cuts
  for(Int_t ii=0;ii<kNCuts;ii++){
      if(!SetCut(cutIds(ii),fCuts[ii]))return kFALSE;
  }

  //PrintCuts();

    return kTRUE;
}
///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetCut(cutIds cutID, const Int_t value) {
  ///Set individual cut ID

  //cout << "Updating cut  " << fgkCutNames[cutID] << " (" << cutID << ") to " << value << endl;

  switch (cutID) {
 
    case kMaxChi2TPCConstrainedGlobal:
	if( SetMaxChi2TPCConstrainedGlobal( value ) ){
	  fCuts[kMaxChi2TPCConstrainedGlobal] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;
	
  case kededxSigmaITSCut:
	if( SetITSdEdxCutElectronLine(value)) { //NOTE SetITSdEdxCutElectronLine: To be implemented 
	  fCuts[kededxSigmaITSCut] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

	case kededxSigmaTPCCut:
		if( SetTPCdEdxCutElectronLine(value)) { //NOTE SetITSdEdxCutElectronLine: To be implemented 
			fCuts[kededxSigmaTPCCut] = value;
			UpdateCutString(cutID, value);
			return kTRUE;
		} else return kFALSE;
		
  case kpidedxSigmaTPCCut:
		if( SetTPCdEdxCutPionLine(value)) { //NOTE SetITSdEdxCutPionLine: To be implemented
			fCuts[kpidedxSigmaTPCCut] = value;
			UpdateCutString(cutID, value);
			return kTRUE;
		} else return kFALSE;
		
  case kpiMinMomdedxSigmaTPCCut:
		if( SetMinMomPiondEdxTPCCut(value)) {
			fCuts[kpiMinMomdedxSigmaTPCCut] = value;
			UpdateCutString(cutID, value);
			return kTRUE;
		} else return kFALSE;
		
  case kpiMaxMomdedxSigmaTPCCut:
		if( SetMaxMomPiondEdxTPCCut(value)) {
			fCuts[kpiMaxMomdedxSigmaTPCCut] = value;
			UpdateCutString(cutID, value);
			return kTRUE;
		} else return kFALSE;
		
  case kLowPRejectionSigmaCut:
		if( SetLowPRejectionCuts(value) ) {
			fCuts[kLowPRejectionSigmaCut] = value;
			UpdateCutString(cutID, value);
			return kTRUE;
		} else return kFALSE;
		
		
  case kTOFelectronPID:
		if( SetTOFElectronPIDCut(value)) {
			fCuts[kTOFelectronPID] = value;
			UpdateCutString(cutID, value);
			return kTRUE;
		} else return kFALSE;
  case kclsITSCut:
		if( SetITSClusterCut(value) ) {
			fCuts[kclsITSCut] = value;
			UpdateCutString(cutID, value);
			return kTRUE;		 	
		} else return kFALSE;
  case kclsTPCCut:
		if( SetTPCClusterCut(value)) {
			fCuts[kclsTPCCut] = value;
			UpdateCutString(cutID, value);
			return kTRUE;
		} else return kFALSE;
		
  case ketaCut:
		if( SetEtaCut(value)) {
			fCuts[ketaCut] = value;
			UpdateCutString(cutID, value);
			return kTRUE;
		} else return kFALSE;
  case kptCut: 	
		if( SetPtCut(value)) {
			fCuts[kptCut] = value;
			UpdateCutString(cutID, value);
			return kTRUE;
		} else return kFALSE;
    
  case kDCACut:
		if( SetDCACut(value)) {
			fCuts[kDCACut] = value;
			UpdateCutString(cutID, value);
			return kTRUE;
		} else return kFALSE;
	      
	      
	case kPsiPair:
		if( SetPsiPairCut(value)) {
			fCuts[kPsiPair] = value;
			UpdateCutString(cutID, value);
			return kTRUE;
		} else return kFALSE;
		
  case kRejectSharedElecGamma:
		if( SetRejectSharedElecGamma(value)) {
			fCuts[kRejectSharedElecGamma] = value;
			UpdateCutString(cutID, value);
          return kTRUE;
		} else return kFALSE;
		
  case kMaxChi2PerClusterTPC:
		if( SetMaxChi2PerClusterTPC(value)) {
			fCuts[kMaxChi2PerClusterTPC] = value;
			UpdateCutString(cutID, value);
			return kTRUE;
		} else return kFALSE;

  case kMaxChi2PerClusterITS:
		if( SetMaxChi2PerClusterITS(value)) {
			fCuts[kMaxChi2PerClusterITS] = value;
			UpdateCutString(cutID, value);
			return kTRUE;
		} else return kFALSE;
		
  case kmassCut:
		if( SetMassCut(value)) {
			fCuts[kmassCut] = value;
			UpdateCutString(cutID, value);
			return kTRUE;
		} else return kFALSE;
  case kWeights:
                if( SetDoWeights(value)) {
                        fCuts[kWeights] = value;
                        UpdateCutString(cutID, value);
                        return kTRUE;
                } else return kFALSE;
               
  case kuseVPhotonMCPSmearing:
    
		  if( SetUseVPhotonMCPmearing(value)) {
                        fCuts[kuseVPhotonMCPSmearing] = value;
                        UpdateCutString(cutID, value);
                        return kTRUE;
                  } else return kFALSE;
		  
  case kNCuts:
		cout << "Error:: Cut id out of range"<< endl;
		return kFALSE;
  }
	
  cout << "Error:: Cut id " << cutID << " not recognized "<< endl;
  return kFALSE;

  //PrintCuts();
  
}

///________________________________________________________________________

void AliDalitzElectronCuts::PrintCuts() {
    // Print out current Cut Selection
  for(Int_t ic = 0; ic < kNCuts; ic++) {
	printf("%-30s : %d \n", fgkCutNames[ic], fCuts[ic]);
  }

}
///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetMaxChi2TPCConstrainedGlobal(Int_t maxChi2)
{
    if( !fesdTrackCuts ) {

	       cout<<"Warning: AliESDtrackCut is not initialized "<<endl;
	       return kFALSE;
    }    
    
    switch( maxChi2 ){
      
      case 0: fesdTrackCuts->SetMaxChi2TPCConstrainedGlobal(1e10);
	       break;
      case 1: fesdTrackCuts->SetMaxChi2TPCConstrainedGlobal(25.);
	       break;
      case 2: fesdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36.);
	       break;
      case 3: fesdTrackCuts->SetMaxChi2TPCConstrainedGlobal(49.);
	       break;
      case 4: fesdTrackCuts->SetMaxChi2TPCConstrainedGlobal(100.);
	       break;
		
      default:  cout<<"Warning: maxChi2 is not defined"<<maxChi2<<endl;
		return kFALSE;
		  
    }
  
    return kTRUE;
}
///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetITSdEdxCutElectronLine(Int_t ededxSigmaCut)
{   // Set Cut

	switch(ededxSigmaCut){

	case 0: 
		fDodEdxSigmaITSCut = kFALSE;
		fPIDnSigmaBelowElectronLineITS=-100;
		fPIDnSigmaAboveElectronLineITS= 100;
		break;
	case 1: // -10,10
		fDodEdxSigmaITSCut = kTRUE;
		fPIDnSigmaBelowElectronLineITS=-10;
		fPIDnSigmaAboveElectronLineITS=10;
		break;
	case 2: // -6,7
		fDodEdxSigmaITSCut = kTRUE;
		fPIDnSigmaBelowElectronLineITS=-6;
		fPIDnSigmaAboveElectronLineITS=7;
		break;
	case 3: // -5,5
		fDodEdxSigmaITSCut = kTRUE;
		fPIDnSigmaBelowElectronLineITS=-5;
		fPIDnSigmaAboveElectronLineITS=5;
		break;
	case 4: // -4,5
		fDodEdxSigmaITSCut = kTRUE;
		fPIDnSigmaBelowElectronLineITS=-4;
		fPIDnSigmaAboveElectronLineITS=5;
		break;
	case 5: // -3,5
		fDodEdxSigmaITSCut = kTRUE;
		fPIDnSigmaBelowElectronLineITS=-3;
		fPIDnSigmaAboveElectronLineITS=5;
		break;
	case 6: // -4,4
		fDodEdxSigmaITSCut = kTRUE;
		fPIDnSigmaBelowElectronLineITS=-4;
		fPIDnSigmaAboveElectronLineITS=4;
		break;
	case 7: // -2.5,4
		fDodEdxSigmaITSCut = kTRUE;
		fPIDnSigmaBelowElectronLineITS=-2.5;
		fPIDnSigmaAboveElectronLineITS=4;
		break;
	case 8: // -2,3.5
		fDodEdxSigmaITSCut = kTRUE;
		fPIDnSigmaBelowElectronLineITS=-2;
		fPIDnSigmaAboveElectronLineITS=3.5;
		break;
	default:
		cout<<"Warning: ITSdEdxCutElectronLine not defined"<<ededxSigmaCut<<endl;
		return kFALSE;
    
	}
	return kTRUE;
}

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetTPCdEdxCutElectronLine(Int_t ededxSigmaCut)
{   // Set Cut
	switch(ededxSigmaCut){

	case 0: fDodEdxSigmaTPCCut = kFALSE;
		fPIDnSigmaBelowElectronLineTPC=-10;
		fPIDnSigmaAboveElectronLineTPC=10;
		break;
	case 1: // -10,10
		fDodEdxSigmaTPCCut = kTRUE;
		fPIDnSigmaBelowElectronLineTPC=-10;
		fPIDnSigmaAboveElectronLineTPC=10;
		break;
	case 2: // -6,7
		fDodEdxSigmaTPCCut = kTRUE;
		fPIDnSigmaBelowElectronLineTPC=-6;
		fPIDnSigmaAboveElectronLineTPC=7;
		break;
	case 3: // -5,5
		fDodEdxSigmaTPCCut = kTRUE;
		fPIDnSigmaBelowElectronLineTPC=-5;
		fPIDnSigmaAboveElectronLineTPC=5;
		break;
	case 4: // -4,5
		fDodEdxSigmaTPCCut = kTRUE;
		fPIDnSigmaBelowElectronLineTPC=-4;
		fPIDnSigmaAboveElectronLineTPC=5;
		break;	
	case 5: // -3,5
		fDodEdxSigmaTPCCut = kTRUE;
		fPIDnSigmaBelowElectronLineTPC=-3;
		fPIDnSigmaAboveElectronLineTPC=5;
		break;
	case 6: // -4,4
		fDodEdxSigmaTPCCut = kTRUE;
		fPIDnSigmaBelowElectronLineTPC=-4;
		fPIDnSigmaAboveElectronLineTPC=4;
		break;
	case 7: // -2.5,4
		fDodEdxSigmaTPCCut = kTRUE;
		fPIDnSigmaBelowElectronLineTPC=-2.5;
		fPIDnSigmaAboveElectronLineTPC=4;
		break;
	case 8: // -2,3.5
		fDodEdxSigmaTPCCut = kTRUE;
		fPIDnSigmaBelowElectronLineTPC=-2;
		fPIDnSigmaAboveElectronLineTPC=3.5;
		break;
	default:
		cout<<"Warning: TPCdEdxCutElectronLine not defined"<<ededxSigmaCut<<endl;
		return kFALSE;
    
	}
	return kTRUE;
}

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetTPCdEdxCutPionLine(Int_t pidedxSigmaCut)
{   // Set Cut

	switch(pidedxSigmaCut){
                
	case 0: fPIDnSigmaAbovePionLineTPC= 0;
		fPIDnSigmaAbovePionLineTPCHighPt=-100;
		break;
	case 1:  // -10
		fPIDnSigmaAbovePionLineTPC=3.0;            //Update Sep-05-2013 from -10 to 3
		fPIDnSigmaAbovePionLineTPCHighPt=-10;
		break;
	case 2:  // 1
		fPIDnSigmaAbovePionLineTPC=2;              //Update Sep-09-2013 from -1  to  2
		fPIDnSigmaAbovePionLineTPCHighPt=-1;       //Update Sep-09-2013 from -10 to -1
		break;
	case 3:   // 0
		fPIDnSigmaAbovePionLineTPC=2;              //Update Sep-09-2013 from   0 to   2
		fPIDnSigmaAbovePionLineTPCHighPt=0;        //Update Sep-09-2013 from -10 to   0
		break;
	case 4:  // 1
		fPIDnSigmaAbovePionLineTPC=1;
		fPIDnSigmaAbovePionLineTPCHighPt=-10;
		break;
	case 5:  // 1
		fPIDnSigmaAbovePionLineTPC=2.;
		fPIDnSigmaAbovePionLineTPCHighPt=-10;
		break;
	case 6:  // 1
		fPIDnSigmaAbovePionLineTPC=2.5;
		fPIDnSigmaAbovePionLineTPCHighPt=-10;
		break;
	case 7:
		fPIDnSigmaAbovePionLineTPC = 2.0; // We need a bit less tight cut on dE/dx //Updated from 3.0 and -10 to +2.0 , +2.0 
		fPIDnSigmaAbovePionLineTPCHighPt = 2.0;
		break;
	case 8:  // 1
		fPIDnSigmaAbovePionLineTPC = 1.5;   // Updated May-16-2013 from 3.5 and -10 to +1.5, +1
		fPIDnSigmaAbovePionLineTPCHighPt = 1.0;
		break;
	case 9:  // 1
		fPIDnSigmaAbovePionLineTPC=1.5;
		fPIDnSigmaAbovePionLineTPCHighPt=-1.0;
		break;
	default:
		cout<<"Warning: pidedxSigmaCut not defined "<<pidedxSigmaCut<<endl;
		return kFALSE;
	}
	return kTRUE;
}

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetMinMomPiondEdxTPCCut(Int_t piMomdedxSigmaCut)
{   // Set Cut
	switch(piMomdedxSigmaCut){
                
	case 0: fPIDMinPnSigmaAbovePionLineTPC=0.;
		break;
	case 1:  // 50.0 GeV
		fPIDMinPnSigmaAbovePionLineTPC=50.;
		break;
	case 2:  // 20.0 GeV
		fPIDMinPnSigmaAbovePionLineTPC=20.;
		break;
	case 3:  // 1.5 GeV
		fPIDMinPnSigmaAbovePionLineTPC=1.5;
		break;
	case 4:  // 1. GeV
		fPIDMinPnSigmaAbovePionLineTPC=1.;
		break;	
	case 5:  // 0.5 GeV
		fPIDMinPnSigmaAbovePionLineTPC=0.5;
		break;
	case 6:  // 0.4 GeV
		fPIDMinPnSigmaAbovePionLineTPC=0.4;
		break;    
	case 7:  // 0.3 GeV
		fPIDMinPnSigmaAbovePionLineTPC=0.3;
		break;
	case 8:  // 0.25 GeV
		fPIDMinPnSigmaAbovePionLineTPC=0.25;
		break;
	default:
		cout<<"Warning: piMomdedxSigmaCut not defined "<<piMomdedxSigmaCut<<endl;
		return kFALSE;
	}
	return kTRUE;
}
///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetITSClusterCut(Int_t clsITSCut){

    
	if( !fesdTrackCuts ) {

		cout<<"Warning: AliESDtrackCut is not initialized "<<endl;
		return kFALSE;
	}

	switch(clsITSCut){

	case 0: fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
		break;
	case 1: fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
		break;  //1 hit first layer of SPD
	case 2: fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
		break; //1 hit in any layer of SPD
	case 3: fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
		fesdTrackCuts->SetMinNClustersITS(4);
		// 4 hits in total in the ITS. At least 1 hit in the first layer of SPD  
		break;
	case 4: fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
		fesdTrackCuts->SetMinNClustersITS(3);
		// 3 hits in total in the ITS. At least 1 hit in any layer of SPD
		break;
        case 5: fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
                fesdTrackCuts->SetMinNClustersITS(4);
                // 4 hits in total in the ITS. At least 1 hit in any layer of SPD
                break;
        case 6: fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
                fesdTrackCuts->SetMinNClustersITS(5);
                // 5 hits in total in the ITS. At least 1 hit in any layer of SPD
                break;
	case 7: fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
                fesdTrackCuts->SetMinNClustersITS(4);
		break;
	case 8: fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kBoth);
                break;
        case 9: fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kBoth);
		fesdTrackCuts->SetMinNClustersITS(4);
                break;
		
	default:
		cout<<"Warning: clsITSCut not defined "<<clsITSCut<<endl;
		return kFALSE;
	}
	
return kTRUE;
}

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetTPCClusterCut(Int_t clsTPCCut)
{   // Set Cut

	//Update function for systematics 2015-10-08
	
	switch(clsTPCCut){
	  
	case 0: // 0
		fMinClsTPC= 0.;
		fMinClsTPCToF = 0.;
		fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		fesdTrackCuts->SetMinNCrossedRowsTPC(fMinClsTPC);
		break;
		
	case 1:  // 70
		fMinClsTPC= 70.;
		fMinClsTPCToF = 0.7;
		 
		if( fUseCrossedRows ){

		  fesdTrackCuts->SetMinNCrossedRowsTPC(fMinClsTPC);
		  fesdTrackCuts->SetMinNClustersTPC(0);
		
		} else {
		  
		  fesdTrackCuts->SetMinNCrossedRowsTPC(0);
		  fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		  fUseCorrectedTPCClsInfo=0;
		}
		break;
	case 2:  // 80
		fMinClsTPC = 70.;
		fMinClsTPCToF = 0.9;
		
		if( fUseCrossedRows ){

		  fesdTrackCuts->SetMinNCrossedRowsTPC(fMinClsTPC);
		  fesdTrackCuts->SetMinNClustersTPC(0);
		
		} else {
		  
		  fesdTrackCuts->SetMinNCrossedRowsTPC(0);
		  fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		  fUseCorrectedTPCClsInfo=0;
		  
		}
		break;
	case 3:  // Changed 2014-02-04  before fMinClsTPC = 50.;
		fMinClsTPC = 70;
		fMinClsTPCToF = 0.8;
		 
		if( fUseCrossedRows ){
		  
		  fesdTrackCuts->SetMinNCrossedRowsTPC(fMinClsTPC);
		  fesdTrackCuts->SetMinNClustersTPC(0);
		  
		} else {
		  
		  fesdTrackCuts->SetMinNCrossedRowsTPC(0);
		  fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		  fUseCorrectedTPCClsInfo=0;
		  
		}
		//fUseCrossedRows = kTRUE;
		break;
	case 4:  // 0% of findable clusters
		fMinClsTPC = 90;
		fMinClsTPCToF = 0.8;
		
	        if( fUseCrossedRows ){
		  
		  fesdTrackCuts->SetMinNCrossedRowsTPC(fMinClsTPC);
		  fesdTrackCuts->SetMinNClustersTPC(0);
		  
		} else {
		  
		  fesdTrackCuts->SetMinNCrossedRowsTPC(0);
		  fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		  fUseCorrectedTPCClsInfo=0;
		  
		}
		break;
	case 5:  // 35% of findable clusters
		fMinClsTPC = 70;
		fMinClsTPCToF = 0.35;
		
		if( fUseCrossedRows ){
		  
		  fesdTrackCuts->SetMinNCrossedRowsTPC(fMinClsTPC);
		  fesdTrackCuts->SetMinNClustersTPC(0);
		  
		} else {
		  
		  fesdTrackCuts->SetMinNCrossedRowsTPC(0);
		  fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		  fUseCorrectedTPCClsInfo=0;
		  
		}
		break;
	case 6:  // 60% of findable clusters
		fMinClsTPC = 70;
		fMinClsTPCToF = 0.60;

		if( fUseCrossedRows ){
		  
		  fesdTrackCuts->SetMinNCrossedRowsTPC(fMinClsTPC);
		  fesdTrackCuts->SetMinNClustersTPC(0);
		  
		} else {
		  
		  fesdTrackCuts->SetMinNCrossedRowsTPC(0);
		  fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		  fUseCorrectedTPCClsInfo=0;
		  
		}
		break;
	case 7:  // 60% Changed 2014-02-04 before fMinClsTPC = 0.7 fUseCorrectedTPCClsInfo = 0
		 // Changed 2014-02-04  before fMinClsTPC = 50.;
		fMinClsTPC = 90;
		fMinClsTPCToF = 0.35;

		if( fUseCrossedRows ){
		  
		  fesdTrackCuts->SetMinNCrossedRowsTPC(fMinClsTPC);
		  fesdTrackCuts->SetMinNClustersTPC(0);
		  
		} else {
		  
		  fesdTrackCuts->SetMinNCrossedRowsTPC(0);
		  fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		  fUseCorrectedTPCClsInfo=0;
		  
		}
		break;
		
	case 8: fMinClsTPC = 0;
		fMinClsTPCToF = 0.35;

		if( fUseCrossedRows ){
		  
		  fesdTrackCuts->SetMinNCrossedRowsTPC(fMinClsTPC);
		  fesdTrackCuts->SetMinNClustersTPC(0);
		  
		} else {
		  
		  fesdTrackCuts->SetMinNCrossedRowsTPC(0);
		  fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		  fUseCorrectedTPCClsInfo=0;
		  
		}
		break;
		
	case 9:  // 35% of findable clusters
		fMinClsTPC = 70;
		fMinClsTPCToF = 0.35;
		
		if( fUseCrossedRows ){
		  
		  fesdTrackCuts->SetMinNCrossedRowsTPC(fMinClsTPC);
		  fesdTrackCuts->SetMinNClustersTPC(0);
		  
		} else {
		  
		  fesdTrackCuts->SetMinNCrossedRowsTPC(0);
		  fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		  fUseCorrectedTPCClsInfo=1;
		  
		}
		break;
	  
	default:
		cout<<"Warning: clsTPCCut not defined "<<clsTPCCut<<endl;
		return kFALSE;
	}
	return kTRUE;
}

/*
///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetTPCClusterCut(Int_t clsTPCCut)
{   // Set Cut
	switch(clsTPCCut){
	case 0: // 0
		fMinClsTPC= 0.;
		fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		break;
	case 1:  // 70
		fMinClsTPC= 70.;
		fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		break;
	case 2:  // 80
		fMinClsTPC= 80.;
		fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		break;
	case 3:  // Changed 2014-02-04  before fMinClsTPC = 50.;
		fMinClsTPCToF = 0.8;
		fesdTrackCuts->SetMinNCrossedRowsTPC(70);
		fesdTrackCuts->SetMinNClustersTPC(0);
		fUseCrossedRows = kTRUE;
		break;
	case 4:  // 0% of findable clusters
	        fMinClsTPC= 70.;  
		fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		fMinClsTPCToF= 0.0;
		fUseCorrectedTPCClsInfo=0;
		break;
	case 5:  // 35% of findable clusters
		fMinClsTPC = 70.;  
		fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		fMinClsTPCToF= 0.35;
		fUseCorrectedTPCClsInfo=0;
		break;
	case 6:  // 60% of findable clusters
		fMinClsTPC= 70.;  
		fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		fMinClsTPCToF= 0.6;
		fUseCorrectedTPCClsInfo=0;
		break;
	case 7:  // 60% Changed 2014-02-04 before fMinClsTPC = 0.7 fUseCorrectedTPCClsInfo = 0
		 // Changed 2014-02-04  before fMinClsTPC = 50.;
		fMinClsTPCToF = 0.6;
		fesdTrackCuts->SetMinNCrossedRowsTPC(70);
		fesdTrackCuts->SetMinNClustersTPC(0);
		fUseCrossedRows = kTRUE;
		break;
	case 8: fMinClsTPC = 0.;  
		fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		fMinClsTPCToF= 0.35;
		fUseCorrectedTPCClsInfo=0;
		break;
	case 9:  // 35% of findable clusters
		fMinClsTPC = 70.;  
		fesdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
		fMinClsTPCToF= 0.35;
		fUseCorrectedTPCClsInfo=1;
		break;
	  
	default:
		cout<<"Warning: clsTPCCut not defined "<<clsTPCCut<<endl;
		return kFALSE;
	}
	return kTRUE;
}*/

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetEtaCut(Int_t etaCut)
{ 
  // Set eta Cut
	switch(etaCut){
	case 0: 
		fEtaCut = 100.;
		fDoEtaCut = kFALSE;
		break;
	case 1:	// 1.4
		fEtaCut	= 1.4;
		fDoEtaCut = kTRUE;
		break;
	case 2:	// 1.2
		fEtaCut	= 1.2;
		fDoEtaCut = kTRUE;
		break;
	case 3: // 0.9
		fEtaCut	= 0.9;
		fDoEtaCut = kTRUE;
		break;
	case 4: // 0.8
		fEtaCut	= 0.8;
		fDoEtaCut = kTRUE;
		break;
	case 5: // 0.75
		fEtaCut	= 0.75;
		fDoEtaCut = kTRUE;
		break;
	case 6: //0.6
		fEtaCut = 0.6; //changed from 0.4 to 0.6 2013.06.10
		fDoEtaCut = kTRUE;
		break;
	case 7: //0.5
		fEtaCut = 0.5; //changed from 0.3 to 0.5 2013.06.10
		fDoEtaCut = kTRUE;
		break;
	case 8: fEtaCut = 0.4;
		fDoEtaCut = kTRUE;
		break;
	case 9: fEtaCut = 0.65;
		fDoEtaCut = kTRUE; 
		break;
	default:
		cout<<"Warning: EtaCut not defined "<<etaCut<<endl;
		return kFALSE;
	}
	return kTRUE;
}

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetPtCut(Int_t ptCut)
{ 
	// Set Pt Cut
	 //0.1GeV, 0.125 GeV, 0.15 GeV
  
	switch(ptCut){
	  
	case 0: fPtMinCut = 0.075;
		fPtMaxCut = 9999;
		break;
	case 1:	 // 0.1
		fPtMinCut  = 0.1; 	
		fPtMaxCut  = 9999;
		break;
	case 2:	 // 0.125 GeV
		fPtMinCut = 0.125;		
		fPtMaxCut = 9999;
		break;
	case 3: // 0.15 GeV
		fPtMinCut = 0.15;
		fPtMaxCut = 9999;
		break;
		// 0.5 - 0.7 
	case 4: fPtMinCut = 0.5;
		fPtMaxCut = 0.7;
		break;
	case 5: // 0.175 GeV
		fPtMinCut = 0.175;
		fPtMaxCut = 9999;
		break;
	default:
		cout<<"Warning: PtCut not defined "<<ptCut<<endl;
		return kFALSE;
	}
	return kTRUE;
}


///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetDCACut(Int_t dcaCut)
{ 
  // Set DCA Cut
  
	if( !fesdTrackCuts ) {

		cout<<"Warning: AliESDtrackCut is not initialized "<<endl;
		return kFALSE;
	}
  
	switch(dcaCut){
	  
	case 0: //Open cuts//
		fesdTrackCuts->SetMaxDCAToVertexZ(1000);
		fesdTrackCuts->SetMaxDCAToVertexXY(1000);
		break;
	case 1: 
		fesdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01"); //Standard 2010
		fesdTrackCuts->SetMaxDCAToVertexZ(2);
		break;
	case 2: fesdTrackCuts->SetMaxDCAToVertexZ(2);
		fesdTrackCuts->SetMaxDCAToVertexXY(1);
		break; 
	case 3: fesdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1"); //Standard 2011
	        fesdTrackCuts->SetMaxDCAToVertexZ(2);
		break;
	case 4: fesdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0525+0.175/pt^1.1");
		fesdTrackCuts->SetMaxDCAToVertexZ(2);
		break;
        case 5: fesdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
		fesdTrackCuts->SetMaxDCAToVertexZ(1);
                break;
        case 6: fesdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
	        fesdTrackCuts->SetMaxDCAToVertexZ(5);
                break;
	case 7: fesdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1"); //Standard 2011
	        fesdTrackCuts->SetMaxDCAToVertexZ(1);
		break;
	case 8:	fesdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1"); //Standard 2011
	        fesdTrackCuts->SetMaxDCAToVertexZ(5);
		break;
	default:
		cout<<"Warning: dcaCut not defined "<<dcaCut<<endl;
		return kFALSE;
	}
	return kTRUE;
}




///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetMaxMomPiondEdxTPCCut(Int_t piMaxMomdedxSigmaCut)
{   // Set Cut
   switch(piMaxMomdedxSigmaCut){

	 case 0: 
		 fPIDMaxPnSigmaAbovePionLineTPC=0.;
		 break;
	 case 1:  // 100. GeV
		 fPIDMaxPnSigmaAbovePionLineTPC=100.;
		 break;
	 case 2:  // 5. GeV
		 fPIDMaxPnSigmaAbovePionLineTPC=5.;
			break;
	 case 3:  // 4. GeV
		 fPIDMaxPnSigmaAbovePionLineTPC=4.;
		 break;
	 case 4:  // 3.5 GeV
		 fPIDMaxPnSigmaAbovePionLineTPC=3.5;
		 break;
	 case 5:  // 3. GeV
		 fPIDMaxPnSigmaAbovePionLineTPC=3.;
		 break;
	 default:
		 cout<<"Warning: piMaxMomdedxSigmaCut not defined "<<piMaxMomdedxSigmaCut<<endl;
		 return kFALSE;
	 }
	 return kTRUE;
}

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetLowPRejectionCuts(Int_t LowPRejectionSigmaCut)
{   // Set Cut
	switch(LowPRejectionSigmaCut){
	case 0:  //
		fDoKaonRejectionLowP=kFALSE;
		fDoProtonRejectionLowP=kFALSE;
		fDoPionRejectionLowP=kFALSE;
		fPIDnSigmaAtLowPAroundKaonLine=0;
		fPIDnSigmaAtLowPAroundProtonLine=0;
		fPIDnSigmaAtLowPAroundPionLine=0;
		fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLineTPC;
		break;
	case 1:  //
		fDoKaonRejectionLowP=kTRUE;
		fDoProtonRejectionLowP=kTRUE;
		fDoPionRejectionLowP=kTRUE;
		fPIDnSigmaAtLowPAroundKaonLine=0.5;
		fPIDnSigmaAtLowPAroundProtonLine=0.5;
		fPIDnSigmaAtLowPAroundPionLine=0.5;
		fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLineTPC;
		break;
	case 2:  //
		fDoKaonRejectionLowP=kTRUE;
		fDoProtonRejectionLowP=kTRUE;
		fDoPionRejectionLowP=kTRUE;
		fPIDnSigmaAtLowPAroundKaonLine=1.0;
		fPIDnSigmaAtLowPAroundProtonLine=1.0;
		fPIDnSigmaAtLowPAroundPionLine=1.0;
		fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLineTPC;
		break;
	case 3:  //
		fDoKaonRejectionLowP=kTRUE;
		fDoProtonRejectionLowP=kTRUE;
		fDoPionRejectionLowP=kTRUE;
		fPIDnSigmaAtLowPAroundKaonLine=1.5;
		fPIDnSigmaAtLowPAroundProtonLine=1.5;
		fPIDnSigmaAtLowPAroundPionLine=1.5;
		fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLineTPC;
		break;
	case 4:  //
		fDoKaonRejectionLowP=kTRUE;
		fDoProtonRejectionLowP=kTRUE;
		fDoPionRejectionLowP=kTRUE;
		fPIDnSigmaAtLowPAroundKaonLine=2.0;
		fPIDnSigmaAtLowPAroundProtonLine=2.0;
		fPIDnSigmaAtLowPAroundPionLine=2.0;
		fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLineTPC;
		break;
	case 5:  //
		fDoKaonRejectionLowP=kTRUE;
		fDoProtonRejectionLowP=kTRUE;
		fDoPionRejectionLowP=kTRUE;
		fPIDnSigmaAtLowPAroundKaonLine=2.0;
		fPIDnSigmaAtLowPAroundProtonLine=2.0;
		fPIDnSigmaAtLowPAroundPionLine=2.5;
		fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLineTPC;
		break;
	case 6:  //
		fDoKaonRejectionLowP=kTRUE;
		fDoProtonRejectionLowP=kTRUE;
		fDoPionRejectionLowP=kTRUE;
		fPIDnSigmaAtLowPAroundKaonLine=0.;
		fPIDnSigmaAtLowPAroundProtonLine=0.;
		fPIDnSigmaAtLowPAroundPionLine=2.;
		fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLineTPC;
		break;
	case 7: //
		fDoKaonRejectionLowP=kFALSE;
		fDoProtonRejectionLowP=kFALSE;
		fDoPionRejectionLowP=kTRUE;
		fPIDnSigmaAtLowPAroundKaonLine=0.0;
		fPIDnSigmaAtLowPAroundProtonLine=0.0;
		fPIDnSigmaAtLowPAroundPionLine=1.0;
		fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLineTPC;
		break;
	case 8:
		fDoKaonRejectionLowP=kFALSE;
		fDoProtonRejectionLowP=kFALSE;
		fDoPionRejectionLowP=kTRUE;
		fPIDnSigmaAtLowPAroundKaonLine=0.;
		fPIDnSigmaAtLowPAroundProtonLine=0.;
		fPIDnSigmaAtLowPAroundPionLine=0.5; 
		fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLineTPC;
		break;	
	default:
		cout<<"Warning: LowPRejectionSigmaCut not defined "<<LowPRejectionSigmaCut<<endl;
		return kFALSE;
	}
	return kTRUE;
}

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetTOFElectronPIDCut(Int_t TOFelectronPID){
    // Set Cut
	switch(TOFelectronPID){ // RRnewTOF start //////////////////////////////////////////////////////////////////////////
	case 0: // no cut
		fRequireTOF = kFALSE;
		fUseTOFpid = kFALSE;
		fTofPIDnSigmaBelowElectronLine=-100;
		fTofPIDnSigmaAboveElectronLine=100;
		break;
	case 1: // -7,7
		fRequireTOF = kFALSE;
		fUseTOFpid = kTRUE;
		fTofPIDnSigmaBelowElectronLine=-7;
		fTofPIDnSigmaAboveElectronLine=7;
		break;
	case 2: // -5,5
		fRequireTOF = kFALSE;
		fUseTOFpid = kTRUE;
		fTofPIDnSigmaBelowElectronLine=-5;
		fTofPIDnSigmaAboveElectronLine=5;
		break;
	case 3: // -3,5
		fRequireTOF = kFALSE;
		fUseTOFpid = kTRUE;
		fTofPIDnSigmaBelowElectronLine=-3;
		fTofPIDnSigmaAboveElectronLine=5;
		break;
	case 4: // -2,3
		fRequireTOF = kFALSE;
		fUseTOFpid = kTRUE;
		fTofPIDnSigmaBelowElectronLine=-2;
		fTofPIDnSigmaAboveElectronLine=3;
		break;
	case 5: // -3, 3 TOF mandatory
		fRequireTOF = kTRUE;
		fUseTOFpid  = kTRUE;
		fTofPIDnSigmaBelowElectronLine= -3;
		fTofPIDnSigmaAboveElectronLine=  3;
		break;
	default:
        cout<<"Warning: TOFElectronCut not defined "<<TOFelectronPID<<endl;
	return kFALSE;
    } //////////////////////// RRnewTOF end //////////////////////////////////////////////////////////////////////////
    return kTRUE;
}
///_______________________________________________________________________________

Bool_t AliDalitzElectronCuts::SetPsiPairCut(Int_t psiCut) {
  

  switch(psiCut) {
  case 0:
        fDoPsiPairCut = kFALSE;
        fPsiPairCut = 10000.; //
        fDeltaPhiCutMin = -1000.;
        fDeltaPhiCutMax =  1000.;
        
        break;
  case 1:
        fDoPsiPairCut = kTRUE;
        fPsiPairCut = 0.45; // Standard
        fDeltaPhiCutMin = 0.0;
        fDeltaPhiCutMax = 0.12;
        break;
  case 2:
	fDoPsiPairCut = kTRUE;
        fPsiPairCut = 0.60; 
        fDeltaPhiCutMin = 0.0;
        fDeltaPhiCutMax = 0.12;
        break;
  case 3:
        fDoPsiPairCut = kTRUE;
        fPsiPairCut = 0.52;
        fDeltaPhiCutMin = 0.0;
        fDeltaPhiCutMax = 0.12;
	break;
  case 4:
        fDoPsiPairCut = kTRUE;
        fPsiPairCut = 0.30;
        fDeltaPhiCutMin = 0.0;
        fDeltaPhiCutMax = 0.12;
        break;
  case 5:
	fDoPsiPairCut = kTRUE;
        fPsiPairCut = 0.60;
        fDeltaPhiCutMin = 0.0;
        fDeltaPhiCutMax = 0.06;
	break;
  case 6:
	fDoPsiPairCut = kTRUE;
        fPsiPairCut = 0.65;
        fDeltaPhiCutMin = 0.0;
        fDeltaPhiCutMax = 0.14;
	break;
	
    
  default:
      cout<<"Warning: PsiPairCut not defined "<<fPsiPairCut<<endl;
      return kFALSE;
  }

  return kTRUE;
}

///_______________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetRejectSharedElecGamma(Int_t RCut) {
  

  switch(RCut) {
  case 0:
        fDoRejectSharedElecGamma = kFALSE;
        fRadiusCut = 10000; // 
        break;
  case 1:
        fDoRejectSharedElecGamma = kTRUE;
        fRadiusCut = 2.0; // cm 
        break;
  case 2:
        fDoRejectSharedElecGamma = kTRUE;
        fRadiusCut = 3.0; // Standard
        break;
  case 3:
        fDoRejectSharedElecGamma = kTRUE;
        fRadiusCut = 4.0; // 
        break;
  case 4:
        fDoRejectSharedElecGamma = kTRUE;
        fRadiusCut = 5.0; // 
        break;
  case 5:
        fDoRejectSharedElecGamma = kTRUE;
        fRadiusCut = 10.0; // 
        break;
  case 6:
        fDoRejectSharedElecGamma = kTRUE;
        fRadiusCut = 15.0; // 
        break;
  default:
      cout<<"Warning: PsiPairCut not defined "<<fDoRejectSharedElecGamma<<endl;
      return kFALSE;
  }

  return kTRUE;
}
///__________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetMaxChi2PerClusterTPC(Int_t maxChi2){
  
  if( !fesdTrackCuts ) {

	       cout<<"Warning: AliESDtrackCuts is not initialized "<<endl;
	       return kFALSE;
  }    
    

    // Set Cut
    switch(maxChi2){

    case 0:   fesdTrackCuts->SetMaxChi2PerClusterTPC(1e10);
	      break;
    case 1:   fesdTrackCuts->SetMaxChi2PerClusterTPC(3.);
	      break;
    case 2:   fesdTrackCuts->SetMaxChi2PerClusterTPC(4.);
	      break;
    case 3:   fesdTrackCuts->SetMaxChi2PerClusterTPC(5.);
	      break;
    default:
        cout<<"Warning: SetMaxChi2PerClusterTPC not defined "<<maxChi2<<endl;
        return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetMaxChi2PerClusterITS(Int_t maxChi2)
{   // Set Cut

    if( !fesdTrackCuts ) {

	       cout<<"Warning: AliESDtrackCuts is not initialized "<<endl;
	       return kFALSE;
    }    

    switch(maxChi2){
    case 0:
        fesdTrackCuts->SetMaxChi2PerClusterITS(1e10);
        break;
    case 1:
        fesdTrackCuts->SetMaxChi2PerClusterITS(25.);
        break;
    case 2:
         fesdTrackCuts->SetMaxChi2PerClusterITS(36.);
        break;
    case 3:
        fesdTrackCuts->SetMaxChi2PerClusterITS(49.);
        break;
    default:
        cout<<"Warning: SetMaxChi2PerClusterITS not defined "<<maxChi2<<endl;
        return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetDoWeights(Int_t opc)
{   // Set Cut
    switch(opc){
      
    case 0:             fDoWeights = kFALSE;
		        break;
    case 1:             fDoWeights = kTRUE;
		        break; 
    default:
			cout<<"Warning: Weights option not defined "<<opc<<endl;
			return kFALSE;
    }
    return kTRUE;
}
///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetMassCut(Int_t massCut)
{   // Set Cut
    switch(massCut){
      
    case 0:
                        
                        fMassCutPtMin  = -999.; //GeV
                        fMassCutLowPt  =  999.; //GeV/c^2
                        fMassCutHighPt =  999.; //GeV/c^2
                        fMassMinCut = -999;
                        fDoMassCut = kFALSE;   
			fDoMassMinCut = kFALSE;
                        break;
    case 1:
                        //fMassCut = 0.135;             //GeV/c^2
                        fMassCutPtMin  = -999.; //GeV
                        fMassCutLowPt  = 0.135; //GeV/c^2
                        fMassCutHighPt = 0.135; //GeV/c^2
                        fMassMinCut = -999;
                        fDoMassCut = kTRUE;
			fDoMassMinCut = kFALSE;
                        break; 
    case 2:
                        //fMassCut = 0.100;     //GeV/c^2
                        fMassCutPtMin  = -999.; //GeV
                        fMassCutLowPt  = 0.100; //GeV/c^2
                        fMassCutHighPt = 0.100; //GeV/c^2
                        fMassMinCut = -999;
                        fDoMassCut = kTRUE;
			fDoMassMinCut = kFALSE;
                        break;
    case 3:
                        /*fMassCut = 0.075;     //GeV/c^2 Changed from Feb 25
                        fMassCutPtMin  = -999.; //GeV
                        fMassCutLowPt  = 0.075; //GeV/c^2
                        fMassCutHighPt = 0.075; //GeV/c^2
                        fDoMassCut = kTRUE;
			fDoMassMinCut = kFALSE;*/
			fMassCutPtMin  = 1.0;   //GeV
                        fMassCutLowPt  = 0.015; //GeV/c^2
                        fMassCutHighPt = 0.035; //GeV/c^2
                        fMassMinCut    = 0.002;
                        fDoMassCut = kTRUE;
			fDoMassMinCut = kTRUE;
                        break;
    case 4:
                        //fMassCut = 0.050;     //GeV/c^2
                        fMassCutPtMin  = -999.; //GeV
                        fMassCutLowPt  = 0.050; //GeV/c^2
                        fMassCutHighPt = 0.050; //GeV/c^2
                        fMassMinCut = -999;
                        fDoMassCut = kTRUE;
			fDoMassMinCut = kFALSE;
                        break;
    case 5:
                        
                        fMassCutPtMin  = -999.; //GeV
                        fMassCutLowPt  = 0.035; //GeV/c^2
                        fMassCutHighPt = 0.035; //GeV/c^2
                        fMassMinCut = -999;
                        fDoMassCut = kTRUE;
			fDoMassMinCut = kFALSE;
                        break;
    case 6:
                        fMassCutPtMin  = -999.; //GeV
                        fMassCutLowPt  = 0.015; //GeV/c^2
                        fMassCutHighPt = 0.015; //GeV/c^2
                        fMassMinCut = -999;
                        fDoMassCut = kTRUE;
			fDoMassMinCut = kFALSE;
                        break;
    case 7:             fMassCutPtMin  = 1.0;   //GeV
                        fMassCutLowPt  = 0.015; //GeV/c^2
                        fMassCutHighPt = 0.035; //GeV/c^2
                        fMassMinCut = -999;
                        fDoMassCut = kTRUE;
			fDoMassMinCut = kFALSE;
                        break;
    case 8:             fMassCutPtMin  = 1.0;   //GeV
                        fMassCutLowPt  = 0.015; //GeV/c^2
                        fMassCutHighPt = 0.050; //GeV/c^2
                        fMassMinCut = -999;
                        fDoMassCut = kTRUE;
			fDoMassMinCut = kFALSE;
			break;
    case 9:             fMassCutPtMin  = 1.0;   //GeV
                        fMassCutLowPt  = 0.025; //GeV/c^2
                        fMassCutHighPt = 0.035; //GeV/c^2
                        fMassMinCut = -999;
                        fDoMassCut = kTRUE;
			fDoMassMinCut = kFALSE;
			break;
    default:
                        cout<<"Warning: MassCut not defined "<<massCut<<endl;
                        return kFALSE;
    }
    return kTRUE;
}
Bool_t AliDalitzElectronCuts::SetUseVPhotonMCPmearing(Int_t useMCPSmearing)
{// Set Cut
   switch(useMCPSmearing){
     
   case 0:
      fUseVPhotonMCPSmearing=kFALSE;
      fUseElectronMCPSmearing=kFALSE;
      break;
   case 1:
      fUseVPhotonMCPSmearing=kTRUE;
      fUseElectronMCPSmearing=kFALSE;
      break;
   case 2:
      fUseVPhotonMCPSmearing=kFALSE;
      fUseElectronMCPSmearing=kTRUE;
      break;
       
      
   default: cout<<"Warning: Virtual Photon SMearing not defined "<<useMCPSmearing<<endl;
	    return kFALSE;
      
   }
   
   return kTRUE;
}


///________________________________________________________________________
TString AliDalitzElectronCuts::GetCutNumber(){
    // returns TString with current cut number
  return fCutStringRead;
}


///________________________________________________________________________
AliDalitzElectronCuts* AliDalitzElectronCuts::GetStandardCuts2010PbPb(){
    //Create and return standard 2010 PbPb cuts
    AliDalitzElectronCuts *cuts=new AliDalitzElectronCuts("StandardCuts2010PbPb","StandardCuts2010PbPb");
    if(!cuts->InitializeCutsFromCutString("9069640364102")){
	cout<<"Warning: Initialization of Standardcuts2010PbPb failed"<<endl;}
    return cuts;
}

///________________________________________________________________________
AliDalitzElectronCuts* AliDalitzElectronCuts::GetStandardCuts2010pp(){
    //Create and return standard 2010 PbPb cuts
    AliDalitzElectronCuts *cuts=new AliDalitzElectronCuts("StandardCuts2010pp","StandardCuts2010pp");
                                          
    if(!cuts->InitializeCutsFromCutString("9069640364102")){
	cout<<"Warning: Initialization of Standardcuts2010pp failed"<<endl;}
     return cuts;
}

