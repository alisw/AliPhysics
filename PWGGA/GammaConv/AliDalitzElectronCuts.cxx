
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
#include "AliStack.h"
#include "TObjString.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "TList.h"
class iostream;

using namespace std;

ClassImp(AliDalitzElectronCuts)


const char* AliDalitzElectronCuts::fgkCutNames[AliDalitzElectronCuts::kNCuts] = {
"GoodId",
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
"BackgroundScheme",
"NumberOfRotations",
};

//________________________________________________________________________
AliDalitzElectronCuts::AliDalitzElectronCuts(const char *name,const char *title) : AliAnalysisCuts(name,title),
    fHistograms(NULL),
    fPIDResponse(NULL),
    fesdTrackCuts(NULL),
    fEtaCut(0.9),
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
    fUseTOFpid(kFALSE),
    fRequireTOF(kFALSE),
    fUseTrackMultiplicityForBG(kFALSE),
    fBKGMethod(0),
    fnumberOfRotationEventsForBG(0),
    fCutString(NULL),
    hCutIndex(NULL),
    hdEdxCuts(NULL),
    hITSdEdxbefore(NULL),
    hITSdEdxafter(NULL),
    hTPCdEdxbefore(NULL),
    hTPCdEdxafter(NULL),
    hTPCdEdxSignalafter(NULL),
    hTOFbefore(NULL),
    hTOFafter(NULL)
   {
    InitPIDResponse();
    for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=0;}
    fCutString=new TObjString((GetCutNumber()).Data());


    //fesdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");

   // Using standard function for setting Cuts
    Bool_t selectPrimaries=kTRUE;
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
    


    TAxis *AxisBeforeITS = NULL;
    TAxis *AxisBeforedEdx = NULL;
    TAxis *AxisBeforeTOF = NULL;

    if(preCut){


       hITSdEdxbefore=new TH2F(Form("Electron_ITS_before %s",cutName.Data()),"ITS dEdx electron before" ,150,0.05,20,400,-10,10);
       fHistograms->Add(hITSdEdxbefore);
       AxisBeforeITS = hITSdEdxbefore->GetXaxis();

       hTPCdEdxbefore=new TH2F(Form("Electron_dEdx_before %s",cutName.Data()),"dEdx electron before" ,150,0.05,20,400,-10,10);
       fHistograms->Add(hTPCdEdxbefore);
       AxisBeforedEdx = hTPCdEdxbefore->GetXaxis();

       hTOFbefore=new TH2F(Form("Electron_TOF_before %s",cutName.Data()),"TOF electron before" ,150,0.05,20,400,-6,10);
       fHistograms->Add(hTOFbefore);
       AxisBeforeTOF = hTOFbefore->GetXaxis();

    }


    hITSdEdxafter=new TH2F(Form("Electron_ITS_after %s",cutName.Data()),"ITS dEdx electron after" ,150,0.05,20,400, -10,10);
    fHistograms->Add(hITSdEdxafter);

    hTPCdEdxafter=new TH2F(Form("Electron_dEdx_after %s",cutName.Data()),"dEdx electron after" ,150,0.05,20,400, -10,10);
    fHistograms->Add(hTPCdEdxafter);

    hTPCdEdxSignalafter=new TH2F(Form("Electron_dEdxSignal_after %s",cutName.Data()),"dEdx electron signal after" ,150,0.0,3.0,200,0.0,200);
    fHistograms->Add(hTPCdEdxSignalafter);

    hTOFafter=new TH2F(Form("Electron_TOF_after %s",cutName.Data()),"TOF electron after" ,150,0.05,20,400,-6,10);
    fHistograms->Add(hTOFafter);

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
    if(preCut){
       AxisBeforeITS->Set(bins, newBins);
       AxisBeforedEdx->Set(bins, newBins);
       AxisBeforeTOF->Set(bins, newBins);
    }
    delete [] newBins;

        
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
Bool_t AliDalitzElectronCuts::ElectronIsSelected(AliESDtrack* lTrack)
{
    //Selection of Reconstructed electrons

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


    // dEdx Cuts
    if( ! dEdxCuts( track ) ) {
         if(hCutIndex)hCutIndex->Fill(kdEdxCuts);
         return kFALSE;

    }

    //Electron passed the cuts
    if(hCutIndex)hCutIndex->Fill(kElectronOut);


    return kTRUE;
}

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::TrackIsSelected(AliESDtrack* lTrack) {
    // Track Selection for Photon Reconstruction



    if( ! fesdTrackCuts->AcceptTrack(lTrack) ){

        return kFALSE;
    }

   if(  TMath::Abs( lTrack->Eta()) > fEtaCut ) {

        return kFALSE;
    }

    if( lTrack->GetNcls(1) < fMinClsTPC ) {

        return kFALSE;
    }

    //Findable clusters

    Double_t clsToF=0;


    if (!fUseCorrectedTPCClsInfo ){
        if(lTrack->GetTPCNclsF()!=0){

              clsToF = (Double_t)lTrack->GetNcls(1)/(Double_t)lTrack->GetTPCNclsF();
        }// Ncluster/Nfindablecluster
    }
    else {

              //clsToF = lTrack->GetTPCClusterInfo(2,0,GetFirstTPCRow(photon->GetConversionRadius()));
              clsToF = lTrack->GetTPCClusterInfo(2,1); //NOTE ask friederike
                
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
   // Initialize Cuts from a given Cut string

  cout<<"Set Cut Number: "<<analysisCutSelection.Data()<<endl;
  if(analysisCutSelection.Length()!=kNCuts) {
	AliError(Form("Cut selection has the wrong length! size is %d, number of cuts is %d", analysisCutSelection.Length(), kNCuts));
	return kFALSE;
  }
  if(!analysisCutSelection.IsDigit()){
	AliError("Cut selection contains characters");
	return kFALSE;
  }
  
  const char *cutSelection = analysisCutSelection.Data();
  #define ASSIGNARRAY(i)	fCuts[i] = cutSelection[i] - '0'
  for(Int_t ii=0;ii<kNCuts;ii++){
      ASSIGNARRAY(ii);
  }

  // TestFlag
  if(fCuts[0] !=9){
    AliError("Analysis Cut Selection does not start with 9");
	PrintCuts();
    return kFALSE;
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
  case kgoodId:
	fCuts[kgoodId] = value;
	if(value != 9) {
	  cout << "Error:: First value of cut string is wrong, aborting!!" << endl;
	  return kFALSE;
	} else {
	  return kTRUE;
	}

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

  case kBackgroundScheme:
        if( SetBackgroundScheme(value)) {
          fCuts[kBackgroundScheme] = value;
          UpdateCutString(cutID, value);
          return kTRUE;
        } else return kFALSE;

  case kNumberOfRotations:
        if( SetNumberOfRotations(value)) {
          fCuts[kNumberOfRotations] = value;
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
			fPIDnSigmaAbovePionLineTPC=-10;
			fPIDnSigmaAbovePionLineTPCHighPt=-10;
		        break;
		case 2:  // 1
			fPIDnSigmaAbovePionLineTPC=-1;
			fPIDnSigmaAbovePionLineTPCHighPt=-10;
			break;
		case 3:   // 0
			fPIDnSigmaAbovePionLineTPC=0;
			fPIDnSigmaAbovePionLineTPCHighPt=-10;
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
			fPIDnSigmaAbovePionLineTPC=3.0; // We need a bit less tight cut on dE/dx
			fPIDnSigmaAbovePionLineTPCHighPt=-10;
			break;
		case 8:  // 1
			fPIDnSigmaAbovePionLineTPC=3.5;
			fPIDnSigmaAbovePionLineTPCHighPt=-10;
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
	default:
        cout<<"Warning: clsITSCut not defined "<<clsITSCut<<endl;
        return kFALSE;
    }

return kTRUE;
}

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetTPCClusterCut(Int_t clsTPCCut)
{   // Set Cut
    switch(clsTPCCut){
    case 0: // 0
	fMinClsTPC= 0.;
	break;
    case 1:  // 70
	fMinClsTPC= 70.;
	break;
    case 2:  // 80
	fMinClsTPC= 80.;
	break;
    case 3:  // 100
	fMinClsTPC= 100.;
	break;
    case 4:  // 0% of findable clusters
	fMinClsTPCToF= 0.0;
	fUseCorrectedTPCClsInfo=0;
	break;
    case 5:  // 35% of findable clusters
	fMinClsTPCToF= 0.35;
	fUseCorrectedTPCClsInfo=0;
	break;
    case 6:  // 60% of findable clusters
	fMinClsTPCToF= 0.6;
	fUseCorrectedTPCClsInfo=0;
	break;
    case 7:  // 70% of findable clusters
	fMinClsTPCToF= 0.7;
	fUseCorrectedTPCClsInfo=0;
	break;
    default:
	cout<<"Warning: clsTPCCut not defined "<<clsTPCCut<<endl;
	return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetEtaCut(Int_t etaCut)
{ 
  // Set eta Cut
	switch(etaCut){
                case 0: fEtaCut         = 100.;
                        break;
		case 1:	// 1.4
			fEtaCut		= 1.4;
			break;
		case 2:	// 1.2
			fEtaCut		= 1.2;
			break;
		case 3: // 0.9
			fEtaCut		= 0.9;
			break;
		case 4: // 0.8
			fEtaCut		= 0.8;
			break;
		case 5: // 0.75
			fEtaCut		= 0.75;
			break;
		default:
			cout<<"Warning: EtaCut not defined "<<etaCut<<endl;
			return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetMaxMomPiondEdxTPCCut(Int_t piMaxMomdedxSigmaCut)
{   // Set Cut
   switch(piMaxMomdedxSigmaCut){

                case 0: fPIDMinPnSigmaAbovePionLineTPC=0.;
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
	break;
    case 1:  //
        fDoKaonRejectionLowP=kTRUE;
        fDoProtonRejectionLowP=kTRUE;
        fDoPionRejectionLowP=kTRUE;
	fPIDnSigmaAtLowPAroundKaonLine=0.5;
	fPIDnSigmaAtLowPAroundProtonLine=0.5;
	fPIDnSigmaAtLowPAroundPionLine=0.5;
	break;
    case 2:  //
        fDoKaonRejectionLowP=kTRUE;
        fDoProtonRejectionLowP=kTRUE;
        fDoPionRejectionLowP=kTRUE;
	fPIDnSigmaAtLowPAroundKaonLine=1;
	fPIDnSigmaAtLowPAroundProtonLine=1;
	fPIDnSigmaAtLowPAroundPionLine=1;
	break;
    case 3:  //
        fDoKaonRejectionLowP=kTRUE;
        fDoProtonRejectionLowP=kTRUE;
        fDoPionRejectionLowP=kTRUE;
	fPIDnSigmaAtLowPAroundKaonLine=1.5;
	fPIDnSigmaAtLowPAroundProtonLine=1.5;
	fPIDnSigmaAtLowPAroundPionLine=1.5;
	break;
    case 4:  //
        fDoKaonRejectionLowP=kTRUE;
        fDoProtonRejectionLowP=kTRUE;
        fDoPionRejectionLowP=kTRUE;
	fPIDnSigmaAtLowPAroundKaonLine=2.0;
	fPIDnSigmaAtLowPAroundProtonLine=2.0;
	fPIDnSigmaAtLowPAroundPionLine=2.0;
	break;
    case 5:  //
        fDoKaonRejectionLowP=kTRUE;
        fDoProtonRejectionLowP=kTRUE;
        fDoPionRejectionLowP=kTRUE;
	fPIDnSigmaAtLowPAroundKaonLine=2.0;
	fPIDnSigmaAtLowPAroundProtonLine=2.0;
	fPIDnSigmaAtLowPAroundPionLine=2.5;
	break;
    case 6:  //
        fDoKaonRejectionLowP=kTRUE;
        fDoProtonRejectionLowP=kTRUE;
        fDoPionRejectionLowP=kTRUE;
	fPIDnSigmaAtLowPAroundKaonLine=0.;
	fPIDnSigmaAtLowPAroundProtonLine=0.;
	fPIDnSigmaAtLowPAroundPionLine=2.;
	break;
    case 7: //
        fDoKaonRejectionLowP=kFALSE;
        fDoProtonRejectionLowP=kFALSE;
        fDoPionRejectionLowP=kTRUE;
        fPIDnSigmaAtLowPAroundKaonLine=0.0;
        fPIDnSigmaAtLowPAroundProtonLine=0.0;
        fPIDnSigmaAtLowPAroundPionLine=1.0;
        break;
    case 8:
        fDoKaonRejectionLowP=kFALSE;
        fDoProtonRejectionLowP=kFALSE;
        fDoPionRejectionLowP=kTRUE;
        fPIDnSigmaAtLowPAroundKaonLine=0.;
        fPIDnSigmaAtLowPAroundProtonLine=0.;
        fPIDnSigmaAtLowPAroundPionLine=0.5; 
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
	fUseTOFpid = kFALSE;
	fTofPIDnSigmaBelowElectronLine=-100;
	fTofPIDnSigmaAboveElectronLine=100;
	break;
    case 1: // -7,7
	fUseTOFpid = kTRUE;
	fTofPIDnSigmaBelowElectronLine=-7;
	fTofPIDnSigmaAboveElectronLine=7;
	break;
    case 2: // -5,5
	fUseTOFpid = kTRUE;
	fTofPIDnSigmaBelowElectronLine=-5;
	fTofPIDnSigmaAboveElectronLine=5;
	break;
    case 3: // -3,5
	fUseTOFpid = kTRUE;
	fTofPIDnSigmaBelowElectronLine=-3;
	fTofPIDnSigmaAboveElectronLine=5;
	break;
    case 4: // -2,3
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
        fDeltaPhiCutMin = 0.;
        fDeltaPhiCutMax = 0.12;
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
Bool_t AliDalitzElectronCuts::SetBackgroundScheme(Int_t BackgroundScheme){

    // Set Cut
    switch(BackgroundScheme){

    case 0: //Rotation
            fBKGMethod = 0;
            fUseTrackMultiplicityForBG = kFALSE;
        break;
    case 1: // mixed event with V0 multiplicity
            fBKGMethod  = 1;
            fUseTrackMultiplicityForBG = kFALSE;
        break;
    case 2: // mixed event with track multiplicity
            fUseTrackMultiplicityForBG = kTRUE;
            fBKGMethod  = 1;
        break;
    case 3: //Rotation
           fUseTrackMultiplicityForBG = kFALSE;
            fBKGMethod  = 2;
        break;
    case 4: //Rotation
            fUseTrackMultiplicityForBG = kTRUE;
            fBKGMethod  = 2;
        break;
    case 5: fUseTrackMultiplicityForBG = kTRUE;
            fBKGMethod  = 3;
        break;

    default:
        cout<<"Warning: BackgroundScheme not defined "<<BackgroundScheme<<endl;
        return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliDalitzElectronCuts::SetNumberOfRotations(Int_t NumberOfRotations)
{   // Set Cut
    switch(NumberOfRotations){
    case 0:
        fnumberOfRotationEventsForBG = 5;
        break;
    case 1:
        fnumberOfRotationEventsForBG = 10;
        break;
    case 2:
        fnumberOfRotationEventsForBG = 15;
        break;
    case 3:
        fnumberOfRotationEventsForBG = 20;
        break;
    case 4:
        fnumberOfRotationEventsForBG = 2;
        break;
    case 5:
        fnumberOfRotationEventsForBG = 50;
        break;
    case 6:
        fnumberOfRotationEventsForBG = 80;
        break;
    case 7:
        fnumberOfRotationEventsForBG = 100;
        break;
    default:
        cout<<"Warning: NumberOfRotations not defined "<<NumberOfRotations<<endl;
        return kFALSE;
    }
    return kTRUE;
}
///________________________________________________________________________
TString AliDalitzElectronCuts::GetCutNumber(){
    // returns TString with current cut number
  TString a(kNCuts);
  for(Int_t ii=0;ii<kNCuts;ii++){
	a.Append(Form("%d",fCuts[ii]));
  }
  return a;
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

