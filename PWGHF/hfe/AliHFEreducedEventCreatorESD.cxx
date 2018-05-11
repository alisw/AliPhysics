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
//
// Debug tree task
// the tree is represented as reduced events
// 
// Authors:
//   M.Fasel <M.Fasel@gsi.de>
//
//
#include <TArrayI.h>
#include <TBits.h>
#include <TFile.h>
#include <TParticle.h>
#include <TString.h>
#include <TTree.h>

#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisUtils.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliInputEventHandler.h"
#include "AliHFEcuts.h"
#include "AliHFEextraCuts.h"
#include "AliHFEmcQA.h"
#include "AliHFEpidTPC.h"
#include "AliHFEreducedEvent.h"
#include "AliHFEreducedTrack.h"
#include "AliHFEreducedMCParticle.h"
#include "AliHFEsignalCuts.h"
#include "AliHFEV0taginfo.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMultiplicity.h"
#include "AliPIDResponse.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliVVZERO.h"
#include "AliVZDC.h"
#include "AliTRDTriggerAnalysis.h"
#include "TTreeStream.h"
#include "AliEventplane.h"
#include <AliOADBContainer.h>


#include "AliHFEreducedEventCreatorESD.h"

ClassImp(AliHFEreducedEventCreatorESD)

AliHFEreducedEventCreatorESD::AliHFEreducedEventCreatorESD():
  AliAnalysisTaskSE(),
  fHFEtree(NULL),
  fAnalysisUtils(NULL),
  fHFEevent(NULL),
  fTrackCuts(NULL),
  fExtraCuts(NULL),
  fSignalCuts(NULL),
  fTPCpid(NULL),
  fV0Tagger(NULL),
  fTRDTriggerAnalysis(NULL),
  fEventNumber(0),
  fNclustersTPC(70),
  fNclustersTPCPID(0),
  fNclustersITS(2),
  fRemoveFirstEvent(kFALSE),
  fFlagPileupEvents(kFALSE),
  fSelectSignalOnly(kFALSE),
  fRunNumber(-1), fRunNumberCaliInfo(-1), fVZEROgainEqualization(0x0), fVZEROApol(0), fVZEROCpol(0), fChi2A(0x0), fChi2C(0x0), fChi3A(0x0), fChi3C(0x0), fSigma2A(0x0), fSigma2C(0x0), fSigma3A(0x0), fSigma3C(0x0), fWeightForVZERO(kChi), fOADB(0x0), fCentralityClasses(0), fInCentralitySelection(-1)
{
  //
  // Default constructor
  //
}

AliHFEreducedEventCreatorESD::AliHFEreducedEventCreatorESD(const char *name):
  AliAnalysisTaskSE(name),
  fHFEtree(NULL),
  fAnalysisUtils(NULL),
  fHFEevent(NULL),
  fTrackCuts(NULL),
  fExtraCuts(NULL),
  fSignalCuts(NULL),
  fTPCpid(NULL),
  fV0Tagger(NULL),
  fTRDTriggerAnalysis(NULL),
  fEventNumber(0),
  fNclustersTPC(70),
  fNclustersTPCPID(0),
  fNclustersITS(2),
  fRemoveFirstEvent(kFALSE),
  fFlagPileupEvents(kFALSE),
  fSelectSignalOnly(kFALSE),
  fRunNumber(-1), fRunNumberCaliInfo(-1), fVZEROgainEqualization(0x0), fVZEROApol(0), fVZEROCpol(0), fChi2A(0x0), fChi2C(0x0), fChi3A(0x0), fChi3C(0x0), fSigma2A(0x0), fSigma2C(0x0), fSigma3A(0x0), fSigma3C(0x0), fWeightForVZERO(kChi), fOADB(0x0), fCentralityClasses(0), fInCentralitySelection(-1)
{
  //
  // Default constructor
  //
  fTPCpid = new AliHFEpidTPC("QAtpcPID");
  fAnalysisUtils = new AliAnalysisUtils;
  fTRDTriggerAnalysis = new AliTRDTriggerAnalysis();
  fTRDTriggerAnalysis->SetRequireMatchElectron(kTRUE);
  fV0Tagger = new AliHFEV0taginfo("Tagger");
  DefineOutput(1, TTree::Class());
}

AliHFEreducedEventCreatorESD::~AliHFEreducedEventCreatorESD(){
  //
  // Default destructor
  //
  if(fAnalysisUtils) delete fAnalysisUtils;
  if(fTPCpid) delete fTPCpid;
  if(fV0Tagger) delete fV0Tagger;
  if(fTRDTriggerAnalysis) delete fTRDTriggerAnalysis;
  if(fHFEevent) delete fHFEevent;
  if(fSignalCuts) delete fSignalCuts;
  if(fTrackCuts) delete fTrackCuts;

}


Int_t AliHFEreducedEventCreatorESD::GetVZEROCentralityBin() const
{
    // return cache index number corresponding to the event centrality
    #ifdef ALIANALYSISTASKJETV2_DEBUG_FLAG_1
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
    Float_t v0Centr(InputEvent()->GetCentrality()->GetCentralityPercentile("V0M"));
    if(v0Centr < 5) return 0;
    else if(v0Centr < 10) return 1;
    else if(v0Centr < 20) return  2;
    else if(v0Centr < 30) return  3;
    else if(v0Centr < 40) return  4;
    else if(v0Centr < 50) return  5;
    else if(v0Centr < 60) return  6;
    else if(v0Centr < 70) return  7;
    else return 8;
}

void AliHFEreducedEventCreatorESD::ReadVZEROCalibration2010h() // From jetv2 task
{
    // necessary for calibration of 10h vzero event plane. code copied from flow package 
    // (duplicate, but i didn't want to introduce an ulgy dependency )
    // this function is only called when the runnumber changes 
    #ifdef ALIANALYSISTASKJETV2_DEBUG_FLAG_1
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif

    // 1) check if the proper chi weights for merging vzero a and vzero c ep are present
    // if not, use sane defaults. centrality binning is equal to that given in the fVZEROcentralityBin snippet
    //
    // when the user wants to, set the weights to 1 (effectively disabling them)
    // chi values can be calculated using the static helper function 
    // AliAnalysisTaskJetV2::CalculateEventPlaneChi(Double_t res) where res is the event plane
    // resolution in a given centrality bin
    // the resolutions that were used for these defaults are
    Double_t chiC2[] = {0.771423, 1.10236, 1.38116, 1.48077, 1.31964, 1.10236, 0.674622, 0.600403, 0.273865};
    Double_t chiA2[] = {0.582214, 0.674622, 0.832214, 0.873962, 0.832214, 0.771423, 0.637146, 0.424255, 0.257385};
    Double_t chiC3[] = {0.493347, 0.493347, 0.458557, 0.407166, 0.356628, 0.273865, 0.176208, 6.10352e-05, 6.10352e-05};
    Double_t chiA3[] = {0.356628, 0.373474, 0.356628, 0.306702, 0.24115, 0.192322, 0.127869, 6.10352e-05, 6.10352e-05};

    if(!fChi2A) fChi2A = new TArrayD(9, chiA2);
    if(!fChi2C) fChi2C = new TArrayD(9, chiC2);
    if(!fChi3A) fChi3A = new TArrayD(9, chiA3);
    if(!fChi3C) fChi3C = new TArrayD(9, chiC3);
   
    Double_t sigmaC2[] = {0.000210563,0.000554248,0.00126934,0.00138031,0.00124522,0.000948494,0.00115442,0.000626186,0.000161246};
    Double_t sigmaA2[] =  {0.000195393,0.000509235,0.00112734,0.00121416,0.00110601,0.00086572,0.0010805,0.000579927,0.00013517};
    Double_t sigmaC3[] = {0.000131573,0.000317261,0.000783971,0.000885244,0.000763271,0.000542612,0.000647701,0.000524767,0};
    Double_t sigmaA3[] = {0.000123304,0.000293338,0.000714463,0.000798547,0.00069079,0.000503398,0.000615878,0.000489984,0};

    if(!fSigma2A) fSigma2A = new TArrayD(9, sigmaA2);
    if(!fSigma2C) fSigma2C = new TArrayD(9, sigmaC2);
    if(!fSigma3A) fSigma3A = new TArrayD(9, sigmaA3);
    if(!fSigma3C) fSigma3C = new TArrayD(9, sigmaC3);

    // 2) check if the database file is open, if not, open it
    if(!fOADB || fOADB->IsZombie()) fOADB = new TFile("$ALICE_PHYSICS/OADB/PWGCF/VZERO/VZEROcalibEP.root"); // changed from TFile::Open()
    if(fOADB->IsZombie()) {
	printf("OADB file $ALICE_PHYSICS/OADB/PWGCF/VZERO/VZEROcalibEP.root cannot be opened, CALIBRATION FAILED !");
	return;
    }

    AliOADBContainer *cont = (AliOADBContainer*) fOADB->Get("hMultV0BefCorr");
    if(!cont){
        // see if database is readable
	printf("OADB object hMultV0BefCorr is not available in the file\n");
	return;	
    }
    Int_t run(fRunNumber);
    if(!(cont->GetObject(run))){
        // if the run isn't recognized fall back to a default run
	printf("OADB object hMultVZEROBefCorr is not available for run %i (used default run 137366)\n",run);
	run = 137366;
    }
    // step 3) get the proper multiplicity weights from the vzero signal
    fVZEROgainEqualization = (TH1*)((TH2F*)cont->GetObject(run))->ProfileX();
    if(!fVZEROgainEqualization) {
        AliFatal(Form("%s: Fatal error, couldn't read fVZEROgainEqualization from OADB object < \n", GetName()));
        return;
    }

    TF1* fpol0 = new TF1("fpol0","pol0");
    fVZEROgainEqualization->Fit(fpol0, "N0", "", 0, 31);
    fVZEROCpol = fpol0->GetParameter(0);
    fVZEROgainEqualization->Fit(fpol0, "N0", "", 32, 64);
    fVZEROApol = fpol0->GetParameter(0);

    // step 4) extract the information to re-weight the q-vectors 
    for(Int_t iside=0;iside<2;iside++){
	for(Int_t icoord=0;icoord<2;icoord++){
	    for(Int_t i=0;i  < 9;i++){
		char namecont[100];
  		if(iside==0 && icoord==0)
		  snprintf(namecont,100,"hQxc2_%i",i);
		else if(iside==1 && icoord==0)
		  snprintf(namecont,100,"hQxa2_%i",i);
		else if(iside==0 && icoord==1)
		  snprintf(namecont,100,"hQyc2_%i",i);
		else if(iside==1 && icoord==1)
		  snprintf(namecont,100,"hQya2_%i",i);

		cont = (AliOADBContainer*) fOADB->Get(namecont);
		if(!cont){
		    printf("OADB object %s is not available in the file\n",namecont);
		    return;	
		}
	
		if(!(cont->GetObject(run))){
		    printf("OADB object %s is not available for run %i (used run 137366)\n",namecont,run);
		    run = 137366;
		}

                // store info for all centralities to cache
                fMeanQ[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetMean();
		fWidthQ[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetRMS();

		//for v3
		if(iside==0 && icoord==0)
		  snprintf(namecont,100,"hQxc3_%i",i);
		else if(iside==1 && icoord==0)
		  snprintf(namecont,100,"hQxa3_%i",i);
		else if(iside==0 && icoord==1)
		  snprintf(namecont,100,"hQyc3_%i",i);
		else if(iside==1 && icoord==1)
		  snprintf(namecont,100,"hQya3_%i",i);

		cont = (AliOADBContainer*) fOADB->Get(namecont);
		if(!cont){
		    printf("OADB object %s is not available in the file\n",namecont);
		    return;	
		}
		
		if(!(cont->GetObject(run))){
		    printf("OADB object %s is not available for run %i (used run 137366)\n",namecont,run);
		    run = 137366;
		}
                // store info for all centralities to cache
		fMeanQv3[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetMean();
		fWidthQv3[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetRMS();
     	    }
	}
    }
    // cleanup. the opened file is closed in the destructor, otherwise fVZEROgainEqualization is no longer available
    delete fpol0;
    // for qa store the runnumber that is currently used for calibration purposes
    fRunNumberCaliInfo = run;
}

void AliHFEreducedEventCreatorESD::CalculateQvectorVZERO(Double_t Qa2[2], Double_t Qc2[2], Double_t Qa3[2], Double_t Qc3[2]) const
{
    // return the calibrated 2nd and 3rd order q-vectors for vzeroa and vzeroc
    // function takes arrays as arguments, which correspond to vzero info in the following way
    // 
    // Qa2[0] = Qx2 for vzero A         Qa2[1] = Qy2 for vzero A (etc)
    
    #ifdef ALIANALYSISTASKJETV2_DEBUG_FLAG_1
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
    // placeholders 
    Double_t phi(-999.), mult(-999.); 
    // reset placeholders for Q-vector components
    Qa2[0] = 0.;    Qc2[0] = 0.;    Qa3[0] = 0.;    Qc3[0] = 0.;
    Qa2[1] = 0.;    Qc2[1] = 0.;    Qa3[1] = 0.;    Qc3[1] = 0.;
    // for qa purposes, save also raw signal
    //Double_t QaX(0), QaY(0), QcX(0), QcY(0);
    for(Int_t i(0); i < 64; i++) {
        // loop over all scintillators, construct Q-vectors in the same loop
        phi     = TMath::PiOver4()*(0.5+i%8);
        mult    = InputEvent()->GetVZEROData()->GetMultiplicity(i);
        //if(fFillQAHistograms) fHistMultVsCellBC->Fill(i, mult);
        // note that disabled rings have already been excluded in ReadVZEROCalibration2010h
        if(i < 32) {    // v0c side
            // fill Q-vectors for v0c side
            Qc2[0] += mult*TMath::Cos(2.*phi)*fVZEROCpol/fVZEROgainEqualization->GetBinContent(1+i);
            Qc3[0] += mult*TMath::Cos(3.*phi)*fVZEROCpol/fVZEROgainEqualization->GetBinContent(1+i);
            Qc2[1] += mult*TMath::Sin(2.*phi)*fVZEROCpol/fVZEROgainEqualization->GetBinContent(1+i);
            Qc3[1] += mult*TMath::Sin(3.*phi)*fVZEROCpol/fVZEROgainEqualization->GetBinContent(1+i);
//             if(fFillQAHistograms) {
//                 fHistMultVsCell->Fill(i, mult*fVZEROCpol/fVZEROgainEqualization->GetBinContent(1+i));
//                 QcX += mult*TMath::Cos(2.*phi);
//                 QcY += mult*TMath::Sin(2.*phi);
//             }
        } else {       // v0a side
            // fill Q-vectors for v0a side
            Qa2[0] += mult*TMath::Cos(2.*phi)*fVZEROApol/fVZEROgainEqualization->GetBinContent(1+i);
            Qa3[0] += mult*TMath::Cos(3.*phi)*fVZEROApol/fVZEROgainEqualization->GetBinContent(1+i);
            Qa2[1] += mult*TMath::Sin(2.*phi)*fVZEROApol/fVZEROgainEqualization->GetBinContent(1+i);
            Qa3[1] += mult*TMath::Sin(3.*phi)*fVZEROApol/fVZEROgainEqualization->GetBinContent(1+i);
//             if(fFillQAHistograms) {
//                 fHistMultVsCell->Fill(i, mult*fVZEROApol/fVZEROgainEqualization->GetBinContent(1+i));
//                 QaX += mult*TMath::Cos(2.*phi);
//                 QaY += mult*TMath::Sin(2.*phi);
//             }
        }
    }
    // get the cache index and read the correction terms from the cache
    Int_t VZEROcentralityBin(GetVZEROCentralityBin());

//     if(fFillQAHistograms) {
//         // recentering qa
//         fHistQxV0aBC->Fill(Qa2[0], VZEROcentralityBin);
//         fHistQyV0aBC->Fill(Qa2[1], VZEROcentralityBin);
//         fHistQxV0cBC->Fill(Qc2[0], VZEROcentralityBin);
//         fHistQyV0cBC->Fill(Qc2[0], VZEROcentralityBin);
//         fHistEPBC->Fill(.5*TMath::ATan2(QaY+QcY, QaX+QcX));
//     }

    Double_t Qx2amean = fMeanQ[VZEROcentralityBin][1][0];
    Double_t Qx2arms  = fWidthQ[VZEROcentralityBin][1][0];
    Double_t Qy2amean = fMeanQ[VZEROcentralityBin][1][1];
    Double_t Qy2arms  = fWidthQ[VZEROcentralityBin][1][1];

    Double_t Qx2cmean = fMeanQ[VZEROcentralityBin][0][0];
    Double_t Qx2crms  = fWidthQ[VZEROcentralityBin][0][0];
    Double_t Qy2cmean = fMeanQ[VZEROcentralityBin][0][1];
    Double_t Qy2crms  = fWidthQ[VZEROcentralityBin][0][1];	

    Double_t Qx3amean = fMeanQv3[VZEROcentralityBin][1][0];
    Double_t Qx3arms  = fWidthQv3[VZEROcentralityBin][1][0];
    Double_t Qy3amean = fMeanQv3[VZEROcentralityBin][1][1];
    Double_t Qy3arms  = fWidthQv3[VZEROcentralityBin][1][1];

    Double_t Qx3cmean = fMeanQv3[VZEROcentralityBin][0][0];
    Double_t Qx3crms  = fWidthQv3[VZEROcentralityBin][0][0];
    Double_t Qy3cmean = fMeanQv3[VZEROcentralityBin][0][1];
    Double_t Qy3crms  = fWidthQv3[VZEROcentralityBin][0][1];	

    // update the weighted q-vectors with the re-centered values
    Qa2[0] = (Qa2[0] - Qx2amean)/Qx2arms;
    Qa2[1] = (Qa2[1] - Qy2amean)/Qy2arms;
    Qc2[0] = (Qc2[0] - Qx2cmean)/Qx2crms;
    Qc2[1] = (Qc2[1] - Qy2cmean)/Qy2crms;

    Qa3[0] = (Qa3[0] - Qx3amean)/Qx3arms;
    Qa3[1] = (Qa3[1] - Qy3amean)/Qy3arms;
    Qc3[0] = (Qc3[0] - Qx3cmean)/Qx3crms;
    Qc3[1] = (Qc3[1] - Qy3cmean)/Qy3crms;

//     if(fFillQAHistograms) {
//         // recentering qa
//         fHistQxV0a->Fill(Qa2[0], VZEROcentralityBin);
//         fHistQyV0a->Fill(Qa2[1], VZEROcentralityBin);
//         fHistQxV0c->Fill(Qc2[0], VZEROcentralityBin);
//         fHistQyV0c->Fill(Qc2[0], VZEROcentralityBin);
//         fHistEP->Fill(.5*TMath::ATan2(Qa2[1]+Qc2[1], Qa2[0]+Qc2[0]));
//     }
}
//_____________________________________________________________________________
void AliHFEreducedEventCreatorESD::CalculateQvectorCombinedVZERO(Double_t Q2[2], Double_t Q3[2]) const
{
    // calculate calibrated q-vector of the combined vzeroa, vzeroc system
    // this is somewhat ugly as CalculateQvectorCombinedVZERO is called more than once per event
    // but for now it will have to do ...
    #ifdef ALIANALYSISTASKJETV2_DEBUG_FLAG_1
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif

    // first step: retrieve the q-vectors component-wise per vzero detector
    Double_t QA2[] = {-999., -999.};
    Double_t QA3[] = {-999., -999.};
    Double_t QC2[] = {-999., -999.};
    Double_t QC3[] = {-999., -999.};
    CalculateQvectorVZERO(QA2, QC2, QA3, QC3);

    // get cache index and retrieve the chi weights for this centrality
    Int_t VZEROcentralityBin(GetVZEROCentralityBin());
    Double_t chi2A(1);
    Double_t chi2C(1);
    Double_t chi3A(1);
    Double_t chi3C(1);

    switch (fWeightForVZERO) {
        case kChi : {
            chi2A = fChi2A->At(VZEROcentralityBin);
            chi2C = fChi2C->At(VZEROcentralityBin);
            chi3A = fChi3A->At(VZEROcentralityBin);
            chi3C = fChi3C->At(VZEROcentralityBin);
        } break;
        case kSigmaSquared : {
            chi2A = fSigma2A->At(VZEROcentralityBin);
            chi2C = fSigma2C->At(VZEROcentralityBin);
            chi3A = fSigma3A->At(VZEROcentralityBin);
            chi3C = fSigma3C->At(VZEROcentralityBin);
            chi2A = (chi2A > 0) ? 1./chi2A : 1.;
            chi2C = (chi2C > 0) ? 1./chi2C : 1.;
            chi3A = (chi3A > 0) ? 1./chi3A : 1.;
            chi3C = (chi3C > 0) ? 1./chi3C : 1.;
        } break;
        default : break;
    }

    // bookkkeep these guys
    //Double_t qx2a(QA2[0]), qy2a(QA2[1]), qx2c(QC2[0]), qy2c(QC2[1]);  
    // combine the vzera and vzeroc signal
    Q2[0] = chi2A*chi2A*QA2[0]+chi2C*chi2C*QC2[0];
    Q2[1] = chi2A*chi2A*QA2[1]+chi2C*chi2C*QC2[1];
    Q3[0] = chi3A*chi3A*QA3[0]+chi3C*chi3C*QC3[0];
    Q3[1] = chi3A*chi3A*QA3[1]+chi3C*chi3C*QC3[1];

    /*Double_t _chi(0), _sigma(0), _none(0);
    // if requested do the EP correlation histos
    if(fHistEPCorrelations[fInCentralitySelection]) {
        switch (fWeightForVZERO) {
            case kNone : {
                chi2A = fChi2A->At(VZEROcentralityBin);
                chi2C = fChi2C->At(VZEROcentralityBin);
                _chi = .5*TMath::ATan2(chi2A*chi2A*qy2a+chi2C*chi2C*qy2c, chi2A*chi2A*qx2a+chi2C*chi2C*qx2c);
                chi2A = fSigma2A->At(VZEROcentralityBin);
                chi2C = fSigma2C->At(VZEROcentralityBin);
                chi2A = (chi2A > 0) ? 1./chi2A : 1.;
                chi2C = (chi2C > 0) ? 1./chi2C : 1.;
                _sigma = .5*TMath::ATan2(chi2A*chi2A*qy2a+chi2C*chi2C*qy2c, chi2A*chi2A*qx2a+chi2C*chi2C*qx2c);
                fHistEPCorrelations[fInCentralitySelection]->Fill(.5*TMath::ATan2(qy2a+qy2c,qx2a+qx2c), _chi, _sigma);
            } break;
            case kChi : {
                _chi = .5*TMath::ATan2(chi2A*chi2A*qy2a+chi2C*chi2C*qy2c, chi2A*chi2A*qx2a+chi2C*chi2C*qx2c);
                chi2A = fSigma2A->At(VZEROcentralityBin);
                chi2C = fSigma2C->At(VZEROcentralityBin);
                chi2A = (chi2A > 0) ? 1./chi2A : 1.;
                chi2C = (chi2C > 0) ? 1./chi2C : 1.;
                _sigma = .5*TMath::ATan2(chi2A*chi2A*qy2a+chi2C*chi2C*qy2c, chi2A*chi2A*qx2a+chi2C*chi2C*qx2c);
                fHistEPCorrelations[fInCentralitySelection]->Fill(.5*TMath::ATan2(qy2a+qy2c,qx2a+qx2c), _chi, _sigma);
            } break;
            case kSigmaSquared : {
                _sigma = .5*TMath::ATan2(chi2A*chi2A*qy2a+chi2C*chi2C*qy2c, chi2A*chi2A*qx2a+chi2C*chi2C*qx2c);
                chi2A = fChi2A->At(VZEROcentralityBin);
                chi2C = fChi2C->At(VZEROcentralityBin);
                _chi = .5*TMath::ATan2(chi2A*chi2A*qy2a+chi2C*chi2C*qy2c, chi2A*chi2A*qx2a+chi2C*chi2C*qx2c);
                fHistEPCorrelations[fInCentralitySelection]->Fill(.5*TMath::ATan2(qy2a+qy2c,qx2a+qx2c), _chi, _sigma);
             } break;
            default : break;
        }
        _none = .5*TMath::ATan2(qy2a+qy2c,qx2a+qx2c);
        fHistEPCorrAvChi[fInCentralitySelection]->Fill(_none, _chi);
        fHistEPCorrAvSigma[fInCentralitySelection]->Fill(_none, _sigma); 
        fHistEPCorrChiSigma[fInCentralitySelection]->Fill(_chi, _sigma);
    }*/ // this part also seems not to be important for the correction
}


void AliHFEreducedEventCreatorESD::CalculateEventPlaneVZERO(Double_t vzero[2][2]) const 
{
    // get the vzero event plane (a and c separately)
    #ifdef ALIANALYSISTASKJETV2_DEBUG_FLAG_1
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
   // switch (fCollisionType) {
   //     case kPbPb10h : {
            // for 10h data, get the calibrated q-vector from the database
            Double_t QA2[] = {-999., -999.};
            Double_t QA3[] = {-999., -999.};
            Double_t QC2[] = {-999., -999.};
            Double_t QC3[] = {-999., -999.};
            CalculateQvectorVZERO(QA2, QC2, QA3, QC3);
            vzero[0][0] = .5*TMath::ATan2(QA2[1], QA2[0]);
            vzero[1][0] = .5*TMath::ATan2(QC2[1], QC2[0]);
            vzero[0][1] = (1./3.)*TMath::ATan2(QA3[1], QA3[0]);
            vzero[1][1] = (1./3.)*TMath::ATan2(QC3[1], QC3[0]);
            return;     // paranoid return
    /*    } break;
        default: {
            // by default use the ep from the event header (make sure EP selection task is enabeled!)
            Double_t a(0), b(0), c(0), d(0), e(0), f(0), g(0), h(0);
            vzero[0][0] = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 8, 2, a, b);
            vzero[1][0] = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 9, 2, c, d);
            vzero[0][1] = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 8, 3, e, f);
            vzero[1][1] = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 9, 3, g, h);
            return;
        }
    }*/
}


void AliHFEreducedEventCreatorESD::CalculateEventPlaneCombinedVZERO(Double_t* comb) const
{
    // return the combined vzero event plane
    #ifdef ALIANALYSISTASKJETV2_DEBUG_FLAG_1
        printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    #endif
 
    // define some placeholders
    Double_t Q2[] = {-999., -999.};            
    Double_t Q3[] = {-999., -999.};
      
            CalculateQvectorCombinedVZERO(Q2, Q3);
            comb[0] = .5*TMath::ATan2(Q2[1], Q2[0]);
            comb[1] = (1./3.)*TMath::ATan2(Q3[1], Q3[0]);
      
}


void AliHFEreducedEventCreatorESD::UserCreateOutputObjects(){
  //
  // Create debug tree, signal cuts and track cuts
  //
    
    // this initialization necessary for EP correction, which uses the centrality classes
    if(!fCentralityClasses) {   // classes must be defined at this point
        Double_t c[] = {0., 20., 40., 60., 80., 100.};
        fCentralityClasses = new TArrayD(sizeof(c)/sizeof(c[0]), c);
        }

  fSignalCuts = new AliHFEsignalCuts("HFEsignalCuts", "HFE MC Signal definition");
  
  fTrackCuts = new AliHFEcuts("fTrackCuts", "Basic HFE track cuts");
  fTrackCuts->CreateStandardCuts();
  // Track cuts
  fTrackCuts->SetMaxChi2perClusterTPC(1000.);   // switch off this cut (for the moment);
  fTrackCuts->SetMinNClustersTPC(fNclustersTPC);
  fTrackCuts->SetMinRatioTPCclusters(0);
  fTrackCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable); 
  fTrackCuts->SetMinNClustersTPCPID(fNclustersTPCPID);
  fTrackCuts->SetMinNClustersITS(fNclustersITS);
  // Event cuts
  fTrackCuts->SetUseMixedVertex(kTRUE);
  fTrackCuts->SetVertexRange(10.);
  fTrackCuts->Initialize();

  fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");

  fHFEevent = new AliHFEreducedEvent;
  OpenFile(1);
  fHFEtree = new TTree("HFEtree", "HFE event tree");
  fHFEtree->Branch("HFEevent", "AliHFEreducedEvent", fHFEevent,128000,0);
  PostData(1, fHFEtree);
}

void AliHFEreducedEventCreatorESD::UserExec(Option_t *){
  //
  // User Exec: Fill debug Tree
  // 

  // Get PID response
  AliPIDResponse *pid = NULL;
  AliInputEventHandler *handler = dynamic_cast<AliInputEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(handler){
    pid = handler->GetPIDResponse();
  } else {
    AliError("No Handler");
  }
  if(!pid){
    AliError("No PID response");
    return;
  }
  if(!fInputEvent) {
    AliError("No Input event");
    return;
  }

  if(fRemoveFirstEvent){
    if(fAnalysisUtils->IsFirstEventInChunk(fInputEvent)) return;
  }

  AliDebug(1, "Event Selected");

  AliESDtrack copyTrack;
  
  if(fRunNumber != fInputEvent->GetRunNumber()) { // This is used to read in the new calibration whenever the run changes
        fRunNumber = fInputEvent->GetRunNumber(); 
    // Calibration file read-in goes here later on! 
        ReadVZEROCalibration2010h();
    }

  // MC info
  Bool_t mcthere = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()) != NULL;
  if(mcthere){ 
    fTrackCuts->SetMCEvent(fMCEvent);
    fSignalCuts->SetMCEvent(fMCEvent);
  }  

  fTrackCuts->SetRecEvent(fInputEvent);

  if(!fTrackCuts->CheckEventCuts("fCutsEvRec", fInputEvent)){
    AliDebug(1, "Event rejected by the event cuts\n");
    return;
  }
  
  // reject pile up in case of pp
  AliESDEvent *event = dynamic_cast<AliESDEvent *>(fInputEvent);
  if(event) {
    TString beamtype = event->GetBeamType();
    //printf("beamtype %s\n",(const char*)beamtype);
    if (strstr(beamtype,"p-p")) {
      //printf("Reject\n");
      if(fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5)){
	      AliDebug(1, "Event flagged as pileup\n");
	      return;
      }
    }
  }

  if(!fExtraCuts){
    fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
  }
  fExtraCuts->SetRecEventInfo(fInputEvent);

  // Tag all v0s in current event
  if(fV0Tagger){
    fV0Tagger->Reset();
    fV0Tagger->TagV0Tracks(fInputEvent);
  }

  // Make Reduced Event 
  //AliHFEreducedEvent hfeevent;
  fHFEevent->~AliHFEreducedEvent();
  new(fHFEevent)AliHFEreducedEvent();

  // Get run number
  fHFEevent->SetRunNumber(fInputEvent->GetRunNumber());

  // Derive trigger 
  UInt_t trigger = fInputHandler->IsEventSelected();
  if(trigger & AliVEvent::kMB) fHFEevent->SetMBTrigger();
  if((trigger & AliVEvent::kINT7)||(trigger & AliVEvent::kINT8)) fHFEevent->SetINTTrigger();
  if(trigger & AliVEvent::kCentral) fHFEevent->SetCentralTrigger();
  if(trigger & AliVEvent::kSemiCentral) fHFEevent->SetCentralTrigger();
  if(trigger & AliVEvent::kEMCEJE) fHFEevent->SetEMCALTrigger();

  /*if(fTRDTriggerAnalysis){
    fTRDTriggerAnalysis->CalcTriggers(event);
  if(fTRDTriggerAnalysis->HasTriggeredConfirmed(AliTRDTriggerAnalysis::kHSE)) fHFEevent->SetTRDSETrigger();
  if(fTRDTriggerAnalysis->HasTriggeredConfirmed(AliTRDTriggerAnalysis::kHQU)) fHFEevent->SetTRDDQTrigger();
  }*/
  // Get Primary Vertex
  const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();
  Double_t vtx[3];
  vertex->GetXYZ(vtx);
  fHFEevent->SetVX(vtx[0]);
  fHFEevent->SetVY(vtx[1]);
  fHFEevent->SetVZ(vtx[2]);
  Int_t ncontrib(vertex->GetNContributors());
  fHFEevent->SetNContribVertex(ncontrib);
  Double_t vcov[6];
  vertex->GetCovarianceMatrix(vcov);
  fHFEevent->SetVertexResolution(TMath::Sqrt(vcov[5]));
  fHFEevent->SetVertexDispersion(static_cast<const AliESDVertex *>(vertex)->GetDispersion());
  // Get Primary Vertex from SPD
  const AliVVertex *vertexSPD = event->GetPrimaryVertexSPD();
  if(vertexSPD){
    memset(vtx, 0, sizeof(Double_t) *3);
    vertexSPD->GetXYZ(vtx);
    fHFEevent->SetVXSPD(vtx[0]);
    fHFEevent->SetVYSPD(vtx[1]);
    fHFEevent->SetVZSPD(vtx[2]);
    fHFEevent->SetNContribVertexSPD(vertexSPD->GetNContributors());
    memset(vcov, 0, sizeof(Double_t)*6);
    vertex->GetCovarianceMatrix(vcov);
    fHFEevent->SetVertexResolutionSPD(TMath::Sqrt(vcov[5]));
    fHFEevent->SetVertexDispersionSPD(static_cast<const AliESDVertex *>(vertex)->GetDispersion());
  }

  // Get centrality
  AliCentrality *hicent = fInputEvent->GetCentrality();
  fHFEevent->SetCentrality(
    hicent->GetCentralityPercentile("V0M"),
    hicent->GetCentralityPercentile("V0A"),
    hicent->GetCentralityPercentile("V0C"),
    hicent->GetCentralityPercentile("TKL"),
    hicent->GetCentralityPercentile("TRK"),
    hicent->GetCentralityPercentile("ZNA"),
    hicent->GetCentralityPercentile("ZNC"),
    hicent->GetCentralityPercentile("CL0"),
    hicent->GetCentralityPercentile("CL1"),
    hicent->GetCentralityPercentile("CND")
  );
  
  // Get Magnetic Field
  fHFEevent->SetMagneticField(fInputEvent->GetMagneticField());
  
  // Event Plane Calculations
  
  AliEventplane* vEPa = fInputEvent->GetEventplane();
  Double_t qVx, qVy;
  TVector2 *qTPC = 0x0;
  Float_t V0PlanePhi = TVector2::Phi_0_2pi(vEPa->CalculateVZEROEventPlane(fInputEvent,10,2,qVx,qVy));
  if(V0PlanePhi > TMath::Pi()) V0PlanePhi = V0PlanePhi - TMath::Pi();
  Float_t V0APlanePhi = TVector2::Phi_0_2pi(vEPa->CalculateVZEROEventPlane(fInputEvent,8,2,qVx,qVy));
  if(V0APlanePhi > TMath::Pi()) V0APlanePhi = V0APlanePhi - TMath::Pi();
  Float_t V0CPlanePhi = TVector2::Phi_0_2pi(vEPa->CalculateVZEROEventPlane(fInputEvent,9,2,qVx,qVy));
  if(V0CPlanePhi > TMath::Pi()) V0CPlanePhi = V0CPlanePhi - TMath::Pi();
  
  qTPC = vEPa->GetQVector(); 
  Double_t qx = -1.0;
  Double_t qy = -1.0;
  if(qTPC) {
    qx = qTPC->X();
    qy = qTPC->Y();
  }  
  TVector2 qVectorfortrack;
  qVectorfortrack.Set(qx,qy);
  
  Float_t TPCPlanePhi = TVector2::Phi_0_2pi(qVectorfortrack.Phi())/2.;
  
  fHFEevent->SetV0PlanePhi(V0PlanePhi);
  fHFEevent->SetV0APlanePhi(V0APlanePhi);
  fHFEevent->SetV0CPlanePhi(V0CPlanePhi);
  fHFEevent->SetTPCPlanePhi(TPCPlanePhi);
  
  
  // Find centrality class for EP corrections
  Float_t centrality = hicent->GetCentralityPercentile("V0M");
  fInCentralitySelection = -1;
      for(Int_t i(0); i < fCentralityClasses->GetSize()-1; i++) {
        if(centrality >= fCentralityClasses->At(i) && centrality <= fCentralityClasses->At(1+i)) {
            fInCentralitySelection = i;
            break;
          }
      }
  
  double epCorrArray[2];
  double epCorr;
  CalculateEventPlaneCombinedVZERO(epCorrArray);
  epCorr = epCorrArray[0];
  if(epCorr<0.) epCorr += TMath::Pi();
  Double_t vzero[2][2];
  vzero[0][0]=vzero[1][0]=vzero[0][1]=vzero[1][1]=10.2;
  CalculateEventPlaneVZERO(vzero);
  if(vzero[0][0]<0.) vzero[0][0] += TMath::Pi();
  if(vzero[1][0]<0.) vzero[1][0] += TMath::Pi();
  // Now fill histograms for the 2010 EP correction
  fHFEevent->SetV0PlanePhiCorrected(epCorr);
  fHFEevent->SetV0APlanePhiCorrected(vzero[0][0]);
  fHFEevent->SetV0CPlanePhiCorrected(vzero[1][0]);
  
  
  // End of Event Plane Calculations
  
  
  // Get VZERO Information
  AliVVZERO *vzeroinfo = fInputEvent->GetVZEROData();
  if(vzeroinfo) fHFEevent->SetV0Multiplicity(vzeroinfo->GetMTotV0A(), vzeroinfo->GetMTotV0C());

  // Get ZDC Information
  AliVZDC *zdcinfo = fInputEvent->GetZDCData();
  if(zdcinfo) fHFEevent->SetZDCEnergy(zdcinfo->GetZNAEnergy(), zdcinfo->GetZNCEnergy(), zdcinfo->GetZPAEnergy(), zdcinfo->GetZPCEnergy()); 

  // Set SPD multiplicity
  const AliMultiplicity *mult = event->GetMultiplicity();
  if(mult) fHFEevent->SetSPDMultiplicity(mult->GetNumberOfTracklets());

  // Flag Pileup
  if(fFlagPileupEvents){
    if(fAnalysisUtils->IsPileUpEvent(fInputEvent)) fHFEevent->SetPileupFlag();
  }

  //
  // Loop on MC tracks only
  //
  AliMCParticle *mctrack(NULL);
  // Monte-Carlo info
  //Int_t source(5);
  if(mcthere){
    const AliVVertex *mcvtx = fMCEvent->GetPrimaryVertex();
    fHFEevent->SetVMC(mcvtx->GetX(), mcvtx->GetY(), mcvtx->GetX());
    for(Int_t itrack = 0; itrack < fMCEvent->GetNumberOfTracks(); itrack++) {
      mctrack = (AliMCParticle *)(fMCEvent->GetTrack(itrack));
      if(!mctrack) continue;
      if(TMath::Abs(mctrack->PdgCode()) != 11) continue;  // in MC truth list store only electrons
      if(fSelectSignalOnly && !fTrackCuts->CheckParticleCuts(static_cast<UInt_t>(AliHFEcuts::kStepMCGenerated), mctrack))  continue;        
      AliHFEreducedMCParticle hfemcpart;
      if(fTrackCuts->CheckParticleCuts(static_cast<UInt_t>(AliHFEcuts::kStepMCGenerated), mctrack)) hfemcpart.SetSignal();
      // Kinematics
      hfemcpart.SetSignedPt(mctrack->Pt(), mctrack->Charge() > 0.);
      hfemcpart.SetP(mctrack->P());
      hfemcpart.SetEta(mctrack->Eta());
      hfemcpart.SetPhi(mctrack->Phi());
      hfemcpart.SetPdg(mctrack->PdgCode());
      
      // Get Production Vertex in radial direction
      hfemcpart.SetProductionVertex(mctrack->Xv(),mctrack->Yv(),mctrack->Zv());

      // Get Mother PDG code of the particle
      Int_t motherlabel = TMath::Abs(mctrack->GetMother());
      if(motherlabel >= 0 && motherlabel < fMCEvent->GetNumberOfTracks()){
        AliMCParticle *mother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(motherlabel));
        if(mother){
          hfemcpart.SetMotherPdg(mother->PdgCode());
          hfemcpart.SetMotherProductionVertex(mother->Xv(),mother->Yv(),mother->Zv());
	  Int_t grmotherlabel = TMath::Abs(mother->GetMother());
	  if(grmotherlabel >= 0 && grmotherlabel < fMCEvent->GetNumberOfTracks()){
	    AliMCParticle *grmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(grmotherlabel));
	    if(grmother){
	      hfemcpart.SetMotherPdg(grmother->PdgCode());
	    }
	  }
        }
      }
      
      // derive source
      /*
      source = 7;	      
      if(fSignalCuts->IsCharmElectron(mctrack)) source = 0;
	    else if(fSignalCuts->IsBeautyElectron(mctrack)) source = 1;
	    else if(fSignalCuts->IsGammaElectron(mctrack)) source = 2;
      else if(fSignalCuts->IsNonHFElectron(mctrack)) source = 3;
      else if(fSignalCuts->IsJpsiElectron(mctrack)) source = 4;
      else if(fSignalCuts->IsB2JpsiElectron(mctrack)) source = 5;
      else if(fSignalCuts->IsKe3Electron(mctrack)) source = 6;
	    else source = 7;
      hfemcpart.SetSource(source);
      */
      hfemcpart.SetSource(static_cast<Int_t>(fSignalCuts->GetSignalSource(mctrack)));
      Double_t mpt = -1;
      Int_t electronSource = fSignalCuts->GetMCQAObject()->GetElecSource(mctrack,kTRUE, mpt);
      hfemcpart.SetElectronSource(electronSource);
      hfemcpart.SetElectronSourcePt(mpt);
      fHFEevent->AddMCParticle(&hfemcpart);
    }
  }
  
  //
  // Loop on reconstructed tracks
  //
  TArrayI arraytrack(fInputEvent->GetNumberOfTracks());
  Int_t counterdc=0;
  
  AliESDtrack *track = 0x0;
  for(Int_t itrack = 0; itrack < fInputEvent->GetNumberOfTracks(); itrack++){
    // Run track loop
    track = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(itrack));
    if(!track) continue;
    // Cut track (Only basic track cuts)
    if(!fTrackCuts->CheckParticleCuts(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepRecKineITSTPC, track)) continue;
    //

    // Kinematics
    AliHFEreducedTrack hfetrack;
    hfetrack.SetSignedPt(track->Pt(), track->Charge() > 0);
    hfetrack.SetP(track->P());
    hfetrack.SetEta(track->Eta());
    hfetrack.SetPhi(track->Phi());
    hfetrack.SetTPCmomentum(track->GetTPCmomentum());

    // Track ID
    hfetrack.SetTrackID(track->GetID());

    // status
    ULong_t status = track->GetStatus();
    if((status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit) hfetrack.SetITSrefit();
    if((status & AliVTrack::kTPCrefit) == AliVTrack::kTPCrefit) hfetrack.SetTPCrefit();
    if((status & AliVTrack::kTOFpid) == AliVTrack::kTOFpid) hfetrack.SetTOFpid();
    //if((status & AliVTrack::kTOFmismatch) == AliVTrack::kTOFmismatch) hfetrack.SetTOFmismatch();
    if(IsTOFmismatch(track, pid)) hfetrack.SetTOFmismatch(); // New version suggested by Pietro Antonioli
    Bool_t isEMCAL(kFALSE);
    Int_t fClsId = track->GetEMCALcluster();
    if(fClsId >= 0) isEMCAL = kTRUE;
    AliDebug(2, Form("cluster ID: %d, EMCAL: %s", fClsId, isEMCAL ? "yes" : "no"));
    if(isEMCAL) hfetrack.SetEMCALpid();
    // no filter bits available for ESDs

    // fill counts of v0-identified particles
    AliPID::EParticleType myv0pid = fV0Tagger->GetV0Info(track->GetID());
    AliHFEreducedTrack::EV0PID_t v0pid = AliHFEreducedTrack::kV0undef;
    if(myv0pid == AliPID::kElectron) v0pid = AliHFEreducedTrack::kV0electron;
    else if(myv0pid == AliPID::kPion) v0pid = AliHFEreducedTrack::kV0pion;
    else if(myv0pid == AliPID::kProton) v0pid = AliHFEreducedTrack::kV0proton;
    hfetrack.SetV0PID(v0pid);

    Double_t v0prodR = fV0Tagger->GetV0ProdR(track->GetID());
    hfetrack.SetV0prodR(v0prodR);

    if(mcthere){
      // Fill Monte-Carlo Information
      Int_t label = TMath::Abs(track->GetLabel());
      if(label < fMCEvent->GetNumberOfTracks())
        mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(label));
      if(mctrack){
        AliDebug(2, "Associated MC particle found");
        if(fTrackCuts->CheckParticleCuts(static_cast<UInt_t>(AliHFEcuts::kStepMCGenerated), mctrack)) hfetrack.SetMCSignal();
        // Kinematics
        hfetrack.SetMCSignedPt(mctrack->Pt(),mctrack->Charge() > 0.);
        hfetrack.SetMCP(mctrack->P());
        hfetrack.SetMCEta(mctrack->Eta());
        hfetrack.SetMCPhi(mctrack->Phi());
        hfetrack.SetMCPDG(mctrack->PdgCode());
      
        // Get Production Vertex in radial direction
        hfetrack.SetMCProdVtx(mctrack->Xv(),mctrack->Yv(),mctrack->Zv());
      
        // Get Mother PDG code of the particle
        Int_t motherlabel = TMath::Abs(mctrack->GetMother());
        if(motherlabel >= 0 && motherlabel < fMCEvent->GetNumberOfTracks()){
          AliMCParticle *mother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(motherlabel));
          if(mother){ 
            hfetrack.SetMCMotherPdg(mother->PdgCode());
            hfetrack.SetMCMotherProdVtx(mother->Xv(),mother->Yv(),mother->Zv());
          }
        }
        
        // derive source
        /*
        source = 7;
        if(fSignalCuts->IsCharmElectron(track)) source = 0;
	      else if(fSignalCuts->IsBeautyElectron(track)) source = 1;
	      else if(fSignalCuts->IsGammaElectron(track)) source = 2;
        else if(fSignalCuts->IsNonHFElectron(track)) source = 3;
        else if(fSignalCuts->IsJpsiElectron(track)) source = 4;
        else if(fSignalCuts->IsB2JpsiElectron(track)) source = 5;
        else if(fSignalCuts->IsKe3Electron(track)) source = 6;
        else source = 7;
        hfetrack.SetMCSource(source); 
        */
        hfetrack.SetMCSource(static_cast<Int_t>(fSignalCuts->GetSignalSource(track))); 
        Double_t mcmpt = -1;
        Int_t mcelectronSource = fSignalCuts->GetMCQAObject()->GetElecSource(mctrack,kTRUE, mcmpt);
        hfetrack.SetMCElectronSource(mcelectronSource);
        hfetrack.SetMCElectronSourcePt(mcmpt);
      } else {
        AliDebug(2, "Associated MC particle not found");
      }
    }

    // HFE DCA
    Float_t dcaxy = -999.,
            dcaz = -999.;
    Double_t dcaErr, dcaxyD;
    fExtraCuts->GetImpactParameters((AliVTrack *)track,dcaxy,dcaz);
    fExtraCuts->GetHFEImpactParameters((AliVTrack *)track,dcaxyD,dcaErr);
    hfetrack.SetDCA(dcaxyD, dcaz);
    hfetrack.SetDCAerr(dcaErr);
    Double_t hfeImpactParam(-999.), hfeImpactParamResol(-999.);
    fExtraCuts->GetHFEImpactParameters((AliVTrack *)track,hfeImpactParam,hfeImpactParamResol);
    hfetrack.SetHFEImpactParam(hfeImpactParam,hfeImpactParamResol);

    // Different number of clusters definitions
    Int_t nclustersITS(track->GetITSclusters(NULL)),
          nclustersTPC(track->GetTPCNcls()),
          nclustersTPCall(track->GetTPCClusterMap().CountBits()),
          nclustersTPCshared(0);
    UChar_t nfindableTPC = track->GetTPCNclsF();
    const TBits &sharedTPC = track->GetTPCSharedMap();
    for(Int_t ibit = 0; ibit < 160; ibit++) if(sharedTPC.TestBitNumber(ibit)) nclustersTPCshared++;
    hfetrack.SetChi2PerTPCcluster(track->GetTPCchi2()/Double_t(nclustersTPC));
    hfetrack.SetITSnclusters(nclustersITS);
    hfetrack.SetTPCnclusters(nclustersTPC);
    hfetrack.SetITSchi2(track->GetITSchi2());
    hfetrack.SetITSsharedClusterMap(track->GetITSSharedMap());
    hfetrack.SetTRDnclusters(track->GetTRDncls());
    hfetrack.SetTPCnclustersPID(track->GetTPCsignalN());
    hfetrack.SetTPCcrossedRows(track->GetTPCCrossedRows());
    hfetrack.SetTPCnclustersAll(nclustersTPCall);
    hfetrack.SetTPCsharedClusters(nclustersTPCshared);
    hfetrack.SetTPCclusterRatio(nfindableTPC ? static_cast<Float_t>(nclustersTPC)/static_cast<Float_t>(nfindableTPC) : 0);
    hfetrack.SetTPCclusterRatioAll(nfindableTPC ? static_cast<Float_t>(nclustersTPCall)/static_cast<Float_t>(nfindableTPC) : 0);
    UChar_t itsPixel = track->GetITSClusterMap();
    for(int ily = 0; ily < 6; ily++) 
            if(TESTBIT(itsPixel, ily)) hfetrack.SetITScluster(ily);
   
    // TRD related quantities (Yvonne)
    Int_t nslices = track->GetNumberOfTRDslices();
    hfetrack.SetTRDntrackletsPID(track->GetTRDntrackletsPID());
    hfetrack.SetTRDnslices(nslices);
    hfetrack.SetTRDchi2(track->GetTRDchi2());
    Int_t nslicetemp=0;
    for(Int_t iplane = 0; iplane < 6; iplane++){
	    nslicetemp=0;
	    for(Int_t isl = 0; isl < nslices; isl++){
	      if(track->GetTRDntrackletsPID()>0){
		      if(track->GetTRDslice(iplane, isl)>0.001) nslicetemp++;
	      }
	    }
	    if(nslicetemp > 0) hfetrack.SetTRDstatus(iplane);
    }


    //test for kink tracks
    if(fExtraCuts->IsKinkMother(track)) hfetrack.SetIsKinkMother();
    else if(fExtraCuts->IsKinkDaughter(track)) hfetrack.SetIsKinkDaughter();
    
    // Double counted
    Int_t id(track->GetID());
    for(Int_t l=0; l < counterdc; l++){
      Int_t iTrack2 = arraytrack.At(l);
      if(iTrack2==id){
         hfetrack.SetDoubleCounted();
         break;
      }
    }
    // Add the id at this place
    arraytrack.AddAt(id,counterdc);
    counterdc++;

    // PID
    hfetrack.SetTPCdEdx(track->GetTPCsignal());
    hfetrack.SetTPCsigmaEl(pid->NumberOfSigmasTPC(track, AliPID::kElectron));
    hfetrack.SetTOFsigmaEl(pid->NumberOfSigmasTOF(track, AliPID::kElectron));
    hfetrack.SetTOFsigmaP(pid->NumberOfSigmasTOF(track, AliPID::kProton));
    if(TMath::Abs(pid->NumberOfSigmasTOF(track, AliPID::kDeuteron)) < 40.) hfetrack.SetTOFsigmaDeuteron(pid->NumberOfSigmasTOF(track, AliPID::kDeuteron));   else hfetrack.SetTOFsigmaDeuteron(100);
    hfetrack.SetTOFmismatchProbability(pid->GetTOFMismatchProbability(track));
    hfetrack.SetITSsigmaEl(pid->NumberOfSigmasITS(track, AliPID::kElectron));
    hfetrack.SetITSsigmaP(pid->NumberOfSigmasITS(track, AliPID::kProton));
    // Eta correction
    copyTrack.~AliESDtrack();
    new(&copyTrack) AliESDtrack(*track);
    if(fTPCpid->HasCentralityCorrection()) fTPCpid->ApplyCentralityCorrection(&copyTrack, static_cast<Double_t>(ncontrib),AliHFEpidObject::kESDanalysis);
    if(fTPCpid->HasEtaCorrection()) fTPCpid->ApplyEtaCorrection(&copyTrack, AliHFEpidObject::kESDanalysis);
    hfetrack.SetTPCsigmaElCorrected(pid->NumberOfSigmasTPC(&copyTrack, AliPID::kElectron));
    hfetrack.SetTPCdEdxCorrected(copyTrack.GetTPCsignal());
    if(isEMCAL){
      AliDebug(2, "Adding EMCAL PID information");
      // EMCAL cluster
      Double_t emcalEnergyOverP = -1.,
               showershape[4] = {0.,0.,0.,0.};
      hfetrack.SetEMCALSigmaEl(pid->NumberOfSigmasEMCAL(track, AliPID::kElectron, emcalEnergyOverP, &showershape[0]));
      hfetrack.SetEMCALEoverP(emcalEnergyOverP);
      hfetrack.SetEMCALShowerShape(showershape);
    }

    // Track finished, add NOW to the Event
    fHFEevent->AddTrack(&hfetrack);
  }
  
  // Fill the debug tree
  //AliInfo(Form("Number of tracks: %d\n", fHFEevent->GetNumberOfTracks()));
  //AliInfo(Form("Number of MC particles: %d\n", fHFEevent->GetNumberOfMCParticles()));
  fHFEtree->Fill();

  fEventNumber++;
  PostData(1, fHFEtree);
}

void AliHFEreducedEventCreatorESD::Terminate(Option_t *){
  //
  // Terminate
  //
  AliInfo("terminating...\n");

}

Bool_t AliHFEreducedEventCreatorESD::IsTOFmismatch(const AliVTrack *const track, const AliPIDResponse *const pid) const {
  //
  // Is TOF mismatch
  //
  Double_t probs[AliPID::kSPECIESC];
  AliPIDResponse::EDetPidStatus status = pid->ComputeTOFProbability(track, AliPID::kSPECIESC, probs);
  return status == AliPIDResponse::kDetMismatch;
}

