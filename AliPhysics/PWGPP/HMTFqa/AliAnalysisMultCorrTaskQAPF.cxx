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

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// This is a simple task designed to test AliPPVsMultUtils functionality.
// --- david.dobrigkeit.chinellato@cern.ch
// Also has been added analysis for dN/dEta and Multiplicity + Multiplicity correlations with V0 estimators+ V0 sectors
//  THIS CODE includes High mult trigger and Past-Future protection.
// ---- hector.bello.martinez@cern.ch
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+


#include <Riostream.h>
#include "TList.h"
#include <TChain.h>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
//#define LOG_NO_INFO
//#define LOG_NO_DEBUG
//#define LOG_NO_WARNING
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include <AliHeader.h>
#include "AliCentrality.h"

#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliPPVsMultUtils.h"
#include "AliAnalysisMultCorrTaskQAPF.h"

#include <AliESDVertex.h>
#include <AliMultiplicity.h>

#include <TTree.h>
#include <TDirectory.h>
#include <TBits.h>

#include "AliESDInputHandler.h"
#include <AliAnalysisFilter.h>

using std::cout;
using std::endl;

ClassImp(AliAnalysisMultCorrTaskQAPF)

AliAnalysisMultCorrTaskQAPF::AliAnalysisMultCorrTaskQAPF()
    : AliAnalysisTaskSE(), 
    fListHist(0), 
    fCuts(0),
    fTrackFilterSKE(0x0),
    fHistEventCounter(0), 
    fHistRefMult08(0), 
    fHistRefMult05(0), 
    fHistV0M(0), 
    fHistV0A(0), 
    fHistV0C(0), 
    fdNdeta(0), 
    fcorrRef05Ref08(0), 
    fcorrV0ARef08(0), 
    fcorrV0CRef08(0), 
    fcorrV0MRef08(0), 
    fHistV0Aamp(0), 
    fHistV0Camp(0), 
    fHistV0Mamp(0), 
    fcorrV0AampRef08(0), 
    fcorrV0CampRef08(0), 
    fcorrV0MampRef08(0),
    fModulesV0(0), 
    fTrackvsclust(0),
     fHistRefMult08PF(0),
     fHistRefMult05PF(0),
     fdNdetaPF(0),
     fcorrRef05Ref08PF(0),
     fHistV0AampPF(0),
     fHistV0CampPF(0),
     fHistV0MampPF(0),
     fcorrV0AampRef08PF(0),
     fcorrV0CampRef08PF(0),
     fcorrV0MampRef08PF(0),
     fModulesV0PF(0),
     fTrackvsclustPF(0),
     fHistRefMult08BG(0),
     fHistRefMult05BG(0),
     fdNdetaBG(0),
     fcorrRef05Ref08BG(0),
     fHistV0AampBG(0),
     fHistV0CampBG(0),
     fHistV0MampBG(0),
     fcorrV0AampRef08BG(0),
     fcorrV0CampRef08BG(0),
     fcorrV0MampRef08BG(0),
     fModulesV0BG(0),
     fTrackvsclustBG(0),
    fHistRefMult08HM(0),
    fHistRefMult05HM(0),
    fdNdetaHM(0),
    fcorrRef05Ref08HM(0),
    fHistV0AampHM(0),
    fHistV0CampHM(0),
    fHistV0MampHM(0),
    fcorrV0AampRef08HM(0),
    fcorrV0CampRef08HM(0),
    fcorrV0MampRef08HM(0),
    fModulesV0HM(0),
    fTrackvsclustHM(0),
    fTrackvsclustHMBG(0),
    fTrackvsclustHMPF(0),
    fPPVsMultUtils(0)
{

}

AliAnalysisMultCorrTaskQAPF::AliAnalysisMultCorrTaskQAPF(const char *name)
    : AliAnalysisTaskSE(name), 
    fListHist(0), 
    fCuts(0), 
    fTrackFilterSKE(0x0),
    fHistEventCounter(0), 
    fHistRefMult08(0), 
    fHistRefMult05(0), 
    fHistV0M(0), 
    fHistV0A(0), 
    fHistV0C(0), 
    fdNdeta(0), 
    fcorrRef05Ref08(0), 
    fcorrV0ARef08(0), 
    fcorrV0CRef08(0),
    fcorrV0MRef08(0), 
    fHistV0Aamp(0), 
    fHistV0Camp(0), 
    fHistV0Mamp(0), 
    fcorrV0AampRef08(0),
    fcorrV0CampRef08(0), 
    fcorrV0MampRef08(0),
    fModulesV0(0),
    fTrackvsclust(0),
     fHistRefMult08PF(0),
     fHistRefMult05PF(0),  
     fdNdetaPF(0),
     fcorrRef05Ref08PF(0),
     fHistV0AampPF(0),
     fHistV0CampPF(0),
     fHistV0MampPF(0),
     fcorrV0AampRef08PF(0),
     fcorrV0CampRef08PF(0),
     fcorrV0MampRef08PF(0),
     fModulesV0PF(0),
     fTrackvsclustPF(0),
     fHistRefMult08BG(0),
     fHistRefMult05BG(0),
     fdNdetaBG(0),
     fcorrRef05Ref08BG(0),
     fHistV0AampBG(0),
     fHistV0CampBG(0),
     fHistV0MampBG(0),
     fcorrV0AampRef08BG(0),
     fcorrV0CampRef08BG(0),
     fcorrV0MampRef08BG(0),
     fModulesV0BG(0),
     fTrackvsclustBG(0),
    fHistRefMult08HM(0),
    fHistRefMult05HM(0),
    fdNdetaHM(0),
    fcorrRef05Ref08HM(0),
    fHistV0AampHM(0),
    fHistV0CampHM(0),
    fHistV0MampHM(0),
    fcorrV0AampRef08HM(0),
    fcorrV0CampRef08HM(0),
    fcorrV0MampRef08HM(0),
    fModulesV0HM(0),
    fTrackvsclustHM(0),
    fTrackvsclustHMBG(0),
    fTrackvsclustHMPF(0),
    fPPVsMultUtils(0)
{
    DefineOutput(1, TList::Class()); // List of Histograms
}


AliAnalysisMultCorrTaskQAPF::~AliAnalysisMultCorrTaskQAPF()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------

    if (fListHist) {
        delete fListHist;
        fListHist = 0x0;
    }
    if (fPPVsMultUtils) {
        delete fPPVsMultUtils;
        fPPVsMultUtils = 0x0;
    }
}

//________________________________________________________________________
void AliAnalysisMultCorrTaskQAPF::UserCreateOutputObjects()
{
    //Helper
    if(! fPPVsMultUtils ) {
        fPPVsMultUtils = new AliPPVsMultUtils();
    }
    //... cuts .....
    fCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    
    inputHandler->SetNeedField();

    // Create histograms
    fListHist = new TList();
    fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner
    if(! fHistEventCounter ){
        fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",13,0,13);
        fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
        fHistEventCounter->GetXaxis()->SetBinLabel(2, "Selected by Analysis");
        fHistEventCounter->GetXaxis()->SetBinLabel(3, "SPDPILEUP"); //PileupSPDInMultBins");
        fHistEventCounter->GetXaxis()->SetBinLabel(4, "notINEL>0"); //NotINEL>0");
        fHistEventCounter->GetXaxis()->SetBinLabel(5, "NotVTXcut");//NotinVertexcut");
        fHistEventCounter->GetXaxis()->SetBinLabel(6, "SPDvsTrackINC");//Inconsistent SPD&TrackVertx");
//        fHistEventCounter->GetXaxis()->SetBinLabel(7, "Not MinimumBias");
	fHistEventCounter->GetXaxis()->SetBinLabel(8, "CINT7 MB"); //CINT7-B-NOPF-CENT");
	fHistEventCounter->GetXaxis()->SetBinLabel(9, "VOHM"); //CVHMV0M-B-");
	fHistEventCounter->GetXaxis()->SetBinLabel(10,"SPDHM"); // CVHMSH2-B-");
	fHistEventCounter->GetXaxis()->SetBinLabel(11,"VOHM+smu");// CVHMV0MMSL-B-");
	fHistEventCounter->GetXaxis()->SetBinLabel(12,"SPDHM+smu");// CVHMSH2MSL-B-"); 
	fHistEventCounter->GetXaxis()->SetBinLabel(13,"SelectedHighMult");// CVHMSH2MSL-B-"); 
        fListHist->Add(fHistEventCounter);
    }
    //Histogram Output: Event-by-Event
    if(! fHistRefMult08 ) {
        fHistRefMult08 = new TH1D( "fHistRefMult08", "Multiplicity |#eta| < 0.8 ;Ref. Mult. |#eta| < 0.8;Count",200,0,200);
        fListHist->Add(fHistRefMult08);
    }
    if(! fHistRefMult05 ) {
        fHistRefMult05 = new TH1D( "fHistRefMult05", "Multiplicity |#eta| < 0.5 ;Ref. Mult. |#eta| < 0.5;Count",200,0,200);
        fListHist->Add(fHistRefMult05);
    }
    if(! fHistV0M ) {
        fHistV0M = new TH1D( "fHistV0M", "Multiplicity V0M;V0M Percentile;Count",100,0,100);
       // fListHist->Add(fHistV0M);
    }
    if(! fHistV0A ) {
        fHistV0A = new TH1D( "fHistV0A", "Multiplicity V0A;V0A Percentile;Count",100,0,100);
       // fListHist->Add(fHistV0A);
    }
    if(! fHistV0C ) {
        fHistV0C = new TH1D( "fHistV0C", "Multiplicity V0C;V0C Percentile;Count",100,0,100);
       // fListHist->Add(fHistV0C);
    }
    if(! fHistV0Aamp ) {
        fHistV0Aamp = new TH1D( "fHistV0Aamp", "Multiplicity V0A;V0A Amplitud;Count",700,0,700);
        fListHist->Add(fHistV0Aamp);
    }
    if(! fHistV0Camp ) {
        fHistV0Camp = new TH1D( "fHistV0Camp", "Multiplicity V0C;V0C Amplitud;Count",700,0,700);
        fListHist->Add(fHistV0Camp);
    }
    if(! fHistV0Mamp ) {
        fHistV0Mamp = new TH1D( "fHistV0Mamp", "Multiplicity V0M;V0M Amplitud;Count",700,0,700);
        fListHist->Add(fHistV0Mamp);
    }
    if(!fdNdeta){
        fdNdeta= new TH1D( "fdNdeta","dN/d#eta (|#eta| < 0.8);#eta (rads);Count",100,-1,1);
        fListHist->Add(fdNdeta);
    }

    if(!fcorrRef05Ref08){
        fcorrRef05Ref08= new TH2D( "fcorrRef05Ref08"," Multiplicity Correlation (|#eta| < 0.5 and |#eta| < 0.8); Nch (|#eta| < 0.5); Nch (|#eta| < 0.8)",200,0,200,200,0,200);
        fListHist->Add(fcorrRef05Ref08);
    }

    if(!fcorrV0ARef08){
	    fcorrV0ARef08= new TH2D( "fcorrV0ARef08"," Multiplicity Correlation (V0A percentile and |#eta| < 0.8); V0A Percentile; Nch (|#eta| < 0.8)",100,0,100,200,0,200);
	    fListHist->Add(fcorrV0ARef08);
    }
    if(!fcorrV0CRef08){
	    fcorrV0CRef08= new TH2D( "fcorrV0CRef08"," Multiplicity Correlation (V0C percentile and |#eta| < 0.8); V0C Percentile; Nch (|#eta| < 0.8)",100,0,100,200,0,200);
	    fListHist->Add(fcorrV0CRef08);
    }
    if(!fcorrV0AampRef08){
	    fcorrV0AampRef08= new TH2D( "fcorrV0AampRef08"," Multiplicity Correlation (V0A and |#eta| < 0.8); V0A Amplitud; Nch (|#eta| < 0.8)",700,0,700,200,0,200);
	    fListHist->Add(fcorrV0AampRef08);
    }
    if(!fcorrV0CampRef08){
	    fcorrV0CampRef08= new TH2D( "fcorrV0CampRef08"," Multiplicity Correlation (V0C and |#eta| < 0.8); V0C Amplitud; Nch (|#eta| < 0.8)",700,0,700,200,0,200);
	    fListHist->Add(fcorrV0CampRef08);
    }
    if(!fcorrV0MRef08){
	    fcorrV0MRef08= new TH2D( "fcorrV0MRef08"," Multiplicity Correlation (V0M and |#eta| < 0.8); V0M Percentile; Nch (|#eta| < 0.8)",100,0,100,200,0,200);
	    fListHist->Add(fcorrV0MRef08);
    }
    if(!fcorrV0MampRef08){
	    fcorrV0MampRef08= new TH2D( "fcorrV0MampRef08"," Multiplicity Correlation (V0M and |#eta| < 0.8); V0M Amplitud; Nch (|#eta| < 0.8)",700,0,700,200,0,200);
	    fListHist->Add(fcorrV0MampRef08);
    }
    if(!fModulesV0){
	    fModulesV0= new TH2D( "fModulesV0"," Multiplicity vs cell; V0 Sector; counts ",64,0,64,1000,0,1000);
	    fListHist->Add(fModulesV0);
    }
    if(!fTrackvsclust){
            fTrackvsclust= new TH2D("fTrackvsclust"," Tracklets vs Clust; N_{tracklets}; N_{SPDclusters} ",100,0,100,1000,0,1000);
	    fListHist->Add(fTrackvsclust);
    } 


    if(! fHistRefMult08PF ) {
     fHistRefMult08PF = new TH1D( "fHistRefMult08PF", "Multiplicity |#eta| < 0.8 ;Ref. Mult. |#eta| < 0.8;Count",200,0,200);
     fListHist->Add(fHistRefMult08PF);
    }
    if(! fHistRefMult05PF ) {
     fHistRefMult05PF = new TH1D( "fHistRefMult05PF", "Multiplicity |#eta| < 0.5 ;Ref. Mult. |#eta| < 0.5;Count",200,0,200);
     fListHist->Add(fHistRefMult05PF);
    }
    if(! fHistV0AampPF ) {
     fHistV0AampPF = new TH1D( "fHistV0AampPF", "Multiplicity V0A;V0A Amplitud;Count",700,0,700);
     fListHist->Add(fHistV0AampPF);
    }
    if(! fHistV0CampPF ) {
     fHistV0CampPF = new TH1D( "fHistV0CampPF", "Multiplicity V0C;V0C Amplitud;Count",700,0,700);
     fListHist->Add(fHistV0CampPF);
    }
    if(! fHistV0MampPF ) {
     fHistV0MampPF = new TH1D( "fHistV0MampPF", "Multiplicity V0M;V0M Amplitud;Count",700,0,700);
     fListHist->Add(fHistV0MampPF);
    }
    if(!fdNdetaPF){
     fdNdetaPF= new TH1D( "fdNdetaPF","dN/d#eta (|#eta| < 0.8);#eta (rads);Count",100,-1,1);
     fListHist->Add(fdNdetaPF);
    } 
    if(!fcorrRef05Ref08PF){
     fcorrRef05Ref08PF= new TH2D( "fcorrRef05Ref08PF"," Multiplicity Correlation (|#eta| < 0.5 and |#eta| < 0.8); Nch (|#eta| < 0.5); Nch (|#eta| <0.8)",200,0,200,200,0,200);
     fListHist->Add(fcorrRef05Ref08PF);
    }
    if(!fcorrV0AampRef08PF){
      fcorrV0AampRef08PF= new TH2D( "fcorrV0AampRef08PF"," Multiplicity Correlation (V0A and |#eta| < 0.8); V0A Amplitud; Nch (|#eta| < 0.8)",700,0,700,200,0,200);
      fListHist->Add(fcorrV0AampRef08PF);
    }
    if(!fcorrV0CampRef08PF){
      fcorrV0CampRef08PF= new TH2D( "fcorrV0CampRef08PF"," Multiplicity Correlation (V0C and |#eta| < 0.8); V0C Amplitud; Nch (|#eta| < 0.8)",700,0,700,200,0,200);
      fListHist->Add(fcorrV0CampRef08PF);
    }
    if(!fcorrV0MampRef08PF){
      fcorrV0MampRef08PF= new TH2D( "fcorrV0MampRef08PF"," Multiplicity Correlation (V0M and |#eta| < 0.8); V0M Amplitud; Nch (|#eta| < 0.8)",700,0,700,200,0,200);
      fListHist->Add(fcorrV0MampRef08PF);
    }
    if(!fModulesV0PF){
      fModulesV0PF= new TH2D( "fModulesV0PF"," Multiplicity vs cell; V0 Sector; counts ",64,0,64,1000,0,1000);
      fListHist->Add(fModulesV0PF);
    }
    if(!fTrackvsclustPF){
      fTrackvsclustPF= new TH2D("fTrackvsclustPF"," Tracklets vs Clust with PF and BGrej; N_{tracklets}; N_{SPDclusters} ",100,0,100,1000,0,1000);
      fListHist->Add(fTrackvsclustPF);
    }

    if(! fHistRefMult08BG ) {
       fHistRefMult08BG = new TH1D( "fHistRefMult08BG", "Multiplicity |#eta| < 0.8 ;Ref. Mult. |#eta| < 0.8;Count",200,0,200);
       fListHist->Add(fHistRefMult08BG);
    }
    if(! fHistRefMult05BG ) {
       fHistRefMult05BG = new TH1D( "fHistRefMult05BG", "Multiplicity |#eta| < 0.5 ;Ref. Mult. |#eta| < 0.5;Count",200,0,200);
       fListHist->Add(fHistRefMult05BG);
    }
    if(! fHistV0AampBG ) {
       fHistV0AampBG = new TH1D( "fHistV0AampBG", "Multiplicity V0A;V0A Amplitud;Count",700,0,700);
       fListHist->Add(fHistV0AampBG);
    } 
    if(! fHistV0CampBG ) {
       fHistV0CampBG = new TH1D( "fHistV0CampBG", "Multiplicity V0C;V0C Amplitud;Count",700,0,700);
       fListHist->Add(fHistV0CampBG);
    }
    if(! fHistV0MampBG ) {
       fHistV0MampBG = new TH1D( "fHistV0MampBG", "Multiplicity V0M;V0M Amplitud;Count",700,0,700);
       fListHist->Add(fHistV0MampBG);
    }
    if(!fdNdetaBG){
       fdNdetaBG= new TH1D( "fdNdetaBG","dN/d#eta (|#eta| < 0.8);#eta (rads);Count",100,-1,1);
       fListHist->Add(fdNdetaBG);
    }
    if(!fcorrRef05Ref08BG){
       fcorrRef05Ref08BG= new TH2D( "fcorrRef05Ref08BG"," Multiplicity Correlation (|#eta| < 0.5 and |#eta| < 0.8); Nch (|#eta| < 0.5); Nch (|#eta| <0.8)",200,0,200,200,0,200);
       fListHist->Add(fcorrRef05Ref08BG); 
    }
    if(!fcorrV0AampRef08BG){
	 fcorrV0AampRef08BG= new TH2D( "fcorrV0AampRef08BG"," Multiplicity Correlation (V0A and |#eta| < 0.8); V0A Amplitud; Nch (|#eta| < 0.8)",700,0,700,200,0,200);
	 fListHist->Add(fcorrV0AampRef08BG);
    }
    if(!fcorrV0CampRef08BG){
	 fcorrV0CampRef08BG= new TH2D( "fcorrV0CampRef08BG"," Multiplicity Correlation (V0C and |#eta| < 0.8); V0C Amplitud; Nch (|#eta| < 0.8)",700,0,700,200,0,200);
	 fListHist->Add(fcorrV0CampRef08BG);
    }
    if(!fcorrV0MampRef08BG){
	 fcorrV0MampRef08BG= new TH2D( "fcorrV0MampRef08BG"," Multiplicity Correlation (V0M and |#eta| < 0.8); V0M Amplitud; Nch (|#eta| < 0.8)",700,0,700,200,0,200);
	 fListHist->Add(fcorrV0MampRef08BG);
    }
    if(!fModulesV0BG){
	 fModulesV0BG= new TH2D( "fModulesV0BG"," Multiplicity vs cell; V0 Sector; counts ",64,0,64,1000,0,1000);
	 fListHist->Add(fModulesV0BG);
    }
    if(!fTrackvsclustBG){
      fTrackvsclustBG= new TH2D("fTrackvsclustBG"," Tracklets vs Clust with BGrej; N_{tracklets}; N_{SPDclusters} ",100,0,100,1000,0,1000);
      fListHist->Add(fTrackvsclustBG);
    }

    if(! fHistRefMult08HM ) {
	fHistRefMult08HM = new TH1D( "fHistRefMult08HM", "Multiplicity |#eta| < 0.8 ;Ref. Mult. |#eta| < 0.8;Count",200,0,200);
	fListHist->Add(fHistRefMult08HM);
    }
    if(! fHistRefMult05HM ) {
        fHistRefMult05HM = new TH1D( "fHistRefMult05HM", "Multiplicity |#eta| < 0.5 ;Ref. Mult. |#eta| < 0.5;Count",200,0,200);
	fListHist->Add(fHistRefMult05HM);
    }
    if(! fHistV0AampHM ) {
        fHistV0AampHM = new TH1D( "fHistV0AampHM", "Multiplicity V0A;V0A Amplitud;Count",700,0,700);
	fListHist->Add(fHistV0AampHM);
    }
    if(! fHistV0CampHM ) {
	fHistV0CampHM = new TH1D( "fHistV0CampHM", "Multiplicity V0C;V0C Amplitud;Count",700,0,700);
	fListHist->Add(fHistV0CampHM);
    }
    if(! fHistV0MampHM ) {
	fHistV0MampHM = new TH1D( "fHistV0MampHM", "Multiplicity V0M;V0M Amplitud;Count",700,0,700);
	fListHist->Add(fHistV0MampHM);
    }
    if(!fdNdetaHM){
	fdNdetaHM= new TH1D( "fdNdetaHM","dN/d#eta (|#eta| < 0.8);#eta (rads);Count",100,-1,1);
	fListHist->Add(fdNdetaHM);
    }
    if(!fcorrRef05Ref08HM){
	fcorrRef05Ref08HM= new TH2D( "fcorrRef05Ref08HM"," Multiplicity Correlation (|#eta| < 0.5 and |#eta| < 0.8); Nch (|#eta| < 0.5); Nch (|#eta| <0.8)",200,0,200,200,0,200);
	fListHist->Add(fcorrRef05Ref08HM);
    }
    if(!fcorrV0AampRef08HM){
	fcorrV0AampRef08HM= new TH2D( "fcorrV0AampRef08HM"," Multiplicity Correlation (V0A and |#eta| < 0.8); V0A Amplitud; Nch (|#eta| < 0.8)",700,0,700,200,0,200);
	 fListHist->Add(fcorrV0AampRef08HM);
    }
    if(!fcorrV0CampRef08HM){
	fcorrV0CampRef08HM= new TH2D( "fcorrV0CampRef08HM"," Multiplicity Correlation (V0C and |#eta| < 0.8); V0C Amplitud; Nch (|#eta| < 0.8)",700,0,700,200,0,200);
        fListHist->Add(fcorrV0CampRef08HM);
    }
    if(!fcorrV0MampRef08HM){
	fcorrV0MampRef08HM= new TH2D( "fcorrV0MampRef08HM"," Multiplicity Correlation (V0M and |#eta| < 0.8); V0M Amplitud; Nch (|#eta| < 0.8)",700,0,700,200,0,200);
        fListHist->Add(fcorrV0MampRef08HM);
    }
    if(!fModulesV0HM){
	fModulesV0HM= new TH2D( "fModulesV0HM"," Multiplicity vs cell; V0 Sector; counts ",64,0,64,1000,0,1000);
        fListHist->Add(fModulesV0HM);
    }
    if(!fTrackvsclustHM){
       fTrackvsclustHM= new TH2D("fTrackvsclustHM"," Tracklets vs Clust HM withBG+PF; N_{tracklets}; N_{SPDclusters} ",100,0,100,1000,0,1000);    
       fListHist->Add(fTrackvsclustHM);
    }
    if(!fTrackvsclustHMBG){
      fTrackvsclustHMBG= new TH2D("fTrackvsclustHMBG"," Tracklets vs Clust HM withBG+PF; N_{tracklets}; N_{SPDclusters} ",100,0,100,1000,0,1000);
      fListHist->Add(fTrackvsclustHMBG);
    }
    if(!fTrackvsclustHMPF){
      fTrackvsclustHMPF= new TH2D("fTrackvsclustHMPF"," Tracklets vs Clust HM withBG+PF; N_{tracklets}; N_{SPDclusters} ",100,0,100,1000,0,1000);
      fListHist->Add(fTrackvsclustHMPF);
    }

    PostData(1, fListHist);
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisMultCorrTaskQAPF::UserExec(Option_t *)
{
	// Main loop
	// Called for each event
	AliESDEvent *lESDevent = 0x0;

	lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
	if (!lESDevent) {
		AliWarning("ERROR: lESDevent not available \n");
		return;
	}
	//------------------------------------------------
	// Selection Investigation with AliPPVsMultUtils
	//------------------------------------------------
	fHistEventCounter->Fill(0.5);

	if( AliPPVsMultUtils::IsNotPileupSPDInMultBins(lESDevent) == kFALSE){fHistEventCounter->Fill(2.5);}
	if( AliPPVsMultUtils::IsINELgtZERO( lESDevent ) == kFALSE){fHistEventCounter->Fill(3.5);}
	if( AliPPVsMultUtils::IsAcceptedVertexPosition( lESDevent ) == kFALSE){fHistEventCounter->Fill(4.5);}
	if( AliPPVsMultUtils::HasNoInconsistentSPDandTrackVertices( lESDevent ) == kFALSE){fHistEventCounter->Fill(5.5);}
//	if( AliPPVsMultUtils::IsMinimumBias( lESDevent ) == kFALSE){fHistEventCounter->Fill(6.5);}

        TString classes = lESDevent->GetFiredTriggerClasses();
	if (classes.Contains("CINT7-B-NOPF-CENT")) { fHistEventCounter->Fill(7.5); /* minimum bias */ }
	if (classes.Contains("CVHMV0M-B-")) { fHistEventCounter->Fill(8.5); /* V0 HM */ }
	if (classes.Contains("CVHMSH2-B-")) { fHistEventCounter->Fill(9.5); /* SPD HM */ }
	if (classes.Contains("CVHMV0MMSL-B-")) { fHistEventCounter->Fill(10.5); /* V0 HM + single muon */ }
	if (classes.Contains("CVHMSH2MSL-B-")) { fHistEventCounter->Fill(11.5);  /* SPD HM + single muon */} 

        // ----- For V0 High Mult for LHC15i 

   if (classes.Contains("CINT7-B-NOPF-CENT")){ // MINIMUM BIAS//
	

	fHistEventCounter->Fill(1.5);

	//Reference Multiplicity   
	fHistRefMult08->Fill(AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );
	fHistRefMult05->Fill(AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.5) );
	//V0 Percentiles
	fHistV0M->Fill( fPPVsMultUtils->GetMultiplicityPercentile( lESDevent, "V0M" ) );
	fHistV0A->Fill( fPPVsMultUtils->GetMultiplicityPercentile( lESDevent, "V0A" ) );
	fHistV0C->Fill( fPPVsMultUtils->GetMultiplicityPercentile( lESDevent, "V0C" ) );
	LoopESD(lESDevent);
	AliVVZERO* esdV0 = lESDevent->GetVZEROData();

	Float_t multV0A= esdV0->GetMTotV0A();
	Float_t multV0C= esdV0->GetMTotV0C();
	Float_t multV0M= multV0A+multV0C;
	
        fHistV0Aamp->Fill( esdV0->GetMTotV0A() );
	fHistV0Camp->Fill( esdV0->GetMTotV0C() );


	fHistV0Mamp->Fill( multV0M );
	//Correlations  
	fcorrRef05Ref08->Fill(AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.5),AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );
	fcorrV0ARef08->Fill( fPPVsMultUtils->GetMultiplicityPercentile( lESDevent, "V0A" ), AliPPVsMultUtils::GetStandardReferenceMultiplicity( lESDevent ));
	fcorrV0CRef08->Fill( fPPVsMultUtils->GetMultiplicityPercentile( lESDevent, "V0C" ), AliPPVsMultUtils::GetStandardReferenceMultiplicity( lESDevent ));
	fcorrV0AampRef08->Fill( esdV0->GetMTotV0A(), AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );
	fcorrV0CampRef08->Fill( esdV0->GetMTotV0C(), AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );
	fcorrV0MampRef08->Fill( esdV0->GetMTotV0A()+esdV0->GetMTotV0C(), AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );
	fcorrV0MRef08->Fill( fPPVsMultUtils->GetMultiplicityPercentile( lESDevent, "V0M" ), AliPPVsMultUtils::GetStandardReferenceMultiplicity( lESDevent ));



	for(Int_t ich =0; ich < 64; ich++){
		fModulesV0->Fill(ich,esdV0->GetMultiplicity(ich));
	}
          Int_t nClustersLayer0 = lESDevent->GetNumberOfITSClusters(0);
	  Int_t nClustersLayer1 = lESDevent->GetNumberOfITSClusters(1);
	  Int_t nTracklets      = lESDevent->GetMultiplicity()->GetNumberOfTracklets();
        fTrackvsclust->Fill(nTracklets, nClustersLayer0+nClustersLayer1);
       //---- BG rejection 
        if ( nClustersLayer0+ nClustersLayer1 < 65+4*nTracklets ) {
            fHistRefMult08BG->Fill(AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );
	    fHistRefMult05BG->Fill(AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.5) );
	    LoopESDBG(lESDevent);
	    AliVVZERO* esdV0BG = lESDevent->GetVZEROData();
            Float_t multV0ABG= esdV0BG->GetMTotV0A();
	    Float_t multV0CBG= esdV0BG->GetMTotV0C();
	    Float_t multV0MBG= multV0ABG+multV0CBG;
	    fHistV0AampBG->Fill( esdV0BG->GetMTotV0A() );
	    fHistV0CampBG->Fill( esdV0BG->GetMTotV0C() );
            fHistV0MampBG->Fill( multV0MBG );
	    fcorrRef05Ref08BG->Fill(AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.5), AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );
            fcorrV0AampRef08BG->Fill( esdV0BG->GetMTotV0A(), AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );
	    fcorrV0CampRef08BG->Fill( esdV0BG->GetMTotV0C(), AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );
            fcorrV0MampRef08BG->Fill( esdV0BG->GetMTotV0A()+esdV0BG->GetMTotV0C(), AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );
            for(Int_t ich =0; ich < 64; ich++){
	        fModulesV0BG->Fill(ich,esdV0BG->GetMultiplicity(ich));
	    }
	    fTrackvsclustBG->Fill(nTracklets, nClustersLayer0+nClustersLayer1);

         //----- Past Future protection ---------
            TBits fIR1 =  lESDevent->GetHeader()->GetIRInt1InteractionMap();
       	    Bool_t isOutOfBunchPileup11BC = 0;
	    for (Int_t i=1;i<=11;i++){ isOutOfBunchPileup11BC|=fIR1.TestBitNumber(90-i);
	        printf("isoOut=%d \n", isOutOfBunchPileup11BC );
	    }
            if( isOutOfBunchPileup11BC==1) return;


            fHistRefMult08PF->Fill(AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );
            fHistRefMult05PF->Fill(AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.5) );
            LoopESDPF(lESDevent);
            AliVVZERO* esdV0PF = lESDevent->GetVZEROData();
            Float_t multV0APF= esdV0PF->GetMTotV0A();
            Float_t multV0CPF= esdV0PF->GetMTotV0C();
            Float_t multV0MPF= multV0APF+multV0CPF;
            fHistV0AampPF->Fill( esdV0PF->GetMTotV0A() );
            fHistV0CampPF->Fill( esdV0PF->GetMTotV0C() );
            fHistV0MampPF->Fill( multV0MPF );
            fcorrRef05Ref08PF->Fill(AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.5), AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );
            fcorrV0AampRef08PF->Fill( esdV0PF->GetMTotV0A(), AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );   
            fcorrV0CampRef08PF->Fill( esdV0PF->GetMTotV0C(), AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );
            fcorrV0MampRef08PF->Fill( esdV0PF->GetMTotV0A()+esdV0PF->GetMTotV0C(), AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );
            for(Int_t ich =0; ich < 64; ich++){
               fModulesV0PF->Fill(ich,esdV0PF->GetMultiplicity(ich));
            }
            Int_t nClustersLayer0PF = lESDevent->GetNumberOfITSClusters(0);
            Int_t nClustersLayer1PF = lESDevent->GetNumberOfITSClusters(1);
	    Int_t nTrackletsPF      = lESDevent->GetMultiplicity()->GetNumberOfTracklets();
	    fTrackvsclustPF->Fill(nTrackletsPF, nClustersLayer0PF+nClustersLayer1PF);

        }

   }
	//----  trigger for V0 High Mult

   if (classes.Contains("CVHMV0M-B-")){ // High Mult //


        Int_t nClustersLayer0HM = lESDevent->GetNumberOfITSClusters(0);
        Int_t nClustersLayer1HM = lESDevent->GetNumberOfITSClusters(1);
        Int_t nTrackletsHM      = lESDevent->GetMultiplicity()->GetNumberOfTracklets();
	 fTrackvsclustHM->Fill(nTrackletsHM, nClustersLayer0HM+nClustersLayer1HM);
        //----- Background rejection ----- 
        if ( nClustersLayer0HM+ nClustersLayer1HM < 65+4*nTrackletsHM ) {
           fTrackvsclustHMBG->Fill(nTrackletsHM, nClustersLayer0HM+nClustersLayer1HM);
	//-----Past Future Protection-------
       	   TBits fIR1 =  lESDevent->GetHeader()->GetIRInt1InteractionMap();
	          Bool_t isOutOfBunchPileup11BC = 0;
	          for (Int_t i=1;i<=11;i++){ isOutOfBunchPileup11BC|=fIR1.TestBitNumber(90-i);
	            printf("isoOut=%d \n", isOutOfBunchPileup11BC );
                  }
           if( isOutOfBunchPileup11BC==1) return;

	   fHistEventCounter->Fill(12.5);
	   fHistRefMult08HM->Fill(AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );
 	   fHistRefMult05HM->Fill(AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.5) );
	   LoopESDHM(lESDevent);
	   AliVVZERO* esdV0HM = lESDevent->GetVZEROData();
           Float_t multV0AHM= esdV0HM->GetMTotV0A();
           Float_t multV0CHM= esdV0HM->GetMTotV0C();
           Float_t multV0MHM= multV0AHM+multV0CHM;
	   fHistV0AampHM->Fill( esdV0HM->GetMTotV0A() );
           fHistV0CampHM->Fill( esdV0HM->GetMTotV0C() );
           fHistV0MampHM->Fill( multV0MHM );
           fcorrRef05Ref08HM->Fill(AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.5), AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8));
	   fcorrV0AampRef08HM->Fill( esdV0HM->GetMTotV0A(), AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );   
	   fcorrV0CampRef08HM->Fill( esdV0HM->GetMTotV0C(), AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8) );
	   fcorrV0MampRef08HM->Fill( esdV0HM->GetMTotV0A()+esdV0HM->GetMTotV0C(), AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8));
	   for(Int_t ich =0; ich < 64; ich++){
	          fModulesV0HM->Fill(ich,esdV0HM->GetMultiplicity(ich));
	         }
	   fTrackvsclustHMPF->Fill(nTrackletsHM, nClustersLayer0HM+nClustersLayer1HM);
        } 
   }
   else{
       PostData(1, fListHist);// Event isn't selected, post output data, done here
       return;
   }


	PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisMultCorrTaskQAPF::Terminate(Option_t *)
{
	// Draw result to the screen
	/*
	   TList *cRetrievedList = 0x0;
	   cRetrievedList = (TList*)GetOutputData(1);
	   if(!cRetrievedList) {
	   Printf("ERROR - AliAnalysisMultCorrTaskQA : ouput data container list not available\n");
	   return;
	   }

	   fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
	   if (!fHistEventCounter) {
	   Printf("ERROR - AliAnalysisMultCorrTaskQA : fHistEventCounter not available");
	   return;
	   }

	   TCanvas *canCheck = new TCanvas("AliAnalysisMultCorrTaskQA","Control Histo",10,10,510,510);
	   canCheck->cd(1)->SetLogy();

	   fHistEventCounter->SetMarkerStyle(22);
	   fHistEventCounter->DrawCopy("E");*/
}

//----------------------------------------------
void AliAnalysisMultCorrTaskQAPF::LoopESD( AliESDEvent *lESDevent)
{

	Int_t ntracks = lESDevent->GetNumberOfTracks();
	if(ntracks==0){
		return;
	}
	Int_t nesdtracks=0;
	TObjArray* acceptedtracks = new TObjArray();
	for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
		AliESDtrack *track = (AliESDtrack *)lESDevent->GetTrack(iTracks);

		if (!track) {
			Error("UserExec", "Could not receive track %d", iTracks);
			continue;
		}

		if(!fCuts->AcceptTrack(track))continue;

		nesdtracks++;
		fdNdeta->Fill(track->Eta());
		acceptedtracks->Add(track);

	} 
	delete acceptedtracks;

}

void AliAnalysisMultCorrTaskQAPF::LoopESDPF( AliESDEvent *lESDevent)
{
         Int_t ntracks = lESDevent->GetNumberOfTracks();
         if(ntracks==0){
                 return;
         }
         Int_t nesdtracks=0;
	 TObjArray* acceptedtracksPF = new TObjArray();
         for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
                 AliESDtrack *track = (AliESDtrack *)lESDevent->GetTrack(iTracks); 
                if (!track) {
                         Error("UserExec", "Could not receive track %d", iTracks);
                         continue;
                 }
                if(!fCuts->AcceptTrack(track)) continue;
                 nesdtracks++;
                 fdNdetaPF->Fill(track->Eta());
                 acceptedtracksPF->Add(track);
	}
	delete acceptedtracksPF;
}

 void AliAnalysisMultCorrTaskQAPF::LoopESDBG( AliESDEvent *lESDevent)
 {
        Int_t ntracks = lESDevent->GetNumberOfTracks();
        if(ntracks==0){
                   return;
        }
        Int_t nesdtracks=0;
        TObjArray* acceptedtracksBG = new TObjArray();
        for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
             AliESDtrack *track = (AliESDtrack *)lESDevent->GetTrack(iTracks);
             if (!track) {
                    Error("UserExec", "Could not receive track %d", iTracks);
                    continue;
             }
            if(!fCuts->AcceptTrack(track)) continue;
                nesdtracks++;
                fdNdetaBG->Fill(track->Eta());
		acceptedtracksBG->Add(track);
        }
        delete acceptedtracksBG;
}


void AliAnalysisMultCorrTaskQAPF::LoopESDHM( AliESDEvent *lESDevent)
{
          Int_t ntracks = lESDevent->GetNumberOfTracks();
          if(ntracks==0){ 
                  return;
          }       
          Int_t nesdtracks=0;
          TObjArray* acceptedtracksHM = new TObjArray();
          for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
                  AliESDtrack *track = (AliESDtrack *)lESDevent->GetTrack(iTracks);
                 if (!track) {
			Error("UserExec", "Could not receive track %d", iTracks);
                        continue;
                  }       
                 if(!fCuts->AcceptTrack(track))continue;
                  nesdtracks++;
                  fdNdetaHM->Fill(track->Eta());
                  acceptedtracksHM->Add(track);
      }        
	delete acceptedtracksHM;
}       

// - david.dobrigkeit.chinellato@cern.ch
// --hector.bello.martinez@cern.ch
