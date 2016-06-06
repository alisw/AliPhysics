/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

// iaaAnalysis main class
// used in local and grid execution

#include <TH1D.h>
#include "AliJIaaAnalysis.h"

#include <TClonesArray.h>

#include "AliJIaaHistos.h"
#include "AliJIaaCorrelations.h"
#include "../AliJTrackCounter.h"
#include "../AliJCard.h"
#include "../AliJEventPool.h"
#include "../AliJDataManager.h"
#include "../AliJEventHeader.h"
#include "../AliJRunHeader.h"
#include "../AliJTrack.h"
#include "../AliJPhoton.h"
#include "../AliJMCTrack.h"
#include "../AliJConst.h"
#include "../AliJAcceptanceCorrection.h"
#include "../AliJEfficiency.h"
#include <iostream>

ClassImp(AliJIaaAnalysis)

AliJIaaAnalysis::AliJIaaAnalysis() :
    TObject(),
    fExecLocal(kTRUE),
    fFirstEvent(kTRUE),
    fjtrigg((particleType)-100),
    fjassoc((particleType)-100),
    fcard(0),
    finputFile(0),
    fInclusiveFile(""),
    fevt(0),
    fhistos(0),
    fcorrelations(0),
    fAcceptanceCorrection(0x0),
    fassocPool(0),
    fphotonList(0),
    fchargedHadronList(0),
    fpizeroList(0),
    ftriggList(0),
    fassocList(0),
    finputList(0),
    fdmg(0),
    feventHeader(0),
    frunHeader(0),
    fcent(0),
    fbTriggCorrel(0),
    fbLPCorrel(0),
    fMinimumPt(0),
    fEventBC(0),
    fEfficiency(0),
    fRunTable(0),
    fHadronSelectionCut(0)
{
    // constructor
}

AliJIaaAnalysis::AliJIaaAnalysis(Bool_t execLocal) :
    TObject(),
    fExecLocal(execLocal),
    fFirstEvent(kTRUE),
    fjtrigg((particleType)-100),
    fjassoc((particleType)-100),
    fcard(0),
    finputFile(0),
    fInclusiveFile(""),
    fevt(0),
    fhistos(0),
    fcorrelations(0),
    fAcceptanceCorrection(0x0),
    fassocPool(0),
    fphotonList(0),
    fchargedHadronList(0),
    fpizeroList(0),
    ftriggList(0),
    fassocList(0),
    finputList(0),
    fdmg(0),
    feventHeader(0),
    frunHeader(0),
    fcent(0),
    fbTriggCorrel(0),
    fbLPCorrel(0),
    fMinimumPt(0),
    fEventBC(0),
    fEfficiency(0),
    fRunTable(0),
    fHadronSelectionCut(0)
{
    // constructor
}

AliJIaaAnalysis::~AliJIaaAnalysis(){
    // destructor
}

AliJIaaAnalysis::AliJIaaAnalysis(const AliJIaaAnalysis& obj) :
    TObject(),
    fExecLocal(obj.fExecLocal),
    fFirstEvent(obj.fFirstEvent),
    fjtrigg(obj.fjtrigg),
    fjassoc(obj.fjassoc),
    fcard(obj.fcard),
    finputFile(obj.finputFile), 
    fInclusiveFile(obj.fInclusiveFile),
    fevt(obj.fevt),
    fhistos(obj.fhistos), 
    fcorrelations(obj.fcorrelations),
    fAcceptanceCorrection(obj.fAcceptanceCorrection),
    fassocPool(obj.fassocPool),  
    fphotonList(obj.fphotonList),  
    fchargedHadronList(obj.fchargedHadronList), 
    fpizeroList(obj.fpizeroList), 
    ftriggList(obj.ftriggList),  
    fassocList(obj.fassocList), 
    finputList(obj.finputList), 
    fdmg(obj.fdmg), 
    feventHeader(obj.feventHeader), 
    frunHeader(obj.frunHeader), 
    fcent(obj.fcent), 
    fbTriggCorrel(obj.fbTriggCorrel), 
    fbLPCorrel(obj.fbLPCorrel), 
    fMinimumPt(obj.fMinimumPt),
    fEventBC(obj.fEventBC),
    fEfficiency(obj.fEfficiency),
    fRunTable(obj.fRunTable),
    fHadronSelectionCut(obj.fHadronSelectionCut)
{
    // copy constructor
    JUNUSED(obj);
}

AliJIaaAnalysis& AliJIaaAnalysis::operator=(const AliJIaaAnalysis& obj){
    // copy constructor
    JUNUSED(obj);
    return *this;
}

void AliJIaaAnalysis::Initialize() const{
    // init

}

void AliJIaaAnalysis::UserCreateOutputObjects(){
    // local init

    cout << "jcorran user create output objects ----------------" << endl;

    fHadronSelectionCut =int ( fcard->Get("HadronSelectionCut"));

    // Initialize the histograms needed to store the output
    fhistos = new AliJIaaHistos( fcard );

    if(fcard->Get("QualityControlLevel")>1) fhistos->Set2DHistoCreate(true);
    if(fcard->Get("QualityControlLevel")>0) fhistos->SetAcceptanceCorrectionQA(true);

    fhistos->CreateEventTrackHistos();
    fhistos->CreateCorrelationHistos();

    fhistos->fHMG->Print();

    fEventBC = (Int_t)(fcard->Get( "eventBC" ));

    // Create a class for acceptance correction
    fAcceptanceCorrection = new AliJAcceptanceCorrection(fcard);

    // Set the number of hits per bin required in the acceptance correction histograms
    int hitsPerBin = fcard->Get("HitsPerBinAcceptance");
    fAcceptanceCorrection->SetMinCountsPerBinInclusive(hitsPerBin);

    // Create the class doing correlation analysis
    fcorrelations = new AliJIaaCorrelations( fcard, fhistos);
    cout<<endl<< " -----" <<endl;

    // If inclusive file is specified, set inclusive sampling for correlation analysis and
    // read the inclusive histograms for the acceptance correction and histogram class
    if( fInclusiveFile.Length() ) {
        fhistos->ReadInclusiveHistos(fInclusiveFile);
        fcorrelations->SetSamplingInclusive(); //kperp background and triangle. Default is flat
        fAcceptanceCorrection->ReadMixedEventHistograms(fInclusiveFile);


        cout<<"Background and acceptance sampling from " << fInclusiveFile <<endl;
    } else {
        cout << "Background and acceptance sampled from flat distributions." <<endl;
    }
    cout<< " -----" <<endl <<endl;

    // Tell the correlation analysis to use the defined acceptance correction
    fcorrelations->SetAcceptanceCorrection(fAcceptanceCorrection);

    //==================================

    // EventPool for Mixing
    fassocPool   = new AliJEventPool( fcard, fhistos, fcorrelations, fjassoc);

    fphotonList = new TClonesArray(kParticleProtoType[kJPhoton],1500);
    fchargedHadronList  = new TClonesArray(kParticleProtoType[kJHadron],1500);
    fpizeroList = new TClonesArray(kParticleProtoType[kJPizero],1500);
    ftriggList  = new TClonesArray(kParticleProtoType[fjtrigg],1500);
    fassocList  = new TClonesArray(kParticleProtoType[fjassoc],1500);
    finputList = NULL;

    fdmg = new AliJDataManager(fcard, fhistos, fcorrelations, fExecLocal);
    fdmg->SetExecLocal( fExecLocal );

    //==== Read the Data files =====
    if( fExecLocal ){
        fdmg->ChainInputStream(finputFile);

        // for grid running, numberEvents is filled by the encapsulating
        // grid task, which has access to the input handlers and can
        // extract event numbers out of it
        int nEvents = fdmg->GetNEvents();
        frunHeader = fdmg->GetRunHeader();
        cout<<"RunID = "<<frunHeader->GetRunNumber()<< " Looping over "<<nEvents<<" events"<<endl;

    } else {
        fdmg->SetRunHeader( frunHeader );
        frunHeader = fdmg->GetRunHeader();
    }

    //==== Efficiency ====
    fEfficiency = new AliJEfficiency;
    fEfficiency->SetMode( fcard->Get("EfficiencyMode") ); // 0:NoEff, 1:Period 2:RunNum 3:Auto
    if(fExecLocal) {
        fEfficiency->SetDataPath("/mnt/flustre/alice/taxi_jcorran/2013/EFF/data"); // Efficiency root file location local or alien
    } else {
        fEfficiency->SetDataPath("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data"); // Efficiency root file location local or alien
    }

    //-------------------------------------------------------------------

    fbTriggCorrel  = fcard->Get("CorrelationType")==0;
    fbLPCorrel     = fcard->Get("CorrelationType")==1;
    fMinimumPt     = fcard->GetBinBorder(kAssocType, 0);

    // Initialize counters
    fcent = -1;

    fhistos->fHMG->WriteConfig();
    fFirstEvent = kTRUE;
    fevt = -1;

    cout << "end of jcorran user create output objects ----------------" << endl;
}

void AliJIaaAnalysis::UserExec(){
    // event loop
    fevt++;

    if( fFirstEvent ) {
        //==== RunTable
        fRunTable = & AliJRunTable::GetSpecialInstance();
        fRunTable->SetRunNumber( frunHeader->GetRunNumber() );
        double sqrtS = fRunTable->GetBeamEnergy(fRunTable->GetPeriod());
        cout << "sqrts = "<< sqrtS << ",runnumber = "<< frunHeader->GetRunNumber() << endl;
        fEfficiency->SetRunNumber( frunHeader->GetRunNumber() );
        fEfficiency->Load();
        fFirstEvent = kFALSE;
    }

    fdmg->LoadEvent(fevt);
    fhistos->fhEvents->Fill( 0 );

    if(!fdmg->IsGoodEvent()) return;  // Vertex cut applied in IsGoodEvent and histo saved there too

    feventHeader  = fdmg->GetEventHeader();
    double zVert    = feventHeader->GetZVertex();
    //----------------------------------------------------------

    UInt_t triggermaskJCorran = feventHeader->GetTriggerMaskJCorran();

    if( fdmg->IsSelectedTrigger((int) triggermaskJCorran))
        fhistos->fhEvents->Fill( 5 );

    // select only some BC%4
    if( feventHeader->GetBunchCrossNumber() % 4 != fEventBC && fEventBC > -1 )
        return;

    if( fdmg->IsSelectedTrigger((int) triggermaskJCorran))
        fhistos->fhEvents->Fill( 6 );

    if(fRunTable->IsHeavyIon() || fRunTable->IsPA()){
        fcent = feventHeader->GetCentrality();
    }else  {
        fcent = 1; //ntracks;
    }

    int cBin = fcard->GetBin(kCentrType, fcent);
    if(cBin<0) return;

    if( fdmg->IsSelectedTrigger((int) triggermaskJCorran))
        fhistos->fhEvents->Fill( 7 );

    int zBin = fcard->GetBin(kZVertType, zVert); //should be alway >0; checked in fdmg->IsGoodEvent()

    // do not fill MB in case of MB mixing
    if( fdmg->IsSelectedTrigger((int) triggermaskJCorran))
        fhistos->fhZVert[cBin]->Fill(zVert);

    //------------------------------------------------------------------
    // Triggers and associated
    //----------------------ooooo---------------------------------------


    if(fjtrigg==kJHadron || fjassoc==kJHadron){
        fchargedHadronList->Clear();
        fdmg->RegisterList(fchargedHadronList, NULL, cBin, zBin, kJHadron);
        // apply efficiencies

        for( int i = 0; i < fchargedHadronList->GetEntries(); i++ )
        {
            AliJBaseTrack *triggTr = (AliJBaseTrack*)fchargedHadronList->At(i);
            double ptt = triggTr->Pt();

            double effCorr = 1./fEfficiency->GetCorrection(ptt, fHadronSelectionCut, fcent);  // here you generate warning if ptt>30
            fhistos->fhTrackingEfficiency[cBin]->Fill( ptt, 1./effCorr );
            triggTr->SetTrackEff( 1./effCorr );
        }
    }

    //---- assign input list ----
    if(fjtrigg==kJPizero)      finputList = fpizeroList;
    else if(fjtrigg==kJHadron) finputList = fchargedHadronList;
    else if(fjtrigg==kJPhoton) finputList = fphotonList;
    int noAllTriggTracks = finputList->GetEntries();
    int noAllChargedTracks = fchargedHadronList->GetEntries();
    fhistos->fhChargedMult[cBin]->Fill(noAllChargedTracks);
    fhistos->fhChargedMultCent->Fill(fcent, noAllChargedTracks>0 ? log(noAllChargedTracks) : 0);

    //---------------------------------------------
    //--    Check fcentrality outliers           ---
    //--    and count enevents only here        ---
    //---------------------------------------------
    // Need to cooperate with different parameters for different analysis cut !!!!! TODO
    //if(fRunTable->IsHeavyIon() && fCentMultLow->GetParameter(0) >0 ){
    //		double logMult = noAllChargedTracks>0 ? log(noAllChargedTracks) : 0 ;
    //		if( fcent<fCentMultLow->GetParameter(1) && logMult <fCentMultLow->Eval(fcent) ) return;
    //		if( fCentMultHigh->Eval(fcent) < logMult) return;
    //	}
    fhistos->fhChargedMultCut[cBin]->Fill(noAllChargedTracks);
    fhistos->fhCentr->Fill(fcent);
    fhistos->fhiCentr->Fill(cBin);

    // ------------------------------
    // --- calculate e-b-e vn
    // ------------------------------
    if(fRunTable->IsHeavyIon()){
        double vdelta[2] = {0};
        int    vdeltaNorm = 0;
        for(int it1=0; it1<noAllTriggTracks-1; it1++){
            AliJBaseTrack *ftk1 = (AliJBaseTrack*)finputList->At(it1);
            if(ftk1->Pt()<fMinimumPt) continue;
            for(int it2=it1+1; it2<noAllTriggTracks; it2++){
                AliJBaseTrack *ftk2 = (AliJBaseTrack*)finputList->At(it2);
                if(ftk2->Pt()<fMinimumPt) continue;
                if(fabs(ftk1->Eta() - ftk2->Eta())<1.0) continue;
                double fdphi = ftk1->DeltaPhi(*ftk2);
                fhistos->fhVN[cBin]->Fill(fdphi);
                vdelta[0] += cos(2*fdphi);
                vdelta[1] += cos(3*fdphi);
                //cout<< ftk1->Pt() <<" "<< ftk2->Pt() <<" "<< fabs(ftk1->Eta() - ftk2->Eta()) <<" "<< 2*fdphi <<" "<< 3*fdphi <<endl;
                vdeltaNorm++;
            }
        }
        vdelta[0] = vdeltaNorm>0 ? vdelta[0]/vdeltaNorm : 0;
        vdelta[1] = vdeltaNorm>0 ? vdelta[1]/vdeltaNorm : 0;
        fhistos->fhVdelta2[cBin]->Fill(vdelta[0]*100);
        fhistos->fhVdelta3[cBin]->Fill(vdelta[1]*100);
        if(vdelta[0]>0) fhistos->fpV2->Fill( fcent, sqrt(vdelta[0])*100 );
        if(vdelta[1]>0) {
            fhistos->fpV3->Fill( fcent, sqrt(vdelta[1])*100 );
            fcard->SetEventV3kv(vdelta[1]);
        }
        fhistos->fpVdeltaNorm->Fill( fcent, vdeltaNorm );
    }

    //----------------------------------------------------
    //----- Generate trigg list and find LP             --
    //----------------------------------------------------
    AliJTrackCounter *lpTrackCounter = new AliJTrackCounter();
    AliJBaseTrack *lPTr = NULL;
    int noTriggs=0;
    ftriggList->Clear();

    for(int itrack=0; itrack<noAllTriggTracks; itrack++)
    {
        AliJBaseTrack *triggTr = (AliJBaseTrack*)finputList->At(itrack);
        triggTr->SetTriggBin( fcard->GetBin(kTriggType, triggTr->Pt()) );

        double ptt = triggTr->Pt();
        double etat = triggTr->Eta();

        double effCorr = 1.0/triggTr->GetTrackEff();

        if( ptt>fMinimumPt ){
            fhistos->fhChargedPt[cBin]->Fill(ptt, effCorr );
            fhistos->fhChargedPtNoCorr[cBin]->Fill( ptt );
            fhistos->fhChargedEta->Fill(triggTr->Eta(), effCorr);
            fhistos->fhChargedPtJacek[cBin]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0); //One CANNOT do 1/ptt here!! First unfold.
        }

        if( !triggTr->IsInTriggerBin() ) continue;
        int iptt = triggTr->GetTriggBin();
        fhistos->fhIphiTrigg[cBin][iptt]->Fill( triggTr->Phi(), effCorr);
        fhistos->fhIetaTrigg[cBin][iptt]->Fill( triggTr->Eta(), effCorr);

        if( ptt > lpTrackCounter->Pt() ) {
            lpTrackCounter->Store(noTriggs, ptt, iptt);
            lPTr = triggTr;
        }
        new ((*ftriggList)[noTriggs++]) AliJBaseTrack(*triggTr);
    }


    //--------------------------------------------------
    //---   Generate assoc list and pool             ---
    //--------------------------------------------------
    fassocList->Clear();
    int noAssocs=0;
    if(fjassoc==kJPizero) finputList = fpizeroList;
    else if(fjassoc==kJHadron) finputList = fchargedHadronList;
    else if(fjassoc==kJPhoton) finputList = fphotonList;

    int noAllAssocTracks = finputList->GetEntries();

    for(int itrack=0;itrack<noAllAssocTracks; itrack++){

        AliJBaseTrack *assocTr = (AliJBaseTrack*)finputList->At(itrack);
        assocTr->SetAssocBin( fcard->GetBin(kAssocType, assocTr->Pt()) );

        if(assocTr->IsInAssocBin()){

            int ipta  = assocTr->GetAssocBin();
            double effCorr = 1.0/assocTr->GetTrackEff();
            fhistos->fhIphiAssoc[cBin][ipta]->Fill( assocTr->Phi(), effCorr);
            fhistos->fhIetaAssoc[cBin][ipta]->Fill( assocTr->Eta(), effCorr);
            new ((*fassocList)[noAssocs++]) AliJBaseTrack(*assocTr);
        }
    }


    //-----------------------------------------------
    // Leading particle pT and eta
    //-----------------------------------------------
    if( lpTrackCounter->Exists() ){
        double effCorr = 1./fEfficiency->GetCorrection(lpTrackCounter->Pt(), fHadronSelectionCut, fcent );
        fhistos->fhLPpt->Fill(lpTrackCounter->Pt(), effCorr);
        fhistos->fhLPeta->Fill(lPTr->Eta(), effCorr);
    }

    fhistos->fhAssocMult->Fill(noAssocs);

    if(noAssocs>0 ) fassocPool->AcceptList(fassocList, fcent, zVert, noAllChargedTracks, fevt);

    //------------------------------------------------------------------
    // Do the Correlation
    //----------------------ooooo---------------------------------------
    int nTriggerTracks=-1;
    nTriggerTracks = fbTriggCorrel ? noTriggs : 1;
    if(fbLPCorrel && !lpTrackCounter->Exists()) return;
    AliJBaseTrack *triggTr = NULL;

    for(int ii=0;ii<nTriggerTracks;ii++){ // trigger loop
        if (fbTriggCorrel)  triggTr = (AliJBaseTrack*)ftriggList->At(ii);
        if (fbLPCorrel)     triggTr = (AliJBaseTrack*)ftriggList->At(lpTrackCounter->GetIndex());

        double ptt = triggTr->Pt();
        int iptt   = triggTr->GetTriggBin();
        if(iptt<0) {
            cout<<"Not registered trigger ! I better stop here." <<endl;
            exit(-1);
        }
        double effCorr = 1.0/triggTr->GetTrackEff();
        fhistos->fhTriggPtBin[cBin][zBin][iptt]->Fill(ptt, effCorr);//inclusive

        for(int jj=0;jj<noAssocs;jj++){ // assoc loop
            AliJBaseTrack  *assocTr = (AliJBaseTrack*)fassocList->At(jj);
            //-------------------------------------------------------------
            fcorrelations->FillCorrelationHistos(kReal, cBin, zBin, triggTr, assocTr);
            //-------------------------------------------------------------
        } // end assoc loop
    } // end trigg loop

    // ==== Event Mixing ====
    fassocPool->Mix(ftriggList, kAzimuthFill, fcent, zVert, noAllChargedTracks, fevt, fbLPCorrel);

    //--------------------------------------------------------------
    // End of the Correlation
    //--------------------------------------------------------------

}

void AliJIaaAnalysis::Terminate() {
    // termination
    fcorrelations->PrintOut();
    fassocPool->PrintOut();
    fEfficiency->Write();
    fhistos->fHMG->Print();
}

particleType  AliJIaaAnalysis::GetParticleType(char *inchar){
    // part type
    for(int i=0;i<kNumberOfParticleTypes;i++) {
        if(strcmp(inchar,kParticleTypeStrName[i])==0) return (particleType)i;
    }
    std::cout<<"Particle type <"<<inchar<<"> not recognized"<<std::endl;
    exit(1);
}

double AliJIaaAnalysis::DeltaPhi(double phi1, double phi2) {
    // dphi
    double res =  atan2(sin(phi1-phi2), cos(phi1-phi2));
    return res > -kJPi*9./20. ? res : kJTwoPi+res ;
}
