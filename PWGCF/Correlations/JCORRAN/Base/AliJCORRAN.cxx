/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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

// jcorran main class
// used in local and grid execution

#include <TH1D.h>
#include "AliJCORRAN.h"

#include "AliJTrackCounter.h"
#include <TClonesArray.h>

#include "AliJCard.h"
#include "AliJHistos.h"
#include "AliJCorrelations.h"
#include "AliJEventPool.h"
#include "AliJDataManager.h"

#include "AliJEventHeader.h"
#include "AliJRunHeader.h"
#include "AliJTrack.h"
#include "AliJPhoton.h"
#include "AliJMCTrack.h"
#include "AliJConst.h"
#include "AliJAcceptanceCorrection.h"



#include "AliJEfficiency.h"
#include <iostream>

ClassImp(AliJCORRAN)

    AliJCORRAN::AliJCORRAN() :
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
        fphotonPool(0), 
        fassocPool(0),  
        fphotonList(0),  
        fchargedHadronList(0), 
        fpizeroList(0), 
        ftriggList(0),  
        fassocList(0), 
        fpairList(0), 
        fpairCounterList(0), 
        finputList(0), 
        fdmg(0), 
        feventHeader(0), 
        frunHeader(0), 
        fnumberEvents(0), 
        fieout(0), 
        fEventCounter(0), 
        fcent(0), 
        fncBin(0), 
        fnPttBin(0), 
        fbTriggCorrel(0), 
        fbLPCorrel(0), 
        fbLPpairCorrel(0), 
        fTrackEtaRange(0), 
        flowerPtAssocBoarder(0), 
        fCentMultLow(0),  
        fCentMultHigh(0),
        fEventBC(0),
        fSQRTS(0),
        fEfficiency(0),
        fRunTable(0),
        fIsolationR(0),
        fHadronSelectionCut(0)
{
    // constructor
}

AliJCORRAN::AliJCORRAN(Bool_t execLocal) :
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
    fphotonPool(0),
    fassocPool(0),  
    fphotonList(0),  
    fchargedHadronList(0), 
    fpizeroList(0), 
    ftriggList(0),  
    fassocList(0), 
    fpairList(0), 
    fpairCounterList(0), 
    finputList(0), 
    fdmg(0), 
    feventHeader(0), 
    frunHeader(0), 
    fnumberEvents(0), 
    fieout(0), 
    fEventCounter(0), 
    fcent(0), 
    fncBin(0), 
    fnPttBin(0), 
    fbTriggCorrel(0), 
    fbLPCorrel(0), 
    fbLPpairCorrel(0), 
    fTrackEtaRange(0), 
    flowerPtAssocBoarder(0), 
    fCentMultLow(0),  
    fCentMultHigh(0),
    fEventBC(0),
    fSQRTS(0),
    fEfficiency(0),
    fRunTable(0),
    fIsolationR(0),
    fHadronSelectionCut(0)
{
    // constructor
}

AliJCORRAN::~AliJCORRAN(){
    // destructor
}

AliJCORRAN::AliJCORRAN(const AliJCORRAN& obj) : 
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
    fphotonPool(obj.fphotonPool), 
    fassocPool(obj.fassocPool),  
    fphotonList(obj.fphotonList),  
    fchargedHadronList(obj.fchargedHadronList), 
    fpizeroList(obj.fpizeroList), 
    ftriggList(obj.ftriggList),  
    fassocList(obj.fassocList), 
    fpairList(obj.fpairList), 
    fpairCounterList(obj.fpairCounterList), 
    finputList(obj.finputList), 
    fdmg(obj.fdmg), 
    feventHeader(obj.feventHeader), 
    frunHeader(obj.frunHeader), 
    fnumberEvents(obj.fnumberEvents), 
    fieout(obj.fieout), 
    fEventCounter(obj.fEventCounter), 
    fcent(obj.fcent), 
    fncBin(obj.fncBin), 
    fnPttBin(obj.fnPttBin), 
    fbTriggCorrel(obj.fbTriggCorrel), 
    fbLPCorrel(obj.fbLPCorrel), 
    fbLPpairCorrel(obj.fbLPpairCorrel), 
    fTrackEtaRange(obj.fTrackEtaRange), 
    flowerPtAssocBoarder(obj.flowerPtAssocBoarder), 
    fCentMultLow(obj.fCentMultLow),  
    fCentMultHigh(obj.fCentMultHigh),
    fEventBC(obj.fEventBC),
    fSQRTS(obj.fSQRTS),
    fEfficiency(obj.fEfficiency),
    fRunTable(obj.fRunTable),
    fIsolationR(obj.fIsolationR),
    fHadronSelectionCut(obj.fHadronSelectionCut)
{
    // copy constructor
    JUNUSED(obj);
}

AliJCORRAN& AliJCORRAN::operator=(const AliJCORRAN& obj){
    // copy constructor
    JUNUSED(obj);
    return *this;
}


void AliJCORRAN::Initialize() const{
    // init

}

void AliJCORRAN::UserCreateOutputObjects(){
  // local init
  
  
  cout << "jcorran user create output objects ----------------" << endl;
  
  fHadronSelectionCut =int ( fcard->Get("HadronSelectionCut"));
  fIsolationR = fcard->Get("IsolationR");
  
  fhistos = new AliJHistos( fcard );
  if(fcard->Get("QualityControlLevel")>1) fhistos->Set2DHistoCreate(true);
  fhistos->CreateEventTrackHistos();
  fhistos->CreateAzimuthCorrHistos();
  //fhistos->CreateIAAMoons();
  fhistos->CreateXEHistos();
  fhistos->CreateXtHistos();
  //fhistos->CreatePairPtCosThetaStar();
  
  fhistos->fHMG->Print();
  
  fEventBC = (Int_t)(fcard->Get( "eventBC" ));
  fSQRTS = 0.;
  
  // Create a class for acceptance correction
  fAcceptanceCorrection = new AliJAcceptanceCorrection(fcard);
  
  // Set the number of hits per bin required in the acceptance correction histograms
  int hitsPerBin = fcard->Get("HitsPerBinAcceptance");
  fAcceptanceCorrection->SetMinCountsPerBinInclusive(hitsPerBin);
  
  // Create the class doing correlation analysis
  fcorrelations = new AliJCorrelations( fcard, fhistos);
  cout<<endl<< " -----" <<endl;
  
  // If inclusive file is specified, set inclusive sampling for correlation analysis and
  // read the inclusive histograms for the acceptance correction and histogram class
  if( fInclusiveFile.Length() ) {
    fhistos->ReadInclusiveHistos(fInclusiveFile);
    fcorrelations->SetSampligInclusive(); //kperp background and triangle. Default is flat
    fAcceptanceCorrection->ReadMixedEventHistograms(fInclusiveFile);
    cout<<"Sampling kperp and triangle from " << fInclusiveFile <<endl;
  } else cout << "Sampling kperp and triangle from flat" <<endl;
  cout<< " -----" <<endl <<endl;
  
  // Tell the correlation analysis to use the defined acceptance correction
  fcorrelations->SetAcceptanceCorrection(fAcceptanceCorrection);
  
   //==================================
  
  //cout<<kParticleTypeStrName[kPhoton]<<" "<<kParticleTypeStrName[fjtrigg]<<endl;
  // EventPool for Mixing
  fphotonPool  = new AliJEventPool( fcard, fhistos, fcorrelations, kJPhoton);  // for pi0 mass
  fassocPool   = new AliJEventPool( fcard, fhistos, fcorrelations, fjassoc);
  
  fphotonList = new TClonesArray(kParticleProtoType[kJPhoton],1500);
  //     TClonesArray *cellList = new TClonesArray("AliJCaloCell",1500);
  fchargedHadronList  = new TClonesArray(kParticleProtoType[kJHadron],1500);
  fpizeroList = new TClonesArray(kParticleProtoType[kJPizero],1500);
  ftriggList  = new TClonesArray(kParticleProtoType[fjtrigg],1500);
  fassocList  = new TClonesArray(kParticleProtoType[fjassoc],1500);
  fpairList     = new TClonesArray(kParticleProtoType[fjtrigg],1500);
  fpairCounterList  = new TClonesArray("AliJTrackCounter",1500);
  finputList = NULL;
  //TClonesArray *isolPizeroList  = new TClonesArray("AliPhJPiZero",1500);
  
  fdmg = new AliJDataManager(fcard, fhistos, fcorrelations, fExecLocal);
  fdmg->SetExecLocal( fExecLocal );
  
  //==== Read the Data files =====
  if( fExecLocal ){
    fdmg->ChainInputStream(finputFile);
    // TODO: run header is not supposed to be here
    // doesn't work fSQRTS = 2.*frunHeader->GetBeamEnergy();
    
    // for grid running, numberEvents is filled by the encapsulating
    // grid task, which has access to the input handlers and can
    // extract event numbers out of it
    fnumberEvents = fdmg->GetNEvents();
    fieout = fnumberEvents/20;
    frunHeader = fdmg->GetRunHeader();
    cout<<"RunID = "<<frunHeader->GetRunNumber()<< " Looping over "<<fnumberEvents<<" events"<<endl;
    
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
	fEventCounter=0;

	fcent = -1;
	fncBin = fcard->GetNoOfBins( kCentrType );
	fnPttBin = fcard->GetNoOfBins( kTriggType );

	fbTriggCorrel  = fcard->Get("CorrelationType")==0;
	fbLPCorrel     = fcard->Get("CorrelationType")==1;
	fbLPpairCorrel = fcard->Get("CorrelationType")==2;
	fTrackEtaRange = fcard->Get("EtaRange");
	flowerPtAssocBoarder = fcard->GetBinBorder(kAssocType, 0);

	//==== Limit the no of track for each fcentrality ======
	fCentMultLow = new TF1("fCentMultLow","[0]*sqrt([1]-x)+[2]", 0, 75);
	fCentMultHigh = new TF1("fCentMultHigh","[0]*sqrt([1]-x)+[2]", 0, 90);

	// fast parameter load


	fhistos->fHMG->WriteConfig();
	fFirstEvent = kTRUE;
	fevt = -1;

	cout << "end of jcorran user create output objects ----------------" << endl;

}

void AliJCORRAN::UserExec(){
	// event loop
	fevt++;

	// TODO Add train mode if( !fExecLocal && fFirstEvent ){
	if( 0 && !fExecLocal && fFirstEvent ){
		fdmg->ChainInputStream("");
		fieout = fnumberEvents/20;
		if (fieout<1) fieout=1;
		cout<<"RunID = "<<frunHeader->GetRunNumber()<< " Looping over "<<fnumberEvents<<" events"<<endl;

		//==== Efficiency ====
		// TODO run can be different in a job?
		fhistos->fhEventPerRun->Fill(fnumberEvents>0?log(fnumberEvents):1);
		fFirstEvent = kFALSE;
	}

	// TODO 
	if( fExecLocal ) {
		if(fevt % fieout == 0) cout << fevt << "\t" << int(float(fevt)/fnumberEvents*100) << "%" << endl ;
	}

	if( fFirstEvent ) {
		//==== RunTable
		fRunTable = & AliJRunTable::GetSpecialInstance();
		fRunTable->SetRunNumber( frunHeader->GetRunNumber() );
		fSQRTS = fRunTable->GetBeamEnergy(fRunTable->GetPeriod());
		cout << "sqrts = "<< fSQRTS << ",runnumber = "<< frunHeader->GetRunNumber() << endl;
		fEfficiency->SetRunNumber( frunHeader->GetRunNumber() );
		fEfficiency->Load();
		fFirstEvent = kFALSE;
	}

	if(fRunTable->IsHeavyIon()){
		fCentMultLow->SetParameters( fcard->Get("CentMultCutLow",0),  fcard->Get("CentMultCutLow",1),  fcard->Get("CentMultCutLow",2) );
		fCentMultHigh->SetParameters(fcard->Get("CentMultCutHigh",0), fcard->Get("CentMultCutHigh",1), fcard->Get("CentMultCutHigh",2) );
		//fCentMultLow->Print();
		//fCentMultHigh->Print();
	}

	fdmg->LoadEvent(fevt); // to be here to load EventHeader
	fhistos->fhEvents->Fill( 0 );

	if(!fdmg->IsGoodEvent()) return;  // Vertex cut applied in IsGoodEvent and histo saved there too

	feventHeader  = fdmg->GetEventHeader();
	double zVert    = feventHeader->GetZVertex();
	//----------------------------------------------------------

	UInt_t triggermaskJCorran = feventHeader->GetTriggerMaskJCorran();
	//cout << fevt <<"\t"<< zVert << "\t"<< triggermaskJCorran <<  endl;

	if( fdmg->IsSelectedTrigger((int) triggermaskJCorran))
		fhistos->fhEvents->Fill( 5 );

	// select only some BC%4
	if( feventHeader->GetBunchCrossNumber() % 4 != fEventBC &&
			fEventBC > -1 )
		return;

	if( fdmg->IsSelectedTrigger((int) triggermaskJCorran))
		fhistos->fhEvents->Fill( 6 );

	if(fRunTable->IsHeavyIon() || fRunTable->IsPA()){
		fcent = feventHeader->GetCentrality();
	}else  {
		fcent = 1; //ntracks;
	}

	//cout<<"evt="<<fevt <<" c="<<  fcent <<endl;
	int cBin        = fcard->GetBin(kCentrType, fcent);
	if(cBin<0) return;

	if( fdmg->IsSelectedTrigger((int) triggermaskJCorran))
		fhistos->fhEvents->Fill( 7 );

	int zBin        = fcard->GetBin(kZVertType, zVert); //should be alway >0; checked in fdmg->IsGoodEvent()

	// do not fill MB in case of MB mixing
	if( fdmg->IsSelectedTrigger((int) triggermaskJCorran))
		fhistos->fhZVert[cBin]->Fill(zVert);

	//temporary fevent skip
	//Trigger to be selected from the JCorran trigger mask is specified in the fCard
	//         if(! fdmg->IsSelectedTrigger((int) triggermaskJCorran))
	//           continue;

	fEventCounter++;

	//==== QA Event
	fhistos->fhV0AMult[cBin]->Fill( feventHeader->GetV0AMult() );

	//------------------------------------------------------------------
	// Triggers and associated
	//----------------------ooooo---------------------------------------

	if(fjtrigg==kJPizero || fjassoc==kJPizero || fjtrigg==kJPhoton || fjassoc==kJPhoton){
	} // pizero || photon
	if(fjtrigg==kJHadron || fjassoc==kJHadron){
		fchargedHadronList->Clear();
		fdmg->RegisterList(fchargedHadronList, NULL, cBin, zBin, kJHadron);
		// apply efficiencies

		for( int i = 0; i < fchargedHadronList->GetEntries(); i++ ){

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
			if(ftk1->Pt()<flowerPtAssocBoarder) continue;
			for(int it2=it1+1; it2<noAllTriggTracks; it2++){
				AliJBaseTrack *ftk2 = (AliJBaseTrack*)finputList->At(it2);
				if(ftk2->Pt()<flowerPtAssocBoarder) continue;
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
	//----- Fiducial cut should be used in AliJCorrelations--
	//----------------------------------------------------
	AliJTrackCounter *lpTrackCounter = new AliJTrackCounter(), *lpPairCounter = new AliJTrackCounter();
	AliJBaseTrack *lPTr = NULL;
	int noTriggs=0;
	ftriggList->Clear();
	for(int itrack=0; itrack<noAllTriggTracks; itrack++){
		AliJBaseTrack *triggTr = (AliJBaseTrack*)finputList->At(itrack);
		triggTr->SetTriggBin( fcard->GetBin(kTriggType, triggTr->Pt()) );

		double ptt = triggTr->Pt();
		double etat = triggTr->Eta();
		//TODO iCut == 0;

		double effCorr = 1.0/triggTr->GetTrackEff();

		if( ptt>flowerPtAssocBoarder ){
			//FK//double effCorr = 1./fcard->TrackEfficiency(ptt, fcent);  // here you generate warning if ptt>30
			//double effCorr = 1./fcard->TrackEfficiency(ptt, etat, cBin);  // here you generate warning if ptt>30
			fhistos->fhChargedPt[cBin]->Fill(ptt, effCorr );
			fhistos->fhChargedPtNoCorr[cBin]->Fill( ptt );
			fhistos->fhChargedEta->Fill(triggTr->Eta(), effCorr);
			//fhistos->fhChargedPtJacek[cBin]->Fill(ptt, effCorr );
			fhistos->fhChargedPtJacek[cBin]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0); //One CANNOT do 1/ptt here!! First unfold.
			if( -0.8<etat && etat<-0.2) fhistos->fhChargedPtJacekEta[cBin][0]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0);
			if( -0.2<etat && etat< 0.3) fhistos->fhChargedPtJacekEta[cBin][1]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0);
			if(  0.3<etat && etat< 0.8) fhistos->fhChargedPtJacekEta[cBin][2]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0);
			fhistos->fhChargedPtFiete->Fill(ptt, effCorr );
		}

		if( !triggTr->IsInTriggerBin() ) continue;
		//FK//triggTr->SetTrackEff( fcard->TrackEfficiency(triggTr->Pt(), fcent) );
		int iptt = triggTr->GetTriggBin();
		fhistos->fhIphiTrigg[cBin][iptt]->Fill( triggTr->Phi(), effCorr);
		fhistos->fhIetaTrigg[cBin][iptt]->Fill( triggTr->Eta(), effCorr);

		if( ptt > lpTrackCounter->Pt() ) {
			lpTrackCounter->Store(noTriggs, ptt, iptt);
			lPTr = triggTr;
		}
		//cout <<"1 ptt="<< triggTr->Pt() <<" bin="<< triggTr->GetTriggBin() << " st=" << triggTr->GetStatus() << " trigg eff=" << triggTr->GetTrackEff() << endl; 
		new ((*ftriggList)[noTriggs++]) AliJBaseTrack(*triggTr);
		fhistos->fhTriggMult[cBin][iptt]->Fill(noAllTriggTracks);
	}

	//----------------------------------------------------
	//----- Find sum of two leading particles ------------
	//----------------------------------------------------
	if(fbLPpairCorrel){
		int noPairs=0;
		fpairList->Clear();
		fpairCounterList->Clear();
		for(int ii=0;ii<noTriggs-1;ii++){
			AliJBaseTrack *triggTr1 = (AliJBaseTrack*)ftriggList->At(ii);
			for(int jj=ii+1;jj<noTriggs;jj++){
				AliJBaseTrack *triggTr2 = (AliJBaseTrack*)ftriggList->At(jj);
				TLorentzVector lVPair = triggTr1->GetLorentzVector()+triggTr2->GetLorentzVector();
				double fdphi = DeltaPhi(triggTr1->Phi(), triggTr2->Phi());
				new ((*fpairList)[noPairs]) AliJBaseTrack(lVPair);
				new ((*fpairCounterList)[noPairs++]) 
					AliJTrackCounter( triggTr1->GetID(), triggTr2->GetID(), ii, jj,  fdphi );
			}
		}

		// ----- Find the Leading Pair --------------------------
		AliJBaseTrack *pairTr = NULL;
		for(int ii=0;ii<noPairs;ii++){
			pairTr = (AliJBaseTrack*)fpairList->At(ii);
			AliJTrackCounter   *pairCounter = (AliJTrackCounter*)fpairCounterList->At(ii);
			//cout<<pairTr->Pt()<<endl;    pairCounter->Print();
			if(pairTr->Pt() > lpPairCounter->Pt() && fabs(pairCounter->GetPairDPhi())<1.0) {
				int iptt  = fcard->GetBin(kTriggType, pairTr->Pt());
				lpPairCounter = pairCounter;
				lpPairCounter->Store(ii, pairTr->Pt(), iptt); 
			}
		}
		if(lpPairCounter->Exists()){
			pairTr->SetParticleType(kJHadron);
			//double effCorr = 1./fcard->TrackEfficiency(lpTrackCounter->GetLPpt());
			fhistos->fhLPpairPt->Fill(pairTr->Pt());
		}
	}

	//--------------------------------------------------
	//---   Generate assoc list and pool             ---
	//--------------------------------------------------
	fassocList->Clear();
	int noAssocs=0;
	double  sumPtAroundLP = 0;
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
		// determine an activity in the cone around trigger 
		if(lpTrackCounter->Exists() && (lPTr->GetID()!=assocTr->GetID()) && (0.5 < assocTr->Pt())){
			double dPhi = DeltaPhi( assocTr->Phi(), lPTr->Phi() );
			double dEta = assocTr->Eta() - lPTr->Eta();
			if(fIsolationR > sqrt(dPhi*dPhi+dEta*dEta))  sumPtAroundLP += assocTr->Pt();//FK//
		}
	}

	FillXtHistos(finputList, lpTrackCounter);

	//-----------------------------------------------
	// Set the isolation flag
	//-----------------------------------------------
	if( lpTrackCounter->Exists() ){
		//double effCorr = 1./fcard->TrackEfficiency(lpTrackCounter->Pt(), fcent);
		//TODO ??? AGAIN AGAIN???
		double effCorr = 1./fEfficiency->GetCorrection(lpTrackCounter->Pt(), fHadronSelectionCut, fcent );
		fhistos->fhLPpt->Fill(lpTrackCounter->Pt(), effCorr);
		fhistos->fhLPeta->Fill(lPTr->Eta(), effCorr);
		fhistos->fhBkgActivity[lpTrackCounter->GetPtBin()]->Fill(sumPtAroundLP/lpTrackCounter->Pt());

		if( sumPtAroundLP/lpTrackCounter->Pt()< fcard->GetCutOnBkgActivity() &&     ///relative activity
				( fabs(lPTr->Eta()) < (fTrackEtaRange - fIsolationR) )   ){ //fiducial cut

			fhistos->fhIsolatedLPpt->Fill(lpTrackCounter->Pt(), effCorr );

			AliJBaseTrack* lpTrigger = (AliJBaseTrack*) ftriggList->At(lpTrackCounter->GetIndex());
			lpTrigger->SetIsIsolated(1);
		}
	}

	fhistos->fhAssocMult->Fill(noAssocs);
	//cout<"Triggs = "<<<noTriggs<<" Assoc = "<<noAssocs<<" all = "<<noAllTracks<<endl;

	//============================================================
	//there is obviously a problem when you do "fixed" correlation.
	//In this case the noTriggs==0 condition is never fullfilled
	//============================================================
	//if(noTriggs==0 || ((fjtrigg==kPizero) && (fjassoc==kPizero)) )
	//if(leadingPt<1.5) //should be fixed
	//if(noTriggs==0 && noAssocs>0 ){}
	if(noAssocs>0 ) fassocPool->AcceptList(fassocList, fcent, zVert, noAllChargedTracks, fevt);

	//------------------------------------------------------------------
	// Do the Correlation 
	//----------------------ooooo---------------------------------------
	int noTriggTracs=-1;
	noTriggTracs = fbTriggCorrel ? noTriggs : 1;
	if(fbLPCorrel && !lpTrackCounter->Exists()) return;
	if(fbLPpairCorrel && !lpPairCounter->Exists()) return;
	AliJBaseTrack *triggTr = NULL;

	for(int ii=0;ii<noTriggTracs;ii++){ // trigger loop 
		if (fbTriggCorrel)  triggTr = (AliJBaseTrack*)ftriggList->At(ii);
		if (fbLPCorrel)     triggTr = (AliJBaseTrack*)ftriggList->At(lpTrackCounter->GetIndex());
		if (fbLPpairCorrel) triggTr = (AliJBaseTrack*)fpairList->At(lpPairCounter->GetIndex());
		//for(int ii=0;ii<noIsolPizero;ii++){ // trigger loop }
		//AliJBaseTrack *triggTr = (AliJBaseTrack*)isolPizeroList->At(ii);
		double ptt = triggTr->Pt();
		int iptt   = triggTr->GetTriggBin(); 
		if(iptt<0) {
			cout<<"Not registered trigger ! I better stop here." <<endl; 
			exit(-1);
		}
		double effCorr = 1.0/triggTr->GetTrackEff();
		fhistos->fhTriggPtBin[cBin][zBin][iptt]->Fill(ptt, effCorr);//inclusive

		if(triggTr->GetIsIsolated()>0) fhistos->fhTriggPtBinIsolTrigg[kReal][cBin][iptt]->Fill(ptt, effCorr);

		for(int jj=0;jj<noAssocs;jj++){ // assoc loop
			AliJBaseTrack  *assocTr = (AliJBaseTrack*)fassocList->At(jj);
			//assocTr->PrintOut("assoc track");
			if(fbLPpairCorrel && 
					(assocTr->GetID()==lpPairCounter->GetPairTrackID(0) || 
					 assocTr->GetID()==lpPairCounter->GetPairTrackID(1)) ) continue;
			//-------------------------------------------------------------
			fcorrelations->FillAzimuthHistos(kReal, cBin, zBin, triggTr, assocTr);
			//-------------------------------------------------------------
		} // end assoc loop
	} // end trigg loop
	// == mix trigg with assoc
	if (fbLPpairCorrel) 
		fassocPool->Mix(fpairList,  kAzimuthFill, fcent, zVert, noAllChargedTracks, fevt);
	else
		fassocPool->Mix(ftriggList, kAzimuthFill, fcent, zVert, noAllChargedTracks, fevt, fbLPCorrel);

	//--------------------------------------------------------------
	// End of the Correlation
	//--------------------------------------------------------------



}

void AliJCORRAN::Terminate() {
	// termination

	/*  TODO
		for (int hic = 0;hic < fcard->GetNoOfBins(kCentrType);hic++){
		ScaleNotEquidistantHisto( fhistos->fhChargedPt[hic], 1);
		ScaleNotEquidistantHisto( fhistos->fhChargedPtNoCorr[hic], 1);
		ScaleNotEquidistantHisto( fhistos->fhChargedPtJacek[hic], 1);
		}
		ScaleNotEquidistantHisto( fhistos->fhLPpt, 1);
		ScaleNotEquidistantHisto( fhistos->fhLPpairPt, 1);
		ScaleNotEquidistantHisto( fhistos->fhIsolatedLPpt, 1);
		ScaleNotEquidistantHisto( fhistos->fhChargedPtFiete, 1);
		*/

	//    cout<<"MB's="<<noMB<<" "<<" ERT's="<<noERT<<endl;
	fcorrelations->PrintOut();
	fassocPool->PrintOut();
	fEfficiency->Write();
	fhistos->fHMG->Print();
	//fhistos->fHMG->Write();
}

particleType  AliJCORRAN::GetParticleType(char *inchar){
	// part type
	for(int i=0;i<kNumberOfParticleTypes;i++) {
		//cout<<"<"<<inchar<<"> <"<<particleTypeStr[i]<<"> "<<strcmp(inchar,particleTypeStr[i])<<" "<<(particleType)i<<endl;
		if(strcmp(inchar,kParticleTypeStrName[i])==0) return (particleType)i;
	}
	std::cout<<"Particle type <"<<inchar<<"> not recognized"<<std::endl;
	exit(1);
}

double AliJCORRAN::DeltaPhi(double phi1, double phi2) {
	// dphi
	double res =  atan2(sin(phi1-phi2), cos(phi1-phi2));
	//return res>-kJPi/3.0 ? res : kJTwoPi+res ; 
	return res>-kJPi*9./20. ? res : kJTwoPi+res ; 
}

void AliJCORRAN::ScaleNotEquidistantHisto(TH1D *hid, const double sc=1){
	// histo scaler
	for(int i=1;i<= hid->GetNbinsX();i++){
		hid->SetBinContent(i,hid->GetBinContent(i)*sc/hid->GetBinWidth(i));
		hid->SetBinError(i,hid->GetBinError(i)*sc/hid->GetBinWidth(i));
	}   
}

// Fill xT histograms (Esko)
void AliJCORRAN::FillXtHistos(TClonesArray *inputList, AliJTrackCounter *lpTrackCounter){
	// Here should be a comment describing this method 
	JUNUSED(lpTrackCounter);

	enum xTtype { kInclusive=0, kIsolated=1, kIsolatedLP=2} ;
	double lowerPtCut = 0.2;
	double isolationR = 0.4;
	double fCutOnBkgActivity = 0.10;
	int cBin  = fcard->GetBin(kCentrType, fcent);
	double sqrts = fSQRTS;
	int noAllTriggTracks = inputList->GetEntries();
	//cout << "Entering Esko's analysis loop" << endl;
	//cout << "Sqrts = " << sqrts << " pT = " << ptt << " xT = " << xtt << endl;
	for(int itrack=0; itrack<noAllTriggTracks; itrack++){
		AliJBaseTrack *triggTr = (AliJBaseTrack*)inputList->At(itrack);
		double  sumPtAroundTrigger = 0;
		double ptt = triggTr->Pt();
		double xtt = 2.0 * ptt / (1.0 * sqrts);
		double effCorr = 1.0/triggTr->GetTrackEff();
		fhistos->fhPtForXt[kInclusive][cBin]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0);
		fhistos->fhXt[kInclusive][cBin]->Fill(xtt, effCorr );
		fhistos->fhXtWeighted[kInclusive][cBin]->Fill(xtt, effCorr*1.0/xtt );

		for(int jj=0;jj<inputList->GetEntriesFast();jj++){ // assoc loop
			AliJBaseTrack *assocTr = (AliJBaseTrack*)inputList->At(jj);
			if(triggTr->GetID()==assocTr->GetID()) continue;
			double pta = assocTr->Pt();
			// determine the activity in the cone around trigger (Esko)
			if( lowerPtCut < pta ){
				double dPhi = DeltaPhi( assocTr->Phi(), triggTr->Phi() );
				double dEta = assocTr->Eta() - triggTr->Eta();
				if(isolationR > sqrt(dPhi*dPhi+dEta*dEta))  sumPtAroundTrigger += assocTr->Pt();//FK//
			}
			//cout << "  Assoc number " << assocCounter++ << endl;
		}
		// If pT sum is below the limit, fill to the isolated histos
		if( sumPtAroundTrigger/ptt < fCutOnBkgActivity &&     ///relative activity
				( fabs(triggTr->Eta()) < (fTrackEtaRange - isolationR) )   ){ //fiducial cut
			// pT and xT
			// kInclusive kIsolated, kLPIsolated
			fhistos->fhPtForXt[kIsolated][cBin]->Fill(ptt,effCorr*1.0/ptt);
			fhistos->fhXt[kIsolated][cBin]->Fill(xtt, effCorr );
			fhistos->fhXtWeighted[kIsolated][cBin]->Fill(xtt,effCorr*1.0/xtt);
		}
		//cout<<"pT sum around trigger = " << sumPtAroundTrigger << endl;
	}

} // end FillXtHistos
