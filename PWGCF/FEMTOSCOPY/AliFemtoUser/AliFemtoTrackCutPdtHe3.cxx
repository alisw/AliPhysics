//==============================================================\\
// dowang track cut for p-d/t/He3 analysis                      \\
// deuteron part refer to AliFemtoWRzTrackCut.h                 \\
//==============================================================\\

#include "AliFemtoTrackCutPdtHe3.h"
#include "TLorentzVector.h"
AliFemtoTrackCutPdtHe3::AliFemtoTrackCutPdtHe3():
AliFemtoESDTrackCut()
{
    fNsigmaP    = 3.;
    fNsigmaD    = 3.;
    fNsigmaT    = 3.;
    fNsigmaHe3  = 3.;
    fNsigmaPi  = 3.;
    fNsigmaRejection = 3.;
    
    SwitchMom_p = 0.5;
    SwitchMom_d = 0.8;
    SwitchMom_t = 999;
    SwitchMom_He3 = 999;
    SwitchMom_Pi = 0.5;

    fdEdxcut = 0;
    fOtherNsigmacut = 0;
    fMinTPCNCrossRows = 0;
    fMinTPCFoundFraction = 0.;
   
    fUseTOFMassCut = 0; 
    TOFMassLowLimit = 0.;
    TOFMassUpLimit = 10.;

    fUsePtCut = 0;
    // for Nsigma reject
    fOnlyTPCreject = 1;
    
    fPionHe3cut = 0;
    fUseDCAvsPt_cut = 0;
    fInversePID = 0;

    MinHe3TPCSignal = 10.;
    MaxHe3TPCSignal = 5000.;

    fRejectPions  = 0;
    fTPCThreshold = 0.7;

    fdNSigmaVspTcut = 0;
    fdcutline_k = 10/1.1;
    fdcutline_b = -15.;

    fUse2pT = 0;
    AlldEdxmode = 0;

    wiolaCrossCheck = 0;
    fUsePtotalCut = 0;
    MinPtotal = 0.;
    MaxPtotal = 100.;
    pionrejectcut = 2.;

}

AliFemtoTrackCutPdtHe3::AliFemtoTrackCutPdtHe3(const AliFemtoTrackCutPdtHe3 &aCut) : 
AliFemtoESDTrackCut(aCut)
{
    //copy constructor 
    fNsigmaP = aCut.fNsigmaP;
    fNsigmaD = aCut.fNsigmaD;
    fNsigmaT = aCut.fNsigmaT;
    fNsigmaHe3 = aCut.fNsigmaHe3;
    fNsigmaPi = aCut.fNsigmaPi;
    fNsigmaRejection = aCut.fNsigmaRejection;

    SwitchMom_p = aCut.SwitchMom_p;
    SwitchMom_d = aCut.SwitchMom_d;
    SwitchMom_t = aCut.SwitchMom_t;
    SwitchMom_He3 = aCut.SwitchMom_He3;
    SwitchMom_Pi = aCut.SwitchMom_Pi;
    fdEdxcut = aCut.fdEdxcut;
    fOtherNsigmacut = aCut.fOtherNsigmacut;
    fMinTPCNCrossRows = aCut.fMinTPCNCrossRows;
    fMinTPCFoundFraction = aCut.fMinTPCFoundFraction;

    fUseTOFMassCut = aCut.fUseTOFMassCut;
    TOFMassLowLimit = aCut.TOFMassLowLimit;
    TOFMassUpLimit = aCut.TOFMassUpLimit;

    fUsePtCut = aCut.fUsePtCut;
    fOnlyTPCreject = aCut.fOnlyTPCreject;
    fPionHe3cut = aCut.fPionHe3cut;
    fUseDCAvsPt_cut = aCut.fUseDCAvsPt_cut;
    fInversePID = aCut.fInversePID;

    MinHe3TPCSignal = aCut.MinHe3TPCSignal;
    MaxHe3TPCSignal = aCut.MaxHe3TPCSignal;

    fRejectPions  = aCut.fRejectPions;
    fTPCThreshold = aCut.fTPCThreshold;


    fdNSigmaVspTcut = aCut.fdNSigmaVspTcut;
    fdcutline_k = aCut.fdcutline_k;
    fdcutline_b = aCut.fdcutline_b;
   
    fUse2pT = aCut.fUse2pT;
    
    AlldEdxmode = aCut.AlldEdxmode;

    wiolaCrossCheck = aCut.wiolaCrossCheck;
    
    fUsePtotalCut = aCut.fUsePtotalCut;
       MinPtotal = aCut.MinPtotal;
    MaxPtotal = aCut.MaxPtotal;
  pionrejectcut = aCut.pionrejectcut; 
}

AliFemtoTrackCutPdtHe3::~AliFemtoTrackCutPdtHe3()
{
  // Destructor
}


AliFemtoTrackCutPdtHe3& AliFemtoTrackCutPdtHe3::operator=(const AliFemtoTrackCutPdtHe3& aCut)
{
    // assignment operator
    if (this == &aCut)
        return *this;

    AliFemtoESDTrackCut::operator=(aCut);

    fNsigmaP = aCut.fNsigmaP;
    fNsigmaD = aCut.fNsigmaD;
    fNsigmaT = aCut.fNsigmaT;
    fNsigmaHe3 = aCut.fNsigmaHe3;
    fNsigmaPi = aCut.fNsigmaPi;
    
    fNsigmaRejection = aCut.fNsigmaRejection;
    SwitchMom_p = aCut.SwitchMom_p;
    SwitchMom_d = aCut.SwitchMom_d;
    SwitchMom_t = aCut.SwitchMom_t;
    SwitchMom_He3 = aCut.SwitchMom_He3;
    SwitchMom_Pi = aCut.SwitchMom_Pi;

    fdEdxcut = aCut.fdEdxcut;
    fOtherNsigmacut = aCut.fOtherNsigmacut;
    fMinTPCNCrossRows = aCut.fMinTPCNCrossRows;
    fMinTPCFoundFraction = aCut.fMinTPCFoundFraction;

    fUseTOFMassCut = aCut.fUseTOFMassCut;
    TOFMassLowLimit = aCut.TOFMassLowLimit;
    TOFMassUpLimit = aCut.TOFMassUpLimit;

    fUsePtCut = aCut.fUsePtCut;
    fPionHe3cut = aCut.fPionHe3cut;
    fUseDCAvsPt_cut = aCut.fUseDCAvsPt_cut;
    fInversePID = aCut.fInversePID;
    MinHe3TPCSignal = aCut.MinHe3TPCSignal;
    MaxHe3TPCSignal = aCut.MaxHe3TPCSignal;

    fRejectPions  = aCut.fRejectPions;
    fTPCThreshold = aCut.fTPCThreshold;

    fdNSigmaVspTcut = aCut.fdNSigmaVspTcut;
    fdcutline_k = aCut.fdcutline_k;
    fdcutline_b = aCut.fdcutline_b;

    fUse2pT = aCut.fUse2pT;

    AlldEdxmode = aCut.AlldEdxmode;

    wiolaCrossCheck = aCut.wiolaCrossCheck;
    fUsePtotalCut = aCut.fUsePtotalCut;
MinPtotal = MinPtotal;
MaxPtotal = MaxPtotal;
pionrejectcut = aCut.pionrejectcut;
    return *this;
}

bool AliFemtoTrackCutPdtHe3::Pass(const AliFemtoTrack* track){
      if(AlldEdxmode){
	if(fCharge==0){
		fNTracksFailed++;
                return false;
	}
	  if (fCharge != 0 && (track->Charge() != fCharge)) {
        	fNTracksFailed++;
        	return false;
	    }

        fNTracksPassed++ ;
        return true;
     
    }

 
    if (fStatus && (track->Flags() & fStatus) != fStatus) {
        return false;
    }
    if (fRemoveKinks && (track->KinkIndex(0) || track->KinkIndex(1) || track->KinkIndex(2))) {
        return false;
    }
    if (fRemoveITSFake && track->ITSncls() < 0) {
        return false;
    }
    if (fminTPCclsF > track->TPCnclsF()) {
        return false;
    }
    if (fminTPCncls > track->TPCncls()) {
        return false;
    }
    if (fminITScls > track->ITSncls()) {
        return false;
    }

    if (fMaxImpactXY < TMath::Abs(track->ImpactD())) {
        return false;
    }
    if (fMinImpactXY > TMath::Abs(track->ImpactD())) {
        return false;
    }
    if (fMaxImpactZ < TMath::Abs(track->ImpactZ())) {
        return false;
    }
    if (fMaxSigmaToVertex < track->SigmaToVertex()) {
        return false;
    }

    if (track->ITSncls() > 0 && (track->ITSchi2() / track->ITSncls()) > fMaxITSchiNdof) {
        return false;
    }

    if (track->TPCchi2perNDF() > fMaxTPCchiNdof) {
        return false;
    }
    if ( fMinTPCNCrossRows && ( fMinTPCNCrossRows > track->TPCNCrossedRows()) ){
   	return false;
    }

    if ( fMinTPCFoundFraction!=0. ){
	float tmpFraction = track->TPCnclsF() > 0 ? ((float)track->TPCNCrossedRows())/track->TPCnclsF() : 0;
	if( fMinTPCFoundFraction > tmpFraction ){
        	return false;
	}
    }

    // ITS cluster requirenments
    for (Int_t i = 0; i < 3; i++) {
        if (!CheckITSClusterRequirement(fCutClusterRequirementITS[i], track->HasPointOnITSLayer(i * 2), track->HasPointOnITSLayer(i*2+1))) {
        return false;
        }
    }

    if (fLabel) {
        if (track->Label() < 0) {
        fNTracksFailed++;
        return false;
        }
    }
    if (fCharge != 0 && (track->Charge() != fCharge)) {
        fNTracksFailed++;
        return false;
    }


    Bool_t tTPCPidIn = (track->Flags() & AliFemtoTrack::kTPCpid) > 0;
    Bool_t tITSPidIn = (track->Flags() & AliFemtoTrack::kITSpid) > 0;
    Bool_t tTOFPidIn = (track->Flags() & AliFemtoTrack::kTOFpid) > 0;
    
    // no use!
    //const double momentum = track->P().Mag();
    
    
//    if (fMinPforTOFpid > 0
//        && fMinPforTOFpid < momentum && momentum < fMaxPforTOFpid
//        && !tTOFPidIn) {
//        fNTracksFailed++;
//        return false;
//    }
//
//    if (fMinPforTPCpid > 0
//        && fMinPforTPCpid < momentum && momentum < fMaxPforTPCpid
//        && !tTPCPidIn) {
//        fNTracksFailed++;
//        return false;
//    }
//
//    if (fMinPforITSpid > 0
//        && fMinPforITSpid < momentum && momentum < fMaxPforITSpid
//        && !tITSPidIn) {
//        fNTracksFailed++;
//        return false;
//    }

    // dowang for He3, it needs change
    float tEnergy = 0.;
    float tRapidity = 0.;
    float tPt = 0.;
    float tEta = 0.;
    float tTotalP = 0.;

    if(fUse2pT == 1){
        TLorentzVector thisTrackMom;
        tEnergy = ::sqrt(track->P().Mag2() * 4. + fMass * fMass);
        thisTrackMom.SetPxPyPzE(2.*track->P().x(),2.*track->P().y(),2.*track->P().z(),tEnergy);
        if (tEnergy-2.*track->P().z() != 0 && (tEnergy + 2. * track->P().z()) / (tEnergy- 2. * track->P().z()) > 0)
            //tRapidity = 0.5 * ::log((tEnergy + 2. * track->P().z())/(tEnergy- 2. * track->P().z()));
            tRapidity = thisTrackMom.Rapidity();
        tPt = thisTrackMom.Pt();
        tTotalP = thisTrackMom.P();
        tEta = thisTrackMom.PseudoRapidity();
    }
    else{
        tEnergy = ::sqrt(track->P().Mag2() + fMass * fMass);
        if (tEnergy-track->P().z() != 0 && (tEnergy + track->P().z()) / (tEnergy-track->P().z()) > 0)
            tRapidity = 0.5 * ::log((tEnergy + track->P().z())/(tEnergy-track->P().z()));
        tPt = track->P().Perp();
	tTotalP = track->P().Mag();
        tEta = track->P().PseudoRapidity();
    }

    // not use
    if (fMaxImpactXYPtOff < 999.0) {
        if ((fMaxImpactXYPtOff + fMaxImpactXYPtNrm*TMath::Power(tPt, fMaxImpactXYPtPow)) < TMath::Abs(track->ImpactD())) {
        fNTracksFailed++;
        return false;
        }
    }

    if ((tRapidity < fRapidity[0]) || (tRapidity > fRapidity[1])) {
        fNTracksFailed++;
        return false;
    }
    if ((tEta < fEta[0]) || (tEta > fEta[1])) {
        fNTracksFailed++;
        return false;
    }

    // dowang 2.22
    if(fUsePtCut){
	    if ((tPt < fPt[0]) || (tPt > fPt[1])) {
		fNTracksFailed++;
		return false;
	    }
    }
    else{
	    if ((tTotalP < fPt[0]) || (tTotalP > fPt[1])) {
		fNTracksFailed++;
		return false;
	    }
    }

    if(fUsePtotalCut){
        if((tTotalP < MinPtotal) || ( tTotalP > MaxPtotal)) { 
		fNTracksFailed++;
		return false;
	}
    }

    // dowang for pion+he3 cut
    if(fPionHe3cut){
	if( track->TPCsignal() < MinHe3TPCSignal || track->TPCsignal() > MaxHe3TPCSignal){
		fNTracksFailed++;
		return false;	
	}
    }

    if ((track->PidProbElectron() < fPidProbElectron[0]) || (track->PidProbElectron() > fPidProbElectron[1])) {
        fNTracksFailed++;
        return false;
    }
    if ((track->PidProbPion() < fPidProbPion[0]) || (track->PidProbPion() > fPidProbPion[1])) {
        fNTracksFailed++;
        return false;
    }
    if ((track->PidProbKaon() < fPidProbKaon[0]) || (track->PidProbKaon() > fPidProbKaon[1])) {
        fNTracksFailed++;
        return false;
    }
    if ((track->PidProbProton() < fPidProbProton[0]) || (track->PidProbProton() > fPidProbProton[1])) {
        fNTracksFailed++;
        return false;
    }
    if ((track->PidProbMuon() < fPidProbMuon[0]) || (track->PidProbMuon() > fPidProbMuon[1])) {
        fNTracksFailed++;
        return false;
    }
    //****N Sigma Method -- electron rejection****
    if (fElectronRejection){
        if (!IsElectron(track->NSigmaTPCE(),track->NSigmaTPCPi(),track->NSigmaTPCK(), track->NSigmaTPCP())){
            return false;
	}
    }
    if (fMostProbable>0) {
        int imost=0;
	int loose_imost = 0;
        //****N Sigma Method****
        if (fPIDMethod==0) {

	float FillMom = 0.;
	if(fUsePtCut){
		FillMom = track->Pt();
	}
	else{
		FillMom = track->P().Mag();
	} 
	
	if(wiolaCrossCheck==1){
		FillMom = track->P().Mag();
	}
	    //\ for proton PID
            if (fMostProbable == 4) { // proton nsigma-PID required contour adjusting (in LHC10h)
                if (IsProtonNSigma(FillMom, track->NSigmaTPCP(), track->NSigmaTOFP(),SwitchMom_p) ) {
                    	imost = 4;
                }
		
		//\ reject
		if( fOtherNsigmacut==1 && RejectFakeP(track,FillMom, SwitchMom_p)){
			imost = 0;		
		}
		if( fOtherNsigmacut==2){
			if(Reject_commom(track,FillMom)
		|| IsDeuteronNSigma(FillMom, track->MassTOF(), fNsigmaMass, track->NSigmaTPCD(), track->NSigmaTOFD(), 1.4) ){
				imost = 0;
			}		
		}
		if ( fdEdxcut && !IsProtonTPCdEdx(track->Pt(), track->TPCsignal()) ){
                        imost = 0;
                }

		if(fUseDCAvsPt_cut){
			float tmpDCAr = TMath::Abs(track->ImpactD());
			float tmpCut = Return_DCAvsPt_cut_p(track->Pt(),fCharge);  
			if(  tmpDCAr > tmpCut ) imost = 0;		
		}
		// if all TOF cut(pt>0.3 add TOF) did not pass, also not pass this one
	if(fInversePID){
		if (IsProtonNSigma(FillMom, track->NSigmaTPCP(), track->NSigmaTOFP(),0.7)){
			loose_imost = 999;
			if( fOtherNsigmacut==1 && RejectFakeP(track,FillMom, 0.7)){
				loose_imost = 0;		
			}
		}
		// no pass strict PID but pass loose PID
		if(imost==0 && loose_imost==999){
			imost = 4;
		}else{
			imost = 0;
		}		
	}//end inversPID

		if(wiolaCrossCheck==1 && WiolaRejectPion(track->P().Mag(),track->NSigmaTPCPi(), track->NSigmaTOFPi()) ){
			imost = 0;
		}
            }
	    //\ for deuteron PID
            else if (fMostProbable == 13 && wiolaCrossCheck==0){   //cut on Nsigma deuteron

                if ( IsDeuteronNSigma(FillMom, track->MassTOF(), fNsigmaMass, track->NSigmaTPCD(), track->NSigmaTOFD(), SwitchMom_d)){
			imost = 13;
		}
		//\ reject 
                if( fOtherNsigmacut==1 && RejectFakeD(track,FillMom) ){
			imost = 0;			
		}
		if( fOtherNsigmacut==2){
			if(Reject_commom(track,track->P().Mag())
			|| IsProtonNSigma(FillMom, track->NSigmaTPCP(), track->NSigmaTOFP(),0.7)){
				imost = 0;
			}
		}
            
		if(fUseDCAvsPt_cut && fCharge>0){
			float tmpDCAr = TMath::Abs(track->ImpactD());
			float tmpCut = Return_DCAvsPt_cut_d(track->Pt(),fCharge);
			if(  tmpDCAr > tmpCut ) imost = 0;
			
		}
	if(fInversePID){
		if ( IsDeuteronNSigma(FillMom, track->MassTOF(), fNsigmaMass, track->NSigmaTPCD(), track->NSigmaTOFD(), 1.4)){
			loose_imost = 999;
		}
		// no pass strict PID but pass loose PID
		if(imost==0 && loose_imost==999){
			imost = 13;
		}else{
			imost = 0;
		}		
	}	
		
            }
	// 2023.2.8
	else if(fMostProbable == 13 && wiolaCrossCheck==1){

		if (WiolaDCut(track->P().Mag(),track->NSigmaTPCD(), track->NSigmaTOFD())) imost = 13;
	}
	    //\ for triton PID
            else if (fMostProbable == 14){   //cut on Nsigma triton

                if (IsTritonNSigma(FillMom, track->MassTOF(), fNsigmaMass, track->NSigmaTPCT(), track->NSigmaTOFT()) ){
                    imost = 14;
                }

                if ( fdEdxcut && !IsTritonTPCdEdx(track->Pt(), track->TPCsignal()) ){
                    imost = 0;
                }   
            }
 	    //\ for He3 PID
            else if (fMostProbable == 15){
                if (IsHe3NSigma(FillMom, track->MassTOF(), fNsigmaMass, track->NSigmaTPCH(), track->NSigmaTOFH()) ){

		        imost = 15;
                }
		if(fPionHe3cut && track->NSigmaTPCT() <= 0) imost = 0;
            }
	    else if (fMostProbable == 2){
		 
                if (IsPionNSigma(FillMom, track->MassTOF(), fNsigmaMass, track->NSigmaTPCPi(), track->NSigmaTOFPi()) ){
			imost = 2;
		}
	    }
            else if(fMostProbable == 16){
                if (IsElectronNSigma(FillMom, track->NSigmaTPCE(), track->NSigmaTOFE()) ){
                    imost = 16;
                }
            }
	    if(fRejectPions==1){
		if(FillMom  < fTPCThreshold){
			if(abs(track->NSigmaTOFPi()) < 3) imost = 0;
		}
		}
  
                if ( fdEdxcut &&fMostProbable == 13 && !IsDeuteronTPCdEdx(track->P().Mag(), track->TPCsignal()) ){
                    	imost = 0;
                }
//cout<<"xxx "<<imost<<" "<<fMostProbable<<endl;
	    if (imost != fMostProbable) return false;
	    if(fUseTOFMassCut){
		//Mass square!
		float TmpTOFMass = ReturnTOFMass(track,imost);
		if(TmpTOFMass == -999) return false;	
		//if(TmpTOFMass == -999){cout<<"imost"<<imost<<endl; return false;}
	    	if(TmpTOFMass < TOFMassLowLimit || TmpTOFMass > TOFMassUpLimit){
			return false;
		}
	    }
            
            
        }
        
        fNTracksPassed++ ;
        return true;
    
    }

    return false;

    
}
bool AliFemtoTrackCutPdtHe3::IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP, float tmp_switch){

    if (mom > tmp_switch) {
	        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < fNsigmaP)
	            return true;	
	}
    else {
        if (TMath::Abs(nsigmaTPCP) < fNsigmaP)
            return true;
    }

    return false;

}
bool AliFemtoTrackCutPdtHe3::IsDeuteronNSigma(float mom, float massTOFPDG, float sigmaMass, float nsigmaTPCD, float nsigmaTOFD,  float tmp_switch){
    //double massPDGD=1.8756;
	if(fdNSigmaVspTcut==0){
	    if (fNsigmaTPCTOF) {
		//Identyfication with only TPC for mom<1.4 and TPC&TOF for mom>1.4
		if (mom > tmp_switch){
		    if (TMath::Hypot( nsigmaTPCD, nsigmaTOFD ) < fNsigmaD)
		    //if ((TMath::Abs(nsigmaTPCD) < fNsigma) && (TMath::Abs(nsigmaTOFD) < fNsigma))
		        return true;
		}
		else{
		    if (TMath::Abs(nsigmaTPCD) < fNsigmaD)
		        return true;
		}
		
	    }
	    else{
		if (TMath::Abs(nsigmaTPCD) < fNsigmaD)
		    return true;
	    } 
	
	}
	else{
		float tmpcut = fdcutline_k * mom + fdcutline_b;

		if(tmp_switch > 4.0){	// only TPC
			if(nsigmaTPCD >= tmpcut && abs(nsigmaTPCD) < fNsigmaD){
				return true;
			}
		}
		else{	
			if(nsigmaTPCD >= tmpcut && abs(nsigmaTPCD) < fNsigmaD && TMath::Hypot( nsigmaTPCD, nsigmaTOFD ) < fNsigmaD){
				return true;
			}
		}
	}
    return false;
}
bool AliFemtoTrackCutPdtHe3::IsTritonNSigma(float mom, float massTOFPDG, float sigmaMass, float nsigmaTPCT, float nsigmaTOFT){
    //double massPDGD=2.8089;
	if(fdNSigmaVspTcut==0){
	    if (fNsigmaTPCTOF) {
		//Identyfication with only TPC for mom<1.4 and TPC&TOF for mom>1.4
		if (mom > SwitchMom_t){
		if(TMath::Hypot(nsigmaTPCT,nsigmaTOFT) < fNsigmaT )    
		//if ((TMath::Abs(nsigmaTPCT) < fNsigmaT) && (TMath::Abs(nsigmaTOFT) < fNsigmaT))
		        return true;
		}
		else{
		    if (TMath::Abs(nsigmaTPCT) < fNsigmaT)
		        return true;
		}
		
	    }
	    else{
		if (TMath::Abs(nsigmaTPCT) < fNsigmaT)
		    return true;
	    } 

	}else{
		float tmpcut = fdcutline_k * mom + fdcutline_b;
		if(SwitchMom_t > 4.0){	// only TPC
			//cout<<"sss "<<nsigmaTPCT<<" "<<fNsigmaT<<" "<<fdcutline_k<<" "<<fdcutline_b<<endl;
			if(nsigmaTPCT >= tmpcut && abs(nsigmaTPCT) < fNsigmaT){
				return true;
			}
		}
		else{	
			if(nsigmaTPCT >= tmpcut && abs(nsigmaTPCT) < fNsigmaT && TMath::Hypot( nsigmaTPCT, nsigmaTOFT ) < fNsigmaT){
				return true;
			}
		}
	}
    return false;
}
bool AliFemtoTrackCutPdtHe3::IsHe3NSigma(float mom, float massTOFPDG, float sigmaMass, float nsigmaTPCHe3, float nsigmaTOFHe3){
//	cout<<"IsHe3NSigma "<<nsigmaTPCHe3<<" "<<nsigmaTOFHe3<<endl;
    //double massPDGD=2.8089;
    if (fNsigmaTPCTOF) {
        if (mom > SwitchMom_He3){
            if ((TMath::Abs(nsigmaTPCHe3) < fNsigmaHe3) && (TMath::Abs(nsigmaTOFHe3) < fNsigmaHe3))
                return true;
        }
        else{
            if (TMath::Abs(nsigmaTPCHe3) < fNsigmaHe3)
                return true;
        }
    }
    else{
        if (TMath::Abs(nsigmaTPCHe3) < fNsigmaHe3)
            return true;
    }
    return false;
  
}
bool AliFemtoTrackCutPdtHe3::IsPionNSigma(float mom, float massTOFPDG, float sigmaMass, float nsigmaTPCPi, float nsigmaTOFPi){
    
    if (fNsigmaTPCTOF) {
        if (mom > SwitchMom_Pi){
	    if (TMath::Hypot( nsigmaTPCPi, nsigmaTOFPi ) < fNsigmaPi)
                return true;
        }
        else{
            if (TMath::Abs(nsigmaTPCPi) < fNsigmaPi)
                return true;
        }
    }
    else{
        if (TMath::Abs(nsigmaTPCPi) < fNsigmaPi)
            return true;
    }
    return false;
  
}
//rejection methods
bool AliFemtoTrackCutPdtHe3::IsElectronNSigmaRejection(float mom, float nsigmaTPCE){
  if(TMath::Abs(nsigmaTPCE) < fNsigmaRejection)
    return true;
 

  return false;
}

bool AliFemtoTrackCutPdtHe3::IsPionNSigmaRejection(float mom, float nsigmaTPCPi, float nsigmaTOFPi){
  if(fOnlyTPCreject){
    if(TMath::Abs(nsigmaTPCPi) < fNsigmaRejection)
      return true;
  }
  else{
    if(TMath::Hypot( nsigmaTOFPi, nsigmaTPCPi ) < fNsigmaRejection)
      return true;
  }

  return false;
}


bool AliFemtoTrackCutPdtHe3::IsKaonNSigmaRejection(float mom, float nsigmaTPCK, float nsigmaTOFK){
  if(fOnlyTPCreject){
    if(TMath::Abs(nsigmaTPCK) < fNsigmaRejection)
      return true;
  }
  else{
    if(TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < fNsigmaRejection)
      return true;
  }

  return false;
}

bool AliFemtoTrackCutPdtHe3::IsProtonNSigmaRejection(float mom, float nsigmaTPCP, float nsigmaTOFP){
  if(fOnlyTPCreject){
    if(TMath::Abs(nsigmaTPCP) < fNsigmaRejection)
      return true;
  }
  else{
    if(TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < fNsigmaRejection)
      return true;
  }
 
  return false;
}
void AliFemtoTrackCutPdtHe3::SetProtonSwitchMom(float SwitchMom){
	SwitchMom_p = SwitchMom;
}
void AliFemtoTrackCutPdtHe3::SetDeuteronSwitchMom(float SwitchMom){
	SwitchMom_d = SwitchMom;
}
void AliFemtoTrackCutPdtHe3::SetTritonSwitchMom(float SwitchMom){
	SwitchMom_t = SwitchMom;
}
void AliFemtoTrackCutPdtHe3::SetHe3SwitchMom(float SwitchMom){
	SwitchMom_He3 = SwitchMom;
}
void AliFemtoTrackCutPdtHe3::SetPionSwitchMom(float SwitchMom){
	SwitchMom_Pi = SwitchMom;
}
//\ for e+e femto
bool AliFemtoTrackCutPdtHe3::IsElectronNSigma(float mom, float nsigmaTPCE, float nsigmaTOFE){

    if(TMath::Abs(nsigmaTPCE) < 3.) return true;
    return false;
}
void AliFemtoTrackCutPdtHe3::SetPionHe3Cut(int aPionHe3cut){
	fPionHe3cut = aPionHe3cut;
}
bool AliFemtoTrackCutPdtHe3::IsProtonTPCdEdx(float mom, float dEdx){
    double a1 = -250.,  b1 = 400.;
    if (mom < 1.) {
        if (dEdx > a1*mom+b1) return false;
    }
    return true;

}
//\ follow wiola
bool AliFemtoTrackCutPdtHe3::IsDeuteronTPCdEdx(float mom, float dEdx){
	double a1 = -250.0,  b1 = 400.0;
  	double a2 = -80,   b2 = 190.0;
  	if (mom < 1.1) {
    		if (dEdx < a1 * mom+b1) return false;
  	}
  	else if (mom < 2) {
    		if (dEdx < a2 *mom+b2) return false;
  	}
    	return true;

}
bool AliFemtoTrackCutPdtHe3::IsTritonTPCdEdx(float mom, float dEdx){
    
        if (dEdx < MinHe3TPCSignal) return false;
    
 
    	return true;

}

void AliFemtoTrackCutPdtHe3::SetdEdxcutLabel(int dEdxcutLabel){
    fdEdxcut = dEdxcutLabel;
}
void AliFemtoTrackCutPdtHe3::SetProtonNsigma(float Nsigma){
    fNsigmaP = Nsigma;
}
void AliFemtoTrackCutPdtHe3::SetDeuteronNsigma(float Nsigma){
    fNsigmaD = Nsigma;
}
void AliFemtoTrackCutPdtHe3::SetTritonNsigma(float Nsigma){
    fNsigmaT = Nsigma;
}
void AliFemtoTrackCutPdtHe3::SetHe3Nsigma(float Nsigma){
    fNsigmaHe3 = Nsigma;
}
void AliFemtoTrackCutPdtHe3::SetPionNsigma(float Nsigma){
    fNsigmaPi = Nsigma;
}
void AliFemtoTrackCutPdtHe3::SetRejectionNsigma(float Nsigma){
    fNsigmaRejection = Nsigma;
}
void AliFemtoTrackCutPdtHe3::SetOtherNsigmacutLabel(int OtherNsigmaLabel){
    fOtherNsigmacut = OtherNsigmaLabel;
}
void AliFemtoTrackCutPdtHe3::SetMinTPCNCrossRows(int MinTPCNCrossRows){
    fMinTPCNCrossRows = MinTPCNCrossRows;
}
void AliFemtoTrackCutPdtHe3::SetMinTPCFoundFraction(float MinTPCFoundFraction){
    fMinTPCFoundFraction = MinTPCFoundFraction;
}
void AliFemtoTrackCutPdtHe3::SetUseTOFMassCut(int UseTOFMassCut){
    fUseTOFMassCut = UseTOFMassCut;
}
void AliFemtoTrackCutPdtHe3::SetTOFMassLimit(float LowMass,float UpMass){
    TOFMassLowLimit = LowMass;
    TOFMassUpLimit  = UpMass;  
}
float AliFemtoTrackCutPdtHe3::ReturnTOFMass(const AliFemtoTrack* track,int imost){
	float tMom = track->P().Mag();
	float c=1.;
	float beta = track->VTOF();
	float massTOF = -1.;//Mass square

	float imostMass = 0.;
	if(imost==4) 	imostMass = 0.938272;
	if(imost==13)	imostMass = 1.8756;
	if(imost==14)	imostMass = 2.8089;
	if(imost==15)	imostMass = 2.8084;
	if(imost==2)	imostMass = 0.13957018;
	if(beta!=0){
		massTOF = tMom*tMom/c/c*(1/(beta*beta)-1);  
		return  massTOF - imostMass*imostMass;
	}
	else{
		return -999;	
	}
	

}
void AliFemtoTrackCutPdtHe3::SetfUsePtCut(int aUsePtCut){
	fUsePtCut = aUsePtCut;
}
void AliFemtoTrackCutPdtHe3::SetfOnlyTPCreject(int aOnlyTPCreject){
	fOnlyTPCreject = aOnlyTPCreject;
}
bool AliFemtoTrackCutPdtHe3::RejectFakeP(const AliFemtoTrack* track, float mom, float tmp_switch){

	bool rejected = true;
	float p_NsigmaCombine = 0.;
	float k_NsigmaCombine = 0.;
	float e_NsigmaCombine = 0.;
	float pi_NsigmaCombine = 0.;

	if (mom > tmp_switch) {
		p_NsigmaCombine = TMath::Hypot(track->NSigmaTPCP(), track->NSigmaTOFP());
		k_NsigmaCombine = TMath::Hypot(track->NSigmaTPCK(), track->NSigmaTOFK());
		e_NsigmaCombine = TMath::Hypot(track->NSigmaTPCE(), track->NSigmaTOFE());
		pi_NsigmaCombine = TMath::Hypot(track->NSigmaTPCPi(), track->NSigmaTOFPi());
	        if ( 	(k_NsigmaCombine < p_NsigmaCombine) || 
			(e_NsigmaCombine < p_NsigmaCombine) || 
			(pi_NsigmaCombine < p_NsigmaCombine) )
	            return rejected;	
	}
    	else {
		p_NsigmaCombine = abs(track->NSigmaTPCP());
		k_NsigmaCombine = abs(track->NSigmaTPCK());
		e_NsigmaCombine = abs(track->NSigmaTPCE());
		pi_NsigmaCombine = abs(track->NSigmaTPCPi());
	        if ( 	(k_NsigmaCombine < p_NsigmaCombine) || 
			(e_NsigmaCombine < p_NsigmaCombine) || 
			(pi_NsigmaCombine < p_NsigmaCombine) )
	            return rejected;
    	}
	return false;

}
bool AliFemtoTrackCutPdtHe3::RejectFakeD(const AliFemtoTrack* track, float mom){

	bool rejected = true;
	float d_NsigmaCombine = 0.;
	float p_NsigmaCombine = 0.;
	float k_NsigmaCombine = 0.;
	//float e_NsigmaCombine = 0.;
	float pi_NsigmaCombine = 0.;
	if (mom > SwitchMom_d) {
		d_NsigmaCombine = TMath::Hypot(track->NSigmaTPCD(), track->NSigmaTOFD());
		p_NsigmaCombine = TMath::Hypot(track->NSigmaTPCP(), track->NSigmaTOFP());
		k_NsigmaCombine = TMath::Hypot(track->NSigmaTPCK(), track->NSigmaTOFK());
		//e_NsigmaCombine = TMath::Hypot(track->NSigmaTPCE(), track->NSigmaTOFE());
		pi_NsigmaCombine = TMath::Hypot(track->NSigmaTPCPi(), track->NSigmaTOFPi());
	        if ( 	(k_NsigmaCombine 	< d_NsigmaCombine) || 
			//(e_NsigmaCombine 	< d_NsigmaCombine) || 
			(pi_NsigmaCombine 	< d_NsigmaCombine) ||
			(p_NsigmaCombine 	< d_NsigmaCombine) )
	            return rejected;	
	}
    	else {
		d_NsigmaCombine = abs(track->NSigmaTPCD());      	
		p_NsigmaCombine = abs(track->NSigmaTPCP());
		k_NsigmaCombine = abs(track->NSigmaTPCK());
		//e_NsigmaCombine = abs(track->NSigmaTPCE());
		pi_NsigmaCombine = abs(track->NSigmaTPCPi());
	        if ( 	(k_NsigmaCombine 	< d_NsigmaCombine) || 
			//(e_NsigmaCombine 	< d_NsigmaCombine) || 
			(pi_NsigmaCombine 	< d_NsigmaCombine) ||
			(p_NsigmaCombine 	< d_NsigmaCombine) )
	            return rejected;
    	}
	return false;

}
bool AliFemtoTrackCutPdtHe3::IsPionNSigma(float mom,float nsigmaTPCpi,float nsigmaTOFpi){
	
	if (mom > 0.5) {
	        if (TMath::Hypot( nsigmaTOFpi, nsigmaTPCpi ) < 3)
	            return true;	
	}
    else {
        if (TMath::Abs(nsigmaTPCpi) < 3.)
            return true;
    }
    return false;
}
bool AliFemtoTrackCutPdtHe3::IsKaonNSigma(float mom,float nsigmaTPCk,float nsigmaTOFk){
	
	if (mom > 0.5) {
	        if (TMath::Hypot( nsigmaTPCk, nsigmaTOFk ) < 3)
	            return true;	
	}
    else {
        if (TMath::Abs(nsigmaTPCk) < 3.)
            return true;
    }
    return false;
}
bool AliFemtoTrackCutPdtHe3::Reject_commom(const AliFemtoTrack* track, float mom){
	bool rejected = true;
	bool passpi = IsPionNSigma(mom,track->NSigmaTPCPi(), track->NSigmaTOFPi());
	bool passK = IsKaonNSigma(mom,track->NSigmaTPCK(), track->NSigmaTOFK());
	if(passpi || passK){ 
		return rejected;
	}else{
		return false;
	}

}
float AliFemtoTrackCutPdtHe3::SetfUseDCAvsPt_cut(int aUseDCAvsPt_cut){
	fUseDCAvsPt_cut = aUseDCAvsPt_cut; 
}
float AliFemtoTrackCutPdtHe3::Return_DCAvsPt_cut_p(float pt,int charge){
	int WhichBin = (pt - 0.4)/0.2;
	if(charge>0){
		return v_DCAvspTcut_p[WhichBin];		
	} 
	if(charge<0){
		return v_DCAvspTcut_antip[WhichBin];	
	}

}
float AliFemtoTrackCutPdtHe3::Return_DCAvsPt_cut_d(float pt,int charge){
	int WhichBin = (pt - 0.4)/0.2;
	return v_DCAvspTcut_d[WhichBin];	

}
void AliFemtoTrackCutPdtHe3::Set_DCAvsPt_cut(float *input_v,int label){
	if(label==0){
		for(int i=0;i<18;i++) v_DCAvspTcut_p[i] = input_v[i];
	}
	if(label==1){
		for(int i=0;i<18;i++) v_DCAvspTcut_antip[i] = input_v[i];			
	}
	if(label==2){
		for(int i=0;i<18;i++) v_DCAvspTcut_d[i] = input_v[i];			
	}
}
void AliFemtoTrackCutPdtHe3::SetInversePID(int aInversePID){
	fInversePID = aInversePID;
}

void AliFemtoTrackCutPdtHe3::SetHe3TPCSignal(float aMin,float aMax){
    MinHe3TPCSignal = aMin;
    MaxHe3TPCSignal = aMax;
}
void AliFemtoTrackCutPdtHe3::SetfUseRejectLowpTPion(int aUse){
	fRejectPions = aUse;
}
void AliFemtoTrackCutPdtHe3::SetfTPCThreshold(float aTPCThreshold){ 
	fTPCThreshold = aTPCThreshold;
}
void AliFemtoTrackCutPdtHe3::SetfdNSigmaVspTcut(int aUse){
	fdNSigmaVspTcut = aUse;
}
void AliFemtoTrackCutPdtHe3::Setfdline(float ak,float ab){
	fdcutline_k = ak;
	fdcutline_b = ab;
}
void AliFemtoTrackCutPdtHe3::Setf2pT(int aUse){
	fUse2pT = aUse;
}
void AliFemtoTrackCutPdtHe3::SetAlldEdxMode(int aUse){
	AlldEdxmode = aUse;
}
void AliFemtoTrackCutPdtHe3::SetwiolaCrossCheck(int aUse){
	wiolaCrossCheck = aUse;
}
bool AliFemtoTrackCutPdtHe3::WiolaDCut(float mom, float nsigmaTPCD, float nsigmaTOFD){

	if(mom < 1.3){
		if(abs(nsigmaTPCD) < 2)	return true;
	}
	else{
		if((abs(nsigmaTPCD) < 2) && (abs(nsigmaTOFD) < 2)) return true;
	}
	return false;
}
bool AliFemtoTrackCutPdtHe3::WiolaRejectPion(float mom,float nsigmaTPCpi,float nsigmaTOFpi){
if (mom > 0.5) {
	        if (TMath::Hypot( nsigmaTOFpi, nsigmaTPCpi ) < pionrejectcut) return true;	
	}
    	else {
        	if (TMath::Abs(nsigmaTPCpi) < pionrejectcut) return true;
    	}
    	return false;

}
void AliFemtoTrackCutPdtHe3::SetUsePtotal(int aUse){
fUsePtotalCut = aUse;
}
void AliFemtoTrackCutPdtHe3::SetPtotalRange(float aMin,float aMax){
	 MinPtotal = aMin;
         MaxPtotal = aMax;
}
void AliFemtoTrackCutPdtHe3::Setpionrejectcut(float aRejectCut){
pionrejectcut = aRejectCut;
}
