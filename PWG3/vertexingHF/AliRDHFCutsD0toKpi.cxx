/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// Class for cuts on AOD reconstructed D0->Kpi
//
// Author: A.Dainese, andrea.dainese@pd.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <Riostream.h>

#include "AliRDHFCutsD0toKpi.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAODPid.h"
#include "AliTPCPIDResponse.h"
#include "AliAODVertex.h"

ClassImp(AliRDHFCutsD0toKpi)

//--------------------------------------------------------------------------
AliRDHFCutsD0toKpi::AliRDHFCutsD0toKpi(const char* name) : 
AliRDHFCuts(name),
fUseSpecialCuts(kFALSE),
fLowPt(kTRUE),
fDefaultPID(kTRUE)
{
  //
  // Default Constructor
  //
  Int_t nvars=9;
  SetNVars(nvars);
  TString varNames[9]={"inv. mass [GeV]",   
		       "dca [cm]",
		       "cosThetaStar", 
		       "pTK [GeV/c]",
		       "pTPi [GeV/c]",
		       "d0K [cm]",
		       "d0Pi [cm]",
		       "d0d0 [cm^2]",
		       "cosThetaPoint"};
  Bool_t isUpperCut[9]={kTRUE,
			kTRUE,
			kTRUE,
			kFALSE,
			kFALSE,
			kTRUE,
			kTRUE,
			kTRUE,
			kFALSE};
  SetVarNames(nvars,varNames,isUpperCut);
  Bool_t forOpt[9]={kFALSE,
		    kTRUE,
		    kTRUE,
		    kFALSE,
		    kFALSE,
		    kFALSE,
		    kFALSE,
		    kTRUE,
		    kTRUE};
  SetVarsForOpt(4,forOpt);
  Float_t limits[2]={0,999999999.};
  SetPtBins(2,limits);

}
//--------------------------------------------------------------------------
AliRDHFCutsD0toKpi::AliRDHFCutsD0toKpi(const AliRDHFCutsD0toKpi &source) :
  AliRDHFCuts(source),   
  fUseSpecialCuts(source.fUseSpecialCuts),
  fLowPt(source.fLowPt),
  fDefaultPID(source.fDefaultPID)
{
  //
  // Copy constructor
  //

}
//--------------------------------------------------------------------------
AliRDHFCutsD0toKpi &AliRDHFCutsD0toKpi::operator=(const AliRDHFCutsD0toKpi &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliRDHFCuts::operator=(source); 
  fUseSpecialCuts=source.fUseSpecialCuts;
  fLowPt=source.fLowPt;
  fDefaultPID=source.fDefaultPID;

  return *this;
}


//---------------------------------------------------------------------------
void AliRDHFCutsD0toKpi::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters) {
  // 
  // Fills in vars the values of the variables 
  //

  if(nvars!=fnVarsForOpt) {
    printf("AliRDHFCutsD0toKpi::GetCutsVarsForOpt: wrong number of variables\n");
    return;
  }

  AliAODRecoDecayHF2Prong *dd = (AliAODRecoDecayHF2Prong*)d;
 
  Int_t iter=-1;
  if(fVarsForOpt[0]){
    iter++;
    if(TMath::Abs(pdgdaughters[0])==211) {
      vars[iter]=dd->InvMassD0();
    } else {
      vars[iter]=dd->InvMassD0bar();
    }
  }
  if(fVarsForOpt[1]){
    iter++;
    vars[iter]=dd->GetDCA();
  }
  if(fVarsForOpt[2]){
    iter++;
    if(TMath::Abs(pdgdaughters[0])==211) {
      vars[iter] = dd->CosThetaStarD0();
    } else {
      vars[iter] = dd->CosThetaStarD0bar();
    }
  }
  if(fVarsForOpt[3]){
    iter++;
   if(TMath::Abs(pdgdaughters[0])==321) {
     vars[iter]=dd->PtProng(0);
   }
   else{
     vars[iter]=dd->PtProng(1);
   }
  }
  if(fVarsForOpt[4]){
    iter++;
   if(TMath::Abs(pdgdaughters[0])==211) {
     vars[iter]=dd->PtProng(0);
   }
   else{
     vars[iter]=dd->PtProng(1);
   }
  }
  if(fVarsForOpt[5]){
    iter++;
    if(TMath::Abs(pdgdaughters[0])==321) {
     vars[iter]=dd->Getd0Prong(0);
   }
   else{
     vars[iter]=dd->Getd0Prong(1);
   }
  }
  if(fVarsForOpt[6]){
    iter++;
     if(TMath::Abs(pdgdaughters[0])==211) {
     vars[iter]=dd->Getd0Prong(0);
   }
   else{
     vars[iter]=dd->Getd0Prong(1);
   }
  }
  if(fVarsForOpt[7]){
    iter++;
    vars[iter]= dd->Prodd0d0();
  }
  if(fVarsForOpt[8]){
    iter++;
    vars[iter]=dd->CosPointingAngle();
  }
  
  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpi::IsSelected(TObject* obj,Int_t selectionLevel,AliAODEvent* aod) {
  //
  // Apply selection
  //

  if(!fCutsRD){
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }
  //PrintAll();
  AliAODRecoDecayHF2Prong* d=(AliAODRecoDecayHF2Prong*)obj;

  if(!d){
    cout<<"AliAODRecoDecayHF2Prong null"<<endl;
    return 0;
  }

  // selection on daughter tracks 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kTracks) {
    if(!AreDaughtersSelected(d)) return 0;
  }


 
  // returnvalue: 0 not sel, 1 only D0, 2 only D0bar, 3 both
  Int_t returnvaluePID=3;
  Int_t returnvalueCuts=3;

  // selection on PID 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate ||
     selectionLevel==AliRDHFCuts::kPID) {
    returnvaluePID = IsSelectedPID(d);
  }



  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {

    //recalculate vertex w/o daughters
    AliAODVertex *origownvtx=0x0;
    AliAODVertex *recvtx=0x0;
  
    if(fRemoveDaughtersFromPrimary) {
      if(!aod) {
	AliError("Can not remove daughters from vertex without AOD event");
	return 0;
      }   
      if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());
      recvtx=d->RemoveDaughtersFromPrimaryVtx(aod);
      if(!recvtx){
	AliDebug(2,"Removal of daughter tracks failed");
	//recvtx=d->GetPrimaryVtx();
	if(origownvtx){
	  delete origownvtx;
	  origownvtx=NULL;
	}
	return 0;
      }
      //set recalculed primary vertex
      d->SetOwnPrimaryVtx(recvtx);
      delete recvtx; recvtx=NULL;
    }

    
    Double_t pt=d->Pt();
   
    Int_t okD0=0,okD0bar=0;
 
    Int_t ptbin=PtBin(pt);
    if (ptbin==-1) {
      if(origownvtx){
	d->SetOwnPrimaryVtx(origownvtx);
	delete origownvtx;
	origownvtx=NULL;
      }
      else d->UnsetOwnPrimaryVtx();
      return 0;
    }
    Double_t mD0,mD0bar,ctsD0,ctsD0bar;
    okD0=1; okD0bar=1;

    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();

    if(d->PtProng(1) < fCutsRD[GetGlobalIndex(3,ptbin)] || d->PtProng(0) < fCutsRD[GetGlobalIndex(4,ptbin)]) okD0 = 0;
    if(d->PtProng(0) < fCutsRD[GetGlobalIndex(3,ptbin)] || d->PtProng(1) < fCutsRD[GetGlobalIndex(4,ptbin)]) okD0bar = 0;
    if(!okD0 && !okD0bar) returnvalueCuts=0;
  
    
    if(TMath::Abs(d->Getd0Prong(1)) > fCutsRD[GetGlobalIndex(5,ptbin)] || 
       TMath::Abs(d->Getd0Prong(0)) > fCutsRD[GetGlobalIndex(6,ptbin)]) okD0 = 0;
    if(TMath::Abs(d->Getd0Prong(0)) > fCutsRD[GetGlobalIndex(6,ptbin)] ||
       TMath::Abs(d->Getd0Prong(1)) > fCutsRD[GetGlobalIndex(5,ptbin)]) okD0bar = 0;
    if(!okD0 && !okD0bar)  returnvalueCuts=0;
    
    if(d->GetDCA() > fCutsRD[GetGlobalIndex(1,ptbin)])  returnvalueCuts=0;
    
    d->InvMassD0(mD0,mD0bar);
    if(TMath::Abs(mD0-mD0PDG) > fCutsRD[GetGlobalIndex(0,ptbin)]) okD0 = 0;
    if(TMath::Abs(mD0bar-mD0PDG) > fCutsRD[GetGlobalIndex(0,ptbin)])  okD0bar = 0;
    if(!okD0 && !okD0bar)  returnvalueCuts=0;
    
    d->CosThetaStarD0(ctsD0,ctsD0bar);
    if(TMath::Abs(ctsD0) > fCutsRD[GetGlobalIndex(2,ptbin)]) okD0 = 0; 
    if(TMath::Abs(ctsD0bar) > fCutsRD[GetGlobalIndex(2,ptbin)]) okD0bar = 0;
    if(!okD0 && !okD0bar)   returnvalueCuts=0;
    
    if(d->Prodd0d0() > fCutsRD[GetGlobalIndex(7,ptbin)])  returnvalueCuts=0;
    
    if(d->CosPointingAngle() < fCutsRD[GetGlobalIndex(8,ptbin)])  returnvalueCuts=0;
    
    if (returnvalueCuts!=0) {
      if (okD0) returnvalueCuts=1; //cuts passed as D0
      if (okD0bar) returnvalueCuts=2; //cuts passed as D0bar
      if (okD0 && okD0bar) returnvalueCuts=3; //cuts passed as D0 and D0bar
    }

    // call special cuts
    Int_t special=1;
    if(fUseSpecialCuts) special=IsSelectedSpecialCuts(d);
    if(!special) returnvalueCuts=0;

    // unset recalculated primary vertex when not needed any more
    if(origownvtx) {
      d->SetOwnPrimaryVtx(origownvtx);
      delete origownvtx;
      origownvtx=NULL;
    } else if(fRemoveDaughtersFromPrimary) {
      d->UnsetOwnPrimaryVtx();
      AliDebug(3,"delete new vertex\n");
    }

  }

 

  //  cout<<"Pid = "<<returnvaluePID<<endl;
  return CombineSelectionLevels(3,returnvalueCuts,returnvaluePID);
}
//---------------------------------------------------------------------------
Bool_t AliRDHFCutsD0toKpi::IsInFiducialAcceptance(Double_t pt, Double_t y) const
{
  //
  // Checking if D0 is in fiducial acceptance region 
  //

  if(pt > 5.) {
    // applying cut for pt > 5 GeV
    AliDebug(2,Form("pt of D0 = %f (> 5), cutting at |y| < 0.8\n",pt)); 
    if (TMath::Abs(y) > 0.8){
      return kFALSE;
    }
  } else {
    // appliying smooth cut for pt < 5 GeV
    Double_t maxFiducialY = -0.2/15*pt*pt+1.9/15*pt+0.5; 
    Double_t minFiducialY = 0.2/15*pt*pt-1.9/15*pt-0.5;		
    AliDebug(2,Form("pt of D0 = %f (< 5), cutting  according to the fiducial zone [%f, %f]\n",pt,minFiducialY,maxFiducialY)); 
    if (y < minFiducialY || y > maxFiducialY){
      return kFALSE;
    }
  }

  return kTRUE;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpi::IsSelectedPID(AliAODRecoDecayHF* d) 
{
  // ############################################################
  //
  // Apply PID selection
  //
  //
  // ############################################################

  if(!fUsePID) return 3;
  if(fDefaultPID) return IsSelectedPIDdefault(d);
  fWhyRejection=0;
  Int_t isD0D0barPID[2]={1,2};
  Int_t combinedPID[2][2];// CONVENTION: [daught][isK,IsPi]; [0][0]=(prong 1, isK)=value [0][1]=(prong 1, isPi)=value; 
  //                                                                                                 same for prong 2
  //                                               values convention -1 = discarded 
  //                                                                  0 = not identified (but compatible) || No PID (->hasPID flag)
  //                                                                  1 = identified
  // PID search:   pion (TPC) or not K (TOF), Kaon hypothesis for both 
  // Initial hypothesis: unknwon (but compatible) 
  combinedPID[0][0]=0;  // prima figlia, Kaon
  combinedPID[0][1]=0;  // prima figlia, pione
  combinedPID[1][0]=0;  // seconda figlia, Kaon
  combinedPID[1][1]=0;  // seconda figlia, pion
  
  Bool_t checkPIDInfo[2]={kTRUE,kTRUE};
  Double_t sigma_tmp[3]={fPidHF->GetSigma(0),fPidHF->GetSigma(1),fPidHF->GetSigma(2)};
  for(Int_t daught=0;daught<2;daught++){
    //Loop con prongs
    AliAODTrack *aodtrack=(AliAODTrack*)d->GetDaughter(daught);
    
    if(!(fPidHF->CheckStatus(aodtrack,"TPC")) && !(fPidHF->CheckStatus(aodtrack,"TOF"))) {
    checkPIDInfo[daught]=kFALSE; 
    continue;
    }

    // identify kaon
    combinedPID[daught][0]=fPidHF->MakeRawPid(aodtrack,3);

    // identify pion

    if(!(fPidHF->CheckStatus(aodtrack,"TPC"))) {
     combinedPID[daught][1]=0;
    }else{
      fPidHF->SetTOF(kFALSE);
      combinedPID[daught][1]=fPidHF->MakeRawPid(aodtrack,2);
      fPidHF->SetTOF(kTRUE);
      fPidHF->SetCompat(kTRUE);
     }


   if(combinedPID[daught][0]<=-1&&combinedPID[daught][1]<=-1){ // if not a K- and not a pi- both D0 and D0bar excluded
    isD0D0barPID[0]=0;
    isD0D0barPID[1]=0;
   }
   else if(combinedPID[daught][0]==2&&combinedPID[daught][1]>=1){
    if(aodtrack->Charge()==-1)isD0D0barPID[1]=0;//if K- D0bar excluded
    else isD0D0barPID[0]=0;// if K+ D0 excluded
   }
   else if(combinedPID[daught][0]==1&&combinedPID[daught][1]>=1){
    isD0D0barPID[0]=0;
    isD0D0barPID[1]=0;
   }
   else if(combinedPID[daught][0]>=1||combinedPID[daught][1]<=-1){ 
   if(aodtrack->Charge()==-1)isD0D0barPID[1]=0;// not a D0bar if K- or if pi- excluded
   else isD0D0barPID[0]=0;//  not a D0 if K+ or if pi+ excluded
        }
   else if(combinedPID[daught][0]<=-1||combinedPID[daught][1]>=1){
    if(aodtrack->Charge()==-1)isD0D0barPID[0]=0;// not a D0 if pi- or if K- excluded
    else isD0D0barPID[1]=0;// not a D0bar if pi+ or if K+ excluded
   }

    if(fLowPt && d->Pt()<2.){
     Double_t sigmaTPC[3]={3.,2.,0.};
     fPidHF->SetSigmaForTPC(sigmaTPC);
    // identify kaon
    combinedPID[daught][0]=fPidHF->MakeRawPid(aodtrack,3);

    Double_t ptProng=aodtrack->P();

    if(ptProng<0.6){
     fPidHF->SetCompat(kFALSE);
     combinedPID[daught][0]=fPidHF->MakeRawPid(aodtrack,3);
     fPidHF->SetCompat(kTRUE);
    }

    if(!(fPidHF->CheckStatus(aodtrack,"TPC"))) {
     combinedPID[daught][1]=0;
    }else{
      fPidHF->SetTOF(kFALSE);
      Double_t sigmaTPCpi[3]={3.,3.,0.};
      fPidHF->SetSigmaForTPC(sigmaTPCpi);
      combinedPID[daught][1]=fPidHF->MakeRawPid(aodtrack,2);
      fPidHF->SetTOF(kTRUE);
       if(ptProng<0.8){
        Bool_t isTPCpion=fPidHF->IsPionRaw(aodtrack,"TPC");
        if(isTPCpion){
         combinedPID[daught][1]=1;
        }else{
         combinedPID[daught][1]=-1;
        }
      }
    }

   }
   fPidHF->SetSigmaForTPC(sigma_tmp);
  }// END OF LOOP ON DAUGHTERS

   if(!checkPIDInfo[0] && !checkPIDInfo[1]) {
    if(fLowPt) fPidHF->SetSigmaForTPC(sigma_tmp);
    return 0;
   }


  // FURTHER PID REQUEST (both daughter info is needed)
  if(combinedPID[0][0]<=-1&&combinedPID[1][0]<=-1){
    fWhyRejection=31;// reject cases in which no kaon-compatible tracks are found
    if(fLowPt) fPidHF->SetSigmaForTPC(sigma_tmp);
    return 0;
  }

  if(d->Pt()<2.){
    if(fLowPt) fPidHF->SetSigmaForTPC(sigma_tmp);
    if(combinedPID[0][0]<=0&&combinedPID[1][0]<=0){
      fWhyRejection=32;// reject cases where the Kaon is not identified
      return 0;
    }
  }
    if(fLowPt) fPidHF->SetSigmaForTPC(sigma_tmp);

  //  cout<<"Why? "<<fWhyRejection<<endl;  
  return isD0D0barPID[0]+isD0D0barPID[1];
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpi::IsSelectedPIDdefault(AliAODRecoDecayHF* d) 
{
  // ############################################################
  //
  // Apply PID selection
  //
  //
  // temporary selection: PID AS USED FOR D0 by Andrea Rossi (up to 28/06/2010)
  //
  // d must be a AliAODRecoDecayHF2Prong object
  // returns 0 if both D0 and D0bar are rejectecd
  //         1 if D0 is accepted while D0bar is rejected
  //         2 if D0bar is accepted while D0 is rejected
  //         3 if both are accepted
  // fWhyRejection variable (not returned for the moment, print it if needed)
  //               keeps some information on why a candidate has been 
  //               rejected according to the following (unfriendly?) scheme 
  //             if more rejection cases are considered interesting, just add numbers
  //
  //      TO BE CONSIDERED WITH A GRAIN OF SALT (the order in which cut are applied is relevant) 
  //              from 20 to 30: "detector" selection (PID acceptance)                                             
  //                                                  26: TPC refit
  //                                                  27: ITS refit
  //                                                  28: no (TOF||TPC) pid information (no kTOFpid,kTOFout,kTIME,kTPCpid,...)
  //
  //              from 30 to 40: PID selection
  //                                                  31: no Kaon compatible tracks found between daughters
  //                                                  32: no Kaon identified tracks found (strong sel. at low momenta)
  //                                                  33: both mass hypotheses are rejected 
  //                  
  // ############################################################

  if(!fUsePID) return 3;
  fWhyRejection=0;
  Int_t isD0D0barPID[2]={1,2};
  Double_t nsigmaTPCpi=-1., nsigmaTPCK=-1.; //used for TPC pid
  Double_t tofSig,times[5];// used fot TOF pid
  Int_t hasPID[2]={2,2};// flag to count how many detectors give PID info for the daughters
  Int_t isKaonPionTOF[2][2],isKaonPionTPC[2][2];
  Int_t combinedPID[2][2];// CONVENTION: [daught][isK,IsPi]; [0][0]=(prong 1, isK)=value [0][1]=(prong 1, isPi)=value; 
  //                                                                                                 same for prong 2
  //                                               values convention -1 = discarded 
  //                                                                  0 = not identified (but compatible) || No PID (->hasPID flag)
  //                                                                  1 = identified
  // PID search:   pion (TPC) or not K (TOF), Kaon hypothesis for both 
  // Initial hypothesis: unknwon (but compatible) 
  isKaonPionTOF[0][0]=0;
  isKaonPionTOF[0][1]=0;
  isKaonPionTOF[1][0]=0;
  isKaonPionTOF[1][1]=0;
  
  isKaonPionTPC[0][0]=0;
  isKaonPionTPC[0][1]=0;
  isKaonPionTPC[1][0]=0;
  isKaonPionTPC[1][1]=0;
  
  combinedPID[0][0]=0;
  combinedPID[0][1]=0;
  combinedPID[1][0]=0;
  combinedPID[1][1]=0;
  
  
 
  for(Int_t daught=0;daught<2;daught++){
    //Loop con prongs
    
    // ########### Step 0- CHECKING minimal PID "ACCEPTANCE" ####################

    AliAODTrack *aodtrack=(AliAODTrack*)d->GetDaughter(daught); 
   
    if(!(aodtrack->GetStatus()&AliESDtrack::kTPCrefit)){
      fWhyRejection=26;
      return 0;
    } 
    if(!(aodtrack->GetStatus()&AliESDtrack::kITSrefit)){
      fWhyRejection=27;
      return 0;
    } 
    
    AliAODPid *pid=aodtrack->GetDetPid();
    if(!pid) {
      //delete esdtrack;
      hasPID[daught]--;
      continue;
    }
  
    // ########### Step 1- Check of TPC and TOF response ####################

    Double_t ptrack=aodtrack->P();
    //#################### TPC PID #######################
     if (!(aodtrack->GetStatus()&AliESDtrack::kTPCpid )){
       // NO TPC PID INFO FOR THIS TRACK 
       hasPID[daught]--;
     }
     else {
       static AliTPCPIDResponse theTPCpid;
       AliAODPid *pidObj = aodtrack->GetDetPid();
       Double_t ptProng=pidObj->GetTPCmomentum();
       nsigmaTPCpi = theTPCpid.GetNumberOfSigmas(ptProng,(Float_t)pid->GetTPCsignal(),(Int_t)aodtrack->GetTPCClusterMap().CountBits(),AliPID::kPion);
       nsigmaTPCK =  theTPCpid.GetNumberOfSigmas(ptProng,(Float_t)pid->GetTPCsignal(),(Int_t)aodtrack->GetTPCClusterMap().CountBits(),AliPID::kKaon);
       //if(ptrack<0.6){
       if(ptProng<0.6){
	 if(TMath::Abs(nsigmaTPCK)<2.)isKaonPionTPC[daught][0]=1;
	 else if(TMath::Abs(nsigmaTPCK)>3.)isKaonPionTPC[daught][0]=-1;
	 if(TMath::Abs(nsigmaTPCpi)<2.)isKaonPionTPC[daught][1]=1;
	 else if(TMath::Abs(nsigmaTPCpi)>3.)isKaonPionTPC[daught][1]=-1;
       }
       //else if(ptrack<.8){
       else if(ptProng<.8){
	 if(TMath::Abs(nsigmaTPCK)<1.)isKaonPionTPC[daught][0]=1;
	 else if(TMath::Abs(nsigmaTPCK)>3.)isKaonPionTPC[daught][0]=-1;
	 if(TMath::Abs(nsigmaTPCpi)<1.)isKaonPionTPC[daught][1]=1;
	 else if(TMath::Abs(nsigmaTPCpi)>3.)isKaonPionTPC[daught][1]=-1;
       }     
       else {
	 //	if(nsigmaTPCK>-2.&&nsigmaTPCK<1.)isKaonPionTPC[daught][0]=1;
	 if(TMath::Abs(nsigmaTPCK)>3.)isKaonPionTPC[daught][0]=-1;
	 //if(nsigmaTPCpi>-1.&&nsigmaTPCpi<2.)isKaonPionTPC[daught][1]=1;
	 if(TMath::Abs(nsigmaTPCpi)>3.)isKaonPionTPC[daught][1]=-1;
       }
     }
    
    
    // ##### TOF PID: do not ask nothing for pion/protons ############
     if(!((aodtrack->GetStatus()&AliESDtrack::kTOFpid)&&(aodtrack->GetStatus()&AliESDtrack::kTOFout)&&(aodtrack->GetStatus()&AliESDtrack::kTIME))){
       // NO TOF PID INFO FOR THIS TRACK      
       hasPID[daught]--;
     }
     else{
       tofSig=pid->GetTOFsignal(); 
       pid->GetIntegratedTimes(times);
       if(TMath::Abs(tofSig-times[3])>3.*160.){
	 isKaonPionTOF[daught][0]=-1;
       }
       else {	 
	 if(ptrack<1.5){
	   isKaonPionTOF[daught][0]=1;
	 }
       }
     }
     
     //######### Step 2: COMBINE TOF and TPC PID ###############
     // we apply the following convention: if TPC and TOF disagree (discarded Vs identified) -> unknown
     combinedPID[daught][0]=isKaonPionTOF[daught][0]+isKaonPionTPC[daught][0];
     combinedPID[daught][1]=isKaonPionTOF[daught][1]+isKaonPionTPC[daught][1];
     
     
     //######### Step 3:   USE PID INFO     
     
     if(combinedPID[daught][0]<=-1&&combinedPID[daught][1]<=-1){// if not a K- and not a pi- both D0 and D0bar excluded
       isD0D0barPID[0]=0;
       isD0D0barPID[1]=0;
     }
     else if(combinedPID[daught][0]==2&&combinedPID[daught][1]>=1){// if in conflict (both pi- and K-), if k for both TPC and TOF -> is K
       if(aodtrack->Charge()==-1)isD0D0barPID[1]=0;//if K- D0bar excluded
       else isD0D0barPID[0]=0;// if K+ D0 excluded
     }
     else if(combinedPID[daught][0]==1&&combinedPID[daught][1]>=1){// if in conflict (both pi- and K-) and k- only for TPC or TOF -> reject
       isD0D0barPID[0]=0;
       isD0D0barPID[1]=0;
     }
     else if(combinedPID[daught][0]>=1||combinedPID[daught][1]<=-1){
       if(aodtrack->Charge()==-1)isD0D0barPID[1]=0;// not a D0bar if K- or if pi- excluded
       else isD0D0barPID[0]=0;//  not a D0 if K+ or if pi+ excluded
     }
     else if(combinedPID[daught][0]<=-1||combinedPID[daught][1]>=1){
       if(aodtrack->Charge()==-1)isD0D0barPID[0]=0;// not a D0 if pi- or if K- excluded
      else isD0D0barPID[1]=0;// not a D0bar if pi+ or if K+ excluded
     }
     
     // ##########  ALSO DIFFERENT TPC PID REQUEST FOR LOW pt D0: request of K identification      ###############################
     // ########## more tolerant criteria for single particle ID-> more selective criteria for D0   ##############################
     // ###############                     NOT OPTIMIZED YET                                  ###################################
     if(d->Pt()<2.){
       isKaonPionTPC[daught][0]=0;
       isKaonPionTPC[daught][1]=0;
       AliAODPid *pidObj = aodtrack->GetDetPid();
       Double_t ptProng=pidObj->GetTPCmomentum();
       //if(ptrack<0.6){
       if(ptProng<0.6){
	 if(TMath::Abs(nsigmaTPCK)<3.)isKaonPionTPC[daught][0]=1;
	 else if(TMath::Abs(nsigmaTPCK)>3.)isKaonPionTPC[daught][0]=-1;
	 if(TMath::Abs(nsigmaTPCpi)<3.)isKaonPionTPC[daught][1]=1;
	 else if(TMath::Abs(nsigmaTPCpi)>3.)isKaonPionTPC[daught][1]=-1;
     }
       //else if(ptrack<.8){
       else if(ptProng<.8){
	 if(TMath::Abs(nsigmaTPCK)<2.)isKaonPionTPC[daught][0]=1;
	 else if(TMath::Abs(nsigmaTPCK)>3.)isKaonPionTPC[daught][0]=-1;
	 if(TMath::Abs(nsigmaTPCpi)<3.)isKaonPionTPC[daught][1]=1;
	 else if(TMath::Abs(nsigmaTPCpi)>3.)isKaonPionTPC[daught][1]=-1;
       }     
       else {
	 if(TMath::Abs(nsigmaTPCK)>3.)isKaonPionTPC[daught][0]=-1;
	 if(TMath::Abs(nsigmaTPCpi)>3.)isKaonPionTPC[daught][1]=-1;
       }
     }
     
  }// END OF LOOP ON DAUGHTERS
  
  // FURTHER PID REQUEST (both daughter info is needed)
  if(combinedPID[0][0]<=-1&&combinedPID[1][0]<=-1){
    fWhyRejection=31;// reject cases in which no kaon-compatible tracks are found
    return 0;
  }
  else if(hasPID[0]==0&&hasPID[1]==0){
    fWhyRejection=28;// reject cases in which no PID info is available  
    return 0;
  }
  if(d->Pt()<2.){
    // request of K identification at low D0 pt
    combinedPID[0][0]=0;
    combinedPID[0][1]=0;
    combinedPID[1][0]=0;
    combinedPID[1][1]=0;
    
    combinedPID[0][0]=isKaonPionTOF[0][0]+isKaonPionTPC[0][0];
    combinedPID[0][1]=isKaonPionTOF[0][1]+isKaonPionTPC[0][1];
    combinedPID[1][0]=isKaonPionTOF[1][0]+isKaonPionTPC[1][0];
    combinedPID[1][1]=isKaonPionTOF[1][1]+isKaonPionTPC[1][1];
    
    if(combinedPID[0][0]<=0&&combinedPID[1][0]<=0){
      fWhyRejection=32;// reject cases where the Kaon is not identified
      return 0;
    }
  }

  //  cout<<"Why? "<<fWhyRejection<<endl;  
  return isD0D0barPID[0]+isD0D0barPID[1];
}



//---------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpi::CombineSelectionLevels(Int_t selectionvalTrack,
						 Int_t selectionvalCand,
						 Int_t selectionvalPID) const
{
  //
  // This method combines the tracks, PID and cuts selection results
  //
  if(selectionvalTrack==0) return 0;

  Int_t returnvalue;

  switch(selectionvalPID) {
  case 0:
    returnvalue=0;
    break;
  case 1:
    returnvalue=((selectionvalCand==1 || selectionvalCand==3) ? 1 : 0);
    break;
  case 2:
    returnvalue=((selectionvalCand==2 || selectionvalCand==3) ? 2 : 0);
    break;
  case 3:
    returnvalue=selectionvalCand;
    break;
  default:
    returnvalue=0;
    break;
  }

  return returnvalue;
}
//----------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpi::IsSelectedSpecialCuts(AliAODRecoDecayHF *d) const
{
  //
  // Note: this method is temporary
  // Additional cuts on decay lenght and lower cut for d0 norm are applied using vertex without candidate's daughters
  //

  //apply cuts

  Float_t normDecLengthCut=1.,decLengthCut=TMath::Min(d->P()*0.0066+0.01,0.06/*cm*/), normd0Cut=0.5;
  // "decay length" expo law with tau' = beta*gamma*ctau= p/m*ctau =p*0.0123/1.864~p*0.0066
  // decay lenght > ctau' implies to retain (1-1/e) (for signal without considering detector resolution), 

  Int_t returnvalue=3; //cut passed
  for(Int_t i=0;i<2/*prongs*/;i++){
    if(TMath::Abs(d->Normalizedd0Prong(i))<normd0Cut) return 0; //normd0Cut not passed
  }
  if(d->DecayLength()<decLengthCut)  return 0; //decLengthCut not passed
  if(d->NormalizedDecayLength()<normDecLengthCut)  return 0; //decLengthCut not passed
    

  return returnvalue;
}

