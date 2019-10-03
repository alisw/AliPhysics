#include <THnSparse.h>
#include <TTreeStream.h>
#include <AliLog.h>

#include <AliESDEvent.h>
#include <AliESDtrack.h>

#include "AliESDUtils.h"

#include "AliMESpidTask.h"
#include "AliMESeventInfo.h"
#include "AliMEStrackInfo.h"

ClassImp(AliMESpidTask)

//________________________________________________________________________
AliMESpidTask::AliMESpidTask()
  : AliMESbaseTask()
{
  //
  // Constructor
  //
}

//________________________________________________________________________
AliMESpidTask::AliMESpidTask(const char *name)
  : AliMESbaseTask(name)
{
  //
  // Constructor
  //
}


//________________________________________________________________________
AliMESpidTask::~AliMESpidTask()
{
  //
  // Destructor
  //
}

//________________________________________________________________________
void AliMESpidTask::UserCreateOutputObjects()
{
  //define user data containers
  AliMESbaseTask::UserCreateOutputObjects();

  //define extra user containers
}

//________________________________________________________________________
void AliMESpidTask::UserExec(Option_t *opt)
{
  // Run user analysis. The following objects are allocated after calling AliMESbaseTask::UserExec(opt)
  // fEvInfo  -  reconstructed event information (class AliMESeventInfo)
  // fTracks  -  reconstructed array of tracks (class TObjArray of AliMEStrackInfo)
  // fMCevInfo-  MC event information (class AliMESeventInfo)
  // fMCtracks-  MC array of tracks (class TObjArray of AliMEStrackInfo)

  AliMESbaseTask::UserExec(opt);


  if( !fEvInfo->HasTriggerMB() ) return;
  // if( !fEvInfo->HasTriggerHM() ) return;

 /*

 // !!!!!!!!!!
 // These are meaningless as long as AliPPVsMultUtils:IsSelected() is used in AliMEStender
 // !!!!!!!!!!

   if( !fEvInfo->HasVertex() ) return;
   vec_hNoEvts[0] = 0.;
   hNoEvts->Fill(vec_hNoEvts);

   if( fEvInfo->IsPileUp() ) return;
   vec_hNoEvts[0] = 1.;
   hNoEvts->Fill(vec_hNoEvts);

   if (TMath::Abs(fEvInfo->GetVertexZ()) > 10.) return;
   vec_hNoEvts[0] = 2.;
   hNoEvts->Fill(vec_hNoEvts);

   if(mult_comb08 < 0. || mult_V0M < 0.) return;
 //   if(mult_comb08 < 0.){
     return;
   }

   // !!!!!!!!!!
 */

  // number of events counter
  Double_t vec_hNoEvts[7];   // vector used to fill hNoEvts
  THnSparseD *hNoEvts = (THnSparseD*)fHistosQA->At(slot_NoEvts);

  Double_t mult_comb08 = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb);  	// combined multiplicity with |eta| < 0.8
  Double_t mult_V0M = fEvInfo->GetMultiplicity(AliMESeventInfo::kV0M);   		// V0M percentile

  // Double_t mult_comb0408 = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb0408);		// combined multiplicity with 0.4 < |eta| < 0.8

  // event shape for data (from ESD)
  Double_t directivity_plus = fEvInfo->GetEventShape()->GetDirectivity(1);
  Double_t directivity_minus = fEvInfo->GetEventShape()->GetDirectivity(0);

  vec_hNoEvts[0] = 0.;
  // hNoEvts->Fill(vec_hNoEvts);


  Double_t directivity = -2;

  // select events with both dirs in the same interval
  const Int_t lenght = 4;
  Double_t intervals[lenght] = {0., 0.3, 0.6, 1.0};

  // NOTE: the intervals are considered half-closed: (a,b]
  if( (directivity_plus >= intervals[0]) && (directivity_plus <= intervals[lenght-1]) ){

      Int_t first = -1;
      for(Int_t i=1; i<lenght; i++){
          if(directivity_plus <= intervals[i]){
    	    first = i;
    		break;
          }
      }

      if( (directivity_minus > intervals[first-1]) && (directivity_minus <= intervals[first]) ){
          directivity =  (directivity_plus + directivity_minus) / 2.0;
      }
  }
/*
  // select events with both dirs close
  if(TMath::Abs(directivity_plus - directivity_minus) < 0.2){
    directivity = (directivity_plus + directivity_minus) / 2.0;
  }

  // select events using only dir plus
  Double_t directivity = directivity_plus;
*/

  // event shape for MC (from MC event)
  Double_t MC_directivity_plus = 0;
  Double_t MC_directivity_minus = 0;
  Double_t MC_directivity = 0;
  if( HasMCdata() ){ // run only on MC
      MC_directivity_plus = fMCevInfo->GetEventShape()->GetDirectivity(1);
      MC_directivity_minus = fMCevInfo->GetEventShape()->GetDirectivity(0);
      MC_directivity =  (MC_directivity_plus + MC_directivity_minus) / 2.0;
      MC_directivity = MC_directivity_plus;
  }

/*
//   ((TH2*)fHistosQA->At(3))->Fill(0., ev_mult);
  // ESD event level cuts
//   AliInfo(Form("\n\n\nHas trigger MB: %d",fEvInfo->HasTriggerMB()));
//   AliInfo(Form("Has trigger HM: %d",fEvInfo->HasTriggerHM()));
  if( !fEvInfo->HasTriggerMB() ) return;
//   ((TH2*)fHistosQA->At(3))->Fill(1, ev_mult);
 //   AliInfo(Form("Has vertex: %d",fEvInfo->HasVertex()));
//   AliInfo (Form("vertex Z = %g",fEvInfo->GetVertexZ()));
  if( !fEvInfo->HasVertex() ) return;
  ((TH2*)fHistosQA->At(3))->Fill(0., ESDmult);
  //   AliInfo(Form("IsPileUp: %d",fEvInfo->IsPileUp()));
  if( fEvInfo->IsPileUp() ) return;
  ((TH2*)fHistosQA->At(3))->Fill(1, ESDmult);
  if (TMath::Abs(fEvInfo->GetVertexZ()) > 10.) return;
  ((TH2*)fHistosQA->At(3))->Fill(2, ESDmult);
  //   if(fEvInfo->GetMultiplicity(AliMESeventInfo::kComb) < 0.) return;
  if(ESDmult < 0.) return;
  ((TH2*)fHistosQA->At(3))->Fill(3, ESDmult);
*/

//   vec_hNoEvts[1] = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb);				// combined multiplicity with |eta| < 0.8
  vec_hNoEvts[1] = mult_comb08;				// combined multiplicity with |eta| < 0.8
//   vec_hNoEvts[2] = fEvInfo->GetMultiplicity(AliMESeventInfo::kV0M);				// V0M percentile
  vec_hNoEvts[2] = mult_V0M;				// V0M percentile
//   vec_hNoEvts[3] = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb0408);		// combined multiplicity with 0.4 < |eta| < 0.8
// vec_hNoEvts[3] = mult_comb0408;		// combined multiplicity with 0.4 < |eta| < 0.8
  vec_hNoEvts[3] = directivity;		// combined multiplicity with 0.4 < |eta| < 0.8

  vec_hNoEvts[4] = 0;
  vec_hNoEvts[5] = 0;
  vec_hNoEvts[6] = 0;
  if( HasMCdata() ){
    vec_hNoEvts[4] = fMCevInfo->GetMultiplicity(AliMESeventInfo::kGlob08);
	  vec_hNoEvts[5] = fMCevInfo->GetMultiplicity(AliMESeventInfo::kV0M);
    vec_hNoEvts[6] = MC_directivity;
  }


  vec_hNoEvts[0] = 3.;
  hNoEvts->Fill(vec_hNoEvts);


  // used ONLY for systematic studies (see line 377)
  AliESDEvent* fESD = NULL;
  // if(DebugLevel()>0){
	  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
	  if (!fESD) {
		  AliError("ESD event not available");
		  return;
	  }
  // }

// 	AliInfo("\n\nStarting the tracks loop:");
// 	printf("tracks = %i\n",fTracks->GetEntries());

  // AliInfo(Form("\n\n VOA signal maybe: %f\n\n", AliESDUtils::GetCorrV0A(fESD->GetVZEROData()->GetMTotV0A(), fESD->GetPrimaryVertexSPD()->GetZ())));
  Double_t V0Asignal = AliESDUtils::GetCorrV0A(fESD->GetVZEROData()->GetMTotV0A(), fESD->GetPrimaryVertexSPD()->GetZ());
  Double_t V0Csignal = AliESDUtils::GetCorrV0C(fESD->GetVZEROData()->GetMTotV0C(), fESD->GetPrimaryVertexSPD()->GetZ());
  // Double_t V0Asignal = 0;
  // AliInfo(Form("V0A signal = %f\t V0C signal = %f\n", V0Asignal, V0Csignal));

  THnSparseD *hMultEst = (THnSparseD*)fHistosQA->At(slot_MultEst);
  Double_t vec_hMultEst[12]; // vector used to fill hMultEst
//   vec_hMultEst[0] = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb);
//   vec_hMultEst[1] = fEvInfo->GetMultiplicity(AliMESeventInfo::kV0M);
//   vec_hMultEst[2] = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb0408);
  vec_hMultEst[0] = mult_comb08;
  vec_hMultEst[1] = mult_V0M;
  // vec_hMultEst[2] = mult_comb0408;
  vec_hMultEst[2] = V0Asignal + V0Csignal;
  vec_hMultEst[3] = directivity;
  if( HasMCdata() ){
    vec_hMultEst[6] = fMCevInfo->GetMultiplicity(AliMESeventInfo::kGlob08);
 		vec_hMultEst[7] = fMCevInfo->GetMultiplicity(AliMESeventInfo::kV0M);
  	vec_hMultEst[8] = fMCevInfo->GetMultiplicity(AliMESeventInfo::kComb0408);
    vec_hMultEst[9] = MC_directivity;
  }

  // get the leading particle direction
  Double_t px_LP=0., py_LP=0., phi_LP=0.;
  px_LP = fEvInfo->GetEventShape()->GetMomLeading(kTRUE);
  py_LP = fEvInfo->GetEventShape()->GetMomLeading(kFALSE);
  phi_LP = TMath::ATan2(py_LP, px_LP);
  phi_LP = (phi_LP>0) ? phi_LP : (phi_LP+TMath::TwoPi());  // if negative add 2*pi

  Double_t px_LP_MC=0., py_LP_MC=0., phi_LP_MC=0.;
  if( HasMCdata() ){
      // get the leading particle direction
      px_LP_MC = fMCevInfo->GetEventShape()->GetMomLeading(kTRUE);
      py_LP_MC = fMCevInfo->GetEventShape()->GetMomLeading(kFALSE);
      phi_LP_MC = TMath::ATan2(py_LP_MC, px_LP_MC);
      phi_LP_MC = (phi_LP_MC>0) ? phi_LP_MC : (phi_LP_MC+TMath::TwoPi());  // if negative add 2*pi
  }


  Double_t vec_hPIDQA[8];			//  vector used to fill hPIDQA
/*
  TTree *testTree = (TTree*)fHistosQA->At(slot_testTree);
  testTree->Branch("p", &vec_hPIDQA[0], "p/D");
  testTree->Branch("charge", &vec_hPIDQA[1], "charge/D");
  testTree->Branch("TPC", &vec_hPIDQA[5], "TPC/D");
  testTree->Branch("TOF", &vec_hPIDQA[6], "TPC/D");
  testTree->Branch("MCpid", &vec_hPIDQA[7], "MCpid/D");
*/

  Double_t y_LP_ESD = -9999.;
  Double_t y_LP_MC = -9999.;
  Double_t pT_LP_ESD = -9999.;
  Double_t pT_LP_MC = -9999.;
  
  // ESD track loop
  AliMEStrackInfo *t(NULL), *tMC(NULL);
  for(Int_t it(0); it<fTracks->GetEntries(); it++){
    if(!(t = (AliMEStrackInfo*)fTracks->At(it))) continue;
// 	if(!HasMCdata() || !(tMC= (AliMEStrackInfo*)fMCtracks->At(t->GetLabel()))) continue;
	if( HasMCdata() ){
		if( !(tMC= (AliMEStrackInfo*)fMCtracks->At(t->GetLabel())) ) continue;

		// 	AliInfo(Form("%2d p_t[GeV/c] = %f  p_t[GeV/c] = %f", it, t->Pt(), tMC->Pt()));
		// 	AliInfo(Form("%2d p[GeV/c] = %f  p[GeV/c] = %f", it, t->P(), tMC->P()));
		// 	AliInfo(Form("ESD charge: %i", t->Charge()));
		// 	AliInfo(Form("MC charge: %i\n", tMC->Charge()));
	}


	Double_t vec_hAllESD[14];    	// vector used to fill hAllESD
  Double_t vec_hDeltaPhi[9];		// vector used to fill hDeltaPhi


	THnSparseD *hAllESD = (THnSparseD*)fHistosQA->At(slot_AllESD);
  // enum axis_hAllESD {l_comb08, l_V0M, l_comb0408, l_pT, l_charge, l_pidTPC, l_pidTOF, l_rapidity, l_TOFmatching, l_MCPID, l_yMCPID, l_MCprimary};  // labels for the hAllESD axis
	enum axis_hAllESD {l_comb08, l_V0M, l_directivity, l_pT, l_charge, l_pidTPC, l_pidTOF, l_rapidity, l_TOFmatching, l_delta_phi, l_delta_y, l_MCPID, l_yMCPID, l_MCprimary};  // labels for the hAllESD axis

	// ---------------------------
	// get ESD multiplicity
// 	vec_hAllESD[0] = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb);
// 	vec_hAllESD[0] = fEvInfo->GetMultiplicity(AliMESeventInfo::kV0M);
// 	vec_hAllESD[0] = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb0408);
// 	vec_hAllESD[0] = mult_comb08;
	vec_hAllESD[l_comb08] = mult_comb08;
// 	AliInfo(Form("mult ESD = %g",vec_hAllESD[0]));
	vec_hAllESD[l_V0M] = mult_V0M;
    // vec_hAllESD[l_comb0408] = mult_comb0408;
	vec_hAllESD[l_directivity] = directivity;

	// ---------------------------
	// get pT
// 	vec_hAllESD[1] = t->Pt();
	vec_hAllESD[l_pT] = t->Pt();

	// ---------------------------
	// get charge
// 	vec_hAllESD[2] = t->Charge();
	vec_hAllESD[l_charge] = t->Charge();
// 	vec_hPIDQA[1] = vec_hAllESD[2];
	vec_hPIDQA[1] = vec_hAllESD[l_charge];

	// ---------------------------
	// make PID
	// TPC
	Int_t maxIndex = -9999;
	const Double_t *prob_TPC = t->GetPID()->GetProb(AliMEStrackInfo::kITS);
	Double_t maxProb = 0; // max probability
	for(Int_t i = 0; i<AliPID::kSPECIES; i++){
// 		AliInfo(Form("probTPC[%i]= %g", i, prob_TPC[i]));
		if(prob_TPC[i] > maxProb){
			maxProb = prob_TPC[i];
			maxIndex = i;
		}
	}
// 	vec_hAllESD[3] = maxIndex;    // particle identification
	vec_hAllESD[l_pidTPC] = maxIndex;    // particle identification
// 	AliInfo(Form("maxIndex TPC = %i\n\n", maxIndex));
// vec_hPIDQA[2] = vec_hAllESD[3];
	vec_hPIDQA[2] = vec_hAllESD[l_pidTPC];

	// TOF
	maxIndex = -9999;
	const Double_t *prob_TOF = t->GetPID()->GetProb(AliMEStrackInfo::kTOF);
	maxProb = 0; // max probability
	for(Int_t i = 0; i<AliPID::kSPECIES; i++){
//  		AliInfo(Form("probTOF[%i]= %f", i, prob_TOF[i]));
		if(prob_TOF[i] > maxProb){
			maxProb = prob_TOF[i];
			maxIndex = i;
		}
	}
// 	vec_hAllESD[4] = maxIndex;    // particle identification with TPC - TOF
	vec_hAllESD[l_pidTOF] = maxIndex;    // particle identification with TPC - TOF
// 		AliInfo(Form("maxIndex TOF = %i\n\n", maxIndex));
// 	vec_hPIDQA[3] = vec_hAllESD[4];
	vec_hPIDQA[3] = vec_hAllESD[l_pidTOF];


	// ---------------------------
	// compute y after PID
	Double_t  mass[AliPID::kSPECIES] = {0.00051, 0.10565, 0.13957, 0.49368, 0.93827};
	Double_t e = TMath::Sqrt(t->P()*t->P() + mass[maxIndex]*mass[maxIndex]);
// 	Double_t pz = t->P() * TMath::Cos(t->Theta());
// 	if( TMath::Abs(t->Pz()) != e ) vec_hAllESD[5] = 0.5*TMath::Log((e + t->Pz())/(e - t->Pz()));
// 	else vec_hAllESD[5] = -9999;
// 	if(TMath::Abs(vec_hAllESD[5]) > 1.0) continue;
	if( TMath::Abs(t->Pz()) != e ) vec_hAllESD[l_rapidity] = 0.5*TMath::Log((e + t->Pz())/(e - t->Pz()));
	else vec_hAllESD[l_rapidity] = -9999;
	if(TMath::Abs(vec_hAllESD[l_rapidity]) > 1.0) continue;
// 	AliInfo(Form("Pz = %g \te = %g \t e+Pz = %g \t e-Pz = %g\n", t->Pz(), e, (e + t->Pz()), (e - t->Pz())));

  // set the y for the LP (NOTE: the LP is always the first one in the list)
  if(it == 0){
    pT_LP_ESD = vec_hAllESD[l_pT];
    y_LP_ESD = vec_hAllESD[l_rapidity];
  }
	// ---------------------------
	// get the TOF mistmatch
// 	vec_hAllESD[6] = (((t->GetPID()->GetTOFmisProb())<0.01)?1:0);
// 	vec_hAllESD[6] = (t->GetPID()->GetTOFmisProb());
	vec_hAllESD[l_TOFmatching] = (t->GetPID()->GetTOFmisProb());
// 	printf("misProb = %g\n",t->GetPID()->GetTOFmisProb());
// 	vec_hPIDQA[4] = vec_hAllESD[6];
	vec_hPIDQA[4] = vec_hAllESD[l_TOFmatching];

  // ---------------------------
  // get the delta phi angle
  vec_hAllESD[l_delta_phi] = ComputeDeltaPhi(t->Phi(), phi_LP);

  // ---------------------------
  // compute delta y
  vec_hAllESD[l_delta_y] = y_LP_ESD - vec_hAllESD[l_rapidity];


  // fill the deltaPhi sparse
  THnSparseD *hDeltaPhi = (THnSparseD*)fHistosQA->At(slot_DeltaPhi);
  // vec_hDeltaPhi[0] = t->Pt();
  vec_hDeltaPhi[0] = mult_comb08;
  vec_hDeltaPhi[1] = directivity;
  vec_hDeltaPhi[2] = vec_hAllESD[l_pT]; // pT ESD
  vec_hDeltaPhi[3] = vec_hAllESD[l_delta_y];  // delta_y ESD
  vec_hDeltaPhi[4] = vec_hAllESD[l_delta_phi];  // delta_phi ESD
  // generated:
  vec_hDeltaPhi[5] = 0; // pt MC
  vec_hDeltaPhi[6] = 0; // delta_y MC
  vec_hDeltaPhi[7] = 0; // delta_phi MC
  vec_hDeltaPhi[8] = 0; // deltaPhi ESD vs MC LP


	if( HasMCdata() ){ // run only on MC
		// ---------------------------
		// get the MC PDG code
		Int_t MC_identity  = -9999;

		if(tMC){
			const Double_t *MCpdg = tMC->GetPID()->GetProb(AliMEStrackInfo::kITS);
			for(Int_t i = 0; i<AliPID::kSPECIES; i++){
				// 			AliInfo(Form("MC: probITS[%i]= %g", i, MCpdg[i]));
				if( MCpdg[i] > 0 ){
// 					vec_hAllESD[7] = i;
					vec_hAllESD[l_MCPID] = i;
					MC_identity = i;
				}
			}
			// 		AliInfo(Form("MC PID = %g", vec_hAllESD[4]));
// 			vec_hPIDQA[7] = vec_hAllESD[7];
			vec_hPIDQA[7] = vec_hAllESD[l_MCPID];

// 			if( (tMC->HasOrigin(AliMEStrackInfo::kPrimary)) ) vec_hAllESD[9] = 1.;
// 			else vec_hAllESD[9] = 0.;
			if( (tMC->HasOrigin(AliMEStrackInfo::kPrimary)) ) vec_hAllESD[l_MCprimary] = 1.;
			else vec_hAllESD[l_MCprimary] = 0.;
		}

		// ---------------------------
		// compute y using true PID
		Double_t eMC = TMath::Sqrt(t->P()*t->P() + mass[MC_identity]*mass[MC_identity]);
		// 	Double_t pz = t->P() * TMath::Cos(t->Theta());
// 		if( TMath::Abs(t->Pz()) != eMC ) vec_hAllESD[8] = 0.5*TMath::Log((eMC + t->Pz())/(eMC - t->Pz()));
// 		else vec_hAllESD[8] = -9999;
// 		if(TMath::Abs(vec_hAllESD[8]) > 1.0) continue;
		if( TMath::Abs(t->Pz()) != eMC ) vec_hAllESD[l_yMCPID] = 0.5*TMath::Log((eMC + t->Pz())/(eMC - t->Pz()));
		else vec_hAllESD[l_yMCPID] = -9999;
		if(TMath::Abs(vec_hAllESD[l_yMCPID]) > 1.0) continue;

    vec_hDeltaPhi[5] = tMC->Pt();
    
    // delta_y MC
    y_LP_MC = ((AliMEStrackInfo*)fMCtracks->At(0))->Y();
    pT_LP_MC = ((AliMEStrackInfo*)fMCtracks->At(0))->Pt();
    vec_hDeltaPhi[6] = y_LP_MC - tMC->Y();
    
    vec_hDeltaPhi[7] = ComputeDeltaPhi(tMC->Phi(), phi_LP_MC);  // gen info
    vec_hDeltaPhi[8] = ComputeDeltaPhi(t->Phi(), phi_LP_MC);    // rec tracks vs gen LP
	}

	// ---------------------------
	// fill the hSparse
    hDeltaPhi->Fill(vec_hDeltaPhi);
	hAllESD->Fill(vec_hAllESD);


	// fill the PID QA sparse
	THnSparseD *hPIDQA = (THnSparseD*)fHistosQA->At(slot_PIDQA);
	vec_hPIDQA[0] = t->P();
	vec_hPIDQA[5] = t->GetdEdx();
	vec_hPIDQA[6] = t->Getbeta();
	hPIDQA->Fill(vec_hPIDQA);

    // fill also the testTree
    // testTree->Fill();


/*
	// THIS IS USED ONLY WHEN RUNNING FOR DCA SHAPES (change AliMEStender.cxx line 163 kTRUE->kFALSE)

	Double_t valDCA[12];

	THnSparseD *hDCA = (THnSparseD*)fHistosQA->At(slot_DCA);

	valDCA[0] = vec_hAllESD[l_comb08];		// multiplicity
    valDCA[1] = vec_hAllESD[l_directivity];	// directivity
	valDCA[2] = vec_hAllESD[l_pT];		    // pT
	valDCA[3] = vec_hAllESD[l_charge];		// charge
	valDCA[4] = vec_hAllESD[l_pidTPC];		// PID_TPC
	valDCA[5] = vec_hAllESD[l_pidTOF];		// PID_TPCTOF
	valDCA[6] = vec_hAllESD[l_rapidity];	// y
    valDCA[7] = vec_hAllESD[l_TOFmatching];	// TOFmatching
	valDCA[8] = vec_hAllESD[l_delta_phi];	// deltaPhi

	// fill DCA
	Double_t dca[2];  // 0 = xy; 1 = z
	t->GetDCA(dca);
	valDCA[10] = dca[0];    // DCA
// 	AliInfo(Form("DCAxy() = %g", valDCA[8]));

	// check if it passes the DCA cut
	if( TMath::Abs(valDCA[10]) < (0.0182+0.0350/ TMath::Power(valDCA[2],1.01)) ){  // valDCA[2] = pT; valDCA[10] = DCAxy
		valDCA[11] = 1.;   // pass_DCAcut
	}
	else{
		valDCA[11] = 0.;   // pass_DCAcut
	}


	if( HasMCdata() ){ // run only on MC
		if(tMC){
			// get origin
			if( (tMC->HasOrigin(AliMEStrackInfo::kPrimary)) ){
				valDCA[9] = 0.;     // MCorigin
			}
			else{
				if( (tMC->HasOrigin(AliMEStrackInfo::kSecondary)) ){
					valDCA[9] = 1.;     // MCorigin
				}
				else{
					if( (tMC->HasOrigin(AliMEStrackInfo::kMaterial)) ){
						valDCA[9] = 2.;     // MCorigin
					}
					else{
						valDCA[9] = -1.;     // MCorigin
					}
				}
			}
		}
	}

	hDCA->Fill(valDCA);
*/

// 	if( HasMCdata() ){ // run only on MC
	if(DebugLevel()>0){ // used ONLY for systematic studies
    Double_t vec_ESDtree[2];
    vec_ESDtree[0] = t->Eta();
		vec_ESDtree[1] = t->Phi();
// 		vec_ESDtree[2] = t->GetPID()->GetTOFmisProb();
		// 		vec_ESDtree[3] = tMC->P();
		// 		vec_ESDtree[4] = t->Charge();
		// 		vec_ESDtree[5] = tMC->Charge();

		// get the ESD track
		if(!fESD) continue;
		AliESDtrack* track = fESD->GetTrack(it);
		if (!track) {
			Printf("ERROR: Could not receive track %d", it);
			continue;
		}

		// very ugly hack allowing for the offline use of esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
		// get the Chi2TPCConstrainedGlobal here and put it in the tree (copied from AliESDtrackCuts::AcceptTrack)
		enum VertexType { kVertexTracks = 0x1, kVertexSPD = 0x2, kVertexTPC = 0x4 };
		Int_t fCutMaxChi2TPCConstrainedVsGlobalVertexType;
		fCutMaxChi2TPCConstrainedVsGlobalVertexType = kVertexTracks | kVertexSPD;

		const AliESDVertex* vertex = 0;
		if (fCutMaxChi2TPCConstrainedVsGlobalVertexType & kVertexTracks)
			vertex = fESD->GetPrimaryVertexTracks();

		if ((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & kVertexSPD)
			vertex = fESD->GetPrimaryVertexSPD();

		if ((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & kVertexTPC)
			vertex = fESD->GetPrimaryVertexTPC();

		Double_t chi2TPCConstrainedVsGlobal = -1.0;
		if (vertex->GetStatus())
			chi2TPCConstrainedVsGlobal = track->GetChi2TPCConstrainedVsGlobal(vertex);

        // fill debug
        (*AliMESbaseTask::DebugStream()) << "pidTrk"
		<<"mult08=" << vec_hAllESD[l_comb08]
		<<"multV0M=" << vec_hAllESD[l_V0M]
        // <<"mult0408=" << vec_hAllESD[l_comb0408]
		<<"directivity=" << vec_hAllESD[l_directivity]
        <<"pT=" << vec_hAllESD[l_pT]
        <<"charge=" << vec_hAllESD[l_charge]
        <<"TPC_PID=" << vec_hAllESD[l_pidTPC]
        <<"TOF_PID=" << vec_hAllESD[l_pidTOF]
        <<"y=" << vec_hAllESD[l_rapidity]
        <<"mismatch=" << vec_hAllESD[l_TOFmatching]
        <<"MC_PID=" << vec_hAllESD[l_MCPID]
        <<"MC_y=" << vec_hAllESD[l_yMCPID]
        <<"MC_primary=" << vec_hAllESD[l_MCprimary]
        // ----
        <<"eta=" << vec_ESDtree[0]
        <<"phi=" << vec_ESDtree[1]
        // ----
        // 			<< "dEta=" << dEta
        // 			<< "dPhi=" << dPhi
        <<"track.=" << track
//         <<"event.=" << fESD
		<< "chi2TPCConstrainedVsGlobal=" << chi2TPCConstrainedVsGlobal
        << "\n";
    }
// 	}

  } // end ESD track loop



  if( HasMCdata() ){ // run only on MC

// 	AliInfo("\n\nNew event");

    enum axis_hGen {l_MC_comb08, l_MC_directivity, l_MC_pT, l_MC_charge, l_MC_PID, l_MC_rapidity, l_MC_delta_phi, l_MC_ESDmult, l_MC_ESDdir};  // labels for the hAllESD axis


  	for(Int_t it(0); it<fMCtracks->GetEntries(); it++){
    	if(!(tMC = (AliMEStrackInfo*)fMCtracks->At(it))) continue;
    	//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	// review this for DCA fits shapes
		  if( !(tMC->HasOrigin(AliMEStrackInfo::kPrimary)) ) continue;
		  // add this as a extra dimension on the sparse
      // if(TMath::Abs(tMC->Eta()) > 0.9) continue;
		  if(TMath::Abs(tMC->Y()) > 1.) continue;
    	//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    	Double_t vec_hGen[9];  // vector used to fill hGen

		  THnSparseD *hGen = (THnSparseD*)fHistosQA->At(slot_Gen);

  		// ---------------------------
  		// get generated multiplicity
  		vec_hGen[l_MC_comb08] = fMCevInfo->GetMultiplicity(AliMESeventInfo::kGlob08);
      // 		vec_hGen[0] = fMCevInfo->GetMultiplicity(AliMESeventInfo::kComb0408);

      // ---------------------------
      // get generated directivity
      vec_hGen[l_MC_directivity] = MC_directivity;

  		// ---------------------------
  		// get pT
  		vec_hGen[l_MC_pT] = tMC->Pt();
      //  AliInfo(Form("pT = %g", vec_hGen[0]));

  		// ---------------------------
  		// get charge
  		vec_hGen[l_MC_charge] = tMC->Charge();

  		// ---------------------------
  		// get the MC PDG code
  		const Double_t *MCpdg = tMC->GetPID()->GetProb(AliMEStrackInfo::kITS);
  		for(Int_t i = 0; i<AliPID::kSPECIES; i++){
      // 			AliInfo(Form("MC: probITS[%i]= %g", i, MCpdg[i]));
  			if( MCpdg[i] > 0 )	vec_hGen[l_MC_PID] = i;
  		}

  		// ---------------------------
  		// get y
  		vec_hGen[l_MC_rapidity] = tMC->Y();

      // ---------------------------
      // get the delta phi angle
      vec_hGen[l_MC_delta_phi] = ComputeDeltaPhi(tMC->Phi(), phi_LP_MC);

  		// ---------------------------
  		// get the ESD multiplicity
  		// 		vec_hGen[5] = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb);
      vec_hGen[l_MC_ESDmult] = mult_comb08;

      // ---------------------------
  		// get the ESD directivity
  		vec_hGen[l_MC_ESDdir] = directivity;

  		// ---------------------------
  		// fill the hSparse
  		hGen->Fill(vec_hGen);
      
      if(it == 0){
        pT_LP_MC = vec_hGen[l_MC_pT];
        y_LP_MC = vec_hGen[l_MC_rapidity];
      }

/*
		TH1D *testCounter = (TH1D*)fHistosQA->At(3);
		if(vec_hGen[2] > 0){
			testCounter->Fill(1.);
			if(TMath::Abs(vec_hGen[4]) < 0.5){PWGLF/SPECTRA/MultEvShape/AliMESpidTask.cxx
				testCounter->Fill(2.);
				if(vec_hGen[3] == 2){
					AliInfo(Form("pT = %g \t charge = %g \t yMC = %g \t MC PID = %g", vec_hGen[1], vec_hGen[2], vec_hGen[4], vec_hGen[3]));
					testCounter->Fill(3.);
				}
			}
		}
*/

  	} // end MC track loop

  } // end HasMCdata IF

  vec_hMultEst[4] = pT_LP_ESD;
  vec_hMultEst[5] = y_LP_ESD;
  vec_hMultEst[10] = pT_LP_MC;
  vec_hMultEst[11] = y_LP_MC;
  hMultEst->Fill(vec_hMultEst);
  

}  // end UserExec

//________________________________________________________________________
Bool_t AliMESpidTask::PostProcess()
{
  return kTRUE;
}

//________________________________________________________
Bool_t AliMESpidTask::BuildQAHistos()
{
  // Make QA sparse histos for
  fHistosQA = new TList(); fHistosQA->SetOwner(kTRUE);

  // multiplicity estimators correlations && generated multiplicity
  Double_t binLimitsV0M[] = {0.0,0.01,0.1,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.0,31.0,32.0,33.0,34.0,35.0,36.0,37.0,38.0,39.0,40.0,41.0,42.0,43.0,44.0,45.0,46.0,47.0,48.0,49.0,50.0,51.0,52.0,53.0,54.0,55.0,56.0,57.0,58.0,59.0,60.0,61.0,62.0,63.0,64.0,65.0,66.0,67.0,68.0,69.0,70.0,71.0,72.0,73.0,74.0,75.0,76.0,77.0,78.0,79.0,80.0,81.0,82.0,83.0,84.0,85.0,86.0,87.0,88.0,89.0,90.0,91.0,92.0,93.0,94.0,95.0,96.0,97.0,98.0,99.0,100.0};
  
  // use for matching, PID and contaminations efficiency
  Double_t binLimits[] = {0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0};
  
  const Int_t ndim(12);
  const Int_t cldNbins[ndim]   = {105, 102, 50, 21, 52, 20, 102, 102, 20, 21, 52, 20};
  const Double_t cldMin[ndim]  = {-5, 0., 0., 0., 0., -1., -1.5, -1.5, -2.5, 0., 0., -1.},
  cldMax[ndim]  = {100., 100., 500., 1.05, 5., 1., 100.5, 100.5, 1.5, 1.05, 5, 1.};
  // THnSparseD *hMultEst = new THnSparseD("hMultEst","hMultEst;combined 0.8;V0M;combined 0.4-0.8; directivity; generated 0.8;generated V0M;generated 0.4-0.8;generated directivity;",ndim, cldNbins, cldMin, cldMax);
  THnSparseD *hMultEst = new THnSparseD("hMultEst","hMultEst;combined 0.8;V0M;V0A signal;directivity;LP pT; LP y;generated 0.8;generated V0M;generated 0.4-0.8;generated directivity;generated LP pT;generated LP y",ndim, cldNbins, cldMin, cldMax);
  hMultEst->GetAxis(1)->Set(102, binLimitsV0M);  // custom made V0M binning (to incorporate the 3 bins below 1)
  // hMultEst->GetAxis(5)->Set(102, binLimitsV0M);  // custom made V0M binning (to incorporate the 3 bins below 1)
  hMultEst->GetAxis(4)->Set(52, binLimits);
  hMultEst->GetAxis(10)->Set(52, binLimits);
  fHistosQA->AddAt(hMultEst, slot_MultEst);

  // used for raw spectra and a lot of corrections
  const Int_t ndimAllESD(14);
  const Int_t cldNbinsAllESD[ndimAllESD]   = {150, 102, 21, 52, 2, 5, 5, 20, 2, 80, 20, 5, 20, 2};
  const Double_t cldMinAllESD[ndimAllESD]  = {0.5, 0., 0., 0., -2., -0.5, -0.5, -1., -0.5, -TMath::PiOver2(), -2., -0.5, -1., -0.5},
  cldMaxAllESD[ndimAllESD]  = {150.5, 100., 1.05, 5., 2., 4.5, 4.5, 1., 1.5, (3.*TMath::PiOver2()), 2., 4.5, 1.,1.5};
  // THnSparseD *hAllESD = new THnSparseD("AllESD","AllESD;combined08;V0M;combined0408;p_{T};charge;PID_TPC;PID_TPCTOF;y;TOFmatching;MCPID;yMCPID;MCprimary;",ndimAllESD, cldNbinsAllESD, cldMinAllESD, cldMaxAllESD);
  THnSparseD *hAllESD = new THnSparseD("AllESD","AllESD;combined08;V0M;directivity;p_{T};charge;PID_TPC;PID_TPCTOF;y;TOFmatching;delta_phi;delta_y;MCPID;yMCPID;MCprimary;",ndimAllESD, cldNbinsAllESD, cldMinAllESD, cldMaxAllESD);
  hAllESD->GetAxis(1)->Set(102, binLimitsV0M);
  hAllESD->GetAxis(3)->Set(52, binLimits);
  fHistosQA->AddAt(hAllESD, slot_AllESD);

  // used for tracking efficiency
  const Int_t ndimGen(9);
  const Int_t cldNbinsGen[ndimGen]   = {150, 21, 52, 2, 5, 20, 80, 150, 20};
  const Double_t cldMinGen[ndimGen]  = {0.5, 0., 0., -2., -0.5, -1., -TMath::PiOver2(), 0.5, 0.},
  cldMaxGen[ndimGen]  = {150.5, 1.05, 5., 2., 4.5, 1., (3.*TMath::PiOver2()), 150.5, 1.};
  THnSparseD *hGen = new THnSparseD("Gen","Gen;MCmultiplicity;MCdirectivity;MCp_{T};MCcharge;MCPID;MCy;MCdelta_phi;ESDmultiplicity;ESDdirectivity;",ndimGen, cldNbinsGen, cldMinGen, cldMaxGen);
  hGen->GetAxis(2)->Set(52, binLimits);
  fHistosQA->AddAt(hGen, slot_Gen);

  // 	TH1D *testCounter = new TH1D("testCounter","testCounter", 4, 0.5, 4.5);
  // 	fHistosQA->AddAt(testCounter, 3);

  // used for scaling
/*
  TH2D *fNoEvt = new TH2D("fNoEvt", "Number of processed events", 4, -0.5, 3.5, 150, 0.5, 150.5);
  // 	fNoEvtMB->GetXaxis()->SetBinLabel(1, "Calls UserExec");
  // 	fNoEvtMB->GetXaxis()->SetBinLabel(2, "Physics Selection");
  // 	fNoEvtMB->GetXaxis()->SetBinLabel(3, "Vertex");
  fNoEvt->GetXaxis()->SetBinLabel(1, "Tender OK");
  fNoEvt->GetXaxis()->SetBinLabel(2, "Pile-up Rejection");
  fNoEvt->GetXaxis()->SetBinLabel(3, "Vertex Cut");
  fNoEvt->GetXaxis()->SetBinLabel(4, "Analyzed");
  fHistosQA->AddAt(fNoEvt, 3);
*/
  const Int_t ndimNoEvts(7);
  const Int_t cldNbinsNoEvts[ndimNoEvts]   = {4, 150, 102, 21, 150, 102, 20};
  const Double_t cldMinNoEvts[ndimNoEvts]  = {-0.5, 0.5, 0., 0., 0.5, 0., 0.},
  cldMaxNoEvts[ndimNoEvts]  = {3.5, 150.5, 100., 1.05, 150.5, 100., 1.};
  THnSparseD *hNoEvts = new THnSparseD("NoEvts","NoEvts;step;combined 0.8;V0M;directivity;MCmultiplicity;MCV0M;MCdirectivity;",ndimNoEvts, cldNbinsNoEvts, cldMinNoEvts, cldMaxNoEvts);
  hNoEvts->GetAxis(0)->SetBinLabel(1, "Tender OK");
  hNoEvts->GetAxis(0)->SetBinLabel(2, "Pile-up Rejection");
  hNoEvts->GetAxis(0)->SetBinLabel(3, "Vertex Cut");
  hNoEvts->GetAxis(0)->SetBinLabel(4, "Analyzed");
  hNoEvts->GetAxis(2)->Set(102, binLimitsV0M);  // custom made V0M binning (to incorporate the 3 bins below 1)
  fHistosQA->AddAt(hNoEvts, slot_NoEvts);


  // PID "QA"
  const Int_t ndimPID(8);
  const Int_t cldNbinsPID[ndimPID]   = {52, 2, 5, 5, 2, 100, 100, 5};
  const Double_t cldMinPID[ndimPID]  = {0., -2., -0.5, -0.5, -0.5, 0., 0., -0.5},
  cldMaxPID[ndimPID]  = {5., 2., 4.5, 4.5, 1.5, 200., 1.1, 4.5};
  THnSparseD *hPIDQA = new THnSparseD("PIDQA","PIDQA;p;charge;PID_TPC;PID_TPCTOF;TOFmatching;dEdx;beta;MCPID",ndimPID, cldNbinsPID, cldMinPID, cldMaxPID);
  hPIDQA->GetAxis(0)->Set(52, binLimits);
  fHistosQA->AddAt(hPIDQA, slot_PIDQA);


  // deltaPhi studies
  const Int_t ndimPhi(9);
  const Int_t cldNbinsPhi[ndimPhi]   = {150, 21, 52, 32, 80, 52, 32, 80, 80};
  const Double_t cldMinPhi[ndimPhi]  = {-0.5, 0., 0., -1.6, -TMath::PiOver2(), 0., -1.6, -TMath::PiOver2(), -TMath::PiOver2()},
  cldMaxPhi[ndimPhi]  = {150.5, 1.05, 5., 1.6, (3.*TMath::PiOver2()), 5., 1.6, (3.*TMath::PiOver2()), (3.*TMath::PiOver2())};
  THnSparseD *hDeltaPhi = new THnSparseD("DeltaPhi","deltaPhi;combined 0.8;directivity;ptESD; delta_y_ESD;deltaPhiESD;ptMC;delta_y_MC;deltaPhiMC;deltaPhiESD_LPMC",ndimPhi, cldNbinsPhi, cldMinPhi, cldMaxPhi);
  // hDeltaPhi->GetAxis(0)->Set(52, binLimits);
  fHistosQA->AddAt(hDeltaPhi, slot_DeltaPhi);

/*
	// used for DCA corrections
	const Int_t ndimDCA(12);
	const Int_t cldNbinsDCA[ndimDCA]   = {150, 21, 52, 2, 5, 5, 20, 2, 80, 3, 121, 2};
	const Double_t cldMinDCA[ndimDCA]  = {0.5, 0., 0.,-2., -0.5, -0.5, -1., -0.5, -TMath::PiOver2(), -0.5, -3., -0.5},
	cldMaxDCA[ndimDCA]  = {150.5, 1.05, 5., 2., 4.5, 4.5, 1., 1.5, (3.*TMath::PiOver2()), 2.5, 3.,1.5};
	THnSparseD *hDCA = new THnSparseD("DCA","DCA;multiplicity;directivity;p_{T};charge;PID_TPC;PID_TPCTOF;y;TOFmatching;delta_phi;MCorigin;DCA;pass_DCAcut;",ndimDCA, cldNbinsDCA, cldMinDCA, cldMaxDCA);
	hDCA->GetAxis(2)->Set(52, binLimits);
	fHistosQA->AddAt(hDCA, slot_DCA);
*/

    // test a tree output
    // TTree *testTree = new TTree("testTree", "test Tree");
    // fHistosQA->AddAt(testTree, slot_testTree);

    fHistosQA->ls();

    return kTRUE;
}


Double_t AliMESpidTask::ComputeDeltaPhi(Double_t phi, Double_t phi_LP)
{
    Double_t delta_phi = phi - phi_LP;

    Double_t result = -9999.;

    if(TMath::Abs(delta_phi) > 0.0001){  // avoid LP
        if(delta_phi < -TMath::PiOver2()){
            result = delta_phi + TMath::TwoPi();
        }
        else if(delta_phi > (3*TMath::PiOver2())){
            result = delta_phi - TMath::TwoPi();
        }
        else{
            result = delta_phi;
        }
    }

    return result;

}
