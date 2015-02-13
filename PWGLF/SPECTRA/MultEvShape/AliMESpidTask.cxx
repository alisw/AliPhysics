#include "THnSparse.h"
#include "TTreeStream.h"
#include "AliLog.h"

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


  Double_t ESDmult = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb);
/*
  Double_t ESDmult = fEvInfo->GetMultiplicity(AliMESeventInfo::kV0M);
  if(ESDmult > 0.){
	  if(ESDmult < 0.01) ESDmult = 1;
	  else if(ESDmult < 0.1) ESDmult = 2;
	  else if(ESDmult < 1.) ESDmult = 3;
	  else if(ESDmult < 5.) ESDmult = 4;   // 1-5 bin
	  else if(ESDmult < 10.) ESDmult = 8;   // 5-10 bin
	  else if(ESDmult < 15.) ESDmult = 12;   // 10-15 bin
	  else if(ESDmult < 20.) ESDmult = 18;   // 15-20 bin
	  else if(ESDmult < 30.) ESDmult = 25;   // 20-30 bin
	  else if(ESDmult < 40.) ESDmult = 35;   // 30-40 bin
	  else if(ESDmult < 50.) ESDmult = 45;   // 40-50 bin
	  else if(ESDmult < 70.) ESDmult = 60;   // 50-70 bin
	  else ESDmult = 90;   // 70-100 bin
  }
*/

//     Double_t ESDmult = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb0408);


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

// 	AliInfo("Starting the tracks loop:");
// 	printf("tracks = %i\n",fTracks->GetEntries());

  THnSparseD *multEst = (THnSparseD*)fHistosQA->At(0);
  Double_t correl_val[3];
  correl_val[0] = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb);
  correl_val[1] = fEvInfo->GetMultiplicity(AliMESeventInfo::kV0M);
  correl_val[2] = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb0408);
  multEst->Fill(correl_val);

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


	Double_t val[10];

	THnSparseD *hAllESD = (THnSparseD*)fHistosQA->At(1);

	// ---------------------------
	// get ESD multiplicity
// 	val[0] = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb);
// 	val[0] = fEvInfo->GetMultiplicity(AliMESeventInfo::kV0M);
// 	val[0] = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb0408);
	val[0] = ESDmult;
// 	AliInfo(Form("mult ESD = %g",val[0]));

	// ---------------------------
	// get pT
	val[1] = t->Pt();

	// ---------------------------
	// get charge
	val[2] = t->Charge();


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
	val[3] = maxIndex;    // particle identification
// 	AliInfo(Form("maxIndex TPC = %i\n\n", maxIndex));

	// TOF
	maxIndex = -9999;
	const Double_t *prob_TOF = t->GetPID()->GetProb(AliMEStrackInfo::kTOF);
	maxProb = 0; // max probability
	for(Int_t i = 0; i<AliPID::kSPECIES; i++){
// 		AliInfo(Form("probTOF[%i]= %g", i, prob_TOF[i]));
		if(prob_TOF[i] > maxProb){
			maxProb = prob_TOF[i];
			maxIndex = i;
		}
	}
	val[4] = maxIndex;    // particle identification
// 		AliInfo(Form("maxIndex TOF = %i\n\n", maxIndex));


	// ---------------------------
	// compute y after PID
	Double_t  mass[AliPID::kSPECIES] = {0.00051, 0.10565, 0.13957, 0.49368, 0.93827};
	Double_t e = TMath::Sqrt(t->P()*t->P() + mass[maxIndex]*mass[maxIndex]);
// 	Double_t pz = t->P() * TMath::Cos(t->Theta());
	if( TMath::Abs(t->Pz()) != e ) val[5] = 0.5*TMath::Log((e + t->Pz())/(e - t->Pz()));
	else val[5] = -9999;
	if(TMath::Abs(val[5]) > 1.0) continue;

// 	AliInfo(Form("Pz = %g \te = %g \t e+Pz = %g \t e-Pz = %g\n", t->Pz(), e, (e + t->Pz()), (e - t->Pz())));

	// ---------------------------
	// get the TOF mistmatch
// 	val[6] = (((t->GetPID()->GetTOFmisProb())<0.01)?1:0);
	val[6] = (t->GetPID()->GetTOFmisProb());
// 	printf("misProb = %g\n",t->GetPID()->GetTOFmisProb());


	if( HasMCdata() ){ // run only on MC
		// ---------------------------
		// get the MC PDG code
		Int_t MC_identity  = -9999;

		if(tMC){
			const Double_t *MCpdg = tMC->GetPID()->GetProb(AliMEStrackInfo::kITS);
			for(Int_t i = 0; i<AliPID::kSPECIES; i++){
				// 			AliInfo(Form("MC: probITS[%i]= %g", i, MCpdg[i]));
				if( MCpdg[i] > 0 ){
					val[7] = i;
					MC_identity = i;
				}
			}
			// 		AliInfo(Form("MC PID = %g", val[4]));

			if( (tMC->HasOrigin(AliMEStrackInfo::kPrimary)) ) val[9] = 1.;
			else val[9] = 0.;
		}

		// ---------------------------
		// compute y using true PID
		Double_t eMC = TMath::Sqrt(t->P()*t->P() + mass[MC_identity]*mass[MC_identity]);
		// 	Double_t pz = t->P() * TMath::Cos(t->Theta());
		if( TMath::Abs(t->Pz()) != eMC ) val[8] = 0.5*TMath::Log((eMC + t->Pz())/(eMC - t->Pz()));
		else val[8] = -9999;
		if(TMath::Abs(val[8]) > 1.0) continue;
	}

	// ---------------------------
	// fill the hSparse
	hAllESD->Fill(val);


	// THIS IS USED ONLY WHEN RUNNING FOR DCA SHAPES (change AliMEStender.cxx line 146 kTRUE->kFALSE)
/*
	Double_t valDCA[10];

	THnSparseD *hDCA = (THnSparseD*)fHistosQA->At(4);

	valDCA[0] = val[0];		// multiplicity
	valDCA[1] = val[1];		// pT
	valDCA[2] = val[2];		// charge
	valDCA[3] = val[3];		// PID_TPC
	valDCA[4] = val[4];		// PID_TPCTOF
	valDCA[5] = val[5];		// y
	valDCA[6] = val[6];		// TOFmatching

	// fill DCA
	Double_t dca[2];  // 0 = xy; 1 = z
	t->GetDCA(dca);
	valDCA[8] = dca[0];
// 	AliInfo(Form("DCAxy() = %g", valDCA[8]));

	// check if it passes the DCA cut
	if( TMath::Abs(valDCA[8]) < (0.0182+0.0350/ TMath::Power(valDCA[1],1.01)) ){  // valDCA[1] = pT; valDCA[8] = DCAxy
				valDCA[9] = 1.;
	}
	else{
		valDCA[9] = 0.;
	}


	if( HasMCdata() ){ // run only on MC
		if(tMC){
			// get origin
			if( (tMC->HasOrigin(AliMEStrackInfo::kPrimary)) ){
				valDCA[7] = 0.;
			}
			else{
				if( (tMC->HasOrigin(AliMEStrackInfo::kSecondary)) ){
					valDCA[7] = 1.;
				}
				else{
					if( (tMC->HasOrigin(AliMEStrackInfo::kMaterial)) ){
						valDCA[7] = 2.;
					}
					else{
						valDCA[7] = -1.;
					}
				}
			}
		}
	}


	hDCA->Fill(valDCA);
	*/

	if( HasMCdata() ){ // run only on MC
		Double_t testVal[6];
		testVal[0] = t->Eta();
		testVal[1] = t->Phi();
		testVal[2] = t->GetPID()->GetTOFmisProb();
		// 		testVal[3] = tMC->P();
		// 		testVal[4] = t->Charge();
		// 		testVal[5] = tMC->Charge();

		// fill debug
		if(DebugLevel()>1){
			(*AliMESbaseTask::DebugStream()) << "pidTrk"
			<<"mult=" << val[0]
			<<"pT=" << val[1]
			<<"charge=" << val[2]
			<<"TPC_PID=" << val[3]
			<<"TOF_PID=" << val[4]
			<<"y=" << val[5]
			<<"mismatch=" << testVal[2]
			<<"MC_PID=" << val[7]
			<<"MC_primary=" << val[9]
			<<"MC_y=" << val[8]
			// ----
			<<"eta=" << testVal[0]
			<<"phi=" << testVal[1]
			// ----
// 			<< "dEta=" << dEta
// 			<< "dPhi=" << dPhi
			// 		<<"t.=" << t
			// 		<<"tMC.=" << tMC
			<< "\n";
		}
	}

  } // end ESD track loop



  if( HasMCdata() ){ // run only on MC

// 	AliInfo("\n\nNew event");

  	for(Int_t it(0); it<fMCtracks->GetEntries(); it++){
    	if(!(tMC = (AliMEStrackInfo*)fMCtracks->At(it))) continue;

    	//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	// review this for DCA fits shapes
		if( !(tMC->HasOrigin(AliMEStrackInfo::kPrimary)) ) continue;

		// add this as a extra dimension on the sparse
// 		if(TMath::Abs(tMC->Eta()) > 0.9) continue;
		if(TMath::Abs(tMC->Y()) > 1.) continue;
    	//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    	Double_t valMC[6];

		THnSparseD *hGen = (THnSparseD*)fHistosQA->At(2);

		// ---------------------------
		// get generated multiplicity
		valMC[0] = fMCevInfo->GetMultiplicity(AliMESeventInfo::kGlob08);

		// ---------------------------
		// get pT
		valMC[1] = tMC->Pt();
//  AliInfo(Form("pT = %g", valMC[0]));

		// ---------------------------
		// get charge
		valMC[2] = tMC->Charge();

		// ---------------------------
		// get the MC PDG code
		const Double_t *MCpdg = tMC->GetPID()->GetProb(AliMEStrackInfo::kITS);
		for(Int_t i = 0; i<AliPID::kSPECIES; i++){
// 			AliInfo(Form("MC: probITS[%i]= %g", i, MCpdg[i]));
			if( MCpdg[i] > 0 )	valMC[3] = i;
		}

		// ---------------------------
		// get y
		valMC[4] = tMC->Y();

		// ---------------------------
		// get the ESD multiplicity
		// 		valMC[5] = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb);
		valMC[5] = ESDmult;

		// ---------------------------
		// fill the hSparse
		hGen->Fill(valMC);

/*
		TH1D *testCounter = (TH1D*)fHistosQA->At(3);
		if(valMC[2] > 0){
			testCounter->Fill(1.);
			if(TMath::Abs(valMC[4]) < 0.5){
				testCounter->Fill(2.);
				if(valMC[3] == 2){
					AliInfo(Form("pT = %g \t charge = %g \t yMC = %g \t MC PID = %g", valMC[1], valMC[2], valMC[4], valMC[3]));
					testCounter->Fill(3.);
				}
			}
		}
*/

	} // end MC track loop

  } // end HasMCdata IF



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

  // example
  const Int_t ndim(3);
  const Int_t cldNbins[ndim]   = {150, 100, 150};
  const Double_t cldMin[ndim]  = {0.5, 0., 0.5},
  cldMax[ndim]  = {150.5, 100., 150.5};
  THnSparseD *multEst = new THnSparseD("multEst","multEst;combined 0.8;V0M;combined 0.4-0.8",ndim, cldNbins, cldMin, cldMax);
  fHistosQA->AddAt(multEst, 0);


	// use for matching, PID and contaminations efficiency
	Double_t binLimits[] = {0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0};

	// used for raw spectra and a lot of corrections
	const Int_t ndimAllESD(10);
	const Int_t cldNbinsAllESD[ndimAllESD]   = {150, 52, 2, 5, 5, 20, 2, 5, 20, 2};
	const Double_t cldMinAllESD[ndimAllESD]  = {0.5, 0.,-2., -0.5, -0.5, -1., -0.5, -0.5, -1., -0.5},
	cldMaxAllESD[ndimAllESD]  = {150.5, 5., 2., 4.5, 4.5, 1., 1.5, 4.5, 1.,1.5};
	THnSparseD *hAllESD = new THnSparseD("AllESD","AllESD;multiplicity;p_{T};charge;PID_TPC;PID_TPCTOF;y;TOFmatching;MCPID;yMCPID;MCprimary;",ndimAllESD, cldNbinsAllESD, cldMinAllESD, cldMaxAllESD);
	hAllESD->GetAxis(1)->Set(52, binLimits);
	fHistosQA->AddAt(hAllESD, 1);

	// used for tracking efficiency
	const Int_t ndimGen(6);
	const Int_t cldNbinsGen[ndimGen]   = {150, 52, 2, 5, 20, 150};
	const Double_t cldMinGen[ndimGen]  = {0.5, 0., -2., -0.5, -1., 0.5},
	cldMaxGen[ndimGen]  = {150.5, 5., 2., 4.5, 1.,150.5};
	THnSparseD *hGen = new THnSparseD("Gen","Gen;MCmultiplicity;MCp_{T};MCcharge;MCPID;MCy;ESDmultiplicity;",ndimGen, cldNbinsGen, cldMinGen, cldMaxGen);
	hGen->GetAxis(1)->Set(52, binLimits);
	fHistosQA->AddAt(hGen, 2);

// 	TH1D *testCounter = new TH1D("testCounter","testCounter", 4, 0.5, 4.5);
// 	fHistosQA->AddAt(testCounter, 3);

	// used for scaling
	TH2D *fNoEvt = new TH2D("fNoEvt", "Number of processed events", 4, -0.5, 3.5, 150, 0.5, 150.5);
// 	fNoEvtMB->GetXaxis()->SetBinLabel(1, "Calls UserExec");
// 	fNoEvtMB->GetXaxis()->SetBinLabel(2, "Physics Selection");
	// 	fNoEvtMB->GetXaxis()->SetBinLabel(3, "Vertex");
	fNoEvt->GetXaxis()->SetBinLabel(1, "Tender OK");
	fNoEvt->GetXaxis()->SetBinLabel(2, "Pile-up Rejection");
	fNoEvt->GetXaxis()->SetBinLabel(3, "Vertex Cut");
	fNoEvt->GetXaxis()->SetBinLabel(4, "Analyzed");
	fHistosQA->AddAt(fNoEvt, 3);

/*
	// used for DCA corrections
	const Int_t ndimDCA(10);
	const Int_t cldNbinsDCA[ndimDCA]   = {100, 52, 2, 5, 5, 20, 2, 3, 121, 2};
	const Double_t cldMinDCA[ndimDCA]  = {0.5, 0.,-2., -0.5, -0.5, -1., -0.5, -0.5, -3., -0.5},
	cldMaxDCA[ndimDCA]  = {100.5, 5., 2., 4.5, 4.5, 1., 1.5, 2.5, 3.,1.5};
	THnSparseD *hDCA = new THnSparseD("DCA","DCA;multiplicity;p_{T};charge;PID_TPC;PID_TPCTOF;y;TOFmatching;MCorigin;DCA;pass_DCAcut;",ndimDCA, cldNbinsDCA, cldMinDCA, cldMaxDCA);
	hDCA->GetAxis(1)->Set(52, binLimits);
	fHistosQA->AddAt(hDCA, 4);
	*/

  return kTRUE;
}
