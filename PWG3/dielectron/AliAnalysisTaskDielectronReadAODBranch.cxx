/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

////////////////////////////////////////////////////////////////////////////
//  Class to retrieve the branch of the dielectron candidates stored in   //
//  filtered AODs file (AliAOD.Dielectron.root). It is possible to        //
//  apply tighter cuts to candidates stored in the branch and do the      //
//  matching with MC truth.                                               //
////////////////////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TDatabasePDG.h"
#include "TF1.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskDielectronReadAODBranch.h"
#include "AliDielectronPair.h"

ClassImp(AliAnalysisTaskDielectronReadAODBranch)
	//________________________________________________________________________
	AliAnalysisTaskDielectronReadAODBranch::AliAnalysisTaskDielectronReadAODBranch():
		AliAnalysisTaskSE(),
		fOutput(0),
                fNtupleJPSI(0),
                fNentries(0),
		fInvMass(0),
		fInvMassNoCuts(0),
                fpsproperSignal(0),
		fpsproperSidebands(0),            
		fpsproperAll(0),                  
		fpsproperUnder(0),
                fpsproperUpper(0),
                fLxyVsPtleg1(0),
		fLxyVsPtleg2(0),
		fLxyVsPt(0),
		fMeeVsPt(0),
		fMeeVsLxy(0),
		fprimvtxZ(0),  
		fsecvtxZ(0),  
		fprimvtxX(0),  
		fsecvtxX(0),  
		fprimvtxY(0),  
		fsecvtxY(0),  
		fPt(0),
		fPtLeg1(0),
		fPtLeg2(0),
		fdEdxP(0),
		fHasMC(0),
		fobj(0),
		fobjMC(0),
                fPtCut(1.),
                fSpdFirstRequired(kFALSE),
                fClsTPC(90),
                fPairType(1),
                fPtJpsi(1.3),
                fInvMassSignalLimits(0),
                fInvMassSideBandsLimits(0)
{
	// Default constructor
}

//________________________________________________________________________
AliAnalysisTaskDielectronReadAODBranch::AliAnalysisTaskDielectronReadAODBranch(const char *name):
	AliAnalysisTaskSE(name),
	fOutput(0),
        fNtupleJPSI(0),
        fNentries(0),
	fInvMass(0),
        fInvMassNoCuts(0),
        fpsproperSignal(0),
	fpsproperSidebands(0),            
	fpsproperAll(0),                  
	fpsproperUnder(0),
        fpsproperUpper(0),
        fLxyVsPtleg1(0),
	fLxyVsPtleg2(0),
	fLxyVsPt(0),
	fMeeVsPt(0),
	fMeeVsLxy(0),
	fprimvtxZ(0), 
	fsecvtxZ(0),  
	fprimvtxX(0), 
	fsecvtxX(0),  
	fprimvtxY(0), 
	fsecvtxY(0),  
	fPt(0),
	fPtLeg1(0),
	fPtLeg2(0),
	fdEdxP(0),
        fHasMC(0),
	fobj(0),
	fobjMC(0),
        fPtCut(1.),
        fSpdFirstRequired(kFALSE),
        fClsTPC(90),
        fPairType(1),
        fPtJpsi(1.3),
        fInvMassSignalLimits(0),
        fInvMassSideBandsLimits(0) 
{
	// Default constructor
	DefineInput(0,TChain::Class());
	DefineOutput(1,TList::Class());  //My private output

        fInvMassSignalLimits = new Double_t[2]; 
        fInvMassSignalLimits[0] = 0.;  fInvMassSignalLimits[1] = 0.; 
        fInvMassSideBandsLimits = new Double_t[2]; 
        fInvMassSideBandsLimits[0] = 0.;  fInvMassSideBandsLimits[1] = 0.;
}

//________________________________________________________________________
AliAnalysisTaskDielectronReadAODBranch::~AliAnalysisTaskDielectronReadAODBranch()
{
	if(fOutput)   delete fOutput;
	if(fobj)      delete fobj;
	if(fobjMC)    delete fobjMC;
}

//________________________________________________________________________
void AliAnalysisTaskDielectronReadAODBranch::Init()
{
	// Initialization
	if(fDebug > 1) printf("AnalysisTaskReadAOD::Init() \n");
	return;
}

//________________________________________________________________________
void AliAnalysisTaskDielectronReadAODBranch::UserCreateOutputObjects()
{
	//
        // Create the output container
	//
	if(fDebug > 1) printf("AnalysisTaskReadAOD::UserCreateOutputObjects() \n");
	fOutput = new TList();
	fOutput->SetOwner();
	//// Histogram booking
        
        // invariant mass
        fInvMass = new TH1F("fInvMass","fInvMass",300,2.0,2.0+300*.04); // step 40MeV
        fInvMassNoCuts = new TH1F("fInvMass_no_cuts","fInvMass_no_cuts",125,0,125*0.04); //step 40MeV
        fNentries = new TH1F("numberOfevent","numbOfEvent",1,-0.5,0.5);
	
        //pseudoproper 
	fpsproperSignal = new TH1F("psproper_decay_length",Form("psproper_decay_length_distrib(AODcuts+#clustTPC>%d, pT(leg)>%fGeV, pT(J/#psi)>%fGeV/c, %f < M < %f);X [#mu m];Entries/40#mu m",fClsTPC,fPtCut,fPtJpsi,fInvMassSignalLimits[0],fInvMassSignalLimits[1]),150,-3000.,3000.);
	fpsproperSidebands = new TH1F("psproper_decay_length_sidebands",Form("psproper_decay_length_distrib_sidebands(AODcuts+#clustTPC>%d, pT(e+e-)>%f GeV,pT(J/#psi)>%f GeV/c,M> %f && M < %f );  X [#mu m]; Entries/40#mu m",fClsTPC,fPtCut,fPtJpsi,fInvMassSideBandsLimits[1],fInvMassSideBandsLimits[0]),150,-3000.,3000.);
	fpsproperAll = new TH1F("psproper_decay_length_all",Form("psproper_decay_length_distrib_all(AODcuts+#clustTPC>%d, pT(e+e-)>%f GeV,pT(J/#psi)>%f GeV/c);X [#mu m];Entries/15#mu m",fClsTPC,fPtCut,fPtJpsi),400,-3000.,3000.);
	fpsproperUnder = new TH1F("psproper_decay_length_under",Form("psproper_decay_length_distrib_under(AODcuts+#clustTPC>%d,pT(e+e-)>%f GeV,pT(J/#psi)>%f GeV/c);X [#mu m];Entries/15#mu m",fClsTPC,fPtCut,fPtJpsi),400,-3000.,3000.); 
        fpsproperUpper = new TH1F("psproper_decay_length_upper",Form("psproper_decay_length_distrib_upper(AODcuts+#clustTPC>%d,pT(e+e-)>%f GeV,pT(J/#psi)>%f GeV/c);X [#mu m];Entries/15#mu m",fClsTPC,fPtCut,fPtJpsi),400,-3000.,3000.);  
        //
	fLxyVsPtleg1 = new TH2F("Lxy_vs_Pt_leg1",Form("Lxy_vs_Pt_leg1_distrib(AODcuts+#clustTPC>%d,pT(e+e-)>%f GeV,pT(J/#psi)>%f GeV/c,%f<M<%f);p_{T}[GeV/c]; X",fClsTPC,fPtCut,fPtJpsi,fInvMassSignalLimits[0],fInvMassSignalLimits[1]), 700, 0., 7.,400, -3000.,3000.);
	fLxyVsPtleg2 = new TH2F("Lxy_vs_Pt_leg2",Form("Lxy_vs_Pt_leg2_distrib(AODcuts+#clustTPC>%d,pT(e+e-)>%f GeV,pT(J/#psi)>%f GeV/c,%f<M<%f);p_{T}[GeV/c]; X",fClsTPC,fPtCut,fPtJpsi,fInvMassSignalLimits[0],fInvMassSignalLimits[1]), 700, 0., 7.,400, -3000.,3000.);
	fLxyVsPt = new TH2F("Lxy_vs_Pt_jpsi",Form("Lxy_vs_Pt_jpsi_distrib(AODcuts+#clustTPC>%d,pT(e+e-)>%f GeV,pT(J/#psi)>%f GeV/c,%f<M<%f); p_{T}[GeV/c]; X",fClsTPC,fPtCut,fPtJpsi,fInvMassSignalLimits[0],fInvMassSignalLimits[1]), 700, 0., 7.,400, -3000.,3000.);
	fMeeVsPt = new TH2F("Mee_vs_Pt_jpsi",Form("Mee_vs_Pt_jpsi_distrib(AODcuts+#clustTPC>%d,pT(e+e-)>%f GeV,pT(J/#psi)>%f GeV/c,%f<M<%f); p_{T}[GeV/c]; M_{ee}",fClsTPC,fPtCut,fPtJpsi,fInvMassSignalLimits[0],fInvMassSignalLimits[1]), 700, 0., 7.,200, 1.99,4.1);
	fMeeVsLxy = new TH2F("Mee_vs_Lxy_jpsi",Form("Mee_vs_Lxy_jpsi_distrib(AODcuts+#clustTPC>%d,pT(e+e-)>%f GeV),pT(J/#psi)>%f GeV/c; X; M_{ee}",fClsTPC,fPtCut,fPtJpsi), 400, -3000, 3000.,200, 1.99,4.1);

	// QA plots
	fprimvtxZ = new TH1F("prim_vtx_Z",Form("prim_vtx_Z_distrib(AODcuts + #clustTPC>%d, pT(e+e-)>%f GeV, pT(J/#psi)>%f);  Z[cm]; Entries/20 mm",fClsTPC,fPtCut,fPtJpsi),1000,-10.,10.);
	fsecvtxZ = new TH1F("sec_vtx_Z",Form("sec_vtx_Z_distrib(AODcuts + #clustTPC>%d, pT(e+e-)>%f GeV, pT(J/#psi)>%f);  Z[cm]; Entries/20 mm",fClsTPC,fPtCut,fPtJpsi),1000,-10.,10.);
	fprimvtxX = new TH1F("prim_vtx_X",Form("prim_vtx_X_distrib(AODcuts + #clustTPC>%d, pT(e+e-)>%f GeV,pT(J/#psi)>%f);  X[cm]; Entries/mm",fClsTPC,fPtCut,fPtJpsi),1000,-0.5,0.5);
	fsecvtxX = new TH1F("sec_vtx_X",Form("sec_vtx_X_distrib(AODcuts + #clustTPC>%d, pT(e+e-)>%f GeV,pT(J/#psi)>%f);  X[cm]; Entries/20 mm",fClsTPC,fPtCut,fPtJpsi),100,-1.,1.);
	fprimvtxY = new TH1F("prim_vtx_Y",Form("prim_vtx_Y_distrib(AODcuts + #clustTPC>%d, pT(e+e-)>%f GeV,pT(J/#psi)>%f);  Y[cm]; Entries/mm",fClsTPC,fPtCut,fPtJpsi),1000,-0.5,0.5);
	fsecvtxY = new TH1F("sec_vtx_Y",Form("sec_vtx_Y_distrib(AODcuts + #clustTPC>%d, pT(e+e-)>%f GeV,pT(J/#psi)>%f);  Y[cm]; Entries/20 mm",fClsTPC,fPtCut,fPtJpsi),100,-1.,1.);
	//
	fPt = new TH1F("Pt(J/psi)",Form("Pt_Jpsi_distrib(AODcuts + #clustTPC>%d, pT(e+e-)>%f GeV, pT(J/#psi)>%f); p_{T}(J/#psi)[GeV/c]; Entries/25 MeV",fClsTPC,fPtCut,fPtJpsi),400,0.,10.);
	fPtLeg1 = new TH1F("Pt_leg1",Form("Pt_leg1_distrib(AODcuts + #clustTPC>%d, pT(e+e-) > %f GeV, pT(J/#psi)>%f);  p_{T}[GeV/c]; Entries/25 MeV",fClsTPC,fPtCut,fPtJpsi),400,0.,10.);
	fPtLeg2 = new TH1F("Pt_leg2",Form("Pt_leg2_distrib(AODcuts + #clustTPC>%d, pT(e+e-) > %f GeV, pT(J/#psi)>%f);  p_{T}[GeV/c]; Entries/25 MeV",fClsTPC,fPtCut,fPtJpsi),400,0.,10.);
	//dE/dx plots
	fdEdxP = new TH2F("dE/dx_vs_Ptpc",Form("dE/dx_vs_Ptpc(AODcuts + #clustTPC>%d, pT(e+e-) > %f GeV, pT(J/#psi)>%f)",fClsTPC,fPtCut,fPtJpsi),400,0.2,20.,200,0.,200.);
        
        //J/psi n-tuple (X and M values)
        fNtupleJPSI = new TNtuple("fNtupleJPSI","J/#psi pseudo-proper decay time & invariant mass","Xdecaytime:Mass");
 
	// add histograms to list
	fOutput->Add(fNentries);
	fOutput->Add(fInvMass);
        fOutput->Add(fInvMassNoCuts); 
	fOutput->Add(fpsproperSignal);
	fOutput->Add(fpsproperSidebands);
	fOutput->Add(fpsproperAll);
	fOutput->Add(fpsproperUnder);
        fOutput->Add(fpsproperUpper);
        fOutput->Add(fLxyVsPtleg1);
	fOutput->Add(fLxyVsPtleg2);
	fOutput->Add(fLxyVsPt);
	fOutput->Add(fMeeVsPt);
	fOutput->Add(fMeeVsLxy);
	fOutput->Add(fprimvtxZ);
	fOutput->Add(fsecvtxZ);
	fOutput->Add(fprimvtxX);
	fOutput->Add(fsecvtxX);
	fOutput->Add(fprimvtxY);
	fOutput->Add(fsecvtxY);
	fOutput->Add(fPt);
	fOutput->Add(fPtLeg1);
	fOutput->Add(fPtLeg2);
        fOutput->Add(fNtupleJPSI);
        fOutput->Add(fdEdxP);
	return;
}

//________________________________________________________________________
void AliAnalysisTaskDielectronReadAODBranch::UserExec(Option_t */*option*/)
{
	// Execute analysis for current event:
	AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
	Double_t vtxPrim[3] = {0.,0.,0.};

	AliAODVertex* primvtx = aod->GetPrimaryVertex();
	vtxPrim[0] = primvtx->GetX();
	vtxPrim[1] = primvtx->GetY();
	vtxPrim[2] = primvtx->GetZ();

	AliAODHandler* aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());

	TTree *aodTree = aodHandler->GetTree();
	if(!fobj) aodTree->SetBranchAddress("dielectrons",&fobj);

	if(fHasMC){ aodTree->SetBranchAddress("mcparticles",&fobjMC);
	if(!fobjMC) printf("AliAnalysisTaskDielectronReadAODBranch::UserExec: MC particles branch not found!\n"); }

	fNentries->Fill(0);
	aodTree->GetEvent(Entry());
	Int_t dgLabels[2] = {0,0};
	Int_t pdgDgJPSItoEE[2] = {11,11};

	// loop over candidates
	if(fobj) {
		TObjArray *objArr = (TObjArray*)fobj->UncheckedAt(fPairType);
		for(int j=0;j<objArr->GetEntriesFast();j++)
		{
			AliDielectronPair *pairObj = (AliDielectronPair*)objArr->UncheckedAt(j);
			fInvMassNoCuts->Fill(pairObj->M()); 
                        Double_t vtxSec[3] = {0.,0.,0.};
			fprimvtxX->Fill(vtxPrim[0]);
			fprimvtxY->Fill(vtxPrim[1]);
			fprimvtxZ->Fill(vtxPrim[2]);
			pairObj->XvYvZv(vtxSec);
			fsecvtxX->Fill(vtxSec[0]);
			fsecvtxY->Fill(vtxSec[1]);
			fsecvtxZ->Fill(vtxSec[2]);

			Double_t lxy = ((vtxSec[0]-vtxPrim[0])*(pairObj->Px()) + (vtxSec[1]-vtxPrim[1])*(pairObj->Py()))/pairObj->Pt();
			Double_t psProperDecayLength = lxy*(TDatabasePDG::Instance()->GetParticle(443)->Mass())/pairObj->Pt();

			AliAODTrack *trk = (AliAODTrack*)pairObj->GetFirstDaughter();
			AliAODTrack *trk1 = (AliAODTrack*)pairObj->GetSecondDaughter();
			if(!trk || !trk1) {printf("ERROR: daughter tracks not available\n"); continue;}
			dgLabels[0] = trk->GetLabel(); 
			dgLabels[1] = trk1->GetLabel(); 

			//check in case of MC analysis if candidate is a true J/psi->ee
			if(fHasMC && !MatchToMC(443,fobjMC,dgLabels,2,2,pdgDgJPSItoEE)) continue;
                           AliAODPid *pid  = trk->GetDetPid();
                                AliAODPid *pid1 = trk1->GetDetPid();
                                if(!pid || !pid1){
                                        if(fDebug>1) printf("No AliAODPid found\n");
                                        continue;
                                }
                                 // ptLeg cut
                                if((trk->Pt()<fPtCut) || (trk1->Pt()<fPtCut)) continue;
                           
                                // pt jpsi cut
                                if((pairObj->Pt())<fPtJpsi) continue; 
 
                              //spd first required
                              if(fSpdFirstRequired)
                              {
                               if((!trk->HasPointOnITSLayer(0)) || (!trk1->HasPointOnITSLayer(0))) continue;
                              }
                           //
                          if((trk->GetTPCNcls()) <= fClsTPC || (trk1->GetTPCNcls()) <= fClsTPC) continue;
                         	if(trk->Charge()==-1){
					fPtLeg1->Fill(trk->Pt());
				} else {fPtLeg2->Fill(trk->Pt());}
				if(trk1->Charge()==-1){
					fPtLeg1->Fill(trk1->Pt());
				} else {fPtLeg2->Fill(trk1->Pt());}
				//Fill dE/dx related plots (before PID cuts)
				fdEdxP->Fill(pid->GetTPCmomentum(),pid->GetTPCsignal());
				fdEdxP->Fill(pid1->GetTPCmomentum(),pid1->GetTPCsignal());
					fInvMass->Fill(pairObj->M());
				        fPt->Fill(pairObj->Pt());
                                 	fMeeVsLxy->Fill(10000*psProperDecayLength,pairObj->M());
					fMeeVsPt->Fill(pairObj->Pt(),pairObj->M());
					fpsproperAll->Fill(10000*psProperDecayLength);
                                           
                                           //psproper in the signal region 
					   if((pairObj->M())<fInvMassSignalLimits[1] && (pairObj->M())>fInvMassSignalLimits[0]){
                                         	fpsproperSignal->Fill(10000*psProperDecayLength);
						fLxyVsPt->Fill(pairObj->Pt(),10000*psProperDecayLength);  //jpsi 
						fLxyVsPtleg1->Fill(trk->Pt(),10000*psProperDecayLength);  //leg1 (warning: mixture of pos and neg tracks)
						fLxyVsPtleg2->Fill(trk1->Pt(),10000*psProperDecayLength); //leg2 (warning: mixture of pos and neg tracks)
					        fNtupleJPSI->Fill(10000*psProperDecayLength,pairObj->M()); // fill the N-tuple (imput of the minimization algorithm)
                                                 } 

                                        // psproper in the sidebands 
					if((pairObj->M())> fInvMassSideBandsLimits[1] || (pairObj->M())< fInvMassSideBandsLimits[0])
                                        fpsproperSidebands->Fill(10000*psProperDecayLength);
                                          
                                        // psproper in the lower sideband
				        if(pairObj->M()< fInvMassSideBandsLimits[0]) fpsproperUnder->Fill(10000*psProperDecayLength);
                                          
                                        // psproper in the upper sideband
                                        if(pairObj->M()> fInvMassSideBandsLimits[1]) fpsproperUpper->Fill(10000*psProperDecayLength);
                        } // end loop over candidates
		}// end loop over pair types
        // Post the data
	PostData(1,fOutput);
	return;
} 

//________________________________________________________________________
void AliAnalysisTaskDielectronReadAODBranch::Terminate(Option_t */*option*/)
{
	if(fDebug > 1) printf("AliAnalysisTaskDielectronReadAODBranch::Terminate() \n");
	return;
}

//___________________________________________________________________
Bool_t AliAnalysisTaskDielectronReadAODBranch::MatchToMC(Int_t pdgabs,const TObjArray *mcArray,
		const Int_t *dgLabels,Int_t ndg,
		Int_t ndgCk, const Int_t *pdgDg) const
{
	//
	// jpsi candidates association to MC truth
        // Check if this candidate is matched to a MC signal
	// If no, return kFALSE
	// If yes, return kTRUE
	// 
	Int_t labMom[2]={0,0};
	Int_t i,j,lab,labMother,pdgMother,pdgPart,labJPSIMother,pdgJPSIMother;
	AliAODMCParticle *part=0;
	AliAODMCParticle *mother=0;
	AliAODMCParticle *jPSImother=0;
	Double_t pxSumDgs=0.,pySumDgs=0.,pzSumDgs=0.;
	Bool_t pdgUsed[2]={kFALSE,kFALSE};

	// loop on daughter labels
	for(i=0; i<ndg; i++) {
		labMom[i]=-1;
		lab = dgLabels[i];
		if(lab<0) {
			printf("daughter with negative label %d\n",lab);
			return kFALSE;
		}
		part = (AliAODMCParticle*)mcArray->At(lab);
		if(!part) {
			printf("no MC particle\n");
			return kFALSE;
		}

		// check the PDG of the daughter, if requested
		if(ndgCk>0) {
			pdgPart=TMath::Abs(part->GetPdgCode());
			printf("pdg code of daughter %d == %d\n",i,pdgPart);
			for(j=0; j<ndg; j++) {
				if(!pdgUsed[j] && pdgPart==pdgDg[j]) {
					pdgUsed[j]=kTRUE;
					break;
				}
			}
		}

		// for the J/psi, check that the daughters are electrons
		if(pdgabs==443) {
			if(TMath::Abs(part->GetPdgCode())!=11) return kFALSE;
		}

		mother = part;
		while(mother->GetMother()>=0) {
			labMother=mother->GetMother();
			mother = (AliAODMCParticle*)mcArray->At(labMother);//get particle mother
			if(!mother) {
				printf("no MC mother particle\n");
				break;
			}
			pdgMother = TMath::Abs(mother->GetPdgCode());//get particle mother pdg code
			printf("pdg code of mother track == %d\n",pdgMother);
			if(pdgMother==pdgabs) {
				labJPSIMother=mother->GetMother();
				jPSImother = (AliAODMCParticle*)mcArray->At(labJPSIMother);//get J/psi mother
				if(jPSImother) {
					pdgJPSIMother = TMath::Abs(jPSImother->GetPdgCode()); //get J/psi mother pdg code
					if( pdgJPSIMother==511   || pdgJPSIMother==521   ||
							pdgJPSIMother==10511 || pdgJPSIMother==10521 ||
							pdgJPSIMother==513   || pdgJPSIMother==523   ||
							pdgJPSIMother==10513 || pdgJPSIMother==10523 ||
							pdgJPSIMother==20513 || pdgJPSIMother==20523 ||
							pdgJPSIMother==515   || pdgJPSIMother==525   ||
							pdgJPSIMother==531   || pdgJPSIMother==10531 ||
							pdgJPSIMother==533   || pdgJPSIMother==10533 ||
							pdgJPSIMother==20533 || pdgJPSIMother==535   ||
							pdgJPSIMother==541   || pdgJPSIMother==10541 ||
							pdgJPSIMother==543   || pdgJPSIMother==10543 ||
							pdgJPSIMother==20543 || pdgJPSIMother==545) return kFALSE;} //check if J/psi comes from B hadron

					labMom[i]=labMother;
					// keep sum of daughters' momenta, to check for mom conservation
					pxSumDgs += part->Px();
					pySumDgs += part->Py();
					pzSumDgs += part->Pz();
					break;
			} else if(pdgMother>pdgabs || pdgMother<10) {
				break;
			}
		}
		if(labMom[i]==-1) {printf("mother PDG not ok for this daughter\n"); return kFALSE;} // mother PDG not ok for this daughter
	} // end loop on daughters

	// check if the candidate is signal
	labMother=labMom[0];
	// all labels have to be the same and !=-1
	for(i=0; i<ndg; i++) {
		if(labMom[i]==-1)        return kFALSE;
		if(labMom[i]!=labMother) return kFALSE;
	}

	// check that all daughter PDGs are matched
	if(ndgCk>0) {
		for(i=0; i<ndg; i++) {
			if(pdgUsed[i]==kFALSE) return kFALSE;
		}
	}

	mother = (AliAODMCParticle*)mcArray->At(labMother);
	Double_t pxMother = mother->Px();
	Double_t pyMother = mother->Py();
	Double_t pzMother = mother->Pz();
	// within 0.1%
	if((TMath::Abs(pxMother-pxSumDgs)/(TMath::Abs(pxMother)+1.e-13)) > 0.00001 &&
			(TMath::Abs(pyMother-pySumDgs)/(TMath::Abs(pyMother)+1.e-13)) > 0.00001 &&
			(TMath::Abs(pzMother-pzSumDgs)/(TMath::Abs(pzMother)+1.e-13)) > 0.00001) return kFALSE;

	return kTRUE;
}
