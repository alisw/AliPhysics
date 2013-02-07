/*************************************************************************************
 *	Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.	*
 *											*
 *	Author: The ALICE Off-line Project.						*
 *	Contributors are mentioned in the code where appropriate.			*
 *											*
 *	Permission to use, copy, modify and distribute this software and its		*
 *	documentation strictly for non-commercial purposes is hereby granted		*
 *	without fee, provided that the above copyright notice appears in all		*
 *	copies and that both the copyright notice and this permission notice		*
 *	appear in the supporting documentation. The authors make no claims		*
 *	about the suitability of this software for any purpose. It is			*
 *	provided "as is" without express or implied warranty.				*
 *											*
 *************************************************************************************/

/*************************************************************************************
 *											*
 *	Class for the Selection of Non-Heavy-Flavour-Electrons trought		*
 *	the invariant mass method. The selection can be done from two			*
 *	different algorithms, which can be choosed calling the function		*
 *		"SetAlgorithm(TString Algorithm)".	  				*
 *											*
 *		Authors: R.Bailhache, C.A.Schmidt					*
 *											*
 *************************************************************************************/

#include "TVector2.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TList.h"
#include "TDatabasePDG.h"

#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliPID.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"

#include "AliHFENonPhotonicElectron.h"

ClassImp(AliHFENonPhotonicElectron)
//________________________________________________________________________
AliHFENonPhotonicElectron::AliHFENonPhotonicElectron(const char *name, const Char_t *title)
  :TNamed		(name, title)
  ,fIsAOD		(kFALSE)
  ,fMCEvent		(NULL)
  ,fAODArrayMCInfo	(NULL)
  ,fHFEBackgroundCuts	(NULL)
  ,fPIDBackground	(0x0)
  ,fPIDBackgroundQA	(0)
  ,fkPIDRespons		(NULL)
  ,fAlgorithmMA		(kTRUE)
  ,fUseFilterAOD	(kTRUE)
  ,fFilter		(-1)
  ,fChi2OverNDFCut	(3.0)
  ,fMaxDCA		(3.0)
//  ,fMaxOpeningTheta	(0.02)
//  ,fMaxOpeningPhi	(0.1)
  ,fMaxOpening3D	(TMath::Pi())
  ,fMaxInvMass		(1000)
  ,fSetMassConstraint	(kFALSE)
  ,fArraytrack		(NULL)
  ,fCounterPoolBackground	(0)
  ,fnumberfound			(0)
  ,fListOutput		(NULL)
  ,fAssElectron		(NULL)
  ,fIncElectron		(NULL)
  ,fUSign		(NULL)
  ,fLSign		(NULL)
//  ,fUSignAngle	(NULL)
//  ,fLSignAngle	(NULL)
{
  //
  // Constructor
  //
  fPIDBackground   = new AliHFEpid("hfePidBackground");
  fPIDBackgroundQA = new AliHFEpidQAmanager;
}

//________________________________________________________________________
AliHFENonPhotonicElectron::AliHFENonPhotonicElectron()
  :TNamed		()
  ,fIsAOD		(kFALSE)
  ,fMCEvent		(NULL)
  ,fAODArrayMCInfo	(NULL)
  ,fHFEBackgroundCuts	(NULL)
  ,fPIDBackground	(0x0)
  ,fPIDBackgroundQA	(0)
  ,fkPIDRespons		(NULL)
  ,fAlgorithmMA		(kTRUE)
  ,fUseFilterAOD	(kTRUE)
  ,fFilter		(-1)
  ,fChi2OverNDFCut	(3.0)
  ,fMaxDCA		(3.0)
//  ,fMaxOpeningTheta	(0.02)
//  ,fMaxOpeningPhi	(0.1)
  ,fMaxOpening3D	(TMath::TwoPi())
  ,fMaxInvMass		(1000)
  ,fSetMassConstraint	(kFALSE)
  ,fArraytrack		(NULL)
  ,fCounterPoolBackground	(0)
  ,fnumberfound			(0)
  ,fListOutput		(NULL)
  ,fAssElectron		(NULL)
  ,fIncElectron		(NULL)
  ,fUSign		(NULL)
  ,fLSign		(NULL)
//  ,fUSignAngle	(NULL)
//  ,fLSignAngle	(NULL)
{
  //
  // Constructor
  //
  fPIDBackground   = new AliHFEpid("hfePidBackground");
  fPIDBackgroundQA = new AliHFEpidQAmanager;
}

//________________________________________________________________________
AliHFENonPhotonicElectron::AliHFENonPhotonicElectron(const AliHFENonPhotonicElectron &ref)
  :TNamed(ref)
  ,fIsAOD		(ref.fIsAOD)
  ,fMCEvent		(NULL)
  ,fAODArrayMCInfo	(NULL)
  ,fHFEBackgroundCuts	(ref.fHFEBackgroundCuts)
  ,fPIDBackground	(ref.fPIDBackground)
  ,fPIDBackgroundQA	(ref.fPIDBackgroundQA)
  ,fkPIDRespons		(ref.fkPIDRespons)
  ,fAlgorithmMA		(ref.fAlgorithmMA)
  ,fUseFilterAOD	(ref.fUseFilterAOD)
  ,fFilter		(ref.fFilter)
  ,fChi2OverNDFCut	(ref.fChi2OverNDFCut)
  ,fMaxDCA		(ref.fMaxDCA)
//  ,fMaxOpeningTheta	(ref.fMaxOpeningTheta)
//  ,fMaxOpeningPhi	(ref.fMaxOpeningPhi)
  ,fMaxOpening3D	(ref.fMaxOpening3D)
  ,fMaxInvMass		(ref.fMaxInvMass)
  ,fSetMassConstraint	(ref.fSetMassConstraint)
  ,fArraytrack		(NULL)
  ,fCounterPoolBackground	(0)
  ,fnumberfound			(0)
  ,fListOutput		(ref.fListOutput)
  ,fAssElectron		(ref.fAssElectron)
  ,fIncElectron		(ref.fIncElectron)
  ,fUSign		(ref.fUSign)
  ,fLSign		(ref.fLSign)
//  ,fUSignAngle	(ref.fUSignAngle)
//  ,fLSignAngle	(ref.fLSignAngle)
{
  //
  // Copy Constructor
  //
  ref.Copy(*this);
}

//____________________________________________________________
AliHFENonPhotonicElectron &AliHFENonPhotonicElectron::operator=(const AliHFENonPhotonicElectron &ref){
  //
  // Assignment operator
  //
  if(this == &ref) ref.Copy(*this);
  return *this;
}

//_________________________________________
AliHFENonPhotonicElectron::~AliHFENonPhotonicElectron()
{
  //
  // Destructor
  //
  if(fArraytrack)		delete fArraytrack;
  //if(fHFEBackgroundCuts)	delete fHFEBackgroundCuts;
  if(fPIDBackground)		delete fPIDBackground;
  if(fPIDBackgroundQA)		delete fPIDBackgroundQA;
}

//_____________________________________________________________________________________________
void AliHFENonPhotonicElectron::Init()
{
  //
  // Init
  //

  //printf("Analysis Mode for AliHFENonPhotonicElectron: %s Analysis\n", fIsAOD ? "AOD" : "ESD");

  if(!fListOutput) fListOutput = new TList;
  fListOutput->SetName("HFENonPhotonicElectron");
  fListOutput->SetOwner();

  if(!fHFEBackgroundCuts) fHFEBackgroundCuts = new AliHFEcuts();
  if(fIsAOD) fHFEBackgroundCuts->SetAOD();
  fHFEBackgroundCuts->Initialize();
  if(fHFEBackgroundCuts->IsQAOn()) {
    fListOutput->Add(fHFEBackgroundCuts->GetQAhistograms());
  }

  // Initialize PID
  if(!fPIDBackground) fPIDBackground = new AliHFEpid("default pid");
  if(fMCEvent || fAODArrayMCInfo) fPIDBackground->SetHasMCData(kTRUE); // does nothing since the fMCEvent are set afterwards at the moment
  if(!fPIDBackground->GetNumberOfPIDdetectors())
  {
    //fPIDBackground->AddDetector("TOF", 0);
    fPIDBackground->AddDetector("TPC", 0);
  }
  AliInfo("PID Background QA switched on");
  fPIDBackgroundQA->Initialize(fPIDBackground);
  fListOutput->Add(fPIDBackgroundQA->MakeList("HFENP_PID_Background"));
  fPIDBackground->SortDetectors();

  Int_t nBinsPt = 35;
  //Double_t binLimPt[25] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 3., 3.5, 4., 5., 6.};
  Double_t binLimPt[36] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};

  Int_t nBinsP = 400;
  Double_t minP = 0.0;
  Double_t maxP = 20.0;
  Double_t binLimP[nBinsP+1];
  for(Int_t i=0; i<=nBinsP; i++) binLimP[i]=(Double_t)minP + (maxP-minP)/nBinsP*(Double_t)i ;

  Int_t nBinsC = 11;
  Double_t minC = 0.0;
  Double_t maxC = 11.0;
  Double_t binLimC[nBinsC+1];
  for(Int_t i=0; i<=nBinsC; i++) binLimC[i]=(Double_t)minC + (maxC-minC)/nBinsC*(Double_t)i ;

  Int_t nBinsSource = 10;
  Double_t minSource = 0.;
  Double_t maxSource = 10.;
  Double_t binLimSource[nBinsSource+1];
  for(Int_t i=0; i<=nBinsSource; i++) binLimSource[i]=(Double_t)minSource + (maxSource-minSource)/nBinsSource*(Double_t)i ;

  Int_t nBinsInvMass = 1000;
  Double_t minInvMass = 0.;
  Double_t maxInvMass = 10.;
  Double_t binLimInvMass[nBinsInvMass+1];
  for(Int_t i=0; i<=nBinsInvMass; i++) binLimInvMass[i]=(Double_t)minInvMass + (maxInvMass-minInvMass)/nBinsInvMass*(Double_t)i ;

  Int_t nBinsPhi = 180;
  Double_t minPhi = 0.0;
  Double_t maxPhi = TMath::Pi();
  Double_t binLimPhi[nBinsPhi+1];
  for(Int_t i=0; i<=nBinsPhi; i++)
  {
    binLimPhi[i]=(Double_t)minPhi + (maxPhi-minPhi)/nBinsPhi*(Double_t)i ;
    AliDebug(2,Form("bin phi is %f for %d",binLimPhi[i],i));
  }

  Int_t nBinsAngle = 180;
  Double_t minAngle = 0.0;
  Double_t maxAngle = TMath::Pi();
  Double_t binLimAngle[nBinsAngle+1];
  for(Int_t i=0; i<=nBinsAngle; i++)
  {
    binLimAngle[i]=(Double_t)minAngle + (maxAngle-minAngle)/nBinsAngle*(Double_t)i ;
    AliDebug(2,Form("bin phi is %f for %d",binLimAngle[i],i));
  }

  Int_t nBinsTPC = 400;
  Double_t minTPC = 20;
  Double_t maxTPC = 220;
  Double_t binLimTPC[nBinsTPC+1];
  for(Int_t i=0; i<=nBinsTPC; i++)
  {
    binLimTPC[i]=(Double_t)minTPC + (maxTPC-minTPC)/nBinsTPC*(Double_t)i ;
    AliDebug(2,Form("bin TPC is %f for %d",binLimTPC[i],i));
  }

  Int_t nBinsTPCSigma = 240;
  Double_t minTPCSigma = -12.0;
  Double_t maxTPCSigma =  12.0;
  Double_t binLimTPCSigma[nBinsTPCSigma+1];
  for(Int_t i=0; i<=nBinsTPCSigma; i++) binLimTPCSigma[i]=(Double_t)minTPCSigma + (maxTPCSigma-minTPCSigma)/nBinsTPCSigma*(Double_t)i ;

  Int_t nBinsTOFSigma = 240;
  Double_t minTOFSigma = -12.0;
  Double_t maxTOFSigma =  12.0;
  Double_t binLimTOFSigma[nBinsTOFSigma+1];
  for(Int_t i=0; i<=nBinsTOFSigma; i++) binLimTOFSigma[i]=(Double_t)minTOFSigma + (maxTOFSigma-minTOFSigma)/nBinsTOFSigma*(Double_t)i ;

  // Associated Electron
  const Int_t nDimAssElectron=7;
  Int_t nBinAssElectron[nDimAssElectron] = {nBinsC,nBinsPt,nBinsSource,nBinsP,nBinsTPC,nBinsTPCSigma,nBinsTOFSigma};
  fAssElectron = new THnSparseF("fAssElectron","fAssElectron",nDimAssElectron,nBinAssElectron);
  fAssElectron->SetBinEdges(0,binLimC);
  fAssElectron->SetBinEdges(1,binLimPt);
  fAssElectron->SetBinEdges(2,binLimSource);
  fAssElectron->SetBinEdges(3,binLimP);
  fAssElectron->SetBinEdges(4,binLimTPC);
  fAssElectron->SetBinEdges(5,binLimTPCSigma);
  fAssElectron->SetBinEdges(6,binLimTOFSigma);
  fAssElectron->Sumw2();
  AliDebug(2,"AliHFENonPhotonicElectron: fAssElectron");

  // Inclusive Electron
  const Int_t nDimIncElectron=7;
  Int_t nBinIncElectron[nDimIncElectron] = {nBinsC,nBinsPt,nBinsSource,nBinsP,nBinsTPC,nBinsTPCSigma,nBinsTOFSigma};
  fIncElectron = new THnSparseF("fIncElectron","fIncElectron",nDimIncElectron,nBinIncElectron);
  fIncElectron->SetBinEdges(0,binLimC);
  fIncElectron->SetBinEdges(1,binLimPt);
  fIncElectron->SetBinEdges(2,binLimSource);
  fIncElectron->SetBinEdges(3,binLimP);
  fIncElectron->SetBinEdges(4,binLimTPC);
  fIncElectron->SetBinEdges(5,binLimTPCSigma);
  fIncElectron->SetBinEdges(6,binLimTOFSigma);
  fIncElectron->Sumw2();
  AliDebug(2,"AliHFENonPhotonicElectron: fIncElectron");

  // ee invariant mass Unlike Sign
  const Int_t nDimUSign=7;
  Int_t nBinUSign[nDimUSign] = {nBinsPhi,nBinsC,nBinsPt,nBinsInvMass,nBinsSource,nBinsAngle,nBinsPt};
  fUSign = new THnSparseF("fUSign","fUSign",nDimUSign,nBinUSign);
  fUSign->SetBinEdges(0,binLimPhi);
  fUSign->SetBinEdges(1,binLimC);
  fUSign->SetBinEdges(2,binLimPt);
  fUSign->SetBinEdges(3,binLimInvMass);
  fUSign->SetBinEdges(4,binLimSource);
  fUSign->SetBinEdges(5,binLimAngle);
  fUSign->SetBinEdges(6,binLimPt);
  fUSign->Sumw2();
  AliDebug(2,"AliHFENonPhotonicElectron: fUSign");

  // ee invariant mass Like Sign
  const Int_t nDimLSign=7;
  Int_t nBinLSign[nDimLSign] = {nBinsPhi,nBinsC,nBinsPt,nBinsInvMass,nBinsSource,nBinsAngle,nBinsPt};
  fLSign = new THnSparseF("fLSign","fLSign",nDimLSign,nBinLSign);
  fLSign->SetBinEdges(0,binLimPhi);
  fLSign->SetBinEdges(1,binLimC);
  fLSign->SetBinEdges(2,binLimPt);
  fLSign->SetBinEdges(3,binLimInvMass);
  fLSign->SetBinEdges(4,binLimSource);
  fLSign->SetBinEdges(5,binLimAngle);
  fLSign->SetBinEdges(6,binLimPt);
  fLSign->Sumw2();
  AliDebug(2,"AliHFENonPhotonicElectron: fLSign");

/*
  // ee angle Unlike Sign
  const Int_t nDimUSignAngle=3;
  Int_t nBinUSignAngle[nDimUSignAngle] = {nBinsAngle,nBinsC,nBinsSource};
  fUSignAngle = new THnSparseF("fUSignAngle","fUSignAngle",nDimUSignAngle,nBinUSignAngle);
  fUSignAngle->SetBinEdges(0,binLimAngle);
  fUSignAngle->SetBinEdges(1,binLimC);
  fUSignAngle->SetBinEdges(2,binLimSource);
  fUSignAngle->Sumw2();
  AliDebug(2,"AliHFENonPhotonicElectron: fUSignAngle");

  // ee angle Like Sign
  const Int_t nDimLSignAngle=3;
  Int_t nBinLSignAngle[nDimLSignAngle] = {nBinsAngle,nBinsC,nBinsSource};
  fLSignAngle = new THnSparseF("fLSignAngle","fLSignAngle",nDimLSignAngle,nBinLSignAngle);
  fLSignAngle->SetBinEdges(0,binLimAngle);
  fLSignAngle->SetBinEdges(1,binLimC);
  fLSignAngle->SetBinEdges(2,binLimSource);
  fLSignAngle->Sumw2();
  AliDebug(2,"AliHFENonPhotonicElectron: fLSignAngle");
*/

  fListOutput->Add(fAssElectron);
  fListOutput->Add(fIncElectron);
  fListOutput->Add(fUSign);
  fListOutput->Add(fLSign);
//  fListOutput->Add(fUSignAngle);
//  fListOutput->Add(fLSignAngle);

}

//_____________________________________________________________________________________________
void AliHFENonPhotonicElectron::InitRun(const AliVEvent *inputEvent,const AliPIDResponse *pidResponse)
{
  //
  // Init run
  //

  if(!pidResponse)
  {
    AliDebug(1, "Using default PID Response");
    Bool_t hasmc = kFALSE;
    if(fMCEvent || fAODArrayMCInfo) hasmc=kTRUE;
    pidResponse = AliHFEtools::GetDefaultPID(hasmc, inputEvent->IsA() == AliESDEvent::Class());
  }

  if(!fPIDBackground) return;
  fPIDBackground->SetPIDResponse(pidResponse);

  if(!fPIDBackground->IsInitialized())
  {
    // Initialize PID with the given run number
    fPIDBackground->InitializePID(inputEvent->GetRunNumber());
  }

}

//_____________________________________________________________________________________________
Int_t AliHFENonPhotonicElectron::FillPoolAssociatedTracks(AliVEvent *inputEvent, Int_t binct)
{
  //
  // Fill the pool of associated tracks
  // Return the number of associated tracks
  //

  fnumberfound = 0;

  Int_t nbtracks = inputEvent->GetNumberOfTracks();

  if( fArraytrack )
  {
    fArraytrack->~TArrayI();
    new(fArraytrack) TArrayI(nbtracks);
  }
  else
  {
    fArraytrack = new TArrayI(nbtracks);
  }

  fCounterPoolBackground = 0;

  for(Int_t k = 0; k < nbtracks; k++)
  {
    AliVTrack *track = (AliVTrack *) inputEvent->GetTrack(k);
    if(!track) continue;

    // Track cuts
    Bool_t survivedbackground = kTRUE;
    AliAODEvent *aodeventu = dynamic_cast<AliAODEvent *>(inputEvent);

    if(aodeventu)
    {
      /**				**
       *	AOD Analysis		 *
       **				**/

      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack *>(track);
      if(aodtrack)
      {
	// filter
	if(fUseFilterAOD)
	{
	  if(!(aodtrack->TestFilterBit(fFilter))) survivedbackground = kFALSE;
	}

      }
    }
   
    if(!fHFEBackgroundCuts->CheckParticleCuts(AliHFEcuts::kStepRecKineITSTPC + AliHFEcuts::kNcutStepsMCTrack, (TObject *) track)) survivedbackground = kFALSE;
    if(!fHFEBackgroundCuts->CheckParticleCuts(AliHFEcuts::kStepRecPrim       + AliHFEcuts::kNcutStepsMCTrack, (TObject *) track)) survivedbackground = kFALSE;

    // PID
    if(survivedbackground)
    {
      // PID track cuts
      AliHFEpidObject hfetrack2;

      if(!aodeventu)	hfetrack2.SetAnalysisType(AliHFEpidObject::kESDanalysis);
      else 		hfetrack2.SetAnalysisType(AliHFEpidObject::kAODanalysis);

      hfetrack2.SetRecTrack(track);
      if(binct>-1)
      {
	hfetrack2.SetCentrality((Int_t)binct);
	AliDebug(2,Form("centrality %d and %d",binct,hfetrack2.GetCentrality()));
	hfetrack2.SetPbPb();
      }

      if(fPIDBackground->IsSelected(&hfetrack2,0x0,"recTrackCont",fPIDBackgroundQA))
      {
	fArraytrack->AddAt(k,fCounterPoolBackground);
	fCounterPoolBackground++;
	AliDebug(2,Form("fCounterPoolBackground %d, track %d",fCounterPoolBackground,k));
      }
    }
  } // loop tracks

  //printf(Form("Associated Pool: Tracks %d, fCounterPoolBackground %d \n", nbtracks, fCounterPoolBackground));

  return fCounterPoolBackground;

}

//_____________________________________________________________________________________________
Int_t AliHFENonPhotonicElectron::CountPoolAssociated(AliVEvent *inputEvent, Int_t binct)
{
  //
  // Count the pool of assiocated tracks
  //


  if(fnumberfound > 0) //!count only events with an inclusive electron
  {
    Double_t valueAssElectron[7] = { binct, -1, -1, -1, -1, -20, -20};		//Centrality	Pt	Source	P	TPCsignal	TPCsigma	TOFsigma
    Int_t iTrack2 = 0;
    Int_t indexmother2 = -1;
    AliVTrack *track2 = 0x0;

    for(Int_t ii = 0; ii < fCounterPoolBackground; ii++)
    {
      iTrack2 = fArraytrack->At(ii);
      AliDebug(2,Form("track %d",iTrack2));
      track2 = (AliVTrack *)inputEvent->GetTrack(iTrack2);

      if(!track2)
      {
	//printf("ERROR: Could not receive track %d", iTrack2);
	continue;
      }

      // if MC look
      if(fMCEvent || fAODArrayMCInfo)
      {
	valueAssElectron[2] = FindMother(TMath::Abs(track2->GetLabel()), indexmother2) ;
      }

      fkPIDRespons = fPIDBackground->GetPIDResponse();

      valueAssElectron[1] = track2->Pt() ;
      valueAssElectron[3] = track2->P() ;
      valueAssElectron[4] = track2->GetTPCsignal() ;
      valueAssElectron[5] = fkPIDRespons->NumberOfSigmasTPC( track2, AliPID::kElectron) ;
      valueAssElectron[6] = fkPIDRespons->NumberOfSigmasTOF( track2, AliPID::kElectron) ;

      fAssElectron->Fill( valueAssElectron) ;
    }
  //printf(Form("Associated Pool: fCounterPoolBackground %d \n", fCounterPoolBackground));
  }
  return fnumberfound;
}

//_____________________________________________________________________________________________
Int_t AliHFENonPhotonicElectron::LookAtNonHFE(Int_t iTrack1, AliVTrack *track1, AliVEvent *vEvent, Double_t weight, Int_t binct, Double_t deltaphi, Int_t source, Int_t indexmother)
{
  //
  // Look At Non HFE
  //

  /***********************************************************************************
   *											*
   *	iTrack1:	index of the tagged electrons in AliVEvent			*
   *	track1:		tagged electron							*
   *	vEvent:		event								*
   *	weight:		weight in pt if not realistic					*
   *	binct:		centrality bin							*
   *	deltaphi:	phi-phi event plane for v2					*
   *	source:		MC sources							*
   *	indexmother:	MC index mother							*
   *											*
   *											*
   *	return -1  if  nothing								*
   *	return  2  if  opposite		charge		within the mass range		*
   *	return  4  if      like		charge		within the mass range		*
   *	return  6  if  opposite & like charge		within the mass range		*
   *											*
   ***********************************************************************************/

  AliAODEvent *aodeventu = dynamic_cast<AliAODEvent*>(vEvent);
  Int_t taggedphotonic = -1;

  AliDebug(2,Form("fCounterPoolBackground %d in LookAtNonHFE!!!",fCounterPoolBackground));
  if(!fArraytrack) return taggedphotonic;
  AliDebug(2,Form("process track %d",iTrack1));

  fkPIDRespons = fPIDBackground->GetPIDResponse();

  //Set Fill-Arrays for THnSparse
  Double_t valueIncElectron[7]	= { binct, track1->Pt(), source, track1->P(), track1->GetTPCsignal(), fkPIDRespons->NumberOfSigmasTPC( track1, AliPID::kElectron), fkPIDRespons->NumberOfSigmasTOF( track1, AliPID::kElectron)};	//Centrality	Pt	Source	P	TPCsignal	TPCsigma	TOFsigma
  Double_t valueSign[7]		= { deltaphi, binct, track1->Pt(), -1, source, -1, -1};			//DeltaPhi	Centrality	Pt	InvariantMass	Source	Angle	Pt
  //Double_t valueAngle[3]	= { -1, binct, source};								//Angle		Centrality	Source

  Int_t pdg1 = CheckPdg(TMath::Abs(track1->GetLabel()));
  Double_t eMass = TDatabasePDG::Instance()->GetParticle(11)->Mass(); //Electron mass in GeV
  Double_t bfield = vEvent->GetMagneticField();

  AliVTrack *track2 = 0x0;
  Int_t iTrack2 = 0;
  Int_t indexmother2 = -1;
  Int_t pdg2 = -100;
  Int_t source2 = -1;
  Int_t fPDGtrack2 = 0;
  Float_t fCharge2 = 0;

  Double_t dca12 = 0;

  TLorentzVector electron1;
  TLorentzVector electron2;
  TLorentzVector mother;

  Double_t xt1 = 0; //radial position track 1 at the DCA point
  Double_t xt2 = 0; //radial position track 2 at the DCA point
  Double_t p1[3] = {0,0,0};
  Double_t p2[3] = {0,0,0};
  Double_t angleESD = -1;
  Double_t invmassESD = -1;

  Double_t chi2OverNDF = -1;
  Double_t width = 0;
  Double_t angleAOD = -1;
  Double_t invmassAOD = -1;

  AliKFVertex primV(*(vEvent->GetPrimaryVertex()));

  Float_t fCharge1 = track1->Charge();							//Charge from track1
  Int_t fPDGtrack1 = 11;
  if(fCharge1>0) fPDGtrack1 = -11;
  AliKFParticle ktrack1(*track1, fPDGtrack1);
  AliESDtrack *esdtrack1 = dynamic_cast<AliESDtrack *>(track1);			//ESD-track1

  AliESDtrack *esdtrack2 = 0x0;

  Bool_t kUSignPhotonic = kFALSE;
  Bool_t kLSignPhotonic = kFALSE;
  Bool_t kHasdcaT1 = kFALSE;
  Bool_t kHasdcaT2 = kFALSE;

  //! FILL Inclusive Electron
  fIncElectron->Fill(valueIncElectron,weight);
  fnumberfound++;
  //printf(Form("Inclusive Pool: TrackNr. %d, fnumberfound %d \n", iTrack1, fnumberfound));

  for(Int_t idex = 0; idex < fCounterPoolBackground; idex++)
  {
    iTrack2 = fArraytrack->At(idex);
    AliDebug(2,Form("track %d",iTrack2));
    track2 = (AliVTrack *)vEvent->GetTrack(iTrack2);

    if(!track2)
    {
      //printf("ERROR: Could not receive track %d", iTrack2);
      continue;
    }

    fCharge2 = track2->Charge();		//Charge from track2

    // Reset the MC info
    //valueAngle[2] = source;
    valueSign[4] = source;
    valueSign[6] = track2->Pt();

    // track cuts and PID already done

    // Checking if it is the same Track!
    if(iTrack2==iTrack1) continue;
    AliDebug(2,"Different");

    // if MC look
    if(fMCEvent || fAODArrayMCInfo)
    {
      source2	 = FindMother(TMath::Abs(track2->GetLabel()), indexmother2);
      pdg2	 = CheckPdg(TMath::Abs(track2->GetLabel()));

      if(source2 >=0 )
      {
	if((indexmother2 == indexmother) && (source == source2) && ((pdg1*pdg2)<0.0))
	{
	  if(source == kElectronfromconversion)
	  {
	    //valueAngle[2] = kElectronfromconversionboth;
	    valueSign[4] = kElectronfromconversionboth;
	  }

	  if(source == kElectronfrompi0)
	  {
	    //valueAngle[2] = kElectronfrompi0both;
	    valueSign[4] = kElectronfrompi0both;
	  }

	  if(source == kElectronfrometa)
	  {
	    //valueAngle[2] = kElectronfrometaboth;
	    valueSign[4] = kElectronfrometaboth;
	  }
	}
      }
    }


    if(fAlgorithmMA && (!aodeventu)) 
    {
      /**				*
       *	ESD-Analysis		*
       **				*/

      esdtrack2 = dynamic_cast<AliESDtrack *>(track2);			//ESD-track2
      if((!esdtrack1) || (!esdtrack2)) continue;

      dca12 = esdtrack2->GetDCA(esdtrack1,bfield,xt2,xt1);		//DCA track1-track2

      if(dca12 > fMaxDCA) continue;					//! Cut on DCA

      //Momento of the track extrapolated to DCA track-track
      kHasdcaT1 = esdtrack1->GetPxPyPzAt(xt1,bfield,p1);		//Track1
      kHasdcaT2 = esdtrack2->GetPxPyPzAt(xt2,bfield,p2);		//Track2
      if(!kHasdcaT1 || !kHasdcaT2) AliWarning("It could be a problem in the extrapolation");

      electron1.SetXYZM(p1[0], p1[1], p1[2], eMass);
      electron2.SetXYZM(p2[0], p2[1], p2[2], eMass);

//      electron1.SetXYZM(esdtrack1->Px(), esdtrack1->Py(), esdtrack1->Pz(), eMass);
//      electron2.SetXYZM(esdtrack2->Px(), esdtrack2->Py(), esdtrack2->Pz(), eMass);

      mother      = electron1 + electron2;
      invmassESD  = mother.M();
      angleESD    = TVector2::Phi_0_2pi(electron1.Angle(electron2.Vect()));

      //valueAngle[0] = angleESD;
      valueSign[3] = invmassESD;
      valueSign[5] = angleESD;

      //if((fCharge1*fCharge2)>0.0)	fLSignAngle->Fill(&valueAngle[0],weight);
      //else				fUSignAngle->Fill(&valueAngle[0],weight);

      if(angleESD > fMaxOpening3D) continue;				 //! Cut on Opening Angle
      if(invmassESD > fMaxInvMass) continue;				//! Cut on Invariant Mass

      if((fCharge1*fCharge2)>0.0)	fLSign->Fill( valueSign, weight);
      else				fUSign->Fill( valueSign, weight);

      if((fCharge1*fCharge2)>0.0)	kLSignPhotonic=kTRUE;
      else				kUSignPhotonic=kTRUE;
    }
    else
    {
      /**				*
       *	AOD-AliKF-Analysis	*
       **				*/

      //printf("AOD HFE non photonic\n");

      fPDGtrack2 = 11;
      if(fCharge2>0) fPDGtrack2 = -11;

      AliKFParticle ktrack2(*track2, fPDGtrack2);
      AliKFParticle recoGamma(ktrack1,ktrack2);

      if(recoGamma.GetNDF()<1) continue;				//! Cut on Reconstruction

      chi2OverNDF = recoGamma.GetChi2()/recoGamma.GetNDF();
      if(TMath::Sqrt(TMath::Abs(chi2OverNDF))>fChi2OverNDFCut) continue;

      // DCA
      //Double_t dca12 = ktrack1.GetDistanceFromParticle(ktrack2);
      //if(dca12 > fMaxDCA) continue;

      // if set mass constraint
      if(fSetMassConstraint) //&& pVtx)
      {
	primV += recoGamma;
	primV -= ktrack1;
	primV -= ktrack2;
	recoGamma.SetProductionVertex(primV);
	recoGamma.SetMassConstraint(0,0.0001);
      }

      recoGamma.GetMass(invmassAOD,width);
      angleAOD = ktrack1.GetAngle(ktrack2);

      //valueAngle[0] = angleAOD;
      valueSign[3] = invmassAOD;
      valueSign[5] = angleAOD;

      //if((fCharge1*fCharge2)>0.0)	fLSignAngle->Fill(&valueAngle[0],weight);
      //else				fUSignAngle->Fill(&valueAngle[0],weight);

      if(angleAOD > fMaxOpening3D) continue;				//! Cut on Opening Angle
      if(invmassAOD > fMaxInvMass) continue;				//! Cut on Invariant Mass

      if((fCharge1*fCharge2)>0.0)	fLSign->Fill( valueSign, weight);
      else				fUSign->Fill( valueSign, weight);

      if((fCharge1*fCharge2)>0.0)	kLSignPhotonic=kTRUE;
      else				kUSignPhotonic=kTRUE;
    }
  }

  if( kUSignPhotonic &&  kLSignPhotonic) taggedphotonic = 6;
  if(!kUSignPhotonic &&  kLSignPhotonic) taggedphotonic = 4;
  if( kUSignPhotonic && !kLSignPhotonic) taggedphotonic = 2;

  return taggedphotonic;
}

//_________________________________________________________________________
Int_t AliHFENonPhotonicElectron::FindMother(Int_t tr, Int_t &indexmother){
  //
  // Find the mother if MC
  //

  if(!fMCEvent && !fAODArrayMCInfo) return 0;

  Int_t pdg = CheckPdg(tr);
  if(TMath::Abs(pdg)!= 11)
  {
    indexmother = -1;
    return kNoElectron;
  }

  indexmother = IsMotherGamma(tr);
  if(indexmother > 0) return kElectronfromconversion;
  indexmother = IsMotherPi0(tr);
  if(indexmother > 0) return kElectronfrompi0;
  indexmother = IsMotherC(tr);
  if(indexmother > 0) return kElectronfromC;
  indexmother = IsMotherB(tr);
  if(indexmother > 0) return kElectronfromB;
  indexmother = IsMotherEta(tr);
  if(indexmother > 0) return kElectronfrometa;

  return kElectronfromother;
}

//________________________________________________________________________________________________
Int_t AliHFENonPhotonicElectron::CheckPdg(Int_t tr) {

  //
  // Return the pdg of the particle
  //

  Int_t pdgcode = -1;
  if(tr < 0) return pdgcode;

  if(fMCEvent)
  {
    AliVParticle *mctrack = fMCEvent->GetTrack(tr);
    if(!mctrack) return -1;
    AliMCParticle *mctrackesd = NULL;
    if(!(mctrackesd = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(tr))))) return pdgcode;
    pdgcode = mctrackesd->PdgCode();
  }

  if(fAODArrayMCInfo)
  {
    if((tr+1)>fAODArrayMCInfo->GetEntriesFast()) return -1;
    AliAODMCParticle *mctrackaod = (AliAODMCParticle *) fAODArrayMCInfo->At(tr);
    if(!mctrackaod) return pdgcode;
    pdgcode = mctrackaod->GetPdgCode();
  }

  return pdgcode;
}

//_______________________________________________________________________________________________
Int_t AliHFENonPhotonicElectron::IsMotherGamma(Int_t tr) {

  //
  // Return the lab of gamma mother or -1 if not gamma
  //

  if(tr < 0) return -1;

  if(fMCEvent)
  {
    AliVParticle *mctrack = fMCEvent->GetTrack(tr);
    if(!mctrack) return -1;
    AliMCParticle *mctrackesd = NULL;
    if(!(mctrackesd = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(tr))))) return -1;
    TParticle *particle = 0x0;
    particle = mctrackesd->Particle();

    // Take mother
    if(!particle) return -1;
    Int_t imother   = particle->GetFirstMother();
    if(imother < 0) return -1;
    AliMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(imother))))) return -1;
    TParticle * mother = mothertrack->Particle();
    if(!mother) return -1;

    // Check gamma
    Int_t pdg = mother->GetPdgCode();
    if(TMath::Abs(pdg) == 22) return imother;
    if(TMath::Abs(pdg) == 11)
    {
      return IsMotherGamma(imother);
    }

    return -1;
  }

  if(fAODArrayMCInfo)
  {
    if((tr+1)>fAODArrayMCInfo->GetEntriesFast()) return -1;
    AliAODMCParticle *mctrackaod = (AliAODMCParticle *) fAODArrayMCInfo->At(tr);
    if(!mctrackaod) return -1;

    // Take mother
    Int_t imother = mctrackaod->GetMother();
    if(imother < 0 || ((imother+1)>fAODArrayMCInfo->GetEntriesFast())) return -1;
    AliAODMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(imother))))) return -1;

    // Check gamma
    Int_t pdg = mothertrack->GetPdgCode();
    if(TMath::Abs(pdg) == 22) return imother;
    if(TMath::Abs(pdg) == 11)
    {
      return IsMotherGamma(imother);
    }

    return -1;
  }

  return -1;
}

//________________________________________________________________________________________________
Int_t AliHFENonPhotonicElectron::IsMotherPi0(Int_t tr) {

  //
  // Return the lab of pi0 mother or -1 if not pi0
  //

  if(tr < 0) return -1;

  if(fMCEvent)
  {
    AliVParticle *mctrack = fMCEvent->GetTrack(tr);
    if(!mctrack) return -1;
    AliMCParticle *mctrackesd = NULL;
    if(!(mctrackesd = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(tr))))) return -1;
    TParticle *particle = 0x0;
    particle = mctrackesd->Particle();
    // Take mother
    if(!particle) return -1;
    Int_t imother   = particle->GetFirstMother();
    if(imother < 0) return -1;
    AliMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(imother))))) return -1;
    TParticle * mother = mothertrack->Particle();
    if(!mother) return -1;
    // Check pi0
    Int_t pdg = mother->GetPdgCode();
    if(TMath::Abs(pdg) == 111) return imother;
    if(TMath::Abs(pdg) == 11)
    {
      return IsMotherPi0(imother);
    }

    return -1;
  }

  if(fAODArrayMCInfo)  {

    if((tr+1)>fAODArrayMCInfo->GetEntriesFast()) return -1;
    AliAODMCParticle *mctrackaod = (AliAODMCParticle *) fAODArrayMCInfo->At(tr);
    if(!mctrackaod) return -1;

    // Take mother
    Int_t imother = mctrackaod->GetMother();
    if(imother < 0 || ((imother+1)>fAODArrayMCInfo->GetEntriesFast())) return -1;
    AliAODMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(imother))))) return -1;
    // Check pi0
    Int_t pdg = mothertrack->GetPdgCode();
    if(TMath::Abs(pdg) == 111) return imother;
    if(TMath::Abs(pdg) == 11)
    {
      return IsMotherPi0(imother);
    }

    return -1;
  }

  return -1;
}
//________________________________________________________________________________________________
Int_t AliHFENonPhotonicElectron::IsMotherC(Int_t tr) {

  //
  // Return the lab of signal mother or -1 if not from C
  //

  if(tr < 0) return -1;

  if(fMCEvent)
  {
    AliVParticle *mctrack = fMCEvent->GetTrack(tr);
    if(!mctrack) return -1;
    AliMCParticle *mctrackesd = NULL;
    if(!(mctrackesd = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(tr))))) return -1;
    TParticle *particle = 0x0;
    particle = mctrackesd->Particle();

    // Take mother
    if(!particle) return -1;
    Int_t imother   = particle->GetFirstMother();
    if(imother < 0) return -1;
    AliMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(imother))))) return -1;
    TParticle * mother = mothertrack->Particle();
    if(!mother) return -1;

    // Check C
    Int_t pdg = mother->GetPdgCode();
    if((TMath::Abs(pdg)==411) || (TMath::Abs(pdg)==421) || (TMath::Abs(pdg)==431) || (TMath::Abs(pdg)==4122) || (TMath::Abs(pdg)==4132) || (TMath::Abs(pdg)==4232) || (TMath::Abs(pdg)==43320)) return imother;
    if(TMath::Abs(pdg) == 11)
    {
      return IsMotherC(imother);
    }

    return -1;
  }

  if(fAODArrayMCInfo)
  {
    if((tr+1)>fAODArrayMCInfo->GetEntriesFast()) return -1;
    AliAODMCParticle *mctrackaod = (AliAODMCParticle *) fAODArrayMCInfo->At(tr);
    if(!mctrackaod) return -1;

    // Take mother
    Int_t imother = mctrackaod->GetMother();
    if(imother < 0 || ((imother+1)>fAODArrayMCInfo->GetEntriesFast())) return -1;
    AliAODMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(imother))))) return -1;

    // Check C
    Int_t pdg = mothertrack->GetPdgCode();
    if((TMath::Abs(pdg)==411) || (TMath::Abs(pdg)==421) || (TMath::Abs(pdg)==431) || (TMath::Abs(pdg)==4122) || (TMath::Abs(pdg)==4132) || (TMath::Abs(pdg)==4232) || (TMath::Abs(pdg)==43320)) return imother;
    if(TMath::Abs(pdg) == 11)
    {
      return IsMotherC(imother);
    }

    return -1;
  }

  return -1;
}
//_______________________________________________________________________________________________
Int_t AliHFENonPhotonicElectron::IsMotherB(Int_t tr) {

  //
  // Return the lab of signal mother or -1 if not B
  //

  if(tr < 0) return -1;

  if(fMCEvent)
  {
    AliVParticle *mctrack = fMCEvent->GetTrack(tr);
    if(!mctrack) return -1;
    AliMCParticle *mctrackesd = NULL;
    if(!(mctrackesd = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(tr))))) return -1;
    TParticle *particle = 0x0;
    particle = mctrackesd->Particle();

    // Take mother
    if(!particle) return -1;
    Int_t imother   = particle->GetFirstMother();
    if(imother < 0) return -1;
    AliMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(imother))))) return -1;
    TParticle * mother = mothertrack->Particle();
    if(!mother) return -1;

    // Check B
    Int_t pdg = mother->GetPdgCode();
    if((TMath::Abs(pdg)==511) || (TMath::Abs(pdg)==521) || (TMath::Abs(pdg)==531) || (TMath::Abs(pdg)==5122) || (TMath::Abs(pdg)==5132) || (TMath::Abs(pdg)==5232) || (TMath::Abs(pdg)==53320)) return imother;
    if(TMath::Abs(pdg) == 11)
    {
      return IsMotherB(imother);
    }

    return -1;
  }

  if(fAODArrayMCInfo)
  {
    if((tr+1)>fAODArrayMCInfo->GetEntriesFast()) return -1;
    AliAODMCParticle *mctrackaod = (AliAODMCParticle *) fAODArrayMCInfo->At(tr);
    if(!mctrackaod) return -1;

    // Take mother
    Int_t imother = mctrackaod->GetMother();
    if(imother < 0 || ((imother+1)>fAODArrayMCInfo->GetEntriesFast())) return -1;
    AliAODMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(imother))))) return -1;
    // Check B
    Int_t pdg = mothertrack->GetPdgCode();
    if((TMath::Abs(pdg)==511) || (TMath::Abs(pdg)==521) || (TMath::Abs(pdg)==531) || (TMath::Abs(pdg)==5122) || (TMath::Abs(pdg)==5132) || (TMath::Abs(pdg)==5232) || (TMath::Abs(pdg)==53320)) return imother;
    if(TMath::Abs(pdg) == 11)
    {
      return IsMotherB(imother);
    }

    return -1;
  }

  return -1;
}

//_______________________________________________________________________________________________
Int_t AliHFENonPhotonicElectron::IsMotherEta(Int_t tr) {

  //
  // Return the lab of eta mother or -1 if not eta
  //

  if(tr < 0) return -1;

  if(fMCEvent)
  {
    AliVParticle *mctrack = fMCEvent->GetTrack(tr);
    if(!mctrack) return -1;
    AliMCParticle *mctrackesd = NULL;
    if(!(mctrackesd = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(tr))))) return -1;
    TParticle *particle = 0x0;
    particle = mctrackesd->Particle();

    // Take mother
    if(!particle) return -1;
    Int_t imother   = particle->GetFirstMother();
    if(imother < 0) return -1;
    AliMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(imother))))) return -1;
    TParticle * mother = mothertrack->Particle();
    if(!mother) return -1;

    // Check eta
    Int_t pdg = mother->GetPdgCode();
    if(TMath::Abs(pdg) == 221) return imother;
    if(TMath::Abs(pdg) == 11)
    {
      return IsMotherEta(imother);
    }

    return -1;
  }

  if(fAODArrayMCInfo)
  {
    if((tr+1)>fAODArrayMCInfo->GetEntriesFast()) return -1;
    AliAODMCParticle *mctrackaod = (AliAODMCParticle *) fAODArrayMCInfo->At(tr);
    if(!mctrackaod) return -1;

    // Take mother
    Int_t imother = mctrackaod->GetMother();
    if(imother < 0 || ((imother+1)>fAODArrayMCInfo->GetEntriesFast())) return -1;
    AliAODMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(imother))))) return -1;

    // Check eta
    Int_t pdg = mothertrack->GetPdgCode();
    if(TMath::Abs(pdg) == 221) return imother;
    if(TMath::Abs(pdg) == 11)
    {
      return IsMotherEta(imother);
    }

    return -1;
  }

  return -1;
}
