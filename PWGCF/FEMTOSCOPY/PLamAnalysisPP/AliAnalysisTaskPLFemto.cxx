/*
 **************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appeuear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */
// Class to calculate P-V0 (especially P-Lambda) correlation functions
// Author: O. Arnold
// inherited a lot from other analyses, especially from H. Beck (thanks!)
//
//-----------------------------------------------------------------------

#include <TSystem.h>
#include <TParticle.h>
#include "TROOT.h"
#include <TDatabasePDG.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecay.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskPLFemto.h"
#include "AliInputEventHandler.h"
#include "AliFemtoLambdaEventCollection2.h"
#include "AliFemtoLambdaEvent.h"
#include "AliFemtoLambdaParticle.h"
#include "AliFemtoProtonParticle.h"
#include "AliPPVsMultUtils.h"
#include "AliAODTracklets.h"
#include "AliMultSelection.h"

#include "TVector2.h"

ClassImp(AliAnalysisTaskPLFemto)

Float_t TPCradii[9] = {85.,105.,125.,145.,165.,185.,205.,225.,245.};//must be global to be used on grid
//__________________________________________________________________________
AliAnalysisTaskPLFemto::AliAnalysisTaskPLFemto():
  AliAnalysisTaskSE(),
  fAliEventCuts(),
  fTrigger(AliVEvent::kMB),
  fIsRun1(true),
  fIsLightweight(false),
  fV0PercentileMax(100),
  fUseMCInfo(kFALSE),
  fOnlineV0(kFALSE),
  fOutput(0),
  fOutputSP(0),
  fOutputTP(0),
  fOutputPID(0),
  fAODMCEvent(0),
  fCEvents(0),
  fwhichV0("Lambda"),
  fwhichAntiV0("Anti-"+fwhichV0),
  fwhichV0region("signal"),
  fV0Counter(0),
  fAntiV0Counter(0),
  fProtonCounter(0),
  fAntiProtonCounter(0),
  fXiCounter(0),
  fEventNumber(0),
  fTPCradii(0),
  fcutType(AliFemtoCutValues::kDefault),
  fPIDResponse(0),
  fGTI(0),
  fTrackBuffSize(2500),
  fV0cand(new AliFemtoLambdaParticle[kV0TrackLimit]),
  fAntiV0cand(new AliFemtoLambdaParticle[kV0TrackLimit]),
  fProtoncand(new AliFemtoProtonParticle[kProtonTrackLimit]),
  fAntiProtoncand(new AliFemtoProtonParticle[kProtonTrackLimit]),
  fXicand(new AliFemtoXiParticle[kXiTrackLimit]),
  fCuts(new AliFemtoCutValues(AliFemtoCutValues::kDefault)),
  fAnaUtils(new AliAnalysisUtils()),
  fProtonTrackVector(),
  fAntiProtonTrackVector(),
  fLambdaTrackVector(),
  fAntiLambdaTrackVector(),
  fXiTrackVector(),
  fProtonEvtBuffer(),
  fAntiProtonEvtBuffer(),
  fLambdaEvtBuffer(),
  fAntiLambdaEvtBuffer(),
  fXiEvtBuffer()
{
  Info("AliAnalysisTaskPLFemto","Calling default Constructor");
  //
  // Constructor. Initialization of Inputs and Outputs
  //

  fTPCradii = TPCradii;

  for(Int_t i=0; i<3;i++) fPrimVertex[i] = -9999.;

  fCuts->SetCutVariation(fcutType);
  fWhichfilterbit = fCuts->GetFilterBit();
  fAnaUtils->SetMinPlpContribSPD(3);
}
//___________________________________________________________________________
AliAnalysisTaskPLFemto::AliAnalysisTaskPLFemto(const Char_t* name,Bool_t OnlineCase,TString whichV0,TString whichV0region) :
  AliAnalysisTaskSE(name),
  fAliEventCuts(),
  fTrigger(AliVEvent::kMB),
  fIsRun1(true),
  fIsLightweight(false),
  fV0PercentileMax(100),
  fUseMCInfo(kFALSE),
  fOnlineV0(OnlineCase),
  fOutput(0),
  fOutputSP(0),
  fOutputTP(0),
  fOutputPID(0),
  fAODMCEvent(0),
  fCEvents(0),
  fwhichV0(whichV0),
  fwhichAntiV0("Anti-"+whichV0),
  fwhichV0region(whichV0region),
  fV0Counter(0),
  fAntiV0Counter(0),
  fProtonCounter(0),
  fAntiProtonCounter(0),
  fXiCounter(0),
  fEventNumber(0),
  fTPCradii(0),
  fcutType(AliFemtoCutValues::kDefault),
  fPIDResponse(0),
  fGTI(0),
  fTrackBuffSize(2500),
  fV0cand(new AliFemtoLambdaParticle[kV0TrackLimit]),
  fAntiV0cand(new AliFemtoLambdaParticle[kV0TrackLimit]),
  fProtoncand(new AliFemtoProtonParticle[kProtonTrackLimit]),
  fAntiProtoncand(new AliFemtoProtonParticle[kProtonTrackLimit]),
  fXicand(new AliFemtoXiParticle[kXiTrackLimit]),
  fCuts(new AliFemtoCutValues(AliFemtoCutValues::kDefault)),
  fAnaUtils(new AliAnalysisUtils()),
  fProtonTrackVector(),
  fAntiProtonTrackVector(),
  fLambdaTrackVector(),
  fAntiLambdaTrackVector(),
  fXiTrackVector(),
  fProtonEvtBuffer(),
  fAntiProtonEvtBuffer(),
  fLambdaEvtBuffer(),
  fAntiLambdaEvtBuffer(),
  fXiEvtBuffer()
{
  Info("AliAnalysisTaskPLFemto","Calling extended Constructor");

  fTPCradii = TPCradii;

  for(Int_t i=0; i<3;i++) fPrimVertex[i] = -9999.;

  DefineOutput(1,TList::Class());  //conters
  DefineOutput(2,TList::Class());
  DefineOutput(3,TList::Class());
  DefineOutput(4,TList::Class());
  DefineOutput(5,TList::Class());

  fAnaUtils->SetMinPlpContribSPD(3);
  fCuts->SetCutVariation(fcutType);
  fWhichfilterbit = fCuts->GetFilterBit();
}

//___________________________________________________________________________
AliAnalysisTaskPLFemto::~AliAnalysisTaskPLFemto() {
  //
  // destructor
  //
  Info("~AliAnalysisTaskPLFemto","Calling Destructor");

  delete fOutput;
  delete fOutputSP;
  delete fOutputTP;
  delete fOutputPID;
  delete fCEvents;

  if (fGTI) { delete[] fGTI; fGTI=NULL;}

  if(fV0cand) delete[] fV0cand;
  if(fAntiV0cand) delete[] fAntiV0cand;
  if(fProtoncand) delete[] fProtoncand;
  if(fAntiProtoncand) delete[] fAntiProtoncand;
  if(fXicand) delete[] fXicand;
  if(fCuts){ delete fCuts; fCuts = NULL;}
  if(fPIDResponse){ delete fPIDResponse; fPIDResponse = NULL;}

  // TODO: Uncomment this line when the AliAnalysisUtils destructor is fixed
  // if(fAnaUtils) delete fAnaUtils;
}
//_________________________________________________
void AliAnalysisTaskPLFemto::Init()
{
  //
  // Initialization
  //
  return;
}
Bool_t AliAnalysisTaskPLFemto::GoodTPCFitMapSharedMap(const AliAODTrack *pTrack,const AliAODTrack *nTrack)
{
  //This method was inherited form H. Beck analysis

  // Rejects tracks with shared clusters
  // This overload is used for positive and negative daughters from V0s

  // Get the shared maps

  const TBits posSharedMap = pTrack->GetTPCSharedMap();
  const TBits negSharedMap = nTrack->GetTPCSharedMap();

  if( ((posSharedMap.CountBits()) >= 1) || ((negSharedMap.CountBits()) >= 1))
    {
      // Bad tracks, have too many shared clusters!
      return kFALSE;
    }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPLFemto::GoodTPCFitMapSharedMap(const AliAODTrack *track)
{
  //This method was inherited form H. Beck analysis

  // Rejects tracks with shared clusters
  // This overload is used for primaries

  // Get the shared maps
  const TBits sharedMap = track->GetTPCSharedMap();
  if((sharedMap.CountBits()) >= 1)
    {
      // Bad track, has too many shared clusters!
      return kFALSE;
    }
  return kTRUE;
}
//________________________________________________________________________
Double_t AliAnalysisTaskPLFemto::LpCorrfunc(Double_t relk,Double_t source)
{
  //Dummy Lambda-p correlation function for a simple MC study
  const Double_t hbarc = 0.197;

  Double_t qinv = 2. * relk;

  return 1. + TMath::Exp(-pow(qinv*source,2.)/pow(hbarc,2.));
}
//________________________________________________________________________
Double_t AliAnalysisTaskPLFemto::ppCorrfunc(Double_t relk,Double_t source)
{
  //Dummy p-p correlation function for a simple MC study
  const Double_t hbarc = 0.197;

  Double_t qinv = 2. * relk;

  Double_t rectang = 0.;
  if(relk > 0.4 && relk < 0.6) rectang = 2.; //artificial rectangular weighting

  return 1. + TMath::Exp(-qinv*source/hbarc) + rectang; // in this momentum region the first part is zero and the rectang part survives
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPLFemto::AcceptTrack(const AliAODTrack *track,int type)
{
  //This method was inherited from H. Beck analysis

  // Apply additional track cuts

  // In the documents
  // https://alisoft.cern.ch/AliRoot/trunk/TPC/doc/Definitions/Definitions.pdf
  // TPC people describe the cut strategy for the TPC. It is explicitly
  // stated that a cut on the number of crossed rows and a cut on the
  // number of crossed rows over findable clusters is recommended to
  // remove fakes. In the pdf a cut value of .83 on the ratio
  // is stated, no value for the number of crossed rows. Looking at the
  // AliESDtrackCuts.cxx one sees that exactly this cut is used with
  // 0.8 on the ratio and 70 on the crossed rows.

  // Checked the filter task and AliAODTrack and AliESDtrack and
  // AliESDtrackCuts and the Definitions.pdf:
  // The function to get the findable clusters is GetTPCNclsF()

  // For the number fo crossed rows for ESD tracks, the function
  // GetTPCCrossedRows() usually is used. Looking at the AliESDtrack.cxx
  // one sees that it's just an alias (with additional caching) for
  // GetTPCClusterInfo(2, 1); The identical function exists in the
  // AliAODTrack.cxx

  TString fillthis = "";

  Float_t nCrossed = track->GetTPCClusterInfo(2,1);
  UShort_t nFindableCluster = track->GetTPCNclsF();
  Float_t ratio = nCrossed / Float_t(nFindableCluster);



  if(type == 4)
    {
      fillthis = "fFindableClusterProton";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(nFindableCluster);

      fillthis = "fNCrossedRowsProton";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(nCrossed);

      fillthis = "fRatioFindableCrossedProton";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(ratio);
    }
  else if(type == -4)
    {
      fillthis = "fFindableClusterAntiProton";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(nFindableCluster);

      fillthis = "fNCrossedRowsAntiProton";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(nCrossed);

      fillthis = "fRatioFindableCrossedAntiProton";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(ratio);
    }

  if(nCrossed<70)
    return kFALSE;
  if(!track->GetTPCNclsF())
    return kFALSE; // Note that the AliESDtrackCuts would here return kTRUE
  if((nCrossed/track->GetTPCNclsF()) < .83)
    return kFALSE;
  return kTRUE;
}
//_________________________________________________
void AliAnalysisTaskPLFemto::ResetGlobalTrackReference()
{
  //This method was inherited form H. Beck analysis

  // Sets all the pointers to zero. To be called at
  // the beginning or end of an event
  for(UShort_t i=0;i<fTrackBuffSize;i++)
    {
      fGTI[i]=0;
    }
}

//_________________________________________________
void AliAnalysisTaskPLFemto::StoreGlobalTrackReference(AliAODTrack *track)
{
  //This method was inherited form H. Beck analysis

  // Stores the pointer to the global track

  // This was AOD073
  // // Don't use the filter bits 2 (ITS standalone) and 128 TPC only
  // // Remove this return statement and you'll see they don't have
  // // any TPC signal
  // if(track->TestFilterBit(128) || track->TestFilterBit(2))
  //   return;
  // This is AOD086
  // Another set of tracks was introduced: Global constrained.
  // We only want filter bit 1 <-- NO! we also want no
  // filter bit at all, which are the v0 tracks
  //  if(!track->TestFilterBit(1))
  //    return;

  // There are also tracks without any filter bit, i.e. filter map 0,
  // at the beginning of the event: they have ~id 1 to 5, 1 to 12
  // This are tracks that didn't survive the primary track filter but
  // got written cause they are V0 daughters

  // Check whether the track has some info
  // I don't know: there are tracks with filter bit 0
  // and no TPC signal. ITS standalone V0 daughters?
  // if(!track->GetTPCsignal()){
  //   printf("Warning: track has no TPC signal, "
  // 	   //	   "not adding it's info! "
  // 	   "ID: %d FilterMap: %d\n"
  // 	   ,track->GetID(),track->GetFilterMap());
  //   //    return;
  // }

  // Check that the id is positive (global track)
  if(track->GetID()<0)
    {
      //    printf("Warning: track has negative ID: %d\n",track->GetID());
      return;
    }

  // Check id is not too big for buffer
  if(track->GetID()>=fTrackBuffSize)
    {
      printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n"
	     ,track->GetID(),fTrackBuffSize);
      return;
    }

  // Warn if we overwrite a track
  if(fGTI[track->GetID()])
    {
      // Seems like there are FilterMap 0 tracks
      // that have zero TPCNcls, don't store these!
      if( (!track->GetFilterMap()) &&
	  (!track->GetTPCNcls())   )
	return;

      if(!track->GetTPCNcls())
	{
	  std::cout << "TPC cluster" << track->GetTPCNcls() << std::endl;
	}

    // Imagine the other way around,
    // the zero map zero clusters track
    // is stored and the good one wants
    // to be added. We ommit the warning
    // and just overwrite the 'bad' track
      if( fGTI[track->GetID()]->GetFilterMap() ||
	  fGTI[track->GetID()]->GetTPCNcls()   )
	{
	  // If we come here, there's a problem
	  printf("Warning! global track info already there!");
	  printf("         TPCNcls track1 %u track2 %u",
		 (fGTI[track->GetID()])->GetTPCNcls(),track->GetTPCNcls());
	  printf("         FilterMap track1 %u track2 %u\n",
		 (fGTI[track->GetID()])->GetFilterMap(),track->GetFilterMap());
	}
    } // Two tracks same id

  // // There are tracks with filter bit 0,
  // // do they have TPCNcls stored?
  // if(!track->GetFilterMap()){
  //   printf("Filter map is zero, TPCNcls: %u\n"
  // 	   ,track->GetTPCNcls());
  // }

  //if(track->GetTPCNcls() <= 0) std::cout << "no tpc cluster" << std::endl;
  //if(fGTI[track->GetID()] >= 0) std::cout << "two tracks same ID " << std::endl;


  // Assign the pointer (ID of global track)
  (fGTI[track->GetID()]) = track;
}
//_________________________________________________________________
void AliAnalysisTaskPLFemto::DCAxy(const AliAODTrack *track, const AliVEvent *evt,Double_t *dcaxy)
{
  //standard "error" values:
  dcaxy[0] = -9999.;
  dcaxy[1] = -9999.;

  if(!track) return;

  // Create an external parameter from the AODtrack
  AliExternalTrackParam etp;
  etp.CopyFromVTrack(track);

  Double_t covar[3]={0.,0.,0.};
  if(!etp.PropagateToDCA(evt->GetPrimaryVertex(),evt->GetMagneticField(),10.,dcaxy,covar))
    {
      dcaxy[0] = -9999.;
      dcaxy[1] = -9999.;
      return;
    }
}
//_________________________________________________________________
Bool_t AliAnalysisTaskPLFemto::GoodDCA(const AliAODTrack *AODtrack,const AliAODEvent *aodEvent,Int_t type)
{
  TString fillthis = "";
  Double_t dcaxyz[2];

  if(fWhichfilterbit == 128) DCAxy(fGTI[-AODtrack->GetID()-1],fInputEvent,dcaxyz);
  else DCAxy(AODtrack,fInputEvent,dcaxyz);

  if(type == 4) //proton
    {
      fillthis = "fProtonDCAxyDCAz";
      if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaxyz[0],dcaxyz[1]);

      fillthis = "fProtonDCAxy";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaxyz[0]);

      if(fabs(dcaxyz[1])<fCuts->GetProtonDCAzCut())//check the influence of z cut on DCAxy
	{
	  fillthis = "fProtonDCAxyCutz";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaxyz[0]);
	}

      fillthis = "fProtonDCAz";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaxyz[1]);

      //Find the pt bin for pt dependent DCA analysis:
      int ptbin = fCuts->FindProtonPtBinDCA(AODtrack->Pt());
      if(!fUseMCInfo && ptbin != -9999)
	{
    if(!fIsLightweight) fillthis = "fProtonDCAxyDCAzPt";
	  fillthis += ptbin;
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaxyz[0],dcaxyz[1]);//total amount of protons in exp.
	}
      else if(fUseMCInfo && ptbin != -9999)
	{
	  TClonesArray *mcarray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	  if(!mcarray) return kFALSE;

	  AliAODMCParticle* mcproton = 0;
	  if(AODtrack->GetLabel() >=0) mcproton = (AliAODMCParticle*)mcarray->At(AODtrack->GetLabel());
	  else return kFALSE;

	  fillthis = "fProtonDCAxyDCAzMCCase";
	  fillthis += 0;
	  fillthis += "PtBin";
	  fillthis += ptbin;
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaxyz[0],dcaxyz[1]);//total amount of protons

	  //check the DCA of the proton in Monte Carlo
	  //if(mcproton->GetPdgCode() == 2212) //be shure it is a proton
	    {
	      if(mcproton->IsPhysicalPrimary() && !mcproton->IsSecondaryFromWeakDecay())
		{
      if(!fIsLightweight) fillthis = "fProtonDCAxyDCAzMCCase";
		  fillthis += 1;
		  fillthis += "PtBin";
		  fillthis += ptbin;
      if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaxyz[0],dcaxyz[1]);//primary protons
		}
	      else if(mcproton->IsSecondaryFromWeakDecay() && !mcproton->IsSecondaryFromMaterial())
		{
      if(!fIsLightweight) fillthis = "fProtonDCAxyDCAzMCCase";
		  fillthis += 2;
		  fillthis += "PtBin";
		  fillthis += ptbin;
      if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaxyz[0],dcaxyz[1]);//secondary protons


		  AliAODMCParticle* mcprotonMother = (AliAODMCParticle*)mcarray->At(mcproton->GetMother());
		  int PDGcode = mcprotonMother->GetPdgCode();

		  //check DCA of Lambda hyperon:
		  if(PDGcode == 3122)
		    {
          if(!fIsLightweight) fillthis = "fProtonDCAxyDCAzMCPtBinLambda";
		      fillthis += "PtBin";
		      fillthis += ptbin;
          if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaxyz[0],dcaxyz[1]);//secondary protons
		    }
		  if(PDGcode == 3222)
		    {
          if(!fIsLightweight) fillthis = "fProtonDCAxyDCAzMCPtBinSigma";
		      fillthis += "PtBin";
		      fillthis += ptbin;
          if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaxyz[0],dcaxyz[1]);//secondary protons
		    }
		}
	      else if(mcproton->IsSecondaryFromMaterial())
		{
      if(!fIsLightweight) fillthis = "fProtonDCAxyDCAzMCCase";
		  fillthis += 3;
		  fillthis += "PtBin";
		  fillthis += ptbin;
      if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaxyz[0],dcaxyz[1]);//material protons
		}
	    }
	    /*
	  else //background
	    {
	      fillthis = "fProtonDCAxyDCAzMCCase";
	      fillthis += 4;
	      fillthis += "PtBin";
	      fillthis += ptbin;
	      ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaxyz[0],dcaxyz[1]);//total amount of protons
	    }
	    */
	}
    }
  else if(type == -4) //antiproton
    {
      fillthis = "fAntiProtonDCAxyDCAz";
      if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaxyz[0],dcaxyz[1]);

      fillthis = "fAntiProtonDCAxy";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaxyz[0]);

      if(fabs(dcaxyz[1])<fCuts->GetProtonDCAzCut())//check the influence of z cut on DCAxy
	{
	  fillthis = "fAntiProtonDCAxyCutz";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaxyz[0]);
	}

      fillthis = "fAntiProtonDCAz";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaxyz[1]);
    }

  if(fabs(dcaxyz[0]) < fCuts->GetProtonDCAxyCut())//first cut on xy-value
    {
      if(fabs(dcaxyz[1]) < fCuts->GetProtonDCAzCut())//then on z-value
	{
	  return kTRUE;
	}
      else
	{
	  return kFALSE;
	}
    }
  else
    {
      return kFALSE;
    }
}
//________________________________________________________________________________________________
Int_t AliAnalysisTaskPLFemto::GetNumberOfTrackletsInEtaRange(AliAODEvent* ev,Double_t mineta,Double_t maxeta)
{
  // counts tracklets in given eta range
  AliAODTracklets* tracklets=ev->GetTracklets();
  Int_t nTr=tracklets->GetNumberOfTracklets();
  Int_t count=0;
  for(Int_t iTr=0; iTr<nTr; iTr++)
    {
      Double_t theta=tracklets->GetTheta(iTr);
      Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
      if(eta>mineta && eta<maxeta) count++;
    }
  TString fillthis = "fNTracklets";
  if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(count);
  /*
  AliAODHeader *header = dynamic_cast<AliAODHeader*>(ev->GetHeader());
  if(header)
    {
      count = header->GetRefMultiplicityComb08();
    }
  */
  return count;
}
//________________________________________________________________________________________________
void AliAnalysisTaskPLFemto::GetNumberOfTPCtracksInEtaRange(AliAODEvent* aodEvent,Double_t mineta,Double_t maxeta)
{

  Int_t count = 0;
  for(Int_t iTrack=0;iTrack<aodEvent->GetNumberOfTracks();iTrack++)
    {
      AliAODTrack *track = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTrack));
      if(!track)
	{
	  AliFatal("Not a standard AOD");
	  continue;
	}

      if(!track->TestFilterBit(128)) continue;

      Double_t eta=track->Eta();
      if(eta>mineta && eta<maxeta) count++;
    }
  TString fillthis = "fNTPCTracks";
  if(count>0 && !fIsLightweight)((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(count);

}
//________________________________________________________________________________________________
Bool_t AliAnalysisTaskPLFemto::PileUpRejection(AliAODEvent* ev)
{
  Bool_t isPileUp = kFALSE;

  Int_t fNofTracklets = ev->GetMultiplicity()->GetNumberOfTracklets();
  Int_t NofITSClusters1 = ev->GetNumberOfITSClusters(1);


  TString fillthis = "fEventSPDClusterNTracklets";
  ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(fNofTracklets,NofITSClusters1);

  if(!fAnaUtils->IsPileUpEvent(ev))
    {
      fillthis = "fEventSPDClusterNTrackletsPileUpCut";
      ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(fNofTracklets,NofITSClusters1);
    }


  return isPileUp;
}
//________________________________________________________________________________________________
Int_t *AliAnalysisTaskPLFemto::MultiplicityBinFinderzVertexBinFinder(AliAODEvent *aodEvent,Double_t zVertex)
{
  Int_t *MixingBins = new Int_t[2];

  Int_t multipl = GetNumberOfTrackletsInEtaRange(aodEvent,-0.8,0.8);//this is the multiplicty calculated with tracklets, (bug, before it was float)
  GetNumberOfTPCtracksInEtaRange(aodEvent,-0.8,0.8);

  Double_t RefMultiplicity08 = -1.f;
  AliAODHeader *header = dynamic_cast<AliAODHeader*>(aodEvent->GetHeader());
  if(header)
    {
      RefMultiplicity08 = header->GetRefMultiplicityComb08();
      TString fillthis = "fRefMultiplicity08";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(RefMultiplicity08);

      fillthis = "fNTrackletsRefMultiplicity";
      if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(multipl,RefMultiplicity08);
    }

  // Use RefMult08 for Run2
  multipl = (fIsRun1) ? multipl : RefMultiplicity08;

  if(multipl <= 0) //multiplicity estimation fails
    {
      MixingBins[0] = -9999;
      MixingBins[1] = -9999;
      return MixingBins;
    }

  Int_t multBin = 0;

  if(multipl <= 4) multBin = 0;
  else if(multipl <= 8) multBin = 1;
  else if(multipl <= 12) multBin = 2;
  else if(multipl <= 16) multBin = 3;
  else if(multipl <= 20) multBin = 4;
  else if(multipl <= 24) multBin = 5;
  else if(multipl <= 28) multBin = 6;
  else if(multipl <= 32) multBin = 7;
  else if(multipl <= 36) multBin = 8;
  else if(multipl <= 40) multBin = 9;
  else if(multipl <= 60) multBin = 10;
  else if(multipl <= 80) multBin = 11;
  else multBin = 12;

  if(multBin+1 > kMultiplicityBins) std::cout << "Change the Multiplicity maximum" << std::endl;

  Double_t deltaZ = (fCuts->GetEvtzVertexUp() - fCuts->GetEvtzVertexLow())/(Double_t)kZVertexBins;

  Int_t zVertexBin = 0;
  for(Int_t i=0;i<kZVertexBins;i++)
    {
      if((fCuts->GetEvtzVertexLow() + i*deltaZ)<zVertex && zVertex<(fCuts->GetEvtzVertexLow() + (i+1)*deltaZ))
	{
	  zVertexBin = i;
	  break;
	}
    }

  MixingBins[0] = zVertexBin;
  MixingBins[1] = multBin;

  //check filling of reference multiplicity:


  return MixingBins;
}
//_________________________________________________________________
Float_t AliAnalysisTaskPLFemto::CalcSphericityEvent(AliAODEvent *aodEvent)
{
  Float_t ptTot = 0.; //total Pt of all protons and v0s in the event

  Float_t s00 = 0.; //Elements of the sphericity matrix
  Float_t s11 = 0.;
  Float_t s10 = 0.;

  Int_t numOfTracks = aodEvent->GetNumberOfTracks();
  if(numOfTracks<3) return -9999.;//if already at this point not enough tracks are in the event -> return

  Int_t nTracks = 0;
  for(Int_t iTrack=0;iTrack<numOfTracks;iTrack++)
    {
      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTrack));

      if(!aodtrack->TestFilterBit(fWhichfilterbit)) continue;

      Float_t pt = aodtrack->Pt();
      //Double_t Phi = aodtrack->Phi();
      Float_t px = aodtrack->Px();
      Float_t py = aodtrack->Py();
      Float_t eta = aodtrack->Eta();


      if(!(eta>-0.8 && eta<0.8)) continue;
      if(pt<0.5) continue;

      ptTot += pt;

      s00 += px*px/pt;
      s11 += py*py/pt;
      s10 += px*py/pt;
      nTracks++;
    }

  if(nTracks<3) return -9999.;//new flag: check

  //normalize to total Pt to obtain a linear form:
  if(ptTot == 0.) return -9999.;
  s00 /= ptTot;
  s11 /= ptTot;
  s10 /= ptTot;

  //Calculate the trace of the sphericity matrix:
  Float_t T = s00+s11;

  //Calculate the determinant of the sphericity matrix:
  Float_t D = s00*s11 - s10*s10;//S10 = S01

  //Calculate the eigenvalues of the sphericity matrix:
  Float_t lambda1 = 0.5*(T + TMath::Sqrt(T*T - 4.*D));
  Float_t lambda2 = 0.5*(T - TMath::Sqrt(T*T - 4.*D));

  if((lambda1 + lambda2) == 0.) return -9999.;

  Float_t ST = -1.;

  if(lambda2>lambda1)
    {
      ST = 2.*lambda1/(lambda1+lambda2);
    }
  else
    {
      ST = 2.*lambda2/(lambda1+lambda2);
    }

  return ST;
}
//_________________________________________________________________
Bool_t AliAnalysisTaskPLFemto::CheckGlobalTrack(const AliAODTrack *track)
{
  //Checks if to the corresponding track a global track exists
  //This is especially useful if one has TPC only tracks and needs the PID information

  Bool_t isGlobal = kTRUE;

  if(!(fGTI[-track->GetID()-1])) isGlobal = kFALSE;

  return isGlobal;
}
//_________________________________________________________________
Float_t AliAnalysisTaskPLFemto::GetCorrectedTOFSignal(const AliVTrack *track,const AliAODEvent *aodevent)
{
  // Return the corrected TOF signal, see https://twiki.cern.ch/twiki/bin/viewauth/ALICE/TOF

  // Check for the global track
  /*
  if(fWhichfilterbit == 128 && !(fGTI[-track->GetID()-1]))
    {
      printf("Warning: no corresponding global track found!\n");
      return -9999.;
    }
  */


  // Request the TOFpid bit

  //AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,fGTI[-track->GetID()-1]);
  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,track);
  if(statusTOF != AliPIDResponse::kDetPidOk) return -9999.;

  // The expected time
  Float_t expectedTimes[AliPID::kSPECIES];
  //Double_t expectedTimes[AliPID::kProton];

  //(fGTI[-track->GetID()-1])->GetIntegratedTimes(expectedTimes);

  // Check for TOF header
  if(aodevent->GetTOFHeader())
    {
      // New AODs without start time subtraction
      //return ((fGTI[-track->GetID()-1])->GetTOFsignal() - expectedTimes[AliPID::kProton] - fPIDResponse->GetTOFResponse().GetStartTime(track->P()));
      return (track->GetTOFsignal() - expectedTimes[AliPID::kProton] - fPIDResponse->GetTOFResponse().GetStartTime(track->P()));
    }

  // Old AODs with start time already subtracted
  //return ((fGTI[-track->GetID()-1])->GetTOFsignal() - expectedTimes[AliPID::kProton]);
  return (track->GetTOFsignal() - expectedTimes[AliPID::kProton]);
}
//_________________________________________________________________
Bool_t AliAnalysisTaskPLFemto::SingleParticleQualityCuts(const AliAODTrack *aodtrack,const AliAODEvent *aodEvent)
{
  Double_t eta = aodtrack->Eta();
  Double_t pt = aodtrack->Pt();
  Double_t ch = aodtrack->Charge();

  if(ch>0.)//protons
    {
      if(fabs(eta) > fCuts->GetProtonEtaRange()) return kFALSE; //should lie within the ALICE acceptance
      if(pt < fCuts->GetProtonPIDthrPtLow()) return kFALSE;
      if(pt > fCuts->GetProtonPIDthrPtUp()) return kFALSE;
      if(!GoodTPCFitMapSharedMap(aodtrack)) return kFALSE;// Reject tracks with shared clusters
      if(aodtrack->GetTPCNcls()<fCuts->GetProtonTPCCluster()) return kFALSE;
      //Accept track should be here
      if(!AcceptTrack(aodtrack,4)) return kFALSE;//reject tracks with row crossing
      if(!(SelectPID(aodtrack,aodEvent,4))) return kFALSE; //select protons
      if(!GoodDCA(aodtrack,aodEvent,4)) return kFALSE; //to increase the probability for primaries
    }
  else//anti-protons
    {
      if(fabs(eta) > fCuts->GetProtonEtaRange()) return kFALSE; //should lie within the ALICE acceptance
      if(pt < fCuts->GetProtonPIDthrPtLow()) return kFALSE;
      if(pt > fCuts->GetProtonPIDthrPtUp()) return kFALSE;;
      if(!GoodTPCFitMapSharedMap(aodtrack)) return kFALSE;// Reject tracks with shared clusters
      if(aodtrack->GetTPCNcls()<fCuts->GetProtonTPCCluster()) return kFALSE;
      //Accept track should be here
      if(!AcceptTrack(aodtrack,-4)) return kFALSE;//reject tracks with row crossing
      if(!(SelectPID(aodtrack,aodEvent,-4))) return kFALSE; //select protons
      if(!GoodDCA(aodtrack,aodEvent,-4)) return kFALSE; //to increase the probability for primaries
    }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskPLFemto::SingleParticleQualityCuts(const AliAODv0 *v0)
{

  AliAODTrack *posTrack = fGTI[v0->GetPosID()];
  AliAODTrack *negTrack = fGTI[v0->GetNegID()];
  Double_t poseta = posTrack->Eta();
  Double_t negeta = negTrack->Eta();
  Double_t posTPCcls = posTrack->GetTPCNcls();
  Double_t negTPCcls = negTrack->GetTPCNcls();

  //needed?
  //Double_t fracSharedClustersPos = Double_t(daughterTrackPos->GetTPCnclsS()) / Double_t(daughterTrackPos->GetTPCncls());
  //Double_t fracSharedClustersNeg = Double_t(daughterTrackNeg->GetTPCnclsS()) / Double_t(daughterTrackNeg->GetTPCncls());

  if(v0->GetNDaughters() != 2) return kFALSE;
  if(v0->GetNProngs() != 2) return kFALSE;
  if(v0->GetCharge() != 0) return kFALSE;
  if(posTPCcls < fCuts->GetV0TPCCluster()) return kFALSE;
  if(negTPCcls < fCuts->GetV0TPCCluster()) return kFALSE;
  if(fabs(poseta) > fCuts->GetV0EtaRange()) return kFALSE;
  if(fabs(negeta) > fCuts->GetV0EtaRange()) return kFALSE;
  if(v0->Pt() < fCuts->GetV0PIDthrPtLow()) return kFALSE; //below it just enhances background

  // Pile-up rejection cut for RUN2
  if (!fIsRun1 && !posTrack->HasPointOnITSLayer(0) &&
      !posTrack->HasPointOnITSLayer(1) && !posTrack->HasPointOnITSLayer(4) &&
      !posTrack->HasPointOnITSLayer(5) && !(posTrack->GetTOFBunchCrossing() == 0)) return kFALSE;

  if (!fIsRun1 && !negTrack->HasPointOnITSLayer(0) &&
      !negTrack->HasPointOnITSLayer(1) && !negTrack->HasPointOnITSLayer(4) &&
      !negTrack->HasPointOnITSLayer(5) && !(negTrack->GetTOFBunchCrossing() == 0)) return kFALSE;


  return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::ProtonSelector(AliAODEvent *aodEvent)
{
  fProtonCounter = 0;
  fAntiProtonCounter = 0;
  for(Int_t iTrack=0;iTrack<aodEvent->GetNumberOfTracks();iTrack++)
    {
      AliAODTrack *track = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTrack));
      if(!track)
	{
	  AliFatal("Not a standard AOD");
	  continue;
	}

      TString fillthis = "";


      //Several selections
      //if(!(track->GetType() == AliAODTrack::kPrimary)) continue; //should select primary tracks (effect?)

      if(fWhichfilterbit == 128 && (-track->GetID()-1 >= fTrackBuffSize)) continue; //check if the fGTI array is not too small
      if(!track->TestFilterBit(fWhichfilterbit)) continue;
      if(fWhichfilterbit == 128 && !CheckGlobalTrack(track)) continue; //Does a global track exist for the TPC only track?

      //to select only good track conditions (very strict):
      //does only work with filterbit(128) commented out?
      //if((!(track->GetStatus() & AliESDtrack::kITSrefit) || (!(track->GetStatus() & AliESDtrack::kTPCrefit)))) continue;
      //if(!(track->GetStatus()&AliESDtrack::kTPCrefit)) continue;


      Double_t charge = track->Charge();
      if(charge>0.)
	{
	  if(fProtonCounter > kProtonTrackLimit) continue;

	  if(!SingleParticleQualityCuts(track,aodEvent)) continue;//Define additional quality criteria

	  AliAODTrack *globaltrack = fGTI[-track->GetID()-1];

	  fillthis = "fHistTrackFilterMap";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(globaltrack->GetFilterMap());


	  //fill proton parameters:
	  //____________________________________________
	  fProtoncand[fProtonCounter].fPt = track->Pt();
	  fProtoncand[fProtonCounter].fMomentum.SetXYZ(track->Px(),track->Py(),track->Pz());
	  if(fWhichfilterbit == 128) fProtoncand[fProtonCounter].fID = (fGTI[-track->GetID()-1])->GetID();
	  else fProtoncand[fProtonCounter].fID = track->GetID();
	  fProtoncand[fProtonCounter].fEta = track->Eta();
	  fProtoncand[fProtonCounter].fProtonTag = kTRUE;
	  //
	  //GetGlobalPositionAtGlobalRadiiThroughTPC(track,aodEvent->GetMagneticField(),fProtoncand[fProtonCounter].fPrimPosTPC,fPrimVertex);

	  //calculate phi* to check for possible merging/splitting effects:
	  for(int radius=0;radius<9;radius++)
	    {
	      fProtoncand[fProtonCounter].fPhistar[radius] = PhiS(track,aodEvent->GetMagneticField(),fTPCradii[radius]);
	      fProtoncand[fProtonCounter].fPositionTPC[radius] = GetGlobalPositionAtGlobalRadiiThroughTPC(track,aodEvent->GetMagneticField(),fTPCradii[radius],fPrimVertex);
	    }

	  fillthis = "fNProtonsTot";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(1);

	  fillthis = "fProtonPt";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(track->Pt());

	  fillthis = "fProtonPhi";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(track->Phi());

	  fillthis = "fProtonEta";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(track->Eta());

	  //Monte Carlo
	  if(fUseMCInfo)
	    {
	      if(CheckMCPID(aodEvent,track,2212)) fProtoncand[fProtonCounter].fReal = kTRUE;
	      else fProtoncand[fProtonCounter].fReal = kFALSE;

	      //Monte Carlo functions:
	      GetProtonOrigin(track,aodEvent,4);//evaluated with all the cuts applied
	      double protontruth[3] = {-9999.,-9999.,-9999.};
	      double protontruth_mother[3] = {-9999.,-9999.,-9999.};
	      double protontruth_motherParton[3] = {-9999.,-9999.,-9999.};
	      int PDGcodes[3] = {-9999,-9999,-9999};
	      GetMCMomentum(track,protontruth,protontruth_mother,protontruth_motherParton,PDGcodes,aodEvent,4);
	      fProtoncand[fProtonCounter].fMomentumMC.SetXYZ(protontruth[0],protontruth[1],protontruth[2]);
	      fProtoncand[fProtonCounter].fMomentumMCMother.SetXYZ(protontruth_mother[0],protontruth_mother[1],protontruth_mother[2]);
	      fProtoncand[fProtonCounter].fMomentumMCMotherParton.SetXYZ(protontruth_motherParton[0],protontruth_motherParton[1],protontruth_motherParton[2]);
	      if(PDGcodes[0] != -9999) fProtoncand[fProtonCounter].fPDGCodeMother = PDGcodes[0];
	      else fProtoncand[fProtonCounter].fPDGCodeMother = 0;
	      if(PDGcodes[1] != -9999 && PDGcodes[2] != -9999)
		{
		  fProtoncand[fProtonCounter].fPDGCodePartonMother = PDGcodes[1];
		  fProtoncand[fProtonCounter].fPartonMotherLabel = PDGcodes[2];
		}
	      else
		{
		  fProtoncand[fProtonCounter].fPDGCodePartonMother = 0;
		  fProtoncand[fProtonCounter].fPartonMotherLabel = 0;
		}
	    }
	  else
	    {
	      fProtoncand[fProtonCounter].fMomentumMC.SetXYZ(-9999.,-9999.,-9999.);
	      fProtoncand[fProtonCounter].fMomentumMCMother.SetXYZ(-9999.,-9999.,-9999.);
	      fProtoncand[fProtonCounter].fPDGCodePartonMother = -9999;
	      fProtoncand[fProtonCounter].fPartonMotherLabel = -9999;
	    }
	  fProtonCounter++;
	}
      else // Antiprotons
	{
	  if(fAntiProtonCounter > kProtonTrackLimit) continue;
	  if(!SingleParticleQualityCuts(track,aodEvent)) continue;//Define additional quality criteria

	  //fill Anti-proton parameters:
	  //____________________________________________
	  fAntiProtoncand[fAntiProtonCounter].fPt = track->Pt();
	  fAntiProtoncand[fAntiProtonCounter].fMomentum.SetXYZ(track->Px(),track->Py(),track->Pz());
	  if(fWhichfilterbit == 128) fAntiProtoncand[fAntiProtonCounter].fID = (fGTI[-track->GetID()-1])->GetID();
	  else fAntiProtoncand[fAntiProtonCounter].fID = track->GetID();
	  fAntiProtoncand[fAntiProtonCounter].fProtonTag = kTRUE;
	  fAntiProtoncand[fAntiProtonCounter].fEta = track->Eta();

	  //calculate phi* to check for possible merging/splitting effects:
	  for(int radius=0;radius<9;radius++)
	    {
	      fAntiProtoncand[fAntiProtonCounter].fPhistar[radius] = PhiS(track,aodEvent->GetMagneticField(),fTPCradii[radius]);
	      fAntiProtoncand[fAntiProtonCounter].fPositionTPC[radius] = GetGlobalPositionAtGlobalRadiiThroughTPC(track,aodEvent->GetMagneticField(),fTPCradii[radius],fPrimVertex);
	    }


	  if(fUseMCInfo)
	    {
	      //Monte Carlo functions:
	      GetProtonOrigin(track,aodEvent,-4);//evaluated with all the cuts applied
	      double antiprotontruth[3] = {-9999.,-9999.,-9999.};
	      double antiprotontruth_mother[3] = {-9999.,-9999.,-9999.};
	      double antiprotontruth_motherParton[3] = {-9999.,-9999.,-9999.};
	      int PDGcodes[3] = {-9999,-9999,-9999};
	      GetMCMomentum(track,antiprotontruth,antiprotontruth_mother,antiprotontruth_motherParton,PDGcodes,aodEvent,-4);
	      fAntiProtoncand[fAntiProtonCounter].fMomentumMC.SetXYZ(antiprotontruth[0],antiprotontruth[1],antiprotontruth[2]);
	      fAntiProtoncand[fAntiProtonCounter].fMomentumMCMother.SetXYZ(antiprotontruth_mother[0],antiprotontruth_mother[1],antiprotontruth_mother[2]);
	      fAntiProtoncand[fAntiProtonCounter].fMomentumMCMotherParton.SetXYZ(antiprotontruth_motherParton[0],antiprotontruth_motherParton[1],antiprotontruth_motherParton[2]);
	      if(PDGcodes[0] != -9999) fAntiProtoncand[fAntiProtonCounter].fPDGCodeMother = PDGcodes[0];
	      else fProtoncand[fAntiProtonCounter].fPDGCodeMother = 0;
	      if(PDGcodes[1] != -9999 && PDGcodes[2] != -9999)
		{
		  fAntiProtoncand[fAntiProtonCounter].fPDGCodePartonMother = PDGcodes[1];
		  fAntiProtoncand[fAntiProtonCounter].fPartonMotherLabel = PDGcodes[2];
		}
	      else
		{
		  fAntiProtoncand[fAntiProtonCounter].fPDGCodePartonMother = 0;
		  fAntiProtoncand[fAntiProtonCounter].fPartonMotherLabel = 0;
		}
	    }
	  else
	    {
	      fAntiProtoncand[fAntiProtonCounter].fMomentumMC.SetXYZ(-9999.,-9999.,-9999.);
	      fAntiProtoncand[fAntiProtonCounter].fMomentumMCMother.SetXYZ(-9999.,-9999.,-9999.);
	      fAntiProtoncand[fAntiProtonCounter].fPDGCodePartonMother = -9999;
	      fAntiProtoncand[fAntiProtonCounter].fPartonMotherLabel = -9999;
	    }

	  //GetGlobalPositionAtGlobalRadiiThroughTPC(track,aodEvent->GetMagneticField(),fAntiProtoncand[fAntiProtonCounter].fPrimPosTPC,fPrimVertex);
	  fAntiProtonCounter++;
	  fillthis = "fNAntiProtonsTot";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(1);

	  fillthis = "fAntiProtonPt";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(track->Pt());
	}
    }
  //LOOP END (PRIMARY) PROTONS_____________________________________________________________
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::V0Selector(AliAODEvent *aodEvent)
{
  TString fillthis = "";
  fV0Counter = 0;//how many V0 candidates do we have in this event?
  fAntiV0Counter = 0;//how many Anti-V0 candidates do we have in this event?

  //Get a V0 from the event:
  TClonesArray *v01 = (TClonesArray*)aodEvent->GetV0s();
  //number of V0s:
  Int_t entriesV0= v01->GetEntriesFast();


  for(Int_t i=0; i<entriesV0; i++)
    {
      AliAODv0 *v0 = aodEvent->GetV0(i);
      if(!v0) continue;
      if(v0->GetNProngs()>2) continue;

      if(fOnlineV0)
	{
	  if(!(v0->GetOnFlyStatus())) continue; //online v0
	}
      else
	{
	  if((v0->GetOnFlyStatus())) continue; //offline v0
	}

      // Check that the array fGTI isn't too small
      // for the track ids
      if(v0->GetPosID() >= fTrackBuffSize || v0->GetNegID() >= fTrackBuffSize) continue;

      // This is AODs: find the track for given id:
      AliAODTrack *pTrack = fGTI[v0->GetPosID()];
      AliAODTrack *nTrack = fGTI[v0->GetNegID()];

      if((!pTrack) || (!nTrack)) continue;
      if(!SingleParticleQualityCuts(v0)) continue;

      //Probably better to do it in this way (no real difference to the method below):
      //if(!SelectV0PID(v0,fwhichV0) && !SelectV0PID(v0,fwhichAntiV0)) continue;


      //FillV0Selection(v0,aodEvent,fwhichV0);
      //FillV0Selection(v0,aodEvent,fwhichAntiV0);


      if(SelectV0PID(v0,fwhichV0))//V0 selection
	{
	  FillV0Selection(v0,aodEvent,fwhichV0);
	}
      else if(SelectV0PID(v0,fwhichAntiV0))//Anti-V0 selection
	{
	  FillV0Selection(v0,aodEvent,fwhichAntiV0);
	}
    }
}
//_______________________________________________________________________
void AliAnalysisTaskPLFemto::XiSelector(AliAODEvent *aodEvent)
{
  TString fillthis = "";
  //fV0Counter = 0;//how many V0 candidates do we have in this event?
  //fAntiV0Counter = 0;//how many Anti-V0 candidates do we have in this event?

  //Get a V0 from the event:
  TClonesArray *Xicand = (TClonesArray*)aodEvent->GetCascades();
  //number of V0s:
  Int_t entriesXi= Xicand->GetEntriesFast();

  Double_t xvP = aodEvent->GetPrimaryVertex()->GetX();
  Double_t yvP = aodEvent->GetPrimaryVertex()->GetY();
  Double_t zvP = aodEvent->GetPrimaryVertex()->GetZ();


  fXiCounter = 0;
  for (Int_t i=0; i<entriesXi; i++)
    {
      AliAODcascade *Xi = aodEvent->GetCascade(i);
      Double_t pointXi = Xi->CosPointingAngleXi(xvP,yvP,zvP);
      if(!Xi) continue;



      if(Xi->ChargeXi() > 0)//Xip
	{
	}
      else if(Xi->ChargeXi() < 0)//Xim
	{
	  Double_t massXi = Xi->MassXi();

	  AliAODTrack *nTrack2=dynamic_cast<AliAODTrack*>(Xi->GetDaughter(1));
	  AliAODTrack *pTrack=dynamic_cast<AliAODTrack*>(Xi->GetDaughter(0));
	  AliAODTrack *nTrack1=dynamic_cast<AliAODTrack*>(Xi->GetDecayVertexXi()->GetDaughter(0));

	  if(!pTrack || !nTrack1 || !nTrack2 ) continue;

	  if((pTrack->GetID() == nTrack1->GetID()) || (pTrack->GetID() == nTrack2->GetID()) || (nTrack1->GetID() == nTrack2->GetID())) continue;//this would be weird

	  Int_t posTPCcluster = pTrack->GetTPCNcls();
	  Int_t negTPCcluster1 = nTrack1->GetTPCNcls();
	  Int_t negTPCcluster2 = nTrack2->GetTPCNcls();

	  //This selections throw out tracking effects seen in the invariant mass of the Xis
	  if(posTPCcluster < 70) continue;
	  if(negTPCcluster1 < 70) continue;
	  if(negTPCcluster2 < 70) continue;

	  fillthis = "fXiInvMasswoCuts";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(massXi);

	  //Do a rough PID for all daughter tracks:
	  AliPIDResponse::EDetPidStatus statusPosTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,pTrack);
	  if(statusPosTPC != AliPIDResponse::kDetPidOk) continue;
	  AliPIDResponse::EDetPidStatus statusNegTPC1 = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,nTrack1);
	  if(statusNegTPC1 != AliPIDResponse::kDetPidOk) continue;
	  AliPIDResponse::EDetPidStatus statusNegTPC2 = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,nTrack2);
	  if(statusNegTPC2 != AliPIDResponse::kDetPidOk) continue;

	  Double_t nSigmaProtonTPC = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,(AliPID::kProton)));
	  Double_t nSigmaPionTPC1 = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack1,(AliPID::kPion)));
	  Double_t nSigmaPionTPC2 = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack2,(AliPID::kPion)));
	  if(nSigmaProtonTPC > 5.) continue;
	  if(nSigmaPionTPC1 > 5.) continue;
	  if(nSigmaPionTPC2 > 5.) continue;

	  if(pointXi<0.98) continue;
	  Double_t MassV0 = Xi->MassLambda();
	  fillthis = "fXiInvMassLambda";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(MassV0);

	  if(!(MassV0>(fCuts->GetHadronMasses(3122)-fCuts->GetV0PIDWidth()) && MassV0<(fCuts->GetHadronMasses(3122)+fCuts->GetV0PIDWidth()))) continue;
	  fillthis = "fXiInvMasswCuts";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(massXi);

	  if(massXi>(fCuts->GetHadronMasses(3312)-0.005) && massXi<(fCuts->GetHadronMasses(3312)+0.005))
	    {
	      fXicand[fXiCounter].fPt = Xi->Pt();
	      fXicand[fXiCounter].fMass = massXi;
	      fXicand[fXiCounter].fDaughterID1 = Xi->GetPosID();
	      fXicand[fXiCounter].fDaughterID2 = Xi->GetNegID();
	      fXicand[fXiCounter].fBachID = Xi->GetBachID();
	      fXicand[fXiCounter].fXitag = kTRUE;//considered as a good Xi candidate
	      fXicand[fXiCounter].fPointing = pointXi;
	      fXicand[fXiCounter].fMomentum.SetXYZ(Xi->Px(),Xi->Py(),Xi->Pz());
	      fXiCounter++;
	    }
	}

    }
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::TrackCleaner()
{
  TString fillthis = "";
  //Cleaning for V0: There is a probability that two V0s in an event share the same tracks
  Int_t nDoubleShare = 0;
  if(fV0Counter>1)
    {
      //A procedure that controls that we have different V0s:
      for(Int_t i=0;i<fV0Counter;i++)
	{
	  if(!fV0cand[i].fV0tag) continue; //already tagged as bad, must not be tested
	  for(Int_t j=i+1;j<fV0Counter;j++)
	    {
	      if(!fV0cand[j].fV0tag) continue; //already tagged as bad, must not be tested

	      if((fV0cand[i].fDaughterID1 == fV0cand[j].fDaughterID2) || (fV0cand[i].fDaughterID2 == fV0cand[j].fDaughterID1) //crossing, should be excluded already at PID level
		  || (fV0cand[i].fDaughterID1 == fV0cand[j].fDaughterID1) || (fV0cand[i].fDaughterID2 == fV0cand[j].fDaughterID2))//same
		{
		  nDoubleShare++;
		  if(fV0cand[i].fPointing > fV0cand[j].fPointing) fV0cand[j].fV0tag = kFALSE;//V0_i is "more primary" than V0_j
		  else fV0cand[i].fV0tag = kFALSE;
		}
	    }
	}
    }


  if(nDoubleShare>0)
    {
      fillthis = "fNV0TrackSharing";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(nDoubleShare);
    }


  nDoubleShare = 0;
  //Same for Anti V0
  if(fAntiV0Counter>1)
    {
      //A procedure that controls that we have different Anti-V0s:
      for(Int_t i=0;i<fAntiV0Counter;i++)
	{
	  if(!fAntiV0cand[i].fV0tag) continue; //already tagged as bad, must not be tested
	  for(Int_t j=i+1;j<fAntiV0Counter;j++)
	    {
	      if(!fAntiV0cand[j].fV0tag) continue; //already tagged as bad, must not be tested

	      if((fAntiV0cand[i].fDaughterID1 == fAntiV0cand[j].fDaughterID2) || (fAntiV0cand[i].fDaughterID2 == fAntiV0cand[j].fDaughterID1) //crossing, should be excluded already at PID level
		  || (fAntiV0cand[i].fDaughterID1 == fAntiV0cand[j].fDaughterID1) || (fAntiV0cand[i].fDaughterID2 == fAntiV0cand[j].fDaughterID2))//same
		{
		  nDoubleShare++;
		  if(fAntiV0cand[i].fPointing > fAntiV0cand[j].fPointing) fAntiV0cand[j].fV0tag = kFALSE;
		  else fAntiV0cand[i].fV0tag = kFALSE;
		}
	    }
	}
    }

  if(nDoubleShare>0)
    {
      fillthis = "fNAntiV0TrackSharing";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(nDoubleShare);
    }


  //Now check that the proton from V0 and the primary proton are not the same
  for(int i=0; i<fV0Counter;i++)
    {
      if(!fV0cand[i].fV0tag) continue; //already tagged as bad, must not be tested
      for(int j=0; j<fProtonCounter;j++)
	{
	  if((fV0cand[i].fDaughterID1 == fProtoncand[j].fID) || (fV0cand[i].fDaughterID2 == fProtoncand[j].fID))
	    {
	      //Then mark the V0 as bad because the probability is larger that it was a primary proton wrongly paired to a V0
	      fV0cand[i].fV0tag = kFALSE;
	      fillthis = "fNV0protonSharedTracks";
        if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(1);

	      //Test what happens if also the proton is rejected:
	      //fProtoncand[j].fProtonTag = kFALSE;
	    }
	}
    }



  //Check that Anti-V0 and proton don't share tracks
  for(int i=0; i<fAntiV0Counter;i++)
    {
      if(!fAntiV0cand[i].fV0tag) continue; //already tagged as bad, must not be tested
      for(int j=0; j<fProtonCounter;j++)
	{
	  if((fAntiV0cand[i].fDaughterID1 == fProtoncand[j].fID) || (fAntiV0cand[i].fDaughterID2 == fProtoncand[j].fID))
	    {
	      //Then mark the Anti-V0 as bad because the probability is larger that it was a primary proton wrongly paired to a V0
	      fAntiV0cand[i].fV0tag = kFALSE;
	      fillthis = "fNAntiV0protonSharedTracks";
        if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(1);

	      //Test what happens if also the Antiproton is rejected:
	      //fAntiProtoncand[j].fProtonTag = kFALSE;
	    }
	}
    }

  //Check that Anti-Proton and Anti-V0 don't share tracks
  for(int i=0; i<fAntiV0Counter;i++)
    {
      if(!fAntiV0cand[i].fV0tag) continue; //already tagged as bad, must not be tested
      for(int j=0; j<fAntiProtonCounter;j++)
	{
	  if((fAntiV0cand[i].fDaughterID1 == fAntiProtoncand[j].fID) || (fAntiV0cand[i].fDaughterID2 == fAntiProtoncand[j].fID))
	    {
	      //Then mark the Anti-V0 as bad because the probability is larger that it was a primary proton wrongly paired to a V0
	      fAntiV0cand[i].fV0tag = kFALSE;
	      fillthis = "fNAntiV0AntiprotonSharedTracks";
        if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(1);
	    }
	}
    }

  //Check that Xi and proton don't share tracks
  for(int i=0; i<fXiCounter;i++)
    {
      if(!fXicand[i].fXitag) continue; //already tagged as bad, must not be tested
      for(int j=0; j<fProtonCounter;j++)
	{
	  if((fXicand[i].fDaughterID1 == fProtoncand[j].fID) || (fXicand[i].fDaughterID2 == fProtoncand[j].fID) || (fXicand[i].fBachID == fProtoncand[j].fID))
	    {
	      fXicand[i].fXitag = kFALSE;
	    }
	}
    }

  nDoubleShare = 0;
  //Check that Xis don't share tracks
  if(fXiCounter>1)
    {
      //A procedure that controls that we have different V0s:
      for(Int_t i=0;i<fXiCounter;i++)
	{
	  if(!fXicand[i].fXitag) continue; //already tagged as bad, must not be tested
	  for(Int_t j=i+1;j<fXiCounter;j++)
	    {
	      if(!fXicand[j].fXitag) continue; //already tagged as bad, must not be tested

	      if((fXicand[i].fDaughterID1 == fXicand[j].fDaughterID2) || (fXicand[i].fDaughterID2 == fXicand[j].fDaughterID1) //crossing, should be excluded already at PID level
		 || (fXicand[i].fDaughterID1 == fXicand[j].fDaughterID1) || (fXicand[i].fDaughterID2 == fXicand[j].fDaughterID2)//same
		 || (fXicand[i].fBachID == fXicand[j].fBachID)//check that bachelor IDs are different
		 || (fXicand[i].fBachID == fXicand[j].fDaughterID1) || (fXicand[i].fBachID == fXicand[j].fDaughterID2) //check that bacherlor IDs from Xi1 and daughter from Xi2 are different
		 || (fXicand[i].fDaughterID1 == fXicand[j].fBachID) || (fXicand[i].fDaughterID2 == fXicand[j].fBachID))//check that bacherlor IDs from Xi2 and daughter from Xi1 are different
		{
		  nDoubleShare++;
		  if(fXicand[i].fPointing > fXicand[j].fPointing) fXicand[j].fXitag = kFALSE;//V0_i is "more primary" than V0_j
		  else fXicand[i].fXitag = kFALSE;
		}
	    }
	}
    }

  if(nDoubleShare>0)
    {
      fillthis = "fNXimTrackSharing";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(nDoubleShare);
    }


}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::BufferFiller(Int_t zBin,Int_t multBin,Int_t *NumberOfParticles)
{
  fSEPairAnalysisDecider[kProton] = kFALSE;
  fSEPairAnalysisDecider[kAntiProton] = kFALSE;
  fSEPairAnalysisDecider[kV0] = kFALSE;
  fSEPairAnalysisDecider[kAntiV0] = kFALSE;
  fSEPairAnalysisDecider[kXi] = kFALSE;

  TString fillthis = "";

  //Fill V0 candidates in EventCollection:
  Int_t tempV0Counter = 0;
  for(int i=0;i<fV0Counter;i++)
    {
      //one can include additional cuts for V0s at this position
      if(fV0cand[i].fV0tag) //use only V0s which have a good tag (unique V0s without track sharing with other V0s in the sample)
	{
	  fillthis = "fNLambdasTot";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(2); //fill only unique candidates

	  if(fProtonCounter>0)
	    {
	      fillthis = "fNLambdasTot";
        if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(3); //fill only unique candidates
	    }
	  tempV0Counter++;

	  fLambdaTrackVector.push_back(fV0cand[i]);


	  fillthis = "fInvMassLambdaKept";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(fV0cand[i].fMass);
	}
      else
	{
	  fillthis = "fInvMassLambdaRejected";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(fV0cand[i].fMass);
	}
    }


  //If it is not empty put it into Mixing Buffer:
  if(fLambdaTrackVector.size()>0)
    {
      fSEPairAnalysisDecider[kV0] = kTRUE;
      fLambdaEvtBuffer[zBin][multBin].push_front(fLambdaTrackVector);
    }
  //If Buffer is larger than mixing depth, kick out the last element:
  if(fLambdaEvtBuffer[zBin][multBin].size() > kEventsToMix) fLambdaEvtBuffer[zBin][multBin].pop_back();
  fLambdaTrackVector.clear();


  Int_t tempAntiV0Counter = 0;
  for(int i=0;i<fAntiV0Counter;i++)
    {
      //one can include additional cuts for V0s at this position
      if(fAntiV0cand[i].fV0tag) //use only V0s which have a good tag (unique V0s without track sharing with other V0s in the sample)
       {

	  fillthis = "fNAntiLambdasTot";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(2); //fill only unique candidates
	  if(fAntiProtonCounter>0)
	    {
	      fillthis = "fNAntiLambdasTot";
        if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(3); //fill only unique candidates
	    }
	  tempAntiV0Counter++;

	  fAntiLambdaTrackVector.push_back(fAntiV0cand[i]);
	}
    }


  //If it is not empty put it into Mixing Buffer:
  if(fAntiLambdaTrackVector.size()>0)
    {
      fSEPairAnalysisDecider[kAntiV0] = kTRUE;
      fAntiLambdaEvtBuffer[zBin][multBin].push_front(fAntiLambdaTrackVector);
    }
  //If Buffer is larger than mixing depth, kick out the last element:
  if(fAntiLambdaEvtBuffer[zBin][multBin].size() > kEventsToMix) fAntiLambdaEvtBuffer[zBin][multBin].pop_back();
  fAntiLambdaTrackVector.clear();


  //Fill Proton candidates in EventCollection:
  Int_t tempProtonCounter = 0;
  for(Int_t i=0;i<fProtonCounter;i++)
    {
      //at this point you can include some cuts
      //if(fProtoncand[i].fProtonTag)
  {
	  tempProtonCounter++;
	  //Put protons in proton array:
	  fProtonTrackVector.push_back(fProtoncand[i]);
	}
    }

  //If it is not empty put it into Mixing Buffer:
  if(fProtonTrackVector.size()>0)
    {
      fSEPairAnalysisDecider[kProton] = kTRUE;
      fProtonEvtBuffer[zBin][multBin].push_front(fProtonTrackVector);
    }
  //If Buffer is larger than mixing depth, kick out the last element:
  if(fProtonEvtBuffer[zBin][multBin].size() > kEventsToMix) fProtonEvtBuffer[zBin][multBin].pop_back();
  fProtonTrackVector.clear();


  //Fill Anti-Proton candidates in EventCollection:
  Int_t tempAntiProtonCounter = 0;
  for(Int_t i=0;i<fAntiProtonCounter;i++)
    {
      //if(fAntiProtoncand[i].fProtonTag)
  {
	  tempAntiProtonCounter++;
	  //Put protons in proton array:
	  fAntiProtonTrackVector.push_back(fAntiProtoncand[i]);
	}
    }

  //If it is not empty put it into Mixing Buffer:
  if(fAntiProtonTrackVector.size()>0)
    {
      fSEPairAnalysisDecider[kAntiProton] = kTRUE;
      fAntiProtonEvtBuffer[zBin][multBin].push_front(fAntiProtonTrackVector);
    }
  //If Buffer is larger than mixing depth, kick out the last element:
  if(fAntiProtonEvtBuffer[zBin][multBin].size() > kEventsToMix) fAntiProtonEvtBuffer[zBin][multBin].pop_back();
  fAntiProtonTrackVector.clear();


  Int_t tempXiCounter = 0;
  for(int i=0;i<fXiCounter;i++)
    {
      if(fXicand[i].fXitag)
  {
    tempXiCounter++;
    fXiTrackVector.push_back(fXicand[i]);
  }
    }
  //If it is not empty put it into Mixing Buffer:
  if(fXiTrackVector.size()>0)
    {
      fSEPairAnalysisDecider[kXi] = kTRUE;
      fXiEvtBuffer[zBin][multBin].push_front(fXiTrackVector);
    }
  //If Buffer is larger than mixing depth, kick out the last element:
  if(fXiEvtBuffer[zBin][multBin].size() > kEventsToMix) fXiEvtBuffer[zBin][multBin].pop_back();
  fXiTrackVector.clear();

  NumberOfParticles[0] = tempV0Counter;
  NumberOfParticles[1] = tempProtonCounter;
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::ParticlePairer(Int_t zBin,Int_t multBin)
{
  //In this function two different event mixing procedures are performed.
  //The "standard" one which includes sometimes a buffering of emtpy evts
  //A second one where empty events are avoided
  //The one where empty events are avoided fills histograms with using the QA enums

  Int_t NPairs = 0;
  //Loop for p-V0:
  for(unsigned int i=0; i<(fLambdaEvtBuffer[zBin][multBin].front()).size();i++)//get same event size
    {
      if(!fSEPairAnalysisDecider[kV0]) break;

      for(unsigned int evnum=0; evnum<fProtonEvtBuffer[zBin][multBin].size(); evnum++)
	{
	  if(evnum == 0) NPairs++;
	  if(evnum == 0 && !(fSEPairAnalysisDecider[kProton] && fSEPairAnalysisDecider[kV0])) continue;//this ensures that same event correlations is only performed if a proton and V0 is in same event

	  for(unsigned int j=0;j<(fProtonEvtBuffer[zBin][multBin].at(evnum)).size();j++)
	    {
	      PairAnalysis((fLambdaEvtBuffer[zBin][multBin].front()).at(i),(fProtonEvtBuffer[zBin][multBin].at(evnum)).at(j),evnum,kProtonV0);
	    }//ProtonA Loop
	}//event Loop
    }//ProtonB Loop


  TString fillthis = "fNProtonLambdaPairs";
  if(NPairs !=0) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(NPairs);
  /*
  //Loop for p-AntiV0
  for(int i=0; i<(fEvt)->fNumAntiV0s;i++)
    {
      for(int evnum=0; evnum<kEventsToMix+1; evnum++)
	{
	  for(int j=0; j<(fEvt+evnum)->fNumProtons;j++)
	    {
	      PairAnalysis((fEvt)->fAntiLambdaParticle[i],(fEvt+evnum)->fProtonParticle[j],evnum,kProtonAntiV0);
	    }//Proton Loop
	}//event Loop
    }//V0 Loop


  //Loop for Antip-V0
  for(int i=0; i<(fEvt)->fNumV0s;i++)
    {
      for(int evnum=0; evnum<kEventsToMix+1; evnum++)
	{
	  for(int j=0; j<(fEvt+evnum)->fNumAntiProtons;j++)
	    {
	      PairAnalysis((fEvt)->fLambdaParticle[i],(fEvt+evnum)->fAntiProtonParticle[j],evnum,kAntiProtonV0);
	    }//AntiProton Loop
	}//event Loop
    }//V0 Loop
  */

  //Loop for Antip-AntiV0
  for(unsigned int i=0; i<(fAntiLambdaEvtBuffer[zBin][multBin].front()).size();i++)//get same event size
    {
      if(!fSEPairAnalysisDecider[kAntiV0]) break;

      for(unsigned int evnum=0; evnum<fAntiProtonEvtBuffer[zBin][multBin].size(); evnum++)
	{

	  if(evnum == 0 && !(fSEPairAnalysisDecider[kAntiProton] && fSEPairAnalysisDecider[kAntiV0])) continue;//this ensures that same event correlations is only performed if a proton and V0 is in same event

	  for(unsigned int j=0;j<(fAntiProtonEvtBuffer[zBin][multBin].at(evnum)).size();j++)
	    {
	      PairAnalysis((fAntiLambdaEvtBuffer[zBin][multBin].front()).at(i),(fAntiProtonEvtBuffer[zBin][multBin].at(evnum)).at(j),evnum,kAntiProtonAntiV0);
	    }//ProtonA Loop
	}//event Loop
    }//ProtonB Loop


  //Loop for p-p
  if(fSEPairAnalysisDecider[kProton])
    {
      //Check second option
      for(unsigned int i=0; i<(fProtonEvtBuffer[zBin][multBin].front()).size();i++)//get same event size
	{
	  for(unsigned int evnum=0; evnum<fProtonEvtBuffer[zBin][multBin].size(); evnum++)
	    {
	      int startbin = 0;
	      if(evnum==0)
		{
		  startbin=i+1;
		  //if((fProtonEvtBuffer.front()).at(i).fID == (fProtonEvtBuffer.at(evnum)).at(startbin).fID) continue;
		}
	      for(unsigned int j=startbin;j<(fProtonEvtBuffer[zBin][multBin].at(evnum)).size();j++)
		{
		  PairAnalysis((fProtonEvtBuffer[zBin][multBin].front()).at(i),(fProtonEvtBuffer[zBin][multBin].at(evnum)).at(j),evnum,kProtonProton);
		}//ProtonA Loop
	    }//event Loop
	}//ProtonB Loop
    }


  //Loop for Antip-Antip
  if(fSEPairAnalysisDecider[kAntiProton])
    {
      //Check second option
      for(unsigned int i=0; i<(fAntiProtonEvtBuffer[zBin][multBin].front()).size();i++)//get same event size
	{
	  for(unsigned int evnum=0; evnum<fAntiProtonEvtBuffer[zBin][multBin].size(); evnum++)
	    {
	      int startbin = 0;
	      if(evnum==0)
		{
		  startbin=i+1;
		  //if((fProtonEvtBuffer.front()).at(i).fID == (fProtonEvtBuffer.at(evnum)).at(startbin).fID) continue;
		}
	      for(unsigned int j=startbin;j<(fAntiProtonEvtBuffer[zBin][multBin].at(evnum)).size();j++)
		{
		  PairAnalysis((fAntiProtonEvtBuffer[zBin][multBin].front()).at(i),(fAntiProtonEvtBuffer[zBin][multBin].at(evnum)).at(j),evnum,kAntiProtonAntiProton);
		}//ProtonA Loop
	    }//event Loop
	}//ProtonB Loop
    }

  /*
  //Loop for Antip-p
  for(int i=0; i<(fEvt)->fNumProtons;i++)
    {
      for(int evnum=0; evnum<kEventsToMix+1; evnum++)
	{
	  for(int j=0; j<(fEvt+evnum)->fNumAntiProtons;j++)
	    {
	      if(evnum==0)
		{
		  if((fEvt)->fProtonParticle[i].fID == (fEvt+evnum)->fAntiProtonParticle[j].fID) continue; //should only be set for same events
		}
	      PairAnalysis((fEvt)->fProtonParticle[i],(fEvt+evnum)->fAntiProtonParticle[j],evnum,kAntiProtonProton);
	    }//AntiProton Loop
	}//event Loop
    }//Proton Loop

  //Loop for p-Xi
  for(int i=0; i<(fEvt)->fNumXis;i++)
    {
      for(int evnum=0; evnum<kEventsToMix+1; evnum++)
	{
	  for(int j=0; j<(fEvt+evnum)->fNumProtons;j++)
	    {
	      PairAnalysis((fEvt)->fXiParticle[i],(fEvt+evnum)->fProtonParticle[j],evnum,kProtonXi);
	    }//Proton Loop
	}//event Loop
    }//Xi Loop


  //Loop for Antip-Xi
  for(int i=0; i<(fEvt)->fNumXis;i++)
    {
      for(int evnum=0; evnum<kEventsToMix+1; evnum++)
	{
	  for(int j=0; j<(fEvt+evnum)->fNumAntiProtons;j++)
	    {
	      PairAnalysis((fEvt)->fXiParticle[i],(fEvt+evnum)->fAntiProtonParticle[j],evnum,kAntiProtonXi);
	    }//AntiProton Loop
	}//event Loop
    }//Xi Loop

  //Loop for V0-V0
  for(int i=0; i<(fEvt)->fNumV0s;i++)
    {
      for(int evnum=0; evnum<kEventsToMix+1; evnum++)
	{
	  int startbin = 0;
	  if(evnum==0) startbin=i+1;
	  for(int j=startbin; j<(fEvt+evnum)->fNumV0s;j++)
	    {
	      PairAnalysis((fEvt)->fLambdaParticle[i],(fEvt+evnum)->fLambdaParticle[j],evnum,kV0V0);
	    }//Proton Loop
	}//event Loop
    }//V0 Loop
  */

  if(fSEPairAnalysisDecider[kV0])
    {
      //Check second option
      for(unsigned int i=0; i<(fLambdaEvtBuffer[zBin][multBin].front()).size();i++)//get same event size
	{
	  for(unsigned int evnum=0; evnum<fLambdaEvtBuffer[zBin][multBin].size(); evnum++)
	    {
	      int startbin = 0;
	      if(evnum==0) startbin = i+1;
	      for(unsigned int j=startbin;j<(fLambdaEvtBuffer[zBin][multBin].at(evnum)).size();j++)
		{
		  PairAnalysis((fLambdaEvtBuffer[zBin][multBin].front()).at(i),(fLambdaEvtBuffer[zBin][multBin].at(evnum)).at(j),evnum,kV0V0);
		}//V0A Loop
	    }//event Loop
	}//V0B Loop
    }


  //Loop for AntiV0-AntiV0
  if(fSEPairAnalysisDecider[kAntiV0])
    {
      //Check second option
      for(unsigned int i=0; i<(fAntiLambdaEvtBuffer[zBin][multBin].front()).size();i++)//get same event size
	{
	  for(unsigned int evnum=0; evnum<fAntiLambdaEvtBuffer[zBin][multBin].size(); evnum++)
	    {
	      int startbin = 0;
	      if(evnum==0) startbin = i+1;
	      for(unsigned int j=startbin;j<(fAntiLambdaEvtBuffer[zBin][multBin].at(evnum)).size();j++)
		{
		  PairAnalysis((fAntiLambdaEvtBuffer[zBin][multBin].front()).at(i),(fAntiLambdaEvtBuffer[zBin][multBin].at(evnum)).at(j),evnum,kAntiV0AntiV0);
		}//V0A Loop
	    }//event Loop
	}//V0B Loop
    }

  /*
  //Loop for AntiV0-V0
  for(int i=0; i<(fEvt)->fNumAntiV0s;i++)
    {
      for(int evnum=0; evnum<kEventsToMix+1; evnum++)
	{
	  for(int j=0; j<(fEvt+evnum)->fNumV0s;j++)
	    {
	      PairAnalysis((fEvt)->fAntiLambdaParticle[i],(fEvt+evnum)->fLambdaParticle[j],evnum,kAntiV0V0);
	    }//Proton Loop
	}//event Loop
    }//V0 Loop
  */


  //Loop for p-Xi
  if(fSEPairAnalysisDecider[kXi])
    {
      //Check second option
      for(unsigned int i=0; i<(fXiEvtBuffer[zBin][multBin].front()).size();i++)//get same event size
  {
    for(unsigned int evnum=0; evnum<fProtonEvtBuffer[zBin][multBin].size(); evnum++)
      {
        int startbin = 0;
        if(evnum==0) startbin = i+1;
        for(unsigned int j=startbin;j<(fProtonEvtBuffer[zBin][multBin].at(evnum)).size();j++)
    {
      PairAnalysis((fXiEvtBuffer[zBin][multBin].front()).at(i),(fProtonEvtBuffer[zBin][multBin].at(evnum)).at(j),evnum,kProtonXi);
    }//V0A Loop
      }//event Loop
  }//V0B Loop
    }

  //Loop for Antip-Xi
  if(fSEPairAnalysisDecider[kXi])
    {
      //Check second option
      for(unsigned int i=0; i<(fXiEvtBuffer[zBin][multBin].front()).size();i++)//get same event size
  {
    for(unsigned int evnum=0; evnum<fAntiProtonEvtBuffer[zBin][multBin].size(); evnum++)
      {
        int startbin = 0;
        if(evnum==0) startbin = i+1;
        for(unsigned int j=startbin;j<(fAntiProtonEvtBuffer[zBin][multBin].at(evnum)).size();j++)
    {
      PairAnalysis((fXiEvtBuffer[zBin][multBin].front()).at(i),(fAntiProtonEvtBuffer[zBin][multBin].at(evnum)).at(j),evnum,kAntiProtonXi);
    }//V0A Loop
      }//event Loop
  }//V0B Loop
    }

}
//________________________________________________________________________
Float_t AliAnalysisTaskPLFemto::PhiS(AliAODTrack *track,const Float_t bfield,const Float_t Radius,const Float_t DecVtx[2])
{
  const Float_t chg = track->Charge();
  const Float_t Pt = track->Pt();
  const Float_t phi0 = track->Phi();
  const Float_t magCurv = Pt/(chg*0.1*0.3*bfield);
  const Float_t xD = DecVtx[0];
  const Float_t yD = DecVtx[1];

  //now the ugly solution of the system of equations:
  Float_t firstTerm = (magCurv + xD)*(Radius*Radius + 2*magCurv*xD + xD*xD + yD*yD);
  Float_t secondTerm = -yD*yD*((2*magCurv-Radius+xD)*(Radius+xD) + yD*yD)*((xD-Radius)*(2*magCurv+Radius + xD)+yD*yD);

  if(secondTerm<0) return -9999.;

  Float_t nominator = firstTerm - TMath::Sqrt(secondTerm);
  Float_t denominator = 2*((magCurv+xD)*(magCurv+xD) + yD*yD);

  //intersection point:
  Float_t interX = nominator/denominator;

  Float_t phis = phi0 + TMath::ASin(interX/Radius);

  return phis;
}
//________________________________________________________________________
Float_t AliAnalysisTaskPLFemto::PhiS(AliAODTrack *track,const Float_t bfield,const Float_t Radius)
{
  const Float_t phi0 = track->Phi();//angle at primary vertex
  const Float_t pt = track->Pt();
  const Float_t chg = track->Charge();

  //To use the following equation:
  //pt must be given in GeV/c
  //bfield in T
  //chg in units of electric charge
  //radius in m

  //0.3 is a conversion factor for pt and bfield can be plugged in in terms of GeV/c and electric charge, 0.1 converts the magnetic field to Tesla, 0.01 transforms the radius from cm to m
  Float_t phis = phi0 + TMath::ASin(0.1*chg*bfield*0.3*Radius*0.01/(2.*pt));

  return phis;
}
//________________________________________________________________________
Float_t AliAnalysisTaskPLFemto::EtaS(const Float_t zVertex, const Float_t Radius)
{
  Float_t thetaS = TMath::Pi()/2. - TMath::ATan(zVertex/Radius);

  return -TMath::Log( TMath::Tan(thetaS/2.) );
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPLFemto::ClosePairRejecter(Float_t deltaPhiS,Float_t deltaEta)
{
  Bool_t isTooClose = kFALSE;

  if(fabs(deltaPhiS)<0.045 && fabs(deltaEta)<0.01) isTooClose = kTRUE;

  return isTooClose;
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::UserExec(Option_t *)
{
  // user exec
  if (!fInputEvent) {
    Error("UserExec","NO EVENT FOUND!");
    return;
  }

  fCuts->SetCutVariation(fcutType);

  TString fillthis = "";
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  fAODMCEvent = MCEvent(); //return NULL pointer in case of exp event, no problem

  if(!fIsLightweight) fCEvents->Fill(1);

  if(fIsRun1) {
      // the AODs with null vertex pointer didn't pass the PhysSel
      if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return;
      if(!fIsLightweight) fCEvents->Fill(2);

      // AOD primary vertex
      AliAODVertex *vtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();
      if(!vtx1) return;
      if(vtx1->GetNContributors()<2) return;
      if(!fIsLightweight) fCEvents->Fill(3);
      if(!fIsLightweight) fPrimVertex[0] = vtx1->GetX();
      if(!fIsLightweight) fPrimVertex[1] = vtx1->GetY();
      if(!fIsLightweight) fPrimVertex[2] = vtx1->GetZ();

      fillthis = "fVtxX";
      if(!fIsLightweight) ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(fPrimVertex[0]);
      fillthis = "fVtxY";
      if(!fIsLightweight) ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(fPrimVertex[1]);
      fillthis = "fVtxZ";
      if(!fIsLightweight) ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(fPrimVertex[2]);


      if(fPrimVertex[2] < fCuts->GetEvtzVertexLow() || fPrimVertex[2] > fCuts->GetEvtzVertexUp()) return;
      if(!fIsLightweight) fCEvents->Fill(4);
    }

  // AliEventCuts for Run2 analysis
  if(!fIsRun1) {
      bool isSelected = fAliEventCuts.AcceptEvent(fInputEvent);
      if(!isSelected) return;

      if(fTrigger == AliVEvent::kHighMultV0 && fV0PercentileMax < 100.f) {
          Float_t lPercentile = 300;
          AliMultSelection *MultSelection = 0x0;
          MultSelection = (AliMultSelection *)fInputEvent->FindListObject("MultSelection");
          if (!MultSelection) {
              // If you get this warning (and lPercentiles 300) please check that the
              // AliMultSelectionTask actually ran (before your task)
              AliWarning("AliMultSelection object not found!");
            } else {
              lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
            }
          if (lPercentile > fV0PercentileMax) return;
        }
    }

  //find the centrality and zBin for mixed event: 0=zBin, 1=MultBin
  Int_t *MixingBins = MultiplicityBinFinderzVertexBinFinder(aodEvent,fPrimVertex[2]);
  if(MixingBins[0]==-9999 || MixingBins[1]==-9999) return;//Multiplicity estimation fails
  if(!fIsLightweight) fCEvents->Fill(5);
  Int_t zBin = MixingBins[0];
  Int_t MultBin = MixingBins[1];
  fillthis = "fNBinsMultmixing";
  if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(MultBin);
  fillthis = "fNBinsVertexmixing";
  if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(zBin);

  if(fIsRun1) {
    if(!fIsLightweight) if(!aodEvent->IsPileupFromSPD(3,0.8,3.,2.,5.)) fCEvents->Fill(7);

    if(fAnaUtils->IsPileUpEvent(aodEvent)) return;
    if(!fIsLightweight) fCEvents->Fill(6);
  }

  //*******************************(O.Arnold)
  //Set the private pointer array which contains information from global tracks
  ResetGlobalTrackReference();
  for(Int_t iTrack=0;iTrack<aodEvent->GetNumberOfTracks();iTrack++)
    {
      AliAODTrack *track = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTrack));
      if(!track) AliFatal("Not a standard AOD");
      if(!track) continue;

      // Store the reference of the global tracks
      StoreGlobalTrackReference(track);
    }
  //*******************************
  fEventNumber++;

  if(!fIsLightweight) fCEvents->Fill(9);

  //FIFOShifter(zBin,MultBin);//Put the information into the FIFO
  ProtonSelector(aodEvent); //Select primary protons
  V0Selector(aodEvent); //Select V0s
  XiSelector(aodEvent); //Select Xis

  //totally empty event for our purposes? then neglect it:
  if(fProtonCounter == 0 && fV0Counter == 0 && fAntiProtonCounter == 0 && fAntiV0Counter == 0 && fXiCounter == 0) return;

  //if(fProtonCounter == 0) return;

  //if(fProtonCounter>0 && fV0Counter>0) {fCEvents->Fill(8);}
  //if(fAntiProtonCounter>0 && fAntiV0Counter>0) {fCEvents->Fill(9);}

  TrackCleaner(); //Check for unique V0s and no track sharing between V0 daughter tracks and primary tracks
  Int_t NumberOfParticles[2];
  BufferFiller(zBin,MultBin,NumberOfParticles);//Fill the candidates to the buffer: 0=V0,1=Proton
  ParticlePairer(zBin,MultBin); //at this point the relevant distributions for the CF are built and mixing takes place

  fillthis = "fNProtonsPerevent";
  if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(NumberOfParticles[1]);

  fillthis = "fNLambdasPerevent";
  if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(NumberOfParticles[0]); //fill only unique candidates

  PostData(1,fOutput);
  PostData(2,fOutputSP);
  PostData(3,fOutputPID);
  PostData(4,fOutputTP);
  PostData(5,fOutputAliEvent);
}
//________________________________________________________________________
Double_t AliAnalysisTaskPLFemto::Qinv(TLorentzVector trackV0,TLorentzVector trackProton)
{
  // Copied from Hans Beck analysis who copied it from NA49:
  // Copied from NA49. See http://na49info.web.cern.ch/na49info/na49/Software/minidst/ana/html/src/T49Tool.cxx.html#T49Tool:Qinv
  //  Always using lambda mass (no mass difference found yet for lam <-> alam (see PDG))

  TLorentzVector trackQ, trackP;
  trackQ = trackV0 - trackProton;
  trackP = trackV0 + trackProton;

  Double_t qinvL = trackQ.Mag2();
  Double_t qP = trackQ.Dot(trackP);
  Double_t pinv = trackP.Mag2();

  Double_t QinvLP = TMath::Sqrt(qP*qP/pinv - qinvL);

  return QinvLP;
}
//________________________________________________________________________
Double_t AliAnalysisTaskPLFemto::relKcalc(TLorentzVector track1,TLorentzVector track2)
{
  //This function calculates the relative momentum k* between any particle pair

  TLorentzVector trackSum, track1cms, track2cms;
  trackSum = track1 + track2;

  Double_t beta = trackSum.Beta();
  Double_t betax = beta*cos(trackSum.Phi())*sin(trackSum.Theta());
  Double_t betay = beta*sin(trackSum.Phi())*sin(trackSum.Theta());
  Double_t betaz = beta*cos(trackSum.Theta());

  track1cms = track1;
  track2cms = track2;

  track1cms.Boost(-betax,-betay,-betaz);
  track2cms.Boost(-betax,-betay,-betaz);

  TLorentzVector trackRelK;

  trackRelK = track1cms - track2cms;
  Double_t relK = 0.5*trackRelK.P();

  return relK;
}
//________________________________________________________________________
Double_t AliAnalysisTaskPLFemto::GetAverageSeparation(const TVector3 globalPositions1st[],const TVector3 globalPositions2nd[])
{
  double sumSquare = 0.;
  int pointUsed = 0;

  for(int RNumber = 0; RNumber < 9; RNumber++)
    {
      TVector3 diffVec = globalPositions1st[RNumber] - globalPositions2nd[RNumber];
      if(globalPositions1st[RNumber].X() == -9999. || globalPositions2nd[RNumber].X() == -9999.) continue;//This happens if the distances deviate too strongly

      sumSquare += diffVec.Mag();
      pointUsed++;
    }

  if(pointUsed < 1) return 0.;
  else return sumSquare / pointUsed;
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::PairAnalysis(const AliFemtoLambdaParticle &v0,const AliFemtoProtonParticle &proton,Int_t REorME,pairType inputPair)
{
  //This function analyses p-V0 correlations

  TLorentzVector trackV0,trackProton,trackSum,trackV0posdaughter;
  trackV0.SetXYZM(v0.fMomentum.X(),v0.fMomentum.Y(),v0.fMomentum.Z(),fCuts->GetHadronMasses(3122));
  trackV0posdaughter.SetXYZM(v0.fMomentumPosDaughter.X(),v0.fMomentumPosDaughter.Y(),v0.fMomentumPosDaughter.Z(),fCuts->GetHadronMasses(2212));
  trackProton.SetXYZM(proton.fMomentum.X(),proton.fMomentum.Y(),proton.fMomentum.Z(),fCuts->GetHadronMasses(2212));

  trackSum = trackV0 + trackProton;

  Double_t relK = relKcalc(trackV0,trackProton);
  Double_t relKPPdaughter = relKcalc(trackProton,trackV0posdaughter);
  Double_t pairKT = 0.5*trackSum.Pt();
  double averageMass = 0.5*(fCuts->GetHadronMasses(3122) + fCuts->GetHadronMasses(2212));
  double pairMT = TMath::Sqrt(pow(pairKT,2.) + pow(averageMass,2.));

  double qinv = Qinv(trackV0,trackProton);
  Double_t avgSepPos = GetAverageSeparation(proton.fPositionTPC,v0.fPositionPosTPC);
  Double_t avgSepNeg = GetAverageSeparation(proton.fPositionTPC,v0.fPositionNegTPC);

  TString fillthis = "fProtonLambdaDiff";
  if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(fabs(relK - 0.5*qinv));

  Int_t *MixingBins = MultiplicityBinFinderzVertexBinFinder(static_cast<AliAODEvent*>(fInputEvent),fPrimVertex[2]);
  Int_t MultBin = MixingBins[1];

  if(REorME == 0 && inputPair == kProtonV0)//real event
    {
      fillthis = "fPairkTLp";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(pairKT);

      fillthis = "fProtonLambdaMT";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(pairMT);

      fillthis = "fProtonLambdaRelK";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fProtonLambdaRelKMulti_";
      fillthis += MultBin;
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fProtonLambdaPosDaughterRelK";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relKPPdaughter);

      fillthis = "fNPairStatistics";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(1);

      if(relK<0.15)
	{
	  fillthis = "fNPairStatistics";
	  ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(2);

	  fillthis = "fProtonLambdaMTrelKcut";
	  ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(pairMT);
	}

      //check for splitting/merging in angle space
      for(int radius=0; radius<9; radius++) {
        Double_t deltaEtaPos = proton.fEta - v0.fEtaPosdaughter;
        Double_t deltaPhistarPos = proton.fPhistar[radius] - v0.fPhiStarPosdaughter[radius];
        Double_t deltaEtaNeg = proton.fEta - v0.fEtaNegdaughter;
        Double_t deltaPhistarNeg = proton.fPhistar[radius] - v0.fPhiStarNegdaughter[radius];

        fillthis = "fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughProton";
        fillthis += radius;
        if(!fIsLightweight) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(TMath::Abs(deltaEtaPos),TMath::Abs(deltaPhistarPos));

        fillthis = "fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughPion";
        fillthis += radius;
        if(!fIsLightweight) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(TMath::Abs(deltaEtaNeg),TMath::Abs(deltaPhistarNeg));
      }

      fillthis = "fProtonLambdaAvgSeparationPP";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPos);

      fillthis = "fProtonLambdaAvgSeparationPPi";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepNeg);
    }
  else if(REorME!=0 && inputPair == kProtonV0)//mixed event
    {
      fillthis = "fProtonLambdaRelKME";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fProtonLambdaRelKMEMulti_";
      fillthis += MultBin;
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fProtonLambdaPosDaughterRelKME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relKPPdaughter);

      fillthis = "fProtonLambdaAvgSeparationPPME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPos);

      fillthis = "fProtonLambdaAvgSeparationPPiME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepNeg);

      //check for splitting/merging in angle space
      for(int radius=0; radius<9; radius++) {
        Double_t deltaEtaPos = proton.fEta - v0.fEtaPosdaughter;
        Double_t deltaPhistarPos = proton.fPhistar[radius] - v0.fPhiStarPosdaughter[radius];
        Double_t deltaEtaNeg = proton.fEta - v0.fEtaNegdaughter;
        Double_t deltaPhistarNeg = proton.fPhistar[radius] - v0.fPhiStarNegdaughter[radius];

        fillthis = "fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughProtonME";
        fillthis += radius;
        if(!fIsLightweight) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(TMath::Abs(deltaEtaPos),TMath::Abs(deltaPhistarPos));

        fillthis = "fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughPionME";
        fillthis += radius;
        if(!fIsLightweight) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(TMath::Abs(deltaEtaNeg),TMath::Abs(deltaPhistarNeg));
      }

      if(fUseMCInfo) GetMomentumMatrix(v0,proton);//determines the amount of momentum resolution (from mixed event to enlarge statistics)
    }
  else if(REorME == 0 && inputPair == kAntiProtonV0)//real event
    {
      fillthis = "fAntiProtonLambdaRelK";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);
    }
  //else if(REorME != 0 && whichpair_combi=="AntiProton-V0")//mixed event
  else if(REorME != 0 && inputPair == kAntiProtonV0)//mixed event
    {
      fillthis = "fAntiProtonLambdaRelKME";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);
    }
  else if(REorME == 0 && inputPair == kProtonAntiV0)//real event
    {
      fillthis = "fProtonAntiLambdaRelK";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);
    }
  else if(REorME != 0 && inputPair == kProtonAntiV0)//mixed event
    {
      fillthis = "fProtonAntiLambdaRelKME";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);
    }
  else if(REorME == 0 && inputPair == kAntiProtonAntiV0)//real event
    {
      fillthis = "fAntiProtonAntiLambdaRelK";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fAntiProtonAntiLambdaRelKMulti_";
      fillthis += MultBin;
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fNPairStatistics";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(3);

      fillthis = "fAntiProtonAntiLambdaPosDaughterRelK";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relKPPdaughter);

      fillthis = "fAntiProtonAntiLambdaAvgSeparationPP";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPos);

      fillthis = "fAntiProtonAntiLambdaAvgSeparationPPi";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepNeg);


      if(relK<0.15)
	{
	  fillthis = "fNPairStatistics";
	  ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(4);
	}
    }
  else if(REorME != 0 && inputPair == kAntiProtonAntiV0)//mixed event
    {
      fillthis = "fAntiProtonAntiLambdaRelKME";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fAntiProtonAntiLambdaRelKMEMulti_";
      fillthis += MultBin;
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fAntiProtonAntiLambdaPosDaughterRelKME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relKPPdaughter);

      fillthis = "fAntiProtonAntiLambdaAvgSeparationPPME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPos);

      fillthis = "fAntiProtonAntiLambdaAvgSeparationPPiME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepNeg);
    }
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::PairAnalysis(const AliFemtoXiParticle &xi,const AliFemtoProtonParticle &proton,Int_t REorME,pairType inputPair)
{
  //This method analyses p-Xi correlations

  TLorentzVector trackXi,trackProton;
  trackXi.SetXYZM(xi.fMomentum.X(),xi.fMomentum.Y(),xi.fMomentum.Z(),fCuts->GetHadronMasses(3312));
  trackProton.SetXYZM(proton.fMomentum.X(),proton.fMomentum.Y(),proton.fMomentum.Z(),fCuts->GetHadronMasses(2212));

  Double_t relK = relKcalc(trackXi,trackProton);

  Int_t *MixingBins = MultiplicityBinFinderzVertexBinFinder(static_cast<AliAODEvent*>(fInputEvent),fPrimVertex[2]);
  Int_t MultBin = MixingBins[1];

  if(REorME == 0 && inputPair == kProtonXi)//real event
    {
      TString fillthis = "fProtonXiRelK";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fProtonXiRelKMulti_";
      fillthis += MultBin;
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);
    }
  else if(REorME!=0 && inputPair == kProtonXi)//mixed event
    {
      TString fillthis = "fProtonXiRelKME";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fProtonXiRelKMEMulti_";
      fillthis += MultBin;
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);
    }
  else if(REorME == 0 && inputPair == kAntiProtonXi)
    {
      TString fillthis = "fAntiProtonXiRelK";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fAntiProtonXiRelKMulti_";
      fillthis += MultBin;
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);
    }
  else if(REorME != 0 && inputPair == kAntiProtonXi)
    {
      TString fillthis = "fAntiProtonXiRelKME";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fAntiProtonXiRelKMEMulti_";
      fillthis += MultBin;
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);
    }
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::PairAnalysis(const AliFemtoLambdaParticle &v01,const AliFemtoLambdaParticle &v02,Int_t REorME,pairType inputPair)
{
  TLorentzVector trackV01,trackV02;

  trackV01.SetXYZM(v01.fMomentum.X(),v01.fMomentum.Y(),v01.fMomentum.Z(),fCuts->GetHadronMasses(3122));
  trackV02.SetXYZM(v02.fMomentum.X(),v02.fMomentum.Y(),v02.fMomentum.Z(),fCuts->GetHadronMasses(3122));

  Double_t relK = relKcalc(trackV01,trackV02);

  Double_t avgSepPP = GetAverageSeparation(v01.fPositionPosTPC,v02.fPositionPosTPC);
  Double_t avgSepPPi = GetAverageSeparation(v01.fPositionPosTPC,v02.fPositionNegTPC);
  Double_t avgSepPiP = GetAverageSeparation(v01.fPositionNegTPC,v02.fPositionPosTPC);
  Double_t avgSepPiPi = GetAverageSeparation(v01.fPositionNegTPC,v02.fPositionNegTPC);

  Int_t *MixingBins = MultiplicityBinFinderzVertexBinFinder(static_cast<AliAODEvent*>(fInputEvent),fPrimVertex[2]);
  Int_t MultBin = MixingBins[1];

  TString fillthis = "";

  if(REorME == 0 && inputPair == kV0V0)//real event
    {
      fillthis = "fLambdaLambdaRelK";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fLambdaLambdaRelKMulti_";
      fillthis += MultBin;
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fLambdaLambdaAvgSeparationPP";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPP);

      fillthis = "fLambdaLambdaAvgSeparationP1Pi2";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPPi);

      fillthis = "fLambdaLambdaAvgSeparationP2Pi1";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPiP);

      fillthis = "fLambdaLambdaAvgSeparationPiPi";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPiPi);
    }
  else if(REorME != 0 && inputPair == kV0V0)//mixed event
    {
      fillthis = "fLambdaLambdaRelKME";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fLambdaLambdaRelKMEMulti_";
      fillthis += MultBin;
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fLambdaLambdaAvgSeparationPPME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPP);

      fillthis = "fLambdaLambdaAvgSeparationP1Pi2ME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPPi);

      fillthis = "fLambdaLambdaAvgSeparationP2Pi1ME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPiP);

      fillthis = "fLambdaLambdaAvgSeparationPiPiME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPiPi);

      if(fUseMCInfo) GetMomentumMatrix(v01, v02);
    }
  else if(REorME == 0 && inputPair == kAntiV0AntiV0)//real event
    {
      fillthis = "fAntiLambdaAntiLambdaRelK";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fAntiLambdaAntiLambdaRelKMulti_";
      fillthis += MultBin;
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fAntiLambdaAntiLambdaAvgSeparationPP";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPP);

      fillthis = "fAntiLambdaAntiLambdaAvgSeparationP1Pi2";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPPi);

      fillthis = "fAntiLambdaAntiLambdaAvgSeparationP2Pi1";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPiP);

      fillthis = "fAntiLambdaAntiLambdaAvgSeparationPiPi";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPiPi);
    }
  else if(REorME != 0 && inputPair == kAntiV0AntiV0)//mixed event
    {
      fillthis = "fAntiLambdaAntiLambdaRelKME";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fAntiLambdaAntiLambdaRelKMEMulti_";
      fillthis += MultBin;
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fAntiLambdaAntiLambdaAvgSeparationPPME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPP);

      fillthis = "fAntiLambdaAntiLambdaAvgSeparationP1Pi2ME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPPi);

      fillthis = "fAntiLambdaAntiLambdaAvgSeparationP2Pi1ME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPiP);

      fillthis = "fAntiLambdaAntiLambdaAvgSeparationPiPiME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSepPiPi);

      if(fUseMCInfo) GetMomentumMatrix(v01, v02);
    }
  else if(REorME == 0 && inputPair == kAntiV0V0)//real event
    {
      fillthis = "fAntiLambdaLambdaRelK";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);
    }
  else if(REorME != 0 && inputPair == kAntiV0V0)//mixed event
    {
      fillthis = "fAntiLambdaLambdaRelKME";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);
    }
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::PairAnalysis(const AliFemtoProtonParticle &protonA,const AliFemtoProtonParticle &protonB,Int_t REorME,pairType inputPair)
{
  //This method analyses p-p correlations

  TLorentzVector trackProtonA,trackProtonB,trackSum;
  TString fillthis = "";

  trackProtonA.SetXYZM(protonA.fMomentum.X(),protonA.fMomentum.Y(),protonA.fMomentum.Z(),fCuts->GetHadronMasses(2212));
  trackProtonB.SetXYZM(protonB.fMomentum.X(),protonB.fMomentum.Y(),protonB.fMomentum.Z(),fCuts->GetHadronMasses(2212));

  trackSum = trackProtonA + trackProtonB;

  Double_t relK = relKcalc(trackProtonA,trackProtonB);
  Double_t pairKT = 0.5*trackSum.Pt();
  Double_t pairMT = TMath::Sqrt(pow(pairKT,2.) + pow(fCuts->GetHadronMasses(2212),2.));

  Double_t e1Timese2 = trackProtonA.E() * trackProtonB.E();
  Double_t e1Pluse2 = trackProtonA.E() + trackProtonB.E();
  Double_t p1zTimesp2z = trackProtonA.Pz() * trackProtonB.Pz();

  Double_t avgSep = GetAverageSeparation(protonA.fPositionTPC,protonB.fPositionTPC);

  TVector2 tPt1;
  TVector2 tPt2;
  tPt1.Set(trackProtonA.Px(),trackProtonA.Py());
  tPt2.Set(trackProtonB.Px(),trackProtonB.Py());
  double pt1Dotpt2 = tPt1*tPt2;

  Int_t *MixingBins = MultiplicityBinFinderzVertexBinFinder(static_cast<AliAODEvent*>(fInputEvent),fPrimVertex[2]);
  Int_t MultBin = MixingBins[1];

  if(REorME == 0 && inputPair == kProtonProton)//real event
    {
      fillthis = "fProtonProtonRelK";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fProtonProtonRelKMulti_";
      fillthis += MultBin;
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fPairkTpp";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(pairKT);

      fillthis = "fProtonProtonMT";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(pairMT);

      fillthis = "fNPairStatistics";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(6);

      if(relK < 0.15)
	{
	  fillthis = "fProtonProtonMTrelKcut";
	  ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(pairMT);

	  fillthis = "fNPairStatistics";
	  ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(7);
	}


      //check for splitting/merging in angle space

      for(int radius=0; radius<9; radius++)
	{
	  Double_t deltaEta = protonA.fEta - protonB.fEta;

	  Double_t deltaPhistar = protonA.fPhistar[radius] - protonB.fPhistar[radius];


	  fillthis = "fDeltaEtaDeltaPhiTPCradpp";
	  fillthis += radius;
    if(!fIsLightweight) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(TMath::Abs(deltaEta),TMath::Abs(deltaPhistar));
	}

      //Fill only with certain separation
      if(avgSep > 4.)//cm
	{
	  fillthis = "fProtonProtonRelKw2TC";
    if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);
	}


      fillthis = "fProtonProtonAvgSeparation";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSep);

      fillthis = "fProtonProtonAvgSeparationVsRelK";
      if(!fIsLightweight) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSep,relK);

      fillthis = "fProtonPtCorrelationRE";
      if(fUseMCInfo) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(trackProtonA.Pt(),trackProtonB.Pt());

      fillthis = "fProtonEtaCorrelationRE";
      if(fUseMCInfo) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(trackProtonA.PseudoRapidity(),trackProtonB.PseudoRapidity());

      fillthis = "fProtonEnCorrelationRE";
      if(fUseMCInfo) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(trackProtonA.E(),trackProtonB.E());

      fillthis = "fProtonPhiCorrelationRE";
      if(fUseMCInfo) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(trackProtonA.Phi(),trackProtonB.Phi());
    }
  else if(REorME!=0 && inputPair == kProtonProton)//mixed event
    {
      fillthis = "fProtonProtonRelKME";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fProtonProtonRelKMEMulti_";
      fillthis += MultBin;
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      //EMCIC histos:
      fillthis = "fProtonProtonRelKEnMultME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK,e1Timese2);

      fillthis = "fProtonProtonRelKEnSumME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK,e1Pluse2);

      fillthis = "fProtonProtonRelKPzMultME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK,p1zTimesp2z);

      fillthis = "fProtonProtonRelKPtMultME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK,pt1Dotpt2);

      if(fUseMCInfo) GetMomentumMatrix(protonA,protonB);//Determine the momentum resolution (from mixed event to enlarge statistics)


      //check for splitting/merging in angle space
      for(int radius=0; radius<9; radius++)
	{
	  Double_t deltaEta = protonA.fEta - protonB.fEta;
	  Double_t deltaPhistar = protonA.fPhistar[radius] - protonB.fPhistar[radius];

	  fillthis = "fDeltaEtaDeltaPhiTPCradppME";
	  fillthis += radius;
    if(!fIsLightweight) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(TMath::Abs(deltaEta),TMath::Abs(deltaPhistar));
	}

      //Fill p-p correlations with a certain track separation:
      if(avgSep > 4.)//cm
	{
	  fillthis = "fProtonProtonRelKw2TCME";
    if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);
	}

      fillthis = "fProtonProtonAvgSeparationME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSep);

      fillthis = "fProtonProtonAvgSeparationVsRelKME";
      if(!fIsLightweight) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSep,relK);

      fillthis = "fProtonPtCorrelationME";
      if(fUseMCInfo) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(trackProtonA.Pt(),trackProtonB.Pt());

      fillthis = "fProtonEtaCorrelationME";
      if(fUseMCInfo) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(trackProtonA.PseudoRapidity(),trackProtonB.PseudoRapidity());

      fillthis = "fProtonEnCorrelationME";
      if(fUseMCInfo) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(trackProtonA.E(),trackProtonB.E());

      fillthis = "fProtonPhiCorrelationME";
      if(fUseMCInfo) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(trackProtonA.Phi(),trackProtonB.Phi());
    }
  else if(REorME==0 && inputPair == kAntiProtonAntiProton)//same event
    {
      fillthis = "fAntiProtonAntiProtonRelK";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fAntiProtonAntiProtonRelKMulti_";
      fillthis += MultBin;
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fPairkTApAp";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(pairKT);

      fillthis = "fNPairStatistics";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(8);

      if(relK < 0.15)
	{
	  fillthis = "fNPairStatistics";
	  ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(9);
	}

      fillthis = "fAntiProtonAntiProtonAvgSeparation";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSep);

      //Fill only with certain separation
      if(avgSep > 4.)//cm
	{
	  fillthis = "fAntiProtonAntiProtonRelKw2TC";
    if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);
	}

    }
  else if(REorME!=0 && inputPair == kAntiProtonAntiProton)//mixed event
    {
      fillthis = "fAntiProtonAntiProtonRelKME";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fAntiProtonAntiProtonRelKMEMulti_";
      fillthis += MultBin;
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fAntiProtonAntiProtonAvgSeparationME";
      if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(avgSep);

      //Fill only with certain separation
      if(avgSep > 4.)//cm
	{
	  fillthis = "fAntiProtonAntiProtonRelKw2TCME";
    if(!fIsLightweight) ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);
	}
    }
  else if(REorME==0 && inputPair == kAntiProtonProton)//same event
    {
      fillthis = "fAntiProtonProtonRelK";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);

      fillthis = "fPairkTApp";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(pairKT);
    }
  else if(REorME!=0 && inputPair == kAntiProtonProton)//mixed event
    {
      fillthis = "fAntiProtonProtonRelKME";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK);
    }
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::TripleAnalysis(const AliFemtoLambdaParticle &v0,const AliFemtoProtonParticle &proton1,const AliFemtoProtonParticle &proton2,Int_t REorME,pairType inputPair)
{
  TLorentzVector trackV0;
  TLorentzVector trackProton1;
  TLorentzVector trackProton2;

  trackProton1.SetXYZM(proton1.fMomentum.X(),proton1.fMomentum.Y(),proton1.fMomentum.Z(),fCuts->GetHadronMasses(2212));
  trackProton2.SetXYZM(proton2.fMomentum.X(),proton2.fMomentum.Y(),proton2.fMomentum.Z(),fCuts->GetHadronMasses(2212));
  trackV0.SetXYZM(v0.fMomentum.X(),v0.fMomentum.Y(),v0.fMomentum.Z(),fCuts->GetHadronMasses(3122));


  Double_t relK12 = relKcalc(trackProton1,trackProton2);
  Double_t relK13 = relKcalc(trackProton1,trackV0);
  Double_t relK23 = relKcalc(trackProton2,trackV0);

  Double_t relK123 = TMath::Sqrt(pow(relK12,2.) + pow(relK13,2.) + pow(relK23,2.));

  TString fillthis = "";

  if(REorME == 0 && inputPair == kProtonProtonV0)//real event
    {
      fillthis = "fProtonProtonLambdaRelK";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK123);
    }
  else if(REorME != 0 && inputPair == kProtonProtonV0)
    {
      fillthis = "fProtonProtonLambdaRelKME";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK123);
    }
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::TripleAnalysis(const AliFemtoProtonParticle &proton1,const AliFemtoProtonParticle &proton2,const AliFemtoProtonParticle &proton3,Int_t REorME,pairType inputPair)
{
  //This method analyses p-p-p correlations

  TLorentzVector trackProton1;
  TLorentzVector trackProton2;
  TLorentzVector trackProton3;

  trackProton1.SetXYZM(proton1.fMomentum.X(),proton1.fMomentum.Y(),proton1.fMomentum.Z(),fCuts->GetHadronMasses(2212));
  trackProton2.SetXYZM(proton2.fMomentum.X(),proton2.fMomentum.Y(),proton2.fMomentum.Z(),fCuts->GetHadronMasses(2212));
  trackProton3.SetXYZM(proton3.fMomentum.X(),proton3.fMomentum.Y(),proton3.fMomentum.Z(),fCuts->GetHadronMasses(2212));

  Double_t relK12 = relKcalc(trackProton1,trackProton2);
  Double_t relK13 = relKcalc(trackProton1,trackProton3);
  Double_t relK23 = relKcalc(trackProton2,trackProton3);

  Double_t relK123 = TMath::Sqrt(pow(relK12,2.) + pow(relK13,2.) + pow(relK23,2.));

  TString fillthis = "";

  if(REorME == 0 && inputPair == kProtonProtonProton)//real event
    {
      fillthis = "fProtonProtonProtonRelK";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK123);
    }
  else if(REorME != 0 && inputPair == kProtonProtonProton)
    {
      fillthis = "fProtonProtonProtonRelKME";
      ((TH1F*)(fOutputTP->FindObject(fillthis)))->Fill(relK123);
    }
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPLFemto::V0TopologicalSelection(AliAODv0 *v0,AliAODEvent *aodEvent,TString ParOrAPar)
{
  //This function selects the V0s according to given criteria which were studied with Monte Carlo

  TString fillthis = "";

  // Get the coordinates of the primary vertex
  Double_t xvP = aodEvent->GetPrimaryVertex()->GetX();
  Double_t yvP = aodEvent->GetPrimaryVertex()->GetY();
  Double_t zvP = aodEvent->GetPrimaryVertex()->GetZ();
  Double_t vecTarget[3]={xvP,yvP,zvP};
  Double_t vV0[3];

  v0->GetXYZ(vV0);
  Double_t V0pt = v0->Pt();

  //Calculate vertex variables:
  Double_t dcaV0Dau = v0->DcaV0Daughters();
  Double_t dcaPrim = v0->DcaV0ToPrimVertex();
  Double_t dcaPrimPos = v0->DcaPosToPrimVertex();
  Double_t dcaPrimNeg = v0->DcaNegToPrimVertex();
  Double_t lenDecay = v0->DecayLengthV0(vecTarget);
  Double_t point = v0->CosPointingAngle(vecTarget);
  Double_t transverseRadius = v0->DecayLengthXY(vecTarget);


  if(ParOrAPar == fwhichV0)
    {
      fillthis = "fLambdaDCADaughterTracks";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaV0Dau);

      fillthis = "fLambdaDCAPrimVertex";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaPrim);

      fillthis = "fLambdaDCAPosdaughPrimVertex";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaPrimPos);

      fillthis = "fLambdaDCANegdaughPrimVertex";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaPrimNeg);

      fillthis = "fLambdaDecaylength";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(lenDecay);

      fillthis = "fLambdaCosPointingAngle";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(point);

      fillthis = "fLambdaTransverseRadius";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(transverseRadius);
    }
  else if(ParOrAPar == fwhichAntiV0)
    {
      //conditions for AntiV0 must be implemented!
    }


  if(fUseMCInfo && ParOrAPar == fwhichV0)
    {
      Bool_t realV0 = CheckIfRealV0(v0,aodEvent,3122,2212,211);

      if(realV0)
	{
	  fillthis = "fDCAV0PosdaughMC";
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaPrimPos,V0pt);

	  fillthis = "fDCAV0NegdaughMC";
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaPrimNeg,V0pt);

	  fillthis = "fTrRadiusV0MC";
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(transverseRadius,V0pt);

	  fillthis = "fPoinAngleV0MC";
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(point,V0pt);

	  fillthis = "fLambdaV0vertexXMC";
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(vV0[0],V0pt);

	  fillthis = "fLambdaV0vertexYMC";
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(vV0[1],V0pt);

	  fillthis = "fLambdaV0vertexZMC";
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(vV0[2],V0pt);

	  fillthis = "fDecayLengthV0MC";
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(lenDecay,V0pt);
	}
      else
	{
	  fillthis = "fDCAbkgPosdaughMC";
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaPrimPos,V0pt);

	  fillthis = "fDCAbkgNegdaughMC";
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(dcaPrimNeg,V0pt);

	  fillthis = "fTrRadiusbkgMC";
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(transverseRadius,V0pt);

	  fillthis = "fPoinAnglebkgMC";
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(point,V0pt);

	  fillthis = "fLambdabkgvertexXMC";
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(vV0[0],V0pt);

	  fillthis = "fLambdabkgvertexYMC";
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(vV0[1],V0pt);

	  fillthis = "fLambdabkgvertexZMC";
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(vV0[2],V0pt);

	  fillthis = "fDecayLengthbkgMC";
    if(!fIsLightweight) ((TH2F*)(fOutputSP->FindObject(fillthis)))->Fill(lenDecay,V0pt);
	}
    }

  //Check if the V0 fulfills the topological cut criteria:
  if(fabs(vV0[0]) > fCuts->GetV0decayvtxcut()) return kFALSE;
  if(fabs(vV0[1]) > fCuts->GetV0decayvtxcut()) return kFALSE;
  if(fabs(vV0[2]) > fCuts->GetV0decayvtxcut()) return kFALSE;

  if(transverseRadius < fCuts->GetV0rxylow()) return kFALSE;
  if(transverseRadius > fCuts->GetV0rxyup()) return kFALSE;
  if(dcaPrimPos < fCuts->GetV0DCAtrPV()) return kFALSE;
  if(dcaPrimNeg < fCuts->GetV0DCAtrPV())  return kFALSE;
  if(dcaV0Dau > fCuts->GetV0DCAtracksV0decay()) return kFALSE;


  //Fill the pointing angle as a function of pt
  int ptbin = fCuts->FindV0PtBinCPA(V0pt);

  if(ParOrAPar == fwhichV0 && ptbin != -9999)
    {
      Double_t MassV0 = v0->MassLambda();

      fillthis = "fLambdaCPAPtBin";
      fillthis += ptbin;
      if(MassV0>(fCuts->GetHadronMasses(3122)-fCuts->GetV0PIDWidth()) && MassV0<(fCuts->GetHadronMasses(3122)+fCuts->GetV0PIDWidth()) && ParOrAPar == fwhichV0)
	{
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(point);//Fill pointing angle only for Lambda candidates
	}

      if(fUseMCInfo)
	{
	  Bool_t realV0 = CheckIfRealV0(v0,aodEvent,3122,2212,211);

	  if(realV0)//these are true MC V0s
	    {
	      TClonesArray *mcarray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	      if (!mcarray) return kFALSE;
	      Int_t daugh[2] = {2212,211};//daughter IDs of Lambda
	      Int_t label = v0->MatchToMC(3122,mcarray,2,daugh);
	      AliAODMCParticle *MCv0 = (AliAODMCParticle*)mcarray->At(label);
	      if(!MCv0) return kFALSE;;

	      fillthis = "fLambdaCPAAllPtBin";
	      fillthis += ptbin;
        if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(point);

	      if(MCv0->IsPhysicalPrimary() && !MCv0->IsSecondaryFromWeakDecay())
		{
		  fillthis = "fLambdaCPAPrimaryPtBin";
		  fillthis += ptbin;
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(point);
		}
	      else if(MCv0->IsSecondaryFromWeakDecay())
		{
		  fillthis = "fLambdaCPASecondaryPtBin";
		  fillthis += ptbin;
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(point);
		}
	      else if(MCv0->IsSecondaryFromMaterial())
		{
		  fillthis = "fLambdaCPAMaterialPtBin";
		  fillthis += ptbin;
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(point);
		}
	    }
	  else
	    {
	      fillthis = "fLambdaCPABkgPtBin";
	      fillthis += ptbin;
	      //background in the mass range of the lambda is interesting for us:
	      if(MassV0>(fCuts->GetHadronMasses(3122)-fCuts->GetV0PIDWidth()) && MassV0<(fCuts->GetHadronMasses(3122)+fCuts->GetV0PIDWidth()) && ParOrAPar == fwhichV0)
		{
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(point);//Fill background only in the selection region
		}
	    }

	}
    }


  if(point < fCuts->GetV0pointing()) return kFALSE;

  return kTRUE;
}
//________________________________________________________________________
inline void AliAnalysisTaskPLFemto::FillV0Selection(AliAODv0 *v0,AliAODEvent* aodEvent,TString ParOrAPar)
{
  Double_t MassV0 = 0.;

  TString fillthis = "";


  if(ParOrAPar == fwhichV0)
    {
      Double_t MassK0 = v0->MassK0Short();
      MassV0 = v0->MassLambda();
      fillthis = "fInvMassMissIDK0s";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(MassK0);

      if(MassK0 < fCuts->GetK0sAntiCutLow() || fCuts->GetK0sAntiCutUp() < MassK0)
	{
	  //This cuts out the peak of the K0s
	  fillthis = "fInvMassMissIDK0swCuts";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(MassK0);
	}
      else return;

    }
  else if(ParOrAPar == fwhichAntiV0)
    {
      MassV0 = v0->MassAntiLambda();
      Double_t MassK0 = v0->MassK0Short();
      fillthis = "fInvMassMissIDK0s";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(MassK0);

      if(MassK0 < fCuts->GetK0sAntiCutLow() || fCuts->GetK0sAntiCutUp() < MassK0)
	{
	  fillthis = "fInvMassMissIDK0swCuts";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(MassK0);
	}
      else return;
    }
  else if(ParOrAPar == "Kaon" || ParOrAPar == "Anti-Kaon") MassV0 = v0->MassK0Short();


  Bool_t V0Cuts = V0TopologicalSelection(v0,aodEvent,ParOrAPar);//calculates the topological cut conditions
  Double_t xvP = aodEvent->GetPrimaryVertex()->GetX();
  Double_t yvP = aodEvent->GetPrimaryVertex()->GetY();
  Double_t zvP = aodEvent->GetPrimaryVertex()->GetZ();
  Double_t vecTarget[3]={xvP,yvP,zvP};
  Double_t point = v0->CosPointingAngle(vecTarget);

  if(ParOrAPar == fwhichV0)
    {
      fillthis = "fInvMassLambdawoCuts";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(MassV0);
    }
  else if(ParOrAPar == fwhichAntiV0)
    {
      fillthis = "fInvMassAntiLambdawoCuts";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(MassV0);
    }
  else
    {
      std::cout << "Specify which V0 you want to analyse" << std::endl;
      return;
    }

  if(V0Cuts)// apply vertex cuts
    {
      AliAODTrack* pTrack = fGTI[v0->GetPosID()];
      AliAODTrack* nTrack = fGTI[v0->GetNegID()];

      if(!pTrack || !nTrack) return;

      Double_t charge1 = pTrack->Charge();
      Double_t charge2 = nTrack->Charge();

      if(charge1 < 0. && charge2 > 0. && ParOrAPar == fwhichV0)//assign correct charge to tracks from V0
	{
	  pTrack = fGTI[v0->GetNegID()];//proton as positive particle
	  nTrack = fGTI[v0->GetPosID()];//pion as negative particle
	}
      else if(charge1 > 0. && charge2 <0. && ParOrAPar == fwhichAntiV0)//assign correct charge to tracks from Anti-V0
	{
	  pTrack = fGTI[v0->GetNegID()];//anti-proton as negative particle
	  nTrack = fGTI[v0->GetPosID()];//anti-pion as positive particle
	}

      // Lambda or AntiLambda
      if(ParOrAPar == fwhichV0)
	{
	  fillthis = "fInvMassLambdawCuts";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(MassV0);

	  if(fProtonCounter>0)
	    {
	      fillthis = "fInvMassLambdawCutsAfterSelection";
        if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(MassV0);
	    }

	  Int_t ptbin = fCuts->FindV0PtBinInvMass(v0->Pt());

	  if(ParOrAPar == fwhichV0 && ptbin != -9999)
	    {
	      fillthis = "fInvMassLambdawCutsPtBin";
	      fillthis += ptbin;
        if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(MassV0);
	    }
	}
      else if(ParOrAPar == fwhichAntiV0)
	{
	  fillthis = "fInvMassAntiLambdawCuts";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(MassV0);
	}

      //Save all the information of the V0
      //_____________________________________________________
      Bool_t which_V0_region = kFALSE;
      Bool_t which_AntiV0_region = kFALSE;

      if(fwhichV0region == "left")
	{
	  if(MassV0<1.105 && ParOrAPar == fwhichV0) which_V0_region = kTRUE;
	  else if(MassV0<1.105 && ParOrAPar == fwhichAntiV0) which_AntiV0_region = kTRUE;
	}
      else if(fwhichV0region == "right")
	{
	  if(MassV0>1.13 && MassV0<1.18 && ParOrAPar == fwhichV0) which_V0_region = kTRUE;
	  else if(MassV0>1.13 && MassV0<1.18 && ParOrAPar == fwhichAntiV0) which_AntiV0_region = kTRUE;
	}
      else if(fwhichV0region == "signal")
	{
	  if(MassV0>(fCuts->GetHadronMasses(3122)-fCuts->GetV0PIDWidth()) && MassV0<(fCuts->GetHadronMasses(3122)+fCuts->GetV0PIDWidth()) && ParOrAPar == fwhichV0) which_V0_region = kTRUE;
	  else if(MassV0>(fCuts->GetHadronMasses(3122)-fCuts->GetV0PIDWidth()) && MassV0<(fCuts->GetHadronMasses(3122)+fCuts->GetV0PIDWidth()) && ParOrAPar == fwhichAntiV0) which_AntiV0_region = kTRUE;
	}


      if(which_V0_region && ParOrAPar == fwhichV0)
	{
	  fillthis = "fLambdaPt";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(v0->Pt());

	  fillthis = "fLambdaPhi";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(v0->Phi());

	  fillthis = "fLambdaEta";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(v0->Eta());

	  fillthis = "fInvMassLambdawCutsAfterMassCut";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(MassV0);

	  fillthis = "fLambdaPosDaughPt";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(pTrack->Pt());

	  fillthis = "fLambdaNegDaughPt";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(nTrack->Pt());

	  if(fV0Counter > kV0TrackLimit) return;
	  if(pTrack->GetID() == nTrack->GetID()) return; //should never happen
	  fV0cand[fV0Counter].fPt = v0->Pt();
	  fV0cand[fV0Counter].fMomentum.SetXYZ(v0->Px(),v0->Py(),v0->Pz());
	  fV0cand[fV0Counter].fMass = MassV0;
	  fV0cand[fV0Counter].fPointing = point;
	  fV0cand[fV0Counter].fEtaPosdaughter = pTrack->Eta();//eta should be translational invariant due to magnetic field along z (?)
	  fV0cand[fV0Counter].fV0tag = kTRUE;//should be by true by initialization but do it again (however i had problems with this flag, maybe order of initialization went sometimes wrong)
	  fV0cand[fV0Counter].fDaughterID1 = pTrack->GetID();
	  fV0cand[fV0Counter].fDaughterID2 = nTrack->GetID();
	  fV0cand[fV0Counter].fMomentumPosDaughter.SetXYZ(pTrack->Px(),pTrack->Py(),pTrack->Pz());
	  fV0cand[fV0Counter].fMomentumNegDaughter.SetXYZ(nTrack->Px(),nTrack->Py(),nTrack->Pz());
    //GetGlobalPositionAtGlobalRadiiThroughTPC(pTrack,aodEvent->GetMagneticField(),fV0cand[fV0Counter].fPositionPosTPC,fPrimVertex);
    //GetGlobalPositionAtGlobalRadiiThroughTPC(nTrack,aodEvent->GetMagneticField(),fV0cand[fV0Counter].fPositionNegTPC,fPrimVertex);

	  for(int radius=0;radius<9;radius++)
      {
        fV0cand[fV0Counter].fPhiStarPosdaughter[radius] = PhiS(pTrack,aodEvent->GetMagneticField(),fTPCradii[radius]);
        fV0cand[fV0Counter].fPhiStarNegdaughter[radius] = PhiS(nTrack,aodEvent->GetMagneticField(),fTPCradii[radius]);
	      fV0cand[fV0Counter].fPositionPosTPC[radius] = GetGlobalPositionAtGlobalRadiiThroughTPC(pTrack,aodEvent->GetMagneticField(),fTPCradii[radius],fPrimVertex);
	      fV0cand[fV0Counter].fPositionNegTPC[radius] = GetGlobalPositionAtGlobalRadiiThroughTPC(nTrack,aodEvent->GetMagneticField(),fTPCradii[radius],fPrimVertex);
	    }



	  if(fUseMCInfo)// if monte carlo
	    {
	      if(CheckMCPID(v0,aodEvent,3122)) fV0cand[fV0Counter].fReal = kTRUE;
	      else fV0cand[fV0Counter].fReal = kFALSE;

	      /*
	      if(fV0cand[fV0Counter].fReal)
		{
		  fillthis = "fLambdaSignalBkg";
		  ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(1);
		}
	      else
		{
		  fillthis = "fLambdaSignalBkg";
		  ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(2);
		}
	      */
	      GetV0Origin(v0,aodEvent);//evaluated with all the cuts applied
	      double v0MCmom[3] = {-9999.,-9999.,-9999.};
	      double v0MCmom_mother[4] = {-9999.,-9999.,-9999.,-9999.};
	      GetMCMomentum(v0,v0MCmom,v0MCmom_mother,aodEvent);
	      fV0cand[fV0Counter].fMomentumMC.SetXYZ(v0MCmom[0],v0MCmom[1],v0MCmom[2]);
	      fV0cand[fV0Counter].fMomentumMCMother.SetXYZ(v0MCmom_mother[0],v0MCmom_mother[1],v0MCmom_mother[2]);
	      if(v0MCmom_mother[3]>0) fV0cand[fV0Counter].fPDGCodeMother = (int)v0MCmom_mother[3];//it stems coming from a weak decay
	      else fV0cand[fV0Counter].fPDGCodeMother = 0;
	    }
	  else
	    {
	      fV0cand[fV0Counter].fMomentumMC.SetXYZ(-9999.,-9999.,-9999.);
	      fV0cand[fV0Counter].fMomentumMCMother.SetXYZ(-9999.,-9999.,-9999.);
	    }

	  if(!fV0cand[fV0Counter].fV0tag) std::cout << "something is wrong with the V0 tag" << std::endl;

	  fillthis = "fNLambdasTot";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(1); //fill only unique candidates

	  fV0Counter++;//This is a V0 (Lambda) candidate
	}
      else if(which_AntiV0_region && ParOrAPar == fwhichAntiV0)//fill only in a certain mass range
	{
	  fillthis = "fAntiLambdaPt";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(v0->Pt());

	  if(fAntiV0Counter > kV0TrackLimit) return;
	  if(pTrack->GetID() == nTrack->GetID()) return; //should never happen
	  fAntiV0cand[fAntiV0Counter].fPt = v0->Pt();
	  fAntiV0cand[fAntiV0Counter].fMomentum.SetXYZ(v0->Px(),v0->Py(),v0->Pz());
	  fAntiV0cand[fAntiV0Counter].fMass = MassV0;
	  fAntiV0cand[fAntiV0Counter].fPointing = point;

	  fAntiV0cand[fAntiV0Counter].fV0tag = kTRUE;
	  fAntiV0cand[fAntiV0Counter].fDaughterID1 = pTrack->GetID();
	  fAntiV0cand[fAntiV0Counter].fDaughterID2 = nTrack->GetID();
	  fAntiV0cand[fAntiV0Counter].fMomentumPosDaughter.SetXYZ(pTrack->Px(),pTrack->Py(),pTrack->Pz());
	  fAntiV0cand[fAntiV0Counter].fMomentumNegDaughter.SetXYZ(nTrack->Px(),nTrack->Py(),nTrack->Pz());

	  fillthis = "fNAntiLambdasTot";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(1); //fill only unique candidates


	  for(int radius=0;radius<9;radius++)
	    {
	      fAntiV0cand[fAntiV0Counter].fPositionPosTPC[radius] = GetGlobalPositionAtGlobalRadiiThroughTPC(pTrack,aodEvent->GetMagneticField(),fTPCradii[radius],fPrimVertex);
	      fAntiV0cand[fAntiV0Counter].fPositionNegTPC[radius] = GetGlobalPositionAtGlobalRadiiThroughTPC(nTrack,aodEvent->GetMagneticField(),fTPCradii[radius],fPrimVertex);
	    }

	  fAntiV0Counter++;//This is a Anti-V0 (Anti-Lambda) candidate
	}
      //_____________________________________________________
    }
}

//________________________________________ terminate ___________________________
void AliAnalysisTaskPLFemto::Terminate(Option_t*)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  //Info("Terminate","");
  AliAnalysisTaskSE::Terminate();

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }

  //fCEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fCEvents"));

  fOutputSP = dynamic_cast<TList*> (GetOutputData(2)); //1
  if (!fOutputSP) {
    printf("ERROR: fOutputAll not available\n");
    return;
  }
  fOutputPID = dynamic_cast<TList*> (GetOutputData(3));//2
  if (!fOutputPID) {
    printf("ERROR: fOutputPID not available\n");
    return;
  }

  fOutputTP = dynamic_cast<TList*> (GetOutputData(4));//3
  if (!fOutputTP) {
    printf("ERROR: fOutputTP not available\n");
    return;
  }

  return;
}
//___________________________________________________________________________
void AliAnalysisTaskPLFemto::UserCreateOutputObjects()
{
  std::cout << "Create Output Objects" << std::endl;

 // output
  Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());

  //slot #1
  //OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("listEvt");

  if(!fIsRun1) {
      fOutputAliEvent = new TList();
      fOutputAliEvent->SetOwner();
      fOutputAliEvent->SetName("listAliEventCuts");
      if(fTrigger != AliVEvent::kINT7) {
        fAliEventCuts.SetManualMode();
        fAliEventCuts.SetupRun2pp();
        fAliEventCuts.fTriggerMask = fTrigger;
      }
      if(!fIsLightweight) fAliEventCuts.AddQAplotsToList(fOutputAliEvent);
    }

  fOutputPID = new TList();
  fOutputPID->SetOwner();
  fOutputPID->SetName("listPID");

  fOutputSP = new TList();
  fOutputSP->SetOwner();
  fOutputSP->SetName("listSP");

  fOutputTP = new TList();
  fOutputTP->SetOwner();
  fOutputTP->SetName("listTP");

  // define histograms
  DefineHistograms(fwhichV0);


  fGTI = new AliAODTrack *[fTrackBuffSize]; // Array of pointers

  //for particle identification
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if(!man){AliError("Couldn't get the analysis manager!");}
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if(!inputHandler){AliError("Couldn't get the input handler!");}
  fPIDResponse = inputHandler->GetPIDResponse();
  if(!fPIDResponse){AliError("Couldn't get the PID response task!");}


  PostData(1,fOutput);
  PostData(2,fOutputSP);
  PostData(3,fOutputPID);
  PostData(4,fOutputTP);
  PostData(5,fOutputAliEvent);

  return;
}
//___________________________________ hiostograms _______________________________________
void AliAnalysisTaskPLFemto::DefineHistograms(TString whichV0)
{
  // Create histograms
  //Define event related histograms:
  if(!fIsLightweight) {
    fCEvents = new TH1F("fCEvents","conter",11,0,11);
    fCEvents->SetStats(kTRUE);
    fCEvents->GetXaxis()->SetTitle("1");
    fCEvents->GetYaxis()->SetTitle("counts");

    //***************************************************************************
    //Global event histograms
    //***************************************************************************
    TH1F* fVtxX = new TH1F("fVtxX","Primary X vertex",600,-1.,1.);
    TH1F* fVtxY = new TH1F("fVtxY","Primary Y vertex",600,-1.,1.);
    TH1F* fVtxZ = new TH1F("fVtxZ","Primary Z vertex",600,-15.,15.);
    TH1F* fProtonPIDthresholdTPCLow = new TH1F("fProtonPIDthresholdTPCLow","lower momentum threshold which was chosen for PID selection of protons",10,0,10);
    fProtonPIDthresholdTPCLow->SetBinContent(1,fCuts->GetProtonPIDthrPtLow());
    TH1F* fProtonPIDthresholdTPCUp = new TH1F("fProtonPIDthresholdTPCUp","upper momentum threshold which was chosen for PID selection of protons",10,0,10);
    fProtonPIDthresholdTPCUp->SetBinContent(1,fCuts->GetProtonPIDthrPtUp());
    TH1F* fProtonPIDTPCTOFSwitch = new TH1F("fProtonPIDTPCTOFSwitc","Switch between TOF and TPC for proton PID",10,0,10);
    fProtonPIDTPCTOFSwitch->SetBinContent(1,fCuts->GetProtonPIDTOFTPCSwitch());
    TH1F* fLambdaPIDthresholdTPCLow = new TH1F("fLambdaPIDthresholdTPCLow","lower momentum threshold which was chosen for PID selection of Lambdas",10,0,10);
    fLambdaPIDthresholdTPCLow->SetBinContent(1,fCuts->GetV0PIDthrPtLow());
    TH1F* fLambdaPIDthresholdTPCUp = new TH1F("fLambdaPIDthresholdTPCUp","upper momentum threshold which was chosen for PID selection of Lambdas",10,0,10);
    fLambdaPIDthresholdTPCUp->SetBinContent(1,fCuts->GetProtonPIDthrPtUp());
    TH1F* fProtonDCAxyValue = new TH1F("fProtonDCAxyValue","DCAxy value for protons",10,0,10);
    fProtonDCAxyValue->SetBinContent(1,fCuts->GetProtonDCAxyCut());
    TH1F* fProtonDCAzValue = new TH1F("fProtonDCAzValue","DCAxy value for protons",10,0,10);
    fProtonDCAzValue->SetBinContent(1,fCuts->GetProtonDCAzCut());

    TH1F* fTrackselection = new TH1F("fTrackselection","which filterbit was used to select the tracks",10,0,10);
    fTrackselection->SetBinContent(1,fWhichfilterbit);

    //Define single particle histograms event related:
    TH1F* fWhichAnalysis = new TH1F("fWhichAnalysis","Which V0 is correlated, 1:Lambda,2:KOshort",10,0,10);
    TH1F* fWhichV0region = new TH1F("fWhichV0region","Sideband or signal: 1:left sideband,2:signal,3:right sideband",10,0,10);
    TH1F* fOnlineOrOffline = new TH1F("fOnlineOrOffline","Histogram that displays if V0 selection was online or on-the-fly, 1:online, 2:offline",10,0,10);
    TH1F* fMCorExp = new TH1F("fMCorExp","Monte Carlo or experimental data, 1: Exp, 2: MC",10,0,10);
    TH2F* fEventSPDClusterNTracklets = new TH2F("fEventSPDClusterNTracklets","SPDCluster vs Tracklets: Is pile-up seen",200,0,200,200,0,200);
    TH2F* fEventSPDClusterNTrackletsPileUpCut = new TH2F("fEventSPDClusterNTrackletsPileUpCut","SPDCluster vs Tracklets: Is pile-up seen",200,0,200,200,0,200);

    if(fOnlineV0) fOnlineOrOffline->Fill(1);
    else if(!fOnlineV0) fOnlineOrOffline->Fill(3);

    if(whichV0=="Lambda") fWhichAnalysis->Fill(1);
    else if(whichV0=="Kaon") fWhichAnalysis->Fill(3);

    if(fwhichV0region == "left") fWhichV0region->Fill(1);
    else if(fwhichV0region == "signal") fWhichV0region->Fill(2);
    else if(fwhichV0region == "right") fWhichV0region->Fill(3);

    if(fUseMCInfo) fMCorExp->Fill(1);
    else fMCorExp->Fill(2);

    fOutput->Add(fCEvents);
    fOutput->Add(fVtxX);
    fOutput->Add(fVtxY);
    fOutput->Add(fVtxZ);
    fOutput->Add(fProtonPIDthresholdTPCLow);
    fOutput->Add(fProtonPIDthresholdTPCUp);
    fOutput->Add(fProtonPIDTPCTOFSwitch);
    fOutput->Add(fLambdaPIDthresholdTPCLow);
    fOutput->Add(fLambdaPIDthresholdTPCUp);
    fOutput->Add(fWhichAnalysis);
    fOutput->Add(fWhichV0region);
    fOutput->Add(fOnlineOrOffline);
    fOutput->Add(fEventSPDClusterNTracklets);
    fOutput->Add(fEventSPDClusterNTrackletsPileUpCut);
    fOutput->Add(fProtonDCAxyValue);
    fOutput->Add(fProtonDCAzValue);
    fOutput->Add(fTrackselection);
    fOutput->Add(fMCorExp);
  }


  //***************************************************************************
  //***************************************************************************



  //***************************************************************************
  //Single particle histograms
  //***************************************************************************
  if(!fIsLightweight) {
    Float_t invMassLow = 0.;
    Float_t invMassUp = 0.;

    if(whichV0=="Lambda")
      {
        invMassLow = 1.;
        invMassUp = 1.2;
      }
    else if(whichV0=="Kaon")
      {
        invMassLow = 0.4;
        invMassUp = 0.6;
      }

    TH1F* fInvMassLambdawCuts = new TH1F("fInvMassLambdawCuts","Invariant mass of Lambda(p #pi) (GeV/c) with topological cuts",400,invMassLow,invMassUp);
    TH1F* fInvMassLambdawoCuts = new TH1F("fInvMassLambdawoCuts","Invariant mass of Lambda(p #pi) (GeV/c) without topological cuts",400,invMassLow,invMassUp);
    TH1F* fInvMassLambdawCutsAfterSelection = new TH1F("fInvMassLambdawCutsAfterSelection","Invariant mass of Lambda(p #pi) (GeV/c) with topological cuts and after proton selection",400,invMassLow,invMassUp);
    TH1F* fInvMassLambdawCutsAfterMassCut = new TH1F("fInvMassLambdawCutsAfterMassCut","Invariant mass of Lambda(p #pi) (GeV/c) with vertex cuts and Mass cut of 4 MeV/c²",400,invMassLow,invMassUp);

    //TH1F* fLambdaSignalBkg = new TH1F("fLambdaSignalBkg","Count true and fake lambda pairs",10,0,10);

    TH1F* fInvMassLambdaRejected = new TH1F("fInvMassLambdaRejected","Invariant mass of Lambda(p #pi) (GeV/c) which were rejected from the sample",400,invMassLow,invMassUp);
    TH1F* fInvMassLambdaKept = new TH1F("fInvMassLambdaKept","Invariant mass of Lambda(p #pi) (GeV/c) which were kept in the sample",400,invMassLow,invMassUp);

    //Differential Analysis in Pt:
    TH1F* fInvMassLambdawCutsPtBin[fCuts->GetV0PtBinsInvMass()];

    for(Int_t i=0; i<fCuts->GetV0PtBinsInvMass(); i++)
      {
        TString HistName = "fInvMassLambdawCutsPtBin";
        HistName += i;
        fInvMassLambdawCutsPtBin[i] = new TH1F(HistName.Data(),HistName.Data(),200,1.,1.2);
        fOutputSP->Add(fInvMassLambdawCutsPtBin[i]);
      }


    TH1F* fHistTrackFilterMap = new TH1F("fHistTrackFilterMap","Includes different filter maps for all the corresponding global tracks",4000,0,2000);

    TH1F* fInvMassAntiLambdawCuts = new TH1F("fInvMassAntiLambdawCuts","Invariant mass of AntiLambda(a-p a-#pi) (GeV/c) with vertex cuts",400,invMassLow,invMassUp);
    TH1F* fInvMassAntiLambdawoCuts = new TH1F("fInvMassAntiLambdawoCuts","Invariant mass of AntiLambda(a-p a-#pi) (GeV/c) without topological cuts",400,invMassLow,invMassUp);

    TH1F* fInvMassMissIDK0s = new TH1F("fInvMassMissIDK0s","Invariant mass of K0s under the hypothesis we have a lambda",400,0.4,0.6);
    TH1F* fInvMassMissIDK0swCuts = new TH1F("fInvMassMissIDK0swCuts","Invariant mass of K0s after mass cut applied",400,0.4,0.6);


    TH1F* fXiInvMasswoCuts = new TH1F("fXiInvMasswoCuts","Invariant mass of Xi candidate",300,1.2,1.5);
    TH1F* fXiInvMasswCuts = new TH1F("fXiInvMasswCuts","Invariant mass of Xi candidate with topological cuts",300,1.2,1.5);
    TH1F* fXiInvMassLambda = new TH1F("fXiInvMassLambda","Invariant mass of p-pi stemming from Xi",400,invMassLow,invMassUp);

    TH1F* fLambdaDCADaughterTracks = new TH1F("fLambdaDCADaughterTracks","Distance of closest approach of p and #pi track",400,0.,2.);
    TH1F* fLambdaDCAPrimVertex = new TH1F("fLambdaDCAPrimVertex","Distance of closest approach of V0 track to primary vertex",500,0.,10.);
    TH1F* fLambdaDCAPosdaughPrimVertex = new TH1F("fLambdaDCAPosdaughPrimVertex","Distance of closest approach of V0 positive daughter track to primary vertex",500,0.,10.);
    TH1F* fLambdaDCANegdaughPrimVertex = new TH1F("fLambdaDCANegdaughPrimVertex","Distance of closest approach of V0 negative daughter track to primary vertex",500,0.,10.);
    TH1F* fLambdaTransverseRadius = new TH1F("fLambdaTransverseRadius","Transverse distance between primary vertex and V0 decay vertex",400,0,60.);
    TH1F* fLambdaDecaylength = new TH1F("fLambdaDecaylength","Distance between primary and secondary vertex",500,0.,100.);
    TH1F* fLambdaCosPointingAngle = new TH1F("fLambdaCosPointingAngle","Cosinuns of Pointing angle",800,0.8,1.);
    TH1F* fProtonDCAxy = new TH1F("fProtonDCAxy","Distance of closest approach of primary proton to primary vertex in xy",1000,-5.,5.);
    TH1F* fProtonDCAxyCutz = new TH1F("fProtonDCAxyCutz","Distance of closest approach of primary proton to primary vertex in xy with DCAz Cut",1000,-5.,5.);
    TH1F* fProtonDCAz = new TH1F("fProtonDCAz","Distance of closest approach of primary proton to primary vertex in z",1000,-20.,20.);
    TH2F* fProtonDCAxyDCAz = new TH2F("fProtonDCAxyDCAz","Distribution of DCAz vs DCAxy",500,-5,5,1000,-20,20);
    TH1F* fAntiProtonDCAxy = new TH1F("fAntiProtonDCAxy","Distance of closest approach of primary proton to primary vertex in xy",1000,-5.,5.);
    TH1F* fAntiProtonDCAxyCutz = new TH1F("fAntiProtonDCAxyCutz","Distance of closest approach of primary proton to primary vertex in xy with DCAz Cut",1000,-5.,5.);
    TH1F* fAntiProtonDCAz = new TH1F("fAntiProtonDCAz","Distance of closest approach of primary proton to primary vertex in z",1000,-20.,20.);
    TH2F* fAntiProtonDCAxyDCAz = new TH2F("fAntiProtonDCAxyDCAz","Distribution of DCAz vs DCAxy",500,-5,5,1000,-20,20);

    TH2F* fProtonDCAxyDCAzMCPtBin[5][fCuts->GetPtBinsDCA()];
    TH2F* fProtonDCAxyDCAzMCPtBinLambda[fCuts->GetPtBinsDCA()];
    TH2F* fProtonDCAxyDCAzMCPtBinSigma[fCuts->GetPtBinsDCA()];
    TH2F* fProtonDCAxyDCAzPt[fCuts->GetPtBinsDCA()];

    for(int ProtonCases = 0; ProtonCases < 5; ProtonCases++)//primary, secondary, ....
      {
        for(int ptbins = 0; ptbins < fCuts->GetPtBinsDCA(); ptbins++)
      {
        TString HistName = "fProtonDCAxyDCAzMCCase";
        HistName += ProtonCases;
      HistName += "PtBin";
      HistName += ptbins;

      fProtonDCAxyDCAzMCPtBin[ProtonCases][ptbins] = new TH2F(HistName.Data(),"DCAz vs DCAxy for different pt Bins",500,-5,5,500,-5,5);
      if(fUseMCInfo) fOutputSP->Add(fProtonDCAxyDCAzMCPtBin[ProtonCases][ptbins]);

      if(ProtonCases == 0)
        {
          HistName = "fProtonDCAxyDCAzMCPtBinLambda";
          HistName += "PtBin";
          HistName += ptbins;

          fProtonDCAxyDCAzMCPtBinLambda[ptbins] = new TH2F(HistName.Data(),"DCAz vs DCAxy for different pt Bins for Lambda hyperons",500,-5,5,500,-5,5);
          if(fUseMCInfo) fOutputSP->Add(fProtonDCAxyDCAzMCPtBinLambda[ptbins]);

          HistName = "fProtonDCAxyDCAzMCPtBinSigma";
          HistName += "PtBin";
          HistName += ptbins;

          fProtonDCAxyDCAzMCPtBinSigma[ptbins] = new TH2F(HistName.Data(),"DCAz vs DCAxy for different pt Bins for Sigma+ hyperons",500,-5,5,500,-5,5);
          if(fUseMCInfo) fOutputSP->Add(fProtonDCAxyDCAzMCPtBinSigma[ptbins]);
        }

      if(ProtonCases == 0)
        {
          HistName = "fProtonDCAxyDCAzPt";
          HistName += ptbins;
          fProtonDCAxyDCAzPt[ptbins] = new TH2F(HistName.Data(),"DCAz vs DCAxy for different pt Bins",500,-5,5,500,-5,5);
          fOutputSP->Add(fProtonDCAxyDCAzPt[ptbins]);
        }
    }
      }

    TH1F* fLambdaCPAPtBin[fCuts->GetV0PtBinsCPA()];
    TH1F* fLambdaCPAAllPtBin[fCuts->GetV0PtBinsCPA()];
    TH1F* fLambdaCPAPrimaryPtBin[fCuts->GetV0PtBinsCPA()];
    TH1F* fLambdaCPASecondaryPtBin[fCuts->GetV0PtBinsCPA()];
    TH1F* fLambdaCPAMaterialPtBin[fCuts->GetV0PtBinsCPA()];
    TH1F* fLambdaCPABkgPtBin[fCuts->GetV0PtBinsCPA()];

    for(int ptbins = 0; ptbins < fCuts->GetV0PtBinsCPA(); ptbins++)
      {
        TString HistName = "fLambdaCPAPtBin";
        HistName += ptbins;
        fLambdaCPAPtBin[ptbins] = new TH1F(HistName.Data(),"Cosine of pointing angle as function of pt for all V0s",100,0.99,1);
        fOutputSP->Add(fLambdaCPAPtBin[ptbins]);

        if(fUseMCInfo)
      {
      HistName = "fLambdaCPAAllPtBin";
      HistName += ptbins;
      fLambdaCPAAllPtBin[ptbins] = new TH1F(HistName.Data(),"Cosine of pointing angle as function of pt for all real V0s",100,0.99,1);
      fOutputSP->Add(fLambdaCPAAllPtBin[ptbins]);

      HistName = "fLambdaCPAPrimaryPtBin";
      HistName += ptbins;
      fLambdaCPAPrimaryPtBin[ptbins] = new TH1F(HistName.Data(),"Cosine of pointing angle as function of pt for primary real V0s",100,0.99,1);
      fOutputSP->Add(fLambdaCPAPrimaryPtBin[ptbins]);

      HistName = "fLambdaCPASecondaryPtBin";
      HistName += ptbins;
      fLambdaCPASecondaryPtBin[ptbins] = new TH1F(HistName.Data(),"Cosine of pointing angle as function of pt for secondary real V0s",100,0.99,1);
      fOutputSP->Add(fLambdaCPASecondaryPtBin[ptbins]);

      HistName = "fLambdaCPAMaterialPtBin";
      HistName += ptbins;
      fLambdaCPAMaterialPtBin[ptbins] = new TH1F(HistName.Data(),"Cosine of pointing angle as function of pt for material real V0s",100,0.99,1);
      fOutputSP->Add(fLambdaCPAMaterialPtBin[ptbins]);

      HistName = "fLambdaCPABkgPtBin";
      HistName += ptbins;
      fLambdaCPABkgPtBin[ptbins] = new TH1F(HistName.Data(),"Cosine of pointing angle as function of pt for background V0s",100,0.99,1);
      fOutputSP->Add(fLambdaCPABkgPtBin[ptbins]);
    }
    }

    TH1F* fLambdaPt = new TH1F("fLambdaPt","Transverse momentum spectrum of V0s",50,0.,10.);
    TH1F* fLambdaPosDaughPt = new TH1F("fLambdaPosDaughPt","Transverse momentum of positive daughter of V0",500,0.,10.);
    TH1F* fLambdaNegDaughPt = new TH1F("fLambdaNegDaughPt","Transverse momentum of negative daughter of V0",500,0.,10.);
    TH1F* fLambdaPhi = new TH1F("fLambdaPhi","Phi angle spectrum of V0s",300,-1,2.*TMath::Pi());
    TH1F* fLambdaEta = new TH1F("fLambdaEta","Eta spectrum of V0s",400,-2,2);
    TH1F* fAntiLambdaPt = new TH1F("fAntiLambdaPt","Transverse momentum spectrum of Anti-V0s",50,0.,10.);
    TH1F* fProtonPt = new TH1F("fProtonPt","Transverse momentum spectrum of protons",400,0.,10.);
    TH1F* fProtonPhi = new TH1F("fProtonPhi","Phi angle spectrum of protons",300,-1,2.*TMath::Pi());
    TH1F* fProtonEta = new TH1F("fProtonEta","Pseudorapidity spectrum of protons",400,-2,2);

    TH1F* fAntiProtonPt = new TH1F("fAntiProtonPt","Transverse momentum spectrum of Anti-protons",50,0.,10.);

    if(fUseMCInfo) {
      TH2F* fLambdaV0vertexXMC = new TH2F("fLambdaV0vertexXMC","Vertex X of V0 decay vertex",800,-300,300,100,0,10);
      TH2F* fLambdaV0vertexYMC = new TH2F("fLambdaV0vertexYMC","Vertex Y of V0 decay vertex",800,-300,300,100,0,10);
      TH2F* fLambdaV0vertexZMC = new TH2F("fLambdaV0vertexZMC","Vertex Z of V0 decay vertex",800,-300,300,100,0,10);
      TH2F* fLambdabkgvertexXMC = new TH2F("fLambdabkgvertexXMC","Vertex X of V0 decay vertex",800,-300,300,100,0,10);
      TH2F* fLambdabkgvertexYMC = new TH2F("fLambdabkgvertexYMC","Vertex Y of V0 decay vertex",800,-300,300,100,0,10);
      TH2F* fLambdabkgvertexZMC = new TH2F("fLambdabkgvertexZMC","Vertex Z of V0 decay vertex",800,-300,300,100,0,10);
      TH2F* fDCAV0PosdaughMC = new TH2F("fDCAV0PosdaughMC","DCA of positive daughter track from V0 to primary vertex in MC as a function of V0 pt",400,0,2,100,0,10);
      TH2F* fDCAbkgPosdaughMC = new TH2F("fDCAbkgPosdaughMC","DCA of positive daughter track from bkg to primary vertex in MC as a function of bkg pt",400,0,2,100,0,10);
      TH2F* fDCAV0NegdaughMC = new TH2F("fDCAV0NegdaughMC","DCA of negative daughter track from V0 to primary vertex in MC as a function of V0 pt",400,0,2,100,0,10);
      TH2F* fDCAbkgNegdaughMC = new TH2F("fDCAbkgNegdaughMC","DCA of negative daughter track from background to primary vertex in MC as a function of bkg pt",400,0,2,100,0,10);
      TH2F* fTrRadiusV0MC = new TH2F("fTrRadiusV0MC","Transverse distance of V0 between primary and secondary vertex MC as a function of V0 pt",3000,0,300,100,0,10);
      TH2F* fTrRadiusbkgMC = new TH2F("fTrRadiusbkgMC","Transverse distance of bkg between primary and secondary vertex MC as a function of bkg pt",3000,0,300,100,0,10);
      TH2F* fDecayLengthV0MC = new TH2F("fDecayLengthV0MC","Transverse distance of V0 between primary and secondary vertex MC as a function of V0 pt",600,0,60,100,0,10);
      TH2F* fDecayLengthbkgMC = new TH2F("fDecayLengthbkgMC","Transverse distance of bkg between primary and secondary vertex MC as a function of bkg pt",600,0,60,100,0,10);
      TH2F* fPoinAngleV0MC = new TH2F("fPoinAngleV0MC","Pointing angle of V0 as a function of V0 pt",500,0.97,1,100,0,10);
      TH2F* fPoinAnglebkgMC = new TH2F("fPoinAnglebkgMC","Pointing angle of bkg as a function of bkg pt",500,0.97,1,100,0,10);
      fOutputSP->Add(fLambdaV0vertexXMC);
      fOutputSP->Add(fLambdabkgvertexXMC);
      fOutputSP->Add(fLambdaV0vertexYMC);
      fOutputSP->Add(fLambdabkgvertexYMC);
      fOutputSP->Add(fLambdaV0vertexZMC);
      fOutputSP->Add(fLambdabkgvertexZMC);
      fOutputSP->Add(fDCAV0PosdaughMC);
      fOutputSP->Add(fDCAbkgPosdaughMC);
      fOutputSP->Add(fDCAV0NegdaughMC);
      fOutputSP->Add(fDCAbkgNegdaughMC);
      fOutputSP->Add(fTrRadiusV0MC);
      fOutputSP->Add(fTrRadiusbkgMC);
      fOutputSP->Add(fPoinAngleV0MC);
      fOutputSP->Add(fPoinAnglebkgMC);
      fOutputSP->Add(fDecayLengthV0MC);
      fOutputSP->Add(fDecayLengthbkgMC);
    }

    TH1F* fNBinsMultmixing = new TH1F("fNBinsMultmixing","Bins in multiplicity that are used for event mixing, which bin is often occupied etc.",kMultiplicityBins,0,kMultiplicityBins);
    TH1F* fNBinsVertexmixing = new TH1F("fNBinsVertexmixing","Bins in z-Vertex that are used for event mixing, which bin is often occupied etc.",kZVertexBins,0,kZVertexBins);
    TH1F* fNTracklets = new TH1F("fNTracklets","Number of Tracklets",200,0,200);
    TH1F* fRefMultiplicity08 = new TH1F("fRefMultiplicity08","Reference multiplicity",200,0,200);
    TH2F* fNTrackletsRefMultiplicity = new TH2F("fNTrackletsRefMultiplicity","Correlation between tracklets and refmultiplicity",200,0,200,200,0,200);
    TH1F* fNTPCTracks = new TH1F("fNTPCTracks","Number of TPC tracks",200,0,200);
    TH1F* fNProtonsPerevent = new TH1F("fNProtonsPerevent","Number of protons in an event",200,1,100);
    TH1F* fNLambdasPerevent = new TH1F("fNLambdasPerevent","Number of V0s in an event",20,1,10);
    TH1F* fNProtonsTot = new TH1F("fNProtonsTot","Total number of protons",10,1,10);
    TH1F* fNAntiProtonsTot = new TH1F("fNAntiProtonsTot","Total number of anti-protons",10,1,10);
    TH1F* fNLambdasTot = new TH1F("fNLambdasTot","Total number of V0s under certain conditions",10,1,10);
    TH1F* fNAntiLambdasTot = new TH1F("fNAntiLambdasTot","Total number of Anti-V0s under certain conditions",10,1,10);
    TH1F* fNV0TrackSharing = new TH1F("fNV0TrackSharing","Number of V0s which shared tracks in an event",20,0,10);
    TH1F* fNAntiV0TrackSharing = new TH1F("fNAntiV0TrackSharing","Number of Anti-V0s which shared tracks in an event",20,0,10);
    TH1F* fNV0protonSharedTracks = new TH1F("fNV0protonSharedTracks","How often share a V0 and a proton the same track?",10,0,10);
    TH1F* fNAntiV0protonSharedTracks = new TH1F("fNAntiV0protonSharedTracks","How often share a Anti-V0 and a proton the same track?",10,0,10);
    TH1F* fNAntiV0AntiprotonSharedTracks = new TH1F("fNAntiV0AntiprotonSharedTracks","How often share a Anti-V0 and a Anti-proton the same track?",10,0,10);
    TH1F* fNXimTrackSharing = new TH1F("fNXimTrackSharing","Number of Xim which shared tracks in an event",20,0,10);

    if(fUseMCInfo) {
      TH1F* fFeeddownProton = new TH1F("fFeeddownProton","1: total number of protons, 2: Primary, 3: from weak decays, 4: from material",10,0,10);
      TH1F* fFeeddownProtonFromWeakPDG = new TH1F("fFeeddownProtonFromWeakPDG","Protons stemming from weak decays",213,3121,3334);
      TH1F* fFeeddownV0 = new TH1F("fFeeddownV0","1: total number of protons, 2: Primary, 3: from weak decays, 4: from material",10,0,10);
      TH1F* fFeeddownV0FromWeakPDG = new TH1F("fFeeddownV0FromWeakPDG","V0s stemming from weak decays",213,3121,3334);
      fOutputSP->Add(fFeeddownProton);
      fOutputSP->Add(fFeeddownProtonFromWeakPDG);
      fOutputSP->Add(fFeeddownV0);
      fOutputSP->Add(fFeeddownV0FromWeakPDG);
    }

    TH1F* fFindableClusterProton = new TH1F("fFindableClusterProton","Number of findable cluster for protons",200,0,200);
    TH1F* fNCrossedRowsProton = new TH1F("fNCrossedRowsProton","Number of crossed rows for protons",200,0,200);
    TH1F* fRatioFindableCrossedProton = new TH1F("fRatioFindableCrossedProton","Number of crossed rows over findable cluster for protons",200,0,1.5);
    TH1F* fFindableClusterAntiProton = new TH1F("fFindableClusterAntiProton","Number of findable cluster for protons",200,0,200);
    TH1F* fNCrossedRowsAntiProton = new TH1F("fNCrossedRowsAntiProton","Number of crossed rows for protons",200,0,200);
    TH1F* fRatioFindableCrossedAntiProton = new TH1F("fRatioFindableCrossedAntiProton","Number of crossed rows over findable cluster for protons",200,0,1.5);

    fOutputSP->Add(fHistTrackFilterMap);
    fOutputSP->Add(fNBinsMultmixing);
    fOutputSP->Add(fNBinsVertexmixing);
    fOutputSP->Add(fNTracklets);
    fOutputSP->Add(fRefMultiplicity08);
    fOutputSP->Add(fNTrackletsRefMultiplicity);
    fOutputSP->Add(fNTPCTracks);
    fOutputSP->Add(fInvMassLambdawCuts);
    fOutputSP->Add(fInvMassLambdawoCuts);
    fOutputSP->Add(fInvMassLambdawCutsAfterSelection);
    fOutputSP->Add(fInvMassLambdaRejected);
    fOutputSP->Add(fInvMassLambdaKept);

    fOutputSP->Add(fInvMassLambdawCutsAfterMassCut);
    fOutputSP->Add(fInvMassAntiLambdawCuts);
    fOutputSP->Add(fInvMassAntiLambdawoCuts);
    fOutputSP->Add(fInvMassMissIDK0s);
    fOutputSP->Add(fInvMassMissIDK0swCuts);

    fOutputSP->Add(fXiInvMasswoCuts);
    fOutputSP->Add(fXiInvMasswCuts);
    fOutputSP->Add(fXiInvMassLambda);

    fOutputSP->Add(fLambdaPt);
    fOutputSP->Add(fLambdaPosDaughPt);
    fOutputSP->Add(fLambdaNegDaughPt);
    fOutputSP->Add(fLambdaPhi);
    fOutputSP->Add(fLambdaEta);
    fOutputSP->Add(fAntiLambdaPt);
    fOutputSP->Add(fProtonPt);
    fOutputSP->Add(fProtonPhi);
    fOutputSP->Add(fProtonEta);
    fOutputSP->Add(fAntiProtonPt);
    fOutputSP->Add(fNLambdasTot);
    fOutputSP->Add(fNAntiLambdasTot);
    fOutputSP->Add(fNLambdasPerevent);
    fOutputSP->Add(fNProtonsTot);
    fOutputSP->Add(fNAntiProtonsTot);
    fOutputSP->Add(fNProtonsPerevent);
    fOutputSP->Add(fNV0TrackSharing);
    fOutputSP->Add(fNAntiV0TrackSharing);
    fOutputSP->Add(fNV0protonSharedTracks);
    fOutputSP->Add(fNAntiV0protonSharedTracks);
    fOutputSP->Add(fNAntiV0AntiprotonSharedTracks);
    fOutputSP->Add(fNXimTrackSharing);

    fOutputSP->Add(fLambdaDCADaughterTracks);
    fOutputSP->Add(fLambdaDCAPrimVertex);
    fOutputSP->Add(fLambdaDCAPosdaughPrimVertex);
    fOutputSP->Add(fLambdaDCANegdaughPrimVertex);
    fOutputSP->Add(fLambdaDecaylength);
    fOutputSP->Add(fLambdaCosPointingAngle);
    fOutputSP->Add(fLambdaTransverseRadius);
    fOutputSP->Add(fProtonDCAxy);
    fOutputSP->Add(fProtonDCAxyCutz);
    fOutputSP->Add(fProtonDCAz);
    fOutputSP->Add(fProtonDCAxyDCAz);
    fOutputSP->Add(fAntiProtonDCAxy);
    fOutputSP->Add(fAntiProtonDCAxyCutz);
    fOutputSP->Add(fAntiProtonDCAz);
    fOutputSP->Add(fAntiProtonDCAxyDCAz);

    fOutputSP->Add(fFindableClusterProton);
    fOutputSP->Add(fNCrossedRowsProton);
    fOutputSP->Add(fRatioFindableCrossedProton);
    fOutputSP->Add(fFindableClusterAntiProton);
    fOutputSP->Add(fNCrossedRowsAntiProton);
    fOutputSP->Add(fRatioFindableCrossedAntiProton);
  }


  //Define two particle histograms:
  //p-lambda
  TH1F* fProtonLambdaRelK = new TH1F("fProtonLambdaRelK","Relative momentum of V0-p RE",150,0.,3.); //first choice was 300 bins
  TH1F* fProtonLambdaRelKME = new TH1F("fProtonLambdaRelKME","Relative momentum of V0-p ME",150,0.,3.);
  if(!fIsLightweight) {
    TH1F* fProtonLambdaRelKMulti[13];
    TH1F* fProtonLambdaRelKMEMulti[13];
    for(int i=0; i<13; ++i) {
        fProtonLambdaRelKMulti[i] = new TH1F(Form("fProtonLambdaRelKMulti_%i", i),"Relative momentum of V0-p RE",150,0.,3.); //first choice was 300 bins
        fProtonLambdaRelKMEMulti[i] = new TH1F(Form("fProtonLambdaRelKMEMulti_%i", i),"Relative momentum of V0-p ME",150,0.,3.);
        fOutputTP->Add(fProtonLambdaRelKMulti[i]);
        fOutputTP->Add(fProtonLambdaRelKMEMulti[i]);
      }
  }

  if(!fIsLightweight) {
    TH1F* fProtonLambdaPosDaughterRelK = new TH1F("fProtonLambdaPosDaughterRelK","Relative momentum of daughter proton and primary proton RE",150,0.,3.); //first choice was 300 bins
    TH1F* fProtonLambdaPosDaughterRelKME = new TH1F("fProtonLambdaPosDaughterRelKME","Relative momentum of daughter proton and primary proton ME",150,0.,3.);
    TH1F* fAntiProtonAntiLambdaPosDaughterRelK = new TH1F("fAntiProtonAntiLambdaPosDaughterRelK","Relative momentum of daughter anti-proton and primary anti-proton RE",150,0.,3.); //first choice was 300 bins
    TH1F* fAntiProtonAntiLambdaPosDaughterRelKME = new TH1F("fAntiProtonAntiLambdaPosDaughterRelKME","Relative momentum of daughter anti-proton and primary anti-proton ME",150,0.,3.);
    TH1F* fProtonLambdaDiff = new TH1F("fProtonLambdaDiff","diff of qinv relk",1000,-0.1,0.1);
    TH1F* fProtonLambdaAvgSeparationPP = new TH1F("fProtonLambdaAvgSeparationPP","Relative distances of pp_daugh RE",2000,0.,500.);
    TH1F* fProtonLambdaAvgSeparationPPME = new TH1F("fProtonLambdaAvgSeparationPPME","Relative distances of pp_daugh ME",2000,0.,500.);
    TH1F* fProtonLambdaAvgSeparationPPi = new TH1F("fProtonLambdaAvgSeparationPPi","Relative distances of ppi_daugh RE",2000,0.,500.);
    TH1F* fProtonLambdaAvgSeparationPPiME = new TH1F("fProtonLambdaAvgSeparationPPiME","Relative distances of ppi_daugh ME",2000,0.,500.);
    fOutputTP->Add(fProtonLambdaPosDaughterRelK);
    fOutputTP->Add(fProtonLambdaPosDaughterRelKME);
    fOutputTP->Add(fAntiProtonAntiLambdaPosDaughterRelK);
    fOutputTP->Add(fAntiProtonAntiLambdaPosDaughterRelKME);
    fOutputTP->Add(fProtonLambdaDiff);
    fOutputTP->Add(fProtonLambdaAvgSeparationPP);
    fOutputTP->Add(fProtonLambdaAvgSeparationPPME);
    fOutputTP->Add(fProtonLambdaAvgSeparationPPi);
    fOutputTP->Add(fProtonLambdaAvgSeparationPPiME);
  }

  //Antip-Antilambda
  TH1F* fAntiProtonAntiLambdaRelK = new TH1F("fAntiProtonAntiLambdaRelK","Relative momentum of Anti-V0-Antip RE",150,0.,3.); //first choice was 300 bins
  TH1F* fAntiProtonAntiLambdaRelKME = new TH1F("fAntiProtonAntiLambdaRelKME","Relative momentum of Anti-V0-Antip ME",150,0.,3.);

  if(!fIsLightweight) {
    TH1F* fAntiProtonAntiLambdaRelKMulti[13];
    TH1F* fAntiProtonAntiLambdaRelKMEMulti[13];
    for(int i=0; i<13; ++i) {
        fAntiProtonAntiLambdaRelKMulti[i] = new TH1F(Form("fAntiProtonAntiLambdaRelKMulti_%i", i), Form("Relative momentum of AntiV0-Antip RE %i", i),150,0.,3.); //first choice was 300 bins
        fAntiProtonAntiLambdaRelKMEMulti[i] = new TH1F(Form("fAntiProtonAntiLambdaRelKMEMulti_%i", i), Form("Relative momentum of AntiV0-Antip ME %i", i),150,0.,3.);
        fOutputTP->Add(fAntiProtonAntiLambdaRelKMulti[i]);
        fOutputTP->Add(fAntiProtonAntiLambdaRelKMEMulti[i]);
      }
  }

  if(!fIsLightweight) {
    TH1F* fAntiProtonAntiLambdaAvgSeparationPP = new TH1F("fAntiProtonAntiLambdaAvgSeparationPP","Relative distances of pp_daugh RE",2000,0.,500.);
    TH1F* fAntiProtonAntiLambdaAvgSeparationPPME = new TH1F("fAntiProtonAntiLambdaAvgSeparationPPME","Relative distances of pp_daugh ME",2000,0.,500.);
    TH1F* fAntiProtonAntiLambdaAvgSeparationPPi = new TH1F("fAntiProtonAntiLambdaAvgSeparationPPi","Relative distances of ppi_daugh RE",2000,0.,500.);
    TH1F* fAntiProtonAntiLambdaAvgSeparationPPiME = new TH1F("fAntiProtonAntiLambdaAvgSeparationPPiME","Relative distances of ppi_daugh ME",2000,0.,500.);
    fOutputTP->Add(fAntiProtonAntiLambdaAvgSeparationPP);
    fOutputTP->Add(fAntiProtonAntiLambdaAvgSeparationPPME);
    fOutputTP->Add(fAntiProtonAntiLambdaAvgSeparationPPi);
    fOutputTP->Add(fAntiProtonAntiLambdaAvgSeparationPPiME);
  }

  TH1F* fProtonLambdaMT = new TH1F("fProtonLambdaMT","Transverse mass of pL",1000,0.,10.);
  TH1F* fProtonLambdaMTrelKcut = new TH1F("fProtonLambdaMTrelKcut","Transverse mass of pL",1000,0.,10.);

  //p-Antilambda
  TH1F* fAntiProtonLambdaRelK = new TH1F("fAntiProtonLambdaRelK","Relative momentum of V0-Anti-p RE",150,0.,3.); //first choice was 300 bins
  TH1F* fAntiProtonLambdaRelKME = new TH1F("fAntiProtonLambdaRelKME","Relative momentum of V0-Anti-p ME",150,0.,3.);
  TH1F* fProtonAntiLambdaRelK = new TH1F("fProtonAntiLambdaRelK","Relative momentum of Anti-V0-p RE",150,0.,3.);
  TH1F* fProtonAntiLambdaRelKME = new TH1F("fProtonAntiLambdaRelKME","Relative momentum of Anti-V0-p ME",150,0.,3.);

  //p-Xi
  TH1F* fProtonXiRelK = new TH1F("fProtonXiRelK","Relative momentum of Xi-p RE",150,0.,3.);
  TH1F* fProtonXiRelKME = new TH1F("fProtonXiRelKME","Relative momentum of Xi-p ME",150,0.,3.);
  if(!fIsLightweight) {
    TH1F* fProtonXiRelKMulti[13];
    TH1F* fProtonXiRelKMEMulti[13];
    for(int i=0; i<13; ++i) {
        fProtonXiRelKMulti[i] = new TH1F(Form("fProtonXiRelKMulti_%i", i), Form("Relative momentum of Xi-p RE %i", i),150,0.,3.); //first choice was 300 bins
        fProtonXiRelKMEMulti[i] = new TH1F(Form("fProtonXiRelKMEMulti_%i", i), Form("Relative momentum of Xi-p ME %i", i),150,0.,3.);
        fOutputTP->Add(fProtonXiRelKMulti[i]);
        fOutputTP->Add(fProtonXiRelKMEMulti[i]);
    }
  }

  TH1F* fAntiProtonXiRelK = new TH1F("fAntiProtonXiRelK","Relative momentum of Xi-Antip RE",150,0.,3.);
  TH1F* fAntiProtonXiRelKME = new TH1F("fAntiProtonXiRelKME","Relative momentum of Xi-Antip ME",150,0.,3.);
  if(!fIsLightweight) {
    TH1F* fAntiProtonXiRelKMulti[13];
    TH1F* fAntiProtonXiRelKMEMulti[13];
    for(int i=0; i<13; ++i) {
        fAntiProtonXiRelKMulti[i] = new TH1F(Form("fAntiProtonXiRelKMulti_%i", i), Form("Relative momentum of Xi-Antip RE %i", i),150,0.,3.); //first choice was 300 bins
        fAntiProtonXiRelKMEMulti[i] = new TH1F(Form("fAntiProtonXiRelKMEMulti_%i", i), Form("Relative momentum of Xi-Antip ME %i", i),150,0.,3.);
        fOutputTP->Add(fAntiProtonXiRelKMulti[i]);
        fOutputTP->Add(fAntiProtonXiRelKMEMulti[i]);
    }
  }

  //p-p
  TH1F* fProtonProtonRelK = new TH1F("fProtonProtonRelK","Relative momentum of pp RE",750,0.,3.);
  TH1F* fProtonProtonRelKME = new TH1F("fProtonProtonRelKME","Relative momentum of pp ME",750,0.,3.);
  if(!fIsLightweight) {
    TH1F* fProtonProtonRelKMulti[13];
    TH1F* fProtonProtonRelKMEMulti[13];
    for(int i=0; i<13; ++i) {
        fProtonProtonRelKMulti[i] = new TH1F(Form("fProtonProtonRelKMulti_%i", i), Form("Relative momentum of pp RE %i", i),750,0.,3.); //first choice was 300 bins
        fProtonProtonRelKMEMulti[i] = new TH1F(Form("fProtonProtonRelKMEMulti_%i", i), Form("Relative momentum of pp ME %i", i),750,0.,3.);
        fOutputTP->Add(fProtonProtonRelKMulti[i]);
        fOutputTP->Add(fProtonProtonRelKMEMulti[i]);
    }
  }

  if(!fIsLightweight) {
    TH1F* fProtonProtonRelKw2TC = new TH1F("fProtonProtonRelKw2TC","Relative momentum of pp RE with two track cuts",750,0.,3.);
    TH1F* fProtonProtonRelKw2TCME = new TH1F("fProtonProtonRelKw2TCME","Relative momentum of pp with  two track cuts ME",750,0.,3.);
    fOutputTP->Add(fProtonProtonRelKw2TC);
    fOutputTP->Add(fProtonProtonRelKw2TCME);

    TH1F* fProtonProtonAvgSeparation = new TH1F("fProtonProtonAvgSeparation","Relative distances of pp RE",2000,0.,500.);
    TH1F* fProtonProtonAvgSeparationME = new TH1F("fProtonProtonAvgSeparationME","Relative distances of pp ME",2000,0.,500.);
    TH2F* fProtonProtonAvgSeparationVsRelK = new TH2F("fProtonProtonAvgSeparationVsRelK","Average separation vs. relK",500,0,50,500,0,3);
    TH2F* fProtonProtonAvgSeparationVsRelKME = new TH2F("fProtonProtonAvgSeparationVsRelKME","Average separation vs. relK",500,0,50,500,0,3);
    fOutputTP->Add(fProtonProtonAvgSeparation);
    fOutputTP->Add(fProtonProtonAvgSeparationME);
    fOutputTP->Add(fProtonProtonAvgSeparationVsRelK);
    fOutputTP->Add(fProtonProtonAvgSeparationVsRelKME);

    //Check also EMCIC histos:
    TH1F* fProtonProtonRelKEnSumME = new TH1F("fProtonProtonRelKEnSumME","Relative momentum of pp weighted with total energy ME",750,0.,3.);
    TH1F* fProtonProtonRelKEnMultME = new TH1F("fProtonProtonRelKEnMultME","Relative momentum of pp weighted with multiplied energy ME",750,0.,3.);
    TH1F* fProtonProtonRelKPzMultME = new TH1F("fProtonProtonRelKPzMultME","Relative momentum of pp weighted with multiplied z-components of momenta ME",750,0.,3.);
    TH1F* fProtonProtonRelKPtMultME = new TH1F("fProtonProtonRelKPtMultME","Relative momentum of pp weighted with multiplied pt-components of momenta ME",750,0.,3.);
    fOutputTP->Add(fProtonProtonRelKEnSumME);
    fOutputTP->Add(fProtonProtonRelKEnMultME);
    fOutputTP->Add(fProtonProtonRelKPzMultME);
    fOutputTP->Add(fProtonProtonRelKPtMultME);

    TH1F* fAntiProtonAntiProtonRelKw2TC = new TH1F("fAntiProtonAntiProtonRelKw2TC","Relative momentum of antip-antip RE with two track cuts",750,0.,3.);
    TH1F* fAntiProtonAntiProtonRelKw2TCME = new TH1F("fAntiProtonAntiProtonRelKw2TCME","Relative momentum of antip-antip ME with two track cuts",750,0.,3.);
    TH1F* fAntiProtonAntiProtonAvgSeparation = new TH1F("fAntiProtonAntiProtonAvgSeparation","Relative distances of pp RE",2000,0.,500.);
    TH1F* fAntiProtonAntiProtonAvgSeparationME = new TH1F("fAntiProtonAntiProtonAvgSeparationME","Relative distances of pp ME",2000,0.,500.);
    fOutputTP->Add(fAntiProtonAntiProtonRelKw2TC);
    fOutputTP->Add(fAntiProtonAntiProtonRelKw2TCME);
    fOutputTP->Add(fAntiProtonAntiProtonAvgSeparation);
    fOutputTP->Add(fAntiProtonAntiProtonAvgSeparationME);
  }

  TH1F* fProtonProtonMT = new TH1F("fProtonProtonMT","Transverse mass of pp",1000,0.,10.);
  TH1F* fProtonProtonMTrelKcut = new TH1F("fProtonProtonMTrelKcut","Transverse mass of pp with relative momentum cut",1000,0.,10.);

  //Antip-Antip
  TH1F* fAntiProtonAntiProtonRelK = new TH1F("fAntiProtonAntiProtonRelK","Relative momentum of antip-antip RE",750,0.,3.);
  TH1F* fAntiProtonAntiProtonRelKME = new TH1F("fAntiProtonAntiProtonRelKME","Relative momentum of antip-antip ME",750,0.,3.);
  if(!fIsLightweight) {
    TH1F* fAntiProtonAntiProtonRelKMulti[13];
    TH1F* fAntiProtonAntiProtonRelKMEMulti[13];
    for(int i=0; i<13; ++i) {
        fAntiProtonAntiProtonRelKMulti[i] = new TH1F(Form("fAntiProtonAntiProtonRelKMulti_%i", i), Form("Relative momentum of AntipAntip RE %i", i),750,0.,3.); //first choice was 300 bins
        fAntiProtonAntiProtonRelKMEMulti[i] = new TH1F(Form("fAntiProtonAntiProtonRelKMEMulti_%i", i), Form("Relative momentum of AntipAntip ME %i", i),750,0.,3.);
        fOutputTP->Add(fAntiProtonAntiProtonRelKMulti[i]);
        fOutputTP->Add(fAntiProtonAntiProtonRelKMEMulti[i]);
    }
  }


  //Antip-p
  TH1F* fAntiProtonProtonRelK = new TH1F("fAntiProtonProtonRelK","Relative momentum of antip-p RE",300,0.,3.);
  TH1F* fAntiProtonProtonRelKME = new TH1F("fAntiProtonProtonRelKME","Relative momentum of p-antip ME",300,0.,3.);

  //Lambda-Lambda
  TH1F* fLambdaLambdaRelK = new TH1F("fLambdaLambdaRelK","Relative momentum of V0-V0 RE",150,0.,3.);
  TH1F* fLambdaLambdaRelKME = new TH1F("fLambdaLambdaRelKME","Relative momentum of V0-V0 ME",150,0.,3.);
  if(!fIsLightweight) {
    TH1F* fLambdaLambdaRelKMulti[13];
    TH1F* fLambdaLambdaRelKMEMulti[13];
    for(int i=0; i<13; ++i) {
        fLambdaLambdaRelKMulti[i] = new TH1F(Form("fLambdaLambdaRelKMulti_%i", i), Form("Relative momentum of V0-V0 RE %i", i),150,0.,3.); //first choice was 300 bins
        fLambdaLambdaRelKMEMulti[i] = new TH1F(Form("fLambdaLambdaRelKMEMulti_%i", i), Form("Relative momentum of V0-V0 ME %i", i),150,0.,3.);
        fOutputTP->Add(fLambdaLambdaRelKMulti[i]);
        fOutputTP->Add(fLambdaLambdaRelKMEMulti[i]);
    }
  }

  if(!fIsLightweight) {
    TH1F* fLambdaLambdaAvgSeparationPP = new TH1F("fLambdaLambdaAvgSeparationPP","Relative distances of p_daugh-p_daugh RE",1000,0.,500.);
    TH1F* fLambdaLambdaAvgSeparationPPME = new TH1F("fLambdaLambdaAvgSeparationPPME","Relative distances of p_daugh-p_daugh ME",1000,0.,500.);
    TH1F* fLambdaLambdaAvgSeparationP1Pi2 = new TH1F("fLambdaLambdaAvgSeparationP1Pi2","Relative distances of ppi_daugh RE",1000,0.,500.);
    TH1F* fLambdaLambdaAvgSeparationP1Pi2ME = new TH1F("fLambdaLambdaAvgSeparationP1Pi2ME","Relative distances of ppi_daugh ME",1000,0.,500.);
    TH1F* fLambdaLambdaAvgSeparationP2Pi1 = new TH1F("fLambdaLambdaAvgSeparationP2Pi1","Relative distances of ppi_daugh RE",1000,0.,500.);
    TH1F* fLambdaLambdaAvgSeparationP2Pi1ME = new TH1F("fLambdaLambdaAvgSeparationP2Pi1ME","Relative distances of ppi_daugh ME",1000,0.,500.);
    TH1F* fLambdaLambdaAvgSeparationPiPi = new TH1F("fLambdaLambdaAvgSeparationPiPi","Relative distances of p_daugh-p_daugh RE",1000,0.,500.);
    TH1F* fLambdaLambdaAvgSeparationPiPiME = new TH1F("fLambdaLambdaAvgSeparationPiPiME","Relative distances of p_daugh-p_daugh ME",1000,0.,500.);
    fOutputTP->Add(fLambdaLambdaAvgSeparationPP);
    fOutputTP->Add(fLambdaLambdaAvgSeparationPPME);
    fOutputTP->Add(fLambdaLambdaAvgSeparationP1Pi2);
    fOutputTP->Add(fLambdaLambdaAvgSeparationP1Pi2ME);
    fOutputTP->Add(fLambdaLambdaAvgSeparationP2Pi1);
    fOutputTP->Add(fLambdaLambdaAvgSeparationP2Pi1ME);
    fOutputTP->Add(fLambdaLambdaAvgSeparationPiPi);
    fOutputTP->Add(fLambdaLambdaAvgSeparationPiPiME);
  }


  //AntiLambda-Lambda
  TH1F* fAntiLambdaLambdaRelK = new TH1F("fAntiLambdaLambdaRelK","Relative momentum of AntiV0-V0 RE",150,0.,3.);
  TH1F* fAntiLambdaLambdaRelKME = new TH1F("fAntiLambdaLambdaRelKME","Relative momentum of V0-V0 ME",150,0.,3.);
  if(!fIsLightweight) {
    TH1F* fAntiLambdaAntiLambdaRelKMulti[13];
    TH1F* fAntiLambdaAntiLambdaRelKMEMulti[13];
    for(int i=0; i<13; ++i) {
        fAntiLambdaAntiLambdaRelKMulti[i] = new TH1F(Form("fAntiLambdaAntiLambdaRelKMulti_%i", i), Form("Relative momentum of AntiV0-AntiV0 RE %i", i),150,0.,3.); //first choice was 300 bins
        fAntiLambdaAntiLambdaRelKMEMulti[i] = new TH1F(Form("fAntiLambdaAntiLambdaRelKMEMulti_%i", i), Form("Relative momentum of AntiV0-AntiV0 ME %i", i),150,0.,3.);
        fOutputTP->Add(fAntiLambdaAntiLambdaRelKMulti[i]);
        fOutputTP->Add(fAntiLambdaAntiLambdaRelKMEMulti[i]);
    }
  }

  //AntiLambda-AntiLambda
  TH1F* fAntiLambdaAntiLambdaRelK = new TH1F("fAntiLambdaAntiLambdaRelK","Relative momentum of AntiV0-AntiV0 RE",150,0.,3.);
  TH1F* fAntiLambdaAntiLambdaRelKME = new TH1F("fAntiLambdaAntiLambdaRelKME","Relative momentum of AntiV0-AntiV0 ME",150,0.,3.);

  if(!fIsLightweight) {
    TH1F* fAntiLambdaAntiLambdaAvgSeparationPP = new TH1F("fAntiLambdaAntiLambdaAvgSeparationPP","Relative distances of p_daugh-p_daugh RE",1000,0.,500.);
    TH1F* fAntiLambdaAntiLambdaAvgSeparationPPME = new TH1F("fAntiLambdaAntiLambdaAvgSeparationPPME","Relative distances of p_daugh-p_daugh ME",1000,0.,500.);
    TH1F* fAntiLambdaAntiLambdaAvgSeparationP1Pi2 = new TH1F("fAntiLambdaAntiLambdaAvgSeparationP1Pi2","Relative distances of ppi_daugh RE",1000,0.,500.);
    TH1F* fAntiLambdaAntiLambdaAvgSeparationP1Pi2ME = new TH1F("fAntiLambdaAntiLambdaAvgSeparationP1Pi2ME","Relative distances of ppi_daugh ME",1000,0.,500.);
    TH1F* fAntiLambdaAntiLambdaAvgSeparationP2Pi1 = new TH1F("fAntiLambdaAntiLambdaAvgSeparationP2Pi1","Relative distances of ppi_daugh RE",1000,0.,500.);
    TH1F* fAntiLambdaAntiLambdaAvgSeparationP2Pi1ME = new TH1F("fAntiLambdaAntiLambdaAvgSeparationP2Pi1ME","Relative distances of ppi_daugh ME",1000,0.,500.);
    TH1F* fAntiLambdaAntiLambdaAvgSeparationPiPi = new TH1F("fAntiLambdaAntiLambdaAvgSeparationPiPi","Relative distances of p_daugh-p_daugh RE",1000,0.,500.);
    TH1F* fAntiLambdaAntiLambdaAvgSeparationPiPiME = new TH1F("fAntiLambdaAntiLambdaAvgSeparationPiPiME","Relative distances of p_daugh-p_daugh ME",1000,0.,500.);
    fOutputTP->Add(fAntiLambdaAntiLambdaAvgSeparationPP);
    fOutputTP->Add(fAntiLambdaAntiLambdaAvgSeparationPPME);
    fOutputTP->Add(fAntiLambdaAntiLambdaAvgSeparationP1Pi2);
    fOutputTP->Add(fAntiLambdaAntiLambdaAvgSeparationP1Pi2ME);
    fOutputTP->Add(fAntiLambdaAntiLambdaAvgSeparationP2Pi1);
    fOutputTP->Add(fAntiLambdaAntiLambdaAvgSeparationP2Pi1ME);
    fOutputTP->Add(fAntiLambdaAntiLambdaAvgSeparationPiPi);
    fOutputTP->Add(fAntiLambdaAntiLambdaAvgSeparationPiPiME);
  }


  TH1F* fNProtonLambdaPairs = new TH1F("fNProtonLambdaPairs","Number of V0-p pairs in an event",200,0,100);
  TH1F* fNPairStatistics = new TH1F("fNPairStatistics","Number of pairs under certain condition",40,0,10);


  TH1F* fPairkTLp = new TH1F("fPairkTLp","total transverse momentum of V0-p pair",1000,0.,10.);
  TH1F* fPairkTpp = new TH1F("fPairkTpp","total transverse momentum of p-p pair",1000,0.,10.);
  TH1F* fPairkTApAp = new TH1F("fPairkTApAp","total transverse momentum of Ap-Ap pair",1000,0.,10.);
  TH1F* fPairkTApp = new TH1F("fPairkTApp","total transverse momentum of Ap-p pair",1000,0.,10.);
  const Int_t NumberTPCrad = 9;

  if(!fIsLightweight) {
    TH2F* fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughProton[NumberTPCrad];
    TH2F* fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughProtonME[NumberTPCrad];
    TH2F* fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughPion[NumberTPCrad];
    TH2F* fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughPionME[NumberTPCrad];
    TH2F* fDeltaEtaDeltaPhiTPCradpp[NumberTPCrad];
    TH2F* fDeltaEtaDeltaPhiTPCradppME[NumberTPCrad];

    for(Int_t i=0;i<NumberTPCrad;i++) //close track eff. as a function of the TPC radius
      {
        TString HistDesc = "#Delta #eta vs. #Delta #phi for primary proton and decay proton evaluated at TPCradnum = ";
        HistDesc += i;
        TString HistName = "fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughProton";
        HistName += i;
        fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughProton[i] = new TH2F(HistName.Data(),HistDesc.Data(),400,0,4,400,0,1.);
        fOutputTP->Add(fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughProton[i]);

        HistDesc = "#Delta #eta vs. #Delta #phi for primary proton and decay pion evaluated at TPCradnum = ";
        HistDesc += i;
        HistName = "fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughPion";
        HistName += i;
        fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughPion[i] = new TH2F(HistName.Data(),HistDesc.Data(),400,0,4,400,0,1.);
        fOutputTP->Add(fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughPion[i]);

        HistDesc = "#Delta #eta vs. #Delta #phi for primary proton and decay proton in mixed event evaluated at TPCradnum = ";
        HistDesc += i;
        HistName = "fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughProtonME";
        HistName += i;
        fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughProtonME[i] = new TH2F(HistName.Data(),HistDesc.Data(),400,0,4,400,0,1.);
        fOutputTP->Add(fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughProtonME[i]);

        HistDesc = "#Delta #eta vs. #Delta #phi for primary proton and decay pion in mixed event evaluated at TPCradnum = ";
        HistDesc += i;
        HistName = "fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughPionME";
        HistName += i;
        fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughPionME[i] = new TH2F(HistName.Data(),HistDesc.Data(),400,0,4,400,0,1.);
        fOutputTP->Add(fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughPionME[i]);

        HistDesc = "#Delta #eta vs. #Delta #phi for proton-proton pairs evaluated at TPCradnum = ";
        HistDesc += i;
        HistName = "fDeltaEtaDeltaPhiTPCradpp";
        HistName += i;
        fDeltaEtaDeltaPhiTPCradpp[i] = new TH2F(HistName.Data(),HistDesc.Data(),400,0,4,400,0,1.);
        fOutputTP->Add(fDeltaEtaDeltaPhiTPCradpp[i]);

        HistDesc = "#Delta #eta vs. #Delta #phi for proton-proton pairs evaluated at TPCradnum = ";
        HistDesc += i;
        HistName = "fDeltaEtaDeltaPhiTPCradppME";
        HistName += i;
        fDeltaEtaDeltaPhiTPCradppME[i] = new TH2F(HistName.Data(),HistDesc.Data(),400,0,4,400,0,1.);
        fOutputTP->Add(fDeltaEtaDeltaPhiTPCradppME[i]);
     }
  }


  if(fUseMCInfo) {
    TH2F* fV0PtTrueReco = new TH2F("fV0PtTrueReco","V0: true transverse momentum vs. reconstructed momentum pt",200,0,10,200,0,10);
    TH2F* fV0PtResoRecoRotated = new TH2F("fV0PtResoRecoRotated","V0: true transverse momentum vs. reconstructed momentum pt",200,0,10,1000,-1,1);
    TH2F* fV0PhiResoRecoRotated = new TH2F("fV0PhiResoRecoRotated","V0: true phi vs. reconstructed momentum pt",200,0,10,3000,-1,1);
    TH2F* fV0ThetaResoRecoRotated = new TH2F("fV0ThetaResoRecoRotated","V0: true theta vs. reconstructed momentum pt",200,0,10,3000,-1,1);
    TH2F* fProtonPtTrueReco = new TH2F("fProtonPtTrueReco","Proton: true transverse momentum vs. reconstructed momentum pt",600,0,3,600,0,10);
    TH2F* fProtonPtResoRecoRotated = new TH2F("fProtonPtResoRecoRotated","Proton: true momentum vs. reconstructed momentum pt in rotated coordinate system",200,0,10,1000,-1,1);
    TH2F* fProtonPhiResoRecoRotated = new TH2F("fProtonPhiResoRecoRotated","Proton: true phi vs. reconstructed momentum pt",200,0,10,3000,-1,1);
    TH2F* fProtonThetaResoRecoRotated = new TH2F("fProtonThetaResoRecoRotated","Proton: true theta vs. reconstructed momentum pt",200,0,10,3000,-1,1);
    TH2F* fV0pRelKTrueReco = new TH2F("fV0pRelKTrueReco","V0p: true momentum vs. reconstructed momentum relK",600,0.,3.,600,0.,3.);
    TH2F* fV0pRelKTrueRecoRotated = new TH2F("fV0pRelKTrueRecoRotated","V0p: true momentum vs. reconstructed momentum relK",600,0.,3.,1000,-1.,1.);
    TH2F* fV0V0RelKTrueReco = new TH2F("fV0V0RelKTrueReco","V0V0: true momentum vs. reconstructed momentum relK",600,0.,3.,600,0.,3.);
    TH2F* fV0V0RelKTrueRecoRotated = new TH2F("fV0V0RelKTrueRecoRotated","V0p: true momentum vs. reconstructed momentum relK",600,0.,3.,1000,-1.,1.);
    TH2F* fPpRelKTrueReco = new TH2F("fPpRelKTrueReco","pp: true momentum vs. reconstructed momentum relK",750,0.,3.,750,0.,3.);
    TH2F* fPpRelKTrueRecoFB = new TH2F("fPpRelKTrueRecoFB","pp: true momentum vs. reconstructed momentum relK",1500,0.,3.,1500,0.,3.);
    TH2F* fPpRelKTrueRecoRotated = new TH2F("fPpRelKTrueRecoRotated","pp: true momentum vs. reconstructed momentum relK in rotated coordinate system",750,0.,3.,1000,-1.,1.);
    fOutputTP->Add(fV0PtTrueReco);
    fOutputTP->Add(fV0PtResoRecoRotated);
    fOutputTP->Add(fV0PhiResoRecoRotated);
    fOutputTP->Add(fV0ThetaResoRecoRotated);
    fOutputTP->Add(fProtonPtTrueReco);
    fOutputTP->Add(fProtonPtResoRecoRotated);
    fOutputTP->Add(fProtonPhiResoRecoRotated);
    fOutputTP->Add(fProtonThetaResoRecoRotated);
    fOutputTP->Add(fV0pRelKTrueReco);
    fOutputTP->Add(fV0pRelKTrueRecoRotated);
    fOutputTP->Add(fV0V0RelKTrueReco);
    fOutputTP->Add(fV0V0RelKTrueRecoRotated);
    fOutputTP->Add(fPpRelKTrueReco);
    fOutputTP->Add(fPpRelKTrueRecoFB);
    fOutputTP->Add(fPpRelKTrueRecoRotated);

    TH2F *fProtonPtCorrelationRE = new TH2F("fProtonPtCorrelationRE","correlation pt1-pt2 in same event",100,0.,5.,100,0.,5.);
    TH2F *fProtonPtCorrelationME = new TH2F("fProtonPtCorrelationME","correlation pt1-pt2 in mixed event",100,0.,5.,100,0.,5.);
    TH2F *fProtonEtaCorrelationRE = new TH2F("fProtonEtaCorrelationRE","correlation eta1-eta2 in same event",100,-1.,1.,100,-1.,1.);
    TH2F *fProtonEtaCorrelationME = new TH2F("fProtonEtaCorrelationME","correlation eta1-eta2 in mixed event",100,-1.,1.,100,-1.,1.);
    TH2F *fProtonEnCorrelationRE = new TH2F("fProtonEnCorrelationRE","correlation En1-En2 in same event",100,0.,5.,100,0.,5.);
    TH2F *fProtonEnCorrelationME = new TH2F("fProtonEnCorrelationME","correlation En1-En2 in mixed event",100,0.,5.,100,0.,5.);
    TH2F *fProtonPhiCorrelationRE = new TH2F("fProtonPhiCorrelationRE","correlation phi1-phi2 in same event",100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());
    TH2F *fProtonPhiCorrelationME = new TH2F("fProtonPhiCorrelationME","correlation phi1-phi2 in mixed event",100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());

    fOutputTP->Add(fProtonPtCorrelationRE);
    fOutputTP->Add(fProtonPtCorrelationME);
    fOutputTP->Add(fProtonEtaCorrelationRE);
    fOutputTP->Add(fProtonEtaCorrelationME);
    fOutputTP->Add(fProtonEnCorrelationRE);
    fOutputTP->Add(fProtonEnCorrelationME);
    fOutputTP->Add(fProtonPhiCorrelationRE);
    fOutputTP->Add(fProtonPhiCorrelationME);

  }

  fOutputTP->Add(fProtonLambdaMT);
  fOutputTP->Add(fProtonLambdaMTrelKcut);
  fOutputTP->Add(fProtonLambdaRelK);
  fOutputTP->Add(fProtonLambdaRelKME);

  fOutputTP->Add(fAntiProtonAntiLambdaRelK);
  fOutputTP->Add(fAntiProtonAntiLambdaRelKME);

  fOutputTP->Add(fProtonAntiLambdaRelK);
  fOutputTP->Add(fProtonAntiLambdaRelKME);
  fOutputTP->Add(fAntiProtonLambdaRelK);
  fOutputTP->Add(fAntiProtonLambdaRelKME);
  fOutputTP->Add(fProtonXiRelK);
  fOutputTP->Add(fProtonXiRelKME);
  fOutputTP->Add(fAntiProtonXiRelK);
  fOutputTP->Add(fAntiProtonXiRelKME);
  fOutputTP->Add(fProtonProtonRelK);
  fOutputTP->Add(fProtonProtonRelKME);

  fOutputTP->Add(fProtonProtonMT);
  fOutputTP->Add(fProtonProtonMTrelKcut);
  fOutputTP->Add(fAntiProtonAntiProtonRelK);
  fOutputTP->Add(fAntiProtonAntiProtonRelKME);

  fOutputTP->Add(fAntiProtonProtonRelK);
  fOutputTP->Add(fAntiProtonProtonRelKME);
  fOutputTP->Add(fLambdaLambdaRelK);
  fOutputTP->Add(fLambdaLambdaRelKME);

  fOutputTP->Add(fAntiLambdaAntiLambdaRelK);
  fOutputTP->Add(fAntiLambdaAntiLambdaRelKME);

  fOutputTP->Add(fAntiLambdaLambdaRelK);
  fOutputTP->Add(fAntiLambdaLambdaRelKME);

  fOutputTP->Add(fNProtonLambdaPairs);
  fOutputTP->Add(fNPairStatistics);

  fOutputTP->Add(fPairkTLp);
  fOutputTP->Add(fPairkTpp);
  fOutputTP->Add(fPairkTApAp);
  fOutputTP->Add(fPairkTApp);

  if(!fIsLightweight) {
    //Define PID related histograms:
    TH2F* fdEdxVsP = new TH2F("fdEdxVsP","dE/dx of all particles vs momentum",1000,0,7,400,-10,200);
    TH2F* fTOFSignalVsP = new TH2F ("fTOFSignalVsP","tof signal vs p (positives);p [GeV/c];t_{meas} - t_{0} - t_{expected} [ps]",20,0.0,5.0,120,-10000.0,5000.0);
    TH2F* fProtonNSigmaTPC = new TH2F("fProtonNSigmaTPC","nsigma_TPC of protons for P<0.75 GeV/c without TOF signal",100,0,1,100,-1.,5.);
    TH2F* fProtonNSigmaCombined = new TH2F("fProtonNSigmaCombined","nsigma_combined of protons for P>0.75 GeV/c",500,0.5,7,100,-1.,5.);

    TH2F* fNsigmaTOFTPCPtBin[fCuts->GetPtBinsPurity()];
    TH2F* fNsigmaTOFTPCHypothesisIncludedPtBin[fCuts->GetPtBinsPurity()];


    TH1F* fProtonsCorrectlyIdentified = new TH1F("fProtonsCorrectlyIdentified","Particles that were correctly identified as protons",100,fCuts->GetProtonPIDthrPtLow(),fCuts->GetProtonPIDthrPtUp());
    TH1F* fProtonsTotallyIdentified = new TH1F("fProtonsTotallyIdentified","All Particles that were identified as protons",100,fCuts->GetProtonPIDthrPtLow(),fCuts->GetProtonPIDthrPtUp());
    TH1F* fProtonsBkgIdentified = new TH1F("fProtonsBkgIdentified","All Particles that were identified as protons",100,fCuts->GetProtonPIDthrPtLow(),fCuts->GetProtonPIDthrPtUp());


    if(fUseMCInfo)
      {
        fOutputPID->Add(fProtonsCorrectlyIdentified);
        fOutputPID->Add(fProtonsTotallyIdentified);
        fOutputPID->Add(fProtonsBkgIdentified);
      }

    for(Int_t i=0;i<fCuts->GetPtBinsPurity();i++)
      {
        TString HistName = "fNsigmaTOFTPCPtBin";
        HistName += i;
        TString Description = "nsigma of TPC vs nsigma of TOF for momentum bin ";
        Description += i;
        fNsigmaTOFTPCPtBin[i] = new TH2F(HistName.Data(),Description.Data(),200,-10,10,200,-10,10);
        fOutputPID->Add(fNsigmaTOFTPCPtBin[i]);

        HistName = "fNsigmaTOFTPCHypothesisIncludedPtBin";
        HistName += i;
        Description = "nsigma of TPC vs nsigma of TOF for momentum bin ";
        Description += i;
        fNsigmaTOFTPCHypothesisIncludedPtBin[i] = new TH2F(HistName.Data(),Description.Data(),200,-10,10,200,-10,10);
        fOutputPID->Add(fNsigmaTOFTPCHypothesisIncludedPtBin[i]);
      }

    fOutputPID->Add(fdEdxVsP);
    fOutputPID->Add(fTOFSignalVsP);
    fOutputPID->Add(fProtonNSigmaTPC);
    fOutputPID->Add(fProtonNSigmaCombined);
  }

  return;
}

Bool_t AliAnalysisTaskPLFemto::SelectPID(const AliAODTrack *track,const AliAODEvent *aodevent,Int_t type)
{
  //Comment: In the SingleParticleQA function goodDCA is called after SelectPID. To have plots which include already
  //the goodDCA decision it was partly included in the SelectPID function even if the clarity gets a bit lost
  Double_t charge = track->Charge();
  if(type == 4 && charge<0.) return kFALSE;
  else if(type == -4 && charge>0.) return kFALSE;


  const AliAODTrack *globaltrack;
  if(fWhichfilterbit == 128) globaltrack = fGTI[-track->GetID()-1];
  else globaltrack = track;

  Bool_t isparticle = kFALSE;

  TClonesArray *mcarray = 0;
  AliAODMCParticle* mcproton = 0;
  if(fUseMCInfo)
    {
      mcarray = dynamic_cast<TClonesArray*>(aodevent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!mcarray) return kFALSE;

      if(track->GetLabel() >=0) mcproton = (AliAODMCParticle*)mcarray->At(track->GetLabel());
      else return kFALSE;//check this flag out! it rejects in principle all fake tracks in the analysis
    }

  //there must be TPC & TOF signal (TOF for P>0.75 GeV/c)
  Bool_t TPCisthere = kFALSE;
  Bool_t TOFisthere = kFALSE;
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,globaltrack);
  if(statusTPC == AliPIDResponse::kDetPidOk) TPCisthere = kTRUE;
  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,globaltrack);
  if(statusTOF == AliPIDResponse::kDetPidOk) TOFisthere = kTRUE;
  if(!TPCisthere) return kFALSE; //TPC signal must be there

  TString fillthis = "fdEdxVsP";
  if(!fIsLightweight) ((TH2F*)(fOutputPID->FindObject(fillthis)))->Fill(globaltrack->GetTPCmomentum(),globaltrack->GetTPCsignal());

  //Now do particle identification:
  if(globaltrack->GetTPCmomentum() <= fCuts->GetProtonPIDTOFTPCSwitch())//for P<0.75 use TPC
    {
      Double_t nSigmaProtonTPConly = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(globaltrack,(AliPID::kProton)));
      fillthis = "fProtonNSigmaTPC";
      if(type == 4 && !fIsLightweight)((TH2F*)(fOutputPID->FindObject(fillthis)))->Fill(globaltrack->GetTPCmomentum(),nSigmaProtonTPConly);
      if(nSigmaProtonTPConly < fCuts->GetProtonnsigma()) isparticle = kTRUE;
    }
  else if(globaltrack->GetTPCmomentum() > fCuts->GetProtonPIDTOFTPCSwitch())//for P>0.75 GeV/c use TPC&TOF
    {
      if(!TOFisthere) return kFALSE;//makes only sense if TOF PID is available (TPC signal was checked already at the beginning)

      Double_t corrTOFsignal = GetCorrectedTOFSignal(globaltrack,aodevent);

      Double_t nSigmaTPCproton = fPIDResponse->NumberOfSigmasTPC(globaltrack,(AliPID::kProton));
      Double_t nSigmaTOFproton = fPIDResponse->NumberOfSigmasTOF(globaltrack,(AliPID::kProton));
      Double_t nSigmaTPCTOFcombined = TMath::Sqrt(pow(nSigmaTPCproton,2.) + pow(nSigmaTOFproton,2.));

      if(nSigmaTPCTOFcombined < fCuts->GetProtonnsigma() && TestPIDHypothesis(globaltrack,type)) isparticle = kTRUE;

      if(type == 4)
	{
	  fillthis = "fProtonNSigmaCombined";
    if(!fIsLightweight) ((TH2F*)(fOutputPID->FindObject(fillthis)))->Fill(globaltrack->P(),nSigmaTPCTOFcombined);
	  fillthis = "fTOFSignalVsP";
    if(!fIsLightweight) ((TH2F*)(fOutputPID->FindObject(fillthis)))->Fill(globaltrack->GetTPCmomentum(),corrTOFsignal);
	}


      //QA the PID differentially as a function of pt
      Int_t ptbin = fCuts->FindProtonPtBinPurity(globaltrack->GetTPCmomentum());
      //Int_t ptbin = fCuts->FindProtonPtBinPurity(track->Pt());

      if(ptbin != -9999)
	{
	  //investigate this as a function of momentum slices:
	  if(type == 4)
	    {
	      fillthis = "fNsigmaTOFTPCPtBin";
	      fillthis += ptbin;
        if(!fIsLightweight) ((TH2F*)(fOutputPID->FindObject(fillthis)))->Fill(nSigmaTPCproton,nSigmaTOFproton);

	      if(TestPIDHypothesis(globaltrack,type))
		{
		  fillthis = "fNsigmaTOFTPCHypothesisIncludedPtBin";
		  fillthis += ptbin;
      if(!fIsLightweight) ((TH2F*)(fOutputPID->FindObject(fillthis)))->Fill(nSigmaTPCproton,nSigmaTOFproton);
		}
	    }
	}
    }


  //Investigate the purity for the selected particles (in Monte Carlo):
  if(isparticle)
    {
      if(fUseMCInfo)
	{
	  //check with monte carlo the purity:
	  if(type == 4)
	    {
	      fillthis = "fProtonsTotallyIdentified";
        if(!fIsLightweight) ((TH1F*)(fOutputPID->FindObject(fillthis)))->Fill(track->Pt());
	    }
	  if(type == 4 && mcproton->GetPdgCode() == 2212)
	    {
	      fillthis = "fProtonsCorrectlyIdentified";
        if(!fIsLightweight) ((TH1F*)(fOutputPID->FindObject(fillthis)))->Fill(track->Pt());
	    }
	  else if(type == 4 && mcproton->GetPdgCode() != 2212)
	    {
	      fillthis = "fProtonsBkgIdentified";
        if(!fIsLightweight) ((TH1F*)(fOutputPID->FindObject(fillthis)))->Fill(track->Pt());
	    }
	}
    }


  return isparticle;
}


Bool_t AliAnalysisTaskPLFemto::TestPIDHypothesis(const AliAODTrack *track,Int_t type)
{
  //For P>0.75 GeV/c the TPC dE/dx bands start to merge, we apply an additional test to the PID of the particle (proton)

  //const AliAODTrack *globaltrack;
  //if(fWhichfilterbit == 128) globaltrack = fGTI[-track->GetID()-1];
  //else globaltrack = track;

  Double_t charge = track->Charge();
  if(type == 4 && charge<0.) return kFALSE;
  else if(type == -4 && charge>0.) return kFALSE;

  Bool_t keepParticle = kTRUE;
  if(track->GetTPCmomentum() <= fCuts->GetProtonPIDTOFTPCSwitch()) return kTRUE;//don't check for low momentum tracks
  else
    {
      AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,track);
      AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,track);
      if(!statusTPC) return kFALSE; //TPC signal must be present, BUG 30.01.2017: Was kTRUE must be kFALSE but it is anyhow checked before in SelectPID
      if(!statusTOF) return kFALSE; //TOF signal must be present

      Double_t nSigmaProtonTPC = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::kProton)));
      Double_t nSigmaKaonTPC = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::kKaon)));
      Double_t nSigmaPionTPC = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::kPion)));
      Double_t nSigmaElectronTPC = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::kElectron)));

      Double_t nSigmaProtonTOF = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,(AliPID::kProton)));
      Double_t nSigmaKaonTOF = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,(AliPID::kKaon)));
      Double_t nSigmaPionTOF = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,(AliPID::kPion)));
      Double_t nSigmaElectronTOF = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,(AliPID::kElectron)));

      //Form the combination:
      Double_t nSigmaProtonCombined = TMath::Sqrt(pow(nSigmaProtonTPC,2.) + pow(nSigmaProtonTOF,2.));
      Double_t nSigmaKaonCombined = TMath::Sqrt(pow(nSigmaKaonTPC,2.) + pow(nSigmaKaonTOF,2.));
      Double_t nSigmaPionCombined = TMath::Sqrt(pow(nSigmaPionTPC,2.) + pow(nSigmaPionTOF,2.));
      Double_t nSigmaElectronCombined = TMath::Sqrt(pow(nSigmaElectronTPC,2.) + pow(nSigmaElectronTOF,2.));

      //The other hypothesis fits better:
      if((nSigmaKaonCombined < nSigmaProtonCombined) || (nSigmaPionCombined < nSigmaProtonCombined) || (nSigmaElectronCombined < nSigmaProtonCombined)) keepParticle = kFALSE;
    }

  return keepParticle;
}

Bool_t AliAnalysisTaskPLFemto::SelectV0PID(const AliAODv0 *v0,TString ParOrAPar)
{
  AliAODTrack* pTrack = fGTI[v0->GetPosID()];
  AliAODTrack* nTrack = fGTI[v0->GetNegID()];

  Double_t charge1 = pTrack->Charge();
  Double_t charge2 = nTrack->Charge();

  if(charge1 < 0. && charge2 > 0. && ParOrAPar == fwhichV0)//assign correct charge to tracks from V0
    {
      pTrack = fGTI[v0->GetNegID()];//proton as positive particle
      nTrack = fGTI[v0->GetPosID()];//pion as negative particle
    }
  else if(charge1 > 0. && charge2 <0. && ParOrAPar == fwhichAntiV0)//assign correct charge to tracks from Anti-V0
    {
      pTrack = fGTI[v0->GetNegID()];//anti-proton as negative particle
      nTrack = fGTI[v0->GetPosID()];//anti-pion as positive particle
    }
  else if((charge1>0. && charge2>0.) || (charge1<0. && charge2<0.))
    {
      return kFALSE;
    }

  //new method for checking detector status:
  AliPIDResponse::EDetPidStatus statusPosTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,pTrack);
  if(statusPosTPC != AliPIDResponse::kDetPidOk) return kFALSE;
  AliPIDResponse::EDetPidStatus statusNegTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,nTrack);
  if(statusNegTPC != AliPIDResponse::kDetPidOk) return kFALSE;

  Bool_t isV0 = kFALSE;
  Bool_t isProton = kFALSE;
  Bool_t isPion = kFALSE;

  Double_t nSigmaProtonTPC = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,(AliPID::kProton)));

  if(nSigmaProtonTPC < fCuts->GetV0nsigma()) isProton = kTRUE;

  Double_t nSigmaPionTPC = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack,(AliPID::kPion)));
  if(nSigmaPionTPC < fCuts->GetV0nsigma()) isPion = kTRUE;

  if(isProton && isPion) isV0 = kTRUE;

  /*
  AliAODTrack* daughterTrackPos = (AliAODTrack*)v0->GetDaughter(0);
  AliAODTrack* daughterTrackNeg = (AliAODTrack*)v0->GetDaughter(1);

  AliAODMCParticle *posTrackMC  = (AliAODMCParticle*)fAODMCEvent->GetTrack(daughterTrackPos->GetLabel());
  if(!posTrackMC) return kFALSE;

  AliAODMCParticle *negTrackMC  = (AliAODMCParticle*)fAODMCEvent->GetTrack(daughterTrackNeg->GetLabel());
  if(!negTrackMC) return kFALSE;


  if(ParOrAPar == fwhichV0)
    {

      if(posTrackMC->GetPdgCode() == 2212 && negTrackMC->GetPdgCode()==-211)
	{

	  std::cout << "particle: " << fwhichV0 << std::endl;
	  std::cout << "posTrackPDG: "  << posTrackMC->GetPdgCode() << std::endl;
	  std::cout << "negTrackPDG: "  << negTrackMC->GetPdgCode() << std::endl;
	  std::cout << "Px: " << daughterTrackPos->Px() << std::endl;
	  std::cout << "px2: " << pTrack->Px() << std::endl;
	  std::cout << ".................................." << std::endl;
	}
    }
  else if(ParOrAPar == fwhichAntiV0)
    {
      if(posTrackMC->GetPdgCode() == 211 && negTrackMC->GetPdgCode()==-2212)
	{

	  std::cout << "particle: " << fwhichAntiV0 << std::endl;
	  std::cout << "posTrackPDG: "  << posTrackMC->GetPdgCode() << std::endl;
	  std::cout << "negTrackPDG: "  << negTrackMC->GetPdgCode() << std::endl;
	  std::cout << "Px: " << daughterTrackPos->Px() << std::endl;
	  std::cout << "px2: " << pTrack->Px() << std::endl;
	  std::cout << ".................................." << std::endl;
	}
    }
  */

  return isV0;
}

//________________________________________________________________________
TVector3 AliAnalysisTaskPLFemto::GetGlobalPositionAtGlobalRadiiThroughTPC(const AliAODTrack *track, const Float_t bfield,double Rwanted,float PrimaryVertex[3])
{
  // Gets the global position of the track at nine different radii in the TPC
  // track is the track you want to propagate
  // bfield is the magnetic field of your event
  // globalPositionsAtRadii is the array of global positions in the radii and xyz


  //Directly from GitHub:  Add method GetXYZatR for fast estimate of track intersection with given X or R
  //https://gitlab.cern.ch/mkrzewic/AliRoot/commit/c8253dd95447cc9624b1761bef3603d7171d8759

  // This method has 3 modes of behaviour
  // 1) xyz[3] array is provided but alpSect pointer is 0: calculate the position of track intersection
  //    with circle of radius xr and fill it in xyz array
  // 2) alpSect pointer is provided: find alpha of the sector where the track reaches local coordinate xr
  //    Note that in this case xr is NOT the radius but the local coordinate.
  //    If the xyz array is provided, it will be filled by track lab coordinates at local X in this sector
  // 3) Neither alpSect nor xyz pointers are provided: just check if the track reaches radius xr

  //*********************************************
  //How the method (probably) works:

  //helix is obtained at starting point of track e.g. at DCA or first track points
  //this track helix is then intersected with a circle centered around the origin

  //*********************************************

  //How to take into account in mixed event hat one has two separated primary vertices?
  //1) Shift at least along the z-coordinate
  //2) Shift the whole primary vertex

  //Method 1) does not influence the Rwanted

  //For method 2) there is basically only one possibility left, to change the Rwanted value

  //If we shift the primary vertex to be in the origin of the coordinate system, it is the same as shifting the circle to the primary vertex
  //But then the track hits not the position it would hit without offset. Well there are always drawbacks

  //If we have a vector R which points from the origin to the intersection point of the track with the radius, and a vector x_0 which contains the primary vertex (x,y), then we have to calculate the shifted vector a = R - x_0
  //The length of the new vector |a|² is the radius we would obtain if all primary vertices would lie at the origin


  //Try to shift only along z, this is the worst direction and easiest to calculate:

  AliExternalTrackParam exParam;
  exParam.CopyFromVTrack(track);
  double Position[3] = {-9999.,-9999.,-9999.};

  TVector3 globalPositionsAtRadii;
  globalPositionsAtRadii.SetXYZ(Position[0],Position[1],Position[2]);


  if(!exParam.GetXYZatR(Rwanted,bfield,Position,NULL)) return globalPositionsAtRadii;
  else
    {
      //Shift z-direction:
      Position[2] -= PrimaryVertex[2];

      globalPositionsAtRadii.SetXYZ(Position[0],Position[1],Position[2]);//set the position to the correct value if the propagation worked
      return globalPositionsAtRadii;
    }


  /*
    AliExternalTrackParam exParam[9];//For every track an individual copy

  for(Int_t RNumber = 0; RNumber < 9; RNumber++)
    {
    exParam[RNumber].CopyFromVTrack(track);

    double Position[3] = {-9999.,-9999.,-9999.};
      globalPositionsAtRadii[RNumber].SetXYZ(Position[0],Position[1],Position[2]);//set the position onto unnatural large values to see that something goes wrong
      if(!exParam[RNumber].GetXYZatR(TPCradii[RNumber],bfield,Position,NULL)) continue;
      globalPositionsAtRadii[RNumber].SetXYZ(Position[0],Position[1],Position[2]);//set the position to the correct value if the propagation worked
    }
  */

  //OLD method:

  /*
  // Initialize the array to something indicating there was no propagation
  for(Int_t i=0;i<9;i++)
    {
      for(Int_t j=0;j<3;j++)
	{
	  globalPositionsAtRadii[i][j]=-9999.;
	}
    }

  // Make a copy of the track to not change parameters of the track
  AliExternalTrackParam etp; etp.CopyFromVTrack(track);
  //printf("\nAfter CopyFromVTrack\n");
  //etp.Print();

  // The global position of the the track
  Double_t xyz[3]={-9999.,-9999.,-9999.};

  // Counter for which radius we want
  Int_t iR=0;
  // The radii at which we get the global positions
  // IROC (OROC) from 84.1 cm to 132.1 cm (134.6 cm to 246.6 cm)
  //Float_t Rwanted[9]={85.,105.,125.,145.,165.,185.,205.,225.,245.};


  // The global radius we are at
  Float_t globalRadius=0;

  // Propagation is done in local x of the track
  for (Float_t x = etp.GetX();x<247.;x+=1.)
    { 
      // GetX returns local coordinates
      // Starts at the tracks fX and goes outwards. x = 245 is the outer radial limit
      // of the TPC when the track is straight, i.e. has inifinite pt and doesn't get bent.
      // If the track's momentum is smaller than infinite, it will develop a y-component, which
      // adds to the global radius

      // Stop if the propagation was not succesful. This can happen for low pt tracks
      // that don't reach outer radii
      if(!etp.PropagateTo(x,bfield)) break;
      etp.GetXYZ(xyz); // GetXYZ returns global coordinates

      globalRadius = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]); //Idea to speed up: compare squared radii



      // Roughly reached the radius we want
      if(globalRadius > fTPCradii[iR])
	{
	  // Bigger loop has bad precision, we're nearly one centimeter too far, go back in small steps.
	  while (globalRadius>fTPCradii[iR])
	    {
	      x-=.1;
	      //printf("propagating to x %5.2f\n",x);
	      if(!etp.PropagateTo(x,bfield)) break;
	      etp.GetXYZ(xyz); // GetXYZ returns global coordinates

	      globalRadius = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]); //Idea to speed up: compare squared radii
	    }
	  //printf("At Radius:%05.2f (local x %5.2f). Setting position to x %4.1f y %4.1f z %4.1f\n",globalRadius,x,xyz[0],xyz[1],xyz[2]);

	  globalPositionsAtRadii[iR][0]=xyz[0];
	  globalPositionsAtRadii[iR][1]=xyz[1];
	  globalPositionsAtRadii[iR][2]=xyz[2];
	  //subtract primary vertex, "zero" track for correct mixed-event comparison
	  globalPositionsAtRadii[iR][0] -= PrimaryVertex[0];
	  globalPositionsAtRadii[iR][1] -= PrimaryVertex[1];
	  globalPositionsAtRadii[iR][2] -= PrimaryVertex[2];

	  // Indicate we want the next radius
	  iR+=1;
	}

      if(iR>8) return; //TPC edge reached
    }
  */
}
//________________________________________________________________________
Float_t AliAnalysisTaskPLFemto::DEtacalc(const double zprime1, const double zprime2, const double radius)
{
  float thetaS1 = TMath::Pi()/2. - TMath::ATan(zprime1/radius);
  float thetaS2 = TMath::Pi()/2. - TMath::ATan(zprime2/radius);

  float etaS1 = -TMath::Log(TMath::Tan(thetaS1/2.));
  float etaS2 = -TMath::Log(TMath::Tan(thetaS2/2.));

  float detaS = etaS1 - etaS2;

  return detaS;
}
//________________________________________________________________________
Float_t AliAnalysisTaskPLFemto::DPhicalc(const double dxprime, const double dyprime, const double radius)
{
  float dist = TMath::Sqrt(pow(dxprime,2.) + pow(dyprime,2.));
  float dphi = 2.*TMath::ATan(dist/(2.*radius));

  return dphi;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPLFemto::CheckIfRealV0(AliAODv0 *v0,AliAODEvent *aodEvent,Int_t PDGmother,Int_t PDGdaugh1, Int_t PDGdaugh2)
{

  //Checks according to daughter tracks and mother PDG if the V0 is real or fake

  Bool_t realV0 = kFALSE;
  TClonesArray *mcarray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if(!mcarray) return kFALSE;

  //Int_t daugh[2] = {2212,211};//daughter IDs of Lambda
  Int_t daugh[2];
  daugh[0] = PDGdaugh1;
  daugh[1] = PDGdaugh2;
  Int_t label = v0->MatchToMC(PDGmother,mcarray,2,daugh);
  //Int_t label = v0->MatchToMC(PDGmother,mcarray);

  if(label >= 0) realV0 = kTRUE;

  return realV0;
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::GetV0Origin(AliAODv0 *v0,AliAODEvent *aodEvent)
{
  //Checks where the V0 comes from, e.g. feed-down or primary

  TString fillthis = "";

  TClonesArray *mcarray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!mcarray) return;

  Int_t daugh[2] = {2212,211};//daughter IDs of Lambda
  Int_t label = v0->MatchToMC(3122,mcarray,2,daugh);
  if(label<0) return;// only "true" if it is a real v0

  AliAODMCParticle *MCv0 = (AliAODMCParticle*)mcarray->At(label);
  if(!MCv0) return;

  fillthis = "fFeeddownV0";
  if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(1);

  if(MCv0->IsPhysicalPrimary() && !MCv0->IsSecondaryFromWeakDecay())
    {
      fillthis = "fFeeddownV0";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(2);
    }
  else if(MCv0->IsSecondaryFromWeakDecay())
    {
      AliAODMCParticle* MCv0Mother = (AliAODMCParticle*)mcarray->At(MCv0->GetMother());

      fillthis = "fFeeddownV0";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(3);

      fillthis = "fFeeddownV0FromWeakPDG";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(MCv0Mother->PdgCode());
    }
  else if(MCv0->IsSecondaryFromMaterial())
    {
      fillthis = "fFeeddownV0";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(4);
    }
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::GetProtonOrigin(AliAODTrack *AODtrack,AliAODEvent *aodEvent,Int_t type)
{
  //Checks where the proton comes from, e.g. feed-down or primary

  Int_t desiredPDGcode = 0;
  if(type == 4) desiredPDGcode = 2212;
  else if(type == -4) desiredPDGcode = -2212;
  else std::cout << "type is wrong" << std::endl;

  //Momentum Matrix for single particle
  TString fillthis = "";
  TClonesArray *mcarray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!mcarray) return;

  if(AODtrack->GetLabel() < 0) return;

  AliAODMCParticle* mcproton = (AliAODMCParticle*)mcarray->At(AODtrack->GetLabel());
  if(!mcproton) return;

  if(mcproton->GetPdgCode() != desiredPDGcode) return; //be shure it is a proton/antiproton


  //check where the proton comes from in Monte Carlo:
  if(type == 4)
    {
      fillthis = "fFeeddownProton";
      if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(1);//total amount of protons

      if(mcproton->IsPhysicalPrimary() && !mcproton->IsSecondaryFromWeakDecay())
	{
	  fillthis = "fFeeddownProton";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(2);//primary and not from weak decay
	}
      else if(mcproton->IsSecondaryFromWeakDecay() && !mcproton->IsSecondaryFromMaterial())
	{
	  AliAODMCParticle* mcprotonMother = (AliAODMCParticle*)mcarray->At(mcproton->GetMother());

	  fillthis = "fFeeddownProton";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(3);//from weak decay

	  fillthis = "fFeeddownProtonFromWeakPDG";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(mcprotonMother->PdgCode());
	}
      else if(mcproton->IsSecondaryFromMaterial())
	{
	  fillthis = "fFeeddownProton";
    if(!fIsLightweight) ((TH1F*)(fOutputSP->FindObject(fillthis)))->Fill(4);//from material
	}
    }
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::GetMomentumMatrix(const AliFemtoLambdaParticle &v01,const AliFemtoLambdaParticle &v02)
{
  //This function calculates the effect of the finite momentum resolution of ALICE for v0s

  TLorentzVector trackV01,trackV01MC,trackV02,trackV02MC;

  trackV01.SetXYZM(v01.fMomentum.X(),v01.fMomentum.Y(),v01.fMomentum.Z(),fCuts->GetHadronMasses(3122));
  trackV01MC.SetXYZM(v01.fMomentumMC.X(),v01.fMomentumMC.Y(),v01.fMomentumMC.Z(),fCuts->GetHadronMasses(3122));

  trackV02.SetXYZM(v02.fMomentum.X(),v02.fMomentum.Y(),v02.fMomentum.Z(),fCuts->GetHadronMasses(3122));
  trackV02MC.SetXYZM(v02.fMomentumMC.X(),v02.fMomentumMC.Y(),v02.fMomentumMC.Z(),fCuts->GetHadronMasses(3122));

  if(trackV01MC.X() == -9999. || trackV02MC.X() == -9999.) return;
  if(trackV01MC.Y() == -9999. || trackV02MC.Y() == -9999.) return;
  if(trackV01MC.Z() == -9999. || trackV02MC.Z() == -9999.) return;

  Double_t relK = relKcalc(trackV01,trackV02);
  Double_t relKMC = relKcalc(trackV01MC,trackV02MC);

  TString fillthis = "fV0V0RelKTrueReco";
  if(fUseMCInfo) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(relKMC,relK);

  Double_t relKTruePrime = 1./TMath::Sqrt(2.)*(relKMC - relK);
  Double_t relKRecoPrime = 1./TMath::Sqrt(2.)*(relKMC + relK);
  fillthis = "fV0V0RelKTrueRecoRotated";
  if(fUseMCInfo) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(relKRecoPrime,relKTruePrime);
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::GetMomentumMatrix(const AliFemtoLambdaParticle &v0,const AliFemtoProtonParticle &proton)
{
  //This function calculates the effect of the finite momentum resolution of ALICE for v0s

  TLorentzVector trackV0,trackProton,trackV0MC,trackProtonMC;

  trackV0.SetXYZM(v0.fMomentum.X(),v0.fMomentum.Y(),v0.fMomentum.Z(),fCuts->GetHadronMasses(3122));
  trackProton.SetXYZM(proton.fMomentum.X(),proton.fMomentum.Y(),proton.fMomentum.Z(),fCuts->GetHadronMasses(2212));

  trackV0MC.SetXYZM(v0.fMomentumMC.X(),v0.fMomentumMC.Y(),v0.fMomentumMC.Z(),fCuts->GetHadronMasses(3122));
  trackProtonMC.SetXYZM(proton.fMomentumMC.X(),proton.fMomentumMC.Y(),proton.fMomentumMC.Z(),fCuts->GetHadronMasses(2212));

  if(trackV0MC.X() == -9999. || trackProtonMC.X() == -9999.) return;
  if(trackV0MC.Y() == -9999. || trackProtonMC.Y() == -9999.) return;
  if(trackV0MC.Z() == -9999. || trackProtonMC.Z() == -9999.) return;


  Double_t relK = relKcalc(trackV0,trackProton);
  Double_t relKMC = relKcalc(trackV0MC,trackProtonMC);

  TString fillthis = "fV0pRelKTrueReco";
  if(fUseMCInfo) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(relKMC,relK);

  Double_t relKTruePrime = 1./TMath::Sqrt(2.)*(relKMC - relK);
  Double_t relKRecoPrime = 1./TMath::Sqrt(2.)*(relKMC + relK);
  fillthis = "fV0pRelKTrueRecoRotated";
  if(fUseMCInfo) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(relKRecoPrime,relKTruePrime);
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::GetMomentumMatrix(const AliFemtoProtonParticle &proton1,const AliFemtoProtonParticle &proton2)
{
  //This function calculates the effect of the finite momentum resolution of ALICE for protons
  //It takes the input of the proton momentum and the reconstructed momentum

  TLorentzVector trackProton1,trackProton2,trackProton1MC,trackProton2MC;

  trackProton1.SetXYZM(proton1.fMomentum.X(),proton1.fMomentum.Y(),proton1.fMomentum.Z(),fCuts->GetHadronMasses(2212));
  trackProton2.SetXYZM(proton2.fMomentum.X(),proton2.fMomentum.Y(),proton2.fMomentum.Z(),fCuts->GetHadronMasses(2212));

  trackProton1MC.SetXYZM(proton1.fMomentumMC.X(),proton1.fMomentumMC.Y(),proton1.fMomentumMC.Z(),fCuts->GetHadronMasses(2212));
  trackProton2MC.SetXYZM(proton2.fMomentumMC.X(),proton2.fMomentumMC.Y(),proton2.fMomentumMC.Z(),fCuts->GetHadronMasses(2212));

  if(trackProton1MC.X() == -9999. || trackProton2MC.X() == -9999.) return;
  if(trackProton1MC.Y() == -9999. || trackProton2MC.Y() == -9999.) return;
  if(trackProton1MC.Z() == -9999. || trackProton2MC.Z() == -9999.) return;


  Double_t relK = relKcalc(trackProton1,trackProton2);
  Double_t relKMC = relKcalc(trackProton1MC,trackProton2MC);


  TString fillthis = "fPpRelKTrueReco";
  if(fUseMCInfo) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(relKMC,relK);

  fillthis = "fPpRelKTrueRecoFB";
  if(fUseMCInfo) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(relKMC,relK);


  Double_t relKTruePrime = 1./TMath::Sqrt(2.)*(relKMC - relK);
  Double_t relKRecoPrime = 1./TMath::Sqrt(2.)*(relKMC + relK);
  fillthis = "fPpRelKTrueRecoRotated";
  if(fUseMCInfo) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(relKRecoPrime,relKTruePrime);
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPLFemto::CheckMCPID(AliAODEvent *aodEvent, AliAODTrack* aodtrack,int pdg)
{
  //Checks the input track according to given PDG value

  TClonesArray *mcarray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if(!mcarray) return kFALSE;

  if(aodtrack->GetLabel() < 0) return kFALSE;

  AliAODMCParticle* mcPart = (AliAODMCParticle*)mcarray->At(aodtrack->GetLabel());
  if(!mcPart) return kFALSE;
  if(mcPart->GetPdgCode() == pdg) return kTRUE;
  else return kFALSE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPLFemto::CheckMCPID(AliAODv0 *v0,AliAODEvent* aodEvent,int pdg)
{
  //Checks if the input v0 is a real V0

  TClonesArray *mcarray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if(!mcarray) return kFALSE;

  Int_t daugh[2] = {2212,211};//daughter IDs of Lambda
  Int_t label = v0->MatchToMC(3122,mcarray,2,daugh);
  if(label<0) return kFALSE;// only "true" if it is a real v0

  AliAODMCParticle* mcPart = (AliAODMCParticle*)mcarray->At(label);
  if(!mcPart) return kFALSE;
  if(mcPart->GetPdgCode() == pdg) return kTRUE;
  else return kFALSE;
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::GetMCMomentum(AliAODTrack *aodtrack,double *Protontruth,double *ProtontruthMother,double *ProtontruthMotherParton,int *PDGcodes,AliAODEvent *aodEvent, Int_t type)
{
  TString fillthis = "";

  if(aodtrack->GetLabel() < 0) return; //reject fakes

  AliAODMCParticle *mcproton  = (AliAODMCParticle*)fAODMCEvent->GetTrack(aodtrack->GetLabel());
  if(!mcproton) return;
  if(type == 4 && mcproton->GetPdgCode() != 2212) return;
  else if(type == -4 && mcproton->GetPdgCode() != -2212) return;
  if(!mcproton->PxPyPz(Protontruth)) std::cout << "Err copying true momentum" << std::endl;


  fillthis = "fProtonPtTrueReco";
  if(fUseMCInfo && mcproton->Pt() && aodtrack->Pt())((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(mcproton->Pt(),aodtrack->Pt());


  if(Protontruth[0] != -9999.)
    {
      TLorentzVector protonVec;//needed if mass is involved
      protonVec.SetXYZM(aodtrack->Px(),aodtrack->Py(),aodtrack->Pz(),fCuts->GetHadronMasses(2212));

      Double_t ptTrue = mcproton->Pt();
      Double_t ptReco = aodtrack->Pt();
      Double_t ptTruePrime = 1./TMath::Sqrt(2.)*(ptTrue - ptReco);

      Double_t phiTrue = mcproton->Phi();
      Double_t phiReco = aodtrack->Phi();
      Double_t phiTruePrime = phiTrue - phiReco;

      Double_t thetaTrue = mcproton->Theta();
      Double_t thetaReco = aodtrack->Theta();
      Double_t thetaTruePrime = thetaTrue - thetaReco;

      fillthis = "fProtonPtResoRecoRotated";
      if(fUseMCInfo) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(ptReco,ptTruePrime);

      fillthis = "fProtonPhiResoRecoRotated";
      if(fUseMCInfo) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(ptReco,phiTruePrime);

      fillthis = "fProtonThetaResoRecoRotated";
      if(fUseMCInfo) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(ptReco,thetaTruePrime);
    }


  //check if proton is from feed-down
  if(mcproton->IsSecondaryFromWeakDecay())
    {
      //AliAODMCParticle* mcproton_motherWeak = (AliAODMCParticle*)mcarray->At(mcproton->GetMother());
      AliAODMCParticle* mcprotonMotherWeak = (AliAODMCParticle*)fAODMCEvent->GetTrack(mcproton->GetMother());

      if(!mcprotonMotherWeak->PxPyPz(ProtontruthMother)) std::cout << "Err copying true momentum" << std::endl;
      PDGcodes[0] = mcprotonMotherWeak->PdgCode();
    }
}
//________________________________________________________________________
void AliAnalysisTaskPLFemto::GetMCMomentum(AliAODv0 *v0,double *V0momtruth,double *V0momtruthMother,AliAODEvent *aodEvent)
{
  TString fillthis = "";
  //Momentum Matrix for single-particle quantities
  //MCv0 is the v0 from MC (checked for correct daughter tracks and if its not a fake)
  //v0 is the reconstructed v0

  TClonesArray *mcarray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!mcarray) return;

  Int_t daugh[2] = {2212,211};//daughter IDs of Lambda
  Int_t label = v0->MatchToMC(3122,mcarray,2,daugh);
  if(label<0) return;// only "true" if it is a real v0

  AliAODMCParticle *MCv0 = (AliAODMCParticle*)mcarray->At(label);
  if(!MCv0) return;
  if(!MCv0->PxPyPz(V0momtruth)) std::cout << "Err copying true momentum" << std::endl;


  fillthis = "fV0PtTrueReco";
  if(fUseMCInfo && MCv0->Pt() && v0->Pt())((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(MCv0->Pt(),v0->Pt());


  if(V0momtruth[0] != -9999.)
    {
      Double_t ptTrue = MCv0->Pt();
      Double_t ptReco = v0->Pt();
      Double_t ptTruePrime = 1./TMath::Sqrt(2.)*(ptTrue - ptReco);

      Double_t phiTrue = MCv0->Phi();
      Double_t phiReco = v0->Phi();
      Double_t phiTruePrime = phiTrue - phiReco;

      Double_t thetaTrue = MCv0->Theta();
      Double_t thetaReco = v0->Theta();
      Double_t thetaTruePrime = thetaTrue - thetaReco;

      fillthis = "fV0PtResoRecoRotated";
      if(fUseMCInfo) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(ptReco,ptTruePrime);

      fillthis = "fV0PhiResoRecoRotated";
      if(fUseMCInfo) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(ptReco,phiTruePrime);

      fillthis = "fV0ThetaResoRecoRotated";
      if(fUseMCInfo) ((TH2F*)(fOutputTP->FindObject(fillthis)))->Fill(ptReco,thetaTruePrime);
    }


  //check if V0 is from feed-down:
  if(MCv0->IsSecondaryFromWeakDecay())
    {
      AliAODMCParticle* MCv0Mother = (AliAODMCParticle*)mcarray->At(MCv0->GetMother());
      if(!MCv0Mother->PxPyPz(V0momtruthMother)) std::cout << "Err copying true momentum" << std::endl;
      V0momtruthMother[3] = MCv0Mother->PdgCode();
    }
}
