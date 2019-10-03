
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

//-----------------------------------------------------------------
//         AliSpectraBothEventCuts class
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#include "AliSpectraBothEventCuts.h"
#include "AliSpectraBothTrackCuts.h"
#include "AliAnalysisUtils.h"
#include "AliPPVsMultUtils.h"	
#include "AliGenPythiaEventHeader.h"
//#include "AliSpectraBothHistoManager.h"
#include <iostream>
#include "TParameter.h"
#include "AliMultSelection.h"

using namespace std;

ClassImp(AliSpectraBothEventCuts)

AliSpectraBothEventCuts::AliSpectraBothEventCuts(const char *name) : TNamed(name, "AOD Event Cuts"), fAOD(0),fAODEvent(AliSpectraBothTrackCuts::kAODobject), fTrackBits(0),fIsMC(0),fCentEstimator(""), fUseCentPatchAOD049(0), fUseSDDPatchforLHC11a(kDoNotCheckforSDD),fTriggerSettings(AliVEvent::kMB),fTrackCuts(0),
fIsSelected(0), fCentralityCutMin(0), fCentralityCutMax(0), fQVectorCutMin(0), fQVectorCutMax(0), fVertexCutMin(0), fVertexCutMax(0), fMultiplicityCutMin(0), fMultiplicityCutMax(0),fMaxChi2perNDFforVertex(0),
fMinRun(0),fMaxRun(0),fetarangeofmultiplicitycut(0.8),fUseAliPPVsMultUtils(false),
fNMCProcessType(-1),fEventMCProcessType(0),fEventMCProcessTypeIncluded(0),fchecktypeofveretxbytitle(kTRUE),fvertexselection(-1),fDotheeventcutsinmultselection(kFALSE),
fDotheBGRejection(kTRUE),fDothePileUpRejection(kTRUE),
fHistoCuts(0),fHistoVtxBefSel(0),fHistoVtxAftSel(0),fHistoEtaBefSel(0),fHistoEtaAftSel(0),fHistoNChAftSel(0),fHistoQVector(0)
,fHistoEP(0),fHistoVtxAftSelwithoutZvertexCut(0),fHistoVtxalltriggerEventswithMCz(0),fHistoVtxAftSelwithoutZvertexCutusingMCz(0),fHistoRunNumbers(0),
fHistoCentrality(0),fHistoMultiplicty(0),fAnalysisUtils(0),fAliPPVsMultUtils(0)

{
  // Constructori
 // Bool_t oldStatus = TH1::AddDirectoryStatus();
  //TH1::AddDirectory(kFALSE);	
 // fHistoCuts = new TH1I("fEventCuts", "Event Cuts", kNVtxCuts, -0.5, kNVtxCuts - 0.5);
 // fHistoVtxBefSel = new TH1F("fHistoVtxBefSel", "Vtx distr before event selection",300,-15,15);
//  fHistoVtxAftSel = new TH1F("fHistoVtxAftSel", "Vtx distr after event selection",300,-15,15);
//  fHistoVtxAftSelwithoutZvertexCut=new TH1F("fHistoVtxAftSelwithoutZvertexcut", "Vtx distr after event selection without Z vertex cut",300,-15,15);
//  fHistoVtxalltriggerEventswithMCz=new TH1F("fHistoVtxalltriggerEventswithMCz", "generated z vertex position",300,-15,15);
//  fHistoVtxAftSelwithoutZvertexCutusingMCz=new TH1F("fHistoVtxAftSelwithoutZvertexCutusingMCz", "Vtx distr after event selection without Z vertex cut using MC z",300,-15,15);
//  fHistoEtaBefSel = new TH1F("fHistoEtaBefSel", "Eta distr before event selection",500,-2,2);
//  fHistoEtaAftSel = new TH1F("fHistoEtaAftSel", "Eta distr after event selection",500,-2,2);
 // fHistoNChAftSel = new TH1F("fHistoNChAftSel", "NCh distr after event selection",2000,-0.5,1999.5);
  //fHistoQVectorPos = new TH1F("fHistoQVectorPos", "QVectorPos distribution",100,0,10);
  //fHistoQVectorNeg = new TH1F("fHistoQVectorNeg", "QVectorNeg distribution",100,0,10);
 // fHistoQVector = new TH1F("fHistoQVector", "QVector with VZERO distribution",100,0,10);
 // fHistoEP = new TH1F("fHistoEP", "EP with VZERO distribution",100,-10,10);
  fCentralityCutMin = 0.0;      // default value of centrality cut minimum, 0 ~ no cut
  fCentralityCutMax = 10000.0;  // default value of centrality cut maximum,  ~ no cut
  // fQVectorPosCutMin=0.0;
  // fQVectorPosCutMax=10000.0;
  // fQVectorNegCutMin=0.0;
  // fQVectorNegCutMax=10000.0;
  fQVectorCutMin=0.0;
  fQVectorCutMax=10000.0;
  fVertexCutMin=-10.0;
  fVertexCutMax=10.0;
  fMultiplicityCutMin=-1.0;
  fMultiplicityCutMax=-1.0;
  fTrackBits=1;
  fCentEstimator="V0M";
  fMaxChi2perNDFforVertex=-1;	
 // TH1::AddDirectory(oldStatus);	
}
//______________________________________________________

AliSpectraBothEventCuts::~AliSpectraBothEventCuts()
{
	if(fHistoCuts)
		delete fHistoCuts;
	if(fHistoVtxBefSel)
		delete fHistoVtxBefSel;
	if(fHistoVtxAftSel)
		delete fHistoVtxAftSel;
	if(fHistoVtxAftSelwithoutZvertexCut)
		delete fHistoVtxAftSelwithoutZvertexCut;
	if(fHistoVtxalltriggerEventswithMCz)
		delete fHistoVtxalltriggerEventswithMCz;
	if(fHistoVtxAftSelwithoutZvertexCutusingMCz)
		delete fHistoVtxAftSelwithoutZvertexCutusingMCz;
	if(fHistoEtaBefSel)
		delete fHistoEtaBefSel;
	if(fHistoEtaAftSel)
		delete fHistoEtaAftSel ;
	if(fHistoNChAftSel)
		delete fHistoNChAftSel;
	if(fHistoQVector)
		delete fHistoQVector;
	if(fHistoEP)
		delete fHistoEP;
	if(fHistoRunNumbers)
		delete fHistoRunNumbers;
	if(fHistoCentrality)
		delete fHistoCentrality;
  	if(fHistoMultiplicty)
		delete fHistoMultiplicty;
	if(fAnalysisUtils)
		delete fAnalysisUtils;
	if(fAliPPVsMultUtils)
		delete fAliPPVsMultUtils;
	if(fEventMCProcessType)
		delete [] fEventMCProcessType;

}
//______________________________________________________
void AliSpectraBothEventCuts::InitHisto()
{
	Bool_t oldStatus = TH1::AddDirectoryStatus();
  	TH1::AddDirectory(kFALSE);
	if(!fHistoCuts)	
 		fHistoCuts = new TH1I("fEventCuts", "Event Cuts", kNVtxCuts, -0.5, kNVtxCuts - 0.5);
 	if(!fHistoVtxBefSel )
		fHistoVtxBefSel = new TH1F("fHistoVtxBefSel", "Vtx distr before event selection",300,-15,15);
	if(!fHistoVtxAftSel)
  		fHistoVtxAftSel = new TH1F("fHistoVtxAftSel", "Vtx distr after event selection",300,-15,15);
	if(!fHistoVtxAftSelwithoutZvertexCut)
  		fHistoVtxAftSelwithoutZvertexCut=new TH1F("fHistoVtxAftSelwithoutZvertexcut", "Vtx distr after event selection without Z vertex cut",300,-15,15);
  	if(!fHistoVtxalltriggerEventswithMCz)
		fHistoVtxalltriggerEventswithMCz=new TH1F("fHistoVtxalltriggerEventswithMCz", "generated z vertex position",300,-15,15);
  	if(!fHistoVtxAftSelwithoutZvertexCutusingMCz)
		fHistoVtxAftSelwithoutZvertexCutusingMCz=new TH1F("fHistoVtxAftSelwithoutZvertexCutusingMCz", "Vtx distr after event selection without Z vertex cut using MC z",300,-15,15);
  	if(!fHistoEtaBefSel)
		fHistoEtaBefSel = new TH1F("fHistoEtaBefSel", "Eta distr before event selection",500,-2,2);
  	if(!fHistoEtaAftSel)
		fHistoEtaAftSel = new TH1F("fHistoEtaAftSel", "Eta distr after event selection",500,-2,2);
  	if(!fHistoNChAftSel)
		fHistoNChAftSel = new TH1F("fHistoNChAftSel", "NCh distr after event selection",2000,-0.5,1999.5);
  //fHistoQVectorPos = new TH1F("fHistoQVectorPos", "QVectorPos distribution",100,0,10);
  //fHistoQVectorNeg = new TH1F("fHistoQVectorNeg", "QVectorNeg distribution",100,0,10);
	if(!fHistoQVector)
		fHistoQVector = new TH1F("fHistoQVector", "QVector with VZERO distribution",100,0,10);
  	if(!fHistoEP)
		fHistoEP = new TH1F("fHistoEP", "EP with VZERO distribution",100,-10,10);
	if(!fHistoRunNumbers)
	{
		if(fMaxRun>fMinRun&&fMinRun>=0)
			fHistoRunNumbers=new TH1F("fHistoRunNumbers","Run numbers",fMaxRun-fMinRun+1,fMinRun-0.5,fMaxRun+0.5);
		else
			fHistoRunNumbers=new TH1F("fHistoRunNumbers","Run numbers",1001,120000-.5,121000+0.5);

	}
	if(!fHistoCentrality)
		fHistoCentrality = new TH2F("fHistoCentrality", "centrality",2,0,2,100,0.0,100);

	if(!fHistoMultiplicty)
		fHistoMultiplicty= new TH2F("fHistoMultiplicty", "multiplicty estimator",2,0,2,155,-4.5,150.5);

	TH1::AddDirectory(oldStatus);		
}
//______________________________________________________
Bool_t AliSpectraBothEventCuts::IsSelected(AliVEvent * aod,AliSpectraBothTrackCuts* trackcuts,Bool_t isMC,Double_t mcZ,TH1F* managerhisteventcuts)
{
  // Returns true if Event Cuts are selected and applied
	fIsSelected =kFALSE;
  	fAOD = aod;
	if(!CheckifESDorAODEvent())
		return kFALSE;
  	fTrackCuts = trackcuts;
  	fHistoCuts->Fill(kProcessedEvents);
  	if(managerhisteventcuts)
		managerhisteventcuts->Fill(0);
	  fHistoRunNumbers->Fill(aod->GetRunNumber());
  	Bool_t IsPhysSel = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerSettings);//FIXME we can add the trigger mask here
	if(!IsPhysSel)
		return IsPhysSel;
	  if(fAnalysisUtils&&(!fDotheeventcutsinmultselection)) // we check for pile-up
	  {	
		if(fDothePileUpRejection)
		{
		 	IsPhysSel = (!fAnalysisUtils->IsPileUpEvent(fAOD));

		}
		if(fDotheBGRejection)
		{
			IsPhysSel = (!fAnalysisUtils->IsSPDClusterVsTrackletBG(fAOD));

		}
	}
 	 if(!IsPhysSel)
		return IsPhysSel;
 	
 	 if(isMC)	
   		fHistoVtxalltriggerEventswithMCz->Fill(mcZ);
   	 if(fUseSDDPatchforLHC11a!=kDoNotCheckforSDD)
		if(!CheckSDDPatchforLHC11a())
			return kFALSE;

    	fHistoCuts->Fill(kPhysSelEvents);
   	if(managerhisteventcuts)
		managerhisteventcuts->Fill(1);


	if(fDotheeventcutsinmultselection)
	{
		const AliVVertex * vertex = fAOD->GetPrimaryVertex();
 	 	if(vertex)
			fHistoVtxBefSel->Fill(vertex->GetZ());

		if( CheckCentralityCut() && CheckMultiplicityCut() && CheckQVectorCut())
			fIsSelected=kTRUE;
		else
			fIsSelected=kFALSE;	
		
 	 	
 	 	if(fIsSelected&&vertex)
 	 	{
  			fHistoVtxAftSelwithoutZvertexCut->Fill(vertex->GetZ());
      			if(isMC)
        			fHistoVtxAftSelwithoutZvertexCutusingMCz->Fill(mcZ);	
     			if (vertex->GetZ() > fVertexCutMin && vertex->GetZ() < fVertexCutMax)
     			{
				fHistoVtxAftSel->Fill(vertex->GetZ());
     			}
  		}
	}
	else
	{
   		const AliVVertex * vertex = fAOD->GetPrimaryVertex();//FIXME vertex is recreated	
 	 	if(vertex)
			fHistoVtxBefSel->Fill(vertex->GetZ());
 	 	if(CheckVtx())
 	 	{ //selection on vertex
 	 	  	fIsSelected=kTRUE;
 	 	}
 	 	if(fIsSelected&&vertex)
 	 	{
  			fHistoVtxAftSelwithoutZvertexCut->Fill(vertex->GetZ());
      			if(isMC)
        			fHistoVtxAftSelwithoutZvertexCutusingMCz->Fill(mcZ);	
     			if (vertex->GetZ() > fVertexCutMin && vertex->GetZ() < fVertexCutMax)
     			{
				fIsSelected=kTRUE;
				fHistoVtxAftSel->Fill(vertex->GetZ());
     			}
    			else	
    			{
				fIsSelected=kFALSE;
    			}	
  		}
  		if(fIsSelected)
  		{
			if( CheckCentralityCut() && CheckMultiplicityCut() && CheckQVectorCut())
				fIsSelected=kTRUE;
			else
				fIsSelected=kFALSE;	
  		}	
	}


  	Int_t Nch=0;
  	for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) 
  	{
   		 AliVTrack* track =dynamic_cast<AliVTrack*>(fAOD->GetTrack(iTracks));
   	/* if(fAODEvent==AliSpectraBothTrackCuts::kESDobject)
			track=dynamic_cast<AliVTrack*>(esdevent->GetTrack(iTracks));
		     else if (fAODEvent==AliSpectraBothTrackCuts::kAODobject)
			track=dynamic_cast<AliVTrack*>(aodevent->GetTrack(iTracks));
     	else return false;*/
     
   		 if (!fTrackCuts->IsSelected(track,kFALSE)) 
			continue;
  	  	fHistoEtaBefSel->Fill(track->Eta());
  	 	 if(fIsSelected)
		{
      			fHistoEtaAftSel->Fill(track->Eta());
      			Nch++;
    		}
  	}
  	if(fIsSelected)
  	{	
  		fHistoNChAftSel->Fill(Nch);
		fHistoCuts->Fill(kAcceptedEvents);
		if(managerhisteventcuts)
			managerhisteventcuts->Fill(2);

	}	
  	return fIsSelected;
}

//______________________________________________________
Bool_t AliSpectraBothEventCuts::CheckVtx()
{
  // reject events outside of range
  	const AliVVertex * vertex = fAOD->GetPrimaryVertex();
	Int_t vertexflag=0;		
  	if(fAODEvent==AliSpectraBothTrackCuts::kESDobject)
 	{
		AliESDEvent* esdevent=dynamic_cast<AliESDEvent*>(fAOD);
		if(!esdevent)
			return kFALSE;
		const AliESDVertex* vertextracks=esdevent->GetPrimaryVertexTracks();
		const AliESDVertex* vertexspd=esdevent->GetPrimaryVertexSPD();
		if(vertextracks)
		{	
			if(vertextracks->GetNContributors()>0)
				vertexflag=vertexflag+1;
		}
		if(vertexspd)
		{
			if(vertexspd->GetNContributors()>0)
				vertexflag=vertexflag+2;
		}

 	}			
 	else if(fAODEvent==AliSpectraBothTrackCuts::kAODobject)
  	{
		AliAODEvent* aodevent=dynamic_cast<AliAODEvent*>(fAOD);
		if(!aodevent)
			return kFALSE;
		AliAODVertex* vertexspd=aodevent->GetPrimaryVertexSPD();
		if(vertexspd)
		{
			if(vertexspd->GetNContributors()>0)
				vertexflag=vertexflag+2;
			TString tmp(vertex->GetTitle());
			if(tmp.Contains("VertexerTracksWithConstraint"))
				vertexflag=vertexflag+1;
		}
		

  	}
  	else 
		return kFALSE;	 
       			
  	//when moving to 2011 wìone has to add a cut using SPD vertex.
 	 //The point is that for events with |z|>20 the vertexer tracks is not working (only 2011!). One has to put a safety cut using SPD vertex large e.g. 15cm
 	if (!vertex)
 	{
      		fHistoCuts->Fill(kVtxNoEvent);
      		return kFALSE;
   	 }

    	if(vertex->GetNContributors()<1)
   	{
      		fHistoCuts->Fill(kZeroCont);
      		return kFALSE;
	}
	
   	TString tmp(vertex->GetTitle());
	if(fchecktypeofveretxbytitle)
	{
   		if(tmp.Contains("NoConstraint"))
   		{
   		     fHistoCuts->Fill(kTPCasPV);
   		     return kFALSE;
   		}
	}
	else
	{
		if(vertexflag==0)
   		{
   		     fHistoCuts->Fill(kTPCasPV);
   		     return kFALSE;
   		}
	
	}
	if(fvertexselection>0&&((vertexflag&fvertexselection)==0))
		return kFALSE;


	// if (vertex->GetZ() > fVertexCutMin && vertex->GetZ() < fVertexCutMax)
   // {
    //  return kTRUE;
   // }
     if(!CheckVtxChi2perNDF())	  
		return kFALSE;
		
  fHistoCuts->Fill(kGoodVtx);
  //return kFALSE;
   return kTRUE;
}

//______________________________________________________
Bool_t AliSpectraBothEventCuts::CheckCentralityCut()
{
  // Check centrality cut
	 if ( fCentralityCutMax<0.0  &&  fCentralityCutMin<0.0 )  
		return kTRUE;
	  Double_t cent=0;
	  Bool_t validcent=kFALSE;
	  AliMultSelection *multselection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
	  if(multselection)
	  {
 		cent = multselection->GetMultiplicityPercentile(fCentEstimator.Data(),fDotheeventcutsinmultselection);
		validcent=kTRUE;
	  }
          else
	  {	
	  	if(fAODEvent==AliSpectraBothTrackCuts::kESDobject)
	  	{
			TString par1name("cent");
			par1name+=fCentEstimator;
			TParameter<Double_t>* par= dynamic_cast<TParameter<Double_t>*>(fAOD->FindListObject(par1name.Data()));
			if(par)
			{
				validcent=kTRUE;
				cent=par->GetVal();
			}

  		}
  		if(!validcent)
  		{			
  			if(fUseAliPPVsMultUtils)
  			{
				if(!fAliPPVsMultUtils)
					fAliPPVsMultUtils=new AliPPVsMultUtils();
				cent=fAliPPVsMultUtils->GetMultiplicityPercentile(fAOD,fCentEstimator.Data());
  			}	
  			else
  			{
  				if(!fUseCentPatchAOD049)
					cent=fAOD->GetCentrality()->GetCentralityPercentile(fCentEstimator.Data());
  				else 
					cent=ApplyCentralityPatchAOD049();
  			}
  		}
	}	
  	fHistoCentrality->Fill(0.5,cent);	
  	if ( (cent < fCentralityCutMax)  &&  (cent >= fCentralityCutMin) )  
  	{
		 fHistoCentrality->Fill(1.5,cent);	
  		return kTRUE;
  	}	   
  	fHistoCuts->Fill(kVtxCentral);

  	return kFALSE;
}

//______________________________________________________
Bool_t AliSpectraBothEventCuts::CheckMultiplicityCut()
{
  // Check multiplicity cut
	if(fMultiplicityCutMin<0 && fMultiplicityCutMax<0)
		return kTRUE;
	Int_t Ncharged=-1;
	if(fAODEvent==AliSpectraBothTrackCuts::kESDobject)
	{
	
		TParameter<Int_t>* par=0;
		if(TMath::Abs(0.8-fetarangeofmultiplicitycut)<0.1)
			par= dynamic_cast<TParameter<Int_t>*>(fAOD->FindListObject("Ncheta0dot8"));
		else  
			par= dynamic_cast<TParameter<Int_t>*>(fAOD->FindListObject("Ncheta0dot5"));	
		if(par)
		{
			Ncharged=par->GetVal();
		}	
		else
		{
			AliESDEvent* esdevent=dynamic_cast<AliESDEvent*>(fAOD);
			if(!esdevent)
			{
				AliFatal("Not a standard ESD");
				return kFALSE;

			}
			AliESDtrackCuts::MultEstTrackType estType = esdevent->GetPrimaryVertexTracks()->GetStatus() ? AliESDtrackCuts::kTrackletsITSTPC : AliESDtrackCuts::kTracklets;
			Ncharged=AliESDtrackCuts::GetReferenceMultiplicity(esdevent,estType,fetarangeofmultiplicitycut);
		}
	}
	else if(fAODEvent==AliSpectraBothTrackCuts::kAODobject)
	{
		AliAODEvent* aodevent=0x0;
		aodevent=dynamic_cast<AliAODEvent*>(fAOD);
		if(!aodevent)
		{
			AliFatal("Not a standard AOD");
			return kFALSE;
		}
                AliAODHeader * header = dynamic_cast<AliAODHeader*>(aodevent->GetHeader());
                if(!header)
		{ 
			AliFatal("Not a standard AOD");
			return kFALSE;
		}

		if(TMath::Abs(0.8-fetarangeofmultiplicitycut)<0.1)
			Ncharged=header->GetRefMultiplicityComb08();
		else if (TMath::Abs(0.5-fetarangeofmultiplicitycut)<0.1)
			Ncharged=header->GetRefMultiplicityComb05();
		else 
			Ncharged=-1;
	}
	else
		return kFALSE;	 

   	fHistoMultiplicty->Fill(0.5,Ncharged);
   	if(Ncharged>=fMultiplicityCutMin && Ncharged<fMultiplicityCutMax)
   	{ 
		fHistoMultiplicty->Fill(1.5,Ncharged);
	  	return kTRUE;
	 }
	fHistoCuts->Fill(kVtxCentral);
	return kFALSE;
}

//______________________________________________________
Bool_t AliSpectraBothEventCuts::CheckQVectorCut()
{
	 if(fQVectorCutMin<0.0 && fQVectorCutMax<0.0)
		return kTRUE;
  // Check qvector cut
  /// FIXME: Q vector
  // //Selection on QVector, before ANY other selection on the event
  // //Spectra MUST be normalized wrt events AFTER the selection on Qvector
  // Double_t Qx2EtaPos = 0, Qy2EtaPos = 0;
  // Double_t Qx2EtaNeg = 0, Qy2EtaNeg = 0;
 
  // Int_t multPos = 0;
  // Int_t multNeg = 0;
  // for(Int_t iT = 0; iT < fAOD->GetNumberOfTracks(); iT++) {
  //   AliAODTrack* aodTrack = fAOD->GetTrack(iT);
  //   if (!aodTrack->TestFilterBit(fTrackBits)) continue;
  //   if (aodTrack->Eta() >= 0){
  //     multPos++;
  //     Qx2EtaPos += TMath::Cos(2*aodTrack->Phi()); 
  //     Qy2EtaPos += TMath::Sin(2*aodTrack->Phi());
  //   } else {
  //     multNeg++;
  //     Qx2EtaNeg += TMath::Cos(2*aodTrack->Phi()); 
  //     Qy2EtaNeg += TMath::Sin(2*aodTrack->Phi());
  //   }
  // } 
  // Double_t qPos=-999;
  // if(multPos!=0)qPos= TMath::Sqrt((Qx2EtaPos*Qx2EtaPos + Qy2EtaPos*Qy2EtaPos)/multPos);
  // Double_t qNeg=-999;
  // if(multNeg!=0)qNeg= TMath::Sqrt((Qx2EtaNeg*Qx2EtaNeg + Qy2EtaNeg*Qy2EtaNeg)/multNeg);
  //if(qPos<fQVectorPosCutMin || qPos>fQVectorPosCutMax || qNeg<fQVectorNegCutMin || qNeg>fQVectorNegCutMax)return kFALSE;
  
  Double_t qxEPVZERO = 0, qyEPVZERO = 0;
  Double_t qVZERO = -999;
  Double_t psi=fAOD->GetEventplane()->CalculateVZEROEventPlane(fAOD,10,2,qxEPVZERO,qyEPVZERO);
  
  qVZERO= TMath::Sqrt(qxEPVZERO*qxEPVZERO + qyEPVZERO*qyEPVZERO);
  if(qVZERO<fQVectorCutMin || qVZERO>fQVectorCutMax)return kFALSE;
  fHistoQVector->Fill(qVZERO);
  fHistoEP->Fill(psi);
  
  fHistoCuts->Fill(kQVector);
  // fHistoQVectorPos->Fill(qPos);
  // fHistoQVectorNeg->Fill(qNeg);
  return kTRUE;
}
//____________________________________________________________
 Bool_t AliSpectraBothEventCuts::CheckVtxChi2perNDF()
 {
	if(fMaxChi2perNDFforVertex<0)
		return kTRUE;
	 const AliVVertex * vertex = fAOD->GetPrimaryVertex();
	if(TMath::Abs(vertex->GetChi2perNDF())>fMaxChi2perNDFforVertex) 
		return kFALSE;
	return kTRUE;
 }



//______________________________________________________
void AliSpectraBothEventCuts::PrintCuts()
{
  // print info about event cuts
  cout << "Event Stats" << endl;
  if(fHistoCuts)
  {		
 	 cout << " > Number of accepted events: " << fHistoCuts->GetBinContent(kAcceptedEvents + 1) << endl;
 	 cout << " > Number of processed events: " << fHistoCuts->GetBinContent(kProcessedEvents + 1) << endl;
 	 cout << " > Number of PhysSel events: " << fHistoCuts->GetBinContent(kPhysSelEvents + 1) << endl;
 	 cout << " > With good veretx: " << fHistoCuts->GetBinContent(kGoodVtx + 1) << endl;
 	 cout << " > Events cut by centrality: " << fHistoCuts->GetBinContent(kVtxCentral + 1) << endl;
 	 cout << " > Events without vertex: " << fHistoCuts->GetBinContent(kVtxNoEvent + 1) << endl;
	  cout << " > QVector cut: " << fHistoCuts->GetBinContent(kQVector + 1) << endl;
  }	
  cout << " > Track type used for the QVector calculation: " << fTrackBits << endl;
  // cout << " > QPosRange: [" << fQVectorPosCutMin <<"," <<fQVectorPosCutMax<<"]"<< endl;
  // cout << " > QNegRange: [" << fQVectorNegCutMin <<"," <<fQVectorNegCutMax<<"]"<< endl;
  cout << " > QRange: [" << fQVectorCutMin <<"," <<fQVectorCutMax<<"]"<< endl;
  cout << " > Vertex: [" << fVertexCutMin <<"," <<fVertexCutMax<<"]"<< endl;
  cout << " > Multiplicity: [" << fMultiplicityCutMin <<"," <<fMultiplicityCutMax<<"]"<< endl;
  cout << " > Centrality: [" << fCentralityCutMin <<"," <<fCentralityCutMax<<"]"<< endl;
}
//______________________________________________________

Double_t AliSpectraBothEventCuts::ApplyCentralityPatchAOD049()
{
   //
   //Apply centrality patch for AOD049 outliers
   //
  // if (fCentralityType!="V0M") {
  //   AliWarning("Requested patch forAOD049 for wrong value (not centrality from V0).");
  //   return -999.0;
  // }
  AliCentrality *centrality = fAOD->GetCentrality();
   if (!centrality) {
     AliWarning("Cannot get centrality from AOD event.");
     return -999.0;
   }
   
   Float_t cent = (Float_t)(centrality->GetCentralityPercentile("V0M"));
   /*
     Bool_t isSelRun = kFALSE;
     Int_t selRun[5] = {138364, 138826, 138828, 138836, 138871};
     if(cent<0){
     Int_t quality = centrality->GetQuality();
     if(quality<=1){
     cent=(Float_t)centrality->GetCentralityPercentileUnchecked("V0M");
     } else {
     Int_t runnum=aodEvent->GetRunNumber();
     for(Int_t ir=0;ir<5;ir++){
     if(runnum==selRun[ir]){
     isSelRun=kTRUE;
     break;
     }
     }
     if((quality==8||quality==9)&&isSelRun) cent=(Float_t)centrality->GetCentralityPercentileUnchecked("V0M");
     }
     }
   */
   if(cent>=0.0) {
     Float_t v0 = 0.0;
     AliAODEvent *aodEvent = (AliAODEvent *)fAOD;
     AliAODVZERO *aodV0 = (AliAODVZERO *) aodEvent->GetVZEROData();
     v0+=aodV0->GetMTotV0A();
     v0+=aodV0->GetMTotV0C();
     if ( (cent==0) && (v0<19500) ) {
       AliDebug(3, Form("Filtering issue in centrality -> cent = %5.2f",cent));
       return -999.0;
     }
     Float_t tkl = (Float_t)(aodEvent->GetTracklets()->GetNumberOfTracklets());
     Float_t val = 1.30552 +  0.147931 * v0;
     
     Float_t tklSigma[101] = {176.644, 156.401, 153.789, 153.015, 142.476, 137.951, 136.127, 129.852, 127.436, 124.86,
			      120.788, 115.611, 113.172, 110.496, 109.127, 104.421, 102.479, 99.9766, 97.5152, 94.0654,
			      92.4602, 89.3364, 87.1342, 83.3497, 82.6216, 81.1084, 78.0793, 76.1234, 72.9434, 72.1334,
			      68.0056, 68.2755, 66.0376, 62.9666, 62.4274, 59.65, 58.3776, 56.6361, 54.5184, 53.4224,
			      51.932, 50.8922, 48.2848, 47.912, 46.5717, 43.4114, 43.2083, 41.3065, 40.1863, 38.5255,
			      37.2851, 37.5396, 34.4949, 33.8366, 31.8043, 31.7412, 30.8392, 30.0274, 28.8793, 27.6398,
			      26.6488, 25.0183, 25.1489, 24.4185, 22.9107, 21.2002, 21.6977, 20.1242, 20.4963, 19.0235,
			      19.298, 17.4103, 16.868, 15.2939, 15.2939, 16.0295, 14.186, 14.186, 15.2173, 12.9504, 12.9504,
			      12.9504, 15.264, 12.3674, 12.3674, 12.3674, 12.3674, 12.3674, 18.3811, 13.7544, 13.7544,
			      13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544
     };
     
     if ( TMath::Abs(tkl-val) > 6.*tklSigma[(Int_t)cent] )  {
       AliDebug(3, Form("Outlier event in centrality -> cent = %5.2f",cent));
       return -999.0;
     }
   } else {
     //force it to be -999. whatever the negative value was
     cent = -999.;
   }
   return cent;
}

//______________________________________________________


Long64_t AliSpectraBothEventCuts::Merge(TCollection* list)
{
  // Merge a list of AliSpectraBothEventCuts objects with this.
  // Returns the number of merged objects (including this).

  //  AliInfo("Merging");

  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of all histograms
  TList collections;//FIXME we should use only 1 collection
  TList collections_histoVtxBefSel;
  TList collections_histoVtxAftSel;
  TList collections_histoEtaBefSel;
  TList collections_histoEtaAftSel;
  TList collections_histoNChAftSel;
  // TList collections_histoQVectorPos;
  // TList collections_histoQVectorNeg;
  TList collections_histoQVector;
  TList collections_histoEP;
  TList collections_histoVtxAftSelwithoutZvertexCut;
  TList collections_histoVtxalltriggerEventswithMCz;
  TList collections_histoVtxAftSelwithoutZvertexCutusingMCz;			
  TList collections_histoRunNumbers;
  TList collections_histoCentrality;			
  TList collections_histoMultiplicty;

  Int_t count = 0;

  while ((obj = iter->Next())) {
    AliSpectraBothEventCuts* entry = dynamic_cast<AliSpectraBothEventCuts*> (obj);
    if (entry == 0) 
      continue;

    TH1I * histo = entry->GetHistoCuts(); 
    if(histo)    
    	collections.Add(histo);
    TH1F * histo_histoVtxBefSel = entry->GetHistoVtxBefSel();      
    if(histo_histoVtxBefSel)
    	collections_histoVtxBefSel.Add(histo_histoVtxBefSel);
    TH1F * histo_histoVtxAftSel = entry->GetHistoVtxAftSel();      
    if(histo_histoVtxAftSel)
    	collections_histoVtxAftSel.Add(histo_histoVtxAftSel);
    TH1F * histo_histoEtaBefSel = entry->GetHistoEtaBefSel(); 
    if(histo_histoEtaBefSel)	     
    	collections_histoEtaBefSel.Add(histo_histoEtaBefSel);
    TH1F * histo_histoEtaAftSel = entry->GetHistoEtaAftSel();
    if(histo_histoEtaAftSel)	      
    	collections_histoEtaAftSel.Add(histo_histoEtaAftSel);
    TH1F * histo_histoNChAftSel = entry->GetHistoNChAftSel();
    if(histo_histoNChAftSel)	      
    	collections_histoNChAftSel.Add(histo_histoNChAftSel);
    // TH1F * histo_histoQVectorPos = entry->GetHistoQVectorPos();      
    // collections_histoQVectorPos.Add(histo_histoQVectorPos);
    // TH1F * histo_histoQVectorNeg = entry->GetHistoQVectorNeg();      
    // collections_histoQVectorNeg.Add(histo_histoQVectorNeg);	
    TH1F * histo_histoQVector = entry->GetHistoQVector();     
    if(histo_histoQVector)	 
    	collections_histoQVector.Add(histo_histoQVector);
    TH1F * histo_histoEP = entry->GetHistoEP();      
    if(histo_histoEP) 
    	collections_histoEP.Add(histo_histoEP);
    TH1F* histo_histoVtxAftSelwithoutZvertexCut=entry->GetHistoVtxAftSelwithoutZvertexCut();
     if(histo_histoVtxAftSelwithoutZvertexCut)
    	collections_histoVtxAftSelwithoutZvertexCut.Add(histo_histoVtxAftSelwithoutZvertexCut);
    TH1F* histo_histoVtxalltriggerEventswithMCz=entry->GetHistoVtxGenerated();
     if(histo_histoVtxalltriggerEventswithMCz)
    	collections_histoVtxalltriggerEventswithMCz.Add(histo_histoVtxalltriggerEventswithMCz);
    
   TH1F* histo_histoVtxAftSelwithoutZvertexCutusingMCz=entry->GetHistoVtxAftSelwithoutZvertexCutusingMCz();
     if(histo_histoVtxAftSelwithoutZvertexCutusingMCz)	
    	collections_histoVtxAftSelwithoutZvertexCutusingMCz.Add(histo_histoVtxAftSelwithoutZvertexCutusingMCz);	
    
    TH1F* histo_histoRunNumbers=entry->GetHistoRunNumbers();
    if(histo_histoRunNumbers)
	collections_histoRunNumbers.Add(histo_histoRunNumbers);
 
   TH2F* histo_histoCentrality=entry->GetHistoCentrality();
  if(histo_histoCentrality)
 	 collections_histoCentrality.Add(histo_histoCentrality);	 	

TH2F* histo_histoMultiplicty=entry->GetHistoMultiplicty();
  if(histo_histoMultiplicty)
  	collections_histoMultiplicty.Add(histo_histoMultiplicty);


    count++;
  }
  if(fHistoCuts)
  	fHistoCuts->Merge(&collections);
  if(fHistoVtxBefSel)
  	fHistoVtxBefSel->Merge(&collections_histoVtxBefSel);
  if(fHistoVtxAftSel)
  	fHistoVtxAftSel->Merge(&collections_histoVtxAftSel);
  if(fHistoEtaBefSel)
 	 fHistoEtaBefSel->Merge(&collections_histoEtaBefSel);
  if(fHistoEtaAftSel)
  	fHistoEtaAftSel->Merge(&collections_histoEtaAftSel);
  if(fHistoNChAftSel)
  	fHistoNChAftSel->Merge(&collections_histoNChAftSel);
  // fHistoQVectorPos->Merge(&collections_histoQVectorPos);
  // fHistoQVectorNeg->Merge(&collections_histoQVectorNeg);
  if(fHistoQVector)
 	 fHistoQVector->Merge(&collections_histoQVector);
  if(fHistoEP)
  	fHistoEP->Merge(&collections_histoEP);
  if(fHistoVtxAftSelwithoutZvertexCut)
  	fHistoVtxAftSelwithoutZvertexCut->Merge(&collections_histoVtxAftSelwithoutZvertexCut);
  if(fHistoVtxalltriggerEventswithMCz)
 	 fHistoVtxalltriggerEventswithMCz->Merge(&collections_histoVtxalltriggerEventswithMCz);
  if(fHistoVtxAftSelwithoutZvertexCutusingMCz)	
  	fHistoVtxAftSelwithoutZvertexCutusingMCz->Merge(&collections_histoVtxAftSelwithoutZvertexCutusingMCz);
  if(fHistoRunNumbers)
  	fHistoRunNumbers->Merge(&collections_histoRunNumbers);
  if(fHistoCentrality)
 	 fHistoCentrality->Merge(&collections_histoCentrality);
  if(fHistoMultiplicty)
  	fHistoMultiplicty->Merge(&collections_histoMultiplicty);


  delete iter;

  return count+1;
}
//__________________________________________________________________________________________________________
void AliSpectraBothEventCuts::SetRunNumberRange(Int_t min, Int_t max)
{
	if(max>min&&min>=0)
	{
		 fMinRun=min; 			 
  		 fMaxRun=max;
	}
}
//__________________________________________________________________________________________________________
Bool_t AliSpectraBothEventCuts::CheckMCProcessType(AliMCEvent* mcevent)
{
	if(fNMCProcessType<=0)
		return kTRUE;
	if(!mcevent)
		return kFALSE;
	AliHeader* aHeader=mcevent->Header();
	if(!aHeader)
		return kFALSE;
	AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(aHeader->GenEventHeader());
	if(!pythiaGenHeader)
		return kFALSE;	
	Int_t processtype=pythiaGenHeader->ProcessType();
	for(int i=0;i<fNMCProcessType;i++)
	{
		if(fEventMCProcessType[i]<0)
			continue;
		if (processtype==fEventMCProcessType[i])
			return fEventMCProcessTypeIncluded;
	}
	return !fEventMCProcessTypeIncluded;
}
//_____________________________________________________________________________________________________________
void AliSpectraBothEventCuts::SetNMCProcessType(Int_t flag) 
{
	fNMCProcessType=flag;
	if(fEventMCProcessType)
		delete [] fEventMCProcessType;
	
	fEventMCProcessType= new Int_t[fNMCProcessType];
	if(!fEventMCProcessType)
	{
		fNMCProcessType=-1;
		return;
	}
	for(int i=0;i<fNMCProcessType;i++)
		fEventMCProcessType[i]=-1.0;
	return;
}
//________________________________________________________________________________________________________________________
void AliSpectraBothEventCuts::AddMCProcessType(Int_t type,Int_t index)
{
	if(index<fNMCProcessType)
		fEventMCProcessType[index]=type;
}
//__________________________________________________________________________________________________________________________
Bool_t AliSpectraBothEventCuts::CheckifESDorAODEvent()
{
	TString nameoftrack(fAOD->ClassName());  
    	if(!nameoftrack.CompareTo("AliESDEvent"))
    	{
		fAODEvent=AliSpectraBothTrackCuts::kESDobject;	
		return true;
	}
	else if(!nameoftrack.CompareTo("AliAODEvent"))
	{
		fAODEvent=AliSpectraBothTrackCuts::kAODobject;
		return true;
	}
	else
		return false; 
}
//______________________________________________________________________________________________________________________________
Bool_t AliSpectraBothEventCuts::CheckSDDPatchforLHC11a()
{
	Bool_t isSDD=kFALSE;
	if(fAODEvent==AliSpectraBothTrackCuts::kESDobject)
	{
		AliESDEvent* esdevent=0x0;
		esdevent=dynamic_cast<AliESDEvent*>(fAOD);
		if(!esdevent)
			return kFALSE;
		if(esdevent->GetFiredTriggerClasses().Contains("ALLNOTRD"))
			isSDD=kTRUE;
	}		
	else if(fAODEvent==AliSpectraBothTrackCuts::kAODobject)
	{	
		AliAODEvent* aodevent=0x0;
		aodevent=dynamic_cast<AliAODEvent*>(fAOD);
		if(!aodevent)
			return kFALSE;
		if(aodevent->GetFiredTriggerClasses().Contains("ALLNOTRD"))
			isSDD=kTRUE;	
	}	
	else
		return false;
      	if(fUseSDDPatchforLHC11a==kwithSDD&&isSDD==kFALSE)
		return false;
    	if(fUseSDDPatchforLHC11a==kwithoutSDD&&isSDD==kTRUE)
		return false;
	return kTRUE;
}
