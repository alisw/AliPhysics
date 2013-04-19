
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
//#include "AliSpectraBothHistoManager.h"
#include <iostream>

using namespace std;

ClassImp(AliSpectraBothEventCuts)

AliSpectraBothEventCuts::AliSpectraBothEventCuts(const char *name) : TNamed(name, "AOD Event Cuts"), fAOD(0),fAODEvent(AliSpectraBothTrackCuts::kAODobject), fTrackBits(0),fIsMC(0),fCentEstimator(""), fUseCentPatchAOD049(0), fUseSDDPatchforLHC11a(kDoNotCheckforSDD),fTriggerSettings(AliVEvent::kMB),fTrackCuts(0),
fIsSelected(0), fCentralityCutMin(0), fCentralityCutMax(0), fQVectorCutMin(0), fQVectorCutMax(0), fVertexCutMin(0), fVertexCutMax(0), fMultiplicityCutMin(0), fMultiplicityCutMax(0),fMaxChi2perNDFforVertex(0),
fHistoCuts(0),fHistoVtxBefSel(0),fHistoVtxAftSel(0),fHistoEtaBefSel(0),fHistoEtaAftSel(0),fHistoNChAftSel(0),fHistoQVector(0)
,fHistoEP(0),fHistoVtxAftSelwithoutZvertexCut(0),fHistoVtxalltriggerEventswithMCz(0),fHistoVtxAftSelwithoutZvertexCutusingMCz(0)
{
  // Constructor
  fHistoCuts = new TH1I("fEventCuts", "Event Cuts", kNVtxCuts, -0.5, kNVtxCuts - 0.5);
  fHistoVtxBefSel = new TH1F("fHistoVtxBefSel", "Vtx distr before event selection",300,-15,15);
  fHistoVtxAftSel = new TH1F("fHistoVtxAftSel", "Vtx distr after event selection",300,-15,15);
  fHistoVtxAftSelwithoutZvertexCut=new TH1F("fHistoVtxAftSelwithoutZvertexcut", "Vtx distr after event selection without Z vertex cut",300,-15,15);
  fHistoVtxalltriggerEventswithMCz=new TH1F("fHistoVtxalltriggerEventswithMCz", "generated z vertex position",300,-15,15);
  fHistoVtxAftSelwithoutZvertexCutusingMCz=new TH1F("fHistoVtxAftSelwithoutZvertexCutusingMCz", "Vtx distr after event selection without Z vertex cut using MC z",300,-15,15);
  fHistoEtaBefSel = new TH1F("fHistoEtaBefSel", "Eta distr before event selection",500,-2,2);
  fHistoEtaAftSel = new TH1F("fHistoEtaAftSel", "Eta distr after event selection",500,-2,2);
  fHistoNChAftSel = new TH1F("fHistoNChAftSel", "NCh distr after event selection",2000,-0.5,1999.5);
  //fHistoQVectorPos = new TH1F("fHistoQVectorPos", "QVectorPos distribution",100,0,10);
  //fHistoQVectorNeg = new TH1F("fHistoQVectorNeg", "QVectorNeg distribution",100,0,10);
  fHistoQVector = new TH1F("fHistoQVector", "QVector with VZERO distribution",100,0,10);
  fHistoEP = new TH1F("fHistoEP", "EP with VZERO distribution",100,-10,10);
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
}

//______________________________________________________
Bool_t AliSpectraBothEventCuts::IsSelected(AliVEvent * aod,AliSpectraBothTrackCuts* trackcuts,Bool_t isMC,Double_t mcZ)
{
  // Returns true if Event Cuts are selected and applied
  fAOD = aod;
  fTrackCuts = trackcuts;
  fHistoCuts->Fill(kProcessedEvents);

  Bool_t IsPhysSel = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerSettings);//FIXME we can add the trigger mask here
  if(!IsPhysSel)return IsPhysSel;

  if(isMC)	
   	fHistoVtxalltriggerEventswithMCz->Fill(mcZ);

   //loop on tracks, before event selection, filling QA histos
 AliESDEvent* esdevent=0x0;
  AliAODEvent* aodevent=0x0;
  Bool_t isSDD=kFALSE;
     TString nameoftrack(fAOD->ClassName());  
    if(!nameoftrack.CompareTo("AliESDEvent"))
    {
		fAODEvent=AliSpectraBothTrackCuts::kESDobject;
		
		if(fUseSDDPatchforLHC11a!=kDoNotCheckforSDD)
		{
			esdevent=dynamic_cast<AliESDEvent*>(fAOD);
			if(!esdevent)
				return kFALSE;
			if(esdevent->GetFiredTriggerClasses().Contains("ALLNOTRD"))
				isSDD=kTRUE;
		}	
	}
	else if(!nameoftrack.CompareTo("AliAODEvent"))
	{
		fAODEvent=AliSpectraBothTrackCuts::kAODobject;
		if(fUseSDDPatchforLHC11a!=kDoNotCheckforSDD)
		{
			aodevent=dynamic_cast<AliAODEvent*>(fAOD);
			if(!aodevent)
				return kFALSE;
			if(aodevent->GetFiredTriggerClasses().Contains("ALLNOTRD"))
				isSDD=kTRUE;	
		}	
	}
	else
		return false; 
      if(fUseSDDPatchforLHC11a==kwithSDD&&isSDD==kFALSE)
		return false;
     if(fUseSDDPatchforLHC11a==kwithoutSDD&&isSDD==kTRUE)
		return false;


	
    fHistoCuts->Fill(kPhysSelEvents);




   const AliVVertex * vertex = fAOD->GetPrimaryVertex();//FIXME vertex is recreated	

  if(vertex)fHistoVtxBefSel->Fill(vertex->GetZ());
  fIsSelected =kFALSE;
  if(CheckVtx() && CheckCentralityCut() && CheckMultiplicityCut() && CheckVtxChi2perNDF())
   { //selection on vertex and Centrality

    fIsSelected=CheckQVectorCut(); // QVector is calculated only if the centrality and vertex are correct (performance)
  }
  if(fIsSelected&&vertex)
 {
      fHistoVtxAftSelwithoutZvertexCut->Fill(vertex->GetZ());
      if(isMC)
          fHistoVtxAftSelwithoutZvertexCutusingMCz->Fill(mcZ);	
     if (vertex->GetZ() > fVertexCutMin && vertex->GetZ() < fVertexCutMax)
     {
      		fHistoCuts->Fill(kAcceptedEvents);
		fIsSelected=kTRUE;
		fHistoVtxAftSel->Fill(vertex->GetZ());
     }
    else	
    {
		fIsSelected=kFALSE;
    }	
  }
  Int_t Nch=0;
  for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
    AliVTrack* track =dynamic_cast<AliVTrack*>(fAOD->GetTrack(iTracks));
   /* if(fAODEvent==AliSpectraBothTrackCuts::kESDobject)
		track=dynamic_cast<AliVTrack*>(esdevent->GetTrack(iTracks));
     else if (fAODEvent==AliSpectraBothTrackCuts::kAODobject)
		track=dynamic_cast<AliVTrack*>(aodevent->GetTrack(iTracks));
     else return false;*/
     
    if (!fTrackCuts->IsSelected(track,kFALSE)) continue;
    fHistoEtaBefSel->Fill(track->Eta());
    if(fIsSelected){
      fHistoEtaAftSel->Fill(track->Eta());
      Nch++;
    }
  }
  //Printf("NCHARGED_EvSel : %d",Nch);
  if(fIsSelected)fHistoNChAftSel->Fill(Nch);
  return fIsSelected;
}

//______________________________________________________
Bool_t AliSpectraBothEventCuts::CheckVtx()
{
  // reject events outside of range
  const AliVVertex * vertex = fAOD->GetPrimaryVertex();
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
   if(tmp.Contains("NoConstraint"))
   {
        fHistoCuts->Fill(kTPCasPV);
        return kFALSE;
   }


 // if (vertex->GetZ() > fVertexCutMin && vertex->GetZ() < fVertexCutMax)
   // {
    //  return kTRUE;
   // }
  fHistoCuts->Fill(kVtxRange);
  //return kFALSE;
   return kTRUE;
}

//______________________________________________________
Bool_t AliSpectraBothEventCuts::CheckCentralityCut()
{
  // Check centrality cut
  if ( fCentralityCutMax<0.0  &&  fCentralityCutMin<0.0 )  return kTRUE;   
  Double_t cent=0;
  if(!fUseCentPatchAOD049)cent=fAOD->GetCentrality()->GetCentralityPercentile(fCentEstimator.Data());
  else cent=ApplyCentralityPatchAOD049();
  
  if ( (cent <= fCentralityCutMax)  &&  (cent >= fCentralityCutMin) )  return kTRUE;   
  fHistoCuts->Fill(kVtxCentral);

  return kFALSE;
}

//______________________________________________________
Bool_t AliSpectraBothEventCuts::CheckMultiplicityCut()
{
  // Check multiplicity cut
if(fMultiplicityCutMin<0.0 && fMultiplicityCutMax<0.0)
	return kTRUE;
  Int_t Ncharged=0;
  for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++){
    AliVTrack* track = dynamic_cast<AliVTrack*>(fAOD->GetTrack(iTracks));

    if (!fTrackCuts->IsSelected(track,kFALSE)) continue;
	
    Ncharged++;
  }
  //Printf("NCHARGED_cut : %d",Ncharged);
  if(Ncharged>fMultiplicityCutMin && Ncharged<fMultiplicityCutMax)return kTRUE;
  
  return kFALSE;
}

//______________________________________________________
Bool_t AliSpectraBothEventCuts::CheckQVectorCut()
{
	 if(fQVectorCutMin<0.0 && fQVectorCutMax<0.0)return kTRUE;
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
  cout << " > Number of accepted events: " << fHistoCuts->GetBinContent(kAcceptedEvents + 1) << endl;
  cout << " > Number of processed events: " << fHistoCuts->GetBinContent(kProcessedEvents + 1) << endl;
  cout << " > Number of PhysSel events: " << fHistoCuts->GetBinContent(kPhysSelEvents + 1) << endl;
  cout << " > Vertex out of range: " << fHistoCuts->GetBinContent(kVtxRange + 1) << endl;
  cout << " > Events cut by centrality: " << fHistoCuts->GetBinContent(kVtxCentral + 1) << endl;
  cout << " > Events without vertex: " << fHistoCuts->GetBinContent(kVtxNoEvent + 1) << endl;
  cout << " > QVector cut: " << fHistoCuts->GetBinContent(kQVector + 1) << endl;
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

  Int_t count = 0;

  while ((obj = iter->Next())) {
    AliSpectraBothEventCuts* entry = dynamic_cast<AliSpectraBothEventCuts*> (obj);
    if (entry == 0) 
      continue;

    TH1I * histo = entry->GetHistoCuts();      
    collections.Add(histo);
    TH1F * histo_histoVtxBefSel = entry->GetHistoVtxBefSel();      
    collections_histoVtxBefSel.Add(histo_histoVtxBefSel);
    TH1F * histo_histoVtxAftSel = entry->GetHistoVtxAftSel();      
    collections_histoVtxAftSel.Add(histo_histoVtxAftSel);
    TH1F * histo_histoEtaBefSel = entry->GetHistoEtaBefSel();      
    collections_histoEtaBefSel.Add(histo_histoEtaBefSel);
    TH1F * histo_histoEtaAftSel = entry->GetHistoEtaAftSel();      
    collections_histoEtaAftSel.Add(histo_histoEtaAftSel);
    TH1F * histo_histoNChAftSel = entry->GetHistoNChAftSel();      
    collections_histoNChAftSel.Add(histo_histoNChAftSel);
    // TH1F * histo_histoQVectorPos = entry->GetHistoQVectorPos();      
    // collections_histoQVectorPos.Add(histo_histoQVectorPos);
    // TH1F * histo_histoQVectorNeg = entry->GetHistoQVectorNeg();      
    // collections_histoQVectorNeg.Add(histo_histoQVectorNeg);
    TH1F * histo_histoQVector = entry->GetHistoQVector();      
    collections_histoQVector.Add(histo_histoQVector);
    TH1F * histo_histoEP = entry->GetHistoEP();      
    collections_histoEP.Add(histo_histoEP);
    count++;
  }
  
  fHistoCuts->Merge(&collections);
  fHistoVtxBefSel->Merge(&collections_histoVtxBefSel);
  fHistoVtxAftSel->Merge(&collections_histoVtxAftSel);
  fHistoEtaBefSel->Merge(&collections_histoEtaBefSel);
  fHistoEtaAftSel->Merge(&collections_histoEtaAftSel);
  fHistoNChAftSel->Merge(&collections_histoNChAftSel);
  // fHistoQVectorPos->Merge(&collections_histoQVectorPos);
  // fHistoQVectorNeg->Merge(&collections_histoQVectorNeg);
  fHistoQVector->Merge(&collections_histoQVector);
  fHistoEP->Merge(&collections_histoEP);
  

  delete iter;

  return count+1;
}

