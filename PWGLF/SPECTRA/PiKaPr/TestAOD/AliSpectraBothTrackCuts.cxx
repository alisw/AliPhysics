
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
//         AliSpectraBothTrackCuts class
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"
#include "AliExternalTrackParam.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#include "AliSpectraBothTrackCuts.h"
//#include "AliSpectraBothHistoManager.h"
#include <iostream>

using namespace std;

const char * AliSpectraBothTrackCuts::kBinLabel[] ={"TrkBit",
						   "TrkCuts",
						   "TrkEta",
						   "TrkDCA",
						   "TrkP",
						   "TrkPt",
						   "TrkPtTOF",
						   "TOFMatching",
						   "kTOFout",
						   "kTIME",
						   "kTOFpid",
						   "Accepted"};


ClassImp(AliSpectraBothTrackCuts)


AliSpectraBothTrackCuts::AliSpectraBothTrackCuts(const char *name) : TNamed(name, "AOD Track Cuts"), fIsSelected(0), fTrackBits(0), fMinTPCcls(0), fEtaCutMin(0), fEtaCutMax(0), fDCACut(0), fPCut(0), fPtCut(0), fYCutMax(0),fYCutMin(0),
  fPtCutTOFMatching(0),fAODtrack(kotherobject), fHashitinSPD1(0),fusedadditionalcuts(kTRUE),
fPtCutTOFMatchingPion(-1.0),fPtCutTOFMatchingKaon(-1.0),fPtCutTOFMatchingProton(-1.0),fUseTypeDependedTOFCut(kFALSE),fMakeQAhisto(kFALSE),
fHistoCuts(0), fHistoNSelectedPos(0), fHistoNSelectedNeg(0), fHistoNMatchedPos(0), fHistoNMatchedNeg(0), fHistoEtaPhiHighPt(0), fHistoNclustersITS(0),
fHistoDCAzQA(0),fHistoNclustersQA(0),fHistochi2perNDFQA(0),
fTrack(0),fCuts(0)
  
{
/*
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);	
  // Constructor
  fHistoCuts = new TH1I("fTrkCuts", "Track Cuts", kNTrkCuts, -0.5, kNTrkCuts - 0.5);
  for(Int_t ibin=1;ibin<=kNTrkCuts;ibin++)fHistoCuts->GetXaxis()->SetBinLabel(ibin,kBinLabel[ibin-1]);
  //standard histo
  const Double_t templBins[] = {0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0};
  Int_t nbinsTempl=52;
  
  fHistoNSelectedPos=new TH1F("fHistoNSelectedPos","fHistoNSelectedPos",nbinsTempl,templBins);
  fHistoNSelectedPos->GetXaxis()->SetTitle("P_{T} (GeV / c)");
  fHistoNSelectedNeg=new TH1F("fHistoNSelectedNeg","fHistoNSelectedNeg",nbinsTempl,templBins);
  fHistoNSelectedNeg->GetXaxis()->SetTitle("P_{T} (GeV / c)");
  fHistoNMatchedPos=new TH1F("fHistoNMatchedPos","fHistoNMatchedPos",nbinsTempl,templBins);
  fHistoNMatchedPos->GetXaxis()->SetTitle("P_{T} (GeV / c)");
  fHistoNMatchedNeg=new TH1F("fHistoNMatchedNeg","fHistoNMatchedNeg",nbinsTempl,templBins);
  fHistoNMatchedNeg->GetXaxis()->SetTitle("P_{T} (GeV / c)");
  fHistoEtaPhiHighPt=new TH2F("fHistoEtaPhiHighPt","fHistoEtaPhiHighPt",200,-1,1,400,0,7);
  fHistoEtaPhiHighPt->SetXTitle("eta");
  fHistoEtaPhiHighPt->SetYTitle("phi");
  fHistoNclustersITS=new TH1F("fHistoNclustersITS","fHistoNclustersITS;N;ITSLayer",6,-0.5,5.5);
  */
  fEtaCutMin = -100000.0; // default value of eta cut ~ no cut
  fEtaCutMax = 100000.0; // default value of eta cut ~ no cut
  fDCACut = 100000.0; // default value of dca cut ~ no cut
  fPCut = 100000.0; // default value of p cut ~ no cut
  fPtCut = 100000.0; // default value of pt cut ~ no cut 
  fPtCutTOFMatching=0.6; //default value fot matching with TOF
  fYCutMax       = 100000.0; // default value of y cut ~ no cut 
  fYCutMin       = -100000.0; // default value of y cut ~ no cut 
  fMinTPCcls=70; // ncls in TPC
 //  TH1::AddDirectory(oldStatus);
	
}
//__________________________________________________________________
void AliSpectraBothTrackCuts::InitHisto()
{
	Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);	
 
  fHistoCuts = new TH1I("fTrkCuts", "Track Cuts", kNTrkCuts, -0.5, kNTrkCuts - 0.5);
  for(Int_t ibin=1;ibin<=kNTrkCuts;ibin++)fHistoCuts->GetXaxis()->SetBinLabel(ibin,kBinLabel[ibin-1]);
  //standard histo
  const Double_t templBins[] = {0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0};
  Int_t nbinsTempl=52;
  const  Double_t parBins[]={-0.5,0.5,1.5,2.5};
  Int_t nbinsnpar=3;
  const  Double_t dcazBins[] ={-4.0,-3.8,-3.6,-3.4,-3.2,-3.0,-2.8,-2.6,-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0}; 		  	
  Int_t nbinsdcaz=40;
  const  Double_t nclsBins[]={45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160}; 	 
  Int_t nbinsncls=23;
  const  Double_t chindfBins[]={0.0,0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6,4.0,4.4,4.8,5.2,5.6,6.0,6.4,6.8,7.2,7.6,8.0}; 	
  Int_t nbinchindf=20;

  if(!fHistoNSelectedPos)
  	fHistoNSelectedPos=new TH1F("fHistoNSelectedPos","fHistoNSelectedPos",nbinsTempl,templBins);
  fHistoNSelectedPos->GetXaxis()->SetTitle("P_{T} (GeV / c)");
  if(!fHistoNSelectedNeg)
  	fHistoNSelectedNeg=new TH1F("fHistoNSelectedNeg","fHistoNSelectedNeg",nbinsTempl,templBins);
  fHistoNSelectedNeg->GetXaxis()->SetTitle("P_{T} (GeV / c)");
  if(!fHistoNMatchedPos)
  	fHistoNMatchedPos=new TH1F("fHistoNMatchedPos","fHistoNMatchedPos",nbinsTempl,templBins);
  fHistoNMatchedPos->GetXaxis()->SetTitle("P_{T} (GeV / c)");
  if(!fHistoNMatchedNeg)
  	fHistoNMatchedNeg=new TH1F("fHistoNMatchedNeg","fHistoNMatchedNeg",nbinsTempl,templBins);
  fHistoNMatchedNeg->GetXaxis()->SetTitle("P_{T} (GeV / c)");
  if(!fHistoEtaPhiHighPt)
  	fHistoEtaPhiHighPt=new TH2F("fHistoEtaPhiHighPt","fHistoEtaPhiHighPt",200,-1,1,400,0,7);
  fHistoEtaPhiHighPt->SetXTitle("eta");
  fHistoEtaPhiHighPt->SetYTitle("phi");
  if(!fHistoNclustersITS)
  	fHistoNclustersITS=new TH1F("fHistoNclustersITS","fHistoNclustersITS;N;ITSLayer",6,-0.5,5.5);
  if(fMakeQAhisto)
  {		
   	if(!fHistoDCAzQA)
		fHistoDCAzQA=new TH3F("fHistoDCAzQA","QA of DCA z;par type; P_{T} (GeV/C);dcaz",nbinsnpar,parBins,nbinsTempl,templBins,nbinsdcaz,dcazBins);
   	if(!fHistoNclustersQA)
		fHistoNclustersQA=new TH3F("fHistoNclustersQA","QA of Ncls ;par type ; P_{T} (GeV/C);ncls",nbinsnpar,parBins,nbinsTempl,templBins,nbinsncls,nclsBins);
	   if(!fHistochi2perNDFQA) 
		fHistochi2perNDFQA=new TH3F("fHistochi2perNDFQA","QA of chi2/ndf ;par type; P_{T} (GeV/C);chi2ndf",nbinsnpar,parBins,nbinsTempl,templBins,nbinchindf,chindfBins);
 }
 TH1::AddDirectory(oldStatus);

if(fUseTypeDependedTOFCut) 
	fPtCutTOFMatching=TMath::Max(fPtCutTOFMatchingPion,TMath::Max(fPtCutTOFMatchingKaon,fPtCutTOFMatchingProton));
 


}
//_______________________________________________________
AliSpectraBothTrackCuts::~AliSpectraBothTrackCuts()
{
	if(fHistoCuts)
		delete fHistoCuts;
	if(fHistoNSelectedPos)
		delete fHistoNSelectedPos;
	if(fHistoNSelectedNeg)
		delete fHistoNSelectedNeg;
	if(fHistoNMatchedPos)
		delete fHistoNMatchedPos;
	if(fHistoNMatchedNeg)
		delete fHistoNMatchedNeg;
	if(fHistoEtaPhiHighPt)
		delete fHistoEtaPhiHighPt;
	if(fHistoNclustersITS)
		delete fHistoNclustersITS;
	if(fHistoDCAzQA)
		delete fHistoDCAzQA;
	if(fHistoNclustersQA)
		delete fHistoNclustersQA;
	if(fHistochi2perNDFQA)
		delete fHistochi2perNDFQA;



}
//_______________________________________________________
Bool_t AliSpectraBothTrackCuts::IsSelected(AliVTrack * track,Bool_t FillHistStat)
{
// Returns true if Track Cuts are selected and applied
  if (!track)
    {
      printf("ERROR: Could not receive track");
      return kFALSE;
    }
    fTrack = track;
   TString nameoftrack(track->ClassName());  
    if(!nameoftrack.CompareTo("AliESDtrack"))
		fAODtrack=kESDobject;
	else if(!nameoftrack.CompareTo("AliAODTrack"))
		fAODtrack=kAODobject;
	else
		fAODtrack=kotherobject;
  if(!CheckTrackType()){
    return kFALSE;
  }
  if(FillHistStat)fHistoCuts->Fill(kTrkBit);
  if(!CheckTrackCuts()){
    return kFALSE;
  }
  if(FillHistStat)fHistoCuts->Fill(kTrkCuts);
  if(!CheckEtaCut()){
    return kFALSE;
  }
  if(FillHistStat)fHistoCuts->Fill(kTrkEta);
  if(!CheckDCACut()){
    return kFALSE;
  }
  if(FillHistStat)fHistoCuts->Fill(kTrkDCA);
  if(!CheckPCut()){
    return kFALSE;
  }
  if(FillHistStat)fHistoCuts->Fill(kTrkP);
  if(!CheckPtCut()){
    return kFALSE;
  }
  if(FillHistStat)fHistoCuts->Fill(kTrkPt);
  if(!CheckTOFMatching(FillHistStat)){
    return kFALSE;
  }
  if(FillHistStat)fHistoCuts->Fill(kAccepted);
  //Printf("-------- %d,%d",kTOFMatching,kAccepted);
  
  return kTRUE;
}
//_________________________________________________________

Bool_t AliSpectraBothTrackCuts::CheckTrackType()
{
  // Check track Type
  if(fAODtrack==kESDobject)
  {
	AliESDtrack* esdtrack=dynamic_cast<AliESDtrack*>(fTrack);
	if(!esdtrack)
		return kFALSE;	
	if(fCuts->AcceptTrack(esdtrack)) return kTRUE;
		return kFALSE;
 }
  else if(fAODtrack==kAODobject)
  {
	AliAODTrack* aodtrack=dynamic_cast<AliAODTrack*>(fTrack);
	if(!aodtrack)
		return kFALSE;
	if (aodtrack->TestFilterBit(fTrackBits)) return kTRUE;
		return kFALSE;
  }

  else
	return kFALSE;
  
}
//_________________________________________________________

Bool_t AliSpectraBothTrackCuts::CheckTrackCuts()
{
  // Check additional track Cuts
  Bool_t PassTrackCuts=kTRUE;
  if(!fusedadditionalcuts)
	return PassTrackCuts;
  AliAODTrack* aodtrack=0;
  AliESDtrack* esdtrack=0;
  if(fAODtrack==kESDobject)
  {
	esdtrack=dynamic_cast<AliESDtrack*>(fTrack);
	if(!esdtrack)
		return kFALSE;
	if (!esdtrack->HasPointOnITSLayer(0) && !esdtrack->HasPointOnITSLayer(1))PassTrackCuts=kFALSE; //FIXME 1 SPD for the moment
	if (fHashitinSPD1&&!esdtrack->HasPointOnITSLayer(0)) PassTrackCuts=kFALSE; 		
	if (esdtrack->GetTPCNcls()<fMinTPCcls)PassTrackCuts=kFALSE;
	if(!esdtrack->IsOn(AliESDtrack::kTPCrefit))PassTrackCuts=kFALSE;
	if(!esdtrack->IsOn(AliESDtrack::kITSrefit))PassTrackCuts=kFALSE;
	if(PassTrackCuts)
	{
		for(int i=0;i<6;i++)
			if(esdtrack->HasPointOnITSLayer(i))
				fHistoNclustersITS->Fill(i);
	}	
  }
  else if (fAODtrack==kAODobject)
  	{
	aodtrack=dynamic_cast<AliAODTrack*>(fTrack);
	if(!aodtrack)
		return kFALSE;
	if (!aodtrack->HasPointOnITSLayer(0) && !aodtrack->HasPointOnITSLayer(1))PassTrackCuts=kFALSE; //FIXME 1 SPD for the moment
	if (fHashitinSPD1&&!aodtrack->HasPointOnITSLayer(0)) PassTrackCuts=kFALSE; 	
	if (aodtrack->GetTPCNcls()<fMinTPCcls)PassTrackCuts=kFALSE;
	if(!aodtrack->IsOn(AliAODTrack::kTPCrefit))PassTrackCuts=kFALSE;
	if(!aodtrack->IsOn(AliAODTrack::kITSrefit))PassTrackCuts=kFALSE;
	if(PassTrackCuts)
	{
		for(int i=0;i<6;i++)
			if(aodtrack->HasPointOnITSLayer(i))
				fHistoNclustersITS->Fill(i);
	}	
  }
  else
	return kFALSE;
    
  
  
  return PassTrackCuts;
}
//________________________________________________________
Bool_t AliSpectraBothTrackCuts::CheckEtaCut()
{
   // Check eta cut
   if (fTrack->Eta() < fEtaCutMax && fTrack->Eta() > fEtaCutMin) return kTRUE;
    return kFALSE;
}

Bool_t AliSpectraBothTrackCuts::CheckYCut(BothParticleSpecies_t species) 
{
 if(fYCutMin>fYCutMax) 	
	return kTRUE;
  // check if the rapidity is within the set range
  Double_t y;
  
  Double_t pz=fTrack->Pz();
  Double_t p=fTrack->P();
  Double_t mass=-1.0;
  /*
  if (species == kSpProton) { y = fTrack->Y(9.38271999999999995e-01); }
  if ( species == kSpKaon ) { y = fTrack->Y(4.93676999999999977e-01); }
  if ( species == kSpPion)  { y = fTrack->Y(1.39570000000000000e-01); }
  
  */
    if (species == kSpProton) { mass=9.38271999999999995e-01; }
  if ( species == kSpKaon ) { mass=4.93676999999999977e-01; }
  if ( species == kSpPion)  { mass=1.39570000000000000e-01; }
  if(mass<0.0)
	y =-999.0 ;
  else
	y=0.5*TMath::Log((TMath::Sqrt(mass*mass+p*p)+pz)/(TMath::Sqrt(mass*mass+p*p)-pz));
  if (y > fYCutMax || y<fYCutMin||y < -998.) return kFALSE;
	return kTRUE;
}
//_______________________________________________________
Bool_t AliSpectraBothTrackCuts::CheckDCACut()
{
   // Check DCA cut
 // if (TMath::Abs(fTrack->DCA()) < fDCACut) return kTRUE; //FIXME for newest AOD fTrack->DCA() always gives -999
   
    AliAODTrack* aodtrack=0;
  AliESDtrack* esdtrack=0;
  if(fAODtrack==kESDobject)
  {
	esdtrack=dynamic_cast<AliESDtrack*>(fTrack);
	if(!esdtrack)
		return kFALSE;
	Float_t dcaxy=0.0; 
	Float_t dcaz=0.0;
	esdtrack->GetImpactParameters(dcaxy,dcaz);	
	if (TMath::Abs(dcaxy) < fDCACut) 
		return kTRUE;
	else 
		return kFALSE;
   }
  else if (fAODtrack==kAODobject)
  {
	aodtrack=dynamic_cast<AliAODTrack*>(fTrack);
	if(!aodtrack)
		return kFALSE;
	if (TMath::Abs(aodtrack->DCA()) < fDCACut) return kTRUE;
	else 
		return kFALSE;
		
   }
	else
	return kFALSE;
}
//________________________________________________________
Bool_t AliSpectraBothTrackCuts::CheckPCut()
{
   // Check P cut
   if (fTrack->P() < fPCut) return kTRUE;
   return kFALSE;
}
//_______________________________________________________
Bool_t AliSpectraBothTrackCuts::CheckPtCut()
{
    // check Pt cut
//    if ((fTrack->Pt() < fPtCut) && (fTrack->Pt() > 0.3 )) return kTRUE;
   if (fTrack->Pt() < fPtCut) return kTRUE;
    return kFALSE;
}

//_______________________________________________________
Bool_t AliSpectraBothTrackCuts::CheckTOFMatching(Bool_t FillHistStat)
{
  // check Pt cut
  //    if ((fTrack->Pt() < fPtCut) && (fTrack->Pt() > 0.3 )) return kTRUE;
	if (fTrack->Pt() < fPtCutTOFMatching) 
		return kTRUE;
	else
 	{
		if(FillHistStat)
			fHistoCuts->Fill(kTrkPtTOF);
    		if(fTrack->Charge()>0)
			fHistoNSelectedPos->Fill(fTrack->Pt());
    		else 
			fHistoNSelectedNeg->Fill(fTrack->Pt());
    		UInt_t status=fTrack->GetStatus();
    		if((status&AliAODTrack::kTOFout)&&FillHistStat)
			fHistoCuts->Fill(kTrTOFout);
    		if((status&AliAODTrack::kTIME)&&FillHistStat)
			fHistoCuts->Fill(kTrTIME);
    		if((status&AliAODTrack::kTOFpid)&&FillHistStat)
			fHistoCuts->Fill(kTrTOFpid);
    
   		 if((status&AliAODTrack::kTOFout)==0 || (status&AliAODTrack::kTIME)==0)
		{//kTOFout and kTIME
      			return kFALSE; 
    		} 
    		if(FillHistStat)
			fHistoCuts->Fill(kTOFMatching);
    		if(fTrack->Charge()>0)
			fHistoNMatchedPos->Fill(fTrack->Pt());
   		 else 
			fHistoNMatchedNeg->Fill(fTrack->Pt());
    		if(fTrack->Pt()>1.5)
		{
      			//fHistoEtaPhiHighPt->Fill(fTrack->GetOuterParam()->Eta(),fTrack->GetOuterParam()->Phi());
      			//Printf("AliExternalTrackParam * extpar=(AliExternalTrackParam*)fTrack->GetOuterParam();");
      			//AliExternalTrackParam * extpar=(AliExternalTrackParam*)fTrack->GetOuterParam();
      			fHistoEtaPhiHighPt->Fill(fTrack->Eta(),fTrack->Phi());
     			 //Printf("fHistoEtaPhiHighPt->Fill(extpar->Eta(),extpar->Phi());");
      			//fHistoEtaPhiHighPt->Fill(extpar->Eta(),extpar->Phi());
      			//delete extpar;
    		}
    		return kTRUE;
  	}
}
//_______________________________________________________
void AliSpectraBothTrackCuts::PrintCuts() const
{
  // Print cuts
    cout << "Track Cuts" << endl;
    cout << " > TrackBit\t" << fTrackBits << endl;
    cout << " > Eta cut\t" << fEtaCutMin <<","<< fEtaCutMax << endl;
    cout << " > DCA cut\t" << fDCACut << endl;
    cout << " > P cut\t" << fPCut << endl;
    cout << " > Pt cut \t" << fPtCut << endl;
    cout << " > TPC cls \t" << fMinTPCcls << endl;
}
//_______________________________________________________
void AliSpectraBothTrackCuts::SetTrackType(UInt_t bit)
{
   // Set the type of track to be used. The argument should be the bit number. The mask is produced automatically.
   fTrackBits = (0x1 << (bit - 1));
}
//_______________________________________________________

Long64_t AliSpectraBothTrackCuts::Merge(TCollection* list)
{
  // Merge a list of AliSpectraBothTrackCuts objects with this.
  // Returns the number of merged objects (including this).

  //  AliInfo("Merging");


  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of all histograms
  TList collections;//FIXME we should only 1 collection
  TList collections_histoNSelectedPos;
  TList collections_histoNSelectedNeg;
  TList collections_histoNMatchedPos;
  TList collections_histoNMatchedNeg;
  TList collections_histoEtaPhiHighPt;
  TList collections_histoDCAzQA;
  TList collections_histoNclustersQA;	
  TList collections_histochi2perNDFQA;
	
  Int_t count = 0;

  while ((obj = iter->Next())) {
    AliSpectraBothTrackCuts* entry = dynamic_cast<AliSpectraBothTrackCuts*> (obj);
    if (entry == 0) 
      continue;
    
    TH1I * histo = entry->GetHistoCuts();   
    if(histo)	   
    	collections.Add(histo);
    TH1F * histoNSelectedPos = entry->GetHistoNSelectedPos();
    if(histoNSelectedPos)	      
    	collections_histoNSelectedPos.Add(histoNSelectedPos);
    TH1F * histoNSelectedNeg = entry->GetHistoNSelectedNeg();      
    if(histoNSelectedNeg)
    	collections_histoNSelectedNeg.Add(histoNSelectedNeg);
    TH1F * histoNMatchedPos = entry->GetHistoNMatchedPos();  
    if(histoNMatchedPos)	    
    	collections_histoNMatchedPos.Add(histoNMatchedPos);
    TH1F * histoNMatchedNeg = entry->GetHistoNMatchedNeg();      
    if(histoNMatchedNeg)
    	collections_histoNMatchedNeg.Add(histoNMatchedNeg);
    TH2F * histoEtaPhiHighPt = entry->GetHistoEtaPhiHighPt();      
    if(histoEtaPhiHighPt)
    	collections_histoEtaPhiHighPt.Add(histoEtaPhiHighPt);
    TH3F* histoDCAzQA=entry->GetHistoDCAzQA();
    if(histoDCAzQA)	 	
    	collections_histoDCAzQA.Add(histoDCAzQA);
    TH3F* histoNclustersQA=entry->GetHistoNclustersQA(); 	
    if(histoNclustersQA)
   	collections_histoNclustersQA.Add(histoNclustersQA);
    TH3F* histochi2perNDFQA=entry->GetHistochi2perNDFQA(); 	
    if(histochi2perNDFQA)
    	collections_histochi2perNDFQA.Add(histochi2perNDFQA);
	


    count++;
  }
  if(fHistoCuts)	 
  	fHistoCuts->Merge(&collections);
  if(fHistoNSelectedPos)
  	fHistoNSelectedPos->Merge(&collections_histoNSelectedPos);
  if(fHistoNSelectedNeg)
  	fHistoNSelectedNeg->Merge(&collections_histoNSelectedNeg);
  if(fHistoNMatchedPos)
  	fHistoNMatchedPos->Merge(&collections_histoNMatchedPos);
  if(fHistoNMatchedNeg)
  	fHistoNMatchedNeg->Merge(&collections_histoNMatchedNeg);
  if(fHistoEtaPhiHighPt)
  	fHistoEtaPhiHighPt->Merge(&collections_histoEtaPhiHighPt);
  if(fHistoDCAzQA)
 	 fHistoDCAzQA->Merge(&collections_histoDCAzQA);
  if(fHistoNclustersQA)
  	fHistoNclustersQA->Merge(&collections_histoNclustersQA);
  if(fHistochi2perNDFQA)
  fHistochi2perNDFQA->Merge(&collections_histochi2perNDFQA);

  delete iter;

  return count+1;
}
//________________________________________________________________________________________________________________________
  void AliSpectraBothTrackCuts::SetPtTOFMatchingPartDepended(Float_t pion,Float_t kaon,Float_t proton)
{
	// configure the task in case of the use particle depended TOF matching cut 
	// the value fPtCutTOFMatching is set as max pt of all
	fPtCutTOFMatchingPion=pion;
	fPtCutTOFMatchingKaon=kaon;
	fPtCutTOFMatchingProton=proton;
	fUseTypeDependedTOFCut=kTRUE;
	fPtCutTOFMatching=TMath::Max(fPtCutTOFMatchingPion,TMath::Max(fPtCutTOFMatchingKaon,fPtCutTOFMatchingProton));
 
}
//___________________________________________________________________________________________________________________________
Bool_t AliSpectraBothTrackCuts::CheckTOFMatchingParticleType(Int_t type) 
{
	Float_t ptcut=fPtCutTOFMatching;
	switch (type) 
	{	
	    case  kSpPion:
	    ptcut=fPtCutTOFMatchingPion;
	    break;
	    case kSpKaon:
	    ptcut=fPtCutTOFMatchingKaon;
	    break;
	    case  kSpProton:
	    ptcut=fPtCutTOFMatchingProton;
	    break;
	    default: 
	    ptcut=0.0;
	    break;		
		

	}	
	if (fTrack->Pt() < ptcut) 
		return kTRUE;
	else
 	{
		UInt_t status=fTrack->GetStatus();
   		 if((status&AliAODTrack::kTOFout)==0 || (status&AliAODTrack::kTIME)==0)
		{//kTOFout and kTIME
      			return kFALSE; 
    		} 
    		return kTRUE;
  	}

}
