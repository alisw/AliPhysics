#include "AliLog.h"

#include "TROOT.h"
#include "TString.h"
#include "TH1.h"
#include "THnSparse.h"
#include "TTreeStream.h"
#include "TGeoGlobalMagField.h"

#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "TParticlePDG.h"
#include "TParticle.h"

#include "AliMESppColTask.h"
#include "AliMESeventInfo.h"
#include "AliMEStrackInfo.h"
#include "AliMCEvent.h"

#include "AliEventPoolManager.h"
#include "AliBasicParticle.h"

ClassImp(AliMESppColTask)


//________________________________________________________________________
AliMESppColTask::AliMESppColTask()
  : AliMESbaseTask()
  ,fPoolMgr1(0x0)
  ,fPoolMgr2(0x0)
  ,fPoolMgr3(0x0)
  ,fPoolMgrMC1(0x0)
  ,fPoolMgrMC2(0x0)
  ,fPoolMgrMC3(0x0)
{
  //
  // Constructor
  //
}

//________________________________________________________________________
AliMESppColTask::AliMESppColTask(const char *name)
  : AliMESbaseTask(name)
  ,fPoolMgr1(0x0)
  ,fPoolMgr2(0x0)
  ,fPoolMgr3(0x0)
  ,fPoolMgrMC1(0x0)
  ,fPoolMgrMC2(0x0)
  ,fPoolMgrMC3(0x0)
{
  //
  // Constructor
  //
  DefineOutput(AliMESbaseTask::kQA, TList::Class());
}

//________________________________________________________________________
AliMESppColTask::~AliMESppColTask()
{
  //deconstructor
  if(fPoolMgr1) {delete fPoolMgr1; fPoolMgr1=0;}
  if(fPoolMgr2) {delete fPoolMgr2; fPoolMgr2=0;}
  if(fPoolMgr3) {delete fPoolMgr3; fPoolMgr3=0;}
  if(fPoolMgrMC1) {delete fPoolMgrMC1; fPoolMgrMC1=0;}
  if(fPoolMgrMC2) {delete fPoolMgrMC2; fPoolMgrMC2=0;}
  if(fPoolMgrMC3) {delete fPoolMgrMC3; fPoolMgrMC3=0;}
}

//________________________________________________________________________
void AliMESppColTask::UserCreateOutputObjects()
{
  //define user data containers
  
  AliMESbaseTask::UserCreateOutputObjects();

  
	DefineMixedEventPool(0);   
	if(HasMCdata()){
		DefineMixedEventPool(1);
	}
}

//________________________________________________________________________
void AliMESppColTask::UserExec(Option_t *opt)
{
  // Run user analysis. The following objects are allocated after calling AliMESbaseTask::UserExec(opt)
  // fEvInfo  -  reconstructed event information (class AliMESeventInfo)
  // fTracks  -  reconstructed array of tracks (class TObjArray of AliMEStrackInfo)
  // fMCevInfo-  MC event information (class AliMESeventInfo)
  // fMCtracks-  MC array of tracks (class TObjArray of AliMEStrackInfo)
  AliMESbaseTask::UserExec(opt);
  
	//event selectors
  if(!(fEvInfo   = dynamic_cast<AliMESeventInfo*>(GetInputData(AliMESbaseTask::kEventInfo)))){ 
	AliInfo("Ev Info Not defined.");
	return;
	}
	//trigger selectors
	if(!fEvInfo->HasTriggerMB()) return ; //Minimum Bias Trigger
// 	if(!fEvInfo->HasTriggerHM()) return ; //High Multiplicity Trigger
	
	Double_t vec_hNoEvts[7]; // vector used to fill hNoEvts
	THnSparseD *hNoEvts = (THnSparseD*)fHistosQA->At(0);
    
	Double_t mult_comb08 = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb);// combined multiplicity with |eta| < 0.8
	// event shape for data (from ESD)
// 	Double_t directivity_plus = fEvInfo->GetEventShape()->GetDirectivity(1);
	Double_t directivity_minus = fEvInfo->GetEventShape()->GetDirectivity(0);
	Double_t sfer = fEvInfo->GetEventShape()->GetSphericity();
	
	vec_hNoEvts[0] = 0.;
	hNoEvts->Fill(vec_hNoEvts);
	

// 	// select events with both dirs in the same interval
// 	const Int_t lenght = 4;
// 	Double_t intervals[lenght] = {0., 0.3, 0.6, 0.9};
// 	// NOTE: the intervals are considered half-closed: (a,b]
// 	if( (directivity_plus < intervals[0]) || (directivity_plus > intervals[lenght-1]) ) return;
// 	vec_hNoEvts[0] = 1.;
// 	hNoEvts->Fill(vec_hNoEvts);
// 	if( (directivity_minus < intervals[0]) || (directivity_minus > intervals[lenght-1]) ) return;
// 	vec_hNoEvts[0] = 2.;
// 	hNoEvts->Fill(vec_hNoEvts);
//   
// 	Int_t first = -1;
// 	for(Int_t i=1; i<lenght; i++){
// 		if(directivity_plus <= intervals[i]){
// 			first = i;
// 			break;
// 		}
// 	}
// 	if( (directivity_minus <= intervals[first-1]) || (directivity_minus > intervals[first]) ) return;
// 	vec_hNoEvts[0] = 3.;
// 	hNoEvts->Fill(vec_hNoEvts);
	
	Double_t directivity = directivity_minus;
	
	vec_hNoEvts[1] = mult_comb08; // combined multiplicity with |eta| < 0.8
	vec_hNoEvts[2] = directivity;
	vec_hNoEvts[3] = sfer; 
	
	// event multiplicity and shape for MC (from MC event)
	Double_t MC_directivity_plus = -2;
	Double_t MC_directivity_minus = -2;
	Double_t MC_directivity = -2;
	Double_t MC_mult_glob08 = -2;
	Double_t MC_sfer = -2;
	if( HasMCdata() ){ // run only on MC
// 		MC_directivity_plus = fMCevInfo->GetEventShape()->GetDirectivity(1);
		MC_directivity_minus = fMCevInfo->GetEventShape()->GetDirectivity(0);
// 		MC_directivity = (MC_directivity_plus + MC_directivity_minus) / 2.0;
		MC_directivity= MC_directivity_minus;
		MC_mult_glob08 = fMCevInfo->GetMultiplicity(AliMESeventInfo::kGlob08);
		MC_sfer = fMCevInfo->GetEventShape()->GetSphericity();
		vec_hNoEvts[4] = fMCevInfo->GetMultiplicity(AliMESeventInfo::kGlob08);
		vec_hNoEvts[5] = MC_directivity;
		vec_hNoEvts[6] = MC_sfer;
	}
	
	vec_hNoEvts[0] = 4.;
	hNoEvts->Fill(vec_hNoEvts);
	
	//-------------------ESD Loop--------------------
	AliMEStrackInfo *t(NULL), *tMC(NULL);
  for(Int_t it(0); it<fTracks->GetEntries(); it++){
    if(!(t = (AliMEStrackInfo*)fTracks->At(it))) continue;
	 if( !(t->HasOrigin(AliMEStrackInfo::kPrimary)) ) continue;
	 if( TMath::Abs(t->Eta())> 0.8 ) continue;
	 if( t->Pt()< 0.2 ) continue;
	 if( HasMCdata() ){
		if( !(tMC= (AliMEStrackInfo*)fMCtracks->At(t->GetLabel())) ) continue;
	}
	Double_t vec_hTrk[4];
	THnSparseD *hTrk = (THnSparseD*)fHistosQA->At(145);
	vec_hTrk[0]=mult_comb08;
	vec_hTrk[1]=directivity;
	vec_hTrk[2]=sfer;
	vec_hTrk[3]=t->Pt();
	
	hTrk->Fill(vec_hTrk);
  }
	if( HasMCdata() ){
	  for(Int_t it(0); it<fMCtracks->GetEntries(); it++){
    	if(!(tMC = (AliMEStrackInfo*)fMCtracks->At(it))) continue;
		if( !(tMC->HasOrigin(AliMEStrackInfo::kPrimary)) ) continue;
		if( TMath::Abs(tMC->Eta())> 0.8 ) continue;
		if( tMC->Pt()< 0.2 ) continue;
		Double_t vec_hMCTrk[4];
		THnSparseD *hMCTrk = (THnSparseD*)fHistosQA->At(146);
		vec_hMCTrk[0]=MC_mult_glob08;
		vec_hMCTrk[1]=MC_directivity;
		vec_hMCTrk[2]=MC_sfer;
		vec_hMCTrk[3]=tMC->Pt();
		hMCTrk->Fill(vec_hMCTrk);
	  }
	}
	  //---------------end of ESD loop--------------------
	
	
	
	// ESD tracks
	Int_t ESD=1;
    
    do{
      // NOTE: the intervals are considered half-closed: (a,b]
      if(mult_comb08>0 && mult_comb08<=80 && TMath::Abs(fEvInfo->GetVertexZ())<10.0 && sfer>0.0 && sfer<=0.3){
        TObjArray *selectedTracks1=FindLeadingObjects(fTracks, 0);
        if(!selectedTracks1) break;
        selectedTracks1->SetOwner(kTRUE);
        FillCorrelationSE(mult_comb08, selectedTracks1, 3, 0);
        FillCorrelationMixing(mult_comb08, fEvInfo->GetVertexZ(), 80., 0., selectedTracks1, 3, 0);
      }
      if(mult_comb08>0 && mult_comb08<=80 && TMath::Abs(fEvInfo->GetVertexZ())<10.0 && sfer>0.3 && sfer<=0.6){
        TObjArray *selectedTracks2=FindLeadingObjects(fTracks, 0);
        if(!selectedTracks2) break;
        selectedTracks2->SetOwner(kTRUE);
        FillCorrelationSE(mult_comb08, selectedTracks2, 6, 0);
        FillCorrelationMixing(mult_comb08, fEvInfo->GetVertexZ(), 80., 0., selectedTracks2, 6, 0);
      }
      if(mult_comb08>0 && mult_comb08<=80 && TMath::Abs(fEvInfo->GetVertexZ())<10.0 && sfer>0.6 && sfer<=1.0){
        TObjArray *selectedTracks3=FindLeadingObjects(fTracks, 0);
        if(!selectedTracks3) break;
        selectedTracks3->SetOwner(kTRUE);
        FillCorrelationSE(mult_comb08, selectedTracks3, 9, 0);
        FillCorrelationMixing(mult_comb08, fEvInfo->GetVertexZ(), 80., 0., selectedTracks3, 9, 0);
      }
        ESD=0;
    }while(ESD==1);

  
	if( HasMCdata() ){// run only on MC  
      // NOTE: the intervals are considered half-closed: (a,b]
      if(MC_mult_glob08>0 && MC_mult_glob08<=80 && TMath::Abs(fMCevInfo->GetVertexZ())<10.0 && MC_sfer>0.0 && MC_sfer<=0.3){
		TObjArray *selectedTracksMC1=FindLeadingObjects(fMCtracks, 1);
		if(!selectedTracksMC1) return;
		selectedTracksMC1->SetOwner(kTRUE);
        FillCorrelationSE(MC_mult_glob08, selectedTracksMC1, 3, 1);
        FillCorrelationMixing(MC_mult_glob08, fMCevInfo->GetVertexZ(), 80., 0., selectedTracksMC1, 3, 1);
      }
      if(MC_mult_glob08>0 && MC_mult_glob08<=80 && TMath::Abs(fMCevInfo->GetVertexZ())<10.0 && MC_sfer>0.3 && MC_sfer<=0.6){
		TObjArray *selectedTracksMC2=FindLeadingObjects(fMCtracks, 1);
		if(!selectedTracksMC2) return;
		selectedTracksMC2->SetOwner(kTRUE);
        FillCorrelationSE(MC_mult_glob08, selectedTracksMC2, 6, 1);
        FillCorrelationMixing(MC_mult_glob08, fMCevInfo->GetVertexZ(), 80., 0., selectedTracksMC2, 6, 1);
      }
      if(MC_mult_glob08>0 && MC_mult_glob08<=80 && TMath::Abs(fMCevInfo->GetVertexZ())<10.0 && MC_sfer>0.6 && MC_sfer<=1.0){
		TObjArray *selectedTracksMC3=FindLeadingObjects(fMCtracks, 1);
		if(!selectedTracksMC3) return;
		selectedTracksMC3->SetOwner(kTRUE);
        FillCorrelationSE(MC_mult_glob08, selectedTracksMC3, 9, 1);
        FillCorrelationMixing(MC_mult_glob08, fMCevInfo->GetVertexZ(), 80., 0., selectedTracksMC3, 9, 1);
      }
    }
}
  
//________________________________________________________________________
Bool_t AliMESppColTask::PostProcess()
{
  return kTRUE;
}
    
//________________________________________________________________________

TObjArray* AliMESppColTask::CloneTracks(TObjArray* tracks)
{
  // clones a track list by using AliBasicParticle which uses much less memory (used for event mixing)
  
  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);

  for (Int_t i=0; i<tracks->GetEntries(); i++)
  {
    AliMEStrackInfo* particle = (AliMEStrackInfo*) tracks->At(i);
    AliBasicParticle* copy = 0;


	copy = new AliBasicParticle(particle->Eta(), particle->Phi(), particle->Pt(), 1);
	copy->SetUniqueID(particle->GetLabel());

    tracksClone->Add(copy);
  }
  return tracksClone;
}


//________________________________________________________________________

TObjArray*  AliMESppColTask::FindLeadingObjects(TObjArray *obj, Int_t MC)
{
// Returns an array of charged particles ordered according to their pT.
  	Int_t nTracks = obj->GetEntries();
	if( !nTracks ) return 0;
    
        // Define array of AliMEStrackInfo objects
    TObjArray* tracks = new TObjArray(nTracks);
    TObjArray* tracksMC = new TObjArray(nTracks);
    AliMEStrackInfo* partMC(NULL);
    AliMEStrackInfo* part(NULL);
    
  if(MC==0){
	// Loop over tracks
      for (Int_t ipart=0; ipart<nTracks; ++ipart) {
        if(!(part = (AliMEStrackInfo*)obj->At(ipart))) continue;
          if(HasMCdata() && MC==0){
            if(!(partMC=(AliMEStrackInfo*)fMCtracks->At(part->GetLabel()))) continue;
          }
          if( !(part->HasOrigin(AliMEStrackInfo::kPrimary)) ) continue;
      // Accept tracks in a limited pT & rapidity range
          if(part->Pt()<1.0 || part->Pt()>2.0) continue;
          if( TMath::Abs(part->Eta())> 0.8 ) continue;
          tracks->AddLast(part);
    }
    // Order tracks by pT, first track is LeadingParticle
	QSortTracks( *tracks, 0, tracks->GetEntriesFast() );
	nTracks = tracks->GetEntriesFast();
	if( !nTracks ) return 0;
    
	TObjArray* ClonedTracks = CloneTracks(tracks);
	ClonedTracks->SetOwner(kTRUE);
    return ClonedTracks;
    }
	
	if(HasMCdata() && MC==1){
        for (Int_t ipart=0; ipart<nTracks; ++ipart) {
          if(!(partMC = (AliMEStrackInfo*)obj ->At(ipart))) continue;
          if(!(partMC->HasOrigin(AliMEStrackInfo::kPrimary))) continue;
        // Accept tracks in a limited pT & rapidity range
          if(partMC->Pt()<1.0 || partMC->Pt()>2.0) continue;
          if( TMath::Abs(partMC->Eta())> 0.8 ) continue;
          tracksMC->AddLast(partMC);
        }
      QSortTracks( *tracksMC, 0, tracksMC->GetEntriesFast() );
      nTracks = tracksMC->GetEntriesFast();
      if( !nTracks ) return 0;

      TObjArray* ClonedTracksMC = CloneTracks(tracksMC);
      ClonedTracksMC->SetOwner(kTRUE);
      return ClonedTracksMC;
    }

  }

void  AliMESppColTask::QSortTracks(TObjArray &a, Int_t first, Int_t last)
{
	// Sort array of TObjArray of tracks by pT using a quicksort algorithm.
	static TObject *tmp;
	static int i;           // "static" to save stack space
	int j;
  
	while (last - first > 1) {
    i = first;
    j = last;
    for (;;) {
      while (++i < last && ((AliVParticle*)a[i])->Pt() > ((AliVParticle*)a[first])->Pt() )
        ;
      while (--j > first && ((AliVParticle*)a[j])->Pt() < ((AliVParticle*)a[first])->Pt() )
        ;
      if (i >= j)
        break;
      
      tmp  = a[i];
      a[i] = a[j];
      a[j] = tmp;
    }
    if (j == first) {
      ++first;
      continue;
    }
    tmp = a[first];
    a[first] = a[j];
    a[j] = tmp;
    if (j - first < last - (j + 1)) {
      QSortTracks(a, first, j);
      first = j + 1;   // QSortTracks(j + 1, last);
    } else {
      QSortTracks(a, j + 1, last);
      last = j;        // QSortTracks(first, j);
    }
  }
}

//____________________________________________________________________
Double_t AliMESppColTask::RangePhi(Double_t DPhi)
{
	if (DPhi < -TMath::Pi()/2)  DPhi += 2*TMath::Pi();
	if (DPhi > 3*TMath::Pi()/2) DPhi -= 2*TMath::Pi();      
	return DPhi;    
}

Bool_t AliMESppColTask::DefineMixedEventPool(Int_t MC)
{
  Int_t PoolMaxNEvents = 200;
  Int_t PoolMinNTracks = 20000;
  
  Int_t NMultBins = 12;
  Double_t MultBins[] = { 1., 4., 7., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80. };
  Int_t NzVtxBins = 5;
  Double_t ZvtxBins[] = {  -10., -5., -2.5, 2.5, 5.,  10. }; 
  
    if(HasMCdata() && MC==1){
    fPoolMgrMC1 = new AliEventPoolManager(PoolMaxNEvents, PoolMinNTracks, NMultBins, MultBins, NzVtxBins, ZvtxBins);
    fPoolMgrMC1->SetDebug(0);
    fPoolMgrMC2 = new AliEventPoolManager(PoolMaxNEvents, PoolMinNTracks, NMultBins, MultBins, NzVtxBins, ZvtxBins);
    fPoolMgrMC2->SetDebug(0);
    fPoolMgrMC3 = new AliEventPoolManager(PoolMaxNEvents, PoolMinNTracks, NMultBins, MultBins, NzVtxBins, ZvtxBins);
    fPoolMgrMC3->SetDebug(0);
//     if(!fPoolMgrMC) return kFALSE;
  }
  
  if(MC==0){
  fPoolMgr1 = new AliEventPoolManager(PoolMaxNEvents, PoolMinNTracks, NMultBins, MultBins, NzVtxBins, ZvtxBins);
  fPoolMgr1->SetDebug(0);
  fPoolMgr2 = new AliEventPoolManager(PoolMaxNEvents, PoolMinNTracks, NMultBins, MultBins, NzVtxBins, ZvtxBins);
  fPoolMgr2->SetDebug(0);
  fPoolMgr3 = new AliEventPoolManager(PoolMaxNEvents, PoolMinNTracks, NMultBins, MultBins, NzVtxBins, ZvtxBins);
  fPoolMgr3->SetDebug(0);
//   if(!fPoolMgr) return kFALSE;
  }
  
  
  return kTRUE;

}
    
void AliMESppColTask::FillCorrelationSE(Double_t MultipOrCent, TObjArray*selectedArray, Int_t d, Int_t MC)
{
	const Int_t nMult(12);
	Double_t multBin[nMult+1] = {1., 4., 7., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80.};
	Int_t jm(-1);
	for(Int_t im(0); im<nMult; im++){
		if(MultipOrCent>= multBin[im] && MultipOrCent<multBin[im+1]){
			jm = im;  
			break;
      }
	}
	AliMEStrackInfo* trigger = (AliMEStrackInfo*)selectedArray->At(0);
	if(!trigger) return;
	if(trigger->Pt()<1.0 || trigger->Pt()>2.0) return;
	if( TMath::Abs(trigger->Eta())> 0.8 ) return;
                      
	Double_t ptL  = trigger->Pt();
	Double_t phiL = trigger->Phi();
	Double_t etaL = trigger->Eta();


	for (Int_t j=1; j<selectedArray->GetEntriesFast(); j++){
		AliMEStrackInfo* associate = (AliMEStrackInfo*)selectedArray->At(j);
		if(!associate) continue;
		if(associate->Pt()<1.0 || associate->Pt()>2.0) continue;
		if( TMath::Abs(associate->Eta())> 0.8 ) continue;

		Double_t ptAs = associate->Pt();
		Double_t phiAs = associate->Phi();
		Double_t etaAs = associate->Eta();

		Double_t dPhi(-999.), dEta(-999.);
		if(ptL>ptAs && jm>-1 && ptL>=1.0 && ptL<=2.0){                   
			dPhi = RangePhi(phiL-phiAs);
			dEta=etaL-etaAs;
			if(MC==0){
              if(d==3){
                  Int_t bin=(1+jm);
                  ((TH2*)fHistosQA->At(bin))->Fill(dEta, dPhi); 
              }
              if(d==6){
                Int_t bin=(13+jm);
                ((TH2*)fHistosQA->At(bin))->Fill(dEta, dPhi);
              }
              if(d==9){
                Int_t bin=(25+jm);
                ((TH2*)fHistosQA->At(bin))->Fill(dEta, dPhi);
              }
			}
			if(HasMCdata() && MC==1){
              if(d==3){
				Int_t bin=(73+jm);
				((TH2*)fHistosQA->At(bin))->Fill(dEta, dPhi);
              }
              if(d==6){
                Int_t bin=(85+jm);
				((TH2*)fHistosQA->At(bin))->Fill(dEta, dPhi);
              }
              if(d==9){
                Int_t bin=(98+jm);
				((TH2*)fHistosQA->At(bin))->Fill(dEta, dPhi);
              }
			}
        }
      }
}

//---------------------------------------------------------------------------------------
void AliMESppColTask::FillCorrelationMixing(Double_t MultipOrCentMix, Double_t Zvtx, Double_t poolmax, Double_t poolmin, TObjArray*selectedArray, Int_t d, Int_t MC)
{
	const Int_t nMult(12);
	Int_t multBin[nMult+1] = {1, 4, 7, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80};
	Int_t jm(-1);
	for(Int_t im(0); im<nMult; im++){
		if(MultipOrCentMix>=multBin[im] && MultipOrCentMix<multBin[im+1]){
			jm = im; 
			break;
		}
	}

	if(TMath::Abs(Zvtx)>=10. || MultipOrCentMix>poolmax || MultipOrCentMix < poolmin){
      AliInfo(Form("pp Event with Zvertex = %.2f cm and multiplicity = %.0f out of pool bounds, SKIPPING",Zvtx,MultipOrCentMix));
      return;
	}

	if(MC==0)
    {
      if(d==3){
          // AliInfo("Bucla de reconstruite!!!!!!!!");
        AliEventPool* pool1 = fPoolMgr1->GetEventPool(MultipOrCentMix, Zvtx);
        if (!pool1){
          AliInfo(Form("No pool found for mult/cent = %f, zVtx = %f", MultipOrCentMix, Zvtx));
          return;
        }
  //       pool1->PrintInfo();
  // 		pool1->SetTargetEvents(5);

        if (pool1->IsReady() || pool1->GetCurrentNEvents() > 5){
          //	AliInfo("Pool este Ready!!!!!!!!");
          Int_t nMix = pool1->GetCurrentNEvents();
          for (Int_t jMix=0; jMix<nMix; jMix++){
            TObjArray* mixEvents = pool1->GetEvent(jMix);
            AliBasicParticle* trigger = (AliBasicParticle*)selectedArray->At(0);
            if(!trigger) continue;
            Double_t ptL  = trigger->Pt();
            Double_t phiL = trigger->Phi();
            Double_t etaL = trigger->Eta();
            if(trigger->Pt()<1.0 || trigger->Pt()>2.0) continue;
            if( TMath::Abs(trigger->Eta())> 0.8 ) continue;
        
            for (Int_t j=0; j<mixEvents->GetEntriesFast(); j++){
          //		 AliInfo("In a doua bucla!!!!!!");
              AliBasicParticle* associate = (AliBasicParticle*) mixEvents->At(j);
              if(!associate)continue;
              if(associate->Pt()<1.0 || associate->Pt()>2.0) continue;
              if( TMath::Abs(associate->Eta())> 0.8 ) continue;
              Double_t ptAs= associate->Pt();
              Double_t phiAs= associate->Phi();
              Double_t etaAs= associate->Eta();
              Double_t dPhi(-999.), dEta(-999.);
              if(ptL>ptAs && jm>-1 && ptL>=1.0 && ptL<=2.0){                     
                dPhi = RangePhi(phiL-phiAs);
                dEta=etaL-etaAs;
                Int_t bin=(37+jm);
              //	  AliInfo("Inainte de fill!!!!!!!!");
                ((TH2*)fHistosQA->At(bin))->Fill(dEta, dPhi); 
              }
            }
          }
        }
          TObjArray* tracksClone1 = new TObjArray;
          tracksClone1->SetOwner(kTRUE);
          for (Int_t i=0; i<selectedArray->GetEntriesFast(); i++){
            AliBasicParticle* particle1 = (AliBasicParticle*) selectedArray->At(i);
            tracksClone1->AddLast(particle1);
          }
            pool1->UpdatePool(tracksClone1);
      }
      if(d==6){
          // AliInfo("Bucla de reconstruite!!!!!!!!");
        AliEventPool* pool2 = fPoolMgr2->GetEventPool(MultipOrCentMix, Zvtx);
        if (!pool2){
          AliInfo(Form("No pool found for mult/cent = %f, zVtx = %f", MultipOrCentMix, Zvtx));
          return;
        }
  //       pool2->PrintInfo();
  // 		pool2->SetTargetEvents(5);

        if (pool2->IsReady() || pool2->GetCurrentNEvents() > 5){
          //	AliInfo("Pool este Ready!!!!!!!!");
          Int_t nMix = pool2->GetCurrentNEvents();
          for (Int_t jMix=0; jMix<nMix; jMix++){
            TObjArray* mixEvents = pool2->GetEvent(jMix);
            AliBasicParticle* trigger = (AliBasicParticle*)selectedArray->At(0);
            if(!trigger) continue;
            Double_t ptL  = trigger->Pt();
            Double_t phiL = trigger->Phi();
            Double_t etaL = trigger->Eta();
            if(trigger->Pt()<1.0 || trigger->Pt()>2.0) continue;
            if( TMath::Abs(trigger->Eta())> 0.8 ) continue;
        
            for (Int_t j=0; j<mixEvents->GetEntriesFast(); j++){
          //		 AliInfo("In a doua bucla!!!!!!");
              AliBasicParticle* associate = (AliBasicParticle*) mixEvents->At(j);
              if(!associate)continue;
              if(associate->Pt()<1.0 || associate->Pt()>2.0) continue;
              if( TMath::Abs(associate->Eta())> 0.8 ) continue;
              Double_t ptAs= associate->Pt();
              Double_t phiAs= associate->Phi();
              Double_t etaAs= associate->Eta();
              Double_t dPhi(-999.), dEta(-999.);
              if(ptL>ptAs && jm>-1 && ptL>=1.0 && ptL<=2.0){                     
                dPhi = RangePhi(phiL-phiAs);
                dEta=etaL-etaAs;
                Int_t bin=(49+jm);
              //	  AliInfo("Inainte de fill!!!!!!!!");
                ((TH2*)fHistosQA->At(bin))->Fill(dEta, dPhi); 
              }
            }
          }
        }
          TObjArray* tracksClone2 = new TObjArray;
          tracksClone2->SetOwner(kTRUE);
          for (Int_t i=0; i<selectedArray->GetEntriesFast(); i++){
            AliBasicParticle* particle2 = (AliBasicParticle*) selectedArray->At(i);
            tracksClone2->AddLast(particle2);
          }
            pool2->UpdatePool(tracksClone2);
      }
      if(d==9){
          // AliInfo("Bucla de reconstruite!!!!!!!!");
        AliEventPool* pool3 = fPoolMgr3->GetEventPool(MultipOrCentMix, Zvtx);
        if (!pool3){
          AliInfo(Form("No pool found for mult/cent = %f, zVtx = %f", MultipOrCentMix, Zvtx));
          return;
        }
  //       pool3->PrintInfo();
  // 		pool3->SetTargetEvents(5);

        if (pool3->IsReady() || pool3->GetCurrentNEvents() > 5){
          //	AliInfo("Pool este Ready!!!!!!!!");
          Int_t nMix = pool3->GetCurrentNEvents();
          for (Int_t jMix=0; jMix<nMix; jMix++){
            TObjArray* mixEvents = pool3->GetEvent(jMix);
            AliBasicParticle* trigger = (AliBasicParticle*)selectedArray->At(0);
            if(!trigger) continue;
            Double_t ptL  = trigger->Pt();
            Double_t phiL = trigger->Phi();
            Double_t etaL = trigger->Eta();
            if(trigger->Pt()<1.0 || trigger->Pt()>2.0) continue;
            if( TMath::Abs(trigger->Eta())> 0.8 ) continue;
        
            for (Int_t j=0; j<mixEvents->GetEntriesFast(); j++){
          //		 AliInfo("In a doua bucla!!!!!!");
              AliBasicParticle* associate = (AliBasicParticle*) mixEvents->At(j);
              if(!associate)continue;
              if(associate->Pt()<1.0 || associate->Pt()>2.0) continue;
              if( TMath::Abs(associate->Eta())> 0.8 ) continue;
              Double_t ptAs= associate->Pt();
              Double_t phiAs= associate->Phi();
              Double_t etaAs= associate->Eta();
              Double_t dPhi(-999.), dEta(-999.);
              if(ptL>ptAs && jm>-1 && ptL>=1.0 && ptL<=2.0){                     
                dPhi = RangePhi(phiL-phiAs);
                dEta=etaL-etaAs;
                Int_t bin=(61+jm);
              //	  AliInfo("Inainte de fill!!!!!!!!");
                ((TH2*)fHistosQA->At(bin))->Fill(dEta, dPhi); 
              }
            }
          }
        }
          TObjArray* tracksClone3 = new TObjArray;
          tracksClone3->SetOwner(kTRUE);
          for (Int_t i=0; i<selectedArray->GetEntriesFast(); i++){
            AliBasicParticle* particle3 = (AliBasicParticle*) selectedArray->At(i);
            tracksClone3->AddLast(particle3);
          }
            pool3->UpdatePool(tracksClone3);
      }
    }
    
    if(HasMCdata() && MC==1){
      if(d==3){
        AliEventPool* poolMC1 = fPoolMgrMC1->GetEventPool(MultipOrCentMix, Zvtx);
        if (!poolMC1){
          AliInfo(Form("No pool found for mult/cent = %f, zVtx = %f", MultipOrCentMix, Zvtx));
          return;
        }
//         poolMC1->PrintInfo();
  // 		poolMC1->SetTargetEvents(5);

        if (poolMC1->IsReady() || poolMC1->GetCurrentNEvents() > 5){
          Int_t nMix = poolMC1->GetCurrentNEvents();
          for (Int_t jMix=0; jMix<nMix; jMix++){
            TObjArray* mixEvents = poolMC1->GetEvent(jMix);
            AliBasicParticle* trigger = (AliBasicParticle*)selectedArray->At(0);
            if(!trigger) continue;
            Double_t ptL  = trigger->Pt();
            Double_t phiL = trigger->Phi();
            Double_t etaL = trigger->Eta();
            if(trigger->Pt()<1.0 || trigger->Pt()>2.0) continue;
            if( TMath::Abs(trigger->Eta())> 0.8 ) continue;
        
            for (Int_t j=0; j<mixEvents->GetEntriesFast(); j++){
              AliBasicParticle* associate = (AliBasicParticle*) mixEvents->At(j);
              if(!associate)continue;
              if(associate->Pt()<1.0 || associate->Pt()>2.0) continue;
              if( TMath::Abs(associate->Eta())> 0.8 ) continue;
              Double_t ptAs= associate->Pt();
              Double_t phiAs= associate->Phi();
              Double_t etaAs= associate->Eta();
              Double_t dPhi(-999.), dEta(-999.);
              if(ptL>ptAs && jm>-1 && ptL>=1.0 && ptL<=2.0){                     
                dPhi = RangePhi(phiL-phiAs);
                dEta=etaL-etaAs;
                Int_t bin=(109+jm);
                ((TH2*)fHistosQA->At(bin))->Fill(dEta, dPhi); 
              }
            }
          }
        }

          TObjArray* tracksCloneMC1 = new TObjArray;
          tracksCloneMC1->SetOwner(kTRUE);
          for (Int_t i=0; i<selectedArray->GetEntriesFast(); i++){
            AliBasicParticle* particleMC1 = (AliBasicParticle*) selectedArray->At(i);
            tracksCloneMC1->AddLast(particleMC1);
          }
            poolMC1->UpdatePool(tracksCloneMC1);
      }
      if(d==6){
        AliEventPool* poolMC2 = fPoolMgrMC2->GetEventPool(MultipOrCentMix, Zvtx);
        if (!poolMC2){
          AliInfo(Form("No pool found for mult/cent = %f, zVtx = %f", MultipOrCentMix, Zvtx));
          return;
        }
//         poolMC2->PrintInfo();
  // 		poolMC2->SetTargetEvents(5);

        if (poolMC2->IsReady() || poolMC2->GetCurrentNEvents() > 5){
          Int_t nMix = poolMC2->GetCurrentNEvents();
          for (Int_t jMix=0; jMix<nMix; jMix++){
            TObjArray* mixEvents = poolMC2->GetEvent(jMix);
            AliBasicParticle* trigger = (AliBasicParticle*)selectedArray->At(0);
            if(!trigger) continue;
            Double_t ptL  = trigger->Pt();
            Double_t phiL = trigger->Phi();
            Double_t etaL = trigger->Eta();
            if(trigger->Pt()<1.0 || trigger->Pt()>2.0) continue;
            if( TMath::Abs(trigger->Eta())> 0.8 ) continue;
        
            for (Int_t j=0; j<mixEvents->GetEntriesFast(); j++){
              AliBasicParticle* associate = (AliBasicParticle*) mixEvents->At(j);
              if(!associate)continue;
              if(associate->Pt()<1.0 || associate->Pt()>2.0) continue;
              if( TMath::Abs(associate->Eta())> 0.8 ) continue;
              Double_t ptAs= associate->Pt();
              Double_t phiAs= associate->Phi();
              Double_t etaAs= associate->Eta();
              Double_t dPhi(-999.), dEta(-999.);
              if(ptL>ptAs && jm>-1 && ptL>=1.0 && ptL<=2.0){                     
                dPhi = RangePhi(phiL-phiAs);
                dEta=etaL-etaAs;
                Int_t bin=(121+jm);
                ((TH2*)fHistosQA->At(bin))->Fill(dEta, dPhi); 
              }
            }
          }
        }

          TObjArray* tracksCloneMC2 = new TObjArray;
          tracksCloneMC2->SetOwner(kTRUE);
          for (Int_t i=0; i<selectedArray->GetEntriesFast(); i++){
            AliBasicParticle* particleMC2 = (AliBasicParticle*) selectedArray->At(i);
            tracksCloneMC2->AddLast(particleMC2);
          }
            poolMC2->UpdatePool(tracksCloneMC2);
      }
      if(d==9){
        AliEventPool* poolMC3 = fPoolMgrMC3->GetEventPool(MultipOrCentMix, Zvtx);
        if (!poolMC3){
          AliInfo(Form("No pool found for mult/cent = %f, zVtx = %f", MultipOrCentMix, Zvtx));
          return;
        }
//         poolMC3->PrintInfo();
  // 		poolMC3->SetTargetEvents(5);

        if (poolMC3->IsReady() || poolMC3->GetCurrentNEvents() > 5){
          Int_t nMix = poolMC3->GetCurrentNEvents();
          for (Int_t jMix=0; jMix<nMix; jMix++){
            TObjArray* mixEvents = poolMC3->GetEvent(jMix);
            AliBasicParticle* trigger = (AliBasicParticle*)selectedArray->At(0);
            if(!trigger) continue;
            Double_t ptL  = trigger->Pt();
            Double_t phiL = trigger->Phi();
            Double_t etaL = trigger->Eta();
            if(trigger->Pt()<1.0 || trigger->Pt()>2.0) continue;
            if( TMath::Abs(trigger->Eta())> 0.8 ) continue;
        
            for (Int_t j=0; j<mixEvents->GetEntriesFast(); j++){
              AliBasicParticle* associate = (AliBasicParticle*) mixEvents->At(j);
              if(!associate)continue;
              if(associate->Pt()<1.0 || associate->Pt()>2.0) continue;
              if( TMath::Abs(associate->Eta())> 0.8 ) continue;
              Double_t ptAs= associate->Pt();
              Double_t phiAs= associate->Phi();
              Double_t etaAs= associate->Eta();
              Double_t dPhi(-999.), dEta(-999.);
              if(ptL>ptAs && jm>-1 && ptL>=1.0 && ptL<=2.0){                     
                dPhi = RangePhi(phiL-phiAs);
                dEta=etaL-etaAs;
                Int_t bin=(133+jm);
                ((TH2*)fHistosQA->At(bin))->Fill(dEta, dPhi); 
              }
            }
          }
        }

          TObjArray* tracksCloneMC3 = new TObjArray;
          tracksCloneMC3->SetOwner(kTRUE);
          for (Int_t i=0; i<selectedArray->GetEntriesFast(); i++){
            AliBasicParticle* particleMC3 = (AliBasicParticle*) selectedArray->At(i);
            tracksCloneMC3->AddLast(particleMC3);
          }
            poolMC3->UpdatePool(tracksCloneMC3);
      }
    }
}

//________________________________________________________
Bool_t AliMESppColTask::BuildQAHistos()
{
  fHistosQA = new TList(); fHistosQA->SetOwner(kTRUE);
  
  
// used for scaling
  const Int_t ndimNoEvts(7);
  const Int_t cldNbinsNoEvts[ndimNoEvts]   = {5, 150, 30, 30, 150, 30, 30};
  const Double_t cldMinNoEvts[ndimNoEvts]  = {-0.5, 0.5, 0., 0., 0.5, 0., 0.}, cldMaxNoEvts[ndimNoEvts]  = {4.5, 150.5, 1., 1., 150.5, 1., 1.};
  THnSparseD *hNoEvts = new THnSparseD("NoEvts","NoEvts;step;combined 0.8;dirNeg;sfericity; MCmultiplicity;MCdirNeg; MCsfericity;",ndimNoEvts, cldNbinsNoEvts, cldMinNoEvts, cldMaxNoEvts);
  hNoEvts->GetAxis(0)->SetBinLabel(1, "Tender OK");
  hNoEvts->GetAxis(0)->SetBinLabel(2, "Pile-up Rejection");
  hNoEvts->GetAxis(0)->SetBinLabel(3, "Vertex Cut");
  hNoEvts->GetAxis(0)->SetBinLabel(4, "Analyzed");
  fHistosQA->AddAt(hNoEvts, 0);

//histos  
  const Int_t nMult(12);
  Int_t selBin[4]={3, 6, 9};
  Int_t multBin[nMult+1] = {1, 4, 7, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80};
  TH2F *hSE[36] = {NULL};
  TH2F *hME[36] = {NULL};
  for(Int_t sel(0); sel<3; sel++){
    for(Int_t im(0); im<nMult; im++){
        hSE[im+12*sel] = new TH2F(Form("hS%02d", im+12*sel), Form("Mult[%2d-%2d]Sfer[%d]", multBin[im], multBin[im+1]-1, selBin[sel]), 72, -1.5, 1.5, 60, -0.5*TMath::Pi(), 1.5*TMath::Pi());
        fHistosQA->AddAt(hSE[im+12*sel], im+1+12*sel);
        hME[im+12*sel] = new TH2F(Form("hE%02d", im+12*sel), Form("MultE[%2d-%2d]Sfer[%d]", multBin[im], multBin[im+1]-1, selBin[sel]), 72, -1.5, 1.5, 60, -0.5*TMath::Pi(), 1.5*TMath::Pi());
        fHistosQA->AddAt(hME[im+12*sel], im+37+12*sel);
    }
  }
    
  TH2F *hSEMC[36] = {NULL};
  TH2F *hMEMC[36] = {NULL};
  for(Int_t sel(0); sel<3; sel++){
    for(Int_t im(0); im<nMult; im++){
        hSEMC[im+12*sel] = new TH2F(Form("hSMC%02d", im+12*sel), Form("MultSMC[%2d-%2d]Sfer[%d]", multBin[im], multBin[im+1]-1, selBin[sel]), 72, -1.5, 1.5, 60, -0.5*TMath::Pi(), 1.5*TMath::Pi());
        fHistosQA->AddAt(hSEMC[im+12*sel], im+12*sel+73);
        
        hMEMC[im+12*sel] = new TH2F(Form("hEMC%02d", im+12*sel), Form("MultEMC[%2d-%2d]Sfer[%d]", multBin[im], multBin[im+1]-1, selBin[sel]), 72, -1.5, 1.5, 60, -0.5*TMath::Pi(), 1.5*TMath::Pi());
        fHistosQA->AddAt(hMEMC[im+12*sel], im+12*sel+109);
      }
  }
  
  //
  const Int_t ndimTrk(4);
  const Int_t cldNbinsTrk[ndimTrk]   = {150, 30, 30, 50};
  const Double_t cldMinTrk[ndimTrk]  = {0.5, 0., 0.,  0.},
					  cldMaxTrk[ndimTrk]  = {150.5, 1., 1., 5.};
  THnSparseD *hTrk = new THnSparseD("infoTrk","infoTrk;combined0.8;dirNeg;sfericity;p_{T};",ndimTrk, cldNbinsTrk, cldMinTrk, cldMaxTrk);
  fHistosQA->AddAt(hTrk, 145);
  
  THnSparseD *hMCTrk = new THnSparseD("infoMCTrk","infoMCTrk;MCmult;MCdirNeg;MCsfericity;MCp_{T};",ndimTrk, cldNbinsTrk, cldMinTrk, cldMaxTrk);
  fHistosQA->AddAt(hMCTrk, 146);
  return kTRUE;
}