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
  ,fPoolMgr(0x0)
{
  //
  // Constructor
  //
}

//________________________________________________________________________
AliMESppColTask::AliMESppColTask(const char *name)
  : AliMESbaseTask(name)
  ,fPoolMgr(0x0)
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
  if(fPoolMgr) {delete fPoolMgr; fPoolMgr=0;}
}

//________________________________________________________________________
void AliMESppColTask::UserCreateOutputObjects()
{
  //define user data containers
  
  AliMESbaseTask::UserCreateOutputObjects();
  
  Bool_t DefPool = DefineMixedEventPool();
  if(!DefPool){
    AliInfo("UserCreateOutput: Pool is not defined properly");
    return;
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
	
	Double_t vec_hNoEvts[7]; // vector used to fill hNoEvts
	THnSparseD *hNoEvts = (THnSparseD*)fHistosQA->At(0);
    
	Double_t mult_comb08 = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb);// combined multiplicity with |eta| < 0.8
	Double_t mult_V0M = fEvInfo->GetMultiplicity(AliMESeventInfo::kV0M); // V0M percentile
	// event shape for data (from ESD)
	Double_t directivity_plus = fEvInfo->GetEventShape()->GetDirectivity(1);
	Double_t directivity_minus = fEvInfo->GetEventShape()->GetDirectivity(0);
	Double_t sfer = fEvInfo->GetEventShape()->GetSphericity();
// 	  AliInfo(Form("Sfer=%f", sfer));
	// select events with both dirs in the same interval
// 	const Int_t lenght = 4;
// 	Double_t intervals[lenght] = {0., 0.3, 0.6, 0.9};
	// NOTE: the intervals are considered half-closed: (a,b]
/*	if( (directivity_plus < intervals[0]) || (directivity_plus > intervals[lenght-1]) ) return;
	Int_t first = -1;
	for(Int_t i=1; i<lenght; i++){
		if(directivity_plus <= intervals[i]){
			first = i;
			break;
		}
	}*/
	Double_t directivity = (directivity_plus + directivity_minus) / 2.0;
//	if( (directivity_minus <= intervals[first-1]) || (directivity_minus > intervals[first]) ) return;
  if (directivity_minus <= 0.0 || directivity_minus > 0.3) return;
  if (directivity_plus <= 0.0 || directivity_plus > 0.3) return;
  if (mult_comb08<30 || mult_comb08>39) return;
	
	vec_hNoEvts[1] = mult_comb08; // combined multiplicity with |eta| < 0.8
	vec_hNoEvts[2] = directivity;
	vec_hNoEvts[3] = sfer; 

	
	
	// event multiplicity and shape for MC (from MC event)
	Double_t MC_directivity_plus = 0;
	Double_t MC_directivity_minus = 0;
	Double_t MC_directivity = 0;
	Double_t MC_mult_glob08 = 0;
	Double_t MC_mult_V0M = 0;
	Double_t MC_sfer = 0;
	if( HasMCdata() ){ // run only on MC
		MC_directivity_plus = fMCevInfo->GetEventShape()->GetDirectivity(1);
		MC_directivity_minus = fMCevInfo->GetEventShape()->GetDirectivity(0);
		MC_directivity = (MC_directivity_plus + MC_directivity_minus) / 2.0;
		MC_mult_glob08 = fMCevInfo->GetMultiplicity(AliMESeventInfo::kGlob08);
		MC_mult_V0M = fMCevInfo->GetMultiplicity(AliMESeventInfo::kV0M);
		MC_sfer = fMCevInfo->GetEventShape()->GetSphericity();
		vec_hNoEvts[4] = fMCevInfo->GetMultiplicity(AliMESeventInfo::kGlob08);
		vec_hNoEvts[5] = MC_directivity;
		vec_hNoEvts[6] = MC_sfer;
	}
	
	vec_hNoEvts[0] = 3.;
	hNoEvts->Fill(vec_hNoEvts);
	
	
	// ESD track loop
	AliMEStrackInfo *t(NULL), *tMC(NULL);
	const Int_t nMult(12);
	Int_t multBin[nMult+1] = {1, 4, 7, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80};
	
	TObjArray *selectedTracks=FindLeadingObjects(fTracks);
	if(!selectedTracks) return;
	selectedTracks->SetOwner(kTRUE);
	
    
// NOTE: the intervals are considered half-closed: (a,b]
    if(directivity_plus>0.0 && directivity_plus<=0.3 && directivity_minus>0.0 && directivity_minus<=0.3){
      FillCorrelationSE(mult_comb08, selectedTracks);
      FillCorrelationMixing(mult_comb08, fEvInfo->GetVertexZ(), 80., 0., selectedTracks);
	}

  
	if( HasMCdata() ){// run only on MC  
		TObjArray *selectedTracksMC=FindLeadingObjects(fMCtracks);
		if(!selectedTracksMC) return;
		selectedTracksMC->SetOwner(kTRUE);
    // NOTE: the intervals are considered half-closed: (a,b]
		if(MC_directivity_plus>0.0 && MC_directivity_plus<=0.3 && MC_directivity_minus>0.0 && MC_directivity_minus<=0.3){
			FillCorrelationSE(MC_mult_glob08, selectedTracksMC);
			FillCorrelationMixing(MC_mult_glob08, fMCevInfo->GetVertexZ(), 80., 0., selectedTracksMC);
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

TObjArray*  AliMESppColTask::FindLeadingObjects(TObjArray *obj)
{
	// Returns an array of charged particles ordered according to their pT.
	Int_t nTracks = obj->GetEntries();
	if( !nTracks ) return 0;
	// Define array of AliMEStrackInfo objects
	TObjArray* tracks = new TObjArray(nTracks);
	AliMEStrackInfo* part(NULL);
	// Loop over tracks
	for (Int_t ipart=0; ipart<nTracks; ++ipart) {
		if(!(part = (AliMEStrackInfo*)obj ->At(ipart))) continue;
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

Bool_t AliMESppColTask::DefineMixedEventPool()
{
  Int_t PoolMaxNEvents = 6000;
  Int_t PoolMinNTracks = 50000;
  
  Int_t NMultBins = 12;
  Double_t MultBins[] = { 1., 4., 7., 10., 15., 20., 25., 30., 40., 50., 60., 70., 80. };
  Int_t NzVtxBins = 5;
  Double_t ZvtxBins[] = {  -10., -5., -2.5, 2.5, 5.,  10. }; 
  
  fPoolMgr = new AliEventPoolManager(PoolMaxNEvents, PoolMinNTracks, NMultBins, MultBins, NzVtxBins, ZvtxBins);
  fPoolMgr->SetDebug(0);
  if(!fPoolMgr) return kFALSE;
  
  return kTRUE;

}
    
void AliMESppColTask::FillCorrelationSE(Double_t MultipOrCent, TObjArray*selectedArray)
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

//                                 if(ptAs>ptL) continue;
		Double_t dPhi(-999.), dEta(-999.);
		if(ptL>ptAs && jm>-1 && ptL>=1.0 && ptL<=2.0){                   
			dPhi = RangePhi(phiL-phiAs);
			dEta=etaL-etaAs;
			Int_t bin=(1+jm);
			((TH2*)fHistosQA->At(bin))->Fill(dEta, dPhi); 
      }
	}
}

//---------------------------------------------------------------------------------------
void AliMESppColTask::FillCorrelationMixing(Double_t MultipOrCentMix, Double_t Zvtx, Double_t poolmax, Double_t poolmin, TObjArray*selectedArray)
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

	AliEventPool* pool = fPoolMgr->GetEventPool(MultipOrCentMix, Zvtx);
	if (!pool){
		AliInfo(Form("No pool found for mult/cent = %f, zVtx = %f", MultipOrCentMix, Zvtx));
		return;
	}
// 	pool->PrintInfo();

	if (pool->IsReady() || pool->NTracksInPool() > 1000  || pool->GetCurrentNEvents() > 5){
// 		AliInfo("InBucla!!!!!!!");
      Int_t nMix = pool->GetCurrentNEvents();
      for (Int_t jMix=0; jMix<nMix; jMix++){
// 			AliInfo("InBucla2!!!!!!!");
        TObjArray* mixEvents = pool->GetEvent(jMix);
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
          Double_t dPhi(-999.), dEta(-999.), yy(-999.), yx(-999.);
          if(ptL>ptAs && jm>-1 && ptL>=1.0 && ptL<=2.0){                     
            dPhi = RangePhi(phiL-phiAs);
            dEta=etaL-etaAs;
            Int_t bin=(13+jm);
            ((TH2*)fHistosQA->At(bin))->Fill(dEta, dPhi); 
          }
        }
      }
    }

        TObjArray* tracksClone = new TObjArray;
        tracksClone->SetOwner(kTRUE);
        
        for (Int_t i=0; i<selectedArray->GetEntriesFast(); i++){
          AliBasicParticle* particle = (AliBasicParticle*) selectedArray->At(i);
          tracksClone->AddLast(particle);
        }
          pool->UpdatePool(tracksClone);
}

//________________________________________________________
Bool_t AliMESppColTask::BuildQAHistos()
{
  fHistosQA = new TList(); fHistosQA->SetOwner(kTRUE);
  
  
// used for scaling
  const Int_t ndimNoEvts(7);
  const Int_t cldNbinsNoEvts[ndimNoEvts]   = {4, 150, 30, 30, 150, 30, 30};
  const Double_t cldMinNoEvts[ndimNoEvts]  = {-0.5, 0.5, 0., 0., 0.5, 0., 0.}, cldMaxNoEvts[ndimNoEvts]  = {3.5, 150.5, 1., 1., 150.5, 1., 1.};
  THnSparseD *hNoEvts = new THnSparseD("NoEvts","NoEvts;step;combined 0.8;directivity;sfericity; MCmultiplicity;MCdirectivity; MCsfericity;",ndimNoEvts, cldNbinsNoEvts, cldMinNoEvts, cldMaxNoEvts);
  hNoEvts->GetAxis(0)->SetBinLabel(1, "Tender OK");
  hNoEvts->GetAxis(0)->SetBinLabel(2, "Pile-up Rejection");
  hNoEvts->GetAxis(0)->SetBinLabel(3, "Vertex Cut");
  hNoEvts->GetAxis(0)->SetBinLabel(4, "Analyzed");
  fHistosQA->AddAt(hNoEvts, 0);

//histos  
  const Int_t nMult(12);
  Int_t multBin[nMult+1] = {1, 4, 7, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80};
  TH2F *hSE[nMult] = {NULL};
  TH2F *hME[nMult] = {NULL};
  for(Int_t im(0); im<nMult; im++){
       hSE[im] = new TH2F(Form("hS%02d", im), Form("Mult [%2d-%2d]", multBin[im], multBin[im+1]-1), 33, -1.5, 1.5, 30, -0.5*TMath::Pi(), 1.5*TMath::Pi());
       fHistosQA->AddAt(hSE[im], im+1);
       hME[im] = new TH2F(Form("hE%02d", im), Form("MultE [%2d-%2d]", multBin[im], multBin[im+1]-1), 33, -1.5, 1.5, 30, -0.5*TMath::Pi(), 1.5*TMath::Pi());
       fHistosQA->AddAt(hME[im], im+13);
    }
  return kTRUE;
}