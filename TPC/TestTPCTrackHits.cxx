/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/// \class TestTPCTrackHits
/// Macro to compare TClonesArray hits with interpolated hits
/// ConvertHits1 read
///
/// \author MI

#include "alles.h"
#include "AliObjectArray.h"
#include "AliTPCTrackHits.h"
#include "AliArrayBranch.h"
#include "TestTPCTrackHits.h"
 
/// \cond CLASSIMP
ClassImp(AliTPChitD)
/// \endcond
 
void CompareHits(TClonesArray * arr, AliTPCTrackHits * myhits,  Bool_t debug, TClonesArray *arrd=0);
AliTPCTrackHits * MakeTrack(TClonesArray * arr, TClonesArray * arrp, AliTPCTrackHits *myhits);

void ConvertHits(const char * benchmark, Bool_t debug)
{      

  /// convert - not parametrised hits stored in Clones array
  /// to

  TFile f("galice.root","update");  
  TClonesArray *arr = new TClonesArray("AliTPChit",100);
  TTree * treeCl =(TTree*)f.Get("TreeH0");
  TBranch *branch = treeCl->GetBranch("TPC"); 
  AliTPCTrackHits *  myhits = new AliTPCTrackHits;
  branch->SetAddress(&arr);
  //particle  
  TTree * treeP = (TTree*)f.Get("TreeK0");
  TBranch *branchP = treeP->GetBranch("Particles"); 
  TClonesArray *arrp = new TClonesArray("TParticle",100);
  branchP->SetAddress(&arrp);
  branchP->GetEvent(0);
  //
  TFile f2("treeh.root","recreate");
  f2.SetCompressionLevel(10);
  TTree * treeh2 = new TTree("TreeTPCH","TreeTPCH");
  //treeh2->AliBranch("TPC","AliTPCTrackHits", &myhits,4096);
  treeh2->GetListOfBranches()->Add(new AliObjectBranch("TPC","AliTPCTrackHits", 
	 &myhits,treeh2,4096,1));

  TFile f3("treehcl.root","recreate");
  f3.SetCompressionLevel(2);
  TTree * treeh3 = new TTree("TreeTPCH","TreeTPCH");
  treeh3->Branch("TPC", &arr,64000,100);

  gBenchmark->Start(benchmark); 
  Int_t trsize = (Int_t)treeCl->GetEntries();
  for (Int_t i=0;i<300;i++){    
    arr->Clear(); 
    Int_t size = branch->GetEvent(i);     
    //if (size>0){
    if ((size>0)&& arr->GetEntriesFast()>0) {    
      f.cd();
      myhits->Clear();
      myhits =MakeTrack(arr,arrp,myhits);            
      CompareHits(arr,myhits,debug);
      f2.cd();
      treeh2->Fill();
      f3.cd();
      treeh3->Fill();      
      if ((i%10)==0) cout<<i<<"\n"; 
    }    
  }
  f3.cd();
  treeh3->Write();
  f2.cd();
  treeh2->Write();
  f2.Close();
  f.Close();
  delete myhits;
  gBenchmark->Stop(benchmark);
  gBenchmark->Show(benchmark); 
}  



void CompareHits(const char * benchmark, Bool_t debug)
{        
  /// compare persistent hits

  TFile f2("treeh.root");
  TTree * treeh2 = (TTree*)f2.Get("TreeTPCH");
  AliTPCTrackHits *  myhits = new AliTPCTrackHits ; 
  AliObjectBranch *mybranch = (AliObjectBranch*)treeh2->GetBranch("TPC"); 
  mybranch->SetAddress(&myhits);

    
  TClonesArray *arr = new TClonesArray("AliTPChit",100);
  TFile f3("treehcl.root");
  TTree * treeh3 = (TTree*)f3.Get("TreeTPCH");
  TBranch *branch = treeh3->GetBranch("TPC"); 
  branch->SetAddress(&arr);
   
  cout<<"Lets go!\n"; 
  gBenchmark->Start(benchmark);
  Int_t trsize = (Int_t)treeh3->GetEntries();
  for (Int_t i=0;i<300;i++){    
    Int_t size = branch->GetEvent(i);     
    mybranch->GetEvent(i);
    //if (arr){
    if ((arr->GetEntriesFast()>0) && size>0) {
      CompareHits(arr,myhits,debug);
      if ((i%10)==0) cout<<i<<"\n"; 
    }    
  }
  gBenchmark->Stop(benchmark);
  gBenchmark->Show(benchmark); 
}  


void CompareHitsG(const char * benchmark, Bool_t debug)
{        
  /// compare persistent hits

  TFile f2("galice.root");
  TTree * treeh2 = (TTree*)f2.Get("TreeH0");
  AliTPCTrackHits *  myhits = new AliTPCTrackHits ; 
  AliObjectBranch *mybranch = (AliObjectBranch*)treeh2->GetBranch("TPC2"); 
  mybranch->SetAddress(&myhits);

    
  TClonesArray *arr = new TClonesArray("AliTPChit",100);
  //TFile f3("treehcl.root");
  //TTree * treeh3 = (TTree*)f3.Get("TreeTPCH");
  TBranch *branch = treeh2->GetBranch("TPC"); 
  branch->SetAddress(&arr);

  TFile f3("treehdelta.root","recreate");
  f3.SetCompressionLevel(2);
  TTree * treeh3 = new TTree("DelataH","DeltaH");
  TClonesArray *arrd = new TClonesArray("AliTPChitD",100);
  treeh3->Branch("TPC", &arrd,64000,100);
   
  cout<<"Lets go!\n"; 
  gBenchmark->Start(benchmark);
  Int_t trsize = treeh2->GetEntries();
  for (Int_t i=0;i<trsize;i++){    
    Int_t size = branch->GetEvent(i);     
    mybranch->GetEvent(i);
    //if (arr){
    if ((arr->GetEntriesFast()>0) && size>0) {
      CompareHits(arr,myhits,debug,arrd);      
    }  
    if ((i%10)==0) cout<<i<<"\n";
    treeh3->Fill();         
  }
  treeh3->Write();
  f3.Close();
  gBenchmark->Stop(benchmark);
  gBenchmark->Show(benchmark); 
}  


 
   
AliTPCTrackHits * MakeTrack(TClonesArray * arr, TClonesArray * arrp, AliTPCTrackHits *myhits)  
{
  /// make track wit hits  according

  AliTPChit * hit;
  //  AliTPCTrackHits * myhits= new AliTPCTrackHits;
  myhits->SetHitPrecision(0.002);
  myhits->SetStepPrecision(0.003);  
  myhits->SetMaxDistance(100);  
  myhits->FlushHitStack(kTRUE);
  for (Int_t i=0;i<arr->GetEntriesFast();i++){
    hit = (AliTPChit*)arr->At(i);    
    if (hit){
      TParticle *p = (TParticle*)arrp->At(hit->GetTrack());
      Float_t momentum = TMath::Sqrt(p->Px()*p->Px()+p->Py()*p->Py());
      Float_t ran= 100.*gRandom->Rndm();
      if (ran<1.)  myhits->FlushHitStack(kTRUE);
      if (momentum<0.01) {//if small momentum - not precise
	myhits->SetHitPrecision(0.05);
	myhits->AddHitKartez(hit->fSector,hit->GetTrack(), hit->X(), hit->Y(),
			     hit->Z(), hit->fQ+1000);
      }
      else {
	myhits->SetHitPrecision(0.002);
	myhits->AddHitKartez(hit->fSector,hit->GetTrack(), hit->X(), hit->Y(),
			     hit->Z(), hit->fQ);
      }
    }
  }
  myhits->FlushHitStack();
  return myhits;  
}



void CompareHits(TClonesArray * arr, AliTPCTrackHits * myhits, Bool_t debug, TClonesArray *arrd)
{     
  /// if debug option  kTRUE
  /// compare hits and write result  to the stdoutput

  AliTPChit * hit, *hit2;
  if (arrd) arrd->Clear();

  for (Int_t i=0;i<arr->GetEntriesFast();i++){
    hit = (AliTPChit*)arr->At(i);  
    if (hit){
      if (i==0) myhits->First();
      else myhits->Next();
      hit2 = myhits->GetHit();   
      if (!hit2) {
	hit2=0;
	if (hit) cout<<"Error _ hits "<<i<<"didn't find\n";
      }

      if (hit&&hit2){
	// if (hit){
	
	AliTrackHitsParam *param= myhits->GetParam();
	AliHitInfo * info = myhits->GetHitInfo();
	//	
	if (arrd) {
	  TClonesArray &larrd = *arrd;
	  AliTPChitD * h = new(larrd[i]) AliTPChitD;
	  h->SetX(hit->X());
	  h->SetY(hit->Y());
	  h->SetZ(hit->Z());
	  h->SetTrack(hit->GetTrack());
	  h->fQ = hit->fQ;
	  h->fSector = hit->fSector;
	  AliTPChit * hitd = h->GetDelta();
	  hitd->SetX(hit->X()-hit2->X());
	  hitd->SetY(hit->Y()-hit2->Y());
	  hitd->SetZ(hit->Z()-hit2->Z());
	  hitd->SetTrack(hit->GetTrack()-hit2->GetTrack());
	  hitd->fQ = hit->fQ-hit2->fQ;
	  hitd->fSector = hit->fSector-hit2->fSector;
	}

	if (debug){
	  Float_t dd = 
	    TMath::Sqrt(
			(hit->X()-hit2->X())*(hit->X()-hit2->X())+
			(hit->Y()-hit2->Y())*(hit->Y()-hit2->Y())+
			(hit->Z()-hit2->Z())*(hit->Z()-hit2->Z())); 
	  printf("C1\t%d\t%d\t%d\t%d\t",
		 hit->fSector,hit2->fSector,hit->GetTrack(),hit2->GetTrack());
	  printf("%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t",
		 hit->X(),hit2->X(), hit->Y(),hit2->Y(),hit->Z(),hit2->Z());
	  printf("%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t",
		 dd, param->fR,param->fZ,param->fFi,param->fAn,
		 param->fAd,param->fTheta);
	  printf("%d\t%d\t%d\n",
		 (Int_t)info->fHitDistance, (Int_t)hit->fQ, (Int_t)hit2->fQ);
	}
      }        
    }
  }
} 




//Filter for paw visualisation of results
//
//sed filter  
 //sed '/*/d' dd.txt | sed '/C/d'|  cat -n > out.txt
 // sed -n '/C2/p' dd.txt | sed 's/C1//' | cat -n >out2.txt 
 // sed -n '/C3/p' dd.txt | sed 's/C2//' | cat -n >out3.txt 

//filter   
//myfilter C1 dd.txt >out1.txt
//myfilter C2 dd.txt >out2.txt
//myfilter C3 dd.txt >out3.txt
// sed -n 1,50000p
// sed -n 1,50000p dd.txt | myfilter C1  | sed -n 1,50000p >out1.txt

//myfilter C1 dd.txt | sed -n 1,50000p >out1.txt
//myfilter C2 dd.txt | sed -n 1,50000p >out2.txt
//myfilter C3 dd.txt | sed -n 1,50000p >out3.txt 
/*
paw visualisation
 Nt/Create 1 'Test of Match' 21 ! ! event v1 v2 t1 t2  x1 x2 y1  y2 z1 z2 dd r z fi an ad theta dist q1 q2
  Nt/Read 1 out1.txt 
 nt 
 Nt/Create 2 'Test of Match2' 13 ! ! event stacksize i r  z fi alfa theta dr1 ddr ddz ddrfi dd
  Nt/Read 2 out2.txt
 
 Nt/Create 3 'Test of Match2' 2 ! ! event stacksize 
  Nt/Read 3 out3.txt     
*/ 

void   Fit2(Double_t fSumY, Double_t fSumYX, Double_t fSumYX2,
	    Double_t fSumX,  Double_t fSumX2, Double_t fSumX3, 
	    Double_t fSumX4, Int_t n,
	    Double_t &a, Double_t &b, Double_t &c)
{
  /// recalc parameters not fixing origin point

  Double_t det = 
    n* (fSumX2*fSumX4-fSumX3*fSumX3) -
    fSumX*      (fSumX*fSumX4-fSumX3*fSumX2)+
    fSumX2*     (fSumX*fSumX3-fSumX2*fSumX2);
    
  if (TMath::Abs(det)>0){    
    a = 
      (fSumY * (fSumX2*fSumX4-fSumX3*fSumX3)-
       fSumX *(fSumYX*fSumX4-fSumYX2*fSumX3)+
       fSumX2*(fSumYX*fSumX3-fSumYX2*fSumX2))/det; 
    b=
      (n*(fSumYX*fSumX4-fSumX3*fSumYX2)-
      fSumY*(fSumX*fSumX4-fSumX3*fSumX2)+
      fSumX2*(fSumX*fSumYX2-fSumYX*fSumX2))/det;
    c=
      (n*(fSumX2*fSumYX2-fSumYX*fSumX3)-
       fSumX*(fSumX*fSumYX2-fSumYX*fSumX2)+
       fSumY*(fSumX*fSumX3-fSumX2*fSumX2))/det;  
    cout<<a<<"\t"<<b<<"\t"<<c<<"\n";
  }
}

void TestFit(Float_t a, Float_t b, Float_t c, Int_t n)
{
  Double_t fSumY,fSumYX,fSumYX2,fSumX, fSumX2,fSumX3, fSumX4;
  fSumY = fSumYX = fSumYX2 = fSumX =  fSumX2 = fSumX3 =  fSumX4 =0;
  Double_t a2,b2,c2;
  for (Int_t i=0;i<n; i++){
    Float_t  x = gRandom->Rndm();
    Float_t  y = a+b*x+c*x*x;
    fSumY+=y;
    fSumYX+=y*x;
    fSumYX2+=y*x*x;
    fSumX += x;
    fSumX2 += x*x;
    fSumX3 += x*x*x;
    fSumX4 += x*x*x*x;
  }
  Fit2(fSumY,fSumYX,fSumYX2, fSumX,fSumX2,fSumX3,fSumX4,n,a2,b2,c2);  
}
