#include <Riostream.h>
#include "TObject.h"
#include "TMath.h"
#include "Database.h"
#include "TString.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "./Topology.h"
#include "TFile.h"
#include "./TopDatabase.h"
#include "AliITSUClusterPix.h"

ClassImp(TopDatabase)

TopDatabase::TopDatabase():TObject(),fN(0),fTotClusters(0),fThreshold(0),
  fOverThr(0),fNGroups(0),fNmax(1e5),fArrTopologies(),fArrHisto(){
}

TopDatabase::TopDatabase(TopDatabase &ogg):TObject(),fArrHisto(){
  fN=ogg.GetN();
  fArrTopologies.Expand(fN);
  for(int i=0; i<fN; i++){
    fArrTopologies.AddAt(new Topology(*(Topology*)(ogg.GetArrTopologies()->At(i))),i);
  }
  fTotClusters=ogg.GetTotClusters();
  fOverThr=ogg.GetOverThr();
  fNGroups=ogg.GetNGroups();
  fNmax=ogg.GetNmax();
}

TopDatabase::~TopDatabase(){
  fArrTopologies.Delete();
}

void TopDatabase::AccountTopology(const AliITSUClusterPix &cluster, Float_t dX, Float_t dZ, Float_t alpha, Float_t beta){

  TBits Top;
  Top.Clear();
  Int_t rs = cluster.GetPatternRowSpan();
  Int_t cs = cluster.GetPatternColSpan();
  for(Int_t ir=0;ir<rs;ir++)
    for(Int_t ic=0;ic<cs;ic++) 
      if(cluster.TestPixel(ir,ic)) Top.SetBitNumber(ir*cs+ic);
  Top.SetUniqueID((rs<<16)+cs);
  Bool_t newPatt = kTRUE;
  Bool_t Junk = kFALSE;
  fTotClusters++;
  Int_t indTop = -1;
  for(Int_t ip=0;ip<fN;ip++) {
    TBits pattOld = ((Topology*)fArrTopologies.At(ip))->GetPattern();
    if(pattOld==Top && pattOld.GetUniqueID()==Top.GetUniqueID()){
      newPatt = kFALSE;
      indTop = ip;
      break;
    }
  }
  if(newPatt){
    if(fN == fNmax){ //Junk bin
      Junk = kTRUE;
    }
    else {
      TBits* pt = new TBits(Top);
      this->ExpandDB(pt);
      indTop=fN-1;
    }
  }
  if(Junk==kTRUE) return;
  ((Topology*)fArrTopologies.At(indTop))->IncreaseCounts();
  TH2*  h0a = ((Topology*)fArrTopologies.At(indTop))->GetHxA();
  h0a->Fill(alpha, dX);
  TH2*  h1a = ((Topology*)fArrTopologies.At(indTop))->GetHzA();
  h1a->Fill(alpha, dZ);
  TH2*  h2a = ((Topology*)fArrTopologies.At(indTop))->GetHxB();
  h2a->Fill(beta, dX);
  TH2*  h3a = ((Topology*)fArrTopologies.At(indTop))->GetHzB();
  h3a->Fill(beta, dZ);
}

void TopDatabase::ExpandDB(TBits* patt){
  fN++;
  fArrTopologies.Expand(fN);
  Topology* top = new Topology(*patt);
  fArrTopologies.AddAt(top,fN-1); 
}

void TopDatabase::EndAndSort(Int_t mode){ 
  Topology::SetMode(mode);
  Int_t nPatterns = fN;
  fArrTopologies.Sort();
  Bool_t alreadyProcessed[nPatterns];
  TArrayF arrFreq;
  arrFreq.Set(nPatterns);
  TArrayI sortIndex;
  sortIndex.Set(nPatterns);
  for(Int_t i=0; i<nPatterns; i++){
    Topology* top = (Topology*)fArrTopologies.At(i);
    Float_t tempFreq = ((Float_t)(top->GetCounts()))/fTotClusters;
    top->SetFreq(tempFreq);
    arrFreq[i] = tempFreq;
    top->SetFlag(0);
    alreadyProcessed[i]=kFALSE;
  }
  TMath::Sort(nPatterns,arrFreq.GetArray(),sortIndex.GetArray());

  //Solving clashes dtoring the first 4 bytes of the topology
  for(Int_t i=0; i<nPatterns; i++){
    if(alreadyProcessed[i]==kTRUE) continue;
    else alreadyProcessed[i]=kTRUE;
    Topology* topFreq = (Topology*)fArrTopologies.At(sortIndex[i]);
    topFreq->SetPattID(i);
    Topology* top = (Topology*)fArrTopologies.At(i);
    Int_t refHash = top->GetHash();
    for(Int_t j=i+1; j<nPatterns; j++){
      Topology* top2comp = (Topology*)fArrTopologies.At(j);
      Int_t counter = 0;
      if(top2comp->GetHash()==refHash){
	alreadyProcessed[j]=kTRUE;
	if(counter==0){
	  counter++;
	  top->SetFlag(1);
	}
	top2comp->SetFlag(1);
      }
    }
  }
}

void TopDatabase::SetThreshold(Float_t thr){
  fThreshold=thr;
  Int_t nPatterns = fN;
  TArrayF arrFreq;
  arrFreq.Set(nPatterns);
  TArrayI sortIndex;
  sortIndex.Set(nPatterns);
  for(Int_t i=0; i<nPatterns; i++){
    Topology* top = (Topology*)fArrTopologies.At(i);
    Float_t tempFreq = top->GetFreq();
    arrFreq[i] = tempFreq;
  }
  TMath::Sort(nPatterns,arrFreq.GetArray(),sortIndex.GetArray());
  Int_t over=0;
  Float_t fr=0;
  for(Int_t j=0; j<fN; j++){
    fr=arrFreq[sortIndex[j]];
    if(fr<thr) break;
    else over++;
  }
  fOverThr=over;
}

void TopDatabase::SetThresholdCumulative(Float_t cumulative){
  if(cumulative<=0 || cumulative >=1) cumulative = 0.99;
  Float_t totFreq = 0.;
  Int_t nPatterns = fN;
  TArrayF arrFreq;
  arrFreq.Set(nPatterns);
  TArrayI sortIndex;
  sortIndex.Set(nPatterns);
  for(Int_t i=0; i<nPatterns; i++){
    Topology* top = (Topology*)fArrTopologies.At(i);
    Float_t tempFreq = top->GetFreq();
    arrFreq[i] = tempFreq;
  }
  TArrayF provvFreq;
  provvFreq.Set(fN);
  TMath::Sort(fN,arrFreq.GetArray(),sortIndex.GetArray());
  for(Int_t j=0; j<fN; j++){
    provvFreq[j]=arrFreq[sortIndex[j]];
  }
  Int_t over=0;
  while(totFreq < cumulative){
    totFreq+=provvFreq[over++];
  }
  fThreshold=provvFreq[--over];
  while(provvFreq[over]==fThreshold) over++;
  over++;
  fOverThr = over;
}

void TopDatabase::PrintDB(const char* output) const{
  ofstream o(output);
  o << "Number of topologies: " << fN << endl << endl;
  for(Int_t i=0; i<fN; i++){
    Topology* top = (Topology*)fArrTopologies.At(i);
    TBits patt = top->GetPattern();
    o << "               POSITION " << i << endl;
    o << Form("\nPattID %4d R:%d C: %d Freq:%f Hash:%d\n\n", top->GetPattID(),top->GetRowSpan(),top->GetColumnSpan(),top->GetFreq(),top->GetHash());
    Topology::printTop(patt,o);
    o << "\nxCentre: " << top->GetxCOGPix() << " + " << top->GetxCOGshift()
      <<" zCentre: " << top->GetzCOGPix() << " + " << top->GetzCOGshift()
      <<" groupID: " << top->GetGroupID()<< "\n\n...............................................\n";
  }
  o.close();
}

void TopDatabase::Grouping(Int_t NumberofShiftXbins, Int_t NumberofShiftZbins){
  printf("\nGROUPING\n\n");
  Float_t threshold=fThreshold;
  Int_t nPatterns=fN;
  TArrayI shiftXID;
  TArrayI shiftZID;
  Int_t notINgorups=0;
  Float_t shiftXwidth=1./NumberofShiftXbins; //(fraction of the pitch)
  Float_t shiftZwidth=1./NumberofShiftZbins; //(fraction of the pitch)
  Float_t xOffset = shiftXwidth/2;
  Float_t zOffset = shiftZwidth/2;

  shiftXID.Set(nPatterns);
  shiftZID.Set(nPatterns);

  TH1F xShiftDistro("xShiftDistro", "x position", NumberofShiftXbins+1 , -0.5-xOffset, 0.5+xOffset);
  xShiftDistro.SetDirectory(0);

  TH1F zShiftDistro("zShiftDistro", "z position", NumberofShiftZbins+1, -0.5-zOffset, 0.5+zOffset);
  zShiftDistro.SetDirectory(0);
  
  printf("Assigning shift-IDs:\n");
  for(Int_t i=0; i<nPatterns; i++){
    Topology* top = (Topology*)fArrTopologies.At(i);
    if(i%100==0)printf("%d / %d\n ", i, nPatterns);
    xShiftDistro.Fill(top->GetxCOGshift());
    shiftXID[i]=xShiftDistro.FindBin(top->GetxCOGshift());
    zShiftDistro.Fill(top->GetzCOGshift());
    shiftZID[i]=zShiftDistro.FindBin(top->GetzCOGshift());
  }
  Int_t tempgroupID=0;
  printf("Processing patterns over threshold:\n");
  //first assigning groupID to patterns whose frequency is above the treshold 
  for(Int_t t=0; t<nPatterns;t++){
    if(t%100==0)printf("%d / %d\n ", t, nPatterns);
    Topology* top = (Topology*)fArrTopologies.At(t);
    if(top->GetFreq()>=threshold){
      top->SetGroupID(tempgroupID);
      tempgroupID++;
      notINgorups++;
    }
  }
  printf("Processing patterns under threshold:\n");
  for(Int_t i=0; i<nPatterns; i++){
    if(i%100==0)printf("%d / %d\n ", i, nPatterns);
    Topology* top = (Topology*)fArrTopologies.At(i);
    if(top->GetGroupID()!=-1) continue;	
    top->SetGroupID(tempgroupID);
    //if(i%10000==0)printf("Assigning group ID %d / %d\n ", i, nPatterns);
    for(Int_t j=i+1; j<nPatterns; j++){
      Topology* top1 = (Topology*)fArrTopologies.At(j);
      if(top1->GetGroupID()!=-1) continue;
      else {
	if(shiftXID[j]==shiftXID[i] && shiftZID[j]==shiftZID[i] && 
	   top->GetRowSpan()==top1->GetRowSpan() && top->GetColumnSpan()==top1->GetColumnSpan()) top1->SetGroupID(tempgroupID);
      }
    }
    tempgroupID++;
  }
  fNGroups = tempgroupID;
  fOverThr = notINgorups;
  //**********************************************************//
  //                                                          //
  //                 OPERATIONS WITH HISTOGRAMS               //
  //                                                          //
  //**********************************************************//
  fArrHisto.Expand(fNGroups*4);
  fArrHisto.Clear();
  printf("Number of groups: %d\n", fNGroups);
  
  printf("Summing group-histograms:\n");
  for(Int_t iGroup=0; iGroup<fNGroups; iGroup++){
    if(iGroup%1000==0) printf("%d / %d\n", iGroup, fNGroups);
    Bool_t FirstMatch = kTRUE;
    for(Int_t i=0; i<nPatterns; i++){
      // printf("iGroup: %d i: %d\n", iGroup, i);
      Topology* top = (Topology*)fArrTopologies.At(i);
      if(top->GetGroupID()==iGroup){
	if(FirstMatch==kTRUE){
	  FirstMatch==kFALSE;
	  fArrHisto[iGroup*4] = new TH2F(*(top->GetHxA()));
	  fArrHisto[iGroup*4+1] = new TH2F(*(top->GetHzA()));
	  fArrHisto[iGroup*4+2] = new TH2F(*(top->GetHxB()));
	  fArrHisto[iGroup*4+3] = new TH2F(*(top->GetHzB()));
	  top->DeleteHistos();
	  top->SetHxA((TH2F*)fArrHisto[iGroup*4]);
	  top->SetHzA((TH2F*)fArrHisto[iGroup*4+1]);
	  top->SetHxB((TH2F*)fArrHisto[iGroup*4+2]);
	  top->SetHzB((TH2F*)fArrHisto[iGroup*4+3]);
	}
	else{
	  ((TH2F*)fArrHisto[iGroup*4])->Add(top->GetHxA());
	  ((TH2F*)fArrHisto[iGroup*4+1])->Add(top->GetHzA());
	  ((TH2F*)fArrHisto[iGroup*4+2])->Add(top->GetHxB());
	  ((TH2F*)fArrHisto[iGroup*4+3])->Add(top->GetHzB());
	  top->DeleteHistos();
	  top->SetHxA((TH2F*)fArrHisto[iGroup*4]);
	  top->SetHzA((TH2F*)fArrHisto[iGroup*4+1]);
	  top->SetHxB((TH2F*)fArrHisto[iGroup*4+2]);
	  top->SetHzB((TH2F*)fArrHisto[iGroup*4+3]);
	}
      }
    }
  }
  static TF1* gs = new TF1("gs","gaus",-50,50);
  static TF1* gs2 = new TF1("gs2","gaus",-50,50);
  printf("Taking data from histograms:\n");
  for(Int_t iGroup=0; iGroup<fNGroups; iGroup++){
    if(iGroup%1000==0) printf("%d / %d\n", iGroup, fNGroups);
    //X projection
    Int_t fitStatusX=0;
    TH2* hXA = (TH2F*)fArrHisto.At(iGroup*4);
    TH1* hdx = hXA->ProjectionY("hdx");
    gs->SetParameters(hdx->GetMaximum(),hdx->GetMean(),hdx->GetRMS());
    if((hdx->GetEntries())<100) fitStatusX = hdx->Fit("gs","ql");
    else fitStatusX = hdx->Fit("gs","q");
    //Z projection
    Int_t fitStatusZ=0;
    TH2* hZA = (TH2F*)fArrHisto.At(iGroup*4+1);
    TH1* hdz = hZA->ProjectionY("hdz");
    gs2->SetParameters(hdz->GetMaximum(),hdz->GetMean(),hdz->GetRMS());
    if((hdz->GetEntries())<100) fitStatusZ = hdz->Fit("gs2","ql");
    else fitStatusZ = hdz->Fit("gs2","q");
    //*******************************************
    for(Int_t i=0; i<nPatterns; i++){
      Topology* top = (Topology*)fArrTopologies.At(i);
      if(top->GetGroupID()==iGroup){
	//x chunk
	if(fitStatusX==0){
	  top->SetFitStuff(gs->GetParameter(1),Topology::kDeltaXmean);
	  top->SetFitStuff(gs->GetParError(1),Topology::kDeltaXmeanErr);
	  top->SetFitStuff(gs->GetParameter(2),Topology::kDeltaXsigma);
	  top->SetFitStuff(gs->GetParError(2),Topology::kDeltaXsigmaErr);
	  top->SetFitStuff(gs->GetChisquare(),Topology::kChi2x);
	  Int_t varNDFx = gs->GetNDF();
	  if(varNDFx>=0){
	    top->SetFitStuff(varNDFx,Topology::kNDFx);
	  }
	  else{
	    top->SetFitStuff(0.,Topology::kNDFx);
	  }
	}
	else{
	  top->SetFitStuff(0.,Topology::kDeltaXmean);
	  top->SetFitStuff(0.,Topology::kDeltaXmeanErr);
	  top->SetFitStuff(0.,Topology::kDeltaXsigma);
	  top->SetFitStuff(0.,Topology::kDeltaXsigmaErr);
	  top->SetFitStuff(-1,Topology::kChi2x);
	}
	//z chunk
	if(fitStatusZ==0){
	  top->SetFitStuff(gs2->GetParameter(1),Topology::kDeltaZmean);
	  top->SetFitStuff(gs2->GetParError(1),Topology::kDeltaZmeanErr);
	  top->SetFitStuff(gs2->GetParameter(2),Topology::kDeltaZsigma);
	  top->SetFitStuff(gs2->GetParError(2),Topology::kDeltaZsigmaErr);
	  top->SetFitStuff(gs2->GetChisquare(),Topology::kChi2z);
	  Int_t varNDFz = gs2->GetNDF();
	  if(varNDFz>=0){
	    top->SetFitStuff(varNDFz,Topology::kNDFz);
	  }
	  else{
	    top->SetFitStuff(0.,Topology::kNDFz);
	  }
	}
	else{
	  top->SetFitStuff(0.,Topology::kDeltaZmean);
	  top->SetFitStuff(0.,Topology::kDeltaZmeanErr);
	  top->SetFitStuff(0.,Topology::kDeltaZsigma);
	  top->SetFitStuff(0.,Topology::kDeltaZsigmaErr);
	  top->SetFitStuff(-1,Topology::kChi2z);
	}
      }
    }
  }
}

Int_t TopDatabase::FromCluster2GroupID(const AliITSUClusterPix &cl) const{
  
  //Passing from cluster to UChar_t*
  Int_t rs = cl.GetPatternRowSpan();
  Int_t cs = cl.GetPatternColSpan();
  Int_t nBytes = rs*cs/8;
  if((rs*cs)%8 != 0) nBytes++;
  nBytes+=2; //first byte: number of rows; second byte: number of columns;
  const Int_t length = nBytes;
  UChar_t Word[length];
  for(int k=length; k--;) Word[k]=0;
  Word[0]=rs;
  Word[1]=cs;
  UChar_t tempChar=0;
  Int_t index=2;
  Int_t BitCounter=7;
  for(Int_t ir=0; ir<rs; ir++){
    for(Int_t ic=0; ic<cs; ic++){
      if(BitCounter<0) {
	Word[index]=tempChar;
	tempChar=0;
	BitCounter=7;
	index++;
      }	
      if(cl.TestPixel(ir,ic)) tempChar+=(1<<BitCounter);
      BitCounter--;
    }
  }
  Word[index]=tempChar;
  //Creating Hash
  Int_t hashcode = (Int_t)Topology::FuncMurmurHash2(Word,length);
  //Looking up in the database with intepolation search
  Int_t low = 0;
  Int_t high = fN-1;
  Int_t interIndex=-1;
  Int_t min=-2147483648;//minimum Int_t value
  Int_t max=2147483647;
  while(1){
    if(low>high || hashcode<min || hashcode>max){
      interIndex=-1;
      break;
    }
    Int_t guess;
    if(high==low) guess=high;
    else{
      Int_t size = high-low;
      Int_t offset = (Int_t)(((size-1)*((Long_t)hashcode-(Long_t)min))/((Long_t)max-(Long_t)min));
      guess = low+offset;
    }
    if( ( (Topology*)fArrTopologies.At(guess) )->GetHash()==hashcode){
      interIndex=guess;
      break;
    }
    else if( ( (Topology*)fArrTopologies.At(guess) )->GetHash()>hashcode){
      high = guess-1;
      max = ( (Topology*)fArrTopologies.At(guess-1) )->GetHash();
    } 
    else{
      low = guess+1;
      min = ( (Topology*)fArrTopologies.At(guess+1) )->GetHash();
    }  
  }
  if(interIndex==-1){//junk case
    return -1;
  }
  //Solving clashes if necessary
  if(( (Topology*)fArrTopologies.At(interIndex) )->GetFlag()==1){
    //printf("Clash found.\n");
    Int_t Part = 0;
    Int_t BitPosition = 0;
    for(Int_t ir=0; ir<rs;ir++){
      for(Int_t ic=0; ic<cs; ic++){
	if(cl.TestPixel(ir,ic)) Part+=(1<<(31-BitPosition));
	BitPosition++;
	if(BitPosition==32) break;
      }
      if(BitPosition==32) break;
    }
    Bool_t IndexFound = kFALSE;
    Int_t guessIndex = interIndex;
    while( ( (Topology*)fArrTopologies.At(guessIndex) )->GetHash()==hashcode && IndexFound==kFALSE){
      if(Part==( (Topology*)fArrTopologies.At(guessIndex) )->GetPartialTop()){
	IndexFound = kTRUE;
	interIndex = guessIndex;
      }
      guessIndex--;
    }
    guessIndex = interIndex;
    while( ( (Topology*)fArrTopologies.At(guessIndex) )->GetHash()==hashcode && IndexFound==kFALSE){
      if(Part==( (Topology*)fArrTopologies.At(guessIndex) )->GetPartialTop() ){
	IndexFound = kTRUE;
	interIndex = guessIndex;
      }
      guessIndex++;
    }
    if(IndexFound==kFALSE){
      printf("Index not found after clash\n");
      exit(1);
    }
  }
  return ( (Topology*)fArrTopologies.At(interIndex) )->GetGroupID(); 
}
