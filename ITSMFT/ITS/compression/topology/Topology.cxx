#include <Riostream.h>
#include "TObject.h"
#include "TMath.h"
#include "TBits.h"
#include "AliITSUClusterPix.h"
#include "TH2F.h"
#include "TH1F.h"
#include "Topology.h"
#include "AliLog.h"
#include "TString.h"

ClassImp(Topology)

Int_t Topology::fMode = Topology::kHashes;

Topology::Topology():TObject(), fPattern(), fWordLength(0), fWord(0x0), fUID(0), fRs(0), fCs(0), fFiredPixels(0),
		     fxCOGPix(0.), fzCOGPix(0.), fxCOGshift(0.), fzCOGshift(0.), fHash(0), fFreq(0.), fCounts(0),
		     fGroupID(-1), fHxA(0), fHxB(0), fHzA(0), fHzB(0), fFlag(0), fPartialTop(0), fPattID(-1){
  for(Int_t i=0; i<kFitLength; i++) fArrFit[i]=0;
}

Topology::Topology(const AliITSMFTClusterPix &cluster):TObject(){//UniqueID of the argument must already have been set
  fPattern.Clear();
  fRs = cluster.GetPatternRowSpan();
  fCs = cluster.GetPatternColSpan();
  fUID = fRs<<16 + fCs;
  fPattern.SetUniqueID(fUID);
  
  //Defining the topology
  Int_t tempxCOG = 0;
  Int_t tempzCOG = 0;
  Int_t tempFiredPixels = 0;
  for(Int_t ir=0;ir<fRs;ir++){
    for(Int_t ic=0;ic<fCs;ic++){
      if(cluster.TestPixel(ir,ic)){
	fPattern.SetBitNumber(ir*fCs+ic);
	tempFiredPixels++;
	tempxCOG+=ir;
	tempzCOG+=ic;
      }
    }
  }
  Float_t xsh=Float_t((tempxCOG%tempFiredPixels))/tempFiredPixels; //distance between COG end centre of the pixel containing COG
  Float_t zsh=Float_t((tempzCOG%tempFiredPixels))/tempFiredPixels;
  tempxCOG/=tempFiredPixels;
  tempzCOG/=tempFiredPixels;
  if(xsh>0.5){
    tempxCOG+=1;
    xsh-=1;
  }
  if(zsh>0.5){
    tempzCOG+=1;
    zsh-=1;
  }
  fxCOGPix = (Float_t) tempxCOG+0.5;
  fzCOGPix = (Float_t) tempzCOG+0.5;
  fxCOGshift = xsh;
  fzCOGshift = zsh;
  fFiredPixels = tempFiredPixels;

  //Creating hash
  Int_t nBytes = fRs*fCs/8;
  if((fRs*fCs)%8 != 0) nBytes++;
  nBytes+=2; //first byte: number of rows; second byte: number of columns;
  fWordLength = nBytes;
  fWord = new UChar_t[fWordLength];
  Top2Word(fPattern,fWord);
  fHash=(Int_t)FuncMurmurHash2(fWord, fWordLength);
  
  fFreq=0.; //WARNING: it is to set in a second time
  fCounts=0; //WARNING: it is to set in a second time
  fGroupID=-1; //WARNING: it is to set in a second time
  
  //Creating histograms

  fHxA = new TH2F("hXA","#DeltaX vs #alpha",10,0,TMath::Pi()/2,50,-30,30);
  fHxA->SetDirectory(0);
  fHxA->GetXaxis()->SetTitle("#alpha");
  fHxA->GetYaxis()->SetTitle("#DeltaX (#mum)");

  fHxB = new TH2F("hXB","#DeltaX vs #beta",10,0,TMath::Pi()/2,50,-30,30);
  fHxB->SetDirectory(0);
  fHxB->GetXaxis()->SetTitle("#alpha");
  fHxB->GetYaxis()->SetTitle("#DeltaX (#mum)");
  
  fHzA = new TH2F("hZA","#DeltaZ vs #alpha",10,0,TMath::Pi()/2,50,-30,30);
  fHzA->SetDirectory(0);
  fHzA->GetXaxis()->SetTitle("#alpha");
  fHzA->GetYaxis()->SetTitle("#DeltaZ (#mum)");

  fHzB = new TH2F("hZB","#DeltaZ vs beta",10,0,TMath::Pi()/2,50,-30,30);
  fHzB->SetDirectory(0);
  fHzB->GetXaxis()->SetTitle("#beta");
  fHzB->GetYaxis()->SetTitle("#DeltaZ (#mum)");

  for(Int_t i=0; i<kFitLength; i++) fArrFit[i]=0;
  fFlag=0;
  fPartialTop = Top2Int(fPattern);
  fPattID = -1;
}



Topology::Topology(const TBits &top2copy):TObject(),fPattern(top2copy){//UniqueID of the argument must already have been set
  fUID = top2copy.GetUniqueID();
  fRs = fUID>>16;
  fCs = fUID&0xffff;
  Int_t nBytes = fRs*fCs/8;
  if((fRs*fCs)%8 != 0) nBytes++;
  nBytes+=2; //first byte: number of rows; second byte: number of columns;
  fWordLength = nBytes;
  fWord = new UChar_t[fWordLength];
  Top2Word(fPattern,fWord);

  //Defining COG position
  Int_t tempxCOG = 0;
  Int_t tempzCOG = 0;
  Int_t tempFiredPixels = 0;
  for(Int_t ir=0;ir<fRs;ir++){
    for(Int_t ic=0;ic<fCs;ic++){
      if(fPattern.TestBitNumber(ir*fCs+ic)){
	tempFiredPixels++;
	tempxCOG+=ir;
	tempzCOG+=ic;
      }
    }
  }
  Float_t xsh=Float_t((tempxCOG%tempFiredPixels))/tempFiredPixels; //distance between COG end centre of the pixel containing COG
  Float_t zsh=Float_t((tempzCOG%tempFiredPixels))/tempFiredPixels;
  tempxCOG/=tempFiredPixels;
  tempzCOG/=tempFiredPixels;
  if(xsh>0.5){
    tempxCOG+=1;
    xsh-=1;
  }
  if(zsh>0.5){
    tempzCOG+=1;
    zsh-=1;
  }
  fxCOGPix = (Float_t) tempxCOG+0.5;
  fzCOGPix = (Float_t) tempzCOG+0.5;
  fxCOGshift = xsh;
  fzCOGshift = zsh;
  fFiredPixels = tempFiredPixels;
  
  fFreq=0.; //WARNING: it is to set in a second time
  fCounts=0; //WARNING: it is to set in a second time
  fGroupID=-1; //WARNING: it is to set in a second time
  fHash=(Int_t)FuncMurmurHash2(fWord, fWordLength);
  
  //Creating histograms

  fHxA = new TH2F("hXA","#DeltaX vs #alpha",10,0,TMath::Pi()/2,50,-30,30);
  fHxA->SetDirectory(0);
  fHxA->GetXaxis()->SetTitle("#alpha");
  fHxA->GetYaxis()->SetTitle("#DeltaX (#mum)");

  fHxB = new TH2F("hXB","#DeltaX vs #beta",10,0,TMath::Pi()/2,50,-30,30);
  fHxB->SetDirectory(0);
  fHxB->GetXaxis()->SetTitle("#alpha");
  fHxB->GetYaxis()->SetTitle("#DeltaX (#mum)");
  
  fHzA = new TH2F("hXA","#DeltaZ vs #alpha",10,0,TMath::Pi()/2,50,-30,30);
  fHzA->SetDirectory(0);
  fHzA->GetXaxis()->SetTitle("#alpha");
  fHzA->GetYaxis()->SetTitle("#DeltaZ (#mum)");

  fHzB = new TH2F("hXA","#DeltaZ vs beta",10,0,TMath::Pi()/2,50,-30,30);
  fHzB->SetDirectory(0);
  fHzB->GetXaxis()->SetTitle("#beta");
  fHzB->GetYaxis()->SetTitle("#DeltaZ (#mum)");

  for(Int_t i=0; i<kFitLength; i++) fArrFit[i]=0;

  fFlag=0;
  fPartialTop = Top2Int(fPattern);
  fPattID=-1;
}

Topology::Topology(const Topology &topo):TObject(),fPattern(topo.GetPattern()){
  fUID = topo.GetUniqueID();
  fRs = fUID>>16;
  fCs = fUID&0xffff;
  Int_t nBytes = fRs*fCs/8;
  if((fRs*fCs)%8 != 0) nBytes++;
  nBytes+=2; //first byte: number of rows; second byte: number of columns;
  fWordLength = nBytes;
  fWord = new UChar_t[fWordLength];
  Top2Word(fPattern,fWord);
  fFreq = topo.GetFreq();
  fCounts = topo.GetCounts();
  fHash = topo.GetHash();
  fGroupID = topo.GetGroupID();
  fFiredPixels = topo.GetFiredPixels();
  fxCOGPix = topo.GetxCOGPix();
  fzCOGPix = topo.GetzCOGPix();
  fxCOGshift = topo.GetxCOGshift();
  fzCOGshift = topo.GetzCOGshift();
  fHxA = new TH2F(*topo.GetHxA());
  fHxB = new TH2F(*topo.GetHxB());
  fHzA = new TH2F(*topo.GetHzA());
  fHzB = new TH2F(*topo.GetHzB());
  for(Int_t i=0; i<kFitLength; i++) fArrFit[i]=topo.GetFitStuff(i);
  fFlag=0;
  fPartialTop = Top2Int(fPattern);
  fPattID=-1;
}

Topology::~Topology(){
  delete[] fWord;
  delete fHxA;
  delete fHxB;
  delete fHzA;
  delete fHzB;
}

void Topology::Top2Word(const TBits &top, UChar_t* Word){
  UInt_t UID = top.GetUniqueID();
  Int_t rs = UID>>16;
  Int_t cs = UID&0xffff;
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
      if(top.TestBitNumber(ir*cs+ic)) tempChar+=(1<<BitCounter);
      BitCounter--;
    }
  }
  Word[index]=tempChar;
}

void Topology::Word2Top(const UChar_t* Word, TBits &top){
  Int_t rs = Word[0];
  Int_t cs = Word[1];
  Int_t nBytes = cs*rs/8;
  UChar_t tempChar=0;
  Int_t indBit=0;
  Int_t s=0;
  if((cs*rs)%8!=0) nBytes++;
  for(int i=2; i<nBytes+2; i++){
    tempChar=Word[i];
    s=128; //0b10000000
    while(s>0){
      if((tempChar&s)!=0) top.SetBitNumber(indBit);
      s/=2;
      indBit++;
    }
  }
  top.SetUniqueID((rs<<16)+cs);
}

std::ostream& Topology::printTop(TBits top, std::ostream &out){
  UInt_t UID = top.GetUniqueID();
  Int_t rs = UID>>16; 
  Int_t cs = UID&0xffff;
  for (Int_t ir=0;ir<rs;ir++){
    out << "|"; 
    for (Int_t ic=0; ic<cs; ic++) {
      out << Form("%c",top.TestBitNumber(ir*cs+ic) ? '+':' ');
    }
    out << ("|\n");
  }
  out<< endl;
}

UInt_t Topology::FuncMurmurHash2(const void* key, Int_t len){
  // 'm' and 'r' are mixing constants generated offline.
  const UInt_t m =0x5bd1e995;
  const Int_t r = 24;
  // Initialize the hash
  UInt_t h = 0;
  // Mix 4 bytes at a time into the hash
  const UChar_t* data = (const UChar_t *)key;
  //Int_t recIndex=0;
  while(len >= 4){
    UInt_t k = *(UInt_t*)data;
    k *= m;
    k ^= k >> r;
    k *= m;

    h *= m;
    h ^= k;
    data += 4;
    len -= 4;
  }
  // Handle the last few bytes of the input array
  switch(len){
  case 3: h ^= data[2] << 16;
  case 2: h ^= data[1] << 8;
  case 1: h ^= data[0];
    h *= m;
  };
  // Do a few final mixes of the hash to ensure the last few
  // bytes are well-incorporated.
  h ^= h >> 13;
  h *= m;
  h ^= h >> 15;
  return h;
}

Int_t Topology::Compare(const TObject* obj) const{
  //Since Sort method of TObjArray sorts object in ascending order ( if +1 means higher and -1 lower),
  //it is necessary to invert the frequency outputs 
  const Topology* top = (const Topology*)obj;
  if(fMode==kFrequency){
    if(fFreq < top->GetFreq()) return +1;
    else if(fFreq == top->GetFreq()) return 0;
    else return -1;
  }
  if(fMode==kHashes){
    if(fHash < top->GetHash()) return -1;
    else if(fHash == top->GetHash()) return 0;
    else return +1;    
  }
  AliFatal(Form("Unknown mode for sorting: %d",fMode));
  return 0;
}

Bool_t Topology::IsEqual(const TObject* obj) const{
  const Topology* top = (const Topology*)obj;
  if(fMode==kFrequency){
    if(fFreq == top->GetFreq()) return kTRUE;
    else return kFALSE;
  }
  if(fMode==kHashes){
    if(fHash == top->GetHash()) return kTRUE;
    else return kFALSE;    
  }
  AliFatal(Form("Unknown mode for sorting: %d",fMode));
  return kFALSE;
}

Int_t Topology::Top2Int(const TBits &top){
  Int_t output=0;
  for(Int_t i=0; i<32; i++){
    if(top.TestBitNumber(i)) output+=(1<<(31-i));
  }
  return output;
}

void Topology::DeleteHistos(){
  delete fHxA;
  delete fHxB;
  delete fHzA;
  delete fHzB;  
}
