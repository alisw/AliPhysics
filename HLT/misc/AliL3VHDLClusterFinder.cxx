// $Id$

// Author: Constantin Loizides <mailto:loizides@ikf.physik.uni-frankfurt.de>
//*-- Copyright & copy CL

#include "AliL3StandardIncludes.h"

#include "AliL3VHDLClusterFinder.h"

#if GCCVERSION == 3
using namespace std;
#endif

/** \class AliL3VHDLClusterFinder
//<pre>
//____________________________________________________
// AliL3VHDLClusterFinder
//
// The current VHDL cluster finder for HLT
// Based on STAR L3
//
// Most important parameters:
// fThreshold - threshold for noise clusters
// fMatch - length in time for overlapping sequences
//</pre> 
*/

ClassImp(AliL3VHDLClusterFinder)

AliL3VHDLClusterFinder::AliL3VHDLClusterFinder()
{
  fMatch = 4;
  fThreshold = 10;
  fMinMerge = 1;
  fNClusters=0;

  fXYErr = 0.2;
  fZErr = 0.3;
  fDeconvPad = kTRUE;
  fDeconvTime = kTRUE;
  fstdout = kFALSE;
  fcalcerr = kTRUE;

  Clear();
#ifdef DEBUG
  fdeb=fopen("vhdlclusterfinder.debug","w");
  //fdeb=stderr;
#endif
}

AliL3VHDLClusterFinder::~AliL3VHDLClusterFinder()
{
#ifdef DEBUG
  fclose(fdeb);
#endif
}

void AliL3VHDLClusterFinder::ProcessDigits()
{
  //Loop over data like the analyzer of the VHDL code
  const UChar_t n=255;
  UShort_t rrow=0,rtime=0;
  UChar_t rpad=0,i=n;
  UShort_t *charges=new UShort_t[n];
  Int_t tc=0,mp=0,mt=0,sp=0,st=0;

  fNClusters=0;
  fRow=0;
  fPad=0;
  Clear();

  //loop over input data
  while(fAltromem.ReadSequence(rrow,rpad,rtime,i,&charges)){
    tc=0;mp=0;mt=0;sp=0;st=0;
#if 0
    cout << "Padrow " << (int)rrow << " pad " << (int)rpad << " time " <<(int)rtime << " charges ";
    for(UChar_t ii=0;ii<i;ii++) cout << (int)charges[ii] << " ";
    cout << endl;
#endif
#ifdef DEBUG
    fprintf(fdeb,"ProcessDigits: Input Data: %d %d %d charges:",(int)rrow,(int)rpad,(int)rtime);
    for(UChar_t ii=0;ii<i;ii++) fprintf(fdeb," %d",(int)charges[ii]);
    fprintf(fdeb,"\n");
#endif

    fNRow=rrow;
    fNPad=rpad;

    //calculate sequence values 
    //no deconvulution so far
    for(UChar_t ii=0;ii<i;ii++){
      tc+=charges[ii];
      mt+=(rtime-ii)*charges[ii];
      st+=(rtime-ii)*(rtime-ii)*charges[ii];
    }
    mp=rpad*tc;
    sp=rpad*rpad*tc;

    fSeq.fTotalCharge=tc;
    fSeq.fPad=mp;
    fSeq.fTime=mt;
    fSeq.fPad2=sp;
    fSeq.fTime2=st; 
    fSeq.fMean=0;
    //if(tc!=0) fSeq.fMean=mt/tc;
    if(tc!=0) fSeq.fMean=rtime-i/2;    
    fSeq.fMerge=0;
    fSeq.fRow=rrow;
    fSeq.fLastPad=rpad;
    
    //work on this sequence
    ProcessSequence();
    //output one cluster
    OutputMemory();

#ifdef DEBUG
    fflush(fdeb);
#endif
    i=n; //store size of charges array
  } //loop over data

  //flush everything left
  FlushMemory();
  while(fOP!=fFP) OutputMemory();
}

void AliL3VHDLClusterFinder::ProcessSequence()
{
  if(fNRow!=fRow) FlushMemory();
  else if(fNPad==fPad+1) PrepareMemory();
  else if(fNPad!=fPad) FlushMemory();

  //store new row and pad values
  fRow=fNRow;
  fPad=fNPad;

#ifdef DEBUG
  fprintf(fdeb,"ProcessSequence: Mean=%d TC=%d ",fSeq.fMean,fSeq.fTotalCharge);
#endif

  CompareSeq(); //merge or insert
}

void AliL3VHDLClusterFinder::PrepareMemory()
{
#ifdef DEBUG
  fprintf(fdeb,"PrepareMemory %d %d %d\n",fRP,fEP,fWP);
#endif
  fRP=fEP;
  fEP=fWP;
}

void AliL3VHDLClusterFinder::FlushMemory()
{
#ifdef DEBUG
  fprintf(fdeb,"FlushMemory %d %d %d %d\n",fFP,fRP,fEP,fWP);
#endif
  fFP=fWP;
  fRP=fWP;
  fEP=fWP;
}

void AliL3VHDLClusterFinder::CompareSeq()
{

  while(fRP!=fEP){
    Int_t diff=fSeqs[fPList[fRP]].fMean-fSeq.fMean;

    if(diff>fMatch){ //no match
      IncRPointer(); //cluster finished
      continue;
    } else if(diff<-fMatch){ //no match
      break;                 //insert new cluster       
    }
    else { //match found, merge it
      MergeSeq(); 
      return;
    }
  }

  InsertSeq(); //start new cluster
}

void AliL3VHDLClusterFinder::MergeSeq()
{
#ifdef DEBUG
  fprintf(fdeb,"merged with Mean=%d TC=%d (new Merge=%d) %d %d\n",fSeqs[fPList[fRP]].fMean,fSeqs[fPList[fRP]].fTotalCharge,fSeqs[fPList[fRP]].fMerge+1,fRow,fPad);
#endif
  if(fSeqs[fPList[fRP]].fRow==fSeq.fRow){
    LOG(AliL3Log::kWarning,"AliL3VHDLClusterFinder::","Memory Check")
      <<"Sequences can be merged on the same rows only."<<ENDLOG;
  }
  if(fSeqs[fPList[fRP]].fLastPad+1!=fSeq.fLastPad){
    LOG(AliL3Log::kWarning,"AliL3VHDLClusterFinder::","Memory Check")
      <<"Sequences can be merged on consecutive pads only."<<ENDLOG;
  }

  fSeqs[fPList[fRP]].fMean=fSeq.fMean; //take the new mean
  fSeqs[fPList[fRP]].fLastPad=fSeq.fLastPad;
  fSeqs[fPList[fRP]].fTotalCharge+=fSeq.fTotalCharge;
  fSeqs[fPList[fRP]].fPad+=fSeq.fPad;
  fSeqs[fPList[fRP]].fTime+=fSeq.fTime;
  fSeqs[fPList[fRP]].fPad2+=fSeq.fPad2;
  fSeqs[fPList[fRP]].fTime2+=fSeq.fTime2;
  fSeqs[fPList[fRP]].fMerge+=1;

  fPList[fWP]=fPList[fRP];
  fPList[fRP]=N_clmem; //mark it to be free

  IncRPointer();
  IncWPointer();
}

void AliL3VHDLClusterFinder::InsertSeq()
{
#ifdef DEBUG
  fprintf(fdeb,"inserted %d %d\n",fRow,fPad);
#endif
  NextFreeIndex();    //get next index
  fSeqs[fFirst]=fSeq; //store data
  fPList[fWP]=fFirst; //store index
  IncWPointer();
}

void AliL3VHDLClusterFinder::OutputMemory()
{
  Float_t mtime=0,mpad=0;
  Float_t mtime2=0,mpad2=0;
  UInt_t tc,row,mno;

  if(fOP!=fFP){
  //while(fOP!=fFP){

    UInt_t index=fPList[fOP];
    fPList[fOP]=N_clmem;
    //cout << fOP << " - " << fFP << " - " << index << endl;
    IncPointer(fOP);
    if(index>=N_clmem) return; //nothing to do
    //if(index>=N_clmem) continue; //nothing to do
    
    tc=fSeqs[index].fTotalCharge;
    mno=fSeqs[index].fMerge;
    row=fSeqs[index].fRow;
    if(tc!=0){
      mtime=(Float_t)(fSeqs[index].fTime)/tc;
      mtime2=sqrt((Float_t)(fSeqs[index].fTime2)/tc-mtime*mtime);
    }
    if(tc!=0){
      mpad=(Float_t)(fSeqs[index].fPad)/tc;
      mpad2=sqrt((Float_t)(fSeqs[index].fPad2)/tc-mpad*mpad);
    }

    if(mno<fMinMerge){ //noise cluster
#ifdef DEBUG
    fprintf(fdeb,"OutputMemory %d %d: noise cluster (merge): Row %d Pad %.3f Time %.3f TC %d Merge %d\n",fOP,fFP,row,mpad,mtime,tc,mno);
#endif          
      FreeSeq(index);  
      return;
      //continue;
    }

    if(tc<fThreshold){ //total charge below threshold
#ifdef DEBUG
    fprintf(fdeb,"OutputMemory %d %d: noise cluster (threshold): Row %d Pad %.3f Time %.3f TC %d Merge %d\n",fOP,fFP,row,mpad,mtime,tc,mno);
#endif          
      FreeSeq(index);
      return;
      //continue;
    }

#ifdef DEBUG
    fprintf(fdeb,"OutputMemory %d: Row %d Pad %.3f Time %.3f TC %d Merge %d\n",fNClusters,row,mpad,mtime,tc,mno);
#endif    

    if(stdout)
      cout<<"WriteCluster: padrow "<<row<<" pad "<<mpad<< " +- "<<mpad2<<" time "<<mtime<<" +- "<<mtime2<<" charge "<<tc<<endl;
    
    fNClusters++;
    FreeSeq(index);
  }
}

inline void AliL3VHDLClusterFinder::FreeSeq(UShort_t i)
{
  ClearSeq(i);
  IncPointer(fLast,1,N_clmem);
}

void AliL3VHDLClusterFinder::Clear()
{
  fFirst=0; //first list pointer
  fLast=0;  //last  list pointer

  fWP=0; //write pointer
  fRP=0; //read pointer
  fEP=0; //end read pointer
  fFP=0; //end flush pointer
  fOP=0; //output pointer

  fSeq.fRow=0;
  fSeq.fLastPad=0;
  fSeq.fTotalCharge=0;
  fSeq.fPad=0;
  fSeq.fTime=0;
  fSeq.fPad2=0;
  fSeq.fTime2=0; 
  fSeq.fMean=0;
  fSeq.fMerge=0;

  for(UInt_t i=0;i<N_mem;i++) fPList[i]=N_clmem;
  for(UInt_t i=0;i<N_clmem;i++) ClearSeq(i);
}

void AliL3VHDLClusterFinder::ClearSeq(UShort_t i){
  fSeqs[i].fRow=0;
  fSeqs[i].fLastPad=0;
  fSeqs[i].fTotalCharge=0;
  fSeqs[i].fPad=0;
  fSeqs[i].fTime=0;
  fSeqs[i].fPad2=0;
  fSeqs[i].fTime2=0; 
  fSeqs[i].fMean=0;
  fSeqs[i].fMerge=0;
}

void AliL3VHDLClusterFinder::IncPointer(UShort_t &p, Short_t add, UShort_t N){
  Short_t pp=p;
  pp+=add;  
  if(pp>=N) p=UShort_t(pp-N);
  else if(pp<0) p=UShort_t(pp+N);
  else p=UShort_t(pp);
}

inline void AliL3VHDLClusterFinder::IncRPointer(){
  IncPointer(fRP);
}

inline void AliL3VHDLClusterFinder::IncWPointer(){
  IncPointer(fWP);

  if(fWP==fOP){
    LOG(AliL3Log::kWarning,"AliL3VHDLClusterFinder::IncWPointer","Memory Check")
      <<"Write pointer overwrites output pointer."<<ENDLOG;
  }
}

inline void AliL3VHDLClusterFinder::NextFreeIndex(){
  IncPointer(fFirst,1,N_clmem);
  if(fFirst==fLast) {
    LOG(AliL3Log::kFatal,"AliL3VHDLClusterFinder::GetFreeIndex","Memory Check")
      <<"No space left in sequence list: "<<fFirst<<"=="<<fLast<<ENDLOG;
  }
}

#if 0
void AliL3ClustFinderNew::WriteClusters(Int_t n_clusters,ClusterData *list)
{
  Int_t thisrow,thissector;
  UInt_t counter = fNClusters;
  
  for(int j=0; j<n_clusters; j++)
    {
      if(!list[j].fFlags) continue; //discard 1 pad clusters
      if(list[j].fTotalCharge < fThreshold) continue; //noise cluster

      Float_t xyz[3];      
      Float_t fpad =(Float_t)list[j].fPad /(Float_t)list[j].fTotalCharge;
      Float_t fpad2=fXYErr;
      Float_t ftime =(Float_t)list[j].fTime /(Float_t)list[j].fTotalCharge;
      Float_t ftime2=fZErr;

      if(fcalcerr) {
	fpad2=(Float_t)list[j].fPad2/(Float_t)list[j].fTotalCharge - fpad*fpad;
	fpad2 = sqrt(fpad2);
	ftime2=(Float_t)list[j].fTime2/(Float_t)list[j].fTotalCharge - ftime*ftime;
	ftime2 = sqrt(ftime2); 
      }
       
      if(fstdout==kTRUE)
	cout<<"WriteCluster: padrow "<<fCurrentRow<<" pad "<<fpad << "+-"<<fpad2<<" time "<<ftime<<"+-"<<ftime2<<" charge "<<list[j].fTotalCharge<<endl;

      AliL3Transform::Slice2Sector(fCurrentSlice,fCurrentRow,thissector,thisrow);
      AliL3Transform::Raw2Local(xyz,thissector,thisrow,fpad,ftime);
      if(xyz[0]==0) LOG(AliL3Log::kError,"AliL3ClustFinder","Cluster Finder")
	<<AliL3Log::kDec<<"Zero cluster"<<ENDLOG;
      if(fNClusters >= fMaxNClusters)
	{
	  LOG(AliL3Log::kError,"AliL3ClustFinder::WriteClusters","Cluster Finder")
	    <<AliL3Log::kDec<<"Too many clusters"<<ENDLOG;
	  return;
	}  
      fSpacePointData[counter].fCharge = list[j].fTotalCharge;
      fSpacePointData[counter].fX = xyz[0];
      fSpacePointData[counter].fY = xyz[1];
      fSpacePointData[counter].fZ = xyz[2];
      fSpacePointData[counter].fPadRow = fCurrentRow;
      fSpacePointData[counter].fXYErr = fpad2;
      fSpacePointData[counter].fZErr  = ftime2;
      fSpacePointData[counter].fID = counter
	+((fCurrentSlice&0x7f)<<25)+((fCurrentPatch&0x7)<<22);//Uli

      fNClusters++;
      counter++;
    }
}
#endif

