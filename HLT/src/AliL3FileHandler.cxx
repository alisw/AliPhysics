
//Author:        Uli Frankenfeld
//Last Modified: 17.12.2000

#include <math.h>
#include <iostream.h>
#include <TFile.h>
#include <TTree.h>
#include <stdio.h>

#include "AliL3Transform.h"
#include "AliL3Logging.h"
#include "AliL3MemHandler.h"
#include "AliL3FileHandler.h"

#include "AliTPCClustersArray.h"
#include "AliTPCcluster.h"
#include "AliTPCClustersRow.h"
#include "AliTPCParam.h"
#include "AliSimDigits.h"

#include "AliL3DigitData.h"
#include "AliL3TrackSegmentData.h"
#include "AliL3SpacePointData.h"
#include "AliL3TrackArray.h"
//_____________________________________________________________
//
// The L3 Binary File handler 
//

ClassImp(AliL3FileHandler)

AliL3FileHandler::AliL3FileHandler(){
  //Default constructor
  fInAli = 0;
  fParam = 0;
  fTransformer = 0;
  fMC =0;
}


AliL3FileHandler::~AliL3FileHandler(){
  //Destructor
  if(fTransformer) delete fTransformer;
  if(fMC) CloseMCOutput();
}

Bool_t AliL3FileHandler::SetMCOutput(char *name){
  fMC = fopen(name,"w");
  if(!fMC){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::SetMCOutput","File Open")
    <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliL3FileHandler::SetMCOutput(FILE *file){
  fMC = file;
  if(!fMC){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::SetMCOutput","File Open")
    <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

void AliL3FileHandler::CloseMCOutput(){
  if(!fMC){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::CloseMCOutPut","File Close")
    <<"Nothing to Close"<<ENDLOG;
    return;
  }
  fclose(fMC);
  fMC =0;
}

Bool_t AliL3FileHandler::SetAliInput(){
  if(!fInAli->IsOpen()){
    LOG(AliL3Log::kError,"AliL3FileHandler::SetAliInput","File Open")
    <<"Ali File "<<fInAli->GetName()<<" does not exist"<<ENDLOG;
    return kFALSE;
  }
  fParam = (AliTPCParam*)fInAli->Get("75x40_100x60");
  if(!fParam){ 
    LOG(AliL3Log::kError,"AliL3FileHandler::SetAliInput","File Open")
    <<"No AliTPCParam 75x40_100x60 in File "<<fInAli->GetName()<<ENDLOG;
     return kFALSE;
  }
  fTransformer = new AliL3Transform();
  return kTRUE;
}

Bool_t AliL3FileHandler::SetAliInput(char *name){
  fInAli= new TFile(name,"READ");
  if(!fInAli){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::SetAliInput","File Open")
    <<"Pointer to TFile = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return SetAliInput();
}

Bool_t AliL3FileHandler::SetAliInput(TFile *file){
  fInAli=file;
  if(!fInAli){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::SetAliInput","File Open")
    <<"Pointer to TFile = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return SetAliInput();
}

void AliL3FileHandler::CloseAliInput(){
  if(!fInAli){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::CloseAliInput","File Close")
    <<"Nothing to Close"<<ENDLOG;
     return;
  }
  if(fInAli->IsOpen()) fInAli->Close();
  delete fInAli;
  fInAli = 0;
}

Bool_t AliL3FileHandler::IsDigit(){
  if(!fInAli){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::IsDigit","File")
    <<"Pointer to TFile = 0x0 "<<ENDLOG;
    return kTRUE;  //may you are use binary input which is Digits!!
  }
  TTree *t=(TTree*)fInAli->Get("TreeD_75x40_100x60");
  if(t){
    LOG(AliL3Log::kInformational,"AliL3FileHandler::IsDigit","File Type")
    <<"Found Digit Tree -> Use Fast Cluster Finder"<<ENDLOG;
    return kTRUE;
  }
  else{
    LOG(AliL3Log::kInformational,"AliL3FileHandler::IsDigit","File Type")
    <<"No Digit Tree -> Use Cluster Tree"<<ENDLOG;
    return kFALSE;
  }
}

///////////////////////////////////////// Digit IO  
Bool_t AliL3FileHandler::AliDigits2Binary(){
  Bool_t out = kTRUE;
  UInt_t nrow;
  AliL3DigitRowData* data = AliDigits2Memory(nrow);
  out = Memory2Binary(nrow,data);
  Free();
  return out;
}


Bool_t AliL3FileHandler::AliDigits2CompBinary(){
  Bool_t out = kTRUE;
  UInt_t ndigits=0;
  AliL3DigitRowData* digits = AliDigits2Memory(ndigits);
  out = Memory2CompBinary(ndigits,digits);
  Free();
  return out;
}


AliL3DigitRowData * AliL3FileHandler::AliDigits2Memory(UInt_t & nrow){
  AliL3DigitRowData *data = 0;
  nrow=0;
  if(!fInAli){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::AliDigits2Memory","File")
    <<"No Input avalible: no object TFile"<<ENDLOG;
    return 0; 
  }
  if(!fInAli->IsOpen()){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::AliDigits2Memory","File")
    <<"No Input avalible: TFile not opend"<<ENDLOG;
    return 0;
  }

  TDirectory *savedir = gDirectory;
  fInAli->cd();
  TTree *t=(TTree*)fInAli->Get("TreeD_75x40_100x60");
  if(!t){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::AliDigits2Binary","AliRoot")
    <<"No Digit Tree inside!"<<ENDLOG;
    return 0;
  }
  AliSimDigits digarr, *dummy=&digarr;
  t->GetBranch("Segment")->SetAddress(&dummy);
  UShort_t dig;
  Int_t time,pad,sector,row;
  Int_t nrows=0;
  Int_t ndigitcount=0;
  Int_t entries = (Int_t)t->GetEntries();
  Int_t ndigits[entries];
  Int_t lslice,lrow;
  for(Int_t n=0; n<t->GetEntries(); n++)
    {
      t->GetEvent(n);
      fParam->AdjustSectorRow(digarr.GetID(),sector,row);
      fTransformer->Sector2Slice(lslice,lrow,sector,row);
      if(fSlice != lslice || lrow<fRowMin || lrow>fRowMax) continue;
//      if(fSlice != lslice) continue;

      Float_t xyz[3];
      ndigits[lrow] = 0;
      digarr.First();
      do {
        time=digarr.CurrentRow();
        pad=digarr.CurrentColumn();
        dig = digarr.GetDigit(time,pad);
        if(dig<=fParam->GetZeroSup()) continue;
        if(time < fParam->GetMaxTBin()-1 && time > 0)
          if(digarr.GetDigit(time+1,pad) <= fParam->GetZeroSup()
             && digarr.GetDigit(time-1,pad) <= fParam->GetZeroSup())
            continue;

        fTransformer->Raw2Local(xyz,sector,row,pad,time);
        if(fParam->GetPadRowRadii(sector,row)<230./250.*fabs(xyz[2]))
          continue;

        ndigits[lrow]++; //for this row only
        ndigitcount++;  //total number of digits to be published

      } while (digarr.Next());

      nrows++;
    }
  Int_t size = sizeof(AliL3DigitData)*ndigitcount
                                      + nrows*sizeof(AliL3DigitRowData);

  LOG(AliL3Log::kDebug,"AliL3FileHandler::AliDigits2Memory","Digits")
  <<AliL3Log::kDec<<"Found "<<ndigitcount<<" Digits"<<ENDLOG;

  data=(AliL3DigitRowData*) Allocate(size);
  nrow = (UInt_t)nrows;
  AliL3DigitRowData *tempPt = data;
  for(Int_t n=0; n<t->GetEntries(); n++)
    {
      t->GetEvent(n);

      Float_t xyz[3];
      fParam->AdjustSectorRow(digarr.GetID(),sector,row);
      fTransformer->Sector2Slice(lslice,lrow,sector,row);
//      if(fSlice != lslice) continue;
      if(fSlice != lslice || lrow<fRowMin || lrow>fRowMax) continue;
      tempPt->fRow = lrow;
      tempPt->fNDigit = ndigits[lrow];

      Int_t localcount=0;
      digarr.First();
      do {
        dig=digarr.CurrentDigit();
        if (dig<=fParam->GetZeroSup()) continue;
        time=digarr.CurrentRow();
        pad=digarr.CurrentColumn();
        if(time < fParam->GetMaxTBin()-1 && time > 0)
          if(digarr.GetDigit(time-1,pad) <= fParam->GetZeroSup() &&
             digarr.GetDigit(time+1,pad) <= fParam->GetZeroSup()) continue;

        //Exclude data outside cone:
        fTransformer->Raw2Local(xyz,sector,row,pad,time);
        if(fParam->GetPadRowRadii(sector,row)<230./250.*fabs(xyz[2]))
          continue;

        if(localcount >= ndigits[lrow])
          LOG(AliL3Log::kFatal,"AliL3FileHandler::AliDigits2Binary","Memory")
          <<AliL3Log::kDec<<"Mismatch: localcount "<<localcount<<" ndigits "
          <<ndigits[lrow]<<ENDLOG;

        tempPt->fDigitData[localcount].fCharge=dig;
        tempPt->fDigitData[localcount].fPad=pad;
        tempPt->fDigitData[localcount].fTime=time;
        localcount++;
      } while (digarr.Next());

      Byte_t *tmp = (Byte_t*)tempPt;
      Int_t size = sizeof(AliL3DigitRowData)
                                      + ndigits[lrow]*sizeof(AliL3DigitData);
      tmp += size;
      tempPt = (AliL3DigitRowData*)tmp;
    }
  savedir->cd(); 
  return data;
}

///////////////////////////////////////// Point IO  
Bool_t AliL3FileHandler::AliPoints2Binary(){
  Bool_t out = kTRUE;
  UInt_t npoint;
  AliL3SpacePointData *data = AliPoints2Memory(npoint);
  out = Memory2Binary(npoint,data);
  Free();
  return out;
}

AliL3SpacePointData * AliL3FileHandler::AliPoints2Memory(UInt_t & npoint){
  AliL3SpacePointData *data = 0;
  npoint=0;
  if(!fInAli){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::AliPoints2Memory","File")
    <<"No Input avalible: no object TFile"<<ENDLOG;
    return 0;
  }
  if(!fInAli->IsOpen()){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::AliPoints2Memory","File")
    <<"No Input avalible: TFile not opend"<<ENDLOG;
    return 0;
  }
  TDirectory *savedir = gDirectory;
  fInAli->cd();

  AliTPCClustersArray carray;
  carray.Setup(fParam);
  carray.SetClusterType("AliTPCcluster");
  Bool_t clusterok = carray.ConnectTree("Segment Tree");
  if(!clusterok) return 0;

  AliTPCClustersRow ** clusterrow = 
               new AliTPCClustersRow*[ (int)carray.GetTree()->GetEntries()];
  Int_t *rows = new int[ (int)carray.GetTree()->GetEntries()];
  Int_t *sects = new int[  (int)carray.GetTree()->GetEntries()];
  Int_t sum=0;

  Int_t lslice,lrow;
  for(Int_t i=0; i<carray.GetTree()->GetEntries(); i++){
    AliSegmentID *s = carray.LoadEntry(i);
    Int_t sector,row;
    fParam->AdjustSectorRow(s->GetID(),sector,row);
    rows[i] = row;
    sects[i] = sector;
    clusterrow[i] = 0;
    fTransformer->Sector2Slice(lslice,lrow,sector,row);
    if(fSlice != lslice || lrow<fRowMin || lrow>fRowMax) continue;
    clusterrow[i] = carray.GetRow(sector,row);
    if(clusterrow[i])
      sum+=clusterrow[i]->GetArray()->GetEntriesFast();
  }
  UInt_t size = sum*sizeof(AliL3SpacePointData);

  LOG(AliL3Log::kDebug,"AliL3FileHandler::AliPoints2Memory","File")
  <<AliL3Log::kDec<<"Found "<<sum<<" SpacePoints"<<ENDLOG;

  data = (AliL3SpacePointData *) Allocate(size);
  npoint = sum;
  UInt_t n=0; 
  for(Int_t i=0; i<carray.GetTree()->GetEntries(); i++){
    if(!clusterrow[i]) continue;
    Int_t row = rows[i];
    Int_t sector = sects[i];
    fTransformer->Sector2Slice(lslice,lrow,sector,row);
    Int_t entries_in_row = clusterrow[i]->GetArray()->GetEntriesFast();
    for(Int_t j = 0;j<entries_in_row;j++){
      AliTPCcluster *c = (AliTPCcluster*)(*clusterrow[i])[j];
      data[n].fZ = c->GetZ();
      data[n].fY = c->GetY();
      data[n].fX = fParam->GetPadRowRadii(sector,row);
      data[n].fID = n+((fSlice&0x7f)<<25)+((fPatch&0x7)<<22);//uli
      data[n].fPadRow = lrow;
      data[n].fXYErr = c->GetSigmaY2();
      data[n].fZErr = c->GetSigmaZ2();
      if(fMC) fprintf(fMC,"%d %d\n",data[n].fID,c->GetLabel(0));
      n++;
    }
  }
  for(Int_t i=0;i<carray.GetTree()->GetEntries();i++){
    Int_t row = rows[i];
    Int_t sector = sects[i];
    if(carray.GetRow(sector,row))
      carray.ClearRow(sector,row);
  }

  delete [] clusterrow;
  delete [] rows;
  delete [] sects;
  savedir->cd();   

  return data;
}

