//$Id$

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>, Anders Vestbo <mailto:vestbo$fi.uib.no>
//*-- Copyright &copy Uli 

#include "AliL3StandardIncludes.h"
#include <AliTPCDigitsArray.h>
#include <AliTPCClustersArray.h>
#include <AliTPCcluster.h>
#include <AliTPCClustersRow.h>

#include "AliL3Logging.h"
#include "AliL3Transform.h"
#include "AliL3MemHandler.h"
#include "AliL3FileHandler.h"

#include "AliL3DigitData.h"
#include "AliL3TrackSegmentData.h"
#include "AliL3SpacePointData.h"
#include "AliL3TrackArray.h"

//_____________________________________________________________
// AliL3FileHandler
//
// The HLT ROOT <-> binary files handling class
//
// This class provides the interface between AliROOT files,
// and HLT binary files. It should be used for converting 
// TPC data stored in AliROOT format (outputfile from a simulation),
// into the data format currently used by in the HLT framework. 
// This enables the possibility to always use the same data format, 
// whether you are using a binary file as an input, or a AliROOT file.
//
// For example on how to create binary files from a AliROOT simulation,
// see example macro exa/Binary.C.
//
// For reading a AliROOT file into HLT format in memory, do the following:
//
// AliL3FileHandler file;
// file.SetAliInput("galice.root");
// AliL3DigitRowData *dataPt = (AliL3DigitRowData*)file.AliDigits2Memory(nrows,eventnr);
// 
// All the data are then stored in memory and accessible via the pointer dataPt.
// Accesing the data is then identical to the example 1) showed in AliL3MemHandler class.
//
// For converting the data back, and writing it to a new AliROOT file do:
//
// AliL3FileHandler file;
// file.SetAliInput("galice.root");
// file.Init(slice,patch,NumberOfRowsInPatch);
// file.AliDigits2RootFile(dataPt,"new_galice.root");
// file.CloseAliInput();

ClassImp(AliL3FileHandler)

AliL3FileHandler::AliL3FileHandler()
{
  //Default constructor
  fInAli = 0;
  fParam = 0;
  fMC =0;
  fLastIndex=0;
  fDigits=0;
  fDigitsTree=0;
}

AliL3FileHandler::~AliL3FileHandler()
{
  //Destructor
  if(fMC) CloseMCOutput();
  FreeDigitsTree();
  if(fInAli) CloseAliInput();
  
}

void AliL3FileHandler::FreeDigitsTree()
{
  if(!fDigitsTree)
    {
      LOG(AliL3Log::kWarning,"AliL3FileHandler::FreeDigitsTree()","Pointer")
	<<"Cannot free digitstree, it is not present"<<ENDLOG;
      return;
    }
  fDigits=0;
  fDigitsTree->Delete();
  fDigitsTree=0;
}

Bool_t AliL3FileHandler::SetMCOutput(char *name)
{
  fMC = fopen(name,"w");
  if(!fMC){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::SetMCOutput","File Open")
      <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliL3FileHandler::SetMCOutput(FILE *file)
{
  fMC = file;
  if(!fMC){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::SetMCOutput","File Open")
      <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

void AliL3FileHandler::CloseMCOutput()
{
  if(!fMC){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::CloseMCOutPut","File Close")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }
  fclose(fMC);
  fMC =0;
}

Bool_t AliL3FileHandler::SetAliInput()
{
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
  return kTRUE;
}

Bool_t AliL3FileHandler::SetAliInput(char *name)
{
  //Open the AliROOT file with name.
  
  fInAli= new TFile(name,"READ");
  if(!fInAli){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::SetAliInput","File Open")
    <<"Pointer to TFile = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return SetAliInput();
}

Bool_t AliL3FileHandler::SetAliInput(TFile *file)
{
  //Specify already opened AliROOT file to use as an input.
  
  fInAli=file;
  if(!fInAli){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::SetAliInput","File Open")
    <<"Pointer to TFile = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return SetAliInput();
}

void AliL3FileHandler::CloseAliInput()
{
  if(!fInAli){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::CloseAliInput","File Close")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }
  if(fInAli->IsOpen()) fInAli->Close();
  delete fInAli;
  fInAli = 0;
  
}

Bool_t AliL3FileHandler::IsDigit()
{
  //Check if there is a TPC digit tree in the current file.
  //Return kTRUE if tree was found, and kFALSE if not found.
  
  if(!fInAli){
    LOG(AliL3Log::kWarning,"AliL3FileHandler::IsDigit","File")
    <<"Pointer to TFile = 0x0 "<<ENDLOG;
    return kTRUE;  //may you are use binary input which is Digits!!
  }
  TTree *t=(TTree*)fInAli->Get("TreeD_75x40_100x60_0");
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
Bool_t AliL3FileHandler::AliDigits2Binary(Int_t event)
{
  
  Bool_t out = kTRUE;
  UInt_t nrow;
  AliL3DigitRowData* data = AliDigits2Memory(nrow,event);
  out = Memory2Binary(nrow,data);
  Free();
  return out;
}

Bool_t AliL3FileHandler::AliDigits2CompBinary(Int_t event)
{
  //Convert AliROOT TPC data, into HLT data format.
  //event specifies the event you want in the aliroot file.
  
  Bool_t out = kTRUE;
  UInt_t ndigits=0;
  AliL3DigitRowData* digits = AliDigits2Memory(ndigits,event);
  out = Memory2CompBinary(ndigits,digits);
  Free();
  return out;
}

AliL3DigitRowData * AliL3FileHandler::AliDigits2Memory(UInt_t & nrow,Int_t event)
{
  //Read data from AliROOT file into memory, and store it in the HLT data format.
  //Returns a pointer to the data.

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

  if(!fDigitsTree)
    GetDigitsTree(event);
  
  UShort_t dig;
  Int_t time,pad,sector,row;
  Int_t nrows=0;
  Int_t ndigitcount=0;
  Int_t entries = (Int_t)fDigitsTree->GetEntries();
  Int_t ndigits[entries];
  Int_t lslice,lrow;
  Float_t xyz[3];
  
  for(Int_t n=fLastIndex; n<fDigitsTree->GetEntries(); n++)
    {
      fDigitsTree->GetEvent(n);
      fParam->AdjustSectorRow(fDigits->GetID(),sector,row);
      AliL3Transform::Sector2Slice(lslice,lrow,sector,row);
      //if(fSlice != lslice || lrow<fRowMin || lrow>fRowMax) continue;
      if(lslice < fSlice) continue;
      if(lslice != fSlice) break;
      if(lrow < fRowMin) continue;
      if(lrow > fRowMax) break;

      ndigits[lrow] = 0;
      fDigits->First();
      do {
        time=fDigits->CurrentRow();
        pad=fDigits->CurrentColumn();
        dig = fDigits->GetDigit(time,pad);
        if(dig<=fParam->GetZeroSup()) continue;
        
	/*
	  if(time < fParam->GetMaxTBin()-1 && time > 0)
	  if(fDigits->GetDigit(time+1,pad) <= fParam->GetZeroSup()
	  && fDigits->GetDigit(time-1,pad) <= fParam->GetZeroSup())
	  continue;
	*/
	
        AliL3Transform::Raw2Local(xyz,sector,row,pad,time);
        if(fParam->GetPadRowRadii(sector,row)<230./250.*fabs(xyz[2]))
          continue;

        ndigits[lrow]++; //for this row only
        ndigitcount++;   //total number of digits to be published

      } while (fDigits->Next());
      //cout << lrow << " " << ndigits[lrow] << " - " << ndigitcount << endl;
      nrows++;
    }

  Int_t size = sizeof(AliL3DigitData)*ndigitcount
    + nrows*sizeof(AliL3DigitRowData);

  LOG(AliL3Log::kDebug,"AliL3FileHandler::AliDigits2Memory","Digits")
    <<AliL3Log::kDec<<"Found "<<ndigitcount<<" Digits"<<ENDLOG;
  
  data=(AliL3DigitRowData*) Allocate(size);
  nrow = (UInt_t)nrows;
  AliL3DigitRowData *tempPt = data;
  for(Int_t n=fLastIndex; n<fDigitsTree->GetEntries(); n++)
    {
      fDigitsTree->GetEvent(n);
      fParam->AdjustSectorRow(fDigits->GetID(),sector,row);
      AliL3Transform::Sector2Slice(lslice,lrow,sector,row);
      //if(fSlice != lslice || lrow<fRowMin || lrow>fRowMax) continue;
      if(lslice < fSlice) continue;
      if(lslice != fSlice) break;
      if(lrow < fRowMin) continue;
      if(lrow > fRowMax) break;

      tempPt->fRow = lrow;
      tempPt->fNDigit = ndigits[lrow];

      Int_t localcount=0;
      fDigits->First();
      do {
        //dig=fDigits->CurrentDigit();
	//if (dig<=fParam->GetZeroSup()) continue;
        time=fDigits->CurrentRow();
        pad=fDigits->CurrentColumn();
        dig = fDigits->GetDigit(time,pad);
	if (dig <= fParam->GetZeroSup()) continue;
	
	/*
	  if(time < fParam->GetMaxTBin()-1 && time > 0)
          if(fDigits->GetDigit(time-1,pad) <= fParam->GetZeroSup() &&
	  fDigits->GetDigit(time+1,pad) <= fParam->GetZeroSup()) continue;
	*/
        
	//Exclude data outside cone:
        AliL3Transform::Raw2Local(xyz,sector,row,pad,time);
        if(fParam->GetPadRowRadii(sector,row)<230./250.*fabs(xyz[2]))
          continue;

        if(localcount >= ndigits[lrow])
          LOG(AliL3Log::kFatal,"AliL3FileHandler::AliDigits2Binary","Memory")
	    <<AliL3Log::kDec<<"Mismatch: localcount "<<localcount<<" ndigits "
	    <<ndigits[lrow]<<ENDLOG;
	
        tempPt->fDigitData[localcount].fCharge=dig;
        tempPt->fDigitData[localcount].fPad=pad;
        tempPt->fDigitData[localcount].fTime=time;
#ifdef do_mc
	tempPt->fDigitData[localcount].fTrackID[0] = fDigits->GetTrackID(time,pad,0);
	tempPt->fDigitData[localcount].fTrackID[1] = fDigits->GetTrackID(time,pad,1);
	tempPt->fDigitData[localcount].fTrackID[2] = fDigits->GetTrackID(time,pad,2);
#endif
        localcount++;
      } while (fDigits->Next());

      Byte_t *tmp = (Byte_t*)tempPt;
      Int_t size = sizeof(AliL3DigitRowData)
                                      + ndigits[lrow]*sizeof(AliL3DigitData);
      tmp += size;
      tempPt = (AliL3DigitRowData*)tmp;
      //fLastIndex=n;
    }
  //fLastIndex++;
  return data;
}

Bool_t AliL3FileHandler::GetDigitsTree(Int_t event)
{
  //Connects to the TPC digit tree in the AliROOT file.
  
  fInAli->cd();
  Char_t dname[100];
  sprintf(dname,"TreeD_75x40_100x60_%d",event);
  fDigitsTree = (TTree*)fInAli->Get(dname);
  if(!fDigitsTree) 
    {
      LOG(AliL3Log::kError,"AliL3FileHandler::GetDigitsTree","Digits Tree")
	<<AliL3Log::kHex<<"Error getting digitstree "<<(Int_t)fDigitsTree<<ENDLOG;
      return kFALSE;
    }
  fDigitsTree->GetBranch("Segment")->SetAddress(&fDigits);
  return kTRUE;
}

void AliL3FileHandler::AliDigits2RootFile(AliL3DigitRowData *rowPt,Char_t *new_digitsfile)
{
  //Write the data stored in rowPt, into a new AliROOT file.
  //The data is stored in the AliROOT format 
  //This is specially a nice thing if you have modified data, and wants to run it  
  //through the offline reconstruction chain.
  //The arguments is a pointer to the data, and the name of the new AliROOT file.
  //Remember to pass the original AliROOT file (the one that contains the original
  //simulated data) to this object, in order to retrieve the MC id's of the digits.
  
  if(!fInAli)
    {
      printf("AliL3FileHandler::AliDigits2RootFile : No rootfile\n");
      return;
    }
  if(!fParam)
    {
      printf("AliL3FileHandler::AliDigits2RootFile : No parameter object. Run on rootfile\n");
      return;
    }
  
  //Get the original digitstree:
  fInAli->cd();
  AliTPCDigitsArray *old_array = new AliTPCDigitsArray();
  old_array->Setup(fParam);
  old_array->SetClass("AliSimDigits");
  Bool_t ok = old_array->ConnectTree("TreeD_75x40_100x60_0");
  if(!ok)
    {
      printf("AliL3FileHandler::AliDigits2RootFile : No digits tree object\n");
      return;
    }
  
  Bool_t create=kFALSE;
  TFile *digFile;
  
  digFile = TFile::Open(new_digitsfile,"NEW");
  if(digFile->IsOpen())
    {    
      create = kTRUE;
      fParam->Write(fParam->GetTitle());
    }
  else
    {
      LOG(AliL3Log::kDebug,"AliL3FileHandler::AliDigits2RootFile","Rootfile")
	<<"Rootfile did already exist, so I will just open it for updates"<<ENDLOG;
      digFile = TFile::Open(new_digitsfile,"UPDATE");
      create=kFALSE;
    }
  if(!digFile->IsOpen())
    {
      LOG(AliL3Log::kError,"AliL3FileHandler::AliDigits2RootFile","Rootfile")
	<<"Error opening rootfile "<<new_digitsfile<<ENDLOG;
      return;
    }
  
  digFile->cd();
    
  //setup a new one, or connect it to the existing one:
  AliTPCDigitsArray *arr = new AliTPCDigitsArray; 
  arr->SetClass("AliSimDigits");
  arr->Setup(fParam);
  if(create)
    arr->MakeTree();
  else
    {
      Bool_t ok = arr->ConnectTree("TreeD_75x40_100x60_0");
      if(!ok)
	{
	  printf("AliL3FileHandler::AliDigits2RootFile : No digits tree object in existing file\n");
	  return;
	}
    }
  Int_t digcounter=0;

  for(Int_t i=fRowMin; i<=fRowMax; i++)
    {
      
      if((Int_t)rowPt->fRow != i) printf("AliL3FileHandler::AliDigits2RootFile : Mismatching row numbering!!!\n");
            
      Int_t sector,row;
      AliL3Transform::Slice2Sector(fSlice,i,sector,row);
      AliSimDigits * dig = (AliSimDigits*)arr->CreateRow(sector,row);
      AliSimDigits *old_dig = (AliSimDigits*)old_array->LoadRow(sector,row);
      if(!old_dig)
	printf("AliL3FileHandler::AliDigits2RootFile : No padrow %d %d\n",sector,row);
      
      AliL3DigitData *digPt = rowPt->fDigitData;
      digcounter=0;
      for(UInt_t j=0; j<rowPt->fNDigit; j++)
	{
	  Int_t charge = (Int_t)digPt[j].fCharge;
	  Int_t pad = (Int_t)digPt[j].fPad;
	  Int_t time = (Int_t)digPt[j].fTime;
	  
	  if(charge == 0) //Only write the digits that has not been removed
	    continue;
	  digcounter++;
	  dig->SetDigitFast(charge,time,pad);
	  
	  Int_t trackID[3] = {old_dig->GetTrackID(time,pad,0),old_dig->GetTrackID(time,pad,1),old_dig->GetTrackID(time,pad,2)};
	  Int_t s_pad = pad;
	  Int_t s_time = time - 1;
	  while(trackID[0] < 0)
	    {
	      if(s_time >= 0 && s_time < AliL3Transform::GetNTimeBins() && s_pad >= 0 && s_pad < AliL3Transform::GetNPads(i))
		{
		  if(old_dig->GetTrackID(s_time,s_pad,0) > 0)
		    {
		      trackID[0]=old_dig->GetTrackID(s_time,s_pad,0); 
		      trackID[1]=old_dig->GetTrackID(s_time,s_pad,1); 
		      trackID[2]=old_dig->GetTrackID(s_time,s_pad,2); 
		    }
		}
	      if(s_pad == pad && s_time == time - 1)
		s_time = time + 1;
	      else if(s_pad == pad && s_time == time + 1)
		{s_pad = pad - 1; s_time = time;}
	      else if(s_pad == pad - 1 && s_time == time)
		s_time = time - 1;
	      else if(s_pad == pad - 1 && s_time == time - 1)
		s_time = time + 1;
	      else if(s_pad == pad - 1 && s_time == time + 1)
		{s_pad = pad + 1; s_time = time;}
	      else if(s_pad == pad + 1 && s_time == time)
		s_time = time - 1;
	      else if(s_pad == pad + 1 && s_time == time - 1)
		s_time = time + 1;
	      else 
		break;
	    }
	  
	  dig->SetTrackIDFast(trackID[0],time,pad,0);
	  dig->SetTrackIDFast(trackID[1],time,pad,1);
	  dig->SetTrackIDFast(trackID[2],time,pad,2);
	  
	}
      //cout<<"Wrote "<<digcounter<<" on row "<<i<<endl;
      UpdateRowPointer(rowPt);
      arr->StoreRow(sector,row);
      arr->ClearRow(sector,row);  
      old_array->ClearRow(sector,row);
    }
  digFile->cd();
  char treeName[100];
  sprintf(treeName,"TreeD_%s_0",fParam->GetTitle());
  printf("Writing tree to file.....");
  arr->GetTree()->Write(treeName,TObject::kOverwrite);
  printf("done\n");
  digFile->Close();
  //arr->GetTree()->Delete();
  //delete arr;
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
  
  Char_t cname[100];
  Int_t eventn = 0;
  sprintf(cname,"TreeC_TPC_%d",eventn);
  AliTPCClustersArray carray;
  carray.Setup(fParam);
  carray.SetClusterType("AliTPCcluster");
  Bool_t clusterok = carray.ConnectTree(cname);
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
    AliL3Transform::Sector2Slice(lslice,lrow,sector,row);
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
    AliL3Transform::Sector2Slice(lslice,lrow,sector,row);
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

