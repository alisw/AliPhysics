/**************************************************************************
 * Copyright(c) 2008-2010, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

////////////////////////////////////////////////////////////////
// Author: A. Mastroserio                                     // 
// This class is used within the detector algorithm framework //
// to write and read FO scan data.                            //
////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <Riostream.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TArrayI.h>
#include <TIterator.h>
#include <TKey.h>
#include <TMath.h>
#include <TList.h>
#include "AliLog.h"
#include "AliITSOnlineSPDfoChipConfig.h"
#include "AliITSOnlineSPDfoChip.h"
#include "AliITSOnlineSPDfoInfo.h"
#include "AliITSOnlineSPDfo.h"


ClassImp(AliITSOnlineSPDfo)
//-----------------------------------------------
AliITSOnlineSPDfo::AliITSOnlineSPDfo():
fRunNr(0),
fNdacs(0),
fFileName(""),
fFile(0x0),
fInfo(0x0),
fDACnames(0x0),
fArray(0x0),
fCheckIndex(-1),
fIndex(-1),
fInitialConfiguration("")
{
// default constructor
} 

//-------------------------------------------------
AliITSOnlineSPDfo::AliITSOnlineSPDfo(TString inputfile, Int_t runNr, Int_t eqId):
fRunNr(runNr),
fNdacs(0),
fFileName(""),
fFile(0x0),
fInfo(0x0),
fDACnames(0x0),
fArray(0x0),
fCheckIndex(-1),
fIndex(-1),
fInitialConfiguration("")
{
 //
 // constructor 
 //   
    fFileName=Form("%i_%s%02i.root",runNr,inputfile.Data(),eqId);
    fFile = TFile::Open(fFileName.Data());
    fArray = new TObjArray();
    fDACnames = new THashList();
}
//--------------------------------------------------
AliITSOnlineSPDfo::AliITSOnlineSPDfo(const AliITSOnlineSPDfo &c):
fRunNr(c.fRunNr),
fNdacs(c.fNdacs),
fFileName(c.fFileName),
fFile(c.fFile),
fInfo(c.fInfo),
fDACnames(c.fDACnames),
fArray(c.fArray),
fCheckIndex(c.fCheckIndex),
fIndex(c.fIndex),
fInitialConfiguration(c.fInitialConfiguration)
{
  //
  //copy constructor
  //
}    
//--------------------------------------------------
void AliITSOnlineSPDfo::SetFile(TString inputfile)
{
  //
  // open the file where the data are
  //
  
 if(fFile) {
   fInfo = (AliITSOnlineSPDfoInfo *)fFile->Get("generalinfo");
   return; 
 } else {
   fFile = TFile::Open(inputfile.Data());
   if(!fFile) {
     Info("AliITSOnlineSPDfo::SetFile",Form(" %s  not existing.... The scan info are not available....crash is expected \n",inputfile.Data()));
     return;
    } else {
     fInfo = (AliITSOnlineSPDfoInfo *)fFile->Get("generalinfo");
     fNdacs = fInfo->GetNumDACindex();
   }
}

}

//------------------------------------------------------
void AliITSOnlineSPDfo::CreateOutputFile()
{
  //
  // Create the file (needed only in the DA ), will delete the previous one!
  //  
  if(fFile) Info("AliITSOnlineSPDfo::CreateOutputFile","removing previous file....\n");
   fFile = TFile::Open(fFileName.Data(),"RECREATE"); 
  
}
//------------------------------------------------------
void AliITSOnlineSPDfo::AddMeasurement(const TArrayS dac, Short_t measure[4], Int_t hs, Int_t chipId)
{
  /*
  // Here a single measurement is added to the chip container. 
  // A single measurement corresponds to a specific pixel-configuration output
  // in the Fast-OR chip. If N configurations are considered, then
  // the data structure is the following:
  //
  //                                          -> HS0_CHIP0 -> measure0[4]
  //                                         /                measure1[4]
  //                                        /                   ...
  //  fArray->At(i) =  DAC1-DAC2-DAC3-DAC4                    measureN[4]
  //                                        \
  //                                         \ 
  //                                          -> HS0_CHIP1 -> measure0[4]
  //                                                          measure1[4]
  //                                                             ....
  //                                                          measureN[4]
  //
  */

  AliITSOnlineSPDfoChipConfig *counts = new AliITSOnlineSPDfoChipConfig(measure);
  Int_t arrayelement = CheckDACEntry(dac);
  
  if(arrayelement< 0){
    TString dacname = CreateDACEntry(dac);
    TObjString *string = new TObjString(dacname.Data());
    fDACnames->Add(string);
     
    TObjArray *array = new TObjArray(60);
    AliITSOnlineSPDfoChip * chip = new AliITSOnlineSPDfoChip(dac.GetSize());
    chip->SetActiveHS(hs);
    chip->SetChipId(chipId); 
    for(Int_t i=0; i< dac.GetSize() ; i++) chip->SetDACParameter(i,dac.At(i));               
    chip->AddMeasurement(counts);
    array->AddAt(chip,hs*10+chipId);  
    fArray->AddLast(array); 
     
  } else {
    
    TObjArray *arr = (TObjArray*)fArray->At(arrayelement);  
    if(!arr->At(hs*10+chipId)){
      AliITSOnlineSPDfoChip * chip = new AliITSOnlineSPDfoChip(dac.GetSize());
      chip->SetActiveHS(hs);
      chip->SetChipId(chipId); 
      for(Int_t i=0; i< dac.GetSize() ; i++) chip->SetDACParameter(i,dac.At(i));               
      chip->AddMeasurement(counts);
      arr->AddAt(chip,hs*10+chipId);     
      
      } else {
      
      AliITSOnlineSPDfoChip *c = (AliITSOnlineSPDfoChip *)arr->At(hs*10+chipId);
      if(c)c->AddMeasurement(counts);    
    }   
  }
}
//---------------------------------------
Int_t AliITSOnlineSPDfo::CheckDACEntry(const TArrayS dac)
{
  //
  // Check if the set of dacs has been already added to the array
  // 
  
   
  TString name = CreateDACEntry(dac);
  
  if(!fDACnames) {
    Info("AliITSOnlineSPDfo::CheckDACEntry"," NO DAC name array is present, exiting.... \n");
    return -1;
  }
  
   Double_t c = 0;
  for(Int_t i=0; i< fNdacs; i++)  c+=dac.At(i)*TMath::Power(10,3*i);
  if(c==fCheckIndex) return fIndex;
   
  
 
   Int_t idx = -1;
   
  if(fDACnames->FindObject( name.Data() )) {
    idx = fDACnames->IndexOf( fDACnames->FindObject( name.Data() )  );
    fCheckIndex=0;
    for(Int_t i=0; i< fNdacs; i++) fCheckIndex+=dac.At(i)*TMath::Power(10,3*i);
    fIndex=idx;

  }
  return idx;

}
//_________________________________________________________
void AliITSOnlineSPDfo::WriteToFile()
{
  //
  //The array of DACS and all its content is written to file.
  // Here the general info on the FO calibration scan are 
  // written in the same file
  
  if(fDACnames->GetEntries() != fArray->GetEntries()) {
   printf("mismatch names-array. Exiting....");   
    return;
  }
  
  for(Int_t i=0; i< fDACnames->GetEntries(); i++){
  fFile->WriteObject(fArray->At(i),fDACnames->At(i)->GetName());
  }
  fFile->WriteTObject(fInfo,"generalinfo"); 
  fFile->Close();
  
}
//______________________________________________________

TString AliITSOnlineSPDfo::CreateDACEntry(const TArrayS dacs) const
{
  //
  // The string of DACs is build
  // 
  
  TString dacvalues;
  for(Int_t i=0; i<dacs.GetSize(); i++) dacvalues+=Form("-%i",dacs.At(i));
  dacvalues.Remove(0,1);
  return dacvalues;
 
}
//_____________________________________________________
TArrayI AliITSOnlineSPDfo::GetDACscanParams() const
{ 
  //
  // this method retrieves the DAC value range and its steps
  //
  TFile *f =0x0;
  if(fFile->IsOpen()) f = fFile;
  else f = TFile::Open(fFileName.Data());
  
  TArrayI dacs;
  if(f->GetNkeys() < 2) return dacs;
  dacs.Set(fNdacs*3);
  
  TArrayI min(fNdacs), max(fNdacs),step(fNdacs), check(fNdacs), refvalues(fNdacs);
  
  for(Int_t i=0; i< fNdacs ; i++) {
    min.AddAt(9999,i);
    max.AddAt(0,i);
    step.AddAt(0,i);
    check.AddAt(0,i);
  }
  
  
  TKey *key;
  TIter iter(f->GetListOfKeys());  
  while ((key = (TKey*)(iter.Next()))) {
    TString classname = key->GetClassName();
    if(classname.Contains("OnlineSPD")) break;
    
    TString values = key->GetName();
    Int_t *val = GetDACvalues(values,fNdacs); // the user has to delete it!
      
    for(Int_t i=0; i< fNdacs; i++) {
      if(val[i]<=min.At(i)) min.AddAt(val[i],i);
      if(val[i] >= max.At(i)) max.AddAt(val[i],i);
      // procedure to get the step size;
      if(!check.At(i)) { 
	refvalues.AddAt(val[i],i);
	check.AddAt(1,i);
      }
      if(step.At(i) ==0 && check.At(i)){
        if(val[i]!=refvalues.At(i)) step.AddAt(TMath::Abs(refvalues.At(i) - val[i]),i);
      }  
      // end procedure to get the step size 
    }   
    delete [] val;
  }

 
 
  for(Int_t i=0; i<fNdacs; i++) {
    dacs.AddAt(min.At(i),3*i);    
    dacs.AddAt(max.At(i),3*i+1); 
    dacs.AddAt(step.At(i),3*i+2); 
  }
 
  
  return dacs;     
}
//___________________________________________________________________
Int_t* AliITSOnlineSPDfo::GetDACvalues(TString s, const Int_t ndacs) const
{
  //
  // Translates the string of DACS values into an array of integers
  //
  
   Int_t *val = new Int_t[ndacs]; 
   
   s.ReplaceAll("-"," ");
   char *pEnd = Form("%s",s.Data());
   for(Int_t i=0; i< ndacs; i++) {
     val[i] = strtol(pEnd,&pEnd,10); // conversion from string to long
   }  
   return val;   
}

//___________________________________________________________________
Double_t* AliITSOnlineSPDfo::GetDACvaluesD(TString s, const Int_t ndacs) const 
{  
  //
  // Translates the string of DACS values into an array of doubles
  // (needed to fill the thnsparse)
  //
   Double_t *val = new Double_t[ndacs];
   Int_t *values =  GetDACvalues(s,ndacs);
   for(Int_t i=0; i< ndacs; i++) val[i] = (Double_t)values[i];
   delete [] values;
   
   return val;   
}
//-------------------------------------------------------------------
TArrayS AliITSOnlineSPDfo::CreateDACArray(const TArrayS dacs, const TArrayS dacId) const
{
  //
  // method to order the data according to the DAC index
  //
  
  TArrayS dacarray(dacs.GetSize());
  for(Int_t i=0; i<dacs.GetSize(); i++) {
    if(dacId.At(i)==kIdFOPOL) dacarray.AddAt(dacs.At(i),kFOPOL);
    else if(dacId.At(i)==kIdCONVPOL) dacarray.AddAt(dacs.At(i),kCONVPOL);
    else if(dacId.At(i)==kIdCOMPREF) dacarray.AddAt(dacs.At(i),kCOMPREF);
    else if(dacId.At(i)==kIdPreVTH) dacarray.AddAt(dacs.At(i),kPreVTH);
    else if(dacId.At(i)==kIdCGPOL) dacarray.AddAt(dacs.At(i),kCGPOL);
    else printf("new DAC included in the scan??\n");
    
  }
 
  return dacarray;
  
}
  
  
