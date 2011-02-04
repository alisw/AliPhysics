/**************************************************************************
 * Copyright(c) 2009-2011, ALICE Experiment at CERN, All rights reserved. *
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
/*
  $Id$
*/

//-------------------------------------------------------------------------
//                          Class AliITSQASSDRefData
//                     ITS SSD reference values for the QA
//
//         Origin: Panos.Christakoglou@cern.ch, NIKHEF-Utrecht University
//-------------------------------------------------------------------------

#include <Riostream.h>
#include <fstream>
#include <TArray.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>

#include "AliLog.h"
#include "AliITSQASSDRefData.h"

ClassImp(AliITSQASSDRefData)

//___________________________________________________________________________
AliITSQASSDRefData::AliITSQASSDRefData() :
  TObject(),
  fRefList(0),
  fNameList(0) { 
  //Default constructor
}

//___________________________________________________________________________
AliITSQASSDRefData::AliITSQASSDRefData(Int_t specie) :
  TObject(),
  fRefList(0),
  fNameList(0) { 
  //Default constructor
  SetDefault(specie);
}

//___________________________________________________________________________
AliITSQASSDRefData::AliITSQASSDRefData(const char* path) :
  TObject(),
  fRefList(0),
  fNameList(0) {
  //Constructor with the path of the ascii file as an argument
  SetReferenceData(path);
}

//___________________________________________________________________________
AliITSQASSDRefData::AliITSQASSDRefData(const AliITSQASSDRefData& refData):
TObject(),
fRefList(refData.fRefList),
fNameList(refData.fNameList) {
  //Copy constructor
}

//___________________________________________________________________________
AliITSQASSDRefData& AliITSQASSDRefData::operator = (const AliITSQASSDRefData& refData) {
  //assignment operator
  if(&refData != this) {
    fRefList = refData.fRefList;
    fNameList = refData.fNameList;
  }
  return *this ;
}

//___________________________________________________________________________
AliITSQASSDRefData::~AliITSQASSDRefData() { 
  //Destructor
  if(fRefList) delete fRefList;
  if(fNameList) delete fNameList;
}

//___________________________________________________________________________
void AliITSQASSDRefData::AddReference(const char* name="", 
				      Int_t id=-1, 
				      Double_t value=0) {
  //Adding a ref. value to the list
  //Printf("(AliITSQASSDRefData::AddReference) Name: %s - Id: %d - Value: %lf",name,id,value);
  if(id>-1&&id<fRefList->GetSize()) {
    AliError(Form("Reference with id %i already exists. Choose other id or use SetReferenceValue(Int_t, Double_t) to overwrite",id));
    return;
  }
  
  if( (strcmp(name,"")!=0) && GetID(name)!=-1) {
    AliError(Form("Reference with name %s already exists. Choose other name or use SetReferenceValue(const char*, Double_t) to overwrite",name));
    return;
  }
  
  if(id==-1) id=fRefList->GetSize();
  fRefList->Set(id+1);
  fRefList->AddAt(value,id);
  fNameList->AddAt(new TObjString(name),id);
}

//___________________________________________________________________________
Int_t AliITSQASSDRefData::GetID(const char* name) {
  //Get the id of the reference value
  Int_t status = -1;
  TString refName = "";
  TString stringName = name;
  TObjString *dummyString = 0;
  for (Int_t id=0; id<fNameList->GetEntriesFast(); id++){
    dummyString = static_cast <TObjString *>(fNameList->At(id));
    refName = dummyString->GetString();
    if(refName == stringName) {
      status = id;
    }
  }

  return status;
}

//___________________________________________________________________________
Double_t AliITSQASSDRefData::GetReferenceValue(const char* name) {
  //Returns the ref. value based on the given name
  TString refName = "";
  TObjString *dummyString = 0;
  for (Int_t id=0; id<fNameList->GetEntriesFast(); id++){
    dummyString = static_cast <TObjString *>(fNameList->At(id));
    refName = dummyString->GetString();

    if(refName.Data()==name) return fRefList->At(id);
  }
  AliError(Form("Reference name %s unknown",name));
  return -1;
}

//___________________________________________________________________________
Double_t AliITSQASSDRefData::GetReferenceValue(Int_t id) {
  //Returns the ref. value based on the given id
  if (id<0||id>fRefList->GetSize()-1){
    AliError("Reference ID out of range");
    return 0;
  }
  return fRefList->At(id);

}

//___________________________________________________________________________
void AliITSQASSDRefData::PrintTable() {
  //Prints the list of reference values
  Printf("___ SSD REFERENCE DATA ___ ");
  Printf("ID ----- Value ------- Name");
  Int_t id=0;
  TString refName = "";
  TObjString *dummyString = 0;
  for(id=0; id<fRefList->GetSize()-1; id++) {
    dummyString = static_cast <TObjString *>(fNameList->At(id));
    refName = dummyString->GetString();
    Printf("%i ------ %4.3g -------- %s",id,fRefList->At(id),refName.Data());
	   
  }

}

//___________________________________________________________________________
void AliITSQASSDRefData::SetDefault(Int_t specie) {
  //Sets the default values
  if(!fNameList) fNameList = new TObjArray();
  fNameList->Add(new TObjString("minSSDDataSize"));
  fNameList->Add(new TObjString("maxSSDDataSize"));
  fNameList->Add(new TObjString("minDDLDataSize"));
  fNameList->Add(new TObjString("maxDDLDataSize"));
  fNameList->Add(new TObjString("minLDCDataSize"));
  fNameList->Add(new TObjString("maxLDCDataSize"));
  fNameList->Add(new TObjString("minMeanDDLDataSize"));
  fNameList->Add(new TObjString("maxMeanDDLDataSize"));
  fNameList->Add(new TObjString("minMeanLDCDataSize"));
  fNameList->Add(new TObjString("maxMeanLDCDataSize"));
  fNameList->Add(new TObjString("maxOccupancy"));
  fNameList->SetOwner(kTRUE);
  
  //specie == 0 ==> Default
  Double_t refValues[11] = {0,0.0,0,0,0,0,0,0,0,0};
  //specie == 1 ==> Low multiplicity
  if(specie == 1) {
    refValues[0] = 0; refValues[1] = 500; refValues[2] = 0; refValues[3] = 50;
    refValues[4] = 0; refValues[5] = 100; refValues[6] = 0; refValues[7] = 50;
    refValues[8] = 0; refValues[9] = 100; refValues[10] = 5;
  }
  //specie == 2 ==> High multiplicity
  if(specie == 2) {
    refValues[0] = 0; refValues[1] = 500; refValues[2] = 0; refValues[3] = 50;
    refValues[4] = 0; refValues[5] = 100; refValues[6] = 0; refValues[7] = 50;
    refValues[8] = 0; refValues[9] = 100; refValues[10] = 5;
  }
  //specie == 3 ==> Cosmics
  if(specie == 3) {
    refValues[0] = 0; refValues[1] = 500; refValues[2] = 0; refValues[3] = 50;
    refValues[4] = 0; refValues[5] = 100; refValues[6] = 0; refValues[7] = 50;
    refValues[8] = 0; refValues[9] = 100; refValues[10] = 5;
  }
  //specie == 4 ==> Calibration
  if(specie == 4) {
    refValues[0] = 0; refValues[1] = 500; refValues[2] = 0; refValues[3] = 50;
    refValues[4] = 0; refValues[5] = 100; refValues[6] = 0; refValues[7] = 50;
    refValues[8] = 0; refValues[9] = 100; refValues[10] = 5;
  }

  if(!fRefList) fRefList = new TArrayD();
  fRefList->Set(11,refValues);
}

//___________________________________________________________________________
void AliITSQASSDRefData::SetReferenceData(const char* path) {
  //Parses an ascii file with the reference values
  if(!fRefList) fRefList = new TArrayD();
  if(!fNameList) fNameList = new TObjArray();

  ifstream file;
  file.open(path);
  
  if (!file) {
    AliWarning(Form("No file found at %s",path));
    SetDefault(0);
    return;
  }
  if(file.bad()){
    AliWarning("Reference data could not be read: Default values are used.");
    SetDefault(0);
    return;
  }
  
  fRefList->Set(0);
  Int_t id = 0, newid = -1;
  Double_t value = 0;
  TString name = "";
  
  while (!file.eof()) {
    //file >> newid;
    file >> name >> id >> value;
    //Printf("Name: %s - Id: %d - Value: %lf",name.Data(),id,value);
    
    if (newid==id) continue; //skip line if id is the same as previous
    AddReference(name.Data(), id, value);
    newid = id;
  }

  file.close();
}

//___________________________________________________________________________
void AliITSQASSDRefData::SetReferenceValue(Int_t id, Double_t value) {
  //Adding a single reference value by id
  if(id<0||id>fRefList->GetSize()-1) {
    AliWarning(Form("Reference ID %i out of range: value not set",id));
  }
  else fRefList->SetAt(value,id);

}

//___________________________________________________________________________
void AliITSQASSDRefData::SetReferenceValue(const char* name, Double_t value) {
  //Adding a single reference value by name
  Int_t id = GetID(name);
  //Printf("Name: %s - Id: %d",name,id);
  if(id == -1) {
    AliError(Form("Reference name %s unknown: value not set",name));
    return;
  }
   
  fRefList->SetAt(value,id);
}
