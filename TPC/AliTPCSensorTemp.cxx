/**************************************************************************
 * Copyright(c) 2006-07, ALICE Experiment at CERN, All rights reserved. *
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


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Class describing TPC temperature sensors (including pointers to graphs/fits//
// Authors: Marian Ivanov, Haavard Helstrup and Martin Siska                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// Running instructions:
/*
  TClonesArray * arr = AliTPCSensorTemp::ReadList("TempSensor.txt");
  TFile f("TempSensors.root","RECREATE");
  TTree * tree = new TTree("TempSensor", "TempSensor");
  tree->Branch("Temp",&arr);
  tree->Fill();
  tree->Write();
  
 */
//


#include "AliTPCSensorTemp.h"
ClassImp(AliTPCSensorTemp)


AliTPCSensorTemp::AliTPCSensorTemp(): AliDCSSensor(),
  fType(0),
  fSide(0),
  fSector(0),
  fNum(0)
{
  //
  //  Standard constructor
  //
}

AliTPCSensorTemp::AliTPCSensorTemp(const AliTPCSensorTemp& source) :
  AliDCSSensor(source),
   fType(source.fType),
   fSide(source.fSide),
   fSector(source.fSector),
   fNum(source.fNum)

//
//  Copy constructor
//
{ }

AliTPCSensorTemp& AliTPCSensorTemp::operator=(const AliTPCSensorTemp& source){
//
// assignment operator
//
  if (&source == this) return *this;
  new (this) AliTPCSensorTemp(source);
  
  return *this;  
}

TClonesArray * AliTPCSensorTemp::ReadList(const char *fname) {
   
   Int_t firstSensor, lastSensor;
   return ReadListInd(fname,firstSensor,lastSensor);
}  

TClonesArray * AliTPCSensorTemp::ReadListInd(const char *fname, 
                                          Int_t& firstSensor,
					  Int_t& lastSensor) {
  //
  // read values from ascii file
  //
  TTree * tree = new TTree("asci","asci");
  tree->ReadFile(fname,"");
  
  Int_t nentries = tree->GetEntries();
  Int_t sensor=0;
  Int_t sector=0;
  char  type[100];
  char  side[100];
  Int_t num=0;
  Int_t echa=0;
  //Double_t x=0;
  //Double_t y=0;
  //Double_t z=0;
  //String_t namedtp[100];

  tree->SetBranchAddress("Sensor",&sensor);
  tree->SetBranchAddress("Type",&type);
  tree->SetBranchAddress("Side",&side);
  tree->SetBranchAddress("Sec",&sector);
  tree->SetBranchAddress("Num",&num);
  tree->SetBranchAddress("ECha",&echa);
  //tree->SetBranchAddress("X",&x);
  //tree->SetBranchAddress("Y",&y);
  //tree->SetBranchAddress("Z",&z);

  firstSensor = (Int_t)tree->GetMinimum("ECha");
  lastSensor = (Int_t)tree->GetMaximum("ECha");

  TClonesArray * array = new TClonesArray("AliTPCSensorTemp",nentries);

  for (Int_t isensor=0; isensor<nentries; isensor++){
    AliTPCSensorTemp * temp = new ((*array)[isensor])AliTPCSensorTemp;
    tree->GetEntry(isensor);
    temp->SetId(sensor);
    temp->SetIdDCS(echa);
    if (side[0]=='C') temp->SetSide(1);
    temp->SetSector(sector);
    temp->SetNum(num);
    //temp->SetType(type);
    if (bcmp(type,"ROC",3)==0) temp->SetType(0);
    if (bcmp(type,"OFC",3)==0) temp->SetType(1);
    if (bcmp(type,"IFC",3)==0) temp->SetType(2);
    if (bcmp(type,"TPC",3)==0) temp->SetType(3); 
    if (bcmp(type,"ELM",3)==0) temp->SetType(4);
    if (bcmp(type,"TS",2)==0)  temp->SetType(5);
    if (bcmp(type,"COOL",3)==0)temp->SetType(6);
    //temp->SetX(x);

    if (temp->GetType()==0){
	temp->SetX(TMath::Cos((2*sector+1)*0.1745)*(83+(num+1)*30));
    }
    if ((temp->GetType()==1) || (temp->GetType()==4)){
      temp->SetX(TMath::Cos((2*sector+1)*0.1745)*260);
    }
    if ((temp->GetType()==2) || (temp->GetType()==3)){
      temp->SetX(TMath::Cos((2*sector+1)*0.1745)*83);
    }
    if ((temp->GetType()==5) || (temp->GetType()==6)){
      temp->SetX(0);
    }
    
    //temp->SetY(y);
    if (temp->GetType()==0){
	  temp->SetY(TMath::Sin((2*sector+1)*0.1745)*(83+(num+1)*30));
    }
    if ((temp->GetType()==1) || (temp->GetType()==4)){
      temp->SetY(TMath::Sin((2*sector+1)*0.1745)*260);
    }
    if ((temp->GetType()==2) || (temp->GetType()==3)){
      temp->SetY(TMath::Sin((2*sector+1)*0.1745)*83);
    }
    if ((temp->GetType()==5) || (temp->GetType()==6)){
      temp->SetY(0);
    }
    //temp->SetZ(z);
    if ((temp->GetType()==0 || temp->GetType()==3 || temp->GetType()==4 || temp->GetType()==5 || temp->GetType()==6) && temp->GetSide()==0) {
      temp->SetZ(260);
      }
    if ((temp->GetType()==0 || temp->GetType()==3 || temp->GetType()==4 || temp->GetType()==5 || temp->GetType()==6) && temp->GetSide()==1){
      temp->SetZ(-260);
      }
    if((temp->GetType()==1 || temp->GetType()==2) && (num==0)) {
      temp->SetZ(250);
      }
    if((temp->GetType()==1 || temp->GetType()==2) && (num==1)) {
      temp->SetZ(180);
      }
    if((temp->GetType()==1 || temp->GetType()==2) && (num==2)) {
      temp->SetZ(60);
      }
    if((temp->GetType()==1 || temp->GetType()==2) && (num==3)) {
      temp->SetZ(-60);
      }
    if((temp->GetType()==1 || temp->GetType()==2) && (num==4)) {
      temp->SetZ(-180);
      }
    if((temp->GetType()==1 || temp->GetType()==2) && (num==5)) {
      temp->SetZ(-250);
      }
  }
  delete tree;  
  return array;
}
