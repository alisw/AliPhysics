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
  TClonesArray * arr = AliTPCSensorTemp::ReadList("TempSensor.txt","tpc_PT_%d.Temperature");
  TFile f("TempSensors.root","RECREATE");
  TTree * tree = new TTree("TempSensor", "TempSensor");
  tree->Branch("Temp",&arr);
  tree->Fill();
  tree->Write();
  
*/

//

#include <strings.h>
#include "AliTPCSensorTemp.h"
ClassImp(AliTPCSensorTemp)



const Float_t kASideX[18][5]={
        { 99.56,  117.59,  160.82,  186.92,  213.11},
	{ 87.56,  103.4,   141.42,  164.37,  187.41},
	{ 64.99,   76.75,  104.97,  122.00,  139.1},
	{ 34.58,   40.84,   55.85,   64.92,  74.01},
	{    0,    0,    0,    0,    0},
	{-34.58,  -40.84,  -55.85,  -64.92,  -74.01},
	{-64.99,  -76.75, -104.97, -122.0,  -139.1},
	{-87.56, -103.4,  -141.42, -164.37, -187.41},
	{-99.56, -117.59, -160.82, -186.92, -213.11},
	{-99.56, -117.59, -160.82, -186.92, -213.11},
	{-87.56, -103.4,  -141.42, -164.37, -187.41},
	{-64.99,  -76.75, -104.97, -122,    -139.1},
	{-34.58,  -40.84,  -55.85,  -64.92,  -74.01},
	{    0,    0,    0,   0,   0},
	{ 34.58,   40.84,   55.85,   64.92,   74.01},
	{ 64.99,   76.75,  104.97,  122,     139.1},
	{ 87.56,  103.4,   141.42,  164.37,  187.41},
	{ 99.56,  117.59,  160.82,  186.92,  213.11}};
	
const Float_t kASideY[18][5]={
        { 17.56,   20.73,   28.36,   32.96,   37.58},
	{ 50.55,   59.7,    81.65,   94.9,   108.2},
	{ 77.45,   91.47,  125.1,   145.4,   165.77},
	{ 95.0,  112.3,    153.45,  178.35,  203.35},
	{101.1,  119.4,  163.3,  189.8,  216.4},
	{ 95.0,  112.2,  153.45,  178.35,  203.35},
	{ 77.45,   91.47,  125.1,  145.4,  165.77},
	{ 50.55,   59.7,   81.65,   94.9,  108.2},
	{ 17.56,   20.73,  28.36,   32.96,  37.58},
	{-17.56,  -20.73, -28.36,  -32.96, -37.58},
	{-50.55,  -59.7,  -81.65,  -94.9, -108.2},
	{-77.45,  -91.47, -125.1, -145.4, -165.77},
	{-95.0, -112.2, -153.45, -178.35, -203.35},
        {-101.1, -119.4, -163.3,  -189.8, -216.4},
        {-95.0, -112.2, -153.45, -178.35, -203.35},
	{-77.45,  -91.47, -125.1, -145.4, -165.77},
	{-50.55,  -59.7,  -81.65,  -94.9, -108.2},
	{-17.56,  -20.73, -28.36,  -32.96, -37.58}};  
	
const Float_t kCSideX[18][5]={
        { 99.56,  117.59,  160.82,  186.92,  213.11},
	{ 87.56,  103.4,   141.42,  164.37,  187.41},
	{ 64.99,   76.75,  104.97,  122,     139.1},
	{ 34.58,   40.84,   55.85,   64.92,   74.01},
	{    0,    0,    0,   0,   0},
	{-34.58,  -40.84,  -55.85,  -64.92,  -74.01},
	{-64.99,  -76.75, -104.97, -122,    -139.1},
	{-87.56, -103.4,  -141.42, -164.37, -187.41},
	{-99.56, -117.59, -160.82, -186.92, -213.11},
	{-99.56, -117.59, -160.82, -186.92, -213.11},
	{-87.56, -103.4,  -141.42, -164.37, -187.41},
	{-64.99,  -76.75, -104.97, -122,    -139.1},
	{-34.58,  -40.84,  -55.85,  -64.92,  -74.01},
	{    0,    0,    0,    0,    0},
	{ 34.58,   40.84,   55.85,   64.92,   74.01},
	{ 64.99,   76.75,  104.97,  122,     139.1},
	{ 87.56,  103.4,   141.42,  164.37,  187.41},
	{ 99.56,  117.59,  160.82,  186.92,  213.11}};

const Float_t kCSideY[18][5]={
        { 17.56,   20.73,   28.36,   32.96,   37.58},
	{ 50.55,   59.7,    81.65,   94.9,   108.2},
	{ 77.45,   91.47,  125.1,   145.4,   165.77},
	{ 95.0,   112.2,   153.54,  178.35,  203.35},
	{101.1,  119.4,  163.3,  189.8,  216.4},
	{ 95.0,   112.2,   153.45,  178.35,  203.35},
	{ 77.45,   91.47,  125.1,   145.4,   165.77},
	{ 50.55,   59.7,    81.65,   94.9,   108.2},
	{ 17.56,   20.73,   28.36,   32.96,   37.58},
	{-17.56,  -20.73,  -28.36,  -32.96,  -37.58},
	{-50.55,  -59.7,   -81.56,  -94.9,  -108.2},
	{-77.45,  -91.47, -125.1,  -145.4,  -165.77},
	{-95.0,  -112.2,  -153.45, -178.35, -203.35},
        {-101.1, -119.4, -163.3, -189.8, -216.4},
        {-95.0, -112.2, -153.45, -178.35, -203.35},
	{-77.45,  -91.47, -125.1, -145.4, -165.77},
	{-50.55,  -59.7,  -81.65,  -94.9, -108.2},
	{-17.56,  -20.73, -28.36,  -32.96,  -37.58}};  

const Float_t kIFCrad[5] = {67.2, 64.4, 60.7, 64.4, 67.2};

const Float_t kTSrad[4] =  {67.2, 61.5, 67.2, 61.5}; 
const Float_t kTSz[4] =  {240.0, 90.0, 240.0, 90.0}; 

//______________________________________________________________________________________________

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
//______________________________________________________________________________________________

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
//______________________________________________________________________________________________

AliTPCSensorTemp& AliTPCSensorTemp::operator=(const AliTPCSensorTemp& source){
//
// assignment operator
//
  if (&source == this) return *this;
  new (this) AliTPCSensorTemp(source);
  
  return *this;  
}
//______________________________________________________________________________________________

TClonesArray * AliTPCSensorTemp::ReadList(const char *fname,
                                          const TString& amandaString) {
  //
  // read values from ascii file
  //
  TTree * tree = new TTree("asci","asci");
  tree->ReadFile(fname,"");
  TClonesArray *arr = ReadTree(tree, amandaString);
  delete tree;
  return arr;
}
     
//______________________________________________________________________________________________

TClonesArray * AliTPCSensorTemp::ReadTree(TTree *tree, 
                                          const TString& amandaString) {
  
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

  // firstSensor = (Int_t)tree->GetMinimum("ECha");
  // lastSensor = (Int_t)tree->GetMaximum("ECha");

  TClonesArray * array = new TClonesArray("AliTPCSensorTemp",nentries);

  for (Int_t isensor=0; isensor<nentries; isensor++){
    AliTPCSensorTemp * temp = new ((*array)[isensor])AliTPCSensorTemp;
    tree->GetEntry(isensor);
    temp->SetId(sensor);
    temp->SetIdDCS(echa);
    TString stringID = Form (amandaString.Data(),echa);
    temp->SetStringID(stringID);
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
//	temp->SetX(TMath::Cos((2*sector+1)*0.1745)*(83+(num+1)*30));
      if (side[0]=='C') {
          temp->SetX(kCSideX[sector][num]);
      } else {
          temp->SetX(kASideX[sector][num]);
      }      
    }
    if ((temp->GetType()==1) || (temp->GetType()==4)){
      temp->SetX(TMath::Cos((2*sector+1)*0.1745)*278);
    }
    if (temp->GetType()==2) {
      temp->SetX(TMath::Cos((2*sector+1)*0.1745)*60.7);
    }
    if (temp->GetType()==3) {
      if (num==0) {
        temp->SetX(TMath::Cos((2*sector+1)*0.1745)*87.5);
      } else {
        temp->SetX(TMath::Cos((2*sector+1)*0.1745)*241.8);
      }
    } 
    if (temp->GetType()==5){
      temp->SetX(TMath::Cos(sector*0.524+(num+1)*0.131)*kTSrad[num]);
    }
    if (temp->GetType()==6){
      temp->SetX(0);
    }
    
    //temp->SetY(y);
    if (temp->GetType()==0){
//	  temp->SetY(TMath::Sin((2*sector+1)*0.1745)*(83+(num+1)*30));
      if (side[0]=='C') {
          temp->SetY(kCSideY[sector][num]);
      } else {
          temp->SetY(kASideY[sector][num]);
      }      
    }
    if ((temp->GetType()==1) || (temp->GetType()==4)){
      temp->SetY(TMath::Sin((2*sector+1)*0.1745)*278);
    }
    if (temp->GetType()==2){
      temp->SetY(TMath::Sin((2*sector+1)*0.1745)*60.7);
    }
    if (temp->GetType()==3) {
      if (num==0) {
        temp->SetY(TMath::Sin((2*sector+1)*0.1745)*87.5);
      } else {
        temp->SetY(TMath::Sin((2*sector+1)*0.1745)*241.8);
      }
    } 

    if (temp->GetType()==5){
      temp->SetY(TMath::Sin(sector*0.524+(num+1)*0.131)*kTSrad[num]);
    }

    if (temp->GetType()==6){
      temp->SetY(0);
    }

    //temp->SetZ(z);
    if ((temp->GetType()==0 || 
         temp->GetType()==4 || temp->GetType()==6) && 
	 temp->GetSide()==0) {
             temp->SetZ(250);
      }
    if ((temp->GetType()==0 || temp->GetType()==4 ||
         temp->GetType()==6) && temp->GetSide()==1){
             temp->SetZ(-250);
      }
    if(temp->GetType()==2 && temp->GetSide()==0) {
         temp->SetZ(52.4);
      }
    if(temp->GetType()==2 && temp->GetSide()==1) {
         temp->SetZ(-52.4);
      }

    if(temp->GetType()==3 && temp->GetSide()==0) {
         temp->SetZ(247);
      }
    if(temp->GetType()==3 && temp->GetSide()==1) {
         temp->SetZ(-247);
      }

    if((temp->GetType()==1 ) && (num==0)) {
      temp->SetZ(240);
      }
    if((temp->GetType()==1 ) && (num==1)) {
      temp->SetZ(168.4);
      }
    if((temp->GetType()==1 ) && (num==2)) {
      temp->SetZ(51);
      }
    if((temp->GetType()==1 ) && (num==3)) {
      temp->SetZ(-51);
      }
    if((temp->GetType()==1 ) && (num==4)) {
      temp->SetZ(-168.4);
      }
    if((temp->GetType()==1 ) && (num==5)) {
      temp->SetZ(-240);
      }

    if ( num < (Int_t)(sizeof(kTSz)/sizeof(kTSz[0]))) {
      if(temp->GetType()==5 && temp->GetSide()==0) {
         temp->SetZ(kTSz[num]);
      }
      if(temp->GetType()==5 && temp->GetSide()==1) {
         temp->SetZ(-kTSz[num]);
      }
    }

  }
  return array;
}
