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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for T0 calibration                       TM-AC-AM_6-02-2006         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <TCanvas.h>

#include "AliT0CalibData.h"
#include "AliT0LookUpValue.h"
#include "TObjArray.h"
#include "TGraph.h"
#include "TFile.h"
#include "AliLog.h"
#include "TObjString.h"

#include "TAxis.h"
#include "TH2F.h"



ClassImp(AliT0CalibData)

//________________________________________________________________
  AliT0CalibData::AliT0CalibData():   TNamed()

{
  //
}

//________________________________________________________________
AliT0CalibData::AliT0CalibData(const char* name):TNamed(),fTimeDelayTVD(0),fWalk(),fSlewingLED(),fSlewingRec()
{
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());

}

//________________________________________________________________
AliT0CalibData::AliT0CalibData(const AliT0CalibData& calibda) :
  TNamed(calibda),fTimeDelayTVD(0),fWalk(),fSlewingLED(),fSlewingRec()
{
// copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());


}

//________________________________________________________________
AliT0CalibData &AliT0CalibData::operator =(const AliT0CalibData& calibda)
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
 
  return *this;
}

//________________________________________________________________
AliT0CalibData::~AliT0CalibData()
{
  //
}
//________________________________________________________________
void AliT0CalibData::Reset()
{
    memset(fTimeDelayCFD,1,24*sizeof(Float_t));
    memset(fTimeDelayLED,1,24*sizeof(Float_t));
    memset(fGain,1,24*sizeof(Float_t));
}


//________________________________________________________________
void  AliT0CalibData::Print(Option_t*) const
{

  printf("\n	----	PM Arrays	----\n\n");
  printf(" Time delay CFD \n");
  for (Int_t i=0; i<24; i++) printf("  %f",fTimeDelayCFD[i]);
  printf(" \n LED \n");
  for (Int_t i=0; i<24; i++) printf("  %f",fTimeDelayLED[i]);
  printf(" \n Gain \n");
  for (Int_t i=0; i<24; i++) printf("  %f",fGain[i]);
  printf(" \n");
} 

//________________________________________________________________
void  AliT0CalibData::PrintLookup(Option_t*, Int_t iTRM, Int_t iTDC, Int_t iChannel) const
{

   AliT0LookUpKey* lookkey= new AliT0LookUpKey();
   AliT0LookUpValue*  lookvalue= new AliT0LookUpValue();
 
     lookvalue->SetTRM(iTRM);
     lookvalue->SetTDC(iTDC);
     lookvalue->SetChain(0);
     lookvalue->SetChannel(iChannel);

     printf(" AliT0CalibData::PrintLookup ::start GetValue %i %i %i \n",iTRM, iTDC, iChannel);
     lookkey = (AliT0LookUpKey*) fLookup.GetValue((TObject*)lookvalue);
     
     cout<<"  AliT0CalibData::PrintLookup :: lookkey "<< lookkey<<endl;
     if (lookkey)
       {
	 cout<<" lookup KEY!!! "<<lookkey->GetKey()<<" VALUE "<<lookvalue->GetTRM()<<" "
	     <<lookvalue->GetTDC()<<" "
	     << lookvalue->GetChain()<<" "
	     <<lookvalue->GetChannel()<<endl;
       }
     

}

//________________________________________________________________
void AliT0CalibData::SetTimeDelayCFD(Float_t* TimeDelay)
{
  if(TimeDelay) for(int t=0; t<24; t++) fTimeDelayCFD[t] = TimeDelay[t];
  //  else for(int t=0; t<24; t++) fTimeDelay[t] = 0.;
}
//________________________________________________________________
void AliT0CalibData::SetTimeDelayLED(Float_t* TimeDelay)
{
  if(TimeDelay) for(int t=0; t<24; t++) fTimeDelayLED[t] = TimeDelay[t];
  //  else for(int t=0; t<24; t++) fTimeDelay[t] = 0.;
}

//________________________________________________________________
void AliT0CalibData::SetGain(Float_t* Gain)
{
  if(Gain) for(int t=0; t<24; t++) fGain[t] = Gain[t];
  // else for(int t=0; t<24; t++) fGain[t] = 0.;
}


//________________________________________________________________
void AliT0CalibData::SetWalk(Int_t ipmt, const Char_t *filename)
{

  TFile *file = new TFile(filename);
  char funcname[256];
  sprintf(funcname,"CFD%i",ipmt+1);
  TF1* gr = (TF1*)file->Get(funcname);
  fWalk.AddAtAndExpand(gr,ipmt);
  file->Close();
}


//________________________________________________________________

void AliT0CalibData::SetSlewingLED(Int_t ipmt,const Char_t *filename)
{
  Float_t mv, ps; 
  Float_t x[100], y[100];
  string buffer;
  
  ifstream inFile(filename);
  if(!inFile) {AliError(Form("Cannot open file %s !",filename));}
  
  inFile >> mv>>ps;
  Int_t i=0;
  
  while(getline(inFile,buffer)){
    x[i]=mv; y[i]=ps;	
    inFile >> mv >> ps;
    i++;
  }
  inFile.close();
  TGraph* gr = new TGraph(i,x,y);
  fSlewingLED.AddAtAndExpand(gr,ipmt);
   
}

//________________________________________________________________

void AliT0CalibData::SetSlewingRec(Int_t ipmt,const Char_t *filename)
{
  Float_t mv, ps; 
  Float_t x[100], y[100];
  string buffer;
  
  ifstream inFile(filename);
  if(!inFile) {AliError(Form("Cannot open file %s !",filename));}
  
  inFile >> mv>>ps;
  Int_t i=0;
  
  while(getline(inFile,buffer)){
    x[i]=mv; y[i]=ps;	
    inFile >> mv >> ps;
    i++;
  }
  inFile.close();
  Float_t y1[100], x1[100];
  for (Int_t ir=0; ir<i; ir++){
    y1[ir]=y[i-ir]; x1[ir]=x[i-ir];}
  TGraph* gr = new TGraph(i,y1,x1);
  fSlewingRec.AddAtAndExpand(gr,ipmt);
  
}


void AliT0CalibData::ReadAsciiLookup(const Char_t *filename)
{

  if(filename == 0){
    AliError(Form("Please, specify file with database")) ;
    return ;
  }

  //  AliT0LookUpKey * lookkey= new AliT0LookUpKey();
  //AliT0LookUpValue * lookvalue= new AliT0LookUpValue();

  ifstream lookup;
  lookup.open(filename);
  if(!lookup) {AliError(Form("Cannot open file %s !",filename));}
  Char_t varname[11];
  Int_t key, trm, tdc, chain, channel;
  // while(lookup.eof())

 for (Int_t i=0; i<108; i++)
//    for (Int_t i=0;i<3;i++)
    {
       AliT0LookUpKey * lookkey= new AliT0LookUpKey();
       AliT0LookUpValue * lookvalue= new AliT0LookUpValue();

 lookup>>varname>>key>>trm>>chain>>tdc>>channel;
      lookvalue->SetTRM(trm);
      lookvalue->SetTDC(tdc);
      lookvalue->SetChain(chain);
      lookvalue->SetChannel(channel);
      lookkey->SetKey(key);

      //      cout<<"TRM="<<trm<<" TDC="<<tdc<<" chain="<<chain<<" channel="<<channel<<" key="<<key<<endl;

      cout<<"AliT0CalibData:: "<<varname<<" "<<key<<" "<<trm<<" "<<chain<<" "<<tdc<<" "<<channel<<endl;

      //    fLookup.Add((TObject*)lookkey,(TObject*)lookvalue);
      // ar_key.AddAt(lookkey,i);
      // ar_val.AddAt(lookvalue,i);

        fLookup.Add((TObject*)lookvalue,(TObject*)lookkey);

    }

   lookup.close();

}

//________________________________________________________________

Int_t AliT0CalibData::GetChannel(Int_t trm,  Int_t tdc, Int_t chain, Int_t channel)
{

  AliT0LookUpKey * lookkey;//= new AliT0LookUpKey();
  AliT0LookUpValue * lookvalue= new AliT0LookUpValue(trm,tdc,chain,channel);

     lookkey = (AliT0LookUpKey*) fLookup.GetValue((TObject*)lookvalue);
    cout<<"AliT0CalibData:: key "<<lookkey->GetKey()<<endl;
  return lookkey->GetKey();

}

