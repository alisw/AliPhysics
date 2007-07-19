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

#include "AliT0CalibData.h"
#include "AliT0LookUpValue.h"
#include "AliLog.h"

#include <TCanvas.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TGraph.h>
#include <TFile.h>
#include <TAxis.h>
#include <TH2F.h>
#include <TMath.h>
#include <TSystem.h>
#include <Riostream.h>

#include <string>

ClassImp(AliT0CalibData)

//________________________________________________________________
  AliT0CalibData::AliT0CalibData():   TNamed(),
				      fTimeDelayTVD(0),
				      fMeanT0(0),
				      fWalk(0),
				      fAmpLEDRec(0),
				      fLookup(0),
				      fNumberOfTRMs(0)

{
  //
}

//________________________________________________________________
AliT0CalibData::AliT0CalibData(const char* name):TNamed(),
				      fTimeDelayTVD(0),
				      fMeanT0(0),
				      fWalk(0),
				      fAmpLEDRec(0),
				      fLookup(0),
				      fNumberOfTRMs(0)
{
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());

}

//________________________________________________________________
AliT0CalibData::AliT0CalibData(const AliT0CalibData& calibda) :
  TNamed(calibda),		
  fTimeDelayTVD(0),
  fMeanT0(0),
  fWalk(0),
  //  fAmpLED(0),
  fAmpLEDRec(0),
  fLookup(0),
  fNumberOfTRMs(0)

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
    memset(fTimeDelayDA,1,24*sizeof(Float_t));
}


//________________________________________________________________
void  AliT0CalibData::Print(Option_t*) const
{

  printf("\n	----	PM Arrays	----\n\n");
  printf(" Time delay CFD & DA\n");
  for (Int_t i=0; i<24; i++) printf(" CFD  %f DA %f ",fTimeDelayCFD[i], fTimeDelayDA[i]);
} 

//________________________________________________________________
void  AliT0CalibData::PrintLookup(Option_t*, Int_t iTRM, Int_t iTDC, Int_t iChannel) const
{
  
  AliT0LookUpKey* lookkey= new AliT0LookUpKey();
  AliT0LookUpValue*  lookvalue= new AliT0LookUpValue();

  cout<<" Number Of TRMs in setup "<<GetNumberOfTRMs()<<endl;
  iTRM=0; iTDC=0; Int_t chain=0; iChannel=0;

  for (Int_t ik=0; ik<105; ik++){
    lookvalue->SetTRM(iTRM);
    lookvalue->SetTDC(iTDC);
    lookvalue->SetChain(chain);
    lookvalue->SetChannel(iChannel);
    
    if (iChannel<6) iChannel +=2;
    else {iChannel = 0; iTDC++;}
    if(ik==57) { iTDC=0; iChannel=0; iTRM=1;}
   
  printf(" AliT0CalibData::PrintLookup ::start GetValue %i %i %i %i\n",iTRM, iTDC,chain, iChannel);
    lookkey = (AliT0LookUpKey*) fLookup.GetValue((TObject*)lookvalue);
    cout<<"  lookkey "<< lookkey<<endl;
    //    TString name= lookkey->GetChannelName();
    // cout<<name.Data()<<endl;
    if (lookkey)
      {
	TString name= lookkey->GetChannelName();
	cout<<" lookup KEY!!! "<<name.Data()<<" "<<lookkey->GetKey()<<" VALUE "<<lookvalue->GetTRM()<<" "
	    <<lookvalue->GetTDC()<<" "
	    << lookvalue->GetChain()<<" "
	  <<lookvalue->GetChannel()<<endl;
      }
  }
  
}

//________________________________________________________________
void AliT0CalibData::SetTimeDelayCFD(Float_t* TimeDelay)
{
  if(TimeDelay) for(int t=0; t<24; t++) fTimeDelayCFD[t] = TimeDelay[t];
}  
  //________________________________________________________________
  void AliT0CalibData::SetTimeDelayDA(Float_t* TimeDelay)
{
  if(TimeDelay) for(int t=0; t<24; t++) fTimeDelayDA[t] = TimeDelay[t];
}


//________________________________________________________________
void AliT0CalibData::SetWalk(Int_t ipmt)
{

  Int_t mv, ps; 
  Int_t x[70000], y[70000], index[70000];
  Float_t time[10000],amplitude[10000];
  string buffer;
  Bool_t down=false;
  
  const char * filename = gSystem->ExpandPathName("$ALICE_ROOT/T0/data/CFD-Amp.txt");
  ifstream inFile(filename);
  if(!inFile) AliError(Form("Cannot open file %s !",filename));
  
  Int_t i=0;
  while(getline(inFile,buffer)){
    inFile >> ps >> mv;

    x[i]=ps; y[i]=mv;
    i++;
  }
  inFile.close();
  cout<<" number of data "<<i<<endl;
 
  TMath::Sort(i, y, index,down);
  Int_t amp=0, iin=0, isum=0, sum=0;
  Int_t ind=0;
  for (Int_t ii=0; ii<i; ii++)
    {
      ind=index[ii];
      if(y[ind] == amp)
	{
	  sum +=x[ind];
	  iin++;
	}
      else
	{
	  if(iin>0)
	    time[isum] = Float_t (sum/(iin));
	  else
	    time[isum] =Float_t (x[ind]);
	  amplitude[isum] = Float_t (amp);
	  amp=y[ind];
	  iin=0;
	  isum++;
	  sum=0;
	}
    }

  inFile.close();

  TGraph* gr = new TGraph(isum, amplitude, time);
  fWalk.AddAtAndExpand(gr,ipmt);
}


//________________________________________________________________
/*
void AliT0CalibData::SetAmpLED(Int_t ipmt)
{
  Float_t mv, ps; 
  Float_t x[100], y[100];
  string buffer;
  
  const char * filename = gSystem->ExpandPathName("$ALICE_ROOT/T0/data/CFD-LED.txt");
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
  fAmpLED.AddAtAndExpand(gr,ipmt);
   
}
*/
//________________________________________________________________

void AliT0CalibData::SetAmpLEDRec(Int_t ipmt)
{
  Float_t mv, ps; 
  Float_t x[100], y[100];
  string buffer;
  
 const char * filename = gSystem->ExpandPathName("$ALICE_ROOT/T0/data/CFD-LED.txt");
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
  fAmpLEDRec.AddAtAndExpand(gr,ipmt);
  
}

//________________________________________________________________

void AliT0CalibData::ReadAsciiLookup(const Char_t *filename)
{
  Int_t key, trm, tdc, chain, channel;

  if(filename == 0){
    AliError(Form("Please, specify file with database")) ;
    return ;
  }


  ifstream lookup;
  lookup.open(filename);
  if(!lookup)
    {
     AliError(Form("!!!!!!!!!!!!!!No look up table in CDB!" ));
 
    }
  Char_t varname[11];
  Int_t ntrms;
  if(lookup)
    {
      lookup>>ntrms;
      cout<<" !!!!!!! ntrms "<<ntrms<<endl;
      //      fNumberOfTRMs=ntrms;
      SetNumberOfTRMs(ntrms);
       while(!lookup.eof())
	{
	  AliT0LookUpKey * lookkey= new AliT0LookUpKey();
	  AliT0LookUpValue * lookvalue= new AliT0LookUpValue();
	  
	  lookup>>varname>>key>>trm>>chain>>tdc>>channel;
	  lookvalue->SetTRM(trm);
	  lookvalue->SetTDC(tdc);
	  lookvalue->SetChain(chain);
	  lookvalue->SetChannel(channel);
	  lookkey->SetKey(key);
	  lookkey->SetChannelName(varname);
	  cout<<"lookup "<<varname<<" "<<key<<" "<<trm<<" "<<chain<<" "<<tdc<<" "<<channel<<endl;	  
	  
	  fLookup.Add((TObject*)lookvalue,(TObject*)lookkey);
	  
	}
      
      lookup.close();
      
    }
}
//________________________________________________________________

Int_t AliT0CalibData::GetChannel(Int_t trm,  Int_t tdc, Int_t chain, Int_t channel)
{

  AliT0LookUpKey * lookkey;//= new AliT0LookUpKey();
  AliT0LookUpValue * lookvalue= new AliT0LookUpValue(trm,tdc,chain,channel);

  lookkey = (AliT0LookUpKey*) fLookup.GetValue((TObject*)lookvalue);

  return lookkey->GetKey();

}

