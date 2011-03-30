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
// class for T0 calibration                       TM--AM_6-02-2006         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliT0CalibData.h"
#include "AliT0LookUpValue.h"
#include "AliT0LookUpKey.h"
#include "AliLog.h"

#include <Riostream.h>

//#include <string>

ClassImp(AliT0CalibData)

//________________________________________________________________
  AliT0CalibData::AliT0CalibData():   TNamed(),
				      fLookup(0),
				      fNumberOfTRMs(0)

{
  //
}

//________________________________________________________________
AliT0CalibData::AliT0CalibData(const char* name):TNamed(),
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
void  AliT0CalibData::PrintLookup(Option_t*, Int_t iTRM, Int_t iTDC, Int_t iChannel) const
{
  // print lookup table

  AliT0LookUpKey* lookkey; //= new AliT0LookUpKey();
  AliT0LookUpValue*  lookvalue= new AliT0LookUpValue();
  printf("Number Of TRMs in setup %i\n",GetNumberOfTRMs());

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
    //    TString name= lookkey->GetChannelName();
    // cout<<name.Data()<<endl;
    if (lookkey)
      {
	TString name= lookkey->GetChannelName();
	/*	cout<<" lookup KEY!!! "<<name.Data()<<" "<<lookkey->GetKey()<<" VALUE "<<lookvalue->GetTRM()<<" "
	    <<lookvalue->GetTDC()<<" "
	    << lookvalue->GetChain()<<" "
	    <<lookvalue->GetChannel()<<endl;*/
      }
  }
  
}
//________________________________________________________________

void AliT0CalibData::ReadAsciiLookup(const Char_t *filename)
{
  // read lookup table from ascii file

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
	  
	  fLookup.Add((TObject*)lookvalue,(TObject*)lookkey);
	  
	}
      
      lookup.close();
      
    }
}
//________________________________________________________________

Int_t AliT0CalibData::GetChannel(Int_t trm,  Int_t tdc, Int_t chain, Int_t channel)
{
  // read number of channel according physical addres 


  AliT0LookUpKey * lookkey;//= new AliT0LookUpKey();
  AliT0LookUpValue * lookvalue= new AliT0LookUpValue(trm,tdc,chain,channel);

  lookkey = (AliT0LookUpKey*) fLookup.GetValue((TObject*)lookvalue);

  return lookkey->GetKey();

}

