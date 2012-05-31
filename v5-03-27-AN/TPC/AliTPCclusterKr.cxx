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





/* $Id: AliTPCclusterKr.cxx,v 1.7 2008/01/22 17:24:53 matyja Exp $ */





//-----------------------------------------------------------------


//           Implementation of the TPC Kr cluster class


//


// Origin: Adam Matyja, INP PAN, adam.matyja@ifj.edu.pl


//-----------------------------------------------------------------





#include "AliTPCclusterKr.h"


#include "AliCluster.h"


#include "AliTPCvtpr.h"


#include "TObjArray.h"


//#include "TH1F.h"


#include "TMath.h"


#include "TArrayI.h"





ClassImp(AliTPCclusterKr)








AliTPCclusterKr::AliTPCclusterKr()


:AliCluster(),


 fMax(),


 fADCcluster(0),


 fSec(0),


 fNPads(0),


 fNRows(0),


 fTimebins1D(0),


 fPads1D(0),


 fPadRMS(0),


 fRowRMS(0),


 fTimebinRMS(0),


 fSize(0),


 fCenterX(0),


 fCenterY(0),


 fCenterT(0),


 fCluster(0),


 fTimeStamp(0),
 fRun(0)


{


//


// default constructor


//


  fCluster=new TObjArray();


}





AliTPCclusterKr::AliTPCclusterKr(const AliTPCclusterKr &param)


:AliCluster(param),


 fMax(),


 fADCcluster(0),


 fSec(0),


 fNPads(0),


 fNRows(0),


 fTimebins1D(0),


 fPads1D(0),


 fPadRMS(0),


 fRowRMS(0),


 fTimebinRMS(0),


 fSize(0),


 fCenterX(0),


 fCenterY(0),


 fCenterT(0),


 fCluster(0),


 fTimeStamp(0),
 fRun(0)


{


//


// copy constructor


//


  fADCcluster = param.fADCcluster;


  fSec  = param.fSec ;


  fNPads = param.fNPads;


  fNRows = param.fNRows;


  fMax = param.fMax;


  //  fCluster = param.fCluster;


  fCenterX = param.fCenterX;


  fCenterY = param.fCenterY;


  fCenterT = param.fCenterT;


  fCluster=new TObjArray(*(param.fCluster));


  fSize = param.fSize;


  fTimebins1D = param.fTimebins1D;


  fPads1D = param.fPads1D;


  fPadRMS = param.fPadRMS;


  fRowRMS = param.fRowRMS;


  fTimebinRMS = param.fTimebinRMS;


  fTimeStamp = param.fTimeStamp;
  fRun = param.fRun;


} 





AliTPCclusterKr &AliTPCclusterKr::operator = (const AliTPCclusterKr & param)


{


  //


  // assignment operator


  // 
  if (this == &param) return (*this);

  (AliCluster&)(*this) = (AliCluster&)param;


  fADCcluster = param.fADCcluster;


  fSec  = param.fSec ;


  fNPads = param.fNPads;


  fNRows = param.fNRows;


  fMax = param.fMax;


  //  fCluster=param.fCluster;


  fCenterX = param.fCenterX;


  fCenterY = param.fCenterY;


  fCenterT = param.fCenterT;


  delete fCluster;


  fCluster=new TObjArray(*(param.fCluster));


  fSize=param.fSize;


  fTimebins1D = param.fTimebins1D;


  fPads1D = param.fPads1D;


  fPadRMS = param.fPadRMS;


  fRowRMS = param.fRowRMS;


  fTimebinRMS = param.fTimebinRMS;


  fTimeStamp = param.fTimeStamp;
  fRun = param.fRun;


  return (*this);


}





AliTPCclusterKr::~AliTPCclusterKr()


{


  //


  // destructor


  //


  if(fCluster) {


    fCluster->SetOwner(kTRUE);


    fCluster->Delete();


    delete fCluster;


  }


  fCluster=0;


}





////____________________________________________________________________________


void AliTPCclusterKr::SetCenter(){


  //


  // calculate geometrical center of the cluster


  //


  Double_t rX=0;


  Double_t rY=0;


  Double_t rT=0;





  Short_t adc;


  fADCcluster=0;


  for(Int_t iter = 0; iter < fCluster->GetEntriesFast(); ++iter) {


    AliTPCvtpr *iclus=(AliTPCvtpr *)fCluster->At(iter);





    //for( std::vector<AliTPCvtpr*>::iterator iclus  = fCluster.begin();


    //iclus != fCluster.end(); ++iclus ) {


    adc = (iclus)->GetAdc();


    fADCcluster+=adc;


    rX += ((iclus)->GetX() * adc);


    rY += ((iclus)->GetY() * adc);


    rT += ((iclus)->GetT() * adc);


  }
  if(fADCcluster){ 

    fCenterX=rX/fADCcluster;


    fCenterY=rY/fADCcluster;


    fCenterT=rT/fADCcluster;
  }





  return;


}





void AliTPCclusterKr::SetPadRMS(){


  //


  // calculate RMS in pad direction


  //


  //  TH1F *histo= new TH1F("","",200,0,200);


  TArrayI *array= new TArrayI(fCluster->GetEntriesFast());


  for(Int_t i=0;i<fCluster->GetEntriesFast();i++)


    {


      array->SetAt(((AliTPCvtpr *)(fCluster->At(i)))->GetPad(),i);


      //histo->Fill( ((AliTPCvtpr *)(fCluster->At(i)))->GetPad() );


    }


  //  fPadRMS=histo->GetRMS();


  fPadRMS=TMath::RMS(array->GetSize(),array->GetArray());


  //  delete histo;


  delete array;


  return;


}





void AliTPCclusterKr::SetRowRMS(){


  //


  // calculate RMS in row direction


  //


  TArrayI *array= new TArrayI(fCluster->GetEntriesFast());


  //  TH1F *histo= new TH1F("","",120,0,120);


  for(Int_t i=0;i<fCluster->GetEntriesFast();i++)


    {


      array->SetAt(((AliTPCvtpr *)(fCluster->At(i)))->GetRow(),i);


      //      histo->Fill( ((AliTPCvtpr *)(fCluster->At(i)))->GetRow() );


    }


  //  fRowRMS=histo->GetRMS();


  fRowRMS=TMath::RMS(array->GetSize(),array->GetArray());


  //  delete histo;


  delete array;


  return;


}





void AliTPCclusterKr::SetTimebinRMS(){


  //


  // calculate RMS in timebin direction


  //


  TArrayI *array= new TArrayI(fCluster->GetEntriesFast());


  //  TH1F *histo= new TH1F("","",1000,0,1000);


  for(Int_t i=0;i<fCluster->GetEntriesFast();i++)


    {


      array->SetAt(((AliTPCvtpr *)(fCluster->At(i)))->GetTime(),i);


      //      histo->Fill( ((AliTPCvtpr *)(fCluster->At(i)))->GetTime() );


    }


  fTimebinRMS=TMath::RMS(array->GetSize(),array->GetArray());


  //histo->GetRMS();


  //  delete histo;


  delete array;


  return;


}





void AliTPCclusterKr::SetRMS(){


  //


  // calculate RMS in pad,row,timebin direction


  //


  TArrayI *arrayPad = new TArrayI(fCluster->GetEntriesFast());


  TArrayI *arrayRow = new TArrayI(fCluster->GetEntriesFast());


  TArrayI *arrayTime= new TArrayI(fCluster->GetEntriesFast());


  //  TH1F *histoPad= new TH1F("p","p",200,0,200);


  //  TH1F *histoRow= new TH1F("r","r",120,0,120);


  //  TH1F *histoTime= new TH1F("t","t",1000,0,1000);


  for(Int_t i=0;i<fCluster->GetEntriesFast();i++)


    {


      arrayPad->SetAt(((AliTPCvtpr *)(fCluster->At(i)))->GetPad(),i);


      arrayRow->SetAt(((AliTPCvtpr *)(fCluster->At(i)))->GetRow(),i);


      arrayTime->SetAt(((AliTPCvtpr *)(fCluster->At(i)))->GetTime(),i);





      //histoPad->Fill( ((AliTPCvtpr *)(fCluster->At(i)))->GetPad() );


      //histoRow->Fill( ((AliTPCvtpr *)(fCluster->At(i)))->GetRow() );


      //histoTime->Fill( ((AliTPCvtpr *)(fCluster->At(i)))->GetTime() );


    }


  //  fPadRMS=histoPad->GetRMS();


  fPadRMS=TMath::RMS(arrayPad->GetSize(),arrayPad->GetArray());


  fRowRMS=TMath::RMS(arrayRow->GetSize(),arrayRow->GetArray());


    //histoRow->GetRMS();


  fTimebinRMS=TMath::RMS(arrayTime->GetSize(),arrayTime->GetArray());


    //histoTime->GetRMS();





  delete arrayPad;


  delete arrayRow;


  delete arrayTime;


  //  delete histoPad;


  //  delete histoRow;


  //  delete histoTime;





  return;


}








void AliTPCclusterKr::Set1D(){


  //


  //


  //


  Short_t maxTime=0;


  Short_t minTime=1000;


  Short_t maxPad=0;


  Short_t minPad=1000;


 


  for(Int_t i=0;i<fCluster->GetEntriesFast();i++)


    {


      if(((AliTPCvtpr *)(fCluster->At(i)))->GetPad()>maxPad)maxPad   =((AliTPCvtpr *)(fCluster->At(i)))->GetPad();


      if(((AliTPCvtpr *)(fCluster->At(i)))->GetPad()<minPad)minPad   =((AliTPCvtpr *)(fCluster->At(i)))->GetPad();


      if(((AliTPCvtpr *)(fCluster->At(i)))->GetTime()>maxTime)maxTime=((AliTPCvtpr *)(fCluster->At(i)))->GetTime();


      if(((AliTPCvtpr *)(fCluster->At(i)))->GetTime()<minTime)minTime=((AliTPCvtpr *)(fCluster->At(i)))->GetTime();


    }


  fPads1D=maxPad-minPad+1;


  fTimebins1D=maxTime-minTime+1;


  return;


}


