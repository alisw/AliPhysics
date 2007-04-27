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
//  TPC calibration class for parameters which saved per pad                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"
#include <TObjArray.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TH2F.h>
#include "TTreeStream.h"

ClassImp(AliTPCCalPad)

//_____________________________________________________________________________
AliTPCCalPad::AliTPCCalPad():TNamed()
{
  //
  // AliTPCCalPad default constructor
  //

  for (Int_t isec = 0; isec < kNsec; isec++) {
    fROC[isec] = 0;
  }

}

//_____________________________________________________________________________
AliTPCCalPad::AliTPCCalPad(const Text_t *name, const Text_t *title)
                :TNamed(name,title)
{
  //
  // AliTPCCalPad constructor
  //
  for (Int_t isec = 0; isec < kNsec; isec++) {
    fROC[isec] = new AliTPCCalROC(isec);
  }
}


//_____________________________________________________________________________
AliTPCCalPad::AliTPCCalPad(const AliTPCCalPad &c):TNamed(c)
{
  //
  // AliTPCCalPad copy constructor
  //

  for (Int_t isec = 0; isec < kNsec; isec++) {
         fROC[isec] = 0;
     if (c.fROC[isec])
       fROC[isec] = new AliTPCCalROC(*(c.fROC[isec]));
  }
}

//_____________________________________________________________________________
AliTPCCalPad::AliTPCCalPad(TObjArray * array):TNamed()
{
  //
  // AliTPCCalPad default constructor
  //

  for (Int_t isec = 0; isec < kNsec; isec++) {
    fROC[isec] = (AliTPCCalROC *)array->At(isec);
  }

}


///_____________________________________________________________________________
AliTPCCalPad::~AliTPCCalPad()
{
  //
  // AliTPCCalPad destructor
  //

  for (Int_t isec = 0; isec < kNsec; isec++) {
    if (fROC[isec]) {
      delete fROC[isec];
      fROC[isec] = 0;
    }
  }

}

//_____________________________________________________________________________
AliTPCCalPad &AliTPCCalPad::operator=(const AliTPCCalPad &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTPCCalPad &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTPCCalPad::Copy(TObject &c) const
{
  //
  // Copy function
  //

  for (Int_t isec = 0; isec < kNsec; isec++) {
    if (fROC[isec]) {
      fROC[isec]->Copy(*((AliTPCCalPad &) c).fROC[isec]);
    }
  }
  TObject::Copy(c);
}

//_____________________________________________________________________________
void AliTPCCalPad::Add(Float_t c1)
{
    //
    // add constant for all channels of all ROCs
    //

    for (Int_t isec = 0; isec < kNsec; isec++) {
	if (fROC[isec]){
	    fROC[isec]->Add(c1);
	}
    }
}

//_____________________________________________________________________________
void AliTPCCalPad::Multiply(Float_t c1)
{
    //
    // multiply constant for all channels of all ROCs
    //
    for (Int_t isec = 0; isec < kNsec; isec++) {
	if (fROC[isec]){
	    fROC[isec]->Multiply(c1);
	}
    }
}

//_____________________________________________________________________________
void AliTPCCalPad::Add(const AliTPCCalPad * pad, Double_t c1)
{
    //
    // add calpad channel by channel multiplied by c1 - all ROCs
    //
    for (Int_t isec = 0; isec < kNsec; isec++) {
	if (fROC[isec]){
	    fROC[isec]->Add(pad->GetCalROC(isec),c1);
	}
    }
}

//_____________________________________________________________________________
void AliTPCCalPad::Multiply(const AliTPCCalPad * pad)
{
    //
    // multiply calpad channel by channel - all ROCs
    //
    for (Int_t isec = 0; isec < kNsec; isec++) {
	if (fROC[isec]){
	    fROC[isec]->Multiply(pad->GetCalROC(isec));
	}
    }
}

//_____________________________________________________________________________
void AliTPCCalPad::Divide(const AliTPCCalPad * pad)
{
    //
    // divide calpad channel by channel - all ROCs
    //
    for (Int_t isec = 0; isec < kNsec; isec++) {
	if (fROC[isec]){
	    fROC[isec]->Divide(pad->GetCalROC(isec));
	}
    }
}

//_____________________________________________________________________________
TGraph  *  AliTPCCalPad::MakeGraph(Int_t type, Float_t ratio){
  //
  //   type=1 - mean
  //        2 - median
  //        3 - LTM
  Int_t npoints = 0;
  for (Int_t i=0;i<72;i++) if (fROC[i]) npoints++;
  TGraph * graph = new TGraph(npoints);
  npoints=0;   
  for (Int_t isec=0;isec<72;isec++){
    if (!fROC[isec]) continue;
    if (type==0)  graph->SetPoint(npoints,isec,fROC[isec]->GetMean());      
    if (type==1)  graph->SetPoint(npoints,isec,fROC[isec]->GetMedian());
    if (type==2)  graph->SetPoint(npoints,isec,fROC[isec]->GetLTM(0,ratio));    
    npoints++;
  }

  graph->GetXaxis()->SetTitle("Sector"); 
  if (type==0) {
    graph->GetYaxis()->SetTitle("Mean");   
    graph->SetMarkerStyle(22);    
  }
  if (type==1) {
    graph->GetYaxis()->SetTitle("Median");   
    graph->SetMarkerStyle(22);    
  }
  if (type==2) {
      graph->GetYaxis()->SetTitle(Form("Mean%f",ratio));      
      graph->SetMarkerStyle(24);
  }

  return graph;
}

//_____________________________________________________________________________
Double_t AliTPCCalPad::GetMeanRMS(Double_t &rms)
{
    //
    // Calculate mean an RMS of all rocs
    //
    Double_t sum = 0, sum2 = 0, n=0, val=0;
    for (Int_t isec = 0; isec < kNsec; isec++) {
        AliTPCCalROC *calRoc = fROC[isec];
	if ( calRoc ){
	    for (UInt_t irow=0; irow<calRoc->GetNrows(); irow++){
		for (UInt_t ipad=0; ipad<calRoc->GetNPads(irow); ipad++){
		    val = calRoc->GetValue(irow,ipad);
		    sum+=val;
		    sum2+=val*val;
                    n++;
		}
	    }

	}
    }
    Double_t n1 = 1./n;
    Double_t mean = sum*n1;
    rms  = TMath::Sqrt(TMath::Abs(sum2*n1-mean*mean));
    return mean;
}


//_____________________________________________________________________________
Double_t AliTPCCalPad::GetMean()
{
    //
    // return mean of the mean of all ROCs
    //
    Double_t arr[kNsec];
    Int_t n=0;
    for (Int_t isec = 0; isec < kNsec; isec++) {
        AliTPCCalROC *calRoc = fROC[isec];
	if ( calRoc ){
	    arr[n] = calRoc->GetMean();
            n++;
	}
    }
    return TMath::Mean(n,arr);
}

//_____________________________________________________________________________
Double_t AliTPCCalPad::GetRMS()
{
    //
    // return mean of the RMS of all ROCs
    //
    Double_t arr[kNsec];
    Int_t n=0;
    for (Int_t isec = 0; isec < kNsec; isec++) {
        AliTPCCalROC *calRoc = fROC[isec];
	if ( calRoc ){
	    arr[n] = calRoc->GetRMS();
            n++;
	}
    }
    return TMath::Mean(n,arr);
}

//_____________________________________________________________________________
Double_t AliTPCCalPad::GetMedian()
{
    //
    // return mean of the median of all ROCs
    //
    Double_t arr[kNsec];
    Int_t n=0;
    for (Int_t isec = 0; isec < kNsec; isec++) {
        AliTPCCalROC *calRoc = fROC[isec];
	if ( calRoc ){
	    arr[n] = calRoc->GetMedian();
            n++;
	}
    }
    return TMath::Mean(n,arr);
}

//_____________________________________________________________________________
Double_t AliTPCCalPad::GetLTM(Double_t *sigma, Double_t fraction)
{
    //
    // return mean of the LTM and sigma of all ROCs
    //
    Double_t arrm[kNsec];
    Double_t arrs[kNsec];
    Double_t *sTemp=0x0;
    Int_t n=0;

    for (Int_t isec = 0; isec < kNsec; isec++) {
        AliTPCCalROC *calRoc = fROC[isec];
	if ( calRoc ){
	    if ( sigma ) sTemp=arrs+n;
	    arrm[n] = calRoc->GetLTM(sTemp,fraction);
            n++;
	}
    }
    if ( sigma ) *sigma = TMath::Mean(n,arrs);
    return TMath::Mean(n,arrm);
}

//_____________________________________________________________________________
TH1F * AliTPCCalPad::MakeHisto1D(Float_t min, Float_t max,Int_t type){
  //
  // make 1D histo
  // type -1 = user defined range
  //       0 = nsigma cut nsigma=min
  if (type>=0){
    if (type==0){
      // nsigma range
      Float_t mean  = GetMean();
      Float_t sigma = GetRMS();
      Float_t nsigma = TMath::Abs(min);
      min = mean-nsigma*sigma;
      max = mean+nsigma*sigma;
    }
    if (type==1){
      // fixed range
      Float_t mean   = GetMedian();
      Float_t  delta = min;
      min = mean-delta;
      max = mean+delta;
    }
    if (type==2){
      //
      // LTM mean +- nsigma
      //
      Double_t sigma;
      Float_t mean  = GetLTM(&sigma,max);
      sigma*=min;
      min = mean-sigma;
      max = mean+sigma;
    }
  }
  char  name[1000];
  sprintf(name,"%s Pad 1D",GetTitle());
  TH1F * his = new TH1F(name,name,100, min,max);
    for (Int_t isec = 0; isec < kNsec; isec++) {
	if (fROC[isec]){
	    for (UInt_t irow=0; irow<fROC[isec]->GetNrows(); irow++){
		UInt_t npads = (Int_t)fROC[isec]->GetNPads(irow);
		for (UInt_t ipad=0; ipad<npads; ipad++){
		    his->Fill(fROC[isec]->GetValue(irow,ipad));
		}
	    }
	}
    }
  return his;
}

//_____________________________________________________________________________
TH2F *AliTPCCalPad::MakeHisto2D(Int_t side){
  //
  // Make 2D graph
  // side  -  specify the side A = 0 C = 1
  // type  -  used types of determination of boundaries in z
  Float_t kEpsilon = 0.000000000001;
  TH2F * his = new TH2F(GetName(), GetName(), 250,-250,250,250,-250,250);
  AliTPCROC * roc  = AliTPCROC::Instance(); 
  for (Int_t isec=0; isec<72; isec++){
    if (side==0 && isec%36>=18) continue;
    if (side>0 && isec%36<18) continue;
    if (fROC[isec]){
      AliTPCCalROC * calRoc = fROC[isec];
      for (UInt_t irow=0; irow<calRoc->GetNrows(); irow++)
	for (UInt_t ipad=0; ipad<calRoc->GetNPads(irow); ipad++)
	  if (TMath::Abs(calRoc->GetValue(irow,ipad))>kEpsilon){
	    Float_t xyz[3];
	    roc->GetPositionGlobal(isec,irow,ipad,xyz);
	    Int_t binx = 1+TMath::Nint((xyz[0]+250.)*0.5);
	    Int_t biny = 1+TMath::Nint((xyz[1]+250.)*0.5);
	    Float_t value = calRoc->GetValue(irow,ipad);	    
	    his->SetBinContent(binx,biny,value);
	  }
    }
  }
  his->SetXTitle("x (cm)");
  his->SetYTitle("y (cm)");
  return his;
}




void AliTPCCalPad::MakeTree(const char * fileName, TObjArray * array) {
  //
  // Write tree with all available information
  //
   TTreeSRedirector cstream(fileName);
   AliTPCROC* tpcROCinstance = AliTPCROC::Instance();
   Int_t arrayEntries = array->GetEntries();
   
   TString* names = new TString[arrayEntries];
   for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++)
      names[ivalue].Append(((AliTPCCalPad*)array->At(ivalue))->GetName());

   for (UInt_t isector = 0; isector < tpcROCinstance->GetNSectors(); isector++) {
      //
      // get statistic for given sector
      //
      TVectorF median(arrayEntries);
      TVectorF mean(arrayEntries);
      TVectorF rms(arrayEntries);
      TVectorF ltm(arrayEntries);
      
      TVectorF *vectorArray = new TVectorF[arrayEntries];
      for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++)
         vectorArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(isector));
      
      for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
         AliTPCCalPad* calPad = (AliTPCCalPad*) array->At(ivalue);
         AliTPCCalROC* calROC = calPad->GetCalROC(isector);
         if (calROC) {
            median[ivalue] = calROC->GetMedian();
            mean[ivalue] = calROC->GetMean();
            rms[ivalue] = calROC->GetRMS();
            ltm[ivalue] = calROC->GetLTM();
         }
         else {
            median[ivalue] = 0.;
            mean[ivalue] = 0.;
            rms[ivalue] = 0.;
            ltm[ivalue] = 0.;
         }
      }
      
      //
      // fill vectors of variable per pad
      //
      TVectorF *posArray = new TVectorF[6];
      for (Int_t ivalue = 0; ivalue < 6; ivalue++)
         posArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(isector));

      Float_t posG[3] = {0};
      Float_t posL[3] = {0};
      Int_t ichannel = 0;
      for (UInt_t irow = 0; irow < tpcROCinstance->GetNRows(isector); irow++) {
         for (UInt_t ipad = 0; ipad < tpcROCinstance->GetNPads(isector, irow); ipad++) {
            tpcROCinstance->GetPositionLocal(isector, irow, ipad, posL);
            tpcROCinstance->GetPositionGlobal(isector, irow, ipad, posG);
            posArray[0][ichannel] = irow;
            posArray[1][ichannel] = ipad;
            posArray[2][ichannel] = posL[0];
            posArray[3][ichannel] = posL[1];
            posArray[4][ichannel] = posG[0];
            posArray[5][ichannel] = posG[1];
            
            // loop over array containing AliTPCCalPads
            for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
               AliTPCCalPad* calPad = (AliTPCCalPad*) array->At(ivalue);
               AliTPCCalROC* calROC = calPad->GetCalROC(isector);
               if (calROC)
                  (vectorArray[ivalue])[ichannel] = calROC->GetValue(irow, ipad);
               else
                  (vectorArray[ivalue])[ichannel] = 0;
            }
            ichannel++;
         }
      }
      
      cstream << "calPads" <<
         "sector=" << isector;
      
      for  (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
         cstream << "calPads" <<
            (Char_t*)((names[ivalue] + "_Median=").Data()) << median[ivalue] <<
            (Char_t*)((names[ivalue] + "_Mean=").Data()) << mean[ivalue] <<
            (Char_t*)((names[ivalue] + "_RMS=").Data()) << rms[ivalue] <<
            (Char_t*)((names[ivalue] + "_LTM=").Data()) << ltm[ivalue];
      }

      for  (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
         cstream << "calPads" <<
            (Char_t*)((names[ivalue] + ".=").Data()) << &vectorArray[ivalue];
      }

      cstream << "calPads" <<
         "row.=" << &posArray[0] <<
         "pad.=" << &posArray[1] <<
         "lx.=" << &posArray[2] <<
         "ly.=" << &posArray[3] <<
         "gx.=" << &posArray[4] <<
         "gy.=" << &posArray[5];
         
      cstream << "calPads" <<
         "\n";

      delete[] posArray;
      delete[] vectorArray;
   }
   delete[] names;
}


