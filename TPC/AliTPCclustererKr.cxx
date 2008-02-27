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

/* $Id: AliTPCclustererKr.cxx,v 1.7 2008/02/06 17:24:53 matyja Exp $ */

//-----------------------------------------------------------------
//           Implementation of the TPC cluster class
//
// Origin: Adam Matyja, INP PAN, adam.matyja@ifj.edu.pl
//-----------------------------------------------------------------

#include "AliTPCclustererKr.h"
#include "AliTPCclusterKr.h"
#include <vector>
#include <list>
#include "TObject.h"
#include "AliPadMax.h"
#include "AliSimDigits.h"
#include "AliTPCv4.h"
#include "AliTPCParam.h"
#include "AliTPCDigitsArray.h"
#include "AliTPCvtpr.h"
#include "AliTPCClustersRow.h"
#include "TTree.h"

//used in raw data finder
#include "AliTPCROC.h"
#include "AliTPCCalPad.h"
#include "AliTPCAltroMapping.h"
#include "AliTPCcalibDB.h"
#include "AliTPCRawStream.h"
#include "AliTPCRecoParam.h"
#include "AliTPCReconstructor.h"
#include "AliRawReader.h"
#include "AliTPCCalROC.h"

ClassImp(AliTPCclustererKr)


AliTPCclustererKr::AliTPCclustererKr()
  :TObject(),
  fRawData(kFALSE),
  fRowCl(0),
  fInput(0),
  fOutput(0),
  fParam(0),
  fDigarr(0),
  fRecoParam(0),
  fIsOldRCUFormat(kFALSE),
  fZeroSup(2),
  fFirstBin(60),
  fLastBin(950),
  fMaxNoiseAbs(2),
  fMaxNoiseSigma(3),
  fMinAdc(3),
  fMinTimeBins(2),
//  fMaxPadRange(4),
//  fMaxRowRange(3),
  fMaxTimeRange(7),
  fValueToSize(3.5),
  fMaxPadRangeCm(2.5),
  fMaxRowRangeCm(3.5)
{
//
// default constructor
//
}

AliTPCclustererKr::AliTPCclustererKr(const AliTPCclustererKr &param)
  :TObject(),
  fRawData(kFALSE),
  fRowCl(0),
  fInput(0),
  fOutput(0),
  fParam(0),
  fDigarr(0),
  fRecoParam(0),
  fIsOldRCUFormat(kFALSE),
  fZeroSup(2),
  fFirstBin(60),
  fLastBin(950),
  fMaxNoiseAbs(2),
  fMaxNoiseSigma(3),
  fMinAdc(3),
  fMinTimeBins(2),
//  fMaxPadRange(4),
//  fMaxRowRange(3),
  fMaxTimeRange(7),
  fValueToSize(3.5),
  fMaxPadRangeCm(2.5),
  fMaxRowRangeCm(3.5)
{
//
// copy constructor
//
  fParam = param.fParam;
  fRecoParam = param.fRecoParam;
  fIsOldRCUFormat = param.fIsOldRCUFormat;
  fRawData = param.fRawData;
  fRowCl  = param.fRowCl ;
  fInput  = param.fInput ;
  fOutput = param.fOutput;
  fDigarr = param.fDigarr;
  fZeroSup       = param.fZeroSup       ;
  fFirstBin	 = param.fFirstBin	;
  fLastBin	 = param.fLastBin	;
  fMaxNoiseAbs	 = param.fMaxNoiseAbs	;
  fMaxNoiseSigma = param.fMaxNoiseSigma ;
  fMinAdc = param.fMinAdc;
  fMinTimeBins = param.fMinTimeBins;
//  fMaxPadRange  = param.fMaxPadRange ;
//  fMaxRowRange  = param.fMaxRowRange ;
  fMaxTimeRange = param.fMaxTimeRange;
  fValueToSize  = param.fValueToSize;
  fMaxPadRangeCm = param.fMaxPadRangeCm;
  fMaxRowRangeCm = param.fMaxRowRangeCm;
} 

AliTPCclustererKr & AliTPCclustererKr::operator = (const AliTPCclustererKr & param)
{
  fParam = param.fParam;
  fRecoParam = param.fRecoParam;
  fIsOldRCUFormat = param.fIsOldRCUFormat;
  fRawData = param.fRawData;
  fRowCl  = param.fRowCl ;
  fInput  = param.fInput ;
  fOutput = param.fOutput;
  fDigarr = param.fDigarr;
  fZeroSup       = param.fZeroSup       ;
  fFirstBin	 = param.fFirstBin	;
  fLastBin	 = param.fLastBin	;
  fMaxNoiseAbs	 = param.fMaxNoiseAbs	;
  fMaxNoiseSigma = param.fMaxNoiseSigma ;
  fMinAdc = param.fMinAdc;
  fMinTimeBins = param.fMinTimeBins;
//  fMaxPadRange  = param.fMaxPadRange ;
//  fMaxRowRange  = param.fMaxRowRange ;
  fMaxTimeRange = param.fMaxTimeRange;
  fValueToSize  = param.fValueToSize;
  fMaxPadRangeCm = param.fMaxPadRangeCm;
  fMaxRowRangeCm = param.fMaxRowRangeCm;
  return (*this);
}

AliTPCclustererKr::~AliTPCclustererKr()
{
  //
  // destructor
  //
}

void AliTPCclustererKr::SetRecoParam(AliTPCRecoParam *recoParam)
{
  if (recoParam) {
    fRecoParam = recoParam;
  }else{
    //set default parameters if not specified
    fRecoParam = AliTPCReconstructor::GetRecoParam();
    if (!fRecoParam)  fRecoParam = AliTPCRecoParam::GetLowFluxParam();
  }
  return;
}


////____________________________________________________________________________
////       I/O
void AliTPCclustererKr::SetInput(TTree * tree)
{
  //
  // set input tree with digits
  //
  fInput = tree;  
  if  (!fInput->GetBranch("Segment")){
    cerr<<"AliTPCclusterKr::FindClusterKr(): no proper input tree !\n";
    fInput=0;
    return;
  }
}

void AliTPCclustererKr::SetOutput(TTree * tree) 
{
  //
  // set output
  //
  fOutput= tree;  
  AliTPCClustersRow clrow;
  AliTPCClustersRow *pclrow=&clrow;  
  clrow.SetClass("AliTPCclusterKr");
  clrow.SetArray(1); // to make Clones array
  fOutput->Branch("Segment","AliTPCClustersRow",&pclrow,32000,200);    
}

////____________________________________________________________________________
//// with new I/O
Int_t AliTPCclustererKr::FinderIO()
{
  // Krypton cluster finder for simulated events from MC

  if (!fInput) { 
    Error("Digits2Clusters", "input tree not initialised");
    return 10;
  }
  
  if (!fOutput) {
    Error("Digits2Clusters", "output tree not initialised");
    return 11;
  }

  FindClusterKrIO();
  return 0;
}

Int_t AliTPCclustererKr::FinderIO(AliRawReader* rawReader)
{
  // Krypton cluster finder for the TPC raw data
  //
  // fParam must be defined before

  // consider noiceROC or not

  if(rawReader)fRawData=kTRUE; //set flag to data

  if (!fOutput) {
    Error("Digits2Clusters", "output tree not initialised");
    return 11;
  }

  fParam->SetMaxTBin(fRecoParam->GetLastBin());//set number of timebins from reco -> param
  //   used later for memory allocation

  Bool_t isAltro=kFALSE;

  AliTPCROC * roc = AliTPCROC::Instance();
  AliTPCCalPad * noiseTPC = AliTPCcalibDB::Instance()->GetPadNoise();
  AliTPCAltroMapping** mapping =AliTPCcalibDB::Instance()->GetMapping();
  //
  AliTPCRawStream input(rawReader,(AliAltroMapping**)mapping);

  const Int_t kNIS = fParam->GetNInnerSector();//number of inner sectors
  const Int_t kNOS = fParam->GetNOuterSector();//number of outer sectors
  const Int_t kNS = kNIS + kNOS;//all sectors

//  //set cluster finder parameters
//  SetZeroSup(fParam->GetZeroSup());//zero suppression parameter
//  SetFirstBin(60);
//  SetLastBin(950);
//  SetMaxNoiseAbs(2);
//  SetMaxNoiseSigma(3);

  //crate TPC view
  AliTPCDigitsArray *digarr=new AliTPCDigitsArray(kFALSE);//data not sim
  digarr->Setup(fParam);//as usually parameters

  //
  // Loop over sectors
  //
  for(Int_t iSec = 0; iSec < kNS; iSec++) {
    AliTPCCalROC * noiseROC;
    if(noiseTPC==0x0){
      noiseROC =new AliTPCCalROC(iSec);//noise=0
    }
    else{
      noiseROC = noiseTPC->GetCalROC(iSec);  // noise per given sector
    }
    Int_t nRows = 0; //number of rows in sector
    Int_t nDDLs = 0; //number of DDLs
    Int_t indexDDL = 0; //DDL index
    if (iSec < kNIS) {
      nRows = fParam->GetNRowLow();
      nDDLs = 2;
      indexDDL = iSec * 2;
    }else {
      nRows = fParam->GetNRowUp();
      nDDLs = 4;
      indexDDL = (iSec-kNIS) * 4 + kNIS * 2;
    }

    //
    // Load the raw data for corresponding DDLs
    //
    rawReader->Reset();
    input.SetOldRCUFormat(fIsOldRCUFormat);
    rawReader->Select("TPC",indexDDL,indexDDL+nDDLs-1);

    if(input.Next()) {
      isAltro=kTRUE;
      // Allocate memory for rows in sector (pads(depends on row) x timebins)
      for(Int_t iRow = 0; iRow < nRows; iRow++) {
	digarr->CreateRow(iSec,iRow);
      }//end loop over rows
    }
    rawReader->Reset();
    input.SetOldRCUFormat(fIsOldRCUFormat);
    rawReader->Select("TPC",indexDDL,indexDDL+nDDLs-1);

    //
    // Begin loop over altro data
    //
    while (input.Next()) {

      //check sector consistency
      if (input.GetSector() != iSec)
	AliFatal(Form("Sector index mismatch ! Expected (%d), but got (%d) !",iSec,input.GetSector()));
      
      //check row consistency
      Short_t iRow = input.GetRow();
      if (iRow < 0 || iRow >= nRows){
	AliError(Form("Pad-row index (%d) outside the range (%d -> %d) !",
		      iRow, 0, nRows -1));
	continue;
      }

      //check pad consistency
      Short_t iPad = input.GetPad();
      if (iPad < 0 || iPad >= (Short_t)(roc->GetNPads(iSec,iRow))) {
	AliError(Form("Pad index (%d) outside the range (%d -> %d) !",
		      iPad, 0, roc->GetNPads(iSec,iRow) ));
	continue;
      }

      //check time consistency
      Short_t iTimeBin = input.GetTime();
      if ( iTimeBin < fRecoParam->GetFirstBin() || iTimeBin >= fRecoParam->GetLastBin()){
	//cout<<iTimeBin<<endl;
	continue;
	AliFatal(Form("Timebin index (%d) outside the range (%d -> %d) !",
		      iTimeBin, 0, fRecoParam->GetLastBin() -1));
      }

      //signal
      Int_t signal = input.GetSignal();
      if (signal <= fZeroSup || 
	  iTimeBin < fFirstBin ||
	  iTimeBin > fLastBin
	  ) {
	digarr->GetRow(iSec,iRow)->SetDigitFast(0,iTimeBin,iPad);
	continue;
      }
      if (!noiseROC) continue;
      Double_t noiseOnPad = noiseROC->GetValue(iRow,iPad);//noise on given pad and row in sector
      if (noiseOnPad > fMaxNoiseAbs) continue; // consider noisy pad as dead
      
      if(signal <= fMaxNoiseSigma * noiseOnPad){
	digarr->GetRow(iSec,iRow)->SetDigitFast(0,iTimeBin,iPad);
	continue;
      }

      digarr->GetRow(iSec,iRow)->SetDigitFast(signal,iTimeBin,iPad);

    }//end of loop over altro data
  }//end of loop over sectors
  
  SetDigArr(digarr);
  if(isAltro) FindClusterKrIO();
  delete digarr;

  return 0;
}

////____________________________________________________________________________
Int_t AliTPCclustererKr::FindClusterKrIO()
{
  //fParam and  fDigarr must be set to run this method

//  //set parameters 
//  SetMinAdc(3);//usually is 3 
//  SetMinTimeBins(2);//should be 2 - the best result of shape in MC
////  SetMaxPadRange(4);
////  SetMaxRowRange(3);
//  SetMaxTimeRange(7);
//  SetValueToSize(3.5);//3.5
//  SetMaxPadRangeCm(2.5);
//  SetMaxRowRangeCm(3.5);

  Int_t clusterCounter=0;
  const Short_t nTotalSector=fParam->GetNSector();//number of sectors
  for(Short_t iSec=0; iSec<nTotalSector; ++iSec){
    
    //vector of maxima for each sector
    std::vector<AliPadMax*> maximaInSector;
    
    //
    //  looking for the maxima on the pad
    //

    const Short_t nRows=fParam->GetNRow(iSec);//number of rows in sector
    for(Short_t iRow=0; iRow<nRows; ++iRow){
      AliSimDigits *digrow;
      if(fRawData){
	digrow = (AliSimDigits*)fDigarr->GetRow(iSec,iRow);//real data
      }else{
	digrow = (AliSimDigits*)fDigarr->LoadRow(iSec,iRow);//MC
      }
      if(digrow){//if pointer exist
	digrow->ExpandBuffer(); //decrunch
	const Short_t nPads = digrow->GetNCols();  // number of pads
	const Short_t nTime = digrow->GetNRows(); // number of timebins
	for(Short_t iPad=0;iPad<nPads;iPad++){
	  
	  Short_t timeBinMax=-1;//timebin of maximum 
	  Short_t valueMaximum=-1;//value of maximum in adc
	  Short_t increaseBegin=-1;//timebin when increase starts
	  Short_t sumAdc=0;//sum of adc on the pad in maximum surrounding
	  bool ifIncreaseBegin=true;//flag - check if increasing started
	  bool ifMaximum=false;//flag - check if it could be maximum

	  for(Short_t iTimeBin=0;iTimeBin<nTime;iTimeBin++){
	    Short_t adc = digrow->GetDigitFast(iTimeBin,iPad);
	    if(adc<fMinAdc){//standard was 3
	      if(ifMaximum){
		if(iTimeBin-increaseBegin<fMinTimeBins){//at least 2 time bins
		  timeBinMax=-1;
		  valueMaximum=-1;
		  increaseBegin=-1;
		  sumAdc=0;
		  ifIncreaseBegin=true;
		  ifMaximum=false;
		  continue;
		}
		//insert maximum, default values and set flags
		Double_t xCord,yCord;
		GetXY(iSec,iRow,iPad,xCord,yCord);
		AliPadMax *oneMaximum = new AliPadMax(AliTPCvtpr(valueMaximum,
								 timeBinMax,
								 iPad,
								 iRow,
								 xCord,
								 yCord,
								 timeBinMax),
						      increaseBegin,
						      iTimeBin-1,
						      sumAdc);
		maximaInSector.push_back(oneMaximum);
		
		timeBinMax=-1;
		valueMaximum=-1;
		increaseBegin=-1;
		sumAdc=0;
		ifIncreaseBegin=true;
		ifMaximum=false;
	      }
	      continue;
	    }
	    
	    if(ifIncreaseBegin){
	      ifIncreaseBegin=false;
	      increaseBegin=iTimeBin;
	    }
	    
	    if(adc>valueMaximum){
	      timeBinMax=iTimeBin;
	      valueMaximum=adc;
	      ifMaximum=true;
	    }
	    sumAdc+=adc;
	    if(iTimeBin==nTime-1 && ifMaximum && nTime-increaseBegin>fMinTimeBins){//on the edge
	      //at least 3 timebins
	      //insert maximum, default values and set flags
	      Double_t xCord,yCord;
	      GetXY(iSec,iRow,iPad,xCord,yCord);
	      AliPadMax *oneMaximum = new AliPadMax(AliTPCvtpr(valueMaximum,
							       timeBinMax,
							       iPad,
							       iRow,
							       xCord,
							       yCord,
							       timeBinMax),
						    increaseBegin,
						    iTimeBin-1,
						    sumAdc);
	      maximaInSector.push_back(oneMaximum);
	      
	      timeBinMax=-1;
	      valueMaximum=-1;
	      increaseBegin=-1;
	      sumAdc=0;
	      ifIncreaseBegin=true;
	      ifMaximum=false;
	      continue;
	    }
	    
	  }//end loop over timebins
	}//end loop over pads
//      }else{
//	cout<<"Pointer does not exist!!"<<endl;
      }//end if poiner exists
    }//end loop over rows
    
    //
    // Making clusters
    //

    Short_t maxDig=0;
    Short_t maxSumAdc=0;
    Short_t maxTimeBin=0;
    Short_t maxPad=0;
    Short_t maxRow=0;
    Double_t maxX=0;
    Double_t maxY=0;
    
    for( std::vector<AliPadMax*>::iterator mp1  = maximaInSector.begin();
	 mp1 != maximaInSector.end(); ++mp1 ) {
      
      AliTPCclusterKr *tmp=new AliTPCclusterKr();
      
      Short_t nUsedPads=1;
      Int_t clusterValue=0;
      clusterValue+=(*mp1)->GetSum();
      list<Short_t> nUsedRows;
      nUsedRows.push_back((*mp1)->GetRow());

      maxDig      =(*mp1)->GetAdc() ;
      maxSumAdc   =(*mp1)->GetSum() ;
      maxTimeBin  =(*mp1)->GetTime();
      maxPad      =(*mp1)->GetPad() ;
      maxRow      =(*mp1)->GetRow() ;
      maxX        =(*mp1)->GetX();
      maxY        =(*mp1)->GetY();

      AliSimDigits *digrowTmp;
      if(fRawData){
	digrowTmp = (AliSimDigits*)fDigarr->GetRow(iSec,(*mp1)->GetRow());
      }else{
	digrowTmp = (AliSimDigits*)fDigarr->LoadRow(iSec,(*mp1)->GetRow());
      }

      digrowTmp->ExpandBuffer(); //decrunch

      for(Short_t itb=(*mp1)->GetBegin(); itb<((*mp1)->GetEnd())+1; itb++){
	Short_t adcTmp = digrowTmp->GetDigitFast(itb,(*mp1)->GetPad());
	AliTPCvtpr *vtpr=new AliTPCvtpr(adcTmp,itb,(*mp1)->GetPad(),(*mp1)->GetRow(),(*mp1)->GetX(),(*mp1)->GetY(),itb);
	tmp->fCluster.push_back(vtpr);

      }
      tmp->SetCenter();//set centr of the cluster
      
      maximaInSector.erase(mp1);
      mp1--;
      
      for( std::vector<AliPadMax*>::iterator mp2  = maximaInSector.begin();
	   mp2 != maximaInSector.end(); ++mp2 ) {
//	if(abs(maxRow - (*mp2)->GetRow()) < fMaxRowRange && //3
//         abs(maxPad - (*mp2)->GetPad()) < fMaxPadRange && //4
	  
	if(fabs(tmp->GetCenterX() - (*mp2)->GetX()) < fMaxPadRangeCm &&
	   fabs(tmp->GetCenterY() - (*mp2)->GetY()) < fMaxRowRangeCm &&
	   fabs(tmp->GetCenterT() - (*mp2)->GetT()) < fMaxTimeRange){//7
	  
	  clusterValue+=(*mp2)->GetSum();

	  nUsedPads++;
	  nUsedRows.push_back((*mp2)->GetRow());

	  AliSimDigits *digrowTmp1;
	  if(fRawData){
	    digrowTmp1 = (AliSimDigits*)fDigarr->GetRow(iSec,(*mp2)->GetRow());
	  }else{
	    digrowTmp1 = (AliSimDigits*)fDigarr->LoadRow(iSec,(*mp2)->GetRow());
	  }
	  digrowTmp1->ExpandBuffer(); //decrunch
	  
	  for(Short_t itb=(*mp2)->GetBegin(); itb<(*mp2)->GetEnd()+1; itb++){
	    Short_t adcTmp = digrowTmp1->GetDigitFast(itb,(*mp2)->GetPad());
	    AliTPCvtpr *vtpr=new AliTPCvtpr(adcTmp,itb,(*mp2)->GetPad(),(*mp2)->GetRow(),(*mp2)->GetX(),(*mp2)->GetY(),itb);
	    tmp->fCluster.push_back(vtpr);
	  }
	  
	  tmp->SetCenter();//set centr of the cluster

	  //which one is bigger
	  if( (*mp2)->GetAdc() > maxDig ){
	    maxDig      =(*mp2)->GetAdc() ;
	    maxSumAdc   =(*mp2)->GetSum() ;
	    maxTimeBin  =(*mp2)->GetTime();
	    maxPad      =(*mp2)->GetPad() ;
	    maxRow      =(*mp2)->GetRow() ;
	    maxX        =(*mp2)->GetX() ;
	    maxY        =(*mp2)->GetY() ;
	  } else if ( (*mp2)->GetAdc() == maxDig ){
	    if( (*mp2)->GetSum() > maxSumAdc){
	      maxDig      =(*mp2)->GetAdc() ;
	      maxSumAdc   =(*mp2)->GetSum() ;
	      maxTimeBin  =(*mp2)->GetTime();
	      maxPad      =(*mp2)->GetPad() ;
	      maxRow      =(*mp2)->GetRow() ;
	      maxX        =(*mp2)->GetX() ;
	      maxY        =(*mp2)->GetY() ;
	    }
	  }
	  maximaInSector.erase(mp2);
	  mp2--;
	}
      }//inside loop
      
      //through out ADC=6,7 on 1 pad, 2 tb and ADC=12 on 2 pads,2 tb
      //if(nUsedPads==1 && clusterValue/tmp->fCluster.size()<3.6)continue;
      //if(nUsedPads==2 && clusterValue/tmp->fCluster.size()<3.1)continue;

      //through out clusters on the edge and noise
      if(clusterValue/tmp->fCluster.size()<fValueToSize)continue;
      
      tmp->SetADCcluster(clusterValue);
      tmp->SetNPads(nUsedPads);
      tmp->SetMax(AliTPCvtpr(maxDig,maxTimeBin,maxPad,maxRow,maxX,maxY,maxTimeBin));
      tmp->SetSec(iSec);
      tmp->SetSize();

      nUsedRows.sort();
      nUsedRows.unique();
      tmp->SetNRows(nUsedRows.size());
      tmp->SetCenter();

      clusterCounter++;
      
      //save each cluster into file

      AliTPCClustersRow *clrow= new AliTPCClustersRow();
      fRowCl=clrow;
      clrow->SetClass("AliTPCclusterKr");
      clrow->SetArray(1);
      fOutput->GetBranch("Segment")->SetAddress(&clrow);

      Int_t tmpCluster=0;
      TClonesArray * arr = fRowCl->GetArray();
      AliTPCclusterKr * cl = new ((*arr)[tmpCluster]) AliTPCclusterKr(*tmp);
      
      fOutput->Fill(); 
      delete clrow;

      //end of save each cluster into file adc.root
    }//outer loop

  }//end sector for
  cout<<"Number of clusters in event: "<<clusterCounter<<endl;

  return 0;
}

////____________________________________________________________________________


void AliTPCclustererKr::GetXY(Short_t sec,Short_t row,Short_t pad,Double_t& xGlob,Double_t& yGlob){

  Double_t yLocal = fParam->GetPadRowRadii(sec,row);//radius of row in sector in cm

  Short_t padmax = fParam->GetNPads(sec,row);//number of pads in a given row
  Float_t padXSize;
  if(sec<fParam->GetNInnerSector())padXSize=0.4;
  else padXSize=0.6;
  Double_t xLocal=(pad+0.5-padmax)*padXSize;//x-value of the center of pad

  Float_t sin,cos;
  fParam->AdjustCosSin((Int_t)sec,cos,sin);//return sinus and cosinus of the sector

  xGlob =  xLocal * cos + yLocal * sin;
  yGlob = -xLocal * sin + yLocal * cos;

  return;
}
