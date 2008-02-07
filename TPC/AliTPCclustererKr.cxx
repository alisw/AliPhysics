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
  fIsOldRCUFormat(kFALSE)
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
  fIsOldRCUFormat(kFALSE)
{
//
// copy constructor
//
  fParam = param.fParam;
  fRecoParam=param.fRecoParam;
  fIsOldRCUFormat=param.fIsOldRCUFormat;
  fRawData=param.fRawData;
  fRowCl =param.fRowCl ;
  fInput =param.fInput ;
  fOutput=param.fOutput;
  fDigarr=param.fDigarr;
} 

AliTPCclustererKr & AliTPCclustererKr::operator = (const AliTPCclustererKr & param)
{
  fParam = param.fParam;
  fRecoParam=param.fRecoParam;
  fIsOldRCUFormat=param.fIsOldRCUFormat;
  fRawData=param.fRawData;
  fRowCl =param.fRowCl ;
  fInput =param.fInput ;
  fOutput=param.fOutput;
  fDigarr=param.fDigarr;
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
    cerr<<"AliTPCclusterKr::findClusterKr(): no proper input tree !\n";
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
Int_t AliTPCclustererKr::finderIO()
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

  findClusterKrIO();
  return 0;
}

Int_t AliTPCclustererKr::finderIO(AliRawReader* rawReader)
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
  //AliTPCCalPad * noiseTPC = AliTPCcalibDB::Instance()->GetPadNoise();
  AliTPCAltroMapping** mapping =AliTPCcalibDB::Instance()->GetMapping();
  //
  AliTPCRawStream input(rawReader,(AliAltroMapping**)mapping);

  const Int_t kNIS = fParam->GetNInnerSector();//number of inner sectors
  const Int_t kNOS = fParam->GetNOuterSector();//number of outer sectors
  const Int_t kNS = kNIS + kNOS;//all sectors
  Int_t zeroSup = fParam->GetZeroSup();//zero suppression parameter

  //crate TPC view
  AliTPCDigitsArray *digarr=new AliTPCDigitsArray(kFALSE);//data not sim
  digarr->Setup(fParam);//as usually parameters

  //
  // Loop over sectors
  //
  for(Int_t sec = 0; sec < kNS; sec++) {
    //AliTPCCalROC * noiseROC   = noiseTPC->GetCalROC(sec);  // noise per given sector
 
    Int_t nRows = 0; //number of rows in sector
    Int_t nDDLs = 0; //number of DDLs
    Int_t indexDDL = 0; //DDL index
    if (sec < kNIS) {
      nRows = fParam->GetNRowLow();
      nDDLs = 2;
      indexDDL = sec * 2;
    }else {
      nRows = fParam->GetNRowUp();
      nDDLs = 4;
      indexDDL = (sec-kNIS) * 4 + kNIS * 2;
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
      for(Int_t row = 0; row < nRows; row++) {
	digarr->CreateRow(sec,row);
      }//end loop over rows
    }
    rawReader->Select("TPC",indexDDL,indexDDL+nDDLs-1);

    //
    // Begin loop over altro data
    //
    while (input.Next()) {

      //check sector consistency
      if (input.GetSector() != sec)
	AliFatal(Form("Sector index mismatch ! Expected (%d), but got (%d) !",sec,input.GetSector()));
      
      //check row consistency
      Short_t iRow = input.GetRow();
      if (iRow < 0 || iRow >= nRows){
	AliError(Form("Pad-row index (%d) outside the range (%d -> %d) !",
		      iRow, 0, nRows -1));
	continue;
      }

      //check pad consistency
      Short_t iPad = input.GetPad();
      if (iPad < 0 || iPad >= (Short_t)(roc->GetNPads(sec,iRow))) {
	AliError(Form("Pad index (%d) outside the range (%d -> %d) !",
		      iPad, 0, roc->GetNPads(sec,iRow) ));
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
      if (signal <= zeroSup) {
	//continue;
	digarr->GetRow(sec,iRow)->SetDigitFast(0,iTimeBin,iPad);
      }
      (//(AliSimDigits*)
       digarr->GetRow(sec,iRow))->SetDigitFast(signal,iTimeBin,iPad);
    }//end of loop over altro data
  }//end of loop over sectors
  
  SetDigArr(digarr);
  if(isAltro) findClusterKrIO();
  delete digarr;

  return 0;
}

////____________________________________________________________________________
Int_t AliTPCclustererKr::findClusterKrIO()
{
  //fParam and  fDigarr must be set to run this method

  Int_t cluster_counter=0;
  const Short_t NTotalSector=fParam->GetNSector();//number of sectors
  for(Short_t sec=0; sec<NTotalSector; sec++){
    
    //vector of maxima for each sector
    std::vector<AliPadMax*> maxima_in_sector;
    
    const Short_t Nrows=fParam->GetNRow(sec);//number of rows in sector
    for(Short_t row=0; row<Nrows; row++){
      AliSimDigits *digrow;
      if(fRawData){
	digrow = (AliSimDigits*)fDigarr->GetRow(sec,row);//real data
      }else{
	digrow = (AliSimDigits*)fDigarr->LoadRow(sec,row);//MC
      }
      if(digrow){//if pointer exist
	digrow->ExpandBuffer(); //decrunch
	const Short_t npads = digrow->GetNCols();  // number of pads
	const Short_t ntime = digrow->GetNRows(); // number of timebins
	for(Short_t np=0;np<npads;np++){
	  
	  Short_t tb_max=-1;//timebin of maximum 
	  Short_t value_maximum=-1;//value of maximum in adc
	  Short_t increase_begin=-1;//timebin when increase starts
	  Short_t sum_adc=0;//sum of adc on the pad in maximum surrounding
	  bool if_increase_begin=true;//flag - check if increasing start
	  bool if_maximum=false;//flag - check if it could be maximum

	  for(Short_t nt=0;nt<ntime;nt++){
	    Short_t adc = digrow->GetDigitFast(nt,np);
	    if(adc<3){
	      if(if_maximum){
		if(nt-1-increase_begin<1){//at least 2 time bins
		  tb_max=-1;
		  value_maximum=-1;
		  increase_begin=-1;
		  sum_adc=0;
		  if_increase_begin=true;
		  if_maximum=false;
		  continue;
		}
		//insert maximum, default values and set flags
		AliPadMax *one_maximum = new AliPadMax(AliTPCvtpr(value_maximum,
								  tb_max,
								  np,
								  row),
						       increase_begin,
						       nt-1,
						       sum_adc);
		maxima_in_sector.push_back(one_maximum);
		
		tb_max=-1;
		value_maximum=-1;
		increase_begin=-1;
		sum_adc=0;
		if_increase_begin=true;
		if_maximum=false;
	      }
	      continue;
	    }
	    
	    if(if_increase_begin){
	      if_increase_begin=false;
	      increase_begin=nt;
	    }
	    
	    if(adc>value_maximum){
	      tb_max=nt;
	      value_maximum=adc;
	      if_maximum=true;
	    }
	    sum_adc+=adc;
	    if(nt==ntime-1 && if_maximum && ntime-1-increase_begin>2){//on the edge
	      //insert maximum, default values and set flags
	      AliPadMax *one_maximum = new AliPadMax(AliTPCvtpr(value_maximum,
								tb_max,
								np,
								row),
						     increase_begin,
						     nt-1,
						     sum_adc);
	      maxima_in_sector.push_back(one_maximum);
	      
	      tb_max=-1;
	      value_maximum=-1;
	      increase_begin=-1;
	      sum_adc=0;
	      if_increase_begin=true;
	      if_maximum=false;
	      continue;
	    }
	    
	  }//end loop over timebins
	}//end loop over pads
//      }else{
//	cout<<"Pointer does not exist!!"<<endl;
      }//end if poiner exists
    }//end loop over rows
    

    Short_t max_dig=0;
    Short_t max_sum_adc=0;
    Short_t max_nt=0;
    Short_t max_np=0;
    Short_t max_row=0;
    
    for( std::vector<AliPadMax*>::iterator mp1  = maxima_in_sector.begin();
	 mp1 != maxima_in_sector.end(); ++mp1 ) {
      
      AliTPCclusterKr *tmp=new AliTPCclusterKr();
      
      Short_t n_used_pads=1;
      Short_t cluster_value=0;
      cluster_value+=(*mp1)->GetSum();
      
      max_dig      =(*mp1)->GetAdc() ;
      max_sum_adc  =(*mp1)->GetSum() ;
      max_nt       =(*mp1)->GetTime();
      max_np       =(*mp1)->GetPad() ;
      max_row      =(*mp1)->GetRow() ;
      
      AliSimDigits *digrow_tmp;
      if(fRawData){
	digrow_tmp = (AliSimDigits*)fDigarr->GetRow(sec,(*mp1)->GetRow());
      }else{
	digrow_tmp = (AliSimDigits*)fDigarr->LoadRow(sec,(*mp1)->GetRow());
      }

      digrow_tmp->ExpandBuffer(); //decrunch
      for(Short_t tb=(*mp1)->GetBegin(); tb<((*mp1)->GetEnd())+1; tb++){
	Short_t adc_tmp = digrow_tmp->GetDigitFast(tb,(*mp1)->GetPad());
	AliTPCvtpr *vtpr=new AliTPCvtpr(adc_tmp,tb,(*mp1)->GetPad(),(*mp1)->GetRow());
	tmp->fCluster.push_back(vtpr);
      }
      
      maxima_in_sector.erase(mp1);
      mp1--;
      
      for( std::vector<AliPadMax*>::iterator mp2  = maxima_in_sector.begin();
	   mp2 != maxima_in_sector.end(); ++mp2 ) {
	
	Short_t two_pad_max=4;
	Short_t two_row_max=3;
	
	if(abs(max_row - (*mp2)->GetRow())<two_row_max && 
	   abs(max_np - (*mp2)->GetPad())<two_pad_max  &&
	   abs(max_nt - (*mp2)->GetTime())<7){
	  
	  cluster_value+=(*mp2)->GetSum();
	  n_used_pads++;
	  
	  AliSimDigits *digrow_tmp1;
	  if(fRawData){
	    digrow_tmp1 = (AliSimDigits*)fDigarr->GetRow(sec,(*mp2)->GetRow());
	  }else{
	    digrow_tmp1 = (AliSimDigits*)fDigarr->LoadRow(sec,(*mp2)->GetRow());
	  }
	  digrow_tmp1->ExpandBuffer(); //decrunch
	  
	  for(Short_t tb=(*mp2)->GetBegin(); tb<(*mp2)->GetEnd()+1; tb++){
	    Short_t adc_tmp = digrow_tmp1->GetDigitFast(tb,(*mp2)->GetPad());
	    AliTPCvtpr *vtpr=new AliTPCvtpr(adc_tmp,tb,(*mp2)->GetPad(),(*mp2)->GetRow());
	      tmp->fCluster.push_back(vtpr);
	  }
	  
	  //which one is bigger
	  if( (*mp2)->GetAdc() > max_dig ){
	    max_dig      =(*mp2)->GetAdc() ;
	    max_sum_adc  =(*mp2)->GetSum() ;
	    max_nt       =(*mp2)->GetTime();
	    max_np       =(*mp2)->GetPad() ;
	    max_row      =(*mp2)->GetRow() ;
	  } else if ( (*mp2)->GetAdc() == max_dig ){
	    if( (*mp2)->GetSum() > max_sum_adc){
	      max_dig      =(*mp2)->GetAdc() ;
	      max_sum_adc  =(*mp2)->GetSum() ;
	      max_nt       =(*mp2)->GetTime();
	      max_np       =(*mp2)->GetPad() ;
	      max_row      =(*mp2)->GetRow() ;
	    }
	  }
	  maxima_in_sector.erase(mp2);
	  mp2--;
	}
      }//inside loop
      
      //through out ADC=6,7 on 1 pad, 2 tb and ADC=12 on 2 pads,2 tb
      if(n_used_pads==1 && cluster_value/tmp->fCluster.size()<3.6)continue;
      if(n_used_pads==2 && cluster_value/tmp->fCluster.size()<3.1)continue;
      
      tmp->SetADCcluster(cluster_value);
      tmp->SetNpads(n_used_pads);
      tmp->SetMax(AliTPCvtpr(max_dig,max_nt,max_np,max_row));
      tmp->SetSec(sec);
      tmp->SetSize();
      
      cluster_counter++;
      
      //save each cluster into file

      AliTPCClustersRow *clrow= new AliTPCClustersRow();
      fRowCl=clrow;
      clrow->SetClass("AliTPCclusterKr");
      clrow->SetArray(1);
      fOutput->GetBranch("Segment")->SetAddress(&clrow);

      Int_t tmp_cluster=0;
      TClonesArray * arr = fRowCl->GetArray();
      AliTPCclusterKr * cl = new ((*arr)[tmp_cluster]) AliTPCclusterKr(*tmp);
      
      fOutput->Fill(); 
      delete clrow;

      //end of save each cluster into file adc.root
    }//outer loop

  }//end sector for
  cout<<"Number of clusters in event: "<<cluster_counter<<endl;

  return 0;
}

////____________________________________________________________________________
