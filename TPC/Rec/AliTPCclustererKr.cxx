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
//           Implementation of the TPC Kr cluster class
//
// Origin: Adam Matyja, INP PAN, adam.matyja@ifj.edu.pl
//-----------------------------------------------------------------

/*
Instruction - how to use that:
There are two macros prepared. One is for preparing clusters from MC 
samples:  FindKrClusters.C. The output is kept in TPC.RecPoints.root.
The other macro is prepared for data analysis: FindKrClustersRaw.C. 
The output is created for each processed file in root file named adc.root. 
For each data subsample the same named file is created. So be careful 
do not overwrite them. 

Additional selection criteria to select the GOLD cluster
Example:
// open  file with clusters
TFile f("Krypton.root");
TTree * tree = (TTree*)f.Get("Kr")
TCut cutR0("cutR0","fADCcluster/fSize<100");        // adjust it according v seetings - 
TCut cutR1("cutR1","fADCcluster/fSize>7");          // cosmic tracks and noise removal
TCut cutR2("cutR2","fMax.fAdc/fADCcluster<0.2");    // digital noise removal
TCut cutR3("cutR3","fMax.fAdc/fADCcluster>0.01");   // noise removal
TCut cutS1("cutS1","fSize<200");    // adjust it according v seetings - cosmic tracks
TCut cutAll = cutR0+cutR1+cutR2+cutR3+cutS1;
This values are typical values to be applied in selectors


*
**** MC ****
*

To run clusterizaton for MC type:
.x FindKrClusters.C

If you don't want to use the standard selection criteria then you 
have to do following:

// load the standard setup
AliRunLoader* rl = AliRunLoader::Open("galice.root");
AliTPCLoader *tpcl = (AliTPCLoader*)rl->GetLoader("TPCLoader");
tpcl->LoadDigits();
rl->LoadgAlice();
gAlice=rl->GetAliRun();
TDirectory *cwd = gDirectory;
AliTPCv4 *tpc = (AliTPCv4*)gAlice->GetDetector("TPC");
Int_t ver = tpc->IsVersion();
rl->CdGAFile();
AliTPCParam *param=(AliTPCParamSR *)gDirectory->Get("75x40_100x60_150x60");
AliTPCDigitsArray *digarr=new AliTPCDigitsArray;
digarr->Setup(param);
cwd->cd();

//loop over events
Int_t nevmax=rl->GetNumberOfEvents();
for(Int_t nev=0;nev<nevmax ;nev++){
  rl->GetEvent(nev);
  TTree* input_tree= tpcl->TreeD();//load tree with digits
  digarr->ConnectTree(input_tree);
  TTree *output_tree =tpcl->TreeR();//load output tree

  AliTPCclustererKr *clusters = new AliTPCclustererKr();
  clusters->SetParam(param);
  clusters->SetInput(input_tree);
  clusters->SetOutput(output_tree);
  clusters->SetDigArr(digarr);
  
//If you want to change the cluster finder parameters for MC there are 
//several of them:

//1. signal threshold (everything below the given number is treated as 0)
  clusters->SetMinAdc(3);

//2. number of neighbouring timebins to be considered
  clusters->SetMinTimeBins(2);

//3. distance of the cluster center to the center of a pad in pad-padrow plane 
//(in cm). Remenber that this is still quantified by pad size.
  clusters->SetMaxPadRangeCm(2.5);

//4. distance of the cluster center to the center of a padrow in pad-padrow 
//plane (in cm). Remenber that this is still quantified by pad size.
  clusters->SetMaxRowRangeCm(3.5);

//5. distance of the cluster center to the max time bin on a pad (in tackts)
//ie. fabs(centerT - time)<7
  clusters->SetMaxTimeRange(7);

//6. cut reduce peak at 0. There are noises which appear mostly as two 
//timebins on one pad.
  clusters->SetValueToSize(3.5);


  clusters->finderIO();
  tpcl->WriteRecPoints("OVERWRITE");
}
delete rl;//cleans everything

*
********* DATA *********
*

To run clusterizaton for DATA for file named raw_data.root type:
.x FindKrClustersRaw.C("raw_data.root")

If you want to change some criteria do the following:

//
// remove Altro warnings
//
AliLog::SetClassDebugLevel("AliRawReaderDate",-5);
AliLog::SetClassDebugLevel("AliTPCAltroMapping",-5);
AliLog::SetModuleDebugLevel("RAW",-5);

//
// Get database with noises
//
//  char *ocdbpath = gSystem->Getenv("OCDB_PATH");
char *ocdbpath ="local:///afs/cern.ch/alice/tpctest/OCDB";
if (ocdbpath==0){
ocdbpath="alien://folder=/alice/data/2007/LHC07w/OCDB/";
}
AliCDBManager * man = AliCDBManager::Instance();
man->SetDefaultStorage(ocdbpath);
man->SetRun(0);
AliTPCCalPad * noiseTPC = AliTPCcalibDB::Instance()->GetPadNoise();
AliTPCAltroMapping** mapping =AliTPCcalibDB::Instance()->GetMapping();

//define tree
TFile *hfile=new TFile("adc.root","RECREATE","ADC file");
// Create a ROOT Tree
TTree *mytree = new TTree("Kr","Krypton cluster tree");

//define infput file
const char *fileName="data.root";
AliRawReader *reader = new AliRawReaderRoot(fileName);
//AliRawReader *reader = new AliRawReaderDate(fileName);
reader->Reset();
AliAltroRawStreamFast* stream = new AliAltroRawStreamFast(reader);
stream->SelectRawData("TPC");

//one general output
AliTPCclustererKr *clusters = new AliTPCclustererKr();
clusters->SetOutput(mytree);
clusters->SetRecoParam(0);//standard reco parameters
AliTPCParamSR *param=new AliTPCParamSR();
clusters->SetParam(param);//TPC parameters(sectors, timebins, etc.)

//set cluster finder parameters (from data):
//1. zero suppression parameter
  clusters->SetZeroSup(param->GetZeroSup());

//2. first bin
  clusters->SetFirstBin(60);

//3. last bin
  clusters->SetLastBin(950);

//4. maximal noise
  clusters->SetMaxNoiseAbs(2);

//5. maximal amount of sigma of noise
  clusters->SetMaxNoiseSigma(3);

//The remaining parameters are the same paramters as for MC (see MC section 
//points 1-6)
  clusters->SetMinAdc(3);
  clusters->SetMinTimeBins(2);
  clusters->SetMaxPadRangeCm(2.5);
  clusters->SetMaxRowRangeCm(3.5);
  clusters->SetMaxTimeRange(7);
  clusters->SetValueToSize(3.5);

while (reader->NextEvent()) {
  clusters->FinderIO(reader);
}

hfile->Write();
hfile->Close();
delete stream;


*/

#include "AliTPCclustererKr.h"
#include "AliTPCclusterKr.h"
//#include <vector>
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
#include "TH1F.h"
#include "TH2F.h"
#include "TTreeStream.h"

#include "AliTPCTransform.h"

//used in raw data finder
#include "AliTPCROC.h"
#include "AliTPCCalPad.h"
#include "AliTPCAltroMapping.h"
#include "AliTPCcalibDB.h"
#include "AliTPCRawStreamV3.h"
#include "AliTPCRecoParam.h"
#include "AliTPCReconstructor.h"
#include "AliRawReader.h"
#include "AliTPCCalROC.h"
#include "AliRawEventHeaderBase.h"

using std::cerr;
using std::cout;
using std::endl;
using std::list;
ClassImp(AliTPCclustererKr)


AliTPCclustererKr::AliTPCclustererKr()
  :TObject(),
  fRawData(kFALSE),
  fInput(0),
  fOutput(0),
  fParam(0),
  fDigarr(0),
  fRecoParam(0),
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
  fMaxRowRangeCm(3.5),
  fIsolCut(3),
  fDebugLevel(-1),
  fHistoRow(0),
  fHistoPad(0),
  fHistoTime(0),
   fHistoRowPad(0),
   fTimeStamp(0),
  fRun(0)
{
//
// default constructor
//
}

AliTPCclustererKr::AliTPCclustererKr(const AliTPCclustererKr &param)
  :TObject(),
  fRawData(kFALSE),
  fInput(0),
  fOutput(0),
  fParam(0),
  fDigarr(0),
  fRecoParam(0),
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
  fMaxRowRangeCm(3.5),
  fIsolCut(3),
  fDebugLevel(-1),
  fHistoRow(0),
  fHistoPad(0),
  fHistoTime(0),
  fHistoRowPad(0),
   fTimeStamp(0),
   fRun(0)
{
//
// copy constructor
//
  fParam = param.fParam;
  fRecoParam = param.fRecoParam;
  fRawData = param.fRawData;
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
  fIsolCut = param.fIsolCut;
  fDebugLevel = param.fDebugLevel;
  fHistoRow    = param.fHistoRow   ;
  fHistoPad    = param.fHistoPad  ;
  fHistoTime   = param.fHistoTime;
  fHistoRowPad = param.fHistoRowPad;
  fTimeStamp = param.fTimeStamp;
  fRun = param.fRun;

} 

AliTPCclustererKr & AliTPCclustererKr::operator = (const AliTPCclustererKr & param)
{
  //
  // assignment operator
  //
  if (this == &param) return (*this);
  
  fParam = param.fParam;
  fRecoParam = param.fRecoParam;
  fRawData = param.fRawData;
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
  fIsolCut = param.fIsolCut;
  fDebugLevel = param.fDebugLevel;
  fHistoRow    = param.fHistoRow   ;
  fHistoPad    = param.fHistoPad  ;
  fHistoTime   = param.fHistoTime;
  fHistoRowPad = param.fHistoRowPad;
  fTimeStamp = param.fTimeStamp;
  fRun = param.fRun;
  return (*this);
}

AliTPCclustererKr::~AliTPCclustererKr()
{
  //
  // destructor
  //
  delete fOutput;
}

void AliTPCclustererKr::SetRecoParam(AliTPCRecoParam *recoParam)
{
  //
  // set reconstruction parameters
  //
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

void AliTPCclustererKr::SetOutput(TTree * /*tree*/) 
{
  //
  //dummy
  //
  fOutput = new TTreeSRedirector("Krypton.root");
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
  // this method is unsing AliAltroRawStreamV3
  // fParam must be defined before
  if (!rawReader) return 1;
  //
  fRawData=kTRUE; //set flag to data
  
  if (!fOutput) {
    Error("Digits2Clusters", "output tree not initialised");
    return 11;
  }
  
  fParam->SetMaxTBin(fRecoParam->GetLastBin());//set number of timebins from reco -> param
  //   used later for memory allocation

  AliRawEventHeaderBase* eventHeader = (AliRawEventHeaderBase*)rawReader->GetEventHeader();
  if (eventHeader){
    fTimeStamp = eventHeader->Get("Timestamp");
    fRun = rawReader->GetRunNumber();
  }


  Bool_t isAltro=kFALSE;
  
  AliTPCROC * roc = AliTPCROC::Instance();
  AliTPCCalPad * noiseTPC = AliTPCcalibDB::Instance()->GetPadNoise();
  AliTPCAltroMapping** mapping =AliTPCcalibDB::Instance()->GetMapping();
  //
  AliTPCRawStreamV3 input(rawReader,(AliAltroMapping**)mapping);
  
  const Int_t kNIS = fParam->GetNInnerSector();//number of inner sectors
  const Int_t kNOS = fParam->GetNOuterSector();//number of outer sectors
  const Int_t kNS = kNIS + kNOS;//all sectors
  
  
  //crate TPC view
  AliTPCDigitsArray *digarr=new AliTPCDigitsArray(kFALSE);//data not sim
  digarr->Setup(fParam);//as usually parameters
  
  for(Int_t iSec = 0; iSec < kNS; iSec++) {
    AliTPCCalROC * noiseROC;
    AliTPCCalROC noiseDummy(iSec);
    if(noiseTPC==0x0){
      noiseROC = &noiseDummy;//noise=0
    }else{
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
    rawReader->Select("TPC",indexDDL,indexDDL+nDDLs-1);
      
    
    while (input.NextDDL()){
      // Allocate memory for rows in sector (pads(depends on row) x timebins)
      if (!digarr->GetRow(iSec,0)){
        for(Int_t iRow = 0; iRow < nRows; iRow++) {
          digarr->CreateRow(iSec,iRow);
        }//end loop over rows
      }
      //loop over pads
      while ( input.NextChannel() ) {
        Int_t iRow = input.GetRow();
        Int_t iPad = input.GetPad();
        //check row consistency
        if (iRow < 0 ) continue;
        if (iRow < 0 || iRow >= nRows){
          AliError(Form("Pad-row index (%d) outside the range (%d -> %d) !",
                        iRow, 0, nRows -1));
          continue;
        }
        
      //check pad consistency
        if (iPad < 0 || iPad >= (Int_t)(roc->GetNPads(iSec,iRow))) {
          AliError(Form("Pad index (%d) outside the range (%d -> %d) !",
                        iPad, 0, roc->GetNPads(iSec,iRow) ));
          continue;
        }
        
      //loop over bunches
        while ( input.NextBunch() ){
          Int_t  startTbin    = (Int_t)input.GetStartTimeBin();
          Int_t  bunchlength  = (Int_t)input.GetBunchLength();
          const UShort_t *sig = input.GetSignals();
          isAltro=kTRUE;
          for (Int_t iTime = 0; iTime<bunchlength; iTime++){
            Int_t iTimeBin=startTbin-iTime;
            //
            if(fDebugLevel==72){
              fHistoRow->Fill(iRow);
              fHistoPad->Fill(iPad);
              fHistoTime->Fill(iTimeBin);
              fHistoRowPad->Fill(iPad,iRow);
            }else if(fDebugLevel>=0&&fDebugLevel<72){
              if(iSec==fDebugLevel){
                fHistoRow->Fill(iRow);
                fHistoPad->Fill(iPad);
                fHistoTime->Fill(iTimeBin);
                fHistoRowPad->Fill(iPad,iRow);
              }
            }else if(fDebugLevel==73){
              if(iSec<36){
                fHistoRow->Fill(iRow);
                fHistoPad->Fill(iPad);
                fHistoTime->Fill(iTimeBin);
                fHistoRowPad->Fill(iPad,iRow);
              }
            }else if(fDebugLevel==74){
              if(iSec>=36){
                fHistoRow->Fill(iRow);
                fHistoPad->Fill(iPad);
                fHistoTime->Fill(iTimeBin);
                fHistoRowPad->Fill(iPad,iRow);
              }
            }
            
            //check time consistency
            if ( iTimeBin < fRecoParam->GetFirstBin() || iTimeBin >= fRecoParam->GetLastBin()){
              //cout<<iTimeBin<<endl;
              continue;
              AliFatal(Form("Timebin index (%d) outside the range (%d -> %d) !",
                            iTimeBin, 0, fRecoParam->GetLastBin() -1));
            }
            //signal
            Float_t signal=(Float_t)sig[iTime];
            if (signal <= fZeroSup ||
                iTimeBin < fFirstBin ||
                iTimeBin > fLastBin
               ) {
                 digarr->GetRow(iSec,iRow)->SetDigitFast(0,iTimeBin,iPad);
                 continue;
               }
            if (!noiseROC) continue;
            Double_t noiseOnPad = noiseROC->GetValue(iRow,iPad);//noise on given pad and row in sector
            if (noiseOnPad > fMaxNoiseAbs){
              digarr->GetRow(iSec,iRow)->SetDigitFast(0,iTimeBin,iPad);
              continue; // consider noisy pad as dead
            }
            if(signal <= fMaxNoiseSigma * noiseOnPad){
              digarr->GetRow(iSec,iRow)->SetDigitFast(0,iTimeBin,iPad);
              continue;
            }
            digarr->GetRow(iSec,iRow)->SetDigitFast(TMath::Nint(signal),iTimeBin,iPad);
          }// end loop signals in bunch
        }// end loop bunches
      } // end loop pads
    }// end ddl loop
  }// end sector loop
  SetDigArr(digarr);
  if(isAltro) FindClusterKrIO();
  delete digarr;
  
  return 0;
}

void AliTPCclustererKr::CleanSector(Int_t sector){
  //
  // clean isolated digits
  //  
  const Int_t kNRows=fParam->GetNRow(sector);//number of rows in sector
  for(Int_t iRow=0; iRow<kNRows; ++iRow){
    AliSimDigits *digrow;
    if(fRawData){
      digrow = (AliSimDigits*)fDigarr->GetRow(sector,iRow);//real data
    }else{
      digrow = (AliSimDigits*)fDigarr->LoadRow(sector,iRow);//MC
    }
    if(!digrow) continue;
    digrow->ExpandBuffer(); //decrunch
    const Int_t kNPads = digrow->GetNCols();  // number of pads
    const Int_t kNTime = digrow->GetNRows(); // number of timebins
    for(Int_t iPad=1;iPad<kNPads-1;iPad++){
      Short_t*  val = digrow->GetDigitsColumn(iPad);

      for(Int_t iTimeBin=1;iTimeBin<kNTime-1;iTimeBin++){
	if (val[iTimeBin]<=0) continue;
	if (val[iTimeBin-1]+val[iTimeBin+1]<fIsolCut) {val[iTimeBin]=0; continue;}
	if (val[iTimeBin-kNTime]+val[iTimeBin+kNTime]<fIsolCut) {val[iTimeBin]=0; continue;}
	//
	if (val[iTimeBin-1-kNTime]+val[iTimeBin+1+kNTime]<fIsolCut) {val[iTimeBin]=0; continue;}
	if (val[iTimeBin+1-kNTime]+val[iTimeBin-1+kNTime]<fIsolCut) {val[iTimeBin]=0; continue;}

      }
    }
  }
}


////____________________________________________________________________________
Int_t AliTPCclustererKr::FindClusterKrIO()
{

  //
  //fParam and  fDigarr must be set to run this method
  //

  Int_t clusterCounter=0;
  const Int_t nTotalSector=fParam->GetNSector();//number of sectors
  for(Int_t iSec=0; iSec<nTotalSector; ++iSec){
    CleanSector(iSec);

    //vector of maxima for each sector
    //std::vector<AliPadMax*> maximaInSector;
    TObjArray *maximaInSector=new TObjArray();//to store AliPadMax*

    //
    //  looking for the maxima on the pad
    //

    const Int_t kNRows=fParam->GetNRow(iSec);//number of rows in sector
    for(Int_t iRow=0; iRow<kNRows; ++iRow){
      AliSimDigits *digrow;
      if(fRawData){
	digrow = (AliSimDigits*)fDigarr->GetRow(iSec,iRow);//real data
      }else{
	digrow = (AliSimDigits*)fDigarr->LoadRow(iSec,iRow);//MC
      }
      if(digrow){//if pointer exist
	digrow->ExpandBuffer(); //decrunch
	const Int_t kNPads = digrow->GetNCols();  // number of pads
	const Int_t kNTime = digrow->GetNRows(); // number of timebins
	for(Int_t iPad=0;iPad<kNPads;iPad++){
	  
	  Int_t timeBinMax=-1;//timebin of maximum 
	  Int_t valueMaximum=-1;//value of maximum in adc
	  Int_t increaseBegin=-1;//timebin when increase starts
	  Int_t sumAdc=0;//sum of adc on the pad in maximum surrounding
	  bool ifIncreaseBegin=true;//flag - check if increasing started
	  bool ifMaximum=false;//flag - check if it could be maximum
	  Short_t* val = digrow->GetDigitsColumn(iPad);
	  for(Int_t iTimeBin=1;iTimeBin<kNTime-1;iTimeBin++){
	    if (!ifMaximum)  {
	      if (val[iTimeBin]==-1) break;   // 0 until the end
	      for( ; iTimeBin<kNTime-2&&val[iTimeBin]<fMinAdc ;iTimeBin++) {}
	    }
	    //
	    Short_t adc = val[iTimeBin];

	    if(adc<fMinAdc){//standard was 3 for fMinAdc
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
		//Double_t xCord,yCord;
		//GetXY(iSec,iRow,iPad,xCord,yCord);
		Double_t x[]={iRow,iPad,iTimeBin};
		Int_t i[]={iSec};
		AliTPCTransform *transform     = AliTPCcalibDB::Instance()->GetTransform() ;

		transform->Transform(x,i,0,1);
		
		AliPadMax *oneMaximum = new AliPadMax(AliTPCvtpr(valueMaximum,
								 timeBinMax,
								 iPad,
								 iRow,
								 x[0],//xCord,
								 x[1],//yCord,
								 x[2]/*timeBinMax*/),
						      increaseBegin,
						      iTimeBin-1,
						      sumAdc);
		maximaInSector->AddLast(oneMaximum);
		
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
	    if(iTimeBin==kNTime-1 && ifMaximum && kNTime-increaseBegin>fMinTimeBins){//on the edge
	      //at least 3 timebins
	      //insert maximum, default values and set flags
	      //Double_t xCord,yCord;
	      //GetXY(iSec,iRow,iPad,xCord,yCord);
	      Double_t x[]={iRow,iPad,iTimeBin};
	      Int_t i[]={iSec};
	      //AliTPCTransform trafo;
	      //trafo.Transform(x,i,0,1);

		AliTPCTransform *transform     = AliTPCcalibDB::Instance()->GetTransform() ;

		transform->Transform(x,i,0,1);

	      AliPadMax *oneMaximum = new AliPadMax(AliTPCvtpr(valueMaximum,
							       timeBinMax,
							       iPad,
							       iRow,
							       x[0],//xCord,
							       x[1],//yCord,
							       x[2]/*timeBinMax*/),
						    increaseBegin,
						    iTimeBin-1,
						    sumAdc);
	      maximaInSector->AddLast(oneMaximum);
		
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

    MakeClusters(maximaInSector,iSec,clusterCounter);
    //
    maximaInSector->SetOwner(kTRUE);
    maximaInSector->Delete();
    delete maximaInSector;
  }//end sector for
  cout<<"Number of clusters in event: "<<clusterCounter<<endl;
  return 0;
}

void AliTPCclustererKr::MakeClusters(TObjArray * maximaInSector, Int_t iSec, Int_t &clusterCounter){
  //
  // Make clusters
  //

  Int_t maxDig=0;
  Int_t maxSumAdc=0;
  Int_t maxTimeBin=0;
  Int_t maxPad=0;
  Int_t maxRow=0;
  Double_t maxX=0;
  Double_t maxY=0;
  Double_t maxT=0;
  Int_t entriesArr = maximaInSector->GetEntriesFast();
  for(Int_t it1 = 0; it1 < entriesArr; ++it1 ) {
    
    AliPadMax *mp1=(AliPadMax *)maximaInSector->UncheckedAt(it1);
    if (!mp1) continue;
    AliTPCclusterKr clusterKr;
    
    Int_t nUsedPads=1;
    Int_t clusterValue=0;
    clusterValue+=(mp1)->GetSum();
    list<Int_t> nUsedRows;
    nUsedRows.push_back((mp1)->GetRow());
    
    maxDig      =(mp1)->GetAdc() ;
    maxSumAdc   =(mp1)->GetSum() ;
    maxTimeBin  =(mp1)->GetTime();
    maxPad      =(mp1)->GetPad() ;
    maxRow      =(mp1)->GetRow() ;
    maxX        =(mp1)->GetX();
    maxY        =(mp1)->GetY();
    maxT        =(mp1)->GetT();
    
    AliSimDigits *digrowTmp;
    if(fRawData){
      digrowTmp = (AliSimDigits*)fDigarr->GetRow(iSec,(mp1)->GetRow());
    }else{
      digrowTmp = (AliSimDigits*)fDigarr->LoadRow(iSec,(mp1)->GetRow());
    }
    
    digrowTmp->ExpandBuffer(); //decrunch
    
    for(Int_t itb=(mp1)->GetBegin(); itb<((mp1)->GetEnd())+1; itb++){
      Int_t adcTmp = digrowTmp->GetDigitUnchecked(itb,(mp1)->GetPad());
      AliTPCvtpr *vtpr=new AliTPCvtpr(adcTmp,itb,(mp1)->GetPad(),(mp1)->GetRow(),(mp1)->GetX(),(mp1)->GetY(),(mp1)->GetT());
      clusterKr.AddDigitToCluster(vtpr);
    }
    clusterKr.SetCenter();//set centr of the cluster
    
    for(Int_t it2 = it1+1; it2 < entriesArr; ++it2 ) {
      AliPadMax *mp2=(AliPadMax *)maximaInSector->UncheckedAt(it2);
      if (!mp2) continue;
      if (TMath::Abs(clusterKr.GetCenterX() - (mp2)->GetX()) > fMaxPadRangeCm) continue;
      if (TMath::Abs(clusterKr.GetCenterY() - (mp2)->GetY()) > fMaxRowRangeCm) continue;      
      if (TMath::Abs(clusterKr.GetCenterT() - (mp2)->GetT()) > fMaxTimeRange) continue;

      {
	clusterValue+=(mp2)->GetSum();
	
	nUsedPads++;
	nUsedRows.push_back((mp2)->GetRow());
	
	AliSimDigits *digrowTmp1;
	if(fRawData){
	  digrowTmp1 = (AliSimDigits*)fDigarr->GetRow(iSec,(mp2)->GetRow());
	}else{
	  digrowTmp1 = (AliSimDigits*)fDigarr->LoadRow(iSec,(mp2)->GetRow());
	}
	digrowTmp1->ExpandBuffer(); //decrunch
	
	for(Int_t itb=(mp2)->GetBegin(); itb<(mp2)->GetEnd()+1; itb++){
	  Int_t adcTmp = digrowTmp1->GetDigitUnchecked(itb,(mp2)->GetPad());
	  AliTPCvtpr *vtpr=new AliTPCvtpr(adcTmp,itb,(mp2)->GetPad(),(mp2)->GetRow(),(mp2)->GetX(),(mp2)->GetY(),(mp2)->GetT());
	  clusterKr.AddDigitToCluster(vtpr);
	}
	
	clusterKr.SetCenter();//set center of the cluster
	
	//which one is bigger
	if( (mp2)->GetAdc() > maxDig ){
	  maxDig      =(mp2)->GetAdc() ;
	  maxSumAdc   =(mp2)->GetSum() ;
	  maxTimeBin  =(mp2)->GetTime();
	  maxPad      =(mp2)->GetPad() ;
	  maxRow      =(mp2)->GetRow() ;
	  maxX        =(mp2)->GetX() ;
	  maxY        =(mp2)->GetY() ;
	  maxT        =(mp2)->GetT() ;
	} else if ( (mp2)->GetAdc() == maxDig ){
	  if( (mp2)->GetSum() > maxSumAdc){
	    maxDig      =(mp2)->GetAdc() ;
	    maxSumAdc   =(mp2)->GetSum() ;
	    maxTimeBin  =(mp2)->GetTime();
	    maxPad      =(mp2)->GetPad() ;
	    maxRow      =(mp2)->GetRow() ;
	    maxX        =(mp2)->GetX() ;
	    maxY        =(mp2)->GetY() ;
	    maxT        =(mp2)->GetT() ;
	  }
	}
	delete maximaInSector->RemoveAt(it2);
      }
    }//inside loop
    delete maximaInSector->RemoveAt(it1);          
    clusterKr.SetSize();
    //through out clusters on the edge and noise
    //if(clusterValue/clusterKr.fCluster.size()<fValueToSize)continue;
    if(clusterValue/(clusterKr.GetSize())<fValueToSize)continue;
    
    clusterKr.SetADCcluster(clusterValue);
    clusterKr.SetNPads(nUsedPads);
    clusterKr.SetMax(AliTPCvtpr(maxDig,maxTimeBin,maxPad,maxRow,maxX,maxY,maxT));
    clusterKr.SetSec(iSec);
    clusterKr.SetSize();
    
    nUsedRows.sort();
    nUsedRows.unique();
    clusterKr.SetNRows(nUsedRows.size());
    clusterKr.SetCenter();
    
    clusterKr.SetRMS();//Set pad,row,timebin RMS
    clusterKr.Set1D();//Set size in pads and timebins

    clusterKr.SetTimeStamp(fTimeStamp);
    clusterKr.SetRun(fRun);

    clusterCounter++;
    
    
    //save each cluster into file
    if (fOutput){
      (*fOutput)<<"Kr"<<
	"Cl.="<<&clusterKr<<
	"\n";
    }
    //end of save each cluster into file adc.root
  }//outer loop
}



////____________________________________________________________________________


void AliTPCclustererKr::GetXY(Int_t sec,Int_t row,Int_t pad,Double_t& xGlob,Double_t& yGlob){
  //
  //gives global XY coordinate of the pad
  //

  Double_t yLocal = fParam->GetPadRowRadii(sec,row);//radius of row in sector in cm

  Int_t padmax = fParam->GetNPads(sec,row);//number of pads in a given row
  Float_t padXSize;
  if(sec<fParam->GetNInnerSector())padXSize=0.4;
  else padXSize=0.6;
  Double_t xLocal=(pad+0.5-padmax/2.)*padXSize;//x-value of the center of pad

  Float_t sin,cos;
  fParam->AdjustCosSin((Int_t)sec,cos,sin);//return sinus and cosinus of the sector

  Double_t xGlob1 =  xLocal * cos + yLocal * sin;
  Double_t yGlob1 = -xLocal * sin + yLocal * cos;


  Double_t rot=0;
  rot=TMath::Pi()/2.;

  xGlob =  xGlob1 * TMath::Cos(rot) + yGlob1 * TMath::Sin(rot);
  yGlob = -xGlob1 * TMath::Sin(rot) + yGlob1 * TMath::Cos(rot);

   yGlob=-1*yGlob;
   if(sec<18||(sec>=36&&sec<54)) xGlob =-1*xGlob;


  return;
}
