/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes hereby granted      *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TH1C.h>
#include <TList.h>
#include <TCanvas.h>
#include <TGeoManager.h>
#include <TRefArray.h>
#include <TKey.h>

#include "AliLog.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliESDpid.h"
#include "AliTOFcalib.h"
#include "AliCDBManager.h"
#include "AliRunTag.h"

#include "AliTOFT0maker.h"
#include "AliVCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVCaloCells.h"
#include "AliESDCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAODCaloCells.h"
#include "AliEMCALGeometry.h"
#include "AliOADBContainer.h"

#include "AliAnalysisTaskEMCALTimeCalib.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEMCALTimeCalib) ;
/// \endcond

using std::cout;
using std::endl;

//________________________________________________________________________
/// Constructor
AliAnalysisTaskEMCALTimeCalib::AliAnalysisTaskEMCALTimeCalib(const char *name)
: AliAnalysisTaskSE(name),
  fRunNumber(-1),
  fTOFmaker(0),
  fOutputList(0x0),
  fgeom(0),
  fGeometryName(),
  fMinClusterEnergy(0),
  fMaxClusterEnergy(0),
  fMinNcells(0),
  fMaxNcells(0),
  fMinLambda0(0),
  fMaxLambda0(0),
  fMinLambda0LG(0),
  fMaxLambda0LG(0),
  fMaxRtrack(0),
  fMinCellEnergy(0),
  fReferenceFileName(),
  fReferenceRunByRunFileName(),
  fPileupFromSPD(kFALSE),
  fMinTime(0),
  fMaxTime(0),
  fMostEneCellOnly(kFALSE),
  fRawTimeNbins (0),
  fRawTimeMin   (0),
  fRawTimeMax   (0),
  fPassTimeNbins(0),
  fPassTimeMin  (0),
  fPassTimeMax  (0),
  fEnergyNbins  (0),
  fEnergyMin(0),
  fEnergyMax(0),
  fEnergyLGNbins  (0),
  fEnergyLGMin(0),
  fEnergyLGMax(0),
  fFineNbins(0),
  fFineTmin(0),
  fFineTmax(0),
  fL1PhaseList(0),
  fBadReco(kFALSE),
  fFillHeavyHisto(kFALSE),
  fBadChannelMapArray(),
  fBadChannelMapSet(kFALSE),
  fSetBadChannelMapSource(0),
  fBadChannelFileName(),
  fhcalcEvtTime(0),
  fhEvtTimeHeader(0),
  fhEvtTimeDiff(0),
  fhEventType(0),
  fhTOFT0vsEventNumber(0),
  fhTcellvsTOFT0(0),
  fhTcellvsTOFT0HD(0),
  fhTcellvsSM(0),
  fhEneVsAbsIdHG(0),
  fhEneVsAbsIdLG(0),
  fhTimeVsBC(0),
  fhTimeSumSq(),
  fhTimeEnt(),
  fhTimeSum(),
  fhTimeLGSumSq(),
  fhTimeLGEnt(),
  fhTimeLGSum(),
  fhAllAverageBC(),
  fhAllAverageLGBC(),
  fhRefRuns(0),
  fhTimeDsup(),
  fhTimeDsupBC(),
  fhTimeDsupLG(),
  fhTimeDsupLGBC(),
  fhRawTimeVsIdBC(),
  fhRawTimeSumBC(),
  fhRawTimeEntriesBC(),
  fhRawTimeSumSqBC(),
  fhRawTimeVsIdLGBC(),
  fhRawTimeSumLGBC(),
  fhRawTimeEntriesLGBC(),
  fhRawTimeSumSqLGBC(),
  fhRawCorrTimeVsIdBC(),
  fhRawCorrTimeVsIdLGBC(),
  fhTimeVsIdBC(),
  fhTimeVsIdLGBC()
{
  for(Int_t i = 0; i < kNBCmask; i++) 
  {
    fhAllAverageBC[i]=0;
    fhAllAverageLGBC[i]=0;

    fhTimeSumSq[i]=0;
    fhTimeEnt[i]=0;    
    fhTimeSum[i]=0;

    fhTimeLGSumSq[i]=0;
    fhTimeLGEnt[i]=0;
    fhTimeLGSum[i]=0;

    fhRawTimeVsIdBC[i]=0;
    fhRawTimeSumBC[i]=0;
    fhRawTimeEntriesBC[i]=0;
    fhRawTimeSumSqBC[i]=0;

    fhRawTimeVsIdLGBC[i]=0;
    fhRawTimeSumLGBC[i]=0;
    fhRawTimeEntriesLGBC[i]=0;
    fhRawTimeSumSqLGBC[i]=0;
    fhRawCorrTimeVsIdBC[i]=0;
    fhRawCorrTimeVsIdLGBC[i]=0;
    fhTimeVsIdBC[i]=0;
    fhTimeVsIdLGBC[i]=0;

  }

  //set default cuts for calibration and geometry name
  SetDefaultCuts();

  //T0 TOF time
  PrepareTOFT0maker();

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
  
} // End ctor

//_____________________________________________________________________
/// HKD Move from constructor
/// Use aliTOFT0maker to get proper T0
/// Look the proper source to have more information
/// Modified July 2, 2010 - HKD to take into account
/// the changes in ALiTOFT0maker
//void AliAnalysisTaskEMCALTimeCalib::LocalInit()
//{
//  AliDebug(1,"AliAnalysisTaskEMCALTimeCalib::LocalInit()");
//}

/// Load reference Histograms (for one period) from file
//_____________________________________________________________________
void AliAnalysisTaskEMCALTimeCalib::LoadReferenceHistos()
{
  if(fReferenceFileName.Length()!=0){
    TFile *myFile = TFile::Open(fReferenceFileName.Data());
    AliInfo(Form("Reference file: %s, pointer %p",fReferenceFileName.Data(),myFile));
    if(myFile==0x0) {
      AliFatal("*** NO REFERENCE FILE");
    } else {
	AliDebug(1,"*** OK TFILE");
	// connect ref run here
	for(Int_t i = 0; i < kNBCmask; i++)
	  {
	    fhAllAverageBC[i]=(TH1F*) myFile->Get(Form("hAllTimeAvBC%d",i));
	    if(fhAllAverageBC[i]==0x0) AliFatal(Form("Reference histogram for BC%d does not exist",i));
	    if(fhAllAverageBC[i]->GetEntries()==0)AliWarning(Form("fhAllAverageLGBC[%d]->GetEntries() = 0",i));
	    fhAllAverageLGBC[i]=(TH1F*) myFile->Get(Form("hAllTimeAvLGBC%d",i));
	    if(fhAllAverageLGBC[i]==0x0) AliFatal(Form("Reference LG histogram for BC%d does not exist",i));
	    if(fhAllAverageLGBC[i]->GetEntries()==0)AliFatal(Form("fhAllAverageLGBC[%d]->GetEntries() = 0",i));
	  }
	
	AliDebug(1,Form("hAllAverage entries BC0 %d", (Int_t)fhAllAverageBC[0]->GetEntries() ));
	AliDebug(1,Form("hAllAverage entries BC2 %d",(Int_t)fhAllAverageBC[2]->GetEntries() ));
	AliDebug(1,Form("hAllAverageLG entries BC0 %d", (Int_t)fhAllAverageLGBC[0]->GetEntries() ));
	AliDebug(1,Form("hAllAverageLG entries BC2 %d",(Int_t)fhAllAverageLGBC[2]->GetEntries() ));
	
    }
  } else { //end of reference file is provided
    AliFatal("You require to load reference histos from file but FILENAME is not provided");
  }
} // End of AliAnalysisTaskEMCALTimeCalib::LoadReferenceHistos()

/// Load reference Histograms (run-by-run in one period) from file into memory
/// This method should be called at the beginning of processing only once
//_____________________________________________________________________
void AliAnalysisTaskEMCALTimeCalib::LoadReferenceRunByRunHistos()
{
  // connect ref run here
  if(fReferenceRunByRunFileName.Length()!=0){
    TFile *referenceFile = TFile::Open(fReferenceRunByRunFileName.Data());
    if(referenceFile==0x0) {
      AliFatal("*** NO REFERENCE R-B-R FILE");
      return;
    } else {
      AliInfo(Form("Reference R-b-R file: %s, pointer %p",fReferenceRunByRunFileName.Data(),referenceFile));

      //load L1 phases to memory
      fL1PhaseList=new TObjArray(referenceFile->GetNkeys());
      TIter next(referenceFile->GetListOfKeys());
      TKey *key;
      while ((key=(TKey*)next())) {
	fL1PhaseList->AddLast((TH1F*)referenceFile->Get(key->GetName()) );
	//printf("key: %s points to an object of class: %s at %dn",key->GetName(),key->GetClassName(),key->GetSeekKey());
      }
    }
  } else { //reference file is not provided
    AliFatal("You require to load reference run-by-run histos from file but FILENAME is not provided");
    return;
  }
} // End of AliAnalysisTaskEMCALTimeCalib::LoadReferenceRunByRunHistos()

/// Load reference histogram with L1 phases for given run
/// This method should be called per run
////_____________________________________________________________________
void AliAnalysisTaskEMCALTimeCalib::SetL1PhaseReferenceForGivenRun()
{
  fhRefRuns=NULL;
  if(!fL1PhaseList) {
    AliFatal("Array with reference L1 phase histograms do not exist in memory");
    return;
  }
  if(fRunNumber<0) {
    AliFatal("Negative run number");
    return;
  }

  fhRefRuns=(TH1C*)fL1PhaseList->FindObject(Form("h%d",fRunNumber));
  if(fhRefRuns==0x0){
    AliError(Form("Reference histogram for run %d does not exist. Use Default",fRunNumber));
    fhRefRuns=(TH1C*)fL1PhaseList->FindObject("h0");
  }
  if(fhRefRuns==0x0) {
    AliFatal(Form("No default histogram with L1 phases! Add default histogram to file %s!!!",fReferenceRunByRunFileName.Data()));
    return;
  }

  AliDebug(1,Form("Reference R-b-R histo %p, list %p, run number %d",fhRefRuns,fL1PhaseList,fRunNumber));
  if(fhRefRuns->GetEntries()==0)AliWarning("fhRefRuns->GetEntries() = 0");
  AliDebug(1,Form("hRefRuns entries %d", (Int_t)fhRefRuns->GetEntries() )); 
}

//_____________________________________________________________________
/// Connect ESD or AOD here
/// Called when run is changed
void AliAnalysisTaskEMCALTimeCalib::NotifyRun()
{
  AliDebug(1,"AnalysisTaskEMCalTimeCalib::NotifyRun()");
  AliDebug(2,Form("Notify(): EMCal geometry: fgeom = %p, fGeometryName=%s\n ",fgeom,fGeometryName.Data()));

  if (!InputEvent())
  {
    AliFatal("ERROR: InputEvent not set");
    return;
  }
  else AliDebug(1,"Good, InputEvent set");

  //  AliInfo(Form("NotifyRun, fCurrentRunnumber %d",fCurrentRunNumber));
  fRunNumber = InputEvent()->GetRunNumber();
  AliDebug(1,Form("RunNumber %d", fRunNumber));

  // Init EMCAL geometry 
  if (!fgeom) SetEMCalGeometry();
  //Init EMCAL geometry done

  //set L1 phases for current run
  if(fReferenceRunByRunFileName.Length()!=0)
    SetL1PhaseReferenceForGivenRun();

  // set bad channel map
  if(!fBadChannelMapSet && fSetBadChannelMapSource>0) LoadBadChannelMap();

  return;
}

//_____________________________________________________________________
/// Set the EMCal Geometry
Bool_t AliAnalysisTaskEMCALTimeCalib::SetEMCalGeometry()
{
  AliDebug(1,"AliAnalysisTaskEMCALTimeCalib::SetEMCalGeometry()");
  if(fGeometryName.Length()==0){
    fgeom=AliEMCALGeometry::GetInstanceFromRunNumber(fRunNumber);
    AliInfo(Form("Get EMCAL geometry name <%s> for run %d",fgeom->GetName(),fRunNumber));
  } else {
    fgeom = AliEMCALGeometry::GetInstance(fGeometryName.Data());
    AliInfo(Form("Set EMCAL geometry name to <%s>",fGeometryName.Data()));
  }

  if (!fgeom){
    AliWarning("Make sure the EMCal geometry is set properly !");
  } else {
    AliDebug(1,Form("EMCal geometry properly set: fGeom = %p, fGeometryName=%s",fgeom,fGeometryName.Data()));
  }

  return kTRUE;
}

//_____________________________________________________________________
/// Get T0 time from TOF
void AliAnalysisTaskEMCALTimeCalib::PrepareTOFT0maker()
{
  //method under development
  AliInfo(Form("<D> -- Run # = %d", fRunNumber));
  AliInfo("prepare TOFT0maker!!");
  //cout<<"Run "<< fRunNumber<<" in TOFT0maker"<<endl;


  AliCDBManager * cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");
  cdb->SetRun(fRunNumber);
  
  AliESDpid *extPID=new AliESDpid();

  // Wonder if some have to be declared as private variables??
  // AliESDpid *extPID = new AliESDpid();
  // AliTOFcalib * tofCalib = new AliTOFcalib();
  // tofCalib->SetCalibrateTOFsignal(kTRUE);
  // tofCalib->Init();
  
  fTOFmaker = new AliTOFT0maker(extPID);
  fTOFmaker->SetTimeResolution(115.0); // if you want set the TOF res
  // fTOFmaker = new AliTOFT0maker(extPID,tofCalib);
  // fTOFmaker->SetTimeResolution(130.0);

  //cout<<"extPID "<<extPID<<" fTOFmaker "<<fTOFmaker<<endl;
  
}// End PrepareTOFT0maker

//________________________________________________________________________
/// Create histograms
/// Called once
void AliAnalysisTaskEMCALTimeCalib::UserCreateOutputObjects()
{
  AliDebug(1,"AliAnalysisTaskEMCALTimeCalib::UserCreateOutputObjects()");

  const Int_t nChannels = 17664;
  //book histograms
  if(fFillHeavyHisto){
    fhcalcEvtTime = new TH1F("fhcalcEvtTime","calculated event time from T0",fFineNbins, fFineNbins,fFineTmax);
    fhcalcEvtTime->GetXaxis()->SetTitle("T ");
    fhcalcEvtTime->GetYaxis()->SetTitle("Counts (a.u.)");
  
    fhEvtTimeHeader = new TH1F("fhEvtTimeHeader","event time from header",fFineNbins, fFineNbins,fFineTmax);
    fhEvtTimeHeader->GetXaxis()->SetTitle("T ");
    fhEvtTimeHeader->GetYaxis()->SetTitle("Counts (a.u.)");

    fhEvtTimeDiff = new TH1F("fhEvtTimeDiff","event time difference",fFineNbins, fFineNbins,fFineTmax);
    fhEvtTimeDiff->GetXaxis()->SetTitle("#Delta T ");
    fhEvtTimeDiff->GetYaxis()->SetTitle("Counts (a.u.)");
  }

  fhEventType = new TH1F("fhEventType","event type",10, 0.,10.);
  //fhEventType ->GetXaxis()->SetTitle("Type ");
  fhEventType->GetXaxis()->SetBinLabel(1 ,"1=No ESD");
  fhEventType->GetXaxis()->SetBinLabel(2 ,"2=Pileup");
  fhEventType->GetXaxis()->SetBinLabel(3 ,"3=No Trigger");
  fhEventType->GetXaxis()->SetBinLabel(4 ,"4=Evt Type != 7");
  fhEventType->GetXaxis()->SetBinLabel(5 ,"5=INT7,8");
  fhEventType->GetXaxis()->SetBinLabel(6 ,"6=EMC7,8");
  fhEventType->GetXaxis()->SetBinLabel(7 ,"7=L1 EMCal");
  fhEventType->GetXaxis()->SetBinLabel(8 ,"8=DMC7,8");
  fhEventType->GetXaxis()->SetBinLabel(9 ,"9=L1 DCal");

  
  fhEventType ->GetYaxis()->SetTitle("Counts (a.u.)");
  if(fFillHeavyHisto){
    fhTcellvsTOFT0 = new TH2F("hTcellvsTOFT0", " T_cell vs TOFT0", 500,-600.0,+400.0,fRawTimeNbins,fRawTimeMin,fRawTimeMax);
    fhTcellvsTOFT0HD = new TH2F("hTcellvsTOFT0HD", " T_cell vs TOFT0,HighEnergy", 500,-600.0,+400.0,4*fRawTimeNbins,fRawTimeMin,fRawTimeMax);
  }
  fhTcellvsSM = new TH2F("hTcellvsSM", " T_cell vs SM", (Int_t)kNSM,0,(Double_t)kNSM,(Int_t)(fRawTimeNbins/2),fRawTimeMin,fRawTimeMax);

  if(fFillHeavyHisto){
    fhEneVsAbsIdHG = new TH2F("fhEneVsAbsIdHG", "energy vs ID for HG",1000,0,18000,200,0,10);
    fhEneVsAbsIdLG = new TH2F("fhEneVsAbsIdLG", "energy vs ID for LG",1000,0,18000,200,0,40);
  }

  for (Int_t i = 0; i < kNBCmask ;  i++)
  {
    //already after correction
    //high gain
    fhTimeSumSq[i] = new TH1F(Form("hTimeSumSq%d", i),
			      Form("cell Sum Square time HG, BC %d ", i),
			      nChannels,0.,(Double_t)nChannels);
    fhTimeSumSq[i]->SetYTitle("Sum Sq Time ");
    fhTimeSumSq[i]->SetXTitle("AbsId");
    
    fhTimeSum[i] = new TH1F(Form("hTimeSum%d", i),
			    Form("cell Sum  time HG, BC %d ", i),
			    nChannels,0.,(Double_t)nChannels);
    fhTimeSum[i]->SetYTitle("Sum  Time ");
    fhTimeSum[i]->SetXTitle("AbsId");
    
    fhTimeEnt[i] = new TH1F(Form("hTimeEnt%d", i),
			    Form("cell Entries HG, BC %d ", i),
			    nChannels,0.,(Double_t)nChannels);
    fhTimeEnt[i]->SetYTitle("Entries for Time ");
    fhTimeEnt[i]->SetXTitle("AbsId");
    
    //low gain 
    fhTimeLGSumSq[i] = new TH1F(Form("hTimeLGSumSq%d", i),
			      Form("cell Sum Square time LG, BC %d ", i),
			      nChannels,0.,(Double_t)nChannels);
    fhTimeLGSumSq[i]->SetYTitle("Sum Sq Time ");
    fhTimeLGSumSq[i]->SetXTitle("AbsId");
    
    fhTimeLGSum[i] = new TH1F(Form("hTimeLGSum%d", i),
			    Form("cell Sum time LG, BC %d ", i),
			    nChannels,0.,(Double_t)nChannels);
    fhTimeLGSum[i]->SetYTitle("Sum  Time ");
    fhTimeLGSum[i]->SetXTitle("AbsId");
    
    fhTimeLGEnt[i] = new TH1F(Form("hTimeLGEnt%d", i),
			    Form("cell Entries LG, BC %d ", i),
			    nChannels,0.,(Double_t)nChannels);
    fhTimeLGEnt[i]->SetYTitle("Entries for Time ");
    fhTimeLGEnt[i]->SetXTitle("AbsId");

    //raw time histograms
    //high gain
    if(fFillHeavyHisto){
      fhRawTimeVsIdBC[i] = new TH2F(Form("RawTimeVsIdBC%d", i),
				    Form("cell raw time vs ID for high gain BC %d ", i),
				    nChannels,0.,(Double_t)nChannels,fRawTimeNbins,fRawTimeMin,fRawTimeMax);
      fhRawTimeVsIdBC[i]->SetXTitle("AbsId");
      fhRawTimeVsIdBC[i]->SetYTitle("Time");
    }

    fhRawTimeSumBC[i] = new TH1F(Form("RawTimeSumBC%d", i),
				 Form("sum of cell raw time for high gain BC %d ", i),
				 nChannels,0.,(Double_t)nChannels);
    fhRawTimeSumBC[i]->SetXTitle("AbsId");
    fhRawTimeSumBC[i]->SetYTitle("Sum Time");

    fhRawTimeEntriesBC[i] = new TH1F(Form("RawTimeEntriesBC%d", i),
				     Form("No. entries of cells raw time for high gain BC %d ", i),
				     nChannels,0.,(Double_t)nChannels);
    fhRawTimeEntriesBC[i]->SetXTitle("AbsId");
    fhRawTimeEntriesBC[i]->SetYTitle("Entries for Time ");

    fhRawTimeSumSqBC[i] = new TH1F(Form("RawTimeSumSqBC%d", i),
				   Form("sum of (cell raw time)^2 for high gain BC %d ", i),
				   nChannels,0.,(Double_t)nChannels);
    fhRawTimeSumSqBC[i]->SetXTitle("AbsId");
    fhRawTimeSumSqBC[i]->SetYTitle("Sum Sq Time");

    //low gain
    if(fFillHeavyHisto){
      fhRawTimeVsIdLGBC[i] = new TH2F(Form("RawTimeVsIdLGBC%d", i),
				      Form("cell raw time vs ID for low gain BC %d ", i),
				      nChannels,0.,(Double_t)nChannels,fRawTimeNbins,fRawTimeMin,fRawTimeMax);
      fhRawTimeVsIdLGBC[i]->SetXTitle("AbsId");
      fhRawTimeVsIdLGBC[i]->SetYTitle("Time");
    }

    fhRawTimeSumLGBC[i] = new TH1F(Form("RawTimeSumLGBC%d", i),
				 Form("sum of cell raw time for low gain BC %d ", i),
				 nChannels,0.,(Double_t)nChannels);
    fhRawTimeSumLGBC[i]->SetXTitle("AbsId");
    fhRawTimeSumLGBC[i]->SetYTitle("Sum Time");

    fhRawTimeEntriesLGBC[i] = new TH1F(Form("RawTimeEntriesLGBC%d", i),
				     Form("No. entries of cells raw time for low gain BC %d ", i),
				     nChannels,0.,(Double_t)nChannels);
    fhRawTimeEntriesLGBC[i]->SetXTitle("AbsId");
    fhRawTimeEntriesLGBC[i]->SetYTitle("Entries for Time ");

    fhRawTimeSumSqLGBC[i] = new TH1F(Form("RawTimeSumSqLGBC%d", i),
				   Form("sum of (cell raw time)^2 for low gain BC %d ", i),
				     nChannels,0.,(Double_t)nChannels);
    fhRawTimeSumSqLGBC[i]->SetXTitle("AbsId");
    fhRawTimeSumSqLGBC[i]->SetYTitle("Sum Sq Time");

    //histograms with corrected raw time for L1 shift and 100ns
    if(fBadReco && fFillHeavyHisto){
      fhRawCorrTimeVsIdBC[i] = new TH2F(Form("RawCorrTimeVsIdBC%d", i),
					Form("cell L1 shift and 100ns corrected raw time vs ID for high gain BC %d ", i),
					nChannels,0.,(Double_t)nChannels,fPassTimeNbins,fPassTimeMin,fPassTimeMax);
      fhRawCorrTimeVsIdBC[i]->SetXTitle("AbsId");
      fhRawCorrTimeVsIdBC[i]->SetYTitle("Time");
      
      fhRawCorrTimeVsIdLGBC[i] = new TH2F(Form("RawCorrTimeVsIdLGBC%d", i),
					  Form("cell L1 shift and 100ns corrected raw time vs ID for low gain BC %d ", i),
					  nChannels,0.,(Double_t)nChannels,fPassTimeNbins,fPassTimeMin,fPassTimeMax);
      fhRawCorrTimeVsIdLGBC[i]->SetXTitle("AbsId");
      fhRawCorrTimeVsIdLGBC[i]->SetYTitle("Time");
    }

    //histograms with corrected raw time for L1 shift and 100ns + new L1 phase
    if(fReferenceRunByRunFileName.Length()!=0 && fFillHeavyHisto){
      fhTimeVsIdBC[i] = new TH2F(Form("TimeVsIdBC%d", i),
				 Form("cell time corrected for L1 shift, 100ns and L1 phase vs ID for high gain BC %d ", i),
				 nChannels,0.,(Double_t)nChannels,fPassTimeNbins,fPassTimeMin,fPassTimeMax);
      fhTimeVsIdBC[i]->SetXTitle("AbsId");
      fhTimeVsIdBC[i]->SetYTitle("Time");
      
      fhTimeVsIdLGBC[i] = new TH2F(Form("TimeVsIdLGBC%d", i),
				   Form("cell time corrected for L1 shift, 100ns and L1 phase vs ID for low gain BC %d ", i),
				   nChannels,0.,(Double_t)nChannels,fPassTimeNbins,fPassTimeMin,fPassTimeMax);
      fhTimeVsIdLGBC[i]->SetXTitle("AbsId");
      fhTimeVsIdLGBC[i]->SetYTitle("Time");
    }

    for (Int_t j = 0; j < kNSM ;  j++) 
    {
      //High gain
      //fhTimeDsupBC[j][i]= new TH2F(Form("SupMod%dBC%d",j,i), Form("SupMod %d time_vs_E  BC %d",j,i),500,0.0,20.0,2200,-350.0,750.0);
      fhTimeDsupBC[j][i]= new TH2F(Form("SupMod%dBC%d",j,i), Form("SupMod %d time_vs_E, high gain, BC %d",j,i),fEnergyNbins,fEnergyMin,fEnergyMax,fPassTimeNbins,fPassTimeMin,fPassTimeMax);
      fhTimeDsupBC[j][i]->SetYTitle(" Time (ns) "); 
      fhTimeDsupBC[j][i]->SetXTitle(" E (GeV) "); 

      //low gain
      fhTimeDsupLGBC[j][i]= new TH2F(Form("SupMod%dBC%dLG",j,i), Form("SupMod %d time_vs_E, low gain, BC %d",j,i),fEnergyLGNbins,fEnergyLGMin,fEnergyLGMax,fPassTimeNbins,fPassTimeMin,fPassTimeMax);
      fhTimeDsupLGBC[j][i]->SetYTitle(" Time (ns) "); 
      fhTimeDsupLGBC[j][i]->SetXTitle(" E (GeV) "); 
    }
  }

  for (Int_t jj = 0; jj < kNSM ;  jj++) 
  {
    //high gain
    fhTimeDsup[jj] =  new TH2F(Form("SupMod%d",jj), Form("SupMod %d time_vs_E, high gain",jj),fEnergyNbins,fEnergyMin,fEnergyMax,fPassTimeNbins,fPassTimeMin,fPassTimeMax);
    fhTimeDsup[jj]->SetYTitle(" Time (ns) "); 
    fhTimeDsup[jj]->SetXTitle(" E (GeV) "); 

    //low gain
    fhTimeDsupLG[jj] =  new TH2F(Form("SupMod%dLG",jj), Form("SupMod %d time_vs_E, low gain ",jj),fEnergyLGNbins,fEnergyLGMin,fEnergyLGMax,fPassTimeNbins,fPassTimeMin,fPassTimeMax);
    fhTimeDsupLG[jj]->SetYTitle(" Time (ns) "); 
    fhTimeDsupLG[jj]->SetXTitle(" E (GeV) "); 
  }
  
  fhTimeVsBC = new TH2F("TimeVsBC"," SupMod time_vs_BC ", 4001,-0.5,4000.5,(Int_t)(fRawTimeNbins/2.),fRawTimeMin,fRawTimeMax); 
  

  //add histos to list
  fOutputList = new TList();
  fOutputList->Add(fhEventType);
  if(fFillHeavyHisto){
    fOutputList->Add(fhcalcEvtTime);
    fOutputList->Add(fhEvtTimeHeader);
    fOutputList->Add(fhEvtTimeDiff);

    fOutputList->Add(fhTcellvsTOFT0);
    fOutputList->Add(fhTcellvsTOFT0HD);
    fOutputList->Add(fhEneVsAbsIdHG);
    fOutputList->Add(fhEneVsAbsIdLG);
  }
  fOutputList->Add(fhTcellvsSM);

  for (Int_t i = 0; i < kNBCmask ;  i++) 
  {
    fOutputList->Add(fhTimeSumSq[i]);
    fOutputList->Add(fhTimeEnt[i]);
    fOutputList->Add(fhTimeSum[i]);

    fOutputList->Add(fhTimeLGSumSq[i]);
    fOutputList->Add(fhTimeLGEnt[i]);
    fOutputList->Add(fhTimeLGSum[i]);

    if(fFillHeavyHisto) {
      fOutputList->Add(fhRawTimeVsIdBC[i]);
      fOutputList->Add(fhRawTimeVsIdLGBC[i]);
    }

    fOutputList->Add(fhRawTimeSumBC[i]);
    fOutputList->Add(fhRawTimeEntriesBC[i]);
    fOutputList->Add(fhRawTimeSumSqBC[i]);

    fOutputList->Add(fhRawTimeSumLGBC[i]);
    fOutputList->Add(fhRawTimeEntriesLGBC[i]);
    fOutputList->Add(fhRawTimeSumSqLGBC[i]);

    if(fBadReco && fFillHeavyHisto) {
      fOutputList->Add(fhRawCorrTimeVsIdBC[i]);
      fOutputList->Add(fhRawCorrTimeVsIdLGBC[i]);
    }
    if(fReferenceRunByRunFileName.Length()!=0 && fFillHeavyHisto) {
      fOutputList->Add(fhTimeVsIdBC[i]);
      fOutputList->Add(fhTimeVsIdLGBC[i]);
    }

    for (Int_t j = 0; j < kNSM ;  j++){
      fOutputList->Add(fhTimeDsupBC[j][i]);
      fOutputList->Add(fhTimeDsupLGBC[j][i]);
    }
  }
  
  for (Int_t j = 0; j < kNSM ;  j++)
  {
    fOutputList->Add(fhTimeDsup[j]);
    fOutputList->Add(fhTimeDsupLG[j]);
  }
	
  fOutputList->Add(fhTimeVsBC);
  
  fOutputList->SetOwner(kTRUE);
  PostData(1,fOutputList);

  
} // End of AliAnalysisTaskEMCALTimeCalib::UserCreateOuputObjects()

//________________________________________________________________________
/// Main loop executed for each event
void AliAnalysisTaskEMCALTimeCalib::UserExec(Option_t *)
{
  // Called for each event
  AliDebug(2,Form("UserExec: EMCal geometry: fgeom = %p fGeometryName %s",fgeom,fGeometryName.Data()));
  AliVEvent   *event = InputEvent();
  //cout<<"T0TOF "<<event->GetT0TOF()<<endl;//bad idea
  //cout<< fEvent->GetTOFHeader()->GetDefaultEventTimeVal()<<endl;
  AliDebug(2,Form("TOF time from header %f ps",event->GetTOFHeader()->GetDefaultEventTimeVal()));
  if(fFillHeavyHisto) fhEvtTimeHeader->Fill(event->GetTOFHeader()->GetDefaultEventTimeVal());

  //fEvent = dynamic_cast<AliESDEvent*>(event);
  if (!event) {
    AliError("ESD not available, exit");
    fhEventType->Fill(0.5);
    return;
  }
  
  if(fPileupFromSPD==kTRUE){
    if(event->IsPileupFromSPD(3,0.8,3.,2.,5.)){
      AliDebug(1,"Event: PileUp skip.");
      fhEventType->Fill(1.5);
      return;
    }	
  }

  TString triggerclasses = event->GetFiredTriggerClasses();
  if(triggerclasses=="") {
    fhEventType->Fill(2.5);
    return;
  }

  Int_t eventType = ((AliVHeader*)event->GetHeader())->GetEventType();
  // physics events eventType=7, select only those
  AliDebug(1,Form("Triggerclasses %s, eventType %d",triggerclasses.Data(),eventType));
  if(eventType != 7) {
    fhEventType->Fill(3.5);
    return;
  }
  
  // Check trigger
  Bool_t bMB  = kFALSE;
  Bool_t bL0  = kFALSE;
  Bool_t bL1G = kFALSE;
  Bool_t bL1J = kFALSE;
  Bool_t bDL0  = kFALSE;
  Bool_t bDL1G = kFALSE;
  Bool_t bDL1J = kFALSE;
  
  if(triggerclasses.Contains("CINT7-B-NOPF-ALLNOTRD") ||
     triggerclasses.Contains("CINT7-I-NOPF-ALLNOTRD") ||
     triggerclasses.Contains("CINT1-I-NOPF-ALLNOTRD") || 
     triggerclasses.Contains("CINT1-B-NOPF-ALLNOTRD") ||
     triggerclasses.Contains("CINT8") ||
     triggerclasses.Contains("CINT7") ||
     triggerclasses.Contains("CPBI2_B1-B-NOPF-ALLNOTRD") )   bMB  = kTRUE;
  
  if(triggerclasses.Contains("CEMC7-B-NOPF-CENTNOTRD") || 
     triggerclasses.Contains("CEMC1-B-NOPF-CENTNOTRD") ||
     triggerclasses.Contains("CEMC7") ||
     triggerclasses.Contains("CEMC8") ||
     triggerclasses.Contains("CEMC8-B-NOPF-CENTNOTRD")   )   bL0  = kTRUE;

  if(triggerclasses.Contains("CDMC7-B-NOPF-CENTNOTRD") || 
     triggerclasses.Contains("CDMC1-B-NOPF-CENTNOTRD") ||
     triggerclasses.Contains("CDMC7") ||
     triggerclasses.Contains("CDMC8") ||
     triggerclasses.Contains("CDMC8-B-NOPF-CENTNOTRD")   )   bDL0  = kTRUE;
  
  if(triggerclasses.Contains("CEMC7EG1-B-NOPF-CENTNOTRD") ||
     triggerclasses.Contains("CEMC7EG2-B-NOPF-CENTNOTRD") ||
     triggerclasses.Contains("CEMC8EG1-B-NOPF-CENTNOTRD") ||
     triggerclasses.Contains("CEMC8EGA") ||
     triggerclasses.Contains("CEMC7EGA") ||
     triggerclasses.Contains("CEMC7EG1-B") ||
     triggerclasses.Contains("CEMC7EG2-B") ||
     triggerclasses.Contains("CPBI2EGA")                 )   bL1G = kTRUE;
 
  if(triggerclasses.Contains("CDMC7DG1-B-NOPF-CENTNOTRD") ||
     triggerclasses.Contains("CDMC7DG2-B-NOPF-CENTNOTRD") ||
     triggerclasses.Contains("CDMC8DG1-B-NOPF-CENTNOTRD") ||
     triggerclasses.Contains("CDMC8DGA") ||
     triggerclasses.Contains("CDMC7DGA") ||
     triggerclasses.Contains("CDMC7DG1-B") ||
     triggerclasses.Contains("CDMC7DG2-B") ||
     triggerclasses.Contains("CPBI2DGA")                 )   bDL1G = kTRUE;
    
  if(triggerclasses.Contains("CEMC7EJ1-B-NOPF-CENTNOTRD") ||
     triggerclasses.Contains("CEMC7EJ2-B-NOPF-CENTNOTRD") ||
     triggerclasses.Contains("CEMC8EJ1-B-NOPF-CENTNOTRD") ||
     triggerclasses.Contains("CEMC7EJE") ||
     triggerclasses.Contains("CEMC8EJE") ||
     triggerclasses.Contains("CEMC7EJ1-B") ||
     triggerclasses.Contains("CEMC7EJ2-B") ||
     triggerclasses.Contains("CPBI2EJE")                 )   bL1J = kTRUE;

  if(triggerclasses.Contains("CDMC7DJ1-B-NOPF-CENTNOTRD") ||
     triggerclasses.Contains("CDMC7DJ2-B-NOPF-CENTNOTRD") ||
     triggerclasses.Contains("CDMC8DJ1-B-NOPF-CENTNOTRD") ||
     triggerclasses.Contains("CDMC7DJE") ||
     triggerclasses.Contains("CDMC8DJE") ||
     triggerclasses.Contains("CDMC7DJ1-B") ||
     triggerclasses.Contains("CDMC7DJ2-B") ||
     triggerclasses.Contains("CPBI2DJE")                 )   bDL1J = kTRUE;
    
  if( bMB ){ fhEventType->Fill(4.5);}//INT7,8
  if( bL0 ){ fhEventType->Fill(5.5);}//EMC7,EMC8
  if( bL1G || bL1J ){ fhEventType->Fill(6.5);}//L1 EMCal
  if( bDL0 ){ fhEventType->Fill(7.5);}//DMC7,DMC8
  if( bDL1G || bDL1J ){ fhEventType->Fill(8.5);}//L1 DCal
  
  //  if(bL1G || bL1J ||  bL0){

// Prepare TOFT0 maker at the beginning of a run
//  if (event->GetRunNumber() != fRunNumber){
//    AliInfo(Form("Runno per event %d",event->GetRunNumber()));
//    fRunNumber = event->GetRunNumber();
//    //	PrepareTOFT0maker();
//    //  	cout<<"tofT0maker per run"<<fRunNumber<<endl;
//  }// fi Check if run number has changed
  
  // --- Use of AliTOFT0maker
  Double_t calcolot0=0.0;
  if(!AODEvent()){
    Double_t* timeTOFtable;
    timeTOFtable=fTOFmaker->ComputeT0TOF(dynamic_cast<AliESDEvent*>(event));
    AliDebug(2,Form("TOF time %f ps, resolution %f ps, tracks at TOF %f/used %f",timeTOFtable[0],timeTOFtable[1],timeTOFtable[3],timeTOFtable[7]));
    //cout<<"event time "<<timeTOFtable[0]<<" resolution "<<timeTOFtable[1]<<"ps av. ev. time "<<timeTOFtable[2]<<" trks at TOF "<<timeTOFtable[3]<<" calc evnt time "<<timeTOFtable[4]<<" resolution "<<timeTOFtable[5]<<" tracks used "<<timeTOFtable[7]<<endl;
    calcolot0=timeTOFtable[0];
  }

  if(fFillHeavyHisto) {
    fhcalcEvtTime->Fill(calcolot0);
    if(calcolot0 != 0 && event->GetTOFHeader()->GetDefaultEventTimeVal() != 0 )
      fhEvtTimeDiff->Fill(calcolot0-event->GetTOFHeader()->GetDefaultEventTimeVal());
  }

  TRefArray* caloClusters = new TRefArray();
  event->GetEMCALClusters(caloClusters);
  //           	cout << " ###########Bunch Cross nb  = " << event->GetBunchCrossNumber() << endl;
  
  Int_t BunchCrossNumber =event->GetBunchCrossNumber(); 
  
  Float_t offset=0.;
  Float_t offsetPerSM=0.;
  Int_t L1phaseshift=0;
  Int_t L1phase=0;
  Int_t L1shiftOffset=0;

  Int_t nBC = 0;
  nBC = BunchCrossNumber%4;
  //Int_t nTriggerMask =event->GetTriggerMask();
  //	cout << " nBC " << nBC << " nTriggerMask " << nTriggerMask<< endl;
  Float_t timeBCoffset = 0.; //manual offset
  //	if( nBC%4 ==0 || nBC%4==1) timeBCoffset = 100.; // correction was for LHC11 when BC was not corrected
  
  Int_t nclus = caloClusters->GetEntries();
  AliDebug(1,Form("###########Bunch Cross nb = %d nclus = %d",nBC,nclus ));
  //cout << " ###########Bunch Cross nb  = " << nBC <<" nclus= "<< nclus<< endl;
  //Int_t ntracks = event-> GetNumberOfTracks() ; 

  AliVCaloCells &cells= *(event->GetEMCALCells());//it is cluster independent
  //Variables used plenty in loops
  Int_t nSupMod=-1, nModule=-1;
  Int_t iphi=-1, ieta=-1, nIphi=-1, nIeta=-1;
  Int_t absId=-1;
  Float_t hkdtime=0.0;
  Float_t amp=0.0;
  Bool_t isHighGain=kTRUE;
  Int_t mostEneId=-1;
  Float_t mostEneEn=0.;
  
  for (Int_t icl = 0; icl < nclus; icl++) {
    //ESD and AOD CaloCells carries the same information
    AliVCluster* clus = (AliVCluster*)caloClusters->At(icl);
    if(!AcceptCluster(clus)) continue;

    //cout<<"nCells="<< clus->GetNCells();<<endl;
   
    UShort_t * index = clus->GetCellsAbsId() ;

    // find index of the most energetic cell in cluster
    mostEneEn=0.;
    mostEneId=-1;
    if(fMostEneCellOnly) { 
      for(Int_t i = 0; i < clus->GetNCells() ; i++) {
	absId      = index[i];
	amp        = cells.GetCellAmplitude(absId) ;
	if(amp > mostEneEn){
	  mostEneEn = amp;
	  mostEneId = absId;
	}
      }
    }//works only for fMostEneCellOnly=kTRUE
    
    for(Int_t i = 0; i < clus->GetNCells() ; i++) {
      absId      = index[i]; // or clus->GetCellNumber(i) ;
      if(fMostEneCellOnly && absId != mostEneId) {
	//printf("tr.%s.cl.%d.cell.%d.rejected\n",triggerclasses.Data(),icl,i);
	continue;
      }
      //printf("tr.%s.cl.%d.cell.%d.accepted\n",triggerclasses.Data(),icl,i);
      hkdtime    = cells.GetCellTime(absId) * 1.0e09; // to get ns
      amp        = cells.GetCellAmplitude(absId) ;
      isHighGain = cells.GetCellHighGain(absId);
      //cout<<"cell absID: "<<absId<<" cellTime: "<<hkdtime<<" cellaplit: "<< amp<<endl;	
      // GEOMETRY tranformations
      fgeom->GetCellIndex(absId,  nSupMod, nModule, nIphi, nIeta);
      fgeom->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi,ieta);

      //bad channel check. 0: good channel, 1-5: bad channel
      if(fSetBadChannelMapSource==1){
	if(GetEMCALChannelStatus(nSupMod,ieta,iphi)) continue;//printf("bad\n");
      } else if(fSetBadChannelMapSource==2){
	if(GetEMCALChannelStatus(absId)) continue;//printf("bad\n");
      }

      //main histograms with raw time information 
      if(amp>fMinCellEnergy){
	if(isHighGain){
	  if(fFillHeavyHisto) fhRawTimeVsIdBC[nBC]->Fill(absId,hkdtime);
	  fhRawTimeSumBC[nBC]->Fill(absId,hkdtime);
	  fhRawTimeEntriesBC[nBC]->Fill(absId,1.);
	  fhRawTimeSumSqBC[nBC]->Fill(absId,hkdtime*hkdtime);
	}else{
	  if(fFillHeavyHisto) fhRawTimeVsIdLGBC[nBC]->Fill(absId,hkdtime);
	  fhRawTimeSumLGBC[nBC]->Fill(absId,hkdtime);
	  fhRawTimeEntriesLGBC[nBC]->Fill(absId,1.);
	  fhRawTimeSumSqLGBC[nBC]->Fill(absId,hkdtime*hkdtime);
	}
      }
      //fgeom->PrintCellIndexes(absId);
      //fgeom->PrintCellIndexes(absId,1);

      // other histograms for cross-check      
      CheckCellRCU(nSupMod,ieta,iphi);//SM, column, row

      fhTcellvsSM->Fill(nSupMod,hkdtime);
      if(fFillHeavyHisto) {
	if(isHighGain==kTRUE) {fhEneVsAbsIdHG->Fill(absId,amp);}
	else {fhEneVsAbsIdLG->Fill(absId,amp);}
      }
      fhTimeVsBC->Fill(1.*BunchCrossNumber,hkdtime-timeBCoffset);
      //important remark: We use 'Underflow bin' for absid=0 in OADB for time calibration 
      if(isHighGain==kTRUE){
	if(fhAllAverageBC[nBC]!=0) {//comming from file after the first iteration
	  offset = (Float_t)(fhAllAverageBC[nBC]->GetBinContent(absId));//channel absId=0 has histogram bin=0
	} else if(fReferenceFileName.Length()!=0){//protection against missing reference histogram
	  AliFatal(Form("Reference histogram for BC%d not properly loaded",nBC));
	}
      } else {
	if(fhAllAverageLGBC[nBC]!=0) {//comming from file after the first iteration
	  offset = (Float_t)(fhAllAverageLGBC[nBC]->GetBinContent(absId));//channel absId=0 has histogram bin=0
	} else if(fReferenceFileName.Length()!=0){//protection against missing reference histogram
	  AliFatal(Form("Reference LG histogram for BC%d not properly loaded",nBC));
	}
      }
      //if(offset==0)cout<<"offset 0 in SM "<<nSupMod<<endl;

      // Solution for 2015 data where L1 phase and L1 shift is not correct in data. We need to calibrate run by run.
      // The shift and phase are done per SM (0-19).
      // L1 phase is necessary in run 2. 
      // L1 shift is necessary only for bad reconstructed runs (muon_calo_pass1 lhc15f-m)
      if(fhRefRuns!=0) {//comming from file after the first iteration
	L1phaseshift = (Int_t)(fhRefRuns->GetBinContent(nSupMod));//SM0 = bin0
	
	// to correct for L1 phase
	// this part works for both: muon_calo_pass1 of LHC15n (pp@2.76) and later reconstructions
	// wrong reconstruction done before in run2 
	L1phase = L1phaseshift & 3; //bit operation
	if(nBC >= L1phase)
	  offsetPerSM = (nBC - L1phase)*25;
	else
	  offsetPerSM = (nBC - L1phase + 4)*25;
	
	// to correct for L1 shift
	// this part is only for wrongly reconstructed runs before LHC15n in run2
	if(fBadReco){
	  L1shiftOffset = L1phaseshift>>2; //bit operation
	  L1shiftOffset*=25;
	  //(we subtract it here because we subtract the whole wrong offset later --=+)
	  if(nBC==0 || nBC==1) L1shiftOffset-=100.;//additional shift for muon_calo_pass1 up to lhc15f-m 
	}
      } else if(fReferenceRunByRunFileName.Length()!=0){//protection against missing reference histogram
	AliFatal("Reference histogram run-by-run not properly loaded");
      }
      //end of load additional offset 
      
      //fill the raw time with L1 shift correction and 100ns
      if(fBadReco && fFillHeavyHisto && amp>fMinCellEnergy){
        if(isHighGain){
	  fhRawCorrTimeVsIdBC[nBC]->Fill(absId,hkdtime-L1shiftOffset);
	}else{
	  fhRawCorrTimeVsIdLGBC[nBC]->Fill(absId,hkdtime-L1shiftOffset);
	}
      }

      //fill time after L1 shift correction and 100ns and new L1 phase
      if(fReferenceRunByRunFileName.Length()!=0 && fFillHeavyHisto && amp>fMinCellEnergy){
        if(isHighGain){
          fhTimeVsIdBC[nBC]->Fill(absId,hkdtime-L1shiftOffset-offsetPerSM);
        }else{
          fhTimeVsIdLGBC[nBC]->Fill(absId,hkdtime-L1shiftOffset-offsetPerSM);
        }
      }

      //other control histograms
      if(amp>0.5) {
	if(isHighGain){				
	  fhTimeDsup[nSupMod]->Fill(amp,hkdtime-offset-offsetPerSM-L1shiftOffset);
	  fhTimeDsupBC[nSupMod][nBC]->Fill(amp,hkdtime-offset-offsetPerSM-L1shiftOffset);
	}else{
	  fhTimeDsupLG[nSupMod]->Fill(amp,hkdtime-offset-offsetPerSM-L1shiftOffset);
	  fhTimeDsupLGBC[nSupMod][nBC]->Fill(amp,hkdtime-offset-offsetPerSM-L1shiftOffset);
	}
      }
      
      if(fFillHeavyHisto) {
	if(amp>0.9) {
	  fhTcellvsTOFT0HD->Fill(calcolot0, hkdtime);
	}
	fhTcellvsTOFT0->Fill(calcolot0, hkdtime-offset-offsetPerSM-L1shiftOffset);
      }

      hkdtime = hkdtime-timeBCoffset;//time corrected by manual offset (default=0)
      Float_t hkdtimecorr;
      hkdtimecorr= hkdtime-offset-offsetPerSM-L1shiftOffset;//time after first iteration

      //main histograms after the first itereation for calibration constants
      //if(hkdtimecorr>=-20. && hkdtimecorr<=20. && amp>0.9 ) {
      if(hkdtimecorr>=fMinTime && hkdtimecorr<=fMaxTime && amp>fMinCellEnergy ) {
	// per cell
//	Float_t entriesTime=fhTimeEnt[nBC]->GetBinContent(absId)+1;
//	Float_t sumTimeSq=(fhTimeSumSq[nBC]->GetBinContent(absId)+(hkdtime*hkdtime));
//	Float_t sumTime=(fhTimeSum[nBC]->GetBinContent(absId)+hkdtime);
//	
//	fhTimeEnt[nBC]->SetBinContent(absId,entriesTime);
//	fhTimeSumSq[nBC]->SetBinContent(absId,sumTimeSq);
//	fhTimeSum[nBC]->SetBinContent(absId,sumTime);

        //correction in 2015 for wrong L1 phase and L1 shift
	hkdtime = hkdtime - offsetPerSM - L1shiftOffset;

	if(isHighGain){
	  fhTimeEnt[nBC]->Fill(absId,1.);
	  fhTimeSumSq[nBC]->Fill(absId,hkdtime*hkdtime);
	  fhTimeSum[nBC]->Fill(absId,hkdtime);
	}else{
	  fhTimeLGEnt[nBC]->Fill(absId,1.);
	  fhTimeLGSumSq[nBC]->Fill(absId,hkdtime*hkdtime);
	  fhTimeLGSum[nBC]->Fill(absId,hkdtime);
	}


      } // hkdtime:[-20;20]
    } // end icell
  } //end cluster
  
  
  // Post output data.
  //cout<<"Post data and delete caloClusters"<<endl;
  caloClusters->Delete();
  delete caloClusters;
// } // end if trigger type 

  PostData(1, fOutputList);  
} // End of AliAnalysisTaskEMCALTimeCalib::UserExec()

//________________________________________________________________________
/// Draw result to the screen
/// Called once at the end of the query
void AliAnalysisTaskEMCALTimeCalib::Terminate(Option_t *)
{
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  
  if(fTOFmaker) delete fTOFmaker;

  if(fL1PhaseList) {
    fL1PhaseList->SetOwner();
    fL1PhaseList->Clear();
    delete fL1PhaseList;
  }

  if (fBadChannelMapArray) { 
    fBadChannelMapArray->Clear();
    delete fBadChannelMapArray;
  }

  if (!fOutputList) 
  {
    AliDebug(1,"ERROR: Output list not available");
    return;
  }
} // End of AliAnalysisTaskEMCALTimeCalib::Terminate

//________________________________________________________________________
/// Selection criteria of good cluster are set here
Bool_t AliAnalysisTaskEMCALTimeCalib::AcceptCluster(AliVCluster* clus)
{
  //fix with noisy EMCAL fee card
  Int_t nCells = clus->GetNCells();
  
  if(clus->IsEMCAL())
  {    
    if ((clus->E() > fMaxClusterEnergy && nCells > fMaxNcells ) || nCells > fMaxNcells){
      AliDebug(1,"very big cluster with enormous energy - cluster rejected");
      return kFALSE;
    }
  }
  
  // remove other than photonlike
  Double_t lambda0=clus->GetM02();
  if (lambda0>fMaxLambda0LG || lambda0<fMinLambda0LG){
    AliDebug(1,"lambda0 loose cut failed - cluster rejected");
    return kFALSE;
  }

  // remove matched clusters
  Double_t Dx=clus->GetTrackDx();
  Double_t Dz=clus->GetTrackDz();
  Double_t Rtrack = TMath::Sqrt(Dx*Dx+Dz*Dz);
  if (Rtrack <fMaxRtrack) 
  {
    AliDebug(1,"track matched - cluster rejected");
    return kFALSE;
  }

  if (nCells<fMinNcells) 
  {
    AliDebug(1,"single cell cluster - cluster rejected");
    return kFALSE;
  }

  if(clus->E()<fMinClusterEnergy) 
  {
    AliDebug(1,"cluster energy < 1 GeV- cluster rejected");
    return kFALSE;
  }


  if(!IsLowGainCellInCluster(clus)) {//no low gain cell in cluster
    //apply more strict lambda0^2 cut
    if (lambda0>fMaxLambda0 || lambda0<fMinLambda0){
      AliDebug(1,"lambda0 strict cut failed - cluster rejected");
      return kFALSE;
    }
  }




  return kTRUE;
}//End AliAnalysisTaskEMCALTimeCalib::AcceptCluster

//________________________________________________________________________
/// Check if low gain cell is in a cluster
Bool_t  AliAnalysisTaskEMCALTimeCalib::IsLowGainCellInCluster(AliVCluster* clus){
  UShort_t * index = clus->GetCellsAbsId() ;
  AliVCaloCells &cells= *(InputEvent()->GetEMCALCells());
  for(Int_t i = 0; i < clus->GetNCells() ; i++) {
    if(cells.GetCellHighGain(index[i])==kFALSE) return kTRUE;//low gain cell
  }
  return kFALSE;

}

//________________________________________________________________________
/// Check RCU for cell given by Super Module, column index, row index
Bool_t AliAnalysisTaskEMCALTimeCalib::CheckCellRCU(Int_t nSupMod,Int_t icol,Int_t irow)
{
  Int_t iRCU;
  if(nSupMod < 10 || (nSupMod >= 12 && nSupMod <18) ) 
  {
    if      (0<=irow&&irow<8)                       iRCU=0; // first cable row
    else if (8<=irow&&irow<16 &&  0<=icol&&icol<24) iRCU=0; // first half; 
    //second cable row
    //RCU1
    else if (8<=irow&&irow<16 && 24<=icol&&icol<48) iRCU=1; // second half; 
    //second cable row
    else if (16<=irow&&irow<24)                     iRCU=1; // third cable row
    
    if (nSupMod%2==1) iRCU = 1 - iRCU; // swap for odd=C side, to allow us to cable both sides the same
  } 
  else 
  {
	// Last 2 SM have one single SRU, just assign RCU 0
    iRCU = 0 ;
  }

  //cout<<"RCU:"<<iRCU<<endl;
  if (iRCU<0) 
    AliFatal(Form("Wrong EMCAL/DCAL RCU number = %d\n", iRCU));

  return kTRUE;
}//End AliAnalysisTaskEMCALTimeCalib::CheckCellRCU

//________________________________________________________________________
/// Set default cuts for calibration
void AliAnalysisTaskEMCALTimeCalib::SetDefaultCuts()
{
  fMinClusterEnergy=1.0;//0.5//0.7
  fMaxClusterEnergy=500;
  fMinNcells=2;
  fMaxNcells=200;
  fMinLambda0=0.1;
  fMaxLambda0=0.4;
  fMinLambda0LG=0.1;
  fMaxLambda0LG=4.0;
  fMaxRtrack=0.025;
  fMinCellEnergy=0.4;//0.1//0.4
  fReferenceFileName="";//Reference.root
  fReferenceRunByRunFileName="";
  fBadReco=kFALSE;
  fFillHeavyHisto=kFALSE;
  fGeometryName="";//EMCAL_COMPLETE12SMV1_DCAL_8SM
  fPileupFromSPD=kFALSE;
  fMinTime=-20.;
  fMaxTime=20.;
  fMostEneCellOnly=kFALSE;
  
  fBadChannelMapSet=kFALSE;
  fSetBadChannelMapSource=0;
  fBadChannelFileName="";

  //histograms
  fRawTimeNbins  = 400;  // Raw time settings should be like that all the time
  fRawTimeMin    = 400.; // importent in pass1
  fRawTimeMax    = 800.;
  fPassTimeNbins = 1000; // in pass2 should be (400,400,800)
  fPassTimeMin   = -250.;// in pass3 should be (1000,-250,250)
  fPassTimeMax   = 250.;
  fEnergyNbins   = 100;  // default settings was 500
  fEnergyMin     = 0.;
  fEnergyMax     = 20.;
  fEnergyLGNbins = 200;  // default settings
  fEnergyLGMin   = 0.;
  fEnergyLGMax   = 100.;
  fFineNbins     = 90; //was 4500 for T0 time studies
  fFineTmin      = -500;
  fFineTmax      = 400;
}

//________________________________________________________________________
/// Calculate calibration constants
/// input - root file with histograms 
/// output - root file with constants in historams
/// isFinal - flag: kFALSE-first iteration, kTRUE-final iteration
void AliAnalysisTaskEMCALTimeCalib::ProduceCalibConsts(TString inputFile,TString outputFile,Bool_t isFinal)
{
  TFile *file =new TFile(inputFile.Data());
  if(file==0x0) {
    //AliWarning("Input file does not exist!");
    return;
  }

  TList *list=(TList*)file->Get("chistolist");
  if(list==0x0) 
  {
    //AliWarning("List chistolist does not exist in file!");
    return;
  }

  //high gain
  TH1F *h1[4];
  TH1F *h2[4];
  TH1F *h3[4];
  TH1F *hAllTimeAvBC[4];
  TH1F *hAllTimeRMSBC[4];

  //low gain
  TH1F *h4[4];
  TH1F *h5[4];
  TH1F *h6[4];
  TH1F *hAllTimeAvLGBC[4];
  TH1F *hAllTimeRMSLGBC[4];

  if(isFinal==kFALSE){//first itereation
    for(Int_t i=0;i<4;i++){
      h1[i]=(TH1F *)list->FindObject(Form("RawTimeSumBC%d",i));
      h2[i]=(TH1F *)list->FindObject(Form("RawTimeEntriesBC%d",i));
      h3[i]=(TH1F *)list->FindObject(Form("RawTimeSumSqBC%d",i));
      
      h4[i]=(TH1F *)list->FindObject(Form("RawTimeSumLGBC%d",i));
      h5[i]=(TH1F *)list->FindObject(Form("RawTimeEntriesLGBC%d",i));
      h6[i]=(TH1F *)list->FindObject(Form("RawTimeSumSqLGBC%d",i));
    }
  } else {//final iteration
    for(Int_t i=0;i<4;i++){
      h1[i]=(TH1F *)list->FindObject(Form("hTimeSum%d",i));
      h2[i]=(TH1F *)list->FindObject(Form("hTimeEnt%d",i));
      h3[i]=(TH1F *)list->FindObject(Form("hTimeSumSq%d",i));
      
      h4[i]=(TH1F *)list->FindObject(Form("hTimeLGSum%d",i));
      h5[i]=(TH1F *)list->FindObject(Form("hTimeLGEnt%d",i));
      h6[i]=(TH1F *)list->FindObject(Form("hTimeLGSumSq%d",i));
    }
  }
  //AliWarning("Input histograms read.");

  for(Int_t i=0;i<4;i++){
    hAllTimeAvBC[i]=new TH1F(Form("hAllTimeAvBC%d",i),Form("hAllTimeAvBC%d",i),h1[i]->GetNbinsX(),h1[i]->GetXaxis()->GetXmin(),h1[i]->GetXaxis()->GetXmax());
    hAllTimeRMSBC[i]=new TH1F(Form("hAllTimeRMSBC%d",i),Form("hAllTimeRMSBC%d",i),h3[i]->GetNbinsX(),h3[i]->GetXaxis()->GetXmin(),h3[i]->GetXaxis()->GetXmax());

    hAllTimeAvLGBC[i]=new TH1F(Form("hAllTimeAvLGBC%d",i),Form("hAllTimeAvLGBC%d",i),h4[i]->GetNbinsX(),h4[i]->GetXaxis()->GetXmin(),h4[i]->GetXaxis()->GetXmax());
    hAllTimeRMSLGBC[i]=new TH1F(Form("hAllTimeRMSLGBC%d",i),Form("hAllTimeRMSLGBC%d",i),h6[i]->GetNbinsX(),h6[i]->GetXaxis()->GetXmin(),h6[i]->GetXaxis()->GetXmax());
  }
  
  //AliWarning("New histograms booked.");

  //important remark: we use 'underflow bin' for absid=0 in OADB  . That's why there is j-1 below.
  for(Int_t i=0;i<4;i++){
    for(Int_t j=1;j<=h1[i]->GetNbinsX();j++){
      //high gain
      if(h2[i]->GetBinContent(j)!=0){
	hAllTimeAvBC[i]->SetBinContent(j-1,h1[i]->GetBinContent(j)/h2[i]->GetBinContent(j));
	hAllTimeRMSBC[i]->SetBinContent(j-1,TMath::Sqrt(h3[i]->GetBinContent(j)/h2[i]->GetBinContent(j)) );
      } else {
	hAllTimeAvBC[i]->SetBinContent(j-1,0.);
	hAllTimeRMSBC[i]->SetBinContent(j-1,0.);
      }
      //low gain
      if(h5[i]->GetBinContent(j)!=0){
	hAllTimeAvLGBC[i]->SetBinContent(j-1,h4[i]->GetBinContent(j)/h5[i]->GetBinContent(j));
	hAllTimeRMSLGBC[i]->SetBinContent(j-1,TMath::Sqrt(h6[i]->GetBinContent(j)/h5[i]->GetBinContent(j)) );
      } else {
	hAllTimeAvLGBC[i]->SetBinContent(j-1,0.);
	hAllTimeRMSLGBC[i]->SetBinContent(j-1,0.);
      }

    }
  }

  //AliWarning("Average and rms calculated.");
  TFile *fileNew=new TFile(outputFile.Data(),"recreate");
  for(Int_t i=0;i<4;i++){
    hAllTimeAvBC[i]->Write();
    hAllTimeRMSBC[i]->Write();
    hAllTimeAvLGBC[i]->Write();
    hAllTimeRMSLGBC[i]->Write();
  }

  //AliWarning(Form("Histograms saved in %s file.",outputFile.Data()));

  fileNew->Close();
  delete fileNew;

  for(Int_t i=0;i<4;i++){
    delete hAllTimeAvBC[i];
    delete hAllTimeRMSBC[i];
    delete hAllTimeAvLGBC[i];
    delete hAllTimeRMSLGBC[i];

    delete h1[i];
    delete h2[i];
    delete h3[i];
    delete h4[i];
    delete h5[i];
    delete h6[i];
  }
  
  file->Close();
  delete file;

  //AliWarning("Pointers deleted. Memory cleaned.");
}

//________________________________________________________________________
/// Calculate calibration constants per SM (equivalent of L1 phase)
/// input - root file with calibration constants from 1st pass
/// output - root file with histograms for given run offset per SM 
void AliAnalysisTaskEMCALTimeCalib::ProduceOffsetForSMsV2(Int_t runNumber,TString inputFile,TString outputFile, Bool_t offset100, Bool_t justL1phase){

const  Double_t lowerLimit[]={
    0,
    1152,
    2304,
    3456,
    4608,
    5760,
    6912,
    8064,
    9216,
    10368,
    11520,
    11904,
    12288,
    13056,
    13824,
    14592,
    15360,
    16128,
    16896,
    17280};

const  Double_t upperLimit[]={
    1151 ,
    2303 ,
    3455 ,
    4607 ,
    5759 ,
    6911 ,
    8063 ,
    9215 ,
    10367,
    11519,
    11903,
    12287,
    13055,
    13823,
    14591,
    15359,
    16127,
    16895,
    17279,
    17663};

  TFile *file =new TFile(inputFile.Data());
  if(file==0x0) return;

  TH1F *ccBC[4];
  Bool_t shouldBeEmpty[4];
  Int_t emptyCounter;
  for(Int_t i = 0; i < kNBCmask; i++){
    ccBC[i]=(TH1F*) file->Get(Form("hAllTimeAvBC%d",i));
    shouldBeEmpty[i]=kFALSE;
    emptyCounter=0;
    for(Int_t j=0;j<upperLimit[19];j++){
      if(ccBC[i]->GetBinContent(j)>0.) emptyCounter++;
    }
    if(emptyCounter<1500) shouldBeEmpty[i]=kTRUE;
    cout<<"Non-zero channels "<<emptyCounter<<" BC"<<i<<" should be empty: "<<shouldBeEmpty[i]<<endl;
  }

  TH1C *hRun=new TH1C(Form("h%d",runNumber),Form("h%d",runNumber),19,0,19);
  Int_t fitResult=0;
  Double_t minimumValue=10000.;
  Int_t minimumIndex=-1;
  Double_t meanBC[4];

  Double_t fitParameter=0;
  TF1 *f1=new TF1("f1","pol0",0,17664);
  Bool_t orderTest=kTRUE;
  Int_t iorder=0;//order index
  Int_t j=0;//BC index
  Int_t L1shift=0;
  Int_t totalValue=0;

  for(Int_t i=0;i<20;i++){
    minimumValue=10000;
    for(j=0;j<kNBCmask;j++){
      if(shouldBeEmpty[j]) {
	meanBC[j]=-1;
	continue;
      }
      fitResult=ccBC[j]->Fit("f1","CQ","",lowerLimit[i],upperLimit[i]);
      if(fitResult<0){
	//hRun->SetBinContent(i,0);//correct it please
	meanBC[j]=-1;
	printf("Fit failed for SM %d BC%d, integral %f\n",i,j,ccBC[j]->Integral(lowerLimit[i],upperLimit[i]));
	continue;
      } else {
	fitParameter = f1->GetParameter(0);
      }
      if(offset100 && (j==0 || j==1)) {
	//the 100 ns offset was removed in LHC15n muon_calo_pass1 and further reconstructions 
	fitParameter+=100;
      }
      meanBC[j]=fitParameter;
	
      if(fitParameter>0 && fitParameter<minimumValue){
	minimumValue = fitParameter;
	minimumIndex = j;
      }
    }

    if( minimumValue/25-(Int_t)(minimumValue/25)>0.5 ) {
      L1shift=(Int_t)(minimumValue/25.)+1;
    } else {
      L1shift=(Int_t)(minimumValue/25.);
    }

    if(TMath::Abs(minimumValue/25-(Int_t)(minimumValue/25)-0.5)<0.05)
      printf("Run %d, SM %d, min %f, next_min %f, next+1_min %f, next+2_min %f, min/25 %f, min%%25 %d, next_min/25 %f, next+1_min/25 %f, next+2_min/25 %f, SMmin %d\n",runNumber,i,minimumValue,meanBC[(minimumIndex+1)%4],meanBC[(minimumIndex+2)%4],meanBC[(minimumIndex+3)%4],minimumValue/25., (Int_t)((Int_t)minimumValue%25), meanBC[(minimumIndex+1)%4]/25., meanBC[(minimumIndex+2)%4]/25., meanBC[(minimumIndex+3)%4]/25., L1shift*25);

    if(justL1phase) totalValue = minimumIndex;
    else totalValue = L1shift<<2 | minimumIndex ;
    //printf("L1 phase %d, L1 shift %d *25ns= %d, L1p+L1s %d, total %d, L1pback %d, L1sback %d\n",minimumIndex,L1shift,L1shift*25,minimumIndex+L1shift,totalValue,totalValue&3,totalValue>>2);

    hRun->SetBinContent(i,totalValue);
    orderTest=kTRUE;
    for(iorder=minimumIndex;iorder<minimumIndex+4-1;iorder++){
      if( meanBC[(iorder+1)%4] <= meanBC[iorder%4] ) orderTest=kFALSE;
    }
    if(!orderTest)
      printf("run %d, SM %d, min index %d meanBC %f %f %f %f, order ok? %d\n",runNumber,i,minimumIndex,meanBC[0],meanBC[1],meanBC[2],meanBC[3],orderTest);

    //patch for runs with not filled one, two or three BCs
    //manual patch for LHC16q - pPb@5TeV - only BC0 is filled and phase rotate
    if(shouldBeEmpty[0] || shouldBeEmpty[1] || shouldBeEmpty[2] || shouldBeEmpty[3]){
    Double_t newMean = meanBC[minimumIndex]-600;
    if(newMean<=12.5) hRun->SetBinContent(i,minimumIndex);
    else {
      Int_t minIndexTmp=-1;
      if(newMean/25. - (Int_t)(newMean/25.) <0.5)
	minIndexTmp = (Int_t)(newMean/25.);
      else
	minIndexTmp = 1+(Int_t)(newMean/25.);

      hRun->SetBinContent(i,(4-minIndexTmp+minimumIndex)%4);
      //cout<<newMean/25.<<" int "<<(Int_t)(newMean/25.)<<" dif "<< newMean/25.-(Int_t)(newMean/25.)<<endl;
      }
    cout << "run with missing BC; new L1 phase set to " << hRun->GetBinContent(i)<<endl;
    }//end of patch for LHC16q and other runs with not filled BCs
  }//end of loop over SM

  delete f1;
  TFile *fileNew=new TFile(outputFile.Data(),"update");
  hRun->Write();
  fileNew->Close();
  delete fileNew;

  file->Close();
  delete file;
}

//____________________________________________________
void AliAnalysisTaskEMCALTimeCalib::LoadBadChannelMapOADB()
{
  if(fBadChannelMapSet) return;
  AliOADBContainer *contBC=new AliOADBContainer("");
  contBC->InitFromFile(Form("%s/EMCALBadChannels.root","$ALICE_PHYSICS/OADB/EMCAL"),"AliEMCALBadChannels"); 
  printf("contBC %p, ent  %d\n",contBC,contBC->GetNumberOfEntries());
  TObjArray *arrayBC=(TObjArray*)contBC->GetObject(fRunNumber);
  if(arrayBC) {
    AliInfo("Remove EMCAL bad cells");
    fBadChannelMapArray = new TObjArray(kNSM);
    for (Int_t i=0; i<kNSM; ++i) {
      TH2I *hbm = (TH2I*)arrayBC->FindObject(Form("EMCALBadChannelMap_Mod%d",i));
      if (!hbm) {
	AliError(Form("Can not get EMCALBadChannelMap_Mod%d",i));
	continue;
      }
      hbm->SetDirectory(0);
      fBadChannelMapArray->AddAt(hbm,i);
          
    } // loop over SMs
  } else AliInfo("Do NOT remove EMCAL bad channels\n"); // run array
      
  delete contBC;
  fBadChannelMapSet=kTRUE;
}  // Bad channel map loaded

//____________________________________________________
void AliAnalysisTaskEMCALTimeCalib::LoadBadChannelMapFile()
{
  if(fBadChannelMapSet) return;

  TFile *referenceFile = TFile::Open(fBadChannelFileName.Data());
  if(referenceFile==0x0) {
    AliFatal("*** NO bad channel map FILE");
  }

  TH1F *hbm = (TH1F*)referenceFile->Get("h1");
  if (!hbm) {
    AliError("Can not get EMCALBadChannelMap");
  }
  fBadChannelMapArray = new TObjArray(1);
  fBadChannelMapArray->AddAt(hbm,0);
  fBadChannelMapSet=kTRUE;
}  // Bad channel map loaded


//_____________________________________________________________________
/// Load Bad Channel Map from different source
void AliAnalysisTaskEMCALTimeCalib::LoadBadChannelMap(){
  if(fSetBadChannelMapSource==1) LoadBadChannelMapOADB();
  else if(fSetBadChannelMapSource==2) LoadBadChannelMapFile();
}
