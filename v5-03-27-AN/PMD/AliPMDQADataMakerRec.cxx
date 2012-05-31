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


/*
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
  B.K. Nandi
*/

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1D.h>  
#include <TH1F.h> 
#include <TH1I.h> 
#include <TH2F.h> 

// --- Standard library ---

// --- AliRoot header files ---

#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliPMDQADataMakerRec.h"
#include "AliQAChecker.h"
#include "AliPMDdigit.h" 
#include "AliPMDrecpoint1.h" 
#include "AliPMDRawStream.h"
#include "AliPMDddldata.h"
#include "AliPMDUtility.h"
#include "AliESDPmdTrack.h"
//#include "AliPMDRecoParam.h"

ClassImp(AliPMDQADataMakerRec)
           
//____________________________________________________________________________ 
  AliPMDQADataMakerRec::AliPMDQADataMakerRec() : 
  AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kPMD), "PMD Quality Assurance Data Maker")
{
  // ctor
}

//____________________________________________________________________________ 
AliPMDQADataMakerRec::AliPMDQADataMakerRec(const AliPMDQADataMakerRec& qadm) :
  AliQADataMakerRec()
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliPMDQADataMakerRec& AliPMDQADataMakerRec::operator = (const AliPMDQADataMakerRec& qadm )
{
  // Equal operator.
  this->~AliPMDQADataMakerRec();
  new(this) AliPMDQADataMakerRec(qadm);
  return *this;
}
 

//____________________________________________________________________________ 
void AliPMDQADataMakerRec::InitRaws() {
  // create Raws histograms in Raws subdir

  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  

  TH1F *h0 = new TH1F("hDdl2304","Cell Wise Frequency for PRE Plane DDL 2304;DDL2304;Frequency",50000 ,0,50000.0);
  TH1F *h1 = new TH1F("hDdl2305","Cell Wise Frequency for PRE Plane DDL 2305;DDL2305;Frequency",50000 ,0,50000.0);
  TH1F *h2 = new TH1F("hDdl2308","Cell Wise Frequency for CPV Plane DDL 2308;DDL2308;Frequency",50000 ,0,50000.0);
  TH1F *h3 = new TH1F("hDdl2309","Cell Wise Frequency for CPV Plane DDL 2309;DDL2309;Frequency",50000 ,0,50000.0);
    
  for(int i = 0; i < 50; i++ )  {
    h0->GetXaxis()->SetBinLabel(i*1000+300,Form("%2d",i+1));
    h1->GetXaxis()->SetBinLabel(i*1000+300,Form("%2d",i+1));
    h2->GetXaxis()->SetBinLabel(i*1000+300,Form("%2d",i+1));
    h3->GetXaxis()->SetBinLabel(i*1000+300,Form("%2d",i+1));
  }


  Add2RawsList(h0, 0, expert, image, saveCorr);  
  Add2RawsList(h1, 1, expert, image, saveCorr);  
  Add2RawsList(h2, 2, expert, image, saveCorr);  
  Add2RawsList(h3, 3, expert, image, saveCorr);  


  TH2F *hPreXY = new TH2F("hPreXY","PRE plane;X [cm];Y [cm]",200,-100.,100.,200,-100.,100.);
  Add2RawsList(hPreXY, 4, expert, !image, saveCorr);
  
  TH2F *hCpvXY = new TH2F("hCpvXY","CPV plane;X [cm];Y [cm]",200,-100.,100.,200,-100.,100.);
  Add2RawsList(hCpvXY, 5, expert, !image, saveCorr);
  
  TH1F *hQuality1 = new TH1F("hPmdQualityWAdc","Quality Plot SumWAdc/SumWAdc;#eta Bins (PRE   <--->    CPV);Quality",25,0,25);
  TH1F *hQuality2 = new TH1F("hPmdQualityHit","Quality Plot Hit/Hit;#eta Bins (PRE  <--->   CPV);Quality",25,0,25);

  for(int i = 0; i < 10; i++ )  {
    hQuality1->GetXaxis()->SetBinLabel(i+3,Form("%2d",i+1));
    hQuality2->GetXaxis()->SetBinLabel(i+3,Form("%2d",i+1));

    hQuality1->GetXaxis()->SetBinLabel(i+15,Form("%2d",i+1));
    hQuality2->GetXaxis()->SetBinLabel(i+15,Form("%2d",i+1));
  }

  Add2RawsList(hQuality1, 6, !expert, image, saveCorr);
  Add2RawsList(hQuality2, 7, !expert, image, saveCorr);

  TH1F *hMP = new TH1F("hPmdMultiplicityP"," Cell Hit Distribution with adc>40; Multiplicity; Entries",1000,0,7000); 
  TH1F *hMC = new TH1F("hPmdMultiplicityC"," Cell Hit Distribution with adc>40; Multiplicity; Entries",1000,0,7000); 
  Add2RawsList(hMP, 8, expert, !image, saveCorr);
  Add2RawsList(hMC, 9, expert, !image, saveCorr);
  



  ClonePerTrigClass(AliQAv1::kRAWS); // this should be the last line
}

//____________________________________________________________________________
void AliPMDQADataMakerRec::InitDigits()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F *h0 = new TH1F("hPreDigitsEdep","Digits energy distribution in PRE(PMD);Amplitude [ADC counts];Counts", 100, 0., 2000.);
  h0->Sumw2();
  Add2DigitsList(h0, 0, !expert, image);
  
  TH1F *h1 = new TH1F("hCpvDigitsEdep","Digits energy distribution in CPV(PMD);Amplitude [ADC counts];Counts", 100, 0., 2000.); 
  h1->Sumw2();
  Add2DigitsList(h1, 1, !expert, image);
  
  TH1I *h2 = new TH1I("hPreDigitsMult","Digits multiplicity distribution in PRE(PMD);# of Digits;Entries", 500, 0, 1000) ; 
  h2->Sumw2();
  Add2DigitsList(h2, 2, !expert, image);
  
  TH1I *h3 = new TH1I("hCpvDigitsMult","Digits multiplicity distribution in CPV(PMD);# of Digits;Entries", 500, 0, 1000);
  h3->Sumw2();
  Add2DigitsList(h3, 3, !expert, image);  
  //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
}

//____________________________________________________________________________ 
void AliPMDQADataMakerRec::InitRecPoints()
{
  // create Reconstructed Points histograms in RecPoints subdir
  
  /*
    TH2F * h0 = new TH2F("hPreXY","RecPoints Y vs X PRE(PMD)", 100,-100.,100.,100,-100.,100.);
    Add2RecPointsList(h0,0);
    
    TH2F * h1 = new TH2F("hCpvXY","RecPoints Y vs X CPV(PMD)", 100,-100.,100.,100,-100.,100.);
    Add2RecPointsList(h1,1);
  */

    //  Ncell distribution in a cluster
  
    // PRE plane

  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH1F * h0 = new TH1F("hPreDdl0Ncell","PRE: Ddl0 Ncell in a cluster;# of cells;Counts",50,0.,50.);
  h0->Sumw2();
  Add2RecPointsList(h0, 0, !expert, image);


  TH1F * h1 = new TH1F("hPreDdl1Ncell","PRE: Ddl1 Ncell in a cluste;# of cells;Countsr",50,0.,50.);
  h1->Sumw2();
  Add2RecPointsList(h1, 1, !expert, image);


  TH1F * h2 = new TH1F("hPreDdl2Ncell","PRE: Ddl2 Ncell in a cluster;# of cells;Counts",50,0.,50.);
  h2->Sumw2();
  Add2RecPointsList(h2, 2, !expert, image);


  TH1F * h3 = new TH1F("hPreDdl3Ncell","PRE: Ddl3 Ncell in a cluster;# of cells;Counts",50,0.,50.);
  h3->Sumw2();
  Add2RecPointsList(h3, 3, !expert, image);

  // CPV plane

  TH1F * h4 = new TH1F("hCpvDdl4Ncell","CPV: Ddl4 Ncell in a cluster;# of cells;Counts",50,0.,50.);
  h4->Sumw2();
  Add2RecPointsList(h4, 4, !expert, image);

  TH1F * h5 = new TH1F("hCpvDdl5Ncell","CPV: Ddl5 Ncell in a cluster;# of cells;Counts",50,0.,50.);
  h5->Sumw2();
  Add2RecPointsList(h5, 5, !expert, image);

  // Correlation plot

  TH2I *h6 = new TH2I("hPre10","Cluster - DDL1 vs DDL0;DDL0;DDL1", 100,0,200,100,0,200);
  Add2RecPointsList(h6,6, !expert, image);

  TH2I *h7 = new TH2I("hPre32","Cluster - DDL3 vs DDL2;DDL2;DDL3", 100,0,200,100,0,200);
  Add2RecPointsList(h7,7, !expert, image);


  TH2I *h8 = new TH2I("hCpv54","Cluster - DDL5 vs DDL4;DDL4;DDL5", 100,0,200,100,0,200);
  Add2RecPointsList(h8,8, !expert, image);
  //
  ClonePerTrigClass(AliQAv1::kRECPOINTS); // this should be the last line
}

//____________________________________________________________________________ 

void AliPMDQADataMakerRec::InitESDs()
{
  //Create histograms to controll ESD

  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH1F *h0 = new TH1F("hPreClADC","Cluster ADC of PRE plane;# of clusters;Counts",500,0.,10000.);
  h0->Sumw2();
  Add2ESDsList(h0, 0, !expert, image)  ;

  TH1F *h1 = new TH1F("hCpvClADC","Cluster ADC of CPV plane;# of clusters;Counts",500,0.,10000.);
  h1->Sumw2();
  Add2ESDsList(h1, 1, !expert, image)  ;

  TH2I *h2 = new TH2I("hPmdClMult","Cluster Multiplicity: PRE vs. CPVplane;CPV multiplicity;PRE Multiplicity",100,0,1000,100,0,1000);
  h2->Sumw2();
  Add2ESDsList(h2, 2, !expert, image)  ;

  /*
    TH1I * h3 = new TH1I("hCpvClMult","Cluster Multiplicity of CPV plane",100,0.,1000.);
    h3->Sumw2();
    Add2ESDsList(h3, 3, !expert, image)  ;
  */
  //
  ClonePerTrigClass(AliQAv1::kESDS); // this should be the last line
}

//____________________________________________________________________________
void AliPMDQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  //Fill prepared histograms with Raw digit properties

  TObjArray *pmdddlcont = 0x0;
    pmdddlcont = new TObjArray();
    AliPMDRawStream stream(rawReader);
    
    AliPMDddldata *pmdddl = 0x0;

    Int_t   iddl =  1;
    Int_t   sddl = -1;
    Float_t xx = 0., yy = 0.;

    Double_t eta;
    Int_t etaindex = 0;

    Float_t nDdl1[10] = {0,0,0,0,0,0,0,0,0,0};
    Float_t nDdl2[10] = {1,1,1,1,1,1,1,1,1,1};
    Float_t nDdl3[10] = {0,0,0,0,0,0,0,0,0,0};
    Float_t nDdl4[10] = {1,1,1,1,1,1,1,1,1,1};

    Float_t nDdl1a[10] = {1,1,1,1,1,1,1,1,1,1};
    Float_t nDdl2a[10] = {1,1,1,1,1,1,1,1,1,1};
    Float_t nDdl3a[10] = {1,1,1,1,1,1,1,1,1,1};
    Float_t nDdl4a[10] = {1,1,1,1,1,1,1,1,1,1};
    
    Float_t nDdl11[10] = {1,1,1,1,1,1,1,1,1,1};
    Float_t nDdl22[10] = {1,1,1,1,1,1,1,1,1,1};
    Float_t nDdl33[10] = {1,1,1,1,1,1,1,1,1,1};
    Float_t nDdl44[10] = {1,1,1,1,1,1,1,1,1,1};

    Float_t nddlp = 0;
    Float_t nddlc  = 0;

    AliPMDUtility cc;

       

    while ((iddl = stream.DdlData(pmdddlcont)) >=0)
    {
	Int_t ientries = pmdddlcont->GetEntries();
	//printf(" ======= DDLNO = %d ientries = %d \n", iddl, ientries);
	if(iddl == 0 || iddl == 1) sddl = iddl;
	if(iddl == 4 ) sddl = 2;
	if(iddl == 5 ) sddl = 3;
	
	for (Int_t ient = 0; ient < ientries; ient++)
	  {
	    pmdddl = (AliPMDddldata*)pmdddlcont->UncheckedAt(ient);

	    Int_t det  = pmdddl->GetDetector();
	    Int_t smn  = pmdddl->GetSMN();
	    Int_t mcm  = pmdddl->GetMCM();
	    Int_t chno = pmdddl->GetChannel();
	    Int_t pbus = pmdddl->GetPatchBusId();
	    Int_t row  = pmdddl->GetRow();
	    Int_t col  = pmdddl->GetColumn();
	    Int_t sig  = pmdddl->GetSignal();

	    if (mcm == 0) continue;

	    if (det < 0 || det > 1)  continue;
	    if (smn < 0 || smn > 23) continue;
	    if (row < 0 || row > 47) continue;
	    if (col < 0 || col > 95) continue;
	    
	    Int_t rc = pbus*1000 + mcm*64 + chno;

	    if(sddl == 0) FillRawsData(0,rc);
	    else if(sddl == 1) FillRawsData(1,rc);
	    else if(sddl == 2) FillRawsData(2,rc);
	    else if(sddl == 3) FillRawsData(3,rc);

	    cc.GetEtaIndexXY(smn,row,col,xx,yy,eta,etaindex);
	    
	    if(etaindex >= 0 && etaindex <10) {

	      if(sddl == 0) { 
		nDdl1a[etaindex] += 1; 
		nDdl1[etaindex] += sig; 
		if(sig >40) {
		  nDdl11[etaindex] += 1;
		  nddlp += 1;
		}
		  
	      }
	      else if(sddl == 1) { 
		nDdl2a[etaindex] += 1; 
		nDdl2[etaindex] += sig; 
		if(sig >40) {
		  nDdl22[etaindex] += 1;
		  nddlp += 1;
		}
 

	      }
	      else if(sddl == 2) { 
		nDdl3a[etaindex] += 1; 
		nDdl3[etaindex] += sig; 
		if(sig >40) {
		  nDdl33[etaindex] += 1;
		  nddlc += 1;
		} 
	      }
	      else if(sddl == 3) { 
		nDdl4a[etaindex] += 1; 
		nDdl4[etaindex] += sig; 
		if(sig >40) {
		  nDdl44[etaindex] += 1; 
		  nddlc += 1;
		}
	      }

	    }

	    if (det == 0) FillRawsData(4,xx,yy);
	    else if (det == 1) FillRawsData(5,xx,yy);
	    
	  }

	pmdddlcont->Delete();
    }


    delete pmdddlcont;
    pmdddlcont = 0x0;
    
    if(nddlp != 0)   FillRawsData(8,nddlp);
    if(nddlc != 0)   FillRawsData(9,nddlc);
   

    ResetRawsData(6);
    ResetRawsData(7);
    for (Int_t i = 0; i < 10; i++) {
     
      Float_t prerC =  nDdl11[i]/nDdl22[i];
      Float_t cpvrC = nDdl33[i]/nDdl44[i];

      Float_t prera =(nDdl1[i]/nDdl1a[i])/(nDdl2[i]/nDdl2a[i]);
      Float_t cpvra =(nDdl3[i]/nDdl3a[i])/(nDdl4[i]/nDdl4a[i]);
      
     
      FillRawsData(6,i+2,prera);
      FillRawsData(6,i+14,cpvra);
      FillRawsData(7,i+2,prerC);
      FillRawsData(7,i+14,cpvrC);
    }

    IncEvCountCycleRaws();
    IncEvCountTotalRaws();
 
    //
}
//____________________________________________________________________________
void AliPMDQADataMakerRec::MakeDigits()
{
  // makes data from Digits
  
   Int_t cpvmul = 0, premul = 0;
  
  TIter next(fDigitsArray) ; 
  AliPMDdigit * digit ; 
  while ( (digit = dynamic_cast<AliPMDdigit *>(next())) )
    {
    if(digit->GetDetector() == 0)
      {
	    FillDigitsData(0, digit->GetADC()) ;
	    premul++;
      }
    if(digit->GetDetector() == 1)
      {
	    FillDigitsData(1, digit->GetADC());
	    cpvmul++;
      }
    }  
  
  if (premul > 0) FillDigitsData(2,premul);
  if (cpvmul > 0) FillDigitsData(3,cpvmul);
  
  
}

//____________________________________________________________________________
void AliPMDQADataMakerRec::MakeDigits(TTree * digitTree)
{
  // makes data from Digit Tree
  
  if (fDigitsArray) 
    fDigitsArray->Clear() ; 
  else
    fDigitsArray = new TClonesArray("AliPMDdigit", 1000) ; 
  
  TBranch * branch = digitTree->GetBranch("PMDDigit") ;
  if ( ! branch ) {AliWarning("PMD branch in Digit Tree not found"); return;} 
  //
  branch->SetAddress(&fDigitsArray) ;
  for (Int_t ient = 0; ient < branch->GetEntries(); ient++) {
    branch->GetEntry(ient) ; 
    MakeDigits() ; 
  }
  //
  IncEvCountCycleDigits();
  IncEvCountTotalDigits();
  //
}

//____________________________________________________________________________
void AliPMDQADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
    // makes data from RecPoints

  Int_t multDdl0 = 0, multDdl1 = 0, multDdl2 = 0;
  Int_t multDdl3 = 0, multDdl4 = 0, multDdl5 = 0;
  
  AliPMDrecpoint1 * recpoint; 
  
  if (fRecPointsArray) 
    fRecPointsArray->Clear() ; 
  else 
    fRecPointsArray = new TClonesArray("AliPMDrecpoint1", 1000) ; 
  
  TBranch * branch = clustersTree->GetBranch("PMDRecpoint") ;

  if ( ! branch )
    {
      AliWarning("PMD branch in Recpoints Tree not found") ; 
    }
  else
    {
      branch->SetAddress(&fRecPointsArray) ;
      
      for (Int_t imod = 0; imod < branch->GetEntries(); imod++)
	{
	  branch->GetEntry(imod) ;
	  
	  TIter next(fRecPointsArray) ; 
	  
	  while ( (recpoint = dynamic_cast<AliPMDrecpoint1 *>(next())) )
	    {
	      //Float_t xpos = recpoint->GetClusX();
	      //Float_t ypos = recpoint->GetClusY();
	      //Int_t smn = recpoint->GetSMNumber();
	      
	      if(recpoint->GetDetector() == 0)
		{
		  if(recpoint->GetSMNumber() >= 0 && recpoint->GetSMNumber() < 6)
		    {
		      FillRecPointsData(0,recpoint->GetClusCells());
		      multDdl0++;
		    }
		  if(recpoint->GetSMNumber() >= 6 && recpoint->GetSMNumber() < 12)
		    {
		      FillRecPointsData(1,recpoint->GetClusCells());
		      multDdl1++;
		    }
		  if(recpoint->GetSMNumber() >= 12 && recpoint->GetSMNumber() < 18)
		    {
		      FillRecPointsData(2,recpoint->GetClusCells());
		      multDdl2++;
		    }
		  if(recpoint->GetSMNumber() >= 18 && recpoint->GetSMNumber() < 24)
		    {
		      FillRecPointsData(3,recpoint->GetClusCells());
		      multDdl3++;
		    }
		}
	      
	      if(recpoint->GetDetector() == 1)
		{
		  if((recpoint->GetSMNumber() >= 0 && recpoint->GetSMNumber() < 6) || 
		     (recpoint->GetSMNumber() >= 18 && recpoint->GetSMNumber() < 24))
		    {
		      FillRecPointsData(4,recpoint->GetClusCells());
		      multDdl4++;
		    }
		  if(recpoint->GetSMNumber() >= 6 && recpoint->GetSMNumber() < 18 )
		    {
		      FillRecPointsData(5,recpoint->GetClusCells());
		      multDdl5++;
		    }
		}
	    } 
	}
    }
  
  FillRecPointsData(6,multDdl0,multDdl1);
  FillRecPointsData(7,multDdl2,multDdl3);
  FillRecPointsData(8,multDdl4,multDdl5);
  //
  IncEvCountCycleRecPoints();
  IncEvCountTotalRecPoints();  
  //
}

//____________________________________________________________________________

void AliPMDQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs

  Int_t premul = 0, cpvmul = 0;

  for (Int_t icl = 0; icl < esd->GetNumberOfPmdTracks(); icl++)
    {
      AliESDPmdTrack *pmdtr = esd->GetPmdTrack(icl);
      
      //Int_t   det   = pmdtr->GetDetector(); 
      //Float_t clsX  = pmdtr->GetClusterX();
      //Float_t clsY  = pmdtr->GetClusterY();
      //Float_t clsZ  = pmdtr->GetClusterZ();
      //Float_t ncell = pmdtr->GetClusterCells();
      Float_t adc   = pmdtr->GetClusterADC();
      //Float_t pid   = pmdtr->GetClusterPID();

      if (pmdtr->GetDetector() == 0)
	{
	  FillESDsData(0,adc);
	  premul++;
	}
      if (pmdtr->GetDetector() == 1)
	{
	  FillESDsData(1,adc) ;
	  cpvmul++;
	}
    }
  
  FillESDsData(2,cpvmul,premul) ;
  //FillESDsData(3,cpvmul) ;  
  //
  IncEvCountCycleESDs();
  IncEvCountTotalESDs(); 
  //
}

//____________________________________________________________________________ 

void AliPMDQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}
//____________________________________________________________________________ 
void AliPMDQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kPMD, task, list) ;  
}
