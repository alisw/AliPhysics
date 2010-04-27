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
void AliPMDQADataMakerRec::InitRaws()
{
  // create Raws histograms in Raws subdir

  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F *hCellAdcPRE = new TH1F("hCellAdcPRE", "Cell Wise ADC fill for PRE Shower Plane;Cells in SMN-ROW-COL;ADC", 234795, 0, 234795);
 
  Add2RawsList(hCellAdcPRE, 0, !expert, image, !saveCorr);  


  TH1F *hCellAdcCPV = new TH1F("hCellAdcCPV", "Cell Wise ADC fill for CPV Veto Plane;Cells in SMN-ROW-COL;ADC", 234795, 0, 234795);
  
  Add2RawsList(hCellAdcCPV, 1, !expert, image, !saveCorr);  

  
  TH1F *hCalibEntPRE = new TH1F("hCalibEntPRE", "Cell Wise Frequency for PRE Shower Plane;Cells in SMN-ROW-COL;Frequency", 234795, 0, 234795);
  
  Add2RawsList(hCalibEntPRE, 2, !expert, image, !saveCorr);  


  TH1F *hCellEntCPV = new TH1F("hCellEntCPV", "Cell Wise Frequency for CPV Veto Plane;Cells in SMN-ROW-COL;Frequency", 234795, 0, 234795);
  
  Add2RawsList(hCellEntCPV, 3, !expert, image, !saveCorr);  

  
  TH2F *hPreXY = new TH2F("hPreXY","PRE plane;X [cm];Y [cm]",200,-100.,100.,200,-100.,100.);
  Add2RawsList(hPreXY, 4, !expert, !image, saveCorr);
  
  TH2F *hCpvXY = new TH2F("hCpvXY","CPV plane;X [cm];Y [cm]",200,-100.,100.,200,-100.,100.);
  Add2RawsList(hCpvXY, 5, !expert, !image, saveCorr);
  

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

}

//____________________________________________________________________________
void AliPMDQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
    //Fill prepared histograms with Raw digit properties

  TObjArray *pmdddlcont = 0x0;
    pmdddlcont = new TObjArray();
    AliPMDRawStream stream(rawReader);
    
    AliPMDddldata *pmdddl = 0x0;

    Int_t   iddl = -1;
    Int_t   xpad = -1;
    Int_t   ypad = -1;
    Float_t xx, yy;

    AliPMDUtility cc;

    while ((iddl = stream.DdlData(pmdddlcont)) >=0)
    {
	Int_t ientries = pmdddlcont->GetEntries();
	//printf(" ======= DDLNO = %d ientries = %d \n", iddl, ientries);
	
	for (Int_t ient = 0; ient < ientries; ient++)
	  {
	      //AliPMDddldata *pmdddl = (AliPMDddldata*)pmdddlcont->UncheckedAt(ient);
	    pmdddl = (AliPMDddldata*)pmdddlcont->UncheckedAt(ient);

	    Int_t det = pmdddl->GetDetector();
	    Int_t smn = pmdddl->GetSMN();
	    Int_t mcm = pmdddl->GetMCM();
	    //Int_t chno = pmdddl->GetChannel();
	    Int_t row = pmdddl->GetRow();
	    Int_t col = pmdddl->GetColumn();
	    Int_t sig = pmdddl->GetSignal();
	    
	    if (mcm == 0) continue;
	    if (det < 0 || det > 1)  continue;
	    if (smn < 0 || smn > 23) continue;
	    if (row < 0 || row > 47) continue;
	    if (col < 0 || col > 95) continue;

	    Int_t ipp = 10000*smn + 100*row + col;

	    if (det == 0)
	    {
	      GetRawsData(0)->Fill(ipp, sig);
	      GetRawsData(2)->Fill(ipp);

		if(smn < 12)
		{
		    xpad = col;
		    ypad = row;
		}
		else if(smn >= 12 && smn < 24)
		{
		    xpad = row;
		    ypad = col;
		}


		cc.RectGeomCellPos(smn,xpad,ypad,xx,yy);
		GetRawsData(4)->Fill(xx,yy);
	    }
	    if (det == 1)
	    {
	      GetRawsData(1)->Fill(ipp, sig);
	      GetRawsData(3)->Fill(ipp);
	
		if(smn < 12)
		{
		    xpad = col;
		    ypad = row;
		}
		else if(smn >= 12 && smn < 24)
		{
		    xpad = row;
		    ypad = col;
		}

		cc.RectGeomCellPos(smn,xpad,ypad,xx,yy);
		GetRawsData(5)->Fill(xx,yy);

	    }

	  }

	pmdddlcont->Delete();
    }

    delete pmdddlcont;
    pmdddlcont = 0x0;

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
	    GetDigitsData(0)->Fill( digit->GetADC()) ;
	    premul++;
      }
    if(digit->GetDetector() == 1)
      {
	    GetDigitsData(1)->Fill( digit->GetADC());
	    cpvmul++;
      }
    }  
  
  if (premul > 0) GetDigitsData(2)->Fill(premul);
  if (cpvmul > 0) GetDigitsData(3)->Fill(cpvmul);
  
  
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
  branch->SetAddress(&fDigitsArray) ;
  
  if ( ! branch )
    {
    AliWarning("PMD branch in Digit Tree not found") ; 
    }
  else
    {
    for (Int_t ient = 0; ient < branch->GetEntries(); ient++)
      {
	    branch->GetEntry(ient) ; 
	    MakeDigits() ; 
      }
    
    }
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
    branch->SetAddress(&fRecPointsArray) ;

    if ( ! branch )
    {
	AliWarning("PMD branch in Recpoints Tree not found") ; 
    }
    else
    {
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
			GetRecPointsData(0)->Fill(recpoint->GetClusCells());
			multDdl0++;
		      }
		    if(recpoint->GetSMNumber() >= 6 && recpoint->GetSMNumber() < 12)
		      {
			GetRecPointsData(1)->Fill(recpoint->GetClusCells());
			multDdl1++;
		      }
		    if(recpoint->GetSMNumber() >= 12 && recpoint->GetSMNumber() < 18)
		      {
			GetRecPointsData(2)->Fill(recpoint->GetClusCells());
			multDdl2++;
		      }
		    if(recpoint->GetSMNumber() >= 18 && recpoint->GetSMNumber() < 24)
		      {
			GetRecPointsData(3)->Fill(recpoint->GetClusCells());
			multDdl3++;
		      }
		  }

		if(recpoint->GetDetector() == 1)
		  {
		    if((recpoint->GetSMNumber() >= 0 && recpoint->GetSMNumber() < 6) || 
		       (recpoint->GetSMNumber() >= 18 && recpoint->GetSMNumber() < 24))
		      {
			GetRecPointsData(4)->Fill(recpoint->GetClusCells());
			multDdl4++;
		      }
		    if(recpoint->GetSMNumber() >= 6 && recpoint->GetSMNumber() < 18 )
		      {
			GetRecPointsData(5)->Fill(recpoint->GetClusCells());
			multDdl5++;
		      }
		  }
	      } 
	}
    }
    
    GetRecPointsData(6)->Fill(multDdl0,multDdl1);
    GetRecPointsData(7)->Fill(multDdl2,multDdl3);
    GetRecPointsData(8)->Fill(multDdl4,multDdl5);
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
	  GetESDsData(0)->Fill(adc);
	  premul++;
	}
      if (pmdtr->GetDetector() == 1)
	{
	  GetESDsData(1)->Fill(adc) ;
	  cpvmul++;
	}
    }
  
  GetESDsData(2)->Fill(cpvmul,premul) ;
  //GetESDsData(3)->Fill(cpvmul) ;  
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
