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
#include "AliPMDrecpoint1.h" 
#include "AliPMDRawStream.h"
#include "AliPMDddldata.h"
#include "AliPMDUtility.h"
#include "AliESDPmdTrack.h"
//#include "AliPMDRecoParam.h"

ClassImp(AliPMDQADataMakerRec)
           
//____________________________________________________________________________ 
  AliPMDQADataMakerRec::AliPMDQADataMakerRec() : 
  AliQADataMakerRec(AliQA::GetDetName(AliQA::kPMD), "PMD Quality Assurance Data Maker")
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
                                                                                                                                           
  TH1F * h0 = new TH1F("hPreDdl0Edep","DDL0 Energy Distribution PRE(PMD)", 100, 0, 2000);
  h0->Sumw2();                                                                                                                           
  Add2RawsList(h0, 0);                                                                                                             

  TH1F * h1 = new TH1F("hPreDdl1Edep","DDL1 Energy Distribution PRE(PMD)", 100, 0, 2000);
  h1->Sumw2();
  Add2RawsList(h1, 1);                                                                                                             

  TH1F * h2 = new TH1F("hPreDdl2Edep","DDL2 Energy Distribution PRE(PMD)", 100, 0, 2000);
  h2->Sumw2();                                                                                                                           
  Add2RawsList(h2, 2);                                                                                                             

  TH1F * h3 = new TH1F("hPreDdl3Edep","DDL3 Energy Distribution PRE(PMD)", 100, 0, 2000);
  h3->Sumw2();                                                                                                                       
  Add2RawsList(h3, 3);                                                                                                             

  TH1F * h4 = new TH1F("hCpvDdl4Edep","DDL4 Energy Distribution CPV(PMD)", 100, 0, 2000);
  h4->Sumw2();                                                                                                                           
  Add2RawsList(h4, 4);                                                                                                             

  TH1F * h5 = new TH1F("hCpvDdl5Edep","DDL5 Energy Distribution CPV(PMD)", 100, 0, 2000);
  h5->Sumw2();                                                                                                                           
  Add2RawsList(h5, 5);                                                                                                             

  
}
//____________________________________________________________________________ 
void AliPMDQADataMakerRec::InitRecPoints()
{
  // create Reconstructed Points histograms in RecPoints subdir

/*
    TH2F * h0 = new TH2F("hPreXY","RecPoints Y vs X PRE(PMD)", 100,-100.,100.,100,-100.,100.);
  Add2RecPointsList(h0,0) ;
 
    TH2F * h1 = new TH2F("hCpvXY","RecPoints Y vs X CPV(PMD)", 100,-100.,100.,100,-100.,100.);
  Add2RecPointsList(h1,1) ;
*/

    TH1F * h0 = new TH1F("hPreClADC","Cluster ADC of PRE(PMD) plane",500,0.,5000.);
    h0->Sumw2();
    Add2RecPointsList(h0, 0)  ;

    TH1F * h1 = new TH1F("hCpvClADC","Cluster ADC of CPV(PMD) plane",500,0.,5000.);
    h1->Sumw2();
    Add2RecPointsList(h1, 1)  ;
    

    TH1I * h2 = new TH1I("hPreClMult","Cluster Multiplicity of PRE(PMD) plane",100,0.,10000.);
    h2->Sumw2();
    Add2RecPointsList(h2, 2)  ;

    TH1I * h3 = new TH1I("hCpvClMult","Cluster Multiplicity of CPV(PMD) plane",100,0.,10000.);
    h3->Sumw2();
    Add2RecPointsList(h3, 3)  ;


}

//____________________________________________________________________________ 

void AliPMDQADataMakerRec::InitESDs()
{
  //Create histograms to controll ESD
 
  TH1F * h0 = new TH1F("hPreClADC","Cluster ADC of PRE(PMD) plane",500,0.,5000.);
  h0->Sumw2();
  Add2ESDsList(h0, 0)  ;

  TH1F * h1 = new TH1F("hCpvClADC","Cluster ADC of CPV(PMD) plane",500,0.,5000.);
  h1->Sumw2();
  Add2ESDsList(h1, 1)  ;

  TH1I * h2 = new TH1I("hPreClMult","Cluster Multiplicity of PRE(PMD) plane",100,0.,10000.);
  h2->Sumw2();
  Add2ESDsList(h2, 2)  ;

  TH1I * h3 = new TH1I("hCpvClMult","Cluster Multiplicity of CPV(PMD) plane",100,0.,10000.);
  h3->Sumw2();
  Add2ESDsList(h3, 3)  ;

}

//____________________________________________________________________________
void AliPMDQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
    //Fill prepared histograms with Raw digit properties

    TObjArray pmdddlcont;
    AliPMDRawStream stream(rawReader);
    


    Int_t iddl = -1;
    while ((iddl = stream.DdlData(&pmdddlcont)) >=0)
    {
	Int_t ientries = pmdddlcont.GetEntries();
	//printf(" ======= DDLNO = %d ientries = %d \n", iddl, ientries);
	for (Int_t ient = 0; ient < ientries; ient++)
	{
	    AliPMDddldata *pmdddl = (AliPMDddldata*)pmdddlcont.UncheckedAt(ient);
	    //Int_t det = pmdddl->GetDetector();
	    //Int_t smn = pmdddl->GetSMN();
	    //Int_t mcm = pmdddl->GetMCM();
	    //Int_t chno = pmdddl->GetChannel();
	    //Int_t row = pmdddl->GetRow();
	    //Int_t col = pmdddl->GetColumn();
	    Int_t sig = pmdddl->GetSignal();
	    //cout<<sig<<endl;
	    
	    if (iddl == 0) GetRawsData(0)->Fill(sig); 
	    if (iddl == 1) GetRawsData(1)->Fill(sig); 
	    if (iddl == 2) GetRawsData(2)->Fill(sig); 
	    if (iddl == 3) GetRawsData(3)->Fill(sig); 
	    if (iddl == 4) GetRawsData(4)->Fill(sig); 
	    if (iddl == 5) GetRawsData(5)->Fill(sig); 
	    
	    
	}
	
	pmdddlcont.Delete();
    }
}
//____________________________________________________________________________
void AliPMDQADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
    // makes data from RecPoints

    Int_t premul = 0, cpvmul = 0;
    AliPMDrecpoint1 * recpoint; 

    TClonesArray * recpoints = new TClonesArray("AliPMDrecpoint1", 1000) ; 
    
    TBranch * branch = clustersTree->GetBranch("PMDRecpoint") ;
    branch->SetAddress(&recpoints) ;

    if ( ! branch )
    {
	AliWarning("PMD branch in SDigit Tree not found") ; 
    }
    else
    {
	for (Int_t imod = 0; imod < branch->GetEntries(); imod++)
	{
	    branch->GetEntry(imod) ;

	    TIter next(recpoints) ; 

	    while ( (recpoint = dynamic_cast<AliPMDrecpoint1 *>(next())) )
	    {
		//Float_t xpos = recpoint->GetClusX();
		//Float_t ypos = recpoint->GetClusY();
		if(recpoint->GetDetector() == 0)
		{
		    GetRecPointsData(0)->Fill(recpoint->GetClusADC());
		    premul++;
		}
		if(recpoint->GetDetector() == 1)
		{
		    GetRecPointsData(1)->Fill(recpoint->GetClusADC());
		    cpvmul++;
		}
	
	    } 

	}
    }

    GetRecPointsData(2)->Fill(premul);
    GetRecPointsData(3)->Fill(cpvmul);
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

    GetESDsData(2)->Fill(premul) ;
    GetESDsData(3)->Fill(cpvmul) ;
}

//____________________________________________________________________________ 

void AliPMDQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}
//____________________________________________________________________________ 
void AliPMDQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQA::kPMD, task, list) ;  
}
