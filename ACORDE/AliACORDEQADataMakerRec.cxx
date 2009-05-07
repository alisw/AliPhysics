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
//---
//  Produces the data needed to calculate the quality assurance. 
//  All data must be mergeable objects.


//  Authors:
//
//  Luciano Diaz Gonzalez <luciano.diaz@nucleares.unam.mx> (ICN-UNAM)
//  Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch> (FCFM-BUAP)
//  Arturo Fernandez Tellez <afernan@mail.cern.ch (FCFM-BUAP)
//
//  Created: June 13th 2008
//---
// Last Update: Aug. 27th 2008 --> Implementation to declare QA expert histogram


// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TDirectory.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliACORDEdigit.h" 
#include "AliACORDEhit.h"
#include "AliACORDEQADataMakerRec.h"
#include "AliQAChecker.h"
#include "AliACORDERawReader.h"
#include "AliACORDERawStream.h"
ClassImp(AliACORDEQADataMakerRec)
           
//____________________________________________________________________________ 
AliACORDEQADataMakerRec::AliACORDEQADataMakerRec():AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kACORDE), "ACORDE Quality Assurance Data Maker")
{

}
//____________________________________________________________________________ 
AliACORDEQADataMakerRec::AliACORDEQADataMakerRec(const AliACORDEQADataMakerRec& qadm):AliQADataMakerRec() 
{
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}
//__________________________________________________________________
AliACORDEQADataMakerRec& AliACORDEQADataMakerRec::operator = (const AliACORDEQADataMakerRec& qadm )
{
  // Equal operator.
  this->~AliACORDEQADataMakerRec();
  new(this) AliACORDEQADataMakerRec(qadm);
  return *this;
}
//____________________________________________________________________________
void AliACORDEQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kACORDE, task, list) ;
}

//____________________________________________________________________________
void AliACORDEQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle

}
 
//____________________________________________________________________________ 
void AliACORDEQADataMakerRec::InitRaws()
{
  // create Raw histograms in Raw subdir

TH1D *fhACORDEBitPattern[4];
fhACORDEBitPattern[0] = new TH1D("ACORDERawDataSM","ACORDE-SingleMuon",60,1,60);//AcordeSingleMuon BitPattern
fhACORDEBitPattern[1] = new TH1D("ACORDERawDataMM","ACORDE-MultiMuon",60,1,60);//AcordeMultiMuon BitPattern
fhACORDEBitPattern[2] = new TH1D("ACORDERawDataSMM","ACORDE-SingleMuonMultiplicity",60,1,60);//AcordeSingleMuon Multiplicity
fhACORDEBitPattern[3] = new TH1D("ACORDERawDataMMM","ACORDE-MultiMuonMultiplicity",60,1,60);//AcordeMultiMuon Multiplicity
for(Int_t i=0;i<4;i++)
{
	Add2RawsList(fhACORDEBitPattern[i],i,kFALSE);
}
}
//____________________________________________________________________________ 

void AliACORDEQADataMakerRec::InitRecPoints()
{
  // create cluster histograms in RecPoint subdir
  // Not needed for ACORDE by now !!!
}
//____________________________________________________________________________
void AliACORDEQADataMakerRec::InitESDs()
{
  //create ESDs histograms in ESDs subdir

   TH1F *    fhESDsSingle;
   TH1F *    fhESDsMulti;

   TString   name;

   name = "hESDsSingle";
   fhESDsSingle = new TH1F(name.Data(),"hESDsSingle",60,0,60);
   Add2ESDsList(fhESDsSingle,0,kFALSE);

   name = "hESDsMulti";
   fhESDsMulti = new TH1F(name.Data(),"hESDsMulti",60,0,60);
   Add2ESDsList(fhESDsMulti,1,kFALSE);
}
//____________________________________________________________________________
void AliACORDEQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  //fills QA histos for RAW
  rawReader->Reset();
  AliACORDERawStream rawStream(rawReader);
  size_t contSingle=0;
  size_t contMulti=0;
  UInt_t dy[4];

  bool kroSingle[60],kroMulti[60];
  UInt_t tmpDy;

  for(Int_t m=0;m<60;m++) {kroSingle[m]=0;kroMulti[m]=0;}

if(rawStream.Next())
{
        dy[0]=rawStream.GetWord(0);
        dy[1]=rawStream.GetWord(1);
        dy[2]=rawStream.GetWord(2);
        dy[3]=rawStream.GetWord(3);
        tmpDy=dy[0];
        for(Int_t r=0;r<30;++r)
        {
                kroSingle[r] = tmpDy & 1;
                tmpDy>>=1;
        }
        tmpDy=dy[1];
        for(Int_t r=30;r<60;++r)
        {
                kroSingle[r] = tmpDy & 1;
                tmpDy>>=1;
        }
        tmpDy=dy[2];
        for(Int_t r=0;r<30;++r)
        {
                kroMulti[r] = tmpDy & 1;
                tmpDy>>=1;
        }
        tmpDy=dy[3];
        for(Int_t r=30;r<60;++r)
        {
                kroMulti[r] = tmpDy & 1;
                tmpDy>>=1;
        }
        contSingle=0;
	contMulti=0;
        for(Int_t r=0;r<60;++r)
        {
                if(kroSingle[r]==1)
                {
                        GetRawsData(0)->Fill(r+1);
                        contSingle=contSingle+1;
                }
		if(kroMulti[r]==1)
		{
			GetRawsData(1)->Fill(r+1);
			contMulti++;
		}

        }GetRawsData(2)->Fill(contSingle);GetRawsData(3)->Fill(contMulti);
}
}
//____________________________________________________________________________
void AliACORDEQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  //fills QA histos for ESD

        AliESDACORDE * fESDACORDE= esd->GetACORDEData();
        Int_t *fACORDEMultiMuon =fESDACORDE->GetACORDEMultiMuon();
        Int_t *fACORDESingleMuon=fESDACORDE->GetACORDESingleMuon();

       for(int i=0;i<60;i++){
            if(fACORDESingleMuon[i]==1)
               GetESDsData(0) -> Fill(i);
            if(fACORDEMultiMuon[i]==1)
               GetESDsData(1) -> Fill(i);
        }




}

