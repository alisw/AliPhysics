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

// This class reads rawdata for MUON detector and produces  
// reconstructed hit position as output, which is stored in a root 
// file  as tree format. It uses another class HLTMUONHitReconstructor of "../PubSub"
// to perform the hitreconstruction which is the same as it is used 
// in online hitreconstruction. The output data is of type AliMUONRawCluster same as that
// of produced in Dimuon reconstruction (i.e MUON.RecPoints.root) 
 
///////////////////////////////////////////////
//Author : Sukalyan Chattopadhyay, SINP, INDIA
//         Indranil Das, SINP, INDIA
//
//Email :  indra.das@saha.ac.in
//         sukalyan.chattopadhyay@saha.ac.in 
///////////////////////////////////////////////


#include "AliHLTMUONHitReconstructor.h"

#include <TTree.h>
#include <TFile.h>
#include <TClonesArray.h>

#include "AliLog.h"
#include "AliRawReader.h"
//#include "AliDAQ.h"

#include "AliMUONRawCluster.h"



ClassImp(AliHLTMUONHitReconstructor)


AliHLTMUONHitReconstructor::AliHLTMUONHitReconstructor(AliRawReader* rawReader): 
  TObject(),
  fRawReader(rawReader),
  fLutPath(0)
{
  // ctor 
  f1 = new TFile("DHLTRecPionts.root","recreate");
  f1->cd();

  fRecPoints = new TClonesArray("AliMUONRawCluster",1024);

  fRecTimer.Start();

}


AliHLTMUONHitReconstructor:: AliHLTMUONHitReconstructor(const AliHLTMUONHitReconstructor& rhs)
  : TObject(rhs)
{
// Protected copy constructor

  AliFatal("Not implemented.");
}


AliHLTMUONHitReconstructor & 
AliHLTMUONHitReconstructor::operator=(const AliHLTMUONHitReconstructor& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}

AliHLTMUONHitReconstructor::~AliHLTMUONHitReconstructor()
{
  // Deletes all dHLT objects created in the constructor.
  fRecTimer.Stop();
  AliInfo(Form("Total time for  RecPoints : R:%.2fs C:%.2fs",
               fRecTimer.RealTime(),fRecTimer.CpuTime()));
  fRecTimer.Delete();

  if(fRecPoints){
    fRecPoints->Delete();
    delete fRecPoints ;
  }

  for(int i=0; i<8; i++){
    if( fLut[i] != NULL)
      delete [] fLut[i];
  }

  delete f1;
}


Bool_t AliHLTMUONHitReconstructor::Init(const char* lutpath, const char* buspatchmappath)
{
  // Initialisation of this
  if(!ReadLookUpTable(lutpath)){
    return kFALSE;
  }
  
  BusToDetElem busToDetElem;
  if (! ReadBusPatchToDetElemFile(busToDetElem, buspatchmappath))
    {
      return kFALSE;
    }

  if (! fkHLTRec.SetBusToDetMap(busToDetElem))
    {
      AliInfo(Form("Failed to load but patch to detector element mapping."));
      return kFALSE;
    }

  return kTRUE;
}

Bool_t AliHLTMUONHitReconstructor::ReadLookUpTable(const char* lutpath)
{
  // Read all the LookupTables
  for(int iDDL = 0 ; iDDL < HLTMUONHitReconstructor::fgkNofDDL ; iDDL++){
    
    int ddlId = iDDL +  HLTMUONHitReconstructor::fgkDDLOffSet;
    

    int lutLine = fkHLTRec.GetLutLine(ddlId);
    fLut[iDDL] = new DHLTLut[lutLine];  
    
    char* filepath = new char[strlen(lutpath)+10];
    sprintf(filepath,"%s/Lut%d.dat",lutpath,ddlId);
    
    FILE* fin = fopen(filepath, "r");
    if (fin == NULL)
      {
      AliInfo(Form("Failed to open file: %s\n",filepath));
      delete [] filepath;
      delete [] fLut[iDDL];
      return kFALSE;
      }
    
    for(int i=0;i<lutLine;i++)
      {
	fscanf(
	       fin,
	       "%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d",
	       &fLut[iDDL][i].fIdManuChannel,
	       &fLut[iDDL][i].fIX,
	       &fLut[iDDL][i].fIY,
	       &fLut[iDDL][i].fRealX,
	       &fLut[iDDL][i].fRealY,
	       &fLut[iDDL][i].fRealZ,
	       &fLut[iDDL][i].fPcbZone,
	       &fLut[iDDL][i].fPlane
	       );
      }
  
    fclose(fin);
    delete [] filepath;
    
    
  }//iDDL loop


//   if (! fkHLTRec.LoadLookUpTable(lut[iDDL], iDDL))
//     {
//       AliInfo(Form("Failed to load lookup table."));
//       delete [] lut[iDDL];
//       return kFALSE;
//     }
//   delete [] lut[iDDL];
  
  return kTRUE;
}


Bool_t AliHLTMUONHitReconstructor::ReadBusPatchToDetElemFile(BusToDetElem& busToDetElem, 
							     const char* buspatchmappath)
{
  // Read the bupatch to detelem mapping
  char getLine[80];
  char temp;
  int detElem, minBusPatch, maxBusPatch;

  char* filepath = new char[strlen(buspatchmappath)+20];
  sprintf(filepath,"%s/BusToDetElem.dat",buspatchmappath);

  FILE* fin = fopen(filepath, "r");
  if (fin == NULL)
    {
      AliInfo(Form("Failed to open file: %s\n",buspatchmappath));	  
      return kFALSE;
    }

  while (feof(fin)==0){
    
    fgets(getLine,80,fin);
    sscanf(getLine, "%d\t%d %c %d\n", &detElem, &minBusPatch, &temp, &maxBusPatch);
    if (detElem >= 700 && detElem <= 1025)
      {
	for(int i = minBusPatch; i <= maxBusPatch; i++)
	  busToDetElem[i] = detElem;
      } // detElem condn
  } // while loop for file
  
  fclose(fin);
  delete [] filepath;
  return kTRUE;
}


Bool_t AliHLTMUONHitReconstructor::WriteDHLTRecHits(Int_t iEvent) 
{  
//   // main function called by AliHLTReconstructor to perform DHLT Hitreconstruction 
  Char_t branchName[30],eventId[20];
  DHLTRecPoint hits[1024];
  int iDDLId;
  int nofHit = 1024;
  AliMUONRawCluster cnew;


  sprintf(eventId,"Event%d",iEvent);
  f1->cd();
  f1->mkdir(eventId);
  f1->cd(eventId);
  

  fDHLTTree = new TTree("TreeR","TreeR");

  

  //printf("fRawReader0 : %p\n",fRawReader);
  for(int iChamber = 6; iChamber < 10 ; iChamber++){

    sprintf(branchName,"HLTMUONRecPoints%d",iChamber+1);
    fDHLTTree->Branch(branchName,&fRecPoints,32000);

    for(int iDDL = 0; iDDL<2; iDDL++){
      iDDLId = iDDL + 2*iChamber;
      fRawReader->Reset();
      fRawReader->Select("MUONTRK",iDDLId,iDDLId);  //Select the DDL file to be read
      
      fRawReader->ReadHeader();
      Int_t totalDataSize = fRawReader->GetDataSize();
      Int_t* buffer = new Int_t[totalDataSize/4];
      fRawReader->ReadNext((UChar_t*)buffer,totalDataSize);
      
      fkHLTRec.LoadLookUpTable(fLut[(iDDLId - HLTMUONHitReconstructor::fgkDDLOffSet)], iDDLId);
      nofHit = 1024;
      fkHLTRec.Run(buffer,totalDataSize,hits,&nofHit);
      delete [] buffer;
      
      if(nofHit>0){
	TClonesArray &temp = *fRecPoints;  
	for(int i=0;i<nofHit;i++){
	  cnew.SetDetElemId(hits[i].DetElemId);
	  cnew.SetX(0,hits[i].X);
	  cnew.SetY(0,hits[i].Y);
	  cnew.SetZ(0,hits[i].Z);
	  new(temp[temp.GetEntriesFast()]) AliMUONRawCluster(cnew);
	}
      }

    }//
  }// iChamber loop

  fDHLTTree->Fill();
  fDHLTTree->Write();
  fDHLTTree->Delete();
  fRecPoints->Clear("C");
  
  return kTRUE;
}

