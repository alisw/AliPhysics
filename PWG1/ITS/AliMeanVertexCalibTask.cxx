/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//*************************************************************************
// Class AliMeanVertexCalibTask
// AliAnalysisTask to extract from ESD the information on primary vertex
// reconstruction in order to compute the MeanVertex object
//
// Author:  D.Caffarri, davide.caffarri@pd.infn.it  
//          A.Dainese, andrea.dainese@pd.infn.it
//*************************************************************************


#include <TH1F.h>
#include <TH2F.h>

#include "AliMultiplicity.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliVertexerTracks.h"

#include "AliMeanVertexCalibTask.h"


ClassImp(AliMeanVertexCalibTask)

//_____________________________________________________________________
AliMeanVertexCalibTask::AliMeanVertexCalibTask(const char *name):
  AliAnalysisTaskSE(name),
  fESD(0), 
  fOutput(0),
  fOnlyITSTPCTracks(kFALSE),
  fOnlyITSSATracks(kTRUE)
{

  // Constructor
  
  // Define input and output slots here
  // Output slot #0 writes into a TList container
  DefineOutput(1, TList::Class());  //My private output
}

//________________________________________________________________________
AliMeanVertexCalibTask::~AliMeanVertexCalibTask()
{
  // Destructor

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  
}
//________________________________________________________________________
void AliMeanVertexCalibTask::UserCreateOutputObjects()
{

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList;
  fOutput->SetOwner();

  
  TH1F* hSPDVertexX = new TH1F("hSPDVertexX","SPDVertex x; x vertex [cm]; events",200,-1,1);
  fOutput->Add(hSPDVertexX);
  TH1F* hSPDVertexY = new TH1F("hSPDVertexY","SPDVertex y; y vertex [cm]; events",200,-1,1);
  fOutput->Add(hSPDVertexY);
  TH1F* hSPDVertexZ = new TH1F("hSPDVertexZ","SPDVertex z; z vertex [cm]; events",200,-20,20);
  fOutput->Add(hSPDVertexZ);
  TH1F* hTRKVertexX = new TH1F("hTRKVertexX","TRKVertex x; x vertex [cm]; events",200,-1,1);
  fOutput->Add(hTRKVertexX);
  TH1F* hTRKVertexY = new TH1F("hTRKVertexY","TRKVertex y; y vertex [cm]; events",200,-1,1);
  fOutput->Add(hTRKVertexY);
  TH1F* hTRKVertexZ = new TH1F("hTRKVertexZ","TRKVertex z; z vertex [cm]; events",200,-20,20);
  fOutput->Add(hTRKVertexZ);
  
  TH2F *hTRKVertexXZ = new TH2F("hTRKVertexXZ", "TRKVertex XZ corr", 200, -1, 1, 200, -20, 20);
  fOutput->Add(hTRKVertexXZ);

  TH2F *hTRKVertexYZ = new TH2F("hTRKVertexYZ", "TRKVertex YZ corr", 200, -1, 1, 200, -20, 20);
  fOutput->Add(hTRKVertexYZ);
  
  TH1F* hTRKVertexXdefMult = new TH1F("hTRKVertexXdefMult","TRKVertex x Mult; x vertex [cm] 30<Mult<45; events",500,-1,1);  
  fOutput->Add(hTRKVertexXdefMult);
  TH1F* hTRKVertexYdefMult = new TH1F("hTRKVertexYdefMult","TRKVertex y Mult; y vertex [cm] 30<Mult<45; events",500,-1,1);
  fOutput->Add(hTRKVertexYdefMult);
	
	TH1F* hTRKVertexXHighMult = new TH1F("hTRKVertexXHighMult","TRKVertex x High Mult; x vertex [cm] Mult>1500; events",500,-0.5,0.5);  
	fOutput->Add(hTRKVertexXHighMult);
	TH1F* hTRKVertexYHighMult = new TH1F("hTRKVertexYHighMult","TRKVertex y High Mult; y vertex [cm] Mult>1500; events",500,-0.5,0.5);
	fOutput->Add(hTRKVertexYHighMult);	
	

  TH1F* hITSSAVertexX = new TH1F("hITSSAVertexX","ITSSAVertex x; x vertex [cm]; events",200,-1,1);
  fOutput->Add(hITSSAVertexX);
  TH1F* hITSSAVertexY = new TH1F("hITSSAVertexY","ITSSAVertex y; y vertex [cm]; events",200,-1,1);
  fOutput->Add(hITSSAVertexY);
  TH1F* hITSSAVertexZ = new TH1F("hITSSAVertexZ","ITSSAVertex z; z vertex [cm]; events",200,-20,20);
  fOutput->Add(hITSSAVertexZ);
  
  TH2F *hITSSAVertexXZ = new TH2F("hITSSAVertexXZ", "ITSSAVertex XZ corr", 200, -1, 1, 200, -20, 20);
  fOutput->Add(hITSSAVertexXZ);

  TH2F *hITSSAVertexYZ = new TH2F("hITSSAVertexYZ", "ITSSAVertex YZ corr", 200, -1, 1, 200, -20, 20);
  fOutput->Add(hITSSAVertexYZ);
  
  TH1F* hITSSAVertexXdefMult = new TH1F("hITSSAVertexXdefMult","ITSSAVertex x Mult; x vertex [cm] 30<Mult<45; events",500,-1,1);  
  fOutput->Add(hITSSAVertexXdefMult);
  TH1F* hITSSAVertexYdefMult = new TH1F("hITSSAVertexYdefMult","ITSSAVertex y Mult; y vertex [cm] 30<Mult<45; events",500,-1,1);
  fOutput->Add(hITSSAVertexYdefMult);

	
   TH1F* hITSSAVertexXHighMult = new TH1F("hITSSAVertexXHighMult","ITSSAVertex x High Mult; x vertex [cm] Mult>1500; events",500,-0.5,0.5);  
   fOutput->Add(hITSSAVertexXHighMult);
   TH1F* hITSSAVertexYHighMult = new TH1F("hITSSAVertexYHighMult","ITSSAVertex y High Mult; y vertex [cm] Mult>1500; events",500,-0.5,0.5);
   fOutput->Add(hITSSAVertexYHighMult);
	
	

  PostData(1, fOutput);
  
  return;
}


//________________________________________________________________________
void AliMeanVertexCalibTask::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  if (!InputEvent()) {
    Printf("ERROR: fESD not available");
    return;
  }
  
  AliESDEvent* esdE = (AliESDEvent*) InputEvent();

  const AliMultiplicity *alimult = esdE->GetMultiplicity();
  Int_t ntrklets=0;
  if(alimult) ntrklets = alimult->GetNumberOfTracklets();
   

  const AliESDVertex *trkv = esdE->GetPrimaryVertexTracks();
  const AliESDVertex *spdv = esdE->GetPrimaryVertexSPD();

  AliESDVertex *itsSAv = ReconstructPrimaryVertex(kFALSE); 

  if(spdv) {
    if(spdv->GetNContributors()>0) {
      TString title=spdv->GetTitle();
      if(title.Contains("3D")) {
	((TH1F*)fOutput->FindObject("hSPDVertexX"))->Fill(spdv->GetXv());
	((TH1F*)fOutput->FindObject("hSPDVertexY"))->Fill(spdv->GetYv());
      }
      ((TH1F*)fOutput->FindObject("hSPDVertexZ"))->Fill(spdv->GetZv());
    }
  }
  
  
  if(trkv) {
    if(trkv->GetNContributors()>0) {
      ((TH1F*)fOutput->FindObject("hTRKVertexX"))->Fill(trkv->GetXv());
      ((TH1F*)fOutput->FindObject("hTRKVertexY"))->Fill(trkv->GetYv());
      ((TH1F*)fOutput->FindObject("hTRKVertexZ"))->Fill(trkv->GetZv());
      
      if (ntrklets>30 && ntrklets<45){
	((TH1F*)fOutput->FindObject("hTRKVertexXdefMult"))->Fill(trkv->GetXv());
	((TH1F*)fOutput->FindObject("hTRKVertexYdefMult"))->Fill(trkv->GetYv());
      }
		
	  if (ntrklets>1500){
			((TH1F*)fOutput->FindObject("hTRKVertexXHighMult"))->Fill(trkv->GetXv());
			((TH1F*)fOutput->FindObject("hTRKVertexYHighMult"))->Fill(trkv->GetYv());
		}
      
      ((TH2F*)fOutput->FindObject("hTRKVertexXZ"))->Fill(trkv->GetXv(),trkv->GetZv());
      ((TH2F*)fOutput->FindObject("hTRKVertexYZ"))->Fill(trkv->GetYv(),trkv->GetZv());
      
    }
  }

  if (itsSAv){
    if (itsSAv->GetNContributors()>0){
  
      ((TH1F*)fOutput->FindObject("hITSSAVertexX"))->Fill(itsSAv->GetXv());
      ((TH1F*)fOutput->FindObject("hITSSAVertexY"))->Fill(itsSAv->GetYv());
      ((TH1F*)fOutput->FindObject("hITSSAVertexZ"))->Fill(itsSAv->GetZv());
      
		if (ntrklets>30 && ntrklets<45){
			((TH1F*)fOutput->FindObject("hITSSAVertexXdefMult"))->Fill(itsSAv->GetXv());
		    ((TH1F*)fOutput->FindObject("hITSSAVertexYdefMult"))->Fill(itsSAv->GetYv());
		}
		
		if (ntrklets>1500){
			((TH1F*)fOutput->FindObject("hITSSAVertexXHighMult"))->Fill(itsSAv->GetXv());
			((TH1F*)fOutput->FindObject("hITSSAVertexYHighMult"))->Fill(itsSAv->GetYv());
		}
      
      ((TH2F*)fOutput->FindObject("hITSSAVertexXZ"))->Fill(itsSAv->GetXv(),itsSAv->GetZv());
      ((TH2F*)fOutput->FindObject("hITSSAVertexYZ"))->Fill(itsSAv->GetYv(),itsSAv->GetZv());
     
    }
  }

  delete itsSAv;
  
  PostData(1, fOutput);
  
  return;
  
}


//________________________________________________________________________
void AliMeanVertexCalibTask::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    Printf("ERROR: fOutput not available");
    return;
  }


  return;

}


//__________________________________________________________________________
AliESDVertex* AliMeanVertexCalibTask::ReconstructPrimaryVertex(Bool_t constr,Int_t mode) const {
  // On the fly reco of ITS+TPC vertex from ESD
  // mode=0 use all tracks
  // mode=1 use odd-number tracks
  // mode=2 use even-number tracks

  AliESDEvent* evt = (AliESDEvent*) fInputEvent;
  AliVertexerTracks vertexer(evt->GetMagneticField());
  if(evt->GetNumberOfTracks()<500) {
    vertexer.SetITSMode(); // defaults
    vertexer.SetMinClusters(4); // default is 5
  } else { 
    vertexer.SetITSMode(0.1,0.1,0.5,5,1,3.,100.,1000.,3.,30.,1,1);// PbPb
  } 
  if (fOnlyITSSATracks) vertexer.SetITSpureSA();
  
  Float_t diamondcovxy[3]; evt->GetDiamondCovXY(diamondcovxy);
  Double_t pos[3]={evt->GetDiamondX(),evt->GetDiamondY(),0}; 
  Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
  AliESDVertex *initVertex = new AliESDVertex(pos,cov,1.,1);
  vertexer.SetVtxStart(initVertex);
  delete initVertex;
  if(!constr) vertexer.SetConstraintOff();

  if(fOnlyITSTPCTracks || fOnlyITSSATracks || mode>0) {
    Int_t iskip=0;
    Int_t *skip = new Int_t[evt->GetNumberOfTracks()];
    for(Int_t itr=0;itr<evt->GetNumberOfTracks(); itr++) {
      AliESDtrack* track = evt->GetTrack(itr);
      if(fOnlyITSTPCTracks && track->GetNcls(1)==0) { // skip ITSSA
	skip[iskip++]=itr;
	continue;
      }
      if(fOnlyITSSATracks && track->GetNcls(1)>0) { // skip ITSTPC
	skip[iskip++]=itr;
	continue;
      }
      if(mode==1 && itr%2==0) skip[iskip++]=itr;
      if(mode==2 && itr%2==1) skip[iskip++]=itr;
    }
    vertexer.SetSkipTracks(iskip,skip);
    delete [] skip; skip=NULL;
  }
  
  return vertexer.FindPrimaryVertex(evt);
}
