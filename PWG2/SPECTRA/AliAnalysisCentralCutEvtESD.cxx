/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

//  *******************************************
//  * ESD event level cuts for azimuthal isotropic  *
//  * expansion in highly central collisions analysis *
//  * author: Cristian Andrei                                    *
//  *         acristian@niham.nipne.ro                        *
//  * *****************************************

#include "TMath.h"
#include <TObjArray.h>
#include <TList.h>

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
 #include "AliMultiplicity.h"

#include "AliAnalysisCentralCutESD.h"
#include "AliAnalysisCentralCutEvtESD.h"

class TObject;


//____________________________________________________________________
ClassImp(AliAnalysisCentralCutEvtESD)

//____________________________________________________________________
AliAnalysisCentralCutEvtESD::AliAnalysisCentralCutEvtESD(const Char_t* name, const Char_t* title) 
    :AliAnalysisCuts(name,title)
    ,fReqMult(kFALSE)
    ,fReqDir(kFALSE)
    ,fReqSPDMult(kFALSE)
    ,fReqSPDDir(kFALSE)
    ,fMultMin(0)
    ,fMultMax(0)
    ,fDirMin(0)
    ,fDirMax(0)
    ,fSPDMultMin(0)
    ,fSPDMultMax(0)
    ,fSPDDirMin(0)
    ,fSPDDirMax(0)

{
//constructor
    for(Int_t i=0; i<10; i++){
		fCutsList[i] = 0;
    }

}

AliAnalysisCentralCutEvtESD::~AliAnalysisCentralCutEvtESD() {
//Destructor
	if(fCutsList){
		delete [] fCutsList;
	}

}



Bool_t AliAnalysisCentralCutEvtESD::IsSelected(TObject *obj){
// check whether the event passes the cuts
    AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(obj);

    if(!esdEvent){
	printf("AliAnalysisCentralCutEvtESD:IsSelected ->Can't get ESD event!\n");
	exit(1);
    }


    if(fReqMult){
		Int_t mult = CalcMult(esdEvent);
		if((mult<fMultMin)||(mult>fMultMax)){
			return kFALSE;
		}
    }



    if(fReqDir){
		Double_t dir = CalcDir(esdEvent);
		if((dir<fDirMin)||(dir>fDirMax)){
			return kFALSE;
		}
    }

    if(fReqSPDMult){
		Double_t spdMult = CalcSPDMult(esdEvent);
		if((spdMult<fSPDMultMin)||(spdMult>fSPDMultMax)){
			return kFALSE;
		}
    }

    if(fReqSPDDir){
		Double_t spdDir = CalcSPDDir(esdEvent);
		if((spdDir<fSPDDirMin)||(spdDir>fSPDDirMax)){
			return kFALSE;
		}
    }

    return kTRUE;
}

//_________________________________________________________________________
void AliAnalysisCentralCutEvtESD::InitCuts(){
//Initialize internal cuts

////////////////  FULL ALICE ///////////////////
//------------General ESD Cuts------------------
    AliESDtrackCuts *esdCutsGen = new AliESDtrackCuts("AliESDtrackCuts", "CutEvtESDInternal");
    esdCutsGen->SetMinNClustersTPC(50);
    esdCutsGen->SetMaxChi2PerClusterTPC(2.2);
    esdCutsGen->SetMaxCovDiagonalElements(0.5,0.5,0.5,0.5,0.5);
    esdCutsGen->SetRequireTPCRefit(kTRUE);
    esdCutsGen->SetAcceptKinkDaughters(kFALSE);
    esdCutsGen->SetMaxNsigmaToVertex(2.0);
    esdCutsGen->SetRequireSigmaToVertex(kTRUE);

    AliAnalysisCentralCutESD *esdCutsGen1 = new AliAnalysisCentralCutESD("AliAnalysisCentralCutESD","NIHAM");
    esdCutsGen1->SetReqIsCharged();

//-------------Specific ESD Cuts------------------
    AliESDtrackCuts *esdCutsMult = new AliESDtrackCuts("AliESDCutsMult", "midrapidity");
    esdCutsMult->SetEtaRange(-0.5,0.5);
    AliESDtrackCuts *esdCutsDir = new AliESDtrackCuts("AliESDCutsDir", "SPD coverage");
    esdCutsDir->SetEtaRange(0.0,1.9);

///////////////////  SPD ONLY  ////////////////////
    AliESDtrackCuts *esdCutsSPD = new AliESDtrackCuts("AliESDtrackCuts", "CutEvtESDInternal");
    esdCutsSPD->SetAcceptKinkDaughters(kFALSE);
    esdCutsSPD->SetMaxNsigmaToVertex(2.0);
    esdCutsSPD->SetRequireSigmaToVertex(kTRUE);


//--------------set the cuts ----------------------
    TObjArray* esdListMult = new TObjArray();
    esdListMult->AddLast(esdCutsGen);
    esdListMult->AddLast(esdCutsGen1);
    esdListMult->AddLast(esdCutsMult);

    TObjArray* esdListDir = new TObjArray();
    esdListDir->AddLast(esdCutsGen);
    esdListDir->AddLast(esdCutsGen1);
    esdListDir->AddLast(esdCutsDir);

    TObjArray* esdListSPD = new TObjArray();
    esdListSPD->AddLast(esdCutsSPD);

    fCutsList[0]=esdListDir;
    fCutsList[1]=esdListMult;
    fCutsList[2]=esdListSPD;

}


//__________________________________________________________________________
Bool_t AliAnalysisCentralCutEvtESD::CheckIntCuts(Int_t no, TObject *obj) const{
// Check if the particle passes the internal cuts

    if(no > 9){
		printf("\n AliAnalysisCentralCutEvtESD::CheckIntCuts -> Cut number is not ok! \n");
		return kFALSE;
    }

    if(!fCutsList[no]){
		printf("AliAnalysisCentralCutEvtESD::CheckIntCuts -> cuts list problem! \n");
		return kFALSE;
    }


    TObjArrayIter iter(fCutsList[no]);
    AliAnalysisCuts *cut = 0;


    while((cut = (AliAnalysisCuts*)iter.Next())){

		if(!cut->IsSelected(obj)) return kFALSE;
    }

    return kTRUE;
}


//__________________________________________________________________________
Double_t AliAnalysisCentralCutEvtESD::CalcDir(AliESDEvent* const esdEv) {

//Compute the directivity - FULL ALICE

    InitCuts();

    Double_t dir;
    Double_t px,py;
    Double_t sumaPt = 0;
    Double_t sumaPx = 0;
    Double_t sumaPy = 0;

    Double_t pt;

    if (!esdEv){
		printf("NULL tree\n");
		return -1;
    }

    Int_t totTracks=esdEv->GetNumberOfTracks();

    for(Int_t itrack = 0; itrack < totTracks; itrack++){//first loop->compute event directivity

		AliESDtrack* track = esdEv->GetTrack(itrack);
		
		if (!track) {
				Printf("ERROR: Could not receive track %d", itrack);
			continue;
		}
	
		if(!CheckIntCuts(0, track)) continue;
		
		px = track->Px();	
		py = track->Py();	
		pt = track->Pt();
	
		sumaPx = sumaPx + px;
		sumaPy = sumaPy + py;
	
		sumaPt = sumaPt + pt;

    }//end track loop


	if(sumaPt < 0.0000001){
		return -1;
    }

    dir = (sqrt(pow(sumaPx,2)+pow(sumaPy,2)))/sumaPt;

    return dir;
}

//__________________________________________________________________________
Int_t AliAnalysisCentralCutEvtESD::CalcMult(AliESDEvent* const esdEv) {

//Compute multiplicity - FULL ALICE

    InitCuts();

    Int_t charged = 0;

    if (!esdEv){
		printf("NULL tree\n");
		return -1;
    }

    Int_t totTracks=esdEv->GetNumberOfTracks();

    for(Int_t iTrack = 0; iTrack < totTracks; iTrack++){//second track loop -> compute event multiplicity

		AliESDtrack* track = esdEv->GetTrack(iTrack);
		
		if (!track) {
			Printf("ERROR: Could not receive track %d", iTrack);
			continue;
		}
	
		if(!CheckIntCuts(1, track)) continue;
	
		
		charged++; //multiplicity

    }//end second track loop


    return charged;
}


//__________________________________________________________________________
Double_t AliAnalysisCentralCutEvtESD::CalcSPDDir(AliESDEvent* const esdEv) {

//Compute directivity - SPD ONLY

    InitCuts();

    Double_t dirU;
    Double_t pxU,pyU;
    Double_t sumaPxU = 0;
    Double_t sumaPyU = 0;

    Double_t goodtrack = 0;
    Double_t spdEta = 0.;

    const AliMultiplicity *spdMult=esdEv->GetMultiplicity();
    if(!spdMult){
        printf("Unable to get multiplicity! \n");
        return -1;
    }

    Int_t spdTracks = spdMult->GetNumberOfTracklets();  //SPD multiplicity

    for(Int_t iTrack = 0; iTrack< spdTracks; iTrack++){  //SPD track loop -> Directivity
	
		AliESDtrack* track = esdEv->GetTrack(iTrack);
		
		if (!track) {
				Printf("ERROR: Could not receive track %d", iTrack);
				continue;
		}
	
		if(!CheckIntCuts(2, track)) continue; 
		
		spdEta = spdMult->GetEta(iTrack);
		
		if((spdEta<0.0)||(spdEta>1.9)) continue;
	
		Double_t phi = spdMult->GetPhi(iTrack);  
	
		pxU = TMath::Cos(phi);
		pyU = TMath::Sin(phi);
	
		sumaPxU = sumaPxU + pxU;
		sumaPyU = sumaPyU + pyU;	
		
		goodtrack++;
	
    }//end SPD track loop

    if(goodtrack < 1.) return -1;

    dirU = (sqrt(pow(sumaPxU,2)+pow(sumaPyU,2)))/goodtrack;

    return dirU;
}

//__________________________________________________________________________

Int_t AliAnalysisCentralCutEvtESD::CalcSPDMult(AliESDEvent* const esdEv) {

	//Compute multiplicity - SPD ONLY

    InitCuts();

    const AliMultiplicity *spdMult=esdEv->GetMultiplicity();
    if(!spdMult){
        printf("Unable to get multiplicity! \n");
        return -1;
    }

    Int_t spdTracks = spdMult->GetNumberOfTracklets();  //SPD multiplicity

    Double_t spdEta;
    Int_t charged = 0;

    for(Int_t iTrack = 0; iTrack< spdTracks; iTrack++){  //second SPD track loop -> Multiplicity
	
		AliESDtrack* track = esdEv->GetTrack(iTrack);
		
		if (!track) {
				Printf("ERROR: Could not receive track %d", iTrack);
				continue;
		}
		
		if(!CheckIntCuts(2, track)) continue;
		
		spdEta = spdMult->GetEta(iTrack);
		
		if((spdEta<-0.5)||(spdEta>0.5)) continue;
	
		charged++;
	
    }//end second SPD track loop 

    return charged;
}
