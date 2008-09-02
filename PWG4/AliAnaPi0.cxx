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
/* $Id: $ */

//_________________________________________________________________________
// Class to collect two-photon invariant mass distributions for
// extractin raw pi0 yield.
//
//-- Author: Dmitri Peressounko (RRC "KI") 
//_________________________________________________________________________


// --- ROOT system ---
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TH3.h"

//---- AliRoot system ----
#include "AliAnaPi0.h"
#include "AliLog.h"
#include "AliCaloPhoton.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"

ClassImp(AliAnaPi0)
  
//____________________________________________________________________________
  AliAnaPi0::AliAnaPi0() : AliAnalysisTaskSE(),
  fNCentrBin(1),fNZvertBin(1),fNrpBin(1),fNPID(2),fNmaxMixEv(10),
  fCurCentrBin(0),fCurZvertBin(0),fCurRPBin(0),fPtMin(0.),fZvtxCut(40.),
  fEventsList(0x0),fCurrentEvent(0x0),fOutputList(0x0),fhEtalon(0x0),
  fhRe1(0x0),fhMi1(0x0),fhRe2(0x0),fhMi2(0x0),fhRe3(0x0),fhMi3(0x0),fhEvents(0x0) 
{
  //Default Ctor

}
//____________________________________________________________________________
  AliAnaPi0::AliAnaPi0(const char *name): AliAnalysisTaskSE(name),
  fNCentrBin(1),fNZvertBin(1),fNrpBin(1),fNPID(9),fNmaxMixEv(10),
  fCurCentrBin(0),fCurZvertBin(0),fCurRPBin(0),fPtMin(0.),fZvtxCut(40.),
  fEventsList(0x0),fCurrentEvent(0x0),fOutputList(0x0),fhEtalon(0x0),
  fhRe1(0x0),fhMi1(0x0),fhRe2(0x0),fhMi2(0x0),fhRe3(0x0),fhMi3(0x0),fhEvents(0x0) 
{
  //Ctor
  DefineOutput(1,TList::Class());
}

//____________________________________________________________________________
AliAnaPi0::AliAnaPi0(const AliAnaPi0 & ex) : AliAnalysisTaskSE(ex)  
{
  // cpy ctor
  //Do not need it
}
//_________________________________________________________________________
AliAnaPi0 & AliAnaPi0::operator = (const AliAnaPi0 & ex)
{
  // assignment operator

  if(this == &ex)return *this;
  ((AliAnalysisTaskSE *)this)->operator=(ex);
 
  return *this;

}
//____________________________________________________________________________
AliAnaPi0::~AliAnaPi0() {
  // Remove event containers
  if(fEventsList){
    for(Int_t ic=0; ic<fNCentrBin; ic++){
      for(Int_t iz=0; iz<fNZvertBin; iz++){
        for(Int_t irp=0; irp<fNrpBin; irp++){
          fEventsList[ic*fNZvertBin*fNrpBin+iz*fNrpBin+irp]->Delete() ;
          delete fEventsList[ic*fNZvertBin*fNrpBin+iz*fNrpBin+irp] ;
        }
      }
    }
    delete[] fEventsList; 
  }

  if(fhEtalon){
    delete fhEtalon ;
    fhEtalon = 0; 
  }
  if(fhRe1) delete[] fhRe1 ; //Do not delete histograms!
  if(fhRe2) delete[] fhRe2 ; //Do not delete histograms!
  if(fhRe3) delete[] fhRe3 ; //Do not delete histograms!
  if(fhMi1) delete[] fhMi1 ; //Do not delete histograms!
  if(fhMi2) delete[] fhMi2 ; //Do not delete histograms!
  if(fhMi3) delete[] fhMi3 ; //Do not delete histograms!

}
//________________________________________________________________________
void  AliAnaPi0::UserCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in fOutputContainer

  AliDebug(1,"Init inv. mass histograms");
  OpenFile(1) ; 

  fOutputList = new TList() ;
  fOutputList->SetName(GetName()) ;
 
  fhRe1=new TH3D*[fNCentrBin*fNPID] ;
  fhRe2=new TH3D*[fNCentrBin*fNPID] ;
  fhRe3=new TH3D*[fNCentrBin*fNPID] ;
  fhMi1=new TH3D*[fNCentrBin*fNPID] ;
  fhMi2=new TH3D*[fNCentrBin*fNPID] ;
  fhMi3=new TH3D*[fNCentrBin*fNPID] ;

  char key[255] ;
  char title[255] ;
  for(Int_t ic=0; ic<fNCentrBin; ic++){
    for(Int_t ipid=0; ipid<fNPID; ipid++){
      //Distance to bad module 1
      sprintf(key,"hRe_cen%d_pid%d_dist1",ic,ipid) ;
      sprintf(title,"Real m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      fhRe1[ic*fNPID+ipid]=(TH3D*)fhEtalon->Clone(key) ;
      fhRe1[ic*fNPID+ipid]->SetName(key) ;
      fhRe1[ic*fNPID+ipid]->SetTitle(title) ;
      fOutputList->Add(fhRe1[ic*fNPID+ipid]) ;

      sprintf(key,"hMi_cen%d_pid%d_dist1",ic,ipid) ;
      sprintf(title,"Mixed m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      fhMi1[ic*fNPID+ipid]=(TH3D*)fhEtalon->Clone(key) ;
      fhMi1[ic*fNPID+ipid]->SetName(key) ;
      fhMi1[ic*fNPID+ipid]->SetTitle(title) ;
      fOutputList->Add(fhMi1[ic*fNPID+ipid]) ;

      //Distance to bad module 2
      sprintf(key,"hRe_cen%d_pid%d_dist2",ic,ipid) ;
      sprintf(title,"Real m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      fhRe2[ic*fNPID+ipid]=(TH3D*)fhEtalon->Clone(key) ;
      fhRe2[ic*fNPID+ipid]->SetName(key) ;
      fhRe2[ic*fNPID+ipid]->SetTitle(title) ;
      fOutputList->Add(fhRe2[ic*fNPID+ipid]) ;

      sprintf(key,"hMi_cen%d_pid%d_dist2",ic,ipid) ;
      sprintf(title,"Mixed m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      fhMi2[ic*fNPID+ipid]=(TH3D*)fhEtalon->Clone(key) ;
      fhMi2[ic*fNPID+ipid]->SetName(key) ;
      fhMi2[ic*fNPID+ipid]->SetTitle(title) ;
      fOutputList->Add(fhMi2[ic*fNPID+ipid]) ;

      //Distance to bad module 3
      sprintf(key,"hRe_cen%d_pid%d_dist3",ic,ipid) ;
      sprintf(title,"Real m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      fhRe3[ic*fNPID+ipid]=(TH3D*)fhEtalon->Clone(key) ;
      fhRe3[ic*fNPID+ipid]->SetName(key) ; 
      fhRe3[ic*fNPID+ipid]->SetTitle(title) ;
      fOutputList->Add(fhRe3[ic*fNPID+ipid]) ;

      sprintf(key,"hMi_cen%d_pid%d_dist3",ic,ipid) ;
      sprintf(title,"Mixed m_{#gamma#gamma} distr. for centrality=%d and PID=%d",ic,ipid) ;
      fhMi3[ic*fNPID+ipid]=(TH3D*)fhEtalon->Clone(key) ;
      fhMi3[ic*fNPID+ipid]->SetName(key) ;
      fhMi3[ic*fNPID+ipid]->SetTitle(title) ;
      fOutputList->Add(fhMi3[ic*fNPID+ipid]) ;
    }
  }

  fhEvents=new TH3D("hEvents","Number of events",fNCentrBin,0.,1.*fNCentrBin,
                    fNZvertBin,0.,1.*fNZvertBin,fNrpBin,0.,1.*fNrpBin) ;
  fOutputList->Add(fhEvents) ;

  //Save parameters used for analysis
  fOutputList->Add(this);
 
}
//__________________________________________________
void AliAnaPi0::InitParameters()
{ 
  //Initialize the parameters of the analysis.
  
}
 //__________________________________________________
void AliAnaPi0::Init()
{ 
  //Make here all memory allocations
  //create etalon histo for all later histograms
  if(!fhEtalon){                                                   //  p_T      alpha   d m_gg    
    fhEtalon = new TH3D("hEtalon","Histo with binning parameters",20,0.,10.,10,0.,1.,200,0.,1.) ; 
    fhEtalon->SetXTitle("P_{T} (GeV)") ;
    fhEtalon->SetYTitle("#alpha") ;
    fhEtalon->SetZTitle("m_{#gamma#gamma} (GeV)") ;
  }
  
  //create event containers
  fEventsList = new TList*[fNCentrBin*fNZvertBin*fNrpBin] ;
  for(Int_t ic=0; ic<fNCentrBin; ic++){
    for(Int_t iz=0; iz<fNZvertBin; iz++){
      for(Int_t irp=0; irp<fNrpBin; irp++){
        fEventsList[ic*fNZvertBin*fNrpBin+iz*fNrpBin+irp] = new TList() ;
      }
    }
  }
}
//__________________________________________________________________
void AliAnaPi0::Print(const Option_t * /*opt*/) const
{
  //Print some relevant parameters set for the analysis
  printf("Class AliAnaPi0 for gamma-gamma inv.mass construction \n") ;
  printf("Number of bins in Centrality:  %d \n",fNCentrBin) ;
  printf("Number of bins in Z vert. pos: %d \n",fNZvertBin) ;
  printf("Number of bins in Reac. Plain: %d \n",fNrpBin) ;
  printf("Depth of event buffer: %d \n",fNmaxMixEv) ;
  printf("Number of different PID used:  %d \n",fNPID) ;
  printf("Cuts: \n") ;
  printf("Z vertex position: -%f < z < %f \n",fZvtxCut,fZvtxCut) ;
  printf("Minimal P_t: %f \n", fPtMin) ;
  printf("------------------------------------------------------\n") ;
  
} 
//__________________________________________________________________
void AliAnaPi0::UserExec(Option_t *)
{
  //Process one event and extract photons 
  //impose cuts on bad modules and bad runs
  //fill internal storage for subsequent histo filling

  AliVEvent* event = InputEvent();
 
  //Apply some cuts on event: vertex position and centrality range  
  Int_t iRun=event->GetRunNumber() ;
  if(IsBadRun(iRun))
    return ;

  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event) ;
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(event) ;
  if(esd){
     if(!FillFromESD(esd))
      return  ; //Event can not be used (vertex, centrality,... cuts not fulfilled)
   }
   else if(aod){
     if(!FillFromAOD(aod))
       return ; //Event can not be used (vertex, centrality,... cuts not fulfilled)
   }
   else{
     printf("Input neither ESD nor AOD, do nothing \n") ;
     return ;
   }


  fhEvents->Fill(fCurCentrBin+0.5,fCurZvertBin+0.5,fCurRPBin+0.5) ;

  if(fCurrentEvent->GetEntriesFast()>0){
    //Reduce size for storing
    fCurrentEvent->Expand(fCurrentEvent->GetEntriesFast()) ;
    FillHistograms() ;
  }

  PostData(1, fOutputList);
 
} 
//__________________________________________________________________
Bool_t AliAnaPi0::FillFromESD(AliESDEvent * esd){
  //Fill photon list from ESD applying 
  //some cut should be applyed:
  const Float_t kMinDist=2. ; //Minimal distance to bad channel to accept cluster
  const Float_t kMinDist2=4.; //Cuts on Minimal distance 
  const Float_t kMinDist3=5.; //used for acceptance-efficiency study
  const Float_t kDispCut=1.5; //Cut on dispersion, used in PID evaluation
  const Float_t kTOFCut=5.e-9;//Cut on TOF, used in PID evaluation 
  const Float_t kPhotPID=0.6; //Baesian PID for photon

  //Impose cut on vertex and calculate Zvertex bin
  const AliESDVertex *esdV = esd->GetVertex() ;
  esdV->GetXYZ(fVert);
  if(fVert[2]<-fZvtxCut || fVert[2]> fZvtxCut)
    return kFALSE ;
  fCurZvertBin=(Int_t)(0.5*fNZvertBin*(fVert[2]+fZvtxCut)/fZvtxCut) ;
    
  //Get Centrality and calculate centrality bin
  //Does not exist in ESD yet???????
  fCurCentrBin=0 ;
 
  //Get Reaction Plain position and calculate RP bin
  //does not exist in ESD yet????
  fCurRPBin=0 ;
 
  //************************  PHOS *************************************
  TRefArray * caloClustersArr  = new TRefArray();
  esd->GetPHOSClusters(caloClustersArr);
 
  const Int_t kNumberOfPhosClusters   = caloClustersArr->GetEntries() ;
  //if fCurrentEvent!=0 it was not filled in previous event, still empty and can be used
  if(!fCurrentEvent) 
    fCurrentEvent = new TClonesArray("AliCaloPhoton",kNumberOfPhosClusters) ;
  Int_t inList=0 ;
 
  // loop over the PHOS Cluster
  for(Int_t i = 0 ; i < kNumberOfPhosClusters ; i++) {
    AliESDCaloCluster * calo = (AliESDCaloCluster *) caloClustersArr->At(i) ;
 
    //Make some tests
    if(calo->E()<fPtMin)
      continue ;

    Double_t distBad=calo->GetDistanceToBadChannel() ; //Distance to bad in cm
    if(distBad<0.)distBad=9999. ; //workout strange convension dist = -1. ;
    if(distBad<kMinDist) //In bad channel (cristall size 2.2x2.2 cm)
      continue ;

    new((*fCurrentEvent)[inList])AliCaloPhoton() ;
    AliCaloPhoton * ph = static_cast<AliCaloPhoton*>(fCurrentEvent->At(inList)) ;
    inList++ ;

    //Set Momentum 
    calo->GetMomentum(*ph,fVert);

    //Dispersion/lambdas
    Double_t disp=calo->GetClusterDisp()  ;
//    Double_t m20=calo->GetM20() ;
//    Double_t m02=calo->GetM02() ; 
    ph->SetDispBit(disp<kDispCut) ;  

    //TOF
    Double_t tof=calo->GetTOF()  ;
    ph->SetTOFBit(TMath::Abs(tof)<kTOFCut) ; 
 
    //Charged veto
//    Double_t cpvR=calo->GetEmcCpvDistance() ; 
    Int_t itr=calo->GetNTracksMatched();  //number of track matched
    ph->SetCPVBit(itr>0) ;  //Temporary cut, should we evaluate distance?

    //Overall PID
    Double_t *pid=calo->GetPid();
    ph->SetPCAPID(pid[AliPID::kPhoton]>kPhotPID) ;
      
    //Set Distance to Bad channel
    if(distBad>kMinDist3)
      ph->SetDistToBad(2) ;
    else if(distBad>kMinDist2)
        ph->SetDistToBad(1) ; 
    else
       ph->SetDistToBad(0) ;

  }//loop
  return kTRUE ;
}
//__________________________________________________________________
Bool_t AliAnaPi0::FillFromAOD(AliAODEvent * aod){
  //Fill photon list from AOD applying 
  //some cuts:
  const Float_t kMinDist=2. ; //Minimal distance to bad channel to accept cluster
  const Float_t kMinDist2=4.; //Cuts on Minimal distance 
  const Float_t kMinDist3=5.; //used for acceptance-efficiency study
  const Float_t kDispCut=1.5; //Cut on dispersion, used in PID evaluation
  const Float_t kTOFCut=5.e-9;//Cut on TOF, used in PID evaluation 
  const Float_t kPhotPID=0.6 ; //Baesian PID for photon
 
  //Impose cut on vertex and calculate Zvertex bin
  const AliAODVertex *aodV = aod->GetPrimaryVertex() ;
  fVert[0]=aodV->GetX();
  fVert[1]=aodV->GetY();
  fVert[2]=aodV->GetZ();
  if(fVert[2]<-fZvtxCut || fVert[2]> fZvtxCut)
    return kFALSE ;
  fCurZvertBin=(Int_t)(0.5*fNZvertBin*(fVert[2]+fZvtxCut)/fZvtxCut) ;
 
  //Get Centrality and calculate centrality bin
  //Does not exist in ESD yet???????
  fCurCentrBin=0 ;
 
  //Get Reaction Plain position and calculate RP bin
  //does not exist in ESD yet????
  fCurRPBin=0 ;
 
  //************************  PHOS *************************************
  const Int_t kNumberOfPhosClusters   = aod->GetNCaloClusters() ;
  if(!fCurrentEvent)
    fCurrentEvent = new TClonesArray("AliCaloPhoton",kNumberOfPhosClusters) ;
  Int_t inList=0 ;
 
  // loop over the PHOS Cluster
  for(Int_t i = 0 ; i < kNumberOfPhosClusters ; i++) {
    AliAODCaloCluster * calo = aod->GetCaloCluster(i);
      
    //Make some tests
    if(calo->E()<fPtMin)
      continue ;

    Double_t distBad=calo->GetDistToBadChannel() ;
    if(distBad<kMinDist) //In bad channel
      continue ;

    new((*fCurrentEvent)[inList])AliCaloPhoton() ;
    AliCaloPhoton * ph = static_cast<AliCaloPhoton*>(fCurrentEvent->At(inList)) ;
    inList++ ;

    //Set Momentum 
    calo->GetMomentum(*ph,fVert);

    //Dispersion/lambdas
    Double_t disp=calo->GetDispersion()  ;
//    Double_t m20=calo->GetM20() ;
//    Double_t m02=calo->GetM02() ; 
    ph->SetDispBit(disp<kDispCut) ; 

    //TOF
    Double_t tof=calo->GetTOF()  ;
    ph->SetTOFBit(TMath::Abs(tof)<kTOFCut) ;
 
    //Charged veto
//    Double_t cpvR=calo->GetEmcCpvDistance() ; 
    Int_t itr=calo->GetNTracksMatched();  //number of track matched
    ph->SetCPVBit(itr>0) ;  //Temporary cut !!!!

    //Overall PID
    Double_t pid[13];
    calo->GetPID(pid);
    ph->SetPCAPID(pid[AliAODCluster::kPhoton]>kPhotPID) ;
      
    //Distance to Bad
    if(distBad>kMinDist3)
      ph->SetDistToBad(2) ; 
    else if(distBad>kMinDist2)
        ph->SetDistToBad(1) ;
    else
       ph->SetDistToBad(0) ;
 

  }//loop
  return kTRUE ;
}
//__________________________________________________________________
void  AliAnaPi0::FillHistograms() 
{
  //Fill Real and Mixed invariant mass histograms
  //then add current event to the buffer and 
  //remove redundant events from buffer if necessary
  
  //Fill histograms only if list of particles for current event was created
  if(fCurrentEvent==0)
    return ;

  //Fill Real distribution
  Int_t nPhot = fCurrentEvent->GetEntriesFast() ;
  for(Int_t i1=0; i1<nPhot-1; i1++){
    AliCaloPhoton * p1 = static_cast<AliCaloPhoton*>(fCurrentEvent->At(i1)) ;
    for(Int_t i2=i1+1; i2<nPhot; i2++){
      AliCaloPhoton * p2 = static_cast<AliCaloPhoton*>(fCurrentEvent->At(i2)) ;
      Double_t m =(*p1 + *p2).M() ;
      Double_t pt=(*p1 + *p2).Pt();
      Double_t a = TMath::Abs(p1->Energy()-p2->Energy())/(p1->Energy()+p2->Energy()) ;
      for(Int_t ipid=0; ipid<fNPID; ipid++){
        if(p1->IsPIDOK(ipid)&&p2->IsPIDOK(ipid)){   
          fhRe1[fCurCentrBin*fNPID+ipid]->Fill(pt,a,m) ;
          if(p1->DistToBad()>0 && p2->DistToBad()>0){
            fhRe2[fCurCentrBin*fNPID+ipid]->Fill(pt,a,m) ;
            if(p1->DistToBad()>1 && p2->DistToBad()>1){
              fhRe3[fCurCentrBin*fNPID+ipid]->Fill(pt,a,m) ;
            }
          }
        }
      } 
    }
  }

  //Fill mixed
  TList * evMixList=fEventsList[fCurCentrBin*fNZvertBin*fNrpBin+fCurZvertBin*fNrpBin+fCurRPBin] ;
  Int_t nMixed = evMixList->GetSize() ;
  for(Int_t ii=0; ii<nMixed; ii++){  
    TClonesArray* ev2=dynamic_cast<TClonesArray*>(evMixList->At(ii));
    Int_t nPhot2=ev2->GetEntriesFast() ;
    for(Int_t i1=0; i1<nPhot; i1++){
      AliCaloPhoton * p1 = static_cast<AliCaloPhoton*>(fCurrentEvent->At(i1)) ;
      for(Int_t i2=0; i2<nPhot2; i2++){
        AliCaloPhoton * p2 = static_cast<AliCaloPhoton*>(ev2->At(i2)) ;
        Double_t m =(*p1 + *p2).M() ;
        Double_t pt=(*p1 + *p2).Pt();
        Double_t a = TMath::Abs(p1->Energy()-p2->Energy())/(p1->Energy()+p2->Energy()) ;
        for(Int_t ipid=0; ipid<fNPID; ipid++){
          if(p1->IsPIDOK(ipid)&&p2->IsPIDOK(ipid)){
            fhMi1[fCurCentrBin*fNPID+ipid]->Fill(pt,a,m) ;
            if(p1->DistToBad()>0 && p2->DistToBad()>0){
              fhMi2[fCurCentrBin*fNPID+ipid]->Fill(pt,a,m) ;
              if(p1->DistToBad()>1 && p2->DistToBad()>1){
                fhMi3[fCurCentrBin*fNPID+ipid]->Fill(pt,a,m) ;
              }
            }
          }
        }
      }
    }
  }

  //Add current event to buffer and Remove redandant events 
  if(fCurrentEvent->GetEntriesFast()>0){
    evMixList->AddFirst(fCurrentEvent) ;
    fCurrentEvent=0 ; //Now list of particles belongs to buffer and it will be deleted with buffer
    if(evMixList->GetSize()>=fNmaxMixEv){
      TClonesArray * tmp = dynamic_cast<TClonesArray*>(evMixList->Last()) ;
      evMixList->RemoveLast() ;
      delete tmp ;
    }
  } 
  else{ //empty event
   delete fCurrentEvent ;
   fCurrentEvent=0 ; 
  }
}
//______________________________________________________________________________
void AliAnaPi0::Terminate(Option_t *)
{
}
