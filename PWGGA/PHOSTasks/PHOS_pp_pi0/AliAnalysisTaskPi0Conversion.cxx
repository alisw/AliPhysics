#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "TFile.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "THashList.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskPi0Conversion.h"
#include "AliCaloPhoton.h"
#include "AliAODMCParticle.h"
#include "AliPHOSGeometry.h"
#include "AliAODEvent.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAODVertex.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliPHOSCalibData.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBMetaData.h"
//#include "AliPHOSAodCluster1.h"
#include "AliOADBContainer.h"


// Analysis task to fill histograms with PHOS AOD clusters and cells
// Authors: Yulia Demkina
// Date   : 24.03.2016

ClassImp(AliAnalysisTaskPi0Conversion)
//________________________________________________________________________
AliAnalysisTaskPi0Conversion::AliAnalysisTaskPi0Conversion(const char *name) 
: AliAnalysisTaskSE(name),
  fOutputContainer(0),
  fPHOSEvents(0x0),
  fStack(0x0),
  fPHOSEvent(0),
  fPHOSGeo(0)
{
  // Constructor
  
  // Output slots #0 write into a TH1 container
  DefineOutput(1,TList::Class());

}

//________________________________________________________________________
void AliAnalysisTaskPi0Conversion::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  // AOD histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
  }
  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);

  fOutputContainer->Add(new TH1I("hCellMultEvent"  ,"PHOS cell multiplicity per event"    ,2000,0,2000));
  fOutputContainer->Add(new TH1I("hClusterMult"      ,"CaloCluster multiplicity"     ,100,0,100));
  fOutputContainer->Add(new TH1F("hSelEvents","Selected events",5,0.,5.));
  
  for(Int_t mod=1; mod<4; mod++){
    fOutputContainer->Add(new TH2F(Form("hCluLowM%d",mod),Form("Cell (X,Z), M%d",mod) ,64,0.5,64.5, 56,0.5,56.5));

   fOutputContainer->Add(new TH3F(Form("hMggAllM%d",mod),Form("M%d",mod) ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
   fOutputContainer->Add(new TH3F(Form("hMggDispM%d",mod),Form("M%d",mod) ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
    fOutputContainer->Add(new TH3F(Form("hMiMggAllM%d",mod),Form("M%d",mod) ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
    fOutputContainer->Add(new TH3F(Form("hMiMggDispM%d",mod),Form("M%d",mod) ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
    fOutputContainer->Add(new TH2F(Form("hMggPtM%d",mod),Form("M%d",mod),400,0.,1.,40,0.,20.));
    fOutputContainer->Add(new TH2F(Form("hMggPtM%dAll",mod),Form("M%d",mod),400,0.,1.,40,0.,20.));
    fOutputContainer->Add(new TH2F(Form("hMggPtM%dCPV",mod),Form("M%d",mod),400,0.,1.,40,0.,20.));
    fOutputContainer->Add(new TH2F(Form("hMggPtM%dDisp",mod),Form("M%d",mod),400,0.,1.,40,0.,20.));
    fOutputContainer->Add(new TH2F(Form("hMggPtM%dBoth",mod),Form("M%d",mod),400,0.,1.,40,0.,20.));
    fOutputContainer->Add(new TH2F(Form("hMiMggPtM%d",mod),Form("M%d",mod),400,0.,1.,40,0.,20.));
    fOutputContainer->Add(new TH2F(Form("hMiMggPtM%dAll",mod),Form("M%d",mod),400,0.,1.,40,0.,20.));
    fOutputContainer->Add(new TH2F(Form("hMiMggPtM%dCPV",mod),Form("M%d",mod),400,0.,1.,40,0.,20.));
    fOutputContainer->Add(new TH2F(Form("hMiMggPtM%dDisp",mod),Form("M%d",mod),400,0.,1.,40,0.,20.));
    fOutputContainer->Add(new TH2F(Form("hMiMggPtM%dBoth",mod),Form("M%d",mod),400,0.,1.,40,0.,20.));
  }

 //MY CODE 
  fOutputContainer->Add(new TH1F("hVertex_pi0","R of primary particles",150,0.,150.));
  fOutputContainer->Add(new TH1F("hVertex_gamma","R of primary particles",150,0.,150.));
  fOutputContainer->Add(new TH1F("hVertex_el","R of primary particles",150,0.,150.));
  fOutputContainer->Add(new TH1F("hVertex_pos","R of primary particles",150,0.,150.));
//  fOutputContainer->Add(new TH1F("hClusterPrimary","N of primary particles in clusters",20,0.,20.)); 
//  fOutputContainer->Add(new TH1F("hInvMass","Inv mass two photons",500,0.,0.5));
//  fOutputContainer->Add(new TH1F("hInvMassNew","Inv mass two photons with P0",500,0.,0.5));
  fOutputContainer->Add(new TH2F("hInvMassP0","Inv mass two photons with P0",500,0.,0.5,100,0.,10.));   
  fOutputContainer->Add(new TH2F("hInvMassGamma","Inv mass two photons with gamma",500,0.,0.5,100,0.,10.));
  fOutputContainer->Add(new TH2F("hInvMassEl","Inv mass two photons with Electron",500,0.,0.5,100,0.,10.));
  fOutputContainer->Add(new TH2F("hInvMassElse","Inv mass two photons Default",500,0.,0.5,100,0.,10.));
  fOutputContainer->Add(new TH2F("hInvMassEt","Inv mass two photons Eta",500,0.,0.5,100,0.,10.));
  fOutputContainer->Add(new TH2F("hInvMassK0s","Inv mass two photons K0 Short",500,0.,0.5,100,0.,10.));
  fOutputContainer->Add(new TH2F("hInvMassK0l","Inv mass two photons K0 Long",500,0.,0.5,100,0.,10.));
  fOutputContainer->Add(new TH2F("hInvMassJet","Inv mass two photons Jet",500,0.,0.5,100,0.,10.));


  PostData(1, fOutputContainer);

}

//________________________________________________________________________
void AliAnalysisTaskPi0Conversion::UserExec(Option_t *) 
{
  // Main loop, called for each event
  // Analyze AOD
  const Double_t logWeight=4.5 ;  
  const Double_t rcut=1. ;
  
  FillHistogram("hSelEvents",0.5) ;
  AliAODEvent *event = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!event) {
     Printf("ERROR: Could not retrieve event");
     PostData(1, fOutputContainer);
     return;
  }
  
 // Initialize the PHOS geometry
  if(!fPHOSGeo){
    //should be initialized by Tender
    fPHOSGeo = AliPHOSGeometry::GetInstance() ; 
    if(!fPHOSGeo){ //if not yet created, create manually   
      fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;
      AliOADBContainer geomContainer("phosGeo");
      geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
      TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(233715,"PHOSRotationMatrixes");
      for(Int_t mod=0; mod<6; mod++) {
        if(!matrixes->At(mod)) continue;
        fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;
        printf(".........Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo) ;
        ((TGeoHMatrix*)matrixes->At(mod))->Print() ;
      } 
    }
  }
  
  FillHistogram("hSelEvents",1.5) ;

  if(fPHOSEvent){
    fPHOSEvent->Clear() ;
  }
  else{
    fPHOSEvent = new TClonesArray("AliCaloPhoton",50) ;
  }

  // Checks if we have a primary vertex
  // Get primary vertices form AOD

  Double_t vtx0[3] = {0,0,0}; // don't rely on AOD vertex, assume (0,0,0) //WHY IT IS SO????
  Double_t vtx5[3] = {0.,0.,0.};
  TVector3 vertex(vtx0);

  FillHistogram("hSelEvents",2.5) ; //WHAT MEANS THE SECOND NUMBER?

  //Vtx class z-bin
  Int_t zvtx = (Int_t)((vtx5[2]+10.)/2.) ;
  if(zvtx<0)zvtx=0 ;
  if(zvtx>9)zvtx=9 ;


  Int_t centr=0 ;
  //always zero centrality
  if(!fPHOSEvents) fPHOSEvents=new TList() ;

  AliAODCaloCluster *clu1;
  TLorentzVector p1,p2,p12, pv1,pv2,pv12;
  AliAODCaloCells *cells      = event->GetPHOSCells();

  Int_t multClust = event->GetNumberOfCaloClusters();
  Int_t multCells      = cells     ->GetNumberOfCells();
  FillHistogram("hClusterMult",multClust);
  FillHistogram("hCellMultEvent",multCells);

  Float_t  energy;
  Int_t    mod1, relId[4], cellAbsId, cellX, cellZ;

  // Single loop over cells

  Int_t nCellModule[4] = {0,0,0,0};

  // Single loop over clusters fills cluster histograms

  Int_t    digMult;
  Int_t    multPHOSClust[4]  = {0,0,0,0};
  Float_t  position[3];

  fStack = (TClonesArray*)event->FindListObject(AliAODMCParticle::StdBranchName());//DEFINE fStack

  //Select photons for inv mass calculation
  Int_t inPHOS=0 ;
  for (Int_t i1=0; i1<multClust; i1++) {
    clu1 = event->GetCaloCluster(i1);
    if ( !clu1->IsPHOS() || clu1->E()<0.3) 
      continue;
    if(clu1->GetNCells()<3)
      continue ;          
    if(clu1->GetM02()<0.2) 
      continue ;       

    clu1->GetPosition(position);
    TVector3 global1(position) ;
    fPHOSGeo->GlobalPos2RelId(global1,relId) ;
    mod1  = relId[0] ;
    cellX = relId[2];
    cellZ = relId[3] ;
  
    UInt_t inCLU = clu1->GetNLabels();
//    FillHistogram("hClusterPrimary",inCLU);
//

    if(clu1->E()>0.5)
      FillHistogram(Form("hCluLowM%d",mod1),cellX,cellZ,1.);

    clu1->GetMomentum(p1 ,vtx0); 
    clu1->GetMomentum(pv1,vtx5);
    new((*fPHOSEvent)[inPHOS]) AliCaloPhoton(p1.X(),p1.Y (),p1.Z(),p1.E()) ;
    AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSEvent->At(inPHOS) ; 
    ph->SetModule(mod1) ;
    ph->SetMomV2(&pv1) ;
    ph->SetDispBit(clu1->Chi2()<2.5*2.5) ;
    ph->SetCPVBit(clu1->GetEmcCpvDistance()>2.5) ;
    ph->SetEMCx(float(cellX)) ;
    ph->SetEMCz(float(cellZ)) ;

    Int_t primLabel=FindAODLabel(clu1->GetLabelAt(0)) ; 
    if(primLabel>-1){
           AliAODMCParticle * prim = (AliAODMCParticle*)fStack->At(primLabel) ;
           Int_t iparent=primLabel;
           AliAODMCParticle * parent = prim;
           Double_t r2=prim->Xv()*prim->Xv()+prim->Yv()*prim->Yv() ;
           while((r2 > rcut*rcut) && (iparent>-1)){
               iparent=parent->GetMother();
               parent=(AliAODMCParticle*)fStack->At(iparent);
               r2=parent->Xv()*parent->Xv()+parent->Yv()*parent->Yv() ;
           }
          ph->SetPrimary(primLabel) ;
          ph->SetPrimaryAtVertex(iparent) ;

     }
     else{
         ph->SetPrimary(-1); //Primary index    
         ph->SetPrimaryAtVertex(-1) ;
     }
  
    inPHOS++ ;
  } 


  //MY CODE

/*
  for (Int_t i=0; i<inPHOS; i++){ //ALL PHOTONS 
    AliCaloPhoton *p = static_cast<AliCaloPhoton*>(fPHOSEvent->At(i));    
    Int_t label=p->GetPrimary() ;//Index of Primary particle
    if (label < 0) 
      continue;
    AliAODMCParticle * prim = (AliAODMCParticle*)fStack->At(p->GetPrimary()) ; //Ykazatel on Primary Particle in fStack 
    Double_t r2=TMath::Sqrt (prim->Xv()*prim->Xv()+prim->Yv()*prim->Yv());
    FillHistogram(Form("hVertex_%s", partName), r2);
  }
*/
/*
  for (Int_t i=0; i<inPHOS; i++){  
    AliCaloPhoton *p1 = static_cast<AliCaloPhoton*>(fPHOSEvent->At(i));    
    for (Int_t j=0; j<i; j++){ 

      AliCaloPhoton *p2 = static_cast<AliCaloPhoton*>(fPHOSEvent->At(j));   
      Double_t m = (*p1 + *p2).M();
      FillHistogram("hInvMass",m);
   }
  }
*/
 
  if(fStack){ //only for MC case
  for (Int_t i=0; i<inPHOS; i++){  
    AliCaloPhoton *p1 = static_cast<AliCaloPhoton*>(fPHOSEvent->At(i)); 
    Int_t klabel = p1->GetPrimary() ;
    if (klabel < 0)
      continue;
    AliAODMCParticle * kprim = (AliAODMCParticle*)fStack->At(p1->GetPrimary()) ;
    Int_t kpdg = kprim -> GetPdgCode();

    while ((kpdg != 111) && (klabel>-1)) {
      Int_t klabel = kprim -> GetMother();
      if (klabel < 0) break; 
      kprim = (AliAODMCParticle*)fStack->At(klabel) ;
      kpdg = kprim -> GetPdgCode();
    } //while prim for p1

    if (kpdg == 111) {
      for (Int_t j=0; j<i; j++){ 
        AliCaloPhoton *p2 = static_cast<AliCaloPhoton*>(fPHOSEvent->At(j));  
        Int_t mlabel = p2->GetPrimary() ;
        if (mlabel < 0)
	  continue;
        AliAODMCParticle * mprim = (AliAODMCParticle*)fStack->At(p2->GetPrimary()) ;
        Int_t mpdg = mprim -> GetPdgCode();
        while ((mpdg != 111) && (mlabel>-1)) {
          Int_t mlabel = mprim -> GetMother();
          if (mlabel < 0) break; 
          mprim = (AliAODMCParticle*)fStack->At(mlabel) ;
          mpdg = mprim -> GetPdgCode(); 
        }  // while prim for p2      
        if (mpdg == 111){
          if(mprim == kprim){
            Double_t m = (*p1 + *p2).M();
//               FillHistogram("hInvMassNew",m); 
          }
        }   
      } //for j
    } //for pdg=111 for p1
  } //for i
  }

// MY CODE 4
 if(fStack){
 for (Int_t i=0; i<inPHOS; i++){  
   AliCaloPhoton *p1 = static_cast<AliCaloPhoton*>(fPHOSEvent->At(i)); 
   Int_t klabel1 = p1->GetPrimary();
   if (klabel1 < 0) 
     continue;
   AliAODMCParticle * kprim1 = (AliAODMCParticle*)fStack->At(p1->GetPrimary());
   AliAODMCParticle  * kprim = kprim1;
 
   for (Int_t j=0; j<i; j++){ 
     AliCaloPhoton *p2 = static_cast<AliCaloPhoton*>(fPHOSEvent->At(j));  
     Int_t mlabel1 = p2->GetPrimary();
     if (mlabel1 < 0)
       continue;
     AliAODMCParticle * mprim1 = (AliAODMCParticle*)fStack->At(p2->GetPrimary());
     AliAODMCParticle  * mprim = mprim1;

     Int_t klabel = klabel1;
     while (klabel>-1) {
        Int_t klabel = kprim -> GetMother();
        if (klabel < 0)
	  break; 
        kprim = (AliAODMCParticle*)fStack->At(klabel) ;
        Int_t mlabel = mlabel1;       
        while (mlabel>-1) {
            Int_t mlabel = mprim -> GetMother();
            if (mlabel < 0) break; 
            mprim = (AliAODMCParticle*)fStack->At(mlabel) ; 
              if (klabel == mlabel) {
                Double_t m = (*p1 + *p2).M(); 
                Double_t pt = (*p1 + *p2).Pt();
                Int_t pdg = mprim -> GetPdgCode(); 
                switch (pdg) {
                case 111: FillHistogram("hInvMassP0",m,pt); break; 
                case 22: FillHistogram("hInvMassGamma",m,pt); break;
                case 11:
                case -11: 
{ 
                    FillHistogram("hInvMassEl",m,pt);
                    Int_t wpdg = pdg;
                    Int_t wlabel = mlabel;
                    AliAODMCParticle  * wprim = mprim;
                    while ((wpdg != 111) && (wlabel>-1)) {
        
                     Int_t wlabel = wprim -> GetMother();
                      if (wlabel < 0) break; 
                      wprim = (AliAODMCParticle*)fStack->At(wlabel) ;
                      wpdg = wprim -> GetPdgCode(); }      
                    if (wpdg == 111) {
                        Double_t r3=TMath::Sqrt (mprim->Xv()*mprim->Xv()+mprim->Yv()*mprim->Yv());    
                        FillHistogram("hConversion",m,pt,r3); 
                        };
}; break;
                case 221: FillHistogram("hInvMassEt",m,pt); break;
                case 130: FillHistogram("hInvMassK0s",m,pt); break;
                case 310: FillHistogram("hInvMassK0l",m,pt); break;
                default:  if  ((abs(pdg) < 10) || (pdg == 21)) FillHistogram("hInvMassJet",m, pt);
                          else FillHistogram("hInvMassElse",m, pt);

                } //switch pdg
                mlabel = -1;
                klabel = -1;
              } // klabel=mlabel
        }  // while prim for p2      
     } //while prim for p1
   } //for j
 } //for i
 }
 // Fill Real disribution

  char key[55];
  for (Int_t i1=0; i1<inPHOS-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    for (Int_t i2=i1+1; i2<inPHOS; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent->At(i2) ;
      p12  = *ph1  + *ph2;
      if(ph1->Module()==ph2->Module())
	FillHistogram(Form("hMggPtM%d",ph1->Module()),p12.M(),p12.Pt()) ;

      FillHistogram(Form("hMggPtM%dAll",ph1->Module()),p12.M(),p12.Pt()) ;
      if(ph1->IsDispOK() && ph2->IsDispOK())
        FillHistogram(Form("hMggPtM%dDisp",ph1->Module()),p12.M(),p12.Pt()) ;
      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
        FillHistogram(Form("hMggPtM%dCPV",ph1->Module()),p12.M(),p12.Pt()) ;
        if(ph1->IsDispOK() && ph2->IsDispOK())
          FillHistogram(Form("hMggPtM%dBoth",ph1->Module()),p12.M(),p12.Pt()) ;
      }
      
      if(p12.Pt()<1.5)
	continue ;
      FillHistogram(Form("hMggAllM%d",ph1->Module()),ph1->EMCx(),ph1->EMCz(),p12.M() );
      FillHistogram(Form("hMggAllM%d",ph2->Module()),ph2->EMCx(),ph2->EMCz(),p12.M() );

      if(ph1->IsDispOK())
	FillHistogram(Form("hMggDispM%d",ph1->Module()),ph1->EMCx(),ph1->EMCz(),p12.M() );
      if(ph2->IsDispOK())
	FillHistogram(Form("hMggDispM%d",ph2->Module()),ph2->EMCx(),ph2->EMCz(),p12.M() );
    } // end of loop i2
  } // end of loop i1

  
  //now mixed
  for (Int_t i1=0; i1<inPHOS; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    for(Int_t ev=0; ev<fPHOSEvents->GetSize();ev++){
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(fPHOSEvents->At(ev)) ;
      for(Int_t i2=0; i2<mixPHOS->GetEntriesFast();i2++){
	AliCaloPhoton * ph2=(AliCaloPhoton*)mixPHOS->At(i2) ;
	p12  = *ph1  + *ph2;
        if(ph1->Module()==ph2->Module())
	  FillHistogram(Form("hMiMggPtM%d",ph1->Module()),p12.M(),p12.Pt()) ;
 
	FillHistogram(Form("hMiMggPtM%dAll",ph1->Module()),p12.M(),p12.Pt()) ;
        if(ph1->IsDispOK() && ph2->IsDispOK())
          FillHistogram(Form("hMiMggPtM%dDisp",ph1->Module()),p12.M(),p12.Pt()) ;
        if(ph1->IsCPVOK() && ph2->IsCPVOK()){
          FillHistogram(Form("hMiMggPtM%dCPV",ph1->Module()),p12.M(),p12.Pt()) ;
          if(ph1->IsDispOK() && ph2->IsDispOK())
            FillHistogram(Form("hMiMggPtM%dBoth",ph1->Module()),p12.M(),p12.Pt()) ;
        }
 
	
        if(p12.Pt()<1.5)
  	  continue ;
	
	FillHistogram(Form("hMiMggAllM%d",ph1->Module()),ph1->EMCx(),ph1->EMCz(),p12.M() );
	FillHistogram(Form("hMiMggAllM%d",ph2->Module()),ph2->EMCx(),ph2->EMCz(),p12.M() );
	
	if(ph1->IsDispOK())
	  FillHistogram(Form("hMiMggDispM%d",ph1->Module()),ph1->EMCx(),ph1->EMCz(),p12.M() );
	if(ph2->IsDispOK())
	  FillHistogram(Form("hMiMggDispM%d",ph2->Module()),ph2->EMCx(),ph2->EMCz(),p12.M() );	  

      } // end of loop i2
    }
  } // end of loop i1
  
  
  //Now we either add current events to stack or remove
  //If no photons in current event - no need to add it to mixed
  if(fPHOSEvent->GetEntriesFast()>0){
    fPHOSEvents->AddFirst(fPHOSEvent) ;
    fPHOSEvent=0;
    if(fPHOSEvents->GetSize()>100){//Remove redundant events
      TClonesArray * tmp = static_cast<TClonesArray*>(fPHOSEvents->Last()) ;
      fPHOSEvents->RemoveLast() ;
      delete tmp ;
    }
  }
  // Post output data.
  PostData(1, fOutputContainer);
}


//_____________________________________________________________________________
void AliAnalysisTaskPi0Conversion::FillHistogram(const char * key,Double_t x)const{
  //FillHistogram
  TH1I * tmpI = dynamic_cast<TH1I*>(fOutputContainer->FindObject(key)) ;
  if(tmpI){
    tmpI->Fill(x) ;
    return ;
  }
  TH1F * tmpF = dynamic_cast<TH1F*>(fOutputContainer->FindObject(key)) ;
  if(tmpF){
    tmpF->Fill(x) ;
    return ;
  }
  TH1D * tmpD = dynamic_cast<TH1D*>(fOutputContainer->FindObject(key)) ;
  if(tmpD){
    tmpD->Fill(x) ;
    return ;
  }
  AliInfo(Form("can not find histogram <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskPi0Conversion::FillHistogram(const char * key,Double_t x,Double_t y)const{
  //FillHistogram
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH1F")){
    ((TH1F*)tmp)->Fill(x,y) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")){
    ((TH2F*)tmp)->Fill(x,y) ;
    return ;
  }
  AliError(Form("Calling FillHistogram with 2 parameters for histo <%s> of type %s",key,tmp->IsA()->GetName())) ;
}

//_____________________________________________________________________________
void AliAnalysisTaskPi0Conversion::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
  //Fills 1D histograms with key
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")){
    ((TH2F*)tmp)->Fill(x,y,z) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH3F")){
    ((TH3F*)tmp)->Fill(x,y,z) ;
    return ;
  }
}
//___________________________________________________________________________
Int_t AliAnalysisTaskPi0Conversion::FindAODLabel(Int_t esdLabel)const{
   
  if(esdLabel<0)
     return -1 ;
  
  Int_t n = fStack->GetEntriesFast();
  if(esdLabel<n){
      AliAODMCParticle* tmp =  (AliAODMCParticle*)fStack->At(esdLabel) ;
      if(tmp->Label()==esdLabel)
	return esdLabel;
      else{
	Int_t i=esdLabel;
	if(esdLabel>tmp->Label()){
	  i++;
	  while(i<n){
	   tmp =  (AliAODMCParticle*)fStack->At(i) ;
	   if(tmp->Label()==esdLabel)
	     return i;
	   i++;
	  }
	  return -1 ; //failed to find 
	}
	else{
	  i--;
	  while(i>=0){
	   tmp =  (AliAODMCParticle*)fStack->At(i) ;
	   if(tmp->Label()==esdLabel)
	     return i;
	   i--;
	  }
	  return -1 ; //failed to find 
	
	}
      }
      
    }
    else{
       Int_t i=n-1;
       while(i>=0){
	   AliAODMCParticle* tmp =  (AliAODMCParticle*)fStack->At(i) ;
	   if(tmp->Label()==esdLabel)
	     return i;
	   i--;
       }
       return -1 ; //failed to find 
    }
    return -1 ;
  
}
