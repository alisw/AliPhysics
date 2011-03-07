/*
 The class for calculating the global (not detector specific) quality assurance.
 It reuses the following TLists from its base class 
    AliQADataMaker::fRecPointsQAList (for keeping the track residuals)
    AliQADataMaker::fESDsQAList      (for keeping global ESD QA data)
*/

#include <TPDGCode.h>
#include <TH1F.h>
#include <TH2F.h>

#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliDetectorRecoParam.h"
#include "AliQAChecker.h"
#include "AliGlobalQADataMaker.h"
#include "AliGeomManager.h"
#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliRawReader.h"
#include "AliESDVZERO.h"
#include "AliMultiplicity.h" 

ClassImp(AliGlobalQADataMaker)
 
//____________________________________________________________________________ 
void AliGlobalQADataMaker::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kGLOBAL, task, list) ;  
}

//____________________________________________________________________________ 
void AliGlobalQADataMaker::InitRaws()
{
  // create Raws histograms in Raws subdir
}

//____________________________________________________________________________
void AliGlobalQADataMaker::InitRecoParams() 
{
  // Get the recoparam form the OCDB 
  if (!fRecoParam) {
    TString name("GRP") ; 
    AliDebug(AliQAv1::GetQADebugLevel(), Form("Loading reconstruction parameter objects for detector %s", name.Data()));
    AliCDBPath path(name.Data(),"Calib","RecoParam");
    AliCDBEntry *entry=AliCDBManager::Instance()->Get(path.GetPath());
    if(!entry) {
      fRecoParam = NULL ; 
      AliDebug(AliQAv1::GetQADebugLevel(), Form("Couldn't find RecoParam entry in OCDB for detector %s",name.Data()));
    }
    else {
      TObject * recoParamObj = entry->GetObject() ; 
      if ( strcmp(recoParamObj->ClassName(), "TObjArray") == 0 ) {
          // The detector has only one set of reco parameters
        AliDebug(AliQAv1::GetQADebugLevel(), Form("Array of reconstruction parameters found for detector %s",name.Data()));
        TObjArray *recoParamArray = static_cast<TObjArray*>(recoParamObj) ;
        for (Int_t iRP=0; iRP<recoParamArray->GetEntriesFast(); iRP++) {
          fRecoParam = static_cast<AliDetectorRecoParam*>(recoParamArray->At(iRP)) ;
          if (!fRecoParam) 
            break ; 
          else if (fRecoParam->IsDefault()) 
            break ; 
        }
      }
      else if (recoParamObj->InheritsFrom("AliDetectorRecoParam")) {
          // The detector has only one set of reco parameters
          // Registering it in AliRecoParam
        AliDebug(AliQAv1::GetQADebugLevel(), Form("Single set of reconstruction parameters found for detector %s",name.Data()));
        fRecoParam = static_cast<AliDetectorRecoParam*>(recoParamObj) ;
        static_cast<AliDetectorRecoParam*>(recoParamObj)->SetAsDefault();
      } else { 
        AliError(Form("No valid RecoParam object found in the OCDB for detector %s",name.Data()));
      }
    }
  }
}

//____________________________________________________________________________ 
void AliGlobalQADataMaker::InitRecPoints() {
  //------------------------------------------------------
  // This function books the histograms of *track*residuals*
  // as a part of global QA
  //------------------------------------------------------
  static Bool_t first = kTRUE ; 
  if ( ! first ) 
    return ; 
  const Char_t *name[]={
    "hGlobalSPD1ResidualsY","SPD1ResidualsZ",
    "hGlobalSPD2ResidualsY","SPD2ResidualsZ",
    "hGlobalSDD1ResidualsY","SDD1ResidualsZ",
    "hGlobalSDD2ResidualsY","SDD2ResidualsZ",
    "hGlobalSSD1ResidualsY","SSD1ResidualsZ",
    "hGlobalSSD2ResidualsY","SSD2ResidualsZ",
    
    "hGlobalTPC1ResidualsY","TPC1ResidualsZ",
    "hGlobalTPC2ResidualsY","TPC2ResidualsZ",
    
    "hGlobalTRD1ResidualsY","TRD1ResidualsZ",
    "hGlobalTRD2ResidualsY","TRD2ResidualsZ",
    "hGlobalTRD3ResidualsY","TRD3ResidualsZ",
    "hGlobalTRD4ResidualsY","TRD4ResidualsZ",
    "hGlobalTRD5ResidualsY","TRD5ResidualsZ",
    "hGlobalTRD6ResidualsY","TRD6ResidualsZ",
    
    "hGlobalTOFResidualsY","TOFResidualsZ",
    
    "hGlobalPHOS1ResidualsY","PHOS1ResidualsZ",
    "hGlobalPHOS2ResidualsY","PHOS2ResidualsZ",
    
    "hGlobalHMPIDResidualsY","HMPIDResidualsZ",
    
    "hGlobalMUONResidualsY","MUONResidualsZ",
    
    "hGlobalEMCALResidualsY","EMCALResidualsZ"
  };
  const Char_t *title[]={
    "SPD1 residuals Y","SPD1 residuals Z",
    "SPD2 residuals Y","SPD2 residuals Z",
    "SDD1 residuals Y","SDD1 residuals Z",
    "SDD2 residuals Y","SDD2 residuals Z",
    "SSD1 residuals Y","SSD1 residuals Z",
    "SSD2 residuals Y","SSD2 residuals Z",
    
    "TPC1 residuals Y","TPC1 residuals Z",
    "TPC2 residuals Y","TPC2 residuals Z",
    
    "TRD1 residuals Y","TRD1 residuals Z",
    "TRD2 residuals Y","TRD2 residuals Z",
    "TRD3 residuals Y","TRD3 residuals Z",
    "TRD4 residuals Y","TRD4 residuals Z",
    "TRD5 residuals Y","TRD5 residuals Z",
    "TRD6 residuals Y","TRD6 residuals Z",
    
    "TOF residuals Y","TOF residuals Z",
    
    "PHOS1 residuals Y","PHOS1 residuals Z",
    "PHOS2 residuals Y","PHOS2 residuals Z",
    
    "HMPID residuals Y","HMPID residuals Z",
    
    "MUON residuals Y","MUON residuals Z",
    
    "EMCAL residuals Y","EMCAL residuals Z"
  };
  
  for (Int_t m=1; m<AliGeomManager::kLastLayer; m++) {
    Int_t i=2*m-2;
    TH1F *h=new TH1F(name[i],title[i],100,-5.,5.);
    Add2RecPointsList(h,i);    
    h=new TH1F(name[i+1],title[i+1],100,-5.,5.);
    Add2RecPointsList(h,i+1);    
  }

  Add2RecPointsList(
  new TH1F("hGlobalSSD1AbsoluteResidualsYNegZ",
           "SSD1 Absolute Residuals Y Neg Z",100,-2.,2.),40);
  Add2RecPointsList(
  new TH1F("hGlobalSSD1AbsoluteResidualsZNegZ",
           "SSD1 Absolute Residuals Z Neg Z",100,-2.,2.),41);
  Add2RecPointsList(
  new TH1F("hGlobalSSD1AbsoluteResidualsYPosZ",
           "SSD1 Absolute Residuals Y Pos Z",100,-2.,2.),42);
  Add2RecPointsList(
  new TH1F("hGlobalSSD1AbsoluteResidualsZPosZ",
           "SSD1 Absolute Residuals Z Pos Z",100,-2.,2.),43);


  Add2RecPointsList(
  new TH1F("hGlobalSSD2AbsoluteResidualsYNegZ",
           "SSD2 Absolute Residuals Y Neg Z",100,-3.,3.),44);
  Add2RecPointsList(
  new TH1F("hGlobalSSD2AbsoluteResidualsZNegZ",
           "SSD2 Absolute Residuals Z Neg Z",100,-3.,3.),45);
  Add2RecPointsList(
  new TH1F("hGlobalSSD2AbsoluteResidualsYPosZ",
           "SSD2 Absolute Residuals Y Pos Z",100,-3.,3.),46);
  Add2RecPointsList(
  new TH1F("hGlobalSSD2AbsoluteResidualsZPosZ",
           "SSD2Absolute Residuals Z Pos Z",100,-3.,3.),47);
  
  first = kFALSE ; 
}

//____________________________________________________________________________ 
void AliGlobalQADataMaker::InitESDs() {
  //------------------------------------------------------
  // This function books the ESD QA histograms
  // as a part of global QA
  //------------------------------------------------------

  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ;

  {// Event related QA
    const Char_t *name[]={
      "hGlobalPrimaryVertex"
    };
    const Char_t *title[]={
      "Z-distribution of the primary vertex"
    };
    Add2ESDsList(new TH1F(name[0],title[0],100,-20.,20.),kEvt0,!expert,image);
  }
 
  {// Cluster related QA
    const Char_t *name[]={
      "hGlobalFractionAssignedClustersITS",
      "hGlobalFractionAssignedClustersTPC",
      "hGlobalFractionAssignedClustersTRD",
      "hGlobalClustersPerITSModule"
    };
    const Char_t *title[]={
      "Fraction of the assigned clusters in ITS",
      "Fraction of the assigned clusters in TPC",
      "Fraction of the assigned clusters in TRD",
      "Number of clusters per an ITS module"
    };
    Add2ESDsList(new TH1F(name[0],title[0],100,0.,2.),kClr0, !expert, image);
    Add2ESDsList(new TH1F(name[1],title[1],100,0.,2.),kClr1, !expert, image);
    Add2ESDsList(new TH1F(name[2],title[2],100,0.,2.),kClr2, !expert, image);
    Add2ESDsList(new TH1F(name[3],title[3],2201,-0.5,2200.5),kClr3, !expert, image);
  }

  {// Track related QA
    const Char_t *name[]={
      "hGlobalTrackAzimuthe",                               // kTrk0
      "hGlobalTrackEta",                                    // kTrk1
      "hGlobalTPCTrackpT",                                  // kTrk2
      "hGlobalTPCITSMatchedpT",                             // kTrk3
      "hGlobalTPCTOFMatchedpT",                             // kTrk4
      "hGlobalTPCITSMatchingProbability",                   // kTrk5
      "hGlobalTPCTOFMatchingProbability",                   // kTrk6
      "hGlobalTPCsideAposDCA",                              // kTrk7
      "hGlobalTPCsideAnegDCA",                              // kTrk8
      "hGlobalTPCsideCposDCA",                              // kTrk9
      "hGlobalTPCsideCnegDCA"                               // kTrk10
  };
    const Char_t *title[]={
      "Track azimuthal distribution (rad)",                   // kTrk0
      "Track pseudo-rapidity distribution",                   // kTrk1
      "TPC: track momentum distribution (GeV)",               // kTrk2
      "TPC-ITS matched: track momentum distribution (GeV)",   // kTrk3
      "TPC-TOF matched: track momentum distribution (GeV)",   // kTrk4
      "TPC-ITS track-matching probability",                   // kTrk5
      "TPC-TOF track-matching probability",                   // kTrk6
      "TPC side A: DCA for the positive tracks (mm)",         // kTrk7
      "TPC side A: DCA for the negative tracks (mm)",         // kTrk8
      "TPC side C: DCA for the positive tracks (mm)",         // kTrk9
      "TPC side C: DCA for the negative tracks (mm)"          // kTrk10
    };
  Add2ESDsList(new TH1F(name[0],title[0],100, 0.,TMath::TwoPi()),kTrk0, !expert, image);
  Add2ESDsList(new TH1F(name[1],title[1],100,-2.00,2.00),kTrk1, !expert, image);
  Add2ESDsList(new TH1F(name[2],title[2],50,  0.20,5.00),kTrk2, !expert, image);
  Add2ESDsList(new TH1F(name[3],title[3],50,  0.20,5.00),kTrk3, !expert, image);
  Add2ESDsList(new TH1F(name[4],title[4],50,  0.20,5.00),kTrk4, !expert, image);
  Add2ESDsList(new TH1F(name[5],title[5],50,  0.20,5.00),kTrk5, !expert, image);
  Add2ESDsList(new TH1F(name[6],title[6],50,  0.20,5.00),kTrk6, !expert, image);
  Add2ESDsList(new TH1F(name[7],title[7],50, -25.0,25.0),kTrk7, !expert, image);
  Add2ESDsList(new TH1F(name[8],title[8],50, -25.0,25.0),kTrk8, !expert, image);
  Add2ESDsList(new TH1F(name[9],title[9],50, -25.0,25.0),kTrk9, !expert, image);
  Add2ESDsList(new TH1F(name[10],title[10],50, -25.0,25.0),kTrk10, !expert, image);
  }

  {// V0 related QA
    const Char_t *name[]={
      "hGlobalPromptK0sMass",
      "hGlobalOfflineK0sMass",
      "hGlobalPromptLambda0Lambda0BarMass",
      "hGlobalOfflineLambda0Lambda0BarMass"
    };
    const Char_t *title[]={
      "On-the-fly K0s mass (GeV)",
      "Offline K0s mass (GeV)",
      "On-the-fly Lambda0 + Lambda0Bar mass (GeV)",
      "Offline Lambda0 + Lambda0Bar mass (GeV)"
    };
    Add2ESDsList(new TH1F(name[0],title[0],50,  0.4477,0.5477),kK0on, !expert, image);
    Add2ESDsList(new TH1F(name[1],title[1],50,  0.4477,0.5477),kK0off, !expert, image);
    Add2ESDsList(new TH1F(name[2],title[2],50,  1.0657,1.1657),kL0on, !expert, image);
    Add2ESDsList(new TH1F(name[3],title[3],50,  1.0657,1.1657),kL0off, !expert, image);
  }

  {// PID related QA
  const Char_t *name[]={
    "hGlobalITSdEdx",
    "hGlobalTPCdEdx",
    "hGlobalTOFTrackingvsMeasured",
    "hGlobalTPCdEdxvsMomentum"
   };
    const Char_t *title[]={
      "ITS: dEdx (ADC) for particles with momentum 0.4 - 0.5 (GeV)",
      "TPC: dEdx (ADC) for particles with momentum 0.4 - 0.5 (GeV)",
      "TOF: tracking - measured (ps)",
      "TPC: dEdx (A.U.) vs momentum (GeV)"
     };
    Add2ESDsList(new TH1F(name[0],title[0],50,0.00,200.),kPid0, !expert, image);
    Add2ESDsList(new TH1F(name[1],title[1],50,0.00,100.),kPid1, !expert, image);
    Add2ESDsList(new TH1F(name[2],title[2],50,-3500.,3500.),kPid2, !expert, image);
    Add2ESDsList(new TH2F(name[3],title[3],1500,0.05,15.,700,0.,700.),kPid3,!expert,image);
   }
  {// Multiplicity related QA
    const Char_t *name[]={
      "hGlobalV0AvsITS",
      "hGlobalV0CvsITS"
    };
    const Char_t *title[]={
      "Multiplicity: V0A vs ITS",
      "Multiplicity: V0C vs ITS"
    };
    TH2F *h0=new TH2F(name[0],title[0],41,-0.5,40.5, 33,-0.5,32.5);
    Add2ESDsList(h0,kMlt0, !expert, image);
    TH2F *h1=new TH2F(name[1],title[1],41,-0.5,40.5, 33,-0.5,32.5);
    Add2ESDsList(h1,kMlt1, !expert, image);
  }

}

//____________________________________________________________________________
void AliGlobalQADataMaker::MakeRaws(AliRawReader* rawReader)
{
  //Fill prepared histograms with Raw digit properties
  rawReader->Reset() ;

}

//____________________________________________________________________________ 
void AliGlobalQADataMaker::MakeESDs(AliESDEvent * event) {
  //-----------------------------------------------------------
  // This function fills the ESD QA histograms
  // as a part of global QA
  //-----------------------------------------------------------

  const AliESDEvent *esd=event;

  // Event related QA
  const AliESDVertex *vtx=esd->GetPrimaryVertex();
  if (!vtx->GetStatus()) return;

  Double_t xv=vtx->GetXv();
  Double_t yv=vtx->GetYv();
  Double_t zv=vtx->GetZv();
  GetESDsData(kEvt0)->Fill(zv);


  Int_t ntrk=esd->GetNumberOfTracks() ; 
  for (Int_t i=0; i<ntrk; i++) {
    const AliESDtrack *track=esd->GetTrack(i);

    // Cluster related QA
    if (track->IsOn(AliESDtrack::kITSrefit)) {
      Int_t n=track->GetITSclusters(0);
      GetESDsData(kClr0)->Fill(Float_t(n)/6.); //6 is the number of ITS layers
    }

    for (Int_t j=0; j<6; ++j) {
      Int_t idet, sts;
      Float_t xloc,zloc;
      if (!track->GetITSModuleIndexInfo(j,idet,sts,xloc,zloc)) continue;
      if (j>=2) idet+=240;
      if (j>=4) idet+=260;
      if ((sts==1)||(sts==2)||(sts==4)) GetESDsData(kClr3)->Fill(idet);  
    }

    if (track->IsOn(AliESDtrack::kTPCrefit)) {
      Int_t n =track->GetTPCNcls();
      Int_t nf=track->GetTPCNclsF();      // number of crossed TPC pad rows
      if (nf>0) {
        Double_t val = n*1.0/nf; 
        GetESDsData(kClr1)->Fill(val); 
      }
    }

    if (track->IsOn(AliESDtrack::kTRDrefit)) {
      Int_t n=track->GetTRDclusters(0);
      GetESDsData(kClr2)->Fill(Float_t(n)/(6*24));//(6*24) is the number of TRD time bins
    }

    Double_t p=track->GetP();

    // Track related QA
    if (track->IsOn(AliESDtrack::kTPCrefit)) {
      Float_t dz[2]; 
      track->GetDZ(xv,yv,zv,esd->GetMagneticField(),dz); 
      if ((TMath::Abs(dz[0])<3.) && (TMath::Abs(dz[1])<3.)) { // beam pipe
        Double_t phi=track->Phi();
	GetESDsData(kTrk0)->Fill(phi);
	Double_t y=track->Eta();
	GetESDsData(kTrk1)->Fill(y);

        if (TMath::Abs(y)<0.9) {
	   GetESDsData(kTrk2)->Fill(p);
	   if (track->IsOn(AliESDtrack::kITSrefit)) GetESDsData(kTrk3)->Fill(p);
	  //if (track->IsOn(AliESDtrack::kTOFout)) GetESDsData(kTrk4)->Fill(p);
	   if (track->GetTOFsignal()>0) GetESDsData(kTrk4)->Fill(p);
	}
      }
    }
    const AliExternalTrackParam *tpcTrack=track->GetTPCInnerParam();
    const AliExternalTrackParam *innTrack=track->GetInnerParam();
    if (tpcTrack)
    if (innTrack) {
       Float_t dz[2];
       tpcTrack->GetDZ(xv,yv,zv,esd->GetMagneticField(),dz);
       dz[0]*=10.; // in mm
       if (innTrack->GetZ()  > 0)
       if (innTrack->GetTgl()> 0) { // TPC side A
	  if (tpcTrack->GetSign() > 0) GetESDsData(kTrk7)->Fill(dz[0]);
          else                         GetESDsData(kTrk8)->Fill(dz[0]);
       }
       if (innTrack->GetZ()  < 0)
       if (innTrack->GetTgl()< 0) { // TPC side C
	  if (tpcTrack->GetSign() > 0) GetESDsData(kTrk9)->Fill(dz[0]);
          else                         GetESDsData(kTrk10)->Fill(dz[0]);
       }
    }

    // PID related QA
    if ((p>0.4) && (p<0.5)) {
      if (track->IsOn(AliESDtrack::kITSpid)) {
	Double_t dedx=track->GetITSsignal();
        GetESDsData(kPid0)->Fill(dedx);
      }
      if (track->IsOn(AliESDtrack::kTPCpid)) {
	Double_t dedx=track->GetTPCsignal();
        GetESDsData(kPid1)->Fill(dedx);
      }
    }
    if (p>1.0) {
      if (track->IsOn(AliESDtrack::kITSrefit))
      if (track->IsOn(AliESDtrack::kTPCrefit))
      if (track->IsOn(AliESDtrack::kTOFout)) {
         Float_t dz[2];
         track->GetDZ(xv,yv,zv,esd->GetMagneticField(),dz);
         if (dz[1]<3.) {
            Double_t times[10];
            track->GetIntegratedTimes(times);
            Double_t tof=track->GetTOFsignal()/*-847055 -1771207*/;
            GetESDsData(kPid2)->Fill(times[2]-tof);
	 }
      }
    }
    const AliExternalTrackParam *par=track->GetInnerParam();
    if (par) {
      Double_t pp=par->GetP();
      Double_t dedx=track->GetTPCsignal();
      TH2F *h = dynamic_cast<TH2F*>(GetESDsData(kPid3));
      if (h) h->Fill(pp,dedx);
    }
 
  }

  // Multiplicity related QA
  AliESDVZERO     *mltV0 =esd->GetVZEROData();
  const AliMultiplicity *mltITS=esd->GetMultiplicity();
  if (mltV0)
    if (mltITS) {
       Short_t nv0a=mltV0->GetNbPMV0A();
       Short_t nv0c=mltV0->GetNbPMV0C();
       Int_t   nits=mltITS->GetNumberOfTracklets();
       TH2F *h0=dynamic_cast<TH2F*>(GetESDsData(kMlt0));
       if (h0) h0->Fill(nits,nv0a);
       TH2F *h1=dynamic_cast<TH2F*>(GetESDsData(kMlt1));
       if (h1) h1->Fill(nits,nv0c);
    }


  TH1 *tpc=GetESDsData(kTrk2); tpc->Sumw2();
  TH1 *its=GetESDsData(kTrk3); its->Sumw2();
  TH1 *tof=GetESDsData(kTrk4); tof->Sumw2();
  GetESDsData(kTrk5)->Divide(its,tpc,1,1.,"b");
  GetESDsData(kTrk6)->Divide(tof,tpc,1,1.,"b");

  // V0 related QA
  Int_t nV0=esd->GetNumberOfV0s();
  for (Int_t i=0; i<nV0; i++) {
    Double_t mass;
    AliESDv0 v0(*esd->GetV0(i));

    Int_t nidx=TMath::Abs(v0.GetNindex());
    AliESDtrack *ntrack1=esd->GetTrack(nidx);
    if (!ntrack1->IsOn(AliESDtrack::kTPCrefit)) continue;

    Int_t pidx=TMath::Abs(v0.GetPindex());
    AliESDtrack *ptrack1=esd->GetTrack(pidx);
    if (!ptrack1->IsOn(AliESDtrack::kTPCrefit)) continue;

    v0.ChangeMassHypothesis(kK0Short);
    mass=v0.GetEffMass();
    if (v0.GetOnFlyStatus())
       GetESDsData(kK0on)->Fill(mass);
    else
       GetESDsData(kK0off)->Fill(mass);

    v0.ChangeMassHypothesis(kLambda0);
    mass=v0.GetEffMass();
    if (v0.GetOnFlyStatus())
       GetESDsData(kL0on)->Fill(mass);
    else
       GetESDsData(kL0off)->Fill(mass);

    v0.ChangeMassHypothesis(kLambda0Bar);
    mass=v0.GetEffMass();
    if (v0.GetOnFlyStatus())
       GetESDsData(kL0on)->Fill(mass);
    else
       GetESDsData(kL0off)->Fill(mass);
  }

}
