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
//
// Reconstruct mesons 
// in central diffractive event 
// in pp collisions
//
// Authors:
//  Xianguo Lu <lu@physi.uni-heidelberg.de>
//

#include <TH1I.h>
#include <TH2I.h>
#include <TH2D.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <THnSparse.h>

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliGeomManager.h"
#include "AliITSAlignMille2Module.h"
#include "AliITSsegmentationSPD.h"
#include "AliKFParticle.h"
#include "AliMultiplicity.h"
#include "AliPID.h"
#include "AliSPDUtils.h"
#include "AliTriggerAnalysis.h"

#include "AliAnalysisTaskDDMeson.h"

void AliAnalysisTaskDDMeson::IniTask()
{
  //
  //initialize values used in many places
  //
  fnmass = 512;
  fmass1 = 0.01 * fnmass;

  fnptv = 10;
  fptv1 = 0.2 * fnptv;

  fnpta = 51*2;
  fpta1 = 0.05 * fnpta;

  fneta = 64;
  feta = 0.04 * fneta/2.;

  fnnsel = 2;
  fnsel1 = 2+fnnsel;

  fncts = 10;
  fnctlab = 20;

  fCHECKVBA = 2;
  fCHECKVBC = 4;

  //------------------
  fOpt.ToUpper();
}

AliAnalysisTaskDDMeson::AliAnalysisTaskDDMeson(const TString opt):
  AliAnalysisTaskSE("DDMeson")
  , fOpt(opt)
  , fESD(0x0)

  , fnmass(-999)
  , fmass1(-999)
  , fnptv(-999)
  , fptv1(-999)
  , fnpta(-999)
  , fpta1(-999)
  , fneta(-999)
  , feta(-999)
  , fnnsel(-999)
  , fnsel1(-999)
  , fncts(-999)
  , fnctlab(-999)
  , fCHECKVBA(-999)
  , fCHECKVBC(-999)

  , fHBIT(0x0)
  , fBitcg(0)

  , fRun(-999)
  , fat(0x0)
  , fct(0x0)
  , fbt(0x0)
  , fnt(0x0)
  , ftt(0x0)

  , fv0ntrk(0x0)
  , frsntrk(0x0)

  , fhps(0x0)
  , fhfo(0x0)
  , fhspd(0x0)
  , fhv0fmd(0x0)
  , fhpriv(0x0)
  , fhntrk(0x0)

  , fHist(0x0)
  , fThnMass(0x0)
  , fThnDPt(0x0)
  , fThnDEta(0x0)
  , fThnKF(0x0)
{
  //
  // Dummy constructor
  //
  //slot in TaskSE must start from 1
  DefineOutput(1, TList::Class());
  DefineOutput(2, THnSparse::Class());
  DefineOutput(3, THnSparse::Class());
  DefineOutput(4, THnSparse::Class());
  DefineOutput(5, THnSparse::Class());
 
  IniTask();
}

AliAnalysisTaskDDMeson::AliAnalysisTaskDDMeson(const AliAnalysisTaskDDMeson  &p):
  AliAnalysisTaskSE("DDMeson")
  , fOpt(p.fOpt)
  , fESD(p.fESD)

  , fnmass(p.fnmass)
  , fmass1(p.fmass1)
  , fnptv(p.fnptv)
  , fptv1(p.fptv1)
  , fnpta(p.fnpta)
  , fpta1(p.fpta1)
  , fneta(p.fneta)
  , feta(p.feta)
  , fnnsel(p.fnnsel)
  , fnsel1(p.fnsel1)
  , fncts(p.fncts)
  , fnctlab(p.fnctlab)
  , fCHECKVBA(p.fCHECKVBA)
  , fCHECKVBC(p.fCHECKVBC)

  , fHBIT(p.fHBIT)
  , fBitcg(p.fBitcg)

  , fRun(p.fRun)
  , fat(p.fat)
  , fct(p.fct)
  , fbt(p.fbt)
  , fnt(p.fnt)
  , ftt(p.ftt)

  , fv0ntrk(p.fv0ntrk)
  , frsntrk(p.frsntrk)

  , fhps(p.fhps)
  , fhfo(p.fhfo)
  , fhspd(p.fhspd)
  , fhv0fmd(p.fhv0fmd)
  , fhpriv(p.fhpriv)
  , fhntrk(p.fhntrk)

  , fHist(p.fHist)
  , fThnMass(p.fThnMass)
  , fThnDPt(p.fThnDPt)
  , fThnDEta(p.fThnDEta)
  , fThnKF(p.fThnKF)
{
  //
  // Copy constructor
  //

  IniTask();
}

AliAnalysisTaskDDMeson & AliAnalysisTaskDDMeson::operator=(const AliAnalysisTaskDDMeson  &p)
{
  //
  // overload =
  //
  fESD = p.fESD;

  fnmass = p.fnmass;
  fmass1 = p.fmass1;
  fnptv = p.fnptv;
  fptv1 = p.fptv1;
  fnpta = p.fnpta;
  fpta1 = p.fpta1;
  fneta = p.fneta;
  feta = p.feta;
  fnnsel = p.fnnsel;
  fnsel1 = p.fnsel1;
  fncts = p.fncts;
  fnctlab = p.fnctlab;
  fCHECKVBA = p.fCHECKVBA;
  fCHECKVBC = p.fCHECKVBC;
  
  fHBIT = p.fHBIT;
  fBitcg = p.fBitcg;

  fRun = p.fRun;

  fat = p.fat;
  fct = p.fct;
  fbt = p.fbt;
  fnt = p.fnt;
  ftt = p.ftt;

  fv0ntrk = p.fv0ntrk;
  frsntrk = p.frsntrk;

  fhps = p.fhps;
  fhfo = p.fhfo;
  fhspd = p.fhspd;
  fhv0fmd = p.fhv0fmd;
  fhpriv = p.fhpriv;
  fhntrk = p.fhntrk;

  fHist = p.fHist;
  fThnMass = p.fThnMass;
  fThnDPt = p.fThnDPt;
  fThnDEta = p.fThnDEta;
  fThnKF = p.fThnKF;

  IniTask();

  return *this;
}

AliAnalysisTaskDDMeson::~AliAnalysisTaskDDMeson()
{
  //
  //Destructor
  //
  if(fESD) delete fESD;

  delete fHBIT;
  delete fat;
  delete fct;
  delete fbt;
  delete fnt;
  delete ftt;

  delete fv0ntrk;
  delete frsntrk;

  delete fhps;
  delete fhfo;
  delete fhspd;
  delete fhv0fmd;
  delete fhpriv;
  delete fhntrk;

  delete fThnMass;
  delete fThnDPt;
  delete fThnDEta;
  delete fThnKF;

  fHist->Clear();
  delete fHist;  
 
}

void AliAnalysisTaskDDMeson::CheckRange(Double_t &ptv, Double_t &pta, Double_t &etaa
                                        , Double_t &mpi
                                        ) const 
{
  //
  //save over/under flow bins
  //

  const Double_t eps = 1e-6;
  if(ptv>=fptv1) ptv = fptv1 - eps;
  if(pta>=fpta1) pta = fpta1 - eps;

  if(etaa >= feta) etaa = feta - eps;

  if(etaa <= -feta) etaa = -feta + eps;

  if(mpi >= fmass1) mpi = fmass1 - eps;

}

void AliAnalysisTaskDDMeson::UserCreateOutputObjects()
{
  //
  //CreateOutputObjects
  //
  //=======================================
  const Double_t kvz0=1, kvz1=5;
  const Int_t nkvz=(Int_t)(kvz1-kvz0);

  const Double_t charge0=1, charge1=4;
  const Int_t ncharge=(Int_t)(charge1-charge0);
  
  const Double_t mass0 = 0;
  //=======================================
  //total momentum and pt in lab frame:
  const Double_t p0=0;
  //=======================================
  const Double_t nsel0 = 2;

  //=======================================
  const Double_t cts0=0, cts1=1;

  //=======================================
  const Double_t ctlab0 = -1, ctlab1 = 1;

  //=======================================
  //------- Mass
  const Int_t    nparPtMass = 7;
  const Int_t    nbinPtMass[]     ={fnnsel, nkvz, ncharge,  fnmass, fnptv, fncts, fnctlab};
  const Double_t xminPtMass[]     ={nsel0,  kvz0, charge0,   mass0,    p0,  cts0,  ctlab0};
  const Double_t xmaxPtMass[]     ={fnsel1, kvz1, charge1,  fmass1, fptv1,  cts1,  ctlab1};
  fThnMass = new THnSparseD("DDMeson_Mass","nsel, kvz, charge, Mpipi, ptv, cts, ctlab", nparPtMass, nbinPtMass, xminPtMass, xmaxPtMass);
 
  //------- DPt
  const Int_t    nparDPt = 4;
  const Int_t    nbinDPt[] ={fnnsel, nkvz, ncharge,  fnpta};
  const Double_t xminDPt[] ={nsel0, kvz0, charge0,  p0};
  const Double_t xmaxDPt[] ={fnsel1, kvz1, charge1,  fpta1};
  fThnDPt = new THnSparseD("DDMeson_DPt","nsel, kvz, charge, pt1", nparDPt, nbinDPt, xminDPt, xmaxDPt);
  //------- DEta
  const Int_t    nparDEta = 4;
  const Int_t    nbinDEta[] ={fnnsel, nkvz, ncharge,  fneta};
  const Double_t xminDEta[] ={nsel0, kvz0, charge0,  -feta};
  const Double_t xmaxDEta[] ={fnsel1, kvz1, charge1,  feta};
  fThnDEta = new THnSparseD("DDMeson_DEta"," nsel, kvz, charge, eta1", nparDEta, nbinDEta, xminDEta, xmaxDEta);

  //------- KF
  const Double_t mks = 0.494;
  const Double_t dks = 0.024;
  const Int_t nparKF = 4;
  const Int_t nbinKF[] = {nkvz, 200, 200, 200};
  const Double_t xminKF[] ={kvz0, mks-dks, mks-dks, 0};
  const Double_t xmaxKF[] ={kvz1, mks+dks, mks+dks, 10};
  fThnKF = new THnSparseD("DDMeson_KF","kvz, rawks, kfks, chi2", nparKF, nbinKF, xminKF, xmaxKF);
  //=======================================
  //=======================================
  fHBIT = new TH1I("HBIT","",32,1,33);

  const Int_t runb0 = 114650; //runa1=114649, runb1=117630, runc0=117631, runc1=121526, rund0=121527, rund1=126460, rune0=126461, rune1 = 130930; //runf0=130931, all checked on elog
  const Int_t runf1 = 133900;//02OCT2010
  const Int_t run0 = runb0-1, run1 = runf1+1;
  const Int_t nrun = (Int_t)(run1-run0);

  fat = new TH1I("at","",nrun,run0,run1);
  fct = new TH1I("ct","",nrun,run0,run1);
  fbt = new TH1I("bt","",nrun,run0,run1);
  fnt = new TH1I("nt","",nrun,run0,run1);
  ftt = new TH1I("tt","",nrun,run0,run1);
  
  fv0ntrk = new TH2D("c_v0ntrk","",80,0,80, 4, 1, 5);//x: ntrk; y: GetV0
  frsntrk = new TH2D("c_rsntrk","",1000,0,1000, 200, 0,200);
 
  fhps =    new TH1I("a000_ps",    "", 20, 0,20);
  fhfo =    new TH2I("a010_fo",    "", 50, 0, 50, 50, 0, 50);
  fhspd =   new TH1I("a011_spd",   "", 50, 0, 50);
  fhv0fmd = new TH2I("a020_v0fmd", "", 4, 0, 4, 4, 0, 4);
  fhpriv =  new TH1I("a100_priv",  "", 2, 0,2);
  fhntrk =  new TH1I("a110_ntrk",  "", (Int_t)fnsel1, 0, (Int_t)fnsel1);

  //------

  fHist = new TList;

  //------>>> Very important!!!!!!!
  fHist->SetOwner();

  fHist->Add(fat);
  fHist->Add(fct);
  fHist->Add(fbt);
  fHist->Add(fnt);
  fHist->Add(ftt);

  fHist->Add(fv0ntrk);
  fHist->Add(frsntrk);
  fHist->Add(fhps);
  fHist->Add(fhfo);
  fHist->Add(fhspd);
  fHist->Add(fhv0fmd);
  fHist->Add(fhpriv);
  fHist->Add(fhntrk);
}

void AliAnalysisTaskDDMeson::UserExec(Option_t *)
{
  //
  //Execute
  // 
  //----------------------------------- general protection
  if(!CheckESD())
    return;

  //----------------------------------- trigger bits
  if(!CheckBit())
    return;  

  //----------------------------------- track cuts
  const Int_t ntrk0 = fESD->GetNumberOfTracks();
  const AliESDtrack * trks[ntrk0];
  for(Int_t ii=0; ii<ntrk0; ii++){
    trks[ii] = 0x0;
  }
  const Int_t nsel = CutESD(trks);
  if(nsel<2)
    return;

  for(Int_t ii=0; ii<nsel; ii++){
    for(Int_t jj=ii+1; jj<nsel; jj++){
      const AliESDtrack * trkpair[2]={trks[ii], trks[jj]};

      //----------------------------------- tracks initailization
      TRandom3 tmprd(0);
      if(tmprd.Rndm()>0.5){
        SwapTrack(trkpair);
      }
  
      //------------------------------------------------------------------------
      const Int_t pipdg = 211;
      const AliKFParticle daughter0(*(trkpair[0]), pipdg);
      const AliKFParticle daughter1(*(trkpair[1]), pipdg);

      const AliKFParticle parent(daughter0, daughter1);
      const Double_t kfchi2 = parent.GetChi2();
      const Double_t kfmass = parent.GetMass();

      //------------------------------------------------------------------------
      const Int_t kvz = GetV0();

      const AliESDtrack * trk1=trkpair[0];
      const AliESDtrack * trk2=trkpair[1];
  
      const Double_t ch1 = trk1->GetSign();
      const Double_t ch2 = trk2->GetSign();
      const Int_t combch = GetCombCh(ch1, ch2);
      
      Double_t pta  = trk1->Pt();
      Double_t etaa = trk1->Eta();
      
      Double_t tmpp1[3], tmpp2[3];
      trk1->GetPxPyPz(tmpp1);
      trk2->GetPxPyPz(tmpp2);
  
      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      const Double_t ctlab = GetCtlab(tmpp1, tmpp2);

      Double_t cts = -999;
      const Double_t masspi = AliPID::ParticleMass(AliPID::kPion);
      const TLorentzVector sumv = GetKinematics(tmpp1, tmpp2, masspi, masspi, cts);
      Double_t pipi = sumv.Mag();
      Double_t ptv = sumv.Pt();

      CheckRange(ptv, pta, etaa, pipi);
      
      const Double_t varMass[]  ={nsel, kvz, combch
                                       , pipi, ptv, cts
                                       , ctlab
      };
      const Double_t varDPt[]        ={nsel, kvz, combch, pta};
      const Double_t varDEta[]       ={nsel, kvz, combch, etaa};
      fThnMass->Fill(varMass);
      fThnDPt->Fill(varDPt);
      fThnDEta->Fill(varDEta);

      const Double_t varKF[] = {kvz, pipi, kfmass, kfchi2};
      if(nsel==2 && combch==3)
        fThnKF->Fill(varKF);
    }
  }
  //------------------------------------------------------------------------
  //------------------------------------------------------------------------

  PostData(1, fHist);
  PostData(2, fThnMass);
  PostData(3, fThnDPt);
  PostData(4, fThnDEta);
  PostData(5, fThnKF);
}
//=======================================================================
//=======================================================================
    
void AliAnalysisTaskDDMeson::Terminate(Option_t *)
{
  //
  //fillbit before leaving
  //
  FillBit();
}

void AliAnalysisTaskDDMeson::FillBit()
{
  //
  //save the bit count in f*t 
  //
  Double_t tot[5];
  CalcBit(fHBIT, tot);

  const Int_t bin = fat->FindBin(fRun);
  fat->AddBinContent(bin, tot[0]);
  fct->AddBinContent(bin, tot[1]);
  fbt->AddBinContent(bin, tot[2]);
  fnt->AddBinContent(bin, tot[3]);
  ftt->AddBinContent(bin, tot[4]);
}

void AliAnalysisTaskDDMeson::CalcBit(TH1I *hc, Double_t tot[]) const
{
  //
  //calculate the bit count in f*t 
  //
  if(hc->GetBinLowEdge(1)!=1){
    printf("xlulog hc defined incorrectly!! %f\n",hc->GetBinLowEdge(1));
    return;
  }
  if(hc->GetBinContent(0)){
    printf("\nxlulog hc !!!!!!!!!!!!!!!!!!! underflow!!!!!!!!!!!!!!!!! %f\n\n", hc->GetBinContent(0));
    return;
  }
  if(hc->GetBinContent(hc->GetNbinsX())){
    printf("\nxlulog hc !!!!!!!!!!!!!!!!! OVERflow!!!!!!!!!!!!!! %f\n\n", hc->GetBinContent(hc->GetNbinsX()));
    return;
  }


  Double_t atot=0, ctot=0, btot=0, ntot=0;
  for(Int_t ib=1; ib<=hc->GetNbinsX(); ib++){
    const Double_t cont = hc->GetBinContent(ib);
    if(       (ib & fCHECKVBA) && !(ib & fCHECKVBC) ){
      atot += cont;
    }
    else if( !(ib & fCHECKVBA) &&  (ib & fCHECKVBC) ){
      ctot += cont;
    }
    else if(  (ib & fCHECKVBA) &&  (ib & fCHECKVBC) ){
      btot += cont;
    }
    else{
      ntot += cont;
    }
  }

  const Double_t ttot   = atot + ctot + btot + ntot;
  
  tot[0]=atot;
  tot[1]=ctot;
  tot[2]=btot;
  tot[3]=ntot;
  tot[4]=ttot;

  //---
  char fullname[5][20];
  strncpy(fullname[0],"V0A-ONLY",19);
  strncpy(fullname[1],"V0C-ONLY",19);
  strncpy(fullname[2],"V0A & V0C",19);
  strncpy(fullname[3],"!V0A & !V0C",19);
  strncpy(fullname[4],"TOTAL",19);
  for (Int_t ii=0;ii<5;ii++){
    fullname[ii][19] = '\0';
  }

  for(Int_t ii=0;ii<5;ii++){
    printf("xlulog %s CTP-L0: %s\tTOT: %10.0f\n",fOpt.Data(), fullname[ii],tot[ii]);
  }

  //---

  const Int_t nb=hc->GetNbinsX();
  for(Int_t ib=0; ib<=nb+1; ib++){
    hc->SetBinContent(ib, 0);
  }
}

//====================================================================================

Bool_t AliAnalysisTaskDDMeson::CheckESD()
{
  //
  //general protection
  //
  const AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler *>(fInputHandler);
  if(esdH)  
    fESD = esdH->GetEvent();
  if(!fESD){
    AliError("xlulog No ESD Event");
    return 0;
  }
  
  if(fabs(fESD->GetMagneticField())<1){
    printf("xlulog strange Bfield! %f\n", fESD->GetMagneticField());
    return 0;
  }

  const Int_t tmprun = fESD->GetRunNumber();

  if(fRun!=tmprun){
    if(fRun>0){
      printf("xlulog run changed during exec!! %d %d\n", tmprun, fRun);
      FillBit();
    }
    fRun = tmprun;
    if(fOpt.Contains("SPD")){
      SPDLoadGeom();
    }
  }

  return 1;  
}

void AliAnalysisTaskDDMeson::SPDLoadGeom() const 
{
  //method to get the gGeomanager
  // it is called at the CreatedOutputObject stage
  // to comply with the CAF environment
  AliCDBManager *man = AliCDBManager::Instance();
  TString cdbpath("local:///lustre/alice/alien/alice/data/2010/OCDB");
  man->SetDefaultStorage(cdbpath);
  man->SetRun(fRun);  
  man->SetSpecificStorage("ITS/Align/Data",cdbpath);
  man->SetSpecificStorage("GRP/Geometry/Data",cdbpath);
  
  AliCDBEntry* obj = man->Get(AliCDBPath("GRP", "Geometry", "Data"));
  if (obj != NULL) {
    AliGeomManager::SetGeometry((TGeoManager*)obj->GetObject());
    AliGeomManager::ApplyAlignObjsFromCDB("ITS"); 
    printf("xlulog DDMeson %s Geom Loaded for run %d !\n", fOpt.Data(), fRun);
  }
  else {
    printf("xlulog DDMeson %s was unable to load Geom for run %d !\n", fOpt.Data(), fRun);
  }
}

Bool_t AliAnalysisTaskDDMeson::SPDLoc2Glo(const Int_t id, const Double_t *loc, Double_t *glo) const 
{
  //
  //SPDLoc2Glo
  //

  static TGeoHMatrix mat;
  Int_t vid = AliITSAlignMille2Module::GetVolumeIDFromIndex(id);
  if (vid<0) {
    printf("xlulog Did not find module with such ID %d\n",id);
    return kFALSE;
  }
  AliITSAlignMille2Module::SensVolMatrix(vid,&mat);
  mat.LocalToMaster(loc,glo);
  return kTRUE;
}

Bool_t AliAnalysisTaskDDMeson::CheckChipEta(const Int_t chipKey) const 
{
  //
  //CheckChipEta
  //

  //no eta cut 
  if(fOpt.Contains("SPD0"))
    return kTRUE;

  const Double_t etacut = 1.0;

  UInt_t module=999, offchip=999;
  AliSPDUtils::GetOfflineFromOfflineChipKey(chipKey,module,offchip);
  UInt_t hs = AliSPDUtils::GetOnlineHSFromOffline(module);
  if(hs<2) offchip = 4 - offchip; // inversion  in the inner layer... 
  
  const Int_t col[]={
    hs<2? 0 : 31, 
    hs<2? 31 : 0, 
    hs<2? 31 : 0, 
    hs<2? 0 : 31};
  const Int_t aa[]={0, 0, 255, 255};
  const AliITSsegmentationSPD seg;

  for(Int_t ic=0; ic<4; ic++){
    Float_t localchip[3]={0.,0.,0.};
    seg.DetToLocal(aa[ic],col[ic]+32*offchip,localchip[0],localchip[2]); // local coordinate of the chip center
    //printf("local coordinates %d %d: %f %f \n",chipKey, ic, localchip[0],localchip[2]);
    const Double_t local[3] = {localchip[0],localchip[1],localchip[2]};
    Double_t glochip[3]={0.,0.,0.};
    if(!SPDLoc2Glo(module,local,glochip)){
      return kFALSE;
    }

    const TVector3 pos(glochip[0],glochip[1],glochip[2]);
    //pos.Print();

    if(fabs(pos.Eta()) > etacut)
      return kFALSE;
  }

  return kTRUE;
}

void AliAnalysisTaskDDMeson::GetNFO(Int_t &ni, Int_t &no) const 
{
  //
  //GetNFO
  //

  const AliMultiplicity *mult = fESD->GetMultiplicity();
  Int_t nFOinner=0;
  Int_t nFOouter=0;
  for(Int_t iChipKey=0; iChipKey < 1200; iChipKey++){
    if(mult->TestFastOrFiredChips(iChipKey) && CheckChipEta(iChipKey)){ // here you check if the FastOr bit is 1 or 0
      if(iChipKey<400) 
        nFOinner++;  // here you count the FastOr bits in the inner layer
      else 
        nFOouter++;  // here you count the FastOr bits in the outer layer
    }
  }

  ni = nFOinner;
  no = nFOouter;
}

Bool_t AliAnalysisTaskDDMeson::CheckBit()
{
  //
  //get offline trigger
  //
  fhps->Fill(1);

  AliTriggerAnalysis triggerAnalysis;

  //hardware 1: hw 0: offline
  const Bool_t khw = kFALSE;
  const Bool_t v0A       = (triggerAnalysis.V0Trigger(fESD, AliTriggerAnalysis::kASide, khw) == AliTriggerAnalysis::kV0BB);
  const Bool_t v0C       = (triggerAnalysis.V0Trigger(fESD, AliTriggerAnalysis::kCSide, khw) == AliTriggerAnalysis::kV0BB);

  /*
    const Bool_t v0ABG     = (triggerAnalysis.V0Trigger(fESD, AliTriggerAnalysis::kASide, khw) == AliTriggerAnalysis::kV0BG);
    const Bool_t v0CBG     = (triggerAnalysis.V0Trigger(fESD, AliTriggerAnalysis::kCSide, khw) == AliTriggerAnalysis::kV0BG);
    const Bool_t v0BG      = v0ABG || v0CBG;
    const Int_t fastOROffline = triggerAnalysis.SPDFiredChips(fESD, 0); // SPD number of chips from clusters (!)        

    //http://alisoft.cern.ch/viewvc/trunk/ANALYSIS/AliPhysicsSelection.cxx?view=markup&root=AliRoot
    //Bool_t mb1 = (fastOROffline > 0 || v0A || v0C) && (!v0BG);
   */

  const Bool_t fmdA = triggerAnalysis.FMDTrigger(fESD, AliTriggerAnalysis::kASide);
  const Bool_t fmdC = triggerAnalysis.FMDTrigger(fESD, AliTriggerAnalysis::kCSide);

  Bool_t bitA = kFALSE, bitC = kFALSE;
  if(fOpt.Contains("V0")){
    bitA = bitA || v0A;
    bitC = bitC || v0C;
  }
  if(fOpt.Contains("FMD")){
    bitA = bitA || fmdA;
    bitC = bitC || fmdC;
  }
  
  //---------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------
  //considering the efficiency per chip, nfo and fastORhw is required only to have non-0 count!!  

  //simulate SPD _hardware_ trigger in eta range
  if(fOpt.Contains("SPD")){
    Int_t ni=-999, no = -999;
    GetNFO(ni, no);
    if(!bitA && !bitC){
      fhfo->Fill(ni, no);
    }
    if(ni<2 || no<2)
      return 0;
  }

  //simulate MB condition, since for double gap SPD is definitely fired, fastORhw>0 to reduce GA, GC, NG
  //http://alisoft.cern.ch/viewvc/trunk/ANALYSIS/AliTriggerAnalysis.cxx?view=markup&root=AliRoot
  //Int_t AliTriggerAnalysis::SPDFiredChips(const AliESDEvent* aEsd, Int_t origin, Bool_t fillHists, Int_t layer)
  // returns the number of fired chips in the SPD
  // origin = 0 --> aEsd->GetMultiplicity()->GetNumberOfFiredChips() (filled from clusters)
  // origin = 1 --> aEsd->GetMultiplicity()->TestFastOrFiredChips() (from hardware bits)
  const Int_t fastORhw = triggerAnalysis.SPDFiredChips(fESD, 1);
  fhspd->Fill(fastORhw);
  if(fastORhw<1)
    return 0;

  //---------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------
    
  const Int_t iv0 = v0A + 2*v0C;
  const Int_t ifmd = fmdA + 2*fmdC;
  fhv0fmd->Fill(iv0, ifmd);

  //---------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------
  
  //need a base number to distringuish from null entry; fHBit starts from 1; important!!
  fBitcg=1; 

  if(bitA)
    fBitcg += fCHECKVBA;
  if(bitC)
    fBitcg += fCHECKVBC;

  fHBIT->Fill( fBitcg );
  
  return 1;
}

Int_t AliAnalysisTaskDDMeson::CutESD(const AliESDtrack *outtrk[])
{
  //
  //track cuts
  //
  
  //collision vertex cut
  //A cut in XY is implicitly done during the reconstruction by constraining the vertex to the beam diamond.

  // Primary vertex
  Bool_t kpr0 = kTRUE;
  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
  if(vertex->GetNContributors()<1) {
    // SPD vertex
    vertex = fESD->GetPrimaryVertexSPD();
    if(vertex->GetNContributors()<1) {
      // NO GOOD VERTEX, SKIP EVENT 
      kpr0 = kFALSE;
    }
  }
  const Bool_t kpriv = kpr0 && ( fabs(fESD->GetPrimaryVertex()->GetZ())<10 ); 
  fhpriv->Fill(kpriv);
  if(!kpriv)  
    return 0;
  
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;//::GetStandardITSTPCTrackCuts2010(kTRUE);

  //--->
  // TPC  
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  // ITS
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  // 7*(0.0026+0.0050/pt^1.01)
  esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

  esdTrackCuts->SetMaxDCAToVertexZ(2);
  esdTrackCuts->SetDCAToVertex2D(kFALSE);
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //---<

  esdTrackCuts->SetEtaRange(-0.9, 0.9);
  
  const TObjArray* seltrk = esdTrackCuts->GetAcceptedTracks(fESD);
  const Int_t n0=seltrk->GetEntries();
  const AliESDtrack * trks[n0];
  Int_t nsel=0;
  for(Int_t ii=0; ii<n0; ii++){
    trks[ii]=0x0;
    const AliESDtrack *esdtrack = (AliESDtrack *)(seltrk->At(ii));

    /*
    if(!CutTrack(esdtrack))
      continue;
    */

    trks[nsel++]=esdtrack;
  }

  delete seltrk;
  delete esdTrackCuts;

  fv0ntrk->Fill(nsel,GetV0());

  fhntrk->Fill(nsel);

  frsntrk->Fill(fESD->GetNumberOfTracks(), nsel);

  if(nsel>=fnsel1)
    return 0;

  for(Int_t ii=0; ii<nsel; ii++){
    outtrk[ii]=trks[ii];
  }

  return nsel;
}

/*
Bool_t AliAnalysisTaskDDMeson::CutTrack(const AliESDtrack *  esdtrack) const 
{
  //
  //track cut
  //  

  return kTRUE;
}
*/

//=================================================================================

void AliAnalysisTaskDDMeson::SwapTrack(const AliESDtrack * trks[]) const 
{
  //
  //swap two esdtracks
  //
  const AliESDtrack * tmpt = trks[0];
  trks[0]=trks[1];
  trks[1]=tmpt;
}

Int_t AliAnalysisTaskDDMeson::GetV0() const
{
  //
  //V0 bit combination
  //

  const Bool_t ka = (fBitcg & fCHECKVBA);
  const Bool_t kc = (fBitcg & fCHECKVBC);
  

  //a1c0 a0c1 must be adjacent so that a cut for signal gas can be eaily applied
  const Int_t a0c0 = 1;
  const Int_t a1c0 = 2;
  const Int_t a0c1 = 3;
  const Int_t a1c1 = 4;

  if(!ka && !kc)
    return a0c0;
  else{
    if(ka && kc)
      return a1c1;
    else if(ka)
      return a1c0;
    else 
      return a0c1;
  }
}

Int_t AliAnalysisTaskDDMeson::GetCombCh(const Double_t s1, const Double_t s2) const
{
  //
  //get combination of charges
  //
  const Int_t combPP = 1;
  const Int_t combMM = 2;
  const Int_t combPM = 3;

  if(s1*s2<0){
    return combPM;
  }
  else if(s1>0)
    return combPP;
  else
    return combMM;
}

//==========================================================================

TLorentzVector AliAnalysisTaskDDMeson::GetKinematics(const Double_t *pa, const Double_t *pb, const Double_t ma, const Double_t mb, Double_t & cts) const
{
  //
  //get kinematics, cts = cos(theta_{#star})
  //
  TLorentzVector va, vb, sumv;
  va.SetXYZM(pa[0], pa[1], pa[2], ma);
  vb.SetXYZM(pb[0], pb[1], pb[2], mb);
  sumv = va+vb;

  const TVector3 bv = -sumv.BoostVector();

  va.Boost(bv);
  vb.Boost(bv);
  //printf("test a %f %f %f -- %f %f\n", va.X(), va.Y(), va.Z(), va.Theta(), va.Phi());
  //printf("test b %f %f %f -- %f %f\n", vb.X(), vb.Y(), vb.Z(), vb.Theta(), vb.Phi());

  cts = fabs(cos(va.Theta()));

  return sumv;
}

Double_t AliAnalysisTaskDDMeson::GetCtlab(const Double_t *pa, const Double_t *pb) const 
{
  //
  //cos(theat_lab)
  //

  TVector3 va, vb;
  va.SetXYZ(pa[0], pa[1], pa[2]);
  vb.SetXYZ(pb[0], pb[1], pb[2]);

  const TVector3 ua = va.Unit();
  const TVector3 ub = vb.Unit();

  return ua.Dot(ub);
}
