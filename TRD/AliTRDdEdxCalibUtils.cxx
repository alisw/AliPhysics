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
// class to calculate TRD dEdx
// xx
// xx
// xx
// xx
//
//  Xianguo Lu 
//  lu@physi.uni-heidelberg.de
//  Xianguo.Lu@cern.ch
//  
//


#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMinuit.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TVectorD.h"

#include "TTreeStream.h"

#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "AliCDBStorage.h"
#include "AliESDEvent.h"
#include "AliESDfriendTrack.h"
#include "AliESDtrack.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCalROC.h"
#include "AliTRDtrackV1.h"

#include "AliTRDdEdxBaseUtils.h"
#include "AliTRDdEdxReconUtils.h"
#include "AliTRDdEdxCalibHistArray.h"
#include "AliTRDdEdxCalibUtils.h"

#define EPSILON 1e-12

AliTRDdEdxCalibHistArray * AliTRDdEdxCalibUtils::fgHistArray = 0x0;
TObjArray * AliTRDdEdxCalibUtils::fgObjArray = 0x0;

//============================================================
//                           CalibObj
//============================================================
TObjArray * AliTRDdEdxCalibUtils::GetObjArray()
{
  //
  //return fgObjArray, initialized if null
  //

  if(!fgObjArray){
    fgObjArray = new TObjArray(8);
  }

  return fgObjArray;
}

TObjArray * AliTRDdEdxCalibUtils::GetObj(const Bool_t kinvq, const Double_t mag, const Int_t charge)
{
  //
  //return calib obj
  //
  if(!fgObjArray){
    printf("AliTRDdEdxCalibUtils::GetObjArray(kinvq, mag, charge) error fgObjArray null!!\n"); exit(1);
  }

  return (TObjArray*) fgObjArray->At(AliTRDdEdxCalibHistArray::GetIterator(kinvq, mag, charge));
}

void AliTRDdEdxCalibUtils::DeleteObjArray()
{
  //
  //delete calib obj
  //
  if(fgObjArray){
    fgObjArray->SetOwner();
    delete fgObjArray;
    fgObjArray = 0x0;
  }
}

Bool_t AliTRDdEdxCalibUtils::GenerateDefaultOCDB(const TString path)
{
  //
  //generate default OCDB object PHQ, do like
  //AliTRDdEdxCalibUtils::GenerateDefaultPHQOCDB("local://./")
  //

  TObjArray * arr8 = new TObjArray(8);
  arr8->SetOwner();

  for(Int_t ii=0; ii<8; ii++){
    TObjArray * arr1 = new TObjArray(1);
    arr1->SetOwner();
    TString objn(AliTRDdEdxCalibHistArray::GetNameAt(ii));
    objn.ReplaceAll("Hist","Obj");
    arr1->SetName(objn);

    const Int_t nbins = AliTRDdEdxBaseUtils::NTRDtimebin();
    TVectorD * vec = new TVectorD(nbins);
    for(Int_t jj=0; jj<nbins; jj++){
      (*vec)[jj] = 1;
    }
    arr1->Add(vec);
    arr8->Add(arr1);
  }

  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Raphaelle Bailhache and Xianguo Lu");

  AliCDBId id1("TRD/Calib/PHQ", 0, 999999999);
  AliCDBStorage * gStorage = AliCDBManager::Instance()->GetStorage(path);
  gStorage->Put(arr8, id1, metaData);

  delete metaData;
  delete arr8;

  return kTRUE;
}

//============================================================
//                           CalibHist
//============================================================
void AliTRDdEdxCalibUtils::DeleteHistArray()
{
  //
  //delete calib hist
  //
  delete fgHistArray;
  fgHistArray = 0x0;
}

THnBase * AliTRDdEdxCalibUtils::GetHistAt(const Int_t iter)
{
  //
  //
  //
  if(iter<fgHistArray->GetSize())
    return (THnBase*)fgHistArray->At(iter);
  else
    return 0x0;
}

void AliTRDdEdxCalibUtils::IniHistArray(TList *list, const Bool_t kNoInv)
{
  //
  //initialize calib hist, list should not own the hist, or list->Clear/delete hist should not be called
  //

  delete fgHistArray;
  fgHistArray = new AliTRDdEdxCalibHistArray(kNoInv);
  list->Add(fgHistArray);
}

Bool_t AliTRDdEdxCalibUtils::ReadHistArray(const TString filename, const TString listname)
{
  //
  //used in AliTRDPreprocessorOffline
  //read in calib hist from file, only for PHQ
  //

  //maybe already open by others... don't close
  TFile fcalib(filename);
  TObjArray * array = (TObjArray*)fcalib.Get(listname);
  const TString histname = AliTRDdEdxCalibHistArray::GetArrayName();
 
  AliTRDdEdxCalibHistArray * tmph=0x0;
  if(array){
    tmph = (AliTRDdEdxCalibHistArray *) array->FindObject(histname);
  }
  else{
    tmph = (AliTRDdEdxCalibHistArray *) fcalib.Get(histname);
  }
  if(!tmph){
    printf("AliTRDdEdxCalibUtils::ReadCalibHist warning calib hist not found! %s %s %s\n", filename.Data(), listname.Data(), histname.Data());
    fcalib.ls();
    if(array){
      array->ls();
    }
    return kFALSE;
  }
  
  delete fgHistArray; 
  fgHistArray = new AliTRDdEdxCalibHistArray(*tmph);

  return kTRUE;
}

void AliTRDdEdxCalibUtils::FillHist(const Int_t ncls, const TVectorD *arrayQ, const TVectorD *arrayX, THnBase * hcalib, const Double_t scale)
{
  //
  //fill calibration hist
  //
  if(!hcalib){printf("AliTRDdEdxCalibUtils::FillCalibHist errro hcalib null!!\n"); exit(1);}

  for(Int_t ii=0; ii<ncls; ii++){
    const Double_t dq = (*arrayQ)[ii];
    const Double_t xx = (*arrayX)[ii];

    const Double_t qmax = hcalib->GetAxis(1)->GetXmax() -0.5 * hcalib->GetAxis(1)->GetBinWidth(1);
    const Double_t xmin = hcalib->GetAxis(0)->GetXmin();
    const Double_t xmax = hcalib->GetAxis(0)->GetXmax();

    if(xx>=xmax || xx<xmin){
      printf("AliTRDdEdxCalibUtils::FillCalibHist error x overflow or underflow! %s %15f %15f %15f\n", hcalib->GetName(),  xx, xmin, xmax); exit(1);
    }

    const Double_t var[]={xx, TMath::Min(dq, qmax)/scale};
    hcalib->Fill(var);
  }
}

void AliTRDdEdxCalibUtils::FillHist(const AliTRDtrackV1 *trdv1, const Bool_t kinvq, const Double_t mag, const Int_t charge, const Double_t scale) 
{
  //
  //get cluster Q and fill calib hist, if kinvq = kTRUE, 1/Q is filled
  //
  if(!fgHistArray){
    printf("AliTRDdEdxCalibUtils::FillHist fgHistArray not initialized!!\n"); exit(1);
  }

  THnBase * hcalib = fgHistArray->GetHist(kinvq, mag, charge);

  TVectorD arrayQ(200), arrayX(200);
  const Int_t ncls = AliTRDdEdxReconUtils::GetArrayClusterQ(kinvq, &arrayQ, &arrayX, trdv1);
  FillHist(ncls, &arrayQ, &arrayX, hcalib, kinvq ? 1/scale : scale);

  static Int_t kprint = 100;
  if(kprint<0){
    printf("\nAliTRDdEdxCalibUtils::FillHist summary: \n");
    printf("\nkinvq= %d;\n", kinvq);
    for(Int_t iq=0; iq<ncls; iq++){
      printf("arrayX[%3d] = %15.0f; arrayQ[%3d] =  %15f;\n", iq, arrayX[iq], iq, arrayQ[iq]);
    }
    printf("\n");
  }
  kprint++;
}

//============================================================
//
//============================================================

Double_t AliTRDdEdxCalibUtils::GetCalibTPCscale(const Int_t tpcncls, const Double_t tpcsig)
{
  //
  //the scale used in calibration
  //

  if(tpcncls < AliTRDdEdxBaseUtils::CalibTPCnclsCut())
    return -999;

  if(tpcsig<EPSILON)
    return -999;

  return tpcsig/120;

}

void AliTRDdEdxCalibUtils::GetPHCountMeanRMS(const TH1D *hnor, TH1D *&hmean)
{
  //
  //calculate from the ph calib hist the (mean-3sigma) ph-count in the chamber, save in the TH1D output
  //
  const Int_t ndet = 540;
  TObjArray *obj=new TObjArray(ndet);
  obj->SetOwner();
  for(Int_t ii=0; ii<ndet; ii++){
    obj->Add(new TVectorD(AliTRDseedV1::kNtb));
  }

  //ibin = binlowedge of bin(ibin+1) = the number fills this bin
  for(Int_t ibin=0; ibin<hnor->GetNbinsX(); ibin++){
    const Double_t stat = hnor->GetBinContent(ibin+1);
    if(stat<EPSILON){
      continue;
    }

    const Int_t idet = AliTRDdEdxBaseUtils::ToDetector(ibin);
    const Int_t itb  = AliTRDdEdxBaseUtils::ToTimeBin(ibin);
    TVectorD *vec=(TVectorD *)obj->At(idet);
    (*vec)[itb] = stat;
  }

  hmean = new TH1D(Form("%sdetmean", hnor->GetName()), "", hnor->GetNbinsX(), hnor->GetXaxis()->GetXmin(), hnor->GetXaxis()->GetXmax());
  for(Int_t ibin=0; ibin<hnor->GetNbinsX(); ibin++){
    const Int_t idet = AliTRDdEdxBaseUtils::ToDetector(ibin);
    const TVectorD *vec=(TVectorD *)obj->At(idet);

    Int_t nonzero = 0;
    for(Int_t ii=0; ii<vec->GetNrows(); ii++){
      if((*vec)[ii]>EPSILON){
        nonzero++;
      }
    }

    Double_t mean = 0;
    const Double_t lowfrac = 0.6;
    //if there are too many 0's, reject this chamber by settig mean=rms=0
    if(nonzero> (AliTRDseedV1::kNtb*(1-lowfrac)) ){
      //only highest (1-lowfrac)*31 timebins are used to estimate the mean and rms! important! otherwise the 0' will make rms very large!
      mean = AliTRDdEdxBaseUtils::TruncatedMean(AliTRDseedV1::kNtb, vec->GetMatrixArray(), lowfrac, 1);
    }

    hmean->SetBinContent(ibin+1, mean);
  }

  delete obj;
}

void AliTRDdEdxCalibUtils::Output(const TList *lin, Int_t run)
{
  //
  //produce calibration objects
  //

  TString objnames;
  for(Int_t iter=0; iter<8; iter++){
    objnames+= AliTRDdEdxCalibHistArray::GetNameAt(iter)+" ";
  }

  TList * lout = new TList;
  lout->SetOwner();

  TTreeSRedirector *calibStream = new TTreeSRedirector(Form("TRDCalibStream_%010d.root", run));
    
  const Int_t nh=lin->GetEntries();
  for(Int_t ii=0; ii<nh; ii++){
    const THnBase *hh=(THnBase*)lin->At(ii);
    const TString hname = hh->GetName();
    if(!objnames.Contains(hname))
      continue;

    TObjArray * cobj0 = HistToObj(hh, run, lout, calibStream);
    lout->Add(cobj0);
  }

  //lout->ls();

  //=============================================================
  //=============================================================
  
  TFile *fout=new TFile(Form("TRDCalibObj_%010d.root", run),"recreate");
  fout->cd();
  const Int_t nout=lout->GetEntries();
  for(Int_t ii=0; ii<nout; ii++){
    const TString oname = lout->At(ii)->GetName();
    if(oname.Contains("Obj")){
      TObjArray * cobj = (TObjArray*) lout->At(ii);
      cobj->Write(oname, TObjArray::kSingleKey);
    }
  }
  fout->Save();
  fout->Close();
  delete fout;

  fout=new TFile(Form("TRDCalibList_%010d.root", run),"recreate");
  fout->cd();
  lin->Write();
  lout->Write();
  fout->Save();
  fout->Close();
  delete fout;
  
  delete calibStream;

  /*
    http://root.cern.ch/root/html/TH1.html
    When an histogram is created, a reference to it is automatically added to the list of in-memory objects for the current file or directory. This default behaviour can be changed by:
    
    h->SetDirectory(0);          for the current histogram h
    TH1::AddDirectory(kFALSE);   sets a global switch disabling the reference
    
    When the histogram is deleted, the reference to it is removed from the list of objects in memory. When a file is closed, all histograms in memory associated with this file are automatically deleted. 
  */
  delete lout;
}

TObjArray* AliTRDdEdxCalibUtils::HistToObj(const THnBase *hh, Int_t run, TList *lout, TTreeSRedirector *calibStream)
{
  //
  //produce calibration objects
  //

  const TString hname = hh->GetName();
  const Bool_t kinvq = TString(hname(hname.First('Q')+1,1)).Atoi()&4;

  //----------------------------------------
  //               Define nbin, tag, cobj0
  //----------------------------------------
  const Int_t nbin = AliTRDdEdxBaseUtils::NTRDtimebin();

    
  TString tag(hname);
  tag.ReplaceAll("Hist","Obj");
  
  TObjArray * cobj0 = new TObjArray(1);
  cobj0->SetOwner();
  cobj0->SetName(tag);
  cobj0->Add(new TVectorD(nbin));
  
  //----------------------------------------
  //               Define lowFrac, highFrac
  //----------------------------------------
  Double_t lowFrac = -999, highFrac = -999;
  if(!kinvq) {
    lowFrac = 0.01; highFrac = AliTRDdEdxBaseUtils::Q0Frac();
  }
  else{
    lowFrac = AliTRDdEdxBaseUtils::Q1Frac(); highFrac = 0.99;
  }
  
  //----------------------------------------
  //               Get analysis result
  //----------------------------------------
  TH1::AddDirectory(kFALSE);//important!
  TH1D *hnor=0x0, *hmpv=0x0, *hres=0x0, *hwid=0x0, *htrdphmean = 0x0;//if(!lout), these have to be deleted
  TH2D *hpj = hh->Projection(1,0);
  AliTRDdEdxBaseUtils::FitSlicesY(hpj, hnor, hmpv, hwid, hres, 0, lowFrac, highFrac);
  GetPHCountMeanRMS(hnor, htrdphmean);
  if(lout){
    lout->Add(htrdphmean);
  }
  delete hpj;
  
  if(lout){
    lout->Add(hnor);
    lout->Add(hmpv);
    lout->Add(hwid);
    lout->Add(hres);
  }
  
  //----------------------------------------
  //               Define Counter
  //----------------------------------------
  TVectorD *countDet=0x0;
  TObjArray *countSSL=0x0;

  if(!kinvq){
    countDet=new TVectorD(540);
    countSSL=new TObjArray(90);//SectorStackLayer
    countSSL->SetOwner();
    for(Int_t ii=0; ii<90; ii++){
      countSSL->Add(new TVectorD(6));
    }
  }

  //----------------------------------------
  //               Fill cobj0
  //----------------------------------------

  //ibin = binlowedge of bin(ibin+1) = the number fills this bin
  for(Int_t ibin=0; ibin<nbin; ibin++){
    Double_t gnor = hnor->GetBinContent(ibin+1);
    Double_t gmpv = hmpv->GetBinContent(ibin+1);
    Double_t gwid = hwid->GetBinContent(ibin+1);
    Double_t gres = hres->GetBinContent(ibin+1);

    //--- set additional cut by kpass
    Bool_t kpass = kTRUE;
    Double_t gtrdphmean = -999;
    if(htrdphmean){
      gtrdphmean = htrdphmean->GetBinContent(ibin+1);
      //chamber no statistics (e.g. too many 0's), not usual, not seen in run 143237
      if(gtrdphmean<EPSILON){
        kpass = kFALSE;
      }
      if(gnor<AliTRDdEdxBaseUtils::TimeBinCountCut()*gtrdphmean){
        kpass = kFALSE;
      }
    }
    
    //--- set calibration constant p0
    Double_t p0= 0;
    
    //reason for gmpv=0:
    //1)gnor<=3; truncation in hist: (0, 0.6*ntot=1.8 with ntot=3]={1}, in hist entries can pile up so that ntot=2, or 3, and (ntot>nstart && ntot<=nstop) is skipped;
    //2)TruncatedMean(hist) out of range (only for Q0, not Q1).
    
    if(gmpv>EPSILON && kpass){ 
      if(tag.Contains("T0")){
        p0 = gmpv;
      }
      else{
        p0 = 1/gmpv;
      }
      //printf("outcalibobj%s %d %15.6e\n", tag.Data(), ibin, p0);
    }

    (*( (TVectorD*)cobj0->At(0) ))[ibin] = p0;

    //--- save optional record
    if(p0>EPSILON && countDet && countSSL){
      const Int_t idet = AliTRDdEdxBaseUtils::ToDetector(ibin);
      (*countDet)[idet]=1;
      
      const Int_t isector = AliTRDdEdxBaseUtils::ToSector(ibin);
      const Int_t istack = AliTRDdEdxBaseUtils::ToStack(ibin);
      const Int_t ilayer = AliTRDdEdxBaseUtils::ToLayer(ibin);
      TVectorD * vecsectorstack = (TVectorD*)countSSL->At(istack*18+isector);
      (*vecsectorstack)[ilayer]=1;
    }
    
    if(calibStream){
      (*calibStream)<<tag<<
        "run="<<run<<
        "p0="<<p0<<
        
        "nor="<<gnor<<
        "mpv="<<gmpv<<
        "wid="<<gwid<<
        "res="<<gres<<
        "gtrdphmean="<<gtrdphmean<<
        
        "ibin="<<ibin<<
        "\n";
    }
  }
  
  //----------------------------------------
  //               Status Report
  //----------------------------------------
  if(countDet && countSSL){
    TVectorD count2Dstack(90);
    for(Int_t ii=0; ii<90; ii++){
      TVectorD * vecsectorstack = (TVectorD*)countSSL->At(ii);
      const Int_t nlayer = (Int_t)vecsectorstack->Sum();
      if(nlayer==6){
        count2Dstack[ii]=1;
      }
    }

    printf("\nAliTRDdEdxCalibUtils::GetCalibObj Summary run: %d name: %s entries: %.0f ndetector: %03.0f n2dstack %02.0f\n\n", run, hname.Data(), hh->GetEntries(), countDet->Sum(), count2Dstack.Sum());
  }

  //----------------------------------------
  //               Clean Up
  //----------------------------------------
  
  TH1D **hhs[]={&hnor, &hmpv, &hwid, &hres, &htrdphmean};
  const Int_t nhh=sizeof(hhs)/sizeof(TH1D**);
  for(Int_t ihh=0; ihh<nhh; ihh++){
    if(!lout){
      delete (*hhs[ihh]);
    }
  }
  
  delete countDet;
  delete countSSL;

  //----------------------------------------

  return cobj0;
}

