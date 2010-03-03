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

//______________________________________________________
// Class for output (histograms) definition and filling.
//-- Author: Paul Constantin

#include "AliJetCorrelWriter.h"

using namespace std;

ClassImp(AliJetCorrelWriter)

AliJetCorrelWriter::AliJetCorrelWriter() :
  fSelector(NULL), fMaker(NULL), hname(""), htit(""), fRecoTrigg(kFALSE), fRndm(6251093),
  hBinsCentr(NULL), hBinsZVert(NULL), hCentr(NULL), hZVert(NULL), ntuParent(NULL) {
  // constructor
  for(UInt_t i=0; i<kMAXNUMCORREL; i++) 
    for(UInt_t j=0; j<kMAXCENTBIN; j++) 
      {fNumReal[i][j]=0; fNumMix[i][j]=0;}
}

AliJetCorrelWriter::~AliJetCorrelWriter(){
  // destructor
}

void AliJetCorrelWriter::Init(AliJetCorrelSelector * const s, AliJetCorrelMaker * const m){
  // initialization method
  fSelector = s;
  fMaker = m;
  fRecoTrigg = m->RecoTrigger();
}

////////////////////////////////////////////
// METHODS FOR CREATION OF OUTPUT HISTOGRAMS
////////////////////////////////////////////

void AliJetCorrelWriter::CreateGeneric(TList *histosContainer){
  // books generic histograms
  UInt_t nTypeTrigg = fMaker->NoOfTrigg();
  UInt_t nTypeAssoc = fMaker->NoOfAssoc();
  UInt_t nBinsCentr = fSelector->NoOfBins(centr);
  UInt_t nBinsZVert = fSelector->NoOfBins(zvert);
  UInt_t nBinsTrigg = fSelector->NumTriggPt();
  Float_t minTrigg  = fSelector->MinTriggPt();
  Float_t maxTrigg  = fSelector->MaxTriggPt();
  UInt_t nBinsAssoc = fSelector->NumAssocPt();
  Float_t minAssoc  = fSelector->MinAssocPt();
  Float_t maxAssoc  = fSelector->MaxAssocPt();

  for(UInt_t ic=0; ic<nBinsCentr; ic++){ // loop over centrality bins
    for(UInt_t tt=0; tt<nTypeTrigg; tt++){ // loop over trigger types
      htit="type:"; htit+=tt; htit+=" cent:"; htit+=ic;
      hname = "hTriggPt"; hname+=tt; hname+=ic;
      hTriggPt[tt][ic] = new TH1F(hname,  htit, 10*nBinsTrigg, minTrigg, maxTrigg); // for <pT> each bin
      histosContainer->AddLast(hTriggPt[tt][ic]);
    } // loop over trigger types
    for(UInt_t at=0; at<nTypeAssoc; at++){ // loop over trigger types
      htit="type:"; htit+=at; htit+=" cent:"; htit+=ic;
      hname = "hAssocPt"; hname+=at; hname+=ic;
      hAssocPt[at][ic] = new TH1F(hname,  htit, 10*nBinsAssoc, minAssoc, maxAssoc); // for <pT> each bin
      histosContainer->AddLast(hAssocPt[at][ic]);
    } // loop over trigger types
  } // centrality binning loop

  hCentr = new TH1F("hCentr","centrality distribution",100, 0., 100.);
  histosContainer->AddLast(hCentr);
  hZVert = new TH1F("hZVert","vertex distribution",100, -50., 50.);
  histosContainer->AddLast(hZVert);

  hBinsCentr = new TH1F("hBinsCentr","centrality binning", nBinsCentr+1, 0, 1); 
  histosContainer->AddLast(hBinsCentr);
  for(UInt_t i=1;i<=nBinsCentr+1; i++)
    hBinsCentr->SetBinContent(i,fSelector->BinBorder(centr,i-1));
  hBinsZVert = new TH1F("hBinsZVert","centrality binning", nBinsZVert+1, 0, 1); 
  histosContainer->AddLast(hBinsZVert);
  for(UInt_t i=1;i<=nBinsZVert+1; i++)
    hBinsZVert->SetBinContent(i,fSelector->BinBorder(zvert,i-1));
}

void AliJetCorrelWriter::CreateQA(TList *histosContainer){
  // books QA histograms
  TString when[2] = {"before cuts","after cuts"};
  for(UInt_t i=0; i<2; i++){
    hname = "hTrkITSQA"; hname += i;
    htit = "ITS nClust vs Chi2/nClust "; htit += when[i];
    hTrkITSQA[i] = new TH2F(hname,htit,20,0.,20.,50,0.,10.);
    histosContainer->AddLast(hTrkITSQA[i]);
    hname = "hTrkTPCQA"; hname += i;
    htit = "TPC nClust vs Chi2/nClust "; htit += when[i];
    hTrkTPCQA[i] = new TH2F(hname,htit,30,0.,150.,50,0.,10.);
    histosContainer->AddLast(hTrkTPCQA[i]);
    hname = "hTrkVTXQA"; hname += i;
    htit = "VTX KinkIndex vs nSigma "; htit += when[i];
    hTrkVTXQA[i] = new TH2F(hname,htit,21,-10.,10,50,0.,10.);
    histosContainer->AddLast(hTrkVTXQA[i]);
  }

  UInt_t nTypeTrigg = fMaker->NoOfTrigg();
  UInt_t nTypeAssoc = fMaker->NoOfAssoc();
  UInt_t nBinsCentr = fSelector->NoOfBins(centr);
  UInt_t nBinsTrigg = fSelector->NumTriggPt();
  Float_t minTrigg  = fSelector->MinTriggPt();
  Float_t maxTrigg  = fSelector->MaxTriggPt();
  UInt_t nBinsAssoc = fSelector->NumAssocPt();
  Float_t minAssoc  = fSelector->MinAssocPt();
  Float_t maxAssoc  = fSelector->MaxAssocPt();
  for(UInt_t ic=0; ic<nBinsCentr; ic++){ // centrality loop
    for(UInt_t it=0; it<2; it++){ // mixing type loop (real/mixed)
      hname="hTrkProx"; hname+=it; hname+=ic;
      htit="hTrkProx"; htit+=it; htit+=ic;
      hTrkProx[it][ic] = new TH3F(hname,htit,100,0.,10.,
				  nBinsTrigg,minTrigg,maxTrigg,nBinsAssoc,minAssoc,maxAssoc);
      histosContainer->AddLast(hTrkProx[it][ic]);
    } // loop over mixing type
  } // loop over centrality
  for(UInt_t tt=0; tt<nTypeTrigg; tt++){ // loop over trigger types
    hname = "hTriggAcc"; hname+=tt;
    hTriggAcc[tt] = new TH3F(hname,hname,62,0.,2*TMath::Pi(),20,-1.,1.,nBinsTrigg,minTrigg,maxTrigg);
    histosContainer->AddLast(hTriggAcc[tt]);
  } // loop over trigger types
  for(UInt_t ta=0; ta<nTypeAssoc; ta++){ // loop over associated types
    hname="hAssocAcc"; hname+=ta;
    hAssocAcc[ta] = new TH3F(hname,hname,62,0.,2*TMath::Pi(),20,-1.,1.,nBinsAssoc,minAssoc,maxAssoc);
    histosContainer->AddLast(hAssocAcc[ta]);
  } // loop over associated types
}

void AliJetCorrelWriter::CreateCorrelations(TList* histosContainer){
  // books correlation histograms
  const Float_t lr=-TMath::Pi()/3, ur=5*TMath::Pi()/3, bwPout=0.2;
  UInt_t nTypeCorrel = fMaker->NoOfCorrel();
  UInt_t nBinsCentr = fSelector->NoOfBins(centr);
  UInt_t nBinsZVert = fSelector->NoOfBins(zvert);
  UInt_t nBinsTrigg = fSelector->NumTriggPt();
  UInt_t nBinsAssoc = fSelector->NumAssocPt();
  Float_t maxAssoc  = fSelector->MaxAssocPt();
  UInt_t kDPhiNumBins = fSelector->DPhiNumBins();
  UInt_t kDEtaNumBins = fSelector->DEtaNumBins();
  UInt_t nPoutBins = UInt_t(TMath::Ceil(2.2*maxAssoc/bwPout)); // since |p_out|<p_Ta	    

  if(fRecoTrigg) {  // if any correlation has reconstructed trigger, define ntuple; use id to differentiate
    ntuParent = new TNtuple("ntuParent","Reconstructed Parent Ntuple","id:q:m:pT:phi:eta:assym:open");
    histosContainer->AddLast(ntuParent);
  }
  for(UInt_t htc=0; htc<nTypeCorrel; htc++){ // loop over correlation types
    for(UInt_t hic=0; hic<nBinsCentr; hic++){ // centrality loop
      for(UInt_t hiv=0; hiv<nBinsZVert; hiv++){ // vertex loop
	for(UInt_t hit=0; hit<nBinsTrigg; hit++){ // trigger loop
	  for(UInt_t hia=0; hia<nBinsAssoc; hia++){ // associated loop
	    htit = fMaker->Descriptor(htc); htit+=":"; htit+=hic; htit+=hiv; htit+=hit; htit+=hia;
	    hname="hReal"; hname+=htc; hname+=hic; hname+=hiv; hname+=hit; hname+=hia;
	    hReal[htc][hic][hiv][hit][hia] = new TH3F(hname,htit, kDPhiNumBins,lr,ur,
						      kDEtaNumBins,-2.,2.,
						      nPoutBins,-1.1*maxAssoc,1.1*maxAssoc);
	    histosContainer->AddLast(hReal[htc][hic][hiv][hit][hia]);
	    hname="hMix"; hname+=htc; hname+=hic; hname+=hiv; hname+=hit; hname+=hia;
	    hMix[htc][hic][hiv][hit][hia] = new TH3F(hname,htit, kDPhiNumBins,lr,ur,
						      kDEtaNumBins,-2.,2.,
						      nPoutBins,-1.1*maxAssoc,1.1*maxAssoc);
	    histosContainer->AddLast(hMix[htc][hic][hiv][hit][hia]);
	  } // loop over associated bins
	} // loop over trigger bins
      } // loop over vertex bins
    } // loop over centrality bins
  } // loop over correlation types
}

/////////////////////////////////////////////////////////
// METHODS FOR FILLING THE OUTPUT HISTOGRAMS
/////////////////////////////////////////////////////////

void AliJetCorrelWriter::FillGlobal(Float_t cent, Float_t zvert){
  // some global event histos
  hCentr->Fill(cent);
  hZVert->Fill(zvert);
}

void AliJetCorrelWriter::FillSingleHistos(UInt_t cBin, CorrelList_t * const TriggList, UInt_t tIdx,
					  CorrelList_t * const AssocList, UInt_t aIdx){
  // fills single-particle histograms
  if(TriggList->Size()>0){
    CorrelListIter_t trigIter = TriggList->Head();
    while(!trigIter.HasEnded()){
      CorrelParticle_t *particle = trigIter.Data();
      hTriggPt[tIdx][cBin]->Fill(particle->Pt());
      if(fSelector->GenQA()) hTriggAcc[tIdx]->Fill(particle->Phi(),particle->Eta(),particle->Pt());
      trigIter.Move();
    }
  }
  if(AssocList->Size()>0){
    CorrelListIter_t assoIter = AssocList->Head();
    while(!assoIter.HasEnded()){
      CorrelParticle_t *particle = assoIter.Data();
      hAssocPt[aIdx][cBin]->Fill(particle->Pt());
      if(fSelector->GenQA()) hAssocAcc[aIdx]->Fill(particle->Phi(),particle->Eta(),particle->Pt());
      assoIter.Move();
    }
  }
}

void AliJetCorrelWriter::FillTrackQA(AliESDtrack * const track, UInt_t idx){
  // fills single-particle QA
  if(idx>1){std::cerr<<"AliJetCorrelWriter::FillTrackQA: wrong idx!"<<std::endl; exit(-1);}

  UInt_t nClusITS = track->GetITSclusters(0);
  UInt_t nClusTPC = track->GetTPCclusters(0); // or track->GetTPCNcls() ?
  Float_t chi2ITS=-1., chi2TPC=-1.;
  if(nClusITS!=0) chi2ITS = track->GetITSchi2()/Float_t(nClusITS);
  if(nClusTPC!=0) chi2TPC = track->GetTPCchi2()/Float_t(nClusTPC);
  UInt_t kinkIndex = track->GetKinkIndex(0);
  Float_t sigVtx = fSelector->GetSigmaToVertex(track);

  hTrkITSQA[idx]->Fill(Float_t(nClusITS),chi2ITS);
  hTrkTPCQA[idx]->Fill(Float_t(nClusTPC),chi2TPC);
  hTrkVTXQA[idx]->Fill(Float_t(kinkIndex),sigVtx);
}

void AliJetCorrelWriter::FillParentNtuple(CorrelList_t * const ParentList){
  // fills ntuple of triggers when they are reconstructed parents (pi0,Z0,etc.)
  if(!fRecoTrigg)
    {std::cerr<<"AliJetCorrelWriter::FillParentNtuple: you shouldn't be here!"<<std::endl; exit(-1);}
  if(ParentList->Size()<1) return;
  Float_t parVar[8]={0.};
  CorrelListIter_t parIter = ParentList->Head();
  while(!parIter.HasEnded()){
    CorrelRecoParent_t *parent = dynamic_cast<CorrelRecoParent_t*>(parIter.Data());
    if(!parent)
      {std::cerr<<"AliJetCorrelWriter::FillParentNtuple: failed casting!"<<std::endl; exit(-1);}
    parVar[0] = parent->ID();
    parVar[1] = parent->Q();
    parVar[2] = parent->M();
    parVar[3] = parent->Pt();
    parVar[4] = parent->Phi();
    parVar[5] = parent->Eta();
    parVar[6] = parent->Assym();
    parVar[7] = parent->OpenAng();
    ntuParent->Fill(parVar);
    parIter.Move();
  }
}

void AliJetCorrelWriter::FillCorrelations(FillType_t fTyp, UInt_t iCorr, UInt_t cBin, UInt_t vBin, 
					  CorrelParticle_t * const Trigg, CorrelParticle_t * const Assoc){
  // fills the correlation (two-particle) histograms
  // trigger information (this is why the first particle has to be the trigger):
  Float_t ptt  = Trigg->Pt();
  Float_t phit = Trigg->Phi();
  Float_t etat = Trigg->Eta();
  Int_t   tBin = fSelector->TriggBin(ptt);
  // associated information:
  Float_t pta  = Assoc->Pt();
  Float_t phia = Assoc->Phi();
  Float_t etaa = Assoc->Eta();
  Int_t   aBin = fSelector->AssocBin(pta);

  if(tBin<0 || aBin<0) return;  // one of them is not in the required pT range
  if(fabs(ptt-pta)<1.e-6 && fabs(phit-phia)<1.e-6 && fabs(etat-etaa)<1.e-6) return; // don't auto-correlate

  // store track pair proximity
  if(Trigg->ID()==hadron && Assoc->ID()==hadron){
    CorrelTrack_t* trk1 = dynamic_cast<CorrelTrack_t*>(Trigg);
    CorrelTrack_t* trk2 = dynamic_cast<CorrelTrack_t*>(Assoc);
    if(!trk1 || !trk2)
      {std::cerr<<"AliJetCorrelWriter::FillCorrelations: failed casting!"<<std::endl; exit(-1);}
    Float_t pairDist = trk1->Dist(trk2);
    if(fSelector->CloseTrackPair(pairDist)) return; // proximity cut
    if(fSelector->GenQA()) hTrkProx[fTyp][cBin]->Fill(pairDist,ptt,pta);
  }

  // Fill correlation histograms:
  Float_t dphi = DeltaPhi(phit,phia);
  Float_t deta = etat-etaa;
  Float_t pout = pta*TMath::Sin(dphi);
  if(fTyp==real){
    hReal[iCorr][cBin][vBin][tBin][aBin]->Fill(dphi,deta,pout);
    fNumReal[iCorr][cBin]++;
  } else {
    hMix[iCorr][cBin][vBin][tBin][aBin]->Fill(dphi,deta,pout);
    fNumMix[iCorr][cBin]++;
  }
}

void AliJetCorrelWriter::ShowStats(){
  // stats printout method
  UInt_t nTypeCorrel = fMaker->NoOfCorrel();
  UInt_t nBinsCentr = fSelector->NoOfBins(centr);
  for(UInt_t i=0; i<nTypeCorrel; i++){
    std::cout<<"Correlation:"<<fMaker->Descriptor(i)<<std::endl;
    for(UInt_t j=0; j<nBinsCentr; j++){
      Float_t cb1 = fSelector->BinBorder(centr, j);
      Float_t cb2 = fSelector->BinBorder(centr, j+1);
      std::cout<<" Centrality:"<<cb1<<"-"<<cb2
	       <<"\t Real Pairs="<<fNumReal[i][j]<<"\t Mixed Pairs="<<fNumMix[i][j]<<std::endl;
    }
  }
}

Float_t AliJetCorrelWriter::DeltaPhi(Float_t  phi1, Float_t phi2){
// wrapps dphi in (-pi/3,5pi/3)
  Float_t kPi = TMath::Pi();
  Float_t res = phi2-phi1;
  if(fRndm.Uniform()<0.5) res = phi1-phi2;
  if(5.*kPi/3.<res && res<=2.*kPi) res-=2*kPi;
  else if(-2.*kPi<=res && res<-kPi/3.) res+=2*kPi;
//   if(3.*kPi/2.<res && res<=2*kPi) res-=2*kPi;
//   else if(-2.*kPi<=res && res<-kPi/2.) res+=2*kPi;
  return res;  
}

