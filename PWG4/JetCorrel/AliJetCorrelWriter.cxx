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
// Contains also the methods for calculation of correlation parameters.
//-- Author: Paul Constantin

#include "AliJetCorrelWriter.h"

using namespace std;

ClassImp(AliJetCorrelWriter)
  
AliJetCorrelWriter::AliJetCorrelWriter() :
  fSelector(NULL), fMaker(NULL), fHname(""), fHtit(""), fRecoTrigg(kFALSE), fRndm(6251093), fNtuParent(NULL),
  fHBinsCentr(NULL), fHBinsZVert(NULL), fHBinsTrigg(NULL), fHBinsAssoc(NULL), fHCentr(NULL), fHZVert(NULL) {
  // constructor
  for(UInt_t i=0; i<2; i++){
    fHTrkITSQA[i] = NULL; fHTrkTPCQA[i] = NULL; fHTrkVTXQA[i] = NULL;
    for(UInt_t ic=0; ic<AliJetCorrelSelector::kMaxCent; ic++)  fHTrkProx[i][ic] = NULL;
  }
  for(UInt_t ik=0; ik<AliJetCorrelSelector::kMaxCorrel; ik++){
    fHTriggAcc[ik] = NULL; fHAssocAcc[ik] = NULL;
    for(UInt_t ic=0; ic<AliJetCorrelSelector::kMaxCent; ic++){
      fHTriggPt[ik][ic] = NULL; fHAssocPt[ik][ic] = NULL;
      for(UInt_t iv=0; iv<AliJetCorrelSelector::kMaxVert; iv++)
	for(UInt_t it=0; it<AliJetCorrelSelector::kMaxTrig; it++)
	  for(UInt_t ia=0; ia<AliJetCorrelSelector::kMaxAsso; ia++){
	    fHReal[ik][ic][iv][it][ia] = NULL;
	    fHMix[ik][ic][iv][it][ia]  = NULL;
	  }
    }
  }
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
  UInt_t nBinsCentr = fSelector->NoOfBins(t_cent);
  UInt_t nBinsZVert = fSelector->NoOfBins(t_vert);
  UInt_t nBinsTrigg = fSelector->NoOfBins(t_trig);
  UInt_t nBinsAssoc = fSelector->NoOfBins(t_asso);
  Float_t minTrigg  = fSelector->MinLowBin(t_trig);
  Float_t maxTrigg  = fSelector->MaxHighBin(t_trig);
  Float_t minAssoc  = fSelector->MinLowBin(t_asso);
  Float_t maxAssoc  = fSelector->MaxHighBin(t_asso);

  for(UInt_t ic=0; ic<nBinsCentr; ic++){ // loop over centrality bins
    for(UInt_t tt=0; tt<nTypeTrigg; tt++){ // loop over trigger types
      fHtit="type:"; fHtit+=tt; fHtit+=" cent:"; fHtit+=ic;
      fHname = "fHTriggPt"; fHname+=tt; fHname+=ic;
      fHTriggPt[tt][ic] = new TH1F(fHname,  fHtit, 10*nBinsTrigg, minTrigg, maxTrigg); // for <pT> each bin
      histosContainer->AddLast(fHTriggPt[tt][ic]);
    } // loop over trigger types
    for(UInt_t at=0; at<nTypeAssoc; at++){ // loop over trigger types
      fHtit="type:"; fHtit+=at; fHtit+=" cent:"; fHtit+=ic;
      fHname = "fHAssocPt"; fHname+=at; fHname+=ic;
      fHAssocPt[at][ic] = new TH1F(fHname,  fHtit, 10*nBinsAssoc, minAssoc, maxAssoc); // for <pT> each bin
      histosContainer->AddLast(fHAssocPt[at][ic]);
    } // loop over trigger types
  } // centrality binning loop

  fHCentr = new TH1F("fHCentr","centrality distribution",100, 0., 200.);
  histosContainer->AddLast(fHCentr);
  fHZVert = new TH1F("fHZVert","vertex distribution",100, -10., 10.);
  histosContainer->AddLast(fHZVert);

  fHBinsCentr = new TH1F("fHBinsCentr","centrality binning", nBinsCentr+1, 0, 1); 
  histosContainer->AddLast(fHBinsCentr);
  for(UInt_t i=1;i<=nBinsCentr+1; i++)
    fHBinsCentr->SetBinContent(i,fSelector->BinBorder(t_cent,i-1));
  fHBinsZVert = new TH1F("fHBinsZVert","centrality binning", nBinsZVert+1, 0, 1); 
  histosContainer->AddLast(fHBinsZVert);
  for(UInt_t i=1;i<=nBinsZVert+1; i++)
    fHBinsZVert->SetBinContent(i,fSelector->BinBorder(t_vert,i-1));
  fHBinsTrigg = new TH1F("fHBinsTrigg","trigger binning", nBinsTrigg+1, 0, 1); 
  histosContainer->AddLast(fHBinsTrigg);
  for(UInt_t i=1;i<=nBinsTrigg+1; i++)
    fHBinsTrigg->SetBinContent(i,fSelector->BinBorder(t_trig,i-1));
  fHBinsAssoc = new TH1F("fHBinsAssoc","associated binning", nBinsAssoc+1, 0, 1); 
  histosContainer->AddLast(fHBinsAssoc);
  for(UInt_t i=1;i<=nBinsAssoc+1; i++)
    fHBinsAssoc->SetBinContent(i,fSelector->BinBorder(t_asso,i-1));
}

void AliJetCorrelWriter::CreateQA(TList *histosContainer){
  // books QA histograms
  TString when[2] = {"before cuts","after cuts"};
  for(UInt_t i=0; i<2; i++){
    fHname = "fHTrkITSQA"; fHname += i;
    fHtit = "ITS nClust vs Chi2/nClust "; fHtit += when[i];
    fHTrkITSQA[i] = new TH2F(fHname,fHtit,20,0.,20.,50,0.,10.);
    histosContainer->AddLast(fHTrkITSQA[i]);
    fHname = "fHTrkTPCQA"; fHname += i;
    fHtit = "TPC nClust vs Chi2/nClust "; fHtit += when[i];
    fHTrkTPCQA[i] = new TH2F(fHname,fHtit,30,0.,150.,50,0.,10.);
    histosContainer->AddLast(fHTrkTPCQA[i]);
    fHname = "fHTrkVTXQA"; fHname += i;
    fHtit = "VTX KinkIndex vs nSigma "; fHtit += when[i];
    fHTrkVTXQA[i] = new TH2F(fHname,fHtit,21,-10.,10,50,0.,10.);
    histosContainer->AddLast(fHTrkVTXQA[i]);
  }

  UInt_t nTypeTrigg = fMaker->NoOfTrigg();
  UInt_t nTypeAssoc = fMaker->NoOfAssoc();
  UInt_t nBinsCentr = fSelector->NoOfBins(t_cent);
  UInt_t nBinsTrigg = fSelector->NoOfBins(t_trig);
  UInt_t nBinsAssoc = fSelector->NoOfBins(t_asso);
  Float_t minTrigg  = fSelector->MinLowBin(t_trig);
  Float_t maxTrigg  = fSelector->MaxHighBin(t_trig);
  Float_t minAssoc  = fSelector->MinLowBin(t_asso);
  Float_t maxAssoc  = fSelector->MaxHighBin(t_asso);
  for(UInt_t ic=0; ic<nBinsCentr; ic++){ // centrality loop
    for(UInt_t it=0; it<2; it++){ // mixing type loop (real/mixed)
      fHname="fHTrkProx"; fHname+=it; fHname+=ic;
      fHtit="fHTrkProx"; fHtit+=it; fHtit+=ic;
      fHTrkProx[it][ic] = new TH3F(fHname,fHtit,100,0.,10.,
				  nBinsTrigg,minTrigg,maxTrigg,nBinsAssoc,minAssoc,maxAssoc);
      histosContainer->AddLast(fHTrkProx[it][ic]);
    } // loop over mixing type
  } // loop over centrality
  for(UInt_t tt=0; tt<nTypeTrigg; tt++){ // loop over trigger types
    fHname = "fHTriggAcc"; fHname+=tt;
    fHTriggAcc[tt] = new TH3F(fHname,fHname,62,0.,2*TMath::Pi(),20,-1.,1.,nBinsTrigg,minTrigg,maxTrigg);
    histosContainer->AddLast(fHTriggAcc[tt]);
  } // loop over trigger types
  for(UInt_t ta=0; ta<nTypeAssoc; ta++){ // loop over associated types
    fHname="fHAssocAcc"; fHname+=ta;
    fHAssocAcc[ta] = new TH3F(fHname,fHname,62,0.,2*TMath::Pi(),20,-1.,1.,nBinsAssoc,minAssoc,maxAssoc);
    histosContainer->AddLast(fHAssocAcc[ta]);
  } // loop over associated types
}

void AliJetCorrelWriter::CreateCorrelations(TList* histosContainer){
  // books correlation histograms
  UInt_t nTypeCorrel = fMaker->NoOfCorrel();
  UInt_t nBinsCentr = fSelector->NoOfBins(t_cent);
  UInt_t nBinsZVert = fSelector->NoOfBins(t_vert);
  UInt_t nBinsTrigg = fSelector->NoOfBins(t_trig);
  UInt_t nBinsAssoc = fSelector->NoOfBins(t_asso);
  Float_t maxAssoc  = fSelector->MaxHighBin(t_asso);
  UInt_t kDPhiNumBins = fSelector->DPhiNumBins();
  UInt_t kDEtaNumBins = fSelector->DEtaNumBins();
  UInt_t nPoutBins = UInt_t(TMath::Ceil(maxAssoc/fSelector->PoutBW())); // since |p_out|<p_Ta
  if(fRecoTrigg) {  // if any correlation has reconstructed trigger, define ntuple; use id to differentiate
    fNtuParent = new TNtuple("fNtuParent","Reconstructed Parent Ntuple","id:q:m:pT:phi:eta:assym:open");
    histosContainer->AddLast(fNtuParent);
  }
  for(UInt_t htc=0; htc<nTypeCorrel; htc++){ // loop over correlation types
    for(UInt_t hic=0; hic<nBinsCentr; hic++){ // centrality loop
      for(UInt_t hiv=0; hiv<nBinsZVert; hiv++){ // vertex loop
	for(UInt_t hit=0; hit<nBinsTrigg; hit++){ // trigger loop
	  for(UInt_t hia=0; hia<nBinsAssoc; hia++){ // associated loop
	    fHtit = fMaker->Descriptor(htc); fHtit+=":"; fHtit+=hic; fHtit+=hiv; fHtit+=hit; fHtit+=hia;
	    fHname="fHReal"; fHname+=htc; fHname+=hic; fHname+=hiv; fHname+=hit; fHname+=hia;
	    fHReal[htc][hic][hiv][hit][hia] = new TH3F(fHname,fHtit, kDPhiNumBins,-1./3.,5./3.,
						      kDEtaNumBins,-1.6,1.6, nPoutBins,0.,maxAssoc);
	    histosContainer->AddLast(fHReal[htc][hic][hiv][hit][hia]);
	    fHname="fHMix"; fHname+=htc; fHname+=hic; fHname+=hiv; fHname+=hit; fHname+=hia;
	    fHMix[htc][hic][hiv][hit][hia] = new TH3F(fHname,fHtit, kDPhiNumBins,-1./3.,5./3.,
						      kDEtaNumBins,-1.6,1.6, nPoutBins,0.,maxAssoc);
	    histosContainer->AddLast(fHMix[htc][hic][hiv][hit][hia]);
	  } // loop over associated bins
	} // loop over trigger bins
      } // loop over vertex bins
    } // loop over centrality bins
  } // loop over correlation types
}

/////////////////////////////////////////////////////////
// METHODS FOR FILLING THE OUTPUT HISTOGRAMS
/////////////////////////////////////////////////////////

void AliJetCorrelWriter::FillGlobal(Float_t cent, Float_t vert){
  // some global event histos
  fHCentr->Fill(cent);
  fHZVert->Fill(vert);
}

void AliJetCorrelWriter::FillSingleHistos(UInt_t cBin, CorrelList_t * const TriggList, UInt_t tIdx,
					  CorrelList_t * const AssocList, UInt_t aIdx){
  // fills single-particle histograms
  if(TriggList->Size()>0){
    CorrelListIter_t trigIter = TriggList->Head();
    while(!trigIter.HasEnded()){
      CorrelParticle_t *particle = trigIter.Data();
      fHTriggPt[tIdx][cBin]->Fill(particle->Pt());
      if(fSelector->GenQA()) fHTriggAcc[tIdx]->Fill(particle->Phi(),particle->Eta(),particle->Pt());
      trigIter.Move();
    }
  }
  if(AssocList->Size()>0){
    CorrelListIter_t assoIter = AssocList->Head();
    while(!assoIter.HasEnded()){
      CorrelParticle_t *particle = assoIter.Data();
      fHAssocPt[aIdx][cBin]->Fill(particle->Pt());
      if(fSelector->GenQA()) fHAssocAcc[aIdx]->Fill(particle->Phi(),particle->Eta(),particle->Pt());
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

  fHTrkITSQA[idx]->Fill(Float_t(nClusITS),chi2ITS);
  fHTrkTPCQA[idx]->Fill(Float_t(nClusTPC),chi2TPC);
  fHTrkVTXQA[idx]->Fill(Float_t(kinkIndex),sigVtx);
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
    fNtuParent->Fill(parVar);
    parIter.Move();
  }
}

void AliJetCorrelWriter::FillCorrelations(UInt_t fTyp, UInt_t iCorr, UInt_t cBin, UInt_t vBin, 
					  CorrelParticle_t * const Trigg, CorrelParticle_t * const Assoc){
  // fills the correlation (two-particle) histograms
  // trigger information (this is why the first particle has to be the trigger):
  Float_t ptt  = Trigg->Pt();
  Float_t phit = Trigg->Phi();
  Float_t etat = Trigg->Eta();
  Int_t   tBin = fSelector->GetBin(t_trig,ptt);
  // associated information:
  Float_t pta  = Assoc->Pt();
  Float_t phia = Assoc->Phi();
  Float_t etaa = Assoc->Eta();
  Int_t   aBin = fSelector->GetBin(t_asso,pta);
  //  Short_t qprod= Trigg->Q()*Assoc->Q();

  if(tBin<0 || aBin<0) return;  // one of them is not in the required pT range
  if(pta>=ptt) return; // use only associated particles below the trigger
  if(fabs(ptt-pta)<1.e-6 && fabs(phit-phia)<1.e-6 && fabs(etat-etaa)<1.e-6) return; // don't auto-correlate

  // track pair proximity 
  // Before uncommenting, first store CorrelTrack_t in lists see AliJetCorrelReader::FillESDTrackLists
//   if(Trigg->ID()==hadron && Assoc->ID()==hadron){
//     CorrelTrack_t* trk1 = dynamic_cast<CorrelTrack_t*>(Trigg);
//     CorrelTrack_t* trk2 = dynamic_cast<CorrelTrack_t*>(Assoc);
//     if(!trk1 || !trk2)
//       {std::cerr<<"AliJetCorrelWriter::FillCorrelations: failed casting!"<<std::endl; exit(-1);}
//     Float_t pairDist = trk1->Dist(trk2);
//     if(fSelector->CloseTrackPair(pairDist)) return; // proximity cut
//     if(fSelector->GenQA()) fHTrkProx[fTyp][cBin]->Fill(pairDist,ptt,pta);
//   }

  // Fill correlation histograms:
  Float_t dphi = DeltaPhi(phit,phia);
  Float_t deta = etat-etaa;
  Float_t pout = TMath::Abs(pta*TMath::Sin(dphi));
  if(fTyp==0)
    fHReal[iCorr][cBin][vBin][tBin][aBin]->Fill(dphi/TMath::Pi(),deta,pout);
  else
    fHMix[iCorr][cBin][vBin][tBin][aBin]->Fill(dphi/TMath::Pi(),deta,pout);
}

Float_t AliJetCorrelWriter::DeltaPhi(Float_t  phi1, Float_t phi2) {
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

