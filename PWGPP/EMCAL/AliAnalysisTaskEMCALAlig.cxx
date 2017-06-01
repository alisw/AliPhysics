/**************************************************************************
 * Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
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

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskEMCALAlig.h"
#include "AliInputEventHandler.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliExternalTrackParam.h"
#include "TVector3.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEMCALAlig)
/// \endcond

AliAnalysisTaskEMCALAlig::AliAnalysisTaskEMCALAlig() :
AliAnalysisTaskEmcal(),
fEMCALRecoUtils(NULL),
fEMCALGeo(NULL),
fPIDResponse(NULL),
fNPartPt(NULL),
fNPartPhi(NULL),
fNPartEta(NULL),
fPPartPt(NULL),
fPPartPhi(NULL),
fPPartEta(NULL),
fTPCnSgima(NULL),
fEOverP(NULL),
fElectronPhiRes(NULL),
fElectronEtaRes(NULL),
fElectronXRes(NULL),
fElectronYRes(NULL),
fElectronPositronZRes(NULL),
fPositronPhiRes(NULL),
fPositronEtaRes(NULL),
fPositronXRes(NULL),
fPositronYRes(NULL),
fAllMatchedTracksPhiRes(NULL),
fAllMatchedTracksEtaRes(NULL),
fElectronPhiResRU(NULL),
fElectronEtaResRU(NULL),
fElectronXResRU(NULL),
fElectronYResRU(NULL),
fElectronPositronZResRU(NULL),
fPositronPhiResRU(NULL),
fPositronEtaResRU(NULL),
fPositronXResRU(NULL),
fPositronYResRU(NULL)
{
    
}

AliAnalysisTaskEMCALAlig::AliAnalysisTaskEMCALAlig(const char *name) :
AliAnalysisTaskEmcal(name, kTRUE),
fEMCALRecoUtils(NULL),
fEMCALGeo(NULL),
fPIDResponse(NULL),
fNPartPt(NULL),
fNPartPhi(NULL),
fNPartEta(NULL),
fPPartPt(NULL),
fPPartPhi(NULL),
fPPartEta(NULL),
fTPCnSgima(NULL),
fEOverP(NULL),
fElectronPhiRes(NULL),
fElectronEtaRes(NULL),
fElectronXRes(NULL),
fElectronYRes(NULL),
fElectronPositronZRes(NULL),
fPositronPhiRes(NULL),
fPositronEtaRes(NULL),
fPositronXRes(NULL),
fPositronYRes(NULL),
fAllMatchedTracksPhiRes(NULL),
fAllMatchedTracksEtaRes(NULL),
fElectronPhiResRU(NULL),
fElectronEtaResRU(NULL),
fElectronXResRU(NULL),
fElectronYResRU(NULL),
fElectronPositronZResRU(NULL),
fPositronPhiResRU(NULL),
fPositronEtaResRU(NULL),
fPositronXResRU(NULL),
fPositronYResRU(NULL)

{
    SetMakeGeneralHistograms(kTRUE);
}


AliAnalysisTaskEMCALAlig::~AliAnalysisTaskEMCALAlig()
{
    
}


void AliAnalysisTaskEMCALAlig::UserCreateOutputObjects()
{
    AliAnalysisTaskEmcal::UserCreateOutputObjects();
    
    fPIDResponse = fInputHandler->GetPIDResponse();
    
    //Negative Particle Histograms;
    fNPartPt = new TH1F("fNPartPt","p_{T} distribution of Electrons ;p_{T} (GeV/c);counts",200,0,20);
    fOutput->Add(fNPartPt);
    fNPartPhi = new TH1F("fNPartPhi","Electron #phi distribution;#phi;counts",100,0,6.3);
    fOutput->Add(fNPartPhi);
    fNPartEta = new TH1F("fNPartEta","Electron #eta distribution;#eta;counts",100,-1.0,1.0);
    fOutput->Add(fNPartEta);
    
    //Positive Particle Histograms;
    fPPartPt = new TH1F("fPPartPt","p_{T} distribution of Positrons;p_{T} (GeV/c);counts",200,0,20);
    fOutput->Add(fPPartPt);
    fPPartPhi = new TH1F("fPPartPhi","Positron #phi distribution;#phi;counts",100,0,6.3);
    fOutput->Add(fPPartPhi);
    fPPartEta = new TH1F("fPPartEta","Positron #eta distribution;#eta;counts",100,-1.0,1.0);
    fOutput->Add(fPPartEta);
    
    //PID Plots
    fTPCnSgima = new TH2F("fTPCnSgima","All Track TPC Nsigma distribution;p (GeV/c);#sigma_{TPC-dE/dx}",300,0,30,200,-20,20); // TPC Nsigma as function of pT
    fOutput->Add(fTPCnSgima);
    
    fEOverP = new TH2F("fEOverP", "E/p distribution;p_{T} (GeV/c);E/p", 200,0,20,60, 0.0, 3.0); // E/p as function of pT
    fOutput->Add(fEOverP);
    
    
    fElectronPhiRes = new TH2F *[21];
    fElectronEtaRes = new TH2F *[21];
    fElectronXRes = new TH2F *[21];
    fElectronYRes = new TH2F *[21];
    fElectronPositronZRes = new TH2F *[21];
    fPositronPhiRes = new TH2F *[21];
    fPositronEtaRes = new TH2F *[21];
    fPositronXRes = new TH2F *[21];
    fPositronYRes = new TH2F *[21];
    fAllMatchedTracksPhiRes = new TH2F *[21];
    fAllMatchedTracksEtaRes = new TH2F *[21];
    
    fElectronPhiResRU = new TH2F *[21];
    fElectronEtaResRU = new TH2F *[21];
    fElectronXResRU = new TH2F *[21];
    fElectronYResRU = new TH2F *[21];
    fElectronPositronZResRU = new TH2F *[21];
    fPositronPhiResRU = new TH2F *[21];
    fPositronEtaResRU = new TH2F *[21];
    fPositronXResRU = new TH2F *[21];
    fPositronYResRU = new TH2F *[21];
    
    
    for(Int_t i = 0; i<21;i++)
    {
        fElectronPhiRes[i] = new TH2F(Form("fElectronPhiRes%d",i),"Distance (phi) of EMCAL cluster in ID electrons;#phi;#Delta#phi",100,0,2*TMath::Pi(),1000,-0.3,0.3);
        fOutput->Add(fElectronPhiRes[i]);
        
        fElectronEtaRes[i] = new TH2F(Form("fElectronEtaRes%d",i),"Distance (eta) of EMCAL clusterin ID electrons;#phi;#Delta#eta",100,0,2*TMath::Pi(),1000,-0.3,0.3);
        fOutput->Add(fElectronEtaRes[i]);
        
        fElectronXRes[i] = new TH2F(Form("fElectronXRes%d",i),"Distance (x) of EMCAL cluster in ID electrons;#phi;#Delta x",100,0,2*TMath::Pi(),1000,-20,20);
        fOutput->Add(fElectronXRes[i]);
        
        fElectronYRes[i] = new TH2F(Form("fElectronYRes%d",i),"Distance (y) of EMCAL clusterin ID electrons;#phi;#Delta y",100,0,2*TMath::Pi(),1000,-20,20);
        fOutput->Add(fElectronYRes[i]);
        
        fElectronPositronZRes[i] = new TH2F(Form("fElectronPositronZRes%d",i),"Distance (z) of EMCAL clusterin ID electrons and positrons;#phi;#Delta z",100,0,2*TMath::Pi(),1000,-20,20);
        fOutput->Add(fElectronPositronZRes[i]);
        
        //Positrons
        
        fPositronPhiRes[i] = new TH2F(Form("fPositronPhiRes%d",i),"Distance (phi) of EMCAL cluster in ID positrons;#Delta#phi;#Delta#phi",100,0,2*TMath::Pi(),1000,-0.3,0.3);
        fOutput->Add(fPositronPhiRes[i]);
        
        fPositronEtaRes[i] = new TH2F(Form("fPositronEtaRes%d",i),"Distance (eta) of EMCAL clusterin ID positrons;#phi;#Delta#eta",100,0,2*TMath::Pi(),1000,-0.3,0.3);
        fOutput->Add(fPositronEtaRes[i]);
        
        fPositronXRes[i] = new TH2F(Form("fPositronXRes%d",i),"Distance (x) of EMCAL cluster in ID positrons;#phi;#Delta x",100,0,2*TMath::Pi(),1000,-20,20);
        fOutput->Add(fPositronXRes[i]);
        
        fPositronYRes[i] = new TH2F(Form("fPositronYRes%d",i),"Distance (y) of EMCAL clusterin ID positrons;#phi;#Delta y",100,0,2*TMath::Pi(),1000,-20,20);
        fOutput->Add(fPositronYRes[i]);
        
        //All tracks
        
        fAllMatchedTracksPhiRes[i] = new TH2F(Form("fAllMatchedTracksPhiRes%d",i),"Distance (phi) of EMCAL cluster in all matched tracks ;#Delta#phi;#Delta#eta",100,0,2*TMath::Pi(),100,-0.3,0.3);
        fOutput->Add(fAllMatchedTracksPhiRes[i]);
        
        fAllMatchedTracksEtaRes[i] = new TH2F(Form("fAllMatchedTracksEtaRes%d",i),"Distance (phi) of EMCAL cluster in all matched tracks ;#Delta#phi;#Delta#eta",100,0,2*TMath::Pi(),100,-0.3,0.3);
        fOutput->Add(fAllMatchedTracksEtaRes[i]);
        
        fElectronPhiResRU[i] = new TH2F(Form("fElectronPhiResRU%d",i),"Distance (phi) of EMCAL cluster in ID electrons RU;#phi;#Delta#phi",100,0,2*TMath::Pi(),1000,-0.3,0.3);
        fOutput->Add(fElectronPhiResRU[i]);
        
        fElectronEtaResRU[i] = new TH2F(Form("fElectronEtaResRU%d",i),"Distance (eta) of EMCAL clusterin ID electrons RU;#phi;#Delta#eta",100,0,2*TMath::Pi(),1000,-0.3,0.3);
        fOutput->Add(fElectronEtaResRU[i]);
        
        fElectronXResRU[i] = new TH2F(Form("fElectronXResRU%d",i),"Distance (x) of EMCAL cluster in ID electrons RU;#phi;#Delta x",100,0,2*TMath::Pi(),1000,-20,20);
        fOutput->Add(fElectronXResRU[i]);
        
        fElectronYResRU[i] = new TH2F(Form("fElectronYResRU%d",i),"Distance (y) of EMCAL clusterin ID electrons RU;#phi;#Delta y",100,0,2*TMath::Pi(),1000,-20,20);
        fOutput->Add(fElectronYResRU[i]);
        
        fElectronPositronZResRU[i] = new TH2F(Form("fElectronPositronZRes%dRU",i),"Distance (z) of EMCAL clusterin ID electrons and positrons RU;#phi;#Delta z",100,0,2*TMath::Pi(),1000,-20,20);
        fOutput->Add(fElectronPositronZResRU[i]);
        
        //Positrons
        
        fPositronPhiResRU[i] = new TH2F(Form("fPositronPhiResRU%d",i),"Distance (phi) of EMCAL cluster in ID positrons RU;#Delta#phi;#Delta#eta",100,0,2*TMath::Pi(),1000,-0.3,0.3);
        fOutput->Add(fPositronPhiResRU[i]);
        
        fPositronEtaResRU[i] = new TH2F(Form("fPositronEtaResRU%d",i),"Distance (eta) of EMCAL clusterin ID positrons RU;#phi;#Delta#eta",100,0,2*TMath::Pi(),1000,-0.3,0.3);
        fOutput->Add(fPositronEtaResRU[i]);
        
        fPositronXResRU[i] = new TH2F(Form("fPositronXResRU%d",i),"Distance (x) of EMCAL cluster in ID positrons RU;#phi;#Delta x",100,0,2*TMath::Pi(),1000,-20,20);
        fOutput->Add(fPositronXResRU[i]);
        
        fPositronYResRU[i] = new TH2F(Form("fPositronYResRU%d",i),"Distance (y) of EMCAL clusterin ID positrons RU;#phi;#Delta y",100,0,2*TMath::Pi(),1000,-20,20);
        fOutput->Add(fPositronYResRU[i]);
        
        
    }
    

    
    PostData(1, fOutput);
}


Bool_t AliAnalysisTaskEMCALAlig::FillHistograms()
{
    DoTrackLoop();
    return kTRUE;
}

void AliAnalysisTaskEMCALAlig::DoTrackLoop()
{
    AliClusterContainer* clusCont = GetClusterContainer(0);
    AliVCaloCells *cells = InputEvent()->GetEMCALCells();
    
    if (!clusCont)
    {
        AliError("No Cluster Container Available\n");
        return;
    }
    
    AliParticleContainer* partCont = 0;
    
    TIter next(&fParticleCollArray);
    
    while ((partCont = static_cast<AliParticleContainer*>(next()))) {
        
        UInt_t count = 0;
        
        for(auto part : partCont->accepted())
        {
            if (!part)
                continue;
            count++;
            
            const AliVTrack* track = static_cast<const AliVTrack*>(part);
            if (!track)
                continue;
            //Electron PID cuts: -1 to 3 sigma TPC
            Double_t TPCNSgima = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
            fTPCnSgima->Fill(track->Pt(),TPCNSgima);
            
            if (TPCNSgima < 1.0 || TPCNSgima > 3.0)
                continue;
            
            if (clusCont)
            {
                Int_t iCluster = track->GetEMCALcluster();
                if (iCluster >= 0)
                {
                    AliVCluster* cluster = clusCont->GetAcceptCluster(iCluster);
                    
                    if (cluster)
                    {
                        Double_t EoverP = cluster->GetNonLinCorrEnergy()/track->P();
                        fEOverP->Fill(track->Pt(),EoverP);
                        
                        //E/p cut
                        if (EoverP>= 0.8 && EoverP<=1.2)
                        {
                            if (track->Pt() > 8)
                            {
                                //M20cut for high pT
                                Double_t M20 = cluster->GetM20();
                                if ( M20<0.01 && M20>0.4)
                                    continue;
                            }
                            
                            //Cluster properties
                            Float_t  emcx[3];
                            cluster->GetPosition(emcx);
                            TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
                            Double_t emcphi = clustpos.Phi();
                            Double_t emceta = clustpos.Eta();
                            
                            if(emcphi < 0) emcphi = emcphi+(2*TMath::Pi());
                            
                            Int_t nCells=cluster->GetNCells();
                            Int_t iSupMod = -9;
                            Int_t ieta = -9;
                            Int_t iphi = -9;
                            Int_t icell = -9;
                            Bool_t Isshared	= kFALSE;
                            //Get the SM number
                            fEMCALRecoUtils->GetMaxEnergyCell(fEMCALGeo,fCaloCells,cluster,icell,iSupMod,ieta,iphi,Isshared);
                            Int_t SMNumber = iSupMod;
                            
                            TVector3 trackposOnEMCAL;
                            trackposOnEMCAL.SetPtEtaPhi(440,track->GetTrackEtaOnEMCal(),track->GetTrackPhiOnEMCal());
                            
                            Float_t phidiff = TVector2::Phi_mpi_pi(trackposOnEMCAL.Phi()-emcphi);
                            Float_t etadiff = trackposOnEMCAL.Eta() - emceta;
                            
                            Float_t xdiff = trackposOnEMCAL.X() - clustpos.X();
                            Float_t ydiff = trackposOnEMCAL.Y() - clustpos.Y();
                            Float_t zdiff = trackposOnEMCAL.Z() - clustpos.Z();
                            
                            
                            //Save Default residuals (from propagation to the EMCal Surface)
                            if (track->Charge()> 0)
                            {
                                //Track properties
                                fPPartPt->Fill(track->Pt());
                                fPPartPhi->Fill(track->Phi());
                                fPPartEta->Fill(track->Eta());
                                //Residual from matching
                                fPositronEtaRes[SMNumber]->Fill(emcphi,etadiff);
                                fPositronPhiRes[SMNumber]->Fill(emcphi,phidiff);
                                fPositronXRes[SMNumber]->Fill(emcphi,xdiff);
                                fPositronYRes[SMNumber]->Fill(emcphi,ydiff);
                                fElectronPositronZRes[SMNumber]->Fill(emcphi,zdiff);
                            }
                            else
                            {
                                fNPartPt->Fill(track->Pt());
                                fNPartPhi->Fill(track->Phi());
                                fNPartEta->Fill(track->Eta());
                                fElectronEtaRes[SMNumber]->Fill(emcphi,etadiff);
                                fElectronPhiRes[SMNumber]->Fill(emcphi,phidiff);
                                fElectronXRes[SMNumber]->Fill(emcphi,xdiff);
                                fElectronYRes[SMNumber]->Fill(emcphi,ydiff);
                                fElectronPositronZRes[SMNumber]->Fill(emcphi,zdiff);
                                
                            }
                            
                            //propagation using electron mass
                            
                            Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
                            track->PxPyPz(pxpypz);
                            track->XvYvZv(xyz);
                            track->GetCovarianceXYZPxPyPz(cv);
                            AliExternalTrackParam trackParam = AliExternalTrackParam(xyz,pxpypz,cv,track->Charge());
                            
                            Double_t trackPosExt[3] = {0.,0.,0.};
                            Double_t ElectronMass = 0.000510998910; //Electron mass in GeV
                            Float_t EtaResidualsForCrossCheck, PhiResidualsForCrossCheck;
                            
                            fEMCALRecoUtils->ExtrapolateTrackToPosition(&trackParam, emcx, ElectronMass, fEMCALRecoUtils->GetStep(), EtaResidualsForCrossCheck,PhiResidualsForCrossCheck);
                            trackParam.GetXYZ(trackPosExt);
                            
                            TVector3 trackposOnEMCALRU;
                            
                            trackposOnEMCALRU.SetXYZ(trackPosExt[0],trackPosExt[1],trackPosExt[2]);
                            
                            Double_t phidiffRU = TVector2::Phi_mpi_pi(trackposOnEMCALRU.Phi()-emcphi);
                            Double_t etadiffRU = trackposOnEMCALRU.Eta() - emceta;
                            
                            Double_t xdiffRU = trackposOnEMCALRU.X() - clustpos.X();
                            Double_t ydiffRU = trackposOnEMCALRU.Y() - clustpos.Y();
                            Double_t zdiffRU = trackposOnEMCALRU.Z() - clustpos.Z();
                            
                            if (track->Charge()>0)
                            {
                                fPositronEtaResRU[SMNumber]->Fill(emcphi,etadiffRU);
                                fPositronPhiResRU[SMNumber]->Fill(emcphi,phidiffRU);
                                fPositronXResRU[SMNumber]->Fill(emcphi,xdiffRU);
                                fPositronYResRU[SMNumber]->Fill(emcphi,ydiffRU);
                                fElectronPositronZResRU[SMNumber]->Fill(emcphi,zdiffRU);
                            }
                            else
                            {
                                fElectronEtaResRU[SMNumber]->Fill(emcphi,etadiffRU);
                                fElectronPhiResRU[SMNumber]->Fill(emcphi,phidiffRU);
                                fElectronXResRU[SMNumber]->Fill(emcphi,xdiffRU);
                                fElectronYResRU[SMNumber]->Fill(emcphi,ydiffRU);
                                fElectronPositronZResRU[SMNumber]->Fill(emcphi,zdiffRU);
                            }
                            
                        }
                    }
                }
                
            }
        }
        
    }
}

void AliAnalysisTaskEMCALAlig::ExecOnce()
{
    AliAnalysisTaskEmcal::ExecOnce();
    
    //EMCal utilits
    fEMCALGeo = AliEMCALGeometry::GetInstance();
    fEMCALRecoUtils  = new AliEMCALRecoUtils();
    fEMCALRecoUtils->InitParameters();
}

Bool_t AliAnalysisTaskEMCALAlig::Run()
{
    return kTRUE;
}

void AliAnalysisTaskEMCALAlig::Terminate(Option_t *)
{
}
