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
ElectronInformation(ElectronForAlignment()),
ElectronTree(NULL)
{
    
}

AliAnalysisTaskEMCALAlig::AliAnalysisTaskEMCALAlig(const char *name) :
AliAnalysisTaskEmcal(name, kTRUE),
fEMCALRecoUtils(NULL),
fEMCALGeo(NULL),
fPIDResponse(NULL),
ElectronInformation(ElectronForAlignment()),
ElectronTree(NULL)
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
    
    ElectronTree = new TTree("electron_information","electron_information");
    ElectronTree->Branch("electrons", &ElectronInformation);
    fOutput->Add(ElectronTree);
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
    //AliVCaloCells *cells = InputEvent()->GetEMCALCells();
    
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
            
            if (TPCNSgima < -1.5 || TPCNSgima > 3.5)
                continue;
            
            Int_t iCluster = track->GetEMCALcluster();
            
            if (iCluster < 0)
                continue;
            
            AliVCluster* cluster = clusCont->GetAcceptCluster(iCluster);
            
            if (!cluster)
                continue;
            
            
            Double_t EoverP = cluster->GetNonLinCorrEnergy()/track->P();
                //E/p cut
            if (EoverP<0.7 || EoverP>1.3)
                continue;
            
            //Cluster properties
            Float_t  emcx[3];
            cluster->GetPosition(emcx);
            TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
            Double_t emcphi = clustpos.Phi();
            Double_t emceta = clustpos.Eta();
            
            if(emcphi < 0) emcphi = emcphi+(2*TMath::Pi());
            
            //Int_t nCells=cluster->GetNCells();
            Int_t iSupMod = -9;
            Int_t ieta = -9;
            Int_t iphi = -9;
            Int_t icell = -9;
            Bool_t Isshared	= kFALSE;
            //Get the SM number
            fEMCALRecoUtils->GetMaxEnergyCell(fEMCALGeo,fCaloCells,cluster,icell,iSupMod,ieta,iphi,Isshared);
            
            TVector3 trackposOnEMCAL;
            trackposOnEMCAL.SetPtEtaPhi(440,track->GetTrackEtaOnEMCal(),track->GetTrackPhiOnEMCal());
            
            Float_t phidiff = TVector2::Phi_mpi_pi(trackposOnEMCAL.Phi()-emcphi);
            Float_t etadiff = trackposOnEMCAL.Eta() - emceta;
            
            Float_t xdiff = trackposOnEMCAL.X() - clustpos.X();
            Float_t ydiff = trackposOnEMCAL.Y() - clustpos.Y();
            Float_t zdiff = trackposOnEMCAL.Z() - clustpos.Z();
            
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
            
            
            ElectronInformation.charge = track->Charge();
            ElectronInformation.pt = track->Pt();
            ElectronInformation.pz = track->Pz();
            ElectronInformation.eta_track = track->Eta();
            ElectronInformation.phi_track = track->Phi();
            
            //cluster properties
            ElectronInformation.energy = cluster->GetNonLinCorrEnergy();
            ElectronInformation.M20 = cluster->GetM20();
            ElectronInformation.M02 = cluster->GetM02();
            ElectronInformation.eta_cluster = clustpos.Eta();
            ElectronInformation.phi_cluster = clustpos.Phi();
            
            //mathing properties using default matcher
            ElectronInformation.x_resitual_def = xdiff;
            ElectronInformation.y_resitual_def = ydiff;
            ElectronInformation.z_resitual_def = zdiff;
            ElectronInformation.phi_resitual_def = cluster->GetTrackDx();
            ElectronInformation.eta_resitual_def = cluster->GetTrackDz();
                   
            //mathing properties using electron mass
            ElectronInformation.x_resitual_e = xdiffRU;
            ElectronInformation.y_resitual_e = ydiffRU;
            ElectronInformation.z_resitual_e = zdiffRU;
            ElectronInformation.phi_resitual_e = phidiffRU;
            ElectronInformation.eta_resitual_e = etadiffRU;
            
            ElectronInformation.super_module_number = iSupMod;
            //PID properties
            ElectronInformation.TPCNSigmaElectron = TPCNSgima;
            
            ElectronTree->Fill();
            
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
