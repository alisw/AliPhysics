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

#include <AliVCluster.h>
#include <AliVParticle.h>
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPIDResponse.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisTaskEMCALAlig.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliExternalTrackParam.h"
#include "TVector3.h"
#include "AliTrackerBase.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEMCALAlig)
/// \endcond

AliAnalysisTaskEMCALAlig::AliAnalysisTaskEMCALAlig() :
AliAnalysisTaskEmcal(),
fEMCALRecoUtils(NULL),
fEMCALGeo(NULL),
fPIDResponse(NULL),
fElectronInformation(ElectronForAlignment()),
fElectronTree(NULL),
fTreeSuffix("")
{
    
}

AliAnalysisTaskEMCALAlig::AliAnalysisTaskEMCALAlig(const char *name) :
AliAnalysisTaskEmcal(name, kTRUE),
fEMCALRecoUtils(NULL),
fEMCALGeo(NULL),
fPIDResponse(NULL),
fElectronInformation(ElectronForAlignment()),
fElectronTree(NULL),
fTreeSuffix("")
{
    SetMakeGeneralHistograms(kTRUE);
    DefineOutput(2, TTree::Class());
}


AliAnalysisTaskEMCALAlig::~AliAnalysisTaskEMCALAlig()
{
    if (fEMCALRecoUtils)
        delete fEMCALRecoUtils;
    if (fElectronTree)
        delete fElectronTree;
}


void AliAnalysisTaskEMCALAlig::UserCreateOutputObjects()
{
    AliAnalysisTaskEmcal::UserCreateOutputObjects();
    
    fPIDResponse = fInputHandler->GetPIDResponse();
    
    TString name_tree("electron_information");
    name_tree += fTreeSuffix;
    
    fElectronTree = new TTree(name_tree,name_tree);
    fElectronTree->Branch("electrons", &fElectronInformation);
    PostData(2, fElectronTree);
}


Bool_t AliAnalysisTaskEMCALAlig::FillHistograms()
{
    DoTrackLoop();
    return kTRUE;
}

void AliAnalysisTaskEMCALAlig::DoTrackLoop()
{
    AliClusterContainer* clusCont = GetClusterContainer(0);
    
    if (!clusCont)
    {
        AliError("No Cluster Container Available\n");
        return;
    }
    
    AliParticleContainer* partCont = 0;
    
    TIter next(&fParticleCollArray);
    
    while ((partCont = static_cast<AliParticleContainer*>(next()))) {        
        for(auto part : partCont->accepted())
        {
            if (!part)
                continue;
            
            const AliVTrack* track = static_cast<const AliVTrack*>(part);
            if (!track)
                continue;
            
            //-1.5 to 3.5 sigma TPC to reduce the size of the trees
            Double_t n_sigma_electron_TPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
            
            if (n_sigma_electron_TPC < -1.5 || n_sigma_electron_TPC > 3.5)
                continue;
            
            Int_t iCluster = track->GetEMCALcluster();
            
            if (iCluster < 0)
                continue;
            
            AliVCluster* cluster = clusCont->GetAcceptCluster(iCluster);
            
            if (!cluster)
                continue;
            
            Double_t EoverP = cluster->E()/track->P();
            
            //Loose E/p cut to reduce tree
            if (EoverP<0.7 || EoverP>1.3)
                continue;
            
            //Cluster properties
            Float_t  emcx[3];
            cluster->GetPosition(emcx);
            TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
            
            Int_t iSupMod = -9;
            Int_t ieta = -9;
            Int_t iphi = -9;
            Int_t icell = -9;
            Bool_t Isshared	= kFALSE;
            
            //Get the SM number
            fEMCALRecoUtils->GetMaxEnergyCell(fEMCALGeo,fCaloCells,cluster,icell,iSupMod,ieta,iphi,Isshared);
                        
            //Propagation using electron mass
            
            Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
            track->PxPyPz(pxpypz);
            track->XvYvZv(xyz);
            track->GetCovarianceXYZPxPyPz(cv);
            AliExternalTrackParam trackParam = AliExternalTrackParam(xyz,pxpypz,cv,track->Charge());
            
            Double_t trackPosExt[3] = {0.,0.,0.};
            Double_t ElectronMass = 0.000510998910; //Electron mass in GeV

            if (!AliTrackerBase::PropagateTrackToBxByBz(&trackParam, clustpos.Perp(), ElectronMass, fEMCALRecoUtils->GetStep(),kTRUE, 0.8, -1)) 
    			continue;
            
            trackParam.GetXYZ(trackPosExt);
            
            TVector3 trackposOnEMCAL;
            trackposOnEMCAL.SetXYZ(trackPosExt[0],trackPosExt[1],trackPosExt[2]);
            
            //Save to Tree
            //Basic track properties
            fElectronInformation.charge = track->Charge();
            fElectronInformation.pt = track->Pt();
            fElectronInformation.pz = track->Pz();
            fElectronInformation.eta_track = track->Eta();
            fElectronInformation.phi_track = track->Phi();
            fElectronInformation.zvtx = fVertex[2];
            
            //cluster properties
            fElectronInformation.energy = cluster->E();
            fElectronInformation.M20 = cluster->GetM20();
            fElectronInformation.M02 = cluster->GetM02();
            fElectronInformation.eta_cluster = clustpos.Eta();
			fElectronInformation.phi_cluster = clustpos.Phi();
            fElectronInformation.x_cluster = clustpos.X();
            fElectronInformation.y_cluster = clustpos.Y();
            fElectronInformation.z_cluster = clustpos.Z();
            
                   
            //mathing properties using electron mass
            fElectronInformation.x_track = trackposOnEMCAL.X();
            fElectronInformation.y_track = trackposOnEMCAL.Y();
            fElectronInformation.z_track = trackposOnEMCAL.Z();
            
            //aditional information
            fElectronInformation.super_module_number = iSupMod;
            fElectronInformation.distance_bad_channel = cluster->GetDistanceToBadChannel();
            
            AliVCaloCells *cells = InputEvent()->GetEMCALCells();
            fElectronInformation.is_in_fid_region = fEMCALRecoUtils->CheckCellFiducialRegion(fEMCALGeo,cluster,cells);
            
            //PID properties
            fElectronInformation.n_sigma_electron_TPC = n_sigma_electron_TPC;
            
            fElectronTree->Fill();
            
        }
        
        
    }
    
    PostData(2, fElectronTree);
}

void AliAnalysisTaskEMCALAlig::ExecOnce()
{
    AliAnalysisTaskEmcal::ExecOnce();
    
    //EMCal utilits
    fEMCALGeo = AliEMCALGeometry::GetInstance();
    fEMCALRecoUtils  = new AliEMCALRecoUtils();
    fEMCALRecoUtils->InitParameters();
    fEMCALRecoUtils->SetNumberOfCellsFromEMCALBorder(1);
}

Bool_t AliAnalysisTaskEMCALAlig::Run()
{
    return kTRUE;
}

void AliAnalysisTaskEMCALAlig::Terminate(Option_t *)
{
}
