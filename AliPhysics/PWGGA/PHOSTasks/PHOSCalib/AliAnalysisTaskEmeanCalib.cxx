#include <stdio.h>

// ROOT includes
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TList.h"
#include "TStyle.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "AliInputEventHandler.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskEmeanCalib.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSCalibData.h"
#include "AliPHOSEsdCluster.h"
#include "AliESDCaloCells.h"
#include "AliOADBContainer.h"

ClassImp(AliAnalysisTaskEmeanCalib);

//__________________________________________________________________________
AliAnalysisTaskEmeanCalib::AliAnalysisTaskEmeanCalib() :
AliAnalysisTaskSE(),fESDEvent(0x0),fOutput(0x0)
,fPHOSGeo(0),fCalibData(),fRunNumber(-999),hmpt_diff(0x0)
{
    // Default ctor.
    
    for (Int_t imod=0; imod<5; imod++) {
        for (Int_t ix=0; ix<64; ix++) {
            for (Int_t iz=0; iz<56; iz++) {
                hCellAmplitudes[imod][ix][iz]=0;
            }
        }
    }
    
    for (Int_t imod=0; imod<5; imod++) {
        hCellMeanAmps[imod]=0;
    }
    
    for (Int_t imod=0; imod<5; imod++) {
        for (Int_t ix=0; ix<64; ix++) {
            for (Int_t iz=0; iz<56; iz++) {
                fCC[imod][ix][iz] = 1.;
            }
        }
    }
    
}

//__________________________________________________________________________
AliAnalysisTaskEmeanCalib::AliAnalysisTaskEmeanCalib(const char *name) :
AliAnalysisTaskSE(name),fESDEvent(0x0),fOutput(0x0)
,fPHOSGeo(0),fCalibData(),fRunNumber(-999),hmpt_diff(0x0)
{
    // Constructor. Initialization of Inputs and Outputs
    printf("#### TASK CREATED!! ########\n");
    DefineOutput(1, TList::Class());
    
    for (Int_t imod=0; imod<5; imod++) {
        for (Int_t ix=0; ix<64; ix++) {
            for (Int_t iz=0; iz<56; iz++) {
                hCellAmplitudes[imod][ix][iz]=0;
            }
        }
    }
    
    for (Int_t imod=0; imod<5; imod++) {
        hCellMeanAmps[imod]=0;
    }
    
    for (Int_t imod=0; imod<5; imod++) {
        for (Int_t ix=0; ix<64; ix++) {
            for (Int_t iz=0; iz<56; iz++) {
                fCC[imod][ix][iz] = 1.;
            }
        }
    }

}

//___________________________________________________________________________
void AliAnalysisTaskEmeanCalib::UserCreateOutputObjects(){
    printf("Starting UserCreateOutputObjects!\n");
    
    // Output objects creation
    fOutput = new TList();
    fOutput->SetOwner();

    char hnam[80];
    char htit[80];
    
    // X-pT,Y-Minv. in modules
    for (Int_t i=0; i<4; i++) {
        
        sprintf(hnam,"hmpt_SM%d",i+1);
        sprintf(htit,"2-cluster inv. mass vs pT in SM%d",i+1);
        hmpt[i] = new TH2F(hnam,htit,100,0.,10.,200,0.,1.);
        
        sprintf(hnam,"hmpt_orig_SM%d",i+1);
        sprintf(htit,"2-cluster inv. mass vs pT before recalibration in SM%d",i+1);
        hmpt_orig[i] = new TH2F(hnam,htit,100,0.,10.,200,0.,1.);
        
        fOutput->Add(hmpt[i]);
        fOutput->Add(hmpt_orig[i]);
        
    }

    hmpt[4] = new TH2F("hmpt","2-cluster inv. mass vs pT",100,0.,10.,200,0.,1.);
    fOutput->Add(hmpt[4]);
    
    hmpt_orig[4] = new TH2F("hmpt_orig","2-cluster inv. mass vs pT before recalibration",100,0.,10.,200,0.,1.);
    fOutput->Add(hmpt_orig[4]);

    hmpt_diff = new TH2F("hmpt_diff","Mgg vs pT for clusters in different PHOS modules",100,0.,10.,200,0.,1.);
    fOutput->Add(hmpt_diff);

    // Required both here and in UserExec()
    PostData(1, fOutput);
}

//___________________________________________________________________________
void AliAnalysisTaskEmeanCalib::UserExec(Option_t *)
{
    // Execute analysis for current InputEvent
    printf("RUNNING UserExec!!\n");
    
    fESDEvent = dynamic_cast<AliESDEvent*> (InputEvent());
    
    if ( ! fESDEvent ) {
        AliError ("ESD event not found. Nothing done!");
        return;
    }
    
    if(fRunNumber != fESDEvent->GetRunNumber()) { //this should run only at first call of UserExec()
        fRunNumber = fESDEvent->GetRunNumber();
        SetGeometry();
        SetMisalignment();
    }
    
    Int_t module = -1; // If <0 - all modules
    
    Int_t minCells=3; //Ncells>2
    Double_t minEclu=0.3; // min. cluster energy = 0.3 GeV
    
    Int_t nCols = 56;
    Int_t nRows = 64;
    
    Int_t absId=0 ;
    TLorentzVector p1;
    TLorentzVector p2;
    TLorentzVector p12;
    TLorentzVector p1o;
    TLorentzVector p2o;
    TLorentzVector p12o;
    Int_t relid[4] ;
    Int_t phosClu=0;
    Int_t phosMod1 = -1;
    Int_t phosMod2 = -1;
    UShort_t *absIds=0;
    Int_t nCells=0;
    TVector3 vtx(0,0,0); // vertex!
    Float_t logWeight = 4.5;
    
    AliESDCaloCells *phsCells = fESDEvent->GetPHOSCells();
    Int_t multCells = phsCells->GetNumberOfCells();
    
    Int_t multClu = fESDEvent->GetNumberOfCaloClusters();
    phosClu = 0;
    
    for (Int_t iCell=0; iCell<multCells; iCell++) {
        
        Int_t cellAbsId = phsCells->GetCellNumber(iCell);
        fPHOSGeo->AbsToRelNumbering(cellAbsId, relid);
        
        Int_t cellMod = relid[0];
        Int_t cellX   = relid[2];
        Int_t cellZ   = relid[3];
        
        if(!hCellMeanAmps[cellMod-1]) {
            char mname[128], mtitle[128];
            sprintf(mname,"hCellMeanAmp%d",cellMod);
            sprintf(mtitle,"Cell Occupancies in Module %d",cellMod);
            hCellMeanAmps[cellMod-1] = new TH2F(mname,mtitle,64,0.,64,56,0.,56);
            fOutput->Add(hCellMeanAmps[cellMod-1]);
        }
        
        if(phsCells->GetAmplitude(iCell)) {
            
            if(!hCellAmplitudes[cellMod-1][cellX-1][cellZ-1]) {
                char chname[128],chtitle[128];
                sprintf(chname,"cell_m%d_x%d_z%d",cellMod,cellX,cellZ);
                sprintf(chtitle,"Amplitude in cell m=%d, x=%d, z=%d",cellMod,cellX,cellZ);
                hCellAmplitudes[cellMod-1][cellX-1][cellZ-1] = new TH1F(chname,chtitle,2000,0.,2.);
                fOutput->Add(hCellAmplitudes[cellMod-1][cellX-1][cellZ-1]);
            }
            
//            Float_t cc_i = fCalibData->GetADCchannelEmc(cellMod,cellZ,cellX);
//            Float_t cc_i = 1.;
            hCellAmplitudes[cellMod-1][cellX-1][cellZ-1]->Fill(fCC[cellMod-1][cellX-1][cellZ-1]*phsCells->GetAmplitude(iCell));
            hCellMeanAmps[cellMod-1]->Fill(cellX,cellZ);
        }
    }
    
    for (Int_t i=0; i<multClu; i++) {
        
        AliESDCaloCluster *c1 = fESDEvent->GetCaloCluster(i);
        if(!c1->IsPHOS()) continue;
        phosClu++;
        
        if(c1->E()<minEclu) continue;
        
        absIds = c1->GetCellsAbsId();
        fPHOSGeo->AbsToRelNumbering(absIds[0], relid) ;
        phosMod1 = relid[0];
        if(module>0 && phosMod1 != module) continue;
        
        Int_t nCells      = c1->GetNCells();
        if(nCells<minCells) continue;
        
        AliPHOSEsdCluster clu1(*c1);
        clu1.Recalibrate(fCalibData, phsCells);
        clu1.EvalAll(logWeight,vtx);
        clu1.EnergyCorrection() ;
        
        c1->GetMomentum(p1o,0);
        clu1.GetMomentum(p1,0);
        
        for (Int_t j=i; j<multClu; j++) {
            
            AliESDCaloCluster *c2 = fESDEvent->GetCaloCluster(j);
            if(!c2->IsPHOS()) continue; if(c2->IsEqual(c1)) continue;
            
            if(c2->E()<minEclu) continue;
            
            absIds = c2->GetCellsAbsId();
            fPHOSGeo->AbsToRelNumbering(absIds[0], relid) ;
            phosMod2 = relid[0];
            if(module>0 && phosMod2 != module) continue;
            
            nCells      = c2->GetNCells();
            if(nCells<minCells) continue;
            
            
            AliPHOSEsdCluster clu2(*c2);
            clu2.Recalibrate(fCalibData, phsCells);
            clu2.EvalAll(logWeight,vtx);
            clu2.EnergyCorrection() ;
            
            c2->GetMomentum(p2o,0);
            clu2.GetMomentum(p2,0);
            
            p12 = p1+p2;
            
            if(phosMod1==phosMod2)
                hmpt[phosMod1-1]->Fill(p12.Pt(),p12.M());
            
            hmpt[4]->Fill(p12.Pt(),p12.M());
            
            if(phosMod1 != phosMod2 && p1.E()>0.5 && p2.E()>0.5)
                hmpt_diff->Fill(p12.Pt(),p12.M());
            
            p12o = p1o+p2o;
            
            if(phosMod1==phosMod2)
                hmpt_orig[phosMod1-1]->Fill(p12o.Pt(),p12o.M());
            
            hmpt_orig[4]->Fill(p12o.Pt(),p12o.M());
            
        }
        
    }
    
    //    printf("\t--- Run %d, Event %d, %d cluster(s), %d in PHOS. ---\n",
    //           runNum,iEvent,multClu,phosClu);
    
    // Required both here and in UserCreateOutputObjects()
    PostData(1, fOutput);
}

//___________________________________________________________________________
void AliAnalysisTaskEmeanCalib::SetCalibrations(AliPHOSCalibData* c)
{
    if(fCalibData) delete fCalibData;
    fCalibData = c;
    
    for (Int_t imod=0; imod<5; imod++) {
        for (Int_t ix=0; ix<64; ix++) {
            for (Int_t iz=0; iz<56; iz++) {
                fCC[imod][ix][iz] = fCalibData->GetADCchannelEmc(imod+1,iz+1,ix+1);
            }
        }
    }

//    Float_t cc_i = fCalibData->GetADCchannelEmc(2,11,12);
//    printf("CC_i=%.3f\n",cc_i);
}


//___________________________________________________________________________
void AliAnalysisTaskEmeanCalib::Terminate(Option_t *)
{
    // Display the Pt histogram
    
//    fOutput = dynamic_cast<TList*> (GetOutputData(1));
//    if (fOutput) {
//        fHistThePt = dynamic_cast<TH1D *>( fOutput->FindObject("hThePt") );
//        if (fHistThePt) {
//            fHistThePt->Draw();
//        }
//    }
}

//___________________________________________________________________________
void AliAnalysisTaskEmeanCalib::SetGeometry()
{
    fPHOSGeo =  AliPHOSGeometry::GetInstance("IHEP") ;
    
//    if(!fPHOSGeo){
//        
//        AliOADBContainer geomContainer("phosGeo");
//        geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
//        TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(fRunNumber,"PHOSRotationMatrixes");
//        fPHOSGeo =  AliPHOSGeometry::GetInstance("IHEP") ;
//        
//        for(Int_t mod=0; mod<5; mod++) {
//            if(!matrixes->At(mod)) {
//                if( fDebug )
//                    AliInfo(Form("No PHOS Matrix for mod:%d, geo=%p\n", mod, fPHOSGeo));
//                continue;
//            }
//            else {
//                fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;
//                if( fDebug >1 )
//                    AliInfo(Form("Adding PHOS Matrix for mod:%d, geo=%p\n", mod, fPHOSGeo));
//            }
//        }
//    }
    
}

//___________________________________________________________________________
void AliAnalysisTaskEmeanCalib::SetMisalignment()
{
    // sets the misalignment vertex for ESD events.
    
    for(Int_t mod=0; mod<5; mod++) {
        const TGeoHMatrix* modMatrix = fESDEvent->GetPHOSMatrix(mod);
        if( ! modMatrix) {
            if( fDebug )
                AliInfo(Form("no PHOS Geometric Misalignment Matrix for module %d", mod));
            continue    ;
        }
        else    {
            fPHOSGeo->SetMisalMatrix(modMatrix, mod);
            if( fDebug )
                AliInfo(Form("PHOS Geometric Misalignment Matrix set for module %d", mod));
        }
    }
}
