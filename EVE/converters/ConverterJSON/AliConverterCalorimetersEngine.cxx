#include "AliConverterCalorimetersEngine.h"

#include "AliGeomManager.h"
#include "AliCDBManager.h"

#include <TGeoManager.h>

#include <iostream>

using namespace std;


/// In the ctor of the CalorimetersEngine calorimeters geometry is set up according to the event
AliConverterCalorimetersEngine::AliConverterCalorimetersEngine(AliESDEvent *event) :
    fESDEvent(event), fCaloCluster(0), fGeomEM(0), fGeomPH(0), fCellsEM(0), fCellsPH(0)
{
    AssertMagField();
    AssertGeometry();
    
    cout<<"\nGeometry set\n"<<endl;
    
    SetUpEMCALGeometry();
    SetUpPHOSGeometry();
    
    cout<<"\nEMCal and PHOS geometry set\n"<<endl;
    
    fCellsEM = fESDEvent->GetEMCALCells();
    fCellsPH = fESDEvent->GetPHOSCells();
}

AliConverterCalorimetersEngine::~AliConverterCalorimetersEngine()
{
    if(fCaloCluster)
        delete fCaloCluster;
    if(fGeomEM)
        delete fGeomEM;
    if(fGeomPH)
        delete fGeomPH;
    if(fCellsEM)
        delete fCellsEM;
    if(fCellsPH)
        delete fCellsPH;
}


/// This method populates AliMinimalisticEvent with information regarding readouts from calorimeters associated with
/// tracks.
void AliConverterCalorimetersEngine::PopulateEventWithCaloClusters(AliMinimalisticEvent &event)
{
    Int_t absIdEMaxCell = -1;
    Int_t id = -1;
    Float_t eMaxCell = 0.;
    Float_t amp = -1;
    Float_t phi = .0;
    Float_t eta = .0;
    Double_t vertex[] = {0., 0., 0.}; /// Vertex container
    
    for (Int_t cluster_id = 0; cluster_id < fESDEvent->GetNumberOfCaloClusters(); cluster_id++)
    {
        fCaloCluster = fESDEvent->GetCaloCluster(cluster_id);
        fCaloCluster->GetMomentum(fClusterMomentum, vertex);
        
        absIdEMaxCell = -1;
        eMaxCell      = 0.;
        GetMaxEnergyCellAbsId(absIdEMaxCell, eMaxCell) ;
        
        if(IsBadCluster(absIdEMaxCell, eMaxCell))
            continue;
        
        if(fCaloCluster->IsEMCAL()){
            AddEMCALClustersToEvent(event, id, amp, phi, eta);
        } else {
            AddPHOSCalClusterToEvent(event, id, amp);
        }
    }
}

void AliConverterCalorimetersEngine::AddPHOSCalClusterToEvent(AliMinimalisticEvent &event, Int_t &id, Float_t &amp) {
    for(Int_t icell = 0; icell < fCaloCluster->GetNCells(); icell++){
        id  = fCaloCluster->GetCellAbsId(icell);
        amp = fCellsPH->GetCellAmplitude(id); // GeV

        Int_t relId[4];
        Float_t xCell, zCell;

        fGeomPH->AbsToRelNumbering(id, relId);
        fGeomPH->RelPosInModule(relId, xCell, zCell);

        cout << "PHOS\tEta: " << fClusterMomentum.Eta() << "\tPhi: " << GetPhi(fClusterMomentum.Phi()) << "\tE: " << fClusterMomentum.Energy() << endl;
        AliMinimalisticCaloCluster cluster(
                fGeomEM->GetIPDistance(),
                GetPhi(fClusterMomentum.Phi()),
                fClusterMomentum.Eta(),
                TMath::Pi() / 80.,  // Uwazam, ze nalezy zapisac te wartosc do zmiennej, nie bardzo wiadomo skad sie bierze taka wartosc
                3./200.,          //                    -"-
                fClusterMomentum.Energy()
        );
        event.AddCaloCluster(cluster);
    }
}

void AliConverterCalorimetersEngine::AddEMCALClustersToEvent(
        AliMinimalisticEvent &event, Int_t id, Float_t amp, Float_t phi, Float_t eta
)
{
    for (Int_t icell = 0; icell < fCaloCluster->GetNCells(); icell++){
        id  = fCaloCluster->GetCellAbsId(icell);
        amp = fCellsEM->GetCellAmplitude(id); // GeV

        fGeomEM->EtaPhiFromIndex(id, eta, phi);
        cout << "Eta: " << eta << "\tPhi: " << GetPhi(phi) << "\tE: " << amp << endl;
        AliMinimalisticCaloCluster cluster(
                fGeomEM->GetIPDistance(),
                GetPhi(phi),
                eta,
                TMath::Pi() / 80.,
                3./200.,
                amp
        );
        event.AddCaloCluster(cluster);
    }
}

void AliConverterCalorimetersEngine::AssertGeometry()
{
    if (AliGeomManager::GetGeometry() == 0)
    {
        gGeoManager = 0;
        
        AliGeomManager::LoadGeometry();
        if ( ! AliGeomManager::GetGeometry())
        {
            std::cerr<<"\n\nCan not load geometry.\n\n"<<std::endl;
        }
        if ( ! AliGeomManager::ApplyAlignObjsFromCDB("ITS TPC TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO ACORDE"))
        {
            std::cerr<<"\n\nMismatch of alignable volumes. Proceeding.\n\n"<<std::endl;
        }
        AliGeomManager::GetGeometry()->DefaultColors();
    }
    gGeoManager = AliGeomManager::GetGeometry();
}

void AliConverterCalorimetersEngine::AssertMagField()
{
    AliCDBManager* cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage(Form("local://%s/../src/OCDB",gSystem->Getenv("ALICE_ROOT")));
    if (!cdb->IsDefaultStorageSet())
    {
        std::cerr<<"\n\nCould not set CDB.\n\n"<<std::endl;
        return;
    }
    cdb->SetRun(fESDEvent->GetRunNumber());
}


float AliConverterCalorimetersEngine::GetECross(bool isEMCAL, int imod, int icol, int irow)
{

    Float_t  ecell1 =  0, ecell2  = 0, ecell3  = 0, ecell4  = 0;
    Int_t    absId1 = -1, absId2 = -1, absId3 = -1, absId4 = -1;
    
    AliESDCaloCells * cells = 0;
    
    if ( isEMCAL )
    {
        cells = fCellsEM;
        
        Int_t rowMax = AliEMCALGeoParams::fgkEMCALRows;
        Int_t colMax = AliEMCALGeoParams::fgkEMCALCols;
        
        if(imod == 11 || imod == 10 || imod == 18 || imod == 19) rowMax = AliEMCALGeoParams::fgkEMCALRows/3;
        if(imod > 12 && imod < 18) colMax = AliEMCALGeoParams::fgkEMCALCols*2/3;
        if( irow < rowMax) absId1 = fGeomEM->GetAbsCellIdFromCellIndexes(imod, irow + 1, icol);
        if( irow > 0     ) absId2 = fGeomEM->GetAbsCellIdFromCellIndexes(imod, irow - 1, icol);
        if( icol < colMax-1 )   absId3 = fGeomEM->GetAbsCellIdFromCellIndexes(imod, irow, icol + 1);
        if( icol > 0 )          absId4 = fGeomEM->GetAbsCellIdFromCellIndexes(imod, irow, icol - 1);
        
        if(imod > 11 && imod < 18)
        {
            if( icol == colMax - 1 && !(imod % 2) )
            {
                absId3 = fGeomEM->GetAbsCellIdFromCellIndexes(imod+1, irow, 0);
                absId4 = fGeomEM->GetAbsCellIdFromCellIndexes(imod  , irow, icol-1);
            }
            else if( icol == 0 && imod % 2 )
            {
                absId3 = fGeomEM->GetAbsCellIdFromCellIndexes(imod  , irow, icol+1);
                absId4 = fGeomEM->GetAbsCellIdFromCellIndexes(imod-1, irow, colMax-1);
            }
        }
    } else { // PHOS
        cells = fCellsPH;
        
        Int_t relId1[] = { imod+1, 0, irow+1, icol   };
        Int_t relId2[] = { imod+1, 0, irow-1, icol   };
        Int_t relId3[] = { imod+1, 0, irow  , icol+1 };
        Int_t relId4[] = { imod+1, 0, irow  , icol-1 };
        
        fGeomPH->RelToAbsNumbering(relId1, absId1);
        fGeomPH->RelToAbsNumbering(relId2, absId2);
        fGeomPH->RelToAbsNumbering(relId3, absId3);
        fGeomPH->RelToAbsNumbering(relId4, absId4);
    }
    
    if(absId1 > 0 )
        ecell1 = cells->GetCellAmplitude(absId1);
    if(absId2 > 0 )
        ecell2 = cells->GetCellAmplitude(absId2);
    if(absId3 > 0 )
        ecell3 = cells->GetCellAmplitude(absId3);
    if(absId4 > 0 )
        ecell4 = cells->GetCellAmplitude(absId4);
    
    return ecell1+ecell2+ecell3+ecell4;
}

void AliConverterCalorimetersEngine::GetModuleNumberColAndRow(int absId, bool isEMCAL,int &iSM, int &icol, int &irow)
{
    if (absId < 0)
        return;

    Int_t imod = -1;

    if (isEMCAL) {
        Int_t iTower = -1;
        Int_t iIphi = -1;
        Int_t iIeta = -1;
        fGeomEM->GetCellIndex(absId, iSM, iTower, iIphi, iIeta);
        fGeomEM->GetCellPhiEtaIndexInSModule(iSM, iTower, iIphi, iIeta, irow, icol);
    } else {
        Int_t relId[4];
        fGeomPH->AbsToRelNumbering(absId,relId);
        irow = relId[2];
        icol = relId[3];
        iSM  = relId[0] - 1;// Czemu -1?
    }
}

void AliConverterCalorimetersEngine::GetMaxEnergyCellAbsId(int &absId, float &eMax)
{
    Double_t eCell       = -1.;
    Int_t    cellAbsId   = -1 ;
    AliESDCaloCells*  cells = fCaloCluster->IsEMCAL() ? fCellsEM : fCellsPH;
    
    for (Int_t iDig=0; iDig < fCaloCluster->GetNCells(); iDig++)
    {
        cellAbsId = fCaloCluster->GetCellAbsId(iDig);
        eCell  = cells->GetCellAmplitude(cellAbsId);
        if(eCell > eMax)
        {
            eMax  = eCell;
            absId = cellAbsId;
        }
    }
}

bool AliConverterCalorimetersEngine::IsBadCluster(int absId, float eMax)
{
    Int_t   nMinCellsCut[] = { 3   , 2    };    /// Number of cells in cluster must be larger than this value.
    Int_t   nMaxCellsCut[] = { 60  , 30   };    /// Number of cells in cluster must be smaller than this value.
    Float_t    exoCut   = 0.95;                 /// Reject clusters with this exoticity value.
    Float_t    m20lowCut[] = { 0.50, -1.0 };    /// Cluster shower shape lower axis must be larger than this value.
    Float_t    m02higCut[] = { 7.00, 7.00 };    /// Cluster shower shape major axis must be larger than this value.
    Float_t    m02lowCut[] = { 0.70, 0.10 };    /// Cluster shower shape major axis must be larger than this value.
    Float_t    energyCut[] = { 0.30, 0.30 };    /// Cluster energy must be larger than this value.
    
    Bool_t isEMCAL = fCaloCluster->IsEMCAL();
    
    if(fCaloCluster->GetNCells() < nMinCellsCut[isEMCAL])
        return kTRUE ;
    if(fCaloCluster->GetNCells() > nMaxCellsCut[isEMCAL])
        return kTRUE ;
    if(fCaloCluster->E()         < energyCut[isEMCAL])
        return kTRUE ;
    if(fCaloCluster->GetM02()    < m02lowCut[isEMCAL])
        return kTRUE ;
    if(fCaloCluster->GetM02()    > m02higCut[isEMCAL])
        return kTRUE ;
    if(fCaloCluster->GetM20()    < m20lowCut[isEMCAL])
        return kTRUE ;
    if(eMax                      < (energyCut[isEMCAL] / 2.) || absId < 0 )
        return kTRUE;
    
    Int_t ism  =  -1;
    Int_t icol = -1;
    Int_t irow = -1;
    GetModuleNumberColAndRow(absId, isEMCAL, ism, icol, irow);
    Float_t eCross = GetECross(isEMCAL, ism, icol, irow);
    if((1 - (eCross / eMax)) > exoCut)
        return kTRUE;
    
    return kFALSE;
}

void AliConverterCalorimetersEngine::SetUpEMCALGeometry()
{
    if(!fGeomEM) fGeomEM  = AliEMCALGeometry::GetInstance();
    if(!fGeomEM)fGeomEM  = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM");
    
    for(Int_t mod = 0; mod < (fGeomEM->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
        if( fESDEvent->GetEMCALMatrix(mod) )
            fGeomEM->SetMisalMatrix(fESDEvent->GetEMCALMatrix(mod),mod) ;
        else
            fGeomEM->SetMisalMatrix((new TGeoHMatrix),mod) ;
    } // loop over super modules
}

void AliConverterCalorimetersEngine::SetUpPHOSGeometry()
{
    if(!fGeomPH) fGeomPH  = AliPHOSGeometry::GetInstance();
    if (!fGeomPH)fGeomPH  = AliPHOSGeometry::GetInstance("Run2");
    
    for(Int_t mod = 0; mod < 5; mod++)
    {
        if(fESDEvent->GetPHOSMatrix(mod))
            fGeomPH->SetMisalMatrix(fESDEvent->GetPHOSMatrix(mod),mod);
        else
            fGeomPH->SetMisalMatrix((new TGeoHMatrix),mod);
    }// loop over modules
}

double AliConverterCalorimetersEngine::GetPhi(double phi)
{
    if (phi > TMath::Pi())
        phi -= 2*TMath::Pi();
    return phi;
}
