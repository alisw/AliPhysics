// --- Custom header files ---
#include <AliPP13TriggerProperties.h>

// --- ROOT system ---

// --- AliRoot header files ---
#include <AliPHOSGeometry.h>

ClassImp(AliPP13TriggerProperties);

//________________________________________________________________
void AliPP13TriggerProperties::FillTriggerInformation(AliPP13AnalysisCluster * cluster)
{
    // Coordinates of a cluster
    Float_t  position[3];
    cluster->GetPosition(position);
    TVector3 global(position);

    Int_t relid[4];
    AliPHOSGeometry::GetInstance()->GlobalPos2RelId(global, relid);
    Int_t tru = TRU(relid[2], relid[3]);

    Int_t tru_ch_x, tru_ch_z;
    Int_t tru_ch = TRUChannel(relid[2], relid[3], tru_ch_x, tru_ch_z);

    cluster->SetTrigger(kFALSE);
    cluster->SetTRU(tru);
    cluster->SetTRUChannel(tru_ch, tru_ch_x, tru_ch_z);
    cluster->SetModule(relid[0]);

    fTrigger->Reset();
    while (fTrigger->Next())
    {
        if (fTrigger->GetL1TimeSum() != fL1Threshold)
            continue;

        if (Matched(fTrigger, relid))
        {
            cluster->SetTrigger(kTRUE);
            return;
        }
    }
}

//________________________________________________________________
Bool_t AliPP13TriggerProperties::Matched(AliVCaloTrigger * trigger, Int_t * relid)
{
    // "Online" module number, bottom-left 4x4 edge cell absId
    // Get the global position of the fired trigger
    Int_t trigger_module, trigger_absId;
    trigger->GetPosition(trigger_module, trigger_absId);

    // Convert it to local coordinates of PHOS
    Int_t trigger_relid[4] ;
    AliPHOSGeometry::GetInstance()->AbsToRelNumbering(trigger_absId, trigger_relid);

    // Returns kTRUE if cluster position coincides with 4x4 patch.
    //
    if (trigger_relid[0] != relid[0]) // different modules!
        return kFALSE;

    if (TMath::Abs(trigger_relid[2] - relid[2]) > 3) // X-distance too large!
        return kFALSE;

    if (TMath::Abs(trigger_relid[3] - relid[3]) > 3) // Z-distance too large!
        return kFALSE;

    return kTRUE;
}

//________________________________________________________________________
Int_t AliPP13TriggerProperties::TRU(Int_t cellx, Int_t cellz)
{
    Int_t tru = -1;
    if (cellx < 1 || 64 < cellx) 
        return -1;

    if (cellz < 1 || 56 < cellz)
        return -1;

    Int_t XID = (cellx - 1) / 16 + 1;
    Int_t ZID = (cellz - 1) / 28;

    tru = 2 * XID - ZID;
    return tru;
}
//________________________________________________________________________
Int_t AliPP13TriggerProperties::TRUChannel(Int_t cellx, Int_t cellz, Int_t &chX, Int_t &chZ)
{
    // Returns TRU channel in the range 0-111.
    Int_t ch = -1;
    if (cellx < 1 || 64 < cellx)
        return -1;

    if (cellz < 1 || 56 < cellz)
        return -1;

    chX = ((cellx - 1) / 2) %  8;
    chZ = ((cellz - 1) / 2) % 14;

    ch = 8 * chZ + chX;
    return ch;
}
