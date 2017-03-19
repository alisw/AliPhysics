#include "AliMuonCompactManuStatus.h"

#include "AliAnalysisTriggerScalers.h"
#include "AliCDBManager.h"
#include "AliMUONCDB.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONPadStatusMaker.h"
#include "AliMUONRecoParam.h"
#include "AliMUONRejectList.h"
#include "AliMpCDB.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include "AliMpManuIterator.h"
#include "AliMuonCompactMapping.h"
#include <cassert>

/// \ingroup compact
const UInt_t AliMuonCompactManuStatus::MANUBADPEDMASK = ( 1 << 0 );
const UInt_t AliMuonCompactManuStatus::MANUBADHVMASK = ( 1 << 1 );
const UInt_t AliMuonCompactManuStatus::MANUBADLVMASK =  (1 << 2 );
const UInt_t AliMuonCompactManuStatus::MANUBADOCCMASK = ( 1 << 3 );
const UInt_t AliMuonCompactManuStatus::MANUOUTOFCONFIGMASK = ( 1 << 4 );
const UInt_t AliMuonCompactManuStatus::MANUREJECTMASK = ( 1 << 5 );

std::string AliMuonCompactManuStatus::CauseAsString(UInt_t cause)
{
    std::string rv = "";

    if ( cause & MANUBADPEDMASK ) rv += "_ped";
    if ( cause & MANUBADHVMASK ) rv += "_hv";
    if ( cause & MANUBADLVMASK ) rv += "_lv";
    if ( cause & MANUBADOCCMASK ) rv += "_occ";
    if ( cause & MANUOUTOFCONFIGMASK ) rv += "_config";
    if ( cause & MANUREJECTMASK ) rv += "_reject";

    return rv;
}

void AliMuonCompactManuStatus::Print(const std::vector<UInt_t>& manuStatus, bool all)
{
    AliMuonCompactMapping* cm = AliMuonCompactMapping::GetCompactMapping();

    for ( std::vector<UInt_t>::size_type i = 0; i < manuStatus.size(); ++i )
    {
        Int_t absManuId = cm->AbsManuId(i);
        Int_t detElemId = cm->GetDetElemIdFromAbsManuId(absManuId);
        Int_t manuId = cm->GetManuIdFromAbsManuId(absManuId);
        Int_t busPatchId = AliMpDDLStore::Instance()->GetBusPatchId(detElemId,manuId);
        if ( manuStatus[i] || all )
        {
            std::cout << Form("status[%04lu]=%6x (DE %04d BP %04d MANU %04d) %s",i,manuStatus[i],detElemId,busPatchId,manuId,CauseAsString(manuStatus[i]).c_str()) << std::endl;
        }
    }
}

std::vector<UInt_t> AliMuonCompactManuStatus::BuildFromOCDB(Int_t runNumber, const char* ocdbPath)
{
    std::vector<UInt_t> vManuStatus(16828,0);

    AliCDBManager* man = AliCDBManager::Instance();
    man->SetDefaultStorage(ocdbPath);
    man->SetRun(runNumber);

    if (!AliMpDDLStore::Instance())
    {
        AliMpCDB::LoadAll();
    }
   
    AliMuonCompactMapping* cm = AliMuonCompactMapping::GetCompactMapping(ocdbPath,runNumber);

    AliMUONCalibrationData cd(runNumber,true);

    AliMUONRejectList* rl = cd.RejectList();
    assert(rl->IsBinary());

    AliMUONPadStatusMaker statusMaker(cd);

    AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();

    statusMaker.SetLimits(*recoParam);

    AliMpManuIterator it;
    Int_t detElemId, manuId;

    // Int_t pedCheck = (
    //         AliMUONPadStatusMaker::kPedMeanZero |
    //         AliMUONPadStatusMaker::kPedMeanTooLow |
    //         AliMUONPadStatusMaker::kPedMeanTooHigh |
    //         AliMUONPadStatusMaker::kPedSigmaTooLow |
    //         AliMUONPadStatusMaker::kPedSigmaTooHigh );
    //
    // Int_t hvCheck = (
    //         AliMUONPadStatusMaker::kHVError |
    //         AliMUONPadStatusMaker::kHVTooLow |
    //         AliMUONPadStatusMaker::kHVTooHigh |
    //         AliMUONPadStatusMaker::kHVChannelOFF |
    //         AliMUONPadStatusMaker::kHVSwitchOFF );
    //
    // Int_t occCheck = (
    //         AliMUONPadStatusMaker::kManuOccupancyTooHigh
    //         );
    //
    // Int_t lvCheck = ( AliMUONPadStatusMaker::kLVTooLow );

    Int_t pedCheck = 62; // should really use the enum above, but can't do that until there is an AliRoot tag with that mod, otherwise I break AliPhysics...
    Int_t hvCheck = 31;
    Int_t occCheck = 4;
    Int_t lvCheck = 8;

    Int_t ntotal(0);
    Int_t nbad(0);
    Int_t nbadped=0;
    Int_t nbadocc=0;
    Int_t nbadhv=0;
    Int_t nbadlv=0;
    Int_t nmissing=0;
    Int_t nreco=0;
    Int_t nrejected=0;

    while ( it.Next(detElemId,manuId) )
    {
        AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
        Int_t busPatchId = AliMpDDLStore::Instance()->GetBusPatchId(detElemId,manuId);
        
        UInt_t manuStatus = 0;

        Int_t manubadped=0;
        Int_t manubadocc=0;
        Int_t manubadhv=0;
        Int_t manubadlv=0;
        Int_t manumissing=0;
        Int_t manureject=0;

        for ( Int_t manuChannel = 0; manuChannel < AliMpConstants::ManuNofChannels(); ++manuChannel )
        {
            if ( de->IsConnectedChannel(manuId,manuChannel) )
            {
                ++ntotal;

                UInt_t status = statusMaker.PadStatus(detElemId, manuId, manuChannel);

                //if (!status) continue;

                if ( status & AliMUONPadStatusMaker::BuildStatus(pedCheck,0,0,0) )
                {
                    ++manubadped;
                }

                if ( status & AliMUONPadStatusMaker::BuildStatus(0,hvCheck,0,0) )
                {
                    ++manubadhv;
                }

                if ( status & AliMUONPadStatusMaker::BuildStatus(0,0,lvCheck,0) )
                {
                    ++manubadlv;
                }

                if ( status & AliMUONPadStatusMaker::BuildStatus(0,0,0,occCheck) )
                {
                    ++manubadocc;
                }

                if ( status & AliMUONPadStatusMaker::BuildStatus(128 /*AliMUONPadStatusMaker::kMissing*/,0,0,0) )
                {
                    ++manumissing;
                }

                Float_t proba = TMath::Max(rl->DetectionElementProbability(detElemId),rl->BusPatchProbability(busPatchId));
                proba = TMath::Max(proba,rl->ManuProbability(detElemId,manuId));
                proba = TMath::Max(proba,rl->ChannelProbability(detElemId,manuId,manuChannel));

                if ( proba > 0 )
                {
                    ++manureject;
                }
            }

            if ( manubadped>=0.9*de->NofChannelsInManu(manuId) )
            {
                manuStatus |= MANUBADPEDMASK;
            }

            if ( manubadhv )
            {
                manuStatus |= MANUBADHVMASK;
            }

            if ( manubadlv )
            {
                manuStatus |= MANUBADLVMASK;
            }
            
            if ( manubadocc ) 
            {
                manuStatus |= MANUBADOCCMASK;
            }
            
            if ( manumissing) 
            {
                manuStatus |= MANUOUTOFCONFIGMASK;
            }
            
            if ( manureject >= 0.9*de->NofChannelsInManu(manuId) )
            {
                manuStatus |= MANUREJECTMASK;

            }
            
            Int_t manuAbsIndex = cm->FindManuAbsIndex(detElemId,manuId);
            vManuStatus[manuAbsIndex] = manuStatus;
        }
    }

    assert(ntotal==1064008);

    man->ClearCache();

    return vManuStatus;
}

void AliMuonCompactManuStatus::WriteToBinaryFile(const char* runlist, const char* outputfile, const char* ocdbPath)
{
    AliAnalysisTriggerScalers ts(runlist,ocdbPath);

    std::vector<int> vrunlist = ts.GetRunList(); // FIXME: should need to bring in all the AliAnalysisTriggerScalers class just to read the runlist... 

    std::ofstream out(outputfile,std::ios::binary);

    std::vector<int>::size_type nruns = vrunlist.size();

    out.write((char*)&nruns,sizeof(int));
    out.write((char*)&vrunlist[0],nruns*sizeof(int));
    
    for ( std::vector<int>::size_type i = 0; i < vrunlist.size(); ++i )
    {
        Int_t runNumber = vrunlist[i];
        std::vector<UInt_t> manuStatus = BuildFromOCDB(runNumber,ocdbPath);
        out.write((char*)&manuStatus[0],
                manuStatus.size()*sizeof(int));
        assert(manuStatus.size()==16828);

    std::cout << Form("RUN %6d",runNumber) << std::endl;
    // gObjectTable->Print();

    }
    out.close();
}

void AliMuonCompactManuStatus::ReadManuStatus(const char* inputfile,
    std::map<int,std::vector<UInt_t> >& manuStatusForRuns)
{
    std::ifstream in(inputfile,std::ios::binary);

    int nruns;

    in.read((char*)&nruns,sizeof(int));

    std::cout << "nruns=" << nruns << std::endl;

    std::vector<int> vrunlist;

    vrunlist.resize(nruns,0);

    std::vector<int> manuStatus;

    in.read((char*)&vrunlist[0],sizeof(int)*nruns);

    for ( std::vector<int>::size_type i = 0; i < vrunlist.size(); ++i ) 
    {
        Int_t runNumber = vrunlist[i];
        std::cout << runNumber << " ";
        manuStatus.resize(16828,0);
        in.read((char*)&manuStatus[0],sizeof(int)*manuStatus.size());
        for ( std::vector<int>::size_type j = 0; j < manuStatus.size(); ++j )
        {
            manuStatusForRuns[runNumber].push_back(manuStatus[j]);
        }
    }
    std::cout << std::endl;
}

