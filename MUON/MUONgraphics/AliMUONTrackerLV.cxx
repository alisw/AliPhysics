#include "AliCDBManager.h"
#include "AliLog.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONPadStatusMaker.h"
#include "AliMUONTrackerLV.h"
#include "AliMpCathodType.h"
#include "AliMpDCSNamer.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include "AliMpPlaneType.h"
#include "AliMpSegmentation.h"
#include "AliMpStationType.h"
#include "AliMpVSegmentation.h"
#include "Riostream.h"
#include "TMap.h"
#include "TObjString.h"
#include <set>

///\cond CLASSIMP
ClassImp(AliMUONTrackerLV)
///\endcond

AliMUONTrackerLV::AliMUONTrackerLV(const char* runlist, const char* ocdbPath)
 : AliMUONTrackerVoltages(runlist,ocdbPath)
{
  fOCDBObjectPath = "MUON/Calib/LV";
}

AliMUONTrackerLV::AliMUONTrackerLV(Int_t runNumber, const char* ocdbPath)
: AliMUONTrackerVoltages(runNumber,ocdbPath)
{
  fOCDBObjectPath = "MUON/Calib/LV";
}

void AliMUONTrackerLV::CheckLV(Int_t runNumber, Int_t verbose)
{
    // TMap* m = AliMUONCalibrationData::CreateLV(runNumber);
    AliMp::PlaneType planeType;
    Int_t* detectionElementIds(0x0);
    Int_t numberOfDetectionElements;

    AliMpDCSNamer namer("TRACKER");
    TArrayI manuIds;

    AliMUONCalibrationData cd(runNumber);
    TMap* m = cd.LV();
    TIter next(m);
    TObjString* s;
    AliMUONPadStatusMaker pm(cd);
    std::set<std::string> badLV;
    
    while ( ( s = static_cast<TObjString*>(next()) ) )
    {
        namer.DecodeDCSMCHLVAlias(s->String().Data(),
                detectionElementIds , numberOfDetectionElements, planeType);

        // to check the LV status, we get the first detection element
        // powered by this LV
        // (as all those DEs should have the same LV status)
        Int_t detElemId = detectionElementIds[0];

        AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(
                detElemId);
        const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,
                de->GetCathodType(planeType)); 
        seg->GetAllElectronicCardIDs(manuIds);
        // and from this detection element, the first manu 
        // (as all manus should have the same LV status)
        
        TString plane("");

        if ( pm.LVStatus(detElemId,manuIds[0]) ) {
            TString bad;
            if (de->GetStationType()==AliMp::kStation12)
            {
                plane = AliMp::PlaneTypeName(planeType);
                bad = Form("%04d %4s ",detectionElementIds[0],plane.Data());
            }
            else 
            {
                for ( Int_t i = 0; i < numberOfDetectionElements; ++i )
                {
                    bad += Form("%04d ",detectionElementIds[i]);
                }
            }

            badLV.insert(Form("RUN %06d DE with LV issue :  %s",runNumber,
                        bad.Data()));
        } 
        delete[] detectionElementIds;
    }

    if (!badLV.empty())
    {
        std::cout << std::string(40,'-') << std::endl;
    }
    for ( std::set<std::string>::const_iterator it = badLV.begin();
            it != badLV.end(); ++it )
    {
        std::cout << it->c_str() << std::endl;
    }
}

void AliMUONTrackerLV::Scan(Int_t verbose)
{
    /// Loop over all LV channels and report their status

  if ( fRunList.empty() )
  {
    AliError("No runs to process...");
    return;
  }

  AliCDBManager::Instance()->SetDefaultStorage(fOCDBPath.Data());

  for ( std::vector<int>::size_type i = 0; i < fRunList.size(); ++i )
  {
    CheckLV(fRunList[i],verbose);
  }
}
