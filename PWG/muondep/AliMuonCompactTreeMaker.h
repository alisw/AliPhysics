#ifndef ALIMUONCOMPACTTREEMAKER_H
#define ALIMUONCOMPACTTREEMAKER_H

#include <string>
#include "Rtypes.h"
#include "AliMuonCompactEvent.h"

class AliESDEvent;
class AliMUONGeometryTransformer;
class TTree;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

/**

@ingroup pwg_muondep_compact

@class AliMuonCompactTreeMaker

@brief Class to transform regular ESD into a very compact muon one

*/

class AliMuonCompactTreeMaker : public AliAnalysisTaskSE
{
    public:

        virtual void   UserCreateOutputObjects();
        virtual void   UserExec(Option_t *);
        virtual void   NotifyRun();

        AliMuonCompactTreeMaker(const char* ocdbPath);
        AliMuonCompactTreeMaker();

    private:

        void CleanupOCDB();
        Bool_t SetupOCDB(const char* alignSpecificStorage);
        void ConvertEvent(AliESDEvent& esd, AliMuonCompactEvent& compactEvent);
        void GetClusterLocation(Int_t detElemId, Double_t xg, Double_t yg, Double_t zg, Int_t& bendingManuAbsIndex, Int_t& nonBendingManuAbsIndex);

    private:

        AliMuonCompactTreeMaker(const AliMuonCompactTreeMaker&);
        AliMuonCompactTreeMaker& operator=(const AliMuonCompactTreeMaker&);
        std::string fOCDBPath;
        Int_t fRunNumber;
        AliMUONGeometryTransformer* fGeometryTransformer;
        TTree* fOutputTree;
        AliMuonCompactEvent fCompactEvent;

    ClassDef(AliMuonCompactTreeMaker,1)
};

#endif
