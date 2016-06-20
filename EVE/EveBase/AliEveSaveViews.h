//
//  AliEveSaveViews.h
//
//
//  Created by Jeremi Niedziela on 16/03/15.
//
//

#ifndef __AliRoot__AliEveSaveViews__
#define __AliRoot__AliEveSaveViews__

#include <AliESDEvent.h>

#include <TString.h>
#include <TASImage.h>

#include <vector>


class AliEveSaveViews
{
public:
    AliEveSaveViews(int width=1440,int height=900,int heightInfoBar=65,bool showLiveBar=true);
    ~AliEveSaveViews();
    
    void SaveForAmore();
    void Save(bool withDialog=true,char* filename="");
    int SendToAmore();
private:
    void ChangeRun();
    
    void BuildEventInfoString();
    void BuildTriggerClassesStrings();
    void BuildClustersInfoString();
    
    int fRunNumber;
    int fCurrentFileNumber;
    int fMaxFiles;
    bool fShowLiveBar;
    int fNumberOfClusters;
    TString fCompositeImgFileName;
    
    AliESDEvent *fESDEvent;
    
    // images and dimensions
    int fHeightInfoBar;     // height of the Information bar
    int fWidth;             // width of resulting image
    int fHeight;            // height of resulting image
        
    // strings with additional info
    TString fEventInfo;
    TString fTriggerClasses[3];
    TString fClustersInfo;
    
    // info from logbook
    std::vector<int> fCluster;
    std::vector<ULong64_t> fInputDetectorMask;
    
    AliEveSaveViews(const AliEveSaveViews&);
    AliEveSaveViews& operator=(const AliEveSaveViews&);
};

#endif
