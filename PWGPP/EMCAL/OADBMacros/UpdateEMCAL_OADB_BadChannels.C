//!  EMCal bad channel map creation macro
/*!
    With this macro, the 2D histogram bad channel maps can be added to OADB/EMCAL/EMCALBadChannels.root.
*/

#include "TGeoManager.h"

/*******************************************************************
 *  NOTE: Main function which needs to be adjusted for new BC maps *
 *******************************************************************/
void UpdateEMCAL_OADB_BadChannels(const char *fileNameOADBAli="$ALICE_DATA/OADB/EMCAL/EMCALBadChannels.root")
{
    gSystem->Load("libOADB");  
    gSystem->Load("libEMCALbase");
    gSystem->Load("libEMCALUtils");
    gSystem->Load("libEMCALrec");

    // create working copy of OADB file from EOS
    const char *fileNameOADB                ="EMCALBadChannels_temp.root";
    gSystem->Exec(Form("cp %s %s",fileNameOADBAli,fileNameOADB));

    // update OADB file with dead, bad and warm cells
    // last parameter:
    //      0: write new map for given run range
    //      1: update existing map for given run range
    //      2: same as 1 but also save over/underhead parts (special setting, not for normal use!)

    // LHC16f - 20180403
    updateFile(fileNameOADB,"BadChannelsLHC16f","LHC16f/LHC16f_INT7_OADBFile_Dead_V0.txt","LHC16f/LHC16f_INT7_OADBFile_Bad_V0.txt","LHC16f/LHC16f_INT7_OADBFile_Warm_V0.txt",253614,253979,0);

    // LHC16h - 20180405 - Additional bad cells
    updateFile(fileNameOADB,"BadChannelsLHC16h_V1_add","LHC16hijl/empty.txt","LHC16hijl/16h_block1.txt","LHC16hijl/empty.txt",254381,254654,1);
    updateFile(fileNameOADB,"BadChannelsLHC16h_V2_add","LHC16hijl/empty.txt","LHC16hijl/16h_block2.txt","LHC16hijl/empty.txt",254655,255467,1);

    // LHC16i - 20180405 - Additional bad cells
    updateFile(fileNameOADB,"BadChannelsLHC16i_V1_add","LHC16hijl/empty.txt","LHC16hijl/16i_block1.txt","LHC16hijl/empty.txt",255515,255583,1);
    updateFile(fileNameOADB,"BadChannelsLHC16i_V2_add","LHC16hijl/empty.txt","LHC16hijl/16i_block2.txt","LHC16hijl/empty.txt",255591,255618,1);

    // LHC16j - 20180405 - Additional bad cells
    updateFile(fileNameOADB,"BadChannelsLHC16j_V1_add","LHC16hijl/empty.txt","LHC16hijl/16j_block1.txt","LHC16hijl/empty.txt",256147,256281,1);
    updateFile(fileNameOADB,"BadChannelsLHC16j_V2_add","LHC16hijl/empty.txt","LHC16hijl/16j_block2.txt","LHC16hijl/empty.txt",256282,256420,1);

    // LHC16l - 20180405 - Additional bad cells
    updateFile(fileNameOADB,"BadChannelsLHC16l_V1_add","LHC16hijl/empty.txt","LHC16hijl/16l_block1.txt","LHC16hijl/empty.txt",258883,258964,1);
    updateFile(fileNameOADB,"BadChannelsLHC16l_V2_add","LHC16hijl/empty.txt","LHC16hijl/16l_block2.txt","LHC16hijl/empty.txt",259086,259396,1);
    updateFile(fileNameOADB,"BadChannelsLHC16l_V3_add","LHC16hijl/empty.txt","LHC16hijl/16l_block3.txt","LHC16hijl/empty.txt",259469,260014,1);

    // LHC16o - 20180205
    updateFile(fileNameOADB,"BadChannelsLHC16o_V1","LHC16op/LHC16o_INT7_OADBFile_Dead_V1.txt","LHC16op/LHC16o_INT7_OADBFile_Bad_V1.txt","LHC16op/LHC16o_INT7_OADBFile_Warm_V1.txt",262395,262532,0);
    updateFile(fileNameOADB,"BadChannelsLHC16o_V2","LHC16op/LHC16o_INT7_OADBFile_Dead_V2.txt","LHC16op/LHC16o_INT7_OADBFile_Bad_V2.txt","LHC16op/LHC16o_INT7_OADBFile_Warm_V2.txt",262533,262858,0);
    updateFile(fileNameOADB,"BadChannelsLHC16o_V3","LHC16op/LHC16o_INT7_OADBFile_Dead_V3.txt","LHC16op/LHC16o_INT7_OADBFile_Bad_V3.txt","LHC16op/LHC16o_INT7_OADBFile_Warm_V3.txt",263331,263744,0);
    updateFile(fileNameOADB,"BadChannelsLHC16o_V4","LHC16op/LHC16o_INT7_OADBFile_Dead_V4.txt","LHC16op/LHC16o_INT7_OADBFile_Bad_V4.txt","LHC16op/LHC16o_INT7_OADBFile_Warm_V4.txt",263784,264035,0);

    // LHC16p - 20180205
    updateFile(fileNameOADB,"BadChannelsLHC16p_V1","LHC16op/LHC16p_INT7_OADBFile_Dead_V1.txt","LHC16op/LHC16p_INT7_OADBFile_Bad_V1.txt","LHC16op/LHC16p_INT7_OADBFile_Warm_V1.txt",264076,264164,0);
    updateFile(fileNameOADB,"BadChannelsLHC16p_V2","LHC16op/LHC16p_INT7_OADBFile_Dead_V2.txt","LHC16op/LHC16p_INT7_OADBFile_Bad_V2.txt","LHC16op/LHC16p_INT7_OADBFile_Warm_V2.txt",264168,264347,0);

    // LHC16q - 20180405 - Additional bad cells
    updateFile(fileNameOADB,"BadChannelsLHC16q_V1_add","LHC16hijl/empty.txt","LHC16qt/Add_BadCells_LHC16q.txt","LHC16hijl/empty.txt",265305,265525,1);

    // LHC16t - 20180405 - Additional bad cells
    updateFile(fileNameOADB,"BadChannelsLHC16t_V1_add","LHC16hijl/empty.txt","LHC16qt/Add_BadCells_LHC16t.txt","LHC16hijl/empty.txt",267161,267166,1);

    // LHC17g - 20180308
    updateFile(fileNameOADB,"BadChannelsLHC17g_V1","LHC17ghjk/17g_block1/LHC17g_INT7_OADBFile_Dead_V1.txt","LHC17ghjk/17g_block1/LHC17g_INT7_OADBFile_Bad_V1.txt","LHC17ghjk/17g_block1/LHC17g_INT7_OADBFile_Warm_V1.txt",270882,271287,0);
    updateFile(fileNameOADB,"BadChannelsLHC17g_V2","LHC17ghjk/17g_block2/LHC17g_INT7_OADBFile_Dead_V2.txt","LHC17ghjk/17g_block2/LHC17g_INT7_OADBFile_Bad_V2.txt","LHC17ghjk/17g_block2/LHC17g_INT7_OADBFile_Warm_V2.txt",271288,271777,0);

    // LHC17h - 20180308
    updateFile(fileNameOADB,"BadChannelsLHC17h_V1","LHC17ghjk/17h_block1/LHC17h_INT7_OADBFile_Dead_V1.txt","LHC17ghjk/17h_block1/LHC17h_INT7_OADBFile_Bad_V1.txt","LHC17ghjk/17h_block1/LHC17h_INT7_OADBFile_Warm_V1.txt",271868,272074,0);
    updateFile(fileNameOADB,"BadChannelsLHC17h_V2","LHC17ghjk/17h_block2/LHC17h_INT7_OADBFile_Dead_V2.txt","LHC17ghjk/17h_block2/LHC17h_INT7_OADBFile_Bad_V2.txt","LHC17ghjk/17h_block2/LHC17h_INT7_OADBFile_Warm_V2.txt",272075,272122,0);
    updateFile(fileNameOADB,"BadChannelsLHC17h_V3","LHC17ghjk/17h_block3/LHC17h_INT7_OADBFile_Dead_V3.txt","LHC17ghjk/17h_block3/LHC17h_INT7_OADBFile_Bad_V3.txt","LHC17ghjk/17h_block3/LHC17h_INT7_OADBFile_Warm_V3.txt",272123,272123,0);
    updateFile(fileNameOADB,"BadChannelsLHC17h_V4","LHC17ghjk/17h_block4/LHC17h_INT7_OADBFile_Dead_V4.txt","LHC17ghjk/17h_block4/LHC17h_INT7_OADBFile_Bad_V4.txt","LHC17ghjk/17h_block4/LHC17h_INT7_OADBFile_Warm_V4.txt",272124,273101,0);

    // LHC17i - 20180405
    updateFile(fileNameOADB,"BadChannelsLHC17i_V1","LHC17il/17i_block1/LHC17i_INT7_OADBFile_Dead_V1.txt","LHC17il/17i_block1/LHC17i_INT7_OADBFile_Bad_V1.txt","LHC17il/17i_block1/LHC17i_INT7_OADBFile_Warm_V1.txt",273486,273825,0);
    updateFile(fileNameOADB,"BadChannelsLHC17i_V2","LHC17il/17i_block2/LHC17i_INT7_OADBFile_Dead_V2.txt","LHC17il/17i_block2/LHC17i_INT7_OADBFile_Bad_V2.txt","LHC17il/17i_block2/LHC17i_INT7_OADBFile_Warm_V2.txt",273826,274283,0);
    updateFile(fileNameOADB,"BadChannelsLHC17i_V3","LHC17il/17i_block3/LHC17i_INT7_OADBFile_Dead_V3.txt","LHC17il/17i_block3/LHC17i_INT7_OADBFile_Bad_V3.txt","LHC17il/17i_block3/LHC17i_INT7_OADBFile_Warm_V3.txt",274284,274442,0);

    // LHC17j - 20180308
    updateFile(fileNameOADB,"BadChannelsLHC17j_V1","LHC17ghjk/17j_block1/LHC17j_INT7_OADBFile_Dead_V1.txt","LHC17ghjk/17j_block1/LHC17j_INT7_OADBFile_Bad_V1.txt","LHC17ghjk/17j_block1/LHC17j_INT7_OADBFile_Warm_V1.txt",274591,274671,0);

    // LHC17k - 20180308
    updateFile(fileNameOADB,"BadChannelsLHC17k_V1","LHC17ghjk/17k_block1/LHC17k_INT7_OADBFile_Dead_V1.txt","LHC17ghjk/17k_block1/LHC17k_INT7_OADBFile_Bad_V1.txt","LHC17ghjk/17k_block1/LHC17k_INT7_OADBFile_Warm_V1.txt",276290,276508,0);
    updateFile(fileNameOADB,"BadChannelsLHC17k_V2","LHC17ghjk/17k_block2/LHC17k_INT7_OADBFile_Dead_V2.txt","LHC17ghjk/17k_block2/LHC17k_INT7_OADBFile_Bad_V2.txt","LHC17ghjk/17k_block2/LHC17k_INT7_OADBFile_Warm_V2.txt",275515,276289,0);
    updateFile(fileNameOADB,"BadChannelsLHC17k_V3","LHC17ghjk/17k_block3/LHC17k_INT7_OADBFile_Dead_V3.txt","LHC17ghjk/17k_block3/LHC17k_INT7_OADBFile_Bad_V3.txt","LHC17ghjk/17k_block3/LHC17k_INT7_OADBFile_Warm_V3.txt",275057,275514,0);
    updateFile(fileNameOADB,"BadChannelsLHC17k_V4","LHC17ghjk/17k_block4/LHC17k_INT7_OADBFile_Dead_V4.txt","LHC17ghjk/17k_block4/LHC17k_INT7_OADBFile_Bad_V4.txt","LHC17ghjk/17k_block4/LHC17k_INT7_OADBFile_Warm_V4.txt",274690,275056,0);

    // LHC17l - 20180405
    updateFile(fileNameOADB,"BadChannelsLHC17l_V1","LHC17il/17l_block1/LHC17l_INT7_OADBFile_Dead_V1.txt","LHC17il/17l_block1/LHC17l_INT7_OADBFile_Bad_V1.txt","LHC17il/17l_block1/LHC17l_INT7_OADBFile_Warm_V1.txt",276551,277117,0);
    updateFile(fileNameOADB,"BadChannelsLHC17l_V2","LHC17il/17l_block2/LHC17l_INT7_OADBFile_Dead_V2.txt","LHC17il/17l_block2/LHC17l_INT7_OADBFile_Bad_V2.txt","LHC17il/17l_block2/LHC17l_INT7_OADBFile_Warm_V2.txt",277118,277389,0);
    updateFile(fileNameOADB,"BadChannelsLHC17l_V3","LHC17il/17l_block3/LHC17l_INT7_OADBFile_Dead_V3.txt","LHC17il/17l_block3/LHC17l_INT7_OADBFile_Bad_V3.txt","LHC17il/17l_block3/LHC17l_INT7_OADBFile_Warm_V3.txt",277390,277721,0);
    updateFile(fileNameOADB,"BadChannelsLHC17l_V4","LHC17il/17l_block4/LHC17l_INT7_OADBFile_Dead_V4.txt","LHC17il/17l_block4/LHC17l_INT7_OADBFile_Bad_V4.txt","LHC17il/17l_block4/LHC17l_INT7_OADBFile_Warm_V4.txt",277722,278216,0);

    // LHC17m - 20180305
    updateFile(fileNameOADB,"BadChannelsLHC17m_V1","LHC17m/17m_block1/LHC17m_INT7_OADBFile_Dead_V1.txt","LHC17m/17m_block1/LHC17m_INT7_OADBFile_Bad_V1.txt","LHC17m/17m_block1/LHC17m_INT7_OADBFile_Warm_V1.txt",278818, 279355,0);
    updateFile(fileNameOADB,"BadChannelsLHC17m_V2","LHC17m/17m_block2/LHC17m_INT7_OADBFile_Dead_V2.txt","LHC17m/17m_block2/LHC17m_INT7_OADBFile_Bad_V2.txt","LHC17m/17m_block2/LHC17m_INT7_OADBFile_Warm_V2.txt",279356, 280140,0);

    // LHC17o - 20180206
    updateFile(fileNameOADB,"BadChannelsLHC17o_V1","LHC17o/17o_block1/LHC17o_INT7_OADBFile_Dead_V1.txt","LHC17o/17o_block1/LHC17o_INT7_OADBFile_Bad_V1.txt","LHC17o/17o_block1/LHC17o_INT7_OADBFile_Warm_V1.txt",280282, 280403,0);
    updateFile(fileNameOADB,"BadChannelsLHC17o_V2","LHC17o/17o_block2/LHC17o_INT7_OADBFile_Dead_V2.txt","LHC17o/17o_block2/LHC17o_INT7_OADBFile_Bad_V2.txt","LHC17o/17o_block2/LHC17o_INT7_OADBFile_Warm_V2.txt",280404, 280999,0);
    updateFile(fileNameOADB,"BadChannelsLHC17o_V3","LHC17o/17o_block3/LHC17o_INT7_OADBFile_Dead_V3.txt","LHC17o/17o_block3/LHC17o_INT7_OADBFile_Bad_V3.txt","LHC17o/17o_block3/LHC17o_INT7_OADBFile_Warm_V3.txt",281000, 281753,0);
    updateFile(fileNameOADB,"BadChannelsLHC17o_V4","LHC17o/17o_block4/LHC17o_INT7_OADBFile_Dead_V4.txt","LHC17o/17o_block4/LHC17o_INT7_OADBFile_Bad_V4.txt","LHC17o/17o_block4/LHC17o_INT7_OADBFile_Warm_V4.txt",281757, 281961,0);
    updateFile(fileNameOADB,"BadChannelsLHC17o_V5","LHC17o/17o_block5/LHC17o_INT7_OADBFile_Dead_V5.txt","LHC17o/17o_block5/LHC17o_INT7_OADBFile_Bad_V5.txt","LHC17o/17o_block5/LHC17o_INT7_OADBFile_Warm_V5.txt",281754, 281756,0);

    // LHC17p - 20180206
    updateFile(fileNameOADB,"BadChannelsLHC17p_V1","LHC17pq/LHC17p_INT7_Bad_Ampl_block1_dead.txt","LHC17pq/LHC17p_INT7_Bad_Ampl_block1_bad.txt","LHC17pq/LHC17p_INT7_Bad_Ampl_block1_warm.txt",282025,282126,0);
    updateFile(fileNameOADB,"BadChannelsLHC17p_V2","LHC17pq/LHC17p_INT7_Bad_Ampl_block2_dead.txt","LHC17pq/LHC17p_INT7_Bad_Ampl_block2_bad.txt","LHC17pq/LHC17p_INT7_Bad_Ampl_block2_warm.txt",282127,282147,0);
    updateFile(fileNameOADB,"BadChannelsLHC17p_V3","LHC17pq/LHC17p_INT7_Bad_Ampl_block3_dead.txt","LHC17pq/LHC17p_INT7_Bad_Ampl_block3_bad.txt","LHC17pq/LHC17p_INT7_Bad_Ampl_block3_warm.txt",282148,282314,0);
    updateFile(fileNameOADB,"BadChannelsLHC17p_V4","LHC17pq/LHC17p_INT7_Bad_Ampl_block4_dead.txt","LHC17pq/LHC17p_INT7_Bad_Ampl_block4_bad.txt","LHC17pq/LHC17p_INT7_Bad_Ampl_block4_warm.txt",282315,282343,0);

    // LHC17q - 20180206
    updateFile(fileNameOADB,"BadChannelsLHC17q_V1","LHC17pq/LHC17q_INT7_Bad_Ampl_block1_dead.txt","LHC17pq/LHC17q_INT7_Bad_Ampl_block1_bad.txt","LHC17pq/LHC17q_INT7_Bad_Ampl_block1_warm.txt",282365,282366,0);
    updateFile(fileNameOADB,"BadChannelsLHC17q_V2","LHC17pq/LHC17q_INT7_Bad_Ampl_block2_dead.txt","LHC17pq/LHC17q_INT7_Bad_Ampl_block2_bad.txt","LHC17pq/LHC17q_INT7_Bad_Ampl_block2_warm.txt",282367,282440,0);



    // the final output will be sorted by runnumber
    sortOutput(fileNameOADB);
    // the OADB container is additionally checked for ownership
    rebuildBadCellContainer("EMCALBadChannelsNEW.root");
}

/*******************************************************************
 *  NOTE: Update function that removes existing BC maps from old   *
 *          file and adds new BC maps at the end of the file       *
 *  NOTE: The final file is saved and directly renamed as new      *
 *      input file for the next loop (this reduces total file size)*
 *******************************************************************/
void updateFile(const char *fileNameOADB,TString arrName, TString fileNameDeadCells, TString fileNameBadCells, TString fileNameWarmCells,Int_t lowRun, Int_t highRun, Int_t updateExistingMap = 0){
    printf("\n\nAdding %s maps to OADB file:\n",arrName.Data());

    AliEMCALGeometry   *fEMCALGeo               = new AliEMCALGeometry();
    AliEMCALRecoUtils  *fEMCALRecoUtils         = new AliEMCALRecoUtils();

    fEMCALGeo                                   = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM");
    const int nSM=20;

    // load OADB file
    TFile *f                                    = TFile::Open(fileNameOADB);
    f->ls();
    AliOADBContainer *con                       =(AliOADBContainer*)f->Get("AliEMCALBadChannels");
    con->SetName("Old"); 

    // make the output container
    AliOADBContainer *con2                      = new AliOADBContainer("AliEMCALBadChannels");

    // List of brand new arrays (including overhead arrays for updateExistingMap == 2)
    TObjArray arrayAdd(nSM);
    arrayAdd.SetName(arrName.Data());
    TObjArray arrayAddOverheadLow(nSM);
    arrayAddOverheadLow.SetName(Form("%s_OHLow",arrName.Data()));
    TObjArray arrayAddOverheadHigh(nSM);
    arrayAddOverheadHigh.SetName(Form("%s_OHHigh",arrName.Data()));

    // load OADB file to have access to existing maps
    TFile *foadb                                = TFile::Open(fileNameOADB);
    AliOADBContainer *conAP                     = (AliOADBContainer*)foadb->Get("AliEMCALBadChannels");
    TObjArray *arrayBClow                       = (TObjArray*)conAP->GetObject(lowRun);
    TObjArray *arrayBChigh                      = (TObjArray*)conAP->GetObject(highRun);

    TH2I *hBadChannels[nSM];
    TH2I *hBadChannelsOverhead[nSM];
    char nameSM[50];

    // initialize histograms in array - either load existing map or create new, empty map
    for(int i=0;i<nSM;i++){
        sprintf(nameSM,"EMCALBadChannelMap_Mod%d",i);
        hBadChannels[i]                              = NULL;
        hBadChannelsOverhead[i]                      = NULL;
        if(updateExistingMap){
            cout << "Initializing " << nameSM << " from current AliPhysics OADB file ... ";
            if(arrayBClow){
                hBadChannels[i]                      = (TH2I*)arrayBClow->FindObject(nameSM);
                if(hBadChannels[i] && updateExistingMap>1)
                    hBadChannelsOverhead[i]          = (TH2I*)hBadChannels[i]->Clone(nameSM);
            } else if (arrayBChigh){
                hBadChannels[i]                      = (TH2I*)arrayBChigh->FindObject(nameSM);
                if(hBadChannels[i] && updateExistingMap>1)
                    hBadChannelsOverhead[i]          = (TH2I*)hBadChannels[i]->Clone(nameSM);
            }
            if(hBadChannels[i])
                cout << "and found it!" << endl;
        }
        if(!updateExistingMap || !hBadChannels[i]){
            cout << "File not found! Creating new map for " << nameSM << endl;
            hBadChannels[i]                          = new TH2I(nameSM,nameSM,48,0,48,24,0,24);
        }
    }

    // Read in the list of dead channels from the dead cell file
    std::vector<TString> vecDeadChannels;
    if(!readin(fileNameDeadCells, vecDeadChannels,kTRUE,0)) {
        cout << "No dead Channels could be found in file: " << fileNameDeadCells.Data() << endl;
    }
    // Read in the list of bad channels from the bad cell file
    std::vector<TString> vecBadChannels;
    if(!readin(fileNameBadCells, vecBadChannels,kTRUE,1)) {
        cout << "No bad Channels could be found in file: " << fileNameBadCells.Data() << endl;
    }
    // Read in the list of warm channels from the warm cell file
    std::vector<TString> vecWarmChannels;
    if(!readin(fileNameWarmCells, vecWarmChannels,kTRUE,2)) {
        cout << "No warm Channels could be found in file: " << fileNameWarmCells.Data() << endl;
    }

    // Fill the bad channel map with different values depending on dead, bad and warm cells
    Fill_histo(vecDeadChannels,1,hBadChannels,fEMCALGeo,fEMCALRecoUtils);
    Fill_histo(vecBadChannels,2,hBadChannels,fEMCALGeo,fEMCALRecoUtils);
    Fill_histo(vecWarmChannels,3,hBadChannels,fEMCALGeo,fEMCALRecoUtils);

    // add histograms to a new array
    for (Int_t mod=0;mod<nSM;mod++){
        arrayAdd.Add(hBadChannels[mod]);
        if( updateExistingMap>1 && hBadChannelsOverhead[mod]){
            arrayAddOverheadLow.Add(hBadChannelsOverhead[mod]);
            arrayAddOverheadHigh.Add(hBadChannelsOverhead[mod]);
        }
    }

    // check old OADB file for existing BC maps in given run range and remove that old BC map
    Int_t replacedContainerLowLimit             = -1;
    Int_t replacedContainerHighLimit            = -1;
    for(int i=0;i<con->GetNumberOfEntries();i++){
        if (lowRun >= con->LowerLimit(i) && lowRun <= con->UpperLimit(i)){
            printf("\n!!! Not adding index %d for runrange %d--%d as low run limit %d is contained\n", con->GetObjectByIndex(i), con->LowerLimit(i),con->UpperLimit(i),lowRun);
            replacedContainerLowLimit           = con->LowerLimit(i);
            replacedContainerHighLimit          = con->UpperLimit(i);
        } else if (highRun >= con->LowerLimit(i) && highRun <= con->UpperLimit(i)){
            printf("\n!!! Not adding index %d for runrange %d--%d as high run limit %d is contained\n", con->GetObjectByIndex(i), con->LowerLimit(i),con->UpperLimit(i),highRun);
            replacedContainerLowLimit           = con->LowerLimit(i);
            replacedContainerHighLimit          = con->UpperLimit(i);
        }else{
            con2->AddDefaultObject(con->GetObjectByIndex(i));
            con2->AppendObject(con->GetObjectByIndex(i),con->LowerLimit(i),con->UpperLimit(i));
        }
    }
    // add new BC maps at the end of the new OADB file (will be sorted later)
    if( updateExistingMap>1 && replacedContainerLowLimit>0 && lowRun > replacedContainerLowLimit){
        con2->AddDefaultObject(&arrayAddOverheadLow);
        con2->AppendObject(&arrayAddOverheadLow, replacedContainerLowLimit,lowRun-1);
    }
    con2->AddDefaultObject(&arrayAdd);
    con2->AppendObject(&arrayAdd, lowRun,highRun);
    if(updateExistingMap>1 && replacedContainerHighLimit>0 && highRun < replacedContainerHighLimit){
        con2->AddDefaultObject(&arrayAddOverheadHigh);
        con2->AppendObject(&arrayAddOverheadHigh, highRun+1,replacedContainerHighLimit);
    }

    // temporarilty save BC map file and rename as new input file
    con2->WriteToFile("tempBC.root");
    gSystem->Exec(Form("mv tempBC.root %s",fileNameOADB));

    printf("\nMaps have been successfully added!\n\n",arrName.Data());

}

/*******************************************************************
 *  NOTE: Function to fill the 2D bad channel histograms           *
 *             using EMCAL Geo for cellID->Eta/Phi                 *
 *******************************************************************/
void Fill_histo(std::vector<TString> vecChannels,Int_t fill=1,TH2I **hBadChannels,AliEMCALGeometry   *fEMCALGeo,AliEMCALRecoUtils  *fEMCALRecoUtils){
    Int_t nSupMod=-1, nModule=-1, nIphi=-1, nIeta=-1, iphi=-1, ieta=-1;
    for(Int_t i=0; i<(Int_t) vecChannels.size(); i++){
        fEMCALGeo->GetCellIndex(vecChannels.at(i).Atoi(),nSupMod,nModule,nIphi,nIeta);
        fEMCALGeo->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta,iphi,ieta);
        fEMCALRecoUtils->SetEMCALChannelStatus(nSupMod, ieta, iphi);
        //cout<<endl<<"ID:"<<vecChannels.at(i).Atoi()<<" NSupMod: "<<nSupMod<<" ieta: "<<ieta<<" iphi:"<<iphi<<" fill:"<<fill<<endl; 
        hBadChannels[nSupMod]->SetBinContent(ieta,iphi,fill);
    }
}

/*******************************************************************
 *  NOTE: Sorting function to sort the final OADB file             *
 *                  by ascending runnumber                         *
 *******************************************************************/
void sortOutput(const char *fileNameOADB=""){

    TFile *f                                    = TFile::Open(fileNameOADB);
    f->ls();
    AliOADBContainer *con                       =(AliOADBContainer*)f->Get("AliEMCALBadChannels");
    con->SetName("Old"); 

    Int_t indexAdd                              = 0;
    Int_t largerthan                            = 0;
    Int_t currentvalue                          = 0;

    AliOADBContainer *con2                      = new AliOADBContainer("AliEMCALBadChannels");
    con2->SetOwner();
    // First entry needs to be added before sorting loop
    con2->AppendObject(con->GetObjectByIndex(0),con->LowerLimit(0),con->UpperLimit(0));
    TString strTemp                             = "";

    // sorting magic happens here
    for(int i=1;i<con->GetNumberOfEntries();i++){
        largerthan                              = con2->UpperLimit(con2->GetNumberOfEntries()-1);
        currentvalue                            = -1;
        indexAdd                                = 0;
        for(int j=0;j<con->GetNumberOfEntries();j++){
            if(con->UpperLimit(j)<=largerthan) 
                continue;
            if(currentvalue < 0){
                currentvalue                    = con->UpperLimit(j);
                indexAdd                        = j;
            }
            if(con->UpperLimit(j)<currentvalue){
                currentvalue                    = con->UpperLimit(j);
                indexAdd                        = j;
            }
        }
        con2->AppendObject(con->GetObjectByIndex(indexAdd),con->LowerLimit(indexAdd),con->UpperLimit(indexAdd));
    }

    printf("\n\n");
    Int_t nentries2                             = con2->GetNumberOfEntries();
    for(int i=0;i<nentries2;i++){
        printf("\n Entry2 --> %d/%d -->",i,nentries2);
        printf("%d -- %d --> obj = %p , %s", con2->LowerLimit(i),con2->UpperLimit(i),con2->GetObjectByIndex(i),con2->GetObjectByIndex(i)->GetName());
    }
    printf("\n\n");

    con2->WriteToFile("EMCALBadChannelsNEW.root");
}

/*******************************************************************
 *  NOTE: Function to read the bad channels from                   *
 *                  the given input file                           *
 *******************************************************************/
Bool_t readin(TString fileRuns, std::vector<TString> &vec, Bool_t output, Int_t cellType)
{
    if(output) cout << Form("\nReading from %s...", fileRuns.Data()) << endl;
    fstream file;
    TString fVar;
    Int_t totalN=0;
    file.open(fileRuns.Data(), ios::in);
    if(file.good())
    {
        file.seekg(0L, ios::beg);
        if(output){
            switch(cellType){
                case 0:
                    cout << "Read in dead cells: \"";
                    break;
                case 1:
                    cout << "Read in bad cells: \"";
                    break;
                case 2:
                    cout << "Read in warm cells: \"";
                    break;
                default:
                    cout << "Read in cells (no type set): \"";
            }
        }
        while(!file.eof())
        {
            file >> fVar;
            if(fVar.Sizeof()>1)
            {
                if(output) cout << fVar.Data() << ", ";
                vec.push_back(fVar);
                totalN++;
            }
        }
        if(output) cout << "\"" << endl;
    }
    file.close();
    if(totalN>0){
        switch(cellType){
                case 0:
                    cout << "...done!\n\nIn total " << totalN << " dead cells were read in!\n" << endl;
                    break;
                case 1:
                    cout << "...done!\n\nIn total " << totalN << " bad cells were read in!\n" << endl;
                    break;
                case 2:
                    cout << "...done!\n\nIn total " << totalN << " warm cells were read in!\n" << endl;
                    break;
                default:
                    cout << "...done!\n\nIn total " << totalN << " cells were read in!\n" << endl;
            }
        return kTRUE;
    }
    else return kFALSE;
}



TObjArray *CreatePeriodContainer(TObjArray *inputcont){
  TObjArray *newcont = new TObjArray(inputcont->GetEntries());
  newcont->SetName(inputcont->GetName());
  for(int i = 0; i < inputcont->GetEntries(); i++){
    newcont->AddAt(inputcont->At(i)->Clone(), i);
  }
  return newcont;
}

/*******************************************************************
 *                                                                 *
 *        NOTE: Function required to fix OADB ownership            *
 *                                                                 *
 *******************************************************************/
void rebuildBadCellContainer(const char *fileNameOADB=""){
  TFile *reader = TFile::Open(fileNameOADB);
  AliOADBContainer *cont = static_cast<AliOADBContainer *>(reader->Get("AliEMCALBadChannels"));
  delete reader;

  AliOADBContainer *newcont = new AliOADBContainer("AliEMCALBadChannels");
  for(int irun = 0; irun < cont->GetNumberOfEntries(); irun++){
    newcont->AppendObject(CreatePeriodContainer(static_cast<TObjArray *>(cont->GetObjArray()->At(irun))), cont->LowerLimit(irun), cont->UpperLimit(irun));
  }

  newcont->WriteToFile("EMCALBadChannelsOADBfix.root");

  TFile *reader = TFile::Open("EMCALBadChannelsOADBfix.root", "READ");
    AliOADBContainer *cont = static_cast<AliOADBContainer *>(reader->Get("AliEMCALBadChannels"));
    delete reader;
    delete cont;
}
