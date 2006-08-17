/************************************************************************
**
** This file is property of and copyright by the Computer Science/Computer 
** Engineering Group, Kirchhoff Institute for Physics, Ruprecht-Karls-
** University, Heidelberg, Germany, 2005
** This file has been written by Jochen Thaeder, 
** thaeder@kip.uni-heidelberg.de
**
**
** See the file license.txt for details regarding usage, modification,
** distribution and warranty.
** Important: This file is provided without any warranty, including
** fitness for any particular purpose.
**
**
** Newer versions of this file's package will be made available from 
** http://www.kip.uni-heidelberg.de/ti/HLT/
** or the corresponding page of the Heidelberg Alice HLT group.
**
*************************************************************************/

// *****************************************************************
// *        USEAGE in ROOT
// *****************************************************************
// * -->LOAD the macro first:
// * .L ("HLT-OnlineDisplay-HOMER.C");
// *
// * --> use with ".p FUNCTION;" OR here as macro function
// * ODH_Init();
// * ODH_Connect(Char_t *hostname1,Int_t port1,Char_t *hostname2,Int_t port2,Char_t *hostname3,Int_t port3);
// * ODH_DisplayEvent();                                    // Display NextEvent
// *
// * ODH_SetSliceRange();                                   // Sets Slice: ALL
// * ODH_SetSliceRange(Int_t slice);                        // Sets Slice: "slice"
// * ODH_SetSliceRange(Int_t minslice, Int_t maxslice);     // Sets Slice: "minslice" to  "maxslice"
// * ODH_SetSlicePair(Int_t slice);                         // Sets Slice: "slice" and  "slice"+9
// * ODH_SetSlicePair(Int_t minslice, Int_t maxslice);      // Sets Slice: ["minslice" to  "maxslice"] and  ["minslice"+9 to  "maxslice"+9]
// *
// * ODH_SetCluster(Bool_t used, Bool_t unused);            // Sets Cluster: used, unused, all(= used && unused)
// * 
// * ODH_SetInvert();                                       // Invert 3D Display
// * ODH_SetKeepView(Bool_t keepview);                      // Keep 3D View
// * 
// * ODH_SetPadRow(Int_t slice, Int_t padrow, Int_t pad);   // Set PadRow "padrow" within Slice "slice" and Pad "pad"(=0 per default)  
// *
// * ODH_SetSelectTrack(Bool_t switch, Int_t slice, Int_t track); // with switch=true :Select Single Track "track" in Slice "slice"
// *                                                              // turn off: ODH_SetSelectTrack();
// *
// * ODH_SetTrack(Int_t minhits, Float_t ptthreshold);      // Sets Cuts in order to display tracks, both default is 0
// * 
// * ODH_Set3D(Bool_t tracks, Bool_t cluster, Bool_t padrow, Bool_t geometry);  // Switches for the 3D Display
// * ODH_Set3DTracks(Bool_t on);                            // Set 3D tracks, default is off
// * ODH_Set3DCluster(Bool_t on);                           // Set 3D cluster, default is off
// * ODH_Set3DPadRow(Bool_t on);                            // Set 3D padrow, default is off
// * ODH_Set3DGeometry(Bool_t on);                          // Set 3D geometry, default is off
// *
// *****************************************************************
{

    gROOT->LoadMacro("HLT-OnlineDisplay-HOMER.C");

    // --- LOAD DISPLAY FUNCTIONS
    HLT_OnlineDisplay_HOMER();
    
    // --- INITIALIZE DISPLAY, LOAD GEOMETRY FILE, LOAD LIBRARIES (HOMER,HLT,..)
    ODH_Init();

    //  --- CONNECT TO TCP DUMP SUBSCRIBERS
    ODH_Connect("e300",42002,NULL,NULL,"e300",42001);   

    // --- Next Event
    ODH_DisplayEvent();

    barCl = new TControlBar("vertical", "HLT DISPLAY - Cluster");
    barCl->AddButton("All Cluster",".p ODH_SetCluster(true,true)", "All Cluster");
    barCl->AddButton("Used Cluster",".p ODH_SetCluster(true,false)", "Used Cluster");
    barCl->AddButton("Unused Cluster",".p ODH_SetCluster(false,true)", "Unused Cluster");
    barCl->Show();
    
    bar = new TControlBar("vertical", "HLT DISPLAY");
    bar->AddButton("Next Event",".p ODH_DisplayEvent()", "Next Event");
    bar->AddSeparator();
    bar->AddButton("Show all slices",".p ODH_SetSliceRange()", "Show all slices");
    bar->AddSeparator();
    bar->AddButton("Keep 3D View",".p ODH_SetKeepView(true)", "Keep 3D View");
    bar->AddButton("!Keep 3D View",".p ODH_SetKeepView(false)", "Keep 3D View");
    bar->AddButton("Invert 3D View",".p ODH_SetInvert()","Invert 3D View");
    bar->AddSeparator();
    bar->AddButton("Show 3D Tracks",".p ODH_Set3DTracks(true)","Show 3D Tracks");
    bar->AddButton("!Show 3D Tracks",".p ODH_Set3DTracks()","!Show 3D Tracks");
    bar->AddButton("Show 3D Cluster",".p ODH_Set3DCluster(true)","Show 3D Cluster");
    bar->AddButton("!Show 3D Cluster",".p ODH_Set3DCluster()","!Show 3D Cluster");
    bar->AddButton("Show 3D PadRow",".p ODH_Set3DPadRow(true)","Show 3D PadRow");
    bar->AddButton("!Show 3D PadRow",".p ODH_Set3DPadRow()","!Show 3D PadRow");
    bar->AddButton("Show 3D Geometry",".p ODH_Set3DGeometry(true)","Show 3D Geometry");
    bar->AddButton("!Show 3D Geometry",".p ODH_Set3DGeometry()","!Show 3D Geometry");
    bar->AddSeparator();
    bar->AddButton("Show Sector 0",".p ODH_SetSliceRange(0)", "Show Next Event");
    bar->AddButton("Show Sector 9",".p ODH_SetSliceRange(9)", "Show Next Event");
    bar->AddButton("Show Sector 18",".p ODH_SetSliceRange(18)", "Show Next Event");
    bar->AddButton("Show Sector 27",".p ODH_SetSliceRange(27)", "Show Next Event");
    bar->AddButton("Show Sector 0 - 2",".p ODH_SetSliceRange(0,2)", "Show Next Event");
    bar->AddButton("Show Sector 0 - 9",".p ODH_SetSliceRange(0,9)", "Show Next Event");
    bar->AddButton("Show pair 0",".p ODH_SetSlicePair(0)", "Show Next Event");
    bar->AddButton("Show pair 0 - 2",".p ODH_SetSlicePair(0,2)", "Show Next Event");
    
    bar->AddSeparator();

bar->AddButton("Close","gROOT.Reset(\"a\")", "Close");
bar->Show();
gROOT->SaveContext();
}


