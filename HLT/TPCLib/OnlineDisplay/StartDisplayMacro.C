{

//gROOT->Reset();

gROOT->LoadMacro("HLT-OnlineDisplay-HOMER.C");

HLT_OnlineDisplay_HOMER();
//ODH_Init("./GEO","/home/HLT/src/versions/devel/Util/HOMER/reader/lib/Linux-i686/");

ODH_Init("./GEO","/home/jthaeder/HLT/AnalysisChain/lib/Linux-i686/");
//ODH_Init("./GEO","/home/jthaeder/HLT/src/lib_jthaeder/");

//ODH_Connect("e107",42002,NULL,NULL,"e107",42003);     // pprun_1_nB
//ODH_Connect(NULL,NULL,NULL,NULL,"e107",42003);


ODH_Connect(NULL,NULL,"eh000",42002,NULL,NULL); // --   HUGE_TEST with RAW
//ODH_DisplayNextEvent(false,false,true,NULL);

bar = new TControlBar("vertical", "HLT DISPLAY");

bar->AddButton("Next clusters",".p ODH_DisplayNextEvent(true,false,false,NULL)", "Show Next Event");
bar->AddButton("Next tracks",".p ODH_DisplayNextEvent(false,true,false,NULL)", "Show Next Event");
bar->AddButton("Next clusters and tracks",".p ODH_DisplayNextEvent(true,true,false,NULL)", "Show Next Event");

bar->AddSeparator();

bar->AddButton("Show all",".p ODH_SetSliceRange()", "Show Next Event");

bar->AddButton("Show Sector 0",".p ODH_SetSliceRange(0)", "Show Next Event");
bar->AddButton("Show Sector 1",".p ODH_SetSliceRange(1)", "Show Next Event");
bar->AddButton("Show Sector 2",".p ODH_SetSliceRange(2)", "Show Next Event");
bar->AddButton("Show Sector 9",".p ODH_SetSliceRange(9)", "Show Next Event");
bar->AddButton("Show Sector 10",".p ODH_SetSliceRange(10)", "Show Next Event");
bar->AddButton("Show Sector 11",".p ODH_SetSliceRange(11)", "Show Next Event");

bar->AddButton("Show Sector 0 - 2",".p ODH_SetSliceRange(0,2)", "Show Next Event");
bar->AddButton("Show Sector 9 - 11",".p ODH_SetSliceRange(9,11)", "Show Next Event");
bar->AddButton("Show Sector 0 - 11",".p ODH_SetSliceRange(0,11)", "Show Next Event");

bar->AddButton("Show pair 0",".p ODH_SetSlicePair(0)", "Show Next Event");
bar->AddButton("Show pair 1",".p ODH_SetSlicePair(1)", "Show Next Event");
bar->AddButton("Show pair 2",".p ODH_SetSlicePair(2)", "Show Next Event");

bar->AddButton("Show Geometry",".p ODH_SetDrawGeo()", "Show Next Event");

bar->AddButton("Invert",".p ODH_SetInvert()","ccc");

bar->AddSeparator();
// PADROW 0 - 158

bar->AddButton("Setup PadRow 20 with Histogram",".p ODH_SetupPadRow(1,2,20)","Setup PadRow");
bar->AddButton("Setup PadRow 20 with Geometry",".p ODH_SetupPadRow(0,2,20)","Setup PadRow");

bar->AddButton("Display PadRow",".p ODH_DisplayNextEvent(false,false,true,NULL)","Display PadRow");

bar->AddButton("Display PadRow with Clusters",".p ODH_DisplayNextEvent(true,false,true,NULL)","Display PadRow");
bar->AddButton("Display PadRow with Tracks",".p ODH_DisplayNextEvent(false,true,true,NULL)","Display PadRow");
bar->AddButton("Display PadRow with Clusters and Tracks",".p ODH_DisplayNextEvent(true,true,true,NULL)","Display PadRow");

bar->Show();
gROOT->SaveContext();
}


