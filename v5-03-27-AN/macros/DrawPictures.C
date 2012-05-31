void DrawPictures()
{
   TControlBar *menu = new TControlBar("vertical","Pictures menu");
   menu->AddButton("TPC shaded",   ".x DrawTPC.C","Draw a shaded view of TPC");
   menu->AddButton("ITS shaded",   ".x DrawITS.C","Draw a shaded view of ITS");
   menu->AddButton("CASTOR shaded",".x DrawCASTOR.C","Draw a shaded view of CASTOR");
   menu->AddButton("ABSO shaded",  ".x DrawABSO.C","Draw a shaded view of ABSO");
   menu->AddButton("DIPO shaded",  ".x DrawDIPO.C","Draw a shaded view of DIPO");
   menu->AddButton("FMD shaded",   ".x DrawFMD.C","Draw a shaded view of FMD");
   menu->AddButton("FRAME shaded", ".x DrawFRAME.C","Draw a shaded view of FRAME");
   menu->AddButton("HALL shaded",  ".x DrawHALL.C","Draw a shaded view of HALL");
   menu->AddButton("MAG shaded",   ".x DrawMAG.C","Draw a shaded view of MAG");
   menu->AddButton("MUON shaded",  ".x DrawMUON.C","Draw a shaded view of MUON");
   menu->AddButton("PHOS shaded",  ".x DrawPHOS.C","Draw a shaded view of PHOS");
   menu->AddButton("PIPE shaded",  ".x DrawPIPE.C","Draw a shaded view of PIPE");
   menu->AddButton("PMD shaded",   ".x DrawPMD.C","Draw a shaded view of PMD");
   menu->AddButton("HMPID shaded",  ".x DrawHMPID.C","Draw a shaded view of HMPID");
   menu->AddButton("SHIL shaded",  ".x DrawSHIL.C","Draw a shaded view of SHIL");
   menu->AddButton("T0 shaded", ".x DrawT0.C","Draw a shaded view of T0");
   menu->AddButton("TOF shaded",   ".x DrawTOF.C","Draw a shaded view of TOF");
   menu->AddButton("TRD shaded",   ".x DrawTRD.C","Draw a shaded view of TRD");
   menu->AddButton("ZDC shaded",   ".x DrawZDC.C","Draw a shaded view of ZDC");
   menu->AddButton("VZERO shaded", ".x DrawVZERO.C","Draw a shaded view of VZERO"); 
   menu->Show();
}
