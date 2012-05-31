TControlBar *menu;

void menu(const char *config="Config.C")
{
   menu = new TControlBar("vertical","gAlice menu");
   menu->AddButton("    Help to run gAlice    ","AliceHelp()", "Explains how to use gAlice menus");
   menu->AddButton("Run",             "gAlice->Run()","Process an Alice event");
   menu->AddButton("Run Lego",        "gAlice->RunLego()","Special run to generate the radl/absl lego plots");
   menu->AddButton("Top view",        "DrawTopView()","Draw Top view (cut) of Alice");
   menu->AddButton("Side view",       "DrawSideView()","Draw Side view (cut) of Alice");
   menu->AddButton("Front view",      "DrawFrontView()","Draw Front view (cut) of Alice");
   menu->AddButton("Menu Trees",      ".x DrawTrees.C","Menu to display the Alice Geant trees");
   menu->AddButton("Menu Pictures",   ".x DrawPictures.C","Menu to display detectors in shaded mode");
   menu->AddButton("Hide ON",         "SetHide(1)","Activate drawing option HIDE");
   menu->AddButton("Hide OFF",        "SetHide(0)","DeActivate drawing option HIDE");
   menu->AddButton("Shading ON",      "SetShade(1)","Activate drawing option SHAD");
   menu->AddButton("Shading OFF",     "SetShade(0)","DeActivate drawing option SHAD");
//   menu->AddButton("RayTracing ON",   "SetRayTracing(1)","Activate drawing option RAYT");
//   menu->AddButton("RayTracing OFF",  "SetRayTracing(0)","DeActivate drawing option RAYT");
   menu->AddButton("Box Clip ON",     "SetBoxClip(1)","Activate clipping box");
   menu->AddButton("Box Clip OFF",    "SetBoxClip(0)","DeActivate clipping box");
//   menu->AddButton("Parallel view",   "SetViewType(1)","Set Parallel view");
//   menu->AddButton("Perspective view","SetViewType(0)","Set Perspective view");
   gROOT->SaveContext();

   gAlice->Init(config);
   ((TGeant3*)gMC)->InitHIGZ();
   menu->Show();
}

void AliceHelp()
{
   gSystem->Exec("nedit AliceHelp.C &");
}

void SetRange()
{
   THIGZ *higz = (THIGZ*)gROOT->GetListOfCanvases()->FindObject("higz");
   if (higz) higz->Range(0,0,20,20);
}

void DrawTopView()
{
   printf("Generating TOP view of Alice. Be patient one minute!\n");
   SetRange();
   ((TGeant3*)gMC)->Gdrawc("ALIC",2,0,10,10,0.01,0.01);
}

void DrawSideView()
{
   printf("Generating SIDE view of Alice. Be patient one minute!\n");
   SetRange();
   ((TGeant3*)gMC)->Gdrawc("ALIC",1,0,10,10,0.01,0.01);
}

void DrawFrontView()
{
   printf("Generating FRONT view of Alice. Be patient one minute!\n");
   SetRange();
   ((TGeant3*)gMC)->Gdrawc("ALIC",3,0,10,10,0.01,0.01);
}

void SetHide(Int_t opt)
{
   if (opt) {
      printf("Drawing option HIDE is ACTIVE!\n");
      ((TGeant3*)gMC)->Gdopt("HIDE","ON");
   } else {
      printf("Drawing option HIDE is NOT ACTIVE!\n");
      ((TGeant3*)gMC)->Gdopt("HIDE","OFF");
   }
}

void SetShade(Int_t opt)
{
   if (opt) {
      printf("Drawing option SHAD is ACTIVE!\n");
      ((TGeant3*)gMC)->Gdopt("SHAD","ON");
   } else {
      printf("Drawing option SHAD is NOT ACTIVE!\n");
      ((TGeant3*)gMC)->Gdopt("SHAD","OFF");
   }
}

void SetRayTracing(Int_t opt)
{
   if (opt) {
      printf("Drawing option RAYT is ACTIVE!\n");
      ((TGeant3*)gMC)->Gdopt("RAYT","ON");
   } else {
      printf("Drawing option RAYT is NOT ACTIVE!\n");
      ((TGeant3*)gMC)->Gdopt("RAYT","OFF");
   }
}

void SetBoxClip(Int_t opt)
{
   if (opt) {
      printf("Clipping BOX is ACTIVE!\n");
      ((TGeant3*)gMC)->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
   } else {
      printf("Clipping BOX is NOT ACTIVE!\n");
      ((TGeant3*)gMC)->SetClipBox(".");
   }
}

void SetViewType(Int_t opt)
{
   if (opt) {
      printf("Setting PARALLEL view!\n");
      ((TGeant3*)gMC)->Gdopt("PROJ","PARA");
   } else {
      printf("Setting PERSPECTIVE view!\n");
      ((TGeant3*)gMC)->Gdopt("PROJ","PERS");
   }
}
