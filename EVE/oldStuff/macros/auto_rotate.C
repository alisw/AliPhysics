// $Header: /soft/cvsroot/AliRoot/EVE/test-macros/tpc_gui.C,v 1.7 2006/10/26 13:24:33 mtadel Exp $

// Function to spawn a gui for reading rootified raw-data from TPC sector test.
//
// To use:
// a) select TPCLoader entry in the list-tree view;
//    you'll get a dialog to steer the data-loading process in an adjacent window
// b) to select a ROOT file containing the raw-data double-click on 'File:'
//    text entry to spawn a file-dialog or type in the name
// c) click open to actually open the file and load an event

#ifndef __CINT__

#include <TEveUtil.h>
#include "TGLViewer.h"
#include "TTimer.h"

#endif

TTimer   *g_rotate_timer = 0;
Double_t  g_rotate_speed = 1;
Double_t  g_rotate_theta = 0;

void auto_rotate(Long_t time=25, Double_t speed=1)
{
  if (g_rotate_timer == 0)
  {
    g_rotate_timer = new TTimer;
  }
  g_rotate_speed = speed;
  g_rotate_theta = 0;
  g_rotate_timer->SetCommand("auto_rotate_camera()");
  g_rotate_timer->Start(time);
}

void auto_rotate_stop()
{
  if (g_rotate_timer == 0)
  {
    Error("auto_rotate_stop", "timer not initialized.");
    return;
  }
  g_rotate_timer->Stop();
}

void auto_rotate_camera()
{
   static Double_t hRotateStep = 0.005;
   static Double_t vRotateStep = 0.025;

   g_rotate_theta += hRotateStep * g_rotate_speed;
   if (g_rotate_theta >= 0.8 || g_rotate_theta <= -0.8) 
   {
     hRotateStep = -hRotateStep;
   }

   TGLViewer *v   = gEve->GetDefaultGLViewer();
   TGLCamera &cam = v->CurrentCamera();
   cam.RotateRad(hRotateStep * g_rotate_speed, vRotateStep * g_rotate_speed);
   v->RequestDraw(TGLRnrCtx::kLODHigh);
}
