

#ifndef AliTPCCalibViewerGUIAlarms_H
#define AliTPCCalibViewerGUIAlarms_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTPCCalibViewerGUIAlarms.h,v */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  GUI for displaying Alarms of type AliTPCCalibQAChecker                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TGFrame.h>

class TGListTree;
class TGListTreeItem;
class TCanvas;
class TGCanvas;
class TGLabel;

class AliTPCCalibQAChecker;
class AliTPCCalibViewerGUI;
class AliTPCCalibViewerGUItime;

class AliTPCCalibViewerGUIAlarms : public TGCompositeFrame {
public:
  AliTPCCalibViewerGUIAlarms(const TGWindow *p, UInt_t w, UInt_t h);

  virtual ~AliTPCCalibViewerGUIAlarms();

  void SetCalibChecker(AliTPCCalibQAChecker *checker) {fCalibChecker=checker;}
  void SetCalibViewerGUI(AliTPCCalibViewerGUI *gui) {fCalibViewerGUI=gui;}
  void SetCalibViewerGUItime(AliTPCCalibViewerGUItime *gui) {fCalibViewerGUItime=gui;}
  
  void InitBrowser();
  void UpdateBrowser();
  void ResetBrowser();
  void OpenAllItems();

  static AliTPCCalibViewerGUIAlarms* Show();

  void OnDoubleClick(TGListTreeItem* item, Int_t id);
  void OnClick(TGListTreeItem* item, Int_t id);
    
protected:
  AliTPCCalibQAChecker *fCalibChecker;           //Calibration checker
  TGListTree           *fAlarmTree;              //tree representation of alarms
  TCanvas              *fMainCanvas;             //canvas for alarm histogram displaying
  TGCanvas             *fTreeCanvas;             //tree canvas
  TGLabel              *fAlarmText;              //alarm information
  //
  AliTPCCalibViewerGUI *fCalibViewerGUI;         //! pointer to gui
  AliTPCCalibViewerGUItime *fCalibViewerGUItime; //! pointer to gui time
  //
  void DrawGUI(const TGWindow */*p*/, UInt_t w, UInt_t h);
  void AddSubItems(AliTPCCalibQAChecker *fChecker, TGListTreeItem *item);
  void UpdateSubItem(TGListTreeItem *item);
  void OpenSubItems(TGListTreeItem *item);
  
private:
  AliTPCCalibViewerGUIAlarms(const AliTPCCalibViewerGUIAlarms &v);
  AliTPCCalibViewerGUIAlarms &operator = (const AliTPCCalibViewerGUIAlarms &v);         // assignment operator
  

  ClassDef(AliTPCCalibViewerGUIAlarms,0);
};

#endif


