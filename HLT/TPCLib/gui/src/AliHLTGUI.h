
#ifndef _AliHLTGUI_H_
#define _AliHLTGUI_H_

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTGUI.h
    @author Jochen Thaeder
    @date   
    @brief  Qt class for ALICE HLT online Display
*/

#include "Rtypes.h"

#include "AliHLTLogging.h"
#include "AliHLTGUIMainForm.h"
#include <vector>
#include <qtimer.h>



class AliHLTTPCDigitReaderRaw;
// DEL
//#if __GNUC__ >= 3
//using namespace std;
//#endif

/**
 * @class AliHLTGUI
 * The class handels the graphical user interface for the ALICE HLT online 
 * display. The input made in the GUI is set in the main class @ref AliHLTTPCDisplayMain.
 * In order to seperate the Qt code from the ROOT code, all functionality is 
 * implemented in the class @ref AliHLTTPCDisplayMain. This class fuctions just as
 * front-end for the GUI.<br>
 * Several options can alredeay defined using the commandline:
 *   -rawreadermode <mode>              Sets the mode for @ref AliHLTDigitReaderRaw, can bei either digit or string.
 *   -n-time-bins <mode>                Sets number of TimeBins for the TPC, can bei either a digit, 'sim' or 'tpc'.
 *   -adc-threshold <mode>              Sets the ADC Threshold for the Zero-Suppression of raw data, should be ADC counts.
 *   -occupancy-limit <mode>            Sets the Occupancy limits for pads, wherby <mode> is [0,1].
 *   -b-field <mode>                    Sets the B-field.
 *   -tcp-source <host:port> ...        Sets Host and Ports to a TCP port of a @ref TCPDumpSubscriber. Multiple Connections can be given.
 *   -connect                           Connects already to given TCP-Sources.
 *   --help                             Shows the Help Menu.
 *   --version                          Prints version of AliHLTGUI.
 * @ingroup alihlt_display
 */

class AliHLTGUI : public AliHLTGUIMainForm , public AliHLTLogging {
  Q_OBJECT

 public:
  /** standard constructor */
  AliHLTGUI();

  /**
   * Constructor
   * @param argc     Number of commandline arguments
   * @param argv     Array of commandline arguments
   */
  AliHLTGUI( Int_t argc, Char_t **argv );

  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTGUI(const AliHLTGUI&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTGUI& operator=(const AliHLTGUI&);

  /** standard destructor */
  virtual ~AliHLTGUI();
    
  public slots:

  /** Qt slots for the GUI objects<br> */
  /** ---------------------------------*/

  /** Qt slot for connect button */
  void connectDisplay();         
  /** Qt slot for adding host to hostlist */
  void addHost();            
  /** Qt slot for removing host from hostlist */
  void removeHost();
  /** Qt slot for next event button */
  void nextEvent();
  /** Qt slot for event loop button */
  void eventLoop();
  /** Qt slot for update time of evemt loop */
  void updateEventLoop();
  /** Qt slot for redisplay  button */
  void redisplayEvent(){ redisplay(kFALSE); }


  /** Qt slots for the GUI objects and setter for @ref AliHLTTPCDisplayMain */
  /** --------------------------------------------------------------------- */

  /**
   * Qt slot for objects, responsible for the selection<br>
   * of the slices which should be displayed in the 3D view.<br>
   * Selection is than set to @ref ALiHLTTPCDisplayMain.
   */
  void setSector();

  /**
   * Qt slot for zero suppression.<br>
   * Selection is than set to @ref ALiHLTTPCDisplayMain.<br>
   * @ref redisplay() is called afterwards
   */
  void setZeroSuppression(); 

  /**
   * Qt slot for selection of new slice for raw data.<br>
   * Selection is than set to @ref ALiHLTTPCDisplayMain.<br>
   * @ref redisplay(kTRUE) is called afterwards, on order to read
   * the data for the new slice.
   */
  void setRawSlice();

  /**
   * Qt slot for selection of padrow for raw data.<br>
   * Selection is than set to @ref ALiHLTTPCDisplayMain.<br>
   * @ref redisplay() is called afterwards
   */
  void setPadRow();

  /**
   * Qt slot for selection of pad for raw data.<br>
   * Selection is than set to @ref ALiHLTTPCDisplayMain.<br>
   * @ref redisplay() is called afterwards
   */
  void setPad();

  /** 
   * Qt slot for selection of timebins (all timebins, 
   * one timebin, range of timebin)  for raw data.<br>
   * Selection is than set to @ref ALiHLTTPCDisplayMain.<br>
   * @ref redisplay() is called afterwards
   */
  void setTimeBin();

  /**
   * Qt slot for selection, how timebins should be displayed (sum, average, maximum).<br>
   * Selection is than set to @ref ALiHLTTPCDisplayMain.<br>
   * @ref redisplay() is called afterwards
   */
  void selectTimeBinData();     
  
  /**
   * Qt slot for splitting padrowCanvas and show padrow 
   * and pad histograms in the same canvas.<br>
   * Selection is than set to @ref ALiHLTTPCDisplayMain.<br>
   * @ref redisplay() is called afterwards
   */
  void setSplitPadRow();

  /**
   * Qt slot for splitting frontCanvas and show front 
   * and pad histograms in the same canvas.<br>
   * Selection is than set to @ref ALiHLTTPCDisplayMain.<br>
   * @ref redisplay() is called afterwards
   */
  void setSplitFront();

  /**
   * Qt slot for selection of what to display in 3D: raw data, cluster data, tracks data, geometry.<br>
   * Selection is than set to @ref ALiHLTTPCDisplayMain.<br>
   * @ref redisplay() is called afterwards
   */
  void set3D();

  /**
   * Qt slot for selection if only one padrow or all padrow of selected 
   * slice for raw data should be displayed.<br>
   * Selection is than set to @ref ALiHLTTPCDisplayMain.<br>
   * @ref redisplay() is called afterwards
   */
  void set3DRaw();

  /**
   * Qt slot for selection of cuts for tracks.<br>
   * Selection is than set to @ref ALiHLTTPCDisplayMain.<br>
   * @ref redisplay() is called afterwards
   */
  void setTrack(); 

  /**
   * Qt slot for selection of what clusters shold be displayed: 
   * All clusters, used and unused clusters in terms of associated to found tracks.<br>
   * Selection is than set to @ref ALiHLTTPCDisplayMain.<br>
   * @ref redisplay() is called afterwards
   */
  void setCluster();

  /**
   * Qt slot for inversion of the background in 3D canvas.<br>
   * Selection is than set to @ref ALiHLTTPCDisplayMain.<br>
   * @ref redisplay() is called afterwards
   */
  void setInvert();

  /**
   * Qt slot for selection of single tracks or all tracks.
   * Selection is than set to @ref ALiHLTTPCDisplayMain.<br>
   * @ref redisplay() is called afterwards
   */
  void setSelectTrack();

  /**
   * Qt slot for selection of slice for track and residual  objectes.<br>
   * Selection is than set to @ref ALiHLTTPCDisplayMain.<br>
   * @ref redisplay() is called afterwards
   */
  void selectTrackSector();

  /**
   * Qt slot for selection of single track.<br>
   * Selection is than set to @ref ALiHLTTPCDisplayMain.<br>
   * @ref redisplay() is called afterwards
   */
  void selectTrack();

  /**
   * Qt slot for selection, if angles should be kept when redisplaying
   * Selection is than set to @ref ALiHLTTPCDisplayMain.<br>
   * @ref redisplay() is called afterwards
   */
  void setKeepView();
    
  /** Functions<br> */
  /** --------------*/

  /** Shows Track parameter of current selected single Track */
  void displayTrackParam();      

  /**
   * Redisplays event without reading new event, 
   * after change of parameters in the GUI.
   * @param newRawSlice        In the case of changed slice for raw data set
   *                           to kTRUE, than data for new raw slice is read.
   *                           Standard is kFALSE
   */
  void redisplay(Bool_t newRawSlice=kFALSE); // redisplay function

  /** Qt slot for save histogram button, saves all canvases. */
  void saveHistograms();
  
  /** 
   * Enables connect-objects according to connection Status 
   * @param connected          Connection status of the GUI to  
   *                           the TCPDumpSubscriber.<br>
   *                           Either kTRUE or kFALSE.
   */
  void enableDisplay(Bool_t connected);           

  /** Enables data objects according to data availibility  */
  void enableDataObjects();

  /** Sets pad value in GUI out of @ref AliHLTTPCDisplayMain */
  void setGUIPad(Int_t pad); 
  
  /** Callback function fot setGUIPad */
  static void callback_setGUIPad(void* pt2Object, Int_t pad){
    AliHLTGUI* myself = (AliHLTGUI*) pt2Object;
    myself->setGUIPad(pad);
  }

 protected:
  /** Pointer to AliHLTTPCDisplayMain */
  void* fDisplay;
  
  /** Qt vector for hostnames */
  vector<QString> fHostnames;
  /** Qt vector for ports */
  vector<QString> fPorts;
  
  /** QtWidget for ROOT padrow canvas */
  TQtWidget *padrowWidget;
  /** QtWidget for ROOT pad canvas */
  TQtWidget *padWidget;
  /** QtWidget for ROOT residuals canvas */
  TQtWidget *residualsWidget;
  /** QtWidget for ROOT charge canvas */
  TQtWidget *chargeWidget;
  /** QtWidget for ROOT threeD canvas */
  TQtWidget *threeDWidget;
  /** QtWidget for ROOT front canvas */
  TQtWidget *frontWidget;
  /** QtWidget for ROOT hist_s canvas */
  TQtWidget *hits_sWidget;
  /** QtWidget for ROOT q_track canvas */
  TQtWidget *q_trackWidget;
  /** QtWidget for ROOT q_s canvas */
  TQtWidget *q_sWidget;
  /** QtWidget for ROOT padrow_pad canvas */
  TQtWidget *padrow_padWidget;

  /** Qt Timer for event loop */
  QTimer *fEventLoop;
};

#endif /*  _AliHLTGUI_H_ */
