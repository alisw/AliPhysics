/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>                *
 *          for The ALICE Off-line Project.                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTGUI.cxx
    @author Jochen Thaeder
    @date   
    @brief  Qt class for ALICE HLT online Display
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include <Rtypes.h>

#include "AliHLTGUI.h"
#include "AliHLTTPCDisplayMain.h"
#include "AliHLTLogging.h"
#include "AliHLTTPCDigitReaderRaw.h"

#include <TQtWidget.h>

#include <qlabel.h>
#include <qcheckbox.h>
#include <qspinbox.h>
#include <qbuttongroup.h>
#include <qgroupbox.h>
#include <qpushbutton.h>
#include <qlineedit.h>
#include <qlistview.h>
#include <qslider.h>
#include <qtabwidget.h>
#include <qtimer.h>

#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>
#include <cerrno>
#include <iostream>
#include <vector>


// -- standard constructor
//____________________________________________________________________________________________________
AliHLTGUI::AliHLTGUI()
  :
  fDisplay(NULL),
  fHostnames(),
  fPorts(),
  padrowWidget(NULL),
  padWidget(NULL),
  residualsWidget(NULL),
  chargeWidget(NULL),
  threeDWidget(NULL),
  frontWidget(NULL),
  hits_sWidget(NULL),
  q_trackWidget(NULL),
  q_sWidget(NULL),
  padrow_padWidget(NULL),
  fEventLoop(NULL) 
{
}

// -- not a valid copy constructor, defined according to effective C++ style
//____________________________________________________________________________________________________
AliHLTGUI::AliHLTGUI(const AliHLTGUI& srcGUI)
  :
  fDisplay(srcGUI.fDisplay),
  fHostnames(srcGUI.fHostnames),
  fPorts(srcGUI.fPorts),
  padrowWidget(srcGUI.padrowWidget),
  padWidget(srcGUI.padWidget),
  residualsWidget(srcGUI.residualsWidget),
  chargeWidget(srcGUI.chargeWidget),
  threeDWidget(srcGUI.threeDWidget),
  frontWidget(srcGUI.frontWidget),
  hits_sWidget(srcGUI.hits_sWidget),
  q_trackWidget(srcGUI.q_trackWidget),
  q_sWidget(srcGUI.q_sWidget),
  padrow_padWidget(srcGUI.padrow_padWidget),
  fEventLoop(srcGUI.fEventLoop) 
{
  HLTFatal("copy constructor not implemented");
}

// -- not a valid assignment op, but defined according to effective C++ style
//____________________________________________________________________________________________________
AliHLTGUI& AliHLTGUI::operator=(const AliHLTGUI&) {
  HLTFatal("assignment operator not implemented");
  return (*this);
}

// -- constructor
//____________________________________________________________________________________________________
AliHLTGUI::AliHLTGUI( Int_t argc, Char_t **argv )
  :
  fDisplay(NULL),
  fHostnames(),
  fPorts(),
  padrowWidget(NULL),
  padWidget(NULL),
  residualsWidget(NULL),
  chargeWidget(NULL),
  threeDWidget(NULL),
  frontWidget(NULL),
  hits_sWidget(NULL),
  q_trackWidget(NULL),
  q_sWidget(NULL),
  padrow_padWidget(NULL),
  fEventLoop(NULL) 
{
  cout << "Creating display" << endl;

  // Set Widgets 
  padrowWidget = tQtWidgetPadRow; 
  residualsWidget = tQtWidgetResiduals;
  threeDWidget = tQtWidget3D;
  padWidget  = tQtWidgetPad; 
  chargeWidget = tQtWidgetCharge;
  frontWidget = tQtWidgetFront;
  hits_sWidget = tQtWidgetHits_s;
  q_trackWidget = tQtWidgetq_track;
  q_sWidget = tQtWidgetQ_s;
  padrow_padWidget = tQtWidgetPadRow_Pad;
  
  AliHLTTPCDisplayMain* display = new AliHLTTPCDisplayMain((void*) this, &AliHLTGUI::callback_setGUIPad);
  if ( !display ) { 
    HLTFatal( "ERROR no Display!" );
    exit(-1);
  }
   
  fDisplay = (void*)display;


#if STEP1

  display->RegisterCanvas( charge, chargeWidget->GetCanvas() );
  display->RegisterCanvas( padrow, padrowWidget->GetCanvas() );
  display->RegisterCanvas( threeD, threeDWidget->GetCanvas() );
  display->RegisterCanvas( pad, padWidget->GetCanvas() );
  display->RegisterCanvas( residuals, residualWidget->GetCanvas() );
  display->RegisterCanvas( front, frontWidget->GetCanvas() );
  display->RegisterCanvas( hits_s, hits_sWidget->GetCanvas() );
  display->RegisterCanvas( q_track, q_trackWidget->GetCanvas() );
  display->RegisterCanvas( q_s, q_sWidget->GetCanvas() );
  display->RegisterCanvas( padrow_pad, padrow_padWidget->GetCanvas() );

#else
  // Set Canvas  // TODO : REDO with enum
  display->SetCanvasCharge(chargeWidget->GetCanvas());
  display->SetCanvasPadRow(padrowWidget->GetCanvas());
  display->SetCanvas3D(threeDWidget->GetCanvas());
  display->SetCanvasPad(padWidget->GetCanvas());
  display->SetCanvasResiduals(residualsWidget->GetCanvas());
  display->SetCanvasFront(frontWidget->GetCanvas());
  display->SetCanvasHits_S(hits_sWidget->GetCanvas());	
  display->SetCanvasQ_Track(q_trackWidget->GetCanvas());		
  display->SetCanvasQ_S(q_sWidget->GetCanvas());
#if ALIHLTTPCDISPLAY_VERSION >= 2
  display->SetCanvasPadRow_Pad(padrow_padWidget->GetCanvas());
#endif

#endif


  // -- Timer for Event-Loop
  // -----------------------
  fEventLoop = new QTimer( this );
  if ( !fEventLoop ) { 
    HLTFatal( "ERROR no QTimer" );
    exit(-1);
  }

  connect( fEventLoop, SIGNAL(timeout()),this, SLOT(nextEvent()) );

  // -- read commandline arguments and set values in AliHLTTPCDiplayMain
  // -------------------------------------------------------------------
  Int_t ii = 0;
  Char_t* cpErr;                     // error string
  Int_t errocc = 0;                  // check if error occured
  Bool_t connectToSource = kFALSE;   // check if connect to source on startup

#if 0
  while ( ii < argc ) {      

    // -- raw reader mode
    if ( !strcmp( argv[ii], "-rawreadermode" ) ) {
      if ( argc <= ii+1 ) {
	HLTError( "Raw Reader Mode not specified" );
	errocc = ENOTSUP;
	break;
      }
	
      //      Int_t readerMode = AliHLTTPCDigitReaderRaw::DecodeMode( argv[ii+1] );
      //      Int_t readerMode = AliHLTTPCDigitReaderRaw::DecodeMode( 1 );
      Int_t readerMode = 0;
      if (readerMode < 0){
	HLTError ("Cannot convert rawreadermode specifier '%s'", argv[ii+1]);
	errocc = EINVAL;
	break;
      }
      
      //      display->SetRawReaderMode( readerMode );
      
      ii += 2;
      continue;
    }
    
    // -- zero suppression threshold
    if ( !strcmp( argv[ii], "-adc-threshold" ) ) {
      if ( argc <= ii+1 ) {
	HLTError( "adc-threshold not specified" );
	errocc = ENOTSUP;
	break;
      }
      
      Int_t sigThresh = strtoul( argv[ii+1], &cpErr ,0 );
      if ( *cpErr ) {
	HLTError("Cannot convert adc-threshold specifier '%s'.", argv[ii+1]);
	errocc = EINVAL;
	break;
      }

      //      display->SetZeroSuppressionThreshold(sigThresh);
      
      ii+=2;
      continue;
    }

    // -- pad occupancy limit
    if ( !strcmp( argv[ii], "-occupancy-limit" ) ) {
      if ( argc <= ii+1 ) {
	HLTError( "occupancy-limit not specified" );
	errocc = ENOTSUP;
	break;
      }
      
      Float_t occuLimit = strtof( argv[ii+1], &cpErr );
      if ( *cpErr ) {
	HLTError("Cannot convert occupancy-limit  specifier '%s'.", argv[ii+1]);
	errocc = EINVAL;
	break;
      }
      
      //      display->SetOccupancyLimit(occuLimit);
      
      ii+=2;
      continue;
    }
    
    // -- B field
    if ( !strcmp( argv[ii], "-b-field" ) ) {
      if ( argc <= ii+1 ) {
	HLTError( "b-field not specified" );
	errocc = ENOTSUP;
	break;
      }
      
      Float_t bfield = strtof( argv[ii+1], &cpErr );
      if ( *cpErr ) {
	HLTError("Cannot convert b-field specifier '%s'.", argv[ii+1]);
	errocc = EINVAL;
	break;
      }
      
      //      display->SetBField( bfield )
	
	ii+=2;
      continue;
    }
    
    // -- N Time Bins
    if ( !strcmp( argv[ii], "-n-time-bins" ) ) {
      if ( argc <= ii+1 ) {
	HLTError( "n-time-bins not specified" );
	errocc = ENOTSUP;
	break;
      }

      Int_t nTimeBins = strtoul( argv[ii+1], &cpErr ,0);
      // check if already a number, if not check if "tpc" or "sim"
      if ( *cpErr ) {
	if ( !strcmp( argv[ii+1], "sim" ) ) 
	  nTimeBins = 442;	  
	else if  ( !strcmp( argv[ii+1], "tpc") ) 
	  nTimeBins = 1024;
	else {
	  HLTError("Cannot convert n-time-bins specifier '%s'.", argv[ii+1]);
	  errocc = EINVAL;
	  break;
	}
      }
      
      //      display->SetNTimeBins( nTimeBins );
      
      ii+=2;
      continue;
    }

    // -- TCP soucre
    if ( !strcmp( argv[ii], "-tcp-source" ) ) {
      if ( argc <= ii+1 ) {
	HLTError( "-tcp-source not specified" );
	errocc = ENOTSUP;
	break;
      }
      
      Int_t sourceCount = 1;

      while ( 1 ){
	// break if more argc is exceeded
	if ( argc <= ii+sourceCount ) break;

	// check is next is new option and source count != 1 ( because at least 1 host:port is needed)
	if ( strstr( argv[ii + sourceCount], "-" ) ) {
	  if ( sourceCount == 1 ){
	    HLTError( "-tcp-source not specified" );
	    errocc = ENOTSUP;
	  }
	  break;
	}

	// break if argument is no host:port
	if ( !strstr( argv[ii + sourceCount], ":" ) ) break;

	Char_t *host = strtok( argv[ii + sourceCount], ":" );
	Char_t *port = strtok( NULL, ":" );

	if ( !port ) {	  
	  HLTError("Cannot convert tcp-source specifier '%s'.", argv[ii+1]);
	  errocc = EINVAL;
	  break;
	}

	// push in vector
	fHostnames.push_back( host );
	fPorts.push_back( port );
    
	// write in List
	QListViewItem * item1 = new QListViewItem( hostListView, host, port );
	hostListView->insertItem(item1);

	sourceCount++;
      } // while ( 1 ) {

      if  ( errocc != 0 ) break;
      
      ii += sourceCount;
      continue;
    }

    // -- Connect
    if ( !strcmp( argv[ii], "--connect" ) ) {
      
      connectToSource = kTRUE;

      ii++;
      continue;
    }

    // -- Help
    if ( !strcmp( argv[ii], "--help" ) ) {
      cout << "Usage: AliHLTGUI [OPTON VALUE] ..." << endl;
      cout << " OPTIONs:\n --------" << endl;
      cout << " -rawreadermode <mode>              Sets the mode for AliHLTDigitReaderRaw, can bei either digit or string." << endl;
      cout << " -n-time-bins <mode>                Sets number of TimeBins for the TPC, can bei either a digit, 'sim' or 'tpc'." << endl;
      cout << " -adc-threshold <mode>              Sets the ADC Threshold for the Zero-Suppression of raw data, should be ADC counts." << endl;
      cout << " -occupancy-limit <mode>            Sets the Occupancy limits for pads, wherby <mode> is [0,1]." << endl;
      cout << " -b-field <mode>                    Sets the B-field." << endl;
      cout << " -tcp-source <host:port> ...        Sets Host and Ports to a TCP port of a TCPDumpSubscriber. Multiple Connections can be given." << endl;
      cout << " -connect                           Connects already to given TCP-Sources." << endl;
      cout << " --help                             Shows the Help Menu." << endl;
      cout << " --version                          Prints version of AliHLTGUI" << endl;
      
      exit (0);
    }
    
    // -- Version
    if ( !strcmp( argv[ii], "--version" ) ) {
      cout << "Version of AliHLTGUI: " << PACKAGE_VERSION << endl;
      
      exit (0);
    }

    // -- Unknown Option
    HLTError( "Unknown Option", "Unknown option '%s'", argv[ii] );
    errocc = EINVAL;
    break;
  } // end while ( ii < argc ) {      

#endif

  // if error occured: Print usage and exit
  if ( errocc != 0 ) {
    cout << "Usage: AliHLTGUI [OPTON VALUE] ..." << endl;
    cout << " OPTIONs:\n --------" << endl;
    cout << " -rawreadermode <mode>              Sets the mode for AliHLTDigitReaderRaw, can bei either digit or string." << endl;
    cout << " -n-time-bins <mode>                Sets number of TimeBins for the TPC, can bei either a digit, 'sim' or 'tpc'." << endl;
    cout << " -adc-threshold <mode>              Sets the ADC Threshold for the Zero-Suppression of raw data, should be ADC counts." << endl;
    cout << " -occupancy-limit <mode>            Sets the Occupancy limits for pads, wherby <mode> is [0,1]." << endl;
    cout << " -b-field <mode>                    Sets the B-field." << endl;
    cout << " -tcp-source <host:port> ...        Sets Host and Ports to a TCP port of a TCPDumpSubscriber. Multiple Connections can be given." << endl;
    cout << " -connect                           Connects already to given TCP-Sources." << endl;
    cout << " --help                             Shows the Help Menu." << endl;
    cout << " --version                          Prints version of AliHLTGUI" << endl;
    exit( errocc );
  }


  if ( connectToSource ) connectDisplay();

}

//____________________________________________________________________________________________________
AliHLTGUI::~AliHLTGUI() {     
  // -- destructor

  if (fEventLoop)
    delete fEventLoop;
  fEventLoop = NULL;

  AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*) fDisplay;

  if (display)
    delete display;
  fDisplay = NULL;
}

//____________________________________________________________________________________________________
void AliHLTGUI::connectDisplay() {  
    // -- connect to Display class / TCPDumpSubscriber

  AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;
  

#ifdef STEP2
  enableDisplay( display->GetConnectionStatus() );
#else
  // if connected => diconnect
  if ( display->GetConnectionStatus() ) {
    cout << "Disonnecting..."    << endl;
    display->SetConnectionStatus(kFALSE);	
    // -- CONNECT-TAB
    connectPushButton->setText("Connect");
    connectTextLabel->setText("DISCONNECTED");
    connectTextLabel->setPaletteForegroundColor(red);	
    HostListTextLabel->setEnabled(TRUE);
    hostListView->setEnabled(TRUE);
    removePushButton->setEnabled(TRUE);
    addGroupBox->setEnabled(TRUE);
    // -- DISPLAY-TAB
    eventGroupBox->setEnabled(FALSE);
    
    
#if ALIHLTTPCDISPLAY_VERSION >= 2
    display->Disconnect();
#else
    HLTWarning("disconnect not available for AliHLTTPCDisplay version < 2");
#endif
    
    return;
  }
#endif

  cout << "Connecting..."    << endl;


  // -- Set geometry file in AliHLTTPCDisplayMain
  // --------------------------------------------
  Char_t * geometryPath = getenv("ALIHLT_TOPDIR");
  Char_t * geometryFile = "alice.geom";
  
  Char_t geoPath[256];

  strcpy(geoPath,"");

  if (geometryPath) {
    strcat(geoPath,geometryPath);
    strcat(geoPath,"/TPCLib/OnlineDisplay");
  }
  else strcat(geoPath,".");
  
  strcat(geoPath,"/");
  strcat(geoPath,geometryFile);

  // -- Set Hostnames in AliHLTTPCDisplayMain
  // ----------------------------------------
  Int_t cnt = fHostnames.size();
  
  if (cnt == 0){
    HLTError( "Error establishing connection: No Hostnames specified" );
    connectTextLabel->setText("ERROR");  
    connectTextLabel->setPaletteForegroundColor(red);
    return;
  }

  unsigned short* ports = new unsigned short[cnt];
  const char** hostnames = new const char*[cnt];

  for ( Int_t i = 0; i < cnt; i++ ) {
    const char * temp = fHostnames[i].latin1();
    hostnames[i] = temp;
    ports[i] = (unsigned short) atoi(fPorts[i]);
  }

  Int_t connect = display->Connect( cnt, (const char**) hostnames, ports, geoPath);
  
  delete[] ports;
  delete[] hostnames; 
  
  // check return value
  if ( connect  == 0) {
    connectTextLabel->setText("CONNECTED");
    connectTextLabel->setPaletteForegroundColor(green);
  }
  else if (connect == -1){
    connectTextLabel->setText("NO CHAIN");
    connectTextLabel->setPaletteForegroundColor(red);
    return;
    }
  else if (connect == -2){
    connectTextLabel->setText("TIMEOUT");
    connectTextLabel->setPaletteForegroundColor(red);
    return;
  }
  else {
    connectTextLabel->setText("GENERAL ERROR");
    connectTextLabel->setPaletteForegroundColor(red);
    return;
  }
#ifdef STEP3

  enableDataObjects();
#else
  // if raw data present
    if ( display->ExistsRawData() ) {
      rawZeroGroupBox->setEnabled(TRUE);
      rawSliceGroupBox->setEnabled(TRUE);
      padrowGroupBox->setEnabled(TRUE);
      frontDataButtonGroup->setEnabled(TRUE);
      frontDisplayButtonGroup->setEnabled(TRUE);
      padrowCheckBox->setEnabled(TRUE);
      threeDRawButtonGroup->setEnabled(TRUE);
    }
    else {
      rawZeroGroupBox->setEnabled(FALSE);
      rawSliceGroupBox->setEnabled(FALSE);
      padrowGroupBox->setEnabled(FALSE);
      frontDataButtonGroup->setEnabled(FALSE);
      frontDisplayButtonGroup->setEnabled(FALSE);
      padrowCheckBox->setEnabled(FALSE);
      threeDRawButtonGroup->setEnabled(FALSE);
    }


    // if cluster data present
    if ( display->ExistsClusterData() ) {
	clusterCheckBox->setEnabled(TRUE);
	clusterPropButtonGroup->setEnabled(TRUE);
    }
    else {
	clusterCheckBox->setEnabled(FALSE);
	clusterPropButtonGroup->setEnabled(FALSE);
    }

    // if track data present
    if ( display->ExistsTrackData() ) {
	tracksCheckBox->setEnabled(TRUE);
	trackPropGroupBox->setEnabled(TRUE);
	selectTrackGroupBox->setEnabled(TRUE);
	trackParamGroupBox->setEnabled(TRUE);
    }
    else {
	tracksCheckBox->setEnabled(FALSE);
	trackPropGroupBox->setEnabled(FALSE);
	selectTrackGroupBox->setEnabled(FALSE);
	trackParamGroupBox->setEnabled(FALSE);
    }
#endif

#ifdef STEP2
    enableDisplay( display->GetConnectionStatus() );
#else
    // -- CONNECT-TAB
    connectPushButton->setText("Disconnect");
    connectTextLabel->setText("READ DATA");	
    connectTextLabel->setPaletteForegroundColor(green);
    HostListTextLabel->setEnabled(FALSE);
    hostListView->setEnabled(FALSE);
    removePushButton->setEnabled(FALSE);
    addGroupBox->setEnabled(FALSE);
	
    // -- DISPLAY-TAB
    eventGroupBox->setEnabled(TRUE);
    nextEventPushButton->setEnabled(TRUE);
    redisplayPushButton->setEnabled(TRUE);
    saveGroupBox->setEnabled(TRUE);
    sectorButtonGroup->setEnabled(TRUE);
    threeDGroupBox->setEnabled(TRUE);
#endif

       
    display->ReadData(kFALSE);
    display->DisplayEvent();

    Char_t eventID[256] = "";
#if ALIHLTTPCDISPLAY_VERSION >= 2
    sprintf(eventID,"%Lu",display->GetEventID());
#endif
    sprintf(eventID,"%Lu",0);
    eventIDTextLabel->setText(eventID);

    displayTrackParam();
}

//____________________________________________________________________________________________________
void AliHLTGUI::nextEvent(){
  // -- Display next event
  AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;
  Int_t ret = 0;
#if ALIHLTTPCDISPLAY_VERSION >= 2
  ret = display->ReadData();
#else
  display->ReadData();
#endif
  if (ret == -1){
    connectTextLabel->setText("NO CHAIN");
    setPaletteForegroundColor(red);
    
    // Stop of Event Loop
    if (!fEventLoop->isActive()){
      fEventLoop->stop();
      eventLoopPushButton->setText("Start Loop");
    }
    
    return;
  }
  else if (ret == -2) {
    connectTextLabel->setText("TIME OUT");    
    connectTextLabel->setPaletteForegroundColor(red);

    // Stop of Event Loop
    if (!fEventLoop->isActive()){
      fEventLoop->stop();
      eventLoopPushButton->setText("Start Loop");
    }

    return;
  }
  else if (ret == -3) {
    connectTextLabel->setText("GENERAL ERROR");    
    connectTextLabel->setPaletteForegroundColor(red);

    // Stop of Event Loop
    if (!fEventLoop->isActive()){
      fEventLoop->stop();
      eventLoopPushButton->setText("Start Loop");
    }

    return;
  }
  else {
    connectTextLabel->setText("READ DATA");    
    connectTextLabel->setPaletteForegroundColor(green); 
  }
    
  display->DisplayEvent();

  Char_t eventID[256] = "";
#if ALIHLTTPCDISPLAY_VERSION >= 2
  sprintf(eventID,"%Lu",display->GetEventID());
#endif
  sprintf(eventID,"%Lu",0);
  eventIDTextLabel->setText(eventID);

  displayTrackParam();
}

//____________________________________________________________________________________________________
void AliHLTGUI::eventLoop(){
    // -- Start / Stop of Event Loop
    if (!fEventLoop->isActive()){
	fEventLoop->start(1000 * atoi( eventLoopSpinBox->text() ) );
	eventLoopPushButton->setText("Stop Loop");
    }
    else {
	fEventLoop->stop();
	eventLoopPushButton->setText("Start Loop");
    }    
}


//____________________________________________________________________________________________________
void AliHLTGUI::enableDisplay(Bool_t connected){
  AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*) fDisplay;
  
  Bool_t inverseconnected  = ~connected; // TODO check!!
  
  // if already connected => disconnect
  // else disconnect => connect
  
  HostListTextLabel->setEnabled( connected );
  hostListView->setEnabled( connected );
  removePushButton->setEnabled( connected );
  addGroupBox->setEnabled( connected );
  eventGroupBox->setEnabled( inverseconnected );

  if ( connected ){
    cout << "Disonnecting..."    << endl;
    display->SetConnectionStatus( inverseconnected );	
    // -- CONNECT-TAB
    connectPushButton->setText("Connect");
    connectTextLabel->setText("DISCONNECTED");
    connectTextLabel->setPaletteForegroundColor(red);	
    
#if ALIHLTTPCDISPLAY_VERSION >= 2
    display->Disconnect();
#else
    HLTWarning("disconnect not available for AliHLTTPCDisplay version < 2");
#endif

    return;
  }

  else { // connected = false
    // -- CONNECT-TAB
    connectPushButton->setText("Disconnect");
    connectTextLabel->setText("READ DATA");	
    connectTextLabel->setPaletteForegroundColor(green);
	
    // -- DISPLAY-TAB
    nextEventPushButton->setEnabled( inverseconnected );
    redisplayPushButton->setEnabled( inverseconnected );
    saveGroupBox->setEnabled( inverseconnected );
    sectorButtonGroup->setEnabled( inverseconnected );
    threeDGroupBox->setEnabled( inverseconnected );
  }
  
}

//____________________________________________________________________________________________________
void AliHLTGUI::enableDataObjects(){
  AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*) fDisplay;

  Bool_t existsRawData = display->ExistsRawData();
  Bool_t existsClusterData = display->ExistsClusterData();
  Bool_t existsTrackData = display->ExistsTrackData();

  // -- RAW data
  rawZeroGroupBox->setEnabled(existsRawData);
  rawSliceGroupBox->setEnabled(existsRawData);
  padrowGroupBox->setEnabled(existsRawData);
  frontDataButtonGroup->setEnabled(existsRawData);
  frontDisplayButtonGroup->setEnabled(existsRawData);
  padrowCheckBox->setEnabled(existsRawData);
  threeDRawButtonGroup->setEnabled(existsRawData);

  // -- CLUSTER data
  clusterCheckBox->setEnabled(existsClusterData);
  clusterPropButtonGroup->setEnabled(existsClusterData);

  // -- TRACK data
  tracksCheckBox->setEnabled(existsTrackData);
  trackPropGroupBox->setEnabled(existsTrackData);
  selectTrackGroupBox->setEnabled(existsTrackData);
  trackParamGroupBox->setEnabled(existsTrackData);
}

//____________________________________________________________________________________________________
void AliHLTGUI::updateEventLoop(){
    if (fEventLoop->isActive())
	fEventLoop->changeInterval(1000 * atoi( eventLoopSpinBox->text() ) );
}

//____________________________________________________________________________________________________
void AliHLTGUI::redisplay(Bool_t newRawSlice){
    // -- ReDisplay current event
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;
#if ALIHLTTPCDISPLAY_VERSION >= 2
    display->DisplayEvent(newRawSlice);
#else
    display->DisplayEvent();
#endif

    displayTrackParam();
}

//____________________________________________________________________________________________________
void  AliHLTGUI::saveHistograms(){
    // -- save all Histograms
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;
    display->SaveHistograms();
}

//____________________________________________________________________________________________________
void AliHLTGUI::addHost(){    
    // -- add host to hostlist
    // push in vector
    fHostnames.push_back( hostnameLineEdit->text() );
    fPorts.push_back( portSpinBox->text() );
    
    // write in List
    QListViewItem * item1 = new QListViewItem(hostListView,  hostnameLineEdit->text(), portSpinBox->text() );
    hostListView->insertItem(item1);
}

//____________________________________________________________________________________________________
void AliHLTGUI::removeHost(){    
    // -- remove host from hostlist
    if (hostListView->selectedItem()) {
	for ( unsigned i = 0; i < fHostnames.size(); i++ ) {
	    if ( fHostnames[i] == hostListView->selectedItem()->text(0)  && fPorts[i] == hostListView->selectedItem()->text(1)){
		fHostnames.erase( fHostnames.begin()+i );
		fPorts.erase( fPorts.begin()+i );
	    }
	}	
	delete hostListView->selectedItem();   // #########
	//	hostListView->removeItem(hostListView->selectedItem());
    }
}

//----------------------------------------------------------------------------------------------------
//       SETTER - minor functions
//____________________________________________________________________________________________________
void AliHLTGUI::setSector(){  
    // -- set sectors to display
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;   
      
    Int_t selectSectorID = sectorButtonGroup->selectedId();

    // Set ALL Sectors
    if (selectSectorID == 0) {
	display->SetSlices();
	display->SetTheta(90.);
    }
    // Set ONE Sector
    else if (selectSectorID == 4) {
	display->SetSlices( atoi( oneSpinBox->text() ) );
	display->SetTheta(0.);
    }
    // Set RANGE of Sectors
    else if (selectSectorID == 1) {
	display->SetSlices( atoi( rangeMinSpinBox->text() ) ,atoi( rangeMaxSpinBox->text() ) );
	display->SetTheta(90.);
    }
    // Set PAIR of Sectors
    else if (selectSectorID == 3) { 
	display->SetSlicesPair(atoi( pairSpinBox->text() ) );	
	display->SetTheta(90.);
    }
    // Set PAIR RANGE of Sectors
    else if (selectSectorID == 2) {
	display->SetSlicesPair(atoi( pairRangeMinSpinBox->text() ) ,atoi( pairRangeMaxSpinBox->text() ) );
	display->SetTheta(90.);
    }
    // redisplay 3D
    redisplay();
}

//____________________________________________________________________________________________________
void AliHLTGUI::setCluster(){    
    // -- set used/unuse/all cluster
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;   

    //  ALL Cluster = 0  | USED Cluster = 1 | UNUSED Cluster = 2
    display->SetSelectCluster( clusterPropButtonGroup->selectedId()  );
          
    // redisplay 3D
    redisplay();
}
//____________________________________________________________________________________________________
void AliHLTGUI::setTrack(){    
    // -- set track properties
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;  

    display->SetCutHits(atoi (cutHitsSpinBox->text() ) );
    display->SetCutPt( (Float_t) atof (cutPtLineEdit->text() ) );
    display->SetCutPsi( (Float_t) atof (cutPsiLineEdit->text() ) );
    display->SetCutLambda( (Float_t) atof (cutLambdaLineEdit->text() ) );	
    display->SetCutS( (Float_t) atof (cutSLineEdit->text() ) );
    display->SetIncidentPadrow( (Int_t) atoi (cutPadrowLineEdit->text() ) );
    redisplay();
}

//____________________________________________________________________________________________________
void AliHLTGUI::setInvert(){    
    // -- set invert 3D display
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;  

    display->SetInvert();
    redisplay();
}

//____________________________________________________________________________________________________
void AliHLTGUI::setKeepView(){    
    // -- keep angle view 
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;   

    if (keepViewCheckBox->isChecked() ) display->SetKeepView(kTRUE);
    else display->SetKeepView(kFALSE);
}


void AliHLTGUI::setZeroSuppression(){    
    // -- set Zero Suppresion
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;   

#if ALIHLTTPCDISPLAY_VERSION >= 2
    if (rawZeroCheckBox->isChecked() ) display->SetZeroSuppression(kTRUE);
    else display->SetZeroSuppression(kFALSE);
#else
    HLTError("function not available for AliHLTTPCDisplay version < 2");
#endif

    redisplay();
}


//____________________________________________________________________________________________________
void AliHLTGUI::setRawSlice(){    
    // -- set  slice for padrow / pad / timebin
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;   

    display->SetSlicePadRow( atoi (rawSectorSpinBox->text() ) );
    
    redisplay(kTRUE/*newRawSlice*/);
}

//____________________________________________________________________________________________________
void AliHLTGUI::setPadRow(){    
    // -- set padrow
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;   

    display->SetPadRow(atoi (padrowSpinBox->text() ) );

    padSpinBox->setMaxValue(display->GetNPads());
    padSlider->setMaxValue(display->GetNPads());

    setPad(); 
}

//____________________________________________________________________________________________________
void AliHLTGUI::setPad(){    
    // -- set pad, or pad histograms 
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;  
 
    display->SetPad(atoi (padSpinBox->text() ) );

    redisplay();
}

//____________________________________________________________________________________________________
void AliHLTGUI::setTimeBin(){    
    // -- set timebins
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;  

    Int_t selectTimebinID = frontDisplayButtonGroup->selectedId();

#if ALIHLTTPCDISPLAY_VERSION >= 2
    // Set ALL TimeBins
    if ( selectTimebinID == 3) 	display->SetTimeBinMinMax(0,display->GetNTimeBins()-1 );

    // Set ONE TimeBin
    else if ( selectTimebinID == 2) display->SetTimeBinMinMax( atoi( oneTimeBinSpinBox->text() ), atoi( oneTimeBinSpinBox->text() ) );

    // Set RANGE of TimeBins
    else if ( selectTimebinID == 1) {
      Int_t  min = atoi( minTimeBinSpinBox->text() );
      Int_t  max = atoi( maxTimeBinSpinBox->text() );

      if ( max < min) {
	max = min;
	maxTimeBinSpinBox->setValue(max);
      }

      display->SetTimeBinMinMax( min , max );

    }
#else
    HLTError("function not available for AliHLTTPCDisplay version < 2");
#endif

    redisplay();
}

//____________________________________________________________________________________________________
void AliHLTGUI::selectTimeBinData(){    
    // -- select max,sum,average of Data
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;  

    Int_t selectFrontDataID = frontDataButtonGroup->selectedId();

#if ALIHLTTPCDISPLAY_VERSION >= 2
    // Set use sum in data
    if (selectFrontDataID == 0) display->SetFrontDataSwitch(0);

    // Set use average in data
    else if (selectFrontDataID == 1) display->SetFrontDataSwitch(1);

    // Set use maximum in data
    else if (selectFrontDataID == 2) display->SetFrontDataSwitch(2);
#else
    HLTError("function not available for AliHLTTPCDisplay version < 2");
#endif

    redisplay();
}

//____________________________________________________________________________________________________
void AliHLTGUI::setSplitPadRow(){
    // -- set split PadRow Canvas
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;  

    if (padrowSplitCheckBox->isChecked() ) display->SetSplitPadRow(kTRUE);
    else display->SetSplitPadRow(kFALSE);

    redisplay();
}
//____________________________________________________________________________________________________
void AliHLTGUI::setSplitFront(){
    // -- set split Front Canvas
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;  

#if ALIHLTTPCDISPLAY_VERSION >= 2
    if (frontSplitCheckBox->isChecked() ) display->SetSplitFront(kTRUE);
    else display->SetSplitFront(kFALSE);
#else
    HLTError("function not available for AliHLTTPCDisplay version < 2");
#endif

    redisplay();
}

/* ----- OLD DELETE
//____________________________________________________________________________________________________
void AliHLTGUI::setAllTimeBins(){
    // -- set show all Timebins in front view
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;  

    if (allTimebinCheckBox->isChecked() ) display->SetAllTimebins(kTRUE);
    else display->SetAllTimebins(kFALSE);

    redisplay();
}
----- */

//____________________________________________________________________________________________________
void AliHLTGUI::setGUIPad(Int_t pad){    
    // -- set GUI pad value
    padSpinBox->setValue(pad);
}

//____________________________________________________________________________________________________
void AliHLTGUI::set3D(){    
    // -- set 3D Display
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;  
   
    if  (padrowCheckBox->isChecked()) setPadRow();
    
    display->SetSwitches(tracksCheckBox->isChecked(), clusterCheckBox->isChecked(), padrowCheckBox->isChecked(), geometryCheckBox->isChecked() );
    redisplay();
}

//____________________________________________________________________________________________________
void AliHLTGUI::set3DRaw(){    
    // -- set 3D Display Raw Properties
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;  
   
    Int_t selectRawID = threeDRawButtonGroup->selectedId();

#if ALIHLTTPCDISPLAY_VERSION >= 2
    // Set ALL PadRows
    if (selectRawID == 0) display->Set3DRawSwitch(0);

    // Set ONE PadRow
    else if (selectRawID == 1) display->Set3DRawSwitch(1);
#else
    HLTError("function not available for AliHLTTPCDisplay version < 2");
#endif

    // redisplay 3D
    redisplay();
}

//____________________________________________________________________________________________________
void AliHLTGUI::setSelectTrack(){    
    // --  turn on /off of single track
   AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;  

    display->SetSelectTrackSwitch( selectTrackCheckBox->isChecked() );
  
    if (! selectTrackCheckBox->isChecked() ) {
	display->SetSelectTrack(0);
	redisplay();
    }
    else selectTrackSector();
}

//____________________________________________________________________________________________________
void AliHLTGUI::selectTrackSector(){      
    // -- select slice for single track
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;  

    Int_t slice = atoi( selectTrackSectorSpinBox->text());

    if (display->GetTracksPerSlice(slice) == 0)  selectTrackSpinBox->setMaxValue( 0 ); 
    else selectTrackSpinBox->setMaxValue( display->GetTracksPerSlice(slice)  - 1 );

    selectTrackSpinBox->setValue(0);

    display->SetSelectTrackSlice(slice); 
    
    selectTrack();
}

//____________________________________________________________________________________________________
void AliHLTGUI::selectTrack(){    
    // -- select single track in slice
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;  

    display->SetSelectTrack( atoi( selectTrackSpinBox->text() ) );
    redisplay();
}


//____________________________________________________________________________________________________
void AliHLTGUI::displayTrackParam(){    
    // -- display Track Parameters
    AliHLTTPCDisplayMain* display = (AliHLTTPCDisplayMain*)fDisplay;  
   
    Char_t nHits[256] = "";
    Char_t charge[256] = "";
    Char_t kappa[256] = "";
    Char_t radius[256] = "";
    Char_t slice[256] = "";
    Char_t phi0[256] = "";
    Char_t psi[256] = "";
    Char_t lambda[256]= "";
    Char_t pt[256] = "";
    Char_t id[256] = "";

    if ( display->GetSelectTrackSwitch() 
	 && display->ExistsTrackData() 
	 && display->GetDisplaySlice( atoi(selectTrackSectorSpinBox->text()) )  
	 && display->GetTracksPerSlice( atoi(selectTrackSectorSpinBox->text()) ) > 0 ){ 
	sprintf(kappa,"%f",display->fTrackParam.kappa);
	sprintf(nHits,"%d",display->fTrackParam.nHits);
	sprintf(charge,"%d",display->fTrackParam.charge);
	sprintf(radius,"%f",display->fTrackParam.radius);
	sprintf(slice,"%d",display->fTrackParam.slice);
	sprintf(phi0,"%f",display->fTrackParam.phi0);
	sprintf(psi,"%f",display->fTrackParam.psi);
	sprintf(lambda,"%f",display->fTrackParam.lambda);
	sprintf(pt,"%f",display->fTrackParam.pt);
	sprintf(id,"%d",display->fTrackParam.id);
     } 
     else {
	sprintf(kappa,"%f",0.);
	sprintf(nHits,"%d",1);
	sprintf(charge,"%d",0);
	sprintf(radius,"%f",0.);
	sprintf(slice,"%d",0);
	sprintf(phi0,"%f",0.);
	sprintf(psi,"%f",0.);
	sprintf(lambda,"%f",0.);
	sprintf(pt,"%f",0.);
	sprintf(id,"%d",0);
     }

     trackParamSliceTextLabel->setText(slice);
     trackParamIDTextLabel->setText(id);
     trackParamKappaTextLabel->setText(kappa);
     trackParamPtTextLabel->setText(pt);
     trackParamNHitsTextLabel->setText(nHits);
     trackParamChargeTextLabel->setText(charge);
     trackParamRadiusTextLabel->setText(radius);
     trackParamPhi0TextLabel->setText(phi0);
     trackParamPsiTextLabel->setText(psi);
     trackParamLambdaTextLabel->setText(lambda);
}
