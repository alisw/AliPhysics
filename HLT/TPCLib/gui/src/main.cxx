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

#include <qapplication.h>
#include "AliHLTGUI.h"

int main( int argc, char ** argv ) {
  QApplication qApp( argc, argv );
  AliHLTGUI gui( argc,  argv );
  gui.show();
  qApp.connect( &qApp, SIGNAL( lastWindowClosed() ), &qApp, SLOT( quit() ) );
  return qApp.exec();
}
