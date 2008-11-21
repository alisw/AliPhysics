// -*- mode: c++ -*- 
/* Copyright (C) 2007 Christian Holm Christensen <cholm@nbi.dk>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 */
#ifdef __CINT__
/* $Id: FMDflowLinkDef.h 23165 2007-12-19 01:36:20Z cholm $ */
/** @file    FMDbaseLinkDef.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 14:18:46 2006
    @brief   Link specifications for base library 
*/
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;


#pragma link C++ class AliFMDAnaCalibBackgroundCorrection+;
#pragma link C++ class AliFMDAnaCalibEnergyDistribution+;
#pragma link C++ class AliFMDAnaParameters+;
#pragma link C++ class AliFMDAnaCalibEnergyDistribution+;
#pragma link C++ class AliFMDAnaCalibBackgroundCorrection+;
#pragma link C++ class AliFMDAnalysisTaskESDReader+;
#pragma link C++ class AliFMDAnalysisTaskSharing+;
#pragma link C++ class AliFMDAnalysisTaskDensity+;
#pragma link C++ class AliFMDAnalysisTaskBackgroundCorrection+;
#pragma link C++ class AliFMDAnalysisTaskCollector+;
#else
# error Not for compilation 
#endif
//
// EOF
//
