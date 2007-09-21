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
/* $Id$ */
/** @file    FMDbaseLinkDef.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 14:18:46 2006
    @brief   Link specifications for base library 
*/
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliFMDFlowAxis+;
#pragma link C++ class AliFMDFlowBin+;
#pragma link C++ class AliFMDFlowBin;+;
#pragma link C++ class AliFMDFlowBinned1D+;
#pragma link C++ class AliFMDFlowBin;+;
#pragma link C++ class AliFMDFlowBinned2D+;
#pragma link C++ class AliFMDFlowEventPlane+;
#pragma link C++ class AliFMDFlowHarmonic+;
#pragma link C++ class AliFMDFlowResolution+;
#pragma link C++ class AliFMDFlowResolutionStar+;
#pragma link C++ class AliFMDFlowResolutionTDR+;
#pragma link C++ class AliFMDFlowStat+;
#pragma link C++ class AliFMDFlowTrueBin+;
#pragma link C++ class AliFMDFlowTrue1D+;
#pragma link C++ function  NormalizeAngle(Double_t);
  

#pragma link C++ namespace AliFMDFlowBessel;
#pragma link C++ function  AliFMDFlowBessel::I01(Double_t,Double_t&,Double_t&,Double_t&,Double_t&)
#pragma link C++ function  AliFMDFlowBessel::Ihalf(Int_t,Double_t,Double_t*,Double_t*);
#pragma link C++ function  AliFMDFlowBessel::Iwhole(Int_t,Double_t,Double_t*,Double_t*);
#pragma link C++ function  AliFMDFlowBessel::Inu(Double_t,Double_t,Double_t,Double_t*,Double_t*);
#pragma link C++ function  AliFMDFlowBessel::I(Double_t,Double_t);
#pragma link C++ function  AliFMDFlowBessel::DiffI(Double_t,Double_t);
  

#else
# error Not for compilation 
#endif
//
// EOF
//
