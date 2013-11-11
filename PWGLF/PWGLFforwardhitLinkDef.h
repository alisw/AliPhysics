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
#if defined(__CINT__)
/**
 * @file   PWGLFforward2LinkDef.h
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri May 24 09:24:36 2013
 * 
 * @brief  Link specifications
 */
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclasses;

// FMD algorithms 
#pragma link C++ class AliFMDMCHitEnergyFitter+;
#if ROOT_VERSION_CODE < 0x56300 // ROOT_VERSION(5,99,0)
#pragma link C++ class AliFMDMCHitEnergyFitter::RingHistos+;
#endif
#pragma link C++ class AliFMDMCHitEnergyFitterTask+;

#endif
//
// EOF
//
