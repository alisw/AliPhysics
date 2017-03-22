/*
 * $Id$
 *
/************************************************************************
**
** ALICE HLT project
** Copyright (c) 2005
**
** This file is property of and copyright by the Experimental Nuclear 
** Physics Group, Dep. of Physics and Technology
** University of Bergen, Norway, 2004
** This file has been written by Matthias Richter,
** Matthias.Richter@ift.uib.no
**
** Permission to use, copy, modify and distribute this software and its  
** documentation strictly for non-commercial purposes is hereby granted  
** without fee, provided that the above copyright notice appears in all  
** copies and that both the copyright notice and this permission notice  
** appear in the supporting documentation. The authors make no claims    
** about the suitability of this software for any purpose. It is         
** provided "as is" without express or implied warranty.                 
**
*************************************************************************/

/** @file   mainpage.c
    @author Matthias Richter
    @date   
    @brief  Title page documentation. */

/** @mainpage ALICE HLT analysis framework

    @section intro Introduction

    This package contains the analysis code of the ALICE High Level Trigger
    system (HLT). The system entails a very large processing farm, designed
    for an anticipated input data stream of 25 GB/s and allows on-line data
    processing at the full input rate and efficient data rate reduction.

    A generic communication framework has been developed based on
    the publisher-subscriber principle, in general referred to be the 
    Publisher-Subscriber Framework (PubSub), which allows an arbitrary
    connectivity metric of processing elements across the underlying
    network interface.

    The Framwork implies serveral rules on the component implementation.

    @section overview Overview

    Two running modes can be destinguished: 

    - The on-line mode entails the Publisher-Subscriber Framework to run the
      analysis and processing components on the different nodes of the HLT
      cluster farm. 
    - The off-line mode integrates the HLT processing components into the ALICE
      Off-line analysis framework AliRoot.

    A common interface for HLT processing components has been designed to
    run the components from either the on-line or off-line analysis framework
    without changes. The interface adapts the component to the needs of the
    on-line processing within the PubSub and allows the developer at the same
    time to use AliRoot for easy development, debugging, and benchmarking. 
    Results can be compared directly.

    The approach is base on shared libraries and an abstract pure-C interface 
    to integrate the analysis code into the PubSub. The libraries containing the
    analysis code are built within AliRoot. The PubSub is the pure on-line data
    transportation framework which uses exactly the same executables (shared 
    libraries) for analysis components as the off-line framework.

    \image html  HLT-AliRoot-Integration_overview.png "Overview"
    \image latex HLT-AliRoot-Integration_overview.eps "Overview" width=14cm

    @section toc Sections
    
    @subsection toc_developer Developer Section

    - @ref alihlt_component <br>
      Description of the component interface and guidelines for component implementation.
      
    - @ref alihlt_system <br>
      Description of the HLT integration into offline AliRoot. <br>
      (AliRoot simulation, AliRoot reconstruction and analysis.)

    - @ref alihlt_tutorial <br>
      HLT examples and tutorials.

    @subsection toc_experts Expert Section

    - @ref alihlt_wrapper_interface <br>
      Description of the wrapper interface for external utilization of the module.

    @section libraries Detector/Module libraries

    - @ref alihlt_modules

    @section links Related links on the web
    
    - <a class="el" href="http://aliceinfo.cern.ch">
          The ALICE Experiment </a> 
    - <a class="el" href="https://twiki.cern.ch/twiki/bin/viewauth/ALICEHLT/">
          The ALICE High Level Trigger web pages </a> 
    - <a class="el" href="http://aliceinfo.cern.ch/Offline">
          The ALICE Offline web pages </a> 

*/

#error Not for compilation
//
// EOF
//
