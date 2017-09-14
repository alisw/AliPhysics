/**************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/

/**
 *\defgroup EMCALFW The EMCAL analysis framework
 *\brief An analysis framework for EMCAL related studies
 *
 *See \subpage READMEemcfw
 */

/**
 * \defgroup EMCALCOREFW The core EMCAL framework
 * \ingroup EMCALFW
 * \brief The core EMCAL framework
 *
 * See \subpage READMEcorefw
 */

/**
 * \defgroup EMCALTRGFW EMCAL trigger framework
 * \ingroup EMCALFW
 * \brief The EMCAL trigger framework
 *
 * See \subpage READMEtrgfw
 */

/**
 * \defgroup EMCALFWTASKS EMCAL framework tasks
 * \ingroup EMCALFW
 * \brief EMCAL framework tasks
 *
 * 
 */
 
 /**
 * \defgroup EMCALJETFW Jet finding framework
 * \ingroup EMCALFW
 * \brief Framework for jet finding and analysis of jet-related observables
 *
 * See \subpage READMEjetfw
 */
 
  /**
 * \defgroup PWGJETASKS PWG-JE user tasks 
 * \ingroup EMCALJETFW
 * \brief User tasks using the EMCal jet framework.
 * 
 * The EMCal jet framework is described in \subpage READMEjetfw.
 * # How to document YOUR task
 * Documenting your class with doxygen is easy! Follow the instructions in [Dario's page](https://dberzano.github.io/alice/doxygen/#documenting_a_c_class).
 * 
 * To include your class into the _Module_ "PWG-JE user tasks" add to the comment in the header file the keywork \ <strong>`ingroup`</strong> (no space) followed by the name of the module, in this case PWGJETASKS.
 * 
 * The keyword \ <strong>`param`</strong> allows to document parameters, and \ <strong>`return`</strong> the return value of the method
 * 
 * You can link any doxygen documentation page to your class using \ <strong>`subpage`</strong> followed by the name of the file.
 * 
 * Some of the classes have already been documented (see above), check them out!
 *
 * If you want to add a more detailed description, you can use this page: \subpage READMEJEtasks. You can edit it by modifying the file $ALICE_PHYSICS/../src/PWGJE/doc/READMEJEtasks.txt that can be linked to your class documentation by using \ <strong>`subpage` READMEJEtasks </strong>.
 */
