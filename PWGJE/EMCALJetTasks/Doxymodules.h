/**************************************************************************************
 * Copyright (C) 2019, Copyright Holders of the ALICE Collaboration                   *
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
  * @defgroup PWGJEBASE PWGJE jet framework base task
  * @ingroup EMCALJETFW
  * @brief Basic tasks of the PWGJE jet framework 
  * 
  * Base tasks are the core of the jet trains and used in 
  * (nearly) every jet analysis. Among the base tasks are the 
  * jet finder, the jet taggers, the rho tasks and several 
  * helpers used in the lego train. Note that framework classes 
  * without dependency on fastjet can be found in @ref JETFW.
  * 
  * Note: These are tasks of common interests. User tasks should
  * be documented within the group @ref PWGJEUSER.
  * 
  * See \subpage READMEjetfw for more information
  */

 /**
  * @defgroup PWGJEUSER PWGJE user tasks
  * @ingroup EMCALJETFW
  * @brief User tasks build with the PWGJE jet framework
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
  * 
  * Note: This group is for user tasks. Framework classes (without
  * dependency on fastjet classes can be found here: @ref JETFW. Basic
  * tasks of common interests can be found here: @ref PWGJEBASE
  * 
  * See \subpage READMEJEtasks for more information
  */