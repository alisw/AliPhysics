/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <cstdlib>
#include <bitset>
#include <iostream>
#include <sstream>
#include <string>

#include "AliEMCALTriggerSTUDCSConfig.h"

namespace EMCAL {
  namespace Base {
    namespace Tests {
      namespace TestAliEMCALTriggerSTUDCSConfig {

        /**
         * @brief Apply reference configuration
         * 
         * Reference configuration taken from pp 2016
         * 
         * @param testobject Object for test to be configured
         */
        void ConfigureReference(AliEMCALTriggerSTUDCSConfig &testobject){
          testobject.SetG(0,0,0);
          testobject.SetG(1,0,0);
          testobject.SetG(2,0,115);
          testobject.SetG(0,1,0);
          testobject.SetG(1,1,0);
          testobject.SetG(2,1,51);
          testobject.SetJ(0,0,0);
          testobject.SetJ(1,0,0);
          testobject.SetJ(2,0,255);
          testobject.SetJ(0,1,0);
          testobject.SetJ(1,1,0);
          testobject.SetJ(2,1,204);
          testobject.SetPatchSize(2);
          testobject.SetFw(0x2A012);
          testobject.SetMedianMode(0);
          testobject.SetRegion(0xffffffff);
          for(int i = 0; i < 4; i++) testobject.SetPHOSScale(i, 0);
        }

        /**
         * @brief Test for operator== on itself
         * 
         * Tests whether operator== returns true in case the object is tested
         * against itself.
         * 
         * @return true Test passed (operator== returning true)
         * @return false Test failed (operator== returning false)
         */
        bool TestEqualSelf(){
          AliEMCALTriggerSTUDCSConfig test;
          ConfigureReference(test);
          return test == test;
        }

        /**
         * @brief Test for operator== on same object
         * 
         * Tests whether operator== returns true in case both
         * objects have the same content.
         * 
         * @return true Test passed (operator== returning true)
         * @return false Test failed (operator== returning false)
         */
        bool TestEqualTrue(){
          AliEMCALTriggerSTUDCSConfig test1, test2;
          ConfigureReference(test1);
          ConfigureReference(test2);
          return test1 == test2;
        }

        /**
         * @brief Test for operator== on different objects
         * 
         * Tests whether the operator== returns false if at least one setting
         * is different. For this operator== is tested with multiple objects
         * based on a reference setting where only one parameter is changed at
         * the time.
         * 
         * @return true Test passed (operator== returns false in all cases)
         * @return false Test failed (operator== returns true at least once)
         */
        bool TestEqualFalse(){
          AliEMCALTriggerSTUDCSConfig ref, test;
          ConfigureReference(ref);
          std::bitset<16> testresults;

          // Variation Gamma
          ConfigureReference(test);
          test.SetG(2, 0, 77);
          test.SetG(2, 1, 51);
          testresults.set(0, ref == test);

          // Variation Jet
          ConfigureReference(test);
          test.SetJ(2, 0, 191);
          test.SetJ(2, 1, 128);
          testresults.set(1, ref == test);

          // Variation patch size
          ConfigureReference(test);
          test.SetPatchSize(0);
          testresults.set(2, ref == test);

          // Variation region
          ConfigureReference(test);
          test.SetRegion(0xffffff7f);
          testresults.set(3, ref == test);

          // Variation fw
          ConfigureReference(test);
          test.SetFw(0x1A012);
          testresults.set(4, ref == test);

          // Variation median mode
          ConfigureReference(test);
          test.SetMedianMode(1);
          testresults.set(5, ref == test);

          // Variation PHOS scale
          ConfigureReference(test);
          test.SetPHOSScale(0, 1);
          test.SetPHOSScale(1, 2);
          test.SetPHOSScale(2, 1);
          test.SetPHOSScale(3, 0);
          testresults.set(6, ref == test);

          return !testresults.any();
        }

        /**
         * @brief Test for the stream operator
         * 
         * Test if operator<< for a reference configuration produces 
         * the expected reference string. Test is implemented using a streaming
         * operator.
         * 
         * @return true Test passed (operator<< to string produces reference string)
         * @return false Test failed (operator<< to string produces different string)
         */
        bool TestStream() {
          std::string reference = std::string("Gamma High: (0, 0, 115)\nGamma Low:  (0, 0, 51)\nJet High:   (0, 0, 255)\nJet Low:    (0, 0, 204)\n")
                                + std::string("GetRawData: 1, Region: -1, Median: 0Firmware: 2a012, PHOS Scale: (0, 0, 0, 0)\n");
          
          AliEMCALTriggerSTUDCSConfig test;
          ConfigureReference(test);
          std::stringstream testmaker;
          testmaker << test;
          return testmaker.str() == reference;
        }
      }
    }
  }
}

/**
 * @brief Test runner for tests of class AliEMCALTriggerSTUDCSConfig
 * 
 * Running all tests implemented for the TRU DCS Config
 * 
 * @param argc Number of commamnd line arguments (not handled)
 * @param argv Array of command line arguments (not handled)
 * @return int 0 if the test is successful, 1 if not
 */
int main(int argc, char **argv) {
  int npass(0);
  if(EMCAL::Base::Tests::TestAliEMCALTriggerSTUDCSConfig::TestEqualSelf()){
    npass++;
  } else {
    std::cerr << "Failure eq_self" << std::endl;
  }
  if(EMCAL::Base::Tests::TestAliEMCALTriggerSTUDCSConfig::TestEqualTrue()){
    npass++;
  } else {
    std::cerr << "Failure_eq_true" << std::endl;
  }
  if(EMCAL::Base::Tests::TestAliEMCALTriggerSTUDCSConfig::TestEqualFalse()){
    npass++;
  } else {
    std::cerr << "Failure_eq_false" << std::endl;
  }
  if(EMCAL::Base::Tests::TestAliEMCALTriggerSTUDCSConfig::TestStream()){
    npass++;
  } else {
    std::cerr << "Failure_stream" << std::endl;
  }
  std::cout << "Passing " << npass << "/4" << std::endl;
  if(npass == 4) return EXIT_SUCCESS;
  return EXIT_FAILURE;
}