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

#include "AliEMCALTriggerTRUDCSConfig.h"

namespace EMCAL {
  namespace Base {
    namespace Tests {
      namespace TestAliEMCALTriggerTRUDCSConfig {

        /**
         * @brief Apply reference configuration
         * 
         * Reference configuration taken from pp 2016
         * 
         * @param testobject Object for test to be configured
         */
        void ConfigureReference(AliEMCALTriggerTRUDCSConfig &testobject){
          testobject.SetSELPF(7711);
	        testobject.SetL0SEL(1);
	        testobject.SetL0COSM(100);
	        testobject.SetGTHRL0(132);
	        testobject.SetMaskReg(1024, 0);
	        testobject.SetMaskReg(0, 1);
	        testobject.SetMaskReg(512, 2);
	        testobject.SetMaskReg(31985, 3);
	        testobject.SetMaskReg(0, 4);
	        testobject.SetMaskReg(0, 5);
	        testobject.SetRLBKSTU(0);
	        testobject.SetFw(0x21);
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
          AliEMCALTriggerTRUDCSConfig testobject;
          ConfigureReference(testobject);

          // Make test          
          return testobject == testobject;
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
          AliEMCALTriggerTRUDCSConfig test1, test2;
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
          AliEMCALTriggerTRUDCSConfig ref, testobject;
          std::bitset<8> testresult;      // bit set in case test fails (operator returns true)
          testresult.reset();

          // Configure reference          
          ConfigureReference(ref);

          // Variation peak finder
          ConfigureReference(testobject);
          testobject.SetSELPF(7000);
          testresult.set(0, ref == testobject);

          // Variation L0 Algorithm
          ConfigureReference(testobject);
          testobject.SetL0SEL(2);
          testresult.set(1, ref == testobject);

          // Variation cosmic threshold
          ConfigureReference(testobject);
          testobject.SetL0COSM(1000);
          testresult.set(2, ref == testobject);

          // Variation global threshold
          ConfigureReference(testobject);
          testobject.SetGTHRL0(184);
          testresult.set(3, ref == testobject);

          // Variation Rollback
          ConfigureReference(testobject);
          testobject.SetRLBKSTU(1);
          testresult.set(4, ref == testobject);

          // Variation Firmware
          ConfigureReference(testobject);
          testobject.SetFw(0x11);
          testresult.set(5, ref == testobject);

          // Variation Mask
          ConfigureReference(testobject);
          testobject.SetMaskReg(768,0);
          testobject.SetMaskReg(15,1);
          testobject.SetMaskReg(37632,2);
          testobject.SetMaskReg(63,3);
          testobject.SetMaskReg(0,4);
          testobject.SetMaskReg(208,5);
          testresult.set(5, ref == testobject);

          return !testresult.any();
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
        bool TestStream(){
          std::string reference = std::string("SELPF: 1e1f, L0SEL: 1, L0COSM: 100, GTHRL0: 132, RLBKSTU: 0, FW: 21\n")
                                + std::string("Reg0: 00000000000000000000010000000000\nReg1: 00000000000000000000000000000000\n")
                                + std::string("Reg2: 00000000000000000000001000000000\nReg3: 00000000000000000111110011110001\n")
                                + std::string("Reg4: 00000000000000000000000000000000\nReg5: 00000000000000000000000000000000\n");
          
          AliEMCALTriggerTRUDCSConfig test;
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
 * @brief Test runner for tests of class AliEMCALTriggerTRUDCSConfig
 * 
 * Running all tests implemented for the TRU DCS Config
 * 
 * @param argc Number of commamnd line arguments (not handled)
 * @param argv Array of command line arguments (not handled)
 * @return int 0 if the test is successful, 1 if not
 */
int main(int argc, char **argv) {
  int npass(0);
  if(EMCAL::Base::Tests::TestAliEMCALTriggerTRUDCSConfig::TestEqualSelf()){
    npass++;
  } else {
    std::cerr << "Failure eq_self" << std::endl;
  }
  if(EMCAL::Base::Tests::TestAliEMCALTriggerTRUDCSConfig::TestEqualTrue()){
    npass++;
  } else {
    std::cerr << "Failure_eq_true" << std::endl;
  }
  if(EMCAL::Base::Tests::TestAliEMCALTriggerTRUDCSConfig::TestEqualFalse()){
    npass++;
  } else {
    std::cerr << "Failure_eq_false" << std::endl;
  }
  if(EMCAL::Base::Tests::TestAliEMCALTriggerTRUDCSConfig::TestStream()){
    npass++;
  } else {
    std::cerr << "Failure_stream" << std::endl;
  }
  std::cout << "Passing " << npass << "/4" << std::endl;
  if(npass == 4) return EXIT_SUCCESS;
  return EXIT_FAILURE;
}