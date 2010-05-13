#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

class Simulator():
   
   def __init__(self):
      self.pwd = os.getenv("PWD")
      
   def mkDirStructure(self):
	 os.system("mkdir -p simulations/single")
	 os.system("mkdir -p simulations/pi0")
   
   def copyFiles(self):
	 os.system("cp ConfigSingle.C simulations/single/Config.C")
	 os.system("cp ConfigPi0.C simulations/pi0/Config.C")
	 os.system("cp sim.C simulations/single/.")
	 os.system("cp sim.C simulations/pi0/.")
   
   def cleanFiles(self):
	 os.system("rm simulations/single/ -rf")
	 os.system("rm simulations/pi0/ -rf")
   
   def initSimulation(self):

      self.cleanFiles()
      self.mkDirStructure()
      self.copyFiles()

   def runSimulation(self):
      
      os.system("cd simulations/single/ && aliroot -b -q sim.C\(100\)")
      os.system("cd " + self.pwd)
   

class Test(object):
   def __init__(self):
      print "Test"

class PhosTest(Test):
   def __init__(self):
      super(PhosTest, self).__init__()

   def initTest(self):
      
      self.pwd = os.getenv("PWD")
      os.system("rm -rf tests/phos/")
      os.system("mkdir -p tests/phos/single")
      os.system("cp rec_hlt_calo_phos.C tests/phos/single/.")
      os.system("cp runPhos.C tests/phos/single/.")
      os.system("cp read_HLT_ESDs.C tests/phos/single/.")
      
   def run(self):
      
      path = self.pwd + "/simulations/single/."
      os.system("cd tests/phos/single/;  ln -s " + path + "/raw* .; ln -s " + path + "/GRP . ; aliroot runPhos.C;")
      
      os.system("cd " + self.pwd)
      


class EmcalTest(Test):
   def __init__(self):
      super(EmcalTest, self).__init__()
      print "EmcalTest"
      
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-s", "--simulatedata", action="store_true", dest="simulatedata",
				 default=False, help="Simulate data for the tests")
parser.add_option("", "--nophos", action="store_false", dest="runphos", 
				 default=True, help="Don't run PHOS tests")
parser.add_option("", "--noemcal", action="store_false", dest="runemcal",
				 default=True, help="Don't run EMCAL tests")

(options, args) = parser.parse_args()

if options.simulatedata:
   simulator = Simulator()
   simulator.initSimulation()
   simulator.runSimulation()
   
if options.runphos:
   phostest = PhosTest()
   phostest.initTest()
   phostest.run()
   
if options.runemcal:
   emcalTest = EmcalTest()
   #emcalTest.run()
   
   