#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

class Simulator():
   
   def __init__(self):
      self.pwd = os.getenv("PWD")
      print "SIMULATOR"
      
   def mkDirStructure(self):
	 os.system("mkdir -p simulations/single")
	 os.system("mkdir -p simulations/pi0")
   
   def copyFiles(self):
	 os.system("cp ConfigSingle.C simulations/single/Config.C")
	 os.system("cp ConfigPi0.C simulations/pi0/Config.C")
	 os.system("cp sim.C simulations/single/sim.C")
	 os.system("cp sim.C simulations/pi0/sim.C")
   
   def cleanFiles(self):
	 os.system("rm simulations/single/ -rf")
	 os.system("rm simulations/pi0/ -rf")

   def editConfigFiles(self, dophos, doemcal, dotm):

      if dophos == True:
         command = "sed -i \'s/Int_t   iPHOS  =  0/Int_t   iPHOS  =  1/g\' simulations/single/Config.C"
         os.system(command)
      if doemcal == True:
         command = "sed -i \'s/Int_t   iEMCAL =  0/Int_t   iEMCAL =  1/g\' simulations/single/Config.C"
         os.system(command)

      if dotm == True:
         command = "sed -i \'s/Int_t   iTPC   =  0/Int_t   iTPC   =  1/g\' simulations/single/Config.C"
         os.system(command)

   def initSimulation(self):
      self.cleanFiles()
      self.mkDirStructure()
      self.copyFiles()

   def runSimulation(self, nevents, dophos, doemcal, dotm):
      self.editConfigFiles(dophos, doemcal, dotm)
      simargs = str(nevents) + ", " + str(int(dophos)) + ", " + str(int(doemcal)) + ", " + str(int(dotm))
      command = "cd simulations/single/ && aliroot -b -q sim.C\'(" + simargs + ")\'"
      os.system(command)
      os.system("cd " + self.pwd)

class Tester():
   def __init__(self):
      print "TESTER"

   def initTest(self):
      self.pwd = os.getenv("PWD")
      os.system("rm -rf tests/all/")
      os.system("mkdir -p tests/all/single")
      os.system("cp rec_hlt_calo.C tests/all/single/.")
      os.system("cp runAll.C tests/all/single/.")
      os.system("cp read_HLT_ESDs.C tests/all/single/.")

   def run(self, dophos, doemcal, dotm):
      
      runargs = "\"./\", \"./\", " + str(int(dophos)) + ", " + str(int(doemcal)) + ", " + str(int(dotm))
      path = self.pwd + "/simulations/single/."
      command = "cd tests/all/single/;  ln -s " + path + "/raw* .; ln -s " + path + "/GRP . ; aliroot runAll.C\'(" + runargs + ")\'"
      os.system(command)
      os.system("cd " + self.pwd)

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-s", "--simulatedata", action="store_true", dest="simulatedata",
				 default=False, help="Simulate data for the tests")
parser.add_option("-n", "--nevents", dest="nevents",
				 default=100, help="Specify the number of events to simulate. If you change this the reference histograms will be useless")
parser.add_option("", "--no-phos", action="store_false", dest="runphos", 
				 default=True, help="Testing: Don't run PHOS tests. Simulation: Don't simulate for PHOS")
parser.add_option("", "--no-emcal", action="store_false", dest="runemcal",
				 default=True, help="Testing: Don't run EMCAL tests. Simulation: Don't simulate for EMCAL")
parser.add_option("", "--no-trackmatching", action="store_false", dest="runtm",
				 default=True, help="Testing: Don't run the track matcher, Simulation: Don't simulate TPC")

(options, args) = parser.parse_args()

if options.simulatedata:
   simulator = Simulator()
   simulator.initSimulation()
   simulator.runSimulation(options.nevents, options.runphos, options.runemcal, options.runtm)
   exit(0)

else:
   tester = Tester()
   tester.initTest()
   tester.run(options.runphos, options.runemcal, options.runtm)
   exit(0)
    
   
