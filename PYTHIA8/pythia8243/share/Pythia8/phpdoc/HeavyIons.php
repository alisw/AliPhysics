<html>
<head>
<title>Heavy Ions</title>
<link rel="stylesheet" type="text/css" href="pythia.css"/>
<link rel="shortcut icon" href="pythia32.gif"/>
</head>
<body>

<script language=javascript type=text/javascript>
function stopRKey(evt) {
var evt = (evt) ? evt : ((event) ? event : null);
var node = (evt.target) ? evt.target :((evt.srcElement) ? evt.srcElement : null);
if ((evt.keyCode == 13) && (node.type=="text"))
{return false;}
}

document.onkeypress = stopRKey;
</script>
<?php
if($_POST['saved'] == 1) {
if($_POST['filepath'] != "files/") {
echo "<font color='red'>SETTINGS SAVED TO FILE</font><br/><br/>"; }
else {
echo "<font color='red'>NO FILE SELECTED YET.. PLEASE DO SO </font><a href='SaveSettings.php'>HERE</a><br/><br/>"; }
}
?>

<form method='post' action='HeavyIons.php'>
 
<h2>Heavy Ions</h2> 
 
PYTHIA does not normally handle collisions involving heavy ions, but it 
includes a facility where a model for combining several nucleon-nucleon 
collisions into one heavy ion collision can be implemented. One such model, 
called <a href="#Angantyr" target="page">Angantyr</a>, is provided with PYTHIA 
and is inspired by the old Fritiof program from the Lund group 
[<a href="Bibliography.php#refAnd87" target="page">And87</a>] with recent improvements [<a href="Bibliography.php#refBie16a" target="page">Bie16a</a>] (see below). 
 
<p/> 
 
To allow for the generation of collisions with heavy ions, PYTHIA includes a 
handful of nuclei with PDG numbers on the form 100ZZZAAAI: 
<sup>4</sup>He (1000020040), 
<sup>6</sup>Li (1000030060), 
<sup>12</sup>C (1000060120), 
<sup>16</sup>O (1000080160), 
<sup>63</sup>Cu (1000290630), 
<sup>129</sup>Xe (1000541290), 
<sup>197</sup>Au (1000791970), and 
<sup>208</sup>Pb (1000822080), but  can be added using the function 
<code>ParticleData::addParticle</code>. 
 
<br/><br/><table><tr><td><strong>HeavyIon:mode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 2</code>)</td></tr></table>
This is the master switch for heavy ions, and determines the mode 
of operation of the HeavyIon model. 
<br/>
<input type="radio" name="1" value="1" checked="checked"><strong>1 </strong>:   The heavy ion machinery will be used in case of ion beams.  <br/>
<input type="radio" name="1" value="2"><strong>2 </strong>:   Even collisions without ions are treated using  the heavy ion machinery. (Typically only used  for debugging purposes.)  <br/>
 
<p/> 
 
If <code>HeavyIon:mode</code> is on, the normal initialization in 
<code>Pythia::init()</code> is early on diverted to an object with 
the base class <code>HeavyIons</code> which may instantiate 
secondary <code>Pythia</code> objects needed to generate different types of 
nucleon-nucleon collisions that can be merged together into a full 
heavy ion event. This is all done in the virtual 
<code>HeavyIons::init()</code> function. Subsequent calls to 
<code>Pythia::next()</code> will then also be diverted to the virtual 
function <code>HeavyIons::next()</code> which will be responsible for 
building up the heavy ion collision. The final event will be available in the 
primary <code>Pythia</code> object. 
 
<br/><br/><strong>HeavyIon:showInit</strong>  <input type="radio" name="2" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="2" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Output detailed initialization information from the heavy ion 
model. Typically there will be several <code>Pythia</code> object initialised 
in a heavy ion run, this flag can be used to reduce the amount of output. 
If off, only the output from initialisation of the primary <code>Pythia</code> 
object will be shown 
   
 
<p/> 
 
The <code>HeavyIon</code> class is very simple and flexible and basically only 
specifies that the <code>HeavyIons::init()</code> and 
<code>HeavyIons::next()</code> functions are overridden in a subclass. But 
there are a few additional helper classes that should be generic enough to be 
used by any model implemented. 
 
<ul> 
<li> The <code>Nucleon</code> class represents a single nucleon in a nuclei. 
It can be a proton or a neutron (<code>id()</code>), it has a position in 
impact parameter space (<code>Vec2</code>), 
both absolute (<code>bPos()</code>) and relative to the nuclei 
(<code>nPos()</code>), and optionally it can be in a particular state 
represented by a vector of real numbers which are completely model dependent. 
</li> 
<li> The <code>SubCollision</code> class represents a potential 
nucleon-nucleon  collision between a projectile and a target 
<code>Nucleon</code>. 
</li> 
<li> The <code>NucleusModel</code> class is a base class for implementing a 
model for the distribution in impact parameter 
space of nucleons in a nucleus. There are two ready-made subclasses called 
<code>WoodsSaxonModel</code>, implementing a standard Woods-Saxon 
distribution, and <code>GLISSANDOModel</code>, implementing the  
advanced model in [<a href="Bibliography.php#refBro09" target="page">Bro09</a>,<a href="Bibliography.php#refRyb14" target="page">Ryb14</a>]. 
<br/><br/><strong>HeavyIon:WSHardCore</strong>  <input type="radio" name="3" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="3" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
  In the default Woods-Saxon model for nucleon distributions, assume 
  that there is a minimum distance between nucleons defined by a hard 
  core radius in 
  <code>HeavyIon:WSRh</code>. 
   
<br/><br/><table><tr><td><strong>HeavyIon:WSRh </td><td></td><td> <input type="text" name="4" value="0.9" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.9</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
  The hard core radius in units of fermi, defining the minimum 
  distance between nucleons in a nucleus in the default Woods-Saxon 
  model for nucleon distributions. 
   
<br/><br/><strong>HeavyIon:gaussHardCore</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
  Option to use a Gaussian profile of the hard core instead of a sharp 
  cut-off, inspired by [<a href="Bibliography.php#refBay95" target="page">Bay95</a>]. 
   
 
<br/><br/><table><tr><td><strong>HeavyIon:WSR </td><td></td><td> <input type="text" name="6" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
The radius of a nucleon in units of fermi in the default Woods-Saxon 
model for nucleon distributions. If zero, the size is given by the 
formulae [<a href="Bibliography.php#refRyb14" target="page">Ryb14</a>], based on the number of nucleons in the 
nuclei and whether a hard core is used or not. 
   
 
<br/><br/><table><tr><td><strong>HeavyIon:WSa </td><td></td><td> <input type="text" name="7" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
The <i>skin width</i> of a nucleus in units of fermi in the default 
Woods-Saxon model for nucleon distributions.  If zero, the size is 
given by the numbers in [<a href="Bibliography.php#refRyb14" target="page">Ryb14</a>], based on the number of 
nucleons in the nuclus and whether a hard core is used or not. 
   
</li> 
<li> The <code>ImpactParameterGenerator</code> is used to sample the impact 
parameter space for the colliding nuclei. The base class implements a Gaussian 
sampling, which means that the events produced will always be weighted. Other 
distributions can be implemented in subclasses. 
<br/><br/><table><tr><td><strong>HeavyIon:bWidth </td><td></td><td> <input type="text" name="8" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
  The width in fermi of the distribution by which the impact parameter 
  is sampled. If zero, a suitable width must be guessed by the 
  <code>ImpactParamerGenerator</code> 
 itself. 
   
</li> 
<li> The <code>SubCollisionModel</code> is used to generate individual, 
potential nucleon-nucleon <code>SubCollision</code>s. Two subclasses are 
available, one called <code>NaiveSubCollisionModel</code> which treats 
nucleons as simple black disks, and one implementing a more advanced 
model called <code>DoubleStrikman</code> described below. 
</li> 
<li> The <code>HIInfo</code> class contains information related to the 
generated heavy ion events. 
</li> 
<li> The <code>HIUserHooks</code> class is provided to simplify the 
customization of a model implemented as a <code>HeavyIons</code> subclass. It 
can be used to eg. change the <code>ImpactParamerGenerator</code> used, in a 
way similar to how the <code>UserHooks</code> and <code>MergingHooks</code> 
classes are used. 
</li> 
</ul> 
 
<a name="section0"></a> 
<h3><a name="Angantyr">Angantyr</a> - the default heavy ion model</h3> 
 
The default model in PYTHIA is called Angantyr and is inspired by the old 
Fritiof model [<a href="Bibliography.php#refAnd86" target="page">And86</a>] with improvements described in 
[<a href="Bibliography.php#refBie16a" target="page">Bie16a</a>]. The main idea is to stack parton level events, 
corresponding to individual nucleon-nucleon sub-collisions, on top of 
each other and hadronise them together. 
<br/><strong>Please note:</strong> although it is possible to use 
<?php $filepath = $_GET["filepath"];
echo "<a href='Ropewalk.php?filepath=".$filepath."' target='page'>";?>Rope Hadronisation</a> in heavy ion collisions, 
these two modules have not yet been validated to work properly together. 
Also the parameters in the model have not been properly tuned, so the 
results from running must not be taken as definitive predictions. 
 
<p/> 
 
To determine which projectile nucleon interacts with which target nucleon, 
special care is taken to determine in which way the nucleons interact. In 
a standard <i>Glauber</i> calculations one typicaly only cares about if a 
sub-collision is inelastic or not, but in Angantyr this is divided up, so that 
each inelastic sub-collision can either be single-diffractive, 
double-diffractive or absorptive (ie. non-diffractive). To achieve this, 
Angantyr uses a model with fluctuating radii of the nucleons resulting in a 
fluctuating nucleon-nucleon cross section inspired by the model by Strikman et 
al. [<a href="Bibliography.php#refAlv13" target="page">Alv13</a>]. The model for this includes a number of parameters which 
should be fitted to reproduce inclusive nucleon-nucleon cross sections. To be 
consistent, the values used comes from PYTHIA's internal model of 
<?php $filepath = $_GET["filepath"];
echo "<a href='TotalCrossSections.php?filepath=".$filepath."' target='page'>";?>total cross sections</a>. 
 
<p/> 
 
The default model for nucleon fluctuations has three parameters, the general 
fitting machinery in <code>SubCollisionModel</code> allows for up to eight 
parameters. 
 
<br/><br/><table><tr><td><strong>HeavyIon:SigFitDefPar </td><td></td><td> <input type="text" name="9" value="17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0</strong></code>)</td></tr></table>
These are the default values of the parameters of the 
<code>SubCollisionModel</code> in Angantyr. They will be used as starting 
point when fitting to the inclusive nucleon cross sections. 
   
 
<p/> 
 
The fitting procedure in <code>SubCollisionModel</code> is a kind of genetic 
algorith where a population of parameter values are allowed to evolve for a 
number of generations. In the end the the parameter set in the final 
population  which gives the best inclusive cross sections is picked. 
Eight different cross sections may be fitted to but it is possible to select 
only some of them: 
<br/><br/><table><tr><td><strong>HeavyIon:SigFitErr </td><td></td><td> <input type="text" name="10" value="0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
The relative error assumed in the calculation of goodness of fit 
corresponding to the different cross sections fitted to. The cross 
sections are obtained from the 
<?php $filepath = $_GET["filepath"];
echo "<a href='TotalCrossSections.php?filepath=".$filepath."' target='page'>";?><code>SigmaTotal</code></a> and are given 
as (in order) total, non-diffractive, double diffractive, wounded 
target, wounded projectile, central diffractive, and elastic cross 
sections, and in addition the elastic slope parameter. A relative 
error of zero for one of these cross sections means the corresponding 
cross section not be used in the fit. 
   
 
<a name="anchor1"></a>
<p/><code>mode&nbsp; </code><strong> HeavyIon:SigFitNInt &nbsp;</strong> 
 (<code>default = <strong>100000</strong></code>; <code>minimum = 0</code>)<br/>
The number of integration points used for each parameter setting to calculate 
the cross sections. 
   
 
<a name="anchor2"></a>
<p/><code>mode&nbsp; </code><strong> HeavyIon:SigFitNPop &nbsp;</strong> 
 (<code>default = <strong>20</strong></code>; <code>minimum = 1</code>)<br/>
  The number individuals (parameter settings) in a population in each 
  generation. 
   
 
<a name="anchor3"></a>
<p/><code>mode&nbsp; </code><strong> HeavyIon:SigFitNGen &nbsp;</strong> 
 (<code>default = <strong>20</strong></code>; <code>minimum = 0</code>)<br/>
The number of generation used in the genetic algorithm. If set to zero, no 
fitting will be performed and the values in <code>HeavyIon:SigFitDefPar</code> 
will be used. 
   
 
<br/><br/><table><tr><td><strong>HeavyIon:SigFitFuzz </td><td></td><td> <input type="text" name="11" value="0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 0.5</code>)</td></tr></table>
A parameter determining the probability that an individual parameter setting 
will evolves further away from the best parameter set in each generation. 
   
 
<br/><br/><strong>HeavyIon:SigFitPrint</strong>  <input type="radio" name="12" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="12" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
fitting procedure. If on, extensive information about the fitting will be 
printed. 
   
 
<br/><br/><table><tr><td><strong>Angantyr:CollisionModel  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)</td></tr></table>
The Angantyr model has a couple of option for the SubCollisionModel 
<br/>
<input type="radio" name="13" value="0"><strong>0 </strong>:  A simplified model with fixed nucleon radii.  <br/>
<input type="radio" name="13" value="1" checked="checked"><strong>1 </strong>:  The default model with fluctuating radii and cross  sections.  <br/>
<input type="radio" name="13" value="2"><strong>2 </strong>:  Fluctuating radii and cross sections but different  treatment of opacity.  <br/>
<input type="radio" name="13" value="3"><strong>3 </strong>:  Black disks with no fluctuations, ie. no diffraction.  <br/>
 
<br/><br/><strong>Angantyr:GlauberOnly</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
event generation will stop after SubCollisions has been determined, allowing 
the user to read out the nucleon configuration only. 
   
 
<p/> 
 
After all possible nucleon-nucleon sub-collisions has been determined, 
they are ordered in increasing nucleon-nucleon impact parameter. This 
list is then gone through in order several time. First all absorptive 
sub-collisions are treated. One full nucleon-nucleon non-diffractive 
minimum bias event is generated for each possible absorptive 
sub-colision. These are also ordered in impact parameter. Note that 
one target nucleon can interact absorptively with several target 
nucleons, in a first round only those absorptive sub-collisions 
involving nucleons that have not already interacted absorptively are 
are assigned a non-diffractive event. 
 
<p/> 
 
If PYTHIA is not set up to generate minimum bias events, one or more of the 
generated non-diffractive events will be replaced by events generated with the 
selected signal process, and the cross section reported will be modified 
accordingly. 
 
<p /> 
 
In a second round only those potential absorptive sub-collisions are 
considered where one nucleon has already been assinged a full 
non-diffractive event. In the Angantyr model it is then assumed that 
the other nuclean will contribute to the final state as if it had just 
been diffractively excited. Therefore a corresponding 
single-diffractive event is generated, the elastically scattered beam 
particle is discarded and the rest is added to the previous 
non-diffractive event, shuffling a bit with the kinematics so that the 
total emergy and momentum is conserved. 
 
<p/> 
 
To generate these single-diffraction events to emulate multiple 
absorptive sub-collisions a special <code>Pythia</code> object is 
used. To allow flexibility this object need not have exactly the same 
settings as the the one generating events for normal 
single-diffraction sub-collisions. To manipulate this <code>Pythia</code> 
object a special form of settings can be used. All settings available for 
<?php $filepath = $_GET["filepath"];
echo "<a href='Diffraction.php?filepath=".$filepath."' target='page'>";?><code>Diffraction</code></a>, 
<?php $filepath = $_GET["filepath"];
echo "<a href='MultipartonInteractions.php?filepath=".$filepath."' target='page'>";?> 
  <code>MultipartonInteractions</code></a>, 
<?php $filepath = $_GET["filepath"];
echo "<a href='PDFSelection.php?filepath=".$filepath."' target='page'>";?><code>PDF</code></a>, 
<?php $filepath = $_GET["filepath"];
echo "<a href='TotalCrossSections.php?filepath=".$filepath."' target='page'>";?><code>SigmaDiffractive</code></a> and 
<?php $filepath = $_GET["filepath"];
echo "<a href='PDFSelection.php?filepath=".$filepath."' target='page'>";?><code>PDF</code></a> 
can be set separately for this <code>Pythia</code> object by prefixing 
their names with <code>HI</code>. 
 
As an example, setting <code>HISigmaDiffractive:PomFlux</code> and 
<code>HIPDF:PomSet</code> will set the 
<code>SigmaDiffractive:PomFlux</code> and <code>PDF:PomSet</code> options 
for this <code>Pythia</code> object. 
 
<br/><br/><table><tr><td><strong>Angantyr:SASDmode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>4</strong></code>; <code>minimum = 0</code>; <code>maximum = 4</code>)</td></tr></table>
Determines how to generate single-diffraction events as secondary absorptive 
(SASD) sub-collisions. 
<br/>
<input type="radio" name="15" value="0"><strong>0 </strong>:  Standard singel-diffraction events as speicfied by  <code>HIDiffraction</code> settings above.  <br/>
<input type="radio" name="15" value="1"><strong>1 </strong>:  Always use <code>HIPDF:PomSet = 11</code> and use the  same initial <code>HIMultipartonInteractions:pT0Ref</code> as for  non-diffractive events for the total nucleon-nucleon collision energy,  independent of the mass of the diffractive system.  <br/>
<input type="radio" name="15" value="2"><strong>2 </strong>:  (Experimental) As for option <code>1</code> but also  rescale the pomeron  proton non-diffractive cross section to match the pp non-diffractive one.  <br/>
<input type="radio" name="15" value="3"><strong>3 </strong>:  (Experimental) As for option <code>1</code> but use the full    nucleon-nucleon cross section for the non-diffractive nucleon-Pomeron in the    multiple interaction machinery. Also rescale the Pomeron PDF with the log of    the ratio of maximum and minimum Pomeron-nucleon collision energy.  <br/>
<input type="radio" name="15" value="4" checked="checked"><strong>4 </strong>:  As for option <code>3</code> but no rescaling of the Pomeron   PDF.  <br/>
 
<br/><br/><table><tr><td><strong>Angantyr:impactMode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Determines how to bias non-diffractive minimum-bias sub-collisions in pythia to 
be appropriately central. 
<br/>
<input type="radio" name="16" value="0"><strong>0 </strong>:  If we have N pirmary sub-collisions and Na secondary  sub-collisions, generate N+Na non-diffractive events and pick the N most central.  <br/>
<input type="radio" name="16" value="1"><strong>1 </strong>:  Use UserHooks to force Pythia to produce events with a  particular impact parameter for the N primary sub collisions according to the  generated impact parameter in the SubCollisionModel.  <br/>
<input type="radio" name="16" value="2" checked="checked"><strong>2 </strong>:  As for option <code>1</code> but also the secondary  absorptive sub-collisions have their impact parameter set.  <br/>
 
<br/><br/><table><tr><td><strong>Angantyr:impactFudge </td><td></td><td> <input type="text" name="17" value="0.85" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.85</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 4.0</code>)</td></tr></table>
Multiplicative factor used to compensate for the fact that the 
<code>SubColllisionModel</code> in Angantyr may have a different 
impact parameter profile than what is assumed in the MPI overlap 
calculation in Pythia. 
   
 
<a name="anchor4"></a>
<p/><code>mode&nbsp; </code><strong> Angantyr:SDTries &nbsp;</strong> 
 (<code>default = <strong>1</strong></code>; <code>minimum = 1</code>)<br/>
When adding single diffractive sub-collisions to other 
sub-collisions, there might not be enough energy for the 
diffractive mass. One option here is to  say that the diffractive 
sub-event simply fails, but setting this larger than unity allows 
for regenerating the single diffractive sub-event a number of times 
to see if a small enough diffractive system is produced. 
   
 
<br/><br/><table><tr><td><strong>Angantyr:SDRecoil  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 2</code>)</td></tr></table>
Determines which particles in a primary sub-collision will take the 
recoil when adding single diffractive sub-collisions to other 
sub-collisions. The choice may be overridded by a user-defined 
<code>HIUserHooks::findRecoilers</code> function. 
<br/>
<input type="radio" name="18" value="1" checked="checked"><strong>1 </strong>:   Only elastically scattered nucleons and nucleon  remnants will take recoil.  <br/>
<input type="radio" name="18" value="2"><strong>2 </strong>:   All particles outside the added diffractive system's  rapidity range are considered.  <br/>
 
<br/><br/><strong>Angantyr:SDTest</strong>  <input type="radio" name="19" value="on"><strong>On</strong>
<input type="radio" name="19" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Used in conjunction with <code>HeavyIon:mode = 2</code> and proton 
beams to generate single diffractive events that would be used as 
secondary non-diffractive scatterings in the Angantyr heavy ion model 
for the given nucleon energies. Used for tuning special 
<code>HI</code>-prefixed parameters of the secondary absorptive 
sub-collisions. 
   
 
<br/><br/><table><tr><td><strong>Angantyr:SDTestB </td><td></td><td> <input type="text" name="20" value="-1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.0</strong></code>)</td></tr></table>
In conjunction with <code>Angantyr:SDTest = on</code> and 
<code>Angantyr:impactMode = 2</code> only pick diffractive 
events with a particular impact parameter (as defined by the scaled value 
given in <code>Info::bMPI()</code>). If negative, the standard impact 
parameter distribution is used. 
   
 
<p/> 
 
After all absorptive sub-collisions have been dealt with, the 
diffractive and elastic sub-collisions are dealt with in a similar 
way. In the end there will be a number of parton level events which 
are finally stacked together, and then hadronised. Finally nucleus 
remnants constructed from the non-interacting nucleans, are added to 
complete the full nucleaus-nucleus collision. 
 
<input type="hidden" name="saved" value="1"/>

<?php
echo "<input type='hidden' name='filepath' value='".$_GET["filepath"]."'/>"?>

<table width="100%"><tr><td align="right"><input type="submit" value="Save Settings" /></td></tr></table>
</form>

<?php

if($_POST["saved"] == 1)
{
$filepath = $_POST["filepath"];
$handle = fopen($filepath, 'a');

if($_POST["1"] != "1")
{
$data = "HeavyIon:mode = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "on")
{
$data = "HeavyIon:showInit = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "on")
{
$data = "HeavyIon:WSHardCore = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "0.9")
{
$data = "HeavyIon:WSRh = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "HeavyIon:gaussHardCore = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "0.0")
{
$data = "HeavyIon:WSR = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "0.0")
{
$data = "HeavyIon:WSa = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "0.0")
{
$data = "HeavyIon:bWidth = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0")
{
$data = "HeavyIon:SigFitDefPar = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0")
{
$data = "HeavyIon:SigFitErr = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "0.2")
{
$data = "HeavyIon:SigFitFuzz = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "on")
{
$data = "HeavyIon:SigFitPrint = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "1")
{
$data = "Angantyr:CollisionModel = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "Angantyr:GlauberOnly = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "4")
{
$data = "Angantyr:SASDmode = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "2")
{
$data = "Angantyr:impactMode = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "0.85")
{
$data = "Angantyr:impactFudge = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "1")
{
$data = "Angantyr:SDRecoil = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "off")
{
$data = "Angantyr:SDTest = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "-1.0")
{
$data = "Angantyr:SDTestB = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
