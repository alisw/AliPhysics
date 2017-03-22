<html>
<head>
<title>SUSY Les Houches Accord</title>
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

<form method='post' action='SUSYLesHouchesAccord.php'>

<h2>SUSY Les Houches Accord</h2>

The PYTHIA 8 program does not contain an internal spectrum calculator
(a.k.a. RGE package) to provide supersymmetric couplings, mixing angles,
masses and branching ratios. Thus the SUSY Les Houches Accord (SLHA)
[<a href="Bibliography.php" target="page">Ska04</a>][<a href="Bibliography.php" target="page">All08</a>] is the only way of
inputting SUSY models, and SUSY processes (see
the <?php $filepath = $_GET["filepath"];
echo "<a href='SUSYProcesses.php?filepath=".$filepath."' target='page'>";?>SUSYProcesses</a> page) 
cannot be run unless such an input has taken place. 

<p/>
The SLHA input format can also be extended for use with more general BSM
models, beyond SUSY. Information specific to  how to use the SLHA
interface for generic BSM models is collected below,
under <a href="#generic">Using SLHA for generic BSM Models</a>, with
more elaborate explanations and examples in [<a href="Bibliography.php" target="page">Des11</a>]. 

<p/>
Most of the SUSY implementation in PYTHIA 8 is compatible with both the 
SLHA1 [<a href="Bibliography.php" target="page">Ska04</a>] and SLHA2 [<a href="Bibliography.php" target="page">All08</a>] 
conventions (with some limitations for the NMSSM 
in the latter case). Internally, PYTHIA 8 uses the 
SLHA2 conventions and translates SLHA1 input to these when necessary. 
See the section on SUSY Processes and [<a href="Bibliography.php" target="page">Des11</a>] for more
information. Note that PYTHIA assumes that a spectrum is either fully SHLA1 
or fully SLHA2 compliant. Mixing of the two standards is discouraged, as 
this can lead to ambiguities and inconsistencies.  

<p/>
When reading LHEF files, Pythia automatically looks for SLHA information
between <code>&lt;slha&gt;...&lt;/slha&gt;</code> tags in the header of such
files. When running Pythia without LHEF input (or if reading an LHEF
file that does not contain SLHA information in the header), a separate 
file containing SLHA information may be specified using 
<code>SLHA:file</code> (see below). 

<p/>
Normally the LHEF would be in uncompressed format, and thus human-readable
if opened in a text editor. A possibility to read gzipped files has 
been added, based on the Boost and zlib libraries, which therefore
have to be linked appropriately in order for this option to work.
See the <code>README</code> file in the main directory for details 
on how to do this. 

<p/>
Finally, the SLHA input capability can of course also be used to input 
SLHA-formatted <code>MASS</code> and <code>DECAY</code> tables for 
other particles, such as the Higgs boson, furnishing a less 
sophisticated but more universal complement to the
standard PYTHIA 8-specific methods for inputting such information (for the
latter, see the section on <?php $filepath = $_GET["filepath"];
echo "<a href='ParticleData.php?filepath=".$filepath."' target='page'>";?>Particle Data</a>
and the <?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDataScheme.php?filepath=".$filepath."' target='page'>";?>scheme</a> to modify it). This 
may at times not be desirable, so a few options can be used to curb the right 
of SLHA to overwrite particle data. 
Conversely, it is sometimes useful to allow the user to modify 
eg a mass parameter relative to its value in the SLHA spectrum. 
This is normally not permitted (the SLHA spectrum is normally self-consistent 
and should not be modified), but an option for allowing it is provided.  

<p/>
The reading-in of information from SLHA or LHEF files is handled by the
<code>SusyLesHouches</code> class, while the subsequent calculation of 
derived quantities of direct application to SUSY processes is done in the
<code>CoupSUSY</code>, <code>SigmaSUSY</code>,
and <code>SUSYResonanceWidths</code> classes.

<h3>SLHA Switches and Parameters</h3>

<p/><code>mode&nbsp; </code><strong> SLHA:readFrom &nbsp;</strong> 
 (<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)<br/>
Controls from where SLHA information is read.
<br/><code>option </code><strong> 0</strong> : is not read at all. Useful when SUSY is not simulated
and normal particle properties should not be overwritten.  
<br/><code>option </code><strong> 1</strong> : read in from the <code>&lt;slha&gt;...&lt;/slha&gt;</code> 
block of a LHEF, if such a file is read during initialization, and else
from the <code>SLHA:file</code> below.  
<br/><code>option </code><strong> 2</strong> : read in from the <code>SLHA:file</code> below.  
  

<br/><br/><table><tr><td><strong>SLHA:file  </td><td></td><td> <input type="text" name="1" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
Name of an SLHA (or LHEF) file containing the SUSY/BSM model definition, 
spectra, and (optionally) decay tables. Default <code>void</code>
signals that no such file has been assigned.
  

<br/><br/><strong>SLHA:keepSM</strong>  <input type="radio" name="2" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="2" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Some programs write SLHA output also for SM particles where normally
one would not want to have masses and decay modes changed unwittingly. 
Therefore, by default, known SM particles are ignored in SLHA files. 
To be more specific, particle data for identity codes in the ranges 
1 - 24 and 81 - 999,999 are ignored. Notably this includes <i>Z^0</i>, 
<i>W^+-</i> and <i>t</i>. The SM Higgs is modified by the SLHA input, 
as is other codes in the range 25 - 80 and 1,000,000 - . If you 
switch off this flag then also SM particles are modified by SLHA input.
  

<br/><br/><table><tr><td><strong>SLHA:minMassSM </td><td></td><td> <input type="text" name="3" value="100.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>100.0</strong></code>)</td></tr></table>
This parameter provides an alternative possibility to ignore SLHA input 
for all particles with identity codes below 1,000,000 (which mainly
means SM particle, but also includes e.g. the Higgs bosons in 
two-Higgs-doublet scenarios) whose default masses in PYTHIA lie below 
some threshold value, given by this parameter. The default value of 
100.0 allows SLHA input to modify the top quark, but not, e.g., the 
<i>Z^0</i> and <i>W^+-</i> bosons. 
  

<br/><br/><strong>SLHA:allowUserOverride</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Flag to set whether the user is allowed to modify the parameters read 
from an SLHA spectrum. Is normally kept <code>off</code> to preserve the 
internal self-consistency of SLHA spectra. If this flag is switched 
<code>on</code>, the mass values read from the SLHA block MASS are 
allowed to be modified by the user, using PYTHIA's standard 
<code>readString</code> and related methods. 
  

<h3>SLHA DECAY Tables</h3>

<br/><br/><strong>SLHA:useDecayTable</strong>  <input type="radio" name="5" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="5" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Switch to choose whether to read in SLHA <code>DECAY</code> tables or not. 
If this switch is set to off, PYTHIA will ignore any decay tables found 
in the SLHA file, and all decay widths will be calculated internally by
PYTHIA. If switched on, SLHA decay tables will be read in, and will
then supersede PYTHIA's internal calculations, with PYTHIA only
computing the decays for particles for which no SLHA decay table is
found. (To set a particle stable, you may either omit an SLHA 
<code>DECAY</code> table for it and then  
use PYTHIA's internal <code>id:MayDecay</code> switch for that
particle, or you may include an SLHA <code>DECAY</code> table for it, 
with the width set explicitly to zero.)
  

<br/><br/><table><tr><td><strong>SLHA:minDecayDeltaM </td><td></td><td> <input type="text" name="6" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>)</td></tr></table>
This parameter sets the smallest allowed mass difference (in GeV,
between the mass of the mother and the sum of the daughter masses) 
for a decay mode in a DECAY table to be switched on inside PYTHIA. The
default is to require at least 1 GeV of open phase space, but this can
be reduced (at the user's risk) for instance to be able to treat
decays in  models with very small mass splittings.
  

<h3>Internal SLHA Variables</h3>

<p/><code>mode&nbsp; </code><strong> SLHA:verbose &nbsp;</strong> 
 (<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)<br/>
Controls amount of text output written by the SLHA interface, with a
value of 0 corresponding to the most quiet mode.
  

The following variables are used internally by PYTHIA as local copies
of SLHA information. User changes will generally have no effect, since
these variables will be reset by the SLHA reader during initialization.

<br/><br/><strong>SLHA:NMSSM</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Corresponds to SLHA block MODSEL entry 3.
  

<a name="generic"></a>
<h2>Using SLHA for generic BSM Models</h2>

</p>
Using the <code>QNUMBERS</code> extension [<a href="Bibliography.php" target="page">Alw07</a>], the SLHA
can also be used to define new particles, with arbitrary quantum
numbers. This already serves as a useful way to introduce new
particles and can be combined with <code>MASS</code> and 
<code>DECAY</code> tables in the usual
way, to generate isotropically distributed decays or even chains of
such decays. (If you want something better than isotropic, sorry, you'll
have to do some actual work ...)
</p>

</p>
A more advanced further option is to make use of the possibility
in the SLHA to include user-defined blocks with arbitrary
names and contents. Obviously, standalone 
PYTHIA 8 does not know what to do with such information. However, it
does not throw it away either, but instead stores the contents of user
blocks as strings, which can be read back later, with the user
having full control over the format used to read the individual entries. 
</p>

<p>
The contents of both standard and user-defined SLHA blocks can be accessed 
in any class inheriting from PYTHIA 8's <code>SigmaProcess</code>
class (i.e., in particular, from any semi-internal process written by
a user), through its SLHA pointer, <code>slhaPtr</code>, by using the 
following methods: 
<a name="method1"></a>
<p/><strong> &nbsp;</strong> <br/>
  bool slhaPtr->getEntry(string blockName, double& val); 
  
<strong> &nbsp;</strong> <br/>
  bool slhaPtr->getEntry(string blockName, int indx, double& val); 
  
<strong> &nbsp;</strong> <br/>
  bool slhaPtr->getEntry(string blockName, int indx, int jndx, double& val); 
  
<strong> &nbsp;</strong> <br/>
  bool slhaPtr->getEntry(string blockName, int indx, int jndx, int
  kndx, double& val); 
  
</p>

<p>
This particular example assumes that the user wants to read the
entries (without index, indexed, matrix-indexed, or 3-tensor-indexed, 
respectively) in the user-defined block <code>blockName</code>, 
and that it should be interpreted as 
a <code>double</code>. The last argument is templated, and hence if
anything other than a <code>double</code> is desired to be read, the
user has only to give the last argument a different type. 
If anything went wrong (i.e., the block doesn't
exist, or it doesn't have an entry with that index, or that entry
can't be read as a double), the method returns false; true
otherwise. This effectively allows to input completely arbitrary
parameters using the SLHA machinery, with the user having full control
over names and conventions. Of course, it is then the user's
responsibility to ensure complete consistency between the names and
conventions used in the SLHA input, and those assumed in any
user-written semi-internal process code. 
</p>

<p>
Note that PYTHIA 8 always initializes at least 
the SLHA blocks MASS and SMINPUTS, starting from its internal 
SM parameters and particle data table values (updated to take into
account user modifications). These blocks can therefore be accessed 
using the <code>slhaPtr->getEntry()</code> methods even in the absence 
of SLHA input. 
Note: in the SMINPUTS block, PYTHIA outputs physically correct
(i.e., measured) values of <i>GF</i>, <i>m_Z</i>, and 
<i>alpha_EM(m_Z)</i>. However, if one attempts to compute, e.g., 
the W mass, at one loop from these quantities, a value of 79 GeV results, 
with a corresponding value for the weak mixing angle. We advise to 
instead take the physically measured W mass from block MASS, and 
recompute the EW parameters as best suited for the application at hand.
</p>

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

if($_POST["1"] != "void")
{
$data = "SLHA:file = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "on")
{
$data = "SLHA:keepSM = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "100.0")
{
$data = "SLHA:minMassSM = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "SLHA:allowUserOverride = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "on")
{
$data = "SLHA:useDecayTable = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "1.0")
{
$data = "SLHA:minDecayDeltaM = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "SLHA:NMSSM = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->


