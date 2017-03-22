<html>
<head>
<title>POWHEG Merging</title>
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

<form method='post' action='POWHEGMerging.php'>
 
<h2>POWHEG Merging</h2> 
 
POWHEG [<a href="Bibliography.php" target="page">Nas04</a>] in its character is very much like a parton shower, 
with a Sudakov factor arising from the ordering of emissions. Both 
POWHEG-BOX [<a href="Bibliography.php" target="page">Ali10</a>] and PYTHIA are based on a combined evolution 
of ISR and FSR in <i>pT</i>-related "hardness" variables, and thus are 
kindred spirits. The hardness definitions differ, however. Frequently we 
will therefore need to distinguish between POWHEG-hardness and 
PYTHIA-hardness in the following. 
 
<p/> 
The simplest merging solution, of continuing the PYTHIA shower at the LHA 
<code>scale</code> hardness where POWHEG leaves off, is obtained if you 
set <code>SpaceShower:pTmaxMatch = 1</code> and 
<code>TimeShower:pTmaxMatch = 1</code>. But then mismatches are bound to 
happen: some regions may be doublecounted, while others may not be counted 
at all. Depending on the choice of hardness, such mismatches might be small. 
 
<p/> 
There are no guarantees, however, so a (hopefully) more accurate merging 
scheme is coded up in the <code>include/Pythia8Plugins/PowHegHooks.h</code> 
file, with a realistic user example in the <code>examples/main31</code> 
files. Here we would like to discuss the (POWHEG-specific) input settings 
for <code>main31.cc</code>, see <code>main31.cmnd</code>, and attempt to 
give some recommendations on how to use the main program to perform a 
matching of POWHEG-BOX with PYTHIA 8. 
 
<p/> 
POWHEG-BOX inputs contain Born-like events (with no resolved emission) and 
Real-type events (containing an additional parton). The mismatch between 
POWHEG-hardness and PYTHIA-hardness can be minimised if the PYTHIA shower 
knows 
<br/>a) The POWHEG-hardness criterion (through which the separation of Born- 
and Real-like events is defined), and 
<br/>b) The POWHEG-hardness value (which separates Born- and Real-like 
events). 
<br/>If these definitions are known, then 
PYTHIA can fill missing phase space regions through vetoed showering: let 
the shower sweep over the full phase space, using its PYTHIA-hardness 
ordering, and use the POWHEG-hardness to veto those emissions that POWHEG 
should already have covered. This is only possible since the 
POWHEG-hardness criterion and the shower ordering criterion are very 
similar. In the more general case a truncated showering would be needed 
[<a href="Bibliography.php" target="page">Nas04</a>]. 
 
<p/> 
For vetoed showering, it is necessary to define the POWHEG-hardness 
criterion. In the presence of multiple partons, the definition 
quickly becomes complicated, and allows for different choices. Similar 
decisions have already been made in the implementation of POWHEG, one example 
being the choice in defining which "hardness value" is transferred as 
POWHEG-hardness, e.g. by deciding if the "singular regions" of the FKS or the 
CS approach are used. If the POWHEG-hardness definition were to be changed, 
or extended to more objects, the <code>PowhegHooks.h</code> code would need 
to be modified accordingly. 
 
<p/> 
The merging code is designed to be very flexible, and allows access 
to many possible choices. However, this flexibility means that many parameters 
can be changed, potentially leading to confusion. Thus, recommendations might 
prove helpful. All mistakes and inaccuracies rest with the author. 
 
<p/> 
We recommend the usage of vetoed showers. This means using 
<br/> &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; <code>POWHEG:veto = 1</code> 
<br/> 
This means that PYTHIA will sweep over the full phase space, and apply a veto 
on parton shower emissions for which the POWHEG-hardness separation between 
radiator and emission is above the POWHEG-hardness value of the current input 
event. The variation <code>POWHEG:veto = 0</code> can be used to assess 
how much phase space is under- or double-counted. 
 
<p/> 
To define the POWHEG-hardness criterion, use 
<br/> &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; <code>POWHEG:pTdef = 1</code> 
<br/> 
Other values can be used by experts to assess variations. 
 
<p/> 
Both POWHEG-BOX and PYTHIA 8 generate emissions through a parton shower 
step, meaning that both programs have a clear definition of a radiator 
that emits particles, which is very similar (if not identical). 
To fix the ambiguity if the radiator or the emitted particle should be 
called "the emission", use 
<br/> &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; <code>POWHEG:emitted = 0</code> 
<br/> 
More complicated choices can be used by experts. For instance, use 
<code>POWHEG:emitted = 2</code> to check the POWHEG-hardness of both 
radiator and emitted. 
 
<p/> 
To exhaustively fix the criterion by which to veto parton shower 
emissions, it is important to decide which partons/parton pairs 
are used to calculate the POWHEG hardness of a PYTHIA 8 emission. 
The minimal and recommended choice is 
<br/> &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; <code>POWHEG:pTemt = 0</code> 
<br/> 
This means that only the POWHEG hardness with respect to the radiating leg 
is checked, and recoil effects are neglected. This prescription should be 
very similar to how a hardness value is assigned to a Real-type event 
in the POWHEG-BOX, since in the (implementation of FKS in the) POWHEG-BOX, 
initial state splittings only have singular regions with the radiating 
initial state parton, and final state splittings only have singular 
regions with respect to the radiating final state line. Other choices of 
<code>POWHEG:pTemt</code> are available. A warning is that the impact of 
changes can be huge, particularly for inputs with many jets. Other choices 
therefore should only be made by experts, and a high degree of caution 
is advised. 
 
<p/> 
It is furthermore necessary to decide on a value of the hardness criterion. 
POWHEG-BOX transfers this value in the <code>SCALUP</code> member of 
Les Houches Events, and we recommend using this value by setting 
<br/> &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; <code>POWHEG:pThard = 0</code> 
<br/> 
As a variation, in order to estimate the uncertainty due this choice of 
POWHEG-hardness definition, it can be useful to also check 
<code>POWHEG:pThard = 2</code>. This will recalculate the POWHEG-hardness 
value as promoted in [<a href="Bibliography.php" target="page">Ole12</a>]. 
 
<p/> 
Finally, you need to decide how many emissions the vetoed shower should 
check after an allowed emission has been constructed. If the hardness 
definitions in POWHEG-BOX and PYTHIA 8 where identical, all checking could 
be stopped after the first allowed PS emission. To be prudent, we 
recommend setting 
<br/> &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; <code>POWHEG:vetoCount = 3</code> 
<br/> 
which will then check up to three allowed emissions. Higher values of 
<code>POWHEG:vetoCount</code> have not lead to visible differences 
for the processes which have been tested. 
 
<h3>The modes</h3> 
 
Note that the modes have generally been defined with several default values 
below corresponding to the "off" state, and thus do not agree with the 
recommended values described above. 
 
<br/><br/><table><tr><td><strong>POWHEG:nFinal  </td><td></td><td> <input type="text" name="1" value="2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 1</code>)</td></tr></table>
Number of outgoing particles of POWHEG Born level process, 
i.e. not counting additional POWHEG radiation. 
   
 
<br/><br/><table><tr><td><strong>POWHEG:veto  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 1</code>)</td></tr></table>
Master switch to perform vetoing or not. 
<br/>
<input type="radio" name="2" value="0" checked="checked"><strong>0 </strong>: No vetoing is performed (the user hooks is not loaded).  <br/>
<input type="radio" name="2" value="1"><strong>1 </strong>: Showers are started at the kinematical limit.  Emissions are vetoed if <ei>pTemt > pThard</ei>.  See also <code>POWHEG:vetoCount</code> below.<br/>
 
<br/><br/><table><tr><td><strong>POWHEG:vetoCount  </td><td></td><td> <input type="text" name="3" value="3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>3</strong></code>; <code>minimum = 0</code>)</td></tr></table>
After this many accepted emissions in a row, no more emissions 
are checked. Value 0 means that all emissions are checked. 
   
 
<br/><br/><table><tr><td><strong>POWHEG:pThard  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Selection of the <ei>pThard</ei> scale. For events where there is no 
radiation, <ei>pThard</ei> is always set to be the <code>SCALUP</code> 
value of the LHA/LHEF standard. 
<br/>
<input type="radio" name="4" value="0" checked="checked"><strong>0 </strong>: Set <ei>pThard</ei> equal to <code>SCALUP</code>.<br/>
<input type="radio" name="4" value="1"><strong>1 </strong>: The <ei>pT</ei> of the POWHEG emission is tested against  all other incoming and outgoing partons, with the minimal value chosen.  <br/>
<input type="radio" name="4" value="2"><strong>2 </strong>: The <ei>pT</ei> of all final-state partons is tested  against all other incoming and outgoing partons, with the minimal value  chosen.<br/>
 
<br/><br/><table><tr><td><strong>POWHEG:pTemt  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Selection of the <ei>pTemt</ei> scale. 
<br/>
<input type="radio" name="5" value="0" checked="checked"><strong>0 </strong>: It is the <ei>pT</ei> of the emitted parton with respect  to the radiating parton.<br/>
<input type="radio" name="5" value="1"><strong>1 </strong>: The <ei>pT</ei> of the emission is checked against all  incoming and outgoing partons, and then <ei>pTemt</ei> is set to the  minimum of these values.<br/>
<input type="radio" name="5" value="2"><strong>2 </strong>: The <ei>pT</ei> of all final-state partons is tested  against all other incoming and outgoing partons, with the minimal value  chosen.<br/>
<br/><b>Warning:</b> the choice here can give significant variations 
in the final distributions, notably in the tail to large <ei>pT</ei> values. 
 
<br/><br/><table><tr><td><strong>POWHEG:emitted  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)</td></tr></table>
Selection of emitted parton for FSR. 
<br/>
<input type="radio" name="6" value="0" checked="checked"><strong>0 </strong>: The PYTHIA definition of emitted.<br/>
<input type="radio" name="6" value="1"><strong>1 </strong>: The PYTHIA definition of radiator.<br/>
<input type="radio" name="6" value="2"><strong>2 </strong>: A random selection of emitted or radiator.<br/>
<input type="radio" name="6" value="3"><strong>3 </strong>: Both emitted and radiator are tried.<br/>
 
<br/><br/><table><tr><td><strong>POWHEG:pTdef  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Use of <ei>pT</ei> definitions. 
<br/>
<input type="radio" name="7" value="0" checked="checked"><strong>0 </strong>: The POWHEG ISR <ei>pT</ei> definition for  both ISR and FSR.<br/>
<input type="radio" name="7" value="1"><strong>1 </strong>: The POWHEG ISR <ei>pT</ei> and FSR <ei>d_ij</ei>  definitions.<br/>
<input type="radio" name="7" value="2"><strong>2 </strong>: The PYTHIA definitions.<br/>
 
<br/><br/><table><tr><td><strong>POWHEG:MPIveto  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 1</code>)</td></tr></table>
MPI vetoing. 
<br/>
<input type="radio" name="8" value="0" checked="checked"><strong>0 </strong>: No MPI vetoing is done.<br/>
<input type="radio" name="8" value="1"><strong>1 </strong>: When there is no radiation, MPIs with a scale above  <ei>pT_1</ei> are vetoed, else MPIs with a scale above  <ei>sum_i pT_i / 2 = (pT_1 + pT_2 + pT_3) / 2</ei> are vetoed.  This option is intended specifically for POWHEG simulations of  <ei>2 &rarr; 2 + 2 &rarr; 3</ei> QCD processes.  <br/>
 
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

if($_POST["1"] != "2")
{
$data = "POWHEG:nFinal = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "0")
{
$data = "POWHEG:veto = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "3")
{
$data = "POWHEG:vetoCount = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "0")
{
$data = "POWHEG:pThard = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "0")
{
$data = "POWHEG:pTemt = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "0")
{
$data = "POWHEG:emitted = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "0")
{
$data = "POWHEG:pTdef = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "0")
{
$data = "POWHEG:MPIveto = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
