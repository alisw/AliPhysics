<html>
<head>
<title>Vertex Information</title>
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

<form method='post' action='VertexInformation.php'>
 
<h2>Vertex Information</h2> 
 
While the setting of secondary production vertices of unstable 
hadrons and leptons is a standard part of the particle decay 
rotines, no corresponding standardized handling is in place for 
the evolution in the partonic or hadronization phases 
of the event generation. The intention is to provide such methods 
in due course. 
 
<p/> 
There are some cases where such information is needed already now, 
specifically for the 
<?php $filepath = $_GET["filepath"];
echo "<a href='RopeHadronization.php?filepath=".$filepath."' target='page'>";?>Rope Hadronization</a> 
framework. Therefore the beginning of a framework is available, 
that can be used to set vertices for partonic production by MPI, 
FSR and ISR. This is done in the <code>PartonVertex</code> class. 
This is a base class, with a default implementation, but the user 
can replace it with a derived class that does a more sophisticated 
handling. 
 
<p/> 
Note that currently the parton-level vertices are expressed in fm, 
unlike the normal mm scale. This will be fixed as the methods 
evolve. Also other improvements and extensions are likely to come. 
So, while people are welcome to write their own derived classes, 
it is likely that these may need to be modified in later PYTHIA 
versions. 
 
<a name="section0"></a> 
<h3>Rope Hadronization Parameters</h3> 
 
Currently the base class implements two alternative approaches to 
picking a partonic vertex, for use inside the rope hadronization 
framework. There are also some free parameters in the models. 
 
<br/><br/><strong>PartonVertex:setVertex</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Master switch to allow the setting of partonic vertices. 
   
 
<br/><br/><table><tr><td><strong>PartonVertex:modeVertex  </td><td>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 1</code>; <code>maximum = 2</code>)</td></tr></table>
Choice between alternative procedures to select vertex. 
<br/>
<input type="radio" name="2" value="1"><strong>1 </strong>: Proton profile is a uniform black disc.  <br/>
<input type="radio" name="2" value="2" checked="checked"><strong>2 </strong>: Proton profile is a two-dimensional Gaussian.  <br/>
 
<br/><br/><table><tr><td><strong>PartonVertex:ProtonRadius </td><td></td><td> <input type="text" name="3" value="0.7" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.7</strong></code>; <code>minimum = 0.</code>; <code>maximum = 10.</code>)</td></tr></table>
The proton radius and shape depends on collision energy. At LHC collision 
energies, say 14 TeV, the profile corresponds roughly to a Gaussian with 
a with of around 0.7 fm, according to the DIPSY model [<a href="Bibliography.php#refFle11" target="page">Fle11</a>]. 
   
 
<br/><br/><table><tr><td><strong>PartonVertex:EmissionWidth </td><td></td><td> <input type="text" name="4" value="0.1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.1</strong></code>; <code>minimum = 0.</code>; <code>maximum = 10.</code>)</td></tr></table>
Vertices of ISR+FSR partons are smeared relative to their mother by a 
Gaussian distribution with a width of <code>EmissionWidth</code>/<i>pT</i>, 
where <i>pT</i> is the transverse momentum of the emission (in GeV). 
This parameter thus determined the overall strength of the transverse space 
smearing. 
   
 
<br/><br/><table><tr><td><strong>PartonVertex:pTmin </td><td></td><td> <input type="text" name="5" value="0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2</strong></code>; <code>minimum = 0.05</code>; <code>maximum = 1.</code>)</td></tr></table>
The parton transverse smearing is assumed proportional to <i>1 / pT</i>, 
but <i>pT</i> is limited to be above this parameter so as to avoid 
unreasonable values. 
   
 
<a name="section1"></a> 
<h3>External models</h3> 
 
A derived class for setting parton vertex information can be provided 
to PYTHIA with the 
<br/><code> 
bool Pythia::setPartonVertexPtr( PartonVertex* partonVertexPtrIn) 
</code><br/> 
method. The methods in the derived <code>PartonVertex</code> class 
can then be used to add vertex information to produced particles, 
at creation time, in MPI, FSR and ISR. The assigned vertex information 
will afterwards be accessible as properties of the individual particles. 
Particles produced in other types of processes than the ones mentioned 
above will not have vertex information assigned (e.g. hard process, 
beam remnants etc.), neither will particles produced in the weak shower. 
 
<a name="anchor1"></a>
<p/><strong> virtual void init() &nbsp;</strong> <br/>
can be used to initialize various parameters of the model or precalculate 
common numbers. Note that a separate non-virtual method will already 
have provided pointers to the <code>Info</code>, <code>Settings</code> 
and <code>Rndm</code> classes, so that these are available in all derived 
classes. 
   
 
<a name="anchor2"></a>
<p/><strong> virtual void vertexMPI( int iBeg, int nAdd, double bNow, Event& event) &nbsp;</strong> <br/>
Method to assign a production vertex to a particle produced in the MPI 
framework. Should set the vertices <code>vProd</code> of the particles 
concerned. 
<br/><code>argument</code><strong> iBeg </strong>  :  is the index of the first parton of a MPI. 
   
<br/><code>argument</code><strong> nAdd </strong>  :  is the number of partons involved in the MPI, 
currently always four: two in and two out. 
   
<br/><code>argument</code><strong> bNow </strong>  :  is the impact parameter of the event. It is not 
expressed in physical units (like fm), but rescaled such that the average 
is unity for MPI events. See the section on 
<?php $filepath = $_GET["filepath"];
echo "<a href='MultipartonInteractions.php?filepath=".$filepath."' target='page'>";?>Multiparton Interactions</a> for 
a description of choices for the <i>b</i> dependence. 
   
<br/><code>argument</code><strong> event </strong>  :  reference to the whole event, to read information 
from and set the relevant <code>vProd</code> values into. 
   
   
 
<a name="anchor3"></a>
<p/><strong> virtual Vec4 vertexFSR( int iNow, Event& event) &nbsp;</strong> <br/>
Method to assign production vertex to a particle produced in the FSR 
(<code>TimeShower</code>). Should set the vertex <code>vProd</code> 
of the particle concerned. 
<br/><code>argument</code><strong> iNow </strong>  :  is the index of the parton concerned. In a 
branching the daughters automatically inherit the vertex of the mother, 
if it has one, and similarly for the recoiler. This method is called 
specifically for what is considered the emitted parton of the process, 
i.e. the gluon in a <i>q &rarr; q g</i> branching, and allows the 
vertex of this parton to be modified. 
   
<br/><code>argument</code><strong> event </strong>  :  reference to the whole event, to read information 
from and set the relevant <code>vProd</code> values into. 
   
   
 
<a name="anchor4"></a>
<p/><strong> virtual Vec4 vertexISR( int iNow, Event& event) &nbsp;</strong> <br/>
Method to assign production vertex to a particle produced in the ISR 
(<code>SpaceShower</code>). Should set the vertices <code>vProd</code> 
of the particle concerned. 
<br/><code>argument</code><strong> iNow </strong>  :  is the index of the parton concerned. This method 
is called three times for each ISR branching, for the daughter, the 
new recoiler and the sister. 
   
<br/><code>argument</code><strong> event </strong>  :  reference to the whole event, to read information 
from and set the relevant <code>vProd</code> values into. 
   
   
 
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

if($_POST["1"] != "off")
{
$data = "PartonVertex:setVertex = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "2")
{
$data = "PartonVertex:modeVertex = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "0.7")
{
$data = "PartonVertex:ProtonRadius = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "0.1")
{
$data = "PartonVertex:EmissionWidth = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "0.2")
{
$data = "PartonVertex:pTmin = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
