<html>
<head>
<title>Semi-Internal Resonances</title>
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

<form method='post' action='SemiInternalResonances.php'>
 
<h2>Semi-Internal Resonances</h2> 
 
The introduction of a new <?php $filepath = $_GET["filepath"];
echo "<a href='SemiInternalProcesses.php?filepath=".$filepath."' target='page'>";?> 
semi-internal process</a> may also involve a new particle, 
not currently implemented in PYTHIA. Often it is then enough to 
use the <?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDataScheme.php?filepath=".$filepath."' target='page'>";?>standard machinery</a> 
to introduce a new particle (<code>id:all = ...</code>) and new 
decay channels (<code>id:addChannel = ...</code>). By default this 
only allows you to define a fixed total width and fixed branching 
ratios. Using <code><?php $filepath = $_GET["filepath"];
echo "<a href='ResonanceDecays.php?filepath=".$filepath."' target='page'>";?>meMode</a></code> 
values 100 or bigger provides the possibility of a very 
simple threshold behaviour. 
 
<p/> 
If you want to have complete freedom, however, there are two 
ways to go. One is that you make the resonance decay part of the 
hard process itself, either using the 
<?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>Les Houches interface</a> or 
a semi-internal process. The other is for you to create a new 
<code>ResonanceWidths</code> object, where you write the code 
needed for a calculation of the partial width of a particular 
channel. 
 
<p/> 
Here we will explain what is involved in setting up a resonance. 
Should you actually go ahead with this, it is strongly recommended 
to use an existing resonance as a template, to get the correct 
structure. There also exists a sample main program, 
<code>main22.cc</code>, that illustrates how you could combine 
a new process and a new resonance. 
 
<p/> 
There are three steps involved in implementing a new resonance: 
<br/>1) providing the standard particle information, as already 
outlined above (<code>id:all = ...</code>, 
<code>id:addChannel = ...</code>), except that now branching 
ratios need not be specified, since they anyway will be overwritten 
by the dynamically calculated values. 
<br/>2) writing the class that calculates the partial widths. 
<br/>3) handing in a pointer to an instance of this class to PYTHIA. 
<br/>We consider the latter two aspects in turn. 
 
<h3>The ResonanceWidths Class</h3> 
 
The resonance-width calculation has to be encoded in a new class. 
The relevant code could either be put before the main program in the 
same file, or be stored separately, e.g. in a matched pair 
of <code>.h</code> and <code>.cc</code> files. The latter may be more 
convenient, in particular if the calculations are lengthy, or 
likely to be used in many different runs, but of course requires 
that these additional files are correctly compiled and linked. 
 
<p/> 
The class has to be derived  from the <code>ResonanceWidths</code> 
base class. It can implement a number of methods. The constructor 
and the <code>calcWidth</code> ones are always needed, while others 
are for convenience. Much of the administrative machinery is handled 
by methods in the base class. 
 
<p/>Thus, in particular, you must implement expressions for all 
possible final states, whether switched on in the current run or not, 
since all contribute to the total width needed in the denominator of 
the Breit-Wigner expression. Then the methods in the base class take 
care of selecting only allowed channels where that is required, and 
also of including effects of closed channels in secondary decays. 
These methods can be accessed indirectly via the 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='ResonanceDecays.php?filepath=".$filepath."' target='page'>";?>res...</a></code> 
methods of the normal 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDataScheme.php?filepath=".$filepath."' target='page'>";?>particle database</a></code>. 
 
<p/> 
A <b>constructor</b> for the derived class obviously must be available. 
Here you are quite free to allow a list of arguments, to set 
the parameters of your model. The constructor must call the 
base-class <code>initBasic(idResIn)</code> method, where the argument 
<code>idResIn</code> is the PDG-style identity code you have chosen 
for the new resonance. When you create several related resonances 
as instances of the same class you would naturally make 
<code>idResIn</code> an argument of the constructor; for the 
PYTHIA classes this convention is used also in cases when it is 
not needed. 
<br/>The <code>initBasic(...)</code> method will 
hook up the <code>ResonanceWidths</code> object with the corresponding 
entry in the generic particle database, i.e. with the normal particle 
information you set up in point 1) above. It will store, in base-class 
member variables, a number of quantities that you later may find useful: 
<br/><code>idRes</code> : the identity code you provide; 
<br/><code>hasAntiRes</code> : whether there is an antiparticle; 
<br/><code>mRes</code> : resonance mass; 
<br/><code>GammaRes</code> resonance width; 
<br/><code>m2Res</code> : the squared mass; 
<br/><code>GamMRat</code> : the ratio of width to mass. 
 
<p/> 
A <b>destructor</b> is only needed if you plan to delete the resonance 
before the natural end of the run, and require some special behaviour 
at that point. If you call such a destructor you will leave a pointer 
dangling inside the <code>Pythia</code> object you gave it in to, 
if that still exists. 
 
<a name="method1"></a>
<p/><strong>void ResonanceWidths::initConstants() &nbsp;</strong> <br/>
is called once during initialization, and can then be used to set up 
further parameters specific to this particle species, such as couplings, 
and perform calculations that need not be repeated for each new event, 
thereby saving time. This method needs not be implemented. 
   
 
<a name="method2"></a>
<p/><strong>void ResonanceWidths::calcPreFac(bool calledFromInit = false) &nbsp;</strong> <br/>
is called once a mass has been chosen for the resonance, but before 
a specific final state is considered. This routine can therefore 
be used to perform calculations that otherwise might have to be repeated 
over and over again in <code>calcWidth</code> below. It is optional 
whether you want to use this method, however, or put 
everything in <code>calcWidth()</code>. 
<br/>The optional argument will have the value <code>true</code> when 
the resonance is initialized, and then be <code>false</code> throughout 
the event generation, should you wish to make a distinction. 
In PYTHIA such a distinction is made for <i>gamma^*/Z^0</i> and 
<i>gamma^*/Z^0/Z'^0</i>, owing to the necessity of a special 
description of interference effects, but not for other resonances. 
<br/>In addition to the base-class member variables already described 
above, <code>mHat</code> contains the current mass of the resonance. 
At initialization this agrees with the nominal mass <code>mRes</code>, 
but during the run it will not (in general). 
   
 
<a name="method3"></a>
<p/><strong>void ResonanceWidths::calcWidth(bool calledFromInit = false) &nbsp;</strong> <br/>
is the key method for width calculations and returns a partial width 
value, as further described below. It is called for a specific 
final state, typically in a loop over all allowed final states, 
subsequent to the <code>calcPreFac(...)</code> call above. 
Information on the final state is stored in a number of base-class 
variables, for you to use in your calculations: 
<br/><code>iChannel</code> : the channel number in the list of 
possible decay channels; 
<br/><code>mult</code> : the number of decay products; 
<br/><code>id1, id2, id3</code> : the identity code of up to the first 
three decay products, arranged in descending order of the absolute value 
of the identity code; 
<br/><code>id1Abs, id2Abs, id3Abs</code> : the absolute value of the 
above three identity codes; 
<br/><code>mHat</code> : the current resonance mass, which is the same 
as in the latest <code>calcPreFac(...)</code> call; 
<br/><code>mf1, mf2, mf3</code> : masses of the above decay products; 
<br/><code>mr1, mr2, mr3</code> : squared ratio of the product masses 
to the resonance mass; 
<br/><code>ps</code> : is only meaningful for two-body decays, where it 
gives the phase-space factor 
<i>ps = sqrt( (1. - mr1 - mr2)^2 - 4. * mr1 * mr2 )</i>; 
<br/>In two-body decays the third slot is zero for the above properties. 
Should there be more than three particles in the decay, you would have 
to take care of the subsequent products yourself, e.g. using 
<br/><code>particlePtr->decay[iChannel].product(j);</code> 
<br/>to extract the <code>j</code>'th decay products (with 
<code>j = 0</code> for the first, etc.). Currently we are not aware 
of any such examples. 
<br/>The base class also contains methods for <i>alpha_em</i> and 
<i>alpha_strong</i> evaluation, and can access many standard-model 
couplings; see the existing code for examples. 
<br/>The result of your calculation should be stored in 
<br/><code>widNow</code> : the partial width of the current channel, 
expressed in GeV. 
   
 
<a name="method4"></a>
<p/><strong>double ResonanceWidths::widthChan( double mHat, int idAbs1, int idAbs2) &nbsp;</strong> <br/>
is not normally used. In PYTHIA the only exception is Higgs decays, 
where it is used to define the width (except for colour factors) 
associated with a specific incoming/outgoing state. It allows the 
results of some loop expressions to be pretabulated. 
   
 
<a name="method5"></a>
<p/><strong>bool ResonanceWidths::allowCalc() &nbsp;</strong> <br/>
can normally be left dummy (and then always returns <code>true</code>) but 
can optionally be used to determine whether to force dynamical width 
calculation to be switched off (return <code>false</code>). 
An example is provided by the 
<code>SUSYResonanceWidths</code> class, in which the implementation of 
this method checks for the existence of SLHA decay tables for the 
particular resonance in question, and checks if those tables should be 
given precedence over the internal width calculation. 
   
 
<a name="method6"></a>
<p/><strong>bool ResonanceWidths::initBSM() &nbsp;</strong> <br/>
can normally be left dummy, but for advanced implementations it 
provides a possibility to initialize data members of the derived class 
at a very early stage during initialization, before any of the other 
members are called. An example is 
provided by the <code>SUSYResonanceWidths</code> class, in which 
an internal pointer to a derived <code>Couplings</code> class must be 
(re)set before any of the other methods are used. A return value of 
<code>false</code> can be used to signal that this 
initialization step failed. 
   
 
<h3>Access to resonance widths</h3> 
 
Once you have implemented a class, it is straightforward to 
make use of it in a run. Assume you have written a new class 
<code>MyResonance</code>, which inherits from 
<code>ResonanceWidths</code>. You then create an instance of 
this class and hand it in to a <code>pythia</code> object with 
<pre> 
      ResonanceWidths* myResonance = new MyResonance(); 
      pythia.setResonancePtr( myResonance); 
</pre> 
If you have several resonances you can repeat the procedure any number 
of times. When <code>pythia.init()</code> is called these resonances 
are initialized along with all the internal resonances, and treated in 
exactly the same manner. See also the <?php $filepath = $_GET["filepath"];
echo "<a href='ProgramFlow.php?filepath=".$filepath."' target='page'>";?>Program 
Flow</a> 
description. 
 
<p/> 
If the code should be of good quality and general usefulness, 
it would be simple to include it as a permanently available process 
in the standard program distribution. The final step of that integration 
ought to be left for the PYTHIA authors, but basically all that is 
needed is to add one line in 
<code>ParticleData::initResonances</code>, where one creates an 
instance of the resonance in the same way as for the resonances already 
there. In addition, the particle data and decay table for the new 
resonance has to be added to the permanent 
<?php $filepath = $_GET["filepath"];
echo "<a href='ParticleData.php?filepath=".$filepath."' target='page'>";?>particle database</a>, and the code itself 
to <code>include/ResonanceWidths.h</code> and 
<code>src/ResonanceWidths.cc</code>. 
 
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
