<html>
<head>
<title>Beam Shape</title>
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

<form method='post' action='BeamShape.php'>
 
<h2>Beam Shape</h2> 
 
The <?php $filepath = $_GET["filepath"];
echo "<a href='BeamParameters.php?filepath=".$filepath."' target='page'>";?>Beam Parameters</a> 
page explains how you can set a momentum spread of the two incoming 
beams, and a spread and offset for the location of the interaction 
vertex. The spread is based on a simple parametrization in terms of 
independent Gaussians, however, which is likely to be too primitive 
for realistic applications. 
 
<p/> 
It is therefore possible to define your own class, derived from the 
<code>BeamShape</code> base class, and hand it in to Pythia with the 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='ProgramFlow.php?filepath=".$filepath."' target='page'>";?> 
pythia.setBeamShapePtr( BeamShape*)</a></code> method. 
Below we describe what such a class has to do. An explicit toy example 
is shown in <code>main23.cc</code>. 
 
<p/> 
The <code>BeamShape</code> base class has a very simple structure. 
It only has two main virtual methods. The first, <code>init()</code>, is 
used for initialization. The second, <code>pick()</code>, selects 
beam momentum and production vertex in the current event. 
 
<a name="method1"></a>
<p/><strong>BeamShape::BeamShape() &nbsp;</strong> <br/>
   
<strong>virtual BeamShape::~BeamShape() &nbsp;</strong> <br/>
the constructor and destructor do not need to do anything. 
   
 
<a name="method2"></a>
<p/><strong>virtual void BeamShape::init( Settings& settings, Rndm* rndmPtrIn) &nbsp;</strong> <br/>
the base-class method simply reads in the relevant values stored 
in the <code>Settings</code> data base, and saves a pointer to the 
random-number generator. You are free to write your own 
derived initialization routine, or use the existing one. In the 
latter case you can then give your own modified interpretation 
to the beam spread parameters defined there. 
<br/>The two flags <code>Beams:allowMomentumSpread</code> and 
<code>Beams:allowVertexSpread</code> should not be tampered with, 
however. These are checked elsewhere to determine whether the beam 
shape should be set or not, whereas the other momentum-spread 
and vertex-spread parameters are local to this class. 
   
 
<a name="method3"></a>
<p/><strong>virtual void BeamShape::pick() &nbsp;</strong> <br/>
this method is the key one to supply in the derived class. Here you 
are free to pick whatever parametrization you desire for beam momenta 
and vertex position, including correlations between the two. 
At the end of the day, you should set a few protected 
<code>double</code> numbers: 
<br/><code>deltaPxA, deltaPyA, deltaPzA</code> for the three-momentum 
shift of the first incoming beam, relative to the nominal values; 
<br/><code>deltaPxB, deltaPyB, deltaPzB</code> for the three-momentum 
shift of the second incoming beam, relative to the nominal values; 
<br/><code>vertexX, vertexY, vertexZ, vertexT</code> for the 
production-vertex position and time. 
<br/>As usual, momentum is given in GeV, and space and time in mm, 
with <i>c = 1</i>. 
   
 
<a name="method4"></a>
<p/><strong>Vec4 BeamShape::deltaPA() &nbsp;</strong> <br/>
   
<strong>Vec4 BeamShape::deltaPB() &nbsp;</strong> <br/>
read out the three-momentum shifts for beams A and B that were set by 
<code>pick()</code>. The energy components are put to zero at this stage, 
since they are most conveniently calculated after the original and the 
shift three-momenta have been added. 
   
 
<a name="method5"></a>
<p/><strong>Vec4 BeamShape::vertex() &nbsp;</strong> <br/>
read out the production-vertex position and time that were set by 
<code>pick()</code>. 
   
 
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
