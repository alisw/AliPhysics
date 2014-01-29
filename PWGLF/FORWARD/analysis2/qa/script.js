function toggleDisplay(node,name,forceHide)
{
    if (node.style.display == "none" && forceHide == false) {
	node.style.display = "table-row";
	return true;
    }
    else {
	node.style.display = "none";
	toggle(name, true);
    }
    return false;
}
function toggleVis(node,name,forceHide)
{
    if (node.style.visibility == "collapse" && forceHide == false) {
	node.style.visibility = "visible";
	return true;
    }
    else {
	node.style.visibility = "collapse";
	toggle(name, true);
    }
    return false;
}

function hideAll()
{
    var i = 1;
    while (true) { 
	var nam = 'r' + i;
	var row = document.getElementById(nam);
	if (row == null) break;
	// console.log('Closing ' + nam);
	toggle(nam, true);
	i++;
    }
}
function expandAll()
{
    hideAll();
    var i = 1;
    while (true) { 
	var nam = 'r' + i;
	var row = document.getElementById(nam);
	if (row == null) break;
	// console.log('Closing ' + nam);
	toggle(nam, false);
	i++;
    }
}
    
function toggle(node,forceHide)
{
    var i = 1;
    var c = 0;
    // console.log('toggle ' + node);
    while (true) {
	var sub = node + "." + i;
	var row = document.getElementById(sub);
	// console.log('got element ' + sub + ': ' + row);
	if (row) {
	    if (toggleDisplay(row, sub, forceHide)) c++;
	}
	else 
	    break;
	i++;
    } // while 
    var arn = node.replace('r','a');
    var arr = document.getElementById(arn);
    // console.log('Got arrow ' + arn + ': ' + arr + ' c=' + c);
    if (arr) {
	// console.log("inner HTML of " + arn + ': ' + arr.innerHTML);
	if (arr.innerHTML != "&nbsp;") {
	    if (c > 0) arr.innerHTML = "&blacktriangledown;";
	    else       arr.innerHTML = "&blacktriangleright;"
	}
    }   
}
