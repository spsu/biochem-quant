/**
 * This small program calculates the Isoelectric Point of proteins.
 * Obviously not as thorough as actual computational chemistry 
 * algorithms, so don't use this for anything meaningful.
 *
 * Copyright 2011 Brandon Thomas
 * Available under BSD License.
 */

// TODO: Don't install everything in global namespace. 

// Amino Acid factory
var AA = function(name, carboxy, amino, side, charge)
{
	var AminoAcid = function() {
		this.name = name;
		this.carboxy = carboxy;
		this.amino = amino;
		this.side = side || null;
		this.sideCharge = charge || 0;
	}

	return new AminoAcid;
}

// These are the pKa values for Amino Acid functional groups.
// Indexed as: Carboxyl, Amino[, SideChain, StartCharge]
var PKA_VALUES = {
	G: AA('Glycine',		2.4, 9.8),
	A: AA('Alanine',		2.4, 9.9),
	V: AA('Valine',			2.3, 9.7),
	L: AA('Leucine',		2.3, 9.7),
	I: AA('Isoleucine',		2.3, 9.8),
	M: AA('Methionine',		2.1, 9.3),
	P: AA('Proline',		2.0, 10.6),
	F: AA('Phenylalanine',	2.2, 9.3),
	W: AA('Tryptophan',		2.5, 9.4),
	S: AA('Serine',			2.2, 9.2),
	T: AA('Threonine',		2.1, 9.1),
	C: AA('Cysteine',		1.9, 10.7, 8.4, 0),
	Y: AA('Tyrosine',		2.2, 9.2, 10.5, 0),
	N: AA('Asparagine',		2.1, 8.7),
	Q: AA('Glutamine',		2.2, 9.1),
	D: AA('Aspartate',		2.0, 9.9, 3.9, 0),
	E: AA('Glutamate',		2.1, 9.5, 4.1, 0),
	K: AA('Lysine',			2.2, 9.1, 10.5, 1),
	R: AA('Arginine',		1.8, 9.0, 12.5, 1),
	H: AA('Histidine',		1.8, 9.3, 6.0, 1)
}

// Group factory
// Type - 0: carboxy, 1: amino, 2: side
var get_group = function(code, type)
{
	var aa, pka, charge;

	// Group class
	var Group = function(aa, type, pka, charge)
	{
		this.aa = aa;
		this.type = type;
		this.pka = pka;
		this.charge = charge;

		// Function of pH
		this.chargeAt = function(ph) {
			if(ph >= pka) {
				return charge - 1;
			}
			return charge;
		}

		this.getType = function() {
			switch(this.type) {
				case 0:
					return 'c-term';
				case 1:
					return 'n-term';
				case 2:
					return 'side chain';
			}
		}
		// formatted
		this.getName = function() {
			return this.aa.name + ' ' + this.getType();
		}
	}

	if(!code || typeof(type) == 'undefined') {
		return false;
	}

	if(!(code in PKA_VALUES)) {
		return false;
	}

	aa = PKA_VALUES[code];

	switch(type) {
		case 0:
			// carboxy group
			pka = aa.carboxy;
			charge = 0;
			break;
		case 1: 
			// amino group
			pka = aa.amino;
			charge = 1;
			break;
		case 2:
			// R-side chain
			if(!aa.side) {
				return false;
			}
			pka = aa.side;
			charge = aa.sideCharge;
			break;
	}

	return new Group(aa, type, pka, charge);
}

var pi_main = function()
{
	// numerically sort groups array
	var sort_groups = function(arr) {
		arr.sort(function(a, b) { return a.pka - b.pka; });
	}

	// Grab and filter text input
	var grab_input = function() {
		return $("#protein").val().replace(/[\W\d]/g, '').toUpperCase();
	}

	// Is it a valid amino acid letter code?
	var is_valid_code = function(code) {
		var ch;
		if(code.length == 0) {
			return false;
		}
		for(var i=0; i < code.length; i++) {
			ch = code[i];
			if(!(ch in PKA_VALUES)) {
				return false;
			}
		}
		return true;
	}

	// Get a list of all ionizable groups in the polypeptide
	var get_groups = function(code) {
		var c, n, r, pkas, i;
		n = get_group(code[0], 1);
		c = get_group(code[code.length-1], 0);

		// list of ionizable groups
		pkas = [n, c];

		for(i = 0; i < code.length; i++) {
			r = get_group(code[i], 2);
			if(!r) {
				continue;
			}
			pkas.push(r);
		}

		//return sorted_copy_pkas(pkas);
		return pkas;
	}

	var list_groups = function(groups) {
		var g;
		var str = '<h3>Ionizable Groups</h1><ul>\n';
		for(var i = 0; i < groups.length; i++) {
			g = groups[i];
			str += '<li>' + g.getName() + ' ';
			str += '<em>pKa ' + g.pka + '</em></li>\n';
		}
		str += '</ul>';
		$('#groups').html(str);
	}

	var calculate_pi = function(groups) {
		sort_groups(groups);
		var steps = [];

		// Initial state (prior to any deprotonation)
		var ch = 0;
		for(var i = 0; i < groups.length; i++) {
			ch += groups[i].chargeAt(-10);
		}
		steps.push({ph: 'min', charge:ch, groups:[]});

		// Titrate from 0 to maximum ph; stepwise deprotonation.
		var cur = 0;
		while(cur < groups.length) {
			var curPh = groups[cur].pka;	
			var numIonized = 0;
			for(var i = 0; i < groups.length; i++) {
				if(groups[i].pka == curPh) {
					numIonized++;
				}
				else if(groups[i].pka > curPh) {
					break;
				}
			}
			var step = {ph: curPh, charge:0, groups:[]};
			for(var i = 0; i < numIonized; i++) {
				step.groups.push(groups[cur + i]);
			}
			for(var i = 0; i < groups.length; i++) {
				var ch = groups[i].chargeAt(curPh);
				step.charge += ch;
			}
			steps.push(step);
			cur += numIonized;
		}

		// Final state (maximum pH).
		/*ch = 0;
		for(var i = 0; i < groups.length; i++) {
			ch += groups[i].chargeAt(100);
		}
		steps.push({ph:'max', charge:ch, groups:[]});*/

		// Calculate pI -- two closest flanking charges to zero.
		var pos = 10000;
		var neg = -10000;
		var posPh = null;
		var negPh = null;
		for(var i = 0; i < steps.length; i++) {
			var step = steps[i];
			if(step.charge > 0 && step.charge < pos) {
				pos = step.charge;
				posPh = steps[i+1].ph;
			}
			if(step.charge < 0 && step.charge > neg) {
				neg = step.charge;
				negPh = steps[i].ph;
			}
		
		}
		var pI = (negPh + posPh)/2;

		var html = '';
		for(var i = 0; i < steps.length; i++) {
			var step = steps[i];
			var title = '<h3 class="charge' + step.charge +'">';
			title += 'State ' + (i+1) + ' ';
			title += '<small>(charge = ' + step.charge;
			title += ', ph = ' + step.ph + ')</small>'+ '</h3>';
			var li = '';
			if(step.groups.length == 0) {
				li = '<ul><li>Nothing deprotonated.</li></ul>';
			}
			else {
				var li = '<ul>\n';
				for(var j = 0; j < step.groups.length; j++) {
					g = step.groups[j];
					li += '<li>' + g.getName() + ' ';
					li += '<em>(pKa ' + g.pka + ')</em> deprotonated</li>\n';
				}
				li += '</ul>';
			}

			html += title + li;
		}

		html += "<h2>pI = ";
		html += "(" + negPh + " + ";
		html += posPh + ")/2";
		html += " = " + pI + "</h2>";
		$('#groups').append(html);
	}

	var process = function() {
		var code = grab_input();
		if(code.length == 0) {
			$('#error').empty();
			$('#groups').empty();
			return false;
		}
		if(!is_valid_code(code)) {
			$('#error').html('<p>Invalid AA code.</p>');
			return false;
		}
		else {
			$('#error').empty();
		}

		var groups = get_groups(code);

		sort_groups(groups);
		list_groups(groups);
		calculate_pi(groups);
	}
	// Install callback
	//$("#protein").change(function(){ process(); });
	$("#protein").keyup(function(){ process(); });
}

