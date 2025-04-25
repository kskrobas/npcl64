/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * help.cpp
 * Copyright (C) 2019 Kazimierz Skrobas <kskrobas@unipress.waw.pl>
 *
 * npcl is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * npcl is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "help.h"
#include "elements.h"
#include <iostream>
#include <iomanip>

using namespace std;
const string tab3("\t\t\t");


inline void print(const string str)
{
cout<<tab3<<str<<endl;
}


void help()
{
cout<<setfill('*')<<setw(80)<<'*'<<endl;

cout<<"usage:"<<endl;
cout<<"npcl <script file> (-v var_name=var_value)*"<<endl;

cout<<" This is a list of acceptable commands/statements:"<<endl<<endl;

cout<<"\t\t1)general statements: "<<endl;
print("avepdh - start if pdh block");
print("break - terminate and next leave the current loop");
print("continue - stop processing of the current loop iteration  and next continue loop with a new iteration");
print("breakIfSNA  - command stops processing the current loop if \'save number of atoms failure (SNA)' flag is set (obsolete)");
print("clearSNA - clear buffer of SNA");
print("cast2int <${math_var}>{1,} - conversion floating point number(s) to integer");
print("diff - start diffraction block");
print("end - the end of any block");
print("end(if|for) - the end of if/for block; it is an optional statement, if not given the 'end' must be used to finish a block");
print("for <iterVar>=<from>(:<step>):<to> - start the 'for' block; range must be given by integer numbers, step is optional, if not given 1 is assumed");
print("for <iterVar>=${var}(:${var}):${var} - as above but range is given by variables. More explanations below");
print("for <iterVar> in <shell listing command> - loop for files returned by a shell listing command");
print("grain - start of grain block");
print("pdh - start of pdh block");
print("print (string|${var}|__progress|\"string\") - print any number of values or strings or loop progress");
print("strRep ${var} <position> <rep. string> - replace a part of string ${var} at position <position> with new value <rep. string>");
print("system command - execute a command by the operating system");
print("threads <value> - global set number of threads for each block");
print("<numericalVariable>=<math expression> - definition of a numerical(double) type variable <var> or expression");
print("<stringVariable>=\"string of characters\" - definition of a string type variable <var>");
print("${<math/string variable>} - call the math/string variable");
print("# or \% - one line comment");
print("/> multiline text </ - respectively begin/end of block/multiline comment");
cout<<endl;


cout<<"\t\t2)grain block: "<<endl;
print("atom <atomName> - atom name of monolattice");
print("atoms <atomName> <atomName> - atom names of binary lattice");
print("center (geom|catom(atom Type)?|id <value>) - center the grain; catom centers to the closest atom; use the specified atom type; center using # of an atom (starts from 0)");
print("charge <elementName> <value> - set the charge of element");
print("comment");
print("csh - start subblock of core-shell paramters");
print("disperse (*|atomName) <value>");
print("disloc (dumbell|intdef) uniform from to prob (x|y|z)? - insert dislocations/defects inside hollow sphere; if dumbell then (x|y|z) selects an insertion plane, else ignored");
print("geometry (sphere|cube|cylinder|cone) - select the shape of a grain");
print("geometry (cubic|oct|dod) <value> - select the shape of a grain based on supersphere for a given p parameter. It works for hcp only.");
print("geometry poly <value>{3} - select supersphere polyhedral shape of grain for parametrs: p, a, b. It works for hcp only.");
print("hcpabc <abcABC|var> - custom sequence of ABC layers. It works only if struct is set to hcp. More explanations below");
print("hcpabc (uniform|normal|lognormal) <value> <value> (ab|abc|abac|abcacb|random) - random number of layers sequences (ab|abc|...) or layers (a|b|c) ");
print("hcpu <value> - distance between two sublayers for zb structures given in percentage");
print("hcpcs (circle|hex|rhomb) - set cross section: circle or hex or rhomb, default: circle");
print("hcpcs poly <value>{3} - set superpolygon shaped cross section based on supercubic shape where 3 parameters control respectively: elongation along x, y axis and third one is a  measure of polyhedrality ");
print("hcpcs poly <value>{5} - set superpolygon shaped cross section where 5 parameters control respectively: 1st- polyhedrality, 2-3 kind of shape, 4-5 elngation along x, y axis");
print("hcpfillup (AB|ABC|ABAC|ABCACB) - hcp filling of perfect hexagonal structure. It works only if struct is set hcp");
print("hcpsl (no|yes)? - enable sublattice for monoatomic grain; hcpsl is always on for biatomic lattices");
print("hcpsurfA (yes|no|AB|H) - surface termination type (obsolete: see hcpsurftype)");
print("hcpsurftyp (<value>|var)  - surface termination type given by a value in range 0...2 corresponding to AA, AB or BB");
print("insfault random <value> - insert stacking faults for hcp at random place(s) with a number given by <value>");
print("lmpmargin <value> - lammps format file, define distance between the side of a box and  most extended atom(s) for xyz axis");
print("lp (<value>|var)  - lattice parameter in Angs");
print("lp (uniform|normal|lognormal) <arg0> <arg1> - lattice parameter based on respective distribution, where:");
print("                               |  arg0  | arg1   |");
print("                       uniform |  from  |   to   |");
print("                       normal  |  mean  |  std   |");
print("                    lognormal  |  mean  |  std   |");
print("mass <atom name> <value> - set the mass of element");
print("open - read atoms from *.xyz. Number of atom types is not limited");
print("printPrm - print current parametrs or status");
print("radius (<value>|var)(lp)? - radius in Angs; if 'lp' is given the multiplicity of lp is taken");
print("radius (uniform|normal|lognormal) <arg0> <arg1> - radius based on respective distribution, see above 'lp'");
print("rename (<oldName> <newName> <prob>){1,} - rename atoms with probability ");
print("rotate <angle> <A> <B> <C> (<xo> <yo> <zo>)? - rotate against ABC axis; angle is given in degrees, axis default position: [0,0,0]");
print("replicate <value>{3} (+/-)?- replicate unit cell along x, y, z axis");
print("side  <value> - size of LAMMPS box side");
print("remove <intvalue - numOfbonds> <value - lenOfbonds> (<value - probability 0..100>)? - remove selected atoms");
print("rescale <value>{1,3}");
print("save <file name> - ASCII file formats: *.ndl, *.xyz, *.dat. Multisaving is possible, just add another entry");
print("saveHeader (fileName|var) - save header only");
print("saveopt (min|max) <value> - save file if number of atoms greater|less than <value>");
print("saveopt (if [HWL][<>][HWL]) - save file if Height/Width/Lenght size is greater|less than other dimension");
print("saveopt lmp tric - save file in lammps format including tilt parameters");
print("struct (sc|bcc|fcc|hcp|zb|zb110|uo2|nife|uc|tric|cif) - select Bravis lattice or zb(110) or hcp or uo2 or uc");
print("tricp - start subblock of unit cell parameters based on triclinic geometry; see below for description");
print("ucp - start subblock of unit cell parameters; see below for description");
print("voids <from> <to> - remove atoms according to uniform distributin given by <from> <to> parametrs");
print("");
print("* The custom hcp sequence is an expresion which consists of number of layers given by a user. An expression has a scheme of numbers and ABC(abc) letters.");
print("* ABC sequences may be nested. The nested sequences are separated by brackets '('  and ')' . Arbitrary nesting is allowed. Expression allows for variables");
print("* Examples: 3(ABC)=> ABCABCABC; 2(A3(BC))=>  ABCBCBCABCBCBC;  AB2(C3(AB))=> AB2(CABABAB)=>ABCABABABCABABAB; var=2 ->  ${var}(abc)=>2(abc)");
cout<<endl;


cout<<"\t\t3)pdh block: "<<endl;
print("bin <value> - bin width. See comment below");
print("comment");
print("difftime - time of computation");
print("mode");
print("open");
print("printPrm");
print("rangeN <start> <number of steps> <stop> - the range of calculation. See comment below");
print("save <fileName> - ASCII file formats: *.pdhs, *.pdhl. ");
print("threads <value> - number of threads engaged during computation");
print("");
print("* the \'bin\' and \'rangeN\' statement must not use within same pdh block");
cout<<endl;


cout<<"\t\t4)diff block: "<<endl;
print("comment <short text> - information (max. 80 chars) added to ASCII *.diff file");
print("difftime - time of computation");
print("dw <u²> - Debyea Waller factor of thermal vibrations (default: 0)");
print("extrapolate (no|yes) - fill up values above 0 to the start of calculation range");
print("fastsinc - fast sinc function; values less than 0.001 are ignored");
print("filesft <file sft> - path to scattering factors file scFact.sft");
print("lambda <value> - radiation wavelength given in Å");
print("lambda (Ag|Cu|Co|Cr|Fe|Mo|W)- radiation wavelength of popular X-Ray anodes given in Å");
print("mode (laue|deb) - select type of calculations, Laue diffraction or  (deb) powder diffraction");
print("noiseSQ");
print("norm - I, S normalization (total intensity equal to 1)");
print("nosf (yes|no) - on/off scaterring factors");
print("polarization (yes|no) - on/off polarization");
print("printPrm - print current parameters");
print("radiation (xray|neutron|electron|nt) - type of radiation; nt - no type, atomic scattering factors are ignored");
print("range <start> <step> <stop>");
print("save <filename> - ASCII file formats: *.dat *.diff. Multisaving is possible, just add another entry");
print("save <filename> ((S|I) (Q|T))? - BIN format: *.diffb with selection of X-Y values to save; if not given all values are stored. See explanation below");
print("threads <value> - number of threads engaged during computation");
print("");
print("* BIN format consists 23 bytes header and data section. Header constists of: version(3 bytes), size of data (8 bytes), saving mode SM(8 bytes) and DATA marker(4 bytes).");
print("* Now only 4 least significant bits are assigned according to scheme: 0 - 2Theta, 1 - Q, 2 - I ,3 - S. The data sections is a set of 8 bytes values");
cout<<endl;


cout<<"\t\t5)gr block: "<<endl;
print("comment - information (max. 80 chars) added to ASCII *.gr file");
print("difftime - time of computation");
print("plot <gnuplot statement> - plot data under gnuplot;  gnuplot statement must not consists of data source, only data formating prms");
print("printPrm - print current parameters");
print("range <start> <step> <stop>");
print("rangeN <start> <number of steps> <stop>");
print("save <filename> - ASCII file formats: *.dat *.gr");
print("save <filename> - BIN file format: *.grb");
print("threads <value> - number of threads engaged during computation. See explanations below");
print("wf (lorch|box) - window function ");
print("");
print("* BIN format consists 15 bytes header and data section. Header constists of: version(3 bytes), size of data (8 bytes) and DATA marker(4 bytes).");
print("* The data sections is a set of 8 bytes values");
print("* the \'range\' and \'rangeN\' statement must not use within same gr block");
cout<<endl;


cout<<"\t\t6)avepdh block: "<<endl;
print("save");
print("open");
cout<<endl;

cout<<"\t\t7)mergepdh block:"<<endl;
print("open ");
print("flimit <N>");
print("range <start> <step> <stop>");
print("weights (uniform|normal|lognormal) <m> <s>");
cout<<endl;


cout<<"\t\t8)ucp subblock:"<<endl;
print("v[xyz] <value>{3} - define translation vectors along x or y or z axis (always 3 vectors must be defined)");
print("<atomName> (<value>|<var>) {3} - define atom name and atom positions of unit cell (at least one atom must be given)");
cout<<endl;

cout<<"\t\t9)tricp subblock:"<<endl;
print("(lpa|lpb|lpc) <value> - define lattice  a, b, c parameters");
print("(alpha|beta|gamma) <value> - define alpha, beta, gamma angles");
print("<atomName> (<value>|<var>){3} - define atom name and atom positions of unit cell (at least one atom must be given)");
cout<<endl;


cout<<"\t\t10)disloc block:"<<endl;
print("mode (rot|cyl|loop) - select type of dislocation");
print("axis <value>{3} - vector perpendicular to plane(s) of dislocated atomic positions");
print("position <value>{3} - position of axis center");
print("rangeR <value>{2,3} - the range of radius values");
print("rangeA <value>{2,3} - the range of rotation angle  values");
print("projh <value> -  cylinder height or max. projection length of vector between axis' center and atom");
print("save <fileName> - save file (same as  grain block 'save' instruction)");
print("scatter xyz  <value>{3} - random scattering of dislocated atoms");
cout<<endl;



cout<<setfill('*')<<setw(80)<<'*'<<endl;

        cout<<"\t\tX-ray anodes and their wavelengths"<<endl;

        for(auto xcat: Elements::xrad)
            cout<<tab3<<xcat.first<<": "<<xcat.second<<endl;


cout<<setfill('*')<<setw(80)<<'*'<<endl;
}
