/*
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
 
# npcl64
The npcl is a program for nanocrystal models building and diffraction calculations. It has command line interface, with simple diagram plotting based on gnuplot. Code is based on standard c++(11 and higher) libraries hovewer was tested for Linux only.


Requirements:
* Function Parser for C++ v4.5.2 by Juha Nieminen, Joel Yliluoma  http://warp.povusers.org/FunctionParser/
* c++ v.11
* gnuplot for plotting (Linux version only)

Installation instructions for Linux like machines:
1. download, install and compile fparser library:

    g++ -c -fpic fparser.cc -o fparser.l
    g++ -shared -o libfparser.so fparser.l 
    
   Set LD_LIBRARY_PATH, e.g.:  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<path_to_fparser>
   If you want add  path permanently, go to the .bashrc file placed in your home directory, insert:
   
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<path_to_fparser> and save. 
   
   Restart your console/terminal.
2. go to npcl64 folder and set link to fparser folder, e.g.: 
    $ ln -s <path_to_fparser> . 
    (remember about dot at then end of instruction)
    Finally the npcl64 folder should look like this:
        
    example
    fparser
    images
    n++
    nf.txt
    README.md
    scFact.sft
    source

    
3. go to the source folder and type:
    make -f Makefile.rel
    
4. ready to use 'npcl' program is placed in the ./npcl64/bin folder
5. Installation: 
    For global installation type:
    $ sudo make -f Makefile.rel sysinstall
    
    For local (to your /home/<user>/bin, the bin folder must exist and should be known for system) folder type:
    $ make -f Makefile.rel install
    
6. set variable NPCLPATH, e.g.:
    $ export NPCLPATH=/home/<user/bin
    
    and copy to it  the scFact.sft   file
    
    To add variable permanently, go to .bashrc (as shown in p.1) and insert export instruction.
    
7. start using 
   
    


Other:
* syntax highlighting for notepad++, see ./n++/npcl.xml

#################################################################


Debyea diffraction mode
![diff](images/im01.png)

##################################################################



Laue diffraction mode
![diffLaue](images/leed_si111_7x7.png)

