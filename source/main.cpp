/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * main.cpp
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

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <thread>
#include <future>

#include <omp.h>
#include <stdlib.h>
//#include <stdio.h>

#include <cstdio>
#include <unordered_map>

//#ifndef __linux
//#include <conio.h>
//#endif


#ifdef DEBUG
#define DB true
#else
#define DB false
#endif

using namespace std;

#include "scriptanalyser.h"
#include "cgr.h"
#include "cavepdh.h"
#include "cmergepdh.h"
#include "../fparser/fparser.hh"
#include "help.h"
#include "crandom.h"
#include "colormsg.h"


typedef std::string str;




double psqr (const double *p){ return p[0]*p[0]; }
double prand(const double *p){ CrandomUni runi(p[0],p[1]); return runi.randNumber(); }

//=============================================================================

class ClMainTask{
private:
        vcmdlist cmdlist;
        stdumap  uvars; //user variables

        NanoGrain::StNanoGrain grain;
        Cpdh pdh;
        Cdiff diff;
        Cgr gr;
        Cavepdh avepdh;
        Cmergepdh mergepdh;

        int argc,options;
        size_t cmdlistSize;
        std::vector<string> argv;

        void mthelp()
        {

            help();

            cerr<<" date: "<<__DATE__<<endl;
            cerr<<" author: Kazimierz Skrobas, kskrobas@unipress.waw.pl"<<endl;
            cerr<<" NPCLPATH: ";

                    if(const char* env_p = std::getenv("NPCLPATH"))
                        cerr<<env_p;

                    cerr<<endl;

            std::unordered_map<unsigned,std::string> map{{200505,"2.5"},{200805,"3.0"},{201107,"3.1"},{201307,"4.0"},{201511,"4.5"}};

         //   printf(" OpenMP version: %s \n", map.at(_OPENMP).c_str());
        }

public:

        ClMainTask()
        {
            pdh.grain=&grain;
            avepdh.pdh=&pdh;
            diff.pdh=&pdh;
            gr.diff=&diff;
            mergepdh.pdh=&pdh;
        }

        int retvalue;

        //=====================================================================
        bool argProcessing(const int &argc__,char *argv__[])
        {
            if(argc__<2){
                mthelp();
                retvalue=-1;
            return false;
            }



            argc=argc__;
            argv.resize(argc);

            for(int i=0;i<argc;i++){
            str sargv(argv__[i]);
                argv[i]=std::move(sargv);

                if(argv[i]=="-v"){
                    if( (i+1) <argc){
                    str sargvP(argv__[++i]);
                    vector<string> varname(split<string>(sargvP,"="));

                            if(varname.size()==1) {cerr<<"value of '"<<varname[0]<<"'' is not given "<<endl; return false;}

                    ClKeyValues kv;
                                    kv<<"$var_lit"<<varname[0]<<varname[1];
                                    cmdlist.emplace_back(kv);
                    }
                    else
                       { cerr<<"wrong number of parameters"<<endl; return false;}
                }
            }

            options=0;           

        return true;
        }

        //=====================================================================
        void operator() ()
        {
                parsingProc();
        }

        //=====================================================================

        void parsingProc()
        {

           // cout<<"script parsing "<<__LINE__<<endl;

         fstream script(argv[1],ios::in);

                if(!script){
                    cerr<<" couldn't open the script file"<<endl;
                return ;
                }

         size_t cline=0;

                if( auto retVal=Script::scriptParsing(script,cline,&cmdlist,&uvars,options)==Script::Result::OK){
                    if(cmdlist.empty())
                        cout<<" WARNING, empty script"<<endl;
                }
                else{
                    if(retVal==Script::Result::ENDRET)
                        errMsg("misplaced 'end' , line: "+std::to_string(cline));
                    else
                        if(retVal==Script::Result::ENDELSE)
                            errMsg("misplaced 'else' , line: "+std::to_string(cline));
                        else
                            errMsg("unknown line: "+std::to_string(cline));

                    retvalue=-2;
                    script.close();
                    return;
                }

                script.close();
                uvars.clear();

         size_t cmdIndex=0;

                cmdlistSize=cmdlist.size();

                doCommands(cmdIndex);

        }
        //=====================================================================

        bool doCommands(size_t &cmdIndex)
        {
        size_t op;
        const size_t cmdSize=cmdlist.size();
        //const size_t cmdIndexBegin=cmdIndex;
                    	

                    try{

                        for(op=1;cmdIndex<cmdSize;cmdIndex++,op++){

                            if(DB) cout<<":::::::: "<<cmdlist[cmdIndex]<<endl;

                            /////////////
                            if(cmdlist[cmdIndex]=="end" || cmdlist[cmdIndex]=="endfor" )
                            break;

                            if(cmdlist[cmdIndex]=="else"){
                            break;
                            }

                            if( cmdlist[cmdIndex]=="break"){
                                if(grain.testSNA){
                                    doIgnoreCommands(cmdIndex);
                                break;
                                }
                                else
                                continue;
                            }


                            ////////////
                            if(cmdlist[cmdIndex]=="$var_math"){
                            string vartoeval(cmdlist[cmdIndex][2]);
                            FunctionParser fparser;

                                    Script::replaceVars(&uvars,vartoeval);
                                    fparser.AddFunction("sqr",psqr,1);
                                    fparser.AddFunction("rand",prand,2);


                            int res = fparser.Parse(vartoeval, "");

                                    if(res >= 0) {
                                        cerr<<"ERROR: illegal math expression "<<cmdlist[cmdIndex][2]<<endl;
                                        throw Script::ERR_MATH_EXP;
                                    }

                            double val;
                            const double veval=fparser.Eval(&val);

                            auto iterVar=std::find(uvars.begin(),uvars.end(),cmdlist[cmdIndex][1]);

                                        if(iterVar==uvars.end())
                                            uvars.emplace_back(strpair(cmdlist[cmdIndex][1],std::to_string(veval)));
                                        else                                                                                    
                                            iterVar->getValue()=std::to_string(veval);

                            continue;
                            }

                            ////////////
                            if(cmdlist[cmdIndex]=="$var_lit"){
                            string vartoeval(cmdlist[cmdIndex][2]);
                                Script::replaceVars(&uvars,vartoeval);

                            auto iterVar=std::find(uvars.begin(),uvars.end(),cmdlist[cmdIndex][1]);

                                    if(iterVar==uvars.end())
                                        uvars.emplace_back(strpair(cmdlist[cmdIndex][1],vartoeval));
                                    else
                                        iterVar->getValue()=vartoeval;

                            continue;
                            }
                            ////////////
                            if(cmdlist[cmdIndex]=="incdec"){
                            const string var(cmdlist[cmdIndex][1]);
                            auto iterVar=std::find(uvars.begin(),uvars.end(),var);

                                    if(iterVar==uvars.end())
                                        throw Script::ERR_UNK_VAR;
                                    else{
                                    const string oper(cmdlist[cmdIndex][2]);
                                    const string sval(iterVar->getValue());
                                    double val=std::stod(sval);

                                        if(oper=="++") val++;
                                        else val--;

                                         iterVar->getValue()=std::to_string(val);
                                    }

                            continue;
                            }

                            ///
                            if(cmdlist[cmdIndex]=="grain"){
                                cout<<"grain building/reading"<<endl;
                                
                                if(!grain.parseCommands(cmdlist,++cmdIndex,&uvars))
                                    break;
                                    				
                            continue;
                            }

                            ////
                            if(cmdlist[cmdIndex]=="pdh"){
                                cout<<"pdh"<<endl;
                                
                                if(!pdh.parseCommands(cmdlist,++cmdIndex,&uvars))
                                    break;
                                                                       
                                pdh.calc();
                                cout<<"\n";

                                if(pdh.status!=Cpdh::OK){
                                    cerr<<"ERROR: PDH status not OK"<<endl;
                                throw 0;
                                }
                                
                            continue;
                            }

                            ////
                            if(cmdlist[cmdIndex]=="diff"){
                                cout<<"diffraction"<<endl;

                                if(!diff.parseCommands(cmdlist,++cmdIndex,&uvars))
                                    break;

                                diff.calc();
                                cout<<"\n";


                            continue;
                            }

                            ////
                            if(cmdlist[cmdIndex]=="gr"){
                                cout<<"G(r)"<<endl;

                                if(!gr.parseCommands(cmdlist,++cmdIndex,&uvars))
                                    break;

                                gr.calc();
                                cout<<"\n";

                            continue;
                            }

                            ////
                            if(cmdlist[cmdIndex]=="avepdh"){
                                cout<<"avepdh"<<endl;

                                if(!avepdh.parseCommands(cmdlist,++cmdIndex,&uvars))
                                    break;

                                avepdh.calc();
                                cout<<"\n";

                                if(avepdh.status!=Cavepdh::OK)
                                    throw 0;
                            continue;
                            }
                            ////
                            if(cmdlist[cmdIndex]=="mergepdh"){
                                cout<<"mergepdh"<<endl;

                                if(!mergepdh.parseCommands(cmdlist,++cmdIndex,&uvars))
                                    break;

                                mergepdh.calc();
                                cout<<"\n";

                                if(mergepdh.status!=Cmergepdh::OK)
                                    throw 0;
                            continue;
                            }


                            ////
                            if(cmdlist[cmdIndex]=="system"){                            	
                            str scmd;
                                                        		                            
                                for(int i=1;i<cmdlist[cmdIndex].numOfKeyValues();i++)
                            		scmd+=cmdlist[cmdIndex][i]+" "; 

                                Script::replaceVars(&uvars,scmd);
                                system (scmd.c_str());
                            continue;	
							}


                            if(cmdlist[cmdIndex]=="for"){
                                if(!doForCommand(cmdIndex))
                                    throw 0;
                            continue;
                            }

                            #ifdef __linux__
                            if(cmdlist[cmdIndex]=="for_in"){
                                if(!doForInCommand(cmdIndex))
                                    throw 0;
                            continue;
                            }
                            #endif

                            if(cmdlist[cmdIndex]=="if"){
                                if(!doIfCommand(cmdIndex))
                                    throw 0;
                            continue;
                            }



                            if(cmdlist[cmdIndex]=="cast2int"){
                            string varTofind(cmdlist[cmdIndex][1]);

                                    /// wyluskuje rdzen nazwy, czyli wartosc pomiedzy {...}
                                    varTofind.erase(varTofind.begin(),varTofind.begin()+2);
                                    varTofind.erase(varTofind.end()-1,varTofind.end());

                            auto it=std::find(uvars.begin(),uvars.end(),varTofind);

                                    if(it!=uvars.end()){
                                    const int ival=std::stoi(it->getValue());
                                        it->getValue()=std::to_string(ival);
                                    }
                                    else{
                                            cout<<" ERROR, unknown keyvalue "<<cmdlist[cmdIndex][1];
                                            throw Script::ERR_UNK_VAR;
                                    }

                            continue;
                            }

                            if(cmdlist[cmdIndex]=="clearSNA"){
                                NanoGrain::StNanoGrain::savedNumOfAtoms.clear();
                            continue;
                            }


                            if(cmdlist[cmdIndex]=="print"){                                

                                for(int i=1;i<cmdlist[cmdIndex].numOfKeyValues();i++){
                                cout<<" ";

                                    if(cmdlist[cmdIndex][i][0]=='$'){
                                    string varTofind(cmdlist[cmdIndex][i]);

                                            /// wyluskuje rdzen nazwy, czyli wartosc pomierdzy {...}
                                            varTofind.erase(varTofind.begin(),varTofind.begin()+2);
                                            varTofind.erase(varTofind.end()-1,varTofind.end());

                                    auto it=std::find(uvars.begin(),uvars.end(),varTofind);

                                        if(it!=uvars.end())
                                            cout<<it->getValue();
                                        else{
                                            cout<<" ERROR, unknown keyvalue "<<cmdlist[cmdIndex][i];
                                            throw Script::ERR_UNK_VAR;
                                        }
                                    }
                                    else{
                                        cout<<cmdlist[cmdIndex][i];
                                    }
                                }

                                cout<<endl;
                            continue;
                            }


                            if(cmdlist[cmdIndex]=="printvars"){
                                cout<<"variables :"<<endl;

                                for(auto & var: uvars)
                                    cout<<var<<endl;
                            continue;
                            }

                            cerr<<__FILE__<<":"<<__LINE__<<" unknown command "<<cmdlist[cmdIndex]<<endl;
                            throw 0;
                        }
                    }
                    catch(Script::ResultRepVar e){
                        cerr<<"ERROR: unknown variable, code: "<<e<<endl;
                    return false;
                    }
                    catch(Script::Result sr){
                        switch (sr){
                        case Script::Result::ERR_UNK_VAR:  cerr<<" unknown variable "<<endl; break;
                        default : cerr<<"ERROR: script error, code: "<<sr<<endl;
                        }
                    return false;
                    }
                    catch(...){
                        cerr<<"ERROR: unknown code"<<endl;
                    return false;
                    }
        return true;
        }

        //=====================================================================

        bool doForCommand(size_t &cmdIndex)
        {        
        const string cmd(cmdlist[cmdIndex][1]);
        vector<string> forTokens (split<string>( cmd          ,"="));
        vector<string> iterRange (split<string>( forTokens[1] ,":"));
        const size_t currPos=uvars.size();

                    Script::replaceVars(&uvars,iterRange[0]);
                    Script::replaceVars(&uvars,iterRange[1]);
                    if(iterRange.size()==3) Script::replaceVars(&uvars,iterRange[2]);

        const int iBeg=std::stoi(iterRange[0]);
        const int iter=(iterRange.size()==2) ? 1 : std::stoi(iterRange[1]);
        const int iEnd=(iterRange.size()==2) ? std::stoi(iterRange[1]) : std::stoi(iterRange[2]);

                    if(iter==0)             {cerr<<"ERROR: iteration step must not be 0"<<endl;            throw Script::ERR_INV_PAR; }
                    if(iEnd>iBeg && iter<0) {cerr<<"ERROR: invalid 'for' parameters  "<<iBeg<<":"<<iEnd<<endl;throw Script::ERR_INF_LOOP; }
                    if(iEnd<iBeg && iter>0) {cerr<<"ERROR: invalid 'for' parameters   "<<iBeg<<":"<<iEnd<<endl;throw Script::ERR_INF_LOOP; }

                    uvars.emplace_back(strpair(forTokens[0],iterRange[0]));

        size_t forCmdIndex=0;

                for(int it=iBeg;it<=iEnd;it+=iter){
                    auto locIterVar=std::find(uvars.begin(),uvars.end(),forTokens[0]);
                    locIterVar->getValue()=std::to_string(it);

                    forCmdIndex=cmdIndex+1; // always reset loop to the first statement
                    if(!doCommands(forCmdIndex))
                        return false;

                }

                /// probably below 'if' lines should be removed (they don't look good)
                /// after update 03-04-2023  searching for "end" instruction
                /// is done by doIgnoreCommand placed near line 230 ( doCommands-> ... if ("break") ....)
                if(forCmdIndex<cmdlist.size()){
                    if(cmdlist[forCmdIndex]=="break"){
                        for(;forCmdIndex<cmdlistSize;forCmdIndex++){

                            if(DB) cout<<":::  "<<cmdlist[forCmdIndex]<<endl;

                            if(cmdlist[forCmdIndex]=="end" )
                            break;
                        }
                    }
                }


                /// remove local and for's iter  variables
                cmdIndex=forCmdIndex;
                uvars.erase(uvars.begin()+currPos,uvars.end());
        return true;
        }

        //-------------------------------------------------------------------



        //-------------------------------------------------------------------

        #ifdef __linux__


        bool doForInCommand(size_t &cmdIndex)
        {
        const string iterName(cmdlist[cmdIndex][1]);
        string forArg;

                    for(size_t i=3;i<cmdlist[cmdIndex].numOfKeyValues();i++)
                        forArg+=cmdlist[cmdIndex][i]+" ";

        FILE *fid=nullptr;
        char buff[1024];
        string fList;


                /////////////// open pipeline , store results to buffer, close pipeline
                fid=popen(forArg.c_str(),"r");

                if(!fid){
                    cerr<<"ERROR: stream failure"<<endl;
                return  false;
                }


                while(fgets(buff,1024,fid))
                    fList+=string(buff);

        const int status=pclose(fid);

                if(status){
                    cerr<<"ERROR: stream status failure"<<endl;
                return false;
                }
                ///////////////////////////////////////////////////////////////////////////

        std::size_t from=0,to;
        string fileName;
        const size_t currPos=uvars.size();
        size_t forCmdIndex=0;

                uvars.emplace_back(strpair(iterName,""));

                // do operations for each file
                do{
                    to=fList.find_first_of('\n',from);

                    if (to==string::npos) break;

                    fileName=fList.substr(from,to-from);

                    if(DB) cout<<"  "<<fileName<<endl;

                    auto locIterVar=std::find(uvars.begin(),uvars.end(),iterName);
                    locIterVar->getValue()=fileName;

                    forCmdIndex=cmdIndex+1; // always reset loop to the initial statement
                    if(!doCommands(forCmdIndex))
                        return false;

                    from=to+1;
                }while(true);


                // remove local variables if any
                cmdIndex=forCmdIndex;
                uvars.erase(uvars.begin()+currPos,uvars.end());

        return true;
        }

        #endif
        //=====================================================================
        bool isBlock(size_t &cmdIndex){
        return  cmdlist[cmdIndex]=="for"   ||
                cmdlist[cmdIndex]=="if"    ||
                cmdlist[cmdIndex]=="grain" ||
                cmdlist[cmdIndex]=="pdh"   ||
                cmdlist[cmdIndex]=="diff"  ||
                cmdlist[cmdIndex]=="gr" ;
        }

        //=====================================================================
        bool doIgnoreCommands(size_t &cmdIndex)
        {
                cmdIndex++;

                for(;cmdIndex<cmdlistSize;cmdIndex++){
                    if(DB) cout<<":::  "<<cmdlist[cmdIndex]<<endl;

                    if(cmdlist[cmdIndex]=="end")
                        break;
                    else                      
                        if( isBlock(cmdIndex) )
                            doIgnoreCommands(cmdIndex);

                }

        return true;
        }        

        //=====================================================================
        bool subExpr(const string &expr)
        {
        vector<string> tokPrm(split<string>(expr,"!=<>"));

        string var(tokPrm[0]);
        string arg(tokPrm[1]);

             Script::replaceVars(&uvars,var);
             Script::replaceVars(&uvars,arg);

             if(expr.find("==")!=string::npos)
                return std::stod(var)==std::stod(arg);

             if(expr.find("!=")!=string::npos)
                return std::stod(var)!=std::stod(arg);

             if(expr.find(">=")!=string::npos)
                return std::stod(var)>=std::stod(arg);

             if(expr.find("<=")!=string::npos)
                return std::stod(var)<=std::stod(arg);

             if(expr.find(">")!=string::npos)
                return std::stod(var)>std::stod(arg);

             if(expr.find("<")!=string::npos)
                return std::stod(var)<std::stod(arg);


             cerr<<"ERROR, unknown operator"<<endl;
             throw Script::ERR_UNK_VAR;
        }

        //.................................................................

        bool doIfCommand(size_t &cmdIndex)
        {
        const string cmd(cmdlist[cmdIndex][1]);
        vector<string> tokensCmd(split<string>(cmd," "));
        //const size_t tokSize=tokensCmd.size();
        vector<bool> subExprRet;
        size_t i,j;


                //// conditional expression results
                for(i=1;i<cmdlist[cmdIndex].numOfKeyValues();i+=2){
                    subExprRet.push_back(subExpr(cmdlist[cmdIndex][i]));
                }


        bool result=subExprRet[0];

                for(i=2,j=1;i<cmdlist[cmdIndex].numOfKeyValues();i+=2,j++){

                    if(cmdlist[cmdIndex][i]=="&")
                        result = (result && subExprRet[j]);
                    if(cmdlist[cmdIndex][i]=="|")
                        result = (result || subExprRet[j]);

                }

                cmdIndex++;

                if(result){                    
                    if(!doCommands(cmdIndex))
                        return false;

                    if(cmdlist[cmdIndex]=="else"){
                        cmdIndex++;

                        for(;cmdIndex<cmdlistSize;cmdIndex++){

                            if(DB) cout<<":::  "<<cmdlist[cmdIndex]<<endl;

                            if(isBlock(cmdIndex)){
                                doIgnoreCommands(cmdIndex);
                            continue;
                            }

                            if(cmdlist[cmdIndex]=="end" )
                            break;
                        }

                    }
                }
                else{

                    for(;cmdIndex<cmdlistSize;cmdIndex++){

                        if(DB) cout<<":::  "<<cmdlist[cmdIndex]<<endl;

                        if(isBlock(cmdIndex)){
                            doIgnoreCommands(cmdIndex);
                        continue;
                        }


                        if(cmdlist[cmdIndex]=="else"){
                            cmdIndex++;
                            if(!doCommands(cmdIndex))
                                return false;
                            else
                                break;
                        }


                        if(cmdlist[cmdIndex]=="end" )
                        break;
                    }
                }

        return true;
        }

        //=====================================================================
};


//=============================================================================

int main(int argc, char *argv[])
{
ClMainTask abc;

        if(DB)
        system("pwd");

        if(!abc.argProcessing(argc,argv)){
            cerr<<"  "<<endl;
        return -1;
        }

#ifdef __linux__
        if(!DB)
        system("setterm -linewrap off --cursor off");
#endif


void (ClMainTask::*pa)()=&ClMainTask::operator ();
std::thread mainTask(pa,&abc);

        mainTask.join();


#ifdef __linux__
        if(!DB)
        system("setterm -linewrap on --cursor on");
#endif


return abc.retvalue;
}


//=============================================================================

