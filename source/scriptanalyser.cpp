/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * scriptanalyser.cpp
 * Copyright (C) 2019 Kazimierz Skrobas <kskrobas@unipress.waw.pl>
 *
 * NPCL is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * NPCL is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "scriptanalyser.h"

#include "colormsg.h"

#include <regex>

#define RE_NUMBER_0  "[[:s:]]+[+-]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?"
#define RE_NUMBER_1  "[[:s:]]*[+-]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?"
#define PRE_NUMBER_0   "[[:s:]]+[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?"
#define PRE_NUMBER_1   "[[:s:]]*[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?"
#define PRE_NUMBER_2   "[[:s:]]*[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?[[:s:]]*"
#define UINT_NUMBER  "[[:s:]]+[0-9]+"
#define VAR          "\\$\\{\\w+\\}"

const std::string sRE_NUMBER(RE_NUMBER_0);    //real number
const std::string sRE_NUMBER_1(RE_NUMBER_1);    //real number
const std::string sPRE_NUMBER(PRE_NUMBER_0); //positive, real number
const std::string sPRE_NUMBER_1(PRE_NUMBER_1); //positive, real number
const std::string sPRE_NUMBER_2(PRE_NUMBER_2); //positive, real number
const std::string sUINT_NUMBER(UINT_NUMBER);
const std::string sVAR(VAR);


#ifdef DEBUG
#define DB true
#else
#define DB false
#endif

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
ClKeyValues &ClKeyValues::operator <<(const char *ckv)
{
const string kv(ckv);

        if(keyvalues.empty()){
            keyvalues.push_back(kv);
            ptr_frontValue=&keyvalues.back();
         }
        else
            keyvalues.push_back(kv);

        ptr_lastValue=&keyvalues.back();

        return *this;
}
////-----------------------------------------------------------------------------
ClKeyValues ClKeyValues::operator =(const ClKeyValues &b)
{

        keyvalues.clear();
        keyvalues.resize(b.keyvalues.size());

size_t i=0;
        for(auto &val:b.keyvalues)
             keyvalues[i++]=val;

        ptr_frontValue=&keyvalues.front();
        ptr_lastValue=&keyvalues.back();


return *this;
}
//-----------------------------------------------------------------------------
ClKeyValues &ClKeyValues::operator <<(const  vector<string>  & vs)
{
            keyvalues=std::move(vs);
            ptr_frontValue=& keyvalues.front();
            ptr_lastValue= & keyvalues.back();

return *this;
}

ClKeyValues &ClKeyValues::operator <<(string & kv)
{
        if(keyvalues.empty()){
            keyvalues.push_back(kv);
            ptr_frontValue=&keyvalues.back();
         }
        else
            keyvalues.push_back(kv);

        ptr_lastValue=&keyvalues.back();


return *this;
}
//-----------------------------------------------------------------------------
ClKeyValues &ClKeyValues::operator >>(string &kv)
{

        if(ptr_frontValue <= ptr_lastValue ){
            kv=*ptr_frontValue;
            ptr_frontValue++;
         }
        else
            kv.clear();

        return *this;
}
//-----------------------------------------------------------------------------
ClKeyValues &ClKeyValues::operator >>(vector<string> &vs)
{
        vs=keyvalues;
        return *this;
}
//-----------------------------------------------------------------------------
ostream &operator <<(ostream &os, const ClKeyValues &kv)
{

        for(auto &str: kv.keyvalues)
            os<<str<<" ";
return os;
}
//-----------------------------------------------------------------------------
ClKeyValues::ClKeyValues(const ClKeyValues &v)
{
       //cout<<"copy constr"<<endl;
        keyvalues=v.keyvalues;
        ptr_frontValue= & keyvalues.front();
        ptr_lastValue= & keyvalues.back();

}

ClKeyValues::ClKeyValues(ClKeyValues &&v)
{
       //cout<<"copy constr"<<endl;
        keyvalues=std::move(v.keyvalues);
        ptr_frontValue= & keyvalues.front();
        ptr_lastValue= & keyvalues.back();

}

void ClKeyValues::dispValues()
{
        for(auto & str : keyvalues)
            cout<<str<<" ";
        cout<<endl;
}
//-----------------------------------------------------------------------------


void appKeyValues(vcmdlist &ptr_cl,const std::string &cl)
{
vector<string> tokens(split<string> (cl," \t"));
ClKeyValues kv;

            kv<<tokens;
            ptr_cl.emplace_back(kv);
}
//-----------------------------------------------------------------------------
void appKey(vcmdlist &ptr_cl,const std::string &key)
{
ClKeyValues kv(key);
        ptr_cl.emplace_back(kv);
}

//-----------------------------------------------------------------------------
bool fCmpKeyString(ClKeyValues &kv,const std::string &key)
{ return (kv<key);}

//-----------------------------------------------------------------------------
void testDuplicate(vcmdlist &ptr_cl,const std::string key)
{
auto testKey=[key](ClKeyValues &kv){return (kv==key);};
        if(std::any_of(ptr_cl.begin(),ptr_cl.end(),testKey))
            throw Script::Result::DUP_ERR;
}
//-----------------------------------------------------------------------------
void testVariables(const std::string * cmd)
{
const std::string scmd(*cmd);

bool prtOpen=false;

            for(size_t i=0;i<scmd.length();i++){

                    if(scmd[i]=='{'){
                        if(!prtOpen){
                            prtOpen=true;
                        continue;
                        }
                        else{
                            cerr<<"ERROR: unbalanced paranthesis {"<<endl;
                        throw Script::ERR_PRT;
                        }
                    }


                    if(scmd[i]=='}'){
                        if(prtOpen)
                            prtOpen=false;
                        else{
                            cerr<<"ERROR: unbalanced paranthesis }"<<endl;
                        throw Script::ERR_PRT;
                        }
                    }
            }
}





//-----------------------------------------------------------------------------

void avePdhBlock(fstream &script, vcmdlist *ptr_cl , string &cmdline, const size_t options,size_t &cline)
{
vcmdlist apdh_cmdlist;

            apdh_cmdlist.reserve(2);



                while(!script.eof()){

                    std::getline(script,cmdline);
                    cline++;

                    trim(cmdline);

                    if(DB){ cout<<cline<<": "<<cmdline<<endl;}

                    //-- empty line
                    if(cmdline.empty()) continue;
                    //--

                    //-- comments
                    if(cmdline[0]=='#' || cmdline[0]=='%') continue;
                    //--

                    // replace tabs by spaces
                     std::replace(std::begin(cmdline),std::end(cmdline),'\t',' ');

                    //----------------//----------------//----------------//----------------

                    if(regex_match(cmdline,regex("files2ave[[:s:]]+[0-9]+"))){
                        testVariables(&cmdline);
                        appKeyValues(*ptr_cl,cmdline);
                    }


                    //----------------

                    if(regex_match(cmdline,regex("weight[[:s:]]+(uniform|normal|lognormal)("+sPRE_NUMBER+"|[[:s:]]+\\$\\{\\w+\\}){2}"))){
                    ClKeyValues kv;
                    string value(cmdline.substr(7));

                                kv<<"weight"<<value;
                                ptr_cl->emplace_back(kv);
                    continue;
                    }

                    //file operations open filename
                    if(regex_match(cmdline,regex("open[[:s:]]+[[:print:]]+"))){
                            testVariables(&cmdline);
                    ClKeyValues kv;
                    std::string value(cmdline.substr(5));

                            kv<<"open"<<value;
                            ptr_cl->emplace_back(kv);

                    continue;
                    }

                    //file operations   save filename
                    if(regex_match(cmdline,regex("save[[:s:]]+[[:print:]]+"))){
                            testVariables(&cmdline);
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("print[pP]rm"))){
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }


                    //the end of block
                    if(regex_match(cmdline,std::regex("end"))){
                    break;
                    }

                    throw Script::Result::MIS_UNK_ERR;
                }


                //inserting the current list to the main list
                for(auto &kv: apdh_cmdlist)
                        ptr_cl->emplace_back(kv);

}

//-----------------------------------------------------------------------------
void mergePdhBlock(fstream &script, vcmdlist *ptr_cl , string &cmdline, const size_t options,size_t &cline)
{
vcmdlist mpdh_cmdlist;

            mpdh_cmdlist.reserve(5);



                while(!script.eof()){

                    std::getline(script,cmdline);
                    cline++;

                    trim(cmdline);

                    if(DB){ cout<<cline<<": "<<cmdline<<endl;}

                    //-- empty line
                    if(cmdline.empty()) continue;
                    //--

                    //-- comments
                    if(cmdline[0]=='#' || cmdline[0]=='%') continue;
                    //--

                    // replace tabs by spaces
                     std::replace(std::begin(cmdline),std::end(cmdline),'\t',' ');


                    if(regex_match(cmdline,std::regex("flimit[[:s:]]+[0-9]+"))){
                        testDuplicate(mpdh_cmdlist,"flimit");
                        appKeyValues(mpdh_cmdlist,cmdline);
                    continue;
                    }


                    //-- range
                    if(regex_match(cmdline,std::regex("range("+sPRE_NUMBER+"){3}"))){
                            testDuplicate(mpdh_cmdlist,"range");
                            appKeyValues(mpdh_cmdlist,cmdline);

                    continue;
                    }


                    //file operations open filename
                    if(regex_match(cmdline,regex("open[[:s:]]+[[:print:]]+"))){
                            testVariables(&cmdline);
                    ClKeyValues kv;
                    std::string value(cmdline.substr(5));

                            kv<<"open"<<value;
                            mpdh_cmdlist.emplace_back(kv);

                    continue;
                    }

                    //file operations   save filename
                    if(regex_match(cmdline,regex("save[[:s:]]+[[:print:]]+"))){
                            testVariables(&cmdline);
                            appKeyValues(mpdh_cmdlist,cmdline);
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("print[pP]rm"))){
                            appKeyValues(mpdh_cmdlist,cmdline);
                    continue;
                    }

                    //----------------

                    if(regex_match(cmdline,regex("weight[[:s:]]+(uniform|normal|lognormal)("+sPRE_NUMBER+"|[[:s:]]+\\$\\{\\w+\\}){2}"))){
                    ClKeyValues kv;
                    string value(cmdline.substr(7));

                                kv<<"weight"<<value;
                                mpdh_cmdlist.emplace_back(kv);
                    continue;
                    }



                    //the end of block
                    if(regex_match(cmdline,std::regex("end"))){
                    break;
                    }

                    throw Script::Result::MIS_UNK_ERR;
                }


                //inserting the current list to the main list
                for(auto &kv: mpdh_cmdlist)
                        ptr_cl->emplace_back(kv);


}
//-----------------------------------------------------------------------------
void ucBlock(fstream &script, vcmdlist *ptr_cl , string &cmdline, const size_t options,size_t &cline)
{
int vx,vy,vz,atoms;
    vx=vy=vz=atoms=0;

vcmdlist uc_cmdlist;
            uc_cmdlist.reserve(5);

            while(!script.eof()){
                std::getline(script,cmdline);
                cline++;

                if(DB){ cout<<cline<<": "<<cmdline<<endl;}

                trim(cmdline);

                //-- empty line
                if(cmdline.empty()) continue;
                //--

                //-- comments
                if(cmdline[0]=='#' || cmdline[0]=='%') continue;
                //--

                // replace tabs by spaces
                 std::replace(std::begin(cmdline),std::end(cmdline),'\t',' ');


                if(regex_match(cmdline,std::regex("v[xyz]("+sRE_NUMBER+"|[[:s:]]+"+sVAR+"){3}"))){
                    if(cmdline.find("vx")!=string::npos) vx++;
                    if(cmdline.find("vy")!=string::npos) vy++;
                    if(cmdline.find("vz")!=string::npos) vz++;

                    appKeyValues(uc_cmdlist,cmdline);
                continue;
                }

                if(regex_match(cmdline,std::regex("[A-Z][[:w:]]*("+sRE_NUMBER+"|[[:s:]]+"+sVAR+"){3}"))){
                    appKeyValues(uc_cmdlist,cmdline);
                    atoms++;
                continue;
                }

                //the end of block
                if(regex_match(cmdline,std::regex("end"))){
                break;
                }

            throw Script::Result::MIS_UNK_ERR;
            }




            if(vx!=1  || vy!=1  || vz!=1){
               errMsg("wrong number of translation vectors");
            throw Script::Result::ERR_TRVUC;
            }

            if(atoms<1){
               errMsg("wrong number of base atoms");
            throw Script::Result::ERR_NAUC;
            }


            //inserting the current list to the main list
            for(auto &kv: uc_cmdlist)
                    ptr_cl->emplace_back(kv);


}
//-----------------------------------------------------------------------------
void cshBlock(fstream &script, vcmdlist *ptr_cl , string &cmdline, const size_t options,size_t &cline)
{

            while(!script.eof()){

                std::getline(script,cmdline);
                cline++;

                if(DB){ cout<<cline<<": "<<cmdline<<endl;}

                trim(cmdline);

                //-- empty line
                if(cmdline.empty()) continue;
                //--

                //-- comments
                if(cmdline[0]=='#' || cmdline[0]=='%') continue;
                //--


                // replace tabs by spaces
                 std::replace(std::begin(cmdline),std::end(cmdline),'\t',' ');


                if(regex_match(cmdline,std::regex(sRE_NUMBER_1+"("+sPRE_NUMBER+"|[[:s:]]+inf)"))){
                    ptr_cl->push_back(cmdline);
                continue;
                }


                if(regex_match(cmdline,std::regex("rand("+sRE_NUMBER+"){2}[[:s:]]+rand"+sPRE_NUMBER+"("+sPRE_NUMBER+"|[[:s:]]+inf)"))){
                    ptr_cl->push_back(cmdline);
                continue;
                }

                //the end of block
                if(regex_match(cmdline,std::regex("end"))){
                break;
                }

                throw Script::Result::MIS_UNK_ERR;
            }


}
//-----------------------------------------------------------------------------
void grainBlock(fstream &script, vcmdlist *ptr_cl , string &cmdline, const size_t options,size_t &cline,stdumap *ptr_uvar)
{
    ///
    /// duplicate test needs separate gb_cmdlist. If OK gb_cmdlist's entries are copied to cmdlist
    ///
vcmdlist gb_cmdlist;
bool radiusOrside=false;


                gb_cmdlist.reserve(10);


                while(!script.eof()){

                    std::getline(script,cmdline);
                    cline++;

                    if(DB){ cout<<cline<<": "<<cmdline<<endl;}

                    trim(cmdline);

                    //-- empty line
                    if(cmdline.empty()) continue;
                    //--

                    //-- comments
                    if(cmdline[0]=='#' || cmdline[0]=='%') continue;

                    // replace tabs by spaces
                     std::replace(std::begin(cmdline),std::end(cmdline),'\t',' ');
                    //--


                    //grain
                    if(regex_match(cmdline,std::regex("atom[[:s:]]+\\w+"))){
                            testDuplicate(gb_cmdlist,"atom");
                            testDuplicate(gb_cmdlist,"atoms");
                            appKeyValues(gb_cmdlist,cmdline);                            
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("atoms[[:s:]]+\\w+[[:s:]]+\\w+"))){
                            testDuplicate(gb_cmdlist,"atom");
                            testDuplicate(gb_cmdlist,"atoms");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("buildstyle[[:s:]]+(md|brav)"))){
                            testDuplicate(gb_cmdlist,"buildstyle");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("center[[:s:]]+(geom|catom([[:s:]]+\\w+)?|id[[:s:]]+[0-9]+)"))){
                            testDuplicate(gb_cmdlist,"center");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("charge[[:s:]]+\\w+"+sRE_NUMBER))){
                            //testDuplicate(gb_cmdlist,"mass");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("disloc[[:s:]]+(dumbell|intdef)[[:s:]]+uniform"+sRE_NUMBER+sPRE_NUMBER+sPRE_NUMBER+"([[:s:]]+(x|y|z))?"))){
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("disperse[[:s:]]+((\\*|\\w+)"+sPRE_NUMBER+")"))){
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("geometry[[:s:]]+(sphere|cube|cylinder)"))){
                            testDuplicate(gb_cmdlist,"geometry");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("geometry[[:s:]]+(cubic|oct|dod)"+sPRE_NUMBER+"("+sPRE_NUMBER+"){0,3}"))){
                            testDuplicate(gb_cmdlist,"geometry");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,
                            std::regex("geometry[[:s:]]+poly([[:s:]]+rand\\("+
                                       sPRE_NUMBER_1+","+sPRE_NUMBER_1+"\\)|"+sPRE_NUMBER_1+"|[[:s:]]+"+sVAR+"){3}("+
                                       sPRE_NUMBER_1+"|[[:s:]]+"+sVAR+"){3}"))){
                            testDuplicate(gb_cmdlist,"geometry");

                    const size_t prmPos(cmdline.find("poly"));
                    ClKeyValues kv;

                            kv<<"geometry"<<"poly"<<cmdline.substr(prmPos+5).c_str();
                            gb_cmdlist.emplace_back(kv);
                            //appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("hcpu"+sPRE_NUMBER))){
                            testDuplicate(gb_cmdlist,"hcpu");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("hcpcs[[:s:]]+(circle|hex|rhomb)?"))){
                            testDuplicate(gb_cmdlist,"hcpcs");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("hcpcs[[:s:]]+poly("+sPRE_NUMBER_1+"|[[:s:]]+"+sVAR+"){3}"))){
                            testDuplicate(gb_cmdlist,"hcpcs");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("hcpabc[[:s:]]+[0-9A-Ca-c\\(\\)"+sVAR+"]+"))){
                            testDuplicate(gb_cmdlist,"hcpabc");
                            testDuplicate(gb_cmdlist,"hcpfillup");

                    vector<string> tokens(split<string> (cmdline," "));                                                               
                    //const string abcexpr=parseABC(tokens[1]);
                    ClKeyValues kv;

                            //kv<<"hcpabc"<<"abc"<<abcexpr.c_str();
                            kv<<"hcpabc"<<"abc"<<tokens[1];
                            gb_cmdlist.emplace_back(kv);
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("hcpabc[[:s:]]+(uniform|normal|lognormal)("+sPRE_NUMBER+"){2}[[:s:]]+(abc|ab|abac|abcacb|random)"))){
                            testDuplicate(gb_cmdlist,"hcpabc");
                            testDuplicate(gb_cmdlist,"hcpfillup");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("hcpfillup[[:s:]]+(ABC|AB|ABAC|ABCACB)+([[:s:]]+[0-9]+)?"))){
                            testDuplicate(gb_cmdlist,"hcpfillup");
                            testDuplicate(gb_cmdlist,"hcpabc");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("hcpsl[[:s:]]+(yes|no)?"))){
                        appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("hcpsurfA([[:s:]]+(yes|no|H))?"))){
                        testDuplicate(gb_cmdlist,"hcpsurfA");
                        appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    ////////////
                    if(regex_match(cmdline,std::regex("insfault([[:s:]]+random[[:s:]]+[[:d:]]+|([[:s:]][[:d:]]+))"))){

                            appKeyValues(gb_cmdlist,cmdline);
                            if(gb_cmdlist.back()[2]=="0")
                                cerr<<"WARNING: number of stacking faults is equal to 0";
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("insfault[[:s:]]+random[[:s:]]+"+sVAR))){
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }



                    if(regex_match(cmdline,std::regex("lmpstyle[[:s:]]+(atomic|charge)"))){
                            testDuplicate(gb_cmdlist,"lmpstyle");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }
                    
                    
                    if(regex_match(cmdline,std::regex("lmpmargin("+sPRE_NUMBER+"){1,3}"))){
                            testDuplicate(gb_cmdlist,"lmpmargin");
                            gb_cmdlist.emplace_back(ClKeyValues(strpair("lmpmargin",cmdline.substr(10))));
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("lp("+sPRE_NUMBER+"|[[:s:]]+"+sVAR+")"))){
                            testDuplicate(gb_cmdlist,"lp");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }



                    if(regex_match(cmdline,std::regex("lp[[:s:]]+(uniform|normal|lognormal)"+sPRE_NUMBER+sPRE_NUMBER))){
                            testDuplicate(gb_cmdlist,"lp");
                            gb_cmdlist.emplace_back(ClKeyValues(strpair("lp",cmdline.substr(3))));
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("mass[[:s:]]+[0-9]"+sPRE_NUMBER+""))){
                            //testDuplicate(gb_cmdlist,"mass");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("numOfatomsTest[[:s:]]+(no|yes)?"))){
                            //testDuplicate(gb_cmdlist,"mass");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("open[[:s:]]+[[:print:]]+"))){
                            testVariables(&cmdline);
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("print[pP]rm"))){
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("push[[:s:]]+(numAtoms|[xyz]size)[[:s:]]+[[:print:]]+"))){
                        appKeyValues(gb_cmdlist,cmdline);
                    vector<string> tokens(split<string> (cmdline," "));

                        ptr_uvar->emplace_back(strpair(tokens[2],"0"));

                    continue;
                    }



                   if(regex_match(cmdline,std::regex("radius("+sPRE_NUMBER+"(lp)?|[[:s:]]+"+sVAR+"(lp)?)"))){

                           if(radiusOrside){
                               cerr<<"ERROR: the size of grain was defined earlier"<<endl;
                               throw Script::Result::DUP_ERR;
                           }

                            radiusOrside=true;

                            testDuplicate(gb_cmdlist,"radius");
                            appKeyValues(gb_cmdlist,cmdline);

                    continue;
                    }

                   if(regex_match(cmdline,std::regex("radius[[:s:]]+(uniform|normal|lognormal)"+sPRE_NUMBER+sPRE_NUMBER))){
                           testDuplicate(gb_cmdlist,"radius");
                            gb_cmdlist.emplace_back(ClKeyValues(strpair("radius",cmdline.substr(7))));
                   continue;
                   }


                   if(regex_match(cmdline,std::regex("replicate("+sUINT_NUMBER+"|[[:s:]]+"+sVAR+"){3}"))){
                           // testDuplicate(gb_cmdlist,"replace");
                            appKeyValues(gb_cmdlist,cmdline);
                   continue;
                   }


                   if(regex_match(cmdline,std::regex("rename([[:s:]]+[[:print:]]+[[:s:]]+[[:print:]]+("+sPRE_NUMBER+"|"+sVAR+"))"))){

                            gb_cmdlist.emplace_back(ClKeyValues(strpair("rename",cmdline.substr(6))));
                   continue;
                   }


                    if(regex_match(cmdline,std::regex("rescale("+sPRE_NUMBER+")"))){
                            testDuplicate(gb_cmdlist,"rescale");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("rescale("+sPRE_NUMBER+"){3}"))){
                            testDuplicate(gb_cmdlist,"rescale");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("remove"+sUINT_NUMBER+sPRE_NUMBER+"("+sPRE_NUMBER+")?"))){
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    //file operations   save filename
                    if(regex_match(cmdline,regex("save[[:s:]]+[[:print:]]+"))){
                            testVariables(&cmdline);
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,regex("saveHeader[[:s:]]+[[:print:]]+"))){
                            testVariables(&cmdline);
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,regex("saveopt([[:s:]]+(min|max)+"+sPRE_NUMBER+")+"))){
                            testVariables(&cmdline);
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("threads[[:s:]]+[0-9]+"))){
                            testDuplicate(gb_cmdlist,"threads");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("struct[[:s:]]+(sc|bcc|fcc|hcp|zb|uo2|zb110|feni|uc)"))){
                            testDuplicate(gb_cmdlist,"struct");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("threads[[:s:]]+[0-9]+"))){
                            testDuplicate(gb_cmdlist,"threads");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("ucp"))){
                        appKeyValues(gb_cmdlist,cmdline);
                        ucBlock(script,&gb_cmdlist,cmdline,options,cline);
                        appKeyValues(gb_cmdlist,"end");
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("voids"+sPRE_NUMBER+"("+sPRE_NUMBER+")?"))){
                    vector<string> tokens(split<string>(cmdline," "));
                    const double value=std::stod(tokens[1]);

                        if(value<0 || value>1) {throw Script::Result::ERR_OUTOFRANGE;}


                        testDuplicate(gb_cmdlist,"voids");
                        appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("csh"))){
                        appKeyValues(gb_cmdlist,cmdline);
                        cshBlock(script,&gb_cmdlist,cmdline,options,cline);
                        appKeyValues(gb_cmdlist,"end");
                    continue;
                    }



                    if(regex_match(cmdline,std::regex("comment[[:s:]]+[[:print:]]*"))){
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    //the end of a grain block
                    if(regex_match(cmdline,std::regex("end"))){
                    break;
                    }


                    throw Script::Result::MIS_UNK_ERR;
                }


                //inserting the current list to the main list
                for(auto &kv: gb_cmdlist)
                        ptr_cl->emplace_back(kv);

}



//-----------------------------------------------------------------------------

void pdhBlock(fstream &script, vcmdlist *ptr_cl , string &cmdline, const size_t options,size_t &cline)
{
vcmdlist gb_cmdlist;

                gb_cmdlist.reserve(10);


                while(!script.eof()){

                    std::getline(script,cmdline);
                    cline++;

                    if(DB){ cout<<cline<<": "<<cmdline<<endl;}

                    trim(cmdline);

                    //-- empty line
                    if(cmdline.empty()) continue;
                    //--

                    //-- comments
                    if(cmdline[0]=='#' || cmdline[0]=='%') continue;
                    //--

                    // replace tabs by spaces
                     std::replace(std::begin(cmdline),std::end(cmdline),'\t',' ');

                    if(regex_match(cmdline,std::regex("bin"+sPRE_NUMBER+"(A|lp|nm)?"))){
                            testDuplicate(gb_cmdlist,"bin");
                            testDuplicate(gb_cmdlist,"rangeN");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("mode[[:s:]]+(single|double)"))){
                            testDuplicate(gb_cmdlist,"mode");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("threads[[:s:]]+[0-9]+"))){
                            testDuplicate(gb_cmdlist,"threads");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }
					
					if(regex_match(cmdline,std::regex("mthmode[[:s:]]+(openmp|stdthread)"))){
                            testDuplicate(gb_cmdlist,"mthmode");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    //file operations open filename
                    if(regex_match(cmdline,regex("open[[:s:]]+[[:print:]]+"))){
                            testVariables(&cmdline);
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("print[pP]rm"))){
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("plot([[:s:]]+[[:print:]]*)?"))){
                    ClKeyValues kv;
                    string plotPrm(cmdline.substr(4));

                            kv<<"plot"<<plotPrm;
                            gb_cmdlist.emplace_back(kv);

                    continue;
                    }


                    if(regex_match(cmdline,std::regex("rangeN"+sPRE_NUMBER+"[[:s:]]+[0-9]+("+sPRE_NUMBER+")?"))){
                    vector<string> tokens(split<string> (cmdline," "));

                            if(tokens[2]=="0")
                                throw Script::Result::ERR_VAL_0;

                            testDuplicate(gb_cmdlist,"rangeN");
                            testDuplicate(gb_cmdlist,"bin");

                    ClKeyValues kv;

                               kv<<tokens;
                               gb_cmdlist.emplace_back(kv);
                    continue;
                    }


                    //file operations   save filename
                    if(regex_match(cmdline,regex("save[[:s:]]+[[:print:]]+"))){
                            testVariables(&cmdline);
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,regex("saveFigure[[:s:]]+[[:print:]]+"))){
                            testVariables(&cmdline);
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }


                    //difftime
                    if(regex_match(cmdline,regex("difftime"))){
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }


                    //the end of block
                    if(regex_match(cmdline,std::regex("end"))){
                    break;
                    }


                    if(regex_match(cmdline,std::regex("comment[[:s:]]+[[:print:]]*"))){
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    throw Script::Result::MIS_UNK_ERR;
                }


                //inserting the current list to the main list
                for(auto &kv: gb_cmdlist)
                        ptr_cl->emplace_back(kv);
}
//-----------------------------------------------------------------------------
void ifBlock(fstream &script, size_t &cline , vcmdlist *ptr_cl , string &cmdline, const size_t options, stdumap *ptr_uvar)
{
vector<string> tokensCmd(split<string>(cmdline.substr(3)," "));
const size_t tokSize=tokensCmd.size();


            for(size_t i=0;i<tokSize;i+=2){


                if( regex_match(tokensCmd[i],std::regex("[[:print:]]+([<>]|[!=<>]=)[[:print:]]+")) ){
                vector<string> tokPrm(split<string>(tokensCmd[i],"!=><"));
                const string ssvarL(tokPrm[0].substr(2,tokPrm[0].length()-3));

                        // left side of operator
                        auto iter=std::find(ptr_uvar->begin(),ptr_uvar->end(),ssvarL);
                        if(iter==ptr_uvar->end())
                            throw Script::ERR_UNK_VAR;

                        // right side of operator
                        if(tokPrm[1][0]=='$'){
                        const string ssvarR(tokPrm[1].substr(2,tokPrm[1].length()-3));
                            iter=std::find(ptr_uvar->begin(),ptr_uvar->end(),ssvarR);
                            if(iter==ptr_uvar->end())
                            throw Script::ERR_UNK_VAR;
                         }
                        else{
                            if(!regex_match(tokPrm[1],std::regex(sRE_NUMBER_1)))
                            throw Script::ERR_NO_NUMBER;
                        }

                }
                else{
                    cerr<<"Error, "<<cline<<" illegal operator '"<<tokensCmd[i]<<"'"<<endl;
                throw Script::ERR0;
                }
            }


            for(size_t i=1;i<tokSize;i+=2)
                if(!regex_match(tokensCmd[i],std::regex("[&|]"))){
                    cerr<<"Error, "<<cline<<" illegal operator '"<<tokensCmd[i]<<"'"<<endl;
                throw Script::ERR0;
                }


            appKeyValues(*ptr_cl,cmdline);

            auto retValue=Script::scriptParsing(script,cline, ptr_cl,ptr_uvar,options);

            if(retValue==Script::Result::ENDELSE){
            auto retValueElse=Script::scriptParsing(script,cline, ptr_cl,ptr_uvar,options);

                if(retValueElse!=Script::Result::ENDRET){
                    cerr<<" ERROR: misplaced END"<<endl;
                throw Script::ERR0;
                }

            }
            else
                if(retValue!=Script::Result::ENDRET){
                    cerr<<" ERROR: misplaced END"<<endl;
                throw Script::ERR0;
            }

//vector<string> tokensArg(split<string>(tokensCmd[1],"&|"));
}
//-----------------------------------------------------------------------------
void forBlock(fstream &script, size_t &cline , vcmdlist *ptr_cl , string &cmdline, const size_t options, stdumap *ptr_uvar)
{
//vcmdlist gb_cmdlist;
vector<string> tokensCmd(split<string>(cmdline," "));
vector<string> tokensArg(split<string>(tokensCmd[1],"="));


                if(regex_match(cmdline,std::regex(string("([-]?[[:d:]]+|")+VAR+")(:"+"([-]?[[:d:]]+|"+VAR+")){1,2}"))){
                    cerr<<"Error, "<<cline<<" illegal loop parameters  "<<endl;
                throw Script::ERR0;
                }

//vector<string> tokensPrm(split<string>(tokensArg[1],":"));


                /// czy argumenty for loop są legalne
                ///if(!regex_match(tokensArg[1],std::regex("[-]?[0-9]+:[0-9]+"))){
                ///    cerr<<"Error, "<<cline<<" illegal loop parameters  "<<endl;
                ///throw Script::ERR0;
                ///}

                /// czy istnieje juz zmienna o tej samej nazwie
                ///
                const auto iter=std::find(ptr_uvar->begin(),ptr_uvar->end(),tokensArg[0]);

                if(iter!=ptr_uvar->end())
                throw Script::DUP_ERR;


const size_t currPos=ptr_uvar->size();

                ptr_uvar->emplace_back(strpair(tokensArg[0],tokensArg[1]));
                appKeyValues(*ptr_cl,cmdline);


                if(Script::scriptParsing(script,cline, ptr_cl,ptr_uvar,options)!=Script::Result::ENDRET){
                    cerr<<" ERROR: misplaced END"<<endl;
                    throw Script::ERR0;
                }

                ptr_uvar->erase(ptr_uvar->begin()+currPos,ptr_uvar->end()); /// delete local variables of for's loop
}

//-----------------------------------------------------------------------------

void forInBlock(fstream &script, size_t &cline , vcmdlist *ptr_cl , string &cmdline, const size_t options, stdumap *ptr_uvar)
{
vector<string> tokensCmd(split<string>(cmdline," "));

string forArg("for_in ");
size_t i;
            for(i=1;i<tokensCmd.size()-1;i++)
                forArg+=tokensCmd[i]+" ";
            forArg+=tokensCmd[i];

            if(DB)
            cout<<forArg<<endl;


const size_t currPos=ptr_uvar->size();

            appKeyValues(*ptr_cl,forArg);
            ptr_uvar->emplace_back(strpair(tokensCmd[1],forArg));

            if(Script::scriptParsing(script,cline, ptr_cl,ptr_uvar,options)!=Script::Result::ENDRET){
                cerr<<" ERROR: misplaced END"<<endl;
                throw Script::ERR0;
            }

            ptr_uvar->erase(ptr_uvar->begin()+currPos,ptr_uvar->end()); /// delete local variables of for's loop

}
//-----------------------------------------------------------------------------
void diffBlock(fstream &script, vcmdlist *ptr_cl , string &cmdline, const size_t options,size_t &cline)
{
vcmdlist gb_cmdlist;

                gb_cmdlist.reserve(10);


                while(!script.eof()){

                    std::getline(script,cmdline);
                    cline++;

                    trim(cmdline);

                    //-- empty line
                    if(cmdline.empty()) continue;
                    //--

                    //-- comments
                    if(cmdline[0]=='#' || cmdline[0]=='%') continue;

                    // replace tabs by spaces
                     std::replace(std::begin(cmdline),std::end(cmdline),'\t',' ');


                    if(regex_match(cmdline,std::regex("comment[[:s:]]+[[:print:]]*"))){
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    //the end of block
                    if(regex_match(cmdline,std::regex("end"))){
                    break;
                    }

                   //difftime
                    if(regex_match(cmdline,regex("difftime"))){
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }

                    //Debyea-Waller factor
                    if(regex_match(cmdline,regex("dw"+sPRE_NUMBER))){
                            testDuplicate(gb_cmdlist,"dw");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    if(regex_match(cmdline,regex("fastsinc[[:s:]]+(yes|no)"))){
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }


                    //-- extrapolate
                    if(regex_match(cmdline,std::regex("ex(trapolate)?[[:s:]]+(yes|no)"))){
                            testDuplicate(gb_cmdlist,"extrapolate");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    //-- Ki [kx, ky, kz]

                    if(regex_match(cmdline,std::regex("Ki("+sPRE_NUMBER+"){3}"))){
                            testDuplicate(gb_cmdlist,"Ki");
                    ClKeyValues kv;
                    string kiPrm(cmdline.substr(2));

                                    kv<<"Ki"<<kiPrm;
                                    gb_cmdlist.emplace_back(kv);
                    continue;
                    }


                    //-- wavelength
                    if(regex_match(cmdline,std::regex("lambda"+sPRE_NUMBER))){
                            testDuplicate(gb_cmdlist,"lambda");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    //-- wavelength
                    if(regex_match(cmdline,std::regex("lambda[[:s:]]+(Ag|Co|Cr|Cu|Fe|Mo|W)"))){
                            testDuplicate(gb_cmdlist,"lambda");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    //-- mode
                    if(regex_match(cmdline,std::regex("mode[[:s:]]+(deb|laue)"))){
                            testDuplicate(gb_cmdlist,"mode");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }



                    if(regex_match(cmdline,regex("noiseSQ[[:s:]]+(uniform|normal|lognormal)"+sPRE_NUMBER+sPRE_NUMBER ))){
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,regex("norm[[:s:]]+(yes|no)?"))){
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }

                    //if(regex_match(cmdline,regex("nosf[[:s:]]+(yes|no)"))){
                    //        appKeyValues(*ptr_cl,cmdline);
                    //continue;
                    //}

                    //--
                   if(regex_match(cmdline,std::regex("polarization[[:s:]]+(yes|no)"))){
                          // testDuplicate(gb_cmdlist,"polarization");
                           appKeyValues(gb_cmdlist,cmdline);
                   continue;
                   }

                    if(regex_match(cmdline,std::regex("plot([[:s:]]+[[:print:]]*)?"))){
                    ClKeyValues kv;
                    string plotPrm(cmdline.substr(4));

                            kv<<"plot"<<plotPrm;
                            gb_cmdlist.emplace_back(kv);

                    continue;
                    }

                    if(regex_match(cmdline,std::regex("print[pP]rm"))){
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    //-- radiation
                    if(regex_match(cmdline,std::regex("radiation[[:s:]]+(xray|neutron|electron|nt)"))){
                            testDuplicate(gb_cmdlist,"radiation");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    //-- range
                    if(regex_match(cmdline,std::regex("range("+sRE_NUMBER+"){3}([[:s:]]+(2th|Q))?"))){
                            testDuplicate(gb_cmdlist,"range");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    //file operations   save filename
                    if(regex_match(cmdline,regex("save[[:s:]]+[[:print:]]+"))){
                            testVariables(&cmdline);
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }
                    
                    //file operations   save filename
                    if(regex_match(cmdline,regex("saveopt[[:s:]]+[[:print:]]+"))){
                            testVariables(&cmdline);
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }

                     if(regex_match(cmdline,regex("filesft[[:s:]]+[[:print:]]+"))){
                            testVariables(&cmdline);
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }



                   //-- threads
                    if(regex_match(cmdline,std::regex("threads[[:s:]]+[0-9]+"))){
                            testDuplicate(gb_cmdlist,"threads");
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    throw Script::Result::MIS_UNK_ERR;
                }



                //inserting the current list to the main list
                for(auto &kv: gb_cmdlist)
                        ptr_cl->emplace_back(kv);

}
//-----------------------------------------------------------------------------
void grBlock(fstream &script, vcmdlist *ptr_cl , string &cmdline, const size_t options,size_t &cline)
{
vcmdlist gb_cmdlist;

                gb_cmdlist.reserve(10);

                while(!script.eof()){

                    std::getline(script,cmdline);
                    cline++;

                    trim(cmdline);

                    //-- empty line
                    if(cmdline.empty()) continue;
                    //--

                    //-- comments
                    if(cmdline[0]=='#' || cmdline[0]=='%') continue;

                    // replace tabs by spaces
                     std::replace(std::begin(cmdline),std::end(cmdline),'\t',' ');


                    if(regex_match(cmdline,std::regex("comment[[:s:]]+[[:print:]]*"))){
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    //cf
                     if(regex_match(cmdline,regex("cf"+sPRE_NUMBER))){
                             appKeyValues(*ptr_cl,cmdline);
                     continue;
                     }


                    //difftime
                     if(regex_match(cmdline,regex("difftime"))){
                             appKeyValues(*ptr_cl,cmdline);
                     continue;
                     }

                     //the end of block
                     if(regex_match(cmdline,std::regex("end"))){
                     break;
                     }

                    //
                    if(regex_match(cmdline,std::regex("norm"))){
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }

                    //
                    if(regex_match(cmdline,std::regex("print[pP]rm"))){
                            appKeyValues(gb_cmdlist,cmdline);
                    continue;
                    }


                    //
                    if(regex_match(cmdline,std::regex("plot([[:s:]]+[[:print:]]*)?"))){
                    ClKeyValues kv;
                    string plotPrm(cmdline.substr(4));

                            kv<<"plot"<<plotPrm;
                            gb_cmdlist.emplace_back(kv);

                    continue;
                    }

                     //-- range
                     if(regex_match(cmdline,std::regex("range("+sPRE_NUMBER+"){2}("+sPRE_NUMBER+")?"))){
                             testDuplicate(gb_cmdlist,"range");
                             appKeyValues(gb_cmdlist,cmdline);
                     continue;
                     }

                     //-- range
                     if(regex_match(cmdline,std::regex("rangeN"+sPRE_NUMBER+"[[:s:]]+[0-9]+("+sPRE_NUMBER+")?"))){
                     vector<string> tokens(split<string> (cmdline," "));

                             if(tokens[2]=="0"||tokens[2]=="1")
                                 throw Script::Result::ERR_VAL_0;

                             testDuplicate(gb_cmdlist,"rangeN");

                     ClKeyValues kv;

                                kv<<tokens;
                                gb_cmdlist.emplace_back(kv);
                     continue;
                     }



                     // windowing function
                     if(regex_match(cmdline,std::regex("wf[[:s:]]+(box|lorch)"))){
                             testDuplicate(gb_cmdlist,"wf");
                             appKeyValues(gb_cmdlist,cmdline);
                     continue;
                     }


                     //file operations   save filename
                     if(regex_match(cmdline,regex("save[[:s:]]+[[:print:]]+"))){
                             testVariables(&cmdline);
                             appKeyValues(*ptr_cl,cmdline);
                     continue;
                     }

                     if(regex_match(cmdline,regex("saveFigure[[:s:]]+[[:print:]]+"))){
                             testVariables(&cmdline);
                             appKeyValues(*ptr_cl,cmdline);
                     continue;
                     }



                    //-- threads
                     if(regex_match(cmdline,std::regex("threads[[:s:]]+[0-9]+"))){
                             testDuplicate(gb_cmdlist,"threads");
                             appKeyValues(gb_cmdlist,cmdline);
                     continue;
                     }

                     throw Script::Result::MIS_UNK_ERR;
            }


                //inserting the current list to the main list
                for(auto &kv: gb_cmdlist)
                        ptr_cl->emplace_back(kv);


}
//-----------------------------------------------------------------------------
Script::Result Script::scriptParsing(fstream &script, size_t &cline, vcmdlist *ptr_cl, stdumap *ptr_uvar, const size_t options)
{
std::string cmdline;
Script::Result scres=Script::Result::OK;

const size_t currPos=ptr_uvar->size();

        try{

                while(!script.eof()){

                    std::getline(script,cmdline);
                    cline++;

                    if(DB){ cout<<cline<<": "<<cmdline<<endl;}

                    trim(cmdline);

                    //--
                    if(cmdline.empty()) continue;
                    //--

                    //--            
                    // one line comment
                    if(cmdline[0]=='>' || cmdline[0]=='#' || cmdline[0]=='%') continue;

                    // block comment
                    if(cmdline.find("/>")!=string::npos){
                        while(cmdline.rfind("</")==string::npos && !script.eof() ){
                            std::getline(script,cmdline);
                            rtrim(cmdline);
                            cline++;
                        }

                        if(script.eof()){
                            errMsg(" block comment brackets are not balanced");
                            throw Script::Result::ERR_EOF;
                        }


                    continue;
                    }

                    //-

                    //the end of block
                    if(regex_match(cmdline,std::regex("end(for|if)?")) )
                    return Script::Result::ENDRET;

                    //the else of block
                     if(regex_match(cmdline,std::regex("else")) ){
                         appKeyValues(*ptr_cl,"else");
                     return Script::Result::ENDELSE;
                     }

                    // stop script parsing
                    if(regex_match(cmdline,std::regex("exit")))
                    break;

                    // replace tabs by spaces
                     std::replace(std::begin(cmdline),std::end(cmdline),'\t',' ');

                    //average pdh computation
                    if(regex_match(cmdline,std::regex("avepdh"))){
                        appKeyValues(*ptr_cl,"avepdh");
                        avePdhBlock(script,ptr_cl,cmdline,options,cline);
                        appKeyValues(*ptr_cl,"end");
                    continue;
                    }

                    //average pdh computation
                    if(regex_match(cmdline,std::regex("mergepdh"))){
                        appKeyValues(*ptr_cl,"mergepdh");
                        mergePdhBlock(script,ptr_cl,cmdline,options,cline);
                        appKeyValues(*ptr_cl,"end");
                    continue;
                    }



                    //budowa ziarna
                    if(regex_match(cmdline,std::regex("grain"))){
                        appKeyValues(*ptr_cl,"grain");
                        grainBlock(script,ptr_cl,cmdline,options,cline,ptr_uvar);
                        appKeyValues(*ptr_cl,"end");
                    continue;
                    }

                    //policz pdh
                    if(regex_match(cmdline,std::regex("pdh"))){
                        appKeyValues(*ptr_cl,"pdh");
                        pdhBlock(script,ptr_cl,cmdline,options,cline);
                        appKeyValues(*ptr_cl,"end");
                    continue;
                    }

                    //policz diff
                    if(regex_match(cmdline,std::regex("diff"))){
                        appKeyValues(*ptr_cl,"diff");
                        diffBlock(script,ptr_cl,cmdline,options,cline);
                        appKeyValues(*ptr_cl,"end");
                    continue;
                    }

                    //policz gr
                    if(regex_match(cmdline,std::regex("gr"))){
                        appKeyValues(*ptr_cl,"gr");
                        grBlock(script,ptr_cl,cmdline,options,cline);
                        appKeyValues(*ptr_cl,"end");
                    continue;
                    }

                    //system
                    if(regex_match(cmdline,std::regex("system([[:s:]]+[[:print:]]+)+"))){
                        testVariables(&cmdline);
                        appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }


                    if(regex_match(cmdline,std::regex("if[[:s:]]+[[:print:]]+"))){
                        ifBlock(script,cline,ptr_cl,cmdline,options,ptr_uvar);
                        appKeyValues(*ptr_cl,"end");
                    continue;
                    }



                    if(regex_match(cmdline,std::regex("for[[:s:]]+\\w+=[[:print:]]+"))){
                        forBlock(script,cline,ptr_cl,cmdline,options,ptr_uvar);
                        appKeyValues(*ptr_cl,"end");
                    continue;
                    }

                    #ifdef __linux__
                    if(regex_match(cmdline,std::regex("for[[:s:]]+\\w+[[:s:]]+in[[:s:]]+[[:print:]]+"))){
                        forInBlock(script,cline,ptr_cl,cmdline,options,ptr_uvar);
                        appKeyValues(*ptr_cl,"end");
                    continue;
                    }
                    #endif

                    //if(regex_match(cmdline,std::regex("clear[[:s:]]+SNA"))){ //CLEAR  Nanograins::savedNumOfAtoms
                    //    appKeyValues(*ptr_cl,cmdline);
                    //continue;
                    //}

                    if(regex_match(cmdline,std::regex("break"))){
                        appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("breakIfSNA[[:s:]]+(no|yes)?"))){
                    vector<string> cmdtok(split<string>(cmdline," "));
                            if(cmdtok[1]=="yes")
                                appKeyValues(*ptr_cl,"break");
                    continue;
                    }

                    ///cast2int
                    if(regex_match(cmdline,std::regex("cast2int[[:s:]]+[[:print:]]+"))){
                            testVariables(&cmdline);
                             appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }


                    /// print
                    if(regex_match(cmdline,std::regex("print([[:s:]]+[[:print:]]+)+"))){
                    //vector<string> tokens(split<string>(cmdline," "));                                                
                            appKeyValues(*ptr_cl,cmdline);
                    continue;
                    }

                    if(regex_match(cmdline,std::regex("clearSNA"))){
                        appKeyValues(*ptr_cl,"clearSNA");
                    continue;
                    }


                    //user string type variables
                    if(regex_match(cmdline,std::regex("\\w+[[:s:]]*=[[:s:]]*\"[[:print:]]+\""))){
                    vector<string> tokens(split<string> (cmdline,"="));
                    string suvar(tokens[1].substr(1,tokens[1].length()-2));// kasuje znaki "                    
                    //auto iter=std::find(ptr_uvar->begin(),ptr_uvar->end(),tokens[0]);


                            testVariables(&suvar);

                    ClKeyValues kv;
                                    kv<<"$var_lit"<<tokens[0]<<suvar;
                                    ptr_cl->emplace_back(kv);
                                    ptr_uvar->emplace_back(strpair(tokens[0],suvar));
                    continue;
                    }

                    //user math type variables
                    if(regex_match(cmdline,std::regex("\\w+[[:s:]]*=[[:s:]]*[[:print:]]+"))){
                    vector<string> tokens(split<string> (cmdline,"="));
                            testVariables(&tokens[1]);

                    ClKeyValues kv;

                            kv<<"$var_math"<<tokens[0]<<tokens[1];
                            ptr_cl->emplace_back(kv);
                            ptr_uvar->emplace_back(strpair(tokens[0],tokens[1]));

                    continue;
                    }

                    ////
                    if(regex_match(cmdline,std::regex("\\$\\{\\w+\\}(\\+\\+|--)"))){
                    vector<string> tokens(split<string>(cmdline,"{}"));
                    ClKeyValues kv;
                            kv<<"incdec"<<tokens[1]<<tokens[2];
                            ptr_cl->emplace_back(kv);
                    continue;
                    }


                    /////////////
                    if(regex_match(cmdline,std::regex("printvars"))){
                        appKeyValues(*ptr_cl,"printvars");
                    continue;
                    }




//                    if(regex_match(cmdline,regex("load[[:s:]]+positions[[:s:]]+[[:print:]]+"))){
//                            appKeyValues(*ptr_cl,cmdline);
//                    continue;
//                    }


                    throw Script::Result::MIS_UNK_ERR;
                }

        }
        catch(Script::Result scerr){

                cerr<<" Error, line: "<<cline<<", command: "<<cmdline<<" ::: ";

                switch(scerr){
                case Script::Result::BLOCK_ERR:      cerr<<endl;break;
                case Script::Result::MIS_UNK_ERR:    cerr<<"wrong format, misplaced or unknown command"<<endl;break;
                case Script::Result::DUP_ERR  :      cerr<<"duplication or ambiguity has been found"<<endl;break;
                case Script::Result::ERR_EXPR :      cerr<<"incorrect expression "<<endl;break;
                case Script::Result::ERR_VAL_0:      cerr<<"argument must be greater than 0"<<endl;break;
                case Script::Result::ERR_NO_NUMBER:  cerr<<"argument must be a number"<<endl;break;
                case Script::Result::ERR_UNK_VAR:    cerr<<"unknown variable "<<endl;break;
                case Script::Result::ERR_EOF:        cerr<<"unexpected end of file "<<endl;break;

                default: cerr<<"  unknown Script::Result,   code"<< scerr<< endl;
                }

                scres=Script::Result::ERR1;
        }
        catch (std::regex_error& e) {
                cerr<<"ERROR: regular expression failure, reg-exp error code: "<<e.code()<<endl;
                cerr<<"\t see: http://www.cplusplus.com/reference/regex/regex_error/"<<endl;
                scres=Script::Result::ERR1;
        }
        catch(...){
                std::exception_ptr p=std::current_exception();

                cerr<<" Error, line: "<<cline<<
                      ", name: "<<(p ? p.__cxa_exception_type()->name() : "null") <<
                      ", command: "<<cmdline<<":";
                scres=Script::Result::ERR1;
        }

        ptr_uvar->erase(ptr_uvar->begin()+currPos,ptr_uvar->end()); /// delete local variables of for's loop

return scres;
}

//-----------------------------------------------------------------------------
struct StIsBrace{
    bool operator () (const char &ch) const {return (ch=='{') || (ch =='}') ;}
    typedef int argument_type;
};

/*
static inline std::string &ltrimBrace(std::string &s)
{
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(StIsBrace())));
        return s;
}

// trim from end
static inline std::string &rtrimBrace(std::string &s)
{
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(StIsBrace())).base(), s.end());
        return s;
}
*/
//-----------------------------------------------------------------------------
Script::ResultRepVar Script::replaceVars(stdumap *ptr_uvar, string &value)
{
            if(value.find("$")==string::npos)
            return Script::EMPTY;


            for(ClKeyValues & uvar: *ptr_uvar){
            std::string uv("${"+uvar.getKey()+"}");
            string::size_type pos=0;

                do{
                    pos=value.find(uv,pos);

                    if( pos!=string::npos)
                        value.replace(pos,uv.length(),uvar.getValue());

                }while(pos!=string::npos);
            }


            if(value.find("$")!=string::npos){
                if(DB){
                    for(ClKeyValues & uvar: *ptr_uvar)
                        cout<<uvar.getKey()<<":"<<uvar.getValue()<<endl;
                }


                errMsg(value);
            throw Script::FAIL;
            }


return Script::SUCC;
}


