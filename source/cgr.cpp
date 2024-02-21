/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * cgr.cpp
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


#include "cgr.h"

#include <omp.h>
#include <iomanip>
#include <chrono>
#include <thread>
#include <future>
#include <assert.h>

#include "cprogress.h"
#include "colormsg.h"
#include "createdir.h"


#ifdef DEBUG
#define DB true
#else
#define DB false
#endif

//=============================================================================

Cgr::Cgr()
{

}

//=============================================================================
void Cgr::saveResults()
{

        for(string file : fileNameOut){

                if(!createDirsIfDontExist(file)){
                    errMsg("couldn't create nested directories for "+file);
                throw Cgr::ERR_FILEOPEN;
                }

                if(file.rfind(".dat")!=string::npos)
                    saveDatFile(file);
                else
                    if(file.rfind(".grb")!=string::npos)
                        saveBinFile(file);
                    else
                        if(file.rfind(".gr")!=string::npos)
                            saveGrFile(file);
                        else
                            cerr<<"ERROR: Unknown file format"<<endl;
        }
}
//=============================================================================
void Cgr::saveDatFile(const string &fileName)
{
fstream file(fileName,ios::out);

            if(!file){
                cerr<<"ERROR: couldn't open file"<<endl;
            throw Cgr::ERR_FILEOPEN;
            }

position *x=dataX.values;
position *y=dataY.values;
const size_t dataSize=dataX.size;

            for(size_t i=0;i<dataSize;i++,x++,y++)
                file<<(*x)<<"\t"<<(*y)<<endl;

            file.close();
}
//=============================================================================
void Cgr::saveGrFile(const string &fileName)
{
fstream file(fileName,ios::out);

            if(!file){
                cerr<<"ERROR: couldn't open file"<<endl;
            throw Cgr::ERR_FILEOPEN;
            }


position *x=dataX.values;
position *y=dataY.values;
const size_t dataSize=dataX.size;
const size_t numOfAtoms=diff->pdh->grain->atoms.size();

            file<<"#ver: 01"<<endl;
            file<<"#title: gr data"<<endl;
            file<<"#numOfAtoms: "<<numOfAtoms<<endl;
            file<<"#sizeRC: "<<dataSize<<"\t2"<<endl;
            file<<(*this);
            file<<"#r\tG(r)"<<endl;

            for(size_t i=0;i<dataSize;i++,x++,y++)
                file<<(*x)<<"\t"<<(*y)<<endl;

            file.close();
}
//=============================================================================
void Cgr::saveBinFile(const string &fileName)
{
fstream file(fileName,ios::out | ios::binary);


            if(!file){
                cerr<<"ERROR: couldn't open file"<<endl;
            throw Cgr::ERR_FILEOPEN;
            }


const position *x=dataX.values;
const position *y=dataY.values;
const size_t dataSize=dataX.size;


const string version("v00");
const string data("DATA");


            file.write( (char *) version.c_str(), version.length()*sizeof(char));
            file.write( (char *) &dataSize,       sizeof(size_t));
            file.write( (char *) data.c_str(),    data.length()*sizeof(char) );

            for(size_t i=0;i<dataSize;i++,x++,y++){
                file.write( (char *) x, sizeof (position));
                file.write( (char *) y, sizeof (position));
            }

            file.close();
}
//=============================================================================
void Cgr::dispPrm()
{

        cout<<"*** G(r) parameters ***"<<endl;
        cout<<(*this);

}
//=============================================================================
void Cgr::normGr()
{
double sum=0;
position *y=dataY.values;
const size_t dataSize=dataX.size;
//auto sqrd=[](const position &x){return x*x;};

            for(size_t i=0;i<dataSize;i++,y++)
                sum+=std::abs((*y));

const double isum=1.0/sum;


            y=dataY.values;
            for(size_t i=0;i<dataSize;i++,y++)
                (*y)*=isum;
}

//=============================================================================

ostream & operator<<(ostream &o, Cgr &cgr)
{
        o<<"#range: "<<cgr.range<<endl;
        o<<"#wf: "<<cgr.wf<<endl;
        o<<"#threads: "<<cgr.threads<<endl;
return o;
}



//=============================================================================
bool Cgr::parseCommands(vcmdlist &cmd, size_t &index, stdumap *uvar__)
{
const str send("end");

                ptr_uvar=uvar__;
                presetData();



                while(cmd[index]!=send){

                    if(DB) cout<<"::: "<<cmd[index]<<endl;

                    if(cmd[index].isKey("cf")){
                        correctionFactor=cmd[index++][1];
                    continue;
                    }


                    if( cmd[index].isKey("difftime")){
                        diffTime=true;
                        index++;
                    continue;
                    }

                    if( cmd[index].isKey("norm")){
                        norm=true;
                        index++;
                    continue;
                    }

				#ifdef __linux__
                    if ( cmd[index].isKey("plot")){
                        plotPrm=cmd[index][1];

                        // if gnuplot statement is empty add a space
                        if(plotPrm.empty()) plotPrm+=" ";

                        index++;
                    continue;
                    }
				#endif


                    if(cmd[index]=="rangeN"){
                    const double from=std::stod(cmd[index][1]);
                    const int    N=   std::stoi(cmd[index][2]);
                    const double to=  std::stod(cmd[index][3]);
                    const double step=(to-from)/(N-1);                                                

                        range=cmd[index][1]+" "+std::to_string(step)+" "+cmd[index][3];

                        index++;
                    continue;
                    }

                    if(cmd[index]=="range"){
                        range=cmd[index][1]+" "+cmd[index][2];

                        if(cmd[index].numOfKeyValues()==4)
                            range+=" "+cmd[index][3];

                        index++;
                    continue;
                    }


                    if(cmd[index]=="save"){
                    string fileName(cmd[index++][1]);
                        Script::replaceVars(ptr_uvar,fileName);
                        fileNameOut.push_back(fileName);
                    continue;
                    }

                    #ifdef __linux__
                    if(cmd[index]=="saveFigure"){
                    string fileName(cmd[index][1]);
                            Script::replaceVars(ptr_uvar,fileName);
                            figFileName=fileName;

                            if(cmd[index].numOfKeyValues()==4){
                                figWidth=cmd[index][2];
                                figHeight=cmd[index][3];
                            }
                            else{
                                figWidth="800";
                                figHeight="600";
                            }

                            index++;
                    continue;
                    }
					#endif


                    if(cmd[index]=="threads"){
                        threads=cmd[index++][1];
                    continue;
                    }

                    if(cmd[index]=="wf"){
                        wf=(cmd[index++][1]);
                    continue;
                    }


                    if( cmd[index].isKey("printprm") || cmd[index].isKey("printPrm")){
                         printPrm=true;
                         index++;
                     continue;
                     }


                     cerr<<"Error: gr.cpp  unknown command "<< cmd[index]<<" line "<<__LINE__<<endl;
                 return false;
                }

return true;
}
//=============================================================================
void Cgr::calc()
{
            if(diff->t.empty){
                cerr<<"ERROR: diffraction data buffer is empty"<<endl;
            throw Cgr::ERR_EMPTY_DIFF;
            }


            cout<<" start G(r) computations"<<endl;
            if(printPrm)
                dispPrm();

            runCalc();

            if(norm)
                normGr();

            saveResults();

            #ifdef  __linux__
            if(!plotPrm.empty())     plotFigure();
            if(!figFileName.empty()) saveFigure(figFileName);
            #else
            infoMsg("plotting, saving to png not implemented for non-Linux OS");
            #endif

}
//=============================================================================

void Cgr::runCalc()
{
vector<string> ststst(split<string>(range," "));//start step stop
const double start=std::stod(ststst[0]);
const double step =std::stod(ststst[1]);
const double stop =(ststst.size()==2) ? diff->pdh->dataX[diff->pdh->dataSize()-1] : std::stod(ststst[2]);
const int  size=static_cast<int>(1+floor((stop-start)/step));
const double *const Q=diff->q.values;
const double *const S=diff->S.values;
const size_t diffDataSize=diff->t.size;


            if(!diffDataSize || size<0){
                cerr<<endl<<"ERROR: diffDataSize==0  or size<0 "<<endl;
            return;
            }


double *dataStmp=new double[diffDataSize];
const double Qmax=Q[diffDataSize-1];
CWindowFunction *cwf;


            if(wf.find("box")!=string::npos)
                cwf=new CWindowBoxCar();
            else{
                cwf=new CWindowLorch();
                cwf->setPrm(Qmax);
            }

            allocMem(size);

constexpr double k0=0.25*2.0*M_1_PI;

            omp_set_num_threads(std::stoi(threads));
            dataStmp[0]=0;



            for(size_t i=1;i<diffDataSize;i++)
                dataStmp[i]=(Q[i]+Q[i-1])*(Q[i]-Q[i-1])*(S[i]+S[i-1]-2)*(*cwf)(Q[i])*k0;

CProgress progress;
			
			progress.title=(std::string(" G(r) "));
			progress.start(size);

size_t j;
double x,sum;
//double vs,sumP,sumN;
const  double cf=std::stod(correctionFactor);
			
			
            for(int i=0;i<size ;i++,progress++){
                x=start+i*step;
                dataX[i]=x*cf;
                sum=0;
                //sumP=0;sumN=0;

                #pragma omp parallel for reduction(+:sum) private(j)
                for(j=1;j<diffDataSize;j++)
                    sum+=dataStmp[j]*sin(Q[j]*x);

                dataY[i]=sum;

                /*#pragma omp parallel for reduction(+:sumP,sumN) private(j)
                for(j=1;j<diffDataSize;j++){
                    vs=sin(Q[j]*x);

                    if(vs>=0) sumP+=dataStmp[j]*vs;
                    else      sumN+=dataStmp[j]*vs;


                }

                dataY[i]=sumP+sumN;*/

            }

            delete cwf;
            delete [] dataStmp;

}
//=============================================================================

void Cgr::allocMem(const size_t &size)
{
        dataX.allocMem(size);
        dataY.allocMem(size);
}
//=============================================================================

#ifdef __linux__

void Cgr::plotData(FILE *plotHandle)
{
const size_t dataSize=dataX.size;
const std::string sizeType{(sizeof(position)==sizeof(double))? "double" : "float"};
const std::string plotcmd("plot '-' binary format=\"%%2"+sizeType+"\" record="+
                            std::to_string(dataSize)+" endian=little "+plotPrm+" \n");

            if(DB){cout<<__FILE__<<":"<<__LINE__<<"  "<<plotcmd<<endl;}

            fprintf(plotHandle,plotcmd.c_str());

size_t i,j;
double *x=dataX.values;
double *y=dataY.values;

position *plotData=new position[2*dataSize];

            for(i=0,j=0;i<dataSize;i++,x++,y++){
                plotData[j++]=(*x);
                plotData[j++]=(*y);
            }


            fwrite(plotData,sizeof(position),2*dataSize,plotHandle);

            delete [] plotData;
}

//=============================================================================
void Cgr::plotFigure()
{

FILE *plotHandle = NULL;

            if(plotHandle == NULL){
              plotHandle = popen("gnuplot -persist 2>/dev/null", "w");
              assert(plotHandle != NULL);
            }

            fprintf(plotHandle,"set xlabel \"r\"\n");
            fprintf(plotHandle,"set ylabel \"G(r)\"\n");
            fprintf(plotHandle,"set xticks font \",16\"\n");
            fprintf(plotHandle,"set yticks font \",16\"\n");
            fprintf(plotHandle,"set grid\n");


            if(plotPrm.find("title")==string::npos)
                plotPrm+=" notitle";

            plotData(plotHandle);

            //fflush(plotHandle);
            pclose(plotHandle);


}
//=============================================================================
void Cgr::saveFigure(const string &fileName)
{
FILE *plotHandle = NULL;

            if(plotHandle == NULL){
              plotHandle = popen("gnuplot -persist 2>/dev/null", "w");
              assert(plotHandle != NULL);
            }

const std::string cmdSetTerm("set terminal png size "+figWidth+", "+figHeight+"\n");
const std::string cmdSetOutput("set output '"+fileName+"'\n");


            fprintf(plotHandle,cmdSetTerm.c_str());
            fprintf(plotHandle,cmdSetOutput.c_str());
            fprintf(plotHandle,"set xlabel \"r\" font\"Arial,18\"\n");
            fprintf(plotHandle,"set ylabel \"G(r)\" font\"Arial,18\"\n");
            fprintf(plotHandle,"set xticks font \",16\"\n");
            fprintf(plotHandle,"set yticks font \",16\"\n");
            fprintf(plotHandle,"set grid\n");

            plotPrm=" w l lw 2";

            if(plotPrm.find("title")==string::npos)
                plotPrm+=" notitle";


            plotData(plotHandle);
            pclose(plotHandle);
}

#endif
