#include "crdh.h"
#include<fstream>
#include<iostream>
#include<iomanip>
#include"createdir.h"
#include"colormsg.h"
#include"cprogress.h"

inline double sqrd (const position &x)  {return x*x;}

using namespace std;
//===================================================================================

Crdh::Crdh()
{

}
//===================================================================================
void Crdh::clearData()
{
    threads="1";
    fileName="";
    saveopt.clear();

    range.clear();
    fileNameIn.clear();
    dataX.freeMem();
    dataYnn.clear();
}
//===================================================================================
position Crdh::getBinWidth()
{
const position wbin=std::stod(bin);
return (wbin>1)? 1.0/wbin : wbin;
}
//===================================================================================
void Crdh::saveResults()
{

        if(!createDirsIfDontExist(fileName)){
            errMsg("couldn't create nested directories for "+fileName);
        throw Crdh::ERR_FILEOPEN;
        }

        if(fileName.rfind(".dat")!=string::npos){
            saveDatFile();
        return;
        }

        if(fileName.rfind(".rdhl")!=string::npos){
            saveRdhlFile();
        return;
        }

        if(fileName.rfind(".rdhs")!=string::npos){
            saveRdhsFile();
        return;
        }


        cerr<<"Warning:  unknown file extension, RDH results not saved"<<endl;
}


//===================================================================================


bool Crdh::parseCommands(vcmdlist &cmd, size_t &index, stdumap *uvar__)
{
const str send("end");

           clearData();
           ptr_uvar=uvar__;

           while(cmd[index]!=send){

               if(cmd[index]=="bin"){
                   bin=cmd[index++][1];
               continue;
               }

               if(cmd[index]=="range"){
                    range=cmd[index][1]+" "+cmd[index][2]+" "+cmd[index][3];
                    index++;
               continue;
               }

               if(cmd[index]=="save"){
                   fileName=cmd[index++][1];
                   Script::replaceVars(ptr_uvar,fileName);
               continue;
               }

               if(cmd[index]=="saveopt"){
                   saveopt=cmd[index++][1];
               continue;
               }

               if(cmd[index]=="threads"){
                   threads=cmd[index++][1];
               continue;
               }

               cerr<<"Error: cprh.cpp  unknown command "<< cmd[index]<<" line "<<__LINE__<<endl;
           return false;
           }

return true;
}
//===================================================================================
void Crdh::calc()
{
csize numOfAtoms=grain->atoms.size();
csize numOfRdhAtoms=grain->uc.rdhAtoms;
auto atoms=grain->atoms.data();

vector<string> trange{split<string>(range," ")};
cdouble rmin =std::stod(trange[0]);
cdouble rstep=std::stod(trange[1]);
cdouble rmax =std::stod(trange[2]);

                if(rmax-rmin<0 || (rstep>rmax-rmin) ){
                    errMsg("wrong values of range instruction, usage: min step max");
                throw Crdh::ERR_RANGE;
                }

csize size=static_cast<size_t>(std::ceil((rmax-rmin)/rstep/0.95));
cdouble ibin=1.0/rstep;


cdouble rmin2=rmin*rmin;
cdouble rmax2=rmax*rmax;
auto calcR2=[](NanoGrain::StAtom *A, NanoGrain::StAtom *B){
        return sqrd(A->x-B->x) + sqrd(A->y-B->y) + sqrd(A->z-B->z) ; };


                dataYnn.resize(numOfRdhAtoms);
                for(size_t i=0;i<numOfRdhAtoms;i++)
                    dataYnn[i].allocMem(size);

                dataX.allocMem(size);
                for(size_t i=0;i<size;i++){
                    dataX[i]=rmin+i*rstep;
                    for(size_t k=0;k<numOfRdhAtoms;k++)
                        dataYnn[k][i]=0;
                }

                /// optimalization on
                //for(size_t i=0;i<numOfAtoms;i++)
                //    grain->atoms[i]*=ibin;


size_t rdhAtom=0;
size_t k,bin;
CProgress progress;
double r,dr2;

                progress.title=(string(" RDH "));
                progress.start(numOfAtoms);

                for(size_t ia=0;ia<numOfAtoms && rdhAtom<numOfRdhAtoms;ia++,progress++){
                    if(atoms[ia].rdh){

                        for(k=0;k<ia;k++){
                            dr2=calcR2(&atoms[ia],&atoms[k]);

                            if(rmin2<=dr2 && dr2<=rmax2){
                                r=std::sqrt(dr2)-rmin;
                                bin=static_cast<size_t>(r/rstep);
                                    dataYnn[rdhAtom][bin]++;
                            }
                        }

                        for(k=ia+1;k<numOfAtoms;k++){
                            dr2=calcR2(&atoms[ia],&atoms[k]);

                            if(rmin2<=dr2 && dr2<=rmax2){
                                r=std::sqrt(dr2)-rmin;
                                bin=static_cast<size_t>(r/rstep);
                                    dataYnn[rdhAtom][bin]++;
                            }
                        }
                        rdhAtom++;
                    }
                }

        /// optimalization off
        //for(size_t i=0;i<numOfAtoms;i++){
        //    grain->atoms[i]*=wbin;
        //}

        if(!fileName.empty())
            saveResults();

}
//===================================================================================
void Crdh::saveDatFile()
{

}
//===================================================================================

void Crdh::saveRdhlFile()
{
fstream fout(fileName,ios::out);

        if(!fout){
            cerr<<"ERROR: couldn't open file for saving, PDH results will be lost"<<endl;
            //status=Cpdh::ERR_FILEOPEN;
        return;
        }


const size_t dataSize=dataX.size;
csize numOfAtomTypes= grain->atomTypes.size();
csize numOfAtoms=   grain->atoms.size();
csize mrdhSize=grain->uc.rdhAtoms;
std::time_t datetime = std::time(nullptr);
size_t i;

            //--------------------------- HEADER ----------------------------//
            fout<<"#ver: 0"<<endl;
       //     fout<<"#title: "<<title<<endl;
            fout<<"#sizeRC: "<<dataSize<<"\t"<<(1+mrdhSize)<<endl;
            fout<<"#date: "<<std::asctime(std::localtime(&datetime));
            fout<<"#atomTypes: "<<numOfAtomTypes<<"\t";

            for(i=0;i<numOfAtomTypes;i++)
                fout<<"\t"<<grain->atomTypes[i].name;

            fout<<endl;

                fout<<"#atomsNumber: "<<numOfAtoms<<endl;
                fout<<"#atomsXYZ: "<<((grain->fileNameIn.empty()) ? "npcl" : grain->fileNameIn)<<endl;
                fout<<"#binWidth: "<<bin<<endl;
                //fout<<"#comment: "<<comment<<endl;
                fout<<"#";

//csize nprec=11;
csize colwh=12;


            /// wypisuje nazwy column
            fout<<setw(colwh)<<'X';
            for(auto ucAtom: grain->uc.atoms)
                if(ucAtom.rdh)
                    fout<<" "<<setw(colwh)<<ucAtom.name;

            fout<<endl;

            //--------------------------------------------------------------------//


            for(i=0;i<dataSize;i++){
                fout<<" "<<setw(colwh)<<dataX[i];

                for(dataRdh & rdh : dataYnn)
                  fout<<" "<<setw(colwh)<<rdh[i];

                fout<<endl;
            }//for's end


            fout.close();
}
//===================================================================================
size_t sumPartsSave(fstream &fout,csize i, vector<dataRdh> &dataYnn)
{
//csize nprec=11;
csize colwh=12;
size_t sum=0;

        for(dataRdh & rdh : dataYnn){
          sum+=rdh[i];
          fout<<" "<<setw(colwh)<<rdh[i];
        }
return sum;
}
//===================================================================================
size_t sumParts(fstream &fout,csize i, vector<dataRdh> &dataYnn)
{
size_t sum=0;
        for(dataRdh & rdh : dataYnn)
          sum+=rdh[i];
return sum;
}
//===================================================================================
void Crdh::saveRdhsFile()
{
fstream fout(fileName,ios::out);

        if(!fout){
            cerr<<"ERROR: couldn't open file for saving, PDH results will be lost"<<endl;
            //status=Cpdh::ERR_FILEOPEN;
        return;
        }


const size_t dataSize=dataX.size;
csize numOfAtomTypes= grain->atomTypes.size();
csize numOfAtoms=   grain->atoms.size();
csize mrdhSize=grain->uc.rdhAtoms;
std::time_t datetime = std::time(nullptr);
size_t i;
std::streampos sizeRCpos;
size_t nonZeroBins=0;
bool nzb;

bool soptParts,soptTot;



            if(saveopt.empty()){
               soptParts=soptTot=true;
            }
            else{
            vector<string> sopt{split<string>(saveopt," ")};
                soptParts=soptTot=false;

                for(auto &opt: sopt){
                    if(opt=="parts") soptParts=true;
                    if(opt=="tot")   soptTot  =true;
                }
            }

            //--------------------------- HEADER ----------------------------//
            fout<<"#ver: 0"<<endl;
       //     fout<<"#title: "<<title<<endl;
           // fout<<"#sizeRC: "<<dataSize<<"\t"<<(1+mrdhSize)<<endl;
            fout<<"#sizeRC: ";
            sizeRCpos=fout.tellp();
            fout<<setw(20)<<' '<<endl;

            fout<<"#date: "<<std::asctime(std::localtime(&datetime));
            fout<<"#atomTypes: "<<numOfAtomTypes<<"\t";

            for(i=0;i<numOfAtomTypes;i++)
                fout<<"\t"<<grain->atomTypes[i].name;

            fout<<endl;

                fout<<"#atomsNumber: "<<numOfAtoms<<endl;
                fout<<"#atomsXYZ: "<<((grain->fileNameIn.empty()) ? "npcl" : grain->fileNameIn)<<endl;
                fout<<"#binWidth: "<<bin<<endl;
                //fout<<"#comment: "<<comment<<endl;
                fout<<"#";


//csize numOfRdhAtoms=grain->uc.rdhAtoms;
//csize nprec=11;
csize colwh=12;



            /// wypisuje nazwy column
            fout<<setw(colwh)<<'X';

            if(soptParts){
                for(auto ucAtom: grain->uc.atoms)
                    if(ucAtom.rdh)
                        fout<<" "<<setw(colwh)<<ucAtom.name;
            }

            if(soptTot){
                fout<<" "<<setw(colwh)<<"Total";
            }

            fout<<endl;

            //--------------------------------------------------------------------//




            if(soptTot){
            auto fs=(soptParts) ? &sumPartsSave: &sumParts;
            size_t sum;

                for(i=0;i<dataSize;i++){

                    nzb=false;
                    for(dataRdh &rdh : dataYnn){
                        if(rdh[i]>0) {
                            nzb=true;
                        break;
                        }
                    }

                    if(nzb){
                        nonZeroBins++;

                        fout<<" "<<setw(colwh)<<dataX[i];
                        sum=fs(fout,i,dataYnn);
                        fout<<" "<<setw(colwh)<<sum;

                        fout<<endl;
                    }
                }//for's end

            }
            else{
                for(i=0;i<dataSize;i++){

                    nzb=false;
                    for(dataRdh &rdh : dataYnn){
                        if(rdh[i]>0) {
                            nzb=true;
                        break;
                        }
                    }

                    if(nzb){
                        nonZeroBins++;
                        fout<<" "<<setw(colwh)<<dataX[i];
                        for(dataRdh & rdh : dataYnn)
                          fout<<" "<<setw(colwh)<<rdh[i];

                        fout<<endl;
                    }
                }//for's end
            }

            //--------------------------------------------------------------------//
            fout.seekp(sizeRCpos,ios_base::beg);

            fout<<nonZeroBins<<"\t";
                if(soptParts && !soptTot) fout <<(1+mrdhSize);
                else
                    if(soptParts && soptTot) fout<<(2+mrdhSize);
                    else
                        fout<<2;

            fout.close();
}
