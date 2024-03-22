/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * cprogress.cpp
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



#include "cprogress.h"
#include <chrono>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>




using namespace std;
using namespace std::chrono;

volatile float CProgress::pos=0;
volatile bool  CProgress::active=false;
bool CProgress::printShown=false;

const float oneDay=3600*24;
const float oneHour=3600;
const float oneMin=60; 




//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

CProgress::CProgress()
{
//int 𝛼𝛽𝛾𝛿=1;
}


CProgress::CProgress(CProgress &&progress__): pthread(std::move(progress__.pthread))
{

}

//-----------------------------------------------------------------------------
CProgress &CProgress::operator=(CProgress && progress__)
{
    stop();

    pthread = std::move(progress__.pthread);
    return *this;
}


//=============================================================================



std::string CProgress::title(" progress: ");
size_t CProgress::stripeWidth=25;



//=============================================================================



#if TCOLOR

#if UC
#define shade0 "░"
#define shade1 "►"   //"▒"  //⚽
#define shade2 " "   //"█" //"▓"  🞊 ⬤ ◯
#define lbrac  "╟"
#define rbrac  "╢"
std::string CProgress::linehead("◉◎");
#define pon     "◉"
#define poff    "◎"
#else
#define shade0 "-"
#define shade1 "*"
#define shade2 "="
#define lbrac  "["
#define rbrac  "]"
std::string CProgress::linehead("/\\*");
#endif


void CProgress::print()
{
const int terminalCols=linehead.length()+title.length()+6+1+stripeWidth+5;
const float invSW=stripeWidth/100.0;
size_t spos,iter=0;


const time_t startTime = std::time(NULL);
size_t numOfTimeSamples=10;
double Tfin=1;
std::string strTime;

    system_clock::time_point today = system_clock::now();

    while (pos<100 && active ) {

        cout<<"\r"<<setw(terminalCols+strTime.length())<<setfill(' '); ///clear whole line
        cout<<"\r";

		
        //cout<<"\u001b[1;4240;315m";   // colors : 1;background/font
        //cout<<( (iter %2 ) ? "\u001b[1;440;315m" : "\u001b[1;4160;315m" );


        cout<<"\u001b[1;4240;315m";   // colors : 1;background/font
        #if UC
		cout<<( (iter %2 ) ? pon : poff );
		#else
		cout<<( (iter %2 ) ? linehead[1] : linehead[0] );
		#endif
			

        cout<<title<<std::fixed<<setprecision(1)<<setw(5)<<pos<<"%";
		
        cout<<"\033[0m \033[30m";
		
        if(pos<25)     cout<<"\u001b[48;5;1m";//red background
		else
            if(pos<75) cout<<"\u001b[48;5;3m";//dark yellow
            else       cout<<"\u001b[48;5;2m";//green
			
        cout<<lbrac;

        for(spos=0;spos<invSW*pos;spos++)   cout<<shade2;
        if(pos<100){              spos++;   cout<<shade1;
        for(   ;spos<=stripeWidth;spos++)   cout<<shade0; }

        cout<<rbrac<<"\033[0m ";

        iter++;


        if( !(iter%numOfTimeSamples) && pos>0 ){
        const time_t stopTime = std::time(NULL);
					Tfin=1000*std::difftime(stopTime,startTime)*100.0/pos;					// task duration in ms
					numOfTimeSamples+=10;
        }


        if(iter%4==0 && numOfTimeSamples>10){  // about 1s period
        const int msecToStop=(int)(Tfin*(1-pos*0.01));	// ms
        std::chrono::duration<int,ratio<1,1000> > timeToStop ( (int)(Tfin));
        system_clock::time_point timeStop = today + timeToStop;
        std::time_t tt (system_clock::to_time_t ( timeStop));
        std::string cstr(ctime(&tt));
                    cstr.pop_back();
        const float sec=msecToStop*0.001;
        const int  days=(int) (sec/oneDay);
        float      frac=(sec- days*oneDay);
        const int hours=(int) (frac/oneHour);
                    frac=(frac-hours*oneHour);
        const int mins=(int) (frac/oneMin);
        const int isec=(int) (frac-mins*oneMin);

        stringstream sTime;

                    sTime<<"  \u001b[38;5;14m";

                    if(days>0)  sTime<<setw(2)<<setfill('0')<< days<<":";
                    if(hours>0) sTime<<setw(2)<<setfill('0')<<hours<<":";
                    if(mins>0)  sTime<<setw(2)<<setfill('0')<< mins<<":";

                    sTime<<setw(2)<<setfill('0')<<isec<<"s # "<<cstr<<"\033[0m";
                    strTime=sTime.str();

                    //if(!sTime.good()) cerr<<" SSTream not good "<<__LINE__<<" "<<__FILE__;
        }

        cout<<strTime;
        cout.flush();

        this_thread::sleep_for (chrono::milliseconds(250));
    }

	if(pos>100) pos=100;

    cout<<"\r"<<setw(terminalCols+strTime.length())<<setfill(' '); ///clear whole line
    cout<<"\rX"<<title<<std::fixed<<setprecision(1)<<setw(5)<<pos<<"%"<<"\033[0m";
	printShown=true;
}


#else

std::string CProgress::linehead("/\\*");

void CProgress::print()
{
const int terminalCols=linehead.length()+title.length()+6+1+stripeWidth+5;
const float invSW=stripeWidth/100.0;
size_t spos,iter=0;


const time_t startTime = std::time(NULL);
size_t numOfTimeSamples=10;
double Tfin=1;
std::string strTime;

    system_clock::time_point today = system_clock::now();

    while (pos<100 && active ) {

        cout<<"\r"<<setw(terminalCols+strTime.length())<<setfill(' '); ///clear whole line
		
        cout<<"\r";

        cout<<( (iter %2 ) ? linehead[1] : linehead[0] );

        cout<<title<<std::fixed<<setprecision(1)<<setw(5)<<pos<<"%";

        cout<<" [";

        for(spos=0;spos<invSW*pos;spos++)   cout<<"=";
        if(pos<100){              spos++;   cout<<"*"; 
        for(   ;spos<=stripeWidth;spos++)   cout<<"-"; }

        cout<<"] ";

        iter++;


        if(iter==numOfTimeSamples && pos>0 ){
        const time_t stopTime = std::time(NULL);
					Tfin=1000*std::difftime(stopTime,startTime)*100.0/pos;					// task duration in ms
					numOfTimeSamples+=20;
        }

        if(numOfTimeSamples>10 && iter%4==0){  // about 1s period
        const int msecToStop=(int)(Tfin*(1-pos*0.01));	// ms
        std::chrono::duration<int,ratio<1,1000> > timeToStop ( (int)(Tfin));
        system_clock::time_point timeStop = today + timeToStop;
        std::time_t tt (system_clock::to_time_t ( timeStop));
        std::string cstr(ctime(&tt));
                    cstr.pop_back();
        const float sec=msecToStop*0.001;
        const int  days=(int) (sec/oneDay);
        float      frac=(sec- days*oneDay);
        const int hours=(int) (frac/oneHour);
                    frac=(frac-hours*oneHour);
        const int mins=(int) (frac/oneMin);
        const int isec=(int) (frac-mins*oneMin);

        stringstream sTime;

                    if(days>0)  sTime<<setw(2)<<setfill('0')<< days<<":";
                    if(hours>0) sTime<<setw(2)<<setfill('0')<<hours<<":";
                    if(mins>0)  sTime<<setw(2)<<setfill('0')<< mins<<":";

                    sTime<<setw(2)<<setfill('0')<<isec<<"s # "<<cstr;
                    strTime=sTime.str();
        }

        cout<<strTime;
        cout.flush();

        this_thread::sleep_for (chrono::milliseconds(250));
    }


    cout<<"\r"<<setw(terminalCols+strTime.length())<<setfill(' '); ///clear whole line
    cout<<"\rX"<<title<<std::fixed<<setprecision(1)<<setw(5)<<pos<<"%";
    cout.flush();
    printShown=true;

}

#endif


//-----------------------------------------------------------------------------
void CProgress::start(const size_t numOfsteps)
{
    if(numOfsteps<1){
        cerr<<"ERROR: number of steps must be greater than 0";
    return;
    }


    istep=100.0/(numOfsteps);
    pos=0;
    active=true;
    printShown=false;

    pthread=std::thread(print);
}

//-----------------------------------------------------------------------------
void CProgress::stop()
{

    pos=100;
    active=false;

    //cout<<"progress stop"<<endl;

    if(pthread.joinable())
        pthread.join();


    //if(!printShown){
    	cout<<"\r"<<setw(80)<<setfill(' '); ///clear whole line
    	cout<<"\rX"<<title<<std::fixed<<setprecision(1)<<setw(5)<<100<<"%";
    	cout.flush();
    //}
    //else
        //printShown=true;
}
//-----------------------------------------------------------------------------
void CProgress::next()
{
    pos+=istep;
}

void CProgress::operator++()
{
    pos +=istep;
}

void CProgress::operator ++(int)
{
    pos +=istep;
}



//=============================================================================
CProgress::~CProgress()
{
    stop();
}
