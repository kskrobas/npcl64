#include "info.h"
#include <iostream>
using namespace std;
void info()
{
            cerr<<" date: "<<__DATE__<<endl;
            cerr<<" author: Kazimierz Skrobas, kskrobas@unipress.waw.pl"<<endl;
            cerr<<" NPCLPATH: ";

            if(const char* env_p = std::getenv("NPCLPATH"))
               cerr<<env_p;        
}
