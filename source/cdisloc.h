#ifndef CDISLOC_H
#define CDISLOC_H

#include <vector>
#include "nanograin.h"


using namespace std;

class Cdisloc{
private:


    void resetPrms();
public:

    Cdisloc();

vector<NanoGrain::StAtomType> atomTypes;
string axis,position,rangeR,rangeA;
string mindist,scatter;

string fileNameIn,fileNameHeader;
vector<string> fileNameOut;


    bool parseCommands(vcmdlist &cmd, size_t &index, stdumap *uvar__);
    void calc();






};




#endif // CDISLOC_H
