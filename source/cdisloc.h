#ifndef CDISLOC_H
#define CDISLOC_H

#include <vector>
#include "nanograin.h"


using namespace std;

typedef const string cstring;

class Cdisloc{
private:


    void clearData();

    void insertLoop();
    void rotateLoop();
    void validateDisloc(NanoGrain::vatoms &dislocAtoms);

public:
    enum Edisstatus{OK, ERR_NOFPARAMS};

vector<NanoGrain::StAtomType> atomTypes;
string axis,axispos,rangeR,rangeA;
string angle,projh;
string mindist,scatter,mode;
string saveopt;

string fileNameIn,fileNameHeader;
vector<string> fileNameOut;

NanoGrain::StNanoGrain *grain;
stdumap *ptr_uvar;




    Cdisloc();

    int findAtomName(cstring &aname__);
    void addAtomName(cstring &aname__);

    bool parseCommands(vcmdlist &cmd, size_t &index, stdumap *uvar__);
    void calc();






};




#endif // CDISLOC_H
