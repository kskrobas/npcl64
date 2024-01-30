#ifndef CRDH_H
#define CRDH_H


#include "nanograin.h"
#include "stdatagen.h"

#include <ctime>


typedef StDataGeneric<position> dataRdh;
typedef const size_t csize;

class Crdh{
private:
        std::string fileName,fileNameIn;
        stdumap *ptr_uvar;

        position getBinWidth();

        void saveResults();
        void saveDatFile();
        void saveRdhlFile();
        void saveRdhsFile();
public:
        Crdh();

        enum Erdhstatus{OK,RUN,ERR_EMPTYGRAIN,ERR_FILEOPEN} status;

        void clearData();
        bool parseCommands(vcmdlist &cmd,size_t &index, stdumap *uvar__);
        void calc();

        NanoGrain::StNanoGrain *grain;

        std::string threads;
        std::string bin;

        dataRdh dataX;
        vector<dataRdh> dataYnn;

};






#endif // CRDF_H
