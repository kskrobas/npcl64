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
        std::string saveopt;
        std::string bin,range;
        stdumap *ptr_uvar;

        dataRdh dataX;
        vector<dataRdh> dataYnn;

        position getBinWidth();
        void clearData();

        void saveResults();
        void saveDatFile();
        void saveRdhlFile();
        void saveRdhsFile();
public:
        Crdh();

        enum Erdhstatus{OK,RUN,ERR_EMPTYGRAIN,ERR_FILEOPEN,ERR_RANGE} status;

        bool parseCommands(vcmdlist &cmd,size_t &index, stdumap *uvar__);
        void calc();

        NanoGrain::StNanoGrain *grain;
        std::string threads;
};






#endif // CRDF_H
