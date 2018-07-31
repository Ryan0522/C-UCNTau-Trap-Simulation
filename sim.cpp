#include <vector>
#include <stdio.h>
#include "inc/track_gen.hpp"
#include "inc/trackUCN.hpp"
#include "inc/fields_nate.h"
#include "inc/lyap.hpp"
#include <cmath>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <cstring>

#include "setup.h"

//#include "inc/geometry.hpp"

extern "C" {
    #include "inc/xorshift.h"
}

/*typedef struct fixedResult {
    double energy;
    double theta;
    double t;
    double ePerp;
    double x;
    double y;
    double z;
    double zOff;
    int nHit;
    int nHitHouseLow;
    int nHitHouseHigh;
    double eStart;
    double deathTime;
} fixedResult;*/

void writeFixedRes(std::ofstream &binfile, fixedResult res) {
    const size_t buff_len = sizeof(unsigned int) + 9*sizeof(double) + 3*sizeof(int) + 2*sizeof(double) + sizeof(unsigned int);
    char buf[buff_len];
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! file closed\n");
        return;
    }
    *((unsigned int *)(&buf[0])) = buff_len - 2*sizeof(unsigned int);
    *((double *)(&buf[0] + sizeof(unsigned int))) = res.energy;
    *((double *)(&buf[0] + sizeof(unsigned int) + 1*sizeof(double))) = res.theta;
    *((double *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double))) = res.t;
    *((double *)(&buf[0] + sizeof(unsigned int) + 3*sizeof(double))) = res.settlingT;
    *((double *)(&buf[0] + sizeof(unsigned int) + 4*sizeof(double))) = res.ePerp;
    *((double *)(&buf[0] + sizeof(unsigned int) + 5*sizeof(double))) = res.x;
    *((double *)(&buf[0] + sizeof(unsigned int) + 6*sizeof(double))) = res.y;
    *((double *)(&buf[0] + sizeof(unsigned int) + 7*sizeof(double))) = res.z;
    *((double *)(&buf[0] + sizeof(unsigned int) + 8*sizeof(double))) = res.zOff;
    *((int *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double))) = res.nHit;
    *((int *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 1*sizeof(int))) = res.nHitHouseLow;
    *((int *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 2*sizeof(int))) = res.nHitHouseHigh;
    *((double *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 3*sizeof(int))) = res.eStart;
    *((double *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 3*sizeof(int) + 1*sizeof(double))) = res.deathTime;
    *((unsigned int *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 3*sizeof(int) + 2*sizeof(double))) = buff_len - 2*sizeof(unsigned int);
    binfile.write(buf, buff_len);
}

void writeNoabsRes(std::ofstream &binfile, noabsResult res) {
    if(res.energy < 0.0) {
        return;
    }
    const size_t buff_len = sizeof(unsigned int) + 2*sizeof(double) + 3*NRECORDS*sizeof(float) + sizeof(unsigned int);
    char buf[buff_len];
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! file closed\n");
        return;
    }
    *((unsigned int *)(&buf[0])) = buff_len - 2*sizeof(unsigned int);
    *((double *)(&buf[0] + sizeof(unsigned int))) = res.energy;
    *((double *)(&buf[0] + sizeof(unsigned int) + 1*sizeof(double))) = res.theta;
    std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double)), (void *)&(res.times[0]), NRECORDS*sizeof(float));
    std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + NRECORDS*sizeof(float)), (void *)&(res.ePerps[0]), NRECORDS*sizeof(float));
    std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 2*NRECORDS*sizeof(float)), (void *)&(res.zetas[0]), NRECORDS*sizeof(float));
    *((unsigned int *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 3*NRECORDS*sizeof(float))) = buff_len - 2*sizeof(unsigned int);
    binfile.write(buf, buff_len);
}

void writeCleanRes(std::ofstream &binfile, cleanResult res) {
    const size_t buff_len = sizeof(unsigned int) + 6*sizeof(float) + 1*sizeof(int) + sizeof(unsigned int);
    char buf[buff_len];
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! file closed\n");
        return;
    }
    *((unsigned int *)(&buf[0])) = buff_len - 2*sizeof(unsigned int);
    *((float *)(&buf[0] + sizeof(unsigned int))) = res.energy;
    *((float *)(&buf[0] + sizeof(unsigned int) + 1*sizeof(float))) = res.theta;
    *((float *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(float))) = res.t;
    *((float *)(&buf[0] + sizeof(unsigned int) + 3*sizeof(float))) = res.x;
    *((float *)(&buf[0] + sizeof(unsigned int) + 4*sizeof(float))) = res.y;
    *((float *)(&buf[0] + sizeof(unsigned int) + 5*sizeof(float))) = res.z;
    *((int *)(&buf[0] + sizeof(unsigned int) + 6*sizeof(float))) = res.code;
    *((unsigned int *)(&buf[0] + sizeof(unsigned int) + 6*sizeof(float) + 1*sizeof(int))) = buff_len - 2*sizeof(unsigned int);
    binfile.write(buf, buff_len);
}

void writeLyapRes(std::ofstream &binfile, lyapResult res) {
    const size_t buff_len = sizeof(unsigned int) + 4*sizeof(float) + sizeof(unsigned int);
    char buf[buff_len];
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! file closed\n");
        return;
    }
    *((unsigned int *)(&buf[0])) = buff_len - 2*sizeof(unsigned int);
    *((float *)(&buf[0] + sizeof(unsigned int))) = res.eStart;
    *((float *)(&buf[0] + sizeof(unsigned int) + 1*sizeof(float))) = res.eEnd;
    *((float *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(float))) = res.theta;
    *((float *)(&buf[0] + sizeof(unsigned int) + 3*sizeof(float))) = res.lce;
    *((unsigned int *)(&buf[0] + sizeof(unsigned int) + 4*sizeof(float))) = buff_len - 2*sizeof(unsigned int);
    binfile.write(buf, buff_len);
}

trace readTrace(const char *xfile, const char *yfile, const char *zfile) {
    trace t;
    t.x = NULL;
    t.y = NULL;
    t.z = NULL;
    t.num = -1;
    
    const size_t buff_len = 1*8;
    char* buf = new char[buff_len];
    
    std::ifstream binfileX(xfile, std::ios::in | std::ios::binary);
    std::ifstream binfileY(yfile, std::ios::in | std::ios::binary);
    std::ifstream binfileZ(zfile, std::ios::in | std::ios::binary);
    if(!binfileX.is_open() || !binfileY.is_open() || !binfileZ.is_open()) {
        fprintf(stderr, "Error! Could not open files!\n");
        return t;
    }
    
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    
    while(!binfileX.eof()) {
        binfileX.read(buf, buff_len);
        if(binfileX.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
            break;
        }
        x.push_back(*(double *)&buf[0]);
    }
    binfileX.close();
    
    while(!binfileY.eof()) {
        binfileY.read(buf, buff_len);
        if(binfileY.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
            break;
        }
        y.push_back(*(double *)&buf[0]);
    }
    binfileY.close();
    
    while(!binfileZ.eof()) {
        binfileZ.read(buf, buff_len);
        if(binfileZ.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
            break;
        }
        z.push_back(*(double *)&buf[0]);
    }
    binfileZ.close();
    
    if(z.size() != x.size() || z.size() != y.size()) {
        fprintf(stderr, "Error! Sample length mismatch!\n");
        return t;
    }
    
    t.x = new double[x.size()];
    t.y = new double[y.size()];
    t.z = new double[z.size()];
    
    for(int i = 0; i < x.size(); i++) {
        t.x[i] = x[i];
        t.y[i] = y[i];
        t.z[i] = z[i];
    }
    
    t.num = x.size();
    
    return t;
}

int main(int argc, char** argv) {
    int ierr = MPI_Init(&argc, &argv);
    int nproc;
    int rank;
    
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    
    
    
    int c;
    
    char fName[1024];
    fName[0] = '\0';
    double dt = 0.0;
    double nTraj = 0;

    while (1) {
        static struct option long_options[] = {
            {"file", required_argument, 0, 'f'},
            {"dt", required_argument, 0, 'd'},
            {"ntraj", required_argument, 0, 'n'},
            {0, 0, 0, 0},
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv, "f:d:n:", long_options, &option_index);

        /* Detect the end of the options. */
        if(c == -1) {
            break;
        }

        switch(c) {
            case 'f':
                strncpy(fName, optarg, 1024-1);
                fName[1024-1] = '\0';
                break;
            case 'd':
                dt = atof(optarg);
                break;
            case 'n':
                nTraj = atoi(optarg);
                break;
            case '?':
                /* getopt_long already printed an error message. */
                break;
            default:
                exit(1);
        }
    }
    
    if(fName[0] == '\0' || dt == 0 || nTraj == 0) {
        fprintf(stderr, "Error! Usage: ./arrival_time_fixed_eff --file=fName --dt=timestep --ntraj=N\n");
        exit(1);
    }
    
//    trace tr = readTrace("./xvals.bin", "./yvals.bin", "./zvals.bin");
    trace tr = readTrace(XFNAME, YFNAME, ZFNAME);
    if(tr.x == NULL || tr.y == NULL || tr.z == NULL) {
        fprintf(stderr, "Bad trace\n");
        return 2;
    }
    printf("num trace bins: %d\n", tr.num);
    
    char fNameRank[1024];
    snprintf(fNameRank, 1024, "%s%d", fName, rank);
    std::ofstream binfile(fNameRank, std::ios::out | std::ios::binary);
    
    
    initxorshift(INITBLOCK);
    
//    printf("%d - %d\n", rank*(nTraj/nproc), (rank+1)*(nTraj/nproc));
    
    std::vector<std::vector<double>> traj;
    traj.reserve(nTraj/nproc);
    for(int i = 0; i < nTraj; i++) {
        if(i >= rank*(nTraj/nproc) && i < (rank+1)*(nTraj/nproc)) {
            traj.push_back(TRACKGENERATOR(tr));
        }
        else {
            TRACKGENERATOR(tr);
        }
    }
    
    for(int i = 0; i < rank; i++) {
        jump();
    }
    
    if(TRACKANDPRINT) {
        trackAndPrint(traj.back(), dt, tr);
        binfile.close();
        ierr = MPI_Finalize();
        return 0;
    }
        
    for(auto it = traj.begin(); it < traj.end(); it++) {
//        double lyap = calcLyap(*it, dt, tr, 0.0);
//        printf("%f\n", lyap);
        auto res = TRACKER(*it, dt, tr);
//        writeFixedRes(binfile, res);
        WRITER(binfile, res);
    }
    
    binfile.close();
    ierr = MPI_Finalize();
    
    return 0;
}
