//-----------------------------------------------------------------------------
// Copyright 2016 seblemaguer
// Author: https://github.com/seblemaguer
// Last update: 2017/02/01
//
// Summary:
// The example synthesizes a .wav file by three files generated from
// "analysis.cpp".
//
// How to use:
// The format is shown in the line 129.
//-----------------------------------------------------------------------------
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <iostream>
#include <fstream>

#if (defined (__WIN32__) || defined (_WIN32)) && !defined (__MINGW32__)
#include <conio.h>
#include <windows.h>
#pragma comment(lib, "winmm.lib")
#pragma warning(disable : 4996)
#endif
#if (defined (__linux__) || defined(__CYGWIN__) || defined(__APPLE__))
#include <stdint.h>
#include <sys/time.h>
#endif

// For .wav input/output functions.
#include "audioio.h"

// WORLD core functions.
// Note: win.sln uses an option in Additional Include Directories.
// To compile the program, the option "-I $(SolutionDir)..\src" was set.
#include "world/d4c.h"
#include "world/dio.h"
#include "world/matlabfunctions.h"
#include "world/cheaptrick.h"
#include "world/stonemask.h"
#include "world/synthesis.h"
#include "world/synthesisrealtime.h"

#if (defined (__linux__) || defined(__CYGWIN__) || defined(__APPLE__))
// Linux porting section: implement timeGetTime() by gettimeofday(),
#ifndef DWORD
#define DWORD uint32_t
#endif
DWORD timeGetTime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    DWORD ret = static_cast<DWORD>(tv.tv_usec / 1000 + tv.tv_sec * 1000);
    return ret;
}
#endif

#define DEFAULT_FRAME_PERIOD 5.0

//-----------------------------------------------------------------------------
// struct for WORLD
// This struct is an option.
// Users are NOT forced to use this struct.
//-----------------------------------------------------------------------------
typedef struct {
    double frame_period;
    int fs;

    double *f0;
    double *time_axis;
    int f0_length;

    double **spectrogram;
    double **aperiodicity;
    int fft_size;
} WorldParameters;


namespace {


    std::ifstream::pos_type filesize(const char* filename)
    {
        std::ifstream in(filename, std::ifstream::binary | std::ifstream::ate);
        return in.tellg();
    }


    void WaveformSynthesis(WorldParameters *world_parameters, double *y)
    {
        DWORD elapsed_time;

        // compute nb of samples  stored info

        int y_length = static_cast<int>((world_parameters->f0_length - 1) *
                                        world_parameters->frame_period / 1000.0 * world_parameters->fs) + 1;
        // Synthesis by the aperiodicity
        elapsed_time = timeGetTime();
        Synthesis(world_parameters->f0, world_parameters->f0_length,
                  world_parameters->spectrogram, world_parameters->aperiodicity,
                  world_parameters->fft_size, world_parameters->frame_period,
                  world_parameters->fs, y_length, y);

        std::cout << "WORLD: " << timeGetTime() - elapsed_time << " [msec]" << std::endl;
    }

    void DestroyMemory(WorldParameters *world_parameters) {
        delete[] world_parameters->f0;
        for (int i = 0; i < world_parameters->f0_length; ++i) {
            delete[] world_parameters->spectrogram[i];
            delete[] world_parameters->aperiodicity[i];
        }
        delete[] world_parameters->spectrogram;
        delete[] world_parameters->aperiodicity;
    }

}  // namespace


/**
 * Main function
 *
 */
int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        fprintf(stderr, "%s <input_f0_file> <input_spectrum_file> <input_aperiodicity_file> <output_wav_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    // Define a default filled structures
    WorldParameters world_parameters;
    world_parameters.fs = 0; // this value will be filled when reading the spectrogram file
    world_parameters.f0_length = filesize(argv[1]) / sizeof(double);
    // The first bytes in the spectrogram file contains the sampling frequency and
    // the frame period, so we need to subtract those bytes from the spectrogram size
    size_t specSize = (size_t) filesize(argv[2]) - sizeof(world_parameters.fs) -
        sizeof(world_parameters.frame_period);
    // Be careful that .sp contains only first half of the spectrum
    world_parameters.fft_size = ( (specSize / (sizeof(double) * world_parameters.f0_length)) - 1 ) * 2;
    std::cout << "fft size = " << world_parameters.fft_size << std::endl;

    // 5.0 ms is the default value.
    // Generally, the inverse of the lowest F0 of speech is the best.
    // However, the more elapsed time is required.
    // This value will be replaced by the value stored in the spectrogram file
    world_parameters.frame_period = DEFAULT_FRAME_PERIOD;


    //---------------------------------------------------------------------------
    // Prepare memory
    //---------------------------------------------------------------------------
    world_parameters.f0 = new double[world_parameters.f0_length];

    world_parameters.spectrogram = new double*[world_parameters.f0_length];
    for (int i=0;i<world_parameters.f0_length; i++)
        world_parameters.spectrogram[i] = new double[world_parameters.fft_size / 2 + 1];

    world_parameters.aperiodicity = new double*[world_parameters.f0_length];
    for (int i=0;i<world_parameters.f0_length; i++)
        world_parameters.aperiodicity[i] = new double[world_parameters.fft_size / 2 + 1];

    //---------------------------------------------------------------------------
    // Loading
    //---------------------------------------------------------------------------
    // F0 loading
    std::ifstream is_f0(argv[1], std::ios::binary | std::ios::in);
    if ( !is_f0.is_open() )
        return false;

    is_f0.read(reinterpret_cast<char*>(world_parameters.f0),
               std::streamsize(world_parameters.f0_length*sizeof(double)));
    // for (int i=0; i<world_parameters.f0_length; i++)
    //     std::cout << world_parameters.f0[i] << std::endl;
    is_f0.close();

    // Spectrogram loading
    std::ifstream is_spectrogram(argv[2], std::ios::binary | std::ios::in);
    if ( !is_spectrogram.is_open() )
        return false;

    // read the sampling frequency
    is_spectrogram.read(reinterpret_cast<char*>(&world_parameters.fs),
            std::streamsize( sizeof(world_parameters.fs) ) );

    // read the frame period
    is_spectrogram.read(reinterpret_cast<char*>(&world_parameters.frame_period),
            std::streamsize( sizeof(world_parameters.frame_period) ) );

    // read the spectrogram data
    for (int i=0; i<world_parameters.f0_length; i++)
    {
        is_spectrogram.read(reinterpret_cast<char*>(world_parameters.spectrogram[i]),
                            std::streamsize((world_parameters.fft_size / 2 + 1)*sizeof(double)));
    }

    is_spectrogram.close();

    // Aperiodicity loading
    std::ifstream is_aperiodicity(argv[3], std::ios::binary | std::ios::in);
    if ( !is_aperiodicity.is_open() )
        return false;

    for (int i=0; i<world_parameters.f0_length; i++)
    {
        is_aperiodicity.read(reinterpret_cast<char*>(world_parameters.aperiodicity[i]),
                             std::streamsize((world_parameters.fft_size / 2 + 1)*sizeof(double)));
    }

    is_aperiodicity.close();


    //---------------------------------------------------------------------------
    // Synthesis
    //---------------------------------------------------------------------------

    int y_length = static_cast<int>((world_parameters.f0_length - 1) *
                                    world_parameters.frame_period / 1000.0 * world_parameters.fs) + 1;
    double *y = new double[y_length];
    for (int i = 0; i < y_length; ++i) y[i] = 0.0;
    WaveformSynthesis(&world_parameters, y);
    wavwrite(y, y_length, world_parameters.fs, 16, argv[4]);

    //---------------------------------------------------------------------------
    // Cleaning part
    //---------------------------------------------------------------------------
    delete[] y;
    DestroyMemory(&world_parameters);

    std::cout << "complete" << std::endl;
    return EXIT_SUCCESS;
}

/* synthesis_prog.cpp ends here */
