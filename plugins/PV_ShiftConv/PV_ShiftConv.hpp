// PluginTrax.hpp
// Gianluca Elia (elgiano@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"

namespace Trax {

class PV_ShiftConv : public SCUnit {
public:
    PV_ShiftConv();

    ~PV_ShiftConv(){ if(m_tempbuf!=0){ RTFree(mWorld,m_tempbuf); RTFree(mWorld, m_binCounter); }};

private:
    // Calc function
    void next(int nSamples);
    void shift(SndBuf* buf1, SndBuf* buf2);
    void shiftStretch(SndBuf* buf1, SndBuf* buf2);
    void shift_complex(SndBuf* buf1, SndBuf* buf2);
    void shiftStretch_complex(SndBuf* buf1, SndBuf* buf2);

    float* m_tempbuf;
    float* m_binCounter;
    int m_numbins;

    float m_ampmult;
    float m_rAmpmult;

    int fromBinA, toBinA, fromBinB, toBinB;
    float thr;
    int f0Bin, interp;

    void (PV_ShiftConv::*combineFunc)(SndBuf*, SndBuf*);

};

} // namespace Trax
