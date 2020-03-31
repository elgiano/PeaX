// PluginTrax.hpp
// Gianluca Elia (elgiano@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"

namespace Trax {

class HPS : public SCUnit {
public:
    HPS();
    // Destructor
    ~HPS();

private:
    // Calc function
    void next(int nSamples);
    double getTrueFreq(int bin_i, float phaseDiff);
    int hps(SCPolar* spectrum, int numbins, int minbin, int maxbin, int harmonics, float thr, bool correction);

    // Member variables
    bool m_init;
    double m_freqPerBin;
    double m_rFreqPerBin;
    double m_expct;
    float  m_oversampling;
    float outval;
    float* m_prevPhases;
};

} // namespace Trax
