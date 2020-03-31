// PluginTrax.cpp
// Gianluca Elia (elgiano@gmail.com)

#include "SC_PlugIn.hpp"
#include "Trax.hpp"

static InterfaceTable* ft;

namespace Trax {

Trax::Trax() {
    mCalcFunc = make_calc_function<Trax, &Trax::next>();
    next(1);
}

void Trax::next(int nSamples) {
    const float* input = in(0);
    const float* gain = in(0);
    float* outbuf = out(0);

    // simple gain function
    for (int i = 0; i < nSamples; ++i) {
        outbuf[i] = input[i] * gain[i];
    }
}

} // namespace Trax

PluginLoad(TraxUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<Trax::Trax>(ft, "Trax");
}
