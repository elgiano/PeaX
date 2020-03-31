// PluginTrax.hpp
// Gianluca Elia (elgiano@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"

namespace Trax {

class Trax : public SCUnit {
public:
    Trax();

    // Destructor
    // ~Trax();

private:
    // Calc function
    void next(int nSamples);

    // Member variables
};

} // namespace Trax
