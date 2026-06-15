#include "common/Types.hpp"

namespace hvdc {

Real const* gravity_direction() {
    static const Real g_dir[3] = {0.0, 0.0, -1.0};
    return g_dir;
}

} // namespace hvdc
