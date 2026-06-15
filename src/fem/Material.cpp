#include "fem/Material.hpp"

namespace hvdc {
namespace fem {

void Material::compute_lame_from_E_nu() {
    lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    mu = E / (2.0 * (1.0 + nu));
    G = mu;
}

} // namespace fem
} // namespace hvdc
