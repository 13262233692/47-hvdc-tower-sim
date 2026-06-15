#include "fem/Node.hpp"
#include "common/Vector.hpp"

namespace hvdc {
namespace fem {

Node::Node(Index id, Real x, Real y, Real z)
    : id_(id), coords_{x, y, z} {
    bc_types_.fill(BoundaryType::Free);
    bc_values_.fill(0.0);
}

Node::Node(Index id, const Vec3& coords)
    : id_(id), coords_(coords) {
    bc_types_.fill(BoundaryType::Free);
    bc_values_.fill(0.0);
}

Vec3 Node::displaced_coords(const Vec3& disp) const {
    return {coords_[0] + disp[0], coords_[1] + disp[1], coords_[2] + disp[2]};
}

Vec3 Node::displaced_coords(const Vector& displacement) const {
    Vec3 c = coords_;
    for (Index i = 0; i < std::min<Index>(3, ndofs_); ++i) {
        Index gdof = dof_start_ + i;
        if (gdof < displacement.size()) c[i] += displacement[gdof];
    }
    return c;
}

} // namespace fem
} // namespace hvdc
