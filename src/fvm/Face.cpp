#include "fvm/Face.hpp"

namespace hvdc {
namespace fvm {

Face::Face(Index id, Index num_nodes)
    : id_(id), node_ids_(num_nodes, -1)
{
}

} // namespace fvm
} // namespace hvdc
