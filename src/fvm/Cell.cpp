#include "fvm/Cell.hpp"

namespace hvdc {
namespace fvm {

Cell::Cell(Index id, CellType type)
    : id_(id), type_(type)
{
    switch (type) {
        case CellType::Hexa:
            node_ids_.resize(8);
            face_ids_.resize(6);
            neighbor_ids_.resize(6, -1);
            break;
        case CellType::Tetra:
            node_ids_.resize(4);
            face_ids_.resize(4);
            neighbor_ids_.resize(4, -1);
            break;
        case CellType::Prism:
            node_ids_.resize(6);
            face_ids_.resize(5);
            neighbor_ids_.resize(5, -1);
            break;
        case CellType::Pyramid:
            node_ids_.resize(5);
            face_ids_.resize(5);
            neighbor_ids_.resize(5, -1);
            break;
        default:
            break;
    }
}

} // namespace fvm
} // namespace hvdc
