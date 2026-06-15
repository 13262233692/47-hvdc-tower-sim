#ifndef HVDC_FVM_FIELD_HPP
#define HVDC_FVM_FIELD_HPP

#include "common/Types.hpp"
#include "common/Vector.hpp"
#include <vector>
#include <array>
#include <string>

namespace hvdc {
namespace fvm {

enum class FieldLocation : UInt8 {
    CellCenter = 0,
    FaceCenter = 1,
    Node = 2
};

template <typename T>
class Field {
public:
    Field() = default;
    Field(Index size, FieldLocation loc = FieldLocation::CellCenter,
          const std::string& name = "field");
    
    Field(Index size, T value, FieldLocation loc = FieldLocation::CellCenter,
          const std::string& name = "field");
    
    const std::string& name() const { return name_; }
    void set_name(const std::string& n) { name_ = n; }
    
    FieldLocation location() const { return location_; }
    
    Index size() const { return static_cast<Index>(data_.size()); }
    
    T& operator[](Index i) { return data_[i]; }
    const T& operator[](Index i) const { return data_[i]; }
    
    T& at(Index i) { return data_.at(i); }
    const T& at(Index i) const { return data_.at(i); }
    
    const std::vector<T>& raw() const { return data_; }
    std::vector<T>& raw() { return data_; }
    
    void resize(Index n) { data_.resize(n); }
    void resize(Index n, T val) { data_.resize(n, val); }
    
    void fill(T value) { std::fill(data_.begin(), data_.end(), value); }
    
    void assign(Index n, T value) { data_.assign(n, value); }
    
    void swap(Field& other) noexcept {
        std::swap(data_, other.data_);
        std::swap(name_, other.name_);
        std::swap(location_, other.location_);
    }
    
    typename std::vector<T>::iterator begin() { return data_.begin(); }
    typename std::vector<T>::iterator end() { return data_.end(); }
    typename std::vector<T>::const_iterator begin() const { return data_.begin(); }
    typename std::vector<T>::const_iterator end() const { return data_.end(); }

protected:
    std::vector<T> data_;
    std::string name_;
    FieldLocation location_ = FieldLocation::CellCenter;
};

using ScalarField = Field<Real>;
using VectorField = Field<Vec3>;
using TensorField = Field<Mat3x3>;

class VelocityField : public VectorField {
public:
    VelocityField() : VectorField() {}
    VelocityField(Index n_cells, const std::string& name = "U")
        : VectorField(n_cells, Vec3{0,0,0}, FieldLocation::CellCenter, name) {}
    
    ScalarField component(Index i) const;
    void set_component(Index i, const ScalarField& field);
    
    Real speed(Index cell_id) const;
    Real speed_squared(Index cell_id) const;
    
    ScalarField speed_field() const;
    ScalarField kinetic_energy() const;
    
    Vec3 face_interpolated(Index face_id, Real weight,
                            const Vec3& owner_val, const Vec3& neighbor_val) const;
};

class PressureField : public ScalarField {
public:
    PressureField() : ScalarField() {}
    PressureField(Index n_cells, Real p0 = 0.0, const std::string& name = "p")
        : ScalarField(n_cells, p0, FieldLocation::CellCenter, name) {}
    
    Real reference_pressure() const { return p_ref_; }
    void set_reference_pressure(Real p) { p_ref_ = p; }
    
    Real gauge_pressure(Index cell_id) const { return data_[cell_id] - p_ref_; }
    
    Vec3 pressure_gradient(Index cell_id, const class Grid& grid) const;

private:
    Real p_ref_ = 0.0;
};

class TurbulenceField {
public:
    TurbulenceField() = default;
    TurbulenceField(Index n_cells);
    
    ScalarField k;
    ScalarField epsilon;
    ScalarField omega;
    ScalarField nu_t;
    ScalarField mut;
    
    Real nut(Index cell_id) const { return nu_t[cell_id]; }
    Real Reynolds_stress_xx(Index cell_id, Real Ux_grad) const;
    
    void compute_nut_from_k_epsilon(Real C_mu = 0.09);
    void compute_nut_from_k_omega();
};

struct FluidState {
    VelocityField U;
    PressureField p;
    TurbulenceField turb;
    ScalarField T;
    ScalarField rho;
    ScalarField mu;
    ScalarField nu;
    
    Real rho_ref = 1.225;
    Real mu_ref = 1.789e-5;
    Real T_ref = 288.15;
    Real p_ref = 101325.0;
    
    Index size() const { return U.size(); }
    
    void resize(Index n);
    void initialize_defaults(Index n_cells);
    void update_thermophysical_properties();
};

} // namespace fvm
} // namespace hvdc

#endif // HVDC_FVM_FIELD_HPP
