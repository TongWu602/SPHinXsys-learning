#include "surface_indication.hpp"

namespace SPH
{
//=================================================================================================//
FreeSurfaceIndication<Inner>::
    FreeSurfaceIndication(BaseInnerRelation &inner_relation)
    : FreeSurfaceIndication<GeneralDataDelegateInner>(inner_relation),
      smoothing_length_(inner_relation.getSPHBody().sph_adaptation_->ReferenceSmoothingLength()) {}
//=================================================================================================//
void FreeSurfaceIndication<Inner>::interaction(size_t index_i, Real dt)
{
    Real pos_div = 0.0;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        pos_div -= inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.r_ij_[n];
    }
    pos_div_[index_i] = pos_div;
}
//=================================================================================================//
void FreeSurfaceIndication<Inner>::update(size_t index_i, Real dt)
{
    indicator_[index_i] = 1;
    if (pos_div_[index_i] > threshold_by_dimensions_ && !isVeryNearFreeSurface(index_i))
        indicator_[index_i] = 0;
}
//=================================================================================================//
bool FreeSurfaceIndication<Inner>::isVeryNearFreeSurface(size_t index_i)
{
    bool is_near_surface = false;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        /** Two layer particles.*/
        if (pos_div_[inner_neighborhood.j_[n]] < threshold_by_dimensions_ &&
            inner_neighborhood.r_ij_[n] < smoothing_length_)
        {
            is_near_surface = true;
            break;
        }
    }
    return is_near_surface;
}
//=================================================================================================//
void FreeSurfaceIndication<Contact>::interaction(size_t index_i, Real dt)
{
    Real pos_div = 0.0;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            pos_div -= contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.r_ij_[n];
        }
    }
    pos_div_[index_i] += pos_div;
}
//=================================================================================================//
FreeSurfaceIndication<SpatialTemporal<Inner>>::
    FreeSurfaceIndication(BaseInnerRelation &inner_relation)
    : FreeSurfaceIndication<Inner>(inner_relation)
{
    particles_->registerVariable(previous_surface_indicator_, "PreviousSurfaceIndicator", 1);
    particles_->registerSortableVariable<int>("PreviousSurfaceIndicator");
}
//=================================================================================================//
void FreeSurfaceIndication<SpatialTemporal<Inner>>::interaction(size_t index_i, Real dt)
{
    FreeSurfaceIndication<Inner>::interaction(index_i, dt);

    if (pos_div_[index_i] < threshold_by_dimensions_ &&
        previous_surface_indicator_[index_i] != 1 &&
        !isNearPreviousFreeSurface(index_i))
        pos_div_[index_i] = 2.0 * threshold_by_dimensions_;
}
//=================================================================================================//
bool FreeSurfaceIndication<SpatialTemporal<Inner>>::isNearPreviousFreeSurface(size_t index_i)
{
    bool is_near_surface = false;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        if (previous_surface_indicator_[inner_neighborhood.j_[n]] == 1)
        {
            is_near_surface = true;
            break;
        }
    }
    return is_near_surface;
}
//=================================================================================================//
void FreeSurfaceIndication<SpatialTemporal<Inner>>::update(size_t index_i, Real dt)
{
    FreeSurfaceIndication<Inner>::update(index_i, dt);

    previous_surface_indicator_[index_i] = indicator_[index_i];
}
//=================================================================================================//
} // namespace SPH
