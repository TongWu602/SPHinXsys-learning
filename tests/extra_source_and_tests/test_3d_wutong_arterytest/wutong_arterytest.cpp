/**
 * @file 	mixed_poiseuille_flow.cpp
 * @brief 	2D mixed poiseuille flow example
 * @details This is the one of the basic test cases for mixed pressure/velocity in-/outlet boundary conditions.
 * @author 	Shuoguo Zhang and Xiangyu Hu
 */

#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"

using namespace SPH;

std::string full_path_to_stl_file_wall = "./input/wall1long.stl";
std::string full_path_to_stl_file_fluid = "./input/blood1long.stl";
std::string full_path_to_stl_file_constrain = "./input/constrain1long.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real scale = 0.0001;
Real fluid_radius_inlet = 0.0125;
Real fluid_radius_outlet = 0.0155;
Real full_length = 0.81;
Vecd translation(0, 0, 0);
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real Outlet_pressure = 0.0;
Real rho0_f = 1100.0;
Real mu_f = 3.6e-3;
Real U_f = 0.4;
Real c_f = 10.0 * U_f;

Real rho0_s = 1000.0;
Real poisson = 0.45;
Real Ae = 1.4e3;
Real Youngs_modulus = 4.8e6;

Real resolution_ref = 0.0005;
//----------------------------------------------------------------------
//	Geometry parameters for shell.
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(Vec3d(-0.025, 0.0325, -0.025), Vec3d(0.025, 0.0905, 0.025));
Vecd buffer_halfsize_inlet = Vecd(resolution_ref * 2.5, 1.5 * fluid_radius_inlet, 1.5 * fluid_radius_inlet);
Vecd buffer_translation_inlet = Vecd(0.0, 0.0325 + resolution_ref * 2.5, 0.0);
Vecd normal_vector_inlet = Vecd(0, 1, 0);
Vecd vector_inlet = Vecd(0, 0, 1);

Vecd buffer_halfsize_outlet = Vecd(resolution_ref * 2.5, 1.5 * fluid_radius_outlet, 1.5 * fluid_radius_outlet);
Vecd buffer_translation_outlet = Vecd(0, 0.0905 - resolution_ref * 2.5, 0);
Vecd normal_vector_outlet = Vecd(0, 1, 0);
Vecd vector_outlet = Vecd(0, 0, 1);

//----------------------------------------------------------------------
//  Define water shape
//----------------------------------------------------------------------
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_stl_file_fluid, translation, scale);
    }
};

class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_stl_file_wall, translation, scale);
    }
};

class ConstrainBlock : public ComplexShape
{
  public:
    explicit ConstrainBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_stl_file_constrain, translation, scale);
    }
};

struct InflowPressure
{
    template <class BoundaryConditionType>
    InflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        return p_;
    }
};
struct OutflowPressure
{
    template <class BoundaryConditionType>
    OutflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        /*constant pressure*/
        Real pressure = Outlet_pressure;
        return pressure;
    }
};

struct InflowVelocity
{
    Real u_ref_, t_ref_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(0.151), t_ref_(0.5) {}

    Real period = 1.0;

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = Vecd::Zero();
        Real run_time = GlobalStaticVariables::physical_time_;

        Real time_in_period = fmod(run_time, period);
        Real radius_squared = position[1] * position[1] + position[2] * position[2];
        /*  if (radius_squared < 0.0125 * 0.0125)
          {*/
        target_velocity[0] = ((10815643.4252119 * pow(time_in_period, 15) +
                               -81009632.4398252 * pow(time_in_period, 14) +
                               268856614.120370 * pow(time_in_period, 13) +
                               -519061147.630914 * pow(time_in_period, 12) +
                               642001648.473172 * pow(time_in_period, 11) +
                               -525607410.449130 * pow(time_in_period, 10) +
                               280836392.561131 * pow(time_in_period, 9) +
                               -88978528.7318248 * pow(time_in_period, 8) +
                               9648311.68913681 * pow(time_in_period, 7) +
                               4106564.92095778 * pow(time_in_period, 6) +
                               -1911070.76313921 * pow(time_in_period, 5) +
                               325263.426444146 * pow(time_in_period, 4) +
                               -23273.4128507954 * pow(time_in_period, 3) +
                               618.643903074195 * pow(time_in_period, 2) +
                               6.16233963952450 * time_in_period +
                               0.0132187683796114) *
                              0.0339) *
                             2 * (1.0 - radius_squared / (0.0125 * 0.0125));
        /*}
        else
        {


                target_velocity[0] = 0.0;


        }*/

        target_velocity[1] = 0.0;
        target_velocity[2] = 0.0;
        return target_velocity;
    }
};
class FluidObserver;
template <>
class ParticleGenerator<FluidObserver> : public ParticleGenerator<Observer>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body) : ParticleGenerator<Observer>(sph_body)
    {
        /** A line of measuring points at the entrance of the channel. */
        size_t number_observation_points = 51;
        Real range_of_measure = 0.31;
        Real start_of_measure = -0.155;
        /** the measuring locations */
        for (size_t i = 0; i < number_observation_points; ++i)
        {
            Vecd point_coordinate(0, 0.071, range_of_measure * (Real)i / (Real)(number_observation_points - 1) + start_of_measure);
            positions_.push_back(point_coordinate);
        }
    }
};
class FluidObserver2;
template <>
class ParticleGenerator<FluidObserver2> : public ParticleGenerator<Observer>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body) : ParticleGenerator<Observer>(sph_body)
    {
        /** A line of measuring points at the entrance of the channel. */
        size_t number_observation_points = 51;
        Real range_of_measure = 0.31;
        Real start_of_measure = -0.155;
        /** the measuring locations */
        for (size_t i = 0; i < number_observation_points; ++i)
        {
            Vecd point_coordinate(0, 0.079, range_of_measure * (Real)i / (Real)(number_observation_points - 1) + start_of_measure);
            positions_.push_back(point_coordinate);
        }
    }
};
class FluidObserver3;
template <>
class ParticleGenerator<FluidObserver3> : public ParticleGenerator<Observer>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body) : ParticleGenerator<Observer>(sph_body)
    {
        /** A line of measuring points at the entrance of the channel. */
        size_t number_observation_points = 11;
        Real range_of_measure = 0.31;
        Real start_of_measure = -0.155;
        /** the measuring locations */
        for (size_t i = 0; i < number_observation_points; ++i)
        {
            Vecd point_coordinate(0, 0.085, range_of_measure * (Real)i / (Real)(number_observation_points - 1) + start_of_measure);
            positions_.push_back(point_coordinate);
        }
    }
};
int main(int ac, char *av[])
{

    //----------------------------------------------------------------------
    //  Build up -- a SPHSystem --
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(true);
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("InsertedBody"));
    wall_boundary.defineAdaptationRatios(1.15, 2.0);
    wall_boundary.defineBodyLevelSetShape()->correctLevelSetSign();
    wall_boundary.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.getName())
        : wall_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<BaseParticles, FluidObserver>();

    ObserverBody fluid_observer2(sph_system, "FluidObserver2");
    fluid_observer2.generateParticles<BaseParticles, FluidObserver2>();

    ObserverBody fluid_observer3(sph_system, "FluidObserver3");
    fluid_observer3.generateParticles<BaseParticles, FluidObserver3>();

    if (sph_system.RunParticleRelaxation())
    {
        /** body topology only for particle relaxation */
        InnerRelation wall_boundary_inner(wall_boundary);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> wall_boundary_random_inserted_body_particles(wall_boundary);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp wall_boundary_write_inserted_body_to_vtp({&wall_boundary});
        /** Write the particle reload files. */
        ReloadParticleIO wall_boundary_write_particle_reload_files({&wall_boundary});
        /** A  Physics relaxation step. */
        RelaxationStepInner wall_boundary_relaxation_step_inner(wall_boundary_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        wall_boundary_random_inserted_body_particles.exec(0.25);
        wall_boundary_relaxation_step_inner.SurfaceBounding().exec();
        wall_boundary_write_inserted_body_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            
            wall_boundary_relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 1 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";

                if (ite_p % 200 == 0)
                {

                  
                    wall_boundary_write_inserted_body_to_vtp.writeToFile(ite_p);
                }
            }
        }

        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        wall_boundary_write_particle_reload_files.writeToFile(0);
        return 0;
    }

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation wall_boundary_inner(wall_boundary);
    ContactRelation water_block_contact(water_block, {&wall_boundary});
    ContactRelation wall_boundary_contact(wall_boundary, {&water_block});
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    ContactRelation fluid_observer_contact2(fluid_observer2, {&water_block});
    ContactRelation fluid_observer_contact3(fluid_observer3, {&water_block});
    //----------------------------------------------------------------------
    //	Algorithms of solid dynamics.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> wall_boundary_corrected_configuration(wall_boundary_inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> wall_boundary_stress_relaxation_first_half(wall_boundary_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> wall_boundary_stress_relaxation_second_half(wall_boundary_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> wall_boundary_computing_time_step_size(wall_boundary);
    BodyRegionByParticle constraint_wall(wall_boundary, makeShared<ConstrainBlock>("Holder"));
    SimpleDynamics<FixBodyPartConstraint> constraint_wall_base(constraint_wall);
    //----------------------------------------------------------------------
    //	Algorithms of fluid dynamics.
    //----------------------------------------------------------------------
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_indicator(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_block_contact);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);

    /** delete outflow particles */
    BodyAlignedBoxByCell disposer_inlet(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(-normal_vector_inlet)), -vector_inlet),
                                                                                           Vecd(buffer_translation_inlet)),
                                                                                 buffer_halfsize_inlet));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_inlet(disposer_inlet, xAxis);
    BodyAlignedBoxByCell disposer_outlet(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_outlet)), vector_outlet),
                                                                                            Vecd(buffer_translation_outlet)),
                                                                                  buffer_halfsize_outlet));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_outlet(disposer_outlet, xAxis);

    BodyAlignedBoxByCell inflow_emitter(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_inlet)), vector_inlet),
                                                                                           Vecd(buffer_translation_inlet)),
                                                                                 buffer_halfsize_inlet));
    fluid_dynamics::NonPrescribedPressureBidirectionalBuffer inflow_injection(inflow_emitter, in_outlet_particle_buffer, xAxis);
    BodyAlignedBoxByCell outflow_emitter(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(-normal_vector_outlet)), -vector_outlet),
                                                                                            Vecd(buffer_translation_outlet)),
                                                                                  buffer_halfsize_outlet));
    fluid_dynamics::BidirectionalBuffer<OutflowPressure> outflow_injection(outflow_emitter, in_outlet_particle_buffer, xAxis);

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_block_contact);
    SimpleDynamics<fluid_dynamics::PressureCondition<InflowPressure>> inflow_pressure_condition(inflow_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<OutflowPressure>> outflow_pressure_condition(outflow_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(inflow_emitter);

    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(wall_boundary);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> wall_boundary_update_normal(wall_boundary);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(wall_boundary_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(wall_boundary_contact);

    //----------------------------------------------------------------------
    ObservedQuantityRecording<Vecd> write_fluid_velocity("Velocity", fluid_observer_contact);
    ObservedQuantityRecording<Vecd> write_fluid_velocity2("Velocity", fluid_observer_contact2);
    ObservedQuantityRecording<Vecd> write_fluid_velocity3("Velocity", fluid_observer_contact3);

    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<int>("Indicator");
    water_block.addBodyStateForRecording<Real>("Density");
    water_block.addBodyStateForRecording<int>("BufferParticleIndicator");
    wall_boundary.addDerivedBodyStateForRecording<VonMisesStress>();
    BodyStatesRecordingToVtp body_states_recording(sph_system.real_bodies_);

    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    boundary_indicator.exec();
    inflow_injection.tag_buffer_particles.exec();
    outflow_injection.tag_buffer_particles.exec();
    wall_boundary_normal_direction.exec();
    wall_boundary_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 10.0;    /**< End time. */
    Real Output_Time = 0.01; /**< Time stamps for output of body states. */
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    body_states_recording.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            time_instance = TickCount::now();
            Real Dt = get_fluid_advection_time_step_size.exec();

            update_fluid_density.exec();
            viscous_acceleration.exec();
            transport_velocity_correction.exec();

            viscous_force_from_fluid.exec();
            wall_boundary_update_normal.exec();
            interval_computing_time_step += TickCount::now() - time_instance;
            size_t inner_ite_dt = 0;
            size_t inner_ite_dt_s = 0;

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                /** Fluid pressure relaxation */
                pressure_relaxation.exec(dt);
                kernel_summation.exec();
                inflow_pressure_condition.exec(dt);
                inflow_velocity_condition.exec();
                outflow_pressure_condition.exec(dt);
                pressure_force_from_fluid.exec();
                /** Fluid density relaxation */
                density_relaxation.exec(dt);

                /** Solid dynamics. */
                inner_ite_dt_s = 0;
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = SMIN(wall_boundary_computing_time_step_size.exec(), dt - dt_s_sum);
                    wall_boundary_stress_relaxation_first_half.exec(dt_s);
                    constraint_wall_base.exec();
                    wall_boundary_stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                    inner_ite_dt_s++;
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                inner_ite_dt++;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << " dt = " << dt << "	dt / dt_s = " << inner_ite_dt_s << "\n";
            }
            number_of_iterations++;
            time_instance = TickCount::now();
            inflow_injection.injection.exec();
            outflow_injection.injection.exec();
            disposer_outflow_inlet.exec();
            disposer_outflow_outlet.exec();
            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_complex.updateConfiguration();
            wall_boundary.updateCellLinkedList();
            wall_boundary_contact.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
            boundary_indicator.exec();
            inflow_injection.tag_buffer_particles.exec();
            outflow_injection.tag_buffer_particles.exec();
        }
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        write_fluid_velocity.writeToFile(number_of_iterations);
        write_fluid_velocity2.writeToFile(number_of_iterations);
        write_fluid_velocity3.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    return 0;
}
