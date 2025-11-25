/**
 * @file 	mixed_poiseuille_flow.cpp
 * @brief 	2D mixed poiseuille flow example
 * @details This is the one of the basic test cases for mixed pressure/velocity in-/outlet boundary conditions.
 * @author 	Shuoguo Zhang and Xiangyu Hu
 */
//调成脉动速度，出口压力都为零，入口速度给了，流体密度粘性换掉，模量的公式也换掉（查文章）入口压力时间脉动（已经完成），参考脉动泊肃叶流。
#include "sphinxsys.h"
#include <gtest/gtest.h>
#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"



using namespace SPH;
//----------------------------------------------------------------------
//	Domain bounds of the system.比原来的在xz上放大了一倍
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const Real scale = 0.001;
const Real diameter = 6.35 * scale;
const Real fluid_radius = 0.5 * diameter;
const Real full_length = 10 * fluid_radius;

    //----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real Inlet_pressure = 0.2;
Real Outlet_pressure = 0.0;
const Real rho0_f = 1050.0; /**< Reference density of fluid. */
const Real mu_f = 3.6e-3;   /**< Viscosity. */
const Real Re = 100;
const Real U_f = Re * mu_f / rho0_f / diameter;
const Real U_max = 2.0 * U_f;  // parabolic inflow, Thus U_max = 2*U_f
const Real c_f = 10.0 * U_max; /**< Reference sound speed. */

Real rho0_s = 1000.0; /**< Reference density.*/
Real poisson = 0.4; /**< Poisson ratio.*/
Real Ae = 1.4e3;
; /**< Normalized Youngs Modulus. */
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
const int number_of_particles = 10;
const Real resolution_ref = diameter / number_of_particles;
const Real resolution_shell = 0.5 * resolution_ref;
const Real shell_thickness = 0.5 * resolution_shell;
//----------------------------------------------------------------------
//	Geometry parameters for shell.
//----------------------------------------------------------------------

const Real inflow_length = resolution_ref * 10.0; // Inflow region
const Real wall_thickness = resolution_ref * 4.0;
const int SimTK_resolution = 20;
const Vec3d translation_fluid(0., full_length * 0.5, 0.);
const BoundingBox system_domain_bounds(Vec3d(-diameter, 0, -diameter) - Vec3d(shell_thickness, wall_thickness,
                                                                              shell_thickness),
                                       Vec3d(diameter, full_length, diameter) +
                                           Vec3d(shell_thickness, wall_thickness,
                                                 shell_thickness));
//定义buffer几何参数
Vecd buffer_halfsize_inlet = Vecd(fluid_radius, resolution_ref * 2, fluid_radius);
Vecd buffer_translation_inlet = Vecd(0.0, resolution_ref * 2, 0.0);
Vecd normal_vector_inlet = Vecd(0, 1, 0);
Vecd vector_inlet = Vecd(0,0, 1);
// 
Vecd buffer_halfsize_outlet = Vecd(fluid_radius, resolution_ref * 2, fluid_radius);
Vecd buffer_translation_outlet = Vecd(0, full_length - buffer_halfsize_outlet[1], 0);
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
        Vec3d translation_fluid(0., full_length * 0.5, 0.);
        add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0., 1., 0.), fluid_radius,
                                       full_length * 0.5, SimTK_resolution,
                                       translation_fluid);
    }
};
// 定义壳体
class ShellBoundary;
template <>
class ParticleGenerator<ShellBoundary> : public ParticleGenerator<Surface>
{
    Real resolution_shell_;
    Real wall_thickness_;
    Real shell_thickness_;

  public:
    explicit ParticleGenerator(SPHBody &sph_body, Real resolution_shell, Real wall_thickness, Real shell_thickness)
        : ParticleGenerator<Surface>(sph_body),
          resolution_shell_(resolution_shell),
          wall_thickness_(wall_thickness), shell_thickness_(shell_thickness){};
    void initializeGeometricVariables() override
    {
        const Real radius_mid_surface = fluid_radius + resolution_shell_ * 0.5;
        const auto particle_number_mid_surface =
            int(2.0 * radius_mid_surface * Pi / resolution_shell_);
        const auto particle_number_height =
            int((full_length + 2.0 * wall_thickness_) / resolution_shell_);
        for (int i = 0; i < particle_number_mid_surface; i++)
        {
            for (int j = 0; j < particle_number_height; j++)
            {
                Real theta = (i + 0.5) * 2 * Pi / (Real)particle_number_mid_surface;
                Real x = radius_mid_surface * cos(theta);
                Real y = -wall_thickness_ + (full_length + 2 * wall_thickness_) * j / (Real)particle_number_height + 0.5 * resolution_shell_;
                Real z = radius_mid_surface * sin(theta);
                initializePositionAndVolumetricMeasure(Vec3d(x, y, z),
                                                       resolution_shell_ * resolution_shell_);
                Vec3d n_0 = Vec3d(x / radius_mid_surface, 0.0, z / radius_mid_surface);
                initializeSurfaceProperties(n_0, shell_thickness_);
            }
        }
    }
};


////定义出入口条件
// struct LeftInflowPressure // 脉动压力入口
//{
//     template <class BoundaryConditionType>
//     LeftInflowPressure(BoundaryConditionType &boundary_condition) {}
//
//     Real operator()(Real &p_)
//     {
//         /*pulsatile pressure*/
//         Real pressure = Inlet_pressure * cos(3 * Pi * GlobalStaticVariables::physical_time_);
//         /*constant pressure*/
//         //        Real pressure = Inlet_pressure;
//         return pressure;
//     }
// };
//
// struct RightInflowPressure//恒定压力出口
//{
//     template <class BoundaryConditionType>
//     RightInflowPressure(BoundaryConditionType &boundary_condition) {}
//
//     Real operator()(Real &p_)
//     {
//         /*constant pressure*/
//         Real pressure = Outlet_pressure;
//         return pressure;
//     }
// };
// } // namespace SPH
struct InletflowPressure
{
    template <class BoundaryConditionType>
    InletflowPressure(BoundaryConditionType &boundary_condition) {}

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
        Real run_time = GlobalStaticVariables::physical_time_;
        Real pressure = run_time < 0.1 ? 0.5 * Outlet_pressure * (1.0 - cos(Pi * run_time / 0.1)) : Outlet_pressure;
        return pressure;
    }
};

struct InflowVelocity
{
    Real u_ref_, t_ref_, u_ave;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(0.1), t_ref_(0.1), u_ave(0.0) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = Vecd::Zero();
        Real run_time = GlobalStaticVariables::physical_time_;
        if (run_time < t_ref_)
        {
            u_ave = u_ref_;
        }
        else if (((run_time - t_ref_) - int((run_time - t_ref_) / 0.5) * 0.5) >= 0.0 && ((run_time - t_ref_) - int((run_time - t_ref_) / 0.5) * 0.5) <= 0.218)
        {
            u_ave = 0.5 * sin(4.0 * Pi * ((run_time - t_ref_) + 0.0160236));
        }
        else
        {
            u_ave = u_ref_;
        }

        target_velocity[0] = u_ave;
        target_velocity[1] = 0.0;
        target_velocity[2] = 0.0;

        return target_velocity;
    }
};









//----------------------------------------------------------------------
//  Define constrain shape自己写的，几何参数待定
//----------------------------------------------------------------------



class ConstrainBlockShape : public ComplexShape
{
  public:
    explicit ConstrainBlockShape(const std::string &shape_name) : ComplexShape(shape_name)
    {

        add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0., 1., 0.), fluid_radius * 2,
                                       full_length * 0.2, SimTK_resolution,
                                       buffer_translation_inlet);

        add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0., 1., 0.), fluid_radius * 2,
                                       full_length * 0.2, SimTK_resolution,
                                       buffer_translation_outlet);
    }
};
 int main(int ac, char *av[])
 {



    //----------------------------------------------------------------------
    //  Build up -- a SPHSystem --
    //----------------------------------------------------------------------
    SPHSystem system(system_domain_bounds, resolution_ref);
    IOEnvironment io_environment(system);
    system.setGenerateRegressionData(false);
    system.setRunParticleRelaxation(false);
    system.setReloadParticles(true);
    system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);

    SolidBody shell_boundary(system, makeShared<DefaultShape>("Shell"));
    shell_boundary.defineAdaptation<SPH::SPHAdaptation>(1.15, resolution_ref / resolution_shell);
    /*shell_boundary.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);*/
    shell_boundary.defineMaterial<LinearElasticSolid>(1, 1e3, 0.45);
    shell_boundary.generateParticles<SurfaceParticles, ShellBoundary>(resolution_shell, wall_thickness, shell_thickness);

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    
    InnerRelation shell_boundary_inner(shell_boundary);
    ShellInnerRelationWithContactKernel wall_curvature_inner(shell_boundary, water_block);
    // shell normal should point from fluid to shell
    // normal corrector set to false if shell normal is already pointing from fluid to shell
    ContactRelationToShell water_block_contact(water_block, {&shell_boundary}, {true});
    ContactRelationFromShell shell_boundary_contact(shell_boundary, {&water_block}, {true}); // 加了一个来自壁面的作用
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Algorithms of solid dynamics.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> shell_boundary_normal_direction(shell_boundary);
    /** bidrectional buffer */
    /** delete outflow particles */
    BodyAlignedBoxByCell disposer_inlet(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(-normal_vector_inlet)),
                                                                                                      -vector_inlet),
                                                                                           Vecd(buffer_translation_inlet)),
                                                                                 buffer_halfsize_inlet));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_inlet(disposer_inlet, xAxis);
    BodyAlignedBoxByCell disposer_outlet(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_outlet)),
                                                                                                       vector_outlet),
                                                                                            Vecd(buffer_translation_outlet)),
                                                                                  buffer_halfsize_outlet));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_outlet(disposer_outlet, xAxis);
    /** surface particle identification */
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_indicator(water_block_inner, water_block_contact);
    BodyAlignedBoxByCell inflow_emitter(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_inlet)),
                                                                      vector_inlet),
                                                           Vec3d(buffer_translation_inlet)),
                                                 buffer_halfsize_inlet));
    fluid_dynamics::NonPrescribedPressureBidirectionalBuffer inflow_injection(inflow_emitter, in_outlet_particle_buffer, xAxis);

    BodyAlignedBoxByCell outflow_emitter(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(-normal_vector_outlet)),
                                                                      vector_outlet),
                                                           Vec3d(buffer_translation_outlet)),
                                                 buffer_halfsize_outlet));
    fluid_dynamics::BidirectionalBuffer<OutflowPressure> outflow_injection(outflow_emitter, in_outlet_particle_buffer, xAxis); // 这里写的不一样
  

    
    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_block_contact);
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);

    
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(inflow_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<InletflowPressure>> inflow_pressure_condition(inflow_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<OutflowPressure>> outflow_pressure_condition(outflow_emitter);

    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_block_contact);


    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> shell_boundary_corrected_configuration(shell_boundary_inner); // 不太清楚作用
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> shell_boundary_stress_relaxation_first_half(shell_boundary_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> shell_boundary_stress_relaxation_second_half(shell_boundary_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> shell_boundary_computing_time_step_size(shell_boundary);
    ConstrainBlockShape constrained_shape("a");
    BodyRegionByParticle constraint_in_outlet(shell_boundary, constrained_shape);
    SimpleDynamics<FixBodyPartConstraint> constraint_in_outlet_base(constraint_in_outlet);
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> wall_corrected_configuration(shell_boundary_inner);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_curvature(wall_curvature_inner);
    
    //----------------------------------------------------------------------
    //	Algorithms of fluid dynamics.
    //----------------------------------------------------------------------
   


    



    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(shell_boundary);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> wall_boundary_update_normal(shell_boundary);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(shell_boundary_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(shell_boundary_contact);
  
    water_block.addBodyStateForRecording<int>("Indicator");
    water_block.addBodyStateForRecording<Real>("Density");
    water_block.addBodyStateForRecording<int>("BufferParticleIndicator");
    shell_boundary.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
    shell_boundary.addBodyStateForRecording<Real>("Average2ndPrincipleCurvature");

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------

        BodyStatesRecordingToVtp body_states_recording(system.real_bodies_);

    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    boundary_indicator.exec();
    inflow_injection.tag_buffer_particles.exec();
    outflow_injection.tag_buffer_particles.exec();

    shell_boundary_normal_direction.exec();
    wall_corrected_configuration.exec();
    shell_curvature.exec();
    water_block_complex.updateConfiguration();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 2.0;                 /**< End time. */
    Real Output_Time = end_time / 100.0; /**< Time stamps for output of body states. */
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
            boundary_indicator.exec();
            update_fluid_density.exec();
            viscous_acceleration.exec();
            transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;
            /** FSI for viscous force. */
            viscous_force_from_fluid.exec();
            /** Update normal direction on elastic body.*/
            wall_boundary_update_normal.exec();

            size_t inner_ite_dt = 0;
            size_t inner_ite_dt_s = 0;

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                /** Fluid pressure relaxation */
                pressure_relaxation.exec(dt);
                kernel_summation.exec();

                /*inflow_velocity_condition.exec();*/
                /** FSI for pressure force. */
                pressure_force_from_fluid.exec();
                /** Fluid density relaxation */
                inflow_pressure_condition.exec(dt);
                inflow_velocity_condition.exec();
                outflow_pressure_condition.exec(dt);
             
                density_relaxation.exec(dt);

                /** Solid dynamics. */
                inner_ite_dt_s = 0;
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = SMIN(shell_boundary_computing_time_step_size.exec(), dt - dt_s_sum);
                    shell_boundary_stress_relaxation_first_half.exec(dt_s);
                    constraint_in_outlet_base.exec();
                    shell_boundary_stress_relaxation_second_half.exec(dt_s);
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
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "	dt / dt_s = " << inner_ite_dt_s << "\n";
            }
            number_of_iterations++;

            inflow_injection.injection.exec();
            outflow_injection.injection.exec();
            disposer_outflow_inlet.exec();
            disposer_outflow_outlet.exec();

            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_contact.updateConfiguration();
            water_block_complex.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
            boundary_indicator.exec();
            inflow_injection.tag_buffer_particles.exec();
            outflow_injection.tag_buffer_particles.exec();

            shell_boundary.updateCellLinkedList();
            shell_boundary_contact.updateConfiguration();
            boundary_indicator.exec();
        }
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();

        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9)
              << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9)
              << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9)
              << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    return 0;
 }
