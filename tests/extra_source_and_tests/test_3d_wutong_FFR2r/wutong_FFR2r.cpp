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

std::string full_path_to_stl_file_wall = "./input/wall.stl";
std::string full_path_to_stl_file_fluid = "./input/blood.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real scale = 0.001;
Real resolution_ref = 0.000125; // buffer size 1.25mm
Real fluid_radius_inlet = 0.0019384;
Real fluid_radius_outlet0 = 0.0009953;
Real fluid_radius_outlet1 = 0.0007774;
Real fluid_radius_outlet2 = 0.0008444;

//---------------------------------------------------------------------------------------------------------------
Vecd buffer_halfsize_inlet = Vecd(resolution_ref * 2.5, 1.5 * fluid_radius_inlet, 1.5 * fluid_radius_inlet);
Vecd buffer_translation_inlet = Vecd(-0.0133054 + resolution_ref * 2.5, 0.1909940, 1.5850189);
Vecd normal_vector_inlet = Vecd(0.3373, 0.9396, 0.0580);
Vecd vector_inlet = Vecd(0, -0.0616, 0.9981);
//---------------------------------------------------------------------------------------------------------------
Vecd buffer_halfsize_outlet0 = Vecd(resolution_ref * 2.5, 2 * fluid_radius_outlet0, 2 * fluid_radius_outlet0);
Vecd buffer_translation_outlet0 = Vecd(-0.0202521 + resolution_ref * 2.5, 0.1604702, 1.5547743);
Vecd normal_vector_outlet0 = Vecd(0.2513, 0.8882, -0.3846);
Vecd vector_outlet0 = Vecd(0, 0.3974, 0.9177);
//---------------------------------------------------------------------------------------------------------------
Vecd buffer_halfsize_outlet1 = Vecd(resolution_ref * 2.5, 2 * fluid_radius_outlet1, 2 * fluid_radius_outlet1);
Vecd buffer_translation_outlet1 = Vecd(-0.0486019 + resolution_ref * 2.5, 0.1902887, 1.5416359);
Vecd normal_vector_outlet1 = Vecd(0.9284, -0.3714, 0.0125);
Vecd vector_outlet1 = Vecd(0, -0.0336, -0.9994);
//---------------------------------------------------------------------------------------------------------------
Vecd buffer_halfsize_outlet2 = Vecd(resolution_ref * 2.5, 2 * fluid_radius_outlet2, 2 * fluid_radius_outlet2);
Vecd buffer_translation_outlet2 = Vecd(-0.0160551 + resolution_ref * 2.5, 0.2352654, 1.55536043);
Vecd normal_vector_outlet2 = Vecd(0.8494, -0.2059, 0.4859);
Vecd vector_outlet2 = Vecd(0, -0.9207, -0.3902);
//---------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------
Real full_length = 0.81;
Vecd translation(0, 0, 0);

//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real Outlet_pressure = 0;
Real rho0_f = 1060.0;
Real mu_f = 3.6e-3;
Real U_f = 1;
Real c_f = 10.0 * U_f;

//----------------------------------------------------------------------
// Geometry parameters for shell.
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(Vec3d(-0.05, 0.15, 1.53), Vec3d(0.02, 0.24, 1.6));

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

struct InflowPressure
{
    template <class BoundaryConditionType>
    InflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        return p_;
    }
};
struct OutflowPressure0
{
    template <class BoundaryConditionType>
    OutflowPressure0(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        /*constant pressure*/
        Real pressure = Outlet_pressure;
        return pressure;
    }
};

struct OutflowPressure1
{
    template <class BoundaryConditionType>
    OutflowPressure1(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        /*constant pressure*/
        Real pressure = Outlet_pressure;
        return pressure;
    }
};
struct OutflowPressure2
{
    template <class BoundaryConditionType>
    OutflowPressure2(BoundaryConditionType &boundary_condition) {}

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

        target_velocity[0] = 0.457;
        target_velocity[1] = 0.0;
        target_velocity[2] = 0.0;
        return target_velocity;
    }
};

int main(int ac, char *av[])
{

    //----------------------------------------------------------------------
    //  Build up -- a SPHSystem --
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(false);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("InsertedBody"));
    // wall_boundary.defineAdaptationRatios(1.15, 2.0);
    wall_boundary.defineBodyLevelSetShape()->correctLevelSetSign();
    wall_boundary.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.getName())
        : wall_boundary.generateParticles<BaseParticles, Lattice>();

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
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                wall_boundary_write_inserted_body_to_vtp.writeToFile(ite_p);
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
    ContactRelation water_block_contact(water_block, {&wall_boundary});
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Algorithms of solid dynamics.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
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

    BodyAlignedBoxByCell disposer_outlet0(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(-normal_vector_outlet0)), -vector_outlet0),
                                                                                             Vecd(buffer_translation_outlet0)),
                                                                                   buffer_halfsize_outlet0));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_outlet0(disposer_outlet0, xAxis);

    BodyAlignedBoxByCell disposer_outlet1(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(-normal_vector_outlet1)), -vector_outlet1),
                                                                                             Vecd(buffer_translation_outlet1)),
                                                                                   buffer_halfsize_outlet1));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_outlet1(disposer_outlet1, xAxis);

    BodyAlignedBoxByCell disposer_outlet2(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(-normal_vector_outlet2)), -vector_outlet2),
                                                                                             Vecd(buffer_translation_outlet2)),
                                                                                   buffer_halfsize_outlet2));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_outlet2(disposer_outlet2, xAxis);

    
    //--------------------------------------------------------------------------------------------------------------------

    BodyAlignedBoxByCell inflow_emitter(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_inlet)), vector_inlet),
                                                                                           Vecd(buffer_translation_inlet)),
                                                                                 buffer_halfsize_inlet));
    fluid_dynamics::NonPrescribedPressureBidirectionalBuffer inflow_injection(inflow_emitter, in_outlet_particle_buffer, xAxis);

    BodyAlignedBoxByCell outflow_emitter0(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_outlet0)), vector_outlet0),
                                                                                             Vecd(buffer_translation_outlet0)),
                                                                                   buffer_halfsize_outlet0));
    fluid_dynamics::BidirectionalBuffer<OutflowPressure0> outflow_injection0(outflow_emitter0, in_outlet_particle_buffer, xAxis);

    BodyAlignedBoxByCell outflow_emitter1(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_outlet1)), vector_outlet1),
                                                                                             Vecd(buffer_translation_outlet1)),
                                                                                   buffer_halfsize_outlet1));
    fluid_dynamics::BidirectionalBuffer<OutflowPressure1> outflow_injection1(outflow_emitter1, in_outlet_particle_buffer, xAxis);

    BodyAlignedBoxByCell outflow_emitter2(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_outlet2)), vector_outlet2),
                                                                                             Vecd(buffer_translation_outlet2)),
                                                                                   buffer_halfsize_outlet2));
    fluid_dynamics::BidirectionalBuffer<OutflowPressure2> outflow_injection2(outflow_emitter2, in_outlet_particle_buffer, xAxis);

    

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_block_contact);
    SimpleDynamics<fluid_dynamics::PressureCondition<InflowPressure>> inflow_pressure_condition(inflow_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(inflow_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<OutflowPressure0>> outflow_pressure_condition0(outflow_emitter0);
    SimpleDynamics<fluid_dynamics::PressureCondition<OutflowPressure1>> outflow_pressure_condition1(outflow_emitter1);
    SimpleDynamics<fluid_dynamics::PressureCondition<OutflowPressure2>> outflow_pressure_condition2(outflow_emitter2);
   
    //----------------------------------------------------------------------
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<int>("Indicator");
    water_block.addBodyStateForRecording<Real>("Density");
    water_block.addBodyStateForRecording<int>("BufferParticleIndicator");
    BodyStatesRecordingToVtp body_states_recording(sph_system.real_bodies_);

    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    boundary_indicator.exec();
    inflow_injection.tag_buffer_particles.exec();
    outflow_injection0.tag_buffer_particles.exec();
    outflow_injection1.tag_buffer_particles.exec();
    outflow_injection2.tag_buffer_particles.exec();
 
    wall_boundary_normal_direction.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 200;
    Real end_time = 10.0;     /**< End time. */
    Real Output_Time = 0.005; /**< Time stamps for output of body states. */
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

            interval_computing_time_step += TickCount::now() - time_instance;

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                /** Fluid pressure relaxation */
                pressure_relaxation.exec(dt);
                kernel_summation.exec();
                inflow_pressure_condition.exec(dt);
                inflow_velocity_condition.exec();
                outflow_pressure_condition0.exec(dt);
                outflow_pressure_condition1.exec(dt);
                outflow_pressure_condition2.exec(dt);
           
                /** Fluid density relaxation */
                density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << " dt = " << dt << "\n";
            }
            number_of_iterations++;
            time_instance = TickCount::now();

            inflow_injection.injection.exec();
            outflow_injection0.injection.exec();
            outflow_injection1.injection.exec();
            outflow_injection2.injection.exec();
          

            disposer_outflow_inlet.exec();
            disposer_outflow_outlet0.exec();
            disposer_outflow_outlet1.exec();
            disposer_outflow_outlet2.exec();
          

            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_complex.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
            boundary_indicator.exec();

            inflow_injection.tag_buffer_particles.exec();
            outflow_injection0.tag_buffer_particles.exec();
            outflow_injection1.tag_buffer_particles.exec();
            outflow_injection2.tag_buffer_particles.exec();
           
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
    return 0;
}
