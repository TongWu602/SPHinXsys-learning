#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"

using namespace SPH;
Real H = 0.02;
Real N = 180.0;
Real resolution_ref = H / 200.0;
Real BW = resolution_ref * 4.0;
BoundingBox system_domain_bounds(Vec2d(-0.25*H, -0.125*H-BW), Vec2d(5.5 * H, 2.0 * H + BW));

//----------------------------------------------------------------------

Real Outlet_pressure = 0.0;
Real rho0_f = 1001.0;
Real rho0_s = 890.0;
Real Re = 800.0;
Real mu_f = 0.0043;
Real U_f = 0.6;
Real c_f = 10.0 * U_f;

Real poisson = 0.49;
Real Youngs_modulus = 1.5e6;
Real round_circle_radius = 0.004;
Vec2d insert_circle_center = Vec2d(3.0 * H, H);
Vec2d round_circle_center = Vec2d(0.08366, 0.024);
Vec2d round_circle_center2 = Vec2d(0.084, 0.0244);

//----------------------------------------------------------------------
Vec2d bidirectional_buffer_halfsize = Vec2d(2.5 * resolution_ref, 0.5 * 1.125*H);
Vec2d left_bidirectional_translation = bidirectional_buffer_halfsize + Vec2d(-0.25 * H, -0.125 * H);
Vec2d right_bidirectional_translation = Vec2d(5.5 * H, H) - bidirectional_buffer_halfsize;
Vec2d normal = Vec2d(1.0, 0.0);

//----------------------------------------------------------------------
struct LeftInflowPressure
{

    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        return p_;
    }
};

struct RightInflowPressure
{
    template <class BoundaryConditionType>
    RightInflowPressure(BoundaryConditionType &boundary_condition) {}

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
        : u_ref_(0.151), t_ref_(0.5)  {}

    Real period = 1.0;

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = Vecd::Zero();
        Real run_time = GlobalStaticVariables::physical_time_;

        Real time_in_period = fmod(run_time, period);

        target_velocity[0] =
           ( ( - 135376.864230860 * pow(time_in_period, 15) +
            791724.952620352 * pow(time_in_period, 14) +
            -1834811.29360339 * pow(time_in_period, 13) +
            1819481.65158154 * pow(time_in_period, 12) +
            312478.227377350 * pow(time_in_period, 11) +
            -3047690.77358122 * pow(time_in_period, 10) +
            3974321.09208313 * pow(time_in_period, 9) +
            -2869541.97590965 * pow(time_in_period, 8) +
            1315082.10948023 * pow(time_in_period, 7) +
            -392673.761210967 * pow(time_in_period, 6) +
            75403.5571688622 * pow(time_in_period, 5) +
            -9019.40331836968 * pow(time_in_period, 4) +
            646.458486451316 * pow(time_in_period, 3) +
            -24.4981207622299 * pow(time_in_period, 2) +
            0.520910865265311 * time_in_period +
             0.00844833949696557+0.0611)*(1.33 * (1.0 - (position[1]) * (position[1]) / 0.01125 / 0.01125)));
 
        target_velocity[1] = 0.0;
        return target_velocity;
    }
};


std::vector<Vecd> Water_Block_Shape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-0.25*H,-0.125*H ));
    water_block_shape.push_back(Vecd(-0.25 * H, H));
    for (size_t i = 0; i < N - 9.0; ++i)
    {
        water_block_shape.push_back(Vecd(insert_circle_center[0] - H * cos(i * Pi / N),
                                         insert_circle_center[1] + H * sin(i * Pi / N)));
    }
    water_block_shape.push_back(Vecd(0.08366, H));
    water_block_shape.push_back(Vecd(5.5 * H, H));
    water_block_shape.push_back(Vecd(5.5 * H, -0.125 * H));
    water_block_shape.push_back(Vecd(-0.25 * H, -0.125 * H));

    return water_block_shape;
}

std::vector<Vecd> Wall_Block_Shape()
{
    std::vector<Vecd> wall_block_shape;
    wall_block_shape.push_back(Vecd(-0.25 * H, -0.125 * H - BW));
    wall_block_shape.push_back(Vecd(-0.25 * H, H + BW));
    wall_block_shape.push_back(Vecd(2 * H - BW, H + BW));
    for (size_t i = 1; i < N - 10.0; ++i)
    {
        wall_block_shape.push_back(Vecd(insert_circle_center[0] - (H + BW) * cos(i * Pi / N),
                                        insert_circle_center[1] + (H + BW) * sin(i * Pi / N)));
    }
    wall_block_shape.push_back(Vecd(0.08366 + BW, H + BW));
    wall_block_shape.push_back(Vecd(5.5 * H, H + BW));
    wall_block_shape.push_back(Vecd(5.5 * H, -0.125 * H - BW));
    wall_block_shape.push_back(Vecd(-0.25 * H, -0.125 * H - BW));

    return wall_block_shape;
}
std::vector<Vecd> Wall_Block_Shape1()
{
    std::vector<Vecd> wall_block_shape1;
    wall_block_shape1.push_back(Vecd(0.07975, 0.02695));
    wall_block_shape1.push_back(Vecd(0.087, 0.02695));
    wall_block_shape1.push_back(Vecd(0.087, 0.02035));
    wall_block_shape1.push_back(Vecd(0.084, 0.02035));
    wall_block_shape1.push_back(Vecd(0.07995, 0.02405));
    wall_block_shape1.push_back(Vecd(0.07975, 0.02505));
    wall_block_shape1.push_back(Vecd(0.07975, 0.02695));
    return wall_block_shape1;
}

std::vector<Vecd> Valve_Shape()
{
    std::vector<Vecd> valve_shape3;
    valve_shape3.push_back(Vecd(0.03872, 0.0215));
    valve_shape3.push_back(Vecd(0.056367, 0.0038426));
    valve_shape3.push_back(Vecd(0.0561574, 0.0036305));
    valve_shape3.push_back(Vecd(0.0383, 0.0215));
    valve_shape3.push_back(Vecd(0.03872, 0.0215));

    return valve_shape3;
}
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(Water_Block_Shape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(Valve_Shape(), ShapeBooleanOps::sub);
        multi_polygon_.addACircle(round_circle_center, round_circle_radius, 100, ShapeBooleanOps::sub);
    }
};
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(Wall_Block_Shape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(Water_Block_Shape(), ShapeBooleanOps::sub);
        multi_polygon_.addACircle(round_circle_center, round_circle_radius, 100, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(Valve_Shape(), ShapeBooleanOps::sub);
        multi_polygon_.addACircle(round_circle_center2, round_circle_radius, 100, ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(Wall_Block_Shape1(), ShapeBooleanOps::sub);
    }
};

std::vector<Vecd> Block1()
{
    std::vector<Vecd> block1;
    block1.push_back(Vecd(0.038, H+BW));
    block1.push_back(Vecd(0.038, 0.022));
    block1.push_back(Vecd(2*H-BW, 0.022));
    block1.push_back(Vecd(2 * H - BW, H + BW));
    block1.push_back(Vecd(0.038, H + BW));

    return block1;
}

std::vector<Vecd> Block2()
{
    std::vector<Vecd> block2;
    block2.push_back(Vecd(0.0395, 0.020));
    block2.push_back(Vecd(0.0395, 0.0205));
    block2.push_back(Vecd(0.040, 0.0205));
    block2.push_back(Vecd(0.040, 0.020));
    block2.push_back(Vecd(0.0395, 0.020));

    return block2;
}
class Insert : public MultiPolygonShape
{
  public:
    explicit Insert(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(Valve_Shape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(Block1(), ShapeBooleanOps::sub);
    }
};

MultiPolygon createConstrainShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(Block2(), ShapeBooleanOps::add);

    return multi_polygon;
}

int main(int ac, char *av[])
{
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(false);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);

    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.getName())
        : wall_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody insert_body(sph_system, makeShared<Insert>("InsertedBody"));
    /*    insert_body.defineAdaptationRatios(1.15, 2.0);*/
    insert_body.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    insert_body.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? insert_body.generateParticles<BaseParticles, Reload>(insert_body.getName())
        : insert_body.generateParticles<BaseParticles, Lattice>();

    if (sph_system.RunParticleRelaxation())
    {
        InnerRelation insert_body_inner(insert_body);
        InnerRelation wall_boundary_inner(wall_boundary);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(insert_body);
        SimpleDynamics<RandomizeParticlePosition> random_wall_boundary_particles(wall_boundary);
        RelaxationStepInner relaxation_step_inner(insert_body_inner);
        RelaxationStepInner relaxation_step_wall(wall_boundary_inner);
        BodyStatesRecordingToVtp write_insert_body_to_vtp({&insert_body});
        BodyStatesRecordingToVtp write_wall_boundary_to_vtp({&wall_boundary});
        ReloadParticleIO write_particle_reload_files_inner({&insert_body});
        ReloadParticleIO write_particle_reload_files_wall({&wall_boundary});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_insert_body_particles.exec(0.25);
        random_wall_boundary_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        relaxation_step_wall.SurfaceBounding().exec();
        write_insert_body_to_vtp.writeToFile(0);
        write_wall_boundary_to_vtp.writeToFile(0);

        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            relaxation_step_wall.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps  N = " << ite_p << "\n";
                /*write_insert_body_to_vtp.writeToFile(ite_p);*/
                write_wall_boundary_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process  finish !" << std::endl;
        /** Output results. */
        write_particle_reload_files_inner.writeToFile(0);
        write_particle_reload_files_wall.writeToFile(0);
        return 0;
    }

    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation insert_body_inner(insert_body);
    ContactRelation water_block_contact(water_block, RealBodyVector{&wall_boundary, &insert_body});
    ContactRelation insert_body_contact(insert_body, {&water_block});
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);

    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> insert_body_normal_direction(insert_body);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> insert_body_corrected_configuration(insert_body_inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> insert_body_stress_relaxation_first_half(insert_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> insert_body_stress_relaxation_second_half(insert_body_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> insert_body_computing_time_step_size(insert_body);
    BodyRegionByParticle beam_base(insert_body, makeShared<MultiPolygonShape>(createConstrainShape()));
    SimpleDynamics<FixBodyPartConstraint> constraint_beam_base(beam_base);

    //----------------------------------------------------------------------
    // TimeDependentAcceleration time_dependent_acceleration(Vec2d::Zero());
    // SimpleDynamics<GravityForce> apply_gravity_force(water_block, time_dependent_acceleration);
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_indicator(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_block_contact);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);

    BodyAlignedBoxByCell left_disposer(water_block, makeShared<AlignedBoxShape>(Transform(Rotation2d(Pi), Vec2d(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> left_disposer_outflow_deletion(left_disposer, xAxis);
    BodyAlignedBoxByCell right_disposer(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(right_bidirectional_translation)), bidirectional_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> right_disposer_outflow_deletion(right_disposer, xAxis);
    BodyAlignedBoxByCell left_emitter(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::NonPrescribedPressureBidirectionalBuffer left_emitter_inflow_injection(left_emitter, in_outlet_particle_buffer, xAxis);
    BodyAlignedBoxByCell right_emitter(water_block, makeShared<AlignedBoxShape>(Transform(Rotation2d(Pi), Vec2d(right_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<RightInflowPressure> right_emitter_inflow_injection(right_emitter, in_outlet_particle_buffer, xAxis);

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_block_contact);
    SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<RightInflowPressure>> right_inflow_pressure_condition(right_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);

    //----------------------------------------------------------------------
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(insert_body);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> insert_body_update_normal(insert_body);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(insert_body_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(insert_body_contact);

    //----------------------------------------------------------------------
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<int>("Indicator");
    water_block.addBodyStateForRecording<Real>("Density");
    water_block.addBodyStateForRecording<int>("BufferParticleIndicator");
    BodyStatesRecordingToVtp body_states_recording(sph_system.real_bodies_);

    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    boundary_indicator.exec();
    left_emitter_inflow_injection.tag_buffer_particles.exec();
    right_emitter_inflow_injection.tag_buffer_particles.exec();
    wall_boundary_normal_direction.exec();
    insert_body_normal_direction.exec();
    insert_body_corrected_configuration.exec();
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
            //            apply_gravity_force.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_fluid_density.exec();
            viscous_acceleration.exec();
            transport_velocity_correction.exec();

            ///** FSI for viscous force. */
            viscous_force_from_fluid.exec();
            insert_body_update_normal.exec();
            time_instance = TickCount::now();
            size_t inner_ite_dt = 0;
            size_t inner_ite_dt_s = 0;

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                /** Fluid pressure relaxation */
                pressure_relaxation.exec(dt);
                kernel_summation.exec();
                left_inflow_pressure_condition.exec(dt);
                right_inflow_pressure_condition.exec(dt);
                inflow_velocity_condition.exec();
                pressure_force_from_fluid.exec();
                /** Fluid density relaxation */
                density_relaxation.exec(dt);

                inner_ite_dt_s = 0;
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = SMIN(insert_body_computing_time_step_size.exec(), dt);
                    insert_body_stress_relaxation_first_half.exec(dt_s);
                    constraint_beam_base.exec();
                    insert_body_stress_relaxation_second_half.exec(dt_s);
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

            time_instance = TickCount::now();

            left_emitter_inflow_injection.injection.exec();
            right_emitter_inflow_injection.injection.exec();
            left_disposer_outflow_deletion.exec();
            right_disposer_outflow_deletion.exec();
            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_complex.updateConfiguration();
            insert_body.updateCellLinkedList();
            insert_body_contact.updateConfiguration();

            interval_updating_configuration += TickCount::now() - time_instance;

            boundary_indicator.exec();
            left_emitter_inflow_injection.tag_buffer_particles.exec();
            right_emitter_inflow_injection.tag_buffer_particles.exec();
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


