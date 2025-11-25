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
Real resolution_ref = H / 80.0;
Real BW = resolution_ref * 4.0;
BoundingBox system_domain_bounds(Vec2d(0.0, -BW), Vec2d(6 * H, 2.0 * H + 2.0 * BW));
Vec2d gravity(0.0, -1.0);
//----------------------------------------------------------------------
Real gravity_g = 9.81;
Real Outlet_pressure = 0.0;
Real rho0_f = 1090.0;
Real rho0_s = 1200.0;
Real Re = 800.0;
Real mu_f = 0.00436;
Real U_f = 0.35;
Real c_f = 10.0 * U_f;

Real poisson = 0.49;
Real Youngs_modulus = 1.5e6;
Real round_circle_radius = 0.004;
Vec2d insert_circle_center = Vec2d(3.0 * H, H);
Vec2d round_circle_center = Vec2d(0.08366, 0.024);
Vec2d round_circle_center2 = Vec2d(0.084, 0.0244);

//----------------------------------------------------------------------
Vec2d bidirectional_buffer_halfsize = Vec2d(2.5 * resolution_ref, 0.5 * H);
Vec2d left_bidirectional_translation = bidirectional_buffer_halfsize;
Vec2d right_bidirectional_translation = Vec2d(6 * H, H) - bidirectional_buffer_halfsize;
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
    Real u_avg_, u_amp_, t_p_, k;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_avg_(0.04), u_amp_(0.11), t_p_(2.45) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = Vecd::Zero();
        Real run_time = GlobalStaticVariables::physical_time_;

        // 计算 t/Tp
        Real time_in_period = fmod(run_time, t_p_) / t_p_;

        // 根据时间段选择不同的表达式
        if (time_in_period <= 0.37)
        {
            target_velocity[0] = u_avg_ + u_amp_ * sin(2.0 * M_PI * time_in_period / 0.74);
        }
        else if (time_in_period > 0.37 && time_in_period <= 1.0)
        {
            target_velocity[0] = u_avg_ + (u_amp_ / 2.0) * sin(2.0 * M_PI * (time_in_period + 0.26) / 1.26);
        }
        if (run_time <= 2.45)
        {
            k = 0.5 - 0.5 * cos(M_PI * run_time / 2.45);
        }
        else
        {
            k = 1;
        }
        // 假设流体沿X方向流动
        target_velocity[0] *= k;

        // Y方向的速度为0
        target_velocity[1] = 0.0;

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
        size_t number_observation_points = 11;
        Real range_of_measure = H ;
        Real start_of_measure = 0;
        /** the measuring locations */
        for (size_t i = 0; i < number_observation_points; ++i)
        {
            Vec2d point_coordinate(0.10, range_of_measure * (Real)i / (Real)(number_observation_points - 1) + start_of_measure);
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
        size_t number_observation_points = 11;
        Real range_of_measure = H;
        Real start_of_measure = 0;
        /** the measuring locations */
        for (size_t i = 0; i < number_observation_points; ++i)
        {
            Vec2d point_coordinate(0.11, range_of_measure * (Real)i / (Real)(number_observation_points - 1) + start_of_measure);
            positions_.push_back(point_coordinate);
        }
    }
};

//----------------------------------------------------------------------
std::vector<Vecd> Water_Block_Shape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(0.0, 0.0));
    water_block_shape.push_back(Vecd(0.0, H));
    for (size_t i = 0; i < N - 9.0; ++i)
    {
        water_block_shape.push_back(Vecd(insert_circle_center[0] - H * cos(i * Pi / N),
                                         insert_circle_center[1] + H * sin(i * Pi / N)));
    }
    water_block_shape.push_back(Vecd(0.08366, H));
    water_block_shape.push_back(Vecd(6.0 * H, H));
    water_block_shape.push_back(Vecd(6.0 * H, 0.0));
    water_block_shape.push_back(Vecd(0.0, 0.0));

    return water_block_shape;
}

std::vector<Vecd> Wall_Block_Shape()
{
    std::vector<Vecd> wall_block_shape;
    wall_block_shape.push_back(Vecd(0.0, -BW));
    wall_block_shape.push_back(Vecd(0.0, H + BW));
    wall_block_shape.push_back(Vecd(2 * H - BW, H + BW));
    for (size_t i = 1; i < N - 10.0; ++i)
    {
        wall_block_shape.push_back(Vecd(insert_circle_center[0] - (H + BW) * cos(i * Pi / N),
                                        insert_circle_center[1] + (H + BW) * sin(i * Pi / N)));
    }
    wall_block_shape.push_back(Vecd(0.08366 + BW, H + BW));
    wall_block_shape.push_back(Vecd(6.0 * H, H + BW));
    wall_block_shape.push_back(Vecd(6.0 * H, -BW));
    wall_block_shape.push_back(Vecd(0.0, -BW));

    return wall_block_shape;
}

std::vector<Vecd> Valve_Shape()
{
    std::vector<Vecd> valve_shape3;
    valve_shape3.push_back(Vecd(0.0404631, 0.0201885));
    valve_shape3.push_back(Vecd(0.048683, 0.0));
    valve_shape3.push_back(Vecd(0.0476033, 0.0));
    valve_shape3.push_back(Vecd(0.0395369, 0.0198115));
    valve_shape3.push_back(Vecd(0.0404631, 0.0201885));

    return valve_shape3;
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
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(Water_Block_Shape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(Valve_Shape(), ShapeBooleanOps::sub);
        multi_polygon_.addACircle(round_circle_center, round_circle_radius, 100, ShapeBooleanOps::sub);
        multi_polygon_.addACircle(Vecd(0.04, 0.02), 0.0005, 100, ShapeBooleanOps::sub); // 剪掉一个小关节圆
    }
};
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(Wall_Block_Shape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(Water_Block_Shape(), ShapeBooleanOps::sub);
        multi_polygon_.addACircle(round_circle_center, round_circle_radius, 100, ShapeBooleanOps::add); // 修圆角
        multi_polygon_.addAPolygon(Valve_Shape(), ShapeBooleanOps::sub);
        // multi_polygon_.addACircle(round_circle_center2, round_circle_radius, 100, ShapeBooleanOps::sub);//修圆角
        // multi_polygon_.addAPolygon(Wall_Block_Shape1(), ShapeBooleanOps::sub);//修圆角
        multi_polygon_.addACircle(Vecd(0.04, 0.02), 0.0005, 100, ShapeBooleanOps::sub); // 剪掉一个小关节圆
    }
};

class Insert : public MultiPolygonShape
{
  public:
    explicit Insert(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(Valve_Shape(), ShapeBooleanOps::add);
        multi_polygon_.addACircle(Vecd(0.04, 0.02), 0.0005, 100, ShapeBooleanOps::add); // 加上一个小关节圆
    }
};

MultiPolygon createInsertShape(SPHBody &sph_body)
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(Valve_Shape(), ShapeBooleanOps::add);
    multi_polygon.addACircle(Vecd(0.04, 0.02), 0.0005, 100, ShapeBooleanOps::add); // 加上一个小关节圆

    return multi_polygon;
}
MultiPolygon createwallShape(SPHBody &sph_body)
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(Wall_Block_Shape(), ShapeBooleanOps::add);
    multi_polygon.addAPolygon(Water_Block_Shape(), ShapeBooleanOps::sub);
    multi_polygon.addACircle(round_circle_center, round_circle_radius, 100, ShapeBooleanOps::add); // 修圆角
    multi_polygon.addAPolygon(Valve_Shape(), ShapeBooleanOps::sub);

    multi_polygon.addACircle(Vecd(0.04, 0.02), 0.0005, 100, ShapeBooleanOps::sub); // 剪掉一个小关节圆

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
    insert_body.defineAdaptationRatios(1.15, 1.0);
    insert_body.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    insert_body.defineMaterial<Solid>(rho0_s);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? insert_body.generateParticles<BaseParticles, Reload>(insert_body.getName())
        : insert_body.generateParticles<BaseParticles, Lattice>();

    ObserverBody beam_observer(sph_system, "BeamObserver");
    StdVec<Vecd> beam_observation_location = {Vecd(0.048125, 0.000125)};
    beam_observer.generateParticles<BaseParticles, Observer>(beam_observation_location);

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<BaseParticles, FluidObserver>();

    ObserverBody fluid_observer2(sph_system, "FluidObserver2");
    fluid_observer2.generateParticles<BaseParticles, FluidObserver2>();


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
    ContactRelation insert_body_contact1(insert_body, {&water_block});
    /*SurfaceContactRelation insert_body_contact2(insert_body, {&wall_boundary});*/
    ContactRelation beam_observer_contact(beam_observer, {&insert_body});
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    ContactRelation fluid_observer_contact2(fluid_observer2, {&water_block});
    //----------------------------------------------------------------------
    //	Algorithms of Fluid.
    //----------------------------------------------------------------------

    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> insert_body_normal_direction(insert_body);
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact); // 这是干什么的？
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
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    /*InteractionWithUpdate<LinearGradientCorrectionMatrixInner> ball_corrected_configuration(insert_body_inner);*/
    /*SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> insert_body_update_normal(insert_body);*/ /*加上这一行就会出问题说没有梯度*/
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(insert_body_contact1);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(insert_body_contact1);
    Gravity constant_gravity(gravity);
    SimpleDynamics<GravityForce> insert_constant_gravity(insert_body, constant_gravity);
    /*InteractionDynamics<solid_dynamics::ShellContactDensity> insertbody_update_contact_density(insert_body_contact2);
    InteractionWithUpdate<solid_dynamics::ContactForceFromWall> insertbody_compute_solid_contact_forces(insert_body_contact2);*/
    /*InteractionWithUpdate<LinearGradientCorrectionMatrixInner> insert_body_corrected_configuration(insert_body_inner);*/
    // Dynamics1Level<solid_dynamics::Integration1stHalfPK2> insert_body_stress_relaxation_first_half(insert_body_inner);
    // Dynamics1Level<solid_dynamics::Integration2ndHalf> insert_body_stress_relaxation_second_half(insert_body_inner);
    /*ReduceDynamics<solid_dynamics::AcousticTimeStepSize> insert_body_computing_time_step_size(insert_body);*/
    /*   BodyRegionByParticle beam_base(insert_body, makeShared<MultiPolygonShape>(createConstrainShape()));*/
    // SimpleDynamics<FixBodyPartConstraint> constraint_beam_base(beam_base);
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<int>("Indicator");
    water_block.addBodyStateForRecording<Real>("Density");
    water_block.addBodyStateForRecording<int>("BufferParticleIndicator");

    //----------------------------------------------------------------------
    /*solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(insert_body);*/
    /*SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> insert_body_update_normal(insert_body);*/

    //----------------------------------------------------------------------
    //	Algorithms of Simbody.
    //----------------------------------------------------------------------
    SimTK::MultibodySystem MBsystem;
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    SimTK::GeneralForceSubsystem forces(MBsystem);
    SimTK::CableTrackerSubsystem cables(MBsystem);

    SolidBodyPartForSimbody insert_multibody(insert_body, makeShared<MultiPolygonShape>(createInsertShape(insert_body), "InsertShape"));
    SimTK::Body::Rigid pin_spot_info(*insert_multibody.body_part_mass_properties_);
    /*SolidBodyPartForSimbody wall_multibody(wall_boundary, makeShared<MultiPolygonShape>(createwallShape(wall_boundary), "wallShape"));
    SimTK::Body::Rigid fix_spot_info(*wall_multibody.body_part_mass_properties_);*/

    SimTK::MobilizedBody::Pin pin_spot(matter.Ground(), SimTK::Transform(SimTK::Vec3(0.04, 0.02, 0.0)), pin_spot_info, SimTK::Transform(SimTK::Vec3(0.0, 0.0, 0.0)));
    // SimTK::MobilizedBody::Weld fix_spot(matter.Ground(), SimTK::Transform(SimTK::Vec3(0.0, 0.0)), fix_spot_info, SimTK::Transform(SimTK::Vec3(0.0, 0.0)));

    SimTK::MobilizerQIndex whichQ(0);
    /*SimTK::Force::MobilityLinearDamper linear_damper(forces, pin_spot, whichQ, 0.00001);*/
    SimTK::Force::MobilityLinearStop myStop(forces, pin_spot, whichQ, 10, 10, 0, Pi); // 留个缝，该出入口压力，改kd
    /*SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTK::Vec3(0.0, Real(-100000), 0.0), 0.0);*/
    SimTK::Force::DiscreteForces force_on_bodies(forces, matter);

    /*SimTK::Visualizer viz(MBsystem);
    MBsystem.addEventReporter(new SimTK::Visualizer::Reporter(viz, 0.1));*/
    SimTK::State state = MBsystem.realizeTopology();
    /*viz.report(state);*/

    /* std::cout << "Hit ENTER to run a short simulation ...";
     getchar();*/
    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);
    ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody> force_on_spot_flap(insert_multibody, MBsystem, pin_spot, integ);
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody> constraint_spot_flap(insert_multibody, MBsystem, pin_spot, integ);

    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_beam_tip_displacement("Position", beam_observer_contact);
    ObservedQuantityRecording<Vecd> write_fluid_velocity("Velocity", fluid_observer_contact);
    ObservedQuantityRecording<Vecd> write_fluid_velocity2("Velocity", fluid_observer_contact2);
    //----------------------------------------------------------------------
    //	Algorithms of data output.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system.real_bodies_);
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    insert_body_normal_direction.exec();
    boundary_indicator.exec();
    insert_constant_gravity.exec();
    left_emitter_inflow_injection.tag_buffer_particles.exec();
    right_emitter_inflow_injection.tag_buffer_particles.exec();

    /* ball_corrected_configuration.exec();*/
    /*insert_body_corrected_configuration.exec();*/
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 25.0;    /**< End time. */
    Real Output_Time = 0.01; /**< Time stamps for output of body states. */
    Real dt = 0.0;

    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    body_states_recording.writeToFile();
    write_beam_tip_displacement.writeToFile(number_of_iterations);
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
            /*insert_body_update_normal.exec();*/
            time_instance = TickCount::now();
            size_t inner_ite_dt = 0;
            /*size_t inner_ite_dt_s = 0;*/

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                /** Fluid pressure relaxation */
                pressure_relaxation.exec(dt);
                pressure_force_from_fluid.exec();
                /*insertbody_update_contact_density.exec();
                insertbody_compute_solid_contact_forces.exec();*/
                kernel_summation.exec();
                left_inflow_pressure_condition.exec(dt);
                right_inflow_pressure_condition.exec(dt);
                inflow_velocity_condition.exec();

                /** Fluid density relaxation */
                density_relaxation.exec(dt);
                integ.stepBy(dt);
                /*inner_ite_dt_s = 0;*/
                /*Real dt_s_sum = 0.0;*/
                /*average_velocity_and_acceleration.initialize_displacement_.exec();*/
                // while (dt_s_sum < dt)
                //{
                //     Real dt_s = SMIN(insert_body_computing_time_step_size.exec(), dt);
                //     insert_body_stress_relaxation_first_half.exec(dt_s);
                //     /*constraint_beam_base.exec();*/
                //     insert_body_stress_relaxation_second_half.exec(dt_s);
                //     dt_s_sum += dt_s;
                //     inner_ite_dt_s++;
                // }
                // average_velocity_and_acceleration.update_averages_.exec(dt);

                SimTK::State &state_for_update = integ.updAdvancedState();
                force_on_bodies.clearAllBodyForces(state_for_update);
                force_on_bodies.setOneBodyForce(state_for_update, pin_spot, force_on_spot_flap.exec());

                constraint_spot_flap.exec();
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
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";
            }
            number_of_iterations++;

            time_instance = TickCount::now();

            left_emitter_inflow_injection.injection.exec();
            right_emitter_inflow_injection.injection.exec();
            left_disposer_outflow_deletion.exec();
            right_disposer_outflow_deletion.exec();
            water_block.updateCellLinkedListWithParticleSort(100);
            insert_body.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            insert_body_contact1.updateConfiguration();
            /*insert_body_contact2.updateConfiguration();*/
            /*insert_body_contact2.updateConfiguration();*/
            interval_updating_configuration += TickCount::now() - time_instance;

            boundary_indicator.exec();
            left_emitter_inflow_injection.tag_buffer_particles.exec();
            right_emitter_inflow_injection.tag_buffer_particles.exec();
            write_beam_tip_displacement.writeToFile(number_of_iterations);
        }
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        write_fluid_velocity.writeToFile(number_of_iterations);
        write_fluid_velocity2.writeToFile(number_of_iterations);
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