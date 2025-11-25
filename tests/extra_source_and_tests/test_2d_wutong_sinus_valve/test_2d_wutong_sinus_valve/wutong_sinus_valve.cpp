/**
 * @file 	mixed_poiseuille_flow.cpp
 * @brief 	2D mixed poiseuille flow example
 * @details This is the one of the basic test cases for mixed pressure/velocity in-/outlet boundary conditions.
 * @author 	Shuoguo Zhang and Xiangyu Hu
 */
//调成脉动速度，出口压力都为零，入口速度给了，流体密度粘性换掉，模量的公式也换掉（查文章）入口压力时间脉动（已经完成），参考脉动泊肃叶流。

#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"

using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real H = 20.0;                                             /**< Channel length. */
Real N = 180.0;                                                //180度
Real resolution_ref = H/500;                             /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;                                /**< Extending width for BCs. */
BoundingBox system_domain_bounds(Vec2d(-BW, - BW), Vec2d(6*H + BW, 2.0 * H + 2*BW));
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real Inlet_pressure = 0.2;
Real Outlet_pressure = 0.0;
Real rho0_f = 1001.0;
Real rho0_s = 890.0;
Real Re = 800.0;
Real mu_f = 0.0043;
Real U_f = pow(0.5 * H, 2.0) * fabs(Inlet_pressure - Outlet_pressure) / (2.0 * mu_f * 6*H);
Real c_f = 10.0 * U_f;

Real poisson = 0.49; /**< Poisson ratio.*/
Real Ae = 1.4e3;
; /**< Normalized Youngs Modulus. */
Real Youngs_modulus = 1.5e6;

Vec2d insert_circle_center = Vec2d(3*H, H);//Location of the sinus center
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d bidirectional_buffer_halfsize = Vec2d(2.5 * resolution_ref, 0.5 * H);
Vec2d left_bidirectional_translation = bidirectional_buffer_halfsize;
Vec2d right_bidirectional_translation = Vec2d(6*H,H) - bidirectional_buffer_halfsize;
Vec2d normal = Vec2d(1.0, 0.0);
//----------------------------------------------------------------------
//	Pressure boundary definition.
//----------------------------------------------------------------------
//struct LeftInflowPressure
//{
//    template <class BoundaryConditionType>
//    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}
//
//    Real operator()(Real &p_)
//    {
//        return p_;
//    }
//};
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

//----------------------------------------------------------------------
//	inflow velocity definition.
////----------------------------------------------------------------------
//struct InflowVelocity
//{
//    Real u_ave;
//
//    template <class BoundaryConditionType>
//    InflowVelocity(BoundaryConditionType &boundary_condition)
//        : u_ave(0.0) {}
//
//    Vecd operator()(Vecd &position, Vecd &velocity)
//    {
//        Vecd target_velocity = Vecd::Zero();
//        Real run_time = GlobalStaticVariables::physical_time_;
//        //这是稳态流
//       /* u_ave = fabs(Inlet_pressure - Outlet_pressure) * (position[1] + 0.5 * DH) * (position[1] + 0.5 * DH - DH) / (2.0 * mu_f * DL) +
//                (4.0 * fabs(Inlet_pressure - Outlet_pressure) * DH * DH) /
//                    (mu_f * DL * Pi * Pi * Pi) * sin(Pi * (position[1] + 0.5 * DH) / DH) * exp(-(Pi * Pi * mu_f * run_time) / (DH * DH));*/
//        
//
//
//        target_velocity[0] = u_ave;
//        target_velocity[1] = 0.0;
//
//        return target_velocity;
//    }
//};
struct InflowVelocity
{
    Real u_ref_; // 参考速度

    // 构造函数
    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(0.1) {}

    Real period = 1.0; // 设置周期为1秒

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = Vecd::Zero();
        Real run_time = GlobalStaticVariables::physical_time_;

        // 计算周期性平均速度
        Real time_in_period = fmod(run_time, period);

        // 根据多项式形式计算平均速度
        target_velocity[0] = 
            -6976.27178067763 * pow(time_in_period, 12) + 
            39718.2651305609 * pow(time_in_period, 11) -
            97231.4260051807 * pow(time_in_period, 10) + 
            133568.920608098 * pow(time_in_period, 9) - 
            112749.851450734 * pow(time_in_period, 8) + 
            60318.9917716997 * pow(time_in_period, 7) - 
            20468.105911115 * pow(time_in_period, 6) + 
            4332.72581256502 * pow(time_in_period, 5) - 
            536.171648078612 * pow(time_in_period, 4) + 
            17.972203531116 * pow(time_in_period, 3) + 
            4.77916722403468 * pow(time_in_period, 2) + 
            0.171368158347297 * time_in_period + 
            0.000354679355236548;

        target_velocity[1] = 0.0; // y 方向速度为零
        return target_velocity;
    }
};

//----------------------------------------------------------------------
//	Fluid body definition.
//----------------------------------------------------------------------
//小长方体
std::vector<Vecd> Water_Block_Shape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(0.0, 0.0));
    water_block_shape.push_back(Vecd(0.0, H));
    for (size_t i = 0; i < N + 1; ++i)
    {
        water_block_shape.push_back(Vecd(insert_circle_center[0] - H * cos(i * Pi / N),
                                          insert_circle_center[1] + H * sin(i * Pi / N)));
    }
    water_block_shape.push_back(Vecd(6 * H, H));
    water_block_shape.push_back(Vecd(6 * H, 0.0));
    water_block_shape.push_back(Vecd(0.0, 0.0));
    return water_block_shape;
}




//现在定义壁面

std::vector<Vecd> Wall_Block_Shape()
{
    std::vector<Vecd> wall_block_shape;
    wall_block_shape.push_back(Vecd(0.0, -BW));
    wall_block_shape.push_back(Vecd(0.0, H + BW));
    wall_block_shape.push_back(Vecd(2*H-BW, H + BW));
    for (size_t i = 1; i < N; ++i)
    {
        wall_block_shape.push_back(Vecd(insert_circle_center[0] - (H + BW) * cos(i * Pi / N),
                                         insert_circle_center[1] + (H + BW) * sin(i * Pi / N)));
    }
    wall_block_shape.push_back(Vecd(4 * H + BW, H+BW));
    wall_block_shape.push_back(Vecd(6 * H, H + BW));
    wall_block_shape.push_back(Vecd(6 * H, -BW));
    wall_block_shape.push_back(Vecd(0.0, -BW));
    return wall_block_shape;
}


// 小弹性瓣

std::vector<Vecd> Valve_Shape()
{
    std::vector<Vecd> valve_shape3;
    valve_shape3.push_back(Vecd(39.7708, 20.3422)); // 点1

    valve_shape3.push_back(Vecd(58.4386, 1.6746)); // 点2
    valve_shape3.push_back(Vecd(58.3254, 1.5559));        // 3
    valve_shape3.push_back(Vecd(39.6576, 20.2291));                                 // 4
    valve_shape3.push_back(Vecd(39.7708, 20.3422));       // 5
                                 
    return valve_shape3;
}
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {

        multi_polygon_.addAPolygon(Water_Block_Shape(), ShapeBooleanOps::add);
    
        multi_polygon_.addAPolygon(Valve_Shape(), ShapeBooleanOps::sub);
    }
};
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(Wall_Block_Shape(), ShapeBooleanOps::add);
        
        multi_polygon_.addAPolygon(Water_Block_Shape(), ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(Valve_Shape(), ShapeBooleanOps::sub);


    }
};


//块1
std::vector<Vecd> Block1()
{
    std::vector<Vecd> block1;
    block1.push_back(Vecd(34.8406, 20.16)); // 点1

    block1.push_back(Vecd(34.8406, 25.16)); // 点2
    block1.push_back(Vecd(39.8406, 25.16)); // 3
    block1.push_back(Vecd(39.8406, 20.16));    
    block1.push_back(Vecd(34.8406, 20.16));
  
    return block1;

}

// 块2
std::vector<Vecd> Block2()
{
    std::vector<Vecd> block2;
    block2.push_back(Vecd(39.5, 20)); // 点1

    block2.push_back(Vecd(39.5, 20.5)); // 点2
    block2.push_back(Vecd(40, 20.5)); // 3
    block2.push_back(Vecd(40, 20));
    block2.push_back(Vecd(39.5, 20));

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
 

 /** create constrain shape. */
 MultiPolygon createConstrainShape()
 {
     MultiPolygon multi_polygon;
     multi_polygon.addAPolygon(Block2(), ShapeBooleanOps::add);

     return multi_polygon;
 }

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false); // Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(false);       // Tag for computation with save particles distribution
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
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
    insert_body.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? insert_body.generateParticles<BaseParticles, Reload>(insert_body.getName())
        : insert_body.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        // ----------------------------------------------------------------------
        // Define body relation map used for particle relaxation.
        // ----------------------------------------------------------------------
        InnerRelation insert_body_inner(insert_body);
        InnerRelation wall_boundary_inner(wall_boundary);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(insert_body);
        RelaxationStepInner relaxation_step_inner(insert_body_inner);
        BodyStatesRecordingToVtp write_insert_body_to_vtp({&insert_body});
        ReloadParticleIO write_particle_reload_files_inner({&insert_body});
        SimpleDynamics<RandomizeParticlePosition> random_wall_boundary_particles(wall_boundary);
        RelaxationStepInner relaxation_step_wall(wall_boundary_inner);
        BodyStatesRecordingToVtp write_wall_boundary_to_vtp({&wall_boundary});
        ReloadParticleIO write_particle_reload_files_wall({&wall_boundary});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_insert_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_insert_body_to_vtp.writeToFile(0);
        random_wall_boundary_particles.exec(0.25);
        relaxation_step_wall.SurfaceBounding().exec();
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
                write_insert_body_to_vtp.writeToFile(ite_p);
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
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation insert_body_inner(insert_body);
    ContactRelation water_block_contact(water_block, RealBodyVector{&wall_boundary, &insert_body});
    ContactRelation insert_body_contact(insert_body, {&water_block});
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Algorithms of solid dynamics.
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
    //SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_disposer); // 这里可能要改成emitter

    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(insert_body);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> insert_body_update_normal(insert_body);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(insert_body_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(insert_body_contact);


    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
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
    left_emitter_inflow_injection.tag_buffer_particles.exec();
    right_emitter_inflow_injection.tag_buffer_particles.exec();
    wall_boundary_normal_direction.exec();
    insert_body_normal_direction.exec();
    insert_body_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 100.0;   /**< End time. */
    Real Output_Time = end_time/1000.0; /**< Time stamps for output of body states. */

    TickCount t1 = TickCount::now();
    TimeInterval interval;
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
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_fluid_density.exec();
            viscous_acceleration.exec();
            transport_velocity_correction.exec();

            /** FSI for viscous force. */
            viscous_force_from_fluid.exec();
            /** Update normal direction on elastic body.*/
            insert_body_update_normal.exec();


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
                /*inflow_velocity_condition.exec();*/
                /** FSI for pressure force. */
                pressure_force_from_fluid.exec();
                /** Fluid density relaxation */
                density_relaxation.exec(dt);

                
                

                /** Solid dynamics. */
                inner_ite_dt_s = 0;
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = SMIN(insert_body_computing_time_step_size.exec(), dt - dt_s_sum);
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

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "	dt / dt_s = " << inner_ite_dt_s << "\n";
            }
            number_of_iterations++;


            left_emitter_inflow_injection.injection.exec();//对应periodic_condition.bounding_.exec();
            right_emitter_inflow_injection.injection.exec();
            left_disposer_outflow_deletion.exec();
            right_disposer_outflow_deletion.exec();
            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_complex.updateConfiguration();
            insert_body.updateCellLinkedList();
            insert_body_contact.updateConfiguration();
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
