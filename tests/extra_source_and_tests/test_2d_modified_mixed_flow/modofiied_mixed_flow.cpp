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
Real DL = 4.0;                                             /**< Channel length. */
Real DH = 1.0;                                             /**< Channel height. */
Real resolution_ref = DH / 20.0;                             /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;                                /**< Extending width for BCs. */
BoundingBox system_domain_bounds(Vec2d(-BW, -DH - BW), Vec2d(DL + BW, 2.0 * DH + BW));
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real Inlet_pressure = 0.2;
Real Outlet_pressure = 0.0;
Real rho0_f = 1060.0;
Real Re = 50.0;
Real mu_f = 0.0035;
Real U_f = pow(0.5 * DH, 2.0) * fabs(Inlet_pressure - Outlet_pressure) / (2.0 * mu_f * DL);
Real c_f = 10.0 * U_f;
Real rho0_s = 1000.0; /**< Reference density.*/
Real poisson = 0.4; /**< Poisson ratio.*/
Real Ae = 1.4e3;
; /**< Normalized Youngs Modulus. */
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d bidirectional_buffer_halfsize = Vec2d(2.5 * resolution_ref, 0.5 * DH);
Vec2d left_bidirectional_translation = bidirectional_buffer_halfsize;
Vec2d right_bidirectional_translation = Vec2d(DL, DH) - bidirectional_buffer_halfsize;
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
struct LeftInflowPressure//脉动压力
{
    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        /*pulsatile pressure*/
        Real pressure = Inlet_pressure * cos(3*Pi * GlobalStaticVariables::physical_time_);
        /*constant pressure*/
        //        Real pressure = Inlet_pressure;
        return pressure;
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
//struct InflowVelocity
//{
//    Real u_ref_, t_ref_, u_ave;
//
//    template <class BoundaryConditionType>
//    InflowVelocity(BoundaryConditionType &boundary_condition)
//        : u_ref_(0.1), t_ref_(0.1), u_ave(0.0) {}
//    Real period = 2.5;
//    Vecd operator()(Vecd &position, Vecd &velocity)
//    {
//        Vecd target_velocity = Vecd::Zero();
//        Real run_time = GlobalStaticVariables::physical_time_;
//        if (run_time < t_ref_)
//        {
//            u_ave = u_ref_;
//        }
//        else if (((run_time - t_ref_) - int((run_time - t_ref_) / period) * period) >= 0.0 && ((run_time - t_ref_) - int((run_time - t_ref_) / period) * period) <= 2*0.218)
//        {
//            u_ave = 0.5 * sin(4.0 * Pi * ((run_time - t_ref_) + 0.0160236));
//        }
//        else
//        {
//            u_ave = u_ref_;
//        }
//
//        target_velocity[0] = u_ave;
//        target_velocity[1] = 0.0;
//        return target_velocity;
//    }
//};
//----------------------------------------------------------------------
//	Fluid body definition.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> water_block_shape;
        water_block_shape.push_back(Vecd(0.0, 0.0));
        water_block_shape.push_back(Vecd(0.0, DH));
        water_block_shape.push_back(Vecd(DL, DH));
        water_block_shape.push_back(Vecd(DL, 0.0));
        water_block_shape.push_back(Vecd(0.0, 0.0));
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
    }
};

//----------------------------------------------------------------------
//	Wall boundary body definition.
//----------------------------------------------------------------------
std::vector<Vecd> outer_wall_shape
{
    Vecd(0.0, -BW),Vecd(0.0, DH + BW),Vecd(DL, DH + BW),Vecd(DL, -BW),Vecd(0.0, -BW)
};
 
 std::vector<Vecd> inner_wall_shape
 {
    Vecd(-BW, 0.0),Vecd(-BW, DH),Vecd(DL + BW, DH),Vecd(DL + BW, 0.0),Vecd(-BW, 0.0)
};
 //约束部分几何
 std::vector<Vecd> wall1_shape{
     Vecd(0.0, -2.0 * BW), Vecd(0.0, DH + 2 * BW), Vecd(5 * resolution_ref, DH + 2 * BW), Vecd(5 * resolution_ref, -2 * BW), Vecd(0.0, -2.0 * BW)};
 std::vector<Vecd> wall2_shape{
     Vecd(DL - 5 * resolution_ref, -2 * BW), Vecd(DL - 5 * resolution_ref, DH + 2 * BW), Vecd(DL, DH + 2 * BW), Vecd(DL, -2 * BW), Vecd(DL - 5 * resolution_ref, -2 * BW)};

 // 弹性部分与连接体
 class WallBoundary : public MultiPolygonShape
 {
   public:
     explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
     {
         multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
         multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
     }
 };
 /** create 头尾 as constrain shape. */
 MultiPolygon createConstrainShape()
 {
     MultiPolygon multi_polygon;
     multi_polygon.addAPolygon(wall1_shape, ShapeBooleanOps::add);
     multi_polygon.addAPolygon(wall2_shape, ShapeBooleanOps::add);
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
    wall_boundary.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    wall_boundary.generateParticles<BaseParticles, Lattice>();
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
    //----------------------------------------------------------------------
    //	Algorithms of solid dynamics.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> wall_boundary_corrected_configuration(wall_boundary_inner); // 不太清楚作用
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> wall_boundary_stress_relaxation_first_half(wall_boundary_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> wall_boundary_stress_relaxation_second_half(wall_boundary_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> wall_boundary_computing_time_step_size(wall_boundary);
    BodyRegionByParticle constraint_in_outlet(wall_boundary, makeShared<MultiPolygonShape>(createConstrainShape()));
    SimpleDynamics<FixBodyPartConstraint> constraint_in_outlet_base(constraint_in_outlet);
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
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_disposer); // 这里可能要改成emitter

    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(wall_boundary);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> wall_boundary_update_normal(wall_boundary);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(wall_boundary_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(wall_boundary_contact);

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
    wall_boundary_corrected_configuration.exec();
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
                    Real dt_s = SMIN(wall_boundary_computing_time_step_size.exec(), dt - dt_s_sum);
                    wall_boundary_stress_relaxation_first_half.exec(dt_s);
                    constraint_in_outlet_base.exec();
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
            wall_boundary.updateCellLinkedList();
            wall_boundary_contact.updateConfiguration();
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
