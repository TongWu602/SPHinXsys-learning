/**
 * @file 	test_2d_bifurcation_flow.cpp
 * @brief 	This is the benchmark test for the bifurcation flow.
 * @author 	Dong Wu
 */

/**
 * @brief 	SPHinXsys Library.
 */
#include "sphinxsys.h"
/**
* @brief Namespace cite here.
*/
using namespace SPH;
/**
 * @brief Basic geometry parameters.
 */
Real DL_inlet = 4.0; 							/**< Inlet channel length. */
Real DL_outlet = 8.0; 							/**< Outlet channel length. */
Real DH = 2.0; 									/**< Inlet and outlet channel height. */
Real R = 4.0; 							        /**< Branch radius. */
Real W_upper = 0.8; 							/**< Upper branch width. */
Real W_under = 1.2; 							/**< Under branch width. */
Real N = 180; 				        			/**< Particle number for creating branch. */
Vec2d insert_circle_center = Vec2d(DL_inlet + R + W_upper, 0.0);			/**< Location of the branch center. */
Real resolution_ref = 0.1; 						/**< Global reference resolution. */
Real DL_sponge = resolution_ref * 20.0;	        /**< Sponge region to impose inflow condition. */
Real BW = resolution_ref * 4.0; 			    /**< Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -R - W_under - BW),
	Vec2d(DL_inlet + DL_outlet + 2.0 * (R + W_upper) + BW, R + W_upper + BW));
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;		/**< Density. */
Real U_f = 1.0;			/**< Characteristic velocity. */
Real c_f = 10.0 * U_f;	/**< Speed of sound. */
Real Re = 100.0;		/**< Reynolds number. */
Real mu_f = rho0_f * U_f * DH / Re;	/**< Dynamics viscosity. */

/** A line of measuring points at the entrance of the channel. */
size_t number_observation_pionts = 21;
Real range_of_measure = DH - resolution_ref * 4.0;
Real start_of_measure = resolution_ref * 2.0 - 0.5 * DH;
StdVec<Vec2d> observation_locations()
{
	StdVec<Vec2d> observation_locations;
	for (size_t i = 0; i < number_observation_pionts; ++i)
	{
		observation_locations.push_back(Vec2d(0.0, range_of_measure * Real(i) / Real(number_observation_pionts - 1) + start_of_measure));
	}
	return observation_locations;
}
/**
* @brief define geometry of SPH bodies
*/
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape1()
{
	//geometry
	std::vector<Vecd> water_block_shape1;
	water_block_shape1.push_back(Vecd(-DL_sponge, -0.5 * DH));
	water_block_shape1.push_back(Vecd(-DL_sponge, 0.5 * DH));
	water_block_shape1.push_back(Vecd(DL_inlet + DL_outlet + 2 * (R + W_upper), 0.5 * DH));
	water_block_shape1.push_back(Vecd(DL_inlet + DL_outlet + 2 * (R + W_upper), -0.5 * DH));
	water_block_shape1.push_back(Vecd(-DL_sponge, -0.5 * DH));

	return water_block_shape1;
}
std::vector<Vecd> createWaterBlockShape2()
{
	//geometry
	std::vector<Vecd> water_block_shape2;
	for (int i = 0; i < N + 1; ++i)
	{
		water_block_shape2.push_back(Vecd(insert_circle_center[0] - (R + W_upper) * cos(i * Pi / N),
			insert_circle_center[1] + (R + W_upper) * sin(i * Pi / N)));
	}
	water_block_shape2.push_back(Vecd(DL_inlet, 0.0));

	return water_block_shape2;
}
std::vector<Vecd> createWaterBlockShape3()
{
	//geometry
	std::vector<Vecd> water_block_shape3;
	for (int i = 0; i < N + 1; ++i)
	{
		water_block_shape3.push_back(Vecd(insert_circle_center[0] + (R + W_under) * cos(i * Pi / N),
			insert_circle_center[1] - (R + W_under) * sin(i * Pi / N)));
	}
	water_block_shape3.push_back(Vecd(insert_circle_center[0] + R + W_under, 0.0));

	return water_block_shape3;
}
/** inflow buffer parameters */
Vec2d buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d buffer_translation = Vec2d(-DL_sponge, -0.5 * DH) + buffer_halfsize;
/** create outer wall shape */
std::vector<Vecd> createOuterWallShape1()
{
	//geometry
	std::vector<Vecd> outer_wall_shape1;
	outer_wall_shape1.push_back(Vecd(-DL_sponge - BW, -0.5 * DH - BW));
	outer_wall_shape1.push_back(Vecd(-DL_sponge - BW, 0.5 * DH + BW));
	outer_wall_shape1.push_back(Vecd(DL_inlet + DL_outlet + 2 * (R + W_upper) + BW, 0.5 * DH + BW));
	outer_wall_shape1.push_back(Vecd(DL_inlet + DL_outlet + 2 * (R + W_upper) + BW, -0.5 * DH - BW));
	outer_wall_shape1.push_back(Vecd(-DL_sponge - BW, -0.5 * DH - BW));

	return outer_wall_shape1;
}
std::vector<Vecd> createOuterWallShape2()
{
	//geometry
	std::vector<Vecd> outer_wall_shape2;
	for (int i = 0; i < N + 1; ++i)
	{
		outer_wall_shape2.push_back(Vecd(insert_circle_center[0] - (R + W_upper + BW) * cos(i * Pi / N),
			insert_circle_center[1] + (R + W_upper + BW) * sin(i * Pi / N)));
	}
	outer_wall_shape2.push_back(Vecd(insert_circle_center[0] - (R + W_upper + BW), 0.0));

	return outer_wall_shape2;
}
std::vector<Vecd> createOuterWallShape3()
{
	//geometry
	std::vector<Vecd> outer_wall_shape3;
	for (int i = 0; i < N + 1; ++i)
	{
		outer_wall_shape3.push_back(Vecd(insert_circle_center[0] + (R + W_under + BW) * cos(i * Pi / N),
			insert_circle_center[1] - (R + W_under + BW) * sin(i * Pi / N)));
	}
	outer_wall_shape3.push_back(Vecd(insert_circle_center[0] + R + W_under + BW, 0.0));

	return outer_wall_shape3;
}
/** create inner wall shape */
std::vector<Vecd> createInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH));
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW, 0.5 * DH));
	inner_wall_shape.push_back(Vecd(DL_inlet + DL_outlet + 2 * (R + W_upper) + BW, 0.5 * DH));
	inner_wall_shape.push_back(Vecd(DL_inlet + DL_outlet + 2 * (R + W_upper) + BW, -0.5 * DH));
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH ));

	return inner_wall_shape;
}
 /** Fluid shape definition */
class WaterBlock : public ComplexShape
{
public:
	explicit WaterBlock(const std::string& shape_name) : ComplexShape(shape_name)
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createWaterBlockShape1(), ShapeBooleanOps::add);
		multi_polygon.addAPolygon(createWaterBlockShape2(), ShapeBooleanOps::add);
		multi_polygon.addAPolygon(createWaterBlockShape3(), ShapeBooleanOps::add);
		multi_polygon.addACircle(insert_circle_center, R, 100, ShapeBooleanOps::sub);
		add<MultiPolygonShape>(multi_polygon);
	}
};
/* Definition of the solid body. */
class WallBoundary : public ComplexShape
{
public:
	explicit WallBoundary(const std::string& shape_name) : ComplexShape(shape_name)
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createOuterWallShape1(), ShapeBooleanOps::add);
		multi_polygon.addAPolygon(createOuterWallShape2(), ShapeBooleanOps::add);
		multi_polygon.addAPolygon(createOuterWallShape3(), ShapeBooleanOps::add);
		multi_polygon.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
		multi_polygon.addAPolygon(createWaterBlockShape2(), ShapeBooleanOps::sub);
		multi_polygon.addAPolygon(createWaterBlockShape3(), ShapeBooleanOps::sub);
		multi_polygon.addACircle(insert_circle_center, R, 100, ShapeBooleanOps::add);
		multi_polygon.addACircle(insert_circle_center, R - BW, 100, ShapeBooleanOps::sub);
		add<MultiPolygonShape>(multi_polygon);
	}
};
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
	Real u_ref_, t_ref_;
	AlignedBoxShape& aligned_box_;
	Vecd halfsize_;

	template <class BoundaryConditionType>
	InflowVelocity(BoundaryConditionType& boundary_condition)
		: u_ref_(U_f), t_ref_(20.0),
		aligned_box_(boundary_condition.getAlignedBox()),
		halfsize_(aligned_box_.HalfSize()) {}

	Vecd operator()(Vecd& position, Vecd& velocity)
	{
		Vecd target_velocity = velocity;
		Real run_time = GlobalStaticVariables::physical_time_;
		Real u_ave = u_ref_ * 0.5 * (1.0 + sin(Pi * run_time / t_ref_ - 0.5 * Pi));
		if (aligned_box_.checkInBounds(0, position))
		{
			target_velocity[0] = (-6.0 * position[1] * position[1] / DH / DH + 1.5) * u_ave;
		}
		return target_velocity;
	}
};

/** Main program starts here. */
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up SPHSystem and IO environment.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false);  // Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(true);        // Tag for computation with save particles distribution
	sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<Lattice>();

	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
	wall_boundary.defineAdaptationRatios(1.15, 1.0);
	wall_boundary.defineBodyLevelSetShape()->writeLevelSet(sph_system);
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? wall_boundary.generateParticles<Reload>(wall_boundary.getName())
		: wall_boundary.generateParticles<Lattice>();

	ObserverBody fluid_observer(sph_system, "FluidObserver");
	fluid_observer.generateParticles<Observer>(observation_locations());
	/**
	 * @brief Define body relation map.
	 * The contact map gives the topological connections between the bodies.
	 * Basically the the range of bodies to build neighbor particle lists.
	 */
	InnerRelation water_block_inner(water_block);
	ContactRelation water_block_contact(water_block, { &wall_boundary });
	ContactRelation fluid_observer_contact(fluid_observer, { &water_block });

	/** Run particle relaxation for body-fitted distribution if chosen. */
	if (sph_system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Define body relation map used for particle relaxation.
		//----------------------------------------------------------------------
		InnerRelation wall_boundary_inner(wall_boundary);
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		using namespace relax_dynamics;
		SimpleDynamics<RandomizeParticlePosition> random_wall_particles(wall_boundary);
		RelaxationStepInner relaxation_step_inner(wall_boundary_inner);
		BodyStatesRecordingToVtp write_wall_to_vtp({ &wall_boundary });
		ReloadParticleIO write_particle_reload_files({ &wall_boundary });
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_wall_particles.exec(0.25);
		relaxation_step_inner.SurfaceBounding().exec();
		write_wall_to_vtp.writeToFile(0);
		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the wall N = " << ite_p << "\n";
				write_wall_to_vtp.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of wall finish !" << std::endl;
		/** Output results. */
		write_particle_reload_files.writeToFile(0);
		return 0;
	}
	/**
	 * @brief 	Methods used for time stepping.
	 */
	/** Evaluation of density by summation approach. */
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(water_block_inner, water_block_contact);
	/** Time step size without considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
	/** Pressure relaxation using verlet time stepping. */
	/** Here, we do not use Riemann solver for pressure as the flow is viscous. */
	Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);
	/** viscous acceleration and transport velocity correction can be combined because they are independent dynamics. */
	InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>> transport_correction(water_block_inner, water_block_contact);
	InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);
	/** Computing vorticity in the flow for visualization. */
	InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);
	/** Inflow boundary condition. */
	BodyAlignedBoxByCell inflow_buffer(
		water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(buffer_translation)), buffer_halfsize));
	SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> parabolic_inflow(inflow_buffer);
	/** Periodic BCs in x direction. */
	PeriodicAlongAxis periodic_along_x(water_block.getSPHBodyBounds(), xAxis);
	PeriodicConditionUsingCellLinkedList periodic_condition(water_block, periodic_along_x);

	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	/**
	* @brief Define the methods for I/O operations and observations of the simulation.
	*/
	BodyStatesRecordingToVtp write_real_body_states(sph_system.real_bodies_);
	ObservedQuantityRecording<Vecd> write_fluid_velocity("Velocity", fluid_observer_contact);
	/**
	 * @brief Prepare the simulation with cell linked list, configuration
	 * and case specified initial condition if necessary.
	 */
	/** initialize cell linked lists for all bodies. */
	sph_system.initializeSystemCellLinkedLists();
	/** periodic condition applied after the mesh cell linked list build up
	  * but before the configuration build up. */
	periodic_condition.update_cell_linked_list_.exec();
	/** initialize configurations for all bodies. */
	sph_system.initializeSystemConfigurations();
	/** computing surface normal direction for the wall. */
	wall_boundary_normal_direction.exec();
	/**
	 * @brief Setup for time-stepping control
	 */
	size_t number_of_iterations = 0;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 100.0;			/**< End time. */
	Real D_Time = End_Time / 100.0;	/**< time stamps for output. */
	Real Dt = 0.0;					/**< Default advection time step sizes for fluid. */
	Real dt = 0.0; 					/**< Default acoustic time step sizes for fluid. */
	size_t inner_ite_dt = 0;
	/** Statistics for computing time. */
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	/** First output before the main loop. */
	write_real_body_states.writeToFile();
	write_fluid_velocity.writeToFile(number_of_iterations);
	/**
	 * @brief Main loop starts here.
	 */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			Dt = get_fluid_advection_time_step_size.exec();
			update_density_by_summation.exec();
			viscous_force.exec();
			transport_correction.exec();

			inner_ite_dt = 0;
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				dt = SMIN(get_fluid_time_step_size.exec(), Dt);
				/** Fluid pressure relaxation */
				pressure_relaxation.exec(dt);
				/** Fluid density relaxation */
				density_relaxation.exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				parabolic_inflow.exec();
				inner_ite_dt++;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";
			}
			number_of_iterations++;

			/** Water block configuration and periodic condition. */
			periodic_condition.bounding_.exec();
			water_block.updateCellLinkedList();
			periodic_condition.update_cell_linked_list_.exec();
			water_block_inner.updateConfiguration();
			water_block_contact.updateConfiguration();
		}
		TickCount t2 = TickCount::now();
		/** write run-time observation into file */
		compute_vorticity.exec();
		write_real_body_states.writeToFile();
		fluid_observer_contact.updateConfiguration();
		write_fluid_velocity.writeToFile(number_of_iterations);
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}