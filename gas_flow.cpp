// Gas Flow Numerical Model Based on Molecular Flow Conditions

// Blake Leonard
// 2012


#include <iostream>
using namespace::std;

#include <fstream>
using namespace::std;

#include <cmath>
using namespace::std;

#include <ctime>
using namespace::std;


void main()
{


	const double tube_diameter = 0.000025;              // Meters

	const double tube_length = 0.381;                   // Meters

	const double inlet_pressure = 1;               // Pascals ( 1 atm = 101325 )

	const int inlet_temperature = 300;                  // Kelvin

	const int molecular_mass = 4;                       // AMU

	const double molecular_diameter = 0.000000000062;   // Meters

	const double time_step = 0.00001;                  // Seconds


 
	const unsigned int global_steps = 4000000000;

	const unsigned int max_particles = 40000000;	    // Max particles allowed in tube

	const _int64 particle_weight = 1000;           // # of molecules represented by virtual particles

	const int local_pressure_resolution = 20;		    // # of local pressure cells within tube

	const int sampling_resolution = 100;              // # of global_steps used to sample flow rate

	int particles = 0;								    // # of particles currently in tube

	int old_particles = 0;

	const double k = 0.0000000000000000000000138;       // Boltzmann's Constant in Joules/Kelvin



	const double molecular_mass_kg = molecular_mass * 0.0000000000000000000000000016605;

	const double mean_molecular_speed = sqrt ( ( 8 * k * inlet_temperature ) / ( 3.141592 * molecular_mass_kg ) ); 

	const double local_space_size = tube_length / local_pressure_resolution;

	const double local_area = ( ( 3.141592 * tube_diameter * tube_diameter ) / 4 );

	const double local_volume = local_area * local_space_size;

	double inlet_conductance;        // = ( mean_molecular_speed * local_area ) / 4;           // Molecular Flow Conductance, Change?

	double pressure_ratio;

	double run_time = 0;


	double inlet_flow_rate;

	double inlet_volume_rate;

	double inlet_molecule_rate;

	double inlet_particle_rate;


	int long_outlet_particle_rate = 0;

	int long_inlet_particle_rate = 0;

	double long_outlet_molecule_rate;

	double long_inlet_molecule_rate;

	double long_outlet_volume_rate;

	double long_inlet_volume_rate;


	double particle_r_position;

	double particle_theta_position;

	double absolute_velocity;
	
	
	static double particle_x_position[max_particles];

	static double particle_y_position[max_particles];

	static double particle_z_position[max_particles];


	static double particle_x_velocity[max_particles];

	static double particle_y_velocity[max_particles];

	static double particle_z_velocity[max_particles];



	double local_pressure[local_pressure_resolution];     // Pascals

	double average_first_pressure = 0;

	int local_particle_count[local_pressure_resolution];

	static double local_real_molecule_count[local_pressure_resolution];

	double collision_probability[local_pressure_resolution];


	srand((unsigned)time(0));



	// Global Iterations


	for ( int istep = 0; istep < global_steps; istep ++ )
	{

		old_particles = particles;


		// Initialize local particle counts


		for ( int kstep = 0; kstep < local_pressure_resolution; kstep ++ )
		{

			local_particle_count[kstep] = 0;
	
		}

		
		// Count Particles in local pressures regions


		for ( int jstep = 0; jstep < particles; jstep ++ )
		{

			for ( int lstep = 0; lstep < local_pressure_resolution; lstep ++ )
			{

				if ( ( particle_z_position[jstep] >= ( local_space_size * lstep ) ) && ( particle_z_position[jstep] < ( local_space_size * ( lstep + 1 ) ) ) )
				{

					local_particle_count[lstep] ++;

					break;

				}
				
			}

		}


		// Calculate local pressures & collision probabilities based on particle densities


		for ( jstep = 0; jstep < local_pressure_resolution; jstep ++ )
		{

			local_real_molecule_count[jstep] = local_particle_count[jstep] * particle_weight;

			local_pressure[jstep] = ( local_real_molecule_count[jstep] * k * inlet_temperature ) / local_volume;

			collision_probability[jstep] = ( 4.442883 * molecular_diameter * molecular_diameter * local_pressure[jstep] * mean_molecular_speed * time_step ) / ( k * inlet_temperature );

		}


		average_first_pressure = average_first_pressure + local_pressure[0];


		// Move particles, handle collisions and exiting particles


		for ( jstep = 0; jstep < particles; jstep ++ )
		{

			particle_x_position[jstep] = particle_x_position[jstep] + particle_x_velocity[jstep] * time_step;

			particle_y_position[jstep] = particle_y_position[jstep] + particle_y_velocity[jstep] * time_step;

			particle_z_position[jstep] = particle_z_position[jstep] + particle_z_velocity[jstep] * time_step;


			// Collision with wall

			particle_r_position = sqrt( ( particle_x_position[jstep] * particle_x_position[jstep] ) + ( particle_y_position[jstep] * particle_y_position[jstep] ) );


			if ( particle_r_position > ( tube_diameter / 2 ) )
			{

				particle_x_position[jstep] = ( particle_x_position[jstep] / particle_r_position ) * ( tube_diameter / 2 );

				particle_y_position[jstep] = ( particle_y_position[jstep] / particle_r_position ) * ( tube_diameter / 2 );


				particle_x_velocity[jstep] = ( 2 * ((float)rand()/RAND_MAX)) - 1;

				particle_y_velocity[jstep] = ( 2 * ((float)rand()/RAND_MAX)) - 1;

				particle_z_velocity[jstep] = ( 2 * ((float)rand()/RAND_MAX)) - 1;


				absolute_velocity = sqrt( ( particle_x_velocity[jstep] * particle_x_velocity[jstep] ) + ( particle_y_velocity[jstep] * particle_y_velocity[jstep] ) + ( particle_z_velocity[jstep] * particle_z_velocity [jstep]) );

				particle_x_velocity[jstep] = particle_x_velocity[jstep] / absolute_velocity;
	
				particle_y_velocity[jstep] = particle_y_velocity[jstep] / absolute_velocity;

				particle_z_velocity[jstep] = particle_z_velocity[jstep] / absolute_velocity;


				particle_x_velocity[jstep] = particle_x_velocity[jstep] * mean_molecular_speed;

				particle_y_velocity[jstep] = particle_y_velocity[jstep] * mean_molecular_speed;

				particle_z_velocity[jstep] = particle_z_velocity[jstep] * mean_molecular_speed;
				

			}


			// Inter-particle Collisions


			for ( int lstep = 0; lstep < local_pressure_resolution; lstep ++ )
			{

				if ( ( particle_z_position[jstep] >= ( local_space_size * lstep ) ) && ( particle_z_position[jstep] < ( local_space_size * ( lstep + 1 ) ) ) )
				{

					if ( (float)rand()/RAND_MAX < collision_probability[lstep] )
					{

						particle_x_velocity[jstep] = ( 2 * ((float)rand()/RAND_MAX)) - 1;

						particle_y_velocity[jstep] = ( 2 * ((float)rand()/RAND_MAX)) - 1;

						particle_z_velocity[jstep] = ( 2 * ((float)rand()/RAND_MAX)) - 1;


						absolute_velocity = sqrt( ( particle_x_velocity[jstep] * particle_x_velocity[jstep] ) + ( particle_y_velocity[jstep] * particle_y_velocity[jstep] ) + ( particle_z_velocity[jstep] * particle_z_velocity [jstep]) );

						particle_x_velocity[jstep] = particle_x_velocity[jstep] / absolute_velocity;
	
						particle_y_velocity[jstep] = particle_y_velocity[jstep] / absolute_velocity;

						particle_z_velocity[jstep] = particle_z_velocity[jstep] / absolute_velocity;


						particle_x_velocity[jstep] = particle_x_velocity[jstep] * mean_molecular_speed;

						particle_y_velocity[jstep] = particle_y_velocity[jstep] * mean_molecular_speed;

						particle_z_velocity[jstep] = particle_z_velocity[jstep] * mean_molecular_speed;

					}

					break;

				}
				
			}


			// "Collision" with front of tube


			if ( particle_z_position[jstep] <= 0 )
			{

				particle_z_position[jstep] = 0;


				particle_x_velocity[jstep] = ( 2 * ((float)rand()/RAND_MAX)) - 1;

				particle_y_velocity[jstep] = ( 2 * ((float)rand()/RAND_MAX)) - 1;

				particle_z_velocity[jstep] = (float)rand()/RAND_MAX;


				absolute_velocity = sqrt( ( particle_x_velocity[jstep] * particle_x_velocity[jstep] ) + ( particle_y_velocity[jstep] * particle_y_velocity[jstep] ) + ( particle_z_velocity[jstep] * particle_z_velocity [jstep]) );

				particle_x_velocity[jstep] = particle_x_velocity[jstep] / absolute_velocity;
	
				particle_y_velocity[jstep] = particle_y_velocity[jstep] / absolute_velocity;

				particle_z_velocity[jstep] = particle_z_velocity[jstep] / absolute_velocity;


				particle_x_velocity[jstep] = particle_x_velocity[jstep] * mean_molecular_speed;

				particle_y_velocity[jstep] = particle_y_velocity[jstep] * mean_molecular_speed;

				particle_z_velocity[jstep] = particle_z_velocity[jstep] * mean_molecular_speed;

			}


			// Handle Exiting Particles


			if ( particle_z_position[jstep] >= tube_length )
			{

				for ( lstep = jstep; lstep < (particles - 1); lstep ++ )
				{

					particle_x_position[lstep] = particle_x_position[lstep+1];

					particle_y_position[lstep] = particle_y_position[lstep+1];

					particle_z_position[lstep] = particle_z_position[lstep+1];


					particle_x_velocity[lstep] = particle_x_velocity[lstep+1];

					particle_y_velocity[lstep] = particle_y_velocity[lstep+1];

					particle_z_velocity[lstep] = particle_z_velocity[lstep+1];


					

				}

				particles --;

				jstep --;

				long_outlet_particle_rate ++;

			}

		}


		// Particles enter tube based on inlet pressure and 1st local pressure


		pressure_ratio = local_pressure[0] / inlet_pressure;


		if ( pressure_ratio >= 1 )
		{
		
			inlet_conductance = 0;

		}
		else if ( ( pressure_ratio < 1 ) & ( pressure_ratio >= 0.52 ) )
		{

			inlet_conductance = ( (766 * local_area) / ( 1 - pressure_ratio ) ) * pow( pressure_ratio, 0.712 ) * sqrt ( 1 - ( pow( pressure_ratio, 0.288 ) ) );

		}
		else if ( pressure_ratio < 0.52 )
		{
			inlet_conductance = ( 200 * local_area ) / ( 1 - pressure_ratio );

		}


		inlet_flow_rate = inlet_conductance * ( inlet_pressure - local_pressure[0] );

		inlet_volume_rate = inlet_flow_rate / inlet_pressure;

		inlet_molecule_rate = ( inlet_pressure * inlet_volume_rate ) / ( k * inlet_temperature );

		inlet_particle_rate = ( inlet_molecule_rate / particle_weight ) * time_step;

		inlet_particle_rate = floor( inlet_particle_rate + 0.5 );

		

		// Give entering particles coordinates & velocities
	

		for ( jstep = 0; jstep < inlet_particle_rate; jstep ++ )
		{

			particle_r_position = ( ( tube_diameter / 2 ) * ((float)rand()/RAND_MAX));

			particle_theta_position = ( ( 2 * 3.141592 ) * ((float)rand()/RAND_MAX));

			particle_z_position[particles + jstep] = 0;


			particle_x_position[particles + jstep] = particle_r_position * cos(particle_theta_position);

			particle_y_position[particles + jstep] = particle_r_position * sin(particle_theta_position);


			particle_x_velocity[particles + jstep] = ( 2 * ((float)rand()/RAND_MAX)) - 1;

			particle_y_velocity[particles + jstep] = ( 2 * ((float)rand()/RAND_MAX)) - 1;

			particle_z_velocity[particles + jstep] = (float)rand()/RAND_MAX;


			absolute_velocity = sqrt( ( particle_x_velocity[particles + jstep] * particle_x_velocity[particles + jstep] ) + ( particle_y_velocity[particles + jstep] * particle_y_velocity[particles + jstep] ) + ( particle_z_velocity[particles + jstep] * particle_z_velocity [particles + jstep]) );

			particle_x_velocity[particles + jstep] = particle_x_velocity[particles + jstep] / absolute_velocity;

			particle_y_velocity[particles + jstep] = particle_y_velocity[particles + jstep] / absolute_velocity;

			particle_z_velocity[particles + jstep] = particle_z_velocity[particles + jstep] / absolute_velocity;


			particle_x_velocity[particles + jstep] = particle_x_velocity[particles + jstep] * mean_molecular_speed;

			particle_y_velocity[particles + jstep] = particle_y_velocity[particles + jstep] * mean_molecular_speed;

			particle_z_velocity[particles + jstep] = particle_z_velocity[particles + jstep] * mean_molecular_speed;


		}


		particles = particles + inlet_particle_rate;

		long_inlet_particle_rate = long_inlet_particle_rate + inlet_particle_rate;



		// If Equilibrium or max_particles is reached, break
		
		

		if ( particles >= ( 0.99 * max_particles ) )
		{

			break;

		}


		

		if ( ( ( istep % sampling_resolution ) == 0 ) & ( istep > 0 ) )
		{

			long_inlet_molecule_rate = ( long_inlet_particle_rate * particle_weight ) / ( sampling_resolution * time_step );
			
			long_inlet_volume_rate = ( long_inlet_molecule_rate * k * inlet_temperature ) / inlet_pressure;

			long_inlet_volume_rate = ( long_inlet_volume_rate * 1000 ) * 60;  // convert ( m^3 / s ) to ( Liters / min )
			

			long_outlet_molecule_rate = ( long_outlet_particle_rate * particle_weight ) / ( sampling_resolution * time_step );

			long_outlet_volume_rate = ( long_outlet_molecule_rate * k * inlet_temperature ) / inlet_pressure;  // local_pressure[local_pressure_resolution - 1];

			long_outlet_volume_rate = ( long_outlet_volume_rate * 1000 ) * 60;  // convert ( m^3 / s ) to ( Liters / min )


			run_time = run_time + ( sampling_resolution * time_step );

			average_first_pressure = average_first_pressure / sampling_resolution;


			cout << "Physical Run Time: " << run_time << endl;

			cout << "Sampling Interval: " << sampling_resolution * time_step << " seconds " << endl; 

			cout << "Incoming: " << long_inlet_particle_rate << endl << "Outgoing: " << long_outlet_particle_rate << endl << "Total: " << particles << endl;
			
			cout << "Flow Rate In: " << long_inlet_volume_rate << " L/min" << endl << "Flow Rate Out: " << long_outlet_volume_rate << " L/min" << endl;
			
			cout << "Inlet Pressure: " << inlet_pressure << " Pascals" << endl << "1st Average Pressure: " << average_first_pressure << " Pascals" << endl << endl;

 
			cout << "Particle Distribution: " << endl;


			for ( kstep = 0; kstep < local_pressure_resolution; kstep ++ )
			{

				cout << local_particle_count[kstep] << endl;

			}


			cout << endl;


			long_inlet_particle_rate = 0;

			long_outlet_particle_rate = 0;

			average_first_pressure = 0;

		}

	}


	return;

}
