//////////////////////////
// output.hpp
//////////////////////////
// 
// Output of snapshots and spectra
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: June 2018
//
//////////////////////////

#ifndef OUTPUT_HEADER
#define OUTPUT_HEADER

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

using namespace std;


//////////////////////////
// writeSnapshots
//////////////////////////
// Description:
//   output of snapshots
// 
// Arguments:
//   sim            simulation metadata structure
//   cosmo          cosmological parameter structure
//   fourpiG        4 pi G (in code units)
//   a              scale factor
//   snapcount      snapshot index
//   h5filename     base name for HDF5 output file
//   pcls_cdm       pointer to (uninitialized) particle handler for CDM
//   pcls_b         pointer to (uninitialized) particle handler for baryons
//   pcls_ncdm      array of (uninitialized) particle handlers for
//                  non-cold DM (may be set to NULL)
//   phi            pointer to allocated field
//   chi            pointer to allocated field
//   Bi             pointer to allocated field
//   source         pointer to allocated field
//   Sij            pointer to allocated field
//   scalarFT       pointer to allocated field
//   BiFT           pointer to allocated field
//   SijFT          pointer to allocated field
//   plan_phi       pointer to FFT planner
//   plan_chi       pointer to FFT planner
//   plan_Bi        pointer to FFT planner
//   plan_source    pointer to FFT planner
//   plan_Sij       pointer to FFT planner
//   Bi_check       pointer to allocated field (or NULL)
//   BiFT_check     pointer to allocated field (or NULL)
//   plan_Bi_check  pointer to FFT planner (or NULL)
//
// Returns:
// 
//////////////////////////

void writeSnapshots(metadata & sim, cosmology & cosmo, const double fourpiG, const double a, const int done_hij, const int snapcount, string h5filename, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, Field<Real> * phi, Field<Real> * chi, Field<Real> * Bi, Field<Real> * source, Field<Real> * Sij, Field<Cplx> * scalarFT, Field<Cplx> * BiFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_phi, PlanFFT<Cplx> * plan_chi, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_source, PlanFFT<Cplx> * plan_Sij, Field<Real> * Bi_check = NULL, Field<Cplx> * BiFT_check = NULL, PlanFFT<Cplx> * plan_Bi_check = NULL)
{
	char filename[2*PARAM_MAX_LENGTH+24];
	char buffer[64];
	int i;
	gadget2_header hdr;
	Site x(phi->lattice());
	Real divB, curlB, divh, traceh, normh;

	sprintf(filename, "%03d", snapcount);
			
#ifdef EXTERNAL_IO
	while (ioserver.openOstream()== OSTREAM_FAIL);
	
	if (sim.out_snapshot & MASK_PCLS)
	{
		pcls_cdm->saveHDF5_server_open(h5filename + filename + "_cdm");
		if (sim.baryon_flag)
			pcls_b->saveHDF5_server_open(h5filename + filename + "_b");
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			sprintf(buffer, "_ncdm%d", i);
			pcls_ncdm[i].saveHDF5_server_open(h5filename + filename + buffer);
		}
	}
	
	if (sim.out_snapshot & MASK_T00)
		source->saveHDF5_server_open(h5filename + filename + "_T00");
				
	if (sim.out_snapshot & MASK_B)
		Bi->saveHDF5_server_open(h5filename + filename + "_B");
	
	if (sim.out_snapshot & MASK_PHI)
		phi->saveHDF5_server_open(h5filename + filename + "_phi");
				
	if (sim.out_snapshot & MASK_CHI)
		chi->saveHDF5_server_open(h5filename + filename + "_chi");
	
	if (sim.out_snapshot & MASK_HIJ)
		Sij->saveHDF5_server_open(h5filename + filename + "_hij");
				
#ifdef CHECK_B
	if (sim.out_snapshot & MASK_B)
		Bi_check->saveHDF5_server_open(h5filename + filename + "_B_check");
#endif
#endif		
			
	if (sim.out_snapshot & MASK_RBARE || sim.out_snapshot & MASK_POT)
	{
		projection_init(source);
		scalarProjectionCIC_project(pcls_cdm, source);
		if (sim.baryon_flag)
			scalarProjectionCIC_project(pcls_b, source);
		for (i = 0; i < cosmo.num_ncdm; i++)
			scalarProjectionCIC_project(pcls_ncdm+i, source);
		scalarProjectionCIC_comm(source);
	}

	if (sim.out_snapshot & MASK_RBARE)
	{
		if (sim.downgrade_factor > 1)
			source->saveHDF5_coarseGrain3D(h5filename + filename + "_rhoN.h5", sim.downgrade_factor);
		else
			source->saveHDF5(h5filename + filename + "_rhoN.h5");
	}
			
	if (sim.out_snapshot & MASK_POT)
	{
		plan_source->execute(FFT_FORWARD);				
		solveModifiedPoissonFT(*scalarFT, *scalarFT, fourpiG / a);
		plan_source->execute(FFT_BACKWARD);
		if (sim.downgrade_factor > 1)
			source->saveHDF5_coarseGrain3D(h5filename + filename + "_psiN.h5", sim.downgrade_factor);
		else
			source->saveHDF5(h5filename + filename + "_psiN.h5");
	}
				
	if (sim.out_snapshot & MASK_T00)
	{
		projection_init(source);
		if (sim.gr_flag > 0)
		{
			projection_T00_project(pcls_cdm, source, a, phi);
			if (sim.baryon_flag)
				projection_T00_project(pcls_b, source, a, phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
				projection_T00_project(pcls_ncdm+i, source, a, phi);
		}
		else
		{
			scalarProjectionCIC_project(pcls_cdm, source);
			if (sim.baryon_flag)
				scalarProjectionCIC_project(pcls_b, source);
			for (i = 0; i < cosmo.num_ncdm; i++)
				scalarProjectionCIC_project(pcls_ncdm+i, source);
		}
		projection_T00_comm(source);
#ifdef EXTERNAL_IO
		source->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (sim.downgrade_factor > 1)
			source->saveHDF5_coarseGrain3D(h5filename + filename + "_T00.h5", sim.downgrade_factor);
		else
			source->saveHDF5(h5filename + filename + "_T00.h5");
#endif
	}
				
	if (sim.out_snapshot & MASK_B)
	{
		if (sim.gr_flag == 0)
		{
			plan_Bi->execute(FFT_BACKWARD);
		}
		for (x.first(); x.test(); x.next())
		{
			(*Bi)(x,0) /= a * a * sim.numpts;
			(*Bi)(x,1) /= a * a * sim.numpts;
			(*Bi)(x,2) /= a * a * sim.numpts;
		}
		Bi->updateHalo();
				
		computeVectorDiagnostics(*Bi, divB, curlB);			
		COUT << " B diagnostics: max |divB| = " << divB << ", max |curlB| = " << curlB << endl;

#ifdef EXTERNAL_IO
		Bi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (sim.downgrade_factor > 1)
			Bi->saveHDF5_coarseGrain3D(h5filename + filename + "_B.h5", sim.downgrade_factor);
		else				
			Bi->saveHDF5(h5filename + filename + "_B.h5");
#endif
				
		if (sim.gr_flag > 0)
		{
			plan_Bi->execute(FFT_BACKWARD);
			Bi->updateHalo();
		}
	}
			
	if (sim.out_snapshot & MASK_PHI)
#ifdef EXTERNAL_IO
		phi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (sim.downgrade_factor > 1)
			phi->saveHDF5_coarseGrain3D(h5filename + filename + "_phi.h5", sim.downgrade_factor);
		else
			phi->saveHDF5(h5filename + filename + "_phi.h5");
#endif
				
	if (sim.out_snapshot & MASK_CHI)
#ifdef EXTERNAL_IO
		chi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else	
		if (sim.downgrade_factor > 1)
			chi->saveHDF5_coarseGrain3D(h5filename + filename + "_chi.h5", sim.downgrade_factor);
		else
			chi->saveHDF5(h5filename + filename + "_chi.h5");
#endif
				
	if (sim.out_snapshot & MASK_HIJ)
	{
		if (done_hij == 0)
		{
			projectFTtensor(*SijFT, *SijFT);
			plan_Sij->execute(FFT_BACKWARD);
			Sij->updateHalo();
		}
				
		computeTensorDiagnostics(*Sij, divh, traceh, normh);
		COUT << " GW diagnostics: max |divh| = " << divh << ", max |traceh| = " << traceh << ", max |h| = " << normh << endl;

#ifdef EXTERNAL_IO
		Sij->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else	
		if (sim.downgrade_factor > 1)
			Sij->saveHDF5_coarseGrain3D(h5filename + filename + "_hij.h5", sim.downgrade_factor);
		else
			Sij->saveHDF5(h5filename + filename + "_hij.h5");
#endif
	}

	if (sim.out_snapshot & MASK_TIJ)
	{						
		projection_init(Sij);
		projection_Tij_project(pcls_cdm, Sij, a, phi);
		if (sim.baryon_flag)
			projection_Tij_project(pcls_b, Sij, a, phi);
		for (i = 0; i < cosmo.num_ncdm; i++)
			projection_Tij_project(pcls_ncdm+i, Sij, a, phi);
		projection_Tij_comm(Sij);

		if (sim.downgrade_factor > 1)
			Sij->saveHDF5_coarseGrain3D(h5filename + filename + "_Tij.h5", sim.downgrade_factor);
		else
			Sij->saveHDF5(h5filename + filename + "_Tij.h5");
	}
			
	if (sim.out_snapshot & MASK_P)
	{
		projection_init(Bi);
		projection_T0i_project(pcls_cdm, Bi, phi);
		if (sim.baryon_flag)
			projection_T0i_project(pcls_b, Bi, phi);
		for (i = 0; i < cosmo.num_ncdm; i++)
			projection_T0i_project(pcls_ncdm+i, Bi, phi);
		projection_T0i_comm(Bi);
		if (sim.downgrade_factor > 1)
			Bi->saveHDF5_coarseGrain3D(h5filename + filename + "_p.h5", sim.downgrade_factor);
		else
			Bi->saveHDF5(h5filename + filename + "_p.h5");
		if (sim.gr_flag > 0)
		{
			plan_Bi->execute(FFT_BACKWARD);
			Bi->updateHalo();
		}
	}
				
#ifdef CHECK_B
	if (sim.out_snapshot & MASK_B)
	{
		if (sim.vector_flag == VECTOR_PARABOLIC)
		{
			projection_init(Bi_check);
			projection_T0i_project(pcls_cdm, Bi_check, phi);
			if (sim.baryon_flag)
				projection_T0i_project(pcls_b, Bi_check, phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
				projection_T0i_project(pcls_ncdm+i, Bi_check, phi);
			projection_T0i_comm(Bi_check);
			plan_Bi_check->execute(FFT_FORWARD);
			projectFTvector(*BiFT_check, *BiFT_check, fourpiG / (double) sim.numpts / (double) sim.numpts);
		}
		plan_Bi_check->execute(FFT_BACKWARD);
			
		for (x.first(); x.test(); x.next())
		{
			(*Bi_check)(x,0) /= a * a * sim.numpts;
			(*Bi_check)(x,1) /= a * a * sim.numpts;
			(*Bi_check)(x,2) /= a * a * sim.numpts;
		}
#ifdef EXTERNAL_IO
		Bi_check->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (sim.downgrade_factor > 1)
			Bi_check->saveHDF5_coarseGrain3D(h5filename + filename + "_B_check.h5", sim.downgrade_factor);
		else
			Bi_check->saveHDF5(h5filename + filename + "_B_check.h5");
#endif
	}
#endif

	if (sim.out_snapshot & MASK_GADGET)
	{
		if (sim.out_snapshot & MASK_MULTI)
			hdr.num_files = parallel.grid_size()[1];
		else
			hdr.num_files = 1;
		hdr.Omega0 = cosmo.Omega_m;
		hdr.OmegaLambda = cosmo.Omega_Lambda;
		hdr.HubbleParam = cosmo.h;
		hdr.BoxSize = sim.boxsize / GADGET_LENGTH_CONVERSION;
		hdr.flag_sfr = 0;
		hdr.flag_cooling = 0;
		hdr.flag_feedback = 0;
		hdr.flag_age = 0;
		hdr.flag_metals = 0;
		for (i = 0; i < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8 - 2 * 4 - 6 * 4; i++)
			hdr.fill[i] = 0;
		for (i = 0; i < 6; i++)
		{
			hdr.npart[i] = 0;
			hdr.npartTotal[i] = 0;
			hdr.npartTotalHW[i] = 0;
			hdr.mass[i] = 0.;
		}

		hdr.time = a;
		hdr.redshift = (1./a) - 1.;
				
		hdr.npart[1] = (uint32_t) ((sim.numpcl[0] / sim.tracer_factor[0]) % (1l << 32));
		hdr.npartTotal[1] = hdr.npart[1];
		hdr.npartTotalHW[1] = (uint32_t) ((sim.numpcl[0] / sim.tracer_factor[0]) >> 32);
		if (sim.baryon_flag)
			hdr.mass[1] = (double) sim.tracer_factor[0] * C_RHO_CRIT * cosmo.Omega_cdm * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;
		else
			hdr.mass[1] = (double) sim.tracer_factor[0] * C_RHO_CRIT * (cosmo.Omega_cdm + cosmo.Omega_b) * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;

		if (hdr.npartTotalHW[1] != 0 && hdr.num_files == 1)
		{
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": number of particles (" << (sim.numpcl[0] / sim.tracer_factor[0]) << ") in Gadget2 file exceeds limit (4294967295). Try using multi-Gadget2 output format." << endl;
		}
		else
			pcls_cdm->saveGadget2(h5filename + filename + "_cdm", hdr, sim.tracer_factor[0]);
				
		if (sim.baryon_flag)
		{
			hdr.npart[1] = (uint32_t) ((sim.numpcl[1] / sim.tracer_factor[1]) % (1l << 32));
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.npartTotalHW[1] = (uint32_t) ((sim.numpcl[1] / sim.tracer_factor[1]) >> 32);
			hdr.mass[1] = (double) sim.tracer_factor[1] * C_RHO_CRIT * cosmo.Omega_b * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[1] / GADGET_MASS_CONVERSION;
			if (hdr.npartTotalHW[1] != 0 && hdr.num_files == 1)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": number of particles (" << (sim.numpcl[1] / sim.tracer_factor[1]) << ") in Gadget2 file exceeds limit (4294967295). Try using multi-Gadget2 output format." << endl;
			}
			else
				pcls_b->saveGadget2(h5filename + filename + "_b", hdr, sim.tracer_factor[1]);
		}
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			sprintf(buffer, "_ncdm%d", i);
			hdr.npart[1] = (uint32_t) ((sim.numpcl[i+1+sim.baryon_flag] / sim.tracer_factor[i+1+sim.baryon_flag]) % (1l << 32));
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.npartTotalHW[1] = (uint32_t) ((sim.numpcl[i+1+sim.baryon_flag] / sim.tracer_factor[i+1+sim.baryon_flag]) >> 32);
			hdr.mass[1] = (double) sim.tracer_factor[i+1+sim.baryon_flag] * C_RHO_CRIT * cosmo.Omega_ncdm[i] * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[i+1+sim.baryon_flag] / GADGET_MASS_CONVERSION;
			if (hdr.npartTotalHW[1] != 0 && hdr.num_files == 1)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": number of particles (" << (sim.numpcl[i+1+sim.baryon_flag] / sim.tracer_factor[i+1+sim.baryon_flag]) << ") in Gadget2 file exceeds limit (4294967295). Try using multi-Gadget2 output format." << endl;
			}
			else
				pcls_ncdm[i].saveGadget2(h5filename + filename + buffer, hdr, sim.tracer_factor[i+1+sim.baryon_flag]);
		}
	}
			
	if (sim.out_snapshot & MASK_PCLS)
	{
#ifdef EXTERNAL_IO
		pcls_cdm->saveHDF5_server_write();
		if (sim.baryon_flag)
			pcls_b->saveHDF5_server_write();
		for (i = 0; i < cosmo.num_ncdm; i++)
			pcls_ncdm[i].saveHDF5_server_write();
#else
		pcls_cdm->saveHDF5(h5filename + filename + "_cdm", 1);
		if (sim.baryon_flag)
			pcls_b->saveHDF5(h5filename + filename + "_b", 1);
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			sprintf(buffer, "_ncdm%d", i);
			pcls_ncdm[i].saveHDF5(h5filename + filename + buffer, 1);
		}
#endif
	}
			
#ifdef EXTERNAL_IO
	ioserver.closeOstream();
#endif
}


#if LIGHTCONE_THICKNESS > 3
#error This implementation of writeLightcones does not support LIGHTCONE_THICKNESS > 3
#endif

void writeLightcones(metadata & sim, cosmology & cosmo, const double fourpiG, const double a, const double tau, const double dtau, const double dtau_old, const double dtau_older, const double maxvel, const int cycle, string h5filename, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, Field<Real> * phi, Field<Real> * chi, Field<Real> * Bi, Field<Real> * Sij, Field<Cplx> * BiFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_Sij, Field<Real> * lcbuffer[LIGHTCONE_MAX_FIELDS*LIGHTCONE_THICKNESS], Lattice * lclat[LIGHTCONE_MAX_FIELDS], int & done_hij)
{
	int i, j, n, p;
	double d;
	double vertex[MAX_INTERSECTS][3];
	double domain[6];
	double pos[3];
	double s[4];
	char filename[2*PARAM_MAX_LENGTH+24];
	char buffer[64];
	FILE * outfile;
	gadget2_header hdr;
	Site xsim;
	Site xlc;
	int done_B = 0;
	
	done_hij = 0;
	
	domain[0] = -0.5;
	domain[1] = phi->lattice().coordSkip()[1] - 0.5;
	domain[2] = phi->lattice().coordSkip()[0] - 0.5;
	for (i = 0; i < 3; i++)
		domain[i+3] = domain[i] + phi->lattice().sizeLocal(i) + 1.;

	for (i = 0; i < 6; i++)
		domain[i] /= (double) sim.numpts;

	for (i = 0; i < sim.num_lightcone; i++)
	{
		if (parallel.isRoot())
		{
			if (sim.num_lightcone > 1)
				sprintf(filename, "%s%s%d_info.dat", sim.output_path, sim.basename_lightcone, i);
			else
				sprintf(filename, "%s%s_info.dat", sim.output_path, sim.basename_lightcone);

			outfile = fopen(filename, "a");
			if (outfile == NULL)
			{
				cout << " error opening file for lightcone info!" << endl;
			}
			else if (cycle == 0)
			{
				if (sim.num_lightcone > 1)
					fprintf(outfile, "# information file for lightcone %d\n# geometric parameters:\n# vertex = (%f, %f, %f) Mpc/h\n# redshift = %f\n# distance = (%f - %f) Mpc/h\n# opening half-angle = %f degrees\n# direction = (%f, %f, %f)\n# cycle   tau/boxsize    a              pcl_inner        pcl_outer        metric_inner", i, sim.lightcone[i].vertex[0]*sim.boxsize, sim.lightcone[i].vertex[1]*sim.boxsize, sim.lightcone[i].vertex[2]*sim.boxsize, sim.lightcone[i].z, sim.lightcone[i].distance[0]*sim.boxsize, sim.lightcone[i].distance[1]*sim.boxsize, sim.lightcone[i].opening * 180. / M_PI, sim.lightcone[i].direction[0], sim.lightcone[i].direction[1], sim.lightcone[i].direction[2]);
				else
					fprintf(outfile, "# information file for lightcone\n# geometric parameters:\n# vertex = (%f, %f, %f) Mpc/h\n# redshift = %f\n# distance = (%f - %f) Mpc/h\n# opening half-angle = %f degrees\n# direction = (%f, %f, %f)\n# cycle   tau/boxsize    a              pcl_inner        pcl_outer        metric_inner", sim.lightcone[i].vertex[0]*sim.boxsize, sim.lightcone[i].vertex[1]*sim.boxsize, sim.lightcone[i].vertex[2]*sim.boxsize, sim.lightcone[i].z, sim.lightcone[i].distance[0]*sim.boxsize, sim.lightcone[i].distance[1]*sim.boxsize, sim.lightcone[i].opening * 180. / M_PI, sim.lightcone[i].direction[0], sim.lightcone[i].direction[1], sim.lightcone[i].direction[2]);

				for (j = 1; j < LIGHTCONE_THICKNESS; j++)
					fprintf(outfile, "     metric_(%d/%d)", j-1, j);

				fprintf(outfile, "     metric_outer\n");
			}
		}

		d = particleHorizon(1. / (1. + sim.lightcone[i].z), fourpiG, cosmo);

#if LIGHTCONE_THICKNESS == 2
		s[0] = d - tau - dtau;
		s[1] = d - tau;
		s[2] = d - tau + dtau_old;
#elif LIGHTCONE_THICKNESS == 3
		s[0] = d - tau - 2. * dtau;
		s[1] = d - tau - 0.5 * dtau;
		s[2] = d - tau + 0.5 * dtau_old;
		s[3] = d - tau + dtau_old + 0.5 * dtau_older;
#else
		s[0] = d - tau - 0.5 * dtau;
		s[1] = d - tau + 0.5 * dtau_old;
#endif

		if (sim.lightcone[i].distance[0] > s[1] && sim.lightcone[i].distance[1] <= s[LIGHTCONE_THICKNESS] && s[LIGHTCONE_THICKNESS] > 0.)
		{
			if (parallel.isRoot() && outfile != NULL)
			{
				fprintf(outfile, "%6d   %e   %e   %2.12f   %2.12f   %2.12f", cycle, tau, a, d - tau - 0.5 * dtau, d - tau + 0.5 * dtau_old, s[0]);
				for (j = 1; j < LIGHTCONE_THICKNESS; j++)
					fprintf(outfile, "   %2.12f", s[j]);
				fprintf(outfile, "   %2.12f\n", s[LIGHTCONE_THICKNESS]);
				fclose(outfile);

				if (sim.num_lightcone > 1)
					sprintf(filename, "%s%s%d_info.bin", sim.output_path, sim.basename_lightcone, i);
				else
					sprintf(filename, "%s%s_info.bin", sim.output_path, sim.basename_lightcone);

				outfile = fopen(filename, "a");
				if (outfile == NULL)
				{
					cout << " error opening file for lightcone info!" << endl;
				}
				else
				{
					((double *) buffer)[0] = tau;
					((double *) buffer)[1] = a;
					((double *) buffer)[2] = d - tau - 0.5 * dtau;
					((double *) buffer)[3] = d - tau + 0.5 * dtau_old;

					fwrite((const void *) &cycle, sizeof(int), 1, outfile);
					fwrite((const void *) buffer, sizeof(double), 4, outfile);
					fwrite((const void *) s, sizeof(double), LIGHTCONE_THICKNESS+1, outfile);

					fclose(outfile);
				}
			}

			for (j = 0; j < LIGHTCONE_THICKNESS; j++)
			{
				if (sim.lightcone[i].distance[0] <= s[j+1] || sim.lightcone[i].distance[1] > s[j+1] || s[j+1] <= 0) continue;

				n = findIntersectingLightcones(sim.lightcone[i], s[j+1], s[j], domain, vertex);

				if (n > 0 && sim.out_lightcone[i] & MASK_PHI)
				{
					xsim.initialize(phi->lattice());
					xlc.initialize((*lclat[LIGHTCONE_PHI_OFFSET]));

					for (xsim.first(), xlc.first(); xsim.test(); xsim.next(), xlc.next())
					{
						for (p = 0; p < 3; p++)
							pos[p] = (double) xsim.coord(p) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_PHI_OFFSET+j])(xlc) = (*phi)(xsim);
								break;
							}
						}
					}
				}

				if (n > 0 && sim.out_lightcone[i] & MASK_CHI)
				{
					xsim.initialize(chi->lattice());
					xlc.initialize((*lclat[LIGHTCONE_CHI_OFFSET]));

					for (xsim.first(), xlc.first(); xsim.test(); xsim.next(), xlc.next())
					{
						for (p = 0; p < 3; p++)
							pos[p] = (double) xsim.coord(p) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_CHI_OFFSET+j])(xlc) = (*chi)(xsim);
								break;
							}
						}
					}
				}

				if (sim.gr_flag == 0 && sim.out_lightcone[i] & MASK_B && done_B == 0)
				{
					plan_Bi->execute(FFT_BACKWARD);
					Bi->updateHalo();
					done_B = 1;
				}

				if (n > 0 && sim.out_lightcone[i] & MASK_B)
				{
					xsim.initialize(Bi->lattice());
					xlc.initialize((*lclat[LIGHTCONE_B_OFFSET]));

					for (xsim.first(), xlc.first(); xsim.test(); xsim.next(), xlc.next())
					{
#ifdef LIGHTCONE_INTERPOLATE
						for (p = 0; p < 3; p++)
							pos[p] = (double) xsim.coord(p) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_B_OFFSET+3*j])(xlc) = ((*Bi)(xsim, 0) + (*Bi)(xsim-0, 0)) / (2. * a * a * sim.numpts);
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_B_OFFSET+3*j+1])(xlc) = ((*Bi)(xsim, 1) + (*Bi)(xsim-1, 1)) / (2. * a * a * sim.numpts);
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_B_OFFSET+3*j+2])(xlc) = ((*Bi)(xsim, 2) + (*Bi)(xsim-2, 2)) / (2. * a * a * sim.numpts);
								break;
							}
						}
#else
						pos[0] = (0.5 + (double) xsim.coord(0)) / (double) sim.numpts;
						pos[1] = (double) xsim.coord(1) / (double) sim.numpts;
						pos[2] = (double) xsim.coord(2) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_B_OFFSET+3*j])(xlc) = (*Bi)(xsim, 0) / (a * a * sim.numpts);
								break;
							}
						}

						pos[0] = (double) xsim.coord(0) / (double) sim.numpts;
						pos[1] = (0.5 + (double) xsim.coord(1)) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_B_OFFSET+3*j+1])(xlc) = (*Bi)(xsim, 1) / (a * a * sim.numpts);
								break;
							}
						}

						pos[1] = (double) xsim.coord(1) / (double) sim.numpts;
						pos[2] = (0.5 + (double) xsim.coord(2)) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_B_OFFSET+3*j+2])(xlc) = (*Bi)(xsim, 2) / (a * a * sim.numpts);
								break;
							}
						}
#endif
					}
				}

				if (sim.out_lightcone[i] & MASK_HIJ && done_hij == 0)
				{
					projectFTtensor(*SijFT, *SijFT);
					plan_Sij->execute(FFT_BACKWARD);
					Sij->updateHalo();
					done_hij = 1;
				}

				if (n > 0 && sim.out_lightcone[i] & MASK_HIJ)
				{
					xsim.initialize(Sij->lattice());
					xlc.initialize((*lclat[LIGHTCONE_HIJ_OFFSET]));

					for (xsim.first(), xlc.first(); xsim.test(); xsim.next())
					{
#ifdef LIGHTCONE_DOWNGRADE
						if ((xsim.coord(0) % LIGHTCONE_DOWNGRADE) > 0 || (xsim.coord(1) % LIGHTCONE_DOWNGRADE) > 0 || (xsim.coord(2) % LIGHTCONE_DOWNGRADE) > 0)
							continue;
#endif

						for (p = 0; p < 3; p++)
							pos[p] = (double) xsim.coord(p) / (double) sim.numpts;

#ifdef LIGHTCONE_INTERPOLATE
						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j])(xlc) = (*Sij)(xsim, 0, 0);
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+1])(xlc) = ((*Sij)(xsim, 0, 1) + (*Sij)(xsim-0, 0, 1) + (*Sij)(xsim-1, 0, 1) + (*Sij)(xsim-0-1, 0, 1)) / 4.;
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+2])(xlc) = ((*Sij)(xsim, 0, 2) + (*Sij)(xsim-0, 0, 2) + (*Sij)(xsim-2, 0, 2) + (*Sij)(xsim-0-2, 0, 2)) / 4.;
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+3])(xlc) = (*Sij)(xsim, 1, 1);
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+4])(xlc) = ((*Sij)(xsim, 1, 2) + (*Sij)(xsim-1, 1, 2) + (*Sij)(xsim-2, 1, 2) + (*Sij)(xsim-1-2, 1, 2)) / 4.;
								break;
							}
						}
#else
#ifdef LIGHTCONE_DOWNGRADE
						for (p = 0; p < 3; p++)
							pos[p] += 0.5 / (double) sim.numpts;
							
						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j])(xlc) = ((*Sij)(xsim, 0, 0) + (*Sij)(xsim+0, 0, 0) + (*Sij)(xsim+1, 0, 0) + (*Sij)(xsim+1+0, 0, 0) + (*Sij)(xsim+2, 0, 0) + (*Sij)(xsim+2+0, 0, 0) + (*Sij)(xsim+2+1, 0, 0) + (*Sij)(xsim+2+1+0, 0, 0)) / 8.;
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+1])(xlc) = ((*Sij)(xsim, 0, 1) + (*Sij)(xsim+2, 0, 1)) / 2.;
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+2])(xlc) = ((*Sij)(xsim, 0, 2) + (*Sij)(xsim+1, 0, 2)) / 2.;
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+3])(xlc) = ((*Sij)(xsim, 1, 1) + (*Sij)(xsim+0, 1, 1) + (*Sij)(xsim+1, 1, 1) + (*Sij)(xsim+1+0, 1, 1) + (*Sij)(xsim+2, 1, 1) + (*Sij)(xsim+2+0, 1, 1) + (*Sij)(xsim+2+1, 1, 1) + (*Sij)(xsim+2+1+0, 1, 1)) / 8.;
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+4])(xlc) = ((*Sij)(xsim, 1, 2) + (*Sij)(xsim+0, 1, 2)) / 2.;
								break;
							}
						}
#else
						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j])(xlc) = (*Sij)(xsim, 0, 0);
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+3])(xlc) = (*Sij)(xsim, 1, 1);
								break;
							}
						}

						pos[0] = (0.5 + (double) xsim.coord(0)) / (double) sim.numpts;
						pos[1] = (0.5 + (double) xsim.coord(1)) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+1])(xlc) = (*Sij)(xsim, 0, 1);
								break;
							}
						}

						pos[1] = (double) xsim.coord(1) / (double) sim.numpts;
						pos[2] = (0.5 + (double) xsim.coord(2)) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+2])(xlc) = (*Sij)(xsim, 0, 2);
								break;
							}
						}

						pos[0] = (double) xsim.coord(0) / (double) sim.numpts;
						pos[1] = (0.5 + (double) xsim.coord(1)) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+4])(xlc) = (*Sij)(xsim, 1, 2);
								break;
							}
						}
#endif
#endif
						xlc.next();
					}
				}
			}

			n = findIntersectingLightcones(sim.lightcone[i], d - tau + (0.5 + maxvel) * dtau_old, d - tau - 0.5 * dtau, domain, vertex);

			if (sim.out_lightcone[i] & MASK_GADGET && sim.lightcone[i].distance[0] > d - tau + 0.5 * dtau_old && sim.lightcone[i].distance[1] <= d - tau + 0.5 * dtau_old && d - tau + 0.5 * dtau_old > 0.)
			{
				hdr.num_files = 1;
				hdr.Omega0 = cosmo.Omega_m;
				hdr.OmegaLambda = cosmo.Omega_Lambda;
				hdr.HubbleParam = cosmo.h;
				hdr.BoxSize = sim.boxsize / GADGET_LENGTH_CONVERSION;
				hdr.flag_sfr = 0;
				hdr.flag_cooling = 0;
				hdr.flag_feedback = 0;
				hdr.flag_age = 0;
				hdr.flag_metals = 0;
				for (p = 0; p < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8 - 2 * 4 - 6 * 4; p++)
					hdr.fill[p] = 0;
				for (p = 0; p < 6; p++)
				{
					hdr.npart[p] = 0;
					hdr.npartTotal[p] = 0;
					hdr.npartTotalHW[p] = 0;
					hdr.mass[p] = 0.;
				}

				hdr.time = a;
				hdr.redshift = (1./a) - 1.;
				
				if (sim.baryon_flag)
					hdr.mass[1] = (double) sim.tracer_factor[0] * C_RHO_CRIT * cosmo.Omega_cdm * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;
				else
					hdr.mass[1] = (double) sim.tracer_factor[0] * C_RHO_CRIT * (cosmo.Omega_cdm + cosmo.Omega_b) * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;

				if (sim.num_lightcone > 1)
					sprintf(filename, "%d_%04d", i, cycle);
				else
					sprintf(filename, "_%04d", cycle);

				pcls_cdm->saveGadget2(h5filename + filename + "_cdm", hdr, sim.lightcone[i], d - tau + (0.5 + maxvel) * dtau_old, d - tau - 0.5 * dtau, vertex, n, sim.tracer_factor[0]);

				if (sim.baryon_flag)
				{
					hdr.mass[1] = (double) sim.tracer_factor[1] * C_RHO_CRIT * cosmo.Omega_b * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[1] / GADGET_MASS_CONVERSION;
					pcls_b->saveGadget2(h5filename + filename + "_b", hdr, sim.lightcone[i], d - tau + (0.5 + maxvel) * dtau_old, d - tau - 0.5 * dtau, vertex, n, sim.tracer_factor[1]);
				}
				for (p = 0; p < cosmo.num_ncdm; p++)
				{
					sprintf(buffer, "_ncdm%d", p);
					hdr.mass[1] = (double) sim.tracer_factor[i+1+sim.baryon_flag] * C_RHO_CRIT * cosmo.Omega_ncdm[i] * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[i+1+sim.baryon_flag] / GADGET_MASS_CONVERSION;
					pcls_ncdm[p].saveGadget2(h5filename + filename + buffer, hdr, sim.lightcone[i], d - tau + (0.5 + maxvel) * dtau_old, d - tau - 0.5 * dtau, vertex, n, sim.tracer_factor[i+1+sim.baryon_flag]);
				}
			}
		}
		else if (parallel.isRoot() && outfile != NULL)
			fclose(outfile);
	}
}


//////////////////////////
// writeSpectra
//////////////////////////
// Description:
//   output of spectra
// 
// Arguments:
//   sim            simulation metadata structure
//   cosmo          cosmological parameter structure
//   fourpiG        4 pi G (in code units)
//   a              scale factor
//   pkcount        spectrum output index
//   pcls_cdm       pointer to (uninitialized) particle handler for CDM
//   pcls_b         pointer to (uninitialized) particle handler for baryons
//   pcls_ncdm      array of (uninitialized) particle handlers for
//                  non-cold DM (may be set to NULL)
//   phi            pointer to allocated field
//   chi            pointer to allocated field
//   Bi             pointer to allocated field
//   source         pointer to allocated field
//   Sij            pointer to allocated field
//   scalarFT       pointer to allocated field
//   BiFT           pointer to allocated field
//   SijFT          pointer to allocated field
//   plan_phi       pointer to FFT planner
//   plan_chi       pointer to FFT planner
//   plan_Bi        pointer to FFT planner
//   plan_source    pointer to FFT planner
//   plan_Sij       pointer to FFT planner
//   Bi_check       pointer to allocated field (or NULL)
//   BiFT_check     pointer to allocated field (or NULL)
//   plan_Bi_check  pointer to FFT planner (or NULL)
//
// Returns:
// 
//////////////////////////

void writeSpectra(metadata & sim, cosmology & cosmo, const double fourpiG, const double a, const int pkcount,
#ifdef HAVE_CLASS
background & class_background, perturbs & class_perturbs, spectra & class_spectra, icsettings & ic,
#endif
Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, Field<Real> * phi, Field<Real> * chi, Field<Real> * Bi, Field<Real> * source, Field<Real> * Sij, Field<Cplx> * scalarFT, Field<Cplx> * BiFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_phi, PlanFFT<Cplx> * plan_chi, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_source, PlanFFT<Cplx> * plan_Sij, Field<Real> * Bi_check = NULL, Field<Cplx> * BiFT_check = NULL, PlanFFT<Cplx> * plan_Bi_check = NULL)
{
	char filename[2*PARAM_MAX_LENGTH+24];
	char buffer[64];
	int i, j;
	Site x(phi->lattice());
	rKSite kFT(scalarFT->lattice());
	long numpts3d = (long) sim.numpts * (long) sim.numpts * (long) sim.numpts;
	Cplx tempk;
	double Omega_ncdm;

	Real * kbin;
	Real * power;
	Real * kscatter;
	Real * pscatter;
	int * occupation;

	kbin = (Real *) malloc(sim.numbins * sizeof(Real));
	power = (Real *) malloc(sim.numbins * sizeof(Real));
	kscatter = (Real *) malloc(sim.numbins * sizeof(Real));
	pscatter = (Real *) malloc(sim.numbins * sizeof(Real));
	occupation = (int *) malloc(sim.numbins * sizeof(int));

	if (sim.out_pk & MASK_RBARE || sim.out_pk & MASK_DBARE || sim.out_pk & MASK_POT || ((sim.out_pk & MASK_T00 || sim.out_pk & MASK_DELTA) && sim.gr_flag == 0))
	{
		projection_init(source);
		scalarProjectionCIC_project(pcls_cdm, source);
		if (sim.baryon_flag)
			scalarProjectionCIC_project(pcls_b, source);
		for (i = 0; i < cosmo.num_ncdm; i++)
			scalarProjectionCIC_project(pcls_ncdm+i, source);
		scalarProjectionCIC_comm(source);
		plan_source->execute(FFT_FORWARD);
				
		if (sim.out_pk & MASK_RBARE || sim.out_pk & MASK_DBARE || ((sim.out_pk & MASK_T00 || sim.out_pk & MASK_DELTA) && sim.gr_flag == 0))
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				
		if (sim.out_pk & MASK_RBARE)
		{
			sprintf(filename, "%s%s%03d_rhoN.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of rho_N", a, sim.z_pk[pkcount]);
		}

		if (sim.out_pk & MASK_DBARE)
		{
			sprintf(filename, "%s%s%03d_deltaN.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_m * cosmo.Omega_m, filename, "power spectrum of delta_N", a, sim.z_pk[pkcount]);
		}
				
		if (sim.out_pk & MASK_T00 && sim.gr_flag == 0)
		{
			sprintf(filename, "%s%s%03d_T00.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00", a, sim.z_pk[pkcount]);
		}

		if (sim.out_pk & MASK_DELTA && sim.gr_flag == 0)
		{
			sprintf(filename, "%s%s%03d_delta.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_m * cosmo.Omega_m, filename, "power spectrum of delta", a, sim.z_pk[pkcount]);
		}
				
		if (sim.out_pk & MASK_POT)
		{
			solveModifiedPoissonFT(*scalarFT, *scalarFT, fourpiG / a);
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
			sprintf(filename, "%s%s%03d_psiN.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of psi_N", a, sim.z_pk[pkcount]);
		}
				
		if ((cosmo.num_ncdm > 0 || sim.baryon_flag) && (sim.out_pk & MASK_DBARE || (sim.out_pk & MASK_DELTA && sim.gr_flag == 0)))
		{
			projection_init(source);
			scalarProjectionCIC_project(pcls_cdm, source);
			scalarProjectionCIC_comm(source);
			plan_source->execute(FFT_FORWARD);
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
			sprintf(filename, "%s%s%03d_cdm.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (sim.baryon_flag ? (cosmo.Omega_cdm * cosmo.Omega_cdm) : ((cosmo.Omega_cdm + cosmo.Omega_b) * (cosmo.Omega_cdm + cosmo.Omega_b))), filename, "power spectrum of delta_N for cdm", a, sim.z_pk[pkcount]);
			if (sim.baryon_flag)
			{
				// store k-space information for cross-spectra using SijFT as temporary array
				if (sim.out_pk & MASK_XSPEC)
				{
					for (kFT.first(); kFT.test(); kFT.next())
						(*SijFT)(kFT, 0) = (*scalarFT)(kFT);
				}
				projection_init(source);
				scalarProjectionCIC_project(pcls_b, source);
				scalarProjectionCIC_comm(source);
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				sprintf(filename, "%s%s%03d_b.dat", sim.output_path, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_b * cosmo.Omega_b, filename, "power spectrum of delta_N for baryons", a, sim.z_pk[pkcount]);
				if (sim.out_pk & MASK_XSPEC)
				{
					extractCrossSpectrum(*scalarFT, *SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
					sprintf(filename, "%s%s%03d_cdmxb.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_cdm * cosmo.Omega_b, filename, "cross power spectrum of delta_N for cdm x baryons", a, sim.z_pk[pkcount]);
				}
			}
			Omega_ncdm = 0.;
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				projection_init(source);
				scalarProjectionCIC_project(pcls_ncdm+i, source);
				scalarProjectionCIC_comm(source);
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				sprintf(filename, "%s%s%03d_ncdm%d.dat", sim.output_path, sim.basename_pk, pkcount, i);
				sprintf(buffer, "power spectrum of delta_N for ncdm %d", i);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_ncdm[i] * cosmo.Omega_ncdm[i], filename, buffer, a, sim.z_pk[pkcount]);
				Omega_ncdm += cosmo.Omega_ncdm[i];
				// store k-space information for cross-spectra using SijFT as temporary array
				if (cosmo.num_ncdm > 1 && i < 6)
				{
					for (kFT.first(); kFT.test(); kFT.next())
						(*SijFT)(kFT, i) = (*scalarFT)(kFT);
				}						
			}
			if (cosmo.num_ncdm > 1 && cosmo.num_ncdm <= 7)
			{
				for (kFT.first(); kFT.test(); kFT.next())
				{
					for (i = 0; i < cosmo.num_ncdm-1; i++)
						(*scalarFT)(kFT) += (*SijFT)(kFT, i);
				}
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				sprintf(filename, "%s%s%03d_ncdm.dat", sim.output_path, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * Omega_ncdm * Omega_ncdm, filename, "power spectrum of delta_N for total ncdm", a, sim.z_pk[pkcount]);
			}
			if (cosmo.num_ncdm > 1)
			{
				for (i = 0; i < cosmo.num_ncdm-1 && i < 5; i++)
				{
					for (j = i+1; j < cosmo.num_ncdm && j < 6; j++)
					{
						if (sim.out_pk & MASK_XSPEC || (i == 0 && j == 1) || (i == 2 && j == 3) || (i == 4 && j == 5))
						{
							extractCrossSpectrum(*SijFT, *SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR, i, j);
							sprintf(filename, "%s%s%03d_ncdm%dx%d.dat", sim.output_path, sim.basename_pk, pkcount, i, j);
							sprintf(buffer, "cross power spectrum of delta_N for ncdm %d x %d", i, j);
							writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_ncdm[i] * cosmo.Omega_ncdm[j], filename, buffer, a, sim.z_pk[pkcount]);
						}
					}
				}
			}
		}
	}
	
	if (sim.out_pk & MASK_PHI)
	{
		plan_phi->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_phi.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of phi", a, sim.z_pk[pkcount]);
	}
			
	if (sim.out_pk & MASK_CHI)
	{
		plan_chi->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_chi.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of chi", a, sim.z_pk[pkcount]);
	}
			
	if (sim.out_pk & MASK_HIJ)
	{
		projection_init(Sij);
		projection_Tij_project(pcls_cdm, Sij, a, phi);
		if (sim.baryon_flag)
			projection_Tij_project(pcls_b, Sij, a, phi);
		for (i = 0; i < cosmo.num_ncdm; i++)
			projection_Tij_project(pcls_ncdm+i, Sij, a, phi);
		projection_Tij_comm(Sij);

		prepareFTsource<Real>(*phi, *Sij, *Sij, 2. * fourpiG / (double) sim.numpts / (double) sim.numpts / a);
		plan_Sij->execute(FFT_FORWARD);
		projectFTtensor(*SijFT, *SijFT);

		extractPowerSpectrum(*SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_hij.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, 2. * M_PI * M_PI, filename, "power spectrum of hij", a, sim.z_pk[pkcount]);
	}
			
	if ((sim.out_pk & MASK_T00 || sim.out_pk & MASK_DELTA) && sim.gr_flag > 0)
	{
		projection_init(source);
#ifdef HAVE_CLASS
		if (sim.radiation_flag > 0 || sim.fluid_flag > 0)
		{
			projection_T00_project(class_background, class_perturbs, class_spectra, *source, *scalarFT, plan_source, sim, ic, cosmo, fourpiG, a);
			if (sim.out_pk & MASK_DELTA)
			{
				Omega_ncdm = 0;
				for (i = 0; i < cosmo.num_ncdm; i++)
				{
					if (a < 1. / (sim.z_switch_deltancdm[i] + 1.) && cosmo.Omega_ncdm[i] > 0)
						Omega_ncdm += bg_ncdm(a, cosmo, i);
				}
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				sprintf(filename, "%s%s%03d_deltaclass.dat", sim.output_path, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * ((a < 1. / (sim.z_switch_deltarad + 1.) ? sim.radiation_flag : 0) * cosmo.Omega_rad / a + sim.fluid_flag * cosmo.Omega_fld / pow(a, 3. * cosmo.w0_fld) + Omega_ncdm) * ((a < 1. / (sim.z_switch_deltarad + 1.) ? sim.radiation_flag : 0) * cosmo.Omega_rad / a + sim.fluid_flag * cosmo.Omega_fld / pow(a, 3. * cosmo.w0_fld) + Omega_ncdm), filename, "power spectrum of delta for linear fields (CLASS)", a, sim.z_pk[pkcount]);
			}
		}
#endif
		projection_T00_project(pcls_cdm, source, a, phi);
		if (sim.baryon_flag)
			projection_T00_project(pcls_b, source, a, phi);
		for (i = 0; i < cosmo.num_ncdm; i++)
			projection_T00_project(pcls_ncdm+i, source, a, phi);
		projection_T00_comm(source);

		plan_source->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
		
		if (sim.out_pk & MASK_T00)
		{
			sprintf(filename, "%s%s%03d_T00.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00", a, sim.z_pk[pkcount]);
		}

		if (sim.out_pk & MASK_DELTA)
		{
			sprintf(filename, "%s%s%03d_delta.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)), filename, "power spectrum of delta", a, sim.z_pk[pkcount]);
		}
				
		if (cosmo.num_ncdm > 0 || sim.baryon_flag || sim.radiation_flag > 0 || sim.fluid_flag > 0)
		{
			projection_init(source);
			projection_T00_project(pcls_cdm, source, a, phi);
			projection_T00_comm(source);
			plan_source->execute(FFT_FORWARD);
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
			if (sim.out_pk & MASK_T00)
			{
				sprintf(filename, "%s%s%03d_T00cdm.dat", sim.output_path, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00 for cdm", a, sim.z_pk[pkcount]);
			}
			if (sim.out_pk & MASK_DELTA)
			{
				sprintf(filename, "%s%s%03d_deltacdm.dat", sim.output_path, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (sim.baryon_flag ? (cosmo.Omega_cdm * cosmo.Omega_cdm) : ((cosmo.Omega_cdm + cosmo.Omega_b) * (cosmo.Omega_cdm + cosmo.Omega_b))), filename, "power spectrum of delta for cdm", a, sim.z_pk[pkcount]);
			}
			if (sim.baryon_flag)
			{
				// store k-space information for cross-spectra using SijFT as temporary array
				if (sim.out_pk & MASK_XSPEC)
				{
					for (kFT.first(); kFT.test(); kFT.next())
						(*SijFT)(kFT, 0) = (*scalarFT)(kFT);
				}
				projection_init(source);
				projection_T00_project(pcls_b, source, a, phi);
				projection_T00_comm(source);
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				if (sim.out_pk & MASK_T00)
				{
					sprintf(filename, "%s%s%03d_T00b.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00 for baryons", a, sim.z_pk[pkcount]);
				}
				if (sim.out_pk & MASK_DELTA)
				{
					sprintf(filename, "%s%s%03d_deltab.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_b * cosmo.Omega_b, filename, "power spectrum of delta for baryons", a, sim.z_pk[pkcount]);
				}
				if (sim.out_pk & MASK_XSPEC)
				{
					extractCrossSpectrum(*scalarFT, *SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
					sprintf(filename, "%s%s%03d_deltacdmxb.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_b * cosmo.Omega_cdm, filename, "cross power spectrum of delta for cdm x baryons", a, sim.z_pk[pkcount]);
				}
			}
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				projection_init(source);
				projection_T00_project(pcls_ncdm+i, source, a, phi);
				projection_T00_comm(source);
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				if (sim.out_pk & MASK_T00)
				{
					sprintf(filename, "%s%s%03d_T00ncdm%d.dat", sim.output_path, sim.basename_pk, pkcount, i);
					sprintf(buffer, "power spectrum of T00 for ncdm %d", i);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, buffer, a, sim.z_pk[pkcount]);
				}
				if (sim.out_pk & MASK_DELTA)
				{
					sprintf(filename, "%s%s%03d_deltancdm%d.dat", sim.output_path, sim.basename_pk, pkcount, i);
					sprintf(buffer, "power spectrum of delta for ncdm %d", i);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * bg_ncdm(a, cosmo, i) * bg_ncdm(a, cosmo, i), filename, buffer, a, sim.z_pk[pkcount]);
				}					
				// store k-space information for cross-spectra using SijFT as temporary array
				if (cosmo.num_ncdm > 1 && i < 6)
				{
					for (kFT.first(); kFT.test(); kFT.next())
						(*SijFT)(kFT, i) = (*scalarFT)(kFT);
				}
			}
			if (cosmo.num_ncdm > 1 && cosmo.num_ncdm <= 7)
			{
				for (kFT.first(); kFT.test(); kFT.next())
				{
					for (i = 0; i < cosmo.num_ncdm-1; i++)
						(*scalarFT)(kFT) += (*SijFT)(kFT, i);
				}
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				if (sim.out_pk & MASK_T00)
				{
					sprintf(filename, "%s%s%03d_T00ncdm.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00 for total ncdm", a, sim.z_pk[pkcount]);
				}
				if (sim.out_pk & MASK_DELTA)
				{
					sprintf(filename, "%s%s%03d_deltancdm.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * bg_ncdm(a, cosmo) * bg_ncdm(a, cosmo), filename, "power spectrum of delta for total ncdm", a, sim.z_pk[pkcount]);
				}
			}
			if (cosmo.num_ncdm > 1)
			{
				for (i = 0; i < cosmo.num_ncdm-1 && i < 5; i++)
				{
					for (j = i+1; j < cosmo.num_ncdm && j < 6; j++)
					{
						if (sim.out_pk & MASK_XSPEC || (i == 0 && j == 1) || (i == 2 && j == 3) || (i == 4 && j == 5))
						{
							extractCrossSpectrum(*SijFT, *SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR, i, j);
							if (sim.out_pk & MASK_T00)
							{
								sprintf(filename, "%s%s%03d_T00ncdm%dx%d.dat", sim.output_path, sim.basename_pk, pkcount, i, j);
								sprintf(buffer, "cross power spectrum of T00 for ncdm %d x %d", i, j);
								writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, buffer, a, sim.z_pk[pkcount]);
							}
							if (sim.out_pk & MASK_DELTA)
							{
								sprintf(filename, "%s%s%03d_deltancdm%dx%d.dat", sim.output_path, sim.basename_pk, pkcount, i, j);
								sprintf(buffer, "cross power spectrum of delta for ncdm %d x %d", i, j);
								writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * bg_ncdm(a, cosmo, i) * bg_ncdm(a, cosmo, j), filename, buffer, a, sim.z_pk[pkcount]);
							}
						}
					}
				}
			}
		}
	}
			
	if (sim.out_pk & MASK_B)
	{
		extractPowerSpectrum(*BiFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_B.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, a * a * a * a * sim.numpts * sim.numpts * 2. * M_PI * M_PI, filename, "power spectrum of B", a, sim.z_pk[pkcount]);
			
#ifdef CHECK_B
		if (sim.vector_flag == VECTOR_PARABOLIC)
		{
			projection_init(Bi_check);
			projection_T0i_project(pcls_cdm, Bi_check, phi);
			if (sim.baryon_flag)
				projection_T0i_project(pcls_b, Bi_check, phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
				projection_T0i_project(pcls_ncdm+i, Bi_check, phi);
			projection_T0i_comm(Bi_check);
			plan_Bi_check->execute(FFT_FORWARD);
			projectFTvector(*BiFT_check, *BiFT_check, fourpiG / (double) sim.numpts / (double) sim.numpts);
		}
		extractPowerSpectrum(*BiFT_check, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_B_check.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, a * a * a * a * sim.numpts * sim.numpts * 2. * M_PI * M_PI, filename, "power spectrum of B", a, sim.z_pk[pkcount]);
#endif
	}

	free(kbin);
	free(power);
	free(kscatter);
	free(pscatter);
	free(occupation);
}

#endif

