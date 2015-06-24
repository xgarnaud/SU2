/*!
 * \file numerics_structure.inl
 * \brief In-Line subroutines of the <i>numerics_structure.hpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 4.0.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

inline double CNumerics::Determinant_3x3(double A00, double A01, double A02, double A10, double A11, double A12, double A20, double A21, double A22) {
	return A00*(A11*A22-A12*A21) - A01*(A10*A22-A12*A20) + A02*(A10*A21-A11*A20);
}

inline void CNumerics::SetFEA_StiffMatrix2D(double **StiffMatrix_Elem, double CoordCorners[8][3], unsigned short nNodes, unsigned short form2d) { }

inline void CNumerics::SetFEA_StiffMatrix3D(double **StiffMatrix_Elem, double CoordCorners[8][3], unsigned short nNodes) { }

inline void CNumerics::SetFEA_StiffMassMatrix2D(double **StiffMatrix_Elem, double **MassMatrix_Elem, double CoordCorners[8][3], unsigned short nNodes, unsigned short form2d) { }

inline void CNumerics::SetFEA_StiffMassMatrix3D(double **StiffMatrix_Elem, double **MassMatrix_Elem, double CoordCorners[8][3], unsigned short nNodes) { }

inline void CNumerics::GetFEA_StressNodal2D(double StressVector[8][3], double DispElement[8], double CoordCorners[8][3], unsigned short nNodes, unsigned short form2d) { }

inline void CNumerics::GetFEA_StressNodal3D(double StressVector[8][6], double DispElement[24], double CoordCorners[8][3], unsigned short nNodes) { }

inline void CNumerics::SetFEA_DeadLoad2D(double *DeadLoadVector_Elem, double CoordCorners[8][3], unsigned short nNodes, double matDensity) { }

inline void CNumerics::SetFEA_DeadLoad3D(double *DeadLoadVector_Elem, double CoordCorners[8][3], unsigned short nNodes, double matDensity) { }

inline void CNumerics::PressInt_Linear(double CoordCorners[4][3], double *tn_e, double *Fnodal) { }

inline void CNumerics::ViscTermInt_Linear(double CoordCorners[2][2], double Tau_0[3][3], double Tau_1[3][3], double FviscNodal[4]) { }

inline void CNumerics::ComputeResidual(double *val_residual, CConfig *config) { }

inline void CNumerics::ComputeResidual(double *val_residual_i, double *val_residual_j) { }

inline void CNumerics::ComputeResidual(double *val_residual_i, double *val_residual_j, CConfig *config) { }

inline void CNumerics::ComputeResidual(double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) { }

inline void CNumerics::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, 
                                   CConfig *config) { }

inline void CNumerics::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j,
                                   double **val_JacobianMeanFlow_i, double **val_JacobianMeanFlow_j, CConfig *config) { }

inline void CNumerics::ComputeResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i, 
								   double **val_Jacobian_j, CConfig *config) { }

inline void CNumerics::ComputeResidual(double *val_residual_i, double *val_residual_j, 
								   double **val_Jacobian_ii, double **val_Jacobian_ij, 
								   double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config) { }
							
inline void CNumerics::ComputeResidual(double *val_resconv_i, double *val_resvisc_i, double *val_resconv_j, 
								   double *val_resvisc_j, double **val_Jacobian_ii, double **val_Jacobian_ij, 
								   double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config) { }
							
inline void CNumerics::ComputeResidual(double **val_stiffmatrix_elem, CConfig *config) { }

inline void CNumerics::GetEq_Rxn_Coefficients(double **EqnRxnConstants, CConfig *config) { };
														
inline void CNumerics::ComputeResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) { }

inline void CNumerics::ComputeResidual_TransLM(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config, double &gamma_sep) {}

inline void CNumerics::ComputeResidual_Axisymmetric(double *val_residual, CConfig *config) { }

inline void CNumerics::ComputeResidual_Axisymmetric_ad(double *val_residual, double *val_residuald, CConfig *config) { }

inline void CNumerics::SetJacobian_Axisymmetric(double **val_Jacobian_i, CConfig *config) { }

inline void CNumerics::ComputeVibRelaxation(double *val_residual, double **val_Jacobian_i, CConfig *config) { }

inline void CNumerics::ComputeChemistry(double *val_residual, double **val_Jacobian_i, CConfig *config) { }

inline void CNumerics::GetKeqConstants(double *A, unsigned short val_reaction, CConfig *config) { }

inline void CNumerics::ComputeSourceViscous(double *val_residual, CConfig *config) { }

inline double CNumerics::GetPrecond_Beta() { return 0; }

inline void CNumerics::SetRhosIndex(unsigned short val_Index) { RHOS_INDEX = val_Index; }

inline void CNumerics::SetRhoIndex(unsigned short val_Index) { RHO_INDEX = val_Index; }

inline void CNumerics::SetPIndex(unsigned short val_Index) { P_INDEX = val_Index; }

inline void CNumerics::SetTIndex(unsigned short val_Index) { T_INDEX = val_Index; }
  
inline void CNumerics::SetTveIndex(unsigned short val_Index) { TVE_INDEX = val_Index; }

inline void CNumerics::SetVelIndex(unsigned short val_Index) { VEL_INDEX = val_Index; }
  
inline void CNumerics::SetHIndex(unsigned short val_Index) { H_INDEX = val_Index; }
  
inline void CNumerics::SetAIndex(unsigned short val_Index) { A_INDEX = val_Index; }
  
inline void CNumerics::SetRhoCvtrIndex(unsigned short val_Index) { RHOCVTR_INDEX = val_Index; }

inline void CNumerics::SetRhoCvveIndex(unsigned short val_Index) { RHOCVVE_INDEX = val_Index; }

inline void CNumerics::SetdPdU(double *val_dPdU_i, double *val_dPdU_j) { dPdU_i = val_dPdU_i; dPdU_j = val_dPdU_j; }

inline void CNumerics::SetdTdU(double *val_dTdU_i, double *val_dTdU_j) { dTdU_i = val_dTdU_i; dTdU_j = val_dTdU_j; }

inline void CNumerics::SetdTvedU(double *val_dTvedU_i, double *val_dTvedU_j) { dTvedU_i = val_dTvedU_i; dTvedU_j = val_dTvedU_j; }
	
inline void CNumerics::SetUndivided_Laplacian(double *val_und_lapl_i, double *val_und_lapl_j) {
	Und_Lapl_i = val_und_lapl_i; 
	Und_Lapl_j = val_und_lapl_j; 
}

inline void CNumerics::SetSensor( double val_sensor_i, double val_sensor_j) {
	Sensor_i = val_sensor_i;
	Sensor_j = val_sensor_j;
}

inline void CNumerics::SetConservative(double *val_u_i, double *val_u_j) {
	U_i = val_u_i;
	U_j = val_u_j;
}

inline void CNumerics::SetConservative_ZeroOrder(double *val_u_i, double *val_u_j) {
	UZeroOrder_i = val_u_i;
	UZeroOrder_j = val_u_j;
}

inline void CNumerics::SetPrimitive(double *val_v_i, double *val_v_j) {
	V_i = val_v_i;
	V_j = val_v_j;
}

inline void CNumerics::SetSecondary(double *val_s_i, double *val_s_j) {
  	S_i = val_s_i;
  	S_j = val_s_j;
}

inline void CNumerics::SetConservative(double *val_u_0, double *val_u_1, double *val_u_2) {
	U_0 = val_u_0;
	U_1 = val_u_1;
	U_2 = val_u_2;
}

inline void CNumerics::SetConservative(double *val_u_0, double *val_u_1, double *val_u_2, double *val_u_3) {
	U_0 = val_u_0;
	U_1 = val_u_1;
	U_2 = val_u_2;
	U_3 = val_u_3;
}

inline void CNumerics::SetVelocity2_Inf(double velocity2) {
	vel2_inf = velocity2;
}

inline void CNumerics::SetVorticity(double *val_vorticity_i, double *val_vorticity_j) {
  Vorticity_i = val_vorticity_i;
  Vorticity_j = val_vorticity_j;
}

inline void CNumerics::SetStrainMag(double val_strainmag_i, double val_strainmag_j) {
  StrainMag_i = val_strainmag_i;
  StrainMag_j = val_strainmag_j;
}

inline void CNumerics::SetTimeStep(double val_timestep) {TimeStep = val_timestep;}

inline void CNumerics::SetLaminarViscosity(double val_lam_viscosity_i, double val_lam_viscosity_j) {
	Laminar_Viscosity_i = val_lam_viscosity_i;
	Laminar_Viscosity_j = val_lam_viscosity_j;
}

inline void CNumerics::SetThermalConductivity(double val_therm_conductivity_i, double val_therm_conductivity_j) {
	Thermal_Conductivity_i = val_therm_conductivity_i;
	Thermal_Conductivity_j = val_therm_conductivity_j;
}

inline void CNumerics::SetThermalConductivity_ve(double val_therm_conductivity_ve_i, double val_therm_conductivity_ve_j) {
	Thermal_Conductivity_ve_i = val_therm_conductivity_ve_i;
	Thermal_Conductivity_ve_j = val_therm_conductivity_ve_j;
}

inline void CNumerics::SetDiffusionCoeff(double* val_diffusioncoeff_i, double* val_diffusioncoeff_j) {
	Diffusion_Coeff_i = val_diffusioncoeff_i;
	Diffusion_Coeff_j = val_diffusioncoeff_j;
}

inline void CNumerics::SetEddyViscosity(double val_eddy_viscosity_i, double val_eddy_viscosity_j) {
	Eddy_Viscosity_i = val_eddy_viscosity_i;
	Eddy_Viscosity_j = val_eddy_viscosity_j;
}

inline void CNumerics::SetIntermittency(double intermittency_in) { }

inline void CNumerics::SetProduction(double val_production) { }

inline void CNumerics::SetDestruction(double val_destruction) { }

inline void CNumerics::SetCrossProduction(double val_crossproduction) { }

inline double CNumerics::GetProduction(void) { return 0; }

inline double CNumerics::GetDestruction(void) { return 0; }

inline double CNumerics::GetCrossProduction(void) { return 0; }

inline void CNumerics::SetTurbKineticEnergy(double val_turb_ke_i, double val_turb_ke_j) {
	turb_ke_i = val_turb_ke_i;
	turb_ke_j = val_turb_ke_j;
}

inline void CNumerics::SetDistance(double val_dist_i, double val_dist_j) {
	dist_i = val_dist_i;
	dist_j = val_dist_j;
}

inline void CNumerics::SetAdjointVar(double *val_psi_i, double *val_psi_j) {
	Psi_i = val_psi_i;
	Psi_j = val_psi_j;
}

inline void CNumerics::SetLinearizedVar(double *val_deltau_i, double *val_deltau_j) {
	DeltaU_i = val_deltau_i;
	DeltaU_j = val_deltau_j;
}

inline void CNumerics::SetAdjointVarGradient(double **val_psivar_grad_i, double **val_psivar_grad_j) {
	PsiVar_Grad_i = val_psivar_grad_i;
	PsiVar_Grad_j = val_psivar_grad_j;
}

inline void CNumerics::SetAdjointVarLimiter(double *val_psivar_lim_i, double *val_psivar_lim_j) {
	PsiVar_Lim_i = val_psivar_lim_i;
	PsiVar_Lim_j = val_psivar_lim_j;
}

inline void CNumerics::SetTurbVar(double *val_turbvar_i, double *val_turbvar_j) {
	TurbVar_i = val_turbvar_i;
	TurbVar_j = val_turbvar_j;
}

inline void CNumerics::SetTransVar(double *val_transvar_i, double *val_transvar_j) {
	TransVar_i = val_transvar_i;
	TransVar_j = val_transvar_j;
}

inline void CNumerics::SetTurbVarGradient(double **val_turbvar_grad_i, double **val_turbvar_grad_j) {
	TurbVar_Grad_i = val_turbvar_grad_i;
	TurbVar_Grad_j = val_turbvar_grad_j;
}

inline void CNumerics::SetTransVarGradient(double **val_transvar_grad_i, double **val_transvar_grad_j) {
	TransVar_Grad_i = val_transvar_grad_i;
	TransVar_Grad_j = val_transvar_grad_j;
}

inline void CNumerics::SetLevelSetVar(double *val_levelsetvar_i, double *val_levelsetvar_j) {
	LevelSetVar_i = val_levelsetvar_i;
	LevelSetVar_j = val_levelsetvar_j;
}

inline void CNumerics::SetLevelSetVarGradient(double **val_levelsetvar_grad_i, double **val_levelsetvar_grad_j) {
	LevelSetVar_Grad_i = val_levelsetvar_grad_i;
	LevelSetVar_Grad_j = val_levelsetvar_grad_j;
}

inline void CNumerics::SetPrimVarGradient(double **val_primvar_grad_i, double **val_primvar_grad_j) {
	PrimVar_Grad_i = val_primvar_grad_i;
	PrimVar_Grad_j = val_primvar_grad_j;
}

inline void CNumerics::SetPrimVarLimiter(double *val_primvar_lim_i, double *val_primvar_lim_j) {
  PrimVar_Lim_i = val_primvar_lim_i;
  PrimVar_Lim_j = val_primvar_lim_j;
}

inline void CNumerics::SetConsVarGradient(double **val_consvar_grad_i, double **val_consvar_grad_j) {
	ConsVar_Grad_i = val_consvar_grad_i;
	ConsVar_Grad_j = val_consvar_grad_j;
}

inline void CNumerics::SetConsVarGradient(double **val_consvar_grad_0, double **val_consvar_grad_1, double **val_consvar_grad_2) {
	ConsVar_Grad_0 = val_consvar_grad_0;
	ConsVar_Grad_1 = val_consvar_grad_1;
	ConsVar_Grad_2 = val_consvar_grad_2;
}

inline void CNumerics::SetConsVarGradient(double **val_consvar_grad_0, double **val_consvar_grad_1, double **val_consvar_grad_2, double **val_consvar_grad_3) {
	ConsVar_Grad_0 = val_consvar_grad_0;
	ConsVar_Grad_1 = val_consvar_grad_1;
	ConsVar_Grad_2 = val_consvar_grad_2;
	ConsVar_Grad_3 = val_consvar_grad_3;
}

inline void CNumerics::SetConsVarGradient(double **val_consvar_grad) {
	ConsVar_Grad = val_consvar_grad;
}

inline void CNumerics::SetCoord(double *val_coord_i, double *val_coord_j) {
	Coord_i = val_coord_i;
	Coord_j = val_coord_j;
}

inline void CNumerics::SetCoord(double *val_coord_0, double *val_coord_1, 
									 double *val_coord_2) {
	Coord_0 = val_coord_0;
	Coord_1 = val_coord_1;
	Coord_2 = val_coord_2;
}

inline void CNumerics::SetCoord(double *val_coord_0, double *val_coord_1, 
									 double *val_coord_2, double *val_coord_3) {
	Coord_0 = val_coord_0;
	Coord_1 = val_coord_1;
	Coord_2 = val_coord_2;
	Coord_3 = val_coord_3;	
}

inline void CNumerics::SetGridVel(double *val_gridvel_i, double *val_gridvel_j) {
	GridVel_i = val_gridvel_i;
	GridVel_j = val_gridvel_j;
}

inline void CNumerics::SetWindGust(double *val_windgust_i, double *val_windgust_j) {
	WindGust_i = val_windgust_i;
	WindGust_j = val_windgust_j;
}

inline void CNumerics::SetWindGustDer(double *val_windgustder_i, double *val_windgustder_j) {
	WindGustDer_i = val_windgustder_i;
	WindGustDer_j = val_windgustder_j;
}

inline void CNumerics::SetPressure(double val_pressure_i, double val_pressure_j) {
	Pressure_i = val_pressure_i;
	Pressure_j = val_pressure_j;
}

inline void CNumerics::SetDensityInc(double val_densityinc_i, double val_densityinc_j) {
	DensityInc_i = val_densityinc_i;
	DensityInc_j = val_densityinc_j;
}

inline void CNumerics::SetBetaInc2(double val_betainc2_i, double val_betainc2_j) {
	BetaInc2_i = val_betainc2_i;
	BetaInc2_j = val_betainc2_j;
}

inline void CNumerics::SetSoundSpeed(double val_soundspeed_i, double val_soundspeed_j) {
	SoundSpeed_i = val_soundspeed_i;
	SoundSpeed_j = val_soundspeed_j;
}

inline void CNumerics::SetEnthalpy(double val_enthalpy_i, double val_enthalpy_j) {
	Enthalpy_i = val_enthalpy_i;
	Enthalpy_j = val_enthalpy_j;
}

inline void CNumerics::SetLambda(double val_lambda_i, double val_lambda_j) {
	Lambda_i = val_lambda_i;
	Lambda_j = val_lambda_j;
}

inline void CNumerics::SetNeighbor(unsigned short val_neighbor_i, unsigned short val_neighbor_j) {
	Neighbor_i = val_neighbor_i;
	Neighbor_j = val_neighbor_j;
}

inline void CNumerics::SetTurbAdjointVar(double *val_turbpsivar_i, double *val_turbpsivar_j) {
	TurbPsi_i = val_turbpsivar_i;
	TurbPsi_j = val_turbpsivar_j;
}

inline void CNumerics::SetTurbAdjointGradient(double **val_turbpsivar_grad_i, double **val_turbpsivar_grad_j) {
	TurbPsi_Grad_i = val_turbpsivar_grad_i;
	TurbPsi_Grad_j = val_turbpsivar_grad_j;
}

inline void CNumerics::SetTemperature(double val_temp_i, double val_temp_j) {
	Temp_i = val_temp_i;
	Temp_j = val_temp_j;
}

inline void CNumerics::SetAuxVarGrad(double *val_auxvargrad_i, double *val_auxvargrad_j) {
	AuxVar_Grad_i = val_auxvargrad_i;
	AuxVar_Grad_j = val_auxvargrad_j;
}

inline void CNumerics::SetNormal(double *val_normal) { Normal = val_normal; }

inline void CNumerics::SetVolume(double val_volume) { Volume = val_volume; }

inline void CSourcePieceWise_TurbSST::SetF1blending(double val_F1_i, double val_F1_j) { 
	F1_i = val_F1_i; 
	F1_j = val_F1_j;
}

inline void CSourcePieceWise_TurbSST::SetF2blending(double val_F2_i, double val_F2_j) { 
	F2_i = val_F2_i; 
	F2_j = val_F2_j;
}

inline void CSourcePieceWise_TurbSST::SetCrossDiff(double val_CDkw_i, double val_CDkw_j) {
	CDkw_i = val_CDkw_i;
  CDkw_j = val_CDkw_j;
}			

inline void CSourcePieceWise_TurbSA::SetIntermittency(double intermittency_in) { intermittency = intermittency_in; }

inline void CSourcePieceWise_TurbSA::SetProduction(double val_production) { Production = val_production; }

inline void CSourcePieceWise_TurbSA::SetDestruction(double val_destruction) { Destruction = val_destruction; }

inline void CSourcePieceWise_TurbSA::SetCrossProduction(double val_crossproduction) { CrossProduction = val_crossproduction; }

inline double CSourcePieceWise_TurbSA::GetProduction(void) { return Production; }

inline double CSourcePieceWise_TurbSA::GetDestruction(void) { return Destruction; }

inline double CSourcePieceWise_TurbSA::GetCrossProduction(void) { return CrossProduction; }

inline void CSourcePieceWise_TurbSA_Neg::SetIntermittency(double intermittency_in) { intermittency = intermittency_in; }

inline void CSourcePieceWise_TurbSA_Neg::SetProduction(double val_production) { Production = val_production; }

inline void CSourcePieceWise_TurbSA_Neg::SetDestruction(double val_destruction) { Destruction = val_destruction; }

inline void CSourcePieceWise_TurbSA_Neg::SetCrossProduction(double val_crossproduction) { CrossProduction = val_crossproduction; }

inline double CSourcePieceWise_TurbSA_Neg::GetProduction(void) { return Production; }

inline double CSourcePieceWise_TurbSA_Neg::GetDestruction(void) { return Destruction; }

inline double CSourcePieceWise_TurbSA_Neg::GetCrossProduction(void) { return CrossProduction; }

inline void CSourcePieceWise_TurbML::SetIntermittency(double intermittency_in) { intermittency = intermittency_in; }

inline void CSourcePieceWise_TurbML::SetProduction(double val_production) { Production = val_production; }

inline void CSourcePieceWise_TurbML::SetDestruction(double val_destruction) { Destruction = val_destruction; }

inline void CSourcePieceWise_TurbML::SetCrossProduction(double val_crossproduction) { CrossProduction = val_crossproduction; }

inline double CSourcePieceWise_TurbML::GetProduction(void) { return Production; }

inline double CSourcePieceWise_TurbML::GetDestruction(void) { return Destruction; }

inline double CSourcePieceWise_TurbML::GetCrossProduction(void) { return CrossProduction; }

inline double CUpwTurkel_Flow::GetPrecond_Beta() { return Beta; }

inline void CNumerics::ComputeResidual(double **val_Jacobian_i, double *val_Jacobian_mui, double ***val_Jacobian_gradi, CConfig *config) { }

inline void CNumerics::ComputeResidual(double **val_Jacobian_i, double *val_Jacobian_mui, double ***val_Jacobian_gradi, 
									double **val_Jacobian_j, double *val_Jacobian_muj, double ***val_Jacobian_gradj, CConfig *config) { }
