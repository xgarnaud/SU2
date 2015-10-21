/*!
 * \file variable_structure.inl
 * \brief In-Line subroutines of the <i>variable_structure.hpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 4.0.1 "Cardinal"
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

inline bool CVariable::SetDensity(void) { return 0; }

inline void CVariable::SetVelSolutionOldDVector(void) { }

inline void CVariable::SetVelSolutionDVector(void) { }

inline void CVariable::SetTraction(unsigned short iVar, unsigned short jVar, su2double val_traction) { }

inline void CVariable::AddTraction(unsigned short iVar, unsigned short jVar, su2double val_traction) { }

inline su2double **CVariable::GetTraction(void) { return NULL; }

inline void CVariable::SetStress(unsigned short iVar, unsigned short jVar, su2double val_stress) { }

inline void CVariable::AddStress(unsigned short iVar, unsigned short jVar, su2double val_stress) { }

inline su2double **CVariable::GetStress(void) { return 0; }

inline void CVariable::SetVonMises_Stress(su2double val_stress) { }

inline su2double CVariable::GetVonMises_Stress(void) { return 0; }

inline void CVariable::SetFlow_Pressure(su2double val_pressure) { }

inline su2double CVariable::GetFlow_Pressure(void) { return 0; }

inline void CVariable::Initialize_Connectivity(void) { }

inline void CVariable::Upgrade_Connectivity(void) { }

inline unsigned short CVariable::Get_Connectivity(void) { return 0; }

inline su2double CVariable::GetBetaInc2(void) { return 0; }

inline su2double CVariable::GetDiffLevelSet(void) { return 0; }

inline su2double CVariable::GetDensityInc(void) { return 0; }

inline su2double CVariable::GetLevelSet(void) { return 0; }

inline su2double CVariable::GetDistance(void) { return 0; }

inline su2double CVariable::GetMassFraction(unsigned short val_Species) { return 0; }

inline void CVariable::SetNon_Physical(bool val_value) { Non_Physical = !val_value; }

inline su2double CVariable::GetNon_Physical(void) { return su2double(Non_Physical); }

inline void CVariable::SetSolution(unsigned short val_var, su2double val_solution) { Solution[val_var] = val_solution; }

inline void CVariable::SetUndivided_Laplacian(unsigned short val_var, su2double val_undivided_laplacian) { Undivided_Laplacian[val_var] = val_undivided_laplacian; }

inline void CVariable::SetAuxVar(su2double val_auxvar) { AuxVar = val_auxvar; }

inline void CVariable::SetSolution_Old(unsigned short val_var, su2double val_solution_old) { Solution_Old[val_var] = val_solution_old; }

inline void CVariable::SetLimiter(unsigned short val_var, su2double val_limiter) { Limiter[val_var] = val_limiter; }

inline void CVariable::SetLimiterPrimitive(unsigned short val_species, unsigned short val_var, su2double val_limiter) { }

inline su2double CVariable::GetLimiterPrimitive(unsigned short val_species, unsigned short val_var) { return 0.0; }

inline void CVariable::SetSolution_Max(unsigned short val_var, su2double val_solution) { Solution_Max[val_var] = val_solution; }

inline void CVariable::SetSolution_Min(unsigned short val_var, su2double val_solution) { Solution_Min[val_var] = val_solution; }

inline void CVariable::SetAuxVarGradient(unsigned short iDim, su2double val_gradient) { Grad_AuxVar[iDim] = val_gradient; }

inline su2double *CVariable::GetSolution(void) { return Solution; }

inline su2double *CVariable::GetSolution_Old(void) { return Solution_Old; }

inline su2double *CVariable::GetSolution_time_n(void) { return Solution_time_n; }

inline su2double *CVariable::GetSolution_time_n1(void) { return Solution_time_n1; }

inline su2double CVariable::GetAuxVar(void) { return AuxVar; }

inline su2double *CVariable::GetUndivided_Laplacian(void) { return Undivided_Laplacian; }

inline su2double CVariable::GetUndivided_Laplacian(unsigned short val_var) { return Undivided_Laplacian[val_var]; }

inline su2double CVariable::GetSolution(unsigned short val_var) { return Solution[val_var]; }

inline su2double CVariable::GetSolution_Old(unsigned short val_var) { return Solution_Old[val_var]; }

inline su2double *CVariable::GetResidual_Sum(void) { return Residual_Sum; }

inline su2double *CVariable::GetResidual_Old(void) { return Residual_Old; }

inline void CVariable::SetGradient(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient[val_var][val_dim] = val_value; }

inline void CVariable::AddGradient(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient[val_var][val_dim] += val_value; }

inline void CVariable::SubtractGradient(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient[val_var][val_dim] -= val_value; }

inline void CVariable::AddAuxVarGradient(unsigned short val_dim, su2double val_value) { Grad_AuxVar[val_dim] += val_value; }

inline void CVariable::SubtractAuxVarGradient(unsigned short val_dim, su2double val_value) { Grad_AuxVar[val_dim] -= val_value; }

inline su2double CVariable::GetGradient(unsigned short val_var, unsigned short val_dim) { return Gradient[val_var][val_dim]; }

inline su2double CVariable::GetLimiter(unsigned short val_var) { return Limiter[val_var]; }

inline su2double CVariable::GetSolution_Max(unsigned short val_var) { return Solution_Max[val_var]; }

inline su2double CVariable::GetSolution_Min(unsigned short val_var) { return Solution_Min[val_var]; }

inline su2double CVariable::GetPreconditioner_Beta() { return 0; }

inline void CVariable::SetPreconditioner_Beta( su2double val_Beta) { }

inline su2double* CVariable::GetWindGust() { return 0; }

inline void CVariable::SetWindGust( su2double* val_WindGust) {}

inline su2double* CVariable::GetWindGustDer() { return 0; }

inline void CVariable::SetWindGustDer( su2double* val_WindGustDer) {}

inline su2double **CVariable::GetGradient(void) { return Gradient; }

inline su2double *CVariable::GetLimiter(void) { return Limiter; }

inline su2double *CVariable::GetAuxVarGradient(void) { return Grad_AuxVar; }

inline su2double CVariable::GetAuxVarGradient(unsigned short val_dim) { return Grad_AuxVar[val_dim]; }

inline su2double *CVariable::GetResTruncError(void) { return Res_TruncError; }

inline void CVariable::SetDelta_Time(su2double val_delta_time) { Delta_Time = val_delta_time; }

inline void CVariable::SetDelta_Time(su2double val_delta_time, unsigned short iSpecies) {  }

inline su2double CVariable::GetDelta_Time(void) { return Delta_Time; }

inline su2double CVariable::GetDelta_Time(unsigned short iSpecies) { return 0;}

inline void CVariable::SetMax_Lambda(su2double val_max_lambda) { Max_Lambda = val_max_lambda; }

inline void CVariable::SetMax_Lambda_Inv(su2double val_max_lambda) { Max_Lambda_Inv = val_max_lambda; }

inline void CVariable::SetMax_Lambda_Inv(su2double val_max_lambda, unsigned short val_species) { }

inline void CVariable::SetMax_Lambda_Visc(su2double val_max_lambda) { Max_Lambda_Visc = val_max_lambda; }

inline void CVariable::SetMax_Lambda_Visc(su2double val_max_lambda, unsigned short val_species) { }

inline void CVariable::SetLambda(su2double val_lambda) { Lambda = val_lambda; }

inline void CVariable::SetLambda(su2double val_lambda, unsigned short iSpecies) {}

inline void CVariable::AddMax_Lambda(su2double val_max_lambda) { Max_Lambda += val_max_lambda; }

inline void CVariable::AddMax_Lambda_Inv(su2double val_max_lambda) { Max_Lambda_Inv += val_max_lambda; }

inline void CVariable::AddMax_Lambda_Visc(su2double val_max_lambda) { Max_Lambda_Visc += val_max_lambda; }

inline void CVariable::AddLambda(su2double val_lambda) { Lambda += val_lambda; }

inline void CVariable::AddLambda(su2double val_lambda, unsigned short iSpecies) {}

inline su2double CVariable::GetMax_Lambda(void) { return Max_Lambda; }

inline su2double CVariable::GetMax_Lambda_Inv(void) { return Max_Lambda_Inv; }

inline su2double CVariable::GetMax_Lambda_Visc(void) { return Max_Lambda_Visc; }

inline su2double CVariable::GetLambda(void) { return Lambda; }

inline su2double CVariable::GetLambda(unsigned short iSpecies) { return 0; }

inline su2double CVariable::GetSensor(void) { return Sensor; }

inline su2double CVariable::GetSensor(unsigned short iSpecies) { return 0;}

inline void CVariable::AddMax_Lambda_Inv(su2double val_max_lambda, unsigned short iSpecies) { }

inline void CVariable::AddMax_Lambda_Visc(su2double val_max_lambda, unsigned short iSpecies) { }

inline void CVariable::SetSensor(su2double val_sensor) { Sensor = val_sensor; }

inline void CVariable::SetSensor(su2double val_sensor, unsigned short val_iSpecies) {}

inline su2double CVariable::GetDensity(void) {	return 0; }

inline su2double CVariable::GetDensity(unsigned short val_iSpecies) {	return 0; }

inline su2double CVariable::GetEnergy(void) { return 0; }

inline su2double *CVariable::GetForceProj_Vector(void) { return NULL; }

inline su2double *CVariable::GetObjFuncSource(void) { return NULL; }

inline su2double *CVariable::GetIntBoundary_Jump(void) { return NULL; }

inline su2double CVariable::GetEddyViscosity(void) { return 0; }

inline su2double CVariable::GetEddyViscosityInc(void) { return 0; }

inline void CVariable::SetGammaEff(void) { }

inline void CVariable::SetGammaSep(su2double gamma_sep) { }

inline su2double CVariable::GetIntermittency(void) { return 0; }

inline su2double CVariable::GetEnthalpy(void) { return 0; }

inline su2double CVariable::GetPressure(void) { return 0; }

inline su2double CVariable::GetPressureInc(void) { return 0; }

inline su2double CVariable::GetDeltaPressure(void) { return 0; }

inline su2double CVariable::GetProjVel(su2double *val_vector) { return 0; }

inline su2double CVariable::GetProjVel(su2double *val_vector, unsigned short val_species) { return 0; }

inline su2double CVariable::GetSoundSpeed(void) { return 0; }

inline su2double CVariable::GetTemperature(void) { return 0; }

inline su2double CVariable::GetTemperature_ve(void) { return 0; }

inline su2double CVariable::GetRhoCv_tr(void) { return 0; }

inline su2double CVariable::GetRhoCv_ve(void) { return 0; }

inline su2double CVariable::GetVelocity(unsigned short val_dim) { return 0; }

inline su2double CVariable::GetVelocity2(void) { return 0; }

inline su2double CVariable::GetVelocity2(unsigned short val_species) { return 0;}

inline su2double CVariable::GetLaminarViscosity(void) { return 0; }

inline su2double CVariable::GetLaminarViscosityInc(void) { return 0; }

inline su2double CVariable::GetLaminarViscosity(unsigned short iSpecies) { return 0; }

inline su2double* CVariable::GetDiffusionCoeff(void) { return NULL; }

inline su2double CVariable::GetThermalConductivity(void) { return 0; }

inline su2double CVariable::GetSpecificHeatCp(void) { return 0; }

inline su2double CVariable::GetThermalConductivity_ve(void) { return 0; }

inline su2double* CVariable::GetVorticity(void) { return 0; }

inline su2double CVariable::GetStrainMag(void) { return 0; }

inline void CVariable::SetForceProj_Vector(su2double *val_ForceProj_Vector) { }

inline void CVariable::SetObjFuncSource(su2double *val_ObjFuncSource) { }

inline void CVariable::SetIntBoundary_Jump(su2double *val_IntBoundary_Jump) { }

inline void CVariable::SetEnthalpy(void) { }

inline bool CVariable::SetPrimVar_Compressible(su2double SharpEdge_Distance, bool check, CConfig *config) { return true; }

inline bool CVariable::SetPrimVar_Incompressible(su2double SharpEdge_Distance, bool check, CConfig *config) { return true; }

inline bool CVariable::SetPrimVar_FreeSurface(su2double SharpEdge_Distance, bool check, CConfig *config) { return true; }

inline bool CVariable::SetPrimVar_Compressible(CConfig *config) { return true; }

inline bool CVariable::SetPrimVar_Compressible(CFluidModel *FluidModel) { return true; }

inline void CVariable::SetSecondaryVar_Compressible(CFluidModel *FluidModel) { }

inline bool CVariable::SetPrimVar_Compressible(su2double eddy_visc, su2double turb_ke, CConfig *config) { return true; }

inline bool CVariable::SetPrimVar_Compressible(su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel) { return true; }

inline bool CVariable::SetPrimVar_Incompressible(su2double Density_Inf, CConfig *config) { return true; }

inline bool CVariable::SetPrimVar_FreeSurface(CConfig *config) { return true; }

inline bool CVariable::SetPrimVar_Incompressible(su2double Density_Inf, su2double Viscosity_Inf, su2double eddy_visc, su2double turb_ke, CConfig *config) { return true; }

inline bool CVariable::SetPrimVar_FreeSurface(su2double eddy_visc, su2double turb_ke, CConfig *config) { return true; }

inline su2double CVariable::GetPrimitive(unsigned short val_var) { return 0; }

inline su2double *CVariable::GetPrimitive(void) { return NULL; }

inline void CVariable::SetPrimitive(unsigned short val_var, su2double val_prim) { }

inline void CVariable::SetPrimitive(su2double *val_prim) { }

inline su2double CVariable::GetSecondary(unsigned short val_var) { return 0; }

inline su2double *CVariable::GetSecondary(void) { return NULL; }

inline void CVariable::SetSecondary(unsigned short val_var, su2double val_secondary) { }

inline void CVariable::SetSecondary(su2double *val_prim) { }

inline bool CVariable::Cons2PrimVar(CConfig *config, su2double *U, su2double *V,
                                    su2double *val_dPdU, su2double *val_dTdU,
                                    su2double *val_dTvedU) { return false; }

inline void CVariable::Prim2ConsVar(CConfig *config, su2double *V, su2double *U) { return; }

inline void CVariable::SetBetaInc2(su2double val_betainc2) { }

inline void CVariable::SetDensityInc(su2double val_density) { }

inline void CVariable::SetPressureInc(void) { }

inline void CVariable::SetVelocityInc(void) { }

inline void CVariable::SetPhi_Old(su2double *val_phi) { }

inline void CVariable::SetDiffLevelSet(su2double val_difflevelset) { }

inline void CVariable::SetdPdrho_e(su2double dPdrho_e) { }

inline void CVariable::SetdPde_rho(su2double dPde_rho) { }

inline void CVariable::SetdTdrho_e(su2double dTdrho_e) { }

inline void CVariable::SetdTde_rho(su2double dTde_rho) { }

inline void CVariable::Setdmudrho_T(su2double dmudrho_T) { }

inline void CVariable::SetdmudT_rho(su2double dmudT_rho) { }

inline void CVariable::Setdktdrho_T(su2double dktdrho_T) { }

inline void CVariable::SetdktdT_rho(su2double dktdT_rho) { }

inline bool CVariable::SetPressure(su2double Gamma) { return false; }

inline bool CVariable::SetPressure(CConfig *config) { return false; }

inline bool CVariable::SetPressure(su2double Gamma, su2double turb_ke) { return false; }

inline void CVariable::SetPressure() { }

inline su2double *CVariable::GetdPdU() { return NULL; }

inline su2double *CVariable::GetdTdU() { return NULL; }

inline su2double *CVariable::GetdTvedU() { return NULL; }

inline su2double CVariable::CalcEve(su2double *V, CConfig *config, unsigned short val_Species) { return 0; }

inline su2double CVariable::CalcHs(su2double *V, CConfig *config, unsigned short val_Species) { return 0; }

inline su2double CVariable::CalcCvve(su2double val_Tve, CConfig *config, unsigned short val_Species) { return 0; }

inline void CVariable::CalcdPdU(su2double *V, CConfig *config, su2double *dPdU) { }

inline void CVariable::CalcdTdU(su2double *V, CConfig *config, su2double *dTdU) { }

inline void CVariable::CalcdTvedU(su2double *V, CConfig *config, su2double *dTvedU) { }

inline void CVariable::SetDeltaPressure(su2double *val_velocity, su2double Gamma) { }

inline bool CVariable::SetSoundSpeed(CConfig *config) { return false; }

inline bool CVariable::SetSoundSpeed() { return false; }

inline bool CVariable::SetSoundSpeed(su2double Gamma) { return false; }

inline bool CVariable::SetTemperature(su2double Gas_Constant) { return false; }

inline bool CVariable::SetTemperature_ve(su2double val_Tve) {return false; }

inline bool CVariable::SetTemperature(CConfig *config) { return false; }

inline void CVariable::SetPrimitive(CConfig *config) { }

inline void CVariable::SetPrimitive(CConfig *config, su2double *Coord) { }

inline void CVariable::SetWallTemperature(su2double Temperature_Wall) { }

inline void CVariable::SetWallTemperature(su2double* Temperature_Wall) { }

inline void CVariable::SetThermalCoeff(CConfig *config) { }

inline void CVariable::SetVelocity(void) { }

inline void CVariable::SetVelocity2(void) { }

inline void CVariable::SetVelocity_Old(su2double *val_velocity) { }

inline void CVariable::SetVelocityInc_Old(su2double *val_velocity) { }

inline void CVariable::SetVel_ResTruncError_Zero(unsigned short iSpecies) { }

inline void CVariable::SetLaminarViscosity(su2double laminarViscosity) { }

inline void CVariable::SetLaminarViscosity(CConfig *config) { }

inline void CVariable::SetLaminarViscosityInc(su2double val_laminar_viscosity_inc) { }

inline void CVariable::SetEddyViscosity(su2double eddy_visc) { }

inline void CVariable::SetEddyViscosityInc(su2double eddy_visc) { }

inline void CVariable::SetThermalConductivity(su2double thermalConductivity) { }

inline void CVariable::SetThermalConductivity(CConfig *config) { }

inline void CVariable::SetSpecificHeatCp(su2double Cp) { }

inline bool CVariable::SetVorticity(bool val_limiter) { return false; }

inline bool CVariable::SetStrainMag(bool val_limiter) { return false; }

inline void CVariable::SetGradient_PrimitiveZero(unsigned short val_primvar) { }

inline void CVariable::AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { }

inline void CVariable::SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { }

inline su2double CVariable::GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) { return 0; }

inline su2double CVariable::GetLimiter_Primitive(unsigned short val_var) { return 0; }

inline void CVariable::SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { }

inline void CVariable::SetLimiter_Primitive(unsigned short val_var, su2double val_value) { }

inline su2double **CVariable::GetGradient_Primitive(void) { return NULL; }

inline su2double *CVariable::GetLimiter_Primitive(void) { return NULL; }

inline void CVariable::SetGradient_SecondaryZero(unsigned short val_secondaryvar) { }

inline void CVariable::AddGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) { }

inline void CVariable::SubtractGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) { }

inline su2double CVariable::GetGradient_Secondary(unsigned short val_var, unsigned short val_dim) { return 0; }

inline su2double CVariable::GetLimiter_Secondary(unsigned short val_var) { return 0; }

inline void CVariable::SetGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) { }

inline void CVariable::SetLimiter_Secondary(unsigned short val_var, su2double val_value) { }

inline su2double **CVariable::GetGradient_Secondary(void) { return NULL; }

inline su2double *CVariable::GetLimiter_Secondary(void) { return NULL; }

inline void CVariable::SetBlendingFunc(su2double val_viscosity, su2double val_dist, su2double val_density) { }

inline su2double CVariable::GetF1blending(void) { return 0; }

inline su2double CVariable::GetF2blending(void) { return 0; }

inline su2double CVariable::GetmuT() { return 0;}

inline void CVariable::SetmuT(su2double val_muT) { }

inline su2double* CVariable::GetSolution_Direct() { return NULL; }

inline void CVariable::SetSolution_Direct(su2double *val_solution_direct) { }

inline void CVariable::SetTimeSpectral_Source(unsigned short val_var, su2double val_source) { }

inline unsigned short CVariable::GetRhosIndex(void) { return 0; }

inline unsigned short CVariable::GetRhoIndex(void) { return 0; }

inline unsigned short CVariable::GetPIndex(void) { return 0; }

inline unsigned short CVariable::GetTIndex(void) { return 0; }

inline unsigned short CVariable::GetTveIndex(void) { return 0; }

inline unsigned short CVariable::GetVelIndex(void) { return 0; }

inline unsigned short CVariable::GetHIndex(void) { return 0; }

inline unsigned short CVariable::GetAIndex(void) { return 0; }

inline unsigned short CVariable::GetRhoCvtrIndex(void) { return 0; }

inline unsigned short CVariable::GetRhoCvveIndex(void) { return 0; }

inline su2double CVariable::GetTimeSpectral_Source(unsigned short val_var) { return 0; }

inline void CVariable::SetEddyViscSens(su2double *val_EddyViscSens, unsigned short numTotalVar) { }

inline su2double *CVariable::GetEddyViscSens(void) { return NULL; }

inline void CVariable::SetSolution_time_n(void) { }

inline void CVariable::SetSolution_Vel(su2double *val_solution_vel) { }

inline void CVariable::SetSolution_Vel(unsigned short val_var, su2double val_solution_vel) { }

inline void CVariable::SetSolution_Vel_time_n(su2double *val_solution_vel_time_n) { }

inline void CVariable::SetSolution_Vel_time_n(void) { }

inline void CVariable::SetSolution_Vel_time_n(unsigned short val_var, su2double val_solution_vel_time_n) { }

inline su2double CVariable::GetSolution_time_n(unsigned short val_var) { return 0; }

inline su2double CVariable::GetSolution_Vel(unsigned short val_var) { return 0; }

inline su2double *CVariable::GetSolution_Vel(void) { return NULL; }

inline su2double CVariable::GetSolution_Vel_time_n(unsigned short val_var) { return 0; }

inline su2double *CVariable::GetSolution_Vel_time_n(void) { return NULL; }

inline void CVariable::SetSolution_Accel(su2double *val_solution_accel) { }

inline void CVariable::SetSolution_Accel(unsigned short val_var, su2double val_solution_accel) { }

inline void CVariable::SetSolution_Accel_time_n(su2double *val_solution_accel_time_n) { }

inline void CVariable::SetSolution_Accel_time_n(void) { }

inline void CVariable::SetSolution_Accel_time_n(unsigned short val_var, su2double val_solution_accel_time_n) { }

inline su2double CVariable::GetSolution_Accel(unsigned short val_var) { return 0; }

inline su2double *CVariable::GetSolution_Accel(void) { return NULL; }

inline su2double CVariable::GetSolution_Accel_time_n(unsigned short val_var) { return 0; }

inline su2double *CVariable::GetSolution_Accel_time_n(void) { return NULL; }

inline void CVariable::SetSolution_Pred(su2double *val_solution_pred){  }

inline void CVariable::SetSolution_Pred(void){ }

inline su2double CVariable::GetSolution_Pred(unsigned short val_var){ return 0.0; }

inline su2double *CVariable::GetSolution_Pred(void){ return NULL; }

inline void CVariable::SetSolution_Pred_Old(su2double *val_solution_pred_Old){  }

inline void CVariable::SetSolution_Pred_Old(void){ }

inline su2double CVariable::GetSolution_Pred_Old(unsigned short val_var){ return 0.0; }

inline su2double *CVariable::GetSolution_Pred_Old(void){ return NULL; }

inline su2double CEulerVariable::GetDensity(void) { return Solution[0]; }

inline su2double CEulerVariable::GetDensityInc(void) { return Primitive[nDim+1]; }

inline su2double CEulerVariable::GetLevelSet(void) { return Primitive[nDim+5]; }

inline su2double CEulerVariable::GetDistance(void) { return Primitive[nDim+6]; }

inline su2double CEulerVariable::GetBetaInc2(void) { return Primitive[nDim+2]; }

inline su2double CEulerVariable::GetEnergy(void) { return Solution[nVar-1]/Solution[0]; };

inline su2double CEulerVariable::GetEnthalpy(void) { return Primitive[nDim+3]; }

inline su2double CEulerVariable::GetPressure(void) { return Primitive[nDim+1]; }

inline su2double CEulerVariable::GetPressureInc(void) { return Primitive[0]; }

inline su2double CEulerVariable::GetSoundSpeed(void) { return Primitive[nDim+4]; }

inline su2double CEulerVariable::GetTemperature(void) { return Primitive[0]; }

inline su2double CEulerVariable::GetVelocity(unsigned short val_dim) { return Primitive[val_dim+1]; }

inline su2double CEulerVariable::GetVelocity2(void) { return Velocity2; }

inline bool CEulerVariable::SetDensity(void) {
  Primitive[nDim+2] = Solution[0];
  if (Primitive[nDim+2] > 0.0) return false;
  else return true;
}

inline void CEulerVariable::SetDensityInc(su2double val_density) { Primitive[nDim+1] = val_density; }

inline bool CEulerVariable::SetPressure(su2double pressure) {
  Primitive[nDim+1] = pressure;
  if (Primitive[nDim+1] > 0.0) return false;
  else return true;
}

inline void CEulerVariable::SetPressureInc(void) { Primitive[0] = Solution[0]; }

inline void CEulerVariable::SetVelocity(void) {
  Velocity2 = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    Primitive[iDim+1] = Solution[iDim+1] / Solution[0];
    Velocity2 += Primitive[iDim+1]*Primitive[iDim+1];
  }
}

inline void CEulerVariable::SetVelocityInc(void) {
  Velocity2 = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    Primitive[iDim+1] = Solution[iDim+1] / Primitive[nDim+1];
    Velocity2 += Primitive[iDim+1]*Primitive[iDim+1];
  }
}

inline void CEulerVariable::SetEnthalpy(void) { Primitive[nDim+3] = (Solution[nVar-1] + Primitive[nDim+1]) / Solution[0]; }

inline void CEulerVariable::SetBetaInc2(su2double val_betainc2) { Primitive[nDim+2] = val_betainc2; }

inline bool CEulerVariable::SetSoundSpeed(su2double soundspeed2) {
  su2double radical = soundspeed2;
  if (radical < 0.0) return true;
  else {
    Primitive[nDim+4] = sqrt(radical);
    return false;
  }
}

inline bool CEulerVariable::SetTemperature(su2double temperature) {
  Primitive[0] = temperature;
  if (Primitive[0] > 0.0) return false;
  else return true;
}

inline void CEulerVariable::SetdPdrho_e(su2double dPdrho_e) {
  Secondary[0] = dPdrho_e;
}

inline void CEulerVariable::SetdPde_rho(su2double dPde_rho) {
  Secondary[1] = dPde_rho;
}

inline su2double CEulerVariable::GetPrimitive(unsigned short val_var) { return Primitive[val_var]; }

inline void CEulerVariable::SetPrimitive(unsigned short val_var, su2double val_prim) { Primitive[val_var] = val_prim; }

inline void CEulerVariable::SetPrimitive(su2double *val_prim) {
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++)
    Primitive[iVar] = val_prim[iVar];
}

inline su2double *CEulerVariable::GetPrimitive(void) { return Primitive; }

inline su2double CEulerVariable::GetSecondary(unsigned short val_var) { return Secondary[val_var]; }

inline void CEulerVariable::SetSecondary(unsigned short val_var, su2double val_secondary) { Secondary[val_var] = val_secondary; }

inline void CEulerVariable::SetSecondary(su2double *val_secondary) {
  for (unsigned short iVar = 0; iVar < nSecondaryVar; iVar++)
    Secondary[iVar] = val_secondary[iVar];
}

inline su2double *CEulerVariable::GetSecondary(void) { return Secondary; }

inline void CEulerVariable::SetVelocity_Old(su2double *val_velocity) {
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Solution_Old[iDim+1] = val_velocity[iDim]*Solution[0];
}

inline void CEulerVariable::SetVelocityInc_Old(su2double *val_velocity) {
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Solution_Old[iDim+1] = val_velocity[iDim]*Primitive[nDim+1];
}

inline void CEulerVariable::AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Primitive[val_var][val_dim] += val_value; }

inline void CEulerVariable::SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Primitive[val_var][val_dim] -= val_value; }

inline su2double CEulerVariable::GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) { return Gradient_Primitive[val_var][val_dim]; }

inline su2double CEulerVariable::GetLimiter_Primitive(unsigned short val_var) { return Limiter_Primitive[val_var]; }

inline void CEulerVariable::SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Primitive[val_var][val_dim] = val_value; }

inline void CEulerVariable::SetLimiter_Primitive(unsigned short val_var, su2double val_value) { Limiter_Primitive[val_var] = val_value; }

inline su2double **CEulerVariable::GetGradient_Primitive(void) { return Gradient_Primitive; }

inline su2double *CEulerVariable::GetLimiter_Primitive(void) { return Limiter_Primitive; }

inline void CEulerVariable::AddGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Secondary[val_var][val_dim] += val_value; }

inline void CEulerVariable::SubtractGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Secondary[val_var][val_dim] -= val_value; }

inline su2double CEulerVariable::GetGradient_Secondary(unsigned short val_var, unsigned short val_dim) { return Gradient_Secondary[val_var][val_dim]; }

inline su2double CEulerVariable::GetLimiter_Secondary(unsigned short val_var) { return Limiter_Secondary[val_var]; }

inline void CEulerVariable::SetGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Secondary[val_var][val_dim] = val_value; }

inline void CEulerVariable::SetLimiter_Secondary(unsigned short val_var, su2double val_value) { Limiter_Secondary[val_var] = val_value; }

inline su2double **CEulerVariable::GetGradient_Secondary(void) { return Gradient_Secondary; }

inline su2double *CEulerVariable::GetLimiter_Secondary(void) { return Limiter_Secondary; }

inline void CEulerVariable::SetTimeSpectral_Source(unsigned short val_var, su2double val_source) { TS_Source[val_var] = val_source; }

inline su2double CEulerVariable::GetTimeSpectral_Source(unsigned short val_var) { return TS_Source[val_var]; }

inline su2double CEulerVariable::GetPreconditioner_Beta() { return Precond_Beta; }

inline void CEulerVariable::SetPreconditioner_Beta(su2double val_Beta) { Precond_Beta = val_Beta; }

inline void CEulerVariable::SetWindGust( su2double* val_WindGust) {
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    WindGust[iDim] = val_WindGust[iDim];}

inline su2double* CEulerVariable::GetWindGust() { return WindGust;}

inline void CEulerVariable::SetWindGustDer( su2double* val_WindGustDer) {
  for (unsigned short iDim = 0; iDim < nDim+1; iDim++)
    WindGustDer[iDim] = val_WindGustDer[iDim];}

inline su2double* CEulerVariable::GetWindGustDer() { return WindGustDer;}

inline su2double CNSVariable::GetEddyViscosity(void) { return Primitive[nDim+6]; }

inline su2double CNSVariable::GetEddyViscosityInc(void) { return Primitive[nDim+4]; }

inline su2double CNSVariable::GetLaminarViscosity(void) { return Primitive[nDim+5]; }

inline su2double CNSVariable::GetLaminarViscosityInc(void) { return Primitive[nDim+3]; }

inline su2double CNSVariable::GetThermalConductivity(void) { return Primitive[nDim+7]; }

inline su2double CNSVariable::GetSpecificHeatCp(void) { return Primitive[nDim+8]; }

inline su2double* CNSVariable::GetVorticity(void) { return Vorticity; }

inline su2double CNSVariable::GetStrainMag(void) { return StrainMag; }

inline void CNSVariable::SetLaminarViscosity(su2double laminarViscosity) {
  Primitive[nDim+5] = laminarViscosity;
}

inline void CNSVariable::SetThermalConductivity(su2double thermalConductivity) {
  Primitive[nDim+7] = thermalConductivity;
}

inline void CNSVariable::SetSpecificHeatCp(su2double Cp) {
  Primitive[nDim+8] = Cp;
}

inline void CNSVariable::SetdTdrho_e(su2double dTdrho_e) {
  Secondary[2] = dTdrho_e;
}

inline void CNSVariable::SetdTde_rho(su2double dTde_rho) {
  Secondary[3] = dTde_rho;
}

inline void CNSVariable::Setdmudrho_T(su2double dmudrho_T) {
  Secondary[4] = dmudrho_T;
}

inline void CNSVariable::SetdmudT_rho(su2double dmudT_rho) {
  Secondary[5] = dmudT_rho;
}

inline void CNSVariable::Setdktdrho_T(su2double dktdrho_T) {
  Secondary[6] = dktdrho_T;
}

inline void CNSVariable::SetdktdT_rho(su2double dktdT_rho) {
  Secondary[7] = dktdT_rho;
}

inline void CNSVariable::SetLaminarViscosityInc(su2double val_laminar_viscosity_inc) { Primitive[nDim+3] = val_laminar_viscosity_inc; }

inline void CNSVariable::SetEddyViscosity(su2double eddy_visc) { Primitive[nDim+6] = eddy_visc; }

inline void CNSVariable::SetEddyViscosityInc(su2double eddy_visc) { Primitive[nDim+4] = eddy_visc; }

inline void CNSVariable::SetWallTemperature(su2double Temperature_Wall ) { Primitive[0] = Temperature_Wall; }

inline su2double CTransLMVariable::GetIntermittency() { return Solution[0]; }

inline void CTransLMVariable::SetGammaSep(su2double gamma_sep_in) {gamma_sep = gamma_sep_in;}

inline su2double *CAdjEulerVariable::GetForceProj_Vector(void) { return ForceProj_Vector; }

inline su2double *CAdjEulerVariable::GetObjFuncSource(void) { return ObjFuncSource; }

inline su2double *CAdjEulerVariable::GetIntBoundary_Jump(void) { return IntBoundary_Jump; }

inline void CAdjEulerVariable::SetForceProj_Vector(su2double *val_ForceProj_Vector) { for (unsigned short iDim = 0; iDim < nDim; iDim++) ForceProj_Vector[iDim] = val_ForceProj_Vector[iDim]; }

inline void CAdjEulerVariable::SetObjFuncSource(su2double *val_ObjFuncSource) { for (unsigned short iVar = 0; iVar < nVar; iVar++) ObjFuncSource[iVar] = val_ObjFuncSource[iVar]; }

inline void CAdjEulerVariable::SetIntBoundary_Jump(su2double *val_IntBoundary_Jump) { for (unsigned short iVar = 0; iVar < nVar; iVar++) IntBoundary_Jump[iVar] = val_IntBoundary_Jump[iVar]; }

inline void CAdjEulerVariable::SetPhi_Old(su2double *val_phi) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1]=val_phi[iDim]; };

inline void CAdjEulerVariable::SetTimeSpectral_Source(unsigned short val_var, su2double val_source) { TS_Source[val_var] = val_source; }

inline su2double CAdjEulerVariable::GetTimeSpectral_Source(unsigned short val_var) { return TS_Source[val_var]; }

inline su2double *CAdjNSVariable::GetForceProj_Vector(void) { return ForceProj_Vector; }

inline void CAdjNSVariable::SetForceProj_Vector(su2double *val_ForceProj_Vector) {	for (unsigned short iDim = 0; iDim < nDim; iDim++) ForceProj_Vector[iDim] = val_ForceProj_Vector[iDim]; }

inline void CAdjNSVariable::SetPhi_Old(su2double *val_phi) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1] = val_phi[iDim]; };

inline void CAdjNSVariable::SetVelSolutionOldDVector(void) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1] = ForceProj_Vector[iDim]; };

inline void CAdjNSVariable::SetVelSolutionDVector(void) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution[iDim+1] = ForceProj_Vector[iDim]; };

inline su2double *CLinEulerVariable::GetForceProj_Vector(void) { return ForceProj_Vector; }

inline void CLinEulerVariable::SetForceProj_Vector(su2double *val_ForceProj_Vector) { for (unsigned short iDim = 0; iDim < nDim; iDim++) ForceProj_Vector[iDim] = val_ForceProj_Vector[iDim]; }

inline void CLinEulerVariable::SetDeltaVel_Old(su2double *val_deltavel) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1]=val_deltavel[iDim]; };

inline su2double CLinEulerVariable::GetDeltaPressure(void) { return DeltaPressure; }

inline void CFEAVariable::SetStress(unsigned short iVar, unsigned short jVar, su2double val_stress) { Stress[iVar][jVar] = val_stress; }

inline void CFEAVariable::AddStress(unsigned short iVar, unsigned short jVar, su2double val_stress) { Stress[iVar][jVar] += val_stress; }

inline su2double **CFEAVariable::GetStress(void) { return Stress; }

inline void CFEAVariable::SetVonMises_Stress(su2double val_stress) { VonMises_Stress = val_stress; }

inline su2double CFEAVariable::GetVonMises_Stress(void) { return VonMises_Stress; }

inline void CFEAVariable::SetFlow_Pressure(su2double val_pressure) { Flow_Pressure = val_pressure; }

inline su2double CFEAVariable::GetFlow_Pressure(void) { return Flow_Pressure; }

inline void CFEAVariable::Initialize_Connectivity(void) { nAttachedElements = 0; }

inline void CFEAVariable::Upgrade_Connectivity(void) { nAttachedElements += 1; }

inline unsigned short CFEAVariable::Get_Connectivity(void) { return nAttachedElements; }

inline void CFEABoundVariable::SetTraction(unsigned short iVar, unsigned short jVar, su2double val_traction) { Traction[iVar][jVar] = val_traction; }

inline void CFEABoundVariable::AddTraction(unsigned short iVar, unsigned short jVar, su2double val_traction) { Traction[iVar][jVar] += val_traction; }

inline su2double **CFEABoundVariable::GetTraction(void) { return Traction; }

inline void CFEAVariable::SetSolution_Vel(unsigned short val_var, su2double val_solution_vel) { Solution_Vel[val_var] = val_solution_vel; }

inline void CFEAVariable::SetSolution_Vel_time_n(unsigned short val_var, su2double val_solution_vel_time_n) { Solution_Vel_time_n[val_var] = val_solution_vel_time_n; }

inline su2double CFEAVariable::GetSolution_time_n(unsigned short val_var) { return Solution_time_n[val_var]; }

inline su2double CFEAVariable::GetSolution_Vel(unsigned short val_var) { return Solution_Vel[val_var]; }

inline su2double *CFEAVariable::GetSolution_Vel(void) { return Solution_Vel; }

inline su2double CFEAVariable::GetSolution_Vel_time_n(unsigned short val_var) { return Solution_Vel_time_n[val_var]; }

inline su2double *CFEAVariable::GetSolution_Vel_time_n(void) { return Solution_Vel_time_n; }

inline void CFEAVariable::SetSolution_Accel(unsigned short val_var, su2double val_solution_accel) { Solution_Accel[val_var] = val_solution_accel;  }

inline void CFEAVariable::SetSolution_Accel_time_n(unsigned short val_var, su2double val_solution_accel_time_n) { Solution_Accel_time_n[val_var] = val_solution_accel_time_n; }

inline su2double CFEAVariable::GetSolution_Accel(unsigned short val_var) { return Solution_Accel[val_var]; }

inline su2double *CFEAVariable::GetSolution_Accel(void) { return Solution_Accel; }

inline su2double CFEAVariable::GetSolution_Accel_time_n(unsigned short val_var) { return Solution_Accel_time_n[val_var]; }

inline su2double *CFEAVariable::GetSolution_Accel_time_n(void) { return Solution_Accel_time_n; }

inline void CFEAVariable::SetSolution_Pred(su2double *val_solution_pred){ Solution_Pred = val_solution_pred;  }

inline su2double CFEAVariable::GetSolution_Pred(unsigned short val_var){ return Solution_Pred[val_var]; }

inline su2double *CFEAVariable::GetSolution_Pred(void){ return Solution_Pred; }

inline void CFEAVariable::SetSolution_Pred_Old(su2double *val_solution_pred_Old){ Solution_Pred_Old = val_solution_pred_Old;  }

inline su2double CFEAVariable::GetSolution_Pred_Old(unsigned short val_var){ return Solution_Pred_Old[val_var]; }

inline su2double *CFEAVariable::GetSolution_Pred_Old(void){ return Solution_Pred_Old; }

inline su2double* CWaveVariable::GetSolution_Direct() { return Solution_Direct;}

inline void CWaveVariable::SetSolution_Direct(su2double *val_solution_direct) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Direct[iVar] += val_solution_direct[iVar];}

inline su2double* CPotentialVariable::GetChargeDensity() { return Charge_Density;}

inline void CPotentialVariable::SetChargeDensity(su2double positive_charge, su2double negative_charge) {Charge_Density[0] = positive_charge; Charge_Density[1] = negative_charge;}

inline su2double* CHeatVariable::GetSolution_Direct() { return Solution_Direct;}

inline void CHeatVariable::SetSolution_Direct(su2double *val_solution_direct) { for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Direct[iVar] += val_solution_direct[iVar];}

inline void CTurbSAVariable::SetTimeSpectral_Source(unsigned short val_var, su2double val_source) { TS_Source[val_var] = val_source; }

inline su2double CTurbSAVariable::GetTimeSpectral_Source(unsigned short val_var) { return TS_Source[val_var]; }


inline void CTurbMLVariable::SetTimeSpectral_Source(unsigned short val_var, su2double val_source) { TS_Source[val_var] = val_source; }

inline su2double CTurbMLVariable::GetTimeSpectral_Source(unsigned short val_var) { return TS_Source[val_var]; }

inline su2double CTurbSSTVariable::GetF1blending() { return F1; }

inline su2double CTurbSSTVariable::GetF2blending() { return F2; }

inline su2double CTurbSSTVariable::GetCrossDiff() { return CDkw; }

inline void CAdjTurbVariable::SetEddyViscSens(su2double *val_EddyViscSens, unsigned short numTotalVar) {
  for (unsigned short iVar = 0; iVar < numTotalVar; iVar++) {
    EddyViscSens[iVar] = val_EddyViscSens[iVar];}
}

inline su2double *CAdjTurbVariable::GetEddyViscSens(void) { return EddyViscSens; }

inline su2double CTNE2EulerVariable::GetDensity(void) { return Primitive[RHO_INDEX]; }

inline su2double CTNE2EulerVariable::GetMassFraction(unsigned short val_Species) {
  return Primitive[RHOS_INDEX+val_Species] / Primitive[RHO_INDEX];
}

inline su2double CTNE2EulerVariable::GetEnergy(void) { return Solution[nSpecies+nDim]/Primitive[RHO_INDEX]; };

inline su2double CTNE2EulerVariable::GetEnthalpy(void) { return Primitive[H_INDEX]; }

inline su2double CTNE2EulerVariable::GetPressure(void) { return Primitive[P_INDEX]; }

inline su2double CTNE2EulerVariable::GetSoundSpeed(void) { return Primitive[A_INDEX]; }

inline su2double CTNE2EulerVariable::GetTemperature(void) { return Primitive[T_INDEX]; }

inline bool CTNE2EulerVariable::SetTemperature(su2double val_T) { Primitive[T_INDEX] = val_T; return false; }

inline su2double CTNE2EulerVariable::GetTemperature_ve(void) { return Primitive[TVE_INDEX]; }

inline bool CTNE2EulerVariable::SetTemperature_ve(su2double val_Tve) { Primitive[TVE_INDEX] = val_Tve; return false; }

inline su2double CTNE2EulerVariable::GetRhoCv_tr(void) { return Primitive[RHOCVTR_INDEX]; }

inline su2double CTNE2EulerVariable::GetRhoCv_ve(void) { return Primitive[RHOCVVE_INDEX]; }

inline su2double* CTNE2EulerVariable::GetdPdU(void) { return dPdU; }

inline su2double* CTNE2EulerVariable::GetdTdU(void) { return dTdU; }

inline su2double* CTNE2EulerVariable::GetdTvedU(void) { return dTvedU; }

inline su2double CTNE2EulerVariable::GetVelocity(unsigned short val_dim) { return Primitive[VEL_INDEX+val_dim]; }

inline su2double CTNE2EulerVariable::GetVelocity2(void) { return Velocity2; }

inline void CTNE2EulerVariable::SetEnthalpy(void) { Primitive[H_INDEX] = (Solution[nSpecies+nDim] + Primitive[P_INDEX]) / Primitive[RHO_INDEX]; }

inline su2double CTNE2EulerVariable::GetPrimitive(unsigned short val_var) { return Primitive[val_var]; }

inline void CTNE2EulerVariable::SetPrimitive(unsigned short val_var, su2double val_prim) { Primitive[val_var] = val_prim; }

inline void CTNE2EulerVariable::SetPrimitive(su2double *val_prim) {
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++)
    Primitive[iVar] = val_prim[iVar];
}

inline su2double *CTNE2EulerVariable::GetPrimitive(void) { return Primitive; }

inline void CTNE2EulerVariable::SetVelocity_Old(su2double *val_velocity) {
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Solution_Old[nSpecies+iDim] = val_velocity[iDim]*Primitive[RHO_INDEX];
}

inline void CTNE2EulerVariable::AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Primitive[val_var][val_dim] += val_value; }

inline void CTNE2EulerVariable::SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Primitive[val_var][val_dim] -= val_value; }

inline su2double CTNE2EulerVariable::GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) { return Gradient_Primitive[val_var][val_dim]; }

inline void CTNE2EulerVariable::SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Primitive[val_var][val_dim] = val_value; }

inline su2double **CTNE2EulerVariable::GetGradient_Primitive(void) { return Gradient_Primitive; }

inline su2double *CTNE2EulerVariable::GetLimiter_Primitive(void) { return Limiter_Primitive; }

inline su2double CTNE2EulerVariable::GetLimiter_Primitive(unsigned short val_var) { return Limiter_Primitive[val_var]; }

inline void CTNE2EulerVariable::SetLimiter_Primitive(unsigned short val_var, su2double val_value) { Limiter_Primitive[val_var] = val_value; }

inline void CTNE2EulerVariable::SetLimiter(unsigned short val_var, su2double val_value) { Limiter[val_var] = val_value; }

inline su2double CTNE2EulerVariable::GetPreconditioner_Beta() { return Precond_Beta; }

inline void CTNE2EulerVariable::SetPreconditioner_Beta(su2double val_Beta) { Precond_Beta = val_Beta; }

inline unsigned short CTNE2EulerVariable::GetRhosIndex(void) { return RHOS_INDEX; }

inline unsigned short CTNE2EulerVariable::GetRhoIndex(void) { return RHO_INDEX; }

inline unsigned short CTNE2EulerVariable::GetPIndex(void) { return P_INDEX; }

inline unsigned short CTNE2EulerVariable::GetTIndex(void) { return T_INDEX; }

inline unsigned short CTNE2EulerVariable::GetTveIndex(void) { return TVE_INDEX; }

inline unsigned short CTNE2EulerVariable::GetVelIndex(void) { return VEL_INDEX; }

inline unsigned short CTNE2EulerVariable::GetHIndex(void) { return H_INDEX; }

inline unsigned short CTNE2EulerVariable::GetAIndex(void) { return A_INDEX; }

inline unsigned short CTNE2EulerVariable::GetRhoCvtrIndex(void) { return RHOCVTR_INDEX; }

inline unsigned short CTNE2EulerVariable::GetRhoCvveIndex(void) { return RHOCVVE_INDEX; }

inline su2double* CTNE2NSVariable::GetDiffusionCoeff(void) { return DiffusionCoeff; }

inline su2double  CTNE2NSVariable::GetLaminarViscosity(void) { return LaminarViscosity; }

inline su2double  CTNE2NSVariable::GetThermalConductivity(void) { return ThermalCond; }

inline su2double  CTNE2NSVariable::GetThermalConductivity_ve(void) { return ThermalCond_ve; }

inline su2double*  CTNE2NSVariable::GetVorticity(void) { return Vorticity; }

inline void    CTNE2NSVariable::SetWallTemperature(su2double Temperature_Wall ) { Primitive[T_INDEX] = Temperature_Wall; }

inline void    CAdjTNE2EulerVariable::SetPhi_Old(su2double *val_phi) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1]=val_phi[iDim]; };

inline su2double *CAdjTNE2EulerVariable::GetObjFuncSource(void) { return ObjFuncSource; }

inline void    CAdjTNE2EulerVariable::SetObjFuncSource(su2double *val_ObjFuncSource) { for (unsigned short iVar = 0; iVar < nVar; iVar++) ObjFuncSource[iVar] = val_ObjFuncSource[iVar]; }

inline su2double *CAdjTNE2EulerVariable::GetForceProj_Vector(void) { return ForceProj_Vector; }

inline void    CAdjTNE2EulerVariable::SetForceProj_Vector(su2double *val_ForceProj_Vector) { for (unsigned short iDim = 0; iDim < nDim; iDim++) ForceProj_Vector[iDim] = val_ForceProj_Vector[iDim]; }

inline su2double  CAdjTNE2EulerVariable::GetTheta(void) { return Theta; }

inline void    CAdjTNE2NSVariable::SetPhi_Old(su2double *val_phi) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1] = val_phi[iDim]; };

inline su2double* CAdjTNE2NSVariable::GetForceProj_Vector(void) { return ForceProj_Vector; }

inline void    CAdjTNE2NSVariable::SetForceProj_Vector(su2double *val_ForceProj_Vector) {	for (unsigned short iDim = 0; iDim < nDim; iDim++) ForceProj_Vector[iDim] = val_ForceProj_Vector[iDim]; }

inline void    CAdjTNE2NSVariable::SetVelSolutionDVector(void) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution[iDim+1] = ForceProj_Vector[iDim]; };

inline void    CAdjTNE2NSVariable::SetVelSolutionOldDVector(void) { for (unsigned short iDim = 0; iDim < nDim; iDim++) Solution_Old[iDim+1] = ForceProj_Vector[iDim]; };

inline su2double  CAdjTNE2NSVariable::GetTheta(void) { return Theta; }

inline void CVariable::RegisterSolution(bool input) {
  if (input) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      AD::RegisterInput(Solution[iVar]);
  }
  else { for (unsigned short iVar = 0; iVar < nVar; iVar++)
      AD::RegisterOutput(Solution[iVar]);}
}

inline void CVariable::RegisterSolution_time_n(){
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    AD::RegisterInput(Solution_time_n[iVar]);
}

inline void CVariable::RegisterSolution_time_n1(){
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    AD::RegisterInput(Solution_time_n1[iVar]);
}

inline void CVariable::SetAdjointSolution(su2double *adj_sol){
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
        SU2_TYPE::SetDerivative(Solution[iVar], SU2_TYPE::GetPrimary(adj_sol[iVar]));
}


inline void CVariable::GetAdjointSolution(su2double *adj_sol){
    for (unsigned short iVar = 0; iVar < nVar; iVar++){
        adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution[iVar]);
    }
}

inline void CVariable::SetAdjointSolution_time_n(su2double *adj_sol){
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
      SU2_TYPE::SetDerivative(Solution_time_n[iVar], SU2_TYPE::GetPrimary(adj_sol[iVar]));
}


inline void CVariable::GetAdjointSolution_time_n(su2double *adj_sol){
  for (unsigned short iVar = 0; iVar < nVar; iVar++){
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_time_n[iVar]);
  }
}

inline void CVariable::SetAdjointSolution_time_n1(su2double *adj_sol){
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
      SU2_TYPE::SetDerivative(Solution_time_n1[iVar], SU2_TYPE::GetPrimary(adj_sol[iVar]));
}


inline void CVariable::GetAdjointSolution_time_n1(su2double *adj_sol){
  for (unsigned short iVar = 0; iVar < nVar; iVar++){
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_time_n1[iVar]);
  }
}
inline void CVariable::SetDual_Time_Derivative(unsigned short iVar, su2double der){}

inline void CDiscAdjVariable::SetDual_Time_Derivative(unsigned short iVar, su2double der){DualTime_Derivative[iVar] = der;}

inline void CVariable::SetDual_Time_Derivative_n(unsigned short iVar, su2double der){}

inline void CDiscAdjVariable::SetDual_Time_Derivative_n(unsigned short iVar, su2double der){DualTime_Derivative_n[iVar] = der;}

inline su2double CVariable::GetDual_Time_Derivative(unsigned short iVar){return 0.0;}

inline su2double CDiscAdjVariable::GetDual_Time_Derivative(unsigned short iVar){return DualTime_Derivative[iVar];}

inline su2double CVariable::GetDual_Time_Derivative_n(unsigned short iVar){return 0.0;}

inline su2double CDiscAdjVariable::GetDual_Time_Derivative_n(unsigned short iVar){return DualTime_Derivative_n[iVar];}

inline void CVariable::SetSensitivity(unsigned short iDim, su2double val){}

inline su2double CVariable::GetSensitivity(unsigned short iDim){ return 0.0; }

inline void CDiscAdjVariable::SetSensitivity(unsigned short iDim, su2double val){Sensitivity[iDim] = val;}

inline su2double CDiscAdjVariable::GetSensitivity(unsigned short iDim){return Sensitivity[iDim];}
