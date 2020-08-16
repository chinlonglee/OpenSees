/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// ReSharper disable CppClangTidyClangDiagnosticSignConversion
#include "LeeSparse.h"
#include <ID.h>
#include <iostream>
#include <Matrix.h>
#include <slu_ddefs.h>
#include <sp_sparse/csc_form.hpp>

int LeeSparse::setSize(Graph& theGraph) {
	const auto flag = SparseGenColLinSOE::setSize(theGraph);

	if(0 != flag) return flag;

	current_stiffness = triplet_form<double>(size, size, 1000);
	stiffness = triplet_form<double>(size, size, 1000);
	damping = triplet_form<double>(size, size, 1000);
	mass = triplet_form<double>(size, size, 1000);

	return 0;
}

int LeeSparse::solve() {
	SuperMatrix A, L, U, B;
	int *perm_r, *perm_c;
	int info;
	superlu_options_t options;
	SuperLUStat_t stat;

	const csc_form<double> csc_mat(global_stiffness);

	increment = residual;

	dCreate_CompCol_Matrix(&A, csc_mat.n_rows, csc_mat.n_cols, csc_mat.c_size, csc_mat.val_idx, csc_mat.row_idx, csc_mat.col_ptr, SLU_NC, SLU_D, SLU_GE);

	dCreate_Dense_Matrix(&B, increment.Size(), 1, &increment(0), increment.Size(), SLU_DN, SLU_D, SLU_GE);

	if(!((perm_r = intMalloc(csc_mat.n_rows))) || !((perm_c = intMalloc(csc_mat.n_cols)))) return -1;

	set_default_options(&options);
	options.ColPerm = NATURAL;
	StatInit(&stat);
	dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

	SUPERLU_FREE(perm_r);
	SUPERLU_FREE(perm_c);
	// Destroy_CompCol_Matrix(&A);
	// Destroy_SuperMatrix_Store(&B);
	Destroy_SuperNode_Matrix(&L);
	Destroy_CompCol_Matrix(&U);
	StatFree(&stat);

	for(auto I = 0; I < size; ++I) *(X + I) = increment(I);

	if(0 != info) throw;

	return 0;
}

int LeeSparse::add_current_stiffness(const Matrix& EK, const ID& EI, const double S) {
	if(1. == S) { for(auto I = 0; I < EI.Size(); ++I) if(EI(I) >= 0) for(auto J = 0; J < EI.Size(); ++J) if(EI(J) >= 0) current_stiffness.at(EI(I), EI(J)) = EK(I, J); } else for(auto I = 0; I < EI.Size(); ++I) if(EI(I) >= 0) for(auto J = 0; J < EI.Size(); ++J) if(EI(J) >= 0) current_stiffness.at(EI(I), EI(J)) = S * EK(I, J);

	return 0;
}

int LeeSparse::add_stiffness(const Matrix& EK, const ID& EI, const double S) {
	if(1. == S) { for(auto I = 0; I < EI.Size(); ++I) if(EI(I) >= 0) for(auto J = 0; J < EI.Size(); ++J) if(EI(J) >= 0) stiffness.at(EI(I), EI(J)) = EK(I, J); } else for(auto I = 0; I < EI.Size(); ++I) if(EI(I) >= 0) for(auto J = 0; J < EI.Size(); ++J) if(EI(J) >= 0) stiffness.at(EI(I), EI(J)) = S * EK(I, J);

	return 0;
}

int LeeSparse::add_damping(const Matrix& EC, const ID& EI, const double S) {
	if(1. == S) { for(auto I = 0; I < EI.Size(); ++I) if(EI(I) >= 0) for(auto J = 0; J < EI.Size(); ++J) if(EI(J) >= 0) damping.at(EI(I), EI(J)) = EC(I, J); } else for(auto I = 0; I < EI.Size(); ++I) if(EI(I) >= 0) for(auto J = 0; J < EI.Size(); ++J) if(EI(J) >= 0) damping.at(EI(I), EI(J)) = S * EC(I, J);

	return 0;
}

int LeeSparse::add_mass(const Matrix& EM, const ID& EI, const double S) {
	if(1. == S) { for(auto I = 0; I < EI.Size(); ++I) if(EI(I) >= 0) for(auto J = 0; J < EI.Size(); ++J) if(EI(J) >= 0) mass.at(EI(I), EI(J)) = EM(I, J); } else for(auto I = 0; I < EI.Size(); ++I) if(EI(I) >= 0) for(auto J = 0; J < EI.Size(); ++J) if(EI(J) >= 0) mass.at(EI(I), EI(J)) = S * EM(I, J);

	return 0;
}

void LeeSparse::zero_current_stiffness() const { current_stiffness.zeros(); }

void LeeSparse::zero_stiffness() const { stiffness.zeros(); }

void LeeSparse::zero_damping() const { damping.zeros(); }

void LeeSparse::zero_mass() const { mass.zeros(); }
