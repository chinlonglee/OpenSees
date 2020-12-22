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

// ReSharper disable CppClangTidyClangDiagnosticShorten64To32
// ReSharper disable CppClangTidyClangDiagnosticSignConversion
#include "LeeNewmark.h"
#include <AnalysisModel.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <elementAPI.h>
#include <FE_EleIter.h>
#include <FE_Element.h>
#include <iostream>
#include <LinearSOE.h>
#include <sparseGEN/LeeSparse.h>

void* OPS_LeeNewmark(const bool use_initial) {
	const auto argc = OPS_GetNumRemainingInputArgs();
	if(argc % 2 != 0) {
		opserr << "WARNING - incorrect number of args want LeeNewmark\n";
		return nullptr;
	}

	double gamma, beta;
	auto num_argument = 1;

	OPS_GetDouble(&num_argument, &gamma);
	OPS_GetDouble(&num_argument, &beta);

	num_argument = argc / 2 - 1;
	std::vector<double> X(num_argument), F(num_argument);

	num_argument = 1;

	for(auto I = 0llu; I < X.size(); ++I) {
		OPS_GetDouble(&num_argument, &X[I]);
		OPS_GetDouble(&num_argument, &F[I]);
	}

	return new LeeNewmark(gamma, beta, X, F, use_initial);
}

LeeNewmark::LeeNewmark(const double _gamma, const double _beta, const std::vector<double>& X, const std::vector<double>& F, const bool _use_initial)
	: Newmark(_gamma, _beta)
	, use_initial(_use_initial) {
	if(X.size() != F.size()) throw;

	n_damping = X.size();

	mass_coef.clear();
	stiffness_coef.clear();
	mass_coef.reserve(n_damping);
	stiffness_coef.reserve(n_damping);

	for(index_t I = 0; I < n_damping; ++I) {
		mass_coef.emplace_back(4. * X[I] * F[I]);
		stiffness_coef.emplace_back(4. * X[I] / F[I]);
		omega += mass_coef.back();
	}
}

int LeeNewmark::formTangent(const int F) {
	statusFlag = F;

	auto result = 0;

	auto t_soe = dynamic_cast<LeeSparse*>(getLinearSOE());
	auto t_model = getAnalysisModel();

	if(t_soe == nullptr || t_model == nullptr) {
		opserr << "WARNING LeeNewmark::formTangent() no LeeSparse or AnalysisModel has been set\n";
		return -1;
	}

	if(first_iteration) {
		t_soe->zero_current_stiffness();
		t_soe->zero_mass();
	}
	t_soe->zero_stiffness();
	t_soe->zero_damping();

	auto& t_dof = t_model->getDOFs();
	DOF_Group* t_dofptr;

	while((t_dofptr = t_dof()) != nullptr) {
		if(first_iteration) {
			which_matrix = use_initial ? MatType::InitialStiffness : MatType::CurrentStiffness;
			if(t_soe->add_current_stiffness(t_dofptr->getTangent(this), t_dofptr->getID()) < 0) {
				opserr << "LeeNewmark::formTangent() - failed to add_stiffness:dof\n";
				result = -1;
			}
			which_matrix = MatType::Mass;
			if(t_soe->add_mass(t_dofptr->getTangent(this), t_dofptr->getID()) < 0) {
				opserr << "LeeNewmark::formTangent() - failed to add_mass:dof\n";
				result = -1;
			}
		}
		which_matrix = MatType::Stiffness;
		if(t_soe->add_stiffness(t_dofptr->getTangent(this), t_dofptr->getID()) < 0) {
			opserr << "LeeNewmark::formTangent() - failed to add_stiffness:dof\n";
			result = -1;
		}
		which_matrix = MatType::Damping;
		if(t_soe->add_damping(t_dofptr->getTangent(this), t_dofptr->getID()) < 0) {
			opserr << "LeeNewmark::formTangent() - failed to add_damping:dof\n";
			result = -1;
		}
	}

	auto& t_ele = t_model->getFEs();
	FE_Element* t_eleptr;
	while((t_eleptr = t_ele()) != nullptr) {
		if(first_iteration) {
			which_matrix = use_initial ? MatType::InitialStiffness : MatType::CurrentStiffness;
			if(t_soe->add_current_stiffness(t_eleptr->getTangent(this), t_eleptr->getID()) < 0) {
				opserr << "LeeNewmark::formTangent() - failed to add_stiffness:ele\n";
				result = -2;
			}
			which_matrix = MatType::Mass;
			if(t_soe->add_mass(t_eleptr->getTangent(this), t_eleptr->getID()) < 0) {
				opserr << "LeeNewmark::formTangent() - failed to add_mass:ele\n";
				result = -2;
			}
		}
		which_matrix = MatType::Stiffness;
		if(t_soe->add_stiffness(t_eleptr->getTangent(this), t_eleptr->getID()) < 0) {
			opserr << "LeeNewmark::formTangent() - failed to add_stiffness:ele\n";
			result = -2;
		}
		which_matrix = MatType::Damping;
		if(t_soe->add_damping(t_eleptr->getTangent(this), t_eleptr->getID()) < 0) {
			opserr << "LeeNewmark::formTangent() - failed to add_damping:ele\n";
			result = -2;
		}
	}

	const auto t_size = n_damping * n_block + n_block;

	auto& stiffness = t_soe->global_stiffness;

	if(first_iteration) stiffness = triplet_form<double>(t_size, t_size, (4 * n_damping + 2) * t_soe->stiffness.n_elem);
	else {
		if(!stiffness.csc_sort()) return -1;

		const auto& row = stiffness.row_idx;
		const auto& col = stiffness.col_idx;
		const auto& val = stiffness.val_idx;

		for(index_t I = 0; I < stiffness.c_size; ++I) {
			// quit if current column is beyond the original size of matrix
			if(col[I] >= n_block) break;
			// erase existing entries if fall in intact stiffness matrix
			if(row[I] < n_block) val[I] = 0.;
		}
	}

	t_soe->stiffness += t_soe->mass * C6 + t_soe->damping * C7;

	// check in tangent stiffness
	if(first_iteration) {
		t_soe->current_stiffness.csc_condense();
		t_soe->mass.csc_condense();
	}
	t_soe->stiffness.csc_condense();
	t_soe->damping.csc_condense();

	stiffness += t_soe->stiffness;

	if(first_iteration) {
		// for the first iteration of each substep
		// store current stiffness to be used in the whole substep
		// check in constant terms that does not change in the substep

		first_iteration = false;

		for(index_t I = 0, J = n_block; I < n_damping; ++I, J += n_block) {
			auto row = t_soe->mass.row_idx;
			auto col = t_soe->mass.col_idx;
			auto val = t_soe->mass.val_idx;
			for(index_t O = 0; O < t_soe->mass.c_size; ++O) {
				const auto K = row[O], L = col[O], M = K + J, N = L + J;
				stiffness.at(M, L) = C7 * (stiffness.at(K, N) = -(stiffness.at(M, N) = mass_coef[I] * val[O]));
			}
			row = t_soe->current_stiffness.row_idx;
			col = t_soe->current_stiffness.col_idx;
			val = t_soe->current_stiffness.val_idx;
			for(index_t K = 0; K < t_soe->current_stiffness.c_size; ++K) stiffness.at(row[K] + J, col[K] + J) = stiffness_coef[I] * val[K];
		}
	}

	return result;
}

int LeeNewmark::formEleTangent(FE_Element* theEle) {
	theEle->zeroTangent();

	if(MatType::InitialStiffness == which_matrix) theEle->addKiToTang();
	else if(MatType::Stiffness == which_matrix || MatType::CurrentStiffness == which_matrix)
		if(statusFlag == CURRENT_TANGENT) theEle->addKtToTang();
		else if(statusFlag == INITIAL_TANGENT) theEle->addKiToTang();
		else if(statusFlag == HALL_TANGENT) {
			theEle->addKtToTang(cFactor);
			theEle->addKiToTang(iFactor);
		} else { opserr << "Newmark::formEleTangent - unknown FLAG\n"; }
	else if(MatType::Mass == which_matrix) theEle->addMtoTang();
	else if(MatType::Damping == which_matrix) theEle->addCtoTang();

	return 0;
}

int LeeNewmark::formNodTangent(DOF_Group* theDof) {
	theDof->zeroTangent();

	if(MatType::Mass == which_matrix) theDof->addMtoTang();
	else if(MatType::Damping == which_matrix) theDof->addCtoTang();

	return 0;
}

int LeeNewmark::formUnbalance() {
	const auto flag = TransientIntegrator::formUnbalance();
	if(0 != flag) return flag;

	const auto t_soe = dynamic_cast<LeeSparse*>(getLinearSOE());

	if(nullptr == t_soe) {
		opserr << "WARNING LeeNewmark::formUnbalance() need a LeeSparse system\n";
		return -1;
	}

	auto& residual = t_soe->residual;

	if(0 == residual.Size()) residual.resize(n_damping * n_block + n_block);

	residual.Zero();
	Vector t_residual;
	t_residual.setData(&residual(0), n_block);
	t_residual = getLinearSOE()->getB();

	Vector internal_velocity = omega * *Udot;
	for(index_t I = 0, J = n_block; I < n_damping; ++I, J += n_block) {
		Vector n_internal;
		n_internal.setData(&trial_internal(J), n_block);
		internal_velocity -= mass_coef[I] * n_internal;
		t_residual.setData(&residual(J), n_block);
		t_residual = t_soe->mass * (*Udot - n_internal) * mass_coef[I] - t_soe->current_stiffness * n_internal * stiffness_coef[I];
	}

	t_residual.setData(&residual(0), n_block);
	t_residual -= t_soe->mass * internal_velocity;

	return 0;
}

int LeeNewmark::commit() {
	first_iteration = true;

	current_internal = trial_internal;

	return IncrementalIntegrator::commit();
}

int LeeNewmark::revertToLastStep() {
	first_iteration = true;

	trial_internal = current_internal;

	return Newmark::revertToLastStep();
}

int LeeNewmark::revertToStart() {
	first_iteration = true;

	trial_internal.Zero();
	current_internal.Zero();

	return Newmark::revertToStart();
}

int LeeNewmark::domainChanged() {
	const auto flag = Newmark::domainChanged();

	if(flag != 0) return flag;

	n_block = getLinearSOE()->getX().Size();

	const auto n_size = n_damping * n_block + n_block;

	first_iteration = true;

	current_internal = Vector(n_size);
	trial_internal = Vector(n_size);

	return 0;
}

int LeeNewmark::update(const Vector& DU) {
	const auto t_soe = dynamic_cast<LeeSparse*>(getLinearSOE());

	trial_internal += t_soe->increment;

	return Newmark::update(DU);
}

int LeeNewmark::newStep(const double DT) {
	const auto flag = Newmark::newStep(DT);

	if(0 != flag) return flag;

	C7 = gamma / (beta * DT);
	C6 = 1. / (beta * DT * DT) + omega * C7;

	return 0;
}
