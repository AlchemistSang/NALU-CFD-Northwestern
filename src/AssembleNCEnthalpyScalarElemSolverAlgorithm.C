/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
//#include <AssembleNCScalarElemSolverAlgorithm.h>
#include <AssembleNCEnthalpyScalarElemSolverAlgorithm.h>
#include <EquationSystem.h>
#include <SolverAlgorithm.h>

#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <PecletFunction.h>
#include <Realm.h>
#include <SupplementalAlgorithm.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra {
namespace nalu {

        //==========================================================================
        // Class Definition
        //==========================================================================
        // AssembleNCEnthalpyScalarElemSolverAlgorithm - add LHS/RHS for scalar
        //==========================================================================
        //--------------------------------------------------------------------------
        //-------- constructor -----------------------------------------------------
        //--------------------------------------------------------------------------
        AssembleNCEnthalpyScalarElemSolverAlgorithm::AssembleNCEnthalpyScalarElemSolverAlgorithm(
            Realm& realm,
            stk::mesh::Part* part,
            EquationSystem* eqSystem,
            ScalarFieldType* scalarQ,
            ScalarFieldType* scalarCp,
            ScalarFieldType* scalarYS,
            ScalarFieldType* scalarYL,
            ScalarFieldType* scalarTS,
            ScalarFieldType* scalarTL,
            ScalarFieldType* scalarFL,
            ScalarFieldType* scalarTemp,
            ScalarFieldType* scalarHl,
            VectorFieldType* dqdx,
            ScalarFieldType* diffFluxCoeff)
            : SolverAlgorithm(realm, part, eqSystem),
            meshMotion_(realm_.does_mesh_move()),
            scalarQ_(scalarQ),
            scalarCp_(scalarCp),
            scalarYS_(scalarYS),
            scalarYL_(scalarYL),
            scalarTS_(scalarTS),
            scalarTL_(scalarTL),
            scalarFL_(scalarFL),
            scalarTemp_(scalarTemp),
            scalarHl_(scalarHl),
            dqdx_(dqdx),
            diffFluxCoeff_(diffFluxCoeff),
            velocityRTM_(NULL),
            coordinates_(NULL),
            density_(NULL),
            massFlowRate_(NULL),
            pecletFunction_(NULL)
        {
            // save off fields
            stk::mesh::MetaData& meta_data = realm_.meta_data();
            if (meshMotion_)
                velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_rtm");
            else
                velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
            coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
            density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
            massFlowRate_ = meta_data.get_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "mass_flow_rate_scs");

            // create the peclet blending function
            pecletFunction_ = eqSystem->create_peclet_function(scalarQ_->name());

            /* Notes:

            Matrix layout is in row major. For a npe = 4 (quad) and nDof = 1:

            RHS = (resQ0, resQ1, resQ2, resQ3)

            The LHS is, therefore,

            row 0: d/dQ0(ResQ0), ., ., ., .,  d/dQ3(ResQ0)
            row 1: d/dQ0(ResQ1), ., ., ., .,  d/dQ3(ResQ1)

            */
        }

        //--------------------------------------------------------------------------
        //-------- destructor ------------------------------------------------------
        //--------------------------------------------------------------------------
        AssembleNCEnthalpyScalarElemSolverAlgorithm::~AssembleNCEnthalpyScalarElemSolverAlgorithm()
        {
            delete pecletFunction_;
        }

        //--------------------------------------------------------------------------
        //-------- initialize_connectivity -----------------------------------------
        //--------------------------------------------------------------------------
        void
            AssembleNCEnthalpyScalarElemSolverAlgorithm::initialize_connectivity()
        {
            eqSystem_->linsys_->buildElemToNodeGraph(partVec_);
        }

        //--------------------------------------------------------------------------
        //-------- execute ---------------------------------------------------------
        //--------------------------------------------------------------------------
        void
            AssembleNCEnthalpyScalarElemSolverAlgorithm::execute()
        {

            stk::mesh::BulkData& bulk_data = realm_.bulk_data();
            stk::mesh::MetaData& meta_data = realm_.meta_data();

            const int nDim = meta_data.spatial_dimension();
            const double small = 1.0e-16;

            double LA = 397e3;
            double LB = 185e3;
            const double rhoA = 2590.5;
            const double rhoB = 6510.0;
            const double cpA = 1176.8;
            const double cpB = 436.2;

            // extract user advection options (allow to potentially change over time)
            const std::string dofName = scalarQ_->name();
            const double alpha = realm_.get_alpha_factor(dofName);
            const double alphaUpw = realm_.get_alpha_upw_factor(dofName);
            const double hoUpwind = realm_.get_upw_factor(dofName);
            const bool useLimiter = realm_.primitive_uses_limiter(dofName);

            // one minus flavor..
            const double om_alpha = 1.0 - alpha;
            const double om_alphaUpw = 1.0 - alphaUpw;

            // space for LHS/RHS; nodesPerElem*nodesPerElem* and nodesPerElem
            std::vector<double> lhs;
            std::vector<double> rhs;
            std::vector<int> scratchIds;
            std::vector<double> scratchVals;
            std::vector<stk::mesh::Entity> connected_nodes;

            // supplemental algorithm setup
            const size_t supplementalAlgSize = supplementalAlg_.size();
            for (size_t i = 0; i < supplementalAlgSize; ++i)
                supplementalAlg_[i]->setup();

            // nodal fields to gather
            std::vector<double> ws_vrtm;
            std::vector<double> ws_coordinates;
            std::vector<double> ws_scalarQNp1;
            std::vector<double> ws_scalarCp;
            std::vector<double> ws_scalarYS;
            std::vector<double> ws_scalarYL;
            std::vector<double> ws_scalarTS;
            std::vector<double> ws_scalarTL;
            std::vector<double> ws_scalarFL;
            std::vector<double> ws_scalarTemp;
            std::vector<double> ws_scalarHl;
            std::vector<double> ws_dqdx;
            std::vector<double> ws_density;
            std::vector<double> ws_diffFluxCoeff;

            // geometry related to populate
            std::vector<double> ws_scs_areav;
            std::vector<double> ws_dndx;
            std::vector<double> ws_deriv;
            std::vector<double> ws_det_j;
            std::vector<double> ws_shape_function;

            // ip values
            std::vector<double>coordIp(nDim);

            // pointers
            double* p_coordIp = &coordIp[0];

            // deal with state
            ScalarFieldType& scalarQNp1 = scalarQ_->field_of_state(stk::mesh::StateNP1);
            ScalarFieldType& densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

            // define some common selectors
            stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
                & stk::mesh::selectUnion(partVec_)
                & !(realm_.get_inactive_selector());

            stk::mesh::BucketVector const& elem_buckets =
                realm_.get_buckets(stk::topology::ELEMENT_RANK, s_locally_owned_union);
            for (stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
                ib != elem_buckets.end(); ++ib) {
                stk::mesh::Bucket& b = **ib;
                const stk::mesh::Bucket::size_type length = b.size();

                // extract master element
                MasterElement* meSCS = realm_.get_surface_master_element(b.topology());
                MasterElement* meSCV = realm_.get_volume_master_element(b.topology());

                // extract master element specifics
                const int nodesPerElement = meSCS->nodesPerElement_;
                const int numScsIp = meSCS->numIntPoints_;
                const int* lrscv = meSCS->adjacentNodes();

                // resize some things; matrix related
                const int lhsSize = nodesPerElement * nodesPerElement;
                const int rhsSize = nodesPerElement;
                lhs.resize(lhsSize);
                rhs.resize(rhsSize);
                scratchIds.resize(rhsSize);
                scratchVals.resize(rhsSize);
                connected_nodes.resize(nodesPerElement);

                // algorithm related
                ws_vrtm.resize(nodesPerElement * nDim);
                ws_coordinates.resize(nodesPerElement * nDim);
                ws_dqdx.resize(nodesPerElement * nDim);
                ws_scalarQNp1.resize(nodesPerElement);
                ws_scalarCp.resize(nodesPerElement);
                ws_scalarYS.resize(nodesPerElement);
                ws_scalarYL.resize(nodesPerElement);
                ws_scalarTS.resize(nodesPerElement);
                ws_scalarTL.resize(nodesPerElement);
                ws_scalarFL.resize(nodesPerElement);
                ws_scalarTemp.resize(nodesPerElement);
                ws_scalarHl.resize(nodesPerElement);
                ws_density.resize(nodesPerElement);
                ws_diffFluxCoeff.resize(nodesPerElement);
                ws_scs_areav.resize(numScsIp * nDim);
                ws_dndx.resize(nDim * numScsIp * nodesPerElement);
                ws_deriv.resize(nDim * numScsIp * nodesPerElement);
                ws_det_j.resize(numScsIp);
                ws_shape_function.resize(numScsIp * nodesPerElement);

                // pointer to lhs/rhs
                double* p_lhs = &lhs[0];
                double* p_rhs = &rhs[0];
                double* p_vrtm = &ws_vrtm[0];
                double* p_coordinates = &ws_coordinates[0];
                double* p_dqdx = &ws_dqdx[0];
                double* p_scalarQNp1 = &ws_scalarQNp1[0];
                double* p_scalarCp = &ws_scalarCp[0];
                double* p_scalarYS = &ws_scalarYS[0];
                double* p_scalarYL = &ws_scalarYL[0];
                double* p_scalarTS = &ws_scalarTS[0];
                double* p_scalarTL = &ws_scalarTL[0];
                double* p_scalarFL = &ws_scalarFL[0];
                double* p_scalarTemp = &ws_scalarTemp[0];
                double* p_scalarHl = &ws_scalarHl[0];
                double* p_density = &ws_density[0];
                double* p_diffFluxCoeff = &ws_diffFluxCoeff[0];
                double* p_scs_areav = &ws_scs_areav[0];
                double* p_dndx = &ws_dndx[0];
                double* p_shape_function = &ws_shape_function[0];
                double p_ubar[3];

                // extract shape function
                meSCS->shape_fcn(&p_shape_function[0]);

                // resize possible supplemental element alg
                for (size_t i = 0; i < supplementalAlgSize; ++i)
                    supplementalAlg_[i]->elem_resize(meSCS, meSCV);

                for (stk::mesh::Bucket::size_type k = 0; k < length; ++k) {

                    // get elem
                    stk::mesh::Entity elem = b[k];

                    // zero lhs/rhs
                    for (int p = 0; p < lhsSize; ++p)
                        p_lhs[p] = 0.0;
                    for (int p = 0; p < rhsSize; ++p)
                        p_rhs[p] = 0.0;


                    // ip data for this element; scs and scv
                    const double* mdot = stk::mesh::field_data(*massFlowRate_, elem);

                    //===============================================
                    // gather nodal data; this is how we do it now..
                    //===============================================
                    stk::mesh::Entity const* node_rels = bulk_data.begin_nodes(elem);
                    int num_nodes = bulk_data.num_nodes(elem);

                    // sanity check on num nodes
                    ThrowAssert(num_nodes == nodesPerElement);

                    for (int ni = 0; ni < num_nodes; ++ni) {
                        stk::mesh::Entity node = node_rels[ni];

                        // set connected nodes
                        connected_nodes[ni] = node;

                        // pointers to real data
                        const double* vrtm = stk::mesh::field_data(*velocityRTM_, node);
                        const double* coords = stk::mesh::field_data(*coordinates_, node);
                        const double* dq = stk::mesh::field_data(*dqdx_, node);

                        // gather scalars
                        p_scalarQNp1[ni] = *stk::mesh::field_data(scalarQNp1, node);
                        p_density[ni] = *stk::mesh::field_data(densityNp1, node);
                        p_diffFluxCoeff[ni] = *stk::mesh::field_data(*diffFluxCoeff_, node);
                        p_scalarCp[ni] = *stk::mesh::field_data(*scalarCp_, node);
                        p_scalarYS[ni] = *stk::mesh::field_data(*scalarYS_, node);
                        p_scalarYL[ni] = *stk::mesh::field_data(*scalarYL_, node);
                        p_scalarTS[ni] = *stk::mesh::field_data(*scalarTS_, node);
                        p_scalarTL[ni] = *stk::mesh::field_data(*scalarTL_, node);
                        p_scalarFL[ni] = *stk::mesh::field_data(*scalarFL_, node);
                        p_scalarTemp[ni] = *stk::mesh::field_data(*scalarTemp_, node);
                        p_scalarHl[ni] = *stk::mesh::field_data(*scalarHl_, node);

                        // gather vectors
                        const int niNdim = ni * nDim;
                        for (int i = 0; i < nDim; ++i) {
                            p_vrtm[niNdim + i] = vrtm[i];
                            p_coordinates[niNdim + i] = coords[i];
                            p_dqdx[niNdim + i] = dq[i];
                        }
                    }

                    // compute geometry
                    double scs_error = 0.0;
                    meSCS->determinant(1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

                    // compute dndx
                    meSCS->grad_op(1, &p_coordinates[0], &p_dndx[0], &ws_deriv[0], &ws_det_j[0], &scs_error);

                    for (int ip = 0; ip < numScsIp; ++ip) {

                        // left and right nodes for this ip
                        const int il = lrscv[2 * ip];
                        const int ir = lrscv[2 * ip + 1];

                        // Non-conservative formulation
                        const double rhoL = p_density[il];
                        const double rhoR = p_density[ir];

                        // corresponding matrix rows
                        const int rowL = il * nodesPerElement;
                        const int rowR = ir * nodesPerElement;

                        // save off mdot FIXME
                        const double tmdot = mdot[ip];

                        // zero out values of interest for this ip
                        for (int j = 0; j < nDim; ++j) {
                            p_coordIp[j] = 0.0;
                        }

                        // save off ip values; offset to Shape Function
                        double rhoIp = 0.0;
                        double muIp = 0.0;
                        double qIp = 0.0;
                        const int offSetSF = ip * nodesPerElement;
                        for (int ic = 0; ic < nodesPerElement; ++ic) {
                            const double r = p_shape_function[offSetSF + ic];
                            rhoIp += r * p_density[ic];
                            muIp += r * p_diffFluxCoeff[ic];
                            qIp += r * p_scalarHl[ic];
                            //qIp += r * p_density[ic] * p_scalarHl[ic];
                        
                            // compute scs point values
                            for (int i = 0; i < nDim; ++i) {
                                p_coordIp[i] += r * p_coordinates[ic * nDim + i];
                            }
                        }

                        // FIXME: Hack in mass flow rate
                        /* p_ubar[0] = 2.0 * M_PI * p_coordIp[1];
                        p_ubar[1] = -2.0 * M_PI * p_coordIp[0];
                        p_ubar[2] = 0.0;
                        double tmdot = 0.0;
                        for (int j = 0; j < nDim; j++)
                        {
                          tmdot += p_ubar[j] * p_scs_areav[ip*nDim+j];
                        }
                        // END FIX HERE*/

                        // Peclet factor; along the edge
                        const double diffIp = 0.5 * (p_diffFluxCoeff[il] / p_density[il]
                            + p_diffFluxCoeff[ir] / p_density[ir]);
                        //const double diffIp = 0.5 * (p_diffFluxCoeff[il] + p_diffFluxCoeff[ir]);
                        double udotx = 0.0;
                        for (int j = 0; j < nDim; ++j) {
                            const double dxj = p_coordinates[ir * nDim + j] - p_coordinates[il * nDim + j];
                            const double uj = 0.5 * (p_vrtm[il * nDim + j] + p_vrtm[ir * nDim + j]);
                            udotx += uj * dxj;
                        }

                        // FIXME: HACK IN pecfac
                        //const double pecfac = 1.0;
                        const double pecfac = pecletFunction_->execute(std::abs(udotx) / (diffIp + small));
                        const double om_pecfac = 1.0 - pecfac;

                        // left and right extrapolation
                        double dqL = 0.0;
                        double dqR = 0.0;
                        for (int j = 0; j < nDim; ++j) {
                            const double dxjL = p_coordIp[j] - p_coordinates[il * nDim + j];
                            const double dxjR = p_coordinates[ir * nDim + j] - p_coordIp[j];
                            dqL += dxjL * p_dqdx[nDim * il + j];
                            dqR += dxjR * p_dqdx[nDim * ir + j];
                            //dqL += dxjL * rhoL * p_dqdx[nDim * il + j];
                            //dqR += dxjR * rhoR * p_dqdx[nDim * ir + j];
                        }

                        // add limiter if appropriate
                        double limitL = 1.0;
                        double limitR = 1.0;
                        if (useLimiter) {
                            const double dq = p_scalarHl[ir] - p_scalarHl[il];
                            //const double dq = rhoR * p_scalarHl[ir] - rhoL * p_scalarHl[il];
                            const double dqMl = 2.0 * 2.0 * dqL - dq;
                            const double dqMr = 2.0 * 2.0 * dqR - dq;
                            limitL = van_leer(dqMl, dq, small);
                            limitR = van_leer(dqMr, dq, small);
                        }

                        // extrapolated; for now limit (along edge is fine)
                        const double qIpL = p_scalarHl[il] + dqL * hoUpwind * limitL;
                        const double qIpR = p_scalarHl[ir] - dqR * hoUpwind * limitR;
                        //const double qIpL = rhoL * p_scalarHl[il] + dqL * hoUpwind * limitL;
                        //const double qIpR = rhoR * p_scalarHl[ir] - dqR * hoUpwind * limitR;
                        
                        // assemble advection; rhs and upwind contributions

                        // 2nd order central; simply qIp from above

                        // upwind
                        const double qUpwind = (tmdot > 0) ? alphaUpw * qIpL + om_alphaUpw * qIp
                            : alphaUpw * qIpR + om_alphaUpw * qIp;

                        // generalized central (2nd and 4th order)
                        const double qHatL = alpha * qIpL + om_alpha * qIp;
                        const double qHatR = alpha * qIpR + om_alpha * qIp;
                        const double qCds = 0.5 * (qHatL + qHatR);

                        // total advection
                        const double aflux = tmdot * (pecfac * qUpwind + om_pecfac * qCds);

                        // right hand side; L and R
                        p_rhs[il] -= rhoL * aflux;
                        p_rhs[ir] += rhoR * aflux; 
                        //p_rhs[il] -= aflux;
                        //p_rhs[ir] += aflux;

                        // advection operator sens; all but central

                        // upwind advection (includes 4th); left node
                        const double alhsfacL = 0.5 * (tmdot + std::abs(tmdot)) * pecfac * alphaUpw
                            + 0.5 * alpha * om_pecfac * tmdot;

                        p_lhs[rowL + il] += alhsfacL * rhoL;
                        p_lhs[rowR + il] -= alhsfacL * rhoR;
                        //p_lhs[rowL + il] += alhsfacL;
                        //p_lhs[rowR + il] -= alhsfacL;

                        if (p_scalarTemp[il] > p_scalarTS[il] && p_scalarTemp[il] < p_scalarTL[il]) {
                            //Go back
                            p_lhs[rowL + il] -= alhsfacL * rhoL;
                            p_lhs[rowR + il] += alhsfacL * rhoR;
                            //p_lhs[rowL + il] -= alhsfacL;
                            //p_lhs[rowR + il] += alhsfacL;

                            //p_lhs[rowL + il] += alhsfacL * rhoL / (1.0 + (rhoA * LA * (1.0 - p_scalarYL[il]) + rhoB * LB * p_scalarYL[il]) * p_scalarFL[il] / (rhoA * cpA * (1.0 - p_scalarYL[il]) + rhoB * cpB * p_scalarYL[il]));
                            //p_lhs[rowR + il] -= alhsfacL * rhoR / (1.0 + (rhoA * LA * (1.0 - p_scalarYL[il]) + rhoB * LB * p_scalarYL[il]) * p_scalarFL[il] / (rhoA * cpA * (1.0 - p_scalarYL[il]) + rhoB * cpB * p_scalarYL[il]));
                            p_lhs[rowL + il] += alhsfacL * rhoL / (1.0 + (LA * (1.0 - p_scalarYL[il]) + LB * p_scalarYL[il]) * p_scalarFL[il] / (cpA * (1.0 - p_scalarYL[il]) + cpB * p_scalarYL[il]));
                            p_lhs[rowR + il] -= alhsfacL * rhoR / (1.0 + (LA * (1.0 - p_scalarYL[il]) + LB * p_scalarYL[il]) * p_scalarFL[il] / (cpA * (1.0 - p_scalarYL[il]) + cpB * p_scalarYL[il]));
                        }

                        // upwind advection; right node
                        const double alhsfacR = 0.5 * (tmdot - std::abs(tmdot)) * pecfac * alphaUpw
                            + 0.5 * alpha * om_pecfac * tmdot;

                        p_lhs[rowR + ir] -= alhsfacR * rhoR;
                        p_lhs[rowL + ir] += alhsfacR * rhoL;
                        //p_lhs[rowR + ir] -= alhsfacR;
                        //p_lhs[rowL + ir] += alhsfacR;

                        if (p_scalarTemp[ir] > p_scalarTS[ir] && p_scalarTemp[ir] < p_scalarTL[ir]) {
                            //GO back
                            p_lhs[rowR + ir] += alhsfacR * rhoR;
                            p_lhs[rowL + ir] -= alhsfacR * rhoL;
                            //p_lhs[rowR + ir] += alhsfacR;
                            //p_lhs[rowL + ir] -= alhsfacR;

                            //p_lhs[rowR + ir] -= alhsfacR * rhoR / (1.0 + (rhoA * LA * (1.0 - p_scalarYL[ir]) + rhoB * LB * p_scalarYL[ir]) * p_scalarFL[ir] / (rhoA * cpA * (1.0 - p_scalarYL[ir]) + rhoB * cpB * p_scalarYL[ir]));
                            //p_lhs[rowL + ir] += alhsfacR * rhoL / (1.0 + (rhoA * LA * (1.0 - p_scalarYL[ir]) + rhoB * LB * p_scalarYL[ir]) * p_scalarFL[ir] / (rhoA * cpA * (1.0 - p_scalarYL[ir]) + rhoB * cpB * p_scalarYL[ir]));
                            p_lhs[rowR + ir] -= alhsfacR * rhoR / (1.0 + (LA * (1.0 - p_scalarYL[ir]) + LB * p_scalarYL[ir]) * p_scalarFL[ir] / (cpA * (1.0 - p_scalarYL[ir]) + cpB * p_scalarYL[ir]));
                            p_lhs[rowL + ir] += alhsfacR * rhoL / (1.0 + (LA * (1.0 - p_scalarYL[ir]) + LB * p_scalarYL[ir]) * p_scalarFL[ir] / (cpA * (1.0 - p_scalarYL[ir]) + cpB * p_scalarYL[ir]));
                        }

                        double qDiff = 0.0;
                        for (int ic = 0; ic < nodesPerElement; ++ic) {

                            // shape function
                            const double r = p_shape_function[offSetSF + ic];

                            // upwind (il/ir) handled above; collect terms on alpha and alphaUpw
                            const double lhsfacAdv = r * tmdot * (pecfac * om_alphaUpw + om_pecfac * om_alpha);

                            // advection operator lhs; rhs handled above
                            // lhs; il then ir
                            p_lhs[rowL + ic] += lhsfacAdv * rhoL;
                            p_lhs[rowR + ic] -= lhsfacAdv * rhoR;
                            //p_lhs[rowL + ic] += lhsfacAdv;
                            //p_lhs[rowR + ic] -= lhsfacAdv;

                            if (p_scalarTemp[ic] > p_scalarTS[ic] && p_scalarTemp[ic] < p_scalarTL[ic]) {
                                // Go back
                                p_lhs[rowL + ic] -= lhsfacAdv * rhoL;
                                p_lhs[rowR + ic] += lhsfacAdv * rhoR;
                                //p_lhs[rowL + ic] -= lhsfacAdv;
                                //p_lhs[rowR + ic] += lhsfacAdv;

                                //p_lhs[rowL + ic] += lhsfacAdv * rhoL / (1.0 + (rhoA * LA * (1.0 - p_scalarYL[ic]) + rhoB * LB * p_scalarYL[ic]) * p_scalarFL[ic] / (rhoA * cpA * (1.0 - p_scalarYL[ic]) + rhoB * cpB * p_scalarYL[ic]));
                                //p_lhs[rowR + ic] -= lhsfacAdv * rhoR / (1.0 + (rhoA * LA * (1.0 - p_scalarYL[ic]) + rhoB * LB * p_scalarYL[ic]) * p_scalarFL[ic] / (rhoA * cpA * (1.0 - p_scalarYL[ic]) + rhoB * cpB * p_scalarYL[ic]));
                                p_lhs[rowL + ic] += lhsfacAdv * rhoL / (1.0 + (LA * (1.0 - p_scalarYL[ic]) + LB * p_scalarYL[ic]) * p_scalarFL[ic] / (cpA * (1.0 - p_scalarYL[ic]) + cpB * p_scalarYL[ic]));
                                p_lhs[rowR + ic] -= lhsfacAdv * rhoR / (1.0 + (LA * (1.0 - p_scalarYL[ic]) + LB * p_scalarYL[ic]) * p_scalarFL[ic] / (cpA * (1.0 - p_scalarYL[ic]) + cpB * p_scalarYL[ic]));
                            }

                            // diffusion
                            double lhsfacDiff = 0.0;

                            const int offSetDnDx = nDim * nodesPerElement * ip + ic * nDim;
                            for (int j = 0; j < nDim; ++j) {
                                lhsfacDiff += -muIp * p_dndx[offSetDnDx + j] * p_scs_areav[ip * nDim + j];
                            }

                            qDiff += lhsfacDiff * p_scalarTemp[ic];
                            //qDiff += lhsfacDiff * p_scalarQNp1[ic];

                            p_lhs[rowL + ic] += lhsfacDiff / (cpA * (1.0 - p_scalarYL[ic]) + cpB * p_scalarYL[ic]);
                            p_lhs[rowR + ic] -= lhsfacDiff / (cpA * (1.0 - p_scalarYL[ic]) + cpB * p_scalarYL[ic]);
                            //// lhs; il then ir
                            if (p_scalarTemp[ic] > p_scalarTS[ic] && p_scalarTemp[ic] < p_scalarTL[ic]) {
                                // Go back first
                                p_lhs[rowL + ic] -= lhsfacDiff / (cpA * (1.0 - p_scalarYL[ic]) + cpB * p_scalarYL[ic]);
                                p_lhs[rowR + ic] += lhsfacDiff / (cpA * (1.0 - p_scalarYL[ic]) + cpB * p_scalarYL[ic]);

                                p_lhs[rowL + ic] += lhsfacDiff / ((cpA * (1.0 - p_scalarYS[ic]) + cpB * p_scalarYS[ic]) + (LA * (1.0 - p_scalarYL[ic]) + LB * p_scalarYL[ic]) / (p_scalarTL[ic] - p_scalarTS[ic]));
                                p_lhs[rowR + ic] -= lhsfacDiff / ((cpA * (1.0 - p_scalarYS[ic]) + cpB * p_scalarYS[ic]) + (LA * (1.0 - p_scalarYL[ic]) + LB * p_scalarYL[ic]) / (p_scalarTL[ic] - p_scalarTS[ic]));
                            }
                        }

                        // rhs; il then ir
                        p_rhs[il] -= qDiff;
                        p_rhs[ir] += qDiff;

                    }

                    // call supplemental
                    for (size_t i = 0; i < supplementalAlgSize; ++i)
                        supplementalAlg_[i]->elem_execute(&lhs[0], &rhs[0], elem, meSCS, meSCV);

                    apply_coeff(connected_nodes, scratchIds, scratchVals, rhs, lhs, __FILE__);

                }
            }
        }

        //--------------------------------------------------------------------------
        //-------- van_leer ---------------------------------------------------------
        //--------------------------------------------------------------------------
        double
            AssembleNCEnthalpyScalarElemSolverAlgorithm::van_leer(
                const double& dqm,
                const double& dqp,
                const double& small)
        {
            double limit = (2.0 * (dqm * dqp + std::abs(dqm * dqp))) /
                ((dqm + dqp) * (dqm + dqp) + small);
            return limit;
        }

    } // namespace nalu
} // namespace Sierra
