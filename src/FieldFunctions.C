/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <FieldFunctions.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Field.hpp>

#include <algorithm>

namespace sierra {
namespace nalu {
 
void field_axpby(
  const stk::mesh::MetaData & metaData,
  const stk::mesh::BulkData & bulkData,
  const double alpha,
  const stk::mesh::FieldBase & xField,
  const double beta,
  const stk::mesh::FieldBase & yField,
  const bool auraIsActive,
  const stk::topology::rank_t entityRankValue)
{
  // decide on selector
  const stk::mesh::Selector selector = auraIsActive 
    ? metaData.universal_part() &
    stk::mesh::selectField(xField) &
    stk::mesh::selectField(yField)
    : (metaData.locally_owned_part() | metaData.globally_shared_part()) &
    stk::mesh::selectField(xField) &
    stk::mesh::selectField(yField);
 
  stk::mesh::BucketVector const& buckets = bulkData.get_buckets( entityRankValue, selector );

  for(size_t i=0; i < buckets.size(); ++i) {
    stk::mesh::Bucket & b = *buckets[i];
    const stk::mesh::Bucket::size_type length = b.size();
    const size_t fieldSize = field_bytes_per_entity(xField, b) / sizeof(double);
    ThrowAssert(fieldSize == field_bytes_per_entity(yField, b) / sizeof(double));
    const unsigned kmax = length * fieldSize;
    const double * x = (double*)stk::mesh::field_data(xField, b);
    double * y = (double*)stk::mesh::field_data(yField, b);
    for(unsigned k = 0 ; k < kmax ; ++k) {
      y[k] = alpha * x[k] + beta*y[k];
    }
  }
}

void field_dTdh(
    const stk::mesh::MetaData& metaData,
    const stk::mesh::BulkData& bulkData,
    const double alpha,
    const double rhoA,
    const double rhoB,
    const double cpA,
    const double cpB,
    const double LA,
    const double LB,
    const stk::mesh::FieldBase& dTField,
    const stk::mesh::FieldBase& dhField,
    const stk::mesh::FieldBase& YSField,
    const stk::mesh::FieldBase& YLField,
    const stk::mesh::FieldBase& TSField,
    const stk::mesh::FieldBase& TLField,
    const stk::mesh::FieldBase& fLField,
    const stk::mesh::FieldBase& rhoField,
    const stk::mesh::FieldBase& TempField,
    const bool auraIsActive,
    const stk::topology::rank_t entityRankValue)
{
    // decide on selector
    const stk::mesh::Selector selector = auraIsActive
        ? metaData.universal_part() &
        stk::mesh::selectField(dTField) &
        stk::mesh::selectField(dhField) &
        stk::mesh::selectField(YSField) &
        stk::mesh::selectField(YLField) &
        stk::mesh::selectField(TSField) &
        stk::mesh::selectField(TLField) &
        stk::mesh::selectField(fLField) &
        stk::mesh::selectField(rhoField) &
        stk::mesh::selectField(TempField)
        : (metaData.locally_owned_part() | metaData.globally_shared_part()) &
        stk::mesh::selectField(dTField) &
        stk::mesh::selectField(dhField) &
        stk::mesh::selectField(YSField) &
        stk::mesh::selectField(YLField) &
        stk::mesh::selectField(TSField) &
        stk::mesh::selectField(TLField) &
        stk::mesh::selectField(fLField) &
        stk::mesh::selectField(rhoField) &
        stk::mesh::selectField(TempField);

    stk::mesh::BucketVector const& buckets = bulkData.get_buckets(entityRankValue, selector);

    for (size_t i = 0; i < buckets.size(); ++i) {
        stk::mesh::Bucket& b = *buckets[i];
        const stk::mesh::Bucket::size_type length = b.size();
        const size_t fieldSize = field_bytes_per_entity(dTField, b) / sizeof(double);
        ThrowAssert(fieldSize == field_bytes_per_entity(dhField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(YSField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(YLField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(TSField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(TLField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(fLField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(rhoField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(TempField, b) / sizeof(double));
        const unsigned kmax = length * fieldSize;
        double* dTdT = (double*)stk::mesh::field_data(dTField, b);
        const double *dhdh = (double*)stk::mesh::field_data(dhField, b);
        const double* YSYS = (double*)stk::mesh::field_data(YSField, b);
        const double* YLYL = (double*)stk::mesh::field_data(YLField, b);
        const double* TSTS = (double*)stk::mesh::field_data(TSField, b);
        const double* TLTL = (double*)stk::mesh::field_data(TLField, b);
        const double* fLfL = (double*)stk::mesh::field_data(fLField, b);
        const double* rhorho = (double*)stk::mesh::field_data(rhoField, b);
        const double* TempTemp = (double*)stk::mesh::field_data(TempField, b);
        for (unsigned k = 0; k < kmax; ++k) {
            //dTdT[k] = ( dhdh[k] / (cpA * ((alpha - YSYS[k]) * (alpha - fLfL[k]) + (alpha - YLYL[k]) * fLfL[k]) + cpB * ((YSYS[k]) * (alpha - fLfL[k]) + (YLYL[k]) * fLfL[k]) ) );
            dTdT[k] = dhdh[k] / (cpA * (alpha - YSYS[k]) + cpB * YSYS[k]);
            if (TempTemp[k] > TSTS[k] && TempTemp[k] < TLTL[k]) {
                //dTdT[k] = dhdh[k] / ((rhoA * cpA * (alpha - YSYS[k]) + rhoB * cpB * YSYS[k]) / rhorho[k] + (rhoA * LA * (1 - YLYL[k]) + rhoB * LB * YLYL[k]) / (TLTL[k] - TSTS[k]) / rhorho[k]);
                dTdT[k] = dhdh[k] / ((cpA * (alpha - YSYS[k]) + cpB * YSYS[k]) + (LA * (1 - YLYL[k]) + LB * YLYL[k]) / (TLTL[k] - TSTS[k]));
            }

            //y[k] = alpha * x[k] + beta * y[k];
        }
    }
}

void field_dfLdh(
    const stk::mesh::MetaData& metaData,
    const stk::mesh::BulkData& bulkData,
    const double alpha,
    const double rhoA,
    const double rhoB,
    const double cpA,
    const double cpB,
    const double LA,
    const double LB,
    const stk::mesh::FieldBase& dfLField,
    const stk::mesh::FieldBase& dhField,
    const stk::mesh::FieldBase& YSField,
    const stk::mesh::FieldBase& YLField,
    const stk::mesh::FieldBase& TSField,
    const stk::mesh::FieldBase& TLField,
    const stk::mesh::FieldBase& fLField,
    const stk::mesh::FieldBase& rhoField,
    const stk::mesh::FieldBase& TempField,
    const bool auraIsActive,
    const stk::topology::rank_t entityRankValue)
{
    // decide on selector
    const stk::mesh::Selector selector = auraIsActive
        ? metaData.universal_part() &
        stk::mesh::selectField(dfLField) &
        stk::mesh::selectField(dhField) &
        stk::mesh::selectField(YSField) &
        stk::mesh::selectField(YLField) &
        stk::mesh::selectField(TSField) &
        stk::mesh::selectField(TLField) &
        stk::mesh::selectField(fLField) &
        stk::mesh::selectField(rhoField) &
        stk::mesh::selectField(TempField)
        : (metaData.locally_owned_part() | metaData.globally_shared_part()) &
        stk::mesh::selectField(dfLField) &
        stk::mesh::selectField(dhField) &
        stk::mesh::selectField(YSField) &
        stk::mesh::selectField(YLField) &
        stk::mesh::selectField(TSField) &
        stk::mesh::selectField(TLField) &
        stk::mesh::selectField(fLField) &
        stk::mesh::selectField(rhoField) &
        stk::mesh::selectField(TempField);

    stk::mesh::BucketVector const& buckets = bulkData.get_buckets(entityRankValue, selector);

    for (size_t i = 0; i < buckets.size(); ++i) {
        stk::mesh::Bucket& b = *buckets[i];
        const stk::mesh::Bucket::size_type length = b.size();
        const size_t fieldSize = field_bytes_per_entity(dfLField, b) / sizeof(double);
        ThrowAssert(fieldSize == field_bytes_per_entity(dhField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(YSField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(YLField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(TSField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(TLField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(fLField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(rhoField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(TempField, b) / sizeof(double));
        const unsigned kmax = length * fieldSize;
        double* dfLdfL = (double*)stk::mesh::field_data(dfLField, b);
        const double* dhdh = (double*)stk::mesh::field_data(dhField, b);
        const double* YSYS = (double*)stk::mesh::field_data(YSField, b);
        const double* YLYL = (double*)stk::mesh::field_data(YLField, b);
        const double* TSTS = (double*)stk::mesh::field_data(TSField, b);
        const double* TLTL = (double*)stk::mesh::field_data(TLField, b);
        const double* fLfL = (double*)stk::mesh::field_data(fLField, b);
        const double* rhorho = (double*)stk::mesh::field_data(rhoField, b);
        const double* TempTemp = (double*)stk::mesh::field_data(TempField, b);
        for (unsigned k = 0; k < kmax; ++k) {
            dfLdfL[k] = 0.0;
            if (TempTemp[k] > TSTS[k] && TempTemp[k] < TLTL[k]) {
                //dfLdfL[k] = dhdh[k] / ((rhoA * LA * (alpha - YLYL[k]) + rhoB * LB * YLYL[k]) / rhorho[k]);
                dfLdfL[k] = dhdh[k] / ((LA * (alpha - YLYL[k]) + LB * YLYL[k]));
            }    
            //y[k] = alpha * x[k] + beta * y[k];
        }
    }
}

void field_dHldh(
    const stk::mesh::MetaData& metaData,
    const stk::mesh::BulkData& bulkData,
    const double alpha,
    const double rhoA,
    const double rhoB,
    const double cpA,
    const double cpB,
    const double LA,
    const double LB,
    const stk::mesh::FieldBase& dHlField,
    const stk::mesh::FieldBase& dhField,
    const stk::mesh::FieldBase& YSField,
    const stk::mesh::FieldBase& YLField,
    const stk::mesh::FieldBase& TSField,
    const stk::mesh::FieldBase& TLField,
    const stk::mesh::FieldBase& fLField,
    const stk::mesh::FieldBase& rhoField,
    const stk::mesh::FieldBase& TempField,
    const bool auraIsActive,
    const stk::topology::rank_t entityRankValue)
{
    // decide on selector
    const stk::mesh::Selector selector = auraIsActive
        ? metaData.universal_part() &
        stk::mesh::selectField(dHlField) &
        stk::mesh::selectField(dhField) &
        stk::mesh::selectField(YSField) &
        stk::mesh::selectField(YLField) &
        stk::mesh::selectField(TSField) &
        stk::mesh::selectField(TLField) &
        stk::mesh::selectField(fLField) &
        stk::mesh::selectField(rhoField) &
        stk::mesh::selectField(TempField)
        : (metaData.locally_owned_part() | metaData.globally_shared_part()) &
        stk::mesh::selectField(dHlField) &
        stk::mesh::selectField(dhField) &
        stk::mesh::selectField(YSField) &
        stk::mesh::selectField(YLField) &
        stk::mesh::selectField(TSField) &
        stk::mesh::selectField(TLField) &
        stk::mesh::selectField(fLField) &
        stk::mesh::selectField(rhoField) &
        stk::mesh::selectField(TempField);

    stk::mesh::BucketVector const& buckets = bulkData.get_buckets(entityRankValue, selector);

    for (size_t i = 0; i < buckets.size(); ++i) {
        stk::mesh::Bucket& b = *buckets[i];
        const stk::mesh::Bucket::size_type length = b.size();
        const size_t fieldSize = field_bytes_per_entity(dHlField, b) / sizeof(double);
        ThrowAssert(fieldSize == field_bytes_per_entity(dhField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(YSField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(YLField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(TSField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(TLField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(fLField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(rhoField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(TempField, b) / sizeof(double));
        const unsigned kmax = length * fieldSize;
        double* dHldHl = (double*)stk::mesh::field_data(dHlField, b);
        const double* dhdh = (double*)stk::mesh::field_data(dhField, b);
        const double* YSYS = (double*)stk::mesh::field_data(YSField, b);
        const double* YLYL = (double*)stk::mesh::field_data(YLField, b);
        const double* TSTS = (double*)stk::mesh::field_data(TSField, b);
        const double* TLTL = (double*)stk::mesh::field_data(TLField, b);
        const double* fLfL = (double*)stk::mesh::field_data(fLField, b);
        const double* rhorho = (double*)stk::mesh::field_data(rhoField, b);
        const double* TempTemp = (double*)stk::mesh::field_data(TempField, b);
        for (unsigned k = 0; k < kmax; ++k) {
            dHldHl[k] = 0.0;
            if (TempTemp[k] > TSTS[k] && TempTemp[k] < TLTL[k]) {
                //dHldHl[k] = dhdh[k] / (alpha + (rhoA * LA * (alpha - YLYL[k]) + rhoB * LB * YLYL[k]) / (rhoA * cpA * (alpha - YLYL[k]) + rhoB * cpB * YLYL[k]) / (TLTL[k] - TSTS[k]));
                dHldHl[k] = fLfL[k];
            }
            if (TempTemp[k] >= TLTL[k]) {
                dHldHl[k] = alpha;
            }
            //y[k] = alpha * x[k] + beta * y[k];
        }
    }
}


void field_update_fL(
    const stk::mesh::MetaData& metaData,
    const stk::mesh::BulkData& bulkData,
    const double alpha,
    const stk::mesh::FieldBase& xField,
    const double beta,
    const stk::mesh::FieldBase& yField,
    const bool auraIsActive,
    const stk::topology::rank_t entityRankValue)
{
    // decide on selector
    const stk::mesh::Selector selector = auraIsActive
        ? metaData.universal_part() &
        stk::mesh::selectField(xField) &
        stk::mesh::selectField(yField)
        : (metaData.locally_owned_part() | metaData.globally_shared_part()) &
        stk::mesh::selectField(xField) &
        stk::mesh::selectField(yField);

    stk::mesh::BucketVector const& buckets = bulkData.get_buckets(entityRankValue, selector);

    for (size_t i = 0; i < buckets.size(); ++i) {
        stk::mesh::Bucket& b = *buckets[i];
        const stk::mesh::Bucket::size_type length = b.size();
        const size_t fieldSize = field_bytes_per_entity(xField, b) / sizeof(double);
        ThrowAssert(fieldSize == field_bytes_per_entity(yField, b) / sizeof(double));
        const unsigned kmax = length * fieldSize;
        const double* x = (double*)stk::mesh::field_data(xField, b);
        double* y = (double*)stk::mesh::field_data(yField, b);
        for (unsigned k = 0; k < kmax; ++k) {
            y[k] = alpha * x[k] + beta * y[k];
            y[k] = std::min(std::max(0.0, y[k]), 1.0);
        }
    }
}

void field_update_hl(
    const stk::mesh::MetaData& metaData,
    const stk::mesh::BulkData& bulkData,
    const double alpha,
    const double rhoA,
    const double rhoB,
    const double cpA,
    const double cpB,
    const double LA,
    const double LB,
    const stk::mesh::FieldBase& YSField,
    const stk::mesh::FieldBase& YLField,
    const stk::mesh::FieldBase& fLField,
    const stk::mesh::FieldBase& rhoField,
    const stk::mesh::FieldBase& TempField,
    const stk::mesh::FieldBase& HlField,
    const bool auraIsActive,
    const stk::topology::rank_t entityRankValue)
{
    // decide on selector
    const stk::mesh::Selector selector = auraIsActive
        ? metaData.universal_part() &
        stk::mesh::selectField(YSField) &
        stk::mesh::selectField(YLField) &
        stk::mesh::selectField(fLField) &
        stk::mesh::selectField(rhoField) &
        stk::mesh::selectField(TempField)
        : (metaData.locally_owned_part() | metaData.globally_shared_part()) &
        stk::mesh::selectField(YSField) &
        stk::mesh::selectField(YLField) &
        stk::mesh::selectField(fLField) &
        stk::mesh::selectField(rhoField) &
        stk::mesh::selectField(TempField) &
        stk::mesh::selectField(HlField);

    stk::mesh::BucketVector const& buckets = bulkData.get_buckets(entityRankValue, selector);

    for (size_t i = 0; i < buckets.size(); ++i) {
        stk::mesh::Bucket& b = *buckets[i];
        const stk::mesh::Bucket::size_type length = b.size();
        const size_t fieldSize = field_bytes_per_entity(fLField, b) / sizeof(double);
        ThrowAssert(fieldSize == field_bytes_per_entity(YSField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(YLField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(rhoField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(TempField, b) / sizeof(double));
        ThrowAssert(fieldSize == field_bytes_per_entity(HlField, b) / sizeof(double));
        const unsigned kmax = length * fieldSize;
        const double* fLfL = (double*)stk::mesh::field_data(fLField, b);
        const double* YSYS = (double*)stk::mesh::field_data(YSField, b);
        const double* YLYL = (double*)stk::mesh::field_data(YLField, b);
        const double* rhorho = (double*)stk::mesh::field_data(rhoField, b);
        const double* TempTemp = (double*)stk::mesh::field_data(TempField, b);
        double* HlHl = (double*)stk::mesh::field_data(HlField, b);
        for (unsigned k = 0; k < kmax; ++k) {
            HlHl[k] = fLfL[k] * ((rhoA * LA * (alpha - YLYL[k]) + rhoB * LB * YLYL[k]) + rhoA * cpA * TempTemp[k] * (alpha - YLYL[k]) + rhoB * cpB * TempTemp[k] * YLYL[k]) / rhorho[k];
            //y[k] = alpha * x[k] + beta * y[k];
        }
    }
}

void field_fill(
  const stk::mesh::MetaData & metaData,
  const stk::mesh::BulkData & bulkData,
  const double alpha,
  const stk::mesh::FieldBase & xField,
  const bool auraIsActive,
  const stk::topology::rank_t entityRankValue)
{
  // decide on selector
  const stk::mesh::Selector selector = auraIsActive 
    ? metaData.universal_part() &
    stk::mesh::selectField(xField)
    : (metaData.locally_owned_part() | metaData.globally_shared_part()) &
    stk::mesh::selectField(xField);

  stk::mesh::BucketVector const& buckets = bulkData.get_buckets( entityRankValue, selector );

  for(size_t i=0; i < buckets.size(); ++i) {
    stk::mesh::Bucket & b = *buckets[i];
    const stk::mesh::Bucket::size_type length = b.size();
    const unsigned fieldSize = field_bytes_per_entity(xField, b) / sizeof(double);
    const unsigned kmax = length * fieldSize;
    double * x = (double*)stk::mesh::field_data(xField, b);
    std::fill(x, x + kmax, alpha);
  }
}

void field_scale(
  const stk::mesh::MetaData & metaData,
  const stk::mesh::BulkData & bulkData,
  const double alpha,
  const stk::mesh::FieldBase & xField,
  const bool auraIsActive,
  const stk::topology::rank_t entityRankValue)
{
  // decide on selector
  const stk::mesh::Selector selector = auraIsActive 
    ? metaData.universal_part() &
    stk::mesh::selectField(xField)
    : (metaData.locally_owned_part() | metaData.globally_shared_part()) &
    stk::mesh::selectField(xField);

  stk::mesh::BucketVector const& buckets = bulkData.get_buckets( entityRankValue, selector );

  for(size_t i=0; i < buckets.size(); ++i) {
    stk::mesh::Bucket & b = *buckets[i];
    const stk::mesh::Bucket::size_type length = b.size();
    const unsigned fieldSize = field_bytes_per_entity(xField, b) / sizeof(double);
    const unsigned kmax = length * fieldSize;
    double * x = (double*)stk::mesh::field_data(xField, b);
    for(unsigned k = 0 ; k < kmax ; ++k) {
      x[k] = alpha * x[k];
    }
  }
}

void field_copy(
  const stk::mesh::MetaData & metaData,
  const stk::mesh::BulkData & bulkData,
  const stk::mesh::FieldBase & xField,
  const stk::mesh::FieldBase & yField,
  const bool auraIsActive,
  const stk::topology::rank_t entityRankValue)
{
  // decide on selector
  const stk::mesh::Selector selector = auraIsActive 
    ? metaData.universal_part() &
    stk::mesh::selectField(xField) &
    stk::mesh::selectField(yField)
    : (metaData.locally_owned_part() | metaData.globally_shared_part()) &
    stk::mesh::selectField(xField) &
    stk::mesh::selectField(yField);

  stk::mesh::BucketVector const& buckets = bulkData.get_buckets( entityRankValue, selector );

  for(size_t i=0; i < buckets.size(); ++i) {
    stk::mesh::Bucket & b = *buckets[i];
    const stk::mesh::Bucket::size_type length = b.size();
    const size_t fieldSize = field_bytes_per_entity(xField, b) / sizeof(double);
    ThrowAssert(fieldSize == field_bytes_per_entity(yField, b) / sizeof(double));
    const unsigned kmax = length * fieldSize;
    const double * x = (double*)stk::mesh::field_data(xField, b);
    double * y = (double*)stk::mesh::field_data(yField, b);
    for(unsigned k = 0 ; k < kmax ; ++k) {
      y[k] = x[k];
    }
  }

}

void field_index_copy(
  const stk::mesh::MetaData & metaData,
  const stk::mesh::BulkData & bulkData,
  const stk::mesh::FieldBase & xField,
  const int xFieldIndex,
  const stk::mesh::FieldBase & yField,
  const int yFieldIndex,
  const bool auraIsActive,
  const stk::topology::rank_t entityRankValue)
{
  // decide on selector
  const stk::mesh::Selector selector = auraIsActive 
    ? metaData.universal_part() &
    stk::mesh::selectField(xField) &
    stk::mesh::selectField(yField)
    : (metaData.locally_owned_part() | metaData.globally_shared_part()) &
    stk::mesh::selectField(xField) &
    stk::mesh::selectField(yField);

  stk::mesh::BucketVector const& buckets = bulkData.get_buckets( entityRankValue, selector );

  for(size_t i=0; i < buckets.size(); ++i) {
    stk::mesh::Bucket & b = *buckets[i];
    const stk::mesh::Bucket::size_type length = b.size();
    const size_t xFieldSize = field_bytes_per_entity(xField, b) / sizeof(double);
    const size_t yFieldSize = field_bytes_per_entity(yField, b) / sizeof(double);
    const double * x = (double*)stk::mesh::field_data(xField, b);
    double * y = (double*)stk::mesh::field_data(yField, b);
    for(unsigned k = 0 ; k < length ; ++k) {
      y[k*yFieldSize+yFieldIndex] = x[k*xFieldSize+xFieldIndex];
    }
  }

}

} // namespace nalu
} // namespace Sierra
