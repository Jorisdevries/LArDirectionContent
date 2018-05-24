/**
 *  @file   lardirectioncontent/LArContent.cc
 *
 *  @brief  Factory implementations for content intended for use with particle flow reconstruction at liquid argon time projection chambers
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"
#include "Pandora/Pandora.h"

#include "LArDirection/ExampleDirectionAlgorithm.h"
#include "LArDirection/DirectionAnalysisAlgorithm.h"
#include "LArDirection/DirectionClusterSplittingAlgorithm.h"

#include "LArDirection/TrackDirectionTool.h"
#include "LArDirection/DirectionFlowProbabilityTool.h"

#include "LArDirectionContent.h"

#define LAR_DIRECTION_ALGORITHM_LIST(d)                                                                                                   \
    d("LArExampleDirection",                  ExampleDirectionAlgorithm)                                                                  \
    d("LArDirectionAnalysis",                 DirectionAnalysisAlgorithm)                                                                 \
    d("LArDirectionClusterSplitting",         DirectionClusterSplittingAlgorithm)

#define LAR_DIRECTION_ALGORITHM_TOOL_LIST(d)                                                                                              \
    d("LArTrackDirectionTool",                TrackDirectionTool)                                                                         \
    d("LArDirectionFlowProbabilityTool",      DirectionFlowProbabilityTool)

#define FACTORY Factory

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

#define LAR_DIRECTION_CONTENT_CREATE_ALGORITHM_FACTORY(a, b)                                                                    \
class b##FACTORY : public pandora::AlgorithmFactory                                                                             \
{                                                                                                                               \
public:                                                                                                                         \
    pandora::Algorithm *CreateAlgorithm() const {return new b;};                                                                \
};

LAR_DIRECTION_ALGORITHM_LIST(LAR_DIRECTION_CONTENT_CREATE_ALGORITHM_FACTORY)

//------------------------------------------------------------------------------------------------------------------------------------------

#define LAR_DIRECTION_CONTENT_CREATE_ALGORITHM_TOOL_FACTORY(a, b)                                                                         \
class b##FACTORY : public pandora::AlgorithmToolFactory                                                                         \
{                                                                                                                               \
public:                                                                                                                         \
    pandora::AlgorithmTool *CreateAlgorithmTool() const {return new b;};                                                        \
};

LAR_DIRECTION_ALGORITHM_TOOL_LIST(LAR_DIRECTION_CONTENT_CREATE_ALGORITHM_TOOL_FACTORY)

} // namespace lar_content

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

#define LAR_DIRECTION_CONTENT_REGISTER_ALGORITHM(a, b)                                                                                    \
{                                                                                                                               \
    const pandora::StatusCode statusCode(PandoraApi::RegisterAlgorithmFactory(pandora, a, new lar_content::b##FACTORY));        \
    if (pandora::STATUS_CODE_SUCCESS != statusCode)                                                                             \
        return statusCode;                                                                                                      \
}

#define LAR_DIRECTION_CONTENT_REGISTER_ALGORITHM_TOOL(a, b)                                                                               \
{                                                                                                                               \
    const pandora::StatusCode statusCode(PandoraApi::RegisterAlgorithmToolFactory(pandora, a, new lar_content::b##FACTORY));    \
    if (pandora::STATUS_CODE_SUCCESS != statusCode)                                                                             \
        return statusCode;                                                                                                      \
}

pandora::StatusCode LArDirectionContent::RegisterAlgorithms(const pandora::Pandora &pandora)
{
    LAR_DIRECTION_ALGORITHM_LIST(LAR_DIRECTION_CONTENT_REGISTER_ALGORITHM);
    LAR_DIRECTION_ALGORITHM_TOOL_LIST(LAR_DIRECTION_CONTENT_REGISTER_ALGORITHM_TOOL);
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
