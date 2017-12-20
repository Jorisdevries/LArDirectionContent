/**
 *  @file   ExampleContent/src/ExampleAlgorithms/ExampleDirectionAlgorithm.cc
 * 
 *  @brief  Implementation of the access lists algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArDirection/ExampleDirectionAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_content
{

StatusCode ExampleDirectionAlgorithm::Run()
{
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListName, pClusterList));

    ClusterVector clusterVector(pClusterList->begin(), pClusterList->end());

    for (const Cluster *const pCluster : clusterVector)
    {    
        if (LArClusterHelper::GetClusterHitType(pCluster) == TPC_VIEW_W)
        {
            TrackDirectionTool::FitResult fitResult = m_pTrackDirectionTool->Run(this, pCluster);
            std::cout << "Probability: " << fitResult.GetProbability() << std::endl;
            std::cout << "Vertex x coordinate: " << fitResult.GetBeginpoint().GetX() << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ExampleDirectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterListName", m_clusterListName));

    AlgorithmTool *pAlgorithmTool(nullptr);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "TrackDirection", pAlgorithmTool));

    if (!(this->m_pTrackDirectionTool = dynamic_cast<TrackDirectionTool *>(pAlgorithmTool)))
        throw STATUS_CODE_FAILURE;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

