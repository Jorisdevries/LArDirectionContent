/**
 *  @file   ExampleContent/src/ExampleAlgorithms/ExampleDirectionAlgorithm.cc
 * 
 *  @brief  Implementation of the access lists algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "lardirectioncontent/LArDirection/ExampleDirectionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ExampleDirectionAlgorithm::Run()
{
    const ClusterList *pCurrentClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCurrentClusterList));

    std::cout << "Test." << std::endl;

    // Access to a named list is demonstrated below. This access only proceeds if a list name has been specified within the
    // algorithm xml tags. This list name is an optional parameter, rather than mandatory.
    if (!m_requestedCaloHitListName.empty())
    {
        const CaloHitList *pNamedCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_requestedCaloHitListName, pNamedCaloHitList));

        // Use the named calo hit list...
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ExampleDirectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RequestedCaloHitListName", m_requestedCaloHitListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

