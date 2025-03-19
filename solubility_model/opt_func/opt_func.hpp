#include "../../solubility_model/src/model.hpp"
#include "../../src/limbo/tools/macros.hpp"

#pragma once

namespace OPT_FUNCTION
{   
    using namespace MODEL;
    
    int const numberOfParams = 4;
    int const numberOfOpt = 1;

    struct OptimFunc {
        BO_PARAM(size_t, dim_in, numberOfParams);
        BO_PARAM(size_t, dim_out, numberOfOpt);
    
        Eigen::VectorXd operator()(const Eigen::VectorXd& inputParams) const
        {   
            /// Initial set-up data for the solubility model
            String unitT = "K";
            String unitP = "Bar";
            double m_s = 0.;
            Tensor1DFloat64 dataT = {298.15, 323.15, 348.15, 373.15, 398.15, 423.15};
            Tensor1DFloat64 dataP;
            double initP = 2.;
            while ( initP <= 200. )
            {
                dataP.push_back(initP);
                if ( initP < 30 ) { initP += 1.; }
                else { initP += 5.; };
            };
            std::ifstream f("/home/ver2king/Desktop/limbo/solubility_model/h2_h2o_yh2o.txt");
            int row = dataP.size(); int col = dataT.size();
            Tensor2DFloat64 waterVaporPhaseData = ReadMatrixFile(f, row, col);
            assert( waterVaporPhaseData.size() == row);
            assert( waterVaporPhaseData[0].size() == col);

            /// Define the lower and higher bounds for all input coeffients
            double lowerRate = 0.95; double higherRate = 1.05; 
            Tensor1DFloat64 h2ACoeffBounds = { 1441753.379 * lowerRate, 1441753.379 * higherRate };
            Tensor1DFloat64 h2BCoeffBounds = { 18.417 * lowerRate, 18.417 * higherRate };
            Tensor1DFloat64 h2oACoeffBounds = { 142666655.8 * lowerRate, 142666655.8 * higherRate };
            Tensor1DFloat64 h2oBCoeffBounds = { 21.127 * lowerRate, 21.127 * higherRate };

            /// Parse the coefficients 
            ModelParams modelParams = { LinearScaler(inputParams(0), h2ACoeffBounds[0], h2ACoeffBounds[1]),
            LinearScaler(inputParams(1), h2BCoeffBounds[0], h2BCoeffBounds[1]),
            LinearScaler(inputParams(2), h2oACoeffBounds[0], h2oACoeffBounds[1]),
            LinearScaler(inputParams(3), h2oBCoeffBounds[0], h2oBCoeffBounds[1]) };

            /// Compute the objective value
            double avgRelativeErr = SolubilityModelRelativeError(dataT, dataP, modelParams, m_s, unitT, unitP, waterVaporPhaseData);

            /// Define the final objective value
            Eigen::VectorXd res(1);
            res(0) = - avgRelativeErr;
            return res;
        };
    };
}