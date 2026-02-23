#ifndef SOL_MODEL_OPT_FUNC_HPP
#define SOL_MODEL_OPT_FUNC_HPP

#include "../../solubility_model/src/model.hpp"
#include "../../src/limbo/tools/macros.hpp"

#pragma once

namespace OPT_FUNCTION
{   
    using namespace MODEL;
    using namespace LOG;
    using namespace DATA;
    
    enum class PROBLEM_UNIT { T, P };
    enum class OBJ_FUNC_TYPE { MSE, RMSE, RE };

    class ObjectiveFunction
    { 
        public:
        ObjectiveFunction(std::map< DataIdentifier, Data > & modelData,
        std::map< DataIdentifier, Data > & experimentalData );
        
        void setObjectiveFunctionType( OBJ_FUNC_TYPE & objFuncType );

        bool verifyDimension();
        double computeSingleProperty( DataIdentifier & dataInd );
        double computeMultipleProperties( std::vector< DataIdentifier > & dataInds );

        public:
        std::map< DataIdentifier, Data > _modelData;
        std::map< DataIdentifier, Data > _experimentalData;
        std::vector< PROPERTY > _modelProps;
        std::vector< PROPERTY > _experimentalProps;
        std::vector< String > _modelPropNames;
        std::vector< String > _experimentalPropNames;
        //
        OBJ_FUNC_TYPE _objFuncType = OBJ_FUNC_TYPE::RE;
        double const _objFuncTol = 0.05;
    };

    struct ProblemParams 
    {
        std::map< PROBLEM_UNIT, String > problemUnits;
        Tensor1DFloat64 temperatureData;
        Tensor1DFloat64 pressureData;
        double saltMolarity = 0.;
        //
        Tensor1DString paramNames;
        std::map< String, Tensor1DFloat64 > allParamsBounds;
        //
        std::vector< PROPERTY > allProps;
        std::vector< String > allPropNames;
        //
        Tensor1DFloat64 getTemperatureData();
        Tensor1DFloat64 getPressureData();
        String getTemperatureUnit();
        String getPressureUnit();
        //
        ProblemParams(  std::map< PROBLEM_UNIT, String > _problemUnits,
        Tensor1DFloat64 _temperatureData,
        Tensor1DFloat64 _pressureData, 
        Tensor1DString _paramNames,
        std::map< String, Tensor1DFloat64 > _allParamsBounds,
        std::vector< PROPERTY > _allProps,
        std::vector < String > _allPropNames );

        Tensor1DFloat64 getParamBounds( String & paramName );
    };

    struct ModelPrediction
    {
        ProblemParams problemParams;
        PROPERTY_DIMENSION pressureDim = PROPERTY_DIMENSION::ROW;
        PROPERTY_DIMENSION temperatureDim = PROPERTY_DIMENSION::COL;

        ModelPrediction( ProblemParams & _problemParams );

        std::map< DataIdentifier, Data > operator()( const Eigen::VectorXd& inputParams )
        {
            std::map< DataIdentifier, Data > modelData;
            //
            Tensor1DFloat64 temperatureData = problemParams.temperatureData;
            Tensor1DFloat64 pressureData = problemParams.pressureData;
            double saltMolarity = problemParams.saltMolarity;
            String temperatureUnit = problemParams.getTemperatureUnit();
            String pressureUnit = problemParams.getPressureUnit();
            //
            Tensor1DFloat64 param0Bounds = problemParams.getParamBounds( problemParams.paramNames[0] );
            Tensor1DFloat64 param1Bounds = problemParams.getParamBounds( problemParams.paramNames[1] );
            Tensor1DFloat64 param2Bounds = problemParams.getParamBounds( problemParams.paramNames[2] );
            Tensor1DFloat64 param3Bounds = problemParams.getParamBounds( problemParams.paramNames[3] );
            //
            ModelParams modelParams = { LinearScaler(inputParams(0), param0Bounds[0], param0Bounds[1]),
            LinearScaler(inputParams(1), param1Bounds[0], param1Bounds[1]),
            LinearScaler(inputParams(2), param2Bounds[0], param2Bounds[1]),
            LinearScaler(inputParams(3), param3Bounds[0], param3Bounds[1]) };
            //
            int numOfProps = problemParams.allPropNames.size();
            for ( int i = 0; i < numOfProps; ++ i)
            {   
                PROPERTY Property = problemParams.allProps[i];
                String propertyName = problemParams.allPropNames[i];
                // 
                Tensor2DFloat64 propertyDataMatrix = SolubilityModelWrapper(temperatureData, pressureData, 
                modelParams, saltMolarity, temperatureUnit, pressureUnit, 
                Property, pressureDim, temperatureDim);
                DATA::Data propertyData( pressureData, temperatureData, propertyDataMatrix );
                DataIdentifier dataInd = { propertyName, Property };
                modelData.insert( { dataInd, propertyData } );
            };
            return modelData;
        };
    };

}

#endif