#include "../../solubility_model/src/model.hpp"
#include "../../solubility_model/data/data.hpp"

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
        ObjectiveFunction(std::map< PROPERTY, Data > & modelData,
        std::map< PROPERTY, Data > & experimentalData );
        
        void setObjectiveFunctionType( OBJ_FUNC_TYPE & objFuncType );

        bool verifyDimension();
        double computeSingleProperty( PROPERTY & propName );
        double computeMultipleProperties( std::vector< PROPERTY > & propNames );

        public:
        std::map< PROPERTY, Data > _modelData;
        std::map< PROPERTY, Data > _experimentalData;
        std::vector< PROPERTY > _modelProps;
        std::vector< PROPERTY > _experimentalProps;
        //
        OBJ_FUNC_TYPE _objFuncType = OBJ_FUNC_TYPE::MSE;
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
        std::vector< PROPERTY > propNames;
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
        std::vector< PROPERTY > _propNames );

        Tensor1DFloat64 getParamBounds( String & paramName );
    };

    struct ModelPrediction
    {
        ProblemParams problemParams;
        PROPERTY_DIMENSION pressureDim = PROPERTY_DIMENSION::ROW;
        PROPERTY_DIMENSION temperatureDim = PROPERTY_DIMENSION::COL;

        ModelPrediction( ProblemParams & _problemParams );

        std::map< PROPERTY, Data > operator()( const Eigen::VectorXd& inputParams )
        {
            std::map< PROPERTY, Data > modelData;
            //
            __CHECK_POINT_WITH_MSG__("GET DATA");

            Tensor1DFloat64 temperatureData = problemParams.temperatureData;
            Tensor1DFloat64 pressureData = problemParams.pressureData;
            double saltMolarity = problemParams.saltMolarity;
            String temperatureUnit = problemParams.getTemperatureUnit();
            String pressureUnit = problemParams.getPressureUnit();
            //
            __CHECK_POINT_WITH_MSG__("GET PARAM BOUNDS");

            Tensor1DFloat64 param0Bounds = problemParams.getParamBounds( problemParams.paramNames[0] );
            Tensor1DFloat64 param1Bounds = problemParams.getParamBounds( problemParams.paramNames[1] );
            Tensor1DFloat64 param2Bounds = problemParams.getParamBounds( problemParams.paramNames[2] );
            Tensor1DFloat64 param3Bounds = problemParams.getParamBounds( problemParams.paramNames[3] );
            //
            __CHECK_POINT_WITH_MSG__("SET MODEL PARAMS");
            
            ModelParams modelParams = { LinearScaler(inputParams(0), param0Bounds[0], param0Bounds[1]),
            LinearScaler(inputParams(1), param1Bounds[0], param1Bounds[1]),
            LinearScaler(inputParams(2), param2Bounds[0], param2Bounds[1]),
            LinearScaler(inputParams(3), param3Bounds[0], param3Bounds[1]) };
            //
            __CHECK_POINT_WITH_MSG__("COMPUTE OBJECTIVE FUNCTION");

            for ( int i = 0; i < problemParams.propNames.size(); ++ i)
            {   
                __CHECK_POINT_WITH_MSG__("GET PROPERTY");

                PROPERTY Property = problemParams.propNames[i];
                
                __CHECK_POINT_WITH_MSG__("COMPUTE MODEL PRED");
                Tensor2DFloat64 propertyDataMatrix = SolubilityModelWrapper(temperatureData, pressureData, 
                modelParams, saltMolarity, temperatureUnit, pressureUnit, 
                Property, pressureDim, temperatureDim);

                __CHECK_POINT_WITH_MSG__("SET MODEL DATA");
                DATA::Data propertyData( pressureData, temperatureData, propertyDataMatrix );
                std::cout << propertyDataMatrix.size() << " , " << propertyDataMatrix[0].size() << "\n";
                
                __CHECK_POINT_WITH_MSG__("INSERT MODEL DATA");
                modelData.insert( { Property, propertyData } );
            };
            return modelData;
        };
    };

}