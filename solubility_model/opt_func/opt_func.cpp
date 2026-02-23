#include "opt_func.hpp"

namespace OPT_FUNCTION
{
    ObjectiveFunction::ObjectiveFunction(std::map< DataIdentifier, Data > & modelData,
    std::map< DataIdentifier, Data > & experimentalData ):
    _modelData( modelData ), _experimentalData ( experimentalData ) 
    { 
        std::vector< PROPERTY > modelProps;
        std::vector< PROPERTY > experimentalProps;
        std::vector< String > modelPropNames;
        std::vector< String > experimentalPropNames;
        for ( auto & i : modelData ) { 
            modelProps.push_back( i.first.second ); 
            modelPropNames.push_back( i.first.first );
        };
        for ( auto & j : experimentalData ) { 
            experimentalProps.push_back( j.first.second ); 
            experimentalPropNames.push_back( j.first.first );
        };
        _modelProps = modelProps;
        _experimentalProps = experimentalProps;
        _modelPropNames = modelPropNames;
        _experimentalPropNames = experimentalPropNames;
    };

    void ObjectiveFunction::setObjectiveFunctionType( OBJ_FUNC_TYPE & objFuncType )
    { _objFuncType = objFuncType; };

    bool ObjectiveFunction::verifyDimension()
    {       
        if ( _modelProps.size() != _experimentalProps.size() ) 
        { return false; }
        else {
            DataIdentifier modelInd = { _modelPropNames[0], _modelProps[0] };
            DataIdentifier experimentalInd = { _experimentalPropNames[0], _experimentalProps[0] };
            Tensor2DFloat64 modelPropData = _modelData.at( modelInd )._propertyData;
            Tensor2DFloat64 expPropData = _experimentalData.at( experimentalInd )._propertyData;
            if ( modelPropData.size() == expPropData.size() && modelPropData[0].size() == expPropData[0].size() )
            { return true; }
            else { return false; };
        };
    };

    double ObjectiveFunction::computeSingleProperty( DataIdentifier & dataInd )
    {    
        String propName = dataInd.first;
        if ( FindElementInTensor1D<String>( propName, _modelPropNames ) == false ||
        FindElementInTensor1D<String>( propName, _experimentalPropNames ) == false )
        { 
            printf("Can not find data identifider. \n");
            return doubleInf;  
        }
        else 
        {   
            Tensor2DFloat64 modelPropData = _modelData.at( dataInd )._propertyData;
            Tensor2DFloat64 expPropData = _experimentalData.at( dataInd )._propertyData;
            //
            int rowDim = modelPropData.size();
            int colDim = modelPropData[0].size();
            //
            int totalValid = 0;
            int totalExceedTol = 0;
            double totalVal = 0.;
            for ( int i = 0; i < rowDim; ++ i )
            {
            for ( int  j = 0; j < colDim; ++ j )
            {
                double modelVal = modelPropData[i][j];
                double expVal = expPropData[i][j];
                std::cout << "Model and experiment values are: " << "\n";
                std::cout << modelVal << " , " << expVal << "\n";
                if ( modelVal < 0 || expVal < 0 ) { 
                    totalValid += 0; 
                    totalVal += 0.; 
                } else {
                    double objVal = 0.;
                    totalValid ++;
                    if ( _objFuncType == OBJ_FUNC_TYPE::MSE )
                    { objVal = std::pow( ( modelVal - expVal ), 2 ); }
                    else if ( _objFuncType == OBJ_FUNC_TYPE::RMSE )
                    { objVal = std::abs( modelVal - expVal ); }
                    else if ( _objFuncType == OBJ_FUNC_TYPE::RE )
                    { objVal = std::abs( ( modelVal - expVal ) / expVal ); }
                    else { objVal = doubleNaN; }
                    //
                    totalVal += objVal;
                    if ( objVal > _objFuncTol ) { 
                        totalExceedTol ++;
                    };
                };
            };
            };
            //
            //__VAR_WITH_EXCEPTION__(totalValid, "TOTAL NUMBER OF OBJECTIVE POINTS.");
            //__VAR_WITH_EXCEPTION__(totalExceedTol, "NUMBER OF OBJECTIVE VALUES THAT EXCEED TOLERANCE.");
            //
            if ( totalVal == doubleNaN )
            { 
                return doubleNaN;
            } else if ( totalValid == 0 )
            { 
                return doubleNaN; 
            }
            else { 
                return totalVal / totalValid; 
            }; 
        };
    };

    double ObjectiveFunction::computeMultipleProperties( std::vector< DataIdentifier > & dataInds )
    {
        double objVal = 0;
        int numOfDataInds = dataInds.size();
        for ( int i = 0; i < numOfDataInds; ++ i )
        {
            objVal += computeSingleProperty( dataInds[i] );
        };
        return objVal;
    }

    ProblemParams::ProblemParams(  std::map< PROBLEM_UNIT, String > _problemUnits,
    Tensor1DFloat64 _temperatureData, Tensor1DFloat64 _pressureData, 
    Tensor1DString _paramNames, std::map< String, Tensor1DFloat64 > _allParamsBounds,
    std::vector< PROPERTY > _allProps, std::vector < String > _allPropNames )
    {
        problemUnits = _problemUnits;
        temperatureData = _temperatureData;
        pressureData = _pressureData;
        paramNames = _paramNames;
        allParamsBounds = _allParamsBounds;
        allProps = _allProps;
        allPropNames =  _allPropNames;
    };

    Tensor1DFloat64 ProblemParams::getTemperatureData()
    { return temperatureData; };

    Tensor1DFloat64 ProblemParams::getPressureData()
    { return pressureData; };

    String ProblemParams::getTemperatureUnit()
    { return problemUnits.at( PROBLEM_UNIT::T ); };

    String ProblemParams::getPressureUnit()
    { return problemUnits.at( PROBLEM_UNIT::P ); };

    Tensor1DFloat64 ProblemParams::getParamBounds(String & paramName )
    { return allParamsBounds.at( paramName ); };

    ModelPrediction::ModelPrediction( ProblemParams & _problemParams ):
    problemParams( _problemParams ) { };
}