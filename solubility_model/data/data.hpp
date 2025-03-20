#ifndef SOL_MODEL_DATA_HPP
#define SOL_MODEL_DATA_HPP 

#include "../../solubility_model/utils/utils.hpp"

#pragma once 

namespace DATA
{   
    using namespace DATATYPE;
    using namespace LOG;
    using namespace CONST;

    enum class PROPERTY { H2LIQUID, H2VAPOR, H2OLIQUID, H2OVAPOR }; 
    enum class PROPERTY_DIMENSION { ROW, COL };

    class Data
    {   
        public:
        Data(Tensor1DFloat64 & pressureData, Tensor1DFloat64 & temperatureData, 
        Tensor2DFloat64 & propertyData);
        
        bool verifyDimension();
        void setPressureDimension( PROPERTY_DIMENSION pressureDim );
        void setTemperatureDimension( PROPERTY_DIMENSION temperatureDim );
        
        Tensor1DFloat64 getPressureData();
        Tensor1DFloat64 getTemperatureData();
        double getPropertyData( double T, double P );

        public:
        Tensor1DFloat64 _pressureData;
        Tensor1DFloat64 _temperatureData;
        Tensor2DFloat64 _propertyData;
        //
        PROPERTY_DIMENSION _pressureDim = PROPERTY_DIMENSION::ROW;
        PROPERTY_DIMENSION _temperatureDim = PROPERTY_DIMENSION::COL;
    };
}

#endif