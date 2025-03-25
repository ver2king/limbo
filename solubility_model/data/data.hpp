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
    typedef std::pair< String, PROPERTY > DataIdentifier;

    class Data
    {   
        public:
        Data() = default;
        Data(Tensor1DFloat64 & pressureData, Tensor1DFloat64 & temperatureData, 
        Tensor2DFloat64 & propertyData);
        
        bool verifyDimension();
        void setPressureDimension( PROPERTY_DIMENSION pressureDim );
        void setTemperatureDimension( PROPERTY_DIMENSION temperatureDim );
        
        void setPressureData( Tensor1DFloat64 & pressureData );
        void setTemperatureData( Tensor1DFloat64 & temperatureData );
        void setPropertyData( Tensor2DFloat64 & propertyData );

        Tensor1DFloat64 getPressureData();
        Tensor1DFloat64 getTemperatureData();
        Tensor2DFloat64 getPropertyData();

        double getPropertyData( double T, double P );

        public:
        Tensor1DFloat64 _pressureData;
        Tensor1DFloat64 _temperatureData;
        Tensor2DFloat64 _propertyData;
        //
        PROPERTY_DIMENSION _pressureDim = PROPERTY_DIMENSION::ROW;
        PROPERTY_DIMENSION _temperatureDim = PROPERTY_DIMENSION::COL;
    };

    class ExperimentalSolubilityDataBlock
    {
        public:
        ExperimentalSolubilityDataBlock( Tensor2DString solBlockRawData );
        void getSolubilityData();

        public:
        Tensor2DString _solBlockRawData = { { "" } };
        //
        Tensor1DFloat64 _temperatureData = std::vector<double>( 1, doubleNaN );
        Tensor1DFloat64 _pressureData = std::vector<double>( 1, doubleNaN );
        Tensor2DFloat64 _propertyData = { { doubleNaN } };
        Data _experimentalData;
    };

    class ExperimentalSolubilityData
    {
        /// TODO: Implement method(s) to extract solubility data 
        /// and reform as matrix to feed into Data class
        public:
        ExperimentalSolubilityData( std::ifstream & dataFile, PROPERTY Property );

        void setPropertyNamePrefix( String propNamePrefix );

        Tensor1DString createPropertyNames( int numOfDataBlocks );
        void getPropertyData();
        void getDataIdentifiers();
    
        public:
        std::ifstream _dataFile;
        PROPERTY _Property;
        String _propNamePrefix;
        std::map < DataIdentifier, Data > _experimentalSolubilityData;
        std::vector < DataIdentifier > _dataIdentifiers;
    };
}

#endif