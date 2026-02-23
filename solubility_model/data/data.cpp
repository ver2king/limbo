#include "../../solubility_model/data/data.hpp"

namespace DATA
{
    Data::Data(Tensor1DFloat64 & pressureData, Tensor1DFloat64 & temperatureData, 
    Tensor2DFloat64 & propertyData): 
    _pressureData(pressureData), _temperatureData(temperatureData),
    _propertyData(propertyData) { };
    
    void Data::setPressureDimension( PROPERTY_DIMENSION pressureDim )
    { _pressureDim = pressureDim; };

    void Data::setTemperatureDimension( PROPERTY_DIMENSION temperatureDim )
    { _temperatureDim = temperatureDim; };

    bool Data::verifyDimension()
    { 
        int rowDim = _propertyData.size();
        int colDim = _propertyData[0].size();
        //
        int pressureDim = _pressureData.size();
        int temperatureDim = _temperatureData.size();
        bool rowDimVerify = false;
        bool colDimVerify = false;
        //
        if ( _pressureDim == PROPERTY_DIMENSION::ROW ) { 
            rowDimVerify = ( pressureDim == rowDim );
        } else { 
            colDimVerify = (pressureDim == colDim );
        };
        if ( _temperatureDim == PROPERTY_DIMENSION::ROW ) { 
            colDimVerify = ( temperatureDim == rowDim );
        } else {
            colDimVerify = ( temperatureDim == colDim );
        };
        //
        if ( rowDimVerify == true && colDimVerify == true ) { 
            return true;
        } else { 
            return false;
        };
    };

    void Data::setPressureData( Tensor1DFloat64 & pressureData )
    { _pressureData = pressureData; };

    void Data::setTemperatureData( Tensor1DFloat64 & temperatureData )
    { _temperatureData = temperatureData; };

    void Data::setPropertyData( Tensor2DFloat64 & propertyData )
    { _propertyData = propertyData; };

    Tensor1DFloat64 Data::getPressureData() { return _pressureData; };
    Tensor1DFloat64 Data::getTemperatureData() { return _temperatureData; };
    Tensor2DFloat64 Data::getPropertyData() { return _propertyData; };

    double Data::getPropertyData( double T, double P ) 
    {   
        //
        double propertyVal = doubleNaN;
        auto temperatureIt = std::find( _temperatureData.begin(), _temperatureData.end(), T );
        auto pressureIt = std::find( _pressureData.begin(), _pressureData.end(), P );
        if (temperatureIt != _temperatureData.end() && pressureIt != _pressureData.end() ) {
            int temperatureIdx = std::distance( _temperatureData.begin(), temperatureIt );
            int pressureIdx = std::distance( _pressureData.begin(), pressureIt );
            if ( _pressureDim == PROPERTY_DIMENSION::ROW && _temperatureDim == PROPERTY_DIMENSION::COL ) 
            { propertyVal = _propertyData[pressureIdx][temperatureIdx]; } 
            else if ( _pressureDim == PROPERTY_DIMENSION::COL && _temperatureDim == PROPERTY_DIMENSION::ROW ) 
            { propertyVal = _propertyData[temperatureIdx][pressureIdx]; }
            else 
            { propertyVal = doubleNaN; };
        } else { 
            propertyVal = doubleNaN;
        };
        return propertyVal;
    };

    ExperimentalSolubilityDataBlock::ExperimentalSolubilityDataBlock( Tensor2DString solBlockRawData ):
    _solBlockRawData( solBlockRawData ) { };

    void ExperimentalSolubilityDataBlock::getSolubilityData()
    {   
        //printf("Raw data of block. \n");
        //PrintTensorString2D( _solBlockRawData );
        int const temperatureDim = 1;
        int const pressureDim = _solBlockRawData.size();
        //
        _temperatureData.resize( temperatureDim );
        _pressureData.resize( pressureDim );
        _propertyData.resize( pressureDim );
        //
        double temperatureVal = std::stod( _solBlockRawData[0][0] );
        _temperatureData[0] = temperatureVal;
        //
        for ( int ip = 0; ip < pressureDim; ++ ip )
        {   
            _propertyData[ip].resize( temperatureDim );
            double pressureVal = std::stod( _solBlockRawData[ip][1] );
            _pressureData[ip] = pressureVal;
            // 
            for ( int it = 0; it < temperatureDim; ++ it )
            {   
                char* endPtr;
                _propertyData[ip][it] = std::stod( _solBlockRawData[ip][2] ) * std::strtod( _solBlockRawData[ip][3].c_str(),
                & endPtr );
            };
        };
        _experimentalData.setPressureData( _pressureData );
        _experimentalData.setTemperatureData( _temperatureData );
        _experimentalData.setPropertyData( _propertyData );
    };

    ExperimentalSolubilityData::ExperimentalSolubilityData( std::ifstream & dataFile, PROPERTY Property ):
    _dataFile( std::move( dataFile ) ), _Property( Property ) { };

    void ExperimentalSolubilityData::setPropertyNamePrefix( String propNamePrefix )
    { _propNamePrefix = propNamePrefix; };

    Tensor1DString ExperimentalSolubilityData::createPropertyNames( int numOfDataBlocks )
    {
        Tensor1DString propNames;
        for ( int i = 0 ; i< numOfDataBlocks; ++ i )
        {
            String propName = std::to_string(i);
            propName = _propNamePrefix + "_" + propName;
            propNames.push_back( propName );
        };
        return propNames;
    };

    void ExperimentalSolubilityData::getPropertyData()
    {
        int numOfDataBlocks = 1;
        int i = 0;
        Tensor1DInt locOfDataBlocks; 
        locOfDataBlocks.push_back( -1 );
        //
        std::regex re("(\\t)");
        String data;
        Tensor2DString rawData;
        //
        if (_dataFile.is_open())
        {
        while ( std::getline(_dataFile, data) )
        {   
            Tensor1DString rawLineData;
            if ( data.size() > 1 )
            {   
                /// FIXME: Some parsing is incorrect here !
                std::sregex_token_iterator first{data.begin(), data.end(), re, -1}; 
                std::sregex_token_iterator last; 
                std::vector<std::string> tokens{first, last};
                for ( std::string str : tokens )
                {   
                    rawLineData.push_back( str );
                };
            } else 
            { 
                numOfDataBlocks ++;
                locOfDataBlocks.push_back( i );
                rawLineData.push_back("\n");
            };
            rawData.push_back( rawLineData );
            i++;
        };
        locOfDataBlocks.push_back( i );
        _dataFile.close();
        } else
        {
        printf("File not opened.");
        };       
        //
        // PrintTensorString2D( rawData );
        std::vector< Tensor2DString > rawDataOfAllBlocks;
        for ( int j = 0; j < numOfDataBlocks; ++ j )
        {
            Tensor2DString rawDataOfOneBlock;
            int dataBlockStart = locOfDataBlocks[j] + 1;
            int dataBlockEnd = locOfDataBlocks[j+1];
            for ( int k = dataBlockStart; k < dataBlockEnd; ++ k )
            { rawDataOfOneBlock.push_back( rawData[k] ); };
            //
            rawDataOfAllBlocks.push_back( rawDataOfOneBlock );
        };
        //
        Tensor1DString allPropNames = createPropertyNames( numOfDataBlocks );
        for ( int j = 0; j < numOfDataBlocks; ++ j )
        {   
            Tensor2DString rawDataOfOneBlock = rawDataOfAllBlocks[j];
            // PrintTensorString2D( rawDataOfOneBlock );
            ExperimentalSolubilityDataBlock expSolDataBlock( rawDataOfOneBlock );
            expSolDataBlock.getSolubilityData();
            //
            //PrintTensor1D<double>(expSolDataBlock._pressureData);
            //PrintTensor1D<double>(expSolDataBlock._temperatureData);
            //PrintTensor2D<double>(expSolDataBlock._propertyData);
            DataIdentifier dataIdentifier = { allPropNames[j], _Property };
            _experimentalSolubilityData.insert( { dataIdentifier, expSolDataBlock._experimentalData } );
        };
    };

    void ExperimentalSolubilityData::getDataIdentifiers()
    {
        for ( auto & e : _experimentalSolubilityData )
        { _dataIdentifiers.push_back( e.first ); };
    };
}